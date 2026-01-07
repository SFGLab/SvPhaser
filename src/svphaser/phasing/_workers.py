"""svphaser.phasing._workers
=========================
Worker-process code.

Step B update (biological correctness):
- `n1/n2` are now *ALT-supporting read counts* per haplotype,
  not raw overlap/coverage.
- Evidence is SV-type-aware:
  - DEL: large CIGAR 'D' spanning breakpoints (and optional split-read via SA)
  - INS: large CIGAR 'I' near POS
  - BND: split-read via SA linking to partner chrom:pos
  - INV: split-read via SA to the END breakpoint with strand flip

This is designed to match what IGV "SV-support" reads typically show,
so counts will be closer to the 5/8 style numbers you observed (instead of 27/30).
"""

from __future__ import annotations

import re
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pandas as pd
import pysam
from cyvcf2 import Reader, Variant  # type: ignore

from .algorithms import classify_haplotype
from .types import WorkerOpts

__all__ = ["_phase_chrom_worker"]

# Default breakpoint window (bp) used for evidence gathering.
DEFAULT_BP_WINDOW = 100

# Require event size in the read CIGAR to be at least this fraction of SVLEN.
MIN_CIGAR_FRACTION = 0.30

# And at least this many bp (avoid counting tiny indels / alignment noise).
MIN_CIGAR_BP = 30

# ALT for BND often contains partner like "]chr3:198172833]N" or "N]chr5:181462057]".
_BND_RE = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")


def _has_tabix_index(vcf_path: Path) -> bool:
    """Return True if <file>.tbi or <file>.csi exists."""
    return (
        vcf_path.with_suffix(vcf_path.suffix + ".tbi").exists()
        or vcf_path.with_suffix(vcf_path.suffix + ".csi").exists()
    )


def _coerce_int(x: Any) -> int | None:
    """Convert cyvcf2 INFO values to int if possible."""
    if x is None:
        return None
    if isinstance(x, (list, tuple)):
        return _coerce_int(x[0]) if x else None
    try:
        return int(x)
    except Exception:
        return None


def _svlen_from_record(rec: Variant, pos: int, end: int) -> int:
    svlen = _coerce_int(rec.INFO.get("SVLEN"))
    if svlen is None:
        # VCF END is 1-based inclusive for DEL; length approx END-POS+1
        return abs(end - pos) + 1
    return abs(svlen)


def _parse_bnd_partner(alt: str, rec: Variant) -> tuple[str | None, int | None]:
    """Return (chr2, pos2) for BND if possible."""
    m = _BND_RE.search(alt)
    if m:
        return m.group(1), int(m.group(2))
    chr2 = rec.INFO.get("CHR2")
    # Some callers store partner position in INFO; Sniffles2 uses ALT string primarily.
    return (str(chr2), None) if chr2 else (None, None)


def _parse_sa_tag(read: pysam.AlignedSegment) -> list[tuple[str, int, str]]:
    """Parse SA tag into list of (rname, pos1, strand). pos1 is 1-based."""
    if not read.has_tag("SA"):
        return []
    sa_raw = read.get_tag("SA")
    out: list[tuple[str, int, str]] = []
    for entry in str(sa_raw).split(";"):
        if not entry:
            continue
        # rname,pos,strand,cigar,mapq,nm
        parts = entry.split(",")
        if len(parts) < 3:
            continue
        rname = parts[0]
        try:
            pos1 = int(parts[1])
        except ValueError:
            continue
        strand = parts[2]
        out.append((rname, pos1, strand))
    return out


def _iter_candidate_reads(
    bam: pysam.AlignmentFile,
    chrom: str,
    regions_1based: list[tuple[int, int]],
) -> Iterable[pysam.AlignedSegment]:
    """Yield reads from multiple 1-based-inclusive regions, allowing overlap."""
    for start1, end1 in regions_1based:
        # pysam fetch: 0-based start, end-exclusive
        start0 = max(0, start1 - 1)
        end0 = max(0, end1)  # end1 is inclusive 1-based, so end-exclusive 0-based is end1
        for read in bam.fetch(chrom, start0, end0):
            if read.is_unmapped or read.is_secondary:
                continue
            yield read


def _supports_del(
    read: pysam.AlignedSegment,
    *,
    pos0: int,
    end_excl0: int,
    svlen: int,
    bp_window: int,
) -> bool:
    """DEL support: a large CIGAR deletion spanning POS..END (within window).

    pos0: 0-based start breakpoint
    end_excl0: 0-based exclusive end breakpoint (VCF END maps naturally to this)
    """
    if read.cigartuples is None:
        return False

    min_len = max(MIN_CIGAR_BP, int(MIN_CIGAR_FRACTION * svlen))

    ref = read.reference_start
    for op, length in read.cigartuples:
        if op in (0, 7, 8):
            ref += length
        elif op in (2, 3):
            # deletion / ref skip
            if op == 2 and length >= min_len:
                del_start = ref
                del_end = ref + length
                if abs(del_start - pos0) <= bp_window and abs(del_end - end_excl0) <= bp_window:
                    return True
            ref += length
        elif op in (1, 4, 5, 6):
            # insertion or clipping: does not consume reference
            continue

    # Fallback: split-read evidence (SA) hitting END region
    for rname, sa_pos1, _strand in _parse_sa_tag(read):
        if rname != read.reference_name:
            continue
        sa_pos0 = sa_pos1 - 1
        if abs(sa_pos0 - (end_excl0 - 1)) <= bp_window:
            return True

    return False


def _supports_ins(
    read: pysam.AlignedSegment,
    *,
    pos0: int,
    svlen: int,
    bp_window: int,
) -> bool:
    """INS support: a large CIGAR insertion near POS (within window)."""
    if read.cigartuples is None:
        return False

    min_len = max(MIN_CIGAR_BP, int(MIN_CIGAR_FRACTION * svlen))

    ref = read.reference_start
    for op, length in read.cigartuples:
        if op in (0, 7, 8):
            ref += length
        elif op == 1:
            if length >= min_len and abs(ref - pos0) <= bp_window:
                return True
        elif op in (2, 3):
            ref += length
        else:
            continue
    return False


def _supports_bnd(
    read: pysam.AlignedSegment,
    *,
    pos0: int,
    chr2: str,
    pos2_1based: int,
    bp_window: int,
) -> bool:
    """BND support: SA tag links to partner chrom:pos."""
    # Must overlap POS window (caller already fetched around POS, but keep conservative)
    if abs(read.reference_start - pos0) > 10 * bp_window:
        return False

    pos2_0 = pos2_1based - 1
    for rname, sa_pos1, _strand in _parse_sa_tag(read):
        if rname != chr2:
            continue
        if abs((sa_pos1 - 1) - pos2_0) <= bp_window:
            return True
    return False


def _supports_inv(
    read: pysam.AlignedSegment,
    *,
    pos0: int,
    end0: int,
    bp_window: int,
) -> bool:
    """INV support: SA to the other breakpoint on same chrom with strand flip."""
    strand_primary = "-" if read.is_reverse else "+"

    for rname, sa_pos1, sa_strand in _parse_sa_tag(read):
        if rname != read.reference_name:
            continue
        sa_pos0 = sa_pos1 - 1
        if abs(sa_pos0 - end0) <= bp_window and sa_strand != strand_primary:
            return True
        if abs(sa_pos0 - pos0) <= bp_window and sa_strand != strand_primary:
            return True
    return False


def _read_supports_variant(
    read: pysam.AlignedSegment,
    svtype: str,
    *,
    pos0: int,
    end_excl0: int,
    svlen: int,
    bp_window: int,
    chr2: str | None = None,
    pos2: int | None = None,
) -> bool:
    """Return True if the read provides evidence for the given SV type.

    This is a thin wrapper around the type-specific support helpers so the
    main counting routine remains small and easier to read/test.
    """
    if svtype == "DEL":
        return _supports_del(read, pos0=pos0, end_excl0=end_excl0, svlen=svlen, bp_window=bp_window)
    if svtype == "INS":
        return _supports_ins(read, pos0=pos0, svlen=svlen, bp_window=bp_window)
    if svtype == "BND" and chr2 and pos2:
        return _supports_bnd(
            read,
            pos0=pos0,
            chr2=str(chr2),
            pos2_1based=int(pos2),
            bp_window=bp_window,
        )
    if svtype == "INV":
        return _supports_inv(read, pos0=pos0, end0=(end_excl0 - 1), bp_window=bp_window)

    # Unknown / other SVs: conservative fallback – look for large indel near POS.
    return _supports_ins(read, pos0=pos0, svlen=svlen, bp_window=bp_window) or _supports_del(
        read, pos0=pos0, end_excl0=end_excl0, svlen=svlen, bp_window=bp_window
    )


def _count_hp_sv_support(
    bam: pysam.AlignmentFile,
    chrom: str,
    rec: Variant,
    *,
    bp_window: int = DEFAULT_BP_WINDOW,
) -> tuple[int, int, int, str]:  # noqa: C901
    """Return (n1, n2, sv_end, alt_str) where n1/n2 are ALT-support counts."""

    pos1 = int(rec.POS)
    sv_end = int(rec.end) if getattr(rec, "end", None) is not None else pos1
    alt = ",".join(rec.ALT) if rec.ALT else "<N>"
    svtype = str(rec.INFO.get("SVTYPE", "NA"))

    svlen = _svlen_from_record(rec, pos1, sv_end)

    # Build regions (1-based inclusive) to fetch candidates.
    regions: list[tuple[int, int]] = [(max(1, pos1 - bp_window), pos1 + bp_window)]
    if svtype in {"DEL", "INV"} and sv_end != pos1:
        regions.append((max(1, sv_end - bp_window), sv_end + bp_window))

    # For BND we need the partner locus
    chr2 = None
    pos2 = None
    if svtype == "BND":
        chr2, pos2 = _parse_bnd_partner(alt, rec)

    # Deduplicate by query name: count each read once if it supports ALT.
    # Track hp per read name (in case different segments appear).
    state: dict[str, dict[str, Any]] = {}

    pos0 = pos1 - 1
    end_excl0 = sv_end  # see note in _supports_del: END maps to 0-based exclusive

    for read in _iter_candidate_reads(bam, chrom, regions):
        qn = read.query_name
        if qn is None:
            continue
        st = state.setdefault(qn, {"hp": None, "support": False})

        if st["hp"] is None and read.has_tag("HP"):
            st["hp"] = read.get_tag("HP")

        # If we already know this read supports, we can skip extra work
        if st["support"] is True:
            continue

        if _read_supports_variant(
            read,
            svtype,
            pos0=pos0,
            end_excl0=end_excl0,
            svlen=svlen,
            bp_window=bp_window,
            chr2=chr2,
            pos2=pos2,
        ):
            st["support"] = True

    n1 = n2 = 0
    for st in state.values():
        if not st["support"]:
            continue
        hp = st["hp"]
        if hp == 1:
            n1 += 1
        elif hp == 2:
            n2 += 1

    return n1, n2, sv_end, alt


def _phase_chrom_worker(
    chrom: str,
    vcf_path: Path,
    bam_path: Path,
    opts: WorkerOpts,
) -> pd.DataFrame:
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    rdr = Reader(str(vcf_path))

    rows: list[dict[str, object]] = []

    # Try fast random access first, fall back to linear scan if needed
    use_region_iter = _has_tabix_index(vcf_path)
    records_iter = (
        rdr(f"{chrom}") if use_region_iter else (rec for rec in rdr if rec.CHROM == chrom)
    )

    for rec in records_iter:  # type: ignore[arg-type]
        assert isinstance(rec, Variant)

        n1, n2, sv_end, alt = _count_hp_sv_support(bam, chrom, rec)

        # Respect caller genotype for homozygous ALT: do not infer 1|1 from balance.
        is_hom_alt = False
        if rec.genotypes:
            # cyvcf2: [a1, a2, phased]
            a1, a2 = rec.genotypes[0][0], rec.genotypes[0][1]
            is_hom_alt = a1 == 1 and a2 == 1

        if is_hom_alt:
            gt, gq = "1|1", 0
        else:
            gt, gq = classify_haplotype(
                n1,
                n2,
                min_support=opts.min_support,
                major_delta=opts.major_delta,
                equal_delta=opts.equal_delta,
            )

        row: dict[str, object] = dict(
            chrom=chrom,
            pos=int(rec.POS),
            end=int(sv_end),
            id=rec.ID or ".",
            alt=alt,
            svtype=rec.INFO.get("SVTYPE", "NA"),
            n1=n1,
            n2=n2,
            gt=gt,
            gq=gq,
            gq_label=None,
        )

        if opts.gq_bins:
            for thr, label in opts.gq_bins:
                if gq >= thr:
                    row["gq_label"] = label
                    break

        rows.append(row)

    rdr.close()
    bam.close()
    return pd.DataFrame(rows)
