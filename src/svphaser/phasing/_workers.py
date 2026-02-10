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

from .algorithms import classify_haplotype_v211
from .types import WorkerOpts

__all__ = ["_phase_chrom_worker"]

DEFAULT_BP_WINDOW = 100
MIN_CIGAR_FRACTION = 0.20
MIN_CIGAR_BP = 30

_BND_RE = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")


def _has_tabix_index(vcf_path: Path) -> bool:
    return (
        vcf_path.with_suffix(vcf_path.suffix + ".tbi").exists()
        or vcf_path.with_suffix(vcf_path.suffix + ".csi").exists()
    )


def _coerce_int(x: Any) -> int | None:
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
        return abs(end - pos) + 1
    return abs(svlen)


def _max_abs_int(x: Any) -> int:
    if x is None:
        return 0
    if isinstance(x, (list, tuple)):
        vals = [v for v in (_coerce_int(i) for i in x) if v is not None]
        return max((abs(v) for v in vals), default=0)
    v = _coerce_int(x)
    return abs(v) if v is not None else 0


def _ins_support_min_len(svlen: int) -> int:
    return max(MIN_CIGAR_BP, min(int(MIN_CIGAR_FRACTION * max(svlen, 1)), 200))


def _parse_rnames(val: Any) -> set[str]:
    if val is None:
        return set()
    if isinstance(val, (list, tuple)):
        out: set[str] = set()
        for x in val:
            if x is None:
                continue
            s = str(x)
            if "," in s:
                out.update({p for p in s.split(",") if p})
            else:
                out.add(s)
        return out
    s = str(val)
    if not s:
        return set()
    if "," in s:
        return {p for p in s.split(",") if p}
    return {s}


def _compute_windows(
    rec: Variant, pos1: int, end1: int, svlen: int, *, opts: WorkerOpts
) -> tuple[int, int]:
    if not opts.dynamic_window:
        return max(1, int(opts.bp_window)), max(1, int(opts.bp_window))

    stdev_pos = _coerce_int(rec.INFO.get("STDEV_POS")) or 0
    cipos = _max_abs_int(rec.INFO.get("CIPOS"))
    ciend = _max_abs_int(rec.INFO.get("CIEND"))
    ci = max(cipos, ciend)

    bp_tol = max(int(opts.bp_window), 3 * stdev_pos, ci)
    bp_tol = max(50, min(bp_tol, 1000))

    base = max(200, 2 * bp_tol)
    len_term = int(0.5 * max(svlen, 1))
    fetch_w = max(base, min(len_term, 5000))

    fetch_w = int(max(200, min(fetch_w, 5000)))
    return fetch_w, bp_tol


def _parse_bnd_partner(alt: str, rec: Variant) -> tuple[str | None, int | None]:
    m = _BND_RE.search(alt)
    if m:
        return m.group(1), int(m.group(2))
    chr2 = rec.INFO.get("CHR2")
    return (str(chr2), None) if chr2 else (None, None)


def _parse_sa_tag(read: pysam.AlignedSegment) -> list[tuple[str, int, str]]:
    if not read.has_tag("SA"):
        return []
    sa_raw = read.get_tag("SA")
    out: list[tuple[str, int, str]] = []
    for entry in str(sa_raw).split(";"):
        if not entry:
            continue
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
    for start1, end1 in regions_1based:
        start0 = max(0, start1 - 1)
        end0 = max(0, end1)
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
    if read.cigartuples is None:
        return False

    min_len = max(MIN_CIGAR_BP, int(MIN_CIGAR_FRACTION * svlen))

    ref = read.reference_start
    for op, length in read.cigartuples:
        if op in (0, 7, 8):
            ref += length
        elif op in (2, 3):
            if op == 2 and length >= min_len:
                del_start = ref
                del_end = ref + length
                if abs(del_start - pos0) <= bp_window and abs(del_end - end_excl0) <= bp_window:
                    return True
            ref += length
        elif op in (1, 4, 5, 6):
            continue

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
    if read.cigartuples is None:
        return False

    min_len = _ins_support_min_len(svlen)

    ref = read.reference_start
    for op, length in read.cigartuples:
        if op in (0, 7, 8):
            ref += length
        elif op == 1:
            if length >= min_len and abs(ref - pos0) <= bp_window:
                return True
        elif op in (4, 5):
            if length >= min_len and abs(ref - pos0) <= bp_window:
                return True
        elif op in (2, 3):
            ref += length
    return False


def _supports_bnd(
    read: pysam.AlignedSegment,
    *,
    pos0: int,
    chr2: str,
    pos2_1based: int,
    bp_window: int,
) -> bool:
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
    if svtype == "DEL":
        return _supports_del(read, pos0=pos0, end_excl0=end_excl0, svlen=svlen, bp_window=bp_window)
    if svtype == "INS":
        return _supports_ins(read, pos0=pos0, svlen=svlen, bp_window=bp_window)
    if svtype == "BND" and chr2 and pos2:
        return _supports_bnd(
            read, pos0=pos0, chr2=str(chr2), pos2_1based=int(pos2), bp_window=bp_window
        )
    if svtype == "INV":
        return _supports_inv(read, pos0=pos0, end0=(end_excl0 - 1), bp_window=bp_window)

    return _supports_ins(read, pos0=pos0, svlen=svlen, bp_window=bp_window) or _supports_del(
        read, pos0=pos0, end_excl0=end_excl0, svlen=svlen, bp_window=bp_window
    )


def _count_hp_sv_support(  # noqa: C901
    bam: pysam.AlignmentFile,
    chrom: str,
    rec: Variant,
    *,
    opts: WorkerOpts,
) -> dict[str, Any]:
    pos1 = int(rec.POS)
    sv_end = int(rec.end) if getattr(rec, "end", None) is not None else pos1
    alt = ",".join(rec.ALT) if rec.ALT else "<N>"
    svtype = str(rec.INFO.get("SVTYPE", "NA"))

    in_gt: str | None = None
    if rec.genotypes:
        a1, a2 = rec.genotypes[0][0], rec.genotypes[0][1]
        in_gt = f"{a1}/{a2}"

    svlen = _svlen_from_record(rec, pos1, sv_end)
    fetch_w, bp_tol = _compute_windows(rec, pos1, sv_end, svlen, opts=opts)

    rset: set[str] = set()
    if opts.support_mode in {"hybrid", "rnames"}:
        rset = _parse_rnames(rec.INFO.get("RNAMES"))

    pos0 = pos1 - 1
    end_excl0 = sv_end

    chr2 = None
    pos2 = None
    if svtype == "BND":
        chr2, pos2 = _parse_bnd_partner(alt, rec)

    if rset:
        state: dict[str, dict[str, Any]] = {}

        left0 = min(pos0, max(0, end_excl0 - 1))
        right0 = max(pos0, max(0, end_excl0 - 1))
        start0 = max(0, left0 - fetch_w)
        stop0 = right0 + 1 + fetch_w

        for read in bam.fetch(chrom, start0, stop0):
            if read.is_unmapped or read.is_secondary:
                continue
            qn = read.query_name
            if not qn or qn not in rset:
                continue
            st = state.setdefault(qn, {"hp": None})
            if st["hp"] is None and read.has_tag("HP"):
                st["hp"] = read.get_tag("HP")

        hp1 = hp2 = nohp = 0
        for _qn, st in state.items():
            hp = st.get("hp")
            if hp == 1:
                hp1 += 1
            elif hp == 2:
                hp2 += 1
            else:
                nohp += 1

        tagged_total = hp1 + hp2
        support_total = hp1 + hp2 + nohp

        return {
            "hp1": hp1,
            "hp2": hp2,
            "nohp": nohp,
            "tagged_total": tagged_total,
            "support_total": support_total,
            "sv_end": sv_end,
            "alt": alt,
            "svtype": svtype,
            "svlen": svlen,
            "mode": "RNAMES",
            "fetch_w": fetch_w,
            "bp_tol": bp_tol,
            "rnames_total": len(rset),
            "rnames_found": len(state),
            "in_gt": in_gt,
        }

    if opts.support_mode == "rnames":
        return {
            "hp1": 0,
            "hp2": 0,
            "nohp": 0,
            "tagged_total": 0,
            "support_total": 0,
            "sv_end": sv_end,
            "alt": alt,
            "svtype": svtype,
            "svlen": svlen,
            "mode": "RNAMES",
            "fetch_w": fetch_w,
            "bp_tol": bp_tol,
            "rnames_total": 0,
            "rnames_found": 0,
            "in_gt": in_gt,
        }

    regions: list[tuple[int, int]] = [(max(1, pos1 - fetch_w), pos1 + fetch_w)]
    if svtype in {"DEL", "INV"} and sv_end != pos1:
        regions.append((max(1, sv_end - fetch_w), sv_end + fetch_w))

    state: dict[str, dict[str, Any]] = {}

    for read in _iter_candidate_reads(bam, chrom, regions):
        qn = read.query_name
        if qn is None:
            continue
        st = state.setdefault(qn, {"hp": None, "support": False})

        if st["hp"] is None and read.has_tag("HP"):
            st["hp"] = read.get_tag("HP")

        if st["support"] is True:
            continue

        if _read_supports_variant(
            read,
            svtype,
            pos0=pos0,
            end_excl0=end_excl0,
            svlen=svlen,
            bp_window=bp_tol,
            chr2=chr2,
            pos2=pos2,
        ):
            st["support"] = True

    hp1 = hp2 = nohp = 0
    for st in state.values():
        if not st["support"]:
            continue
        hp = st["hp"]
        if hp == 1:
            hp1 += 1
        elif hp == 2:
            hp2 += 1
        else:
            nohp += 1

    tagged_total = hp1 + hp2
    support_total = hp1 + hp2 + nohp

    return {
        "hp1": hp1,
        "hp2": hp2,
        "nohp": nohp,
        "tagged_total": tagged_total,
        "support_total": support_total,
        "sv_end": sv_end,
        "alt": alt,
        "svtype": svtype,
        "svlen": svlen,
        "mode": "HEURISTIC",
        "fetch_w": fetch_w,
        "bp_tol": bp_tol,
        "rnames_total": 0,
        "rnames_found": 0,
        "in_gt": in_gt,
    }


def _phase_chrom_worker(
    chrom: str,
    vcf_path: Path,
    bam_path: Path,
    opts: WorkerOpts,
) -> pd.DataFrame:
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    rdr = Reader(str(vcf_path))

    rows: list[dict[str, object]] = []

    use_region_iter = _has_tabix_index(vcf_path)
    records_iter = rdr(chrom) if use_region_iter else iter(rdr)

    for rec in records_iter:
        if rec.CHROM != chrom:
            continue

        sup = _count_hp_sv_support(bam, chrom, rec, opts=opts)
        hp1 = int(sup["hp1"])
        hp2 = int(sup["hp2"])
        nohp = int(sup["nohp"])
        tagged_total = int(sup["tagged_total"])
        support_total = int(sup["support_total"])

        # classify returns: (gt, gq, reason, delta)
        gt, gq, reason, delta = classify_haplotype_v211(
            n1=hp1,
            n2=hp2,
            min_support=opts.min_support,
            min_tagged_support=opts.min_tagged_support,
            major_delta=opts.major_delta,
            equal_delta=opts.equal_delta,
            support_total=support_total,
            tie_to_hom_alt=opts.tie_to_hom_alt,
        )

        tag_frac = (tagged_total / support_total) if support_total else 0.0

        row = {
            "chrom": chrom,
            "pos": int(rec.POS),
            "id": rec.ID,
            "svtype": str(sup.get("svtype", "NA")),
            "svlen": int(sup.get("svlen") or 0),
            "end": int(sup.get("sv_end") or rec.POS),
            # carry ALT + input genotype for downstream/debugging
            "alt": sup.get("alt"),
            "in_gt": sup.get("in_gt"),
            # evidence counts
            "hp1": hp1,
            "hp2": hp2,
            "nohp": nohp,
            "tagged_total": tagged_total,
            "support_total": support_total,
            # aliases expected by VCF writer (n1/n2)
            "n1": hp1,
            "n2": hp2,
            # classification outputs (previously missing -> caused GT=./.)
            "gt": gt,
            "gq": int(gq),
            "reason": reason,
            "delta": float(delta),
            "tag_frac": float(tag_frac),
            # provenance
            "mode": sup.get("mode"),
            "fetch_w": sup.get("fetch_w"),
            "bp_tol": sup.get("bp_tol"),
            "rnames_total": sup.get("rnames_total"),
            "rnames_found": sup.get("rnames_found"),
        }
        rows.append(row)

    return pd.DataFrame(rows)
