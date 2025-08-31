"""
svphaser.phasing._workers
=========================
Worker-process code.  Each worker:

1. Opens the (possibly un-indexed) SV VCF.
2. Scans only the records for *its* chromosome.
3. Counts HP-tagged reads in the long-read BAM/CRAM.
4. Classifies the haplotype + GQ, adds optional GQ-bin label.
5. Returns a DataFrame to the parent.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pysam
from cyvcf2 import Reader, Variant  # type: ignore

from .algorithms import classify_haplotype
from .types import WorkerOpts

__all__ = ["_phase_chrom_worker"]


def _count_hp_reads(
    bam: pysam.AlignmentFile,
    chrom: str,
    start: int,
    end: int,
) -> tuple[int, int]:
    hp1 = hp2 = 0
    for read in bam.fetch(chrom, max(0, start - 1), end + 1):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if not read.has_tag("HP"):
            continue
        tag = read.get_tag("HP")
        if tag == 1:
            hp1 += 1
        elif tag == 2:
            hp2 += 1
    return hp1, hp2


def _has_tabix_index(vcf_path: Path) -> bool:
    """Return True if <file>.tbi or <file>.csi exists (supports .vcf.gz.{tbi,csi})."""
    return (
        vcf_path.with_suffix(vcf_path.suffix + ".tbi").exists()
        or vcf_path.with_suffix(vcf_path.suffix + ".csi").exists()
    )


def _phase_chrom_worker(
    chrom: str,
    vcf_path: Path,
    bam_path: Path,
    opts: WorkerOpts,
) -> pd.DataFrame:
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    rdr = Reader(str(vcf_path))

    rows: list[dict[str, object]] = []

    # Try fast random access first, fall back to linear scan if that fails
    use_region_iter = _has_tabix_index(vcf_path)
    records_iter = (
        rdr(f"{chrom}") if use_region_iter else (rec for rec in rdr if rec.CHROM == chrom)
    )

    for rec in records_iter:  # type: ignore[arg-type]
        assert isinstance(rec, Variant)
        sv_end = rec.end if getattr(rec, "end", None) is not None else rec.POS
        n1, n2 = _count_hp_reads(bam, chrom, rec.POS, sv_end)

        gt, gq = classify_haplotype(
            n1,
            n2,
            min_support=opts.min_support,
            major_delta=opts.major_delta,
            equal_delta=opts.equal_delta,
        )

        row = dict(
            chrom=chrom,
            pos=rec.POS,  # cyvcf2 POS is already 1-based
            id=rec.ID or ".",
            svtype=rec.INFO.get("SVTYPE", "NA"),
            n1=n1,
            n2=n2,
            gt=gt,
            gq=gq,
        )

        if opts.gq_bins:
            for thr, label in opts.gq_bins:
                if gq >= thr:
                    row["gq_label"] = label
                    break

        rows.append(row)  # type: ignore[assignment]

    rdr.close()
    bam.close()
    return pd.DataFrame(rows)
