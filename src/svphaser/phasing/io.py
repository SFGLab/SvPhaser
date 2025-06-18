"""
svphaser.phasing.io
===================
High-level “engine” – orchestrates per-chromosome workers, merges their
results, applies global depth filters, then writes the CSV + VCF.

Only *simple* (picklable) objects are sent to worker processes; each worker
re-opens the VCF/BAM for its chromosome.  That avoids the earlier
`TypeError: no default __reduce__ due to non-trivial __cinit__` from trying
to ship `cyvcf2.Variant` instances across process boundaries.
"""
from __future__ import annotations

import logging
import multiprocessing as mp
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
from cyvcf2 import Reader  # only used in the parent for header/seqnames

from ._workers import _phase_chrom_worker
from .types import GQBin, SVKey, WorkerOpts

__all__ = ["phase_vcf"]

logger = logging.getLogger("svphaser.io")


# ────────────────────────────────────────────────────────────────────────
#  Public entry-point
# ────────────────────────────────────────────────────────────────────────
def phase_vcf(
    sv_vcf: Path,
    bam: Path,
    *,
    out_dir: Path,
    min_support: int,
    major_delta: float,
    equal_delta: float,
    gq_bins: str,
    threads: int | None,
) -> None:
    """
    Phase *sv_vcf* against *bam* and write ``*_phased.{vcf,csv}`` to *out_dir*.
    All numeric knobs are documented in `svphaser --help`.
    """

    # 1 ─ Parse --gq-bins → List[Tuple[int,label]]
    bins: List[GQBin] = []
    if gq_bins.strip():
        for part in gq_bins.split(","):
            thr_s, lbl = part.strip().split(":")
            bins.append((int(thr_s), lbl))
        # make sure highest threshold comes first so “first match wins”
        bins.sort(key=lambda x: x[0], reverse=True)

    # 2 ─ Build immutable options holder for workers
    opts = WorkerOpts(
        min_support=min_support,
        major_delta=major_delta,
        equal_delta=equal_delta,
        gq_bins=bins,
    )

    # 3 ─ Discover chromosomes (cheap – no variants parsed yet)
    rdr = Reader(str(sv_vcf))
    chroms: Tuple[str, ...] = tuple(rdr.seqnames)
    rdr.close()

    # 4 ─ Launch one worker per chromosome (or ≤threads)
    worker_args: List[Tuple[str, Path, Path, WorkerOpts]] = [
        (chrom, sv_vcf, bam, opts) for chrom in chroms
    ]

    threads = threads or mp.cpu_count() or 1
    logger.info("SvPhaser ▶ workers: %d", threads)

    dataframes: List[pd.DataFrame] = []
    ctx = mp.get_context("fork")
    with ctx.Pool(processes=threads) as pool:
        # starmap expands each tuple into positional args
        for df in pool.starmap(_phase_chrom_worker, worker_args, chunksize=1):
            dataframes.append(df)
            chrom = df.iloc[0]["chrom"] if not df.empty else "?"
            logger.info("chr %-6s ✔ phased %5d SVs", chrom, len(df))

    # 5 ─ Merge & apply *global* depth filter
    merged = pd.concat(dataframes, ignore_index=True)
    pre = len(merged)
    keep = (merged["n1"] >= min_support) | (merged["n2"] >= min_support)
    merged = merged[keep].reset_index(drop=True)
    if (dropped := pre - len(merged)):
        logger.info("Depth filter removed %d SVs", dropped)

    # 6 ─ Write CSV
    stem = sv_vcf.name.removesuffix(".vcf.gz").removesuffix(".vcf")
    out_csv = out_dir / f"{stem}_phased.csv"
    out_vcf = out_dir / f"{stem}_phased.vcf"
    merged.to_csv(out_csv, index=False)
    logger.info("CSV → %s  (%d SVs)", out_csv, len(merged))

    # 7 ─ Write VCF
    _write_phased_vcf(out_vcf, sv_vcf, merged)
    logger.info("VCF → %s", out_vcf)


# ────────────────────────────────────────────────────────────────────────
#  Helper – copy input VCF, injecting SvPhaser INFO fields
# ────────────────────────────────────────────────────────────────────────
from cyvcf2 import Reader, Writer, Variant   # ← Writer added

def _write_phased_vcf(out_vcf: Path, in_vcf: Path, df: pd.DataFrame) -> None:
    """Copy *in_vcf* → *out_vcf* while adding HP_GT / HP_GQ / HP_GQBIN INFO."""

    key2call: Dict[SVKey, Tuple[str, int, str | None]] = {
        (str(r.chrom), int(r.pos), str(r.id)): ( #type: ignore[assignment]
            str(r.gt),
            int(r.gq), #type: ignore[assignment]
            getattr(r, "gq_label", None),   # safe even if column absent
        )
        for r in df.itertuples()
    }

    rdr = Reader(str(in_vcf))

    def _ensure(iid: str, desc: str, typ: str) -> None:
        if f"ID={iid}," not in rdr.raw_header:
            rdr.add_info_to_header(dict(ID=iid, Description=desc, Type=typ, Number=1))

    _ensure("HP_GT", "Haplotype genotype (0|1,1|0,1|1,./.) from SvPhaser", "String")
    _ensure("HP_GQ", "Phred-scaled genotype quality (0–99) from SvPhaser", "Integer")
    if "gq_label" in df.columns:
        _ensure("HP_GQBIN", "Confidence bin label from SvPhaser", "String")

    # Portable header cloning (works on all cyvcf2 versions)
    wtr = Writer(str(out_vcf), rdr)

    for rec in rdr:                       
        key = (rec.CHROM, rec.POS, rec.ID or ".")
        if key in key2call:
            gt, gq, lbl = key2call[key]
            rec.INFO["HP_GT"] = gt
            rec.INFO["HP_GQ"] = gq
            if lbl is not None:
                rec.INFO["HP_GQBIN"] = lbl
        wtr.write_record(rec)

    wtr.close()
    rdr.close()
