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
from cyvcf2 import Reader

from ._workers import _phase_chrom_worker
from .types import GQBin, SVKey, WorkerOpts

__all__ = ["phase_vcf"]

logger = logging.getLogger("svphaser.io")


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
        for df in pool.starmap(_phase_chrom_worker, worker_args, chunksize=1):
            dataframes.append(df)
            chrom = df.iloc[0]["chrom"] if not df.empty else "?"
            logger.info("chr %-6s ✔ phased %5d SVs", chrom, len(df))

    # 5 ─ Merge & apply *global* depth filter
    merged = pd.concat(dataframes, ignore_index=True)
    pre = len(merged)
    keep = ~((merged["n1"] < min_support) & (merged["n2"] < min_support))

    # --- Save dropped SVs here ---
    dropped_svs = merged[~keep]
    stem = sv_vcf.name.removesuffix(".vcf.gz").removesuffix(".vcf")
    dropped_csv = out_dir / f"{stem}_dropped_svs.csv"
    dropped_svs.to_csv(dropped_csv, index=False)
    logger.info("Dropped SVs saved to %s (%d SVs)", dropped_csv, len(dropped_svs))
    # ----------------------------

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


def _write_phased_vcf(out_vcf: Path, in_vcf: Path, df: pd.DataFrame) -> None:
    """
    Write a phased VCF using the merged DataFrame, reconstructing INFO and all columns,
    tab-delimited, compliant with VCF tools.
    """
    input_vcf = Reader(str(in_vcf))

    # Build info lookup: key = (CHROM, POS, ID) for fast ref/alt/qual/FILTER/INFO lookup
    vcf_info_lookup = {}
    for rec in input_vcf:
        key = (rec.CHROM, rec.POS, rec.ID or ".")
        info_dict = {}
        for k in rec.INFO:
            info_key = k[0] if isinstance(k, tuple) else k
            v = rec.INFO.get(info_key)
            if v is not None:
                info_dict[info_key] = v
        vcf_info_lookup[key] = {
            'REF': rec.REF,
            'ALT': rec.ALT[0] if rec.ALT else "<N>",
            'QUAL': rec.QUAL if rec.QUAL is not None else '.',
            'FILTER': rec.FILTER if rec.FILTER else 'PASS',
            'INFO': info_dict,
        }

    # Re-open to get header from the start
    input_vcf = Reader(str(in_vcf))

    with open(out_vcf, "w") as out:
        # Write all header lines (those starting with "##")
        for line in input_vcf.raw_header.strip().splitlines():
            if line.startswith("##"):
                out.write(line.rstrip() + "\n")

        # Write a SINGLE column header (the only one!)
        sample_name = input_vcf.samples[0] if input_vcf.samples else "SAMPLE"
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_name + "\n")

        # Now write all records
        for row in df.itertuples(index=False):
            chrom = getattr(row, "chrom", ".")
            pos = int(getattr(row, "pos", 0))
            id = getattr(row, "id", ".")
            gt = getattr(row, "gt", "./.")
            gq = getattr(row, "gq", "0")
            svtype = getattr(row, "svtype", None)
            gq_label = getattr(row, "gq_label", None)

            key = (str(chrom), pos, str(id))
            info = vcf_info_lookup.get(key)
            if info is None:
                key_alt = (str(chrom), pos - 1, str(id))
                info = vcf_info_lookup.get(key_alt)
            if info is None:
                logger.warning(f"Could not find VCF info for {chrom}:{pos} {id}")
                continue

            ref = info['REF']
            alt = info['ALT']
            qual = info['QUAL']
            filt = info['FILTER']

            # Compose INFO: always put SVTYPE first
            orig_info = [f"{k}={v}" for k, v in info['INFO'].items() if k != "SVTYPE"]
            if svtype:
                orig_info.insert(0, f"SVTYPE={svtype}")
            if gq_label and pd.notnull(gq_label):
                orig_info.append(f"GQBIN={gq_label}")
            info_str = ";".join(orig_info) if orig_info else "."

            format_str = "GT:GQ"
            sample_str = f"{gt}:{gq}"

            fields = [
                str(chrom), str(pos), str(id), str(ref), str(alt),
                str(qual), str(filt), info_str, format_str, sample_str
            ]
            out.write('\t'.join(fields) + '\n')