"""
svphaser.phasing.io
===================
High-level “engine” – orchestrates per-chromosome workers, merges results,
applies the global depth filter, then writes CSV + VCF.

Workers receive only simple (pickle-safe) arguments; each worker opens its
own BAM/VCF to avoid sharing handles between processes.
"""

from __future__ import annotations

import logging
import multiprocessing as mp
from pathlib import Path

import pandas as pd
from cyvcf2 import Reader

from ._workers import _phase_chrom_worker
from .types import GQBin, WorkerOpts

__all__ = ["phase_vcf"]

logger = logging.getLogger(__name__)


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
    """Phase *sv_vcf* against *bam* and write outputs to *out_dir*.

    Files:
      - *_phased.vcf
      - *_phased.csv
      - *_dropped_svs.csv
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1 ─ Parse --gq-bins → list[(int,label)]
    bins: list[GQBin] = []
    if gq_bins.strip():
        for part in gq_bins.split(","):
            thr_lbl = part.strip()
            if not thr_lbl:
                continue
            try:
                thr_s, lbl = thr_lbl.split(":")
            except ValueError as err:
                raise ValueError(
                    f"Invalid gq-bin specifier: '{thr_lbl}'. " "Use '30:High,10:Moderate'."
                ) from err
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
    chroms: tuple[str, ...] = tuple(rdr.seqnames)
    rdr.close()

    # 4 ─ Launch one worker per chromosome (or ≤threads)
    worker_args: list[tuple[str, Path, Path, WorkerOpts]] = [(c, sv_vcf, bam, opts) for c in chroms]

    threads = threads or mp.cpu_count() or 1
    logger.info("SvPhaser ▶ workers: %d", threads)

    dataframes: list[pd.DataFrame] = []

    # Use 'fork' when available (fast on Linux); fall back to 'spawn' elsewhere.
    try:
        ctx = mp.get_context("fork")
    except ValueError:
        ctx = mp.get_context("spawn")

    if threads == 1:
        # Serial path is handy for debugging
        for args in worker_args:
            df = _phase_chrom_worker(*args)
            dataframes.append(df)
            chrom = df.iloc[0]["chrom"] if not df.empty else "?"
            logger.info("chr %-6s ✔ phased %5d SVs", chrom, len(df))
    else:
        with ctx.Pool(processes=threads) as pool:
            for df in pool.starmap(_phase_chrom_worker, worker_args, chunksize=1):
                dataframes.append(df)
                chrom = df.iloc[0]["chrom"] if not df.empty else "?"
                logger.info("chr %-6s ✔ phased %5d SVs", chrom, len(df))

    # 5 ─ Merge & apply *global* depth filter
    if dataframes:
        merged = pd.concat(dataframes, ignore_index=True)
    else:
        merged = pd.DataFrame(
            columns=["chrom", "pos", "id", "svtype", "n1", "n2", "gt", "gq", "gq_label"]
        )

    pre = len(merged)
    keep = ~((merged["n1"] < min_support) & (merged["n2"] < min_support))

    stem = sv_vcf.name.removesuffix(".vcf.gz").removesuffix(".vcf")

    # Save dropped SVs for transparency
    dropped_csv = out_dir / f"{stem}_dropped_svs.csv"
    merged.loc[~keep].to_csv(dropped_csv, index=False)
    logger.info("Dropped SVs → %s (%d SVs)", dropped_csv, int((~keep).sum()))

    kept = merged.loc[keep].reset_index(drop=True)
    if dropped := pre - len(kept):
        logger.info("Depth filter removed %d SVs", dropped)

    # 6 ─ Write CSV
    out_csv = out_dir / f"{stem}_phased.csv"
    kept.to_csv(out_csv, index=False)
    logger.info("CSV → %s  (%d SVs)", out_csv, len(kept))

    # 7 ─ Write VCF
    out_vcf = out_dir / f"{stem}_phased.vcf"
    _write_phased_vcf(out_vcf, sv_vcf, kept, gqbin_in_header=bool(bins))
    logger.info("VCF → %s", out_vcf)


# ──────────────────────────────────────────────────────────────────────
#  Small helpers to keep complexity down
# ──────────────────────────────────────────────────────────────────────
def _vcf_info_lookup(
    in_vcf: Path,
) -> tuple[dict[tuple[str, int, str], dict[str, object]], list[str], str]:
    """Scan input VCF once: return (lookup, raw_header_lines, sample_name)."""
    rdr = Reader(str(in_vcf))
    raw_header_lines = rdr.raw_header.strip().splitlines()
    sample_name = rdr.samples[0] if rdr.samples else "SAMPLE"

    lookup: dict[tuple[str, int, str], dict[str, object]] = {}
    for rec in rdr:
        key = (rec.CHROM, rec.POS, rec.ID or ".")
        info_dict: dict[str, object] = {}
        for k in rec.INFO:
            info_key = k[0] if isinstance(k, tuple) else k
            v = rec.INFO.get(info_key)
            if v is not None:
                info_dict[info_key] = v
        lookup[key] = {
            "REF": rec.REF,
            "ALT": ",".join(rec.ALT) if rec.ALT else "<N>",
            "QUAL": rec.QUAL if rec.QUAL is not None else ".",
            "FILTER": rec.FILTER if rec.FILTER else "PASS",
            "INFO": info_dict,
        }
    rdr.close()
    return lookup, raw_header_lines, sample_name


def _write_headers(
    out,
    raw_header_lines: list[str],
    sample_name: str,
    *,
    gqbin_in_header: bool,
) -> None:
    """Write preserved meta headers + ensure GT/GQ/GQBIN, then the column header."""
    have_gt = any("##FORMAT=<ID=GT" in ln for ln in raw_header_lines)
    have_gq = any("##FORMAT=<ID=GQ" in ln for ln in raw_header_lines)
    have_gqbin = any("##INFO=<ID=GQBIN" in ln for ln in raw_header_lines)

    for line in raw_header_lines:
        if line.startswith("##"):
            out.write(line.rstrip() + "\n")

    if not have_gt:
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotype">\n')
    if not have_gq:
        out.write(
            "##FORMAT=<ID=GQ,Number=1,Type=Integer," 'Description="Genotype Quality (Phred)">\n'
        )
    if gqbin_in_header and not have_gqbin:
        out.write(
            "##INFO=<ID=GQBIN,Number=1,Type=String," 'Description="GQ bin label from SvPhaser">\n'
        )

    out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_name + "\n")


def _compose_info_str(orig_info: dict[str, object], svtype: object, gq_label: object) -> str:
    """Compose the INFO string with SVTYPE first, original keys (no duplicate), then GQBIN."""
    items: list[str] = []
    for k, v in orig_info.items():
        if k == "SVTYPE":
            continue
        items.append(f"{k}={v}")
        # treat boolean True as a FLAG (bare key). Keep everything else as k=v.
        if v is True:
            items.append(k)
        else:
            items.append(f"{k}={v}")
    if svtype:
        items.insert(0, f"SVTYPE={svtype}")
    if gq_label is not None and pd.notnull(gq_label):
        items.append(f"GQBIN={gq_label}")
    return ";".join(items) if items else "."


def _write_phased_vcf(
    out_vcf: Path,
    in_vcf: Path,
    df: pd.DataFrame,
    *,
    gqbin_in_header: bool,
) -> None:
    """Write a phased VCF: tab-delimited, compliant, with ensured GT/GQ (and GQBIN if used)."""
    lookup, raw_header_lines, sample_name = _vcf_info_lookup(in_vcf)

    with open(out_vcf, "w", newline="") as out:
        _write_headers(out, raw_header_lines, sample_name, gqbin_in_header=gqbin_in_header)

        for row in df.itertuples(index=False):
            chrom = str(getattr(row, "chrom", "."))
            pos = int(getattr(row, "pos", 0))
            vid = str(getattr(row, "id", "."))
            gt = str(getattr(row, "gt", "./."))
            gq = str(getattr(row, "gq", "0"))
            svtype = getattr(row, "svtype", None)
            gq_label = getattr(row, "gq_label", None)

            info = lookup.get((chrom, pos, vid))
            if info is None:
                logger.warning("Could not find VCF info for %s:%s %s", chrom, pos, vid)
                continue

            info_str = _compose_info_str(info["INFO"], svtype, gq_label)

            fields = [
                chrom,
                str(pos),
                vid,
                str(info["REF"]),
                str(info["ALT"]),
                str(info["QUAL"]),
                str(info["FILTER"]),
                info_str,
                "GT:GQ",
                f"{gt}:{gq}",
            ]
            out.write("\t".join(fields) + "\n")
