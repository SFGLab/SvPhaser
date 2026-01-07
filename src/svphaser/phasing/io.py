"""svphaser.phasing.io
===================
High-level “engine” – orchestrates per-chromosome workers, merges results,
applies the global support filter, then writes CSV + VCF.

Step B update (biological correctness):
- `min_support` is interpreted as a *total ALT-support* threshold (n1+n2).
  The filter drops an SV only if (n1+n2) < min_support.

Step A fixes retained:
- collision-resistant VCF record matching using (CHROM, POS, ID, END, ALT)
- correct INFO composition (no duplicated keys; proper FLAG handling)
- typing fixes to satisfy Pylance/Mypy
"""

from __future__ import annotations

import logging
import math
import multiprocessing as mp
from pathlib import Path
from typing import Any, TypedDict

import pandas as pd
from cyvcf2 import Reader

from ._workers import _phase_chrom_worker
from .types import GQBin, SVKey, SVKeyLegacy, WorkerOpts

__all__ = ["phase_vcf"]

logger = logging.getLogger(__name__)


class VcfRec(TypedDict):
    REF: str
    ALT: str
    QUAL: object
    FILTER: str
    INFO: dict[str, Any]


def _is_missing_scalar(x: Any) -> bool:
    """True for None / NaN / empty string."""
    if x is None:
        return True
    if isinstance(x, float) and math.isnan(x):
        return True
    if isinstance(x, str) and x.strip() == "":
        return True
    return False


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
                    f"Invalid gq-bin specifier: '{thr_lbl}'. Use '30:High,10:Moderate'."
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

    # 5 ─ Merge & apply *global* support filter (Step B: total ALT-support)
    if dataframes:
        merged = pd.concat(dataframes, ignore_index=True)
    else:
        merged = pd.DataFrame(
            columns=[
                "chrom",
                "pos",
                "end",
                "id",
                "alt",
                "svtype",
                "n1",
                "n2",
                "gt",
                "gq",
                "gq_label",
            ]
        )

    pre = len(merged)
    total_support = merged["n1"].astype(int) + merged["n2"].astype(int)
    keep = total_support >= int(min_support)

    stem = sv_vcf.name.removesuffix(".vcf.gz").removesuffix(".vcf")

    # Save dropped SVs for transparency
    dropped_csv = out_dir / f"{stem}_dropped_svs.csv"
    merged.loc[~keep].to_csv(dropped_csv, index=False)
    logger.info("Dropped SVs → %s (%d SVs)", dropped_csv, int((~keep).sum()))

    kept = merged.loc[keep].reset_index(drop=True)
    if dropped := pre - len(kept):
        logger.info("Support filter removed %d SVs", dropped)

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
) -> tuple[dict[SVKey, VcfRec], dict[SVKeyLegacy, list[SVKey]], list[str], str]:
    """Scan input VCF once.

    Returns:
      - full_lookup: maps (CHROM, POS, ID, END, ALT) -> record components
      - legacy_index: maps (CHROM, POS, ID) -> list of full keys (fallback)
      - raw_header_lines
      - sample_name
    """
    rdr = Reader(str(in_vcf))
    raw_header_lines = rdr.raw_header.strip().splitlines()
    sample_name = rdr.samples[0] if rdr.samples else "SAMPLE"

    full_lookup: dict[SVKey, VcfRec] = {}
    legacy_index: dict[SVKeyLegacy, list[SVKey]] = {}

    for rec in rdr:
        chrom = rec.CHROM
        pos = int(rec.POS)
        vid = rec.ID or "."
        end = int(rec.end) if getattr(rec, "end", None) is not None else int(pos)
        alt = ",".join(rec.ALT) if rec.ALT else "<N>"

        info_dict: dict[str, Any] = {}
        for k in rec.INFO:
            info_key = k[0] if isinstance(k, tuple) else k
            v = rec.INFO.get(info_key)
            if v is not None:
                info_dict[info_key] = v

        fkey: SVKey = (chrom, pos, vid, end, alt)
        lkey: SVKeyLegacy = (chrom, pos, vid)

        full_lookup[fkey] = {
            "REF": rec.REF,
            "ALT": alt,
            "QUAL": rec.QUAL if rec.QUAL is not None else ".",
            "FILTER": rec.FILTER if rec.FILTER else "PASS",
            "INFO": info_dict,
        }
        legacy_index.setdefault(lkey, []).append(fkey)

    rdr.close()
    return full_lookup, legacy_index, raw_header_lines, sample_name


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


def _compose_info_str(orig_info: dict[str, Any], svtype: Any, gq_label: Any) -> str:
    """Compose INFO with SVTYPE first, proper FLAG handling, then GQBIN."""
    items: list[str] = []

    if svtype:
        items.append(f"SVTYPE={svtype}")

    for k, v in orig_info.items():
        if k == "SVTYPE":
            continue
        if v is None:
            continue
        # cyvcf2 represents INFO flags as boolean True
        if v is True:
            items.append(str(k))
        else:
            items.append(f"{k}={v}")

    if not _is_missing_scalar(gq_label):
        items.append(f"GQBIN={gq_label}")

    return ";".join(items) if items else "."


def _select_info_record(
    full_lookup: dict[SVKey, VcfRec],
    legacy_index: dict[SVKeyLegacy, list[SVKey]],
    *,
    chrom: str,
    pos: int,
    vid: str,
    end: int | None,
    alt: str | None,
) -> VcfRec | None:
    """Pick the best matching input VCF record for this phased row."""
    if end is not None and alt is not None:
        hit = full_lookup.get((chrom, pos, vid, int(end), str(alt)))
        if hit is not None:
            return hit

    candidates = legacy_index.get((chrom, pos, vid), [])
    if not candidates:
        return None

    if len(candidates) == 1:
        return full_lookup[candidates[0]]

    if end is not None:
        end_matches = [k for k in candidates if k[3] == int(end)]
        if len(end_matches) == 1:
            return full_lookup[end_matches[0]]

    # Still ambiguous: refuse to guess
    return None


def _write_phased_vcf(
    out_vcf: Path,
    in_vcf: Path,
    df: pd.DataFrame,
    *,
    gqbin_in_header: bool,
) -> None:
    """Write a phased VCF: tab-delimited, compliant, with ensured GT/GQ (and GQBIN if used)."""
    full_lookup, legacy_index, raw_header_lines, sample_name = _vcf_info_lookup(in_vcf)

    with open(out_vcf, "w", newline="") as out:
        _write_headers(out, raw_header_lines, sample_name, gqbin_in_header=gqbin_in_header)

        for row in df.itertuples(index=False):
            chrom = str(getattr(row, "chrom", "."))
            pos = int(getattr(row, "pos", 0))
            vid = str(getattr(row, "id", "."))

            end = getattr(row, "end", None)
            alt = getattr(row, "alt", None)

            gt = str(getattr(row, "gt", "./."))
            gq = str(getattr(row, "gq", "0"))
            svtype = getattr(row, "svtype", None)
            gq_label = getattr(row, "gq_label", None)

            info = _select_info_record(
                full_lookup,
                legacy_index,
                chrom=chrom,
                pos=pos,
                vid=vid,
                end=int(end) if end is not None else None,
                alt=str(alt) if alt is not None else None,
            )
            if info is None:
                logger.warning("Could not uniquely match VCF record for %s:%s %s", chrom, pos, vid)
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
