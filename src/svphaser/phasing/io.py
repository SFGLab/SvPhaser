"""svphaser.phasing.io
===================

High-level “engine” – orchestrates per-chromosome workers, merges results,
applies the global support filter, then writes CSV + VCF.

Hardening patch (v2.1.3+):
- FAIL LOUD if worker output is missing required phasing columns (gt/gq).
- Backfill expected aliases (n1/n2, tag_frac) from hp1/hp2 + totals when possible.
- Compute gq_label (GQBIN) from --gq-bins if missing.
- Warn if output looks suspicious (e.g., 99% './.' despite high tagged support).
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


def _gq_label_from_bins(gq: Any, bins: list[GQBin]) -> str | None:
    """Return label for first threshold satisfied (bins sorted desc)."""
    if not bins:
        return None
    try:
        gq_i = int(gq)
    except Exception:
        return None
    for thr, lbl in bins:
        if gq_i >= int(thr):
            return str(lbl)
    return None


# ──────────────────────────────────────────────────────────────────────
#  Hardening: normalize worker output before writing
# ──────────────────────────────────────────────────────────────────────


def _backfill_n1_n2(out: pd.DataFrame) -> None:
    # n1/n2 are what the VCF writer uses (SVP_HP1/HP2).
    if "n1" not in out.columns and "hp1" in out.columns:
        out["n1"] = out["hp1"]
    if "n2" not in out.columns and "hp2" in out.columns:
        out["n2"] = out["hp2"]


def _backfill_tag_frac(out: pd.DataFrame) -> None:
    # tag_frac is used in SVP_TAGFRAC; derive if totals exist.
    if "tag_frac" in out.columns:
        return
    if "tagged_total" not in out.columns or "support_total" not in out.columns:
        return

    denom = pd.to_numeric(out["support_total"], errors="coerce").astype(float)
    denom = denom.where(denom != 0.0, other=float("nan"))
    numer = pd.to_numeric(out["tagged_total"], errors="coerce").astype(float)
    out["tag_frac"] = numer / denom


def _backfill_gq_label(out: pd.DataFrame, *, bins: list[GQBin]) -> None:
    # gq_label is optional but CLI advertises it; compute if missing and bins present.
    if not bins:
        return
    if "gq_label" in out.columns:
        return
    if "gq" not in out.columns:
        return
    out["gq_label"] = out["gq"].apply(lambda x: _gq_label_from_bins(x, bins))


def _validate_required_columns(out: pd.DataFrame) -> None:
    required = ["chrom", "pos", "id", "end", "svtype", "gt", "gq"]
    missing = [c for c in required if c not in out.columns]
    if missing:
        raise RuntimeError(
            "SvPhaser internal error: worker output is missing required columns "
            f"{missing}. Refusing to write phased VCF/CSV because it would silently "
            "degrade to './.' genotypes. Fix worker to emit gt/gq "
            "(and ideally reason/delta)."
        )


def _normalize_gt_gq(out: pd.DataFrame) -> None:
    # Avoid pandas NA weirdness: make sure gt is present and string-like.
    out["gt"] = out["gt"].fillna("./.").astype(str)

    # Make sure gq is int-ish.
    out["gq"] = pd.to_numeric(out["gq"], errors="coerce").fillna(0).astype(int)


def _warn_if_suspicious(out: pd.DataFrame) -> None:
    # If almost all GT are './.' but we clearly have tagged support, warn loudly.
    if "tagged_total" not in out.columns:
        return
    try:
        frac_ambig = float((out["gt"] == "./.").mean())
        med_tagged = float(pd.to_numeric(out["tagged_total"], errors="coerce").fillna(0).median())
        if frac_ambig >= 0.95 and med_tagged >= 5:
            logger.warning(
                "Suspicious output: %.1f%% genotypes are './.' despite "
                "median tagged_total=%.1f. This usually indicates a gt/gq "
                "propagation bug or thresholds too strict.",
                100.0 * frac_ambig,
                med_tagged,
            )
    except Exception:
        return


def _ensure_required_columns(df: pd.DataFrame, *, bins: list[GQBin]) -> pd.DataFrame:
    """Backfill columns the writer expects; then validate."""
    if df.empty:
        return df

    out = df.copy()

    _backfill_n1_n2(out)
    _backfill_tag_frac(out)
    _backfill_gq_label(out, bins=bins)

    _validate_required_columns(out)
    _normalize_gt_gq(out)
    _warn_if_suspicious(out)

    return out


# ──────────────────────────────────────────────────────────────────────
#  Public entrypoint
# ──────────────────────────────────────────────────────────────────────


def phase_vcf(
    sv_vcf: Path,
    bam: Path,
    *,
    out_dir: Path,
    min_support: int = 10,
    min_tagged_support: int = 3,
    major_delta: float = 0.60,
    equal_delta: float = 0.10,
    gq_bins: str = "0:LOW,20:MED,50:HIGH",
    support_mode: str = "hybrid",
    bp_window: int = 100,
    dynamic_window: bool = True,
    tie_to_hom_alt: bool = True,
    svp_info: bool = True,
    threads: int | None = None,
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
        min_tagged_support=min_tagged_support,
        major_delta=major_delta,
        equal_delta=equal_delta,
        gq_bins=bins,
        support_mode=support_mode,
        bp_window=bp_window,
        dynamic_window=dynamic_window,
        tie_to_hom_alt=tie_to_hom_alt,
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

    # Backfill + validate before any writing.
    merged = _ensure_required_columns(merged, bins=bins)

    pre = len(merged)
    if "support_total" in merged.columns:
        total_support = (
            pd.to_numeric(merged["support_total"], errors="coerce").fillna(0).astype(int)
        )
    else:
        # Fallback: if support_total not present, treat tagged as total (legacy)
        total_support = pd.to_numeric(merged["n1"], errors="coerce").fillna(0).astype(
            int
        ) + pd.to_numeric(merged["n2"], errors="coerce").fillna(0).astype(int)

    keep = total_support >= int(min_support)

    stem = sv_vcf.name.removesuffix(".vcf.gz").removesuffix(".vcf")

    # Save dropped SVs for transparency
    dropped_csv = out_dir / f"{stem}_dropped_svs.csv"
    merged.loc[~keep].to_csv(dropped_csv, index=False)
    logger.info("Dropped SVs → %s (%d SVs)", dropped_csv, int((~keep).sum()))

    kept = merged.loc[keep].reset_index(drop=True)
    dropped = pre - len(kept)
    if dropped:
        logger.info("Support filter removed %d SVs", dropped)

    # 6 ─ Write CSV
    out_csv = out_dir / f"{stem}_phased.csv"
    kept.to_csv(out_csv, index=False)
    logger.info("CSV → %s  (%d SVs)", out_csv, len(kept))

    # 7 ─ Write VCF
    out_vcf = out_dir / f"{stem}_phased.vcf"
    _write_phased_vcf(
        out_vcf,
        sv_vcf,
        kept,
        gqbin_in_header=bool(bins),
        svp_info_in_header=svp_info,
        svp_info=svp_info,
    )
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
    svp_info_in_header: bool,
) -> None:
    """Write preserved meta headers + ensure GT/GQ/GQBIN, then the column header."""
    have_gt = any("##FORMAT=<ID=GT" in ln for ln in raw_header_lines)
    have_gq = any("##FORMAT=<ID=GQ" in ln for ln in raw_header_lines)
    have_gqbin = any("##INFO=<ID=GQBIN" in ln for ln in raw_header_lines)
    have_svp_hp1 = any("##INFO=<ID=SVP_HP1" in ln for ln in raw_header_lines)

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

    if svp_info_in_header and not have_svp_hp1:
        out.write(
            "##INFO=<ID=SVP_MODE,Number=1,Type=String,"
            'Description="SvPhaser: evidence mode used (RNAMES/HEURISTIC)">\n'
        )
        out.write(
            "##INFO=<ID=SVP_HP1,Number=1,Type=Integer,"
            'Description="SvPhaser: HP=1 supporting reads">\n'
        )
        out.write(
            "##INFO=<ID=SVP_HP2,Number=1,Type=Integer,"
            'Description="SvPhaser: HP=2 supporting reads">\n'
        )
        out.write(
            "##INFO=<ID=SVP_NOHP,Number=1,Type=Integer,"
            'Description="SvPhaser: supporting reads without HP tag">\n'
        )
        out.write(
            "##INFO=<ID=SVP_SUPPORT,Number=1,Type=Integer,"
            'Description="SvPhaser: total supporting reads (HP1+HP2+NO_HP)">\n'
        )
        out.write(
            "##INFO=<ID=SVP_TAGGED,Number=1,Type=Integer,"
            'Description="SvPhaser: supporting reads with HP tag (HP1+HP2)">\n'
        )
        out.write(
            "##INFO=<ID=SVP_TAGFRAC,Number=1,Type=Float,"
            'Description="SvPhaser: fraction of supporting reads that are HP-tagged">\n'
        )
        out.write(
            "##INFO=<ID=SVP_DELTA,Number=1,Type=Float,"
            'Description="SvPhaser: |HP1-HP2|/(HP1+HP2) on tagged support">\n'
        )
        out.write(
            "##INFO=<ID=SVP_REASON,Number=1,Type=String,"
            'Description="SvPhaser: decision reason code">\n'
        )
        out.write(
            "##INFO=<ID=SVP_FETCHW,Number=1,Type=Integer,"
            'Description="SvPhaser: fetch window (bp) used around breakpoints">\n'
        )
        out.write(
            "##INFO=<ID=SVP_BPWIN,Number=1,Type=Integer,"
            'Description="SvPhaser: breakpoint tolerance window (bp) used in evidence checks">\n'
        )
        out.write(
            "##INFO=<ID=SVP_RNAMES_TOTAL,Number=1,Type=Integer,"
            'Description="SvPhaser: RNAMES count in input VCF (if present)">\n'
        )
        out.write(
            "##INFO=<ID=SVP_RNAMES_FOUND,Number=1,Type=Integer,"
            'Description="SvPhaser: RNAMES found in BAM fetch window">\n'
        )

    out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_name + "\n")


def _compose_info_str(
    orig_info: dict[str, Any],
    svtype: Any,
    gq_label: Any,
    svp_fields: dict[str, Any] | None,
) -> str:
    """Compose INFO with SVTYPE first, proper FLAG handling, then SvPhaser fields."""
    items: list[str] = []

    if svtype:
        items.append(f"SVTYPE={svtype}")

    for k, v in orig_info.items():
        if k == "SVTYPE":
            continue
        if v is None:
            continue
        if v is True:
            items.append(str(k))
        else:
            items.append(f"{k}={v}")

    if not _is_missing_scalar(gq_label):
        items.append(f"GQBIN={gq_label}")

    if svp_fields:
        order = [
            "SVP_MODE",
            "SVP_HP1",
            "SVP_HP2",
            "SVP_NOHP",
            "SVP_SUPPORT",
            "SVP_TAGGED",
            "SVP_TAGFRAC",
            "SVP_DELTA",
            "SVP_REASON",
            "SVP_FETCHW",
            "SVP_BPWIN",
            "SVP_RNAMES_TOTAL",
            "SVP_RNAMES_FOUND",
        ]
        for k in order:
            if k not in svp_fields:
                continue
            v = svp_fields[k]
            if v is None:
                continue
            if isinstance(v, float):
                items.append(f"{k}={v:.6g}")
            else:
                items.append(f"{k}={v}")

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

    return None


def _write_phased_vcf(
    out_vcf: Path,
    in_vcf: Path,
    df: pd.DataFrame,
    *,
    gqbin_in_header: bool,
    svp_info_in_header: bool,
    svp_info: bool,
) -> None:
    """Write a phased VCF: tab-delimited, compliant, with ensured GT/GQ (and GQBIN if used)."""
    full_lookup, legacy_index, raw_header_lines, sample_name = _vcf_info_lookup(in_vcf)

    with open(out_vcf, "w", newline="") as out:
        _write_headers(
            out,
            raw_header_lines,
            sample_name,
            gqbin_in_header=gqbin_in_header,
            svp_info_in_header=svp_info_in_header,
        )

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

            svp_fields = None
            if svp_info:
                svp_fields = {
                    "SVP_MODE": getattr(row, "mode", None),
                    "SVP_HP1": int(getattr(row, "n1", 0) or 0),
                    "SVP_HP2": int(getattr(row, "n2", 0) or 0),
                    "SVP_NOHP": int(getattr(row, "nohp", 0) or 0),
                    "SVP_SUPPORT": int(getattr(row, "support_total", 0) or 0),
                    "SVP_TAGGED": int(getattr(row, "tagged_total", 0) or 0),
                    "SVP_TAGFRAC": float(getattr(row, "tag_frac", 0.0) or 0.0),
                    "SVP_DELTA": float(getattr(row, "delta", 0.0) or 0.0),
                    "SVP_REASON": getattr(row, "reason", None),
                    "SVP_FETCHW": int(getattr(row, "fetch_w", 0) or 0),
                    "SVP_BPWIN": int(getattr(row, "bp_tol", getattr(row, "bp_window", 0)) or 0),
                    "SVP_RNAMES_TOTAL": int(getattr(row, "rnames_total", 0) or 0),
                    "SVP_RNAMES_FOUND": int(getattr(row, "rnames_found", 0) or 0),
                }

            info_str = _compose_info_str(info["INFO"], svtype, gq_label, svp_fields)

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
