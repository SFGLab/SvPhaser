#!/usr/bin/env python3
"""
svphaser.cli
============
Command-line interface for **SvPhaser**.

The program writes two files inside **--out-dir** (or the CWD):

* ``<stem>_phased.vcf``   (uncompressed; HP_GT / HP_GQ / HP_GQBIN injected)
* ``<stem>_phased.csv``   (tabular summary including gq_label column)

"""

from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

app = typer.Typer(add_completion=False, rich_markup_mode="rich")

# ──────────────────────────────────────────────────────────────────────────
#  phase command
# ──────────────────────────────────────────────────────────────────────────
@app.command("phase")
def phase_cmd(
    sv_vcf: Annotated[
        Path,
        typer.Argument(
            exists=True,
            help="Input *un-phased* SV VCF (.vcf or .vcf.gz)",
        ),
    ],
    bam: Annotated[
        Path,
        typer.Argument(
            exists=True,
            help="Long-read BAM/CRAM with HP tags",
        ),
    ],
    out_dir: Annotated[
        Path,
        typer.Option(
            "--out-dir",
            "-o",
            exists=False,
            file_okay=False,
            dir_okay=True,
            writable=True,
            help=(
                "Directory in which to write <stem>_phased.vcf & .csv "
                "(created if missing; defaults to current dir)."
            ),
        ),
    ] = Path("."),
    # ---------- thresholds ------------------------------------------------
    min_support: Annotated[
        int,
        typer.Option(
            help=(
                "Minimum HP-tagged reads per haplotype. "
                "SVs where *both* n1 AND n2 fall below this "
                "are now dropped entirely."
            ),
        ),
    ] = 10,
    major_delta: Annotated[
        float,
        typer.Option(
            help="Δ/N > this ⇒ strong majority ⇒ GT 1|0 or 0|1",
        ),
    ] = 0.70,
    equal_delta: Annotated[  # ▲ tweaked wording
        float,
        typer.Option(
            help="Δ/N ≤ this ⇒ near-tie ⇒ GT 1|1",
        ),
    ] = 0.25,
    # ---------- confidence bins ------------------------------------------
    gq_bins: Annotated[
        str,
        typer.Option(
            help=(
                "Comma-separated GQ≥threshold:Label definitions "
                "(e.g. '30:High,10:Moderate').  Labels appear in the CSV "
                "column [gq_label] and the VCF INFO field HP_GQBIN."
            ),
        ),
    ] = "30:High,10:Moderate",
    # ---------- multiprocessing ------------------------------------------
    threads: Annotated[
        int | None,
        typer.Option(
            "-t",
            "--threads",
            help="Worker processes to use (defaults to all CPU cores).",
        ),
    ] = None,
) -> None:
    """Phase structural variants using HP-tagged read evidence."""
    # ── Initialise logging BEFORE we import anything that might log ──
    from svphaser.logging import init as _init_logging

    _init_logging("INFO")       # or "DEBUG" if you want more detail

    # ── Resolve output paths --------------------------------------------
    if not out_dir.exists():
        out_dir.mkdir(parents=True)

    stem = sv_vcf.name
    if stem.endswith(".vcf.gz"):
        stem = stem[:-7]
    elif stem.endswith(".vcf"):
        stem = stem[:-4]

    out_vcf = out_dir / f"{stem}_phased.vcf"
    out_csv = out_dir / f"{stem}_phased.csv"

    # Lazy import so `svphaser --help` works without heavy deps
    from svphaser.phasing.io import phase_vcf

    try:
        phase_vcf(
            sv_vcf,
            bam,
            out_dir=out_dir,           # type: ignore[arg-type]
            min_support=min_support,
            major_delta=major_delta,
            equal_delta=equal_delta,
            gq_bins=gq_bins,           # type: ignore[arg-type]
            threads=threads,
        )
        typer.secho(f"✔ Phased VCF → {out_vcf}", fg=typer.colors.GREEN)
        typer.secho(f"✔ Phased CSV → {out_csv}", fg=typer.colors.GREEN)
    except Exception as exc:  # pragma: no cover
        typer.secho(
            "[SvPhaser] 💥  Unhandled error during phasing",
            fg=typer.colors.RED,
        )
        raise exc
