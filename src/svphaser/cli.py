#!/usr/bin/env python3
"""svphaser.cli
============
Command-line interface for **SvPhaser**.

The program writes two files inside **--out-dir** (or the CWD):

* ``<stem>_phased.vcf``   (uncompressed; GT/GQ injected; optional INFO=GQBIN)
* ``<stem>_phased.csv``   (tabular summary incl. n1/n2/gt/gq and optional gq_label)
"""

from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from svphaser import (
    DEFAULT_EQUAL_DELTA,
    DEFAULT_GQ_BINS,
    DEFAULT_MAJOR_DELTA,
    DEFAULT_MIN_SUPPORT,
    __version__,
)

app = typer.Typer(add_completion=False, rich_markup_mode="rich")


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(__version__)
        raise typer.Exit()


@app.callback()
def main(
    version: Annotated[
        bool | None,
        typer.Option(
            "--version",
            help="Show SvPhaser version and exit.",
            is_flag=True,
            callback=_version_callback,
        ),
    ] = None,
) -> None:
    """SvPhaser – Structural-variant phasing from HP-tagged long-read BAMs."""
    return


@app.command("phase")
def phase_cmd(
    sv_vcf: Annotated[
        Path,
        typer.Argument(exists=True, help="Input *un-phased* SV VCF (.vcf or .vcf.gz)"),
    ],
    bam: Annotated[
        Path,
        typer.Argument(exists=True, help="Long-read BAM/CRAM with HP tags"),
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
            show_default=True,
        ),
    ] = Path("."),
    # ---------- thresholds ------------------------------------------------
    min_support: Annotated[
        int,
        typer.Option(
            help=(
                "Minimum TOTAL ALT-supporting reads required to keep an SV. "
                "Support is counted as HP1+HP2+NO_HP among SV-supporting reads. "
                "If support < min_support the SV is dropped (written to *_dropped_svs.csv)."
            ),
            show_default=True,
        ),
    ] = DEFAULT_MIN_SUPPORT,
    min_tagged_support: Annotated[
        int,
        typer.Option(
            help=(
                "Minimum number of HP-tagged supporting reads required before attempting a "
                "directional call (1|0 or 0|1). This is separate from --min-support."
            ),
            show_default=True,
        ),
    ] = 3,
    major_delta: Annotated[
        float,
        typer.Option(
            help="max(n1,n2)/N >= this ⇒ strong majority ⇒ GT 1|0 or 0|1",
            show_default=True,
        ),
    ] = DEFAULT_MAJOR_DELTA,
    equal_delta: Annotated[
        float,
        typer.Option(
            help=(
                "|n1−n2|/N <= this ⇒ near-tie. "
                "With --tie-to-hom-alt this emits 1|1; otherwise it emits ./."
            ),
            show_default=True,
        ),
    ] = DEFAULT_EQUAL_DELTA,
    # ---------- support selection / windows -----------------------------
    support_mode: Annotated[
        str,
        typer.Option(
            "--support-mode",
            help=(
                "How to define the SV-supporting read set: "
                "'hybrid' (use RNAMES if present, else heuristics), "
                "'rnames' (require RNAMES), or 'heuristic' (ignore RNAMES)."
            ),
            show_default=True,
        ),
    ] = "hybrid",
    bp_window: Annotated[
        int,
        typer.Option(
            "--bp-window",
            help=(
                "Base breakpoint tolerance (bp) used for evidence matching. Dynamic mode "
                "may expand it."
            ),
            show_default=True,
        ),
    ] = 100,
    dynamic_window: Annotated[
        bool,
        typer.Option(
            "--dynamic-window/--fixed-window",
            help="Use dynamic fetch windows based on SVLEN and positional uncertainty.",
            show_default=True,
        ),
    ] = True,
    tie_to_hom_alt: Annotated[
        bool,
        typer.Option(
            "--tie-to-hom-alt/--tie-to-ambig",
            help="In the near-tie zone, emit 1|1 instead of ./. .",
            show_default=True,
        ),
    ] = True,
    svp_info: Annotated[
        bool,
        typer.Option(
            "--svp-info/--no-svp-info",
            help="Write SvPhaser INFO annotations (SVP_*) into the phased VCF.",
            show_default=True,
        ),
    ] = True,
    # ---------- confidence bins ------------------------------------------
    gq_bins: Annotated[
        str,
        typer.Option(
            help=(
                "Comma-separated GQ≥threshold:Label definitions (e.g. '30:High,10:Moderate'). "
                "Labels appear in CSV column [gq_label] and in the VCF INFO field GQBIN."
            ),
            show_default=True,
        ),
    ] = DEFAULT_GQ_BINS,
    # ---------- multiprocessing ------------------------------------------
    threads: Annotated[
        int | None,
        typer.Option(
            "-t",
            "--threads",
            help="Worker processes to use (defaults to all CPU cores).",
            show_default=True,
        ),
    ] = None,
) -> None:
    """Phase structural variants using SV-type-aware ALT-support evidence."""
    from svphaser.logging import init as _init_logging

    _init_logging("INFO")

    if not out_dir.exists():
        out_dir.mkdir(parents=True)

    stem = sv_vcf.name
    if stem.endswith(".vcf.gz"):
        stem = stem[:-7]
    elif stem.endswith(".vcf"):
        stem = stem[:-4]

    out_vcf = out_dir / f"{stem}_phased.vcf"
    out_csv = out_dir / f"{stem}_phased.csv"

    from svphaser.phasing.io import phase_vcf

    try:
        phase_vcf(
            sv_vcf,
            bam,
            out_dir=out_dir,
            min_support=min_support,
            min_tagged_support=min_tagged_support,
            major_delta=major_delta,
            equal_delta=equal_delta,
            gq_bins=gq_bins,
            support_mode=support_mode,
            bp_window=bp_window,
            dynamic_window=dynamic_window,
            tie_to_hom_alt=tie_to_hom_alt,
            svp_info=svp_info,
            threads=threads,
        )
        typer.secho(f"✔ Phased VCF → {out_vcf}", fg=typer.colors.GREEN)
        typer.secho(f"✔ Phased CSV → {out_csv}", fg=typer.colors.GREEN)
    except Exception:
        typer.secho("[SvPhaser] 💥  Unhandled error during phasing", fg=typer.colors.RED)
        raise
