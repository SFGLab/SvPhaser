"""Top-level SvPhaser package.

Public surface kept tiny:
- __version__
- a convenience `phase()` wrapper around svphaser.phasing.io.phase_vcf()

Defaults are chosen to match the recommended SvPhaser settings for long-read SV phasing.
"""

from __future__ import annotations

from pathlib import Path

# --------------------------------------------------------------------
# Robust version lookup:
# - Prefer installed package metadata (works for wheels and PEP 660 editables)
# - Fall back to _version.py for raw-source/dev use
# --------------------------------------------------------------------
try:
    from importlib.metadata import version as _pkg_version  # Python 3.8+

    __version__ = _pkg_version("svphaser")
except Exception:
    try:
        from ._version import __version__  # overwritten in builds when using setuptools-scm
    except Exception:  # highly defensive
        __version__ = "0+unknown"

# Centralized defaults (keep CLI in sync)
DEFAULT_MIN_SUPPORT: int = 10
DEFAULT_MAJOR_DELTA: float = 0.60
DEFAULT_EQUAL_DELTA: float = 0.10
DEFAULT_GQ_BINS: str = "30:High,10:Moderate"


def phase(
    sv_vcf: Path | str,
    bam: Path | str,
    /,
    *,
    out_dir: Path | str = ".",
    min_support: int = DEFAULT_MIN_SUPPORT,
    major_delta: float = DEFAULT_MAJOR_DELTA,
    equal_delta: float = DEFAULT_EQUAL_DELTA,
    gq_bins: str = DEFAULT_GQ_BINS,
    threads: int | None = None,
) -> tuple[Path, Path]:
    """Phase *sv_vcf* using HP-tagged *bam*, writing outputs into *out_dir*.

    Notes
    -----
    - Step B semantics: `min_support` is applied to TOTAL ALT-supporting reads (n1+n2).
    - Near-ties (<= equal_delta) are treated as ambiguous (./.), not homozygous ALT.

    Returns
    -------
    (out_vcf_path, out_csv_path)
    """
    from .phasing.io import phase_vcf  # local import avoids heavy deps at import-time

    out_dir_p = Path(out_dir)
    out_dir_p.mkdir(parents=True, exist_ok=True)

    stem = Path(sv_vcf).name
    if stem.endswith(".vcf.gz"):
        stem = stem[:-7]
    elif stem.endswith(".vcf"):
        stem = stem[:-4]

    out_vcf = out_dir_p / f"{stem}_phased.vcf"
    out_csv = out_dir_p / f"{stem}_phased.csv"

    phase_vcf(
        Path(sv_vcf),
        Path(bam),
        out_dir=out_dir_p,
        min_support=min_support,
        major_delta=major_delta,
        equal_delta=equal_delta,
        gq_bins=gq_bins,
        threads=threads,
    )
    return out_vcf, out_csv


__all__ = [
    "phase",
    "__version__",
    "DEFAULT_MIN_SUPPORT",
    "DEFAULT_MAJOR_DELTA",
    "DEFAULT_EQUAL_DELTA",
    "DEFAULT_GQ_BINS",
]
