"""Top-level SvPhaser package.

Public surface kept small:
- __version__
- a convenience `phase()` wrapper around svphaser.phasing.io.phase_vcf()

Versioning:
- Primary: installed package metadata (works for wheels and PEP 660 editables).
- Fallback: "0+unknown" when running from an uninstalled source tree.

Note: SvPhaser uses hatch-vcs for versioning. Do not rely on a generated _version.py
file at runtime.
"""

from __future__ import annotations

from pathlib import Path

try:
    # Python 3.8+: preferred source for installed distributions
    from importlib.metadata import PackageNotFoundError
    from importlib.metadata import version as _pkg_version
except Exception:  # pragma: no cover
    PackageNotFoundError = Exception  # type: ignore[assignment]
    _pkg_version = None  # type: ignore[assignment]

try:
    if _pkg_version is None:
        raise PackageNotFoundError
    __version__ = _pkg_version("svphaser")
except Exception:
    # Source-tree fallback (e.g., running without installing the package)
    __version__ = "0+unknown"


# Centralized defaults (keep CLI in sync)
DEFAULT_MIN_SUPPORT: int = 10
DEFAULT_MIN_TAGGED_SUPPORT: int = 3
DEFAULT_MAJOR_DELTA: float = 0.60
DEFAULT_EQUAL_DELTA: float = 0.10
DEFAULT_GQ_BINS: str = "30:High,10:Moderate"
DEFAULT_SUPPORT_MODE: str = "hybrid"
DEFAULT_BP_WINDOW: int = 100
DEFAULT_DYNAMIC_WINDOW: bool = True
DEFAULT_TIE_TO_HOM_ALT: bool = True
DEFAULT_SVP_INFO: bool = True


def phase(
    sv_vcf: Path | str,
    bam: Path | str,
    /,
    *,
    out_dir: Path | str = ".",
    min_support: int = DEFAULT_MIN_SUPPORT,
    min_tagged_support: int = DEFAULT_MIN_TAGGED_SUPPORT,
    major_delta: float = DEFAULT_MAJOR_DELTA,
    equal_delta: float = DEFAULT_EQUAL_DELTA,
    gq_bins: str = DEFAULT_GQ_BINS,
    support_mode: str = DEFAULT_SUPPORT_MODE,
    bp_window: int = DEFAULT_BP_WINDOW,
    dynamic_window: bool = DEFAULT_DYNAMIC_WINDOW,
    tie_to_hom_alt: bool = DEFAULT_TIE_TO_HOM_ALT,
    svp_info: bool = DEFAULT_SVP_INFO,
    threads: int | None = None,
    size_match_required: bool = True,
    size_tol_abs: int = 10,
    size_tol_frac: float = 0.0,
) -> tuple[Path, Path]:
    """Phase *sv_vcf* using HP-tagged *bam*, writing outputs into *out_dir*.

    Semantics (matches current SvPhaser behavior)
    ---------------------------------------------
    - Support is ALT-support evidence counted as:
        hp1 (HP=1), hp2 (HP=2), nohp (missing HP tag)
      and:
        tagged_total = hp1 + hp2
        support_total = hp1 + hp2 + nohp
    - `min_support` is applied to support_total (HP + NOHP). Below threshold, SVs are dropped.
    - Directional phasing (1|0 / 0|1) requires at least `min_tagged_support` HP-tagged reads.
    - Strong majority uses `major_delta = max(hp1,hp2)/tagged_total`.
    - Near-ties use `equal_delta = |hp1-hp2|/tagged_total`:
        * if tie_to_hom_alt=True: emit 1|1
        * else: emit ./.

    Returns
    -------
    (out_vcf_path, out_csv_path)
    """
    # Local import avoids heavy deps at import-time
    from .phasing.io import phase_vcf

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
        size_match_required=size_match_required,
        size_tol_abs=size_tol_abs,
        size_tol_frac=size_tol_frac,
    )
    return out_vcf, out_csv


__all__ = [
    "phase",
    "__version__",
    "DEFAULT_MIN_SUPPORT",
    "DEFAULT_MIN_TAGGED_SUPPORT",
    "DEFAULT_MAJOR_DELTA",
    "DEFAULT_EQUAL_DELTA",
    "DEFAULT_GQ_BINS",
    "DEFAULT_SUPPORT_MODE",
    "DEFAULT_BP_WINDOW",
    "DEFAULT_DYNAMIC_WINDOW",
    "DEFAULT_TIE_TO_HOM_ALT",
    "DEFAULT_SVP_INFO",
]
