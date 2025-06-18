"""Top-level SvPhaser package.

This file keeps the public surface tiny: a version string and a single
convenience helper that calls the library’s main phasing routine.
"""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

try:  # when installed via pip / wheel
    __version__: str = version("svphaser")
except PackageNotFoundError:  # editable install without a built dist
    __version__ = "0.0.0+dev"

# --------------------------------------------------------------------
#  Convenience façade – lets users do `svphaser.phase(…)` in notebooks
# --------------------------------------------------------------------
from pathlib import Path
from typing import Optional

def phase(
    sv_vcf: Path | str,
    bam: Path | str,
    /,
    *,
    out_vcf: Path | str = "phased.vcf",
    min_support: int = 10,
    major_delta: float = 0.70,
    equal_delta: float = 0.20,
    gq_bins: str = "30:High,10:Moderate",
) -> Optional[int]:
    """Phase *sv_vcf* using HP-tagged *bam*.

    This is a very thin wrapper around :pyfunc:`svphaser.phasing.io.phase_vcf`
    so that casual users (or tests) can skip importing sub-modules.
    """
    from .phasing.io import phase_vcf  # local import avoids heavy deps at import-time

    return phase_vcf(
        sv_vcf,
        bam,
        out_vcf=out_vcf, # type: ignore
        min_support=min_support,
        major_delta=major_delta,
        equal_delta=equal_delta,
        gq_bins=gq_bins, # type: ignore
    )

__all__ = ["phase", "__version__"]
