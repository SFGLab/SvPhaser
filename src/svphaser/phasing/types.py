"""svphaser.phasing.types
========================
Central place for common type aliases & lightweight data classes.
Keeping them here avoids circular imports and MyPy noise.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import NamedTuple

SVKey = tuple[str, int, str]  # (chrom, POS, ID) – ID is "." if empty
GQBin = tuple[int, str]  # (threshold, label), e.g. (30, "High")


@dataclass(slots=True, frozen=True)
class WorkerOpts:
    """Non-changing knobs passed into every worker."""

    min_support: int
    major_delta: float
    equal_delta: float
    gq_bins: list[GQBin]  # already parsed by cli → phase_vcf


class CallTuple(NamedTuple):
    """Return type per-variant from algorithms.classify_haplotype()."""

    gt: str
    gq: int
    gq_label: str | None
