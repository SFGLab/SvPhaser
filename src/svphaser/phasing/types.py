"""svphaser.phasing.types
========================
Common type aliases & small data structures.

We keep this module light to avoid circular imports.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import NamedTuple

# Legacy key (older writer used this; can collide when ID='.' or same POS repeats)
SVKeyLegacy = tuple[str, int, str]  # (CHROM, POS, ID)

# Collision-resistant key for matching phased rows back to original VCF records
SVKey = tuple[str, int, str, int, str]  # (CHROM, POS, ID, END, ALT)

# GQ bin spec: (threshold, label)
GQBin = tuple[int, str]  # e.g. (30, "High")


@dataclass(slots=True, frozen=True)
class WorkerOpts:
    """Non-changing knobs passed into every worker."""

    min_support: int
    major_delta: float
    equal_delta: float
    gq_bins: list[GQBin]


class CallTuple(NamedTuple):
    """Return type per-variant from algorithms.classify_haplotype()."""

    gt: str
    gq: int
    gq_label: str | None
