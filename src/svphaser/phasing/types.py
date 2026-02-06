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
    # Minimum number of HP-tagged supporting reads required before we attempt
    # a directional call (1|0 or 0|1). This does *not* replace `min_support`.
    min_tagged_support: int
    major_delta: float
    equal_delta: float
    # If True, near-ties are output as 1|1 ("both haplotypes support") instead of ./.
    tie_to_hom_alt: bool
    # Read-support selection strategy.
    #  - "hybrid": use INFO/RNAMES when present, else fall back to heuristics
    #  - "rnames": require INFO/RNAMES; otherwise treat as no-support
    #  - "heuristic": ignore RNAMES even if present
    support_mode: str
    # Breakpoint tolerance (bp) used inside evidence checks; can be dynamically enlarged.
    bp_window: int
    # Fetch window strategy: if True, use SVLEN/CIs/STDEV to widen read fetching.
    dynamic_window: bool
    gq_bins: list[GQBin]


class CallTuple(NamedTuple):
    """Return type per-variant from algorithms.classify_haplotype()."""

    gt: str
    gq: int
    gq_label: str | None
