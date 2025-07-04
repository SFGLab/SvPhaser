# src/svphaser/phasing/algorithms.py
"""Pure maths for SvPhaser – **overflow‑safe** revision

*Fix:* the earlier exact binomial tail used ``math.comb`` which explodes
for deep coverage (> ~1 500 reads).  We now:

1.  Use the **normal approximation** (continuity‑corrected) once *N* > 200.
2.  Cap any resulting GQ at 99 (VCF convention).
3.  Keep the exact combinatorial sum for small *N* so unit tests are
    unchanged.
"""
from __future__ import annotations

import math
from typing import Tuple

__all__ = [
    "classify_haplotype",
    "phasing_gq",
]

# ---------------------------------------------------------------------
#  Public helpers
# ---------------------------------------------------------------------

NORMAL_THRESHOLD = 200  # switch to Gaussian SF above this depth
MAX_GQ = 99


def phasing_gq(n1: int, n2: int) -> int:
    """Return Phred‑scaled Genotype Quality, overflow‑safe for deep data."""

    total = n1 + n2
    if total == 0:
        return 0

    k = max(n1, n2)

    if total > NORMAL_THRESHOLD:
        # Normal approximation to Binom(n,0.5)
        mu = total / 2.0
        sigma = math.sqrt(total * 0.25)
        # continuity correction: use k‑0.5 for SF
        z = (k - 0.5 - mu) / sigma
        # survival function of standard normal via erfc
        p_err = 0.5 * math.erfc(z / math.sqrt(2.0))
    else:
        # Exact tail sum – small n so safe
        p = 0.5
        tail = 0.0
        for i in range(k, total + 1):
            tail += math.comb(total, i) * (p**i) * ((1 - p) ** (total - i))
        p_err = tail

    p_err = max(p_err, 1e-300)  # guard log(0)
    gq = int(round(-10.0 * math.log10(p_err)))
    return min(gq, MAX_GQ)


def classify_haplotype(
    n1: int,
    n2: int,
    *,
    min_support: int = 10,
    major_delta: float = 0.70,
    equal_delta: float = 0.20,
) -> Tuple[str, int]:
    total = n1 + n2
    if total < min_support:
        return "./.", 0

    gq = phasing_gq(n1, n2)
    n1_ratio = n1 / total
    n2_ratio = n2 / total

    if n1_ratio >= major_delta:
        gt = "1|0"
    elif n2_ratio >= major_delta:
        gt = "0|1"
    elif abs(n1 - n2) / total <= equal_delta:
        gt = "1|1"
    else:
        gt = "./."
    return gt, gq

