"""Pure maths for SvPhaser – overflow-safe revision.

1) Exact binomial tail for small depth (N ≤ 200).
2) Continuity-corrected normal approximation for deep coverage (N > 200).
3) Phred GQ capped at 99.
"""

from __future__ import annotations

import math

__all__ = ["classify_haplotype", "phasing_gq"]

NORMAL_THRESHOLD = 200  # switch to Gaussian SF above this depth
MAX_GQ = 99


def phasing_gq(n1: int, n2: int) -> int:
    """Return Phred-scaled Genotype Quality, overflow-safe for deep data."""
    total = n1 + n2
    if total == 0:
        return 0

    k = max(n1, n2)

    if total > NORMAL_THRESHOLD:
        mu = total / 2.0
        sigma = math.sqrt(total * 0.25)
        z = (k - 0.5 - mu) / sigma
        p_err = 0.5 * math.erfc(z / math.sqrt(2.0))  # survival function
    else:
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
    equal_delta: float = 0.25,
) -> tuple[str, int]:
    """Return (GT, GQ) using ratio thresholds and an overflow-safe GQ."""
    total = n1 + n2

    if n1 < min_support and n2 < min_support:
        return "./.", 0
    if total == 0:
        return "./.", 0

    gq = phasing_gq(n1, n2)
    r1 = n1 / total
    r2 = n2 / total

    if r1 >= major_delta:
        gt = "1|0"
    elif r2 >= major_delta:
        gt = "0|1"
    elif abs(n1 - n2) / total <= equal_delta:
        gt = "1|1"
    else:
        gt = "./."
    return gt, gq
