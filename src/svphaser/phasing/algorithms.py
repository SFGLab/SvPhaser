"""Pure maths for SvPhaser – overflow-safe revision.

1) Exact binomial tail for small depth (N ≤ 200).
2) Continuity-corrected normal approximation for deep coverage (N > 200).
3) Phred GQ capped at 99.

Step B semantics (V2.1.1):
- `min_support` is interpreted as a *total ALT-support* threshold
  (HP1 + HP2 + NO_HP among supporting reads).
- Directional phasing (1|0 / 0|1) uses only HP-tagged supporting reads.
- Near-ties (|HP1-HP2|/tagged_total <= equal_delta) are interpreted as
  "both haplotypes support" and emitted as 1|1 by default.
"""

from __future__ import annotations

import math

__all__ = ["classify_haplotype", "classify_haplotype_v211", "phasing_gq"]

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
    major_delta: float = 0.60,
    equal_delta: float = 0.10,
    support_total: int | None = None,
    min_tagged_support: int = 0,
    tie_to_hom_alt: bool = True,
) -> tuple[str, int]:
    """Classify which haplotype carries the ALT allele.

    Returns:
      - GT: "1|0" (ALT on hap1) or "0|1" (ALT on hap2) or "./." (ambiguous)
      - GQ: phred-scaled confidence based on haplotype imbalance

    Notes:
      - `min_support` is applied to total ALT-support reads.
      - Near-ties (<= equal_delta) are treated as ambiguous (./.).
    """

    tagged_total = n1 + n2
    if tagged_total <= 0:
        return "./.", 0
    # Step B: `min_support` is a *total ALT-support* gate (tagged + untagged).
    # If `support_total` wasn't provided, fall back to tagged depth (legacy behavior).
    total_support = tagged_total if support_total is None else int(support_total)
    if total_support < min_support:
        return "./.", 0

    if tagged_total < min_tagged_support:
        return "./.", 0

    gq = phasing_gq(n1, n2)

    # 1) near-tie → ambiguous
    if abs(n1 - n2) / tagged_total <= equal_delta:
        # V2.1.1 default: interpret near-tie as "both haplotypes support".
        if tie_to_hom_alt and n1 > 0 and n2 > 0:
            return "1|1", gq
        return "./.", gq

    # 2) strong majority → heterozygous phased
    r1 = n1 / tagged_total
    r2 = n2 / tagged_total
    if r1 >= major_delta:
        return "1|0", gq
    if r2 >= major_delta:
        return "0|1", gq

    # 3) otherwise ambiguous
    return "./.", gq


def classify_haplotype_v211(
    n1: int,
    n2: int,
    *,
    min_support: int,
    min_tagged_support: int,
    major_delta: float,
    equal_delta: float,
    support_total: int,
    tie_to_hom_alt: bool,
) -> tuple[str, int, str, float]:
    """V2.1.1 classification with reason codes.

    Returns:
      (gt, gq, reason, delta)

    Where:
      - reason is one of: NO_SUPPORT, LOW_SUPPORT, NO_TAGGED, LOW_TAGGED,
        TIE, MAJOR_HP1, MAJOR_HP2, AMBIGUOUS
      - delta = |n1-n2|/tagged_total (0..1), or 1.0 if tagged_total==0
    """
    tagged_total = n1 + n2
    if support_total <= 0:
        return "./.", 0, "NO_SUPPORT", 1.0
    if support_total < min_support:
        return "./.", 0, "LOW_SUPPORT", 1.0
    if tagged_total <= 0:
        return "./.", 0, "NO_TAGGED", 1.0
    if tagged_total < min_tagged_support:
        return "./.", 0, "LOW_TAGGED", 1.0

    gq = phasing_gq(n1, n2)
    delta = abs(n1 - n2) / tagged_total

    if delta <= equal_delta:
        if tie_to_hom_alt and n1 > 0 and n2 > 0:
            return "1|1", gq, "TIE", delta
        return "./.", gq, "TIE", delta

    r1 = n1 / tagged_total
    r2 = n2 / tagged_total
    if r1 >= major_delta:
        return "1|0", gq, "MAJOR_HP1", delta
    if r2 >= major_delta:
        return "0|1", gq, "MAJOR_HP2", delta

    return "./.", gq, "AMBIGUOUS", delta
