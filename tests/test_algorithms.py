"""Tests for svphaser.phasing.algorithms — phasing_gq and both classifiers."""

from svphaser.phasing.algorithms import (
    classify_haplotype,
    classify_haplotype_v211,
    phasing_gq,
)

# ─── phasing_gq ────────────────────────────────────────────────────────


def test_gq_zero_and_cap():
    assert phasing_gq(0, 0) == 0
    assert phasing_gq(1000, 0) == 99


def test_gq_symmetry_and_monotonicity():
    assert phasing_gq(37, 13) == phasing_gq(13, 37)
    # monotonic with |n1-n2| at fixed N
    N = 20
    deltas = [0, 2, 4, 6, 8, 10]
    prev = None
    for d in deltas:
        n1 = (N + d) // 2
        n2 = N - n1
        gq = phasing_gq(n1, n2)
        if prev is not None:
            assert gq >= prev
        prev = gq


# ─── classify_haplotype (legacy) ───────────────────────────────────────


def test_classify_strong_majority_hp1():
    gt, gq = classify_haplotype(70, 30, min_support=1, major_delta=0.6, equal_delta=0.1)
    assert gt == "1|0" and gq >= 0


def test_classify_strong_majority_hp2():
    gt, _ = classify_haplotype(30, 70, min_support=1, major_delta=0.6, equal_delta=0.1)
    assert gt == "0|1"


def test_classify_near_tie_hom_alt():
    """Default tie_to_hom_alt=True: near-tie with both haplotypes → 1|1."""
    gt, _ = classify_haplotype(51, 49, min_support=1, major_delta=0.6, equal_delta=0.05)
    assert gt == "1|1"


def test_classify_near_tie_ambiguous():
    """With tie_to_hom_alt=False: near-tie → ./."""
    gt, _ = classify_haplotype(
        51,
        49,
        min_support=1,
        major_delta=0.6,
        equal_delta=0.05,
        tie_to_hom_alt=False,
    )
    assert gt == "./."


def test_classify_no_support():
    gt, _ = classify_haplotype(0, 0, min_support=1, major_delta=0.6, equal_delta=0.1)
    assert gt == "./."


def test_classify_below_min_support():
    gt, _ = classify_haplotype(3, 2, min_support=10, major_delta=0.6, equal_delta=0.1)
    assert gt == "./."


def test_classify_below_min_tagged_support():
    gt, _ = classify_haplotype(
        2,
        1,
        min_support=1,
        major_delta=0.6,
        equal_delta=0.1,
        min_tagged_support=5,
    )
    assert gt == "./."


def test_classify_ambiguous_zone():
    """Counts in the gap between equal_delta and major_delta → ./."""
    # 60/40 split: ratio=0.6 but |diff|/total = 0.2 which is > equal_delta(0.05)
    # and max_ratio = 0.6 which is NOT >= major_delta(0.7)
    gt, _ = classify_haplotype(60, 40, min_support=1, major_delta=0.7, equal_delta=0.05)
    assert gt == "./."


# ─── classify_haplotype_v211 ──────────────────────────────────────────


def test_v211_strong_majority_hp1():
    gt, gq, reason, delta = classify_haplotype_v211(
        n1=70,
        n2=10,
        min_support=10,
        min_tagged_support=3,
        major_delta=0.6,
        equal_delta=0.1,
        support_total=80,
        tie_to_hom_alt=True,
    )
    assert gt == "1|0"
    assert reason == "MAJOR_HP1"
    assert gq > 0


def test_v211_strong_majority_hp2():
    gt, _, reason, _ = classify_haplotype_v211(
        n1=10,
        n2=70,
        min_support=10,
        min_tagged_support=3,
        major_delta=0.6,
        equal_delta=0.1,
        support_total=80,
        tie_to_hom_alt=True,
    )
    assert gt == "0|1"
    assert reason == "MAJOR_HP2"


def test_v211_tie_hom_alt():
    gt, _, reason, delta = classify_haplotype_v211(
        n1=25,
        n2=25,
        min_support=10,
        min_tagged_support=3,
        major_delta=0.6,
        equal_delta=0.1,
        support_total=50,
        tie_to_hom_alt=True,
    )
    assert gt == "1|1"
    assert reason == "TIE"
    assert delta == 0.0


def test_v211_tie_ambiguous():
    gt, _, reason, _ = classify_haplotype_v211(
        n1=25,
        n2=25,
        min_support=10,
        min_tagged_support=3,
        major_delta=0.6,
        equal_delta=0.1,
        support_total=50,
        tie_to_hom_alt=False,
    )
    assert gt == "./."
    assert reason == "TIE"


def test_v211_no_support():
    gt, gq, reason, _ = classify_haplotype_v211(
        n1=0,
        n2=0,
        min_support=10,
        min_tagged_support=3,
        major_delta=0.6,
        equal_delta=0.1,
        support_total=0,
        tie_to_hom_alt=True,
    )
    assert gt == "./."
    assert reason == "NO_SUPPORT"
    assert gq == 0


def test_v211_low_support():
    gt, _, reason, _ = classify_haplotype_v211(
        n1=3,
        n2=2,
        min_support=10,
        min_tagged_support=3,
        major_delta=0.6,
        equal_delta=0.1,
        support_total=5,
        tie_to_hom_alt=True,
    )
    assert gt == "./."
    assert reason == "LOW_SUPPORT"


def test_v211_no_tagged():
    gt, _, reason, _ = classify_haplotype_v211(
        n1=0,
        n2=0,
        min_support=10,
        min_tagged_support=3,
        major_delta=0.6,
        equal_delta=0.1,
        support_total=15,
        tie_to_hom_alt=True,
    )
    assert gt == "./."
    assert reason == "NO_TAGGED"


def test_v211_low_tagged():
    gt, _, reason, _ = classify_haplotype_v211(
        n1=1,
        n2=1,
        min_support=10,
        min_tagged_support=5,
        major_delta=0.6,
        equal_delta=0.1,
        support_total=15,
        tie_to_hom_alt=True,
    )
    assert gt == "./."
    assert reason == "LOW_TAGGED"


def test_v211_ambiguous():
    """Counts in gap between equal_delta and major_delta → AMBIGUOUS."""
    gt, _, reason, _ = classify_haplotype_v211(
        n1=35,
        n2=25,
        min_support=10,
        min_tagged_support=3,
        major_delta=0.7,
        equal_delta=0.05,
        support_total=60,
        tie_to_hom_alt=True,
    )
    assert gt == "./."
    assert reason == "AMBIGUOUS"
