from svphaser.phasing.algorithms import classify_haplotype, phasing_gq


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


def test_classify_thresholds_basic():
    # strong majority to HP1
    gt, gq = classify_haplotype(70, 30, min_support=1, major_delta=0.6, equal_delta=0.1)
    assert gt == "1|0" and gq >= 0
    # strong majority to HP2
    gt, _ = classify_haplotype(30, 70, min_support=1, major_delta=0.6, equal_delta=0.1)
    assert gt == "0|1"
    # near tie → 1|1
    gt, _ = classify_haplotype(51, 49, min_support=1, major_delta=0.6, equal_delta=0.05)
    assert gt == "1|1"
    # below support on both → drop
    gt, _ = classify_haplotype(0, 0, min_support=1, major_delta=0.6, equal_delta=0.1)
    assert gt == "./."
