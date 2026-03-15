"""Tests for _supports_ins and _supports_del edge cases in svphaser.phasing._workers."""

from unittest.mock import MagicMock

from svphaser.phasing._workers import _supports_del, _supports_ins


def _make_read(cigartuples, reference_start=0, query_name="test_read"):
    """Create a minimal mock pysam.AlignedSegment."""
    read = MagicMock()
    read.cigartuples = cigartuples
    read.reference_start = reference_start
    read.query_name = query_name
    read.reference_name = "chr1"
    read.is_unmapped = False
    read.is_secondary = False
    read.is_reverse = False
    read.has_tag.return_value = False
    read.get_tag.return_value = None
    return read


class TestSupportsIns:
    """Test _supports_ins with edge cases around CIGAR operations."""

    def test_insertion_op_matches(self):
        """Standard insertion CIGAR op (1) near breakpoint should be accepted."""
        # 100bp match, then 150bp insertion at ref pos 100
        read = _make_read([(0, 100), (1, 150), (0, 50)], reference_start=0)
        assert _supports_ins(
            read,
            pos0=100,
            svlen=150,
            bp_window=100,
            size_match_required=True,
            size_tol_abs=10,
            size_tol_frac=0.0,
        )

    def test_insertion_op_too_far(self):
        """Insertion CIGAR op outside bp_window should be rejected."""
        # insertion at ref pos 100 but breakpoint at 500
        read = _make_read([(0, 100), (1, 150), (0, 50)], reference_start=0)
        assert not _supports_ins(
            read,
            pos0=500,
            svlen=150,
            bp_window=100,
            size_match_required=True,
            size_tol_abs=10,
            size_tol_frac=0.0,
        )

    def test_insertion_size_mismatch(self):
        """Insertion with wrong size should be rejected when size_match_required."""
        # 50bp insertion but expecting 200bp
        read = _make_read([(0, 100), (1, 50), (0, 50)], reference_start=0)
        assert not _supports_ins(
            read,
            pos0=100,
            svlen=200,
            bp_window=100,
            size_match_required=True,
            size_tol_abs=10,
            size_tol_frac=0.0,
        )

    def test_soft_clip_accepted(self):
        """Soft-clip (CIGAR op 4) near breakpoint should count as INS evidence."""
        # 100bp soft-clip at the start (ref pos = 0)
        read = _make_read([(4, 150), (0, 100)], reference_start=0)
        assert _supports_ins(
            read,
            pos0=0,
            svlen=150,
            bp_window=100,
            size_match_required=True,
            size_tol_abs=10,
            size_tol_frac=0.0,
        )

    def test_hard_clip_rejected(self):
        """Hard-clip (CIGAR op 5) must NOT count as insertion evidence.

        This is the core fix for Issue #2: hard-clips don't carry sequence
        data, so their length is meaningless for insertion support.
        """
        # 150bp hard-clip at the start (ref pos = 0)
        read = _make_read([(5, 150), (0, 100)], reference_start=0)
        assert not _supports_ins(
            read,
            pos0=0,
            svlen=150,
            bp_window=100,
            size_match_required=True,
            size_tol_abs=10,
            size_tol_frac=0.0,
        )

    def test_hard_clip_only_read(self):
        """Read with only hard-clip and match should not be counted."""
        read = _make_read([(5, 200), (0, 100), (5, 200)], reference_start=0)
        assert not _supports_ins(
            read,
            pos0=0,
            svlen=200,
            bp_window=100,
            size_match_required=True,
            size_tol_abs=10,
            size_tol_frac=0.0,
        )

    def test_no_cigar(self):
        """Read with no CIGAR should return False."""
        read = _make_read(None, reference_start=0)
        assert not _supports_ins(
            read,
            pos0=0,
            svlen=100,
            bp_window=100,
            size_match_required=True,
            size_tol_abs=10,
            size_tol_frac=0.0,
        )

    def test_size_match_not_required(self):
        """When size_match_required=False, size should not matter."""
        # 50bp insertion but expecting 200bp — should still pass
        read = _make_read([(0, 100), (1, 50), (0, 50)], reference_start=0)
        assert _supports_ins(
            read,
            pos0=100,
            svlen=200,
            bp_window=100,
            size_match_required=False,
            size_tol_abs=10,
            size_tol_frac=0.0,
        )

    def test_insertion_too_small(self):
        """Insertion smaller than min_len threshold should be rejected."""
        # 5bp insertion (below MIN_CIGAR_BP=30)
        read = _make_read([(0, 100), (1, 5), (0, 50)], reference_start=0)
        assert not _supports_ins(
            read,
            pos0=100,
            svlen=150,
            bp_window=100,
            size_match_required=False,
            size_tol_abs=10,
            size_tol_frac=0.0,
        )


class TestSupportsDel:
    """Test _supports_del with edge cases."""

    def test_deletion_op_matches(self):
        """DEL CIGAR op near expected breakpoints should be accepted."""
        # 100bp match then 200bp deletion then 100bp match
        read = _make_read([(0, 100), (2, 200), (0, 100)], reference_start=0)
        assert _supports_del(
            read,
            pos0=100,
            end_excl0=300,
            svlen=200,
            bp_window=100,
            size_match_required=True,
            size_tol_abs=10,
            size_tol_frac=0.0,
        )

    def test_deletion_wrong_position(self):
        """DEL far from expected breakpoints should be rejected."""
        read = _make_read([(0, 100), (2, 200), (0, 100)], reference_start=0)
        assert not _supports_del(
            read,
            pos0=500,
            end_excl0=700,
            svlen=200,
            bp_window=50,
            size_match_required=True,
            size_tol_abs=10,
            size_tol_frac=0.0,
        )

    def test_deletion_size_mismatch(self):
        """DEL with wrong size should fail with size_match_required."""
        # 50bp deletion but expecting 200bp
        read = _make_read([(0, 100), (2, 50), (0, 100)], reference_start=0)
        assert not _supports_del(
            read,
            pos0=100,
            end_excl0=300,
            svlen=200,
            bp_window=100,
            size_match_required=True,
            size_tol_abs=10,
            size_tol_frac=0.0,
        )

    def test_no_cigar(self):
        """Read with no CIGAR should return False."""
        read = _make_read(None)
        assert not _supports_del(
            read,
            pos0=100,
            end_excl0=300,
            svlen=200,
            bp_window=100,
            size_match_required=True,
            size_tol_abs=10,
            size_tol_frac=0.0,
        )
