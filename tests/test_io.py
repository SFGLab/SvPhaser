"""Tests for _compose_info_str and _vcf_safe_value in svphaser.phasing.io."""

from svphaser.phasing.io import _compose_info_str, _vcf_safe_value


class TestVcfSafeValue:
    def test_tuple_of_floats(self):
        """COVERAGE-like tuple must become comma-separated, no parens/spaces."""
        val = (22.0, 20.0, 19.0, 20.0, 17.0)
        result = _vcf_safe_value(val)
        assert result == "22,20,19,20,17"
        assert "(" not in result
        assert " " not in result

    def test_list_of_ints(self):
        val = [1, 2, 3]
        assert _vcf_safe_value(val) == "1,2,3"

    def test_tuple_of_mixed(self):
        val = (0.0, 0.0, 73.0, 108.0, 108.0)
        result = _vcf_safe_value(val)
        assert result == "0,0,73,108,108"

    def test_float_no_trailing_zeros(self):
        assert _vcf_safe_value(0.0) == "0"
        assert _vcf_safe_value(1.5) == "1.5"
        assert _vcf_safe_value(0.699999988079071) == "0.7"

    def test_integer(self):
        assert _vcf_safe_value(42) == "42"

    def test_string(self):
        assert _vcf_safe_value("DEL") == "DEL"

    def test_empty_tuple(self):
        assert _vcf_safe_value(()) == ""

    def test_singleton_tuple(self):
        assert _vcf_safe_value((5.0,)) == "5"


class TestComposeInfoStr:
    def test_basic_scalar(self):
        info = {"END": 100, "SVLEN": -200}
        result = _compose_info_str(info, "DEL", None, None)
        assert result == "SVTYPE=DEL;END=100;SVLEN=-200"

    def test_flag_value(self):
        info = {"PRECISE": True, "END": 500}
        result = _compose_info_str(info, "INS", None, None)
        assert "PRECISE" in result
        assert "PRECISE=" not in result  # flag must not have value

    def test_tuple_coverage(self):
        """Reproduction of Issue #1: COVERAGE tuple must not have parens/spaces."""
        info = {"COVERAGE": (22.0, 20.0, 19.0, 20.0, 17.0), "END": 372751}
        result = _compose_info_str(info, "DEL", None, None)
        assert "COVERAGE=22,20,19,20,17" in result
        assert "(" not in result
        assert ")" not in result
        assert ", " not in result  # no Python-style space after comma

    def test_none_values_skipped(self):
        info = {"END": 100, "MISSING": None}
        result = _compose_info_str(info, "DEL", None, None)
        assert "MISSING" not in result

    def test_svtype_comes_first(self):
        info = {"END": 100, "SUPPORT": 14}
        result = _compose_info_str(info, "DEL", None, None)
        assert result.startswith("SVTYPE=DEL")

    def test_svtype_not_duplicated(self):
        info = {"SVTYPE": "DEL", "END": 100}
        result = _compose_info_str(info, "DEL", None, None)
        assert result.count("SVTYPE=DEL") == 1

    def test_gq_label(self):
        info = {"END": 100}
        result = _compose_info_str(info, "DEL", "High", None)
        assert "GQBIN=High" in result

    def test_svp_fields(self):
        info = {"END": 100}
        svp = {"SVP_HP1": 10, "SVP_HP2": 5, "SVP_DELTA": 0.333333}
        result = _compose_info_str(info, "DEL", None, svp)
        assert "SVP_HP1=10" in result
        assert "SVP_HP2=5" in result
        assert "SVP_DELTA=" in result

    def test_empty_info(self):
        result = _compose_info_str({}, None, None, None)
        assert result == "."
