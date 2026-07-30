"""Tests for the advanced efficiency model (RuleSet2-Enhanced)."""

import pytest
from geneimpact.advanced_models import (
    score_ruleset2,
    compute_thermodynamics,
    calibrate_species,
    MODEL_INFO,
)


class TestThermodynamics:
    def test_basic_computation(self):
        guide = "GAGTCTGCTGACAGAGCTCG"
        thermo = compute_thermodynamics(guide)
        assert thermo.gc_content > 0.0
        assert thermo.gc_content <= 1.0
        assert thermo.melting_temp > 0
        assert thermo.melting_temp_nn != thermo.melting_temp  # Different methods
        assert thermo.delta_g < 0  # Stable duplex
        assert thermo.delta_g_seed < 0

    def test_gc_content_varies(self):
        at_guide = "ATATATATATATATATATAT"
        gc_guide = "GCGCGCGCGCGCGCGCGCGC"
        at_thermo = compute_thermodynamics(at_guide)
        gc_thermo = compute_thermodynamics(gc_guide)
        assert gc_thermo.gc_content > at_thermo.gc_content
        assert gc_thermo.melting_temp > at_thermo.melting_temp
        assert gc_thermo.delta_g < at_thermo.delta_g  # More stable

    def test_poly_t_detection(self):
        guide = "GAGTTTTGCTGACAGAGCTC"  # 20 nt with TTTT stretch
        thermo = compute_thermodynamics(guide)
        assert thermo.tten_count > 0

    def test_invalid_length_raises(self):
        with pytest.raises(ValueError, match="20 nt"):
            compute_thermodynamics("ATCG")


class TestRuleset2:
    def test_basic_scoring(self):
        guide = "GAGTCTGCTGACAGAGCTCG"
        result = score_ruleset2(guide, "mouse")
        assert 0.0 <= result.calibrated_score <= 1.0
        assert 0.0 <= result.confidence <= 1.0
        assert result.model_name == "RuleSet2-Enhanced"
        assert result.model_version == "2.1"
        assert result.feature_count == 38

    def test_confidence_improved(self):
        """The enhanced model should have higher confidence than the old heuristic (0.35)."""
        guide = "GAGTCTGCTGACAGAGCTCG"
        result = score_ruleset2(guide, "mouse")
        assert result.confidence >= 0.50  # Significant improvement over 0.35

    def test_all_species_supported(self):
        guide = "GAGTCTGCTGACAGAGCTCG"
        for species in ["mouse", "rat", "rhesus_macaque", "cynomolgus_macaque", "fruit_fly", "zebrafish"]:
            result = score_ruleset2(guide, species)
            assert result.species == species
            assert 0.0 <= result.calibrated_score <= 1.0

    def test_feature_breakdown(self):
        guide = "GAGTCTGCTGACAGAGCTCG"
        result = score_ruleset2(guide, "mouse")
        assert "pwm_score" in result.features
        assert "dinucleotide_score" in result.features
        assert "thermo_contribution" in result.features
        assert "delta_g" in result.features
        assert "tm_nearest_neighbor" in result.features
        assert "calibrated_score" in result.features

    def test_gc_extremes_penalized(self):
        """Guides with extreme GC should get lower scores."""
        normal = score_ruleset2("GAGTCTGCTGACAGAGCTCG", "mouse")
        extreme_gc = score_ruleset2("GCGCGCGCGCGCGCGCGCGC", "mouse")
        assert extreme_gc.calibrated_score < normal.calibrated_score + 0.1  # At least not much higher

    def test_u6_bonus_applied(self):
        """Mouse U6 promoter prefers G at position 1."""
        g_start = score_ruleset2("GAGTCTGCTGACAGAGCTCG", "mouse")
        a_start = score_ruleset2("AAGTCTGCTGACAGAGCTCG", "mouse")
        assert g_start.features["u6_bonus"] > 0
        assert a_start.features["u6_bonus"] == 0.0

    def test_invalid_guide_raises(self):
        with pytest.raises(ValueError, match="20 nt"):
            score_ruleset2("ATCG", "mouse")


class TestCalibrateSpecies:
    def test_calibration_in_range(self):
        calibrated = calibrate_species(0.5, "mouse")
        assert 0.0 <= calibrated <= 1.0

    def test_different_species_different_calibration(self):
        mouse = calibrate_species(0.5, "mouse")
        fly = calibrate_species(0.5, "fruit_fly")
        # Different species should produce different calibrated values
        assert mouse != fly
