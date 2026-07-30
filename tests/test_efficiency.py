"""Tests for editing efficiency prediction (efficiency module)."""

from __future__ import annotations

import pytest

from geneimpact.efficiency import (
    EfficiencyPrediction,
    IndelOutcome,
    predict_efficiency,
    predict_indel_outcomes,
    compute_evidence_scores,
    SPECIES_EFFICIENCY_MODELS,
)
from geneimpact.sgrna_design import NucleaseType, SgrnaCandidate, compute_guide_features
from geneimpact.offtarget import OffTargetReport, OffTargetSite, MismatchPattern


@pytest.fixture
def good_candidate():
    """A well-designed guide with good features."""
    return SgrnaCandidate(
        guide_id="test_guide_1",
        guide_sequence="GAGTCTGCTGACAGAGCTCG",
        pam="AGG",
        pam_strand="+",
        chrom="chr1",
        start=100,
        end=119,
        strand="+",
        context_30nt="AAAAGAGTCTGCTGACAGAGCTCGAGGAAA",
        context_35nt="AAAAGAGAGTCTGCTGACAGAGCTCGAGGAAAAA",
        gc_content=0.55,
        nuclease=NucleaseType.SPCAS9,
        features=compute_guide_features("GAGTCTGCTGACAGAGCTCG"),
    )


@pytest.fixture
def poor_candidate():
    """A poorly designed guide with extreme GC and poly-T."""
    return SgrnaCandidate(
        guide_id="test_guide_2",
        guide_sequence="GCGCGCGCGCGCGCGCGCGC",
        pam="AGG",
        pam_strand="+",
        chrom="chr1",
        start=200,
        end=219,
        strand="+",
        context_30nt="",
        context_35nt="",
        gc_content=1.0,
        nuclease=NucleaseType.SPCAS9,
        features=compute_guide_features("GCGCGCGCGCGCGCGCGCGC"),
    )


class TestPredictEfficiency:
    def test_mouse_prediction(self, good_candidate):
        pred = predict_efficiency(good_candidate, species_key="mouse")
        assert 0.0 <= pred.efficiency_score <= 1.0
        assert pred.model_name == "RuleSet2-Transfer"
        assert pred.species == "mouse"

    def test_zebrafish_with_context(self, good_candidate):
        pred = predict_efficiency(good_candidate, species_key="zebrafish")
        assert pred.model_name == "CRISPRscan"
        assert pred.confidence > 0.5

    def test_zebrafish_without_context(self, poor_candidate):
        pred = predict_efficiency(poor_candidate, species_key="zebrafish")
        assert pred.model_name == "CRISPRscan"
        # Should warn about insufficient context
        assert any("Insufficient" in w for w in pred.warnings)

    def test_unsupported_species_raises(self, good_candidate):
        with pytest.raises(ValueError, match="No efficiency model"):
            predict_efficiency(good_candidate, species_key="alien")

    def test_poor_guide_lower_efficiency(self, good_candidate, poor_candidate):
        good = predict_efficiency(good_candidate, species_key="mouse")
        poor = predict_efficiency(poor_candidate, species_key="mouse")
        # The poor guide has 100% GC which should be penalized
        assert good.efficiency_score >= poor.efficiency_score

    def test_poly_t_warning(self):
        candidate = SgrnaCandidate(
            guide_id="test", guide_sequence="ATCGTTTTATCGATCGATCG",
            pam="AGG", pam_strand="+", chrom="chr1", start=1, end=20,
            strand="+", context_30nt="", context_35nt="",
            gc_content=0.35, nuclease=NucleaseType.SPCAS9,
            features=compute_guide_features("ATCGTTTTATCGATCGATCG"),
        )
        pred = predict_efficiency(candidate, species_key="mouse")
        assert any("Poly-T" in w for w in pred.warnings)


class TestPredictIndelOutcomes:
    def test_basic_prediction(self):
        outcome = predict_indel_outcomes("GAGTCTGCTGACAGAGCTCG", "mouse")
        assert 0.0 <= outcome.insertion_rate <= 1.0
        assert 0.0 <= outcome.deletion_rate <= 1.0
        assert 0.0 <= outcome.no_edit_rate <= 1.0
        assert outcome.insertion_rate + outcome.deletion_rate + outcome.no_edit_rate <= 1.01

    def test_most_likely_outcome(self):
        outcome = predict_indel_outcomes("GAGTCTGCTGACAGAGCTCG", "mouse")
        assert outcome.most_likely_outcome in ("deletion", "insertion", "no_edit")

    def test_predicted_indel_size_negative_for_deletion(self):
        outcome = predict_indel_outcomes("GAGTCTGCTGACAGAGCTCG", "mouse")
        if outcome.most_likely_outcome == "deletion":
            assert outcome.predicted_indel_size < 0


class TestComputeEvidenceScores:
    def test_basic_scores(self, good_candidate):
        eff = predict_efficiency(good_candidate, species_key="mouse")
        scores = compute_evidence_scores(eff)
        assert "on_target_uncertainty" in scores
        assert "off_target_evidence" in scores
        assert "network_impact_evidence" in scores
        assert "welfare_relevance" in scores
        for key in ("on_target_uncertainty", "off_target_evidence",
                     "network_impact_evidence", "welfare_relevance"):
            assert 0.0 <= scores[key] <= 1.0

    def test_with_offtarget_report(self, good_candidate):
        eff = predict_efficiency(good_candidate, species_key="mouse")
        ot_report = OffTargetReport(
            guide_sequence=good_candidate.guide_sequence,
            nuclease=NucleaseType.SPCAS9,
            total_sites_scanned=1,
            off_targets=[
                OffTargetSite(
                    chrom="chr2", start=50, end=69, strand="+",
                    off_target_sequence="GAGTCTGCTGACAGAGCTCA",
                    pam="AGG", mismatch_count=1,
                    mismatch_positions=(19,),
                    mismatch_pattern=MismatchPattern(positions=(19,), count=1),
                    score=0.7, risk_level="high",
                )
            ],
            max_mismatches_searched=3,
        )
        scores = compute_evidence_scores(eff, offtarget_report=ot_report)
        assert scores["off_target_evidence"] > 0.1

    def test_gene_essentiality_increases_network_impact(self, good_candidate):
        eff = predict_efficiency(good_candidate, species_key="mouse")
        low = compute_evidence_scores(eff, gene_essentiality=0.0)
        high = compute_evidence_scores(eff, gene_essentiality=1.0)
        assert high["network_impact_evidence"] > low["network_impact_evidence"]

    def test_phenotype_severity_increases_welfare(self, good_candidate):
        eff = predict_efficiency(good_candidate, species_key="mouse")
        low = compute_evidence_scores(eff, phenotype_severity=0.0)
        high = compute_evidence_scores(eff, phenotype_severity=1.0)
        assert high["welfare_relevance"] > low["welfare_relevance"]
