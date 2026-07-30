"""Tests for the end-to-end prediction pipeline (pipeline module)."""

from __future__ import annotations

import json
import pytest

from geneimpact.pipeline import (
    PipelineConfig,
    GuideResult,
    PipelineReport,
    run_pipeline,
    run_pipeline_from_fasta,
)
from geneimpact.provenance import StudyContext
from geneimpact.sgrna_design import NucleaseType


@pytest.fixture
def mouse_context():
    return StudyContext(
        species="mouse",
        strain_or_breed="C57BL/6J",
        genome_build="GRCm39",
        edit_class="knockout",
        evidence_snapshot="test_snapshot_v1",
    )


@pytest.fixture
def test_sequence():
    """A 200 bp sequence with multiple SpCas9 PAM sites."""
    return (
        "ATCGATCGATCGATCGATCGGGATCGATCGATCGATCGATCGG"
        "ATCGATCGATCGATCGATCGGGATCGATCGATCGATCGATCGG"
        "GAGTCTGCTGACAGAGCTCGGGATCGATCGATCGATCGATCGG"
        "ATCGATCGATCGATCGATCGGGGAGTCTGCTGACAGAGCTCGG"
        "ATCGATCGATCGATCGATCGGGATCGATCGATCGATCGATCGG"
    )


@pytest.fixture
def small_reference():
    """Small reference genome for off-target search."""
    return {
        "chr1": (
            "ATCGATCGATCGATCGATCGGGATCGATCGATCGATCGATCGG"
            "GAGTCTGCTGACAGAGCTAGGGCATCGATCGATCGATCGATCG"
        ),
        "chr2": (
            "TTTTAAAATTTTAAAATTTTAAAACCCGGGATCGATCGATCGG"
            "GAGTCTGCTGACAGAGCTAGGGCATCGATCGATCGATCGATCG"
        ),
    }


class TestRunPipeline:
    def test_basic_pipeline_run(self, test_sequence, mouse_context, small_reference):
        config = PipelineConfig(
            species="mouse",
            nuclease=NucleaseType.SPCAS9,
            max_candidates=10,
            max_offtargets=3,
            top_k=5,
        )
        report = run_pipeline(
            sequence=test_sequence,
            config=config,
            study_context=mouse_context,
            reference_sequences=small_reference,
        )
        assert isinstance(report, PipelineReport)
        assert len(report.guides) > 0
        for guide in report.guides:
            assert isinstance(guide, GuideResult)
            assert guide.rank > 0
            assert guide.assessment is not None
            assert 0.0 <= guide.efficiency.efficiency_score <= 1.0

    def test_pipeline_to_json(self, test_sequence, mouse_context, small_reference):
        config = PipelineConfig(species="mouse", max_candidates=5, top_k=3)
        report = run_pipeline(
            sequence=test_sequence,
            config=config,
            study_context=mouse_context,
            reference_sequences=small_reference,
        )
        json_str = report.to_json()
        data = json.loads(json_str)
        assert "guides" in data
        assert "study_context" in data
        assert "report_notice" in data
        assert len(data["guides"]) <= 3

    def test_pipeline_to_json_file(self, test_sequence, mouse_context, small_reference, tmp_path):
        config = PipelineConfig(species="mouse", max_candidates=5, top_k=3)
        report = run_pipeline(
            sequence=test_sequence,
            config=config,
            study_context=mouse_context,
            reference_sequences=small_reference,
        )
        output_path = tmp_path / "report.json"
        report.to_json(output_path)
        assert output_path.exists()
        data = json.loads(output_path.read_text(encoding="utf-8"))
        assert "guides" in data

    def test_invalid_species_raises(self, test_sequence):
        config = PipelineConfig(species="mouse")
        bad_context = StudyContext(
            species="alien",
            strain_or_breed="unknown",
            genome_build="unknown",
            edit_class="knockout",
            evidence_snapshot="test",
        )
        with pytest.raises(ValueError, match="Invalid study context"):
            run_pipeline(
                sequence=test_sequence,
                config=config,
                study_context=bad_context,
            )

    def test_no_pam_sites(self, mouse_context):
        # ATAT repeat has no NGG PAMs
        seq = "ATATATATATATATATATATATATATATATATATATATAT"
        config = PipelineConfig(species="mouse", max_candidates=5)
        report = run_pipeline(
            sequence=seq,
            config=config,
            study_context=mouse_context,
            reference_sequences={"chr1": seq},
        )
        assert len(report.guides) == 0
        assert any("No sgRNA" in w for w in report.warnings)

    def test_guide_ranking(self, test_sequence, mouse_context, small_reference):
        config = PipelineConfig(
            species="mouse", max_candidates=10, top_k=5,
            min_efficiency=0.0, min_specificity=0.0,
        )
        report = run_pipeline(
            sequence=test_sequence,
            config=config,
            study_context=mouse_context,
            reference_sequences=small_reference,
        )
        # Guides should be ranked
        ranks = [g.rank for g in report.guides]
        assert ranks == sorted(ranks)
        assert ranks[0] == 1

    def test_evidence_scores_present(self, test_sequence, mouse_context, small_reference):
        config = PipelineConfig(species="mouse", max_candidates=5, top_k=3)
        report = run_pipeline(
            sequence=test_sequence,
            config=config,
            study_context=mouse_context,
            reference_sequences=small_reference,
        )
        for guide in report.guides:
            assert "on_target_uncertainty" in guide.evidence_scores
            assert "off_target_evidence" in guide.evidence_scores
            assert "network_impact_evidence" in guide.evidence_scores
            assert "welfare_relevance" in guide.evidence_scores

    def test_recommendation_generated(self, test_sequence, mouse_context, small_reference):
        config = PipelineConfig(species="mouse", max_candidates=5, top_k=3)
        report = run_pipeline(
            sequence=test_sequence,
            config=config,
            study_context=mouse_context,
            reference_sequences=small_reference,
        )
        for guide in report.guides:
            assert len(guide.recommendation) > 0

    def test_offtarget_analysis_included(self, test_sequence, mouse_context, small_reference):
        config = PipelineConfig(species="mouse", max_candidates=5, top_k=3)
        report = run_pipeline(
            sequence=test_sequence,
            config=config,
            study_context=mouse_context,
            reference_sequences=small_reference,
        )
        for guide in report.guides:
            if guide.offtarget_report is not None:
                assert guide.offtarget_report.specificity_score >= 0.0


class TestRunPipelineFromFasta:
    def test_fasta_pipeline(self, tmp_path, mouse_context):
        # Create a test FASTA file
        fasta_path = tmp_path / "target.fa"
        seq = (
            "ATCGATCGATCGATCGATCGGGATCGATCGATCGATCGATCGG"
            "GAGTCTGCTGACAGAGCTCGGGATCGATCGATCGATCGATCGG"
        )
        fasta_path.write_text(f">chr1\n{seq}\n", encoding="utf-8")

        config = PipelineConfig(species="mouse", max_candidates=5, top_k=3)
        report = run_pipeline_from_fasta(
            fasta_path=fasta_path,
            target_chrom="chr1",
            target_start=1,
            target_end=len(seq),
            config=config,
            study_context=mouse_context,
        )
        assert isinstance(report, PipelineReport)
        assert len(report.guides) > 0
