"""Tests for off-target detection and scoring (offtarget module)."""

from __future__ import annotations

import pytest

from geneimpact.offtarget import (
    OffTargetSite,
    OffTargetReport,
    MismatchPattern,
    find_offtargets,
    score_offtarget,
    compute_offtarget_risk,
)
from geneimpact.sgrna_design import NucleaseType


class TestMismatchPattern:
    def test_seed_mismatch_detection(self):
        pattern = MismatchPattern(positions=(15, 18), count=2)
        assert pattern.has_seed_mismatch
        assert pattern.seed_mismatch_count == 2

    def test_no_seed_mismatch(self):
        pattern = MismatchPattern(positions=(1, 3, 5), count=3)
        assert not pattern.has_seed_mismatch
        assert pattern.seed_mismatch_count == 0


class TestScoreOffTarget:
    def test_zero_mismatch_is_perfect(self):
        pattern = MismatchPattern(positions=(), count=0)
        guide = "ATCGATCGATCGATCGATCG"
        score = score_offtarget(guide, guide, pattern, NucleaseType.SPCAS9)
        assert score == 1.0

    def test_one_mismatch_higher_than_two(self):
        guide = "ATCGATCGATCGATCGATCG"
        target1 = "ATCGATCGATCGATCGATCA"  # 1 mismatch at pos 20
        target2 = "ATCGATCGATCGATCGATAC"  # 2 mismatches
        p1 = MismatchPattern(positions=(19,), count=1)
        p2 = MismatchPattern(positions=(18, 19), count=2)
        score1 = score_offtarget(guide, target1, p1, NucleaseType.SPCAS9)
        score2 = score_offtarget(guide, target2, p2, NucleaseType.SPCAS9)
        assert score1 > score2

    def test_seed_mismatch_lower_score(self):
        guide = "ATCGATCGATCGATCGATCG"
        # Mismatch in PAM-distal (position 0) vs PAM-proximal (position 19)
        target_distal = "TTCGATCGATCGATCGATCG"
        target_proximal = "ATCGATCGATCGATCGATCA"
        p_distal = MismatchPattern(positions=(0,), count=1)
        p_proximal = MismatchPattern(positions=(19,), count=1)
        score_distal = score_offtarget(guide, target_distal, p_distal, NucleaseType.SPCAS9)
        score_proximal = score_offtarget(guide, target_proximal, p_proximal, NucleaseType.SPCAS9)
        # Proximal mismatch should be more penalizing
        assert score_proximal < score_distal

    def test_score_in_range(self):
        guide = "ATCGATCGATCGATCGATCG"
        for pos in range(20):
            target = list(guide)
            target[pos] = "A" if target[pos] != "A" else "T"
            target = "".join(target)
            pattern = MismatchPattern(positions=(pos,), count=1)
            score = score_offtarget(guide, target, pattern, NucleaseType.SPCAS9)
            assert 0.0 <= score <= 1.0


class TestFindOffTargets:
    def test_finds_off_targets_in_reference(self):
        guide = "GAGTCTGCTGACAGAGCTCG"
        # Create a reference with a 1-mismatch site followed by a valid NGG PAM
        # guide:   GAGTCTGCTGACAGAGCTCG
        # target:  GAGTCTGCTGACAGAGCTAG  (1 mismatch at pos 19: C->A)
        # PAM:     AGG (matches NGG)
        reference = {
            "chr1": "ATCGGAGTCTGCTGACAGAGCTAGAGGCATCGATCG"
        }
        report = find_offtargets(
            guide_sequence=guide,
            reference_sequences=reference,
            nuclease=NucleaseType.SPCAS9,
            max_mismatches=3,
        )
        assert len(report.off_targets) > 0
        assert any(ot.mismatch_count == 1 for ot in report.off_targets)

    def test_excludes_exact_match(self):
        guide = "GAGTCTGCTGACAGAGCTCG"
        # Reference contains the exact on-target site
        reference = {
            "chr1": "ATC" + guide + "AGG" + "ATCGATCG"
        }
        report = find_offtargets(
            guide_sequence=guide,
            reference_sequences=reference,
            nuclease=NucleaseType.SPCAS9,
            max_mismatches=3,
        )
        for ot in report.off_targets:
            assert ot.mismatch_count > 0

    def test_max_mismatches_filter(self):
        guide = "GAGTCTGCTGACAGAGCTCG"
        reference = {
            "chr1": "ATCGAAAAAAAAAAAAAAAAAAAAGGATCG",  # 19 mismatches
        }
        report = find_offtargets(
            guide_sequence=guide,
            reference_sequences=reference,
            nuclease=NucleaseType.SPCAS9,
            max_mismatches=2,
        )
        assert all(ot.mismatch_count <= 2 for ot in report.off_targets)

    def test_specificity_score_range(self):
        guide = "GAGTCTGCTGACAGAGCTCG"
        reference = {"chr1": "ATCGGAGTCTGCTGACAGAGCTAGGGCATCG"}
        report = find_offtargets(
            guide_sequence=guide,
            reference_sequences=reference,
            nuclease=NucleaseType.SPCAS9,
            max_mismatches=3,
        )
        assert 0.0 <= report.specificity_score <= 1.0

    def test_no_reference_raises(self):
        guide = "GAGTCTGCTGACAGAGCTCG"
        with pytest.raises(ValueError):
            find_offtargets(guide_sequence=guide)

    def test_risk_level_classification(self):
        guide = "GAGTCTGCTGACAGAGCTCG"
        # 1-mismatch site with valid NGG PAM → should be high risk
        reference = {"chr1": "ATCGGAGTCTGCTGACAGAGCTAGAGGCATCGATCG"}
        report = find_offtargets(
            guide_sequence=guide,
            reference_sequences=reference,
            nuclease=NucleaseType.SPCAS9,
            max_mismatches=3,
        )
        risk_levels = {ot.risk_level for ot in report.off_targets}
        assert "high" in risk_levels or "moderate" in risk_levels

    def test_invalid_mismatch_count(self):
        guide = "GAGTCTGCTGACAGAGCTCG"
        with pytest.raises(ValueError):
            find_offtargets(
                guide_sequence=guide,
                reference_sequences={"chr1": "A" * 30},
                max_mismatches=10,
            )

    def test_invalid_guide_length(self):
        with pytest.raises(ValueError):
            find_offtargets(
                guide_sequence="ATCG",
                reference_sequences={"chr1": "A" * 30},
            )


class TestComputeOffTargetRisk:
    def test_no_off_targets(self):
        report = OffTargetReport(
            guide_sequence="GAGTCTGCTGACAGAGCTCG",
            nuclease=NucleaseType.SPCAS9,
            total_sites_scanned=0,
            off_targets=[],
            max_mismatches_searched=4,
        )
        risk = compute_offtarget_risk(report)
        assert risk["off_target_evidence"] == 0.0
        assert risk["specificity_score"] == 1.0

    def test_high_risk_increases_concern(self):
        from geneimpact.offtarget import OffTargetSite
        guide = "GAGTCTGCTGACAGAGCTCG"
        ot = OffTargetSite(
            chrom="chr1", start=1, end=20, strand="+",
            off_target_sequence="GAGTCTGCTGACAGAGCTCA",
            pam="AGG", mismatch_count=1,
            mismatch_positions=(19,),
            mismatch_pattern=MismatchPattern(positions=(19,), count=1),
            score=0.7, risk_level="high",
        )
        report = OffTargetReport(
            guide_sequence=guide,
            nuclease=NucleaseType.SPCAS9,
            total_sites_scanned=1,
            off_targets=[ot],
            max_mismatches_searched=3,
        )
        risk = compute_offtarget_risk(report)
        assert risk["off_target_evidence"] > 0.3
