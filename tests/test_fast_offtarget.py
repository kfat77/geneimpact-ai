"""Tests for the fast off-target search (k-mer seed-and-extend)."""

import pytest
from geneimpact.fast_offtarget import (
    fast_find_offtargets,
    build_seed_index,
    SeedIndex,
    DEFAULT_SEED_LENGTH,
)
from geneimpact.sgrna_design import NucleaseType
from geneimpact.offtarget import find_offtargets


class TestSeedIndex:
    def test_build_index(self):
        ref = {"chr1": "ATCGATCGATCGATCGATCGGGATCGATCGATCGATCGATCGG"}
        idx = build_seed_index(ref, seed_length=10)
        assert idx.kmer_length == 10
        assert idx.total_positions > 0
        assert idx.reference_length == len(ref["chr1"])

    def test_lookup(self):
        ref = {"chr1": "ATCGATCGATCGATCGATCGGGATCGATCGATCGATCGATCGG"}
        idx = build_seed_index(ref, seed_length=10)
        # The sequence contains "ATCGATCGAT" which should be in the index
        positions = idx.lookup("ATCGATCGAT")
        assert len(positions) > 0

    def test_lookup_missing(self):
        ref = {"chr1": "ATCGATCGATCGATCGATCGGG"}
        idx = build_seed_index(ref, seed_length=10)
        positions = idx.lookup("GGGGGGGGGG")
        assert positions == []

    def test_handles_n_bases(self):
        ref = {"chr1": "ATCGNATCGATCGATCGATCGGG"}
        idx = build_seed_index(ref, seed_length=10)
        # K-mers with N should be skipped
        assert idx.total_positions > 0


class TestFastFindOffTargets:
    def test_finds_off_targets(self):
        guide = "GAGTCTGCTGACAGAGCTCG"
        # Reference with a 1-mismatch off-target site (PAM must be NGG)
        reference = {
            "chr1": "ATCGGAGTCTGCTGACAGAGCTAGGGGATCG" + "A" * 30
        }
        report = fast_find_offtargets(
            guide_sequence=guide,
            reference_sequences=reference,
            nuclease=NucleaseType.SPCAS9,
            max_mismatches=3,
        )
        assert len(report.off_targets) > 0
        assert any(ot.mismatch_count == 1 for ot in report.off_targets)

    def test_consistent_with_brute_force(self):
        """Fast algorithm should find the same off-targets as brute-force."""
        guide = "GAGTCTGCTGACAGAGCTCG"
        reference = {
            "chr1": "ATCGGAGTCTGCTGACAGAGCTAGGGGATCG" + "A" * 50,
            "chr2": "GAGTCTGCTGACAGAGCTCGGG" + "T" * 50,  # On-target + near matches
        }

        fast_report = fast_find_offtargets(
            guide_sequence=guide,
            reference_sequences=reference,
            nuclease=NucleaseType.SPCAS9,
            max_mismatches=3,
        )
        brute_report = find_offtargets(
            guide_sequence=guide,
            reference_sequences=reference,
            nuclease=NucleaseType.SPCAS9,
            max_mismatches=3,
        )

        # Both should find off-target sites
        assert len(fast_report.off_targets) > 0
        assert len(brute_report.off_targets) > 0
        # The fast algorithm should find a subset (it uses seed matching)
        # but should find the same high-confidence sites
        fast_seqs = {ot.off_target_sequence for ot in fast_report.off_targets}
        brute_seqs = {ot.off_target_sequence for ot in brute_report.off_targets}
        # There should be overlap
        assert fast_seqs & brute_seqs

    def test_no_off_targets_for_unique_guide(self):
        guide = "GAGTCTGCTGACAGAGCTCG"
        reference = {"chr1": "A" * 100}
        report = fast_find_offtargets(
            guide_sequence=guide,
            reference_sequences=reference,
            nuclease=NucleaseType.SPCAS9,
            max_mismatches=3,
        )
        assert len(report.off_targets) == 0

    def test_requires_reference_or_index(self):
        with pytest.raises(ValueError, match="reference_sequences or seed_index"):
            fast_find_offtargets("GAGTCTGCTGACAGAGCTCG")

    def test_invalid_guide_length(self):
        with pytest.raises(ValueError, match="20 nt"):
            fast_find_offtargets(
                "ATCG",
                reference_sequences={"chr1": "ATCGATCG"},
            )

    def test_returns_offtarget_report(self):
        guide = "GAGTCTGCTGACAGAGCTCG"
        reference = {"chr1": "GAGTCTGCTGACAGAGCTAGGG" + "A" * 30}
        report = fast_find_offtargets(
            guide_sequence=guide,
            reference_sequences=reference,
            nuclease=NucleaseType.SPCAS9,
            max_mismatches=3,
        )
        assert hasattr(report, "off_targets")
        assert hasattr(report, "specificity_score")
        assert hasattr(report, "high_risk_count")
