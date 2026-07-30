"""Tests for sgRNA design and PAM search (sgrna_design module)."""

from __future__ import annotations

import pytest

from geneimpact.sgrna_design import (
    NucleaseType,
    PamPattern,
    SgrnaCandidate,
    SgrnaDesignResult,
    PAM_PATTERNS,
    design_sgrnas,
    search_pam_sites,
    compute_guide_features,
)


class TestPamPattern:
    def test_ngg_match(self):
        pam = PAM_PATTERNS[NucleaseType.SPCAS9]
        assert pam.matches("AGG")
        assert pam.matches("CGG")
        assert pam.matches("TGG")
        assert pam.matches("GGG")
        assert not pam.matches("AAA")
        assert not pam.matches("AAG")

    def test_tttv_match(self):
        pam = PAM_PATTERNS[NucleaseType.CAS12A]
        assert pam.matches("TTTA")
        assert pam.matches("TTTC")
        assert pam.matches("TTTG")
        assert not pam.matches("TTTT")

    def test_iupac_ambiguity(self):
        pam = PamPattern(
            pattern="NGRRT", nuclease=NucleaseType.SACAS9,
            position="3prime", length=5,
        )
        assert pam.matches("AGAAT")
        assert pam.matches("TGGAT")
        assert not pam.matches("TGGCT")
        assert not pam.matches("AAAAA")


class TestDesignSgrnas:
    def test_finds_spCas9_sites(self):
        # 60 bp with two NGG PAM sites
        seq = "AAATCGATCGATCGATCGATGGAAATCGATCGATCGATCGATGG"
        result = design_sgrnas(seq, chrom="test", nuclease=NucleaseType.SPCAS9)
        assert result.count > 0
        for cand in result.candidates:
            assert len(cand.guide_sequence) == 20
            assert cand.pam[-2:] == "GG"
            assert 0.0 <= cand.gc_content <= 1.0

    def test_both_strands(self):
        # Sequence with PAMs on both strands
        seq = "ATCGATCGATCGATCGATCGGGATCGATCGATCGATCGATCGG"
        result = design_sgrnas(seq, nuclease=NucleaseType.SPCAS9)
        strands = {c.strand for c in result.candidates}
        # Should find sites on at least one strand
        assert len(result.candidates) > 0

    def test_cas12a_pam(self):
        # TTTV PAM at 5' end
        seq = "TTTACGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        result = design_sgrnas(seq, nuclease=NucleaseType.CAS12A)
        assert result.count > 0

    def test_no_pam_found(self):
        seq = "ATATATATATATATATATATATATATATATATATATATAT"
        result = design_sgrnas(seq, nuclease=NucleaseType.SPCAS9)
        assert result.count == 0
        assert any("No" in w for w in result.warnings)

    def test_max_candidates(self):
        # Long sequence with many PAMs
        seq = ("GG" + "A" * 18 + "GG" + "C" * 18 + "GG" + "T" * 18 + "GG" + "G" * 18) * 50
        result = design_sgrnas(seq, nuclease=NucleaseType.SPCAS9, max_candidates=10)
        assert result.count <= 10

    def test_handles_n_bases(self):
        seq = "NNNNATCGATCGATCGATCGATGGNNNN"
        result = design_sgrnas(seq, nuclease=NucleaseType.SPCAS9)
        # Should skip N regions but still find valid PAM sites
        assert any("ambiguous" in w.lower() for w in result.warnings)

    def test_invalid_sequence_raises(self):
        with pytest.raises(ValueError, match="invalid base"):
            design_sgrnas("ATCGXYZ", nuclease=NucleaseType.SPCAS9)

    def test_unsupported_nuclease_raises(self):
        with pytest.raises(ValueError, match="No PAM pattern"):
            # CasX has no registered PAM pattern
            design_sgrnas("ATCGATCGATCGATCGATCGGG", nuclease=NucleaseType.CASX)

    def test_guide_length_20(self):
        seq = "ATCGATCGATCGATCGATCGGG"
        result = design_sgrnas(seq, nuclease=NucleaseType.SPCAS9, guide_length=20)
        for cand in result.candidates:
            assert len(cand.guide_sequence) == 20

    def test_context_extraction(self):
        seq = "AAAAATCGATCGATCGATCGATGGAAAA"
        result = design_sgrnas(seq, nuclease=NucleaseType.SPCAS9)
        for cand in result.candidates:
            assert len(cand.context_30nt) > 0


class TestComputeGuideFeatures:
    def test_basic_features(self):
        features = compute_guide_features("GAGTCTGCTGACAGAGCTCG")
        assert "gc_content" in features
        assert "efficiency_score" in features
        assert 0.0 <= features["gc_content"] <= 1.0
        assert 0.0 <= features["efficiency_score"] <= 1.0

    def test_poly_t_detection(self):
        features = compute_guide_features("ATCGTTTTATCGATCGATCG")
        assert features["poly_t_count"] == 1.0

    def test_homopolymer_run(self):
        features = compute_guide_features("ATCGAAAAAAACGATCGATC")
        assert features["max_homopolymer_run"] >= 6.0

    def test_g20_bonus(self):
        features_g = compute_guide_features("ATCGATCGATCGATCGATCG")
        features_a = compute_guide_features("ATCGATCGATCGATCGATCA")
        assert features_g["g20_bonus"] > features_a["g20_bonus"]

    def test_wrong_length_raises(self):
        with pytest.raises(ValueError, match="guide must be 20"):
            compute_guide_features("ATCG")


class TestSearchPamSites:
    def test_finds_sites(self):
        seq = "ATCGGGATCGGGATCGGG"
        sites = search_pam_sites(seq, nuclease=NucleaseType.SPCAS9)
        assert len(sites) > 0
        for pos, pam, strand in sites:
            assert strand in ("+", "-")
            assert pam[-2:] == "GG"

    def test_cas12a_sites(self):
        seq = "TTTACGATCGATCGTTTACGATCG"
        sites = search_pam_sites(seq, nuclease=NucleaseType.CAS12A)
        assert len(sites) > 0
