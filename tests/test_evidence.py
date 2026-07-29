from geneimpact.evidence import EvidenceLevel, permitted_wording


def test_evidence_level_sets_claim_language():
    assert permitted_wording(EvidenceLevel.EXPLORATORY) == "candidate concern"
    assert "strongly supported" in permitted_wording(EvidenceLevel.CAUSAL_SUPPORT)
