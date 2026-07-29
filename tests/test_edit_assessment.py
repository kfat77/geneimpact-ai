import pytest

from geneimpact.edit_assessment import EditEvidence, ReviewTier, assess_edit


def test_highest_concern_signal_controls_review_tier():
    result = assess_edit(EditEvidence(0.1, 0.2, 0.3, 0.8))

    assert result.tier is ReviewTier.HIGH_CONCERN_REVIEW
    assert result.concern_score == 0.8
    assert result.rationale == ("welfare relevance",)


def test_invalid_evidence_is_rejected():
    with pytest.raises(ValueError, match="between 0 and 1"):
        assess_edit(EditEvidence(-0.1, 0.2, 0.3, 0.4))
