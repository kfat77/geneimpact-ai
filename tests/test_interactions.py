import pytest

from geneimpact.interactions import rank_interactions


def test_interactions_are_ranked_by_signal_and_evidence():
    results = rank_interactions(
        {"A": 1.0, "B": 0.81, "C": 0.25},
        {frozenset(("A", "B")): 0.8, frozenset(("A", "C")): 1.0},
    )

    assert results[0].genes == ("A", "B")
    assert results[0].priority == pytest.approx(0.72)
    assert results[-1].priority == 0.0


def test_invalid_scores_are_rejected():
    with pytest.raises(ValueError, match="between 0 and 1"):
        rank_interactions({"A": 1.1}, {})
