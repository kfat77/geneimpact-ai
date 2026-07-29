"""Dependency-free calibration checks for held-out animal-edit outcome data."""

from __future__ import annotations

from collections.abc import Sequence


def brier_score(predictions: Sequence[float], outcomes: Sequence[int]) -> float:
    """Return mean squared probability error for binary held-out outcomes."""
    _validate_pairs(predictions, outcomes)
    return sum((prediction - outcome) ** 2 for prediction, outcome in zip(predictions, outcomes)) / len(predictions)


def expected_calibration_error(
    predictions: Sequence[float], outcomes: Sequence[int], bins: int = 10
) -> float:
    """Return an equal-width-bin calibration error for reported probabilities."""
    _validate_pairs(predictions, outcomes)
    if bins < 2:
        raise ValueError("bins must be at least 2.")

    weighted_error = 0.0
    for index in range(bins):
        lower, upper = index / bins, (index + 1) / bins
        members = [
            (prediction, outcome)
            for prediction, outcome in zip(predictions, outcomes)
            if lower <= prediction < upper or (index == bins - 1 and prediction == upper)
        ]
        if members:
            average_prediction = sum(item[0] for item in members) / len(members)
            observed_frequency = sum(item[1] for item in members) / len(members)
            weighted_error += len(members) / len(predictions) * abs(average_prediction - observed_frequency)
    return weighted_error


def _validate_pairs(predictions: Sequence[float], outcomes: Sequence[int]) -> None:
    if not predictions or len(predictions) != len(outcomes):
        raise ValueError("predictions and outcomes must be non-empty and have equal length.")
    if any(not 0 <= prediction <= 1 for prediction in predictions):
        raise ValueError("predictions must be between 0 and 1.")
    if any(outcome not in (0, 1) for outcome in outcomes):
        raise ValueError("outcomes must contain only 0 or 1.")
