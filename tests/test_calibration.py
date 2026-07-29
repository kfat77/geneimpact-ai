import pytest

from geneimpact.calibration import brier_score, expected_calibration_error


def test_perfect_predictions_have_zero_calibration_error():
    assert brier_score([0.0, 1.0], [0, 1]) == 0.0
    assert expected_calibration_error([0.0, 1.0], [0, 1]) == 0.0


def test_calibration_rejects_invalid_inputs():
    with pytest.raises(ValueError, match="equal length"):
        brier_score([0.1], [])
