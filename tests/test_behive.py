from dataclasses import asdict
from math import exp

import pytest

from geneimpact.behive import (
    BEHIVE_EFFICIENCY_COMMIT,
    integrate_behive_efficiency,
    normalize_behive_efficiency,
)


SEQUENCE = "TATCAGCGGGAATTCAAGCGCACCAGCCAGAGGTGTACCGTGGACGTGAG"


def _document(**request_overrides):
    request = {
        "sequence": SEQUENCE,
        "base_editor": "BE4",
        "cell_type": "mES",
        "model_commit": BEHIVE_EFFICIENCY_COMMIT,
    }
    request.update(request_overrides)
    return {
        "request": request,
        "raw_output": {"Predicted logit score": 0.4},
    }


def test_normalizes_unscaled_logit_without_retaining_sequence():
    prediction = normalize_behive_efficiency(_document())
    rendered = asdict(prediction)

    assert prediction.predicted_logit_score == 0.4
    assert prediction.calibrated_fraction is None
    assert prediction.biological_context == "mouse embryonic stem cells (mES)"
    assert SEQUENCE not in str(rendered)
    assert len(prediction.sequence_sha256) == 64


def test_recomputes_and_checks_calibrated_fraction():
    document = _document(calibration_mean=-0.2, calibration_std=1.5)
    expected = 1 / (1 + exp(-(0.4 * 1.5 - 0.2)))
    document["raw_output"][
        "Predicted fraction of sequenced reads with base editing activity"
    ] = expected

    prediction = normalize_behive_efficiency(document)

    assert prediction.calibrated_fraction == pytest.approx(expected)


def test_omits_fraction_without_declared_calibration():
    document = _document()
    document["raw_output"]["predicted_fraction"] = 0.8

    prediction = normalize_behive_efficiency(document)

    assert prediction.calibrated_fraction is None
    assert any("omitted" in warning for warning in prediction.warnings)


def test_rejects_fraction_inconsistent_with_calibration():
    document = _document(calibration_mean=-0.2, calibration_std=1.5)
    document["raw_output"]["predicted_fraction"] = 0.99

    with pytest.raises(ValueError, match="does not match"):
        normalize_behive_efficiency(document)


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"sequence": "A" * 50}, "NGG PAM"),
        ({"cell_type": "HEK293T"}, "only accepts"),
        ({"base_editor": "unknown"}, "not declared"),
        ({"model_commit": "abc"}, "40-character"),
    ],
)
def test_rejects_inputs_outside_verified_mouse_scope(overrides, message):
    with pytest.raises(ValueError, match=message):
        normalize_behive_efficiency(_document(**overrides))


def test_applicability_requires_mouse_base_editing_context():
    prediction = normalize_behive_efficiency(_document())

    matched = integrate_behive_efficiency((prediction,), "mouse", "base_editing")[0]
    mismatched = integrate_behive_efficiency((prediction,), "mouse", "knockout")[0]

    assert matched.applicability == "declared_match"
    assert mismatched.applicability == "out_of_scope"
