import pytest

from geneimpact.predictors import Applicability, PredictionTask, PredictorOutput, integrate_outputs


def test_output_is_only_applicable_inside_declared_scope():
    output = PredictorOutput(
        "repair-model", "1.0", PredictionTask.REPAIR_OUTCOME, 0.5, 0.8,
        ("mouse",), ("knockout",), "validation-record"
    )

    applicable, out_of_scope = integrate_outputs(
        (output, output), "mouse", "knockout"
    )[0], integrate_outputs((output,), "rat", "knockout")[0]

    assert applicable.applicability is Applicability.DECLARED_MATCH
    assert out_of_scope.applicability is Applicability.OUT_OF_SCOPE


def test_predictor_requires_declared_scope():
    with pytest.raises(ValueError, match="supported species"):
        integrate_outputs(
            (PredictorOutput("x", "1", PredictionTask.OFF_TARGET, 0.1, 0.1, (), ("knockout",), "ref"),),
            "mouse", "knockout"
        )
