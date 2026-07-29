import pytest

from geneimpact.behive import BEHIVE_EFFICIENCY_COMMIT
from geneimpact.behive_bystander import BEHIVE_BYSTANDER_COMMIT
from geneimpact.workflow import assess_request


def test_request_generates_auditable_report():
    report = assess_request(
        {
            "study_context": {
                "species": "mouse",
                "strain_or_breed": "C57BL/6",
                "genome_build": "GRCm39",
                "edit_class": "knockout",
                "evidence_snapshot": "snapshot-1",
            },
            "evidence": {
                "on_target_uncertainty": 0.1,
                "off_target_evidence": 0.2,
                "network_impact_evidence": 0.3,
                "welfare_relevance": 0.8,
            },
        }
    )

    assert report["assessment"]["tier"] == "high_concern_review"
    assert report["context"]["species"] == "mouse"
    assert report["species_validation"]["supported"]
    assert report["predictor_outputs"] == []
    assert report["model_predictions"] == []
    assert "does not establish safety" in report["report_notice"]


def test_request_rejects_missing_context():
    with pytest.raises(ValueError, match="study_context"):
        assess_request({"evidence": {}})


def test_request_integrates_behive_without_treating_efficiency_as_concern():
    request = {
        "study_context": {
            "species": "mouse",
            "strain_or_breed": "C57BL/6",
            "genome_build": "GRCm39",
            "edit_class": "base_editing",
            "evidence_snapshot": "snapshot-1",
        },
        "evidence": {
            "on_target_uncertainty": 0.1,
            "off_target_evidence": 0.2,
            "network_impact_evidence": 0.3,
            "welfare_relevance": 0.4,
        },
        "behive_efficiency_outputs": [
            {
                "request": {
                    "sequence": "TATCAGCGGGAATTCAAGCGCACCAGCCAGAGGTGTACCGTGGACGTGAG",
                    "base_editor": "BE4",
                    "cell_type": "mES",
                    "model_commit": BEHIVE_EFFICIENCY_COMMIT,
                },
                "raw_output": {"Predicted logit score": 0.4},
            }
        ],
    }

    report = assess_request(request)

    prediction = report["model_predictions"][0]
    assert prediction["task"] == "base_editing_efficiency"
    assert prediction["applicability"] == "declared_match"
    assert "concern_score" not in prediction


def test_request_integrates_bystander_distribution_as_task_specific_output():
    sequence = "TATCAGCGGGAATTCAAGCGCACCAGCCAGAGGTGTACCGTGGACGTGAG"
    edited = "TATCAGCGGGAATTCAAGCGCATCAGCCAGAGGTGTACCGTGGACGTGAG"
    request = {
        "study_context": {
            "species": "mouse",
            "strain_or_breed": "C57BL/6",
            "genome_build": "GRCm39",
            "edit_class": "base_editing",
            "evidence_snapshot": "snapshot-1",
        },
        "evidence": {
            "on_target_uncertainty": 0.1,
            "off_target_evidence": 0.2,
            "network_impact_evidence": 0.3,
            "welfare_relevance": 0.4,
        },
        "behive_bystander_outputs": [
            {
                "request": {
                    "sequence": sequence,
                    "base_editor": "BE4",
                    "cell_type": "mES",
                    "model_commit": BEHIVE_BYSTANDER_COMMIT,
                },
                "raw_output": {
                    "stats": {
                        "Total predicted probability": 1.0,
                        "50-nt target sequence": sequence,
                        "Assumed protospacer sequence": sequence[20:40],
                        "Celltype": "mES",
                        "Base editor": "BE4",
                    },
                    "outcomes": [
                        {"Genotype": edited, "Predicted frequency": 1.0}
                    ],
                },
            }
        ],
    }

    report = assess_request(request)

    prediction = report["model_predictions"][0]
    assert prediction["task"] == "base_editing_bystander_outcomes"
    assert prediction["applicability"] == "declared_match"
    assert prediction["top_outcomes"][0]["edits"] == ("C3>T",)
    assert sequence not in str(prediction)
