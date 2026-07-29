from dataclasses import asdict

import pytest

from geneimpact.behive_bystander import (
    BEHIVE_BYSTANDER_COMMIT,
    normalize_behive_bystander,
)


SEQUENCE = "TATCAGCGGGAATTCAAGCGCACCAGCCAGAGGTGTACCGTGGACGTGAG"


def _document():
    edited_once = list(SEQUENCE)
    edited_once[22] = "T"
    edited_twice = list(SEQUENCE)
    edited_twice[22] = "T"
    edited_twice[26] = "T"
    return {
        "request": {
            "sequence": SEQUENCE,
            "base_editor": "BE4",
            "cell_type": "mES",
            "model_commit": BEHIVE_BYSTANDER_COMMIT,
        },
        "raw_output": {
            "stats": {
                "Total predicted probability": 0.95,
                "50-nt target sequence": SEQUENCE,
                "Assumed protospacer sequence": SEQUENCE[20:40],
                "Celltype": "mES",
                "Base editor": "BE4",
            },
            "outcomes": [
                {"Genotype": "".join(edited_once), "Predicted frequency": 0.7},
                {"Genotype": "".join(edited_twice), "Predicted frequency": 0.25},
            ],
        },
    }


def test_normalizes_conditional_outcomes_without_raw_sequences():
    prediction = normalize_behive_bystander(_document())
    rendered = asdict(prediction)

    assert prediction.total_predicted_probability == 0.95
    assert prediction.unreported_probability_mass == pytest.approx(0.05)
    assert prediction.outcome_count == 2
    assert prediction.top_outcomes[0].predicted_frequency == 0.7
    assert prediction.top_outcomes[0].edits == ("C3>T",)
    assert SEQUENCE not in str(rendered)
    assert len(prediction.raw_output_sha256) == 64


def test_rejects_frequency_sum_that_disagrees_with_stats():
    document = _document()
    document["raw_output"]["stats"]["Total predicted probability"] = 0.8

    with pytest.raises(ValueError, match="do not sum"):
        normalize_behive_bystander(document)


def test_rejects_duplicate_genotypes():
    document = _document()
    document["raw_output"]["outcomes"][1]["Genotype"] = document["raw_output"]["outcomes"][0][
        "Genotype"
    ]

    with pytest.raises(ValueError, match="duplicates"):
        normalize_behive_bystander(document)


def test_rejects_unedited_genotype():
    document = _document()
    document["raw_output"]["stats"]["Total predicted probability"] = 0.25
    document["raw_output"]["outcomes"] = [
        {"Genotype": SEQUENCE, "Predicted frequency": 0.25}
    ]

    with pytest.raises(ValueError, match="unedited"):
        normalize_behive_bystander(document)


@pytest.mark.parametrize(
    ("section", "field", "value", "message"),
    [
        ("request", "cell_type", "HEK293T", "only accepts"),
        ("request", "base_editor", "unknown", "not declared"),
        ("request", "model_commit", "bad", "verified only"),
        ("stats", "Base editor", "ABE", "does not match"),
        ("stats", "Assumed protospacer sequence", "A" * 20, "does not match"),
    ],
)
def test_rejects_out_of_scope_or_inconsistent_metadata(section, field, value, message):
    document = _document()
    target = document["request"] if section == "request" else document["raw_output"]["stats"]
    target[field] = value

    with pytest.raises(ValueError, match=message):
        normalize_behive_bystander(document)
