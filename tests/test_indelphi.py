import copy

import pytest

from geneimpact.indelphi import (
    INDELPHI_COMMIT,
    INDELPHI_MODEL_ARTIFACT_FAMILY,
    normalize_indelphi,
)
from geneimpact.indelphi_validation import (
    INDELPHI_MOUSE_SUPPLEMENT_SHA256,
    indelphi_mouse_transfer_evidence,
)


SEQUENCE = "ACGT" * 30


def indelphi_document(*, target_context="coding_sequence"):
    return {
        "request": {
            "target_id": "Tyr-guide-01",
            "species_profile": "mouse",
            "genome_build": "GRCm39",
            "assembly_accession": "GCF_000001635.27",
            "sequence_source_strain_or_isolate": "C57BL/6J",
            "sequence": SEQUENCE,
            "cutsite": 60,
            "cell_type": "mESC",
            "model_commit": INDELPHI_COMMIT,
            "model_artifact_family": INDELPHI_MODEL_ARTIFACT_FAMILY,
            "nuclease": "SpCas9",
            "delivery_context": "spcas9_ribonucleoprotein",
            "developmental_context": "one_cell_embryo",
            "target_context": target_context,
        },
        "execution": {
            "repository_commit": INDELPHI_COMMIT,
            "model_artifact_family": INDELPHI_MODEL_ARTIFACT_FAMILY,
            "model_bundle_sha256": "a" * 64,
            "cell_type": "mESC",
            "python_version": "3.6.15",
            "sklearn_version": "0.18.1",
        },
        "raw_output": {
            "stats": {
                "Phi": 8.0,
                "Precision": 0.06276943678387026,
                "1-bp ins frequency": 20.0,
                "MH del frequency": 50.0,
                "MHless del frequency": 30.0,
                "Frameshift frequency": 70.0,
                "Frame +0 frequency": 30.0,
                "Frame +1 frequency": 20.0,
                "Frame +2 frequency": 50.0,
                "Highest outcome frequency": 50.0,
                "Highest del frequency": 50.0,
                "Highest ins frequency": 20.0,
                "Expected indel length": 3.1,
                "Reference sequence": SEQUENCE,
                "Cutsite": 60,
                "Celltype": "mESC",
            },
            "outcomes": [
                {
                    "Category": "ins",
                    "Length": 1,
                    "Inserted Bases": "A",
                    "Predicted frequency": 20.0,
                },
                {
                    "Category": "del",
                    "Length": 3,
                    "Genotype position": "e",
                    "Predicted frequency": 30.0,
                },
                {
                    "Category": "del",
                    "Length": 4,
                    "Genotype position": 2,
                    "Predicted frequency": 50.0,
                },
            ],
        },
    }


def test_normalizes_version_locked_mes_prediction_without_retaining_sequence():
    prediction = normalize_indelphi(
        indelphi_document(),
        source_document_sha256="b" * 64,
    )

    assert prediction.predictor == "inDelphi"
    assert prediction.species_profile == "mouse"
    assert prediction.frameshift_frequency_percent == 70.0
    assert prediction.total_predicted_frequency_percent == 100.0
    assert prediction.reported_outcome_count == 3
    assert prediction.top_outcomes[0].predicted_frequency_percent == 50.0
    assert prediction.sequence_sha256
    assert SEQUENCE not in repr(prediction)
    assert prediction.source_document_sha256 == "b" * 64
    assert prediction.external_validation.overall_pearson_r == pytest.approx(
        0.6375914615783601
    )


def test_hides_frameshift_semantics_for_noncoding_or_unknown_target():
    prediction = normalize_indelphi(
        indelphi_document(target_context="noncoding_or_unknown")
    )

    assert prediction.frameshift_frequency_percent is None
    assert prediction.frame_0_frequency_percent is None


def test_rejects_tampered_summary_statistics():
    document = indelphi_document()
    document["raw_output"]["stats"]["Frameshift frequency"] = 90

    with pytest.raises(ValueError, match="Frameshift"):
        normalize_indelphi(document)


def test_rejects_truncated_sequence_context():
    document = indelphi_document()
    document["request"]["sequence"] = "ACGT" * 20
    document["request"]["cutsite"] = 40
    document["raw_output"]["stats"]["Reference sequence"] = "ACGT" * 20
    document["raw_output"]["stats"]["Cutsite"] = 40

    with pytest.raises(ValueError, match="sequence"):
        normalize_indelphi(document)


def test_rejects_duplicate_outcomes():
    document = indelphi_document()
    document["raw_output"]["outcomes"].append(
        copy.deepcopy(document["raw_output"]["outcomes"][0])
    )
    document["raw_output"]["outcomes"][0]["Predicted frequency"] = 10
    document["raw_output"]["outcomes"][-1]["Predicted frequency"] = 10

    with pytest.raises(ValueError, match="duplicate"):
        normalize_indelphi(document)


def test_external_validation_record_binds_the_source_without_copying_rows():
    evidence = indelphi_mouse_transfer_evidence()

    assert evidence.source_supplement_sha256 == INDELPHI_MOUSE_SUPPLEMENT_SHA256
    assert evidence.guide_count == 14
    assert evidence.compared_outcome_pairs == 1182
    assert (
        evidence.very_strong_guide_count
        + evidence.strong_guide_count
        + evidence.moderate_guide_count
        + evidence.weak_guide_count
        == evidence.guide_count
    )
    assert not evidence.source_rows_redistributed
