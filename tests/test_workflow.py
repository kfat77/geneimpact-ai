import pytest

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
    assert report["predictor_outputs"] == []
    assert "does not establish safety" in report["report_notice"]


def test_request_rejects_missing_context():
    with pytest.raises(ValueError, match="study_context"):
        assess_request({"evidence": {}})
