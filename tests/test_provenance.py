import pytest

from geneimpact.edit_assessment import EditEvidence, assess_edit
from geneimpact.provenance import StudyContext, create_record


def test_record_captures_applicability_boundary():
    record = create_record(
        StudyContext("mouse", "C57BL/6", "GRCm39", "knockout", "evidence-2026-07"),
        assess_edit(EditEvidence(0.1, 0.2, 0.3, 0.4)),
        "0.2.0",
    )

    assert record.context.species == "mouse"


def test_record_rejects_missing_provenance():
    with pytest.raises(ValueError, match="genome_build"):
        create_record(
            StudyContext("mouse", "C57BL/6", "", "knockout", "snapshot"),
            assess_edit(EditEvidence(0.1, 0.2, 0.3, 0.4)),
            "0.2.0",
        )
