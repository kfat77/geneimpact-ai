import copy
from hashlib import sha256
from io import BytesIO
import json
import sys

import pytest
import xlwt

from geneimpact.cli import main
from geneimpact.housden import (
    HOUSDEN_SERVICE_URL,
    normalize_housden,
)


PROTOSPACER = "ATCTGACCTCCCGGCTAATT"
SERVICE_HEADERS = (
    "Identifier",
    "Input Sequence",
    "Input Sequence Length",
    "Analyzed Sequence",
    "Analyzed Sequence Length",
    "Number of Bases Analyzed",
    "Start Index",
    "End Index",
    "Score",
    "Comments",
    "U6 Terminator",
)


def housden_service_response(*, score="6.93243"):
    workbook = xlwt.Workbook()
    sheet = workbook.add_sheet("EfficiencyScoreList")
    for column, value in enumerate(SERVICE_HEADERS):
        sheet.write(0, column, value)
    values = (
        "1",
        PROTOSPACER,
        20,
        PROTOSPACER,
        20,
        20,
        1,
        20,
        score,
        "None",
        "None",
    )
    for column, value in enumerate(values):
        sheet.write(1, column, value)
    buffer = BytesIO()
    workbook.save(buffer)
    return buffer.getvalue()


SERVICE_RESPONSE = housden_service_response()


def housden_document(*, developmental_context="drosophila_s2r_plus_cell_culture"):
    return {
        "request": {
            "guide_id": "drsc-guide-01",
            "species_profile": "fruit_fly",
            "genome_build": "Release 6 plus ISO1 MT",
            "assembly_accession": "GCF_000001215.4",
            "sequence_source_strain_or_isolate": "ISO-1",
            "protospacer": PROTOSPACER,
            "nuclease": "SpCas9",
            "guide_expression": "u6_sgrna",
            "developmental_context": developmental_context,
        },
        "execution": {
            "source": "flyrnai_evaluate_crispr",
            "source_url": HOUSDEN_SERVICE_URL,
            "retrieved_at": "2026-07-29T09:30:00+00:00",
            "source_response_sha256": sha256(SERVICE_RESPONSE).hexdigest(),
        },
        "raw_output": {
            "housden_score": 6.93243,
        },
    }


def test_normalizes_official_housden_result_without_retaining_sequence():
    prediction = normalize_housden(
        housden_document(),
        source_response=SERVICE_RESPONSE,
        source_document_sha256="b" * 64,
    )

    assert prediction.predictor == "Housden"
    assert prediction.species_profile == "fruit_fly"
    assert prediction.housden_score == 6.93243
    assert prediction.score_semantics == "ranking_score_not_probability"
    assert prediction.sequence_sha256
    assert PROTOSPACER not in repr(prediction)
    assert prediction.source_document_sha256 == "b" * 64
    assert prediction.service_version_status == "live_service_unversioned"
    assert prediction.published_high_efficiency_threshold == 7.5
    assert prediction.current_recommended_threshold == 5.0


def test_rejects_housden_result_outside_s2r_cell_domain():
    with pytest.raises(ValueError, match="S2R"):
        normalize_housden(
            housden_document(developmental_context="fruit_fly_embryo"),
            source_response=SERVICE_RESPONSE,
        )


def test_rejects_unofficial_source_or_invalid_digest():
    document = housden_document()
    document["execution"]["source_url"] = "https://example.invalid/scorer"

    with pytest.raises(ValueError, match="source_url"):
        normalize_housden(document, source_response=SERVICE_RESPONSE)

    document = copy.deepcopy(housden_document())
    document["execution"]["source_response_sha256"] = "not-a-digest"

    with pytest.raises(ValueError, match="source_response_sha256"):
        normalize_housden(document, source_response=SERVICE_RESPONSE)


def test_rejects_score_not_present_in_retained_service_response():
    document = housden_document()
    document["raw_output"]["housden_score"] = 9.9

    with pytest.raises(ValueError, match="does not match"):
        normalize_housden(document, source_response=SERVICE_RESPONSE)


def test_housden_cli_writes_sequence_redacted_audit_record(
    tmp_path,
    monkeypatch,
):
    input_path = tmp_path / "housden-input.json"
    response_path = tmp_path / "official-response.xls"
    output_path = tmp_path / "housden-audit.json"
    input_path.write_text(json.dumps(housden_document()), encoding="utf-8")
    response_path.write_bytes(SERVICE_RESPONSE)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "geneimpact",
            "import-housden",
            "--input",
            str(input_path),
            "--source-response",
            str(response_path),
            "--output",
            str(output_path),
        ],
    )

    main()

    audit = json.loads(output_path.read_text(encoding="utf-8"))
    assert audit["predictor"] == "Housden"
    assert audit["housden_score"] == 6.93243
    assert PROTOSPACER not in output_path.read_text(encoding="utf-8")
