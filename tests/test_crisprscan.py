from dataclasses import asdict
import json

import pytest

from geneimpact.crisprscan import (
    CRISPRSCAN_COEFFICIENTS_SHA256,
    CRISPRSCORE_COMMIT,
    score_crisprscan,
)


EXAMPLE_CONTEXT = "ACCTGGATCGATGCTGATGCTAGATAAGGTTGAGC"


def _request(**overrides):
    request = {
        "species_profile": "zebrafish",
        "genome_build": "GRCz12tu",
        "assembly_accession": "GCF_049306965.2",
        "reference_strain_or_isolate": "Tuebingen",
        "nuclease": "SpCas9",
        "guide_expression": "t7_in_vitro_transcription",
        "developmental_context": "zebrafish_embryo",
        "guides": [
            {
                "guide_id": "zf-guide-001",
                "context_35nt": EXAMPLE_CONTEXT,
            }
        ],
    }
    request.update(overrides)
    return request


def test_matches_version_locked_crisprscore_oracle():
    report = score_crisprscan(_request())
    prediction = report.predictions[0]

    assert prediction.score == pytest.approx(0.5691911213797, abs=1e-12)
    assert prediction.published_threshold_label == (
        "above_published_efficient_threshold"
    )
    assert report.implementation_commit == CRISPRSCORE_COMMIT
    assert report.coefficients_sha256 == CRISPRSCAN_COEFFICIENTS_SHA256


def test_output_hashes_sequence_and_does_not_repeat_raw_context():
    report = score_crisprscan(_request())
    rendered = json.dumps(asdict(report))

    assert len(report.predictions[0].context_sha256) == 64
    assert EXAMPLE_CONTEXT not in rendered
    assert len(report.request_sha256) == 64


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"species_profile": "mouse"}, "only"),
        ({"genome_build": "GRCz11"}, "genome_build"),
        ({"assembly_accession": "GCF_old"}, "assembly_accession"),
        ({"reference_strain_or_isolate": "AB"}, "reference_strain"),
        ({"nuclease": "Cas12a"}, "nuclease"),
        ({"guide_expression": "u6"}, "guide_expression"),
        ({"developmental_context": "adult_liver"}, "developmental_context"),
    ],
)
def test_rejects_out_of_domain_requests(overrides, message):
    with pytest.raises(ValueError, match=message):
        score_crisprscan(_request(**overrides))


@pytest.mark.parametrize(
    "context",
    [
        "A" * 35,
        EXAMPLE_CONTEXT[:10] + "N" + EXAMPLE_CONTEXT[11:],
        EXAMPLE_CONTEXT[:-1],
    ],
)
def test_rejects_invalid_or_noncanonical_contexts(context):
    with pytest.raises(ValueError, match="context_35nt|canonical NGG"):
        score_crisprscan(
            _request(guides=[{"guide_id": "zf-guide-001", "context_35nt": context}])
        )


def test_rejects_duplicate_guide_identifiers():
    guide = {"guide_id": "duplicate", "context_35nt": EXAMPLE_CONTEXT}

    with pytest.raises(ValueError, match="duplicate"):
        score_crisprscan(_request(guides=[guide, guide]))
