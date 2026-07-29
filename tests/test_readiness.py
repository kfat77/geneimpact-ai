import json
import sys

from geneimpact.cli import main
from geneimpact.readiness import (
    EvidenceUseStatus,
    readiness_for_species,
)


def test_hazard_evidence_cannot_be_promoted_to_predictive_adapter():
    for species_profile in ("rhesus_macaque", "cynomolgus_macaque"):
        readiness = readiness_for_species(species_profile)

        assert readiness.predictive_adapter_available is False
        assert readiness.evidence_records
        assert all(
            record.status is EvidenceUseStatus.HAZARD_EVIDENCE_ONLY
            for record in readiness.evidence_records
        )
        assert not any(
            record.eligible_for_predictive_capability
            for record in readiness.evidence_records
        )


def test_fruit_fly_readiness_is_limited_to_the_declared_cell_domain():
    readiness = readiness_for_species("fruit_fly")
    housden = next(
        record
        for record in readiness.evidence_records
        if record.predictor_or_method == "Housden"
    )

    assert readiness.predictive_adapter_available is True
    assert housden.status is EvidenceUseStatus.USABLE_ADAPTER
    assert housden.eligible_for_predictive_capability is True
    assert housden.biological_domain == "Drosophila S2R+ cell culture"
    assert "in vivo" in housden.limitations


def test_rat_readiness_reports_the_public_label_gap():
    readiness = readiness_for_species("rat")

    assert readiness.predictive_adapter_available is False
    assert readiness.evidence_records[0].status is (
        EvidenceUseStatus.INSUFFICIENT_PUBLIC_DATA
    )
    assert readiness.evidence_records[0].labels_public is False


def test_readiness_cli_emits_machine_readable_species_report(
    monkeypatch,
    capsys,
):
    monkeypatch.setattr(
        sys,
        "argv",
        ["geneimpact", "readiness", "--species", "cynomolgus_macaque"],
    )

    main()

    report = json.loads(capsys.readouterr().out)
    assert report["species_profile"] == "cynomolgus_macaque"
    assert report["predictive_adapter_available"] is False
    assert report["evidence_records"][0]["status"] == "hazard_evidence_only"
