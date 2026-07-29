import json
import sys

from geneimpact.cli import main
from geneimpact.readiness import (
    EvidenceUseStatus,
    readiness_for_species,
)


def test_hazard_evidence_cannot_be_promoted_to_predictive_adapter():
    readiness = readiness_for_species("rhesus_macaque")

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


def test_cynomolgus_base_editing_is_bounded_benchmark_not_predictor():
    readiness = readiness_for_species("cynomolgus_macaque")
    benchmark = next(
        record
        for record in readiness.evidence_records
        if record.task == "base_editing_embryo_transfer_validation"
    )

    assert readiness.predictive_adapter_available is False
    assert benchmark.status is EvidenceUseStatus.USABLE_BOUNDED_BENCHMARK
    assert benchmark.eligible_for_predictive_capability is False
    assert benchmark.labels_public is True
    assert benchmark.target_count == 11
    assert benchmark.sample_count == 273
    assert "not independent" in benchmark.limitations


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

    cas12a = next(
        record
        for record in readiness.evidence_records
        if record.task == "in_vivo_cas12a_array_loh_evidence"
    )
    assert cas12a.status is EvidenceUseStatus.USABLE_BOUNDED_BENCHMARK
    assert cas12a.eligible_for_predictive_capability is False
    assert cas12a.labels_public is True
    assert cas12a.target_count == 845
    assert "array" in cas12a.limitations


def test_rat_readiness_reports_transfer_evidence_without_predictor_promotion():
    readiness = readiness_for_species("rat")

    assert readiness.predictive_adapter_available is False
    assert readiness.evidence_records[0].status is (
        EvidenceUseStatus.TRANSFER_EVIDENCE_ONLY
    )
    assert readiness.evidence_records[0].labels_public is True
    assert readiness.evidence_records[0].target_count == 14
    assert readiness.evidence_records[0].sample_count == 186


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
