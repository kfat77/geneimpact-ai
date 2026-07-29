"""Unified, integrity-checkable research dossier for animal genome-edit studies."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from hashlib import sha256
import json
import math
from pathlib import Path
import re
from typing import Any, Mapping, Sequence

from .capabilities import CapabilityStatus, capabilities_for_species
from .crispritz import import_crispritz_targets
from .crisprscan import score_crisprscan
from .housden import (
    HOUSDEN_EDIT_CLASS,
    HOUSDEN_PREDICTOR,
    HOUSDEN_SPECIES_PROFILE,
    normalize_housden,
)
from .indelphi import normalize_indelphi
from .interactions import rank_interactions
from .species import PROFILES
from .workflow import DEFAULT_MODEL_VERSION, assess_request


DOSSIER_SCHEMA_VERSION = "1.0"
MAX_TARGET_GENES = 50
MAX_DECLARED_ENDPOINTS = 20
_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
_IDENTIFIER_PATTERN = re.compile(r"^[A-Za-z0-9_.:-]{1,80}$")
_EDIT_CLASSES = {"knockout", "base_editing", "prime_editing"}
_INTENDED_CHANGES = {
    "loss_of_function",
    "gain_of_function",
    "base_substitution",
    "regulatory_change",
    "other_declared",
}
_AVAILABLE_STATUSES = {
    CapabilityStatus.AVAILABLE_DECLARED_DOMAIN,
    CapabilityStatus.AVAILABLE_REFERENCE_SEARCH,
}


@dataclass(frozen=True)
class DossierIntegrityResult:
    matches: bool
    algorithm: str
    expected_content_sha256: str
    actual_content_sha256: str
    signature_status: str
    note: str


def build_research_dossier(
    request: Mapping[str, Any],
    *,
    attachment_base_dir: Path,
    source_request_sha256: str | None = None,
    model_version: str = DEFAULT_MODEL_VERSION,
) -> dict[str, Any]:
    """Build one auditable dossier without converting predictions into safety claims."""
    if request.get("dossier_schema_version") != DOSSIER_SCHEMA_VERSION:
        raise ValueError(
            f"dossier_schema_version must be {DOSSIER_SCHEMA_VERSION!r}."
        )
    context = _mapping(request, "study_context")
    profile_key = str(context.get("species_profile", "")).strip()
    if profile_key not in PROFILES:
        raise ValueError("study_context.species_profile must be an exact registered key.")
    profile = PROFILES[profile_key]
    _validate_context(context, profile_key)
    targets = _target_genes(request.get("target_genes"))
    intended_outcomes = _text_list(
        request.get("intended_outcomes"),
        "intended_outcomes",
    )
    welfare_endpoints = _text_list(
        request.get("welfare_endpoints"),
        "welfare_endpoints",
    )
    interaction_evidence = _interaction_evidence(
        request.get("interaction_evidence", []),
        {target["gene_id"] for target in targets},
    )
    interaction_results = rank_interactions(
        {target["gene_id"]: target["evidence_signal"] for target in targets},
        {
            frozenset(item["genes"]): item["evidence_weight"]
            for item in interaction_evidence
        },
    )
    evidence = _mapping(request, "evidence")
    predictors = request.get("predictors", {})
    if not isinstance(predictors, Mapping):
        raise ValueError("predictors must be an object.")

    assessment_input = {
        "study_context": {
            "species": profile_key,
            "strain_or_breed": context["strain_or_breed"],
            "genome_build": context["genome_build"],
            "edit_class": context["edit_class"],
            "evidence_snapshot": context["evidence_snapshot_sha256"],
        },
        "evidence": evidence,
        "predictor_outputs": predictors.get("external_concern_outputs", []),
        "behive_efficiency_outputs": predictors.get(
            "behive_efficiency_outputs", []
        ),
        "behive_bystander_outputs": predictors.get(
            "behive_bystander_outputs", []
        ),
    }
    assessment = assess_request(assessment_input, model_version=model_version)
    model_predictions = list(assessment.pop("model_predictions"))
    executed_predictors = {
        item["predictor"] for item in model_predictions
    }

    if "crisprscan" in predictors:
        config = _mapping(predictors, "crisprscan")
        if config.get("guide_expression") != context["delivery_context"]:
            raise ValueError(
                "predictors.crisprscan.guide_expression must match "
                "study_context.delivery_context."
            )
        if config.get("developmental_context") != context["developmental_context"]:
            raise ValueError(
                "predictors.crisprscan.developmental_context must match "
                "study_context.developmental_context."
            )
        crisprscan_request = {
            "species_profile": profile_key,
            "genome_build": context["genome_build"],
            "assembly_accession": context["assembly_accession"],
            "reference_strain_or_isolate": profile.reference_strain,
            "nuclease": config.get("nuclease"),
            "guide_expression": config.get("guide_expression"),
            "developmental_context": config.get("developmental_context"),
            "guides": config.get("guides"),
        }
        model_predictions.append(asdict(score_crisprscan(crisprscan_request)))
        executed_predictors.add("CRISPRscan")

    if "crispritz" in predictors:
        config = _mapping(predictors, "crispritz")
        targets_path = _resolve_attachment(
            config.get("targets_file"),
            attachment_base_dir,
        )
        crispritz_metadata = {
            "species_profile": profile_key,
            "genome_build": context["genome_build"],
            "assembly_accession": context["assembly_accession"],
            "reference_strain_or_isolate": profile.reference_strain,
            "edit_class": context["edit_class"],
            "pam_definition": config.get("pam_definition"),
            "max_mismatches": config.get("max_mismatches"),
            "max_dna_bulge": config.get("max_dna_bulge"),
            "max_rna_bulge": config.get("max_rna_bulge"),
            "crispritz_commit": config.get("crispritz_commit"),
            "reference_fasta_sha256": config.get("reference_fasta_sha256"),
            "variant_aware": config.get("variant_aware"),
            "variant_snapshot_sha256": config.get("variant_snapshot_sha256"),
        }
        model_predictions.append(
            asdict(import_crispritz_targets(crispritz_metadata, targets_path))
        )
        executed_predictors.add("CRISPRitz")

    if "indelphi" in predictors:
        if profile_key != "mouse" or context["edit_class"] != "knockout":
            raise ValueError(
                "predictors.indelphi is available only for mouse knockout studies."
            )
        config = _mapping(predictors, "indelphi")
        result_files = config.get("result_files")
        if (
            not isinstance(result_files, Sequence)
            or isinstance(result_files, (str, bytes))
            or not 1 <= len(result_files) <= 100
        ):
            raise ValueError(
                "predictors.indelphi.result_files must contain 1-100 relative paths."
            )
        indelphi_predictions = []
        target_ids: set[str] = set()
        for raw_result_path in result_files:
            result_document, result_digest = _load_json_attachment(
                raw_result_path,
                attachment_base_dir,
                "predictors.indelphi.result_files",
                "inDelphi",
            )
            prediction = normalize_indelphi(
                result_document,
                source_document_sha256=result_digest,
            )
            if prediction.genome_build.casefold() != str(
                context["genome_build"]
            ).casefold():
                raise ValueError(
                    "inDelphi result genome_build must match study_context."
                )
            if prediction.assembly_accession != context["assembly_accession"]:
                raise ValueError(
                    "inDelphi result assembly_accession must match study_context."
                )
            if prediction.delivery_context != context["delivery_context"]:
                raise ValueError(
                    "inDelphi result delivery_context must match study_context."
                )
            if prediction.developmental_context != context[
                "developmental_context"
            ]:
                raise ValueError(
                    "inDelphi result developmental_context must match study_context."
                )
            if prediction.target_id in target_ids:
                raise ValueError("inDelphi result target_id values must be unique.")
            target_ids.add(prediction.target_id)
            indelphi_predictions.append(asdict(prediction))
        model_predictions.extend(indelphi_predictions)
        executed_predictors.add("inDelphi")

    if "housden" in predictors:
        if (
            profile_key != HOUSDEN_SPECIES_PROFILE
            or context["edit_class"] != HOUSDEN_EDIT_CLASS
        ):
            raise ValueError(
                "predictors.housden is available only for fruit_fly knockout studies."
            )
        config = _mapping(predictors, "housden")
        result_files = config.get("result_files")
        source_response_files = config.get("source_response_files")
        if (
            not isinstance(result_files, Sequence)
            or isinstance(result_files, (str, bytes))
            or not 1 <= len(result_files) <= 100
        ):
            raise ValueError(
                "predictors.housden.result_files must contain 1-100 relative paths."
            )
        if (
            not isinstance(source_response_files, Sequence)
            or isinstance(source_response_files, (str, bytes))
            or len(source_response_files) != len(result_files)
        ):
            raise ValueError(
                "predictors.housden.source_response_files must contain one "
                "relative XLS path per result file."
            )
        housden_predictions = []
        guide_ids: set[str] = set()
        for raw_result_path, raw_response_path in zip(
            result_files,
            source_response_files,
            strict=True,
        ):
            result_document, result_digest = _load_json_attachment(
                raw_result_path,
                attachment_base_dir,
                "predictors.housden.result_files",
                HOUSDEN_PREDICTOR,
            )
            prediction = normalize_housden(
                result_document,
                source_response=_resolve_attachment(
                    raw_response_path,
                    attachment_base_dir,
                    "predictors.housden.source_response_files",
                ).read_bytes(),
                source_document_sha256=result_digest,
            )
            if prediction.genome_build.casefold() != str(
                context["genome_build"]
            ).casefold():
                raise ValueError(
                    "Housden result genome_build must match study_context."
                )
            if prediction.assembly_accession != context["assembly_accession"]:
                raise ValueError(
                    "Housden result assembly_accession must match study_context."
                )
            if prediction.guide_expression != context["delivery_context"]:
                raise ValueError(
                    "Housden result guide_expression must match "
                    "study_context.delivery_context."
                )
            if prediction.developmental_context != context[
                "developmental_context"
            ]:
                raise ValueError(
                    "Housden result developmental_context must match study_context."
                )
            if prediction.guide_id in guide_ids:
                raise ValueError("Housden result guide_id values must be unique.")
            guide_ids.add(prediction.guide_id)
            housden_predictions.append(asdict(prediction))
        model_predictions.extend(housden_predictions)
        executed_predictors.add(HOUSDEN_PREDICTOR)

    capability_coverage = _capability_coverage(
        profile_key,
        context["edit_class"],
        executed_predictors,
    )
    available_not_run = [
        item["predictor"]
        for item in capability_coverage
        if item["execution_state"] == "available_not_run"
    ]
    review_flags = _review_flags(
        assessment,
        targets,
        interaction_evidence,
        available_not_run,
        model_predictions,
        context,
        profile.reference_strain,
    )
    canonical_request_sha256 = _canonical_sha256(request)
    if source_request_sha256 is not None and not _SHA256_PATTERN.fullmatch(
        source_request_sha256
    ):
        raise ValueError("source_request_sha256 must be a lowercase SHA-256 digest.")

    dossier = {
        "dossier_schema_version": DOSSIER_SCHEMA_VERSION,
        "dossier_id": f"geneimpact-{canonical_request_sha256[:16]}",
        "model_version": model_version,
        "input_provenance": {
            "canonical_request_sha256": canonical_request_sha256,
            "source_request_file_sha256": source_request_sha256,
            "attachment_policy": "relative_paths_confined_to_request_directory",
        },
        "study": {
            "study_id": context["study_id"],
            "species_profile": profile_key,
            "scientific_name": profile.scientific_name,
            "taxon_id": profile.taxon_id,
            "strain_or_breed": context["strain_or_breed"],
            "reference_strain_or_isolate": profile.reference_strain,
            "genome_build": context["genome_build"],
            "assembly_accession": context["assembly_accession"],
            "edit_class": context["edit_class"],
            "delivery_context": context["delivery_context"],
            "developmental_context": context["developmental_context"],
            "evidence_snapshot_sha256": context["evidence_snapshot_sha256"],
            "target_genes": targets,
            "intended_outcomes": intended_outcomes,
            "welfare_endpoints": welfare_endpoints,
        },
        "declared_evidence_inputs": {
            key: float(evidence[key])
            for key in (
                "on_target_uncertainty",
                "off_target_evidence",
                "network_impact_evidence",
                "welfare_relevance",
            )
        },
        "assessment": assessment,
        "interaction_hypotheses": {
            "total_possible_pairs": len(targets) * (len(targets) - 1) // 2,
            "evidence_supported_pairs": len(interaction_evidence),
            "declared_pair_evidence": interaction_evidence,
            "ranked_pairs": [asdict(item) for item in interaction_results],
            "interpretation": (
                "Hypothesis prioritization from declared gene signals and curated "
                "pair evidence; zero weight means missing interaction evidence."
            ),
        },
        "model_predictions": model_predictions,
        "capability_coverage": capability_coverage,
        "evidence_completeness": {
            "status": (
                "partial_predictor_coverage"
                if available_not_run
                else "all_available_relevant_predictors_included"
            ),
            "available_predictors_not_run": available_not_run,
            "prospective_empirical_validation_supplied": False,
        },
        "review_flags": review_flags,
        "report_notice": (
            "Research decision-support only. Model outputs rank hypotheses and "
            "candidate sites; they do not establish phenotype, safety, or authorization "
            "and do not replace ethics, biosafety, veterinary, or experimental review."
        ),
    }
    dossier["integrity"] = {
        "algorithm": "sha256",
        "canonicalization": "json_sorted_keys_compact_utf8_excluding_integrity",
        "content_sha256": _canonical_sha256(dossier),
        "signature_status": "unsigned",
    }
    return dossier


def verify_dossier_integrity(report: Mapping[str, Any]) -> DossierIntegrityResult:
    """Verify the dossier's self-contained content hash; this is not a signature."""
    integrity = report.get("integrity")
    if not isinstance(integrity, Mapping):
        raise ValueError("dossier integrity must be an object.")
    if integrity.get("algorithm") != "sha256":
        raise ValueError("dossier integrity algorithm must be sha256.")
    expected = integrity.get("content_sha256")
    if not isinstance(expected, str) or not _SHA256_PATTERN.fullmatch(expected):
        raise ValueError("dossier integrity content_sha256 is invalid.")
    unsigned_content = dict(report)
    unsigned_content.pop("integrity", None)
    actual = _canonical_sha256(unsigned_content)
    return DossierIntegrityResult(
        matches=actual == expected,
        algorithm="sha256",
        expected_content_sha256=expected,
        actual_content_sha256=actual,
        signature_status=str(integrity.get("signature_status", "unknown")),
        note=(
            "The content hash detects changes but is not an identity or authenticity signature."
        ),
    )


def _validate_context(context: Mapping[str, Any], profile_key: str) -> None:
    required = (
        "study_id",
        "species_profile",
        "strain_or_breed",
        "genome_build",
        "assembly_accession",
        "edit_class",
        "delivery_context",
        "developmental_context",
        "evidence_snapshot_sha256",
    )
    missing = [key for key in required if key not in context]
    if missing:
        raise ValueError(
            "study_context is missing required fields: " + ", ".join(missing)
        )
    if not _IDENTIFIER_PATTERN.fullmatch(str(context["study_id"])):
        raise ValueError("study_context.study_id must be a safe 1-80 character identifier.")
    profile = PROFILES[profile_key]
    accepted_builds = {
        value.casefold()
        for value in (profile.genome_build, *profile.accepted_build_names)
    }
    if str(context["genome_build"]).casefold() not in accepted_builds:
        raise ValueError(
            f"study_context.genome_build must match registered {profile.genome_build!r}."
        )
    if context["assembly_accession"] != profile.assembly_accession:
        raise ValueError(
            "study_context.assembly_accession must match registered "
            f"{profile.assembly_accession!r}."
        )
    if context["edit_class"] not in _EDIT_CLASSES:
        raise ValueError(
            "study_context.edit_class must be knockout, base_editing, or prime_editing."
        )
    for field in ("strain_or_breed", "delivery_context", "developmental_context"):
        value = context[field]
        if not isinstance(value, str) or not value.strip() or len(value) > 160:
            raise ValueError(f"study_context.{field} must be a non-empty short string.")
    if not _SHA256_PATTERN.fullmatch(str(context["evidence_snapshot_sha256"])):
        raise ValueError(
            "study_context.evidence_snapshot_sha256 must be a lowercase SHA-256 digest."
        )


def _target_genes(raw: Any) -> list[dict[str, Any]]:
    if (
        not isinstance(raw, Sequence)
        or isinstance(raw, (str, bytes))
        or not 1 <= len(raw) <= MAX_TARGET_GENES
    ):
        raise ValueError(
            f"target_genes must contain between 1 and {MAX_TARGET_GENES} entries."
        )
    targets: list[dict[str, Any]] = []
    seen: set[str] = set()
    for index, item in enumerate(raw, start=1):
        if not isinstance(item, Mapping):
            raise ValueError(f"target gene {index} must be an object.")
        gene_id = str(item.get("gene_id", "")).strip()
        gene_symbol = str(item.get("gene_symbol", "")).strip()
        intended_change = item.get("intended_change")
        reference = str(item.get("evidence_reference", "")).strip()
        if not _IDENTIFIER_PATTERN.fullmatch(gene_id) or gene_id in seen:
            raise ValueError(f"target gene {index} has an invalid or duplicate gene_id.")
        if not _IDENTIFIER_PATTERN.fullmatch(gene_symbol):
            raise ValueError(f"target gene {index} gene_symbol is invalid.")
        if intended_change not in _INTENDED_CHANGES:
            raise ValueError(f"target gene {index} intended_change is unsupported.")
        if not reference or len(reference) > 512:
            raise ValueError(f"target gene {index} evidence_reference is required.")
        signal = _bounded(item.get("evidence_signal"), f"target gene {index} evidence_signal")
        seen.add(gene_id)
        targets.append(
            {
                "gene_id": gene_id,
                "gene_symbol": gene_symbol,
                "intended_change": intended_change,
                "evidence_signal": signal,
                "evidence_reference": reference,
            }
        )
    return targets


def _interaction_evidence(
    raw: Any,
    target_ids: set[str],
) -> list[dict[str, Any]]:
    if not isinstance(raw, list):
        raise ValueError("interaction_evidence must be a list.")
    maximum = len(target_ids) * (len(target_ids) - 1) // 2
    if len(raw) > maximum:
        raise ValueError("interaction_evidence contains more rows than possible pairs.")
    results: list[dict[str, Any]] = []
    seen: set[frozenset[str]] = set()
    for index, item in enumerate(raw, start=1):
        if not isinstance(item, Mapping):
            raise ValueError(f"interaction evidence {index} must be an object.")
        genes = item.get("genes")
        if (
            not isinstance(genes, list)
            or len(genes) != 2
            or genes[0] == genes[1]
            or any(gene not in target_ids for gene in genes)
        ):
            raise ValueError(
                f"interaction evidence {index} must reference two distinct target genes."
            )
        pair = frozenset(genes)
        if pair in seen:
            raise ValueError(f"interaction evidence {index} duplicates a gene pair.")
        reference = str(item.get("evidence_reference", "")).strip()
        if not reference or len(reference) > 512:
            raise ValueError(
                f"interaction evidence {index} evidence_reference is required."
            )
        seen.add(pair)
        results.append(
            {
                "genes": tuple(sorted(genes)),
                "evidence_weight": _bounded(
                    item.get("evidence_weight"),
                    f"interaction evidence {index} evidence_weight",
                ),
                "evidence_reference": reference,
            }
        )
    return results


def _text_list(raw: Any, label: str) -> list[str]:
    if (
        not isinstance(raw, Sequence)
        or isinstance(raw, (str, bytes))
        or not 1 <= len(raw) <= MAX_DECLARED_ENDPOINTS
    ):
        raise ValueError(
            f"{label} must contain between 1 and {MAX_DECLARED_ENDPOINTS} entries."
        )
    values = []
    for index, item in enumerate(raw, start=1):
        if not isinstance(item, str) or not item.strip() or len(item) > 300:
            raise ValueError(f"{label} item {index} must be a non-empty short string.")
        values.append(item.strip())
    if len(set(values)) != len(values):
        raise ValueError(f"{label} must not contain duplicates.")
    return values


def _capability_coverage(
    profile_key: str,
    edit_class: str,
    executed_predictors: set[str],
) -> list[dict[str, Any]]:
    coverage = []
    for capability in capabilities_for_species(profile_key):
        relevant = edit_class in capability.edit_classes
        if capability.predictor in executed_predictors:
            execution_state = "included"
        elif not relevant:
            execution_state = "irrelevant_edit_class"
        elif capability.status in _AVAILABLE_STATUSES:
            execution_state = "available_not_run"
        elif capability.status is CapabilityStatus.OUT_OF_DOMAIN_ONLY:
            execution_state = "out_of_domain"
        else:
            execution_state = "not_integrated"
        row = asdict(capability)
        row["status"] = capability.status.value
        row["relevant_to_edit_class"] = relevant
        row["execution_state"] = execution_state
        coverage.append(row)
    return coverage


def _review_flags(
    assessment: Mapping[str, Any],
    targets: Sequence[Mapping[str, Any]],
    interaction_evidence: Sequence[Mapping[str, Any]],
    available_not_run: Sequence[str],
    model_predictions: Sequence[Mapping[str, Any]],
    context: Mapping[str, Any],
    reference_strain: str,
) -> list[str]:
    flags = list(assessment["species_validation"]["warnings"])
    if assessment["assessment"]["tier"] != "standard_review":
        flags.append(
            f"Evidence triage requires {assessment['assessment']['tier']}."
        )
    possible_pairs = len(targets) * (len(targets) - 1) // 2
    if possible_pairs > len(interaction_evidence):
        flags.append(
            f"Curated interaction evidence covers {len(interaction_evidence)} of "
            f"{possible_pairs} target-gene pairs."
        )
    if available_not_run:
        flags.append(
            "Available relevant predictors not run: " + ", ".join(available_not_run)
        )
    if (
        context["strain_or_breed"].casefold() != reference_strain.casefold()
        and any(
            item.get("predictor") == "CRISPRitz" and not item.get("variant_aware")
            for item in model_predictions
        )
    ):
        flags.append(
            "Reference-only CRISPRitz search does not include declared study-strain variants."
        )
    flags.append(
        "Target gene identifiers are researcher-declared and have not yet been resolved against a versioned species annotation."
    )
    flags.append(
        "No prospective empirical validation result is attached to this dossier."
    )
    return flags


def _resolve_attachment(
    raw_path: Any,
    base_dir: Path,
    label: str = "predictors.crispritz.targets_file",
) -> Path:
    if not isinstance(raw_path, str) or not raw_path.strip():
        raise ValueError(f"{label} requires a non-empty relative path.")
    relative = Path(raw_path)
    if relative.is_absolute():
        raise ValueError("attachment paths must be relative to the request directory.")
    resolved_base = base_dir.resolve()
    resolved = (resolved_base / relative).resolve()
    if not resolved.is_relative_to(resolved_base):
        raise ValueError("attachment path escapes the request directory.")
    if not resolved.is_file():
        raise ValueError(f"attachment file does not exist: {raw_path}")
    return resolved


def _load_json_attachment(
    raw_path: Any,
    base_dir: Path,
    label: str,
    predictor: str,
) -> tuple[Mapping[str, Any], str]:
    path = _resolve_attachment(raw_path, base_dir, label)
    content = path.read_bytes()
    try:
        document = json.loads(content.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(
            f"{predictor} result is not valid UTF-8 JSON: {raw_path}"
        ) from error
    if not isinstance(document, Mapping):
        raise ValueError(f"each {predictor} result must be a JSON object.")
    return document, sha256(content).hexdigest()


def _mapping(parent: Mapping[str, Any], key: str) -> Mapping[str, Any]:
    value = parent.get(key)
    if not isinstance(value, Mapping):
        raise ValueError(f"{key} must be an object.")
    return value


def _bounded(value: Any, label: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{label} must be numeric.")
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be numeric.") from error
    if not math.isfinite(number) or not 0 <= number <= 1:
        raise ValueError(f"{label} must be between 0 and 1.")
    return number


def _canonical_sha256(value: Mapping[str, Any]) -> str:
    try:
        canonical = json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
        )
    except (TypeError, ValueError) as error:
        raise ValueError("dossier content must be JSON serializable.") from error
    return sha256(canonical.encode("utf-8")).hexdigest()
