"""Bounded audit import for externally executed CRISPRitz target searches."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
import re
from typing import Any, Mapping

from .species import PROFILES


CRISPRITZ_VERSION = "v2.7.0"
CRISPRITZ_COMMIT = "24b893ecb0c2354d5c76697e116d2febe1ee6265"
CRISPRITZ_REFERENCE = (
    "https://github.com/pinellolab/CRISPRitz/"
    f"tree/{CRISPRITZ_COMMIT}"
)
MAX_REPORTED_HITS = 100
MAX_GUIDES = 10_000
_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
_SEQUENCE_PATTERN = re.compile(r"^[ACGTNRYSWKMBDHV-]+$")


@dataclass(frozen=True)
class CrispritzHit:
    guide_sha256: str
    target_sha256: str
    chromosome: str
    position: int
    cluster_position: int
    direction: str
    mismatches: int
    bulge_type: str
    bulge_size: int
    total_differences: int


@dataclass(frozen=True)
class CrispritzDifferenceCount:
    total_differences: int
    hit_count: int


@dataclass(frozen=True)
class CrispritzAuditReport:
    predictor: str
    predictor_version: str
    predictor_commit: str
    task: str
    species_profile: str
    genome_build: str
    assembly_accession: str
    reference_strain_or_isolate: str
    edit_class: str
    pam_definition: str
    max_mismatches: int
    max_dna_bulge: int
    max_rna_bulge: int
    variant_aware: bool
    reference_fasta_sha256: str
    variant_snapshot_sha256: str | None
    targets_file_sha256: str
    observed_guide_count: int
    chromosome_count: int
    candidate_site_count: int
    exact_sequence_match_count: int
    counts_by_total_difference: tuple[CrispritzDifferenceCount, ...]
    reported_hit_count: int
    top_candidate_hits: tuple[CrispritzHit, ...]
    evidence_reference: str
    execution_status: str
    license_notice: str
    warnings: tuple[str, ...]


def import_crispritz_targets(
    metadata: Mapping[str, Any], targets_path: Path
) -> CrispritzAuditReport:
    """Validate one CRISPRitz targets file and return a bounded audit record."""
    declaration = _validate_metadata(metadata)
    targets_digest = _sha256_file(targets_path)
    counts: dict[int, int] = {}
    guide_hashes: set[str] = set()
    chromosomes: set[str] = set()
    top_hits: list[CrispritzHit] = []
    total_hits = 0
    exact_matches = 0

    with targets_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        _validate_header(reader.fieldnames)
        for row_number, row in enumerate(reader, start=2):
            hit = _parse_hit(row, row_number, declaration)
            total_hits += 1
            exact_matches += int(hit.total_differences == 0)
            counts[hit.total_differences] = counts.get(hit.total_differences, 0) + 1
            guide_hashes.add(hit.guide_sha256)
            if len(guide_hashes) > MAX_GUIDES:
                raise ValueError(f"targets file may contain at most {MAX_GUIDES} unique guides.")
            chromosomes.add(hit.chromosome)
            top_hits.append(hit)
            if len(top_hits) > MAX_REPORTED_HITS * 2:
                top_hits.sort(key=_hit_sort_key)
                del top_hits[MAX_REPORTED_HITS:]

    if total_hits == 0:
        raise ValueError("CRISPRitz targets file contains no result rows.")
    top_hits.sort(key=_hit_sort_key)
    del top_hits[MAX_REPORTED_HITS:]

    return CrispritzAuditReport(
        predictor="CRISPRitz",
        predictor_version=CRISPRITZ_VERSION,
        predictor_commit=CRISPRITZ_COMMIT,
        task="reference_genome_off_target_enumeration",
        species_profile=declaration["species_profile"],
        genome_build=declaration["genome_build"],
        assembly_accession=declaration["assembly_accession"],
        reference_strain_or_isolate=declaration["reference_strain_or_isolate"],
        edit_class=declaration["edit_class"],
        pam_definition=declaration["pam_definition"],
        max_mismatches=declaration["max_mismatches"],
        max_dna_bulge=declaration["max_dna_bulge"],
        max_rna_bulge=declaration["max_rna_bulge"],
        variant_aware=declaration["variant_aware"],
        reference_fasta_sha256=declaration["reference_fasta_sha256"],
        variant_snapshot_sha256=declaration["variant_snapshot_sha256"],
        targets_file_sha256=targets_digest,
        observed_guide_count=len(guide_hashes),
        chromosome_count=len(chromosomes),
        candidate_site_count=total_hits,
        exact_sequence_match_count=exact_matches,
        counts_by_total_difference=tuple(
            CrispritzDifferenceCount(difference, counts[difference])
            for difference in sorted(counts)
        ),
        reported_hit_count=len(top_hits),
        top_candidate_hits=tuple(top_hits),
        evidence_reference=CRISPRITZ_REFERENCE,
        execution_status="externally_executed_output_validated_not_recomputed",
        license_notice=(
            "CRISPRitz declares AGPL availability for academic research and a "
            "separate commercial license requirement; GeneImpact AI does not redistribute it."
        ),
        warnings=(
            "Enumerated sites are sequence-similarity candidates, not measured cleavage events or calibrated probabilities.",
            "An exact sequence match is not assumed to be the intended locus; intended and unintended loci require explicit coordinates.",
            "A reference assembly does not represent every animal, colony, stock, or population variant.",
            "Scores and hit counts do not establish edit safety and require empirical off-target validation.",
        ),
    )


def _validate_metadata(metadata: Mapping[str, Any]) -> dict[str, Any]:
    required = (
        "species_profile",
        "genome_build",
        "assembly_accession",
        "reference_strain_or_isolate",
        "edit_class",
        "pam_definition",
        "max_mismatches",
        "max_dna_bulge",
        "max_rna_bulge",
        "crispritz_commit",
        "reference_fasta_sha256",
        "variant_aware",
    )
    missing = [key for key in required if key not in metadata]
    if missing:
        raise ValueError(f"CRISPRitz metadata is missing fields: {', '.join(missing)}")
    profile_key = str(metadata["species_profile"])
    if profile_key not in PROFILES:
        raise ValueError(f"unknown species_profile {profile_key!r}.")
    profile = PROFILES[profile_key]
    accepted_builds = {
        value.casefold()
        for value in (profile.genome_build, *profile.accepted_build_names)
    }
    genome_build = str(metadata["genome_build"])
    if genome_build.casefold() not in accepted_builds:
        raise ValueError(
            f"genome_build must match registered {profile.genome_build!r}."
        )
    if metadata["assembly_accession"] != profile.assembly_accession:
        raise ValueError(
            f"assembly_accession must match registered {profile.assembly_accession!r}."
        )
    reference_context = str(metadata["reference_strain_or_isolate"]).strip()
    if reference_context.casefold() != profile.reference_strain.casefold():
        raise ValueError(
            "reference_strain_or_isolate must match registered "
            f"{profile.reference_strain!r}; study animals belong in a separate study context."
        )
    if metadata["crispritz_commit"] != CRISPRITZ_COMMIT:
        raise ValueError(
            f"adapter is verified only for CRISPRitz commit {CRISPRITZ_COMMIT}."
        )
    edit_class = str(metadata["edit_class"])
    if edit_class not in {"knockout", "base_editing", "prime_editing"}:
        raise ValueError("edit_class must be knockout, base_editing, or prime_editing.")
    pam_definition = str(metadata["pam_definition"]).strip().upper()
    if not pam_definition or " " not in pam_definition:
        raise ValueError("pam_definition must retain the CRISPRitz sequence and PAM-length format.")
    thresholds = {
        key: _bounded_integer(metadata[key], key, maximum=10)
        for key in ("max_mismatches", "max_dna_bulge", "max_rna_bulge")
    }
    reference_digest = str(metadata["reference_fasta_sha256"])
    if not _SHA256_PATTERN.fullmatch(reference_digest):
        raise ValueError("reference_fasta_sha256 must be a lowercase SHA-256 digest.")
    if metadata["variant_aware"] not in {True, False}:
        raise ValueError("variant_aware must be a boolean.")
    variant_digest = metadata.get("variant_snapshot_sha256")
    if metadata["variant_aware"]:
        if not isinstance(variant_digest, str) or not _SHA256_PATTERN.fullmatch(variant_digest):
            raise ValueError(
                "variant_snapshot_sha256 is required for a variant-aware search."
            )
    elif variant_digest is not None:
        raise ValueError("variant_snapshot_sha256 must be null when variant_aware is false.")
    return {
        "species_profile": profile_key,
        "genome_build": genome_build,
        "assembly_accession": profile.assembly_accession,
        "reference_strain_or_isolate": profile.reference_strain,
        "edit_class": edit_class,
        "pam_definition": pam_definition,
        "variant_aware": metadata["variant_aware"],
        "reference_fasta_sha256": reference_digest,
        "variant_snapshot_sha256": variant_digest,
        **thresholds,
    }


def _validate_header(fieldnames: list[str] | None) -> None:
    if fieldnames is None:
        raise ValueError("CRISPRitz targets file is missing its header.")
    normalized = {_normalize_header(name) for name in fieldnames}
    required = {
        "bulge_type",
        "crrna",
        "dna",
        "chromosome",
        "position",
        "cluster_position",
        "direction",
        "mismatches",
        "bulge_size",
        "total",
    }
    missing = sorted(required - normalized)
    if missing:
        raise ValueError(f"CRISPRitz targets header is missing: {', '.join(missing)}")


def _parse_hit(
    row: Mapping[str, str], row_number: int, declaration: Mapping[str, Any]
) -> CrispritzHit:
    if None in row or any(value is None for value in row.values()):
        raise ValueError(f"row {row_number} does not match the CRISPRitz targets header.")
    values = {_normalize_header(key): value for key, value in row.items()}
    guide = values["crrna"].strip().upper()
    target = values["dna"].strip().upper()
    if not guide or not _SEQUENCE_PATTERN.fullmatch(guide):
        raise ValueError(f"row {row_number} has an invalid crRNA sequence.")
    if not target or not _SEQUENCE_PATTERN.fullmatch(target):
        raise ValueError(f"row {row_number} has an invalid DNA sequence.")
    mismatch = _bounded_integer(
        values["mismatches"], f"row {row_number} mismatches", declaration["max_mismatches"]
    )
    bulge_size = _bounded_integer(values["bulge_size"], f"row {row_number} bulge_size", 10)
    bulge_type = values["bulge_type"].strip().upper()
    if bulge_type not in {"X", "DNA", "RNA", "RNA,DNA"}:
        raise ValueError(f"row {row_number} has an unsupported bulge type.")
    if bulge_type == "X" and bulge_size != 0:
        raise ValueError(f"row {row_number} has a bulge size for type X.")
    if bulge_type == "DNA" and bulge_size > declaration["max_dna_bulge"]:
        raise ValueError(f"row {row_number} exceeds max_dna_bulge.")
    if bulge_type == "RNA" and bulge_size > declaration["max_rna_bulge"]:
        raise ValueError(f"row {row_number} exceeds max_rna_bulge.")
    if bulge_type == "RNA,DNA" and bulge_size > (
        declaration["max_rna_bulge"] + declaration["max_dna_bulge"]
    ):
        raise ValueError(f"row {row_number} exceeds combined bulge thresholds.")
    total = _bounded_integer(
        values["total"],
        f"row {row_number} total",
        declaration["max_mismatches"]
        + declaration["max_dna_bulge"]
        + declaration["max_rna_bulge"],
    )
    if total != mismatch + bulge_size:
        raise ValueError(f"row {row_number} total must equal mismatches plus bulge_size.")
    direction = values["direction"].strip()
    if direction not in {"+", "-"}:
        raise ValueError(f"row {row_number} direction must be + or -.")
    chromosome = values["chromosome"].strip()
    if not chromosome:
        raise ValueError(f"row {row_number} chromosome is required.")
    return CrispritzHit(
        guide_sha256=sha256(guide.encode("ascii")).hexdigest(),
        target_sha256=sha256(target.encode("ascii")).hexdigest(),
        chromosome=chromosome,
        position=_bounded_integer(values["position"], f"row {row_number} position", 10**12),
        cluster_position=_bounded_integer(
            values["cluster_position"], f"row {row_number} cluster_position", 10**12
        ),
        direction=direction,
        mismatches=mismatch,
        bulge_type=bulge_type,
        bulge_size=bulge_size,
        total_differences=total,
    )


def _normalize_header(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.lstrip("#").strip().casefold()).strip("_")


def _bounded_integer(value: Any, label: str, maximum: int) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{label} must be an integer.")
    try:
        number = int(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be an integer.") from error
    if str(value).strip() != str(number) or not 0 <= number <= maximum:
        raise ValueError(f"{label} must be between 0 and {maximum}.")
    return number


def _hit_sort_key(hit: CrispritzHit) -> tuple[int, int, int, str, int]:
    return (
        hit.total_differences,
        hit.mismatches,
        hit.bulge_size,
        hit.chromosome,
        hit.position,
    )


def _sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
