"""Version-aware metadata connectors for authoritative research data sources."""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any, Callable, Mapping
from urllib.request import Request, urlopen

from .species import SpeciesProfile


ENSEMBL_REST_URL = "https://rest.ensembl.org"


@dataclass(frozen=True)
class EnsemblSpeciesMetadata:
    """The subset of Ensembl metadata needed to verify a species profile."""

    name: str
    taxon_id: str
    assembly: str
    accession: str
    release: str


@dataclass(frozen=True)
class SourceCheck:
    """Result of comparing live source metadata with a local profile."""

    source: str
    matches: bool
    checked_release: str
    errors: tuple[str, ...]


class EnsemblMetadataClient:
    """Small Ensembl REST client with an injectable reader for deterministic tests."""

    def __init__(
        self,
        base_url: str = ENSEMBL_REST_URL,
        reader: Callable[[str], Mapping[str, Any]] | None = None,
    ) -> None:
        self.base_url = base_url.rstrip("/")
        self.reader = reader or _read_json

    def species_metadata(self, species_name: str) -> EnsemblSpeciesMetadata:
        payload = self.reader(f"{self.base_url}/info/species?content-type=application/json")
        species = payload.get("species")
        if not isinstance(species, list):
            raise ValueError("Ensembl response is missing the species list.")
        match = next(
            (
                item
                for item in species
                if isinstance(item, Mapping)
                and str(item.get("name", "")).casefold() == species_name.casefold()
            ),
            None,
        )
        if match is None:
            raise ValueError(f"Ensembl response does not contain {species_name!r}.")
        return EnsemblSpeciesMetadata(
            name=str(match["name"]),
            taxon_id=str(match["taxon_id"]),
            assembly=str(match["assembly"]),
            accession=str(match["accession"]),
            release=str(match["release"]),
        )


def check_ensembl_profile(
    profile: SpeciesProfile, client: EnsemblMetadataClient | None = None
) -> SourceCheck:
    """Verify that live Ensembl metadata still matches the registered profile."""
    metadata = (client or EnsemblMetadataClient()).species_metadata(profile.scientific_name)
    expected = {
        "taxon_id": profile.taxon_id,
        "assembly": profile.genome_build,
        "accession": profile.assembly_accession,
    }
    observed = {
        "taxon_id": metadata.taxon_id,
        "assembly": metadata.assembly,
        "accession": metadata.accession,
    }
    errors = tuple(
        f"{field} changed: expected {expected[field]!r}, observed {observed[field]!r}."
        for field in expected
        if expected[field] != observed[field]
    )
    return SourceCheck(
        source="Ensembl REST",
        matches=not errors,
        checked_release=metadata.release,
        errors=errors,
    )


def _read_json(url: str) -> Mapping[str, Any]:
    request = Request(url, headers={"Accept": "application/json", "User-Agent": "GeneImpact-AI/0.3"})
    with urlopen(request, timeout=20) as response:
        payload = json.load(response)
    if not isinstance(payload, Mapping):
        raise ValueError("Data source returned a non-object JSON response.")
    return payload
