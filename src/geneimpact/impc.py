"""Read-only access to significant IMPC gene-phenotype results."""

from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any, Callable, Mapping
from urllib.parse import urlencode
from urllib.request import Request, urlopen


IMPC_SOLR_URL = "https://www.ebi.ac.uk/mi/impc/solr"
STATISTICAL_FIELDS = (
    "marker_symbol",
    "mp_term_id",
    "mp_term_name",
    "top_level_mp_term_name",
    "effect_size",
    "p_value",
    "significant",
    "procedure_name",
    "parameter_name",
    "sex",
    "zygosity",
)


@dataclass(frozen=True)
class ImpcGenePhenotype:
    """A significant IMPC statistical result with experimental context."""

    marker_symbol: str
    mp_term_id: str | None
    mp_term_name: str | None
    top_level_mp_terms: tuple[str, ...]
    effect_size: float | None
    p_value: float | None
    procedure_name: str | None
    parameter_name: str | None
    sex: str | None
    zygosity: str | None


@dataclass(frozen=True)
class ImpcGeneEvidence:
    """Auditable result of one bounded IMPC query."""

    source_url: str
    retrieved_at: str
    marker_symbol: str
    num_found: int
    results: tuple[ImpcGenePhenotype, ...]


class ImpcClient:
    """Bounded client for the IMPC statistical-result Solr core."""

    def __init__(
        self,
        base_url: str = IMPC_SOLR_URL,
        reader: Callable[[str], Mapping[str, Any]] | None = None,
    ) -> None:
        self.base_url = base_url.rstrip("/")
        self.reader = reader or _read_json

    def significant_gene_phenotypes(
        self, marker_symbol: str, *, rows: int = 1000
    ) -> ImpcGeneEvidence:
        marker_symbol = marker_symbol.strip()
        if not marker_symbol:
            raise ValueError("marker_symbol is required.")
        if not 1 <= rows <= 1000:
            raise ValueError("rows must be between 1 and 1000.")
        query = f"marker_symbol:{marker_symbol} AND significant:true"
        parameters = urlencode(
            {
                "q": query,
                "fl": ",".join(STATISTICAL_FIELDS),
                "rows": rows,
                "wt": "json",
            }
        )
        source_url = f"{self.base_url}/statistical-result/select?{parameters}"
        payload = self.reader(source_url)
        response = payload.get("response")
        if not isinstance(response, Mapping) or not isinstance(response.get("docs"), list):
            raise ValueError("IMPC response is missing response.docs.")
        results = tuple(_normalize(document) for document in response["docs"])
        return ImpcGeneEvidence(
            source_url=source_url,
            retrieved_at=datetime.now(timezone.utc).isoformat(),
            marker_symbol=marker_symbol,
            num_found=int(response.get("numFound", len(results))),
            results=results,
        )


def _normalize(document: Any) -> ImpcGenePhenotype:
    if not isinstance(document, Mapping):
        raise ValueError("IMPC result document must be an object.")
    return ImpcGenePhenotype(
        marker_symbol=str(document["marker_symbol"]),
        mp_term_id=_optional_string(document.get("mp_term_id")),
        mp_term_name=_optional_string(document.get("mp_term_name")),
        top_level_mp_terms=tuple(
            str(value) for value in document.get("top_level_mp_term_name", [])
        ),
        effect_size=_optional_float(document.get("effect_size")),
        p_value=_optional_float(document.get("p_value")),
        procedure_name=_optional_string(document.get("procedure_name")),
        parameter_name=_optional_string(document.get("parameter_name")),
        sex=_optional_string(document.get("sex")),
        zygosity=_optional_string(document.get("zygosity")),
    )


def _optional_string(value: Any) -> str | None:
    return None if value is None else str(value)


def _optional_float(value: Any) -> float | None:
    return None if value is None else float(value)


def _read_json(url: str) -> Mapping[str, Any]:
    request = Request(url, headers={"Accept": "application/json", "User-Agent": "GeneImpact-AI/0.3"})
    with urlopen(request, timeout=30) as response:
        payload = json.load(response)
    if not isinstance(payload, Mapping):
        raise ValueError("IMPC returned a non-object JSON response.")
    return payload
