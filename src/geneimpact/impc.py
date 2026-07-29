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
    significant: bool
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
    pages: int
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
        return self.gene_phenotypes(marker_symbol, significant=True, rows=rows)

    def gene_phenotypes(
        self,
        marker_symbol: str,
        *,
        significant: bool | None = None,
        rows: int = 1000,
        max_documents: int = 5000,
    ) -> ImpcGeneEvidence:
        """Fetch bounded statistical results, optionally filtered by significance."""
        marker_symbol = marker_symbol.strip()
        if not marker_symbol:
            raise ValueError("marker_symbol is required.")
        if not 1 <= rows <= 1000:
            raise ValueError("rows must be between 1 and 1000.")
        query = f"marker_symbol:{marker_symbol}"
        if significant is not None:
            query += f" AND significant:{str(significant).lower()}"
        if max_documents < rows:
            raise ValueError("max_documents must be at least as large as rows.")
        documents: list[Any] = []
        start = 0
        pages = 0
        source_url = ""
        num_found = 0
        while True:
            parameters = urlencode(
                {
                    "q": query,
                    "fl": ",".join(STATISTICAL_FIELDS),
                    "rows": rows,
                    "start": start,
                    "wt": "json",
                }
            )
            page_url = f"{self.base_url}/statistical-result/select?{parameters}"
            source_url = source_url or page_url
            payload = self.reader(page_url)
            response = payload.get("response")
            if not isinstance(response, Mapping) or not isinstance(response.get("docs"), list):
                raise ValueError("IMPC response is missing response.docs.")
            num_found = int(response.get("numFound", len(response["docs"])))
            if num_found > max_documents:
                raise ValueError(
                    f"IMPC returned {num_found} results for {marker_symbol!r}, exceeding "
                    f"the run limit {max_documents}; use a narrower query."
                )
            page = response["docs"]
            documents.extend(page)
            pages += 1
            if len(documents) >= num_found:
                break
            if not page:
                raise ValueError("IMPC pagination ended before all reported results were returned.")
            start = len(documents)

        if len(documents) != num_found:
            raise ValueError(
                f"IMPC returned {len(documents)} documents but reported {num_found}."
            )
        results = tuple(_normalize(document) for document in documents)
        return ImpcGeneEvidence(
            source_url=source_url,
            retrieved_at=datetime.now(timezone.utc).isoformat(),
            marker_symbol=marker_symbol,
            num_found=num_found,
            pages=pages,
            results=results,
        )


def _normalize(document: Any) -> ImpcGenePhenotype:
    if not isinstance(document, Mapping):
        raise ValueError("IMPC result document must be an object.")
    return ImpcGenePhenotype(
        marker_symbol=str(document["marker_symbol"]),
        significant=bool(document.get("significant", False)),
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
