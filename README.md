<div align="center">

# 🧬 GeneImpact AI

[![Tests](https://github.com/kfat77/geneimpact-ai/actions/workflows/test.yml/badge.svg)](https://github.com/kfat77/geneimpact-ai/actions/workflows/test.yml)
[![License: MIT](https://img.shields.io/github/license/kfat77/geneimpact-ai?color=10b981)](LICENSE)
[![Version](https://img.shields.io/badge/version-0.15.0-2563eb)](https://github.com/kfat77/geneimpact-ai)
[![Python](https://img.shields.io/badge/python-3.11%2B-0ea5e9?logo=python&logoColor=white)](pyproject.toml)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-10b981.svg)](https://github.com/kfat77/geneimpact-ai/pulls)

**Evidence-aware genome-editing analysis that turns declared study context, predictor outputs, and biological evidence into reproducible, uncertainty-bounded audit reports.**

</div>

## 🌟 Live Interactive Demo

> ### Explore GeneImpact AI before installing anything
>
> The light-theme, bilingual demo runs a synthetic analysis entirely in your browser, with real-time biostatistical feedback, species validation, multi-omics context checks, and an audit-style result view.

[![GeneImpact AI interactive demo preview](docs/assets/demo-preview.gif)](https://kfat77.github.io/geneimpact-ai/)

🚀 **[Try the Live Demo Here](https://kfat77.github.io/geneimpact-ai/)**

> [!IMPORTANT]
> The Web Demo uses synthetic data and illustrative outputs. It is a product preview—not experimental evidence, an edit-design service, or a safety determination.

## Key Features

- 🧬 **Evidence-aware multi-omics integration** — Combine versioned reference assemblies, functional annotations, phenotype evidence, target-gene context, and independently produced predictor outputs without erasing source provenance.
- 📊 **Real-time statistical feedback** — Surface bounded evidence scores, applicability status, calibration metrics, uncertainty, and animal-welfare relevance instead of presenting an opaque pass/fail result.
- 🔌 **Developer-friendly API and CLI** — Use typed Python interfaces or composable command-line workflows to generate JSON reports, verify dossier integrity, inspect species readiness, and integrate external predictors.
- 🧪 **Strict biological-domain gating** — Bind every result to a species, strain or isolate, genome build, edit class, delivery context, evidence snapshot, and model version.
- 🔍 **Audit-ready provenance** — Preserve request hashes, model and source versions, reference accessions, evidence references, limitations, and machine-readable report notices.
- 🐁 **Multi-species research contexts** — Validate registered contexts for mouse, rat, zebrafish, fruit fly, rhesus macaque, and cynomolgus macaque while reporting predictor maturity separately for each domain.

## Quick Start

GeneImpact AI requires **Python 3.11 or later**.

```bash
git clone https://github.com/kfat77/geneimpact-ai.git
cd geneimpact-ai

python -m venv .venv
source .venv/bin/activate

python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
pytest
```

Run the included synthetic assessment and write a machine-readable report:

```bash
geneimpact assess examples/assessment-request.json \
  --output assessment-report.json
```

Build and verify the preferred unified research dossier:

```bash
geneimpact dossier examples/dossier-zebrafish-request.json \
  --output research-dossier.json

geneimpact verify-dossier research-dossier.json
```

Preview the bilingual Web Demo locally:

```bash
python -m http.server 8000 --directory demo
```

Then open `http://localhost:8000`.

## Code Example / Usage

The example below evaluates a **synthetic 35-nt SpCas9 context** with the version-locked CRISPRscan implementation and returns a provenance-bearing report. This predictor is deliberately limited to its declared zebrafish embryo domain.

```python
import json
from dataclasses import asdict

from geneimpact.crisprscan import score_crisprscan

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
            "context_35nt": "ACCTGGATCGATGCTGATGCTAGATAAGGTTGAGC",
        }
    ],
}

report = score_crisprscan(request)
print(json.dumps(asdict(report), indent=2))
```

Equivalent CLI workflow:

```bash
geneimpact score-crisprscan \
  --input examples/crisprscan-zebrafish-request.json \
  --output crisprscan-report.json
```

The report retains the request checksum, implementation version and commit, coefficient checksum, declared biological context, score labels, and interpretation warnings. A guide-activity score is **not** a phenotype prediction, editing probability, off-target assessment, or safety claim.

## Architecture / Research Philosophy

GeneImpact AI is designed around a simple principle: **a prediction is useful only when its evidence, applicability, uncertainty, and provenance remain visible**.

```mermaid
flowchart LR
    A["Declared study context"] --> B["Species & assembly validation"]
    B --> C["Versioned evidence adapters"]
    C --> D["Predictor applicability gates"]
    D --> E["Biostatistical evaluation"]
    E --> F["Audit-ready report"]
    F --> G["Human, welfare & institutional review"]
```

### Evidence-aware by construction

The system keeps empirical observations, model outputs, and biological hypotheses distinct. External predictors are accepted only through explicit adapter contracts, then checked against their supported species, edit class, experimental context, implementation version, and evidence reference. Out-of-domain outputs may remain visible for audit purposes, but they are not promoted to applicable evidence.

### Biostatistical rigor

GeneImpact AI favors bounded scores, held-out evaluation, gene-disjoint splits, Brier score, expected calibration error, Recall@K, and explicit transfer-evidence labels over unqualified accuracy claims. High-consequence or high-uncertainty signals cannot be averaged away by reassuring values elsewhere.

### Reproducibility and bounded uncertainty

Assessments are tied to exact reference assemblies, strains or isolates, evidence snapshots, model versions, and content hashes. The project treats uncertainty as a reportable result—not a defect to conceal—and requires prospective experimental validation before any operational interpretation.

> [!CAUTION]
> GeneImpact AI is research decision-support software. It does not authorize an animal genome edit, establish safety, replace experimental validation, or substitute for ethics, biosafety, veterinary, animal-welfare, or regulatory review.

For methodology and governance details, see the [research protocol](docs/research-protocol.md), [model card](docs/model-card.md), [multi-species registry](docs/multispecies-registry.md), [validation plan](docs/validation-plan.md), and [data-governance policy](docs/data-governance.md).

## Contributing

Contributions that improve scientific validity, reproducibility, documentation, test coverage, accessibility, or bounded species-specific validation are welcome.

1. Fork the repository and create a focused branch.
2. Add or update tests for behavioral changes.
3. Run `pytest` locally.
4. Document the biological domain, evidence source, assumptions, and limitations of any new predictor or adapter.
5. Open a pull request with a concise scientific and technical rationale.

Please do not commit raw sequencing files, restricted study data, facility records, animal identifiers, or third-party model assets without confirmed redistribution rights.

## License

GeneImpact AI is released under the [MIT License](LICENSE). Third-party datasets, publications, services, and model integrations may have separate terms; review the [third-party notices](docs/third-party-notices.md) before redistribution or operational use.
