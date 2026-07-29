# GeneImpact AI

> An evidence-first research framework for predicting and documenting the plausible consequences of proposed animal genome edits.

## Purpose

GeneImpact AI is a research-use-only decision-support framework for assessing **on-target disruption**, **candidate off-target effects**, **biological-network context**, and **animal-welfare relevance** before an animal genome-editing experiment. It prioritizes questions for validation; it does not approve an edit, design an experiment, or guarantee an outcome.

The research target is deliberately narrow and falsifiable: for each declared species, edit class, target locus, and intended trait, determine whether the system ranks observed unintended consequences and welfare-relevant outcomes better than a documented baseline. See [research aims](docs/research-aims.md) and the [pre-registration template](docs/preregistration-template.md).

The first release focuses on a transparent pipeline:

1. Register the animal species, genome build, intended edit outcome, and welfare endpoints.
2. Ingest versioned reference-genome, functional-annotation, and prior-outcome evidence.
3. Generate an evidence trace for plausible on-target, off-target, and network-level consequences.
4. Produce a calibrated risk *tier* and uncertainty record, not a pass/fail answer.
5. Evaluate against held-out historical outcomes and an independent laboratory or dataset.
6. Require human and animal-welfare review for every high-uncertainty or high-consequence result.

Every assessment is also bound to a species, strain/breed, genome build, edit class, evidence snapshot, and model version. This prevents a result from being reused outside its demonstrated setting.

Registered contexts now cover mouse, rat, zebrafish, fruit fly, rhesus
macaque, and cynomolgus macaque. Registration provides strict context and
assembly auditing; predictor validation remains separately reported for each
species. See the [multi-species registry](docs/multispecies-registry.md).

## Scientific guardrails

- **Predictions are not guarantees.** Outputs distinguish biological plausibility, replicated empirical evidence, and uncertainty.
- **No deterministic language.** Results must report effect sizes, confidence intervals, cohort context, and limitations.
- **Pre-registration before evaluation.** Species, edit class, endpoints, data split, and success criteria are locked before model assessment.
- **Species-aware validity.** Results are reported separately by species, strain/breed, genome build, and edit class; no cross-species generalization is assumed.
- **Animal welfare is a first-class endpoint.** A high predicted welfare consequence or high uncertainty cannot be collapsed into a low-risk recommendation.
- **No autonomous use.** Any future operational use needs institution-specific ethics, biosafety, veterinary, and regulatory review.
- **Data stewardship.** Never commit raw sequencing files, facility records, animal identifiers, or restricted study data.

## Repository layout

```
docs/                  Research protocol, data governance, and model card
src/geneimpact/        Feature and interaction-scoring library
tests/                 Unit tests
examples/              Runnable assessment request example
```

## Quick start

```bash
python -m pip install -e ".[dev]"
pytest
```

## Run an assessment

```bash
python -m geneimpact dossier examples/dossier-zebrafish-request.json --output research-dossier.json
python -m geneimpact verify-dossier research-dossier.json
python -m geneimpact assess examples/assessment-request.json --output assessment-report.json
python -m geneimpact source-check --species mouse
python -m geneimpact snapshot-mgi --report all-phenotypes --output-dir data/mgi
python -m geneimpact normalize-mgi --input data/mgi/MGI_PhenotypicAllele.rpt --output data/mgi/endonuclease-alleles.jsonl
python -m geneimpact impc-gene --gene Prkdc --output data/impc/Prkdc-significant.json
python -m geneimpact benchmark-mgi --input data/mgi/endonuclease-alleles.jsonl --output-dir data/benchmarks/mgi-v1
python -m geneimpact evaluate-baseline --benchmark-dir data/benchmarks/mgi-v1 --k 5
python -m geneimpact benchmark-impc --gene Prkdc --gene Kit --output data/benchmarks/impc-validation.jsonl
python -m geneimpact calibrate-impc --calibration data/benchmarks/impc-calibration-v1.jsonl --test data/benchmarks/impc-test-v1.jsonl --output data/benchmarks/impc-calibration-report-v1.json
python -m geneimpact import-behive-efficiency --input examples/behive-efficiency-import.json --output behive-audit.json
python -m geneimpact import-behive-bystander --input examples/behive-bystander-import.json --output behive-bystander-audit.json
python -m geneimpact import-indelphi --input examples/indelphi-mouse-result.json --output indelphi-audit.json
python -m geneimpact dossier examples/dossier-mouse-indelphi-request.json --output mouse-research-dossier.json
python -m geneimpact import-housden --input examples/housden-fruit-fly-result.json --source-response downloaded-flyrnai-response.xls --output housden-audit.json
python -m geneimpact capabilities --species zebrafish
python -m geneimpact readiness --all
python -m geneimpact import-crispritz --metadata examples/crispritz-rat-metadata.json --targets examples/crispritz-synthetic.targets.txt --output crispritz-audit.json
python -m geneimpact score-crisprscan --input examples/crisprscan-zebrafish-request.json --output crisprscan-report.json
python -m geneimpact validate-crisprscan-transfer --input data/benchmarks/crisprscan-nhgri1-2022.json --output crisprscan-transfer-report.json
```

Read the [researcher guide](docs/researcher-guide.md) before using the tool with
study evidence. The [unified dossier guide](docs/research-dossier.md) documents
the preferred one-request workflow. The
[BE-Hive adapter guide](docs/behive-adapter.md) documents
the first real model integration and its deliberately narrow mES scope. The
[CRISPRitz adapter guide](docs/crispritz-adapter.md) documents cross-species
reference-search auditing and its interpretation limits. The
[CRISPRscan adapter guide](docs/crisprscan-adapter.md) defines the narrower
zebrafish embryo activity-scoring domain. The
[inDelphi adapter guide](docs/indelphi-adapter.md) documents the licensed
external-result workflow and mouse-embryo transfer evidence.
The [Housden adapter guide](docs/housden-adapter.md) documents the
external-result workflow for fruit-fly S2R+ cell guide ranking and its strict
no-in-vivo-extrapolation rule.

## Initial API concept

```python
from geneimpact.edit_assessment import EditEvidence, assess_edit

assess_edit(
    EditEvidence(
        on_target_uncertainty=0.2,
        off_target_evidence=0.7,
        network_impact_evidence=0.5,
        welfare_relevance=0.8,
    )
)
```

The result is an evidence trace and review tier for downstream research. It is not an authorization to edit an animal.

## Roadmap

- [x] Research scope, governance, and model-card templates
- [x] Transparent evidence-to-review-tier baseline
- [x] Animal-edit pre-registration and evidence-grading templates
- [x] Audit-ready study context and held-out calibration metrics
- [x] Runnable JSON-to-report researcher workflow
- [x] Registered mouse/GRCm39 profile and live authoritative metadata adapter
- [x] Current NCBI reference profiles for rat, zebrafish, fruit fly, rhesus macaque, and cynomolgus macaque
- [x] Ambiguity rejection and strain/isolate warnings across registered species
- [x] Code-level species × predictor evidence-status matrix
- [x] Version-locked CRISPRitz cross-species reference-search audit adapter
- [x] Version-locked CRISPRscan zebrafish embryo activity scorer
- [x] Independent 50-guide zebrafish RNP transfer benchmark with honest domain labeling
- [x] Unified multi-species research dossier with target-gene interactions and integrity verification
- [x] Version-locked inDelphi mESC repair-outcome import with external mouse-embryo transfer evidence
- [x] Machine-readable evidence qualification for every registered species
- [x] Official-service Housden result import for fruit-fly S2R+ cell culture
- [x] Checksum-bearing MGI phenotype snapshot adapter
- [x] Streaming MGI genome-edit evidence normalization
- [x] Bounded IMPC significant gene-phenotype query adapter
- [x] Leakage-aware, gene-grouped MGI positive-association benchmark
- [x] Manifest-bound Recall@K baseline for unseen genes
- [x] Bounded independent IMPC validation builder with tested outcomes
- [x] Gene-disjoint IMPC Brier/ECE calibration baseline
- [x] Version-locked BE-Hive mES efficiency import and audit adapter
- [x] Version-locked BE-Hive mES bystander-outcome import and audit adapter
- [x] Leakage- and domain-gated BE-Hive independent validation evaluator
- [ ] IMPC phenotype snapshot adapter
- [ ] Isolated, license-reviewed BE-Hive execution environment
- [ ] Held-out mES efficiency comparison against an independent experimental dataset
- [ ] Reproducible species/edit-class split and calibration workflow
- [ ] Independent laboratory replication benchmarks
- [ ] Prospective locus- and laboratory-specific mESC-to-embryo validation
- [ ] Restricted animal-study data access layer

## License

MIT. See [LICENSE](LICENSE) and [third-party notices](docs/third-party-notices.md).
