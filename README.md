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
```

## Quick start

```bash
python -m pip install -e ".[dev]"
pytest
```

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
- [ ] Versioned reference-genome and annotation adapters
- [ ] Reproducible species/edit-class split and calibration workflow
- [ ] Independent laboratory replication benchmarks
- [ ] Restricted animal-study data access layer

## License

MIT. See [LICENSE](LICENSE).
