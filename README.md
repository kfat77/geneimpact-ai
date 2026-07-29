# GeneImpact AI

> An evidence-first research framework for testing how single genes and pre-specified gene combinations are associated with human phenotypes.

## Purpose

GeneImpact AI is a research-use-only framework for testing **phenotype associations**, **gene-gene interactions**, and **uncertainty** from curated genomic and biomedical evidence. It is not a diagnostic system, a clinical decision-support tool, or a substitute for genetic counselling.

The research target is deliberately narrow and falsifiable: for each pre-registered phenotype, determine whether a gene-level signal or a biologically motivated gene pair improves replicated, calibrated association estimates beyond a covariate-adjusted baseline. See [research aims](docs/research-aims.md) and the [pre-registration template](docs/preregistration-template.md).

The first release focuses on a transparent pipeline:

1. Ingest versioned, consented genotype/variant and phenotype data.
2. Normalize variants to a declared genome build and annotation version.
3. Create single-gene and *pre-specified* multi-gene interaction features.
4. Fit a covariate-adjusted baseline, then a calibrated interaction model.
5. Validate in held-out and independent cohorts.
6. Return an effect estimate with uncertainty, provenance, and an evidence grade.

## Scientific guardrails

- **Association is not causation.** Outputs distinguish statistical association, replicated evidence, and causal evidence.
- **No deterministic language.** Results must report effect sizes, confidence intervals, cohort context, and limitations.
- **Pre-registration before discovery.** Phenotype, interaction search space, covariates, and success criteria are locked before access to outcome results.
- **Multiplicity-aware inference.** Interaction claims must use the pre-declared correction method and report the full number of tested hypotheses.
- **Population-aware evaluation.** Performance is reported by ancestry/cohort strata; models are not released if material disparities are unresolved.
- **No clinical use.** Any future clinical workflow needs prospective validation, regulatory review, and qualified human oversight.
- **Privacy by design.** Never commit raw genotype files, identifiers, or re-identifiable phenotype records.

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
from geneimpact.interactions import rank_interactions

rank_interactions(
    gene_scores={"BRCA1": 0.8, "BRCA2": 0.6, "TP53": 0.4},
    evidence={frozenset({"BRCA1", "BRCA2"}): 0.9},
)
```

The returned scores are prioritization signals for downstream research. They are not health-risk predictions.

## Roadmap

- [x] Research scope, governance, and model-card templates
- [x] Transparent, evidence-gated interaction-prioritization baseline
- [x] Pre-registration and evidence-grading templates
- [ ] Curated public-data connectors (ClinVar, GWAS Catalog, Ensembl)
- [ ] Reproducible cohort-split, confounder adjustment, and calibration workflow
- [ ] External replication benchmarks
- [ ] Secure, consent-aware data access layer

## License

MIT. See [LICENSE](LICENSE).
