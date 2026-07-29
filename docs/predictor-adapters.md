# Predictor adapter contract

GeneImpact AI does not treat a named prediction tool as universally valid. Each upstream result is submitted with its model version, source record, declared species, declared edit class, concern score, and confidence.

## Supported task labels

- `guide_activity`
- `off_target`
- `repair_outcome`
- `base_editing`
- `prime_editing`

## Scope rule

An upstream output is included in the report only as a **declared match** when both its supported species and edit class match the study context. Otherwise it is retained for audit but is explicitly labelled `out_of_scope` and excluded from scoring.

## Why this matters

Tools such as inDelphi, FORECasT, BE-Hive, and guide-design predictors target different biological questions and have different training domains. The adapter layer exposes those boundaries; it does not silently convert an external model output into a cross-species safety claim.
