# Validation plan

## Objective

Evaluate whether concern scores help teams identify observed unintended consequences and welfare-relevant outcomes earlier than a documented evidence-only baseline.

## Dataset partitioning

Split by study or laboratory first, never by individual record alone. Hold out entire edit classes where practical. Keep species, strain/breed, genome build, and edit-class labels intact so transfer limits can be measured rather than hidden.

## Required metrics

- Recall of observed high-concern outcomes at the chosen review threshold.
- False-positive rate and review burden at that threshold.
- Brier score and expected calibration error on held-out outcomes.
- All metrics by species, strain/breed, genome build, edit class, and laboratory/dataset.

## Release gate

Do not release an operational recommendation unless its held-out calibration, review burden, and subgroup results meet pre-registered criteria and have been independently replicated. A missing subgroup result is reported as an evidence gap, not inferred from the aggregate.

## Audit artifact

Every evaluation result must retain: study context, evidence snapshot identifier, model version, code revision, data-split identifier, metrics, and an explanation of exclusions.
