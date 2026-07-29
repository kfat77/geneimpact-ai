# Research aims and success criteria

## Overall objective

Build a reproducible framework to assess whether observed associations between genetic variation—at one gene or across a pre-specified gene pair—and a human phenotype are robust, calibrated, and independently replicable.

This project does **not** attempt to claim that a gene “determines” a human outcome. The expected result is a bounded estimate of association in a defined cohort and phenotype context.

## Aim 1 — Single-gene association baseline

For each pre-registered phenotype, fit an interpretable gene-level baseline adjusted for pre-declared confounders (at minimum: age or time scale where applicable, sex where applicable, cohort/study site, and population-structure variables). Report the effect estimate, 95% confidence interval, calibration, and a multiplicity-adjusted result.

**Success criterion:** the estimate is estimable, calibrated on a held-out split, and reported with all quality-control and cohort details. Statistical significance alone is not success.

## Aim 2 — Evidence-gated gene-pair interaction test

Test only gene pairs supported before outcome analysis by a curated evidence source or pre-declared biological network. Compare the interaction model with the single-gene baseline using an appropriate held-out metric and report the entire tested interaction set.

**Success criterion:** an improvement over the baseline is observed in a held-out cohort, survives the declared multiplicity control, and has a directionally consistent estimate in an independent cohort.

## Aim 3 — Generalizability and failure analysis

Evaluate calibration and uncertainty in every adequately represented cohort/ancestry stratum. Document missingness, selection bias, transfer failures, and negative results.

**Success criterion:** subgroup results and limitations are released alongside the aggregate result; unresolved material performance gaps block promotion beyond exploratory use.

## Decision rule

| Evidence level | Required evidence | Permitted wording |
| --- | --- | --- |
| Exploratory | One discovery analysis | “Candidate association” |
| Replicated | Held-out result plus independent cohort consistency | “Replicated association” |
| Causal support | Replication plus a separately documented causal design and assumptions | “Causal evidence consistent with…” |

No output may use “causes”, “determines”, or an individual-risk interpretation without a study-specific governance and validation review.
