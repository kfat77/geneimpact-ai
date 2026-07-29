# Researcher guide: running an assessment

## What this MVP does

The tool turns an approved, structured evidence summary into a consistent review tier and a JSON report that preserves the study context and model version. It is designed for research triage and recordkeeping—not for designing edits or approving experiments.

## Install and run

```bash
python -m pip install -e ".[dev]"
python -m geneimpact assess examples/assessment-request.json --output assessment-report.json
```

The report includes the declared applicability boundary, the controlling concern signal, the review tier, and a non-authorization notice.

## Request fields

`study_context` records where a result applies: `species`, `strain_or_breed`, `genome_build`, `edit_class`, and `evidence_snapshot`.

`evidence` uses four bounded values from 0 (no concern signal) to 1 (strong concern signal): `on_target_uncertainty`, `off_target_evidence`, `network_impact_evidence`, and `welfare_relevance`.

Scores must come from the team's documented evidence-curation procedure. The source data and scoring rationale belong in the evidence snapshot outside this repository when restricted.

## Review rules

- 0.00–0.39: standard review
- 0.40–0.69: enhanced review
- 0.70–1.00: high-concern review

The highest signal controls the tier. This is intentional: a potential welfare concern or large uncertainty must not be offset by lower scores elsewhere.

## Before operational deployment

Validate the score definitions and thresholds with held-out outcomes, an independent laboratory or dataset, and the relevant ethics, biosafety, veterinary, and regulatory processes. See [validation plan](validation-plan.md).
