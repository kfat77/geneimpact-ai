# BE-Hive efficiency adapter

## What this integration does

GeneImpact AI can validate and normalize an efficiency prediction produced by
the official BE-Hive efficiency implementation. The adapter records the exact
model commit, editor, cell type, sequence hash, raw logit score, calibration
parameters, declared biological scope, and limitations.

This is an **import adapter**, not a reimplementation of BE-Hive. It does not
design a guide or run an edit.

## Verified upstream scope

The adapter is locked to:

- repository: <https://github.com/maxwshen/be_predict_efficiency>
- reviewed commit:
  [`fbd495910d6c95b24081649015d6257a8badc9d7`](https://github.com/maxwshen/be_predict_efficiency/tree/fbd495910d6c95b24081649015d6257a8badc9d7)
- cell type: `mES` (mouse embryonic stem cells)
- task: base-editing efficiency
- input: 50 DNA nucleotides spanning positions -19 through +30 relative to
  the guide, with an NGG PAM at positions +21 through +23

The model's declared mES scope is not evidence of validity for a mouse strain,
tissue, developmental stage, germline transmission, mosaicism, welfare
outcome, or whole-animal phenotype.

## Why execution is isolated

The reviewed upstream implementation declares Python 3.7-era dependencies,
including scikit-learn 0.20.3 and Biopython 1.73. Its GitHub repository does
not expose a standard software license. GeneImpact AI therefore does not copy,
vendor, install, or redistribute the upstream code or model parameters.

Researchers must execute an appropriately authorized copy in a separately
managed environment and import its machine-readable result. The audit record
states that the result was externally reported and was not independently
recomputed.

## Import format

Use [`examples/behive-efficiency-import.json`](../examples/behive-efficiency-import.json)
as the input contract:

```json
{
  "request": {
    "sequence": "50-nt DNA sequence",
    "base_editor": "BE4",
    "cell_type": "mES",
    "model_commit": "fbd495910d6c95b24081649015d6257a8badc9d7"
  },
  "raw_output": {
    "Predicted logit score": 0.4
  }
}
```

Normalize it with:

```bash
python -m geneimpact import-behive-efficiency \
  --input examples/behive-efficiency-import.json \
  --output behive-audit.json
```

The output contains a SHA-256 digest of the sequence, not the raw sequence.

## Calibration rule

BE-Hive returns a model logit. A fraction of edited reads is accepted only
when both `calibration_mean` and `calibration_std` are explicitly declared in
the request. GeneImpact AI recomputes:

```text
fraction = logistic(logit * calibration_std + calibration_mean)
```

If the upstream result also reports a fraction, the importer checks that it
matches the declared transformation. Without both calibration parameters, a
reported fraction is omitted and the logit is retained.

Calibration parameters should come from a pre-specified, representative
experimental domain. They do not convert an in-vitro efficiency estimate into
a probability of organism-level benefit, harm, or phenotypic effect.

## Unified assessment report

An assessment request may include a list named `behive_efficiency_outputs`,
where each item follows the import format above. Results appear under
`model_predictions`, separate from `predictor_outputs`.

This separation is intentional: editing efficiency is not a concern score.
The transparent concern tier remains based on the submitted evidence record,
and the BE-Hive result is retained as task-specific supporting evidence with
its applicability label.
