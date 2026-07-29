# Cynomolgus embryo base-editing transfer benchmark

GeneImpact AI v0.16.1 provides a checksum-pinned workflow for testing an
external base-editing model against Zhang et al. 2020 cynomolgus embryo data.
It is a retrospective transfer evaluator, not a macaque predictor.

## Qualified scope

- species: *Macaca fascicularis*;
- population: source colony; geographic ancestry was not reported;
- stage: zygotes collected after ICSI and cultured for at least three days;
- delivery: base-editor mRNA and T7 sgRNA cytoplasmic microinjection;
- editors: BE3, ABE7.10, and SaKKH-BE3;
- multiplex level: one, two, or three co-injected guides;
- assay: pooled Sanger-clone counts for intended conversion at named target
  bases;
- source assembly: legacy `GCF_000364345.1`, not the registered
  T2T-MFA8v1.1 assembly.

The qualified source contains 66 candidate-site rows, including 11 on-target
sites, and 30 target-base/context records. Those records represent 273
embryo-by-target-base observations with 8,296 repeated clone denominators.
They are not 273 independent biological samples.

## Download the publisher files

Download:

1. Supplementary Data 1, `41467_2020_16173_MOESM1_ESM.xlsx`;
2. Source Data, `41467_2020_16173_MOESM11_ESM.xlsx`.

The workflow rejects altered bytes. The repository does not redistribute
either workbook.

## Prepare a prediction submission

```bash
python -m geneimpact prepare-cynomolgus-base-editing-transfer \
  --target-sites 41467_2020_16173_MOESM1_ESM.xlsx \
  --source-data 41467_2020_16173_MOESM11_ESM.xlsx \
  --output cynomolgus-predictions.json
```

The template identifies every record by editor, multiplex context, target
base, target-site ID, target-sequence SHA-256, and sequence length. It does not
emit the source target sequence or observed editing label.

Run the external model using the publisher target sequence, then fill:

- immutable model name and version or commit;
- `submitted_code_revision`: the caller-supplied model-run code revision;
- `score_semantics`: `ranking_score` or `expected_edit_fraction`;
- training-overlap status;
- model/run evidence reference;
- every `predicted_score`.

Known training overlap is rejected. An `unknown` overlap declaration remains
visible and does not become an independently verified test.
The report marks `code_revision_verified=false`: the caller-supplied revision
is retained for lineage but is not cryptographically checked by this workflow.
It also records the running GeneImpact AI package version and an automatic
SHA-256 of the evaluator module (`evaluator_code_revision_status=
module_source_sha256`), which are the authoritative evaluator provenance fields.

## Evaluate the submission

```bash
python -m geneimpact validate-cynomolgus-base-editing-transfer \
  --target-sites 41467_2020_16173_MOESM1_ESM.xlsx \
  --source-data 41467_2020_16173_MOESM11_ESM.xlsx \
  --predictions cynomolgus-predictions-filled.json \
  --output cynomolgus-transfer-report.json
```

For any score type, the report computes pairwise ranking agreement only among
records sharing both an editor and one injection context. The six injection
contexts therefore form seven comparison strata because one context combines
SaKKH-BE3 and ABE7.10. This prevents an arbitrary BE3 score from being ranked
directly against ABE7.10 or against a different multiplex delivery.
Prediction ties receive half credit, observation ties are excluded, and the
report exposes candidate-pair count, eligible-pair count, weighted concordance,
strict-pair count, and prediction coverage so ties cannot be hidden by a point
estimate.

If the submitted values explicitly represent expected edit fractions, the
report additionally computes MAE and RMSE against the pooled intended
base-conversion fraction. These are descriptive transfer errors, not
calibration metrics.

The report records the fixed external-transfer split identifier, source and
target assembly accessions, `liftover_status=not_performed`, and an explicit
confidence-interval status. The publisher target-sequence record is checked,
but the sequence is not independently reconstructed from the source assembly;
the report therefore keeps `target_sequence_verified_on_source_assembly=false`.
No interval is reported here because the source observations are clustered
within embryo batches and injection contexts.

## Interpretation limits

- Eleven target sites are insufficient for species-level training or
  calibration.
- Multiple target bases, clones, and targets from one embryo batch are
  dependent observations.
- The source result cannot be transferred automatically to another editor,
  delivery, stage, population, tissue, laboratory, or assay.
- No coordinate lift-over to T2T-MFA8v1.1 is performed.
- Candidate-site or single-fetus WGS evidence does not establish genome-wide
  off-target safety.
- The benchmark does not predict repair spectrum, embryo development,
  phenotype, welfare, live birth, or safety.

The machine-readable capability status is `usable_bounded_benchmark`, while
`predictive_adapter_available` remains false.
