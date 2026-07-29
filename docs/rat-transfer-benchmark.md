# Rat guide-activity transfer benchmark

GeneImpact AI can now evaluate an externally generated SpCas9 guide-activity
score against a small, pinned rat in-vivo dataset. This is a retrospective
**transfer benchmark**, not a rat predictor, a probability calibration, or an
edit-safety test.

The source is Anderson et al. (2018), which reports mouse and rat embryo-editing
projects. GeneImpact uses the 14 rat guide labels that can be mapped uniquely
between Supplementary Tables 1 and 5. Those records represent 186 reported
animal-or-embryo observations. Two additional `Usp30` labels are excluded
because the guide identity is ambiguous.

## Source files

Download the source workbooks directly from the publisher:

- [Supplementary Table 1](https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-018-0011-5/MediaObjects/41592_2018_11_MOESM3_ESM.xlsx)
- [Supplementary Table 5](https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-018-0011-5/MediaObjects/41592_2018_11_MOESM7_ESM.xlsx)
- [Anderson et al. article](https://doi.org/10.1038/s41592-018-0011-5)

The importer requires these exact SHA-256 values:

| File | Expected SHA-256 |
|---|---|
| Supplementary Table 1 | `1efe7e4bb49a3feb9b5757b3c42989d0e95cd1da621ae4fd12c8f0aa3f10fdc8` |
| Supplementary Table 5 | `7dfd554af5b9677723970a892c39ac2529a09894cb178f3252121fd18bd2e0c8` |

The workbooks are not copied into this repository. Download access does not
establish broad reuse or redistribution rights; no explicit dataset reuse
licence was located during the evidence review.

## Researcher workflow

First generate a sequence-redacted submission template:

```bash
python -m geneimpact prepare-rat-guide-transfer \
  --table1 41592_2018_11_MOESM3_ESM.xlsx \
  --table5 41592_2018_11_MOESM7_ESM.xlsx \
  --output rat-predictions.json
```

The command verifies both source files and emits 14 target names with SHA-256
identifiers and lengths for both the published design spacer and the actual
5′-G guide. The source contains 19/20nt design spacers and 20/21nt actual
guides, so those two sequence meanings must not be conflated. The template
does not copy raw guide sequences. Edit only:

- predictor `name` and immutable `version` or commit;
- `score_semantics`: `ranking_score` or `expected_edit_fraction`;
- `prediction_target`: keep the fixed value
  `mean_on_target_edit_fraction`;
- `sequence_basis`: `design_sequence` or `actual_guide_sequence`, matching the
  sequence actually passed to the model;
- `training_overlap_status`: `declared_no_overlap` or `unknown`; a known
  overlap is rejected because it is not an external transfer test;
- the model or run `evidence_reference`;
- every `predicted_score`.

Do not alter target names, sequence lengths, or hashes. Scores must increase
with predicted activity. An `expected_edit_fraction` must be between 0 and 1
and must mean the expected mean on-target edit fraction—not a probability of
crossing a threshold or another event.

Then run the evaluation:

```bash
python -m geneimpact validate-rat-guide-transfer \
  --table1 41592_2018_11_MOESM3_ESM.xlsx \
  --table5 41592_2018_11_MOESM7_ESM.xlsx \
  --predictions rat-predictions.json \
  --output rat-transfer-report.json
```

The report includes source checksums, source and current assembly identities,
training-overlap status, Pearson correlation with a Fisher 95% interval,
Spearman rank correlation, and—for explicit expected-edit-fraction inputs
only—descriptive MAE and RMSE. The report stores aggregate metrics and a
checksum of the prediction submission; it does not serialize per-guide source
outcomes or counts. The training-overlap field is retained as a submitter
declaration, but the report does not mark independence as verified. That
requires a separate, reproducible audit of the model's actual training
sequences and projects.

## Interpretation limits

- The source uses legacy `rn5`; the current rat registry uses GRCr8.
- The project-level rat strains are not fully disclosed.
- Fourteen selected, mostly high-activity guides are too few for training,
  stable calibration, or broad sequence generalization.
- Animal-or-embryo counts are not independent guide-level replicates.
- Candidate off-target sites in the source were selected by an existing
  algorithm; this benchmark does not measure genome-wide off-target recall.
- Correlation does not predict phenotype, multi-gene interactions, structural
  variants, animal welfare, safety, or another strain, delivery system,
  developmental stage, nuclease, or edit class.

The broader eight-source evidence comparison and data-acquisition requirements
are documented in [rat predictor evidence](rat-predictor-evidence.md).
