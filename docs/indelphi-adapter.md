# inDelphi mouse repair-outcome adapter

GeneImpact AI validates and normalizes results produced by the official
inDelphi mESC model. It does not redistribute or silently execute the upstream
code or model files.

## Exact supported scope

The adapter is locked to:

- official repository:
  <https://github.com/maxwshen/inDelphi-model>
- commit:
  `9ab67ca53ebb91e49aeb4530ec1e999ee9827ca1`
- model artifact directory: `model-sklearn-0.18.1`
- scikit-learn: `0.18.1`
- cell model: `mESC`
- species profile: `mouse`
- nuclease: canonical SpCas9
- task: distribution of 1-bp insertions and 1–59-bp deletions after a
  double-strand break

The source is under a limited copyright license whose terms name eligible
academic users and US government research institutions for non-commercial
research use. Commercial or industrially sponsored use requires a separate
upstream agreement. Researchers must review their own eligibility. GeneImpact
AI therefore provides an external-result adapter and does not vendor the
source, models, or a derived execution image.

## Import

An authorized external environment must emit the contract demonstrated by
`examples/indelphi-mouse-result.json`. The file declares the exact source
commit, model folder, model-bundle checksum, runtime versions, sequence,
cutsite, study contexts, official summary statistics, and the full modeled
outcome table.

```bash
python -m geneimpact import-indelphi \
  --input examples/indelphi-mouse-result.json \
  --output indelphi-audit.json
```

The repository example is synthetic and exists only to test the contract. Its
outcomes are not an inDelphi prediction.

The importer:

- requires at least 60 bases on both sides of the cut;
- recomputes the outcome sum, insertion/deletion class sums, reading-frame
  sums, entropy-derived precision, highest outcomes, and expected indel
  length;
- rejects duplicate outcomes, altered statistics, mismatched sequences,
  mismatched cuts, and unpinned runtimes;
- retains DNA and optional genotype SHA-256 values rather than raw sequences;
- reports only the 25 highest-frequency outcomes while binding the full
  external result by checksum.

Frameshift fields are exposed only when the request declares
`target_context: coding_sequence`. They are conditional repair-product
frequencies, not probabilities of knockout, phenotype, or successful animal
production.

## Unified dossier

One mouse dossier can include up to 100 independently generated guide results:

```json
{
  "predictors": {
    "indelphi": {
      "result_files": [
        "guide-01.json",
        "guide-02.json"
      ]
    }
  }
}
```

Every path must remain beneath the dossier request directory. Genome build,
assembly accession, delivery context, and developmental context must match the
study declaration.

## External mouse-embryo evidence

Lkhagvadorj et al. (2026) compared the mESC-trained inDelphi predictions with
HiFi-Cas9 RNP editing outcomes in C57BL/6JJmsSlc embryos. GeneImpact checked
the article's supplementary workbook, bound it by SHA-256, and reproduced the
aggregate Tyr comparison:

- 14 guides;
- 1,182 prediction/observation outcome pairs;
- overall Pearson correlation `0.6375914615783601`;
- 1 very-strong, 5 strong, 3 moderate, and 5 weak per-guide correlations.

The source workbook is not redistributed because it is marked CC BY-NC-ND
4.0. Only source identifiers, its checksum, and aggregate facts are retained.

This is useful but insufficient transfer evidence. The same study found
experimentally measured mESC outcome profiles were more concordant with
blastocysts than inDelphi alone. GeneImpact therefore requires the audit record
to say that prospective mESC validation remains necessary before embryo work.

## What it cannot predict

inDelphi does not predict:

- cutting or editing efficiency;
- wild-type read frequency;
- large deletions, complex rearrangements, integrations, or chromosome loss;
- endogenous-chromatin, strain, delivery, mosaic, germline, or laboratory
  effects;
- organism phenotype, welfare consequence, or safety.

Use repair-outcome prioritization together with off-target search, empirical
mESC assays, short- and long-read genotyping plans, and institutional animal,
biosafety, and veterinary review.
