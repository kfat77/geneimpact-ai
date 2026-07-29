# CRISPRscan zebrafish activity adapter

GeneImpact AI implements the published CRISPRscan linear score for ranking
canonical SpCas9 guides in its declared zebrafish embryo domain.

## Declared domain

The current adapter accepts only:

- species profile: `zebrafish`
- registered reference: GRCz12tu / `GCF_049306965.2` / Tuebingen
- nuclease: canonical-PAM SpCas9
- guide expression: T7 in-vitro transcription
- developmental context: zebrafish embryo
- input: 35 bases in `[6 upstream][20 spacer][NGG PAM][6 downstream]` order

The model was trained on in-vivo zebrafish embryo mutagenesis data, principally
in the TU background, using the older Zv9 assembly context. The score itself is
sequence based, but every submitted context must be rechecked against the
current registered assembly and the study animals' variants.

The public CRISPRscan service exposes additional species and assemblies. That
does not establish species-specific model validity, so this adapter rejects
mouse, rat, fruit fly, and macaque requests.

## Run

```bash
python -m geneimpact score-crisprscan \
  --input examples/crisprscan-zebrafish-request.json \
  --output crisprscan-report.json
```

Raw 35-nt contexts are accepted in the request but are not repeated in the
report; the report stores their SHA-256 hashes and researcher-supplied safe
identifiers.

## Reproducibility

The 91 published features and intercept are locked to the coefficient file in
`crisprScore` 1.15.3 commit
`cbd6f9f60dc7fb50d14b90485b9561d582caf21e`. The source coefficient-file
SHA-256 is
`6e3f1bbfd58e5426651a15cfd0db6ac2094e0a93158dc51639b5929fc9ced5a4`.

The documented `crisprScore` example is used as a numerical oracle in the test
suite. The GeneImpact implementation returns `0.5691911213797` for that
35-base context.

## Interpretation limits

The original paper reported score thresholds above 0.55 for efficient guides
and above 0.70 for highly efficient guides in its validation setting.
GeneImpact records these as published threshold labels, not universal
categories.

A CRISPRscan score is not a calibrated editing probability. It does not predict
repair products, tissue mosaicism, phenotype, off-target cutting, welfare
effects, or safety. Results outside the declared expression, nuclease,
developmental, strain, or species context require separate validation.
