# Unified research dossier

The `dossier` command is the preferred researcher-facing entry point. It
combines study context, target-gene evidence, candidate gene-pair interactions,
predictor outputs, applicability checks, evidence gaps, and integrity metadata
in one JSON report.

It supports all registered profiles:

- mouse
- rat
- zebrafish
- fruit fly
- rhesus macaque
- cynomolgus macaque

Generic `monkey` or `macaque` values are not accepted. The exact species,
reference assembly, and study strain must be declared.

## Run

```bash
python -m geneimpact dossier \
  examples/dossier-zebrafish-request.json \
  --output research-dossier.json

python -m geneimpact verify-dossier research-dossier.json

python -m geneimpact dossier \
  examples/dossier-mouse-indelphi-request.json \
  --output mouse-research-dossier.json
```

The example is a synthetic format demonstration. Its placeholder evidence
references, target identifiers, and evidence snapshot hash must be replaced
before research use.

## Required study declaration

The request binds the result to:

- a safe study identifier;
- one exact species profile;
- study strain or breed;
- current registered genome build and assembly accession;
- edit class;
- delivery and developmental contexts;
- SHA-256 of the curated evidence snapshot;
- one or more target genes;
- intended outcomes and animal-welfare endpoints.

Each target gene carries a bounded evidence signal and a versioned evidence
reference. Gene identifiers are currently researcher-declared; the report
flags that automatic resolution against a versioned annotation is not yet
implemented.

## Multiple-gene hypotheses

For two or more target genes, the request may declare independently curated
pair evidence. The report ranks every possible pair using the geometric mean
of the two declared gene signals multiplied by the pair evidence weight.
Missing pair evidence produces a zero priority rather than an invented
interaction.

These are hypothesis priorities, not predictions of organism-level phenotype
or epistasis.

## Predictor integration

The dossier can currently include:

- locally computed CRISPRscan scores in the declared zebrafish domain;
- externally generated, version-locked CRISPRitz target-search results for
  every registered reference genome;
- externally generated BE-Hive efficiency and bystander results in mouse mES;
- externally executed, version-locked inDelphi mESC repair-outcome results for
  mouse knockout studies;
- official FlyRNAi Housden results for fruit-fly S2R+ cell-culture knockout
  studies;
- validated generic concern outputs under the predictor adapter contract.

For CRISPRitz, `targets_file` must be a relative path beneath the request file's
directory. Absolute paths and directory traversal are rejected. The imported
file checksum is retained in the report.

For inDelphi, `result_files` contains 1–100 relative JSON paths beneath the
request directory. Every result is checked against the dossier's assembly,
delivery context, and developmental context. The adapter recalculates the
reported distribution statistics and includes pinned external mouse-embryo
transfer evidence. See the [inDelphi adapter](indelphi-adapter.md).

For Housden, `result_files` contains 1–100 relative JSON envelopes and
`source_response_files` contains the corresponding retained official XLS
responses. Every response checksum, guide row, score, sequence, fruit-fly
assembly, U6 guide-expression context, and exact S2R+ cell-culture context is
checked. An embryo or in-vivo dossier cannot include the score. See the
[Housden adapter](housden-adapter.md).

Every capability relevant to the selected species and edit class is listed as
`included`, `available_not_run`, `not_integrated`, `out_of_domain`, or
`irrelevant_edit_class`. A missing available predictor cannot disappear from
the report.

## Integrity and interpretation

The report stores:

- canonical request SHA-256;
- original request-file SHA-256;
- predictor and attachment checksums;
- model and source commits;
- a report content SHA-256.

`verify-dossier` detects content changes after generation. The dossier is
currently unsigned, so the hash is not proof of author identity or
authenticity.

The output also records incomplete predictor coverage, missing interaction
evidence, strain/reference differences, and the absence of prospective
empirical validation. It never converts those gaps into an edit-safety
approval.
