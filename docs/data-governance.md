# Data governance

## Allowed inputs

Only data with a documented legal basis, participant consent (where applicable), and an approved data-use agreement may enter the pipeline.

## Prohibited repository content

Raw VCF/BAM/CRAM files, participant identifiers, dates of birth, addresses, free-text clinical notes, and re-identifiable phenotype records must never be committed.

## Operational controls

- Store raw data in an approved controlled-access environment.
- Use pseudonymous study identifiers and least-privilege access.
- Keep an immutable lineage record for dataset, genome build, annotation, feature set, code revision, and model artifact.
- Complete security and privacy review before any data connector is enabled.
