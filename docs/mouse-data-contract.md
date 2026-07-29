# Mouse data contract

Mouse is the first registered benchmark species. The initial profile is intentionally narrow:

- species: `mouse` / `mus_musculus`
- taxonomy ID: `10090`
- reference assembly: `GRCm39`
- INSDC assembly accession: `GCA_000001635.9`
- reference strain boundary: `C57BL/6J` (`C57BL/6` accepted as an input alias)

## Authoritative sources

The platform uses [Ensembl REST species metadata](https://rest.ensembl.org/documentation/info/species) to detect assembly or accession drift. Functional and positional lookups will use versioned Ensembl endpoints.

Mouse phenotype and allele evidence will be sourced from [Mouse Genome Informatics data reports](https://www.informatics.jax.org/downloads/reports/index.html) and its documented programmatic access. MGI exposes phenotype, allele, genotype, strain, ontology, and orthology reports. Each imported snapshot must retain its source URL, retrieval timestamp, release identifier when available, file checksum, and licensing/provenance note.

IMPC results are planned as an independent high-throughput phenotype source. They are not treated as interchangeable with MGI annotations; dataset and procedure identifiers must remain attached to every observation.

## Applicability rule

An assessment is rejected if its genome build does not match the registered profile. A non-reference strain is allowed only with an explicit warning that strain-specific validation is still required. No coordinate or annotation is silently lifted between assemblies.

## Live metadata check

```bash
python -m geneimpact source-check --species mouse
```

The command verifies the locally registered taxonomy ID, assembly name, and assembly accession against the current Ensembl REST response.

## Reproducible MGI snapshot

```bash
python -m geneimpact snapshot-mgi \
  --report all-phenotypes \
  --output-dir data/mgi
```

The download is accompanied by a JSON manifest containing the exact source URL, UTC retrieval time, byte count, and SHA-256 checksum. The `data/` directory is ignored by Git because public source files may be large and restricted study data must never be committed.
