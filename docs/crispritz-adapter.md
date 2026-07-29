# CRISPRitz reference-search adapter

GeneImpact AI can validate and summarize a CRISPRitz target-search result for
every registered reference assembly. This is a sequence-search capability, not
a species-calibrated cleavage predictor.

## Supported registered references

The metadata must select one exact profile from the multi-species registry:

- `mouse` — GRCm39 / C57BL/6J
- `rat` — GRCr8 / BN/NHsdMcwi
- `zebrafish` — GRCz12tu / Tuebingen
- `fruit_fly` — Release 6 plus ISO1 MT / ISO-1
- `rhesus_macaque` — T2T-MMU8v2.0 / MMU2019108-1
- `cynomolgus_macaque` — T2T-MFA8v1.1 / 582-1

The assembly accession, reference strain or isolate, FASTA checksum, PAM
definition, mismatch and bulge limits, and variant snapshot status are bound
into the audit report. A generic `monkey` profile is intentionally rejected
because rhesus and cynomolgus macaques use different assemblies.

## External execution boundary

The adapter is verified against CRISPRitz v2.7.0 commit
`24b893ecb0c2354d5c76697e116d2febe1ee6265`. Run CRISPRitz separately, after
reviewing its installation instructions and license, then import the resulting
`.targets.txt` file:

```bash
python -m geneimpact import-crispritz \
  --metadata examples/crispritz-rat-metadata.json \
  --targets examples/crispritz-synthetic.targets.txt \
  --output crispritz-audit.json
```

The examples are synthetic format fixtures. Their placeholder FASTA checksum
must be replaced with the checksum of the reference FASTA actually searched
before a research run is auditable.

GeneImpact AI does not vendor or execute CRISPRitz. The upstream project
declares AGPL availability for academic researchers and a separate commercial
license requirement; users must review the upstream terms for their use.

## Interpretation limits

The report counts all candidate genomic sites and exact sequence matches,
summarizes total sequence differences, and retains at most 100
sequence-hashed hits. An exact match is not automatically labelled as the
intended locus because the targets file alone does not prove locus identity.
The report does not retain raw guide or target sequences.

Candidate sites are not measured cleavage events, calibrated cleavage
probabilities, phenotype predictions, or evidence that an edit is safe. A
reference-only search can miss colony, stock, breed, or individual variation.
Set `variant_aware` only when the search actually used a versioned variant
snapshot, and include that snapshot's SHA-256 digest.

Empirical off-target assays and species-, tissue-, editor-, delivery-, and
developmental-context validation remain necessary.
