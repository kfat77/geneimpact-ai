# Multi-species registry

GeneImpact AI treats each species, reference assembly, strain or isolate, and
edit class as a separate applicability domain. Registration means the system
can validate and audit that context. It does **not** mean that every predictor
has been validated for that species.

## Registered reference contexts

The following profiles were checked against the NCBI Datasets assembly reports
on 2026-07-29:

| Profile key | Scientific name | Current reference assembly | RefSeq accession | Reference strain or isolate |
|---|---|---|---|---|
| `mouse` | *Mus musculus* | GRCm39 | GCF_000001635.27 | C57BL/6J |
| `rat` | *Rattus norvegicus* | GRCr8 | GCF_036323735.1 | BN/NHsdMcwi |
| `zebrafish` | *Danio rerio* | GRCz12tu | GCF_049306965.2 | Tuebingen |
| `fruit_fly` | *Drosophila melanogaster* | Release 6 plus ISO1 MT | GCF_000001215.4 | ISO-1 |
| `rhesus_macaque` | *Macaca mulatta* | T2T-MMU8v2.0 | GCF_049350105.2 | MMU2019108-1 |
| `cynomolgus_macaque` | *Macaca fascicularis* | T2T-MFA8v1.1 | GCF_037993035.2 | 582-1 |

`dm6` and `BDGP6.54` are accepted aliases for the registered fruit-fly
coordinate system. Species names may be supplied using the profile key,
scientific name, selected English common names, or the Chinese aliases 小鼠,
大鼠, 斑马鱼, 果蝇, 恒河猴, and 食蟹猴.

The generic labels `monkey`, `macaque`, `猴`, and `猕猴` are intentionally
rejected as ambiguous. A primate assessment must identify the exact species.

## Live source check

Run one profile check at a time:

```bash
python -m geneimpact source-check --species rat
python -m geneimpact source-check --species zebrafish
python -m geneimpact source-check --species fruit_fly
python -m geneimpact source-check --species rhesus_macaque
python -m geneimpact source-check --species cynomolgus_macaque
```

The command checks organism, taxonomy ID, assembly name, exact accession,
current/suppressed status, and reference-genome category. A source mismatch is
reported as an error rather than silently accepting coordinate drift.

## Predictor maturity

| Capability | Mouse | Rat | Zebrafish | Fruit fly | Rhesus macaque | Cynomolgus macaque |
|---|---:|---:|---:|---:|---:|---:|
| Context and assembly validation | Ready | Ready | Ready | Ready | Ready | Ready |
| Generic external-output applicability audit | Ready | Ready | Ready | Ready | Ready | Ready |
| MGI/IMPC phenotype evidence baseline | Research baseline | Not applicable | Not applicable | Not applicable | Not applicable | Not applicable |
| BE-Hive mES efficiency and bystander import | Narrow mES scope | Out of scope | Out of scope | Out of scope | Out of scope | Out of scope |
| inDelphi repair-outcome import | Narrow mESC scope | Out of scope | Out of scope | Out of scope | Out of scope | Out of scope |
| CRISPRscan guide-activity scoring | Out of scope | Out of scope | Narrow embryo scope | Out of scope | Out of scope | Out of scope |
| Housden guide-activity ranking | Out of scope | Out of scope | Out of scope | Narrow S2R+ cell scope | Out of scope | Out of scope |
| Port 2026 Cas12a array LOH evidence | Out of scope | Out of scope | Out of scope | Bounded in-vivo array audit | Out of scope | Out of scope |
| Independent species-specific transfer evidence | Mouse embryo retrospective | 14-guide rat embryo retrospective | Zebrafish RNP retrospective | Qualitative published comparisons | Hazard evidence only | Hazard evidence only |

An output from a human-cell or mouse-cell model is not promoted to another
species. It may remain in an audit report as `out_of_scope`, but it must not
contribute applicable evidence or a species-level performance claim.

The same matrix is available as machine-readable JSON:

```bash
python -m geneimpact capabilities --species mouse
python -m geneimpact capabilities --species zebrafish
python -m geneimpact capabilities --species rat
python -m geneimpact capabilities --species rhesus_macaque
python -m geneimpact readiness --species fruit_fly
python -m geneimpact readiness --all
```

The `readiness` command additionally qualifies public datasets and papers as
`usable_adapter`, `usable_bounded_benchmark`, `transfer_evidence_only`,
`hazard_evidence_only`, or `insufficient_public_data`. A bounded benchmark,
hazard observation, or method paper is never promoted to executable prediction
capability.

For rat, `validation_candidate` means that a pinned external-transfer evaluator
exists. It does not mean that GeneImpact ships or endorses a rat guide-activity
predictor. See the [rat transfer benchmark](rat-transfer-benchmark.md).

For fruit fly, `usable_bounded_benchmark` identifies checksum-pinned Port 2026
Cas12a evidence at the three- or four-guide array level. It does not mean that
the component guides, repair outcomes, phenotypes, or safety are predictable.
See the [fruit-fly Cas12a evidence guide](fruit-fly-cas12a-evidence.md).

Status meanings:

- `available_declared_domain`: a version-locked adapter exists, but only for
  the stated biological domain;
- `validation_candidate`: the published model is relevant enough to evaluate,
  but the adapter or independent validation is not complete;
- `reference_search_candidate`: the method can operate against a supplied
  species reference, but the GeneImpact adapter and empirical calibration are
  pending;
- `out_of_domain_only`: a service may technically produce a score, but the
  training evidence does not support treating it as validated for that species.

## Strain and population rule

The registered strain or isolate identifies the origin of the reference
assembly; it is not a claim that all experiments use that background. A
different declared strain, stock, colony, or geographic population is accepted
only with a visible warning that strain-specific validation is required.

For macaques, population structure and individual variation are especially
important. The reference isolate must never be interpreted as a universal
sequence for all rhesus or cynomolgus macaques.
