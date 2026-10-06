# VH/VHH data collection QC report

Spec: `.claude/skills/vh-vhh-data-collection/SKILL.md`. Config and versions: `run_config.json`.

## Counts

| stage | n |
|---|---|
| candidate PDB entries | 1320 |
| downloaded (mmCIF + FASTA) | 1320 |
| heavy-domain samples before strict cleaning (all chain copies) | 2462 |
| removed (rows in excluded.tsv, incl. entry-level failures) | 46 |
| **clean samples kept** | 2420 |
| clean representative samples (1 per entry x heavy sequence) | 1295 |
| entries with >=1 clean sample | 1295 |

## Why records were removed (`excluded.tsv`)

| reason (a sample can have several) | n |
|---|---|
| non-standard Chothia insertions | 22 |
| manual review | 13 |
| pipeline_review | 7 |
| processing | 4 |

## Domain type (clean samples)

| domain_type | samples | entries |
|---|---|---|
| VH | 1831 | 979 |
| VHH | 580 | 312 |
| unknown | 9 | 4 |

VHH is assigned only from an explicit entity annotation (nanobody/VHH/sybody/single-domain antibody). AbNumber H type alone is not used as VHH evidence.

## Missing-coordinate statistics (clean samples)

| reason | atom records |
|---|---|
| missing_residue | 11684 |
| zero_occupancy | 94 |
| missing_atom | 14 |

Residues: 297415; with coordinate record: 294494 (99.02%); full backbone: 294403 (98.99%). Samples with >=1 missing residue: 904.

Chothia alignment: 157 common columns over 2420 clean samples.

## Acceptance checks

| # | check | result |
|---|---|---|
| 1 | L == len(variable_sequence) == mapping rows == coords.shape[0] == atom_mask.shape[0] | PASS |
| 2 | seq_idx1 is 1..L; variable sequence restorable from full submitted sequence via recorded boundaries | PASS |
| 3 | coordinate residues map uniquely; chain/model identity retained | PASS |
| 4 | formal samples have no unexplained sequence/structure mismatch | PASS |
| 5 | 1-based export: coordinates unchanged, residue IDs = seq_idx1, atom count = valid atoms | PASS |
| 6 | atom_mask=True coords finite, masked coords NaN, backbone_mask == all(atom_mask) | PASS |
| 7 | Chothia alignment minus gaps restores sequence; structure-missing residues not turned into gaps | PASS |
| 8 | every candidate entry has a final status | PASS |
