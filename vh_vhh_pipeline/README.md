# VH/VHH data collection pipeline

Implements `.claude/skills/vh-vhh-data-collection/SKILL.md` for the 1,320 candidate entries in
`pdb_ab_nb_after_2024-09-09_novel_vs_train_val_fulllength.csv`.

```
pip install abnumber gemmi pandas numpy requests     # plus HMMER (hmmscan) for ANARCI
python vh_vhh_pipeline/pipeline.py --limit 5          # trial run
python vh_vhh_pipeline/pipeline.py                    # full batch -> vh_vhh_data/
python -m pytest vh_vhh_pipeline/test_mapping.py      # synthetic mapping regression tests
```

| file | role |
|---|---|
| `config.json` | sources, numbering scheme, altloc/model/representative rules, explicitly unrestricted filters |
| `pipeline.py` | download -> full sequence -> AbNumber Chothia -> label_seq_id mapping -> 1-based export -> backbone + masks -> Chothia alignment |
| `mapping.py` | per-chain sequence/structure mapping (no I/O) |
| `qc.py` | the 8 acceptance checks and `qc_report.md` |
| `test_mapping.py` | synthetic regression cases (N-term gap, internal gap, missing O, OXT, insertion codes, compressed author numbering, altloc, zero occupancy, microheterogeneity, mismatch, repeats) |

Outputs in `vh_vhh_data/` (see `qc_report.md` there): `metadata.tsv` (one row per heavy-domain sample),
`entries.tsv` (one row per candidate entry), `residue_mapping.tsv.gz`, `missing_atoms.tsv.gz`,
`full_chain.fasta`, `variable.fasta`, `samples/<sample_id>/backbone.npz` + `backbone_1based.cif`,
`chothia_alignment.{tsv,npz}`, `chothia_columns.tsv`, `failed_or_review.tsv`, `download_log.tsv`,
`run_config.json`.

`sample_id` = `{PDB}_e{entity_id}_{label_asym_id}_m{model}_d{domain}`. Use `status == formal`
(and `is_representative` for one copy per entry/heavy sequence).
