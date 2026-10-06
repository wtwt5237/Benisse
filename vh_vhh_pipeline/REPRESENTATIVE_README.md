# VH/VHH representative set (one record per unique heavy-chain V sequence)

Antibody/nanobody PDB entries released after 2024-09-09 (after AF3, Chai-1, ESM3, IgFold and RFAA),
heavy chain >10% different (normalized Levenshtein) from the train/val heavy chains, processed with
`.claude/skills/vh-vhh-data-collection` and strictly cleaned (see `vh_vhh_data/qc_report.md` and
`vh_vhh_data/excluded.tsv` in the repository). For each entry, the chain copy with the most complete
backbone is kept. Then:

1. only entries with overall resolution <= 3.5 A (NMR entries, which have no resolution, are excluded);
2. records with unknown (X/UNK) residues in the heavy V sequence, with <90% of residues having a complete
   N/CA/C/O backbone, or with >5 residues missing individual backbone atoms are removed;
3. **one record per unique heavy V sequence**: among entries with an identical sequence, the one with the
   best (lowest) resolution is kept (ties: more complete backbone). `n_entries_same_sequence` and
   `other_pdb_ids_same_sequence` list the entries that were collapsed.

Every removed record and the reason is in `removed_records.csv`.

## Files

| file | content |
|---|---|
| `master.csv` | one row per entry: identifiers, method, **resolution_A**, dates, heavy V sequence, Chothia CDR/FR sequences, full submitted chain sequence, missing-residue summary, novelty vs train/val, file paths |
| `structures_cif/<sample_id>.cif` | backbone (N, CA, C, O) of the heavy V domain, residues numbered 1..L along the V sequence (exact original coordinates) |
| `structures_pdb/<sample_id>.pdb` | same in PDB format (coordinates rounded to 0.001 A by the format) |
| `backbone_npz/<sample_id>.npz` | `coords [L,4,3]` (NaN = missing), `atom_mask [L,4]`, `backbone_mask [L]`, `sequence`, `chothia_positions`, `regions`, original numbering |
| `residue_mapping.csv.gz` | one row per V-domain residue: seq_idx1, Chothia position, region, original label/author numbering, coordinate presence |
| `missing_atoms.csv.gz` | every missing / zero-occupancy backbone atom with reason |
| `heavy_v.fasta` | heavy V sequences |
| `removed_records.csv` | records dropped from this release and why |

## Key conventions

- Residue numbers in the structure files = position in `heavy_v_sequence` (1-based). Residues without
  coordinates are absent from the structure files but keep their number, so numbering never shifts
  (e.g. residues 27-31 missing -> file jumps from 26 to 32). They are listed in
  `missing_residues_seq_idx1` / `missing_residues_chothia`.
- `resolution_A`: overall entry resolution in Angstrom (X-ray: refinement high-resolution limit;
  cryo-EM: reconstruction resolution), identical to RCSB `resolution_combined`. Empty for NMR.
- `domain_type`: VHH only with explicit nanobody/VHH annotation; VH when the entry has a light chain;
  unknown otherwise.
- Numbering: AbNumber, Chothia scheme and CDR definition (ANARCI human+mouse HMMs).
