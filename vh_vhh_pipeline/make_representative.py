"""Representative-only release: one heavy-chain V domain per PDB entry.

python vh_vhh_pipeline/make_representative.py
-> vh_vhh_representative/ (+ vh_vhh_representative.zip)
"""
import argparse, os, shutil, zipfile

import gemmi
import numpy as np
import pandas as pd
import requests

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "vh_vhh_data")
OUT = os.path.join(ROOT, "vh_vhh_representative")
CAND = os.path.join(ROOT, "pdb_ab_nb_after_2024-09-09_novel_vs_train_val_fulllength.csv")


def rcsb_resolution(ids):
    q = "query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id rcsb_entry_info{resolution_combined}}}"
    res = {}
    for i in range(0, len(ids), 200):
        j = requests.post("https://data.rcsb.org/graphql", json={"query": q, "variables": {"ids": ids[i:i + 200]}},
                          timeout=120).json()
        for e in j["data"]["entries"]:
            rc = e["rcsb_entry_info"]["resolution_combined"]
            res[e["rcsb_id"]] = rc[0] if rc else None
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-resolution", type=float, default=3.5,
                    help="keep entries with resolution_A <= this; entries without a resolution (NMR) are dropped")
    ap.add_argument("--min-backbone-fraction", type=float, default=0.9,
                    help="drop records whose fraction of V residues with complete N/CA/C/O backbone is below this")
    ap.add_argument("--max-partial-residues", type=int, default=5,
                    help="drop records with more residues missing individual backbone atoms than this")
    args = ap.parse_args()
    md = pd.read_csv(os.path.join(SRC, "metadata.tsv"), sep="\t", keep_default_na=False)
    rep = md[md.is_representative.astype(str) == "True"].copy()
    assert rep.pdb_id.is_unique
    mp = pd.read_csv(os.path.join(SRC, "residue_mapping.tsv.gz"), sep="\t", keep_default_na=False)
    mp = mp[mp.sample_id.isin(rep.sample_id)]
    miss = pd.read_csv(os.path.join(SRC, "missing_atoms.tsv.gz"), sep="\t", keep_default_na=False)
    miss = miss[miss.sample_id.isin(rep.sample_id)]
    cand = pd.read_csv(CAND)
    cand["pdb_id"] = cand.pdb_id.str.upper()
    copies = md.groupby("pdb_id").size()
    full = {}
    with open(os.path.join(SRC, "full_chain.fasta")) as f:
        lines = f.read().split("\n")
    for h, q in zip(lines[0::2], lines[1::2]):
        if h.startswith(">"):
            full[h[1:].split("|")[0]] = q

    # resolution: entry-level value from the mmCIF, checked against RCSB resolution_combined
    rc = rcsb_resolution(sorted(rep.pdb_id))
    mine = pd.to_numeric(rep.resolution, errors="coerce")
    for p, v in zip(rep.pdb_id, mine):
        r = rc.get(p)
        assert (pd.isna(v) and r is None) or (r is not None and abs(v - r) < 1e-6), (p, v, r)

    if os.path.exists(OUT):
        shutil.rmtree(OUT)
    for d in ("structures_cif", "structures_pdb", "backbone_npz"):
        os.makedirs(os.path.join(OUT, d))

    rows = []
    for s in rep.itertuples():
        sid = s.sample_id
        g = mp[mp.sample_id == sid].sort_values("seq_idx1")
        assert list(g.seq_idx1) == list(range(1, s.L + 1)) and "".join(g.aa) == s.variable_sequence
        fs = full[f"{s.pdb_id}_{s.entity_id}"]
        assert fs[s.domain_start_full_idx1 - 1:s.domain_end_full_idx1] == s.variable_sequence
        reg = {k: "".join(gg.aa) for k, gg in g.groupby("region", sort=False)}
        nocoord = g[g.coordinate_record_present.astype(str) == "False"]
        partial = g[(g.coordinate_record_present.astype(str) == "True") & (g.backbone_mask.astype(str) == "False")]
        c = cand[cand.pdb_id == s.pdb_id].iloc[0]
        src_dir = os.path.join(SRC, "samples", sid)
        shutil.copy(os.path.join(src_dir, "backbone_1based.cif"), os.path.join(OUT, "structures_cif", f"{sid}.cif"))
        shutil.copy(os.path.join(src_dir, "backbone.npz"), os.path.join(OUT, "backbone_npz", f"{sid}.npz"))
        st = gemmi.read_structure(os.path.join(src_dir, "backbone_1based.cif"))
        st.setup_entities()
        st.write_pdb(os.path.join(OUT, "structures_pdb", f"{sid}.pdb"))
        # PDB format keeps 3 decimals; check residue numbering and coordinates (to 0.0005 A)
        st2 = gemmi.read_structure(os.path.join(OUT, "structures_pdb", f"{sid}.pdb"))
        a1 = {(r.seqid.num, a.name): a.pos for r in st[0][0] for a in r}
        a2 = {(r.seqid.num, a.name): a.pos for r in st2[0][0] for a in r}
        assert a1.keys() == a2.keys() and all(a1[k].dist(a2[k]) < 1e-3 for k in a1), sid
        res = mine[s.Index]
        rows.append(dict(
            sample_id=sid, pdb_id=s.pdb_id, entity_id=s.entity_id, label_asym_id=s.label_asym_id,
            auth_asym_id=s.auth_asym_id, model_id=s.model_id, domain_type=s.domain_type,
            domain_type_evidence=s.domain_type_evidence, entity_description=s.entity_description, title=s.title,
            method=s.method, resolution_A="" if pd.isna(res) else res, deposition_date=s.deposition_date,
            initial_release_date=s.initial_release_date, n_copies_in_entry=int(copies[s.pdb_id]),
            heavy_v_sequence=s.variable_sequence, heavy_v_length=s.L,
            cdr_h1_chothia=reg.get("CDR1", ""), cdr_h2_chothia=reg.get("CDR2", ""), cdr_h3_chothia=reg.get("CDR3", ""),
            fr1=reg.get("FR1", ""), fr2=reg.get("FR2", ""), fr3=reg.get("FR3", ""), fr4=reg.get("FR4", ""),
            chothia_positions=";".join(g.chothia_pos),
            full_chain_sequence=full[f"{s.pdb_id}_{s.entity_id}"], v_start_in_full_chain=s.domain_start_full_idx1,
            v_end_in_full_chain=s.domain_end_full_idx1,
            n_residues_with_coords=s.n_residues_with_coords, n_missing_residues=len(nocoord),
            missing_residues_seq_idx1=";".join(map(str, nocoord.seq_idx1)),
            missing_residues_chothia=";".join(nocoord.chothia_pos),
            n_residues_partial_backbone=len(partial), backbone_complete_fraction=round(s.n_full_backbone / s.L, 4),
            min_norm_lv_dist_to_train_val=c.min_norm_lv_dist_to_train_val, closest_train_val_pdb_id=c.closest_ref_pdb_id,
            structure_cif=f"structures_cif/{sid}.cif", structure_pdb=f"structures_pdb/{sid}.pdb",
            backbone_npz=f"backbone_npz/{sid}.npz", source_mmcif_url=s.mmcif_url, source_mmcif_sha256=s.mmcif_sha256))
    master = pd.DataFrame(rows).sort_values("pdb_id")
    rep_paths = master[["structure_cif", "structure_pdb", "backbone_npz"]].copy()
    n0 = len(master)
    r = pd.to_numeric(master.resolution_A, errors="coerce")
    dropped = master[~(r <= args.max_resolution)]
    master = master[r <= args.max_resolution]
    for col, d in (("structure_cif", "structures_cif"), ("structure_pdb", "structures_pdb"), ("backbone_npz", "backbone_npz")):
        for p in dropped[col]:
            os.remove(os.path.join(OUT, p))
    # quality: unknown residues / largely unmodeled / badly built backbones
    removed = [dict(sample_id=x, reason="resolution > cutoff or no resolution (NMR)") for x in dropped.sample_id]
    bad = {}
    for r_ in master.itertuples():
        why = []
        if "X" in r_.heavy_v_sequence:
            why.append(f"{r_.heavy_v_sequence.count('X')} unknown (X/UNK) residues in heavy V sequence")
        if r_.backbone_complete_fraction < args.min_backbone_fraction:
            why.append(f"only {r_.backbone_complete_fraction:.0%} of residues have a complete backbone")
        if r_.n_residues_partial_backbone > args.max_partial_residues:
            why.append(f"{r_.n_residues_partial_backbone} residues missing individual backbone atoms")
        if why:
            bad[r_.sample_id] = "; ".join(why)
    removed += [dict(sample_id=k, reason=v) for k, v in bad.items()]
    master = master[~master.sample_id.isin(bad)]
    # one record per identical heavy V sequence: best (lowest) resolution, then most complete backbone
    master = master.assign(_res=pd.to_numeric(master.resolution_A))
    master = master.sort_values(["heavy_v_sequence", "_res", "backbone_complete_fraction", "n_missing_residues", "pdb_id"],
                                ascending=[True, True, False, True, True])
    grp = master.groupby("heavy_v_sequence")
    master["n_entries_same_sequence"] = grp.pdb_id.transform("size")
    master["other_pdb_ids_same_sequence"] = grp.pdb_id.transform(lambda x: ";".join(x)).combine(
        master.pdb_id, lambda allp, p: ";".join(q for q in allp.split(";") if q != p))
    dup = master[master.duplicated("heavy_v_sequence", keep="first")]
    removed += [dict(sample_id=x, reason=f"duplicate sequence; kept {k} (better or equal resolution)")
                for x, k in zip(dup.sample_id, dup.heavy_v_sequence.map(master.drop_duplicates("heavy_v_sequence").set_index("heavy_v_sequence").sample_id))]
    master = master.drop_duplicates("heavy_v_sequence", keep="first").drop(columns="_res").sort_values("pdb_id")
    keep = set(master.sample_id)
    for col in ("structure_cif", "structure_pdb", "backbone_npz"):
        for p in rep_paths[col]:
            if os.path.basename(p).rsplit(".", 1)[0] not in keep and os.path.exists(os.path.join(OUT, p)):
                os.remove(os.path.join(OUT, p))
    pd.DataFrame(removed).to_csv(os.path.join(OUT, "removed_records.csv"), index=False)
    print(f"quality removals: {len(bad)}; duplicate-sequence removals: {len(dup)}")
    mp = mp[mp.sample_id.isin(master.sample_id)]
    miss = miss[miss.sample_id.isin(master.sample_id)]
    print(f"resolution filter <= {args.max_resolution} A: dropped {len(dropped)} of {n0} "
          f"({(r > args.max_resolution).sum()} worse than cutoff, {r.isna().sum()} without resolution)")
    master.to_csv(os.path.join(OUT, "master.csv"), index=False)
    mp.to_csv(os.path.join(OUT, "residue_mapping.csv.gz"), index=False)
    miss.to_csv(os.path.join(OUT, "missing_atoms.csv.gz"), index=False)
    with open(os.path.join(OUT, "heavy_v.fasta"), "w") as f:
        for r in master.itertuples():
            f.write(f">{r.sample_id}|{r.domain_type}|res={r.resolution_A}\n{r.heavy_v_sequence}\n")
    shutil.copy(os.path.join(ROOT, "vh_vhh_pipeline", "REPRESENTATIVE_README.md"), os.path.join(OUT, "README.md"))
    zp = os.path.join(ROOT, "vh_vhh_representative.zip")
    with zipfile.ZipFile(zp, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for d, _, fs in os.walk(OUT):
            for fn in sorted(fs):
                full = os.path.join(d, fn)
                z.write(full, os.path.relpath(full, ROOT))
    print(len(master), "representatives;", master.method.value_counts().to_dict())
    print(master.domain_type.value_counts().to_dict())
    print("zip MB", round(os.path.getsize(zp) / 1e6, 1))


if __name__ == "__main__":
    main()
