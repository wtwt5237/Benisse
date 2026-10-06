"""Acceptance checks (spec section 10) and qc_report.md."""
import collections, os, re
import numpy as np
import pandas as pd

BB = ("N", "CA", "C", "O")


def run(out, samples, formal, entry_rows, excluded, ids, aln_rows, col_names, status, n_all):
    checks = collections.OrderedDict()
    fails = collections.defaultdict(list)
    mp = pd.read_csv(os.path.join(out, "residue_mapping.tsv.gz"), sep="\t", keep_default_na=False)
    rows_per = mp.groupby("sample_id").size().to_dict()
    for s in samples:
        sid = s["sample_id"]
        z = np.load(os.path.join(out, "samples", sid, "backbone.npz"))
        L = s["L"]
        # 1 lengths
        if not (L == len(s["variable_sequence"]) == rows_per.get(sid) == z["coords"].shape[0] == z["atom_mask"].shape[0]):
            fails[1].append(sid)
        # 2 seq_idx1 1..L, sequence restorable from full sequence via recorded boundary
        sub = mp[mp.sample_id == sid]
        if list(sub.seq_idx1) != list(range(1, L + 1)) or \
                s["full_chain_sequence"][s["domain_start_full_idx1"] - 1:s["domain_end_full_idx1"]] != s["variable_sequence"] or \
                "".join(sub.aa) != s["variable_sequence"]:
            fails[2].append(sid)
        # 3 unique mapping of coordinate residues; identities kept
        present = sub[sub.coordinate_record_present.astype(str) == "True"]
        if present.label_seq_id.duplicated().any() or (present.auth_asym_id == "").any() or \
                (sub.model_id.astype(str) == "").any():
            fails[3].append(sid)
        # 4 residue identity (formal samples must have no unexplained mismatch)
        if s["status"] == "formal" and any(not str(m).startswith("ok") for m in sub.mapping_status):
            fails[4].append(sid)
        # 5 export keeps coordinates and seq_idx1 numbering
        if s["export_check"] not in ("ok", "not_written_no_valid_backbone_atoms"):
            fails[5].append(sid)
        # 6 masks
        c, m = z["coords"], z["atom_mask"]
        if not np.isfinite(c[m]).all() or np.isfinite(c[~m]).any() or not (z["backbone_mask"] == m.all(1)).all():
            fails[6].append(sid)
    # 7 alignment restores sequences, structure-missing kept as residues
    seq_of = {s["sample_id"]: s for s in formal}
    for r in aln_rows:
        s = seq_of[r["sample_id"]]
        if r["aligned_sequence"].replace("-", "") != s["variable_sequence"]:
            fails[7].append(r["sample_id"])
        nres = sum(1 for ch in r["structure_status"] if ch != "-")
        if nres != s["L"]:
            fails[7].append(r["sample_id"])
    # 8 every candidate has final status
    have = {e["pdb_id"] for e in entry_rows}
    missing = [p for p in ids if p not in have]
    if missing:
        fails[8] += missing
    names = {1: "L == len(variable_sequence) == mapping rows == coords.shape[0] == atom_mask.shape[0]",
             2: "seq_idx1 is 1..L; variable sequence restorable from full submitted sequence via recorded boundaries",
             3: "coordinate residues map uniquely; chain/model identity retained",
             4: "formal samples have no unexplained sequence/structure mismatch",
             5: "1-based export: coordinates unchanged, residue IDs = seq_idx1, atom count = valid atoms",
             6: "atom_mask=True coords finite, masked coords NaN, backbone_mask == all(atom_mask)",
             7: "Chothia alignment minus gaps restores sequence; structure-missing residues not turned into gaps",
             8: "every candidate entry has a final status"}
    for k in range(1, 9):
        checks[k] = (names[k], "PASS" if not fails[k] else f"FAIL ({len(set(fails[k]))}: {sorted(set(fails[k]))[:5]})")

    md = pd.DataFrame([{k: v for k, v in s.items() if k not in ("rows",)} for s in samples])
    miss = pd.read_csv(os.path.join(out, "missing_atoms.tsv.gz"), sep="\t") if os.path.getsize(os.path.join(out, "missing_atoms.tsv.gz")) > 30 else pd.DataFrame(columns=["reason", "sample_id"])
    ent = pd.DataFrame(entry_rows)
    L = []
    L.append("# VH/VHH data collection QC report\n")
    L.append("Spec: `.claude/skills/vh-vhh-data-collection/SKILL.md`. Config and versions: `run_config.json`.\n")
    L.append("## Counts\n")
    n_dl = sum(1 for p in ids if status.get(p) in ("downloaded",))
    L.append("| stage | n |\n|---|---|")
    L.append(f"| candidate PDB entries | {len(ids)} |")
    L.append(f"| downloaded (mmCIF + FASTA) | {n_dl} |")
    L.append(f"| heavy-domain samples before strict cleaning (all chain copies) | {n_all} |")
    L.append(f"| removed (rows in excluded.tsv, incl. entry-level failures) | {len(excluded)} |")
    L.append(f"| **clean samples kept** | {len(md)} |")
    if len(md):
        L.append(f"| clean representative samples (1 per entry x heavy sequence) | {int(md.is_representative.sum())} |")
        L.append(f"| entries with >=1 clean sample | {md.pdb_id.nunique()} |")
    L.append("\n## Why records were removed (`excluded.tsv`)\n")
    if len(excluded):
        reasons = collections.Counter()
        for r in excluded:
            for part in str(r["reason"]).split(" | "):
                key = part.split(": ")[0] if ": " in part else re.sub(r"=\S+", "=x", part)
                reasons[re.sub(r" [A-Z0-9,]+$", "", key)[:80]] += 1
        L.append("| reason (a sample can have several) | n |\n|---|---|")
        for rs, n in reasons.most_common():
            L.append(f"| {rs} | {n} |")
    else:
        L.append("none")
    if len(md):
        L.append("\n## Domain type (clean samples)\n")
        f = md
        L.append("| domain_type | samples | entries |\n|---|---|---|")
        for t, g in f.groupby("domain_type"):
            L.append(f"| {t} | {len(g)} | {g.pdb_id.nunique()} |")
        L.append("\nVHH is assigned only from an explicit entity annotation (nanobody/VHH/sybody/single-domain antibody). AbNumber H type alone is not used as VHH evidence.\n")
        L.append("## Missing-coordinate statistics (clean samples)\n")
        fm = miss[miss.sample_id.isin(set(f.sample_id))] if len(miss) else miss
        L.append("| reason | atom records |\n|---|---|")
        for rs, n in fm.reason.value_counts().items():
            L.append(f"| {rs} | {n} |")
        tot = f.L.sum()
        L.append(f"\nResidues: {tot}; with coordinate record: {f.n_residues_with_coords.sum()} "
                 f"({100*f.n_residues_with_coords.sum()/max(tot,1):.2f}%); full backbone: {f.n_full_backbone.sum()} "
                 f"({100*f.n_full_backbone.sum()/max(tot,1):.2f}%). Samples with >=1 missing residue: "
                 f"{(f.n_residues_with_coords < f.L).sum()}.\n")
        L.append(f"Chothia alignment: {len(col_names)} common columns over {len(aln_rows)} clean samples.\n")
    L.append("## Acceptance checks\n")
    L.append("| # | check | result |\n|---|---|---|")
    for k, (n, r) in checks.items():
        L.append(f"| {k} | {n} | {r} |")
    open(os.path.join(out, "qc_report.md"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))
