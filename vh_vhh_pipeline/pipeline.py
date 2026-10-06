"""VH/VHH data collection following .claude/skills/vh-vhh-data-collection.

Usage:
  python vh_vhh_pipeline/pipeline.py --limit 5            # trial run on first 5 candidates
  python vh_vhh_pipeline/pipeline.py --ids 9ABC,9XYZ      # specific entries
  python vh_vhh_pipeline/pipeline.py                      # full batch from config candidate table

Steps: download original mmCIF + entry FASTA -> full submitted sequences (FASTA, checked
against _entity_poly / _entity_poly_seq) -> AbNumber (Chothia) domain parse on the full
sequence -> explicit label_seq_id mapping -> 1-based export -> backbone [N, CA, C, O]
with masks -> dataset Chothia alignment -> QC.
"""
import argparse, datetime, gzip, hashlib, io, json, math, os, re, sys, time, traceback, warnings
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
import requests
import gemmi
from anarci import anarci as run_anarci
import abnumber
from abnumber import Chain

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mapping import map_chain, BACKBONE  # noqa: E402

warnings.filterwarnings("ignore", category=UserWarning)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CFG = json.load(open(os.path.join(ROOT, "vh_vhh_pipeline", "config.json")))
VHH_KW = re.compile(r"nanobod|\bvhh\b|sybod|single[- ]domain antibod|\bsdab\b|heavy[- ]chain[- ]only", re.I)


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# ---------------------------------------------------------------- download
def fetch(url, dest, log, tries=5):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return "cached"
    err = ""
    for k in range(tries):
        try:
            r = requests.get(url, timeout=120)
            if r.status_code == 200 and r.content:
                tmp = dest + ".part"
                open(tmp, "wb").write(r.content)
                os.replace(tmp, dest)
                log.append(dict(url=url, time=now(), attempt=k + 1, status="ok"))
                return "downloaded"
            err = f"HTTP {r.status_code}"
        except Exception as e:  # network errors are retried, then recorded
            err = repr(e)
        log.append(dict(url=url, time=now(), attempt=k + 1, status=err))
        time.sleep(2 ** k)
    raise RuntimeError(f"download failed after {tries} attempts: {url}: {err}")


def download_entry(pid, raw_dir, log):
    P = pid.upper()
    cif = os.path.join(raw_dir, f"{P}.cif.gz")
    fa = os.path.join(raw_dir, f"{P}.fasta")
    s1 = fetch(CFG["sources"]["mmcif"].format(PDB_ID=P), cif, log)
    s2 = fetch(CFG["sources"]["fasta_entry"].format(PDB_ID=P), fa, log)
    return cif, fa, s1, s2


STANDARD = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())


def ccd_parent(comp, raw_dir, log):
    """Parent residue and one-letter code of a modified residue from the wwPDB CCD."""
    d = os.path.join(raw_dir, "ccd")
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, f"{comp}.cif")
    fetch(f"https://files.rcsb.org/ligands/view/{comp}.cif", path, log)
    b = gemmi.cif.read(path).sole_block()
    return (val(b.find_value("_chem_comp.mon_nstd_parent_comp_id")) or "",
            val(b.find_value("_chem_comp.one_letter_code")) or "",
            gemmi.cif.as_string(b.find_value("_chem_comp.name") or ""))


# ---------------------------------------------------------------- parsing
def cat(block, name):
    c = block.get_mmcif_category(name)
    return c if c else {}


def val(v):
    return None if v in (None, False, "?", ".") else v


def parse_fasta(path, P):
    out = {}
    header, seq = None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if header:
                out[header] = "".join(seq)
            header, seq = line[1:], []
        elif line:
            seq.append(line.strip())
    if header:
        out[header] = "".join(seq)
    ent = {}
    for h, s in out.items():
        m = re.match(rf"{P}_(\d+)\|", h)
        if m:
            ent[m.group(1)] = (h, s)
    return ent


def entry_info(block):
    g = lambda c, k: (cat(block, c).get(k) or [None])
    res = None
    for c, k in (("_refine.", "ls_d_res_high"), ("_em_3d_reconstruction.", "resolution"), ("_reflns.", "d_resolution_high")):
        v = [val(x) for x in g(c, k) if val(x)]
        if v:
            res = v[0]
            break
    revs = [val(x) for x in g("_pdbx_audit_revision_history.", "revision_date") if val(x)]
    return dict(
        title=(val(g("_struct.", "title")[0]) or ""),
        method=";".join(sorted({x for x in g("_exptl.", "method") if val(x)})),
        resolution=res or "",
        deposition_date=val(g("_pdbx_database_status.", "recvd_initial_deposition_date")[0]) or "",
        initial_release_date=min(revs) if revs else "",
    )


def parse_cif(path, P, fasta):
    doc = gemmi.cif.read(path)
    block = doc.sole_block()
    info = entry_info(block)
    ep = cat(block, "_entity_poly.")
    eps = cat(block, "_entity_poly_seq.")
    ent = cat(block, "_entity.")
    desc = dict(zip(ent.get("id", []), ent.get("pdbx_description", [])))
    asym = cat(block, "_struct_asym.")
    asym_entity = dict(zip(asym.get("id", []), asym.get("entity_id", [])))
    entities = {}
    for i, eid in enumerate(ep.get("entity_id", [])):
        typ = val(ep["type"][i]) or ""
        can = re.sub(r"\s+", "", val(ep["pdbx_seq_one_letter_code_can"][i]) or "")
        entities[eid] = dict(entity_id=eid, type=typ, can=can, description=val(desc.get(eid)) or "",
                             asyms=[a for a, e in asym_entity.items() if e == eid], poly_seq={}, issues=[])
    for eid, num, mon in zip(eps.get("entity_id", []), eps.get("num", []), eps.get("mon_id", [])):
        if eid in entities:
            entities[eid]["poly_seq"].setdefault(int(num), []).append(mon)
    for eid, e in entities.items():
        nums = sorted(e["poly_seq"])
        if nums != list(range(1, len(nums) + 1)):
            e["issues"].append("entity_poly_seq_num_not_contiguous")
        if len(e["can"]) != len(nums):
            e["issues"].append(f"can_length_{len(e['can'])}_vs_entity_poly_seq_{len(nums)}")
        fh = fasta.get(eid)
        e["fasta_header"], e["fasta_seq"] = (fh if fh else ("", ""))
        if not fh:
            e["issues"].append("no_fasta_record")
        elif fh[1] != e["can"]:
            e["issues"].append("fasta_vs_mmcif_sequence_conflict")
    # poly seq scheme
    pss = cat(block, "_pdbx_poly_seq_scheme.")
    scheme = {}
    for i in range(len(pss.get("asym_id", []))):
        scheme.setdefault(pss["asym_id"][i], {})[int(pss["seq_id"][i])] = dict(
            auth_seq_num=val(pss["pdb_seq_num"][i]), ins_code=val(pss["pdb_ins_code"][i]),
            strand_id=val(pss["pdb_strand_id"][i]), mon_id=pss["mon_id"][i])
    # unobserved / zero occupancy annotations
    ur = cat(block, "_pdbx_unobs_or_zero_occ_residues.")
    unobs_res = {}
    for i in range(len(ur.get("label_asym_id", []))):
        if val(ur["label_seq_id"][i]):
            unobs_res.setdefault((ur["label_asym_id"][i], int(ur["PDB_model_num"][i])), set()).add(int(ur["label_seq_id"][i]))
    ua = cat(block, "_pdbx_unobs_or_zero_occ_atoms.")
    unobs_atoms = {}
    for i in range(len(ua.get("label_asym_id", []))):
        if val(ua["label_seq_id"][i]):
            unobs_atoms.setdefault((ua["label_asym_id"][i], int(ua["PDB_model_num"][i])), set()).add(
                (int(ua["label_seq_id"][i]), ua["label_atom_id"][i]))
    return block, info, entities, scheme, unobs_res, unobs_atoms


def atom_rows_for(block, asyms):
    t = block.find("_atom_site.", ["group_PDB", "label_atom_id", "label_alt_id", "label_comp_id", "label_asym_id",
                                    "label_seq_id", "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z",
                                    "occupancy", "B_iso_or_equiv", "auth_seq_id", "auth_asym_id",
                                    "pdbx_PDB_model_num", "type_symbol"])
    rows = {}
    for r in t:
        a = r[4]
        if a not in asyms:
            continue
        ls = r[5]
        alt = r[2]
        rows.setdefault(a, []).append(dict(
            group=r[0], atom=gemmi.cif.as_string(r[1]), alt=None if alt in (".", "?") else alt, comp=r[3],
            label_seq_id=None if ls in (".", "?") else int(ls), ins_code=None if r[6] in (".", "?") else r[6],
            x=float(r[7]), y=float(r[8]), z=float(r[9]),
            occ=None if r[10] in (".", "?") else float(r[10]), b=None if r[11] in (".", "?") else float(r[11]),
            auth_seq_id=r[12], auth_asym_id=r[13], model=int(r[14]), element=r[15]))
    return rows


# ---------------------------------------------------------------- numbering
def parse_domains(seqs):
    """seqs: dict key -> full sequence. Returns key -> list of domain dicts (all chain types)."""
    keys = list(seqs)
    out = {k: [] for k in keys}
    errs = {}
    for i in range(0, len(keys), 300):
        chunk = keys[i:i + 300]
        numbered, ali, _ = run_anarci([(k, seqs[k]) for k in chunk], scheme="chothia",
                                      allowed_species=CFG["numbering"]["allowed_species"], ncpu=4)
        chains, cerr = Chain.batch({k: seqs[k] for k in chunk}, scheme="chothia", cdr_definition="chothia",
                                   allowed_species=CFG["numbering"]["allowed_species"], multiple_domains=True)
        errs.update(cerr)
        for k, num, al in zip(chunk, numbered, ali):
            if not num:
                continue
            abn = chains.get(k, [])
            if len(abn) != len(num):
                errs[k] = f"AbNumber found {len(abn)} domains, ANARCI boundary parse found {len(num)}"
                continue
            for d, ((positions, start, end), a, ch) in enumerate(zip(num, al, abn), 1):
                dom_seq = seqs[k][start:end + 1]
                resid = "".join(aa for _, aa in positions if aa != "-")
                ok = (ch.seq == dom_seq == resid) and ch.chain_type == a["chain_type"]
                out[k].append(dict(domain_index=d, chain_type=ch.chain_type, start0=start, end0=end, chain=ch,
                                   boundary_check="ok" if ok else
                                   f"boundary_mismatch abnumber={ch.seq} slice={dom_seq} anarci={resid}"))
    return out, errs


# ---------------------------------------------------------------- export
def write_backbone_cif(path, sample, chain_name):
    st = gemmi.Structure()
    st.name = sample["sample_id"]
    try:
        model = gemmi.Model(1)
    except TypeError:
        model = gemmi.Model("1")
    ch = gemmi.Chain(chain_name)
    for i, r in enumerate(sample["rows"]):
        valid = [n for n in BACKBONE if r["atoms"][n]["valid"]]
        if not valid:
            continue
        res = gemmi.Residue()
        res.name = r["original_resname"]
        res.seqid = gemmi.SeqId(i + 1, " ")
        res.label_seq = i + 1
        res.subchain = "A"
        res.entity_type = gemmi.EntityType.Polymer
        res.het_flag = "A" if gemmi.find_tabulated_residue(res.name) and gemmi.find_tabulated_residue(res.name).is_standard() else "H"
        for n in valid:
            a = r["atoms"][n]
            at = gemmi.Atom()
            at.name = n
            at.element = gemmi.Element(r["elements"][n])
            at.pos = gemmi.Position(a["x"], a["y"], a["z"])
            at.occ = a["occ"]
            at.b_iso = r["b"][n]
            res.add_atom(at)
        ch.add_residue(res)
    model.add_chain(ch)
    st.add_model(model)
    ent = gemmi.Entity("1")
    ent.entity_type = gemmi.EntityType.Polymer
    ent.polymer_type = gemmi.PolymerType.PeptideL
    ent.subchains = ["A"]
    ent.full_sequence = [r["seq_mon_ids"] for r in sample["rows"]]
    st.entities.append(ent)
    doc = st.make_mmcif_document()
    doc.write_file(path)
    # verify: coordinates unchanged, numbering == seq_idx1
    st2 = gemmi.read_structure(path)
    seen = 0
    for res in st2[0][0]:
        r = sample["rows"][res.seqid.num - 1]
        if res.name != r["original_resname"]:
            return f"exported residue name mismatch at {res.seqid.num}"
        for at in res:
            a = r["atoms"][at.name]
            if max(abs(at.pos.x - a["x"]), abs(at.pos.y - a["y"]), abs(at.pos.z - a["z"])) > 1e-6:
                return f"coordinate changed at {res.seqid.num} {at.name}"
            seen += 1
    expect = sum(r["atoms"][n]["valid"] for r in sample["rows"] for n in BACKBONE)
    return "ok" if seen == expect else f"exported atom count {seen} != valid atoms {expect}"


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int)
    ap.add_argument("--ids")
    ap.add_argument("--out", default=os.path.join(ROOT, CFG["output_dir"]))
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--download-only", action="store_true")
    args = ap.parse_args()
    out = args.out
    raw = os.path.join(out, "raw")
    sdir = os.path.join(out, "samples")
    os.makedirs(raw, exist_ok=True)
    os.makedirs(sdir, exist_ok=True)

    cand = pd.read_csv(os.path.join(ROOT, CFG["candidate_table"]))
    if args.ids:
        want = {x.strip().lower() for x in args.ids.split(",")}
        cand = cand[cand.pdb_id.str.lower().isin(want)]
    if args.limit:
        cand = cand.head(args.limit)
    ids = [p.upper() for p in cand.pdb_id]
    prior_seq = dict(zip(cand.pdb_id.str.upper(), cand.heavy_chain_variable_seq))
    run = dict(started=now(), command=" ".join(sys.argv), n_candidates=len(ids),
               versions=dict(abnumber=getattr(abnumber, "__version__", "?"), gemmi=gemmi.__version__,
                             numpy=np.__version__, pandas=pd.__version__, python=sys.version.split()[0]),
               config=CFG)

    # 1. download (cached, retried, every failure recorded)
    dl_log, status, failed = [], {}, []
    def _dl(p):
        log = []
        try:
            r = download_entry(p, raw, log)
            return p, r, log, None
        except Exception as e:
            return p, None, log, str(e)
    with ThreadPoolExecutor(args.workers) as ex:
        for p, r, log, err in ex.map(_dl, ids):
            dl_log += [dict(pdb_id=p, **l) for l in log]
            if err:
                status[p] = "download_failed"
                failed.append(dict(stage="download", error_type="download_failed", reason=err, pdb_id=p, sample_id="",
                                   attempted_fix="5 retries with backoff"))
            else:
                status[p] = "downloaded"

    if args.download_only:
        pd.DataFrame(dl_log).to_csv(os.path.join(out, "download_log.tsv"), sep="\t", index=False)
        print({k: sum(1 for v in status.values() if v == k) for k in set(status.values())})
        return

    # 2. parse + full sequences
    parsed = {}
    for p in ids:
        if status[p] != "downloaded":
            continue
        cif, fa = os.path.join(raw, f"{p}.cif.gz"), os.path.join(raw, f"{p}.fasta")
        try:
            fasta = parse_fasta(fa, p)
            block, info, ents, scheme, ur, ua = parse_cif(cif, p, fasta)
            info.update(mmcif_sha256=sha256(cif), fasta_sha256=sha256(fa),
                        mmcif_url=CFG["sources"]["mmcif"].format(PDB_ID=p),
                        fasta_url=CFG["sources"]["fasta_entry"].format(PDB_ID=p),
                        retrieved=datetime.datetime.fromtimestamp(os.path.getmtime(cif), datetime.timezone.utc).isoformat(timespec="seconds"))
            del block  # keep only extracted metadata; the file is re-read per entry later (memory)
            parsed[p] = dict(cif=cif, info=info, ents=ents, scheme=scheme, ur=ur, ua=ua)
        except Exception as e:
            status[p] = "parse_failed"
            failed.append(dict(stage="parse", error_type="parse_failed", reason=repr(e), pdb_id=p, sample_id="",
                               attempted_fix=""))

    # 3. AbNumber Chothia domain parse on full submitted sequences (all protein entities)
    seqs = {f"{p}|{eid}": e["can"] for p, d in parsed.items() for eid, e in d["ents"].items()
            if e["type"].startswith("polypeptide") and len(e["can"]) >= 70}
    doms, derr = parse_domains(seqs)

    samples, entry_rows = [], []
    full_fa, var_fa = [], []
    for p in ids:
        if p not in parsed:
            entry_rows.append(dict(pdb_id=p, entry_status=status[p]))
            continue
        d = parsed[p]
        has_light = any(dm["chain_type"] in "KL" for k, v in doms.items() if k.startswith(p + "|") for dm in v)
        heavy = [(eid, dm) for eid in d["ents"] for dm in doms.get(f"{p}|{eid}", []) if dm["chain_type"] == "H"]
        for eid in d["ents"]:
            err = derr.get(f"{p}|{eid}")
            if err and err.startswith("Variable chain sequence not recognized"):
                d["ents"][eid]["antibody_domain"] = "none (not an antibody V domain)"
            elif err:
                failed.append(dict(stage="numbering", error_type="numbering_error", reason=err,
                                   pdb_id=p, sample_id=f"{p}_{eid}", attempted_fix=""))
        heavy_asyms = {a for eid, _ in heavy for a in d["ents"][eid]["asyms"]}
        entry_atoms = atom_rows_for(gemmi.cif.read(d["cif"]).sole_block(), heavy_asyms) if heavy_asyms else {}
        distinct = sorted({dm["chain"].seq for _, dm in heavy})
        entry_issue = []
        if not heavy:
            entry_issue.append("no_heavy_domain_found")
        if len(distinct) > 1:
            entry_issue.append(f"{len(distinct)}_distinct_heavy_domains")
        if prior_seq.get(p) not in distinct:
            entry_issue.append("chothia_domain_differs_from_screening_sequence")
        for eid, dm in heavy:
            e = d["ents"][eid]
            ch = dm["chain"]
            positions = list(ch)
            assert "".join(aa for _, aa in positions) == ch.seq
            first, last = positions[0][0], positions[-1][0]
            complete = first.number == 1 and last.number == 113 and not first.letter and not last.letter
            hallmark = {f"H{n}": ch[f"H{n}"] if f"H{n}" in [pp.format() for pp, _ in positions] else "-" for n in (37, 44, 45, 47)}
            if VHH_KW.search(e["description"]):
                dtype, evidence = "VHH", f"entity description: {e['description']}"
            elif has_light:
                dtype, evidence = "VH", "light-chain V domain present in entry (pairing not verified)"
            else:
                dtype, evidence = "unknown", "no VHH annotation on entity and no light-chain V domain in entry"
            for asym in e["asyms"]:
                rows_all = entry_atoms.get(asym, [])
                models = sorted({r["model"] for r in rows_all})
                model = models[0] if models else 1
                rows = [r for r in rows_all if r["model"] == model]
                mapped, issues = map_chain(e["poly_seq"], d["scheme"].get(asym, {}), rows,
                                           d["ur"].get((asym, model)), d["ua"].get((asym, model)))
                auth_asym = sorted({r["auth_asym_id"] for r in rows}) or [
                    (d["scheme"].get(asym, {}).get(1) or {}).get("strand_id") or ""]
                sid = f"{p}_e{eid}_{asym}_m{model}_d{dm['domain_index']}"
                elem = {}
                bfac = {}
                for r in rows:
                    if r["label_seq_id"] is not None and r["atom"] in BACKBONE:
                        elem.setdefault(r["label_seq_id"], {})[r["atom"]] = r["element"]
                        bfac.setdefault(r["label_seq_id"], {})[r["atom"]] = r["b"]
                srows = []
                for i, (pos, aa) in enumerate(positions):
                    num = dm["start0"] + i + 1  # 1-based position in full submitted sequence == entity_poly_seq.num
                    m = mapped[num]
                    can_aa = e["can"][num - 1]
                    st = m["mapping_status"]
                    norm = ""
                    nonstd = [c for c in m["seq_mon_ids"].split(",") if c not in STANDARD]
                    for comp in nonstd:
                        try:
                            par, one, nm = ccd_parent(comp, raw, dl_log)
                        except Exception as ex:
                            par, one, nm = "", "", f"CCD lookup failed: {ex}"
                        norm += f"{comp}({nm})->{can_aa} via _entity_poly.pdbx_seq_one_letter_code_can; CCD parent={par or '?'} one_letter={one or '?'};"
                        if one.strip().upper() != can_aa and st.startswith("ok"):
                            st = "modified_residue_parent_conflict"
                        elif st == "ok":
                            st = "ok_modified_residue"
                    if can_aa != aa:
                        st = "mismatch_sequence_vs_numbering"
                    srows.append(dict(m, sample_id=sid, full_seq_idx1=num, seq_idx1=i + 1, aa=aa,
                                      chothia_pos=pos.format(), region=pos.get_region(), entity_id=eid,
                                      label_asym_id=asym, auth_asym_id=",".join(x for x in auth_asym if x),
                                      model_id=model, mapping_status=st, residue_normalization=norm,
                                      elements={n: elem.get(num, {}).get(n, n[0]) for n in BACKBONE},
                                      b={n: bfac.get(num, {}).get(n) or 0.0 for n in BACKBONE}))
                L = len(srows)
                bad = sorted({r["mapping_status"] for r in srows if not r["mapping_status"].startswith("ok")})
                n_valid = sum(r["atoms"][n]["valid"] for r in srows for n in BACKBONE)
                sample = dict(sample_id=sid, pdb_id=p, entity_id=eid, label_asym_id=asym,
                              auth_asym_id=",".join(x for x in auth_asym if x), model_id=model,
                              n_models=len(models), domain_id=f"d{dm['domain_index']}",
                              domain_type=dtype, domain_type_evidence=evidence,
                              vhh_hallmark_chothia=";".join(f"{k}:{v}" for k, v in hallmark.items()),
                              entity_description=e["description"], full_chain_length=len(e["can"]),
                              domain_start_full_idx1=dm["start0"] + 1, domain_end_full_idx1=dm["end0"] + 1,
                              n_term_prefix_removed=dm["start0"], c_term_removed=len(e["can"]) - dm["end0"] - 1,
                              L=L, variable_sequence=ch.seq, full_chain_sequence=e["can"],
                              sequence_source="RCSB entry FASTA; verified vs _entity_poly.pdbx_seq_one_letter_code_can and _entity_poly_seq",
                              fasta_header=e["fasta_header"], sequence_issues=";".join(e["issues"]),
                              boundary_check=dm["boundary_check"], chothia_complete_H1_H113=complete,
                              n_residues_with_coords=sum(r["coordinate_record_present"] for r in srows),
                              n_full_backbone=sum(all(r["atoms"][n]["valid"] for n in BACKBONE) for r in srows),
                              n_valid_backbone_atoms=n_valid, mapping_problems=";".join(bad),
                              chain_issues=";".join(f"{a}:{b}" for a, b in issues[:5]),
                              entry_issues=";".join(entry_issue), rows=srows, **d["info"])
                reasons = []
                if e["issues"]:
                    reasons.append("sequence_source:" + ";".join(e["issues"]))
                if dm["boundary_check"] != "ok":
                    reasons.append("domain_boundary")
                if bad:
                    reasons.append("mapping:" + ";".join(bad))
                if issues:
                    reasons.append("chain_mapping_issues")
                if not complete:
                    reasons.append("v_domain_not_complete_H1_to_H113")
                if entry_issue:
                    reasons.append("entry:" + ";".join(entry_issue))
                if sample["n_residues_with_coords"] == 0:
                    reasons.append("no_coordinates")
                sample["status"] = "formal" if not reasons else "review"
                sample["status_reason"] = " | ".join(reasons)
                samples.append(sample)
        entry_rows.append(dict(pdb_id=p, entry_status="processed", n_heavy_domain_samples=sum(1 for s in samples if s["pdb_id"] == p),
                               entry_issues=";".join(entry_issue), **d["info"]))
        for eid, e in d["ents"].items():
            if e["type"].startswith("polypeptide"):
                full_fa.append((f"{p}_{eid}|asyms={','.join(e['asyms'])}|{e['fasta_header']}", e["can"]))

    # representative copy per (entry, heavy V sequence)
    for key in {(s["pdb_id"], s["variable_sequence"]) for s in samples}:
        grp = [s for s in samples if (s["pdb_id"], s["variable_sequence"]) == key]
        best = sorted(grp, key=lambda s: (-s["n_valid_backbone_atoms"], s["label_asym_id"]))[0]
        for s in grp:
            s["is_representative"] = s is best
    # duplicate sequence groups across entries
    groups = {}
    for s in samples:
        groups.setdefault(s["variable_sequence"], []).append(s["sample_id"])
    gid = {seq: f"seqgrp{i+1:04d}" for i, seq in enumerate(sorted(groups))}
    for s in samples:
        s["sequence_group"] = gid[s["variable_sequence"]]
        s["sequence_group_size"] = len(groups[s["variable_sequence"]])

    # 4. per-sample outputs
    map_rows, miss_rows = [], []
    for s in samples:
        L = s["L"]
        sd = os.path.join(sdir, s["sample_id"])
        os.makedirs(sd, exist_ok=True)
        coords = np.full((L, 4, 3), np.nan)
        occ = np.full((L, 4), np.nan)
        mask = np.zeros((L, 4), bool)
        for i, r in enumerate(s["rows"]):
            for j, n in enumerate(BACKBONE):
                a = r["atoms"][n]
                occ[i, j] = a["occ"]
                if a["valid"]:
                    coords[i, j] = (a["x"], a["y"], a["z"])
                    mask[i, j] = True
                else:
                    miss_rows.append(dict(sample_id=s["sample_id"], seq_idx1=i + 1, chothia_pos=r["chothia_pos"],
                                          label_seq_id=r["label_seq_id"], auth_seq_id=r["auth_seq_id"],
                                          insertion_code=r["insertion_code"], aa=r["aa"], atom=n, reason=a["reason"],
                                          occupancy="" if math.isnan(a["occ"]) else a["occ"],
                                          annotated_in_pdbx_unobs_or_zero_occ=a["annotated"]))
            map_rows.append(dict(sample_id=s["sample_id"], full_seq_idx1=r["full_seq_idx1"], seq_idx1=r["seq_idx1"],
                                 aa=r["aa"], chothia_pos=r["chothia_pos"], region=r["region"], entity_id=r["entity_id"],
                                 label_asym_id=r["label_asym_id"], auth_asym_id=r["auth_asym_id"],
                                 label_seq_id=r["label_seq_id"], auth_seq_id=r["auth_seq_id"],
                                 insertion_code=r["insertion_code"], model_id=r["model_id"],
                                 original_resname=r["original_resname"], seq_mon_ids=r["seq_mon_ids"],
                                 residue_normalization=r["residue_normalization"],
                                 renumbered_res_id=r["seq_idx1"], export_chain_id="H",
                                 coordinate_record_present=r["coordinate_record_present"],
                                 atom_mask="".join("1" if r["atoms"][n]["valid"] else "0" for n in BACKBONE),
                                 backbone_mask=all(r["atoms"][n]["valid"] for n in BACKBONE),
                                 altloc_choice=r["altloc_choice"], mapping_method=r["mapping_method"],
                                 mapping_status=r["mapping_status"]))
        np.savez_compressed(
            os.path.join(sd, "backbone.npz"), sequence=np.array(s["variable_sequence"]), coords=coords,
            atom_mask=mask, backbone_mask=mask.all(1),
            coordinate_record_present=np.array([r["coordinate_record_present"] for r in s["rows"]]),
            chothia_positions=np.array([r["chothia_pos"] for r in s["rows"]]),
            regions=np.array([r["region"] for r in s["rows"]]), seq_idx1=np.arange(1, L + 1),
            full_seq_idx1=np.array([r["full_seq_idx1"] for r in s["rows"]]),
            label_seq_id=np.array([r["label_seq_id"] for r in s["rows"]]),
            auth_seq_id=np.array([r["auth_seq_id"] for r in s["rows"]]),
            auth_insertion_code=np.array([r["insertion_code"] for r in s["rows"]]),
            occupancy=occ, atom_names=np.array(BACKBONE))
        s["export_check"] = write_backbone_cif(os.path.join(sd, "backbone_1based.cif"), s, "H")
        if s["export_check"] != "ok":
            s["status"], s["status_reason"] = "review", (s["status_reason"] + " | export:" + s["export_check"]).strip(" |")
        var_fa.append((f"{s['sample_id']}|{s['domain_type']}|full_idx1={s['domain_start_full_idx1']}-{s['domain_end_full_idx1']}",
                       s["variable_sequence"]))

    # 5. dataset Chothia alignment (formal samples)
    formal = [s for s in samples if s["status"] == "formal"]
    from abnumber.position import Position
    cols = sorted({Position.from_string(r["chothia_pos"], "H", "chothia") for s in formal for r in s["rows"]})
    col_names = [c.format() for c in cols]
    cidx = {c: i for i, c in enumerate(col_names)}
    M = len(cols)
    aln_rows = []
    seq_mask = np.zeros((len(formal), M), bool)
    col2idx = np.zeros((len(formal), M), int)  # 0 = sequence_gap, else seq_idx1
    for k, s in enumerate(formal):
        aa = ["-"] * M
        stt = ["-"] * M  # '-' sequence_gap, 'M' missing_residue, 'A' missing_atom/zero-occupancy, 'C' complete backbone
        for r in s["rows"]:
            j = cidx[r["chothia_pos"]]
            aa[j] = r["aa"]
            seq_mask[k, j] = True
            col2idx[k, j] = r["seq_idx1"]
            nv = sum(r["atoms"][n]["valid"] for n in BACKBONE)
            stt[j] = "C" if nv == 4 else ("M" if not r["coordinate_record_present"] else "A")
        aln_rows.append(dict(sample_id=s["sample_id"], aligned_sequence="".join(aa), structure_status="".join(stt)))
    if formal:
        np.savez_compressed(os.path.join(out, "chothia_alignment.npz"), columns=np.array(col_names),
                            sample_ids=np.array([s["sample_id"] for s in formal]), sequence_mask=seq_mask,
                            col_to_seq_idx1=col2idx)
    pd.DataFrame(aln_rows).to_csv(os.path.join(out, "chothia_alignment.tsv"), sep="\t", index=False)
    pd.DataFrame(dict(column_index=range(1, M + 1), chothia_pos=col_names,
                      region=[c.get_region() if hasattr(c, "get_region") else "" for c in cols])).to_csv(
        os.path.join(out, "chothia_columns.tsv"), sep="\t", index=False)

    # 6. tables
    meta_cols = [k for k in samples[0] if k not in ("rows", "full_chain_sequence")] if samples else []
    pd.DataFrame([{k: s[k] for k in meta_cols} for s in samples]).to_csv(os.path.join(out, "metadata.tsv"), sep="\t", index=False)
    pd.DataFrame(entry_rows).to_csv(os.path.join(out, "entries.tsv"), sep="\t", index=False)
    pd.DataFrame(map_rows).to_csv(os.path.join(out, "residue_mapping.tsv.gz"), sep="\t", index=False)
    pd.DataFrame(miss_rows).to_csv(os.path.join(out, "missing_atoms.tsv.gz"), sep="\t", index=False)
    for s in samples:
        if s["status"] != "formal":
            failed.append(dict(stage="sample_qc", error_type="review", reason=s["status_reason"], pdb_id=s["pdb_id"],
                               sample_id=s["sample_id"], attempted_fix=""))
    pd.DataFrame(failed, columns=["stage", "error_type", "reason", "pdb_id", "sample_id", "attempted_fix"]).to_csv(
        os.path.join(out, "failed_or_review.tsv"), sep="\t", index=False)
    with open(os.path.join(out, "full_chain.fasta"), "w") as f:
        f.writelines(f">{h}\n{s}\n" for h, s in full_fa)
    with open(os.path.join(out, "variable.fasta"), "w") as f:
        f.writelines(f">{h}\n{s}\n" for h, s in var_fa)
    pd.DataFrame(dl_log).to_csv(os.path.join(out, "download_log.tsv"), sep="\t", index=False)
    run["finished"] = now()
    json.dump(run, open(os.path.join(out, "run_config.json"), "w"), indent=1)

    # 7. QC
    import qc
    qc.run(out, samples, formal, entry_rows, failed, ids, aln_rows, col_names, status)


if __name__ == "__main__":
    main()
