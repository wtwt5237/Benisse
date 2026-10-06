"""Core sequence-to-structure mapping for one polymer chain (no I/O).

Everything here is keyed on the explicit mmCIF correspondence
(_atom_site.label_seq_id <-> _entity_poly_seq.num), never on list positions,
author numbering offsets or zip(sequence, residues).
"""
import math

BACKBONE = ("N", "CA", "C", "O")
BLANK_ALT = {None, False, "", ".", "?"}


def is_blank(v):
    return v in BLANK_ALT


def select_conformer(atoms):
    """Residue-level altloc choice.

    atoms: list of dicts with keys atom, alt, occ.
    Returns (label or None, reason). Shared (blank-altloc) atoms are always kept;
    exactly one altloc label is added. Never mixes labels atom by atom.
    """
    labels = sorted({a["alt"] for a in atoms if not is_blank(a["alt"])})
    if not labels:
        return None, "no_altloc"

    def score(label):
        kept = [a for a in atoms if is_blank(a["alt"]) or a["alt"] == label]
        bb = {a["atom"]: a["occ"] for a in kept if a["atom"] in BACKBONE}
        n_valid = sum(1 for o in bb.values() if o is not None and o > 0)
        occ_sum = sum(o for o in bb.values() if o is not None)
        return (n_valid, occ_sum)

    best = max(labels, key=lambda l: (score(l), [-ord(c) for c in l]))
    return best, f"altloc {best} of {labels}: max(backbone atoms with occ>0, backbone occupancy), ties alphabetical"


def map_chain(poly_seq, scheme_rows, atom_rows, unobs_res=None, unobs_atoms=None):
    """Map every polymer sequence position of one chain (label_asym, one model).

    poly_seq: dict num(int) -> list of allowed mon_ids (>1 = microheterogeneity),
              from _entity_poly_seq.
    scheme_rows: dict num -> dict(auth_seq_num, ins_code, strand_id, mon_id)
              from _pdbx_poly_seq_scheme (may be empty).
    atom_rows: list of dicts with keys label_seq_id(int|None), comp, atom, alt, occ,
              x, y, z, auth_seq_id, ins_code, group.
    unobs_res: set of label_seq_id annotated as unobserved/zero-occupancy residues.
    unobs_atoms: set of (label_seq_id, atom_name) annotated as unobserved/zero-occ atoms.

    Returns dict num -> residue record (one per poly_seq num) and list of chain-level issues.
    """
    unobs_res = unobs_res or set()
    unobs_atoms = unobs_atoms or set()
    issues = []
    by_num = {}
    for a in atom_rows:
        n = a["label_seq_id"]
        if n is None:
            issues.append(("atom_without_label_seq_id", f"{a['comp']} {a['atom']} auth {a['auth_seq_id']}"))
            continue
        if n not in poly_seq:
            issues.append(("label_seq_id_not_in_entity_poly_seq", str(n)))
            continue
        by_num.setdefault(n, []).append(a)

    out = {}
    for num in sorted(poly_seq):
        allowed = poly_seq[num]
        rec = {
            "label_seq_id": num,
            "seq_mon_ids": ",".join(allowed),
            "coordinate_record_present": False,
            "original_resname": "",
            "auth_seq_id": "",
            "insertion_code": "",
            "altloc_choice": "",
            "atoms": {},           # atom -> dict(x,y,z,occ,valid,reason)
            "mapping_method": "mmcif_label_seq_id",
            "mapping_status": "ok",
            "annotated_unobserved": num in unobs_res,
        }
        sch = scheme_rows.get(num)
        atoms = by_num.get(num, [])
        if atoms:
            rec["coordinate_record_present"] = True
            label, why = select_conformer(atoms)
            kept = [a for a in atoms if is_blank(a["alt"]) or a["alt"] == label]
            rec["altloc_choice"] = "" if label is None else why
            comps = sorted({a["comp"] for a in kept})
            auths = sorted({(str(a["auth_seq_id"]), "" if is_blank(a["ins_code"]) else a["ins_code"]) for a in kept})
            if len(comps) != 1:
                rec["mapping_status"] = "ambiguous_mapping"
                rec["original_resname"] = ",".join(comps)
            else:
                rec["original_resname"] = comps[0]
                if comps[0] not in allowed:
                    rec["mapping_status"] = "mismatch"
                elif len(allowed) > 1:
                    rec["mapping_status"] = "ok_microheterogeneity"
            if len(auths) != 1:
                rec["mapping_status"] = "ambiguous_mapping"
            else:
                rec["auth_seq_id"], rec["insertion_code"] = auths[0]
                if sch is not None and sch.get("auth_seq_num") not in (None, False, "?", "."):
                    s_ins = "" if is_blank(sch.get("ins_code")) else sch["ins_code"]
                    if (str(sch["auth_seq_num"]), s_ins) != auths[0]:
                        rec["mapping_status"] = "auth_scheme_mismatch"
            names_seen = {}
            for a in kept:
                if a["atom"] in BACKBONE:
                    names_seen.setdefault(a["atom"], []).append(a)
            for name in BACKBONE:
                lst = names_seen.get(name, [])
                if len(lst) > 1:
                    rec["mapping_status"] = "ambiguous_mapping"
                if not lst:
                    rec["atoms"][name] = dict(x=math.nan, y=math.nan, z=math.nan, occ=math.nan, valid=False,
                                              reason="missing_atom",
                                              annotated=(num, name) in unobs_atoms or num in unobs_res)
                    continue
                a = lst[0]
                finite = all(math.isfinite(v) for v in (a["x"], a["y"], a["z"]))
                occ = a["occ"] if a["occ"] is not None else math.nan
                if not finite:
                    valid, reason = False, "non_finite_coordinate"
                elif not (occ > 0):
                    valid, reason = False, "zero_occupancy"
                else:
                    valid, reason = True, ""
                rec["atoms"][name] = dict(x=a["x"], y=a["y"], z=a["z"], occ=occ, valid=valid, reason=reason,
                                          annotated=(num, name) in unobs_atoms or num in unobs_res)
        else:
            if sch is not None:
                # author numbering of an unobserved residue, as given by the file (not guessed)
                if sch.get("auth_seq_num") not in (None, False, "?", "."):
                    rec["auth_seq_id"] = str(sch["auth_seq_num"])
                if not is_blank(sch.get("ins_code")):
                    rec["insertion_code"] = sch["ins_code"]
            for name in BACKBONE:
                rec["atoms"][name] = dict(x=math.nan, y=math.nan, z=math.nan, occ=math.nan, valid=False,
                                          reason="missing_residue", annotated=num in unobs_res)
        out[num] = rec
    return out, issues
