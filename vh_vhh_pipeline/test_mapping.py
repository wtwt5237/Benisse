"""Regression tests for the residue mapping core.

All data here is SYNTHETIC (minimal hand-made residues), used only to check mapping logic.
Run: python -m pytest vh_vhh_pipeline/test_mapping.py -q
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mapping import map_chain, select_conformer, BACKBONE  # noqa: E402

SEQ = ["GLN", "VAL", "GLN", "LEU", "SER", "GLY", "TYR"]  # synthetic 7-residue polymer


def poly(seq=SEQ):
    return {i + 1: [m] for i, m in enumerate(seq)}


def atoms(num, comp, auth=None, ins=None, names=BACKBONE, alt=None, occ=1.0, shift=0.0):
    auth = num if auth is None else auth
    return [dict(label_seq_id=num, comp=comp, atom=n, alt=alt, occ=occ, x=num + shift, y=k * 1.0, z=0.0,
                 auth_seq_id=str(auth), ins_code=ins, group="ATOM") for k, n in enumerate(names)]


def build(nums, seq=SEQ, **kw):
    rows = []
    for n in nums:
        rows += atoms(n, seq[n - 1], **kw)
    return rows


def test_n_terminal_missing():
    out, issues = map_chain(poly(), {}, build([3, 4, 5, 6, 7]))
    assert not issues
    assert [out[n]["coordinate_record_present"] for n in range(1, 8)] == [False, False, True, True, True, True, True]
    assert out[1]["atoms"]["CA"]["reason"] == "missing_residue" and math.isnan(out[1]["atoms"]["CA"]["x"])
    assert out[3]["atoms"]["CA"]["x"] == 3  # residue 3 keeps sequence position 3, not 1


def test_internal_consecutive_missing():
    out, _ = map_chain(poly(), {}, build([1, 2, 6, 7]))
    assert [n for n in range(1, 8) if not out[n]["coordinate_record_present"]] == [3, 4, 5]
    assert out[6]["original_resname"] == "GLY" and out[6]["atoms"]["N"]["x"] == 6


def test_single_backbone_atom_missing():
    rows = build([1, 2, 3, 5, 6, 7]) + atoms(4, "LEU", names=("N", "CA", "C"))
    out, _ = map_chain(poly(), {}, rows)
    a = out[4]["atoms"]
    assert a["N"]["valid"] and a["CA"]["valid"] and a["C"]["valid"]
    assert not a["O"]["valid"] and a["O"]["reason"] == "missing_atom"
    assert out[4]["coordinate_record_present"]


def test_oxt_does_not_replace_o():
    rows = build([1, 2, 3, 4, 5, 6]) + atoms(7, "TYR", names=("N", "CA", "C", "OXT"))
    out, _ = map_chain(poly(), {}, rows)
    assert out[7]["atoms"]["O"]["reason"] == "missing_atom"


def test_author_insertion_codes_and_compressed_author_numbering():
    # author numbering 100,100A,100B,... ; mapping must still follow label_seq_id
    rows = []
    for n in range(1, 8):
        rows += atoms(n, SEQ[n - 1], auth=100, ins=None if n == 1 else chr(64 + n - 1))
    scheme = {n: dict(auth_seq_num="100", ins_code=None if n == 1 else chr(64 + n - 1)) for n in range(1, 8)}
    out, _ = map_chain(poly(), scheme, rows)
    assert [out[n]["insertion_code"] for n in range(1, 8)] == ["", "A", "B", "C", "D", "E", "F"]
    assert all(out[n]["mapping_status"] == "ok" for n in out)


def test_structure_renumbered_from_one_does_not_shift_sequence():
    # observed residues 3..7 carry author numbers 1..5 (compressed); label_seq_id is authoritative
    rows = []
    for k, n in enumerate([3, 4, 5, 6, 7]):
        rows += atoms(n, SEQ[n - 1], auth=k + 1)
    out, _ = map_chain(poly(), {}, rows)
    assert out[3]["auth_seq_id"] == "1" and out[3]["original_resname"] == "GLN"
    assert not out[1]["coordinate_record_present"]


def test_author_numbering_conflict_with_scheme_is_flagged():
    scheme = {n: dict(auth_seq_num=str(n), ins_code=None) for n in range(1, 8)}
    rows = build([1, 2, 3, 4, 5, 6]) + atoms(7, "TYR", auth=99)
    out, _ = map_chain(poly(), scheme, rows)
    assert out[7]["mapping_status"] == "auth_scheme_mismatch"


def test_residue_mismatch_flagged():
    rows = build([1, 2, 3, 4, 5, 6]) + atoms(7, "TRP")
    out, _ = map_chain(poly(), {}, rows)
    assert out[7]["mapping_status"] == "mismatch"


def test_zero_occupancy_is_not_valid():
    rows = build([1, 2, 3, 4, 5, 6]) + atoms(7, "TYR", occ=0.0)
    out, _ = map_chain(poly(), {}, rows)
    assert all(out[7]["atoms"][n]["reason"] == "zero_occupancy" for n in BACKBONE)
    assert out[7]["coordinate_record_present"]


def test_altloc_whole_residue_choice_no_mixing():
    a = atoms(4, "LEU", alt="A", occ=0.4, shift=0.0)
    b = atoms(4, "LEU", alt="B", occ=0.6, shift=0.5)
    label, _ = select_conformer(a + b)
    assert label == "B"
    out, _ = map_chain(poly(), {}, build([1, 2, 3, 5, 6, 7]) + a + b)
    xs = {out[4]["atoms"][n]["x"] for n in BACKBONE}
    assert xs == {4.5}  # every backbone atom from conformer B


def test_altloc_tie_is_alphabetical():
    a = atoms(4, "LEU", alt="A", occ=0.5)
    b = atoms(4, "LEU", alt="B", occ=0.5, shift=0.5)
    assert select_conformer(b + a)[0] == "A"


def test_microheterogeneity():
    p = poly()
    p[5] = ["SER", "THR"]
    rows = build([1, 2, 3, 4, 6, 7]) + atoms(5, "THR")
    out, _ = map_chain(p, {}, rows)
    assert out[5]["mapping_status"] == "ok_microheterogeneity" and out[5]["original_resname"] == "THR"


def test_repeated_segment_does_not_confuse_explicit_mapping():
    # sequence with a repeat GSGSGS: positions come from label_seq_id, so the repeat cannot cause ambiguity
    seq = ["GLY", "SER", "GLY", "SER", "GLY", "SER"]
    out, _ = map_chain(poly(seq), {}, build([1, 2, 5, 6], seq=seq))
    assert [n for n in out if not out[n]["coordinate_record_present"]] == [3, 4]
