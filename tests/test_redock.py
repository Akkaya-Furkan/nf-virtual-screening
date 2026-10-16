import json

import pytest
import redock_summary as rs

rdkit = pytest.importorskip("rdkit")
from rdkit import Chem  # noqa: E402
from rdkit.Chem import AllChem  # noqa: E402

import rmsd_check as rc  # noqa: E402
from conftest import PDB  # noqa: E402

AZM = "CC(=O)NC1=NN=C(S1)S(=O)(=O)N"


@pytest.fixture
def ref(tmp_path):
    import prep_receptor as pr
    pr.main([str(PDB), "--out", str(tmp_path)])
    return rc.load_reference(str(tmp_path / "ref_ligand.pdb"), AZM), tmp_path


def test_reference_has_the_expected_graph(ref):
    mol, _ = ref
    assert mol.GetNumHeavyAtoms() == 13
    assert Chem.MolToSmiles(mol) == Chem.MolToSmiles(Chem.MolFromSmiles(AZM))


def test_identical_pose_has_zero_rmsd_and_matching_zn_distance(ref):
    mol, tmp = ref
    pose = Chem.Mol(mol)
    pose.SetProp("minimizedAffinity", "-7.0")
    rows, ref_zn = rc.evaluate(mol, [pose], rc.zn_coords(str(tmp / "receptor_clean.pdb")))
    assert rows[0]["rmsd"] == 0.0 and rows[0]["zn_N_distance"] == ref_zn
    assert 1.8 < ref_zn < 2.1  # sulfonamide N is bound to the zinc in the crystal


def test_translated_pose_rmsd_equals_shift(ref):
    mol, tmp = ref
    pose = Chem.Mol(mol)
    conf = pose.GetConformer()
    for i in range(pose.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (p.x + 3.0, p.y, p.z))
    rows, _ = rc.evaluate(mol, [pose], rc.zn_coords(str(tmp / "receptor_clean.pdb")))
    assert rows[0]["rmsd"] == pytest.approx(3.0, abs=0.01)


def test_symmetric_atoms_swapped_does_not_inflate_rmsd(ref):
    mol, _ = ref
    pose = Chem.Mol(mol)
    patt = Chem.MolFromSmarts("S(=O)(=O)N")
    s, o1, o2, _n = pose.GetSubstructMatch(patt)
    conf = pose.GetConformer()
    p1, p2 = conf.GetAtomPosition(o1), conf.GetAtomPosition(o2)
    p1, p2 = (p1.x, p1.y, p1.z), (p2.x, p2.y, p2.z)
    conf.SetAtomPosition(o1, p2)
    conf.SetAtomPosition(o2, p1)
    assert rc.rdMolAlign.CalcRMS(pose, mol) == pytest.approx(0.0, abs=1e-6)


def test_template_mismatch_raises(ref):
    _, tmp = ref
    with pytest.raises(ValueError):
        rc.load_reference(str(tmp / "ref_ligand.pdb"), "CCO")


def test_summary_counts_seeds():
    mk = lambda top, best: {"threshold": 2.0, "top_pose_rmsd": top, "best_of_modes_rmsd": best,
                            "crystal_zn_N_distance": 1.94}
    s = rs.summarize({1: mk(1.5, 1.5), 2: mk(5.0, 1.8), 3: mk(2.5, 2.2), 4: mk(1.9, 1.9)})
    assert s["top_pose_pass_fraction"] == 0.5 and s["best_of_modes_pass_fraction"] == 0.75
    assert s["top_pose_rmsd_min"] == 1.5 and s["top_pose_rmsd_max"] == 5.0
