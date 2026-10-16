import json

import prep_receptor as pr
import pytest
from conftest import PDB


def lines():
    return PDB.read_text().splitlines()


def test_split_keeps_active_site_ligand_zinc_and_protein_only():
    receptor, lig = pr.split_structure(lines(), "A", "AZM", 701)
    names = {l[17:20].strip() for l in receptor}
    assert "ZN" in names and "HOH" not in names and "GOL" not in names and "AZM" not in names
    assert len(lig) == 13 and {l[22:26].strip() for l in lig} == {"701"}


def test_only_first_altloc_is_kept():
    receptor, _ = pr.split_structure(lines(), "A", "AZM", 701)
    assert {l[16] for l in receptor} <= {" ", "A"}
    assert not any(l[16] == "B" for l in receptor)


def test_other_azm_copies_are_not_the_reference():
    _, lig = pr.split_structure(lines(), "A", "AZM", 702)
    xs = [float(l[30:38]) for l in lig]
    assert max(xs) < -15  # copy 702 sits on the protein surface, far from the zinc


def test_box_centre_and_size():
    center, size = pr.box_from_ligand(
        ["ATOM      1  C   AZM A   1       0.000   0.000   0.000",
         "ATOM      2  C   AZM A   1       4.000   2.000   6.000"], padding=1.0)
    assert center == [2.0, 1.0, 3.0] and size == [6.0, 4.0, 8.0]


def test_box_needs_atoms():
    with pytest.raises(ValueError):
        pr.box_from_ligand([], 4.0)


def test_cli_writes_files_and_zinc_is_inside_the_box(tmp_path):
    pr.main([str(PDB), "--out", str(tmp_path)])
    box = json.loads((tmp_path / "box.json").read_text())
    zn = [l for l in (tmp_path / "receptor_clean.pdb").read_text().splitlines() if " ZN " in l][0]
    z = [float(zn[30:38]), float(zn[38:46]), float(zn[46:54])]
    assert all(abs(z[i] - box["center"][i]) < box["size"][i] / 2 for i in range(3))


def test_missing_ligand_exits(tmp_path):
    with pytest.raises(SystemExit):
        pr.main([str(PDB), "--resnum", "999", "--out", str(tmp_path)])
