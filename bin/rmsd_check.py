#!/usr/bin/env python3
"""Re-docking check: heavy-atom RMSD of docked poses against the co-crystal ligand.

RMSD is symmetry-aware and computed in place (no superposition), because the docked pose lives
in the receptor frame. Also reports the distance from the zinc to the sulfonamide nitrogen, which
shows whether the pose reproduces the zinc binding that docking scoring does not model well.
"""
from __future__ import annotations

import argparse
import json
import sys

from rdkit import Chem
from rdkit.Chem import rdMolAlign
from rdkit.Chem.AllChem import AssignBondOrdersFromTemplate

SULFONAMIDE = Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(N-C=O)]S(=O)(=O)")


def load_reference(ref_pdb, template_smiles):
    ref = Chem.MolFromPDBFile(ref_pdb, removeHs=True, sanitize=True)
    tmpl = Chem.MolFromSmiles(template_smiles)
    if ref is None or tmpl is None:
        raise ValueError("cannot read reference ligand or template SMILES")
    if tmpl.GetNumHeavyAtoms() != ref.GetNumHeavyAtoms():
        raise ValueError(f"template has {tmpl.GetNumHeavyAtoms()} heavy atoms, "
                         f"reference ligand {ref.GetNumHeavyAtoms()}")
    return AssignBondOrdersFromTemplate(tmpl, ref)


def sulfonamide_n_index(mol):
    hit = mol.GetSubstructMatch(Chem.MolFromSmarts("[NX3]S(=O)(=O)c"))
    if not hit:
        hit = mol.GetSubstructMatch(Chem.MolFromSmarts("[NX3]S(=O)(=O)"))
    return hit[0] if hit else None


def zn_distance(mol, conf_id, zn_xyz, n_idx):
    if n_idx is None:
        return None
    p = mol.GetConformer(conf_id).GetAtomPosition(n_idx)
    return round(((p.x - zn_xyz[0]) ** 2 + (p.y - zn_xyz[1]) ** 2 + (p.z - zn_xyz[2]) ** 2) ** 0.5, 2)


def evaluate(ref, poses, zn_xyz):
    """poses: list of RDKit mols (one per docking mode, rank order). Returns list of dicts."""
    rows = []
    ref_n = sulfonamide_n_index(ref)
    for i, pose in enumerate(poses, 1):
        pose = Chem.RemoveHs(pose)
        rmsd = rdMolAlign.CalcRMS(pose, ref)
        rows.append({
            "mode": i,
            "affinity": float(pose.GetProp("minimizedAffinity")) if pose.HasProp("minimizedAffinity") else None,
            "rmsd": round(rmsd, 2),
            "zn_N_distance": zn_distance(pose, 0, zn_xyz, sulfonamide_n_index(pose)),
        })
    return rows, zn_distance(ref, 0, zn_xyz, ref_n)


def zn_coords(receptor_pdb):
    for line in open(receptor_pdb):
        if line.startswith(("ATOM", "HETATM")) and line[17:20].strip() == "ZN":
            return float(line[30:38]), float(line[38:46]), float(line[46:54])
    raise ValueError("no zinc in receptor")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ref", required=True, help="co-crystal ligand PDB")
    ap.add_argument("--smiles", required=True, help="SMILES of the ligand (bond-order template)")
    ap.add_argument("--poses", required=True, help="smina SDF output")
    ap.add_argument("--receptor", required=True, help="cleaned receptor PDB (for the zinc position)")
    ap.add_argument("--threshold", type=float, default=2.0)
    ap.add_argument("--out", default="redock.json")
    a = ap.parse_args(argv)

    ref = load_reference(a.ref, a.smiles)
    poses = [m for m in Chem.SDMolSupplier(a.poses, removeHs=False) if m is not None]
    if not poses:
        sys.exit("no poses in " + a.poses)
    rows, ref_zn = evaluate(ref, poses, zn_coords(a.receptor))
    top = rows[0]["rmsd"]
    best = min(r["rmsd"] for r in rows)
    result = {
        "threshold": a.threshold,
        "top_pose_rmsd": top,
        "best_of_modes_rmsd": best,
        "top_pose_pass": top < a.threshold,
        "best_of_modes_pass": best < a.threshold,
        "crystal_zn_N_distance": ref_zn,
        "modes": rows,
    }
    json.dump(result, open(a.out, "w"), indent=2)
    print(f"top-ranked pose RMSD {top} Å, best of {len(rows)} modes {best} Å "
          f"(threshold {a.threshold}); crystal Zn–N {ref_zn} Å")


if __name__ == "__main__":
    main()
