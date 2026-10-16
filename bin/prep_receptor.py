#!/usr/bin/env python3
"""Clean a PDB entry into a docking receptor, a reference ligand and a search box.

Keeps the chosen protein chain (first alternate location only) and the catalytic zinc;
drops waters, glycerol, other ligands and extra copies of the co-crystal ligand.
The box is centred on the reference ligand with `padding` Å around its extent.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

WATER = {"HOH", "WAT", "DOD"}


def _altloc_ok(line: str) -> bool:
    return line[16] in (" ", "A")


def split_structure(lines, chain, ligand, resnum, metal="ZN"):
    """Return (receptor_lines, ligand_lines) as PDB text lines."""
    receptor, lig = [], []
    for line in lines:
        rec = line[:6]
        if rec not in ("ATOM  ", "HETATM") or line[21] != chain or not _altloc_ok(line):
            continue
        resname = line[17:20].strip()
        element = line[76:78].strip().upper()
        if rec == "ATOM  ":
            receptor.append(line.rstrip("\n"))
        elif resname == ligand and int(line[22:26]) == resnum:
            lig.append(line.rstrip("\n"))
        elif resname == metal or (element == metal and resname not in WATER):
            receptor.append(line.rstrip("\n"))
    return receptor, lig


def box_from_ligand(ligand_lines, padding):
    xyz = [(float(l[30:38]), float(l[38:46]), float(l[46:54])) for l in ligand_lines]
    if not xyz:
        raise ValueError("reference ligand has no atoms")
    lo = [min(c[i] for c in xyz) for i in range(3)]
    hi = [max(c[i] for c in xyz) for i in range(3)]
    center = [round((a + b) / 2, 3) for a, b in zip(lo, hi)]
    size = [round(b - a + 2 * padding, 3) for a, b in zip(lo, hi)]
    return center, size


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("pdb")
    ap.add_argument("--chain", default="A")
    ap.add_argument("--ligand", default="AZM", help="co-crystal ligand residue name")
    ap.add_argument("--resnum", type=int, default=701, help="residue number of the copy in the active site")
    ap.add_argument("--padding", type=float, default=8.0)
    ap.add_argument("--out", default=".")
    a = ap.parse_args(argv)

    lines = Path(a.pdb).read_text().splitlines()
    receptor, lig = split_structure(lines, a.chain, a.ligand, a.resnum)
    if not lig:
        raise SystemExit(f"ligand {a.ligand} {a.resnum} chain {a.chain} not found")
    if not any(l[17:20].strip() == "ZN" for l in receptor):
        raise SystemExit("no zinc found in chain " + a.chain)
    center, size = box_from_ligand(lig, a.padding)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "receptor_clean.pdb").write_text("\n".join(receptor) + "\nEND\n")
    (out / "ref_ligand.pdb").write_text("\n".join(lig) + "\nEND\n")
    (out / "box.json").write_text(json.dumps({"center": center, "size": size}, indent=2) + "\n")
    print(f"receptor atoms {len(receptor)}, ligand atoms {len(lig)}, center {center}, size {size}")


if __name__ == "__main__":
    main()
