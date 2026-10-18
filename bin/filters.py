#!/usr/bin/env python3
"""Annotate a ranking with Lipinski, PAINS and (optionally) hERG flags and pick the top candidates.

SMILES come from the ligprep manifest. The hERG flag uses the P3 random-forest model
(Morgan r2 2048 bits + MW, cLogP, TPSA; label pIC50 >= 5): probability >= 0.5 is flagged.
"""
from __future__ import annotations

import argparse
import csv
import re

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams


def ligprep_label(entry):
    """Mirror of ligprep.inputs.Entry.label for the manifest `input` column."""
    text = entry.strip()
    prefix, sep, rest = text.partition(":")
    value = rest.strip() if sep and prefix.lower() in ("name", "cid", "smiles") else text
    base = "smiles" if sep and prefix.lower() == "smiles" else value
    return re.sub(r"[^A-Za-z0-9._-]+", "_", base).strip("_")[:40] or "ligand"


def lipinski_violations(mol):
    return sum([Descriptors.MolWt(mol) > 500, Crippen.MolLogP(mol) > 5,
                Lipinski.NumHDonors(mol) > 5, Lipinski.NumHAcceptors(mol) > 10])


def pains_catalog():
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    return FilterCatalog(params)


def pains_hits(catalog, mol):
    return [e.GetDescription() for e in catalog.GetMatches(mol)]


def herg_features(mol):
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    bits = list(gen.GetFingerprint(mol))
    return bits + [Descriptors.MolWt(mol), Crippen.MolLogP(mol), rdMolDescriptors.CalcTPSA(mol)]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ranked", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--herg-model", help="joblib file {'model','features'} from the P3 pipeline")
    ap.add_argument("--top", type=int, default=20)
    ap.add_argument("--out", default=".")
    a = ap.parse_args(argv)

    smiles = {ligprep_label(r["input"]): r["smiles"] for r in csv.DictReader(open(a.manifest)) if r["smiles"]}
    catalog = pains_catalog()
    model = None
    if a.herg_model:
        import joblib
        import numpy as np
        model = joblib.load(a.herg_model)["model"]
    rows = list(csv.DictReader(open(a.ranked)))
    mols = {}
    for r in rows:
        mol = Chem.MolFromSmiles(smiles.get(r["ligand"], ""))
        mols[r["ligand"]] = mol
    if model is not None:
        ok = [r["ligand"] for r in rows if mols[r["ligand"]] is not None]
        probs = model.predict_proba(np.array([herg_features(mols[k]) for k in ok], dtype=np.float32))[:, 1]
        herg = dict(zip(ok, probs))
    for r in rows:
        mol = mols[r["ligand"]]
        r["smiles"] = smiles.get(r["ligand"], "")
        if mol is None:
            r.update(lipinski_violations="", pains="", herg_prob="", herg_flag="")
            continue
        hits = pains_hits(catalog, mol)
        r["lipinski_violations"] = lipinski_violations(mol)
        r["pains"] = ";".join(hits)
        if model is not None:
            r["herg_prob"] = round(float(herg[r["ligand"]]), 3)
            r["herg_flag"] = int(herg[r["ligand"]] >= 0.5)
        else:
            r["herg_prob"] = r["herg_flag"] = ""
    fields = list(rows[0])
    with open(f"{a.out}/annotated.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    # top-N by score among molecules that pass Lipinski (<= 1 violation) and have no PAINS alert
    passing = [r for r in rows if r["lipinski_violations"] != "" and r["lipinski_violations"] <= 1 and not r["pains"]]
    for name, subset in (("top_candidates.csv", passing[: a.top]),
                         ("top_unlabelled_or_decoy.csv", [r for r in passing if r["activity"] != "1"][: a.top])):
        with open(f"{a.out}/{name}", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            w.writerows(subset)
    print(f"{len(rows)} annotated; {len(passing)} pass Lipinski(<=1) and PAINS-free")


if __name__ == "__main__":
    main()
