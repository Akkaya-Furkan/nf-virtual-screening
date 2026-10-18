#!/usr/bin/env python3
"""Compare several pipeline result directories (one per scoring variant) and two trivial baselines.

Usage: compare_runs.py results_vinardo results_vina ... > comparison.md
Baselines are computed on the molecules docked in the first directory: heavy-atom count and
molecular weight as 'scores' (bigger = more active), to show how much of an AUC is plain size.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from rank import auc_roc, bootstrap_ci, enrichment_factor  # noqa: E402
from rdkit import Chem  # noqa: E402
from rdkit.Chem import Descriptors  # noqa: E402


def main(dirs):
    print("| run | docked (act/decoy) | ROC-AUC (95 % CI) | EF1 % (95 % CI) | EF10 % | re-dock top pose < 2 Å | median RMSD | acetazolamide rank |")
    print("| --- | --- | --- | --- | --- | --- | --- | --- |")
    for d in dirs:
        d = Path(d)
        m = json.load(open(d / "metrics.json"))
        s = json.load(open(d / "redock_summary.json"))
        rows = list(csv.DictReader(open(d / "ranked.csv")))
        az = next((r["rank"] for r in rows if r["ligand"] == "1986"), "-")
        ci = m["roc_auc_ci95"]
        e1, e10 = m["ef"]["0.01"], m["ef"]["0.1"]
        c1 = m["ef_ci95"]["0.01"]
        print(f"| {d.name.replace('results_', '')} | {m['n_actives']}/{m['n_decoys']} | {m['roc_auc']:.3f} "
              f"({ci[0]}–{ci[1]}) | {e1:.2f} ({c1[0]}–{c1[1]}) | {e10:.2f} | "
              f"{s['top_pose_pass_fraction']:.0%} of {s['n_seeds']} seeds | {s['top_pose_rmsd_median']} | {az}/{len(rows)} |")
    rows = list(csv.DictReader(open(Path(dirs[0]) / "annotated.csv")))
    y = [int(r["activity"]) for r in rows]
    print("\nSize-only baselines on the same molecules:\n")
    print("| score | ROC-AUC (95 % CI) | EF1 % | EF10 % |")
    print("| --- | --- | --- | --- |")
    mols = [Chem.MolFromSmiles(r["smiles"]) for r in rows]
    for name, f in (("heavy atoms", lambda m: m.GetNumHeavyAtoms()), ("molecular weight", Descriptors.MolWt)):
        sc = [f(m) for m in mols]
        order = sorted(range(len(sc)), key=lambda i: -sc[i])
        yy = [y[i] for i in order]
        ci = bootstrap_ci(sc, y, None)
        print(f"| {name} | {auc_roc(sc, y):.3f} ({ci[0]}–{ci[1]}) | {enrichment_factor(yy, .01):.2f} | {enrichment_factor(yy, .1):.2f} |")
    herg = [r for r in rows if r["herg_flag"] == "1"]
    print(f"\nhERG flag (P3 model, p >= 0.5): {len(herg)} of {len(rows)} docked molecules.")


if __name__ == "__main__":
    main(sys.argv[1:])
