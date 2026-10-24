#!/usr/bin/env python3
"""Collect metrics.json / redock_summary.json of several result directories into results/summary.{csv,md}.

Usage: make_summary.py OUTDIR results_vina results_vinardo ...   (variant name = directory name minus 'results_')
Every number is read from the files; nothing is typed by hand.
"""
import csv
import json
import sys
from pathlib import Path

COLS = ["variant", "n_actives", "n_decoys", "prep_failed", "roc_auc", "roc_auc_ci_lo", "roc_auc_ci_hi",
        "ef1", "ef1_ci_lo", "ef1_ci_hi", "ef10", "ef10_ci_lo", "ef10_ci_hi", "ef_max_possible",
        "redock_top_pose_pass_fraction", "redock_best_of_modes_pass_fraction", "redock_n_seeds",
        "redock_rmsd_median", "redock_rmsd_max"]


def row(d):
    d = Path(d)
    m = json.load(open(d / "metrics.json"))
    s = json.load(open(d / "redock_summary.json"))
    return {
        "variant": d.name.replace("results_", ""), "n_actives": m["n_actives"], "n_decoys": m["n_decoys"],
        "prep_failed": len(m["prep_failed"]), "roc_auc": round(m["roc_auc"], 3),
        "roc_auc_ci_lo": m["roc_auc_ci95"][0], "roc_auc_ci_hi": m["roc_auc_ci95"][1],
        "ef1": round(m["ef"]["0.01"], 2), "ef1_ci_lo": m["ef_ci95"]["0.01"][0], "ef1_ci_hi": m["ef_ci95"]["0.01"][1],
        "ef10": round(m["ef"]["0.1"], 2), "ef10_ci_lo": m["ef_ci95"]["0.1"][0], "ef10_ci_hi": m["ef_ci95"]["0.1"][1],
        "ef_max_possible": m["ef_max_possible"]["0.01"],
        "redock_top_pose_pass_fraction": s["top_pose_pass_fraction"],
        "redock_best_of_modes_pass_fraction": s["best_of_modes_pass_fraction"], "redock_n_seeds": s["n_seeds"],
        "redock_rmsd_median": s["top_pose_rmsd_median"], "redock_rmsd_max": s["top_pose_rmsd_max"],
    }


def main(out, dirs):
    rows = [row(d) for d in dirs]
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, COLS)
        w.writeheader()
        w.writerows(rows)
    r0 = rows[0]
    lines = [f"Docked {r0['n_actives']} actives / {r0['n_decoys']} decoys; {r0['prep_failed']} ligands failed preparation.", "",
             "| variant | ROC-AUC (95 % CI) | EF1 % (95 % CI) | EF10 % (95 % CI) | re-dock top pose < 2 Å | best of modes < 2 Å | median top-pose RMSD (Å) |",
             "| --- | --- | --- | --- | --- | --- | --- |"]
    for r in rows:
        lines.append(f"| {r['variant']} | {r['roc_auc']:.3f} ({r['roc_auc_ci_lo']}–{r['roc_auc_ci_hi']}) "
                     f"| {r['ef1']:.2f} ({r['ef1_ci_lo']}–{r['ef1_ci_hi']}) | {r['ef10']:.2f} ({r['ef10_ci_lo']}–{r['ef10_ci_hi']}) "
                     f"| {r['redock_top_pose_pass_fraction']:.0%} of {r['redock_n_seeds']} seeds "
                     f"| {r['redock_best_of_modes_pass_fraction']:.0%} | {r['redock_rmsd_median']} |")
    lines += ["", f"Maximum possible EF1 % with this class balance: {r0['ef_max_possible']}."]
    (out / "summary.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
