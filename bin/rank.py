#!/usr/bin/env python3
"""Rank docked ligands and score the ranking against known activity labels."""
from __future__ import annotations

import argparse
import csv
import json
import math


def auc_roc(scores, labels):
    """ROC-AUC by the rank-sum (Mann-Whitney) formula; ties count 0.5. scores: higher = more active."""
    pos = [s for s, y in zip(scores, labels) if y == 1]
    neg = [s for s, y in zip(scores, labels) if y == 0]
    if not pos or not neg:
        return None
    wins = sum((p > n) + 0.5 * (p == n) for p in pos for n in neg)
    return wins / (len(pos) * len(neg))


def enrichment_factor(ranked_labels, fraction):
    """EF at a fraction of the ranked list: hit rate in the top fraction / overall hit rate.
    The top slice holds ceil(n * fraction) molecules (at least 1)."""
    n = len(ranked_labels)
    n_act = sum(ranked_labels)
    if n == 0 or n_act == 0:
        return None
    k = max(1, math.ceil(n * fraction))
    return (sum(ranked_labels[:k]) / k) / (n_act / n)


def bootstrap_ci(scores, labels, fraction, n_boot=1000, seed=0):
    """95 % percentile bootstrap (resampling molecules) for ROC-AUC (fraction=None) or EF at `fraction`."""
    import random
    rng = random.Random(seed)
    idx = range(len(scores))
    vals = []
    for _ in range(n_boot):
        s = [rng.choice(idx) for _ in idx]
        sc, y = [scores[i] for i in s], [labels[i] for i in s]
        if fraction is None:
            v = auc_roc(sc, y)
        else:
            order = sorted(range(len(sc)), key=lambda i: -sc[i])
            v = enrichment_factor([y[i] for i in order], fraction)
        if v is not None:
            vals.append(v)
    vals.sort()
    if not vals:
        return None
    return [round(vals[int(0.025 * len(vals))], 3), round(vals[int(0.975 * len(vals)) - 1], 3)]


def rank(score_rows, labels):
    """score_rows: dicts with ligand, affinity (kcal/mol, more negative = better). Returns ranked dicts."""
    rows = []
    for r in score_rows:
        rows.append({"ligand": r["ligand"], "affinity": float(r["affinity"]),
                     "activity": labels.get(r["ligand"])})
    rows.sort(key=lambda r: r["affinity"])
    for i, r in enumerate(rows, 1):
        r["rank"] = i
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scores", required=True, help="CSV: ligand,affinity")
    ap.add_argument("--labels", required=True, help="CSV: label,activity (1 active, 0 decoy)")
    ap.add_argument("--manifest", help="ligprep manifest.csv, to list ligands that failed preparation")
    ap.add_argument("--ef", default="0.01,0.1", help="comma-separated EF fractions")
    ap.add_argument("--out", default=".")
    a = ap.parse_args(argv)

    labels = {r["label"]: int(r["activity"]) for r in csv.DictReader(open(a.labels))}
    scores = list(csv.DictReader(open(a.scores)))
    ranked = rank(scores, labels)
    unlabeled = [r["ligand"] for r in ranked if r["activity"] is None]
    scored = [r for r in ranked if r["activity"] is not None]
    y = [r["activity"] for r in scored]
    metrics = {
        "n_docked": len(ranked),
        "n_scored_with_label": len(scored),
        "n_actives": sum(y),
        "n_decoys": len(y) - sum(y),
        "roc_auc": auc_roc([-r["affinity"] for r in scored], y),
        "ef": {a_: enrichment_factor(y, float(a_)) for a_ in a.ef.split(",")},
        "roc_auc_ci95": bootstrap_ci([-r["affinity"] for r in scored], y, None) if len(scored) > 20 else None,
        "ef_ci95": ({a_: bootstrap_ci([-r["affinity"] for r in scored], y, float(a_))
                     for a_ in a.ef.split(",")} if len(scored) > 20 else None),
        "ef_max_possible": {a_: round(min(1.0, sum(y) / max(1, math.ceil(len(y) * float(a_))))
                                      / (sum(y) / len(y)), 3) if sum(y) else None
                            for a_ in a.ef.split(",")},
        "unlabeled_ligands": unlabeled,
        "labelled_but_not_docked": sorted(set(labels) - {r["ligand"] for r in ranked}),
    }
    if a.manifest:
        metrics["prep_failed"] = [r["input"] for r in csv.DictReader(open(a.manifest))
                                  if r["status"] != "ok"]
    with open(f"{a.out}/ranked.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["rank", "ligand", "affinity", "activity"])
        w.writeheader()
        w.writerows(ranked)
    json.dump(metrics, open(f"{a.out}/metrics.json", "w"), indent=2)
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
