#!/usr/bin/env python3
"""Two figures (ROC curve, score distributions) and a self-contained HTML report."""
from __future__ import annotations

import argparse
import base64
import csv
import html
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def roc_points(rows):
    """rows sorted best-first with activity 0/1. Returns (fpr, tpr)."""
    n_pos = sum(r["activity"] == 1 for r in rows)
    n_neg = len(rows) - n_pos
    tp = fp = 0
    fpr, tpr = [0.0], [0.0]
    for r in rows:
        tp += r["activity"] == 1
        fp += r["activity"] == 0
        fpr.append(fp / n_neg)
        tpr.append(tp / n_pos)
    return fpr, tpr


def load(path):
    out = []
    for r in csv.DictReader(open(path)):
        if r["activity"] in ("0", "1"):
            r["activity"] = int(r["activity"])
            r["affinity"] = float(r["affinity"])
            out.append(r)
    return sorted(out, key=lambda r: r["affinity"])


def figures(rows, metrics, outdir):
    fpr, tpr = roc_points(rows)
    fig, ax = plt.subplots(figsize=(4.6, 4.4), dpi=150)
    ax.plot(fpr, tpr, color="#1f6f8b", lw=2, label=f"docking score (AUC {metrics['roc_auc']:.2f})")
    ax.plot([0, 1], [0, 1], "--", color="#999", lw=1, label="random")
    ax.set(xlabel="false positive rate", ylabel="true positive rate",
           title=f"CA-II: {metrics['n_actives']} actives vs {metrics['n_decoys']} decoys")
    ax.legend(loc="lower right", frameon=False)
    fig.tight_layout()
    fig.savefig(Path(outdir) / "roc_curve.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.4, 4.0), dpi=150)
    bins = [x / 2 for x in range(-30, 0)]
    ax.hist([r["affinity"] for r in rows if r["activity"] == 0], bins=bins, alpha=.65, color="#999",
            label="property-matched decoys")
    ax.hist([r["affinity"] for r in rows if r["activity"] == 1], bins=bins, alpha=.65, color="#c8553d",
            label="known CA-II actives")
    ax.set(xlabel="best smina score (kcal/mol, lower = better)", ylabel="molecules",
           title="Score distributions")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(Path(outdir) / "score_distribution.png")
    plt.close(fig)


def table(rows, cols):
    head = "".join(f"<th>{html.escape(c)}</th>" for c in cols)
    body = "".join("<tr>" + "".join(f"<td>{html.escape(str(r.get(c, '')))}</td>" for c in cols) + "</tr>"
                   for r in rows)
    return f"<table><tr>{head}</tr>{body}</table>"


def img(path):
    return "data:image/png;base64," + base64.b64encode(Path(path).read_bytes()).decode()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ranked", required=True)
    ap.add_argument("--metrics", required=True)
    ap.add_argument("--redock", required=True)
    ap.add_argument("--annotated")
    ap.add_argument("--top")
    ap.add_argument("--top-nonactive")
    ap.add_argument("--note", default="")
    ap.add_argument("--out", default=".")
    a = ap.parse_args(argv)
    metrics = json.load(open(a.metrics))
    redock = json.load(open(a.redock))
    rows = load(a.ranked)
    figures(rows, metrics, a.out)
    cols = ["rank", "ligand", "affinity", "activity", "lipinski_violations", "pains", "herg_prob"]
    rd = redock
    ef = ", ".join(f"EF{float(k):.0%} = {v:.2f} (95 % CI {metrics['ef_ci95'][k][0]}–{metrics['ef_ci95'][k][1]}, "
                   f"max {metrics['ef_max_possible'][k]})" for k, v in metrics["ef"].items()) \
        if metrics.get("ef_ci95") else ", ".join(f"EF{k}={v}" for k, v in metrics["ef"].items())
    ci = metrics.get("roc_auc_ci95")
    doc = f"""<!doctype html><meta charset="utf-8"><title>CA-II virtual screening report</title>
<style>body{{font:15px/1.5 system-ui,sans-serif;max-width:980px;margin:2rem auto;padding:0 1rem;color:#222}}
table{{border-collapse:collapse;font-size:13px}}td,th{{border:1px solid #ccc;padding:2px 6px;text-align:left}}
img{{max-width:48%}}code{{background:#f2f2f2;padding:0 3px}}</style>
<h1>Carbonic anhydrase II — docking-based virtual screening</h1>
<p><b>Question.</b> Does the pipeline rank known CA-II binders above property-matched decoys?</p>
<h2>Result</h2>
<p>{metrics['n_actives']} actives and {metrics['n_decoys']} decoys docked
({metrics['n_scored_with_label']} scored). ROC-AUC <b>{metrics['roc_auc']:.3f}</b>{f" (95 % CI {ci[0]}–{ci[1]})" if ci else ""}.
{ef}.</p>
<p>Acetazolamide re-docking (top-ranked pose, {rd['n_seeds']} seeds): RMSD &lt; {rd['threshold']} Å in
{rd['top_pose_pass_fraction']:.0%} of seeds (min {rd['top_pose_rmsd_min']}, median {rd['top_pose_rmsd_median']},
max {rd['top_pose_rmsd_max']} Å). The scoring function and box were tuned on this same crystal ligand, so this
figure is optimistic.</p>
<p>{html.escape(a.note)}</p>
<p><img src="{img(Path(a.out)/'roc_curve.png')}"> <img src="{img(Path(a.out)/'score_distribution.png')}"></p>
<h2>Top 20 (Lipinski ≤ 1 violation, no PAINS alert)</h2>{table(list(csv.DictReader(open(a.top)))[:20], cols) if a.top else ''}
<h2>Top 20 that are not known actives</h2>{table(list(csv.DictReader(open(a.top_nonactive)))[:20], cols) if a.top_nonactive else ''}
<h2>Preparation failures</h2><p>{html.escape(', '.join(metrics.get('prep_failed', [])) or 'none')}</p>
<p>Nextflow <a href="pipeline_info/report.html">report</a> · <a href="pipeline_info/timeline.html">timeline</a></p>"""
    Path(a.out, "report.html").write_text(doc)
    print("wrote report.html, roc_curve.png, score_distribution.png")


if __name__ == "__main__":
    main()
