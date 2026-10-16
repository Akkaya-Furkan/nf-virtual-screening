#!/usr/bin/env python3
"""Aggregate re-docking JSON files (one per seed) into one table and a summary."""
from __future__ import annotations

import argparse
import glob
import json
import re


def summarize(results):
    """results: {seed: result dict from rmsd_check}. Returns summary dict."""
    seeds = sorted(results)
    top = [results[s]["top_pose_rmsd"] for s in seeds]
    best = [results[s]["best_of_modes_rmsd"] for s in seeds]
    thr = results[seeds[0]]["threshold"]
    return {
        "threshold": thr,
        "n_seeds": len(seeds),
        "seeds": seeds,
        "top_pose_rmsd_per_seed": dict(zip(seeds, top)),
        "top_pose_pass_fraction": sum(r < thr for r in top) / len(top),
        "best_of_modes_pass_fraction": sum(r < thr for r in best) / len(best),
        "top_pose_rmsd_min": min(top),
        "top_pose_rmsd_median": sorted(top)[len(top) // 2],
        "top_pose_rmsd_max": max(top),
        "crystal_zn_N_distance": results[seeds[0]]["crystal_zn_N_distance"],
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--glob", default="redock_seed*.json")
    ap.add_argument("--out", default="redock_summary.json")
    a = ap.parse_args(argv)
    results = {}
    for f in glob.glob(a.glob):
        results[int(re.search(r"seed(\d+)", f).group(1))] = json.load(open(f))
    if not results:
        raise SystemExit("no re-docking results found")
    s = summarize(results)
    json.dump(s, open(a.out, "w"), indent=2)
    print(f"re-docking over {s['n_seeds']} seeds: top-ranked pose RMSD < {s['threshold']} Å in "
          f"{s['top_pose_pass_fraction']:.0%} of seeds (min {s['top_pose_rmsd_min']}, median "
          f"{s['top_pose_rmsd_median']}, max {s['top_pose_rmsd_max']}); best-of-modes pass "
          f"{s['best_of_modes_pass_fraction']:.0%}")


if __name__ == "__main__":
    main()
