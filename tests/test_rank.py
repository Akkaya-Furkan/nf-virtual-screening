import csv

import pytest
import rank as rk


def test_auc_perfect_inverted_and_tied():
    assert rk.auc_roc([3, 2, 1, 0], [1, 1, 0, 0]) == 1.0
    assert rk.auc_roc([0, 1, 2, 3], [1, 1, 0, 0]) == 0.0
    assert rk.auc_roc([1, 1, 1, 1], [1, 1, 0, 0]) == 0.5
    assert rk.auc_roc([1, 2], [1, 1]) is None


def test_auc_hand_computed():
    # actives 0.9, 0.4; decoys 0.8, 0.3 -> pairs won: (0.9>0.8)(0.9>0.3)(0.4>0.3) = 3 of 4
    assert rk.auc_roc([0.9, 0.4, 0.8, 0.3], [1, 1, 0, 0]) == 0.75


def test_ef():
    labels = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # 1 active in 10
    assert rk.enrichment_factor(labels, 0.1) == 10.0
    assert rk.enrichment_factor([0, 1] + [0] * 8, 0.1) == 0.0
    assert rk.enrichment_factor([0, 0], 0.5) is None


def test_rank_orders_most_negative_first():
    rows = rk.rank([{"ligand": "a", "affinity": "-5"}, {"ligand": "b", "affinity": "-8"}], {"a": 0, "b": 1})
    assert [r["ligand"] for r in rows] == ["b", "a"] and rows[0]["rank"] == 1


def test_cli_reports_unlabeled_missing_and_failed(tmp_path):
    (tmp_path / "s.csv").write_text("ligand,affinity\nx,-7\ny,-5\nz,-6\n")
    (tmp_path / "l.csv").write_text("label,activity\nx,1\ny,0\nw,1\n")
    (tmp_path / "m.csv").write_text("input,status\nx,ok\nq,failed\n")
    rk.main(["--scores", str(tmp_path / "s.csv"), "--labels", str(tmp_path / "l.csv"),
             "--manifest", str(tmp_path / "m.csv"), "--out", str(tmp_path)])
    import json
    m = json.loads((tmp_path / "metrics.json").read_text())
    assert m["unlabeled_ligands"] == ["z"] and m["labelled_but_not_docked"] == ["w"]
    assert m["prep_failed"] == ["q"] and m["roc_auc"] == 1.0
    assert [r["ligand"] for r in csv.DictReader(open(tmp_path / "ranked.csv"))] == ["x", "z", "y"]
