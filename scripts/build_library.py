#!/usr/bin/env python3
"""Build the CA-II screening library from live ChEMBL and PubChem responses.

Actives   : ChEMBL human CA-II (CHEMBL205) Ki/IC50 records, '=' relation, nM, pChEMBL >= 6; one row per
            molecule (median pChEMBL); ECFP4 Butina clusters (Tanimoto distance 0.6), the best-pChEMBL member
            of the largest clusters first, until N actives (acetazolamide is always kept).
Pool      : ChEMBL approved (max_phase 4) small molecules, minus any molecule with a recorded activity on
            ANY human carbonic anhydrase and minus primary sulfonamides.
Decoys    : per active, the DECOYS_PER_ACTIVE nearest pool molecules in standardised
            (MW, cLogP, HBD, HBA, rotatable bonds, formal charge, heavy atoms) space, without reuse, with
            max ECFP4 Tanimoto to every active < MAX_SIM. Only molecules with a PubChem CID are kept
            (ligprep-cli needs one); an unresolvable decoy is replaced by the next nearest.
Every API response is cached under data/raw/ so the build is repeatable offline once cached.
"""
from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import requests
from rdkit import Chem, DataStructs
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdFingerprintGenerator
from rdkit.ML.Cluster import Butina

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
SESSION = requests.Session()
SESSION.headers["User-Agent"] = "nf-virtual-screening/0.2 (research)"
RAW = Path("data/raw")
PRIMARY_SULFONAMIDE = Chem.MolFromSmarts("[SX4](=O)(=O)[NX3;H2]")


def get_json(url, params=None, cache=None, retries=4):
    path = RAW / cache if cache else None
    if path and path.exists():
        return json.loads(path.read_text())
    for attempt in range(retries):
        r = SESSION.get(url, params=params, timeout=60)
        if r.status_code == 200:
            data = r.json()
            if path:
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(json.dumps(data))
            return data
        if r.status_code == 404:
            return None
        time.sleep(2 ** attempt)
    raise RuntimeError(f"{url} failed: HTTP {r.status_code}")


def chembl_all(resource, params, cache):
    """Follow ChEMBL pagination; cache the concatenated list."""
    path = RAW / cache
    if path.exists():
        return json.loads(path.read_text())
    out, url, p = [], f"{CHEMBL}/{resource}.json", dict(params, limit=1000)
    while url:
        data = get_json(url, p)
        key = [k for k in data if k != "page_meta"][0]
        out += data[key]
        nxt = data["page_meta"]["next"]
        url, p = (f"https://www.ebi.ac.uk{nxt}" if nxt else None), None
        print(f"  {resource}: {len(out)}/{data['page_meta']['total_count']}", file=sys.stderr)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out))
    return out


def ca_targets():
    data = get_json(f"{CHEMBL}/target.json", {"pref_name__icontains": "carbonic anhydrase",
                    "organism": "Homo sapiens", "limit": 100}, cache="chembl_ca_targets.json")
    return sorted(t["target_chembl_id"] for t in data["targets"] if t["target_type"] == "SINGLE PROTEIN")


def props(mol):
    return [Descriptors.MolWt(mol), Crippen.MolLogP(mol), Lipinski.NumHDonors(mol),
            Lipinski.NumHAcceptors(mol), Lipinski.NumRotatableBonds(mol),
            Chem.GetFormalCharge(mol), mol.GetNumHeavyAtoms()]


PROP_NAMES = ["mw", "clogp", "hbd", "hba", "rotb", "charge", "heavy"]


def clean_mol(smiles):
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        return None
    frags = Chem.GetMolFrags(mol, asMols=True)  # drop counter-ions: keep the largest fragment
    return max(frags, key=lambda m: m.GetNumHeavyAtoms())


def select_actives(acts, n, seed_chembl_id="CHEMBL20"):
    by_mol = defaultdict(list)
    smiles = {}
    for a in acts:
        if (a["standard_relation"] == "=" and a["standard_units"] == "nM"
                and a["standard_type"] in ("Ki", "IC50") and a.get("pchembl_value")
                and float(a["pchembl_value"]) >= 6 and not a.get("data_validity_comment")
                and a.get("canonical_smiles")):
            by_mol[a["molecule_chembl_id"]].append(float(a["pchembl_value"]))
            smiles[a["molecule_chembl_id"]] = a["canonical_smiles"]
    mols = {}
    for cid, s in smiles.items():
        m = clean_mol(s)
        if m is not None and 8 <= m.GetNumHeavyAtoms() <= 60:
            mols[cid] = m
    ids = sorted(mols)
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = [gen.GetFingerprint(mols[i]) for i in ids]
    dists = []
    for i in range(1, len(fps)):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend(1 - s for s in sims)
    clusters = sorted(Butina.ClusterData(dists, len(fps), 0.6, isDistData=True), key=len, reverse=True)
    chosen = []
    rank = 0
    while len(chosen) < n and any(len(c) > rank for c in clusters):
        # round-robin over clusters, best pChEMBL first inside each cluster
        for c in clusters:
            if len(c) > rank and len(chosen) < n:
                members = sorted(c, key=lambda k: -statistics.median(by_mol[ids[k]]))
                chosen.append(ids[members[rank]])
        rank += 1
    if seed_chembl_id in mols and seed_chembl_id not in chosen:
        chosen[-1] = seed_chembl_id
    return [(i, mols[i], statistics.median(by_mol[i]), len(by_mol[i]), len(clusters)) for i in chosen]


def pubchem_cid(inchikey, cache={}):
    if inchikey in cache:
        return cache[inchikey]
    time.sleep(0.25)
    d = get_json(f"{PUBCHEM}/compound/inchikey/{inchikey}/cids/JSON", cache=f"pubchem/{inchikey}.json")
    cids = (d or {}).get("IdentifierList", {}).get("CID", [])
    cache[inchikey] = cids[0] if cids and cids[0] > 0 else None
    return cache[inchikey]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--actives", type=int, default=150)
    ap.add_argument("--decoys-per-active", type=int, default=3)
    ap.add_argument("--max-sim", type=float, default=0.4)
    ap.add_argument("--out", default="data")
    a = ap.parse_args()

    targets = ca_targets()
    print("CA targets:", targets, file=sys.stderr)
    acts = chembl_all("activity", {"target_chembl_id": "CHEMBL205", "standard_type__in": "Ki,IC50",
                                   "pchembl_value__isnull": "false"}, "chembl_ca2_activities.json")
    any_ca = set()
    for t in targets:
        any_ca |= {x["molecule_chembl_id"] for x in chembl_all(
            "activity", {"target_chembl_id": t, "pchembl_value__isnull": "false"}, f"chembl_act_{t}.json")}
    drugs = chembl_all("molecule", {"max_phase": 4, "molecule_type": "Small molecule"}, "chembl_approved.json")
    print(f"CA-II activity rows {len(acts)}; molecules with any CA activity {len(any_ca)}; "
          f"approved small molecules {len(drugs)}", file=sys.stderr)

    actives = select_actives(acts, a.actives)
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    active_fps = [gen.GetFingerprint(m) for _, m, *_ in actives]
    ids_act = {i for i, *_ in actives}

    pool = []
    seen_blocks = {Chem.MolToInchiKey(m)[:14] for _, m, *_ in actives}  # no decoy may be an active's parent
    for d in drugs:
        cid = d["molecule_chembl_id"]
        st = d.get("molecule_structures") or {}
        if cid in any_ca or cid in ids_act or not st.get("canonical_smiles"):
            continue
        m = clean_mol(st["canonical_smiles"])
        if m is None or m.HasSubstructMatch(PRIMARY_SULFONAMIDE) or not 8 <= m.GetNumHeavyAtoms() <= 60:
            continue
        fp = gen.GetFingerprint(m)
        if max(DataStructs.BulkTanimotoSimilarity(fp, active_fps)) >= a.max_sim:
            continue
        key = Chem.MolToInchiKey(m)  # parent structure, not the salt ChEMBL may list
        if key[:14] in seen_blocks:
            continue
        seen_blocks.add(key[:14])
        pool.append((cid, d.get("pref_name") or cid, m, key, props(m)))
    print(f"decoy pool after exclusions: {len(pool)}", file=sys.stderr)

    P = np.array([p[4] for p in pool])
    sd = P.std(axis=0)
    sd[sd == 0] = 1
    used, rows, stats = set(), [], {"active_no_cid": 0, "decoy_unresolved": 0}
    for cid, m, pchembl, nrec, ncl in actives:
        key = Chem.MolToInchiKey(m)
        pcid = pubchem_cid(key)
        if pcid is None:
            stats["active_no_cid"] += 1
            continue
        pr = props(m)
        rows.append({"label": str(pcid), "activity": 1, "chembl_id": cid, "name": cid,
                     "pchembl": round(pchembl, 2), "matched_to": "", "smiles": Chem.MolToSmiles(m),
                     **dict(zip(PROP_NAMES, np.round(pr, 2)))})
        dist = np.sqrt((((P - np.array(pr)) / sd) ** 2).sum(axis=1))
        got = 0
        for j in np.argsort(dist):
            if got == a.decoys_per_active:
                break
            if j in used:
                continue
            dcid, name, dm, dkey, dp = pool[j]
            pc = pubchem_cid(dkey)
            if pc is None:
                stats["decoy_unresolved"] += 1
                used.add(j)
                continue
            used.add(j)
            got += 1
            rows.append({"label": str(pc), "activity": 0, "chembl_id": dcid, "name": name, "pchembl": "",
                         "matched_to": cid, "smiles": Chem.MolToSmiles(dm),
                         **dict(zip(PROP_NAMES, np.round(dp, 2)))})
    labels = {}
    for r in rows:  # a PubChem CID must not appear twice
        labels.setdefault(r["label"], r)
    rows = list(labels.values())

    out = Path(a.out)
    with open(out / "ca2_labels.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with open(out / "ca2_library.txt", "w") as fh:
        fh.write("# CA-II library built by scripts/build_library.py from ChEMBL + PubChem (cid: lines)\n")
        fh.writelines(f"cid:{r['label']}\n" for r in rows)
    n1 = sum(r["activity"] for r in rows)
    summary = {"chembl_ca2_activity_rows": len(acts), "ca_targets": targets,
               "approved_small_molecules": len(drugs), "decoy_pool": len(pool), "actives": n1,
               "decoys": len(rows) - n1, **stats}
    json.dump(summary, open(out / "ca2_library_summary.json", "w"), indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
