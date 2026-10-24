# Nextflow virtual screening — carbonic anhydrase II

**Question.** Can a reproducible docking pipeline rank known carbonic anhydrase II (CA-II) binders above
property-matched decoys? Acetazolamide is the positive control.

**Status.** Complete. The pipeline docked a 568-molecule CA-II library (142 ChEMBL actives, 426 property-matched
decoys; 538 docked) under five scoring variants. **In this setup the ranking is close to chance** (best ROC-AUC 0.622,
vinardo); details in Results and Limitations.

## Structure
PDB **3HS4**, human CA-II with acetazolamide, X-ray, **1.10 Å** (checked against the RCSB API on 2026-10-16:
`resolution_combined` 1.1, bound components AZM and ZN). Chain A, first alternate location only, zinc kept,
waters/glycerol removed. The file has three AZM copies; only residue 701 binds the zinc (sulfonamide N to Zn
1.94 Å in the crystal), 702 and 703 are partial-occupancy surface copies and are ignored.

## Pipeline
```
PREP_RECEPTOR  PDB -> clean chain + Zn, reference ligand, box, receptor.pdbqt (Open Babel, pH 7.4)
PREP_LIGANDS   names/CIDs/SMILES -> PDBQT with ligprep-cli (PubChem 3D, validated; manifest.csv)
DOCK           smina, one task per ligand            -> best affinity per ligand
RANK           ranked.csv, ROC-AUC, EF               -> results/
REDOCK         acetazolamide through the same files, once per seed; RMSD to the crystal pose
REDOCK_SUMMARY per-seed RMSDs and the pass fraction  -> results/redock_summary.json
```
Defaults live in `nextflow.config` (scoring `vinardo`, exhaustiveness 8, 9 modes, box = crystal ligand extent + 4 Å).

## Run
Needs Nextflow >= 25.04 and Java. Pick one profile:
```bash
nextflow run . -profile conda                        # builds envs/nfvs.yml on first use
nextflow run . -profile docker                       # image built from docker/Dockerfile (not tested here)
nextflow run . -profile local                        # smina, obabel, ligprep, RDKit already on PATH
nextflow run . -profile conda --ligands my.txt --labels my_labels.csv --scoring vina
```
Inputs: a ligprep-cli list (`data/trial10.txt`) and a labels CSV `label,activity` (1 active, 0 decoy; the label is
the name as written in the list). `results/` gets `ranked.csv`, `metrics.json`, `redock_summary.json` and
Nextflow's `pipeline_info/` report, timeline and trace.

`ligprep-cli` is a private repository, so the conda environment installs it over git and needs credentials
(`gh auth setup-git`). Open Babel is pinned to 3.1.1: conda-forge's 3.2.1 aborted with "double free" on exit
(2026-10-16), which made ligprep report every ligand as failed. Fedora sets `BABEL_DATADIR`; the config unsets it
so the conda Open Babel finds its own data.

## Tests
```bash
pytest -q        # 18 tests, no network: box/receptor cleaning on 3HS4, RMSD, ROC-AUC, EF, seed summary
```
Not covered by tests: the Nextflow workflow itself (checked by running it) and smina.

## What was measured on Day 1 (smoke test and re-docking tuning)
All numbers below are from runs on 2026-10-16 (smina Nov 2017 build from bioconda; desktop for the final run).

**Smoke test** (`data/trial10.txt`: 5 CA sulfonamides — acetazolamide, methazolamide, dorzolamide,
brinzolamide, ethoxzolamide — and 5 unrelated drugs — aspirin, caffeine, ibuprofen, metformin, paracetamol).
The decoys are **not** property-matched and n = 10, so these numbers show that the pipeline works, not how well
docking separates actives. Final run, defaults: 24 of 24 tasks succeeded in 1 m 43 s; ROC-AUC 0.66; the ranking was
brinzolamide, ibuprofen, aspirin, ethoxzolamide, methazolamide, dorzolamide, paracetamol, **acetazolamide (8th of
10)**, caffeine, metformin. Scoring `vina` gave AUC 0.9 in the same 10 molecules; with n = 10 neither value is
distinguishable from noise and the scoring function is not chosen on this basis.

**Acetazolamide re-docking** (ligprep-cli PDBQT, flexible ligand, rigid receptor; heavy-atom, symmetry-aware RMSD
against the crystal ligand, no superposition). Criterion: RMSD < 2 Å for the top-ranked pose. Final default
configuration (vinardo, box +4 Å), 10 seeds: top-ranked pose **1.76–1.89 Å in 9 of 10 seeds** (1.76, 1.76, 1.89
×7) and **5.81 Å in one seed**; best of the 9 modes < 2 Å in 10 of 10 seeds. This meets the criterion, with the
caveats in the next paragraph. Settings explored before that (same crystal ligand, so this is tuning on the test
case):

| scoring | box padding | top-pose RMSD < 2 Å |
| --- | --- | --- |
| vina | 8 Å | 0 of 1 seed (5.39 Å; exhaustiveness 8 and 32 alike; best of 9 modes 1.95 Å) |
| vina | 4 Å | 0 of 10 seeds (5.38–5.39 Å) |
| ad4_scoring | 4 and 8 Å | 0 of 1 seed each (5.33, 5.31 Å) |
| vinardo | 8 Å | 2 of 10 seeds (1.76 ×2, 2.47 ×5, 5.81 ×3); exhaustiveness 64 did not fix it |
| **vinardo** | **4 Å** | **9 of 10 seeds** |

**Caveats on that result.** (1) The box is centred on the crystal ligand and kept tight; this is the standard
re-docking protocol, but a real screen has no such hint and pose quality there will not be better. (2) The padding
and scoring were chosen because they re-dock this one ligand better; that is optimistic by construction. (3) With
the 8 Å box the top three vinardo poses differ by 0.005 kcal/mol (−5.773, −5.772, −5.768), so which one is
ranked first is decided by the random seed; the correct pose is always among them. Re-docking one ligand says
nothing about the other CA-II ligands.

## Results (full screen, 5 scoring variants)
Library: `scripts/build_library.py` (ChEMBL human CA-II actives, pChEMBL >= 6, Butina-clustered; decoys are approved
drugs with no recorded CA activity and no primary sulfonamide, matched to each active on MW, cLogP, HBD, HBA,
rotatable bonds, charge and heavy atoms, max ECFP4 Tanimoto to any active < threshold). 142 actives + 426 decoys were
built; ligands that failed preparation are excluded from the metrics. Tables below are generated, not typed:
`results/summary.csv` / `results/summary.md` by `scripts/make_summary.py` and `results/comparison.md` by
`scripts/compare_runs.py`, both from each run's `metrics.json` and `redock_summary.json`. Confidence intervals are
bootstrap 95 %.

Docked 129 actives / 409 decoys; 30 ligands failed preparation.

| variant | ROC-AUC (95 % CI) | EF1 % (95 % CI) | EF10 % (95 % CI) | re-dock top pose < 2 Å | best of modes < 2 Å | median top-pose RMSD (Å) |
| --- | --- | --- | --- | --- | --- | --- |
| vina | 0.548 (0.491–0.606) | 1.39 (0.0–3.297) | 1.39 (0.879–1.873) | 0% of 10 seeds | 90% | 5.38 |
| vinardo | 0.622 (0.569–0.675) | 1.39 (0.0–2.514) | 1.85 (1.266–2.317) | 90% of 10 seeds | 100% | 1.89 |
| vinardo_like | 0.588 (0.535–0.642) | 2.09 (0.0–3.371) | 1.31 (0.754–1.736) | 0% of 10 seeds | 100% | 5.41 |
| zn1 | 0.537 (0.485–0.591) | 1.39 (0.0–3.449) | 1.00 (0.553–1.445) | 0% of 10 seeds | 100% | 5.41 |
| zn2 | 0.465 (0.413–0.517) | 1.39 (0.0–3.321) | 0.62 (0.241–1.09) | 0% of 10 seeds | 20% | 6.19 |

Maximum possible EF1 % with this class balance: 4.171.

Size-only baselines (no docking, same 538 molecules):

| score | ROC-AUC (95 % CI) | EF1 % | EF10 % |
| --- | --- | --- | --- |
| heavy atoms | 0.475 (0.419–0.531) | 0.70 | 0.93 |
| molecular weight | 0.559 (0.504–0.615) | 2.09 | 1.54 |

Reading: ROC-AUC 0.5 is chance. The best variant (vinardo, 0.622, CI 0.569–0.675) is only weakly above it, and the
scoring function that re-docks acetazolamide best is the one with the highest AUC only because it was chosen on
that re-docking (see Caveats). EF1 % has CIs that start at 0 for every variant (1 % of 538 is ~5 molecules).
zn1/zn2 are the experimental zinc-term variants (`scoring/vinardo_like_zn*.txt`, a heuristic Zn–N reward); they did not help (zn2 AUC 0.465).
This is not a pipeline failure: the pipeline ran 538/538 ligands to a score and the metrics reproduce from the files.
It is a finding that **docking scores in this setup barely separate CA-II binders from property-matched decoys.**

## Limitations
- **Zinc.** The catalytic Zn²⁺ is bound by His94, His96, His119 and, in the crystal, the sulfonamide nitrogen
  (1.94 Å). Vina-type scoring functions have no metal-coordination term, and the sulfonamide is docked neutral
  although the zinc-bound form is the anion. There is no physical zinc-coordination term. The experimental variants zn1/zn2 add a
  heuristic Gaussian reward for a donor nitrogen 2.0 Å from the zinc (`scoring/README.md`; weights arbitrary, not
  fitted) and did not improve the ranking.
  No Zn–N-restrained pose selection was done, so a sulfonamide pose that does not coordinate the zinc can still
  rank first.
- **Close to chance.** Best ROC-AUC 0.622 (vinardo); vina 0.548, zn2 0.465. EF1 % CIs include 0. In this setup the
  discrimination is near chance, not a demonstrated screening capability.
- **Size confound.** Heavy-atom count alone gives ROC-AUC 0.475 and molecular weight 0.559, i.e. size explains
  little here, but the gap between 0.559 and the best docking AUC is also small and the intervals overlap.
- **Preparation failures are not random.** 30 of 568 ligands failed preparation (13 of 142 actives = 9.2 %, 17 of
  426 decoys = 4.0 %; from `prep_failed` in `metrics.json` and `data/ca2_labels.csv`). Actives are lost more often,
  so the 538 docked molecules are not exactly the designed set.
- **Re-docking is circular.** The docking box and the re-docking reference are both the crystal acetazolamide, and
  box padding and scoring function were chosen because they re-dock that ligand (vinardo, 4 Å: 9 of 10 seeds < 2 Å;
  vina: 0 of 10). A real screen has no such hint. Acetazolamide itself ranks 455/538 with vinardo.
- Histidine protonation and hydrogens come from Open Babel at pH 7.4; the zinc-bound water is removed.
- Decoys are approved drugs, not experimentally confirmed non-binders; property-matching limits but does not remove
  bias. Docking scores (kcal/mol) are not affinities.
- Docker profile: not verified (no Docker on the development machines).
- CI runs only the pytest suite (no Nextflow, smina or Open Babel); the workflow itself is checked by running it.
