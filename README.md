# Nextflow virtual screening — carbonic anhydrase II

**Question.** Can a reproducible docking pipeline rank known carbonic anhydrase II (CA-II) binders above
property-matched decoys? Acetazolamide is the positive control.

**Status (Day 1 of 2).** The pipeline runs end to end on a 10-molecule smoke test. The full library, enrichment
metrics (EF 1 %, ROC-AUC on a real decoy set), Lipinski/PAINS/hERG filters and the HTML report are Day 2 and
are **not done**. Nothing below answers the question yet.

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

## What was measured on Day 1
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

## Limitations
- **Zinc.** The catalytic Zn²⁺ is bound by His94, His96, His119 and, in the crystal, the sulfonamide nitrogen
  (1.94 Å). Vina-type scoring functions have no metal-coordination term: the zinc is an atom with a radius, and
  the sulfonamide is docked neutral although the zinc-bound form is the deprotonated anion. The zinc contribution
  to the score is therefore unreliable, which can distort the ranking of sulfonamides against other
  chemotypes. In the 8 Å vinardo runs the correctly placed pose (1.76 Å) had a Zn–N distance of 2.3 Å against
  1.94 Å in the crystal; poses 2.5–5.8 Å away sat 3.3–6.9 Å from the zinc. No restraint or special zinc handling is applied.
- Histidine protonation and hydrogens come from Open Babel at pH 7.4; the crystal structure has no hydrogens.
  The zinc-bound water/hydroxide is removed with all other waters.
- Decoy selection directly shapes the enrichment metrics; the Day 2 method will be written here before the run.
- Docking scores (kcal/mol) are not affinities; small molecules score lower by construction (acetazolamide,
  13 heavy atoms, ranks below ibuprofen and aspirin in the smoke test).
- Docker profile is written but was not run (no Docker on the development machines).
