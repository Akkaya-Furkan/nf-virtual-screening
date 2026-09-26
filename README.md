# Nextflow virtual screening — carbonic anhydrase II

**Question.** Can a reproducible docking pipeline rank known carbonic anhydrase II binders above property-matched decoys? (acetazolamide as positive control)

## Data
A high-resolution (< 2 Å) sulfonamide-bound CA-II structure from the PDB; ChEMBL approved drugs, known CA-II actives and property-matched decoys.

## Approach
1. Nextflow DSL2 with conda and Docker profiles
2. Processes: PREP_RECEPTOR, PREP_LIGANDS (ligprep-cli), DOCK (smina), RANK
3. Acetazolamide redocking (RMSD < 2 Å); enrichment (EF 1 %, ROC-AUC)
4. Lipinski + PAINS filters and a hERG flag from herg-cardiotox-classifier
5. HTML report with Nextflow report/timeline

## Deliverables
End-to-end run, enrichment curve, score distribution, report.

## Limitations
- CA-II has a catalytic zinc; docking scoring functions model metal coordination poorly, which can distort sulfonamide ranking.
- Decoy selection directly shapes enrichment metrics; the method is documented.

## Status
Planned — part of a 30-project computational biology and scientific software portfolio.
