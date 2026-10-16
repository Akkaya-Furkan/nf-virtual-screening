nextflow.enable.dsl = 2

/*
 * Docking-based virtual screening against human carbonic anhydrase II.
 *   PREP_RECEPTOR -> receptor PDBQT + reference ligand + search box
 *   PREP_LIGANDS  -> ligprep-cli: names/CIDs/SMILES -> PDBQT
 *   DOCK          -> smina, one task per ligand
 *   RANK          -> ranking, ROC-AUC, EF
 *   REDOCK        -> co-crystal ligand re-docking RMSD over several seeds
 */

process PREP_RECEPTOR {
    input:
    path pdb

    output:
    path 'receptor.pdbqt',       emit: receptor
    path 'receptor_clean.pdb',   emit: clean_pdb
    path 'ref_ligand.pdb',       emit: ref_ligand
    path 'box.json',             emit: box

    script:
    """
    prep_receptor.py ${pdb} --chain ${params.chain} --ligand ${params.ref_ligand} \\
        --resnum ${params.ref_resnum} --padding ${params.box_padding} --out .
    # polar hydrogens + AD4 atom types; the crystal structure carries no hydrogens
    obabel receptor_clean.pdb -xr -h -p 7.4 -O receptor.pdbqt
    # Open Babel leaves the zinc uncharged; set its formal +2 so the file is self-consistent
    # (Vina-type scoring ignores partial charges, so this does not change the scores)
    sed -i 's/^\\(ATOM .* ZN .*\\)+0.000 Zn/\\1+2.000 Zn/' receptor.pdbqt
    grep -q '+2.000 Zn' receptor.pdbqt
    """
}

process PREP_LIGANDS {
    input:
    path ligand_list

    output:
    path 'lig/*.pdbqt',      emit: pdbqt, optional: true
    path 'lig/manifest.csv', emit: manifest

    script:
    // ligprep exits 1 when some ligands fail; the manifest records which, RANK reports them,
    // and the pipeline stops only if nothing was prepared (exit 2 = bad input stays fatal)
    """
    ligprep --input ${ligand_list} --out lig || [ \$? -eq 1 ]
    ls lig/*.pdbqt > /dev/null
    """
}

process DOCK {
    tag "${pdbqt.baseName}"

    input:
    tuple path(pdbqt), path(receptor), path(box)

    output:
    tuple val(ligand), path("${ligand}.sdf"), emit: poses
    path "${ligand}.score.csv",               emit: score

    script:
    ligand = pdbqt.baseName.replaceFirst(/_CID\d+$/, '')
    """
    read cx cy cz sx sy sz < <(python3 -c "import json;b=json.load(open('${box}'));print(*b['center'],*b['size'])")
    smina -r ${receptor} -l ${pdbqt} \\
        --center_x \$cx --center_y \$cy --center_z \$cz --size_x \$sx --size_y \$sy --size_z \$sz \\
        --scoring ${params.scoring} --exhaustiveness ${params.exhaustiveness} \\
        --num_modes ${params.num_modes} --seed ${params.seed} --cpu ${task.cpus} \\
        -o ${ligand}.sdf --log ${ligand}.log > /dev/null
    echo "ligand,affinity" > ${ligand}.score.csv
    awk -v l=${ligand} '\$1 == 1 { print l "," \$2; exit }' ${ligand}.log >> ${ligand}.score.csv
    test \$(wc -l < ${ligand}.score.csv) -eq 2
    """
}

process RANK {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path scores
    path labels
    path manifest

    output:
    path 'ranked.csv'
    path 'metrics.json'

    script:
    """
    (echo "ligand,affinity"; tail -q -n +2 ${scores}) > all_scores.csv
    rank.py --scores all_scores.csv --labels ${labels} --manifest ${manifest} --out .
    """
}

process REDOCK {
    tag "seed ${seed}"

    input:
    tuple val(seed), path(pdbqt), path(receptor), path(clean_pdb), path(ref_ligand), path(box)

    output:
    path "redock_seed${seed}.json"

    script:
    """
    read cx cy cz sx sy sz < <(python3 -c "import json;b=json.load(open('${box}'));print(*b['center'],*b['size'])")
    smina -r ${receptor} -l ${pdbqt} \\
        --center_x \$cx --center_y \$cy --center_z \$cz --size_x \$sx --size_y \$sy --size_z \$sz \\
        --scoring ${params.scoring} --exhaustiveness ${params.exhaustiveness} \\
        --num_modes ${params.num_modes} --seed ${seed} --cpu ${task.cpus} \\
        -o poses.sdf > /dev/null
    rmsd_check.py --ref ${ref_ligand} --smiles '${params.redock_smiles}' --poses poses.sdf \\
        --receptor ${clean_pdb} --threshold ${params.rmsd_threshold} --out redock_seed${seed}.json
    """
}

process REDOCK_SUMMARY {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path jsons

    output:
    path 'redock_summary.json'

    script:
    """
    redock_summary.py --glob 'redock_seed*.json' --out redock_summary.json
    """
}

workflow {
    receptor = PREP_RECEPTOR(file(params.pdb))
    lig      = PREP_LIGANDS(file(params.ligands))

    dock_in = lig.pdbqt.flatten()
        .combine(receptor.receptor)
        .combine(receptor.box)
    docked = DOCK(dock_in)

    RANK(docked.score.collect(), file(params.labels), lig.manifest)

    // re-docking of the co-crystal ligand through the same prepared files
    azm = lig.pdbqt.flatten().filter { it.name.startsWith("${params.redock_label}_CID") }
    seeds = Channel.of(params.redock_seeds.toString().tokenize(',')*.toInteger()).flatten()
    redock_in = seeds.combine(azm)
        .combine(receptor.receptor).combine(receptor.clean_pdb)
        .combine(receptor.ref_ligand).combine(receptor.box)
    REDOCK_SUMMARY(REDOCK(redock_in).collect())
}
