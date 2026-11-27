"""Small synthetic fixtures for configurable targets; no external data requests."""
import importlib.util
import json
from pathlib import Path

import pytest
import prep_receptor as pr
import rmsd_check as rc
import redock_summary as rs
import report
from rdkit import Chem
from rdkit.Chem import AllChem


def tiny_pdb(path):
    mol = Chem.MolFromSmiles('CCO')
    AllChem.Compute2DCoords(mol)
    text = Chem.MolToPDBBlock(mol)
    lines = [l[:17] + 'LIG' + ' A' + '   1' + l[26:] for l in text.splitlines()
             if l.startswith('HETATM')]
    protein = 'ATOM      9  CA  ALA A   2       9.000   9.000   9.000  1.00 20.00           C'
    path.write_text(protein + '\n' + '\n'.join(lines) + '\nEND\n')
    return mol


def test_prep_no_metal_and_default_error(tmp_path):
    pdb = tmp_path / 'tiny.pdb'
    tiny_pdb(pdb)
    args = [str(pdb), '--ligand', 'LIG', '--resnum', '1', '--out', str(tmp_path / 'out')]
    with pytest.raises(SystemExit, match='no zinc found in chain A'):
        pr.main(args)
    pr.main(args + ['--metal', 'none'])
    assert len((tmp_path / 'out/receptor_clean.pdb').read_text().splitlines()) == 2
    assert json.loads((tmp_path / 'out/box.json').read_text())['size']
    with pytest.raises(SystemExit, match='no MG found'):
        pr.main(args + ['--metal', 'MG'])


def test_other_metal_required_and_preserved(tmp_path):
    pdb = tmp_path / 'tiny.pdb'
    tiny_pdb(pdb)
    metal = 'HETATM   10 MG    MG A   3       1.000   2.000   3.000  1.00 20.00          MG'
    # Use fixed PDB residue columns rather than relying on atom-name spacing.
    metal = metal[:17] + ' MG' + metal[20:]
    pdb.write_text(pdb.read_text() + metal + '\n')
    pr.main([str(pdb), '--ligand', 'LIG', '--resnum', '1', '--metal', 'MG', '--out', str(tmp_path / 'out')])
    assert ' MG ' in (tmp_path / 'out/receptor_clean.pdb').read_text()


@pytest.mark.parametrize('metal', ['none', 'ZN'])
def test_rmsd_without_zinc_keeps_schema(tmp_path, metal):
    pdb = tmp_path / 'tiny.pdb'
    mol = tiny_pdb(pdb)
    poses = tmp_path / 'poses.sdf'
    with Chem.SDWriter(str(poses)) as w:
        w.write(mol)
    ref = tmp_path / 'ref.pdb'
    Chem.MolToPDBFile(mol, str(ref))
    out = tmp_path / 'redock.json'
    rc.main(['--ref', str(ref), '--smiles', 'CCO', '--poses', str(poses),
             '--receptor', str(pdb), '--metal', metal, '--out', str(out)])
    result = json.loads(out.read_text())
    assert result['top_pose_rmsd'] == 0.0
    assert result['crystal_zn_N_distance'] is None
    assert result['modes'][0]['zn_N_distance'] is None
    assert rs.summarize({1: result})['crystal_zn_N_distance'] is None


def test_report_parameters_and_default_decoy_label(tmp_path, monkeypatch):
    ranked = tmp_path / 'ranked.csv'
    ranked.write_text('rank,ligand,affinity,activity\n1,a,-7,1\n2,b,-6,0\n')
    metrics = tmp_path / 'metrics.json'
    metrics.write_text(json.dumps(dict(n_actives=1, n_decoys=1, n_scored_with_label=2,
                                     roc_auc=1., ef={'0.1': 2.})))
    labels = []
    original = report.plt.Axes.hist
    def hist(ax, *args, **kwargs):
        labels.append(kwargs.get('label'))
        return original(ax, *args, **kwargs)
    monkeypatch.setattr(report.plt.Axes, 'hist', hist)
    args = ['--ranked', str(ranked), '--metrics', str(metrics), '--out', str(tmp_path),
            '--report-target', 'Example <target>', '--report-question', 'Example question?']
    report.main(args)
    assert labels == ['decoys', 'known Example <target> actives']
    doc = (tmp_path / 'report.html').read_text()
    assert 'Example &lt;target&gt; —' in doc and 'Example question?' in doc
    assert 'Re-docking skipped' in doc and 'acetazolamide' not in doc.lower()
    labels.clear()
    rd = tmp_path / 'redock.json'
    rd.write_text(json.dumps(rs.summarize({1: dict(threshold=2., top_pose_rmsd=1.,
                      best_of_modes_rmsd=1., crystal_zn_N_distance=None)})))
    report.main(args + ['--redock', str(rd), '--redock-label', 'reference', '--decoy-label', 'matched controls'])
    assert labels[0] == 'matched controls'
    assert 'Reference re-docking' in (tmp_path / 'report.html').read_text()


def test_build_library_parser():
    spec = importlib.util.spec_from_file_location('build_library', Path(__file__).parents[1] / 'scripts/build_library.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    args = module.parse_args([])
    assert vars(args) == dict(actives=150, decoys_per_active=3, max_sim=.4, out='data', target_chembl_id='CHEMBL205')
    assert module.parse_args(['--target-chembl-id', 'CHEMBL203']).target_chembl_id == 'CHEMBL203'


def test_rmsd_restores_conversion_lost_bond_orders(capsys):
    ref = Chem.MolFromSmiles('c1ccncc1')
    AllChem.Compute2DCoords(ref)
    pose = Chem.Mol(ref)
    for atom in pose.GetAtoms():
        atom.SetIsAromatic(False)
    for bond in pose.GetBonds():
        bond.SetIsAromatic(False)
        bond.SetBondType(Chem.BondType.SINGLE)
    Chem.SanitizeMol(pose)
    pose.SetProp('minimizedAffinity', '-7.0')
    rows, distance = rc.evaluate(ref, [pose], None)
    assert rows[0]['rmsd'] == 0.0 and rows[0]['affinity'] == -7.0
    assert distance is None
    assert 'restored pose bond orders/charges' in capsys.readouterr().err
