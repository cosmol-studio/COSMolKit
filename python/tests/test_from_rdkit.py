from pathlib import Path

import numpy as np
import pytest

import cosmolkit

Chem = pytest.importorskip("rdkit.Chem")
Point3D = pytest.importorskip("rdkit.Geometry").Point3D


def _load_smiles_cases():
    corpus = Path(__file__).resolve().parents[2] / "tests" / "smiles.smi"
    smiles = [
        line.strip()
        for line in corpus.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    assert smiles, f"no SMILES found in {corpus}"
    return smiles


SMILES_CASES = _load_smiles_cases()


def _rdkit_mol_or_skip(smiles):
    rd_mol = Chem.MolFromSmiles(smiles)
    if rd_mol is None:
        pytest.skip(f"RDKit cannot parse corpus SMILES: {smiles}")
    return rd_mol


def _topology_signature(mol):
    atoms = [
        (
            atom.idx(),
            atom.atomic_num(),
            atom.formal_charge(),
            atom.chiral_tag(),
            atom.isotope(),
            atom.atom_map_num(),
        )
        for atom in mol.atoms()
    ]
    bonds = [
        (
            bond.idx(),
            min(bond.begin_atom_idx(), bond.end_atom_idx()),
            max(bond.begin_atom_idx(), bond.end_atom_idx()),
            bond.bond_type(),
        )
        for bond in mol.bonds()
    ]
    return atoms, sorted(bonds)


def _feature_signature(mol):
    atoms = [
        (
            atom.idx(),
            atom.atomic_num(),
            atom.formal_charge(),
            atom.chiral_tag(),
            atom.isotope(),
            atom.atom_map_num(),
            atom.is_aromatic(),
            atom.explicit_hydrogens(),
            atom.no_implicit(),
            atom.num_radical_electrons(),
            atom.degree(),
            atom.explicit_valence(),
            atom.implicit_hydrogens(),
            atom.total_num_hs(),
            atom.total_valence(),
        )
        for atom in mol.atoms()
    ]
    bonds = [
        (
            bond.idx(),
            min(bond.begin_atom_idx(), bond.end_atom_idx()),
            max(bond.begin_atom_idx(), bond.end_atom_idx()),
            bond.bond_type(),
            bond.bond_dir(),
            bond.stereo(),
            tuple(bond.stereo_atoms()),
            bond.is_aromatic(),
        )
        for bond in mol.bonds()
    ]
    return atoms, sorted(bonds)


def _rdkit_signature(rd_mol):
    atoms = [
        (
            atom.GetIdx(),
            atom.GetAtomicNum(),
            atom.GetFormalCharge(),
            cosmolkit.CHIRAL_TAG_MAP[str(atom.GetChiralTag())],
            atom.GetIsotope() or None,
            atom.GetAtomMapNum() or None,
            atom.GetIsAromatic(),
            atom.GetNumExplicitHs(),
            atom.GetNoImplicit(),
            atom.GetNumRadicalElectrons(),
            atom.GetDegree(),
            atom.GetExplicitValence(),
            atom.GetNumImplicitHs(),
            atom.GetTotalNumHs(),
            atom.GetTotalValence(),
        )
        for atom in rd_mol.GetAtoms()
    ]
    bonds = [
        (
            bond.GetIdx(),
            min(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
            max(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
            cosmolkit.BOND_ORDER_MAP[
                "AROMATIC" if bond.GetIsAromatic() else str(bond.GetBondType())
            ],
            cosmolkit.BOND_DIRECTION_MAP[str(bond.GetBondDir())],
            cosmolkit.BOND_STEREO_MAP[str(bond.GetStereo())],
            tuple(bond.GetStereoAtoms()),
            bond.GetIsAromatic(),
        )
        for bond in rd_mol.GetBonds()
    ]
    return atoms, sorted(bonds)


@pytest.mark.parametrize("smiles", SMILES_CASES)
def test_from_rdkit_copies_basic_graph_features(smiles):
    rd_mol = _rdkit_mol_or_skip(smiles)

    direct = cosmolkit.Molecule.from_smiles(smiles)
    bridged = cosmolkit.Molecule.from_rdkit(rd_mol)

    assert _topology_signature(bridged) == _topology_signature(direct)


def test_atom_and_bond_feature_enums_are_intenum_values():
    mol = cosmolkit.Molecule.from_smiles("C=C")
    atom = mol.atoms()[0]
    bond = mol.bonds()[0]

    assert atom.chiral_tag() == cosmolkit.ChiralTag.CHI_UNSPECIFIED
    assert atom.chiral_tag_code() == int(cosmolkit.ChiralTag.CHI_UNSPECIFIED)
    assert bond.bond_type() == cosmolkit.BondOrder.DOUBLE
    assert bond.bond_type_code() == int(cosmolkit.BondOrder.DOUBLE)
    assert bond.bond_dir() == cosmolkit.BondDirection.NONE
    assert bond.stereo() == cosmolkit.BondStereo.NONE
    assert cosmolkit.BOND_ORDER_MAP["DOUBLE"] == cosmolkit.BondOrder.DOUBLE
    assert not hasattr(bond, "order")


def test_public_bond_enums_include_hydrogen_and_unknown_members():
    assert cosmolkit.BondOrder.HYDROGEN == 18
    assert cosmolkit.BOND_ORDER_MAP["HYDROGEN"] == cosmolkit.BondOrder.HYDROGEN
    assert cosmolkit.BondDirection.UNKNOWN == 6
    assert (
        cosmolkit.BOND_DIRECTION_MAP["UNKNOWN"] == cosmolkit.BondDirection.UNKNOWN
    )


@pytest.mark.parametrize("smiles", SMILES_CASES)
def test_from_rdkit_exposes_rdkit_basic_atom_and_bond_features(smiles):
    rd_mol = _rdkit_mol_or_skip(smiles)

    bridged = cosmolkit.Molecule.from_rdkit(rd_mol)

    assert _feature_signature(bridged) == _rdkit_signature(rd_mol)


@pytest.mark.parametrize("smiles", SMILES_CASES)
def test_from_rdkit_matches_direct_cosmolkit_smiles(smiles):
    rd_mol = _rdkit_mol_or_skip(smiles)

    direct = cosmolkit.Molecule.from_smiles(smiles)
    bridged = cosmolkit.Molecule.from_rdkit(rd_mol)

    assert _feature_signature(bridged) == _feature_signature(direct)


def _add_conformer(rd_mol, coords, is_3d):
    conf = Chem.Conformer(rd_mol.GetNumAtoms())
    conf.Set3D(is_3d)
    for idx, (x, y, z) in enumerate(coords):
        conf.SetAtomPosition(idx, Point3D(float(x), float(y), float(z)))
    rd_mol.AddConformer(conf, assignId=True)


def test_from_rdkit_copies_3d_conformers():
    rd_mol = Chem.MolFromSmiles("CCO")
    coords = np.array(
        [
            [0.1, 0.2, 0.3],
            [1.1, 1.2, 1.3],
            [2.1, 2.2, 2.3],
        ]
    )
    _add_conformer(rd_mol, coords, is_3d=True)

    bridged = cosmolkit.Molecule.from_rdkit(rd_mol)

    assert bridged.num_conformers() == 1
    assert np.allclose(bridged.coordinates_3d(), coords)


def test_from_rdkit_copies_multiple_3d_conformers_and_skips_2d():
    rd_mol = Chem.MolFromSmiles("CO")
    coordinates_2d = np.array([[10.0, 11.0, 0.0], [12.0, 13.0, 0.0]])
    coordinates_3d_a = np.array([[0.0, 0.1, 0.2], [1.0, 1.1, 1.2]])
    coordinates_3d_b = np.array([[2.0, 2.1, 2.2], [3.0, 3.1, 3.2]])
    _add_conformer(rd_mol, coordinates_2d, is_3d=False)
    _add_conformer(rd_mol, coordinates_3d_a, is_3d=True)
    _add_conformer(rd_mol, coordinates_3d_b, is_3d=True)

    bridged = cosmolkit.Molecule.from_rdkit(rd_mol)

    assert bridged.num_conformers() == 2
    assert np.allclose(bridged.coordinates_3d(0), coordinates_3d_a)
    assert np.allclose(bridged.coordinates_3d(1), coordinates_3d_b)


def test_from_rdkit_does_not_copy_2d_conformer():
    rd_mol = Chem.MolFromSmiles("CO")
    _add_conformer(rd_mol, [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]], is_3d=False)

    bridged = cosmolkit.Molecule.from_rdkit(rd_mol)

    assert bridged.num_conformers() == 0
    with pytest.raises(ValueError, match="no 3D conformer"):
        bridged.coordinates_3d()


def test_from_rdkit_rejects_non_object():
    with pytest.raises(ValueError, match="from_rdkit failed calling GetNumAtoms"):
        cosmolkit.Molecule.from_rdkit(object())
