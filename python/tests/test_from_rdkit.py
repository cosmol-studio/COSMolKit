from pathlib import Path

import pytest

import cosmolkit

Chem = pytest.importorskip("rdkit.Chem")


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


def _normalized_rdkit_bond_stereo(stereo):
    name = str(stereo)
    if name == "STEREOE":
        return "STEREOTRANS"
    if name == "STEREOZ":
        return "STEREOCIS"
    return name


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
            str(atom.GetChiralTag()),
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
            "AROMATIC" if bond.GetIsAromatic() else str(bond.GetBondType()),
            str(bond.GetBondDir()),
            _normalized_rdkit_bond_stereo(bond.GetStereo()),
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

    assert bridged.to_smiles(True) == direct.to_smiles(True)
    assert bridged.to_smiles(False) == direct.to_smiles(False)


def test_from_rdkit_rejects_non_object():
    with pytest.raises(ValueError, match="from_rdkit failed calling GetNumAtoms"):
        cosmolkit.Molecule.from_rdkit(object())
