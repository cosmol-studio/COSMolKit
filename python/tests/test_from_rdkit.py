from pathlib import Path

import numpy as np
import pytest

import cosmolkit

Chem = pytest.importorskip("rdkit.Chem")
Point3D = pytest.importorskip("rdkit.Geometry").Point3D


def _load_smiles_cases():
    corpus = (
        Path(__file__).resolve().parents[2]
        / "testdata"
        / "smiles"
        / "corpus"
        / "smiles_small.smi"
    )
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
            atom.id(),
            atom.atomic_number(),
            atom.formal_charge(),
            atom.chiral_tag(),
            atom.isotope(),
            atom.atom_map(),
        )
        for atom in mol.atoms()
    ]
    bonds = [
        (
            bond.id(),
            min(bond.begin(), bond.end()),
            max(bond.begin(), bond.end()),
            bond.order(),
        )
        for bond in mol.bonds()
    ]
    return atoms, sorted(bonds)


def _feature_signature(mol):
    atoms = [
        (
            atom.id(),
            atom.atomic_number(),
            atom.formal_charge(),
            atom.chiral_tag(),
            atom.isotope(),
            atom.atom_map(),
            atom.is_aromatic(),
            atom.explicit_hydrogens(),
            atom.no_implicit(),
            atom.radical_electrons(),
            atom.hybridization(),
            atom.degree(),
            atom.explicit_valence(),
            atom.implicit_hydrogens(),
            atom.total_hydrogens(),
            atom.total_valence(),
        )
        for atom in mol.atoms()
    ]
    bonds = [
        (
            bond.id(),
            min(bond.begin(), bond.end()),
            max(bond.begin(), bond.end()),
            bond.order(),
            bond.direction(),
            bond.stereo(),
            tuple(bond.stereo_atoms() or ()),
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
            cosmolkit.ChiralTag(int(atom.GetChiralTag())),
            atom.GetIsotope() or None,
            atom.GetAtomMapNum() or None,
            atom.GetIsAromatic(),
            atom.GetNumExplicitHs(),
            atom.GetNoImplicit(),
            atom.GetNumRadicalElectrons(),
            cosmolkit.Hybridization(int(atom.GetHybridization())),
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
            cosmolkit.BondOrder(int(bond.GetBondType())),
            cosmolkit.BondDirection(int(bond.GetBondDir())),
            cosmolkit.BondStereo(int(bond.GetStereo())),
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
    assert bond.order() == cosmolkit.BondOrder.DOUBLE
    assert bond.order_code() == int(cosmolkit.BondOrder.DOUBLE)
    assert bond.direction() == cosmolkit.BondDirection.NONE
    assert bond.stereo() == cosmolkit.BondStereo.STEREONONE
    assert cosmolkit.BondOrder.DOUBLE == int(Chem.BondType.DOUBLE)
    assert not hasattr(bond, "bond_type")


def test_public_bond_enums_include_hydrogen_and_unknown_members():
    # Bond.h enum order and the pinned native enum both give HYDROGEN=14,
    # not the historical CK-only ordinal 18. Compare the source directly.
    assert cosmolkit.BondOrder.HYDROGEN == int(Chem.BondType.HYDROGEN) == 14
    assert cosmolkit.BondDirection.UNKNOWN == int(Chem.BondDir.UNKNOWN) == 6


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

    assert bridged.num_3d_conformers() == 1
    assert np.allclose(bridged.coordinates_3d(0), coords)
    assert np.array_equal(np.asarray(bridged.coordinates_3d(0)).view(np.uint64), coords.view(np.uint64))


def test_from_rdkit_defaults_to_prepared_graph_for_3d_atom_pair_fingerprint():
    rd_mol = Chem.MolFromSmiles("C=C")
    _add_conformer(
        rd_mol,
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        is_3d=True,
    )

    default = cosmolkit.Molecule.from_rdkit(rd_mol)
    explicit = cosmolkit.Molecule.from_rdkit(rd_mol, sanitize=True)

    assert default.fingerprint_atom_pair_with_params(cosmolkit.AtomPairFingerprintParams(generator=cosmolkit.AtomPairParams(use_2d=False)), None).on_bits() == [1432]
    assert (
        default.fingerprint_atom_pair_with_params(cosmolkit.AtomPairFingerprintParams(generator=cosmolkit.AtomPairParams(use_2d=False)), None).on_bits()
        == explicit.fingerprint_atom_pair_with_params(cosmolkit.AtomPairFingerprintParams(generator=cosmolkit.AtomPairParams(use_2d=False)), None).on_bits()
    )


def test_from_rdkit_sanitize_false_preserves_unprepared_graph_state():
    rd_mol = Chem.MolFromSmiles("C=C")
    _add_conformer(
        rd_mol,
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        is_3d=True,
    )

    raw = cosmolkit.Molecule.from_rdkit(rd_mol, sanitize=False)

    # Retain the unprepared-cache assertion independently of the current
    # fingerprint facade's earlier prepared-assignment precondition.
    with pytest.raises(cosmolkit.ValenceError) as error:
        raw.atom_metadata(recalculate=False)
    assert error.value.kind == "ExplicitValenceCacheNotInitialized"
    assert error.value.atom == 0
    with pytest.raises(cosmolkit.AtomPairReadError, match="Fingerprint preparation requires a valid prepared valence assignment"):
        raw.fingerprint_atom_pair_with_params(cosmolkit.AtomPairFingerprintParams(generator=cosmolkit.AtomPairParams(use_2d=False)), None)


def test_from_rdkit_copies_multiple_3d_conformers_and_skips_2d():
    rd_mol = Chem.MolFromSmiles("CO")
    coordinates_2d = np.array([[10.0, 11.0, 0.0], [12.0, 13.0, 0.0]])
    coordinates_3d_a = np.array([[0.0, 0.1, 0.2], [1.0, 1.1, 1.2]])
    coordinates_3d_b = np.array([[2.0, 2.1, 2.2], [3.0, 3.1, 3.2]])
    _add_conformer(rd_mol, coordinates_2d, is_3d=False)
    _add_conformer(rd_mol, coordinates_3d_a, is_3d=True)
    _add_conformer(rd_mol, coordinates_3d_b, is_3d=True)

    bridged = cosmolkit.Molecule.from_rdkit(rd_mol)

    assert bridged.num_3d_conformers() == 2
    assert np.allclose(bridged.coordinates_3d(0), coordinates_3d_a)
    assert np.allclose(bridged.coordinates_3d(1), coordinates_3d_b)


def test_from_rdkit_does_not_copy_2d_conformer():
    rd_mol = Chem.MolFromSmiles("CO")
    _add_conformer(rd_mol, [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]], is_3d=False)

    bridged = cosmolkit.Molecule.from_rdkit(rd_mol)

    assert bridged.num_3d_conformers() == 0
    with pytest.raises(ValueError, match="no 3D conformer"):
        bridged.coordinates_3d(0)


def test_from_rdkit_rejects_non_object():
    with pytest.raises(ValueError, match="from_rdkit failed calling GetNumAtoms"):
        cosmolkit.Molecule.from_rdkit(object())


@pytest.mark.parametrize("sanitize", [None, False, True])
def test_from_rdkit_preserves_source_and_returns_independent_storage(sanitize):
    source = Chem.MolFromSmiles("[13CH3:7][C@H](F)Cl")
    positions = np.array([[-0., 0.1, -0.2], [1., 2., 3.], [4., 5., 6.], [7., 8., 9.]])
    _add_conformer(source, positions, is_3d=True)
    before = source.ToBinary()
    imported = cosmolkit.Molecule.from_rdkit(source, sanitize=sanitize)
    assert source.ToBinary() == before
    assert np.array_equal(np.asarray(imported.coordinates_3d(0)).view(np.uint64), positions.view(np.uint64))
    source.GetAtomWithIdx(0).SetIsotope(12)
    source.GetConformer().SetAtomPosition(0, Point3D(100., 200., 300.))
    assert imported.atoms()[0].isotope() == 13
    assert np.array_equal(np.asarray(imported.coordinates_3d(0)).view(np.uint64), positions.view(np.uint64))


def test_from_rdkit_default_prepares_valence_without_sanitizing():
    source = Chem.MolFromSmiles("CC", sanitize=False)
    before = source.ToBinary()
    prepared = cosmolkit.Molecule.from_rdkit(source)
    sanitized = cosmolkit.Molecule.from_rdkit(source, sanitize=True)
    reference = Chem.Mol(source)
    reference.UpdatePropertyCache(strict=True)
    assert _feature_signature(prepared) == _rdkit_signature(reference)
    Chem.SanitizeMol(reference)
    assert _feature_signature(sanitized) == _rdkit_signature(reference)
    assert prepared.atoms()[0].hybridization() == cosmolkit.Hybridization.UNSPECIFIED
    assert sanitized.atoms()[0].hybridization() == cosmolkit.Hybridization.SP3
    assert source.ToBinary() == before


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -float("inf")])
def test_from_rdkit_rejects_invalid_coordinates_without_changing_source(value):
    source = Chem.MolFromSmiles("C")
    _add_conformer(source, [[value, 0., 0.]], is_3d=True)
    before = source.ToBinary()
    with pytest.raises(cosmolkit.OperationError):
        cosmolkit.Molecule.from_rdkit(source)
    assert source.ToBinary() == before
