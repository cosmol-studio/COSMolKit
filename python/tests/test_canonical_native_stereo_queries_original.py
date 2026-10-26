"""Preserved project-native stereo query conditions and original tuple protocols."""
import pickle
import pytest
import cosmolkit as ck


@pytest.mark.parametrize('text, expected', [
    ('F[C@H](Cl)Br', [(1, [0, 2, 3, None])]),
    ('F[C@@H](Cl)Br', [(1, [0, 3, 2, None])]),
    ('F[C@](Cl)(Br)I', [(1, [0, 2, 3, 4])]),
    ('F[C@@](Cl)(Br)I', [(1, [0, 2, 4, 3])]),
    ('[13CH3:7][C@H](F)Cl', [(1, [0, 2, 3, None])]),
])
def test_original_ordered_ligand_values_and_typed_record(text, expected):
    source = ck.Molecule.from_smiles(text)
    before = source.to_smiles()
    result = source.tetrahedral_stereo()
    assert [(row.center, [ligand.atom_id for ligand in row.ligands]) for row in result] == expected
    row = result[0]
    assert isinstance(row, ck.TetrahedralStereo)
    center, ligands = row.center, list(row.ligands)
    assert center == expected[0][0] and len(ligands) == 4
    assert source.perceive_stereochemistry() is None
    assert source.to_smiles() == before
    ligands.append(99)
    assert [(row.center, [ligand.atom_id for ligand in row.ligands]) for row in source.tetrahedral_stereo()] == expected
    with pytest.raises(AttributeError): row.center = 0


def test_modern_chiral_center_filter_defaults_and_exact_cip_labels():
    source = ck.Molecule.from_smiles('F[C@H](Cl)Br')
    labels = source.find_chiral_centers()
    # Pinned RDKit 2026.03.1 modern FindMolChiralCenters.
    assert labels == [(1, 'R')]
    assert source.find_chiral_centers(include_unassigned=False) == labels
    assert ck.Molecule.from_smiles('CCO').find_chiral_centers(False) == []
    assert ck.Molecule.from_smiles('CCO').find_chiral_centers() == []
    potential = ck.Molecule.from_smiles('CC(O)Cl')
    assert potential.find_chiral_centers() == []
    assert potential.find_chiral_centers(include_unassigned=True) == [(1, '?')]


def test_original_empty_read_and_absent_valence_fallback_are_source_defined():
    assert ck.Molecule.new().tetrahedral_stereo() == []
    assert ck.Molecule.new().perceive_stereochemistry() is None
    assert ck.Molecule.new().find_chiral_centers() == []
    raw = ck.Molecule.from_smiles_with_params('F[C@](Cl)Br', ck.SmilesParseParams(sanitize=False))
    assert raw.tetrahedral_stereo() == []
    assert raw.perceive_stereochemistry() is None


def test_canonical_ligand_values_are_immutable_and_isolate_identifier_state():
    atom = ck.LigandRef(3)
    hydrogen = ck.LigandRef()
    assert atom.atom_id == 3 and atom != ck.LigandRef.ImplicitHydrogen
    assert hydrogen.atom_id is None and hydrogen == ck.LigandRef.ImplicitHydrogen
    with pytest.raises(AttributeError): atom.atom_id = 1
    with pytest.raises(OverflowError): ck.LigandRef(-1)
