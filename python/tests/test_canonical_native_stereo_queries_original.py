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
def test_original_ordered_ligand_values_and_tuple_protocol(text, expected):
    source = ck.Molecule.from_smiles(text)
    before = source.to_smiles()
    result = source.tetrahedral_stereo()
    assert result == expected
    row = result[0]
    assert isinstance(row, ck.TetrahedralStereo)
    assert isinstance(row, tuple)
    center, ligands = row
    assert len(row) == 2
    assert row.center == row[0] == center
    assert row.ligands == row[1] == ligands
    assert pickle.loads(pickle.dumps(row)) == row
    assert source.perceive_stereochemistry() is None
    assert source.to_smiles() == before
    ligands.append(99)
    assert source.tetrahedral_stereo() == expected
    with pytest.raises(AttributeError): row.center = 0


def test_original_chiral_center_filter_defaults_and_exact_tag_labels():
    source = ck.Molecule.from_smiles('F[C@H](Cl)Br')
    labels = source.find_chiral_centers()
    assert labels == [(0, '?'), (1, 'CHI_TETRAHEDRAL_CCW'), (2, '?'), (3, '?')]
    assert source.find_chiral_centers(include_unassigned=False) == [labels[1]]
    assert ck.Molecule.from_smiles('CCO').find_chiral_centers(False) == []
    assert ck.Molecule.from_smiles('CCO').find_chiral_centers() == [(0, '?'), (1, '?'), (2, '?')]


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
    assert atom.atom == 3 and atom.is_implicit_hydrogen is False
    assert hydrogen.atom is None and hydrogen.is_implicit_hydrogen is True
    with pytest.raises(AttributeError): atom.atom = 1
    with pytest.raises(OverflowError): ck.LigandRef(-1)
