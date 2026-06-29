import pytest

import cosmolkit

Chem = pytest.importorskip("rdkit.Chem")


CASES = [
    ("CCO", "CO"),
    ("CCCO", "CC"),
    ("c1ccccc1O", "c1ccccc1"),
    ("O=[N+]([O-])c1ccccc1", "N(=O)[O-]"),
    ("N[C@@H](C)C(=O)O", "C(=O)O"),
]


def _rdkit_matches(smiles, query):
    mol = Chem.MolFromSmiles(smiles)
    qmol = Chem.MolFromSmiles(query)
    assert mol is not None
    assert qmol is not None
    return sorted(tuple(match) for match in mol.GetSubstructMatches(qmol, uniquify=True))


def _ck_matches(smiles, query):
    mol = cosmolkit.Molecule.from_smiles(smiles)
    qmol = cosmolkit.Molecule.from_smiles(query)
    return sorted(
        tuple(match.atom_mapping()) for match in cosmolkit.get_substruct_matches(mol, qmol)
    )


@pytest.mark.parametrize(("smiles", "query"), CASES)
def test_molecule_query_substructure_matches_rdkit_exactly(smiles, query):
    rdkit_matches = _rdkit_matches(smiles, query)
    ck_matches = _ck_matches(smiles, query)

    assert bool(ck_matches) == bool(rdkit_matches)
    assert ck_matches == rdkit_matches
