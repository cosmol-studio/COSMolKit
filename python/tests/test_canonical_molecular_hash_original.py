import pickle

import cosmolkit
import pytest


def test_original_mol_hash_benzene():
    assert cosmolkit.Molecule.from_smiles("c1ccccc1").molecular_hash() != 0


def test_original_mol_hash_deterministic():
    molecule = cosmolkit.Molecule.from_smiles("c1ccccc1")
    before = molecule.to_binary()
    assert molecule.molecular_hash() == molecule.molecular_hash()
    assert molecule.to_binary() == before


def test_original_mol_hash_different_molecules():
    assert cosmolkit.Molecule.from_smiles("c1ccccc1").molecular_hash() != cosmolkit.Molecule.from_smiles("CCO").molecular_hash()


def test_original_mol_hash_empty_error():
    with pytest.raises(cosmolkit.MoleculeHashError) as caught:
        cosmolkit.Molecule.from_smiles("").molecular_hash()
    assert caught.value.domain == "molecular_hash"
    assert caught.value.kind == "EmptyMolecule"


def test_original_mol_hash_with_ranks():
    # Benzene symmetry assigns the same source CIP rank to every atom.
    molecule = cosmolkit.Molecule.from_smiles("c1ccccc1")
    assert molecule.molecular_hash_with_ranks([0] * 6) != 0
    assert molecule.molecular_hash_with_ranks([0] * 6) == molecule.molecular_hash()


def test_full_original_chembl_hash_morgan_2d_and_binary_pickle_condition():
    molecule = cosmolkit.Molecule.from_smiles("CNC(=O)[C@H](CCCNC(=O)OC(C)(C)C)NC(=O)[C@H](CCCc1ccccc1)[C@@](C)(O)C(=O)NO").with_2d_coordinates()
    before = molecule.to_binary()
    expected_hash = molecule.molecular_hash()
    expected_morgan = molecule.fingerprint_morgan()
    restored = pickle.loads(pickle.dumps(molecule))
    assert restored.to_binary() == before
    assert restored.coordinates_2d() == molecule.coordinates_2d()
    assert restored.molecular_hash() == expected_hash
    actual_morgan = restored.fingerprint_morgan()
    assert actual_morgan.n_bits() == expected_morgan.n_bits()
    assert actual_morgan.on_bits() == expected_morgan.on_bits()
    assert molecule.to_binary() == before


def test_rank_row_and_absent_valence_errors_keep_carriers_and_source():
    molecule = cosmolkit.Molecule.from_smiles("CCO")
    before = molecule.to_binary()
    for ranks in ([], [0, 0], [0, 0, 0, 0]):
        with pytest.raises(cosmolkit.MoleculeHashError) as caught:
            molecule.molecular_hash_with_ranks(ranks)
        assert caught.value.kind == "RankCount"
        assert caught.value.actual == len(ranks)
        assert caught.value.atom_count == 3
        assert molecule.to_binary() == before
    builder = cosmolkit.MoleculeBuilder.new()
    builder.add_atom(cosmolkit.AtomSpec(cosmolkit.Element.C))
    raw = builder.build()
    before = raw.to_binary()
    with pytest.raises(cosmolkit.MoleculeHashError) as caught:
        raw.molecular_hash()
    assert caught.value.kind == "MissingPreparedValence"
    assert raw.molecular_hash_with_ranks([0]) != 0
    assert raw.to_binary() == before


def test_original_cip_rank_map_limit_is_typed_and_failure_atomic():
    builder = cosmolkit.MoleculeBuilder.new()
    builder.add_atom(cosmolkit.AtomSpec(cosmolkit.Element.C).with_atom_map(4294967295))
    molecule = builder.build().sanitize()
    before = molecule.to_binary()
    with pytest.raises(cosmolkit.MoleculeHashError) as caught:
        molecule.molecular_hash()
    assert caught.value.kind == "CipRanks"
    cause = caught.value.__cause__
    assert isinstance(cause, cosmolkit.CipRankError)
    assert cause.domain == "cip_ranking"
    assert cause.kind == "AtomMapOutOfRange"
    assert cause.atom == 0
    assert cause.map_number == 4294967295
    assert molecule.to_binary() == before
