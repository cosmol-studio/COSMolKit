import cosmolkit
import pytest


def test_atom_pair_default_keywords_convert_all_four_exact_result_forms():
    molecule = cosmolkit.Molecule.from_smiles("CCCO")

    explicit = molecule.fingerprint_atom_pair()
    assert isinstance(explicit, cosmolkit.Fingerprint)
    assert explicit.n_bits() == 2048
    assert explicit.on_bits() == [624, 1144, 1336, 1337, 1404, 1596]

    sparse_count = molecule.fingerprint_atom_pair_sparse_count()
    assert isinstance(sparse_count, cosmolkit.AtomPairSparseCountFingerprint)
    assert sparse_count.size() == 1 << 23
    assert sparse_count.nonzero_elements() == {
        558113: 1,
        558114: 1,
        558145: 1,
        1590307: 1,
        1590337: 1,
        1590338: 1,
    }
    assert sparse_count.value(558113) == 1
    assert sparse_count.value(0) == 0

    count = molecule.fingerprint_atom_pair_count()
    assert isinstance(count, cosmolkit.AtomPairSparseCountFingerprint)
    assert count.size() == 2048
    assert count.nonzero_elements() == {
        1310: 1,
        1358: 2,
        1375: 1,
        1423: 1,
        1692: 1,
    }

    sparse_bits = molecule.fingerprint_atom_pair_sparse_bits()
    assert isinstance(sparse_bits, cosmolkit.AtomPairSparseBitFingerprint)
    assert sparse_bits.size() == 1 << 23
    assert sparse_bits.on_bits() == [
        7918712,
        7918972,
        7920240,
        8066360,
        8066361,
        8066620,
    ]


def test_atom_pair_provenance_matches_exact_count_simulation_projection():
    result = cosmolkit.Molecule.from_smiles("CCCO").fingerprint_atom_pair_with_output()
    assert isinstance(result, cosmolkit.AtomPairFingerprintResult)
    assert result.fingerprint().on_bits() == [624, 1144, 1336, 1337, 1404, 1596]

    output = result.additional_output()
    assert isinstance(output, cosmolkit.AtomPairAdditionalOutput)
    assert output.atom_counts() == [3, 3, 3, 3]
    assert output.atom_to_bits() == [
        [624, 1144, 1404],
        [1336, 1337, 1404],
        [1144, 1336, 1337, 1596],
        [624, 1336, 1337, 1596],
    ]
    expected_pairs = {
        624: [(0, 3)],
        1144: [(0, 2)],
        1336: [(1, 2), (1, 3)],
        1337: [(1, 2), (1, 3)],
        1404: [(0, 1)],
        1596: [(2, 3)],
    }
    assert output.bit_info_map() == expected_pairs
    assert output.atoms_per_bit() == {
        bit: [list(pair) for pair in pairs] for bit, pairs in expected_pairs.items()
    }


def test_atom_pair_option_interactions_and_repeated_calls_are_exact_and_immutable():
    molecule = cosmolkit.Molecule.from_smiles("C[C@H](O)F")
    before = molecule.to_smiles()
    kwargs = {
        "n_bits": 64,
        "min_distance": 1,
        "max_distance": 2,
        "include_chirality": True,
        "count_simulation": False,
        "from_atoms": [1],
        "ignore_atoms": [3],
        "custom_atom_invariants": [11, 22, 33, 44],
    }
    first = molecule.fingerprint_atom_pair(**kwargs)
    second = molecule.fingerprint_atom_pair(**kwargs)
    assert first.on_bits() == second.on_bits() == [0, 37]
    assert molecule.fingerprint_atom_pair_sparse_count(
        **kwargs
    ).nonzero_elements() == {1442145: 1, 2163393: 1}
    assert molecule.to_smiles() == before


def test_atom_pair_binding_returns_typed_value_errors():
    molecule = cosmolkit.Molecule.from_smiles("CCO")
    with pytest.raises(ValueError, match="n_bits > 0"):
        molecule.fingerprint_atom_pair(n_bits=0)
    with pytest.raises(ValueError, match="minDistance"):
        molecule.fingerprint_atom_pair(min_distance=4, max_distance=3)
    with pytest.raises(ValueError, match="(?i)conformer"):
        molecule.fingerprint_atom_pair(use_2d=False)


def test_atom_pair_python_surface_has_no_rdkit_style_duplicate_names():
    forbidden = {
        "GetAtomPairGenerator",
        "GetAtomPairFingerprint",
        "GetHashedAtomPairFingerprint",
        "get_atom_pair_generator",
        "get_atom_pair_fingerprint",
        "get_hashed_atom_pair_fingerprint",
    }
    assert forbidden.isdisjoint(dir(cosmolkit))
    assert forbidden.isdisjoint(dir(cosmolkit.Molecule))
