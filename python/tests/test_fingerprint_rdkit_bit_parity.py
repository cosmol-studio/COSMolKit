import pytest

import cosmolkit

def _ck_bits(fp):
    return set(fp.on_bits())


def _ck_additional_output_record(additional_output):
    return {
        "atom_counts": additional_output.atom_counts(),
        "atom_to_bits": additional_output.atom_to_bits(),
        "bit_info_map": additional_output.bit_info_map(),
        "atoms_per_bit": additional_output.atoms_per_bit(),
    }


def test_morgan_additional_output_python_container_projection():
    # Fixed RDKit 2026.03.6 result: Morgan(radius=2, fpSize=2048), CCO.
    # Full chemistry comparisons belong to parity-tests_fixed. Keep one
    # representative list/tuple/dict projection without executing an oracle.
    ck_mol = cosmolkit.Molecule.from_smiles("CCO")
    ck_generator = cosmolkit.MorganFingerprintGenerator(
        params=cosmolkit.MorganParams(radius=2, fp_size=2048)
    )
    ck_output = cosmolkit.FingerprintAdditionalOutput()
    ck_output.allocate_atom_counts()
    ck_output.allocate_atom_to_bits()
    ck_output.allocate_bit_info_map()
    ck_output.allocate_atoms_per_bit()
    ck_fp = ck_mol.fingerprint_morgan_with_generator(ck_generator, output=ck_output)

    assert _ck_bits(ck_fp) == {80, 222, 294, 807, 1057, 1410}
    assert _ck_additional_output_record(ck_output) == {
        "atom_counts": [2, 2, 2],
        "atom_to_bits": [[1057, 294], [80, 1410], [807, 222]],
        "bit_info_map": {
            80: [(1, 0)], 222: [(2, 1)], 294: [(0, 1)],
            807: [(2, 0)], 1057: [(0, 0)], 1410: [(1, 1)],
        },
        "atoms_per_bit": {
            80: [[1]], 222: [[2, 1]], 294: [[0, 1]],
            807: [[2]], 1057: [[0]], 1410: [[1, 0, 2]],
        },
    }


def test_maccs_fingerprint_rejects_non_rdkit_bit_length():
    ck_mol = cosmolkit.Molecule.from_smiles("NCCO")

    with pytest.raises(ValueError, match="MaccsFingerprintParams.n_bits"):
        ck_mol.fingerprint_maccs(n_bits=64)


def test_topological_fingerprint_matches_rdkit_exact_bits_and_is_deterministic():
    ck_mol = cosmolkit.Molecule.from_smiles("CCO")

    # Fixed RDKit 2026.03.6 RDKFingerprint(fpSize=64, nBitsPerHash=1).
    ck_fp = ck_mol.fingerprint_topological(fp_size=64, num_bits_per_feature=1)
    assert _ck_bits(ck_fp) == {0, 28, 59}
    assert _ck_bits(ck_mol.fingerprint_topological(fp_size=64, num_bits_per_feature=1)) == _ck_bits(ck_fp)


def test_topological_fingerprint_with_output_matches_rdkit_provenance():
    ck_mol = cosmolkit.Molecule.from_smiles("CCO")
    result = ck_mol.fingerprint_topological_with_output(
        fp_size=64,
        num_bits_per_feature=1,
        atom_bits=True,
        bit_info=True,
    )
    assert _ck_bits(result.fingerprint()) == {0, 28, 59}
    assert result.atom_bits() == [[28, 0], [28, 59, 0], [59, 0]]
    assert result.bit_info() == {0: [[0, 1]], 28: [[0]], 59: [[1]]}


def test_topological_fingerprint_rejects_source_precondition_ranges():
    ck_mol = cosmolkit.Molecule.from_smiles("CCO")
    with pytest.raises(ValueError, match="minPath==0"):
        ck_mol.fingerprint_topological(min_path=0)


def test_avalon_fingerprint_returns_source_backed_bits_without_mutating_molecule():
    ck_mol = cosmolkit.Molecule.from_smiles("CCO")
    before = ck_mol.to_smiles()
    fp = ck_mol.fingerprint_avalon(n_bits=64, bit_flags=0x007FFF)

    assert fp.n_bits() == 64
    assert _ck_bits(fp) == {6, 14, 30, 31, 42}
    assert ck_mol.to_smiles() == before


def test_avalon_python_default_profile_matches_explicit_python_defaults():
    mol = cosmolkit.Molecule.from_smiles("CCO")

    default = mol.fingerprint_avalon(n_bits=64)
    explicit = mol.fingerprint_avalon(
        n_bits=64,
        is_query=False,
        bit_flags=0xF07FFF,
    )
    cpp_profile = mol.fingerprint_avalon(n_bits=64, bit_flags=0x007FFF)

    assert _ck_bits(default) == _ck_bits(explicit)
    assert _ck_bits(cpp_profile) == {6, 14, 30, 31, 42}
    assert _ck_bits(default) == {3, 6, 14, 30, 31, 42}


@pytest.mark.parametrize(
    "bit_flags",
    [0x000001, 0x000020, 0x000800, 0x004000, 0xF00000, 0xF07FFF],
)
def test_avalon_flag_profiles_are_typed_and_deterministic(bit_flags):
    mol = cosmolkit.Molecule.from_smiles("c1ccccc1O")

    first = mol.fingerprint_avalon(n_bits=128, bit_flags=bit_flags)
    second = mol.fingerprint_avalon(n_bits=128, bit_flags=bit_flags)

    assert first.n_bits() == 128
    assert _ck_bits(first) == _ck_bits(second)


def test_avalon_query_profile_uses_query_preprocessing_and_is_repeatable():
    mol = cosmolkit.Molecule.from_smiles("C[NH2+]C")

    first = mol.fingerprint_avalon(n_bits=64, is_query=True, bit_flags=0x007FFF)
    second = mol.fingerprint_avalon(n_bits=64, is_query=True, bit_flags=0x007FFF)

    assert first.n_bits() == 64
    assert _ck_bits(first) == set()
    assert _ck_bits(first) == _ck_bits(second)


@pytest.mark.parametrize("n_bits", [8, 9, 31, 32, 33, 511, 512, 513])
def test_avalon_size_boundaries_preserve_requested_public_length(n_bits):
    mol = cosmolkit.Molecule.from_smiles("CCO")

    fingerprint = mol.fingerprint_avalon(n_bits=n_bits, bit_flags=0x007FFF)

    assert fingerprint.n_bits() == n_bits


@pytest.mark.parametrize("n_bits", [0, 1, 7])
def test_avalon_rejects_sub_byte_sizes(n_bits):
    mol = cosmolkit.Molecule.from_smiles("CCO")

    with pytest.raises(ValueError, match="Avalon n_bits must be at least 8"):
        mol.fingerprint_avalon(n_bits=n_bits)


def test_avalon_rejects_unknown_flag_bits():
    mol = cosmolkit.Molecule.from_smiles("CCO")

    with pytest.raises(ValueError, match="Avalon bit_flags contains undefined source bits"):
        mol.fingerprint_avalon(n_bits=64, bit_flags=0x80000000)
