import pytest

import cosmolkit

Chem = pytest.importorskip("rdkit.Chem")
AllChem = pytest.importorskip("rdkit.Chem.AllChem")
MACCSkeys = pytest.importorskip("rdkit.Chem.MACCSkeys")
rdFingerprintGenerator = pytest.importorskip("rdkit.Chem.rdFingerprintGenerator")


SMOKE_SMILES = [
    "CCO",
    "c1ccccc1",
    "c1ccccc1O",
    "O=[N+]([O-])c1ccccc1",
    "N[C@@H](C)C(=O)O",
]


def _ck_bits(fp):
    return set(fp.on_bits())


def _rdkit_bits(fp):
    return set(fp.GetOnBits())


def _rdkit_additional_output_record(additional_output):
    return {
        "atom_counts": list(additional_output.GetAtomCounts()),
        "atom_to_bits": [list(bits) for bits in additional_output.GetAtomToBits()],
        "bit_info_map": {
            int(bit): [tuple(pair) for pair in pairs]
            for bit, pairs in additional_output.GetBitInfoMap().items()
        },
        "atoms_per_bit": {
            int(bit): [list(atoms) for atoms in atoms_per_bit]
            for bit, atoms_per_bit in additional_output.GetAtomsPerBit().items()
        },
    }


def _ck_additional_output_record(additional_output):
    return {
        "atom_counts": additional_output.atom_counts(),
        "atom_to_bits": additional_output.atom_to_bits(),
        "bit_info_map": additional_output.bit_info_map(),
        "atoms_per_bit": additional_output.atoms_per_bit(),
    }


@pytest.mark.parametrize("smiles", SMOKE_SMILES)
def test_morgan_fingerprint_is_rdkit_bit_identical(smiles):
    ck_mol = cosmolkit.Molecule.from_smiles(smiles)
    rd_mol = Chem.MolFromSmiles(smiles)
    generator = AllChem.GetMorganGenerator(radius=2, fpSize=2048)

    assert _ck_bits(ck_mol.fingerprint_morgan(radius=2, n_bits=2048)) == _rdkit_bits(
        generator.GetFingerprint(rd_mol)
    )


@pytest.mark.parametrize("smiles", SMOKE_SMILES)
def test_morgan_fingerprint_with_output_is_rdkit_bit_identical(smiles):
    ck_mol = cosmolkit.Molecule.from_smiles(smiles)
    rd_mol = Chem.MolFromSmiles(smiles)
    generator = AllChem.GetMorganGenerator(radius=2, fpSize=2048)
    additional_output = rdFingerprintGenerator.AdditionalOutput()
    additional_output.AllocateAtomCounts()
    additional_output.AllocateAtomToBits()
    additional_output.AllocateBitInfoMap()
    additional_output.AllocateAtomsPerBit()

    rd_fp = generator.GetFingerprint(rd_mol, additionalOutput=additional_output)
    ck_result = ck_mol.fingerprint_morgan_with_output(radius=2, n_bits=2048)

    assert _ck_bits(ck_result.fingerprint()) == _rdkit_bits(rd_fp)
    assert _ck_additional_output_record(
        ck_result.additional_output()
    ) == _rdkit_additional_output_record(additional_output)


@pytest.mark.parametrize("smiles", SMOKE_SMILES)
def test_maccs_fingerprint_is_rdkit_bit_identical(smiles):
    ck_mol = cosmolkit.Molecule.from_smiles(smiles)
    rd_mol = Chem.MolFromSmiles(smiles)
    rd_bits = {bit - 1 for bit in MACCSkeys.GenMACCSKeys(rd_mol).GetOnBits() if bit > 0}

    assert _ck_bits(ck_mol.maccs_fingerprint(n_bits=166)) == rd_bits


def test_maccs_fingerprint_rejects_non_rdkit_bit_length():
    ck_mol = cosmolkit.Molecule.from_smiles("NCCO")

    with pytest.raises(ValueError, match="MaccsFingerprintParams.n_bits"):
        ck_mol.maccs_fingerprint(n_bits=64)


def test_topological_fingerprint_matches_rdkit_exact_bits_and_is_deterministic():
    ck_mol = cosmolkit.Molecule.from_smiles("CCO")
    rd_mol = Chem.MolFromSmiles("CCO")

    ck_fp = ck_mol.topological_fingerprint(fp_size=64, num_bits_per_feature=1)
    rd_fp = Chem.RDKFingerprint(rd_mol, fpSize=64, nBitsPerHash=1)
    assert _ck_bits(ck_fp) == _rdkit_bits(rd_fp) == {0, 28, 59}
    assert _ck_bits(ck_mol.topological_fingerprint(fp_size=64, num_bits_per_feature=1)) == _ck_bits(ck_fp)


def test_topological_fingerprint_with_output_matches_rdkit_provenance():
    ck_mol = cosmolkit.Molecule.from_smiles("CCO")
    rd_mol = Chem.MolFromSmiles("CCO")
    atom_bits = []
    bit_info = {}
    rd_fp = Chem.RDKFingerprint(
        rd_mol,
        fpSize=64,
        nBitsPerHash=1,
        atomBits=atom_bits,
        bitInfo=bit_info,
    )
    result = ck_mol.topological_fingerprint_with_output(
        fp_size=64,
        num_bits_per_feature=1,
        atom_bits=True,
        bit_info=True,
    )
    assert _ck_bits(result.fingerprint()) == _rdkit_bits(rd_fp)
    assert result.atom_bits() == [[28, 0], [28, 59, 0], [59, 0]]
    assert result.bit_info() == {int(bit): [list(path) for path in paths] for bit, paths in bit_info.items()}


def test_topological_fingerprint_rejects_source_precondition_ranges():
    ck_mol = cosmolkit.Molecule.from_smiles("CCO")
    with pytest.raises(ValueError, match="minPath==0"):
        ck_mol.topological_fingerprint(min_path=0)


def test_avalon_fingerprint_returns_source_backed_bits_without_mutating_molecule():
    ck_mol = cosmolkit.Molecule.from_smiles("CCO")
    before = ck_mol.to_smiles()
    fp = ck_mol.avalon_fingerprint(n_bits=64, bit_flags=0x007FFF)

    assert fp.n_bits() == 64
    assert _ck_bits(fp) == {6, 14, 30, 31, 42}
    assert ck_mol.to_smiles() == before


def test_avalon_python_default_profile_matches_explicit_python_defaults():
    mol = cosmolkit.Molecule.from_smiles("CCO")

    default = mol.avalon_fingerprint(n_bits=64)
    explicit = mol.avalon_fingerprint(
        n_bits=64,
        is_query=False,
        bit_flags=0xF07FFF,
    )
    cpp_profile = mol.avalon_fingerprint(n_bits=64, bit_flags=0x007FFF)

    assert _ck_bits(default) == _ck_bits(explicit)
    assert _ck_bits(cpp_profile) == {6, 14, 30, 31, 42}
    assert _ck_bits(default) == {3, 6, 14, 30, 31, 42}


@pytest.mark.parametrize(
    "bit_flags",
    [0x000001, 0x000020, 0x000800, 0x004000, 0xF00000, 0xF07FFF],
)
def test_avalon_flag_profiles_are_typed_and_deterministic(bit_flags):
    mol = cosmolkit.Molecule.from_smiles("c1ccccc1O")

    first = mol.avalon_fingerprint(n_bits=128, bit_flags=bit_flags)
    second = mol.avalon_fingerprint(n_bits=128, bit_flags=bit_flags)

    assert first.n_bits() == 128
    assert _ck_bits(first) == _ck_bits(second)


def test_avalon_query_profile_uses_query_preprocessing_and_is_repeatable():
    mol = cosmolkit.Molecule.from_smiles("C[NH2+]C")

    first = mol.avalon_fingerprint(n_bits=64, is_query=True, bit_flags=0x007FFF)
    second = mol.avalon_fingerprint(n_bits=64, is_query=True, bit_flags=0x007FFF)

    assert first.n_bits() == 64
    assert _ck_bits(first) == set()
    assert _ck_bits(first) == _ck_bits(second)


@pytest.mark.parametrize("n_bits", [8, 9, 31, 32, 33, 511, 512, 513])
def test_avalon_size_boundaries_preserve_requested_public_length(n_bits):
    mol = cosmolkit.Molecule.from_smiles("CCO")

    fingerprint = mol.avalon_fingerprint(n_bits=n_bits, bit_flags=0x007FFF)

    assert fingerprint.n_bits() == n_bits


@pytest.mark.parametrize("n_bits", [0, 1, 7])
def test_avalon_rejects_sub_byte_sizes(n_bits):
    mol = cosmolkit.Molecule.from_smiles("CCO")

    with pytest.raises(ValueError, match="Avalon n_bits must be at least 8"):
        mol.avalon_fingerprint(n_bits=n_bits)


def test_avalon_rejects_unknown_flag_bits():
    mol = cosmolkit.Molecule.from_smiles("CCO")

    with pytest.raises(ValueError, match="Avalon bit_flags contains undefined source bits"):
        mol.avalon_fingerprint(n_bits=64, bit_flags=0x80000000)
