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


@pytest.mark.xfail(reason="topological fingerprint is not source-ported to RDKit exact parity yet")
@pytest.mark.parametrize("smiles", SMOKE_SMILES)
def test_topological_fingerprint_is_rdkit_bit_identical(smiles):
    ck_mol = cosmolkit.Molecule.from_smiles(smiles)
    rd_mol = Chem.MolFromSmiles(smiles)
    rd_fp = Chem.RDKFingerprint(
        rd_mol,
        minPath=1,
        maxPath=7,
        fpSize=2048,
        nBitsPerHash=2,
        useHs=False,
        tgtDensity=0.0,
        minSize=2048,
    )

    assert _ck_bits(ck_mol.topological_fingerprint()) == _rdkit_bits(rd_fp)


@pytest.mark.xfail(reason="Avalon fingerprint is not source-ported to RDKit exact parity yet")
@pytest.mark.parametrize("smiles", SMOKE_SMILES)
def test_avalon_fingerprint_is_rdkit_bit_identical(smiles):
    py_avalon_tools = pytest.importorskip("rdkit.Avalon.pyAvalonTools")
    ck_mol = cosmolkit.Molecule.from_smiles(smiles)
    rd_mol = Chem.MolFromSmiles(smiles)

    assert _ck_bits(ck_mol.avalon_fingerprint(n_bits=2048)) == _rdkit_bits(
        py_avalon_tools.GetAvalonFP(rd_mol, nBits=2048)
    )
