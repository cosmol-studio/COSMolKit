"""Original Python5 MACCS smoke inputs/source expectations; canonical parameter transport."""
import pytest
import cosmolkit

ORIGINAL_SMOKE_CASES = [('CCO', [81, 108, 113, 138, 152, 154, 156, 159, 163]), ('c1ccccc1', [161, 162, 164]), ('c1ccccc1O', [112, 126, 138, 142, 151, 156, 161, 162, 163, 164]), ('O=[N+]([O-])c1ccccc1', [23, 48, 55, 62, 69, 70, 93, 101, 118, 121, 123, 129, 132, 134, 147, 155, 157, 158, 160, 161, 162, 163, 164]), ('N[C@@H](C)C(=O)O', [53, 83, 94, 122, 130, 138, 150, 153, 155, 156, 157, 158, 159, 160, 163])]

@pytest.mark.parametrize("smiles,expected", ORIGINAL_SMOKE_CASES)
def test_maccs_fingerprint_is_rdkit_bit_identical(smiles, expected):
    molecule = cosmolkit.Molecule.from_smiles(smiles)
    actual = molecule.fingerprint_maccs_with_params(cosmolkit.MaccsFingerprintParams(n_bits=166))
    assert set(actual.on_bits()) == set(expected)

def test_maccs_fingerprint_rejects_non_rdkit_bit_length():
    molecule = cosmolkit.Molecule.from_smiles("NCCO")
    with pytest.raises(ValueError, match="MaccsFingerprintParams.n_bits"):
        molecule.fingerprint_maccs_with_params(cosmolkit.MaccsFingerprintParams(n_bits=64))
