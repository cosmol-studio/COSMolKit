"""Exercise every frozen corpus row through the four canonical result forms.

These projection consistency checks are delivery proposals. They do not replace
independent source-native chemistry parity or the original fixed regressions.
"""

import hashlib
from pathlib import Path

import pytest

import cosmolkit as ck


CORPUS = Path(__file__).resolve().parents[2] / "testdata/smiles/corpus/smiles_5000.smi"
RAW = CORPUS.read_bytes()
assert hashlib.sha256(RAW).hexdigest() == "a4d579cd72621af27772256bb23ba796452276bb924fd20aac83625ffa67d849"
SMILES = [line.strip() for line in RAW.decode().splitlines() if line.strip() and not line.lstrip().startswith("#")]
assert len(SMILES) == 5000


@pytest.mark.parametrize("smiles", SMILES, ids=[str(i) for i in range(5000)])
def test_current_canonical_morgan_four_forms_and_populated_output(smiles):
    molecule = ck.Molecule.from_smiles(smiles)
    atom_count = molecule.num_atoms()
    params = ck.MorganFingerprintParams()
    output = ck.AdditionalOutput()
    output.allocate_atom_counts()
    output.allocate_atom_to_bits()
    output.allocate_bit_info_map()
    counts = molecule.morgan_sparse_count_fingerprint_with_params(params, output).nonzero_elements()
    presence = molecule.morgan_sparse_fingerprint().on_bits()
    folded = molecule.morgan_count_fingerprint()
    dense = molecule.morgan_fingerprint()

    assert dense.n_bits() == 2048
    assert folded.length() == 2048
    assert {value & 0xFFFFFFFF for value in presence} == set(counts)
    assert set(dense.on_bits()) == {value % 2048 for value in counts}
    assert sum(folded.nonzero_elements().values()) == sum(counts.values())
    assert len(output.atom_counts()) == atom_count
    assert len(output.atom_to_bits()) == atom_count
    assert set(output.bit_info_map()) == set(counts)
    assert output.bit_paths() is None
    assert output.atoms_per_bit() is None
    assert molecule.num_atoms() == atom_count
