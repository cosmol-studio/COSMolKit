from typing import assert_type
import numpy as np
import fixed_projection as ck

mol = ck.mol_from_smiles("CCO")
assert_type(mol, ck.Molecule)
assert_type(mol.num_atoms(), int)
params = ck.MaccsFingerprintParams(n_bits=166)
assert_type(params.n_bits, int)
params.n_bits = 166
assert_type(mol.fingerprint_maccs(), ck.Fingerprint)
assert_type(mol.fingerprint_maccs(params), ck.Fingerprint)
assert_type(mol.fingerprint_maccs(n_bits=166), ck.Fingerprint)
fp = ck.Fingerprint.from_on_bits(8, [1, 4])
batch = ck.FingerprintBatch([fp, None])
assert_type(batch[0], ck.Fingerprint | None)
assert_type(batch[:1], ck.FingerprintBatch)
assert_type(batch[np.int64(0)], ck.Fingerprint | None)
params.n_bits = np.int64(166)
assert_type(fp.to_numpy(), np.ndarray[tuple[int], np.dtype[np.uint8]])
assert_type(batch.to_numpy(), np.ndarray[tuple[int, int], np.dtype[np.uint8]])

try:
    fp.get_bit(9)
except ck.FingerprintError as error:
    assert_type(error.kind, ck.FingerprintErrorKind)
    assert_type(error.index, int | None)
