"""Finite framework-category cases; no chemistry corpus or per-API stub inventory."""
import json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(sys.argv[1]).resolve()))
import fixed_projection as ck
from fixed_projection import _native

count = 0


def require(value):
    global count
    assert value
    count += 1


def raises(kind, call):
    global count
    try:
        call()
    except kind as error:
        count += 1
        return error
    raise AssertionError(f"expected {kind}")


mol = ck.mol_from_smiles("CCO")
require(mol.num_atoms() == 3)
default = mol.fingerprint_maccs().on_bits()
require(default == mol.fingerprint_maccs(ck.MaccsFingerprintParams()).on_bits())
require(default == mol.fingerprint_maccs(n_bits=166).on_bits())
require(default == mol._fingerprint_maccs_params().on_bits())
require(default == mol._fingerprint_maccs_keywords(n_bits=166).on_bits())
params = ck.MaccsFingerprintParams()
params.n_bits = 167
native_maccs = _native.Molecule.fingerprint_maccs
require(ck.Molecule.fingerprint_maccs is native_maccs)
require(ck.mol_from_smiles_configured is _native.mol_from_smiles_configured)
native_error = raises(ck.MaccsFingerprintError, lambda: native_maccs(mol, params))
for call in (lambda: mol.fingerprint_maccs(params), lambda: mol.fingerprint_maccs(n_bits=167)):
    error = raises(ck.MaccsFingerprintError, call)
    require(str(error) == str(native_error) and error.kind == native_error.kind)
    require(error.option == native_error.option and error.reason == native_error.reason)
raises(TypeError, lambda: mol.fingerprint_maccs(params, n_bits=166))
raises(TypeError, lambda: mol.fingerprint_maccs(unknown=True))
raises(TypeError, lambda: setattr(params, "n_bits", "wrong"))
require(params.n_bits == 167)
raises(AttributeError, lambda: setattr(params, "unknown", 1))

fp = ck.Fingerprint.from_on_bits(8, [1, 4])
array = fp.to_numpy()
require(array.dtype == np.uint8 and array.shape == (8,))
require(np.flatnonzero(array).tolist() == [1, 4])
array[:] = 0
require(fp.on_bits() == [1, 4])
batch = ck.FingerprintBatch([fp, None, fp])
require(batch[0] is fp and batch[1] is None and batch[-1] is fp)
class IntegerIndex:
    def __index__(self):
        return 0
require(batch[IntegerIndex()] is fp)
require(batch[np.int64(0)] is fp)
require(ck.MaccsFingerprintParams(n_bits=np.int64(166)).n_bits == 166)
params.n_bits = np.int64(166)
require(params.n_bits == 166)
require(isinstance(batch[:2], ck.FingerprintBatch) and len(batch[:2]) == 2)
require(batch[::-1][0] is fp)
raises(IndexError, lambda: batch[9])
raises(TypeError, lambda: batch["wrong"])
raises(ValueError, lambda: batch[::0])
raises(ValueError, batch.to_numpy)
matrix = ck.FingerprintBatch([fp, fp]).to_numpy()
require(matrix.shape == (2, 8) and matrix.dtype == np.uint8)
require(matrix.tolist() == [[0, 1, 0, 0, 1, 0, 0, 0]] * 2)
error = raises(ck.FingerprintError, lambda: fp.get_bit(9))
require(error.index == 9 and error.size == 8)
require(error.kind == ck.FingerprintErrorKind.SparseIndexOutOfRange)
raises(AttributeError, lambda: setattr(error, "index", 0))
error = raises(ck.SmilesError, lambda: ck.mol_from_smiles("C(C"))
require(error.source_message == str(error.__cause__))

p = ck.SmilesParseParams(remove_hs=False)
require(ck.mol_from_smiles_configured("[H]O[H]", p).num_atoms() == 3)
require(ck.mol_from_smiles_configured("[H]O[H]", remove_hs=False).num_atoms() == 3)
require(ck.mol_from_smiles_configured("[H]O[H]").num_atoms() == 1)
require(_native._mol_from_smiles_configured_params("[H]O[H]", p).num_atoms() == 3)
require(_native._mol_from_smiles_configured_keywords("[H]O[H]", remove_hs=False).num_atoms() == 3)
p.replacements = {"{A}": "O"}
require(ck.mol_from_smiles_configured("C{A}", p).num_atoms() == 2)
require(ck.mol_from_smiles_configured("C{A}", replacements={"{A}": "O"}).num_atoms() == 2)
raises(TypeError, lambda: ck.mol_from_smiles_configured("CCO", p, sanitize=False))
raises(TypeError, lambda: setattr(p, "replacements", {1: "O"}))
require(p.replacements == {"{A}": "O"})
mols = ck.MoleculeBatch([mol, None, mol])
require(mols[0] is mol and mols[1] is None)
require(isinstance(mols[::-1], ck.MoleculeBatch) and mols[::-1][0] is mol)
require(ck.SmilesParseParams().replacements == {})
# Two optional-configuration calls must not share mutable defaults.
p.replacements["unused"] = "N"
require(ck.SmilesParseParams().replacements == {})

print(json.dumps({"checks": count, "public_exports": len(ck.__all__)}))
