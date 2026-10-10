"""The registry must enforce the NumPy collection projection, not just docs."""

import ast
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import cosmolkit as ck


_path = Path(__file__).resolve().parents[2] / "dev/tools/check_python_stub_contract.py"
_spec = importlib.util.spec_from_file_location("fingerprint_numpy_checker", _path)
_checker = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_checker)

DOCUMENT = {
    "entries": [
        {"item": "type", "semantic_id": "types.MoleculeBatch", "python_name": "MoleculeBatch"},
        {"item": "callable", "semantic_id": "MoleculeBatch.fingerprint_morgan_list",
         "python_name": "fingerprint_morgan_list", "output": "Vec<Option<crate::Fingerprint>>"},
    ],
    "python_collections": [{"name": "FingerprintBatch", "rust_output": "Vec<Option<crate::Fingerprint>>"}],
}

STUB = """
class FingerprintBatch:
    def __new__(cls, fingerprints): ...
    def __len__(self): ...
    def __iter__(self): ...
    def __repr__(self): ...
    @overload
    def __getitem__(self, key: int) -> Fingerprint | None: ...
    @overload
    def __getitem__(self, key: slice) -> FingerprintBatch: ...
    def to_numpy(self) -> numpy.typing.NDArray[numpy.uint8]: ...
class MoleculeBatch:
    def fingerprint_morgan_list(self) -> FingerprintBatch: ...
"""


def test_numpy_collection_gate_rejects_missing_class_list_return_and_wrong_dtype():
    check = _checker.check_collection_declarations
    assert check(ast.parse(STUB), DOCUMENT) == []
    missing = "class MoleculeBatch:\n    def fingerprint_morgan_list(self) -> list: ..."
    assert any("collection missing" in error for error in check(ast.parse(missing), DOCUMENT))
    wrong_return = STUB.replace("fingerprint_morgan_list(self) -> FingerprintBatch", "fingerprint_morgan_list(self) -> list")
    assert any("expected registered FingerprintBatch" in error for error in check(ast.parse(wrong_return), DOCUMENT))
    wrong_dtype = STUB.replace("numpy.uint8", "numpy.float64")
    assert any("uint8" in error for error in check(ast.parse(wrong_dtype), DOCUMENT))


def test_numpy_collection_gate_requires_a_real_runtime_implementation():
    check = _checker.check_collection_runtime
    assert check(SimpleNamespace(), DOCUMENT)
    assert check(ck, DOCUMENT) == []
