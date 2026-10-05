"""Complete unchanged original MACCS records through actual native Python projections."""
import hashlib
import json
import os
from pathlib import Path

import pytest
import cosmolkit

REFERENCE_PATH = Path(os.environ["COSMOLKIT_MACCS_ORIGINAL5000_GOLDEN"])
REFERENCE_BYTES = REFERENCE_PATH.read_bytes()
assert hashlib.sha256(REFERENCE_BYTES).hexdigest() == "5e94d3b504cbdfab3a0c6b65f8d37c696f444b58625d1b894a69e0358f675687"
RECORDS = [json.loads(line) for line in REFERENCE_BYTES.splitlines()]
assert len(RECORDS) == 5093
assert sum(row["record_type"] == "fixture" for row in RECORDS) == 93
assert sum(row["record_type"] == "corpus" for row in RECORDS) == 5000
assert all(row["rdkit_ok"] and row["error"] is None for row in RECORDS)

@pytest.fixture(scope="session")
def actual_observations():
    with Path(os.environ["COSMOLKIT_MACCS_OBSERVATIONS"]).open("x") as stream:
        yield stream

@pytest.mark.parametrize("index,reference", list(enumerate(RECORDS)))
def test_original_raw_public_and_option_error_native_record(index, reference, actual_observations):
    molecule = cosmolkit.Molecule.from_smiles(reference["smiles"])
    params = cosmolkit.MaccsFingerprintParams()
    raw = molecule.maccs_fingerprint_raw()
    public = molecule.maccs_fingerprint()
    configured = molecule.maccs_fingerprint_with_params(params)
    with pytest.raises(cosmolkit.MaccsFingerprintError) as failure:
        molecule.maccs_fingerprint_with_params(cosmolkit.MaccsFingerprintParams(n_bits=64))
    error = failure.value
    observation = {
        "index": index, "record_type": reference["record_type"], "label": reference["label"],
        "smiles": reference["smiles"], "raw_n_bits": raw.n_bits(), "raw_on_bits": raw.on_bits(),
        "public_n_bits": public.n_bits(), "public_on_bits": public.on_bits(),
        "configured_n_bits": configured.n_bits(), "configured_on_bits": configured.on_bits(),
        "params_n_bits": params.n_bits, "unsupported_width": 64,
        "error": {"type": type(error).__name__, "message": str(error),
                  "domain": error.domain, "kind": error.kind, "option": error.option, "reason": error.reason},
    }
    actual_observations.write(json.dumps(observation, sort_keys=True) + "\n")
    assert observation["raw_n_bits"] == reference["raw_n_bits"] == 167
    assert 0 not in observation["raw_on_bits"]
    assert observation["raw_on_bits"] == reference["raw_on_bits"]
    assert observation["public_n_bits"] == reference["public_n_bits"] == 166
    assert observation["public_on_bits"] == reference["public_on_bits"]
    assert observation["configured_n_bits"] == reference["public_n_bits"]
    assert observation["configured_on_bits"] == reference["public_on_bits"]
    assert params.n_bits == 166
    assert error.domain == "Fingerprint"
    assert error.kind == "UnsupportedOption"
    assert error.option == "MaccsFingerprintParams.n_bits"
    assert error.reason == "RDKit MACCS exposes a fixed 167-bit raw vector with bit 0 unused; COSMolKit only exposes the exact 166-bit public projection"
    assert str(error) == "unsupported fingerprint option " + error.option + ": " + error.reason
