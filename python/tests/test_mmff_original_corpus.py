"""Initial MMFF energy/gradient references, not covered by parity's optimizer tasks."""
import hashlib
import json
import os
from pathlib import Path
import cosmolkit as ck

def reference_path(relative):
    receipt_root = os.environ.get("COSMOLKIT_MMFF_RECEIPT_ROOT")
    if receipt_root:
        return Path(receipt_root) / relative
    # Ordinary repository runs use the same original fixtures directly;
    # immutable receipt runs can point at their frozen copies explicitly.
    return Path(__file__).resolve().parents[2] / relative.split("/", 1)[1]


def original_rows(relative, expected_sha, expected_count):
    data = reference_path(relative).read_bytes()
    assert hashlib.sha256(data).hexdigest() == expected_sha
    rows = [json.loads(line) for line in data.splitlines()]
    assert len(rows) == expected_count
    return rows


def assert_vector(actual, expected, tolerance):
    assert len(actual) == len(expected)
    for index, (observed, oracle) in enumerate(zip(actual, expected)):
        assert abs(observed - oracle) <= tolerance, (index, observed, oracle, tolerance)


def test_mmff_original_initial_energy_gradient():
    records = original_rows("numerical-reference-inputs/testdata/forcefield/expected/rdkit/smiles_small/forcefield_params.jsonl", "7cf9e84bea3a66ee4d29f198639ed62b6caf0f5b0ad369d7812b5417e91519be", 152)
    counts = {"energy": 0, "gradient": 0}
    for row, record in enumerate(records, 1):
        embedded = record.get("embedded")
        if not embedded or not embedded["ok"]:
            if embedded:
                assert embedded["error"] is not None
            continue
        try:
            molecule = ck.Molecule.from_smiles(embedded["cxsmiles"])
            assert len(molecule.conformers_3d()) == 1
            assert molecule.conformers_3d()[0].coordinates() == embedded["coords"]
            expected = embedded["mmff"]
            if expected["ok"]:
                actual = molecule.with_mmff_optimized_confs_with_params(ck.MmffConformerOptimizationParams(max_iterations=0, non_bonded_threshold=100.0)).conformer_results()[0]
                assert actual.status_code() == expected["needs_more"]
                assert abs(actual.energy()-expected["energy"]) <= 1.0e-6
                counts["energy"] += 1
                if expected.get("gradient") is not None:
                    evaluated = molecule.mmff_energy_gradient_with_params(ck.MmffEvaluationParams(non_bonded_threshold=100.0))
                    assert evaluated is not None
                    assert abs(evaluated.energy()-expected["energy"]) <= 1.0e-6
                    assert_vector(evaluated.gradient(), expected["gradient"], 1.0e-6)
                    counts["gradient"] += 1
            else:
                assert expected["error"] is not None
            assert molecule.conformers_3d()[0].coordinates() == embedded["coords"]
        except Exception as error:
            raise AssertionError((row, record["smiles"])) from error
    assert counts == {"energy": 137, "gradient": 122}
    print(f"MMFF initial numerical original rows: {counts}; energy/gradient1e-6")
