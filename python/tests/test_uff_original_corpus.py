"""All 137 initial energy/gradient source rows via actual native UFF APIs."""
import cosmolkit as ck
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

# Resolve the existing frozen fixture helper by its exact sibling path, also
# when this file is imported as a module. No chemistry or fixture copy occurs.
_fixture_path = Path(__file__).with_name("test_mmff_original_corpus.py")
_fixture_spec = spec_from_file_location("_forcefield_original_fixture", _fixture_path)
assert _fixture_spec is not None and _fixture_spec.loader is not None
_fixture_module = module_from_spec(_fixture_spec)
_fixture_spec.loader.exec_module(_fixture_module)
original_rows = _fixture_module.original_rows
assert_vector = _fixture_module.assert_vector


def test_uff_original_152_inputs_all_137_initial_energy_gradient():
    records = original_rows("numerical-reference-inputs/testdata/forcefield/expected/rdkit/smiles_small/forcefield_params.jsonl", "7cf9e84bea3a66ee4d29f198639ed62b6caf0f5b0ad369d7812b5417e91519be", 152)
    counts = {"energy": 0, "gradient": 0}
    for row, record in enumerate(records, 1):
        embedded = record.get("embedded")
        if not embedded or not embedded["ok"]:
            if embedded:
                assert embedded["error"] is not None
            continue
        try:
            molecule = ck.Molecule.from_smiles(embedded["cxsmiles"]).with_assigned_valence()
            assert len(molecule.conformers_3d()) == 1
            assert molecule.conformers_3d()[0].coordinates() == embedded["coords"]
            expected = embedded["uff"]
            assert expected["ok"] is True and expected["error"] is None
            result = molecule.with_uff_optimized_with_params(ck.UffOptimizationParams(max_iterations=0, vdw_threshold=100.0))
            assert result.status_code() == expected["needs_more"]
            assert abs(result.energy()-expected["energy"]) <= 1.0e-6
            counts["energy"] += 1
            evaluated = molecule.uff_energy_gradient_with_params(ck.UffEvaluationParams(vdw_threshold=100.0))
            assert abs(evaluated.energy()-expected["energy"]) <= 1.0e-6
            assert_vector(evaluated.gradient(), expected["gradient"], 1.0e-6)
            counts["gradient"] += 1
            assert molecule.conformers_3d()[0].coordinates() == embedded["coords"]
        except Exception as error:
            raise AssertionError((row, record["smiles"])) from error
    assert counts == {"energy": 137, "gradient": 137}
    print(f"UFF initial numerical original rows: {counts}; original energy/gradient1e-6")
