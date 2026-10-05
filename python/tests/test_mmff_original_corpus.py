"""Original source-bound MMFF oracles through canonical native public APIs."""
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
    for index, (observed, oracle) in enumerate(zip(actual, expected, strict=True)):
        assert abs(observed - oracle) <= tolerance, (index, observed, oracle, tolerance)


def assert_coordinates(actual, expected):
    assert len(actual) == len(expected)
    for observed, oracle in zip(actual, expected, strict=True):
        assert_vector(observed, oracle, 1.0e-6)


def assert_properties(molecule, expected):
    assert expected["ok"], expected["error"]
    assert expected["has_all"] is not None
    assert molecule.mmff_has_all_molecule_params() == expected["has_all"]
    if expected["atom_types"] is None:
        assert expected["has_all"] is False
        assert expected["formal_charges"] is None and expected["partial_charges"] is None
        return
    props = molecule.mmff_properties()
    assert [props.atom_type(i) for i in range(molecule.num_atoms())] == expected["atom_types"]
    assert_vector([props.formal_charge(i) for i in range(molecule.num_atoms())], expected["formal_charges"], 1.0e-12)
    assert_vector([props.partial_charge(i) for i in range(molecule.num_atoms())], expected["partial_charges"], 1.0e-12)
    atoms = props.atoms()
    assert [a.atom_type() for a in atoms] == expected["atom_types"]
    assert_vector([a.formal_charge() for a in atoms], expected["formal_charges"], 1.0e-12)
    assert_vector([a.partial_charge() for a in atoms], expected["partial_charges"], 1.0e-12)


def test_mmff_properties_original_5000_implicit_and_explicit_h():
    data = reference_path("reference-inputs/testdata/smiles/corpus/smiles_5000.smi").read_bytes()
    assert hashlib.sha256(data).hexdigest() == "a4d579cd72621af27772256bb23ba796452276bb924fd20aac83625ffa67d849"
    corpus = [line.strip() for line in data.decode().splitlines() if line.strip() and not line.startswith("#")]
    assert len(corpus) == 5000
    records = original_rows("reference-inputs/testdata/forcefield_coverage/expected/rdkit/smiles_5000/forcefield_coverage.jsonl", "0c61a5dd6d39c2503346c518c0f4d4d1b7d9eb9e048797d7ea837ec08de1f8e1", 5000)
    compared = 0
    source_errors = 0
    for row, record in enumerate(records, 1):
        assert record["smiles"] == corpus[row-1]
        if not record["rdkit_ok"]:
            assert record["error"] is not None
            source_errors += 1
            continue
        try:
            molecule = ck.Molecule.from_smiles(record["smiles"])
            original = molecule.to_smiles()
            assert_properties(molecule, record["mmff"])
            expanded = molecule.with_hydrogens()
            assert_properties(expanded, record["mmff_explicit_h"])
            assert molecule.to_smiles() == original
            compared += 1
        except Exception as error:
            raise AssertionError((row, record["smiles"])) from error
    assert compared + source_errors == 5000
    assert compared > 0
    print(f"MMFF5000: rows=5000 compared={compared} source_errors={source_errors}, implicit+explicit H, exact types, charges1e-12")


def test_mmff_original_numerical_energy_gradient_single_multi():
    records = original_rows("numerical-reference-inputs/testdata/forcefield/expected/rdkit/smiles_small/forcefield_params.jsonl", "7cf9e84bea3a66ee4d29f198639ed62b6caf0f5b0ad369d7812b5417e91519be", 152)
    counts = {"energy": 0, "gradient": 0, "single": 0, "multi": 0}
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
            expected = embedded["mmff_optimized"]
            if expected["ok"]:
                assert expected["error"] is None
                single = molecule.with_mmff_optimized_with_params(ck.MmffOptimizationParams(max_iterations=200,non_bonded_threshold=100.0))
                energy_result = molecule.with_mmff_optimized_confs_with_params(ck.MmffConformerOptimizationParams(max_iterations=200,non_bonded_threshold=100.0)).conformer_results()[0]
                assert single.status_code() == energy_result.status_code() == expected["needs_more"]
                assert abs(energy_result.energy()-expected["energy"]) <= 1.0e-6
                if expected["coords"] is None:
                    assert expected["needs_more"] == -1 and expected["energy"] == -1.0
                    assert_coordinates(single.molecule().conformers_3d()[0].coordinates(), embedded["coords"])
                else:
                    assert_coordinates(single.molecule().conformers_3d()[0].coordinates(), expected["coords"])
                counts["single"] += 1
            expected = embedded["mmff_multi_optimized"]
            if expected["ok"]:
                assert expected["error"] is None
                # Source initial rows use original SMILES ordering. Drop the
                # CX coordinate suffix before rebuilding exactly those rows.
                builder = ck.Molecule.from_smiles(record["smiles"]).to_builder()
                for rows in expected["initial_coords"]:
                    builder.add_3d_conformer(rows)
                multi_molecule = builder.build()
                actual = multi_molecule.with_mmff_optimized_confs_with_params(ck.MmffConformerOptimizationParams(max_iterations=200,non_bonded_threshold=100.0))
                results = actual.conformer_results()
                conformers = actual.molecule().conformers_3d()
                assert len(results) == len(conformers) == len(expected["conformer_results"])
                for result, conformer, oracle in zip(results, conformers, expected["conformer_results"], strict=True):
                    assert result.status_code() == oracle["needs_more"]
                    assert abs(result.energy()-oracle["energy"]) <= 1.0e-6
                    assert_coordinates(conformer.coordinates(),oracle["coords"])
                counts["multi"] += 1
            assert molecule.conformers_3d()[0].coordinates() == embedded["coords"]
        except Exception as error:
            raise AssertionError((row, record["smiles"])) from error
    assert counts == {"energy": 137, "gradient": 122, "single": 137, "multi": 137}
    print(f"MMFF numerical original rows: {counts}; energy/gradient/coordinates1e-6")
