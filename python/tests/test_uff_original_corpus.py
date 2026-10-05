"""Original UFF5000 availability and all 137 numeric source rows via actual native APIs."""
import hashlib
import json
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
reference_path = _fixture_module.reference_path
original_rows = _fixture_module.original_rows
assert_vector = _fixture_module.assert_vector
assert_coordinates = _fixture_module.assert_coordinates


def test_uff_original_5000_implicit_and_explicit_h_source_availability():
    data = reference_path("reference-inputs/testdata/smiles/corpus/smiles_5000.smi").read_bytes()
    assert hashlib.sha256(data).hexdigest() == "a4d579cd72621af27772256bb23ba796452276bb924fd20aac83625ffa67d849"
    corpus = [line.strip() for line in data.decode().splitlines() if line.strip() and not line.startswith("#")]
    assert len(corpus) == 5000
    records = original_rows("reference-inputs/testdata/forcefield_coverage/expected/rdkit/smiles_5000/forcefield_coverage.jsonl", "0c61a5dd6d39c2503346c518c0f4d4d1b7d9eb9e048797d7ea837ec08de1f8e1", 5000)
    compared = 0
    for row, record in enumerate(records, 1):
        assert record["smiles"] == corpus[row-1]
        assert record["rdkit_ok"] is True
        try:
            source = ck.Molecule.from_smiles(record["smiles"])
            original = source.to_smiles()
            molecule = source.with_assigned_valence()
            expected = record["uff"]
            assert expected["ok"] is True and expected["error"] is None
            assert molecule.uff_has_all_molecule_params() == expected["has_all"]
            expanded = source.with_hydrogens().with_assigned_valence()
            expected = record["uff_explicit_h"]
            assert expected["ok"] is True and expected["error"] is None
            assert expanded.uff_has_all_molecule_params() == expected["has_all"]
            assert source.to_smiles() == original
            compared += 1
        except Exception as error:
            raise AssertionError((row, record["smiles"])) from error
    assert compared == 5000
    print(f"UFF5000: rows=5000 compared={compared}, original implicit+explicit H availability")


def test_uff_original_152_inputs_all_137_energy_gradient_single_multi():
    records = original_rows("numerical-reference-inputs/testdata/forcefield/expected/rdkit/smiles_small/forcefield_params.jsonl", "7cf9e84bea3a66ee4d29f198639ed62b6caf0f5b0ad369d7812b5417e91519be", 152)
    counts = {"energy": 0, "gradient": 0, "single": 0, "multi": 0}
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
            expected = embedded["uff_optimized"]
            assert expected["ok"] is True and expected["error"] is None
            single = molecule.with_uff_optimized_with_params(ck.UffOptimizationParams(max_iterations=200,vdw_threshold=100.0))
            assert single.status_code() == expected["needs_more"]
            assert abs(single.energy()-expected["energy"]) <= 1.0e-6
            assert_coordinates(single.molecule().conformers_3d()[0].coordinates(), expected["coords"])
            counts["single"] += 1
            expected = embedded["uff_multi_optimized"]
            assert expected["ok"] is True and expected["error"] is None
            # Original generator uses canonicalCX topology and source-ordered
            # binary64 initial rows. Preserve every dative/CX annotation, replacing
            # only rounded coordinate text; builder appends original later rows.
            prefix, separator, tail = embedded["cxsmiles"].partition(" |(")
            assert separator
            _, closing, suffix = tail.partition(")")
            assert closing
            first_rows = expected["initial_coords"][0]
            coordinates_text = ";".join(",".join(repr(value) for value in point) for point in first_rows)
            initial = ck.Molecule.from_smiles(prefix + " |(" + coordinates_text + ")" + suffix)
            assert initial.conformers_3d()[0].coordinates() == first_rows
            builder = initial.to_builder()
            for rows in expected["initial_coords"][1:]:
                builder.add_3d_conformer(rows)
            multi_molecule = builder.build().with_assigned_valence()
            actual = multi_molecule.with_uff_optimized_confs_with_params(ck.UffConformerOptimizationParams(max_iterations=200,vdw_threshold=100.0))
            results = actual.conformer_results()
            conformers = actual.molecule().conformers_3d()
            assert len(results) == len(conformers) == len(expected["conformer_results"])
            for result, conformer, oracle in zip(results, conformers, expected["conformer_results"]):
                assert result.status_code() == oracle["needs_more"]
                assert abs(result.energy()-oracle["energy"]) <= 1.0e-6
                assert_coordinates(conformer.coordinates(), oracle["coords"])
            counts["multi"] += 1
            assert molecule.conformers_3d()[0].coordinates() == embedded["coords"]
        except Exception as error:
            raise AssertionError((row, record["smiles"])) from error
    assert counts == {"energy": 137, "gradient": 137, "single": 137, "multi": 137}
    print(f"UFF numerical original rows: {counts}; original energy/gradient/coordinates1e-6")
