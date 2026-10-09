"""Runtime regressions for registered Python accessors and thin bindings."""

from __future__ import annotations

import ast
from pathlib import Path
from typing import cast

import cosmolkit as ck
import pytest


def _record(smiles: str) -> str:
    return ck.Molecule.from_smiles(smiles).to_sdf()


def test_stereoisomer_constructor_properties_and_readonly_iterator() -> None:
    options = ck.StereoisomerOptions()
    assert (options.try_embedding, options.only_unassigned, options.max_isomers,
            options.random_source, options.unique, options.only_stereo_groups) == (
                False, True, 1024, None, True, False)
    options.try_embedding = True
    options.only_unassigned = False
    options.max_isomers = 2
    options.random_source = 123
    options.unique = False
    options.only_stereo_groups = True
    assert (options.try_embedding, options.only_unassigned, options.max_isomers,
            options.random_source, options.unique, options.only_stereo_groups) == (
                True, False, 2, 123, False, True)
    molecule = ck.Molecule.from_smiles("FC(Cl)Br")
    iterator = molecule.enumerate_stereoisomers()
    assert iterator.yielded_count == 0
    assert iterator.next() is not None
    assert iterator.yielded_count == 1
    assert len(list(iterator)) == 1
    assert iterator.yielded_count == 2
    with pytest.raises(AttributeError):
        setattr(iterator, "yielded_count", 0)


def test_big_integer_count_and_random_seed_do_not_truncate() -> None:
    molecule = ck.Molecule.from_smiles("Br" + "[CH](Cl)" * 70 + "F")
    assert molecule.stereoisomer_count() == 1 << 70
    assert molecule.stereoisomer_count_with_options(ck.StereoisomerOptions()) == 1 << 70
    small = ck.Molecule.from_smiles("FC(Cl)Br")
    for seed in [1 << 100, -(1 << 100)]:
        source = ck.StereoisomerRandomSource.from_integer_seed(seed)
        options = ck.StereoisomerOptions(max_isomers=1, random_source=source)
        first = [m.to_smiles() for m in small.enumerate_stereoisomers_with_options(options)]
        repeat = ck.StereoisomerOptions(max_isomers=1, random_source=ck.StereoisomerRandomSource.from_integer_seed(seed))
        assert first == [m.to_smiles() for m in small.enumerate_stereoisomers_with_options(repeat)]
        assert len(first) == 1


def test_property_text_and_exception_type_are_real_python_values() -> None:
    molecule = ck.Molecule.from_smiles("F[C@](Cl)(Br)I")
    assert molecule.atom_property_string(1, "_CIPCode") == "S"
    assert molecule.atom_property_string(1, "absent") is None
    assert molecule.atom_property_string(999, "_CIPCode") is None
    assert molecule.bond_property_string(0, "absent") is None
    assert molecule.bond_property_string(999, "absent") is None
    value = molecule.property("__computedProps")
    assert value is not None
    # Pinned RDKit GetProp("__computedProps") for this molecule yields this text.
    # Every currently modeled property kind is supported by the Rust formatter;
    # do not invent an UnsupportedKind input that cannot exist in this model.
    assert ck.property_value_to_text(value) == "[numArom,_StereochemDone]"
    assert value.kind() == ck.PropertyValueKind.StringVector
    assert issubclass(ck.PropertyStringError, ValueError)
    assert callable(ck.PropertyStringError.kind)


def test_batch_record_projection_keeps_error_records_and_input_values() -> None:
    original = ck.Molecule.from_smiles("CCO")
    before = original.to_smiles()
    failed = ck.MoleculeBatch.from_smiles_list_with_params(["not-a-smiles"], ck.SmilesParseParams(), ck.BatchParams(errors="keep"))
    error = failed.errors()[0]
    batch = ck.MoleculeBatch.from_records([original, error], ck.BatchErrorMode.KEEP)
    records = batch.records()
    assert len(records) == 2
    assert isinstance(records[0], ck.Molecule)
    assert records[0].to_smiles() == before
    assert isinstance(records[1], ck.BatchError)
    assert (records[1].index(), records[1].operation(), records[1].message()) == (error.index(), error.operation(), error.message())
    assert isinstance(batch.get(1), ck.BatchError)
    assert batch.get(2) is None
    assert batch.valid_mask() == [True, False]
    assert original.to_smiles() == before
    with pytest.raises(ck.BatchValidationError):
        _ = ck.MoleculeBatch.from_records([original, error], ck.BatchErrorMode.RAISE)
    with pytest.raises(TypeError):
        _ = ck.MoleculeBatch.from_records([object()], ck.BatchErrorMode.KEEP)  # pyright: ignore[reportArgumentType]


def test_typed_sdf_factories_exports_and_report(tmp_path: Path) -> None:
    text = _record("CCO") + _record("CCN")
    read = ck.SdfReadParams()
    batch = ck.MoleculeBatch.from_sdf_records_with_params(text, read, ck.BatchErrorMode.RAISE, 2)
    assert len(batch) == 2
    path = tmp_path / "molecules.sdf"
    _ = path.write_text(text)
    loaded = ck.MoleculeBatch.read_sdf_with_params(str(path), read, ck.BatchErrorMode.RAISE, 2, False)
    assert loaded.to_smiles_list() == batch.to_smiles_list()
    dataset = ck.SdfDataset.open(str(path))
    reordered = ck.MoleculeBatch.from_dataset_indices(dataset, [1, 0, 1], ck.BatchErrorMode.RAISE)
    assert reordered.to_smiles_list() == ["CCN", "CCO", "CCN"]
    invalid_text = _record("CCO") + "invalid molecule\n$$$$\n"
    kept = ck.MoleculeBatch.from_sdf_records_with_params(invalid_text, read, ck.BatchErrorMode.KEEP, 2)
    assert kept.valid_mask() == [True, False]
    with pytest.raises(ck.BatchValidationError):
        _ = ck.MoleculeBatch.from_sdf_records_with_params(invalid_text, read, ck.BatchErrorMode.RAISE, 2)
    # Fixed RDKit MolToMolBlock(AddHs(MolFromSmiles("CO"))) input. Test
    # read-parameter forwarding independently of CK's hydrogen-addition owner.
    explicit_text = "\n     RDKit          2D\n\n  6  5  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.5000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    1.5000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0000   -1.5000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.2500   -1.2990    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  1  3  1  0\n  1  4  1  0\n  1  5  1  0\n  2  6  1  0\nM  END\n$$$$\n"
    retained = ck.MoleculeBatch.from_sdf_records_with_params(explicit_text, ck.SdfReadParams(remove_hs=False), ck.BatchErrorMode.RAISE, 2).get(0)
    removed = ck.MoleculeBatch.from_sdf_records_with_params(explicit_text, ck.SdfReadParams(remove_hs=True), ck.BatchErrorMode.RAISE, 2).get(0)
    assert isinstance(retained, ck.Molecule) and isinstance(removed, ck.Molecule)
    assert retained.num_atoms() == 6
    assert removed.num_atoms() == 2
    options = ck.BatchExportParams(format="v3000", errors="raise", n_jobs=2, progress_bar=False)
    assert (options.format, cast(object, options.errors), options.n_jobs, options.progress_bar) == ("v3000", ck.BatchErrorMode.RAISE, 2, False)
    output = tmp_path / "output.sdf"
    report_file = tmp_path / "report.json"
    report = batch.write_sdf_with_params(str(output), options, str(report_file))
    assert (report.total(), report.success(), report.failed()) == (2, 2, 0)
    assert "V3000" in output.read_text()
    copied_report = tmp_path / "copied-report.json"
    report.write_report(str(copied_report))
    assert copied_report.read_bytes() == report_file.read_bytes()
    directory = tmp_path / "separate"
    directory.mkdir()
    files = batch.write_sdf_files_with_params(str(directory), options, ["one.sdf", "two.sdf"], None)
    assert files.success() == 2
    assert sorted(p.name for p in directory.iterdir()) == ["one.sdf", "two.sdf"]
    with pytest.raises(ValueError):
        _ = ck.BatchExportParams(n_jobs=0)
    with pytest.raises(ck.BatchValidationError):
        _ = batch.write_sdf_with_params(str(tmp_path), options, None)
    with pytest.raises(ck.BatchValidationError):
        _ = ck.MoleculeBatch.read_sdf_with_params(str(tmp_path / "absent.sdf"), read, ck.BatchErrorMode.RAISE, None, False)
    # Typed entrypoints must preserve the Rust error, not replace it with a
    # Python-only prevalidation error before reaching the facade.
    with pytest.raises(ck.BatchValidationError):
        _ = ck.MoleculeBatch.read_sdf_with_params(str(path), read, ck.BatchErrorMode.RAISE, 0, False)


def test_batch_iterators_and_stream_transfer_continue_at_current_record(tmp_path: Path) -> None:
    path = tmp_path / "molecules.sdf"
    _ = path.write_text(_record("C") + _record("CC") + _record("CCC"))
    dataset = ck.SdfDataset.open(str(path))
    indexed = dataset.batches(size=2)
    first = indexed.next_batch()
    assert first is not None and first.to_smiles_list() == ["C", "CC"]
    second = next(indexed)
    assert second.to_smiles_list() == ["CCC"]
    assert indexed.next_batch() is None
    forward = ck.SdfReader.open(str(path)).batches(size=2)
    assert forward.next_batch() is not None
    assert next(forward).to_smiles_list() == ["CCC"]
    assert forward.next_batch() is None
    stream = ck.SdfRecordStream.open(str(path))
    assert stream.next_record() is not None
    assert stream.records_consumed() == 1
    batches = stream.batches(2, ck.BatchErrorMode.RAISE, 2)
    remaining = batches.next_batch()
    assert remaining is not None and remaining.to_smiles_list() == ["CC", "CCC"]
    assert batches.next_batch() is None
    for call in [stream.next_record, stream.is_end, stream.records_consumed, stream.bytes_consumed, stream.lines_consumed]:
        with pytest.raises(RuntimeError, match="transferred"):
            _ = call()
    with pytest.raises(RuntimeError, match="transferred"):
        _ = stream.batches(1, ck.BatchErrorMode.RAISE, None)
    invalid = ck.SdfRecordStream.open(str(path))
    with pytest.raises(ck.BatchValidationError):
        _ = invalid.batches(0, ck.BatchErrorMode.RAISE, None)


def test_torsion_batch_thin_bindings_match_scalar_and_propagate_callback_errors() -> None:
    molecules = [ck.Molecule.from_smiles(s) for s in ["CCCC", "CCCO", "CC"]]
    before = [m.to_smiles() for m in molecules]
    batch = ck.MoleculeBatch.from_records(molecules, ck.BatchErrorMode.RAISE)
    generator = ck.TopologicalTorsionFingerprintGenerator()
    expected = [m.fingerprint_topological_torsion_with_generator(generator).on_bits() for m in molecules]
    assert [v.on_bits() for v in batch.fingerprint_topological_torsion_list() if v is not None] == expected
    generator_params = ck.TopologicalTorsionParams(fp_size=256, count_simulation=False)
    options = ck.TopologicalTorsionFingerprintParams(generator=generator_params)
    calls: list[None] = []
    values = batch.fingerprint_topological_torsion_list_with_params(options, ck.BatchQueryParams(n_jobs=2, progress_callback=lambda: calls.append(None)))
    expected_custom = [m.fingerprint_topological_torsion_with_params(options, None).on_bits() for m in molecules]
    assert [v.on_bits() for v in values if v is not None] == expected_custom
    assert len(values) == 3 and len(calls) == 3
    def fail() -> None:
        raise LookupError("callback sentinel")
    with pytest.raises(LookupError, match="callback sentinel"):
        _ = batch.fingerprint_topological_torsion_list_with_params(options, ck.BatchQueryParams(progress_callback=fail))
    assert [m.to_smiles() for m in molecules] == before


def test_generated_big_integer_annotations_are_python_int() -> None:
    tree = ast.parse((Path(__file__).resolve().parents[1] / "cosmolkit.pyi").read_text())
    classes = {node.name: node for node in tree.body if isinstance(node, ast.ClassDef)}
    for name in ["stereoisomer_count", "stereoisomer_count_with_options"]:
        method = next(node for node in classes["Molecule"].body if isinstance(node, ast.FunctionDef) and node.name == name)
        assert method.returns is not None
        assert ast.unparse(method.returns) == "builtins.int"
    factory = next(node for node in classes["StereoisomerRandomSource"].body if isinstance(node, ast.FunctionDef) and node.name == "from_integer_seed")
    annotation = factory.args.args[0].annotation
    assert annotation is not None
    assert ast.unparse(annotation) == "builtins.int"
    for owner, result in [("StereoisomerIterator", "Molecule"), ("SdfRecordStream", "SdfRecord")]:
        method = next(node for node in classes[owner].body if isinstance(node, ast.FunctionDef) and node.name == "__next__")
        assert method.returns is not None and ast.unparse(method.returns) == result


def test_dynamic_batch_types_have_stubs_and_preserve_error_fields() -> None:
    tree = ast.parse((Path(__file__).resolve().parents[1] / "cosmolkit.pyi").read_text())
    for name in ["BatchErrorMode", "BatchValidationError", "PropertyStringError"]:
        assert len([node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == name]) == 1
    assert int(ck.BatchErrorMode.RAISE) == 1 and int(ck.BatchErrorMode.KEEP) == 2
    failed = ck.MoleculeBatch.from_smiles_list_with_params(["invalid"], ck.SmilesParseParams(), ck.BatchParams(errors="keep"))
    record = failed.errors()[0]
    error = ck.BatchValidationError("message", 1, "reason", [record])
    assert str(error) == "message"
    assert error.error_count == 1 and error.reason == "reason"
    assert error.errors()[0].message() == record.message()


@pytest.mark.parametrize("factory,changes", [
    ("MorganFingerprintGenerator", {"radius": 4, "only_nonzero_invariants": True,
        "include_redundant_environments": True, "include_chirality": True,
        "count_simulation": False, "fp_size": 256, "bits_per_feature": 2, "count_bounds": [2, 4]}),
    ("TopologicalTorsionFingerprintGenerator", {"torsion_atom_count": 5,
        "only_shortest_paths": True, "include_chirality": True,
        "count_simulation": False, "fp_size": 256, "bits_per_feature": 2, "count_bounds": [2, 4]}),
])
def test_fingerprint_settings_properties_and_explicit_setters_are_distinct(factory: str, changes: dict[str, object]) -> None:
    settings = ck.MorganFingerprintGenerator().settings() if factory == "MorganFingerprintGenerator" else ck.TopologicalTorsionFingerprintGenerator().settings()
    for name, value in changes.items():
        before = cast(object, getattr(settings, name))
        method = cast(object, getattr(settings, "set_" + name))
        assert callable(method)
        _ = method(value)
        assert cast(object, getattr(settings, name)) == value
        setattr(settings, name, before)
        assert cast(object, getattr(settings, name)) == before
