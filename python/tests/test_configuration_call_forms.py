"""Real calls through the registry-installed normalization adapter."""
import pytest
import cosmolkit as ck


@pytest.mark.parametrize("method", ["fingerprint_morgan", "fingerprint_morgan_sparse", "fingerprint_morgan_count", "fingerprint_morgan_sparse_count"])
def test_morgan_short_methods_accept_generator_keywords(method):
    molecule = ck.Molecule.from_smiles("CC(C)O")
    function = getattr(molecule, method)
    params = ck.MorganFingerprintParams(generator=ck.MorganParams(radius=2, fp_size=256, include_chirality=True), from_atoms=[0, 1])
    expected = getattr(molecule, method + "_with_params")(params, None)
    actual = function(radius=2, fp_size=256, include_chirality=True, from_atoms=[0, 1])
    project = lambda value: value.on_bits() if hasattr(value, "on_bits") else value.nonzero_elements()
    assert project(actual) == project(expected) == project(function(params))
    assert project(function(generator=params.generator, from_atoms=[0, 1])) == project(expected)
    if hasattr(actual, "n_bits") and "sparse" not in method:
        assert actual.n_bits() == 256
    with pytest.raises(TypeError, match="mutually exclusive"):
        function(params, radius=2)
    with pytest.raises(TypeError, match="mutually exclusive"):
        function(generator=ck.MorganParams(), radius=2)
    with pytest.raises((TypeError, OverflowError)):
        function(radius=-1)
    with pytest.raises(TypeError):
        function(unknown_morgan_option=True)


def test_search_configuration_forms_are_effective_and_mutually_exclusive():
    molecule = ck.Molecule.from_smiles("CCC")
    query = ck.parse_smarts("C")
    params = ck.SubstructMatchParams(max_matches=1, uniquify=False)
    assert len(molecule.substruct_matches(query)) == 3
    assert len(molecule.substruct_matches(query, params)) == 1
    assert len(molecule.substruct_matches(query, max_matches=1, uniquify=False)) == 1
    with pytest.raises(TypeError, match="mutually exclusive"):
        molecule.substruct_matches(query, params, max_matches=1)
    with pytest.raises(TypeError):
        molecule.substruct_matches(query, not_a_field=True)
    with pytest.raises((TypeError, OverflowError)):
        molecule.substruct_matches(query, max_matches=-1)
    assert molecule.to_smiles() == "CCC"
    assert params.max_matches == 1


def test_smiles_factories_writers_and_required_fragment_inputs():
    params = ck.SmilesParseParams(sanitize=False)
    molecule = ck.Molecule.from_smiles("CCO", params)
    assert molecule.num_atoms() == ck.Molecule.from_smiles("CCO", sanitize=False).num_atoms() == 3
    writer = ck.SmilesWriteParams(canonical=False, all_bonds_explicit=True)
    assert molecule.to_smiles(writer) == molecule.to_smiles(canonical=False, all_bonds_explicit=True)
    fragment = ck.FragmentSmilesWriteParams([0, 1])
    assert molecule.to_fragment_smiles(fragment) == molecule.to_fragment_smiles(atoms=[0, 1]) == "CC"
    with pytest.raises(TypeError):
        molecule.to_fragment_smiles()


def test_two_configuration_objects_have_independent_keyword_groups():
    batch = ck.MoleculeBatch.from_smiles_list(["CCO", "CCN"], ck.SmilesParseParams(), ck.BatchParams(n_jobs=2))
    assert batch.to_smiles_list(n_jobs=2, canonical=False) == ["CCO", "CCN"]
    assert ck.MoleculeBatch.from_smiles_list(["CCO"], sanitize=False, n_jobs=2).to_smiles_list() == ["CCO"]
    with pytest.raises(TypeError, match="mutually exclusive"):
        ck.MoleculeBatch.from_smiles_list(["C"], parse=ck.SmilesParseParams(), sanitize=False)


def test_class_factories_keep_defaults_outside_the_configuration_record(tmp_path):
    text = ck.Molecule.from_smiles("C").to_sdf()
    assert len(ck.MoleculeBatch.from_sdf_records(text, ck.SdfReadParams())) == 1
    assert len(ck.MoleculeBatch.from_sdf_records(text, sanitize=False)) == 1
    path = tmp_path / "one.sdf"
    path.write_text(text)
    assert len(ck.SdfDataset.open(str(path), ck.SdfReadParams())) == 1
    assert len(ck.SdfDataset.open(str(path), sanitize=False)) == 1


def test_sdf_text_conversion_and_batch_file_writers_are_distinct(tmp_path):
    molecule = ck.Molecule.from_smiles("CO")
    text = molecule.to_sdf()
    assert isinstance(text, str) and text.endswith("$$$$\n")
    batch = ck.MoleculeBatch.from_sdf_records(text)
    original = batch.to_smiles_list()
    params = ck.BatchExportParams(format="v3000", progress_bar=False)

    for label, args, kwargs in (
        ("default", (), {}),
        ("params", (params,), {}),
        ("keywords", (), {"format": "v3000", "progress_bar": False}),
    ):
        path = tmp_path / f"{label}.sdf"
        report = batch.write_sdf(str(path), *args, **kwargs)
        assert isinstance(report, ck.BatchExportReport)
        assert (report.total(), report.success(), report.failed()) == (1, 1, 0)
        assert ck.MoleculeBatch.read_sdf(str(path)).to_smiles_list() == original
        if label != "default":
            assert "V3000" in path.read_text()

        directory = tmp_path / label
        directory.mkdir()
        report = batch.write_sdf_files(str(directory), *args, **kwargs)
        assert report.success() == 1
        files = list(directory.glob("*.sdf"))
        assert len(files) == 1
        assert ck.MoleculeBatch.read_sdf(str(files[0])).to_smiles_list() == original
        if label != "default":
            assert "V3000" in files[0].read_text()

    for name in ("to_sdf", "to_sdf_with_params", "to_sdf_files", "to_sdf_files_with_params"):
        assert not hasattr(batch, name)
    assert batch.to_smiles_list() == original
    assert molecule.to_sdf() == text


def test_registered_value_equality_keeps_python_reflected_comparison():
    class Reflected:
        def __eq__(self, other):
            return "reflected"
    value = ck.TautomerScoreTerm("carbon", "[#6]", 3)
    assert value == ck.TautomerScoreTerm("carbon", "[#6]", 3)
    assert value != ck.TautomerScoreTerm("carbon", "[#6]", 4)
    assert value.__eq__(Reflected()) is NotImplemented
    assert (value == Reflected()) == "reflected"


def test_pickle_protocol_accepts_the_declared_index_protocol():
    class Protocol:
        def __index__(self):
            return 5
    molecule = ck.Molecule.from_smiles("CCO")
    constructor, arguments = molecule.__reduce_ex__(Protocol())
    assert constructor(*arguments).to_smiles() == "CCO"
    with pytest.raises(TypeError):
        molecule.__reduce_ex__(1.5)


def test_native_batch_record_union_keeps_both_real_variants():
    batch = ck.MoleculeBatch.from_smiles_list(["C", "C#"], errors=ck.BatchErrorMode.KEEP)
    records = batch.records()
    assert len(records) == 2
    assert isinstance(records[0], ck.Molecule)
    assert isinstance(records[1], ck.BatchError)
    assert isinstance(batch[0], ck.Molecule)
    assert batch[1] is None
