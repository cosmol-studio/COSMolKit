"""UFF21 original row mapping and source options through the canonical native API."""
import math
import pytest
import cosmolkit as ck


def make_molecule():
    builder = ck.Molecule.from_smiles("CC").to_builder()
    builder.add_3d_conformer([[0., 0., 0.], [2., 0., 0.]])
    return builder.build().with_assigned_valence()


def coordinates(molecule):
    return [row.coordinates() for row in molecule.conformers_3d()]


BASELINE_ROWS = [
    "Molecule.has_uff_params", "Molecule.with_uff_optimized", "Molecule.with_uff_optimized_confs",
    "UffOptimizeMoleculeConfResult", "UffOptimizeMoleculeConfResult.__repr__",
    "UffOptimizeMoleculeConfResult.energy", "UffOptimizeMoleculeConfResult.needs_more", "UffOptimizeMoleculeConfResult.status_code",
    "UffOptimizeMoleculeConfsResult", "UffOptimizeMoleculeConfsResult.__repr__", "UffOptimizeMoleculeConfsResult.conformer_results", "UffOptimizeMoleculeConfsResult.molecule",
    "UffOptimizeMoleculeResult", "UffOptimizeMoleculeResult.__repr__", "UffOptimizeMoleculeResult.energy", "UffOptimizeMoleculeResult.molecule", "UffOptimizeMoleculeResult.needs_more", "UffOptimizeMoleculeResult.status_code",
    "uff_has_all_molecule_params", "uff_optimize_molecule", "uff_optimize_molecule_confs",
]


@pytest.mark.parametrize("original_row", BASELINE_ROWS, ids=BASELINE_ROWS)
def test_uff_original_21_rows_canonical_native_behavior(original_row):
    molecule = make_molecule()
    before = coordinates(molecule)
    if original_row in ["Molecule.has_uff_params", "uff_has_all_molecule_params"]:
        assert molecule.uff_has_all_molecule_params() is True
        assert ck.Molecule.from_smiles("[He]").with_assigned_valence().uff_has_all_molecule_params() is False
    elif original_row in ["Molecule.with_uff_optimized", "uff_optimize_molecule"] or original_row.startswith("UffOptimizeMoleculeResult"):
        result = molecule.with_uff_optimized()
        assert isinstance(result, ck.UffOptimizationResult)
        assert result.status_code() == 0 and result.needs_more() is False
        assert isinstance(result.energy(), float) and math.isfinite(result.energy())
        assert result.molecule().num_atoms() == 2
        assert coordinates(result.molecule()) != before
        assert "UffOptimizationResult(needs_more=0, energy=" in repr(result)
    else:
        multi = molecule.with_uff_optimized_confs()
        assert isinstance(multi, ck.UffConformerOptimizationResult)
        assert multi.molecule().num_atoms() == 2
        assert coordinates(multi.molecule()) != before
        assert repr(multi) == "UffConformerOptimizationResult(conformers=1)"
        assert len(multi.conformer_results()) == 1
        result = multi.conformer_results()[0]
        assert isinstance(result, ck.UffConformerResult)
        assert result.status_code() == 0 and result.needs_more() is False
        assert isinstance(result.energy(), float) and math.isfinite(result.energy())
        assert result.conformer_id() == molecule.conformers_3d()[0].id()
        assert "UffConformerResult(needs_more=0, energy=" in repr(result)
    assert coordinates(molecule) == before


def test_uff_source_defaults_immutable_parameter_values():
    single = ck.UffOptimizationParams()
    multi = ck.UffConformerOptimizationParams()
    evaluation = ck.UffEvaluationParams()
    assert (single.max_iterations, single.vdw_threshold, single.conformer_id, single.ignore_interfragment_interactions) == (1000, 10., None, True)
    assert (multi.num_threads, multi.max_iterations, multi.vdw_threshold, multi.ignore_interfragment_interactions) == (1, 1000, 10., True)
    assert (evaluation.vdw_threshold, evaluation.conformer_id, evaluation.ignore_interfragment_interactions) == (10., None, True)
    for value in [single, multi, evaluation]:
        with pytest.raises(AttributeError):
            setattr(value, "vdw_threshold", 1.)


@pytest.mark.parametrize("selector", [-1, -2, -3, -2**31])
def test_uff_source_signed_selectors_choose_first_stored_3d(selector):
    molecule = make_molecule()
    before = coordinates(molecule)
    selected = molecule.with_uff_optimized_with_params(ck.UffOptimizationParams(max_iterations=0, conformer_id=selector))
    default = molecule.with_uff_optimized_with_params(ck.UffOptimizationParams(max_iterations=0))
    assert selected.status_code() == default.status_code() == 1
    assert selected.needs_more() is True
    assert selected.energy() == default.energy()
    assert coordinates(selected.molecule()) == before
    actual = molecule.uff_energy_gradient_with_params(ck.UffEvaluationParams(conformer_id=selector))
    first = molecule.uff_energy_gradient()
    assert actual.energy() == first.energy() == selected.energy()
    assert actual.gradient() == first.gradient()
    assert coordinates(molecule) == before


@pytest.mark.parametrize("selector", [-2**31-1, 2**31])
@pytest.mark.parametrize("parameter", [ck.UffOptimizationParams, ck.UffEvaluationParams])
def test_uff_source_selector_i32_conversion_bounds(selector, parameter):
    with pytest.raises(OverflowError):
        parameter(conformer_id=selector)


@pytest.mark.parametrize("threads", [1, 2, 0, -1, -100000])
def test_uff_existing_source_dispatch_preserves_order_rows_and_source(threads):
    molecule = make_molecule()
    builder = molecule.to_builder()
    builder.add_3d_conformer([[0., 0., 0.], [2.4, 0., 0.]])
    molecule = builder.build().with_assigned_valence()
    before = coordinates(molecule)
    serial = molecule.with_uff_optimized_confs_with_params(ck.UffConformerOptimizationParams(num_threads=1, max_iterations=25))
    actual = molecule.with_uff_optimized_confs_with_params(ck.UffConformerOptimizationParams(num_threads=threads, max_iterations=25))
    assert [r.conformer_id() for r in actual.conformer_results()] == [r.id() for r in molecule.conformers_3d()]
    assert [(r.status_code(), r.energy()) for r in actual.conformer_results()] == [(r.status_code(), r.energy()) for r in serial.conformer_results()]
    assert coordinates(actual.molecule()) == coordinates(serial.molecule())
    assert coordinates(molecule) == before


def test_uff_source_undefined_thread_signed_negation_is_typed_atomic():
    molecule = make_molecule()
    before = coordinates(molecule)
    with pytest.raises(ck.OperationError) as error:
        molecule.with_uff_optimized_confs_with_params(ck.UffConformerOptimizationParams(num_threads=-2**31))
    assert error.value.kind == "UffOptimization"
    assert isinstance(error.value.__cause__, ck.UffOptimizationError)
    assert "UndefinedSignedNegation" in str(error.value)
    assert coordinates(molecule) == before


@pytest.mark.parametrize("query", [False, True])
def test_uff_missing_selected_geometry_has_canonical_source_cause(query):
    molecule = make_molecule()
    before = coordinates(molecule)
    with pytest.raises(ck.OperationError) as error:
        if query:
            molecule.uff_energy_gradient_with_params(ck.UffEvaluationParams(conformer_id=999))
        else:
            molecule.with_uff_optimized_with_params(ck.UffOptimizationParams(conformer_id=999))
    assert error.value.kind == "UffOptimization"
    assert isinstance(error.value.__cause__, ck.UffOptimizationError)
    assert error.value.__cause__.kind == "MissingConformer"
    assert error.value.__cause__.requested == 999
    assert coordinates(molecule) == before


def test_uff_parameter_missing_cache_has_original_typed_borrowed_cause():
    molecule = ck.Molecule.new().to_builder().build()
    assert molecule.uff_has_all_molecule_params() is True
    molecule = ck.Molecule.from_smiles("CC").with_hydrogens()
    with pytest.raises(ck.UffParameterQueryError) as error:
        molecule.uff_has_all_molecule_params()
    assert error.value.kind == "Cache"
    assert isinstance(error.value.__cause__, ck.OperationError)
    assert molecule.with_assigned_valence().uff_has_all_molecule_params() is True


def test_uff_default_and_configured_single_multi_protocols_equal():
    molecule = make_molecule()
    default = molecule.with_uff_optimized()
    explicit = molecule.with_uff_optimized_with_params(ck.UffOptimizationParams())
    assert default.energy() == explicit.energy()
    assert default.status_code() == explicit.status_code()
    assert coordinates(default.molecule()) == coordinates(explicit.molecule())
    default = molecule.with_uff_optimized_confs()
    explicit = molecule.with_uff_optimized_confs_with_params(ck.UffConformerOptimizationParams())
    assert [(r.status_code(),r.energy()) for r in default.conformer_results()] == [(r.status_code(),r.energy()) for r in explicit.conformer_results()]
    assert coordinates(default.molecule()) == coordinates(explicit.molecule())
