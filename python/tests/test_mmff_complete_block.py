import json
from pathlib import Path
import pytest
import cosmolkit as ck


def make_molecule(smiles="CC", coords=((0.0, 0.0, 0.0), (2.0, 0.0, 0.0))):
    base = ck.Molecule.from_smiles(smiles)
    builder = base.to_builder()
    builder.add_3d_conformer(list(coords))
    return builder.build()


def coordinates(molecule):
    return [row.coordinates() for row in molecule.conformers_3d()]


def test_mmff_source_default_parameters():
    single = ck.MmffOptimizationParams()
    assert (single.mmff_variant, single.max_iterations, single.non_bonded_threshold, single.conformer_id, single.ignore_interfragment_interactions) == ("MMFF94", 200, 100.0, None, True)
    multi = ck.MmffConformerOptimizationParams()
    assert (multi.num_threads, multi.max_iterations, multi.mmff_variant, multi.non_bonded_threshold, multi.ignore_interfragment_interactions) == (1, 1000, "MMFF94", 10.0, True)


def test_mmff_value_results_and_original_protocols():
    molecule = make_molecule()
    original = coordinates(molecule)
    result = molecule.with_mmff_optimized()
    assert isinstance(result, ck.MmffOptimizeMoleculeResult)
    assert result.status_code() == 0
    assert result.needs_more() is False
    assert result.molecule().num_atoms() == molecule.num_atoms()
    assert "MmffOptimizeMoleculeResult" in repr(result)
    assert coordinates(molecule) == original
    multi = molecule.with_mmff_optimized_confs()
    assert isinstance(multi, ck.MmffOptimizeMoleculeConfsResult)
    assert len(multi.conformer_results()) == 1
    assert "MmffOptimizeMoleculeConfsResult" in repr(multi)
    row = multi.conformer_results()[0]
    assert isinstance(row, ck.MmffOptimizeMoleculeConfResult)
    assert row.status_code() == 0 and row.needs_more() is False
    assert isinstance(row.energy(), float)
    assert "MmffOptimizeMoleculeConfResult" in repr(row)
    assert multi.molecule().num_atoms() == molecule.num_atoms()
    assert coordinates(molecule) == original


def test_mmff_zero_iterations_and_configured_calls():
    molecule = make_molecule()
    original = coordinates(molecule)
    result = molecule.with_mmff_optimized_with_params(ck.MmffOptimizationParams(max_iterations=0))
    configured = molecule.with_mmff_optimized_with_params(ck.MmffOptimizationParams(max_iterations=0))
    assert result.status_code() == configured.status_code() == 1
    assert result.needs_more() and configured.needs_more()
    assert coordinates(result.molecule()) == coordinates(configured.molecule()) == original
    multi = molecule.with_mmff_optimized_confs_with_params(ck.MmffConformerOptimizationParams(max_iterations=0))
    configured_multi = molecule.with_mmff_optimized_confs_with_params(ck.MmffConformerOptimizationParams(max_iterations=0))
    assert multi.conformer_results()[0].status_code() == 1
    assert configured_multi.conformer_results()[0].status_code() == 1
    assert multi.conformer_results()[0].energy() == configured_multi.conformer_results()[0].energy()
    assert coordinates(molecule) == original


def test_mmff_source_invalid_typing_sentinels():
    molecule = make_molecule("[He]", ((0.0, 0.0, 0.0),))
    assert molecule.mmff_has_all_molecule_params() is False
    single = molecule.with_mmff_optimized()
    assert single.status_code() == -1 and single.needs_more() is False
    multi = molecule.with_mmff_optimized_confs()
    row = multi.conformer_results()[0]
    assert row.status_code() == -1 and row.energy() == -1.0
    assert coordinates(single.molecule()) == coordinates(multi.molecule()) == coordinates(molecule)
    assert molecule.mmff_energy_gradient() is None


@pytest.mark.parametrize("selector", [-1, -2, -3, -(2**31)])
def test_mmff_source_negative_selectors_choose_first(selector):
    molecule = make_molecule()
    before = coordinates(molecule)
    params = ck.MmffOptimizationParams(max_iterations=0, conformer_id=selector)
    assert params.conformer_id is None
    selected = molecule.with_mmff_optimized_with_params(params)
    default = molecule.with_mmff_optimized_with_params(ck.MmffOptimizationParams(max_iterations=0))
    assert selected.status_code() == default.status_code() == 1
    assert coordinates(selected.molecule()) == coordinates(default.molecule()) == before
    evaluated = molecule.mmff_energy_gradient_with_params(ck.MmffEvaluationParams(conformer_id=selector))
    original = molecule.mmff_energy_gradient()
    assert evaluated is not None and original is not None
    assert evaluated.energy() == original.energy()
    assert evaluated.gradient() == original.gradient()


@pytest.mark.parametrize("selector", [-2**31 - 1, 2**31])
@pytest.mark.parametrize("params_type", [ck.MmffOptimizationParams, ck.MmffEvaluationParams])
def test_mmff_selector_source_signed_int32_boundary(selector, params_type):
    with pytest.raises(OverflowError):
        params_type(conformer_id=selector)


def test_mmff_original_empty_conformer_negative_fixture():
    builder = ck.Molecule.new().to_builder()
    builder.add_3d_conformer([])
    molecule = builder.build()
    failures = []
    for selector in [-2, None]:
        with pytest.raises(ck.OperationError) as error:
            molecule.with_mmff_optimized_with_params(ck.MmffOptimizationParams(max_iterations=25,conformer_id=selector))
        failures.append(error.value)
        assert isinstance(error.value.__cause__, ck.MmffOptimizationError)
        assert "NoPoints" in str(error.value)
    assert str(failures[0]) == str(failures[1])
    assert coordinates(molecule) == [[]]


def test_mmff_absent_selector_is_typed_and_input_is_preserved():
    molecule = make_molecule()
    before = coordinates(molecule)
    with pytest.raises(ck.OperationError) as error:
        molecule.with_mmff_optimized_with_params(ck.MmffOptimizationParams(conformer_id=999))
    assert error.value.kind == "MmffOptimization"
    assert isinstance(error.value.__cause__, ck.MmffOptimizationError)
    assert error.value.__cause__.domain == "mmff_optimization"
    assert error.value.__cause__.__cause__ is not None
    assert coordinates(molecule) == before


@pytest.mark.parametrize("variant, expected", [("MMFF94s", ck.MmffVariant.Mmff94s), ("MMFF94S", ck.MmffVariant.Mmff94), ("invalid", ck.MmffVariant.Mmff94)])
def test_mmff_variant_parser_matches_original_source(variant, expected):
    molecule = make_molecule()
    props = molecule.mmff_properties_with_params(ck.MmffPropertiesParams(variant))
    assert props.variant() == expected
    result = molecule.with_mmff_optimized_with_params(ck.MmffOptimizationParams(mmff_variant=variant,max_iterations=0))
    assert result.status_code() == 1


def test_mmff_property_lookup_bounds_and_canonical_types():
    molecule = ck.Molecule.from_smiles("CCO")
    props = molecule.mmff_properties()
    assert props.is_valid()
    assert [props.atom_type(i) for i in range(3)] == [1, 1, 6]
    atoms = props.atoms()
    assert len(atoms) == 3
    assert [atom.atom_type() for atom in atoms] == [1, 1, 6]
    with pytest.raises(ck.MmffMolPropertiesError) as error:
        props.atom_type(3)
    assert error.value.kind == "AtomIndexOutOfRange"
    assert molecule.mmff_has_all_molecule_params() is True
    assert ck.Molecule.new().mmff_has_all_molecule_params() is True


def test_mmff_thread_dispatch_order_matches_serial():
    base = make_molecule()
    builder = base.to_builder()
    builder.add_3d_conformer([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
    molecule = builder.build()
    before = coordinates(molecule)
    serial = molecule.with_mmff_optimized_confs_with_params(ck.MmffConformerOptimizationParams(num_threads=1,max_iterations=25))
    threaded = molecule.with_mmff_optimized_confs_with_params(ck.MmffConformerOptimizationParams(num_threads=2,max_iterations=25))
    assert [(r.status_code(), r.energy()) for r in serial.conformer_results()] == [(r.status_code(), r.energy()) for r in threaded.conformer_results()]
    assert coordinates(serial.molecule()) == coordinates(threaded.molecule())
    assert [r.id() for r in serial.molecule().conformers_3d()] == [0, 1]
    assert coordinates(molecule) == before


def test_mmff_parameter_values_are_immutable_and_short_methods_are_default_only():
    params = ck.MmffOptimizationParams()
    with pytest.raises(AttributeError):
        setattr(params, "max_iterations", 0)
    molecule = make_molecule()
    with pytest.raises(TypeError):
        getattr(molecule, "with_mmff_optimized")(max_iters=0)
    with pytest.raises(TypeError):
        getattr(molecule, "with_mmff_optimized_confs")(num_threads=2)
