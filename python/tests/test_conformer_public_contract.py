"""Delivery proposal for native conformer projections, checked after native build."""
import json
import math
import numpy as np
import pytest
import cosmolkit as ck
from rdkit.Chem import rdDistGeom


def flags(is_3d):
    return ck.Coordinate3DInputParams(is_3d=is_3d)


def replacement(identifier):
    return ck.Replace3DCoordinatesParams(conformer_id=identifier)


def water():
    return ck.Molecule.from_smiles("O").with_hydrogens()


def seeded():
    params = ck.EmbedParams.etkdg_v3()
    params = configured(params, random_seed=42)
    params = configured(params, max_iterations=3)
    params = configured(params, num_threads=1)
    params = configured(params, track_failures=True)
    return params


def test_manual_coordinate_forms_ids_flags_and_atomic_errors():
    source = water()
    coords = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    appended = source.with_added_3d_conformer(coords)
    assert source.num_3d_conformers() == 0
    assert appended.num_3d_conformers() == 1
    assert appended.conformers_3d()[0].id() == 0
    assert source.add_3d_conformer_(coords) == 0
    assert source.add_3d_conformer_with_params_(coords, flags(False)) == 1
    assert not source.conformers_3d()[1].is_3d()
    changed = [[2.0, 3.0, 4.0]] * 3
    value = source.with_3d_coordinates_with_params(changed, replacement(1))
    assert source.conformers_3d()[1].coordinates() == coords
    assert value.conformers_3d()[1].coordinates() == changed
    assert not value.conformers_3d()[1].is_3d()
    source.set_3d_coordinates_with_params_(changed, replacement(1))
    snapshot = [(c.id(), c.coordinates(), c.is_3d()) for c in source.conformers_3d()]
    for coords, identifier in [([[0.0, 0.0, 0.0]], 1), ([[math.nan, 0.0, 0.0]] * 3, 1), ([[0.0, 0.0, 0.0]] * 3, 99)]:
        with pytest.raises(ValueError):
            source.set_3d_coordinates_with_params_(coords, replacement(identifier))
        assert [(c.id(), c.coordinates(), c.is_3d()) for c in source.conformers_3d()] == snapshot
    only = source.with_only_3d_conformer(changed)
    assert [c.id() for c in only.conformers_3d()] == [0]
    assert source.num_3d_conformers() == 2
    assert source.set_only_3d_conformer_(changed) == 0
    assert source.num_3d_conformers() == 1
    cleared = source.with_cleared_3d_conformers()
    assert cleared.num_3d_conformers() == 0
    assert source.num_3d_conformers() == 1
    source.clear_3d_conformers_()
    assert source.num_3d_conformers() == 0


def test_results_snapshots_all_generation_forms_and_atomic_typed_error():
    original = water()
    params = seeded()
    result = original.with_3d_conformer_result(params)
    assert result.ok() and result.conf_id() == 0
    assert result.molecule().num_3d_conformers() == 1
    assert original.num_3d_conformers() == 0
    assert result.params().random_seed == 42
    assert params.failures == []
    assert len(result.params().failures) == 12
    snapshot = result.params().to_json()
    params = configured(params, random_seed=123)
    assert result.params().to_json() == snapshot
    params = seeded()
    inplace = water()
    committed = inplace.embed_3d_conformer_result_(params)
    assert committed.conf_id() == 0 and committed.molecule().num_3d_conformers() == 1
    assert inplace.conformers_3d()[0].coordinates() == result.molecule().conformers_3d()[0].coordinates()
    params = configured(params, clear_confs=False)
    multiple = inplace.with_3d_conformers_result(3, params)
    assert multiple.conf_ids() == [1, 2, 3]
    assert multiple.requested_num_confs() == 3 and multiple.generated_count() == 3
    assert multiple.molecule().num_3d_conformers() == 4
    committed = inplace.embed_3d_conformers_result_(3, params)
    assert committed.conf_ids() == [1, 2, 3]
    assert inplace.num_3d_conformers() == 4
    assert params.failures == []
    assert len(committed.params().failures) == 12
    before = [c.coordinates() for c in inplace.conformers_3d()]
    params = configured(params, et_version=99)
    with pytest.raises(ValueError, match="Only version 1 and 2"):
        inplace.embed_3d_conformer_result_(params)
    assert [c.coordinates() for c in inplace.conformers_3d()] == before
    for value, inplace_name, count in [("with_3d_conformer", "embed_3d_conformer_", 1), ("with_3d_conformers", "embed_3d_conformers_", 2)]:
        arguments = [] if count == 1 else [count]
        assert getattr(original, value)(*arguments).num_3d_conformers() == count
        actual = water()
        assert getattr(actual, inplace_name)(*arguments) is None
        assert actual.num_3d_conformers() == count
        params = seeded()
        assert getattr(original, value)(*arguments, params=params).num_3d_conformers() == count
        assert getattr(actual, inplace_name)(*arguments, params=params) is None
        assert actual.num_3d_conformers() == count


@pytest.mark.parametrize("factory,source", [("etkdg", "ETKDG"), ("etkdg_v2", "ETKDGv2"), ("etkdg_v3", "ETKDGv3"), ("sr_etkdg_v3", "srETKDGv3")])
def test_preset_scalar_fields_and_defaults_match_native_reference(factory, source):
    actual = getattr(ck.EmbedParams, factory)()
    expected = getattr(rdDistGeom, source)()
    fields = {
        "max_iterations":"maxIterations", "num_threads":"numThreads", "random_seed":"randomSeed", "clear_confs":"clearConfs", "use_random_coords":"useRandomCoords", "box_size_mult":"boxSizeMult", "rand_neg_eig":"randNegEig", "num_zero_fail":"numZeroFail", "optimizer_force_tol":"optimizerForceTol", "ignore_smoothing_failures":"ignoreSmoothingFailures", "enforce_chirality":"enforceChirality", "use_exp_torsion_angle_prefs":"useExpTorsionAnglePrefs", "use_basic_knowledge":"useBasicKnowledge", "verbose":"verbose", "basin_thresh":"basinThresh", "prune_rms_thresh":"pruneRmsThresh", "only_heavy_atoms_for_rms":"onlyHeavyAtomsForRMS", "et_version":"ETversion", "embed_fragments_separately":"embedFragmentsSeparately", "use_small_ring_torsions":"useSmallRingTorsions", "use_macrocycle_torsions":"useMacrocycleTorsions", "use_macrocycle14config":"useMacrocycle14config", "timeout":"timeout", "force_trans_amides":"forceTransAmides", "use_symmetry_for_pruning":"useSymmetryForPruning", "bounds_mat_force_scaling":"boundsMatForceScaling", "track_failures":"trackFailures", "enable_sequential_random_seeds":"enableSequentialRandomSeeds", "symmetrize_conjugated_terminal_groups_for_pruning":"symmetrizeConjugatedTerminalGroupsForPruning",
    }
    for field, native_field in fields.items():
        assert getattr(actual, field) == getattr(expected, native_field), field
        with pytest.raises(AttributeError):
            setattr(actual, field, getattr(actual, field))
    assert actual.coord_map is None and actual.cpci is None and actual.failures == []
    with pytest.raises(AttributeError):
        actual.failures = []


def test_params_eight_factories_json_maps_and_roundtrip():
    for factory in ["dg", "kdg", "etdg", "etdg_v2", "etkdg", "etkdg_v2", "etkdg_v3", "sr_etkdg_v3"]:
        params = getattr(ck.EmbedParams, factory)()
        assert isinstance(params, ck.EmbedParams)
        assert isinstance(repr(params), str)
        restored = params.with_json(params.to_json())
        assert json.loads(restored.to_json()) == json.loads(params.to_json())
    source = ck.EmbedParams()
    changed = source.with_json('{"randomSeed":42,"clearConfs":false,"maxIterations":3}')
    assert source.random_seed == -1 and source.clear_confs
    assert changed.random_seed == 42 and not changed.clear_confs and changed.max_iterations == 3
    changed = configured(changed, coord_map={0: [1.25, 2.5, 3.75]})
    changed = configured(changed, cpci={(0, 1): -0.5})
    restored = changed.with_json(changed.to_json())
    assert restored.coord_map == changed.coord_map and restored.cpci == changed.cpci
    changed = configured(changed, coord_map=None)
    changed = configured(changed, cpci=None)
    assert changed.coord_map is None and changed.cpci is None
    with pytest.raises(ValueError):
        source.with_json("{")


def test_bounds_array_is_numeric_readonly_query():
    source = water()
    actual = source.dg_bounds_matrix()
    assert isinstance(actual, np.ndarray) and actual.shape == (3, 3)
    assert actual.dtype == np.dtype("float64")
    assert np.isfinite(actual).all() and np.diag(actual).tolist() == [0.0, 0.0, 0.0]
    assert source.num_3d_conformers() == 0


def configured(params, **changes):
    values = {field: getattr(params, field) for field in ['max_iterations', 'num_threads', 'random_seed', 'clear_confs', 'use_random_coords', 'box_size_mult', 'rand_neg_eig', 'num_zero_fail', 'coord_map', 'optimizer_force_tol', 'ignore_smoothing_failures', 'enforce_chirality', 'use_exp_torsion_angle_prefs', 'use_basic_knowledge', 'verbose', 'basin_thresh', 'prune_rms_thresh', 'only_heavy_atoms_for_rms', 'et_version', 'embed_fragments_separately', 'use_small_ring_torsions', 'use_macrocycle_torsions', 'use_macrocycle14config', 'timeout', 'cpci', 'force_trans_amides', 'use_symmetry_for_pruning', 'bounds_mat_force_scaling', 'track_failures', 'enable_sequential_random_seeds', 'symmetrize_conjugated_terminal_groups_for_pruning']}
    values.update(changes)
    return ck.EmbedParams(**values)
