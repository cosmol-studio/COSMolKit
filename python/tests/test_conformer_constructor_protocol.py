"""Full immutable constructor/value protocol proposal under ROOT decision001."""
import inspect
import json
import cosmolkit as ck
import pytest

DEFAULTS = {
    'max_iterations':0, 'num_threads':1, 'random_seed':-1,
    'clear_confs':True, 'use_random_coords':False, 'box_size_mult':2.0,
    'rand_neg_eig':True, 'num_zero_fail':1, 'coord_map':None,
    'optimizer_force_tol':1e-3, 'ignore_smoothing_failures':False,
    'enforce_chirality':True, 'use_exp_torsion_angle_prefs':False,
    'use_basic_knowledge':False, 'verbose':False, 'basin_thresh':5.0,
    'prune_rms_thresh':-1.0, 'only_heavy_atoms_for_rms':True,
    'et_version':2, 'embed_fragments_separately':True,
    'use_small_ring_torsions':False, 'use_macrocycle_torsions':False,
    'use_macrocycle14config':False, 'timeout':0, 'cpci':None,
    'force_trans_amides':True, 'use_symmetry_for_pruning':True,
    'bounds_mat_force_scaling':1.0, 'track_failures':False,
    'enable_sequential_random_seeds':False,
    'symmetrize_conjugated_terminal_groups_for_pruning':True,
}
DISTINCT = {
    'max_iterations':19, 'num_threads':4, 'random_seed':123,
    'clear_confs':False, 'use_random_coords':True, 'box_size_mult':3.25,
    'rand_neg_eig':False, 'num_zero_fail':2, 'coord_map':{0:[1.,2.,3.]},
    'optimizer_force_tol':0.125, 'ignore_smoothing_failures':True,
    'enforce_chirality':False, 'use_exp_torsion_angle_prefs':True,
    'use_basic_knowledge':True, 'verbose':True, 'basin_thresh':6.5,
    'prune_rms_thresh':0.25, 'only_heavy_atoms_for_rms':False,
    'et_version':1, 'embed_fragments_separately':False,
    'use_small_ring_torsions':True, 'use_macrocycle_torsions':True,
    'use_macrocycle14config':True, 'timeout':7, 'cpci':{(0,1):0.125},
    'force_trans_amides':False, 'use_symmetry_for_pruning':False,
    'bounds_mat_force_scaling':2.5, 'track_failures':True,
    'enable_sequential_random_seeds':True,
    'symmetrize_conjugated_terminal_groups_for_pruning':False,
}


def test_constructor_declares_all_fields_explicit_source_defaults():
    signature = inspect.signature(ck.EmbedParams)
    assert list(signature.parameters) == list(DEFAULTS)
    params = ck.EmbedParams()
    assert params.failures == []
    for field, expected in DEFAULTS.items():
        assert getattr(params, field) == expected, field
        assert signature.parameters[field].default == expected, field
        assert signature.parameters[field].kind == inspect.Parameter.KEYWORD_ONLY
    with pytest.raises(TypeError):
        ck.EmbedParams(failures=[1])


@pytest.mark.parametrize('field', DEFAULTS)
def test_each_constructor_field_preserves_other_defaults(field):
    params = ck.EmbedParams(**{field:DISTINCT[field]})
    assert getattr(params, field) == DISTINCT[field]
    for other, expected in DEFAULTS.items():
        if other != field:
            assert getattr(params, other) == expected, other


def test_constructor_maps_and_value_replacement_are_detached():
    maps={0:[1.,2.,3.]};cpci={(0,1):0.125}
    params=ck.EmbedParams(coord_map=maps,cpci=cpci)
    maps[0][0]=9.;cpci[(0,1)]=9.
    assert params.coord_map == {0:[1.,2.,3.]}
    assert params.cpci == {(0,1):0.125}
    changed=params.with_json('{"randomSeed":42}')
    assert params.random_seed == -1 and changed.random_seed == 42
    assert changed.coord_map == params.coord_map and changed.cpci == params.cpci
    restored=changed.with_json(changed.to_json())
    assert json.loads(restored.to_json()) == json.loads(changed.to_json())
    returned=restored.coord_map;returned[0][0]=99.
    assert restored.coord_map == {0:[1.,2.,3.]}
    failures=restored.failures;failures.append(1)
    assert restored.failures == []
