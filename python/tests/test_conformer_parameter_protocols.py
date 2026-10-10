"""Writable configuration fields and independent result snapshots."""
import cosmolkit as ck
import pytest


def test_all_parameter_fields_construct_distinct_values_and_preserve_source_snapshot():
    params = ck.EmbedParams()
    snapshot = params.to_json()
    scalar = {
        "max_iterations": 19, "num_threads": 4, "random_seed": 123,
        "clear_conformers": False, "use_random_coords": True, "box_size_mult": 3.25,
        "rand_neg_eig": False, "num_zero_fail": 2, "optimizer_force_tol": 0.125,
        "ignore_smoothing_failures": True, "enforce_chirality": False,
        "use_exp_torsion_angle_prefs": True, "use_basic_knowledge": True,
        "verbose": True, "basin_thresh": 6.5, "prune_rms_thresh": 0.25,
        "only_heavy_atoms_for_rms": False, "et_version": 1,
        "embed_fragments_separately": False, "use_small_ring_torsions": True,
        "use_macrocycle_torsions": True, "use_macrocycle14config": True,
        "timeout": 7, "force_trans_amides": False, "use_symmetry_for_pruning": False,
        "bounds_mat_force_scaling": 2.5, "track_failures": True,
        "enable_sequential_random_seeds": True,
        "symmetrize_conjugated_terminal_groups_for_pruning": False,
    }
    for field, value in scalar.items():
        # et_version's default may already be1. Its separate v3 state is2.
        if field == "et_version":
            params = configured(params, et_version=2)
        assert getattr(params, field) != value, field
        previous = params
        params = configured(params, **{field: value})
        assert getattr(previous, field) != value
        assert getattr(params, field) == value, field
    params = configured(params, coord_map={0: [1.0, 2.0, 3.0]})
    params = configured(params, cpci={(0, 1): 0.125})
    mapping = params.coord_map
    mapping[0] = [9.0, 9.0, 9.0]
    assert params.coord_map == {0: [9.0, 9.0, 9.0]}
    cpci = params.cpci
    cpci[(0, 1)] = -8.0
    assert params.cpci == {(0, 1): -8.0}
    assert ck.EmbedParams().to_json() == snapshot
    assert repr(params) == "EmbedParams(random_seed=123, num_threads=4, prune_rms_thresh=0.25, clear_conformers=false)"
    assert params.failures == []
    with pytest.raises(AttributeError):
        params.failures = [1]


def test_result_types_reprs_and_returned_molecule_params_are_detached_snapshots():
    source = ck.Molecule.from_smiles("O").with_hydrogens()
    params = ck.EmbedParams.etkdg_v3()
    params = configured(params, random_seed=42)
    single = source.with_3d_conformer_result(params)
    multiple = source.with_3d_conformers_result(2, params)
    assert isinstance(single, ck.EmbedMoleculeResult)
    assert isinstance(multiple, ck.EmbedMultipleConfsResult)
    assert repr(single) == "EmbedMoleculeResult(conf_id=0, ok=true)"
    assert repr(multiple) == "EmbedMultipleConfsResult(requested_num_confs=2, generated_count=2)"
    assert source.num_3d_conformers() == 0
    for result, count in [(single, 1), (multiple, 2)]:
        returned = result.molecule()
        returned.clear_3d_conformers_()
        assert result.molecule().num_3d_conformers() == count
        copied = result.params()
        copied = configured(copied, random_seed=987)
        assert copied.random_seed == 987
        assert result.params().random_seed == 42
    ids = multiple.conf_ids()
    ids.append(99)
    assert multiple.conf_ids() == [0, 1]


def configured(params, **changes):
    values = {field: getattr(params, field) for field in ['max_iterations', 'num_threads', 'random_seed', 'clear_conformers', 'use_random_coords', 'box_size_mult', 'rand_neg_eig', 'num_zero_fail', 'coord_map', 'optimizer_force_tol', 'ignore_smoothing_failures', 'enforce_chirality', 'use_exp_torsion_angle_prefs', 'use_basic_knowledge', 'verbose', 'basin_thresh', 'prune_rms_thresh', 'only_heavy_atoms_for_rms', 'et_version', 'embed_fragments_separately', 'use_small_ring_torsions', 'use_macrocycle_torsions', 'use_macrocycle14config', 'timeout', 'cpci', 'force_trans_amides', 'use_symmetry_for_pruning', 'bounds_mat_force_scaling', 'track_failures', 'enable_sequential_random_seeds', 'symmetrize_conjugated_terminal_groups_for_pruning']}
    values.update(changes)
    return ck.EmbedParams(**values)
