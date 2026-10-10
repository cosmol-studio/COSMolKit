"""Writable input configuration and immutable generation result boundaries."""
import pytest
import cosmolkit as ck

FIELDS = [('max_iterations', 19), ('num_threads', 4), ('random_seed', 123), ('clear_conformers', False), ('use_random_coords', True), ('box_size_mult', 3.25), ('rand_neg_eig', False), ('num_zero_fail', 2), ('optimizer_force_tol', 0.125), ('ignore_smoothing_failures', True), ('enforce_chirality', False), ('use_exp_torsion_angle_prefs', True), ('use_basic_knowledge', True), ('verbose', True), ('basin_thresh', 6.5), ('prune_rms_thresh', 0.25), ('only_heavy_atoms_for_rms', False), ('et_version', 1), ('embed_fragments_separately', False), ('use_small_ring_torsions', True), ('use_macrocycle_torsions', True), ('use_macrocycle14config', True), ('timeout', 7), ('force_trans_amides', False), ('use_symmetry_for_pruning', False), ('bounds_mat_force_scaling', 2.5), ('track_failures', True), ('enable_sequential_random_seeds', True), ('symmetrize_conjugated_terminal_groups_for_pruning', False), ('coord_map', {0: [1.0, 2.0, 3.0]}), ('cpci', {(0, 1): 0.125}), ('failures', [1])]

@pytest.mark.parametrize("field,value", [(field, value) for field, value in FIELDS if field != "failures"])
def test_embed_configuration_assignment_matches_constructor_input(field, value):
    params = ck.EmbedParams()
    setattr(params, field, value)
    expected = ck.EmbedParams(**{field: value})
    assert getattr(params, field) == getattr(expected, field)
    assert params.to_json() == expected.to_json()
    assert list(params.failures) == []


def test_computed_failures_remain_readonly():
    params = ck.EmbedParams()
    before = (params.to_json(), list(params.failures))
    with pytest.raises(AttributeError):
        params.failures = [1]
    assert (params.to_json(), list(params.failures)) == before


@pytest.mark.parametrize("method,count", [("with_3d_conformer",None),("embed_3d_conformer_",None),("with_3d_conformer_result",None),("embed_3d_conformer_result_",None),("with_3d_conformers",2),("embed_3d_conformers_",2),("with_3d_conformers_result",2),("embed_3d_conformers_result_",2)])
def test_all_generation_forms_borrow_configuration_without_failure_progress_writeback(method,count):
    mol = ck.Molecule.from_smiles("O").with_hydrogens()
    params = ck.EmbedParams.etkdg_v3().with_json('{"randomSeed":42,"maxIterations":3,"numThreads":1,"trackFailures":true}')
    before = (params.to_json(), list(params.failures))
    result = getattr(mol, method)(*([] if count is None else [count]), params=params)
    assert (params.to_json(), list(params.failures)) == before
    if "result" in method:
        # RDKit 2026.03.6 Embedder.h: EmbedFailureCauses::END_OF_ENUM.
        assert len(result.params().failures) == 15
        assert result.params().random_seed == 42
