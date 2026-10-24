//! Representation-only embedding keyword configuration.
use crate::embedding_values::EmbedParams;
use cosmolkit as ck;
use std::collections::BTreeMap;
pub(crate) fn keywords(
    max_iterations: Option<u32>,
    num_threads: Option<i32>,
    random_seed: Option<i32>,
    clear_conformers: Option<bool>,
    use_random_coords: Option<bool>,
    box_size_mult: Option<f64>,
    rand_neg_eig: Option<bool>,
    num_zero_fail: Option<u32>,
    coord_map: Option<BTreeMap<i32, [f64; 3]>>,
    optimizer_force_tol: Option<f64>,
    ignore_smoothing_failures: Option<bool>,
    enforce_chirality: Option<bool>,
    use_exp_torsion_angle_prefs: Option<bool>,
    use_basic_knowledge: Option<bool>,
    verbose: Option<bool>,
    basin_thresh: Option<f64>,
    prune_rms_thresh: Option<f64>,
    only_heavy_atoms_for_rms: Option<bool>,
    et_version: Option<u32>,
    embed_fragments_separately: Option<bool>,
    use_small_ring_torsions: Option<bool>,
    use_macrocycle_torsions: Option<bool>,
    use_macrocycle14config: Option<bool>,
    timeout: Option<u32>,
    cpci: Option<BTreeMap<(u32, u32), f64>>,
    force_trans_amides: Option<bool>,
    use_symmetry_for_pruning: Option<bool>,
    bounds_mat_force_scaling: Option<f64>,
    track_failures: Option<bool>,
    enable_sequential_random_seeds: Option<bool>,
    symmetrize_conjugated_terminal_groups_for_pruning: Option<bool>,
) -> ck::EmbedParams {
    // All eight registered embedding operations declare etkdg_v3 as their
    // default_args; partial keywords extend that operation default, not the
    // distinct generic EmbedParams::new constructor.
    let mut inner = ck::EmbedParams::etkdg_v3();
    if let Some(value) = max_iterations {
        inner.max_iterations = value;
    }
    if let Some(value) = num_threads {
        inner.num_threads = value;
    }
    if let Some(value) = random_seed {
        inner.random_seed = value;
    }
    if let Some(value) = clear_conformers {
        inner.clear_conformers = value;
    }
    if let Some(value) = use_random_coords {
        inner.use_random_coords = value;
    }
    if let Some(value) = box_size_mult {
        inner.box_size_mult = value;
    }
    if let Some(value) = rand_neg_eig {
        inner.rand_neg_eig = value;
    }
    if let Some(value) = num_zero_fail {
        inner.num_zero_fail = value;
    }
    if let Some(value) = coord_map {
        inner.coord_map = Some(value);
    }
    if let Some(value) = optimizer_force_tol {
        inner.optimizer_force_tol = value;
    }
    if let Some(value) = ignore_smoothing_failures {
        inner.ignore_smoothing_failures = value;
    }
    if let Some(value) = enforce_chirality {
        inner.enforce_chirality = value;
    }
    if let Some(value) = use_exp_torsion_angle_prefs {
        inner.use_exp_torsion_angle_prefs = value;
    }
    if let Some(value) = use_basic_knowledge {
        inner.use_basic_knowledge = value;
    }
    if let Some(value) = verbose {
        inner.verbose = value;
    }
    if let Some(value) = basin_thresh {
        inner.basin_thresh = value;
    }
    if let Some(value) = prune_rms_thresh {
        inner.prune_rms_thresh = value;
    }
    if let Some(value) = only_heavy_atoms_for_rms {
        inner.only_heavy_atoms_for_rms = value;
    }
    if let Some(value) = et_version {
        inner.et_version = value;
    }
    if let Some(value) = embed_fragments_separately {
        inner.embed_fragments_separately = value;
    }
    if let Some(value) = use_small_ring_torsions {
        inner.use_small_ring_torsions = value;
    }
    if let Some(value) = use_macrocycle_torsions {
        inner.use_macrocycle_torsions = value;
    }
    if let Some(value) = use_macrocycle14config {
        inner.use_macrocycle14config = value;
    }
    if let Some(value) = timeout {
        inner.timeout = value;
    }
    if let Some(value) = cpci {
        inner.cpci = Some(value);
    }
    if let Some(value) = force_trans_amides {
        inner.force_trans_amides = value;
    }
    if let Some(value) = use_symmetry_for_pruning {
        inner.use_symmetry_for_pruning = value;
    }
    if let Some(value) = bounds_mat_force_scaling {
        inner.bounds_mat_force_scaling = value;
    }
    if let Some(value) = track_failures {
        inner.track_failures = value;
    }
    if let Some(value) = enable_sequential_random_seeds {
        inner.enable_sequential_random_seeds = value;
    }
    if let Some(value) = symmetrize_conjugated_terminal_groups_for_pruning {
        inner.symmetrize_conjugated_terminal_groups_for_pruning = value;
    }
    inner
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::native::Molecule;
    use crate::operation_errors::OperationError;
    use pyo3::ffi::c_str;
    use pyo3::prelude::*;
    use pyo3::types::{PyDict, PyModule};

    #[test]
    fn all_embedding_value_inplace_report_forms_match_native_coordinates_and_failure_atomicity() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native=ck::Molecule::from_smiles("CCO").unwrap();
            let before=native.to_binary().unwrap();
            let mut params=ck::EmbedParams::etkdg_v3();
            params.random_seed=17;
            params.track_failures=true;
            let one=native.with_3d_conformer_result_with_params(&params).unwrap();
            let two=native.with_3d_conformers_result_with_params(2,&params).unwrap();
            let module=PyModule::new(py,"embedding_call_test")?;
            module.add_class::<Molecule>()?;
            module.add_class::<OperationError>()?;
            let locals=PyDict::new(py);
            locals.set_item("api",module)?;
            locals.set_item("params",Py::new(py,EmbedParams::from_rust(py, params)?)?)?;
            locals.set_item("before",before)?;
            locals.set_item("expected_one",one.molecule().to_binary().unwrap())?;
            locals.set_item("expected_two",two.molecule().to_binary().unwrap())?;
            locals.set_item("one_id",one.conf_id())?;
            locals.set_item("two_ids",two.conf_ids())?;
            locals.set_item("one_failures",one.params().failures().clone())?;
            locals.set_item("two_failures",two.params().failures().clone())?;
            py.run(c_str!(r#"
def fresh(): return api.Molecule.from_smiles('CCO')
for name in ['with_3d_conformer', 'embed_3d_conformer_',
             'with_3d_conformer_result', 'embed_3d_conformer_result_',
             'with_3d_conformers', 'embed_3d_conformers_',
             'with_3d_conformers_result', 'embed_3d_conformers_result_']:
    multiple = 'conformers' in name
    inplace = name.endswith('_')
    report = 'result' in name
    args = (2,) if multiple else ()
    expected = expected_two if multiple else expected_one
    alias = name[:-1] + '_with_params_' if inplace else name + '_with_params'
    for form in ['params', 'keywords', 'alias']:
        source = fresh()
        if form == 'params': output = getattr(source, name)(*args, params)
        elif form == 'keywords': output = getattr(source, name)(*args, random_seed=17, track_failures=True)
        else: output = getattr(source, alias)(*args, params)
        if inplace:
            produced = source
            if not report: assert output is None
        else:
            assert bytes(source.to_binary()) == bytes(before)
            produced = output.molecule() if report else output
        assert bytes(produced.to_binary()) == bytes(expected), (name, form)
        if report:
            assert bytes(output.molecule().to_binary()) == bytes(expected)
            assert output.params().random_seed == 17 and output.params().track_failures
            assert output.params().failures == (two_failures if multiple else one_failures)
            if multiple:
                assert output.conf_ids() == two_ids
                assert output.requested_num_confs() == 2
                assert output.generated_count() == len(two_ids)
            else: assert output.conf_id() == one_id and output.ok() == (one_id >= 0)
    source = fresh()
    output = getattr(source, name)(*args)  # real native default branch
    produced = source if inplace else (output.molecule() if report else output)
    assert produced.num_3d_conformers() == (2 if multiple else 1)
    if report:
        assert output.params().use_exp_torsion_angle_prefs
        assert output.params().use_basic_knowledge
        assert output.params().use_macrocycle_torsions
    try: getattr(source, name)(*args, params, random_seed=17)
    except TypeError: pass
    else: raise AssertionError('mixed embedding configuration accepted')
source = fresh()
try: source.embed_3d_conformer_(coord_map={999: [0.0, 0.0, 0.0]})
except api.OperationError as error:
    assert error.conformer is not None
    assert error.cause is error.__cause__
else: raise AssertionError('invalid native coordinate map accepted')
assert bytes(source.to_binary()) == bytes(before)
assert params.failures == []
"#),Some(&locals),None)?;
            Ok(())
        }).unwrap();
    }
}
