//! Native force-field evaluation and optimization projection only.
use cosmolkit as ck;
use pyo3::prelude::*;
pub(crate) fn uff_energy_gradient_params(
    vdw_threshold: Option<f64>,
    conformer_id: Option<usize>,
    ignore_interfragment_interactions: Option<bool>,
) -> crate::params::UffEvaluationParams {
    let mut inner = ck::UffEvaluationParams::default();
    if let Some(value) = vdw_threshold {
        inner.vdw_threshold = value;
    }
    if let Some(value) = conformer_id {
        inner.conformer_id = Some(value);
    }
    if let Some(value) = ignore_interfragment_interactions {
        inner.ignore_interfragment_interactions = value;
    }
    crate::params::UffEvaluationParams { inner }
}
pub(crate) fn mmff_energy_gradient_params(
    mmff_variant: Option<String>,
    non_bonded_threshold: Option<f64>,
    conformer_id: Option<usize>,
    ignore_interfragment_interactions: Option<bool>,
) -> crate::params::MmffEvaluationParams {
    let mut inner = ck::MmffEvaluationParams::default();
    if let Some(value) = mmff_variant {
        inner.mmff_variant = value;
    }
    if let Some(value) = non_bonded_threshold {
        inner.non_bonded_threshold = value;
    }
    if let Some(value) = conformer_id {
        inner.conformer_id = Some(value);
    }
    if let Some(value) = ignore_interfragment_interactions {
        inner.ignore_interfragment_interactions = value;
    }
    crate::params::MmffEvaluationParams { inner }
}
pub(crate) fn with_uff_optimized_params(
    max_iterations: Option<i32>,
    vdw_threshold: Option<f64>,
    conformer_id: Option<usize>,
    ignore_interfragment_interactions: Option<bool>,
) -> crate::params::UffOptimizationParams {
    let mut inner = ck::UffOptimizationParams::default();
    if let Some(value) = max_iterations {
        inner.max_iterations = value;
    }
    if let Some(value) = vdw_threshold {
        inner.vdw_threshold = value;
    }
    if let Some(value) = conformer_id {
        inner.conformer_id = Some(value);
    }
    if let Some(value) = ignore_interfragment_interactions {
        inner.ignore_interfragment_interactions = value;
    }
    crate::params::UffOptimizationParams { inner }
}
pub(crate) fn with_uff_optimized_conformers_params(
    num_threads: Option<i32>,
    max_iterations: Option<i32>,
    vdw_threshold: Option<f64>,
    ignore_interfragment_interactions: Option<bool>,
) -> crate::params::UffConformerOptimizationParams {
    let mut inner = ck::UffConformerOptimizationParams::default();
    if let Some(value) = num_threads {
        inner.num_threads = value;
    }
    if let Some(value) = max_iterations {
        inner.max_iterations = value;
    }
    if let Some(value) = vdw_threshold {
        inner.vdw_threshold = value;
    }
    if let Some(value) = ignore_interfragment_interactions {
        inner.ignore_interfragment_interactions = value;
    }
    crate::params::UffConformerOptimizationParams { inner }
}
pub(crate) fn with_mmff_optimized_params(
    mmff_variant: Option<String>,
    max_iterations: Option<i32>,
    non_bonded_threshold: Option<f64>,
    conformer_id: Option<usize>,
    ignore_interfragment_interactions: Option<bool>,
) -> crate::params::MmffOptimizationParams {
    let mut inner = ck::MmffOptimizationParams::default();
    if let Some(value) = mmff_variant {
        inner.mmff_variant = value;
    }
    if let Some(value) = max_iterations {
        inner.max_iterations = value;
    }
    if let Some(value) = non_bonded_threshold {
        inner.non_bonded_threshold = value;
    }
    if let Some(value) = conformer_id {
        inner.conformer_id = Some(value);
    }
    if let Some(value) = ignore_interfragment_interactions {
        inner.ignore_interfragment_interactions = value;
    }
    crate::params::MmffOptimizationParams { inner }
}
pub(crate) fn with_mmff_optimized_conformers_params(
    num_threads: Option<i32>,
    max_iterations: Option<i32>,
    mmff_variant: Option<String>,
    non_bonded_threshold: Option<f64>,
    ignore_interfragment_interactions: Option<bool>,
) -> crate::params::MmffConformerOptimizationParams {
    let mut inner = ck::MmffConformerOptimizationParams::default();
    if let Some(value) = num_threads {
        inner.num_threads = value;
    }
    if let Some(value) = max_iterations {
        inner.max_iterations = value;
    }
    if let Some(value) = mmff_variant {
        inner.mmff_variant = value;
    }
    if let Some(value) = non_bonded_threshold {
        inner.non_bonded_threshold = value;
    }
    if let Some(value) = ignore_interfragment_interactions {
        inner.ignore_interfragment_interactions = value;
    }
    crate::params::MmffConformerOptimizationParams { inner }
}

use crate::error_framework::{Source, error_projection};
use crate::operation_errors::OperationError;
use crate::retained_error::RetainedError;
use crate::uff_parameter_errors::UffParameterError;
use pyo3::exceptions::PyValueError;
use std::error::Error;

error_projection!(retained UffParameterQueryError,UffParameterQueryErrorKind,ck::UffParameterQueryError,"forcefield",{
Cache(_),Parameters(_),
},cause=|inner:&RetainedError<ck::UffParameterQueryError>|{
 let node=inner.node().source().expect("both native parameter-query variants have an actual source");
 Some(match inner.as_ref(){ck::UffParameterQueryError::Cache(_)=>OperationError::from_ref(node),ck::UffParameterQueryError::Parameters(_)=>UffParameterError::from_ref(node)})
},source=Option<Source<Py<OperationError>,Py<UffParameterError>>>);
#[cfg(test)]
mod tests {
    use super::*;
    use crate::native::Molecule;
    use pyo3::types::PyDict;
    #[test]
    fn evaluation_optimization_and_coverage_match_native_defaults_options_results_and_failures() {
        fn check_archive(
            actual: &ck::Molecule,
            native: &ck::Molecule,
            family: &str,
            form: usize,
            failures: &mut Vec<String>,
        ) {
            let bits = |molecule: &ck::Molecule| {
                molecule
                    .conformers_3d()
                    .iter()
                    .map(|row| {
                        (
                            row.id(),
                            row.coordinates()
                                .iter()
                                .map(|xyz| xyz.map(f64::to_bits))
                                .collect::<Vec<_>>(),
                        )
                    })
                    .collect::<Vec<_>>()
            };
            assert_eq!(bits(actual), bits(native));
            assert_eq!(actual.topology(), native.topology());
            assert_eq!(actual.properties(), native.properties());
            match (actual.to_binary(), native.to_binary()) {
                (Ok(actual), Ok(native)) => assert_eq!(actual, native),
                (Err(actual), Err(native)) => {
                    assert_eq!(actual.to_string(), native.to_string());
                    let evidence = format!(
                        "{family} form={form} actual_archive={actual:?} native_archive={native:?}"
                    );
                    eprintln!("{evidence}");
                    failures.push(evidence);
                }
                (actual, native) => panic!(
                    "archive outcome diverged: {family} form={form} actual={actual:?} native={native:?}"
                ),
            }
        }
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let mut archive_failures = Vec::new();
            let mut embed = ck::EmbedParams::etkdg_v3();
            embed.random_seed = 17;
            let source = ck::Molecule::from_smiles("CCO")
                .unwrap()
                .with_hydrogens()
                .unwrap()
                .with_3d_conformers_with_params(2, &embed)
                .unwrap();
            let before = source.to_binary().unwrap();
            let projected = Py::new(
                py,
                Molecule {
                    inner: source.clone(),
                },
            )?;
            assert_eq!(
                projected
                    .bind(py)
                    .call_method0("uff_has_all_molecule_params")?
                    .extract::<bool>()?,
                source.uff_has_all_molecule_params().unwrap()
            );
            {
                let mut config = ck::UffEvaluationParams::default();
                config.vdw_threshold = 7.5;
                let default = source.uff_energy_gradient().unwrap();
                let configured = source.uff_energy_gradient_with_params(&config).unwrap();
                let params = Py::new(py, crate::params::UffEvaluationParams { inner: config })?;
                let keywords = PyDict::new(py);
                keywords.set_item("vdw_threshold", 7.5)?;
                for form in 0..4 {
                    let result = match form {
                        0 => projected.bind(py).call_method0("uff_energy_gradient")?,
                        1 => projected
                            .bind(py)
                            .call_method1("uff_energy_gradient", (params.clone_ref(py),))?,
                        2 => projected.bind(py).call_method(
                            "uff_energy_gradient",
                            (),
                            Some(&keywords),
                        )?,
                        _ => projected.bind(py).call_method1(
                            "uff_energy_gradient_with_params",
                            (params.clone_ref(py),),
                        )?,
                    };
                    let want = if form == 0 { &default } else { &configured };
                    let actual =
                        result.extract::<Py<crate::numerical_results::UffEnergyGradient>>()?;
                    assert_eq!(&actual.borrow(py).inner, want);
                    let actual = actual.borrow(py);
                    assert_eq!(actual.inner.energy().to_bits(), want.energy().to_bits());
                    assert_eq!(
                        actual
                            .inner
                            .gradient()
                            .iter()
                            .map(|v| v.to_bits())
                            .collect::<Vec<_>>(),
                        want.gradient()
                            .iter()
                            .map(|v| v.to_bits())
                            .collect::<Vec<_>>()
                    );
                    assert_eq!(projected.borrow(py).inner.to_binary().unwrap(), before);
                }
                let mixed = projected
                    .bind(py)
                    .call_method("uff_energy_gradient", (params,), Some(&keywords))
                    .unwrap_err();
                assert!(mixed.is_instance_of::<pyo3::exceptions::PyTypeError>(py));
            }
            {
                let mut config = ck::MmffEvaluationParams::default();
                config.non_bonded_threshold = 7.5;
                let default = source.mmff_energy_gradient().unwrap();
                let configured = source.mmff_energy_gradient_with_params(&config).unwrap();
                let params = Py::new(py, crate::params::MmffEvaluationParams { inner: config })?;
                let keywords = PyDict::new(py);
                keywords.set_item("non_bonded_threshold", 7.5)?;
                for form in 0..4 {
                    let result = match form {
                        0 => projected.bind(py).call_method0("mmff_energy_gradient")?,
                        1 => projected
                            .bind(py)
                            .call_method1("mmff_energy_gradient", (params.clone_ref(py),))?,
                        2 => projected.bind(py).call_method(
                            "mmff_energy_gradient",
                            (),
                            Some(&keywords),
                        )?,
                        _ => projected.bind(py).call_method1(
                            "mmff_energy_gradient_with_params",
                            (params.clone_ref(py),),
                        )?,
                    };
                    let want = if form == 0 { &default } else { &configured };
                    let actual = result
                        .extract::<Option<Py<crate::numerical_results::MmffEnergyGradient>>>()?;
                    match (actual, want) {
                        (Some(actual), Some(want)) => {
                            let actual = actual.borrow(py);
                            assert_eq!(&actual.inner, want);
                            assert_eq!(actual.inner.energy().to_bits(), want.energy().to_bits());
                            assert_eq!(
                                actual
                                    .inner
                                    .gradient()
                                    .iter()
                                    .map(|v| v.to_bits())
                                    .collect::<Vec<_>>(),
                                want.gradient()
                                    .iter()
                                    .map(|v| v.to_bits())
                                    .collect::<Vec<_>>()
                            );
                        }
                        (None, None) => {}
                        _ => panic!("native optional MMFF result diverged"),
                    }
                    assert_eq!(projected.borrow(py).inner.to_binary().unwrap(), before);
                }
                let mixed = projected
                    .bind(py)
                    .call_method("mmff_energy_gradient", (params,), Some(&keywords))
                    .unwrap_err();
                assert!(mixed.is_instance_of::<pyo3::exceptions::PyTypeError>(py));
            }
            {
                let mut config = ck::UffOptimizationParams::default();
                config.max_iterations = 3;
                let default = source.with_uff_optimized().unwrap();
                let configured = source.with_uff_optimized_with_params(&config).unwrap();
                let params = Py::new(py, crate::params::UffOptimizationParams { inner: config })?;
                let keywords = PyDict::new(py);
                keywords.set_item("max_iterations", 3)?;
                for form in 0..4 {
                    let result = match form {
                        0 => projected.bind(py).call_method0("with_uff_optimized")?,
                        1 => projected
                            .bind(py)
                            .call_method1("with_uff_optimized", (params.clone_ref(py),))?,
                        2 => projected.bind(py).call_method(
                            "with_uff_optimized",
                            (),
                            Some(&keywords),
                        )?,
                        _ => projected.bind(py).call_method1(
                            "with_uff_optimized_with_params",
                            (params.clone_ref(py),),
                        )?,
                    };
                    let want = if form == 0 { &default } else { &configured };
                    let actual =
                        result.extract::<Py<crate::numerical_results::UffOptimizationResult>>()?;
                    assert_eq!(&actual.borrow(py).inner, want);
                    let actual = actual.borrow(py);
                    assert_eq!(
                        actual.inner.molecule().to_binary().unwrap(),
                        want.molecule().to_binary().unwrap()
                    );
                    assert_eq!(actual.inner.energy().to_bits(), want.energy().to_bits());
                    assert_eq!(projected.borrow(py).inner.to_binary().unwrap(), before);
                }
                let mixed = projected
                    .bind(py)
                    .call_method("with_uff_optimized", (params,), Some(&keywords))
                    .unwrap_err();
                assert!(mixed.is_instance_of::<pyo3::exceptions::PyTypeError>(py));
            }
            {
                let mut config = ck::UffConformerOptimizationParams::default();
                config.max_iterations = 3;
                let default = source.with_uff_optimized_conformers().unwrap();
                let configured = source
                    .with_uff_optimized_conformers_with_params(&config)
                    .unwrap();
                let params = Py::new(
                    py,
                    crate::params::UffConformerOptimizationParams { inner: config },
                )?;
                let keywords = PyDict::new(py);
                keywords.set_item("max_iterations", 3)?;
                for form in 0..4 {
                    let result = match form {
                        0 => projected
                            .bind(py)
                            .call_method0("with_uff_optimized_conformers")?,
                        1 => projected.bind(py).call_method1(
                            "with_uff_optimized_conformers",
                            (params.clone_ref(py),),
                        )?,
                        2 => projected.bind(py).call_method(
                            "with_uff_optimized_conformers",
                            (),
                            Some(&keywords),
                        )?,
                        _ => projected.bind(py).call_method1(
                            "with_uff_optimized_conformers_with_params",
                            (params.clone_ref(py),),
                        )?,
                    };
                    let want = if form == 0 { &default } else { &configured };
                    let actual = result
                        .extract::<Py<crate::numerical_results::UffConformerOptimizationResult>>(
                        )?;
                    assert_eq!(&actual.borrow(py).inner, want);
                    let actual = actual.borrow(py);
                    assert_eq!(
                        actual.inner.molecule().to_binary().unwrap(),
                        want.molecule().to_binary().unwrap()
                    );
                    assert_eq!(
                        actual
                            .inner
                            .conformer_results()
                            .iter()
                            .map(|v| v.energy().to_bits())
                            .collect::<Vec<_>>(),
                        want.conformer_results()
                            .iter()
                            .map(|v| v.energy().to_bits())
                            .collect::<Vec<_>>()
                    );
                    assert_eq!(projected.borrow(py).inner.to_binary().unwrap(), before);
                }
                let mixed = projected
                    .bind(py)
                    .call_method("with_uff_optimized_conformers", (params,), Some(&keywords))
                    .unwrap_err();
                assert!(mixed.is_instance_of::<pyo3::exceptions::PyTypeError>(py));
            }
            {
                let mut config = ck::MmffOptimizationParams::default();
                config.max_iterations = 3;
                let default = source.with_mmff_optimized().unwrap();
                let configured = source.with_mmff_optimized_with_params(&config).unwrap();
                let params = Py::new(py, crate::params::MmffOptimizationParams { inner: config })?;
                let keywords = PyDict::new(py);
                keywords.set_item("max_iterations", 3)?;
                for form in 0..4 {
                    let result = match form {
                        0 => projected.bind(py).call_method0("with_mmff_optimized")?,
                        1 => projected
                            .bind(py)
                            .call_method1("with_mmff_optimized", (params.clone_ref(py),))?,
                        2 => projected.bind(py).call_method(
                            "with_mmff_optimized",
                            (),
                            Some(&keywords),
                        )?,
                        _ => projected.bind(py).call_method1(
                            "with_mmff_optimized_with_params",
                            (params.clone_ref(py),),
                        )?,
                    };
                    let want = if form == 0 { &default } else { &configured };
                    let actual = result
                        .extract::<Py<crate::numerical_results::MmffOptimizeMoleculeResult>>()?;
                    assert_eq!(&actual.borrow(py).inner, want);
                    let actual = actual.borrow(py);
                    let actual_binary = actual.inner.molecule().to_binary();
                    let expected_binary = want.molecule().to_binary();
                    let bits = |molecule: &ck::Molecule| molecule.conformers_3d().iter()
                        .map(|row| (row.id(), row.coordinates().iter()
                            .map(|xyz| xyz.map(f64::to_bits)).collect::<Vec<_>>()))
                        .collect::<Vec<_>>();
                    let actual_bits = bits(actual.inner.molecule());
                    let expected_bits = bits(want.molecule());
                    if actual_binary.is_err() || expected_binary.is_err() {
                        eprintln!("MMFF form={form} actual_archive={actual_binary:?} native_archive={expected_binary:?} actual_status={} native_status={} actual_coordinate_bits={actual_bits:?} native_coordinate_bits={expected_bits:?}", actual.inner.status_code(), want.status_code());
                    }
                    assert_eq!(actual.inner.status_code(), want.status_code());
                    assert_eq!(actual_bits, expected_bits);
                    assert_eq!(actual.inner.molecule().topology(), want.molecule().topology());
                    assert_eq!(actual.inner.molecule().properties(), want.molecule().properties());
                    check_archive(actual.inner.molecule(), want.molecule(), "with_mmff_optimized", form, &mut archive_failures);
                    assert_eq!(projected.borrow(py).inner.to_binary().unwrap(), before);
                }
                let mixed = projected
                    .bind(py)
                    .call_method("with_mmff_optimized", (params,), Some(&keywords))
                    .unwrap_err();
                assert!(mixed.is_instance_of::<pyo3::exceptions::PyTypeError>(py));
            }
            {
                let mut config = ck::MmffConformerOptimizationParams::default();
                config.max_iterations = 3;
                let default = source.with_mmff_optimized_conformers().unwrap();
                let configured = source
                    .with_mmff_optimized_conformers_with_params(&config)
                    .unwrap();
                let params = Py::new(
                    py,
                    crate::params::MmffConformerOptimizationParams { inner: config },
                )?;
                let keywords = PyDict::new(py);
                keywords.set_item("max_iterations", 3)?;
                for form in 0..4 {
                    let result = match form {
                        0 => projected
                            .bind(py)
                            .call_method0("with_mmff_optimized_conformers")?,
                        1 => projected.bind(py).call_method1(
                            "with_mmff_optimized_conformers",
                            (params.clone_ref(py),),
                        )?,
                        2 => projected.bind(py).call_method(
                            "with_mmff_optimized_conformers",
                            (),
                            Some(&keywords),
                        )?,
                        _ => projected.bind(py).call_method1(
                            "with_mmff_optimized_conformers_with_params",
                            (params.clone_ref(py),),
                        )?,
                    };
                    let want = if form == 0 { &default } else { &configured };
                    let actual = result
                        .extract::<Py<crate::numerical_results::MmffOptimizeMoleculeConfsResult>>(
                        )?;
                    assert_eq!(&actual.borrow(py).inner, want);
                    let actual = actual.borrow(py);
                    check_archive(actual.inner.molecule(), want.molecule(), "with_mmff_optimized_conformers", form, &mut archive_failures);
                    assert_eq!(
                        actual
                            .inner
                            .conformer_results()
                            .iter()
                            .map(|v| v.energy().to_bits())
                            .collect::<Vec<_>>(),
                        want.conformer_results()
                            .iter()
                            .map(|v| v.energy().to_bits())
                            .collect::<Vec<_>>()
                    );
                    assert_eq!(projected.borrow(py).inner.to_binary().unwrap(), before);
                }
                let mixed = projected
                    .bind(py)
                    .call_method("with_mmff_optimized_conformers", (params,), Some(&keywords))
                    .unwrap_err();
                assert!(mixed.is_instance_of::<pyo3::exceptions::PyTypeError>(py));
            }

            let empty = ck::Molecule::new();
            let actual = Py::new(
                py,
                Molecule {
                    inner: empty.clone(),
                },
            )?;
            assert_eq!(
                actual
                    .bind(py)
                    .call_method0("uff_has_all_molecule_params")?
                    .extract::<bool>()?,
                empty.uff_has_all_molecule_params().unwrap()
            );
            let detached = ck::Molecule::from_parts(
                source.topology().clone(),
                source.to_builder().coordinates().clone(),
                source.properties().clone(),
            )
            .unwrap();
            let expected = detached.uff_has_all_molecule_params().unwrap_err();
            let actual = Py::new(py, Molecule { inner: detached })?;
            let error = actual
                .bind(py)
                .call_method0("uff_has_all_molecule_params")
                .unwrap_err();
            assert!(error.is_instance_of::<UffParameterQueryError>(py));
            assert_eq!(error.value(py).str()?.to_str()?, expected.to_string());
            assert!(
                error
                    .cause(py)
                    .unwrap()
                    .is_instance_of::<crate::operation_errors::OperationError>(py)
            );
            assert!(
                error
                    .value(py)
                    .getattr("source")?
                    .is(error.cause(py).unwrap().value(py))
            );
            assert!(error.value(py).setattr("variant", 0).is_err());
            let missing = PyDict::new(py);
            missing.set_item("conformer_id", usize::MAX)?;
            let expected = source
                .uff_energy_gradient_with_params(&ck::UffEvaluationParams {
                    conformer_id: Some(usize::MAX),
                    ..Default::default()
                })
                .unwrap_err();
            let error = projected
                .bind(py)
                .call_method("uff_energy_gradient", (), Some(&missing))
                .unwrap_err();
            assert!(error.is_instance_of::<crate::operation_errors::OperationError>(py));
            assert_eq!(error.value(py).str()?.to_str()?, expected.to_string());
            assert_eq!(projected.borrow(py).inner.to_binary().unwrap(), before);
            assert!(archive_failures.is_empty(), "strict successful-archive assertion remains failed: {archive_failures:?}");
            Ok(())
        })
        .unwrap();
    }
}
