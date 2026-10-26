//! RDKFingerprint projection. Native code owns hashing, folding and provenance.
use crate::batch_fingerprint_values::TopologicalFingerprintParams;
use crate::leaf_values::TopologicalFingerprintOutput;
use crate::native::Fingerprint;
use crate::params::TopologicalFingerprintOutputRequest;
use crate::query_compiler::QueryGraph;
use crate::topological_fingerprint_errors::TopologicalFingerprintError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::collections::BTreeMap;
use std::sync::Arc;

pub(crate) fn request_keywords(
    atom_bits: Option<bool>,
    bit_info: Option<bool>,
) -> ck::TopologicalFingerprintOutputRequest {
    let defaults = ck::TopologicalFingerprintOutputRequest::default();
    ck::TopologicalFingerprintOutputRequest {
        atom_bits: atom_bits.unwrap_or(defaults.atom_bits),
        bit_info: bit_info.unwrap_or(defaults.bit_info),
    }
}

pub(crate) fn keyword_params(
    min_path: Option<u32>,
    max_path: Option<u32>,
    fp_size: Option<u32>,
    num_bits_per_feature: Option<u32>,
    use_hs: Option<bool>,
    target_density: Option<f64>,
    min_size: Option<u32>,
    branched_paths: Option<bool>,
    use_bond_order: Option<bool>,
    atom_invariants: Option<Vec<u32>>,
    from_atoms: Option<Vec<u32>>,
    ignore_atoms: Option<Vec<u32>>,
) -> ck::TopologicalFingerprintParams {
    let defaults = ck::TopologicalFingerprintParams::default();
    ck::TopologicalFingerprintParams {
        min_path: min_path.unwrap_or(defaults.min_path),
        max_path: max_path.unwrap_or(defaults.max_path),
        fp_size: fp_size.unwrap_or(defaults.fp_size),
        num_bits_per_feature: num_bits_per_feature.unwrap_or(defaults.num_bits_per_feature),
        use_hs: use_hs.unwrap_or(defaults.use_hs),
        target_density: target_density.unwrap_or(defaults.target_density),
        min_size: min_size.unwrap_or(defaults.min_size),
        branched_paths: branched_paths.unwrap_or(defaults.branched_paths),
        use_bond_order: use_bond_order.unwrap_or(defaults.use_bond_order),
        atom_invariants,
        from_atoms,
        ignore_atoms,
    }
}

#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
/// Topological fingerprint together with explicitly requested contribution metadata.
pub struct TopologicalFingerprintResult {
    pub(crate) inner: Arc<ck::TopologicalFingerprintResult>,
}
impl TopologicalFingerprintResult {
    pub(crate) fn from_rust(inner: ck::TopologicalFingerprintResult) -> Self {
        Self {
            inner: Arc::new(inner),
        }
    }
}
#[pymethods]
impl TopologicalFingerprintResult {
    #[new]
    fn new(fingerprint: &Fingerprint, output: &TopologicalFingerprintOutput) -> Self {
        Self::from_rust(ck::TopologicalFingerprintResult {
            fingerprint: fingerprint.inner.clone(),
            output: output.inner.clone(),
        })
    }
    #[getter]
    /// Generated topological bit fingerprint.
    fn fingerprint(&self) -> Fingerprint {
        Fingerprint {
            inner: self.inner.fingerprint().clone(),
        }
    }
    #[getter]
    /// Bit identifiers attributed to each atom; raises the native error if this output was not collected.
    fn atom_bits(&self) -> PyResult<Vec<Vec<u32>>> {
        self.inner
            .atom_bits()
            .map(<[Vec<u32>]>::to_vec)
            .map_err(TopologicalFingerprintError::from_rust)
    }
    #[getter]
    /// Contributing bond paths keyed by fingerprint bit; raises the native error if this output was not collected.
    fn bit_info(&self) -> PyResult<BTreeMap<u32, Vec<Vec<i32>>>> {
        self.inner
            .bit_info()
            .cloned()
            .map_err(TopologicalFingerprintError::from_rust)
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}

#[pyfunction]
/// Calculate a path-based topological fingerprint for a QueryGraph using explicit parameters.
pub fn fingerprint_topological_query_with_params(
    py: Python<'_>,
    query: &QueryGraph,
    params: &TopologicalFingerprintParams,
) -> PyResult<Fingerprint> {
    ck::fingerprint_topological_query_with_params(query.inner(), &params.to_rust(py)?)
        .map(|inner| Fingerprint { inner })
        .map_err(TopologicalFingerprintError::from_rust)
}
#[pyfunction]
/// Calculate a QueryGraph topological fingerprint together with the requested atom/bit provenance.
pub fn fingerprint_topological_query_with_output_with_params(
    py: Python<'_>,
    query: &QueryGraph,
    params: &TopologicalFingerprintParams,
    request: &TopologicalFingerprintOutputRequest,
) -> PyResult<TopologicalFingerprintResult> {
    let request = pyo3::Python::attach(|py| request.to_rust(py))?;

    ck::fingerprint_topological_query_with_output_with_params(
        query.inner(),
        &params.to_rust(py)?,
        request,
    )
    .map(TopologicalFingerprintResult::from_rust)
    .map_err(TopologicalFingerprintError::from_rust)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::native::Molecule;
    use pyo3::types::{PyDict, PyModule};

    fn exact(py: Python<'_>, actual: Py<Fingerprint>, expected: &ck::Fingerprint) {
        assert_eq!(actual.borrow(py).inner.n_bits(), expected.n_bits());
        assert_eq!(actual.borrow(py).inner.on_bits(), expected.on_bits());
    }

    #[test]
    fn actual_topological_producers_preserve_folding_options_provenance_and_native_failures() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let module = PyModule::new(py, "topological_probe")?;
            module.add_function(wrap_pyfunction!(
                fingerprint_topological_query_with_params,
                &module
            )?)?;
            module.add_function(wrap_pyfunction!(
                fingerprint_topological_query_with_output_with_params,
                &module
            )?)?;
            let query = Py::new(py, QueryGraph::from_rust(ck::parse_smarts("C~N").unwrap()))?;
            let original_query = query.borrow(py).inner().clone();
            for text in ["CCO", "c1ccccc1"] {
                let native = ck::Molecule::from_smiles(text).unwrap();
                let before = native.to_binary().unwrap();
                let count = native.num_atoms();
                let molecule = Py::new(py, Molecule { inner: native })?;
                let baseline = molecule.borrow(py).inner.fingerprint_topological().unwrap();
                exact(
                    py,
                    molecule
                        .bind(py)
                        .call_method0("fingerprint_topological")?
                        .extract()?,
                    &baseline,
                );
                let default_result = molecule
                    .bind(py)
                    .call_method0("fingerprint_topological_with_output")?
                    .extract::<Py<TopologicalFingerprintResult>>()?;
                assert_eq!(
                    default_result.borrow(py).inner.as_ref(),
                    &molecule
                        .borrow(py)
                        .inner
                        .fingerprint_topological_with_output()
                        .unwrap()
                );
                for native_params in [
                    ck::TopologicalFingerprintParams::default(),
                    ck::TopologicalFingerprintParams {
                        fp_size: 128,
                        max_path: 3,
                        num_bits_per_feature: 4,
                        use_hs: false,
                        ..Default::default()
                    },
                    ck::TopologicalFingerprintParams {
                        fp_size: 512,
                        target_density: 0.5,
                        min_size: 64,
                        branched_paths: false,
                        ..Default::default()
                    },
                    ck::TopologicalFingerprintParams {
                        atom_invariants: Some(vec![7; count]),
                        use_bond_order: false,
                        ..Default::default()
                    },
                    ck::TopologicalFingerprintParams {
                        from_atoms: Some(vec![]),
                        ignore_atoms: Some(vec![]),
                        ..Default::default()
                    },
                    ck::TopologicalFingerprintParams {
                        from_atoms: Some(vec![0, 0]),
                        ignore_atoms: Some(vec![1, 1]),
                        ..Default::default()
                    },
                ] {
                    let original_params = native_params.clone();
                    let expected = molecule
                        .borrow(py)
                        .inner
                        .fingerprint_topological_with_params(&native_params)
                        .unwrap();
                    let params = Py::new(
                        py,
                        TopologicalFingerprintParams::from_rust(py, native_params)?,
                    )?;
                    let kwargs = PyDict::new(py);
                    {
                        let value = params.borrow(py);
                        let p = value.to_rust(py)?;
                        kwargs.set_item("min_path", p.min_path)?;
                        kwargs.set_item("max_path", p.max_path)?;
                        kwargs.set_item("fp_size", p.fp_size)?;
                        kwargs.set_item("num_bits_per_feature", p.num_bits_per_feature)?;
                        kwargs.set_item("use_hs", p.use_hs)?;
                        kwargs.set_item("target_density", p.target_density)?;
                        kwargs.set_item("min_size", p.min_size)?;
                        kwargs.set_item("branched_paths", p.branched_paths)?;
                        kwargs.set_item("use_bond_order", p.use_bond_order)?;
                        kwargs.set_item("atom_invariants", p.atom_invariants.clone())?;
                        kwargs.set_item("from_atoms", p.from_atoms.clone())?;
                        kwargs.set_item("ignore_atoms", p.ignore_atoms.clone())?;
                    }
                    for actual in [
                        molecule
                            .bind(py)
                            .call_method1("fingerprint_topological", (params.clone_ref(py),))?,
                        molecule.bind(py).call_method(
                            "fingerprint_topological",
                            (),
                            Some(&kwargs),
                        )?,
                        molecule.bind(py).call_method1(
                            "fingerprint_topological_with_params",
                            (params.clone_ref(py),),
                        )?,
                    ] {
                        exact(py, actual.extract()?, &expected);
                    }
                    assert!(
                        molecule
                            .bind(py)
                            .call_method(
                                "fingerprint_topological",
                                (params.clone_ref(py),),
                                Some(&kwargs)
                            )
                            .is_err()
                    );
                    let expected_query = ck::fingerprint_topological_query_with_params(
                        query.borrow(py).inner(),
                        &original_params,
                    )
                    .unwrap();
                    exact(
                        py,
                        module
                            .getattr("fingerprint_topological_query_with_params")?
                            .call1((query.clone_ref(py), params.clone_ref(py)))?
                            .extract()?,
                        &expected_query,
                    );
                    for (atom_bits, bit_info) in
                        [(false, false), (true, false), (false, true), (true, true)]
                    {
                        let request = ck::TopologicalFingerprintOutputRequest {
                            atom_bits,
                            bit_info,
                        };
                        let expected = molecule
                            .borrow(py)
                            .inner
                            .fingerprint_topological_with_output_with_params(
                                &original_params,
                                request,
                            )
                            .unwrap();
                        let expected_query =
                            ck::fingerprint_topological_query_with_output_with_params(
                                query.borrow(py).inner(),
                                &original_params,
                                request,
                            )
                            .unwrap();
                        let request = Py::new(
                            py,
                            TopologicalFingerprintOutputRequest {
                                inner: crate::parameters::NativeField::Owned(request),
                            },
                        )?;
                        kwargs.set_item("request", request.clone_ref(py))?;
                        let mut retained = Vec::new();
                        for actual in [
                            molecule.bind(py).call_method1(
                                "fingerprint_topological_with_output",
                                (params.clone_ref(py), request.clone_ref(py)),
                            )?,
                            molecule.bind(py).call_method(
                                "fingerprint_topological_with_output",
                                (),
                                Some(&kwargs),
                            )?,
                            molecule.bind(py).call_method1(
                                "fingerprint_topological_with_output_with_params",
                                (params.clone_ref(py), request.clone_ref(py)),
                            )?,
                        ] {
                            let actual = actual.extract::<Py<TopologicalFingerprintResult>>()?;
                            assert_eq!(actual.borrow(py).inner.as_ref(), &expected);
                            exact(
                                py,
                                actual.bind(py).getattr("fingerprint")?.extract()?,
                                &expected.fingerprint,
                            );
                            if atom_bits {
                                assert_eq!(
                                    &actual
                                        .bind(py)
                                        .getattr("atom_bits")?
                                        .extract::<Vec<Vec<u32>>>()?,
                                    expected.atom_bits().unwrap()
                                );
                            } else {
                                let error = actual.bind(py).getattr("atom_bits").unwrap_err();
                                assert!(error.is_instance(py, <TopologicalFingerprintError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
                                assert_eq!(
                                    error.value(py).getattr("field")?.extract::<String>()?,
                                    "atom_bits"
                                );
                            }
                            if bit_info {
                                assert_eq!(
                                    &actual
                                        .bind(py)
                                        .getattr("bit_info")?
                                        .extract::<BTreeMap<u32, Vec<Vec<i32>>>>()?,
                                    expected.bit_info().unwrap()
                                );
                            } else {
                                let error = actual.bind(py).getattr("bit_info").unwrap_err();
                                assert!(error.is_instance(py, <TopologicalFingerprintError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
                                assert_eq!(
                                    error.value(py).getattr("field")?.extract::<String>()?,
                                    "bit_info"
                                );
                            }
                            retained.push(actual);
                        }
                        let actual_query = module
                            .getattr("fingerprint_topological_query_with_output_with_params")?
                            .call1((
                                query.clone_ref(py),
                                params.clone_ref(py),
                                request.clone_ref(py),
                            ))?
                            .extract::<Py<TopologicalFingerprintResult>>()?;
                        assert_eq!(actual_query.borrow(py).inner.as_ref(), &expected_query);
                        assert!(
                            molecule
                                .bind(py)
                                .call_method(
                                    "fingerprint_topological_with_output",
                                    (params.clone_ref(py),),
                                    Some(&kwargs)
                                )
                                .is_err()
                        );
                        drop(request);
                        for result in retained {
                            assert_eq!(result.borrow(py).inner.as_ref(), &expected);
                        }
                    }
                    assert_eq!(params.borrow(py).to_rust(py)?, original_params);
                }
                assert_eq!(molecule.borrow(py).inner.to_binary().unwrap(), before);
                drop(molecule);
                exact(
                    py,
                    default_result
                        .bind(py)
                        .getattr("fingerprint")?
                        .extract()?,
                    &baseline,
                );
            }
            assert_eq!(query.borrow(py).inner(), &original_query);
            let molecule = Py::new(
                py,
                Molecule {
                    inner: ck::Molecule::from_smiles("CCO").unwrap(),
                },
            )?;
            for native in [
                ck::TopologicalFingerprintParams {
                    fp_size: 0,
                    ..Default::default()
                },
                ck::TopologicalFingerprintParams {
                    atom_invariants: Some(vec![]),
                    ..Default::default()
                },
                ck::TopologicalFingerprintParams {
                    ignore_atoms: Some(vec![u32::MAX]),
                    ..Default::default()
                },
                ck::TopologicalFingerprintParams {
                    target_density: f64::NAN,
                    ..Default::default()
                },
            ] {
                let expected = molecule
                    .borrow(py)
                    .inner
                    .fingerprint_topological_with_params(&native)
                    .unwrap_err();
                let params = Py::new(py, TopologicalFingerprintParams::from_rust(py, native)?)?;
                for method in [
                    "fingerprint_topological",
                    "fingerprint_topological_with_output",
                ] {
                    let error = molecule
                        .bind(py)
                        .call_method1(method, (params.clone_ref(py),))
                        .unwrap_err();
                    assert!(error.is_instance(py, <TopologicalFingerprintError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
                    assert_eq!(error.value(py).str()?.to_str()?, expected.to_string());
                }
            }
            Ok(())
        })
        .unwrap();
    }
}
