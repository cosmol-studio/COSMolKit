//! Layered projection: options and detached results, never path computation.
use crate::collection_params::LayeredFingerprintParams;
use crate::layered_fingerprint_errors::LayeredFingerprintError;
use crate::native::Fingerprint;
use crate::ordered_results::LayeredFingerprintResult;
use crate::query_compiler::QueryGraph;
use cosmolkit as ck;
use pyo3::prelude::*;

pub(crate) fn keyword_params(
    layers: Option<u32>,
    min_path: Option<u32>,
    max_path: Option<u32>,
    fp_size: Option<u32>,
    atom_counts: Option<Vec<u32>>,
    set_only_bits: Option<&Fingerprint>,
    branched_paths: Option<bool>,
    from_atoms: Option<Vec<u32>>,
) -> ck::LayeredFingerprintParams {
    let defaults = ck::LayeredFingerprintParams::default();
    ck::LayeredFingerprintParams {
        layers: layers.map_or(
            defaults.layers,
            ck::LayeredFingerprintLayers::from_bits_retain,
        ),
        min_path: min_path.unwrap_or(defaults.min_path),
        max_path: max_path.unwrap_or(defaults.max_path),
        fp_size: fp_size.unwrap_or(defaults.fp_size),
        atom_counts,
        set_only_bits: set_only_bits.map(|value| value.inner.clone()),
        branched_paths: branched_paths.unwrap_or(defaults.branched_paths),
        from_atoms,
    }
}

#[pyfunction]
pub fn fingerprint_layered_query_with_params(
    py: Python<'_>,
    query: &QueryGraph,
    params: &LayeredFingerprintParams,
) -> PyResult<Fingerprint> {
    ck::fingerprint_layered_query_with_params(query.inner(), &params.to_rust(py)?)
        .map(|inner| Fingerprint { inner })
        .map_err(LayeredFingerprintError::from_rust)
}

#[pyfunction]
pub fn fingerprint_layered_query_with_output_with_params(
    py: Python<'_>,
    query: &QueryGraph,
    params: &LayeredFingerprintParams,
) -> PyResult<LayeredFingerprintResult> {
    ck::fingerprint_layered_query_with_output_with_params(query.inner(), &params.to_rust(py)?)
        .map(|inner| LayeredFingerprintResult { inner })
        .map_err(LayeredFingerprintError::from_rust)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::native::Molecule;
    use pyo3::types::{PyDict, PyModule};

    fn exact(py: Python<'_>, actual: Py<Fingerprint>, expected: &ck::Fingerprint) {
        let actual = actual.borrow(py);
        assert_eq!(actual.inner.n_bits(), expected.n_bits());
        assert_eq!(actual.inner.on_bits(), expected.on_bits());
    }

    #[test]
    fn actual_layered_producers_preserve_masks_seeds_roots_and_owned_results() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let module = PyModule::new(py, "layered_probe")?;
            module.add_function(wrap_pyfunction!(
                fingerprint_layered_query_with_params,
                &module
            )?)?;
            module.add_function(wrap_pyfunction!(
                fingerprint_layered_query_with_output_with_params,
                &module
            )?)?;
            let query = Py::new(py, QueryGraph::from_rust(ck::parse_smarts("C~N").unwrap()))?;
            let query_before = query.borrow(py).inner().clone();
            for text in ["CCO", "c1ccccc1"] {
                let native = ck::Molecule::from_smiles(text).unwrap();
                let before = native.to_binary().unwrap();
                let count = native.num_atoms();
                let molecule = Py::new(py, Molecule { inner: native })?;
                let expected = molecule.borrow(py).inner.fingerprint_layered().unwrap();
                exact(
                    py,
                    molecule
                        .bind(py)
                        .call_method0("fingerprint_layered")?
                        .extract()?,
                    &expected,
                );
                let expected_result = molecule
                    .borrow(py)
                    .inner
                    .fingerprint_layered_with_output()
                    .unwrap();
                let default_result = molecule
                    .bind(py)
                    .call_method0("fingerprint_layered_with_output")?
                    .extract::<Py<LayeredFingerprintResult>>()?;
                assert_eq!(default_result.borrow(py).inner, expected_result);
                for native_params in [
                    ck::LayeredFingerprintParams::default(),
                    ck::LayeredFingerprintParams {
                        layers: ck::LayeredFingerprintLayers::SUBSTRUCTURE,
                        fp_size: 128,
                        max_path: 3,
                        ..Default::default()
                    },
                    ck::LayeredFingerprintParams {
                        layers: ck::LayeredFingerprintLayers::from_bits_retain(0x8000_003f),
                        atom_counts: Some(vec![u32::MAX; count]),
                        from_atoms: Some(vec![0, 0]),
                        ..Default::default()
                    },
                    ck::LayeredFingerprintParams {
                        atom_counts: Some(vec![7; count]),
                        from_atoms: Some(vec![]),
                        ..Default::default()
                    },
                    ck::LayeredFingerprintParams {
                        atom_counts: Some(vec![10; count]),
                        set_only_bits: Some(ck::Fingerprint::from_on_bits(2048, [0, 674]).unwrap()),
                        ..Default::default()
                    },
                    ck::LayeredFingerprintParams {
                        set_only_bits: Some(ck::Fingerprint::new(2048)),
                        branched_paths: false,
                        from_atoms: Some(vec![0]),
                        ..Default::default()
                    },
                ] {
                    let original_params = native_params.clone();
                    let expected = molecule
                        .borrow(py)
                        .inner
                        .fingerprint_layered_with_output_with_params(&native_params)
                        .unwrap();
                    let expected_query = ck::fingerprint_layered_query_with_output_with_params(
                        query.borrow(py).inner(),
                        &ck::LayeredFingerprintParams {
                            atom_counts: None,
                            ..native_params.clone()
                        },
                    )
                    .unwrap();
                    let params =
                        Py::new(py, LayeredFingerprintParams::from_rust(py, native_params)?)?;
                    let kwargs = PyDict::new(py);
                    {
                        let borrowed = params.borrow(py);
                        let p = borrowed.to_rust(py)?;
                        kwargs.set_item("layers", p.layers.bits())?;
                        kwargs.set_item("min_path", p.min_path)?;
                        kwargs.set_item("max_path", p.max_path)?;
                        kwargs.set_item("fp_size", p.fp_size)?;
                        kwargs.set_item("atom_counts", p.atom_counts.clone())?;
                        kwargs.set_item(
                            "set_only_bits",
                            p.set_only_bits.clone().map(|inner| Fingerprint { inner }),
                        )?;
                        kwargs.set_item("branched_paths", p.branched_paths)?;
                        kwargs.set_item("from_atoms", p.from_atoms.clone())?;
                    }
                    for actual in [
                        molecule
                            .bind(py)
                            .call_method1("fingerprint_layered", (params.clone_ref(py),))?,
                        molecule
                            .bind(py)
                            .call_method("fingerprint_layered", (), Some(&kwargs))?,
                        molecule.bind(py).call_method1(
                            "fingerprint_layered_with_params",
                            (params.clone_ref(py),),
                        )?,
                    ] {
                        exact(py, actual.extract()?, &expected.fingerprint);
                    }
                    let mut retained = Vec::new();
                    for actual in [
                        molecule.bind(py).call_method1(
                            "fingerprint_layered_with_output",
                            (params.clone_ref(py),),
                        )?,
                        molecule.bind(py).call_method(
                            "fingerprint_layered_with_output",
                            (),
                            Some(&kwargs),
                        )?,
                        molecule.bind(py).call_method1(
                            "fingerprint_layered_with_output_with_params",
                            (params.clone_ref(py),),
                        )?,
                    ] {
                        let actual = actual.extract::<Py<LayeredFingerprintResult>>()?;
                        assert_eq!(actual.borrow(py).inner, expected);
                        retained.push(actual);
                    }
                    assert_eq!(params.borrow(py).to_rust(py)?, original_params);
                    assert!(
                        molecule
                            .bind(py)
                            .call_method(
                                "fingerprint_layered",
                                (params.clone_ref(py),),
                                Some(&kwargs)
                            )
                            .is_err()
                    );
                    assert!(
                        molecule
                            .bind(py)
                            .call_method(
                                "fingerprint_layered_with_output",
                                (params.clone_ref(py),),
                                Some(&kwargs)
                            )
                            .is_err()
                    );
                    let query_params = Py::new(
                        py,
                        LayeredFingerprintParams::from_rust(
                            py,
                            ck::LayeredFingerprintParams {
                                atom_counts: None,
                                ..original_params
                            },
                        )?,
                    )?;
                    exact(
                        py,
                        module
                            .getattr("fingerprint_layered_query_with_params")?
                            .call1((query.clone_ref(py), query_params.clone_ref(py)))?
                            .extract()?,
                        &expected_query.fingerprint,
                    );
                    let actual_query = module
                        .getattr("fingerprint_layered_query_with_output_with_params")?
                        .call1((query.clone_ref(py), query_params))?
                        .extract::<Py<LayeredFingerprintResult>>()?;
                    assert_eq!(actual_query.borrow(py).inner, expected_query);
                    drop(params);
                    for actual in retained {
                        let counts = actual
                            .bind(py)
                            .call_method0("atom_counts")?
                            .extract::<Option<Vec<u32>>>()?;
                        assert_eq!(counts, expected.atom_counts);
                    }
                }
                assert_eq!(molecule.borrow(py).inner.to_binary().unwrap(), before);
                drop(molecule);
                assert_eq!(default_result.borrow(py).inner, expected_result);
            }
            assert_eq!(query.borrow(py).inner(), &query_before);
            let molecule = Py::new(
                py,
                Molecule {
                    inner: ck::Molecule::from_smiles("CCO").unwrap(),
                },
            )?;
            // The native unrooted linear branch uses atom selection as bond
            // indices (layered.rs LayeredFingerprintParams.branched_paths).
            // Keep its actual failure as well as the rooted success above.
            let expected = molecule
                .borrow(py)
                .inner
                .fingerprint_layered_with_params(&ck::LayeredFingerprintParams {
                    branched_paths: false,
                    ..Default::default()
                })
                .unwrap_err();
            let linear = PyDict::new(py);
            linear.set_item("branched_paths", false)?;
            for name in ["fingerprint_layered", "fingerprint_layered_with_output"] {
                let actual = molecule
                    .bind(py)
                    .call_method(name, (), Some(&linear))
                    .unwrap_err();
                assert!(actual.is_instance(py, <LayeredFingerprintError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
                assert_eq!(actual.value(py).str()?.to_str()?, expected.to_string());
            }
            for kwargs in [
                {
                    let kw = PyDict::new(py);
                    kw.set_item("fp_size", 0)?;
                    kw
                },
                {
                    let kw = PyDict::new(py);
                    kw.set_item("atom_counts", Vec::<u32>::new())?;
                    kw
                },
                {
                    let kw = PyDict::new(py);
                    kw.set_item(
                        "set_only_bits",
                        Py::new(
                            py,
                            Fingerprint {
                                inner: ck::Fingerprint::new(64),
                            },
                        )?,
                    )?;
                    kw
                },
            ] {
                let actual = molecule
                    .bind(py)
                    .call_method("fingerprint_layered_with_output", (), Some(&kwargs))
                    .unwrap_err();
                assert!(actual.is_instance(py, <LayeredFingerprintError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
            }
            Ok(())
        })
        .unwrap();
    }
}
