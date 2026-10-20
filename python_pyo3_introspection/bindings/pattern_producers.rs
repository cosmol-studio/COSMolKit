//! Pattern call-form projection; molecule and query computation stay in the facade.
use crate::native::Fingerprint;
use crate::params::PatternFingerprintParams;
use crate::pattern_fingerprint_errors::PatternFingerprintError;
use crate::query_compiler::QueryGraph;
use cosmolkit as ck;
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;

pub(crate) fn keyword_params(
    n_bits: Option<usize>,
    tautomeric: Option<bool>,
) -> PatternFingerprintParams {
    let defaults = ck::PatternFingerprintParams::default();
    PatternFingerprintParams {
        inner: crate::parameters::NativeField::Owned(ck::PatternFingerprintParams {
            n_bits: n_bits.unwrap_or(defaults.n_bits),
            tautomeric: tautomeric.unwrap_or(defaults.tautomeric),
        }),
    }
}

#[pyfunction]
#[pyo3(signature = (query, params=crate::arguments::Argument::Missing, *, n_bits=crate::arguments::Argument::Missing, tautomeric=crate::arguments::Argument::Missing))]
pub fn fingerprint_pattern_query(
    query: &QueryGraph,
    params: crate::arguments::Argument<Option<PyRef<'_, PatternFingerprintParams>>>,
    n_bits: crate::arguments::Argument<Option<usize>>,
    tautomeric: crate::arguments::Argument<Option<bool>>,
) -> PyResult<Fingerprint> {
    let keyword_form = n_bits.is_supplied() || tautomeric.is_supplied();
    crate::arguments::configuration_arguments!(params; n_bits, tautomeric);
    let params = params.as_deref();
    if keyword_form {
        _fingerprint_pattern_query_keywords(query, n_bits, tautomeric)
    } else {
        _fingerprint_pattern_query_params(query, params)
    }
}
#[pyfunction]
#[doc = r#"{"role":"overload","target":"fingerprint_pattern_query"}"#]
#[pyo3(signature = (query, params=None))]
pub fn _fingerprint_pattern_query_params(
    query: &QueryGraph,
    params: Option<&PatternFingerprintParams>,
) -> PyResult<Fingerprint> {
    let params = pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

    match params {
        Some(params) => ck::fingerprint_pattern_query_with_params(query.inner(), &params),
        None => ck::fingerprint_pattern_query(query.inner()),
    }
    .map(|inner| Fingerprint { inner })
    .map_err(PatternFingerprintError::from_rust)
}
#[pyfunction]
#[doc = r#"{"role":"overload","target":"fingerprint_pattern_query"}"#]
#[pyo3(signature = (query, *, n_bits=None, tautomeric=None))]
pub fn _fingerprint_pattern_query_keywords(
    query: &QueryGraph,
    n_bits: Option<usize>,
    tautomeric: Option<bool>,
) -> PyResult<Fingerprint> {
    let params = keyword_params(n_bits, tautomeric);
    _fingerprint_pattern_query_params(query, Some(&params))
}
#[pyfunction]
pub fn fingerprint_pattern_query_with_params(
    query: &QueryGraph,
    params: &PatternFingerprintParams,
) -> PyResult<Fingerprint> {
    let params = pyo3::Python::attach(|py| params.to_rust(py))?;

    ck::fingerprint_pattern_query_with_params(query.inner(), &params)
        .map(|inner| Fingerprint { inner })
        .map_err(PatternFingerprintError::from_rust)
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
    fn actual_molecule_and_query_pattern_producers_keep_bits_call_forms_and_errors() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let module = PyModule::new(py, "pattern_probe")?;
            module.add_function(wrap_pyfunction!(fingerprint_pattern_query, &module)?)?;
            module.add_function(wrap_pyfunction!(
                fingerprint_pattern_query_with_params,
                &module
            )?)?;
            let query = Py::new(py, QueryGraph::from_rust(ck::parse_smarts("C~N").unwrap()))?;
            let query_before = query.borrow(py).inner().clone();
            let expected_query = ck::fingerprint_pattern_query(query.borrow(py).inner()).unwrap();
            exact(
                py,
                module
                    .getattr("fingerprint_pattern_query")?
                    .call1((query.clone_ref(py),))?
                    .extract()?,
                &expected_query,
            );
            for text in ["CCO", "c1ccccc1", "N=C=O"] {
                for prepared in [true, false] {
                    let native = ck::Molecule::from_smiles_with_params(
                        text,
                        &ck::SmilesParseParams {
                            sanitize: prepared,
                            remove_hs: prepared,
                            ..Default::default()
                        },
                    )
                    .unwrap();
                    let before = native.to_binary().unwrap();
                    let molecule = Py::new(py, Molecule { inner: native })?;
                    let expected = molecule.borrow(py).inner.fingerprint_pattern().unwrap();
                    exact(
                        py,
                        molecule
                            .bind(py)
                            .call_method0("fingerprint_pattern")?
                            .extract()?,
                        &expected,
                    );
                    for n_bits in [128, 2048] {
                        for tautomeric in [false, true] {
                            let params = Py::new(
                                py,
                                PatternFingerprintParams {
                                    inner: crate::parameters::NativeField::Owned(
                                        ck::PatternFingerprintParams { n_bits, tautomeric },
                                    ),
                                },
                            )?;
                            let kwargs = PyDict::new(py);
                            kwargs.set_item("n_bits", n_bits)?;
                            kwargs.set_item("tautomeric", tautomeric)?;
                            let expected = molecule
                                .borrow(py)
                                .inner
                                .fingerprint_pattern_with_params(&params.borrow(py).to_rust(py)?)
                                .unwrap();
                            for actual in [
                                molecule
                                    .bind(py)
                                    .call_method1("fingerprint_pattern", (params.clone_ref(py),))?,
                                molecule.bind(py).call_method(
                                    "fingerprint_pattern",
                                    (),
                                    Some(&kwargs),
                                )?,
                                molecule.bind(py).call_method1(
                                    "fingerprint_pattern_with_params",
                                    (params.clone_ref(py),),
                                )?,
                            ] {
                                exact(py, actual.extract()?, &expected);
                            }
                            let expected = ck::fingerprint_pattern_query_with_params(
                                query.borrow(py).inner(),
                                &params.borrow(py).to_rust(py)?,
                            )
                            .unwrap();
                            for actual in [
                                module
                                    .getattr("fingerprint_pattern_query")?
                                    .call1((query.clone_ref(py), params.clone_ref(py)))?,
                                module
                                    .getattr("fingerprint_pattern_query")?
                                    .call((query.clone_ref(py),), Some(&kwargs))?,
                                module
                                    .getattr("fingerprint_pattern_query_with_params")?
                                    .call1((query.clone_ref(py), params.clone_ref(py)))?,
                            ] {
                                exact(py, actual.extract()?, &expected);
                            }
                            assert!(
                                molecule
                                    .bind(py)
                                    .call_method(
                                        "fingerprint_pattern",
                                        (params.clone_ref(py),),
                                        Some(&kwargs)
                                    )
                                    .is_err()
                            );
                            assert!(
                                module
                                    .getattr("fingerprint_pattern_query")?
                                    .call((query.clone_ref(py), params), Some(&kwargs))
                                    .is_err()
                            );
                        }
                    }
                    let kwargs = PyDict::new(py);
                    kwargs.set_item("n_bits", 0)?;
                    let expected = molecule
                        .borrow(py)
                        .inner
                        .fingerprint_pattern_with_params(&ck::PatternFingerprintParams {
                            n_bits: 0,
                            ..Default::default()
                        })
                        .unwrap_err()
                        .to_string();
                    let error = molecule
                        .bind(py)
                        .call_method("fingerprint_pattern", (), Some(&kwargs))
                        .unwrap_err();
                    assert!(error.is_instance_of::<PatternFingerprintError>(py));
                    assert_eq!(error.value(py).str()?.to_str()?, expected);
                    assert_eq!(molecule.borrow(py).inner.to_binary().unwrap(), before);
                }
            }
            let kwargs = PyDict::new(py);
            kwargs.set_item("n_bits", 0)?;
            let error = module
                .getattr("fingerprint_pattern_query")?
                .call((query.clone_ref(py),), Some(&kwargs))
                .unwrap_err();
            assert!(error.is_instance_of::<PatternFingerprintError>(py));
            assert_eq!(query.borrow(py).inner(), &query_before);
            Ok(())
        })
        .unwrap();
    }
}
