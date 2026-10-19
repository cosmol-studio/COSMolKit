//! Native detector results; the binding neither sanitizes nor detects chemistry.
use crate::error_framework::error_projection;
use crate::kekulize_errors::{KekulizeError, KekulizeSource};
use crate::sanitize_values::SanitizeStage;
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::error::Error;

error_projection!(ChemistryProblemError, ChemistryProblemErrorKind, ck::ChemistryProblemError, "sanitize", {
    Valence(_), Kekulize(_),
}, cause = |inner: &ck::ChemistryProblemError| Python::attach(|py| match inner {
    ck::ChemistryProblemError::Valence(value) => ValenceError::from_rust(value.clone()).cause(py),
    ck::ChemistryProblemError::Kekulize(value) => KekulizeError::from_rust(value.clone()).cause(py),
}), source = Option<KekulizeSource>, methods = {
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn valence(&self, py: Python<'_>) -> PyResult<Option<Py<ValenceError>>> { match &self.inner {
        ck::ChemistryProblemError::Valence(value) => ValenceError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn kekulize(&self, py: Python<'_>) -> PyResult<Option<Py<KekulizeError>>> { match &self.inner {
        ck::ChemistryProblemError::Kekulize(value) => KekulizeError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
});

#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub struct ChemistryProblem {
    pub(crate) inner: ck::ChemistryProblem,
}
#[pymethods]
impl ChemistryProblem {
    #[new]
    fn new(operation: SanitizeStage, error: PyRef<'_, ChemistryProblemError>) -> Self {
        Self {
            inner: ck::ChemistryProblem {
                operation: operation.into(),
                error: error.inner.clone(),
            },
        }
    }
    #[getter]
    fn operation(&self) -> SanitizeStage {
        self.inner.operation.into()
    }
    #[getter]
    fn error(&self, py: Python<'_>) -> PyResult<Py<ChemistryProblemError>> {
        ChemistryProblemError::from_rust(self.inner.error.clone())
            .value(py)
            .extract()
            .map_err(Into::into)
    }
}

#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub struct ChemistryProblemReport {
    pub(crate) inner: ck::ChemistryProblemReport,
}
#[pymethods]
impl ChemistryProblemReport {
    #[new]
    #[pyo3(signature = (problems=None))]
    fn new(problems: Option<Vec<ChemistryProblem>>) -> Self {
        Self {
            inner: match problems {
                None => ck::ChemistryProblemReport::default(),
                Some(rows) => ck::ChemistryProblemReport {
                    problems: rows.into_iter().map(|row| row.inner).collect(),
                },
            },
        }
    }
    #[getter]
    fn problems(&self) -> Vec<ChemistryProblem> {
        self.inner
            .problems
            .iter()
            .cloned()
            .map(|inner| ChemistryProblem { inner })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sanitize_values::SanitizeParams;

    #[test]
    fn real_detector_reports_keep_native_order_stages_errors_and_call_forms() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            for text in ["CCO", "C(C)(C)(C)(C)C", "c1cccc1"] {
                let native = ck::Molecule::from_smiles_with_params(
                    text,
                    &ck::SmilesParseParams {
                        sanitize: false,
                        remove_hs: false,
                        ..Default::default()
                    },
                )
                .unwrap();
                let before = native.to_binary().unwrap();
                let molecule = Py::new(py, crate::native::Molecule { inner: native })?;
                let expected = molecule
                    .borrow(py)
                    .inner
                    .detect_chemistry_problems()
                    .unwrap();
                let report = molecule
                    .bind(py)
                    .call_method0("detect_chemistry_problems")?
                    .extract::<Py<ChemistryProblemReport>>()?;
                assert_eq!(report.borrow(py).inner, expected);
                let rows = report
                    .bind(py)
                    .getattr("problems")?
                    .extract::<Vec<ChemistryProblem>>()?;
                assert_eq!(rows.len(), expected.problems.len());
                for (row, expected) in rows.iter().zip(&expected.problems) {
                    assert_eq!(row.inner, *expected);
                    let operation = row.operation();
                    assert_eq!(ck::SanitizeStage::from(operation), expected.operation);
                    let error = row.error(py)?;
                    assert_eq!(error.bind(py).str()?.to_str()?, expected.error.to_string());
                    match &expected.error {
                        ck::ChemistryProblemError::Valence(value) => {
                            assert!(
                                error
                                    .bind(py)
                                    .getattr("valence")?
                                    .is_instance_of::<ValenceError>()
                            );
                            assert_eq!(
                                error.bind(py).getattr("valence")?.str()?.to_str()?,
                                value.to_string()
                            );
                        }
                        ck::ChemistryProblemError::Kekulize(value) => {
                            assert!(
                                error
                                    .bind(py)
                                    .getattr("kekulize")?
                                    .is_instance_of::<KekulizeError>()
                            );
                            assert_eq!(
                                error.bind(py).getattr("kekulize")?.str()?.to_str()?,
                                value.to_string()
                            );
                        }
                    }
                    assert!(error.bind(py).setattr("variant", py.None()).is_err());
                }
                assert!(report.bind(py).setattr("problems", py.None()).is_err());
                for operations in [
                    ck::SanitizeOperations::NONE,
                    ck::SanitizeOperations::PROPERTIES,
                    ck::SanitizeOperations::KEKULIZE,
                ] {
                    let params = Py::new(
                        py,
                        SanitizeParams {
                            inner: ck::SanitizeParams { operations },
                        },
                    )?;
                    let expected = molecule
                        .borrow(py)
                        .inner
                        .detect_chemistry_problems_with_params(&params.borrow(py).inner)
                        .unwrap();
                    let kwargs = pyo3::types::PyDict::new(py);
                    kwargs.set_item(
                        "operations",
                        crate::sanitize_values::SanitizeOperations { inner: operations },
                    )?;
                    for actual in [
                        molecule
                            .bind(py)
                            .call_method1("detect_chemistry_problems", (params.clone_ref(py),))?,
                        molecule.bind(py).call_method(
                            "detect_chemistry_problems",
                            (),
                            Some(&kwargs),
                        )?,
                        molecule.bind(py).call_method1(
                            "detect_chemistry_problems_with_params",
                            (params.clone_ref(py),),
                        )?,
                    ] {
                        assert_eq!(
                            actual.extract::<PyRef<'_, ChemistryProblemReport>>()?.inner,
                            expected
                        );
                    }
                    assert!(
                        molecule
                            .bind(py)
                            .call_method("detect_chemistry_problems", (params,), Some(&kwargs))
                            .is_err()
                    );
                }
                assert_eq!(molecule.borrow(py).inner.to_binary().unwrap(), before);
                drop(molecule);
                assert_eq!(report.borrow(py).inner, expected);
            }
            Ok(())
        })
        .unwrap();
    }
}
