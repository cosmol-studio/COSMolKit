//! MACCS values and error transport over the sole canonical facade.
use ::cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
use std::error::Error;
pyo3::create_exception!(cosmolkit, MaccsFingerprintError, PyValueError);
pub(crate) fn maccs_pyerr(py: Python<'_>, source: ck::MaccsFingerprintError) -> PyErr {
    use ck::MaccsFingerprintError as E;
    let err = MaccsFingerprintError::new_err(source.to_string());
    let attrs = || -> PyResult<()> {
        let object = err.value(py);
        object.setattr("domain", "Fingerprint")?;
        object.setattr(
            "kind",
            match &source {
                E::UnsupportedOption { .. } => "UnsupportedOption",
                E::MissingPattern { .. } => "MissingPattern",
                E::Topology(_) => "Topology",
                E::Rings(_) => "Rings",
                E::Valence(_) => "Valence",
                E::Paths(_) => "Paths",
                E::Smarts(_) => "Smarts",
                E::QueryCompile(_) => "QueryCompile",
                E::QueryContext(_) => "QueryContext",
                E::Match(_) => "Match",
                E::Value(_) => "Value",
            },
        )?;
        match &source {
            E::UnsupportedOption { option, reason } => {
                object.setattr("option", option)?;
                object.setattr("reason", reason)?;
            }
            E::MissingPattern { bit } => object.setattr("bit", bit)?,
            _ => {}
        }
        Ok(())
    };
    if let Err(e) = attrs() {
        return e;
    }
    err.set_cause(
        py,
        source
            .source()
            .map(|e| crate::canonical_values::source_pyerr(py, e)),
    );
    err
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct MaccsFingerprintParams {
    pub(crate) inner: ck::MaccsFingerprintParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl MaccsFingerprintParams {
    #[new]
    #[pyo3(signature=(*,n_bits=166))]
    fn new(n_bits: usize) -> Self {
        Self {
            inner: ck::MaccsFingerprintParams { n_bits },
        }
    }
    #[getter]
    fn n_bits(&self) -> usize {
        self.inner.n_bits
    }
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<MaccsFingerprintParams>()?;
    module.add(
        "MaccsFingerprintError",
        module.py().get_type::<MaccsFingerprintError>(),
    )?;
    Ok(())
}
