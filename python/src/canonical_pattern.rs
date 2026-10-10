//! Pattern values and errors project the sole canonical Rust facade.
use ::cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
use std::error::Error;
pyo3::create_exception!(
    cosmolkit,
    PatternFingerprintError,
    PyValueError,
    "Pattern fingerprint generation failed for the molecule or query."
);
pub(crate) fn pattern_pyerr(
    py: Python<'_>,
    source: impl std::borrow::Borrow<ck::PatternFingerprintError>,
) -> PyErr {
    let source = source.borrow();
    use ck::PatternFingerprintError as E;
    let err = PatternFingerprintError::new_err(source.to_string());
    let attrs = || -> PyResult<()> {
        let object = err.value(py);
        object.setattr("domain", "Fingerprint")?;
        object.setattr(
            "kind",
            match source {
                E::EmptyFingerprint => "EmptyFingerprint",
                E::InvalidArguments { .. } => "InvalidArguments",
                E::BitLengthMismatch { .. } => "BitLengthMismatch",
                E::Topology(_) => "Topology",
                E::Query(_) => "Query",
                E::QueryCarrier(_) => "QueryCarrier",
                E::Rings(_) => "Rings",
                E::Smarts(_) => "Smarts",
                E::QueryCompile(_) => "QueryCompile",
                E::QueryContext(_) => "QueryContext",
                E::Match(_) => "Match",
                E::Value(_) => "Value",
            },
        )?;
        match source {
            E::InvalidArguments { reason } => object.setattr("reason", reason)?,
            E::BitLengthMismatch { left, right } => {
                object.setattr("left", left)?;
                object.setattr("right", right)?;
            }
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
/// Writable configuration for Pattern fingerprint generation.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct PatternFingerprintParams {
    pub(crate) inner: ck::PatternFingerprintParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl PatternFingerprintParams {
    /// Configure Pattern fingerprint generation; omitted fields use the defaults shown in the signature.
    #[new]
    #[pyo3(signature=(*,n_bits=2048,tautomeric=false))]
    fn new(n_bits: usize, tautomeric: bool) -> Self {
        Self {
            inner: ck::PatternFingerprintParams { n_bits, tautomeric },
        }
    }
    /// Logical width of the fingerprint in bits.
    #[getter]
    fn n_bits(&self) -> usize {
        self.inner.n_bits
    }
    /// Whether single, double, and aromatic bonds use tautomer-aware hashing.
    #[getter]
    fn tautomeric(&self) -> bool {
        self.inner.tautomeric
    }
}
/// Compute pattern fixed-width bit fingerprints for a QueryGraph without discarding its predicates.
#[cfg_attr(feature = "stubgen", pyo3_stub_gen::derive::gen_stub_pyfunction)]
#[pyfunction]
fn fingerprint_pattern_query(
    py: Python<'_>,
    query: &crate::canonical_search::QueryGraph,
) -> PyResult<crate::canonical_values::Fingerprint> {
    ck::fingerprint_pattern_query(&query.inner)
        .map(|inner| crate::canonical_values::Fingerprint { inner })
        .map_err(|e| pattern_pyerr(py, e))
}
/// Compute pattern fixed-width bit fingerprints for a QueryGraph without discarding its predicates.
#[cfg_attr(feature = "stubgen", pyo3_stub_gen::derive::gen_stub_pyfunction)]
#[pyfunction]
fn fingerprint_pattern_query_with_params(
    py: Python<'_>,
    query: &crate::canonical_search::QueryGraph,
    params: &PatternFingerprintParams,
) -> PyResult<crate::canonical_values::Fingerprint> {
    ck::fingerprint_pattern_query_with_params(&query.inner, &params.inner)
        .map(|inner| crate::canonical_values::Fingerprint { inner })
        .map_err(|e| pattern_pyerr(py, e))
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PatternFingerprintParams>()?;
    module.add(
        "PatternFingerprintError",
        module.py().get_type::<PatternFingerprintError>(),
    )?;
    module.add_function(wrap_pyfunction!(fingerprint_pattern_query, module)?)?;
    module.add_function(wrap_pyfunction!(
        fingerprint_pattern_query_with_params,
        module
    )?)?;
    Ok(())
}
