//! Thin path-score/explanation projections through the sole canonical facade.
use ::cosmolkit as ck;
use pyo3::exceptions::{PyIndexError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use std::error::Error;
pyo3::create_exception!(
    cosmolkit,
    TopologicalTorsionPathScoreError,
    PyValueError,
    "A topological-torsion atom path could not be encoded or scored."
);

pub(crate) fn score_pyerr(py: Python<'_>, source: ck::TopologicalTorsionPathScoreError) -> PyErr {
    use ck::TopologicalTorsionPathScoreError as E;
    let err = if matches!(source, E::AtomIndexOutOfRange { .. }) {
        PyIndexError::new_err(source.to_string())
    } else {
        TopologicalTorsionPathScoreError::new_err(source.to_string())
    };
    let kind = match &source {
        E::ZeroSize => "ZeroSize",
        E::ShortPath { .. } => "ShortPath",
        E::ShortAtomCodes { .. } => "ShortAtomCodes",
        E::AtomIndexOutOfRange { .. } => "AtomIndexOutOfRange",
        E::AtomCodeUnderflow { .. } => "AtomCodeUnderflow",
        E::InvalidTopology(_) => "InvalidTopology",
        E::AtomCode(_) => "AtomCode",
        E::PackedCode(_) => "PackedCode",
    };
    let object = err.value(py);
    let assign = || -> PyResult<()> {
        object.setattr("domain", "Fingerprint")?;
        object.setattr("kind", kind)?;
        match &source {
            E::ShortPath { actual, required } | E::ShortAtomCodes { actual, required } => {
                object.setattr("actual", actual)?;
                object.setattr("required", required)?;
            }
            E::AtomIndexOutOfRange { index, atom_count } => {
                object.setattr("index", index)?;
                object.setattr("atom_count", atom_count)?;
            }
            E::AtomCodeUnderflow {
                index,
                code,
                subtract,
            } => {
                object.setattr("index", index)?;
                object.setattr("code", code)?;
                object.setattr("subtract", subtract)?;
            }
            _ => {}
        }
        Ok(())
    };
    if let Err(attribute_error) = assign() {
        return attribute_error;
    }
    err.set_cause(
        py,
        source
            .source()
            .map(|e| crate::canonical_values::source_pyerr(py, e)),
    );
    err
}
/// Complete source decoding, including size zero and zero chunks after u64 ends.
#[cfg_attr(feature = "stubgen", pyo3_stub_gen::derive::gen_stub_pyfunction)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pyfunction]
#[gen_stub(override_return_type(type_repr="tuple[tuple[builtins.str, builtins.int, builtins.int], ...]",imports=("builtins")))]
#[pyo3(signature=(score,size=4))]
fn explain_path_score<'py>(
    py: Python<'py>,
    score: u64,
    size: usize,
) -> PyResult<Bound<'py, PyTuple>> {
    PyTuple::new(py, ck::explain_path_score(score, size))
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add(
        "TopologicalTorsionPathScoreError",
        module.py().get_type::<TopologicalTorsionPathScoreError>(),
    )?;
    module.add_function(wrap_pyfunction!(explain_path_score, module)?)?;
    Ok(())
}
