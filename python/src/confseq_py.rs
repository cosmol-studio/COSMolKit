use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyModule;

use crate::Molecule;

fn confseq_pyerr(error: cosmolkit_core::ConfSeqDecodeError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

fn parse_errors_mode(errors: &str) -> PyResult<bool> {
    match errors.to_ascii_lowercase().as_str() {
        "raise" | "strict" => Ok(false),
        "keep" => Ok(true),
        _ => Err(PyValueError::new_err("errors must be 'raise' or 'keep'")),
    }
}

fn parse_template_backend(value: &str) -> PyResult<cosmolkit_core::ConfSeqTemplateBackend> {
    match value.to_ascii_lowercase().as_str() {
        "distance_geometry" | "distance-geometry" | "dg" => {
            Ok(cosmolkit_core::ConfSeqTemplateBackend::DistanceGeometry)
        }
        "fast_geometry" | "fast-geometry" => Ok(cosmolkit_core::ConfSeqTemplateBackend::FastGeometry),
        _ => Err(PyValueError::new_err(
            "template_backend must be 'distance_geometry' or 'fast_geometry'",
        )),
    }
}

#[pyfunction]
#[pyo3(signature = (confseq, *, optimize_with_uff = true, template_backend = "distance_geometry"))]
/// Decode one ConfSeq record.
///
/// `template_backend` accepts "distance_geometry" for the reference backend or
/// "fast_geometry" for the fast initial-coordinate backend.
fn decode(confseq: &str, optimize_with_uff: bool, template_backend: &str) -> PyResult<Molecule> {
    let options = cosmolkit_core::ConfSeqDecodeOptions {
        optimize_with_uff,
        template_backend: parse_template_backend(template_backend)?,
        ..cosmolkit_core::ConfSeqDecodeOptions::default()
    };
    let inner = cosmolkit_core::decode_confseq_record_with_options(confseq, &options).map_err(confseq_pyerr)?;
    Ok(Molecule { inner })
}

#[pyfunction]
#[pyo3(signature = (in_smiles, confseq, *, optimize_with_uff = true, template_backend = "distance_geometry"))]
/// Decode one ConfSeq corpus row with its explicit tokenized input SMILES.
///
/// This is intended for corpus parity and diagnostics. Normal callers should
/// use `decode(confseq, ...)`.
fn decode_with_input_smiles(
    in_smiles: &str,
    confseq: &str,
    optimize_with_uff: bool,
    template_backend: &str,
) -> PyResult<Molecule> {
    let options = cosmolkit_core::ConfSeqDecodeOptions {
        optimize_with_uff,
        template_backend: parse_template_backend(template_backend)?,
        ..cosmolkit_core::ConfSeqDecodeOptions::default()
    };
    let inner = cosmolkit_core::decode_confseq_with_options(in_smiles, confseq, &options).map_err(confseq_pyerr)?;
    Ok(Molecule { inner })
}

#[pyfunction]
#[pyo3(signature = (confseq_list, *, errors = "raise", n_jobs = None, optimize_with_uff = true, template_backend = "distance_geometry"))]
/// Decode ConfSeq records in input order.
///
/// `template_backend` accepts "distance_geometry" for the reference backend or
/// "fast_geometry" for the fast initial-coordinate backend.
fn decode_batch(
    confseq_list: Vec<String>,
    errors: &str,
    n_jobs: Option<usize>,
    optimize_with_uff: bool,
    template_backend: &str,
) -> PyResult<Vec<Option<Molecule>>> {
    if matches!(n_jobs, Some(0)) {
        return Err(PyValueError::new_err("n_jobs must be >= 1"));
    }
    let keep_errors = parse_errors_mode(errors)?;
    let mut options = cosmolkit_core::ConfSeqDecodeOptions::default();
    options.num_threads = n_jobs;
    options.optimize_with_uff = optimize_with_uff;
    options.template_backend = parse_template_backend(template_backend)?;
    let result = cosmolkit_core::decode_confseq_record_batch_with_options(&confseq_list, &options, keep_errors)
        .map_err(confseq_pyerr)?;
    Ok(result
        .molecules
        .into_iter()
        .map(|molecule| molecule.map(|inner| Molecule { inner }))
        .collect())
}

#[pyfunction]
#[pyo3(signature = (in_smiles_list, confseq_list, *, errors = "raise", n_jobs = None, optimize_with_uff = true, template_backend = "distance_geometry"))]
/// Decode ConfSeq corpus rows with explicit tokenized input SMILES strings.
///
/// This is intended for corpus parity and diagnostics. Normal callers should
/// use `decode_batch(confseq_list, ...)`.
fn decode_batch_with_input_smiles(
    in_smiles_list: Vec<String>,
    confseq_list: Vec<String>,
    errors: &str,
    n_jobs: Option<usize>,
    optimize_with_uff: bool,
    template_backend: &str,
) -> PyResult<Vec<Option<Molecule>>> {
    if matches!(n_jobs, Some(0)) {
        return Err(PyValueError::new_err("n_jobs must be >= 1"));
    }
    let keep_errors = parse_errors_mode(errors)?;
    let mut options = cosmolkit_core::ConfSeqDecodeOptions::default();
    options.num_threads = n_jobs;
    options.optimize_with_uff = optimize_with_uff;
    options.template_backend = parse_template_backend(template_backend)?;
    let result =
        cosmolkit_core::decode_confseq_batch_with_options(&in_smiles_list, &confseq_list, &options, keep_errors)
            .map_err(confseq_pyerr)?;
    Ok(result
        .molecules
        .into_iter()
        .map(|molecule| molecule.map(|inner| Molecule { inner }))
        .collect())
}

pub(crate) fn add_confseq_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = m.py();
    let submodule = PyModule::new(py, "confseq")?;
    submodule.add_function(wrap_pyfunction!(decode, &submodule)?)?;
    submodule.add_function(wrap_pyfunction!(decode_with_input_smiles, &submodule)?)?;
    submodule.add_function(wrap_pyfunction!(decode_batch, &submodule)?)?;
    submodule.add_function(wrap_pyfunction!(decode_batch_with_input_smiles, &submodule)?)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item("cosmolkit.confseq", &submodule)?;
    m.add_submodule(&submodule)?;
    m.add("confseq", submodule)?;
    Ok(())
}
