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

#[pyfunction]
#[pyo3(signature = (in_smiles, confseq, *, optimize_with_uff = true))]
fn decode(in_smiles: &str, confseq: &str, optimize_with_uff: bool) -> PyResult<Molecule> {
    let options = cosmolkit_core::ConfSeqDecodeOptions {
        optimize_with_uff,
        ..cosmolkit_core::ConfSeqDecodeOptions::default()
    };
    let inner = cosmolkit_core::decode_confseq_with_options(in_smiles, confseq, &options)
        .map_err(confseq_pyerr)?;
    Ok(Molecule { inner })
}

#[pyfunction]
#[pyo3(signature = (in_smiles_list, confseq_list, *, errors = "raise", n_jobs = None, optimize_with_uff = true))]
fn decode_batch(
    in_smiles_list: Vec<String>,
    confseq_list: Vec<String>,
    errors: &str,
    n_jobs: Option<usize>,
    optimize_with_uff: bool,
) -> PyResult<Vec<Option<Molecule>>> {
    if matches!(n_jobs, Some(0)) {
        return Err(PyValueError::new_err("n_jobs must be >= 1"));
    }
    let keep_errors = parse_errors_mode(errors)?;
    let mut options = cosmolkit_core::ConfSeqDecodeOptions::default();
    options.num_threads = n_jobs;
    options.optimize_with_uff = optimize_with_uff;
    let result = cosmolkit_core::decode_confseq_batch_with_options(
        &in_smiles_list,
        &confseq_list,
        &options,
        keep_errors,
    )
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
    submodule.add_function(wrap_pyfunction!(decode_batch, &submodule)?)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item("cosmolkit.confseq", &submodule)?;
    m.add_submodule(&submodule)?;
    m.add("confseq", submodule)?;
    Ok(())
}
