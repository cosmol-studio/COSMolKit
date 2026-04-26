use pyo3::prelude::*;

#[pyfunction]
fn placeholder() -> &'static str {
    "COSMolKit PyO3 placeholder module"
}

#[pyfunction]
fn rust_version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[pymodule]
fn cosmolkit(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_function(wrap_pyfunction!(placeholder, m)?)?;
    m.add_function(wrap_pyfunction!(rust_version, m)?)?;
    Ok(())
}
