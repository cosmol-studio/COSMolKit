//! Registry-driven argument normalization and configuration display; no chemistry.
use pyo3::prelude::*;
use pyo3::types::PyModule;

#[path = "../../crates/cosmolkit/examples/support/binding_contract_manifest.rs"]
mod contract;

pub(crate) fn repr(value: &Bound<'_, PyAny>) -> PyResult<String> {
    value.getattr("_configuration_repr")?.call0()?.extract()
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = module.py();
    let code = std::ffi::CString::new(include_str!("configuration_projection.py"))
        .expect("embedded Python contains no NUL");
    let adapter = PyModule::from_code(
        py,
        &code,
        c"_configuration_projection.py",
        c"_cosmolkit_configuration",
    )?;
    let document = py
        .import("json")?
        .getattr("loads")?
        .call1((contract::manifest().to_string(),))?;
    adapter.getattr("install")?.call1((module, &document))?;
    // Stub generation consumes the very adapter installed on the actual module.
    module.add(
        "_configuration_declarations",
        adapter.getattr("declarations")?,
    )?;
    Ok(())
}
