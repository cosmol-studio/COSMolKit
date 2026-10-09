//! Configuration and typed errors; Avalon chemistry stays in the Rust facade.
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

pyo3::create_exception!(cosmolkit, AvalonFingerprintError, PyValueError);
pyo3::create_exception!(cosmolkit, AvalonEngineError, PyValueError);

pub(crate) fn error(py: Python<'_>, source: ck::AvalonFingerprintError) -> PyErr {
    let err = AvalonFingerprintError::new_err(source.to_string());
    let attrs = || -> PyResult<PyErr> {
        let (kind, cause) = match source {
            ck::AvalonFingerprintError::Input(e) => {
                ("Input", crate::canonical_molecular_io::error_pyerr(py, e))
            }
            ck::AvalonFingerprintError::Engine(ref e) => {
                let (kind, reason) = match e {
                    ck::AvalonEngineError::InvalidArguments { reason } => {
                        ("InvalidArguments", *reason)
                    }
                    ck::AvalonEngineError::AvalonConversion { reason } => {
                        ("Conversion", reason.as_str())
                    }
                };
                let cause = AvalonEngineError::new_err(e.to_string());
                cause.value(py).setattr("domain", "Avalon")?;
                cause.value(py).setattr("kind", kind)?;
                cause.value(py).setattr("reason", reason)?;
                err.value(py).setattr("reason", reason)?;
                (kind, cause)
            }
        };
        err.value(py).setattr("domain", "Fingerprint")?;
        err.value(py).setattr("kind", kind)?;
        Ok(cause)
    };
    match attrs() {
        Ok(cause) => {
            err.set_cause(py, Some(cause));
            err
        }
        Err(e) => e,
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct AvalonFingerprintFlags {
    inner: ck::AvalonFingerprintFlags,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl AvalonFingerprintFlags {
    #[staticmethod]
    fn from_bits_retain(bits: u32) -> Self {
        Self {
            inner: ck::AvalonFingerprintFlags::from_bits_retain(bits),
        }
    }
    fn bits(&self) -> u32 {
        self.inner.bits()
    }
    fn __repr__(&self) -> String {
        format!("AvalonFingerprintFlags({:#x})", self.inner.bits())
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct AvalonFingerprintParams {
    pub(crate) inner: ck::AvalonFingerprintParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl AvalonFingerprintParams {
    #[new]
    #[pyo3(signature = (*, n_bits=512, is_query=false, bit_flags=32767))]
    fn new(n_bits: u32, is_query: bool, bit_flags: u32) -> Self {
        Self {
            inner: ck::AvalonFingerprintParams {
                n_bits,
                is_query,
                bit_flags: ck::AvalonFingerprintFlags::from_bits_retain(bit_flags),
            },
        }
    }
    #[getter]
    fn n_bits(&self) -> u32 {
        self.inner.n_bits
    }
    #[getter]
    fn is_query(&self) -> bool {
        self.inner.is_query
    }
    #[getter]
    fn bit_flags(&self) -> u32 {
        self.inner.bit_flags.bits()
    }
    fn __repr__(&self) -> String {
        format!(
            "AvalonFingerprintParams(n_bits={}, is_query={}, bit_flags={:#x})",
            self.inner.n_bits,
            self.inner.is_query,
            self.inner.bit_flags.bits()
        )
    }
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<AvalonFingerprintFlags>()?;
    module.add_class::<AvalonFingerprintParams>()?;
    module.add(
        "AvalonFingerprintError",
        module.py().get_type::<AvalonFingerprintError>(),
    )?;
    module.add(
        "AvalonEngineError",
        module.py().get_type::<AvalonEngineError>(),
    )?;
    Ok(())
}
