//! Shared exception projection; domain errors stay owned and typed in Rust.
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use std::error::Error;

pub(crate) fn cause_chain(py: Python<'_>, source: &(dyn Error + 'static)) -> PyErr {
    let error = PyRuntimeError::new_err(source.to_string());
    error.set_cause(py, source.source().map(|next| cause_chain(py, next)));
    error
}

macro_rules! variant_pattern {
    ($module:ident::$rust:ident, $variant:ident {$($field:ident),+}) => {
        $module::$rust::$variant { $($field: _),+ }
    };
    ($module:ident::$rust:ident, $variant:ident ($source:tt)) => {
        $module::$rust::$variant($source)
    };
    ($module:ident::$rust:ident, $variant:ident) => { $module::$rust::$variant };
}

macro_rules! error_projection {
    ($name:ident, $kind:ident, $module:ident::$rust:ident, $domain:literal, {
        $($variant:ident $({$($field:ident: $ty:ty),+ $(,)?})? $(($source:tt))?),+ $(,)?
    }) => {
        #[pyclass(module = "fixed_projection._native", eq, eq_int, from_py_object)]
        #[derive(Clone, Copy, PartialEq)]
        pub enum $kind { $($variant),+ }

        #[pyclass(module = "fixed_projection._native", extends=PyValueError)]
        pub struct $name { inner: $module::$rust }

        #[pymethods]
        impl $name {
            #[getter]
            fn domain(&self) -> &'static str { $domain }
            #[getter]
            fn kind(&self) -> $kind {
                // No wildcard variants or omitted struct fields: additions fail compilation.
                match &self.inner {
                    $(variant_pattern!($module::$rust, $variant $({$($field),+})? $(($source))?) => $kind::$variant),+
                }
            }
            #[getter]
            fn source_message(&self) -> Option<String> {
                self.inner.source().map(ToString::to_string)
            }
            $($($(#[getter]
            fn $field(&self) -> Option<$ty> {
                match &self.inner {
                    $module::$rust::$variant { $field, .. } => Some(*$field),
                    _ => None,
                }
            })+)?) +
            fn __str__(&self) -> String { self.inner.to_string() }
        }

        impl $name {
            pub fn from_rust(inner: $module::$rust) -> PyErr {
                Python::attach(|py| {
                    let cause = inner.source().map(|source| $crate::error_framework::cause_chain(py, source));
                    match Py::new(py, Self { inner }) {
                        Ok(value) => {
                            let error = PyErr::from_value(value.into_bound(py).into_any());
                            error.set_cause(py, cause);
                            error
                        }
                        Err(error) => error,
                    }
                })
            }
        }
    }
}

pub(crate) use {error_projection, variant_pattern};
