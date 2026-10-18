//! Shared exception projection; domain errors stay owned and typed in Rust.
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use pyo3::type_hint_union;
use std::error::Error;

/// A concrete native source union, retaining official annotations without
/// exposing an opaque exception base as the declared source type.
pub enum Source<A, B> {
    First(A),
    Second(B),
}

impl<'a, 'py, A: FromPyObject<'a, 'py>, B: FromPyObject<'a, 'py>> FromPyObject<'a, 'py>
    for Source<A, B>
{
    type Error = PyErr;
    const INPUT_TYPE: pyo3::inspect::PyStaticExpr =
        pyo3::type_hint_union!(A::INPUT_TYPE, B::INPUT_TYPE);

    fn extract(value: pyo3::Borrowed<'a, 'py, pyo3::PyAny>) -> PyResult<Self> {
        use pyo3::exceptions::PyTypeError;
        use pyo3::impl_::frompyobject::failed_to_extract_enum;

        // Match pinned PyO3's FromPyObject derive, including branch order,
        // tuple-newtype field context and combined extraction errors. Only its
        // generic INPUT_TYPE fallback (Incomplete) is replaced above.
        let field_error = |name, cause| {
            let error = PyTypeError::new_err(format!("failed to extract field Source::{name}.0"));
            error.set_cause(value.py(), Some(cause));
            error
        };
        let first = match A::extract(value) {
            Ok(value) => return Ok(Self::First(value)),
            Err(error) => field_error("First", error.into()),
        };
        let second = match B::extract(value) {
            Ok(value) => return Ok(Self::Second(value)),
            Err(error) => field_error("Second", error.into()),
        };
        Err(failed_to_extract_enum(
            value.py(),
            "Source",
            &["First", "Second"],
            &["First", "Second"],
            &[first, second],
        ))
    }
}

impl<'py, A: IntoPyObject<'py>, B: IntoPyObject<'py>> IntoPyObject<'py> for Source<A, B> {
    type Target = pyo3::types::PyAny;
    type Output = Bound<'py, Self::Target>;
    type Error = PyErr;
    const OUTPUT_TYPE: pyo3::inspect::PyStaticExpr =
        pyo3::type_hint_union!(A::OUTPUT_TYPE, B::OUTPUT_TYPE);

    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        use pyo3::BoundObject;
        match self {
            Self::First(value) => value
                .into_pyobject(py)
                .map(BoundObject::into_bound)
                .map(Bound::into_any)
                .map_err(Into::into),
            Self::Second(value) => value
                .into_pyobject(py)
                .map(BoundObject::into_bound)
                .map(Bound::into_any)
                .map_err(Into::into),
        }
    }
}

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
    (retained $name:ident, $($rest:tt)*) => {
        $crate::error_framework::error_projection!(@emit [retained] $name, $($rest)*);
    };
    ($name:ident, $($rest:tt)*) => {
        $crate::error_framework::error_projection!(@emit [] $name, $($rest)*);
    };
    (@emit $mode:tt $name:ident, $kind:ident, $module:ident::$rust:ident, $domain:literal, {
        $($variant:ident $({$($field:ident $(: $ty:ty $(=> $convert:expr)?)?),+ $(,)?})? $(($source:tt))?),+ $(,)?
    } $(, cause = $cause:expr)? $(, source = $source_type:ty)? $(, methods = {$($method:item)*})? ) => {
        $crate::error_framework::require_retained_cause!($mode $(, $cause)?);
        #[pyclass(module = "fixed_projection._native", eq, eq_int, from_py_object)]
        #[derive(Clone, Copy, PartialEq)]
        pub enum $kind { $($variant),+ }

        #[pyclass(module = "fixed_projection._native", extends=PyValueError)]
        pub struct $name {
            inner: $crate::error_framework::error_storage_type!($mode, $module::$rust)
        }

        #[pymethods]
        impl $name {
            #[getter]
            fn domain(&self) -> &'static str { $domain }
            #[getter]
            fn variant(&self) -> $kind {
                // No wildcard variants or omitted struct fields: additions fail compilation.
                match $crate::error_framework::native_error_ref!($mode, self.inner) {
                    $($crate::error_framework::variant_pattern!($module::$rust, $variant $({$($field),+})? $(($source))?) => $kind::$variant),+
                }
            }
            #[getter]
            fn source_message(&self) -> Option<String> {
                self.inner.source().map(ToString::to_string)
            }
            $(#[getter]
            #[doc = r#"{"role":"error_source"}"#]
            fn source<'py>(slf: Bound<'py, Self>) -> PyResult<$source_type> {
                // Read the actual CPython cause, not a separately reconstructed
                // wrapper whose identity/state could differ from __cause__.
                Ok(slf.getattr("__cause__")?.extract()?)
            })?
            $($($($(#[getter]
            fn $field(&self) -> Option<$ty> {
                match $crate::error_framework::native_error_ref!($mode, self.inner) {
                    $module::$rust::$variant { $field, .. } => Some($crate::error_framework::payload_value!($field $(, $convert)?)),
                    _ => None,
                }
            })?)+)?) +
            $($($method)*)?
            fn __str__(&self) -> String { self.inner.to_string() }
        }

        impl $name {
            pub fn from_rust(inner: $module::$rust) -> PyErr {
                Self::from_stored($crate::error_framework::error_store!($mode, inner))
            }
            $crate::error_framework::error_from_ref!($mode);
            fn from_stored(inner: $crate::error_framework::error_storage_type!($mode, $module::$rust)) -> PyErr {
                Python::attach(|py| {
                    let cause = $crate::error_framework::project_cause!(inner, py $(, $cause)?);
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

// Storage is the only optional difference. Ordinary projections still own T;
// retained projections reuse the one actual Arc-backed ErrorRef mechanism.
macro_rules! error_storage_type {
    ([], $rust:ty) => { $rust };
    ([retained], $rust:ty) => { $crate::retained_error::RetainedError<$rust> };
}
macro_rules! error_store {
    ([], $inner:ident) => {
        $inner
    };
    ([retained], $inner:ident) => {
        $crate::retained_error::RetainedError::new($inner)
    };
}
macro_rules! native_error_ref {
    ([], $inner:expr) => {
        &($inner)
    };
    ([retained], $inner:expr) => {
        ($inner).as_ref()
    };
}
macro_rules! error_from_ref {
    ([]) => {};
    ([retained]) => {
        pub(crate) fn from_ref(node: $crate::io_source_errors::ErrorRef) -> pyo3::PyErr {
            Self::from_stored($crate::retained_error::RetainedError::from_ref(node))
        }
    };
}
macro_rules! require_retained_cause {
    ([] $(, $cause:expr)?) => {};
    ([retained], $cause:expr) => {};
    ([retained]) => {
        compile_error!(
            "retained error projections require explicit typed native source conversion"
        );
    };
}

macro_rules! payload_value {
    ($field:ident) => {
        $field.clone()
    };
    ($field:ident, $convert:expr) => {
        ($convert)($field)
    };
}
macro_rules! project_cause {
    ($inner:ident, $py:ident) => {
        $inner
            .source()
            .map(|source| $crate::error_framework::cause_chain($py, source))
    };
    ($inner:ident, $py:ident, $convert:expr) => {
        ($convert)(&$inner)
    };
}
pub(crate) use {
    error_from_ref, error_projection, error_storage_type, error_store, native_error_ref,
    payload_value, project_cause, require_retained_cause, variant_pattern,
};
