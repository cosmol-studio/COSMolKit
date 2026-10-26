//! Shared exception projection; domain errors stay owned and typed in Rust.
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use pyo3::type_hint_union;
use std::error::Error;

use pyo3::exceptions::PyBaseException;
use pyo3::sync::PyOnceLock;
use pyo3::types::PyType;

/// The declaration associates a native payload with its exception backend once.
pub(crate) trait ExceptionPayload: pyo3::PyClass {
    fn exception_type(py: Python<'_>) -> PyResult<Bound<'_, PyType>>;
}

pub(crate) fn install_backend<T: ExceptionPayload>(
    py: Python<'_>,
    exception: &Bound<'_, PyType>,
) -> PyResult<()> {
    static BACKEND: PyOnceLock<Py<PyModule>> = PyOnceLock::new();
    let backend = BACKEND.get_or_try_init(py, || {
        PyModule::from_code(
            py,
            c_str_backend(),
            c"error_backend.py",
            c"_ck_error_backend",
        )
        .map(Bound::unbind)
    })?;
    backend
        .bind(py)
        .getattr("install")?
        .call1((exception, py.get_type::<T>()))?;
    Ok(())
}

fn c_str_backend() -> &'static std::ffi::CStr {
    // The included framework source is fixed; no runtime code generation.
    std::ffi::CStr::from_bytes_with_nul(concat!(include_str!("error_backend.py"), "\0").as_bytes())
        .expect("error backend source contains no NUL")
}

pub(crate) fn from_payload<T: ExceptionPayload>(value: Py<T>) -> PyErr {
    Python::attach(|py| {
        let result = (|| -> PyResult<PyErr> {
            // Only native projection can construct a payload-bearing error,
            // just as the original pyclass had no public #[new] constructor.
            let exception = py
                .get_type::<pyo3::exceptions::PyValueError>()
                .getattr("__new__")?
                .call1((T::exception_type(py)?,))?;
            exception.setattr("_ck_error_payload", value)?;
            Ok(PyErr::from_value(exception))
        })();
        match result {
            Ok(error) | Err(error) => error,
        }
    })
}

/// A typed reference to the public exception, not to its hidden Rust payload.
pub(crate) struct ErrorObject<T: ExceptionPayload> {
    value: Py<PyBaseException>,
    payload: Py<T>,
}
impl<T: ExceptionPayload> ErrorObject<T> {
    pub(crate) fn bind<'py>(&self, py: Python<'py>) -> &Bound<'py, PyBaseException> {
        self.value.bind(py)
    }
    pub(crate) fn borrow<'py>(&'py self, py: Python<'py>) -> PyRef<'py, T> {
        self.payload.borrow(py)
    }
}
impl<'a, 'py, T: ExceptionPayload> FromPyObject<'a, 'py> for ErrorObject<T> {
    type Error = PyErr;
    const INPUT_TYPE: pyo3::inspect::PyStaticExpr = <Py<T> as FromPyObject>::INPUT_TYPE;
    fn extract(value: pyo3::Borrowed<'a, 'py, PyAny>) -> PyResult<Self> {
        if !value.is_instance(T::exception_type(value.py())?.as_any())? {
            return Err(pyo3::exceptions::PyTypeError::new_err(format!(
                "expected {}",
                T::NAME
            )));
        }
        Ok(Self {
            value: value.cast::<PyBaseException>()?.to_owned().unbind(),
            payload: value.getattr("_ck_error_payload")?.extract()?,
        })
    }
}
impl<'py, T: ExceptionPayload> IntoPyObject<'py> for ErrorObject<T> {
    type Target = PyBaseException;
    type Output = Bound<'py, PyBaseException>;
    type Error = std::convert::Infallible;
    const OUTPUT_TYPE: pyo3::inspect::PyStaticExpr = <Py<T> as IntoPyObject>::OUTPUT_TYPE;
    fn into_pyobject(self, py: Python<'py>) -> Result<Self::Output, Self::Error> {
        Ok(self.value.into_bound(py))
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    // Only explicitly marked declarations participate; no error-name inference.
    for (name, value) in module.dict().iter() {
        if value.hasattr("_exception_type")? {
            let factory = value.getattr("_exception_type")?;
            if factory.getattr("__doc__")?.extract::<String>()? == r#"{"role":"error_projection"}"# {
                module.dict().set_item(name, factory.call0()?)?;
            }
        }
    }
    Ok(())
}

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
    (record $(#[$attr:meta])* $vis:vis struct $name:ident { $($fields:tt)* }) => {
        paste::paste! {
            pyo3::create_exception!(cosmolkit, [<$name Exception>], pyo3::exceptions::PyValueError);
            $(#[$attr])*
            #[doc = concat!("\n\n", stringify!($name), " is a structured native exception. ",
                "Its readonly fields retain the Rust failure data. Use str(error) for the ",
                "message and __cause__ for the underlying failure, when present. ",
                "Instances are raised by operations, not constructed directly.")]
            #[pyclass(module = "cosmolkit")]
            $vis struct $name { $($fields)* }
            impl $crate::error_framework::ExceptionPayload for $name {
                fn exception_type(py: Python<'_>) -> PyResult<Bound<'_, pyo3::types::PyType>> {
                    static READY: pyo3::sync::PyOnceLock<()> = pyo3::sync::PyOnceLock::new();
                    let class = py.get_type::<[<$name Exception>]>();
                    READY.get_or_try_init(py, || $crate::error_framework::install_backend::<Self>(py, &class))?;
                    Ok(class)
                }
            }
        }
    };
    (methods impl $name:ident { $($methods:item)* }) => {
        #[pymethods]
        impl $name {
            #[staticmethod]
            #[doc = r#"{"role":"error_projection"}"#]
            fn _exception_type(py: Python<'_>) -> PyResult<Bound<'_, pyo3::types::PyType>> {
                <Self as $crate::error_framework::ExceptionPayload>::exception_type(py)
            }
            $($methods)*
        }
    };
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
        $crate::enums::binding_enum!($kind, { $($variant),+ });

        $crate::error_framework::error_projection! {record
        #[doc = concat!("Failure in the ", $domain, " domain.\n\n",
            "The readonly variant identifies the Rust error case. Payload fields are present only ",
            "for their documented variant and otherwise return None. str(error) preserves the ",
            "native message; __cause__ retains the underlying exception when one exists.")]
        pub struct $name {
            inner: $crate::error_framework::error_storage_type!($mode, $module::$rust),
        }
        }

        $crate::error_framework::error_projection! {methods impl $name {
            /// Domain that reported this failure; this does not replace the typed variant.
            #[getter]
            fn domain(&self) -> &'static str { $domain }
            /// Exact Rust error variant, independent of the human-readable message.
            #[getter]
            fn variant(&self) -> $kind {
                // No wildcard variants or omitted struct fields: additions fail compilation.
                match $crate::error_framework::native_error_ref!($mode, self.inner) {
                    $($crate::error_framework::variant_pattern!($module::$rust, $variant $({$($field),+})? $(($source))?) => $kind::$variant),+
                }
            }
            /// Message of the immediate underlying cause, or None when there is no cause.
            #[getter]
            fn source_message(&self) -> Option<String> {
                self.inner.source().map(ToString::to_string)
            }
            $(#[staticmethod]
            #[doc = r#"{"role":"error_source"}"#]
            /// The actual underlying exception, identical to __cause__, or None.
            fn source<'py>(slf: Bound<'py, pyo3::exceptions::PyBaseException>) -> PyResult<$source_type> {
                // Read the actual CPython cause, not a separately reconstructed
                // wrapper whose identity/state could differ from __cause__.
                Ok(slf.getattr("__cause__")?.extract()?)
            })?
            $($($($(#[doc = concat!("Readonly `", stringify!($field), "` payload of the `",
                stringify!($variant), "` error variant; None for other variants.")]
            #[getter]
            fn $field(&self) -> Option<$ty> {
                match $crate::error_framework::native_error_ref!($mode, self.inner) {
                    $module::$rust::$variant { $field, .. } => Some($crate::error_framework::payload_value!($field $(, $convert)?)),
                    _ => None,
                }
            })?)+)?) +
            $($($method)*)?
            fn __str__(&self) -> String { self.inner.to_string() }
        }
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
                            let error = $crate::error_framework::from_payload(value);
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
