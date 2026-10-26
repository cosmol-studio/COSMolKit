//! Owned standard/external source chains. Known CK errors are projected by owners.
use crate::error_framework::Source;
use crate::utf8_errors::{FromUtf8Error, Utf8Error};
use pyo3::prelude::*;
use std::error::Error;
use std::sync::Arc;

/// Share the actual immutable native root, including non-Clone IO errors.
/// Source access follows native Error::source, never reconstructed messages.
#[derive(Clone)]
pub(crate) struct ErrorRef {
    root: Arc<dyn Error + Send + Sync>,
    projection: Option<Arc<ErrorProjection>>,
    depth: usize,
}
type ErrorSelector = for<'a> fn(&'a (dyn Error + 'static)) -> &'a (dyn Error + 'static);
struct ErrorProjection {
    parent: ErrorRef,
    select: ErrorSelector,
}
impl ErrorRef {
    pub(crate) fn new(root: Arc<dyn Error + Send + Sync>) -> Self {
        Self {
            root,
            projection: None,
            depth: 0,
        }
    }
    pub(crate) fn error(&self) -> &(dyn Error + 'static) {
        let mut value: &(dyn Error + 'static) = match &self.projection {
            Some(projection) => (projection.select)(projection.parent.error()),
            None => self.root.as_ref(),
        };
        for _ in 0..self.depth {
            value = value.source().expect("retained native source edge");
        }
        value
    }
    pub(crate) fn source(&self) -> Option<Self> {
        self.error().source().map(|_| Self {
            root: Arc::clone(&self.root),
            projection: self.projection.clone(),
            depth: self.depth + 1,
        })
    }
    // A transparent error's payload is not necessarily in Error::source().
    // Retain a readonly native field path from the same original owner instead
    // of inventing a cause edge or cloning/reconstructing a non-Clone error.
    // The higher-ranked selector cannot return a borrow outliving its parent.
    // Ordinary source traversal and ErrorRef cloning still copy only Arcs.
    pub(crate) fn project(&self, select: ErrorSelector) -> Self {
        Self {
            root: Arc::clone(&self.root),
            projection: Some(Arc::new(ErrorProjection {
                parent: self.clone(),
                select,
            })),
            depth: 0,
        }
    }
}

// Native ErrorKind is non-exhaustive and includes OS kinds unavailable as
// stable Rust constants. Keep the actual native value even for those kinds;
// named Final constants are convenience values, not a lossy closed enum map.
#[pyclass(module = "cosmolkit", frozen, from_py_object, eq)]
#[derive(Clone, Copy, PartialEq)]
/// Native IO failure category; named constants compare by category rather than object identity.
pub struct IoErrorKind {
    inner: std::io::ErrorKind,
}
impl IoErrorKind {
    pub(crate) fn from_native(inner: std::io::ErrorKind) -> Self {
        Self { inner }
    }
}
macro_rules! io_kinds {
    ($($kind:ident),+ $(,)?) => {
        #[pymethods]
        impl IoErrorKind {
            $(#[classattr]
            #[allow(non_upper_case_globals)]
            const $kind: Self = Self { inner: std::io::ErrorKind::$kind };)+
            fn __repr__(&self) -> String { format!("IoErrorKind.{:?}", self.inner) }
            fn __str__(&self) -> String { self.inner.to_string() }
            fn __hash__(&self) -> u64 { self.inner as u64 }
        }
    };
}
io_kinds! {
    NotFound, PermissionDenied, ConnectionRefused, ConnectionReset, HostUnreachable,
    NetworkUnreachable, ConnectionAborted, NotConnected, AddrInUse, AddrNotAvailable,
    NetworkDown, BrokenPipe, AlreadyExists, WouldBlock, NotADirectory, IsADirectory,
    DirectoryNotEmpty, ReadOnlyFilesystem, StaleNetworkFileHandle, InvalidInput,
    InvalidData, TimedOut, WriteZero, StorageFull, NotSeekable, QuotaExceeded,
    FileTooLarge, ResourceBusy, ExecutableFileBusy, Deadlock, CrossesDevices,
    TooManyLinks, InvalidFilename, ArgumentListTooLong, Interrupted, Unsupported,
    UnexpectedEof, OutOfMemory, Other,
}

pub(crate) type StandardSource = Source<
    crate::error_framework::ErrorObject<IoError>,
    Source<
        crate::error_framework::ErrorObject<Utf8Error>,
        Source<
            crate::error_framework::ErrorObject<FromUtf8Error>,
            crate::error_framework::ErrorObject<ExternalError>,
        >,
    >,
>;

/// Only standard/third-party source branches enter here. A CK owner must use
/// its concrete CK projection before handing off a standard/external branch.
pub(crate) fn standard_source(node: ErrorRef) -> PyErr {
    if node.error().is::<std::io::Error>() {
        return IoError::from_ref(node);
    }
    if let Some(value) = node.error().downcast_ref::<std::str::Utf8Error>() {
        return Utf8Error::from_rust(*value);
    }
    if let Some(value) = node.error().downcast_ref::<std::string::FromUtf8Error>() {
        return FromUtf8Error::from_rust(value.clone());
    }
    ExternalError::from_ref(node)
}

crate::error_framework::error_projection! {record
#[doc = r#"{"role":"error_record"}"#]
pub struct IoError {
    node: ErrorRef,
}
}
impl IoError {
    fn native(&self) -> &std::io::Error {
        self.node.error().downcast_ref().expect("native IO source")
    }
    pub(crate) fn from_ref(node: ErrorRef) -> PyErr {
        Python::attach(|py| {
            let cause = node.source().map(standard_source);
            match Py::new(py, Self { node }) {
                Ok(value) => {
                    let error = crate::error_framework::from_payload(value);
                    error.set_cause(py, cause);
                    error
                }
                Err(error) => error,
            }
        })
    }
}
crate::error_framework::error_projection! {methods impl IoError {
    #[getter]
    /// Native subsystem reporting this failure, preserved independently of the human-readable message.
    fn domain(&self) -> &'static str {
        "std"
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Native IO error category."] fn kind(&self) -> IoErrorKind {
        IoErrorKind {
            inner: self.native().kind(),
        }
    }
    #[getter]
    #[doc = "Original operating-system error number, or None for a non-OS error."] fn raw_os_error(&self) -> Option<i32> {
        self.native().raw_os_error()
    }
    #[staticmethod]
    #[doc = r#"{"role":"error_source"}"#]
    /// Immediate native cause, also available as __cause__; None when the native error has no cause.
    fn source(slf: Bound<'_, pyo3::exceptions::PyBaseException>) -> PyResult<Option<StandardSource>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.native().to_string()
    }
}
}

// An actual open-ended external error, not a replacement for a known CK type.
crate::error_framework::error_projection! {record
#[doc = r#"{"role":"error_record"}"#]
pub struct ExternalError {
    node: ErrorRef,
}
}
impl ExternalError {
    fn from_ref(node: ErrorRef) -> PyErr {
        Python::attach(|py| {
            let cause = node.source().map(standard_source);
            match Py::new(py, Self { node }) {
                Ok(value) => {
                    let error = crate::error_framework::from_payload(value);
                    error.set_cause(py, cause);
                    error
                }
                Err(error) => error,
            }
        })
    }
}
crate::error_framework::error_projection! {methods impl ExternalError {
    #[getter]
    /// Native subsystem reporting this failure, preserved independently of the human-readable message.
    fn domain(&self) -> &'static str {
        "external"
    }
    #[staticmethod]
    #[doc = r#"{"role":"error_source"}"#]
    /// Immediate native cause, also available as __cause__; None when the native error has no cause.
    fn source(slf: Bound<'_, pyo3::exceptions::PyBaseException>) -> PyResult<Option<StandardSource>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.node.error().to_string()
    }
}
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Debug)]
    struct Leaf;
    impl std::fmt::Display for Leaf {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            f.write_str("external leaf\0context")
        }
    }
    impl Error for Leaf {}
    #[derive(Debug)]
    struct Parent(Leaf);
    impl std::fmt::Display for Parent {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            f.write_str("external parent")
        }
    }
    impl Error for Parent {
        fn source(&self) -> Option<&(dyn Error + 'static)> {
            Some(&self.0)
        }
    }

    #[test]
    fn actual_io_kinds_and_sources_keep_nonclone_owners_without_invented_edges() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            // Even an OS kind without a stable Rust constant remains its exact
            // native value, not a guessed Other or a string discriminator.
            let native = std::io::Error::from_raw_os_error(i32::MAX);
            let kind = native.kind();
            let display = native.to_string();
            let root: Arc<dyn Error + Send + Sync> = Arc::new(native);
            let error = IoError::from_ref(ErrorRef::new(root));
            assert_eq!(
                error
                    .value(py)
                    .getattr("kind")?
                    .extract::<IoErrorKind>()?
                    .inner,
                kind
            );
            assert_eq!(
                error.value(py).getattr("raw_os_error")?.extract::<i32>()?,
                i32::MAX
            );
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert!(error.value(py).getattr("source")?.is_none());
            assert!(error.value(py).setattr("kind", py.None()).is_err());

            let native = std::io::Error::new(std::io::ErrorKind::Other, Parent(Leaf));
            // std::io::Error forwards its custom payload's source, not the
            // payload itself. Follow this real edge rather than adding Parent.
            assert!(native.source().unwrap().is::<Leaf>());
            let root: Arc<dyn Error + Send + Sync> = Arc::new(native);
            let owner = Arc::downgrade(&root);
            let error = IoError::from_ref(ErrorRef::new(root));
            let cause = error.value(py).getattr("source")?;
            assert!(
                cause.is_instance(
                    <ExternalError as crate::error_framework::ExceptionPayload>::exception_type(py)
                        .expect("error backend")
                        .as_any()
                )?
            );
            assert!(cause.is(&error.value(py).getattr("__cause__")?));
            assert_eq!(cause.str()?.to_str()?, "external leaf\0context");
            assert!(cause.getattr("source")?.is_none());
            let retained = cause.unbind();
            drop(error);
            assert!(owner.upgrade().is_some());
            assert_eq!(retained.bind(py).str()?.to_str()?, "external leaf\0context");
            drop(retained);
            Ok(())
        })
        .unwrap();
    }
}
