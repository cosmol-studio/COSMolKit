//! Typed Python callables with invocation-local exception preservation.
use pyo3::prelude::*;
use std::sync::Mutex;

macro_rules! typed_callback {
    ($name:ident, [$($argument:expr),+ $(,)?]) => {
        typed_callback!($name, [$($argument),+], pyo3::type_hint_identifier!("builtins", "bool"));
    };
    ($name:ident, [$($argument:expr),* $(,)?], $output:expr) => {
        pub struct $name { handle: Py<pyo3::PyAny> }
        impl FromPyObject<'_, '_> for $name {
            type Error = PyErr;
            const INPUT_TYPE: pyo3::inspect::PyStaticExpr = pyo3::type_hint_subscript!(
                pyo3::type_hint_identifier!("collections.abc", "Callable"),
                pyo3::inspect::PyStaticExpr::List { elts: &[$($argument),*] },
                $output
            );
            fn extract(value: pyo3::Borrowed<'_, '_, pyo3::PyAny>) -> PyResult<Self> {
                if !value.is_callable() {
                    return Err(pyo3::exceptions::PyTypeError::new_err("expected a callable"));
                }
                Ok(Self { handle: value.to_owned().unbind() })
            }
        }
        impl<'py> IntoPyObject<'py> for $name {
            type Target = pyo3::PyAny;
            type Output = Bound<'py, pyo3::PyAny>;
            type Error = PyErr;
            const OUTPUT_TYPE: pyo3::inspect::PyStaticExpr = <Self as FromPyObject>::INPUT_TYPE;
            fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
                Ok(self.handle.into_bound(py))
            }
        }
        impl $name {
            pub(crate) fn clone_ref(&self, py: Python<'_>) -> Self {
                Self { handle: self.handle.clone_ref(py) }
            }
            pub(crate) fn bind<'py>(&self, py: Python<'py>) -> &Bound<'py, pyo3::PyAny> {
                self.handle.bind(py)
            }
            pub(crate) fn traverse(&self, visit: &pyo3::gc::PyVisit<'_>) -> Result<(), pyo3::PyTraverseError> {
                visit.call(&self.handle)
            }
        }
    };
}
pub(crate) use typed_callback;

// Rust callbacks return bool, not Result. Keep the original exception in a
// per-call slot, return false to the owner, then restore that very PyErr after
// the owner's call. Reusable configuration has no shared error state.
pub(crate) struct CallbackErrors(Mutex<Option<PyErr>>);
impl CallbackErrors {
    pub(crate) fn new() -> Self {
        Self(Mutex::new(None))
    }
    pub(crate) fn call(&self, run: impl FnOnce(Python<'_>) -> PyResult<bool>) -> bool {
        if self
            .0
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
            .is_some()
        {
            return false;
        }
        match Python::attach(run) {
            Ok(value) => value,
            Err(error) => {
                self.record(error);
                false
            }
        }
    }
    pub(crate) fn notify(&self, run: impl FnOnce(Python<'_>) -> PyResult<()>) {
        // A void progress notification is not a filtering decision: preserve
        // the owner's one-call-per-row policy even after a previous failure.
        if let Err(error) = Python::attach(run) {
            self.record(error);
        }
    }
    pub(crate) fn record(&self, error: PyErr) {
        let mut slot = self
            .0
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        if slot.is_none() {
            *slot = Some(error);
        }
    }
    pub(crate) fn take(&self) -> Option<PyErr> {
        self.0
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
            .take()
    }
}
