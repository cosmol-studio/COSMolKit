//! Preserve call-site presence independently of an optional argument's value.
use pyo3::prelude::*;

pub(crate) enum Argument<T> {
    Missing,
    Value(T),
}

impl<T> Default for Argument<T> {
    fn default() -> Self {
        Self::Missing
    }
}

impl<T> Argument<T> {
    pub(crate) fn is_supplied(&self) -> bool {
        matches!(self, Self::Value(_))
    }
}

impl<T: Default> Argument<T> {
    pub(crate) fn value_or_default(self) -> T {
        match self {
            Self::Missing => T::default(),
            Self::Value(value) => value,
        }
    }
}

// Check presence before unwrapping optional values. Independent call arguments
// (for example, fingerprint output buffers) do not belong in this keyword list.
macro_rules! configuration_arguments {
    ($params:ident; $($keyword:ident),+ $(,)?) => {
        if $params.is_supplied() && (false $(|| $keyword.is_supplied())+) {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "parameter object and configuration keywords are mutually exclusive",
            ));
        }
        let $params = $params.value_or_default();
        $(let $keyword = $keyword.value_or_default();)+
    };
}
pub(crate) use configuration_arguments;

impl<'a, 'py, T: FromPyObject<'a, 'py>> FromPyObject<'a, 'py> for Argument<T> {
    type Error = PyErr;
    const INPUT_TYPE: pyo3::inspect::PyStaticExpr = T::INPUT_TYPE;

    fn extract(value: pyo3::Borrowed<'a, 'py, pyo3::PyAny>) -> PyResult<Self> {
        T::extract(value).map(Self::Value).map_err(Into::into)
    }
}
