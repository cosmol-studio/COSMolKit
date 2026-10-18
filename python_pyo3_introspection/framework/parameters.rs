//! A parameter declaration generates construction and writable fields.
//! It deliberately does not generate domain calls or a second default policy.
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyDictMethods, PyMapping, PyMappingProxy};
use std::collections::BTreeMap;

const STRING_MAPPING: pyo3::inspect::PyStaticExpr = pyo3::type_hint_subscript!(
    pyo3::type_hint_identifier!("collections.abc", "Mapping"),
    pyo3::type_hint_identifier!("builtins", "str"),
    pyo3::type_hint_identifier!("builtins", "str")
);

/// Language conversion only: owned Rust strings in, builtin readonly snapshot out.
pub(crate) struct StringMap(pub(crate) BTreeMap<String, String>);

impl<'a, 'py> FromPyObject<'a, 'py> for StringMap {
    type Error = PyErr;
    const INPUT_TYPE: pyo3::inspect::PyStaticExpr = STRING_MAPPING;

    fn extract(value: pyo3::Borrowed<'a, 'py, pyo3::PyAny>) -> PyResult<Self> {
        let mapping = value.cast::<PyMapping>()?;
        let dictionary = PyDict::new(value.py());
        dictionary.update(&mapping)?;
        dictionary.extract().map(Self)
    }
}

impl<'py> IntoPyObject<'py> for StringMap {
    type Target = PyMappingProxy;
    type Output = Bound<'py, PyMappingProxy>;
    type Error = PyErr;
    const OUTPUT_TYPE: pyo3::inspect::PyStaticExpr = STRING_MAPPING;

    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        let dictionary = self.0.into_pyobject(py)?;
        Ok(PyMappingProxy::new(py, dictionary.as_mapping()))
    }
}

macro_rules! parameters {
    ($name:ident, $module:ident::$rust:ident, {$($field:ident: $ty:ty),+ $(,)?}
        $(, string_maps = [$($map:ident),+ $(,)?])?
        $(, defaults = [$($unexposed:ident),+ $(,)?])?
        $(, methods = {$($method:item)*})?) => {
        #[pyo3::pyclass(module = "fixed_projection._native")]
        pub struct $name { pub(crate) inner: $module::$rust }

        paste::paste! {
            #[pyo3::pymethods]
            impl $name {
                #[new]
                #[pyo3(signature = (*, $($field=$module::$rust::default().$field),+ $(, $($map=None),+)?))]
                fn new($($field: $ty),+ $(, $($map: Option<$crate::parameters::StringMap>),+)?) -> Self {
                    // The exhaustive literal rejects new domain fields at compile time.
                    Self { inner: $module::$rust { $($field,)+
                        $($($map: $map.map(|value| value.0).unwrap_or_else(|| $module::$rust::default().$map),)+)?
                        $($($unexposed: $module::$rust::default().$unexposed,)+)?
                    } }
                }
                $(
                    #[getter]
                    fn $field(&self) -> $ty { self.inner.$field.clone() }
                    #[setter]
                    fn [<set_ $field>](&mut self, value: $ty) { self.inner.$field = value; }
                )+
                $($(
                    #[getter]
                    fn $map(&self) -> $crate::parameters::StringMap {
                        $crate::parameters::StringMap(self.inner.$map.clone())
                    }
                    #[setter]
                    fn [<set_ $map>](&mut self, value: Option<$crate::parameters::StringMap>) {
                        self.inner.$map = value.map(|value| value.0).unwrap_or_else(|| $module::$rust::default().$map);
                    }
                )+)?
                $($($method)*)?
            }
        }
    };
}
pub(crate) use parameters;
