//! Native Python transport for the existing detached PropertyValue vocabulary.
use cosmolkit as ck;
use pyo3::{
    exceptions::{PyOverflowError, PyTypeError},
    prelude::*,
    types::{PyBool, PyFloat, PyInt, PyList, PyString},
};

pub(crate) struct NativeProperty(pub ck::PropertyValue);

impl<'a, 'py> FromPyObject<'a, 'py> for NativeProperty {
    type Error = PyErr;
    fn extract(value: Borrowed<'a, 'py, PyAny>) -> PyResult<Self> {
        let property = if value.is_instance_of::<PyBool>() {
            ck::PropertyValue::Bool(value.extract()?)
        } else if value.is_instance_of::<PyInt>() {
            if let Ok(v) = value.extract::<i32>() {
                ck::PropertyValue::Int(v)
            } else {
                ck::PropertyValue::UInt(value.extract::<u32>().map_err(|_| {
                    PyOverflowError::new_err("atom property integer must fit i32 or u32")
                })?)
            }
        } else if value.is_instance_of::<PyFloat>() {
            ck::PropertyValue::Double(value.extract()?)
        } else if value.is_instance_of::<PyString>() {
            ck::PropertyValue::from(value.extract::<String>()?)
        } else if let Ok(list) = value.cast::<PyList>() {
            if list.iter().all(|v| v.is_instance_of::<PyString>()) {
                ck::PropertyValue::StringVector(
                    list.extract::<Vec<String>>()?
                        .into_iter()
                        .map(Into::into)
                        .collect(),
                )
            } else if list
                .iter()
                .all(|v| v.is_instance_of::<PyInt>() && !v.is_instance_of::<PyBool>())
            {
                ck::PropertyValue::IntVector(list.extract()?)
            } else {
                return Err(PyTypeError::new_err(
                    "atom property lists must contain only strings or only i32 integers",
                ));
            }
        } else {
            return Err(PyTypeError::new_err(
                "atom property value must be bool, int, float, str, list[int] or list[str]",
            ));
        };
        Ok(Self(property))
    }
}

impl NativeProperty {
    pub(crate) fn to_python(value: &ck::PropertyValue, py: Python<'_>) -> PyResult<Py<PyAny>> {
        Ok(match value {
            ck::PropertyValue::Bool(v) => v.into_pyobject(py)?.to_owned().into_any().unbind(),
            ck::PropertyValue::Int(v) => v.into_pyobject(py)?.into_any().unbind(),
            ck::PropertyValue::UInt(v) => v.into_pyobject(py)?.into_any().unbind(),
            ck::PropertyValue::Double(v) => v.into_pyobject(py)?.into_any().unbind(),
            ck::PropertyValue::String(v) => crate::canonical_sdf::decode_source_text(py, v)?
                .into_pyobject(py)?
                .into_any()
                .unbind(),
            ck::PropertyValue::IntVector(v) => v.into_pyobject(py)?.into_any().unbind(),
            ck::PropertyValue::StringVector(v) => v
                .iter()
                .map(|s| crate::canonical_sdf::decode_source_text(py, s))
                .collect::<PyResult<Vec<_>>>()?
                .into_pyobject(py)?
                .into_any()
                .unbind(),
        })
    }
}

#[cfg(feature = "stubgen")]
impl pyo3_stub_gen::PyStubType for NativeProperty {
    fn type_output() -> pyo3_stub_gen::TypeInfo {
        use pyo3_stub_gen::{PyStubType, TypeInfo};
        TypeInfo::builtin("bool")
            | TypeInfo::builtin("int")
            | TypeInfo::builtin("float")
            | TypeInfo::builtin("str")
            | Vec::<i32>::type_output()
            | Vec::<String>::type_output()
    }
}
