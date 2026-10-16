//! One implementation of integer/slice correspondence for object collections.
use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;
use pyo3::types::PySlice;
use pyo3::{BoundObject, PyClass, type_hint_union};

#[derive(FromPyObject)]
pub enum Index<'py> {
    Integer(isize),
    Slice(Bound<'py, PySlice>),
}

pub enum Indexed<T, B> {
    Item(T),
    Batch(B),
}

// PyO3 generic derives do not retain this branch relationship. Carry the
// actual branch types through a shared conversion instead of annotating a stub.
impl<'py, T: IntoPyObject<'py>, B: IntoPyObject<'py>> IntoPyObject<'py> for Indexed<T, B> {
    type Target = pyo3::types::PyAny;
    type Output = Bound<'py, Self::Target>;
    type Error = PyErr;
    const OUTPUT_TYPE: pyo3::inspect::PyStaticExpr =
        pyo3::type_hint_union!(T::OUTPUT_TYPE, B::OUTPUT_TYPE);

    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        match self {
            Self::Item(value) => value
                .into_pyobject(py)
                .map(BoundObject::into_bound)
                .map(Bound::into_any)
                .map_err(Into::into),
            Self::Batch(value) => value
                .into_pyobject(py)
                .map(BoundObject::into_bound)
                .map(Bound::into_any)
                .map_err(Into::into),
        }
    }
}

pub fn get<T: PyClass, B: PyClass>(
    py: Python<'_>,
    index: Index<'_>,
    values: &[Option<Py<T>>],
    wrap: impl FnOnce(Vec<Option<Py<T>>>) -> B,
) -> PyResult<Indexed<Option<Py<T>>, Py<B>>>
where
    B: Into<PyClassInitializer<B>>,
{
    match index {
        Index::Integer(index) => {
            let len = values.len() as isize;
            let normalized = if index < 0 { index + len } else { index };
            let item = usize::try_from(normalized)
                .ok()
                .and_then(|index| values.get(index))
                .ok_or_else(|| PyIndexError::new_err(index))?;
            Ok(Indexed::Item(
                item.as_ref().map(|value| value.clone_ref(py)),
            ))
        }
        Index::Slice(slice) => {
            let indices = slice.indices(values.len() as isize)?;
            let mut at = indices.start;
            let selected = (0..indices.slicelength)
                .map(|_| {
                    let item = values[at as usize]
                        .as_ref()
                        .map(|value| value.clone_ref(py));
                    at += indices.step;
                    item
                })
                .collect();
            Ok(Indexed::Batch(Py::new(py, wrap(selected))?))
        }
    }
}
