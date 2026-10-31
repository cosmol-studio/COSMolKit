//! One implementation of integer/slice correspondence for object collections.
use pyo3::BoundObject;
use pyo3::exceptions::{PyIndexError, PyOverflowError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::type_hint_union;
use pyo3::types::{PyBool, PyBytes, PySlice, PyString};

pub struct Integer(pub isize);

impl<'a, 'py> FromPyObject<'a, 'py> for Integer {
    type Error = PyErr;
    const INPUT_TYPE: pyo3::inspect::PyStaticExpr = <isize as FromPyObject>::INPUT_TYPE;
    fn extract(value: Borrowed<'a, 'py, PyAny>) -> PyResult<Self> {
        if value.is_instance_of::<PyBool>() {
            return Err(PyTypeError::new_err("boolean is not an integer index"));
        }
        value.extract::<isize>().map(Self).map_err(|error| {
            if error.is_instance_of::<PyOverflowError>(value.py()) {
                PyIndexError::new_err("index exceeds addressable range")
            } else {
                error
            }
        })
    }
}

pub enum Index<'py> {
    Integer(Integer),
    Slice(Bound<'py, PySlice>),
    Indices(Vec<Integer>),
    Mask(Vec<bool>),
}

impl<'a, 'py> FromPyObject<'a, 'py> for Index<'py> {
    type Error = PyErr;
    const INPUT_TYPE: pyo3::inspect::PyStaticExpr = pyo3::type_hint_union!(
        Integer::INPUT_TYPE,
        <Bound<'py, PySlice> as FromPyObject>::INPUT_TYPE,
        <Vec<Integer> as FromPyObject>::INPUT_TYPE,
        <Vec<bool> as FromPyObject>::INPUT_TYPE
    );
    fn extract(value: Borrowed<'a, 'py, PyAny>) -> PyResult<Self> {
        if let Ok(slice) = value.cast::<PySlice>() {
            return Ok(Self::Slice(slice.to_owned().into_bound()));
        }
        if value.is_instance_of::<PyString>() || value.is_instance_of::<PyBytes>() {
            return Err(PyTypeError::new_err(
                "expected integer, slice, integer sequence or boolean mask",
            ));
        }
        // Match PyO3's Vec extraction protocol, not Sequence ABC membership:
        // NumPy arrays and custom __getitem__ sequences need not register there.
        // SAFETY: value holds a live Python object while attached to Python.
        if unsafe { pyo3::ffi::PySequence_Check(value.as_ptr()) } != 0 {
            let values = value.extract::<Vec<Bound<'py, PyAny>>>()?;
            if !values.is_empty() && values.iter().all(|v| v.is_instance_of::<PyBool>()) {
                return Ok(Self::Mask(
                    values
                        .iter()
                        .map(|v| v.extract())
                        .collect::<PyResult<_>>()?,
                ));
            }
            return Ok(Self::Indices(
                values
                    .iter()
                    .map(|v| v.extract())
                    .collect::<PyResult<_>>()?,
            ));
        }
        Integer::extract(value).map(Self::Integer)
    }
}

pub fn normalize(index: Integer, len: usize) -> PyResult<usize> {
    let len = isize::try_from(len).map_err(|_| PyOverflowError::new_err("collection too large"))?;
    let at = if index.0 < 0 { index.0 + len } else { index.0 };
    usize::try_from(at)
        .ok()
        .filter(|&i| i < len as usize)
        .ok_or_else(|| PyIndexError::new_err(index.0))
}

pub enum Selection {
    One(usize),
    Many(Vec<usize>),
}

impl Index<'_> {
    pub fn resolve(self, len: usize) -> PyResult<Selection> {
        match self {
            Self::Integer(index) => normalize(index, len).map(Selection::One),
            Self::Indices(indices) => indices
                .into_iter()
                .map(|i| normalize(i, len))
                .collect::<PyResult<Vec<_>>>()
                .map(Selection::Many),
            Self::Mask(mask) => {
                if mask.len() != len {
                    return Err(PyValueError::new_err(
                        "boolean mask length must equal collection length",
                    ));
                }
                Ok(Selection::Many(
                    mask.into_iter()
                        .enumerate()
                        .filter_map(|(i, selected)| selected.then_some(i))
                        .collect(),
                ))
            }
            Self::Slice(slice) => {
                let len = isize::try_from(len)
                    .map_err(|_| PyOverflowError::new_err("collection too large"))?;
                let indices = slice.indices(len)?;
                let mut at = indices.start;
                let selected = (0..indices.slicelength)
                    .map(|position| {
                        let current = at as usize;
                        if position + 1 < indices.slicelength {
                            at += indices.step;
                        }
                        current
                    })
                    .collect();
                Ok(Selection::Many(selected))
            }
        }
    }
}

/// Native storage stays in the owner. Only selected items are accessed.
pub trait Collection {
    type Item;
    type Selected;
    fn length(&self) -> usize;
    fn item(&self, py: Python<'_>, index: usize) -> PyResult<Self::Item>;
    fn select(&self, py: Python<'_>, indices: &[usize]) -> PyResult<Self::Selected>;
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

pub fn get<C: Collection>(
    collection: &C,
    py: Python<'_>,
    index: Index<'_>,
) -> PyResult<Indexed<C::Item, C::Selected>> {
    match index.resolve(collection.length())? {
        Selection::One(index) => collection.item(py, index).map(Indexed::Item),
        Selection::Many(indices) => collection.select(py, &indices).map(Indexed::Batch),
    }
}
