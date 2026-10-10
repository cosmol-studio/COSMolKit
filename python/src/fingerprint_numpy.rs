//! Python-only numeric projections of Rust-owned fingerprint values.

use crate::canonical_values::Fingerprint;
use ::cosmolkit as ck;
use numpy::IntoPyArray;
use pyo3::exceptions::{PyIndexError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PySlice;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

fn fill_row(fingerprint: &ck::Fingerprint, row: &mut [u8]) {
    // Only representation changes: preserve logical bit order, without an
    // intermediate on-bit list or a temporary dense array for each row.
    for (bit, value) in row.iter_mut().enumerate() {
        *value = u8::from(
            fingerprint
                .get_bit(bit as u32)
                .expect("row length equals the validated fingerprint width"),
        );
    }
}

pub(crate) fn to_numpy<'py>(
    py: Python<'py>,
    fingerprint: &ck::Fingerprint,
) -> Bound<'py, numpy::PyArray1<u8>> {
    let mut array = numpy::ndarray::Array1::zeros(fingerprint.n_bits() as usize);
    fill_row(fingerprint, array.as_slice_mut().unwrap());
    array.into_pyarray(py)
}

/// Ordered Rust-backed bit fingerprints, including failed input slots.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct FingerprintBatch {
    values: Vec<Option<Py<Fingerprint>>>,
    n_bits: usize,
}

impl FingerprintBatch {
    pub(crate) fn from_values(
        py: Python<'_>,
        values: Vec<Option<ck::Fingerprint>>,
    ) -> PyResult<Self> {
        let values = values
            .into_iter()
            .map(|value| {
                value
                    .map(|inner| Py::new(py, Fingerprint { inner }))
                    .transpose()
            })
            .collect::<PyResult<Vec<_>>>()?;
        Ok(Self::new(py, values))
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl FingerprintBatch {
    #[new]
    fn new(py: Python<'_>, fingerprints: Vec<Option<Py<Fingerprint>>>) -> Self {
        let n_bits = fingerprints
            .iter()
            .flatten()
            .next()
            .map_or(0, |value| value.borrow(py).inner.n_bits() as usize);
        Self {
            values: fingerprints,
            n_bits,
        }
    }

    /// Return an independent, C-contiguous uint8 matrix (rows, bits).
    /// None rows and unequal widths raise ValueError; rows are never discarded
    /// or silently replaced with zero fingerprints. An empty batch is (0, 0);
    /// empty slices retain the original batch's width.
    #[gen_stub(override_return_type(type_repr = "numpy.typing.NDArray[numpy.uint8]", imports = ("numpy", "numpy.typing")))]
    fn to_numpy<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, numpy::PyArray2<u8>>> {
        // Validate every row before allocating the output matrix.
        for (index, value) in self.values.iter().enumerate() {
            let value = value
                .as_ref()
                .ok_or_else(|| PyValueError::new_err(format!("fingerprint row {index} is None")))?;
            let actual = value.borrow(py).inner.n_bits() as usize;
            if actual != self.n_bits {
                return Err(PyValueError::new_err(format!(
                    "fingerprint row {index} has {actual} bits; expected {}",
                    self.n_bits
                )));
            }
        }
        let mut array = numpy::ndarray::Array2::zeros((self.values.len(), self.n_bits));
        for (value, mut row) in self.values.iter().zip(array.rows_mut()) {
            let value = value.as_ref().unwrap().borrow(py);
            fill_row(&value.inner, row.as_slice_mut().unwrap());
        }
        Ok(array.into_pyarray(py))
    }

    fn __len__(&self) -> usize {
        self.values.len()
    }

    #[gen_stub(skip)]
    fn __getitem__(&self, py: Python<'_>, key: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        if let Ok(slice) = key.cast::<PySlice>() {
            let indices = slice.indices(self.values.len() as isize)?;
            let values = (0..indices.slicelength)
                .map(|offset| {
                    let index = indices.start + offset as isize * indices.step;
                    let value = self.values[index as usize]
                        .as_ref()
                        .map(|value| value.clone_ref(py));
                    value
                })
                .collect();
            return Ok(Py::new(
                py,
                Self {
                    values,
                    n_bits: self.n_bits,
                },
            )?
            .into_any());
        }
        let mut index = key.extract::<isize>().map_err(|_| {
            PyTypeError::new_err("FingerprintBatch indices must be integers or slices")
        })?;
        if index < 0 {
            index += self.values.len() as isize;
        }
        let value = usize::try_from(index)
            .ok()
            .and_then(|index| self.values.get(index))
            .ok_or_else(|| PyIndexError::new_err("FingerprintBatch index out of range"))?;
        Ok(value
            .as_ref()
            .map_or_else(|| py.None(), |value| value.clone_ref(py).into_any()))
    }

    #[gen_stub(override_return_type(type_repr = "typing.Iterator[Fingerprint | None]", imports = ("typing")))]
    fn __iter__(slf: PyRef<'_, Self>) -> FingerprintBatchIterator {
        FingerprintBatchIterator {
            batch: slf.into(),
            index: 0,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "FingerprintBatch(rows={}, n_bits={})",
            self.values.len(),
            self.n_bits
        )
    }
}

#[pyclass]
struct FingerprintBatchIterator {
    batch: Py<FingerprintBatch>,
    index: usize,
}

#[pymethods]
impl FingerprintBatchIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> Option<Option<Py<Fingerprint>>> {
        let batch = self.batch.borrow(py);
        let value = batch.values.get(self.index)?;
        self.index += 1;
        Some(value.as_ref().map(|value| value.clone_ref(py)))
    }
}

#[cfg(feature = "stubgen")]
pyo3_stub_gen::inventory::submit! {
    pyo3_stub_gen::derive::gen_methods_from_python! {
        r#"
        class FingerprintBatch:
            @overload
            def __getitem__(self, key: int) -> Fingerprint | None: ...
            @overload
            def __getitem__(self, key: slice) -> FingerprintBatch: ...
        "#
    }
}
