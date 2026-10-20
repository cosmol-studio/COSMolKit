//! NumPy is the only fingerprint-specific collection adaptation.
use crate::native::Fingerprint;
use cosmolkit as ck;
use numpy::PyArrayMethods;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::error::Error;

crate::error_framework::error_projection!(FingerprintBatchError, FingerprintBatchErrorKind, ck::FingerprintBatchError, "fingerprint", {
    MissingRow { index },
    WidthMismatch { index, actual: u32, expected: u32 },
}, methods = {
    #[getter]
    fn index(&self) -> usize {
        match self.inner {
            ck::FingerprintBatchError::MissingRow { index }
            | ck::FingerprintBatchError::WidthMismatch { index, .. } => index,
        }
    }
});

crate::sequence::project_sequence!(FingerprintBatch, ck::FingerprintBatch, Fingerprint, {
    #[doc = r#"{"role":"return_type","imports":["numpy"],"type":"numpy.ndarray[tuple[int, int], numpy.dtype[numpy.uint8]]"}"#]
    fn to_numpy<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, numpy::PyArray2<u8>>> {
        self.inner
            .validate_matrix()
            .map_err(FingerprintBatchError::from_rust)?;
        let data = numpy::PyArray2::<u8>::zeros(
            py,
            [self.inner.len(), self.inner.n_bits() as usize],
            false,
        );
        {
            let mut storage = data.readwrite();
            let mut array = storage.as_array_mut();
            for (fingerprint, mut row) in self.inner.iter().zip(array.rows_mut()) {
                let fingerprint = fingerprint.as_ref().expect("validated complete matrix");
                for (bit, value) in row.iter_mut().enumerate() {
                    *value = u8::from(
                        fingerprint
                            .get_bit(bit as u32)
                            .expect("validated matrix width"),
                    );
                }
            }
        }
        Ok(data)
    }
});
