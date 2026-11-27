//! Graph information content with a borrowed matrix or numeric rows.
use crate::error_framework::error_projection;
use crate::native::Molecule;
use crate::native_numerical_values::DenseMatrix;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

error_projection!(IpcError, IpcErrorKind, ck::IpcError, "descriptors", {
    InvalidDistanceMatrixShape { row, actual: usize, expected: usize },
    InvalidBondEndpoint { bond: usize, endpoint: &'static str, atom: usize, atom_count: usize },
}, cause = |_: &ck::IpcError| None::<PyErr>, methods = {
    #[getter]
    /// Index of a row with the wrong width; None for an overall dimension mismatch.
    fn row(&self) -> Option<usize> {
        match self.inner {
            ck::IpcError::InvalidDistanceMatrixShape { row, .. } => row,
            ck::IpcError::InvalidBondEndpoint { .. } => None,
        }
    }
});

#[derive(FromPyObject)]
pub enum IpcMatrixInput {
    Dense(Py<DenseMatrix>),
    Rows(Vec<Vec<f64>>),
}

impl IpcMatrixInput {
    fn with_matrix<R>(
        &self,
        py: Python<'_>,
        apply: impl FnOnce(ck::IpcDistanceMatrix<'_>) -> R,
    ) -> R {
        match self {
            Self::Dense(matrix) => apply((&matrix.borrow(py).inner).into()),
            Self::Rows(rows) => apply(rows.as_slice().into()),
        }
    }
}

#[pymethods]
impl Molecule {
    /// Total information content from exact characteristic-polynomial coefficients.
    #[pyo3(signature = (distance_matrix=None))]
    fn ipc(&self, py: Python<'_>, distance_matrix: Option<IpcMatrixInput>) -> PyResult<f64> {
        match distance_matrix {
            Some(matrix) => self.ipc_with_distance_matrix(py, matrix),
            None => py.detach(|| self.inner.ipc()).map_err(IpcError::from_rust),
        }
    }

    /// Use entries equal to one in the lower triangle, including the diagonal.
    fn ipc_with_distance_matrix(&self, py: Python<'_>, matrix: IpcMatrixInput) -> PyResult<f64> {
        matrix
            .with_matrix(py, |matrix| {
                py.detach(|| self.inner.ipc_with_distance_matrix(matrix))
            })
            .map_err(IpcError::from_rust)
    }

    /// Entropy in bits of the normalized coefficient magnitudes.
    #[pyo3(signature = (distance_matrix=None))]
    fn avg_ipc(&self, py: Python<'_>, distance_matrix: Option<IpcMatrixInput>) -> PyResult<f64> {
        match distance_matrix {
            Some(matrix) => self.avg_ipc_with_distance_matrix(py, matrix),
            None => py
                .detach(|| self.inner.avg_ipc())
                .map_err(IpcError::from_rust),
        }
    }

    /// Average information content with a square graph override.
    fn avg_ipc_with_distance_matrix(
        &self,
        py: Python<'_>,
        matrix: IpcMatrixInput,
    ) -> PyResult<f64> {
        matrix
            .with_matrix(py, |matrix| {
                py.detach(|| self.inner.avg_ipc_with_distance_matrix(matrix))
            })
            .map_err(IpcError::from_rust)
    }
}
