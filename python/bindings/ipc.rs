//! Graph information content with a borrowed matrix or numeric rows.
use crate::chemical_values::enum_value;
use crate::error_framework::error_projection;
use crate::native::Molecule;
use crate::native_numerical_values::DenseMatrix;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

error_projection!(IpcError, IpcErrorKind, ck::IpcError, "descriptors", {
    InvalidDistanceMatrixShape { row, actual, expected: usize },
    InvalidBondEndpoint { bond: usize, endpoint: &'static str, atom: usize, atom_count: usize },
    EigenvaluesDidNotConverge { routine, info: usize },
    InvalidNumericalArgument { routine, argument: u8 },
    InsufficientNumericalStorage { routine, matrix: &'static str, required: usize, actual },
    NumericalDimensionOverflow,
    NumericalSortingStackExhausted,
}, cause = |_: &ck::IpcError| None::<PyErr>, methods = {
    #[getter]
    /// Index of a row with the wrong width; None for an overall dimension mismatch.
    fn row(&self) -> Option<usize> {
        match self.inner {
            ck::IpcError::InvalidDistanceMatrixShape { row, .. } => row,
            ck::IpcError::InvalidBondEndpoint { .. }
            | ck::IpcError::EigenvaluesDidNotConverge { .. }
            | ck::IpcError::InvalidNumericalArgument { .. }
            | ck::IpcError::InsufficientNumericalStorage { .. }
            | ck::IpcError::NumericalDimensionOverflow
            | ck::IpcError::NumericalSortingStackExhausted => None,
        }
    }
    #[getter]
    /// Actual matrix dimension or supplied numerical backing length; None when absent.
    fn actual(&self) -> Option<usize> {
        match self.inner {
            ck::IpcError::InvalidDistanceMatrixShape { actual, .. }
            | ck::IpcError::InsufficientNumericalStorage { actual, .. } => Some(actual),
            _ => None,
        }
    }
    #[getter]
    /// Source numerical routine; None for errors outside a named numerical routine.
    fn routine(&self) -> Option<&'static str> {
        match self.inner {
            ck::IpcError::EigenvaluesDidNotConverge { routine, .. }
            | ck::IpcError::InvalidNumericalArgument { routine, .. }
            | ck::IpcError::InsufficientNumericalStorage { routine, .. } => Some(routine),
            _ => None,
        }
    }
});

enum_value!(IpcMethod, { RdkitCompatible, ExactPolynomial });

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
    /// Total information content using the selected polynomial algorithm.
    /// ExactPolynomial uses exact integer coefficients and floating-point entropy.
    fn ipc_with_method(&self, py: Python<'_>, method: IpcMethod) -> PyResult<f64> {
        let method = method.into();
        py.detach(|| self.inner.ipc_with_method(method))
            .map_err(IpcError::from_rust)
    }

    /// Average information content using the selected polynomial algorithm.
    /// ExactPolynomial uses exact integer coefficients and floating-point entropy.
    fn avg_ipc_with_method(&self, py: Python<'_>, method: IpcMethod) -> PyResult<f64> {
        let method = method.into();
        py.detach(|| self.inner.avg_ipc_with_method(method))
            .map_err(IpcError::from_rust)
    }

    /// Total information content from floating-point characteristic-polynomial coefficients.
    #[pyo3(signature = (distance_matrix=None))]
    fn ipc(&self, py: Python<'_>, distance_matrix: Option<IpcMatrixInput>) -> PyResult<f64> {
        match distance_matrix {
            Some(matrix) => self.ipc_with_distance_matrix(py, matrix),
            None => py.detach(|| self.inner.ipc()).map_err(IpcError::from_rust),
        }
    }

    /// Use every matrix entry equal to one, including the diagonal.
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
