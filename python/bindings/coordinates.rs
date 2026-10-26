//! Coordinate values copy rows into native NumPy arrays without generating geometry.
use cosmolkit as ck;
use numpy::IntoPyArray;
use pyo3::prelude::*;

pub(crate) fn array<'py, const DIM: usize>(
    py: Python<'py>,
    rows: &[[f64; DIM]],
) -> Bound<'py, numpy::PyArray2<f64>> {
    numpy::ndarray::Array2::from_shape_fn((rows.len(), DIM), |(row, col)| rows[row][col])
        .into_pyarray(py)
}

crate::configuration::type_projection! {
    /// Owned 2D coordinate rows with a conformer identifier. Reading coordinates
    /// returns an independent NumPy array, not mutable molecular storage.
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    pub struct Conformer2D {
        pub(crate) inner: ck::Conformer2D,
    }
}
#[pymethods]
impl Conformer2D {
    #[new]
    fn new(id: usize, coordinates: Vec<[f64; 2]>) -> Self {
        Self {
            inner: ck::Conformer2D::new(id, coordinates),
        }
    }
    /// Identifier associated with these coordinate rows.
    #[getter]
    fn id(&self) -> usize {
        self.inner.id()
    }
    #[doc = r#"{"role": "return_type", "imports": ["numpy"], "type": "numpy.ndarray[tuple[int, int], numpy.dtype[numpy.float64]]"}"#]
    /// Return an independent float64 array of shape (num_atoms, 2), in atom-index order.
    #[getter]
    fn coordinates<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray2<f64>> {
        array(py, self.inner.coordinates())
    }
}

crate::configuration::type_projection! {
    /// Owned three-component coordinate rows and their explicit dimensionality flag.
    /// A false is_3d flag describes planar data stored as three-component rows.
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    pub struct Conformer3D {
        pub(crate) inner: ck::Conformer3D,
    }
}
#[pymethods]
impl Conformer3D {
    #[new]
    fn new(id: usize, coordinates: Vec<[f64; 3]>, is_3d: bool) -> Self {
        Self {
            inner: ck::Conformer3D::new(id, coordinates, is_3d),
        }
    }
    /// Identifier associated with these coordinate rows.
    #[getter]
    fn id(&self) -> usize {
        self.inner.id()
    }
    /// Whether these rows are marked as three-dimensional rather than planar.
    #[getter]
    fn is_3d(&self) -> bool {
        self.inner.is_3d()
    }
    #[doc = r#"{"role": "return_type", "imports": ["numpy"], "type": "numpy.ndarray[tuple[int, int], numpy.dtype[numpy.float64]]"}"#]
    /// Return an independent float64 array of shape (num_atoms, 3), in atom-index order.
    #[getter]
    fn coordinates<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray2<f64>> {
        array(py, self.inner.coordinates())
    }
}
