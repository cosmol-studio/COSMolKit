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

#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub struct Conformer2D {
    pub(crate) inner: ck::Conformer2D,
}
#[pymethods]
impl Conformer2D {
    #[new]
    fn new(id: usize, coordinates: Vec<[f64; 2]>) -> Self {
        Self {
            inner: ck::Conformer2D::new(id, coordinates),
        }
    }
    #[getter]
    fn id(&self) -> usize {
        self.inner.id()
    }
    #[doc = r#"{"role": "return_type", "imports": ["numpy"], "type": "numpy.ndarray[tuple[int, int], numpy.dtype[numpy.float64]]"}"#]
    fn coordinates<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray2<f64>> {
        array(py, self.inner.coordinates())
    }
}

#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub struct Conformer3D {
    pub(crate) inner: ck::Conformer3D,
}
#[pymethods]
impl Conformer3D {
    #[new]
    fn new(id: usize, coordinates: Vec<[f64; 3]>, is_3d: bool) -> Self {
        Self {
            inner: ck::Conformer3D::new(id, coordinates, is_3d),
        }
    }
    #[getter]
    fn id(&self) -> usize {
        self.inner.id()
    }
    #[getter]
    fn is_3d(&self) -> bool {
        self.inner.is_3d()
    }
    #[doc = r#"{"role": "return_type", "imports": ["numpy"], "type": "numpy.ndarray[tuple[int, int], numpy.dtype[numpy.float64]]"}"#]
    fn coordinates<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray2<f64>> {
        array(py, self.inner.coordinates())
    }
}
