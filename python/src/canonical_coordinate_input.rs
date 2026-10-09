//! Typed numeric ingress and immutable parameters for canonical coordinate APIs.
//! Chemical interpretation and detached edits are forwarded to cosmolkit.
use ::cosmolkit as ck;
use numpy::ndarray::Ix2;
use numpy::{AllowTypeChange, PyArrayLike, PyUntypedArray, PyUntypedArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pymethods};
pyo3::create_exception!(cosmolkit, CoordinateInputError, PyValueError);
pyo3::create_exception!(cosmolkit, Coordinate3DReadError, PyValueError);

#[cosmolkit_macros::python_enum(existing_methods)]
#[cfg_attr(feature = "stubgen", gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", frozen, eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum CoordinateZPolicy {
    Ignore = 0,
    RequireZero = 1,
    Error = 2,
}
impl CoordinateZPolicy {
    fn core(self) -> ck::CoordinateZPolicy {
        match self {
            Self::Ignore => ck::CoordinateZPolicy::Ignore,
            Self::RequireZero => ck::CoordinateZPolicy::RequireZero,
            Self::Error => ck::CoordinateZPolicy::Error,
        }
    }
    fn from_core(value: ck::CoordinateZPolicy) -> Self {
        match value {
            ck::CoordinateZPolicy::Ignore => Self::Ignore,
            ck::CoordinateZPolicy::RequireZero => Self::RequireZero,
            ck::CoordinateZPolicy::Error => Self::Error,
        }
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl CoordinateZPolicy {
    #[classattr]
    fn _enum_string_values() -> Vec<(&'static str, CoordinateZPolicy)> {
        Self::enum_string_values()
    }
    #[staticmethod]
    fn from_name(py: Python<'_>, value: &str) -> PyResult<CoordinateZPolicy> {
        ck::CoordinateZPolicy::from_name(value)
            .map(Self::from_core)
            .map_err(|e| error_pyerr(py, &e))
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct Coordinate2DInputParams {
    pub(crate) inner: ck::Coordinate2DInputParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl Coordinate2DInputParams {
    #[new]
    #[pyo3(signature=(z_policy=CoordinateZPolicy::Ignore))]
    fn new(z_policy: CoordinateZPolicy) -> Self {
        Self {
            inner: ck::Coordinate2DInputParams {
                z_policy: z_policy.core(),
            },
        }
    }
    #[getter]
    fn z_policy(&self) -> CoordinateZPolicy {
        CoordinateZPolicy::from_core(self.inner.z_policy)
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct Coordinate3DInputParams {
    pub(crate) inner: ck::Coordinate3DInputParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl Coordinate3DInputParams {
    #[new]
    #[pyo3(signature=(is_3d=true))]
    fn new(is_3d: bool) -> Self {
        Self {
            inner: ck::Coordinate3DInputParams { is_3d },
        }
    }
    #[getter]
    fn is_3d(&self) -> bool {
        self.inner.is_3d
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct Replace3DCoordinatesParams {
    pub(crate) inner: ck::Replace3DCoordinatesParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl Replace3DCoordinatesParams {
    #[new]
    #[pyo3(signature=(conformer_id=0))]
    fn new(conformer_id: usize) -> Self {
        Self {
            inner: ck::Replace3DCoordinatesParams { conformer_id },
        }
    }
    #[getter]
    fn conformer_id(&self) -> usize {
        self.inner.conformer_id
    }
}

pub(crate) fn error_pyerr(py: Python<'_>, source: &ck::CoordinateInputError) -> PyErr {
    use ck::CoordinateInputError as E;
    let kind = match source {
        E::UnknownZPolicy { .. } => "UnknownZPolicy",
        E::RowCount { .. } => "RowCount",
        E::Shape { .. } => "Shape",
        E::NonFinite { .. } => "NonFinite",
        E::NonZeroZ { .. } => "NonZeroZ",
        E::ThreeColumnsForbidden => "ThreeColumnsForbidden",
        E::ConformerNotFound { .. } => "ConformerNotFound",
        E::ConformerIdOverflow { .. } => "ConformerIdOverflow",
        E::InvalidCoordinates(..) => "InvalidCoordinates",
    };
    let error = crate::canonical_values::annotate(
        py,
        CoordinateInputError::new_err(source.to_string()),
        "coordinate_input",
        kind,
        source,
    );
    let value = error.value(py);
    let annotate = (|| -> PyResult<()> {
        match source {
            E::UnknownZPolicy { value: input } => value.setattr("value", input)?,
            E::RowCount {
                dimension,
                expected,
                actual,
            } => {
                value.setattr("dimension", *dimension)?;
                value.setattr("expected", *expected)?;
                value.setattr("actual", *actual)?;
            }
            E::Shape {
                dimension,
                row,
                columns,
                expected,
            } => {
                value.setattr("dimension", *dimension)?;
                value.setattr("row", *row)?;
                value.setattr("columns", *columns)?;
                value.setattr("expected_columns", *expected)?;
            }
            E::NonFinite {
                dimension,
                row,
                column,
                value: coordinate,
            } => {
                value.setattr("dimension", *dimension)?;
                value.setattr("row", *row)?;
                value.setattr("column", *column)?;
                value.setattr("value", *coordinate)?;
            }
            E::NonZeroZ { row, z } => {
                value.setattr("row", *row)?;
                value.setattr("z", *z)?;
            }
            E::ConformerNotFound {
                conformer_id,
                count,
            } => {
                value.setattr("conformer_id", *conformer_id)?;
                value.setattr("count", *count)?;
            }
            E::ConformerIdOverflow { max_id } => value.setattr("max_id", *max_id)?,
            E::InvalidCoordinates(cause) => match cause {
                ck::CoordinateValidationError::AtomIndexOverflow { atom } => {
                    value.setattr("cause_kind", "AtomIndexOverflow")?;
                    value.setattr("atom", *atom)?;
                }
                ck::CoordinateValidationError::MissingSourceConformerOrder => {
                    value.setattr("cause_kind", "MissingSourceConformerOrder")?;
                }
                ck::CoordinateValidationError::SourceConformerOrder {
                    two_d,
                    three_d,
                    expected_two_d,
                    expected_three_d,
                } => {
                    value.setattr("cause_kind", "SourceConformerOrder")?;
                    value.setattr("two_d", *two_d)?;
                    value.setattr("three_d", *three_d)?;
                    value.setattr("expected_two_d", *expected_two_d)?;
                    value.setattr("expected_three_d", *expected_three_d)?;
                }
                ck::CoordinateValidationError::Missing3DConformer { id } => {
                    value.setattr("conformer_id", *id)?;
                }
                ck::CoordinateValidationError::ConformerIdOverflow { max_id } => {
                    value.setattr("max_id", *max_id)?;
                }
                ck::CoordinateValidationError::RowCount {
                    dimension,
                    conformer,
                    rows,
                    atom_count,
                } => {
                    value.setattr("dimension", *dimension)?;
                    value.setattr("conformer", *conformer)?;
                    value.setattr("actual", *rows)?;
                    value.setattr("expected", *atom_count)?;
                }
                ck::CoordinateValidationError::DuplicateConformerId { dimension, id } => {
                    value.setattr("dimension", *dimension)?;
                    value.setattr("conformer", *id)?;
                }
                ck::CoordinateValidationError::NonFiniteCoordinate {
                    dimension,
                    conformer,
                    atom,
                    axis,
                } => {
                    value.setattr("dimension", *dimension)?;
                    value.setattr("conformer", *conformer)?;
                    value.setattr("row", *atom)?;
                    value.setattr("axis", *axis)?;
                }
            },
            E::ThreeColumnsForbidden => {}
        }
        Ok(())
    })();
    if let Err(e) = annotate {
        return e;
    }
    error
}

/// Original ndarray row preflight occurs before dtype conversion/copying.
/// Numeric protocol conversion only; finite/z interpretation stays in core.
pub(crate) fn matrix(
    py: Python<'_>,
    coordinates: &Bound<'_, PyAny>,
    expected_rows: usize,
    two_d: bool,
) -> PyResult<Vec<Vec<f64>>> {
    let dimension = if two_d { "2D" } else { "3D" };
    if let Ok(array) = coordinates.cast::<PyUntypedArray>()
        && let [rows, _] = array.shape()
        && *rows != expected_rows
    {
        return Err(error_pyerr(
            py,
            &ck::CoordinateInputError::RowCount {
                dimension,
                expected: expected_rows,
                actual: *rows,
            },
        ));
    }
    let converted = coordinates
        .extract::<PyArrayLike<'_, f64, Ix2, AllowTypeChange>>()
        .map_err(|e| {
            pyo3::exceptions::PyTypeError::new_err(format!(
                "{dimension} coordinates must be a 2D numeric array: {e}"
            ))
        })?;
    let array = converted.as_array();
    let shape = array.shape();
    if shape[0] != expected_rows {
        return Err(error_pyerr(
            py,
            &ck::CoordinateInputError::RowCount {
                dimension,
                expected: expected_rows,
                actual: shape[0],
            },
        ));
    }
    if !(shape[1] == 3 || two_d && shape[1] == 2) {
        return Err(error_pyerr(
            py,
            &ck::CoordinateInputError::Shape {
                dimension,
                row: 0,
                columns: shape[1],
                expected: if two_d { "2 or 3" } else { "3" },
            },
        ));
    }
    Ok(array.outer_iter().map(|row| row.to_vec()).collect())
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add(
        "Coordinate3DReadError",
        module.py().get_type::<Coordinate3DReadError>(),
    )?;
    module.add_class::<CoordinateZPolicy>()?;
    module.add_class::<Coordinate2DInputParams>()?;
    module.add_class::<Coordinate3DInputParams>()?;
    module.add_class::<Replace3DCoordinatesParams>()?;
    module.add(
        "CoordinateInputError",
        module.py().get_type::<CoordinateInputError>(),
    )?;
    Ok(())
}

pub(crate) fn read_pyerr(py: Python<'_>, source: &ck::Coordinate3DReadError) -> PyErr {
    let error = crate::canonical_values::annotate(
        py,
        Coordinate3DReadError::new_err(source.to_string()),
        "coordinate_read",
        "ConformerNotFound",
        source,
    );
    let ck::Coordinate3DReadError::ConformerNotFound {
        conformer_id,
        count,
    } = source;
    if let Err(error) = error.value(py).setattr("conformer_id", *conformer_id) {
        return error;
    }
    if let Err(error) = error.value(py).setattr("count", *count) {
        return error;
    }
    error
}
