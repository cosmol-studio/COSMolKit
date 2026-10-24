//! Numeric transport only: shape errors are projected before scalar conversion.
use crate::coordinate_input_errors::CoordinateInputError;
use cosmolkit as ck;
use numpy::ndarray::Ix2;
use numpy::{AllowTypeChange, PyArrayLike, PyUntypedArray, PyUntypedArrayMethods};
use pyo3::BoundObject;
use pyo3::prelude::*;

pub(crate) struct Coordinates3D<'py>(Bound<'py, PyAny>);

impl<'a, 'py> FromPyObject<'a, 'py> for Coordinates3D<'py> {
    type Error = PyErr;
    const INPUT_TYPE: pyo3::inspect::PyStaticExpr = <Vec<[f64; 3]> as FromPyObject>::INPUT_TYPE;

    fn extract(value: Borrowed<'a, 'py, PyAny>) -> PyResult<Self> {
        Ok(Self(value.to_owned().into_bound()))
    }
}

impl Coordinates3D<'_> {
    pub(crate) fn into_rows(self, expected_rows: usize) -> PyResult<Vec<[f64; 3]>> {
        let check_rows = |actual| {
            if actual == expected_rows {
                Ok(())
            } else {
                Err(CoordinateInputError::from_rust(
                    ck::CoordinateInputError::RowCount {
                        dimension: "3D",
                        expected: expected_rows,
                        actual,
                    },
                ))
            }
        };
        let check_columns = |row, columns| {
            if columns == 3 {
                Ok(())
            } else {
                Err(CoordinateInputError::from_rust(
                    ck::CoordinateInputError::Shape {
                        dimension: "3D",
                        row,
                        columns,
                        expected: "3",
                    },
                ))
            }
        };
        if let Ok(array) = self.0.cast::<PyUntypedArray>() {
            let [rows, columns] = array.shape() else {
                return Err(pyo3::exceptions::PyTypeError::new_err(
                    "3D coordinates must be a 2D numeric array",
                ));
            };
            check_rows(*rows)?;
            check_columns(0, *columns)?;
            let converted = self
                .0
                .extract::<PyArrayLike<'_, f64, Ix2, AllowTypeChange>>()
                .map_err(|error| {
                    pyo3::exceptions::PyTypeError::new_err(format!(
                        "3D coordinates must be a 2D numeric array: {error}"
                    ))
                })?;
            // ndarray indexing supports C/F layout and noncontiguous/reversed strides.
            return Ok(converted
                .as_array()
                .outer_iter()
                .map(|row| [row[0], row[1], row[2]])
                .collect());
        }
        check_rows(self.0.len()?)?;
        let rows = self.0.extract::<Vec<Bound<'_, PyAny>>>()?;
        for (index, row) in rows.iter().enumerate() {
            check_columns(index, row.len()?)?;
        }
        rows.iter().map(|row| row.extract::<[f64; 3]>()).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numeric_transport_and_shape_errors() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let globals = pyo3::types::PyDict::new(py);
            globals.set_item("np", py.import("numpy")?)?;
            // All representations have the same values, including strided arrays.
            for expression in [
                "[[0, 1, 2], [3, 4, 5]]",
                "((0, 1, 2), (3, 4, 5))",
                "np.arange(6, dtype=np.float64).reshape(2, 3)",
                "np.asfortranarray(np.arange(6).reshape(2, 3))",
                "np.arange(12, dtype=np.float32).reshape(2, 6)[:, ::2] / 2",
                "np.arange(5, -1, -1).reshape(2, 3)[::-1, ::-1]",
            ] {
                let expression = std::ffi::CString::new(expression).unwrap();
                let input = py.eval(&expression, Some(&globals), None)?;
                assert_eq!(
                    input.extract::<Coordinates3D>()?.into_rows(2)?,
                    vec![[0., 1., 2.], [3., 4., 5.]]
                );
            }
            for (expression, variant) in [
                ("np.full((1, 2), 'bad')", "RowCount"),
                ("np.full((2, 2), 'bad')", "Shape"),
                ("[[0, 1, 2], [3, 4]]", "Shape"),
            ] {
                let expression = std::ffi::CString::new(expression).unwrap();
                let input = py.eval(&expression, Some(&globals), None)?;
                let error = input.extract::<Coordinates3D>()?.into_rows(2).unwrap_err();
                assert!(error.is_instance(py, <CoordinateInputError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
                assert_eq!(
                    error
                        .value(py)
                        .getattr("variant")?
                        .getattr("name")?
                        .extract::<String>()?,
                    variant
                );
            }
            let empty = py.eval(c"[]", None, None)?;
            assert!(empty.extract::<Coordinates3D>()?.into_rows(0)?.is_empty());
            let invalid = py.eval(c"np.full((2, 3), 'bad')", Some(&globals), None)?;
            let error = invalid
                .extract::<Coordinates3D>()?
                .into_rows(2)
                .unwrap_err();
            assert!(error.is_instance_of::<pyo3::exceptions::PyTypeError>(py));
            Ok(())
        })
        .unwrap();
    }
}
