//! Exact native coordinate input errors and selector parsing.
use crate::error_framework::error_projection;
use crate::model_errors::CoordinateValidationError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

#[derive(IntoPyObject)]
enum InputValue {
    Name(String),
    Coordinate(f64),
}

#[derive(IntoPyObject)]
enum Expected {
    Count(usize),
    Shape(&'static str),
}

error_projection!(CoordinateInputError, CoordinateInputErrorKind, ck::CoordinateInputError, "coordinates", {
    UnknownZPolicy { value },
    RowCount { dimension, expected, actual: usize },
    Shape { dimension, row, columns: usize, expected },
    NonFinite { dimension, row, column: usize, value },
    NonZeroZ { row, z: f64 }, ThreeColumnsForbidden,
    ConformerNotFound { conformer_id: usize, count: usize },
    ConformerIdOverflow { max_id: usize }, InvalidCoordinates(_),
}, cause = |_: &ck::CoordinateInputError| None::<PyErr>, methods = {
    #[getter]
    #[doc = "Coordinate or matrix dimension involved in the failure. Applies to CoordinateInputError.RowCount, CoordinateInputError.Shape, CoordinateInputError.NonFinite. None when this payload is absent."] fn dimension(&self) -> Option<&'static str> { match self.inner {
        ck::CoordinateInputError::RowCount { dimension, .. }
        | ck::CoordinateInputError::Shape { dimension, .. }
        | ck::CoordinateInputError::NonFinite { dimension, .. } => Some(dimension), _ => None,
    } }
    #[getter]
    #[doc = "Row reported by the native failure. Applies to CoordinateInputError.Shape, CoordinateInputError.NonFinite, CoordinateInputError.NonZeroZ. None when this payload is absent."] fn row(&self) -> Option<usize> { match self.inner {
        ck::CoordinateInputError::Shape { row, .. }
        | ck::CoordinateInputError::NonFinite { row, .. }
        | ck::CoordinateInputError::NonZeroZ { row, .. } => Some(row), _ => None,
    } }
    #[getter]
    #[doc = "Value rejected or reported by the native operation. Applies to CoordinateInputError.UnknownZPolicy, CoordinateInputError.NonFinite. None when this payload is absent."] fn value(&self) -> Option<InputValue> { match &self.inner {
        ck::CoordinateInputError::UnknownZPolicy { value } => Some(InputValue::Name(value.clone())),
        ck::CoordinateInputError::NonFinite { value, .. } => Some(InputValue::Coordinate(*value)), _ => None,
    } }
    #[getter]
    #[doc = "Required value reported by the native check. Applies to CoordinateInputError.RowCount, CoordinateInputError.Shape. None when this payload is absent."] fn expected(&self) -> Option<Expected> { match self.inner {
        ck::CoordinateInputError::RowCount { expected, .. } => Some(Expected::Count(expected)),
        ck::CoordinateInputError::Shape { expected, .. } => Some(Expected::Shape(expected)), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured CoordinateValidationError payload of CoordinateInputError.InvalidCoordinates; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn coordinates(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<CoordinateValidationError>>> {
        match &self.inner {
            ck::CoordinateInputError::InvalidCoordinates(value) => CoordinateValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),
            _ => Ok(None),
        }
    }
});

#[cfg(test)]
mod tests {
    use super::*;
    use crate::selectors::CoordinateZPolicy;

    #[test]
    fn native_selector_and_coordinate_payloads_preserve_exact_error_shapes() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let cls = py.get_type::<CoordinateZPolicy>();
            for text in ["ignore", "IGNORE", "Require_Zero", "error"] {
                let actual = cls
                    .call_method1("from_name", (text,))?
                    .extract::<CoordinateZPolicy>()?;
                assert!(
                    ck::CoordinateZPolicy::from(actual)
                        == ck::CoordinateZPolicy::from_name(text).unwrap()
                );
            }
            let native = ck::CoordinateZPolicy::from_name(" not a policy ").unwrap_err();
            let display = native.to_string();
            let actual = cls
                .call_method1("from_name", (" not a policy ",))
                .unwrap_err();
            assert!(actual.is_instance(py, <CoordinateInputError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
            assert_eq!(actual.value(py).str()?.to_str()?, display);
            assert_eq!(
                actual.value(py).getattr("value")?.extract::<String>()?,
                " not a policy "
            );
            assert!(actual.cause(py).is_none());
            let nan = f64::from_bits(0x7ff8_0000_0000_0019);
            let error = CoordinateInputError::from_rust(ck::CoordinateInputError::NonFinite {
                dimension: "3D",
                row: 8,
                column: 2,
                value: nan,
            });
            assert_eq!(
                error
                    .value(py)
                    .getattr("value")?
                    .extract::<f64>()?
                    .to_bits(),
                nan.to_bits()
            );
            let error = CoordinateInputError::from_rust(ck::CoordinateInputError::Shape {
                dimension: "2D",
                row: 3,
                columns: 1,
                expected: "2 or 3",
            });
            assert_eq!(
                error.value(py).getattr("expected")?.extract::<String>()?,
                "2 or 3"
            );
            let native = ck::CoordinateInputError::InvalidCoordinates(
                ck::CoordinateValidationError::RowCount {
                    dimension: "3D",
                    conformer: 7,
                    rows: 2,
                    atom_count: 5,
                },
            );
            assert!(native.source().is_none());
            let error = CoordinateInputError::from_rust(native);
            assert!(error.cause(py).is_none());
            assert!(
                error
                    .value(py)
                    .getattr("coordinates")?
                    .is_instance(<CoordinateValidationError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?
            );
            Ok(())
        })
        .unwrap();
    }
}
