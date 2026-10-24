//! Complete detached MOL2 read errors; no postprocessing or parser duplication.
use crate::canonical_errors::PropertyStringError;
use crate::error_framework::{Source, error_projection};
use crate::model_edit_errors::{AtomPropertyError, MoleculePropertyError};
use crate::model_errors::{CoordinateValidationError, TopologyValidationError};
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

type Mol2Source = Source<
    crate::error_framework::ErrorObject<TopologyValidationError>,
    Source<
        crate::error_framework::ErrorObject<CoordinateValidationError>,
        Source<
            crate::error_framework::ErrorObject<AtomPropertyError>,
            crate::error_framework::ErrorObject<MoleculePropertyError>,
        >,
    >,
>;

error_projection!(Mol2ReadError, Mol2ReadErrorKind, ck::Mol2ReadError, "io", {
    Parse(_), Unsupported { feature: &'static str }, Topology(_), Coordinates(_),
    PropertyString(_), AtomProperty(_), MoleculeProperty(_),
}, cause = |inner: &ck::Mol2ReadError| match inner {
    ck::Mol2ReadError::Topology(value) => Some(TopologyValidationError::from_rust(value.clone())),
    ck::Mol2ReadError::Coordinates(value) => Some(CoordinateValidationError::from_rust(value.clone())),
    ck::Mol2ReadError::AtomProperty(value) => Some(AtomPropertyError::from_rust(value.clone())),
    ck::Mol2ReadError::MoleculeProperty(value) => Some(MoleculePropertyError::from_rust(value.clone())),
    // Native PropertyString is transparent; its source-free child is a payload,
    // not an extra Error::source level introduced by the projection.
    ck::Mol2ReadError::Parse(_) | ck::Mol2ReadError::Unsupported { .. }
    | ck::Mol2ReadError::PropertyString(_) => None,
}, source = Option<Mol2Source>, methods = {
    #[getter]
    fn message(&self) -> Option<&str> {
        match &self.inner { ck::Mol2ReadError::Parse(value) => Some(value), _ => None }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn property_error(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<PropertyStringError>>> {
        match self.inner {
            ck::Mol2ReadError::PropertyString(value) => PropertyStringError::from_rust(value)
                .value(py).extract().map(Some).map_err(Into::into),
            _ => Ok(None),
        }
    }
});

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn actual_mol2_parse_failure_and_all_native_cause_edges_remain_typed() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native = ck::Molecule::from_mol2("not a MOL2 block\n").unwrap_err();
            let ck::MolecularIoError::Mol2Read(native) = native else {
                panic!("expected the actual detached MOL2 parse error")
            };
            let ck::Mol2ReadError::Parse(ref message) = native else {
                panic!("expected actual missing-section parse failure")
            };
            let message = message.clone();
            let display = native.to_string();
            assert!(native.source().is_none());
            let error = Mol2ReadError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert_eq!(
                error.value(py).getattr("message")?.extract::<String>()?,
                message
            );
            assert!(
                error
                    .value(py)
                    .getattr("variant")?
                    .extract::<Mol2ReadErrorKind>()?
                    == Mol2ReadErrorKind::Parse
            );
            assert!(error.value(py).getattr("source")?.is_none());
            assert!(error.value(py).setattr("message", "different").is_err());

            for (native, kind) in [
                (
                    ck::Mol2ReadError::Topology(ck::TopologyValidationError::AtomIdMismatch {
                        position: usize::MAX,
                        id: ck::AtomId::new(7),
                    }),
                    Mol2ReadErrorKind::Topology,
                ),
                (
                    ck::Mol2ReadError::Coordinates(
                        ck::CoordinateValidationError::Missing3DConformer { id: usize::MAX },
                    ),
                    Mol2ReadErrorKind::Coordinates,
                ),
                (
                    ck::Mol2ReadError::AtomProperty(ck::AtomPropertyError::EmptyKey),
                    Mol2ReadErrorKind::AtomProperty,
                ),
                (
                    ck::Mol2ReadError::MoleculeProperty(ck::MoleculePropertyError::EmptyKey),
                    Mol2ReadErrorKind::MoleculeProperty,
                ),
            ] {
                let source_display = native.source().unwrap().to_string();
                let display = native.to_string();
                let error = Mol2ReadError::from_rust(native);
                assert_eq!(error.value(py).str()?.to_str()?, display);
                assert!(
                    error
                        .value(py)
                        .getattr("variant")?
                        .extract::<Mol2ReadErrorKind>()?
                        == kind
                );
                let source = error.value(py).getattr("source")?;
                assert!(source.is(&error.value(py).getattr("__cause__")?));
                assert_eq!(source.str()?.to_str()?, source_display);
                assert!(source.getattr("__cause__")?.is_none());
                assert!(error.value(py).getattr("message")?.is_none());
                assert!(error.value(py).getattr("property_error")?.is_none());
                assert!(error.value(py).setattr("source", py.None()).is_err());
                match kind {
                    Mol2ReadErrorKind::Topology => {
                        assert!(source.is_instance(<TopologyValidationError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
                        assert_eq!(source.getattr("position")?.extract::<usize>()?, usize::MAX);
                        assert_eq!(source.getattr("id")?.extract::<usize>()?, 7);
                    }
                    Mol2ReadErrorKind::Coordinates => {
                        assert!(source.is_instance(<CoordinateValidationError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
                        assert_eq!(source.getattr("id")?.extract::<usize>()?, usize::MAX);
                    }
                    Mol2ReadErrorKind::AtomProperty => {
                        assert!(source.is_instance(<AtomPropertyError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?)
                    }
                    Mol2ReadErrorKind::MoleculeProperty => {
                        assert!(source.is_instance(<MoleculePropertyError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?)
                    }
                    _ => unreachable!(),
                }
            }

            let native =
                ck::Mol2ReadError::PropertyString(ck::PropertyStringError::UnsupportedKind {
                    kind: ck::PropertyValueKind::IntVector,
                });
            assert!(native.source().is_none());
            let display = native.to_string();
            let error = Mol2ReadError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert!(error.value(py).getattr("source")?.is_none());
            let payload = error.value(py).getattr("property_error")?;
            assert!(payload.is_instance(<PropertyStringError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert!(
                payload
                    .getattr("kind")?
                    .extract::<crate::property_values::PropertyValueKind>()?
                    == ck::PropertyValueKind::IntVector.into()
            );
            assert!(
                error
                    .value(py)
                    .setattr("property_error", py.None())
                    .is_err()
            );

            let error = Mol2ReadError::from_rust(ck::Mol2ReadError::Unsupported {
                feature: "scope\0é",
            });
            assert_eq!(
                error.value(py).getattr("feature")?.extract::<String>()?,
                "scope\0é"
            );
            assert!(error.value(py).getattr("source")?.is_none());
            assert!(error.value(py).setattr("feature", "other").is_err());
            let error = Mol2ReadError::from_rust(ck::Mol2ReadError::Parse("raw\0é".into()));
            assert_eq!(
                error.value(py).getattr("message")?.extract::<String>()?,
                "raw\0é"
            );
            Ok(())
        })
        .unwrap();
    }
}
