//! Complete detached SDF read failures, including recursive record positions.
use crate::error_framework::{Source, error_projection};
use crate::model_edit_errors::{AtomPropertyError, BondValueError, MoleculePropertyError};
use crate::model_errors::{CoordinateValidationError, StereoGroupError, TopologyValidationError};
use crate::mol_post_errors::MolPostError;
use crate::query_model_errors::QueryGraphError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

type ModelSource = Source<
    crate::error_framework::ErrorObject<TopologyValidationError>,
    Source<
        crate::error_framework::ErrorObject<CoordinateValidationError>,
        Source<
            crate::error_framework::ErrorObject<QueryGraphError>,
            Source<
                crate::error_framework::ErrorObject<AtomPropertyError>,
                Source<
                    crate::error_framework::ErrorObject<BondValueError>,
                    crate::error_framework::ErrorObject<MoleculePropertyError>,
                >,
            >,
        >,
    >,
>;
type ReadSource = Source<
    Source<
        crate::error_framework::ErrorObject<StereoGroupError>,
        crate::error_framework::ErrorObject<MolPostError>,
    >,
    Source<ModelSource, crate::error_framework::ErrorObject<SdfReadError>>,
>;

error_projection!(SdfReadError, SdfReadErrorKind, ck::SdfReadError, "io", {
    StereoGroup(_), MissingCapability(_), MolPost(_), QueryRecord, Empty, Counts,
    Field { kind: &'static str, line: usize, value: String }, Unsupported(_), Parse(_),
    Topology(_), Coordinates(_), QueryGraph(_), AtomProperty(_), BondValue(_), MoleculeProperty(_),
    PropertyListCount { target: &'static str, name: String, actual: usize, expected: usize },
    RecordIndexOutOfRange { index, record_count: usize },
    Record { index, byte_offset: u64, line_offset: usize, source },
}, cause = |inner: &ck::SdfReadError| match inner {
    ck::SdfReadError::StereoGroup(value) => Some(StereoGroupError::from_rust(*value)),
    ck::SdfReadError::MolPost(value) => Some(MolPostError::from_rust(value.clone())),
    ck::SdfReadError::Topology(value) => Some(TopologyValidationError::from_rust(value.clone())),
    ck::SdfReadError::Coordinates(value) => Some(CoordinateValidationError::from_rust(value.clone())),
    ck::SdfReadError::QueryGraph(value) => Some(QueryGraphError::from_rust(value.clone())),
    ck::SdfReadError::AtomProperty(value) => Some(AtomPropertyError::from_rust(value.clone())),
    ck::SdfReadError::BondValue(value) => Some(BondValueError::from_rust(value.clone())),
    ck::SdfReadError::MoleculeProperty(value) => Some(MoleculePropertyError::from_rust(value.clone())),
    ck::SdfReadError::Record { source, .. } => Some(SdfReadError::from_rust(source.as_ref().clone())),
    ck::SdfReadError::MissingCapability(_) | ck::SdfReadError::QueryRecord | ck::SdfReadError::Empty
    | ck::SdfReadError::Counts | ck::SdfReadError::Field { .. } | ck::SdfReadError::Unsupported(_)
    | ck::SdfReadError::Parse(_) | ck::SdfReadError::PropertyListCount { .. }
    | ck::SdfReadError::RecordIndexOutOfRange { .. } => None,
}, source = Option<ReadSource>, methods = {
    #[getter]
    #[doc = "Index associated with the failure. Present for `SdfReadError.RecordIndexOutOfRange`, `SdfReadError.Record`; None when unavailable."] fn index(&self) -> Option<usize> { match self.inner {
        ck::SdfReadError::RecordIndexOutOfRange { index, .. } | ck::SdfReadError::Record { index, .. } => Some(index), _ => None,
    } }
    #[getter]
    #[doc = "Capability required by this input. Present for `SdfReadError.MissingCapability`; None when unavailable."] fn capability(&self) -> Option<&'static str> { match self.inner {
        ck::SdfReadError::MissingCapability(value) => Some(value), _ => None,
    } }
    #[getter]
    #[doc = "Format associated with the failure. Present for `SdfReadError.Unsupported`; None when unavailable."] fn format(&self) -> Option<&'static str> { match self.inner {
        ck::SdfReadError::Unsupported(value) => Some(value), _ => None,
    } }
    #[getter]
    #[doc = "Native diagnostic message. Present for `SdfReadError.Parse`; None when unavailable."] fn message(&self) -> Option<&str> { match &self.inner {
        ck::SdfReadError::Parse(value) => Some(value), _ => None,
    } }
});

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mol_post_errors::MolProcessingError;
    use crate::mol_post_support_errors::BondDirectionStereoError;

    #[test]
    fn real_mol_failure_and_recursive_record_context_keep_exact_native_sources() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let ck::SdfError::Read(native) = ck::Molecule::from_mol("").unwrap_err() else {
                panic!("actual detached parser failure");
            };
            assert!(matches!(native, ck::SdfReadError::Empty));
            let display = native.to_string();
            let error = SdfReadError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert!(
                error
                    .value(py)
                    .getattr("variant")?
                    .extract::<SdfReadErrorKind>()?
                    == SdfReadErrorKind::Empty
            );
            assert!(error.value(py).getattr("source")?.is_none());

            let native = ck::SdfReadError::Record {
                index: usize::MAX,
                byte_offset: u64::MAX,
                line_offset: usize::MAX,
                source: Box::new(ck::SdfReadError::Record {
                    index: 7,
                    byte_offset: 9,
                    line_offset: 11,
                    source: Box::new(ck::SdfReadError::MolPost(ck::MolPostError::Processing(
                        ck::MolProcessingError::BondDirectionStereo(
                            ck::BondDirectionStereoError::InvalidState("state\0text".into()),
                        ),
                    ))),
                }),
            };
            assert!(native.source().unwrap().is::<ck::SdfReadError>());
            let display = native.to_string();
            let error = SdfReadError::from_rust(native);
            let value = error.value(py);
            assert_eq!(value.str()?.to_str()?, display);
            assert_eq!(value.getattr("index")?.extract::<usize>()?, usize::MAX);
            assert_eq!(value.getattr("byte_offset")?.extract::<u64>()?, u64::MAX);
            assert_eq!(
                value.getattr("line_offset")?.extract::<usize>()?,
                usize::MAX
            );
            let second = value.getattr("source")?;
            assert!(second.is_instance(<SdfReadError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert!(second.is(&value.getattr("__cause__")?));
            assert_eq!(second.getattr("index")?.extract::<usize>()?, 7);
            let molpost = second.getattr("source")?.getattr("source")?;
            assert!(molpost.is_instance(<MolPostError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            let processing = molpost.getattr("source")?;
            assert!(processing.is_instance(<MolProcessingError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            let direction = processing.getattr("source")?;
            assert!(direction.is_instance(<BondDirectionStereoError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert_eq!(
                direction.getattr("message")?.extract::<String>()?,
                "state\0text"
            );
            assert!(direction.getattr("source")?.is_none());
            for field in ["index", "byte_offset", "line_offset", "source"] {
                assert!(value.setattr(field, py.None()).is_err());
            }

            for native in [
                ck::SdfReadError::Field {
                    kind: "kind\0text",
                    line: usize::MAX,
                    value: "value\0text".into(),
                },
                ck::SdfReadError::PropertyListCount {
                    target: "target\0text",
                    name: "name\0text".into(),
                    actual: usize::MAX,
                    expected: 0,
                },
                ck::SdfReadError::RecordIndexOutOfRange {
                    index: usize::MAX,
                    record_count: 0,
                },
            ] {
                assert!(native.source().is_none());
                let display = native.to_string();
                let error = SdfReadError::from_rust(native);
                let value = error.value(py);
                assert_eq!(value.str()?.to_str()?, display);
                assert!(value.getattr("source")?.is_none());
                match value.getattr("variant")?.extract::<SdfReadErrorKind>()? {
                    SdfReadErrorKind::Field => {
                        assert_eq!(value.getattr("kind")?.extract::<String>()?, "kind\0text");
                        assert_eq!(value.getattr("line")?.extract::<usize>()?, usize::MAX);
                        assert_eq!(value.getattr("value")?.extract::<String>()?, "value\0text");
                    }
                    SdfReadErrorKind::PropertyListCount => {
                        assert_eq!(value.getattr("name")?.extract::<String>()?, "name\0text");
                        assert_eq!(value.getattr("actual")?.extract::<usize>()?, usize::MAX);
                    }
                    SdfReadErrorKind::RecordIndexOutOfRange => {
                        assert_eq!(value.getattr("index")?.extract::<usize>()?, usize::MAX)
                    }
                    _ => unreachable!(),
                }
            }
            Ok(())
        })
        .unwrap();
    }
}
