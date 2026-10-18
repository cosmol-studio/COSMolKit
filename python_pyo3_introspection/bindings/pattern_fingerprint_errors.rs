//! Pattern fingerprint failures preserve every explicit native source edge.
use crate::error_framework::{Source, error_projection};
use crate::errors::FingerprintError;
use crate::match_errors::{QueryMatchContextError, SubstructMatchError};
use crate::model_edit_errors::QueryAtomConversionError;
use crate::model_errors::TopologyValidationError;
use crate::query_model_errors::QueryGraphError;
use crate::ring_errors::RingFindingError;
use crate::smarts_errors::SmartsParseError;
use crate::string_errors::QueryCompileError;
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::error::Error;

type PatternSource = Source<
    Py<TopologyValidationError>,
    Source<
        Py<QueryGraphError>,
        Source<
            Py<QueryAtomConversionError>,
            Source<
                Py<RingFindingError>,
                Source<
                    Py<SmartsParseError>,
                    Source<
                        Py<QueryCompileError>,
                        Source<
                            Py<QueryMatchContextError>,
                            Source<Py<SubstructMatchError>, Py<FingerprintError>>,
                        >,
                    >,
                >,
            >,
        >,
    >,
>;

error_projection!(PatternFingerprintError, PatternFingerprintErrorKind, ck::PatternFingerprintError, "fingerprint", {
    EmptyFingerprint, InvalidArguments { reason: &'static str },
    BitLengthMismatch { left: usize, right: usize },
    Topology(_), Query(_), QueryCarrier(_), Rings(_), Smarts(_),
    QueryCompile(_), QueryContext(_), Match(_), Value(_),
}, cause = |inner: &ck::PatternFingerprintError| match inner {
    ck::PatternFingerprintError::Topology(value) => Some(TopologyValidationError::from_rust(value.clone())),
    ck::PatternFingerprintError::Query(value) => Some(QueryGraphError::from_rust(value.clone())),
    ck::PatternFingerprintError::QueryCarrier(value) => Some(QueryAtomConversionError::from_rust(value.clone())),
    ck::PatternFingerprintError::Rings(value) => Some(RingFindingError::from_rust(value.clone())),
    ck::PatternFingerprintError::Smarts(value) => Some(SmartsParseError::from_rust(value.clone())),
    ck::PatternFingerprintError::QueryCompile(value) => Some(QueryCompileError::from_rust(value.clone())),
    ck::PatternFingerprintError::QueryContext(value) => Some(QueryMatchContextError::from_rust(value.clone())),
    ck::PatternFingerprintError::Match(value) => Some(SubstructMatchError::from_rust(value.clone())),
    ck::PatternFingerprintError::Value(value) => Some(FingerprintError::from_rust(value.clone())),
    ck::PatternFingerprintError::EmptyFingerprint | ck::PatternFingerprintError::InvalidArguments { .. }
    | ck::PatternFingerprintError::BitLengthMismatch { .. } => None,
}, source = Option<PatternSource>);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn actual_pattern_failure_and_all_source_categories_keep_native_fields_and_levels() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let molecule = ck::Molecule::from_smiles("CCO").unwrap();
            let before = molecule.to_binary().unwrap();
            let native = molecule.fingerprint_pattern_with_params(&ck::PatternFingerprintParams {
                n_bits: 0,
                ..Default::default()
            }).unwrap_err();
            assert!(matches!(native, ck::PatternFingerprintError::EmptyFingerprint));
            let display = native.to_string();
            let error = PatternFingerprintError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert!(error.value(py).getattr("source")?.is_none());
            assert!(error.value(py).getattr("variant")?.extract::<PatternFingerprintErrorKind>()?
                == PatternFingerprintErrorKind::EmptyFingerprint);
            assert_eq!(molecule.to_binary().unwrap(), before);

            let cases = [
                (ck::PatternFingerprintError::Topology(ck::TopologyValidationError::AdjacencyMismatch), "TopologyValidationError"),
                (ck::PatternFingerprintError::Query(ck::QueryGraphError::AdjacencyMismatch), "QueryGraphError"),
                (ck::PatternFingerprintError::QueryCarrier(ck::QueryAtomConversionError::NonElementAtomicNumber { atom: ck::AtomId::new(usize::MAX), atomic_number: u8::MAX }), "QueryAtomConversionError"),
                (ck::PatternFingerprintError::Rings(ck::RingFindingError::Value { message: "ring context" }), "RingFindingError"),
                (ck::PatternFingerprintError::Smarts(ck::parse_smarts("[").unwrap_err()), "SmartsParseError"),
                (ck::PatternFingerprintError::QueryCompile(ck::QueryCompileError::InvalidGraph("native\0detail".to_owned())), "QueryCompileError"),
                (ck::PatternFingerprintError::QueryContext(ck::QueryMatchContextError::UninitializedRings), "QueryMatchContextError"),
                (ck::PatternFingerprintError::Match(ck::SubstructMatchError::PropertyInteger { property: "count", source: ck::PropertyIntReadError::InvalidKind { kind: ck::PropertyValueKind::Bool } }), "SubstructMatchError"),
                (ck::PatternFingerprintError::Value(ck::FingerprintError::SparseIndexOutOfRange { index: u64::MAX, size: 1 }), "FingerprintError"),
            ];
            for (native, class) in cases {
                let display = native.to_string();
                let source_display = native.source().unwrap().to_string();
                let error = PatternFingerprintError::from_rust(native);
                assert_eq!(error.value(py).str()?.to_str()?, display);
                let source = error.value(py).getattr("source")?;
                assert!(source.is(&error.value(py).getattr("__cause__")?));
                assert_eq!(source.str()?.to_str()?, source_display);
                assert_eq!(source.getattr("__class__")?.getattr("__name__")?.extract::<String>()?, class);
                assert!(error.value(py).setattr("source", py.None()).is_err());
                if class == "QueryAtomConversionError" {
                    assert_eq!(source.getattr("atom")?.extract::<usize>()?, usize::MAX);
                    assert_eq!(source.getattr("atomic_number")?.extract::<u8>()?, u8::MAX);
                }
                if class == "SubstructMatchError" {
                    assert!(source.getattr("source")?.is_instance_of::<crate::property_numeric_errors::PropertyIntReadError>());
                }
                if class == "FingerprintError" {
                    assert_eq!(source.getattr("index")?.extract::<u64>()?, u64::MAX);
                }
            }
            let error = PatternFingerprintError::from_rust(ck::PatternFingerprintError::InvalidArguments { reason: "native\0reason" });
            assert_eq!(error.value(py).getattr("reason")?.extract::<String>()?, "native\0reason");
            assert!(error.value(py).setattr("reason", "changed").is_err());
            assert!(error.value(py).getattr("source")?.is_none());
            let error = PatternFingerprintError::from_rust(ck::PatternFingerprintError::BitLengthMismatch { left: usize::MAX, right: 1 });
            assert_eq!(error.value(py).getattr("left")?.extract::<usize>()?, usize::MAX);
            assert_eq!(error.value(py).getattr("right")?.extract::<usize>()?, 1);
            assert!(error.value(py).getattr("source")?.is_none());
            Ok(())
        }).unwrap();
    }
}
