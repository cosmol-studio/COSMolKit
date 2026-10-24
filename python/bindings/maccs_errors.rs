//! The native MACCS error keeps each explicit typed source level.
use crate::error_framework::{Source, error_projection};
use crate::errors::FingerprintError;
use crate::match_errors::{QueryMatchContextError, SubstructMatchError};
use crate::model_errors::TopologyValidationError;
use crate::path_errors::PathError;
use crate::ring_errors::RingFindingError;
use crate::smarts_errors::SmartsParseError;
use crate::string_errors::QueryCompileError;
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

type MaccsSource = Source<
    crate::error_framework::ErrorObject<TopologyValidationError>,
    Source<
        crate::error_framework::ErrorObject<RingFindingError>,
        Source<
            crate::error_framework::ErrorObject<ValenceError>,
            Source<
                crate::error_framework::ErrorObject<PathError>,
                Source<
                    crate::error_framework::ErrorObject<SmartsParseError>,
                    Source<
                        crate::error_framework::ErrorObject<QueryCompileError>,
                        Source<
                            crate::error_framework::ErrorObject<QueryMatchContextError>,
                            Source<
                                crate::error_framework::ErrorObject<SubstructMatchError>,
                                crate::error_framework::ErrorObject<FingerprintError>,
                            >,
                        >,
                    >,
                >,
            >,
        >,
    >,
>;

error_projection!(MaccsFingerprintError, MaccsFingerprintErrorKind, ck::MaccsFingerprintError, "fingerprint", {
    UnsupportedOption { option: &'static str, reason: &'static str },
    MissingPattern { bit: usize },
    Topology(_), Rings(_), Valence(_), Paths(_), Smarts(_), QueryCompile(_),
    QueryContext(_), Match(_), Value(_),
}, cause = |inner: &ck::MaccsFingerprintError| match inner {
    ck::MaccsFingerprintError::UnsupportedOption { .. }
    | ck::MaccsFingerprintError::MissingPattern { .. } => None,
    ck::MaccsFingerprintError::Topology(value) => Some(TopologyValidationError::from_rust(value.clone())),
    ck::MaccsFingerprintError::Rings(value) => Some(RingFindingError::from_rust(value.clone())),
    ck::MaccsFingerprintError::Valence(value) => Some(ValenceError::from_rust(value.clone())),
    ck::MaccsFingerprintError::Paths(value) => Some(PathError::from_rust(value.clone())),
    ck::MaccsFingerprintError::Smarts(value) => Some(SmartsParseError::from_rust(value.clone())),
    ck::MaccsFingerprintError::QueryCompile(value) => Some(QueryCompileError::from_rust(value.clone())),
    ck::MaccsFingerprintError::QueryContext(value) => Some(QueryMatchContextError::from_rust(value.clone())),
    ck::MaccsFingerprintError::Match(value) => Some(SubstructMatchError::from_rust(value.clone())),
    ck::MaccsFingerprintError::Value(value) => Some(FingerprintError::from_rust(value.clone())),
}, source = Option<MaccsSource>);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn public_maccs_failure_and_typed_sources_keep_native_levels_and_payloads() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native = ck::Molecule::from_smiles("CCO").unwrap();
            let before = native.to_binary().unwrap();
            let expected = native.fingerprint_maccs_with_params(&ck::MaccsFingerprintParams { n_bits: 200 }).unwrap_err();
            let message = expected.to_string();
            let ck::MaccsFingerprintError::UnsupportedOption { option, reason } = expected else { panic!("native unsupported width expected") };
            let molecule = Py::new(py, crate::native::Molecule { inner: native })?;
            let keywords = pyo3::types::PyDict::new(py);
            keywords.set_item("n_bits", 200)?;
            let error = molecule.bind(py).call_method("fingerprint_maccs", (), Some(&keywords)).unwrap_err();
            assert!(error.is_instance(py, <MaccsFingerprintError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
            assert_eq!(error.value(py).str()?.to_str()?, message);
            assert_eq!(error.value(py).getattr("option")?.extract::<String>()?, option);
            assert_eq!(error.value(py).getattr("reason")?.extract::<String>()?, reason);
            assert!(error.value(py).getattr("source")?.is_none());
            assert!(error.value(py).setattr("reason", "different").is_err());
            assert_eq!(molecule.borrow(py).inner.to_binary().unwrap(), before);

            // Each wrapped category is an explicit native source, even when
            // that child's own source is absent. It must not be skipped as if
            // MACCS used transparent Error::source delegation.
            let cases = [
                (ck::MaccsFingerprintError::Topology(ck::TopologyValidationError::AdjacencyMismatch), "TopologyValidationError"),
                (ck::MaccsFingerprintError::Rings(ck::RingFindingError::Value { message: "ring context" }), "RingFindingError"),
                (ck::MaccsFingerprintError::Valence(ck::ValenceError::ExplicitValenceCacheNotInitialized { atom: ck::AtomId::new(usize::MAX) }), "ValenceError"),
                (ck::MaccsFingerprintError::Paths(ck::PathError::InvalidLengthRange { lower: 3, upper: 2 }), "PathError"),
                (ck::MaccsFingerprintError::Smarts(ck::parse_smarts("[").unwrap_err()), "SmartsParseError"),
                (ck::MaccsFingerprintError::QueryCompile(ck::QueryCompileError::InvalidGraph("native\0detail".to_owned())), "QueryCompileError"),
                (ck::MaccsFingerprintError::QueryContext(ck::QueryMatchContextError::UninitializedRings), "QueryMatchContextError"),
                (ck::MaccsFingerprintError::Match(ck::SubstructMatchError::PropertyInteger { property: "count", source: ck::PropertyIntReadError::InvalidKind { kind: ck::PropertyValueKind::Bool } }), "SubstructMatchError"),
                (ck::MaccsFingerprintError::Value(ck::FingerprintError::SparseIndexOutOfRange { index: u64::MAX, size: 166 }), "FingerprintError"),
            ];
            for (native, class) in cases {
                let expected = native.to_string();
                assert!(native.source().is_some());
                let error = MaccsFingerprintError::from_rust(native);
                assert_eq!(error.value(py).str()?.to_str()?, expected);
                let cause = error.value(py).getattr("source")?;
                assert!(cause.is(&error.value(py).getattr("__cause__")?));
                assert_eq!(cause.getattr("__class__")?.getattr("__name__")?.extract::<String>()?, class);
                assert!(error.value(py).setattr("source", py.None()).is_err());
                if class == "ValenceError" { assert_eq!(cause.getattr("atom")?.extract::<usize>()?, usize::MAX); }
                if class == "FingerprintError" { assert_eq!(cause.getattr("index")?.extract::<u64>()?, u64::MAX); }
                if class == "QueryCompileError" { assert_eq!(cause.getattr("message")?.extract::<String>()?, "native\0detail"); }
                if class == "SubstructMatchError" {
                    assert_eq!(cause.getattr("property")?.extract::<String>()?, "count");
                    assert!(cause.getattr("source")?.is_instance(<crate::property_numeric_errors::PropertyIntReadError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
                }
            }
            let error = MaccsFingerprintError::from_rust(ck::MaccsFingerprintError::MissingPattern { bit: usize::MAX });
            assert_eq!(error.value(py).getattr("bit")?.extract::<usize>()?, usize::MAX);
            assert!(error.value(py).getattr("source")?.is_none());
            Ok(())
        }).unwrap();
    }
}
