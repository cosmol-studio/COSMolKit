//! Exact RDKFingerprint argument/provenance errors and explicit native causes.
use crate::error_framework::{Source, error_projection};
use crate::errors::FingerprintError;
use crate::layered_fingerprint_errors::LayeredFingerprintError;
use crate::model_errors::TopologyValidationError;
use crate::query_model_errors::QueryGraphError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

type TopologicalSource = Source<
    crate::error_framework::ErrorObject<TopologyValidationError>,
    Source<
        crate::error_framework::ErrorObject<QueryGraphError>,
        Source<
            crate::error_framework::ErrorObject<LayeredFingerprintError>,
            crate::error_framework::ErrorObject<FingerprintError>,
        >,
    >,
>;

error_projection!(TopologicalFingerprintError, TopologicalFingerprintErrorKind, ck::TopologicalFingerprintError, "fingerprint", {
    InvalidArguments { reason: &'static str }, OutputNotRequested { field: &'static str },
    Topology(_), Query(_), Paths(_), Value(_),
}, cause = |inner: &ck::TopologicalFingerprintError| match inner {
    ck::TopologicalFingerprintError::InvalidArguments { .. } | ck::TopologicalFingerprintError::OutputNotRequested { .. } => None,
    ck::TopologicalFingerprintError::Topology(value) => Some(TopologyValidationError::from_rust(value.clone())),
    ck::TopologicalFingerprintError::Query(value) => Some(QueryGraphError::from_rust(value.clone())),
    ck::TopologicalFingerprintError::Paths(value) => Some(LayeredFingerprintError::from_borrowed(value)),
    ck::TopologicalFingerprintError::Value(value) => Some(FingerprintError::from_rust(value.clone())),
}, source = Option<TopologicalSource>);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn actual_topological_failures_keep_provenance_fields_and_all_native_cause_levels() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let molecule = ck::Molecule::from_smiles("CCO").unwrap();
            let before = molecule.to_binary().unwrap();
            let native = molecule
                .fingerprint_topological_with_params(&ck::TopologicalFingerprintParams {
                    fp_size: 0,
                    ..Default::default()
                })
                .unwrap_err();
            let message = native.to_string();
            let error = TopologicalFingerprintError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, message);
            assert_eq!(
                error.value(py).getattr("reason")?.extract::<String>()?,
                "fpSize==0"
            );
            assert!(error.value(py).getattr("source")?.is_none());
            assert!(error.value(py).setattr("reason", "other").is_err());
            let result = molecule.fingerprint_topological_with_output().unwrap();
            for native in [
                result.atom_bits().unwrap_err(),
                result.bit_info().unwrap_err(),
            ] {
                let ck::TopologicalFingerprintError::OutputNotRequested { field } = &native else {
                    panic!("actual unrequested output expected")
                };
                let field = *field;
                let display = native.to_string();
                let error = TopologicalFingerprintError::from_rust(native);
                assert_eq!(
                    error.value(py).getattr("field")?.extract::<String>()?,
                    field
                );
                assert_eq!(error.value(py).str()?.to_str()?, display);
                assert!(error.value(py).getattr("source")?.is_none());
                assert!(error.value(py).setattr("field", "different").is_err());
            }
            for (native, class) in [
                (
                    ck::TopologicalFingerprintError::Topology(
                        ck::TopologyValidationError::AdjacencyMismatch,
                    ),
                    "TopologyValidationError",
                ),
                (
                    ck::TopologicalFingerprintError::Query(ck::QueryGraphError::AdjacencyMismatch),
                    "QueryGraphError",
                ),
                (
                    ck::TopologicalFingerprintError::Paths(
                        ck::LayeredFingerprintError::InvalidArguments {
                            reason: "native path",
                        },
                    ),
                    "LayeredFingerprintError",
                ),
                (
                    ck::TopologicalFingerprintError::Value(
                        ck::FingerprintError::SparseIndexOutOfRange {
                            index: u64::MAX,
                            size: 1,
                        },
                    ),
                    "FingerprintError",
                ),
            ] {
                let display = native.to_string();
                let source = native.source().unwrap().to_string();
                let error = TopologicalFingerprintError::from_rust(native);
                assert_eq!(error.value(py).str()?.to_str()?, display);
                let actual = error.cause(py).unwrap();
                assert_eq!(actual.value(py).get_type().name()?.to_str()?, class);
                assert_eq!(actual.value(py).str()?.to_str()?, source);
                assert!(error.value(py).getattr("source")?.is(actual.value(py)));
            }
            assert_eq!(molecule.to_binary().unwrap(), before);
            Ok(())
        })
        .unwrap();
    }
}
