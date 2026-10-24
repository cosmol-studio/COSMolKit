//! Layered failures retain the native error and each explicit cause level.
use crate::error_framework::{Source, error_projection};
use crate::errors::FingerprintError;
use crate::model_errors::TopologyValidationError;
use crate::path_errors::PathError;
use crate::query_model_errors::QueryGraphError;
use crate::ring_errors::RingFindingError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

type LayeredSource = Source<
    crate::error_framework::ErrorObject<TopologyValidationError>,
    Source<
        crate::error_framework::ErrorObject<QueryGraphError>,
        Source<
            crate::error_framework::ErrorObject<RingFindingError>,
            Source<
                crate::error_framework::ErrorObject<PathError>,
                crate::error_framework::ErrorObject<FingerprintError>,
            >,
        >,
    >,
>;

error_projection!(LayeredFingerprintError, LayeredFingerprintErrorKind, ck::LayeredFingerprintError, "fingerprint", {
    InvalidArguments { reason: &'static str }, Topology(_), Query(_), Rings(_), Paths(_), Value(_),
}, cause = |inner: &ck::LayeredFingerprintError| match inner {
    ck::LayeredFingerprintError::InvalidArguments { .. } => None,
    ck::LayeredFingerprintError::Topology(value) => Some(TopologyValidationError::from_rust(value.clone())),
    ck::LayeredFingerprintError::Query(value) => Some(QueryGraphError::from_rust(value.clone())),
    ck::LayeredFingerprintError::Rings(value) => Some(RingFindingError::from_rust(value.clone())),
    ck::LayeredFingerprintError::Paths(value) => Some(PathError::from_rust(value.clone())),
    ck::LayeredFingerprintError::Value(value) => Some(FingerprintError::from_rust(value.clone())),
}, source = Option<LayeredSource>);

impl LayeredFingerprintError {
    pub(crate) fn from_borrowed(inner: &ck::LayeredFingerprintError) -> PyErr {
        // The native outer enum is not Clone, but all its actual public
        // payloads are. Exhaustive value copying preserves this borrowed
        // cause without changing native types or inferring any message.
        let owned = match inner {
            ck::LayeredFingerprintError::InvalidArguments { reason } => {
                ck::LayeredFingerprintError::InvalidArguments { reason }
            }
            ck::LayeredFingerprintError::Topology(value) => {
                ck::LayeredFingerprintError::Topology(value.clone())
            }
            ck::LayeredFingerprintError::Query(value) => {
                ck::LayeredFingerprintError::Query(value.clone())
            }
            ck::LayeredFingerprintError::Rings(value) => {
                ck::LayeredFingerprintError::Rings(value.clone())
            }
            ck::LayeredFingerprintError::Paths(value) => {
                ck::LayeredFingerprintError::Paths(value.clone())
            }
            ck::LayeredFingerprintError::Value(value) => {
                ck::LayeredFingerprintError::Value(value.clone())
            }
        };
        Self::from_rust(owned)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn actual_layered_failure_and_all_sources_keep_native_levels_and_readonly_payloads() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let molecule = ck::Molecule::from_smiles("CCO").unwrap();
            let before = molecule.to_binary().unwrap();
            let native = molecule
                .fingerprint_layered_with_params(&ck::LayeredFingerprintParams {
                    fp_size: 0,
                    ..Default::default()
                })
                .unwrap_err();
            let ck::LayeredFingerprintError::InvalidArguments { reason } = &native else {
                panic!("native zero-width precondition expected")
            };
            let expected_reason = *reason;
            let display = native.to_string();
            let error = LayeredFingerprintError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert_eq!(
                error.value(py).getattr("reason")?.extract::<String>()?,
                expected_reason
            );
            assert!(error.value(py).getattr("source")?.is_none());
            assert!(
                error
                    .value(py)
                    .getattr("variant")?
                    .extract::<LayeredFingerprintErrorKind>()?
                    == LayeredFingerprintErrorKind::InvalidArguments
            );
            assert!(error.value(py).setattr("reason", "different").is_err());
            assert_eq!(molecule.to_binary().unwrap(), before);

            let cases = [
                (
                    ck::LayeredFingerprintError::Topology(
                        ck::TopologyValidationError::AdjacencyMismatch,
                    ),
                    "TopologyValidationError",
                ),
                (
                    ck::LayeredFingerprintError::Query(ck::QueryGraphError::AdjacencyMismatch),
                    "QueryGraphError",
                ),
                (
                    ck::LayeredFingerprintError::Rings(ck::RingFindingError::Value {
                        message: "ring context",
                    }),
                    "RingFindingError",
                ),
                (
                    ck::LayeredFingerprintError::Paths(ck::PathError::AtomPathOutOfRange {
                        position: usize::MAX,
                        atom: ck::AtomId::new(usize::MAX),
                        atom_count: 3,
                    }),
                    "PathError",
                ),
                (
                    ck::LayeredFingerprintError::Value(
                        ck::FingerprintError::SparseIndexOutOfRange {
                            index: u64::MAX,
                            size: 1,
                        },
                    ),
                    "FingerprintError",
                ),
            ];
            for (native, class) in cases {
                let display = native.to_string();
                let source_display = native.source().unwrap().to_string();
                let error = LayeredFingerprintError::from_rust(native);
                assert_eq!(error.value(py).str()?.to_str()?, display);
                let source = error.value(py).getattr("source")?;
                assert!(source.is(&error.value(py).getattr("__cause__")?));
                assert_eq!(source.str()?.to_str()?, source_display);
                assert_eq!(
                    source
                        .getattr("__class__")?
                        .getattr("__name__")?
                        .extract::<String>()?,
                    class
                );
                assert!(error.value(py).setattr("source", py.None()).is_err());
                if class == "PathError" {
                    assert_eq!(source.getattr("atom")?.extract::<usize>()?, usize::MAX);
                    assert_eq!(source.getattr("position")?.extract::<usize>()?, usize::MAX);
                }
                if class == "FingerprintError" {
                    assert_eq!(source.getattr("index")?.extract::<u64>()?, u64::MAX);
                }
            }
            Ok(())
        })
        .unwrap();
    }
}
