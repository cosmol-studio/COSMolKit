//! Exact stage-tagged native sanitize errors, without a second failure policy.
use crate::error_framework::{Source, error_projection};
use crate::kekulize_errors::KekulizeError;
use crate::model_edit_errors::{
    AtomPropertyError, BondValueError, MoleculePropertyError, QueryStateError,
};
use crate::model_errors::TopologyValidationError;
use crate::ring_errors::RingFindingError;
use crate::sanitize_leaf_errors::{AdjustHsError, HybridizationError, RadicalError};
use crate::sanitize_support_errors::{
    AromaticityError, CleanupError, ConjugationError, PropertyCacheError,
};
use crate::sanitize_values::SanitizeStage;
use crate::stereo_errors::{AtropisomerError, StereoError};
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::error::Error;

pub(crate) type SanitizeSource = Source<
    Py<MoleculePropertyError>,
    Source<
        Py<BondValueError>,
        Source<
            Py<AtomPropertyError>,
            Source<
                Py<TopologyValidationError>,
                Source<
                    Py<QueryStateError>,
                    Source<
                        Py<CleanupError>,
                        Source<
                            Py<PropertyCacheError>,
                            Source<
                                Py<RingFindingError>,
                                Source<
                                    Py<KekulizeError>,
                                    Source<
                                        Py<RadicalError>,
                                        Source<
                                            Py<AromaticityError>,
                                            Source<
                                                Py<ConjugationError>,
                                                Source<
                                                    Py<HybridizationError>,
                                                    Source<
                                                        Py<AtropisomerError>,
                                                        Source<Py<StereoError>, Py<AdjustHsError>>,
                                                    >,
                                                >,
                                            >,
                                        >,
                                    >,
                                >,
                            >,
                        >,
                    >,
                >,
            >,
        >,
    >,
>;

error_projection!(SanitizeError, SanitizeErrorKind, ck::SanitizeError, "sanitize", {
    MoleculeProperty(_), BondProperty(_), AtomProperty(_),
    InvalidOperations { bits: u32, unknown_bits: u32 },
    InvalidTopology { stage, source }, InvalidQueryState { stage, source },
    Cleanup { stage, source }, Properties { stage, source }, Rings { stage, source },
    Kekulize { stage, source }, Radicals { stage, source }, Aromaticity { stage, source },
    Conjugation { stage, source }, Hybridization { stage, source },
    Atropisomers { stage, source }, Chirality { stage, source }, AdjustHs { stage, source },
}, cause = |inner: &ck::SanitizeError| match inner {
    ck::SanitizeError::MoleculeProperty(value) => Some(MoleculePropertyError::from_rust(value.clone())),
    ck::SanitizeError::BondProperty(value) => Some(BondValueError::from_rust(value.clone())),
    ck::SanitizeError::AtomProperty(value) => Some(AtomPropertyError::from_rust(value.clone())),
    ck::SanitizeError::InvalidTopology { source, .. } => Some(TopologyValidationError::from_rust(source.clone())),
    ck::SanitizeError::InvalidQueryState { source, .. } => Some(QueryStateError::from_rust(source.clone())),
    ck::SanitizeError::Cleanup { source, .. } => Some(CleanupError::from_rust(source.clone())),
    ck::SanitizeError::Properties { source, .. } => Some(PropertyCacheError::from_rust(source.clone())),
    ck::SanitizeError::Rings { source, .. } => Some(RingFindingError::from_rust(source.clone())),
    ck::SanitizeError::Kekulize { source, .. } => Some(KekulizeError::from_rust(source.clone())),
    ck::SanitizeError::Radicals { source, .. } => Some(RadicalError::from_rust(source.clone())),
    ck::SanitizeError::Aromaticity { source, .. } => Some(AromaticityError::from_rust(source.clone())),
    ck::SanitizeError::Conjugation { source, .. } => Some(ConjugationError::from_rust(source.clone())),
    ck::SanitizeError::Hybridization { source, .. } => Some(HybridizationError::from_rust(source.clone())),
    ck::SanitizeError::Atropisomers { source, .. } => Some(AtropisomerError::from_rust(source.clone())),
    ck::SanitizeError::Chirality { source, .. } => Some(StereoError::from_rust(source.clone())),
    ck::SanitizeError::AdjustHs { source, .. } => Some(AdjustHsError::from_rust(source.clone())),
    ck::SanitizeError::InvalidOperations { .. } => None,
}, source = Option<SanitizeSource>, methods = {
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn stage(&self) -> Option<SanitizeStage> { match self.inner {
        ck::SanitizeError::InvalidTopology { stage, .. }
        | ck::SanitizeError::InvalidQueryState { stage, .. }
        | ck::SanitizeError::Cleanup { stage, .. }
        | ck::SanitizeError::Properties { stage, .. }
        | ck::SanitizeError::Rings { stage, .. }
        | ck::SanitizeError::Kekulize { stage, .. }
        | ck::SanitizeError::Radicals { stage, .. }
        | ck::SanitizeError::Aromaticity { stage, .. }
        | ck::SanitizeError::Conjugation { stage, .. }
        | ck::SanitizeError::Hybridization { stage, .. }
        | ck::SanitizeError::Atropisomers { stage, .. }
        | ck::SanitizeError::Chirality { stage, .. }
        | ck::SanitizeError::AdjustHs { stage, .. } => Some(stage.into()), _ => None,
    } }
});

#[cfg(test)]
mod tests {
    use super::*;
    use crate::valence_errors::ValenceError;

    #[test]
    fn sanitize_errors_preserve_stage_and_each_real_cause_level() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native = ck::SanitizeOperations::from_bits(1 << 20).unwrap_err();
            let message = native.to_string();
            let error = SanitizeError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, message);
            assert_eq!(error.value(py).getattr("bits")?.extract::<u32>()?, 1 << 20);
            assert_eq!(
                error.value(py).getattr("unknown_bits")?.extract::<u32>()?,
                1 << 20
            );
            assert!(error.value(py).getattr("stage")?.is_none());
            assert!(error.cause(py).is_none());
            assert!(error.value(py).setattr("bits", 0).is_err());

            let native = ck::SanitizeError::Properties {
                stage: ck::SanitizeStage::Properties,
                source: ck::PropertyCacheError::Valence(ck::ValenceError::InvalidTopology {
                    source: ck::TopologyValidationError::AdjacencyMismatch,
                }),
            };
            assert!(native.source().unwrap().is::<ck::PropertyCacheError>());
            let error = SanitizeError::from_rust(native);
            assert!(
                error
                    .value(py)
                    .getattr("stage")?
                    .extract::<SanitizeStage>()?
                    == SanitizeStage::Properties
            );
            let source = error.value(py).getattr("source")?;
            assert!(source.is(&error.value(py).getattr("__cause__")?));
            assert!(source.is_instance_of::<PropertyCacheError>());
            let valence = source.getattr("source")?;
            assert!(valence.is_instance_of::<ValenceError>());
            assert!(
                valence
                    .getattr("source")?
                    .is_instance_of::<TopologyValidationError>()
            );
            Ok(())
        })
        .unwrap();
    }
}
