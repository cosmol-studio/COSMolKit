//! Concrete ring failures retained by tautomer inspection and score queries.
use crate::attachment_values::TemplateAttachmentOrderError;
use crate::error_framework::{Source, error_projection};
use crate::model_edit_errors::{AdjacencyError, MoleculePropertyError};
use crate::model_errors::{StereoGroupError, TopologyValidationError};
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

error_projection!(RingDecomposerError, RingDecomposerErrorKind, ck::RingDecomposerError, "rings", {
    NodeOutOfRange { node, node_count: usize }, SelfLoop { node },
    DuplicateEdge { from, to }, EmptyGraph, EdgeNotFound { from, to },
    RelevantCycleCountOverflow, InvalidShortestPathDag { root: usize, target: usize },
}, cause = |_: &ck::RingDecomposerError| None::<PyErr>, methods = {
    #[getter]
    #[doc = "Graph node involved in the decomposition failure. Present for `RingDecomposerError.NodeOutOfRange`, `RingDecomposerError.SelfLoop`; None when unavailable."] fn node(&self) -> Option<usize> { match self.inner {
        ck::RingDecomposerError::NodeOutOfRange { node, .. }
        | ck::RingDecomposerError::SelfLoop { node } => Some(node), _ => None,
    } }
    /// Rust `from` is escaped only because it is a Python keyword.
    #[getter]
    fn from_(&self) -> Option<usize> { match self.inner {
        ck::RingDecomposerError::DuplicateEdge { from, .. }
        | ck::RingDecomposerError::EdgeNotFound { from, .. } => Some(from), _ => None,
    } }
    #[getter]
    #[doc = "Destination graph node involved in the failure. Present for `RingDecomposerError.DuplicateEdge`, `RingDecomposerError.EdgeNotFound`; None when unavailable."] fn to(&self) -> Option<usize> { match self.inner {
        ck::RingDecomposerError::DuplicateEdge { to, .. }
        | ck::RingDecomposerError::EdgeNotFound { to, .. } => Some(to), _ => None,
    } }
});

pub(crate) type RingSource = Source<
    crate::error_framework::ErrorObject<MoleculePropertyError>,
    Source<
        crate::error_framework::ErrorObject<StereoGroupError>,
        crate::error_framework::ErrorObject<TemplateAttachmentOrderError>,
    >,
>;

error_projection!(RingFindingError, RingFindingErrorKind, ck::RingFindingError, "rings", {
    SourceRingFindType { tag: u32 }, MoleculeProperty(_), Value { message: &'static str },
    ExpectedBondNotFound { begin: usize => |value: &ck::AtomId| value.index(), end: usize => |value: &ck::AtomId| value.index() },
    RingAtomOutOfRange { atom: usize, atom_count: usize },
    RingBondOutOfRange { bond: usize, bond_count: usize },
    RingFamilyEdgeMappingMissing { edge: usize }, RelevantCycleCountOutOfRange,
    UnsupportedBranch { reason: &'static str }, Adjacency(_), InvalidTopology(_), RingDecomposer(_),
}, cause = |inner: &ck::RingFindingError| Python::attach(|py| match inner {
    ck::RingFindingError::MoleculeProperty(value) => Some(MoleculePropertyError::from_rust(value.clone())),
    ck::RingFindingError::InvalidTopology(value) => TopologyValidationError::from_rust(value.clone()).cause(py),
    ck::RingFindingError::SourceRingFindType { .. } | ck::RingFindingError::Value { .. }
    | ck::RingFindingError::ExpectedBondNotFound { .. } | ck::RingFindingError::RingAtomOutOfRange { .. }
    | ck::RingFindingError::RingBondOutOfRange { .. } | ck::RingFindingError::RingFamilyEdgeMappingMissing { .. }
    | ck::RingFindingError::RelevantCycleCountOutOfRange | ck::RingFindingError::UnsupportedBranch { .. }
    | ck::RingFindingError::Adjacency(_) | ck::RingFindingError::RingDecomposer(_) => None,
}), source = Option<RingSource>, methods = {
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured AdjacencyError payload of RingFindingError.Adjacency; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn adjacency(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<AdjacencyError>>> {
        match &self.inner { ck::RingFindingError::Adjacency(value) =>
            AdjacencyError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None) }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured TopologyValidationError payload of RingFindingError.InvalidTopology; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn topology(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<TopologyValidationError>>> {
        match &self.inner { ck::RingFindingError::InvalidTopology(value) =>
            TopologyValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None) }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured RingDecomposerError payload of RingFindingError.RingDecomposer; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn ring_decomposer(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<RingDecomposerError>>> {
        match &self.inner { ck::RingFindingError::RingDecomposer(value) =>
            RingDecomposerError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None) }
    }
});

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn ring_transparent_payload_preserves_exact_native_variant_and_keyword_field() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native =
                ck::RingFindingError::RingDecomposer(ck::RingDecomposerError::DuplicateEdge {
                    from: 2,
                    to: 7,
                });
            assert!(native.source().is_none());
            let display = native.to_string();
            let error = RingFindingError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert!(error.cause(py).is_none());
            let payload = error.value(py).getattr("ring_decomposer")?;
            assert!(payload.is_instance(<RingDecomposerError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert_eq!(payload.getattr("from_")?.extract::<usize>()?, 2);
            assert_eq!(payload.getattr("to")?.extract::<usize>()?, 7);
            assert!(
                payload
                    .getattr("variant")?
                    .extract::<RingDecomposerErrorKind>()?
                    == RingDecomposerErrorKind::DuplicateEdge
            );
            assert!(payload.getattr("__cause__")?.is_none());
            assert!(payload.setattr("from_", 9).is_err());
            Ok(())
        })
        .unwrap();
    }
}
