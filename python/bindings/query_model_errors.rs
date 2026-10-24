use crate::attachment_values::TemplateAttachmentOrderError;
use crate::error_framework::{Source, error_projection};
use crate::model_errors::{CoordinateValidationError, StereoGroupError};
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

error_projection!(QueryGraphError, QueryGraphErrorKind, ck::QueryGraphError, "model", {
    StereoGroup(_), AtomIdMismatch { position, id },
    TemplateAttachmentOrder { atom, source }, BondIdMismatch { position, id },
    InvalidBondEndpoint(_), StereoAtomOutOfRange { bond, begin, end, atom_count },
    StereoGroupAtomOutOfRange { atom, atom_count },
    StereoGroupBondOutOfRange { bond, bond_count },
    SubstanceGroupIdMismatch { position, id },
    SubstanceGroupAtomOutOfRange { sgroup, atom, atom_count },
    SubstanceGroupBondOutOfRange { sgroup, bond, bond_count },
    SubstanceGroupParentOutOfRange { sgroup, parent },
    CoordinateValidation(_), AdjacencyMismatch,
}, cause = |inner: &ck::QueryGraphError| match inner {
    ck::QueryGraphError::StereoGroup(source) => Some(StereoGroupError::from_rust(*source)),
    ck::QueryGraphError::TemplateAttachmentOrder { source, .. } => Some(TemplateAttachmentOrderError::from_rust(source.clone())),
    ck::QueryGraphError::AtomIdMismatch { .. } | ck::QueryGraphError::BondIdMismatch { .. }
    | ck::QueryGraphError::InvalidBondEndpoint(_) | ck::QueryGraphError::StereoAtomOutOfRange { .. }
    | ck::QueryGraphError::StereoGroupAtomOutOfRange { .. } | ck::QueryGraphError::StereoGroupBondOutOfRange { .. }
    | ck::QueryGraphError::SubstanceGroupIdMismatch { .. } | ck::QueryGraphError::SubstanceGroupAtomOutOfRange { .. }
    | ck::QueryGraphError::SubstanceGroupBondOutOfRange { .. } | ck::QueryGraphError::SubstanceGroupParentOutOfRange { .. }
    | ck::QueryGraphError::CoordinateValidation(_) | ck::QueryGraphError::AdjacencyMismatch => None,
}, source = Option<Source<crate::error_framework::ErrorObject<StereoGroupError>, crate::error_framework::ErrorObject<TemplateAttachmentOrderError>>>, methods = {
    #[getter]
    fn position(&self) -> Option<usize> {
        match self.inner {
            ck::QueryGraphError::AtomIdMismatch { position, .. }
            | ck::QueryGraphError::BondIdMismatch { position, .. }
            | ck::QueryGraphError::SubstanceGroupIdMismatch { position, .. } => Some(position),
            _ => None,
        }
    }
    #[getter]
    fn id(&self) -> Option<usize> {
        match self.inner {
            ck::QueryGraphError::AtomIdMismatch { id, .. } => Some(id.index()),
            ck::QueryGraphError::BondIdMismatch { id, .. } => Some(id.index()),
            ck::QueryGraphError::SubstanceGroupIdMismatch { id, .. } => Some(id.index()),
            _ => None,
        }
    }
    #[getter]
    fn atom(&self) -> Option<usize> {
        match self.inner {
            ck::QueryGraphError::TemplateAttachmentOrder { atom, .. }
            | ck::QueryGraphError::StereoGroupAtomOutOfRange { atom, .. }
            | ck::QueryGraphError::SubstanceGroupAtomOutOfRange { atom, .. } => Some(atom.index()),
            _ => None,
        }
    }
    #[getter]
    fn bond(&self) -> Option<usize> {
        match self.inner {
            ck::QueryGraphError::InvalidBondEndpoint(bond) => Some(bond),
            ck::QueryGraphError::StereoAtomOutOfRange { bond, .. }
            | ck::QueryGraphError::StereoGroupBondOutOfRange { bond, .. }
            | ck::QueryGraphError::SubstanceGroupBondOutOfRange { bond, .. } => Some(bond.index()),
            _ => None,
        }
    }
    #[getter]
    fn begin(&self) -> Option<usize> {
        match self.inner { ck::QueryGraphError::StereoAtomOutOfRange { begin, .. } => Some(begin.index()), _ => None }
    }
    #[getter]
    fn end(&self) -> Option<usize> {
        match self.inner { ck::QueryGraphError::StereoAtomOutOfRange { end, .. } => Some(end.index()), _ => None }
    }
    #[getter]
    fn atom_count(&self) -> Option<usize> {
        match self.inner {
            ck::QueryGraphError::StereoAtomOutOfRange { atom_count, .. }
            | ck::QueryGraphError::StereoGroupAtomOutOfRange { atom_count, .. }
            | ck::QueryGraphError::SubstanceGroupAtomOutOfRange { atom_count, .. } => Some(atom_count),
            _ => None,
        }
    }
    #[getter]
    fn bond_count(&self) -> Option<usize> {
        match self.inner {
            ck::QueryGraphError::StereoGroupBondOutOfRange { bond_count, .. }
            | ck::QueryGraphError::SubstanceGroupBondOutOfRange { bond_count, .. } => Some(bond_count),
            _ => None,
        }
    }
    #[getter]
    fn sgroup(&self) -> Option<usize> {
        match self.inner {
            ck::QueryGraphError::SubstanceGroupAtomOutOfRange { sgroup, .. }
            | ck::QueryGraphError::SubstanceGroupBondOutOfRange { sgroup, .. }
            | ck::QueryGraphError::SubstanceGroupParentOutOfRange { sgroup, .. } => Some(sgroup.index()),
            _ => None,
        }
    }
    #[getter]
    fn parent(&self) -> Option<usize> {
        match self.inner { ck::QueryGraphError::SubstanceGroupParentOutOfRange { parent, .. } => Some(parent.index()), _ => None }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn coordinate_validation(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<CoordinateValidationError>>> {
        match &self.inner {
            // This native tuple is a payload, not Error::source. Do not invent a cause.
            ck::QueryGraphError::CoordinateValidation(value) => CoordinateValidationError::from_rust(value.clone()).value(py).extract::<crate::error_framework::ErrorObject<CoordinateValidationError>>().map(Some).map_err(Into::into),
            _ => Ok(None),
        }
    }
});
