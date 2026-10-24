//! Concrete model error dependencies. They are exported only when a registered
//! parent's actual typed projection reaches them; no second API inventory.
use crate::attachment_values::TemplateAttachmentOrderError;
use crate::chemical_values::BondStereo;
use crate::error_framework::{Source, error_projection};
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

error_projection!(StereoGroupError, StereoGroupErrorKind, ck::StereoGroupError, "model", {
    DuplicateAtom, DuplicateBond,
}, cause = |_: &ck::StereoGroupError| None::<PyErr>);

error_projection!(CoordinateValidationError, CoordinateValidationErrorKind, ck::CoordinateValidationError, "coordinates", {
    AtomIndexOverflow { atom },
    Missing3DConformer { id },
    ConformerIdOverflow { max_id: usize },
    RowCount { dimension, conformer, rows: usize, atom_count: usize },
    MissingSourceConformerOrder,
    SourceConformerOrder { two_d: usize, three_d: usize, expected_two_d: usize, expected_three_d: usize },
    DuplicateConformerId { dimension, id },
    NonFiniteCoordinate { dimension, conformer, atom, axis: &'static str },
}, cause = |_: &ck::CoordinateValidationError| None::<PyErr>, methods = {
    #[getter]
    fn atom(&self) -> Option<usize> {
        match self.inner {
            ck::CoordinateValidationError::AtomIndexOverflow { atom }
            | ck::CoordinateValidationError::NonFiniteCoordinate { atom, .. } => Some(atom),
            _ => None,
        }
    }
    #[getter]
    fn id(&self) -> Option<usize> {
        match self.inner {
            ck::CoordinateValidationError::Missing3DConformer { id }
            | ck::CoordinateValidationError::DuplicateConformerId { id, .. } => Some(id),
            _ => None,
        }
    }
    #[getter]
    fn dimension(&self) -> Option<&'static str> {
        match self.inner {
            ck::CoordinateValidationError::RowCount { dimension, .. }
            | ck::CoordinateValidationError::DuplicateConformerId { dimension, .. }
            | ck::CoordinateValidationError::NonFiniteCoordinate { dimension, .. } => Some(dimension),
            _ => None,
        }
    }
    #[getter]
    fn conformer(&self) -> Option<usize> {
        match self.inner {
            ck::CoordinateValidationError::RowCount { conformer, .. }
            | ck::CoordinateValidationError::NonFiniteCoordinate { conformer, .. } => Some(conformer),
            _ => None,
        }
    }
});

error_projection!(TopologyValidationError, TopologyValidationErrorKind, ck::TopologyValidationError, "model", {
    StereoGroup(_),
    AtomIdMismatch { position, id },
    TemplateAttachmentOrder { atom, source },
    BondIdMismatch { position, id },
    BondEndpointOutOfRange { bond, endpoint: &'static str, atom, atom_count },
    SelfLoopBond { bond, atom },
    StereoAtomOutOfRange { bond, begin: usize => |id: &ck::AtomId| id.index(), end: usize => |id: &ck::AtomId| id.index(), atom_count },
    StereoReferenceOutOfRange { bond, atom, atom_count },
    StereoAtomsRequired { bond, stereo: BondStereo => |value: &ck::BondStereo| (*value).into() },
    SubstanceGroupIdMismatch { position, id },
    SubstanceGroupAtomOutOfRange { sgroup, atom, atom_count },
    SubstanceGroupBondOutOfRange { sgroup, bond, bond_count },
    SubstanceGroupParentOutOfRange { sgroup, parent: usize => |id: &ck::SubstanceGroupId| id.index() },
    StereoGroupAtomOutOfRange { atom, atom_count },
    StereoGroupBondOutOfRange { bond, bond_count },
    AdjacencyMismatch,
}, cause = |inner: &ck::TopologyValidationError| match inner {
    ck::TopologyValidationError::StereoGroup(source) => Some(StereoGroupError::from_rust(*source)),
    ck::TopologyValidationError::TemplateAttachmentOrder { source, .. } => Some(TemplateAttachmentOrderError::from_rust(source.clone())),
    ck::TopologyValidationError::AtomIdMismatch { .. }
    | ck::TopologyValidationError::BondIdMismatch { .. }
    | ck::TopologyValidationError::BondEndpointOutOfRange { .. }
    | ck::TopologyValidationError::SelfLoopBond { .. }
    | ck::TopologyValidationError::StereoAtomOutOfRange { .. }
    | ck::TopologyValidationError::StereoReferenceOutOfRange { .. }
    | ck::TopologyValidationError::StereoAtomsRequired { .. }
    | ck::TopologyValidationError::SubstanceGroupIdMismatch { .. }
    | ck::TopologyValidationError::SubstanceGroupAtomOutOfRange { .. }
    | ck::TopologyValidationError::SubstanceGroupBondOutOfRange { .. }
    | ck::TopologyValidationError::SubstanceGroupParentOutOfRange { .. }
    | ck::TopologyValidationError::StereoGroupAtomOutOfRange { .. }
    | ck::TopologyValidationError::StereoGroupBondOutOfRange { .. }
    | ck::TopologyValidationError::AdjacencyMismatch => None,
}, source = Option<Source<crate::error_framework::ErrorObject<StereoGroupError>, crate::error_framework::ErrorObject<TemplateAttachmentOrderError>>>, methods = {
    #[getter]
    fn position(&self) -> Option<usize> {
        match self.inner {
            ck::TopologyValidationError::AtomIdMismatch { position, .. }
            | ck::TopologyValidationError::BondIdMismatch { position, .. }
            | ck::TopologyValidationError::SubstanceGroupIdMismatch { position, .. } => Some(position),
            _ => None,
        }
    }
    #[getter]
    fn id(&self) -> Option<usize> {
        match self.inner {
            ck::TopologyValidationError::AtomIdMismatch { id, .. } => Some(id.index()),
            ck::TopologyValidationError::BondIdMismatch { id, .. } => Some(id.index()),
            ck::TopologyValidationError::SubstanceGroupIdMismatch { id, .. } => Some(id.index()),
            _ => None,
        }
    }
    #[getter]
    fn atom(&self) -> Option<usize> {
        match self.inner {
            ck::TopologyValidationError::TemplateAttachmentOrder { atom, .. }
            | ck::TopologyValidationError::BondEndpointOutOfRange { atom, .. }
            | ck::TopologyValidationError::SelfLoopBond { atom, .. }
            | ck::TopologyValidationError::StereoReferenceOutOfRange { atom, .. }
            | ck::TopologyValidationError::SubstanceGroupAtomOutOfRange { atom, .. }
            | ck::TopologyValidationError::StereoGroupAtomOutOfRange { atom, .. } => Some(atom.index()),
            _ => None,
        }
    }
    #[getter]
    fn bond(&self) -> Option<usize> {
        match self.inner {
            ck::TopologyValidationError::BondEndpointOutOfRange { bond, .. }
            | ck::TopologyValidationError::SelfLoopBond { bond, .. }
            | ck::TopologyValidationError::StereoAtomOutOfRange { bond, .. }
            | ck::TopologyValidationError::StereoReferenceOutOfRange { bond, .. }
            | ck::TopologyValidationError::StereoAtomsRequired { bond, .. }
            | ck::TopologyValidationError::SubstanceGroupBondOutOfRange { bond, .. }
            | ck::TopologyValidationError::StereoGroupBondOutOfRange { bond, .. } => Some(bond.index()),
            _ => None,
        }
    }
    #[getter]
    fn atom_count(&self) -> Option<usize> {
        match self.inner {
            ck::TopologyValidationError::BondEndpointOutOfRange { atom_count, .. }
            | ck::TopologyValidationError::StereoAtomOutOfRange { atom_count, .. }
            | ck::TopologyValidationError::StereoReferenceOutOfRange { atom_count, .. }
            | ck::TopologyValidationError::SubstanceGroupAtomOutOfRange { atom_count, .. }
            | ck::TopologyValidationError::StereoGroupAtomOutOfRange { atom_count, .. } => Some(atom_count),
            _ => None,
        }
    }
    #[getter]
    fn bond_count(&self) -> Option<usize> {
        match self.inner {
            ck::TopologyValidationError::SubstanceGroupBondOutOfRange { bond_count, .. }
            | ck::TopologyValidationError::StereoGroupBondOutOfRange { bond_count, .. } => Some(bond_count),
            _ => None,
        }
    }
    #[getter]
    fn sgroup(&self) -> Option<usize> {
        match self.inner {
            ck::TopologyValidationError::SubstanceGroupAtomOutOfRange { sgroup, .. }
            | ck::TopologyValidationError::SubstanceGroupBondOutOfRange { sgroup, .. }
            | ck::TopologyValidationError::SubstanceGroupParentOutOfRange { sgroup, .. } => Some(sgroup.index()),
            _ => None,
        }
    }
});

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn model_errors_keep_actual_validation_fields_and_sources() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let carbon = ck::Element::from_atomic_number(6).unwrap();
            let bad_atom = ck::Atom::from_spec(ck::AtomId::new(3), ck::AtomSpec::new(carbon));
            let native = ck::TopologyBlock::try_from_parts(vec![bad_atom], vec![], vec![], vec![])
                .unwrap_err();
            let display = native.to_string();
            let error = TopologyValidationError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert_eq!(error.value(py).getattr("id")?.extract::<usize>()?, 3);
            assert_eq!(error.value(py).getattr("position")?.extract::<usize>()?, 0);
            assert!(error.cause(py).is_none());
            let error = TopologyValidationError::from_rust(
                ck::TopologyValidationError::TemplateAttachmentOrder {
                    atom: ck::AtomId::new(7),
                    source: ck::TemplateAttachmentOrderError::TargetRemoved {
                        position: 4,
                        target: ck::AtomId::new(9),
                    },
                },
            );
            let cause = error.cause(py).unwrap();
            assert!(cause.is_instance(py, <TemplateAttachmentOrderError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
            assert!(error.value(py).getattr("source")?.is(cause.value(py)));
            assert_eq!(cause.value(py).getattr("position")?.extract::<usize>()?, 4);
            assert_eq!(cause.value(py).getattr("target")?.extract::<usize>()?, 9);
            let error = TopologyValidationError::from_rust(
                ck::TopologyValidationError::StereoGroup(ck::StereoGroupError::DuplicateBond),
            );
            assert!(
                error
                    .cause(py)
                    .unwrap()
                    .is_instance(py, <StereoGroupError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())
            );
            let error =
                CoordinateValidationError::from_rust(ck::CoordinateValidationError::RowCount {
                    dimension: "3D",
                    conformer: 11,
                    rows: 2,
                    atom_count: 5,
                });
            assert_eq!(
                error.value(py).getattr("conformer")?.extract::<usize>()?,
                11
            );
            assert_eq!(error.value(py).getattr("rows")?.extract::<usize>()?, 2);
            assert_eq!(
                error.value(py).getattr("atom_count")?.extract::<usize>()?,
                5
            );
            assert!(error.cause(py).is_none());
            Ok(())
        })
        .unwrap();
    }
}
