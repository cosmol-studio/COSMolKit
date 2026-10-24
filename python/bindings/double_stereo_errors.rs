//! Exact double-bond stereo errors and their native, nontransparent causes.
use crate::chemical_values::{BondDirection, BondOrder, BondStereo};
use crate::error_framework::{Source, error_projection};
use crate::model_edit_errors::BondValueError;
use crate::model_errors::{CoordinateValidationError, TopologyValidationError};
use crate::property_numeric_errors::PropertyIntReadError;
use crate::property_values::PropertyValueKind;
use crate::ring_errors::RingFindingError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

pub(crate) type DoubleStereoSource = Source<
    crate::error_framework::ErrorObject<RingFindingError>,
    Source<
        crate::error_framework::ErrorObject<PropertyIntReadError>,
        Source<
            crate::error_framework::ErrorObject<TopologyValidationError>,
            Source<
                crate::error_framework::ErrorObject<BondValueError>,
                crate::error_framework::ErrorObject<CoordinateValidationError>,
            >,
        >,
    >,
>;

error_projection!(DoubleBondStereoError, DoubleBondStereoErrorKind, ck::DoubleBondStereoError, "stereo", {
    Rings(_), UnsignedPropertyOverflow { atom, bond, property, value: u32 },
    InvalidPropertyKind { atom, bond, property, kind }, NumericPropertyRead { atom, bond, property, source },
    InvalidTopology(_), BondValue(_), BondOutOfRange { bond, bond_count: usize },
    AtomOutOfRange { atom, atom_count }, NotDoubleBond { bond, order },
    UndefinedStereo { bond, stereo }, UnsupportedStereo { bond, stereo },
    InvalidEndpointDegree { bond, endpoint, degree },
    InvalidDirection { direction: BondDirection => |value: &ck::BondDirection| (*value).into() },
    RankCount { actual, atom_count }, RingInfoNotInitialized,
    RingRowCount { dimension: &'static str, actual, expected: usize }, InvalidConformer(_),
    StereoReferenceRequired { bond, stereo }, StereoReferenceNotNeighbor { bond, endpoint, reference },
    StereoReferenceIsOppositeEndpoint { bond, endpoint, reference },
    ControllingAtomMismatch { bond, endpoint, reference }, ControllingSlotOverflow { bond, endpoint, degree },
}, cause = |inner: &ck::DoubleBondStereoError| match inner {
    ck::DoubleBondStereoError::Rings(value) => Some(RingFindingError::from_rust(value.clone())),
    ck::DoubleBondStereoError::NumericPropertyRead { source, .. } => Some(PropertyIntReadError::from_rust(source.clone())),
    ck::DoubleBondStereoError::InvalidTopology(value) => Some(TopologyValidationError::from_rust(value.clone())),
    ck::DoubleBondStereoError::BondValue(value) => Some(BondValueError::from_rust(value.clone())),
    ck::DoubleBondStereoError::InvalidConformer(value) => Some(CoordinateValidationError::from_rust(value.clone())),
    ck::DoubleBondStereoError::UnsignedPropertyOverflow { .. } | ck::DoubleBondStereoError::InvalidPropertyKind { .. }
    | ck::DoubleBondStereoError::BondOutOfRange { .. } | ck::DoubleBondStereoError::AtomOutOfRange { .. }
    | ck::DoubleBondStereoError::NotDoubleBond { .. } | ck::DoubleBondStereoError::UndefinedStereo { .. }
    | ck::DoubleBondStereoError::UnsupportedStereo { .. } | ck::DoubleBondStereoError::InvalidEndpointDegree { .. }
    | ck::DoubleBondStereoError::InvalidDirection { .. } | ck::DoubleBondStereoError::RankCount { .. }
    | ck::DoubleBondStereoError::RingInfoNotInitialized | ck::DoubleBondStereoError::RingRowCount { .. }
    | ck::DoubleBondStereoError::StereoReferenceRequired { .. } | ck::DoubleBondStereoError::StereoReferenceNotNeighbor { .. }
    | ck::DoubleBondStereoError::StereoReferenceIsOppositeEndpoint { .. }
    | ck::DoubleBondStereoError::ControllingAtomMismatch { .. } | ck::DoubleBondStereoError::ControllingSlotOverflow { .. } => None,
}, source = Option<DoubleStereoSource>, methods = {
    #[getter]
    fn atom(&self) -> Option<usize> { match self.inner {
        ck::DoubleBondStereoError::UnsignedPropertyOverflow { atom, .. }
        | ck::DoubleBondStereoError::InvalidPropertyKind { atom, .. }
        | ck::DoubleBondStereoError::NumericPropertyRead { atom, .. } => atom.map(ck::AtomId::index),
        ck::DoubleBondStereoError::AtomOutOfRange { atom, .. } => Some(atom.index()), _ => None,
    } }
    #[getter]
    fn bond(&self) -> Option<usize> { match self.inner {
        ck::DoubleBondStereoError::UnsignedPropertyOverflow { bond, .. }
        | ck::DoubleBondStereoError::InvalidPropertyKind { bond, .. }
        | ck::DoubleBondStereoError::NumericPropertyRead { bond, .. } => bond.map(ck::BondId::index),
        ck::DoubleBondStereoError::BondOutOfRange { bond, .. } | ck::DoubleBondStereoError::NotDoubleBond { bond, .. }
        | ck::DoubleBondStereoError::UndefinedStereo { bond, .. } | ck::DoubleBondStereoError::UnsupportedStereo { bond, .. }
        | ck::DoubleBondStereoError::InvalidEndpointDegree { bond, .. }
        | ck::DoubleBondStereoError::StereoReferenceRequired { bond, .. }
        | ck::DoubleBondStereoError::StereoReferenceNotNeighbor { bond, .. }
        | ck::DoubleBondStereoError::StereoReferenceIsOppositeEndpoint { bond, .. }
        | ck::DoubleBondStereoError::ControllingAtomMismatch { bond, .. }
        | ck::DoubleBondStereoError::ControllingSlotOverflow { bond, .. } => Some(bond.index()), _ => None,
    } }
    #[getter]
    fn property(&self) -> Option<&'static str> { match self.inner {
        ck::DoubleBondStereoError::UnsignedPropertyOverflow { property, .. }
        | ck::DoubleBondStereoError::InvalidPropertyKind { property, .. }
        | ck::DoubleBondStereoError::NumericPropertyRead { property, .. } => Some(property), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn kind(&self) -> Option<PropertyValueKind> { match self.inner {
        ck::DoubleBondStereoError::InvalidPropertyKind { kind, .. } => Some(kind.into()), _ => None,
    } }
    #[getter]
    fn atom_count(&self) -> Option<usize> { match self.inner {
        ck::DoubleBondStereoError::AtomOutOfRange { atom_count, .. }
        | ck::DoubleBondStereoError::RankCount { atom_count, .. } => Some(atom_count), _ => None,
    } }
    #[getter]
    fn actual(&self) -> Option<usize> { match self.inner {
        ck::DoubleBondStereoError::RankCount { actual, .. }
        | ck::DoubleBondStereoError::RingRowCount { actual, .. } => Some(actual), _ => None,
    } }
    #[getter]
    fn order(&self) -> Option<BondOrder> { match self.inner {
        ck::DoubleBondStereoError::NotDoubleBond { order, .. } => Some(order.into()), _ => None,
    } }
    #[getter]
    fn stereo(&self) -> Option<BondStereo> { match self.inner {
        ck::DoubleBondStereoError::UndefinedStereo { stereo, .. }
        | ck::DoubleBondStereoError::UnsupportedStereo { stereo, .. }
        | ck::DoubleBondStereoError::StereoReferenceRequired { stereo, .. } => Some(stereo.into()), _ => None,
    } }
    #[getter]
    fn endpoint(&self) -> Option<&'static str> { match self.inner {
        ck::DoubleBondStereoError::InvalidEndpointDegree { endpoint, .. }
        | ck::DoubleBondStereoError::StereoReferenceNotNeighbor { endpoint, .. }
        | ck::DoubleBondStereoError::StereoReferenceIsOppositeEndpoint { endpoint, .. }
        | ck::DoubleBondStereoError::ControllingAtomMismatch { endpoint, .. }
        | ck::DoubleBondStereoError::ControllingSlotOverflow { endpoint, .. } => Some(endpoint), _ => None,
    } }
    #[getter]
    fn degree(&self) -> Option<usize> { match self.inner {
        ck::DoubleBondStereoError::InvalidEndpointDegree { degree, .. }
        | ck::DoubleBondStereoError::ControllingSlotOverflow { degree, .. } => Some(degree), _ => None,
    } }
    #[getter]
    fn reference(&self) -> Option<usize> { match self.inner {
        ck::DoubleBondStereoError::StereoReferenceNotNeighbor { reference, .. }
        | ck::DoubleBondStereoError::StereoReferenceIsOppositeEndpoint { reference, .. }
        | ck::DoubleBondStereoError::ControllingAtomMismatch { reference, .. } => Some(reference.index()), _ => None,
    } }
});
