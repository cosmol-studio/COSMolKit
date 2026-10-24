//! CIP labeler errors retain the original property kind and concrete cause.
use crate::{
    atropisomer_errors, double_stereo_errors,
    enums::{BondOrder, ChiralTag},
    errors::{enum_error, exception},
    kekulize_errors,
    model_errors::{
        AtomPropertyError, BondValueError, MoleculePropertyError, PropertyValueKind,
        TopologyValidationError,
    },
    periodic_errors::PeriodicTableError,
    property_errors::PropertyUIntReadError,
    ring_errors, valence_errors,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "StereoGroupError | ValenceError | ConjugationError | HybridizationError | RingFindingError | TopologyValidationError | CoordinateValidationError | BondValueError | AtomPropertyError | PropertyIntReadError | TemplateAttachmentOrderError | MappingValidationError | PropertyUIntReadError | MoleculePropertyError"
    )]
    #[derive(Clone)]
    pub type CipLabelerCause;
    #[wasm_bindgen(typescript_type = "PropertyValueKind | string")]
    pub type CipPayloadKind;
    #[wasm_bindgen(typescript_type = "AtropisomerError")]
    pub type CipAtropisomer;
    #[wasm_bindgen(typescript_type = "DoubleBondStereoError")]
    pub type CipDoubleStereo;
    #[wasm_bindgen(typescript_type = "RingFindingError")]
    pub type CipRings;
    #[wasm_bindgen(typescript_type = "ValenceError")]
    pub type CipValence;
    #[wasm_bindgen(typescript_type = "KekulizeError")]
    pub type CipKekulize;
    #[wasm_bindgen(typescript_type = "PeriodicTableError")]
    pub type CipPeriodicTable;
}
pub(crate) fn cause(value: &ck::CipLabelerError) -> Option<CipLabelerCause> {
    match value {
        ck::CipLabelerError::UnsignedPropertyRead { source, .. } => {
            Some(PropertyUIntReadError::from_rust(source.clone()).unchecked_into())
        }
        ck::CipLabelerError::InvalidTopology(value) => {
            Some(TopologyValidationError::from_rust(value.clone()).unchecked_into())
        }
        ck::CipLabelerError::AtomProperty { source, .. } => {
            Some(AtomPropertyError::from_rust(source.clone()).unchecked_into())
        }
        ck::CipLabelerError::BondValue { source, .. } => {
            Some(BondValueError::from_rust(source.clone()).unchecked_into())
        }
        ck::CipLabelerError::MoleculeProperty(value) => {
            Some(MoleculePropertyError::from_rust(value.clone()).unchecked_into())
        }
        ck::CipLabelerError::Atropisomer(value) => {
            atropisomer_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::CipLabelerError::DoubleStereo(value) => {
            double_stereo_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::CipLabelerError::RingFinding(value) => {
            ring_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::CipLabelerError::Valence(value) => {
            valence_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::CipLabelerError::Kekulize(value) => {
            kekulize_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::CipLabelerError::PeriodicTable(_)
        | ck::CipLabelerError::InvalidPropertyKind { .. }
        | ck::CipLabelerError::AtomIndexOutOfRange { .. }
        | ck::CipLabelerError::BondIndexOutOfRange { .. }
        | ck::CipLabelerError::SourceIndexWidthExceeded { .. }
        | ck::CipLabelerError::BondNotIncident { .. }
        | ck::CipLabelerError::NonIntegerBondOrder { .. }
        | ck::CipLabelerError::EdgeEndpointMismatch { .. }
        | ck::CipLabelerError::TooManyNodes { .. }
        | ck::CipLabelerError::InvalidInternalState { .. }
        | ck::CipLabelerError::UnsupportedConfiguration { .. }
        | ck::CipLabelerError::MaxIterationsExceeded
        | ck::CipLabelerError::ControlCCaught
        | ck::CipLabelerError::UnexpectedUpEdgeOrdering
        | ck::CipLabelerError::NoSequenceRuleProvided
        | ck::CipLabelerError::EmptyPairListReference
        | ck::CipLabelerError::DescriptorListLengthMismatch
        | ck::CipLabelerError::InvalidStereoDescriptor
        | ck::CipLabelerError::SubstituentsShouldBeTopologicallyEquivalent
        | ck::CipLabelerError::SomethingUnexpected
        | ck::CipLabelerError::Rule4bInstanceNotInRuleSet
        | ck::CipLabelerError::Rule5NewInstanceNotInRuleSet
        | ck::CipLabelerError::ParityVectorsMustHaveSize4
        | ck::CipLabelerError::EmptyConfigurationFoci
        | ck::CipLabelerError::BadTetrahedralConfig
        | ck::CipLabelerError::TetrahedralConfigurationMustHave4Carriers
        | ck::CipLabelerError::DescriptorNotSupportedForAtoms
        | ck::CipLabelerError::InvalidAtomDescriptor
        | ck::CipLabelerError::CarrierMismatch
        | ck::CipLabelerError::BadSp2BondFoci
        | ck::CipLabelerError::BadSp2BondConfig
        | ck::CipLabelerError::IncorrectNumberOfStereoAtoms
        | ck::CipLabelerError::DescriptorNotSupportedForDoubleBonds
        | ck::CipLabelerError::InvalidBondDescriptor
        | ck::CipLabelerError::BadAtropisomerBondFoci
        | ck::CipLabelerError::BadAtropisomerBondConfig
        | ck::CipLabelerError::DescriptorNotSupportedForAtropisomerBonds => None,
    }
}
enum_error!(CipLabelerError, CipLabelerErrorKind, ck::CipLabelerError, {
    UnsignedPropertyRead => ck::CipLabelerError::UnsignedPropertyRead { .. }, InvalidPropertyKind => ck::CipLabelerError::InvalidPropertyKind { .. },
    AtomIndexOutOfRange => ck::CipLabelerError::AtomIndexOutOfRange { .. }, BondIndexOutOfRange => ck::CipLabelerError::BondIndexOutOfRange { .. },
    SourceIndexWidthExceeded => ck::CipLabelerError::SourceIndexWidthExceeded { .. }, BondNotIncident => ck::CipLabelerError::BondNotIncident { .. },
    NonIntegerBondOrder => ck::CipLabelerError::NonIntegerBondOrder { .. }, EdgeEndpointMismatch => ck::CipLabelerError::EdgeEndpointMismatch { .. }, TooManyNodes => ck::CipLabelerError::TooManyNodes { .. },
    MaxIterationsExceeded => ck::CipLabelerError::MaxIterationsExceeded, ControlCCaught => ck::CipLabelerError::ControlCCaught, UnexpectedUpEdgeOrdering => ck::CipLabelerError::UnexpectedUpEdgeOrdering,
    NoSequenceRuleProvided => ck::CipLabelerError::NoSequenceRuleProvided, EmptyPairListReference => ck::CipLabelerError::EmptyPairListReference,
    DescriptorListLengthMismatch => ck::CipLabelerError::DescriptorListLengthMismatch, InvalidStereoDescriptor => ck::CipLabelerError::InvalidStereoDescriptor,
    SubstituentsShouldBeTopologicallyEquivalent => ck::CipLabelerError::SubstituentsShouldBeTopologicallyEquivalent, SomethingUnexpected => ck::CipLabelerError::SomethingUnexpected,
    Rule4bInstanceNotInRuleSet => ck::CipLabelerError::Rule4bInstanceNotInRuleSet, Rule5NewInstanceNotInRuleSet => ck::CipLabelerError::Rule5NewInstanceNotInRuleSet,
    ParityVectorsMustHaveSize4 => ck::CipLabelerError::ParityVectorsMustHaveSize4, EmptyConfigurationFoci => ck::CipLabelerError::EmptyConfigurationFoci,
    BadTetrahedralConfig => ck::CipLabelerError::BadTetrahedralConfig, TetrahedralConfigurationMustHave4Carriers => ck::CipLabelerError::TetrahedralConfigurationMustHave4Carriers,
    DescriptorNotSupportedForAtoms => ck::CipLabelerError::DescriptorNotSupportedForAtoms, InvalidAtomDescriptor => ck::CipLabelerError::InvalidAtomDescriptor,
    CarrierMismatch => ck::CipLabelerError::CarrierMismatch, BadSp2BondFoci => ck::CipLabelerError::BadSp2BondFoci, BadSp2BondConfig => ck::CipLabelerError::BadSp2BondConfig,
    IncorrectNumberOfStereoAtoms => ck::CipLabelerError::IncorrectNumberOfStereoAtoms, DescriptorNotSupportedForDoubleBonds => ck::CipLabelerError::DescriptorNotSupportedForDoubleBonds,
    InvalidBondDescriptor => ck::CipLabelerError::InvalidBondDescriptor, BadAtropisomerBondFoci => ck::CipLabelerError::BadAtropisomerBondFoci,
    BadAtropisomerBondConfig => ck::CipLabelerError::BadAtropisomerBondConfig, DescriptorNotSupportedForAtropisomerBonds => ck::CipLabelerError::DescriptorNotSupportedForAtropisomerBonds,
    InvalidInternalState => ck::CipLabelerError::InvalidInternalState { .. }, UnsupportedConfiguration => ck::CipLabelerError::UnsupportedConfiguration { .. },
    InvalidTopology => ck::CipLabelerError::InvalidTopology(_), AtomProperty => ck::CipLabelerError::AtomProperty { .. }, BondValue => ck::CipLabelerError::BondValue { .. },
    MoleculeProperty => ck::CipLabelerError::MoleculeProperty(_), Atropisomer => ck::CipLabelerError::Atropisomer(_), DoubleStereo => ck::CipLabelerError::DoubleStereo(_),
    RingFinding => ck::CipLabelerError::RingFinding(_), Valence => ck::CipLabelerError::Valence(_), Kekulize => ck::CipLabelerError::Kekulize(_), PeriodicTable => ck::CipLabelerError::PeriodicTable(_),
}, {
    #[wasm_bindgen(getter)] pub fn atom(&self) -> Option<usize> { match *self.inner { ck::CipLabelerError::UnsignedPropertyRead { atom, .. } | ck::CipLabelerError::InvalidPropertyKind { atom, .. } | ck::CipLabelerError::BondNotIncident { atom, .. } | ck::CipLabelerError::UnsupportedConfiguration { atom, .. } | ck::CipLabelerError::AtomProperty { atom, .. } => Some(atom), _ => None } }
    #[wasm_bindgen(getter)] pub fn bond(&self) -> Option<usize> { match *self.inner { ck::CipLabelerError::BondNotIncident { bond, .. } | ck::CipLabelerError::BondValue { bond, .. } => Some(bond), _ => None } }
    #[wasm_bindgen(getter)] pub fn property(&self) -> Option<String> { match *self.inner { ck::CipLabelerError::UnsignedPropertyRead { property, .. } | ck::CipLabelerError::InvalidPropertyKind { property, .. } => Some(property.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn kind(&self) -> Option<CipPayloadKind> { match *self.inner { ck::CipLabelerError::InvalidPropertyKind { kind, .. } => Some(JsValue::from(PropertyValueKind::from(kind)).unchecked_into()), ck::CipLabelerError::SourceIndexWidthExceeded { kind, .. } => Some(JsValue::from_str(kind).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn index(&self) -> Option<usize> { match *self.inner { ck::CipLabelerError::AtomIndexOutOfRange { index, .. } | ck::CipLabelerError::BondIndexOutOfRange { index, .. } | ck::CipLabelerError::SourceIndexWidthExceeded { index, .. } => Some(index), _ => None } }
    #[wasm_bindgen(getter, js_name = atomCount)] pub fn atom_count(&self) -> Option<usize> { match *self.inner { ck::CipLabelerError::AtomIndexOutOfRange { atom_count, .. } => Some(atom_count), _ => None } }
    #[wasm_bindgen(getter, js_name = bondCount)] pub fn bond_count(&self) -> Option<usize> { match *self.inner { ck::CipLabelerError::BondIndexOutOfRange { bond_count, .. } => Some(bond_count), _ => None } }
    #[wasm_bindgen(getter)] pub fn order(&self) -> Option<BondOrder> { match *self.inner { ck::CipLabelerError::NonIntegerBondOrder { order } => Some(order.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn edge(&self) -> Option<usize> { match *self.inner { ck::CipLabelerError::EdgeEndpointMismatch { edge, .. } => Some(edge), _ => None } }
    #[wasm_bindgen(getter)] pub fn node(&self) -> Option<usize> { match *self.inner { ck::CipLabelerError::EdgeEndpointMismatch { node, .. } => Some(node), _ => None } }
    #[wasm_bindgen(getter)] pub fn limit(&self) -> Option<usize> { match *self.inner { ck::CipLabelerError::TooManyNodes { limit } => Some(limit), _ => None } }
    #[wasm_bindgen(getter)] pub fn detail(&self) -> Option<String> { match *self.inner { ck::CipLabelerError::InvalidInternalState { detail } => Some(detail.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn tag(&self) -> Option<ChiralTag> { match *self.inner { ck::CipLabelerError::UnsupportedConfiguration { tag, .. } => Some(tag.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn atropisomer(&self) -> Option<CipAtropisomer> { match self.inner.as_ref() { ck::CipLabelerError::Atropisomer(value) => Some(atropisomer_errors::AtropisomerError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter, js_name = doubleStereo)] pub fn double_stereo(&self) -> Option<CipDoubleStereo> { match self.inner.as_ref() { ck::CipLabelerError::DoubleStereo(value) => Some(double_stereo_errors::DoubleBondStereoError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter, js_name = ringFinding)] pub fn ring_finding(&self) -> Option<CipRings> { match self.inner.as_ref() { ck::CipLabelerError::RingFinding(value) => Some(ring_errors::RingFindingError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn valence(&self) -> Option<CipValence> { match self.inner.as_ref() { ck::CipLabelerError::Valence(value) => Some(valence_errors::ValenceError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn kekulize(&self) -> Option<CipKekulize> { match self.inner.as_ref() { ck::CipLabelerError::Kekulize(value) => Some(kekulize_errors::KekulizeError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter, js_name = periodicTable)] pub fn periodic_table(&self) -> Option<CipPeriodicTable> { match self.inner.as_ref() { ck::CipLabelerError::PeriodicTable(value) => Some(PeriodicTableError::from_rust(*value).unchecked_into()), _ => None } }
}, source = CipLabelerCause, cause = cause);
