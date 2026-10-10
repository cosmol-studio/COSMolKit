//! Concrete stereo error categories, contexts and cause chains.
use crate::alignment_values::{set, source_error};
use cosmolkit_wasm::rust as ck;
use std::error::Error as RustError;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct StereoError {
    inner: ck::StereoError,
}
#[wasm_bindgen]
impl StereoError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "stereo".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        stereo_error_kind(&self.inner).into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.to_string()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> Result<JsValue, JsValue> {
        self.inner.source().map_or(Ok(JsValue::NULL), source_error)
    }
}
fn stereo_error_kind(source: &ck::StereoError) -> &'static str {
    use ck::StereoError as E;
    match source {
        E::StereoGroup(..) => "StereoGroup",
        E::InvalidTopology(..) => "InvalidTopology",
        E::DuplicateConformerId { .. } => "DuplicateConformerId",
        E::ConformerNotFound { .. } => "ConformerNotFound",
        E::ConformerAtomCountMismatch { .. } => "ConformerAtomCountMismatch",
        E::InvalidValence { .. } => "InvalidValence",
        E::InvalidValenceValue { .. } => "InvalidValenceValue",
        E::ZeroLengthVector { .. } => "ZeroLengthVector",
        E::StereoOrder(..) => "StereoOrder",
        E::AtomProperty(..) => "AtomProperty",
        E::HydrogenCount(..) => "HydrogenCount",
    }
}
pub(crate) fn stereo_error(source: &ck::StereoError) -> Result<JsValue, JsValue> {
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("StereoError");
    let error: JsValue = error.into();
    set(&error, "domain", "stereo".into())?;
    set(&error, "kind", stereo_error_kind(source).into())?;
    use ck::StereoError as E;
    match source {
        E::DuplicateConformerId { id, .. } => {
            set(&error, "id", (*id as f64).into())?;
        }
        E::ConformerNotFound { requested, .. } => {
            set(&error, "requested", (*requested as f64).into())?;
        }
        E::ConformerAtomCountMismatch {
            conformer,
            rows,
            atom_count,
            ..
        } => {
            set(&error, "conformer", (*conformer as f64).into())?;
            set(&error, "rows", (*rows as f64).into())?;
            set(&error, "atomCount", (*atom_count as f64).into())?;
        }
        E::InvalidValence {
            field,
            actual,
            atom_count,
            ..
        } => {
            set(&error, "field", (*field).into())?;
            set(&error, "actual", (*actual as f64).into())?;
            set(&error, "atomCount", (*atom_count as f64).into())?;
        }
        E::InvalidValenceValue {
            field, atom, value, ..
        } => {
            set(&error, "field", (*field).into())?;
            set(&error, "atom", (atom.index() as f64).into())?;
            set(&error, "value", (*value as f64).into())?;
        }
        E::ZeroLengthVector {
            center, neighbor, ..
        } => {
            set(&error, "center", (center.index() as f64).into())?;
            set(&error, "neighbor", (neighbor.index() as f64).into())?;
        }
        _ => {}
    }
    if let Some(cause) = source.source() {
        set(&error, "cause", source_error(cause)?)?;
    }
    set(
        &error,
        "detail",
        StereoError {
            inner: source.clone(),
        }
        .into(),
    )?;
    Ok(error)
}
#[wasm_bindgen]
pub struct PotentialStereoError {
    inner: ck::PotentialStereoError,
}
#[wasm_bindgen]
impl PotentialStereoError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "stereo".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        potential_error_kind(&self.inner).into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.to_string()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> Result<JsValue, JsValue> {
        self.inner.source().map_or(Ok(JsValue::NULL), source_error)
    }
}
fn potential_error_kind(source: &ck::PotentialStereoError) -> &'static str {
    use ck::PotentialStereoError as E;
    match source {
        E::InvalidPropertyKind { .. } => "InvalidPropertyKind",
        E::InvalidRingStereoReference { .. } => "InvalidRingStereoReference",
        E::EmptyRingStereoReferences { .. } => "EmptyRingStereoReferences",
        E::Ring(..) => "Ring",
        E::InvalidTopology(..) => "InvalidTopology",
        E::InvalidValence { .. } => "InvalidValence",
        E::InvalidValenceValue { .. } => "InvalidValenceValue",
        E::AtomOutOfRange { .. } => "AtomOutOfRange",
        E::InvalidRingInfo { .. } => "InvalidRingInfo",
        E::InvalidAtomDegree { .. } => "InvalidAtomDegree",
        E::InvalidBondDegree { .. } => "InvalidBondDegree",
        E::InvalidStereoReferences { .. } => "InvalidStereoReferences",
        E::UnsupportedBondOrder { .. } => "UnsupportedBondOrder",
        E::InvalidChiralPermutation { .. } => "InvalidChiralPermutation",
        E::AtropisomerDependencyUnavailable { .. } => "AtropisomerDependencyUnavailable",
        E::RefinementDidNotConverge { .. } => "RefinementDidNotConverge",
        E::Valence(..) => "Valence",
        E::StereoOrder(..) => "StereoOrder",
        E::DoubleStereo(..) => "DoubleStereo",
        E::BondValue(..) => "BondValue",
    }
}
pub(crate) fn potential_error(source: &ck::PotentialStereoError) -> Result<JsValue, JsValue> {
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("PotentialStereoError");
    let error: JsValue = error.into();
    set(&error, "domain", "stereo".into())?;
    set(&error, "kind", potential_error_kind(source).into())?;
    use ck::PotentialStereoError as E;
    match source {
        E::InvalidPropertyKind {
            atom,
            property,
            kind,
            ..
        } => {
            set(&error, "atom", (atom.index() as f64).into())?;
            set(&error, "property", (*property).into())?;
            set(
                &error,
                "propertyKind",
                crate::property_values::kind(*kind).into(),
            )?;
        }
        E::InvalidRingStereoReference { atom, value, .. } => {
            set(&error, "atom", (atom.index() as f64).into())?;
            set(&error, "value", (*value as f64).into())?;
        }
        E::EmptyRingStereoReferences { atom, .. } => {
            set(&error, "atom", (atom.index() as f64).into())?;
        }
        E::InvalidValence {
            field,
            actual,
            atom_count,
            ..
        } => {
            set(&error, "field", (*field).into())?;
            set(&error, "actual", (*actual as f64).into())?;
            set(&error, "atomCount", (*atom_count as f64).into())?;
        }
        E::InvalidValenceValue {
            field, atom, value, ..
        } => {
            set(&error, "field", (*field).into())?;
            set(&error, "atom", (atom.index() as f64).into())?;
            set(&error, "value", (*value as f64).into())?;
        }
        E::AtomOutOfRange {
            atom, atom_count, ..
        } => {
            set(&error, "atom", (atom.index() as f64).into())?;
            set(&error, "atomCount", (*atom_count as f64).into())?;
        }
        E::InvalidRingInfo {
            reason,
            row,
            value,
            limit,
            ..
        } => {
            set(&error, "reason", (*reason).into())?;
            set(&error, "row", (*row as f64).into())?;
            set(&error, "value", (*value as f64).into())?;
            set(&error, "limit", (*limit as f64).into())?;
        }
        E::InvalidAtomDegree { atom, degree, .. } => {
            set(&error, "atom", (atom.index() as f64).into())?;
            set(&error, "degree", (*degree as f64).into())?;
        }
        E::InvalidBondDegree {
            bond,
            endpoint,
            degree,
            ..
        } => {
            set(&error, "bond", (bond.index() as f64).into())?;
            set(&error, "endpoint", (*endpoint).into())?;
            set(&error, "degree", (*degree as f64).into())?;
        }
        E::InvalidStereoReferences { bond, reason, .. } => {
            set(&error, "bond", (bond.index() as f64).into())?;
            set(&error, "reason", (*reason).into())?;
        }
        E::UnsupportedBondOrder { bond, order, .. } => {
            set(&error, "bond", (bond.index() as f64).into())?;
            set(&error, "order", (order.rdkit_code() as f64).into())?;
        }
        E::InvalidChiralPermutation {
            atom, permutation, ..
        } => {
            set(&error, "atom", (atom.index() as f64).into())?;
            set(&error, "permutation", (*permutation as f64).into())?;
        }
        E::AtropisomerDependencyUnavailable { bond, .. } => {
            set(&error, "bond", (bond.index() as f64).into())?;
        }
        E::RefinementDidNotConverge { iterations, .. } => {
            set(&error, "iterations", (*iterations as f64).into())?;
        }
        _ => {}
    }
    if let Some(cause) = source.source() {
        set(&error, "cause", source_error(cause)?)?;
    }
    set(
        &error,
        "detail",
        PotentialStereoError {
            inner: source.clone(),
        }
        .into(),
    )?;
    Ok(error)
}
#[wasm_bindgen]
pub struct CipLabelerError {
    inner: ck::CipLabelerError,
}
#[wasm_bindgen]
impl CipLabelerError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "stereo".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        cip_error_kind(&self.inner).into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.to_string()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> Result<JsValue, JsValue> {
        self.inner.source().map_or(Ok(JsValue::NULL), source_error)
    }
}
fn cip_error_kind(source: &ck::CipLabelerError) -> &'static str {
    use ck::CipLabelerError as E;
    match source {
        E::UnsignedPropertyRead { .. } => "UnsignedPropertyRead",
        E::InvalidPropertyKind { .. } => "InvalidPropertyKind",
        E::AtomIndexOutOfRange { .. } => "AtomIndexOutOfRange",
        E::BondIndexOutOfRange { .. } => "BondIndexOutOfRange",
        E::SourceIndexWidthExceeded { .. } => "SourceIndexWidthExceeded",
        E::BondNotIncident { .. } => "BondNotIncident",
        E::NonIntegerBondOrder { .. } => "NonIntegerBondOrder",
        E::EdgeEndpointMismatch { .. } => "EdgeEndpointMismatch",
        E::TooManyNodes { .. } => "TooManyNodes",
        E::MaxIterationsExceeded => "MaxIterationsExceeded",
        E::UnexpectedUpEdgeOrdering => "UnexpectedUpEdgeOrdering",
        E::NoSequenceRuleProvided => "NoSequenceRuleProvided",
        E::EmptyPairListReference => "EmptyPairListReference",
        E::DescriptorListLengthMismatch => "DescriptorListLengthMismatch",
        E::InvalidStereoDescriptor => "InvalidStereoDescriptor",
        E::SubstituentsShouldBeTopologicallyEquivalent => {
            "SubstituentsShouldBeTopologicallyEquivalent"
        }
        E::SomethingUnexpected => "SomethingUnexpected",
        E::Rule4bInstanceNotInRuleSet => "Rule4bInstanceNotInRuleSet",
        E::Rule5NewInstanceNotInRuleSet => "Rule5NewInstanceNotInRuleSet",
        E::InvalidInternalState { .. } => "InvalidInternalState",
        E::ParityVectorsMustHaveSize4 => "ParityVectorsMustHaveSize4",
        E::EmptyConfigurationFoci => "EmptyConfigurationFoci",
        E::BadTetrahedralConfig => "BadTetrahedralConfig",
        E::TetrahedralConfigurationMustHave4Carriers => "TetrahedralConfigurationMustHave4Carriers",
        E::DescriptorNotSupportedForAtoms => "DescriptorNotSupportedForAtoms",
        E::InvalidAtomDescriptor => "InvalidAtomDescriptor",
        E::CarrierMismatch => "CarrierMismatch",
        E::BadSp2BondFoci => "BadSp2BondFoci",
        E::BadSp2BondConfig => "BadSp2BondConfig",
        E::IncorrectNumberOfStereoAtoms => "IncorrectNumberOfStereoAtoms",
        E::DescriptorNotSupportedForDoubleBonds => "DescriptorNotSupportedForDoubleBonds",
        E::InvalidBondDescriptor => "InvalidBondDescriptor",
        E::BadAtropisomerBondFoci => "BadAtropisomerBondFoci",
        E::BadAtropisomerBondConfig => "BadAtropisomerBondConfig",
        E::DescriptorNotSupportedForAtropisomerBonds => "DescriptorNotSupportedForAtropisomerBonds",
        E::UnsupportedConfiguration { .. } => "UnsupportedConfiguration",
        E::InvalidTopology(..) => "InvalidTopology",
        E::AtomProperty { .. } => "AtomProperty",
        E::BondValue { .. } => "BondValue",
        E::MoleculeProperty(..) => "MoleculeProperty",
        E::Atropisomer(..) => "Atropisomer",
        E::DoubleStereo(..) => "DoubleStereo",
        E::RingFinding(..) => "RingFinding",
        E::Valence(..) => "Valence",
        E::Kekulize(..) => "Kekulize",
        E::PeriodicTable(..) => "PeriodicTable",
    }
}
pub(crate) fn cip_error(source: &ck::CipLabelerError) -> Result<JsValue, JsValue> {
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("CipLabelerError");
    let error: JsValue = error.into();
    set(&error, "domain", "stereo".into())?;
    set(&error, "kind", cip_error_kind(source).into())?;
    use ck::CipLabelerError as E;
    match source {
        E::UnsignedPropertyRead { atom, property, .. } => {
            set(&error, "atom", (*atom as f64).into())?;
            set(&error, "property", (*property).into())?;
        }
        E::InvalidPropertyKind {
            atom,
            property,
            kind,
            ..
        } => {
            set(&error, "atom", (*atom as f64).into())?;
            set(&error, "property", (*property).into())?;
            set(
                &error,
                "propertyKind",
                crate::property_values::kind(*kind).into(),
            )?;
        }
        E::AtomIndexOutOfRange {
            index, atom_count, ..
        } => {
            set(&error, "index", (*index as f64).into())?;
            set(&error, "atomCount", (*atom_count as f64).into())?;
        }
        E::BondIndexOutOfRange {
            index, bond_count, ..
        } => {
            set(&error, "index", (*index as f64).into())?;
            set(&error, "bondCount", (*bond_count as f64).into())?;
        }
        E::SourceIndexWidthExceeded { kind, index, .. } => {
            set(&error, "propertyKind", (*kind).into())?;
            set(&error, "index", (*index as f64).into())?;
        }
        E::BondNotIncident { bond, atom, .. } => {
            set(&error, "bond", (*bond as f64).into())?;
            set(&error, "atom", (*atom as f64).into())?;
        }
        E::NonIntegerBondOrder { order, .. } => {
            set(&error, "order", (order.rdkit_code() as f64).into())?;
        }
        E::EdgeEndpointMismatch { edge, node, .. } => {
            set(&error, "edge", (*edge as f64).into())?;
            set(&error, "node", (*node as f64).into())?;
        }
        E::TooManyNodes { limit, .. } => {
            set(&error, "limit", (*limit as f64).into())?;
        }
        E::InvalidInternalState { detail, .. } => {
            set(&error, "detail", (*detail).into())?;
        }
        E::UnsupportedConfiguration { atom, tag, .. } => {
            set(&error, "atom", (*atom as f64).into())?;
            set(&error, "tag", (tag.rdkit_code() as f64).into())?;
        }
        E::AtomProperty { atom, .. } => {
            set(&error, "atom", (*atom as f64).into())?;
        }
        E::BondValue { bond, .. } => {
            set(&error, "bond", (*bond as f64).into())?;
        }
        _ => {}
    }
    if let Some(cause) = source.source() {
        set(&error, "cause", source_error(cause)?)?;
    }
    set(
        &error,
        "detail",
        CipLabelerError {
            inner: source.clone(),
        }
        .into(),
    )?;
    Ok(error)
}
