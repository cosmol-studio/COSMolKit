//! Operation failures retain the canonical outer type and its actual source.
use crate::{
    errors::{enum_error, exception},
    model_errors::*,
    ring_errors::RingFindingError,
    topology_errors::TopologyEditError,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};

#[wasm_bindgen]
pub struct FunctionStatus(ck::FunctionStatus);
#[wasm_bindgen]
impl FunctionStatus {
    #[wasm_bindgen(getter,js_name=Native)]
    pub fn native() -> Self {
        Self(ck::FunctionStatus::Native)
    }
    #[wasm_bindgen(getter,js_name=Experimental)]
    pub fn experimental() -> Self {
        Self(ck::FunctionStatus::Experimental)
    }
    #[wasm_bindgen(getter)]
    pub fn reference(&self) -> Option<String> {
        match self.0 {
            ck::FunctionStatus::Parity { reference }
            | ck::FunctionStatus::ParityWithDifferences { reference, .. } => Some(reference.into()),
            ck::FunctionStatus::Native | ck::FunctionStatus::Experimental => None,
        }
    }
    #[wasm_bindgen(getter)]
    pub fn explanation(&self) -> Option<String> {
        match self.0 {
            ck::FunctionStatus::ParityWithDifferences { explanation, .. } => {
                Some(explanation.into())
            }
            _ => None,
        }
    }
    #[wasm_bindgen(js_name = toString)]
    pub fn display(&self) -> String {
        format!("{:?}", self.0)
    }
}
#[wasm_bindgen]
pub struct MoleculeOpSpec(pub(crate) &'static ck::MoleculeOpSpec);
#[wasm_bindgen]
impl MoleculeOpSpec {
    #[wasm_bindgen(getter)]
    pub fn method(&self) -> String {
        self.0.method.into()
    }
    #[wasm_bindgen(getter, js_name = implFn)]
    pub fn impl_fn(&self) -> String {
        self.0.impl_fn.into()
    }
    #[wasm_bindgen(getter, js_name = resultType)]
    pub fn result_type(&self) -> String {
        self.0.result_type.into()
    }
    #[wasm_bindgen(getter)]
    pub fn status(&self) -> FunctionStatus {
        FunctionStatus(self.0.status)
    }
    #[wasm_bindgen(getter)]
    pub fn parity(&self) -> crate::enums::ParityPolicy {
        self.0.parity.into()
    }
    #[wasm_bindgen(getter, js_name = ioRoundtrip)]
    pub fn io_roundtrip(&self) -> bool {
        self.0.io_roundtrip
    }
}
#[wasm_bindgen]
pub struct UnsupportedFeatureError(ck::UnsupportedFeatureError);
#[wasm_bindgen]
impl UnsupportedFeatureError {
    #[wasm_bindgen(getter)]
    pub fn feature(&self) -> String {
        self.0.feature.into()
    }
    #[wasm_bindgen(getter)]
    pub fn reason(&self) -> String {
        self.0.reason.into()
    }
}
impl UnsupportedFeatureError {
    fn from_rust(value: ck::UnsupportedFeatureError) -> JsValue {
        Self(value).into()
    }
}

macro_rules! context_value {
    ($name:ident, $getter:ident, $js:literal, $ty:ty, $read:expr, $display:expr) => {
        #[wasm_bindgen]
        pub struct $name(crate::storage::ReadonlyValue<ck::OperationError>);
        #[wasm_bindgen]
        impl $name {
            #[wasm_bindgen(getter, js_name = $js)]
            pub fn $getter(&self) -> $ty {
                ($read)(self.0.as_ref())
            }
            #[wasm_bindgen(js_name = toString)]
            pub fn display(&self) -> String {
                ($display)(self.0.as_ref())
            }
        }
    };
}
context_value!(
    OperationDerivedState,
    bits,
    "bits",
    u16,
    |error: &ck::OperationError| match error {
        ck::OperationError::DerivedEffectContract { states, .. }
        | ck::OperationError::DerivedPreservationMismatch { states, .. } => states.bits() as u16,
        _ => unreachable!("constructor selects derived-state variants"),
    },
    |error: &ck::OperationError| match error {
        ck::OperationError::DerivedEffectContract { states, .. }
        | ck::OperationError::DerivedPreservationMismatch { states, .. } => format!("{states:?}"),
        _ => unreachable!("constructor selects derived-state variants"),
    }
);
context_value!(
    OperationSemanticPreconditions,
    bits,
    "bits",
    u16,
    |error: &ck::OperationError| match error {
        ck::OperationError::SemanticPreconditionContract { missing, .. } => missing.bits() as u16,
        _ => unreachable!("constructor selects semantic-precondition variant"),
    },
    |error: &ck::OperationError| match error {
        ck::OperationError::SemanticPreconditionContract { missing, .. } => format!("{missing:?}"),
        _ => unreachable!("constructor selects semantic-precondition variant"),
    }
);
context_value!(
    OperationCipStatePolicy,
    discriminant,
    "discriminant",
    i32,
    |error: &ck::OperationError| match error {
        ck::OperationError::CipStateContract { policy, .. } => *policy as i32,
        _ => unreachable!("constructor selects CIP-policy variant"),
    },
    |error: &ck::OperationError| match error {
        ck::OperationError::CipStateContract { policy, .. } => format!("{policy:?}"),
        _ => unreachable!("constructor selects CIP-policy variant"),
    }
);
context_value!(
    OperationMappingRequirement,
    discriminant,
    "discriminant",
    i32,
    |error: &ck::OperationError| match error {
        ck::OperationError::MappingContract { requirement, .. } => *requirement as i32,
        _ => unreachable!("constructor selects mapping variant"),
    },
    |error: &ck::OperationError| match error {
        ck::OperationError::MappingContract { requirement, .. } => format!("{requirement:?}"),
        _ => unreachable!("constructor selects mapping variant"),
    }
);

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "TopologyValidationError | TopologyEditError | CoordinateValidationError | MoleculePropertyError | AtomPropertyError | BondValueError | RingFindingError | KekulizeError | AromaticityError | ScaffoldError | ValenceError | PotentialStereoError | StereoError | CipLabelerError | SanitizeError | HydrogenError | AtomCodeError | Coordinate2DError | AlignmentError | TransformError | CoordinateInputError | MoleculeFragmentsError | RadicalError | UffOptimizationError | MmffOptimizationError | ReactionRunError | ReactionApplyError | TautomerRunError | ConformerRunError | EnumerationRunError"
    )]
    #[derive(Clone)]
    pub type OperationCause;
    #[wasm_bindgen(typescript_type = "string | MoleculeOpSpec")]
    pub type OperationName;
    #[wasm_bindgen(typescript_type = "number | string | MoleculeOpOutput")]
    pub type OperationScalar;
    #[wasm_bindgen(typescript_type = "UnsupportedFeatureError | MappingValidationError")]
    pub type OperationSource;
    #[wasm_bindgen(typescript_type = "TopologyValidationError")]
    pub type InvalidTopologyPayload;
    #[wasm_bindgen(typescript_type = "TopologyEditError")]
    pub type InvalidTopologyEditPayload;
    #[wasm_bindgen(typescript_type = "CoordinateValidationError")]
    pub type CoordinatesPayload;
    #[wasm_bindgen(typescript_type = "MoleculePropertyError")]
    pub type PropertyPayload;
    #[wasm_bindgen(typescript_type = "AtomPropertyError")]
    pub type AtomPropertyPayload;
    #[wasm_bindgen(typescript_type = "BondValueError")]
    pub type BondPropertyPayload;
    #[wasm_bindgen(typescript_type = "RingFindingError")]
    pub type RingsPayload;
    #[wasm_bindgen(typescript_type = "KekulizeError")]
    pub type KekulizePayload;
    #[wasm_bindgen(typescript_type = "AromaticityError")]
    pub type AromaticityPayload;
    #[wasm_bindgen(typescript_type = "ScaffoldError")]
    pub type ScaffoldPayload;
    #[wasm_bindgen(typescript_type = "ValenceError")]
    pub type ValencePayload;
    #[wasm_bindgen(typescript_type = "PotentialStereoError")]
    pub type PotentialStereoPayload;
    #[wasm_bindgen(typescript_type = "StereoError")]
    pub type StereoPayload;
    #[wasm_bindgen(typescript_type = "CipLabelerError")]
    pub type CipLabelerPayload;
    #[wasm_bindgen(typescript_type = "SanitizeError")]
    pub type SanitizePayload;
    #[wasm_bindgen(typescript_type = "HydrogenError")]
    pub type HydrogenPayload;
    #[wasm_bindgen(typescript_type = "AtomCodeError")]
    pub type AtomCodePayload;
    #[wasm_bindgen(typescript_type = "Coordinate2DError")]
    pub type Coordinate2DPayload;
    #[wasm_bindgen(typescript_type = "AlignmentError")]
    pub type AlignmentPayload;
    #[wasm_bindgen(
        typescript_type = "TransformError | CoordinateInputError | MoleculeFragmentsError"
    )]
    pub type TransformPayload;
    #[wasm_bindgen(typescript_type = "RadicalError")]
    pub type RadicalPayload;
    #[wasm_bindgen(typescript_type = "ConformerRunError")]
    pub type ConformerPayload;
    #[wasm_bindgen(typescript_type = "EnumerationRunError")]
    pub type EnumerationPayload;
    #[wasm_bindgen(typescript_type = "TautomerRunError")]
    pub type TautomerPayload;
}
fn cause(inner: &ck::OperationError) -> Option<OperationCause> {
    match inner {
        #[cfg(feature = "cap-stereoisomers")]
        ck::OperationError::Enumeration(v) => {
            Some(crate::stereoisomers::EnumerationRunError::from_rust(v.clone()).unchecked_into())
        }
        #[cfg(feature = "cap-conformer")]
        ck::OperationError::Conformer(v) => {
            Some(crate::conformer_errors::ConformerRunError::from_rust(v.clone()).unchecked_into())
        }
        #[cfg(feature = "cap-tautomer")]
        ck::OperationError::Tautomer(v) => {
            Some(crate::tautomer_errors::TautomerRunError::from_rust(v.clone()).unchecked_into())
        }
        #[cfg(feature = "cap-reaction")]
        ck::OperationError::ReactionRun(v) => Some(
            crate::reaction_execution_errors::ReactionRunError::from_rust(v.clone())
                .unchecked_into(),
        ),
        #[cfg(feature = "cap-reaction")]
        ck::OperationError::ReactionApply(v) => Some(
            crate::reaction_execution_errors::ReactionApplyError::from_rust(v.clone())
                .unchecked_into(),
        ),
        #[cfg(feature = "cap-forcefields")]
        ck::OperationError::UffOptimization(v) => Some(
            crate::forcefield_execution_errors::UffOptimizationError::from_rust(v.clone())
                .unchecked_into(),
        ),
        #[cfg(feature = "cap-forcefields")]
        ck::OperationError::MmffOptimization(v) => Some(
            crate::forcefield_execution_errors::MmffOptimizationError::from_rust(v.clone())
                .unchecked_into(),
        ),
        #[cfg(feature = "cap-radicals")]
        ck::OperationError::Radical(v) => {
            Some(crate::sanitize_leaf_errors::RadicalError::from_rust(v.clone()).unchecked_into())
        }
        #[cfg(feature = "cap-transforms")]
        ck::OperationError::Transform(v) => {
            Some(crate::transform_errors::TransformError::from_rust(v.clone()).unchecked_into())
        }
        #[cfg(feature = "cap-transforms")]
        ck::OperationError::CoordinateInput(v) => Some(
            crate::transform_errors::CoordinateInputError::from_rust(v.clone()).unchecked_into(),
        ),
        #[cfg(feature = "cap-transforms")]
        ck::OperationError::Fragments(v) => Some(
            crate::fragment_errors::MoleculeFragmentsError::from_rust(v.clone()).unchecked_into(),
        ),
        #[cfg(feature = "cap-alignment")]
        ck::OperationError::Alignment(v) => {
            Some(crate::alignment_errors::AlignmentError::from_rust(v.clone()).unchecked_into())
        }
        #[cfg(feature = "cap-fingerprints")]
        ck::OperationError::AtomCode(v) => {
            Some(crate::fingerprint_errors::AtomCodeError::from_rust(v.clone()).unchecked_into())
        }
        #[cfg(feature = "cap-depict")]
        ck::OperationError::Coordinate2D(v) => {
            Some(crate::depict_errors::Coordinate2DError::from_rust(v.clone()).unchecked_into())
        }
        #[cfg(feature = "cap-hashing")]
        ck::OperationError::Scaffold(value) => {
            Some(crate::scaffold_errors::ScaffoldError::from_rust(value.clone()).unchecked_into())
        }
        #[cfg(any(feature = "cap-valence", feature = "cap-stereo"))]
        ck::OperationError::Valence(value) => {
            Some(crate::valence_errors::ValenceError::from_rust(value.clone()).unchecked_into())
        }
        #[cfg(feature = "cap-stereo")]
        ck::OperationError::PotentialStereo(value) => Some(
            crate::potential_stereo_errors::PotentialStereoError::from_rust(value.clone())
                .unchecked_into(),
        ),
        #[cfg(feature = "cap-stereo")]
        ck::OperationError::Stereo(value) => {
            Some(crate::stereo_errors::StereoError::from_rust(value.clone()).unchecked_into())
        }
        #[cfg(feature = "cap-stereo")]
        ck::OperationError::CipLabeler(value) => Some(
            crate::cip_labeler_errors::CipLabelerError::from_rust(value.clone()).unchecked_into(),
        ),
        #[cfg(feature = "cap-sanitize")]
        ck::OperationError::Sanitize(value) => {
            Some(crate::sanitize_errors::SanitizeError::from_rust(value.clone()).unchecked_into())
        }
        #[cfg(feature = "cap-hydrogens")]
        ck::OperationError::Hydrogen(value) => {
            Some(crate::hydrogen_errors::HydrogenError::from_rust(value.clone()).unchecked_into())
        }
        #[cfg(feature = "cap-aromaticity")]
        ck::OperationError::Aromaticity(value) => Some(
            crate::aromaticity_errors::AromaticityError::from_rust(value.clone()).unchecked_into(),
        ),
        #[cfg(feature = "cap-kekulize")]
        ck::OperationError::Kekulize(value) => {
            Some(crate::kekulize_errors::KekulizeError::from_rust(value.clone()).unchecked_into())
        }
        ck::OperationError::InvalidTopology(
            value @ ck::TopologyValidationError::StereoGroup(_),
        ) => Some(TopologyValidationError::from_rust(value.clone()).unchecked_into()),
        ck::OperationError::InvalidTopologyEdit(value @ ck::TopologyEditError::StereoGroup(_)) => {
            Some(TopologyEditError::from_rust(value.clone()).unchecked_into())
        }
        ck::OperationError::InvalidCoordinates(value) => {
            Some(CoordinateValidationError::from_rust(value.clone()).unchecked_into())
        }
        ck::OperationError::InvalidProperty(value) => {
            Some(MoleculePropertyError::from_rust(value.clone()).unchecked_into())
        }
        ck::OperationError::AtomProperty(value) => {
            Some(AtomPropertyError::from_rust(value.clone()).unchecked_into())
        }
        ck::OperationError::BondProperty(value) => {
            Some(BondValueError::from_rust(value.clone()).unchecked_into())
        }
        #[cfg(any(
            feature = "cap-inchi",
            feature = "cap-forcefields",
            feature = "cap-rings",
            feature = "cap-stereo",
            feature = "cap-aromaticity",
            feature = "cap-fingerprints",
            feature = "cap-serialization",
            feature = "cap-io"
        ))]
        ck::OperationError::Rings(value) => {
            Some(RingFindingError::from_rust(value.clone()).unchecked_into())
        }
        ck::OperationError::AtomPropertyIndex { .. }
        | ck::OperationError::ReservedAtomPropertyKey { .. }
        | ck::OperationError::EmptyFragments
        | ck::OperationError::UnsupportedFeature { .. }
        | ck::OperationError::Unsupported { .. }
        | ck::OperationError::OutputMismatch { .. }
        | ck::OperationError::AccessDenied { .. }
        | ck::OperationError::BlockCheckedOut { .. }
        | ck::OperationError::BlockNotCheckedOut { .. }
        | ck::OperationError::IncompleteCommit { .. }
        | ck::OperationError::MappingContract { .. }
        | ck::OperationError::InvalidTopologyMapping { .. }
        | ck::OperationError::AutoRemapContract { .. }
        | ck::OperationError::OperationContract { .. }
        | ck::OperationError::SemanticPreconditionContract { .. }
        | ck::OperationError::CoordinateAppendRequiresValues { .. }
        | ck::OperationError::DerivedEffectContract { .. }
        | ck::OperationError::DerivedPreservationMismatch { .. }
        | ck::OperationError::CipStateContract { .. }
        | ck::OperationError::InvalidTopology(_)
        | ck::OperationError::InvalidTopologyEdit(_)
        | ck::OperationError::InvalidPropertyList { .. }
        | ck::OperationError::InvalidDerivedCache { .. }
        | ck::OperationError::InvalidAlgorithmResult { .. }
        | ck::OperationError::InvalidReconstructionOrigin { .. }
        | ck::OperationError::Algorithm { .. } => None,
    }
}
enum_error!(OperationError, OperationErrorKind, ck::OperationError, {
    #[cfg(feature = "cap-reaction")]
    ReactionRun=>ck::OperationError::ReactionRun(_),
    #[cfg(feature = "cap-tautomer")]
    Tautomer=>ck::OperationError::Tautomer(_),
    #[cfg(feature = "cap-conformer")]
    Conformer=>ck::OperationError::Conformer(_),
    #[cfg(feature = "cap-stereoisomers")]
    Enumeration=>ck::OperationError::Enumeration(_),
    #[cfg(feature = "cap-reaction")]
    ReactionApply=>ck::OperationError::ReactionApply(_),
    #[cfg(feature = "cap-forcefields")]
    UffOptimization=>ck::OperationError::UffOptimization(_),
    #[cfg(feature = "cap-forcefields")]
    MmffOptimization=>ck::OperationError::MmffOptimization(_),
    #[cfg(feature = "cap-radicals")]
    Radical=>ck::OperationError::Radical(_),
    #[cfg(feature = "cap-transforms")]
    Transform=>ck::OperationError::Transform(_),#[cfg(feature = "cap-transforms")]
    CoordinateInput=>ck::OperationError::CoordinateInput(_),#[cfg(feature = "cap-transforms")]
    Fragments=>ck::OperationError::Fragments(_),
    #[cfg(feature = "cap-alignment")]
    Alignment=>ck::OperationError::Alignment(_),
    #[cfg(feature = "cap-fingerprints")]
    AtomCode=>ck::OperationError::AtomCode(_),#[cfg(feature = "cap-depict")]
    Coordinate2D=>ck::OperationError::Coordinate2D(_),
    #[cfg(feature = "cap-aromaticity")]
    Aromaticity => ck::OperationError::Aromaticity(_),
    #[cfg(feature = "cap-hashing")]
    Scaffold => ck::OperationError::Scaffold(_),
    #[cfg(any(feature = "cap-valence", feature = "cap-stereo"))]
    Valence => ck::OperationError::Valence(_),
    #[cfg(feature = "cap-stereo")]
    PotentialStereo => ck::OperationError::PotentialStereo(_),
    #[cfg(feature = "cap-stereo")]
    Stereo => ck::OperationError::Stereo(_),
    #[cfg(feature = "cap-stereo")]
    CipLabeler => ck::OperationError::CipLabeler(_),
    #[cfg(feature = "cap-sanitize")]
    Sanitize => ck::OperationError::Sanitize(_),
    #[cfg(feature = "cap-hydrogens")]
    Hydrogen => ck::OperationError::Hydrogen(_),
    #[cfg(feature = "cap-kekulize")]
    Kekulize => ck::OperationError::Kekulize(_),
    AtomPropertyIndex => ck::OperationError::AtomPropertyIndex { .. }, ReservedAtomPropertyKey => ck::OperationError::ReservedAtomPropertyKey { .. },
    EmptyFragments => ck::OperationError::EmptyFragments, UnsupportedFeature => ck::OperationError::UnsupportedFeature { .. }, Unsupported => ck::OperationError::Unsupported { .. },
    OutputMismatch => ck::OperationError::OutputMismatch { .. }, AccessDenied => ck::OperationError::AccessDenied { .. }, BlockCheckedOut => ck::OperationError::BlockCheckedOut { .. },
    BlockNotCheckedOut => ck::OperationError::BlockNotCheckedOut { .. }, IncompleteCommit => ck::OperationError::IncompleteCommit { .. }, MappingContract => ck::OperationError::MappingContract { .. },
    InvalidTopologyMapping => ck::OperationError::InvalidTopologyMapping { .. }, AutoRemapContract => ck::OperationError::AutoRemapContract { .. }, OperationContract => ck::OperationError::OperationContract { .. },
    SemanticPreconditionContract => ck::OperationError::SemanticPreconditionContract { .. }, CoordinateAppendRequiresValues => ck::OperationError::CoordinateAppendRequiresValues { .. },
    DerivedEffectContract => ck::OperationError::DerivedEffectContract { .. }, DerivedPreservationMismatch => ck::OperationError::DerivedPreservationMismatch { .. }, CipStateContract => ck::OperationError::CipStateContract { .. },
    InvalidTopology => ck::OperationError::InvalidTopology(_), InvalidTopologyEdit => ck::OperationError::InvalidTopologyEdit(_), InvalidCoordinates => ck::OperationError::InvalidCoordinates(_),
    InvalidProperty => ck::OperationError::InvalidProperty(_), AtomProperty => ck::OperationError::AtomProperty(_), BondProperty => ck::OperationError::BondProperty(_), #[cfg(any(
        feature = "cap-inchi",
        feature = "cap-forcefields",
        feature = "cap-rings",
        feature = "cap-stereo",
        feature = "cap-aromaticity",
        feature = "cap-fingerprints",
        feature = "cap-serialization",
        feature = "cap-io"
    ))]
    Rings => ck::OperationError::Rings(_),
    InvalidPropertyList => ck::OperationError::InvalidPropertyList { .. }, InvalidDerivedCache => ck::OperationError::InvalidDerivedCache { .. }, InvalidAlgorithmResult => ck::OperationError::InvalidAlgorithmResult { .. },
    InvalidReconstructionOrigin => ck::OperationError::InvalidReconstructionOrigin { .. }, Algorithm => ck::OperationError::Algorithm { .. },
}, {
    #[cfg(feature = "cap-radicals")]
    #[wasm_bindgen(getter)]pub fn radical(&self)->Option<RadicalPayload>{match self.inner.as_ref(){ck::OperationError::Radical(v)=>Some(crate::sanitize_leaf_errors::RadicalError::from_rust(v.clone()).unchecked_into()),_=>None}}
    #[cfg(feature = "cap-transforms")]
    #[wasm_bindgen(getter)] pub fn transform(&self)->Option<TransformPayload>{match self.inner.as_ref(){ck::OperationError::Transform(v)=>Some(crate::transform_errors::TransformError::from_rust(v.clone()).unchecked_into()),_=>None}}
    #[cfg(feature = "cap-transforms")]
    #[wasm_bindgen(getter,js_name=coordinateInput)] pub fn coordinate_input(&self)->Option<TransformPayload>{match self.inner.as_ref(){ck::OperationError::CoordinateInput(v)=>Some(crate::transform_errors::CoordinateInputError::from_rust(v.clone()).unchecked_into()),_=>None}}
    #[cfg(feature = "cap-transforms")]
    #[wasm_bindgen(getter)] pub fn fragments(&self)->Option<TransformPayload>{match self.inner.as_ref(){ck::OperationError::Fragments(v)=>Some(crate::fragment_errors::MoleculeFragmentsError::from_rust(v.clone()).unchecked_into()),_=>None}}
    #[cfg(feature = "cap-fingerprints")]
    #[wasm_bindgen(getter, js_name = atomCode)] pub fn atom_code(&self) -> Option<AtomCodePayload> { match self.inner.as_ref() { ck::OperationError::AtomCode(value) => Some(crate::fingerprint_errors::AtomCodeError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[cfg(feature = "cap-depict")]
    #[wasm_bindgen(getter, js_name = coordinate2D)] pub fn coordinate_2d(&self) -> Option<Coordinate2DPayload> { match self.inner.as_ref() { ck::OperationError::Coordinate2D(value) => Some(crate::depict_errors::Coordinate2DError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[cfg(feature = "cap-alignment")]
    #[wasm_bindgen(getter)] pub fn alignment(&self) -> Option<AlignmentPayload> { match self.inner.as_ref() { ck::OperationError::Alignment(value) => Some(crate::alignment_errors::AlignmentError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[cfg(feature = "cap-hashing")]
    #[wasm_bindgen(getter, js_name = scaffold)] pub fn scaffold(&self) -> Option<ScaffoldPayload> { match self.inner.as_ref() { ck::OperationError::Scaffold(value) => Some(crate::scaffold_errors::ScaffoldError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[cfg(any(feature = "cap-valence", feature = "cap-stereo"))]
    #[wasm_bindgen(getter, js_name = valence)] pub fn valence(&self) -> Option<ValencePayload> { match self.inner.as_ref() { ck::OperationError::Valence(value) => Some(crate::valence_errors::ValenceError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[cfg(feature = "cap-stereo")]
    #[wasm_bindgen(getter, js_name = potentialStereo)] pub fn potential_stereo(&self) -> Option<PotentialStereoPayload> { match self.inner.as_ref() { ck::OperationError::PotentialStereo(value) => Some(crate::potential_stereo_errors::PotentialStereoError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[cfg(feature = "cap-stereo")]
    #[wasm_bindgen(getter, js_name = stereo)] pub fn stereo(&self) -> Option<StereoPayload> { match self.inner.as_ref() { ck::OperationError::Stereo(value) => Some(crate::stereo_errors::StereoError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[cfg(feature = "cap-stereo")]
    #[wasm_bindgen(getter, js_name = cipLabeler)] pub fn cip_labeler(&self) -> Option<CipLabelerPayload> { match self.inner.as_ref() { ck::OperationError::CipLabeler(value) => Some(crate::cip_labeler_errors::CipLabelerError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[cfg(feature = "cap-sanitize")]
    #[wasm_bindgen(getter, js_name = sanitize)] pub fn sanitize(&self) -> Option<SanitizePayload> { match self.inner.as_ref() { ck::OperationError::Sanitize(value) => Some(crate::sanitize_errors::SanitizeError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[cfg(feature = "cap-hydrogens")]
    #[wasm_bindgen(getter, js_name = hydrogen)] pub fn hydrogen(&self) -> Option<HydrogenPayload> { match self.inner.as_ref() { ck::OperationError::Hydrogen(value) => Some(crate::hydrogen_errors::HydrogenError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[cfg(feature = "cap-aromaticity")]
    #[wasm_bindgen(getter)] pub fn aromaticity(&self) -> Option<AromaticityPayload> { match self.inner.as_ref() { ck::OperationError::Aromaticity(value) => Some(crate::aromaticity_errors::AromaticityError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[cfg(feature = "cap-conformer")]
    #[wasm_bindgen(getter)]pub fn conformer(&self)->Option<ConformerPayload>{match self.inner.as_ref(){ck::OperationError::Conformer(v)=>Some(crate::conformer_errors::ConformerRunError::from_rust(v.clone()).unchecked_into()),_=>None}}
    #[cfg(feature = "cap-stereoisomers")]
    #[wasm_bindgen(getter)]pub fn enumeration(&self)->Option<EnumerationPayload>{match self.inner.as_ref(){ck::OperationError::Enumeration(v)=>Some(crate::stereoisomers::EnumerationRunError::from_rust(v.clone()).unchecked_into()),_=>None}}
    #[cfg(feature = "cap-tautomer")]
    #[wasm_bindgen(getter)]pub fn tautomer(&self)->Option<TautomerPayload>{match self.inner.as_ref(){ck::OperationError::Tautomer(v)=>Some(crate::tautomer_errors::TautomerRunError::from_rust(v.clone()).unchecked_into()),_=>None}}
    #[cfg(feature = "cap-kekulize")]
    #[wasm_bindgen(getter)] pub fn kekulize(&self) -> Option<KekulizePayload> { match self.inner.as_ref() { ck::OperationError::Kekulize(value) => Some(crate::kekulize_errors::KekulizeError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter, js_name = invalidTopology)] pub fn invalid_topology(&self) -> Option<InvalidTopologyPayload> { match self.inner.as_ref() { ck::OperationError::InvalidTopology(value) => Some(TopologyValidationError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter, js_name = invalidTopologyEdit)] pub fn invalid_topology_edit(&self) -> Option<InvalidTopologyEditPayload> { match self.inner.as_ref() { ck::OperationError::InvalidTopologyEdit(value) => Some(TopologyEditError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter, js_name = invalidCoordinates)] pub fn invalid_coordinates(&self) -> Option<CoordinatesPayload> { match self.inner.as_ref() { ck::OperationError::InvalidCoordinates(value) => Some(CoordinateValidationError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter, js_name = invalidProperty)] pub fn invalid_property(&self) -> Option<PropertyPayload> { match self.inner.as_ref() { ck::OperationError::InvalidProperty(value) => Some(MoleculePropertyError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter, js_name = atomProperty)] pub fn atom_property(&self) -> Option<AtomPropertyPayload> { match self.inner.as_ref() { ck::OperationError::AtomProperty(value) => Some(AtomPropertyError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter, js_name = bondProperty)] pub fn bond_property(&self) -> Option<BondPropertyPayload> { match self.inner.as_ref() { ck::OperationError::BondProperty(value) => Some(BondValueError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[cfg(any(
        feature = "cap-inchi",
        feature = "cap-forcefields",
        feature = "cap-rings",
        feature = "cap-stereo",
        feature = "cap-aromaticity",
        feature = "cap-fingerprints",
        feature = "cap-serialization",
        feature = "cap-io"
    ))]
    #[wasm_bindgen(getter)] pub fn rings(&self) -> Option<RingsPayload> { match self.inner.as_ref() { ck::OperationError::Rings(value) => Some(RingFindingError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn atom(&self) -> Option<usize> { match *self.inner { ck::OperationError::AtomPropertyIndex { atom, .. } => Some(atom.index()), _ => None } }
    #[wasm_bindgen(getter, js_name = atomCount)] pub fn atom_count(&self) -> Option<usize> { match *self.inner { ck::OperationError::AtomPropertyIndex { atom_count, .. } => Some(atom_count), _ => None } }
    #[wasm_bindgen(getter)] pub fn key(&self) -> Option<String> { match self.inner.as_ref() { ck::OperationError::ReservedAtomPropertyKey { key } => Some(key.clone()), _ => None } }
    #[wasm_bindgen(getter)] pub fn operation(&self) -> Option<OperationName> { match *self.inner {
        ck::OperationError::UnsupportedFeature { operation, .. } => Some(JsValue::from(MoleculeOpSpec(operation)).unchecked_into()),
        ck::OperationError::Unsupported { operation } | ck::OperationError::OutputMismatch { operation, .. } | ck::OperationError::AccessDenied { operation, .. }
        | ck::OperationError::BlockCheckedOut { operation, .. } | ck::OperationError::BlockNotCheckedOut { operation, .. } | ck::OperationError::IncompleteCommit { operation, .. }
        | ck::OperationError::MappingContract { operation, .. } | ck::OperationError::InvalidTopologyMapping { operation, .. } | ck::OperationError::AutoRemapContract { operation, .. }
        | ck::OperationError::OperationContract { operation, .. } | ck::OperationError::SemanticPreconditionContract { operation, .. } | ck::OperationError::CoordinateAppendRequiresValues { operation }
        | ck::OperationError::DerivedEffectContract { operation, .. } | ck::OperationError::DerivedPreservationMismatch { operation, .. } | ck::OperationError::CipStateContract { operation, .. }
        | ck::OperationError::InvalidAlgorithmResult { operation, .. } | ck::OperationError::InvalidReconstructionOrigin { operation, .. } | ck::OperationError::Algorithm { operation, .. } => Some(JsValue::from_str(operation).unchecked_into()), _ => None,
    } }
    #[wasm_bindgen(getter)] pub fn block(&self) -> Option<String> { match *self.inner {
        ck::OperationError::AccessDenied { block, .. } | ck::OperationError::BlockCheckedOut { block, .. } | ck::OperationError::BlockNotCheckedOut { block, .. }
        | ck::OperationError::IncompleteCommit { block, .. } | ck::OperationError::AutoRemapContract { block, .. } => Some(block.into()), _ => None,
    } }
    #[wasm_bindgen(getter)] pub fn issue(&self) -> Option<String> { match *self.inner {
        ck::OperationError::MappingContract { issue, .. } | ck::OperationError::AutoRemapContract { issue, .. } | ck::OperationError::OperationContract { issue, .. }
        | ck::OperationError::SemanticPreconditionContract { issue, .. } | ck::OperationError::DerivedEffectContract { issue, .. } | ck::OperationError::CipStateContract { issue, .. } => Some(issue.into()), _ => None,
    } }
    #[wasm_bindgen(getter)] pub fn field(&self) -> Option<String> { match self.inner.as_ref() {
        ck::OperationError::DerivedPreservationMismatch { field, .. } => Some(field.clone()),
        ck::OperationError::OperationContract { field, .. } | ck::OperationError::InvalidDerivedCache { field, .. } | ck::OperationError::InvalidAlgorithmResult { field, .. } => Some((*field).into()), _ => None,
    } }
    #[wasm_bindgen(getter)] pub fn expected(&self) -> Option<OperationScalar> { match self.inner.as_ref() {
        ck::OperationError::InvalidDerivedCache { expected, .. } | ck::OperationError::InvalidAlgorithmResult { expected, .. } | ck::OperationError::InvalidPropertyList { expected, .. } => Some(JsValue::from(*expected).unchecked_into()),
        ck::OperationError::OperationContract { expected, .. } => Some(JsValue::from(*expected).unchecked_into()),
        ck::OperationError::DerivedPreservationMismatch { expected, .. } => Some(JsValue::from_str(expected).unchecked_into()),
        ck::OperationError::OutputMismatch { expected, .. } => Some(JsValue::from(crate::enums::MoleculeOpOutput::from(*expected)).unchecked_into()), _ => None,
    } }
    #[wasm_bindgen(getter)] pub fn actual(&self) -> Option<OperationScalar> { match self.inner.as_ref() {
        ck::OperationError::InvalidDerivedCache { actual, .. } | ck::OperationError::InvalidAlgorithmResult { actual, .. } => Some(JsValue::from(*actual).unchecked_into()),
        ck::OperationError::OperationContract { actual, .. } => Some(JsValue::from(*actual).unchecked_into()),
        ck::OperationError::DerivedPreservationMismatch { actual, .. } => Some(JsValue::from_str(actual).unchecked_into()),
        ck::OperationError::OutputMismatch { actual, .. } => Some(JsValue::from(crate::enums::MoleculeOpOutput::from(*actual)).unchecked_into()), _ => None,
    } }
    #[wasm_bindgen(getter)] pub fn name(&self) -> Result<Option<String>, JsValue> { match self.inner.as_ref() {
        ck::OperationError::InvalidPropertyList { name, .. } => std::str::from_utf8(name.as_bytes()).map(|value| Some(value.into())).map_err(|error| js_sys::TypeError::new(&error.to_string()).into()), _ => Ok(None),
    } }
    #[wasm_bindgen(getter)] pub fn source(&self) -> Option<OperationSource> { match self.inner.as_ref() {
        ck::OperationError::UnsupportedFeature { source, .. } => Some(UnsupportedFeatureError::from_rust(source.clone()).unchecked_into()),
        ck::OperationError::InvalidTopologyMapping { source, .. } => Some(MappingValidationError::from_rust(source.clone()).unchecked_into()), _ => None,
    } }
    #[wasm_bindgen(getter)] pub fn states(&self) -> Option<OperationDerivedState> { match *self.inner {
        ck::OperationError::DerivedEffectContract { .. } | ck::OperationError::DerivedPreservationMismatch { .. } => Some(OperationDerivedState(self.inner.clone())), _ => None,
    } }
    #[wasm_bindgen(getter)] pub fn missing(&self) -> Option<OperationSemanticPreconditions> { match *self.inner { ck::OperationError::SemanticPreconditionContract { .. } => Some(OperationSemanticPreconditions(self.inner.clone())), _ => None } }
    #[wasm_bindgen(getter)] pub fn policy(&self) -> Option<OperationCipStatePolicy> { match *self.inner { ck::OperationError::CipStateContract { .. } => Some(OperationCipStatePolicy(self.inner.clone())), _ => None } }
    #[wasm_bindgen(getter)] pub fn requirement(&self) -> Option<OperationMappingRequirement> { match *self.inner { ck::OperationError::MappingContract { .. } => Some(OperationMappingRequirement(self.inner.clone())), _ => None } }
    #[wasm_bindgen(getter)] pub fn target(&self) -> Option<String> { match *self.inner { ck::OperationError::InvalidPropertyList { target, .. } => Some(target.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn values(&self) -> Option<usize> { match *self.inner { ck::OperationError::InvalidPropertyList { values, .. } => Some(values), _ => None } }
    #[wasm_bindgen(getter)] pub fn state(&self) -> Option<String> { match *self.inner { ck::OperationError::InvalidDerivedCache { state, .. } => Some(state.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn action(&self) -> Option<String> { match *self.inner { ck::OperationError::DerivedEffectContract { action, .. } => Some(action.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn entity(&self) -> Option<String> { match *self.inner { ck::OperationError::InvalidReconstructionOrigin { entity, .. } => Some(entity.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn destination(&self) -> Option<usize> { match *self.inner { ck::OperationError::InvalidReconstructionOrigin { destination, .. } => Some(destination), _ => None } }
    #[wasm_bindgen(getter)] pub fn input(&self) -> Option<usize> { match *self.inner { ck::OperationError::InvalidReconstructionOrigin { input, .. } => Some(input), _ => None } }
    #[wasm_bindgen(getter)] pub fn row(&self) -> Option<usize> { match *self.inner { ck::OperationError::InvalidReconstructionOrigin { row, .. } => Some(row), _ => None } }
    #[wasm_bindgen(getter, js_name = inputCount)] pub fn input_count(&self) -> Option<usize> { match *self.inner { ck::OperationError::InvalidReconstructionOrigin { input_count, .. } => Some(input_count), _ => None } }
    #[wasm_bindgen(getter, js_name = rowCount)] pub fn row_count(&self) -> Option<usize> { match *self.inner { ck::OperationError::InvalidReconstructionOrigin { row_count, .. } => row_count, _ => None } }
    #[wasm_bindgen(getter)] pub fn detail(&self) -> Option<String> { match self.inner.as_ref() { ck::OperationError::Algorithm { detail, .. } => Some(detail.clone()), _ => None } }
}, source = OperationCause, cause = cause);
