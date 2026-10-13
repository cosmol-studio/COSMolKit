//! CK-owned archive 2.0 wire schema, independent of runtime/model layout.
//!
//! New writes use codec 3 / ck-wire-v1 (musli 0.1.9 wire, Binary mode).
//! Codec 2 / ck-storage-v1 archives written by 0.5.0 remain readable.
//! The implementation version is not part of the persistent codec identity.
//! Field and variant identities below are permanent; retired identities must not be reused.
//! No legacy raw or canonical companion is emitted by this writer.

use super::*;
use cosmolkit_core::{
    PersistedRingDecomposition, PersistedRingEdgeId, PersistedRingGraph,
    ring_persistence::{ComponentState, ComponentView, CycleFamilyState},
};
use musli::mode::Binary;
use musli::{Decode, Encode};

const MAJOR: u16 = 2;
pub(super) const MAGIC: &[u8; 8] = b"COSMOL\0\0";
const MINOR: u16 = 0;
const CODEC: u8 = 3;
const CODEC_CONTRACT: &str = "ck-wire-v1";
const STORAGE_CODEC: u8 = 2;
const STORAGE_CODEC_CONTRACT: &str = "ck-storage-v1";
const METADATA_SCHEMA: u16 = 1;
const MOLECULE_SCHEMA: u16 = 1;
const DERIVED_SCHEMA: u16 = 1;
const MAX_ARCHIVE_BYTES: usize = 256 * 1024 * 1024;
const MAX_ROWS: usize = 1_000_000;

#[derive(Debug, Clone, Copy, PartialEq, Encode, Decode)]
#[musli(Binary, name(type = u8))]
enum ChiralRecord {
    #[musli(Binary, name = 0)]
    Unspecified,
    #[musli(Binary, name = 1)]
    TetrahedralCw,
    #[musli(Binary, name = 2)]
    TetrahedralCcw,
    #[musli(Binary, name = 3)]
    Other,
    #[musli(Binary, name = 4)]
    Tetrahedral,
    #[musli(Binary, name = 5)]
    Allene,
    #[musli(Binary, name = 6)]
    SquarePlanar,
    #[musli(Binary, name = 7)]
    TrigonalBipyramidal,
    #[musli(Binary, name = 8)]
    Octahedral,
}
impl From<ChiralTag> for ChiralRecord {
    fn from(value: ChiralTag) -> Self {
        match value {
            ChiralTag::Unspecified => Self::Unspecified,
            ChiralTag::TetrahedralCw => Self::TetrahedralCw,
            ChiralTag::TetrahedralCcw => Self::TetrahedralCcw,
            ChiralTag::Other => Self::Other,
            ChiralTag::Tetrahedral => Self::Tetrahedral,
            ChiralTag::Allene => Self::Allene,
            ChiralTag::SquarePlanar => Self::SquarePlanar,
            ChiralTag::TrigonalBipyramidal => Self::TrigonalBipyramidal,
            ChiralTag::Octahedral => Self::Octahedral,
        }
    }
}
impl From<ChiralRecord> for ChiralTag {
    fn from(value: ChiralRecord) -> Self {
        match value {
            ChiralRecord::Unspecified => Self::Unspecified,
            ChiralRecord::TetrahedralCw => Self::TetrahedralCw,
            ChiralRecord::TetrahedralCcw => Self::TetrahedralCcw,
            ChiralRecord::Other => Self::Other,
            ChiralRecord::Tetrahedral => Self::Tetrahedral,
            ChiralRecord::Allene => Self::Allene,
            ChiralRecord::SquarePlanar => Self::SquarePlanar,
            ChiralRecord::TrigonalBipyramidal => Self::TrigonalBipyramidal,
            ChiralRecord::Octahedral => Self::Octahedral,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Encode, Decode)]
#[musli(Binary, name(type = u8))]
enum HybridRecord {
    #[musli(Binary, name = 0)]
    Unspecified,
    #[musli(Binary, name = 1)]
    S,
    #[musli(Binary, name = 2)]
    Sp,
    #[musli(Binary, name = 3)]
    Sp2,
    #[musli(Binary, name = 4)]
    Sp3,
    #[musli(Binary, name = 5)]
    Sp2d,
    #[musli(Binary, name = 6)]
    Sp3d,
    #[musli(Binary, name = 7)]
    Sp3d2,
    #[musli(Binary, name = 8)]
    Other,
}
impl From<Hybridization> for HybridRecord {
    fn from(value: Hybridization) -> Self {
        match value {
            Hybridization::Unspecified => Self::Unspecified,
            Hybridization::S => Self::S,
            Hybridization::Sp => Self::Sp,
            Hybridization::Sp2 => Self::Sp2,
            Hybridization::Sp3 => Self::Sp3,
            Hybridization::Sp2d => Self::Sp2d,
            Hybridization::Sp3d => Self::Sp3d,
            Hybridization::Sp3d2 => Self::Sp3d2,
            Hybridization::Other => Self::Other,
        }
    }
}
impl From<HybridRecord> for Hybridization {
    fn from(value: HybridRecord) -> Self {
        match value {
            HybridRecord::Unspecified => Self::Unspecified,
            HybridRecord::S => Self::S,
            HybridRecord::Sp => Self::Sp,
            HybridRecord::Sp2 => Self::Sp2,
            HybridRecord::Sp3 => Self::Sp3,
            HybridRecord::Sp2d => Self::Sp2d,
            HybridRecord::Sp3d => Self::Sp3d,
            HybridRecord::Sp3d2 => Self::Sp3d2,
            HybridRecord::Other => Self::Other,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Encode, Decode)]
#[musli(Binary, name(type = u8))]
enum BondOrderRecord {
    #[musli(Binary, name = 0)]
    Unspecified,
    #[musli(Binary, name = 1)]
    Single,
    #[musli(Binary, name = 2)]
    Double,
    #[musli(Binary, name = 3)]
    Triple,
    #[musli(Binary, name = 4)]
    Quadruple,
    #[musli(Binary, name = 5)]
    Quintuple,
    #[musli(Binary, name = 6)]
    Hextuple,
    #[musli(Binary, name = 7)]
    OneAndHalf,
    #[musli(Binary, name = 8)]
    TwoAndHalf,
    #[musli(Binary, name = 9)]
    ThreeAndHalf,
    #[musli(Binary, name = 10)]
    FourAndHalf,
    #[musli(Binary, name = 11)]
    FiveAndHalf,
    #[musli(Binary, name = 12)]
    Aromatic,
    #[musli(Binary, name = 13)]
    Ionic,
    #[musli(Binary, name = 14)]
    Dative,
    #[musli(Binary, name = 15)]
    DativeOne,
    #[musli(Binary, name = 16)]
    DativeLeft,
    #[musli(Binary, name = 17)]
    DativeRight,
    #[musli(Binary, name = 18)]
    Hydrogen,
    #[musli(Binary, name = 19)]
    ThreeCenter,
    #[musli(Binary, name = 20)]
    Other,
    #[musli(Binary, name = 21)]
    Zero,
}
impl From<BondOrder> for BondOrderRecord {
    fn from(value: BondOrder) -> Self {
        match value {
            BondOrder::Unspecified => Self::Unspecified,
            BondOrder::Single => Self::Single,
            BondOrder::Double => Self::Double,
            BondOrder::Triple => Self::Triple,
            BondOrder::Quadruple => Self::Quadruple,
            BondOrder::Quintuple => Self::Quintuple,
            BondOrder::Hextuple => Self::Hextuple,
            BondOrder::OneAndHalf => Self::OneAndHalf,
            BondOrder::TwoAndHalf => Self::TwoAndHalf,
            BondOrder::ThreeAndHalf => Self::ThreeAndHalf,
            BondOrder::FourAndHalf => Self::FourAndHalf,
            BondOrder::FiveAndHalf => Self::FiveAndHalf,
            BondOrder::Aromatic => Self::Aromatic,
            BondOrder::Ionic => Self::Ionic,
            BondOrder::Dative => Self::Dative,
            BondOrder::DativeOne => Self::DativeOne,
            BondOrder::DativeLeft => Self::DativeLeft,
            BondOrder::DativeRight => Self::DativeRight,
            BondOrder::Hydrogen => Self::Hydrogen,
            BondOrder::ThreeCenter => Self::ThreeCenter,
            BondOrder::Other => Self::Other,
            BondOrder::Zero => Self::Zero,
        }
    }
}
impl From<BondOrderRecord> for BondOrder {
    fn from(value: BondOrderRecord) -> Self {
        match value {
            BondOrderRecord::Unspecified => Self::Unspecified,
            BondOrderRecord::Single => Self::Single,
            BondOrderRecord::Double => Self::Double,
            BondOrderRecord::Triple => Self::Triple,
            BondOrderRecord::Quadruple => Self::Quadruple,
            BondOrderRecord::Quintuple => Self::Quintuple,
            BondOrderRecord::Hextuple => Self::Hextuple,
            BondOrderRecord::OneAndHalf => Self::OneAndHalf,
            BondOrderRecord::TwoAndHalf => Self::TwoAndHalf,
            BondOrderRecord::ThreeAndHalf => Self::ThreeAndHalf,
            BondOrderRecord::FourAndHalf => Self::FourAndHalf,
            BondOrderRecord::FiveAndHalf => Self::FiveAndHalf,
            BondOrderRecord::Aromatic => Self::Aromatic,
            BondOrderRecord::Ionic => Self::Ionic,
            BondOrderRecord::Dative => Self::Dative,
            BondOrderRecord::DativeOne => Self::DativeOne,
            BondOrderRecord::DativeLeft => Self::DativeLeft,
            BondOrderRecord::DativeRight => Self::DativeRight,
            BondOrderRecord::Hydrogen => Self::Hydrogen,
            BondOrderRecord::ThreeCenter => Self::ThreeCenter,
            BondOrderRecord::Other => Self::Other,
            BondOrderRecord::Zero => Self::Zero,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Encode, Decode)]
#[musli(Binary, name(type = u8))]
enum BondDirectionRecord {
    #[musli(Binary, name = 0)]
    None,
    #[musli(Binary, name = 1)]
    BeginWedge,
    #[musli(Binary, name = 2)]
    BeginDash,
    #[musli(Binary, name = 3)]
    EndUpRight,
    #[musli(Binary, name = 4)]
    EndDownRight,
    #[musli(Binary, name = 5)]
    EitherDouble,
    #[musli(Binary, name = 6)]
    Unknown,
}
impl From<BondDirection> for BondDirectionRecord {
    fn from(value: BondDirection) -> Self {
        match value {
            BondDirection::None => Self::None,
            BondDirection::BeginWedge => Self::BeginWedge,
            BondDirection::BeginDash => Self::BeginDash,
            BondDirection::EndUpRight => Self::EndUpRight,
            BondDirection::EndDownRight => Self::EndDownRight,
            BondDirection::EitherDouble => Self::EitherDouble,
            BondDirection::Unknown => Self::Unknown,
        }
    }
}
impl From<BondDirectionRecord> for BondDirection {
    fn from(value: BondDirectionRecord) -> Self {
        match value {
            BondDirectionRecord::None => Self::None,
            BondDirectionRecord::BeginWedge => Self::BeginWedge,
            BondDirectionRecord::BeginDash => Self::BeginDash,
            BondDirectionRecord::EndUpRight => Self::EndUpRight,
            BondDirectionRecord::EndDownRight => Self::EndDownRight,
            BondDirectionRecord::EitherDouble => Self::EitherDouble,
            BondDirectionRecord::Unknown => Self::Unknown,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Encode, Decode)]
#[musli(Binary, name(type = u8))]
enum BondStereoRecord {
    #[musli(Binary, name = 0)]
    None,
    #[musli(Binary, name = 1)]
    Any,
    #[musli(Binary, name = 2)]
    Z,
    #[musli(Binary, name = 3)]
    E,
    #[musli(Binary, name = 4)]
    Cis,
    #[musli(Binary, name = 5)]
    Trans,
    #[musli(Binary, name = 6)]
    AtropCw,
    #[musli(Binary, name = 7)]
    AtropCcw,
}
impl From<BondStereo> for BondStereoRecord {
    fn from(value: BondStereo) -> Self {
        match value {
            BondStereo::None => Self::None,
            BondStereo::Any => Self::Any,
            BondStereo::Z => Self::Z,
            BondStereo::E => Self::E,
            BondStereo::Cis => Self::Cis,
            BondStereo::Trans => Self::Trans,
            BondStereo::AtropCw => Self::AtropCw,
            BondStereo::AtropCcw => Self::AtropCcw,
        }
    }
}
impl From<BondStereoRecord> for BondStereo {
    fn from(value: BondStereoRecord) -> Self {
        match value {
            BondStereoRecord::None => Self::None,
            BondStereoRecord::Any => Self::Any,
            BondStereoRecord::Z => Self::Z,
            BondStereoRecord::E => Self::E,
            BondStereoRecord::Cis => Self::Cis,
            BondStereoRecord::Trans => Self::Trans,
            BondStereoRecord::AtropCw => Self::AtropCw,
            BondStereoRecord::AtropCcw => Self::AtropCcw,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Encode, Decode)]
#[musli(Binary, name(type = u8))]
enum RingFindRecord {
    #[musli(Binary, name = 0)]
    OtherOrUnknown,
    #[musli(Binary, name = 1)]
    Fast,
    #[musli(Binary, name = 2)]
    Sssr,
    #[musli(Binary, name = 3)]
    SymmSssr,
}
impl From<RingFindType> for RingFindRecord {
    fn from(value: RingFindType) -> Self {
        match value {
            RingFindType::OtherOrUnknown => Self::OtherOrUnknown,
            RingFindType::Fast => Self::Fast,
            RingFindType::Sssr => Self::Sssr,
            RingFindType::SymmSssr => Self::SymmSssr,
        }
    }
}
impl From<RingFindRecord> for RingFindType {
    fn from(value: RingFindRecord) -> Self {
        match value {
            RingFindRecord::OtherOrUnknown => Self::OtherOrUnknown,
            RingFindRecord::Fast => Self::Fast,
            RingFindRecord::Sssr => Self::Sssr,
            RingFindRecord::SymmSssr => Self::SymmSssr,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Encode, Decode)]
#[musli(Binary, name(type = u8))]
enum StereoGroupKindRecord {
    #[musli(Binary, name = 0)]
    Absolute,
    #[musli(Binary, name = 1)]
    Or,
    #[musli(Binary, name = 2)]
    And,
}
impl From<StereoGroupKind> for StereoGroupKindRecord {
    fn from(value: StereoGroupKind) -> Self {
        match value {
            StereoGroupKind::Absolute => Self::Absolute,
            StereoGroupKind::Or => Self::Or,
            StereoGroupKind::And => Self::And,
        }
    }
}
impl From<StereoGroupKindRecord> for StereoGroupKind {
    fn from(value: StereoGroupKindRecord) -> Self {
        match value {
            StereoGroupKindRecord::Absolute => Self::Absolute,
            StereoGroupKindRecord::Or => Self::Or,
            StereoGroupKindRecord::And => Self::And,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Encode, Decode)]
#[musli(Binary, name(type = u8))]
enum DimensionRecord {
    #[musli(Binary, name = 2)]
    TwoD,
    #[musli(Binary, name = 3)]
    ThreeD,
}
impl From<CoordinateDimension> for DimensionRecord {
    fn from(value: CoordinateDimension) -> Self {
        match value {
            CoordinateDimension::TwoD => Self::TwoD,
            CoordinateDimension::ThreeD => Self::ThreeD,
        }
    }
}
impl From<DimensionRecord> for CoordinateDimension {
    fn from(value: DimensionRecord) -> Self {
        match value {
            DimensionRecord::TwoD => Self::TwoD,
            DimensionRecord::ThreeD => Self::ThreeD,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Encode, Decode)]
#[musli(Binary, name(type = u8))]
enum BondRoleRecord {
    #[musli(Binary, name = 0)]
    Crossing,
    #[musli(Binary, name = 1)]
    Contained,
}
impl From<SGroupBondRole> for BondRoleRecord {
    fn from(value: SGroupBondRole) -> Self {
        match value {
            SGroupBondRole::Crossing => Self::Crossing,
            SGroupBondRole::Contained => Self::Contained,
        }
    }
}
impl From<BondRoleRecord> for SGroupBondRole {
    fn from(value: BondRoleRecord) -> Self {
        match value {
            BondRoleRecord::Crossing => Self::Crossing,
            BondRoleRecord::Contained => Self::Contained,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Encode, Decode)]
#[musli(Binary, name(type = u8))]
enum PropertyTargetRecord {
    #[musli(Binary, name = 0)]
    Atom,
    #[musli(Binary, name = 1)]
    Bond,
}
impl From<SdfPropertyListTarget> for PropertyTargetRecord {
    fn from(value: SdfPropertyListTarget) -> Self {
        match value {
            SdfPropertyListTarget::Atom => Self::Atom,
            SdfPropertyListTarget::Bond => Self::Bond,
        }
    }
}
impl From<PropertyTargetRecord> for SdfPropertyListTarget {
    fn from(value: PropertyTargetRecord) -> Self {
        match value {
            PropertyTargetRecord::Atom => Self::Atom,
            PropertyTargetRecord::Bond => Self::Bond,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Encode, Decode)]
#[musli(Binary, name(type = u8))]
enum SGroupKindRecord {
    #[musli(Binary, name = 0)]
    Data,
    #[musli(Binary, name = 1)]
    Superatom,
    #[musli(Binary, name = 2)]
    MultipleGroup,
    #[musli(Binary, name = 3)]
    StructuralRepeatUnit,
    #[musli(Binary, name = 4)]
    Monomer,
    #[musli(Binary, name = 5)]
    Copolymer,
    #[musli(Binary, name = 6)]
    Crosslink,
    #[musli(Binary, name = 7)]
    Graft,
    #[musli(Binary, name = 8)]
    Modification,
    #[musli(Binary, name = 9)]
    Mer,
    #[musli(Binary, name = 10)]
    AnyPolymer,
    #[musli(Binary, name = 11)]
    MixtureComponent,
    #[musli(Binary, name = 12)]
    Mixture,
    #[musli(Binary, name = 13)]
    Formulation,
    #[musli(Binary, name = 14, packed)]
    Generic(Vec<u8>),
}
impl From<&SubstanceGroupKind> for SGroupKindRecord {
    fn from(value: &SubstanceGroupKind) -> Self {
        match value {
            SubstanceGroupKind::Data => Self::Data,
            SubstanceGroupKind::Superatom => Self::Superatom,
            SubstanceGroupKind::MultipleGroup => Self::MultipleGroup,
            SubstanceGroupKind::StructuralRepeatUnit => Self::StructuralRepeatUnit,
            SubstanceGroupKind::Monomer => Self::Monomer,
            SubstanceGroupKind::Copolymer => Self::Copolymer,
            SubstanceGroupKind::Crosslink => Self::Crosslink,
            SubstanceGroupKind::Graft => Self::Graft,
            SubstanceGroupKind::Modification => Self::Modification,
            SubstanceGroupKind::Mer => Self::Mer,
            SubstanceGroupKind::AnyPolymer => Self::AnyPolymer,
            SubstanceGroupKind::MixtureComponent => Self::MixtureComponent,
            SubstanceGroupKind::Mixture => Self::Mixture,
            SubstanceGroupKind::Formulation => Self::Formulation,
            SubstanceGroupKind::Generic(value) => Self::Generic(value.as_bytes().to_vec()),
        }
    }
}
impl From<SGroupKindRecord> for SubstanceGroupKind {
    fn from(value: SGroupKindRecord) -> Self {
        match value {
            SGroupKindRecord::Data => Self::Data,
            SGroupKindRecord::Superatom => Self::Superatom,
            SGroupKindRecord::MultipleGroup => Self::MultipleGroup,
            SGroupKindRecord::StructuralRepeatUnit => Self::StructuralRepeatUnit,
            SGroupKindRecord::Monomer => Self::Monomer,
            SGroupKindRecord::Copolymer => Self::Copolymer,
            SGroupKindRecord::Crosslink => Self::Crosslink,
            SGroupKindRecord::Graft => Self::Graft,
            SGroupKindRecord::Modification => Self::Modification,
            SGroupKindRecord::Mer => Self::Mer,
            SGroupKindRecord::AnyPolymer => Self::AnyPolymer,
            SGroupKindRecord::MixtureComponent => Self::MixtureComponent,
            SGroupKindRecord::Mixture => Self::Mixture,
            SGroupKindRecord::Formulation => Self::Formulation,
            SGroupKindRecord::Generic(value) => Self::Generic(value.into()),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Encode, Decode)]
#[musli(Binary, name(type = u8))]
enum ConnectionRecord {
    #[musli(Binary, name = 0)]
    HeadToHead,
    #[musli(Binary, name = 1)]
    HeadToTail,
    #[musli(Binary, name = 2)]
    Either,
    #[musli(Binary, name = 3, packed)]
    Unknown(Vec<u8>),
}
impl From<&SGroupConnection> for ConnectionRecord {
    fn from(value: &SGroupConnection) -> Self {
        match value {
            SGroupConnection::HeadToHead => Self::HeadToHead,
            SGroupConnection::HeadToTail => Self::HeadToTail,
            SGroupConnection::Either => Self::Either,
            SGroupConnection::Unknown(value) => Self::Unknown(value.as_bytes().to_vec()),
        }
    }
}
impl From<ConnectionRecord> for SGroupConnection {
    fn from(value: ConnectionRecord) -> Self {
        match value {
            ConnectionRecord::HeadToHead => Self::HeadToHead,
            ConnectionRecord::HeadToTail => Self::HeadToTail,
            ConnectionRecord::Either => Self::Either,
            ConnectionRecord::Unknown(value) => Self::Unknown(value.into()),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Encode, Decode)]
#[musli(Binary, name(type = u8))]
enum BracketRecord {
    #[musli(Binary, name = 0)]
    Bracket,
    #[musli(Binary, name = 1)]
    Parenthesis,
    #[musli(Binary, name = 2)]
    None,
    #[musli(Binary, name = 3, packed)]
    Unknown(Vec<u8>),
}
impl From<&SGroupBracketStyle> for BracketRecord {
    fn from(value: &SGroupBracketStyle) -> Self {
        match value {
            SGroupBracketStyle::Bracket => Self::Bracket,
            SGroupBracketStyle::Parenthesis => Self::Parenthesis,
            SGroupBracketStyle::None => Self::None,
            SGroupBracketStyle::Unknown(value) => Self::Unknown(value.as_bytes().to_vec()),
        }
    }
}
impl From<BracketRecord> for SGroupBracketStyle {
    fn from(value: BracketRecord) -> Self {
        match value {
            BracketRecord::Bracket => Self::Bracket,
            BracketRecord::Parenthesis => Self::Parenthesis,
            BracketRecord::None => Self::None,
            BracketRecord::Unknown(value) => Self::Unknown(value.into()),
        }
    }
}

#[derive(Debug, Clone, Encode, Decode)]
#[musli(Binary, name(type = u8))]
enum Value {
    #[musli(Binary, name = 0, packed)]
    String(#[musli(bytes)] Vec<u8>),
    #[musli(Binary, name = 1, packed)]
    Int(i32),
    #[musli(Binary, name = 2, packed)]
    UInt(u32),
    #[musli(Binary, name = 3, packed)]
    IntVector(Vec<i32>),
    // Integers carry float bits, including signaling NaNs and signed zero.
    #[musli(Binary, name = 4, packed)]
    Double(u64),
    #[musli(Binary, name = 5, packed)]
    Bool(bool),
    #[musli(Binary, name = 6, packed)]
    StringVector(Vec<Vec<u8>>),
}
impl From<&PropertyValue> for Value {
    fn from(value: &PropertyValue) -> Self {
        match value {
            PropertyValue::String(s) => Self::String(s.as_bytes().to_vec()),
            PropertyValue::Int(v) => Self::Int(*v),
            PropertyValue::UInt(v) => Self::UInt(*v),
            PropertyValue::IntVector(v) => Self::IntVector(v.clone()),
            PropertyValue::Double(v) => Self::Double(v.to_bits()),
            PropertyValue::Bool(v) => Self::Bool(*v),
            PropertyValue::StringVector(v) => {
                Self::StringVector(v.iter().map(|s| s.as_bytes().to_vec()).collect())
            }
        }
    }
}
impl Value {
    fn into_model(self) -> Result<PropertyValue, PickleError> {
        Ok(match self {
            Self::String(v) => PropertyValue::String(v.into()),
            Self::Int(v) => PropertyValue::Int(v),
            Self::UInt(v) => PropertyValue::UInt(v),
            Self::IntVector(v) => PropertyValue::IntVector(v),
            Self::Double(v) => PropertyValue::Double(f64::from_bits(v)),
            Self::Bool(v) => PropertyValue::Bool(v),
            Self::StringVector(v) => {
                PropertyValue::StringVector(v.into_iter().map(Into::into).collect())
            }
        })
    }
}

#[derive(Debug, Clone, Encode, Decode)]
#[musli(Binary, name(type = u32))]
struct Metadata {
    #[musli(Binary, name = 0)]
    producer: String,
    #[musli(Binary, name = 1)]
    codec_contract: String,
    #[musli(Binary, name = 2)]
    molecule_schema: u16,
    #[musli(Binary, name = 3)]
    derived_schema: u16,
}
#[derive(Debug, Clone, Encode, Decode)]
#[musli(Binary, name(type = u32))]
struct Property {
    #[musli(Binary, name = 0)]
    key: Vec<u8>,
    #[musli(Binary, name = 1)]
    value: Value,
    #[musli(Binary, name = 2)]
    computed: bool,
}
#[derive(Debug, Clone, Encode, Decode)]
#[musli(Binary, name(type = u32))]
struct PdbInfo {
    #[musli(Binary, name = 0)]
    atom_name: Vec<u8>,
    #[musli(Binary, name = 1)]
    serial_number: i32,
    #[musli(Binary, name = 2)]
    alt_loc: Vec<u8>,
    #[musli(Binary, name = 3)]
    residue_name: Vec<u8>,
    #[musli(Binary, name = 4)]
    residue_number: i32,
    #[musli(Binary, name = 5)]
    chain_id: Vec<u8>,
    #[musli(Binary, name = 6)]
    insertion_code: Vec<u8>,
    #[musli(Binary, name = 7)]
    occupancy_bits: u64,
    #[musli(Binary, name = 8)]
    temp_factor_bits: u64,
    #[musli(Binary, name = 9)]
    hetero: bool,
    #[musli(Binary, name = 10)]
    secondary_structure: u32,
    #[musli(Binary, name = 11)]
    segment_number: u32,
    #[musli(Binary, name = 12)]
    monomer_class: Vec<u8>,
}
#[derive(Debug, Clone, Encode, Decode)]
#[musli(Binary, name(type = u32))]
struct Attachment {
    #[musli(Binary, name = 0)]
    target: u64,
    #[musli(Binary, name = 1)]
    label: String,
}
#[derive(Debug, Clone, Encode, Decode)]
#[musli(Binary, name(type = u32))]
struct AtomRecord {
    #[musli(Binary, name = 0)]
    atomic_number: u8,
    #[musli(Binary, name = 1)]
    formal_charge: i8,
    #[musli(Binary, name = 2)]
    isotope: Option<u16>,
    #[musli(Binary, name = 3)]
    chiral: ChiralRecord,
    #[musli(Binary, name = 4)]
    chiral_permutation: Option<u32>,
    #[musli(Binary, name = 5)]
    unknown_stereo: bool,
    #[musli(Binary, name = 6)]
    mol_parity: Option<i32>,
    #[musli(Binary, name = 7)]
    mol_inversion: Option<i32>,
    #[musli(Binary, name = 8)]
    radical_electrons: u8,
    #[musli(Binary, name = 9)]
    aromatic: bool,
    #[musli(Binary, name = 10)]
    hybridization: HybridRecord,
    #[musli(Binary, name = 11)]
    atom_map: Option<u32>,
    #[musli(Binary, name = 12)]
    no_implicit: bool,
    #[musli(Binary, name = 13)]
    implicit_hydrogen: bool,
    #[musli(Binary, name = 14)]
    explicit_hydrogens: u8,
    #[musli(Binary, name = 15)]
    tracked_isotopes: Vec<u16>,
    #[musli(Binary, name = 16)]
    properties: Vec<Property>,
    #[musli(Binary, name = 17)]
    temporary_flags: u64,
    #[musli(Binary, name = 18)]
    pdb: Option<PdbInfo>,
    #[musli(Binary, name = 19)]
    attachment_order: Option<Vec<Attachment>>,
}
#[derive(Debug, Clone, Encode, Decode)]
#[musli(Binary, name(type = u32))]
struct BondRecord {
    #[musli(Binary, name = 0)]
    begin: u64,
    #[musli(Binary, name = 1)]
    end: u64,
    #[musli(Binary, name = 2)]
    order: BondOrderRecord,
    #[musli(Binary, name = 3)]
    stereo: BondStereoRecord,
    #[musli(Binary, name = 4)]
    direction: BondDirectionRecord,
    #[musli(Binary, name = 5)]
    aromatic: bool,
    #[musli(Binary, name = 6)]
    conjugated: bool,
    #[musli(Binary, name = 7)]
    stereo_atoms: Option<[u64; 2]>,
    #[musli(Binary, name = 8)]
    unknown_stereo: bool,
    #[musli(Binary, name = 9)]
    properties: Vec<Property>,
    #[musli(Binary, name = 10)]
    temporary_flags: u64,
    // Keep tag 7 readable for existing V2 pairs. Source vectors with another
    // cardinality use a new defaulted field instead of losing their entries.
    #[musli(Binary, name = 11, default, skip_encoding_if = Option::is_none)]
    source_stereo_atoms: Option<Vec<u64>>,
}
#[derive(Debug, Clone, Encode, Decode)]
#[musli(Binary, name(type = u32))]
struct Conformer2 {
    #[musli(Binary, name = 0)]
    id: u64,
    #[musli(Binary, name = 1)]
    coordinates: Vec<[u64; 2]>,
    #[musli(Binary, name = 2)]
    properties: Vec<(Vec<u8>, Vec<u8>)>,
}
#[derive(Debug, Clone, Encode, Decode)]
#[musli(Binary, name(type = u32))]
struct Conformer3 {
    #[musli(Binary, name = 0)]
    id: u64,
    #[musli(Binary, name = 1)]
    coordinates: Vec<[u64; 3]>,
    #[musli(Binary, name = 2)]
    is_3d: bool,
    #[musli(Binary, name = 3)]
    properties: Vec<(Vec<u8>, Vec<u8>)>,
}
#[derive(Debug, Clone, Encode, Decode)]
#[musli(Binary, name(type = u32))]
struct GroupDisplay {
    #[musli(Binary, name = 0)]
    brackets: Vec<[[u64; 3]; 3]>,
    #[musli(Binary, name = 1)]
    field_position: Option<[u64; 2]>,
    #[musli(Binary, name = 2)]
    tag: Option<Vec<u8>>,
}
#[derive(Debug, Clone, Encode, Decode)]
#[musli(Binary, name(type = u32))]
struct GroupData {
    #[musli(Binary, name = 0)]
    field_name: Option<Vec<u8>>,
    #[musli(Binary, name = 1)]
    field_type: Option<Vec<u8>>,
    #[musli(Binary, name = 2)]
    field_info: Option<Vec<u8>>,
    #[musli(Binary, name = 3)]
    field_display: Option<Vec<u8>>,
    #[musli(Binary, name = 4)]
    units: Option<Vec<u8>>,
    #[musli(Binary, name = 5)]
    query_type: Option<Vec<u8>>,
    #[musli(Binary, name = 6)]
    query_op: Option<Vec<u8>>,
    #[musli(Binary, name = 7)]
    values: Vec<Vec<u8>>,
}
#[derive(Debug, Clone, Encode, Decode)]
#[musli(Binary, name(type = u32))]
struct AttachPoint {
    #[musli(Binary, name = 0)]
    atom: u64,
    #[musli(Binary, name = 1)]
    leaving: Option<u64>,
    #[musli(Binary, name = 2)]
    label: Option<Vec<u8>>,
    #[musli(Binary, name = 3)]
    order: Option<u32>,
}
#[derive(Debug, Clone, Encode, Decode)]
#[musli(Binary, name(type = u32))]
struct CState {
    #[musli(Binary, name = 0)]
    bond: u64,
    #[musli(Binary, name = 1)]
    vector: [u64; 3],
}
#[derive(Debug, Clone, Encode, Decode)]
#[musli(Binary, name(type = u32))]
struct SGroup {
    #[musli(Binary, name = 0)]
    id: u64,
    #[musli(Binary, name = 1)]
    sequence_id: Option<u32>,
    #[musli(Binary, name = 2)]
    external_id: Option<u32>,
    #[musli(Binary, name = 3)]
    kind: SGroupKindRecord,
    #[musli(Binary, name = 5)]
    atoms: Vec<u64>,
    #[musli(Binary, name = 6)]
    bonds: Vec<u64>,
    #[musli(Binary, name = 7)]
    roles: Vec<(u64, BondRoleRecord)>,
    #[musli(Binary, name = 8)]
    parent_atoms: Vec<u64>,
    #[musli(Binary, name = 9)]
    parent: Option<u64>,
    #[musli(Binary, name = 10)]
    label: Option<Vec<u8>>,
    #[musli(Binary, name = 11)]
    connection: Option<ConnectionRecord>,
    #[musli(Binary, name = 12)]
    subtype: Option<Vec<u8>>,
    #[musli(Binary, name = 13)]
    bracket_style: Option<BracketRecord>,
    #[musli(Binary, name = 14)]
    expansion: Option<Vec<u8>>,
    #[musli(Binary, name = 15)]
    class: Option<Vec<u8>>,
    #[musli(Binary, name = 16)]
    component_number: Option<u32>,
    #[musli(Binary, name = 17)]
    display: Option<GroupDisplay>,
    #[musli(Binary, name = 18)]
    data: Option<GroupData>,
    #[musli(Binary, name = 19)]
    attachments: Vec<AttachPoint>,
    #[musli(Binary, name = 20)]
    cstates: Vec<CState>,
    #[musli(Binary, name = 21)]
    properties: Vec<(Vec<u8>, Value)>,
    #[musli(Binary, name = 22)]
    data_fields: Vec<Vec<u8>>,
    #[musli(Binary, name = 23)]
    head_crossing: Vec<u64>,
    #[musli(Binary, name = 24)]
    correspondence: Vec<u64>,
}
#[derive(Debug, Clone, Encode, Decode)]
#[musli(Binary, name(type = u32))]
struct StereoRecord {
    #[musli(Binary, name = 0)]
    id: Option<u32>,
    #[musli(Binary, name = 1)]
    write_id: u32,
    #[musli(Binary, name = 2)]
    kind: StereoGroupKindRecord,
    #[musli(Binary, name = 3)]
    atoms: Vec<u64>,
    #[musli(Binary, name = 4)]
    bonds: Vec<u64>,
}
#[derive(Debug, Clone, Encode, Decode)]
#[musli(Binary, name(type = u32))]
struct PropertyList {
    #[musli(Binary, name = 0)]
    target: PropertyTargetRecord,
    #[musli(Binary, name = 1)]
    name: Vec<u8>,
    #[musli(Binary, name = 2)]
    values: Vec<Option<Value>>,
}
#[derive(Debug, Clone, Encode, Decode)]
#[musli(Binary, name(type = u32))]
struct MoleculeState {
    #[musli(Binary, name = 0)]
    atoms: Vec<AtomRecord>,
    #[musli(Binary, name = 1)]
    bonds: Vec<BondRecord>,
    #[musli(Binary, name = 2)]
    conformers_2d: Vec<Conformer2>,
    #[musli(Binary, name = 3)]
    conformers_3d: Vec<Conformer3>,
    #[musli(Binary, name = 4)]
    source_dimension: Option<DimensionRecord>,
    #[musli(Binary, name = 5)]
    #[musli(default)]
    source_order: Option<Vec<DimensionRecord>>,
    #[musli(Binary, name = 6)]
    sgroups: Vec<SGroup>,
    #[musli(Binary, name = 7)]
    stereo_groups: Vec<StereoRecord>,
    #[musli(Binary, name = 8)]
    name: Option<Vec<u8>>,
    #[musli(Binary, name = 9)]
    properties: Vec<Property>,
    #[musli(Binary, name = 10)]
    sdf_fields: Vec<(Vec<u8>, Vec<u8>)>,
    #[musli(Binary, name = 11)]
    sdf_lists: Vec<PropertyList>,
}
#[derive(Debug, Clone, Encode, Decode)]
#[musli(Binary, name(type = u32))]
struct RingRecord {
    #[musli(Binary, name = 0)]
    initialized: bool,
    #[musli(Binary, name = 1)]
    find_type: RingFindRecord,
    #[musli(Binary, name = 2)]
    atom_extent: u64,
    #[musli(Binary, name = 3)]
    bond_extent: u64,
    #[musli(Binary, name = 4)]
    atom_rings: Vec<Vec<u64>>,
    #[musli(Binary, name = 5)]
    bond_rings: Vec<Vec<u64>>,
    #[musli(Binary, name = 6)]
    atom_families: Vec<Vec<u64>>,
    #[musli(Binary, name = 7)]
    bond_families: Vec<Vec<u64>>,
    #[musli(Binary, name = 8)]
    relevant_cycles: Option<u64>,
    #[musli(Binary, name = 9)]
    fused_rings: Vec<Vec<bool>>,
    #[musli(Binary, name = 10)]
    fused_bonds: Vec<u64>,
    #[musli(Binary, name = 11, default)]
    decomposition: Option<RdlRecord>,
    #[musli(Binary, name = 12, default)]
    atom_members: Option<Vec<Vec<u64>>>,
    #[musli(Binary, name = 13, default)]
    bond_members: Option<Vec<Vec<u64>>>,
}

#[derive(Debug, Clone, Encode, Decode)]
#[musli(Binary, name(type = u32))]
struct RdlGraph {
    #[musli(Binary, name = 0)]
    nodes: u64,
    #[musli(Binary, name = 1)]
    edges: Vec<[u64; 2]>,
}

#[derive(Debug, Clone, Encode, Decode)]
#[musli(Binary, name(type = u32))]
struct RdlFamily {
    #[musli(Binary, name = 0)]
    weight: u64,
    #[musli(Binary, name = 1)]
    root: u64,
    #[musli(Binary, name = 2)]
    p: u64,
    #[musli(Binary, name = 3)]
    q: u64,
    #[musli(Binary, name = 4)]
    x: Option<u64>,
    #[musli(Binary, name = 5)]
    prototype: Vec<bool>,
    #[musli(Binary, name = 6)]
    relevant: bool,
}

#[derive(Debug, Clone, Encode, Decode)]
#[musli(Binary, name(type = u32))]
struct RdlComponent {
    #[musli(Binary, name = 0)]
    graph: RdlGraph,
    #[musli(Binary, name = 1)]
    original_nodes: Vec<u64>,
    #[musli(Binary, name = 2)]
    original_edges: Vec<u64>,
    #[musli(Binary, name = 3)]
    predecessor: Vec<Vec<Option<u64>>>,
    #[musli(Binary, name = 4)]
    distance: Vec<Vec<Option<u64>>>,
    #[musli(Binary, name = 5)]
    reachable: Vec<Vec<bool>>,
    #[musli(Binary, name = 6)]
    directed_paths: Vec<Vec<Vec<u64>>>,
    #[musli(Binary, name = 7)]
    families: Vec<RdlFamily>,
    #[musli(Binary, name = 8)]
    prototype_counts: Vec<u64>,
    #[musli(Binary, name = 9)]
    relations: Vec<Vec<Vec<bool>>>,
    #[musli(Binary, name = 10)]
    urfs: Vec<Vec<u64>>,
}

#[derive(Debug, Clone, Encode, Decode)]
#[musli(Binary, name(type = u32))]
struct RdlRecord {
    #[musli(Binary, name = 0)]
    graph: RdlGraph,
    #[musli(Binary, name = 1)]
    family_nodes: Vec<Vec<u64>>,
    #[musli(Binary, name = 2)]
    family_edges: Vec<Vec<u64>>,
    #[musli(Binary, name = 3)]
    relevant_cycle_count_bits: u64,
    #[musli(Binary, name = 4)]
    components: Vec<RdlComponent>,
}

fn index_rows(values: Vec<Vec<u64>>) -> Result<Vec<Vec<usize>>, PickleError> {
    rows(values.len())?;
    values
        .into_iter()
        .map(|row| {
            rows(row.len())?;
            row.into_iter().map(index).collect()
        })
        .collect()
}

fn wire_rows(values: &[Vec<usize>]) -> Vec<Vec<u64>> {
    values
        .iter()
        .map(|row| row.iter().map(|&v| v as u64).collect())
        .collect()
}

impl From<&PersistedRingGraph> for RdlGraph {
    fn from(graph: &PersistedRingGraph) -> Self {
        Self {
            nodes: graph.node_count() as u64,
            edges: graph
                .edges()
                .iter()
                .map(|e| [e.from() as u64, e.to() as u64])
                .collect(),
        }
    }
}

impl RdlGraph {
    fn into_model(self) -> Result<PersistedRingGraph, PickleError> {
        let nodes = index(self.nodes)?;
        rows(nodes)?;
        rows(self.edges.len())?;
        let mut graph = PersistedRingGraph::new(nodes);
        for [from, to] in self.edges {
            graph
                .add_undirected_edge(index(from)?, index(to)?)
                .map_err(invalid)?;
        }
        Ok(graph)
    }
}

impl From<ComponentView<'_>> for RdlComponent {
    fn from(c: ComponentView<'_>) -> Self {
        fn optional_rows(values: &[Vec<Option<usize>>]) -> Vec<Vec<Option<u64>>> {
            values
                .iter()
                .map(|row| row.iter().map(|v| v.map(|v| v as u64)).collect())
                .collect()
        }
        Self {
            graph: c.graph.into(),
            original_nodes: c.original_nodes.iter().map(|&v| v as u64).collect(),
            original_edges: c.original_edges.iter().map(|v| v.index() as u64).collect(),
            predecessor: optional_rows(c.predecessor),
            distance: optional_rows(c.distance),
            reachable: c.reachable.to_vec(),
            directed_paths: c.directed_paths().map(wire_rows).collect(),
            families: c
                .cycle_families()
                .map(|f| RdlFamily {
                    weight: f.weight as u64,
                    root: f.root as u64,
                    p: f.p as u64,
                    q: f.q as u64,
                    x: f.x.map(|v| v as u64),
                    prototype: f.prototype.to_vec(),
                    relevant: f.relevant,
                })
                .collect(),
            prototype_counts: c.prototype_counts.iter().map(|&v| v as u64).collect(),
            relations: c.relations.to_vec(),
            urfs: wire_rows(c.urfs),
        }
    }
}

impl RdlComponent {
    fn into_model(self) -> Result<ComponentState, PickleError> {
        fn optional_rows(
            values: Vec<Vec<Option<u64>>>,
        ) -> Result<Vec<Vec<Option<usize>>>, PickleError> {
            rows(values.len())?;
            values
                .into_iter()
                .map(|row| {
                    rows(row.len())?;
                    row.into_iter().map(|v| v.map(index).transpose()).collect()
                })
                .collect()
        }
        for count in [
            self.original_nodes.len(),
            self.original_edges.len(),
            self.directed_paths.len(),
            self.families.len(),
            self.prototype_counts.len(),
            self.reachable.len(),
            self.relations.len(),
        ] {
            rows(count)?;
        }
        for row in &self.reachable {
            rows(row.len())?;
        }
        for matrix in &self.relations {
            rows(matrix.len())?;
            for row in matrix {
                rows(row.len())?;
            }
        }
        let families = self
            .families
            .into_iter()
            .map(|f| {
                rows(f.prototype.len())?;
                Ok(CycleFamilyState {
                    weight: index(f.weight)?,
                    root: index(f.root)?,
                    p: index(f.p)?,
                    q: index(f.q)?,
                    x: f.x.map(index).transpose()?,
                    prototype: f.prototype,
                    relevant: f.relevant,
                })
            })
            .collect::<Result<_, PickleError>>()?;
        Ok(ComponentState {
            graph: self.graph.into_model()?,
            original_nodes: self
                .original_nodes
                .into_iter()
                .map(index)
                .collect::<Result<_, _>>()?,
            original_edges: self
                .original_edges
                .into_iter()
                .map(|v| index(v).map(PersistedRingEdgeId::new))
                .collect::<Result<_, _>>()?,
            predecessor: optional_rows(self.predecessor)?,
            distance: optional_rows(self.distance)?,
            reachable: self.reachable,
            directed_paths: self
                .directed_paths
                .into_iter()
                .map(index_rows)
                .collect::<Result<_, _>>()?,
            cycle_families: families,
            prototype_counts: self
                .prototype_counts
                .into_iter()
                .map(index)
                .collect::<Result<_, _>>()?,
            relations: self.relations,
            urfs: index_rows(self.urfs)?,
        })
    }
}

impl From<&PersistedRingDecomposition> for RdlRecord {
    fn from(d: &PersistedRingDecomposition) -> Self {
        Self {
            graph: d.graph().into(),
            family_nodes: d
                .urfs()
                .iter()
                .map(|u| u.nodes().iter().map(|&v| v as u64).collect())
                .collect(),
            family_edges: d
                .urfs()
                .iter()
                .map(|u| u.edges().iter().map(|v| v.index() as u64).collect())
                .collect(),
            relevant_cycle_count_bits: d.relevant_cycle_count().to_bits(),
            components: d.persisted_components().map(Into::into).collect(),
        }
    }
}

impl RdlRecord {
    fn into_model(self, atoms: usize) -> Result<PersistedRingDecomposition, PickleError> {
        if index(self.graph.nodes)? > atoms.max(1) {
            return Err(invalid("RDL graph extent exceeds topology"));
        }
        rows(self.components.len())?;
        let nodes = index_rows(self.family_nodes)?;
        let edges = index_rows(self.family_edges)?;
        if nodes.len() != edges.len() {
            return Err(invalid("RDL family table length mismatch"));
        }
        let urfs = nodes
            .into_iter()
            .zip(edges)
            .map(|(nodes, edges)| {
                (
                    nodes,
                    edges.into_iter().map(PersistedRingEdgeId::new).collect(),
                )
            })
            .collect();
        PersistedRingDecomposition::from_persisted_state(
            self.graph.into_model()?,
            urfs,
            f64::from_bits(self.relevant_cycle_count_bits),
            self.components
                .into_iter()
                .map(RdlComponent::into_model)
                .collect::<Result<_, _>>()?,
        )
        .map_err(invalid)
    }
}
#[derive(Debug, Clone, Encode, Decode)]
#[musli(Binary, name(type = u32))]
struct ValenceRecord {
    #[musli(Binary, name = 0)]
    explicit: Vec<i32>,
    #[musli(Binary, name = 1)]
    implicit: Vec<i32>,
}
#[derive(Debug, Clone, Encode, Decode)]
#[musli(Binary, name(type = u32))]
struct DerivedState {
    #[musli(Binary, name = 0)]
    valid_bits: u16,
    #[musli(Binary, name = 1)]
    rings: Option<RingRecord>,
    #[musli(Binary, name = 2)]
    ring_families: Option<RingRecord>,
    #[musli(Binary, name = 3)]
    valence: Option<ValenceRecord>,
}

fn rows(count: usize) -> Result<(), PickleError> {
    if count > MAX_ROWS {
        return Err(PickleError::InvalidArchive(
            "archive 2 row limit exceeded".into(),
        ));
    }
    Ok(())
}
fn text(bytes: Vec<u8>) -> Result<String, PickleError> {
    String::from_utf8(bytes).map_err(invalid)
}
fn index(value: u64) -> Result<usize, PickleError> {
    usize::try_from(value).map_err(invalid)
}
fn ids<T: Copy>(values: &[T], f: impl Fn(T) -> usize) -> Vec<u64> {
    values.iter().map(|v| f(*v) as u64).collect()
}
fn atom_ids(values: Vec<u64>) -> Result<Vec<AtomId>, PickleError> {
    values
        .into_iter()
        .map(|v| index(v).map(AtomId::new))
        .collect()
}
fn bond_ids(values: Vec<u64>) -> Result<Vec<BondId>, PickleError> {
    values
        .into_iter()
        .map(|v| index(v).map(BondId::new))
        .collect()
}
fn props<'a>(values: impl Iterator<Item = (&'a PropertyText, &'a PropertyValue)>) -> Vec<Property> {
    values
        .map(|(key, value)| Property {
            key: key.as_bytes().to_vec(),
            value: value.into(),
            computed: false,
        })
        .collect()
}
fn checked_props(
    values: Vec<Property>,
) -> Result<Vec<(PropertyText, PropertyValue, bool)>, PickleError> {
    let values = values
        .into_iter()
        .map(|p| {
            let key = PropertyText::from(p.key);
            Ok((key, p.value.into_model()?, p.computed))
        })
        .collect::<Result<Vec<_>, PickleError>>()?;
    let mut seen = BTreeSet::new();
    for (key, _, _) in &values {
        if !seen.insert(key) {
            return Err(PickleError::InvalidArchive("duplicate property key".into()));
        }
    }
    Ok(values)
}
fn checked_text_props(
    values: Vec<(Vec<u8>, Vec<u8>)>,
) -> Result<Vec<(Vec<u8>, Vec<u8>)>, PickleError> {
    let mut seen = BTreeSet::new();
    for (key, _) in &values {
        if !seen.insert(key) {
            return Err(PickleError::InvalidArchive(
                "duplicate text property key".into(),
            ));
        }
    }
    Ok(values)
}

impl From<&AtomPdbResidueInfo> for PdbInfo {
    fn from(p: &AtomPdbResidueInfo) -> Self {
        Self {
            atom_name: p.atom_name().as_bytes().to_vec(),
            serial_number: p.serial_number(),
            alt_loc: p.alt_loc().as_bytes().to_vec(),
            residue_name: p.residue_name().as_bytes().to_vec(),
            residue_number: p.residue_number(),
            chain_id: p.chain_id().as_bytes().to_vec(),
            insertion_code: p.insertion_code().as_bytes().to_vec(),
            occupancy_bits: p.occupancy().to_bits(),
            temp_factor_bits: p.temp_factor().to_bits(),
            hetero: p.is_hetero_atom(),
            secondary_structure: p.secondary_structure(),
            segment_number: p.segment_number(),
            monomer_class: p.monomer_class().as_bytes().to_vec(),
        }
    }
}
impl PdbInfo {
    fn into_model(self) -> AtomPdbResidueInfo {
        AtomPdbResidueInfo::new(
            self.atom_name,
            self.serial_number,
            self.residue_name,
            self.residue_number,
            self.chain_id,
            self.hetero,
        )
        .with_alt_loc(self.alt_loc)
        .with_insertion_code(self.insertion_code)
        .with_occupancy(f64::from_bits(self.occupancy_bits))
        .with_temp_factor(f64::from_bits(self.temp_factor_bits))
        .with_secondary_structure(self.secondary_structure)
        .with_segment_number(self.segment_number)
        .with_monomer_class(self.monomer_class)
    }
}
impl From<&Atom> for AtomRecord {
    fn from(a: &Atom) -> Self {
        Self {
            atomic_number: a.atomic_number(),
            formal_charge: a.formal_charge(),
            isotope: a.isotope(),
            chiral: a.chiral_tag().into(),
            chiral_permutation: a.chiral_permutation(),
            unknown_stereo: a.unknown_stereo(),
            mol_parity: a.mol_parity(),
            mol_inversion: a.mol_inversion_flag(),
            radical_electrons: a.radical_electrons(),
            aromatic: a.is_aromatic(),
            hybridization: a.hybridization().into(),
            atom_map: a.atom_map(),
            no_implicit: a.no_implicit(),
            implicit_hydrogen: a.implicit_hydrogen(),
            explicit_hydrogens: a.explicit_hydrogens(),
            tracked_isotopes: a.tracked_isotopic_hydrogens().to_vec(),
            properties: props(ordered_atom_properties(a)),
            temporary_flags: a.temporary_flags(),
            pdb: a.pdb_residue_info().map(Into::into),
            attachment_order: a.template_attachment_order().map(|order| {
                order
                    .entries()
                    .iter()
                    .map(|a| Attachment {
                        target: a.target().index() as u64,
                        label: a.label().into(),
                    })
                    .collect()
            }),
        }
    }
}
impl AtomRecord {
    fn into_model(self, id: usize) -> Result<Atom, PickleError> {
        let element = Element::from_atomic_number(self.atomic_number).ok_or(
            PickleError::InvalidEnumValue {
                value: self.atomic_number,
                type_name: "Element",
            },
        )?;
        let mut a = AtomSpec::new(element)
            .with_formal_charge(self.formal_charge)
            .with_chiral_tag(self.chiral.into())
            .with_unknown_stereo(self.unknown_stereo)
            .with_radical_electrons(self.radical_electrons)
            .with_aromatic(self.aromatic)
            .with_hybridization(self.hybridization.into())
            .with_no_implicit(self.no_implicit)
            .with_implicit_hydrogen(self.implicit_hydrogen)
            .with_explicit_hydrogens(self.explicit_hydrogens)
            .with_tracked_isotopic_hydrogens(self.tracked_isotopes);
        if let Some(v) = self.isotope {
            a = a.with_isotope(v);
        }
        if let Some(v) = self.chiral_permutation {
            a = a.with_chiral_permutation(v);
        }
        if let Some(v) = self.mol_parity {
            a = a.with_mol_parity(v);
        }
        if let Some(v) = self.mol_inversion {
            a = a.with_mol_inversion_flag(v);
        }
        if let Some(v) = self.atom_map {
            a = a.with_atom_map(v);
        }
        if let Some(v) = self.pdb {
            a = a.with_pdb_residue_info(v.into_model());
        }
        if let Some(values) = self.attachment_order {
            let entries = values
                .into_iter()
                .map(|v| {
                    Ok(TemplateAttachment::new(
                        AtomId::new(index(v.target)?),
                        v.label,
                    ))
                })
                .collect::<Result<_, PickleError>>()?;
            a = a.with_template_attachment_order(
                TemplateAttachmentOrder::new(entries).map_err(invalid)?,
            );
        }
        for (key, value, computed) in checked_props(self.properties)? {
            a = if computed {
                a.with_computed_prop(key, value)
            } else {
                a.with_prop(key, value)
            }
            .map_err(invalid)?;
        }
        let mut a = Atom::from_spec(AtomId::new(id), a);
        a.set_temporary_flags(self.temporary_flags);
        Ok(a)
    }
}
impl From<&Bond> for BondRecord {
    fn from(b: &Bond) -> Self {
        Self {
            begin: b.begin().index() as u64,
            end: b.end().index() as u64,
            order: b.order().into(),
            stereo: b.stereo().into(),
            direction: b.direction().into(),
            aromatic: b.is_aromatic(),
            conjugated: b.is_conjugated(),
            stereo_atoms: b.stereo_atoms().map(|v| v.map(|id| id.index() as u64)),
            unknown_stereo: b.unknown_stereo(),
            properties: props(ordered_bond_properties(b)),
            temporary_flags: b.temporary_flags(),
            source_stereo_atoms: (!b.stereo_atom_references().is_empty()
                && b.stereo_atoms().is_none())
            .then(|| ids(b.stereo_atom_references(), AtomId::index)),
        }
    }
}
impl BondRecord {
    fn into_model(self, id: usize) -> Result<Bond, PickleError> {
        let mut b = BondSpec::new(
            AtomId::new(index(self.begin)?),
            AtomId::new(index(self.end)?),
            self.order.into(),
        )
        .with_stereo(self.stereo.into())
        .with_direction(self.direction.into())
        .with_aromatic(self.aromatic)
        .with_conjugated(self.conjugated)
        .with_unknown_stereo(self.unknown_stereo);
        if let Some([x, y]) = self.stereo_atoms {
            b = b.with_stereo_atoms(AtomId::new(index(x)?), AtomId::new(index(y)?));
        }
        for (key, value, computed) in checked_props(self.properties)? {
            b = if computed {
                b.with_computed_prop(key, value)
            } else {
                b.with_prop(key, value)
            }
            .map_err(invalid)?;
        }
        let mut b = Bond::from_spec(BondId::new(id), b);
        b.set_temporary_flags(self.temporary_flags);
        if let Some(references) = self.source_stereo_atoms {
            if self.stereo_atoms.is_some() {
                return Err(invalid("conflicting bond stereo reference representations"));
            }
            b.set_source_stereo_atom_references(atom_ids(references)?);
        }
        Ok(b)
    }
}

impl From<&SubstanceGroup> for SGroup {
    fn from(g: &SubstanceGroup) -> Self {
        Self {
            id: g.id().index() as u64,
            sequence_id: g.rdkit_sequence_id(),
            external_id: g.external_id(),
            kind: g.kind().into(),
            atoms: ids(g.atoms(), AtomId::index),
            bonds: ids(g.bonds(), BondId::index),
            roles: g
                .stored_bond_roles()
                .iter()
                .map(|(b, role)| (b.index() as u64, (*role).into()))
                .collect(),
            parent_atoms: ids(g.parent_atoms(), AtomId::index),
            parent: g.parent().map(|p| p.index() as u64),
            label: g.label().map(|v| v.as_bytes().to_vec()),
            connection: g.connection().map(Into::into),
            subtype: g.subtype().map(|v| v.as_bytes().to_vec()),
            bracket_style: g.bracket_style().map(Into::into),
            expansion: g.expansion_state().map(|v| v.as_bytes().to_vec()),
            class: g.class().map(|v| v.as_bytes().to_vec()),
            component_number: g.component_number(),
            display: g.display().map(|d| GroupDisplay {
                brackets: d
                    .brackets
                    .iter()
                    .map(|b| b.points.map(|p| p.map(f64::to_bits)))
                    .collect(),
                field_position: d.field_position.map(|p| p.map(f64::to_bits)),
                tag: d.display_tag.as_ref().map(|v| v.as_bytes().to_vec()),
            }),
            data: g.data().map(|d| GroupData {
                field_name: d.field_name.as_ref().map(|v| v.as_bytes().to_vec()),
                field_type: d.field_type.as_ref().map(|v| v.as_bytes().to_vec()),
                field_info: d.field_info.as_ref().map(|v| v.as_bytes().to_vec()),
                field_display: d.field_display.as_ref().map(|v| v.as_bytes().to_vec()),
                units: d.units.as_ref().map(|v| v.as_bytes().to_vec()),
                query_type: d.query_type.as_ref().map(|v| v.as_bytes().to_vec()),
                query_op: d.query_op.as_ref().map(|v| v.as_bytes().to_vec()),
                values: d.values.iter().map(|v| v.as_bytes().to_vec()).collect(),
            }),
            attachments: g
                .attach_points()
                .iter()
                .map(|a| AttachPoint {
                    atom: a.atom.index() as u64,
                    leaving: a.leaving_atom.map(|v| v.index() as u64),
                    label: a.label.as_ref().map(|v| v.as_bytes().to_vec()),
                    order: a.order,
                })
                .collect(),
            cstates: g
                .cstates()
                .iter()
                .map(|c| CState {
                    bond: c.bond.index() as u64,
                    vector: c.vector.map(f64::to_bits),
                })
                .collect(),
            properties: g
                .property_records()
                .map(|(k, v)| (k.as_bytes().to_vec(), v.into()))
                .collect(),
            data_fields: g
                .data_fields()
                .iter()
                .map(|v| v.as_bytes().to_vec())
                .collect(),
            head_crossing: ids(g.head_crossing_bonds(), BondId::index),
            correspondence: ids(g.crossing_bond_correspondence(), BondId::index),
        }
    }
}
impl SGroup {
    fn into_model(self) -> Result<SubstanceGroup, PickleError> {
        let mut g = SubstanceGroup::new(SubstanceGroupId::new(index(self.id)?), self.kind.into())
            .with_atoms(atom_ids(self.atoms)?)
            .with_bonds(bond_ids(self.bonds)?)
            .with_parent_atoms(atom_ids(self.parent_atoms)?)
            .with_head_crossing_bonds(bond_ids(self.head_crossing)?)
            .with_crossing_bond_correspondence(bond_ids(self.correspondence)?);
        if let Some(v) = self.sequence_id {
            g = g.with_rdkit_sequence_id(v);
        }
        if let Some(v) = self.external_id {
            g = g.with_external_id(v);
        }
        if let Some(v) = self.parent {
            g = g.with_parent(SubstanceGroupId::new(index(v)?));
        }
        if let Some(v) = self.label {
            g = g.with_label(v);
        }
        if let Some(v) = self.subtype {
            g = g.with_subtype(v);
        }
        if let Some(v) = self.expansion {
            g = g.with_expansion_state(v);
        }
        if let Some(v) = self.class {
            g = g.with_class(v);
        }
        if let Some(v) = self.component_number {
            g = g.with_component_number(v);
        }
        if let Some(value) = self.connection {
            g = g.with_connection(value.into());
        }
        if let Some(value) = self.bracket_style {
            g = g.with_bracket_style(value.into());
        }
        if let Some(d) = self.display {
            g = g.with_display(SGroupDisplay {
                brackets: d
                    .brackets
                    .into_iter()
                    .map(|p| SGroupBracket::new(p.map(|p| p.map(f64::from_bits))))
                    .collect(),
                field_position: d.field_position.map(|p| p.map(f64::from_bits)),
                display_tag: d.tag.map(Into::into),
            });
        }
        if let Some(d) = self.data {
            g = g.with_data(SGroupData {
                field_name: d.field_name.map(Into::into),
                field_type: d.field_type.map(Into::into),
                field_info: d.field_info.map(Into::into),
                field_display: d.field_display.map(Into::into),
                units: d.units.map(Into::into),
                query_type: d.query_type.map(Into::into),
                query_op: d.query_op.map(Into::into),
                values: d.values.into_iter().map(Into::into).collect(),
            });
        }

        g = g.with_attach_points(
            self.attachments
                .into_iter()
                .map(|a| {
                    Ok(SGroupAttachPoint {
                        atom: AtomId::new(index(a.atom)?),
                        leaving_atom: a.leaving.map(index).transpose()?.map(AtomId::new),
                        label: a.label.map(Into::into),
                        order: a.order,
                    })
                })
                .collect::<Result<_, PickleError>>()?,
        );
        g = g.with_cstates(
            self.cstates
                .into_iter()
                .map(|c| {
                    Ok(SGroupCState::new(
                        BondId::new(index(c.bond)?),
                        c.vector.map(f64::from_bits),
                    ))
                })
                .collect::<Result<_, PickleError>>()?,
        );

        let mut seen = BTreeSet::new();
        for (bond, tag) in self.roles {
            let bond = BondId::new(index(bond)?);
            if !seen.insert(bond) || !g.bonds().contains(&bond) {
                return Err(PickleError::InvalidArchive(
                    "invalid or duplicate SGroup role".into(),
                ));
            }
            g = g.with_bond_role(bond, tag.into());
        }
        for (k, v) in self.properties {
            g = g
                .with_prop(PropertyText::from(k), v.into_model()?)
                .map_err(invalid)?;
        }
        for v in self.data_fields {
            g = g.with_data_field(v);
        }
        Ok(g)
    }
}

impl From<&RingInfo> for RingRecord {
    fn from(r: &RingInfo) -> Self {
        Self {
            initialized: r.is_initialized(),
            find_type: r.persisted_find_type().into(),
            atom_extent: r.atom_row_count() as u64,
            bond_extent: r.bond_row_count() as u64,
            atom_rings: r
                .atom_rings()
                .iter()
                .map(|v| ids(v, AtomId::index))
                .collect(),
            bond_rings: r
                .bond_rings()
                .iter()
                .map(|v| ids(v, BondId::index))
                .collect(),
            atom_families: r
                .atom_ring_families()
                .iter()
                .map(|v| ids(v, AtomId::index))
                .collect(),
            bond_families: r
                .bond_ring_families()
                .iter()
                .map(|v| ids(v, BondId::index))
                .collect(),
            relevant_cycles: r.persisted_relevant_cycle_count().map(|v| v as u64),
            fused_rings: r.persisted_fused_rings().to_vec(),
            fused_bonds: r
                .persisted_num_fused_bonds()
                .iter()
                .map(|v| *v as u64)
                .collect(),
            decomposition: r.persisted_decomposition().map(Into::into),
            atom_members: Some(wire_rows(r.persisted_atom_members())),
            bond_members: Some(wire_rows(r.persisted_bond_members())),
        }
    }
}
impl RingRecord {
    fn into_model(self, atoms: usize, bonds: usize) -> Result<RingInfo, PickleError> {
        let ae = index(self.atom_extent)?;
        let be = index(self.bond_extent)?;
        if ae > atoms || be > bonds {
            return Err(PickleError::InvalidArchive(
                "ring cache extents exceed topology".into(),
            ));
        }
        for count in [
            self.atom_rings.len(),
            self.bond_rings.len(),
            self.atom_families.len(),
            self.bond_families.len(),
            self.fused_rings.len(),
            self.fused_bonds.len(),
        ] {
            rows(count)?;
        }
        for row in &self.fused_rings {
            rows(row.len())?;
        }
        let mut info = RingInfo::from_persisted_components(
            self.initialized,
            self.find_type.into(),
            ae,
            be,
            self.atom_rings
                .into_iter()
                .map(atom_ids)
                .collect::<Result<_, _>>()?,
            self.bond_rings
                .into_iter()
                .map(bond_ids)
                .collect::<Result<_, _>>()?,
            self.atom_families
                .into_iter()
                .map(atom_ids)
                .collect::<Result<_, _>>()?,
            self.bond_families
                .into_iter()
                .map(bond_ids)
                .collect::<Result<_, _>>()?,
            self.relevant_cycles.map(index).transpose()?,
            self.fused_rings,
            self.fused_bonds
                .into_iter()
                .map(index)
                .collect::<Result<_, _>>()?,
        )
        .map_err(invalid)?;
        info.restore_persisted_state(
            self.decomposition
                // Source addHs preserves RingInfo/RDL while appending atoms.
                // Membership rows and RDL have independently stored extents.
                // AddHs can subsequently enlarge the molecule without
                // enlarging RDL. Bound its graph by the
                // live topology here; restore_persisted_state separately
                // checks that it covers the saved membership rows.
                .map(|d| d.into_model(atoms))
                .transpose()?,
            self.atom_members.map(index_rows).transpose()?,
            self.bond_members.map(index_rows).transpose()?,
        )
        .map_err(invalid)?;
        Ok(info)
    }
}

impl MoleculeState {
    fn from_input(m: &BinaryInput<'_>) -> Result<Self, PickleError> {
        // A concrete molecule archive is not a QueryGraph archive. Do not
        // silently drop the query flag as the historical raw writer did.
        if m.topology.bonds.iter().any(|b| b.query().is_some()) {
            return Err(PickleError::InvalidMolecule(
                "query state in concrete molecule archive".into(),
            ));
        }
        Ok(Self {
            atoms: m.atoms().iter().map(Into::into).collect(),
            bonds: m.bonds().iter().map(Into::into).collect(),
            conformers_2d: m
                .coordinates
                .conformers_2d
                .iter()
                .map(|c| Conformer2 {
                    id: c.id() as u64,
                    coordinates: c
                        .coordinates()
                        .iter()
                        .map(|p| p.map(f64::to_bits))
                        .collect(),
                    properties: c
                        .props()
                        .iter()
                        .map(|(k, v)| (k.as_bytes().to_vec(), v.as_bytes().to_vec()))
                        .collect(),
                })
                .collect(),
            conformers_3d: m
                .coordinates
                .conformers_3d
                .iter()
                .map(|c| Conformer3 {
                    id: c.id() as u64,
                    coordinates: c
                        .coordinates()
                        .iter()
                        .map(|p| p.map(f64::to_bits))
                        .collect(),
                    is_3d: c.is_3d(),
                    properties: c
                        .props()
                        .iter()
                        .map(|(k, v)| (k.as_bytes().to_vec(), v.as_bytes().to_vec()))
                        .collect(),
                })
                .collect(),
            source_dimension: m.coordinates.source_coordinate_dim.map(Into::into),
            source_order: m
                .coordinates
                .source_conformer_order
                .as_ref()
                .map(|v| v.iter().copied().map(Into::into).collect()),
            sgroups: m.topology.substance_groups.iter().map(Into::into).collect(),
            stereo_groups: m
                .topology
                .stereo_groups
                .iter()
                .map(|s| StereoRecord {
                    id: s.id(),
                    write_id: s.write_id(),
                    kind: s.kind().into(),
                    atoms: ids(s.atoms(), AtomId::index),
                    bonds: ids(s.bonds(), BondId::index),
                })
                .collect(),
            name: m.properties.name().map(|v| v.as_bytes().to_vec()),
            properties: props(m.properties.ordered_props()),
            sdf_fields: m
                .properties
                .sdf_data_fields()
                .iter()
                .map(|(k, v)| (k.as_bytes().to_vec(), v.as_bytes().to_vec()))
                .collect(),
            sdf_lists: m
                .properties
                .sdf_property_lists()
                .iter()
                .map(|p| PropertyList {
                    target: p.target().into(),
                    name: p.name().as_bytes().to_vec(),
                    values: p
                        .values()
                        .iter()
                        .map(|v| v.as_ref().map(Into::into))
                        .collect(),
                })
                .collect(),
        })
    }
    fn into_record(self, derived: DerivedState) -> Result<BinaryRecord, PickleError> {
        if self.atoms.len() > MAX_ROWS {
            return Err(PickleError::TooManyAtoms(self.atoms.len()));
        }
        if self.bonds.len() > MAX_ROWS {
            return Err(PickleError::TooManyBonds(self.bonds.len()));
        }
        let atoms = self
            .atoms
            .into_iter()
            .enumerate()
            .map(|(i, a)| a.into_model(i))
            .collect::<Result<Vec<_>, _>>()?;
        let bonds = self
            .bonds
            .into_iter()
            .enumerate()
            .map(|(i, b)| b.into_model(i))
            .collect::<Result<Vec<_>, _>>()?;
        let sgroups = self
            .sgroups
            .into_iter()
            .map(SGroup::into_model)
            .collect::<Result<_, _>>()?;
        let stereo_groups = self
            .stereo_groups
            .into_iter()
            .map(|s| {
                let mut group =
                    StereoGroup::new(s.kind.into(), atom_ids(s.atoms)?, bond_ids(s.bonds)?)?
                        .with_write_id(s.write_id);
                if let Some(id) = s.id {
                    group = group.with_id(id).with_write_id(s.write_id);
                }
                Ok(group)
            })
            .collect::<Result<_, PickleError>>()?;
        let topology =
            TopologyBlock::try_from_parts(atoms, bonds, sgroups, stereo_groups).map_err(invalid)?;
        let atom_count = topology.atoms.len();
        let mut two = Vec::new();
        let mut three = Vec::new();
        for c in self.conformers_2d {
            if c.coordinates.len() != atom_count {
                return Err(PickleError::DataLengthMismatch {
                    expected: atom_count,
                    actual: c.coordinates.len(),
                });
            }
            let mut conf = Conformer2D::new(
                index(c.id)?,
                c.coordinates
                    .into_iter()
                    .map(|p| p.map(f64::from_bits))
                    .collect(),
            );
            for (k, v) in checked_text_props(c.properties)? {
                conf = conf.with_prop(k, v);
            }
            two.push(conf);
        }
        for c in self.conformers_3d {
            if c.coordinates.len() != atom_count {
                return Err(PickleError::DataLengthMismatch {
                    expected: atom_count,
                    actual: c.coordinates.len(),
                });
            }
            let mut conf = Conformer3D::new(
                index(c.id)?,
                c.coordinates
                    .into_iter()
                    .map(|p| p.map(f64::from_bits))
                    .collect(),
                c.is_3d,
            );
            for (k, v) in checked_text_props(c.properties)? {
                conf = conf.with_prop(k, v);
            }
            three.push(conf);
        }
        let coordinates = CoordinateBlock {
            conformers_2d: two,
            conformers_3d: three,
            source_coordinate_dim: self.source_dimension.map(Into::into),
            source_conformer_order: self
                .source_order
                .map(|v| v.into_iter().map(Into::into).collect()),
        };
        let mut properties = MoleculeProperties::default();
        if let Some(v) = self.name {
            properties = properties.with_name(v);
        }
        let mut seen = BTreeSet::new();
        for p in &self.properties {
            if !seen.insert(&p.key) {
                return Err(PickleError::InvalidArchive(
                    "duplicate molecule property".into(),
                ));
            }
        }
        for p in self.properties {
            properties = if p.computed {
                properties.with_computed_prop(PropertyText::from(p.key), p.value.into_model()?)
            } else {
                properties.with_prop(PropertyText::from(p.key), p.value.into_model()?)
            }
            .map_err(invalid)?;
        }
        for (k, v) in self.sdf_fields {
            properties = properties.with_sdf_data_field(k, v);
        }
        for p in self.sdf_lists {
            let target = p.target.into();
            properties = properties.with_sdf_property_list(SdfPropertyList::new(
                target,
                p.name,
                p.values
                    .into_iter()
                    .map(|v| v.map(Value::into_model).transpose())
                    .collect::<Result<_, _>>()?,
            ));
        }
        let rings = derived
            .rings
            .map(|r| r.into_model(atom_count, topology.bonds.len()))
            .transpose()?;
        let ring_families = derived
            .ring_families
            .map(|r| r.into_model(atom_count, topology.bonds.len()))
            .transpose()?;
        let valence = derived
            .valence
            .map(|v| {
                if v.explicit.len() != atom_count || v.implicit.len() != atom_count {
                    return Err(PickleError::InvalidArchive(
                        "valence cache row count differs from topology".into(),
                    ));
                }
                Ok(ValenceAssignment {
                    explicit_valence: v.explicit,
                    implicit_hydrogens: v.implicit,
                })
            })
            .transpose()?;
        let record = BinaryRecord {
            topology,
            coordinates,
            properties,
            derived: BinaryDerivedState {
                rings,
                ring_families,
                valence,
                aromaticity_valid: derived.valid_bits & 8 != 0,
                stereo_valid: derived.valid_bits & 16 != 0,
                valid_bits: Some(derived.valid_bits),
            },
        };
        record.validate()?;
        Ok(record)
    }
}

pub(super) fn encode(input: &BinaryInput<'_>) -> Result<Vec<u8>, PickleError> {
    let molecule = MoleculeState::from_input(input)?;
    let derived = DerivedState {
        valid_bits: input.derived.valid_bits,
        rings: input.derived.rings.map(Into::into),
        ring_families: input.derived.ring_families.map(Into::into),
        valence: input.derived.valence.map(|v| ValenceRecord {
            explicit: v.explicit_valence.clone(),
            implicit: v.implicit_hydrogens.clone(),
        }),
    };
    // Use the same structural admission on writes and reads. This does not
    // perceive chemistry, update caches, or acquire runtime commit authority.
    rows(input.num_atoms())?;
    rows(input.topology.bonds.len())?;
    if let Some(valence) = input.derived.valence {
        if valence.explicit_valence.len() != input.num_atoms()
            || valence.implicit_hydrogens.len() != input.num_atoms()
        {
            return Err(PickleError::InvalidArchive(
                "valence cache row count differs from topology".into(),
            ));
        }
    }
    // Encode by reference, then consume the DTO in the existing model
    // validator instead of deep-cloning its ring tables just for validation.
    // Invalid state still returns an error before any bytes leave this call.
    let derived_payload = musli::wire::to_vec(&derived);
    for ring in [derived.rings, derived.ring_families].into_iter().flatten() {
        ring.into_model(input.num_atoms(), input.topology.bonds.len())?;
    }
    let metadata = Metadata {
        producer: env!("CARGO_PKG_VERSION").into(),
        codec_contract: CODEC_CONTRACT.into(),
        molecule_schema: MOLECULE_SCHEMA,
        derived_schema: DERIVED_SCHEMA,
    };
    let mut output = MAGIC.to_vec();
    write_u16_le(&mut output, MAJOR);
    write_u16_le(&mut output, MINOR);
    write_u16_le(&mut output, 3);
    for (id, schema, payload) in [
        (
            SECTION_MANIFEST,
            METADATA_SCHEMA,
            musli::wire::to_vec(&metadata),
        ),
        (
            SECTION_MOLECULE_STATE,
            MOLECULE_SCHEMA,
            musli::wire::to_vec(&molecule),
        ),
        (SECTION_DERIVED_STATE, DERIVED_SCHEMA, derived_payload),
    ] {
        let payload =
            payload.map_err(|e| PickleError::InvalidArchive(format!("archive 2 encode: {e}")))?;
        write_archive_section(
            &mut output,
            id,
            schema,
            SECTION_FLAG_REQUIRED,
            CODEC,
            &payload,
        )?;
    }
    if output.len() > MAX_ARCHIVE_BYTES {
        return Err(PickleError::InvalidArchive(
            "archive 2 byte limit exceeded".into(),
        ));
    }
    Ok(output)
}

pub(super) fn decode(data: &[u8]) -> Result<BinaryRecord, PickleError> {
    if data.len() > MAX_ARCHIVE_BYTES {
        return Err(PickleError::InvalidArchive(
            "archive 2 byte limit exceeded".into(),
        ));
    }
    let (major, minor, sections) = read_archive_envelope(data, MAGIC)?;
    if (major, minor) != (MAJOR, MINOR) {
        return Err(PickleError::UnsupportedArchiveVersion { major, minor });
    }
    let mut seen = BTreeSet::new();
    let mut metadata = None;
    let mut molecule = None;
    let mut derived = None;
    let mut codec = None;
    for section in sections {
        if !seen.insert(section.id) {
            return Err(PickleError::DuplicateSection(section.id));
        }
        if ![
            SECTION_MANIFEST,
            SECTION_MOLECULE_STATE,
            SECTION_DERIVED_STATE,
        ]
        .contains(&section.id)
        {
            if section.is_required() {
                return Err(PickleError::UnknownRequiredSection(section.id));
            }
            continue;
        }
        if !section.is_required() {
            return Err(PickleError::InvalidArchive(
                "archive 2 required block flag absent".into(),
            ));
        }
        if !matches!(section.codec, STORAGE_CODEC | CODEC)
            || codec
                .replace(section.codec)
                .is_some_and(|c| c != section.codec)
        {
            return Err(PickleError::InvalidArchive(
                "archive 2 block codec mismatch".into(),
            ));
        }
        match (section.id, section.version) {
            (SECTION_MANIFEST, METADATA_SCHEMA) => {
                metadata = Some(decode_payload::<Metadata>(section.payload, section.codec)?)
            }
            (SECTION_MOLECULE_STATE, MOLECULE_SCHEMA) => {
                molecule = Some(decode_payload::<MoleculeState>(
                    section.payload,
                    section.codec,
                )?)
            }
            (SECTION_DERIVED_STATE, DERIVED_SCHEMA) => {
                derived = Some(decode_payload::<DerivedState>(
                    section.payload,
                    section.codec,
                )?)
            }
            (section, version) => {
                return Err(PickleError::UnsupportedSectionVersion { section, version });
            }
        }
    }
    let metadata = metadata.ok_or(PickleError::MissingRequiredSection(SECTION_MANIFEST))?;
    let expected_contract = if codec == Some(STORAGE_CODEC) {
        STORAGE_CODEC_CONTRACT
    } else {
        CODEC_CONTRACT
    };
    if metadata.codec_contract != expected_contract
        || metadata.molecule_schema != MOLECULE_SCHEMA
        || metadata.derived_schema != DERIVED_SCHEMA
    {
        return Err(PickleError::InvalidArchive(
            "archive 2 metadata disagrees with block contract".into(),
        ));
    }
    molecule
        .ok_or(PickleError::MissingRequiredSection(SECTION_MOLECULE_STATE))?
        .into_record(derived.ok_or(PickleError::MissingRequiredSection(SECTION_DERIVED_STATE))?)
}
fn decode_payload<'de, T: Decode<'de, Binary, musli::alloc::Global>>(
    data: &'de [u8],
    codec: u8,
) -> Result<T, PickleError> {
    match codec {
        STORAGE_CODEC => musli::storage::from_slice(data)
            .map_err(|e| PickleError::InvalidArchive(format!("archive 2 storage decode: {e}"))),
        CODEC => musli::wire::from_slice(data)
            .map_err(|e| PickleError::InvalidArchive(format!("archive 2 wire decode: {e}"))),
        _ => Err(PickleError::InvalidArchive(
            "archive 2 unknown codec".into(),
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rdl_fixture() -> (TopologyBlock, RingInfo) {
        // Three independent paths yield two basis cycles and three relevant
        // cycles, so restoring only the public family rows cannot pass.
        let edges = [(0, 2), (2, 1), (0, 3), (3, 1), (0, 4), (4, 1)];
        let topology = TopologyBlock::try_from_parts(
            (0..5)
                .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
                .collect(),
            edges
                .into_iter()
                .enumerate()
                .map(|(i, (a, b))| {
                    Bond::from_spec(
                        BondId::new(i),
                        BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
                    )
                })
                .collect(),
            vec![],
            vec![],
        )
        .unwrap();
        let ring = cosmolkit_core::find_ring_families(&topology, &Default::default()).unwrap();
        (topology, ring)
    }

    #[test]
    fn archive20_rdl_full_state_roundtrip_and_legacy_missing_fields() {
        for (nodes, edges) in [
            (1, vec![]),
            (3, vec![(0, 1), (1, 2)]),
            (3, vec![(0, 1), (1, 2), (2, 0)]),
            (6, vec![(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3)]),
        ] {
            let mut graph = PersistedRingGraph::new(nodes);
            for (from, to) in edges {
                graph.add_undirected_edge(from, to).unwrap();
            }
            let original = PersistedRingDecomposition::calculate(graph).unwrap();
            let record = RdlRecord::from(&original);
            let restored = record.into_model(nodes).unwrap();
            assert_eq!(restored, original);
            assert_eq!(
                restored.smallest_ring_basis(),
                original.smallest_ring_basis()
            );
            assert_eq!(
                restored.relevant_cycles().collect::<Vec<_>>(),
                original.relevant_cycles().collect::<Vec<_>>()
            );
        }
        let (topology, original) = rdl_fixture();
        let payload = musli::wire::to_vec(&RingRecord::from(&original)).unwrap();
        let record: RingRecord = decode_payload(&payload, CODEC).unwrap();
        let restored = record.into_model(5, 6).unwrap();
        assert_eq!(restored, original);
        let before = original.persisted_decomposition().unwrap();
        let after = restored.persisted_decomposition().unwrap();
        assert_eq!(after.smallest_ring_basis(), before.smallest_ring_basis());
        assert_eq!(
            after.relevant_cycles().collect::<Vec<_>>(),
            before.relevant_cycles().collect::<Vec<_>>()
        );
        assert_eq!(
            after.relevant_cycle_count().to_bits(),
            before.relevant_cycle_count().to_bits()
        );
        assert_eq!(
            musli::wire::to_vec(&RingRecord::from(&restored)).unwrap(),
            payload
        );

        // Write the real pre-extension field set: no field 11, 12 or 13.
        #[derive(Encode)]
        #[musli(Binary, name(type = u32))]
        struct LegacyRing<'a> {
            #[musli(Binary, name = 0)]
            initialized: bool,
            #[musli(Binary, name = 1)]
            find_type: &'a RingFindRecord,
            #[musli(Binary, name = 2)]
            atom_extent: u64,
            #[musli(Binary, name = 3)]
            bond_extent: u64,
            #[musli(Binary, name = 4)]
            atom_rings: &'a [Vec<u64>],
            #[musli(Binary, name = 5)]
            bond_rings: &'a [Vec<u64>],
            #[musli(Binary, name = 6)]
            atom_families: &'a [Vec<u64>],
            #[musli(Binary, name = 7)]
            bond_families: &'a [Vec<u64>],
            #[musli(Binary, name = 8)]
            relevant_cycles: Option<u64>,
            #[musli(Binary, name = 9)]
            fused_rings: &'a [Vec<bool>],
            #[musli(Binary, name = 10)]
            fused_bonds: &'a [u64],
        }
        let record = RingRecord::from(&original);
        let legacy = LegacyRing {
            initialized: record.initialized,
            find_type: &record.find_type,
            atom_extent: record.atom_extent,
            bond_extent: record.bond_extent,
            atom_rings: &record.atom_rings,
            bond_rings: &record.bond_rings,
            atom_families: &record.atom_families,
            bond_families: &record.bond_families,
            relevant_cycles: record.relevant_cycles,
            fused_rings: &record.fused_rings,
            fused_bonds: &record.fused_bonds,
        };
        for codec in [STORAGE_CODEC, CODEC] {
            let bytes = if codec == STORAGE_CODEC {
                musli::storage::to_vec(&legacy).unwrap()
            } else {
                musli::wire::to_vec(&legacy).unwrap()
            };
            let decoded: RingRecord = decode_payload(&bytes, codec).unwrap();
            let mut old = decoded.into_model(5, 6).unwrap();
            assert!(!old.are_ring_families_initialized());
            assert_eq!(old.atom_ring_families(), original.atom_ring_families());
            assert_eq!(old.bond_ring_families(), original.bond_ring_families());
            assert_eq!(
                old.num_relevant_cycles().unwrap(),
                original.num_relevant_cycles().unwrap()
            );
            cosmolkit_core::symmetrized_sssr_with_existing(
                &topology,
                &mut old,
                false,
                &Default::default(),
            )
            .unwrap();
            assert!(old.are_ring_families_initialized());
            assert_eq!(old.atom_ring_families(), original.atom_ring_families());
            assert_eq!(
                old.persisted_decomposition(),
                original.persisted_decomposition()
            );
            assert_eq!(old.atom_rings().len(), 3);
        }
    }

    #[test]
    fn archive20_rdl_preserves_cache_extent_after_appending_atoms() {
        // AddHs.cpp::addHs: mol.clearComputedProps(false) retains RingInfo;
        // appended hydrogens do not expand the existing RDL graph.
        let (mut cyclic, _) = rdl_fixture();
        for i in 5..7 {
            cyclic
                .atoms
                .push(Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)));
        }
        let cyclic =
            TopologyBlock::try_from_parts(cyclic.atoms, cyclic.bonds, vec![], vec![]).unwrap();
        let acyclic = TopologyBlock::try_from_parts(
            (0..3)
                .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
                .collect(),
            vec![],
            vec![],
            vec![],
        )
        .unwrap();
        for topology in [cyclic, acyclic] {
            let original =
                cosmolkit_core::find_ring_families(&topology, &Default::default()).unwrap();
            assert_eq!(original.atom_row_count(), topology.atoms.len());
            assert_eq!(
                original
                    .persisted_decomposition()
                    .unwrap()
                    .graph()
                    .node_count(),
                topology.atoms.len()
            );
            let payload = musli::wire::to_vec(&RingRecord::from(&original)).unwrap();
            let record: RingRecord = decode_payload(&payload, CODEC).unwrap();
            assert_eq!(record.into_model(9, 10).unwrap(), original);
        }
        let (_, original) = rdl_fixture();
        let mut invalid_record = RingRecord::from(&original);
        invalid_record.decomposition.as_mut().unwrap().graph.nodes = 10;
        assert!(matches!(
            invalid_record.into_model(9, 10),
            Err(PickleError::InvalidMolecule(message))
                if message == "RDL graph extent exceeds topology"
        ));
        assert!(RingRecord::from(&original).into_model(4, 6).is_err());
    }

    #[test]
    fn archive20_rdl_rejects_invalid_internal_indices_and_paths() {
        let (_, original) = rdl_fixture();
        for mutation in 0..4 {
            let mut record = RingRecord::from(&original);
            let rdl = record.decomposition.as_mut().unwrap();
            match mutation {
                0 => rdl.components[0].original_edges[0] = 99,
                1 => rdl.components[0].directed_paths[0][0].push(99),
                2 => rdl.components[0].urfs[0][0] = 99,
                _ => rdl.relevant_cycle_count_bits = f64::NAN.to_bits(),
            }
            let expected = [
                "RDL component edge out of range",
                "invalid RDL shortest-path DAG edge",
                "RDL family index out of range",
                "invalid RDL graph or relevant-cycle count",
            ][mutation];
            assert!(
                matches!(record.into_model(5, 6), Err(PickleError::InvalidMolecule(message)) if message == expected),
                "mutation {mutation}"
            );
        }
    }

    #[test]
    fn archive20_complete_source_bond_reference_vectors_are_not_pair_truncated() {
        for references in [
            vec![],
            vec![AtomId::new(1)],
            vec![AtomId::new(0), AtomId::new(1)],
            vec![AtomId::new(1), AtomId::new(0), AtomId::new(1)],
        ] {
            let mut record = fixture();
            record.topology.bonds[0].set_source_stereo_atom_references(references.clone());
            record.topology.validate().unwrap();
            let bytes = encode_molecule_binary(&input(&record)).unwrap();
            let restored = decode_molecule_binary(&bytes).unwrap();
            assert_eq!(
                restored.topology.bonds[0].stereo_atom_references(),
                references
            );
            assert_eq!(restored.topology, record.topology);
            assert_eq!(encode_molecule_binary(&input(&restored)).unwrap(), bytes);
        }
        let record = fixture();
        let bytes = encode_molecule_binary(&input(&record)).unwrap();
        // Absent tag 11 exercises this schema's declared default.
        let restored = decode_molecule_binary(&bytes).unwrap();
        assert!(
            restored.topology.bonds[0]
                .stereo_atom_references()
                .is_empty()
        );
        let bad = mutate_molecule(&bytes, |m| m.bonds[0].source_stereo_atoms = Some(vec![2]));
        assert!(decode_molecule_binary(&bad).is_err());
    }

    #[test]
    fn archive20_pdb_text_preserves_bytes_and_utf8_sequence_encoding() {
        // Codec 2 encoded strings and byte vectors identically. Preserve this
        // old-reader guarantee separately from the new wire PDB round trip.
        for value in ["", " CA ", "é\0链"] {
            let text = musli::storage::to_vec(&value.to_owned()).unwrap();
            let bytes = musli::storage::to_vec(&value.as_bytes().to_vec()).unwrap();
            assert_eq!(text, bytes);
            assert_eq!(
                decode_payload::<Vec<u8>>(&text, STORAGE_CODEC).unwrap(),
                value.as_bytes()
            );
        }
        let mut record = fixture();
        let info = AtomPdbResidueInfo::new(
            vec![b'C', 0xff, 0],
            17,
            vec![0xfe, b'R'],
            3,
            vec![0xfd],
            true,
        )
        .with_alt_loc(vec![0xfc])
        .with_insertion_code(vec![0xfb])
        .with_monomer_class(vec![0xfa]);
        record.topology.atoms[0].set_pdb_residue_info(Some(info));
        let bytes = encode_molecule_binary(&input(&record)).unwrap();
        let restored = decode_molecule_binary(&bytes).unwrap();
        assert_eq!(restored.topology, record.topology);
        assert_eq!(encode_molecule_binary(&input(&restored)).unwrap(), bytes);
    }

    fn fixture() -> BinaryRecord {
        BinaryRecord {
            topology: TopologyBlock::try_from_parts(
                vec![
                    Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
                    Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::O)),
                ],
                vec![Bond::from_spec(
                    BondId::new(0),
                    BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
                )],
                vec![],
                vec![],
            )
            .unwrap(),
            coordinates: CoordinateBlock::default(),
            properties: MoleculeProperties::default(),
            derived: BinaryDerivedState::default(),
        }
    }
    fn input(r: &BinaryRecord) -> BinaryInput<'_> {
        BinaryInput {
            topology: &r.topology,
            coordinates: &r.coordinates,
            properties: &r.properties,
            derived: BinaryDerivedView::default(),
        }
    }
    fn pack(sections: &[ArchiveSection<'_>]) -> Vec<u8> {
        let mut data = MAGIC.to_vec();
        write_u16_le(&mut data, 2);
        write_u16_le(&mut data, 0);
        write_u16_le(&mut data, sections.len() as u16);
        for s in sections {
            write_archive_section(&mut data, s.id, s.version, s.flags, s.codec, s.payload).unwrap();
        }
        data
    }
    fn replace(data: &[u8], id: u16, payload: &[u8]) -> Vec<u8> {
        let (_, _, mut sections) = read_archive_envelope(data, MAGIC).unwrap();
        for s in &mut sections {
            if s.id == id {
                s.payload = payload;
            }
        }
        pack(&sections)
    }
    fn mutate_molecule(data: &[u8], f: impl FnOnce(&mut MoleculeState)) -> Vec<u8> {
        let (_, _, sections) = read_archive_envelope(data, MAGIC).unwrap();
        let mut state: MoleculeState =
            decode_payload(sections.iter().find(|s| s.id == 2).unwrap().payload, CODEC).unwrap();
        f(&mut state);
        replace(data, 2, &musli::wire::to_vec(&state).unwrap())
    }

    #[test]
    fn archive20_has_exactly_one_molecule_and_no_raw_or_canonical_companion() {
        let r = fixture();
        let data = encode_molecule_binary(&input(&r)).unwrap();
        assert_eq!(&data[..12], b"COSMOL\0\0\x02\0\0\0");
        let (major, minor, sections) = read_archive_envelope(&data, MAGIC).unwrap();
        assert_eq!((major, minor), (2, 0));
        assert_eq!(sections.iter().map(|s| s.id).collect::<Vec<_>>(), [1, 2, 3]);
        assert!(
            sections
                .iter()
                .all(|s| s.version == 1 && s.codec == 3 && s.is_required())
        );
        let metadata: Metadata = decode_payload(sections[0].payload, CODEC).unwrap();
        assert_eq!(metadata.codec_contract, "ck-wire-v1");
        assert_eq!((metadata.molecule_schema, metadata.derived_schema), (1, 1));
        let restored = decode_molecule_binary(&data).unwrap();
        assert_eq!(restored.topology, r.topology);
        assert_eq!(encode_molecule_binary(&input(&restored)).unwrap(), data);
        assert!(matches!(
            decode_sectioned_archive(&data),
            Err(PickleError::InvalidArchive(_))
        ));
        // Magic selects exactly one format family, not a fallback decoder.
        let mut legacy_magic = data.clone();
        legacy_magic[..8].copy_from_slice(ARCHIVE_MAGIC);
        assert_eq!(
            decode_molecule_binary(&legacy_magic).unwrap_err(),
            PickleError::UnsupportedArchiveVersion { major: 2, minor: 0 }
        );
        let mut legacy_version = data.clone();
        legacy_version[8..10].copy_from_slice(&1u16.to_le_bytes());
        assert_eq!(
            decode_molecule_binary(&legacy_version).unwrap_err(),
            PickleError::UnsupportedArchiveVersion { major: 1, minor: 0 }
        );
    }

    #[test]
    fn archive20_schema_dispatch_and_metadata_are_checked_per_block() {
        let data = encode_molecule_binary(&input(&fixture())).unwrap();
        let (_, _, sections) = read_archive_envelope(&data, MAGIC).unwrap();
        for id in [
            SECTION_MANIFEST,
            SECTION_MOLECULE_STATE,
            SECTION_DERIVED_STATE,
        ] {
            let mut unsupported = sections.clone();
            let block = unsupported.iter_mut().find(|s| s.id == id).unwrap();
            block.version = 2;
            block.payload = b"\xff";
            // An unsupported schema is identified before decoding its payload;
            // the other two block schemas remain at their own current versions.
            assert_eq!(
                decode(&pack(&unsupported)).unwrap_err(),
                PickleError::UnsupportedSectionVersion {
                    section: id,
                    version: 2
                }
            );
        }
        for which in 0..3 {
            let mut metadata: Metadata = decode_payload(sections[0].payload, CODEC).unwrap();
            match which {
                0 => metadata.molecule_schema = 2,
                1 => metadata.derived_schema = 2,
                _ => metadata.codec_contract = "unknown-ck-codec".into(),
            }
            assert_eq!(
                decode(&replace(
                    &data,
                    SECTION_MANIFEST,
                    &musli::wire::to_vec(&metadata).unwrap()
                ))
                .unwrap_err(),
                PickleError::InvalidArchive(
                    "archive 2 metadata disagrees with block contract".into()
                )
            );
        }
    }

    #[test]
    fn archive20_sparse_roles_and_reserved_property_states_are_lossless() {
        for explicit in [
            None,
            Some(SGroupBondRole::Crossing),
            Some(SGroupBondRole::Contained),
        ] {
            let mut r = fixture();
            let mut g =
                SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Superatom)
                    .with_atoms(vec![AtomId::new(0)])
                    .with_bonds(vec![BondId::new(0)]);
            if let Some(role) = explicit {
                g = g.with_bond_role(BondId::new(0), role);
            }
            r.topology.substance_groups.push(g);
            r.properties = r
                .properties
                .with_prop("__computedProps", "opaque\0value")
                .unwrap();
            r.topology.atoms[0]
                .set_prop("__computedProps", "opaque")
                .unwrap();
            assert!(r.topology.atoms[0].set_computed_prop("rank", 12).is_err());
            let data = encode_molecule_binary(&input(&r)).unwrap();
            let restored = decode_molecule_binary(&data).unwrap();
            assert_eq!(restored.topology, r.topology, "{explicit:?}");
            assert_eq!(restored.properties, r.properties);
        }
    }

    #[test]
    fn archive20_optional_unknown_blocks_skip_but_required_and_duplicates_fail() {
        let data = encode_molecule_binary(&input(&fixture())).unwrap();
        let (_, _, sections) = read_archive_envelope(&data, MAGIC).unwrap();
        let extra = ArchiveSection {
            id: 99,
            version: 999,
            flags: 0,
            codec: 255,
            payload: b"future",
        };
        let mut mutated = sections.clone();
        mutated.push(extra);
        assert!(decode(&pack(&mutated)).is_ok());
        mutated.last_mut().unwrap().flags = 1;
        assert_eq!(
            decode(&pack(&mutated)).unwrap_err(),
            PickleError::UnknownRequiredSection(99)
        );
        mutated = sections.clone();
        mutated.push(sections[0]);
        assert_eq!(
            decode(&pack(&mutated)).unwrap_err(),
            PickleError::DuplicateSection(1)
        );
        for id in 1..=3 {
            let mut missing = sections.clone();
            missing.retain(|s| s.id != id);
            assert_eq!(
                decode(&pack(&missing)).unwrap_err(),
                PickleError::MissingRequiredSection(id)
            );
            for which in 0..3 {
                let mut bad = sections.clone();
                let s = bad.iter_mut().find(|s| s.id == id).unwrap();
                match which {
                    0 => s.version = 2,
                    1 => s.codec = 1,
                    _ => s.flags = 0,
                }
                assert!(decode(&pack(&bad)).is_err(), "block {id}, mutation {which}");
            }
        }
    }

    #[test]
    fn archive20_truncation_trailing_payload_and_future_versions_rejected() {
        let data = encode_molecule_binary(&input(&fixture())).unwrap();
        for n in 0..data.len() {
            assert!(decode_molecule_binary(&data[..n]).is_err(), "truncated {n}");
        }
        let mut trailing = data.clone();
        trailing.push(0);
        assert!(decode(&trailing).is_err());
        let (_, _, sections) = read_archive_envelope(&data, MAGIC).unwrap();
        for s in sections {
            let mut payload = s.payload.to_vec();
            payload.push(0);
            assert!(
                decode(&replace(&data, s.id, &payload)).is_err(),
                "trailing block {}",
                s.id
            );
        }
        for (major, minor) in [(2_u16, 1_u16), (3, 0), (1, 3)] {
            let mut future = data.clone();
            future[8..10].copy_from_slice(&major.to_le_bytes());
            future[10..12].copy_from_slice(&minor.to_le_bytes());
            assert!(matches!(
                decode_molecule_binary(&future),
                Err(PickleError::UnsupportedArchiveVersion { .. })
            ));
        }
    }

    #[test]
    fn archive20_storage_golden_field_order_defaults_and_missing_fields() {
        let r = BinaryRecord {
            topology: TopologyBlock::default(),
            coordinates: CoordinateBlock::default(),
            properties: MoleculeProperties::default(),
            derived: BinaryDerivedState::default(),
        };
        let state = MoleculeState::from_input(&input(&r)).unwrap();
        // Frozen ck-storage-v1/molecule-schema-1 empty record, verified with musli 0.1.9.
        let golden = [
            12, 0, 0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8, 0, 9, 0, 10, 0, 11, 0,
        ];
        assert_eq!(musli::storage::to_vec(&state).unwrap(), golden);
        let mut reordered = vec![12];
        for pair in golden[1..].chunks_exact(2).rev() {
            reordered.extend_from_slice(pair);
        }
        assert!(decode_payload::<MoleculeState>(&reordered, STORAGE_CODEC).is_ok());
        // A missing optional field uses the declared schema default.
        let mut without_source_order = golden.to_vec();
        without_source_order[0] = 11;
        without_source_order.drain(11..13);
        assert!(
            decode_payload::<MoleculeState>(&without_source_order, STORAGE_CODEC)
                .unwrap()
                .source_order
                .is_none()
        );
        let mut missing = golden.to_vec();
        missing[0] = 11;
        missing.drain(1..3);
        assert!(decode_payload::<MoleculeState>(&missing, STORAGE_CODEC).is_err());
        let mut unknown = golden.to_vec();
        unknown[0] = 13;
        unknown.extend_from_slice(&[99, 0]);
        assert!(decode_payload::<MoleculeState>(&unknown, STORAGE_CODEC).is_err());
    }

    #[test]
    fn archive20_reads_frozen_storage_archive_and_rewrites_wire() {
        // Fixed codec-2/schema-1 payloads, not generated by the new writer.
        // Metadata: producer 0.5.0, ck-storage-v1, molecule/derived schema 1.
        let metadata = b"\x04\0\x050.5.0\x01\x0dck-storage-v1\x02\x01\x03\x01";
        let molecule = [
            12, 0, 0, 1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8, 0, 9, 0, 10, 0, 11, 0,
        ];
        let derived = [4, 0, 0, 1, 0, 2, 0, 3, 0];
        let sections = [
            ArchiveSection {
                id: 1,
                version: 1,
                flags: 1,
                codec: 2,
                payload: metadata,
            },
            ArchiveSection {
                id: 2,
                version: 1,
                flags: 1,
                codec: 2,
                payload: &molecule,
            },
            ArchiveSection {
                id: 3,
                version: 1,
                flags: 1,
                codec: 2,
                payload: &derived,
            },
        ];
        let old_bytes = pack(&sections);
        let restored = decode_molecule_binary(&old_bytes).unwrap();
        assert_eq!(restored.topology, TopologyBlock::default());
        assert_eq!(restored.coordinates, CoordinateBlock::default());
        assert_eq!(restored.properties, MoleculeProperties::default());
        assert_eq!(restored.derived.valid_bits, Some(0));
        let new_bytes = encode_molecule_binary(&input(&restored)).unwrap();
        let (_, _, new_sections) = read_archive_envelope(&new_bytes, MAGIC).unwrap();
        assert!(new_sections.iter().all(|s| s.codec == 3));
        let rewritten = decode_molecule_binary(&new_bytes).unwrap();
        assert_eq!(rewritten.topology, restored.topology);
        assert_eq!(rewritten.coordinates, restored.coordinates);
        assert_eq!(rewritten.properties, restored.properties);
        assert_eq!(rewritten.derived.valid_bits, restored.derived.valid_bits);
        assert!(rewritten.derived.rings.is_none());
        assert!(rewritten.derived.ring_families.is_none());
        assert!(rewritten.derived.valence.is_none());

        // Known codec IDs must agree across required blocks and with metadata.
        // Relabeling a payload must never cause a fallback to another decoder.
        for id in 1..=3 {
            let mut mislabeled = sections;
            mislabeled[(id - 1) as usize].codec = 3;
            assert!(decode(&pack(&mislabeled)).is_err());
        }
        let mut bad_contract: Metadata = decode_payload(metadata, STORAGE_CODEC).unwrap();
        bad_contract.codec_contract = CODEC_CONTRACT.into();
        let bad_payload = musli::storage::to_vec(&bad_contract).unwrap();
        assert!(matches!(
            decode(&replace(&old_bytes, 1, &bad_payload)),
            Err(PickleError::InvalidArchive(message))
                if message == "archive 2 metadata disagrees with block contract"
        ));
    }

    #[test]
    fn archive20_wire_golden_field_order_defaults_and_missing_fields() {
        // Wire uses a sequence of 24 map items, continuation-tagged numeric
        // keys and empty sequence values. IDs and Binary options are frozen.
        let golden = [
            0x98, 0xc0, 0x80, 0xc1, 0x80, 0xc2, 0x80, 0xc3, 0x80, 0xc4, 0x80, 0xc5, 0x80, 0xc6,
            0x80, 0xc7, 0x80, 0xc8, 0x80, 0xc9, 0x80, 0xca, 0x80, 0xcb, 0x80,
        ];
        let state: MoleculeState = decode_payload(&golden, CODEC).unwrap();
        assert!(state.atoms.is_empty());
        assert!(state.bonds.is_empty());
        assert_eq!(musli::wire::to_vec(&state).unwrap(), golden);
        let mut reordered = vec![0x98];
        for pair in golden[1..].chunks_exact(2).rev() {
            reordered.extend_from_slice(pair);
        }
        assert!(decode_payload::<MoleculeState>(&reordered, CODEC).is_ok());
        let mut without_source_order = golden.to_vec();
        without_source_order[0] = 0x96;
        without_source_order.drain(11..13);
        assert!(
            decode_payload::<MoleculeState>(&without_source_order, CODEC)
                .unwrap()
                .source_order
                .is_none()
        );
        let mut missing_required = golden.to_vec();
        missing_required[0] = 0x96;
        missing_required.drain(1..3);
        assert!(decode_payload::<MoleculeState>(&missing_required, CODEC).is_err());
    }

    #[test]
    fn archive20_wire_skips_unknown_fields_in_each_block_without_retaining_them() {
        let data = encode_molecule_binary(&input(&fixture())).unwrap();
        let expected = decode(&data).unwrap();
        let (_, _, sections) = read_archive_envelope(&data, MAGIC).unwrap();
        for section in &sections {
            // Extend the fixed schema's small map with a nested unknown field.
            // Müsli, not CK's decoder, must skip its complete value.
            let mut payload = section.payload.to_vec();
            assert!((0x80..0xbd).contains(&payload[0]));
            payload[0] += 2;
            payload.extend(musli::wire::to_vec(&99_u32).unwrap());
            payload
                .extend(musli::wire::to_vec(&vec![("future", vec![Some(1_u64), None])]).unwrap());
            let extended = replace(&data, section.id, &payload);
            let restored = decode(&extended).unwrap();
            assert_eq!(restored.topology, expected.topology, "block {}", section.id);
            assert_eq!(restored.coordinates, expected.coordinates);
            assert_eq!(restored.properties, expected.properties);
            assert_eq!(restored.derived.valid_bits, expected.derived.valid_bits);
            assert!(restored.derived.rings.is_none());
            assert!(restored.derived.ring_families.is_none());
            assert!(restored.derived.valence.is_none());
            assert_eq!(encode_molecule_binary(&input(&restored)).unwrap(), data);
            // Truncation inside an unknown value still returns a structured error.
            payload.pop();
            assert!(matches!(
                decode(&replace(&data, section.id, &payload)),
                Err(PickleError::InvalidArchive(_))
            ));
        }
    }

    #[test]
    fn archive20_invalid_indices_coordinates_and_duplicate_properties_rejected() {
        let data = encode_molecule_binary(&input(&fixture())).unwrap();
        for which in 0..6 {
            let bad = mutate_molecule(&data, |s| match which {
                0 => s.bonds[0].end = 99,
                1 => s.source_order = Some(vec![DimensionRecord::TwoD]),
                2 => s.conformers_2d.push(Conformer2 {
                    id: 0,
                    coordinates: vec![],
                    properties: vec![],
                }),
                3 => {
                    s.atoms[0].properties = vec![
                        Property {
                            key: vec![b'x'],
                            value: Value::Int(1),
                            computed: false
                        };
                        2
                    ]
                }
                4 => {
                    s.atoms[0].properties = vec![Property {
                        key: vec![],
                        value: Value::Bool(true),
                        computed: false,
                    }]
                }
                _ => {
                    s.properties = vec![
                        Property {
                            key: "duplicate".into(),
                            value: Value::String(b"value".to_vec()),
                            computed: false,
                        };
                        2
                    ]
                }
            });
            assert!(decode(&bad).is_err(), "mutation {which}");
        }
    }

    #[test]
    fn archive20_incomplete_sequences_are_rejected_by_musli() {
        // Müsli owns sequence framing, including a declared count without data.
        assert!(
            decode_payload::<MoleculeState>(&[12, 0, 0xc1, 0x84, 0x3d], STORAGE_CODEC).is_err()
        );
        assert!(decode_payload::<RingRecord>(&[1, 4, 1, 0xc1, 0x84, 0x3d], STORAGE_CODEC).is_err());
        assert!(decode_payload::<MoleculeState>(&[0x82, 0xc0, 0x81], CODEC).is_err());
        assert!(decode_payload::<RingRecord>(&[0x82, 0xc4, 0x81], CODEC).is_err());
    }

    #[test]
    fn archive20_cache_payload_validity_and_metadata_are_checked_independently() {
        let r = fixture();
        // Both cache slots must keep the same write admission after removing
        // the validation-only clone. Source caches remain borrowed throughout.
        for families in [false, true] {
            for atom_extent in [r.num_atoms(), r.num_atoms() + 1] {
                let ring = RingInfo::new(RingFindType::Fast, atom_extent, r.num_bonds());
                let mut borrowed = input(&r);
                if families {
                    borrowed.derived.ring_families = Some(&ring);
                } else {
                    borrowed.derived.rings = Some(&ring);
                }
                let result = encode_molecule_binary(&borrowed);
                if atom_extent == r.num_atoms() {
                    let bytes = result.unwrap();
                    let restored = decode(&bytes).unwrap();
                    let cache = if families {
                        restored.derived.ring_families.as_ref().unwrap()
                    } else {
                        restored.derived.rings.as_ref().unwrap()
                    };
                    assert_eq!(cache.atom_row_count(), atom_extent);
                    assert_eq!(cache.bond_row_count(), r.num_bonds());
                } else {
                    assert!(result.is_err());
                }
                assert_eq!(ring.atom_row_count(), atom_extent);
                assert_eq!(ring.bond_row_count(), r.num_bonds());
            }
        }
        for bits in [0, 4] {
            let valence = ValenceAssignment {
                explicit_valence: vec![1, 1],
                implicit_hydrogens: vec![3, 1],
            };
            let mut borrowed = input(&r);
            borrowed.derived.valence = Some(&valence);
            borrowed.derived.valid_bits = bits;
            let data = encode_molecule_binary(&borrowed).unwrap();
            let restored = decode(&data).unwrap();
            assert_eq!(restored.derived.valid_bits, Some(bits));
            assert_eq!(
                restored.derived.valence.as_ref().unwrap().explicit_valence,
                [1, 1]
            );
            assert_eq!(
                restored
                    .derived
                    .valence
                    .as_ref()
                    .unwrap()
                    .implicit_hydrogens,
                [3, 1]
            );
        }
        let data = encode_molecule_binary(&input(&r)).unwrap();
        let (_, _, sections) = read_archive_envelope(&data, MAGIC).unwrap();
        let original: DerivedState =
            decode_payload(sections.iter().find(|s| s.id == 3).unwrap().payload, CODEC).unwrap();
        for which in 0..3 {
            let mut bad = original.clone();
            match which {
                0 => bad.valid_bits = 0x100,
                1 => {
                    bad.valence = Some(ValenceRecord {
                        explicit: vec![1],
                        implicit: vec![1, 1],
                    })
                }
                2 => {
                    let mut ring = RingRecord::from(&RingInfo::new(RingFindType::Fast, 0, 0));
                    ring.atom_extent = 3;
                    bad.rings = Some(ring);
                }
                _ => unreachable!(),
            }
            assert!(
                decode(&replace(&data, 3, &musli::wire::to_vec(&bad).unwrap())).is_err(),
                "derived {which}"
            );
        }
        let mut metadata: Metadata =
            decode_payload(sections.iter().find(|s| s.id == 1).unwrap().payload, CODEC).unwrap();
        metadata.molecule_schema = 2;
        assert!(decode(&replace(&data, 1, &musli::wire::to_vec(&metadata).unwrap())).is_err());
    }

    #[test]
    fn archive20_enums_use_musli_fixed_tags_and_reject_unknown_variants() {
        fn verify<T>(values: &[(u8, T)], invalid: u8)
        where
            T: std::fmt::Debug
                + PartialEq
                + Encode<Binary>
                + for<'de> Decode<'de, Binary, musli::alloc::Global>,
        {
            for (tag, value) in values {
                let bytes = musli::storage::to_vec(value).unwrap();
                assert_eq!(bytes[0], *tag);
                assert_eq!(&musli::storage::from_slice::<T>(&bytes).unwrap(), value);
                let wire = musli::wire::to_vec(value).unwrap();
                // Derive encodes unit variants as their bare numeric tag;
                // payload variants encode a (tag, payload) sequence.
                if wire.len() == 1 {
                    assert_eq!(wire, [0xc0 | tag]);
                } else {
                    assert_eq!(&wire[..2], &[0x82, 0xc0 | tag]);
                }
                assert_eq!(&musli::wire::from_slice::<T>(&wire).unwrap(), value);
            }
            assert!(musli::storage::from_slice::<T>(&[invalid, 0]).is_err());
            assert!(musli::wire::from_slice::<T>(&[0x82, 0xc0 | invalid, 0x80]).is_err());
            assert!(musli::wire::from_slice::<T>(&[0xc0 | invalid]).is_err());
        }
        verify(
            &[
                (0, ChiralRecord::Unspecified),
                (1, ChiralRecord::TetrahedralCw),
                (2, ChiralRecord::TetrahedralCcw),
                (3, ChiralRecord::Other),
                (4, ChiralRecord::Tetrahedral),
                (5, ChiralRecord::Allene),
                (6, ChiralRecord::SquarePlanar),
                (7, ChiralRecord::TrigonalBipyramidal),
                (8, ChiralRecord::Octahedral),
            ],
            9,
        );
        verify(
            &[
                (0, HybridRecord::Unspecified),
                (1, HybridRecord::S),
                (2, HybridRecord::Sp),
                (3, HybridRecord::Sp2),
                (4, HybridRecord::Sp3),
                (5, HybridRecord::Sp2d),
                (6, HybridRecord::Sp3d),
                (7, HybridRecord::Sp3d2),
                (8, HybridRecord::Other),
            ],
            9,
        );
        verify(
            &[
                (0, BondOrderRecord::Unspecified),
                (1, BondOrderRecord::Single),
                (2, BondOrderRecord::Double),
                (3, BondOrderRecord::Triple),
                (4, BondOrderRecord::Quadruple),
                (5, BondOrderRecord::Quintuple),
                (6, BondOrderRecord::Hextuple),
                (7, BondOrderRecord::OneAndHalf),
                (8, BondOrderRecord::TwoAndHalf),
                (9, BondOrderRecord::ThreeAndHalf),
                (10, BondOrderRecord::FourAndHalf),
                (11, BondOrderRecord::FiveAndHalf),
                (12, BondOrderRecord::Aromatic),
                (13, BondOrderRecord::Ionic),
                (14, BondOrderRecord::Dative),
                (15, BondOrderRecord::DativeOne),
                (16, BondOrderRecord::DativeLeft),
                (17, BondOrderRecord::DativeRight),
                (18, BondOrderRecord::Hydrogen),
                (19, BondOrderRecord::ThreeCenter),
                (20, BondOrderRecord::Other),
                (21, BondOrderRecord::Zero),
            ],
            22,
        );
        verify(
            &[
                (0, BondDirectionRecord::None),
                (1, BondDirectionRecord::BeginWedge),
                (2, BondDirectionRecord::BeginDash),
                (3, BondDirectionRecord::EndUpRight),
                (4, BondDirectionRecord::EndDownRight),
                (5, BondDirectionRecord::EitherDouble),
                (6, BondDirectionRecord::Unknown),
            ],
            7,
        );
        verify(
            &[
                (0, BondStereoRecord::None),
                (1, BondStereoRecord::Any),
                (2, BondStereoRecord::Z),
                (3, BondStereoRecord::E),
                (4, BondStereoRecord::Cis),
                (5, BondStereoRecord::Trans),
                (6, BondStereoRecord::AtropCw),
                (7, BondStereoRecord::AtropCcw),
            ],
            8,
        );
        verify(
            &[
                (0, RingFindRecord::OtherOrUnknown),
                (1, RingFindRecord::Fast),
                (2, RingFindRecord::Sssr),
                (3, RingFindRecord::SymmSssr),
            ],
            4,
        );
        verify(
            &[
                (0, StereoGroupKindRecord::Absolute),
                (1, StereoGroupKindRecord::Or),
                (2, StereoGroupKindRecord::And),
            ],
            3,
        );
        verify(
            &[(2, DimensionRecord::TwoD), (3, DimensionRecord::ThreeD)],
            4,
        );
        verify(
            &[
                (0, BondRoleRecord::Crossing),
                (1, BondRoleRecord::Contained),
            ],
            2,
        );
        verify(
            &[
                (0, PropertyTargetRecord::Atom),
                (1, PropertyTargetRecord::Bond),
            ],
            2,
        );
        verify(
            &[
                (0, SGroupKindRecord::Data),
                (1, SGroupKindRecord::Superatom),
                (2, SGroupKindRecord::MultipleGroup),
                (3, SGroupKindRecord::StructuralRepeatUnit),
                (4, SGroupKindRecord::Monomer),
                (5, SGroupKindRecord::Copolymer),
                (6, SGroupKindRecord::Crosslink),
                (7, SGroupKindRecord::Graft),
                (8, SGroupKindRecord::Modification),
                (9, SGroupKindRecord::Mer),
                (10, SGroupKindRecord::AnyPolymer),
                (11, SGroupKindRecord::MixtureComponent),
                (12, SGroupKindRecord::Mixture),
                (13, SGroupKindRecord::Formulation),
                (14, SGroupKindRecord::Generic("unknown\0text".into())),
            ],
            15,
        );
        verify(
            &[
                (0, ConnectionRecord::HeadToHead),
                (1, ConnectionRecord::HeadToTail),
                (2, ConnectionRecord::Either),
                (3, ConnectionRecord::Unknown("unknown\0text".into())),
            ],
            4,
        );
        verify(
            &[
                (0, BracketRecord::Bracket),
                (1, BracketRecord::Parenthesis),
                (2, BracketRecord::None),
                (3, BracketRecord::Unknown("unknown\0text".into())),
            ],
            4,
        );
        assert!(musli::storage::from_slice::<Value>(&[7, 0]).is_err());
        assert!(musli::wire::from_slice::<Value>(&[0x82, 0xc7, 0x80]).is_err());
    }

    #[test]
    fn archive20_property_variants_float_bits_and_embedded_nul_survive() {
        let mut r = fixture();
        let values = [
            PropertyValue::String("a\0中".into()),
            PropertyValue::StringVector(vec![vec![0xff, 0].into(), "rank".into()]),
            PropertyValue::Int(i32::MIN),
            PropertyValue::UInt(u32::MAX),
            PropertyValue::IntVector(vec![i32::MIN, 0, i32::MAX]),
            PropertyValue::Double(f64::from_bits(0x7ff0_1234_5678_9abc)),
            PropertyValue::Bool(true),
            PropertyValue::Double(-0.0),
        ];
        for (i, v) in values.into_iter().enumerate() {
            r.topology.atoms[0]
                .set_prop(PropertyText::from(vec![0xff, 0, i as u8]), v)
                .unwrap();
        }
        let data = encode_molecule_binary(&input(&r)).unwrap();
        let restored = decode(&data).unwrap();
        assert_eq!(r.topology, restored.topology);
        assert_eq!(data, encode_molecule_binary(&input(&restored)).unwrap());
    }

    #[test]
    fn recovery_chem29_archive20_malformed_wire_members_reject_with_typed_cause_without_schema_change()
     {
        use cosmolkit_model::StereoGroupError;
        let record = fixture();
        let good = encode_molecule_binary(&input(&record)).unwrap();
        for (atoms, bonds, expected) in [
            (vec![0, 0], vec![0, 0], StereoGroupError::DuplicateAtom),
            (vec![0], vec![0, 0], StereoGroupError::DuplicateBond),
        ] {
            let malformed = mutate_molecule(&good, |m| {
                m.stereo_groups.push(StereoRecord {
                    id: Some(7),
                    write_id: 9,
                    kind: StereoGroupKindRecord::Or,
                    atoms,
                    bonds,
                })
            });
            assert_eq!(&malformed[..12], &good[..12]); // Same actual archive 2.0 envelope.
            let error = decode_molecule_binary(&malformed).unwrap_err();
            assert_eq!(error, PickleError::StereoGroup(expected));
            assert_eq!(
                std::error::Error::source(&error)
                    .unwrap()
                    .downcast_ref::<StereoGroupError>(),
                Some(&expected)
            );
        }
    }
}
