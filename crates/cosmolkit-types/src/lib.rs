//! Public, dependency-light chemical vocabulary shared by COSMolKit crates.
//!
//! This crate intentionally has no dependency on a molecule owner, parser,
//! operation runtime, or chemistry algorithm implementation.

use std::fmt;

mod element;
pub mod periodic_table;

pub use element::{ELEMENTS, ELEMENTS_WITH_DUMMY, Element, ElementInfo, ElementParseError};

/// Source-compatible molecular bond order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BondOrder {
    Null,
    Single,
    Double,
    Triple,
    Quadruple,
    Quintuple,
    Hextuple,
    OneAndHalf,
    TwoAndHalf,
    ThreeAndHalf,
    FourAndHalf,
    FiveAndHalf,
    Aromatic,
    Ionic,
    Dative,
    DativeOne,
    DativeLeft,
    DativeRight,
    Hydrogen,
    ThreeCenter,
    Other,
    Zero,
    Unspecified,
}

/// RDKit-compatible atom chirality tag.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ChiralTag {
    Unspecified,
    TetrahedralCw,
    TetrahedralCcw,
    Other,
    Tetrahedral,
    Allene,
    SquarePlanar,
    TrigonalBipyramidal,
    Octahedral,
}

impl ChiralTag {
    #[must_use]
    pub const fn python_code(self) -> i64 {
        match self {
            Self::Unspecified => 0,
            Self::TetrahedralCw => 1,
            Self::TetrahedralCcw => 2,
            Self::Other => 3,
            Self::Tetrahedral => 4,
            Self::Allene => 5,
            Self::SquarePlanar => 6,
            Self::TrigonalBipyramidal => 7,
            Self::Octahedral => 8,
        }
    }

    #[must_use]
    pub const fn rdkit_name(self) -> &'static str {
        match self {
            Self::Unspecified => "CHI_UNSPECIFIED",
            Self::TetrahedralCw => "CHI_TETRAHEDRAL_CW",
            Self::TetrahedralCcw => "CHI_TETRAHEDRAL_CCW",
            Self::Other => "CHI_OTHER",
            Self::Tetrahedral => "CHI_TETRAHEDRAL",
            Self::Allene => "CHI_ALLENE",
            Self::SquarePlanar => "CHI_SQUAREPLANAR",
            Self::TrigonalBipyramidal => "CHI_TRIGONALBIPYRAMIDAL",
            Self::Octahedral => "CHI_OCTAHEDRAL",
        }
    }

    #[must_use]
    pub fn from_rdkit_name(name: &str) -> Option<Self> {
        match name {
            "CHI_UNSPECIFIED" => Some(Self::Unspecified),
            "CHI_TETRAHEDRAL_CW" => Some(Self::TetrahedralCw),
            "CHI_TETRAHEDRAL_CCW" => Some(Self::TetrahedralCcw),
            "CHI_OTHER" => Some(Self::Other),
            "CHI_TETRAHEDRAL" => Some(Self::Tetrahedral),
            "CHI_ALLENE" => Some(Self::Allene),
            "CHI_SQUAREPLANAR" => Some(Self::SquarePlanar),
            "CHI_TRIGONALBIPYRAMIDAL" => Some(Self::TrigonalBipyramidal),
            "CHI_OCTAHEDRAL" => Some(Self::Octahedral),
            _ => None,
        }
    }
}

/// RDKit-compatible atom hybridization classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Hybridization {
    Unspecified,
    S,
    Sp,
    Sp2,
    Sp3,
    Sp2d,
    Sp3d,
    Sp3d2,
    Other,
}

/// Direction annotation used by 2D stereochemical bonds.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BondDirection {
    None,
    BeginWedge,
    BeginDash,
    EndUpRight,
    EndDownRight,
    EitherDouble,
    Unknown,
}

impl BondDirection {
    #[must_use]
    pub const fn python_code(self) -> i64 {
        match self {
            Self::None => 0,
            Self::BeginWedge => 1,
            Self::BeginDash => 2,
            Self::EndUpRight => 3,
            Self::EndDownRight => 4,
            Self::EitherDouble => 5,
            Self::Unknown => 6,
        }
    }

    #[must_use]
    pub const fn rdkit_name(self) -> &'static str {
        match self {
            Self::None => "NONE",
            Self::BeginWedge => "BEGINWEDGE",
            Self::BeginDash => "BEGINDASH",
            Self::EndUpRight => "ENDUPRIGHT",
            Self::EndDownRight => "ENDDOWNRIGHT",
            Self::EitherDouble => "EITHERDOUBLE",
            Self::Unknown => "UNKNOWN",
        }
    }
}

/// Double-bond and axial stereochemistry annotation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BondStereo {
    None,
    Any,
    Z,
    E,
    Cis,
    Trans,
    AtropCw,
    AtropCcw,
}

impl BondStereo {
    #[must_use]
    pub const fn python_code(self) -> i64 {
        match self {
            Self::None => 0,
            Self::Any => 1,
            Self::Z => 2,
            Self::E => 3,
            Self::Cis => 4,
            Self::Trans => 5,
            Self::AtropCw => 6,
            Self::AtropCcw => 7,
        }
    }

    #[must_use]
    pub const fn rdkit_name(self) -> &'static str {
        match self {
            Self::None => "NONE",
            Self::Any => "ANY",
            Self::Z => "Z",
            Self::E => "E",
            Self::Cis => "CIS",
            Self::Trans => "TRANS",
            Self::AtropCw => "ATROP_CW",
            Self::AtropCcw => "ATROP_CCW",
        }
    }
}

impl BondOrder {
    #[must_use]
    pub const fn python_code(self) -> i64 {
        match self {
            Self::Null | Self::Unspecified => 0,
            Self::Single => 1,
            Self::Double => 2,
            Self::Triple => 3,
            Self::Quadruple => 4,
            Self::Quintuple => 5,
            Self::Hextuple => 6,
            Self::OneAndHalf => 7,
            Self::TwoAndHalf => 8,
            Self::ThreeAndHalf => 9,
            Self::FourAndHalf => 10,
            Self::FiveAndHalf => 11,
            Self::Aromatic => 12,
            Self::Ionic => 13,
            Self::Dative => 14,
            Self::DativeOne => 15,
            Self::DativeLeft => 16,
            Self::DativeRight => 17,
            Self::Hydrogen => 18,
            Self::ThreeCenter => 19,
            Self::Other => 20,
            Self::Zero => 21,
        }
    }

    #[must_use]
    pub const fn rdkit_name(self) -> &'static str {
        match self {
            Self::Null | Self::Unspecified => "UNSPECIFIED",
            Self::Single => "SINGLE",
            Self::Double => "DOUBLE",
            Self::Triple => "TRIPLE",
            Self::Quadruple => "QUADRUPLE",
            Self::Quintuple => "QUINTUPLE",
            Self::Hextuple => "HEXTUPLE",
            Self::OneAndHalf => "ONEANDAHALF",
            Self::TwoAndHalf => "TWOANDAHALF",
            Self::ThreeAndHalf => "THREEANDAHALF",
            Self::FourAndHalf => "FOURANDAHALF",
            Self::FiveAndHalf => "FIVEANDAHALF",
            Self::Aromatic => "AROMATIC",
            Self::Ionic => "IONIC",
            Self::Dative => "DATIVE",
            Self::DativeOne => "DATIVEONE",
            Self::DativeLeft => "DATIVEL",
            Self::DativeRight => "DATIVER",
            Self::Hydrogen => "HYDROGEN",
            Self::ThreeCenter => "THREECENTER",
            Self::Other => "OTHER",
            Self::Zero => "ZERO",
        }
    }

    #[must_use]
    pub fn from_rdkit_name(name: &str) -> Option<Self> {
        match name {
            "UNSPECIFIED" => Some(Self::Unspecified),
            "ZERO" => Some(Self::Zero),
            "SINGLE" => Some(Self::Single),
            "DOUBLE" => Some(Self::Double),
            "TRIPLE" => Some(Self::Triple),
            "QUADRUPLE" => Some(Self::Quadruple),
            "QUINTUPLE" => Some(Self::Quintuple),
            "HEXTUPLE" => Some(Self::Hextuple),
            "ONEANDAHALF" => Some(Self::OneAndHalf),
            "TWOANDAHALF" => Some(Self::TwoAndHalf),
            "THREEANDAHALF" => Some(Self::ThreeAndHalf),
            "FOURANDAHALF" => Some(Self::FourAndHalf),
            "FIVEANDAHALF" => Some(Self::FiveAndHalf),
            "AROMATIC" => Some(Self::Aromatic),
            "IONIC" => Some(Self::Ionic),
            "DATIVE" => Some(Self::Dative),
            "DATIVEONE" => Some(Self::DativeOne),
            "DATIVEL" => Some(Self::DativeLeft),
            "DATIVER" => Some(Self::DativeRight),
            "HYDROGEN" => Some(Self::Hydrogen),
            "THREECENTER" => Some(Self::ThreeCenter),
            "OTHER" => Some(Self::Other),
            _ => None,
        }
    }
}

impl fmt::Display for BondOrder {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.rdkit_name())
    }
}
