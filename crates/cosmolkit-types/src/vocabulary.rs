use std::{fmt, str::FromStr};

use serde::{Deserialize, Deserializer, Serialize, Serializer};
use thiserror::Error;

/// Stable chemical-element identity.
///
/// The private atomic number is always in the inclusive range `0..=118`.
/// Zero is the source-compatible dummy atom (`*`); `1..=118` are H through
/// Og.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Element {
    atomic_number: u8,
}

/// Error returned when a string is not a source-recognized element symbol.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
#[error("unknown element symbol '{input}'")]
pub struct ElementParseError {
    input: String,
}

impl ElementParseError {
    #[must_use]
    pub fn input(&self) -> &str {
        &self.input
    }
}

/// Source-aligned periodic-table metadata returned by the algorithm layer.
///
/// This crate owns only the dependency-light result schema. The numerical
/// table and the function that populates this record belong to `CORE-tables`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct ElementInfo {
    pub element: Element,
    pub symbol: &'static str,
    pub atomic_number: u8,
    pub period: u8,
    pub outer_electrons: i32,
    pub valences: &'static [i32],
    /// RDKit's source-defined `Rb0` bond radius in angstroms.
    pub rb0: f64,
    pub atomic_weight: f64,
}

macro_rules! define_elements {
    ($( $constant:ident = $number:literal => $symbol:literal ),+ $(,)?) => {
        impl Element {
            $(pub const $constant: Self = Self { atomic_number: $number };)+

            /// Construct an element from a source-recognized symbol.
            #[must_use]
            pub fn from_symbol(symbol: &str) -> Option<Self> {
                // BEGIN RDKIT CPP FUNCTION PeriodicTable::getAtomicNumber /
                // PeriodicTable::PeriodicTable
                // RDKit✔️✔️: int getAtomicNumber(const std::string &elementSymbol) const {
                // RDKit✔️✔️:   int anum = -1;
                // RDKit✔️✔️:   if (elementSymbol == "C") {
                // RDKit✔️✔️:     anum = 6;
                // RDKit✔️✔️:   } else if (elementSymbol == "N") {
                // RDKit✔️✔️:     anum = 7;
                // RDKit✔️✔️:   } else if (elementSymbol == "O") {
                // RDKit✔️✔️:     anum = 8;
                // RDKit✔️✔️:   } else {
                // RDKit✔️✔️:     STR_UINT_MAP::const_iterator iter = byname.find(elementSymbol);
                // RDKit✔️✔️:     if (iter != byname.end()) {
                // RDKit✔️✔️:       anum = iter->second;
                // RDKit✔️✔️:     }
                // RDKit✔️✔️:   }
                // RDKit✔️✔️:   POSTCONDITION(anum > -1, "Element '" + elementSymbol + "' not found");
                // RDKit✔️✔️:   return anum;
                // RDKit✔️✔️: }
                // RDKit✔️✔️: std::string enam = adata.Symbol();
                // RDKit✔️✔️: byname[enam] = adata.AtomicNum();
                // END RDKIT CPP FUNCTION PeriodicTable::getAtomicNumber /
                // PeriodicTable::PeriodicTable
                match symbol {
                    $($symbol => Some(Self::$constant),)+
                    "Uut" => Some(Self::NH),
                    "Uup" => Some(Self::MC),
                    _ => None,
                }
            }
        }

        const SYMBOLS: [&str; 119] = [$($symbol),+];
    };
}

define_elements! {
    DUMMY = 0 => "*", H = 1 => "H", HE = 2 => "He", LI = 3 => "Li",
    BE = 4 => "Be", B = 5 => "B", C = 6 => "C", N = 7 => "N",
    O = 8 => "O", F = 9 => "F", NE = 10 => "Ne", NA = 11 => "Na",
    MG = 12 => "Mg", AL = 13 => "Al", SI = 14 => "Si", P = 15 => "P",
    S = 16 => "S", CL = 17 => "Cl", AR = 18 => "Ar", K = 19 => "K",
    CA = 20 => "Ca", SC = 21 => "Sc", TI = 22 => "Ti", V = 23 => "V",
    CR = 24 => "Cr", MN = 25 => "Mn", FE = 26 => "Fe", CO = 27 => "Co",
    NI = 28 => "Ni", CU = 29 => "Cu", ZN = 30 => "Zn", GA = 31 => "Ga",
    GE = 32 => "Ge", AS = 33 => "As", SE = 34 => "Se", BR = 35 => "Br",
    KR = 36 => "Kr", RB = 37 => "Rb", SR = 38 => "Sr", Y = 39 => "Y",
    ZR = 40 => "Zr", NB = 41 => "Nb", MO = 42 => "Mo", TC = 43 => "Tc",
    RU = 44 => "Ru", RH = 45 => "Rh", PD = 46 => "Pd", AG = 47 => "Ag",
    CD = 48 => "Cd", IN = 49 => "In", SN = 50 => "Sn", SB = 51 => "Sb",
    TE = 52 => "Te", I = 53 => "I", XE = 54 => "Xe", CS = 55 => "Cs",
    BA = 56 => "Ba", LA = 57 => "La", CE = 58 => "Ce", PR = 59 => "Pr",
    ND = 60 => "Nd", PM = 61 => "Pm", SM = 62 => "Sm", EU = 63 => "Eu",
    GD = 64 => "Gd", TB = 65 => "Tb", DY = 66 => "Dy", HO = 67 => "Ho",
    ER = 68 => "Er", TM = 69 => "Tm", YB = 70 => "Yb", LU = 71 => "Lu",
    HF = 72 => "Hf", TA = 73 => "Ta", W = 74 => "W", RE = 75 => "Re",
    OS = 76 => "Os", IR = 77 => "Ir", PT = 78 => "Pt", AU = 79 => "Au",
    HG = 80 => "Hg", TL = 81 => "Tl", PB = 82 => "Pb", BI = 83 => "Bi",
    PO = 84 => "Po", AT = 85 => "At", RN = 86 => "Rn", FR = 87 => "Fr",
    RA = 88 => "Ra", AC = 89 => "Ac", TH = 90 => "Th", PA = 91 => "Pa",
    U = 92 => "U", NP = 93 => "Np", PU = 94 => "Pu", AM = 95 => "Am",
    CM = 96 => "Cm", BK = 97 => "Bk", CF = 98 => "Cf", ES = 99 => "Es",
    FM = 100 => "Fm", MD = 101 => "Md", NO = 102 => "No", LR = 103 => "Lr",
    RF = 104 => "Rf", DB = 105 => "Db", SG = 106 => "Sg", BH = 107 => "Bh",
    HS = 108 => "Hs", MT = 109 => "Mt", DS = 110 => "Ds", RG = 111 => "Rg",
    CN = 112 => "Cn", NH = 113 => "Nh", FL = 114 => "Fl", MC = 115 => "Mc",
    LV = 116 => "Lv", TS = 117 => "Ts", OG = 118 => "Og",
}

impl Element {
    /// Construct an element from a checked atomic number.
    #[must_use]
    pub const fn from_atomic_number(atomic_number: u8) -> Option<Self> {
        // BEGIN RDKIT CPP FUNCTION PeriodicTable::getMaxAtomicNumber /
        // PeriodicTable::getElementSymbol
        // RDKit✔️✔️: UINT getMaxAtomicNumber() const { return byanum.size() - 1; }
        // RDKit✔️✔️: std::string getElementSymbol(UINT atomicNumber) const {
        // RDKit✔️✔️:   PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
        // RDKit✔️✔️:   return byanum[atomicNumber].Symbol();
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION PeriodicTable::getMaxAtomicNumber /
        // PeriodicTable::getElementSymbol
        if atomic_number <= 118 {
            Some(Self { atomic_number })
        } else {
            None
        }
    }

    #[must_use]
    pub const fn atomic_number(self) -> u8 {
        self.atomic_number
    }

    /// Return the canonical element symbol (`*`, `H` through `Og`).
    #[must_use]
    pub fn symbol(self) -> &'static str {
        // BEGIN RDKIT CPP FUNCTION PeriodicTable::getElementSymbol
        // RDKit✔️✔️: std::string getElementSymbol(UINT atomicNumber) const {
        // RDKit✔️✔️:   PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
        // RDKit✔️✔️:   return byanum[atomicNumber].Symbol();
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION PeriodicTable::getElementSymbol
        SYMBOLS[usize::from(self.atomic_number)]
    }

    pub fn iter() -> impl ExactSizeIterator<Item = Self> + DoubleEndedIterator {
        ELEMENTS.iter().copied()
    }

    pub fn iter_with_dummy() -> impl ExactSizeIterator<Item = Self> + DoubleEndedIterator {
        ELEMENTS_WITH_DUMMY.iter().copied()
    }
}

impl fmt::Display for Element {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.symbol())
    }
}

impl FromStr for Element {
    type Err = ElementParseError;

    fn from_str(input: &str) -> Result<Self, Self::Err> {
        Self::from_symbol(input).ok_or_else(|| ElementParseError {
            input: input.to_string(),
        })
    }
}

impl Serialize for Element {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_str(self.symbol())
    }
}

impl<'de> Deserialize<'de> for Element {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let symbol = String::deserialize(deserializer)?;
        symbol.parse().map_err(serde::de::Error::custom)
    }
}

const fn element_array<const N: usize>(first_atomic_number: u8) -> [Element; N] {
    let mut elements = [Element::DUMMY; N];
    let mut index = 0;
    while index < N {
        elements[index] = Element {
            atomic_number: first_atomic_number + index as u8,
        };
        index += 1;
    }
    elements
}

/// All 118 real elements in ascending atomic-number order (H through Og).
pub static ELEMENTS: [Element; 118] = element_array::<118>(1);

/// The dummy atom followed by all real elements in atomic-number order.
pub static ELEMENTS_WITH_DUMMY: [Element; 119] = element_array::<119>(0);

/// Source-compatible molecular bond order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i64)]
pub enum BondOrder {
    Unspecified = 0,
    Single = 1,
    Double = 2,
    Triple = 3,
    Quadruple = 4,
    Quintuple = 5,
    Hextuple = 6,
    OneAndHalf = 7,
    TwoAndHalf = 8,
    ThreeAndHalf = 9,
    FourAndHalf = 10,
    FiveAndHalf = 11,
    Aromatic = 12,
    Ionic = 13,
    Hydrogen = 14,
    ThreeCenter = 15,
    DativeOne = 16,
    Dative = 17,
    DativeLeft = 18,
    DativeRight = 19,
    Other = 20,
    Zero = 21,
}

impl BondOrder {
    #[must_use]
    pub const fn rdkit_code(self) -> i64 {
        // RDKit✔️✔️: typedef enum {
        // RDKit✔️✔️:   UNSPECIFIED = 0,
        // RDKit✔️✔️:   SINGLE,
        // RDKit✔️✔️:   DOUBLE,
        // RDKit✔️✔️:   TRIPLE,
        // RDKit✔️✔️:   QUADRUPLE,
        // RDKit✔️✔️:   QUINTUPLE,
        // RDKit✔️✔️:   HEXTUPLE,
        // RDKit✔️✔️:   ONEANDAHALF,
        // RDKit✔️✔️:   TWOANDAHALF,
        // RDKit✔️✔️:   THREEANDAHALF,
        // RDKit✔️✔️:   FOURANDAHALF,
        // RDKit✔️✔️:   FIVEANDAHALF,
        // RDKit✔️✔️:   AROMATIC,
        // RDKit✔️✔️:   IONIC,
        // RDKit✔️✔️:   HYDROGEN,
        // RDKit✔️✔️:   THREECENTER,
        // RDKit✔️✔️:   DATIVEONE,  //!< one-electron dative (e.g. from a C in a Cp ring to a metal)
        // RDKit✔️✔️:   DATIVE,     //!< standard two-electron dative
        // RDKit✔️✔️:   DATIVEL,    //!< standard two-electron dative
        // RDKit✔️✔️:   DATIVER,    //!< standard two-electron dative
        // RDKit✔️✔️:   OTHER,
        // RDKit✔️✔️:   ZERO  //!< Zero-order bond (from
        // RDKit✔️✔️:   // http://pubs.acs.org/doi/abs/10.1021/ci200488k)
        // RDKit✔️✔️: } BondType;
        self as i64
    }

    #[must_use]
    pub const fn rdkit_name(self) -> &'static str {
        // RDKit✔️✔️: python::enum_<Bond::BondType>("BondType")
        // RDKit✔️✔️:     .value("UNSPECIFIED", Bond::UNSPECIFIED)
        // RDKit✔️✔️:     .value("SINGLE", Bond::SINGLE)
        // RDKit✔️✔️:     .value("DOUBLE", Bond::DOUBLE)
        // RDKit✔️✔️:     .value("TRIPLE", Bond::TRIPLE)
        // RDKit✔️✔️:     .value("QUADRUPLE", Bond::QUADRUPLE)
        // RDKit✔️✔️:     .value("QUINTUPLE", Bond::QUINTUPLE)
        // RDKit✔️✔️:     .value("HEXTUPLE", Bond::HEXTUPLE)
        // RDKit✔️✔️:     .value("ONEANDAHALF", Bond::ONEANDAHALF)
        // RDKit✔️✔️:     .value("TWOANDAHALF", Bond::TWOANDAHALF)
        // RDKit✔️✔️:     .value("THREEANDAHALF", Bond::THREEANDAHALF)
        // RDKit✔️✔️:     .value("FOURANDAHALF", Bond::FOURANDAHALF)
        // RDKit✔️✔️:     .value("FIVEANDAHALF", Bond::FIVEANDAHALF)
        // RDKit✔️✔️:     .value("AROMATIC", Bond::AROMATIC)
        // RDKit✔️✔️:     .value("IONIC", Bond::IONIC)
        // RDKit✔️✔️:     .value("HYDROGEN", Bond::HYDROGEN)
        // RDKit✔️✔️:     .value("THREECENTER", Bond::THREECENTER)
        // RDKit✔️✔️:     .value("DATIVEONE", Bond::DATIVEONE)
        // RDKit✔️✔️:     .value("DATIVE", Bond::DATIVE)
        // RDKit✔️✔️:     .value("DATIVEL", Bond::DATIVEL)
        // RDKit✔️✔️:     .value("DATIVER", Bond::DATIVER)
        // RDKit✔️✔️:     .value("OTHER", Bond::OTHER)
        // RDKit✔️✔️:     .value("ZERO", Bond::ZERO);
        match self {
            Self::Unspecified => "UNSPECIFIED",
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
            Self::Hydrogen => "HYDROGEN",
            Self::ThreeCenter => "THREECENTER",
            Self::DativeOne => "DATIVEONE",
            Self::Dative => "DATIVE",
            Self::DativeLeft => "DATIVEL",
            Self::DativeRight => "DATIVER",
            Self::Other => "OTHER",
            Self::Zero => "ZERO",
        }
    }

    #[must_use]
    pub const fn from_rdkit_code(code: i64) -> Option<Self> {
        // RDKit✔️✔️: typedef enum {
        // RDKit✔️✔️:   UNSPECIFIED = 0,
        // RDKit✔️✔️:   SINGLE,
        // RDKit✔️✔️:   DOUBLE,
        // RDKit✔️✔️:   TRIPLE,
        // RDKit✔️✔️:   QUADRUPLE,
        // RDKit✔️✔️:   QUINTUPLE,
        // RDKit✔️✔️:   HEXTUPLE,
        // RDKit✔️✔️:   ONEANDAHALF,
        // RDKit✔️✔️:   TWOANDAHALF,
        // RDKit✔️✔️:   THREEANDAHALF,
        // RDKit✔️✔️:   FOURANDAHALF,
        // RDKit✔️✔️:   FIVEANDAHALF,
        // RDKit✔️✔️:   AROMATIC,
        // RDKit✔️✔️:   IONIC,
        // RDKit✔️✔️:   HYDROGEN,
        // RDKit✔️✔️:   THREECENTER,
        // RDKit✔️✔️:   DATIVEONE,  //!< one-electron dative (e.g. from a C in a Cp ring to a metal)
        // RDKit✔️✔️:   DATIVE,     //!< standard two-electron dative
        // RDKit✔️✔️:   DATIVEL,    //!< standard two-electron dative
        // RDKit✔️✔️:   DATIVER,    //!< standard two-electron dative
        // RDKit✔️✔️:   OTHER,
        // RDKit✔️✔️:   ZERO  //!< Zero-order bond (from
        // RDKit✔️✔️:   // http://pubs.acs.org/doi/abs/10.1021/ci200488k)
        // RDKit✔️✔️: } BondType;
        match code {
            0 => Some(Self::Unspecified),
            1 => Some(Self::Single),
            2 => Some(Self::Double),
            3 => Some(Self::Triple),
            4 => Some(Self::Quadruple),
            5 => Some(Self::Quintuple),
            6 => Some(Self::Hextuple),
            7 => Some(Self::OneAndHalf),
            8 => Some(Self::TwoAndHalf),
            9 => Some(Self::ThreeAndHalf),
            10 => Some(Self::FourAndHalf),
            11 => Some(Self::FiveAndHalf),
            12 => Some(Self::Aromatic),
            13 => Some(Self::Ionic),
            14 => Some(Self::Hydrogen),
            15 => Some(Self::ThreeCenter),
            16 => Some(Self::DativeOne),
            17 => Some(Self::Dative),
            18 => Some(Self::DativeLeft),
            19 => Some(Self::DativeRight),
            20 => Some(Self::Other),
            21 => Some(Self::Zero),
            _ => None,
        }
    }

    #[must_use]
    pub fn from_rdkit_name(name: &str) -> Option<Self> {
        // This is the checked inverse of the exact Wrap/Bond.cpp table copied
        // into `rdkit_name` above.
        match name {
            "UNSPECIFIED" => Some(Self::Unspecified),
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
            "HYDROGEN" => Some(Self::Hydrogen),
            "THREECENTER" => Some(Self::ThreeCenter),
            "DATIVEONE" => Some(Self::DativeOne),
            "DATIVE" => Some(Self::Dative),
            "DATIVEL" => Some(Self::DativeLeft),
            "DATIVER" => Some(Self::DativeRight),
            "OTHER" => Some(Self::Other),
            "ZERO" => Some(Self::Zero),
            _ => None,
        }
    }
}

/// RDKit-compatible atom chirality tag.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i64)]
pub enum ChiralTag {
    Unspecified = 0,
    TetrahedralCw = 1,
    TetrahedralCcw = 2,
    Other = 3,
    Tetrahedral = 4,
    Allene = 5,
    SquarePlanar = 6,
    TrigonalBipyramidal = 7,
    Octahedral = 8,
}

impl ChiralTag {
    #[must_use]
    pub const fn rdkit_code(self) -> i64 {
        // RDKit✔️✔️: typedef enum {
        // RDKit✔️✔️:   CHI_UNSPECIFIED = 0,  //!< chirality that hasn't been specified
        // RDKit✔️✔️:   CHI_TETRAHEDRAL_CW,   //!< tetrahedral: clockwise rotation (SMILES \@\@)
        // RDKit✔️✔️:   CHI_TETRAHEDRAL_CCW,  //!< tetrahedral: counter-clockwise rotation (SMILES
        // RDKit✔️✔️:                           //\@)
        // RDKit✔️✔️:   CHI_OTHER,            //!< some unrecognized type of chirality
        // RDKit✔️✔️:   CHI_TETRAHEDRAL,      //!< tetrahedral, use permutation flag
        // RDKit✔️✔️:   CHI_ALLENE,           //!< allene, use permutation flag
        // RDKit✔️✔️:   CHI_SQUAREPLANAR,     //!< square planar, use permutation flag
        // RDKit✔️✔️:   CHI_TRIGONALBIPYRAMIDAL,  //!< trigonal bipyramidal, use permutation flag
        // RDKit✔️✔️:   CHI_OCTAHEDRAL            //!< octahedral, use permutation flag
        // RDKit✔️✔️: } ChiralType;
        self as i64
    }

    #[must_use]
    pub const fn rdkit_name(self) -> &'static str {
        // RDKit✔️✔️: python::enum_<Atom::ChiralType>("ChiralType")
        // RDKit✔️✔️:     .value("CHI_UNSPECIFIED", Atom::CHI_UNSPECIFIED)
        // RDKit✔️✔️:     .value("CHI_TETRAHEDRAL_CW", Atom::CHI_TETRAHEDRAL_CW)
        // RDKit✔️✔️:     .value("CHI_TETRAHEDRAL_CCW", Atom::CHI_TETRAHEDRAL_CCW)
        // RDKit✔️✔️:     .value("CHI_OTHER", Atom::CHI_OTHER)
        // RDKit✔️✔️:     .value("CHI_TETRAHEDRAL", Atom::CHI_TETRAHEDRAL)
        // RDKit✔️✔️:     .value("CHI_ALLENE", Atom::CHI_ALLENE)
        // RDKit✔️✔️:     .value("CHI_SQUAREPLANAR", Atom::CHI_SQUAREPLANAR)
        // RDKit✔️✔️:     .value("CHI_TRIGONALBIPYRAMIDAL", Atom::CHI_TRIGONALBIPYRAMIDAL)
        // RDKit✔️✔️:     .value("CHI_OCTAHEDRAL", Atom::CHI_OCTAHEDRAL)
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
    pub const fn from_rdkit_code(code: i64) -> Option<Self> {
        // RDKit✔️✔️: typedef enum {
        // RDKit✔️✔️:   CHI_UNSPECIFIED = 0,  //!< chirality that hasn't been specified
        // RDKit✔️✔️:   CHI_TETRAHEDRAL_CW,   //!< tetrahedral: clockwise rotation (SMILES \@\@)
        // RDKit✔️✔️:   CHI_TETRAHEDRAL_CCW,  //!< tetrahedral: counter-clockwise rotation (SMILES
        // RDKit✔️✔️:                           //\@)
        // RDKit✔️✔️:   CHI_OTHER,            //!< some unrecognized type of chirality
        // RDKit✔️✔️:   CHI_TETRAHEDRAL,      //!< tetrahedral, use permutation flag
        // RDKit✔️✔️:   CHI_ALLENE,           //!< allene, use permutation flag
        // RDKit✔️✔️:   CHI_SQUAREPLANAR,     //!< square planar, use permutation flag
        // RDKit✔️✔️:   CHI_TRIGONALBIPYRAMIDAL,  //!< trigonal bipyramidal, use permutation flag
        // RDKit✔️✔️:   CHI_OCTAHEDRAL            //!< octahedral, use permutation flag
        // RDKit✔️✔️: } ChiralType;
        match code {
            0 => Some(Self::Unspecified),
            1 => Some(Self::TetrahedralCw),
            2 => Some(Self::TetrahedralCcw),
            3 => Some(Self::Other),
            4 => Some(Self::Tetrahedral),
            5 => Some(Self::Allene),
            6 => Some(Self::SquarePlanar),
            7 => Some(Self::TrigonalBipyramidal),
            8 => Some(Self::Octahedral),
            _ => None,
        }
    }

    #[must_use]
    pub fn from_rdkit_name(name: &str) -> Option<Self> {
        // This is the checked inverse of the exact Wrap/Atom.cpp table copied
        // into `rdkit_name` above.
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
#[repr(i64)]
pub enum Hybridization {
    Unspecified = 0,
    S = 1,
    Sp = 2,
    Sp2 = 3,
    Sp3 = 4,
    Sp2d = 5,
    Sp3d = 6,
    Sp3d2 = 7,
    Other = 8,
}

impl Hybridization {
    #[must_use]
    pub const fn rdkit_code(self) -> i64 {
        // RDKit✔️✔️: typedef enum {
        // RDKit✔️✔️:   UNSPECIFIED = 0,  //!< hybridization that hasn't been specified
        // RDKit✔️✔️:   S,
        // RDKit✔️✔️:   SP,
        // RDKit✔️✔️:   SP2,
        // RDKit✔️✔️:   SP3,
        // RDKit✔️✔️:   SP2D,
        // RDKit✔️✔️:   SP3D,
        // RDKit✔️✔️:   SP3D2,
        // RDKit✔️✔️:   OTHER  //!< unrecognized hybridization
        // RDKit✔️✔️: } HybridizationType;
        self as i64
    }

    #[must_use]
    pub const fn rdkit_name(self) -> &'static str {
        // RDKit✔️✔️: python::enum_<Atom::HybridizationType>("HybridizationType")
        // RDKit✔️✔️:     .value("UNSPECIFIED", Atom::UNSPECIFIED)
        // RDKit✔️✔️:     .value("S", Atom::S)
        // RDKit✔️✔️:     .value("SP", Atom::SP)
        // RDKit✔️✔️:     .value("SP2", Atom::SP2)
        // RDKit✔️✔️:     .value("SP3", Atom::SP3)
        // RDKit✔️✔️:     .value("SP2D", Atom::SP2D)
        // RDKit✔️✔️:     .value("SP3D", Atom::SP3D)
        // RDKit✔️✔️:     .value("SP3D2", Atom::SP3D2)
        // RDKit✔️✔️:     .value("OTHER", Atom::OTHER);
        match self {
            Self::Unspecified => "UNSPECIFIED",
            Self::S => "S",
            Self::Sp => "SP",
            Self::Sp2 => "SP2",
            Self::Sp3 => "SP3",
            Self::Sp2d => "SP2D",
            Self::Sp3d => "SP3D",
            Self::Sp3d2 => "SP3D2",
            Self::Other => "OTHER",
        }
    }

    #[must_use]
    pub const fn from_rdkit_code(code: i64) -> Option<Self> {
        // RDKit✔️✔️: typedef enum {
        // RDKit✔️✔️:   UNSPECIFIED = 0,  //!< hybridization that hasn't been specified
        // RDKit✔️✔️:   S,
        // RDKit✔️✔️:   SP,
        // RDKit✔️✔️:   SP2,
        // RDKit✔️✔️:   SP3,
        // RDKit✔️✔️:   SP2D,
        // RDKit✔️✔️:   SP3D,
        // RDKit✔️✔️:   SP3D2,
        // RDKit✔️✔️:   OTHER  //!< unrecognized hybridization
        // RDKit✔️✔️: } HybridizationType;
        match code {
            0 => Some(Self::Unspecified),
            1 => Some(Self::S),
            2 => Some(Self::Sp),
            3 => Some(Self::Sp2),
            4 => Some(Self::Sp3),
            5 => Some(Self::Sp2d),
            6 => Some(Self::Sp3d),
            7 => Some(Self::Sp3d2),
            8 => Some(Self::Other),
            _ => None,
        }
    }

    #[must_use]
    pub fn from_rdkit_name(name: &str) -> Option<Self> {
        // This is the checked inverse of the exact Wrap/Atom.cpp table copied
        // into `rdkit_name` above.
        match name {
            "UNSPECIFIED" => Some(Self::Unspecified),
            "S" => Some(Self::S),
            "SP" => Some(Self::Sp),
            "SP2" => Some(Self::Sp2),
            "SP3" => Some(Self::Sp3),
            "SP2D" => Some(Self::Sp2d),
            "SP3D" => Some(Self::Sp3d),
            "SP3D2" => Some(Self::Sp3d2),
            "OTHER" => Some(Self::Other),
            _ => None,
        }
    }
}

/// Direction annotation used by 2D stereochemical bonds.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i64)]
pub enum BondDirection {
    None = 0,
    BeginWedge = 1,
    BeginDash = 2,
    EndDownRight = 3,
    EndUpRight = 4,
    EitherDouble = 5,
    Unknown = 6,
}

impl BondDirection {
    #[must_use]
    pub const fn rdkit_code(self) -> i64 {
        // RDKit✔️✔️: typedef enum {
        // RDKit✔️✔️:   NONE = 0,    //!< no special style
        // RDKit✔️✔️:   BEGINWEDGE,  //!< wedged: narrow at begin
        // RDKit✔️✔️:   BEGINDASH,   //!< dashed: narrow at begin
        // RDKit✔️✔️:   ENDDOWNRIGHT,  //!< for cis/trans
        // RDKit✔️✔️:   ENDUPRIGHT,    //!<  ditto
        // RDKit✔️✔️:   EITHERDOUBLE,  //!< a "crossed" double bond
        // RDKit✔️✔️:   UNKNOWN,       //!< intentionally unspecified stereochemistry
        // RDKit✔️✔️: } BondDir;
        self as i64
    }

    #[must_use]
    pub const fn rdkit_name(self) -> &'static str {
        // RDKit✔️✔️: python::enum_<Bond::BondDir>("BondDir")
        // RDKit✔️✔️:     .value("NONE", Bond::NONE)
        // RDKit✔️✔️:     .value("BEGINWEDGE", Bond::BEGINWEDGE)
        // RDKit✔️✔️:     .value("BEGINDASH", Bond::BEGINDASH)
        // RDKit✔️✔️:     .value("ENDDOWNRIGHT", Bond::ENDDOWNRIGHT)
        // RDKit✔️✔️:     .value("ENDUPRIGHT", Bond::ENDUPRIGHT)
        // RDKit✔️✔️:     .value("EITHERDOUBLE", Bond::EITHERDOUBLE)
        // RDKit✔️✔️:     .value("UNKNOWN", Bond::UNKNOWN);
        match self {
            Self::None => "NONE",
            Self::BeginWedge => "BEGINWEDGE",
            Self::BeginDash => "BEGINDASH",
            Self::EndDownRight => "ENDDOWNRIGHT",
            Self::EndUpRight => "ENDUPRIGHT",
            Self::EitherDouble => "EITHERDOUBLE",
            Self::Unknown => "UNKNOWN",
        }
    }

    #[must_use]
    pub const fn from_rdkit_code(code: i64) -> Option<Self> {
        // RDKit✔️✔️: typedef enum {
        // RDKit✔️✔️:   NONE = 0,    //!< no special style
        // RDKit✔️✔️:   BEGINWEDGE,  //!< wedged: narrow at begin
        // RDKit✔️✔️:   BEGINDASH,   //!< dashed: narrow at begin
        // RDKit✔️✔️:   ENDDOWNRIGHT,  //!< for cis/trans
        // RDKit✔️✔️:   ENDUPRIGHT,    //!<  ditto
        // RDKit✔️✔️:   EITHERDOUBLE,  //!< a "crossed" double bond
        // RDKit✔️✔️:   UNKNOWN,       //!< intentionally unspecified stereochemistry
        // RDKit✔️✔️: } BondDir;
        match code {
            0 => Some(Self::None),
            1 => Some(Self::BeginWedge),
            2 => Some(Self::BeginDash),
            3 => Some(Self::EndDownRight),
            4 => Some(Self::EndUpRight),
            5 => Some(Self::EitherDouble),
            6 => Some(Self::Unknown),
            _ => None,
        }
    }

    #[must_use]
    pub fn from_rdkit_name(name: &str) -> Option<Self> {
        // This is the checked inverse of the exact Wrap/Bond.cpp table copied
        // into `rdkit_name` above.
        match name {
            "NONE" => Some(Self::None),
            "BEGINWEDGE" => Some(Self::BeginWedge),
            "BEGINDASH" => Some(Self::BeginDash),
            "ENDDOWNRIGHT" => Some(Self::EndDownRight),
            "ENDUPRIGHT" => Some(Self::EndUpRight),
            "EITHERDOUBLE" => Some(Self::EitherDouble),
            "UNKNOWN" => Some(Self::Unknown),
            _ => None,
        }
    }
}

/// Double-bond and axial stereochemistry annotation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i64)]
pub enum BondStereo {
    None = 0,
    Any = 1,
    Z = 2,
    E = 3,
    Cis = 4,
    Trans = 5,
    AtropCw = 6,
    AtropCcw = 7,
}

impl BondStereo {
    #[must_use]
    pub const fn rdkit_code(self) -> i64 {
        // RDKit✔️✔️: typedef enum {
        // RDKit✔️✔️:   STEREONONE = 0,  // no special style
        // RDKit✔️✔️:   STEREOANY,       // intentionally unspecified
        // RDKit✔️✔️:   STEREOZ,         // Z double bond
        // RDKit✔️✔️:   STEREOE,         // E double bond
        // RDKit✔️✔️:   STEREOCIS,       // cis double bond
        // RDKit✔️✔️:   STEREOTRANS,     // trans double bond
        // RDKit✔️✔️:   STEREOATROPCW,   //  atropisomer clockwise rotation
        // RDKit✔️✔️:   STEREOATROPCCW,  //  atropisomer counter clockwise rotation
        // RDKit✔️✔️: } BondStereo;
        self as i64
    }

    #[must_use]
    pub const fn rdkit_name(self) -> &'static str {
        // RDKit✔️✔️: python::enum_<Bond::BondStereo>("BondStereo")
        // RDKit✔️✔️:     .value("STEREONONE", Bond::STEREONONE)
        // RDKit✔️✔️:     .value("STEREOANY", Bond::STEREOANY)
        // RDKit✔️✔️:     .value("STEREOZ", Bond::STEREOZ)
        // RDKit✔️✔️:     .value("STEREOE", Bond::STEREOE)
        // RDKit✔️✔️:     .value("STEREOCIS", Bond::STEREOCIS)
        // RDKit✔️✔️:     .value("STEREOTRANS", Bond::STEREOTRANS)
        // RDKit✔️✔️:     .value("STEREOATROPCW", Bond::STEREOATROPCW)
        // RDKit✔️✔️:     .value("STEREOATROPCCW", Bond::STEREOATROPCCW);
        match self {
            Self::None => "STEREONONE",
            Self::Any => "STEREOANY",
            Self::Z => "STEREOZ",
            Self::E => "STEREOE",
            Self::Cis => "STEREOCIS",
            Self::Trans => "STEREOTRANS",
            Self::AtropCw => "STEREOATROPCW",
            Self::AtropCcw => "STEREOATROPCCW",
        }
    }

    #[must_use]
    pub const fn from_rdkit_code(code: i64) -> Option<Self> {
        // RDKit✔️✔️: typedef enum {
        // RDKit✔️✔️:   STEREONONE = 0,  // no special style
        // RDKit✔️✔️:   STEREOANY,       // intentionally unspecified
        // RDKit✔️✔️:   STEREOZ,         // Z double bond
        // RDKit✔️✔️:   STEREOE,         // E double bond
        // RDKit✔️✔️:   STEREOCIS,       // cis double bond
        // RDKit✔️✔️:   STEREOTRANS,     // trans double bond
        // RDKit✔️✔️:   STEREOATROPCW,   //  atropisomer clockwise rotation
        // RDKit✔️✔️:   STEREOATROPCCW,  //  atropisomer counter clockwise rotation
        // RDKit✔️✔️: } BondStereo;
        match code {
            0 => Some(Self::None),
            1 => Some(Self::Any),
            2 => Some(Self::Z),
            3 => Some(Self::E),
            4 => Some(Self::Cis),
            5 => Some(Self::Trans),
            6 => Some(Self::AtropCw),
            7 => Some(Self::AtropCcw),
            _ => None,
        }
    }

    #[must_use]
    pub fn from_rdkit_name(name: &str) -> Option<Self> {
        // This is the checked inverse of the exact Wrap/Bond.cpp table copied
        // into `rdkit_name` above.
        match name {
            "STEREONONE" => Some(Self::None),
            "STEREOANY" => Some(Self::Any),
            "STEREOZ" => Some(Self::Z),
            "STEREOE" => Some(Self::E),
            "STEREOCIS" => Some(Self::Cis),
            "STEREOTRANS" => Some(Self::Trans),
            "STEREOATROPCW" => Some(Self::AtropCw),
            "STEREOATROPCCW" => Some(Self::AtropCcw),
            _ => None,
        }
    }
}

macro_rules! impl_display_as_rdkit_name {
    ($($type:ty),+ $(,)?) => {
        $(
            impl fmt::Display for $type {
                fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
                    formatter.write_str(self.rdkit_name())
                }
            }
        )+
    };
}

impl_display_as_rdkit_name!(
    BondOrder,
    ChiralTag,
    Hybridization,
    BondDirection,
    BondStereo
);
