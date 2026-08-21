use std::{collections::BTreeMap, fmt, str::FromStr};

use serde::{Deserialize, Deserializer, Serialize, Serializer};
use thiserror::Error;

use crate::{AtomQueryPredicate, QueryNode};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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

/// Typed atom-level PDB residue metadata.
///
/// This models the RDKit `AtomPDBResidueInfo` subset needed by hydrogen
/// addition. It is atom state, not generic string props, so topology operations
/// remap it with the atom row.
#[derive(Debug, Clone)]
pub struct AtomPdbResidueInfo {
    atom_name: String,
    serial_number: i32,
    alt_loc: String,
    residue_name: String,
    residue_number: i32,
    chain_id: String,
    insertion_code: String,
    occupancy: f64,
    temp_factor: f64,
    is_hetero_atom: bool,
}

impl PartialEq for AtomPdbResidueInfo {
    fn eq(&self, other: &Self) -> bool {
        self.atom_name == other.atom_name
            && self.serial_number == other.serial_number
            && self.alt_loc == other.alt_loc
            && self.residue_name == other.residue_name
            && self.residue_number == other.residue_number
            && self.chain_id == other.chain_id
            && self.insertion_code == other.insertion_code
            && self.occupancy.to_bits() == other.occupancy.to_bits()
            && self.temp_factor.to_bits() == other.temp_factor.to_bits()
            && self.is_hetero_atom == other.is_hetero_atom
    }
}

impl Eq for AtomPdbResidueInfo {}

impl AtomPdbResidueInfo {
    #[must_use]
    pub fn new(
        atom_name: impl Into<String>,
        serial_number: i32,
        residue_name: impl Into<String>,
        residue_number: i32,
        chain_id: impl Into<String>,
        is_hetero_atom: bool,
    ) -> Self {
        Self {
            atom_name: atom_name.into(),
            serial_number,
            alt_loc: String::new(),
            residue_name: residue_name.into(),
            residue_number,
            chain_id: chain_id.into(),
            insertion_code: String::new(),
            occupancy: 1.0,
            temp_factor: 0.0,
            is_hetero_atom,
        }
    }

    #[must_use]
    pub fn with_alt_loc(mut self, alt_loc: impl Into<String>) -> Self {
        self.alt_loc = alt_loc.into();
        self
    }

    #[must_use]
    pub fn with_insertion_code(mut self, insertion_code: impl Into<String>) -> Self {
        self.insertion_code = insertion_code.into();
        self
    }

    #[must_use]
    pub const fn with_occupancy(mut self, occupancy: f64) -> Self {
        self.occupancy = occupancy;
        self
    }

    #[must_use]
    pub const fn with_temp_factor(mut self, temp_factor: f64) -> Self {
        self.temp_factor = temp_factor;
        self
    }

    #[must_use]
    pub fn atom_name(&self) -> &str {
        &self.atom_name
    }

    #[must_use]
    pub const fn serial_number(&self) -> i32 {
        self.serial_number
    }

    #[must_use]
    pub fn alt_loc(&self) -> &str {
        &self.alt_loc
    }

    #[must_use]
    pub fn residue_name(&self) -> &str {
        &self.residue_name
    }

    #[must_use]
    pub const fn residue_number(&self) -> i32 {
        self.residue_number
    }

    #[must_use]
    pub fn chain_id(&self) -> &str {
        &self.chain_id
    }

    #[must_use]
    pub fn insertion_code(&self) -> &str {
        &self.insertion_code
    }

    #[must_use]
    pub const fn occupancy(&self) -> f64 {
        self.occupancy
    }

    #[must_use]
    pub const fn temp_factor(&self) -> f64 {
        self.temp_factor
    }

    #[must_use]
    pub const fn is_hetero_atom(&self) -> bool {
        self.is_hetero_atom
    }
}

/// Stable atom table index.
///
/// This is intentionally a newtype instead of a bare `usize` so APIs make index
/// spaces explicit. Do not replace this with raw indices in public APIs without
/// human-author approval.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AtomId(usize);

impl AtomId {
    #[must_use]
    pub const fn new(index: usize) -> Self {
        Self(index)
    }

    #[must_use]
    pub const fn index(self) -> usize {
        self.0
    }
}

impl fmt::Display for AtomId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// Stable chemical-element identity.
///
/// The private atomic number is always in the inclusive range `0..=118`.
/// Zero is the source-compatible dummy atom (`*`); `1..=118` are H through
/// Og. Use [`Element::from_atomic_number`] or [`FromStr`] for checked external
/// construction.
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

/// Source-aligned periodic-table metadata for an [`Element`].
///
/// This is a read-only projection of the periodic-table data used by the
/// chemistry core. `valences` preserves the source sentinel values, including
/// `-1`; callers must not reinterpret those values as an exhaustive IUPAC
/// oxidation-state table.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
#[non_exhaustive]
pub struct ElementInfo {
    pub element: Element,
    pub symbol: &'static str,
    pub atomic_number: u8,
    pub period: u8,
    pub outer_electrons: i32,
    pub valences: &'static [i32],
    pub atomic_weight: f64,
}

impl Element {
    pub const DUMMY: Self = Self { atomic_number: 0 };
    // Period 1
    pub const H: Self = Self { atomic_number: 1 };
    pub const HE: Self = Self { atomic_number: 2 };
    // Period 2
    pub const LI: Self = Self { atomic_number: 3 };
    pub const BE: Self = Self { atomic_number: 4 };
    pub const B: Self = Self { atomic_number: 5 };
    pub const C: Self = Self { atomic_number: 6 };
    pub const N: Self = Self { atomic_number: 7 };
    pub const O: Self = Self { atomic_number: 8 };
    pub const F: Self = Self { atomic_number: 9 };
    pub const NE: Self = Self { atomic_number: 10 };
    // Period 3
    pub const NA: Self = Self { atomic_number: 11 };
    pub const MG: Self = Self { atomic_number: 12 };
    pub const AL: Self = Self { atomic_number: 13 };
    pub const SI: Self = Self { atomic_number: 14 };
    pub const P: Self = Self { atomic_number: 15 };
    pub const S: Self = Self { atomic_number: 16 };
    pub const CL: Self = Self { atomic_number: 17 };
    pub const AR: Self = Self { atomic_number: 18 };
    // Period 4
    pub const K: Self = Self { atomic_number: 19 };
    pub const CA: Self = Self { atomic_number: 20 };
    pub const SC: Self = Self { atomic_number: 21 };
    pub const TI: Self = Self { atomic_number: 22 };
    pub const V: Self = Self { atomic_number: 23 };
    pub const CR: Self = Self { atomic_number: 24 };
    pub const MN: Self = Self { atomic_number: 25 };
    pub const FE: Self = Self { atomic_number: 26 };
    pub const CO: Self = Self { atomic_number: 27 };
    pub const NI: Self = Self { atomic_number: 28 };
    pub const CU: Self = Self { atomic_number: 29 };
    pub const ZN: Self = Self { atomic_number: 30 };
    pub const GA: Self = Self { atomic_number: 31 };
    pub const GE: Self = Self { atomic_number: 32 };
    pub const AS: Self = Self { atomic_number: 33 };
    pub const SE: Self = Self { atomic_number: 34 };
    pub const BR: Self = Self { atomic_number: 35 };
    pub const KR: Self = Self { atomic_number: 36 };
    // Period 5
    pub const RB: Self = Self { atomic_number: 37 };
    pub const SR: Self = Self { atomic_number: 38 };
    pub const Y: Self = Self { atomic_number: 39 };
    pub const ZR: Self = Self { atomic_number: 40 };
    pub const NB: Self = Self { atomic_number: 41 };
    pub const MO: Self = Self { atomic_number: 42 };
    pub const TC: Self = Self { atomic_number: 43 };
    pub const RU: Self = Self { atomic_number: 44 };
    pub const RH: Self = Self { atomic_number: 45 };
    pub const PD: Self = Self { atomic_number: 46 };
    pub const AG: Self = Self { atomic_number: 47 };
    pub const CD: Self = Self { atomic_number: 48 };
    pub const IN: Self = Self { atomic_number: 49 };
    pub const SN: Self = Self { atomic_number: 50 };
    pub const SB: Self = Self { atomic_number: 51 };
    pub const TE: Self = Self { atomic_number: 52 };
    pub const I: Self = Self { atomic_number: 53 };
    pub const XE: Self = Self { atomic_number: 54 };
    // Period 6
    pub const CS: Self = Self { atomic_number: 55 };
    pub const BA: Self = Self { atomic_number: 56 };
    // Lanthanides
    pub const LA: Self = Self { atomic_number: 57 };
    pub const CE: Self = Self { atomic_number: 58 };
    pub const PR: Self = Self { atomic_number: 59 };
    pub const ND: Self = Self { atomic_number: 60 };
    pub const PM: Self = Self { atomic_number: 61 };
    pub const SM: Self = Self { atomic_number: 62 };
    pub const EU: Self = Self { atomic_number: 63 };
    pub const GD: Self = Self { atomic_number: 64 };
    pub const TB: Self = Self { atomic_number: 65 };
    pub const DY: Self = Self { atomic_number: 66 };
    pub const HO: Self = Self { atomic_number: 67 };
    pub const ER: Self = Self { atomic_number: 68 };
    pub const TM: Self = Self { atomic_number: 69 };
    pub const YB: Self = Self { atomic_number: 70 };
    pub const LU: Self = Self { atomic_number: 71 };
    // Period 6 continued
    pub const HF: Self = Self { atomic_number: 72 };
    pub const TA: Self = Self { atomic_number: 73 };
    pub const W: Self = Self { atomic_number: 74 };
    pub const RE: Self = Self { atomic_number: 75 };
    pub const OS: Self = Self { atomic_number: 76 };
    pub const IR: Self = Self { atomic_number: 77 };
    pub const PT: Self = Self { atomic_number: 78 };
    pub const AU: Self = Self { atomic_number: 79 };
    pub const HG: Self = Self { atomic_number: 80 };
    pub const TL: Self = Self { atomic_number: 81 };
    pub const PB: Self = Self { atomic_number: 82 };
    pub const BI: Self = Self { atomic_number: 83 };
    pub const PO: Self = Self { atomic_number: 84 };
    pub const AT: Self = Self { atomic_number: 85 };
    pub const RN: Self = Self { atomic_number: 86 };
    // Period 7
    pub const FR: Self = Self { atomic_number: 87 };
    pub const RA: Self = Self { atomic_number: 88 };
    // Actinides
    pub const AC: Self = Self { atomic_number: 89 };
    pub const TH: Self = Self { atomic_number: 90 };
    pub const PA: Self = Self { atomic_number: 91 };
    pub const U: Self = Self { atomic_number: 92 };
    pub const NP: Self = Self { atomic_number: 93 };
    pub const PU: Self = Self { atomic_number: 94 };
    pub const AM: Self = Self { atomic_number: 95 };
    pub const CM: Self = Self { atomic_number: 96 };
    pub const BK: Self = Self { atomic_number: 97 };
    pub const CF: Self = Self { atomic_number: 98 };
    pub const ES: Self = Self { atomic_number: 99 };
    pub const FM: Self = Self { atomic_number: 100 };
    pub const MD: Self = Self { atomic_number: 101 };
    pub const NO: Self = Self { atomic_number: 102 };
    pub const LR: Self = Self { atomic_number: 103 };
    // Period 7 continued
    pub const RF: Self = Self { atomic_number: 104 };
    pub const DB: Self = Self { atomic_number: 105 };
    pub const SG: Self = Self { atomic_number: 106 };
    pub const BH: Self = Self { atomic_number: 107 };
    pub const HS: Self = Self { atomic_number: 108 };
    pub const MT: Self = Self { atomic_number: 109 };
    pub const DS: Self = Self { atomic_number: 110 };
    pub const RG: Self = Self { atomic_number: 111 };
    pub const CN: Self = Self { atomic_number: 112 };
    pub const NH: Self = Self { atomic_number: 113 };
    pub const FL: Self = Self { atomic_number: 114 };
    pub const MC: Self = Self { atomic_number: 115 };
    pub const LV: Self = Self { atomic_number: 116 };
    pub const TS: Self = Self { atomic_number: 117 };
    pub const OG: Self = Self { atomic_number: 118 };

    /// Construct an element from a checked atomic number.
    ///
    /// `0` is the dummy atom and `1..=118` are the real elements. Values above
    /// 118 are rejected rather than creating an invalid element identity.
    #[must_use]
    pub const fn from_atomic_number(atomic_number: u8) -> Option<Self> {
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

    /// Construct an element from a source-recognized symbol.
    ///
    /// Parsing is case-sensitive, matching the source periodic-table API.
    /// Legacy source aliases such as `Uut` and `Uup` are accepted and
    /// canonicalize to `Nh` and `Mc` respectively.
    #[must_use]
    pub fn from_symbol(symbol: &str) -> Option<Self> {
        crate::chemistry::valence::rdkit_atomic_number_from_symbol(symbol)
            .and_then(Self::from_atomic_number)
    }

    /// Return the canonical element symbol (`*`, `H` through `Og`).
    #[must_use]
    pub fn symbol(self) -> &'static str {
        crate::chemistry::valence::rdkit_element_symbol(self.atomic_number)
            .expect("Element's private atomic number is always in 0..=118")
    }

    /// Return the source-aligned periodic-table metadata used by COSMolKit.
    #[must_use]
    pub fn info(self) -> ElementInfo {
        let atomic_number = self.atomic_number;
        ElementInfo {
            element: self,
            symbol: self.symbol(),
            atomic_number,
            period: crate::chemistry::valence::periodic_table_row(atomic_number)
                .expect("Element's private atomic number is always in 0..=118"),
            outer_electrons: crate::chemistry::valence::periodic_table_outer_electrons(
                atomic_number,
            )
            .expect("Element's private atomic number is always in 0..=118"),
            valences: crate::chemistry::valence::rdkit_valence_list(atomic_number)
                .expect("Element's private atomic number is always in 0..=118")
                .expect("every valid Element has a periodic-table row"),
            atomic_weight: crate::chemistry::valence::rdkit_atomic_mass(atomic_number, None),
        }
    }

    /// Iterate over the 118 real elements, excluding the dummy atom.
    pub fn iter() -> impl ExactSizeIterator<Item = Self> + DoubleEndedIterator {
        ELEMENTS.iter().copied()
    }

    /// Iterate over the dummy atom followed by all 118 real elements.
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

/// Atom construction payload.
///
/// `AtomSpec` is deliberately separate from `Atom`: callers provide facts, and
/// builders assign indices. Future agents must not add an `index` field here.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomSpec {
    element: Element,
    formal_charge: i8,
    explicit_hydrogens: u8,
    chiral_tag: ChiralTag,
    chiral_permutation: Option<u32>,
    unknown_stereo: bool,
    mol_parity: Option<i32>,
    mol_inversion_flag: Option<i32>,
    implicit_hydrogen: bool,
    tracked_isotopic_hydrogens: Vec<u16>,
    is_aromatic: bool,
    isotope: Option<u16>,
    atom_map: Option<u32>,
    no_implicit: bool,
    radical_electrons: u8,
    hybridization: Hybridization,
    query: Option<QueryNode<AtomQueryPredicate>>,
    props: BTreeMap<String, String>,
    pdb_residue_info: Option<AtomPdbResidueInfo>,
}

impl AtomSpec {
    #[must_use]
    pub const fn new(element: Element) -> Self {
        Self {
            element,
            formal_charge: 0,
            explicit_hydrogens: 0,
            chiral_tag: ChiralTag::Unspecified,
            chiral_permutation: None,
            unknown_stereo: false,
            mol_parity: None,
            mol_inversion_flag: None,
            implicit_hydrogen: false,
            tracked_isotopic_hydrogens: Vec::new(),
            is_aromatic: false,
            isotope: None,
            atom_map: None,
            no_implicit: false,
            radical_electrons: 0,
            hybridization: Hybridization::Unspecified,
            query: None,
            props: BTreeMap::new(),
            pdb_residue_info: None,
        }
    }

    #[must_use]
    pub const fn with_element(mut self, element: Element) -> Self {
        self.element = element;
        self
    }

    #[must_use]
    pub const fn with_formal_charge(mut self, formal_charge: i8) -> Self {
        self.formal_charge = formal_charge;
        self
    }

    #[must_use]
    pub const fn with_explicit_hydrogens(mut self, explicit_hydrogens: u8) -> Self {
        self.explicit_hydrogens = explicit_hydrogens;
        self
    }

    #[must_use]
    pub const fn with_chiral_tag(mut self, chiral_tag: ChiralTag) -> Self {
        self.chiral_tag = chiral_tag;
        self
    }

    #[must_use]
    pub const fn with_chiral_permutation(mut self, chiral_permutation: u32) -> Self {
        self.chiral_permutation = Some(chiral_permutation);
        self
    }

    #[must_use]
    pub const fn without_chiral_permutation(mut self) -> Self {
        self.chiral_permutation = None;
        self
    }

    #[must_use]
    pub const fn with_unknown_stereo(mut self, unknown_stereo: bool) -> Self {
        self.unknown_stereo = unknown_stereo;
        self
    }

    #[must_use]
    pub const fn with_mol_parity(mut self, mol_parity: i32) -> Self {
        self.mol_parity = Some(mol_parity);
        self
    }

    #[must_use]
    pub const fn without_mol_parity(mut self) -> Self {
        self.mol_parity = None;
        self
    }

    #[must_use]
    pub const fn with_mol_inversion_flag(mut self, mol_inversion_flag: i32) -> Self {
        self.mol_inversion_flag = Some(mol_inversion_flag);
        self
    }

    #[must_use]
    pub const fn without_mol_inversion_flag(mut self) -> Self {
        self.mol_inversion_flag = None;
        self
    }

    #[must_use]
    pub const fn with_implicit_hydrogen(mut self, implicit_hydrogen: bool) -> Self {
        self.implicit_hydrogen = implicit_hydrogen;
        self
    }

    #[must_use]
    pub fn with_tracked_isotopic_hydrogens(mut self, isotopes: Vec<u16>) -> Self {
        self.tracked_isotopic_hydrogens = isotopes;
        self
    }

    #[must_use]
    pub fn without_tracked_isotopic_hydrogens(mut self) -> Self {
        self.tracked_isotopic_hydrogens.clear();
        self
    }

    #[must_use]
    pub const fn with_aromatic(mut self, is_aromatic: bool) -> Self {
        self.is_aromatic = is_aromatic;
        self
    }

    #[must_use]
    pub const fn with_isotope(mut self, isotope: u16) -> Self {
        self.isotope = Some(isotope);
        self
    }

    #[must_use]
    pub const fn without_isotope(mut self) -> Self {
        self.isotope = None;
        self
    }

    #[must_use]
    pub const fn with_atom_map(mut self, atom_map: u32) -> Self {
        self.atom_map = Some(atom_map);
        self
    }

    #[must_use]
    pub const fn without_atom_map(mut self) -> Self {
        self.atom_map = None;
        self
    }

    #[must_use]
    pub const fn with_no_implicit(mut self, no_implicit: bool) -> Self {
        self.no_implicit = no_implicit;
        self
    }

    #[must_use]
    pub const fn with_radical_electrons(mut self, radical_electrons: u8) -> Self {
        self.radical_electrons = radical_electrons;
        self
    }

    #[must_use]
    pub const fn with_hybridization(mut self, hybridization: Hybridization) -> Self {
        self.hybridization = hybridization;
        self
    }

    #[must_use]
    pub fn with_query(mut self, query: QueryNode<AtomQueryPredicate>) -> Self {
        self.query = Some(query);
        self
    }

    #[must_use]
    pub fn with_prop(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        self.props.insert(key.into(), value.into());
        self
    }

    #[must_use]
    pub fn with_pdb_residue_info(mut self, info: AtomPdbResidueInfo) -> Self {
        self.pdb_residue_info = Some(info);
        self
    }

    #[must_use]
    pub fn without_query(mut self) -> Self {
        self.query = None;
        self
    }

    #[must_use]
    pub const fn element(&self) -> Element {
        self.element
    }

    #[must_use]
    pub const fn formal_charge(&self) -> i8 {
        self.formal_charge
    }

    #[must_use]
    pub const fn chiral_tag(&self) -> ChiralTag {
        self.chiral_tag
    }

    #[must_use]
    pub const fn chiral_permutation(&self) -> Option<u32> {
        self.chiral_permutation
    }

    #[must_use]
    pub const fn unknown_stereo(&self) -> bool {
        self.unknown_stereo
    }

    #[must_use]
    pub const fn mol_parity(&self) -> Option<i32> {
        self.mol_parity
    }

    #[must_use]
    pub const fn mol_inversion_flag(&self) -> Option<i32> {
        self.mol_inversion_flag
    }

    #[must_use]
    pub const fn implicit_hydrogen(&self) -> bool {
        self.implicit_hydrogen
    }

    #[must_use]
    pub fn tracked_isotopic_hydrogens(&self) -> &[u16] {
        &self.tracked_isotopic_hydrogens
    }

    #[must_use]
    pub const fn is_aromatic(&self) -> bool {
        self.is_aromatic
    }

    #[must_use]
    pub const fn isotope(&self) -> Option<u16> {
        self.isotope
    }

    #[must_use]
    pub const fn atom_map(&self) -> Option<u32> {
        self.atom_map
    }

    #[must_use]
    pub const fn no_implicit(&self) -> bool {
        self.no_implicit
    }

    #[must_use]
    pub const fn radical_electrons(&self) -> u8 {
        self.radical_electrons
    }

    #[must_use]
    pub const fn query(&self) -> Option<&QueryNode<AtomQueryPredicate>> {
        self.query.as_ref()
    }

    #[must_use]
    pub(crate) fn query_mut(&mut self) -> Option<&mut QueryNode<AtomQueryPredicate>> {
        self.query.as_mut()
    }

    #[must_use]
    pub fn props(&self) -> &BTreeMap<String, String> {
        &self.props
    }

    #[must_use]
    pub fn prop(&self, key: &str) -> Option<&str> {
        self.props.get(key).map(String::as_str)
    }

    #[must_use]
    pub const fn pdb_residue_info(&self) -> Option<&AtomPdbResidueInfo> {
        self.pdb_residue_info.as_ref()
    }
}

/// Immutable atom record owned by `Molecule`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Atom {
    id: AtomId,
    element: Element,
    formal_charge: i8,
    explicit_hydrogens: u8,
    chiral_tag: ChiralTag,
    chiral_permutation: Option<u32>,
    unknown_stereo: bool,
    mol_parity: Option<i32>,
    mol_inversion_flag: Option<i32>,
    implicit_hydrogen: bool,
    tracked_isotopic_hydrogens: Vec<u16>,
    is_aromatic: bool,
    isotope: Option<u16>,
    atom_map: Option<u32>,
    no_implicit: bool,
    radical_electrons: u8,
    hybridization: Hybridization,
    query: Option<QueryNode<AtomQueryPredicate>>,
    props: BTreeMap<String, String>,
    pdb_residue_info: Option<AtomPdbResidueInfo>,
}

impl Atom {
    pub(crate) fn from_spec(id: AtomId, spec: AtomSpec) -> Self {
        Self {
            id,
            element: spec.element,
            formal_charge: spec.formal_charge,
            explicit_hydrogens: spec.explicit_hydrogens,
            chiral_tag: spec.chiral_tag,
            chiral_permutation: spec.chiral_permutation,
            unknown_stereo: spec.unknown_stereo,
            mol_parity: spec.mol_parity,
            mol_inversion_flag: spec.mol_inversion_flag,
            implicit_hydrogen: spec.implicit_hydrogen,
            tracked_isotopic_hydrogens: spec.tracked_isotopic_hydrogens,
            is_aromatic: spec.is_aromatic,
            isotope: spec.isotope,
            atom_map: spec.atom_map,
            no_implicit: spec.no_implicit,
            radical_electrons: spec.radical_electrons,
            hybridization: spec.hybridization,
            query: spec.query,
            props: spec.props,
            pdb_residue_info: spec.pdb_residue_info,
        }
    }

    #[allow(dead_code)]
    pub(crate) fn with_id(mut self, id: AtomId) -> Self {
        self.id = id;
        self
    }

    #[must_use]
    pub const fn id(&self) -> AtomId {
        self.id
    }

    #[must_use]
    pub const fn element(&self) -> Element {
        self.element
    }

    #[must_use]
    pub const fn atomic_number(&self) -> u8 {
        self.element.atomic_number()
    }

    #[must_use]
    pub const fn formal_charge(&self) -> i8 {
        self.formal_charge
    }

    #[must_use]
    pub const fn explicit_hydrogens(&self) -> u8 {
        self.explicit_hydrogens
    }

    #[must_use]
    pub const fn chiral_tag(&self) -> ChiralTag {
        self.chiral_tag
    }

    #[must_use]
    pub const fn chiral_permutation(&self) -> Option<u32> {
        self.chiral_permutation
    }

    #[must_use]
    pub const fn unknown_stereo(&self) -> bool {
        self.unknown_stereo
    }

    #[must_use]
    pub const fn mol_parity(&self) -> Option<i32> {
        self.mol_parity
    }

    #[must_use]
    pub const fn mol_inversion_flag(&self) -> Option<i32> {
        self.mol_inversion_flag
    }

    #[must_use]
    pub const fn implicit_hydrogen(&self) -> bool {
        self.implicit_hydrogen
    }

    #[must_use]
    pub fn tracked_isotopic_hydrogens(&self) -> &[u16] {
        &self.tracked_isotopic_hydrogens
    }

    #[must_use]
    pub const fn is_aromatic(&self) -> bool {
        self.is_aromatic
    }

    #[must_use]
    pub const fn isotope(&self) -> Option<u16> {
        self.isotope
    }

    #[must_use]
    pub const fn atom_map(&self) -> Option<u32> {
        self.atom_map
    }

    #[must_use]
    pub const fn no_implicit(&self) -> bool {
        self.no_implicit
    }

    #[must_use]
    pub const fn radical_electrons(&self) -> u8 {
        self.radical_electrons
    }

    #[must_use]
    pub const fn hybridization(&self) -> Hybridization {
        self.hybridization
    }

    #[must_use]
    pub const fn query(&self) -> Option<&QueryNode<AtomQueryPredicate>> {
        self.query.as_ref()
    }

    #[must_use]
    pub(crate) fn query_mut(&mut self) -> Option<&mut QueryNode<AtomQueryPredicate>> {
        self.query.as_mut()
    }

    #[must_use]
    pub fn props(&self) -> &BTreeMap<String, String> {
        &self.props
    }

    #[must_use]
    pub fn prop(&self, key: &str) -> Option<&str> {
        self.props.get(key).map(String::as_str)
    }

    #[must_use]
    pub const fn pdb_residue_info(&self) -> Option<&AtomPdbResidueInfo> {
        self.pdb_residue_info.as_ref()
    }

    #[allow(dead_code)]
    pub(crate) fn set_chiral_tag(&mut self, chiral_tag: ChiralTag) {
        self.chiral_tag = chiral_tag;
    }

    #[allow(dead_code)]
    pub(crate) fn set_chiral_permutation(&mut self, chiral_permutation: Option<u32>) {
        self.chiral_permutation = chiral_permutation;
    }

    #[allow(dead_code)]
    pub(crate) fn set_unknown_stereo(&mut self, unknown_stereo: bool) {
        self.unknown_stereo = unknown_stereo;
    }

    #[allow(dead_code)]
    pub(crate) fn set_mol_parity(&mut self, mol_parity: Option<i32>) {
        self.mol_parity = mol_parity;
    }

    #[allow(dead_code)]
    pub(crate) fn set_mol_inversion_flag(&mut self, mol_inversion_flag: Option<i32>) {
        self.mol_inversion_flag = mol_inversion_flag;
    }

    #[allow(dead_code)]
    pub(crate) fn set_implicit_hydrogen(&mut self, implicit_hydrogen: bool) {
        self.implicit_hydrogen = implicit_hydrogen;
    }

    #[allow(dead_code)]
    pub(crate) fn set_tracked_isotopic_hydrogens(&mut self, isotopes: Vec<u16>) {
        self.tracked_isotopic_hydrogens = isotopes;
    }

    #[allow(dead_code)]
    pub(crate) fn set_aromatic(&mut self, is_aromatic: bool) {
        self.is_aromatic = is_aromatic;
    }

    #[allow(dead_code)]
    pub(crate) fn set_formal_charge(&mut self, formal_charge: i8) {
        self.formal_charge = formal_charge;
    }

    #[allow(dead_code)]
    pub(crate) fn set_explicit_hydrogens(&mut self, explicit_hydrogens: u8) {
        self.explicit_hydrogens = explicit_hydrogens;
    }

    #[allow(dead_code)]
    pub(crate) fn set_isotope(&mut self, isotope: Option<u16>) {
        self.isotope = isotope;
    }

    #[allow(dead_code)]
    pub(crate) fn set_atom_map(&mut self, atom_map: Option<u32>) {
        self.atom_map = atom_map;
    }

    #[allow(dead_code)]
    pub(crate) fn set_no_implicit(&mut self, no_implicit: bool) {
        self.no_implicit = no_implicit;
    }

    #[allow(dead_code)]
    pub(crate) fn set_radical_electrons(&mut self, radical_electrons: u8) {
        self.radical_electrons = radical_electrons;
    }

    #[allow(dead_code)]
    pub(crate) fn set_query(&mut self, query: Option<QueryNode<AtomQueryPredicate>>) {
        self.query = query;
    }

    #[allow(dead_code)]
    pub(crate) fn set_prop(&mut self, key: impl Into<String>, value: impl Into<String>) {
        self.props.insert(key.into(), value.into());
    }

    #[allow(dead_code)]
    pub(crate) fn clear_prop(&mut self, key: &str) {
        self.props.remove(key);
    }

    #[allow(dead_code)]
    pub(crate) fn set_hybridization(&mut self, hybridization: Hybridization) {
        self.hybridization = hybridization;
    }

    #[allow(dead_code)]
    pub(crate) fn set_pdb_residue_info(&mut self, info: Option<AtomPdbResidueInfo>) {
        self.pdb_residue_info = info;
    }
}

#[cfg(test)]
mod tests {
    use super::{ELEMENTS, ELEMENTS_WITH_DUMMY, Element};

    #[test]
    fn element_domain_and_public_tables_cover_exactly_dummy_through_oganesson() {
        assert_eq!(ELEMENTS.len(), 118);
        assert_eq!(ELEMENTS_WITH_DUMMY.len(), 119);
        assert_eq!(Element::iter().count(), 118);
        assert_eq!(Element::iter_with_dummy().count(), 119);

        for (atomic_number, element) in ELEMENTS_WITH_DUMMY.iter().copied().enumerate() {
            let atomic_number = u8::try_from(atomic_number).unwrap();
            assert_eq!(element.atomic_number(), atomic_number);
            assert_eq!(Element::from_atomic_number(atomic_number), Some(element));
        }
        assert_eq!(ELEMENTS[0], Element::H);
        assert_eq!(ELEMENTS[117], Element::OG);
        assert_eq!(Element::from_atomic_number(119), None);
        assert_eq!(Element::from_atomic_number(u8::MAX), None);
    }

    #[test]
    fn every_element_symbol_roundtrips_through_display_parse_and_serde() {
        for element in Element::iter_with_dummy() {
            let symbol = element.symbol();
            assert_eq!(element.to_string(), symbol);
            assert_eq!(symbol.parse(), Ok(element));
            assert_eq!(Element::from_symbol(symbol), Some(element));

            let json = serde_json::to_string(&element).unwrap();
            assert_eq!(serde_json::from_str::<Element>(&json).unwrap(), element);
            assert_eq!(json, format!(r#""{symbol}""#));
        }
    }

    #[test]
    fn element_parsing_preserves_source_case_and_legacy_alias_behavior() {
        assert_eq!("Uut".parse(), Ok(Element::NH));
        assert_eq!("Uup".parse(), Ok(Element::MC));
        assert_eq!(
            serde_json::to_string(&"Uut".parse::<Element>().unwrap()).unwrap(),
            r#""Nh""#
        );

        for invalid in ["", "cl", "CL", "Carbon", "Xx"] {
            let error = invalid.parse::<Element>().unwrap_err();
            assert_eq!(error.input(), invalid);
        }
    }

    #[test]
    fn element_info_is_self_consistent_for_every_periodic_table_row() {
        for element in Element::iter_with_dummy() {
            let info = element.info();
            assert_eq!(info.element, element);
            assert_eq!(info.atomic_number, element.atomic_number());
            assert_eq!(info.symbol, element.symbol());
            assert!(!info.valences.is_empty());

            let serialized = serde_json::to_value(info).unwrap();
            assert_eq!(serialized["element"], element.symbol());
            assert_eq!(serialized["symbol"], element.symbol());
            assert_eq!(serialized["atomic_number"], element.atomic_number());
        }

        let carbon = Element::C.info();
        assert_eq!(carbon.symbol, "C");
        assert_eq!(carbon.atomic_number, 6);
        assert_eq!(carbon.period, 2);
        assert_eq!(carbon.outer_electrons, 4);
        assert_eq!(carbon.valences, &[4]);
        assert_eq!(carbon.atomic_weight, 12.011);
    }
}
