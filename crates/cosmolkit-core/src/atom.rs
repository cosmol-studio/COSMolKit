use std::{collections::BTreeMap, fmt};

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

/// Minimal element identity for the redesigned core.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Element {
    atomic_number: u8,
}

impl Element {
    pub const DUMMY: Self = Self { atomic_number: 0 };
    pub const H: Self = Self { atomic_number: 1 };
    pub const C: Self = Self { atomic_number: 6 };
    pub const N: Self = Self { atomic_number: 7 };
    pub const O: Self = Self { atomic_number: 8 };

    pub const fn from_atomic_number(atomic_number: u8) -> Option<Self> {
        Some(Self { atomic_number })
    }

    #[must_use]
    pub const fn atomic_number(self) -> u8 {
        self.atomic_number
    }
}

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
    is_aromatic: bool,
    isotope: Option<u16>,
    atom_map: Option<u32>,
    no_implicit: bool,
    radical_electrons: u8,
    hybridization: Hybridization,
    query: Option<QueryNode<AtomQueryPredicate>>,
    props: BTreeMap<String, String>,
}

impl AtomSpec {
    #[must_use]
    pub const fn new(element: Element) -> Self {
        Self {
            element,
            formal_charge: 0,
            explicit_hydrogens: 0,
            chiral_tag: ChiralTag::Unspecified,
            is_aromatic: false,
            isotope: None,
            atom_map: None,
            no_implicit: false,
            radical_electrons: 0,
            hybridization: Hybridization::Unspecified,
            query: None,
            props: BTreeMap::new(),
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
    pub fn props(&self) -> &BTreeMap<String, String> {
        &self.props
    }

    #[must_use]
    pub fn prop(&self, key: &str) -> Option<&str> {
        self.props.get(key).map(String::as_str)
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
    is_aromatic: bool,
    isotope: Option<u16>,
    atom_map: Option<u32>,
    no_implicit: bool,
    radical_electrons: u8,
    hybridization: Hybridization,
    query: Option<QueryNode<AtomQueryPredicate>>,
    props: BTreeMap<String, String>,
}

impl Atom {
    pub(crate) fn from_spec(id: AtomId, spec: AtomSpec) -> Self {
        Self {
            id,
            element: spec.element,
            formal_charge: spec.formal_charge,
            explicit_hydrogens: spec.explicit_hydrogens,
            chiral_tag: spec.chiral_tag,
            is_aromatic: spec.is_aromatic,
            isotope: spec.isotope,
            atom_map: spec.atom_map,
            no_implicit: spec.no_implicit,
            radical_electrons: spec.radical_electrons,
            hybridization: spec.hybridization,
            query: spec.query,
            props: spec.props,
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
    pub fn props(&self) -> &BTreeMap<String, String> {
        &self.props
    }

    #[must_use]
    pub fn prop(&self, key: &str) -> Option<&str> {
        self.props.get(key).map(String::as_str)
    }

    #[allow(dead_code)]
    pub(crate) fn set_chiral_tag(&mut self, chiral_tag: ChiralTag) {
        self.chiral_tag = chiral_tag;
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
}
