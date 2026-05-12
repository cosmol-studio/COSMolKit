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
