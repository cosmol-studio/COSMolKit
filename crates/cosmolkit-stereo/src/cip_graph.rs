// Source-backed modern RDKit CIP digraph and sequence-rule machinery.
//
// The assignment/configuration layer consumes this crate-private unit; this
// module neither accepts a live molecule nor owns runtime state.

use std::collections::VecDeque;

use cosmolkit_core::{
    AtropisomerError, DoubleBondStereoError, KekulizeAttempt, KekulizeError, KekulizeParams,
    PeriodicTableError, RingFindingError, RingInfo, ValenceAssignment, ValenceError, ValenceModel,
    assign_valence_for_topology, atomic_mass, fast_find_rings, kekulize_if_possible,
};
use cosmolkit_model::{
    Atom, AtomPropertyError, Bond, BondId, BondValueError, CipDescriptor, MoleculePropertyError,
    TopologyBlock, TopologyValidationError,
};
use cosmolkit_types::{BondOrder, ChiralTag, Element};

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum CipLabelerError {
    #[error("CIPLabeler atom index {index} is out of range for {atom_count} atoms")]
    AtomIndexOutOfRange { index: usize, atom_count: usize },
    #[error("CIPLabeler bond index {index} is out of range for {bond_count} bonds")]
    BondIndexOutOfRange { index: usize, bond_count: usize },
    #[error("CIPLabeler {kind} index {index} exceeds the source unsigned-int width")]
    SourceIndexWidthExceeded { kind: &'static str, index: usize },
    #[error("CIPLabeler bond {bond} is not incident to atom {atom}")]
    BondNotIncident { bond: usize, atom: usize },
    #[error("CIPLabeler non integer-order bond is not allowed: {order:?}")]
    NonIntegerBondOrder { order: BondOrder },
    #[error("CIPLabeler node {node} is not an endpoint of edge {edge}")]
    EdgeEndpointMismatch { edge: usize, node: usize },
    #[error("Digraph generation failed: more than {limit} nodes found.")]
    TooManyNodes { limit: usize },
    #[error("Max Iterations Exceeded in CIP label calculation")]
    MaxIterationsExceeded,
    #[error("CIPLabeler unexpected up-edge ordering")]
    UnexpectedUpEdgeOrdering,
    #[error("No sequence rule provided")]
    NoSequenceRuleProvided,
    #[error("Descriptor lists should be the same length!")]
    DescriptorListLengthMismatch,
    #[error("Invalid stereo descriptor")]
    InvalidStereoDescriptor,
    #[error("Substituents should be topologically equivalent!")]
    SubstituentsShouldBeTopologicallyEquivalent,
    #[error("Something unexpected!")]
    SomethingUnexpected,
    #[error("Rule4b instance not in rule set")]
    Rule4bInstanceNotInRuleSet,
    #[error("Rule5New instance not in rule set")]
    Rule5NewInstanceNotInRuleSet,
    #[error("invalid CIPLabeler internal state: {detail}")]
    InvalidInternalState { detail: &'static str },
    #[error("Parity vectors must have size 4.")]
    ParityVectorsMustHaveSize4,
    #[error("CIPLabeler Configuration requires at least one focus atom")]
    EmptyConfigurationFoci,
    #[error("CIPLabeler Tetrahedral received a bad config")]
    BadTetrahedralConfig,
    #[error("CIPLabeler Tetrahedral configuration must have 4 carriers")]
    TetrahedralConfigurationMustHave4Carriers,
    #[error("Received a Descriptor that is not supported for atoms")]
    DescriptorNotSupportedForAtoms,
    #[error("Received an invalid Atom Descriptor")]
    InvalidAtomDescriptor,
    #[error("Could not calculate parity! Carrier mismatch")]
    CarrierMismatch,
    #[error("CIPLabeler Sp2Bond received bad foci")]
    BadSp2BondFoci,
    #[error("CIPLabeler Sp2Bond received bad config")]
    BadSp2BondConfig,
    #[error("CIPLabeler Sp2Bond has incorrect number of stereo atoms")]
    IncorrectNumberOfStereoAtoms,
    #[error("Received a Descriptor that is not supported for double bonds")]
    DescriptorNotSupportedForDoubleBonds,
    #[error("Received an invalid Bond Descriptor")]
    InvalidBondDescriptor,
    #[error("CIPLabeler AtropisomerBond received bad foci")]
    BadAtropisomerBondFoci,
    #[error("CIPLabeler AtropisomerBond received bad config")]
    BadAtropisomerBondConfig,
    #[error("Received a Descriptor that is not supported for atropisomer bonds")]
    DescriptorNotSupportedForAtropisomerBonds,
    #[error("unsupported CIP configuration on atom {atom}: {tag:?}")]
    UnsupportedConfiguration { atom: usize, tag: ChiralTag },
    #[error("invalid detached topology: {0}")]
    InvalidTopology(#[from] TopologyValidationError),
    #[error("failed to update CIP atom {atom} property: {source}")]
    AtomProperty {
        atom: usize,
        source: AtomPropertyError,
    },
    #[error("failed to update CIP bond {bond}: {source}")]
    BondValue { bond: usize, source: BondValueError },
    #[error("failed to update CIP molecule property: {0}")]
    MoleculeProperty(#[from] MoleculePropertyError),
    #[error(transparent)]
    Atropisomer(#[from] AtropisomerError),
    #[error(transparent)]
    DoubleStereo(#[from] DoubleBondStereoError),
    #[error(transparent)]
    RingFinding(#[from] RingFindingError),
    #[error(transparent)]
    Valence(#[from] ValenceError),
    #[error(transparent)]
    Kekulize(#[from] KekulizeError),
    #[error(transparent)]
    PeriodicTable(#[from] PeriodicTableError),
}

trait CipTopologyView {
    fn num_atoms(&self) -> usize;
    fn num_bonds(&self) -> usize;
    fn atoms(&self) -> &[Atom];
    fn bonds(&self) -> &[Bond];
    fn topology_block(&self) -> &TopologyBlock;
}

impl CipTopologyView for TopologyBlock {
    fn num_atoms(&self) -> usize {
        self.atoms.len()
    }
    fn num_bonds(&self) -> usize {
        self.bonds.len()
    }
    fn atoms(&self) -> &[Atom] {
        &self.atoms
    }
    fn bonds(&self) -> &[Bond] {
        &self.bonds
    }
    fn topology_block(&self) -> &TopologyBlock {
        self
    }
}

fn rdkit_atomic_mass(atomic_number: u8, isotope: Option<u16>) -> Result<f64, CipLabelerError> {
    // BEGIN RDKIT CPP FUNCTION Atom::getMass
    // RDKit✔️✔️: double Atom::getMass() const {
    // RDKit✔️✔️:   if (d_isotope) {
    // RDKit✔️✔️:     double res =
    // RDKit✔️✔️:         PeriodicTable::getTable()->getMassForIsotope(d_atomicNum, d_isotope);
    // RDKit✔️✔️:     if (d_atomicNum != 0 && res == 0.0) {
    // RDKit✔️✔️:       res = d_isotope;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return res;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return PeriodicTable::getTable()->getAtomicWeight(d_atomicNum);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Atom::getMass
    let element = Element::from_atomic_number(atomic_number).ok_or(
        CipLabelerError::InvalidInternalState {
            detail: "atomic number exceeds periodic table",
        },
    )?;
    match atomic_mass(element, isotope) {
        Ok(mass) => Ok(mass),
        Err(PeriodicTableError::UnknownIsotope { isotope, .. }) if atomic_number != 0 => {
            Ok(f64::from(isotope))
        }
        Err(PeriodicTableError::UnknownIsotope { .. }) => Ok(0.0),
        Err(error) => Err(error.into()),
    }
}

type CipSourceIndex = u32;
type CipPairing = u64;

const CIP_NO_ATOM: CipSourceIndex = CipSourceIndex::MAX;

// BEGIN RDKIT CPP ENUM Descriptor (CIPLabeler/Descriptor.h)
// RDKit✔️✔️: enum class Descriptor {
// RDKit✔️✔️:   NONE,  // Unspecified
// RDKit✔️✔️:   UNKNOWN,
// RDKit✔️✔️:   ns,  // Other
// RDKit✔️✔️:   /**
// RDKit✔️✔️:    * Tetrahedral
// RDKit✔️✔️:    */
// RDKit✔️✔️:   R,
// RDKit✔️✔️:   S,
// RDKit✔️✔️:   r,
// RDKit✔️✔️:   s,
// RDKit✔️✔️:   /**
// RDKit✔️✔️:    * Cis/Trans
// RDKit✔️✔️:    */
// RDKit✔️✔️:   seqTrans,
// RDKit✔️✔️:   seqCis,
// RDKit✔️✔️:   E,
// RDKit✔️✔️:   Z,
// RDKit✔️✔️:   /* Axial */
// RDKit✔️✔️:   M,
// RDKit✔️✔️:   P,
// RDKit✔️✔️:   m,
// RDKit✔️✔️:   p,
// RDKit✔️✔️:
// RDKit✔️✔️:   SP_4,
// RDKit✔️✔️:   TBPY_5,
// RDKit✔️✔️:   OC_6
// RDKit✔️✔️: };
// END RDKIT CPP ENUM Descriptor
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[allow(non_camel_case_types)]
pub(crate) enum Descriptor {
    None,
    Unknown,
    ns,
    R,
    S,
    r,
    s,
    seqTrans,
    seqCis,
    E,
    Z,
    M,
    P,
    m,
    p,
    SP_4,
    TBPY_5,
    OC_6,
}

impl Descriptor {
    pub(crate) const ALL_IN_RDKIT_ORDER: [Self; 18] = [
        Self::None,
        Self::Unknown,
        Self::ns,
        Self::R,
        Self::S,
        Self::r,
        Self::s,
        Self::seqTrans,
        Self::seqCis,
        Self::E,
        Self::Z,
        Self::M,
        Self::P,
        Self::m,
        Self::p,
        Self::SP_4,
        Self::TBPY_5,
        Self::OC_6,
    ];
}

// BEGIN RDKIT CPP FUNCTION to_string (CIPLabeler/Descriptor.h)
// RDKit✔️✔️: static std::string to_string(const Descriptor &desc) {
// RDKit✔️✔️:   switch (desc) {
// RDKit✔️✔️:     case Descriptor::NONE:
// RDKit✔️✔️:       return "NONE";
// RDKit✔️✔️:     case Descriptor::UNKNOWN:
// RDKit✔️✔️:       return "UNKNOWN";
// RDKit✔️✔️:     case Descriptor::ns:
// RDKit✔️✔️:       return "ns";
// RDKit✔️✔️:     case Descriptor::R:
// RDKit✔️✔️:       return "R";
// RDKit✔️✔️:     case Descriptor::S:
// RDKit✔️✔️:       return "S";
// RDKit✔️✔️:     case Descriptor::r:
// RDKit✔️✔️:       return "r";
// RDKit✔️✔️:     case Descriptor::s:
// RDKit✔️✔️:       return "s";
// RDKit✔️✔️:     case Descriptor::seqTrans:
// RDKit✔️✔️:       return "e";
// RDKit✔️✔️:     case Descriptor::seqCis:
// RDKit✔️✔️:       return "z";
// RDKit✔️✔️:     case Descriptor::E:
// RDKit✔️✔️:       return "E";
// RDKit✔️✔️:     case Descriptor::Z:
// RDKit✔️✔️:       return "Z";
// RDKit✔️✔️:     case Descriptor::M:
// RDKit✔️✔️:       return "M";
// RDKit✔️✔️:     case Descriptor::P:
// RDKit✔️✔️:       return "P";
// RDKit✔️✔️:     case Descriptor::m:
// RDKit✔️✔️:       return "m";
// RDKit✔️✔️:     case Descriptor::p:
// RDKit✔️✔️:       return "p";
// RDKit✔️✔️:     case Descriptor::SP_4:
// RDKit✔️✔️:       return "SP_4";
// RDKit✔️✔️:     case Descriptor::TBPY_5:
// RDKit✔️✔️:       return "TBPY_5";
// RDKit✔️✔️:     case Descriptor::OC_6:
// RDKit✔️✔️:       return "OC_6";
// RDKit✔️✔️:   }
// RDKit✔️✔️:   throw std::runtime_error("Unknown descriptor");
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION to_string
pub(crate) fn descriptor_to_string(desc: Descriptor) -> &'static str {
    if let Some(public_descriptor) = public_descriptor_from_source(desc) {
        return public_descriptor.as_str();
    }
    match desc {
        Descriptor::None => "NONE",
        Descriptor::Unknown => "UNKNOWN",
        Descriptor::ns => "ns",
        Descriptor::SP_4 => "SP_4",
        Descriptor::TBPY_5 => "TBPY_5",
        Descriptor::OC_6 => "OC_6",
        Descriptor::R
        | Descriptor::S
        | Descriptor::r
        | Descriptor::s
        | Descriptor::seqTrans
        | Descriptor::seqCis
        | Descriptor::E
        | Descriptor::Z
        | Descriptor::M
        | Descriptor::P
        | Descriptor::m
        | Descriptor::p => unreachable!("public emitted descriptors returned above"),
    }
}

fn public_descriptor_from_source(desc: Descriptor) -> Option<CipDescriptor> {
    match desc {
        Descriptor::R => Some(CipDescriptor::R),
        Descriptor::S => Some(CipDescriptor::S),
        Descriptor::r => Some(CipDescriptor::LowerR),
        Descriptor::s => Some(CipDescriptor::LowerS),
        Descriptor::seqTrans => Some(CipDescriptor::LowerE),
        Descriptor::seqCis => Some(CipDescriptor::LowerZ),
        Descriptor::E => Some(CipDescriptor::E),
        Descriptor::Z => Some(CipDescriptor::Z),
        Descriptor::M => Some(CipDescriptor::M),
        Descriptor::P => Some(CipDescriptor::P),
        Descriptor::m => Some(CipDescriptor::LowerM),
        Descriptor::p => Some(CipDescriptor::LowerP),
        _ => None,
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) struct RationalI32 {
    numerator: i32,
    denominator: i32,
}

impl PartialOrd for RationalI32 {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for RationalI32 {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // boost::rational compares normalized values, not numerator/denominator
        // tuples.  The products fit in i64 for the source int32 domain.
        (i64::from(self.numerator) * i64::from(other.denominator))
            .cmp(&(i64::from(other.numerator) * i64::from(self.denominator)))
    }
}

impl RationalI32 {
    fn new(numerator: i32, denominator: i32) -> Self {
        debug_assert_ne!(
            denominator, 0,
            "boost::rational does not accept 0 as denominator"
        );
        let mut numerator = numerator;
        let mut denominator = denominator;
        if denominator < 0 {
            numerator = -numerator;
            denominator = -denominator;
        }
        let divisor = gcd_i32(numerator, denominator);
        Self {
            numerator: numerator / divisor,
            denominator: denominator / divisor,
        }
    }

    fn assign(&mut self, numerator: i32, denominator: i32) {
        *self = Self::new(numerator, denominator);
    }

    #[cfg(test)]
    fn tuple(self) -> (i32, i32) {
        (self.numerator, self.denominator)
    }
}

fn gcd_i32(a: i32, b: i32) -> i32 {
    let mut a = i64::from(a).abs();
    let mut b = i64::from(b).abs();
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    if a == 0 { 1 } else { a as i32 }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MancudeType {
    Cv4D3,
    Nv3D2,
    Nv4D3Plus,
    Nv2D2Minus,
    Cv3D3Minus,
    Ov3D2Plus,
    Other,
}

pub(crate) struct CipMol<'a> {
    molecule: &'a TopologyBlock,
    rings: Option<RingInfo>,
    kekulized_bond_orders: Option<Vec<BondOrder>>,
    fractional_atomic_numbers: Option<Vec<RationalI32>>,
    valence: Option<ValenceAssignment>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct CipNodeId(usize);

impl CipNodeId {
    fn new(index: usize) -> Self {
        Self(index)
    }

    fn index(self) -> usize {
        self.0
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct CipEdgeId(usize);

impl CipEdgeId {
    fn new(index: usize) -> Self {
        Self(index)
    }

    fn index(self) -> usize {
        self.0
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CipEdge {
    beg: CipNodeId,
    end: CipNodeId,
    bond_idx: Option<usize>,
    aux: Descriptor,
}

impl CipEdge {
    // BEGIN RDKIT CPP FUNCTION Edge::Edge (CIPLabeler/Edge.cpp)
    // RDKit✔️✔️: Edge::Edge(Node *beg, Node *end, Bond *bond)
    // RDKit✔️✔️:     : dp_beg{beg}, dp_end{end}, dp_bond{bond} {}
    // END RDKIT CPP FUNCTION Edge::Edge
    pub(crate) fn new(beg: CipNodeId, end: CipNodeId, bond_idx: Option<usize>) -> Self {
        Self {
            beg,
            end,
            bond_idx,
            aux: Descriptor::None,
        }
    }

    // BEGIN RDKIT CPP FUNCTION Edge::getOther (CIPLabeler/Edge.cpp)
    // RDKit✔️✔️: Node *Edge::getOther(const Node *node) const {
    // RDKit✔️✔️:   PRECONDITION(node, "bad node")
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (isBeg(node)) {
    // RDKit✔️✔️:     return getEnd();
    // RDKit✔️✔️:   } else if (isEnd(node)) {
    // RDKit✔️✔️:     return getBeg();
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     throw std::runtime_error("Not an end-point of this edge!");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Edge::getOther
    pub(crate) fn get_other(
        &self,
        self_id: CipEdgeId,
        node: CipNodeId,
    ) -> Result<CipNodeId, CipLabelerError> {
        if self.is_beg(node) {
            Ok(self.get_end())
        } else if self.is_end(node) {
            Ok(self.get_beg())
        } else {
            Err(CipLabelerError::EdgeEndpointMismatch {
                edge: self_id.index(),
                node: node.index(),
            })
        }
    }

    // BEGIN RDKIT CPP FUNCTION Edge::getBeg (CIPLabeler/Edge.cpp)
    // RDKit✔️✔️: Node *Edge::getBeg() const { return dp_beg; }
    // END RDKIT CPP FUNCTION Edge::getBeg
    pub(crate) fn get_beg(&self) -> CipNodeId {
        self.beg
    }

    // BEGIN RDKIT CPP FUNCTION Edge::getEnd (CIPLabeler/Edge.cpp)
    // RDKit✔️✔️: Node *Edge::getEnd() const { return dp_end; }
    // END RDKIT CPP FUNCTION Edge::getEnd
    pub(crate) fn get_end(&self) -> CipNodeId {
        self.end
    }

    // BEGIN RDKIT CPP FUNCTION Edge::getBond (CIPLabeler/Edge.cpp)
    // RDKit✔️✔️: Bond *Edge::getBond() const { return dp_bond; }
    // END RDKIT CPP FUNCTION Edge::getBond
    pub(crate) fn get_bond_idx(&self) -> Option<usize> {
        self.bond_idx
    }

    // BEGIN RDKIT CPP FUNCTION Edge::getAux (CIPLabeler/Edge.cpp)
    // RDKit✔️✔️: Descriptor Edge::getAux() const { return d_aux; }
    // END RDKIT CPP FUNCTION Edge::getAux
    pub(crate) fn get_aux(&self) -> Descriptor {
        self.aux
    }

    // BEGIN RDKIT CPP FUNCTION Edge::isBeg (CIPLabeler/Edge.cpp)
    // RDKit✔️✔️: bool Edge::isBeg(const Node *node) const { return node == dp_beg; }
    // END RDKIT CPP FUNCTION Edge::isBeg
    pub(crate) fn is_beg(&self, node: CipNodeId) -> bool {
        node == self.beg
    }

    // BEGIN RDKIT CPP FUNCTION Edge::isEnd (CIPLabeler/Edge.cpp)
    // RDKit✔️✔️: bool Edge::isEnd(const Node *node) const { return node == dp_end; }
    // END RDKIT CPP FUNCTION Edge::isEnd
    pub(crate) fn is_end(&self, node: CipNodeId) -> bool {
        node == self.end
    }

    // BEGIN RDKIT CPP FUNCTION Edge::setAux (CIPLabeler/Edge.cpp)
    // RDKit✔️✔️: void Edge::setAux(Descriptor aux) { d_aux = std::move(aux); }
    // END RDKIT CPP FUNCTION Edge::setAux
    pub(crate) fn set_aux(&mut self, aux: Descriptor) {
        self.aux = aux;
    }

    // BEGIN RDKIT CPP FUNCTION Edge::flip (CIPLabeler/Edge.cpp)
    // RDKit✔️✔️: void Edge::flip() { std::swap(dp_beg, dp_end); }
    // END RDKIT CPP FUNCTION Edge::flip
    pub(crate) fn flip(&mut self) {
        std::mem::swap(&mut self.beg, &mut self.end);
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CipNode {
    digraph: usize,
    atom_idx: Option<usize>,
    distance: i32,
    atomic_num_fraction: RationalI32,
    atomic_mass: f64,
    aux: Descriptor,
    flags: i32,
    edges: Vec<CipEdgeId>,
    visit: Vec<i8>,
}

pub(crate) struct CipDigraph<'a> {
    mol: CipMol<'a>,
    origin: CipNodeId,
    root: CipNodeId,
    rule6_ref: Option<usize>,
    atropisomer_mode: bool,
    nodes: Vec<CipNode>,
    edges: Vec<CipEdge>,
}

pub(crate) struct CipLabelerContext {
    remaining_call_count: u32,
}

impl CipLabelerContext {
    pub(crate) const CONSTITUTIONAL_RULE_TIMEOUT: u32 = 2_000;

    pub(crate) fn new(max_recursive_iterations: u32) -> Self {
        Self {
            remaining_call_count: if max_recursive_iterations == 0 {
                u32::MAX
            } else {
                max_recursive_iterations
            },
        }
    }

    pub(crate) fn with_remaining_call_count(remaining_call_count: u32) -> Self {
        Self {
            remaining_call_count,
        }
    }

    // BEGIN RDKIT CPP FUNCTION decrementRemainingCallCountAndCheck (CIPLabeler.cpp)
    // RDKit✔️✔️: bool decrementRemainingCallCountAndCheck() {
    // RDKit✔️✔️:   return (--CIPLabeler::remainingCallCount) > 0;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION decrementRemainingCallCountAndCheck
    fn decrement_remaining_call_count_and_check(&mut self) -> bool {
        self.remaining_call_count = self.remaining_call_count.wrapping_sub(1);
        self.remaining_call_count > 0
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct CipPriority {
    unique: bool,
    pseudo_asym: bool,
}

impl CipPriority {
    // BEGIN RDKIT CPP CLASS Priority (CIPLabeler/Priority.h)
    // RDKit✔️✔️: class Priority {
    // RDKit✔️✔️:  public:
    // RDKit✔️✔️:   Priority() = delete;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   Priority(bool unique, bool pseudoAsym)
    // RDKit✔️✔️:       : d_unique{unique}, d_pseudoAsym{pseudoAsym} {}
    // RDKit✔️✔️:
    // RDKit✔️✔️:   bool isUnique() const { return d_unique; }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   bool isPseudoAsymetric() const { return d_pseudoAsym; }
    // RDKit✔️✔️:
    // RDKit✔️✔️:  private:
    // RDKit✔️✔️:   bool d_unique;
    // RDKit✔️✔️:   bool d_pseudoAsym;
    // RDKit✔️✔️: };
    // END RDKIT CPP CLASS Priority
    pub(crate) fn new(unique: bool, pseudo_asym: bool) -> Self {
        Self {
            unique,
            pseudo_asym,
        }
    }

    pub(crate) fn is_unique(self) -> bool {
        self.unique
    }

    pub(crate) fn is_pseudo_asymetric(self) -> bool {
        self.pseudo_asym
    }
}

pub(crate) trait CipSequenceRule {
    // BEGIN RDKIT CPP FUNCTION SequenceRule::getBondLabel (rules/SequenceRule.cpp)
    // RDKit✔️✔️: Descriptor SequenceRule::getBondLabel(const Edge *edge) const {
    // RDKit✔️✔️:   Bond *bond = edge->getBond();
    // RDKit✔️✔️:   if (bond == nullptr) {
    // RDKit✔️✔️:     return Descriptor::NONE;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   Descriptor label = edge->getAux();
    // RDKit✔️✔️:   if (label != Descriptor::NONE) {
    // RDKit✔️✔️:     return label;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return label;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION SequenceRule::getBondLabel
    fn get_bond_label(&self, edge: &CipEdge) -> Descriptor {
        if edge.get_bond_idx().is_none() {
            return Descriptor::None;
        }
        edge.get_aux()
    }

    // BEGIN RDKIT CPP FUNCTION SequenceRule::getComparision (rules/SequenceRule.cpp)
    // RDKit✔️✔️: int SequenceRule::getComparision(const Edge *a, const Edge *b) const {
    // RDKit✔️✔️:   return getComparision(a, b, true);
    // RDKit✔️✔️: }
    // RDKit✔️✔️:
    // RDKit✔️✔️: int SequenceRule::getComparision(const Edge *a, const Edge *b,
    // RDKit✔️✔️:                                  bool deep) const {
    // RDKit✔️✔️:   return deep ? recursiveCompare(a, b) : compare(a, b);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION SequenceRule::getComparision
    fn get_comparison(
        &self,
        digraph: &mut CipDigraph<'_>,
        context: &mut CipLabelerContext,
        a: CipEdgeId,
        b: CipEdgeId,
        deep: bool,
    ) -> Result<i32, CipLabelerError> {
        self.get_comparison_with_sort_rules(None, digraph, context, a, b, deep)
    }

    fn get_comparison_with_sort_rules(
        &self,
        sort_rules: Option<&[&dyn CipSequenceRule]>,
        digraph: &mut CipDigraph<'_>,
        context: &mut CipLabelerContext,
        a: CipEdgeId,
        b: CipEdgeId,
        deep: bool,
    ) -> Result<i32, CipLabelerError> {
        if deep {
            recursive_compare_sequence_rule(self, sort_rules, digraph, context, a, b)
        } else {
            self.compare_with_sort_rules(sort_rules, digraph, context, a, b)
        }
    }

    fn compare(
        &self,
        digraph: &mut CipDigraph<'_>,
        context: &mut CipLabelerContext,
        a: CipEdgeId,
        b: CipEdgeId,
    ) -> Result<i32, CipLabelerError>;

    fn compare_with_sort_rules(
        &self,
        _sort_rules: Option<&[&dyn CipSequenceRule]>,
        digraph: &mut CipDigraph<'_>,
        context: &mut CipLabelerContext,
        a: CipEdgeId,
        b: CipEdgeId,
    ) -> Result<i32, CipLabelerError> {
        self.compare(digraph, context, a, b)
    }

    // BEGIN RDKIT CPP FUNCTION SequenceRule::recursiveCompare (rules/SequenceRule.cpp)
    // RDKit✔️✔️: int SequenceRule::recursiveCompare(const Edge *a, const Edge *b) const {
    // RDKit✔️✔️:   if (!CIPLabeler_detail::decrementRemainingCallCountAndCheck()) {
    // RDKit✔️✔️:     throw MaxIterationsExceeded();
    // RDKit✔️✔️:   }
    // RDKit❌❌:   if (ControlCHandler::getGotSignal()) {
    // RDKit❌❌:     throw ControlCCaught();
    // RDKit❌❌:   }
    // COSMolKit has no source-equivalent process-global ControlCHandler.
    // Explicit requests for that integration fail before mutation through
    // CipLabelerError::CancellationUnsupported.
    // RDKit✔️✔️:
    // RDKit✔️✔️:   int cmp = compare(a, b);
    // RDKit✔️✔️:   if (cmp != 0) {
    // RDKit✔️✔️:     return cmp;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto aQueue = std::vector<const Edge *>({a});
    // RDKit✔️✔️:   auto bQueue = std::vector<const Edge *>({b});
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (auto pos = 0u; pos < aQueue.size() && pos < bQueue.size(); ++pos) {
    // RDKit✔️✔️:     a = aQueue[pos];
    // RDKit✔️✔️:     b = bQueue[pos];
    // RDKit✔️✔️:     auto as = a->getEnd()->getEdges();
    // RDKit✔️✔️:     auto bs = b->getEnd()->getEdges();
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // shallow sort first of all
    // RDKit✔️✔️:     sort(a->getEnd(), as, false);
    // RDKit✔️✔️:     sort(b->getEnd(), bs, false);
    // RDKit✔️✔️:
    // RDKit✔️✔️:     int sizediff = three_way_comparison(static_cast<int>(as.size()),
    // RDKit✔️✔️:                                         static_cast<int>(bs.size()));
    // RDKit✔️✔️:
    // RDKit✔️✔️:     {
    // RDKit✔️✔️:       auto aIt = as.begin();
    // RDKit✔️✔️:       auto bIt = bs.begin();
    // RDKit✔️✔️:       for (; aIt != as.end() && bIt != bs.end(); ++aIt, ++bIt) {
    // RDKit✔️✔️:         Node *aNode = a->getEnd();
    // RDKit✔️✔️:         Node *bNode = b->getEnd();
    // RDKit✔️✔️:         Edge *aEdge = *aIt;
    // RDKit✔️✔️:         Edge *bEdge = *bIt;
    // RDKit✔️✔️:
    // RDKit✔️✔️:         if (areUpEdges(aNode, bNode, aEdge, bEdge)) {
    // RDKit✔️✔️:           continue;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:
    // RDKit✔️✔️:         cmp = compare(aEdge, bEdge);
    // RDKit✔️✔️:         if (cmp != 0) {
    // RDKit✔️✔️:           return cmp;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (sizediff != 0) {
    // RDKit✔️✔️:       return sizediff;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     sort(a->getEnd(), as);
    // RDKit✔️✔️:     sort(b->getEnd(), bs);
    // RDKit✔️✔️:
    // RDKit✔️✔️:     {
    // RDKit✔️✔️:       auto aIt = as.begin();
    // RDKit✔️✔️:       auto bIt = bs.begin();
    // RDKit✔️✔️:       for (; aIt != as.end() && bIt != bs.end(); ++aIt, ++bIt) {
    // RDKit✔️✔️:         Node *aNode = a->getEnd();
    // RDKit✔️✔️:         Node *bNode = b->getEnd();
    // RDKit✔️✔️:         Edge *aEdge = *aIt;
    // RDKit✔️✔️:         Edge *bEdge = *bIt;
    // RDKit✔️✔️:
    // RDKit✔️✔️:         if (areUpEdges(aNode, bNode, aEdge, bEdge)) {
    // RDKit✔️✔️:           continue;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:
    // RDKit✔️✔️:         cmp = compare(aEdge, bEdge);
    // RDKit✔️✔️:         if (cmp != 0) {
    // RDKit✔️✔️:           return cmp;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:
    // RDKit✔️✔️:         aQueue.push_back(aEdge);
    // RDKit✔️✔️:         bQueue.push_back(bEdge);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return 0;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION SequenceRule::recursiveCompare
    fn recursive_compare(
        &self,
        digraph: &mut CipDigraph<'_>,
        context: &mut CipLabelerContext,
        a: CipEdgeId,
        b: CipEdgeId,
    ) -> Result<i32, CipLabelerError> {
        recursive_compare_sequence_rule(self, None, digraph, context, a, b)
    }

    fn recursive_compare_with_sort_rules(
        &self,
        sort_rules: &[&dyn CipSequenceRule],
        digraph: &mut CipDigraph<'_>,
        context: &mut CipLabelerContext,
        a: CipEdgeId,
        b: CipEdgeId,
    ) -> Result<i32, CipLabelerError> {
        recursive_compare_sequence_rule(self, Some(sort_rules), digraph, context, a, b)
    }

    // BEGIN RDKIT CPP FUNCTION SequenceRule::sort (rules/SequenceRule.cpp)
    // RDKit✔️✔️: Priority SequenceRule::sort(const Node *node, std::vector<Edge *> &edges,
    // RDKit✔️✔️:                             bool deep) const {
    // RDKit✔️✔️:   return getSorter()->prioritize(node, edges, deep);
    // RDKit✔️✔️: }
    // RDKit✔️✔️:
    // RDKit✔️✔️: Priority SequenceRule::sort(const Node *node,
    // RDKit✔️✔️:                             std::vector<Edge *> &edges) const {
    // RDKit✔️✔️:   return sort(node, edges, true);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION SequenceRule::sort
    fn sort(
        &self,
        digraph: &mut CipDigraph<'_>,
        context: &mut CipLabelerContext,
        node: CipNodeId,
        edges: &mut [CipEdgeId],
        deep: bool,
    ) -> Result<CipPriority, CipLabelerError> {
        prioritize_single_sequence_rule(self, digraph, context, node, edges, deep)
    }

    // BEGIN RDKIT CPP FUNCTION SequenceRule::areUpEdges (rules/SequenceRule.cpp)
    // RDKit✔️✔️: bool SequenceRule::areUpEdges(Node *aNode, Node *bNode, Edge *aEdge,
    // RDKit✔️✔️:                               Edge *bEdge) const {
    // RDKit✔️✔️:   // step over 'up' edges
    // RDKit✔️✔️:   if (aEdge->isEnd(aNode)) {
    // RDKit✔️✔️:     // if b is 'down' something's not right!
    // RDKit✔️✔️:     if (!bEdge->isEnd(bNode)) {
    // RDKit✔️✔️:       throw std::runtime_error("Something unexpected!");
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION SequenceRule::areUpEdges
    fn are_up_edges(
        &self,
        digraph: &CipDigraph<'_>,
        a_node: CipNodeId,
        b_node: CipNodeId,
        a_edge: CipEdgeId,
        b_edge: CipEdgeId,
    ) -> Result<bool, CipLabelerError> {
        if digraph.edge(a_edge).is_end(a_node) {
            if !digraph.edge(b_edge).is_end(b_node) {
                return Err(CipLabelerError::UnexpectedUpEdgeOrdering);
            }
            return Ok(true);
        }
        Ok(false)
    }
}

fn recursive_compare_sequence_rule<R: CipSequenceRule + ?Sized>(
    rule: &R,
    sort_rules: Option<&[&dyn CipSequenceRule]>,
    digraph: &mut CipDigraph<'_>,
    context: &mut CipLabelerContext,
    mut a: CipEdgeId,
    mut b: CipEdgeId,
) -> Result<i32, CipLabelerError> {
    if !context.decrement_remaining_call_count_and_check() {
        return Err(CipLabelerError::MaxIterationsExceeded);
    }

    let mut cmp = rule.compare_with_sort_rules(sort_rules, digraph, context, a, b)?;
    if cmp != 0 {
        return Ok(cmp);
    }

    let mut a_queue = vec![a];
    let mut b_queue = vec![b];
    let mut pos = 0_usize;
    while pos < a_queue.len() && pos < b_queue.len() {
        a = a_queue[pos];
        b = b_queue[pos];
        let a_end = digraph.edge(a).get_end();
        let b_end = digraph.edge(b).get_end();
        let mut as_edges = digraph.node_edges(a_end)?;
        let mut bs_edges = digraph.node_edges(b_end)?;

        sort_edges_for_sequence_rule(
            rule,
            sort_rules,
            digraph,
            context,
            a_end,
            &mut as_edges,
            false,
        )?;
        sort_edges_for_sequence_rule(
            rule,
            sort_rules,
            digraph,
            context,
            b_end,
            &mut bs_edges,
            false,
        )?;

        let sizediff = three_way_comparison_i32(as_edges.len() as i32, bs_edges.len() as i32);

        for (a_edge, b_edge) in as_edges.iter().zip(bs_edges.iter()) {
            if rule.are_up_edges(digraph, a_end, b_end, *a_edge, *b_edge)? {
                continue;
            }

            cmp = rule.compare_with_sort_rules(sort_rules, digraph, context, *a_edge, *b_edge)?;
            if cmp != 0 {
                return Ok(cmp);
            }
        }

        if sizediff != 0 {
            return Ok(sizediff);
        }

        sort_edges_for_sequence_rule(
            rule,
            sort_rules,
            digraph,
            context,
            a_end,
            &mut as_edges,
            true,
        )?;
        sort_edges_for_sequence_rule(
            rule,
            sort_rules,
            digraph,
            context,
            b_end,
            &mut bs_edges,
            true,
        )?;

        for (a_edge, b_edge) in as_edges.iter().zip(bs_edges.iter()) {
            if rule.are_up_edges(digraph, a_end, b_end, *a_edge, *b_edge)? {
                continue;
            }

            cmp = rule.compare_with_sort_rules(sort_rules, digraph, context, *a_edge, *b_edge)?;
            if cmp != 0 {
                return Ok(cmp);
            }

            a_queue.push(*a_edge);
            b_queue.push(*b_edge);
        }
        pos += 1;
    }
    Ok(0)
}

fn sort_edges_for_sequence_rule<R: CipSequenceRule + ?Sized>(
    rule: &R,
    sort_rules: Option<&[&dyn CipSequenceRule]>,
    digraph: &mut CipDigraph<'_>,
    context: &mut CipLabelerContext,
    node: CipNodeId,
    edges: &mut [CipEdgeId],
    deep: bool,
) -> Result<CipPriority, CipLabelerError> {
    if let Some(sort_rules) = sort_rules {
        CipSort::from_rules(sort_rules.to_vec()).prioritize(digraph, context, node, edges, deep)
    } else {
        rule.sort(digraph, context, node, edges, deep)
    }
}

fn prioritize_single_sequence_rule<R: CipSequenceRule + ?Sized>(
    rule: &R,
    digraph: &mut CipDigraph<'_>,
    context: &mut CipLabelerContext,
    node: CipNodeId,
    edges: &mut [CipEdgeId],
    deep: bool,
) -> Result<CipPriority, CipLabelerError> {
    let mut unique = true;
    let mut num_pseudo_asym = 0_i32;

    for i in 0..edges.len() {
        let mut j = i;
        while j > 0 {
            let cmp = compare_substituents_with_rule(
                rule,
                digraph,
                context,
                node,
                edges[j - 1],
                edges[j],
                deep,
            )?;

            if !(-1..=1).contains(&cmp) {
                num_pseudo_asym += 1;
            }

            if cmp < 0 {
                edges.swap(j, j - 1);
            } else {
                if cmp == 0 {
                    unique = false;
                }
                break;
            }
            j -= 1;
        }
    }

    Ok(CipPriority::new(unique, num_pseudo_asym == 1))
}

fn compare_substituents_with_rule<R: CipSequenceRule + ?Sized>(
    rule: &R,
    digraph: &mut CipDigraph<'_>,
    context: &mut CipLabelerContext,
    node: CipNodeId,
    a: CipEdgeId,
    b: CipEdgeId,
    deep: bool,
) -> Result<i32, CipLabelerError> {
    let a_is_beg = digraph.edge(a).is_beg(node);
    let b_is_beg = digraph.edge(b).is_beg(node);
    if !a_is_beg && b_is_beg {
        return Ok(1);
    } else if a_is_beg && !b_is_beg {
        return Ok(-1);
    }

    rule.get_comparison(digraph, context, a, b, deep)
}

fn three_way_comparison_i32(x: i32, y: i32) -> i32 {
    if x < y {
        -1
    } else if x == y {
        0
    } else {
        1
    }
}

pub(crate) struct CipSort<'r> {
    rules: Vec<&'r dyn CipSequenceRule>,
}

impl<'r> CipSort<'r> {
    // BEGIN RDKIT CPP FUNCTION Sort::Sort (CIPLabeler/Sort.cpp)
    // RDKit✔️✔️: Sort::Sort(const SequenceRule *comparator) : d_rules{comparator} {}
    // RDKit✔️✔️:
    // RDKit✔️✔️: Sort::Sort(std::vector<const SequenceRule *> comparators)
    // RDKit✔️✔️:     : d_rules{std::move(comparators)} {}
    // END RDKIT CPP FUNCTION Sort::Sort
    pub(crate) fn new(rule: &'r dyn CipSequenceRule) -> Self {
        Self { rules: vec![rule] }
    }

    pub(crate) fn from_rules(rules: Vec<&'r dyn CipSequenceRule>) -> Self {
        Self { rules }
    }

    // BEGIN RDKIT CPP FUNCTION Sort::getRules (CIPLabeler/Sort.cpp)
    // RDKit✔️✔️: const std::vector<const SequenceRule *> &Sort::getRules() const {
    // RDKit✔️✔️:   return d_rules;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Sort::getRules
    pub(crate) fn get_rules(&self) -> &[&'r dyn CipSequenceRule] {
        &self.rules
    }

    // BEGIN RDKIT CPP FUNCTION Sort::prioritize (CIPLabeler/Sort.cpp)
    // RDKit✔️✔️: Priority Sort::prioritize(const Node *node, std::vector<Edge *> &edges,
    // RDKit✔️✔️:                           bool deep) const {
    // RDKit✔️✔️:   bool unique = true;
    // RDKit✔️✔️:   int numPseudoAsym = 0;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (auto i = 0u; i < edges.size(); ++i) {
    // RDKit✔️✔️:     for (auto j = i; j > 0; --j) {
    // RDKit✔️✔️:       int cmp = compareSubstituents(node, edges[j - 1], edges[j], deep);
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (cmp < -1 || cmp > +1) {
    // RDKit✔️✔️:         ++numPseudoAsym;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (cmp < 0) {
    // RDKit✔️✔️:         std::swap(edges[j], edges[j - 1]);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         if (cmp == 0) {
    // RDKit✔️✔️:           unique = false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return {unique, numPseudoAsym == 1};
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Sort::prioritize
    pub(crate) fn prioritize(
        &self,
        digraph: &mut CipDigraph<'_>,
        context: &mut CipLabelerContext,
        node: CipNodeId,
        edges: &mut [CipEdgeId],
        deep: bool,
    ) -> Result<CipPriority, CipLabelerError> {
        let mut unique = true;
        let mut num_pseudo_asym = 0_i32;

        for i in 0..edges.len() {
            let mut j = i;
            while j > 0 {
                let cmp = self.compare_substituents(
                    digraph,
                    context,
                    node,
                    edges[j - 1],
                    edges[j],
                    deep,
                )?;

                if !(-1..=1).contains(&cmp) {
                    num_pseudo_asym += 1;
                }

                if cmp < 0 {
                    edges.swap(j, j - 1);
                } else {
                    if cmp == 0 {
                        unique = false;
                    }
                    break;
                }
                j -= 1;
            }
        }

        Ok(CipPriority::new(unique, num_pseudo_asym == 1))
    }

    // BEGIN RDKIT CPP FUNCTION Sort::getGroups (CIPLabeler/Sort.cpp)
    // RDKit✔️✔️: std::vector<std::vector<Edge *>> Sort::getGroups(
    // RDKit✔️✔️:     const std::vector<Edge *> &sorted) const {
    // RDKit✔️✔️:   // would be nice to have this integrated whilst sorting - may provide a
    // RDKit✔️✔️:   // small speed increase but as most of our lists are small we take use
    // RDKit✔️✔️:   // ugly sort then group approach
    // RDKit✔️✔️:   std::vector<std::vector<Edge *>> groups;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   Edge *prev = nullptr;
    // RDKit✔️✔️:   for (auto *edge : sorted) {
    // RDKit✔️✔️:     if (prev == nullptr ||
    // RDKit✔️✔️:         compareSubstituents(prev->getBeg(), prev, edge, true) != 0) {
    // RDKit✔️✔️:       groups.emplace_back();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     prev = edge;
    // RDKit✔️✔️:     groups.back().push_back(edge);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return groups;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Sort::getGroups
    pub(crate) fn get_groups(
        &self,
        digraph: &mut CipDigraph<'_>,
        context: &mut CipLabelerContext,
        sorted: &[CipEdgeId],
    ) -> Result<Vec<Vec<CipEdgeId>>, CipLabelerError> {
        let mut groups = Vec::<Vec<CipEdgeId>>::new();
        let mut prev = None;

        for edge in sorted {
            if prev.is_none()
                || self.compare_substituents(
                    digraph,
                    context,
                    digraph.edge(prev.expect("checked")).get_beg(),
                    prev.expect("checked"),
                    *edge,
                    true,
                )? != 0
            {
                groups.push(Vec::new());
            }
            prev = Some(*edge);
            groups
                .last_mut()
                .expect("RDKit Sort::getGroups creates a group before push")
                .push(*edge);
        }

        Ok(groups)
    }

    // BEGIN RDKIT CPP FUNCTION Sort::compareSubstituents (CIPLabeler/Sort.cpp)
    // RDKit✔️✔️: int Sort::compareSubstituents(const Node *node, const Edge *a, const Edge *b,
    // RDKit✔️✔️:                               bool deep) const {
    // RDKit✔️✔️:   // ensure 'out' edges are moved to the front
    // RDKit✔️✔️:   if (!a->isBeg(node) && b->isBeg(node)) {
    // RDKit✔️✔️:     return +1;
    // RDKit✔️✔️:   } else if (a->isBeg(node) && !b->isBeg(node)) {
    // RDKit✔️✔️:     return -1;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (const auto &rule : d_rules) {
    // RDKit✔️✔️:     int cmp = rule->getComparision(a, b, deep);
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (cmp != 0) {
    // RDKit✔️✔️:       return cmp;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return 0;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Sort::compareSubstituents
    fn compare_substituents(
        &self,
        digraph: &mut CipDigraph<'_>,
        context: &mut CipLabelerContext,
        node: CipNodeId,
        a: CipEdgeId,
        b: CipEdgeId,
        deep: bool,
    ) -> Result<i32, CipLabelerError> {
        let a_is_beg = digraph.edge(a).is_beg(node);
        let b_is_beg = digraph.edge(b).is_beg(node);
        if !a_is_beg && b_is_beg {
            return Ok(1);
        } else if a_is_beg && !b_is_beg {
            return Ok(-1);
        }

        for (idx, rule) in self.rules.iter().enumerate() {
            // RDKit assigns each SequenceRule the sorter snapshot that existed
            // when Rules::add() installed it.  Keep that prefix when the
            // aggregate sorter dispatches into the rule; passing the final
            // sorter here lets Rule4b/Rule5New see rules that did not exist in
            // their source snapshot and can cause unbounded recursive work.
            let rule_snapshot = &self.rules[..=idx];
            let cmp = rule.get_comparison_with_sort_rules(
                Some(rule_snapshot),
                digraph,
                context,
                a,
                b,
                deep,
            )?;
            if cmp != 0 {
                return Ok(cmp);
            }
        }
        Ok(0)
    }
}

pub(crate) struct CipRules {
    rules: Vec<Box<dyn CipSequenceRule>>,
}

impl CipRules {
    // BEGIN RDKIT CPP CLASS Rules constructor/add/destructor/getNumSubRules/getSorter (rules/Rules.h)
    // RDKit✔️❌: Rules() = delete;
    // RDKit✔️❌:
    // RDKit✔️❌: Rules(std::initializer_list<SequenceRule *> rules) {
    // RDKit✔️❌:   for (auto &rule : rules) {
    // RDKit✔️❌:     add(rule);
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // RDKit✔️❌:
    // RDKit✔️❌: ~Rules() override {
    // RDKit✔️❌:   for (auto &rule : d_rules) {
    // RDKit✔️❌:     delete rule;
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // RDKit✔️❌:
    // RDKit✔️❌: void add(SequenceRule *rule) {
    // RDKit✔️❌:   if (rule == nullptr) {
    // RDKit✔️❌:     throw std::runtime_error("No sequence rule provided");
    // RDKit✔️❌:   }
    // RDKit✔️❌:   d_rules.push_back(rule);
    // RDKit✔️❌:   rule->setSorter(new Sort(d_rules));
    // RDKit✔️❌: }
    // RDKit✔️❌:
    // RDKit✔️❌: int getNumSubRules() const { return d_rules.size(); }
    // RDKit✔️❌:
    // RDKit✔️❌: const Sort *getSorter() const override {
    // RDKit✔️❌:   if (dp_sorter == nullptr) {
    // RDKit✔️❌:     const_cast<Rules *>(this)->setSorter(new Sort(this));
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return dp_sorter.get();
    // RDKit✔️❌: }
    // END RDKIT CPP CLASS Rules constructor/add/destructor/getNumSubRules/getSorter
    // Performance note: unlike RDKit's sorter snapshot stored on each rule,
    // the detached Rust representation rebuilds a short borrowed-rule vector
    // when a comparison starts. The ordering is identical, but the extra
    // allocation is a known hot-path cost.
    pub(crate) fn new(rules: Vec<Box<dyn CipSequenceRule>>) -> Result<Self, CipLabelerError> {
        let mut result = Self { rules: Vec::new() };
        for rule in rules {
            result.add(Some(rule))?;
        }
        Ok(result)
    }

    fn add(&mut self, rule: Option<Box<dyn CipSequenceRule>>) -> Result<(), CipLabelerError> {
        let rule = rule.ok_or(CipLabelerError::NoSequenceRuleProvided)?;
        self.rules.push(rule);
        Ok(())
    }

    pub(crate) fn get_num_sub_rules(&self) -> usize {
        self.rules.len()
    }

    pub(crate) fn get_sorter(&self) -> CipSort<'_> {
        CipSort::new(self)
    }

    fn rule_refs(&self) -> Vec<&dyn CipSequenceRule> {
        self.rules
            .iter()
            .map(|rule| rule.as_ref() as &dyn CipSequenceRule)
            .collect()
    }
}

pub(crate) fn cip_all_rules() -> Result<CipRules, CipLabelerError> {
    // BEGIN RDKIT CPP CONSTANT all_rules (CIPLabeler.cpp)
    // RDKit✔️✔️: const Rules all_rules({new Rule1a, new Rule1b, new Rule2, new Rule3, new Rule4a,
    // RDKit✔️✔️:                        new Rule4b, new Rule4c, new Rule5New, new Rule6});
    // END RDKIT CPP CONSTANT all_rules
    // The pinned modern dispatcher instantiates Rule5New only. The separate
    // legacy Rule5 source is therefore intentionally unreachable here.
    CipRules::new(vec![
        Box::new(CipRule1a),
        Box::new(CipRule1b),
        Box::new(CipRule2),
        Box::new(CipRule3),
        Box::new(CipRule4a),
        Box::new(CipRule4b::new()),
        Box::new(CipRule4c),
        Box::new(CipRule5New::new()),
        Box::new(CipRule6),
    ])
}

pub(crate) fn cip_constitutional_rules() -> Result<CipRules, CipLabelerError> {
    // BEGIN RDKIT CPP CONSTANT constitutional_rules (CIPLabeler.cpp)
    // RDKit✔️✔️: const Rules constitutional_rules({new Rule1a, new Rule1b, new Rule2});
    // END RDKIT CPP CONSTANT constitutional_rules
    CipRules::new(vec![
        Box::new(CipRule1a),
        Box::new(CipRule1b),
        Box::new(CipRule2),
    ])
}

impl CipSequenceRule for CipRules {
    // BEGIN RDKIT CPP FUNCTION Rules::compare (rules/Rules.h)
    // RDKit✔️❌: int compare(const Edge *o1, const Edge *o2) const override {
    // RDKit✔️❌:   // Try using each rules. The rules will expand the search exhaustively
    // RDKit✔️❌:   // to all child substituents
    // RDKit✔️❌:   for (const auto &rule : d_rules) {
    // RDKit✔️❌:     // compare expands exhaustively across the whole graph
    // RDKit✔️❌:     int value = rule->recursiveCompare(o1, o2);
    // RDKit✔️❌:     if (value != 0) {
    // RDKit✔️❌:       return value;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return 0;
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION Rules::compare
    // `rule_refs()` materializes the source pointer view for Rust's borrow
    // checker; this preserves behavior but adds an allocation per dispatch.
    fn compare(
        &self,
        digraph: &mut CipDigraph<'_>,
        context: &mut CipLabelerContext,
        a: CipEdgeId,
        b: CipEdgeId,
    ) -> Result<i32, CipLabelerError> {
        let sort_rules = self.rule_refs();
        for (idx, rule) in sort_rules.iter().enumerate() {
            let rule_snapshot = &sort_rules[..=idx];
            let value =
                rule.recursive_compare_with_sort_rules(rule_snapshot, digraph, context, a, b)?;
            if value != 0 {
                return Ok(value);
            }
        }
        Ok(0)
    }

    // BEGIN RDKIT CPP FUNCTION Rules::getComparision (rules/Rules.h)
    // RDKit✔️❌: int getComparision(const Edge *a, const Edge *b, bool deep) const override {
    // RDKit✔️❌:   (void)deep;
    // RDKit✔️❌:
    // RDKit✔️❌:   // Try using each rules. The rules will expand the search exhaustively
    // RDKit✔️❌:   // to all child substituents
    // RDKit✔️❌:   for (const auto &rule : d_rules) {
    // RDKit✔️❌:     // compare expands exhaustively across the whole graph
    // RDKit✔️❌:     int value = rule->recursiveCompare(a, b);
    // RDKit✔️❌:
    // RDKit✔️❌:     if (value != 0) {
    // RDKit✔️❌:       return value;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   return 0;
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION Rules::getComparision
    fn get_comparison(
        &self,
        digraph: &mut CipDigraph<'_>,
        context: &mut CipLabelerContext,
        a: CipEdgeId,
        b: CipEdgeId,
        _deep: bool,
    ) -> Result<i32, CipLabelerError> {
        let sort_rules = self.rule_refs();
        for (idx, rule) in sort_rules.iter().enumerate() {
            let rule_snapshot = &sort_rules[..=idx];
            let value =
                rule.recursive_compare_with_sort_rules(rule_snapshot, digraph, context, a, b)?;
            if value != 0 {
                return Ok(value);
            }
        }
        Ok(0)
    }

    // BEGIN RDKIT CPP FUNCTION Rules::getComparision (rules/Rules.h)
    // RDKit✔️❌: int getComparision(const Edge *a, const Edge *b, bool deep) const override {
    // RDKit✔️❌:   (void)deep;
    // RDKit✔️❌:
    // RDKit✔️❌:   // Try using each rules. The rules will expand the search exhaustively
    // RDKit✔️❌:   // to all child substituents
    // RDKit✔️❌:   for (const auto &rule : d_rules) {
    // RDKit✔️❌:     // compare expands exhaustively across the whole graph
    // RDKit✔️❌:     int value = rule->recursiveCompare(a, b);
    // RDKit✔️❌:
    // RDKit✔️❌:     if (value != 0) {
    // RDKit✔️❌:       return value;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   return 0;
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION Rules::getComparision
    fn get_comparison_with_sort_rules(
        &self,
        _sort_rules: Option<&[&dyn CipSequenceRule]>,
        digraph: &mut CipDigraph<'_>,
        context: &mut CipLabelerContext,
        a: CipEdgeId,
        b: CipEdgeId,
        _deep: bool,
    ) -> Result<i32, CipLabelerError> {
        let sort_rules = self.rule_refs();
        for (idx, rule) in sort_rules.iter().enumerate() {
            let rule_snapshot = &sort_rules[..=idx];
            let value =
                rule.recursive_compare_with_sort_rules(rule_snapshot, digraph, context, a, b)?;
            if value != 0 {
                return Ok(value);
            }
        }
        Ok(0)
    }

    fn sort(
        &self,
        digraph: &mut CipDigraph<'_>,
        context: &mut CipLabelerContext,
        node: CipNodeId,
        edges: &mut [CipEdgeId],
        deep: bool,
    ) -> Result<CipPriority, CipLabelerError> {
        self.get_sorter()
            .prioritize(digraph, context, node, edges, deep)
    }
}

#[derive(Debug, Default, Clone, Copy)]
pub(crate) struct CipRule1a;

impl CipSequenceRule for CipRule1a {
    // BEGIN RDKIT CPP FUNCTION Rule1a::compare (rules/Rule1a.cpp)
    // RDKit✔️✔️: Rule1a::Rule1a() = default;
    // RDKit✔️✔️:
    // RDKit✔️✔️: // CIP Rule 1a: Higher atomic number precedes lower.
    // RDKit✔️✔️: int Rule1a::compare(const Edge *a, const Edge *b) const {
    // RDKit✔️✔️:   const auto afrac = a->getEnd()->getAtomicNumFraction();
    // RDKit✔️✔️:   const auto bfrac = b->getEnd()->getAtomicNumFraction();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return three_way_comparison(afrac, bfrac);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Rule1a::compare
    fn compare(
        &self,
        digraph: &mut CipDigraph<'_>,
        _context: &mut CipLabelerContext,
        a: CipEdgeId,
        b: CipEdgeId,
    ) -> Result<i32, CipLabelerError> {
        let afrac = digraph
            .node(digraph.edge(a).get_end())
            .get_atomic_num_fraction();
        let bfrac = digraph
            .node(digraph.edge(b).get_end())
            .get_atomic_num_fraction();
        Ok(match afrac.cmp(&bfrac) {
            std::cmp::Ordering::Less => -1,
            std::cmp::Ordering::Equal => 0,
            std::cmp::Ordering::Greater => 1,
        })
    }
}

#[derive(Debug, Default, Clone, Copy)]
pub(crate) struct CipRule1b;

impl CipRule1b {
    // BEGIN RDKIT CPP CONSTANT Rule1b::IUPAC_2013 (rules/Rule1b.h)
    // RDKit✔️✔️: static const bool IUPAC_2013 = false;
    // END RDKIT CPP CONSTANT Rule1b::IUPAC_2013
    const IUPAC_2013: bool = false;
}

impl CipSequenceRule for CipRule1b {
    // BEGIN RDKIT CPP FUNCTION Rule1b::compare (rules/Rule1b.cpp)
    // RDKit✔️✔️: Rule1b::Rule1b() = default;
    // RDKit✔️✔️:
    // RDKit✔️✔️: int Rule1b::compare(const Edge *a, const Edge *b) const {
    // RDKit✔️✔️:   if (IUPAC_2013) {
    // RDKit✔️✔️:     return -three_way_comparison(a->getEnd()->getDistance(),
    // RDKit✔️✔️:                                  b->getEnd()->getDistance());
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     if (a->getEnd()->isSet(Node::RING_DUPLICATE) &&
    // RDKit✔️✔️:         b->getEnd()->isSet(Node::RING_DUPLICATE)) {
    // RDKit✔️✔️:       return -three_way_comparison(a->getEnd()->getDistance(),
    // RDKit✔️✔️:                                    b->getEnd()->getDistance());
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       if (a->getEnd()->isSet(Node::RING_DUPLICATE) &&
    // RDKit✔️✔️:           !b->getEnd()->isSet(Node::RING_DUPLICATE)) {
    // RDKit✔️✔️:         return +1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (!a->getEnd()->isSet(Node::RING_DUPLICATE) &&
    // RDKit✔️✔️:           b->getEnd()->isSet(Node::RING_DUPLICATE)) {
    // RDKit✔️✔️:         return -1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Rule1b::compare
    fn compare(
        &self,
        digraph: &mut CipDigraph<'_>,
        _context: &mut CipLabelerContext,
        a: CipEdgeId,
        b: CipEdgeId,
    ) -> Result<i32, CipLabelerError> {
        let a_end = digraph.edge(a).get_end();
        let b_end = digraph.edge(b).get_end();
        let a_node = digraph.node(a_end);
        let b_node = digraph.node(b_end);
        if Self::IUPAC_2013 {
            return Ok(-three_way_comparison_i32(
                a_node.get_distance(),
                b_node.get_distance(),
            ));
        }
        if a_node.is_set(CipNode::RING_DUPLICATE) && b_node.is_set(CipNode::RING_DUPLICATE) {
            return Ok(-three_way_comparison_i32(
                a_node.get_distance(),
                b_node.get_distance(),
            ));
        }
        if a_node.is_set(CipNode::RING_DUPLICATE) && !b_node.is_set(CipNode::RING_DUPLICATE) {
            return Ok(1);
        }
        if !a_node.is_set(CipNode::RING_DUPLICATE) && b_node.is_set(CipNode::RING_DUPLICATE) {
            return Ok(-1);
        }
        Ok(0)
    }
}

#[derive(Debug, Default, Clone, Copy)]
pub(crate) struct CipRule2;

impl CipSequenceRule for CipRule2 {
    // BEGIN RDKIT CPP FUNCTION Rule2::compare (rules/Rule2.cpp)
    // RDKit✔️✔️: Rule2::Rule2() = default;
    // RDKit✔️✔️:
    // RDKit✔️✔️: int Rule2::compare(const Edge *a, const Edge *b) const {
    // RDKit✔️✔️:   auto a_end = a->getEnd();
    // RDKit✔️✔️:   auto b_end = b->getEnd();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto aAtomNum = a_end->getAtomicNum();
    // RDKit✔️✔️:   auto bAtomNum = b_end->getAtomicNum();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (aAtomNum == 0 && bAtomNum == 0) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   } else if (aAtomNum == 0 || bAtomNum == 0) {
    // RDKit✔️✔️:     // This should be caught by Rule 1a, but just in case
    // RDKit✔️✔️:     return three_way_comparison(aAtomNum, bAtomNum);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto aMassNum = a_end->getMassNum();
    // RDKit✔️✔️:   auto bMassNum = b_end->getMassNum();
    // RDKit✔️✔️:   if (aMassNum == 0u && bMassNum == 0u) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto aweight = a_end->getAtomicMass();
    // RDKit✔️✔️:   auto bweight = b_end->getAtomicMass();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return three_way_comparison(aweight, bweight);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Rule2::compare
    fn compare(
        &self,
        digraph: &mut CipDigraph<'_>,
        _context: &mut CipLabelerContext,
        a: CipEdgeId,
        b: CipEdgeId,
    ) -> Result<i32, CipLabelerError> {
        let a_end = digraph.edge(a).get_end();
        let b_end = digraph.edge(b).get_end();
        let a_node = digraph.node(a_end);
        let b_node = digraph.node(b_end);

        let a_atom_num = a_node.get_atomic_num(digraph.mol())?;
        let b_atom_num = b_node.get_atomic_num(digraph.mol())?;
        if a_atom_num == 0 && b_atom_num == 0 {
            return Ok(0);
        } else if a_atom_num == 0 || b_atom_num == 0 {
            return Ok(three_way_comparison_i32(
                i32::from(a_atom_num),
                i32::from(b_atom_num),
            ));
        }

        let a_mass_num = a_node.get_mass_num(digraph.mol())?;
        let b_mass_num = b_node.get_mass_num(digraph.mol())?;
        if a_mass_num == 0 && b_mass_num == 0 {
            return Ok(0);
        }

        Ok(
            match a_node
                .get_atomic_mass()
                .partial_cmp(&b_node.get_atomic_mass())
                .expect("RDKit CIP atomic masses are finite")
            {
                std::cmp::Ordering::Less => -1,
                std::cmp::Ordering::Equal => 0,
                std::cmp::Ordering::Greater => 1,
            },
        )
    }
}

#[derive(Debug, Default, Clone, Copy)]
pub(crate) struct CipRule3;

impl CipRule3 {
    // BEGIN RDKIT CPP FUNCTION ord (rules/Rule3.cpp)
    // RDKit✔️✔️: int ord(Descriptor lab) {
    // RDKit✔️✔️:   switch (lab) {
    // RDKit✔️✔️:     case Descriptor::E:
    // RDKit✔️✔️:       return 1;
    // RDKit✔️✔️:     case Descriptor::Z:
    // RDKit✔️✔️:       return 2;
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION ord
    fn ord(lab: Descriptor) -> i32 {
        match lab {
            Descriptor::E => 1,
            Descriptor::Z => 2,
            _ => 0,
        }
    }
}

impl CipSequenceRule for CipRule3 {
    // BEGIN RDKIT CPP FUNCTION Rule3::compare (rules/Rule3.cpp)
    // RDKit✔️✔️: Rule3::Rule3() = default;
    // RDKit✔️✔️:
    // RDKit✔️✔️: int Rule3::compare(const Edge *a, const Edge *b) const {
    // RDKit✔️✔️:   return three_way_comparison(ord(a->getEnd()->getAux()),
    // RDKit✔️✔️:                               ord(b->getEnd()->getAux()));
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Rule3::compare
    fn compare(
        &self,
        digraph: &mut CipDigraph<'_>,
        _context: &mut CipLabelerContext,
        a: CipEdgeId,
        b: CipEdgeId,
    ) -> Result<i32, CipLabelerError> {
        let a_ord = Self::ord(digraph.node(digraph.edge(a).get_end()).get_aux());
        let b_ord = Self::ord(digraph.node(digraph.edge(b).get_end()).get_aux());
        Ok(three_way_comparison_i32(a_ord, b_ord))
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct CipPairList {
    descriptors: Vec<Descriptor>,
    pairing: CipPairing,
}

impl Default for CipPairList {
    fn default() -> Self {
        Self::new()
    }
}

impl CipPairList {
    const NUM_PAIRING_BITS: usize = 64;

    // BEGIN RDKIT CPP CLASS PairList (rules/Pairlist.h)
    // RDKit✔️✔️: using pairing_t = std::uint64_t;
    // RDKit✔️✔️: static constexpr int numPairingBits = sizeof(pairing_t) * 8;
    // RDKit✔️✔️:
    // RDKit✔️✔️: static Descriptor ref(Descriptor descriptor) {
    // RDKit✔️✔️:   switch (descriptor) {
    // RDKit✔️✔️:     case Descriptor::R:
    // RDKit✔️✔️:     case Descriptor::M:
    // RDKit✔️✔️:     case Descriptor::seqCis:
    // RDKit✔️✔️:       return Descriptor::R;
    // RDKit✔️✔️:     case Descriptor::S:
    // RDKit✔️✔️:     case Descriptor::P:
    // RDKit✔️✔️:     case Descriptor::seqTrans:
    // RDKit✔️✔️:       return Descriptor::S;
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       return Descriptor::NONE;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP CLASS PairList ref
    pub(crate) fn ref_descriptor(descriptor: Descriptor) -> Descriptor {
        match descriptor {
            Descriptor::R | Descriptor::M | Descriptor::seqCis => Descriptor::R,
            Descriptor::S | Descriptor::P | Descriptor::seqTrans => Descriptor::S,
            _ => Descriptor::None,
        }
    }

    // BEGIN RDKIT CPP FUNCTION PairList constructors (rules/Pairlist.h)
    // RDKit✔️✔️: PairList() = default;
    // RDKit✔️✔️:
    // RDKit✔️✔️: PairList(Descriptor ref) { add(ref); }
    // END RDKIT CPP FUNCTION PairList constructors
    pub(crate) fn new() -> Self {
        Self {
            descriptors: Vec::new(),
            pairing: 0,
        }
    }

    pub(crate) fn with_ref(ref_descriptor: Descriptor) -> Self {
        let mut result = Self::new();
        result.add(ref_descriptor);
        result
    }

    // BEGIN RDKIT CPP FUNCTION PairList::PairList(head, tail) (rules/Pairlist.h)
    // RDKit✔️✔️: PairList(const PairList &head, const PairList &tail) {
    // RDKit✔️✔️:   // add descriptors to the new instance (ignored descriptors not added)
    // RDKit✔️✔️:   addAll(head.d_descriptors);
    // RDKit✔️✔️:   addAll(tail.d_descriptors);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION PairList::PairList(head, tail)
    pub(crate) fn from_head_tail(head: &Self, tail: &Self) -> Self {
        let mut result = Self::new();
        result.add_all(&head.descriptors);
        result.add_all(&tail.descriptors);
        result
    }

    // BEGIN RDKIT CPP FUNCTION PairList::getRefDescriptor (rules/Pairlist.h)
    // RDKit✔️✔️: Descriptor getRefDescriptor() const { return ref(d_descriptors[0]); }
    // END RDKIT CPP FUNCTION PairList::getRefDescriptor
    pub(crate) fn get_ref_descriptor(&self) -> Descriptor {
        Self::ref_descriptor(self.descriptors[0])
    }

    // BEGIN RDKIT CPP FUNCTION PairList::add (rules/Pairlist.h)
    // RDKit✔️✔️: bool add(Descriptor descriptor) {
    // RDKit✔️✔️:   switch (descriptor) {
    // RDKit✔️✔️:     case Descriptor::R:
    // RDKit✔️✔️:     case Descriptor::S:
    // RDKit✔️✔️:     case Descriptor::M:
    // RDKit✔️✔️:     case Descriptor::P:
    // RDKit✔️✔️:     case Descriptor::seqTrans:
    // RDKit✔️✔️:     case Descriptor::seqCis:
    // RDKit✔️✔️:       addAndPair(descriptor);
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION PairList::add
    pub(crate) fn add(&mut self, descriptor: Descriptor) -> bool {
        match descriptor {
            Descriptor::R
            | Descriptor::S
            | Descriptor::M
            | Descriptor::P
            | Descriptor::seqTrans
            | Descriptor::seqCis => {
                self.add_and_pair(descriptor);
                true
            }
            _ => false,
        }
    }

    // BEGIN RDKIT CPP FUNCTION PairList::addAll (rules/Pairlist.h)
    // RDKit✔️✔️: template <typename T>
    // RDKit✔️✔️: void addAll(const T &descriptors) {
    // RDKit✔️✔️:   for (const auto &descriptor : descriptors) {
    // RDKit✔️✔️:     add(descriptor);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION PairList::addAll
    pub(crate) fn add_all(&mut self, descriptors: &[Descriptor]) {
        for descriptor in descriptors {
            self.add(*descriptor);
        }
    }

    // BEGIN RDKIT CPP FUNCTION PairList::getPairing (rules/Pairlist.h)
    // RDKit✔️✔️: pairing_t getPairing() const { return d_pairing; }
    // END RDKIT CPP FUNCTION PairList::getPairing
    pub(crate) fn get_pairing(&self) -> CipPairing {
        self.pairing
    }

    // BEGIN RDKIT CPP FUNCTION PairList::compareTo/operator< (rules/Pairlist.h)
    // RDKit✔️✔️: int compareTo(const PairList &that) const {
    // RDKit✔️✔️:   if (d_descriptors.size() != that.d_descriptors.size()) {
    // RDKit✔️✔️:     throw std::runtime_error("Descriptor lists should be the same length!");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   Descriptor thisRef = d_descriptors[0];
    // RDKit✔️✔️:   Descriptor thatRef = that.d_descriptors[0];
    // RDKit✔️✔️:   for (auto i = 1u; i < d_descriptors.size(); ++i) {
    // RDKit✔️✔️:     if (thisRef == d_descriptors[i] && thatRef != that.d_descriptors[i]) {
    // RDKit✔️✔️:       return +1;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (thisRef != d_descriptors[i] && thatRef == that.d_descriptors[i]) {
    // RDKit✔️✔️:       return -1;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return 0;
    // RDKit✔️✔️: }
    // RDKit✔️✔️:
    // RDKit✔️✔️: bool operator<(const PairList &that) const { return compareTo(that) == -1; }
    // END RDKIT CPP FUNCTION PairList::compareTo/operator<
    pub(crate) fn compare_to(&self, that: &Self) -> Result<i32, CipLabelerError> {
        if self.descriptors.len() != that.descriptors.len() {
            return Err(CipLabelerError::DescriptorListLengthMismatch);
        }
        let this_ref = self.descriptors[0];
        let that_ref = that.descriptors[0];
        for i in 1..self.descriptors.len() {
            if this_ref == self.descriptors[i] && that_ref != that.descriptors[i] {
                return Ok(1);
            }
            if this_ref != self.descriptors[i] && that_ref == that.descriptors[i] {
                return Ok(-1);
            }
        }
        Ok(0)
    }

    pub(crate) fn less_than(&self, that: &Self) -> Result<bool, CipLabelerError> {
        Ok(self.compare_to(that)? == -1)
    }

    fn sort_descending(lists: &mut [Self]) -> Result<(), CipLabelerError> {
        let Some(expected_len) = lists.first().map(|list| list.descriptors.len()) else {
            return Ok(());
        };
        if lists
            .iter()
            .any(|list| list.descriptors.len() != expected_len)
        {
            return Err(CipLabelerError::DescriptorListLengthMismatch);
        }
        lists.sort_unstable_by(|left, right| {
            right
                .compare_to(left)
                .expect("descriptor lengths were validated before sorting")
                .cmp(&0)
        });
        Ok(())
    }

    // BEGIN RDKIT CPP FUNCTION PairList::toString (rules/Pairlist.h)
    // RDKit✔️✔️: std::string toString() const {
    // RDKit✔️✔️:   // handles cases that would break the toString method
    // RDKit✔️✔️:   if (d_descriptors.empty() || d_descriptors[0] == Descriptor::NONE) {
    // RDKit✔️✔️:     return "";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::stringstream ss;
    // RDKit✔️✔️:   auto basis = d_descriptors[0];
    // RDKit✔️✔️:   ss << to_string(basis) << ':';
    // RDKit✔️✔️:
    // RDKit✔️✔️:   basis = ref(basis);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // build like (l) / unlike (u) descriptor pairing
    // RDKit✔️✔️:   for (auto it = d_descriptors.begin() + 1; it != d_descriptors.end(); ++it) {
    // RDKit✔️✔️:     ss << (basis == ref(*it) ? "l" : "u");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return ss.str();
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION PairList::toString
    pub(crate) fn to_rdkit_string(&self) -> String {
        if self.descriptors.is_empty() || self.descriptors[0] == Descriptor::None {
            return String::new();
        }

        let mut result = String::new();
        let mut basis = self.descriptors[0];
        result.push_str(descriptor_to_string(basis));
        result.push(':');

        basis = Self::ref_descriptor(basis);
        for descriptor in self.descriptors.iter().skip(1) {
            result.push(if basis == Self::ref_descriptor(*descriptor) {
                'l'
            } else {
                'u'
            });
        }
        result
    }

    // BEGIN RDKIT CPP FUNCTION PairList::addAndPair (rules/Pairlist.h)
    // RDKit✔️✔️: void addAndPair(Descriptor descriptor) {
    // RDKit✔️✔️:   // if this isn't the first descriptor - check the pairing
    // RDKit✔️✔️:   if (!d_descriptors.empty() && d_descriptors[0] == descriptor) {
    // RDKit✔️✔️:     // set the bit to indicate a pair
    // RDKit✔️✔️:     d_pairing |= static_cast<pairing_t>(1)
    // RDKit✔️✔️:                  << (numPairingBits - 1 - d_descriptors.size());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   d_descriptors.push_back(ref(descriptor));
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION PairList::addAndPair
    fn add_and_pair(&mut self, descriptor: Descriptor) {
        if !self.descriptors.is_empty() && self.descriptors[0] == descriptor {
            self.pairing |= 1_u64 << (Self::NUM_PAIRING_BITS - 1 - self.descriptors.len());
        }
        self.descriptors.push(Self::ref_descriptor(descriptor));
    }
}

#[derive(Debug, Default, Clone, Copy)]
pub(crate) struct CipRule4a;

impl CipRule4a {
    // BEGIN RDKIT CPP FUNCTION ord (rules/Rule4a.cpp)
    // RDKit✔️✔️: int ord(Descriptor lab) {
    // RDKit✔️✔️:   switch (lab) {
    // RDKit✔️✔️:     case Descriptor::UNKNOWN:
    // RDKit✔️✔️:     case Descriptor::ns:
    // RDKit✔️✔️:     case Descriptor::NONE:
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     case Descriptor::r:
    // RDKit✔️✔️:     case Descriptor::s:
    // RDKit✔️✔️:     case Descriptor::m:
    // RDKit✔️✔️:     case Descriptor::p:
    // RDKit✔️✔️:     case Descriptor::E:
    // RDKit✔️✔️:     case Descriptor::Z:
    // RDKit✔️✔️:       return 1;
    // RDKit✔️✔️:     case Descriptor::R:
    // RDKit✔️✔️:     case Descriptor::S:
    // RDKit✔️✔️:     case Descriptor::M:
    // RDKit✔️✔️:     case Descriptor::P:
    // RDKit✔️✔️:     case Descriptor::seqTrans:
    // RDKit✔️✔️:     case Descriptor::seqCis:
    // RDKit✔️✔️:       return 2;
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       throw std::logic_error("Invalid stereo descriptor");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION ord
    fn ord(lab: Descriptor) -> Result<i32, CipLabelerError> {
        match lab {
            Descriptor::Unknown | Descriptor::ns | Descriptor::None => Ok(0),
            Descriptor::r
            | Descriptor::s
            | Descriptor::m
            | Descriptor::p
            | Descriptor::E
            | Descriptor::Z => Ok(1),
            Descriptor::R
            | Descriptor::S
            | Descriptor::M
            | Descriptor::P
            | Descriptor::seqTrans
            | Descriptor::seqCis => Ok(2),
            _ => Err(CipLabelerError::InvalidStereoDescriptor),
        }
    }
}

impl CipSequenceRule for CipRule4a {
    // BEGIN RDKIT CPP FUNCTION Rule4a::compare (rules/Rule4a.cpp)
    // RDKit✔️✔️: Rule4a::Rule4a() = default;
    // RDKit✔️✔️:
    // RDKit✔️✔️: int Rule4a::compare(const Edge *a, const Edge *b) const {
    // RDKit✔️✔️:   int aOrdinal = ord(getBondLabel(a));
    // RDKit✔️✔️:   int bOrdinal = ord(getBondLabel(b));
    // RDKit✔️✔️:   int cmp = three_way_comparison(aOrdinal, bOrdinal);
    // RDKit✔️✔️:   if (cmp != 0) {
    // RDKit✔️✔️:     return cmp;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   aOrdinal = ord(a->getEnd()->getAux());
    // RDKit✔️✔️:   bOrdinal = ord(b->getEnd()->getAux());
    // RDKit✔️✔️:   return three_way_comparison(aOrdinal, bOrdinal);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Rule4a::compare
    fn compare(
        &self,
        digraph: &mut CipDigraph<'_>,
        _context: &mut CipLabelerContext,
        a: CipEdgeId,
        b: CipEdgeId,
    ) -> Result<i32, CipLabelerError> {
        let a_ordinal = Self::ord(self.get_bond_label(digraph.edge(a)))?;
        let b_ordinal = Self::ord(self.get_bond_label(digraph.edge(b)))?;
        let cmp = three_way_comparison_i32(a_ordinal, b_ordinal);
        if cmp != 0 {
            return Ok(cmp);
        }
        let a_ordinal = Self::ord(digraph.node(digraph.edge(a).get_end()).get_aux())?;
        let b_ordinal = Self::ord(digraph.node(digraph.edge(b).get_end()).get_aux())?;
        Ok(three_way_comparison_i32(a_ordinal, b_ordinal))
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct CipRule4b {
    ref_descriptor: Descriptor,
}

impl CipRule4b {
    pub(crate) fn new() -> Self {
        Self {
            ref_descriptor: Descriptor::None,
        }
    }

    pub(crate) fn with_ref(ref_descriptor: Descriptor) -> Self {
        Self { ref_descriptor }
    }

    // BEGIN RDKIT CPP FUNCTION Rule4b::getReferenceDescriptors (rules/Rule4b.cpp)
    // RDKit✔️✔️: std::vector<Descriptor> Rule4b::getReferenceDescriptors(
    // RDKit✔️✔️:     const Node *node) const {
    // RDKit✔️✔️:   std::vector<Descriptor> result;
    // RDKit✔️✔️:   auto prev = initialLevel(node);
    // RDKit✔️✔️:   while (!prev.empty()) {
    // RDKit✔️✔️:     for (const auto &nodes : prev) {
    // RDKit✔️✔️:       if (getReference(nodes, result)) {
    // RDKit✔️✔️:         return result;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     prev = getNextLevel(prev);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return {};
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Rule4b::getReferenceDescriptors
    fn get_reference_descriptors(
        &self,
        sort_rules: Option<&[&dyn CipSequenceRule]>,
        digraph: &mut CipDigraph<'_>,
        context: &mut CipLabelerContext,
        node: CipNodeId,
    ) -> Result<Vec<Descriptor>, CipLabelerError> {
        let mut result = Vec::new();
        let source_sort_rules = if let Some(sort_rules) = sort_rules {
            let self_ptr = self as &dyn CipSequenceRule as *const dyn CipSequenceRule;
            let end = sort_rules
                .iter()
                .position(|rule| std::ptr::addr_eq(*rule as *const dyn CipSequenceRule, self_ptr))
                .map(|index| index + 1)
                .ok_or(CipLabelerError::Rule4bInstanceNotInRuleSet)?;
            Some(sort_rules[..end].to_vec())
        } else {
            None
        };
        let mut prev = self.initial_level(node);
        while !prev.is_empty() {
            for nodes in &prev {
                if self.get_reference(digraph, nodes, &mut result) {
                    return Ok(result);
                }
            }
            prev = self.get_next_level(source_sort_rules.as_deref(), digraph, context, &prev)?;
        }
        Ok(Vec::new())
    }

    // BEGIN RDKIT CPP FUNCTION Rule4b::hasDescriptors (rules/Rule4b.cpp)
    // RDKit✔️✔️: bool Rule4b::hasDescriptors(const Node *node) const {
    // RDKit✔️✔️:   auto queue = std::list<const Node *>({node});
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (const auto &node : queue) {
    // RDKit✔️✔️:     if (node->getAux() != Descriptor::NONE) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (const auto &e : node->getEdges()) {
    // RDKit✔️✔️:       if (e->getEnd() == node) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (getBondLabel(e) != Descriptor::NONE) {
    // RDKit✔️✔️:         return true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       queue.push_back(e->getEnd());
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Rule4b::hasDescriptors
    fn has_descriptors(
        &self,
        digraph: &mut CipDigraph<'_>,
        node: CipNodeId,
    ) -> Result<bool, CipLabelerError> {
        let mut queue = vec![node];
        let mut pos = 0_usize;
        while pos < queue.len() {
            let current = queue[pos];
            if digraph.node(current).get_aux() != Descriptor::None {
                return Ok(true);
            }
            for edge_id in digraph.node_edges(current)? {
                let edge = digraph.edge(edge_id);
                if edge.get_end() == current {
                    continue;
                }
                if self.get_bond_label(edge) != Descriptor::None {
                    return Ok(true);
                }
                queue.push(edge.get_end());
            }
            pos += 1;
        }
        Ok(false)
    }

    // BEGIN RDKIT CPP FUNCTION Rule4b::getReference (rules/Rule4b.cpp)
    // RDKit✔️✔️: bool Rule4b::getReference(const std::vector<const Node *> &nodes,
    // RDKit✔️✔️:                           std::vector<Descriptor> &result) const {
    // RDKit✔️✔️:   int right = 0;
    // RDKit✔️✔️:   int left = 0;
    // RDKit✔️✔️:   for (const auto &node : nodes) {
    // RDKit✔️✔️:     auto desc = node->getAux();
    // RDKit✔️✔️:     switch (desc) {
    // RDKit✔️✔️:       case Descriptor::NONE:
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       case Descriptor::R:
    // RDKit✔️✔️:       case Descriptor::M:
    // RDKit✔️✔️:       case Descriptor::seqCis:
    // RDKit✔️✔️:         ++right;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       case Descriptor::S:
    // RDKit✔️✔️:       case Descriptor::P:
    // RDKit✔️✔️:       case Descriptor::seqTrans:
    // RDKit✔️✔️:         ++left;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       default:
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (right + left == 0) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   } else if (right > left) {
    // RDKit✔️✔️:     result.push_back(Descriptor::R);
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   } else if (right < left) {
    // RDKit✔️✔️:     result.push_back(Descriptor::S);
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     result.push_back(Descriptor::R);
    // RDKit✔️✔️:     result.push_back(Descriptor::S);
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Rule4b::getReference
    fn get_reference(
        &self,
        digraph: &CipDigraph<'_>,
        nodes: &[CipNodeId],
        result: &mut Vec<Descriptor>,
    ) -> bool {
        let mut right = 0_i32;
        let mut left = 0_i32;
        for node in nodes {
            match digraph.node(*node).get_aux() {
                Descriptor::None => continue,
                Descriptor::R | Descriptor::M | Descriptor::seqCis => right += 1,
                Descriptor::S | Descriptor::P | Descriptor::seqTrans => left += 1,
                _ => {}
            }
        }
        if right + left == 0 {
            false
        } else if right > left {
            result.push(Descriptor::R);
            true
        } else if right < left {
            result.push(Descriptor::S);
            true
        } else {
            result.push(Descriptor::R);
            result.push(Descriptor::S);
            true
        }
    }

    // BEGIN RDKIT CPP FUNCTION Rule4b::initialLevel (rules/Rule4b.cpp)
    // RDKit✔️✔️: std::vector<std::vector<const Node *>> Rule4b::initialLevel(
    // RDKit✔️✔️:     const Node *node) const {
    // RDKit✔️✔️:   return {{node}};
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Rule4b::initialLevel
    fn initial_level(&self, node: CipNodeId) -> Vec<Vec<CipNodeId>> {
        vec![vec![node]]
    }

    // BEGIN RDKIT CPP FUNCTION Rule4b::getNextLevel (rules/Rule4b.cpp)
    // RDKit✔️✔️: std::vector<std::vector<const Node *>> Rule4b::getNextLevel(
    // RDKit✔️✔️:     const std::vector<std::vector<const Node *>> &prevLevel) const {
    // RDKit✔️✔️:   std::vector<std::vector<const Node *>> nextLevel;
    // RDKit✔️✔️:   nextLevel.reserve(4 * prevLevel.size());
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (const auto &prev : prevLevel) {
    // RDKit✔️✔️:     std::vector<std::vector<std::vector<Edge *>>> tmp;
    // RDKit✔️✔️:     for (const auto &node : prev) {
    // RDKit✔️✔️:       auto edges = node->getNonTerminalOutEdges();
    // RDKit✔️✔️:       sort(node, edges);
    // RDKit✔️✔️:       tmp.push_back(getSorter()->getGroups(edges));
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // check sizes
    // RDKit✔️✔️:     int size = -1;
    // RDKit✔️✔️:     for (auto i = 0u; i < tmp.size(); ++i) {
    // RDKit✔️✔️:       int localSize = tmp[0].size();
    // RDKit✔️✔️:       if (size < 0) {
    // RDKit✔️✔️:         size = localSize;
    // RDKit✔️✔️:       } else if (size != localSize) {
    // RDKit✔️✔️:         throw std::runtime_error("Something unexpected!");
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     for (int i = 0; i < size; ++i) {
    // RDKit✔️✔️:       std::vector<const Node *> eq;
    // RDKit✔️✔️:       for (const auto &aTmp : tmp) {
    // RDKit✔️✔️:         auto tmpNodes = toNodeList(aTmp[i]);
    // RDKit✔️✔️:         eq.insert(eq.end(), tmpNodes.begin(), tmpNodes.end());
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (!eq.empty()) {
    // RDKit✔️✔️:         nextLevel.push_back(eq);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return nextLevel;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Rule4b::getNextLevel
    fn get_next_level(
        &self,
        sort_rules: Option<&[&dyn CipSequenceRule]>,
        digraph: &mut CipDigraph<'_>,
        context: &mut CipLabelerContext,
        prev_level: &[Vec<CipNodeId>],
    ) -> Result<Vec<Vec<CipNodeId>>, CipLabelerError> {
        let mut next_level = Vec::with_capacity(4 * prev_level.len());
        for prev in prev_level {
            let mut tmp = Vec::<Vec<Vec<CipEdgeId>>>::new();
            for node in prev {
                let mut edges = digraph.non_terminal_out_edges(*node)?;
                sort_edges_for_sequence_rule(
                    self, sort_rules, digraph, context, *node, &mut edges, true,
                )?;
                let groups = if let Some(sort_rules) = sort_rules {
                    CipSort::from_rules(sort_rules.to_vec()).get_groups(digraph, context, &edges)?
                } else {
                    CipSort::new(self).get_groups(digraph, context, &edges)?
                };
                tmp.push(groups);
            }

            let mut size = -1_i32;
            for groups in &tmp {
                let local_size = groups.len() as i32;
                if size < 0 {
                    size = local_size;
                } else if size != local_size {
                    return Err(CipLabelerError::SomethingUnexpected);
                }
            }

            for i in 0..usize::try_from(size.max(0)).expect("nonnegative") {
                let mut eq = Vec::new();
                for a_tmp in &tmp {
                    let tmp_nodes = self.to_node_list(digraph, &a_tmp[i]);
                    eq.extend(tmp_nodes);
                }
                if !eq.is_empty() {
                    next_level.push(eq);
                }
            }
        }
        Ok(next_level)
    }

    // BEGIN RDKIT CPP FUNCTION Rule4b::toNodeList (rules/Rule4b.cpp)
    // RDKit✔️✔️: std::vector<const Node *> Rule4b::toNodeList(
    // RDKit✔️✔️:     const std::vector<Edge *> &eqEdges) const {
    // RDKit✔️✔️:   std::vector<const Node *> eqNodes;
    // RDKit✔️✔️:   eqNodes.reserve(eqEdges.size());
    // RDKit✔️✔️:   for (const auto &edge : eqEdges) {
    // RDKit✔️✔️:     eqNodes.push_back(edge->getEnd());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return eqNodes;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Rule4b::toNodeList
    fn to_node_list(&self, digraph: &CipDigraph<'_>, eq_edges: &[CipEdgeId]) -> Vec<CipNodeId> {
        let mut eq_nodes = Vec::with_capacity(eq_edges.len());
        for edge in eq_edges {
            eq_nodes.push(digraph.edge(*edge).get_end());
        }
        eq_nodes
    }

    // BEGIN RDKIT CPP FUNCTION Rule4b::newPairLists (rules/Rule4b.cpp)
    // RDKit✔️✔️: std::vector<PairList> Rule4b::newPairLists(
    // RDKit✔️✔️:     const std::vector<Descriptor> &descriptors) const {
    // RDKit✔️✔️:   std::vector<PairList> pairs;
    // RDKit✔️✔️:   pairs.reserve(descriptors.size());
    // RDKit✔️✔️:   for (Descriptor descriptor : descriptors) {
    // RDKit✔️✔️:     pairs.emplace_back(descriptor);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return pairs;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Rule4b::newPairLists
    fn new_pair_lists(&self, descriptors: &[Descriptor]) -> Vec<CipPairList> {
        let mut pairs = Vec::with_capacity(descriptors.len());
        for descriptor in descriptors {
            pairs.push(CipPairList::with_ref(*descriptor));
        }
        pairs
    }

    // BEGIN RDKIT CPP FUNCTION Rule4b::fillPairs (rules/Rule4b.cpp)
    // RDKit✔️✔️: void Rule4b::fillPairs(const Node *beg, PairList &plist) const {
    // RDKit✔️✔️:   const Rule4b replacement_rule(plist.getRefDescriptor());
    // RDKit✔️✔️:   const auto &sorter = getRefSorter(&replacement_rule);
    // RDKit✔️✔️:   auto queue = std::list<const Node *>({beg});
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (const auto &node : queue) {
    // RDKit✔️✔️:     plist.add(node->getAux());
    // RDKit✔️✔️:     auto edges = node->getEdges();
    // RDKit✔️✔️:     sorter.prioritize(node, edges);
    // RDKit✔️✔️:     for (const auto &edge : edges) {
    // RDKit✔️✔️:       if (edge->isBeg(node) && !edge->getEnd()->isTerminal()) {
    // RDKit✔️✔️:         queue.push_back(edge->getEnd());
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Rule4b::fillPairs
    fn fill_pairs(
        &self,
        sort_rules: Option<&[&dyn CipSequenceRule]>,
        digraph: &mut CipDigraph<'_>,
        context: &mut CipLabelerContext,
        beg: CipNodeId,
        plist: &mut CipPairList,
    ) -> Result<(), CipLabelerError> {
        let replacement_rule = CipRule4b::with_ref(plist.get_ref_descriptor());
        let ref_sort_rules = self.get_ref_sorter(sort_rules, &replacement_rule)?;
        let sorter = CipSort::from_rules(ref_sort_rules);
        let mut queue = vec![beg];
        let mut pos = 0_usize;
        while pos < queue.len() {
            let node = queue[pos];
            plist.add(digraph.node(node).get_aux());
            let mut edges = digraph.node_edges(node)?;
            sorter.prioritize(digraph, context, node, &mut edges, true)?;
            for edge in edges {
                if digraph.edge(edge).is_beg(node)
                    && !digraph.node(digraph.edge(edge).get_end()).is_terminal()
                {
                    queue.push(digraph.edge(edge).get_end());
                }
            }
            pos += 1;
        }
        Ok(())
    }

    // BEGIN RDKIT CPP FUNCTION Rule4b::comparePairs (rules/Rule4b.cpp)
    // RDKit✔️✔️: int Rule4b::comparePairs(const Node *a, const Node *b, Descriptor refA,
    // RDKit✔️✔️:                          Descriptor refB) const {
    // RDKit✔️✔️:   const Rule4b replacementA(refA);
    // RDKit✔️✔️:   const Rule4b replacementB(refB);
    // RDKit✔️✔️:   const auto &aSorter = getRefSorter(&replacementA);
    // RDKit✔️✔️:   const auto &bSorter = getRefSorter(&replacementB);
    // RDKit✔️✔️:   auto aQueue = std::vector<const Node *>({a});
    // RDKit✔️✔️:   auto bQueue = std::vector<const Node *>({b});
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (auto pos = 0u; pos < aQueue.size() && pos < bQueue.size(); ++pos) {
    // RDKit✔️✔️:     const auto aNode = aQueue[pos];
    // RDKit✔️✔️:     const auto bNode = bQueue[pos];
    // RDKit✔️✔️:
    // RDKit✔️✔️:     const auto &desA = PairList::ref(aNode->getAux());
    // RDKit✔️✔️:     const auto &desB = PairList::ref(bNode->getAux());
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (desA == refA && desB != refB) {
    // RDKit✔️✔️:       return +1;
    // RDKit✔️✔️:     } else if (desA != refA && desB == refB) {
    // RDKit✔️✔️:       return -1;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     auto edges = aNode->getEdges();
    // RDKit✔️✔️:     aSorter.prioritize(aNode, edges);
    // RDKit✔️✔️:     for (const auto &edge : edges) {
    // RDKit✔️✔️:       if (edge->isBeg(aNode) && !edge->getEnd()->isTerminal()) {
    // RDKit✔️✔️:         aQueue.push_back(edge->getEnd());
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     edges = bNode->getEdges();
    // RDKit✔️✔️:     bSorter.prioritize(bNode, edges);
    // RDKit✔️✔️:     for (const auto &edge : edges) {
    // RDKit✔️✔️:       if (edge->isBeg(bNode) && !edge->getEnd()->isTerminal()) {
    // RDKit✔️✔️:         bQueue.push_back(edge->getEnd());
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return 0;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Rule4b::comparePairs
    fn compare_pairs(
        &self,
        sort_rules: Option<&[&dyn CipSequenceRule]>,
        digraph: &mut CipDigraph<'_>,
        context: &mut CipLabelerContext,
        a: CipNodeId,
        b: CipNodeId,
        ref_a: Descriptor,
        ref_b: Descriptor,
    ) -> Result<i32, CipLabelerError> {
        let replacement_a = CipRule4b::with_ref(ref_a);
        let replacement_b = CipRule4b::with_ref(ref_b);
        let a_sorter = CipSort::from_rules(self.get_ref_sorter(sort_rules, &replacement_a)?);
        let b_sorter = CipSort::from_rules(self.get_ref_sorter(sort_rules, &replacement_b)?);
        let mut a_queue = vec![a];
        let mut b_queue = vec![b];
        let mut pos = 0_usize;
        while pos < a_queue.len() && pos < b_queue.len() {
            let a_node = a_queue[pos];
            let b_node = b_queue[pos];
            let des_a = CipPairList::ref_descriptor(digraph.node(a_node).get_aux());
            let des_b = CipPairList::ref_descriptor(digraph.node(b_node).get_aux());

            if des_a == ref_a && des_b != ref_b {
                return Ok(1);
            } else if des_a != ref_a && des_b == ref_b {
                return Ok(-1);
            }

            let mut edges = digraph.node_edges(a_node)?;
            a_sorter.prioritize(digraph, context, a_node, &mut edges, true)?;
            for edge in edges {
                if digraph.edge(edge).is_beg(a_node)
                    && !digraph.node(digraph.edge(edge).get_end()).is_terminal()
                {
                    a_queue.push(digraph.edge(edge).get_end());
                }
            }

            let mut edges = digraph.node_edges(b_node)?;
            b_sorter.prioritize(digraph, context, b_node, &mut edges, true)?;
            for edge in edges {
                if digraph.edge(edge).is_beg(b_node)
                    && !digraph.node(digraph.edge(edge).get_end()).is_terminal()
                {
                    b_queue.push(digraph.edge(edge).get_end());
                }
            }
            pos += 1;
        }
        Ok(0)
    }

    // BEGIN RDKIT CPP FUNCTION Rule4b::getRefSorter (rules/Rule4b.cpp)
    // RDKit✔️✔️: Sort Rule4b::getRefSorter(const SequenceRule *replacement_rule) const {
    // RDKit✔️✔️:   const auto &rules = getSorter()->getRules();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   CHECK_INVARIANT(std::find(rules.begin(), rules.end(), this) != rules.end(),
    // RDKit✔️✔️:                   "Rule4b instance not in rule set");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::vector<const SequenceRule *> new_rules;
    // RDKit✔️✔️:   new_rules.reserve(rules.size());
    // RDKit✔️✔️:   for (const auto &rule : rules) {
    // RDKit✔️✔️:     if (this != rule) {
    // RDKit✔️✔️:       new_rules.push_back(rule);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   new_rules.push_back(replacement_rule);
    // RDKit✔️✔️:   return {new_rules};
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Rule4b::getRefSorter
    fn get_ref_sorter<'r>(
        &'r self,
        sort_rules: Option<&'r [&'r dyn CipSequenceRule]>,
        replacement_rule: &'r dyn CipSequenceRule,
    ) -> Result<Vec<&'r dyn CipSequenceRule>, CipLabelerError> {
        let mut new_rules = Vec::new();
        if let Some(sort_rules) = sort_rules {
            let self_ptr = self as &dyn CipSequenceRule as *const dyn CipSequenceRule;
            let mut found = false;
            for rule in sort_rules {
                let rule_ptr = *rule as *const dyn CipSequenceRule;
                if std::ptr::addr_eq(rule_ptr, self_ptr) {
                    found = true;
                    break;
                } else {
                    new_rules.push(*rule);
                }
            }
            if !found {
                return Err(CipLabelerError::Rule4bInstanceNotInRuleSet);
            }
        }
        new_rules.push(replacement_rule);
        Ok(new_rules)
    }
}

impl CipSequenceRule for CipRule4b {
    fn compare(
        &self,
        digraph: &mut CipDigraph<'_>,
        context: &mut CipLabelerContext,
        a: CipEdgeId,
        b: CipEdgeId,
    ) -> Result<i32, CipLabelerError> {
        self.compare_with_sort_rules(None, digraph, context, a, b)
    }

    // BEGIN RDKIT CPP FUNCTION Rule4b::compare (rules/Rule4b.cpp)
    // RDKit✔️✔️: Rule4b::Rule4b() = default;
    // RDKit✔️✔️:
    // RDKit✔️✔️: Rule4b::Rule4b(Descriptor ref) : d_ref{ref} {}
    // RDKit✔️✔️:
    // RDKit✔️✔️: int Rule4b::compare(const Edge *a, const Edge *b) const {
    // RDKit✔️✔️:   const auto &aBeg = a->getBeg();
    // RDKit✔️✔️:   const auto &aEnd = a->getEnd();
    // RDKit✔️✔️:   const auto &bBeg = b->getBeg();
    // RDKit✔️✔️:   const auto &bEnd = b->getEnd();
    // RDKit✔️✔️:   if (aBeg->getDigraph()->getCurrentRoot() != aBeg ||
    // RDKit✔️✔️:       bBeg->getDigraph()->getCurrentRoot() != bBeg) {
    // RDKit✔️✔️:     if (d_ref == Descriptor::NONE) {
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     Descriptor aDesc = aEnd->getAux();
    // RDKit✔️✔️:     Descriptor bDesc = bEnd->getAux();
    // RDKit✔️✔️:     if (aDesc != Descriptor::NONE && bDesc != Descriptor::NONE &&
    // RDKit✔️✔️:         aDesc != Descriptor::ns && bDesc != Descriptor::ns) {
    // RDKit✔️✔️:       bool alike = PairList::ref(d_ref) == PairList::ref(aDesc);
    // RDKit✔️✔️:       bool blike = PairList::ref(d_ref) == PairList::ref(bDesc);
    // RDKit✔️✔️:       if (alike && !blike) {
    // RDKit✔️✔️:         return +1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (blike && !alike) {
    // RDKit✔️✔️:         return -1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     auto list1 = newPairLists(getReferenceDescriptors(aEnd));
    // RDKit✔️✔️:
    // RDKit✔️✔️:     auto list2 = newPairLists(getReferenceDescriptors(bEnd));
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (list1.empty() != list2.empty()) {
    // RDKit✔️✔️:       throw std::runtime_error(
    // RDKit✔️✔️:           "Substituents should be topologically equivalent!");
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (list1.size() == 1) {
    // RDKit✔️✔️:       return comparePairs(aEnd, bEnd, list1[0].getRefDescriptor(),
    // RDKit✔️✔️:                           list2[0].getRefDescriptor());
    // RDKit✔️✔️:     } else if (list1.size() > 1) {
    // RDKit✔️✔️:       for (auto &plist : list1) {
    // RDKit✔️✔️:         fillPairs(aEnd, plist);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       for (auto &plist : list2) {
    // RDKit✔️✔️:         fillPairs(bEnd, plist);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       std::sort(list1.rbegin(), list1.rend());
    // RDKit✔️✔️:       std::sort(list2.rbegin(), list2.rend());
    // RDKit✔️✔️:
    // RDKit✔️✔️:       for (auto i = 0u; i < list1.size(); ++i) {
    // RDKit✔️✔️:         int cmp = list1[i].compareTo(list2[i]);
    // RDKit✔️✔️:         if (cmp != 0) {
    // RDKit✔️✔️:           return cmp;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Rule4b::compare
    fn compare_with_sort_rules(
        &self,
        sort_rules: Option<&[&dyn CipSequenceRule]>,
        digraph: &mut CipDigraph<'_>,
        context: &mut CipLabelerContext,
        a: CipEdgeId,
        b: CipEdgeId,
    ) -> Result<i32, CipLabelerError> {
        let a_beg = digraph.edge(a).get_beg();
        let a_end = digraph.edge(a).get_end();
        let b_beg = digraph.edge(b).get_beg();
        let b_end = digraph.edge(b).get_end();

        if digraph.get_current_root() != a_beg || digraph.get_current_root() != b_beg {
            if self.ref_descriptor == Descriptor::None {
                return Ok(0);
            }
            let a_desc = digraph.node(a_end).get_aux();
            let b_desc = digraph.node(b_end).get_aux();
            if a_desc != Descriptor::None
                && b_desc != Descriptor::None
                && a_desc != Descriptor::ns
                && b_desc != Descriptor::ns
            {
                let alike = CipPairList::ref_descriptor(self.ref_descriptor)
                    == CipPairList::ref_descriptor(a_desc);
                let blike = CipPairList::ref_descriptor(self.ref_descriptor)
                    == CipPairList::ref_descriptor(b_desc);
                if alike && !blike {
                    return Ok(1);
                }
                if blike && !alike {
                    return Ok(-1);
                }
            }
            return Ok(0);
        }

        let mut list1 = self
            .new_pair_lists(&self.get_reference_descriptors(sort_rules, digraph, context, a_end)?);
        let mut list2 = self
            .new_pair_lists(&self.get_reference_descriptors(sort_rules, digraph, context, b_end)?);

        if list1.is_empty() != list2.is_empty() {
            return Err(CipLabelerError::SubstituentsShouldBeTopologicallyEquivalent);
        }
        if list1.len() == 1 {
            self.compare_pairs(
                sort_rules,
                digraph,
                context,
                a_end,
                b_end,
                list1[0].get_ref_descriptor(),
                list2[0].get_ref_descriptor(),
            )
        } else if list1.len() > 1 {
            for plist in &mut list1 {
                self.fill_pairs(sort_rules, digraph, context, a_end, plist)?;
            }
            for plist in &mut list2 {
                self.fill_pairs(sort_rules, digraph, context, b_end, plist)?;
            }

            CipPairList::sort_descending(&mut list1)?;
            CipPairList::sort_descending(&mut list2)?;

            for (left, right) in list1.iter().zip(list2.iter()) {
                let cmp = left.compare_to(right)?;
                if cmp != 0 {
                    return Ok(cmp);
                }
            }
            Ok(0)
        } else {
            Ok(0)
        }
    }
}

#[derive(Debug, Default, Clone, Copy)]
pub(crate) struct CipRule4c;

impl CipRule4c {
    // BEGIN RDKIT CPP FUNCTION ord (rules/Rule4c.cpp)
    // RDKit✔️✔️: int ord(Descriptor lab) {
    // RDKit✔️✔️:   switch (lab) {
    // RDKit✔️✔️:     case Descriptor::m:
    // RDKit✔️✔️:     case Descriptor::r:
    // RDKit✔️✔️:       return 2;
    // RDKit✔️✔️:     case Descriptor::p:
    // RDKit✔️✔️:     case Descriptor::s:
    // RDKit✔️✔️:       return 1;
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION ord
    fn ord(lab: Descriptor) -> i32 {
        match lab {
            Descriptor::m | Descriptor::r => 2,
            Descriptor::p | Descriptor::s => 1,
            _ => 0,
        }
    }
}

impl CipSequenceRule for CipRule4c {
    // BEGIN RDKIT CPP FUNCTION Rule4c::compare (rules/Rule4c.cpp)
    // RDKit✔️✔️: Rule4c::Rule4c() = default;
    // RDKit✔️✔️:
    // RDKit✔️✔️: int Rule4c::compare(const Edge *a, const Edge *b) const {
    // RDKit✔️✔️:   // m vs p
    // RDKit✔️✔️:   int aOrdinal = ord(getBondLabel(a));
    // RDKit✔️✔️:   int bOrdinal = ord(getBondLabel(b));
    // RDKit✔️✔️:   int cmp = three_way_comparison(aOrdinal, bOrdinal);
    // RDKit✔️✔️:   if (cmp != 0) {
    // RDKit✔️✔️:     return cmp;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // r vs s
    // RDKit✔️✔️:   aOrdinal = ord(a->getEnd()->getAux());
    // RDKit✔️✔️:   bOrdinal = ord(b->getEnd()->getAux());
    // RDKit✔️✔️:   return three_way_comparison(aOrdinal, bOrdinal);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Rule4c::compare
    fn compare(
        &self,
        digraph: &mut CipDigraph<'_>,
        _context: &mut CipLabelerContext,
        a: CipEdgeId,
        b: CipEdgeId,
    ) -> Result<i32, CipLabelerError> {
        let a_ordinal = Self::ord(self.get_bond_label(digraph.edge(a)));
        let b_ordinal = Self::ord(self.get_bond_label(digraph.edge(b)));
        let cmp = three_way_comparison_i32(a_ordinal, b_ordinal);
        if cmp != 0 {
            return Ok(cmp);
        }
        let a_ordinal = Self::ord(digraph.node(digraph.edge(a).get_end()).get_aux());
        let b_ordinal = Self::ord(digraph.node(digraph.edge(b).get_end()).get_aux());
        Ok(three_way_comparison_i32(a_ordinal, b_ordinal))
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct CipRule5New {
    ref_descriptor: Descriptor,
}

impl CipRule5New {
    pub(crate) fn new() -> Self {
        Self {
            ref_descriptor: Descriptor::None,
        }
    }

    pub(crate) fn with_ref(ref_descriptor: Descriptor) -> Self {
        Self { ref_descriptor }
    }

    // BEGIN RDKIT CPP FUNCTION Rule5New::fillPairs (rules/Rule5New.cpp)
    // RDKit✔️✔️: void Rule5New::fillPairs(const Node *beg, PairList &plist) const {
    // RDKit✔️✔️:   const Rule5New replacement_rule(plist.getRefDescriptor());
    // RDKit✔️✔️:   const auto &sorter = getRefSorter(&replacement_rule);
    // RDKit✔️✔️:   auto queue = std::list<const Node *>({beg});
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (const auto &node : queue) {
    // RDKit✔️✔️:     plist.add(node->getAux());
    // RDKit✔️✔️:     auto edges = node->getEdges();
    // RDKit✔️✔️:     sorter.prioritize(node, edges);
    // RDKit✔️✔️:     for (const auto &edge : edges) {
    // RDKit✔️✔️:       if (edge->isBeg(node) && !edge->getEnd()->isTerminal()) {
    // RDKit✔️✔️:         queue.push_back(edge->getEnd());
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Rule5New::fillPairs
    fn fill_pairs(
        &self,
        sort_rules: Option<&[&dyn CipSequenceRule]>,
        digraph: &mut CipDigraph<'_>,
        context: &mut CipLabelerContext,
        beg: CipNodeId,
        plist: &mut CipPairList,
    ) -> Result<(), CipLabelerError> {
        let replacement_rule = CipRule5New::with_ref(plist.get_ref_descriptor());
        let ref_sort_rules = self.get_ref_sorter(sort_rules, &replacement_rule)?;
        let sorter = CipSort::from_rules(ref_sort_rules);
        let mut queue = vec![beg];
        let mut pos = 0_usize;
        while pos < queue.len() {
            let node = queue[pos];
            plist.add(digraph.node(node).get_aux());
            let mut edges = digraph.node_edges(node)?;
            sorter.prioritize(digraph, context, node, &mut edges, true)?;
            for edge in edges {
                if digraph.edge(edge).is_beg(node)
                    && !digraph.node(digraph.edge(edge).get_end()).is_terminal()
                {
                    queue.push(digraph.edge(edge).get_end());
                }
            }
            pos += 1;
        }
        Ok(())
    }

    // BEGIN RDKIT CPP FUNCTION Rule5New::getRefSorter (rules/Rule5New.cpp)
    // RDKit✔️✔️: Sort Rule5New::getRefSorter(const SequenceRule *replacement_rule) const {
    // RDKit✔️✔️:   const auto &rules = getSorter()->getRules();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   CHECK_INVARIANT(std::find(rules.begin(), rules.end(), this) != rules.end(),
    // RDKit✔️✔️:                   "Rule5New instance not in rule set");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::vector<const SequenceRule *> new_rules;
    // RDKit✔️✔️:   new_rules.reserve(rules.size());
    // RDKit✔️✔️:   for (const auto &rule : rules) {
    // RDKit✔️✔️:     if (this != rule) {
    // RDKit✔️✔️:       new_rules.push_back(rule);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   new_rules.push_back(replacement_rule);
    // RDKit✔️✔️:   return {new_rules};
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Rule5New::getRefSorter
    fn get_ref_sorter<'r>(
        &'r self,
        sort_rules: Option<&'r [&'r dyn CipSequenceRule]>,
        replacement_rule: &'r dyn CipSequenceRule,
    ) -> Result<Vec<&'r dyn CipSequenceRule>, CipLabelerError> {
        let mut new_rules = Vec::new();
        if let Some(sort_rules) = sort_rules {
            let self_ptr = self as &dyn CipSequenceRule as *const dyn CipSequenceRule;
            let mut found = false;
            for rule in sort_rules {
                let rule_ptr = *rule as *const dyn CipSequenceRule;
                if std::ptr::addr_eq(rule_ptr, self_ptr) {
                    found = true;
                    break;
                } else {
                    new_rules.push(*rule);
                }
            }
            if !found {
                return Err(CipLabelerError::Rule5NewInstanceNotInRuleSet);
            }
        }
        new_rules.push(replacement_rule);
        Ok(new_rules)
    }
}

impl CipSequenceRule for CipRule5New {
    fn compare(
        &self,
        digraph: &mut CipDigraph<'_>,
        context: &mut CipLabelerContext,
        a: CipEdgeId,
        b: CipEdgeId,
    ) -> Result<i32, CipLabelerError> {
        self.compare_with_sort_rules(None, digraph, context, a, b)
    }

    // BEGIN RDKIT CPP FUNCTION Rule5New::compare (rules/Rule5New.cpp)
    // RDKit✔️✔️: Rule5New::Rule5New() = default;
    // RDKit✔️✔️:
    // RDKit✔️✔️: Rule5New::Rule5New(Descriptor ref) : d_ref{ref} {}
    // RDKit✔️✔️:
    // RDKit✔️✔️: int Rule5New::compare(const Edge *a, const Edge *b) const {
    // RDKit✔️✔️:   const auto &aBeg = a->getBeg();
    // RDKit✔️✔️:   const auto &aEnd = a->getEnd();
    // RDKit✔️✔️:   const auto &bBeg = b->getBeg();
    // RDKit✔️✔️:   const auto &bEnd = b->getEnd();
    // RDKit✔️✔️:   if (aBeg->getDigraph()->getCurrentRoot() != aBeg ||
    // RDKit✔️✔️:       bBeg->getDigraph()->getCurrentRoot() != bBeg) {
    // RDKit✔️✔️:     if (d_ref == Descriptor::NONE) {
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     Descriptor aDesc = aEnd->getAux();
    // RDKit✔️✔️:     Descriptor bDesc = bEnd->getAux();
    // RDKit✔️✔️:     if (aDesc != Descriptor::NONE && bDesc != Descriptor::NONE &&
    // RDKit✔️✔️:         aDesc != Descriptor::ns && bDesc != Descriptor::ns) {
    // RDKit✔️✔️:       bool alike = PairList::ref(d_ref) == PairList::ref(aDesc);
    // RDKit✔️✔️:       bool blike = PairList::ref(d_ref) == PairList::ref(bDesc);
    // RDKit✔️✔️:       if (alike && !blike) {
    // RDKit✔️✔️:         return +1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (blike && !alike) {
    // RDKit✔️✔️:         return -1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     auto listRA = PairList(Descriptor::R);
    // RDKit✔️✔️:     auto listRB = PairList(Descriptor::R);
    // RDKit✔️✔️:     auto listSA = PairList(Descriptor::S);
    // RDKit✔️✔️:     auto listSB = PairList(Descriptor::S);
    // RDKit✔️✔️:     fillPairs(aEnd, listRA);
    // RDKit✔️✔️:     fillPairs(aEnd, listSA);
    // RDKit✔️✔️:     fillPairs(bEnd, listRB);
    // RDKit✔️✔️:     fillPairs(bEnd, listSB);
    // RDKit✔️✔️:     int cmpR = listRA.compareTo(listRB);
    // RDKit✔️✔️:     int cmpS = listSA.compareTo(listSB);
    // RDKit✔️✔️:     // -2/+2 for pseudo-asymetric
    // RDKit✔️✔️:     // -1/+1 if not (e.g. the R > R and S > S lists)
    // RDKit✔️✔️:     if (cmpR < 0) {
    // RDKit✔️✔️:       return cmpS < 0 ? -1 : -2;
    // RDKit✔️✔️:     } else if (cmpR > 0) {
    // RDKit✔️✔️:       return cmpS > 0 ? +1 : +2;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Rule5New::compare
    fn compare_with_sort_rules(
        &self,
        sort_rules: Option<&[&dyn CipSequenceRule]>,
        digraph: &mut CipDigraph<'_>,
        context: &mut CipLabelerContext,
        a: CipEdgeId,
        b: CipEdgeId,
    ) -> Result<i32, CipLabelerError> {
        let a_beg = digraph.edge(a).get_beg();
        let a_end = digraph.edge(a).get_end();
        let b_beg = digraph.edge(b).get_beg();
        let b_end = digraph.edge(b).get_end();

        if digraph.get_current_root() != a_beg || digraph.get_current_root() != b_beg {
            if self.ref_descriptor == Descriptor::None {
                return Ok(0);
            }
            let a_desc = digraph.node(a_end).get_aux();
            let b_desc = digraph.node(b_end).get_aux();
            if a_desc != Descriptor::None
                && b_desc != Descriptor::None
                && a_desc != Descriptor::ns
                && b_desc != Descriptor::ns
            {
                let alike = CipPairList::ref_descriptor(self.ref_descriptor)
                    == CipPairList::ref_descriptor(a_desc);
                let blike = CipPairList::ref_descriptor(self.ref_descriptor)
                    == CipPairList::ref_descriptor(b_desc);
                if alike && !blike {
                    return Ok(1);
                }
                if blike && !alike {
                    return Ok(-1);
                }
            }
            return Ok(0);
        }

        let mut list_ra = CipPairList::with_ref(Descriptor::R);
        let mut list_rb = CipPairList::with_ref(Descriptor::R);
        let mut list_sa = CipPairList::with_ref(Descriptor::S);
        let mut list_sb = CipPairList::with_ref(Descriptor::S);
        self.fill_pairs(sort_rules, digraph, context, a_end, &mut list_ra)?;
        self.fill_pairs(sort_rules, digraph, context, a_end, &mut list_sa)?;
        self.fill_pairs(sort_rules, digraph, context, b_end, &mut list_rb)?;
        self.fill_pairs(sort_rules, digraph, context, b_end, &mut list_sb)?;
        let cmp_r = list_ra.compare_to(&list_rb)?;
        let cmp_s = list_sa.compare_to(&list_sb)?;
        if cmp_r < 0 {
            Ok(if cmp_s < 0 { -1 } else { -2 })
        } else if cmp_r > 0 {
            Ok(if cmp_s > 0 { 1 } else { 2 })
        } else {
            Ok(0)
        }
    }
}

#[derive(Debug, Default, Clone, Copy)]
pub(crate) struct CipRule6;

impl CipSequenceRule for CipRule6 {
    // BEGIN RDKIT CPP FUNCTION Rule6::compare (rules/Rule6.cpp)
    // RDKit✔️✔️: Rule6::Rule6() = default;
    // RDKit✔️✔️:
    // RDKit✔️✔️: int Rule6::compare(const Edge *a, const Edge *b) const {
    // RDKit✔️✔️:   const auto &digraph = a->getBeg()->getDigraph();
    // RDKit✔️✔️:   const auto &ref = digraph->getRule6Ref();
    // RDKit✔️✔️:   if (ref == nullptr) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   const auto &aAtom = a->getEnd()->getAtom();
    // RDKit✔️✔️:   const auto &bAtom = b->getEnd()->getAtom();
    // RDKit✔️✔️:   if (ref == aAtom && ref != bAtom) {
    // RDKit✔️✔️:     return +1;  // a is ref (has priority)
    // RDKit✔️✔️:   } else if (ref != aAtom && ref == bAtom) {
    // RDKit✔️✔️:     return -1;  // b is ref (has priority)
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return 0;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Rule6::compare
    fn compare(
        &self,
        digraph: &mut CipDigraph<'_>,
        _context: &mut CipLabelerContext,
        a: CipEdgeId,
        b: CipEdgeId,
    ) -> Result<i32, CipLabelerError> {
        let Some(ref_atom) = digraph.get_rule6_ref() else {
            return Ok(0);
        };
        let a_atom = digraph.node(digraph.edge(a).get_end()).atom_idx();
        let b_atom = digraph.node(digraph.edge(b).get_end()).atom_idx();
        if Some(ref_atom) == a_atom && Some(ref_atom) != b_atom {
            Ok(1)
        } else if Some(ref_atom) != a_atom && Some(ref_atom) == b_atom {
            Ok(-1)
        } else {
            Ok(0)
        }
    }
}

impl CipNode {
    // BEGIN RDKIT CPP CONSTANTS Node flags (CIPLabeler/Node.h)
    // RDKit✔️✔️: static const int EXPANDED = 0x1;
    // RDKit✔️✔️: static const int RING_DUPLICATE = 0x2;
    // RDKit✔️✔️: static const int BOND_DUPLICATE = 0x4;
    // RDKit✔️✔️: static const int DUPLICATE = RING_DUPLICATE | BOND_DUPLICATE;
    // RDKit✔️✔️: static const int IMPL_HYDROGEN = 0x8;
    // RDKit✔️✔️: static const int DUPLICATE_OR_H =
    // RDKit✔️✔️:     RING_DUPLICATE | BOND_DUPLICATE | IMPL_HYDROGEN;
    // END RDKIT CPP CONSTANTS Node flags
    pub(crate) const EXPANDED: i32 = 0x1;
    pub(crate) const RING_DUPLICATE: i32 = 0x2;
    pub(crate) const BOND_DUPLICATE: i32 = 0x4;
    pub(crate) const DUPLICATE: i32 = Self::RING_DUPLICATE | Self::BOND_DUPLICATE;
    pub(crate) const IMPL_HYDROGEN: i32 = 0x8;
    pub(crate) const DUPLICATE_OR_H: i32 =
        Self::RING_DUPLICATE | Self::BOND_DUPLICATE | Self::IMPL_HYDROGEN;
    pub(crate) const NO_ATOM_INDEX: CipSourceIndex = CIP_NO_ATOM;

    // BEGIN RDKIT CPP FUNCTION Node::Node (CIPLabeler/Node.cpp)
    // RDKit✔️✔️: Node::Node(Digraph *g, std::vector<char> &&visit, Atom *atom,
    // RDKit✔️✔️:            boost::rational<int> &&frac, int dist, int flags)
    // RDKit✔️✔️:     : dp_g{g},
    // RDKit✔️✔️:       dp_atom{atom},
    // RDKit✔️✔️:       d_dist{dist},
    // RDKit✔️✔️:       d_atomic_num{std::move(frac)},
    // RDKit✔️✔️:       d_flags{flags},
    // RDKit✔️✔️:       d_visit{std::move(visit)} {
    // RDKit✔️✔️:   if (d_flags & DUPLICATE) {
    // RDKit✔️✔️:     d_edges.reserve(4);
    // RDKit✔️✔️:     d_atomic_mass = 0.;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     const auto &table = RDKit::PeriodicTable::getTable();
    // RDKit✔️✔️:     auto atomic_number = getAtomicNum();
    // RDKit✔️✔️:     auto isotope = getMassNum();
    // RDKit✔️✔️:     if (isotope == 0u) {
    // RDKit✔️✔️:       d_atomic_mass = table->getAtomicWeight(atomic_number);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       d_atomic_mass = table->getMassForIsotope(atomic_number, isotope);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (d_visit.empty() || d_flags & DUPLICATE) {
    // RDKit✔️✔️:     d_flags |= EXPANDED;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Node::Node
    pub(crate) fn new(
        digraph: usize,
        visit: Vec<i8>,
        atom_idx: Option<usize>,
        frac: RationalI32,
        dist: i32,
        flags: i32,
        mol: &CipMol<'_>,
    ) -> Result<Self, CipLabelerError> {
        let mut flags = flags;
        let (edge_capacity, atomic_mass) = if flags & Self::DUPLICATE != 0 {
            (4, 0.0)
        } else {
            let atomic_number = atom_idx
                .map(|idx| mol.atom(idx).map(Atom::atomic_number))
                .transpose()?
                .unwrap_or(1);
            let isotope = atom_idx
                .map(|idx| mol.atom(idx).map(Atom::isotope))
                .transpose()?
                .flatten();
            (0, rdkit_atomic_mass(atomic_number, isotope)?)
        };
        if visit.is_empty() || flags & Self::DUPLICATE != 0 {
            flags |= Self::EXPANDED;
        }
        Ok(Self {
            digraph,
            atom_idx,
            distance: dist,
            atomic_num_fraction: frac,
            atomic_mass,
            aux: Descriptor::None,
            flags,
            edges: Vec::with_capacity(edge_capacity),
            visit,
        })
    }

    // BEGIN RDKIT CPP FUNCTION Node::newTerminalChild (CIPLabeler/Node.cpp)
    // RDKit✔️✔️: Node *Node::newTerminalChild(int idx, Atom *atom, int flags) const {
    // RDKit✔️✔️:   int new_dist = flags & DUPLICATE ? d_visit[idx] : d_dist + 1;
    // RDKit✔️✔️:   std::vector<char> new_visit;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (flags & BOND_DUPLICATE) {
    // RDKit✔️✔️:     auto frac = dp_g->getMol().getFractionalAtomicNum(dp_atom);
    // RDKit✔️✔️:     if (frac.denominator() > 1) {
    // RDKit✔️✔️:       return &dp_g->addNode(std::move(new_visit), atom, std::move(frac),
    // RDKit✔️✔️:                             new_dist, flags);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto atomic_num = atom ? atom->getAtomicNum() : 1;
    // RDKit✔️✔️:   return &dp_g->addNode(std::move(new_visit), atom, atomic_num, new_dist,
    // RDKit✔️✔️:                         flags);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Node::newTerminalChild
    fn new_terminal_child(
        &self,
        idx: Option<usize>,
        atom_idx: Option<usize>,
        flags: i32,
        mol: &mut CipMol<'_>,
    ) -> Result<Self, CipLabelerError> {
        let new_dist = if flags & Self::DUPLICATE != 0 {
            let idx = idx.ok_or(CipLabelerError::InvalidInternalState {
                detail: "Node::newTerminalChild duplicate without atom index",
            })?;
            i32::from(self.visit[idx])
        } else {
            self.distance + 1
        };
        let new_visit = Vec::new();

        if flags & Self::BOND_DUPLICATE != 0 {
            let current_atom = self.atom_idx.ok_or(CipLabelerError::InvalidInternalState {
                detail: "Node::newTerminalChild bond duplicate from null atom",
            })?;
            let frac = mol.get_fractional_atomic_num(current_atom)?;
            if frac.denominator > 1 {
                return Self::new(
                    self.digraph,
                    new_visit,
                    atom_idx,
                    frac,
                    new_dist,
                    flags,
                    mol,
                );
            }
        }

        let atomic_num = atom_idx
            .map(|idx| mol.atom(idx).map(Atom::atomic_number))
            .transpose()?
            .unwrap_or(1);
        Self::new(
            self.digraph,
            new_visit,
            atom_idx,
            RationalI32::new(i32::from(atomic_num), 1),
            new_dist,
            flags,
            mol,
        )
    }

    // BEGIN RDKIT CPP FUNCTION Node::getDigraph (CIPLabeler/Node.cpp)
    // RDKit✔️✔️: Digraph *Node::getDigraph() const { return dp_g; }
    // END RDKIT CPP FUNCTION Node::getDigraph
    pub(crate) fn get_digraph(&self) -> usize {
        self.digraph
    }

    // BEGIN RDKIT CPP FUNCTION Node::getAtom (CIPLabeler/Node.cpp)
    // RDKit✔️✔️: Atom *Node::getAtom() const { return dp_atom; }
    // END RDKIT CPP FUNCTION Node::getAtom
    pub(crate) fn atom_idx(&self) -> Option<usize> {
        self.atom_idx
    }

    // BEGIN RDKIT CPP FUNCTION Node::getAtomIdx (CIPLabeler/Node.cpp)
    // RDKit✔️✔️: unsigned int Node::getAtomIdx() const {
    // RDKit✔️✔️:   if (isSet(IMPL_HYDROGEN)) {
    // RDKit✔️✔️:     return Atom::NOATOM;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return dp_atom->getIdx();
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Node::getAtomIdx
    pub(crate) fn get_atom_idx(&self) -> Result<CipSourceIndex, CipLabelerError> {
        if self.is_set(Self::IMPL_HYDROGEN) {
            Ok(Self::NO_ATOM_INDEX)
        } else {
            let index = self
                .atom_idx
                .expect("RDKit Node::getAtomIdx requires non-null atom");
            CipSourceIndex::try_from(index).map_err(|_| CipLabelerError::SourceIndexWidthExceeded {
                kind: "atom",
                index,
            })
        }
    }

    // BEGIN RDKIT CPP FUNCTION Node::getDistance (CIPLabeler/Node.cpp)
    // RDKit✔️✔️: int Node::getDistance() const { return d_dist; }
    // END RDKIT CPP FUNCTION Node::getDistance
    pub(crate) fn get_distance(&self) -> i32 {
        self.distance
    }

    // BEGIN RDKIT CPP FUNCTION Node::getAtomicNumFraction (CIPLabeler/Node.cpp)
    // RDKit✔️✔️: boost::rational<int> Node::getAtomicNumFraction() const { return d_atomic_num; }
    // END RDKIT CPP FUNCTION Node::getAtomicNumFraction
    pub(crate) fn get_atomic_num_fraction(&self) -> RationalI32 {
        self.atomic_num_fraction
    }

    // BEGIN RDKIT CPP FUNCTION Node::getAtomicNum (CIPLabeler/Node.cpp)
    // RDKit✔️✔️: int Node::getAtomicNum() const {
    // RDKit✔️✔️:   if (dp_atom == nullptr) {
    // RDKit✔️✔️:     return 1;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return dp_atom->getAtomicNum();
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION Node::getAtomicNum
    pub(crate) fn get_atomic_num(&self, mol: &CipMol<'_>) -> Result<u8, CipLabelerError> {
        self.atom_idx
            .map(|idx| mol.atom(idx).map(Atom::atomic_number))
            .transpose()
            .map(|atomic_number| atomic_number.unwrap_or(1))
    }

    // BEGIN RDKIT CPP FUNCTION Node::getMassNum (CIPLabeler/Node.cpp)
    // RDKit✔️✔️: unsigned Node::getMassNum() const {
    // RDKit✔️✔️:   if (dp_atom == nullptr || isDuplicate()) {
    // RDKit✔️✔️:     return 0u;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return dp_atom->getIsotope();
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Node::getMassNum
    pub(crate) fn get_mass_num(&self, mol: &CipMol<'_>) -> Result<u16, CipLabelerError> {
        if self.atom_idx.is_none() || self.is_duplicate() {
            return Ok(0);
        }
        Ok(mol
            .atom(self.atom_idx.expect("checked"))?
            .isotope()
            .unwrap_or(0))
    }

    // BEGIN RDKIT CPP FUNCTION Node::getAtomicMass (CIPLabeler/Node.cpp)
    // RDKit✔️✔️: double Node::getAtomicMass() const { return d_atomic_mass; }
    // END RDKIT CPP FUNCTION Node::getAtomicMass
    pub(crate) fn get_atomic_mass(&self) -> f64 {
        self.atomic_mass
    }

    // BEGIN RDKIT CPP FUNCTION Node::getAux (CIPLabeler/Node.cpp)
    // RDKit✔️✔️: Descriptor Node::getAux() const { return d_aux; }
    // END RDKIT CPP FUNCTION Node::getAux
    pub(crate) fn get_aux(&self) -> Descriptor {
        self.aux
    }

    // BEGIN RDKIT CPP FUNCTION Node::isSet (CIPLabeler/Node.cpp)
    // RDKit✔️✔️: bool Node::isSet(int mask) const { return mask & d_flags; }
    // END RDKIT CPP FUNCTION Node::isSet
    pub(crate) fn is_set(&self, mask: i32) -> bool {
        mask & self.flags != 0
    }

    // BEGIN RDKIT CPP FUNCTION Node::isDuplicate (CIPLabeler/Node.cpp)
    // RDKit✔️✔️: bool Node::isDuplicate() const { return d_flags & DUPLICATE; }
    // END RDKIT CPP FUNCTION Node::isDuplicate
    pub(crate) fn is_duplicate(&self) -> bool {
        self.flags & Self::DUPLICATE != 0
    }

    // BEGIN RDKIT CPP FUNCTION Node::isDuplicateOrH (CIPLabeler/Node.cpp)
    // RDKit✔️✔️: bool Node::isDuplicateOrH() const { return d_flags & DUPLICATE_OR_H; }
    // END RDKIT CPP FUNCTION Node::isDuplicateOrH
    pub(crate) fn is_duplicate_or_h(&self) -> bool {
        self.flags & Self::DUPLICATE_OR_H != 0
    }

    // BEGIN RDKIT CPP FUNCTION Node::isTerminal (CIPLabeler/Node.cpp)
    // RDKit✔️✔️: bool Node::isTerminal() const {
    // RDKit✔️✔️:   return d_visit.empty() || (isExpanded() && d_edges.size() == 1);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Node::isTerminal
    pub(crate) fn is_terminal(&self) -> bool {
        self.visit.is_empty() || (self.is_expanded() && self.edges.len() == 1)
    }

    // BEGIN RDKIT CPP FUNCTION Node::isExpanded (CIPLabeler/Node.cpp)
    // RDKit✔️✔️: bool Node::isExpanded() const { return d_flags & EXPANDED; }
    // END RDKIT CPP FUNCTION Node::isExpanded
    pub(crate) fn is_expanded(&self) -> bool {
        self.flags & Self::EXPANDED != 0
    }

    // BEGIN RDKIT CPP FUNCTION Node::isVisited (CIPLabeler/Node.cpp)
    // RDKit✔️✔️: bool Node::isVisited(int idx) const { return d_visit[idx] != 0; }
    // END RDKIT CPP FUNCTION Node::isVisited
    pub(crate) fn is_visited(&self, idx: usize) -> bool {
        self.visit[idx] != 0
    }

    // BEGIN RDKIT CPP FUNCTION Node::newChild (CIPLabeler/Node.cpp)
    // RDKit✔️✔️: Node *Node::newChild(int idx, Atom *atom) const {
    // RDKit✔️✔️:   auto new_visit = d_visit;
    // RDKit✔️✔️:   new_visit[idx] = static_cast<char>(d_dist + 1);
    // RDKit✔️✔️:   auto atomic_num = atom ? atom->getAtomicNum() : 1;
    // RDKit✔️✔️:   return &dp_g->addNode(std::move(new_visit), atom, atomic_num, d_dist + 1, 0);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Node::newChild
    pub(crate) fn new_child(
        &self,
        idx: usize,
        atom_idx: Option<usize>,
        mol: &CipMol<'_>,
    ) -> Result<Self, CipLabelerError> {
        let mut new_visit = self.visit.clone();
        new_visit[idx] = (self.distance + 1) as i8;
        let atomic_num = atom_idx
            .map(|idx| mol.atom(idx).map(Atom::atomic_number))
            .transpose()?
            .unwrap_or(1);
        Self::new(
            self.digraph,
            new_visit,
            atom_idx,
            RationalI32::new(i32::from(atomic_num), 1),
            self.distance + 1,
            0,
            mol,
        )
    }

    // BEGIN RDKIT CPP FUNCTION Node::newBondDuplicateChild (CIPLabeler/Node.cpp)
    // RDKit✔️✔️: Node *Node::newBondDuplicateChild(int idx, Atom *atom) const {
    // RDKit✔️✔️:   return newTerminalChild(idx, atom, BOND_DUPLICATE);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Node::newBondDuplicateChild
    pub(crate) fn new_bond_duplicate_child(
        &self,
        idx: usize,
        atom_idx: Option<usize>,
        mol: &mut CipMol<'_>,
    ) -> Result<Self, CipLabelerError> {
        self.new_terminal_child(Some(idx), atom_idx, Self::BOND_DUPLICATE, mol)
    }

    // BEGIN RDKIT CPP FUNCTION Node::newRingDuplicateChild (CIPLabeler/Node.cpp)
    // RDKit✔️✔️: Node *Node::newRingDuplicateChild(int idx, Atom *atom) const {
    // RDKit✔️✔️:   return newTerminalChild(idx, atom, RING_DUPLICATE);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Node::newRingDuplicateChild
    pub(crate) fn new_ring_duplicate_child(
        &self,
        idx: usize,
        atom_idx: Option<usize>,
        mol: &mut CipMol<'_>,
    ) -> Result<Self, CipLabelerError> {
        self.new_terminal_child(Some(idx), atom_idx, Self::RING_DUPLICATE, mol)
    }

    // BEGIN RDKIT CPP FUNCTION Node::newImplicitHydrogenChild (CIPLabeler/Node.cpp)
    // RDKit✔️✔️: Node *Node::newImplicitHydrogenChild() const {
    // RDKit✔️✔️:   return newTerminalChild(-1, nullptr, IMPL_HYDROGEN);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Node::newImplicitHydrogenChild
    pub(crate) fn new_implicit_hydrogen_child(
        &self,
        mol: &mut CipMol<'_>,
    ) -> Result<Self, CipLabelerError> {
        self.new_terminal_child(None, None, Self::IMPL_HYDROGEN, mol)
    }

    // BEGIN RDKIT CPP FUNCTION Node::add (CIPLabeler/Node.cpp)
    // RDKit✔️✔️: void Node::add(Edge *e) { d_edges.push_back(e); }
    // END RDKIT CPP FUNCTION Node::add
    pub(crate) fn add(&mut self, edge: CipEdgeId) {
        self.edges.push(edge);
    }

    // BEGIN RDKIT CPP FUNCTION Node::setAux (CIPLabeler/Node.cpp)
    // RDKit✔️✔️: void Node::setAux(Descriptor desc) { d_aux = desc; }
    // END RDKIT CPP FUNCTION Node::setAux
    pub(crate) fn set_aux(&mut self, desc: Descriptor) {
        self.aux = desc;
    }

    // Rust ownership note: the lazy expansion part of RDKit Node::getEdges is
    // source-reproduced in CipDigraph::node_edges(), where the graph can mutate
    // nodes and edges without self-referential borrowing.
    pub(crate) fn get_edges(&self) -> Result<&[CipEdgeId], CipLabelerError> {
        if !self.is_expanded() {
            return Err(CipLabelerError::InvalidInternalState {
                detail: "use CipDigraph::node_edges for Node::getEdges lazy expansion",
            });
        }
        Ok(&self.edges)
    }

    // BEGIN RDKIT CPP FUNCTION Node::getEdges(Atom *) (CIPLabeler/Node.cpp)
    // RDKit✔️✔️: std::vector<Edge *> Node::getEdges(Atom *end) const {
    // RDKit✔️✔️:   std::vector<Edge *> res;
    // RDKit✔️✔️:   for (auto &edge : getEdges()) {
    // RDKit✔️✔️:     if (edge->getEnd()->isDuplicate()) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     };
    // RDKit✔️✔️:     if (end == edge->getBeg()->getAtom() || end == edge->getEnd()->getAtom()) {
    // RDKit✔️✔️:       res.push_back(edge);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Node::getEdges(Atom *)
    pub(crate) fn get_edges_for_atom(
        &self,
        end_atom_idx: Option<usize>,
        nodes: &[CipNode],
        edges: &[CipEdge],
    ) -> Result<Vec<CipEdgeId>, CipLabelerError> {
        let mut result = Vec::new();
        for edge_id in self.get_edges()? {
            let edge = &edges[edge_id.index()];
            if nodes[edge.end.index()].is_duplicate() {
                continue;
            }
            if end_atom_idx == nodes[edge.beg.index()].atom_idx
                || end_atom_idx == nodes[edge.end.index()].atom_idx
            {
                result.push(*edge_id);
            }
        }
        Ok(result)
    }

    // BEGIN RDKIT CPP FUNCTION Node::getNonTerminalOutEdges (CIPLabeler/Node.cpp)
    // RDKit✔️✔️: std::vector<Edge *> Node::getNonTerminalOutEdges() const {
    // RDKit✔️✔️:   std::vector<Edge *> edges;
    // RDKit✔️✔️:   for (auto &edge : getEdges()) {
    // RDKit✔️✔️:     if (edge->isBeg(this) && !edge->getEnd()->isTerminal()) {
    // RDKit✔️✔️:       edges.push_back(edge);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return edges;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Node::getNonTerminalOutEdges
    pub(crate) fn get_non_terminal_out_edges(
        &self,
        self_id: CipNodeId,
        nodes: &[CipNode],
        edges: &[CipEdge],
    ) -> Result<Vec<CipEdgeId>, CipLabelerError> {
        let mut result = Vec::new();
        for edge_id in self.get_edges()? {
            let edge = &edges[edge_id.index()];
            if edge.beg == self_id && !nodes[edge.end.index()].is_terminal() {
                result.push(*edge_id);
            }
        }
        Ok(result)
    }
}

impl<'a> CipDigraph<'a> {
    const MAX_NODE_COUNT: usize = 100_000;
    const MAX_NODE_DIST: i32 = 0;

    // BEGIN RDKIT CPP FUNCTION Digraph::Digraph (CIPLabeler/Digraph.cpp)
    // RDKit✔️✔️: Digraph::Digraph(const CIPMol &mol, Atom *atom, bool atropisomerMode)
    // RDKit✔️✔️:     : d_mol{mol} {
    // RDKit✔️✔️:   PRECONDITION(atom, "cannot init digraph on a nullptr")
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto visit = std::vector<char>(d_mol.getNumAtoms());
    // RDKit✔️✔️:   visit[atom->getIdx()] = 1;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto dist = 1;
    // RDKit✔️✔️:   auto flags = 0x0;
    // RDKit✔️✔️:   auto atomic_num = atom->getAtomicNum();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   dp_root = &addNode(std::move(visit), atom, atomic_num, dist, flags);
    // RDKit✔️✔️:   dp_origin = dp_root;
    // RDKit✔️✔️:   d_atropisomerMode = atropisomerMode;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Digraph::Digraph
    pub(crate) fn new(
        molecule: &'a TopologyBlock,
        atom_idx: usize,
        atropisomer_mode: bool,
    ) -> Result<Self, CipLabelerError> {
        let mol = CipMol::new(molecule);
        mol.atom(atom_idx)?;
        let mut digraph = Self {
            mol,
            origin: CipNodeId::new(0),
            root: CipNodeId::new(0),
            rule6_ref: None,
            atropisomer_mode,
            nodes: Vec::new(),
            edges: Vec::new(),
        };
        let mut visit = vec![0_i8; digraph.mol.get_num_atoms()];
        visit[atom_idx] = 1;
        let atomic_num = digraph.mol.atom(atom_idx)?.atomic_number();
        let root = digraph.add_node(
            visit,
            Some(atom_idx),
            RationalI32::new(i32::from(atomic_num), 1),
            1,
            0,
        )?;
        digraph.root = root;
        digraph.origin = root;
        Ok(digraph)
    }

    // BEGIN RDKIT CPP FUNCTION Digraph::getMol (CIPLabeler/Digraph.cpp)
    // RDKit✔️✔️: const CIPMol &Digraph::getMol() const { return d_mol; };
    // END RDKIT CPP FUNCTION Digraph::getMol
    pub(crate) fn mol(&self) -> &CipMol<'a> {
        &self.mol
    }

    // BEGIN RDKIT CPP FUNCTION Digraph::getOriginalRoot (CIPLabeler/Digraph.cpp)
    // RDKit✔️✔️: Node *Digraph::getOriginalRoot() const { return dp_origin; };
    // END RDKIT CPP FUNCTION Digraph::getOriginalRoot
    pub(crate) fn get_original_root(&self) -> CipNodeId {
        self.origin
    }

    // BEGIN RDKIT CPP FUNCTION Digraph::getCurrentRoot (CIPLabeler/Digraph.cpp)
    // RDKit✔️✔️: Node *Digraph::getCurrentRoot() const { return dp_root; }
    // END RDKIT CPP FUNCTION Digraph::getCurrentRoot
    pub(crate) fn get_current_root(&self) -> CipNodeId {
        self.root
    }

    // BEGIN RDKIT CPP FUNCTION Digraph::getNumNodes (CIPLabeler/Digraph.cpp)
    // RDKit✔️✔️: int Digraph::getNumNodes() const { return d_nodes.size(); }
    // END RDKIT CPP FUNCTION Digraph::getNumNodes
    pub(crate) fn get_num_nodes(&self) -> usize {
        self.nodes.len()
    }

    pub(crate) fn node(&self, node: CipNodeId) -> &CipNode {
        &self.nodes[node.index()]
    }

    pub(crate) fn set_node_aux(&mut self, node: CipNodeId, descriptor: Descriptor) {
        self.nodes[node.index()].set_aux(descriptor);
    }

    pub(crate) fn edge(&self, edge: CipEdgeId) -> &CipEdge {
        &self.edges[edge.index()]
    }

    // BEGIN RDKIT CPP FUNCTION Digraph::addNode (CIPLabeler/Digraph.cpp)
    // RDKit✔️✔️: Node &Digraph::addNode(std::vector<char> &&visit, Atom *atom,
    // RDKit✔️✔️:                        boost::rational<int> &&frac, int dist, int flags) {
    // RDKit✔️✔️:   d_nodes.emplace_back(this, std::move(visit), atom, std::move(frac), dist,
    // RDKit✔️✔️:                        flags);
    // RDKit✔️✔️:   return d_nodes.back();
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Digraph::addNode
    fn add_node(
        &mut self,
        visit: Vec<i8>,
        atom_idx: Option<usize>,
        frac: RationalI32,
        dist: i32,
        flags: i32,
    ) -> Result<CipNodeId, CipLabelerError> {
        let node = CipNode::new(0, visit, atom_idx, frac, dist, flags, &self.mol)?;
        let id = CipNodeId::new(self.nodes.len());
        self.nodes.push(node);
        Ok(id)
    }

    // BEGIN RDKIT CPP FUNCTION Digraph::addEdge (CIPLabeler/Digraph.cpp)
    // RDKit✔️✔️: void Digraph::addEdge(Node *beg, Bond *bond, Node *end) {
    // RDKit✔️✔️:   d_edges.emplace_back(beg, end, bond);
    // RDKit✔️✔️:   auto &e = d_edges.back();
    // RDKit✔️✔️:   beg->add(&e);
    // RDKit✔️✔️:   end->add(&e);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Digraph::addEdge
    fn add_edge(&mut self, beg: CipNodeId, bond_idx: Option<usize>, end: CipNodeId) {
        let edge_id = CipEdgeId::new(self.edges.len());
        self.edges.push(CipEdge::new(beg, end, bond_idx));
        self.nodes[beg.index()].add(edge_id);
        self.nodes[end.index()].add(edge_id);
    }

    // BEGIN RDKIT CPP FUNCTION Node::getEdges (CIPLabeler/Node.cpp)
    // RDKit✔️❌: const std::vector<Edge *> &Node::getEdges() const {
    // RDKit✔️❌:   if (!isExpanded()) {
    // RDKit✔️❌:     auto non_const_this = const_cast<Node *>(this);
    // RDKit✔️❌:     non_const_this->d_flags |= EXPANDED;
    // RDKit✔️❌:     dp_g->expand(non_const_this);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return d_edges;
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION Node::getEdges
    // Returning owned stable IDs avoids exposing references across lazy graph
    // mutation, but clones the incident-edge vector instead of borrowing it.
    pub(crate) fn node_edges(
        &mut self,
        node: CipNodeId,
    ) -> Result<Vec<CipEdgeId>, CipLabelerError> {
        if !self.nodes[node.index()].is_expanded() {
            self.nodes[node.index()].flags |= CipNode::EXPANDED;
            self.expand(node)?;
        }
        Ok(self.nodes[node.index()].edges.clone())
    }

    pub(crate) fn node_edges_for_atom(
        &mut self,
        node: CipNodeId,
        end_atom_idx: Option<usize>,
    ) -> Result<Vec<CipEdgeId>, CipLabelerError> {
        self.node_edges(node)?;
        self.nodes[node.index()].get_edges_for_atom(end_atom_idx, &self.nodes, &self.edges)
    }

    pub(crate) fn non_terminal_out_edges(
        &mut self,
        node: CipNodeId,
    ) -> Result<Vec<CipEdgeId>, CipLabelerError> {
        self.node_edges(node)?;
        self.nodes[node.index()].get_non_terminal_out_edges(node, &self.nodes, &self.edges)
    }

    // BEGIN RDKIT CPP FUNCTION Digraph::getNodes (CIPLabeler/Digraph.cpp)
    // RDKit✔️❌: std::vector<Node *> Digraph::getNodes(Atom *atom) const {
    // RDKit✔️❌:   std::vector<Node *> result;
    // RDKit✔️❌:   std::vector<Node*> queue = {getCurrentRoot()};
    // RDKit✔️❌:
    // RDKit✔️❌:   for (size_t i=0; i<queue.size(); ++i) {
    // RDKit✔️❌:     auto node = queue[i];
    // RDKit✔️❌:     if (atom == node->getAtom()) {
    // RDKit✔️❌:       result.push_back(node);
    // RDKit✔️❌:     }
    // RDKit✔️❌:     for (const auto &e : node->getEdges()) {
    // RDKit✔️❌:       if (!e->isBeg(node)) {
    // RDKit✔️❌:         continue;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       queue.push_back(e->getEnd());
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return result;
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION Digraph::getNodes
    // Each lazy edge access returns an owned ID vector, adding temporary
    // allocations beyond the source's borrowed Node edge view.
    pub(crate) fn get_nodes(&mut self, atom_idx: usize) -> Result<Vec<CipNodeId>, CipLabelerError> {
        self.mol.atom(atom_idx)?;
        let mut result = Vec::new();
        let mut queue = vec![self.get_current_root()];
        let mut i = 0_usize;
        while i < queue.len() {
            let node = queue[i];
            if self.nodes[node.index()].atom_idx == Some(atom_idx) {
                result.push(node);
            }
            for edge_id in self.node_edges(node)? {
                let edge = &self.edges[edge_id.index()];
                if !edge.is_beg(node) {
                    continue;
                }
                queue.push(edge.get_end());
            }
            i += 1;
        }
        Ok(result)
    }

    // BEGIN RDKIT CPP FUNCTION Digraph::getRule6Ref (CIPLabeler/Digraph.cpp)
    // RDKit✔️✔️: Atom *Digraph::getRule6Ref() const { return dp_rule6Ref; }
    // END RDKIT CPP FUNCTION Digraph::getRule6Ref
    pub(crate) fn get_rule6_ref(&self) -> Option<usize> {
        self.rule6_ref
    }

    // BEGIN RDKIT CPP FUNCTION Digraph::setRule6Ref (CIPLabeler/Digraph.cpp)
    // RDKit✔️✔️: void Digraph::setRule6Ref(Atom *ref) { dp_rule6Ref = ref; }
    // END RDKIT CPP FUNCTION Digraph::setRule6Ref
    pub(crate) fn set_rule6_ref(&mut self, atom_idx: Option<usize>) -> Result<(), CipLabelerError> {
        if let Some(atom_idx) = atom_idx {
            self.mol.atom(atom_idx)?;
        }
        self.rule6_ref = atom_idx;
        Ok(())
    }

    // BEGIN RDKIT CPP FUNCTION Digraph::changeRoot (CIPLabeler/Digraph.cpp)
    // RDKit✔️❌: void Digraph::changeRoot(Node *newroot) {
    // RDKit✔️❌:   std::vector<Edge *> toflip;
    // RDKit✔️❌:   auto queue = std::list<Node *>({newroot});
    // RDKit✔️❌:   for (const auto &node : queue) {
    // RDKit✔️❌:     for (const auto &e : node->getEdges()) {
    // RDKit✔️❌:       if (e->isEnd(node)) {
    // RDKit✔️❌:         toflip.push_back(e);
    // RDKit✔️❌:         queue.push_back(e->getBeg());
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   for (auto &e : toflip) {
    // RDKit✔️❌:     e->flip();
    // RDKit✔️❌:   }
    // RDKit✔️❌:   dp_root = newroot;
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION Digraph::changeRoot
    // The Rust traversal retains source complexity and order but additionally
    // clones each lazily expanded node's edge-ID vector.
    pub(crate) fn change_root(&mut self, new_root: CipNodeId) -> Result<(), CipLabelerError> {
        let mut to_flip = Vec::new();
        let mut queue = vec![new_root];
        let mut i = 0_usize;
        while i < queue.len() {
            let node = queue[i];
            for edge_id in self.node_edges(node)? {
                let edge = &self.edges[edge_id.index()];
                if edge.is_end(node) {
                    to_flip.push(edge_id);
                    queue.push(edge.get_beg());
                }
            }
            i += 1;
        }
        for edge_id in to_flip {
            self.edges[edge_id.index()].flip();
        }
        self.root = new_root;
        Ok(())
    }

    // BEGIN RDKIT CPP FUNCTION Digraph::expand (CIPLabeler/Digraph.cpp)
    // RDKit✔️❌: void Digraph::expand(Node *beg) {
    // RDKit✔️❌:   const auto &atom = beg->getAtom();
    // RDKit✔️❌:   const auto &edges = beg->getEdges();
    // RDKit✔️❌:   const auto &prev =
    // RDKit✔️❌:       edges.size() > 0 && !edges[0]->isBeg(beg) ? edges[0]->getBond() : nullptr;
    // RDKit✔️❌:
    // RDKit✔️❌:   if (MAX_NODE_DIST > 0 && beg->getDistance() > MAX_NODE_DIST) {
    // RDKit✔️❌:     return;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   if (MAX_NODE_COUNT > 0 && d_nodes.size() >= MAX_NODE_COUNT) {
    // RDKit✔️❌:     std::stringstream errmsg;
    // RDKit✔️❌:     errmsg << "Digraph generation failed: more than " << MAX_NODE_COUNT
    // RDKit✔️❌:            << "nodes found.";
    // RDKit✔️❌:     throw TooManyNodesException(errmsg.str());
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // create 'explicit' nodes
    // RDKit✔️❌:   for (const auto &bond : d_mol.getBonds(atom)) {
    // RDKit✔️❌:     const auto &nbr = bond->getOtherAtom(atom);
    // RDKit✔️❌:     const int nbrIdx = nbr->getIdx();
    // RDKit✔️❌:     const int bord = d_mol.getBondOrder(bond);
    // RDKit✔️❌:     const int virtual_nodes = bord - 1;
    // RDKit✔️❌:
    // RDKit✔️❌:     if (!beg->isVisited(nbrIdx)) {
    // RDKit✔️❌:       auto end = beg->newChild(nbrIdx, nbr);
    // RDKit✔️❌:       addEdge(beg, bond, end);
    // RDKit✔️❌:
    // RDKit✔️❌:       // duplicate nodes for bond orders (except for root atoms...)
    // RDKit✔️❌:       // for example >S=O
    // RDKit✔️❌:       if (dp_origin != beg || d_atropisomerMode) {
    // RDKit✔️❌:         if (atom->getFormalCharge() < 0 &&
    // RDKit✔️❌:             d_mol.getFractionalAtomicNum(atom).denominator() > 1) {
    // RDKit✔️❌:           end = beg->newBondDuplicateChild(nbrIdx, nbr);
    // RDKit✔️❌:           addEdge(beg, bond, end);
    // RDKit✔️❌:         } else {
    // RDKit✔️❌:           for (int i = 0; i < virtual_nodes; ++i) {
    // RDKit✔️❌:             end = beg->newBondDuplicateChild(nbrIdx, nbr);
    // RDKit✔️❌:             addEdge(beg, bond, end);
    // RDKit✔️❌:           }
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:     } else if (bond == prev) {  // bond order expansion (backwards)
    // RDKit✔️❌:       if (dp_origin->getAtom() != nbr || d_atropisomerMode) {
    // RDKit✔️❌:         for (int i = 0; i < virtual_nodes; ++i) {
    // RDKit✔️❌:           auto end = beg->newBondDuplicateChild(nbrIdx, nbr);
    // RDKit✔️❌:           addEdge(beg, bond, end);
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:     } else {  // ring closures
    // RDKit✔️❌:       auto end = beg->newRingDuplicateChild(nbrIdx, nbr);
    // RDKit✔️❌:       addEdge(beg, bond, end);
    // RDKit✔️❌:
    // RDKit✔️❌:       if (atom->getFormalCharge() < 0 &&
    // RDKit✔️❌:           d_mol.getFractionalAtomicNum(atom).denominator() > 1) {
    // RDKit✔️❌:         end = beg->newBondDuplicateChild(nbrIdx, nbr);
    // RDKit✔️❌:         addEdge(beg, bond, end);
    // RDKit✔️❌:       } else {
    // RDKit✔️❌:         for (int i = 0; i < virtual_nodes; ++i) {
    // RDKit✔️❌:           end = beg->newBondDuplicateChild(nbrIdx, nbr);
    // RDKit✔️❌:           addEdge(beg, bond, end);
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // Create implicit hydrogen nodes
    // RDKit✔️❌:   const int hcnt = atom->getTotalNumHs();
    // RDKit✔️❌:   for (int i = 0; i < hcnt; ++i) {
    // RDKit✔️❌:     auto end = beg->newImplicitHydrogenChild();
    // RDKit✔️❌:     addEdge(beg, nullptr, end);
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION Digraph::expand
    // `bond_indices_for_atom()` materializes the adjacency span so expansion
    // can mutably grow the graph; this is a known allocation difference.
    fn expand(&mut self, beg: CipNodeId) -> Result<(), CipLabelerError> {
        let atom_idx =
            self.nodes[beg.index()]
                .atom_idx
                .ok_or(CipLabelerError::InvalidInternalState {
                    detail: "Digraph::expand on null atom node",
                })?;
        let prev = self.nodes[beg.index()].edges.first().and_then(|edge_id| {
            let edge = &self.edges[edge_id.index()];
            (!edge.is_beg(beg)).then_some(edge.get_bond_idx()).flatten()
        });

        if Self::MAX_NODE_DIST > 0 && self.nodes[beg.index()].get_distance() > Self::MAX_NODE_DIST {
            return Ok(());
        }
        if Self::MAX_NODE_COUNT > 0 && self.nodes.len() >= Self::MAX_NODE_COUNT {
            return Err(CipLabelerError::TooManyNodes {
                limit: Self::MAX_NODE_COUNT,
            });
        }

        for bond_idx in self.mol.bond_indices_for_atom(atom_idx)? {
            let nbr_idx = self.mol.other_atom_idx(bond_idx, atom_idx)?;
            let bord = self.mol.get_bond_order(bond_idx)?;
            let virtual_nodes = bord - 1;

            if !self.nodes[beg.index()].is_visited(nbr_idx) {
                let end_node =
                    self.nodes[beg.index()].new_child(nbr_idx, Some(nbr_idx), &self.mol)?;
                let mut end = self.add_existing_node(end_node);
                self.add_edge(beg, Some(bond_idx), end);

                if self.origin != beg || self.atropisomer_mode {
                    let atom_formal_charge = self.mol.atom(atom_idx)?.formal_charge();
                    if atom_formal_charge < 0
                        && self.mol.get_fractional_atomic_num(atom_idx)?.denominator > 1
                    {
                        let end_node = self.nodes[beg.index()].new_bond_duplicate_child(
                            nbr_idx,
                            Some(nbr_idx),
                            &mut self.mol,
                        )?;
                        end = self.add_existing_node(end_node);
                        self.add_edge(beg, Some(bond_idx), end);
                    } else {
                        for _ in 0..virtual_nodes {
                            let end_node = self.nodes[beg.index()].new_bond_duplicate_child(
                                nbr_idx,
                                Some(nbr_idx),
                                &mut self.mol,
                            )?;
                            end = self.add_existing_node(end_node);
                            self.add_edge(beg, Some(bond_idx), end);
                        }
                    }
                }
            } else if Some(bond_idx) == prev {
                if self.nodes[self.origin.index()].atom_idx != Some(nbr_idx)
                    || self.atropisomer_mode
                {
                    for _ in 0..virtual_nodes {
                        let end_node = self.nodes[beg.index()].new_bond_duplicate_child(
                            nbr_idx,
                            Some(nbr_idx),
                            &mut self.mol,
                        )?;
                        let end = self.add_existing_node(end_node);
                        self.add_edge(beg, Some(bond_idx), end);
                    }
                }
            } else {
                let end_node = self.nodes[beg.index()].new_ring_duplicate_child(
                    nbr_idx,
                    Some(nbr_idx),
                    &mut self.mol,
                )?;
                let mut end = self.add_existing_node(end_node);
                self.add_edge(beg, Some(bond_idx), end);

                let atom_formal_charge = self.mol.atom(atom_idx)?.formal_charge();
                if atom_formal_charge < 0
                    && self.mol.get_fractional_atomic_num(atom_idx)?.denominator > 1
                {
                    let end_node = self.nodes[beg.index()].new_bond_duplicate_child(
                        nbr_idx,
                        Some(nbr_idx),
                        &mut self.mol,
                    )?;
                    end = self.add_existing_node(end_node);
                    self.add_edge(beg, Some(bond_idx), end);
                } else {
                    for _ in 0..virtual_nodes {
                        let end_node = self.nodes[beg.index()].new_bond_duplicate_child(
                            nbr_idx,
                            Some(nbr_idx),
                            &mut self.mol,
                        )?;
                        end = self.add_existing_node(end_node);
                        self.add_edge(beg, Some(bond_idx), end);
                    }
                }
            }
        }

        let hcnt = self.mol.total_num_hs(atom_idx)?;
        for _ in 0..hcnt {
            let end_node = self.nodes[beg.index()].new_implicit_hydrogen_child(&mut self.mol)?;
            let end = self.add_existing_node(end_node);
            self.add_edge(beg, None, end);
        }
        Ok(())
    }

    fn add_existing_node(&mut self, node: CipNode) -> CipNodeId {
        let id = CipNodeId::new(self.nodes.len());
        self.nodes.push(node);
        id
    }
}

impl<'a> CipMol<'a> {
    // BEGIN RDKIT CPP FUNCTION CIPMol::CIPMol (CIPMol.cpp)
    // RDKit✔️✔️: CIPMol::CIPMol(ROMol &mol) : d_mol{mol} {}
    // END RDKIT CPP FUNCTION CIPMol::CIPMol
    pub(crate) fn new(molecule: &'a TopologyBlock) -> Self {
        Self {
            molecule,
            rings: None,
            kekulized_bond_orders: None,
            fractional_atomic_numbers: None,
            valence: None,
        }
    }

    // BEGIN RDKIT CPP FUNCTION CIPMol::getFractionalAtomicNum (CIPMol.cpp)
    // RDKit✔️✔️: boost::rational<int> CIPMol::getFractionalAtomicNum(Atom *atom) const {
    // RDKit✔️✔️:   PRECONDITION(atom, "bad atom")
    // RDKit✔️✔️:   if (d_atomnums.empty()) {
    // RDKit✔️✔️:     const_cast<CIPMol *>(this)->d_atomnums = calcFracAtomNums(*this);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return d_atomnums[atom->getIdx()];
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION CIPMol::getFractionalAtomicNum
    pub(crate) fn get_fractional_atomic_num(
        &mut self,
        atom_idx: usize,
    ) -> Result<RationalI32, CipLabelerError> {
        self.atom(atom_idx)?;
        if self.fractional_atomic_numbers.is_none() {
            self.fractional_atomic_numbers = Some(calc_frac_atom_nums(self)?);
        }
        Ok(self
            .fractional_atomic_numbers
            .as_ref()
            .expect("initialized")[atom_idx])
    }

    // BEGIN RDKIT CPP FUNCTION CIPMol::getNumAtoms (CIPMol.cpp)
    // RDKit✔️✔️: unsigned CIPMol::getNumAtoms() const { return d_mol.getNumAtoms(); }
    // END RDKIT CPP FUNCTION CIPMol::getNumAtoms
    pub(crate) fn get_num_atoms(&self) -> usize {
        self.molecule.num_atoms()
    }

    // BEGIN RDKIT CPP FUNCTION CIPMol::getNumBonds (CIPMol.cpp)
    // RDKit✔️✔️: unsigned CIPMol::getNumBonds() const { return d_mol.getNumBonds(); };
    // END RDKIT CPP FUNCTION CIPMol::getNumBonds
    pub(crate) fn get_num_bonds(&self) -> usize {
        self.molecule.num_bonds()
    }

    // BEGIN RDKIT CPP FUNCTION CIPMol::getAtom (CIPMol.cpp)
    // RDKit✔️✔️: Atom *CIPMol::getAtom(int idx) const { return d_mol.getAtomWithIdx(idx); };
    // END RDKIT CPP FUNCTION CIPMol::getAtom
    pub(crate) fn atom(&self, atom_idx: usize) -> Result<&'a Atom, CipLabelerError> {
        self.molecule
            .atoms()
            .get(atom_idx)
            .ok_or(CipLabelerError::AtomIndexOutOfRange {
                index: atom_idx,
                atom_count: self.molecule.num_atoms(),
            })
    }

    // BEGIN RDKIT CPP FUNCTION CIPMol::atoms (CIPLabeler/CIPMol.cpp)
    // RDKit✔️✔️: CXXAtomIterator<MolGraph, Atom *> CIPMol::atoms() const {
    // RDKit✔️✔️:   return d_mol.atoms();
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION CIPMol::atoms
    pub(crate) fn atoms(&self) -> &'a [Atom] {
        self.molecule.atoms()
    }

    // BEGIN RDKIT CPP FUNCTION CIPMol::getBond (CIPMol.cpp)
    // RDKit✔️✔️: Bond *CIPMol::getBond(int idx) const { return d_mol.getBondWithIdx(idx); };
    // END RDKIT CPP FUNCTION CIPMol::getBond
    pub(crate) fn bond(&self, bond_idx: usize) -> Result<&'a Bond, CipLabelerError> {
        self.molecule
            .bonds()
            .get(bond_idx)
            .ok_or(CipLabelerError::BondIndexOutOfRange {
                index: bond_idx,
                bond_count: self.molecule.num_bonds(),
            })
    }

    // BEGIN RDKIT CPP FUNCTION CIPMol::getBonds (CIPMol.cpp)
    // RDKit✔️❌: CIPMolSpan<Bond *, ROMol::OEDGE_ITER> CIPMol::getBonds(Atom *atom) const {
    // RDKit✔️❌:   PRECONDITION(atom, "bad atom")
    // RDKit✔️❌:   return {d_mol, d_mol.getAtomBonds(atom)};
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION CIPMol::getBonds
    // The source returns a lazy span; this detached adapter materializes stable
    // bond IDs and therefore is behaviorally exact but allocation-heavier.
    pub(crate) fn bond_indices_for_atom(
        &self,
        atom_idx: usize,
    ) -> Result<Vec<usize>, CipLabelerError> {
        self.atom(atom_idx)?;
        Ok(self
            .molecule
            .topology_block()
            .adjacency
            .neighbors_of(atom_idx)
            .iter()
            .map(|neighbor| neighbor.bond.index())
            .collect())
    }

    // BEGIN RDKIT CPP FUNCTION CIPMol::getNeighbors (CIPMol.cpp)
    // RDKit✔️❌: CIPMolSpan<Atom *, ROMol::ADJ_ITER> CIPMol::getNeighbors(Atom *atom) const {
    // RDKit✔️❌:   PRECONDITION(atom, "bad atom")
    // RDKit✔️❌:   return {d_mol, d_mol.getAtomNeighbors(atom)};
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION CIPMol::getNeighbors
    // The source returns a lazy span; this detached adapter materializes stable
    // atom IDs and therefore is behaviorally exact but allocation-heavier.
    pub(crate) fn neighbor_indices(&self, atom_idx: usize) -> Result<Vec<usize>, CipLabelerError> {
        self.atom(atom_idx)?;
        Ok(self
            .molecule
            .topology_block()
            .adjacency
            .neighbors_of(atom_idx)
            .iter()
            .map(|neighbor| neighbor.atom_index)
            .collect())
    }

    // BEGIN RDKIT CPP FUNCTION CIPMol::isInRing (CIPMol.cpp)
    // RDKit✔️✔️: bool CIPMol::isInRing(Bond *bond) const {
    // RDKit✔️✔️:   PRECONDITION(bond, "bad bond")
    // RDKit✔️✔️:   const auto rings = d_mol.getRingInfo();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (!rings->isFindFastOrBetter()) {
    // RDKit✔️✔️:     MolOps::fastFindRings(d_mol);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return rings->numBondRings(bond->getIdx()) != 0u;
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION CIPMol::isInRing
    pub(crate) fn is_in_ring(&mut self, bond_idx: usize) -> Result<bool, CipLabelerError> {
        self.bond(bond_idx)?;
        if self.rings.is_none() {
            self.rings = Some(fast_find_rings(self.molecule)?);
        }
        Ok(self
            .rings
            .as_ref()
            .expect("initialized")
            .num_bond_rings(BondId::new(bond_idx))
            != 0)
    }

    // BEGIN RDKIT CPP FUNCTION CIPMol::getBondOrder (CIPMol.cpp)
    // RDKit✔️✔️: int CIPMol::getBondOrder(Bond *bond) const {
    // RDKit✔️✔️:   PRECONDITION(bond, "bad bond")
    // RDKit✔️✔️:   if (d_kekulized_bonds.empty()) {
    // RDKit✔️✔️:     RWMol tmp{d_mol};
    // RDKit✔️✔️:     try {
    // RDKit✔️✔️:       MolOps::Kekulize(tmp);
    // RDKit✔️✔️:     } catch (const MolSanitizeException &) {
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     auto& bonds = const_cast<std::vector<RDKit::Bond::BondType>&>(d_kekulized_bonds);
    // RDKit✔️✔️:     bonds.reserve(d_mol.getNumBonds());
    // RDKit✔️✔️:     for (const auto &b : tmp.bonds()) {
    // RDKit✔️✔️:       bonds.push_back(b->getBondType());
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   const auto bond_type = d_kekulized_bonds.at(bond->getIdx());
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // Dative bonds might need to be considered with a different bond order
    // RDKit✔️✔️:   // for the end atom at the end of the bond.
    // RDKit✔️✔️:   switch (bond_type) {
    // RDKit✔️✔️:     case Bond::ZERO:
    // RDKit✔️✔️:     case Bond::HYDROGEN:
    // RDKit✔️✔️:     case Bond::DATIVE:
    // RDKit✔️✔️:     case Bond::DATIVEL:
    // RDKit✔️✔️:     case Bond::DATIVER:
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     case Bond::SINGLE:
    // RDKit✔️✔️:       return 1;
    // RDKit✔️✔️:     case Bond::AROMATIC:
    // RDKit✔️✔️:       BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:           << "non kekulizable aromatic bond being treated as bond order 1"
    // RDKit✔️✔️:           << std::endl;
    // RDKit✔️✔️:       return 1;
    // RDKit✔️✔️:     case Bond::DOUBLE:
    // RDKit✔️✔️:       return 2;
    // RDKit✔️✔️:     case Bond::TRIPLE:
    // RDKit✔️✔️:       return 3;
    // RDKit✔️✔️:     case Bond::QUADRUPLE:
    // RDKit✔️✔️:       return 4;
    // RDKit✔️✔️:     case Bond::QUINTUPLE:
    // RDKit✔️✔️:       return 5;
    // RDKit✔️✔️:     case Bond::HEXTUPLE:
    // RDKit✔️✔️:       return 6;
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       throw std::runtime_error("Non integer-order bonds are not allowed.");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION CIPMol::getBondOrder
    pub(crate) fn get_bond_order(&mut self, bond_idx: usize) -> Result<i32, CipLabelerError> {
        self.bond(bond_idx)?;
        if self.kekulized_bond_orders.is_none() {
            let mut orders = self
                .molecule
                .bonds()
                .iter()
                .map(Bond::order)
                .collect::<Vec<_>>();
            match kekulize_if_possible(self.molecule, &KekulizeParams::default())? {
                KekulizeAttempt::Applied(assignment) => {
                    for (order, bond) in orders.iter_mut().zip(&assignment.topology.bonds) {
                        *order = bond.order();
                    }
                }
                KekulizeAttempt::NotKekulizable { .. } => {}
            }
            self.kekulized_bond_orders = Some(orders);
        }

        match self.kekulized_bond_orders.as_ref().expect("initialized")[bond_idx] {
            BondOrder::Zero
            | BondOrder::Hydrogen
            | BondOrder::Dative
            | BondOrder::DativeLeft
            | BondOrder::DativeRight => Ok(0),
            BondOrder::Single => Ok(1),
            BondOrder::Aromatic => Ok(1),
            BondOrder::Double => Ok(2),
            BondOrder::Triple => Ok(3),
            BondOrder::Quadruple => Ok(4),
            BondOrder::Quintuple => Ok(5),
            BondOrder::Hextuple => Ok(6),
            order => Err(CipLabelerError::NonIntegerBondOrder { order }),
        }
    }

    fn total_num_hs(&mut self, atom_idx: usize) -> Result<i32, CipLabelerError> {
        let explicit = i32::from(self.atom(atom_idx)?.explicit_hydrogens());
        if self.valence.is_none() {
            self.valence = Some(assign_valence_for_topology(
                self.molecule,
                ValenceModel::RdkitLike,
            )?);
        }
        let implicit = self
            .valence
            .as_ref()
            .and_then(|valence| valence.implicit_hydrogens.get(atom_idx))
            .copied()
            .unwrap_or(0)
            .max(0);
        Ok(explicit + implicit)
    }

    fn other_atom_idx(&self, bond_idx: usize, atom_idx: usize) -> Result<usize, CipLabelerError> {
        let bond = self.bond(bond_idx)?;
        if bond.begin().index() == atom_idx {
            Ok(bond.end().index())
        } else if bond.end().index() == atom_idx {
            Ok(bond.begin().index())
        } else {
            Err(CipLabelerError::BondNotIncident {
                bond: bond_idx,
                atom: atom_idx,
            })
        }
    }
}

// BEGIN RDKIT CPP ENUM Type (CIPLabeler/Mancude.h)
// RDKit✔️✔️: enum class Type {
// RDKit✔️✔️:   Cv4D3,       // =C(X)-
// RDKit✔️✔️:   Nv3D2,       // =N-
// RDKit✔️✔️:   Nv4D3Plus,   // =[N+]<
// RDKit✔️✔️:   Nv2D2Minus,  // -[N-]-
// RDKit✔️✔️:   Cv3D3Minus,  // -[C(X)-]-
// RDKit✔️✔️:   Ov3D2Plus,   // -[O+]=
// RDKit✔️✔️:   Other
// RDKit✔️✔️: };
// END RDKIT CPP ENUM Type

fn seed_types(types: &mut [MancudeType], mol: &mut CipMol<'_>) -> Result<bool, CipLabelerError> {
    // BEGIN RDKIT CPP FUNCTION SeedTypes (Mancude.cpp)
    // RDKit✔️✔️: bool SeedTypes(std::vector<Type> &types, const CIPMol &mol) {
    // RDKit✔️✔️:   bool result = false;
    // RDKit✔️✔️:   for (const auto &atom : mol.atoms()) {
    // RDKit✔️✔️:     const int aidx = atom->getIdx();
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // check ring
    // RDKit✔️✔️:     int btypes = atom->getTotalNumHs();
    // RDKit✔️✔️:     bool ring = false;
    // RDKit✔️✔️:     for (const auto &bond : mol.getBonds(atom)) {
    // RDKit✔️✔️:       // Given the possible types we have, we only care
    // RDKit✔️✔️:       // for single and double bonds which are in rings.
    // RDKit✔️✔️:       switch (mol.getBondOrder(bond)) {
    // RDKit✔️✔️:         case 1:
    // RDKit✔️✔️:           btypes += 0x00000001;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case 2:
    // RDKit✔️✔️:           btypes += 0x00000100;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         default:
    // RDKit✔️✔️:           btypes += 0x01000000;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (mol.isInRing(bond)) {
    // RDKit✔️✔️:         ring = true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (ring) {
    // RDKit✔️✔️:       int q = atom->getFormalCharge();
    // RDKit✔️✔️:       switch (atom->getAtomicNum()) {
    // RDKit✔️✔️:         case 6:   // C
    // RDKit✔️✔️:         case 14:  // Si
    // RDKit✔️✔️:         case 32:  // Ge
    // RDKit✔️✔️:           if (q == 0 && btypes == 0x0102) {
    // RDKit✔️✔️:             types[aidx] = Type::Cv4D3;
    // RDKit✔️✔️:           } else if (q == -1 && btypes == 0x0003) {
    // RDKit✔️✔️:             types[aidx] = Type::Cv3D3Minus;
    // RDKit✔️✔️:             result = true;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case 7:   // N
    // RDKit✔️✔️:         case 15:  // P
    // RDKit✔️✔️:         case 33:  // As
    // RDKit✔️✔️:           if (q == 0 && btypes == 0x0101) {
    // RDKit✔️✔️:             types[aidx] = Type::Nv3D2;
    // RDKit✔️✔️:             result = true;
    // RDKit✔️✔️:           } else if (q == -1 && btypes == 0x0002) {
    // RDKit✔️✔️:             types[aidx] = Type::Nv2D2Minus;
    // RDKit✔️✔️:             result = true;
    // RDKit✔️✔️:           } else if (q == +1 && btypes == 0x0102) {
    // RDKit✔️✔️:             types[aidx] = Type::Nv4D3Plus;
    // RDKit✔️✔️:             result = true;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case 8:  // O
    // RDKit✔️✔️:           if (q == 1 && btypes == 0x0101) {
    // RDKit✔️✔️:             types[aidx] = Type::Ov3D2Plus;
    // RDKit✔️✔️:             result = true;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return result;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION SeedTypes
    let mut result = false;
    for atom_idx in 0..mol.get_num_atoms() {
        let mut btypes = mol.total_num_hs(atom_idx)?;
        let mut ring = false;
        for bond_idx in mol.bond_indices_for_atom(atom_idx)? {
            match mol.get_bond_order(bond_idx)? {
                1 => btypes += 0x00000001,
                2 => btypes += 0x00000100,
                _ => btypes += 0x01000000,
            }
            if mol.is_in_ring(bond_idx)? {
                ring = true;
            }
        }
        if !ring {
            continue;
        }
        let atom = mol.atom(atom_idx)?;
        let q = atom.formal_charge();
        match atom.atomic_number() {
            6 | 14 | 32 => {
                if q == 0 && btypes == 0x0102 {
                    types[atom_idx] = MancudeType::Cv4D3;
                } else if q == -1 && btypes == 0x0003 {
                    types[atom_idx] = MancudeType::Cv3D3Minus;
                    result = true;
                }
            }
            7 | 15 | 33 => {
                if q == 0 && btypes == 0x0101 {
                    types[atom_idx] = MancudeType::Nv3D2;
                    result = true;
                } else if q == -1 && btypes == 0x0002 {
                    types[atom_idx] = MancudeType::Nv2D2Minus;
                    result = true;
                } else if q == 1 && btypes == 0x0102 {
                    types[atom_idx] = MancudeType::Nv4D3Plus;
                    result = true;
                }
            }
            8 => {
                if q == 1 && btypes == 0x0101 {
                    types[atom_idx] = MancudeType::Ov3D2Plus;
                    result = true;
                }
            }
            _ => {}
        }
    }
    Ok(result)
}

fn relax_types(types: &mut [MancudeType], mol: &CipMol<'_>) -> Result<(), CipLabelerError> {
    // BEGIN RDKIT CPP FUNCTION RelaxTypes (Mancude.cpp)
    // RDKit✔️✔️: void RelaxTypes(std::vector<Type> &types, const CIPMol &mol) {
    // RDKit✔️✔️:   std::list<Atom *> queue;
    // RDKit✔️✔️:   auto counts = std::vector<int>(mol.getNumAtoms());
    // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:     const auto aidx = atom->getIdx();
    // RDKit✔️✔️:     for (const auto &nbr : mol.getNeighbors(atom)) {
    // RDKit✔️✔️:       if (types[nbr->getIdx()] != Type::Other) {
    // RDKit✔️✔️:         ++counts[aidx];
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (counts[aidx] == 1) {
    // RDKit✔️✔️:       queue.push_back(atom);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (const auto &atom : queue) {
    // RDKit✔️✔️:     const auto aidx = atom->getIdx();
    // RDKit✔️✔️:     if (types[aidx] != Type::Other) {
    // RDKit✔️✔️:       types[aidx] = Type::Other;
    // RDKit✔️✔️:
    // RDKit✔️✔️:       for (auto &nbr : mol.getNeighbors(atom)) {
    // RDKit✔️✔️:         auto nbridx = nbr->getIdx();
    // RDKit✔️✔️:         --counts[nbridx];
    // RDKit✔️✔️:         if (counts[nbridx] == 1) {
    // RDKit✔️✔️:           queue.push_back(nbr);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RelaxTypes
    let mut queue = VecDeque::new();
    let mut counts = vec![0_i32; mol.get_num_atoms()];
    for atom_idx in 0..mol.get_num_atoms() {
        for nbr_idx in mol.neighbor_indices(atom_idx)? {
            if types[nbr_idx] != MancudeType::Other {
                counts[atom_idx] += 1;
            }
        }
        if counts[atom_idx] == 1 {
            queue.push_back(atom_idx);
        }
    }

    while let Some(atom_idx) = queue.pop_front() {
        if types[atom_idx] == MancudeType::Other {
            continue;
        }
        types[atom_idx] = MancudeType::Other;
        for nbr_idx in mol.neighbor_indices(atom_idx)? {
            counts[nbr_idx] -= 1;
            if counts[nbr_idx] == 1 {
                queue.push_back(nbr_idx);
            }
        }
    }
    Ok(())
}

fn visit_part(
    parts: &mut [i32],
    types: &[MancudeType],
    part: i32,
    mut atom_idx: usize,
    mol: &mut CipMol<'_>,
) -> Result<(), CipLabelerError> {
    // BEGIN RDKIT CPP FUNCTION VisitPart (Mancude.cpp)
    // RDKit✔️✔️: void VisitPart(std::vector<int> &parts, const std::vector<Type> &types,
    // RDKit✔️✔️:                int part, Atom *atom, const CIPMol &mol) {
    // RDKit✔️✔️:   Atom *next;
    // RDKit✔️✔️:   do {
    // RDKit✔️✔️:     next = nullptr;
    // RDKit✔️✔️:     for (auto &bond : mol.getBonds(atom)) {
    // RDKit✔️✔️:       if (!mol.isInRing(bond)) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       auto nbr = bond->getOtherAtom(atom);
    // RDKit✔️✔️:       int aidx = nbr->getIdx();
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (parts[aidx] == 0 && types[aidx] != Type::Other) {
    // RDKit✔️✔️:         parts[aidx] = part;
    // RDKit✔️✔️:         if (next != nullptr) {
    // RDKit✔️✔️:           VisitPart(parts, types, part, nbr, mol);
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           next = nbr;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     atom = next;
    // RDKit✔️✔️:   } while (atom != nullptr);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION VisitPart
    loop {
        let mut next = None;
        for bond_idx in mol.bond_indices_for_atom(atom_idx)? {
            if !mol.is_in_ring(bond_idx)? {
                continue;
            }
            let nbr_idx = mol.other_atom_idx(bond_idx, atom_idx)?;
            if parts[nbr_idx] == 0 && types[nbr_idx] != MancudeType::Other {
                parts[nbr_idx] = part;
                if next.is_some() {
                    visit_part(parts, types, part, nbr_idx, mol)?;
                } else {
                    next = Some(nbr_idx);
                }
            }
        }
        if let Some(next_idx) = next {
            atom_idx = next_idx;
        } else {
            break;
        }
    }
    Ok(())
}

fn visit_parts(
    parts: &mut [i32],
    types: &[MancudeType],
    mol: &mut CipMol<'_>,
) -> Result<i32, CipLabelerError> {
    // BEGIN RDKIT CPP FUNCTION VisitParts (Mancude.cpp)
    // RDKit✔️✔️: int VisitParts(std::vector<int> &parts, const std::vector<Type> &types,
    // RDKit✔️✔️:                const CIPMol &mol) {
    // RDKit✔️✔️:   int numparts = 0;
    // RDKit✔️✔️:   for (auto &atom : mol.atoms()) {
    // RDKit✔️✔️:     int aidx = atom->getIdx();
    // RDKit✔️✔️:     if (parts[aidx] == 0 && types[aidx] != Type::Other) {
    // RDKit✔️✔️:       parts[aidx] = ++numparts;
    // RDKit✔️✔️:       VisitPart(parts, types, parts[aidx], atom, mol);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return numparts;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION VisitParts
    let mut numparts = 0_i32;
    for atom_idx in 0..mol.get_num_atoms() {
        if parts[atom_idx] == 0 && types[atom_idx] != MancudeType::Other {
            numparts += 1;
            parts[atom_idx] = numparts;
            visit_part(parts, types, numparts, atom_idx, mol)?;
        }
    }
    Ok(numparts)
}

fn calc_frac_atom_nums(mol: &mut CipMol<'_>) -> Result<Vec<RationalI32>, CipLabelerError> {
    // BEGIN RDKIT CPP FUNCTION calcFracAtomNums (Mancude.cpp)
    // RDKit✔️✔️: std::vector<boost::rational<int>> calcFracAtomNums(const CIPMol &mol) {
    // RDKit✔️✔️:   const auto num_atoms = mol.getNumAtoms();
    // RDKit✔️✔️:   std::vector<boost::rational<int>> fractions;
    // RDKit✔️✔️:   fractions.reserve(num_atoms);
    // RDKit✔️✔️:   for (const auto &atom : mol.atoms()) {
    // RDKit✔️✔️:     fractions.emplace_back(atom->getAtomicNum(), 1);
    // RDKit✔️✔️:   }
    let num_atoms = mol.get_num_atoms();
    let mut fractions = Vec::with_capacity(num_atoms);
    for atom_idx in 0..num_atoms {
        fractions.push(RationalI32::new(
            i32::from(mol.atom(atom_idx)?.atomic_number()),
            1,
        ));
    }

    // RDKit✔️✔️:   // Mark all atoms which are potentially part of a resonance system.
    // RDKit✔️✔️:   auto types = std::vector<Type>(num_atoms, Type::Other);
    // RDKit✔️✔️:   if (SeedTypes(types, mol)) {
    let mut types = vec![MancudeType::Other; num_atoms];
    if seed_types(&mut types, mol)? {
        // RDKit✔️✔️:     // Filter out atoms which cannot be resonant because
        // RDKit✔️✔️:     // of not having the proper environment.
        // RDKit✔️✔️:     RelaxTypes(types, mol);
        relax_types(&mut types, mol)?;

        // RDKit✔️✔️:     // Find resonant systems: parts stores the ids of the
        // RDKit✔️✔️:     // systems each atom is involved in.
        // RDKit✔️✔️:     auto parts = std::vector<int>(num_atoms);
        // RDKit✔️✔️:     int numparts = VisitParts(parts, types, mol);
        let mut parts = vec![0_i32; num_atoms];
        let numparts = visit_parts(&mut parts, &types, mol)?;

        // RDKit✔️✔️:     auto resparts = std::vector<int>(numparts);
        // RDKit✔️✔️:     int numres = 0;
        let mut resparts = vec![0_i32; usize::try_from(numparts).unwrap_or(0)];
        let mut numres = 0_usize;

        // RDKit✔️✔️:     if (numparts > 0) {
        // RDKit✔️✔️:       for (auto i = 0u; i < num_atoms; ++i) {
        // RDKit✔️✔️:         if (parts[i] == 0) {
        // RDKit✔️✔️:           continue;
        // RDKit✔️✔️:         }
        if numparts > 0 {
            for i in 0..num_atoms {
                if parts[i] == 0 {
                    continue;
                }

                // RDKit✔️✔️:         // Find resonant structures caused by relocation of a negative charge.
                // RDKit✔️✔️:         if (types[i] == Type::Cv3D3Minus || types[i] == Type::Nv2D2Minus) {
                // RDKit✔️✔️:           int j = 0;
                // RDKit✔️✔️:           for (; j < numres; ++j) {
                // RDKit✔️✔️:             if (resparts[j] == parts[i]) {
                // RDKit✔️✔️:               break;
                // RDKit✔️✔️:             }
                // RDKit✔️✔️:           }
                // RDKit✔️✔️:           if (j >= numres) {
                // RDKit✔️✔️:             resparts[numres] = parts[i];
                // RDKit✔️✔️:             ++numres;
                // RDKit✔️✔️:           }
                // RDKit✔️✔️:         }
                if matches!(types[i], MancudeType::Cv3D3Minus | MancudeType::Nv2D2Minus) {
                    let mut j = 0_usize;
                    while j < numres {
                        if resparts[j] == parts[i] {
                            break;
                        }
                        j += 1;
                    }
                    if j >= numres {
                        resparts[numres] = parts[i];
                        numres += 1;
                    }
                }

                // RDKit✔️✔️:         int numerator = 0;
                // RDKit✔️✔️:         int denominator = 0;
                // RDKit✔️✔️:         for (const auto &nbr : mol.getNeighbors(atom)) {
                // RDKit✔️✔️:           if (parts[nbr->getIdx()] == parts[i]) {
                // RDKit✔️✔️:             numerator += nbr->getAtomicNum();
                // RDKit✔️✔️:             ++denominator;
                // RDKit✔️✔️:           }
                // RDKit✔️✔️:         }
                let mut numerator = 0_i32;
                let mut denominator = 0_i32;
                for nbr_idx in mol.neighbor_indices(i)? {
                    if parts[nbr_idx] == parts[i] {
                        numerator += i32::from(mol.atom(nbr_idx)?.atomic_number());
                        denominator += 1;
                    }
                }

                // RDKit✔️✔️:         // boost::rational does not accept 0 as denominator.
                // RDKit✔️✔️:         if (denominator == 0) {
                // RDKit✔️✔️:           fractions[i].assign(0, 1);
                // RDKit✔️✔️:         } else {
                // RDKit✔️✔️:           fractions[i].assign(numerator, denominator);
                // RDKit✔️✔️:         }
                if denominator == 0 {
                    fractions[i].assign(0, 1);
                } else {
                    fractions[i].assign(numerator, denominator);
                }
            }
        }

        // RDKit✔️✔️:     // If there are any resonant structures due to negative charges,
        // RDKit✔️✔️:     // recalculate the average atomic number considering relocation
        // RDKit✔️✔️:     // of the charge through higher order bonds.
        // RDKit✔️✔️:     if (numres > 0) {
        // RDKit✔️✔️:       for (int j = 0; j < numres; ++j) {
        // RDKit✔️✔️:         int numerator = 0;
        // RDKit✔️✔️:         int denominator = 0;
        // RDKit✔️✔️:         int part = resparts[j];
        // RDKit✔️✔️:         for (auto i = 0u; i < num_atoms; ++i) {
        // RDKit✔️✔️:           if (parts[i] == part) {
        if numres > 0 {
            for &part in resparts.iter().take(numres) {
                let mut numerator = 0_i32;
                let mut denominator = 0_i32;
                for i in 0..num_atoms {
                    if parts[i] != part {
                        continue;
                    }
                    // RDKit✔️✔️:             // boost::rational does not accept 0 as denominator
                    // RDKit✔️✔️:             if (denominator == 0) {
                    // RDKit✔️✔️:               fractions[i].assign(0, 1);
                    // RDKit✔️✔️:             } else {
                    // RDKit✔️✔️:               fractions[i].assign(numerator, denominator);
                    // RDKit✔️✔️:             }
                    if denominator == 0 {
                        fractions[i].assign(0, 1);
                    } else {
                        fractions[i].assign(numerator, denominator);
                    }

                    // RDKit✔️✔️:             ++denominator;
                    // RDKit✔️✔️:             auto atom = mol.getAtom(i);
                    // RDKit✔️✔️:             for (auto &bond : mol.getBonds(atom)) {
                    // RDKit✔️✔️:               auto nbr = bond->getOtherAtom(atom);
                    // RDKit✔️✔️:               int bord = mol.getBondOrder(bond);
                    // RDKit✔️✔️:               if (bord > 1 && parts[nbr->getIdx()] == part) {
                    // RDKit✔️✔️:                 numerator += (bord - 1) * nbr->getAtomicNum();
                    // RDKit✔️✔️:               }
                    // RDKit✔️✔️:             }
                    denominator += 1;
                    for bond_idx in mol.bond_indices_for_atom(i)? {
                        let nbr_idx = mol.other_atom_idx(bond_idx, i)?;
                        let bord = mol.get_bond_order(bond_idx)?;
                        if bord > 1 && parts[nbr_idx] == part {
                            numerator += (bord - 1) * i32::from(mol.atom(nbr_idx)?.atomic_number());
                        }
                    }
                }
            }
        }
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return fractions;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION calcFracAtomNums
    Ok(fractions)
}
