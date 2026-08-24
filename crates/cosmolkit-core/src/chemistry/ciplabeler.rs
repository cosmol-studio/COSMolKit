//! Single source-backed implementation of modern RDKit CIP labeling.
//!
//! Public molecule operations and fingerprint consumers delegate to this
//! crate-private `CIPLabeler::assignCIPLabels` core.

use std::collections::VecDeque;

use crate::chemistry::atropisomer::atropisomer_atoms_and_bonds;
use crate::chemistry::valence::rdkit_atomic_mass;
use crate::{
    Atom, AtomId, Bond, BondId, BondOrder, BondStereo, ChiralTag, Molecule, RingInfo,
    ValenceAssignment, ValenceModel,
};

pub(crate) use crate::CipLabelerError;

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
    if let Ok(public_descriptor) = crate::CipDescriptor::try_from(desc) {
        return public_descriptor.as_str();
    }
    match desc {
        Descriptor::None => "NONE",
        Descriptor::Unknown => "UNKNOWN",
        Descriptor::ns => "ns",
        Descriptor::seqTrans => "e",
        Descriptor::seqCis => "z",
        Descriptor::SP_4 => "SP_4",
        Descriptor::TBPY_5 => "TBPY_5",
        Descriptor::OC_6 => "OC_6",
        Descriptor::R
        | Descriptor::S
        | Descriptor::r
        | Descriptor::s
        | Descriptor::E
        | Descriptor::Z
        | Descriptor::M
        | Descriptor::P
        | Descriptor::m
        | Descriptor::p => unreachable!("public emitted descriptors returned above"),
    }
}

impl TryFrom<Descriptor> for crate::CipDescriptor {
    type Error = ();

    fn try_from(descriptor: Descriptor) -> Result<Self, Self::Error> {
        match descriptor {
            Descriptor::R => Ok(Self::R),
            Descriptor::S => Ok(Self::S),
            Descriptor::r => Ok(Self::LowerR),
            Descriptor::s => Ok(Self::LowerS),
            Descriptor::E => Ok(Self::E),
            Descriptor::Z => Ok(Self::Z),
            Descriptor::M => Ok(Self::M),
            Descriptor::P => Ok(Self::P),
            Descriptor::m => Ok(Self::LowerM),
            Descriptor::p => Ok(Self::LowerP),
            _ => Err(()),
        }
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
    molecule: &'a Molecule,
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CipCancellationMode {
    NotRequested,
    RdkitControlC,
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

    fn with_remaining_call_count(remaining_call_count: u32) -> Self {
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
    // RDKit✔️✔️: Rules() = delete;
    // RDKit✔️✔️:
    // RDKit✔️✔️: Rules(std::initializer_list<SequenceRule *> rules) {
    // RDKit✔️✔️:   for (auto &rule : rules) {
    // RDKit✔️✔️:     add(rule);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️:
    // RDKit✔️✔️: ~Rules() override {
    // RDKit✔️✔️:   for (auto &rule : d_rules) {
    // RDKit✔️✔️:     delete rule;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️:
    // RDKit✔️✔️: void add(SequenceRule *rule) {
    // RDKit✔️✔️:   if (rule == nullptr) {
    // RDKit✔️✔️:     throw std::runtime_error("No sequence rule provided");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   d_rules.push_back(rule);
    // RDKit✔️✔️:   rule->setSorter(new Sort(d_rules));
    // RDKit✔️✔️: }
    // RDKit✔️✔️:
    // RDKit✔️✔️: int getNumSubRules() const { return d_rules.size(); }
    // RDKit✔️✔️:
    // RDKit✔️✔️: const Sort *getSorter() const override {
    // RDKit✔️✔️:   if (dp_sorter == nullptr) {
    // RDKit✔️✔️:     const_cast<Rules *>(this)->setSorter(new Sort(this));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return dp_sorter.get();
    // RDKit✔️✔️: }
    // END RDKIT CPP CLASS Rules constructor/add/destructor/getNumSubRules/getSorter
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

fn cip_all_rules() -> Result<CipRules, CipLabelerError> {
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

fn cip_constitutional_rules() -> Result<CipRules, CipLabelerError> {
    // BEGIN RDKIT CPP CONSTANT constitutional_rules (CIPLabeler.cpp)
    // RDKit✔️✔️: const Rules constitutional_rules({new Rule1a, new Rule1b, new Rule2});
    // END RDKIT CPP CONSTANT constitutional_rules
    CipRules::new(vec![
        Box::new(CipRule1a),
        Box::new(CipRule1b),
        Box::new(CipRule2),
    ])
}

// BEGIN RDKIT CPP FUNCTION findConfigs (CIPLabeler.cpp)
// RDKit✔️✔️: std::vector<std::unique_ptr<Configuration>> findConfigs(
// RDKit✔️✔️:     CIPMol &mol, const boost::dynamic_bitset<> &atoms,
// RDKit✔️✔️:     const boost::dynamic_bitset<> &bonds) {
// RDKit✔️✔️:   std::vector<std::unique_ptr<Configuration>> configs;
// RDKit✔️✔️:
// RDKit✔️✔️:   for (auto index = atoms.find_first(); index != boost::dynamic_bitset<>::npos;
// RDKit✔️✔️:        index = atoms.find_next(index)) {
// RDKit✔️✔️:     auto atom = mol.getAtom(index);
// RDKit✔️✔️:     auto chiraltag = atom->getChiralTag();
// RDKit✔️✔️:     if (chiraltag == Atom::CHI_TETRAHEDRAL_CW ||
// RDKit✔️✔️:         chiraltag == Atom::CHI_TETRAHEDRAL_CCW) {
// RDKit✔️✔️:       std::unique_ptr<Tetrahedral> cfg{new Tetrahedral(mol, atom)};
// RDKit✔️✔️:       configs.push_back(std::move(cfg));
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   for (auto index = bonds.find_first(); index != boost::dynamic_bitset<>::npos;
// RDKit✔️✔️:        index = bonds.find_next(index)) {
// RDKit✔️✔️:     auto bond = mol.getBond(index);
// RDKit✔️✔️:
// RDKit✔️✔️:     auto bond_cfg = bond->getStereo();
// RDKit✔️✔️:     switch (bond_cfg) {
// RDKit✔️✔️:       case Bond::STEREOE:
// RDKit✔️✔️:         bond_cfg = Bond::STEREOTRANS;
// RDKit✔️✔️:         break;
// RDKit✔️✔️:       case Bond::STEREOZ:
// RDKit✔️✔️:         bond_cfg = Bond::STEREOCIS;
// RDKit✔️✔️:         break;
// RDKit✔️✔️:       default:
// RDKit✔️✔️:         break;
// RDKit✔️✔️:     }
// RDKit✔️✔️:     switch (bond_cfg) {
// RDKit✔️✔️:       case Bond::STEREOTRANS:
// RDKit✔️✔️:       case Bond::STEREOCIS: {
// RDKit✔️✔️:         std::unique_ptr<Sp2Bond> cfg(new Sp2Bond(
// RDKit✔️✔️:             mol, bond, bond->getBeginAtom(), bond->getEndAtom(), bond_cfg));
// RDKit✔️✔️:         configs.push_back(std::move(cfg));
// RDKit✔️✔️:       } break;
// RDKit✔️✔️:
// RDKit✔️✔️:       case Bond::STEREOATROPCCW:
// RDKit✔️✔️:       case Bond::STEREOATROPCW: {
// RDKit✔️✔️:         std::unique_ptr<AtropisomerBond> cfgAtrop(new AtropisomerBond(
// RDKit✔️✔️:             mol, bond, bond->getBeginAtom(), bond->getEndAtom(), bond_cfg));
// RDKit✔️✔️:         configs.push_back(std::move(cfgAtrop));
// RDKit✔️✔️:       } break;
// RDKit✔️✔️:
// RDKit✔️✔️:       default:
// RDKit✔️✔️:         break;
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   return configs;
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION findConfigs
fn cip_find_configs<'a>(
    molecule: &'a Molecule,
    atom_mask: &[bool],
    bond_mask: &[bool],
) -> Result<Vec<CipConfig<'a>>, CipLabelerError> {
    if atom_mask.len() != molecule.num_atoms() {
        return Err(CipLabelerError::AtomMaskLengthMismatch {
            actual: atom_mask.len(),
            expected: molecule.num_atoms(),
        });
    }
    if bond_mask.len() != molecule.num_bonds() {
        return Err(CipLabelerError::BondMaskLengthMismatch {
            actual: bond_mask.len(),
            expected: molecule.num_bonds(),
        });
    }
    let mut configs = Vec::new();

    for (index, selected) in atom_mask.iter().copied().enumerate() {
        if !selected {
            continue;
        }
        let atom = molecule
            .atoms()
            .get(index)
            .ok_or(CipLabelerError::AtomIndexOutOfRange {
                index,
                atom_count: molecule.num_atoms(),
            })?;
        if matches!(
            atom.chiral_tag(),
            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
        ) {
            configs.push(CipConfig::Tetrahedral(CipTetrahedral::new(
                molecule, index,
            )?));
        }
    }

    for (index, selected) in bond_mask.iter().copied().enumerate() {
        if !selected {
            continue;
        }
        let bond = molecule
            .bonds()
            .get(index)
            .ok_or(CipLabelerError::BondIndexOutOfRange {
                index,
                bond_count: molecule.num_bonds(),
            })?;
        let bond_cfg = match bond.stereo() {
            BondStereo::E => BondStereo::Trans,
            BondStereo::Z => BondStereo::Cis,
            other => other,
        };
        match bond_cfg {
            BondStereo::Trans | BondStereo::Cis => {
                configs.push(CipConfig::Sp2Bond(CipSp2Bond::new(
                    molecule,
                    index,
                    bond.begin().index(),
                    bond.end().index(),
                    bond_cfg,
                )?));
            }
            BondStereo::AtropCcw | BondStereo::AtropCw => {
                configs.push(CipConfig::AtropisomerBond(CipAtropisomerBond::new(
                    molecule,
                    index,
                    bond.begin().index(),
                    bond.end().index(),
                    bond_cfg,
                )?));
            }
            _ => {}
        }
    }

    Ok(configs)
}

fn cip_label_with_center_digraph(
    configs: &mut [CipConfig<'_>],
    center_idx: usize,
    target_idx: usize,
    node: CipNodeId,
    rules: &CipRules,
    context: &mut CipLabelerContext,
) -> Result<Descriptor, CipLabelerError> {
    debug_assert_ne!(center_idx, target_idx);
    if center_idx < target_idx {
        let (left, right) = configs.split_at_mut(target_idx);
        let center = &mut left[center_idx];
        let target = &mut right[0];
        let digraph = center.get_digraph_mut();
        target.label_with_external_digraph(node, digraph, rules, context)
    } else {
        let (left, right) = configs.split_at_mut(center_idx);
        let target = &mut left[target_idx];
        let center = &mut right[0];
        let digraph = center.get_digraph_mut();
        target.label_with_external_digraph(node, digraph, rules, context)
    }
}

fn cip_set_center_node_aux(
    configs: &mut [CipConfig<'_>],
    center_idx: usize,
    node: CipNodeId,
    desc: Descriptor,
) {
    let center = &mut configs[center_idx];
    center.get_digraph_mut().nodes[node.index()].set_aux(desc);
}

// BEGIN RDKIT CPP FUNCTION labelAux (CIPLabeler.cpp)
// RDKit✔️✔️: bool labelAux(std::vector<std::unique_ptr<Configuration>> &configs,
// RDKit✔️✔️:               const Rules &rules,
// RDKit✔️✔️:               const std::unique_ptr<Configuration> &center) {
// RDKit✔️✔️:   using Node_Cfg_Pair = std::pair<Node *, Configuration *>;
// RDKit✔️✔️:   std::vector<Node_Cfg_Pair> aux;
// RDKit✔️✔️:
// RDKit✔️✔️:   auto &digraph = center->getDigraph();
// RDKit✔️✔️:   for (const auto &config : configs) {
// RDKit✔️✔️:     if (config == center) {
// RDKit✔️✔️:       continue;
// RDKit✔️✔️:     }
// RDKit✔️✔️:     // FIXME: specific to each descriptor
// RDKit✔️✔️:     const auto &foci = config->getFoci();
// RDKit✔️✔️:     for (const auto &node : digraph.getNodes(foci[0])) {
// RDKit✔️✔️:       if (node->isDuplicate()) {
// RDKit✔️✔️:         continue;
// RDKit✔️✔️:       }
// RDKit✔️✔️:       auto low = node;
// RDKit✔️✔️:       if (foci.size() == 2) {
// RDKit✔️✔️:         for (const auto &edge : node->getEdges(foci[1])) {
// RDKit✔️✔️:           const auto &other_node = edge->getOther(node);
// RDKit✔️✔️:           if (other_node->getDistance() < node->getDistance()) {
// RDKit✔️✔️:             low = other_node;
// RDKit✔️✔️:           }
// RDKit✔️✔️:         }
// RDKit✔️✔️:       }
// RDKit✔️✔️:       if (!low->isDuplicate()) {
// RDKit✔️✔️:         aux.emplace_back(low, config.get());
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   auto farthest = [](const Node_Cfg_Pair &a, const Node_Cfg_Pair &b) {
// RDKit✔️✔️:     return a.first->getDistance() > b.first->getDistance();
// RDKit✔️✔️:   };
// RDKit✔️✔️:   std::sort(aux.begin(), aux.end(), farthest);
// RDKit✔️✔️:
// RDKit✔️✔️:   // Using a boost::unordered_map because it is more performant
// RDKit✔️✔️:   // than the STL version.
// RDKit✔️✔️:   boost::unordered_map<Node *, Descriptor> queue;
// RDKit✔️✔️:   int prev = std::numeric_limits<int>::max();
// RDKit✔️✔️:   for (const auto &e : aux) {
// RDKit✔️✔️:     const auto &node = e.first;
// RDKit✔️✔️:
// RDKit✔️✔️:     if (node->getDistance() < prev) {
// RDKit✔️✔️:       for (const auto &e2 : queue) {
// RDKit✔️✔️:         e2.first->setAux(e2.second);
// RDKit✔️✔️:       }
// RDKit✔️✔️:       queue.clear();
// RDKit✔️✔️:       prev = node->getDistance();
// RDKit✔️✔️:     }
// RDKit✔️✔️:     const auto &config = e.second;
// RDKit✔️✔️:     auto label = config->label(node, digraph, rules);
// RDKit✔️✔️:     queue.emplace(node, label);
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   for (const auto &e : queue) {
// RDKit✔️✔️:     e.first->setAux(e.second);
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   return true;
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION labelAux
fn cip_label_aux(
    configs: &mut [CipConfig<'_>],
    rules: &CipRules,
    center_idx: usize,
    context: &mut CipLabelerContext,
) -> Result<bool, CipLabelerError> {
    let config_foci = configs
        .iter()
        .enumerate()
        .filter(|(idx, _)| *idx != center_idx)
        .map(|(idx, config)| {
            let foci = config.get_foci();
            (idx, foci[0], foci.get(1).copied())
        })
        .collect::<Vec<_>>();
    let mut aux = Vec::<(CipNodeId, usize, i32)>::new();

    {
        let digraph = configs[center_idx].get_digraph_mut();
        for (config_idx, first_focus, second_focus) in config_foci {
            for node in digraph.get_nodes(first_focus)? {
                if digraph.node(node).is_duplicate() {
                    continue;
                }
                let mut low = node;
                if let Some(second_focus) = second_focus {
                    for edge in digraph.node_edges_for_atom(node, Some(second_focus))? {
                        let other_node = digraph.edge(edge).get_other(edge, node)?;
                        if digraph.node(other_node).get_distance()
                            < digraph.node(node).get_distance()
                        {
                            low = other_node;
                        }
                    }
                }
                if !digraph.node(low).is_duplicate() {
                    aux.push((low, config_idx, digraph.node(low).get_distance()));
                }
            }
        }
    }

    aux.sort_by(|left, right| right.2.cmp(&left.2));

    let mut queue = Vec::<(CipNodeId, Descriptor)>::new();
    let mut prev = i32::MAX;
    for (node, config_idx, distance) in aux {
        if distance < prev {
            for (queued_node, desc) in queue.drain(..) {
                cip_set_center_node_aux(configs, center_idx, queued_node, desc);
            }
            prev = distance;
        }
        let label =
            cip_label_with_center_digraph(configs, center_idx, config_idx, node, rules, context)?;
        if !queue.iter().any(|(queued_node, _)| *queued_node == node) {
            queue.push((node, label));
        }
    }

    for (queued_node, desc) in queue {
        cip_set_center_node_aux(configs, center_idx, queued_node, desc);
    }

    Ok(true)
}

// BEGIN RDKIT CPP FUNCTION label (CIPLabeler.cpp)
// RDKit✔️✔️: void label(std::vector<std::unique_ptr<Configuration>> &configs,
// RDKit✔️✔️:            unsigned int maxRecursiveIterations) {
// RDKit✔️✔️:   // First, if the specified number of iterations allows it, run all centers
// RDKit✔️✔️:   // through a fast pass with the constitutional rules allow easy stuff to be
// RDKit✔️✔️:   // resolved.
// RDKit✔️✔️:   for (auto &conf : configs) {
// RDKit✔️✔️:     // Make sure this stereo center has no label
// RDKit✔️✔️:     conf->resetPrimaryLabel();
// RDKit✔️✔️:
// RDKit✔️✔️:     remainingCallCount = constitutionalRuleTimeout;
// RDKit✔️✔️:     try {
// RDKit✔️✔️:       auto desc = conf->label(constitutional_rules);
// RDKit✔️✔️:       if (desc != Descriptor::UNKNOWN) {
// RDKit✔️✔️:         conf->setPrimaryLabel(desc);
// RDKit✔️✔️:       }
// RDKit✔️✔️:     } catch (const MaxIterationsExceeded &) {
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   // Now, retry everything that hasn't been solved with a more generous
// RDKit✔️✔️:   // threshold
// RDKit✔️✔️:   if (maxRecursiveIterations != 0) {
// RDKit✔️✔️:     remainingCallCount = maxRecursiveIterations;
// RDKit✔️✔️:   } else {
// RDKit✔️✔️:     remainingCallCount = UINT_MAX;  // really big - will never be hit
// RDKit✔️✔️:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   // try again on everything that hasn't been resolved yet
// RDKit✔️✔️:   for (const auto &conf : configs) {
// RDKit✔️✔️:     if (conf->hasPrimaryLabel()) {
// RDKit✔️✔️:       // already resolved!
// RDKit✔️✔️:       continue;
// RDKit✔️✔️:     }
// RDKit✔️✔️:
// RDKit✔️✔️:     auto desc = conf->label(constitutional_rules);
// RDKit✔️✔️:     if (desc != Descriptor::UNKNOWN) {
// RDKit✔️✔️:       conf->setPrimaryLabel(desc);
// RDKit✔️✔️:     } else {
// RDKit✔️✔️:       if (labelAux(configs, all_rules, conf)) {
// RDKit✔️✔️:         desc = conf->label(all_rules);
// RDKit✔️✔️:
// RDKit✔️✔️:         if (desc != Descriptor::UNKNOWN) {
// RDKit✔️✔️:           conf->setPrimaryLabel(desc);
// RDKit✔️✔️:         }
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION label
fn cip_label(
    configs: &mut [CipConfig<'_>],
    max_recursive_iterations: u32,
) -> Result<(), CipLabelerError> {
    let constitutional_rules = cip_constitutional_rules()?;
    for conf in configs.iter_mut() {
        conf.reset_primary_label();
        let mut context = CipLabelerContext::with_remaining_call_count(
            CipLabelerContext::CONSTITUTIONAL_RULE_TIMEOUT,
        );
        match conf.label(&constitutional_rules, &mut context) {
            Ok(desc) if desc != Descriptor::Unknown => conf.set_primary_label(desc)?,
            Ok(_) | Err(CipLabelerError::MaxIterationsExceeded) => {}
            Err(err) => return Err(err),
        }
    }

    let constitutional_rules = cip_constitutional_rules()?;
    let all_rules = cip_all_rules()?;
    let mut context = CipLabelerContext::new(max_recursive_iterations);
    for idx in 0..configs.len() {
        if configs[idx].has_primary_label() {
            continue;
        }

        let desc = configs[idx].label(&constitutional_rules, &mut context)?;
        if desc != Descriptor::Unknown {
            configs[idx].set_primary_label(desc)?;
        } else if cip_label_aux(configs, &all_rules, idx, &mut context)? {
            let desc = configs[idx].label(&all_rules, &mut context)?;
            if desc != Descriptor::Unknown {
                configs[idx].set_primary_label(desc)?;
            }
        }
    }
    Ok(())
}

fn cip_neighbor_order_value(indices: &[CipSourceIndex]) -> String {
    let body = indices
        .iter()
        .map(CipSourceIndex::to_string)
        .collect::<Vec<_>>()
        .join(",");
    format!("[{body}]")
}

fn cip_clear_selected_primary_codes(
    molecule: &mut Molecule,
    atom_mask: &[bool],
    bond_mask: &[bool],
) {
    for (idx, selected) in atom_mask.iter().copied().enumerate() {
        if !selected {
            continue;
        }
        if let Some(atom) = molecule.topology_block_mut().atoms.get_mut(idx) {
            if matches!(
                atom.chiral_tag(),
                ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
            ) {
                atom.clear_prop("_CIPCode");
            }
        }
    }
    for (idx, selected) in bond_mask.iter().copied().enumerate() {
        if !selected {
            continue;
        }
        if let Some(bond) = molecule.topology_block_mut().bonds.get_mut(idx) {
            if matches!(
                bond.stereo(),
                BondStereo::E
                    | BondStereo::Z
                    | BondStereo::Cis
                    | BondStereo::Trans
                    | BondStereo::AtropCcw
                    | BondStereo::AtropCw
            ) {
                bond.clear_prop("_CIPCode");
            }
        }
    }
}

fn cip_apply_primary_labels(
    molecule: &mut Molecule,
    labels: Vec<CipPrimaryLabel>,
) -> Result<(), CipLabelerError> {
    let topology = molecule.topology_block_mut();
    for label in labels {
        match label {
            CipPrimaryLabel::Atom(label) => {
                let atom_count = topology.atoms.len();
                let atom = topology.atoms.get_mut(label.atom_idx).ok_or(
                    CipLabelerError::AtomIndexOutOfRange {
                        index: label.atom_idx,
                        atom_count,
                    },
                )?;
                atom.set_prop("_CIPCode", label.cip_code);
                atom.set_computed_prop(
                    "_CIPNeighborOrder",
                    cip_neighbor_order_value(&label.cip_neighbor_order),
                );
            }
            CipPrimaryLabel::Bond(label) => {
                let bond_count = topology.bonds.len();
                let bond = topology.bonds.get_mut(label.bond_idx).ok_or(
                    CipLabelerError::BondIndexOutOfRange {
                        index: label.bond_idx,
                        bond_count,
                    },
                )?;
                bond.set_stereo_atoms(Some([
                    AtomId::new(label.stereo_atoms[0]),
                    AtomId::new(label.stereo_atoms[1]),
                ]));
                bond.set_stereo(label.stereo);
                bond.set_prop("_CIPCode", label.cip_code);
                bond.set_computed_prop(
                    "_CIPNeighborOrder",
                    cip_neighbor_order_value(&label.cip_neighbor_order),
                );
            }
            CipPrimaryLabel::AtropisomerBond(label) => {
                let bond_count = topology.bonds.len();
                let bond = topology.bonds.get_mut(label.bond_idx).ok_or(
                    CipLabelerError::BondIndexOutOfRange {
                        index: label.bond_idx,
                        bond_count,
                    },
                )?;
                bond.set_prop("_CIPCode", label.cip_code);
                bond.set_computed_prop(
                    "_CIPNeighborOrder",
                    cip_neighbor_order_value(&label.cip_neighbor_order),
                );
            }
        }
    }
    Ok(())
}

struct CipAssignmentOutcome {
    molecule: Molecule,
    error: Option<CipLabelerError>,
}

impl CipAssignmentOutcome {
    fn into_result(self) -> Result<Molecule, CipLabelerError> {
        match self.error {
            Some(error) => Err(error),
            None => Ok(self.molecule),
        }
    }
}

// BEGIN RDKIT CPP FUNCTION assignCIPLabels selected overload (CIPLabeler.cpp)
// RDKit✔️✔️: void assignCIPLabels(ROMol &mol, const boost::dynamic_bitset<> &atoms,
// RDKit✔️✔️:                      const boost::dynamic_bitset<> &bonds,
// RDKit✔️✔️:                      unsigned int maxRecursiveIterations) {
// RDKit✔️✔️:   ControlCHandler::reset();
// RDKit✔️✔️:
// RDKit✔️✔️:   // reset the mark, for the case that this fails
// RDKit✔️✔️:   mol.clearProp(common_properties::_CIPComputed);
// RDKit✔️✔️:   CIPMol cipmol{mol};
// RDKit✔️✔️:   auto configs = findConfigs(cipmol, atoms, bonds);
// RDKit✔️✔️:
// RDKit✔️✔️:   try {
// RDKit✔️✔️:     label(configs, maxRecursiveIterations);
// RDKit❌❌:   } catch (const ControlCCaught &) {
// RDKit❌❌:   }
// RDKit❌❌:   if (ControlCHandler::getGotSignal()) {
// RDKit❌❌:     BOOST_LOG(rdWarningLog)
// RDKit❌❌:         << "Interrupted, cancelling CIP label calculation" << std::endl;
// RDKit❌❌:     return;
// RDKit❌❌:   }
// RDKit✔️✔️:
// RDKit✔️✔️:   const bool computed = true;
// RDKit✔️✔️:   mol.setProp(common_properties::_CIPComputed, true, computed);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION assignCIPLabels selected overload
fn assign_cip_labels_for_masks_with_cancellation(
    molecule: &Molecule,
    atom_mask: &[bool],
    bond_mask: &[bool],
    max_recursive_iterations: u32,
    cancellation: CipCancellationMode,
) -> Result<Molecule, CipLabelerError> {
    assign_cip_labels_for_masks_outcome_with_cancellation(
        molecule,
        atom_mask,
        bond_mask,
        max_recursive_iterations,
        cancellation,
    )?
    .into_result()
}

fn assign_cip_labels_for_masks_outcome_with_cancellation(
    molecule: &Molecule,
    atom_mask: &[bool],
    bond_mask: &[bool],
    max_recursive_iterations: u32,
    cancellation: CipCancellationMode,
) -> Result<CipAssignmentOutcome, CipLabelerError> {
    if cancellation == CipCancellationMode::RdkitControlC {
        return Err(CipLabelerError::CancellationUnsupported);
    }

    let mut labeled = molecule.clone();
    labeled.properties_mut().clear_prop("_CIPComputed");

    let mut configs = match cip_find_configs(&labeled, atom_mask, bond_mask) {
        Ok(configs) => configs,
        Err(error) => {
            return Ok(CipAssignmentOutcome {
                molecule: labeled,
                error: Some(error),
            });
        }
    };
    let label_result = cip_label(&mut configs, max_recursive_iterations);
    let labels = configs
        .iter()
        .filter_map(CipConfig::primary_label)
        .collect::<Vec<_>>();
    drop(configs);

    cip_clear_selected_primary_codes(&mut labeled, atom_mask, bond_mask);
    if let Err(error) = cip_apply_primary_labels(&mut labeled, labels) {
        return Ok(CipAssignmentOutcome {
            molecule: labeled,
            error: Some(error),
        });
    }
    if let Err(error) = label_result {
        return Ok(CipAssignmentOutcome {
            molecule: labeled,
            error: Some(error),
        });
    }
    labeled
        .properties_mut()
        .set_computed_prop("_CIPComputed", "1");
    Ok(CipAssignmentOutcome {
        molecule: labeled,
        error: None,
    })
}

fn assign_cip_labels_for_masks(
    molecule: &Molecule,
    atom_mask: &[bool],
    bond_mask: &[bool],
    max_recursive_iterations: u32,
) -> Result<Molecule, CipLabelerError> {
    assign_cip_labels_for_masks_with_cancellation(
        molecule,
        atom_mask,
        bond_mask,
        max_recursive_iterations,
        CipCancellationMode::NotRequested,
    )
}

fn cip_exact_selection_masks(
    molecule: &Molecule,
    atom_indices: &[CipSourceIndex],
    bond_indices: &[CipSourceIndex],
) -> Result<(Vec<bool>, Vec<bool>), CipLabelerError> {
    let mut atom_mask = vec![false; molecule.num_atoms()];
    for &index in atom_indices {
        let index = index as usize;
        let selected = atom_mask
            .get_mut(index)
            .ok_or(CipLabelerError::AtomIndexOutOfRange {
                index,
                atom_count: molecule.num_atoms(),
            })?;
        *selected = true;
    }

    let mut bond_mask = vec![false; molecule.num_bonds()];
    for &index in bond_indices {
        let index = index as usize;
        let selected = bond_mask
            .get_mut(index)
            .ok_or(CipLabelerError::BondIndexOutOfRange {
                index,
                bond_count: molecule.num_bonds(),
            })?;
        *selected = true;
    }
    Ok((atom_mask, bond_mask))
}

pub(crate) fn assign_cip_labels_for_indices(
    molecule: &Molecule,
    atom_indices: &[CipSourceIndex],
    bond_indices: &[CipSourceIndex],
    max_recursive_iterations: u32,
) -> Result<Molecule, CipLabelerError> {
    let (atom_mask, bond_mask) = cip_exact_selection_masks(molecule, atom_indices, bond_indices)?;
    assign_cip_labels_for_masks(molecule, &atom_mask, &bond_mask, max_recursive_iterations)
}

// BEGIN RDKIT CPP FUNCTION assignCIPLabels all-molecule overload (CIPLabeler.cpp)
// RDKit✔️✔️: void assignCIPLabels(ROMol &mol, unsigned int maxRecursiveIterations) {
// RDKit✔️✔️:   boost::dynamic_bitset<> atoms(mol.getNumAtoms());
// RDKit✔️✔️:   boost::dynamic_bitset<> bonds(mol.getNumBonds());
// RDKit✔️✔️:   atoms.set();
// RDKit✔️✔️:   bonds.set();
// RDKit✔️✔️:   assignCIPLabels(mol, atoms, bonds, maxRecursiveIterations);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION assignCIPLabels all-molecule overload
pub(crate) fn assign_cip_labels(
    molecule: &Molecule,
    max_recursive_iterations: u32,
) -> Result<Molecule, CipLabelerError> {
    let atom_mask = vec![true; molecule.num_atoms()];
    let bond_mask = vec![true; molecule.num_bonds()];
    assign_cip_labels_for_masks(molecule, &atom_mask, &bond_mask, max_recursive_iterations)
}

pub(crate) struct CipAssignmentTransition {
    topology: crate::molecule::TopologyBlock,
    properties: crate::MoleculeProperties,
    error: Option<CipLabelerError>,
}

impl CipAssignmentTransition {
    pub(crate) fn into_parts(
        self,
    ) -> (
        crate::molecule::TopologyBlock,
        crate::MoleculeProperties,
        Option<CipLabelerError>,
    ) {
        (self.topology, self.properties, self.error)
    }
}

pub(crate) fn assign_cip_labels_from_read_parts(
    read: crate::ops::MoleculeReadParts<'_>,
    options: &crate::CipLabelOptions,
) -> Result<CipAssignmentTransition, CipLabelerError> {
    let molecule = Molecule::from_operation_blocks(
        read.topology().clone(),
        crate::molecule::CoordinateBlock::default(),
        read.properties().clone(),
        crate::molecule::DerivedCacheBlock::default(),
        read.capabilities(),
    )?;

    // BEGIN RDKIT CPP FUNCTION assignCIPLabelsWrapHelper selection dispatch
    // (CIPLabeler/Wrap/rdCIPLabeler.cpp:31-44; RDBoost/Wrap.cpp:75-90)
    // RDKit✔️✔️: auto atoms = pythonObjectToDynBitset(atomsToLabel, mol.getNumAtoms());
    // RDKit✔️✔️: auto bonds = pythonObjectToDynBitset(bondsToLabel, mol.getNumBonds());
    // RDKit✔️✔️: // If both atoms and bonds are None, assign all the mol.
    // RDKit✔️✔️: if (!atomsToLabel && !bondsToLabel) {
    // RDKit✔️✔️:   atoms.set();
    // RDKit✔️✔️:   bonds.set();
    // RDKit✔️✔️: }
    // RDKit✔️✔️: assignCIPLabels(mol, atoms, bonds, maxRecursiveIterations);
    // END RDKIT CPP FUNCTION assignCIPLabelsWrapHelper selection dispatch
    // Boost.Python truth-tests its objects here, so both omitted and empty
    // collections are false. The lower-level exact-mask function deliberately
    // retains the distinct C++ selected-overload semantics for internal users.
    let atoms = options.atom_indices();
    let bonds = options.bond_indices();
    let has_atom_selection = atoms.is_some_and(|indices| !indices.is_empty());
    let has_bond_selection = bonds.is_some_and(|indices| !indices.is_empty());
    let outcome = if !has_atom_selection && !has_bond_selection {
        let atom_mask = vec![true; molecule.num_atoms()];
        let bond_mask = vec![true; molecule.num_bonds()];
        assign_cip_labels_for_masks_outcome_with_cancellation(
            &molecule,
            &atom_mask,
            &bond_mask,
            options.max_recursive_iterations(),
            CipCancellationMode::NotRequested,
        )?
    } else {
        let atom_indices = atoms
            .unwrap_or_default()
            .iter()
            .map(|atom| {
                CipSourceIndex::try_from(atom.index()).map_err(|_| {
                    CipLabelerError::SourceIndexWidthExceeded {
                        kind: "atom",
                        index: atom.index(),
                    }
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        let bond_indices = bonds
            .unwrap_or_default()
            .iter()
            .map(|bond| {
                CipSourceIndex::try_from(bond.index()).map_err(|_| {
                    CipLabelerError::SourceIndexWidthExceeded {
                        kind: "bond",
                        index: bond.index(),
                    }
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        let (atom_mask, bond_mask) =
            cip_exact_selection_masks(&molecule, &atom_indices, &bond_indices)?;
        assign_cip_labels_for_masks_outcome_with_cancellation(
            &molecule,
            &atom_mask,
            &bond_mask,
            options.max_recursive_iterations(),
            CipCancellationMode::NotRequested,
        )?
    };
    Ok(CipAssignmentTransition {
        topology: outcome.molecule.topology_block().clone(),
        properties: outcome.molecule.properties().clone(),
        error: outcome.error,
    })
}

impl CipSequenceRule for CipRules {
    // BEGIN RDKIT CPP FUNCTION Rules::compare (rules/Rules.h)
    // RDKit✔️✔️: int compare(const Edge *o1, const Edge *o2) const override {
    // RDKit✔️✔️:   // Try using each rules. The rules will expand the search exhaustively
    // RDKit✔️✔️:   // to all child substituents
    // RDKit✔️✔️:   for (const auto &rule : d_rules) {
    // RDKit✔️✔️:     // compare expands exhaustively across the whole graph
    // RDKit✔️✔️:     int value = rule->recursiveCompare(o1, o2);
    // RDKit✔️✔️:     if (value != 0) {
    // RDKit✔️✔️:       return value;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return 0;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Rules::compare
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
    // RDKit✔️✔️: int getComparision(const Edge *a, const Edge *b, bool deep) const override {
    // RDKit✔️✔️:   (void)deep;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // Try using each rules. The rules will expand the search exhaustively
    // RDKit✔️✔️:   // to all child substituents
    // RDKit✔️✔️:   for (const auto &rule : d_rules) {
    // RDKit✔️✔️:     // compare expands exhaustively across the whole graph
    // RDKit✔️✔️:     int value = rule->recursiveCompare(a, b);
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (value != 0) {
    // RDKit✔️✔️:       return value;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return 0;
    // RDKit✔️✔️: }
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
    // RDKit✔️✔️: int getComparision(const Edge *a, const Edge *b, bool deep) const override {
    // RDKit✔️✔️:   (void)deep;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // Try using each rules. The rules will expand the search exhaustively
    // RDKit✔️✔️:   // to all child substituents
    // RDKit✔️✔️:   for (const auto &rule : d_rules) {
    // RDKit✔️✔️:     // compare expands exhaustively across the whole graph
    // RDKit✔️✔️:     int value = rule->recursiveCompare(a, b);
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (value != 0) {
    // RDKit✔️✔️:       return value;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return 0;
    // RDKit✔️✔️: }
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

pub(crate) struct CipConfiguration<'a> {
    foci: Vec<usize>,
    carriers: Vec<Option<usize>>,
    digraph: CipDigraph<'a>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct CipAtomPrimaryLabel {
    atom_idx: usize,
    cip_code: &'static str,
    cip_neighbor_order: Vec<CipSourceIndex>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct CipBondPrimaryLabel {
    bond_idx: usize,
    stereo_atoms: [usize; 2],
    stereo: BondStereo,
    cip_code: &'static str,
    cip_neighbor_order: Vec<CipSourceIndex>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct CipAtropisomerBondPrimaryLabel {
    bond_idx: usize,
    cip_code: &'static str,
    cip_neighbor_order: Vec<CipSourceIndex>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum CipPrimaryLabel {
    Atom(CipAtomPrimaryLabel),
    Bond(CipBondPrimaryLabel),
    AtropisomerBond(CipAtropisomerBondPrimaryLabel),
}

pub(crate) enum CipConfig<'a> {
    Tetrahedral(CipTetrahedral<'a>),
    Sp2Bond(CipSp2Bond<'a>),
    AtropisomerBond(CipAtropisomerBond<'a>),
}

impl<'a> CipConfig<'a> {
    fn get_foci(&self) -> &[usize] {
        match self {
            Self::Tetrahedral(config) => config.get_foci(),
            Self::Sp2Bond(config) => config.get_foci(),
            Self::AtropisomerBond(config) => config.get_foci(),
        }
    }

    fn get_digraph_mut(&mut self) -> &mut CipDigraph<'a> {
        match self {
            Self::Tetrahedral(config) => config.configuration.get_digraph(),
            Self::Sp2Bond(config) => config.configuration.get_digraph(),
            Self::AtropisomerBond(config) => config.configuration.get_digraph(),
        }
    }

    fn reset_primary_label(&mut self) {
        match self {
            Self::Tetrahedral(config) => config.reset_primary_label(),
            Self::Sp2Bond(config) => config.reset_primary_label(),
            Self::AtropisomerBond(config) => config.reset_primary_label(),
        }
    }

    fn has_primary_label(&self) -> bool {
        match self {
            Self::Tetrahedral(config) => config.has_primary_label(),
            Self::Sp2Bond(config) => config.has_primary_label(),
            Self::AtropisomerBond(config) => config.has_primary_label(),
        }
    }

    fn label(
        &mut self,
        rules: &CipRules,
        context: &mut CipLabelerContext,
    ) -> Result<Descriptor, CipLabelerError> {
        match self {
            Self::Tetrahedral(config) => config.label(rules, context),
            Self::Sp2Bond(config) => config.label(rules, context),
            Self::AtropisomerBond(config) => config.label(rules, context),
        }
    }

    fn label_with_external_digraph(
        &mut self,
        node: CipNodeId,
        digraph: &mut CipDigraph<'_>,
        rules: &CipRules,
        context: &mut CipLabelerContext,
    ) -> Result<Descriptor, CipLabelerError> {
        match self {
            Self::Tetrahedral(config) => {
                config.label_with_external_digraph(node, digraph, rules, context)
            }
            Self::Sp2Bond(config) => {
                config.label_with_external_digraph(node, digraph, rules, context)
            }
            Self::AtropisomerBond(config) => {
                config.label_with_external_digraph(node, digraph, rules, context)
            }
        }
    }

    fn set_primary_label(&mut self, desc: Descriptor) -> Result<(), CipLabelerError> {
        match self {
            Self::Tetrahedral(config) => config.set_primary_label(desc),
            Self::Sp2Bond(config) => config.set_primary_label(desc),
            Self::AtropisomerBond(config) => config.set_primary_label(desc),
        }
    }

    fn primary_label(&self) -> Option<CipPrimaryLabel> {
        match self {
            Self::Tetrahedral(config) => config.primary_label().cloned().map(CipPrimaryLabel::Atom),
            Self::Sp2Bond(config) => config.primary_label().cloned().map(CipPrimaryLabel::Bond),
            Self::AtropisomerBond(config) => config
                .primary_label()
                .cloned()
                .map(CipPrimaryLabel::AtropisomerBond),
        }
    }
}

impl<'a> CipConfiguration<'a> {
    pub(crate) const IMPLICIT_H: CipSourceIndex = CIP_NO_ATOM;

    // BEGIN RDKIT CPP FUNCTION Configuration::parity4 (configs/Configuration.h)
    // RDKit✔️✔️: template <typename T>
    // RDKit✔️✔️: static int parity4(const std::vector<T> &trg, const std::vector<T> &ref) {
    // RDKit✔️✔️:   if (ref.size() != 4 || trg.size() != ref.size()) {
    // RDKit✔️✔️:     throw std::runtime_error("Parity vectors must have size 4.");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (ref[0] == trg[0]) {
    // RDKit✔️✔️:     if (ref[1] == trg[1]) {
    // RDKit✔️✔️:       // a,b,c,d -> a,b,c,d
    // RDKit✔️✔️:       if (ref[2] == trg[2] && ref[3] == trg[3]) {
    // RDKit✔️✔️:         return 2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // a,b,c,d -> a,b,d,c
    // RDKit✔️✔️:       if (ref[2] == trg[3] && ref[3] == trg[2]) {
    // RDKit✔️✔️:         return 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (ref[1] == trg[2]) {
    // RDKit✔️✔️:       // a,b,c,d -> a,c,b,d
    // RDKit✔️✔️:       if (ref[2] == trg[1] && ref[3] == trg[3]) {
    // RDKit✔️✔️:         return 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // a,b,c,d -> a,c,d,b
    // RDKit✔️✔️:       if (ref[2] == trg[3] && ref[3] == trg[1]) {
    // RDKit✔️✔️:         return 2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (ref[1] == trg[3]) {
    // RDKit✔️✔️:       // a,b,c,d -> a,d,c,b
    // RDKit✔️✔️:       if (ref[2] == trg[2] && ref[3] == trg[1]) {
    // RDKit✔️✔️:         return 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // a,b,c,d -> a,d,b,c
    // RDKit✔️✔️:       if (ref[2] == trg[1] && ref[3] == trg[2]) {
    // RDKit✔️✔️:         return 2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (ref[0] == trg[1]) {
    // RDKit✔️✔️:     if (ref[1] == trg[0]) {
    // RDKit✔️✔️:       // a,b,c,d -> b,a,c,d
    // RDKit✔️✔️:       if (ref[2] == trg[2] && ref[3] == trg[3]) {
    // RDKit✔️✔️:         return 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // a,b,c,d -> b,a,d,c
    // RDKit✔️✔️:       if (ref[2] == trg[3] && ref[3] == trg[2]) {
    // RDKit✔️✔️:         return 2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (ref[1] == trg[2]) {
    // RDKit✔️✔️:       // a,b,c,d -> b,c,a,d
    // RDKit✔️✔️:       if (ref[2] == trg[0] && ref[3] == trg[3]) {
    // RDKit✔️✔️:         return 2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // a,b,c,d -> b,c,d,a
    // RDKit✔️✔️:       if (ref[2] == trg[3] && ref[3] == trg[0]) {
    // RDKit✔️✔️:         return 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (ref[1] == trg[3]) {
    // RDKit✔️✔️:       // a,b,c,d -> b,d,c,a
    // RDKit✔️✔️:       if (ref[2] == trg[2] && ref[3] == trg[0]) {
    // RDKit✔️✔️:         return 2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // a,b,c,d -> b,d,a,c
    // RDKit✔️✔️:       if (ref[2] == trg[0] && ref[3] == trg[2]) {
    // RDKit✔️✔️:         return 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (ref[0] == trg[2]) {
    // RDKit✔️✔️:     if (ref[1] == trg[1]) {
    // RDKit✔️✔️:       // a,b,c,d -> c,b,a,d
    // RDKit✔️✔️:       if (ref[2] == trg[0] && ref[3] == trg[3]) {
    // RDKit✔️✔️:         return 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // a,b,c,d -> c,b,d,a
    // RDKit✔️✔️:       if (ref[2] == trg[3] && ref[3] == trg[0]) {
    // RDKit✔️✔️:         return 2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (ref[1] == trg[0]) {
    // RDKit✔️✔️:       // a,b,c,d -> c,a,b,d
    // RDKit✔️✔️:       if (ref[2] == trg[1] && ref[3] == trg[3]) {
    // RDKit✔️✔️:         return 2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // a,b,c,d -> c,a,d,b
    // RDKit✔️✔️:       if (ref[2] == trg[3] && ref[3] == trg[1]) {
    // RDKit✔️✔️:         return 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (ref[1] == trg[3]) {
    // RDKit✔️✔️:       // a,b,c,d -> c,d,a,b
    // RDKit✔️✔️:       if (ref[2] == trg[0] && ref[3] == trg[1]) {
    // RDKit✔️✔️:         return 2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // a,b,c,d -> c,d,b,a
    // RDKit✔️✔️:       if (ref[2] == trg[1] && ref[3] == trg[0]) {
    // RDKit✔️✔️:         return 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (ref[0] == trg[3]) {
    // RDKit✔️✔️:     if (ref[1] == trg[1]) {
    // RDKit✔️✔️:       // a,b,c,d -> d,b,c,a
    // RDKit✔️✔️:       if (ref[2] == trg[2] && ref[3] == trg[0]) {
    // RDKit✔️✔️:         return 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // a,b,c,d -> d,b,a,c
    // RDKit✔️✔️:       if (ref[2] == trg[0] && ref[3] == trg[2]) {
    // RDKit✔️✔️:         return 2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (ref[1] == trg[2]) {
    // RDKit✔️✔️:       // a,b,c,d -> d,c,b,a
    // RDKit✔️✔️:       if (ref[2] == trg[1] && ref[3] == trg[0]) {
    // RDKit✔️✔️:         return 2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // a,b,c,d -> d,c,a,b
    // RDKit✔️✔️:       if (ref[2] == trg[0] && ref[3] == trg[1]) {
    // RDKit✔️✔️:         return 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (ref[1] == trg[0]) {
    // RDKit✔️✔️:       // a,b,c,d -> d,a,c,b
    // RDKit✔️✔️:       if (ref[2] == trg[2] && ref[3] == trg[1]) {
    // RDKit✔️✔️:         return 2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // a,b,c,d -> d,a,b,c
    // RDKit✔️✔️:       if (ref[2] == trg[1] && ref[3] == trg[2]) {
    // RDKit✔️✔️:         return 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // We should never hit this, but the compiler still complains
    // RDKit✔️✔️:   // about a missing return statement.
    // RDKit✔️✔️:   return 0;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Configuration::parity4
    pub(crate) fn parity4<T: PartialEq>(
        trg: &[T],
        reference: &[T],
    ) -> Result<i32, CipLabelerError> {
        if reference.len() != 4 || trg.len() != reference.len() {
            return Err(CipLabelerError::ParityVectorsMustHaveSize4);
        }

        let r = reference;
        let t = trg;
        if r[0] == t[0] {
            if r[1] == t[1] {
                if r[2] == t[2] && r[3] == t[3] {
                    return Ok(2);
                }
                if r[2] == t[3] && r[3] == t[2] {
                    return Ok(1);
                }
            } else if r[1] == t[2] {
                if r[2] == t[1] && r[3] == t[3] {
                    return Ok(1);
                }
                if r[2] == t[3] && r[3] == t[1] {
                    return Ok(2);
                }
            } else if r[1] == t[3] {
                if r[2] == t[2] && r[3] == t[1] {
                    return Ok(1);
                }
                if r[2] == t[1] && r[3] == t[2] {
                    return Ok(2);
                }
            }
        } else if r[0] == t[1] {
            if r[1] == t[0] {
                if r[2] == t[2] && r[3] == t[3] {
                    return Ok(1);
                }
                if r[2] == t[3] && r[3] == t[2] {
                    return Ok(2);
                }
            } else if r[1] == t[2] {
                if r[2] == t[0] && r[3] == t[3] {
                    return Ok(2);
                }
                if r[2] == t[3] && r[3] == t[0] {
                    return Ok(1);
                }
            } else if r[1] == t[3] {
                if r[2] == t[2] && r[3] == t[0] {
                    return Ok(2);
                }
                if r[2] == t[0] && r[3] == t[2] {
                    return Ok(1);
                }
            }
        } else if r[0] == t[2] {
            if r[1] == t[1] {
                if r[2] == t[0] && r[3] == t[3] {
                    return Ok(1);
                }
                if r[2] == t[3] && r[3] == t[0] {
                    return Ok(2);
                }
            } else if r[1] == t[0] {
                if r[2] == t[1] && r[3] == t[3] {
                    return Ok(2);
                }
                if r[2] == t[3] && r[3] == t[1] {
                    return Ok(1);
                }
            } else if r[1] == t[3] {
                if r[2] == t[0] && r[3] == t[1] {
                    return Ok(2);
                }
                if r[2] == t[1] && r[3] == t[0] {
                    return Ok(1);
                }
            }
        } else if r[0] == t[3] {
            if r[1] == t[1] {
                if r[2] == t[2] && r[3] == t[0] {
                    return Ok(1);
                }
                if r[2] == t[0] && r[3] == t[2] {
                    return Ok(2);
                }
            } else if r[1] == t[2] {
                if r[2] == t[1] && r[3] == t[0] {
                    return Ok(2);
                }
                if r[2] == t[0] && r[3] == t[1] {
                    return Ok(1);
                }
            } else if r[1] == t[0] {
                if r[2] == t[2] && r[3] == t[1] {
                    return Ok(2);
                }
                if r[2] == t[1] && r[3] == t[2] {
                    return Ok(1);
                }
            }
        }
        Ok(0)
    }

    // BEGIN RDKIT CPP FUNCTION Configuration constructors/getters/setCarriers/label (configs/Configuration.cpp)
    // RDKit✔️✔️: Configuration::Configuration(const CIPMol &mol, Atom *focus)
    // RDKit✔️✔️:     : d_foci{focus}, d_digraph{mol, focus} {};
    // RDKit✔️✔️:
    // RDKit✔️✔️: Configuration::Configuration(const CIPMol &mol, std::vector<Atom *> &&foci,
    // RDKit✔️✔️:                              bool atropisomerMode)
    // RDKit✔️✔️:     : d_foci{std::move(foci)}, d_digraph{mol, d_foci[0], atropisomerMode} {}
    // RDKit✔️✔️:
    // RDKit✔️✔️: Configuration::~Configuration() = default;
    // RDKit✔️✔️:
    // RDKit✔️✔️: Atom *Configuration::getFocus() const { return d_foci[0]; }
    // RDKit✔️✔️:
    // RDKit✔️✔️: const std::vector<Atom *> &Configuration::getFoci() const { return d_foci; }
    // RDKit✔️✔️:
    // RDKit✔️✔️: const std::vector<Atom *> &Configuration::getCarriers() const {
    // RDKit✔️✔️:   return d_carriers;
    // RDKit✔️✔️: }
    // RDKit✔️✔️:
    // RDKit✔️✔️: Digraph &Configuration::getDigraph() { return d_digraph; }
    // RDKit✔️✔️:
    // RDKit✔️✔️: Descriptor Configuration::label(Node *node, Digraph &digraph,
    // RDKit✔️✔️:                                 const Rules &comp) {
    // RDKit✔️✔️:   (void)node;
    // RDKit✔️✔️:   (void)digraph;
    // RDKit✔️✔️:   (void)comp;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return Descriptor::UNKNOWN;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Configuration constructors/getters/setCarriers/label
    pub(crate) fn new(molecule: &'a Molecule, focus: usize) -> Result<Self, CipLabelerError> {
        Ok(Self {
            foci: vec![focus],
            carriers: Vec::new(),
            digraph: CipDigraph::new(molecule, focus, false)?,
        })
    }

    pub(crate) fn with_foci(
        molecule: &'a Molecule,
        foci: Vec<usize>,
        atropisomer_mode: bool,
    ) -> Result<Self, CipLabelerError> {
        let focus = *foci
            .first()
            .ok_or(CipLabelerError::EmptyConfigurationFoci)?;
        Ok(Self {
            foci,
            carriers: Vec::new(),
            digraph: CipDigraph::new(molecule, focus, atropisomer_mode)?,
        })
    }

    pub(crate) fn get_focus(&self) -> usize {
        self.foci[0]
    }

    pub(crate) fn get_foci(&self) -> &[usize] {
        &self.foci
    }

    pub(crate) fn set_carriers(&mut self, carriers: Vec<Option<usize>>) {
        self.carriers = carriers;
    }

    pub(crate) fn get_carriers(&self) -> &[Option<usize>] {
        &self.carriers
    }

    pub(crate) fn get_digraph(&mut self) -> &mut CipDigraph<'a> {
        &mut self.digraph
    }

    pub(crate) fn label(
        &self,
        _node: CipNodeId,
        _digraph: &mut CipDigraph<'_>,
        _comp: &CipRules,
    ) -> Descriptor {
        Descriptor::Unknown
    }

    // BEGIN RDKIT CPP FUNCTION Configuration::findInternalEdge/isInternalEdge/removeInternalEdges (configs/Configuration.cpp)
    // RDKit✔️✔️: Edge *Configuration::findInternalEdge(const std::vector<Edge *> &edges,
    // RDKit✔️✔️:                                       Atom *f1, Atom *f2) {
    // RDKit✔️✔️:   for (const auto &edge : edges) {
    // RDKit✔️✔️:     if (edge->getBeg()->isDuplicate() || edge->getEnd()->isDuplicate()) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (isInternalEdge(edge, f1, f2)) {
    // RDKit✔️✔️:       return edge;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return nullptr;
    // RDKit✔️✔️: }
    // RDKit✔️✔️:
    // RDKit✔️✔️: bool Configuration::isInternalEdge(const Edge *edge, Atom *f1, Atom *f2) {
    // RDKit✔️✔️:   const auto &beg = edge->getBeg();
    // RDKit✔️✔️:   const auto &end = edge->getEnd();
    // RDKit✔️✔️:   if (f1 == beg->getAtom() && f2 == end->getAtom()) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   } else if (f1 == end->getAtom() && f2 == beg->getAtom()) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // RDKit✔️✔️:
    // RDKit✔️✔️: void Configuration::removeInternalEdges(std::vector<Edge *> &edges, Atom *f1,
    // RDKit✔️✔️:                                         Atom *f2) {
    // RDKit✔️✔️:   std::vector<Edge *> new_edges;
    // RDKit✔️✔️:   for (auto &&e : edges) {
    // RDKit✔️✔️:     if (!isInternalEdge(e, f1, f2)) {
    // RDKit✔️✔️:       new_edges.push_back(std::move(e));
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::swap(edges, new_edges);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Configuration::findInternalEdge/isInternalEdge/removeInternalEdges
    fn find_internal_edge(
        digraph: &CipDigraph<'_>,
        edges: &[CipEdgeId],
        f1: usize,
        f2: usize,
    ) -> Option<CipEdgeId> {
        edges.iter().copied().find(|edge_id| {
            let edge = digraph.edge(*edge_id);
            if digraph.node(edge.get_beg()).is_duplicate()
                || digraph.node(edge.get_end()).is_duplicate()
            {
                return false;
            }
            Self::is_internal_edge(digraph, *edge_id, f1, f2)
        })
    }

    fn is_internal_edge(
        digraph: &CipDigraph<'_>,
        edge_id: CipEdgeId,
        f1: usize,
        f2: usize,
    ) -> bool {
        let edge = digraph.edge(edge_id);
        let beg = digraph.node(edge.get_beg()).atom_idx();
        let end = digraph.node(edge.get_end()).atom_idx();
        (beg == Some(f1) && end == Some(f2)) || (beg == Some(f2) && end == Some(f1))
    }

    fn remove_internal_edges(
        digraph: &CipDigraph<'_>,
        edges: &mut Vec<CipEdgeId>,
        f1: usize,
        f2: usize,
    ) {
        edges.retain(|edge_id| !Self::is_internal_edge(digraph, *edge_id, f1, f2));
    }

    // BEGIN RDKIT CPP FUNCTION Configuration::isDuplicateOrHydrogenEdge/removeDuplicatesAndHs (configs/Configuration.cpp)
    // RDKit✔️✔️: bool Configuration::isDuplicateOrHydrogenEdge(const Edge *edge) {
    // RDKit✔️✔️:   return edge->getBeg()->isDuplicateOrH() || edge->getEnd()->isDuplicateOrH();
    // RDKit✔️✔️: }
    // RDKit✔️✔️:
    // RDKit✔️✔️: void Configuration::removeDuplicatesAndHs(std::vector<Edge *> &edges) {
    // RDKit✔️✔️:   std::vector<Edge *> new_edges;
    // RDKit✔️✔️:   for (auto &&e : edges) {
    // RDKit✔️✔️:     if (!isDuplicateOrHydrogenEdge(e)) {
    // RDKit✔️✔️:       new_edges.push_back(std::move(e));
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   edges = std::move(new_edges);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Configuration::isDuplicateOrHydrogenEdge/removeDuplicatesAndHs
    fn is_duplicate_or_hydrogen_edge(digraph: &CipDigraph<'_>, edge_id: CipEdgeId) -> bool {
        let edge = digraph.edge(edge_id);
        digraph.node(edge.get_beg()).is_duplicate_or_h()
            || digraph.node(edge.get_end()).is_duplicate_or_h()
    }

    fn remove_duplicates_and_hs(digraph: &CipDigraph<'_>, edges: &mut Vec<CipEdgeId>) {
        edges.retain(|edge_id| !Self::is_duplicate_or_hydrogen_edge(digraph, *edge_id));
    }
}

pub(crate) struct CipTetrahedral<'a> {
    configuration: CipConfiguration<'a>,
    ranked_anchors: Vec<CipSourceIndex>,
    primary_label: Option<CipAtomPrimaryLabel>,
    source_primary_label_visible: bool,
}

impl<'a> CipTetrahedral<'a> {
    // BEGIN RDKIT CPP FUNCTION Tetrahedral::Tetrahedral (configs/Tetrahedral.cpp)
    // RDKit✔️✔️: Tetrahedral::Tetrahedral(const CIPMol &mol, Atom *focus)
    // RDKit✔️✔️:     : Configuration(mol, focus) {
    // RDKit✔️✔️:   CHECK_INVARIANT(focus, "bad atom")
    // RDKit✔️✔️:   CHECK_INVARIANT(focus->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW ||
    // RDKit✔️✔️:                       focus->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW,
    // RDKit✔️✔️:                   "bad config")
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::vector<Atom *> carriers;
    // RDKit✔️✔️:   carriers.reserve(4);
    // RDKit✔️✔️:   for (auto &nbr : mol.getNeighbors(focus)) {
    // RDKit✔️✔️:     carriers.push_back(nbr);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (carriers.size() < 4) {
    // RDKit✔️✔️:     // Implicit H -- use the central atom instead of a dummy H
    // RDKit✔️✔️:     carriers.push_back(focus);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (carriers.size() < 4) {
    // RDKit✔️✔️:     // Trigonal pyramid centers with an implicit H need a phantom
    // RDKit✔️✔️:     // atom as fourth carrier. This one must be represented differently
    // RDKit✔️✔️:     // than the implicit H.
    // RDKit✔️✔️:     carriers.push_back(nullptr);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   POSTCONDITION(carriers.size() == 4, "configuration must have 4 carriers");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   setCarriers(std::move(carriers));
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION Tetrahedral::Tetrahedral
    pub(crate) fn new(molecule: &'a Molecule, focus: usize) -> Result<Self, CipLabelerError> {
        let mut configuration = CipConfiguration::new(molecule, focus)?;
        let atom = configuration.digraph.mol().atom(focus)?;
        let source_primary_label_visible = atom.prop("_CIPCode").is_some();
        match atom.chiral_tag() {
            ChiralTag::TetrahedralCcw | ChiralTag::TetrahedralCw => {}
            ChiralTag::Unspecified
            | ChiralTag::Other
            | ChiralTag::Tetrahedral
            | ChiralTag::Allene
            | ChiralTag::SquarePlanar
            | ChiralTag::TrigonalBipyramidal
            | ChiralTag::Octahedral => return Err(CipLabelerError::BadTetrahedralConfig),
        }

        let mut carriers = Vec::with_capacity(4);
        for nbr in configuration.digraph.mol().neighbor_indices(focus)? {
            carriers.push(Some(nbr));
        }
        if carriers.len() < 4 {
            carriers.push(Some(focus));
        }
        if carriers.len() < 4 {
            carriers.push(None);
        }
        if carriers.len() != 4 {
            return Err(CipLabelerError::TetrahedralConfigurationMustHave4Carriers);
        }

        configuration.set_carriers(carriers);
        Ok(Self {
            configuration,
            ranked_anchors: Vec::new(),
            primary_label: None,
            source_primary_label_visible,
        })
    }

    pub(crate) fn get_focus(&self) -> usize {
        self.configuration.get_focus()
    }

    pub(crate) fn get_foci(&self) -> &[usize] {
        self.configuration.get_foci()
    }

    pub(crate) fn get_carriers(&self) -> &[Option<usize>] {
        self.configuration.get_carriers()
    }

    pub(crate) fn ranked_anchors(&self) -> &[CipSourceIndex] {
        &self.ranked_anchors
    }

    pub(crate) fn primary_label(&self) -> Option<&CipAtomPrimaryLabel> {
        self.primary_label.as_ref()
    }

    // BEGIN RDKIT CPP FUNCTION Tetrahedral::setPrimaryLabel (configs/Tetrahedral.cpp)
    // RDKit✔️✔️: void Tetrahedral::setPrimaryLabel(Descriptor desc) {
    // RDKit✔️✔️:   switch (desc) {
    // RDKit✔️✔️:     case Descriptor::R:
    // RDKit✔️✔️:     case Descriptor::S:
    // RDKit✔️✔️:     case Descriptor::r:
    // RDKit✔️✔️:     case Descriptor::s: {
    // RDKit✔️✔️:       auto chiralAtom = getFocus();
    // RDKit✔️✔️:       chiralAtom->setProp(common_properties::_CIPCode, to_string(desc));
    // RDKit✔️✔️:       chiralAtom->setProp(common_properties::_CIPNeighborOrder,
    // RDKit✔️✔️:                           d_ranked_anchors, true);
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     case Descriptor::seqTrans:
    // RDKit✔️✔️:     case Descriptor::seqCis:
    // RDKit✔️✔️:     case Descriptor::E:
    // RDKit✔️✔️:     case Descriptor::Z:
    // RDKit✔️✔️:     case Descriptor::M:
    // RDKit✔️✔️:     case Descriptor::P:
    // RDKit✔️✔️:     case Descriptor::m:
    // RDKit✔️✔️:     case Descriptor::p:
    // RDKit✔️✔️:     case Descriptor::SP_4:
    // RDKit✔️✔️:     case Descriptor::TBPY_5:
    // RDKit✔️✔️:     case Descriptor::OC_6:
    // RDKit✔️✔️:       throw std::runtime_error(
    // RDKit✔️✔️:           "Received a Descriptor that is not supported for atoms");
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       throw std::runtime_error("Received an invalid Atom Descriptor");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Tetrahedral::setPrimaryLabel
    pub(crate) fn set_primary_label(&mut self, desc: Descriptor) -> Result<(), CipLabelerError> {
        match desc {
            Descriptor::R | Descriptor::S | Descriptor::r | Descriptor::s => {
                self.primary_label = Some(CipAtomPrimaryLabel {
                    atom_idx: self.configuration.get_focus(),
                    cip_code: descriptor_to_string(desc),
                    cip_neighbor_order: self.ranked_anchors.clone(),
                });
                self.source_primary_label_visible = false;
                Ok(())
            }
            Descriptor::seqTrans
            | Descriptor::seqCis
            | Descriptor::E
            | Descriptor::Z
            | Descriptor::M
            | Descriptor::P
            | Descriptor::m
            | Descriptor::p
            | Descriptor::SP_4
            | Descriptor::TBPY_5
            | Descriptor::OC_6 => Err(CipLabelerError::DescriptorNotSupportedForAtoms),
            Descriptor::None | Descriptor::Unknown | Descriptor::ns => {
                Err(CipLabelerError::InvalidAtomDescriptor)
            }
        }
    }

    // BEGIN RDKIT CPP FUNCTION Tetrahedral::hasPrimaryLabel (configs/Tetrahedral.cpp)
    // RDKit✔️✔️: bool Tetrahedral::hasPrimaryLabel() const {
    // RDKit✔️✔️:   return getFocus()->hasProp(common_properties::_CIPCode);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Tetrahedral::hasPrimaryLabel
    pub(crate) fn has_primary_label(&self) -> bool {
        self.primary_label.is_some() || self.source_primary_label_visible
    }

    // BEGIN RDKIT CPP FUNCTION Tetrahedral::resetPrimaryLabel (configs/Tetrahedral.cpp)
    // RDKit✔️✔️: void Tetrahedral::resetPrimaryLabel() const {
    // RDKit✔️✔️:   getFocus()->clearProp(common_properties::_CIPCode);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Tetrahedral::resetPrimaryLabel
    pub(crate) fn reset_primary_label(&mut self) {
        self.primary_label = None;
        self.source_primary_label_visible = false;
    }

    // BEGIN RDKIT CPP FUNCTION Tetrahedral::label(const Rules &) (configs/Tetrahedral.cpp)
    // RDKit✔️✔️: Descriptor Tetrahedral::label(const Rules &comp) {
    // RDKit✔️✔️:   auto &digraph = getDigraph();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto root = digraph.getOriginalRoot();
    // RDKit✔️✔️:   if (digraph.getCurrentRoot() != root) {
    // RDKit✔️✔️:     digraph.changeRoot(root);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return label(root, comp);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Tetrahedral::label(const Rules &)
    pub(crate) fn label(
        &mut self,
        comp: &CipRules,
        context: &mut CipLabelerContext,
    ) -> Result<Descriptor, CipLabelerError> {
        let root = self.configuration.digraph.get_original_root();
        if self.configuration.digraph.get_current_root() != root {
            self.configuration.digraph.change_root(root)?;
        }
        self.label_node(root, comp, context)
    }

    // BEGIN RDKIT CPP FUNCTION Tetrahedral::label(Node *, Digraph &, const Rules &) (configs/Tetrahedral.cpp)
    // RDKit✔️✔️: Descriptor Tetrahedral::label(Node *node, Digraph &digraph, const Rules &comp) {
    // RDKit✔️✔️:   digraph.changeRoot(node);
    // RDKit✔️✔️:   return label(node, comp);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Tetrahedral::label(Node *, Digraph &, const Rules &)
    pub(crate) fn label_with_external_digraph(
        &mut self,
        node: CipNodeId,
        digraph: &mut CipDigraph<'_>,
        comp: &CipRules,
        context: &mut CipLabelerContext,
    ) -> Result<Descriptor, CipLabelerError> {
        digraph.change_root(node)?;
        self.label_node_in_digraph(node, digraph, comp, context)
    }

    fn label_node(
        &mut self,
        node: CipNodeId,
        comp: &CipRules,
        context: &mut CipLabelerContext,
    ) -> Result<Descriptor, CipLabelerError> {
        let focus = self.configuration.get_focus();
        let carriers = self.configuration.get_carriers().to_vec();
        Self::label_node_impl(
            focus,
            &carriers,
            &mut self.ranked_anchors,
            node,
            &mut self.configuration.digraph,
            comp,
            context,
        )
    }

    // BEGIN RDKIT CPP FUNCTION Tetrahedral::label(Node *, const Rules &) (configs/Tetrahedral.cpp)
    // RDKit✔️✔️: Descriptor Tetrahedral::label(Node *node, const Rules &comp) {
    // RDKit✔️✔️:   auto focus = getFocus();
    // RDKit✔️✔️:   auto edges = node->getEdges();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   d_ranked_anchors.clear();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // something not right!?! bad creation
    // RDKit✔️✔️:   if (edges.size() < 3) {
    // RDKit✔️✔️:     return Descriptor::ns;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto priority = comp.sort(node, edges);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   bool isUnique = priority.isUnique();
    // RDKit✔️✔️:   if (!isUnique && edges.size() == 4) {
    // RDKit✔️✔️:     if (comp.getNumSubRules() == 3) {
    // RDKit✔️✔️:       return Descriptor::UNKNOWN;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     auto partition = comp.getSorter()->getGroups(edges);
    // RDKit✔️✔️:     if (partition.size() == 2) {
    // RDKit✔️✔️:       node->getDigraph()->setRule6Ref(edges[1]->getEnd()->getAtom());
    // RDKit✔️✔️:       priority = comp.sort(node, edges);
    // RDKit✔️✔️:       node->getDigraph()->setRule6Ref(nullptr);
    // RDKit✔️✔️:     } else if (partition.size() == 1) {
    // RDKit✔️✔️:       // S4 symmetric case
    // RDKit✔️✔️:       node->getDigraph()->setRule6Ref(edges[0]->getEnd()->getAtom());
    // RDKit✔️✔️:       comp.sort(node, edges);
    // RDKit✔️✔️:       auto nbrs1 = std::vector<Edge *>(edges.begin(), edges.end());
    // RDKit✔️✔️:
    // RDKit✔️✔️:       node->getDigraph()->setRule6Ref(edges[1]->getEnd()->getAtom());
    // RDKit✔️✔️:       priority = comp.sort(node, edges);
    // RDKit✔️✔️:
    // RDKit✔️✔️:       node->getDigraph()->setRule6Ref(nullptr);
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (parity4(nbrs1, edges) == 1) {
    // RDKit✔️✔️:         return Descriptor::UNKNOWN;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (!priority.isUnique()) {
    // RDKit✔️✔️:       return Descriptor::UNKNOWN;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (!isUnique) {
    // RDKit✔️✔️:     return Descriptor::UNKNOWN;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto ordered = std::vector<Atom *>(4, nullptr);
    // RDKit✔️✔️:   int idx = 0;
    // RDKit✔️✔️:   d_ranked_anchors.reserve(4);
    // RDKit✔️✔️:   for (const auto &edge : edges) {
    // RDKit✔️✔️:     if (edge->getEnd()->isSet(Node::BOND_DUPLICATE) ||
    // RDKit✔️✔️:         edge->getEnd()->isSet(Node::IMPL_HYDROGEN)) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     auto atom = edge->getEnd()->getAtom();
    // RDKit✔️✔️:     ordered[idx] = atom;
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // In this case we don't worry about implicit H (see Sp2Bond
    // RDKit✔️✔️:     // and Atropisomer): chirality is positional, and we don't
    // RDKit✔️✔️:     // know where the implicit H may be ("before" or "after" a
    // RDKit✔️✔️:     // potential 1H with lower priority?), so we just ignore it
    // RDKit✔️✔️:     // in the ranked neighbors list
    // RDKit✔️✔️:     d_ranked_anchors.push_back(atom->getIdx());
    // RDKit✔️✔️:
    // RDKit✔️✔️:     ++idx;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // if we are resolving a trigonal pyramid with an implicit H,
    // RDKit✔️✔️:   // the 4th carrier will be a nullptr: we need to add a phantom
    // RDKit✔️✔️:   // atom, which will always have the lowest priority, so that
    // RDKit✔️✔️:   // it must be different than the representation of the implicit H.
    // RDKit✔️✔️:   if (idx < 4) {
    // RDKit✔️✔️:     ordered[idx] = focus;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   int parity = parity4(ordered, getCarriers());
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (parity == 0) {
    // RDKit✔️✔️:     throw std::runtime_error("Could not calculate parity! Carrier mismatch");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto config = focus->getChiralTag();
    // RDKit✔️✔️:   if (parity == 1) {
    // RDKit✔️✔️:     if (config == Atom::CHI_TETRAHEDRAL_CCW) {
    // RDKit✔️✔️:       config = Atom::CHI_TETRAHEDRAL_CW;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       config = Atom::CHI_TETRAHEDRAL_CCW;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (config == Atom::CHI_TETRAHEDRAL_CCW) {
    // RDKit✔️✔️:     if (priority.isPseudoAsymetric()) {
    // RDKit✔️✔️:       return Descriptor::s;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       return Descriptor::S;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (config == Atom::CHI_TETRAHEDRAL_CW) {
    // RDKit✔️✔️:     if (priority.isPseudoAsymetric()) {
    // RDKit✔️✔️:       return Descriptor::r;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       return Descriptor::R;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return Descriptor::UNKNOWN;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Tetrahedral::label(Node *, const Rules &)
    fn label_node_in_digraph(
        &mut self,
        node: CipNodeId,
        digraph: &mut CipDigraph<'_>,
        comp: &CipRules,
        context: &mut CipLabelerContext,
    ) -> Result<Descriptor, CipLabelerError> {
        let focus = self.configuration.get_focus();
        let carriers = self.configuration.get_carriers().to_vec();
        Self::label_node_impl(
            focus,
            &carriers,
            &mut self.ranked_anchors,
            node,
            digraph,
            comp,
            context,
        )
    }

    fn label_node_impl(
        focus: usize,
        carriers: &[Option<usize>],
        ranked_anchors: &mut Vec<CipSourceIndex>,
        node: CipNodeId,
        digraph: &mut CipDigraph<'_>,
        comp: &CipRules,
        context: &mut CipLabelerContext,
    ) -> Result<Descriptor, CipLabelerError> {
        let mut edges = digraph.node_edges(node)?;
        ranked_anchors.clear();

        if edges.len() < 3 {
            return Ok(Descriptor::ns);
        }

        let mut priority = comp.sort(digraph, context, node, &mut edges, true)?;
        let is_unique = priority.is_unique();
        if !is_unique && edges.len() == 4 {
            if comp.get_num_sub_rules() == 3 {
                return Ok(Descriptor::Unknown);
            }
            let partition = comp.get_sorter().get_groups(digraph, context, &edges)?;
            if partition.len() == 2 {
                let ref_atom = digraph.node(digraph.edge(edges[1]).get_end()).atom_idx();
                digraph.set_rule6_ref(ref_atom)?;
                priority = comp.sort(digraph, context, node, &mut edges, true)?;
                digraph.set_rule6_ref(None)?;
            } else if partition.len() == 1 {
                let ref_atom = digraph.node(digraph.edge(edges[0]).get_end()).atom_idx();
                digraph.set_rule6_ref(ref_atom)?;
                comp.sort(digraph, context, node, &mut edges, true)?;
                let nbrs1 = edges.clone();

                let ref_atom = digraph.node(digraph.edge(edges[1]).get_end()).atom_idx();
                digraph.set_rule6_ref(ref_atom)?;
                priority = comp.sort(digraph, context, node, &mut edges, true)?;

                digraph.set_rule6_ref(None)?;

                if CipConfiguration::parity4(&nbrs1, &edges)? == 1 {
                    return Ok(Descriptor::Unknown);
                }
            }
            if !priority.is_unique() {
                return Ok(Descriptor::Unknown);
            }
        } else if !is_unique {
            return Ok(Descriptor::Unknown);
        }

        let mut ordered = vec![None; 4];
        let mut idx = 0_usize;
        ranked_anchors.reserve(4);
        for edge in &edges {
            let end = digraph.edge(*edge).get_end();
            if digraph.node(end).is_set(CipNode::BOND_DUPLICATE)
                || digraph.node(end).is_set(CipNode::IMPL_HYDROGEN)
            {
                continue;
            }

            let atom = digraph.node(end).atom_idx();
            if idx < 4 {
                ordered[idx] = atom;
            }
            if let Some(atom_idx) = atom {
                ranked_anchors.push(CipSourceIndex::try_from(atom_idx).map_err(|_| {
                    CipLabelerError::SourceIndexWidthExceeded {
                        kind: "atom",
                        index: atom_idx,
                    }
                })?);
            }
            idx += 1;
        }

        if idx < 4 {
            ordered[idx] = Some(focus);
        }

        let parity = CipConfiguration::parity4(&ordered, carriers)?;
        if parity == 0 {
            return Err(CipLabelerError::CarrierMismatch);
        }

        let mut config = digraph.mol().atom(focus)?.chiral_tag();
        if parity == 1 {
            config = match config {
                ChiralTag::TetrahedralCcw => ChiralTag::TetrahedralCw,
                ChiralTag::TetrahedralCw => ChiralTag::TetrahedralCcw,
                _ => config,
            };
        }

        if config == ChiralTag::TetrahedralCcw {
            if priority.is_pseudo_asymetric() {
                Ok(Descriptor::s)
            } else {
                Ok(Descriptor::S)
            }
        } else if config == ChiralTag::TetrahedralCw {
            if priority.is_pseudo_asymetric() {
                Ok(Descriptor::r)
            } else {
                Ok(Descriptor::R)
            }
        } else {
            Ok(Descriptor::Unknown)
        }
    }
}

pub(crate) struct CipSp2Bond<'a> {
    configuration: CipConfiguration<'a>,
    bond_idx: usize,
    cfg: BondStereo,
    ranked_anchors: Vec<CipSourceIndex>,
    primary_label: Option<CipBondPrimaryLabel>,
    source_primary_label_visible: bool,
}

impl<'a> CipSp2Bond<'a> {
    // BEGIN RDKIT CPP FUNCTION Sp2Bond::Sp2Bond (configs/Sp2Bond.cpp)
    // RDKit✔️✔️: Sp2Bond::Sp2Bond(const CIPMol &mol, Bond *bond, Atom *startAtom, Atom *endAtom,
    // RDKit✔️✔️:                  Bond::BondStereo cfg)
    // RDKit✔️✔️:     : Configuration(mol, {startAtom, endAtom}), dp_bond{bond}, d_cfg{cfg} {
    // RDKit✔️✔️:   CHECK_INVARIANT(startAtom && endAtom, "bad foci")
    // RDKit✔️✔️:   CHECK_INVARIANT(d_cfg == Bond::STEREOTRANS || d_cfg == Bond::STEREOCIS,
    // RDKit✔️✔️:                   "bad config")
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto stereo_atoms = Chirality::findStereoAtoms(bond);
    // RDKit✔️✔️:   CHECK_INVARIANT(stereo_atoms.size() == 2, "incorrect number of stereo atoms")
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::vector<Atom *> anchors{
    // RDKit✔️✔️:       {mol.getAtom(stereo_atoms[0]), mol.getAtom(stereo_atoms[1])}};
    // RDKit✔️✔️:
    // RDKit✔️✔️:   setCarriers(std::move(anchors));
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Sp2Bond::Sp2Bond
    // BEGIN RDKIT CPP FUNCTION Chirality::findStereoAtoms (GraphMol/Chirality.cpp)
    // RDKit✔️✔️: INT_VECT findStereoAtoms(const Bond *bond) {
    // RDKit✔️✔️:   PRECONDITION(bond, "bad bond");
    // RDKit✔️✔️:   PRECONDITION(bond->hasOwningMol(), "no mol");
    // RDKit✔️✔️:   PRECONDITION(bond->getBondType() == Bond::DOUBLE, "not double bond");
    // RDKit✔️✔️:   PRECONDITION(bond->getStereo() > Bond::BondStereo::STEREOANY,
    // RDKit✔️✔️:                "no defined stereo");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (!bond->getStereoAtoms().empty()) {
    // RDKit✔️✔️:     return bond->getStereoAtoms();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (bond->getStereo() == Bond::BondStereo::STEREOE ||
    // RDKit✔️✔️:       bond->getStereo() == Bond::BondStereo::STEREOZ) {
    // RDKit✔️✔️:     const Atom *startStereoAtom =
    // RDKit✔️✔️:         findHighestCIPNeighbor(bond->getBeginAtom(), bond->getEndAtom());
    // RDKit✔️✔️:     const Atom *endStereoAtom =
    // RDKit✔️✔️:         findHighestCIPNeighbor(bond->getEndAtom(), bond->getBeginAtom());
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (startStereoAtom == nullptr || endStereoAtom == nullptr) {
    // RDKit✔️✔️:       return {};
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     int startStereoAtomIdx = static_cast<int>(startStereoAtom->getIdx());
    // RDKit✔️✔️:     int endStereoAtomIdx = static_cast<int>(endStereoAtom->getIdx());
    // RDKit✔️✔️:
    // RDKit✔️✔️:     return {startStereoAtomIdx, endStereoAtomIdx};
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     BOOST_LOG(rdWarningLog) << "Unable to assign stereo atoms for bond "
    // RDKit✔️✔️:                             << bond->getIdx() << std::endl;
    // RDKit✔️✔️:     return {};
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Chirality::findStereoAtoms
    pub(crate) fn new(
        molecule: &'a Molecule,
        bond_idx: usize,
        start_atom: usize,
        end_atom: usize,
        cfg: BondStereo,
    ) -> Result<Self, CipLabelerError> {
        let mut configuration =
            CipConfiguration::with_foci(molecule, vec![start_atom, end_atom], false)?;
        configuration.digraph.mol().atom(start_atom)?;
        configuration.digraph.mol().atom(end_atom)?;
        let bond = configuration.digraph.mol().bond(bond_idx)?;
        let source_primary_label_visible = bond.prop("_CIPCode").is_some();
        if !((bond.begin().index() == start_atom && bond.end().index() == end_atom)
            || (bond.begin().index() == end_atom && bond.end().index() == start_atom))
        {
            return Err(CipLabelerError::BadSp2BondFoci);
        }
        if !matches!(cfg, BondStereo::Trans | BondStereo::Cis) {
            return Err(CipLabelerError::BadSp2BondConfig);
        }
        if bond.order() != BondOrder::Double {
            return Err(CipLabelerError::Sp2BondNotDoubleBond);
        }
        if matches!(bond.stereo(), BondStereo::None | BondStereo::Any) {
            return Err(CipLabelerError::Sp2BondHasNoDefinedStereo);
        }

        let stereo_atoms = if let Some([left, right]) = bond.stereo_atoms() {
            [left.index(), right.index()]
        } else if matches!(bond.stereo(), BondStereo::E | BondStereo::Z) {
            let left = Self::find_highest_cip_neighbor_like_rdkit(
                configuration.digraph.mol(),
                start_atom,
                end_atom,
            )?;
            let right = Self::find_highest_cip_neighbor_like_rdkit(
                configuration.digraph.mol(),
                end_atom,
                start_atom,
            )?;
            match (left, right) {
                (Some(left), Some(right)) => [left, right],
                _ => return Err(CipLabelerError::IncorrectNumberOfStereoAtoms),
            }
        } else {
            return Err(CipLabelerError::IncorrectNumberOfStereoAtoms);
        };
        configuration.digraph.mol().atom(stereo_atoms[0])?;
        configuration.digraph.mol().atom(stereo_atoms[1])?;
        configuration.set_carriers(vec![Some(stereo_atoms[0]), Some(stereo_atoms[1])]);

        Ok(Self {
            configuration,
            bond_idx,
            cfg,
            ranked_anchors: Vec::new(),
            primary_label: None,
            source_primary_label_visible,
        })
    }

    // BEGIN RDKIT CPP FUNCTION findHighestCIPNeighbor (GraphMol/Chirality.cpp)
    // RDKit✔️✔️: const Atom *findHighestCIPNeighbor(const Atom *atom, const Atom *skipAtom) {
    // RDKit✔️✔️:   PRECONDITION(atom, "bad atom");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   unsigned bestCipRank = 0;
    // RDKit✔️✔️:   const Atom *bestCipRankedAtom = nullptr;
    // RDKit✔️✔️:   const auto &mol = atom->getOwningMol();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (const auto neighbor : mol.atomNeighbors(atom)) {
    // RDKit✔️✔️:     if (neighbor == skipAtom) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     unsigned cip = 0;
    // RDKit✔️✔️:     if (!neighbor->getPropIfPresent(common_properties::_CIPRank, cip)) {
    // RDKit✔️✔️:       // If at least one of the atoms doesn't have a CIP rank, the highest rank
    // RDKit✔️✔️:       // does not make sense, so return a nullptr.
    // RDKit✔️✔️:       return nullptr;
    // RDKit✔️✔️:     } else if (cip > bestCipRank || bestCipRankedAtom == nullptr) {
    // RDKit✔️✔️:       bestCipRank = cip;
    // RDKit✔️✔️:       bestCipRankedAtom = neighbor;
    // RDKit✔️✔️:     } else if (cip == bestCipRank) {
    // RDKit✔️✔️:       // This also doesn't make sense if there is a tie (if that's possible).
    // RDKit✔️✔️:       // We still keep the best CIP rank in case something better comes around
    // RDKit✔️✔️:       // (also not sure if that's possible).
    // RDKit✔️✔️:       BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:           << "Warning: duplicate CIP ranks found in findHighestCIPNeighbor()"
    // RDKit✔️✔️:           << std::endl;
    // RDKit✔️✔️:       bestCipRankedAtom = nullptr;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return bestCipRankedAtom;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION findHighestCIPNeighbor
    fn find_highest_cip_neighbor_like_rdkit(
        mol: &CipMol<'_>,
        atom_idx: usize,
        skip_atom_idx: usize,
    ) -> Result<Option<usize>, CipLabelerError> {
        mol.atom(atom_idx)?;
        let mut best_cip_rank = 0_u32;
        let mut best_cip_ranked_atom = None;
        for neighbor in mol.neighbor_indices(atom_idx)? {
            if neighbor == skip_atom_idx {
                continue;
            }
            let Some(rank_text) = mol.atom(neighbor)?.prop("_CIPRank") else {
                return Ok(None);
            };
            let Ok(cip) = rank_text.parse::<u32>() else {
                return Ok(None);
            };
            if cip > best_cip_rank || best_cip_ranked_atom.is_none() {
                best_cip_rank = cip;
                best_cip_ranked_atom = Some(neighbor);
            } else if cip == best_cip_rank {
                best_cip_ranked_atom = None;
            }
        }
        Ok(best_cip_ranked_atom)
    }

    pub(crate) fn get_foci(&self) -> &[usize] {
        self.configuration.get_foci()
    }

    pub(crate) fn get_carriers(&self) -> &[Option<usize>] {
        self.configuration.get_carriers()
    }

    pub(crate) fn ranked_anchors(&self) -> &[CipSourceIndex] {
        &self.ranked_anchors
    }

    pub(crate) fn primary_label(&self) -> Option<&CipBondPrimaryLabel> {
        self.primary_label.as_ref()
    }

    // BEGIN RDKIT CPP FUNCTION Sp2Bond::setPrimaryLabel (configs/Sp2Bond.cpp)
    // RDKit✔️✔️: void Sp2Bond::setPrimaryLabel(Descriptor desc) {
    // RDKit✔️✔️:   switch (desc) {
    // RDKit✔️✔️:     case Descriptor::seqTrans:
    // RDKit✔️✔️:     case Descriptor::E:
    // RDKit✔️✔️:     case Descriptor::seqCis:
    // RDKit✔️✔️:     case Descriptor::Z: {
    // RDKit✔️✔️:       auto carriers = getCarriers();
    // RDKit✔️✔️:       dp_bond->setStereoAtoms(carriers[0]->getIdx(), carriers[1]->getIdx());
    // RDKit✔️✔️:       dp_bond->setStereo(d_cfg);
    // RDKit✔️✔️:       dp_bond->setProp(common_properties::_CIPCode, to_string(desc));
    // RDKit✔️✔️:       dp_bond->setProp(common_properties::_CIPNeighborOrder, d_ranked_anchors,
    // RDKit✔️✔️:                        true);
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     case Descriptor::R:
    // RDKit✔️✔️:     case Descriptor::S:
    // RDKit✔️✔️:     case Descriptor::r:
    // RDKit✔️✔️:     case Descriptor::s:
    // RDKit✔️✔️:     case Descriptor::M:
    // RDKit✔️✔️:     case Descriptor::P:
    // RDKit✔️✔️:     case Descriptor::m:
    // RDKit✔️✔️:     case Descriptor::p:
    // RDKit✔️✔️:     case Descriptor::SP_4:
    // RDKit✔️✔️:     case Descriptor::TBPY_5:
    // RDKit✔️✔️:     case Descriptor::OC_6:
    // RDKit✔️✔️:       throw std::runtime_error(
    // RDKit✔️✔️:           "Received a Descriptor that is not supported for double bonds");
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       throw std::runtime_error("Received an invalid Bond Descriptor");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Sp2Bond::setPrimaryLabel
    pub(crate) fn set_primary_label(&mut self, desc: Descriptor) -> Result<(), CipLabelerError> {
        match desc {
            Descriptor::seqTrans | Descriptor::E | Descriptor::seqCis | Descriptor::Z => {
                let carriers = self.configuration.get_carriers();
                let stereo_atoms = [
                    carriers
                        .first()
                        .and_then(|carrier| *carrier)
                        .ok_or(CipLabelerError::IncorrectNumberOfStereoAtoms)?,
                    carriers
                        .get(1)
                        .and_then(|carrier| *carrier)
                        .ok_or(CipLabelerError::IncorrectNumberOfStereoAtoms)?,
                ];
                self.primary_label = Some(CipBondPrimaryLabel {
                    bond_idx: self.bond_idx,
                    stereo_atoms,
                    stereo: self.cfg,
                    cip_code: descriptor_to_string(desc),
                    cip_neighbor_order: self.ranked_anchors.clone(),
                });
                self.source_primary_label_visible = false;
                Ok(())
            }
            Descriptor::R
            | Descriptor::S
            | Descriptor::r
            | Descriptor::s
            | Descriptor::M
            | Descriptor::P
            | Descriptor::m
            | Descriptor::p
            | Descriptor::SP_4
            | Descriptor::TBPY_5
            | Descriptor::OC_6 => Err(CipLabelerError::DescriptorNotSupportedForDoubleBonds),
            Descriptor::None | Descriptor::Unknown | Descriptor::ns => {
                Err(CipLabelerError::InvalidBondDescriptor)
            }
        }
    }

    // BEGIN RDKIT CPP FUNCTION Sp2Bond::hasPrimaryLabel (configs/Sp2Bond.cpp)
    // RDKit✔️✔️: bool Sp2Bond::hasPrimaryLabel() const {
    // RDKit✔️✔️:   return dp_bond->hasProp(common_properties::_CIPCode);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Sp2Bond::hasPrimaryLabel
    pub(crate) fn has_primary_label(&self) -> bool {
        self.primary_label.is_some() || self.source_primary_label_visible
    }

    // BEGIN RDKIT CPP FUNCTION Sp2Bond::resetPrimaryLabel (configs/Sp2Bond.cpp)
    // RDKit✔️✔️: void Sp2Bond::resetPrimaryLabel() const {
    // RDKit✔️✔️:   dp_bond->clearProp(common_properties::_CIPCode);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Sp2Bond::resetPrimaryLabel
    pub(crate) fn reset_primary_label(&mut self) {
        self.primary_label = None;
        self.source_primary_label_visible = false;
    }

    // BEGIN RDKIT CPP FUNCTION Sp2Bond::label(const Rules &) (configs/Sp2Bond.cpp)
    // RDKit✔️✔️: Descriptor Sp2Bond::label(const Rules &comp) {
    // RDKit✔️✔️:   auto &digraph = getDigraph();
    // RDKit✔️✔️:   auto root1 = digraph.getOriginalRoot();
    // RDKit✔️✔️:   if (digraph.getCurrentRoot() != root1) {
    // RDKit✔️✔️:     digraph.changeRoot(root1);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return label(root1, digraph, comp);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Sp2Bond::label(const Rules &)
    pub(crate) fn label(
        &mut self,
        comp: &CipRules,
        context: &mut CipLabelerContext,
    ) -> Result<Descriptor, CipLabelerError> {
        let root1 = self.configuration.digraph.get_original_root();
        if self.configuration.digraph.get_current_root() != root1 {
            self.configuration.digraph.change_root(root1)?;
        }
        self.label_node(root1, comp, context)
    }

    pub(crate) fn label_with_external_digraph(
        &mut self,
        root1: CipNodeId,
        digraph: &mut CipDigraph<'_>,
        comp: &CipRules,
        context: &mut CipLabelerContext,
    ) -> Result<Descriptor, CipLabelerError> {
        let focus1 = self.configuration.foci[0];
        let focus2 = self.configuration.foci[1];
        let carriers = self.configuration.get_carriers().to_vec();
        Self::label_node_impl(
            focus1,
            focus2,
            carriers,
            self.cfg,
            &mut self.ranked_anchors,
            root1,
            digraph,
            comp,
            context,
        )
    }

    // BEGIN RDKIT CPP FUNCTION Sp2Bond::label(Node *, Digraph &, const Rules &) (configs/Sp2Bond.cpp)
    // RDKit✔️✔️: Descriptor Sp2Bond::label(Node *root1, Digraph &digraph, const Rules &comp) {
    // RDKit✔️✔️:   const auto &focus1 = getFoci()[0];
    // RDKit✔️✔️:   const auto &focus2 = getFoci()[1];
    // RDKit✔️✔️:
    // RDKit✔️✔️:   d_ranked_anchors.clear();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   const auto &internal = findInternalEdge(root1->getEdges(), focus1, focus2);
    // RDKit✔️✔️:   if (internal == nullptr) {
    // RDKit✔️✔️:     return Descriptor::UNKNOWN;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   const auto &root2 = internal->getOther(root1);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto edges1 = root1->getEdges();
    // RDKit✔️✔️:   auto edges2 = root2->getEdges();
    // RDKit✔️✔️:   removeInternalEdges(edges1, focus1, focus2);
    // RDKit✔️✔️:   removeInternalEdges(edges2, focus1, focus2);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto carriers = std::vector<Atom *>(getCarriers());
    // RDKit✔️✔️:   auto config = d_cfg;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (root1->getAtom() == focus2) {
    // RDKit✔️✔️:     std::swap(carriers[0], carriers[1]);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   digraph.changeRoot(root1);
    // RDKit✔️✔️:   const auto &priority1 = comp.sort(root1, edges1);
    // RDKit✔️✔️:   if (!priority1.isUnique()) {
    // RDKit✔️✔️:     return Descriptor::UNKNOWN;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // swap
    // RDKit✔️✔️:   if (edges1.size() > 1 && carriers[0] != edges1[0]->getEnd()->getAtom()) {
    // RDKit✔️✔️:     if (config == Bond::STEREOCIS) {
    // RDKit✔️✔️:       config = Bond::STEREOTRANS;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       config = Bond::STEREOCIS;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   digraph.changeRoot(root2);
    // RDKit✔️✔️:   const auto &priority2 = comp.sort(root2, edges2);
    // RDKit✔️✔️:   if (!priority2.isUnique()) {
    // RDKit✔️✔️:     return Descriptor::UNKNOWN;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // swap
    // RDKit✔️✔️:   if (edges2.size() > 1 && carriers[1] != edges2[0]->getEnd()->getAtom()) {
    // RDKit✔️✔️:     if (config == Bond::STEREOCIS) {
    // RDKit✔️✔️:       config = Bond::STEREOTRANS;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       config = Bond::STEREOCIS;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   {
    // RDKit✔️✔️:     // At this point, edges1 and edges2 are sorted by priority starting from
    // RDKit✔️✔️:     // this node. Record that now! - they may be resorted after processing
    // RDKit✔️✔️:     // other nodes.
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // As weird as it seems, these may actually be implicit Hs: Rule 2
    // RDKit✔️✔️:     // in the paper on which this code is based states that,
    // RDKit✔️✔️:     // in CIP ranks, H > 1H, so implicit H actually has a higher
    // RDKit✔️✔️:     // priority than 1H (!!!). getAtomIdx() returns Atom::NOATOM
    // RDKit✔️✔️:     // if that is the case.
    // RDKit✔️✔️:     auto carrier1_idx = edges1[0]->getEnd()->getAtomIdx();
    // RDKit✔️✔️:     auto carrier2_idx = edges2[0]->getEnd()->getAtomIdx();
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // Make sure the stereo atoms are in the right order
    // RDKit✔️✔️:     if (edges1[0]->getBeg()->getAtom() == focus1) {
    // RDKit✔️✔️:       d_ranked_anchors.assign({carrier1_idx, carrier2_idx});
    // RDKit✔️✔️:     } else if (edges2[0]->getBeg()->getAtom() == focus1) {
    // RDKit✔️✔️:       d_ranked_anchors.assign({carrier2_idx, carrier1_idx});
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (config == Bond::STEREOCIS) {
    // RDKit✔️✔️:     if (priority1.isPseudoAsymetric() != priority2.isPseudoAsymetric()) {
    // RDKit✔️✔️:       return Descriptor::seqCis;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       return Descriptor::Z;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (config == Bond::STEREOTRANS) {
    // RDKit✔️✔️:     if (priority1.isPseudoAsymetric() != priority2.isPseudoAsymetric()) {
    // RDKit✔️✔️:       return Descriptor::seqTrans;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       return Descriptor::E;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return Descriptor::UNKNOWN;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Sp2Bond::label(Node *, Digraph &, const Rules &)
    fn label_node(
        &mut self,
        root1: CipNodeId,
        comp: &CipRules,
        context: &mut CipLabelerContext,
    ) -> Result<Descriptor, CipLabelerError> {
        let focus1 = self.configuration.foci[0];
        let focus2 = self.configuration.foci[1];
        let carriers = self.configuration.get_carriers().to_vec();
        Self::label_node_impl(
            focus1,
            focus2,
            carriers,
            self.cfg,
            &mut self.ranked_anchors,
            root1,
            &mut self.configuration.digraph,
            comp,
            context,
        )
    }

    fn label_node_impl(
        focus1: usize,
        focus2: usize,
        mut carriers: Vec<Option<usize>>,
        cfg: BondStereo,
        ranked_anchors: &mut Vec<CipSourceIndex>,
        root1: CipNodeId,
        digraph: &mut CipDigraph<'_>,
        comp: &CipRules,
        context: &mut CipLabelerContext,
    ) -> Result<Descriptor, CipLabelerError> {
        ranked_anchors.clear();

        let root1_edges = digraph.node_edges(root1)?;
        let Some(internal) =
            CipConfiguration::find_internal_edge(digraph, &root1_edges, focus1, focus2)
        else {
            return Ok(Descriptor::Unknown);
        };
        let root2 = digraph.edge(internal).get_other(internal, root1)?;

        let mut edges1 = digraph.node_edges(root1)?;
        let mut edges2 = digraph.node_edges(root2)?;
        CipConfiguration::remove_internal_edges(digraph, &mut edges1, focus1, focus2);
        CipConfiguration::remove_internal_edges(digraph, &mut edges2, focus1, focus2);

        if edges1.is_empty() || edges2.is_empty() {
            return Ok(Descriptor::Unknown);
        }

        let mut config = cfg;
        if digraph.node(root1).atom_idx() == Some(focus2) {
            carriers.swap(0, 1);
        }

        digraph.change_root(root1)?;
        let priority1 = comp.sort(digraph, context, root1, &mut edges1, true)?;
        if !priority1.is_unique() {
            return Ok(Descriptor::Unknown);
        }
        if edges1.len() > 1
            && carriers[0] != digraph.node(digraph.edge(edges1[0]).get_end()).atom_idx()
        {
            config = match config {
                BondStereo::Cis => BondStereo::Trans,
                BondStereo::Trans => BondStereo::Cis,
                _ => config,
            };
        }

        digraph.change_root(root2)?;
        let priority2 = comp.sort(digraph, context, root2, &mut edges2, true)?;
        if !priority2.is_unique() {
            return Ok(Descriptor::Unknown);
        }
        if edges2.len() > 1
            && carriers[1] != digraph.node(digraph.edge(edges2[0]).get_end()).atom_idx()
        {
            config = match config {
                BondStereo::Cis => BondStereo::Trans,
                BondStereo::Trans => BondStereo::Cis,
                _ => config,
            };
        }

        let carrier1_idx = digraph
            .node(digraph.edge(edges1[0]).get_end())
            .get_atom_idx()?;
        let carrier2_idx = digraph
            .node(digraph.edge(edges2[0]).get_end())
            .get_atom_idx()?;
        if digraph.node(digraph.edge(edges1[0]).get_beg()).atom_idx() == Some(focus1) {
            ranked_anchors.extend([carrier1_idx, carrier2_idx]);
        } else if digraph.node(digraph.edge(edges2[0]).get_beg()).atom_idx() == Some(focus1) {
            ranked_anchors.extend([carrier2_idx, carrier1_idx]);
        }

        if config == BondStereo::Cis {
            if priority1.is_pseudo_asymetric() != priority2.is_pseudo_asymetric() {
                Ok(Descriptor::seqCis)
            } else {
                Ok(Descriptor::Z)
            }
        } else if config == BondStereo::Trans {
            if priority1.is_pseudo_asymetric() != priority2.is_pseudo_asymetric() {
                Ok(Descriptor::seqTrans)
            } else {
                Ok(Descriptor::E)
            }
        } else {
            Ok(Descriptor::Unknown)
        }
    }
}

pub(crate) struct CipAtropisomerBond<'a> {
    configuration: CipConfiguration<'a>,
    bond_idx: usize,
    cfg: BondStereo,
    ranked_anchors: Vec<CipSourceIndex>,
    primary_label: Option<CipAtropisomerBondPrimaryLabel>,
    source_primary_label_visible: bool,
}

impl<'a> CipAtropisomerBond<'a> {
    // BEGIN RDKIT CPP FUNCTION AtropisomerBond::AtropisomerBond (configs/AtropisomerBond.cpp)
    // RDKit✔️✔️: AtropisomerBond::AtropisomerBond(const CIPMol &mol, Bond *bond, Atom *startAtom,
    // RDKit✔️✔️:                                  Atom *endAtom, Bond::BondStereo cfg)
    // RDKit✔️✔️:     : Configuration(mol, {startAtom, endAtom}, true),
    // RDKit✔️✔️:       dp_bond{bond},
    // RDKit✔️✔️:       d_cfg{cfg} {
    // RDKit✔️✔️:   CHECK_INVARIANT(startAtom && endAtom, "bad foci")
    // RDKit✔️✔️:   CHECK_INVARIANT(d_cfg == Bond::STEREOATROPCW || d_cfg == Bond::STEREOATROPCCW,
    // RDKit✔️✔️:                   "bad config")
    // RDKit✔️✔️:
    // RDKit✔️✔️:   Atropisomers::AtropAtomAndBondVec atomAndBondVecs[2];
    // RDKit✔️✔️:   if (!Atropisomers::getAtropisomerAtomsAndBonds(bond, atomAndBondVecs,
    // RDKit✔️✔️:                                                  bond->getOwningMol())) {
    // RDKit✔️✔️:     return;  // not an atropisomer
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto atom1 = mol.getAtom(atomAndBondVecs[0].second[0]->getOtherAtomIdx(
    // RDKit✔️✔️:       atomAndBondVecs[0].first->getIdx()));
    // RDKit✔️✔️:   auto atom2 = mol.getAtom(atomAndBondVecs[1].second[0]->getOtherAtomIdx(
    // RDKit✔️✔️:       atomAndBondVecs[1].first->getIdx()));
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::vector<Atom *> anchors{atom1, atom2};
    // RDKit✔️✔️:
    // RDKit✔️✔️:   setCarriers(std::move(anchors));
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION AtropisomerBond::AtropisomerBond
    pub(crate) fn new(
        molecule: &'a Molecule,
        bond_idx: usize,
        start_atom: usize,
        end_atom: usize,
        cfg: BondStereo,
    ) -> Result<Self, CipLabelerError> {
        let mut configuration =
            CipConfiguration::with_foci(molecule, vec![start_atom, end_atom], true)?;
        configuration.digraph.mol().atom(start_atom)?;
        configuration.digraph.mol().atom(end_atom)?;
        let bond = configuration.digraph.mol().bond(bond_idx)?;
        let source_primary_label_visible = bond.prop("_CIPCode").is_some();
        if !((bond.begin().index() == start_atom && bond.end().index() == end_atom)
            || (bond.begin().index() == end_atom && bond.end().index() == start_atom))
        {
            return Err(CipLabelerError::BadAtropisomerBondFoci);
        }
        if !matches!(cfg, BondStereo::AtropCw | BondStereo::AtropCcw) {
            return Err(CipLabelerError::BadAtropisomerBondConfig);
        }

        if let Some(parts) = atropisomer_atoms_and_bonds(molecule, BondId::new(bond_idx))
            && let Some(carriers) = parts.first_neighbor_atoms(molecule)
        {
            configuration.set_carriers(vec![Some(carriers[0].index()), Some(carriers[1].index())]);
        }

        Ok(Self {
            configuration,
            bond_idx,
            cfg,
            ranked_anchors: Vec::new(),
            primary_label: None,
            source_primary_label_visible,
        })
    }

    pub(crate) fn get_foci(&self) -> &[usize] {
        self.configuration.get_foci()
    }

    pub(crate) fn get_carriers(&self) -> &[Option<usize>] {
        self.configuration.get_carriers()
    }

    pub(crate) fn ranked_anchors(&self) -> &[CipSourceIndex] {
        &self.ranked_anchors
    }

    pub(crate) fn primary_label(&self) -> Option<&CipAtropisomerBondPrimaryLabel> {
        self.primary_label.as_ref()
    }

    // BEGIN RDKIT CPP FUNCTION AtropisomerBond::setPrimaryLabel (configs/AtropisomerBond.cpp)
    // RDKit✔️✔️: void AtropisomerBond::setPrimaryLabel(Descriptor desc) {
    // RDKit✔️✔️:   switch (desc) {
    // RDKit✔️✔️:     case Descriptor::M:
    // RDKit✔️✔️:     case Descriptor::P:
    // RDKit✔️✔️:     case Descriptor::m:
    // RDKit✔️✔️:     case Descriptor::p: {
    // RDKit✔️✔️:       dp_bond->setProp(common_properties::_CIPCode, to_string(desc));
    // RDKit✔️✔️:       dp_bond->setProp(common_properties::_CIPNeighborOrder, d_ranked_anchors,
    // RDKit✔️✔️:                        true);
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     case Descriptor::R:
    // RDKit✔️✔️:     case Descriptor::S:
    // RDKit✔️✔️:     case Descriptor::r:
    // RDKit✔️✔️:     case Descriptor::s:
    // RDKit✔️✔️:     case Descriptor::SP_4:
    // RDKit✔️✔️:     case Descriptor::TBPY_5:
    // RDKit✔️✔️:     case Descriptor::OC_6:
    // RDKit✔️✔️:     case Descriptor::seqTrans:
    // RDKit✔️✔️:     case Descriptor::E:
    // RDKit✔️✔️:     case Descriptor::seqCis:
    // RDKit✔️✔️:     case Descriptor::Z:
    // RDKit✔️✔️:       throw std::runtime_error(
    // RDKit✔️✔️:           "Received a Descriptor that is not supported for atropisomer bonds");
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       throw std::runtime_error("Received an invalid Bond Descriptor");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION AtropisomerBond::setPrimaryLabel
    pub(crate) fn set_primary_label(&mut self, desc: Descriptor) -> Result<(), CipLabelerError> {
        match desc {
            Descriptor::M | Descriptor::P | Descriptor::m | Descriptor::p => {
                self.primary_label = Some(CipAtropisomerBondPrimaryLabel {
                    bond_idx: self.bond_idx,
                    cip_code: descriptor_to_string(desc),
                    cip_neighbor_order: self.ranked_anchors.clone(),
                });
                self.source_primary_label_visible = false;
                Ok(())
            }
            Descriptor::R
            | Descriptor::S
            | Descriptor::r
            | Descriptor::s
            | Descriptor::SP_4
            | Descriptor::TBPY_5
            | Descriptor::OC_6
            | Descriptor::seqTrans
            | Descriptor::E
            | Descriptor::seqCis
            | Descriptor::Z => Err(CipLabelerError::DescriptorNotSupportedForAtropisomerBonds),
            Descriptor::None | Descriptor::Unknown | Descriptor::ns => {
                Err(CipLabelerError::InvalidBondDescriptor)
            }
        }
    }

    // BEGIN RDKIT CPP FUNCTION AtropisomerBond::hasPrimaryLabel/resetPrimaryLabel (configs/AtropisomerBond.cpp)
    // RDKit✔️✔️: bool AtropisomerBond::hasPrimaryLabel() const {
    // RDKit✔️✔️:   return dp_bond->hasProp(common_properties::_CIPCode);
    // RDKit✔️✔️: }
    // RDKit✔️✔️:
    // RDKit✔️✔️: void AtropisomerBond::resetPrimaryLabel() const {
    // RDKit✔️✔️:   dp_bond->clearProp(common_properties::_CIPCode);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION AtropisomerBond::hasPrimaryLabel/resetPrimaryLabel
    pub(crate) fn has_primary_label(&self) -> bool {
        self.primary_label.is_some() || self.source_primary_label_visible
    }

    pub(crate) fn reset_primary_label(&mut self) {
        self.primary_label = None;
        self.source_primary_label_visible = false;
    }

    // BEGIN RDKIT CPP FUNCTION AtropisomerBond::label(const Rules &) (configs/AtropisomerBond.cpp)
    // RDKit✔️✔️: Descriptor AtropisomerBond::label(const Rules &comp) {
    // RDKit✔️✔️:   auto &digraph = getDigraph();
    // RDKit✔️✔️:   auto root1 = digraph.getOriginalRoot();
    // RDKit✔️✔️:   if (digraph.getCurrentRoot() != root1) {
    // RDKit✔️✔️:     digraph.changeRoot(root1);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return label(root1, digraph, comp);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION AtropisomerBond::label(const Rules &)
    pub(crate) fn label(
        &mut self,
        comp: &CipRules,
        context: &mut CipLabelerContext,
    ) -> Result<Descriptor, CipLabelerError> {
        let root1 = self.configuration.digraph.get_original_root();
        if self.configuration.digraph.get_current_root() != root1 {
            self.configuration.digraph.change_root(root1)?;
        }
        self.label_with_root(root1, comp, context)
    }

    pub(crate) fn label_with_external_digraph(
        &mut self,
        root1: CipNodeId,
        digraph: &mut CipDigraph<'_>,
        comp: &CipRules,
        context: &mut CipLabelerContext,
    ) -> Result<Descriptor, CipLabelerError> {
        let focus1 = self.configuration.foci[0];
        let focus2 = self.configuration.foci[1];
        let carriers = self.configuration.get_carriers().to_vec();
        Self::label_node_impl(
            focus1,
            focus2,
            carriers,
            self.cfg,
            &mut self.ranked_anchors,
            root1,
            digraph,
            comp,
            context,
        )
    }

    // BEGIN RDKIT CPP FUNCTION AtropisomerBond::label(Node *, Digraph &, const Rules &) (configs/AtropisomerBond.cpp)
    // RDKit✔️✔️: Descriptor AtropisomerBond::label(Node *root1, Digraph &digraph,
    // RDKit✔️✔️:                                   const Rules &comp) {
    // RDKit✔️✔️:   const auto &focus1 = getFoci()[0];
    // RDKit✔️✔️:   const auto &focus2 = getFoci()[1];
    // RDKit✔️✔️:
    // RDKit✔️✔️:   d_ranked_anchors.clear();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   const auto &internal = findInternalEdge(root1->getEdges(), focus1, focus2);
    // RDKit✔️✔️:   if (internal == nullptr) {
    // RDKit✔️✔️:     return Descriptor::UNKNOWN;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   const auto &root2 = internal->getOther(root1);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto edges1 = root1->getEdges();
    // RDKit✔️✔️:   auto edges2 = root2->getEdges();
    // RDKit✔️✔️:   removeInternalEdges(edges1, focus1, focus2);
    // RDKit✔️✔️:   removeInternalEdges(edges2, focus1, focus2);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   removeDuplicatesAndHs(edges1);
    // RDKit✔️✔️:   removeDuplicatesAndHs(edges2);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto carriers = std::vector<Atom *>(getCarriers());
    // RDKit✔️✔️:   auto config = d_cfg;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (root1->getAtom() == focus2) {
    // RDKit✔️✔️:     std::swap(carriers[0], carriers[1]);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   digraph.changeRoot(root1);
    // RDKit✔️✔️:   const auto &priority1 = comp.sort(root1, edges1);
    // RDKit✔️✔️:   if (!priority1.isUnique()) {
    // RDKit✔️✔️:     return Descriptor::UNKNOWN;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // swap
    // RDKit✔️✔️:   if (edges1.size() > 1 && carriers[0] == edges1[1]->getEnd()->getAtom()) {
    // RDKit✔️✔️:     if (config == Bond::STEREOATROPCCW) {
    // RDKit✔️✔️:       config = Bond::STEREOATROPCW;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       config = Bond::STEREOATROPCCW;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   digraph.changeRoot(root2);
    // RDKit✔️✔️:   const auto &priority2 = comp.sort(root2, edges2);
    // RDKit✔️✔️:   if (!priority2.isUnique()) {
    // RDKit✔️✔️:     return Descriptor::UNKNOWN;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // swap
    // RDKit✔️✔️:   if (edges2.size() > 1 && carriers[1] == edges2[1]->getEnd()->getAtom()) {
    // RDKit✔️✔️:     if (config == Bond::STEREOATROPCCW) {
    // RDKit✔️✔️:       config = Bond::STEREOATROPCW;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       config = Bond::STEREOATROPCCW;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   {
    // RDKit✔️✔️:     // This is mostly the same as in Sp2Bonds, but I doubt the anchors will be
    // RDKit✔️✔️:     // implicit Hs in this case.
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // At this point, edges1 and edges2 are sorted by priority starting from
    // RDKit✔️✔️:     // this node. Record that now! - they may be resorted after processing
    // RDKit✔️✔️:     // other nodes.
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // As weird as it seems, these may actually be implicit Hs: Rule 2
    // RDKit✔️✔️:     // in the paper on which this code is based states that,
    // RDKit✔️✔️:     // in CIP ranks, H > 1H, so implicit H actually has a higher
    // RDKit✔️✔️:     // priority than 1H (!!!). getAtomIdx() returns Atom::NOATOM
    // RDKit✔️✔️:     // if that is the case.
    // RDKit✔️✔️:     auto carrier1_idx = edges1[0]->getEnd()->getAtomIdx();
    // RDKit✔️✔️:     auto carrier2_idx = edges2[0]->getEnd()->getAtomIdx();
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // Make sure the stereo atoms are in the right order
    // RDKit✔️✔️:     if (edges1[0]->getBeg()->getAtom() == focus1) {
    // RDKit✔️✔️:       d_ranked_anchors.assign({carrier1_idx, carrier2_idx});
    // RDKit✔️✔️:     } else if (edges2[0]->getBeg()->getAtom() == focus1) {
    // RDKit✔️✔️:       d_ranked_anchors.assign({carrier2_idx, carrier1_idx});
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (config == Bond::STEREOATROPCCW) {
    // RDKit✔️✔️:     if (priority1.isPseudoAsymetric() || priority2.isPseudoAsymetric()) {
    // RDKit✔️✔️:       return Descriptor::m;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       return Descriptor::M;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (config == Bond::STEREOATROPCW) {
    // RDKit✔️✔️:     if (priority1.isPseudoAsymetric() || priority2.isPseudoAsymetric()) {
    // RDKit✔️✔️:       return Descriptor::p;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       return Descriptor::P;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return Descriptor::UNKNOWN;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION AtropisomerBond::label(Node *, Digraph &, const Rules &)
    fn label_with_root(
        &mut self,
        root1: CipNodeId,
        comp: &CipRules,
        context: &mut CipLabelerContext,
    ) -> Result<Descriptor, CipLabelerError> {
        let focus1 = self.configuration.foci[0];
        let focus2 = self.configuration.foci[1];
        let carriers = self.configuration.get_carriers().to_vec();
        Self::label_node_impl(
            focus1,
            focus2,
            carriers,
            self.cfg,
            &mut self.ranked_anchors,
            root1,
            &mut self.configuration.digraph,
            comp,
            context,
        )
    }

    fn label_node_impl(
        focus1: usize,
        focus2: usize,
        mut carriers: Vec<Option<usize>>,
        cfg: BondStereo,
        ranked_anchors: &mut Vec<CipSourceIndex>,
        root1: CipNodeId,
        digraph: &mut CipDigraph<'_>,
        comp: &CipRules,
        context: &mut CipLabelerContext,
    ) -> Result<Descriptor, CipLabelerError> {
        ranked_anchors.clear();
        if carriers.len() < 2 {
            return Ok(Descriptor::Unknown);
        }

        let root1_edges = digraph.node_edges(root1)?;
        let Some(internal) =
            CipConfiguration::find_internal_edge(digraph, &root1_edges, focus1, focus2)
        else {
            return Ok(Descriptor::Unknown);
        };
        let root2 = digraph.edge(internal).get_other(internal, root1)?;

        let mut edges1 = digraph.node_edges(root1)?;
        let mut edges2 = digraph.node_edges(root2)?;
        CipConfiguration::remove_internal_edges(digraph, &mut edges1, focus1, focus2);
        CipConfiguration::remove_internal_edges(digraph, &mut edges2, focus1, focus2);
        CipConfiguration::remove_duplicates_and_hs(digraph, &mut edges1);
        CipConfiguration::remove_duplicates_and_hs(digraph, &mut edges2);
        if edges1.is_empty() || edges2.is_empty() {
            return Ok(Descriptor::Unknown);
        }

        let mut config = cfg;
        if digraph.node(root1).atom_idx() == Some(focus2) {
            carriers.swap(0, 1);
        }

        digraph.change_root(root1)?;
        let priority1 = comp.sort(digraph, context, root1, &mut edges1, true)?;
        if !priority1.is_unique() {
            return Ok(Descriptor::Unknown);
        }
        if edges1.len() > 1
            && carriers[0] == digraph.node(digraph.edge(edges1[1]).get_end()).atom_idx()
        {
            config = if config == BondStereo::AtropCcw {
                BondStereo::AtropCw
            } else {
                BondStereo::AtropCcw
            };
        }

        digraph.change_root(root2)?;
        let priority2 = comp.sort(digraph, context, root2, &mut edges2, true)?;
        if !priority2.is_unique() {
            return Ok(Descriptor::Unknown);
        }
        if edges2.len() > 1
            && carriers[1] == digraph.node(digraph.edge(edges2[1]).get_end()).atom_idx()
        {
            config = if config == BondStereo::AtropCcw {
                BondStereo::AtropCw
            } else {
                BondStereo::AtropCcw
            };
        }

        let carrier1_idx = digraph
            .node(digraph.edge(edges1[0]).get_end())
            .get_atom_idx()?;
        let carrier2_idx = digraph
            .node(digraph.edge(edges2[0]).get_end())
            .get_atom_idx()?;
        if digraph.node(digraph.edge(edges1[0]).get_beg()).atom_idx() == Some(focus1) {
            ranked_anchors.extend([carrier1_idx, carrier2_idx]);
        } else if digraph.node(digraph.edge(edges2[0]).get_beg()).atom_idx() == Some(focus1) {
            ranked_anchors.extend([carrier2_idx, carrier1_idx]);
        }

        if config == BondStereo::AtropCcw {
            if priority1.is_pseudo_asymetric() || priority2.is_pseudo_asymetric() {
                Ok(Descriptor::m)
            } else {
                Ok(Descriptor::M)
            }
        } else if config == BondStereo::AtropCw {
            if priority1.is_pseudo_asymetric() || priority2.is_pseudo_asymetric() {
                Ok(Descriptor::p)
            } else {
                Ok(Descriptor::P)
            }
        } else {
            Ok(Descriptor::Unknown)
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
            (0, rdkit_atomic_mass(atomic_number, isotope))
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
        molecule: &'a Molecule,
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
    // RDKit✔️✔️: const std::vector<Edge *> &Node::getEdges() const {
    // RDKit✔️✔️:   if (!isExpanded()) {
    // RDKit✔️✔️:     auto non_const_this = const_cast<Node *>(this);
    // RDKit✔️✔️:     non_const_this->d_flags |= EXPANDED;
    // RDKit✔️✔️:     dp_g->expand(non_const_this);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return d_edges;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Node::getEdges
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
    // RDKit✔️✔️: std::vector<Node *> Digraph::getNodes(Atom *atom) const {
    // RDKit✔️✔️:   std::vector<Node *> result;
    // RDKit✔️✔️:   std::vector<Node*> queue = {getCurrentRoot()};
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (size_t i=0; i<queue.size(); ++i) {
    // RDKit✔️✔️:     auto node = queue[i];
    // RDKit✔️✔️:     if (atom == node->getAtom()) {
    // RDKit✔️✔️:       result.push_back(node);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (const auto &e : node->getEdges()) {
    // RDKit✔️✔️:       if (!e->isBeg(node)) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       queue.push_back(e->getEnd());
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return result;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Digraph::getNodes
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
    // RDKit✔️✔️: void Digraph::changeRoot(Node *newroot) {
    // RDKit✔️✔️:   std::vector<Edge *> toflip;
    // RDKit✔️✔️:   auto queue = std::list<Node *>({newroot});
    // RDKit✔️✔️:   for (const auto &node : queue) {
    // RDKit✔️✔️:     for (const auto &e : node->getEdges()) {
    // RDKit✔️✔️:       if (e->isEnd(node)) {
    // RDKit✔️✔️:         toflip.push_back(e);
    // RDKit✔️✔️:         queue.push_back(e->getBeg());
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (auto &e : toflip) {
    // RDKit✔️✔️:     e->flip();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   dp_root = newroot;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Digraph::changeRoot
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
    // RDKit✔️✔️: void Digraph::expand(Node *beg) {
    // RDKit✔️✔️:   const auto &atom = beg->getAtom();
    // RDKit✔️✔️:   const auto &edges = beg->getEdges();
    // RDKit✔️✔️:   const auto &prev =
    // RDKit✔️✔️:       edges.size() > 0 && !edges[0]->isBeg(beg) ? edges[0]->getBond() : nullptr;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (MAX_NODE_DIST > 0 && beg->getDistance() > MAX_NODE_DIST) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (MAX_NODE_COUNT > 0 && d_nodes.size() >= MAX_NODE_COUNT) {
    // RDKit✔️✔️:     std::stringstream errmsg;
    // RDKit✔️✔️:     errmsg << "Digraph generation failed: more than " << MAX_NODE_COUNT
    // RDKit✔️✔️:            << "nodes found.";
    // RDKit✔️✔️:     throw TooManyNodesException(errmsg.str());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // create 'explicit' nodes
    // RDKit✔️✔️:   for (const auto &bond : d_mol.getBonds(atom)) {
    // RDKit✔️✔️:     const auto &nbr = bond->getOtherAtom(atom);
    // RDKit✔️✔️:     const int nbrIdx = nbr->getIdx();
    // RDKit✔️✔️:     const int bord = d_mol.getBondOrder(bond);
    // RDKit✔️✔️:     const int virtual_nodes = bord - 1;
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (!beg->isVisited(nbrIdx)) {
    // RDKit✔️✔️:       auto end = beg->newChild(nbrIdx, nbr);
    // RDKit✔️✔️:       addEdge(beg, bond, end);
    // RDKit✔️✔️:
    // RDKit✔️✔️:       // duplicate nodes for bond orders (except for root atoms...)
    // RDKit✔️✔️:       // for example >S=O
    // RDKit✔️✔️:       if (dp_origin != beg || d_atropisomerMode) {
    // RDKit✔️✔️:         if (atom->getFormalCharge() < 0 &&
    // RDKit✔️✔️:             d_mol.getFractionalAtomicNum(atom).denominator() > 1) {
    // RDKit✔️✔️:           end = beg->newBondDuplicateChild(nbrIdx, nbr);
    // RDKit✔️✔️:           addEdge(beg, bond, end);
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           for (int i = 0; i < virtual_nodes; ++i) {
    // RDKit✔️✔️:             end = beg->newBondDuplicateChild(nbrIdx, nbr);
    // RDKit✔️✔️:             addEdge(beg, bond, end);
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (bond == prev) {  // bond order expansion (backwards)
    // RDKit✔️✔️:       if (dp_origin->getAtom() != nbr || d_atropisomerMode) {
    // RDKit✔️✔️:         for (int i = 0; i < virtual_nodes; ++i) {
    // RDKit✔️✔️:           auto end = beg->newBondDuplicateChild(nbrIdx, nbr);
    // RDKit✔️✔️:           addEdge(beg, bond, end);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {  // ring closures
    // RDKit✔️✔️:       auto end = beg->newRingDuplicateChild(nbrIdx, nbr);
    // RDKit✔️✔️:       addEdge(beg, bond, end);
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (atom->getFormalCharge() < 0 &&
    // RDKit✔️✔️:           d_mol.getFractionalAtomicNum(atom).denominator() > 1) {
    // RDKit✔️✔️:         end = beg->newBondDuplicateChild(nbrIdx, nbr);
    // RDKit✔️✔️:         addEdge(beg, bond, end);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         for (int i = 0; i < virtual_nodes; ++i) {
    // RDKit✔️✔️:           end = beg->newBondDuplicateChild(nbrIdx, nbr);
    // RDKit✔️✔️:           addEdge(beg, bond, end);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // Create implicit hydrogen nodes
    // RDKit✔️✔️:   const int hcnt = atom->getTotalNumHs();
    // RDKit✔️✔️:   for (int i = 0; i < hcnt; ++i) {
    // RDKit✔️✔️:     auto end = beg->newImplicitHydrogenChild();
    // RDKit✔️✔️:     addEdge(beg, nullptr, end);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Digraph::expand
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
    pub(crate) fn new(molecule: &'a Molecule) -> Self {
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
    // RDKit✔️✔️: CIPMolSpan<Bond *, ROMol::OEDGE_ITER> CIPMol::getBonds(Atom *atom) const {
    // RDKit✔️✔️:   PRECONDITION(atom, "bad atom")
    // RDKit✔️✔️:   return {d_mol, d_mol.getAtomBonds(atom)};
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION CIPMol::getBonds
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
    // RDKit✔️✔️: CIPMolSpan<Atom *, ROMol::ADJ_ITER> CIPMol::getNeighbors(Atom *atom) const {
    // RDKit✔️✔️:   PRECONDITION(atom, "bad atom")
    // RDKit✔️✔️:   return {d_mol, d_mol.getAtomNeighbors(atom)};
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION CIPMol::getNeighbors
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
            self.rings = Some(crate::rings::fast_find_rings(self.molecule)?);
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
            if let Ok(assignment) = crate::kekulize::kekulize_assignment(
                self.molecule,
                self.rings.as_ref(),
                true,
                true,
                100,
            ) {
                for (idx, order) in orders.iter_mut().enumerate() {
                    if let Some(kekulized) = assignment.bond_order(BondId::new(idx)) {
                        *order = kekulized;
                    }
                }
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
            self.valence = Some(crate::valence::assign_valence(
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

#[cfg(test)]
mod tests {
    use std::{cell::Cell, rc::Rc};

    use super::{
        CipAtropisomerBond, CipCancellationMode, CipConfiguration, CipDigraph, CipEdge, CipEdgeId,
        CipLabelerContext, CipLabelerError, CipMol, CipNode, CipNodeId, CipPairList, CipPriority,
        CipRule1a, CipRule1b, CipRule2, CipRule3, CipRule4a, CipRule4b, CipRule4c, CipRule5New,
        CipRule6, CipRules, CipSequenceRule, CipSort, CipSp2Bond, CipTetrahedral, Descriptor,
        MancudeType, RationalI32, assign_cip_labels, assign_cip_labels_for_indices,
        assign_cip_labels_for_masks, assign_cip_labels_for_masks_with_cancellation, cip_all_rules,
        descriptor_to_string, relax_types, seed_types, three_way_comparison_i32, visit_parts,
    };
    use crate::{
        AtomSpec, BondId, BondOrder, BondSpec, BondStereo, ChiralTag, Element, Molecule,
        MoleculeBuilder,
    };

    #[test]
    fn ciplabeler_rational_order_matches_boost_rational() {
        assert!(RationalI32::new(13, 2) < RationalI32::new(7, 1));
        assert!(RationalI32::new(-3, 2) < RationalI32::new(-1, 3));
        assert_eq!(RationalI32::new(6, 4), RationalI32::new(3, 2));
    }

    #[test]
    fn ciplabeler_chembl_recursive_priority_regressions_match_rdkit() {
        let molecule = Molecule::from_smiles("CCOC(C)(C)C(=O)N[C@@H](COC)c1ccccn1").unwrap();
        let atom_mask = vec![true; molecule.num_atoms()];
        let bond_mask = vec![true; molecule.num_bonds()];
        let labeled = assign_cip_labels_for_masks(&molecule, &atom_mask, &bond_mask, 0).unwrap();
        assert_eq!(labeled.atoms()[9].prop("_CIPCode"), Some("R"));
        assert_eq!(
            labeled.atoms()[9].prop("_CIPNeighborOrder"),
            Some("[8,10,13]")
        );
    }

    #[test]
    fn ciplabeler_sort_priority_descriptor_to_string_matches_rdkit() {
        let expected = [
            (Descriptor::None, "NONE"),
            (Descriptor::Unknown, "UNKNOWN"),
            (Descriptor::ns, "ns"),
            (Descriptor::R, "R"),
            (Descriptor::S, "S"),
            (Descriptor::r, "r"),
            (Descriptor::s, "s"),
            (Descriptor::seqTrans, "e"),
            (Descriptor::seqCis, "z"),
            (Descriptor::E, "E"),
            (Descriptor::Z, "Z"),
            (Descriptor::M, "M"),
            (Descriptor::P, "P"),
            (Descriptor::m, "m"),
            (Descriptor::p, "p"),
            (Descriptor::SP_4, "SP_4"),
            (Descriptor::TBPY_5, "TBPY_5"),
            (Descriptor::OC_6, "OC_6"),
        ];
        for (desc, label) in expected {
            assert_eq!(descriptor_to_string(desc), label);
        }
    }

    #[test]
    fn ciplabeler_sort_priority_descriptor_order_matches_rdkit_enum_order() {
        assert_eq!(
            Descriptor::ALL_IN_RDKIT_ORDER,
            [
                Descriptor::None,
                Descriptor::Unknown,
                Descriptor::ns,
                Descriptor::R,
                Descriptor::S,
                Descriptor::r,
                Descriptor::s,
                Descriptor::seqTrans,
                Descriptor::seqCis,
                Descriptor::E,
                Descriptor::Z,
                Descriptor::M,
                Descriptor::P,
                Descriptor::m,
                Descriptor::p,
                Descriptor::SP_4,
                Descriptor::TBPY_5,
                Descriptor::OC_6,
            ]
        );
        for window in Descriptor::ALL_IN_RDKIT_ORDER.windows(2) {
            assert!(window[0] < window[1]);
        }
    }

    #[test]
    fn ciplabeler_cipmol_accessors_match_rdkit_molecule_view() {
        let molecule = Molecule::from_smiles("CCO").unwrap();
        let cipmol = CipMol::new(&molecule);
        assert_eq!(cipmol.get_num_atoms(), 3);
        assert_eq!(cipmol.get_num_bonds(), 2);
        assert_eq!(cipmol.atom(0).unwrap().atomic_number(), 6);
        assert_eq!(cipmol.atom(2).unwrap().atomic_number(), 8);
        assert_eq!(cipmol.bond(0).unwrap().order(), BondOrder::Single);
        assert_eq!(cipmol.neighbor_indices(1).unwrap(), vec![0, 2]);
        assert_eq!(cipmol.bond_indices_for_atom(1).unwrap(), vec![0, 1]);
    }

    #[test]
    fn ciplabeler_cipmol_ring_membership_uses_fast_ring_info() {
        let molecule = Molecule::from_smiles("C1CCCCC1").unwrap();
        let mut cipmol = CipMol::new(&molecule);
        assert!(cipmol.is_in_ring(0).unwrap());
        assert!(cipmol.is_in_ring(5).unwrap());

        let chain = Molecule::from_smiles("CCCC").unwrap();
        let mut chain_cipmol = CipMol::new(&chain);
        assert!(!chain_cipmol.is_in_ring(0).unwrap());
    }

    #[test]
    fn ciplabeler_cipmol_bond_order_uses_kekulized_aromatic_copy() {
        let molecule = Molecule::from_smiles("c1ccccc1").unwrap();
        let mut cipmol = CipMol::new(&molecule);
        let mut orders = (0..molecule.num_bonds())
            .map(|idx| cipmol.get_bond_order(idx).unwrap())
            .collect::<Vec<_>>();
        orders.sort_unstable();
        assert_eq!(orders, vec![1, 1, 1, 2, 2, 2]);
    }

    #[test]
    fn ciplabeler_cipmol_bond_order_maps_zero_hydrogen_and_dative_to_zero() {
        let mut builder = MoleculeBuilder::new();
        let a = builder.add_atom(AtomSpec::new(Element::N));
        let b = builder.add_atom(AtomSpec::new(Element::B));
        let c = builder.add_atom(AtomSpec::new(Element::H));
        let d = builder.add_atom(AtomSpec::new(Element::H));
        builder
            .add_bond(BondSpec::new(a, b, BondOrder::Dative))
            .unwrap();
        builder
            .add_bond(BondSpec::new(b, c, BondOrder::Zero))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c, d, BondOrder::Hydrogen))
            .unwrap();
        let molecule = builder.build().unwrap();
        let mut cipmol = CipMol::new(&molecule);
        assert_eq!(cipmol.get_bond_order(0).unwrap(), 0);
        assert_eq!(cipmol.get_bond_order(1).unwrap(), 0);
        assert_eq!(cipmol.get_bond_order(2).unwrap(), 0);
    }

    #[test]
    fn ciplabeler_cipmol_fractional_atomic_numbers_handle_mancude_nitrogen() {
        let molecule = Molecule::from_smiles("c1ccncc1").unwrap();
        let mut cipmol = CipMol::new(&molecule);
        let fractions = (0..molecule.num_atoms())
            .map(|idx| cipmol.get_fractional_atomic_num(idx).unwrap().tuple())
            .collect::<Vec<_>>();
        assert_eq!(
            fractions,
            vec![(6, 1), (6, 1), (13, 2), (6, 1), (13, 2), (6, 1)]
        );
    }

    #[test]
    fn ciplabeler_cipmol_fractional_atomic_numbers_recompute_negative_charge_resonance() {
        let molecule = Molecule::from_smiles("C1=C[CH-]C=C1").unwrap();
        let mut cipmol = CipMol::new(&molecule);
        let fractions = (0..molecule.num_atoms())
            .map(|idx| cipmol.get_fractional_atomic_num(idx).unwrap().tuple())
            .collect::<Vec<_>>();
        assert_eq!(fractions, vec![(0, 1), (6, 1), (6, 1), (4, 1), (9, 2)]);
    }

    #[test]
    fn ciplabeler_cipmol_mancude_seed_relax_and_partition_match_rdkit() {
        let molecule = Molecule::from_smiles("c1ccncc1").unwrap();
        let mut cipmol = CipMol::new(&molecule);
        let mut types = vec![MancudeType::Other; molecule.num_atoms()];

        assert!(seed_types(&mut types, &mut cipmol).unwrap());
        assert_eq!(
            types,
            vec![
                MancudeType::Cv4D3,
                MancudeType::Cv4D3,
                MancudeType::Cv4D3,
                MancudeType::Nv3D2,
                MancudeType::Cv4D3,
                MancudeType::Cv4D3,
            ]
        );
        relax_types(&mut types, &cipmol).unwrap();
        assert!(types.iter().all(|kind| *kind != MancudeType::Other));

        let mut parts = vec![0; molecule.num_atoms()];
        assert_eq!(visit_parts(&mut parts, &types, &mut cipmol).unwrap(), 1);
        assert_eq!(parts, vec![1; molecule.num_atoms()]);
    }

    #[test]
    fn ciplabeler_cipmol_mancude_relax_removes_non_resonant_typed_chain() {
        let molecule = Molecule::from_smiles("CCC").unwrap();
        let cipmol = CipMol::new(&molecule);
        let mut types = vec![MancudeType::Cv4D3, MancudeType::Cv4D3, MancudeType::Other];

        relax_types(&mut types, &cipmol).unwrap();
        assert_eq!(types, vec![MancudeType::Other; 3]);
    }

    #[test]
    fn ciplabeler_cipmol_rejects_non_integer_bond_orders() {
        let mut builder = MoleculeBuilder::new();
        let first = builder.add_atom(AtomSpec::new(Element::C));
        let second = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(first, second, BondOrder::OneAndHalf))
            .unwrap();
        let molecule = builder.build().unwrap();
        let mut cipmol = CipMol::new(&molecule);

        assert!(matches!(
            cipmol.get_bond_order(0),
            Err(CipLabelerError::NonIntegerBondOrder {
                order: BondOrder::OneAndHalf
            })
        ));
    }

    #[test]
    fn ciplabeler_node_constructor_flags_mass_and_visits_match_rdkit() {
        let molecule = Molecule::from_smiles("[13CH4]").unwrap();
        let cipmol = CipMol::new(&molecule);
        let node =
            CipNode::new(11, vec![1], Some(0), RationalI32::new(6, 1), 1, 0, &cipmol).unwrap();

        assert_eq!(node.get_digraph(), 11);
        assert_eq!(node.atom_idx(), Some(0));
        assert_eq!(node.get_atom_idx().unwrap(), 0);
        assert_eq!(node.get_distance(), 1);
        assert_eq!(node.get_atomic_num_fraction().tuple(), (6, 1));
        assert_eq!(node.get_atomic_num(&cipmol).unwrap(), 6);
        assert_eq!(node.get_mass_num(&cipmol).unwrap(), 13);
        assert!(node.get_atomic_mass() > 13.0);
        assert_eq!(node.get_aux(), Descriptor::None);
        assert!(!node.is_duplicate());
        assert!(!node.is_duplicate_or_h());
        assert!(!node.is_expanded());
        assert!(node.is_visited(0));
        assert!(!node.is_terminal());

        let duplicate = CipNode::new(
            11,
            vec![1],
            Some(0),
            RationalI32::new(6, 1),
            1,
            CipNode::BOND_DUPLICATE,
            &cipmol,
        )
        .unwrap();
        assert!(duplicate.is_duplicate());
        assert!(duplicate.is_duplicate_or_h());
        assert!(duplicate.is_expanded());
        assert_eq!(duplicate.get_mass_num(&cipmol).unwrap(), 0);
        assert_eq!(duplicate.get_atomic_mass(), 0.0);
    }

    #[test]
    fn ciplabeler_node_child_creation_and_aux_storage_match_rdkit() {
        let molecule = Molecule::from_smiles("CC").unwrap();
        let mut cipmol = CipMol::new(&molecule);
        let root = CipNode::new(
            3,
            vec![1, 0],
            Some(0),
            RationalI32::new(6, 1),
            1,
            0,
            &cipmol,
        )
        .unwrap();

        let child = root.new_child(1, Some(1), &cipmol).unwrap();
        assert_eq!(child.atom_idx(), Some(1));
        assert_eq!(child.get_distance(), 2);
        assert!(child.is_visited(0));
        assert!(child.is_visited(1));
        assert!(!child.is_expanded());

        let ring_duplicate = root
            .new_ring_duplicate_child(0, Some(0), &mut cipmol)
            .unwrap();
        assert!(ring_duplicate.is_duplicate());
        assert!(ring_duplicate.is_expanded());
        assert_eq!(ring_duplicate.get_distance(), 1);
        assert_eq!(ring_duplicate.get_atomic_num_fraction().tuple(), (6, 1));

        let implicit_h = root.new_implicit_hydrogen_child(&mut cipmol).unwrap();
        assert_eq!(implicit_h.atom_idx(), None);
        assert_eq!(implicit_h.get_atom_idx().unwrap(), CipNode::NO_ATOM_INDEX);
        assert_eq!(implicit_h.get_atomic_num(&cipmol).unwrap(), 1);
        assert_eq!(implicit_h.get_mass_num(&cipmol).unwrap(), 0);
        assert!(implicit_h.is_duplicate_or_h());
        assert!(implicit_h.is_expanded());

        let mut aux_node = root.clone();
        aux_node.set_aux(Descriptor::R);
        assert_eq!(aux_node.get_aux(), Descriptor::R);
    }

    #[test]
    fn ciplabeler_node_bond_duplicate_uses_parent_mancude_fraction() {
        let molecule = Molecule::from_smiles("c1ccncc1").unwrap();
        let mut cipmol = CipMol::new(&molecule);
        let parent = CipNode::new(
            5,
            vec![0, 0, 1, 0, 0, 0],
            Some(2),
            RationalI32::new(13, 2),
            1,
            0,
            &cipmol,
        )
        .unwrap();

        let duplicate = parent
            .new_bond_duplicate_child(1, Some(1), &mut cipmol)
            .unwrap();
        assert_eq!(duplicate.atom_idx(), Some(1));
        assert_eq!(duplicate.get_atomic_num_fraction().tuple(), (13, 2));
        assert_eq!(duplicate.get_distance(), 0);
        assert_eq!(duplicate.get_atomic_mass(), 0.0);
        assert!(duplicate.is_set(CipNode::BOND_DUPLICATE));
        assert!(duplicate.is_duplicate());
        assert!(duplicate.is_duplicate_or_h());
        assert!(duplicate.is_expanded());
        assert!(duplicate.is_terminal());
        assert!(duplicate.edges.capacity() >= 4);
    }

    #[test]
    fn ciplabeler_node_terminal_state_tracks_visit_expansion_and_edge_count() {
        let molecule = Molecule::from_smiles("C").unwrap();
        let cipmol = CipMol::new(&molecule);

        let empty_visit = CipNode::new(
            0,
            Vec::new(),
            Some(0),
            RationalI32::new(6, 1),
            1,
            0,
            &cipmol,
        )
        .unwrap();
        assert!(empty_visit.is_expanded());
        assert!(empty_visit.is_terminal());

        let mut expanded = CipNode::new(
            0,
            vec![1],
            Some(0),
            RationalI32::new(6, 1),
            1,
            CipNode::EXPANDED,
            &cipmol,
        )
        .unwrap();
        assert!(!expanded.is_terminal());
        expanded.add(CipEdgeId::new(0));
        assert!(expanded.is_terminal());
        expanded.add(CipEdgeId::new(1));
        assert!(!expanded.is_terminal());
    }

    #[test]
    fn ciplabeler_node_atom_edge_filter_skips_duplicate_end_nodes() {
        let molecule = Molecule::from_smiles("C=C").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 0, true).unwrap();
        let root = digraph.get_current_root();

        let matching = digraph.node_edges_for_atom(root, Some(1)).unwrap();
        assert_eq!(matching.len(), 1);
        let edge = digraph.edge(matching[0]);
        assert_eq!(digraph.node(edge.get_end()).atom_idx(), Some(1));
        assert!(!digraph.node(edge.get_end()).is_duplicate());
    }

    #[test]
    fn ciplabeler_node_non_terminal_out_edges_follow_direction_and_terminal_state() {
        let molecule = Molecule::from_smiles("CCC").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 1, false).unwrap();
        let root = digraph.get_current_root();

        let outgoing = digraph.non_terminal_out_edges(root).unwrap();
        let atom_indices = outgoing
            .iter()
            .map(|edge_id| {
                let edge = digraph.edge(*edge_id);
                assert!(edge.is_beg(root));
                assert!(!digraph.node(edge.get_end()).is_terminal());
                digraph.node(edge.get_end()).atom_idx().unwrap()
            })
            .collect::<Vec<_>>();
        assert_eq!(atom_indices, vec![0, 2]);
    }

    #[test]
    fn ciplabeler_sort_priority_edge_endpoints_aux_and_flip_match_rdkit() {
        let mut edge = CipEdge::new(CipNodeId::new(0), CipNodeId::new(1), Some(7));
        assert_eq!(edge.get_beg(), CipNodeId::new(0));
        assert_eq!(edge.get_end(), CipNodeId::new(1));
        assert_eq!(edge.get_bond_idx(), Some(7));
        assert!(edge.is_beg(CipNodeId::new(0)));
        assert!(edge.is_end(CipNodeId::new(1)));
        assert_eq!(
            edge.get_other(CipEdgeId::new(4), CipNodeId::new(0))
                .unwrap(),
            CipNodeId::new(1)
        );
        assert_eq!(
            edge.get_other(CipEdgeId::new(4), CipNodeId::new(1))
                .unwrap(),
            CipNodeId::new(0)
        );
        assert!(matches!(
            edge.get_other(CipEdgeId::new(4), CipNodeId::new(2)),
            Err(CipLabelerError::EdgeEndpointMismatch { edge: 4, node: 2 })
        ));

        edge.set_aux(Descriptor::seqTrans);
        assert_eq!(edge.get_aux(), Descriptor::seqTrans);
        edge.flip();
        assert_eq!(edge.get_beg(), CipNodeId::new(1));
        assert_eq!(edge.get_end(), CipNodeId::new(0));
    }

    #[test]
    fn ciplabeler_node_raw_edges_fail_closed_until_digraph_lazy_expansion() {
        let molecule = Molecule::from_smiles("CC").unwrap();
        let cipmol = CipMol::new(&molecule);
        let node = CipNode::new(
            3,
            vec![1, 0],
            Some(0),
            RationalI32::new(6, 1),
            1,
            0,
            &cipmol,
        )
        .unwrap();
        assert!(matches!(
            node.get_edges(),
            Err(CipLabelerError::InvalidInternalState {
                detail: "use CipDigraph::node_edges for Node::getEdges lazy expansion"
            })
        ));
    }

    #[test]
    fn ciplabeler_digraph_root_construction_and_rule6_ref_match_rdkit() {
        let molecule = Molecule::from_smiles("CCO").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 1, false).unwrap();

        assert_eq!(digraph.get_original_root(), CipNodeId::new(0));
        assert_eq!(digraph.get_current_root(), CipNodeId::new(0));
        assert_eq!(digraph.get_num_nodes(), 1);
        assert_eq!(digraph.node(digraph.get_current_root()).atom_idx(), Some(1));
        assert_eq!(digraph.node(digraph.get_current_root()).get_distance(), 1);
        assert_eq!(digraph.mol().get_num_atoms(), 3);
        assert_eq!(digraph.get_rule6_ref(), None);

        digraph.set_rule6_ref(Some(2)).unwrap();
        assert_eq!(digraph.get_rule6_ref(), Some(2));
        digraph.set_rule6_ref(None).unwrap();
        assert_eq!(digraph.get_rule6_ref(), None);
    }

    #[test]
    fn ciplabeler_digraph_lazy_expansion_creates_explicit_and_implicit_nodes() {
        let molecule = Molecule::from_smiles("CC").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 0, false).unwrap();
        let root = digraph.get_current_root();

        let root_edges = digraph.node_edges(root).unwrap();
        assert_eq!(root_edges.len(), 4);
        assert_eq!(digraph.get_num_nodes(), 5);
        assert!(digraph.node(root).is_expanded());

        let explicit_edges = root_edges
            .iter()
            .filter(|edge_id| digraph.edge(**edge_id).get_bond_idx().is_some())
            .copied()
            .collect::<Vec<_>>();
        assert_eq!(explicit_edges.len(), 1);
        let explicit_end = digraph.edge(explicit_edges[0]).get_end();
        assert_eq!(digraph.node(explicit_end).atom_idx(), Some(1));
        assert_eq!(digraph.node(explicit_end).get_distance(), 2);

        let implicit_count = root_edges
            .iter()
            .filter(|edge_id| digraph.edge(**edge_id).get_bond_idx().is_none())
            .count();
        assert_eq!(implicit_count, 3);

        let node_count = digraph.get_num_nodes();
        let repeated_edges = digraph.node_edges(root).unwrap();
        assert_eq!(repeated_edges, root_edges);
        assert_eq!(digraph.get_num_nodes(), node_count);
    }

    #[test]
    fn ciplabeler_digraph_expansion_creates_bond_and_ring_duplicates() {
        let normal_ethene = Molecule::from_smiles("C=C").unwrap();
        let mut normal_graph = CipDigraph::new(&normal_ethene, 0, false).unwrap();
        let normal_root = normal_graph.get_current_root();
        let normal_root_edges = normal_graph.node_edges(normal_root).unwrap();
        assert!(!normal_root_edges.iter().any(|edge_id| {
            let end = normal_graph.edge(*edge_id).get_end();
            normal_graph.node(end).is_duplicate()
        }));

        let ethene = Molecule::from_smiles("C=C").unwrap();
        let mut ethene_graph = CipDigraph::new(&ethene, 0, true).unwrap();
        let root = ethene_graph.get_current_root();
        let root_edges = ethene_graph.node_edges(root).unwrap();
        let duplicate_edges = root_edges
            .iter()
            .filter(|edge_id| {
                let end = ethene_graph.edge(**edge_id).get_end();
                ethene_graph.node(end).is_duplicate()
            })
            .count();
        assert_eq!(duplicate_edges, 1);

        let cyclopropane = Molecule::from_smiles("C1CC1").unwrap();
        let mut ring_graph = CipDigraph::new(&cyclopropane, 0, false).unwrap();
        let root = ring_graph.get_current_root();
        let first_layer = ring_graph.node_edges(root).unwrap();
        let non_root_neighbor = first_layer
            .iter()
            .map(|edge_id| ring_graph.edge(*edge_id).get_end())
            .find(|node_id| ring_graph.node(*node_id).atom_idx() == Some(1))
            .unwrap();
        let second_layer = ring_graph.node_edges(non_root_neighbor).unwrap();
        let closing_node = second_layer
            .iter()
            .map(|edge_id| ring_graph.edge(*edge_id).get_end())
            .find(|node_id| ring_graph.node(*node_id).atom_idx() == Some(2))
            .unwrap();
        let third_layer = ring_graph.node_edges(closing_node).unwrap();
        assert!(third_layer.iter().any(|edge_id| {
            let end = ring_graph.edge(*edge_id).get_end();
            ring_graph.node(end).is_set(CipNode::RING_DUPLICATE)
        }));
    }

    #[test]
    fn ciplabeler_digraph_get_nodes_and_change_root_match_rdkit_direction_rules() {
        let molecule = Molecule::from_smiles("CCC").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 1, false).unwrap();
        let root = digraph.get_current_root();
        let root_edges = digraph.node_edges(root).unwrap();
        let left = root_edges
            .iter()
            .map(|edge_id| digraph.edge(*edge_id).get_end())
            .find(|node_id| digraph.node(*node_id).atom_idx() == Some(0))
            .unwrap();

        digraph.node_edges(left).unwrap();
        let center_nodes = digraph.get_nodes(1).unwrap();
        assert!(center_nodes.iter().any(|node| *node == root));

        digraph.change_root(left).unwrap();
        assert_eq!(digraph.get_current_root(), left);
        assert_eq!(digraph.get_original_root(), root);
        let flipped_back_edge = digraph
            .node_edges(left)
            .unwrap()
            .into_iter()
            .find(|edge_id| digraph.edge(*edge_id).get_end() == root)
            .unwrap();
        assert!(digraph.edge(flipped_back_edge).is_beg(left));
    }

    #[test]
    fn ciplabeler_digraph_get_nodes_returns_repeated_atoms_in_distance_order() {
        let molecule = Molecule::from_smiles("C1CC1").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 0, false).unwrap();

        let nodes = digraph.get_nodes(0).unwrap();
        assert_eq!(nodes[0], digraph.get_current_root());
        assert!(nodes.len() > 1);
        let distances = nodes
            .iter()
            .map(|node| digraph.node(*node).get_distance())
            .collect::<Vec<_>>();
        assert!(distances.windows(2).all(|pair| pair[0] <= pair[1]));
        assert!(
            nodes
                .iter()
                .all(|node| digraph.node(*node).atom_idx() == Some(0))
        );
    }

    struct AtomicNumberSortRule;

    impl CipSequenceRule for AtomicNumberSortRule {
        fn compare(
            &self,
            digraph: &mut CipDigraph<'_>,
            _context: &mut CipLabelerContext,
            a: CipEdgeId,
            b: CipEdgeId,
        ) -> Result<i32, CipLabelerError> {
            let a_end = digraph.edge(a).get_end();
            let b_end = digraph.edge(b).get_end();
            let a_atomic_num = i32::from(digraph.node(a_end).get_atomic_num(digraph.mol())?);
            let b_atomic_num = i32::from(digraph.node(b_end).get_atomic_num(digraph.mol())?);
            Ok(match a_atomic_num.cmp(&b_atomic_num) {
                std::cmp::Ordering::Less => -1,
                std::cmp::Ordering::Equal => 0,
                std::cmp::Ordering::Greater => 1,
            })
        }
    }

    struct PseudoAsymSortRule;

    impl CipSequenceRule for PseudoAsymSortRule {
        fn compare(
            &self,
            digraph: &mut CipDigraph<'_>,
            _context: &mut CipLabelerContext,
            a: CipEdgeId,
            b: CipEdgeId,
        ) -> Result<i32, CipLabelerError> {
            let a_end = digraph.edge(a).get_end();
            let b_end = digraph.edge(b).get_end();
            let a_atomic_num = i32::from(digraph.node(a_end).get_atomic_num(digraph.mol())?);
            let b_atomic_num = i32::from(digraph.node(b_end).get_atomic_num(digraph.mol())?);
            Ok(match a_atomic_num.cmp(&b_atomic_num) {
                std::cmp::Ordering::Less => -2,
                std::cmp::Ordering::Equal => 0,
                std::cmp::Ordering::Greater => 2,
            })
        }
    }

    struct SinglePseudoAsymSortRule;

    impl CipSequenceRule for SinglePseudoAsymSortRule {
        fn compare(
            &self,
            digraph: &mut CipDigraph<'_>,
            _context: &mut CipLabelerContext,
            a: CipEdgeId,
            b: CipEdgeId,
        ) -> Result<i32, CipLabelerError> {
            let a_end = digraph.edge(a).get_end();
            let b_end = digraph.edge(b).get_end();
            let a_atomic_num = i32::from(digraph.node(a_end).get_atomic_num(digraph.mol())?);
            let b_atomic_num = i32::from(digraph.node(b_end).get_atomic_num(digraph.mol())?);
            let magnitude = if (a_atomic_num == 35 && b_atomic_num == 53)
                || (a_atomic_num == 53 && b_atomic_num == 35)
            {
                2
            } else {
                1
            };
            Ok(match a_atomic_num.cmp(&b_atomic_num) {
                std::cmp::Ordering::Less => -magnitude,
                std::cmp::Ordering::Equal => 0,
                std::cmp::Ordering::Greater => magnitude,
            })
        }
    }

    struct BondLabelSortRule;

    impl CipSequenceRule for BondLabelSortRule {
        fn compare(
            &self,
            digraph: &mut CipDigraph<'_>,
            _context: &mut CipLabelerContext,
            a: CipEdgeId,
            b: CipEdgeId,
        ) -> Result<i32, CipLabelerError> {
            let a_label = self.get_bond_label(digraph.edge(a));
            let b_label = self.get_bond_label(digraph.edge(b));
            Ok(match a_label.cmp(&b_label) {
                std::cmp::Ordering::Less => -1,
                std::cmp::Ordering::Equal => 0,
                std::cmp::Ordering::Greater => 1,
            })
        }
    }

    struct AlwaysEqualSortRule;

    impl CipSequenceRule for AlwaysEqualSortRule {
        fn compare(
            &self,
            _digraph: &mut CipDigraph<'_>,
            _context: &mut CipLabelerContext,
            _a: CipEdgeId,
            _b: CipEdgeId,
        ) -> Result<i32, CipLabelerError> {
            Ok(0)
        }
    }

    struct CountingAtomicNumberSortRule {
        calls: Rc<Cell<usize>>,
    }

    impl CipSequenceRule for CountingAtomicNumberSortRule {
        fn compare(
            &self,
            digraph: &mut CipDigraph<'_>,
            _context: &mut CipLabelerContext,
            a: CipEdgeId,
            b: CipEdgeId,
        ) -> Result<i32, CipLabelerError> {
            self.calls.set(self.calls.get() + 1);
            let a_end = digraph.edge(a).get_end();
            let b_end = digraph.edge(b).get_end();
            let a_atomic_num = digraph.node(a_end).get_atomic_num(digraph.mol())?;
            let b_atomic_num = digraph.node(b_end).get_atomic_num(digraph.mol())?;
            Ok(three_way_comparison_i32(
                i32::from(a_atomic_num),
                i32::from(b_atomic_num),
            ))
        }
    }

    fn edge_end_atomic_nums(
        digraph: &CipDigraph<'_>,
        edges: &[CipEdgeId],
    ) -> Result<Vec<u8>, CipLabelerError> {
        edges
            .iter()
            .map(|edge_id| {
                let end = digraph.edge(*edge_id).get_end();
                digraph.node(end).get_atomic_num(digraph.mol())
            })
            .collect()
    }

    #[test]
    fn ciplabeler_sort_priority_orders_edges_and_reports_unique_like_rdkit() {
        let priority = CipPriority::new(true, false);
        assert!(priority.is_unique());
        assert!(!priority.is_pseudo_asymetric());

        let molecule = Molecule::from_smiles("C(F)(Cl)Br").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 0, false).unwrap();
        let root = digraph.get_current_root();
        let mut edges = digraph.node_edges(root).unwrap();

        let rule = AtomicNumberSortRule;
        let sorter = CipSort::new(&rule);
        let mut context = CipLabelerContext::new(0);
        assert_eq!(sorter.get_rules().len(), 1);
        let priority = sorter
            .prioritize(&mut digraph, &mut context, root, &mut edges, true)
            .unwrap();

        assert!(priority.is_unique());
        assert!(!priority.is_pseudo_asymetric());
        assert_eq!(
            edge_end_atomic_nums(&digraph, &edges).unwrap(),
            vec![35, 17, 9, 1]
        );
    }

    #[test]
    fn ciplabeler_sort_priority_groups_equal_edges_like_rdkit() {
        let molecule = Molecule::from_smiles("CC(C)N").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 1, false).unwrap();
        let root = digraph.get_current_root();
        let mut edges = digraph.node_edges(root).unwrap();

        let rule = AtomicNumberSortRule;
        let sorter = CipSort::from_rules(vec![&rule]);
        let mut context = CipLabelerContext::new(0);
        let priority = sorter
            .prioritize(&mut digraph, &mut context, root, &mut edges, true)
            .unwrap();
        let groups = sorter
            .get_groups(&mut digraph, &mut context, &edges)
            .unwrap();

        assert!(!priority.is_unique());
        assert!(!priority.is_pseudo_asymetric());
        assert_eq!(
            edge_end_atomic_nums(&digraph, &edges).unwrap(),
            vec![7, 6, 6, 1]
        );
        assert_eq!(
            groups.iter().map(Vec::len).collect::<Vec<_>>(),
            vec![1, 2, 1]
        );
    }

    #[test]
    fn ciplabeler_sort_priority_counts_single_pseudoasym_comparison_like_rdkit() {
        let molecule = Molecule::from_smiles("FON").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 1, false).unwrap();
        let root = digraph.get_current_root();
        let mut edges = digraph.node_edges(root).unwrap();

        let rule = PseudoAsymSortRule;
        let sorter = CipSort::new(&rule);
        let mut context = CipLabelerContext::new(0);
        let priority = sorter
            .prioritize(&mut digraph, &mut context, root, &mut edges, true)
            .unwrap();

        assert!(priority.is_unique());
        assert!(priority.is_pseudo_asymetric());
        assert_eq!(edge_end_atomic_nums(&digraph, &edges).unwrap(), vec![9, 7]);
    }

    #[test]
    fn ciplabeler_sort_priority_forwards_deep_and_uses_first_nonzero_rule() {
        let molecule = Molecule::from_smiles("FON").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 1, false).unwrap();
        let root = digraph.get_current_root();
        let original_edges = digraph.node_edges(root).unwrap();

        let equal = AlwaysEqualSortRule;
        let atomic_number = AtomicNumberSortRule;
        let sorter = CipSort::from_rules(vec![&equal, &atomic_number]);
        assert_eq!(sorter.get_rules().len(), 2);

        let mut shallow_edges = original_edges.clone();
        let mut shallow_context = CipLabelerContext::with_remaining_call_count(1);
        let shallow_priority = sorter
            .prioritize(
                &mut digraph,
                &mut shallow_context,
                root,
                &mut shallow_edges,
                false,
            )
            .unwrap();
        assert!(shallow_priority.is_unique());
        assert_eq!(
            edge_end_atomic_nums(&digraph, &shallow_edges).unwrap(),
            vec![9, 7]
        );

        let mut deep_edges = original_edges;
        let mut deep_context = CipLabelerContext::with_remaining_call_count(1);
        assert_eq!(
            sorter.prioritize(&mut digraph, &mut deep_context, root, &mut deep_edges, true,),
            Err(CipLabelerError::MaxIterationsExceeded)
        );
    }

    #[test]
    fn ciplabeler_sequence_rule_get_bond_label_matches_rdkit() {
        let molecule = Molecule::from_smiles("CC").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 0, false).unwrap();
        let root = digraph.get_current_root();
        let edges = digraph.node_edges(root).unwrap();
        let rule = BondLabelSortRule;

        let real_bond = edges
            .iter()
            .copied()
            .find(|edge_id| digraph.edge(*edge_id).get_bond_idx().is_some())
            .unwrap();
        let implicit_h = edges
            .iter()
            .copied()
            .find(|edge_id| digraph.edge(*edge_id).get_bond_idx().is_none())
            .unwrap();

        assert_eq!(
            rule.get_bond_label(digraph.edge(real_bond)),
            Descriptor::None
        );
        assert_eq!(
            rule.get_bond_label(digraph.edge(implicit_h)),
            Descriptor::None
        );

        digraph.edges[real_bond.index()].set_aux(Descriptor::seqTrans);
        digraph.edges[implicit_h.index()].set_aux(Descriptor::seqCis);

        assert_eq!(
            rule.get_bond_label(digraph.edge(real_bond)),
            Descriptor::seqTrans
        );
        assert_eq!(
            rule.get_bond_label(digraph.edge(implicit_h)),
            Descriptor::None
        );
    }

    #[test]
    fn ciplabeler_sequence_rule_get_comparison_deep_and_shallow_match_rdkit() {
        let molecule = Molecule::from_smiles("C(CF)CCl").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 0, false).unwrap();
        let root = digraph.get_current_root();
        let edges = digraph.node_edges(root).unwrap();
        let carbon_edges = edges
            .iter()
            .copied()
            .filter(|edge_id| {
                let end = digraph.edge(*edge_id).get_end();
                digraph.node(end).get_atomic_num(digraph.mol()).unwrap() == 6
            })
            .collect::<Vec<_>>();
        assert_eq!(carbon_edges.len(), 2);

        let rule = AtomicNumberSortRule;
        let mut context = CipLabelerContext::new(0);
        assert_eq!(
            rule.get_comparison(
                &mut digraph,
                &mut context,
                carbon_edges[0],
                carbon_edges[1],
                false,
            )
            .unwrap(),
            0
        );

        let mut context = CipLabelerContext::new(0);
        assert_eq!(
            rule.get_comparison(
                &mut digraph,
                &mut context,
                carbon_edges[0],
                carbon_edges[1],
                true,
            )
            .unwrap(),
            -1
        );
    }

    #[test]
    fn ciplabeler_sequence_rule_recursive_compare_returns_early_nonzero_like_rdkit() {
        let molecule = Molecule::from_smiles("C(F)Cl").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 0, false).unwrap();
        let root = digraph.get_current_root();
        let edges = digraph.node_edges(root).unwrap();
        let fluorine = edges
            .iter()
            .copied()
            .find(|edge_id| {
                let end = digraph.edge(*edge_id).get_end();
                digraph.node(end).get_atomic_num(digraph.mol()).unwrap() == 9
            })
            .unwrap();
        let chlorine = edges
            .iter()
            .copied()
            .find(|edge_id| {
                let end = digraph.edge(*edge_id).get_end();
                digraph.node(end).get_atomic_num(digraph.mol()).unwrap() == 17
            })
            .unwrap();

        let rule = AtomicNumberSortRule;
        let mut context = CipLabelerContext::new(0);
        assert_eq!(
            rule.recursive_compare(&mut digraph, &mut context, fluorine, chlorine)
                .unwrap(),
            -1
        );
    }

    #[test]
    fn ciplabeler_sequence_rule_recursive_compare_enforces_remaining_call_count_like_rdkit() {
        let molecule = Molecule::from_smiles("C(F)Cl").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 0, false).unwrap();
        let root = digraph.get_current_root();
        let edges = digraph.node_edges(root).unwrap();
        let fluorine = edges
            .iter()
            .copied()
            .find(|edge_id| {
                let end = digraph.edge(*edge_id).get_end();
                digraph.node(end).get_atomic_num(digraph.mol()).unwrap() == 9
            })
            .unwrap();
        let chlorine = edges
            .iter()
            .copied()
            .find(|edge_id| {
                let end = digraph.edge(*edge_id).get_end();
                digraph.node(end).get_atomic_num(digraph.mol()).unwrap() == 17
            })
            .unwrap();

        let rule = AtomicNumberSortRule;
        let mut context = CipLabelerContext::with_remaining_call_count(1);
        assert!(matches!(
            rule.recursive_compare(&mut digraph, &mut context, fluorine, chlorine),
            Err(CipLabelerError::MaxIterationsExceeded)
        ));
    }

    #[test]
    fn ciplabeler_sequence_rule_sort_delegates_to_single_rule_sorter_like_rdkit() {
        let molecule = Molecule::from_smiles("C(F)(Cl)Br").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 0, false).unwrap();
        let root = digraph.get_current_root();
        let mut edges = digraph.node_edges(root).unwrap();

        let rule = AtomicNumberSortRule;
        let mut context = CipLabelerContext::new(0);
        let priority = rule
            .sort(&mut digraph, &mut context, root, &mut edges, true)
            .unwrap();

        assert!(priority.is_unique());
        assert!(!priority.is_pseudo_asymetric());
        assert_eq!(
            edge_end_atomic_nums(&digraph, &edges).unwrap(),
            vec![35, 17, 9, 1]
        );
    }

    #[test]
    fn ciplabeler_sequence_rule_are_up_edges_skip_and_error_match_rdkit() {
        let molecule = Molecule::from_smiles("CCC").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 1, false).unwrap();
        let root = digraph.get_current_root();
        let root_edges = digraph.node_edges(root).unwrap();
        let left = root_edges
            .iter()
            .copied()
            .find(|edge_id| {
                let end = digraph.edge(*edge_id).get_end();
                digraph.node(end).atom_idx() == Some(0)
            })
            .unwrap();
        let right = root_edges
            .iter()
            .copied()
            .find(|edge_id| {
                let end = digraph.edge(*edge_id).get_end();
                digraph.node(end).atom_idx() == Some(2)
            })
            .unwrap();
        let left_node = digraph.edge(left).get_end();
        let right_node = digraph.edge(right).get_end();
        digraph.node_edges(left_node).unwrap();

        let rule = AtomicNumberSortRule;
        assert!(
            rule.are_up_edges(&digraph, left_node, left_node, left, left)
                .unwrap()
        );
        assert!(matches!(
            rule.are_up_edges(&digraph, left_node, root, left, right),
            Err(CipLabelerError::UnexpectedUpEdgeOrdering)
        ));
    }

    #[test]
    fn ciplabeler_rules_constructor_count_and_null_error_match_rdkit() {
        let rules = CipRules::new(vec![
            Box::new(AlwaysEqualSortRule),
            Box::new(AtomicNumberSortRule),
        ])
        .unwrap();
        assert_eq!(rules.get_num_sub_rules(), 2);
        assert_eq!(rules.get_sorter().get_rules().len(), 1);

        let mut empty = CipRules::new(Vec::new()).unwrap();
        assert!(matches!(
            empty.add(None),
            Err(CipLabelerError::NoSequenceRuleProvided)
        ));
    }

    #[test]
    fn ciplabeler_sequence_rule_subrule_sorters_use_addition_prefixes() {
        let molecule = Molecule::from_smiles("C(CF)CCl").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 0, false).unwrap();
        let root = digraph.get_current_root();
        let carbon_edges = digraph
            .node_edges(root)
            .unwrap()
            .into_iter()
            .filter(|edge_id| {
                let end = digraph.edge(*edge_id).get_end();
                digraph.node(end).get_atomic_num(digraph.mol()).unwrap() == 6
            })
            .collect::<Vec<_>>();
        let calls = Rc::new(Cell::new(0));
        let rules = CipRules::new(vec![
            Box::new(AlwaysEqualSortRule),
            Box::new(CountingAtomicNumberSortRule {
                calls: Rc::clone(&calls),
            }),
        ])
        .unwrap();
        let rule_refs = rules.rule_refs();

        let mut context = CipLabelerContext::new(0);
        assert_eq!(
            rule_refs[0]
                .recursive_compare_with_sort_rules(
                    &rule_refs[..=0],
                    &mut digraph,
                    &mut context,
                    carbon_edges[0],
                    carbon_edges[1],
                )
                .unwrap(),
            0
        );
        assert_eq!(calls.get(), 0);

        let mut context = CipLabelerContext::new(0);
        rule_refs[0]
            .recursive_compare_with_sort_rules(
                &rule_refs,
                &mut digraph,
                &mut context,
                carbon_edges[0],
                carbon_edges[1],
            )
            .unwrap();
        assert!(calls.get() > 0);
    }

    #[test]
    fn ciplabeler_rules_compare_tries_subrules_in_order_like_rdkit() {
        let molecule = Molecule::from_smiles("C(F)Cl").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 0, false).unwrap();
        let root = digraph.get_current_root();
        let edges = digraph.node_edges(root).unwrap();
        let fluorine = edges
            .iter()
            .copied()
            .find(|edge_id| {
                let end = digraph.edge(*edge_id).get_end();
                digraph.node(end).get_atomic_num(digraph.mol()).unwrap() == 9
            })
            .unwrap();
        let chlorine = edges
            .iter()
            .copied()
            .find(|edge_id| {
                let end = digraph.edge(*edge_id).get_end();
                digraph.node(end).get_atomic_num(digraph.mol()).unwrap() == 17
            })
            .unwrap();

        let rules = CipRules::new(vec![
            Box::new(AlwaysEqualSortRule),
            Box::new(AtomicNumberSortRule),
        ])
        .unwrap();
        let mut context = CipLabelerContext::new(0);
        assert_eq!(
            rules
                .get_comparison(&mut digraph, &mut context, fluorine, chlorine, false)
                .unwrap(),
            -1
        );
    }

    #[test]
    fn ciplabeler_rules_get_comparison_ignores_deep_flag_like_rdkit() {
        let molecule = Molecule::from_smiles("C(CF)CCl").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 0, false).unwrap();
        let root = digraph.get_current_root();
        let edges = digraph.node_edges(root).unwrap();
        let carbon_edges = edges
            .iter()
            .copied()
            .filter(|edge_id| {
                let end = digraph.edge(*edge_id).get_end();
                digraph.node(end).get_atomic_num(digraph.mol()).unwrap() == 6
            })
            .collect::<Vec<_>>();
        assert_eq!(carbon_edges.len(), 2);

        let rules = CipRules::new(vec![Box::new(AtomicNumberSortRule)]).unwrap();
        let mut shallow_context = CipLabelerContext::new(0);
        let mut deep_context = CipLabelerContext::new(0);
        let shallow = rules
            .get_comparison(
                &mut digraph,
                &mut shallow_context,
                carbon_edges[0],
                carbon_edges[1],
                false,
            )
            .unwrap();
        let deep = rules
            .get_comparison(
                &mut digraph,
                &mut deep_context,
                carbon_edges[0],
                carbon_edges[1],
                true,
            )
            .unwrap();
        assert_eq!(shallow, -1);
        assert_eq!(deep, shallow);
    }

    #[test]
    fn ciplabeler_rules_sort_uses_rules_own_sorter_like_rdkit() {
        let molecule = Molecule::from_smiles("C(F)(Cl)Br").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 0, false).unwrap();
        let root = digraph.get_current_root();
        let mut edges = digraph.node_edges(root).unwrap();

        let rules = CipRules::new(vec![
            Box::new(AlwaysEqualSortRule),
            Box::new(AtomicNumberSortRule),
        ])
        .unwrap();
        let mut context = CipLabelerContext::new(0);
        let priority = rules
            .sort(&mut digraph, &mut context, root, &mut edges, true)
            .unwrap();

        assert!(priority.is_unique());
        assert!(!priority.is_pseudo_asymetric());
        assert_eq!(
            edge_end_atomic_nums(&digraph, &edges).unwrap(),
            vec![35, 17, 9, 1]
        );
    }

    #[test]
    fn ciplabeler_rules_rule1a_compares_fractional_atomic_number_like_rdkit() {
        let molecule = Molecule::from_smiles("CC").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 0, false).unwrap();
        let root = digraph.get_current_root();
        let carbon_node = digraph
            .add_node(vec![], Some(1), RationalI32::new(6, 1), 2, 0)
            .unwrap();
        let fractional_node = digraph
            .add_node(vec![], Some(1), RationalI32::new(13, 2), 2, 0)
            .unwrap();
        digraph.add_edge(root, Some(0), carbon_node);
        let carbon = CipEdgeId::new(digraph.edges.len() - 1);
        digraph.add_edge(root, Some(0), fractional_node);
        let fractional_n = CipEdgeId::new(digraph.edges.len() - 1);

        let mut context = CipLabelerContext::new(0);
        assert_eq!(
            CipRule1a
                .compare(&mut digraph, &mut context, carbon, fractional_n)
                .unwrap(),
            -1
        );
        assert_eq!(
            CipRule1a
                .compare(&mut digraph, &mut context, fractional_n, carbon)
                .unwrap(),
            1
        );
        assert_eq!(
            CipRule1a
                .compare(&mut digraph, &mut context, carbon, carbon)
                .unwrap(),
            0
        );
    }

    #[test]
    fn ciplabeler_rules_rule1b_orders_ring_duplicates_like_rdkit_non_iupac_2013() {
        let molecule = Molecule::from_smiles("CC").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 0, false).unwrap();
        let root = digraph.get_current_root();
        let non_duplicate_node = digraph
            .add_node(vec![], Some(1), RationalI32::new(6, 1), 2, 0)
            .unwrap();
        let ring_duplicate_node = digraph
            .add_node(
                vec![],
                Some(1),
                RationalI32::new(6, 1),
                3,
                CipNode::RING_DUPLICATE,
            )
            .unwrap();
        digraph.add_edge(root, Some(0), non_duplicate_node);
        let non_duplicate = CipEdgeId::new(digraph.edges.len() - 1);
        digraph.add_edge(root, Some(0), ring_duplicate_node);
        let ring_duplicate = CipEdgeId::new(digraph.edges.len() - 1);

        let mut context = CipLabelerContext::new(0);
        assert_eq!(
            CipRule1b
                .compare(&mut digraph, &mut context, ring_duplicate, non_duplicate)
                .unwrap(),
            1
        );
        assert_eq!(
            CipRule1b
                .compare(&mut digraph, &mut context, non_duplicate, ring_duplicate)
                .unwrap(),
            -1
        );

        let nearer_ring_duplicate_node = digraph
            .add_node(
                vec![],
                Some(1),
                RationalI32::new(6, 1),
                2,
                CipNode::RING_DUPLICATE,
            )
            .unwrap();
        digraph.add_edge(root, Some(0), nearer_ring_duplicate_node);
        let nearer_ring_duplicate = CipEdgeId::new(digraph.edges.len() - 1);
        assert_eq!(
            CipRule1b
                .compare(
                    &mut digraph,
                    &mut context,
                    nearer_ring_duplicate,
                    ring_duplicate,
                )
                .unwrap(),
            1
        );
        assert_eq!(
            CipRule1b
                .compare(&mut digraph, &mut context, non_duplicate, non_duplicate)
                .unwrap(),
            0
        );
    }

    #[test]
    fn ciplabeler_rules_rule2_compares_isotopic_mass_like_rdkit() {
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(AtomSpec::new(Element::C));
        let c12 = builder.add_atom(AtomSpec::new(Element::C).with_isotope(12));
        let c13 = builder.add_atom(AtomSpec::new(Element::C).with_isotope(13));
        builder
            .add_bond(BondSpec::new(center, c12, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, c13, BondOrder::Single))
            .unwrap();
        let molecule = builder.build().unwrap();
        let mut digraph = CipDigraph::new(&molecule, center.index(), false).unwrap();
        let root = digraph.get_current_root();
        let root_edges = digraph.node_edges(root).unwrap();
        let c13 = root_edges
            .iter()
            .copied()
            .find(|edge_id| {
                let end = digraph.edge(*edge_id).get_end();
                digraph.node(end).atom_idx() == Some(c13.index())
            })
            .unwrap();
        let c12 = root_edges
            .iter()
            .copied()
            .find(|edge_id| {
                let end = digraph.edge(*edge_id).get_end();
                digraph.node(end).atom_idx() == Some(c12.index())
            })
            .unwrap();

        let mut context = CipLabelerContext::new(0);
        assert_eq!(
            CipRule2
                .compare(&mut digraph, &mut context, c12, c13)
                .unwrap(),
            -1
        );
        let mut context = CipLabelerContext::new(0);
        let c13_mass = digraph.node(digraph.edge(c13).get_end()).get_atomic_mass();
        let c12_mass = digraph.node(digraph.edge(c12).get_end()).get_atomic_mass();
        assert!(c13_mass > c12_mass);
        assert_eq!(
            CipRule2
                .compare(&mut digraph, &mut context, c13, c12)
                .unwrap(),
            1
        );
    }

    #[test]
    fn ciplabeler_rules_rule2_handles_zero_atomic_and_mass_numbers_like_rdkit() {
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(AtomSpec::new(Element::C));
        let dummy = builder.add_atom(AtomSpec::new(Element::DUMMY));
        let carbon = builder.add_atom(AtomSpec::new(Element::C));
        let nitrogen = builder.add_atom(AtomSpec::new(Element::N));
        builder
            .add_bond(BondSpec::new(center, dummy, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, carbon, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, nitrogen, BondOrder::Single))
            .unwrap();
        let molecule = builder.build().unwrap();
        let mut digraph = CipDigraph::new(&molecule, center.index(), false).unwrap();
        let root = digraph.get_current_root();
        let edges = digraph.node_edges(root).unwrap();
        let edge_for = |digraph: &CipDigraph<'_>, atom_idx| {
            edges
                .iter()
                .copied()
                .find(|edge_id| {
                    digraph.node(digraph.edge(*edge_id).get_end()).atom_idx() == Some(atom_idx)
                })
                .unwrap()
        };
        let dummy_edge = edge_for(&digraph, dummy.index());
        let carbon_edge = edge_for(&digraph, carbon.index());
        let nitrogen_edge = edge_for(&digraph, nitrogen.index());

        let mut context = CipLabelerContext::new(0);
        assert_eq!(
            CipRule2
                .compare(&mut digraph, &mut context, dummy_edge, dummy_edge)
                .unwrap(),
            0
        );
        assert_eq!(
            CipRule2
                .compare(&mut digraph, &mut context, dummy_edge, carbon_edge)
                .unwrap(),
            -1
        );
        assert_eq!(
            CipRule2
                .compare(&mut digraph, &mut context, carbon_edge, nitrogen_edge)
                .unwrap(),
            0
        );
    }

    #[test]
    fn ciplabeler_rules_rule3_orders_e_z_aux_labels_like_rdkit() {
        assert_eq!(CipRule3::ord(Descriptor::E), 1);
        assert_eq!(CipRule3::ord(Descriptor::Z), 2);
        assert_eq!(CipRule3::ord(Descriptor::R), 0);
        assert_eq!(CipRule3::ord(Descriptor::None), 0);

        let molecule = Molecule::from_smiles("CC").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 0, false).unwrap();
        let root = digraph.get_current_root();
        let mut edges = digraph.node_edges(root).unwrap();
        assert!(edges.len() >= 2);

        let e_edge = edges[0];
        let z_edge = edges[1];
        let e_node = digraph.edge(e_edge).get_end();
        let z_node = digraph.edge(z_edge).get_end();
        digraph.nodes[e_node.index()].set_aux(Descriptor::E);
        digraph.nodes[z_node.index()].set_aux(Descriptor::Z);

        let mut context = CipLabelerContext::new(0);
        assert_eq!(
            CipRule3
                .compare(&mut digraph, &mut context, e_edge, z_edge)
                .unwrap(),
            -1
        );
        edges.swap(0, 1);
        assert_eq!(
            CipRule3
                .compare(&mut digraph, &mut context, edges[0], edges[1])
                .unwrap(),
            1
        );
        digraph.nodes[e_node.index()].set_aux(Descriptor::R);
        digraph.nodes[z_node.index()].set_aux(Descriptor::S);
        assert_eq!(
            CipRule3
                .compare(&mut digraph, &mut context, e_edge, z_edge)
                .unwrap(),
            0
        );
    }

    #[test]
    fn ciplabeler_pairlist_ref_and_add_match_rdkit() {
        assert_eq!(CipPairList::ref_descriptor(Descriptor::R), Descriptor::R);
        assert_eq!(CipPairList::ref_descriptor(Descriptor::M), Descriptor::R);
        assert_eq!(
            CipPairList::ref_descriptor(Descriptor::seqCis),
            Descriptor::R
        );
        assert_eq!(CipPairList::ref_descriptor(Descriptor::S), Descriptor::S);
        assert_eq!(CipPairList::ref_descriptor(Descriptor::P), Descriptor::S);
        assert_eq!(
            CipPairList::ref_descriptor(Descriptor::seqTrans),
            Descriptor::S
        );
        assert_eq!(
            CipPairList::ref_descriptor(Descriptor::None),
            Descriptor::None
        );
        assert_eq!(CipPairList::ref_descriptor(Descriptor::r), Descriptor::None);

        let mut pairs = CipPairList::new();
        assert_eq!(pairs.to_rdkit_string(), "");
        assert!(!pairs.add(Descriptor::None));
        assert!(!pairs.add(Descriptor::r));
        assert!(pairs.add(Descriptor::M));
        assert_eq!(pairs.get_ref_descriptor(), Descriptor::R);
        assert_eq!(pairs.to_rdkit_string(), "R:");
    }

    #[test]
    fn ciplabeler_pairlist_pairing_bits_and_string_match_rdkit() {
        let mut pairs = CipPairList::with_ref(Descriptor::R);
        pairs.add(Descriptor::R);
        pairs.add(Descriptor::P);
        pairs.add(Descriptor::R);

        assert_eq!(pairs.get_pairing(), (1_u64 << 62) | (1_u64 << 60));
        assert_eq!(pairs.to_rdkit_string(), "R:lul");
    }

    #[test]
    fn ciplabeler_pairlist_head_tail_constructor_filters_descriptors_like_rdkit() {
        let mut head = CipPairList::with_ref(Descriptor::R);
        head.add(Descriptor::S);
        let mut tail = CipPairList::new();
        tail.add(Descriptor::None);
        tail.add(Descriptor::seqCis);
        tail.add(Descriptor::Unknown);

        let combined = CipPairList::from_head_tail(&head, &tail);
        assert_eq!(combined.to_rdkit_string(), "R:ul");
        assert_eq!(combined.get_pairing(), 1_u64 << 61);
    }

    #[test]
    fn ciplabeler_pairlist_compare_to_and_order_match_rdkit() {
        let mut like = CipPairList::with_ref(Descriptor::R);
        like.add(Descriptor::R);
        let mut unlike = CipPairList::with_ref(Descriptor::R);
        unlike.add(Descriptor::S);

        assert_eq!(like.compare_to(&unlike).unwrap(), 1);
        assert_eq!(unlike.compare_to(&like).unwrap(), -1);
        assert!(unlike.less_than(&like).unwrap());
        assert!(!like.less_than(&unlike).unwrap());

        let mut lists = vec![unlike.clone(), like.clone()];
        CipPairList::sort_descending(&mut lists).unwrap();
        assert_eq!(lists[0], like);
        assert_eq!(lists[1], unlike);

        let mut shorter = CipPairList::with_ref(Descriptor::R);
        shorter.add(Descriptor::R);
        let mut longer = shorter.clone();
        longer.add(Descriptor::R);
        assert!(matches!(
            shorter.compare_to(&longer),
            Err(CipLabelerError::DescriptorListLengthMismatch)
        ));
        assert!(matches!(
            CipPairList::sort_descending(&mut [shorter, longer]),
            Err(CipLabelerError::DescriptorListLengthMismatch)
        ));
    }

    #[test]
    fn ciplabeler_rules_rule4a_orders_descriptor_classes_like_rdkit() {
        for descriptor in [Descriptor::Unknown, Descriptor::ns, Descriptor::None] {
            assert_eq!(CipRule4a::ord(descriptor).unwrap(), 0);
        }
        for descriptor in [
            Descriptor::r,
            Descriptor::s,
            Descriptor::m,
            Descriptor::p,
            Descriptor::E,
            Descriptor::Z,
        ] {
            assert_eq!(CipRule4a::ord(descriptor).unwrap(), 1);
        }
        for descriptor in [
            Descriptor::R,
            Descriptor::S,
            Descriptor::M,
            Descriptor::P,
            Descriptor::seqTrans,
            Descriptor::seqCis,
        ] {
            assert_eq!(CipRule4a::ord(descriptor).unwrap(), 2);
        }
        assert!(matches!(
            CipRule4a::ord(Descriptor::SP_4),
            Err(CipLabelerError::InvalidStereoDescriptor)
        ));

        let molecule = Molecule::from_smiles("CC").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 0, false).unwrap();
        let root = digraph.get_current_root();
        let node_r = digraph
            .add_node(vec![], Some(1), RationalI32::new(6, 1), 2, 0)
            .unwrap();
        let node_s = digraph
            .add_node(vec![], Some(1), RationalI32::new(6, 1), 2, 0)
            .unwrap();
        digraph.add_edge(root, Some(0), node_r);
        let edge_r = CipEdgeId::new(digraph.edges.len() - 1);
        digraph.add_edge(root, Some(0), node_s);
        let edge_s = CipEdgeId::new(digraph.edges.len() - 1);

        digraph.edges[edge_r.index()].set_aux(Descriptor::R);
        digraph.edges[edge_s.index()].set_aux(Descriptor::r);
        let mut context = CipLabelerContext::new(0);
        assert_eq!(
            CipRule4a
                .compare(&mut digraph, &mut context, edge_r, edge_s)
                .unwrap(),
            1
        );

        digraph.edges[edge_r.index()].set_aux(Descriptor::None);
        digraph.edges[edge_s.index()].set_aux(Descriptor::None);
        digraph.nodes[node_r.index()].set_aux(Descriptor::R);
        digraph.nodes[node_s.index()].set_aux(Descriptor::r);
        assert_eq!(
            CipRule4a
                .compare(&mut digraph, &mut context, edge_r, edge_s)
                .unwrap(),
            1
        );
        digraph.nodes[node_r.index()].set_aux(Descriptor::SP_4);
        assert!(matches!(
            CipRule4a.compare(&mut digraph, &mut context, edge_r, edge_s),
            Err(CipLabelerError::InvalidStereoDescriptor)
        ));
    }

    #[test]
    fn ciplabeler_rules_rule4c_orders_lowercase_descriptor_classes_like_rdkit() {
        assert_eq!(CipRule4c::ord(Descriptor::m), 2);
        assert_eq!(CipRule4c::ord(Descriptor::r), 2);
        assert_eq!(CipRule4c::ord(Descriptor::p), 1);
        assert_eq!(CipRule4c::ord(Descriptor::s), 1);
        assert_eq!(CipRule4c::ord(Descriptor::R), 0);

        let molecule = Molecule::from_smiles("CC").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 0, false).unwrap();
        let root = digraph.get_current_root();
        let node_m = digraph
            .add_node(vec![], Some(1), RationalI32::new(6, 1), 2, 0)
            .unwrap();
        let node_p = digraph
            .add_node(vec![], Some(1), RationalI32::new(6, 1), 2, 0)
            .unwrap();
        digraph.add_edge(root, Some(0), node_m);
        let edge_m = CipEdgeId::new(digraph.edges.len() - 1);
        digraph.add_edge(root, Some(0), node_p);
        let edge_p = CipEdgeId::new(digraph.edges.len() - 1);

        digraph.edges[edge_m.index()].set_aux(Descriptor::m);
        digraph.edges[edge_p.index()].set_aux(Descriptor::p);
        let mut context = CipLabelerContext::new(0);
        assert_eq!(
            CipRule4c
                .compare(&mut digraph, &mut context, edge_m, edge_p)
                .unwrap(),
            1
        );

        digraph.edges[edge_m.index()].set_aux(Descriptor::None);
        digraph.edges[edge_p.index()].set_aux(Descriptor::None);
        digraph.nodes[node_m.index()].set_aux(Descriptor::r);
        digraph.nodes[node_p.index()].set_aux(Descriptor::s);
        assert_eq!(
            CipRule4c
                .compare(&mut digraph, &mut context, edge_m, edge_p)
                .unwrap(),
            1
        );
    }

    #[test]
    fn ciplabeler_rules_rule4b_nonroot_ref_branch_matches_rdkit() {
        let molecule = Molecule::from_smiles("CCC").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 1, false).unwrap();
        let root = digraph.get_current_root();
        let branch = digraph
            .add_node(vec![1, 1, 0], Some(0), RationalI32::new(6, 1), 2, 0)
            .unwrap();
        digraph.add_edge(root, Some(0), branch);
        let r_child = digraph
            .add_node(vec![], Some(0), RationalI32::new(6, 1), 3, 0)
            .unwrap();
        let s_child = digraph
            .add_node(vec![], Some(0), RationalI32::new(6, 1), 3, 0)
            .unwrap();
        digraph.nodes[r_child.index()].set_aux(Descriptor::R);
        digraph.nodes[s_child.index()].set_aux(Descriptor::S);
        digraph.add_edge(branch, Some(0), r_child);
        let r_edge = CipEdgeId::new(digraph.edges.len() - 1);
        digraph.add_edge(branch, Some(0), s_child);
        let s_edge = CipEdgeId::new(digraph.edges.len() - 1);

        let mut context = CipLabelerContext::new(0);
        assert_eq!(
            CipRule4b::with_ref(Descriptor::R)
                .compare(&mut digraph, &mut context, r_edge, s_edge)
                .unwrap(),
            1
        );
        assert_eq!(
            CipRule4b::with_ref(Descriptor::S)
                .compare(&mut digraph, &mut context, r_edge, s_edge)
                .unwrap(),
            -1
        );
        assert_eq!(
            CipRule4b::new()
                .compare(&mut digraph, &mut context, r_edge, s_edge)
                .unwrap(),
            0
        );
    }

    #[test]
    fn ciplabeler_rules_rule4b_reference_descriptor_search_matches_rdkit() {
        let molecule = Molecule::from_smiles("CCCC").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 1, false).unwrap();
        let root = digraph.get_current_root();
        let root_edges = digraph.node_edges(root).unwrap();
        let left = root_edges
            .iter()
            .copied()
            .find(|edge_id| {
                let end = digraph.edge(*edge_id).get_end();
                digraph.node(end).atom_idx() == Some(0)
            })
            .unwrap();
        let right = root_edges
            .iter()
            .copied()
            .find(|edge_id| {
                let end = digraph.edge(*edge_id).get_end();
                digraph.node(end).atom_idx() == Some(2)
            })
            .unwrap();
        let left_node = digraph.edge(left).get_end();
        let right_node = digraph.edge(right).get_end();
        digraph.nodes[left_node.index()].set_aux(Descriptor::R);
        digraph.nodes[right_node.index()].set_aux(Descriptor::S);

        let rule = CipRule4b::new();
        let mut context = CipLabelerContext::new(0);
        assert_eq!(
            rule.get_reference_descriptors(None, &mut digraph, &mut context, left_node)
                .unwrap(),
            vec![Descriptor::R]
        );
        assert_eq!(
            rule.get_reference_descriptors(None, &mut digraph, &mut context, right_node)
                .unwrap(),
            vec![Descriptor::S]
        );
        assert!(rule.has_descriptors(&mut digraph, left_node).unwrap());
    }

    #[test]
    fn ciplabeler_pairlist_rule4b_reference_counts_and_level_shape_match_rdkit() {
        let molecule = Molecule::from_smiles("C(F)(Cl)Br").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 0, false).unwrap();
        let mut nodes = Vec::new();
        for descriptor in [
            Descriptor::R,
            Descriptor::M,
            Descriptor::seqCis,
            Descriptor::S,
            Descriptor::P,
            Descriptor::seqTrans,
            Descriptor::r,
            Descriptor::None,
        ] {
            let node = digraph
                .add_node(vec![], Some(1), RationalI32::new(9, 1), 2, 0)
                .unwrap();
            digraph.nodes[node.index()].set_aux(descriptor);
            nodes.push(node);
        }

        let rule = CipRule4b::new();
        let mut result = Vec::new();
        assert!(rule.get_reference(&digraph, &nodes[..3], &mut result));
        assert_eq!(result, vec![Descriptor::R]);
        result.clear();
        assert!(rule.get_reference(&digraph, &nodes[3..6], &mut result));
        assert_eq!(result, vec![Descriptor::S]);
        result.clear();
        assert!(rule.get_reference(&digraph, &[nodes[0], nodes[3]], &mut result));
        assert_eq!(result, vec![Descriptor::R, Descriptor::S]);
        result.clear();
        assert!(!rule.get_reference(&digraph, &nodes[6..], &mut result));
        assert!(result.is_empty());

        let parent_a = digraph
            .add_node(
                vec![1, 0, 0, 0],
                Some(0),
                RationalI32::new(6, 1),
                2,
                CipNode::EXPANDED,
            )
            .unwrap();
        let parent_b = digraph
            .add_node(
                vec![1, 0, 0, 0],
                Some(0),
                RationalI32::new(6, 1),
                2,
                CipNode::EXPANDED,
            )
            .unwrap();
        let fluorine = digraph
            .add_node(vec![1, 1, 0, 0], Some(1), RationalI32::new(9, 1), 3, 0)
            .unwrap();
        let chlorine = digraph
            .add_node(vec![1, 0, 1, 0], Some(2), RationalI32::new(17, 1), 3, 0)
            .unwrap();
        digraph.add_edge(parent_a, Some(0), fluorine);
        digraph.add_edge(parent_b, Some(0), fluorine);
        digraph.add_edge(parent_b, Some(1), chlorine);

        let sort_rule = CipRule1a;
        let sort_rules = [&sort_rule as &dyn CipSequenceRule];
        let mut context = CipLabelerContext::new(0);
        assert!(matches!(
            rule.get_next_level(
                Some(&sort_rules),
                &mut digraph,
                &mut context,
                &[vec![parent_a, parent_b]],
            ),
            Err(CipLabelerError::SomethingUnexpected)
        ));
    }

    #[test]
    fn ciplabeler_rules_rule4b_root_pair_comparison_matches_rdkit() {
        let molecule = Molecule::from_smiles("CCCC").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 1, false).unwrap();
        let root = digraph.get_current_root();
        let root_edges = digraph.node_edges(root).unwrap();
        let left = root_edges
            .iter()
            .copied()
            .find(|edge_id| {
                let end = digraph.edge(*edge_id).get_end();
                digraph.node(end).atom_idx() == Some(0)
            })
            .unwrap();
        let right = root_edges
            .iter()
            .copied()
            .find(|edge_id| {
                let end = digraph.edge(*edge_id).get_end();
                digraph.node(end).atom_idx() == Some(2)
            })
            .unwrap();
        let left_node = digraph.edge(left).get_end();
        let right_node = digraph.edge(right).get_end();
        digraph.nodes[left_node.index()].set_aux(Descriptor::R);
        digraph.nodes[right_node.index()].set_aux(Descriptor::S);

        let rule = CipRule4b::new();
        let sort_rules = vec![&rule as &dyn CipSequenceRule];
        let mut context = CipLabelerContext::new(0);
        assert_eq!(
            rule.compare_with_sort_rules(
                Some(&sort_rules),
                &mut digraph,
                &mut context,
                left,
                right,
            )
            .unwrap(),
            0
        );
    }

    #[test]
    fn ciplabeler_rules_rule5new_rule6_nonroot_ref_branch_matches_rdkit() {
        let molecule = Molecule::from_smiles("CCC").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 1, false).unwrap();
        let root = digraph.get_current_root();
        let branch = digraph
            .add_node(vec![1, 1, 0], Some(0), RationalI32::new(6, 1), 2, 0)
            .unwrap();
        digraph.add_edge(root, Some(0), branch);
        let r_child = digraph
            .add_node(vec![], Some(0), RationalI32::new(6, 1), 3, 0)
            .unwrap();
        let s_child = digraph
            .add_node(vec![], Some(0), RationalI32::new(6, 1), 3, 0)
            .unwrap();
        digraph.nodes[r_child.index()].set_aux(Descriptor::R);
        digraph.nodes[s_child.index()].set_aux(Descriptor::S);
        digraph.add_edge(branch, Some(0), r_child);
        let r_edge = CipEdgeId::new(digraph.edges.len() - 1);
        digraph.add_edge(branch, Some(0), s_child);
        let s_edge = CipEdgeId::new(digraph.edges.len() - 1);

        let mut context = CipLabelerContext::new(0);
        assert_eq!(
            CipRule5New::with_ref(Descriptor::R)
                .compare(&mut digraph, &mut context, r_edge, s_edge)
                .unwrap(),
            1
        );
        assert_eq!(
            CipRule5New::with_ref(Descriptor::S)
                .compare(&mut digraph, &mut context, r_edge, s_edge)
                .unwrap(),
            -1
        );
        assert_eq!(
            CipRule5New::new()
                .compare(&mut digraph, &mut context, r_edge, s_edge)
                .unwrap(),
            0
        );
    }

    #[test]
    fn ciplabeler_rules_rule5new_rule6_root_pair_comparison_matches_rdkit() {
        let molecule = Molecule::from_smiles("CCCC").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 1, false).unwrap();
        let root = digraph.get_current_root();
        let root_edges = digraph.node_edges(root).unwrap();
        let left = root_edges
            .iter()
            .copied()
            .find(|edge_id| {
                let end = digraph.edge(*edge_id).get_end();
                digraph.node(end).atom_idx() == Some(0)
            })
            .unwrap();
        let right = root_edges
            .iter()
            .copied()
            .find(|edge_id| {
                let end = digraph.edge(*edge_id).get_end();
                digraph.node(end).atom_idx() == Some(2)
            })
            .unwrap();
        let left_node = digraph.edge(left).get_end();
        let right_node = digraph.edge(right).get_end();
        digraph.nodes[left_node.index()].set_aux(Descriptor::R);
        digraph.nodes[right_node.index()].set_aux(Descriptor::S);

        let rule = CipRule5New::new();
        let sort_rules = vec![&rule as &dyn CipSequenceRule];
        let mut context = CipLabelerContext::new(0);
        assert_eq!(
            rule.compare_with_sort_rules(
                Some(&sort_rules),
                &mut digraph,
                &mut context,
                left,
                right,
            )
            .unwrap(),
            2
        );
        assert_eq!(
            rule.compare_with_sort_rules(
                Some(&sort_rules),
                &mut digraph,
                &mut context,
                right,
                left,
            )
            .unwrap(),
            -2
        );

        let unrelated_rule = CipRule1a;
        let unrelated_rules = [&unrelated_rule as &dyn CipSequenceRule];
        let replacement = CipRule5New::with_ref(Descriptor::R);
        assert!(matches!(
            rule.get_ref_sorter(Some(&unrelated_rules), &replacement),
            Err(CipLabelerError::Rule5NewInstanceNotInRuleSet)
        ));
    }

    #[test]
    fn ciplabeler_rules_rule5new_rule6_rule6_ref_atom_matches_rdkit() {
        let molecule = Molecule::from_smiles("CCC").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 1, false).unwrap();
        let root = digraph.get_current_root();
        let root_edges = digraph.node_edges(root).unwrap();
        let left = root_edges
            .iter()
            .copied()
            .find(|edge_id| {
                let end = digraph.edge(*edge_id).get_end();
                digraph.node(end).atom_idx() == Some(0)
            })
            .unwrap();
        let right = root_edges
            .iter()
            .copied()
            .find(|edge_id| {
                let end = digraph.edge(*edge_id).get_end();
                digraph.node(end).atom_idx() == Some(2)
            })
            .unwrap();

        let mut context = CipLabelerContext::new(0);
        assert_eq!(
            CipRule6
                .compare(&mut digraph, &mut context, left, right)
                .unwrap(),
            0
        );
        digraph.set_rule6_ref(Some(0)).unwrap();
        assert_eq!(
            CipRule6
                .compare(&mut digraph, &mut context, left, right)
                .unwrap(),
            1
        );
        digraph.set_rule6_ref(Some(2)).unwrap();
        assert_eq!(
            CipRule6
                .compare(&mut digraph, &mut context, left, right)
                .unwrap(),
            -1
        );
        assert_eq!(
            CipRule6
                .compare(&mut digraph, &mut context, right, right)
                .unwrap(),
            0
        );
    }

    #[test]
    fn ciplabeler_rules_rule5new_rule6_modern_dispatch_uses_nine_pinned_rules() {
        assert_eq!(cip_all_rules().unwrap().get_num_sub_rules(), 9);
    }

    #[test]
    fn ciplabeler_configuration_constructors_accessors_and_carriers_match_rdkit() {
        let molecule = Molecule::from_smiles("CCO").unwrap();
        let mut config = CipConfiguration::new(&molecule, 1).unwrap();
        assert_eq!(CipConfiguration::IMPLICIT_H, u32::MAX);
        assert_eq!(config.get_focus(), 1);
        assert_eq!(config.get_foci(), &[1]);
        assert!(config.get_carriers().is_empty());
        assert_eq!(config.get_digraph().get_original_root(), CipNodeId::new(0));
        assert_eq!(
            config.get_digraph().node(CipNodeId::new(0)).atom_idx(),
            Some(1)
        );

        config.set_carriers(vec![Some(0), Some(2), None]);
        assert_eq!(config.get_carriers(), &[Some(0), Some(2), None]);

        let mut multi = CipConfiguration::with_foci(&molecule, vec![0, 1], true).unwrap();
        assert_eq!(multi.get_focus(), 0);
        assert_eq!(multi.get_foci(), &[0, 1]);
        assert_eq!(multi.get_digraph().get_original_root(), CipNodeId::new(0));
        assert!(matches!(
            CipConfiguration::with_foci(&molecule, Vec::new(), false),
            Err(CipLabelerError::EmptyConfigurationFoci)
        ));
    }

    #[test]
    fn ciplabeler_configuration_parity4_matches_rdkit_table() {
        let reference = [0, 1, 2, 3];
        let even = [
            [0, 1, 2, 3],
            [0, 2, 3, 1],
            [0, 3, 1, 2],
            [1, 0, 3, 2],
            [1, 2, 0, 3],
            [1, 3, 2, 0],
            [2, 0, 1, 3],
            [2, 1, 3, 0],
            [2, 3, 0, 1],
            [3, 0, 2, 1],
            [3, 1, 0, 2],
            [3, 2, 1, 0],
        ];
        let odd = [
            [0, 1, 3, 2],
            [0, 2, 1, 3],
            [0, 3, 2, 1],
            [1, 0, 2, 3],
            [1, 2, 3, 0],
            [1, 3, 0, 2],
            [2, 0, 3, 1],
            [2, 1, 0, 3],
            [2, 3, 1, 0],
            [3, 0, 1, 2],
            [3, 1, 2, 0],
            [3, 2, 0, 1],
        ];

        for target in even {
            assert_eq!(
                CipConfiguration::parity4(&target, &reference).unwrap(),
                2,
                "target={target:?}"
            );
        }
        for target in odd {
            assert_eq!(
                CipConfiguration::parity4(&target, &reference).unwrap(),
                1,
                "target={target:?}"
            );
        }
        assert_eq!(
            CipConfiguration::parity4(&[9, 8, 7, 6], &reference).unwrap(),
            0
        );
        assert!(matches!(
            CipConfiguration::parity4(&[0, 1, 2], &reference),
            Err(CipLabelerError::ParityVectorsMustHaveSize4)
        ));
        assert!(matches!(
            CipConfiguration::parity4(&reference, &[0, 1, 2]),
            Err(CipLabelerError::ParityVectorsMustHaveSize4)
        ));
    }

    #[test]
    fn ciplabeler_configuration_base_label_returns_unknown_like_rdkit() {
        let molecule = Molecule::from_smiles("CCO").unwrap();
        let config = CipConfiguration::new(&molecule, 1).unwrap();
        let mut digraph = CipDigraph::new(&molecule, 1, false).unwrap();
        let node = digraph.get_current_root();
        let rules = CipRules::new(vec![Box::new(CipRule1a)]).unwrap();
        assert_eq!(
            config.label(node, &mut digraph, &rules),
            Descriptor::Unknown
        );
    }

    fn tetrahedral_molecule(tag: ChiralTag) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(AtomSpec::new(Element::C).with_chiral_tag(tag));
        let fluorine = builder.add_atom(AtomSpec::new(Element::F));
        let chlorine = builder.add_atom(AtomSpec::new(Element::CL));
        let bromine = builder.add_atom(AtomSpec::new(Element::BR));
        let iodine = builder.add_atom(AtomSpec::new(Element::I));
        builder
            .add_bond(BondSpec::new(center, fluorine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, chlorine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, bromine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, iodine, BondOrder::Single))
            .unwrap();
        builder.build().unwrap()
    }

    #[test]
    fn ciplabeler_tetrahedral_constructor_builds_carriers_like_rdkit() {
        let molecule = tetrahedral_molecule(ChiralTag::TetrahedralCcw);
        let config = CipTetrahedral::new(&molecule, 0).unwrap();
        assert_eq!(config.get_focus(), 0);
        assert_eq!(config.get_carriers(), &[Some(1), Some(2), Some(3), Some(4)]);

        let implicit_h = Molecule::from_smiles_with_sanitize("C[C@H](F)Cl", false).unwrap();
        let implicit_config = CipTetrahedral::new(&implicit_h, 1).unwrap();
        assert_eq!(
            implicit_config.get_carriers(),
            &[Some(0), Some(2), Some(3), Some(1)]
        );

        let mut builder = MoleculeBuilder::new();
        let center =
            builder.add_atom(AtomSpec::new(Element::N).with_chiral_tag(ChiralTag::TetrahedralCcw));
        let carbon = builder.add_atom(AtomSpec::new(Element::C));
        let oxygen = builder.add_atom(AtomSpec::new(Element::O));
        builder
            .add_bond(BondSpec::new(center, carbon, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(center, oxygen, BondOrder::Single))
            .unwrap();
        let trigonal_pyramid = builder.build().unwrap();
        let trigonal_config = CipTetrahedral::new(&trigonal_pyramid, center.index()).unwrap();
        assert_eq!(
            trigonal_config.get_carriers(),
            &[
                Some(carbon.index()),
                Some(oxygen.index()),
                Some(center.index()),
                None
            ]
        );

        let methane = Molecule::from_smiles("C").unwrap();
        assert!(matches!(
            CipTetrahedral::new(&methane, 0),
            Err(CipLabelerError::BadTetrahedralConfig)
        ));
    }

    #[test]
    fn ciplabeler_tetrahedral_primary_label_records_atom_props_like_rdkit() {
        let molecule = tetrahedral_molecule(ChiralTag::TetrahedralCcw);
        let mut config = CipTetrahedral::new(&molecule, 0).unwrap();
        assert!(!config.has_primary_label());

        config.ranked_anchors = vec![4, 3, 2, 1];
        for (descriptor, code) in [
            (Descriptor::R, "R"),
            (Descriptor::S, "S"),
            (Descriptor::r, "r"),
            (Descriptor::s, "s"),
        ] {
            config.set_primary_label(descriptor).unwrap();
            assert!(config.has_primary_label());
            assert_eq!(
                config.primary_label(),
                Some(&super::CipAtomPrimaryLabel {
                    atom_idx: 0,
                    cip_code: code,
                    cip_neighbor_order: vec![4, 3, 2, 1],
                })
            );
            config.reset_primary_label();
            assert!(!config.has_primary_label());
        }
        assert!(matches!(
            config.set_primary_label(Descriptor::E),
            Err(CipLabelerError::DescriptorNotSupportedForAtoms)
        ));
        assert!(matches!(
            config.set_primary_label(Descriptor::Unknown),
            Err(CipLabelerError::InvalidAtomDescriptor)
        ));
    }

    #[test]
    fn ciplabeler_tetrahedral_label_assigns_r_s_and_ranked_anchors_like_rdkit() {
        let rules = cip_all_rules().unwrap();

        let molecule = tetrahedral_molecule(ChiralTag::TetrahedralCcw);
        let mut config = CipTetrahedral::new(&molecule, 0).unwrap();
        let mut context = CipLabelerContext::new(0);
        assert_eq!(config.label(&rules, &mut context).unwrap(), Descriptor::S);
        assert_eq!(config.ranked_anchors(), &[4, 3, 2, 1]);

        let molecule = tetrahedral_molecule(ChiralTag::TetrahedralCw);
        let mut config = CipTetrahedral::new(&molecule, 0).unwrap();
        let mut context = CipLabelerContext::new(0);
        assert_eq!(config.label(&rules, &mut context).unwrap(), Descriptor::R);
        assert_eq!(config.ranked_anchors(), &[4, 3, 2, 1]);
    }

    #[test]
    fn ciplabeler_tetrahedral_label_assigns_lowercase_for_pseudoasymmetric_priority() {
        let rules = CipRules::new(vec![Box::new(SinglePseudoAsymSortRule)]).unwrap();

        let molecule = tetrahedral_molecule(ChiralTag::TetrahedralCcw);
        let mut config = CipTetrahedral::new(&molecule, 0).unwrap();
        let mut context = CipLabelerContext::new(0);
        assert_eq!(config.label(&rules, &mut context).unwrap(), Descriptor::s);

        let molecule = tetrahedral_molecule(ChiralTag::TetrahedralCw);
        let mut config = CipTetrahedral::new(&molecule, 0).unwrap();
        let mut context = CipLabelerContext::new(0);
        assert_eq!(config.label(&rules, &mut context).unwrap(), Descriptor::r);
    }

    #[test]
    fn ciplabeler_tetrahedral_configuration_edge_filters_match_rdkit() {
        let molecule = Molecule::from_smiles("CC").unwrap();
        let mut digraph = CipDigraph::new(&molecule, 0, false).unwrap();
        let root = digraph.get_current_root();
        let edges = digraph.node_edges(root).unwrap();
        let internal = CipConfiguration::find_internal_edge(&digraph, &edges, 0, 1).unwrap();
        assert!(CipConfiguration::is_internal_edge(&digraph, internal, 0, 1));
        assert!(CipConfiguration::is_internal_edge(&digraph, internal, 1, 0));

        let mut without_internal = edges.clone();
        CipConfiguration::remove_internal_edges(&digraph, &mut without_internal, 0, 1);
        assert!(!without_internal.contains(&internal));

        let mut without_duplicates_and_hs = edges;
        CipConfiguration::remove_duplicates_and_hs(&digraph, &mut without_duplicates_and_hs);
        assert_eq!(without_duplicates_and_hs, vec![internal]);
    }

    #[test]
    fn ciplabeler_tetrahedral_label_external_digraph_changes_root_like_rdkit() {
        let molecule = tetrahedral_molecule(ChiralTag::TetrahedralCcw);
        let mut config = CipTetrahedral::new(&molecule, 0).unwrap();
        let mut external = CipDigraph::new(&molecule, 0, false).unwrap();
        let root = external.get_original_root();
        let rules = cip_all_rules().unwrap();
        let mut context = CipLabelerContext::new(0);

        assert_eq!(
            config
                .label_with_external_digraph(root, &mut external, &rules, &mut context)
                .unwrap(),
            Descriptor::S
        );
        assert_eq!(external.get_current_root(), root);
    }

    #[test]
    fn ciplabeler_tetrahedral_label_returns_ns_or_unknown_for_unresolved_cases_like_rdkit() {
        let molecule = Molecule::from_smiles_with_sanitize("C[C@H](C)C", false).unwrap();
        let mut config = CipTetrahedral::new(&molecule, 1).unwrap();
        let rules = cip_all_rules().unwrap();
        let mut context = CipLabelerContext::new(0);
        assert_eq!(
            config.label(&rules, &mut context).unwrap(),
            Descriptor::Unknown
        );

        let leaf_molecule = tetrahedral_molecule(ChiralTag::TetrahedralCcw);
        let mut config = CipTetrahedral::new(&leaf_molecule, 0).unwrap();
        let mut external = CipDigraph::new(&leaf_molecule, 0, false).unwrap();
        let root = external.get_original_root();
        let leaf_edge = external.node_edges(root).unwrap()[0];
        let leaf = external.edge(leaf_edge).get_end();
        let mut context = CipLabelerContext::new(0);
        assert_eq!(
            config
                .label_with_external_digraph(leaf, &mut external, &rules, &mut context)
                .unwrap(),
            Descriptor::ns
        );
    }

    fn sp2_bond_molecule(stereo: BondStereo, stereo_atoms: (usize, usize)) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let fluorine = builder.add_atom(AtomSpec::new(Element::F));
        let begin = builder.add_atom(AtomSpec::new(Element::C));
        let end = builder.add_atom(AtomSpec::new(Element::C));
        let chlorine = builder.add_atom(AtomSpec::new(Element::CL));
        let bromine = builder.add_atom(AtomSpec::new(Element::BR));
        let iodine = builder.add_atom(AtomSpec::new(Element::I));
        builder
            .add_bond(
                BondSpec::new(begin, end, BondOrder::Double)
                    .with_stereo(stereo)
                    .with_stereo_atoms(
                        match stereo_atoms.0 {
                            0 => fluorine,
                            3 => chlorine,
                            4 => bromine,
                            5 => iodine,
                            _ => panic!("unexpected stereo atom fixture index"),
                        },
                        match stereo_atoms.1 {
                            0 => fluorine,
                            3 => chlorine,
                            4 => bromine,
                            5 => iodine,
                            _ => panic!("unexpected stereo atom fixture index"),
                        },
                    ),
            )
            .unwrap();
        builder
            .add_bond(BondSpec::new(begin, fluorine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(begin, bromine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(end, chlorine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(end, iodine, BondOrder::Single))
            .unwrap();
        builder.build().unwrap()
    }

    fn sp2_ez_bond_molecule_without_explicit_stereo_atoms() -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let fluorine = builder.add_atom(AtomSpec::new(Element::F).with_prop("_CIPRank", "10"));
        let begin = builder.add_atom(AtomSpec::new(Element::C));
        let end = builder.add_atom(AtomSpec::new(Element::C));
        let chlorine = builder.add_atom(AtomSpec::new(Element::CL).with_prop("_CIPRank", "20"));
        let bromine = builder.add_atom(AtomSpec::new(Element::BR).with_prop("_CIPRank", "30"));
        let iodine = builder.add_atom(AtomSpec::new(Element::I).with_prop("_CIPRank", "40"));
        builder
            .add_bond(BondSpec::new(begin, end, BondOrder::Double).with_stereo(BondStereo::E))
            .unwrap();
        builder
            .add_bond(BondSpec::new(begin, fluorine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(begin, bromine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(end, chlorine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(end, iodine, BondOrder::Single))
            .unwrap();
        builder.build().unwrap()
    }

    fn sp2_pseudoasymmetric_bond_molecule(stereo: BondStereo) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let bromine = builder.add_atom(AtomSpec::new(Element::BR));
        let begin = builder.add_atom(AtomSpec::new(Element::C));
        let end = builder.add_atom(AtomSpec::new(Element::C));
        let fluorine = builder.add_atom(AtomSpec::new(Element::F));
        let iodine = builder.add_atom(AtomSpec::new(Element::I));
        let chlorine = builder.add_atom(AtomSpec::new(Element::CL));
        builder
            .add_bond(
                BondSpec::new(begin, end, BondOrder::Double)
                    .with_stereo(stereo)
                    .with_stereo_atoms(iodine, chlorine),
            )
            .unwrap();
        builder
            .add_bond(BondSpec::new(begin, bromine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(begin, iodine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(end, fluorine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(end, chlorine, BondOrder::Single))
            .unwrap();
        builder.build().unwrap()
    }

    fn cip_rank_neighbor_molecule(ranks: &[Option<&str>]) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(AtomSpec::new(Element::C));
        let skipped = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(center, skipped, BondOrder::Single))
            .unwrap();
        for rank in ranks {
            let mut spec = AtomSpec::new(Element::F);
            if let Some(rank) = rank {
                spec = spec.with_prop("_CIPRank", *rank);
            }
            let neighbor = builder.add_atom(spec);
            builder
                .add_bond(BondSpec::new(center, neighbor, BondOrder::Single))
                .unwrap();
        }
        builder.build().unwrap()
    }

    #[test]
    fn ciplabeler_sp2bond_constructor_builds_carriers_like_rdkit_find_stereo_atoms() {
        let molecule = sp2_bond_molecule(BondStereo::Cis, (4, 5));
        let config = CipSp2Bond::new(&molecule, 0, 1, 2, BondStereo::Cis).unwrap();
        assert_eq!(config.get_foci(), &[1, 2]);
        assert_eq!(config.get_carriers(), &[Some(4), Some(5)]);

        assert!(matches!(
            CipSp2Bond::new(&molecule, 0, 1, 2, BondStereo::E),
            Err(CipLabelerError::BadSp2BondConfig)
        ));
        assert!(matches!(
            CipSp2Bond::new(&molecule, 1, 1, 0, BondStereo::Cis),
            Err(CipLabelerError::Sp2BondNotDoubleBond)
        ));
        assert!(matches!(
            CipSp2Bond::new(&molecule, 1, 1, 2, BondStereo::Cis),
            Err(CipLabelerError::BadSp2BondFoci)
        ));

        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(a0, a1, BondOrder::Double).with_stereo(BondStereo::E))
            .unwrap();
        let missing = builder.build().unwrap();
        assert!(matches!(
            CipSp2Bond::new(&missing, 0, 0, 1, BondStereo::Cis),
            Err(CipLabelerError::IncorrectNumberOfStereoAtoms)
        ));
    }

    #[test]
    fn ciplabeler_sp2bond_constructor_finds_highest_cip_neighbors_like_rdkit() {
        let molecule = sp2_ez_bond_molecule_without_explicit_stereo_atoms();
        let config = CipSp2Bond::new(&molecule, 0, 1, 2, BondStereo::Trans).unwrap();
        assert_eq!(config.get_foci(), &[1, 2]);
        assert_eq!(config.get_carriers(), &[Some(4), Some(5)]);
    }

    #[test]
    fn ciplabeler_sp2bond_find_highest_cip_neighbor_matches_missing_tie_and_recovery() {
        let missing = cip_rank_neighbor_molecule(&[Some("10"), None, Some("20")]);
        assert_eq!(
            CipSp2Bond::find_highest_cip_neighbor_like_rdkit(&CipMol::new(&missing), 0, 1).unwrap(),
            None
        );

        let tied = cip_rank_neighbor_molecule(&[Some("10"), Some("10")]);
        assert_eq!(
            CipSp2Bond::find_highest_cip_neighbor_like_rdkit(&CipMol::new(&tied), 0, 1).unwrap(),
            None
        );

        let recovered = cip_rank_neighbor_molecule(&[Some("10"), Some("10"), Some("20")]);
        assert_eq!(
            CipSp2Bond::find_highest_cip_neighbor_like_rdkit(&CipMol::new(&recovered), 0, 1)
                .unwrap(),
            Some(4)
        );
    }

    #[test]
    fn ciplabeler_sp2bond_rejects_undefined_source_stereo_before_explicit_stereo_atoms() {
        for stereo in [BondStereo::None, BondStereo::Any] {
            let mut builder = MoleculeBuilder::new();
            let left = builder.add_atom(AtomSpec::new(Element::F));
            let begin = builder.add_atom(AtomSpec::new(Element::C));
            let end = builder.add_atom(AtomSpec::new(Element::C));
            let right = builder.add_atom(AtomSpec::new(Element::CL));
            builder
                .add_bond(
                    BondSpec::new(begin, end, BondOrder::Double)
                        .with_stereo(stereo)
                        .with_stereo_atoms(left, right),
                )
                .unwrap();
            builder
                .add_bond(BondSpec::new(begin, left, BondOrder::Single))
                .unwrap();
            builder
                .add_bond(BondSpec::new(end, right, BondOrder::Single))
                .unwrap();
            let molecule = builder.build().unwrap();

            assert!(matches!(
                CipSp2Bond::new(&molecule, 0, 1, 2, BondStereo::Cis),
                Err(CipLabelerError::Sp2BondHasNoDefinedStereo)
            ));
        }
    }

    #[test]
    fn ciplabeler_sp2bond_primary_label_records_bond_props_like_rdkit() {
        let molecule = sp2_bond_molecule(BondStereo::Cis, (4, 5));
        for (descriptor, code) in [
            (Descriptor::seqTrans, "e"),
            (Descriptor::E, "E"),
            (Descriptor::seqCis, "z"),
            (Descriptor::Z, "Z"),
        ] {
            let mut config = CipSp2Bond::new(&molecule, 0, 1, 2, BondStereo::Cis).unwrap();
            assert!(!config.has_primary_label());
            config.ranked_anchors = vec![4, 5];
            config.set_primary_label(descriptor).unwrap();
            assert!(config.has_primary_label());
            assert_eq!(
                config.primary_label(),
                Some(&super::CipBondPrimaryLabel {
                    bond_idx: 0,
                    stereo_atoms: [4, 5],
                    stereo: BondStereo::Cis,
                    cip_code: code,
                    cip_neighbor_order: vec![4, 5],
                })
            );
            config.reset_primary_label();
            assert!(!config.has_primary_label());
        }

        let mut config = CipSp2Bond::new(&molecule, 0, 1, 2, BondStereo::Cis).unwrap();
        assert!(matches!(
            config.set_primary_label(Descriptor::R),
            Err(CipLabelerError::DescriptorNotSupportedForDoubleBonds)
        ));
        assert!(matches!(
            config.set_primary_label(Descriptor::Unknown),
            Err(CipLabelerError::InvalidBondDescriptor)
        ));
    }

    #[test]
    fn ciplabeler_sp2bond_label_assigns_e_z_and_ranked_anchors_like_rdkit() {
        let rules = cip_all_rules().unwrap();

        let molecule = sp2_bond_molecule(BondStereo::Cis, (4, 5));
        let mut config = CipSp2Bond::new(&molecule, 0, 1, 2, BondStereo::Cis).unwrap();
        let mut context = CipLabelerContext::new(0);
        assert_eq!(config.label(&rules, &mut context).unwrap(), Descriptor::Z);
        assert_eq!(config.ranked_anchors(), &[4, 5]);

        let molecule = sp2_bond_molecule(BondStereo::Trans, (4, 5));
        let mut config = CipSp2Bond::new(&molecule, 0, 1, 2, BondStereo::Trans).unwrap();
        let mut context = CipLabelerContext::new(0);
        assert_eq!(config.label(&rules, &mut context).unwrap(), Descriptor::E);
        assert_eq!(config.ranked_anchors(), &[4, 5]);
    }

    #[test]
    fn ciplabeler_sp2bond_label_flips_config_when_carrier_is_not_top_priority_like_rdkit() {
        let rules = cip_all_rules().unwrap();

        let molecule = sp2_bond_molecule(BondStereo::Cis, (0, 5));
        let mut config = CipSp2Bond::new(&molecule, 0, 1, 2, BondStereo::Cis).unwrap();
        let mut context = CipLabelerContext::new(0);
        assert_eq!(config.label(&rules, &mut context).unwrap(), Descriptor::E);
        assert_eq!(config.ranked_anchors(), &[4, 5]);

        let molecule = sp2_bond_molecule(BondStereo::Trans, (0, 5));
        let mut config = CipSp2Bond::new(&molecule, 0, 1, 2, BondStereo::Trans).unwrap();
        let mut context = CipLabelerContext::new(0);
        assert_eq!(config.label(&rules, &mut context).unwrap(), Descriptor::Z);
        assert_eq!(config.ranked_anchors(), &[4, 5]);
    }

    #[test]
    fn ciplabeler_sp2bond_label_assigns_pseudoasymmetric_seq_descriptors_like_rdkit() {
        let rules = CipRules::new(vec![Box::new(SinglePseudoAsymSortRule)]).unwrap();
        for (stereo, expected) in [
            (BondStereo::Cis, Descriptor::seqCis),
            (BondStereo::Trans, Descriptor::seqTrans),
        ] {
            let molecule = sp2_pseudoasymmetric_bond_molecule(stereo);
            let mut config = CipSp2Bond::new(&molecule, 0, 1, 2, stereo).unwrap();
            let mut context = CipLabelerContext::new(0);
            assert_eq!(config.label(&rules, &mut context).unwrap(), expected);
            assert_eq!(config.ranked_anchors(), &[4, 5]);
        }
    }

    #[test]
    fn ciplabeler_sp2bond_label_returns_unknown_for_nonunique_substituents() {
        let molecule = sp2_bond_molecule(BondStereo::Cis, (4, 5));
        let mut config = CipSp2Bond::new(&molecule, 0, 1, 2, BondStereo::Cis).unwrap();
        let rules = CipRules::new(vec![Box::new(AlwaysEqualSortRule)]).unwrap();
        let mut context = CipLabelerContext::new(0);

        assert_eq!(
            config.label(&rules, &mut context).unwrap(),
            Descriptor::Unknown
        );
        assert!(config.ranked_anchors().is_empty());
    }

    #[test]
    fn ciplabeler_sp2bond_label_external_digraph_matches_rdkit_overload() {
        let molecule = sp2_bond_molecule(BondStereo::Cis, (4, 5));
        let mut config = CipSp2Bond::new(&molecule, 0, 1, 2, BondStereo::Cis).unwrap();
        let mut external = CipDigraph::new(&molecule, 1, false).unwrap();
        let root = external.get_original_root();
        let rules = cip_all_rules().unwrap();
        let mut context = CipLabelerContext::new(0);

        assert_eq!(
            config
                .label_with_external_digraph(root, &mut external, &rules, &mut context)
                .unwrap(),
            Descriptor::Z
        );
        assert_eq!(config.ranked_anchors(), &[4, 5]);
    }

    fn atropisomer_bond_molecule(stereo: BondStereo, low_index_carriers: bool) -> Molecule {
        atropisomer_bond_molecule_with_carrier_positions(
            stereo,
            low_index_carriers,
            low_index_carriers,
        )
    }

    fn atropisomer_bond_molecule_with_carrier_positions(
        stereo: BondStereo,
        begin_low_index_carrier: bool,
        end_low_index_carrier: bool,
    ) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let fluorine = builder.add_atom(AtomSpec::new(Element::F));
        let begin = builder.add_atom(AtomSpec::new(Element::C));
        let end = builder.add_atom(AtomSpec::new(Element::C));
        let chlorine = builder.add_atom(AtomSpec::new(Element::CL));
        let bromine = builder.add_atom(AtomSpec::new(Element::BR));
        let iodine = builder.add_atom(AtomSpec::new(Element::I));
        builder
            .add_bond(BondSpec::new(begin, end, BondOrder::Single).with_stereo(stereo))
            .unwrap();
        if begin_low_index_carrier {
            builder
                .add_bond(BondSpec::new(begin, fluorine, BondOrder::Single))
                .unwrap();
            builder
                .add_bond(BondSpec::new(begin, bromine, BondOrder::Single))
                .unwrap();
        } else {
            builder
                .add_bond(BondSpec::new(begin, bromine, BondOrder::Single))
                .unwrap();
            builder
                .add_bond(BondSpec::new(begin, fluorine, BondOrder::Single))
                .unwrap();
        }
        if end_low_index_carrier {
            builder
                .add_bond(BondSpec::new(end, chlorine, BondOrder::Single))
                .unwrap();
            builder
                .add_bond(BondSpec::new(end, iodine, BondOrder::Single))
                .unwrap();
        } else {
            builder
                .add_bond(BondSpec::new(end, iodine, BondOrder::Single))
                .unwrap();
            builder
                .add_bond(BondSpec::new(end, chlorine, BondOrder::Single))
                .unwrap();
        }
        builder.build().unwrap()
    }

    fn atropisomer_bond_molecule_one_side_second_carrier(stereo: BondStereo) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let bromine = builder.add_atom(AtomSpec::new(Element::BR));
        let begin = builder.add_atom(AtomSpec::new(Element::C));
        let end = builder.add_atom(AtomSpec::new(Element::C));
        let chlorine = builder.add_atom(AtomSpec::new(Element::CL));
        let fluorine = builder.add_atom(AtomSpec::new(Element::F));
        let iodine = builder.add_atom(AtomSpec::new(Element::I));
        builder
            .add_bond(BondSpec::new(begin, end, BondOrder::Single).with_stereo(stereo))
            .unwrap();
        builder
            .add_bond(BondSpec::new(begin, bromine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(begin, fluorine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(end, chlorine, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(end, iodine, BondOrder::Single))
            .unwrap();
        builder.build().unwrap()
    }

    #[test]
    fn ciplabeler_atropisomerbond_constructor_builds_carriers_like_rdkit() {
        let molecule = atropisomer_bond_molecule(BondStereo::AtropCcw, true);
        let config = CipAtropisomerBond::new(&molecule, 0, 1, 2, BondStereo::AtropCcw).unwrap();
        assert_eq!(config.get_foci(), &[1, 2]);
        assert_eq!(config.get_carriers(), &[Some(0), Some(3)]);

        assert!(matches!(
            CipAtropisomerBond::new(&molecule, 0, 1, 2, BondStereo::Cis),
            Err(CipLabelerError::BadAtropisomerBondConfig)
        ));
        assert!(matches!(
            CipAtropisomerBond::new(&molecule, 1, 1, 2, BondStereo::AtropCcw),
            Err(CipLabelerError::BadAtropisomerBondFoci)
        ));

        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(a0, a1, BondOrder::Single).with_stereo(BondStereo::AtropCcw))
            .unwrap();
        let missing = builder.build().unwrap();
        let config = CipAtropisomerBond::new(&missing, 0, 0, 1, BondStereo::AtropCcw).unwrap();
        assert!(config.get_carriers().is_empty());
    }

    #[test]
    fn ciplabeler_atropisomerbond_primary_label_records_bond_props_like_rdkit() {
        let molecule = atropisomer_bond_molecule(BondStereo::AtropCcw, true);
        let mut config = CipAtropisomerBond::new(&molecule, 0, 1, 2, BondStereo::AtropCcw).unwrap();
        assert!(!config.has_primary_label());

        for (descriptor, code) in [
            (Descriptor::M, "M"),
            (Descriptor::P, "P"),
            (Descriptor::m, "m"),
            (Descriptor::p, "p"),
        ] {
            config.ranked_anchors = vec![4, 5];
            config.set_primary_label(descriptor).unwrap();
            assert!(config.has_primary_label());
            assert_eq!(
                config.primary_label(),
                Some(&super::CipAtropisomerBondPrimaryLabel {
                    bond_idx: 0,
                    cip_code: code,
                    cip_neighbor_order: vec![4, 5],
                })
            );
            config.reset_primary_label();
            assert!(!config.has_primary_label());
        }

        for descriptor in [
            Descriptor::R,
            Descriptor::S,
            Descriptor::r,
            Descriptor::s,
            Descriptor::SP_4,
            Descriptor::TBPY_5,
            Descriptor::OC_6,
            Descriptor::seqTrans,
            Descriptor::E,
            Descriptor::seqCis,
            Descriptor::Z,
        ] {
            assert!(matches!(
                config.set_primary_label(descriptor),
                Err(CipLabelerError::DescriptorNotSupportedForAtropisomerBonds)
            ));
        }
        for descriptor in [Descriptor::None, Descriptor::Unknown, Descriptor::ns] {
            assert!(matches!(
                config.set_primary_label(descriptor),
                Err(CipLabelerError::InvalidBondDescriptor)
            ));
        }
    }

    #[test]
    fn ciplabeler_atropisomerbond_reset_hides_a_source_primary_label_like_rdkit() {
        let molecule = atropisomer_bond_molecule(BondStereo::AtropCcw, true);
        let mut builder = molecule.to_builder();
        builder
            .bond_mut(BondId::new(0))
            .unwrap()
            .set_prop("_CIPCode", "P");
        let molecule = builder.build().unwrap();
        let mut config = CipAtropisomerBond::new(&molecule, 0, 1, 2, BondStereo::AtropCcw).unwrap();

        assert!(config.has_primary_label());
        config.reset_primary_label();
        assert!(!config.has_primary_label());
    }

    #[test]
    fn ciplabeler_atropisomerbond_label_assigns_m_p_and_ranked_anchors_like_rdkit() {
        let rules = cip_all_rules().unwrap();

        let molecule = atropisomer_bond_molecule(BondStereo::AtropCcw, false);
        let mut config = CipAtropisomerBond::new(&molecule, 0, 1, 2, BondStereo::AtropCcw).unwrap();
        let mut context = CipLabelerContext::new(0);
        assert_eq!(config.label(&rules, &mut context).unwrap(), Descriptor::M);
        assert_eq!(config.ranked_anchors(), &[4, 5]);

        let molecule = atropisomer_bond_molecule(BondStereo::AtropCw, false);
        let mut config = CipAtropisomerBond::new(&molecule, 0, 1, 2, BondStereo::AtropCw).unwrap();
        let mut context = CipLabelerContext::new(0);
        assert_eq!(config.label(&rules, &mut context).unwrap(), Descriptor::P);
        assert_eq!(config.ranked_anchors(), &[4, 5]);
    }

    #[test]
    fn ciplabeler_atropisomerbond_label_flips_config_when_carrier_is_second_like_rdkit() {
        let rules = cip_all_rules().unwrap();

        let molecule = atropisomer_bond_molecule_one_side_second_carrier(BondStereo::AtropCcw);
        let mut config = CipAtropisomerBond::new(&molecule, 0, 1, 2, BondStereo::AtropCcw).unwrap();
        let mut context = CipLabelerContext::new(0);
        assert_eq!(config.label(&rules, &mut context).unwrap(), Descriptor::P);
        assert_eq!(config.ranked_anchors(), &[0, 5]);

        let molecule = atropisomer_bond_molecule_one_side_second_carrier(BondStereo::AtropCw);
        let mut config = CipAtropisomerBond::new(&molecule, 0, 1, 2, BondStereo::AtropCw).unwrap();
        let mut context = CipLabelerContext::new(0);
        assert_eq!(config.label(&rules, &mut context).unwrap(), Descriptor::M);
        assert_eq!(config.ranked_anchors(), &[0, 5]);
    }

    #[test]
    fn ciplabeler_atropisomerbond_label_external_digraph_matches_rdkit_overload() {
        let molecule = atropisomer_bond_molecule(BondStereo::AtropCcw, false);
        let mut config = CipAtropisomerBond::new(&molecule, 0, 1, 2, BondStereo::AtropCcw).unwrap();
        let mut external = CipDigraph::new(&molecule, 1, true).unwrap();
        let root = external.get_original_root();
        let rules = cip_all_rules().unwrap();
        let mut context = CipLabelerContext::new(0);

        assert_eq!(
            config
                .label_with_external_digraph(root, &mut external, &rules, &mut context)
                .unwrap(),
            Descriptor::M
        );
        assert_eq!(config.ranked_anchors(), &[4, 5]);
    }

    #[test]
    fn ciplabeler_atropisomerbond_label_restores_the_original_root_before_labeling() {
        let molecule = atropisomer_bond_molecule(BondStereo::AtropCcw, false);
        let mut config = CipAtropisomerBond::new(&molecule, 0, 1, 2, BondStereo::AtropCcw).unwrap();
        let focus2_node = config.configuration.digraph.get_nodes(2).unwrap()[0];
        config
            .configuration
            .digraph
            .change_root(focus2_node)
            .unwrap();
        let rules = cip_all_rules().unwrap();
        let mut context = CipLabelerContext::new(0);

        assert_eq!(config.label(&rules, &mut context).unwrap(), Descriptor::M);
        assert_eq!(config.ranked_anchors(), &[4, 5]);
    }

    #[test]
    fn ciplabeler_atropisomerbond_label_reversed_root_preserves_source_anchor_order() {
        let molecule = atropisomer_bond_molecule(BondStereo::AtropCcw, false);
        let mut config = CipAtropisomerBond::new(&molecule, 0, 1, 2, BondStereo::AtropCcw).unwrap();
        let mut external = CipDigraph::new(&molecule, 2, true).unwrap();
        let root = external.get_original_root();
        let rules = cip_all_rules().unwrap();
        let mut context = CipLabelerContext::new(0);

        assert_eq!(
            config
                .label_with_external_digraph(root, &mut external, &rules, &mut context)
                .unwrap(),
            Descriptor::M
        );
        assert_eq!(config.ranked_anchors(), &[4, 5]);
    }

    #[test]
    fn ciplabeler_atropisomerbond_label_assigns_pseudoasymmetric_m_and_p_like_rdkit() {
        let rules = CipRules::new(vec![Box::new(PseudoAsymSortRule)]).unwrap();
        for (stereo, expected) in [
            (BondStereo::AtropCcw, Descriptor::m),
            (BondStereo::AtropCw, Descriptor::p),
        ] {
            let molecule = atropisomer_bond_molecule(stereo, false);
            let mut config = CipAtropisomerBond::new(&molecule, 0, 1, 2, stereo).unwrap();
            let mut context = CipLabelerContext::new(0);

            assert_eq!(config.label(&rules, &mut context).unwrap(), expected);
            assert_eq!(config.ranked_anchors(), &[4, 5]);
        }
    }

    #[test]
    fn ciplabeler_atropisomerbond_label_returns_unknown_for_nonunique_substituents() {
        let molecule = atropisomer_bond_molecule(BondStereo::AtropCcw, false);
        let mut config = CipAtropisomerBond::new(&molecule, 0, 1, 2, BondStereo::AtropCcw).unwrap();
        config.ranked_anchors = vec![u32::MAX];
        let rules = CipRules::new(vec![Box::new(AlwaysEqualSortRule)]).unwrap();
        let mut context = CipLabelerContext::new(0);

        assert_eq!(
            config.label(&rules, &mut context).unwrap(),
            Descriptor::Unknown
        );
        assert!(config.ranked_anchors().is_empty());
    }

    #[test]
    fn ciplabeler_atropisomerbond_label_returns_unknown_without_internal_edge() {
        let molecule = atropisomer_bond_molecule(BondStereo::AtropCcw, false);
        let mut config = CipAtropisomerBond::new(&molecule, 0, 1, 2, BondStereo::AtropCcw).unwrap();
        config.ranked_anchors = vec![u32::MAX];
        let mut external = CipDigraph::new(&molecule, 0, true).unwrap();
        let root = external.get_original_root();
        let rules = cip_all_rules().unwrap();
        let mut context = CipLabelerContext::new(0);

        assert_eq!(
            config
                .label_with_external_digraph(root, &mut external, &rules, &mut context)
                .unwrap(),
            Descriptor::Unknown
        );
        assert!(config.ranked_anchors().is_empty());
    }

    #[test]
    fn ciplabeler_configs_golden_configuration_parity_and_base_label_match_rdkit() {
        let reference = [0, 1, 2, 3];
        let rdkit_parity4_golden = [
            ([0, 1, 2, 3], 2),
            ([0, 1, 3, 2], 1),
            ([0, 2, 1, 3], 1),
            ([0, 2, 3, 1], 2),
            ([0, 3, 1, 2], 2),
            ([0, 3, 2, 1], 1),
            ([1, 0, 2, 3], 1),
            ([1, 0, 3, 2], 2),
            ([1, 2, 0, 3], 2),
            ([1, 2, 3, 0], 1),
            ([1, 3, 0, 2], 1),
            ([1, 3, 2, 0], 2),
            ([2, 0, 1, 3], 2),
            ([2, 0, 3, 1], 1),
            ([2, 1, 0, 3], 1),
            ([2, 1, 3, 0], 2),
            ([2, 3, 0, 1], 2),
            ([2, 3, 1, 0], 1),
            ([3, 0, 1, 2], 1),
            ([3, 0, 2, 1], 2),
            ([3, 1, 0, 2], 2),
            ([3, 1, 2, 0], 1),
            ([3, 2, 0, 1], 1),
            ([3, 2, 1, 0], 2),
        ];
        for (target, expected) in rdkit_parity4_golden {
            assert_eq!(
                CipConfiguration::parity4(&target, &reference).unwrap(),
                expected,
                "RDKit parity4 golden mismatch for target {target:?}"
            );
        }

        let molecule = Molecule::from_smiles("CCO").unwrap();
        let config = CipConfiguration::new(&molecule, 1).unwrap();
        let mut digraph = CipDigraph::new(&molecule, 1, false).unwrap();
        let root = digraph.get_original_root();
        let rules = cip_all_rules().unwrap();
        assert_eq!(
            config.label(root, &mut digraph, &rules),
            Descriptor::Unknown
        );
    }

    #[test]
    fn ciplabeler_configs_golden_tetrahedral_labels_reset_and_neighbor_order_match_rdkit() {
        let rules = cip_all_rules().unwrap();
        let pseudo_rules = CipRules::new(vec![Box::new(SinglePseudoAsymSortRule)]).unwrap();

        let tetrahedral_cases = [
            (ChiralTag::TetrahedralCcw, &rules, Descriptor::S, "S"),
            (ChiralTag::TetrahedralCw, &rules, Descriptor::R, "R"),
            (ChiralTag::TetrahedralCcw, &pseudo_rules, Descriptor::s, "s"),
            (ChiralTag::TetrahedralCw, &pseudo_rules, Descriptor::r, "r"),
        ];
        for (tag, rules, expected_descriptor, expected_code) in tetrahedral_cases {
            let molecule = tetrahedral_molecule(tag);
            let mut config = CipTetrahedral::new(&molecule, 0).unwrap();
            let mut context = CipLabelerContext::new(0);
            assert_eq!(
                config.label(rules, &mut context).unwrap(),
                expected_descriptor
            );
            assert_eq!(config.ranked_anchors(), &[4, 3, 2, 1]);

            config.set_primary_label(expected_descriptor).unwrap();
            assert_eq!(
                config.primary_label(),
                Some(&super::CipAtomPrimaryLabel {
                    atom_idx: 0,
                    cip_code: expected_code,
                    cip_neighbor_order: vec![4, 3, 2, 1],
                })
            );
            assert!(config.has_primary_label());

            config.reset_primary_label();
            assert_eq!(config.primary_label(), None);
            assert!(!config.has_primary_label());
        }
    }

    #[test]
    fn ciplabeler_configs_golden_sp2bond_labels_reset_and_neighbor_order_match_rdkit() {
        let rules = cip_all_rules().unwrap();
        let sp2_cases = [
            (BondStereo::Cis, Descriptor::Z, "Z"),
            (BondStereo::Trans, Descriptor::E, "E"),
        ];
        for (stereo, expected_descriptor, expected_code) in sp2_cases {
            let molecule = sp2_bond_molecule(stereo, (4, 5));
            let mut config = CipSp2Bond::new(&molecule, 0, 1, 2, stereo).unwrap();
            let mut context = CipLabelerContext::new(0);
            assert_eq!(
                config.label(&rules, &mut context).unwrap(),
                expected_descriptor
            );
            assert_eq!(config.ranked_anchors(), &[4, 5]);

            config.set_primary_label(expected_descriptor).unwrap();
            assert_eq!(
                config.primary_label(),
                Some(&super::CipBondPrimaryLabel {
                    bond_idx: 0,
                    stereo_atoms: [4, 5],
                    stereo,
                    cip_code: expected_code,
                    cip_neighbor_order: vec![4, 5],
                })
            );
            assert!(config.has_primary_label());

            config.reset_primary_label();
            assert_eq!(config.primary_label(), None);
            assert!(!config.has_primary_label());
        }
    }

    #[test]
    fn ciplabeler_configs_golden_atropisomerbond_labels_reset_and_neighbor_order_match_rdkit() {
        let rules = cip_all_rules().unwrap();
        let atrop_cases = [
            (BondStereo::AtropCcw, Descriptor::M, "M"),
            (BondStereo::AtropCw, Descriptor::P, "P"),
        ];
        for (stereo, expected_descriptor, expected_code) in atrop_cases {
            let molecule = atropisomer_bond_molecule(stereo, false);
            let mut config = CipAtropisomerBond::new(&molecule, 0, 1, 2, stereo).unwrap();
            let mut context = CipLabelerContext::new(0);
            assert_eq!(
                config.label(&rules, &mut context).unwrap(),
                expected_descriptor
            );
            assert_eq!(config.ranked_anchors(), &[4, 5]);

            config.set_primary_label(expected_descriptor).unwrap();
            assert_eq!(
                config.primary_label(),
                Some(&super::CipAtropisomerBondPrimaryLabel {
                    bond_idx: 0,
                    cip_code: expected_code,
                    cip_neighbor_order: vec![4, 5],
                })
            );
            assert!(config.has_primary_label());

            config.reset_primary_label();
            assert_eq!(config.primary_label(), None);
            assert!(!config.has_primary_label());
        }
    }

    #[test]
    fn ciplabeler_assign_writes_tetrahedral_labels_and_molecule_computed_prop_like_rdkit() {
        let molecule = tetrahedral_molecule(ChiralTag::TetrahedralCcw);
        let labeled = assign_cip_labels(&molecule, 0).unwrap();

        assert_eq!(molecule.prop("_CIPComputed"), None);
        assert_eq!(molecule.atoms()[0].prop("_CIPCode"), None);
        assert_eq!(labeled.prop("_CIPComputed"), Some("1"));
        assert_eq!(labeled.atoms()[0].prop("_CIPCode"), Some("S"));
        assert_eq!(
            labeled.atoms()[0].prop("_CIPNeighborOrder"),
            Some("[4,3,2,1]")
        );
    }

    #[test]
    fn ciplabeler_assign_writes_sp2bond_labels_stereo_atoms_and_neighbor_order_like_rdkit() {
        let molecule = sp2_bond_molecule(BondStereo::Cis, (4, 5));
        let labeled = assign_cip_labels(&molecule, 0).unwrap();

        assert_eq!(labeled.prop("_CIPComputed"), Some("1"));
        assert_eq!(labeled.bonds()[0].prop("_CIPCode"), Some("Z"));
        assert_eq!(labeled.bonds()[0].prop("_CIPNeighborOrder"), Some("[4,5]"));
        assert_eq!(labeled.bonds()[0].stereo(), BondStereo::Cis);
        assert_eq!(
            labeled.bonds()[0].stereo_atoms(),
            Some([crate::AtomId::new(4), crate::AtomId::new(5)])
        );
    }

    #[test]
    fn ciplabeler_assign_writes_atropisomer_bond_labels_like_rdkit() {
        let molecule = atropisomer_bond_molecule(BondStereo::AtropCcw, false);
        let labeled = assign_cip_labels(&molecule, 0).unwrap();

        assert_eq!(labeled.prop("_CIPComputed"), Some("1"));
        assert_eq!(labeled.bonds()[0].prop("_CIPCode"), Some("M"));
        assert_eq!(labeled.bonds()[0].prop("_CIPNeighborOrder"), Some("[4,5]"));
    }

    #[test]
    fn ciplabeler_assign_selected_masks_only_label_selected_atoms_and_bonds_like_rdkit() {
        let mut molecule = tetrahedral_molecule(ChiralTag::TetrahedralCcw);
        molecule.topology_block_mut().atoms[0].set_prop("_CIPCode", "old");
        molecule.topology_block_mut().atoms[0].set_prop("_CIPNeighborOrder", "[0]");
        let atom_mask = vec![false; molecule.num_atoms()];
        let bond_mask = vec![false; molecule.num_bonds()];

        let labeled = assign_cip_labels_for_masks(&molecule, &atom_mask, &bond_mask, 0).unwrap();

        assert_eq!(labeled.prop("_CIPComputed"), Some("1"));
        assert_eq!(labeled.atoms()[0].prop("_CIPCode"), Some("old"));
        assert_eq!(labeled.atoms()[0].prop("_CIPNeighborOrder"), Some("[0]"));

        let mut atom_mask = vec![false; molecule.num_atoms()];
        atom_mask[0] = true;
        let labeled = assign_cip_labels_for_masks(&molecule, &atom_mask, &bond_mask, 0).unwrap();
        assert_eq!(labeled.atoms()[0].prop("_CIPCode"), Some("S"));
        assert_eq!(
            labeled.atoms()[0].prop("_CIPNeighborOrder"),
            Some("[4,3,2,1]")
        );
    }

    #[test]
    fn ciplabeler_assign_selected_atom_preserves_unselected_primary_state_like_rdkit() {
        let mut molecule =
            Molecule::from_smiles_with_sanitize("F[C@](Cl)(Br)I.F[C@@](Cl)(Br)I", false).unwrap();
        molecule.topology_block_mut().atoms[1].set_prop("_CIPCode", "selected-old");
        molecule.topology_block_mut().atoms[1].set_prop("_CIPNeighborOrder", "[91]");
        molecule.topology_block_mut().atoms[6].set_prop("_CIPCode", "unselected-old");
        molecule.topology_block_mut().atoms[6].set_prop("_CIPNeighborOrder", "[92]");
        let mut atom_mask = vec![false; molecule.num_atoms()];
        atom_mask[1] = true;
        let bond_mask = vec![false; molecule.num_bonds()];

        let labeled = assign_cip_labels_for_masks(&molecule, &atom_mask, &bond_mask, 0).unwrap();

        assert_ne!(labeled.atoms()[1].prop("_CIPCode"), Some("selected-old"));
        assert_ne!(labeled.atoms()[1].prop("_CIPNeighborOrder"), Some("[91]"));
        assert_eq!(labeled.atoms()[6].prop("_CIPCode"), Some("unselected-old"));
        assert_eq!(labeled.atoms()[6].prop("_CIPNeighborOrder"), Some("[92]"));
    }

    #[test]
    fn ciplabeler_assign_unresolved_selected_center_clears_only_primary_code_like_rdkit() {
        let mut molecule = Molecule::from_smiles_with_sanitize("C[C@H](C)C", false).unwrap();
        molecule.topology_block_mut().atoms[1].set_prop("_CIPCode", "stale");
        molecule.topology_block_mut().atoms[1].set_prop("_CIPNeighborOrder", "[7,8,9]");
        let mut atom_mask = vec![false; molecule.num_atoms()];
        atom_mask[1] = true;
        let bond_mask = vec![false; molecule.num_bonds()];

        let labeled = assign_cip_labels_for_masks(&molecule, &atom_mask, &bond_mask, 0).unwrap();

        assert_eq!(labeled.atoms()[1].prop("_CIPCode"), None);
        assert_eq!(
            labeled.atoms()[1].prop("_CIPNeighborOrder"),
            Some("[7,8,9]")
        );
        assert_eq!(labeled.prop("_CIPComputed"), Some("1"));
    }

    #[test]
    fn ciplabeler_assign_empty_masks_replace_computed_mark_without_touching_center_state() {
        let mut molecule = tetrahedral_molecule(ChiralTag::TetrahedralCcw);
        molecule.properties_mut().set_prop("_CIPComputed", "stale");
        molecule.topology_block_mut().atoms[0].set_prop("_CIPCode", "old");
        molecule.topology_block_mut().atoms[0].set_prop("_CIPNeighborOrder", "[9]");

        let labeled = assign_cip_labels_for_masks(
            &molecule,
            &vec![false; molecule.num_atoms()],
            &vec![false; molecule.num_bonds()],
            0,
        )
        .unwrap();

        assert_eq!(labeled.prop("_CIPComputed"), Some("1"));
        assert_eq!(labeled.atoms()[0].prop("_CIPCode"), Some("old"));
        assert_eq!(labeled.atoms()[0].prop("_CIPNeighborOrder"), Some("[9]"));
    }

    #[test]
    fn ciplabeler_assign_normalizes_e_z_configurations_to_trans_cis_like_rdkit() {
        for (input, normalized) in [
            (BondStereo::E, BondStereo::Trans),
            (BondStereo::Z, BondStereo::Cis),
        ] {
            let molecule = sp2_bond_molecule(input, (4, 5));
            let labeled = assign_cip_labels(&molecule, 0).unwrap();

            assert_eq!(labeled.bonds()[0].stereo(), normalized);
            assert!(matches!(
                labeled.bonds()[0].prop("_CIPCode"),
                Some("E" | "Z")
            ));
        }
    }

    #[test]
    fn ciplabeler_assign_ignores_bond_stereo_outside_find_configs_dispatcher() {
        let mut builder = MoleculeBuilder::new();
        let begin = builder.add_atom(AtomSpec::new(Element::C));
        let end = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(
                BondSpec::new(begin, end, BondOrder::Double)
                    .with_stereo(BondStereo::Any)
                    .with_prop("_CIPCode", "preserved")
                    .with_prop("_CIPNeighborOrder", "[4,5]"),
            )
            .unwrap();
        let molecule = builder.build().unwrap();

        let labeled = assign_cip_labels(&molecule, 0).unwrap();

        assert_eq!(labeled.bonds()[0].prop("_CIPCode"), Some("preserved"));
        assert_eq!(labeled.bonds()[0].prop("_CIPNeighborOrder"), Some("[4,5]"));
        assert_eq!(labeled.prop("_CIPComputed"), Some("1"));
    }

    #[test]
    fn ciplabeler_assign_atom_and_bond_masks_share_one_ordered_assignment_call() {
        let mut builder = MoleculeBuilder::new();
        let center =
            builder.add_atom(AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCcw));
        let fluorine = builder.add_atom(AtomSpec::new(Element::F));
        let chlorine = builder.add_atom(AtomSpec::new(Element::CL));
        let bromine = builder.add_atom(AtomSpec::new(Element::BR));
        let iodine = builder.add_atom(AtomSpec::new(Element::I));
        for neighbor in [fluorine, chlorine, bromine, iodine] {
            builder
                .add_bond(BondSpec::new(center, neighbor, BondOrder::Single))
                .unwrap();
        }
        let alkene_begin = builder.add_atom(AtomSpec::new(Element::C));
        let alkene_end = builder.add_atom(AtomSpec::new(Element::C));
        let left_low = builder.add_atom(AtomSpec::new(Element::F));
        let left_high = builder.add_atom(AtomSpec::new(Element::BR));
        let right_low = builder.add_atom(AtomSpec::new(Element::CL));
        let right_high = builder.add_atom(AtomSpec::new(Element::I));
        let alkene = builder
            .add_bond(
                BondSpec::new(alkene_begin, alkene_end, BondOrder::Double)
                    .with_stereo(BondStereo::Cis)
                    .with_stereo_atoms(left_high, right_high),
            )
            .unwrap();
        for (focus, neighbor) in [
            (alkene_begin, left_low),
            (alkene_begin, left_high),
            (alkene_end, right_low),
            (alkene_end, right_high),
        ] {
            builder
                .add_bond(BondSpec::new(focus, neighbor, BondOrder::Single))
                .unwrap();
        }
        let molecule = builder.build().unwrap();
        let mut atom_mask = vec![false; molecule.num_atoms()];
        atom_mask[center.index()] = true;
        let mut bond_mask = vec![false; molecule.num_bonds()];
        bond_mask[alkene.index()] = true;

        let labeled = assign_cip_labels_for_masks(&molecule, &atom_mask, &bond_mask, 0).unwrap();

        assert_eq!(labeled.atoms()[center.index()].prop("_CIPCode"), Some("S"));
        assert!(matches!(
            labeled.bonds()[alkene.index()].prop("_CIPCode"),
            Some("E" | "Z")
        ));
    }

    #[test]
    fn ciplabeler_assign_second_pass_uses_the_requested_shared_recursion_budget() {
        let molecule =
            Molecule::from_smiles("C/C=C/[C@@H](/C=C/O)[C@H](C)[C@H](/C=C/C)/C=C/O").unwrap();

        assert!(matches!(
            assign_cip_labels(&molecule, 1),
            Err(CipLabelerError::MaxIterationsExceeded)
        ));
        let labeled = assign_cip_labels(&molecule, 0).unwrap();
        assert!(
            labeled
                .atoms()
                .iter()
                .any(|atom| matches!(atom.prop("_CIPCode"), Some("r" | "s")))
        );
    }

    #[test]
    fn ciplabeler_assign_max_recursive_iterations_does_not_disable_constitutional_fast_pass_like_rdkit()
     {
        let molecule = tetrahedral_molecule(ChiralTag::TetrahedralCcw);
        let labeled = assign_cip_labels(&molecule, 1).unwrap();

        assert_eq!(labeled.prop("_CIPComputed"), Some("1"));
        assert_eq!(labeled.atoms()[0].prop("_CIPCode"), Some("S"));
        assert_eq!(
            labeled.atoms()[0].prop("_CIPNeighborOrder"),
            Some("[4,3,2,1]")
        );
    }

    #[test]
    fn ciplabeler_exact_width_source_sentinels_are_unsigned_int_max() {
        assert_eq!(CipNode::NO_ATOM_INDEX, 4_294_967_295);
        assert_eq!(CipConfiguration::IMPLICIT_H, 4_294_967_295);
    }

    #[test]
    fn ciplabeler_exact_width_implicit_h_neighbor_order_serializes_as_u32_like_rdkit() {
        // VS185 from the pinned RDKit CIPLabeler catch tests exercises a real
        // _CIPNeighborOrder containing Atom::NOATOM.
        let mut molecule = Molecule::from_smiles(r"[2H]/C(=C(/[1H])\[H])/[H]").unwrap();
        molecule.topology_block_mut().bonds[1]
            .set_stereo_atoms(Some([crate::AtomId::new(0), crate::AtomId::new(3)]));

        let labeled = assign_cip_labels(&molecule, 100).unwrap();

        assert_eq!(
            labeled.bonds()[1].prop("_CIPNeighborOrder"),
            Some("[0,4294967295]")
        );
        assert!(
            !labeled.bonds()[1]
                .prop("_CIPNeighborOrder")
                .unwrap()
                .contains("18446744073709551615")
        );
    }

    #[test]
    fn ciplabeler_exact_width_real_atom_index_conversion_accepts_u32_max() {
        let node = CipNode {
            digraph: 0,
            atom_idx: Some(u32::MAX as usize),
            distance: 0,
            atomic_num_fraction: RationalI32::new(6, 1),
            atomic_mass: 12.0,
            aux: Descriptor::None,
            flags: 0,
            edges: Vec::new(),
            visit: Vec::new(),
        };

        assert_eq!(node.get_atom_idx().unwrap(), u32::MAX);
    }

    #[cfg(target_pointer_width = "64")]
    #[test]
    fn ciplabeler_exact_width_real_atom_index_conversion_rejects_above_u32_max() {
        let index = u32::MAX as usize + 1;
        let node = CipNode {
            digraph: 0,
            atom_idx: Some(index),
            distance: 0,
            atomic_num_fraction: RationalI32::new(6, 1),
            atomic_mass: 12.0,
            aux: Descriptor::None,
            flags: 0,
            edges: Vec::new(),
            visit: Vec::new(),
        };

        assert_eq!(
            node.get_atom_idx(),
            Err(CipLabelerError::SourceIndexWidthExceeded {
                kind: "atom",
                index,
            })
        );
    }

    #[test]
    fn ciplabeler_exact_width_selection_indexes_equal_to_counts_are_rejected() {
        let molecule = Molecule::from_smiles("CC").unwrap();

        assert_eq!(
            assign_cip_labels_for_indices(&molecule, &[molecule.num_atoms() as u32], &[], 0,),
            Err(CipLabelerError::AtomIndexOutOfRange {
                index: molecule.num_atoms(),
                atom_count: molecule.num_atoms(),
            })
        );
        assert_eq!(
            assign_cip_labels_for_indices(&molecule, &[], &[molecule.num_bonds() as u32], 0,),
            Err(CipLabelerError::BondIndexOutOfRange {
                index: molecule.num_bonds(),
                bond_count: molecule.num_bonds(),
            })
        );
    }

    #[test]
    fn ciplabeler_exact_width_duplicate_selected_indexes_are_idempotent() {
        let molecule = tetrahedral_molecule(ChiralTag::TetrahedralCcw);

        let once = assign_cip_labels_for_indices(&molecule, &[0], &[], 0).unwrap();
        let repeated = assign_cip_labels_for_indices(&molecule, &[0, 0, 0], &[], 0).unwrap();

        assert_eq!(repeated, once);
    }

    #[test]
    fn ciplabeler_exact_width_short_and_long_selection_masks_are_rejected() {
        let molecule = Molecule::from_smiles("CC").unwrap();
        let atom_count = molecule.num_atoms();
        let bond_count = molecule.num_bonds();

        for actual in [atom_count - 1, atom_count + 1] {
            assert_eq!(
                assign_cip_labels_for_masks(
                    &molecule,
                    &vec![false; actual],
                    &vec![false; bond_count],
                    0,
                ),
                Err(CipLabelerError::AtomMaskLengthMismatch {
                    actual,
                    expected: atom_count,
                })
            );
        }
        for actual in [bond_count - 1, bond_count + 1] {
            assert_eq!(
                assign_cip_labels_for_masks(
                    &molecule,
                    &vec![false; atom_count],
                    &vec![false; actual],
                    0,
                ),
                Err(CipLabelerError::BondMaskLengthMismatch {
                    actual,
                    expected: bond_count,
                })
            );
        }
    }

    #[test]
    fn ciplabeler_exact_width_pairlist_uses_highest_and_lowest_relevant_u64_bits() {
        let mut pairs = CipPairList::with_ref(Descriptor::R);
        pairs.add(Descriptor::R);
        for _ in 0..61 {
            pairs.add(Descriptor::S);
        }
        pairs.add(Descriptor::R);

        assert_eq!(pairs.get_pairing(), (1_u64 << 62) | 1_u64);
    }

    #[test]
    fn ciplabeler_exact_width_recursion_counter_reproduces_unsigned_predecrement() {
        let mut zero = CipLabelerContext::with_remaining_call_count(0);
        assert!(zero.decrement_remaining_call_count_and_check());
        assert_eq!(zero.remaining_call_count, u32::MAX);

        let mut one = CipLabelerContext::with_remaining_call_count(1);
        assert!(!one.decrement_remaining_call_count_and_check());
        assert_eq!(one.remaining_call_count, 0);

        let mut two = CipLabelerContext::with_remaining_call_count(2);
        assert!(two.decrement_remaining_call_count_and_check());
        assert_eq!(two.remaining_call_count, 1);
        assert!(!two.decrement_remaining_call_count_and_check());
        assert_eq!(two.remaining_call_count, 0);
    }

    #[test]
    fn ciplabeler_cancellation_rejects_without_mutating_source_state() {
        let mut molecule = tetrahedral_molecule(ChiralTag::TetrahedralCcw);
        molecule.properties_mut().set_prop("_CIPComputed", "stale");
        molecule.topology_block_mut().atoms[0].set_prop("_CIPCode", "old");
        molecule.topology_block_mut().atoms[0].set_prop("_CIPNeighborOrder", "[9]");
        let before = molecule.clone();
        let atom_mask = vec![true; molecule.num_atoms()];
        let bond_mask = vec![true; molecule.num_bonds()];

        assert_eq!(
            assign_cip_labels_for_masks_with_cancellation(
                &molecule,
                &atom_mask,
                &bond_mask,
                0,
                CipCancellationMode::RdkitControlC,
            ),
            Err(CipLabelerError::CancellationUnsupported)
        );
        assert_eq!(molecule, before);
    }
}
