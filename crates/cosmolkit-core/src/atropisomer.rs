//! Detached, source-backed atropisomer perception and wedge assignment.

use crate::stereo_graph::StereoGraphAccess;
use crate::stereo_graph::{StereoAtomAccess, StereoGraphMut};
use std::collections::{BTreeMap, BTreeSet};

use cosmolkit_model::{
    AtomId, Bond, BondId, BondValueError, Conformer2D, Conformer3D, CoordinateValidationError,
    SourceAtomValenceFacts, StereoGroup, TopologyBlock, TopologyValidationError,
};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo, Hybridization};

use crate::{HybridizationAssignment, RingInfo, WedgeAssignments, WedgeInfo};

const REALLY_SMALL_BOND_LEN: f64 = 0.000_000_1;

#[derive(Debug, Clone, Copy)]
pub enum AtropisomerConformer<'a> {
    TwoD(&'a Conformer2D),
    ThreeD(&'a Conformer3D),
}

fn conformer_is_3d(conformer: AtropisomerConformer<'_>) -> bool {
    // BEGIN RDKIT CPP FUNCTION getBondFrameOfReference
    // RDKit✔️✔️:   if (!conf->is3D()) {
    // RDKit✔️✔️:     yAxis = RDGeom::Point3D(-xAxis.y, xAxis.x, 0);
    // RDKit✔️✔️:     yAxis.normalize();
    // RDKit✔️✔️:     zAxis = RDGeom::Point3D(0.0, 0.0, 1.0);
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION
    // Behavior review: a Conformer3D is a lossless XYZ storage carrier, not
    // proof that the source conformer flag is 3D. The independent is_3d bit
    // selects the same source branch while all XYZ bits remain available.
    // Complexity review: one enum match and one flag read are constant-time
    // and introduce no allocation or coordinate projection.
    match conformer {
        AtropisomerConformer::TwoD(_) => false,
        AtropisomerConformer::ThreeD(value) => value.is_3d(),
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AtropisomerRejectionKind {
    MissingCarrier,
    UnknownCarrierDirection,
    InconsistentDirections,
    ZeroLengthAxis,
    CollinearCarrier,
    SameSideCarriers,
    CoplanarCarriers,
    NoUsableWedgeBond,
    DirectionConflict,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AtropisomerDiagnostic {
    pub bond: BondId,
    pub kind: AtropisomerRejectionKind,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AtropisomerBondUpdate {
    pub bond: BondId,
    pub stereo: BondStereo,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct AtropisomerAssignment {
    pub bond_updates: Vec<AtropisomerBondUpdate>,
    pub diagnostics: Vec<AtropisomerDiagnostic>,
    /// Source cache writes, in execution order. No graph or commit authority.
    pub atom_valence_updates: Vec<(AtomId, SourceAtomValenceFacts)>,
    /// Present only when the native whole-molecule prelude executes.
    pub conjugated_bonds: Option<Vec<bool>>,
    pub hybridization: Option<HybridizationAssignment>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AtropisomerWedgeUpdate {
    pub bond: BondId,
    pub begin: AtomId,
    pub end: AtomId,
    pub direction: BondDirection,
    pub atropisomer_bond: BondId,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct AtropisomerWedgeAssignment {
    /// Carrier IDs where the source actually inserts or replaces WedgeInfo.
    pub source_map_writes: Vec<BondId>,
    pub bond_updates: Vec<AtropisomerWedgeUpdate>,
    pub diagnostics: Vec<AtropisomerDiagnostic>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtropisomerCarrierEnd {
    focus: AtomId,
    carrier_bonds: Vec<BondId>,
}

impl AtropisomerCarrierEnd {
    #[must_use]
    pub const fn focus(&self) -> AtomId {
        self.focus
    }

    #[must_use]
    pub fn carrier_bonds(&self) -> &[BondId] {
        &self.carrier_bonds
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct StereoGroupAssignment {
    pub groups: Vec<StereoGroup>,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum AtropisomerError {
    #[error("{0}")]
    StereoGroup(#[from] cosmolkit_model::StereoGroupError),

    #[error("source atropisomer precondition failed: {message}")]
    SourcePrecondition { message: &'static str },
    #[error("candidate traversal bond {bond} is outside {bond_count} source bond rows")]
    CandidateTraversalOutOfRange { bond: BondId, bond_count: usize },
    #[error("candidate traversal repeats source bond {bond}")]
    CandidateTraversalDuplicate { bond: BondId },
    #[error("candidate traversal contains non-candidate source bond {bond}")]
    CandidateTraversalUnexpected { bond: BondId },
    #[error("candidate traversal omits source bond {bond}")]
    CandidateTraversalMissing { bond: BondId },
    #[error("source atropisomer cache/getter failed: {0}")]
    Valence(#[from] crate::ValenceError),
    #[error("source atropisomer conjugation failed: {0}")]
    Conjugation(#[from] crate::ConjugationError),
    #[error("source atropisomer hybridization failed: {0}")]
    Hybridization(#[from] crate::HybridizationError),
    #[error("source Point3D normalization failed: {0}")]
    Normalization(crate::StereoError),
    #[error("source carrier count precondition failed at atom {atom}: {count}")]
    CarrierCount { atom: AtomId, count: usize },
    #[error("source atropisomer ring acquisition failed: {0}")]
    RingFinding(#[from] crate::RingFindingError),
    #[error("invalid topology: {source}")]
    InvalidTopology { source: TopologyValidationError },
    #[error("invalid coordinates: {source}")]
    InvalidCoordinates { source: CoordinateValidationError },
    #[error("3D conformer {conformer} is marked as non-3D")]
    ConformerNotThreeDimensional { conformer: usize },
    #[error("ring information must be SSSR or better")]
    RingInfoNotSssr,
    #[error("ring atom {atom} is out of range for {atom_count} atoms")]
    RingAtomOutOfRange { atom: AtomId, atom_count: usize },
    #[error("ring bond {bond} is out of range for {bond_count} bonds")]
    RingBondOutOfRange { bond: BondId, bond_count: usize },
    #[error("ring atom membership has {actual} rows; expected {expected}")]
    RingAtomRowCount { actual: usize, expected: usize },
    #[error("ring bond membership has {actual} rows; expected {expected}")]
    RingBondRowCount { actual: usize, expected: usize },
    #[error("hybridization assignment has {actual} rows; expected {expected}")]
    HybridizationAssignmentLength { actual: usize, expected: usize },
    #[error("failed to clear atropisomer stereo from bond {bond}: {source}")]
    BondUpdate {
        bond: BondId,
        source: BondValueError,
    },
    #[error("stereo group atom {atom} is out of range for {atom_count} atoms")]
    StereoGroupAtomOutOfRange { atom: AtomId, atom_count: usize },
    #[error("stereo group bond {bond} is out of range for {bond_count} bonds")]
    StereoGroupBondOutOfRange { bond: BondId, bond_count: usize },
    #[error("atropisomer assignment bond {bond} is out of range for {bond_count} bonds")]
    AssignmentBondOutOfRange { bond: BondId, bond_count: usize },
    #[error("atropisomer axial bond {bond} is out of range for {bond_count} bonds")]
    AxialBondOutOfRange { bond: BondId, bond_count: usize },
    #[error("bond {bond} has invalid atropisomer stereo {stereo:?}")]
    InvalidAtropisomerStereo { bond: BondId, stereo: BondStereo },
}

#[derive(Debug, Clone)]
struct AtropEnd {
    atom: AtomId,
    bonds: Vec<BondId>,
}

#[derive(Debug, Clone, Copy)]
struct Frame {
    y: [f64; 3],
    z: [f64; 3],
}

fn can_have_direction(bond: &Bond) -> bool {
    // Complete pinned source: Bond.h::canHaveDirection.
    // RDKit✔️✔️: inline bool canHaveDirection(const Bond &bond) {
    // RDKit✔️✔️:   auto bondType = bond.getBondType();
    // RDKit✔️✔️:   return (bondType == Bond::SINGLE || bondType == Bond::AROMATIC);
    // RDKit✔️✔️: }
    matches!(bond.order(), BondOrder::Single | BondOrder::Aromatic)
}

#[derive(Debug, Clone, Copy)]
struct SourceCarrierError {
    bond: BondId,
    atom: AtomId,
}
impl From<SourceCarrierError> for AtropisomerError {
    fn from(_error: SourceCarrierError) -> Self {
        Self::SourcePrecondition {
            message: "bad index",
        }
    }
}
impl From<SourceCarrierError> for PerceptionError {
    fn from(_error: SourceCarrierError) -> Self {
        Self::SourcePrecondition {
            message: "bad index",
        }
    }
}

fn other_atom(bond: &Bond, atom: AtomId) -> Result<AtomId, SourceCarrierError> {
    // RDKit✔️✔️: Atom *Bond::getOtherAtom(Atom const *what) const {
    // RDKit✔️✔️:   PRECONDITION(dp_mol != nullptr, "no owning molecule for bond");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return dp_mol->getAtomWithIdx(getOtherAtomIdx(what->getIdx()));
    // RDKit✔️✔️: }
    // RDKit✔️✔️: unsigned int Bond::getOtherAtomIdx(const unsigned int thisIdx) const {
    // RDKit✔️✔️:   if (d_beginAtomIdx == thisIdx) {
    // RDKit✔️✔️:     return d_endAtomIdx;
    // RDKit✔️✔️:   } else if (d_endAtomIdx == thisIdx) {
    // RDKit✔️✔️:     return d_beginAtomIdx;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // This "precondition" would check exactly the same that is checked
    // RDKit✔️✔️:   // above, but no need to be redundant, so just throw.
    // RDKit✔️✔️:   POSTCONDITION(false, "bad index");
    // RDKit✔️✔️: }
    // Behavior/cost review: both native endpoint tests precede the actual bad
    // index postcondition; no non-incident-bond endpoint guess. Three scalar
    // branches, no allocation or additional graph scan.
    if bond.begin() == atom {
        Ok(bond.end())
    } else if bond.end() == atom {
        Ok(bond.begin())
    } else {
        Err(SourceCarrierError {
            bond: bond.id(),
            atom,
        })
    }
}

// Native getTotalDegree uses getTotalNumHs, including the exact stored count;
// a legacy boolean is not the source implicit-valence cache.
fn source_total_degree<G: StereoGraphAccess>(
    topology: &G,
    valence: &crate::ValenceAssignment,
    atom: AtomId,
) -> Result<u32, crate::ValenceError> {
    // BEGIN RDKIT CPP FUNCTION Atom::getTotalDegree complete source
    // RDKit✔️✔️: unsigned int Atom::getTotalDegree() const {
    // RDKit✔️✔️:   unsigned int res = this->getTotalNumHs(false) + this->getDegree();
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Atom::getTotalDegree complete source
    // RDKit✔️✔️: unsigned int Atom::getDegree() const {
    // RDKit✔️✔️:   return dp_mol ? getOwningMol().getAtomDegree(this) : 0;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: unsigned int ROMol::getAtomDegree(const Atom *at) const {
    // RDKit✔️✔️:   PRECONDITION(at, "no atom");
    // RDKit✔️✔️:   PRECONDITION(&at->getOwningMol() == this,
    // RDKit✔️✔️:                "atom not associated with this molecule");
    // RDKit✔️✔️:   return rdcast<unsigned int>(boost::out_degree(at->getIdx(), d_graph));
    // RDKit✔️✔️: };
    // RDKit✔️✔️: #ifdef RDDEBUG
    // RDKit✔️✔️: // use rdcast to convert between types
    // RDKit✔️✔️: //  when RDDEBUG is defined, this checks for
    // RDKit✔️✔️: //  validity (overflow, etc)
    // RDKit✔️✔️: //  when RDDEBUG is off, the cast is a no-cost
    // RDKit✔️✔️: //   static_cast
    // RDKit✔️✔️: #define rdcast boost::numeric_cast
    // RDKit✔️✔️: #else
    // RDKit✔️✔️: #define rdcast static_cast
    // RDKit✔️✔️: #endif
    // Normal pinned release build uses the source static_cast specialization.
    // Degree on an actual topology-owned atom is the undirected adjacency
    // row length (constant time), independent of bond order/implicit H flag.
    // Native uint32 addition wraps; totalNumHs(false) never counts H neighbors
    // separately, since they already contribute to this graph degree.
    Ok(
        crate::hcount::total_hydrogen_count_from_validated(topology, valence, atom, false)?
            .wrapping_add(topology.adjacency().neighbors_of(atom.index()).len() as u32),
    )
}

#[derive(Debug)]
enum PerceptionError {
    SourcePrecondition { message: &'static str },
    Rejected(AtropisomerRejectionKind),
    Normalization(crate::StereoError),
    CarrierCount { atom: AtomId, count: usize },
}
impl From<AtropisomerRejectionKind> for PerceptionError {
    fn from(value: AtropisomerRejectionKind) -> Self {
        Self::Rejected(value)
    }
}
impl From<crate::StereoError> for PerceptionError {
    fn from(value: crate::StereoError) -> Self {
        Self::Normalization(value)
    }
}

// Scalar views of actual current source Bond getters. The projected carrier
// retains graph directions/orientations separately; this view borrows those
// real endpoint values without making a synthetic Bond or cloning properties.
trait AtropBondEndpoints {
    fn id(&self) -> BondId;
    fn begin(&self) -> AtomId;
    fn end(&self) -> AtomId;
}
impl AtropBondEndpoints for Bond {
    fn id(&self) -> BondId {
        Bond::id(self)
    }
    fn begin(&self) -> AtomId {
        Bond::begin(self)
    }
    fn end(&self) -> AtomId {
        Bond::end(self)
    }
}
#[derive(Clone, Copy)]
struct CurrentAtropBondEndpoints {
    id: BondId,
    begin: AtomId,
    end: AtomId,
}
impl AtropBondEndpoints for CurrentAtropBondEndpoints {
    fn id(&self) -> BondId {
        self.id
    }
    fn begin(&self) -> AtomId {
        self.begin
    }
    fn end(&self) -> AtomId {
        self.end
    }
}

fn atropisomer_ends<G: StereoGraphAccess, B: AtropBondEndpoints>(
    topology: &G,
    bond: &B,
    ends: &mut [AtropEnd; 2],
) -> Result<bool, SourceCarrierError> {
    // Complete pinned source: getAtropisomerAtomsAndBonds.
    // RDKit✔️✔️: bool getAtropisomerAtomsAndBonds(const Bond *bond,
    // RDKit✔️✔️:                                  AtropAtomAndBondVec atomsAndBondVects[2],
    // RDKit✔️✔️:                                  const ROMol &mol) {
    // RDKit✔️✔️:   PRECONDITION(bond, "no bond");
    // RDKit✔️✔️:   atomsAndBondVects[0].first = bond->getBeginAtom();
    // RDKit✔️✔️:   atomsAndBondVects[1].first = bond->getEndAtom();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // get the one or two bonds on each end
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (int bondAtomIndex = 0; bondAtomIndex < 2; ++bondAtomIndex) {
    // RDKit✔️✔️:     for (const auto nbrBond :
    // RDKit✔️✔️:          mol.atomBonds(atomsAndBondVects[bondAtomIndex].first)) {
    // RDKit✔️✔️:       if (nbrBond == bond) {
    // RDKit✔️✔️:         continue;  // a bond is NOT its own neighbor
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       atomsAndBondVects[bondAtomIndex].second.push_back(nbrBond);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (atomsAndBondVects[bondAtomIndex].second.size() == 0) {
    // RDKit✔️✔️:       return false;  // no neighbor bonds found
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // make sure the bond with this lowest atom is is first
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (atomsAndBondVects[bondAtomIndex].second.size() == 2 &&
    // RDKit✔️✔️:         atomsAndBondVects[bondAtomIndex]
    // RDKit✔️✔️:                 .second[1]
    // RDKit✔️✔️:                 ->getOtherAtom(atomsAndBondVects[bondAtomIndex].first)
    // RDKit✔️✔️:                 ->getIdx() <
    // RDKit✔️✔️:             atomsAndBondVects[bondAtomIndex]
    // RDKit✔️✔️:                 .second[0]
    // RDKit✔️✔️:                 ->getOtherAtom(atomsAndBondVects[bondAtomIndex].first)
    // RDKit✔️✔️:                 ->getIdx()) {
    // RDKit✔️✔️:       std::swap(atomsAndBondVects[bondAtomIndex].second[0],
    // RDKit✔️✔️:                 atomsAndBondVects[bondAtomIndex].second[1]);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // Behavior review: borrow actual caller output; set both focus atoms before
    // scanning, append without clearing, preserve the first-end prefix on false,
    // and swap only a final vector of exactly two carriers. Valid source bond
    // references are represented by stable IDs in this validated topology.
    // Complexity review: two borrowed incident passes, amortized vector append,
    // constant-size conditional swap; no graph clone, sort or replacement buffer.
    ends[0].atom = bond.begin();
    ends[1].atom = bond.end();
    for end in ends {
        for neighbor in topology.adjacency().neighbors_of(end.atom.index()).iter() {
            if neighbor.bond == bond.id() {
                continue;
            }
            end.bonds.push(neighbor.bond);
        }
        if end.bonds.is_empty() {
            return Ok(false);
        }
        if end.bonds.len() == 2
            && other_atom(&topology.bonds()[end.bonds[1].index()], end.atom)?.index()
                < other_atom(&topology.bonds()[end.bonds[0].index()], end.atom)?.index()
        {
            end.bonds.swap(0, 1);
        }
    }
    Ok(true)
}

fn atropisomer_ends_fresh<G: StereoGraphAccess, B: AtropBondEndpoints>(
    topology: &G,
    bond: &B,
) -> Result<Option<[AtropEnd; 2]>, SourceCarrierError> {
    let mut ends = [
        AtropEnd {
            atom: bond.begin(),
            bonds: Vec::new(),
        },
        AtropEnd {
            atom: bond.end(),
            bonds: Vec::new(),
        },
    ];
    Ok(atropisomer_ends(topology, bond, &mut ends)?.then_some(ends))
}

pub fn atropisomer_carriers(
    topology: &TopologyBlock,
    axial_bond: BondId,
) -> Result<Option<[AtropisomerCarrierEnd; 2]>, AtropisomerError> {
    topology
        .validate()
        .map_err(|source| AtropisomerError::InvalidTopology { source })?;
    let bond =
        topology
            .bonds
            .get(axial_bond.index())
            .ok_or(AtropisomerError::AxialBondOutOfRange {
                bond: axial_bond,
                bond_count: topology.bonds.len(),
            })?;
    Ok(atropisomer_ends_fresh(topology, bond)
        .map_err(AtropisomerError::from)?
        .map(|ends| {
            ends.map(|end| AtropisomerCarrierEnd {
                focus: end.atom,
                carrier_bonds: end.bonds,
            })
        }))
}

/// Returns carrier ends for requested axial bonds after validating topology once.
pub fn atropisomer_carriers_for_bonds(
    topology: &TopologyBlock,
    axial_bonds: &[BondId],
) -> Result<Vec<(BondId, Option<[AtropisomerCarrierEnd; 2]>)>, AtropisomerError> {
    // Validate once, then reuse the source-shaped `atropisomer_ends` owner
    // helper for each requested bond instead of repeating a whole-topology
    // scan for every atropisomer neighbor in the CX bond-config writer.
    topology
        .validate()
        .map_err(|source| AtropisomerError::InvalidTopology { source })?;
    axial_bonds
        .iter()
        .copied()
        .map(|axial_bond| {
            let bond = topology.bonds.get(axial_bond.index()).ok_or(
                AtropisomerError::AxialBondOutOfRange {
                    bond: axial_bond,
                    bond_count: topology.bonds.len(),
                },
            )?;
            let carriers = atropisomer_ends_fresh(topology, bond)
                .map_err(AtropisomerError::from)?
                .map(|ends| {
                    ends.map(|end| AtropisomerCarrierEnd {
                        focus: end.atom,
                        carrier_bonds: end.bonds,
                    })
                });
            Ok((axial_bond, carriers))
        })
        .collect()
}

/// Source carrier ends over the actual validated query graph, without Atom coercion.
#[doc(hidden)]
pub fn query_atropisomer_carriers_source(
    query: &cosmolkit_model::QueryGraph,
    axial_bond: BondId,
) -> Result<Option<[AtropisomerCarrierEnd; 2]>, AtropisomerError> {
    let bond = query
        .bonds()
        .get(axial_bond.index())
        .ok_or(AtropisomerError::AxialBondOutOfRange {
            bond: axial_bond,
            bond_count: query.num_bonds(),
        })?
        .bond();
    // QueryGraph construction validates graph references. Reuse the source
    // incident-row kernel directly; no eager all-graph pass or carrier copy.
    Ok(atropisomer_ends_fresh(query, bond)
        .map_err(AtropisomerError::from)?
        .map(|ends| {
            ends.map(|end| AtropisomerCarrierEnd {
                focus: end.atom,
                carrier_bonds: end.bonds,
            })
        }))
}

fn point(conformer: AtropisomerConformer<'_>, atom: AtomId) -> [f64; 3] {
    match conformer {
        AtropisomerConformer::TwoD(value) => {
            let p = value.coordinates()[atom.index()];
            [p[0], p[1], 0.0]
        }
        AtropisomerConformer::ThreeD(value) => value.coordinates()[atom.index()],
    }
}

fn sub(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    // RDKit✔️✔️: Point3D operator-(const Point3D &p1, const Point3D &p2) {
    // RDKit✔️✔️:   Point3D res;
    // RDKit✔️✔️:   res.x = p1.x - p2.x;
    // RDKit✔️✔️:   res.y = p1.y - p2.y;
    // RDKit✔️✔️:   res.z = p1.z - p2.z;
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: constexpr Point3D &operator-=(const Point3D &other) {
    // RDKit✔️✔️:     x -= other.x;
    // RDKit✔️✔️:     y -= other.y;
    // RDKit✔️✔️:     z -= other.z;
    // RDKit✔️✔️:     return *this;
    // RDKit✔️✔️:   }
    // Behavior/cost review: native component arithmetic and evaluation order,
    // constant scalar stack operations without allocations or buffering.
    [left[0] - right[0], left[1] - right[1], left[2] - right[2]]
}

fn neg(value: [f64; 3]) -> [f64; 3] {
    // RDKit✔️✔️: constexpr Point3D operator-() const {
    // RDKit✔️✔️:     Point3D res(x, y, z);
    // RDKit✔️✔️:     res.x *= -1.0;
    // RDKit✔️✔️:     res.y *= -1.0;
    // RDKit✔️✔️:     res.z *= -1.0;
    // RDKit✔️✔️:     return res;
    // RDKit✔️✔️:   }
    // Behavior/cost review: native component arithmetic and evaluation order,
    // constant scalar stack operations without allocations or buffering.
    [value[0] * -1.0, value[1] * -1.0, value[2] * -1.0]
}

fn dot(left: [f64; 3], right: [f64; 3]) -> f64 {
    // RDKit✔️✔️: constexpr double dotProduct(const Point3D &other) const {
    // RDKit✔️✔️:     double res = x * (other.x) + y * (other.y) + z * (other.z);
    // RDKit✔️✔️:     return res;
    // RDKit✔️✔️:   }
    // Behavior/cost review: native component arithmetic and evaluation order,
    // constant scalar stack operations without allocations or buffering.
    left[0] * right[0] + left[1] * right[1] + left[2] * right[2]
}

fn cross(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    // RDKit✔️✔️: constexpr Point3D crossProduct(const Point3D &other) const {
    // RDKit✔️✔️:     Point3D res;
    // RDKit✔️✔️:     res.x = y * (other.z) - z * (other.y);
    // RDKit✔️✔️:     res.y = -x * (other.z) + z * (other.x);
    // RDKit✔️✔️:     res.z = x * (other.y) - y * (other.x);
    // RDKit✔️✔️:     return res;
    // RDKit✔️✔️:   }
    // Behavior/cost review: native component arithmetic and evaluation order,
    // constant scalar stack operations without allocations or buffering.
    [
        left[1] * right[2] - left[2] * right[1],
        -left[0] * right[2] + left[2] * right[0],
        left[0] * right[1] - left[1] * right[0],
    ]
}

fn length(value: [f64; 3]) -> f64 {
    // RDKit✔️✔️: double length() const override {
    // RDKit✔️✔️:     double res = x * x + y * y + z * z;
    // RDKit✔️✔️:     return sqrt(res);
    // RDKit✔️✔️:   }
    // Behavior/cost review: native component arithmetic and evaluation order,
    // constant scalar stack operations without allocations or buffering.
    dot(value, value).sqrt()
}

fn normalized(
    value: [f64; 3],
    center: AtomId,
    neighbor: AtomId,
) -> Result<[f64; 3], crate::StereoError> {
    crate::structure_tags::normalize_vector_components(value, center, neighbor)
}

fn frame_of_reference_source<B: AtropBondEndpoints>(
    bond: &B,
    conformer: AtropisomerConformer<'_>,
    x_axis: &mut [f64; 3],
    y_axis: &mut [f64; 3],
    z_axis: &mut [f64; 3],
) -> Result<bool, crate::StereoError> {
    // Complete pinned source: getBondFrameOfReference.
    // RDKit✔️✔️: bool getBondFrameOfReference(const Bond *bond, const Conformer *conf,
    // RDKit✔️✔️:                              RDGeom::Point3D &xAxis, RDGeom::Point3D &yAxis,
    // RDKit✔️✔️:                              RDGeom::Point3D &zAxis) {
    // RDKit✔️✔️:   // create a frame of reference that has its X-axis along the atrop bond
    // RDKit✔️✔️:   // for 2D confs, the yAxis is in the 2D plane and the zAxis is perpendicular
    // RDKit✔️✔️:   // to that plane) for 3D confs  the yAxis and the zAxis are arbitrary.
    // RDKit✔️✔️:
    // RDKit✔️✔️:   PRECONDITION(bond, "bad bond");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   xAxis = conf->getAtomPos(bond->getEndAtom()->getIdx()) -
    // RDKit✔️✔️:           conf->getAtomPos(bond->getBeginAtom()->getIdx());
    // RDKit✔️✔️:   if (xAxis.length() < REALLY_SMALL_BOND_LEN) {
    // RDKit✔️✔️:     return false;  // bond len is xero
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   xAxis.normalize();
    // RDKit✔️✔️:   if (!conf->is3D()) {
    // RDKit✔️✔️:     yAxis = RDGeom::Point3D(-xAxis.y, xAxis.x, 0);
    // RDKit✔️✔️:     yAxis.normalize();
    // RDKit✔️✔️:     zAxis = RDGeom::Point3D(0.0, 0.0, 1.0);
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // here for 3D conf
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (fabs(xAxis.x) > REALLY_SMALL_BOND_LEN ||
    // RDKit✔️✔️:       fabs(xAxis.y) > REALLY_SMALL_BOND_LEN) {
    // RDKit✔️✔️:     zAxis = RDGeom::Point3D(
    // RDKit✔️✔️:         0, 0, 1);  // since X or Y value of the new x xaxis is NOT 0, this
    // RDKit✔️✔️:                    // new temp z axis cannnot be colinear with it
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     zAxis = RDGeom::Point3D(1, 0, 0);  // since the new x axis is exactly along
    // RDKit✔️✔️:                                        // the (old) z axis, this new temp z axis
    // RDKit✔️✔️:                                        // is NOT colinear with it
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   yAxis = zAxis.crossProduct(xAxis);
    // RDKit✔️✔️:   zAxis = xAxis.crossProduct(yAxis);
    // RDKit✔️✔️:   yAxis.normalize();
    // RDKit✔️✔️:   zAxis.normalize();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // Behavior review: overwrite the actual source outputs in statement order.
    // A short axis preserves the new raw X while caller Y/Z are untouched.
    // The 3D path writes both cross products before either normalization.
    // Complexity review: constant stack vectors and scalar operations only;
    // no coordinate projection, graph clone, heap buffers or repeated scan.
    *x_axis = sub(point(conformer, bond.end()), point(conformer, bond.begin()));
    if length(*x_axis) < REALLY_SMALL_BOND_LEN {
        return Ok(false);
    }
    *x_axis = normalized(*x_axis, bond.begin(), bond.end())?;
    if !conformer_is_3d(conformer) {
        *y_axis = [-x_axis[1], x_axis[0], 0.0];
        *y_axis = normalized(*y_axis, bond.begin(), bond.end())?;
        *z_axis = [0.0, 0.0, 1.0];
        return Ok(true);
    }
    *z_axis = if x_axis[0].abs() > REALLY_SMALL_BOND_LEN || x_axis[1].abs() > REALLY_SMALL_BOND_LEN
    {
        [0.0, 0.0, 1.0]
    } else {
        [1.0, 0.0, 0.0]
    };
    *y_axis = cross(*z_axis, *x_axis);
    *z_axis = cross(*x_axis, *y_axis);
    *y_axis = normalized(*y_axis, bond.begin(), bond.end())?;
    *z_axis = normalized(*z_axis, bond.begin(), bond.end())?;
    Ok(true)
}

fn frame_of_reference<B: AtropBondEndpoints>(
    bond: &B,
    conformer: AtropisomerConformer<'_>,
) -> Result<Option<Frame>, crate::StereoError> {
    let (mut x, mut y, mut z) = ([0.0; 3], [0.0; 3], [0.0; 3]);
    Ok(
        frame_of_reference_source(bond, conformer, &mut x, &mut y, &mut z)?
            .then_some(Frame { y, z }),
    )
}

fn end_vector_source<G: StereoGraphAccess>(
    topology: &G,
    end: &AtropEnd,
    frame: Frame,
    conformer: AtropisomerConformer<'_>,
    normalize_output: bool,
    bond_vector: &mut [f64; 3],
) -> Result<(), PerceptionError> {
    // Complete pinned source: getAtropIsomerEndVect.
    // RDKit✔️✔️: bool getAtropIsomerEndVect(const AtropAtomAndBondVec &atomAndBondVec,
    // RDKit✔️✔️:                            const RDGeom::Point3D &yAxis,
    // RDKit✔️✔️:                            const RDGeom::Point3D &zAxis, const Conformer *conf,
    // RDKit✔️✔️:                            RDGeom::Point3D &bondVec) {
    // RDKit✔️✔️:   PRECONDITION(
    // RDKit✔️✔️:       atomAndBondVec.second.size() > 0 && atomAndBondVec.second.size() < 3,
    // RDKit✔️✔️:       "bad bond size");
    // RDKit✔️✔️:   PRECONDITION(atomAndBondVec.second[0], "bad first bond");
    // RDKit✔️✔️:   PRECONDITION(atomAndBondVec.second.size() == 1 || atomAndBondVec.second[1],
    // RDKit✔️✔️:                "bad second bond");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   bondVec = conf->getAtomPos(atomAndBondVec.second[0]
    // RDKit✔️✔️:                                  ->getOtherAtom(atomAndBondVec.first)
    // RDKit✔️✔️:                                  ->getIdx()) -
    // RDKit✔️✔️:             conf->getAtomPos(
    // RDKit✔️✔️:                 atomAndBondVec.first->getIdx());  // in old frame of reference
    // RDKit✔️✔️:
    // RDKit✔️✔️:   bondVec = RDGeom::Point3D(0.0, bondVec.dotProduct(yAxis),
    // RDKit✔️✔️:                             bondVec.dotProduct(zAxis));  // in new frame
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // make sure the other atom is on the other side
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (atomAndBondVec.second.size() == 2) {
    // RDKit✔️✔️:     RDGeom::Point3D otherVec =
    // RDKit✔️✔️:         conf->getAtomPos(atomAndBondVec.second[1]
    // RDKit✔️✔️:                              ->getOtherAtom(atomAndBondVec.first)
    // RDKit✔️✔️:                              ->getIdx()) -
    // RDKit✔️✔️:         conf->getAtomPos(
    // RDKit✔️✔️:             atomAndBondVec.first->getIdx());  // in old frame of reference
    // RDKit✔️✔️:     otherVec = RDGeom::Point3D(0.0, otherVec.dotProduct(yAxis),
    // RDKit✔️✔️:                                otherVec.dotProduct(zAxis));  // in new frame
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (bondVec.length() < REALLY_SMALL_BOND_LEN) {
    // RDKit✔️✔️:       bondVec = -otherVec;  // put it on the other side of otherVec
    // RDKit✔️✔️:     } else if (bondVec.dotProduct(otherVec) > REALLY_SMALL_BOND_LEN) {
    // RDKit✔️✔️:       // the product of dotproducts (y-values) should be
    // RDKit✔️✔️:       // negative (or at least zero)
    // RDKit✔️✔️:       BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:           << "Both bonds on one end of an atropisomer are on the same side - atoms is : "
    // RDKit✔️✔️:           << atomAndBondVec.first->getIdx() << std::endl;
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (bondVec.length() < REALLY_SMALL_BOND_LEN) {
    // RDKit✔️✔️:     BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:         << "Could not find a bond on one end of an atropisomer that is not co-linear - atoms are : "
    // RDKit✔️✔️:         << atomAndBondVec.first->getIdx() << std::endl;
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   bondVec.normalize();
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // BEGIN RDKIT CPP complete reached 3D projection branch
    // RDKit❗✔️:     } else {  // the conf is 3D
    // RDKit❗✔️:       // to be considered, one or more neighbor bonds must have a wedge or
    // RDKit❗✔️:       // hash
    // RDKit❗✔️:
    // RDKit❗✔️:       // find the projection of the bond(s) on this end in the frame of
    // RDKit❗✔️:       // reference's  x=0  plane
    // RDKit❗✔️:       RDGeom::Point3D tempBondVec =
    // RDKit❗✔️:           conf->getAtomPos(
    // RDKit❗✔️:               atomAndBondVecs[bondAtomIndex]
    // RDKit❗✔️:                   .second[0]
    // RDKit❗✔️:                   ->getOtherAtom(atomAndBondVecs[bondAtomIndex].first)
    // RDKit❗✔️:                   ->getIdx()) -
    // RDKit❗✔️:           conf->getAtomPos(atomAndBondVecs[bondAtomIndex].first->getIdx());
    // RDKit❗✔️:       bondVecs[bondAtomIndex] = RDGeom::Point3D(
    // RDKit❗✔️:           0.0, tempBondVec.dotProduct(yAxis), tempBondVec.dotProduct(zAxis));
    // RDKit❗✔️:
    // RDKit❗✔️:       if (atomAndBondVecs[bondAtomIndex].second.size() == 2) {
    // RDKit❗✔️:         tempBondVec =
    // RDKit❗✔️:             conf->getAtomPos(
    // RDKit❗✔️:                 atomAndBondVecs[bondAtomIndex]
    // RDKit❗✔️:                     .second[1]
    // RDKit❗✔️:                     ->getOtherAtom(atomAndBondVecs[bondAtomIndex].first)
    // RDKit❗✔️:                     ->getIdx()) -
    // RDKit❗✔️:             conf->getAtomPos(atomAndBondVecs[bondAtomIndex].first->getIdx());
    // RDKit❗✔️:
    // RDKit❗✔️:         // get the projection of the 2nd bond on the x=0 plane
    // RDKit❗✔️:
    // RDKit❗✔️:         RDGeom::Point3D otherBondVec = RDGeom::Point3D(
    // RDKit❗✔️:             0.0, tempBondVec.dotProduct(yAxis), tempBondVec.dotProduct(zAxis));
    // RDKit❗✔️:
    // RDKit❗✔️:         // if the first atom is co-linear with the main atrop bond, use
    // RDKit❗✔️:         // the opposite of the 2nd atom
    // RDKit❗✔️:
    // RDKit❗✔️:         if (bondVecs[bondAtomIndex].length() < REALLY_SMALL_BOND_LEN) {
    // RDKit❗✔️:           bondVecs[bondAtomIndex] =
    // RDKit❗✔️:               -otherBondVec;  // note - it might still be co-linear-
    // RDKit❗✔️:                               // this is checked below
    // RDKit❗✔️:         } else if (bondVecs[bondAtomIndex].dotProduct(otherBondVec) >
    // RDKit❗✔️:                    REALLY_SMALL_BOND_LEN) {
    // RDKit❗✔️:           BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:               << "Both bonds on one end of an atropisomer are on the same side - atoms are: "
    // RDKit❗✔️:               << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx()
    // RDKit❗✔️:               << std::endl;
    // RDKit❗✔️:           return;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:
    // RDKit❗✔️:       if (bondVecs[bondAtomIndex].length() < REALLY_SMALL_BOND_LEN) {
    // RDKit❗✔️:         BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:             << "Failed to find a bond on one end of an atropisomer that is NOT co-linear - atoms are: "
    // RDKit❗✔️:             << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx()
    // RDKit❗✔️:             << std::endl;
    // RDKit❗✔️:         return;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // END RDKIT CPP complete reached 3D projection branch
    // getAtropIsomerEndVect (2D/wedging) normalizes; the 3D detection
    // branch projects the same vectors/checks and retains their lengths.
    // Native normalization exceptions are fatal, never diagnostic fallback.
    // Behavior review: mutate the real projected-vector output in native order.
    // Same-side/collinear false branches and normalization errors retain the
    // preceding raw or opposite-second projection, without clearing/restoring.
    // The shared 3D detection projection preserves its original non-normalized
    // branch; the getAtropIsomerEndVect caller enforces the source size guard.
    // Complexity review: at most two borrowed bonds and coordinate rows, fixed
    // scalar projection arithmetic; no heap buffers or cloned geometry.
    if normalize_output && !(1..=2).contains(&end.bonds.len()) {
        return Err(PerceptionError::CarrierCount {
            atom: end.atom,
            count: end.bonds.len(),
        });
    }
    let raw = |bond_id: BondId| -> Result<[f64; 3], SourceCarrierError> {
        let carrier = &topology.bonds()[bond_id.index()];
        Ok(sub(
            point(conformer, other_atom(carrier, end.atom)?),
            point(conformer, end.atom),
        ))
    };
    *bond_vector = raw(end.bonds[0])?;
    *bond_vector = [0.0, dot(*bond_vector, frame.y), dot(*bond_vector, frame.z)];
    if end.bonds.len() == 2 {
        let other = raw(end.bonds[1])?;
        let other = [0.0, dot(other, frame.y), dot(other, frame.z)];
        if length(*bond_vector) < REALLY_SMALL_BOND_LEN {
            *bond_vector = neg(other);
        } else if dot(*bond_vector, other) > REALLY_SMALL_BOND_LEN {
            return Err(AtropisomerRejectionKind::SameSideCarriers.into());
        }
    }
    if length(*bond_vector) < REALLY_SMALL_BOND_LEN {
        return Err(AtropisomerRejectionKind::CollinearCarrier.into());
    }
    if normalize_output {
        let neighbor = other_atom(&topology.bonds()[end.bonds[0].index()], end.atom)?;
        *bond_vector = normalized(*bond_vector, end.atom, neighbor)?;
    }
    Ok(())
}

fn end_vector(
    topology: &TopologyBlock,
    end: &AtropEnd,
    frame: Frame,
    conformer: AtropisomerConformer<'_>,
    normalize_output: bool,
) -> Result<[f64; 3], PerceptionError> {
    let mut result = [0.0; 3];
    end_vector_source(
        topology,
        end,
        frame,
        conformer,
        normalize_output,
        &mut result,
    )?;
    Ok(result)
}

fn effective_direction<G: StereoGraphAccess>(
    topology: &G,
    updates: &BTreeMap<BondId, AtropisomerWedgeUpdate>,
    bond: BondId,
) -> BondDirection {
    updates
        .get(&bond)
        .map_or(topology.bonds()[bond.index()].direction(), |value| {
            value.direction
        })
}

fn interpreted_end_direction(
    topology: &TopologyBlock,
    end: &AtropEnd,
    updates: &BTreeMap<BondId, AtropisomerWedgeUpdate>,
) -> Result<BondDirection, AtropisomerRejectionKind> {
    // Complete pinned source: getBondDir.
    // RDKit✔️✔️: std::pair<bool, Bond::BondDir> getBondDir(
    // RDKit✔️✔️:     const Bond *bond, const AtropAtomAndBondVec &atomAndBondVec) {
    // RDKit✔️✔️:   // get the wedge dir for this end of the bond
    // RDKit✔️✔️:   // if the first bond 1 has a bondDir, use it
    // RDKit✔️✔️:   // if the second bond has a bond dir use the opposite of if
    // RDKit✔️✔️:   // if both bonds have a dir, make sure they are different
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto bond1Dir = atomAndBondVec.second[0]->getBondDir();
    // RDKit✔️✔️:   if (bond1Dir != Bond::BEGINWEDGE && bond1Dir != Bond::BEGINDASH) {
    // RDKit✔️✔️:     bond1Dir = Bond::NONE;  //  we dont care if it any thing else
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto bond2Dir = atomAndBondVec.second.size() == 2
    // RDKit✔️✔️:                       ? atomAndBondVec.second[1]->getBondDir()
    // RDKit✔️✔️:                       : Bond::NONE;
    // RDKit✔️✔️:   if (bond2Dir != Bond::BEGINWEDGE && bond2Dir != Bond::BEGINDASH) {
    // RDKit✔️✔️:     bond2Dir = Bond::NONE;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // if both are set to a direction, they must NOT be the same - one
    // RDKit✔️✔️:   // must be a dash and the other a hash
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (bond1Dir != Bond::NONE && bond2Dir != Bond::NONE &&
    // RDKit✔️✔️:       bond1Dir == bond2Dir) {
    // RDKit✔️✔️:     BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:         << "The bonds on one end of an atropisomer are both UP or both DOWN - atoms are: "
    // RDKit✔️✔️:         << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx() << std::endl;
    // RDKit✔️✔️:     return {false, Bond::BondDir::NONE};
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (bond1Dir == Bond::BEGINWEDGE || bond2Dir == Bond::BEGINDASH) {
    // RDKit✔️✔️:     return {true, Bond::BondDir::BEGINWEDGE};
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (bond1Dir == Bond::BEGINDASH || bond2Dir == Bond::BEGINWEDGE) {
    // RDKit✔️✔️:     return {true, Bond::BondDir::BEGINDASH};
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return {true, Bond::BondDir::NONE};
    // RDKit✔️✔️: }
    let normalize = |value| match value {
        BondDirection::BeginWedge | BondDirection::BeginDash => value,
        _ => BondDirection::None,
    };
    let first = normalize(effective_direction(topology, updates, end.bonds[0]));
    let second = if end.bonds.len() == 2 {
        normalize(effective_direction(topology, updates, end.bonds[1]))
    } else {
        BondDirection::None
    };
    if first != BondDirection::None && first == second {
        return Err(AtropisomerRejectionKind::InconsistentDirections);
    }
    if first == BondDirection::BeginWedge || second == BondDirection::BeginDash {
        Ok(BondDirection::BeginWedge)
    } else if first == BondDirection::BeginDash || second == BondDirection::BeginWedge {
        Ok(BondDirection::BeginDash)
    } else {
        Ok(BondDirection::None)
    }
}

fn validate_conformer(
    conformer: AtropisomerConformer<'_>,
    atom_count: usize,
) -> Result<(), AtropisomerError> {
    match conformer {
        AtropisomerConformer::TwoD(value) => value
            .validate_for_atom_count(atom_count)
            .map_err(|source| AtropisomerError::InvalidCoordinates { source }),
        AtropisomerConformer::ThreeD(value) => value
            .validate_for_atom_count(atom_count)
            .map_err(|source| AtropisomerError::InvalidCoordinates { source }),
    }
}

fn detect_one(
    topology: &TopologyBlock,
    bond: &Bond,
    conformer: Option<AtropisomerConformer<'_>>,
) -> Result<BondStereo, PerceptionError> {
    // Complete pinned source: DetectAtropisomerChiralityOneBond.
    // RDKit✔️✔️: void DetectAtropisomerChiralityOneBond(Bond *bond, ROMol &mol,
    // RDKit✔️✔️:                                        const Conformer *conf) {
    // RDKit✔️✔️:   // the approach is this:
    // RDKit✔️✔️:   // we will view the system along the line from the potential atropisomer
    // RDKit✔️✔️:   // bond, from atom1 to atom 2 and we do a coordinate transformation to
    // RDKit✔️✔️:   // the plane of reference where that vector, from a1 to a2, is the x-AXIS.
    // RDKit✔️✔️:   // For 2D, the y axis is in the 2D plane, and the zaxis is perpendicaul to
    // RDKit✔️✔️:   // the 2D plane For 3D, the Y and Z axes are taken arbitrarily to form
    // RDKit✔️✔️:   // a right-handed system with the X-axis.
    // RDKit✔️✔️:   //  atoms 1 and 2 each have one or two bonds out from the main potential
    // RDKit✔️✔️:   //  atrop bond. for each end of the main bond, we find a vector to reprent
    // RDKit✔️✔️:   //  the neighbor atom with the smallest index as its projection onto the
    // RDKit✔️✔️:   //  x=0 plane.
    // RDKit✔️✔️:   // (In 2d, this projection is on the y-AXIS for the end that does NOT have
    // RDKit✔️✔️:   // a wedge/hash bond, and  on the z axis - out of the plane - for the end
    // RDKit✔️✔️:   // that does have a wedge/hash). The chirality is recorded as the
    // RDKit✔️✔️:   // direction we rotate from, atom 1's projection to atom2's proejection -
    // RDKit✔️✔️:   // either clockwise or counter clockwise
    // RDKit✔️✔️:
    // RDKit✔️✔️:   PRECONDITION(bond, "bad bond");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // one vector for each end - each one - should end up with 1 or 2 entries
    // RDKit✔️✔️:   AtropAtomAndBondVec atomAndBondVecs[2];
    // RDKit✔️✔️:   if (!getAtropisomerAtomsAndBonds(bond, atomAndBondVecs, mol)) {
    // RDKit✔️✔️:     return;  // not an atropisomer
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // make sure we do not have wiggle bonds
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (auto atomAndBondVec : atomAndBondVecs) {
    // RDKit✔️✔️:     for (auto endBond : atomAndBondVec.second) {
    // RDKit✔️✔️:       if (endBond->getBondDir() == Bond::UNKNOWN) {
    // RDKit✔️✔️:         return;  // not an atropisomer
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // the convention is that in the absence of coords, the coordiates are choosen
    // RDKit✔️✔️:   // with the lowest numbered atom of the atrop bond down, and the other atom
    // RDKit✔️✔️:   // straight up.
    // RDKit✔️✔️:   // On each end, the lowest numbered connecting atom is on the left
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   //              a      b
    // RDKit✔️✔️:   //               \   /
    // RDKit✔️✔️:   //                 c
    // RDKit✔️✔️:   //                 |
    // RDKit✔️✔️:   //                 d
    // RDKit✔️✔️:   //               /   \     aaa
    // RDKit✔️✔️:   //              e      f
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   // where  c > d
    // RDKit✔️✔️:   //        a < b
    // RDKit✔️✔️:   //        e < f
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (conf == nullptr) {
    // RDKit✔️✔️:     std::pair<bool, Bond::BondDir> bond1DirResult;
    // RDKit✔️✔️:     bond1DirResult = getBondDir(bond, atomAndBondVecs[0]);
    // RDKit✔️✔️:     if (!bond1DirResult.first) {
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     std::pair<bool, Bond::BondDir> bond2DirResult;
    // RDKit✔️✔️:     bond2DirResult = getBondDir(bond, atomAndBondVecs[1]);
    // RDKit✔️✔️:     if (!bond2DirResult.first) {
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (bond1DirResult.second == bond2DirResult.second) {
    // RDKit✔️✔️:       BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:           << "inconsistent bond wedging for an atropisomer.  Atoms are: "
    // RDKit✔️✔️:           << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx()
    // RDKit✔️✔️:           << std::endl;
    // RDKit✔️✔️:       return;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (bond1DirResult.second == Bond::BEGINWEDGE ||
    // RDKit✔️✔️:         bond2DirResult.second == Bond::BEGINDASH) {
    // RDKit✔️✔️:       bond->setStereo(Bond::BondStereo::STEREOATROPCCW);
    // RDKit✔️✔️:     } else if (bond1DirResult.second == Bond::BEGINDASH ||
    // RDKit✔️✔️:                bond2DirResult.second == Bond::BEGINWEDGE) {
    // RDKit✔️✔️:       bond->setStereo(Bond::BondStereo::STEREOATROPCW);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // create a frame of reference that has its X-axis along the atrop bond
    // RDKit✔️✔️:
    // RDKit✔️✔️:   RDGeom::Point3D xAxis, yAxis, zAxis;
    // RDKit✔️✔️:   if (!getBondFrameOfReference(bond, conf, xAxis, yAxis, zAxis)) {
    // RDKit✔️✔️:     // connot percieve atroisomer
    // RDKit✔️✔️:     BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:         << "Failed to get a frame of reference along an atropisomer bond - atoms are: "
    // RDKit✔️✔️:         << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx() << std::endl;
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   RDGeom::Point3D bondVecs[2];  // one bond vector from each end of the
    // RDKit✔️✔️:                                 // potential atropisomer bond
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (int bondAtomIndex = 0; bondAtomIndex < 2; ++bondAtomIndex) {
    // RDKit✔️✔️:     // if the conf is 2D, we use the wedge bonds to set the coords for the
    // RDKit✔️✔️:     // projected vector onto the xAxis perpendicular plane (looking down
    // RDKit✔️✔️:     // the atrop bond )
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (!conf->is3D()) {
    // RDKit✔️✔️:       // get the wedge dir for this end of the bond
    // RDKit✔️✔️:       // if the first bond 1 has a bondDir, use it
    // RDKit✔️✔️:       // if the second bond has a bond dir use the opposite of if
    // RDKit✔️✔️:       // if both bonds have a dir, make sure they are different
    // RDKit✔️✔️:
    // RDKit✔️✔️:       std::pair<bool, Bond::BondDir> bondDirResult;
    // RDKit✔️✔️:
    // RDKit✔️✔️:       bondDirResult = getBondDir(bond, atomAndBondVecs[bondAtomIndex]);
    // RDKit✔️✔️:       if (!bondDirResult.first) {
    // RDKit✔️✔️:         return;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (!getAtropIsomerEndVect(atomAndBondVecs[bondAtomIndex], yAxis, zAxis,
    // RDKit✔️✔️:                                  conf, bondVecs[bondAtomIndex])) {
    // RDKit✔️✔️:         return;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (bondDirResult.second == Bond::BEGINWEDGE) {
    // RDKit✔️✔️:         bondVecs[bondAtomIndex].y *= 0.707;
    // RDKit✔️✔️:         bondVecs[bondAtomIndex].z = fabs(bondVecs[bondAtomIndex].y);
    // RDKit✔️✔️:       } else if (bondDirResult.second == Bond::BEGINDASH) {
    // RDKit✔️✔️:         bondVecs[bondAtomIndex].y *= 0.707;
    // RDKit✔️✔️:         bondVecs[bondAtomIndex].z = -fabs(bondVecs[bondAtomIndex].y);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {  // the conf is 3D
    // RDKit✔️✔️:       // to be considered, one or more neighbor bonds must have a wedge or
    // RDKit✔️✔️:       // hash
    // RDKit✔️✔️:
    // RDKit✔️✔️:       // find the projection of the bond(s) on this end in the frame of
    // RDKit✔️✔️:       // reference's  x=0  plane
    // RDKit✔️✔️:       RDGeom::Point3D tempBondVec =
    // RDKit✔️✔️:           conf->getAtomPos(
    // RDKit✔️✔️:               atomAndBondVecs[bondAtomIndex]
    // RDKit✔️✔️:                   .second[0]
    // RDKit✔️✔️:                   ->getOtherAtom(atomAndBondVecs[bondAtomIndex].first)
    // RDKit✔️✔️:                   ->getIdx()) -
    // RDKit✔️✔️:           conf->getAtomPos(atomAndBondVecs[bondAtomIndex].first->getIdx());
    // RDKit✔️✔️:       bondVecs[bondAtomIndex] = RDGeom::Point3D(
    // RDKit✔️✔️:           0.0, tempBondVec.dotProduct(yAxis), tempBondVec.dotProduct(zAxis));
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (atomAndBondVecs[bondAtomIndex].second.size() == 2) {
    // RDKit✔️✔️:         tempBondVec =
    // RDKit✔️✔️:             conf->getAtomPos(
    // RDKit✔️✔️:                 atomAndBondVecs[bondAtomIndex]
    // RDKit✔️✔️:                     .second[1]
    // RDKit✔️✔️:                     ->getOtherAtom(atomAndBondVecs[bondAtomIndex].first)
    // RDKit✔️✔️:                     ->getIdx()) -
    // RDKit✔️✔️:             conf->getAtomPos(atomAndBondVecs[bondAtomIndex].first->getIdx());
    // RDKit✔️✔️:
    // RDKit✔️✔️:         // get the projection of the 2nd bond on the x=0 plane
    // RDKit✔️✔️:
    // RDKit✔️✔️:         RDGeom::Point3D otherBondVec = RDGeom::Point3D(
    // RDKit✔️✔️:             0.0, tempBondVec.dotProduct(yAxis), tempBondVec.dotProduct(zAxis));
    // RDKit✔️✔️:
    // RDKit✔️✔️:         // if the first atom is co-linear with the main atrop bond, use
    // RDKit✔️✔️:         // the opposite of the 2nd atom
    // RDKit✔️✔️:
    // RDKit✔️✔️:         if (bondVecs[bondAtomIndex].length() < REALLY_SMALL_BOND_LEN) {
    // RDKit✔️✔️:           bondVecs[bondAtomIndex] =
    // RDKit✔️✔️:               -otherBondVec;  // note - it might still be co-linear-
    // RDKit✔️✔️:                               // this is checked below
    // RDKit✔️✔️:         } else if (bondVecs[bondAtomIndex].dotProduct(otherBondVec) >
    // RDKit✔️✔️:                    REALLY_SMALL_BOND_LEN) {
    // RDKit✔️✔️:           BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:               << "Both bonds on one end of an atropisomer are on the same side - atoms are: "
    // RDKit✔️✔️:               << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx()
    // RDKit✔️✔️:               << std::endl;
    // RDKit✔️✔️:           return;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (bondVecs[bondAtomIndex].length() < REALLY_SMALL_BOND_LEN) {
    // RDKit✔️✔️:         BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:             << "Failed to find a bond on one end of an atropisomer that is NOT co-linear - atoms are: "
    // RDKit✔️✔️:             << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx()
    // RDKit✔️✔️:             << std::endl;
    // RDKit✔️✔️:         return;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto crossProduct = bondVecs[1].crossProduct(bondVecs[0]);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (crossProduct.x > REALLY_SMALL_BOND_LEN) {
    // RDKit✔️✔️:     bond->setStereo(Bond::BondStereo::STEREOATROPCCW);
    // RDKit✔️✔️:   } else if (crossProduct.x < -REALLY_SMALL_BOND_LEN) {
    // RDKit✔️✔️:     bond->setStereo(Bond::BondStereo::STEREOATROPCW);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:         << "The 2 defining bonds for an atropisomer are co-planar - atoms are: "
    // RDKit✔️✔️:         << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx() << std::endl;
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    let ends =
        atropisomer_ends_fresh(topology, bond)?.ok_or(AtropisomerRejectionKind::MissingCarrier)?;
    if ends
        .iter()
        .flat_map(|end| &end.bonds)
        .any(|id| topology.bonds[id.index()].direction() == BondDirection::Unknown)
    {
        return Err(AtropisomerRejectionKind::UnknownCarrierDirection.into());
    }
    let no_updates = BTreeMap::new();
    if conformer.is_none() {
        let first = interpreted_end_direction(topology, &ends[0], &no_updates)?;
        let second = interpreted_end_direction(topology, &ends[1], &no_updates)?;
        if first == second {
            return Err(AtropisomerRejectionKind::InconsistentDirections.into());
        }
        return if first == BondDirection::BeginWedge || second == BondDirection::BeginDash {
            Ok(BondStereo::AtropCcw)
        } else if first == BondDirection::BeginDash || second == BondDirection::BeginWedge {
            Ok(BondStereo::AtropCw)
        } else {
            Err(AtropisomerRejectionKind::InconsistentDirections.into())
        };
    }
    let conformer = conformer.expect("checked above");
    let frame =
        frame_of_reference(bond, conformer)?.ok_or(AtropisomerRejectionKind::ZeroLengthAxis)?;
    let two_d = !conformer_is_3d(conformer);
    let mut vectors = [[0.0; 3]; 2];
    for (index, end) in ends.iter().enumerate() {
        // Native 2D getBondDir failure occurs before this end's geometry read.
        let direction = if two_d {
            Some(interpreted_end_direction(topology, end, &no_updates)?)
        } else {
            None
        };
        vectors[index] = end_vector(topology, end, frame, conformer, two_d)?;
        if let Some(direction) = direction {
            match direction {
                BondDirection::BeginWedge => {
                    vectors[index][1] *= 0.707;
                    vectors[index][2] = vectors[index][1].abs();
                }
                BondDirection::BeginDash => {
                    vectors[index][1] *= 0.707;
                    vectors[index][2] = -vectors[index][1].abs();
                }
                _ => {}
            }
        }
    }
    let x = cross(vectors[1], vectors[0])[0];
    if x > REALLY_SMALL_BOND_LEN {
        Ok(BondStereo::AtropCcw)
    } else if x < -REALLY_SMALL_BOND_LEN {
        Ok(BondStereo::AtropCw)
    } else {
        Err(AtropisomerRejectionKind::CoplanarCarriers.into())
    }
}

pub fn detect_atropisomer_chirality(
    topology: &TopologyBlock,
    conformer: Option<AtropisomerConformer<'_>>,
) -> Result<AtropisomerAssignment, AtropisomerError> {
    topology
        .validate()
        .map_err(|source| AtropisomerError::InvalidTopology { source })?;
    if let Some(value) = conformer {
        validate_conformer(value, topology.atoms.len())?;
    }
    // Complete pinned source: detectAtropisomerChirality.
    // Source candidate storage is std::set<Bond*> (allocation-address order).
    // BondId ordering below is the existing detached identity order, not a
    // reproduction of those native addresses. A partially dirty, otherwise
    // hybridized cache can change the global-update guard and final chemistry
    // with native pointer order. Retain ❗ pending the source-state decision.
    // Fresh parser caller invariants require their own proof; extra
    // detached scalar result/overlay storage makes the cost axis ❌.
    // RDKit❗❌: void detectAtropisomerChirality(ROMol &mol, const Conformer *conf) {
    // RDKit❗❌:   PRECONDITION(conf == nullptr || &(conf->getOwningMol()) == &mol,
    // RDKit❗❌:                "conformer does not belong to molecule");
    // RDKit❗❌:
    // RDKit❗❌:   std::set<Bond *> bondsToTry;
    // RDKit❗❌:
    // RDKit❗❌:   for (auto bond : mol.bonds()) {
    // RDKit❗❌:     if (canHaveDirection(*bond) &&
    // RDKit❗❌:         (bond->getBondDir() == Bond::BondDir::BEGINDASH ||
    // RDKit❗❌:          bond->getBondDir() == Bond::BondDir::BEGINWEDGE)) {
    // RDKit❗❌:       for (const auto &nbrBond : mol.atomBonds(bond->getBeginAtom())) {
    // RDKit❗❌:         if (nbrBond == bond) {
    // RDKit❗❌:           continue;  // a bond is NOT its own neighbor
    // RDKit❗❌:         }
    // RDKit❗❌:         bondsToTry.insert(nbrBond);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (bondsToTry.empty()) {
    // RDKit❗❌:     return;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // First, do a simple check with TotalDegree to see if any bonds might be
    // RDKit❗❌:   // candidates before doing the expensive hybridization calculation.
    // RDKit❗❌:   bool anyBondPassesDegreeCheck = false;
    // RDKit❗❌:   for (auto bondToTry : bondsToTry) {
    // RDKit❗❌:     if (bondToTry->getBeginAtom()->needsUpdatePropertyCache()) {
    // RDKit❗❌:       bondToTry->getBeginAtom()->updatePropertyCache(false);
    // RDKit❗❌:     }
    // RDKit❗❌:     if (bondToTry->getEndAtom()->needsUpdatePropertyCache()) {
    // RDKit❗❌:       bondToTry->getEndAtom()->updatePropertyCache(false);
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     if (bondToTry->getBondType() == Bond::SINGLE &&
    // RDKit❗❌:         bondToTry->getStereo() != Bond::BondStereo::STEREOANY &&
    // RDKit❗❌:         bondToTry->getBeginAtom()->getTotalDegree() >= 2 &&
    // RDKit❗❌:         bondToTry->getBeginAtom()->getTotalDegree() <= 3 &&
    // RDKit❗❌:         bondToTry->getEndAtom()->getTotalDegree() >= 2 &&
    // RDKit❗❌:         bondToTry->getEndAtom()->getTotalDegree() <= 3) {
    // RDKit❗❌:       anyBondPassesDegreeCheck = true;
    // RDKit❗❌:       break;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (!anyBondPassesDegreeCheck) {
    // RDKit❗❌:     return;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // defer cache update on the whole mol unless we actually have bonds to try
    // RDKit❗❌:   // we need to do an update on the whole mol and not just incident atoms
    // RDKit❗❌:   // because we need to calculate hybridization, which is non-local
    // RDKit❗❌:   bool needsUpdate =
    // RDKit❗❌:       mol.needsUpdatePropertyCache() ||
    // RDKit❗❌:       std::any_of(mol.atoms().begin(), mol.atoms().end(), [](const auto atom) {
    // RDKit❗❌:         return atom->getAtomicNum() != 0 &&
    // RDKit❗❌:                atom->getHybridization() == Atom::HybridizationType::UNSPECIFIED;
    // RDKit❗❌:       });
    // RDKit❗❌:   if (needsUpdate) {
    // RDKit❗❌:     mol.updatePropertyCache(false);
    // RDKit❗❌:     MolOps::setConjugation(mol);
    // RDKit❗❌:     MolOps::setHybridization(mol);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   for (auto bondToTry : bondsToTry) {
    // RDKit❗❌:     if (bondToTry->getBondType() != Bond::SINGLE ||
    // RDKit❗❌:         bondToTry->getStereo() == Bond::BondStereo::STEREOANY ||
    // RDKit❗❌:         // before, we checked only on totalDegree = 2 or 3,
    // RDKit❗❌:         // but this causes false positives for something like a chiral sulfoxide
    // RDKit❗❌:         // since the S is tetrahedral (sp3) but has only 3 substituents.
    // RDKit❗❌:         // the hybridization code relies on totalDegree,
    // RDKit❗❌:         // but modified to include and making sure to include conjugation
    // RDKit❗❌:         // so while this is more expensive per molecule, it is closer to intent
    // RDKit❗❌:         bondToTry->getBeginAtom()->getHybridization() != Atom::SP2 ||
    // RDKit❗❌:         bondToTry->getEndAtom()->getHybridization() != Atom::SP2) {
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     DetectAtropisomerChiralityOneBond(bondToTry, mol, conf);
    // RDKit❗❌:   }
    // RDKit❗❌: }
    let candidates: BTreeSet<_> = source_atropisomer_candidate_membership(topology)
        .into_iter()
        .enumerate()
        .filter_map(|(index, present)| present.then_some(BondId::new(index)))
        .collect();
    let ordered: Vec<_> = candidates.into_iter().collect();
    let traversal = SourceAtropisomerCandidateTraversal::try_from_parts(topology, &ordered)?;
    detect_atropisomer_from_candidate_order(traversal, conformer)
}

/// A graph-bound borrow of the complete source candidate-set traversal.
///
/// The producer must supply the actual iteration order of the native
/// `std::set<Bond *>`, mapped through source bond identity to detached BondId.
/// Validation checks identity, membership, uniqueness and completeness; it
/// cannot authenticate the order's origin. In particular, successful
/// validation of BondId order is not evidence of native pointer order.
/// This input never sorts, fills missing rows, or reads Rust object addresses.
#[derive(Debug, Clone, Copy)]
struct SourceAtropisomerCandidateTraversal<'a> {
    topology: &'a TopologyBlock,
    candidates: &'a [BondId],
}

impl<'a> SourceAtropisomerCandidateTraversal<'a> {
    fn try_from_parts(
        topology: &'a TopologyBlock,
        candidates: &'a [BondId],
    ) -> Result<Self, AtropisomerError> {
        topology
            .validate()
            .map_err(|source| AtropisomerError::InvalidTopology { source })?;
        let membership = source_atropisomer_candidate_membership(topology);
        let mut seen = vec![false; membership.len()];
        for &bond in candidates {
            let Some(is_candidate) = membership.get(bond.index()) else {
                return Err(AtropisomerError::CandidateTraversalOutOfRange {
                    bond,
                    bond_count: membership.len(),
                });
            };
            if seen[bond.index()] {
                return Err(AtropisomerError::CandidateTraversalDuplicate { bond });
            }
            if !is_candidate {
                return Err(AtropisomerError::CandidateTraversalUnexpected { bond });
            }
            seen[bond.index()] = true;
        }
        for (index, is_candidate) in membership.into_iter().enumerate() {
            if is_candidate && !seen[index] {
                return Err(AtropisomerError::CandidateTraversalMissing {
                    bond: BondId::new(index),
                });
            }
        }
        Ok(Self {
            topology,
            candidates,
        })
    }
}

// Only called with structurally validated topology.
fn source_atropisomer_candidate_membership(topology: &TopologyBlock) -> Vec<bool> {
    // This is the source-defined membership pass, independent of set
    // iteration order. A bool row retains no guessed pointer comparator.
    // RDKit✔️❌:   std::set<Bond *> bondsToTry;
    // RDKit✔️❌:
    // RDKit✔️❌:   for (auto bond : mol.bonds()) {
    // RDKit✔️❌:     if (canHaveDirection(*bond) &&
    // RDKit✔️❌:         (bond->getBondDir() == Bond::BondDir::BEGINDASH ||
    // RDKit✔️❌:          bond->getBondDir() == Bond::BondDir::BEGINWEDGE)) {
    // RDKit✔️❌:       for (const auto &nbrBond : mol.atomBonds(bond->getBeginAtom())) {
    // RDKit✔️❌:         if (nbrBond == bond) {
    // RDKit✔️❌:           continue;  // a bond is NOT its own neighbor
    // RDKit✔️❌:         }
    // RDKit✔️❌:         bondsToTry.insert(nbrBond);
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // Only membership, not std::set pointer order, is reproduced here.
    // Cost: two O(E) validation buffers plus structural graph validation
    // exceed the native candidate-only set on sparse candidate inputs.
    let mut membership = vec![false; topology.bonds.len()];
    for marker in &topology.bonds {
        if can_have_direction(marker)
            && matches!(
                marker.direction(),
                BondDirection::BeginDash | BondDirection::BeginWedge
            )
        {
            for neighbor in topology.adjacency.neighbors_of(marker.begin().index()) {
                if neighbor.bond != marker.id() {
                    membership[neighbor.bond.index()] = true;
                }
            }
        }
    }
    membership
}

// The sole processing body borrows its graph from the checked input, so a
// candidate list cannot accidentally be paired with another topology here.
// Its default public collector still lacks genuine native order provenance.
fn detect_atropisomer_from_candidate_order(
    traversal: SourceAtropisomerCandidateTraversal<'_>,
    conformer: Option<AtropisomerConformer<'_>>,
) -> Result<AtropisomerAssignment, AtropisomerError> {
    let SourceAtropisomerCandidateTraversal {
        topology,
        candidates,
    } = traversal;
    if let Some(value) = conformer {
        validate_conformer(value, topology.atoms.len())?;
    }
    // BEGIN COMPLETE PINNED detectAtropisomerChirality candidate processing
    // RDKit❗❌:   if (bondsToTry.empty()) {
    // RDKit❗❌:     return;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // First, do a simple check with TotalDegree to see if any bonds might be
    // RDKit❗❌:   // candidates before doing the expensive hybridization calculation.
    // RDKit❗❌:   bool anyBondPassesDegreeCheck = false;
    // RDKit❗❌:   for (auto bondToTry : bondsToTry) {
    // RDKit❗❌:     if (bondToTry->getBeginAtom()->needsUpdatePropertyCache()) {
    // RDKit❗❌:       bondToTry->getBeginAtom()->updatePropertyCache(false);
    // RDKit❗❌:     }
    // RDKit❗❌:     if (bondToTry->getEndAtom()->needsUpdatePropertyCache()) {
    // RDKit❗❌:       bondToTry->getEndAtom()->updatePropertyCache(false);
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     if (bondToTry->getBondType() == Bond::SINGLE &&
    // RDKit❗❌:         bondToTry->getStereo() != Bond::BondStereo::STEREOANY &&
    // RDKit❗❌:         bondToTry->getBeginAtom()->getTotalDegree() >= 2 &&
    // RDKit❗❌:         bondToTry->getBeginAtom()->getTotalDegree() <= 3 &&
    // RDKit❗❌:         bondToTry->getEndAtom()->getTotalDegree() >= 2 &&
    // RDKit❗❌:         bondToTry->getEndAtom()->getTotalDegree() <= 3) {
    // RDKit❗❌:       anyBondPassesDegreeCheck = true;
    // RDKit❗❌:       break;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (!anyBondPassesDegreeCheck) {
    // RDKit❗❌:     return;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // defer cache update on the whole mol unless we actually have bonds to try
    // RDKit❗❌:   // we need to do an update on the whole mol and not just incident atoms
    // RDKit❗❌:   // because we need to calculate hybridization, which is non-local
    // RDKit❗❌:   bool needsUpdate =
    // RDKit❗❌:       mol.needsUpdatePropertyCache() ||
    // RDKit❗❌:       std::any_of(mol.atoms().begin(), mol.atoms().end(), [](const auto atom) {
    // RDKit❗❌:         return atom->getAtomicNum() != 0 &&
    // RDKit❗❌:                atom->getHybridization() == Atom::HybridizationType::UNSPECIFIED;
    // RDKit❗❌:       });
    // RDKit❗❌:   if (needsUpdate) {
    // RDKit❗❌:     mol.updatePropertyCache(false);
    // RDKit❗❌:     MolOps::setConjugation(mol);
    // RDKit❗❌:     MolOps::setHybridization(mol);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   for (auto bondToTry : bondsToTry) {
    // RDKit❗❌:     if (bondToTry->getBondType() != Bond::SINGLE ||
    // RDKit❗❌:         bondToTry->getStereo() == Bond::BondStereo::STEREOANY ||
    // RDKit❗❌:         // before, we checked only on totalDegree = 2 or 3,
    // RDKit❗❌:         // but this causes false positives for something like a chiral sulfoxide
    // RDKit❗❌:         // since the S is tetrahedral (sp3) but has only 3 substituents.
    // RDKit❗❌:         // the hybridization code relies on totalDegree,
    // RDKit❗❌:         // but modified to include and making sure to include conjugation
    // RDKit❗❌:         // so while this is more expensive per molecule, it is closer to intent
    // RDKit❗❌:         bondToTry->getBeginAtom()->getHybridization() != Atom::SP2 ||
    // RDKit❗❌:         bondToTry->getEndAtom()->getHybridization() != Atom::SP2) {
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     DetectAtropisomerChiralityOneBond(bondToTry, mol, conf);
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // END COMPLETE PINNED detectAtropisomerChirality candidate processing
    let mut assignment = AtropisomerAssignment::default();
    if candidates.is_empty() {
        return Ok(assignment);
    }
    // The native cache is scalar attached state. Borrow the graph and project
    // only those scalars: no topology/conformer/property cloning. The O(V)
    // overlay is additional memory versus source in-place access (known cost).
    let mut valence = crate::ValenceAssignment {
        explicit_valence: topology
            .atoms
            .iter()
            .map(|a| i32::from(a.source_valence_facts().explicit_valence))
            .collect(),
        implicit_hydrogens: topology
            .atoms
            .iter()
            .map(|a| i32::from(a.source_valence_facts().implicit_valence))
            .collect(),
    };
    let mut any_degree_candidate = false;
    for &id in candidates {
        let bond = &topology.bonds[id.index()];
        for endpoint in [bond.begin(), bond.end()] {
            let index = endpoint.index();
            let facts = SourceAtomValenceFacts {
                explicit_valence: valence.explicit_valence[index] as i8,
                implicit_valence: valence.implicit_hydrogens[index] as i8,
            };
            if crate::valence::source_atom_needs_cache_update_with_facts(
                &topology.atoms[index],
                facts,
            ) {
                let (explicit, implicit) = crate::assign_valence_state_for_atom_from_parts(
                    &topology.atoms,
                    &topology.bonds,
                    &topology.adjacency,
                    endpoint,
                    false,
                )?;
                valence.explicit_valence[index] = explicit;
                valence.implicit_hydrogens[index] = implicit;
                assignment.atom_valence_updates.push((
                    endpoint,
                    SourceAtomValenceFacts {
                        explicit_valence: explicit as i8,
                        implicit_valence: implicit as i8,
                    },
                ));
            }
        }
        if bond.order() == BondOrder::Single
            && bond.stereo() != BondStereo::Any
            && source_total_degree(topology, &valence, bond.begin())? >= 2
            && source_total_degree(topology, &valence, bond.begin())? <= 3
            && source_total_degree(topology, &valence, bond.end())? >= 2
            && source_total_degree(topology, &valence, bond.end())? <= 3
        {
            any_degree_candidate = true;
            break;
        }
    }
    if !any_degree_candidate {
        return Ok(assignment);
    }
    let needs_update = topology.atoms.iter().enumerate().any(|(index, atom)| {
        crate::valence::source_atom_needs_cache_update_with_facts(
            atom,
            SourceAtomValenceFacts {
                explicit_valence: valence.explicit_valence[index] as i8,
                implicit_valence: valence.implicit_hydrogens[index] as i8,
            },
        )
    }) || topology.atoms.iter().any(|atom| {
        atom.atomic_number() != 0 && atom.hybridization() == Hybridization::Unspecified
    });
    if needs_update {
        valence = crate::assign_valence_with_options_for_topology(
            topology,
            crate::ValenceModel::RdkitLike,
            false,
        )?;
        for (index, atom) in topology.atoms.iter().enumerate() {
            assignment.atom_valence_updates.push((
                atom.id(),
                SourceAtomValenceFacts {
                    explicit_valence: valence.explicit_valence[index] as i8,
                    implicit_valence: valence.implicit_hydrogens[index] as i8,
                },
            ));
        }
        let conjugated = crate::assign_conjugation_flags(topology, &valence)?;
        assignment.hybridization = Some(crate::assign_hybridization_with_conjugation(
            topology,
            &valence,
            &conjugated,
        )?);
        assignment.conjugated_bonds = Some(conjugated);
    }
    for &id in candidates {
        let bond = &topology.bonds[id.index()];
        let hybridization = |atom: AtomId| {
            assignment
                .hybridization
                .as_ref()
                .map_or(topology.atoms[atom.index()].hybridization(), |a| {
                    a.values[atom.index()]
                })
        };
        if bond.order() != BondOrder::Single
            || bond.stereo() == BondStereo::Any
            || hybridization(bond.begin()) != Hybridization::Sp2
            || hybridization(bond.end()) != Hybridization::Sp2
        {
            continue;
        }
        match detect_one(topology, bond, conformer) {
            Ok(stereo) => assignment
                .bond_updates
                .push(AtropisomerBondUpdate { bond: id, stereo }),
            Err(PerceptionError::Rejected(kind)) => assignment
                .diagnostics
                .push(AtropisomerDiagnostic { bond: id, kind }),
            Err(PerceptionError::SourcePrecondition { message }) => {
                return Err(AtropisomerError::SourcePrecondition { message });
            }
            Err(PerceptionError::Normalization(error)) => {
                return Err(AtropisomerError::Normalization(error));
            }
            Err(PerceptionError::CarrierCount { atom, count }) => {
                return Err(AtropisomerError::CarrierCount { atom, count });
            }
        }
    }
    Ok(assignment)
}

fn effective_stereo(
    topology: &TopologyBlock,
    assignment: &AtropisomerAssignment,
    bond: BondId,
) -> BondStereo {
    assignment
        .bond_updates
        .iter()
        .find(|update| update.bond == bond)
        .map_or(topology.bonds[bond.index()].stereo(), |update| {
            update.stereo
        })
}

pub fn cleanup_atropisomer_stereo_groups(
    topology: &TopologyBlock,
    detected: &AtropisomerAssignment,
) -> Result<StereoGroupAssignment, AtropisomerError> {
    topology
        .validate()
        .map_err(|source| AtropisomerError::InvalidTopology { source })?;
    for update in &detected.bond_updates {
        if update.bond.index() >= topology.bonds.len() {
            return Err(AtropisomerError::AssignmentBondOutOfRange {
                bond: update.bond,
                bond_count: topology.bonds.len(),
            });
        }
        if !matches!(update.stereo, BondStereo::AtropCw | BondStereo::AtropCcw) {
            return Err(AtropisomerError::InvalidAtropisomerStereo {
                bond: update.bond,
                stereo: update.stereo,
            });
        }
    }
    // Complete pinned source: cleanupAtropisomerStereoGroups.
    // RDKit✔️✔️: void cleanupAtropisomerStereoGroups(ROMol &mol) {
    // RDKit✔️✔️:   std::vector<StereoGroup> newsgs;
    // RDKit✔️✔️:   for (auto sg : mol.getStereoGroups()) {
    // RDKit✔️✔️:     std::vector<Atom *> okatoms;
    // RDKit✔️✔️:     std::vector<Bond *> okbonds;
    // RDKit✔️✔️:
    // RDKit✔️✔️:     for (auto atom : sg.getAtoms()) {
    // RDKit✔️✔️:       bool foundAtrop = false;
    // RDKit✔️✔️:       for (auto bndI : boost::make_iterator_range(mol.getAtomBonds(atom))) {
    // RDKit✔️✔️:         auto bond = (mol)[bndI];
    // RDKit✔️✔️:         if (bond->getStereo() == Bond::BondStereo::STEREOATROPCCW ||
    // RDKit✔️✔️:             bond->getStereo() == Bond::BondStereo::STEREOATROPCW) {
    // RDKit✔️✔️:           foundAtrop = true;
    // RDKit✔️✔️:           if (std::find(okbonds.begin(), okbonds.end(), bond) ==
    // RDKit✔️✔️:               okbonds.end()) {
    // RDKit✔️✔️:             okbonds.push_back(bond);
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (!foundAtrop) {
    // RDKit✔️✔️:         okatoms.push_back(atom);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (okbonds.empty()) {
    // RDKit✔️✔️:       newsgs.push_back(sg);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       newsgs.emplace_back(sg.getGroupType(), std::move(okatoms),
    // RDKit✔️✔️:                           std::move(okbonds));
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   mol.setStereoGroups(std::move(newsgs));
    // RDKit✔️✔️: }
    let mut groups = Vec::with_capacity(topology.stereo_groups.len());
    for group in &topology.stereo_groups {
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
        for atom in group.atoms() {
            let mut found = false;
            for neighbor in topology.adjacency.neighbors_of(atom.index()) {
                if matches!(
                    effective_stereo(topology, detected, neighbor.bond),
                    BondStereo::AtropCw | BondStereo::AtropCcw
                ) {
                    found = true;
                    if !bonds.contains(&neighbor.bond) {
                        bonds.push(neighbor.bond);
                    }
                }
            }
            if !found {
                atoms.push(*atom);
            }
        }
        if bonds.is_empty() {
            groups.push(group.clone());
        } else {
            // The pinned three-argument constructor propagates neither the
            // source read ID nor its write ID; this fresh value resets both.
            groups.push(StereoGroup::new(group.kind(), atoms, bonds)?);
        }
    }
    Ok(StereoGroupAssignment { groups })
}

pub fn stereo_group_atom_ids(
    topology: &TopologyBlock,
    group: &StereoGroup,
    wedges: &AtropisomerWedgeAssignment,
) -> Result<Vec<AtomId>, AtropisomerError> {
    validate_stereo_group_topology(topology)?;
    validate_stereo_group_members(topology, group)?;
    let source_map = WedgeAssignments::from_atropisomer_wedge_parts_source(
        &wedges.bond_updates,
        &wedges.source_map_writes,
    );
    let mut atom_ids = Vec::new();
    collect_stereo_group_atom_ids_source(topology, group, &mut atom_ids, &source_map)?;
    Ok(atom_ids)
}

pub fn get_all_atom_ids_for_stereo_group(
    topology: &TopologyBlock,
    group: &StereoGroup,
    wedge_bonds: &WedgeAssignments,
) -> Result<Vec<AtomId>, AtropisomerError> {
    validate_stereo_group_topology(topology)?;
    validate_stereo_group_members(topology, group)?;
    let mut atom_ids = Vec::new();
    collect_stereo_group_atom_ids_source(topology, group, &mut atom_ids, wedge_bonds)?;
    Ok(atom_ids)
}

/// Collects source-ordered atom IDs for each checked stereo group.
///
/// The topology is validated once for the batch. Group membership is checked
/// in input order before each source collection pass.
pub fn get_all_atom_ids_for_stereo_groups(
    topology: &TopologyBlock,
    groups: &[StereoGroup],
    wedge_bonds: &WedgeAssignments,
) -> Result<Vec<Vec<AtomId>>, AtropisomerError> {
    validate_stereo_group_topology(topology)?;

    let mut atom_ids_by_group = Vec::with_capacity(groups.len());
    for group in groups {
        validate_stereo_group_members(topology, group)?;
        let mut atom_ids = Vec::new();
        collect_stereo_group_atom_ids_source(topology, group, &mut atom_ids, wedge_bonds)?;
        atom_ids_by_group.push(atom_ids);
    }
    Ok(atom_ids_by_group)
}

#[cfg(test)]
std::thread_local! {
    static STEREO_GROUP_TOPOLOGY_VALIDATION_CALLS: std::cell::Cell<usize> = const {
        std::cell::Cell::new(0)
    };
}

fn validate_stereo_group_topology(topology: &TopologyBlock) -> Result<(), AtropisomerError> {
    #[cfg(test)]
    STEREO_GROUP_TOPOLOGY_VALIDATION_CALLS.with(|calls| calls.set(calls.get() + 1));

    topology
        .validate()
        .map_err(|source| AtropisomerError::InvalidTopology { source })
}

fn validate_stereo_group_members(
    topology: &TopologyBlock,
    group: &StereoGroup,
) -> Result<(), AtropisomerError> {
    for atom in group.atoms() {
        if atom.index() >= topology.atoms.len() {
            return Err(AtropisomerError::StereoGroupAtomOutOfRange {
                atom: *atom,
                atom_count: topology.atoms.len(),
            });
        }
    }
    for bond in group.bonds() {
        if bond.index() >= topology.bonds.len() {
            return Err(AtropisomerError::StereoGroupBondOutOfRange {
                bond: *bond,
                bond_count: topology.bonds.len(),
            });
        }
    }
    Ok(())
}

/// Full native collector over actual detached graph rows and reusable output.
#[doc(hidden)]
pub fn collect_stereo_group_atom_ids_source(
    topology: &TopologyBlock,
    group: &StereoGroup,
    atom_ids: &mut Vec<AtomId>,
    wedge_bonds: &WedgeAssignments,
) -> Result<(), AtropisomerError> {
    collect_stereo_group_atom_ids_impl(topology, group, atom_ids, wedge_bonds)
}

/// The same source collector over actual borrowed canonical query members.
#[doc(hidden)]
pub fn collect_query_stereo_group_atom_ids_source(
    query: &cosmolkit_model::QueryGraph,
    group: &StereoGroup,
    atom_ids: &mut Vec<AtomId>,
    wedge_bonds: &WedgeAssignments,
) -> Result<(), AtropisomerError> {
    collect_stereo_group_atom_ids_impl(query, group, atom_ids, wedge_bonds)
}

// Closed borrowed member access: this source helper reads only atom IDs,
// carrier bond endpoints/directions and actual adjacency. Query predicates,
// element conversion, valence/cache, geometry and lifecycle are not touched.
trait SourceStereoGroupGraph {
    fn source_atom_count(&self) -> usize;
    fn source_bond_count(&self) -> usize;
    fn source_atom_id(&self, index: usize) -> Option<AtomId>;
    fn source_bond(&self, index: usize) -> Option<&Bond>;
    fn source_adjacent_bonds(&self, index: usize) -> Option<impl Iterator<Item = BondId>>;
}

impl SourceStereoGroupGraph for TopologyBlock {
    fn source_atom_count(&self) -> usize {
        self.atoms.len()
    }
    fn source_bond_count(&self) -> usize {
        self.bonds.len()
    }
    fn source_atom_id(&self, index: usize) -> Option<AtomId> {
        self.atoms.get(index).map(|atom| atom.id())
    }
    fn source_bond(&self, index: usize) -> Option<&Bond> {
        self.bonds.get(index)
    }
    fn source_adjacent_bonds(&self, index: usize) -> Option<impl Iterator<Item = BondId>> {
        self.adjacency
            .try_neighbors_of(index)
            .map(|row| row.iter().map(|edge| edge.bond))
    }
}

impl SourceStereoGroupGraph for cosmolkit_model::QueryGraph {
    fn source_atom_count(&self) -> usize {
        self.num_atoms()
    }
    fn source_bond_count(&self) -> usize {
        self.num_bonds()
    }
    fn source_atom_id(&self, index: usize) -> Option<AtomId> {
        self.atom(index).map(|atom| atom.id())
    }
    fn source_bond(&self, index: usize) -> Option<&Bond> {
        self.bond(index).map(|bond| bond.bond())
    }
    fn source_adjacent_bonds(&self, index: usize) -> Option<impl Iterator<Item = BondId>> {
        self.adjacency()
            .get(index)
            .map(|row| row.iter().map(|(_, bond)| BondId::new(*bond)))
    }
}

fn collect_stereo_group_atom_ids_impl<G: SourceStereoGroupGraph>(
    graph: &G,
    group: &StereoGroup,
    atom_ids: &mut Vec<AtomId>,
    wedge_bonds: &WedgeAssignments,
) -> Result<(), AtropisomerError> {
    // RDKit❗✔️: void getAllAtomIdsForStereoGroup(
    // RDKit❗✔️:     const ROMol &mol, const StereoGroup &group,
    // RDKit❗✔️:     std::vector<unsigned int> &atomIds,
    // RDKit❗✔️:     const std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
    // RDKit❗✔️:         &wedgeBonds) {
    // RDKit❗✔️:   atomIds.clear();
    // RDKit❗✔️:   for (auto &&atom : group.getAtoms()) {
    // RDKit❗✔️:     atomIds.push_back(atom->getIdx());
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   for (auto &&bond : group.getBonds()) {
    // RDKit❗✔️:     // figure out which atoms of the bond get wedge/hash indications
    // RDKit❗✔️:     // mark the atom with the wedge/hash
    // RDKit❗✔️:
    // RDKit❗✔️:     for (auto atom : {bond->getBeginAtom(), bond->getEndAtom()}) {
    // RDKit❗✔️:       for (const auto atomBond : mol.atomBonds(atom)) {
    // RDKit❗✔️:         if (atomBond->getIdx() == bond->getIdx()) {
    // RDKit❗✔️:           continue;
    // RDKit❗✔️:         }
    // RDKit❗✔️:
    // RDKit❗✔️:         if (atomBond->getBondDir() == Bond::BEGINWEDGE ||
    // RDKit❗✔️:             atomBond->getBondDir() == Bond::BEGINDASH ||
    // RDKit❗✔️:             (wedgeBonds.find(atomBond->getIdx()) != wedgeBonds.end() &&
    // RDKit❗✔️:              (wedgeBonds.at(atomBond->getIdx())->getType()) ==
    // RDKit❗✔️:                  Chirality::WedgeInfoType::WedgeInfoTypeAtropisomer)) {
    // RDKit❗✔️:           if (std::find(atomIds.begin(), atomIds.end(), atom->getIdx()) ==
    // RDKit❗✔️:               atomIds.end()) {
    // RDKit❗✔️:             atomIds.push_back(atom->getIdx());
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // RDKit❗✔️: const std::vector<Atom *> &StereoGroup::getAtoms() const { return d_atoms; }
    // RDKit❗✔️: const std::vector<Bond *> &StereoGroup::getBonds() const { return d_bonds; }
    // RDKit❗✔️:   unsigned int getIdx() const { return d_index; }
    // RDKit❗✔️:   unsigned int getIdx() const { return d_index; }
    // RDKit❗✔️:   BondDir getBondDir() const { return static_cast<BondDir>(d_dirTag); }
    // RDKit❗✔️: Atom *Bond::getBeginAtom() const {
    // RDKit❗✔️:   PRECONDITION(dp_mol != nullptr, "no owning molecule for bond");
    // RDKit❗✔️:   return dp_mol->getAtomWithIdx(d_beginAtomIdx);
    // RDKit❗✔️: };
    // RDKit❗✔️: Atom *Bond::getEndAtom() const {
    // RDKit❗✔️:   PRECONDITION(dp_mol != nullptr, "no owning molecule for bond");
    // RDKit❗✔️:   return dp_mol->getAtomWithIdx(d_endAtomIdx);
    // RDKit❗✔️: };
    // RDKit❗✔️: Atom *ROMol::getAtomWithIdx(unsigned int idx) {
    // RDKit❗✔️:   URANGE_CHECK(idx, getNumAtoms());
    // RDKit❗✔️:
    // RDKit❗✔️:   auto vd = boost::vertex(idx, d_graph);
    // RDKit❗✔️:   auto res = d_graph[vd];
    // RDKit❗✔️:   POSTCONDITION(res, "");
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // RDKit❗✔️:
    // RDKit❗✔️: const Atom *ROMol::getAtomWithIdx(unsigned int idx) const {
    // RDKit❗✔️:   URANGE_CHECK(idx, getNumAtoms());
    // RDKit❗✔️:
    // RDKit❗✔️:   auto vd = boost::vertex(idx, d_graph);
    // RDKit❗✔️:   const auto res = d_graph[vd];
    // RDKit❗✔️:
    // RDKit❗✔️:   POSTCONDITION(res, "");
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // RDKit❗✔️:
    // RDKit❗✔️: ROMol::ADJ_ITER_PAIR ROMol::getAtomNeighbors(Atom const *at) const {
    // RDKit❗✔️:   PRECONDITION(at, "no atom");
    // RDKit❗✔️:   PRECONDITION(&at->getOwningMol() == this,
    // RDKit❗✔️:                "atom not associated with this molecule");
    // RDKit❗✔️:   return boost::adjacent_vertices(at->getIdx(), d_graph);
    // RDKit❗✔️: };
    // RDKit❗✔️:
    // RDKit❗✔️: ROMol::OBOND_ITER_PAIR ROMol::getAtomBonds(Atom const *at) const {
    // RDKit❗✔️:   PRECONDITION(at, "no atom");
    // RDKit❗✔️:   PRECONDITION(&at->getOwningMol() == this,
    // RDKit❗✔️:                "atom not associated with this molecule");
    // RDKit❗✔️:   return boost::out_edges(at->getIdx(), d_graph);
    // RDKit❗✔️: }
    // RDKit❗✔️:
    // RDKit❗✔️: ROMol::ATOM_ITER_PAIR ROMol::getVertices() { return boost::vertices(d_graph); }
    // RDKit❗✔️: ROMol::BOND_ITER_PAIR ROMol::getEdges() { return boost::edges(d_graph); }
    // RDKit❗✔️: ROMol::ATOM_ITER_PAIR ROMol::getVertices() const {
    // RDKit❗✔️:   return boost::vertices(d_graph);
    // RDKit❗✔️: }
    // RDKit❗✔️: ROMol::BOND_ITER_PAIR ROMol::getEdges() const { return boost::edges(d_graph); }
    // RDKit❗✔️:
    // RDKit❗✔️: unsigned int ROMol::addAtom(Atom *atom_pin, bool updateLabel,
    // RDKit❗✔️:                             bool takeOwnership) {
    // RDKit❗✔️:   PRECONDITION(atom_pin, "null atom passed in");
    // RDKit❗✔️:   PRECONDITION(!takeOwnership || !atom_pin->hasOwningMol() ||
    // RDKit❗✔️:                    &atom_pin->getOwningMol() == this,
    // RDKit❗✔️:                "cannot take ownership of an atom which already has an owner");
    // RDKit❗✔️:   Atom *atom_p;
    // RDKit❗✔️:   if (!takeOwnership) {
    // RDKit❗✔️:     atom_p = atom_pin->copy();
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     atom_p = atom_pin;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   atom_p->setOwningMol(this);
    // RDKit❗✔️:   auto which = boost::add_vertex(d_graph);
    // RDKit❗✔️:   d_graph[which] = atom_p;
    // RDKit❗✔️:   atom_p->setIdx(which);
    // RDKit❗✔️:   if (updateLabel) {
    // RDKit❗✔️:     replaceAtomBookmark(atom_p, ci_RIGHTMOST_ATOM);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   for (auto &conf : d_confs) {
    // RDKit❗✔️:     conf->setAtomPos(which, RDGeom::Point3D(0.0, 0.0, 0.0));
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return rdcast<unsigned int>(which);
    // RDKit❗✔️: };
    // RDKit❗✔️:   WedgeInfoType getType() const override {
    // RDKit❗✔️:     return Chirality::WedgeInfoType::WedgeInfoTypeChiral;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   WedgeInfoType getType() const override {
    // RDKit❗✔️:     return Chirality::WedgeInfoType::WedgeInfoTypeAtropisomer;
    // RDKit❗✔️:   }
    // Source WedgeInfo getType implementations consume only the actual enum
    // tag; direction, associated axial ID and chiral center are not inspected.
    // Behavior: clear before metadata reads, retain group atom order/duplicates,
    // resolve both endpoint pointers before the endpoint loop, visit actual
    // adjacency order, then append only newly marked endpoint IDs. Structural
    // extraction errors preserve the already-cleared/appended output prefix.
    // Cost: source O(A + sum(degree * (log W + current output length))). No
    // whole-graph validation, topology clone, sorted/deduplicated atom set or
    // extra neighbor scan. Existing checked APIs keep their own validation.
    atom_ids.clear();
    for atom in group.atoms() {
        atom_ids.push(*atom);
    }
    for group_bond_id in group.bonds() {
        let group_bond = graph.source_bond(group_bond_id.index()).ok_or(
            AtropisomerError::StereoGroupBondOutOfRange {
                bond: *group_bond_id,
                bond_count: graph.source_bond_count(),
            },
        )?;
        // C++ initializer-list resolves begin and end before visiting either.
        let begin = graph.source_atom_id(group_bond.begin().index()).ok_or(
            AtropisomerError::StereoGroupAtomOutOfRange {
                atom: group_bond.begin(),
                atom_count: graph.source_atom_count(),
            },
        )?;
        let end = graph.source_atom_id(group_bond.end().index()).ok_or(
            AtropisomerError::StereoGroupAtomOutOfRange {
                atom: group_bond.end(),
                atom_count: graph.source_atom_count(),
            },
        )?;
        for atom in [begin, end] {
            let neighbors = graph.source_adjacent_bonds(atom.index()).ok_or(
                AtropisomerError::SourcePrecondition {
                    message: "source stereo-group atom adjacency row is absent",
                },
            )?;
            for adjacent in neighbors {
                let adjacent_bond = graph.source_bond(adjacent.index()).ok_or(
                    AtropisomerError::StereoGroupBondOutOfRange {
                        bond: adjacent,
                        bond_count: graph.source_bond_count(),
                    },
                )?;
                if adjacent_bond.id() == group_bond.id() {
                    continue;
                }
                if matches!(
                    adjacent_bond.direction(),
                    BondDirection::BeginWedge | BondDirection::BeginDash
                ) || matches!(
                    wedge_bonds.get(adjacent_bond.id()),
                    Some(WedgeInfo::Atropisomer { .. })
                ) {
                    if !atom_ids.contains(&atom) {
                        atom_ids.push(atom);
                    }
                }
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod stereo_group_batch_tests {
    use super::*;
    use crate::{WedgeAssignments, pick_bonds_to_wedge};
    use cosmolkit_model::{Atom, AtomSpec, BondSpec, StereoGroupKind, TopologyValidationError};
    use cosmolkit_types::{BondOrder, ChiralTag, Element};

    fn group_topology(directions: [BondDirection; 3]) -> TopologyBlock {
        let atoms = (0..4)
            .map(|index| {
                Atom::from_spec(
                    AtomId::new(index),
                    AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::Unspecified),
                )
            })
            .collect();
        let edges = [(0, 1), (1, 2), (1, 3)];
        let bonds = edges
            .into_iter()
            .zip(directions)
            .enumerate()
            .map(|(index, ((begin, end), direction))| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single)
                        .with_direction(direction),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("valid stereo-group test topology")
    }

    fn tetrahedral_wedge_topology() -> TopologyBlock {
        let atoms = (0..5)
            .map(|index| {
                let tag = if index == 0 {
                    ChiralTag::TetrahedralCw
                } else {
                    ChiralTag::Unspecified
                };
                Atom::from_spec(
                    AtomId::new(index),
                    AtomSpec::new(Element::C).with_chiral_tag(tag),
                )
            })
            .collect();
        let bonds = (0..4)
            .map(|index| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(0), AtomId::new(index + 1), BondOrder::Single),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("valid tetrahedral wedge topology")
    }

    fn reset_validation_count() {
        STEREO_GROUP_TOPOLOGY_VALIDATION_CALLS.with(|calls| calls.set(0));
    }

    fn validation_count() -> usize {
        STEREO_GROUP_TOPOLOGY_VALIDATION_CALLS.with(std::cell::Cell::get)
    }

    #[test]
    fn cf_smi_groups_batch_preserves_member_order_duplicates_and_overlaps() {
        let topology = group_topology([
            BondDirection::None,
            BondDirection::BeginWedge,
            BondDirection::BeginDash,
        ]);
        let groups = vec![
            StereoGroup::new(
                StereoGroupKind::Absolute,
                vec![AtomId::new(2), AtomId::new(0)],
                Vec::new(),
            )
            .expect("valid distinct stereo members"),
            StereoGroup::new(
                StereoGroupKind::Or,
                vec![AtomId::new(0), AtomId::new(1)],
                Vec::new(),
            )
            .expect("valid distinct stereo members"),
            StereoGroup::new(
                StereoGroupKind::And,
                vec![AtomId::new(3)],
                vec![BondId::new(0)],
            )
            .expect("valid distinct stereo members"),
            StereoGroup::new(
                StereoGroupKind::Or,
                vec![AtomId::new(1)],
                vec![BondId::new(0)],
            )
            .expect("valid distinct stereo members"),
        ];

        let expected = vec![
            vec![AtomId::new(2), AtomId::new(0)],
            vec![AtomId::new(0), AtomId::new(1)],
            vec![AtomId::new(3), AtomId::new(1)],
            vec![AtomId::new(1)],
        ];
        let batch =
            get_all_atom_ids_for_stereo_groups(&topology, &groups, &WedgeAssignments::default())
                .expect("collect all valid groups");
        let independent = groups
            .iter()
            .map(|group| {
                get_all_atom_ids_for_stereo_group(&topology, group, &WedgeAssignments::default())
                    .expect("collect one valid group")
            })
            .collect::<Vec<_>>();

        assert_eq!(batch, expected);
        assert_eq!(batch, independent);
    }

    #[test]
    fn cf_smi_groups_batch_handles_empty_and_single_group() {
        let topology = group_topology([BondDirection::None; 3]);
        let empty =
            get_all_atom_ids_for_stereo_groups(&topology, &[], &WedgeAssignments::default())
                .expect("empty checked batch");
        assert!(empty.is_empty());

        let group = StereoGroup::new(
            StereoGroupKind::Absolute,
            vec![AtomId::new(3), AtomId::new(0)],
            Vec::new(),
        )
        .expect("valid distinct stereo members");
        assert_eq!(
            get_all_atom_ids_for_stereo_groups(
                &topology,
                std::slice::from_ref(&group),
                &WedgeAssignments::default(),
            )
            .expect("one checked group"),
            vec![vec![AtomId::new(3), AtomId::new(0)]]
        );
    }

    #[test]
    fn cf_smi_groups_wedge_map_accepts_only_atropisomer_entries() {
        let topology = group_topology([BondDirection::None; 3]);
        let axial_group =
            StereoGroup::new(StereoGroupKind::Absolute, Vec::new(), vec![BondId::new(0)])
                .expect("valid distinct stereo members");
        let update = AtropisomerWedgeUpdate {
            bond: BondId::new(1),
            begin: AtomId::new(1),
            end: AtomId::new(2),
            direction: BondDirection::None,
            atropisomer_bond: BondId::new(0),
        };
        let atrop_map =
            WedgeAssignments::from_atropisomer_wedge_assignment(AtropisomerWedgeAssignment {
                source_map_writes: vec![update.bond],
                bond_updates: vec![update],
                diagnostics: Vec::new(),
            });
        assert_eq!(
            get_all_atom_ids_for_stereo_group(&topology, &axial_group, &atrop_map)
                .expect("Atropisomer map entry expands its endpoint"),
            vec![AtomId::new(1)]
        );

        let tetrahedral = tetrahedral_wedge_topology();
        let chiral_map = pick_bonds_to_wedge(&tetrahedral, None)
            .expect("valid tetrahedral topology receives a Chiral wedge map entry");
        let (chiral_bond, entry) = chiral_map
            .iter()
            .find(|(_, entry)| matches!(entry, WedgeInfo::Chiral { .. }))
            .expect("a tetrahedral center receives a Chiral map entry");
        assert!(matches!(entry, WedgeInfo::Chiral { .. }));
        let group_bond = (0..4)
            .map(BondId::new)
            .find(|bond| *bond != chiral_bond)
            .expect("a second central bond remains for the group");
        let chiral_group =
            StereoGroup::new(StereoGroupKind::Absolute, Vec::new(), vec![group_bond])
                .expect("valid distinct stereo members");
        assert_eq!(
            get_all_atom_ids_for_stereo_group(&tetrahedral, &chiral_group, &chiral_map)
                .expect("Chiral map entries do not expand an enhanced group"),
            Vec::<AtomId>::new()
        );
    }

    #[test]
    fn cf_smi_groups_batch_preserves_typed_errors_and_precedence() {
        let topology = group_topology([BondDirection::None; 3]);
        let invalid_members = StereoGroup::new(
            StereoGroupKind::Absolute,
            vec![AtomId::new(99)],
            vec![BondId::new(99)],
        )
        .expect("valid distinct stereo members");
        assert_eq!(
            get_all_atom_ids_for_stereo_groups(
                &topology,
                std::slice::from_ref(&invalid_members),
                &WedgeAssignments::default(),
            ),
            Err(AtropisomerError::StereoGroupAtomOutOfRange {
                atom: AtomId::new(99),
                atom_count: 4,
            })
        );

        let invalid_bond = StereoGroup::new(StereoGroupKind::Or, Vec::new(), vec![BondId::new(99)])
            .expect("valid distinct stereo members");
        assert_eq!(
            get_all_atom_ids_for_stereo_groups(
                &topology,
                &[
                    StereoGroup::new(StereoGroupKind::Absolute, vec![AtomId::new(0)], Vec::new(),)
                        .expect("valid distinct stereo members"),
                    invalid_bond
                ],
                &WedgeAssignments::default(),
            ),
            Err(AtropisomerError::StereoGroupBondOutOfRange {
                bond: BondId::new(99),
                bond_count: 3,
            })
        );

        let mut invalid_topology = group_topology([BondDirection::None; 3]);
        invalid_topology.bonds[0] = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(99), AtomId::new(1), BondOrder::Single),
        );
        assert!(matches!(
            get_all_atom_ids_for_stereo_groups(
                &invalid_topology,
                &[StereoGroup::new(
                    StereoGroupKind::Absolute,
                    vec![AtomId::new(0)],
                    Vec::new(),
                ).expect("valid distinct stereo members")],
                &WedgeAssignments::default(),
            ),
            Err(AtropisomerError::InvalidTopology {
                source: TopologyValidationError::BondEndpointOutOfRange {
                    bond,
                    ..
                }
            })
                if bond == BondId::new(0)
        ));
    }

    #[test]
    fn cf_smi_groups_batch_validates_topology_once_per_call() {
        let topology = group_topology([BondDirection::None; 3]);
        let groups = (0..12)
            .map(|index| {
                StereoGroup::new(
                    StereoGroupKind::Or,
                    vec![AtomId::new(index % 4)],
                    Vec::new(),
                )
                .expect("valid distinct stereo members")
            })
            .collect::<Vec<_>>();

        reset_validation_count();
        get_all_atom_ids_for_stereo_groups(&topology, &groups, &WedgeAssignments::default())
            .expect("batch of twelve groups");
        assert_eq!(validation_count(), 1);

        reset_validation_count();
        get_all_atom_ids_for_stereo_group(&topology, &groups[0], &WedgeAssignments::default())
            .expect("checked singular group");
        assert_eq!(validation_count(), 1);

        reset_validation_count();
        get_all_atom_ids_for_stereo_groups(&topology, &[], &WedgeAssignments::default())
            .expect("empty checked batch still validates the topology");
        assert_eq!(validation_count(), 1);
    }
}

fn no_conf_direction(
    stereo: BondStereo,
    which_end: usize,
    which_bond: usize,
) -> Result<BondDirection, PerceptionError> {
    // Complete pinned source: getBondDirForAtropisomerNoConf.
    // RDKit✔️✔️: Bond::BondDir getBondDirForAtropisomerNoConf(Bond::BondStereo bondStereo,
    // RDKit✔️✔️:                                              unsigned int whichEnd,
    // RDKit✔️✔️:                                              unsigned int whichBond) {
    // RDKit✔️✔️:   // the convention is that in the absence of coords, the coordiates are choosen
    // RDKit✔️✔️:   // with the lowest numbered atom of the atrop bond down, and the other atom
    // RDKit✔️✔️:   // straight up.
    // RDKit✔️✔️:   // On each end, the lowest numbered connecting atom is on the left
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   //              a      b
    // RDKit✔️✔️:   //               \   /
    // RDKit✔️✔️:   //                 c
    // RDKit✔️✔️:   //                 |
    // RDKit✔️✔️:   //                 d
    // RDKit✔️✔️:   //               /   \     aaa
    // RDKit✔️✔️:   //              e      f
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   // where  c > d
    // RDKit✔️✔️:   //        a < b
    // RDKit✔️✔️:   //        e < f
    // RDKit✔️✔️:
    // RDKit✔️✔️:   PRECONDITION(whichEnd <= 1, "whichEnd must be 0 or 1");
    // RDKit✔️✔️:   PRECONDITION(whichBond <= 1, "whichBond must be 0 or 1");
    // RDKit✔️✔️:   PRECONDITION(bondStereo == Bond::BondStereo::STEREOATROPCW ||
    // RDKit✔️✔️:                    bondStereo == Bond::BondStereo::STEREOATROPCCW,
    // RDKit✔️✔️:                "bondStereo must be BondAtropisomerCW or BondAtropisomerCCW");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   int flips = 0;
    // RDKit✔️✔️:   if (bondStereo == Bond::BondStereo::STEREOATROPCW) {
    // RDKit✔️✔️:     ++flips;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (whichBond == 1) {
    // RDKit✔️✔️:     ++flips;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (whichEnd == 1) {
    // RDKit✔️✔️:     ++flips;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return flips % 2 ? Bond::BEGINDASH : Bond::BEGINWEDGE;
    // RDKit✔️✔️: }
    // Behavior review: preserve source preconditions in their original order;
    // only 0/1 end and carrier indices and CW/CCW stereo reach the parity sum.
    // Complexity review: scalar guards and three additions, O(1), no allocation.
    if which_end > 1 {
        return Err(PerceptionError::SourcePrecondition {
            message: "whichEnd must be 0 or 1",
        });
    }
    if which_bond > 1 {
        return Err(PerceptionError::SourcePrecondition {
            message: "whichBond must be 0 or 1",
        });
    }
    if !matches!(stereo, BondStereo::AtropCw | BondStereo::AtropCcw) {
        return Err(PerceptionError::SourcePrecondition {
            message: "bondStereo must be BondAtropisomerCW or BondAtropisomerCCW",
        });
    }
    let flips = usize::from(stereo == BondStereo::AtropCw)
        + usize::from(which_bond == 1)
        + usize::from(which_end == 1);
    Ok(if flips % 2 == 1 {
        BondDirection::BeginDash
    } else {
        BondDirection::BeginWedge
    })
}

fn two_d_direction(
    vectors: [[f64; 3]; 2],
    stereo: BondStereo,
    which_end: usize,
    which_bond: usize,
) -> Result<BondDirection, PerceptionError> {
    // Complete pinned source: getBondDirForAtropisomer2d.
    // RDKit✔️✔️: Bond::BondDir getBondDirForAtropisomer2d(RDGeom::Point3D bondVecs[2],
    // RDKit✔️✔️:                                          Bond::BondStereo bondStereo,
    // RDKit✔️✔️:                                          unsigned int whichEnd,
    // RDKit✔️✔️:                                          unsigned int whichBond) {
    // RDKit✔️✔️:   PRECONDITION(whichEnd <= 1, "whichEnd must be 0 or 1");
    // RDKit✔️✔️:   PRECONDITION(whichBond <= 1, "whichBond must be 0 or 1");
    // RDKit✔️✔️:   PRECONDITION(bondStereo == Bond::BondStereo::STEREOATROPCW ||
    // RDKit✔️✔️:                    bondStereo == Bond::BondStereo::STEREOATROPCCW,
    // RDKit✔️✔️:                "bondStereo must be BondAtropisomerCW or BondAtropisomerCCW");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   int flips = 0;
    // RDKit✔️✔️:   if (bondStereo == Bond::BondStereo::STEREOATROPCCW) {
    // RDKit✔️✔️:     ++flips;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (whichBond == 1) {
    // RDKit✔️✔️:     ++flips;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (whichEnd == 1) {
    // RDKit✔️✔️:     ++flips;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (bondVecs[1 - whichEnd].y < 0) {
    // RDKit✔️✔️:     ++flips;  // if the OTHER end is negative for the low index bond vec, it
    // RDKit✔️✔️:               // is a flip
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return flips % 2 ? Bond::BEGINWEDGE : Bond::BEGINDASH;
    // RDKit✔️✔️: }
    // Behavior review: source preconditions precede both the stereo parity
    // and the other-end Y read; negative zero and NaN do not satisfy Y < 0.
    // Complexity review: constant scalar checks and one indexed vector read,
    // no allocation, geometry reconstruction or graph traversal.
    if which_end > 1 {
        return Err(PerceptionError::SourcePrecondition {
            message: "whichEnd must be 0 or 1",
        });
    }
    if which_bond > 1 {
        return Err(PerceptionError::SourcePrecondition {
            message: "whichBond must be 0 or 1",
        });
    }
    if !matches!(stereo, BondStereo::AtropCw | BondStereo::AtropCcw) {
        return Err(PerceptionError::SourcePrecondition {
            message: "bondStereo must be BondAtropisomerCW or BondAtropisomerCCW",
        });
    }
    let flips = usize::from(stereo == BondStereo::AtropCcw)
        + usize::from(which_bond == 1)
        + usize::from(which_end == 1)
        + usize::from(vectors[1 - which_end][1] < 0.0);
    Ok(if flips % 2 == 1 {
        BondDirection::BeginWedge
    } else {
        BondDirection::BeginDash
    })
}

fn three_d_direction<B: AtropBondEndpoints>(
    bond: &B,
    conformer: AtropisomerConformer<'_>,
) -> BondDirection {
    // Complete pinned source: getBondDirForAtropisomer3d.
    // RDKit✔️✔️: Bond::BondDir getBondDirForAtropisomer3d(Bond *whichBond,
    // RDKit✔️✔️:                                          const Conformer *conf) {
    // RDKit✔️✔️:   // for 3D we mark it as wedge or hash depending on the z-value of the bond
    // RDKit✔️✔️:   // vector
    // RDKit✔️✔️:   //  IT really doesn't matter since we ignore these except as MARKERS for
    // RDKit✔️✔️:   //  which bonds are atropisomer bonds
    // RDKit✔️✔️:   if ((conf->getAtomPos(whichBond->getEndAtom()->getIdx()).z -
    // RDKit✔️✔️:        conf->getAtomPos(whichBond->getBeginAtom()->getIdx()).z) >
    // RDKit✔️✔️:       REALLY_SMALL_BOND_LEN) {
    // RDKit✔️✔️:     return Bond::BondDir::BEGINWEDGE;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return Bond::BondDir::BEGINDASH;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Behavior review: actual current end Z minus begin Z, strict threshold;
    // no stereo/tag/flag adjustment. NaN and all values <= threshold are dash.
    // Complexity review: two borrowed coordinate reads and scalar arithmetic,
    // constant stack work with no allocation or geometry reconstruction.
    if point(conformer, bond.end())[2] - point(conformer, bond.begin())[2] > REALLY_SMALL_BOND_LEN {
        BondDirection::BeginWedge
    } else {
        BondDirection::BeginDash
    }
}

// Closed source-state access: graph and native wedge map are separate stores.
trait AtropWedgeState {
    type Graph: StereoGraphAccess;
    fn topology(&self) -> &Self::Graph;
    fn begin(&self, bond: BondId) -> AtomId;
    fn end(&self, bond: BondId) -> AtomId;
    fn direction(&self, bond: BondId) -> BondDirection;
    fn occupied(&self, bond: BondId) -> bool;
    fn set_direction(&mut self, bond: BondId, direction: BondDirection, axial: BondId);
    fn orient(&mut self, bond: BondId, begin: AtomId, axial: BondId);
    fn insert_wedge(&mut self, bond: BondId, axial: BondId);
    fn source_warning(&mut self, prefix: &'static str, axial: BondId) {
        // The defining Native stream insertion is anchored at each source
        // branch below; source warning bytes precede its false return.
        eprintln!(
            "{prefix} {} {}",
            self.begin(axial).index(),
            self.end(axial).index()
        );
    }
}

struct MutableAtropWedgeState<'a, G: StereoGraphMut> {
    topology: &'a mut G,
    wedges: &'a mut WedgeAssignments,
}
impl<G: StereoGraphMut> AtropWedgeState for MutableAtropWedgeState<'_, G> {
    type Graph = G;
    fn topology(&self) -> &G {
        self.topology
    }
    fn begin(&self, bond: BondId) -> AtomId {
        self.topology.bonds()[bond.index()].begin()
    }
    fn end(&self, bond: BondId) -> AtomId {
        self.topology.bonds()[bond.index()].end()
    }
    fn direction(&self, bond: BondId) -> BondDirection {
        self.topology.bonds()[bond.index()].direction()
    }
    fn occupied(&self, bond: BondId) -> bool {
        self.wedges.get(bond).is_some()
    }
    fn set_direction(&mut self, bond: BondId, direction: BondDirection, _axial: BondId) {
        self.topology
            .source_bond_mut(bond.index())
            .set_direction(direction);
    }
    fn orient(&mut self, bond: BondId, begin: AtomId, _axial: BondId) {
        let previous_begin = self.begin(bond);
        if previous_begin != begin {
            self.topology
                .source_bond_mut(bond.index())
                .set_endpoints(begin, previous_begin);
        }
    }
    fn insert_wedge(&mut self, bond: BondId, axial: BondId) {
        self.wedges
            .insert_atropisomer_source(AtropisomerWedgeUpdate {
                bond,
                begin: self.begin(bond),
                end: self.end(bond),
                direction: self.direction(bond),
                atropisomer_bond: axial,
            });
    }
}

struct ProjectedAtropWedgeState<'a, G: StereoGraphAccess> {
    topology: &'a G,
    occupied: &'a BTreeSet<BondId>,
    updates: &'a mut BTreeMap<BondId, AtropisomerWedgeUpdate>,
    map_writes: &'a mut BTreeSet<BondId>,
}
impl<G: StereoGraphAccess> ProjectedAtropWedgeState<'_, G> {
    fn update(&self, bond: BondId, axial: BondId) -> AtropisomerWedgeUpdate {
        self.updates
            .get(&bond)
            .copied()
            .unwrap_or_else(|| AtropisomerWedgeUpdate {
                bond,
                begin: self.topology.bonds()[bond.index()].begin(),
                end: self.topology.bonds()[bond.index()].end(),
                direction: self.topology.bonds()[bond.index()].direction(),
                atropisomer_bond: axial,
            })
    }
}
impl<G: StereoGraphAccess> AtropWedgeState for ProjectedAtropWedgeState<'_, G> {
    type Graph = G;
    fn topology(&self) -> &G {
        self.topology
    }
    fn begin(&self, bond: BondId) -> AtomId {
        self.updates.get(&bond).map_or_else(
            || self.topology.bonds()[bond.index()].begin(),
            |update| update.begin,
        )
    }
    fn end(&self, bond: BondId) -> AtomId {
        self.updates.get(&bond).map_or_else(
            || self.topology.bonds()[bond.index()].end(),
            |update| update.end,
        )
    }
    fn direction(&self, bond: BondId) -> BondDirection {
        effective_direction(self.topology, self.updates, bond)
    }
    fn occupied(&self, bond: BondId) -> bool {
        self.occupied.contains(&bond) || self.map_writes.contains(&bond)
    }
    fn set_direction(&mut self, bond: BondId, direction: BondDirection, axial: BondId) {
        let mut update = self.update(bond, axial);
        update.direction = direction;
        self.updates.insert(bond, update);
    }
    fn orient(&mut self, bond: BondId, begin: AtomId, axial: BondId) {
        if self.begin(bond) != begin {
            let mut update = self.update(bond, axial);
            update.end = update.begin;
            update.begin = begin;
            self.updates.insert(bond, update);
        }
    }
    fn insert_wedge(&mut self, bond: BondId, axial: BondId) {
        let mut update = self.update(bond, axial);
        update.atropisomer_bond = axial;
        self.updates.insert(bond, update);
        self.map_writes.insert(bond);
    }
}

fn wedge_one_no_conformer_source<S: AtropWedgeState>(
    state: &mut S,
    rings: &RingInfo,
    axial_id: BondId,
) -> Result<(), PerceptionError> {
    // RDKit❗❌: bool WedgeBondFromAtropisomerOneBondNoConf(
    // RDKit❗❌:     Bond *bond, const ROMol &mol,
    // RDKit❗❌:     std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
    // RDKit❗❌:         &wedgeBonds) {
    // RDKit❗❌:   PRECONDITION(bond, "no bond");
    // RDKit❗❌:
    // RDKit❗❌:   AtropAtomAndBondVec atomAndBondVecs[2];
    // RDKit❗❌:   if (!getAtropisomerAtomsAndBonds(bond, atomAndBondVecs, mol)) {
    // RDKit❗❌:     return false;  // not an atropisomer
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   //  make sure we do not have wiggle bonds
    // RDKit❗❌:
    // RDKit❗❌:   for (auto atomAndBondVec : atomAndBondVecs) {
    // RDKit❗❌:     for (auto endBond : atomAndBondVec.second) {
    // RDKit❗❌:       if (endBond->getBondDir() == Bond::UNKNOWN) {
    // RDKit❗❌:         return false;  // not an atropisomer)
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // first see if any candidate bond is already set to a wedge or hash
    // RDKit❗❌:   // if so, we will use that bond as a wedge or hash
    // RDKit❗❌:
    // RDKit❗❌:   std::vector<int> useBondsAtEnd[2];
    // RDKit❗❌:   bool foundBondDir = false;
    // RDKit❗❌:
    // RDKit❗❌:   for (unsigned int whichEnd = 0; whichEnd < 2; ++whichEnd) {
    // RDKit❗❌:     for (unsigned int whichBond = 0;
    // RDKit❗❌:          whichBond < atomAndBondVecs[whichEnd].second.size(); ++whichBond) {
    // RDKit❗❌:       auto bondDir = atomAndBondVecs[whichEnd].second[whichBond]->getBondDir();
    // RDKit❗❌:
    // RDKit❗❌:       // see if it is a wedge or hash and its origin is the atom in the
    // RDKit❗❌:       // main bond
    // RDKit❗❌:
    // RDKit❗❌:       if ((bondDir == Bond::BEGINWEDGE || bondDir == Bond::BEGINDASH) &&
    // RDKit❗❌:           atomAndBondVecs[whichEnd].second[whichBond]->getBeginAtom() ==
    // RDKit❗❌:               atomAndBondVecs[whichEnd].first &&
    // RDKit❗❌:           canHaveDirection(*bond)) {
    // RDKit❗❌:         useBondsAtEnd[whichEnd].push_back(whichBond);
    // RDKit❗❌:         foundBondDir = true;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (foundBondDir) {
    // RDKit❗❌:     for (unsigned int whichEnd = 0; whichEnd < 2; ++whichEnd) {
    // RDKit❗❌:       for (unsigned int whichBondIndex = 0;
    // RDKit❗❌:            whichBondIndex < useBondsAtEnd[whichEnd].size(); ++whichBondIndex) {
    // RDKit❗❌:         atomAndBondVecs[whichEnd]
    // RDKit❗❌:             .second[useBondsAtEnd[whichEnd][whichBondIndex]]
    // RDKit❗❌:             ->setBondDir(getBondDirForAtropisomerNoConf(
    // RDKit❗❌:                 bond->getStereo(), whichEnd,
    // RDKit❗❌:                 useBondsAtEnd[whichEnd][whichBondIndex]));
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     return true;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // did not find a good bond dir - pick one to use
    // RDKit❗❌:   // we would like to have one that is not in a ring, and will be a wedge
    // RDKit❗❌:
    // RDKit❗❌:   const RingInfo *ri = bond->getOwningMol().getRingInfo();
    // RDKit❗❌:
    // RDKit❗❌:   int bestBondEnd = -1, bestBondNumber = -1;
    // RDKit❗❌:   bool bestBondIsSingle = false;
    // RDKit❗❌:   unsigned int bestRingCount = INT_MAX;
    // RDKit❗❌:   Bond::BondDir bestBondDir = Bond::BondDir::NONE;
    // RDKit❗❌:   for (unsigned int whichEnd = 0; whichEnd < 2; ++whichEnd) {
    // RDKit❗❌:     for (unsigned int whichBond = 0;
    // RDKit❗❌:          whichBond < atomAndBondVecs[whichEnd].second.size(); ++whichBond) {
    // RDKit❗❌:       auto bondToTry = atomAndBondVecs[whichEnd].second[whichBond];
    // RDKit❗❌:
    // RDKit❗❌:       if (!canHaveDirection(*bondToTry) ||
    // RDKit❗❌:           wedgeBonds.find(bondToTry->getIdx()) != wedgeBonds.end()) {
    // RDKit❗❌:         continue;  // must be a single OR aromatic bond and not already
    // RDKit❗❌:                    // spoken for by a chiral center
    // RDKit❗❌:       }
    // RDKit❗❌:
    // RDKit❗❌:       if (bondToTry->getBondDir() != Bond::BondDir::NONE) {
    // RDKit❗❌:         if (bondToTry->getBeginAtom()->getIdx() ==
    // RDKit❗❌:             atomAndBondVecs[whichEnd].first->getIdx()) {
    // RDKit❗❌:           BOOST_LOG(rdWarningLog)
    // RDKit❗❌:               << "Wedge or hash bond found on atropisomer where not expected - atoms are: "
    // RDKit❗❌:               << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx()
    // RDKit❗❌:               << std::endl;
    // RDKit❗❌:           return false;
    // RDKit❗❌:         } else {
    // RDKit❗❌:           continue;  // wedge or hash bond affecting the OTHER atom
    // RDKit❗❌:                      // = perhaps a chiral center
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       auto ringCount = ri->numBondRings(bondToTry->getIdx());
    // RDKit❗❌:       if (ringCount > bestRingCount) {
    // RDKit❗❌:         continue;
    // RDKit❗❌:       }
    // RDKit❗❌:
    // RDKit❗❌:       else if (ringCount < bestRingCount) {
    // RDKit❗❌:         bestBondEnd = whichEnd;
    // RDKit❗❌:         bestBondNumber = whichBond;
    // RDKit❗❌:         bestRingCount = ringCount;
    // RDKit❗❌:         bestBondIsSingle = (bondToTry->getBondType() == Bond::BondType::SINGLE);
    // RDKit❗❌:         bestBondDir = getBondDirForAtropisomerNoConf(bond->getStereo(),
    // RDKit❗❌:                                                      whichEnd, whichBond);
    // RDKit❗❌:       } else if (bestBondIsSingle &&
    // RDKit❗❌:                  bondToTry->getBondType() != Bond::BondType::SINGLE) {
    // RDKit❗❌:         continue;
    // RDKit❗❌:
    // RDKit❗❌:       } else if (!bestBondIsSingle &&
    // RDKit❗❌:                  bondToTry->getBondType() == Bond::BondType::SINGLE) {
    // RDKit❗❌:         bestBondEnd = whichEnd;
    // RDKit❗❌:         bestBondNumber = whichBond;
    // RDKit❗❌:         bestRingCount = ringCount;
    // RDKit❗❌:         bestBondIsSingle = true;
    // RDKit❗❌:         bestBondDir = getBondDirForAtropisomerNoConf(bond->getStereo(),
    // RDKit❗❌:                                                      whichEnd, whichBond);
    // RDKit❗❌:
    // RDKit❗❌:       } else {
    // RDKit❗❌:         auto bondDir = getBondDirForAtropisomerNoConf(bond->getStereo(),
    // RDKit❗❌:                                                       whichEnd, whichBond);
    // RDKit❗❌:         if (bestBondDir == Bond::BondDir::NONE ||
    // RDKit❗❌:             (bestBondDir == Bond::BondDir::BEGINDASH &&
    // RDKit❗❌:              bondDir == Bond::BondDir::BEGINWEDGE)) {
    // RDKit❗❌:           bestBondEnd = whichEnd;
    // RDKit❗❌:           bestBondNumber = whichBond;
    // RDKit❗❌:           bestRingCount = ringCount;
    // RDKit❗❌:           bestBondIsSingle =
    // RDKit❗❌:               (bondToTry->getBondType() == Bond::BondType::SINGLE);
    // RDKit❗❌:           bestBondDir = bondDir;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (bestBondEnd >= 0)  // we found a good one
    // RDKit❗❌:   {
    // RDKit❗❌:     // make sure the atoms on the bond are in the right order for the
    // RDKit❗❌:     // wedge/hash the atom on the end of the main bond must be listed
    // RDKit❗❌:     // first for the wedge/has bond
    // RDKit❗❌:
    // RDKit❗❌:     auto bestBond = atomAndBondVecs[bestBondEnd].second[bestBondNumber];
    // RDKit❗❌:     if (bestBond->getBeginAtom() != atomAndBondVecs[bestBondEnd].first) {
    // RDKit❗❌:       bestBond->setEndAtom(bestBond->getBeginAtom());
    // RDKit❗❌:       bestBond->setBeginAtom(atomAndBondVecs[bestBondEnd].first);
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     bestBond->setBondDir(bestBondDir);
    // RDKit❗❌:
    // RDKit❗❌:     auto newWedgeInfo = std::unique_ptr<RDKit::Chirality::WedgeInfoBase>(
    // RDKit❗❌:         new RDKit::Chirality::WedgeInfoAtropisomer(bond->getIdx(),
    // RDKit❗❌:                                                    bestBondDir));
    // RDKit❗❌:     wedgeBonds[bestBond->getIdx()] = std::move(newWedgeInfo);
    // RDKit❗❌:   } else {
    // RDKit❗❌:     BOOST_LOG(rdWarningLog)
    // RDKit❗❌:         << "Failed to find a good bond to set as UP or DOWN for an atropisomer - atoms are: "
    // RDKit❗❌:         << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx() << std::endl;
    // RDKit❗❌:     return false;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   return true;
    // RDKit❗❌: }
    // Behavior review: graph direction/orientation writes and native map writes
    // have independent carriers. Existing wedges receive only direction writes;
    // candidate selection executes each source guard before direction lookup.
    // Real mutable state and the detached projection run this sole kernel.
    // Complexity review: two degree-sized carrier vectors, two use-index vectors,
    // borrowed bond/ring access and std::map-shaped BTree lookups; no graph clone.
    // The detached projection adds O(log U) overlay reads/writes where native
    // direct Bond pointers have O(1) getters; that carrier cost is marked ❌.
    let axial = &state.topology().bonds()[axial_id.index()];
    let stereo = axial.stereo();
    let axial_can_have_direction = can_have_direction(axial);
    let current_axial = CurrentAtropBondEndpoints {
        id: axial_id,
        begin: state.begin(axial_id),
        end: state.end(axial_id),
    };
    let ends = atropisomer_ends_fresh(state.topology(), &current_axial)?
        .ok_or(AtropisomerRejectionKind::MissingCarrier)?;
    for end in &ends {
        for &carrier in &end.bonds {
            if state.direction(carrier) == BondDirection::Unknown {
                return Err(AtropisomerRejectionKind::UnknownCarrierDirection.into());
            }
        }
    }
    let mut use_bonds = [Vec::new(), Vec::new()];
    let mut found_direction = false;
    for (which_end, end) in ends.iter().enumerate() {
        for (which_bond, &carrier) in end.bonds.iter().enumerate() {
            let direction = state.direction(carrier);
            if matches!(
                direction,
                BondDirection::BeginWedge | BondDirection::BeginDash
            ) && state.begin(carrier) == end.atom
                && axial_can_have_direction
            {
                use_bonds[which_end].push(which_bond);
                found_direction = true;
            }
        }
    }
    if found_direction {
        for which_end in 0..2 {
            for &which_bond in &use_bonds[which_end] {
                let direction = no_conf_direction(stereo, which_end, which_bond)?;
                state.set_direction(ends[which_end].bonds[which_bond], direction, axial_id);
            }
        }
        return Ok(());
    }
    let mut best: Option<(usize, usize)> = None;
    let mut best_single = false;
    let mut best_ring_count = i32::MAX as u32;
    let mut best_direction = BondDirection::None;
    for (which_end, end) in ends.iter().enumerate() {
        for (which_bond, &carrier) in end.bonds.iter().enumerate() {
            let candidate = &state.topology().bonds()[carrier.index()];
            if !can_have_direction(candidate) || state.occupied(carrier) {
                continue;
            }
            if state.direction(carrier) != BondDirection::None {
                if state.begin(carrier) == end.atom {
                    state.source_warning(
                        "Wedge or hash bond found on atropisomer where not expected - atoms are:",
                        axial_id,
                    );
                    return Err(AtropisomerRejectionKind::DirectionConflict.into());
                }
                continue;
            }
            if !rings.is_initialized() {
                return Err(PerceptionError::SourcePrecondition {
                    message: "RingInfo not initialized",
                });
            }
            let ring_count = rings.num_bond_rings(carrier) as u32;
            let single = candidate.order() == BondOrder::Single;
            if ring_count > best_ring_count {
                continue;
            } else if ring_count < best_ring_count {
                best = Some((which_end, which_bond));
                best_ring_count = ring_count;
                best_single = single;
                best_direction = no_conf_direction(stereo, which_end, which_bond)?;
            } else if best_single && !single {
                continue;
            } else if !best_single && single {
                best = Some((which_end, which_bond));
                best_ring_count = ring_count;
                best_single = true;
                best_direction = no_conf_direction(stereo, which_end, which_bond)?;
            } else {
                let direction = no_conf_direction(stereo, which_end, which_bond)?;
                if best_direction == BondDirection::None
                    || (best_direction == BondDirection::BeginDash
                        && direction == BondDirection::BeginWedge)
                {
                    best = Some((which_end, which_bond));
                    best_ring_count = ring_count;
                    best_single = single;
                    best_direction = direction;
                }
            }
        }
    }
    if let Some((which_end, which_bond)) = best {
        let carrier = ends[which_end].bonds[which_bond];
        state.orient(carrier, ends[which_end].atom, axial_id);
        state.set_direction(carrier, best_direction, axial_id);
        state.insert_wedge(carrier, axial_id);
        Ok(())
    } else {
        {
            state.source_warning(
                "Failed to find a good bond to set as UP or DOWN for an atropisomer - atoms are:",
                axial_id,
            );
            Err(AtropisomerRejectionKind::NoUsableWedgeBond.into())
        }
    }
}

/// Applies the native no-conformer single-atropisomer behavior to actual detached state.
#[doc(hidden)]
pub fn wedge_atropisomer_no_conformer_source(
    topology: &mut TopologyBlock,
    rings: &RingInfo,
    axial_bond: BondId,
    wedges: &mut WedgeAssignments,
) -> Result<bool, AtropisomerError> {
    topology
        .validate()
        .map_err(|source| AtropisomerError::InvalidTopology { source })?;
    if axial_bond.index() >= topology.bonds.len() {
        return Err(AtropisomerError::AxialBondOutOfRange {
            bond: axial_bond,
            bond_count: topology.bonds.len(),
        });
    }
    let mut state = MutableAtropWedgeState { topology, wedges };
    match wedge_one_no_conformer_source(&mut state, rings, axial_bond) {
        Ok(()) => Ok(true),
        Err(PerceptionError::Rejected(_)) => Ok(false),
        Err(PerceptionError::SourcePrecondition { message }) => {
            Err(AtropisomerError::SourcePrecondition { message })
        }
        Err(PerceptionError::Normalization(error)) => Err(AtropisomerError::Normalization(error)),
        Err(PerceptionError::CarrierCount { atom, count }) => {
            Err(AtropisomerError::CarrierCount { atom, count })
        }
    }
}

fn wedge_one_two_d_source<S: AtropWedgeState>(
    state: &mut S,
    rings: &RingInfo,
    axial_id: BondId,
    conformer: AtropisomerConformer<'_>,
) -> Result<(), PerceptionError> {
    // RDKit❗❌: bool WedgeBondFromAtropisomerOneBond2d(
    // RDKit❗❌:     Bond *bond, const ROMol &mol, const Conformer *conf,
    // RDKit❗❌:     std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
    // RDKit❗❌:         &wedgeBonds) {
    // RDKit❗❌:   PRECONDITION(bond, "no bond");
    // RDKit❗❌:
    // RDKit❗❌:   AtropAtomAndBondVec atomAndBondVecs[2];
    // RDKit❗❌:   if (!getAtropisomerAtomsAndBonds(bond, atomAndBondVecs, mol)) {
    // RDKit❗❌:     return false;  // not an atropisomer
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   //  make sure we do not have wiggle bonds
    // RDKit❗❌:
    // RDKit❗❌:   for (auto atomAndBondVec : atomAndBondVecs) {
    // RDKit❗❌:     for (auto endBond : atomAndBondVec.second) {
    // RDKit❗❌:       if (endBond->getBondDir() == Bond::UNKNOWN) {
    // RDKit❗❌:         return false;  // not an atropisomer)
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // create a frame of reference that has its X-axis along the atrop bond
    // RDKit❗❌:
    // RDKit❗❌:   RDGeom::Point3D xAxis, yAxis, zAxis;
    // RDKit❗❌:
    // RDKit❗❌:   if (!getBondFrameOfReference(bond, conf, xAxis, yAxis, zAxis)) {
    // RDKit❗❌:     // connot percieve atroisomer bond
    // RDKit❗❌:
    // RDKit❗❌:     BOOST_LOG(rdWarningLog)
    // RDKit❗❌:         << "Cound not get a frame of reference for an atropisomer bond - atoms are: "
    // RDKit❗❌:         << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx() << std::endl;
    // RDKit❗❌:     return false;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   RDGeom::Point3D bondVecs[2];  // one bond vector from each end of the
    // RDKit❗❌:                                 // potential atropisome bond
    // RDKit❗❌:
    // RDKit❗❌:   for (int bondAtomIndex = 0; bondAtomIndex < 2; ++bondAtomIndex) {
    // RDKit❗❌:     // find a vector to represent the lowest numbered atom on each end
    // RDKit❗❌:     // this vector is NOT the bond vector, but is y-value in the bond
    // RDKit❗❌:     // frame or reference
    // RDKit❗❌:
    // RDKit❗❌:     if (!getAtropIsomerEndVect(atomAndBondVecs[bondAtomIndex], yAxis, zAxis,
    // RDKit❗❌:                                conf, bondVecs[bondAtomIndex])) {
    // RDKit❗❌:       return false;
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     if (bondVecs[bondAtomIndex].length() < REALLY_SMALL_BOND_LEN) {
    // RDKit❗❌:       // did not find a non-colinear bond
    // RDKit❗❌:
    // RDKit❗❌:       BOOST_LOG(rdWarningLog)
    // RDKit❗❌:           << "Failed to get a representative vector for the defining bond of an atropisomer - atoms are: "
    // RDKit❗❌:           << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx()
    // RDKit❗❌:           << std::endl;
    // RDKit❗❌:       return false;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // first see if any candidate bond is already set to a wedge or hash
    // RDKit❗❌:   // if so, we will use that bond as a wedge or hash
    // RDKit❗❌:
    // RDKit❗❌:   std::vector<int> useBondsAtEnd[2];
    // RDKit❗❌:   bool foundBondDir = false;
    // RDKit❗❌:
    // RDKit❗❌:   for (unsigned int whichEnd = 0; whichEnd < 2; ++whichEnd) {
    // RDKit❗❌:     for (unsigned int whichBond = 0;
    // RDKit❗❌:          whichBond < atomAndBondVecs[whichEnd].second.size(); ++whichBond) {
    // RDKit❗❌:       auto bondDir = atomAndBondVecs[whichEnd].second[whichBond]->getBondDir();
    // RDKit❗❌:
    // RDKit❗❌:       // see if it is a wedge or hash and its origin is the atom in the
    // RDKit❗❌:       // main bond
    // RDKit❗❌:
    // RDKit❗❌:       if ((bondDir == Bond::BEGINWEDGE || bondDir == Bond::BEGINDASH) &&
    // RDKit❗❌:           atomAndBondVecs[whichEnd].second[whichBond]->getBeginAtom() ==
    // RDKit❗❌:               atomAndBondVecs[whichEnd].first &&
    // RDKit❗❌:           canHaveDirection(*bond)) {
    // RDKit❗❌:         useBondsAtEnd[whichEnd].push_back(whichBond);
    // RDKit❗❌:         foundBondDir = true;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (foundBondDir) {
    // RDKit❗❌:     for (unsigned int whichEnd = 0; whichEnd < 2; ++whichEnd) {
    // RDKit❗❌:       for (unsigned int whichBondIndex = 0;
    // RDKit❗❌:            whichBondIndex < useBondsAtEnd[whichEnd].size(); ++whichBondIndex) {
    // RDKit❗❌:         atomAndBondVecs[whichEnd]
    // RDKit❗❌:             .second[useBondsAtEnd[whichEnd][whichBondIndex]]
    // RDKit❗❌:             ->setBondDir(getBondDirForAtropisomer2d(
    // RDKit❗❌:                 bondVecs, bond->getStereo(), whichEnd,
    // RDKit❗❌:                 useBondsAtEnd[whichEnd][whichBondIndex]));
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     return true;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // did not find a good bond dir - pick one to use
    // RDKit❗❌:   // we would like to have one that is in a ring, and will favor it being a
    // RDKit❗❌:   // wedge
    // RDKit❗❌:
    // RDKit❗❌:   // We favor rings here because wedging non-ring bonds makes it too likely that
    // RDKit❗❌:   // we'll end up accidentally creating new atropisomeric bonds. This was github
    // RDKit❗❌:   // issue 7371
    // RDKit❗❌:
    // RDKit❗❌:   const RingInfo *ri = bond->getOwningMol().getRingInfo();
    // RDKit❗❌:
    // RDKit❗❌:   int bestBondEnd = -1, bestBondNumber = -1;
    // RDKit❗❌:   bool bestBondIsSingle = false;
    // RDKit❗❌:   unsigned int bestRingCount = INT_MAX;
    // RDKit❗❌:   unsigned int largestRingSize = 0;
    // RDKit❗❌:   Bond::BondDir bestBondDir = Bond::BondDir::NONE;
    // RDKit❗❌:   for (unsigned int whichEnd = 0; whichEnd < 2; ++whichEnd) {
    // RDKit❗❌:     for (unsigned int whichBond = 0;
    // RDKit❗❌:          whichBond < atomAndBondVecs[whichEnd].second.size(); ++whichBond) {
    // RDKit❗❌:       auto bondToTry = atomAndBondVecs[whichEnd].second[whichBond];
    // RDKit❗❌:
    // RDKit❗❌:       if (!canHaveDirection(*bondToTry) ||
    // RDKit❗❌:           wedgeBonds.find(bondToTry->getIdx()) != wedgeBonds.end()) {
    // RDKit❗❌:         continue;  // must be a single OR aromatic bond and not already
    // RDKit❗❌:                    // spoken for by a chiral center
    // RDKit❗❌:       }
    // RDKit❗❌:
    // RDKit❗❌:       if (bondToTry->getBondDir() != Bond::BondDir::NONE) {
    // RDKit❗❌:         if (bondToTry->getBeginAtom()->getIdx() ==
    // RDKit❗❌:             atomAndBondVecs[whichEnd].first->getIdx()) {
    // RDKit❗❌:           if (bondToTry->getBondDir() == Bond::BEGINWEDGE ||
    // RDKit❗❌:               bondToTry->getBondDir() == Bond::BEGINDASH) {
    // RDKit❗❌:             BOOST_LOG(rdWarningLog)
    // RDKit❗❌:                 << "Wedge or hash bond found on atropisomer where not expected - atoms are: "
    // RDKit❗❌:                 << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx()
    // RDKit❗❌:                 << std::endl;
    // RDKit❗❌:             return false;
    // RDKit❗❌:           } else {
    // RDKit❗❌:             continue;  // probably a slash up or down for a double bond
    // RDKit❗❌:           }
    // RDKit❗❌:         } else {
    // RDKit❗❌:           continue;  // wedge or hash bond affecting the OTHER atom
    // RDKit❗❌:                      // = perhaps a chiral center
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       auto ringCount = ri->numBondRings(bondToTry->getIdx());
    // RDKit❗❌:       unsigned int ringSize = 0;
    // RDKit❗❌:       if (!ringCount) {
    // RDKit❗❌:         ringCount = 10;
    // RDKit❗❌:       } else {
    // RDKit❗❌:         // we're going to prefer to put wedges in larger rings, but don't want
    // RDKit❗❌:         // to end up wedging macrocyles if it's avoidable.
    // RDKit❗❌:         ringSize = ri->minBondRingSize(bondToTry->getIdx());
    // RDKit❗❌:         if (ringSize > 8) {
    // RDKit❗❌:           ringSize = 0;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       if (ringCount > bestRingCount) {
    // RDKit❗❌:         continue;
    // RDKit❗❌:       } else if (ringCount < bestRingCount || ringSize > largestRingSize) {
    // RDKit❗❌:         bestBondEnd = whichEnd;
    // RDKit❗❌:         bestBondNumber = whichBond;
    // RDKit❗❌:         bestRingCount = ringCount;
    // RDKit❗❌:         largestRingSize = ringSize;
    // RDKit❗❌:         bestBondIsSingle = (bondToTry->getBondType() == Bond::BondType::SINGLE);
    // RDKit❗❌:         bestBondDir = getBondDirForAtropisomer2d(bondVecs, bond->getStereo(),
    // RDKit❗❌:                                                  whichEnd, whichBond);
    // RDKit❗❌:       } else if (bestBondIsSingle &&
    // RDKit❗❌:                  bondToTry->getBondType() != Bond::BondType::SINGLE) {
    // RDKit❗❌:         continue;
    // RDKit❗❌:
    // RDKit❗❌:       } else if (!bestBondIsSingle &&
    // RDKit❗❌:                  bondToTry->getBondType() == Bond::BondType::SINGLE) {
    // RDKit❗❌:         bestBondEnd = whichEnd;
    // RDKit❗❌:         bestBondNumber = whichBond;
    // RDKit❗❌:         bestRingCount = ringCount;
    // RDKit❗❌:         bestBondIsSingle = true;
    // RDKit❗❌:         bestBondDir = getBondDirForAtropisomer2d(bondVecs, bond->getStereo(),
    // RDKit❗❌:                                                  whichEnd, whichBond);
    // RDKit❗❌:
    // RDKit❗❌:       } else {
    // RDKit❗❌:         auto bondDir = getBondDirForAtropisomer2d(bondVecs, bond->getStereo(),
    // RDKit❗❌:                                                   whichEnd, whichBond);
    // RDKit❗❌:         if (bestBondDir == Bond::BondDir::NONE ||
    // RDKit❗❌:             (bestBondDir == Bond::BondDir::BEGINDASH &&
    // RDKit❗❌:              bondDir == Bond::BondDir::BEGINWEDGE)) {
    // RDKit❗❌:           bestBondEnd = whichEnd;
    // RDKit❗❌:           bestBondNumber = whichBond;
    // RDKit❗❌:           bestRingCount = ringCount;
    // RDKit❗❌:           bestBondIsSingle =
    // RDKit❗❌:               (bondToTry->getBondType() == Bond::BondType::SINGLE);
    // RDKit❗❌:           bestBondDir = bondDir;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (bestBondEnd >= 0) {
    // RDKit❗❌:     // we found a good one
    // RDKit❗❌:     // make sure the atoms on the bond are in the right order for the
    // RDKit❗❌:     // wedge/hash the atom on the end of the main bond must be listed
    // RDKit❗❌:     // first for the wedge/has bond
    // RDKit❗❌:
    // RDKit❗❌:     auto bestBond = atomAndBondVecs[bestBondEnd].second[bestBondNumber];
    // RDKit❗❌:     if (bestBond->getBeginAtom() != atomAndBondVecs[bestBondEnd].first) {
    // RDKit❗❌:       bestBond->setEndAtom(bestBond->getBeginAtom());
    // RDKit❗❌:       bestBond->setBeginAtom(atomAndBondVecs[bestBondEnd].first);
    // RDKit❗❌:     }
    // RDKit❗❌:     bestBond->setBondDir(bestBondDir);
    // RDKit❗❌:
    // RDKit❗❌:     auto newWedgeInfo = std::unique_ptr<RDKit::Chirality::WedgeInfoBase>(
    // RDKit❗❌:         new RDKit::Chirality::WedgeInfoAtropisomer(bond->getIdx(),
    // RDKit❗❌:                                                    bestBondDir));
    // RDKit❗❌:     wedgeBonds[bestBond->getIdx()] = std::move(newWedgeInfo);
    // RDKit❗❌:
    // RDKit❗❌:   } else {
    // RDKit❗❌:     BOOST_LOG(rdWarningLog)
    // RDKit❗❌:         << "Failed to find a good bond to set as UP or DOWN for an atropisomer - atoms are: "
    // RDKit❗❌:         << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx() << std::endl;
    // RDKit❗❌:     return false;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   return true;
    // RDKit❗❌: }
    // Behavior review: literal source guard/mutation order and separate graph
    // direction versus wedge-map writes. Axis and frame use actual current
    // endpoint getters, including earlier orientation writes in a projection.
    // Complexity review: native-sized carrier/use-index vectors and scalar
    // geometry. BTree map/set transport matches membership; the projection
    // adds logarithmic overlay queries/writes to native constant-time getters.
    let axial = &state.topology().bonds()[axial_id.index()];
    let stereo = axial.stereo();
    let axial_can_have_direction = can_have_direction(axial);
    let current_axial = CurrentAtropBondEndpoints {
        id: axial_id,
        begin: state.begin(axial_id),
        end: state.end(axial_id),
    };
    let ends = atropisomer_ends_fresh(state.topology(), &current_axial)?
        .ok_or(AtropisomerRejectionKind::MissingCarrier)?;
    for end in &ends {
        for &carrier in &end.bonds {
            if state.direction(carrier) == BondDirection::Unknown {
                return Err(AtropisomerRejectionKind::UnknownCarrierDirection.into());
            }
        }
    }
    let Some(frame) = frame_of_reference(&current_axial, conformer)? else {
        state.source_warning(
            "Cound not get a frame of reference for an atropisomer bond - atoms are:",
            axial_id,
        );
        return Err(AtropisomerRejectionKind::ZeroLengthAxis.into());
    };
    let mut vectors = [[0.0; 3]; 2];
    for index in 0..2 {
        end_vector_source(
            state.topology(),
            &ends[index],
            frame,
            conformer,
            true,
            &mut vectors[index],
        )?;
        if length(vectors[index]) < REALLY_SMALL_BOND_LEN {
            state.source_warning("Failed to get a representative vector for the defining bond of an atropisomer - atoms are:",axial_id);
            return Err(AtropisomerRejectionKind::CollinearCarrier.into());
        }
    }
    let mut use_bonds = [Vec::new(), Vec::new()];
    let mut found_direction = false;
    for (which_end, end) in ends.iter().enumerate() {
        for (which_bond, &carrier) in end.bonds.iter().enumerate() {
            if matches!(
                state.direction(carrier),
                BondDirection::BeginWedge | BondDirection::BeginDash
            ) && state.begin(carrier) == end.atom
                && axial_can_have_direction
            {
                use_bonds[which_end].push(which_bond);
                found_direction = true;
            }
        }
    }
    if found_direction {
        for which_end in 0..2 {
            for &which_bond in &use_bonds[which_end] {
                let direction = two_d_direction(vectors, stereo, which_end, which_bond)?;
                state.set_direction(ends[which_end].bonds[which_bond], direction, axial_id);
            }
        }
        return Ok(());
    }
    let mut best: Option<(usize, usize)> = None;
    let mut best_single = false;
    let mut best_ring_count = i32::MAX as u32;
    let mut largest_ring_size = 0;
    let mut best_direction = BondDirection::None;
    for (which_end, end) in ends.iter().enumerate() {
        for (which_bond, &carrier) in end.bonds.iter().enumerate() {
            let candidate = &state.topology().bonds()[carrier.index()];
            if !can_have_direction(candidate) || state.occupied(carrier) {
                continue;
            }
            let current_direction = state.direction(carrier);
            if current_direction != BondDirection::None {
                if state.begin(carrier) == end.atom
                    && matches!(
                        current_direction,
                        BondDirection::BeginWedge | BondDirection::BeginDash
                    )
                {
                    state.source_warning(
                        "Wedge or hash bond found on atropisomer where not expected - atoms are:",
                        axial_id,
                    );
                    return Err(AtropisomerRejectionKind::DirectionConflict.into());
                }
                continue;
            }
            if !rings.is_initialized() {
                return Err(PerceptionError::SourcePrecondition {
                    message: "RingInfo not initialized",
                });
            }
            let mut ring_count = rings.num_bond_rings(carrier) as u32;
            let mut ring_size = 0;
            if ring_count == 0 {
                ring_count = 10;
            } else {
                ring_size = rings.min_bond_ring_size(carrier) as u32;
                if ring_size > 8 {
                    ring_size = 0;
                }
            }
            let single = candidate.order() == BondOrder::Single;
            if ring_count > best_ring_count {
                continue;
            } else if ring_count < best_ring_count || ring_size > largest_ring_size {
                best = Some((which_end, which_bond));
                best_ring_count = ring_count;
                largest_ring_size = ring_size;
                best_single = single;
                best_direction = two_d_direction(vectors, stereo, which_end, which_bond)?;
            } else if best_single && !single {
                continue;
            } else if !best_single && single {
                best = Some((which_end, which_bond));
                best_ring_count = ring_count;
                best_single = true;
                best_direction = two_d_direction(vectors, stereo, which_end, which_bond)?;
            } else {
                let direction = two_d_direction(vectors, stereo, which_end, which_bond)?;
                if best_direction == BondDirection::None
                    || (best_direction == BondDirection::BeginDash
                        && direction == BondDirection::BeginWedge)
                {
                    best = Some((which_end, which_bond));
                    best_ring_count = ring_count;
                    best_single = single;
                    best_direction = direction;
                }
            }
        }
    }
    if let Some((which_end, which_bond)) = best {
        let carrier = ends[which_end].bonds[which_bond];
        state.orient(carrier, ends[which_end].atom, axial_id);
        state.set_direction(carrier, best_direction, axial_id);
        state.insert_wedge(carrier, axial_id);
        Ok(())
    } else {
        {
            state.source_warning(
                "Failed to find a good bond to set as UP or DOWN for an atropisomer - atoms are:",
                axial_id,
            );
            Err(AtropisomerRejectionKind::NoUsableWedgeBond.into())
        }
    }
}

/// Applies the native 2D single-atropisomer behavior to actual detached state.
#[doc(hidden)]
pub fn wedge_atropisomer_two_d_source(
    topology: &mut TopologyBlock,
    rings: &RingInfo,
    axial_bond: BondId,
    conformer: AtropisomerConformer<'_>,
    wedges: &mut WedgeAssignments,
) -> Result<bool, AtropisomerError> {
    topology
        .validate()
        .map_err(|source| AtropisomerError::InvalidTopology { source })?;
    if axial_bond.index() >= topology.bonds.len() {
        return Err(AtropisomerError::AxialBondOutOfRange {
            bond: axial_bond,
            bond_count: topology.bonds.len(),
        });
    }
    match conformer {
        AtropisomerConformer::TwoD(value) => value.validate_for_atom_count(topology.atoms.len()),
        AtropisomerConformer::ThreeD(value) => value.validate_for_atom_count(topology.atoms.len()),
    }
    .map_err(|source| AtropisomerError::InvalidCoordinates { source })?;
    let mut state = MutableAtropWedgeState { topology, wedges };
    match wedge_one_two_d_source(&mut state, rings, axial_bond, conformer) {
        Ok(()) => Ok(true),
        Err(PerceptionError::Rejected(_)) => Ok(false),
        Err(PerceptionError::SourcePrecondition { message }) => {
            Err(AtropisomerError::SourcePrecondition { message })
        }
        Err(PerceptionError::Normalization(error)) => Err(AtropisomerError::Normalization(error)),
        Err(PerceptionError::CarrierCount { atom, count }) => {
            Err(AtropisomerError::CarrierCount { atom, count })
        }
    }
}

fn current_three_d_direction<S: AtropWedgeState>(
    state: &S,
    bond: BondId,
    conformer: AtropisomerConformer<'_>,
) -> BondDirection {
    // A borrowed scalar getter transport to the sole native marker owner.
    three_d_direction(
        &CurrentAtropBondEndpoints {
            id: bond,
            begin: state.begin(bond),
            end: state.end(bond),
        },
        conformer,
    )
}

fn wedge_one_three_d_source<S: AtropWedgeState>(
    state: &mut S,
    rings: &RingInfo,
    axial_id: BondId,
    conformer: AtropisomerConformer<'_>,
) -> Result<(), PerceptionError> {
    // RDKit❗❌: bool WedgeBondFromAtropisomerOneBond3d(
    // RDKit❗❌:     Bond *bond, const ROMol &mol, const Conformer *conf,
    // RDKit❗❌:     std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
    // RDKit❗❌:         &wedgeBonds) {
    // RDKit❗❌:   PRECONDITION(bond, "bad bond");
    // RDKit❗❌:
    // RDKit❗❌:   AtropAtomAndBondVec atomAndBondVecs[2];
    // RDKit❗❌:   if (!getAtropisomerAtomsAndBonds(bond, atomAndBondVecs, mol)) {
    // RDKit❗❌:     return false;  // not an atropisomer
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   //  make sure we do not have wiggle bonds
    // RDKit❗❌:
    // RDKit❗❌:   for (auto atomAndBondVecs : atomAndBondVecs) {
    // RDKit❗❌:     for (auto endBond : atomAndBondVecs.second) {
    // RDKit❗❌:       if (endBond->getBondDir() == Bond::UNKNOWN) {
    // RDKit❗❌:         return false;  // not an atropisomer)
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // first see if any candidate bond is already set to a wedge or hash
    // RDKit❗❌:   // if so, we will use that bond as a wedge or hash
    // RDKit❗❌:
    // RDKit❗❌:   std::vector<Bond *> useBonds;
    // RDKit❗❌:
    // RDKit❗❌:   for (unsigned int whichEnd = 0; whichEnd < 2; ++whichEnd) {
    // RDKit❗❌:     for (unsigned int whichBond = 0;
    // RDKit❗❌:          whichBond < atomAndBondVecs[whichEnd].second.size(); ++whichBond) {
    // RDKit❗❌:       auto bond = atomAndBondVecs[whichEnd].second[whichBond];
    // RDKit❗❌:       auto bondDir = bond->getBondDir();
    // RDKit❗❌:
    // RDKit❗❌:       // see if it is a wedge or hash and its origin is the atom in the
    // RDKit❗❌:       // main bond
    // RDKit❗❌:
    // RDKit❗❌:       if ((bondDir == Bond::BEGINWEDGE || bondDir == Bond::BEGINDASH) &&
    // RDKit❗❌:           bond->getBeginAtom() == atomAndBondVecs[whichEnd].first &&
    // RDKit❗❌:           canHaveDirection(*bond)) {
    // RDKit❗❌:         useBonds.push_back(bond);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // the following may seem redundant, since we just found the useBonds
    // RDKit❗❌:   // based on their bond dir PRESENCE, but this endures that the values are
    // RDKit❗❌:   // correct.
    // RDKit❗❌:
    // RDKit❗❌:   if (useBonds.size() > 0) {
    // RDKit❗❌:     for (auto useBond : useBonds) {
    // RDKit❗❌:       useBond->setBondDir(getBondDirForAtropisomer3d(useBond, conf));
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     return true;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // did not find a used bond dir - pick one to use
    // RDKit❗❌:   // we would like to have one that is not in a ring, and will be a dash
    // RDKit❗❌:
    // RDKit❗❌:   const RingInfo *ri = bond->getOwningMol().getRingInfo();
    // RDKit❗❌:
    // RDKit❗❌:   Bond *bestBond = nullptr;
    // RDKit❗❌:   int bestBondEnd = -1;
    // RDKit❗❌:   unsigned int bestRingCount = UINT_MAX;
    // RDKit❗❌:   unsigned int largestRingSize = 0;
    // RDKit❗❌:   Bond::BondDir bestBondDir = Bond::BondDir::NONE;
    // RDKit❗❌:   bool bestBondIsSingle = false;
    // RDKit❗❌:   for (unsigned int whichEnd = 0; whichEnd < 2; ++whichEnd) {
    // RDKit❗❌:     for (unsigned int whichBond = 0;
    // RDKit❗❌:          whichBond < atomAndBondVecs[whichEnd].second.size(); ++whichBond) {
    // RDKit❗❌:       auto bondToTry = atomAndBondVecs[whichEnd].second[whichBond];
    // RDKit❗❌:
    // RDKit❗❌:       // cannot use a bond that is not single, nor if it is already slated
    // RDKit❗❌:       // to be used for a chiral center
    // RDKit❗❌:
    // RDKit❗❌:       if (!canHaveDirection(*bondToTry) ||
    // RDKit❗❌:           wedgeBonds.find(bond->getIdx()) != wedgeBonds.end()) {
    // RDKit❗❌:         continue;  // must be a single bond and not already spoken
    // RDKit❗❌:                    // for by a chiral center
    // RDKit❗❌:       }
    // RDKit❗❌:
    // RDKit❗❌:       // make sure the atoms on the bond are in the right order for the
    // RDKit❗❌:       // wedge/hash the atom on the end of the main bond must be listed
    // RDKit❗❌:       // first
    // RDKit❗❌:
    // RDKit❗❌:       if (bondToTry->getBeginAtom() != atomAndBondVecs[whichEnd].first) {
    // RDKit❗❌:         bondToTry->setEndAtom(bondToTry->getBeginAtom());
    // RDKit❗❌:         bondToTry->setBeginAtom(atomAndBondVecs[whichEnd].first);
    // RDKit❗❌:       }
    // RDKit❗❌:
    // RDKit❗❌:       if (bondToTry->getBondDir() != Bond::BondDir::NONE) {
    // RDKit❗❌:         if (bondToTry->getBeginAtom()->getIdx() ==
    // RDKit❗❌:             atomAndBondVecs[whichEnd].first->getIdx()) {
    // RDKit❗❌:           BOOST_LOG(rdWarningLog)
    // RDKit❗❌:               << "Wedge or hash bond found on atropisomer where not expected - atoms are: "
    // RDKit❗❌:               << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx()
    // RDKit❗❌:               << std::endl;
    // RDKit❗❌:           return false;
    // RDKit❗❌:         } else {
    // RDKit❗❌:           continue;  // wedge or hash bond affecting the OTHER atom
    // RDKit❗❌:                      // = perhaps a chiral center
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       auto ringCount = ri->numBondRings(bondToTry->getIdx());
    // RDKit❗❌:       unsigned int ringSize = 0;
    // RDKit❗❌:       if (!ringCount) {
    // RDKit❗❌:         ringCount = 10;
    // RDKit❗❌:       } else {
    // RDKit❗❌:         // we're going to prefer to put wedges in larger rings, but don't want
    // RDKit❗❌:         // to end up wedging macrocyles if it's avoidable.
    // RDKit❗❌:         ringSize = ri->minBondRingSize(bondToTry->getIdx());
    // RDKit❗❌:         if (ringSize > 8) {
    // RDKit❗❌:           ringSize = 0;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       if (ringCount > bestRingCount) {
    // RDKit❗❌:         continue;
    // RDKit❗❌:       } else if (ringCount < bestRingCount || ringSize > largestRingSize) {
    // RDKit❗❌:         bestBond = bondToTry;
    // RDKit❗❌:         bestBondEnd = whichEnd;
    // RDKit❗❌:         bestRingCount = ringCount;
    // RDKit❗❌:         largestRingSize = ringSize;
    // RDKit❗❌:         bestBondIsSingle = (bondToTry->getBondType() == Bond::BondType::SINGLE);
    // RDKit❗❌:         bestBondDir = getBondDirForAtropisomer3d(bondToTry, conf);
    // RDKit❗❌:       } else if (bestBondIsSingle &&
    // RDKit❗❌:                  bondToTry->getBondType() != Bond::BondType::SINGLE) {
    // RDKit❗❌:         continue;
    // RDKit❗❌:       } else if (!bestBondIsSingle &&
    // RDKit❗❌:                  bondToTry->getBondType() == Bond::BondType::SINGLE) {
    // RDKit❗❌:         bestBondEnd = whichEnd;
    // RDKit❗❌:         bestBond = bondToTry;
    // RDKit❗❌:         bestRingCount = ringCount;
    // RDKit❗❌:         bestBondIsSingle = true;
    // RDKit❗❌:         bestBondDir = getBondDirForAtropisomer3d(bondToTry, conf);
    // RDKit❗❌:       } else {
    // RDKit❗❌:         auto bondDir = getBondDirForAtropisomer3d(bondToTry, conf);
    // RDKit❗❌:         if (bestBondDir == Bond::BondDir::NONE ||
    // RDKit❗❌:             (bestBondDir == Bond::BondDir::BEGINDASH &&
    // RDKit❗❌:              bondDir == Bond::BondDir::BEGINWEDGE)) {
    // RDKit❗❌:           bestBond = bondToTry;
    // RDKit❗❌:           bestBondEnd = whichEnd;
    // RDKit❗❌:           bestRingCount = ringCount;
    // RDKit❗❌:           bestBondIsSingle =
    // RDKit❗❌:               (bondToTry->getBondType() == Bond::BondType::SINGLE);
    // RDKit❗❌:
    // RDKit❗❌:           bestBondDir = bondDir;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (bestBond != nullptr) {
    // RDKit❗❌:     // we found a good one
    // RDKit❗❌:
    // RDKit❗❌:     // make sure the atoms on the bond are in the right order for the
    // RDKit❗❌:     // wedge/hash the atom on the end of the main bond must be listed
    // RDKit❗❌:     // first for the wedge/has bond
    // RDKit❗❌:
    // RDKit❗❌:     if (bestBond->getBeginAtom() != atomAndBondVecs[bestBondEnd].first) {
    // RDKit❗❌:       bestBond->setEndAtom(bestBond->getBeginAtom());
    // RDKit❗❌:       bestBond->setBeginAtom(atomAndBondVecs[bestBondEnd].first);
    // RDKit❗❌:     }
    // RDKit❗❌:     bestBond->setBondDir(bestBondDir);
    // RDKit❗❌:     auto newWedgeInfo = std::unique_ptr<RDKit::Chirality::WedgeInfoBase>(
    // RDKit❗❌:         new RDKit::Chirality::WedgeInfoAtropisomer(bond->getIdx(),
    // RDKit❗❌:                                                    bestBondDir));
    // RDKit❗❌:
    // RDKit❗❌:     wedgeBonds[bestBond->getIdx()] = std::move(newWedgeInfo);
    // RDKit❗❌:   } else {
    // RDKit❗❌:     BOOST_LOG(rdWarningLog)
    // RDKit❗❌:         << "Failed to find a good bond to set as UP or DOWN for an atropisomer - atoms are: "
    // RDKit❗❌:         << bond->getBeginAtomIdx() << " " << bond->getEndAtomIdx() << std::endl;
    // RDKit❗❌:     return false;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   return true;
    // RDKit❗❌: }
    // Behavior review: actual current carrier endpoints/directions/map state,
    // literal source guards and trial-orientation prefix before direction/ring
    // errors. Native occupancy checks the AXIAL key, not each carrier key.
    // No 2D frame/representative-vector or stereo precondition is added.
    // Complexity review: native-sized carrier/use vectors and degree scan.
    // Projected transport adds logarithmic getter/write queries; checked entry
    // validates detached shapes once. No whole graph/property/coordinate copy.
    let current_axial = CurrentAtropBondEndpoints {
        id: axial_id,
        begin: state.begin(axial_id),
        end: state.end(axial_id),
    };
    let ends = atropisomer_ends_fresh(state.topology(), &current_axial)?
        .ok_or(AtropisomerRejectionKind::MissingCarrier)?;
    for end in &ends {
        for &carrier in &end.bonds {
            if state.direction(carrier) == BondDirection::Unknown {
                return Err(AtropisomerRejectionKind::UnknownCarrierDirection.into());
            }
        }
    }
    let mut use_bonds = Vec::new();
    for end in &ends {
        for &carrier in &end.bonds {
            if matches!(
                state.direction(carrier),
                BondDirection::BeginWedge | BondDirection::BeginDash
            ) && state.begin(carrier) == end.atom
                && can_have_direction(&state.topology().bonds()[carrier.index()])
            {
                use_bonds.push(carrier);
            }
        }
    }
    if !use_bonds.is_empty() {
        for carrier in use_bonds {
            let direction = current_three_d_direction(state, carrier, conformer);
            state.set_direction(carrier, direction, axial_id);
        }
        return Ok(());
    }
    let mut best: Option<(usize, BondId)> = None;
    let mut best_ring_count = u32::MAX;
    let mut largest_ring_size = 0;
    let mut best_direction = BondDirection::None;
    let mut best_single = false;
    for (which_end, end) in ends.iter().enumerate() {
        for &carrier in &end.bonds {
            if !can_have_direction(&state.topology().bonds()[carrier.index()])
                || state.occupied(axial_id)
            {
                continue;
            }
            state.orient(carrier, end.atom, axial_id);
            if state.direction(carrier) != BondDirection::None {
                if state.begin(carrier) == end.atom {
                    state.source_warning(
                        "Wedge or hash bond found on atropisomer where not expected - atoms are:",
                        axial_id,
                    );
                    return Err(AtropisomerRejectionKind::DirectionConflict.into());
                }
                continue;
            }
            if !rings.is_initialized() {
                return Err(PerceptionError::SourcePrecondition {
                    message: "RingInfo not initialized",
                });
            }
            let mut ring_count = rings.num_bond_rings(carrier) as u32;
            let mut ring_size = 0;
            if ring_count == 0 {
                ring_count = 10;
            } else {
                ring_size = rings.min_bond_ring_size(carrier) as u32;
                if ring_size > 8 {
                    ring_size = 0;
                }
            }
            let single = state.topology().bonds()[carrier.index()].order() == BondOrder::Single;
            if ring_count > best_ring_count {
                continue;
            } else if ring_count < best_ring_count || ring_size > largest_ring_size {
                best = Some((which_end, carrier));
                best_ring_count = ring_count;
                largest_ring_size = ring_size;
                best_single = single;
                best_direction = current_three_d_direction(state, carrier, conformer);
            } else if best_single && !single {
                continue;
            } else if !best_single && single {
                best = Some((which_end, carrier));
                best_ring_count = ring_count;
                best_single = true;
                best_direction = current_three_d_direction(state, carrier, conformer);
            } else {
                let direction = current_three_d_direction(state, carrier, conformer);
                if best_direction == BondDirection::None
                    || (best_direction == BondDirection::BeginDash
                        && direction == BondDirection::BeginWedge)
                {
                    best = Some((which_end, carrier));
                    best_ring_count = ring_count;
                    best_single = single;
                    best_direction = direction;
                }
            }
        }
    }
    if let Some((which_end, carrier)) = best {
        state.orient(carrier, ends[which_end].atom, axial_id);
        state.set_direction(carrier, best_direction, axial_id);
        state.insert_wedge(carrier, axial_id);
        Ok(())
    } else {
        {
            state.source_warning(
                "Failed to find a good bond to set as UP or DOWN for an atropisomer - atoms are:",
                axial_id,
            );
            Err(AtropisomerRejectionKind::NoUsableWedgeBond.into())
        }
    }
}

/// Applies the native 3D single-atropisomer behavior to actual detached state.
#[doc(hidden)]
pub fn wedge_atropisomer_three_d_source(
    topology: &mut TopologyBlock,
    rings: &RingInfo,
    axial_bond: BondId,
    conformer: AtropisomerConformer<'_>,
    wedges: &mut WedgeAssignments,
) -> Result<bool, AtropisomerError> {
    topology
        .validate()
        .map_err(|source| AtropisomerError::InvalidTopology { source })?;
    if axial_bond.index() >= topology.bonds.len() {
        return Err(AtropisomerError::AxialBondOutOfRange {
            bond: axial_bond,
            bond_count: topology.bonds.len(),
        });
    }
    match conformer {
        AtropisomerConformer::TwoD(value) => value.validate_for_atom_count(topology.atoms.len()),
        AtropisomerConformer::ThreeD(value) => value.validate_for_atom_count(topology.atoms.len()),
    }
    .map_err(|source| AtropisomerError::InvalidCoordinates { source })?;
    let mut state = MutableAtropWedgeState { topology, wedges };
    match wedge_one_three_d_source(&mut state, rings, axial_bond, conformer) {
        Ok(()) => Ok(true),
        Err(PerceptionError::Rejected(_)) => Ok(false),
        Err(PerceptionError::SourcePrecondition { message }) => {
            Err(AtropisomerError::SourcePrecondition { message })
        }
        Err(PerceptionError::Normalization(error)) => Err(AtropisomerError::Normalization(error)),
        Err(PerceptionError::CarrierCount { atom, count }) => {
            Err(AtropisomerError::CarrierCount { atom, count })
        }
    }
}

fn validate_rings(topology: &TopologyBlock, rings: &RingInfo) -> Result<(), AtropisomerError> {
    if !rings.is_sssr_or_better() {
        return Err(AtropisomerError::RingInfoNotSssr);
    }
    if rings.atom_row_count() != topology.atoms.len() {
        return Err(AtropisomerError::RingAtomRowCount {
            actual: rings.atom_row_count(),
            expected: topology.atoms.len(),
        });
    }
    if rings.bond_row_count() != topology.bonds.len() {
        return Err(AtropisomerError::RingBondRowCount {
            actual: rings.bond_row_count(),
            expected: topology.bonds.len(),
        });
    }
    for ring in rings.atom_rings() {
        for atom in ring {
            if atom.index() >= topology.atoms.len() {
                return Err(AtropisomerError::RingAtomOutOfRange {
                    atom: *atom,
                    atom_count: topology.atoms.len(),
                });
            }
        }
    }
    for ring in rings.bond_rings() {
        for bond in ring {
            if bond.index() >= topology.bonds.len() {
                return Err(AtropisomerError::RingBondOutOfRange {
                    bond: *bond,
                    bond_count: topology.bonds.len(),
                });
            }
        }
    }
    Ok(())
}

fn check_invalid_atrop_bond(
    bond: &mut Bond,
    hybridization: &HybridizationAssignment,
    rings: &RingInfo,
) -> Result<bool, AtropisomerError> {
    // Complete pinned source: MolOps.cpp anonymous-namespace checkBond.
    // RDKit✔️✔️: bool checkBond(RWMol &mol, Bond *bond, MolOps::Hybridizations &hybs) {
    // RDKit✔️✔️:   if (!mol.getRingInfo()->isSssrOrBetter()) {
    // RDKit✔️✔️:     RDKit::MolOps::findSSSR(mol);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   const RingInfo *ri = mol.getRingInfo();
    // RDKit✔️✔️:   if (hybs[bond->getBeginAtomIdx()] != Atom::SP2 ||
    // RDKit✔️✔️:       hybs[bond->getEndAtomIdx()] != Atom::SP2 ||
    // RDKit✔️✔️:       // do not clear bonds that part of a macrocycle
    // RDKit✔️✔️:       // because they can be linking actual atropisomeric portions
    // RDKit✔️✔️:       (ri->numBondRings(bond->getIdx()) > 0 &&
    // RDKit✔️✔️:        ri->minBondRingSize(bond->getIdx()) < 8)) {
    // RDKit✔️✔️:     bond->setStereo(Bond::BondStereo::STEREONONE);
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // Ring discovery is performed by the canonical ring owner before this
    // detached primitive is called; `validate_rings()` enforces the same
    // SSSR-or-better postcondition before any bond can be changed.
    let invalid = hybridization.values[bond.begin().index()] != Hybridization::Sp2
        || hybridization.values[bond.end().index()] != Hybridization::Sp2
        || (rings.num_bond_rings(bond.id()) > 0 && rings.min_bond_ring_size(bond.id()) < 8);
    if invalid {
        let bond_id = bond.id();
        bond.set_stereo(BondStereo::None)
            .map_err(|source| AtropisomerError::BondUpdate {
                bond: bond_id,
                source,
            })?;
        return Ok(true);
    }
    Ok(false)
}

/// Clear source-invalid atropisomer bond stereo over detached topology values.
///
/// The caller supplies canonical hybridization and ring assignments from the
/// immediately preceding sanitize stages. The borrowed topology is never
/// modified; failure returns no partial result.
pub fn cleanup_invalid_atropisomers(
    topology: &TopologyBlock,
    hybridization: &HybridizationAssignment,
    rings: &RingInfo,
) -> Result<TopologyBlock, AtropisomerError> {
    topology
        .validate()
        .map_err(|source| AtropisomerError::InvalidTopology { source })?;
    validate_rings(topology, rings)?;
    if hybridization.values.len() != topology.atoms.len() {
        return Err(AtropisomerError::HybridizationAssignmentLength {
            actual: hybridization.values.len(),
            expected: topology.atoms.len(),
        });
    }

    // The one private engine below owns the clone/tag-loop/group-cleanup/
    // final-validation body; this wrapper keeps the source wrapper-before-
    // loop prevalidation and supplies the existing detached check primitive.
    cleanup_invalid_atropisomers_engine(topology, |bond| {
        check_invalid_atrop_bond(bond, hybridization, rings)
    })
}

/// One private error-generic cleanup engine shared by the public detached
/// wrapper and the owned-ring-state sanitize wrapper.
///
/// `check` is called ONLY for AtropCw/AtropCcw bonds; it reports whether it
/// cleared the tag. Stereo-group cleanup inspects the RESULT topology after
/// tag updates and runs only if a tag was cleared, exactly as the source
/// does after its bond loop.
fn cleanup_invalid_atropisomers_engine<E>(
    topology: &TopologyBlock,
    mut check: impl FnMut(&mut Bond) -> Result<bool, E>,
) -> Result<TopologyBlock, E>
where
    E: From<AtropisomerError>,
{
    // Complete pinned source: MolOps::cleanupAtropisomers(RWMol &, Hybridizations &).
    // RDKit✔️❌: void cleanupAtropisomers(RWMol &mol, MolOps::Hybridizations &hybs) {
    // RDKit✔️❌:   // make sure that ring info is available
    // RDKit✔️❌:   // (defensive, current calls have it available)
    // RDKit✔️❌:   bool needCleanupAtropisomerStereoGroups = false;
    // RDKit✔️❌:   for (auto bond : mol.bonds()) {
    // RDKit✔️❌:     switch (bond->getStereo()) {
    // RDKit✔️❌:       case Bond::BondStereo::STEREOATROPCW:
    // RDKit✔️❌:       case Bond::BondStereo::STEREOATROPCCW:
    // RDKit✔️❌:         if (checkBond(mol, bond, hybs)) {
    // RDKit✔️❌:           needCleanupAtropisomerStereoGroups = true;
    // RDKit✔️❌:         }
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       default:
    // RDKit✔️❌:         break;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   if (needCleanupAtropisomerStereoGroups) {
    // RDKit✔️❌:     Atropisomers::cleanupAtropisomerStereoGroups(mol);
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // Returning a detached owned value requires one O(atoms+bonds) clone that
    // the source in-place operation does not perform; loop and lookup costs
    // after that clone remain linear/direct-indexed.
    let mut result = topology.clone();
    let mut need_stereo_group_cleanup = false;
    for bond in &mut result.bonds {
        if matches!(bond.stereo(), BondStereo::AtropCw | BondStereo::AtropCcw) && check(bond)? {
            need_stereo_group_cleanup = true;
        }
    }
    if need_stereo_group_cleanup {
        #[cfg(test)]
        cleanup_ring_state_probe::record_group_cleanup();
        result.stereo_groups =
            cleanup_atropisomer_stereo_groups(&result, &AtropisomerAssignment::default())
                .map_err(E::from)?
                .groups;
    }
    result
        .validate()
        .map_err(|source| E::from(AtropisomerError::InvalidTopology { source }))?;
    Ok(result)
}

/// Private two-cause error for the owned-ring-state cleanup wrapper.
///
/// `Rings` carries the canonical finder failure raised while lazily
/// acquiring SSSR for a tagged bond; `Algorithm` carries the existing
/// detached atropisomer causes. `source()` borrows the inner error.
#[derive(Debug, thiserror::Error)]
pub(crate) enum AtropisomerCleanupError {
    #[error("ring finding failed during atropisomer cleanup: {0}")]
    Rings(#[from] crate::RingFindingError),
    #[error("atropisomer cleanup failed: {0}")]
    Algorithm(#[from] AtropisomerError),
}

/// Test-only observation points at the actual wrapper acquisition site.
/// Counters are never reset by production code; production builds contain
/// no counters or hooks.
#[cfg(test)]
pub(crate) mod cleanup_ring_state_probe {
    use std::cell::Cell;

    thread_local! {
        static FIND_CALLS: Cell<u64> = const { Cell::new(0) };
        static CONSUMED_CALLS: Cell<u64> = const { Cell::new(0) };
        static GROUP_CLEANUP_CALLS: Cell<u64> = const { Cell::new(0) };
        static ACQUIRED_ATOM_ROW_PTRS: Cell<Vec<usize>> = const { Cell::new(Vec::new()) };
    }

    pub(crate) fn record_find() {
        FIND_CALLS.with(|calls| calls.set(calls.get() + 1));
    }

    pub(crate) fn record_consumed() {
        CONSUMED_CALLS.with(|calls| calls.set(calls.get() + 1));
    }

    pub(crate) fn record_group_cleanup() {
        GROUP_CLEANUP_CALLS.with(|calls| calls.set(calls.get() + 1));
    }

    pub(crate) fn record_acquired_atom_rows(rows: *const Vec<cosmolkit_model::AtomId>) {
        ACQUIRED_ATOM_ROW_PTRS.with(|buffer| {
            let mut buffer = buffer.take();
            buffer.push(rows as usize);
            ACQUIRED_ATOM_ROW_PTRS.with(|cell| cell.set(buffer));
        });
    }

    pub(crate) fn find_calls() -> u64 {
        FIND_CALLS.with(Cell::get)
    }

    pub(crate) fn consumed_calls() -> u64 {
        CONSUMED_CALLS.with(Cell::get)
    }

    pub(crate) fn group_cleanup_calls() -> u64 {
        GROUP_CLEANUP_CALLS.with(Cell::get)
    }

    pub(crate) fn acquired_atom_row_pointers() -> Vec<usize> {
        ACQUIRED_ATOM_ROW_PTRS.with(Cell::take)
    }
}

/// Owned-ring-state cleanup used by the sanitize CLEANUP_ATROPISOMERS stage.
///
/// The supplied `Option<RingInfo>` is the one state carrier: it is moved in
/// and moved back out; no full RingInfo clone serves this path. With no
/// AtropCw/AtropCcw bonds the callback never runs: NO ring acquisition or
/// validation happens and the supplied state is preserved exactly (a
/// malformed no-tag state remains unconsumed). For a tagged bond the source
/// checkBond guard runs first: a state below SSSR-or-better acquires ONE
/// canonical find_sssr from the ORIGINAL topology (cleanup changes only
/// stereo tags in its clone, never graph/bond orders/endpoints/IDs) and
/// replaces the carrier by move BEFORE any Sp2 endpoint short-circuit; an
/// initialized-empty SSSR/Symm satisfies the guard and is never re-found.
/// The consumed state is then dimension/index validated and the existing
/// detached check primitive decides the tag.
pub(crate) fn cleanup_invalid_atropisomers_with_ring_state(
    topology: &TopologyBlock,
    hybridization: &HybridizationAssignment,
    rings: Option<RingInfo>,
) -> Result<(TopologyBlock, Option<RingInfo>), AtropisomerCleanupError> {
    topology
        .validate()
        .map_err(|source| AtropisomerError::InvalidTopology { source })?;
    if hybridization.values.len() != topology.atoms.len() {
        return Err(AtropisomerError::HybridizationAssignmentLength {
            actual: hybridization.values.len(),
            expected: topology.atoms.len(),
        }
        .into());
    }
    let mut state = rings;
    let output =
        cleanup_invalid_atropisomers_engine::<AtropisomerCleanupError>(topology, |bond| {
            #[cfg(test)]
            cleanup_ring_state_probe::record_consumed();
            // Complete pinned source: MolOps.cpp anonymous-namespace checkBond
            // ring-state guard, applied to the owned carrier.
            // RDKit✔️✔️:   if (!mol.getRingInfo()->isSssrOrBetter()) {
            // RDKit✔️✔️:     RDKit::MolOps::findSSSR(mol);
            // RDKit✔️✔️:   }
            let below_quality = !state
                .as_ref()
                .is_some_and(|carrier| carrier.is_sssr_or_better());
            if below_quality {
                let fresh = crate::rings::find_sssr(topology, &crate::RingSearchParams::default())?;
                #[cfg(test)]
                cleanup_ring_state_probe::record_find();
                #[cfg(test)]
                cleanup_ring_state_probe::record_acquired_atom_rows(fresh.atom_rings().as_ptr());
                state = Some(fresh);
            }
            let carrier = state
                .as_ref()
                .expect("ring state present after the guard branch");
            validate_rings(topology, carrier)?;
            check_invalid_atrop_bond(bond, hybridization, carrier)
                .map_err(AtropisomerCleanupError::from)
        })?;
    Ok((output, state))
}

// Native source cache can be mutable; the existing checked query supplies
// already-validated shared SSSR state and never copies unchanged ring rows.
enum AtropSourceRings<'a> {
    Mutable {
        rings: &'a mut RingInfo,
        properties: Option<&'a mut cosmolkit_model::MoleculeProperties>,
    },
    Validated(&'a RingInfo),
}
impl AtropSourceRings<'_> {
    fn current(&self) -> &RingInfo {
        match self {
            Self::Mutable { rings, .. } => rings,
            Self::Validated(rings) => rings,
        }
    }
    fn acquire<G: StereoGraphAccess>(&mut self, topology: &G) -> Result<(), AtropisomerError> {
        match self {
            Self::Mutable { rings, properties } => {
                crate::rings::find_sssr_with_source_outputs_from_graph(
                    topology,
                    rings,
                    properties.as_deref_mut(),
                    None,
                    false,
                    false,
                )?;
                Ok(())
            }
            // Only the private checked adapter constructs this case, after
            // validate_rings; the native cache-acquisition branch is unreachable.
            Self::Validated(_) => Err(AtropisomerError::RingInfoNotSssr),
        }
    }
}

fn wedge_bonds_from_atropisomers_kernel<S: AtropWedgeState>(
    state: &mut S,
    rings: &mut AtropSourceRings<'_>,
    conformer: Option<AtropisomerConformer<'_>>,
    diagnostics: &mut Vec<AtropisomerDiagnostic>,
) -> Result<(), AtropisomerError> {
    // RDKit❗❌: void wedgeBondsFromAtropisomers(
    // RDKit❗❌:     const ROMol &mol, const Conformer *conf,
    // RDKit❗❌:     std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
    // RDKit❗❌:         &wedgeBonds) {
    // RDKit❗❌:   PRECONDITION(conf == nullptr || &(conf->getOwningMol()) == &mol,
    // RDKit❗❌:                "conformer does not belong to molecule");
    // RDKit❗❌:
    // RDKit❗❌:   // WedgeBondFromAtropisomerOneBond 2d/3d requires ring bond counts
    // RDKit❗❌:   if (!mol.getRingInfo()->isSssrOrBetter()) {
    // RDKit❗❌:     RDKit::MolOps::findSSSR(mol);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   for (auto bond : mol.bonds()) {
    // RDKit❗❌:     auto bondStereo = bond->getStereo();
    // RDKit❗❌:
    // RDKit❗❌:     if (bond->getBondType() != Bond::BondType::SINGLE ||
    // RDKit❗❌:         (bondStereo != Bond::BondStereo::STEREOATROPCW &&
    // RDKit❗❌:          bondStereo != Bond::BondStereo::STEREOATROPCCW) ||
    // RDKit❗❌:         bond->getBeginAtom()->getTotalDegree() < 2 ||
    // RDKit❗❌:         bond->getEndAtom()->getTotalDegree() < 2 ||
    // RDKit❗❌:         bond->getBeginAtom()->getTotalDegree() > 3 ||
    // RDKit❗❌:         bond->getEndAtom()->getTotalDegree() > 3) {
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     if (conf) {
    // RDKit❗❌:       if (conf->is3D()) {
    // RDKit❗❌:         WedgeBondFromAtropisomerOneBond3d(bond, mol, conf, wedgeBonds);
    // RDKit❗❌:       } else {
    // RDKit❗❌:         WedgeBondFromAtropisomerOneBond2d(bond, mol, conf, wedgeBonds);
    // RDKit❗❌:       }
    // RDKit❗❌:     } else {  // no conformer
    // RDKit❗❌:       WedgeBondFromAtropisomerOneBondNoConf(bond, mol, wedgeBonds);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // Behavior review: promote source ring state before bond iteration, trust
    // legal initialized sparse source rows, and read actual cached H counts.
    // Native condition/getter order and ignored single-bond false results are
    // preserved. Current endpoints include earlier source trial orientation.
    // Complexity review: unique single-bond kernels, indexed graph iteration,
    // source ordered wedge map. Detached actual implicit cache projection is
    // O(V) extra storage; no whole graph, properties or ring-state clone.
    if !rings.current().is_sssr_or_better() {
        rings.acquire(state.topology())?;
    }
    let valence = crate::ValenceAssignment {
        explicit_valence: Vec::new(),
        implicit_hydrogens: state
            .topology()
            .atoms()
            .iter()
            .map(|atom| i32::from(atom.source_valence_facts().implicit_valence))
            .collect(),
    };
    for index in 0..state.topology().bonds().len() {
        let axial_id = BondId::new(index);
        let axial = &state.topology().bonds()[index];
        if axial.order() != BondOrder::Single
            || !matches!(axial.stereo(), BondStereo::AtropCw | BondStereo::AtropCcw)
            || source_total_degree(state.topology(), &valence, state.begin(axial_id))? < 2
            || source_total_degree(state.topology(), &valence, state.end(axial_id))? < 2
            || source_total_degree(state.topology(), &valence, state.begin(axial_id))? > 3
            || source_total_degree(state.topology(), &valence, state.end(axial_id))? > 3
        {
            continue;
        }
        let result = match conformer {
            Some(value) if conformer_is_3d(value) => {
                wedge_one_three_d_source(state, rings.current(), axial_id, value)
            }
            Some(value) => wedge_one_two_d_source(state, rings.current(), axial_id, value),
            None => wedge_one_no_conformer_source(state, rings.current(), axial_id),
        };
        match result {
            Ok(()) => {}
            Err(PerceptionError::Rejected(kind)) => diagnostics.push(AtropisomerDiagnostic {
                bond: axial_id,
                kind,
            }),
            Err(PerceptionError::SourcePrecondition { message }) => {
                return Err(AtropisomerError::SourcePrecondition { message });
            }
            Err(PerceptionError::Normalization(error)) => {
                return Err(AtropisomerError::Normalization(error));
            }
            Err(PerceptionError::CarrierCount { atom, count }) => {
                return Err(AtropisomerError::CarrierCount { atom, count });
            }
        }
    }
    Ok(())
}

/// Runs native cache/mutation ordering on actual detached graph and wedge map.
#[doc(hidden)]
pub fn wedge_bonds_from_atropisomers_source(
    topology: &mut TopologyBlock,
    rings: &mut RingInfo,
    properties: Option<&mut cosmolkit_model::MoleculeProperties>,
    conformer: Option<AtropisomerConformer<'_>>,
    wedges: &mut WedgeAssignments,
) -> Result<(), AtropisomerError> {
    topology
        .validate()
        .map_err(|source| AtropisomerError::InvalidTopology { source })?;
    if let Some(value) = conformer {
        validate_conformer(value, topology.atoms.len())?;
    }
    wedge_bonds_from_atropisomers_graph_source(topology, rings, properties, conformer, wedges)
}

pub(crate) fn wedge_bonds_from_atropisomers_graph_source<G: StereoGraphMut>(
    topology: &mut G,
    rings: &mut RingInfo,
    properties: Option<&mut cosmolkit_model::MoleculeProperties>,
    conformer: Option<AtropisomerConformer<'_>>,
    wedges: &mut WedgeAssignments,
) -> Result<(), AtropisomerError> {
    let mut source_rings = AtropSourceRings::Mutable { rings, properties };
    let mut diagnostics = Vec::new();
    let mut state = MutableAtropWedgeState { topology, wedges };
    let result = wedge_bonds_from_atropisomers_kernel(
        &mut state,
        &mut source_rings,
        conformer,
        &mut diagnostics,
    );
    state
        .wedges
        .append_source_atropisomer_diagnostics(diagnostics);
    result
}

fn wedge_atropisomer_projection<G: StereoGraphAccess>(
    topology: &G,
    rings: &mut AtropSourceRings<'_>,
    conformer: Option<AtropisomerConformer<'_>>,
    occupied_bonds: &BTreeSet<BondId>,
) -> Result<AtropisomerWedgeAssignment, AtropisomerError> {
    let mut updates = BTreeMap::new();
    let mut map_writes = BTreeSet::new();
    let mut diagnostics = Vec::new();
    let mut state = ProjectedAtropWedgeState {
        topology,
        occupied: occupied_bonds,
        updates: &mut updates,
        map_writes: &mut map_writes,
    };
    wedge_bonds_from_atropisomers_kernel(&mut state, rings, conformer, &mut diagnostics)?;
    Ok(AtropisomerWedgeAssignment {
        source_map_writes: map_writes.into_iter().collect(),
        bond_updates: updates.into_values().collect(),
        diagnostics,
    })
}

/// Projects native source graph writes while retaining actual mutable cache.
#[doc(hidden)]
pub fn wedge_bonds_from_atropisomers_projected_source(
    topology: &TopologyBlock,
    rings: &mut RingInfo,
    properties: Option<&mut cosmolkit_model::MoleculeProperties>,
    conformer: Option<AtropisomerConformer<'_>>,
    occupied_bonds: &BTreeSet<BondId>,
) -> Result<AtropisomerWedgeAssignment, AtropisomerError> {
    topology
        .validate()
        .map_err(|source| AtropisomerError::InvalidTopology { source })?;
    if let Some(value) = conformer {
        validate_conformer(value, topology.atoms.len())?;
    }
    for &bond in occupied_bonds {
        if bond.index() >= topology.bonds.len() {
            return Err(AtropisomerError::AssignmentBondOutOfRange {
                bond,
                bond_count: topology.bonds.len(),
            });
        }
    }
    wedge_bonds_from_atropisomers_projected_graph_source(
        topology,
        rings,
        properties,
        conformer,
        occupied_bonds,
    )
}

pub(crate) fn wedge_bonds_from_atropisomers_projected_graph_source<G: StereoGraphAccess>(
    topology: &G,
    rings: &mut RingInfo,
    properties: Option<&mut cosmolkit_model::MoleculeProperties>,
    conformer: Option<AtropisomerConformer<'_>>,
    occupied_bonds: &BTreeSet<BondId>,
) -> Result<AtropisomerWedgeAssignment, AtropisomerError> {
    let mut source_rings = AtropSourceRings::Mutable { rings, properties };
    wedge_atropisomer_projection(topology, &mut source_rings, conformer, occupied_bonds)
}

pub fn wedge_bonds_from_atropisomers(
    topology: &TopologyBlock,
    rings: &RingInfo,
    conformer: Option<AtropisomerConformer<'_>>,
    occupied_bonds: &BTreeSet<BondId>,
) -> Result<AtropisomerWedgeAssignment, AtropisomerError> {
    topology
        .validate()
        .map_err(|source| AtropisomerError::InvalidTopology { source })?;
    validate_rings(topology, rings)?;
    if let Some(value) = conformer {
        validate_conformer(value, topology.atoms.len())?;
    }
    for &bond in occupied_bonds {
        if bond.index() >= topology.bonds.len() {
            return Err(AtropisomerError::AssignmentBondOutOfRange {
                bond,
                bond_count: topology.bonds.len(),
            });
        }
    }
    let mut source_rings = AtropSourceRings::Validated(rings);
    wedge_atropisomer_projection(topology, &mut source_rings, conformer, occupied_bonds)
}

// L05: the owned-ring-state cleanup wrapper products. Counters are
// thread-local and never reset; all observations are per-call deltas.
// Probe pointers prove MOVE preservation of the one carrier (supplied
// buffers survive with the same heap address; the acquired SSSR address
// recorded at the real acquisition site is the address returned).
#[cfg(test)]
mod cleanup_ring_state_tests {
    use super::cleanup_invalid_atropisomers_with_ring_state;
    use super::cleanup_ring_state_probe as probe;
    use crate::atropisomer::AtropisomerError;
    use crate::{HybridizationAssignment, RingFindType, RingInfo};
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, StereoGroup, StereoGroupKind, TopologyBlock,
    };
    use cosmolkit_types::{BondOrder, BondStereo, Element, Hybridization};

    fn cycle(n: usize, stereo: BondStereo) -> TopologyBlock {
        let atoms = (0..n)
            .map(|id| Atom::from_spec(AtomId::new(id), AtomSpec::new(Element::C)))
            .collect::<Vec<_>>();
        let bonds = (0..n)
            .map(|id| {
                let mut spec = BondSpec::new(
                    AtomId::new(id),
                    AtomId::new((id + 1) % n),
                    BondOrder::Single,
                );
                if id == 0 {
                    spec = spec.with_stereo(stereo);
                }
                Bond::from_spec(BondId::new(id), spec)
            })
            .collect::<Vec<_>>();
        // Source-shaped sentinel: one ordered stereo group over all cycle
        // atoms; complete ordered group output is asserted every call.
        let sentinel = StereoGroup::new(
            StereoGroupKind::Or,
            (0..n).map(AtomId::new).collect(),
            Vec::new(),
        )
        .expect("valid distinct stereo members")
        .with_id(7);
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), vec![sentinel]).unwrap()
    }

    fn hybs(values: Vec<Hybridization>) -> HybridizationAssignment {
        HybridizationAssignment { values }
    }

    fn all(n: usize, value: Hybridization) -> Vec<Hybridization> {
        vec![value; n]
    }

    fn full_ring(find_type: RingFindType, n: usize) -> RingInfo {
        // Frozen full-row convention: reversed atom/bond order, e.g. the
        // six-cycle rows atoms[0,5,4,3,2,1], bonds[5,4,3,2,1,0].
        let mut info = RingInfo::new(find_type, n, n);
        let rows: Vec<usize> = (0..n).rev().collect();
        info.add_ring(&rows, &rows).unwrap();
        info
    }

    fn reset_state() -> RingInfo {
        let mut info = RingInfo::new(RingFindType::Sssr, 6, 6);
        info.reset();
        info
    }

    fn row_indices(rings: &RingInfo) -> Vec<Vec<usize>> {
        rings
            .atom_rings()
            .iter()
            .map(|row| row.iter().map(|atom| atom.index()).collect())
            .collect()
    }

    fn bond_row_indices(rings: &RingInfo) -> Vec<Vec<usize>> {
        rings
            .bond_rings()
            .iter()
            .map(|row| row.iter().map(|bond| bond.index()).collect())
            .collect()
    }

    #[test]
    fn sanitize_ring_l05_lazy_quality_sixty_call_product() {
        let states: Vec<(&'static str, Option<RingInfo>)> = vec![
            ("none", None),
            ("reset", Some(reset_state())),
            (
                "other-empty",
                Some(RingInfo::new(RingFindType::OtherOrUnknown, 6, 6)),
            ),
            ("fast-empty", Some(RingInfo::new(RingFindType::Fast, 6, 6))),
            ("sssr-empty", Some(RingInfo::new(RingFindType::Sssr, 6, 6))),
            (
                "symm-empty",
                Some(RingInfo::new(RingFindType::SymmSssr, 6, 6)),
            ),
            (
                "other-full6",
                Some(full_ring(RingFindType::OtherOrUnknown, 6)),
            ),
            ("fast-full6", Some(full_ring(RingFindType::Fast, 6))),
            ("sssr-full6", Some(full_ring(RingFindType::Sssr, 6))),
            ("symm-full6", Some(full_ring(RingFindType::SymmSssr, 6))),
        ];
        let mut calls = 0usize;
        for (state_name, state) in states {
            for stereo in [BondStereo::None, BondStereo::AtropCw, BondStereo::AtropCcw] {
                for endpoint in [Hybridization::Sp2, Hybridization::Sp3] {
                    let label = format!("{state_name}/{stereo:?}/{endpoint:?}");
                    let topology = cycle(6, stereo);
                    let groups_snapshot = topology.stereo_groups.clone();
                    let topology_snapshot = topology.clone();
                    let hybridization = hybs(all(6, endpoint));
                    let hybs_snapshot = hybridization.clone();
                    let state_snapshot = state.clone();
                    // The pointer must be captured from the value actually
                    // PASSED (the moved-in carrier), not the loop-local
                    // original, for a real move-preservation proof.
                    let supplied = state.clone();
                    let supplied_ptr = supplied.as_ref().map(|rings| rings.atom_rings().as_ptr());
                    let find_before = probe::find_calls();
                    let consumed_before = probe::consumed_calls();
                    let group_before = probe::group_cleanup_calls();
                    let acquired_before = probe::acquired_atom_row_pointers().len();
                    let (output, returned) = cleanup_invalid_atropisomers_with_ring_state(
                        &topology,
                        &hybridization,
                        supplied,
                    )
                    .unwrap_or_else(|error| panic!("{label}: unexpected error {error:?}"));
                    calls += 1;
                    let find_delta = probe::find_calls() - find_before;
                    let consumed_delta = probe::consumed_calls() - consumed_before;
                    let group_delta = probe::group_cleanup_calls() - group_before;
                    let acquired = probe::acquired_atom_row_pointers();
                    assert_eq!(topology, topology_snapshot, "{label}: input mutated");
                    assert_eq!(hybridization, hybs_snapshot, "{label}: hybs mutated");
                    let tagged = !matches!(stereo, BondStereo::None);
                    let quality = state_snapshot
                        .as_ref()
                        .is_some_and(|rings| rings.is_sssr_or_better());
                    if !tagged {
                        // No tags: NO acquisition, NO validation, NO group
                        // cleanup; the supplied state is preserved exactly.
                        assert_eq!(consumed_delta, 0, "{label}: consumed");
                        assert_eq!(find_delta, 0, "{label}: find");
                        assert_eq!(group_delta, 0, "{label}: group cleanup");
                        assert_eq!(returned, state_snapshot, "{label}: state changed");
                        assert_eq!(output.bonds[0].stereo(), BondStereo::None, "{label}");
                    } else {
                        assert!(consumed_delta >= 1, "{label}: tagged bond not consumed");
                        // Disposition on the six-cycle: Sp3 always clears;
                        // Sp2 clears whenever b0 sits in a size-6 ring
                        // (acquired full6 or supplied full6) and retains only
                        // for supplied initialized-EMPTY SSSR/Symm.
                        let supplied_empty = quality
                            && state_snapshot
                                .as_ref()
                                .is_some_and(|r| r.atom_rings().is_empty());
                        let retains = endpoint == Hybridization::Sp2 && supplied_empty;
                        if quality {
                            assert_eq!(find_delta, 0, "{label}: re-found supplied quality state");
                            assert_eq!(returned, state_snapshot, "{label}: borrowed state changed");
                            assert_eq!(
                                returned.as_ref().map(|r| r.atom_rings().as_ptr()),
                                supplied_ptr,
                                "{label}: supplied buffer not move-preserved"
                            );
                        } else {
                            // Below SSSR: ONE canonical full6 SSSR acquisition
                            // replaces the carrier by move BEFORE the Sp2
                            // short-circuit, whatever the endpoint value.
                            assert_eq!(find_delta, 1, "{label}: acquisition count");
                            let acquired_state = returned
                                .as_ref()
                                .unwrap_or_else(|| panic!("{label}: acquired state missing"));
                            assert_eq!(acquired_state.find_type(), RingFindType::Sssr, "{label}");
                            assert_eq!(
                                row_indices(acquired_state),
                                vec![vec![0, 5, 4, 3, 2, 1]],
                                "{label}"
                            );
                            assert_eq!(
                                bond_row_indices(acquired_state),
                                vec![vec![5, 4, 3, 2, 1, 0]],
                                "{label}"
                            );
                            assert_eq!(acquired.len(), acquired_before + 1, "{label}");
                            assert_eq!(
                                acquired_state.atom_rings().as_ptr() as usize,
                                acquired[acquired_before],
                                "{label}: returned buffer differs from real acquisition site"
                            );
                        }
                        if retains {
                            assert_eq!(output.bonds[0].stereo(), stereo, "{label}");
                            assert_eq!(group_delta, 0, "{label}: group cleanup ran");
                        } else {
                            assert_eq!(output.bonds[0].stereo(), BondStereo::None, "{label}");
                            // Cleared tags alone trigger group cleanup; with
                            // no remaining atrop bonds the ordered group
                            // output is unchanged.
                            assert_eq!(group_delta, 1, "{label}: group cleanup missing");
                            assert_eq!(output.stereo_groups, groups_snapshot, "{label}: groups");
                        }
                    }
                }
            }
        }
        assert_eq!(calls, 60, "exact census");
    }

    #[test]
    fn sanitize_ring_l05_macrocycle_endpoint_forty_eight_call_product() {
        let mk_state = |name: &str, n: usize| -> Option<RingInfo> {
            match name {
                "none" => None,
                "fast-full" => Some(full_ring(RingFindType::Fast, n)),
                "sssr-full" => Some(full_ring(RingFindType::Sssr, n)),
                _ => Some(full_ring(RingFindType::SymmSssr, n)),
            }
        };
        let mut calls = 0usize;
        for (state_name, quality) in [
            ("none", false),
            ("fast-full", false),
            ("sssr-full", true),
            ("symm-full", true),
        ] {
            for stereo in [BondStereo::AtropCw, BondStereo::AtropCcw] {
                for n in [6usize, 8usize] {
                    for endpoint in [0usize, 1usize, 2usize] {
                        let label = format!("{state_name}/{stereo:?}/c{n}/ep{endpoint}");
                        let mut values = all(n, Hybridization::Sp2);
                        if endpoint == 1 {
                            values[0] = Hybridization::Sp3;
                        }
                        if endpoint == 2 {
                            values[1] = Hybridization::Sp3;
                        }
                        let topology = cycle(n, stereo);
                        let groups_snapshot = topology.stereo_groups.clone();
                        let topology_snapshot = topology.clone();
                        let hybridization = hybs(values);
                        let hybs_snapshot = hybridization.clone();
                        let state = mk_state(state_name, n);
                        let state_snapshot = state.clone();
                        let supplied = state.clone();
                        let supplied_ptr = supplied.as_ref().map(|r| r.atom_rings().as_ptr());
                        let find_before = probe::find_calls();
                        let group_before = probe::group_cleanup_calls();
                        let acquired_before = probe::acquired_atom_row_pointers().len();
                        let (output, returned) = cleanup_invalid_atropisomers_with_ring_state(
                            &topology,
                            &hybridization,
                            supplied,
                        )
                        .unwrap_or_else(|error| panic!("{label}: unexpected error {error:?}"));
                        calls += 1;
                        let find_delta = probe::find_calls() - find_before;
                        let group_delta = probe::group_cleanup_calls() - group_before;
                        let acquired = probe::acquired_atom_row_pointers();
                        assert_eq!(topology, topology_snapshot, "{label}: input mutated");
                        assert_eq!(hybridization, hybs_snapshot, "{label}: hybs mutated");
                        if quality {
                            assert_eq!(find_delta, 0, "{label}: re-found");
                            assert_eq!(returned, state_snapshot, "{label}: state changed");
                            assert_eq!(
                                returned.as_ref().map(|r| r.atom_rings().as_ptr()),
                                supplied_ptr,
                                "{label}: buffer not move-preserved"
                            );
                        } else {
                            assert_eq!(find_delta, 1, "{label}: acquisition count");
                            let acquired_state = returned.as_ref().unwrap();
                            assert_eq!(acquired_state.find_type(), RingFindType::Sssr, "{label}");
                            // Frozen row convention: atoms start at 0 and
                            // run backwards to 1; bonds run n-1..0.
                            let expected_atoms: Vec<usize> =
                                std::iter::once(0).chain((1..n).rev()).collect();
                            let expected_bonds: Vec<usize> = (0..n).rev().collect();
                            assert_eq!(
                                row_indices(acquired_state),
                                vec![expected_atoms],
                                "{label}"
                            );
                            assert_eq!(
                                bond_row_indices(acquired_state),
                                vec![expected_bonds],
                                "{label}"
                            );
                            assert_eq!(acquired.len(), acquired_before + 1, "{label}");
                            assert_eq!(
                                acquired_state.atom_rings().as_ptr() as usize,
                                acquired[acquired_before],
                                "{label}: acquired buffer mismatch"
                            );
                        }
                        // Six-cycle always clears; eight-cycle retains only
                        // with BOTH endpoints Sp2.
                        let retains = n == 8 && endpoint == 0;
                        if retains {
                            assert_eq!(output.bonds[0].stereo(), stereo, "{label}");
                            assert_eq!(group_delta, 0, "{label}: group cleanup ran");
                        } else {
                            assert_eq!(output.bonds[0].stereo(), BondStereo::None, "{label}");
                            assert_eq!(group_delta, 1, "{label}: group cleanup missing");
                            assert_eq!(output.stereo_groups, groups_snapshot, "{label}: groups");
                        }
                    }
                }
            }
        }
        assert_eq!(calls, 48, "exact census");
    }

    #[test]
    fn sanitize_ring_l05_malformed_inputs_and_vocabulary_controls() {
        // Bad hybridization length is rejected before the engine runs.
        let tagged = cycle(6, BondStereo::AtropCw);
        let short = hybs(all(5, Hybridization::Sp2));
        assert!(matches!(
            cleanup_invalid_atropisomers_with_ring_state(&tagged, &short, None),
            Err(super::AtropisomerCleanupError::Algorithm(
                AtropisomerError::HybridizationAssignmentLength {
                    actual: 5,
                    expected: 6
                }
            ))
        ));
        // Invalid topology is rejected first in the existing order.
        let invalid = TopologyBlock {
            adjacency: cosmolkit_model::AdjacencyList::from_topology(0, &[]),
            ..tagged.clone()
        };
        assert!(matches!(
            cleanup_invalid_atropisomers_with_ring_state(
                &invalid,
                &hybs(all(6, Hybridization::Sp2)),
                None
            ),
            Err(super::AtropisomerCleanupError::Algorithm(
                AtropisomerError::InvalidTopology { .. }
            ))
        ));
        // A CONSUMED SSSR with wrong membership dimensions is rejected by
        // the existing validate_rings causes.
        let atom_rows_wrong = RingInfo::new(RingFindType::Sssr, 7, 6);
        assert!(matches!(
            cleanup_invalid_atropisomers_with_ring_state(
                &tagged,
                &hybs(all(6, Hybridization::Sp2)),
                Some(atom_rows_wrong)
            ),
            Err(super::AtropisomerCleanupError::Algorithm(
                AtropisomerError::RingAtomRowCount {
                    actual: 7,
                    expected: 6
                }
            ))
        ));
        let bond_rows_wrong = RingInfo::new(RingFindType::Sssr, 6, 5);
        assert!(matches!(
            cleanup_invalid_atropisomers_with_ring_state(
                &tagged,
                &hybs(all(6, Hybridization::Sp2)),
                Some(bond_rows_wrong)
            ),
            Err(super::AtropisomerCleanupError::Algorithm(
                AtropisomerError::RingBondRowCount {
                    actual: 5,
                    expected: 6
                }
            ))
        ));
        // Out-of-range membership cannot reach the wrapper from valid
        // RingInfo values: add_ring GROWS membership rows, and the
        // persisted-components constructor rejects out-of-range rows at
        // construction, so the dimension check above always fires first.
        // That structural fact is recorded here instead of an unreachable
        // assertion.
        // No-tag malformed supplied state stays unconsumed and is returned
        // exactly.
        let untagged = cycle(6, BondStereo::None);
        let malformed = RingInfo::new(RingFindType::Sssr, 7, 6);
        let malformed_snapshot = malformed.clone();
        let find_before = probe::find_calls();
        let (output, returned) = cleanup_invalid_atropisomers_with_ring_state(
            &untagged,
            &hybs(all(6, Hybridization::Sp2)),
            Some(malformed),
        )
        .unwrap();
        assert_eq!(returned, Some(malformed_snapshot), "no-tag state consumed");
        assert_eq!(probe::find_calls() - find_before, 0, "no-tag find");
        assert_eq!(output.bonds[0].stereo(), BondStereo::None);
        // Vocabulary control ONLY: constructed errors prove Display and the
        // borrowed Error::source chain, never a real acquisition failure.
        // A canonical find failure is unreachable from these valid inputs.
        let rings_error = crate::RingFindingError::Value { message: "probe" };
        let cleanup_error = super::AtropisomerCleanupError::Rings(rings_error);
        assert_eq!(
            cleanup_error.to_string(),
            "ring finding failed during atropisomer cleanup: probe"
        );
        assert!(std::error::Error::source(&cleanup_error).is_some());
        let algorithm_error =
            super::AtropisomerCleanupError::Algorithm(AtropisomerError::RingInfoNotSssr);
        assert_eq!(
            algorithm_error.to_string(),
            "atropisomer cleanup failed: ring information must be SSSR or better"
        );
        assert!(std::error::Error::source(&algorithm_error).is_some());
    }
}

#[must_use]
pub fn does_topology_have_atropisomers(topology: &TopologyBlock) -> bool {
    // Complete pinned source: doesMolHaveAtropisomers.
    // RDKit✔️✔️: bool doesMolHaveAtropisomers(const ROMol &mol) {
    // RDKit✔️✔️:   for (auto bond : mol.bonds()) {
    // RDKit✔️✔️:     auto bondStereo = bond->getStereo();
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (bondStereo == Bond::BondStereo::STEREOATROPCW ||
    // RDKit✔️✔️:         bondStereo == Bond::BondStereo::STEREOATROPCCW) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    topology
        .bonds
        .iter()
        .any(|bond| matches!(bond.stereo(), BondStereo::AtropCw | BondStereo::AtropCcw))
}

#[cfg(test)]
mod source_atrop_pointer_order_cases {
    use super::*;
    use cosmolkit_model::{Atom, AtomSpec, BondSpec};
    use cosmolkit_types::Element;

    #[test]
    fn detached_candidate_traversal_checks_every_identity_without_sorting_or_repair() {
        let atoms = (0..4)
            .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
            .collect();
        let edges = [
            (1, 0, BondOrder::Single, BondDirection::BeginWedge),
            (1, 2, BondOrder::Single, BondDirection::None),
            (1, 3, BondOrder::Double, BondDirection::None),
        ];
        let bonds = edges
            .into_iter()
            .enumerate()
            .map(|(i, (a, b, o, d))| {
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(a), AtomId::new(b), o).with_direction(d),
                )
            })
            .collect();
        let graph = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap();
        let before = graph.clone();
        // These are explicit abstract source-state witnesses, not claims that
        // either sequence was captured from a particular native execution.
        let descending = [BondId::new(2), BondId::new(1)];
        let checked =
            SourceAtropisomerCandidateTraversal::try_from_parts(&graph, &descending).unwrap();
        assert_eq!(checked.candidates, &descending);
        assert!(std::ptr::eq(checked.topology, &graph));
        let invalid = [
            (
                vec![BondId::new(1), BondId::new(1), BondId::new(2)],
                AtropisomerError::CandidateTraversalDuplicate {
                    bond: BondId::new(1),
                },
            ),
            (
                vec![BondId::new(1)],
                AtropisomerError::CandidateTraversalMissing {
                    bond: BondId::new(2),
                },
            ),
            (
                vec![BondId::new(1), BondId::new(2), BondId::new(0)],
                AtropisomerError::CandidateTraversalUnexpected {
                    bond: BondId::new(0),
                },
            ),
            (
                vec![BondId::new(1), BondId::new(2), BondId::new(3)],
                AtropisomerError::CandidateTraversalOutOfRange {
                    bond: BondId::new(3),
                    bond_count: 3,
                },
            ),
        ];
        for (order, expected) in invalid {
            assert_eq!(
                SourceAtropisomerCandidateTraversal::try_from_parts(&graph, &order).unwrap_err(),
                expected,
            );
            assert_eq!(graph, before);
        }
        let empty = TopologyBlock::default();
        let checked = SourceAtropisomerCandidateTraversal::try_from_parts(&empty, &[]).unwrap();
        assert!(checked.candidates.is_empty());
        assert_eq!(
            detect_atropisomer_from_candidate_order(checked, None).unwrap(),
            AtropisomerAssignment::default()
        );
    }
    #[test]
    fn same_source_candidate_body_exposes_pointer_prefix_dependence_in_final_chemistry() {
        let elements = [Element::C, Element::C, Element::C, Element::C, Element::O];
        let cache = [
            SourceAtomValenceFacts {
                explicit_valence: 1,
                implicit_valence: 3,
            },
            SourceAtomValenceFacts {
                explicit_valence: 4,
                implicit_valence: 0,
            },
            SourceAtomValenceFacts {
                explicit_valence: 2,
                implicit_valence: 0,
            },
            SourceAtomValenceFacts {
                explicit_valence: 1,
                implicit_valence: 3,
            },
            SourceAtomValenceFacts::UNINITIALIZED,
        ];
        let atoms = elements
            .into_iter()
            .enumerate()
            .map(|(i, e)| {
                let mut atom = Atom::from_spec(
                    AtomId::new(i),
                    AtomSpec::new(e).with_hybridization(if i == 0 || i == 3 {
                        Hybridization::Sp3
                    } else {
                        Hybridization::Sp2
                    }),
                );
                atom.set_source_valence_facts(cache[i]);
                atom
            })
            .collect();
        let edges = [
            (1, 0, BondOrder::Single, BondDirection::BeginWedge),
            (1, 2, BondOrder::Single, BondDirection::None),
            (2, 3, BondOrder::Single, BondDirection::BeginDash),
            (1, 4, BondOrder::Double, BondDirection::None),
        ];
        let bonds = edges
            .into_iter()
            .enumerate()
            .map(|(i, (a, b, o, d))| {
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(a), AtomId::new(b), o).with_direction(d),
                )
            })
            .collect();
        let graph = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap();
        let before = graph.clone();
        // The cached C2 implicit zero is a legal source state: setNoImplicit
        // true, cache(false), restore false (the native setter retains cache).
        // A=double bond3 to the last dirty atom; B=eligible single axis bond1.
        let a_then_b = detect_atropisomer_from_candidate_order(
            SourceAtropisomerCandidateTraversal::try_from_parts(
                &graph,
                &[BondId::new(3), BondId::new(1)],
            )
            .unwrap(),
            None,
        )
        .unwrap();
        let b_then_a = detect_atropisomer_from_candidate_order(
            SourceAtropisomerCandidateTraversal::try_from_parts(
                &graph,
                &[BondId::new(1), BondId::new(3)],
            )
            .unwrap(),
            None,
        )
        .unwrap();
        assert_eq!(a_then_b.hybridization, None);
        assert_eq!(a_then_b.conjugated_bonds, None);
        assert_eq!(
            a_then_b.atom_valence_updates,
            vec![(
                AtomId::new(4),
                SourceAtomValenceFacts {
                    explicit_valence: 2,
                    implicit_valence: 0
                }
            )]
        );
        assert_eq!(
            a_then_b.bond_updates,
            vec![AtropisomerBondUpdate {
                bond: BondId::new(1),
                stereo: BondStereo::AtropCcw
            }]
        );
        assert_eq!(
            b_then_a.hybridization.as_ref().unwrap().values[2],
            Hybridization::Sp3
        );
        assert!(b_then_a.conjugated_bonds.is_some());
        assert!(b_then_a.bond_updates.is_empty());
        assert_eq!(
            b_then_a.atom_valence_updates[2],
            (
                AtomId::new(2),
                SourceAtomValenceFacts {
                    explicit_valence: 2,
                    implicit_valence: 2
                }
            )
        );
        assert_eq!(graph, before);
    }

    #[test]
    fn fresh_all_dummy_parser_state_also_has_pointer_dependent_conjugation_effects() {
        // Fresh native *(*)(:*)*(*)* graph, all source cache rows UNINITIALIZED
        // and all hybridizations UNSPECIFIED. Native atomicNum!=0 condition
        // excludes every atom from the unspecified-hybridization guard.
        let dummy = Element::from_atomic_number(0).unwrap();
        let atoms = (0..6)
            .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(dummy)))
            .collect();
        let edges = [
            (0, 1, BondOrder::Single, BondDirection::BeginWedge),
            (0, 2, BondOrder::Aromatic, BondDirection::None),
            (0, 3, BondOrder::Single, BondDirection::BeginWedge),
            (3, 4, BondOrder::Single, BondDirection::BeginWedge),
            (3, 5, BondOrder::Single, BondDirection::BeginWedge),
        ];
        let bonds = edges
            .into_iter()
            .enumerate()
            .map(|(i, (a, b, o, d))| {
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(a), AtomId::new(b), o)
                        .with_direction(d)
                        .with_aromatic(o == BondOrder::Aromatic),
                )
            })
            .collect();
        let graph = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap();
        let before = graph.clone();
        let leaves_before_axis = detect_atropisomer_from_candidate_order(
            SourceAtropisomerCandidateTraversal::try_from_parts(
                &graph,
                &[
                    BondId::new(0),
                    BondId::new(1),
                    BondId::new(3),
                    BondId::new(4),
                    BondId::new(2),
                ],
            )
            .unwrap(),
            None,
        )
        .unwrap();
        let axis_before_leaves = detect_atropisomer_from_candidate_order(
            SourceAtropisomerCandidateTraversal::try_from_parts(
                &graph,
                &[
                    BondId::new(2),
                    BondId::new(0),
                    BondId::new(1),
                    BondId::new(3),
                    BondId::new(4),
                ],
            )
            .unwrap(),
            None,
        )
        .unwrap();
        assert_eq!(leaves_before_axis.conjugated_bonds, None);
        assert_eq!(leaves_before_axis.hybridization, None);
        assert!(!graph.bonds[1].is_conjugated());
        assert_eq!(
            axis_before_leaves.conjugated_bonds,
            Some(vec![false, true, false, false, false])
        );
        assert_eq!(
            axis_before_leaves.hybridization.unwrap().values,
            vec![Hybridization::Unspecified; 6]
        );
        assert!(
            leaves_before_axis.bond_updates.is_empty()
                && axis_before_leaves.bond_updates.is_empty()
        );
        assert_eq!(graph, before);
    }
}

#[cfg(test)]
mod source610_carrier_output_tests {
    use super::*;
    use cosmolkit_model::{Atom, AtomSpec, BondSpec};
    use cosmolkit_types::Element;

    fn graph(edges: &[(usize, usize)]) -> TopologyBlock {
        let atoms = (0..6)
            .map(|id| Atom::from_spec(AtomId::new(id), AtomSpec::new(Element::C)))
            .collect();
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(id, &(begin, end))| {
                Bond::from_spec(
                    BondId::new(id),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
    }
    fn ends(left: &[usize], right: &[usize]) -> [AtropEnd; 2] {
        [
            AtropEnd {
                atom: AtomId::new(5),
                bonds: left.iter().copied().map(BondId::new).collect(),
            },
            AtropEnd {
                atom: AtomId::new(5),
                bonds: right.iter().copied().map(BondId::new).collect(),
            },
        ]
    }
    #[test]
    fn source610_appends_existing_outputs_and_only_sorts_exactly_two() {
        let topology = graph(&[(0, 1), (0, 4), (0, 2), (1, 5), (1, 3)]);
        let mut output = ends(&[1], &[3]);
        assert!(atropisomer_ends(&topology, &topology.bonds[0], &mut output).unwrap());
        assert_eq!(
            [output[0].atom, output[1].atom],
            [AtomId::new(0), AtomId::new(1)]
        );
        assert_eq!(
            output[0].bonds,
            vec![BondId::new(1), BondId::new(1), BondId::new(2)]
        );
        assert_eq!(
            output[1].bonds,
            vec![BondId::new(3), BondId::new(3), BondId::new(4)]
        );
        let fresh = atropisomer_ends_fresh(&topology, &topology.bonds[0])
            .unwrap()
            .unwrap();
        assert_eq!(fresh[0].bonds, vec![BondId::new(2), BondId::new(1)]);
        assert_eq!(fresh[1].bonds, vec![BondId::new(4), BondId::new(3)]);
    }
    #[test]
    fn source610_false_keeps_both_focus_writes_and_visited_end_prefix() {
        let topology = graph(&[(0, 1), (0, 4), (0, 2)]);
        let mut output = ends(&[], &[]);
        assert!(!atropisomer_ends(&topology, &topology.bonds[0], &mut output).unwrap());
        assert_eq!(
            [output[0].atom, output[1].atom],
            [AtomId::new(0), AtomId::new(1)]
        );
        assert_eq!(output[0].bonds, vec![BondId::new(2), BondId::new(1)]);
        assert!(output[1].bonds.is_empty());
        let mut seeded = ends(&[], &[0]);
        assert!(atropisomer_ends(&topology, &topology.bonds[0], &mut seeded).unwrap());
        assert_eq!(seeded[1].bonds, vec![BondId::new(0)]);
    }
    #[test]
    fn source610_first_empty_end_does_not_visit_or_reorder_second_output() {
        let topology = graph(&[(0, 1), (1, 5), (1, 3)]);
        let mut output = ends(&[], &[1, 2]);
        assert!(!atropisomer_ends(&topology, &topology.bonds[0], &mut output).unwrap());
        assert_eq!(
            [output[0].atom, output[1].atom],
            [AtomId::new(0), AtomId::new(1)]
        );
        assert_eq!(output[1].bonds, vec![BondId::new(1), BondId::new(2)]);
    }
}

#[cfg(test)]
mod source614_no_conf_direction_tests {
    use super::*;
    #[test]
    fn source614_eight_native_parity_cells() {
        let expected = [
            (
                BondStereo::AtropCw,
                [
                    [BondDirection::BeginDash, BondDirection::BeginWedge],
                    [BondDirection::BeginWedge, BondDirection::BeginDash],
                ],
            ),
            (
                BondStereo::AtropCcw,
                [
                    [BondDirection::BeginWedge, BondDirection::BeginDash],
                    [BondDirection::BeginDash, BondDirection::BeginWedge],
                ],
            ),
        ];
        for (stereo, rows) in expected {
            for (end, row) in rows.into_iter().enumerate() {
                for (bond, direction) in row.into_iter().enumerate() {
                    assert_eq!(no_conf_direction(stereo, end, bond).unwrap(), direction);
                }
            }
        }
    }
    #[test]
    fn source614_preconditions_fail_in_source_order_without_parity_fallback() {
        for (stereo, end, bond, expected) in [
            (
                BondStereo::None,
                usize::MAX,
                usize::MAX,
                "whichEnd must be 0 or 1",
            ),
            (BondStereo::None, 1, usize::MAX, "whichBond must be 0 or 1"),
            (
                BondStereo::None,
                1,
                1,
                "bondStereo must be BondAtropisomerCW or BondAtropisomerCCW",
            ),
            (
                BondStereo::Any,
                0,
                0,
                "bondStereo must be BondAtropisomerCW or BondAtropisomerCCW",
            ),
        ] {
            assert!(
                matches!(no_conf_direction(stereo,end,bond), Err(PerceptionError::SourcePrecondition {message}) if message == expected)
            );
        }
    }
}

#[cfg(test)]
mod source618_no_conf_wedge_tests {
    use super::*;
    use crate::RingFindType;
    use cosmolkit_model::{Atom, AtomSpec, BondSpec};
    use cosmolkit_types::Element;
    fn graph(
        edges: &[(usize, usize, BondDirection, BondOrder)],
        stereo: BondStereo,
    ) -> TopologyBlock {
        let atoms = (0..6)
            .map(|id| Atom::from_spec(AtomId::new(id), AtomSpec::new(Element::C)))
            .collect();
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(id, &(begin, end, direction, order))| {
                let spec = BondSpec::new(AtomId::new(begin), AtomId::new(end), order)
                    .with_direction(direction);
                Bond::from_spec(
                    BondId::new(id),
                    if id == 0 {
                        spec.with_stereo(stereo)
                    } else {
                        spec
                    },
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
    }
    fn sparse() -> RingInfo {
        RingInfo::new(RingFindType::Sssr, 0, 0)
    }
    #[test]
    fn source618_existing_wedges_write_directions_without_map_or_ring_reads() {
        let mut topology = graph(
            &[
                (0, 1, BondDirection::None, BondOrder::Single),
                (0, 2, BondDirection::BeginWedge, BondOrder::Single),
                (1, 3, BondDirection::BeginDash, BondOrder::Single),
            ],
            BondStereo::AtropCw,
        );
        let mut wedges = WedgeAssignments::default();
        wedges.insert_atropisomer_source(AtropisomerWedgeUpdate {
            bond: BondId::new(1),
            begin: AtomId::new(0),
            end: AtomId::new(2),
            direction: BondDirection::BeginWedge,
            atropisomer_bond: BondId::new(2),
        });
        let before = wedges.clone();
        let mut rings = sparse();
        rings.reset();
        assert!(
            wedge_atropisomer_no_conformer_source(
                &mut topology,
                &rings,
                BondId::new(0),
                &mut wedges
            )
            .unwrap()
        );
        assert_eq!(topology.bonds[1].direction(), BondDirection::BeginDash);
        assert_eq!(topology.bonds[2].direction(), BondDirection::BeginWedge);
        assert_eq!(wedges, before);
    }
    #[test]
    fn source618_sparse_ring_cache_selects_wedge_orients_and_inserts_real_map() {
        let mut topology = graph(
            &[
                (0, 1, BondDirection::None, BondOrder::Single),
                (2, 0, BondDirection::None, BondOrder::Single),
                (3, 1, BondDirection::None, BondOrder::Single),
            ],
            BondStereo::AtropCw,
        );
        let adjacency = topology.adjacency.neighbors_of(1).as_ptr();
        let mut wedges = WedgeAssignments::default();
        assert!(
            wedge_atropisomer_no_conformer_source(
                &mut topology,
                &sparse(),
                BondId::new(0),
                &mut wedges
            )
            .unwrap()
        );
        assert_eq!(topology.bonds[1].direction(), BondDirection::None);
        assert_eq!(
            (
                topology.bonds[2].begin(),
                topology.bonds[2].end(),
                topology.bonds[2].direction()
            ),
            (AtomId::new(1), AtomId::new(3), BondDirection::BeginWedge)
        );
        assert_eq!(topology.adjacency.neighbors_of(1).as_ptr(), adjacency);
        assert_eq!(wedges.iter().count(), 1);
        assert!(
            matches!(wedges.get(BondId::new(2)),Some(WedgeInfo::Atropisomer {update}) if update.atropisomer_bond==BondId::new(0))
        );
    }
    #[test]
    fn source618_missing_unknown_and_occupied_candidates_return_source_false() {
        for edges in [
            vec![
                (0, 1, BondDirection::None, BondOrder::Single),
                (0, 2, BondDirection::None, BondOrder::Single),
            ],
            vec![
                (0, 1, BondDirection::None, BondOrder::Single),
                (0, 2, BondDirection::Unknown, BondOrder::Single),
                (1, 3, BondDirection::None, BondOrder::Single),
            ],
            vec![
                (0, 1, BondDirection::None, BondOrder::Single),
                (0, 2, BondDirection::None, BondOrder::Double),
                (1, 3, BondDirection::None, BondOrder::Triple),
            ],
        ] {
            let mut topology = graph(&edges, BondStereo::AtropCw);
            let before = topology.clone();
            let mut wedges = WedgeAssignments::default();
            assert!(
                !wedge_atropisomer_no_conformer_source(
                    &mut topology,
                    &sparse(),
                    BondId::new(0),
                    &mut wedges
                )
                .unwrap()
            );
            assert_eq!(topology, before);
            assert_eq!(wedges.iter().count(), 0);
        }
        let mut topology = graph(
            &[
                (0, 1, BondDirection::None, BondOrder::Single),
                (0, 2, BondDirection::None, BondOrder::Single),
                (1, 3, BondDirection::None, BondOrder::Single),
            ],
            BondStereo::AtropCw,
        );
        let mut wedges = WedgeAssignments::default();
        for id in [1, 2] {
            let b = &topology.bonds[id];
            wedges.insert_atropisomer_source(AtropisomerWedgeUpdate {
                bond: b.id(),
                begin: b.begin(),
                end: b.end(),
                direction: b.direction(),
                atropisomer_bond: BondId::new(0),
            });
        }
        let before = wedges.clone();
        let mut rings = sparse();
        rings.reset();
        assert!(
            !wedge_atropisomer_no_conformer_source(
                &mut topology,
                &rings,
                BondId::new(0),
                &mut wedges
            )
            .unwrap()
        );
        assert_eq!(wedges, before);
    }
    #[test]
    fn source618_late_direction_precondition_preserves_earlier_direct_write() {
        let mut topology = graph(
            &[
                (0, 1, BondDirection::None, BondOrder::Single),
                (0, 2, BondDirection::BeginWedge, BondOrder::Single),
                (0, 3, BondDirection::None, BondOrder::Single),
                (0, 4, BondDirection::BeginWedge, BondOrder::Single),
                (1, 5, BondDirection::None, BondOrder::Single),
            ],
            BondStereo::AtropCw,
        );
        let mut wedges = WedgeAssignments::default();
        let mut rings = sparse();
        rings.reset();
        assert!(matches!(
            wedge_atropisomer_no_conformer_source(
                &mut topology,
                &rings,
                BondId::new(0),
                &mut wedges
            ),
            Err(AtropisomerError::SourcePrecondition {
                message: "whichBond must be 0 or 1"
            })
        ));
        assert_eq!(topology.bonds[1].direction(), BondDirection::BeginDash);
        assert_eq!(topology.bonds[3].direction(), BondDirection::BeginWedge);
        assert_eq!(wedges.iter().count(), 0);
    }
    #[test]
    fn source618_single_preference_precedes_wedge_preference_and_map_occupancy() {
        let mut topology = graph(
            &[
                (0, 1, BondDirection::None, BondOrder::Single),
                (0, 2, BondDirection::None, BondOrder::Single),
                (1, 3, BondDirection::None, BondOrder::Aromatic),
            ],
            BondStereo::AtropCw,
        );
        let mut wedges = WedgeAssignments::default();
        assert!(
            wedge_atropisomer_no_conformer_source(
                &mut topology,
                &sparse(),
                BondId::new(0),
                &mut wedges
            )
            .unwrap()
        );
        assert_eq!(topology.bonds[1].direction(), BondDirection::BeginDash);
        assert_eq!(topology.bonds[2].direction(), BondDirection::None);
        assert!(wedges.get(BondId::new(1)).is_some());
    }
    #[test]
    fn source618_projected_existing_wedges_do_not_become_new_native_map_entries() {
        let topology = graph(
            &[
                (0, 1, BondDirection::None, BondOrder::Single),
                (0, 2, BondDirection::BeginWedge, BondOrder::Single),
                (1, 3, BondDirection::BeginDash, BondOrder::Single),
            ],
            BondStereo::AtropCw,
        );
        let mut updates = BTreeMap::new();
        let mut writes = BTreeSet::new();
        let occupied = BTreeSet::new();
        let mut state = ProjectedAtropWedgeState {
            topology: &topology,
            occupied: &occupied,
            updates: &mut updates,
            map_writes: &mut writes,
        };
        wedge_one_no_conformer_source(&mut state, &sparse(), BondId::new(0)).unwrap();
        assert_eq!(updates.len(), 2);
        assert!(writes.is_empty());
        let assignment = AtropisomerWedgeAssignment {
            source_map_writes: writes.into_iter().collect(),
            bond_updates: updates.into_values().collect(),
            diagnostics: Vec::new(),
        };
        assert_eq!(
            WedgeAssignments::from_atropisomer_wedge_assignment(assignment)
                .iter()
                .count(),
            0
        );
    }
}

#[cfg(test)]
mod source622_frame_output_tests {
    use super::*;
    use cosmolkit_model::BondSpec;
    fn axis() -> Bond {
        Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        )
    }
    #[test]
    fn source622_short_axis_writes_raw_x_and_preserves_other_actual_outputs() {
        let conf = Conformer3D::new(0, vec![[0.0; 3], [5e-8, 0.0, 0.0]], true);
        let (mut x, mut y, mut z) = ([9.0; 3], [8.0; 3], [7.0; 3]);
        assert!(
            !frame_of_reference_source(
                &axis(),
                AtropisomerConformer::ThreeD(&conf),
                &mut x,
                &mut y,
                &mut z
            )
            .unwrap()
        );
        assert_eq!(x, [5e-8, 0.0, 0.0]);
        assert_eq!(y, [8.0; 3]);
        assert_eq!(z, [7.0; 3]);
    }
    #[test]
    fn source622_exact_length_boundary_and_flagged_2d_xyz_use_actual_source_axes() {
        let conf = Conformer3D::new(0, vec![[0.0; 3], [1e-7, 0.0, 0.0]], false);
        let (mut x, mut y, mut z) = ([9.0; 3], [8.0; 3], [7.0; 3]);
        assert!(
            frame_of_reference_source(
                &axis(),
                AtropisomerConformer::ThreeD(&conf),
                &mut x,
                &mut y,
                &mut z
            )
            .unwrap()
        );
        assert_eq!(x, [1.0, 0.0, 0.0]);
        assert_eq!(y, [-0.0, 1.0, 0.0]);
        assert_eq!(z, [0.0, 0.0, 1.0]);
        assert_eq!(y[0].to_bits(), (-0.0f64).to_bits());
        let conf = Conformer3D::new(0, vec![[0.0; 3], [3.0, 4.0, 12.0]], false);
        assert!(
            frame_of_reference_source(
                &axis(),
                AtropisomerConformer::ThreeD(&conf),
                &mut x,
                &mut y,
                &mut z
            )
            .unwrap()
        );
        assert_eq!(x, [3.0 / 13.0, 4.0 / 13.0, 12.0 / 13.0]);
        assert!((y[0] + 0.8).abs() < 1e-15 && (y[1] - 0.6).abs() < 1e-15);
        assert_eq!(z, [0.0, 0.0, 1.0]);
    }
    #[test]
    fn source622_three_dimensional_seed_axes_follow_source_z_threshold() {
        for (coordinate, expected_x, expected_y, expected_z) in [
            (
                [2.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ),
            (
                [0.0, 0.0, 2.0],
                [0.0, 0.0, 1.0],
                [0.0, -1.0, 0.0],
                [1.0, 0.0, -0.0],
            ),
        ] {
            let conf = Conformer3D::new(0, vec![[0.0; 3], coordinate], true);
            let (mut x, mut y, mut z) = ([9.0; 3], [8.0; 3], [7.0; 3]);
            assert!(
                frame_of_reference_source(
                    &axis(),
                    AtropisomerConformer::ThreeD(&conf),
                    &mut x,
                    &mut y,
                    &mut z
                )
                .unwrap()
            );
            assert_eq!(x, expected_x);
            assert_eq!(y, expected_y);
            assert_eq!(z, expected_z);
        }
    }
    #[test]
    fn source622_ieee_nan_follows_native_comparisons_and_normalization() {
        let conf = Conformer3D::new(0, vec![[0.0; 3], [f64::NAN, 0.0, 1.0]], true);
        let (mut x, mut y, mut z) = ([9.0; 3], [8.0; 3], [7.0; 3]);
        assert!(
            frame_of_reference_source(
                &axis(),
                AtropisomerConformer::ThreeD(&conf),
                &mut x,
                &mut y,
                &mut z
            )
            .unwrap()
        );
        assert!(x.into_iter().chain(y).chain(z).all(f64::is_nan));
    }
    #[test]
    fn source622_normalization_error_keeps_raw_2d_y_and_previous_z() {
        let conf = Conformer3D::new(0, vec![[0.0; 3], [0.0, 0.0, 2.0]], false);
        let (mut x, mut y, mut z) = ([9.0; 3], [8.0; 3], [7.0; 3]);
        assert!(
            frame_of_reference_source(
                &axis(),
                AtropisomerConformer::ThreeD(&conf),
                &mut x,
                &mut y,
                &mut z
            )
            .is_err()
        );
        assert_eq!(x, [0.0, 0.0, 1.0]);
        assert_eq!(y, [-0.0, 0.0, 0.0]);
        assert_eq!(z, [7.0; 3]);
    }
}

#[cfg(test)]
mod source626_two_d_direction_tests {
    use super::*;
    #[test]
    fn source626_native_parity_other_end_y_and_ieee_cells() {
        for (stereo, rows) in [
            (
                BondStereo::AtropCw,
                [
                    [BondDirection::BeginDash, BondDirection::BeginWedge],
                    [BondDirection::BeginWedge, BondDirection::BeginDash],
                ],
            ),
            (
                BondStereo::AtropCcw,
                [
                    [BondDirection::BeginWedge, BondDirection::BeginDash],
                    [BondDirection::BeginDash, BondDirection::BeginWedge],
                ],
            ),
        ] {
            for (end, row) in rows.into_iter().enumerate() {
                for (bond, base) in row.into_iter().enumerate() {
                    for (other_y, negative) in [
                        (1.0, false),
                        (-1.0, true),
                        (0.0, false),
                        (-0.0, false),
                        (f64::NAN, false),
                        (f64::INFINITY, false),
                        (f64::NEG_INFINITY, true),
                    ] {
                        let mut vectors = [[0.0; 3]; 2];
                        vectors[end][1] = -other_y;
                        vectors[1 - end][1] = other_y;
                        let expected = if negative {
                            match base {
                                BondDirection::BeginWedge => BondDirection::BeginDash,
                                BondDirection::BeginDash => BondDirection::BeginWedge,
                                _ => unreachable!(),
                            }
                        } else {
                            base
                        };
                        assert_eq!(
                            two_d_direction(vectors, stereo, end, bond).unwrap(),
                            expected,
                            "{stereo:?}, end={end}, bond={bond}, otherY={other_y}"
                        );
                    }
                }
            }
        }
    }
    #[test]
    fn source626_preconditions_precede_other_end_index_and_stereo_fallback() {
        for (stereo, end, bond, expected) in [
            (
                BondStereo::None,
                usize::MAX,
                usize::MAX,
                "whichEnd must be 0 or 1",
            ),
            (BondStereo::None, 0, usize::MAX, "whichBond must be 0 or 1"),
            (
                BondStereo::None,
                0,
                0,
                "bondStereo must be BondAtropisomerCW or BondAtropisomerCCW",
            ),
        ] {
            assert!(
                matches!(two_d_direction([[f64::NAN;3];2],stereo,end,bond),Err(PerceptionError::SourcePrecondition {message}) if message==expected)
            );
        }
    }
}

#[cfg(test)]
mod source630_end_vector_output_tests {
    use super::*;
    use cosmolkit_model::{Atom, AtomSpec, BondSpec};
    use cosmolkit_types::Element;
    fn graph() -> TopologyBlock {
        let atoms = (0..3)
            .map(|id| Atom::from_spec(AtomId::new(id), AtomSpec::new(Element::C)))
            .collect();
        let bonds = (1..3)
            .map(|id| {
                Bond::from_spec(
                    BondId::new(id - 1),
                    BondSpec::new(AtomId::new(0), AtomId::new(id), BondOrder::Single),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
    }
    fn frame() -> Frame {
        Frame {
            y: [0.0, 1.0, 0.0],
            z: [0.0, 0.0, 1.0],
        }
    }
    fn end(count: usize) -> AtropEnd {
        AtropEnd {
            atom: AtomId::new(0),
            bonds: (0..count).map(BondId::new).collect(),
        }
    }
    #[test]
    fn source630_native_precondition_preserves_existing_output_before_reads() {
        let topology = graph();
        let conf = Conformer3D::new(0, vec![[0.0; 3]; 3], true);
        for count in [0, 3] {
            let mut output = [7.0; 3];
            assert!(
                matches!(end_vector_source(&topology,&end(count),frame(),AtropisomerConformer::ThreeD(&conf),true,&mut output),Err(PerceptionError::CarrierCount {count: actual,..}) if actual==count)
            );
            assert_eq!(output, [7.0; 3]);
        }
    }
    #[test]
    fn source630_same_side_failure_keeps_first_projection() {
        let topology = graph();
        let conf = Conformer3D::new(0, vec![[0.0; 3], [2.0, 3.0, 4.0], [1.0, 6.0, 8.0]], true);
        let mut output = [7.0; 3];
        assert!(matches!(
            end_vector_source(
                &topology,
                &end(2),
                frame(),
                AtropisomerConformer::ThreeD(&conf),
                true,
                &mut output
            ),
            Err(PerceptionError::Rejected(
                AtropisomerRejectionKind::SameSideCarriers
            ))
        ));
        assert_eq!(output, [0.0, 3.0, 4.0]);
    }
    #[test]
    fn source630_collinear_first_uses_opposite_second_before_normalizing() {
        let topology = graph();
        let conf = Conformer3D::new(0, vec![[0.0; 3], [2.0, 0.0, 0.0], [1.0, 3.0, 4.0]], true);
        let mut output = [7.0; 3];
        end_vector_source(
            &topology,
            &end(2),
            frame(),
            AtropisomerConformer::ThreeD(&conf),
            true,
            &mut output,
        )
        .unwrap();
        assert_eq!(output, [-0.0, -0.6, -0.8]);
        assert_eq!(output[0].to_bits(), (-0.0f64).to_bits());
        end_vector_source(
            &topology,
            &end(2),
            frame(),
            AtropisomerConformer::ThreeD(&conf),
            false,
            &mut output,
        )
        .unwrap();
        assert_eq!(output, [-0.0, -3.0, -4.0]);
    }
    #[test]
    fn source630_collinear_false_retains_opposite_zero_and_single_projection() {
        let topology = graph();
        let conf = Conformer3D::new(0, vec![[0.0; 3], [2.0, 0.0, 0.0], [1.0, 0.0, 0.0]], true);
        let mut output = [7.0; 3];
        assert!(matches!(
            end_vector_source(
                &topology,
                &end(2),
                frame(),
                AtropisomerConformer::ThreeD(&conf),
                true,
                &mut output
            ),
            Err(PerceptionError::Rejected(
                AtropisomerRejectionKind::CollinearCarrier
            ))
        ));
        assert_eq!(output.map(f64::to_bits), [(-0.0f64).to_bits(); 3]);
        assert!(matches!(
            end_vector_source(
                &topology,
                &end(1),
                frame(),
                AtropisomerConformer::ThreeD(&conf),
                true,
                &mut output
            ),
            Err(PerceptionError::Rejected(
                AtropisomerRejectionKind::CollinearCarrier
            ))
        ));
        assert_eq!(output.map(f64::to_bits), [0.0f64.to_bits(); 3]);
    }
    #[test]
    fn source630_exact_projection_length_threshold_and_ieee_nan() {
        let topology = graph();
        let conf = Conformer3D::new(0, vec![[0.0; 3], [0.0, 1e-7, 0.0], [0.0; 3]], true);
        let mut output = [7.0; 3];
        end_vector_source(
            &topology,
            &end(1),
            frame(),
            AtropisomerConformer::ThreeD(&conf),
            true,
            &mut output,
        )
        .unwrap();
        assert_eq!(output, [0.0, 1.0, 0.0]);
        let conf = Conformer3D::new(0, vec![[0.0; 3], [f64::NAN, 1.0, 0.0], [0.0; 3]], true);
        end_vector_source(
            &topology,
            &end(1),
            frame(),
            AtropisomerConformer::ThreeD(&conf),
            true,
            &mut output,
        )
        .unwrap();
        assert!(output.iter().all(|value| value.is_nan()));
    }
}

#[cfg(test)]
mod source630_other_atom_source_error_tests {
    use super::*;
    use cosmolkit_model::{Atom, AtomSpec, BondSpec};
    use cosmolkit_types::Element;
    fn graph() -> TopologyBlock {
        let atoms = (0..4)
            .map(|id| Atom::from_spec(AtomId::new(id), AtomSpec::new(Element::C)))
            .collect();
        let bonds = [(0, 1), (0, 2), (1, 3)]
            .into_iter()
            .enumerate()
            .map(|(id, (begin, end))| {
                Bond::from_spec(
                    BondId::new(id),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
    }
    #[test]
    fn source630_reused_carrier_output_preserves_appends_when_other_atom_throws() {
        let topology = graph();
        let mut ends = [
            AtropEnd {
                atom: AtomId::new(3),
                bonds: vec![BondId::new(2)],
            },
            AtropEnd {
                atom: AtomId::new(3),
                bonds: Vec::new(),
            },
        ];
        let error = atropisomer_ends(&topology, &topology.bonds[0], &mut ends).unwrap_err();
        assert_eq!((error.bond, error.atom), (BondId::new(2), AtomId::new(0)));
        assert_eq!(ends[0].atom, AtomId::new(0));
        assert_eq!(ends[1].atom, AtomId::new(1));
        assert_eq!(ends[0].bonds, vec![BondId::new(2), BondId::new(1)]);
        assert!(ends[1].bonds.is_empty());
    }
    #[test]
    fn source630_second_other_atom_error_preserves_first_projected_vector() {
        let topology = graph();
        let conf = Conformer3D::new(
            0,
            vec![[0.0; 3], [1.0, 0.0, 0.0], [2.0, 3.0, 4.0], [3.0; 3]],
            true,
        );
        let end = AtropEnd {
            atom: AtomId::new(0),
            bonds: vec![BondId::new(1), BondId::new(2)],
        };
        let mut output = [7.0; 3];
        assert!(matches!(
            end_vector_source(
                &topology,
                &end,
                Frame {
                    y: [0.0, 1.0, 0.0],
                    z: [0.0, 0.0, 1.0]
                },
                AtropisomerConformer::ThreeD(&conf),
                true,
                &mut output
            ),
            Err(PerceptionError::SourcePrecondition {
                message: "bad index"
            })
        ));
        assert_eq!(output, [0.0, 3.0, 4.0]);
        let end = AtropEnd {
            atom: AtomId::new(0),
            bonds: vec![BondId::new(2)],
        };
        output = [7.0; 3];
        assert!(matches!(
            end_vector_source(
                &topology,
                &end,
                Frame {
                    y: [0.0, 1.0, 0.0],
                    z: [0.0, 0.0, 1.0]
                },
                AtropisomerConformer::ThreeD(&conf),
                true,
                &mut output
            ),
            Err(PerceptionError::SourcePrecondition {
                message: "bad index"
            })
        ));
        assert_eq!(output, [7.0; 3]);
    }
}

#[cfg(test)]
mod source634_two_d_wedge_tests {
    use super::*;
    use crate::RingFindType;
    use cosmolkit_model::{Atom, AtomSpec, BondSpec};
    use cosmolkit_types::Element;
    fn graph(
        edges: &[(usize, usize, BondDirection, BondOrder)],
        stereo: BondStereo,
    ) -> TopologyBlock {
        let atoms = (0..6)
            .map(|id| Atom::from_spec(AtomId::new(id), AtomSpec::new(Element::C)))
            .collect();
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(id, &(begin, end, direction, order))| {
                let spec = BondSpec::new(AtomId::new(begin), AtomId::new(end), order)
                    .with_direction(direction);
                Bond::from_spec(
                    BondId::new(id),
                    if id == 0 {
                        spec.with_stereo(stereo)
                    } else {
                        spec
                    },
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
    }
    fn sparse() -> RingInfo {
        RingInfo::new(RingFindType::Sssr, 0, 0)
    }
    fn coordinates() -> Conformer2D {
        Conformer2D::new(
            0,
            vec![
                [0.0, 0.0],
                [0.0, 1.0],
                [-1.0, 0.0],
                [-1.0, 1.0],
                [1.0, 0.0],
                [1.0, 1.0],
            ],
        )
    }
    fn edges() -> [(usize, usize, BondDirection, BondOrder); 3] {
        [
            (0, 1, BondDirection::None, BondOrder::Single),
            (2, 0, BondDirection::None, BondOrder::Single),
            (3, 1, BondDirection::None, BondOrder::Single),
        ]
    }
    #[test]
    fn source634_existing_wedges_reuse_after_geometry_without_native_map_or_ring_reads() {
        let mut topology = graph(
            &[
                (0, 1, BondDirection::None, BondOrder::Single),
                (0, 2, BondDirection::BeginWedge, BondOrder::Single),
                (1, 3, BondDirection::BeginDash, BondOrder::Single),
            ],
            BondStereo::AtropCw,
        );
        let conf = coordinates();
        let mut rings = sparse();
        rings.reset();
        let mut map = WedgeAssignments::default();
        assert!(
            wedge_atropisomer_two_d_source(
                &mut topology,
                &rings,
                BondId::new(0),
                AtropisomerConformer::TwoD(&conf),
                &mut map
            )
            .unwrap()
        );
        assert_eq!(topology.bonds[1].direction(), BondDirection::BeginDash);
        assert_eq!(topology.bonds[2].direction(), BondDirection::BeginWedge);
        assert_eq!(map.iter().count(), 0);
    }
    #[test]
    fn source634_sparse_cache_selects_preferred_wedge_and_reorients_only_choice() {
        let mut topology = graph(&edges(), BondStereo::AtropCw);
        let conf = coordinates();
        let mut map = WedgeAssignments::default();
        assert!(
            wedge_atropisomer_two_d_source(
                &mut topology,
                &sparse(),
                BondId::new(0),
                AtropisomerConformer::TwoD(&conf),
                &mut map
            )
            .unwrap()
        );
        assert_eq!(
            (
                topology.bonds[1].begin(),
                topology.bonds[1].end(),
                topology.bonds[1].direction()
            ),
            (AtomId::new(2), AtomId::new(0), BondDirection::None)
        );
        assert_eq!(
            (
                topology.bonds[2].begin(),
                topology.bonds[2].end(),
                topology.bonds[2].direction()
            ),
            (AtomId::new(1), AtomId::new(3), BondDirection::BeginWedge)
        );
        assert!(map.get(BondId::new(2)).is_some());
    }
    #[test]
    fn source634_geometry_false_precedes_existing_direction_writes() {
        for same_side in [false, true] {
            let mut es = edges().to_vec();
            es[1] = (0, 2, BondDirection::BeginWedge, BondOrder::Single);
            if same_side {
                es.push((0, 4, BondDirection::None, BondOrder::Single));
            }
            let mut topology = graph(&es, BondStereo::AtropCw);
            let before = topology.clone();
            let conf = if same_side {
                Conformer2D::new(
                    0,
                    vec![
                        [0.0, 0.0],
                        [0.0, 1.0],
                        [-1.0, 0.0],
                        [-1.0, 1.0],
                        [-2.0, 0.0],
                        [1.0, 1.0],
                    ],
                )
            } else {
                Conformer2D::new(0, vec![[0.0; 2]; 6])
            };
            let mut map = WedgeAssignments::default();
            let mut rings = sparse();
            rings.reset();
            assert!(
                !wedge_atropisomer_two_d_source(
                    &mut topology,
                    &rings,
                    BondId::new(0),
                    AtropisomerConformer::TwoD(&conf),
                    &mut map
                )
                .unwrap()
            );
            assert_eq!(topology, before);
            assert_eq!(map.iter().count(), 0);
        }
    }
    #[test]
    fn source634_slash_candidate_is_skipped_without_no_conformer_conflict_rule() {
        let mut topology = graph(
            &[
                (0, 1, BondDirection::None, BondOrder::Single),
                (0, 2, BondDirection::EndUpRight, BondOrder::Single),
                (3, 1, BondDirection::None, BondOrder::Single),
            ],
            BondStereo::AtropCw,
        );
        let conf = coordinates();
        let mut map = WedgeAssignments::default();
        assert!(
            wedge_atropisomer_two_d_source(
                &mut topology,
                &sparse(),
                BondId::new(0),
                AtropisomerConformer::TwoD(&conf),
                &mut map
            )
            .unwrap()
        );
        assert_eq!(topology.bonds[1].direction(), BondDirection::EndUpRight);
        assert!(map.get(BondId::new(2)).is_some());
    }
    #[test]
    fn source634_single_preference_precedes_equal_ring_count_wedge_preference() {
        let mut es = edges();
        es[2].3 = BondOrder::Aromatic;
        let mut topology = graph(&es, BondStereo::AtropCw);
        let conf = coordinates();
        let mut map = WedgeAssignments::default();
        assert!(
            wedge_atropisomer_two_d_source(
                &mut topology,
                &sparse(),
                BondId::new(0),
                AtropisomerConformer::TwoD(&conf),
                &mut map
            )
            .unwrap()
        );
        assert_eq!(topology.bonds[1].direction(), BondDirection::BeginDash);
        assert_eq!(topology.bonds[2].direction(), BondDirection::None);
        assert!(map.get(BondId::new(1)).is_some());
    }
    #[test]
    fn source634_actual_cycle_cache_ring_candidate_precedes_acyclic_candidate() {
        let mut topology = graph(
            &[
                (0, 1, BondDirection::None, BondOrder::Single),
                (2, 0, BondDirection::None, BondOrder::Single),
                (2, 4, BondDirection::None, BondOrder::Single),
                (4, 1, BondDirection::None, BondOrder::Single),
                (3, 1, BondDirection::None, BondOrder::Single),
            ],
            BondStereo::AtropCcw,
        );
        let rings = crate::find_sssr(&topology, &crate::RingSearchParams::default()).unwrap();
        assert_eq!(rings.num_rings(), 1);
        let conf = Conformer2D::new(
            0,
            vec![
                [0.0, 0.0],
                [0.0, 1.0],
                [-1.0, 0.0],
                [1.0, 1.0],
                [-2.0, 1.0],
                [3.0, 3.0],
            ],
        );
        let mut map = WedgeAssignments::default();
        assert!(
            wedge_atropisomer_two_d_source(
                &mut topology,
                &rings,
                BondId::new(0),
                AtropisomerConformer::TwoD(&conf),
                &mut map
            )
            .unwrap()
        );
        assert_eq!(topology.bonds[4].direction(), BondDirection::None);
        assert_eq!(map.iter().count(), 1);
        assert!(map.get(BondId::new(3)).is_some());
        assert_eq!(topology.bonds[3].direction(), BondDirection::BeginWedge);
    }
    #[test]
    fn source634_projected_earlier_axial_orientation_matches_actual_mutable_source_state() {
        let topology = graph(&edges(), BondStereo::AtropCw);
        let conf = Conformer2D::new(
            0,
            vec![
                [0.0, 0.0],
                [0.0, 1.0],
                [-1.0, 0.0],
                [1.0, 1.0],
                [2.0; 2],
                [3.0; 2],
            ],
        );
        let mut actual = topology.clone();
        actual.bonds[0].set_endpoints(AtomId::new(1), AtomId::new(0));
        let mut actual_map = WedgeAssignments::default();
        assert!(
            wedge_atropisomer_two_d_source(
                &mut actual,
                &sparse(),
                BondId::new(0),
                AtropisomerConformer::TwoD(&conf),
                &mut actual_map
            )
            .unwrap()
        );
        assert!(actual_map.get(BondId::new(2)).is_some());
        let mut updates = BTreeMap::from([(
            BondId::new(0),
            AtropisomerWedgeUpdate {
                bond: BondId::new(0),
                begin: AtomId::new(1),
                end: AtomId::new(0),
                direction: BondDirection::None,
                atropisomer_bond: BondId::new(0),
            },
        )]);
        let occupied = BTreeSet::new();
        let mut writes = BTreeSet::new();
        let mut state = ProjectedAtropWedgeState {
            topology: &topology,
            occupied: &occupied,
            updates: &mut updates,
            map_writes: &mut writes,
        };
        wedge_one_two_d_source(
            &mut state,
            &sparse(),
            BondId::new(0),
            AtropisomerConformer::TwoD(&conf),
        )
        .unwrap();
        assert_eq!(writes, BTreeSet::from([BondId::new(2)]));
        let projected = updates.get(&BondId::new(2)).unwrap();
        let real = &actual.bonds[2];
        assert_eq!(
            (projected.begin, projected.end, projected.direction),
            (real.begin(), real.end(), real.direction())
        );
        // The no-conformer branch also reads the same actual prior orientation.
        let mut actual = topology.clone();
        actual.bonds[0].set_endpoints(AtomId::new(1), AtomId::new(0));
        let mut actual_map = WedgeAssignments::default();
        assert!(
            wedge_atropisomer_no_conformer_source(
                &mut actual,
                &sparse(),
                BondId::new(0),
                &mut actual_map
            )
            .unwrap()
        );
        assert!(actual_map.get(BondId::new(1)).is_some());
        updates.retain(|id, _| *id == BondId::new(0));
        writes.clear();
        let mut state = ProjectedAtropWedgeState {
            topology: &topology,
            occupied: &occupied,
            updates: &mut updates,
            map_writes: &mut writes,
        };
        wedge_one_no_conformer_source(&mut state, &sparse(), BondId::new(0)).unwrap();
        assert_eq!(writes, BTreeSet::from([BondId::new(1)]));
        let projected = updates.get(&BondId::new(1)).unwrap();
        let real = &actual.bonds[1];
        assert_eq!(
            (projected.begin, projected.end, projected.direction),
            (real.begin(), real.end(), real.direction())
        );
    }
}

#[cfg(test)]
mod source638_three_d_direction_tests {
    use super::*;
    #[test]
    fn source638_actual_z_delta_source_threshold_and_ieee_cells() {
        let bond = CurrentAtropBondEndpoints {
            id: BondId::new(0),
            begin: AtomId::new(0),
            end: AtomId::new(1),
        };
        for (begin_z, end_z, expected) in [
            (0.0, 1e-7, BondDirection::BeginDash),
            (
                0.0,
                f64::from_bits(1e-7f64.to_bits() + 1),
                BondDirection::BeginWedge,
            ),
            (
                0.0,
                f64::from_bits(1e-7f64.to_bits() - 1),
                BondDirection::BeginDash,
            ),
            (0.0, -1.0, BondDirection::BeginDash),
            (1.0, 2.0, BondDirection::BeginWedge),
            (0.0, -0.0, BondDirection::BeginDash),
            (0.0, f64::NAN, BondDirection::BeginDash),
            (f64::NAN, 1.0, BondDirection::BeginDash),
            (0.0, f64::INFINITY, BondDirection::BeginWedge),
            (f64::INFINITY, f64::INFINITY, BondDirection::BeginDash),
            (f64::NEG_INFINITY, 0.0, BondDirection::BeginWedge),
            (0.0, f64::NEG_INFINITY, BondDirection::BeginDash),
        ] {
            for is_3d in [false, true] {
                let conf = Conformer3D::new(
                    0,
                    vec![[100.0, 200.0, begin_z], [-100.0, -200.0, end_z]],
                    is_3d,
                );
                assert_eq!(
                    three_d_direction(&bond, AtropisomerConformer::ThreeD(&conf)),
                    expected,
                    "z0={begin_z}, z1={end_z}, flag={is_3d}"
                );
            }
        }
    }
    #[test]
    fn source638_current_source_endpoints_and_two_d_storage_markers() {
        let conf = Conformer3D::new(0, vec![[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]], true);
        let forward = CurrentAtropBondEndpoints {
            id: BondId::new(0),
            begin: AtomId::new(0),
            end: AtomId::new(1),
        };
        let reverse = CurrentAtropBondEndpoints {
            id: BondId::new(0),
            begin: AtomId::new(1),
            end: AtomId::new(0),
        };
        assert_eq!(
            three_d_direction(&forward, AtropisomerConformer::ThreeD(&conf)),
            BondDirection::BeginWedge
        );
        assert_eq!(
            three_d_direction(&reverse, AtropisomerConformer::ThreeD(&conf)),
            BondDirection::BeginDash
        );
        let conf = Conformer2D::new(0, vec![[100.0, 200.0], [-100.0, -200.0]]);
        assert_eq!(
            three_d_direction(&forward, AtropisomerConformer::TwoD(&conf)),
            BondDirection::BeginDash
        );
    }
}
#[cfg(test)]
mod source642_three_d_wedge_tests {
    use super::*;
    use crate::RingFindType;
    use cosmolkit_model::{Atom, AtomSpec, BondSpec};
    use cosmolkit_types::Element;
    fn graph(
        edges: &[(usize, usize, BondDirection, BondOrder)],
        stereo: BondStereo,
    ) -> TopologyBlock {
        let atoms = (0..6)
            .map(|id| Atom::from_spec(AtomId::new(id), AtomSpec::new(Element::C)))
            .collect();
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(id, &(begin, end, direction, order))| {
                let spec = BondSpec::new(AtomId::new(begin), AtomId::new(end), order)
                    .with_direction(direction);
                Bond::from_spec(
                    BondId::new(id),
                    if id == 0 {
                        spec.with_stereo(stereo)
                    } else {
                        spec
                    },
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
    }
    fn sparse() -> RingInfo {
        RingInfo::new(RingFindType::Sssr, 0, 0)
    }
    fn coordinates() -> Conformer2D {
        Conformer2D::new(
            0,
            vec![
                [0.0, 0.0],
                [0.0, 1.0],
                [-1.0, 0.0],
                [-1.0, 1.0],
                [1.0, 0.0],
                [1.0, 1.0],
            ],
        )
    }
    fn edges() -> [(usize, usize, BondDirection, BondOrder); 3] {
        [
            (0, 1, BondDirection::None, BondOrder::Single),
            (2, 0, BondDirection::None, BondOrder::Single),
            (3, 1, BondDirection::None, BondOrder::Single),
        ]
    }

    fn conf(z2: f64, z3: f64) -> Conformer3D {
        Conformer3D::new(
            0,
            vec![
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, z2],
                [0.0, 0.0, z3],
                [0.0, 0.0, -1.0],
                [0.0, 0.0, -1.0],
            ],
            true,
        )
    }
    fn apply(
        t: &mut TopologyBlock,
        r: &RingInfo,
        c: &Conformer3D,
        m: &mut WedgeAssignments,
    ) -> Result<bool, AtropisomerError> {
        wedge_atropisomer_three_d_source(t, r, BondId::new(0), AtropisomerConformer::ThreeD(c), m)
    }
    fn mark(t: &TopologyBlock, m: &mut WedgeAssignments, carrier: usize, axial: usize) {
        m.insert_atropisomer_source(AtropisomerWedgeUpdate {
            bond: BondId::new(carrier),
            begin: t.bonds[carrier].begin(),
            end: t.bonds[carrier].end(),
            direction: BondDirection::BeginDash,
            atropisomer_bond: BondId::new(axial),
        });
    }
    #[test]
    fn source642_existing_carrier_uses_carrier_eligibility_without_axis_geometry_or_ring_reads() {
        let mut t = graph(
            &[
                (0, 1, BondDirection::None, BondOrder::Double),
                (0, 2, BondDirection::BeginDash, BondOrder::Single),
                (1, 3, BondDirection::BeginWedge, BondOrder::Single),
            ],
            BondStereo::None,
        );
        let mut r = sparse();
        r.reset();
        let mut m = WedgeAssignments::default();
        let c = Conformer3D::new(
            0,
            vec![
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, -1.0],
                [0.0; 3],
                [0.0; 3],
            ],
            false,
        );
        assert!(apply(&mut t, &r, &c, &mut m).unwrap());
        assert_eq!(t.bonds[1].direction(), BondDirection::BeginWedge);
        assert_eq!(t.bonds[2].direction(), BondDirection::BeginDash);
        assert_eq!(m.iter().count(), 0);
    }
    #[test]
    fn source642_all_allowed_trials_are_oriented_even_when_not_chosen() {
        let mut t = graph(&edges(), BondStereo::None);
        let mut m = WedgeAssignments::default();
        assert!(apply(&mut t, &sparse(), &conf(-1.0, 1.0), &mut m).unwrap());
        assert_eq!(
            (t.bonds[1].begin(), t.bonds[1].end(), t.bonds[1].direction()),
            (AtomId::new(0), AtomId::new(2), BondDirection::None)
        );
        assert_eq!(
            (t.bonds[2].begin(), t.bonds[2].end(), t.bonds[2].direction()),
            (AtomId::new(1), AtomId::new(3), BondDirection::BeginWedge)
        );
        assert_eq!(m.iter().count(), 1);
        assert!(m.get(BondId::new(2)).is_some());
    }
    #[test]
    fn source642_later_direction_conflict_retains_actual_mutation_prefix() {
        let mut es = edges();
        es[2].2 = BondDirection::EndUpRight;
        let mut t = graph(&es, BondStereo::AtropCw);
        let mut m = WedgeAssignments::default();
        assert!(!apply(&mut t, &sparse(), &conf(1.0, -1.0), &mut m).unwrap());
        assert_eq!(
            (t.bonds[1].begin(), t.bonds[1].end(), t.bonds[1].direction()),
            (AtomId::new(0), AtomId::new(2), BondDirection::None)
        );
        assert_eq!(
            (t.bonds[2].begin(), t.bonds[2].end(), t.bonds[2].direction()),
            (AtomId::new(1), AtomId::new(3), BondDirection::EndUpRight)
        );
        assert_eq!(m.iter().count(), 0);
    }
    #[test]
    fn source642_source_checks_axial_map_key_and_can_overwrite_carrier_map_entry() {
        let mut t = graph(&edges(), BondStereo::AtropCcw);
        let mut m = WedgeAssignments::default();
        mark(&t, &mut m, 1, 2);
        assert!(apply(&mut t, &sparse(), &conf(1.0, -1.0), &mut m).unwrap());
        match m.get(BondId::new(1)).unwrap() {
            WedgeInfo::Atropisomer { update } => {
                assert_eq!(update.atropisomer_bond, BondId::new(0));
                assert_eq!(update.direction, BondDirection::BeginWedge);
            }
            _ => panic!("native atrop info required"),
        }
        let mut t = graph(&edges(), BondStereo::AtropCcw);
        let before = t.clone();
        let mut m = WedgeAssignments::default();
        mark(&t, &mut m, 0, 2);
        let before_m = m.clone();
        let mut r = sparse();
        r.reset();
        assert!(!apply(&mut t, &r, &conf(1.0, -1.0), &mut m).unwrap());
        assert_eq!(t, before);
        assert_eq!(m, before_m);
    }
    #[test]
    fn source642_uninitialized_ring_error_follows_first_actual_trial_orientation() {
        let mut t = graph(&edges(), BondStereo::AtropCw);
        let mut m = WedgeAssignments::default();
        let mut r = sparse();
        r.reset();
        assert!(matches!(
            apply(&mut t, &r, &conf(1.0, -1.0), &mut m),
            Err(AtropisomerError::SourcePrecondition {
                message: "RingInfo not initialized"
            })
        ));
        assert_eq!(
            (t.bonds[1].begin(), t.bonds[1].end()),
            (AtomId::new(0), AtomId::new(2))
        );
        assert_eq!(
            (t.bonds[2].begin(), t.bonds[2].end()),
            (AtomId::new(3), AtomId::new(1))
        );
        assert_eq!(m.iter().count(), 0);
        assert_eq!(t.bonds[1].direction(), BondDirection::None);
    }
    #[test]
    fn source642_single_preference_and_real_ring_cache_precede_wedge_preference() {
        let mut es = edges();
        es[2].3 = BondOrder::Aromatic;
        let mut t = graph(&es, BondStereo::None);
        let mut m = WedgeAssignments::default();
        assert!(apply(&mut t, &sparse(), &conf(-1.0, 1.0), &mut m).unwrap());
        assert_eq!(t.bonds[1].direction(), BondDirection::BeginDash);
        assert_eq!(t.bonds[2].begin(), AtomId::new(1));
        assert!(m.get(BondId::new(1)).is_some());
        let mut t = graph(
            &[
                (0, 1, BondDirection::None, BondOrder::Single),
                (2, 0, BondDirection::None, BondOrder::Single),
                (2, 4, BondDirection::None, BondOrder::Single),
                (4, 1, BondDirection::None, BondOrder::Single),
                (3, 1, BondDirection::None, BondOrder::Single),
            ],
            BondStereo::None,
        );
        let r = crate::find_sssr(&t, &crate::RingSearchParams::default()).unwrap();
        assert_eq!(r.num_rings(), 1);
        let mut m = WedgeAssignments::default();
        assert!(apply(&mut t, &r, &conf(-1.0, 1.0), &mut m).unwrap());
        assert!(m.get(BondId::new(1)).is_some());
        assert_eq!(t.bonds[1].direction(), BondDirection::BeginDash);
        assert_eq!(t.bonds[4].begin(), AtomId::new(1));
        assert_eq!(t.bonds[4].direction(), BondDirection::None);
    }
    #[test]
    fn source642_nan_marker_and_unknown_or_missing_guards_preserve_source_order() {
        let mut t = graph(&edges(), BondStereo::None);
        let mut m = WedgeAssignments::default();
        assert!(apply(&mut t, &sparse(), &conf(f64::NAN, f64::NAN), &mut m).unwrap());
        assert!(m.get(BondId::new(1)).is_some());
        assert_eq!(t.bonds[1].direction(), BondDirection::BeginDash);
        for missing in [false, true] {
            let es = if missing {
                vec![edges()[0], edges()[1]]
            } else {
                let mut es = edges().to_vec();
                es[2].2 = BondDirection::Unknown;
                es
            };
            let mut t = graph(&es, BondStereo::AtropCw);
            let before = t.clone();
            let mut m = WedgeAssignments::default();
            let mut r = sparse();
            r.reset();
            assert!(!apply(&mut t, &r, &conf(1.0, -1.0), &mut m).unwrap());
            assert_eq!(t, before);
            assert_eq!(m.iter().count(), 0);
        }
    }
    #[test]
    fn source642_projected_earlier_axial_orientation_and_failure_prefix_match_real_state() {
        for failure in [false, true] {
            let mut es = edges();
            if failure {
                es[2].2 = BondDirection::EndUpRight;
            }
            let original = graph(&es, BondStereo::None);
            let mut real = original.clone();
            real.bonds[0].set_endpoints(AtomId::new(1), AtomId::new(0));
            let mut map = WedgeAssignments::default();
            let c = conf(-1.0, 1.0);
            let actual = apply(&mut real, &sparse(), &c, &mut map).unwrap();
            let mut updates = BTreeMap::from([(
                BondId::new(0),
                AtropisomerWedgeUpdate {
                    bond: BondId::new(0),
                    begin: AtomId::new(1),
                    end: AtomId::new(0),
                    direction: BondDirection::None,
                    atropisomer_bond: BondId::new(0),
                },
            )]);
            let occupied = BTreeSet::new();
            let mut writes = BTreeSet::new();
            let mut state = ProjectedAtropWedgeState {
                topology: &original,
                occupied: &occupied,
                updates: &mut updates,
                map_writes: &mut writes,
            };
            let result = wedge_one_three_d_source(
                &mut state,
                &sparse(),
                BondId::new(0),
                AtropisomerConformer::ThreeD(&c),
            );
            assert_eq!(result.is_ok(), actual);
            for (id, bond) in real.bonds.iter().enumerate() {
                let id = BondId::new(id);
                assert_eq!(
                    (state.begin(id), state.end(id), state.direction(id)),
                    (bond.begin(), bond.end(), bond.direction())
                );
            }
            assert_eq!(writes, map.iter().map(|(id, _)| id).collect());
        }
    }
}

#[cfg(test)]
mod source654_total_degree_tests {
    use super::*;
    use cosmolkit_model::{Atom, AtomSpec, BondSpec};
    use cosmolkit_types::Element;
    fn graph(flag: bool, no_implicit: bool) -> TopologyBlock {
        let atoms = vec![
            Atom::from_spec(
                AtomId::new(0),
                AtomSpec::new(Element::C)
                    .with_explicit_hydrogens(255)
                    .with_implicit_hydrogen(flag)
                    .with_no_implicit(no_implicit),
            ),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::H)),
            Atom::from_spec(AtomId::new(2), AtomSpec::new(Element::DUMMY)),
        ];
        let bonds = vec![
            Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Zero),
            ),
            Bond::from_spec(
                BondId::new(1),
                BondSpec::new(AtomId::new(2), AtomId::new(0), BondOrder::Dative),
            ),
        ];
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
    }
    #[test]
    fn source654_actual_total_h_and_graph_degree_without_neighbor_double_count() {
        for flag in [false, true] {
            let t = graph(flag, false);
            let v = crate::ValenceAssignment {
                explicit_valence: vec![],
                implicit_hydrogens: vec![127],
            };
            assert_eq!(source_total_degree(&t, &v, AtomId::new(0)), Ok(384));
            assert_eq!(
                crate::hcount::total_hydrogen_count_from_validated(&t, &v, AtomId::new(0), true),
                Ok(383)
            );
        }
    }
    #[test]
    fn source654_source_cache_error_or_no_implicit_short_circuit() {
        for no_implicit in [false, true] {
            let t = graph(false, no_implicit);
            for rows in [vec![], vec![-1]] {
                let v = crate::ValenceAssignment {
                    explicit_valence: vec![],
                    implicit_hydrogens: rows,
                };
                if no_implicit {
                    assert_eq!(source_total_degree(&t, &v, AtomId::new(0)), Ok(257));
                } else {
                    assert_eq!(
                        source_total_degree(&t, &v, AtomId::new(0)),
                        Err(crate::ValenceError::ImplicitValenceCacheNotInitialized {
                            atom: AtomId::new(0)
                        })
                    );
                }
            }
        }
    }
}
#[cfg(test)]
mod source658_whole_atrop_wedge_tests {
    use super::*;
    use crate::RingFindType;
    use cosmolkit_model::{Atom, AtomSpec, BondSpec};
    use cosmolkit_types::Element;
    fn graph(
        edges: &[(usize, usize, BondDirection, BondOrder)],
        stereo: BondStereo,
    ) -> TopologyBlock {
        let atoms = (0..6)
            .map(|id| {
                Atom::from_spec(
                    AtomId::new(id),
                    AtomSpec::new(Element::C).with_no_implicit(true),
                )
            })
            .collect();
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(id, &(begin, end, direction, order))| {
                let spec = BondSpec::new(AtomId::new(begin), AtomId::new(end), order)
                    .with_direction(direction);
                Bond::from_spec(
                    BondId::new(id),
                    if id == 0 {
                        spec.with_stereo(stereo)
                    } else {
                        spec
                    },
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
    }
    fn sparse() -> RingInfo {
        RingInfo::new(RingFindType::Sssr, 0, 0)
    }
    fn coordinates() -> Conformer2D {
        Conformer2D::new(
            0,
            vec![
                [0.0, 0.0],
                [0.0, 1.0],
                [-1.0, 0.0],
                [-1.0, 1.0],
                [1.0, 0.0],
                [1.0, 1.0],
            ],
        )
    }
    fn edges() -> [(usize, usize, BondDirection, BondOrder); 3] {
        [
            (0, 1, BondDirection::None, BondOrder::Single),
            (2, 0, BondDirection::None, BondOrder::Single),
            (3, 1, BondDirection::None, BondOrder::Single),
        ]
    }

    fn weak() -> RingInfo {
        RingInfo::new(RingFindType::Fast, 6, 3)
    }
    fn apply(
        t: &mut TopologyBlock,
        r: &mut RingInfo,
        c: Option<AtropisomerConformer<'_>>,
        m: &mut WedgeAssignments,
    ) -> Result<(), AtropisomerError> {
        wedge_bonds_from_atropisomers_source(t, r, None, c, m)
    }
    fn cache(t: &mut TopologyBlock, id: usize, value: i8) {
        t.atoms[id].set_no_implicit(false);
        t.atoms[id].set_source_valence_facts(SourceAtomValenceFacts {
            explicit_valence: -1,
            implicit_valence: value,
        });
    }
    #[test]
    fn source658_promotes_real_sparse_sssr_before_even_empty_candidate_iteration() {
        let mut t = graph(&edges(), BondStereo::None);
        let before = t.clone();
        let mut r = weak();
        let mut m = WedgeAssignments::default();
        apply(&mut t, &mut r, None, &mut m).unwrap();
        assert!(r.is_sssr_or_better());
        assert_eq!(r.atom_row_count(), 0);
        assert_eq!(r.bond_row_count(), 0);
        assert_eq!(t, before);
        assert_eq!(m.iter().count(), 0);
        let mut t = graph(&edges(), BondStereo::AtropCcw);
        apply(&mut t, &mut r, None, &mut m).unwrap();
        assert!(m.get(BondId::new(1)).is_some());
    }
    #[test]
    fn source658_trusted_sparse_cache_and_mutable_projected_dispatch_agree() {
        let original = graph(&edges(), BondStereo::AtropCcw);
        let c2 = coordinates();
        let c3 = Conformer3D::new(
            0,
            vec![
                [0.0; 3],
                [0.0; 3],
                [0.0, 0.0, -1.0],
                [0.0, 0.0, 1.0],
                [0.0; 3],
                [0.0; 3],
            ],
            true,
        );
        for c in [
            None,
            Some(AtropisomerConformer::TwoD(&c2)),
            Some(AtropisomerConformer::ThreeD(&c3)),
        ] {
            let mut real = original.clone();
            let mut r = sparse();
            let mut m = WedgeAssignments::default();
            apply(&mut real, &mut r, c, &mut m).unwrap();
            let mut rp = sparse();
            let projection = wedge_bonds_from_atropisomers_projected_source(
                &original,
                &mut rp,
                None,
                c,
                &BTreeSet::new(),
            )
            .unwrap();
            assert_eq!(
                projection.source_map_writes,
                m.iter().map(|(id, _)| id).collect::<Vec<_>>()
            );
            assert_eq!(projection.diagnostics, m.diagnostics());
            for update in &projection.bond_updates {
                let actual = &real.bonds[update.bond.index()];
                assert_eq!(
                    (update.begin, update.end, update.direction),
                    (actual.begin(), actual.end(), actual.direction())
                );
            }
            assert_eq!(r, rp);
        }
    }
    #[test]
    fn source658_actual_cache_counts_can_skip_axis_without_legacy_flag_guess() {
        for flag in [false, true] {
            let mut t = graph(&edges(), BondStereo::AtropCcw);
            cache(&mut t, 0, 2);
            t.atoms[0].set_implicit_hydrogen(flag);
            let before = t.clone();
            let mut r = sparse();
            let mut m = WedgeAssignments::default();
            apply(&mut t, &mut r, None, &mut m).unwrap();
            assert_eq!(t, before);
            assert_eq!(m.iter().count(), 0);
            cache(&mut t, 0, 1);
            apply(&mut t, &mut r, None, &mut m).unwrap();
            assert!(m.get(BondId::new(1)).is_some());
        }
    }
    #[test]
    fn source658_condition_order_reads_end_after_begin_lower_bound_before_begin_upper_bound() {
        let mut t = graph(&edges(), BondStereo::AtropCcw);
        cache(&mut t, 0, 2);
        cache(&mut t, 1, -1);
        let before = t.clone();
        let mut r = sparse();
        let mut m = WedgeAssignments::default();
        assert!(
            matches!(apply(&mut t,&mut r,None,&mut m),Err(AtropisomerError::Valence(crate::ValenceError::ImplicitValenceCacheNotInitialized {atom})) if atom==AtomId::new(1))
        );
        assert_eq!(t, before);
        assert_eq!(m.iter().count(), 0);
        let mut t = graph(&[edges()[0], edges()[2]], BondStereo::AtropCcw);
        cache(&mut t, 0, 0);
        cache(&mut t, 1, -1);
        let before = t.clone();
        apply(&mut t, &mut r, None, &mut m).unwrap();
        assert_eq!(t, before);
    }
    #[test]
    fn source658_ring_property_error_preserves_native_cache_prefix_before_any_candidate_getter() {
        let mut t = graph(&edges(), BondStereo::AtropCcw);
        cache(&mut t, 0, -1);
        let before = t.clone();
        let mut r = weak();
        let mut m = WedgeAssignments::default();
        let mut props = cosmolkit_model::MoleculeProperties::default();
        props.set_prop("__computedProps", false).unwrap();
        assert!(matches!(
            wedge_bonds_from_atropisomers_source(&mut t, &mut r, Some(&mut props), None, &mut m),
            Err(AtropisomerError::RingFinding(
                crate::RingFindingError::MoleculeProperty(_)
            ))
        ));
        assert!(r.is_sssr_or_better());
        assert_eq!(r.atom_row_count(), 0);
        assert_eq!(t, before);
        assert_eq!(m.iter().count(), 0);
        assert_eq!(
            props.prop("__computedProps"),
            Some(&cosmolkit_model::PropertyValue::Bool(false))
        );
    }
    #[test]
    fn source658_ignored_false_keeps_three_d_orientation_prefix_and_source_diagnostic() {
        let mut es = edges();
        es[2].2 = BondDirection::EndUpRight;
        let mut t = graph(&es, BondStereo::AtropCcw);
        let c = Conformer3D::new(0, vec![[0.0; 3]; 6], true);
        let mut r = sparse();
        let mut m = WedgeAssignments::default();
        apply(
            &mut t,
            &mut r,
            Some(AtropisomerConformer::ThreeD(&c)),
            &mut m,
        )
        .unwrap();
        assert_eq!(t.bonds[1].begin(), AtomId::new(0));
        assert_eq!(t.bonds[2].begin(), AtomId::new(1));
        assert_eq!(m.iter().count(), 0);
        assert_eq!(
            m.diagnostics(),
            &[AtropisomerDiagnostic {
                bond: BondId::new(0),
                kind: AtropisomerRejectionKind::DirectionConflict
            }]
        );
    }
    #[test]
    fn source658_checked_query_keeps_exact_row_guard_while_source_adapter_accepts_sparse() {
        let t = graph(&edges(), BondStereo::AtropCcw);
        assert!(matches!(
            wedge_bonds_from_atropisomers(&t, &sparse(), None, &BTreeSet::new()),
            Err(AtropisomerError::RingAtomRowCount {
                actual: 0,
                expected: 6
            })
        ));
        let mut r = sparse();
        let a = wedge_bonds_from_atropisomers_projected_source(
            &t,
            &mut r,
            None,
            None,
            &BTreeSet::new(),
        )
        .unwrap();
        assert_eq!(a.source_map_writes, vec![BondId::new(1)]);
        assert_eq!(r.atom_row_count(), 0);
    }
    #[test]
    fn source658_later_actual_cache_error_preserves_earlier_native_graph_and_map_writes() {
        let mut t = graph(
            &[
                (0, 1, BondDirection::None, BondOrder::Single),
                (2, 0, BondDirection::None, BondOrder::Single),
                (3, 1, BondDirection::None, BondOrder::Single),
                (3, 4, BondDirection::None, BondOrder::Single),
                (5, 4, BondDirection::None, BondOrder::Single),
            ],
            BondStereo::AtropCcw,
        );
        t.bonds[3].set_stereo(BondStereo::AtropCcw).unwrap();
        cache(&mut t, 3, -1);
        let mut r = sparse();
        let mut m = WedgeAssignments::default();
        assert!(
            matches!(apply(&mut t,&mut r,None,&mut m),Err(AtropisomerError::Valence(crate::ValenceError::ImplicitValenceCacheNotInitialized {atom})) if atom==AtomId::new(3))
        );
        assert_eq!(
            (t.bonds[1].begin(), t.bonds[1].end(), t.bonds[1].direction()),
            (AtomId::new(0), AtomId::new(2), BondDirection::BeginWedge)
        );
        assert!(m.get(BondId::new(1)).is_some());
        assert_eq!(m.iter().count(), 1);
        assert_eq!(t.bonds[4].begin(), AtomId::new(5));
        assert_eq!(t.bonds[4].direction(), BondDirection::None);
    }
}

#[cfg(test)]
mod complete_stereo_group_source_tests {
    use super::*;
    use cosmolkit_model::{AdjacencyList, Atom, AtomSpec, BondSpec, StereoGroupKind};
    use cosmolkit_types::Element;

    fn graph(count: usize, edges: &[(usize, usize, BondDirection)]) -> TopologyBlock {
        let atoms = (0..count)
            .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
            .collect();
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(i, (a, b, direction))| {
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(
                        AtomId::new(*a),
                        AtomId::new(*b),
                        cosmolkit_types::BondOrder::Single,
                    )
                    .with_direction(*direction),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
    }
    fn group(atoms: &[usize], bonds: &[usize]) -> StereoGroup {
        StereoGroup::new(
            StereoGroupKind::Absolute,
            atoms.iter().copied().map(AtomId::new).collect(),
            bonds.iter().copied().map(BondId::new).collect(),
        )
        .expect("valid distinct stereo members")
    }
    fn indices(output: &[AtomId]) -> Vec<usize> {
        output.iter().map(|id| id.index()).collect()
    }

    #[test]
    fn duplicate_members_are_rejected_and_reused_output_keeps_order_and_capacity() {
        // The original illegal fixture is now rejected at the source constructor.
        assert_eq!(
            StereoGroup::new(
                StereoGroupKind::Absolute,
                [2, 0, 2].into_iter().map(AtomId::new).collect(),
                [0, 0].into_iter().map(BondId::new).collect(),
            ),
            Err(cosmolkit_model::StereoGroupError::DuplicateAtom),
        );
        // The original illegal fixture is now rejected at the source constructor.
        assert_eq!(
            StereoGroup::new(
                StereoGroupKind::Absolute,
                [2, 0].into_iter().map(AtomId::new).collect(),
                [0, 0].into_iter().map(BondId::new).collect(),
            ),
            Err(cosmolkit_model::StereoGroupError::DuplicateBond),
        );
        let topology = graph(
            3,
            &[
                (0, 1, BondDirection::None),
                (1, 2, BondDirection::BeginWedge),
            ],
        );
        let mut output = Vec::with_capacity(40);
        output.push(AtomId::new(999));
        let pointer = output.as_ptr();
        collect_stereo_group_atom_ids_source(
            &topology,
            &group(&[2, 0], &[0]),
            &mut output,
            &WedgeAssignments::default(),
        )
        .unwrap();
        assert_eq!(indices(&output), vec![2, 0, 1]);
        assert_eq!(output.as_ptr(), pointer);
        collect_stereo_group_atom_ids_source(
            &topology,
            &group(&[], &[]),
            &mut output,
            &WedgeAssignments::default(),
        )
        .unwrap();
        assert!(output.is_empty());
        assert_eq!(output.as_ptr(), pointer);
    }

    #[test]
    fn group_atom_metadata_is_copied_without_a_source_graph_dereference() {
        // The original illegal fixture is now rejected at the source constructor.
        assert_eq!(
            StereoGroup::new(
                StereoGroupKind::Absolute,
                [77, 0, 77].into_iter().map(AtomId::new).collect(),
                Vec::new(),
            ),
            Err(cosmolkit_model::StereoGroupError::DuplicateAtom),
        );
        let mut output = vec![AtomId::new(500)];
        collect_stereo_group_atom_ids_source(
            &TopologyBlock::default(),
            &group(&[77, 0], &[]),
            &mut output,
            &WedgeAssignments::default(),
        )
        .unwrap();
        assert_eq!(indices(&output), vec![77, 0]);
    }

    #[test]
    fn begin_then_end_order_and_neighbor_direction_ignore_carrier_begin_orientation() {
        // The original illegal fixture is now rejected at the source constructor.
        assert_eq!(
            StereoGroup::new(
                StereoGroupKind::Absolute,
                [2, 2].into_iter().map(AtomId::new).collect(),
                [0].into_iter().map(BondId::new).collect(),
            ),
            Err(cosmolkit_model::StereoGroupError::DuplicateAtom),
        );
        let topology = graph(
            4,
            &[
                (1, 2, BondDirection::None),
                (0, 1, BondDirection::BeginWedge),
                (3, 2, BondDirection::BeginDash),
            ],
        );
        let mut output = Vec::new();
        collect_stereo_group_atom_ids_source(
            &topology,
            &group(&[], &[0]),
            &mut output,
            &WedgeAssignments::default(),
        )
        .unwrap();
        assert_eq!(indices(&output), vec![1, 2]);
        collect_stereo_group_atom_ids_source(
            &topology,
            &group(&[2], &[0]),
            &mut output,
            &WedgeAssignments::default(),
        )
        .unwrap();
        assert_eq!(indices(&output), vec![2, 1]);
    }

    #[test]
    fn only_wedge_dash_or_actual_atrop_info_type_marks_endpoint() {
        for direction in [
            BondDirection::None,
            BondDirection::BeginWedge,
            BondDirection::BeginDash,
            BondDirection::EndDownRight,
            BondDirection::EndUpRight,
            BondDirection::EitherDouble,
            BondDirection::Unknown,
        ] {
            let topology = graph(3, &[(0, 1, BondDirection::None), (1, 2, direction)]);
            let mut output = Vec::new();
            collect_stereo_group_atom_ids_source(
                &topology,
                &group(&[], &[0]),
                &mut output,
                &WedgeAssignments::default(),
            )
            .unwrap();
            assert_eq!(
                indices(&output),
                if matches!(
                    direction,
                    BondDirection::BeginWedge | BondDirection::BeginDash
                ) {
                    vec![1]
                } else {
                    vec![]
                }
            );
        }
        let topology = graph(
            3,
            &[(0, 1, BondDirection::None), (1, 2, BondDirection::None)],
        );
        let update = AtropisomerWedgeUpdate {
            bond: BondId::new(1),
            begin: AtomId::new(1),
            end: AtomId::new(2),
            direction: BondDirection::Unknown,
            atropisomer_bond: BondId::new(99),
        };
        let map = WedgeAssignments::from_atropisomer_wedge_parts_source(&[update], &[update.bond]);
        let mut output = Vec::new();
        collect_stereo_group_atom_ids_source(&topology, &group(&[], &[0]), &mut output, &map)
            .unwrap();
        assert_eq!(indices(&output), vec![1]);
    }

    #[test]
    fn later_bond_error_preserves_cleared_metadata_and_earlier_appended_prefix() {
        // The original illegal fixture is now rejected at the source constructor.
        assert_eq!(
            StereoGroup::new(
                StereoGroupKind::Absolute,
                [2, 2].into_iter().map(AtomId::new).collect(),
                [0, 99].into_iter().map(BondId::new).collect(),
            ),
            Err(cosmolkit_model::StereoGroupError::DuplicateAtom),
        );
        let topology = graph(
            3,
            &[
                (0, 1, BondDirection::None),
                (0, 2, BondDirection::BeginWedge),
            ],
        );
        let mut output = vec![AtomId::new(88)];
        assert_eq!(
            collect_stereo_group_atom_ids_source(
                &topology,
                &group(&[2], &[0, 99]),
                &mut output,
                &WedgeAssignments::default()
            ),
            Err(AtropisomerError::StereoGroupBondOutOfRange {
                bond: BondId::new(99),
                bond_count: 2
            })
        );
        assert_eq!(indices(&output), vec![2, 0]);
    }

    #[test]
    fn both_endpoint_getters_resolve_before_visiting_first_endpoint_adjacency() {
        let mut topology = graph(
            4,
            &[
                (0, 1, BondDirection::None),
                (0, 2, BondDirection::BeginWedge),
            ],
        );
        topology.bonds[0] = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(
                AtomId::new(0),
                AtomId::new(99),
                cosmolkit_types::BondOrder::Single,
            ),
        );
        let mut output = vec![AtomId::new(88)];
        assert_eq!(
            collect_stereo_group_atom_ids_source(
                &topology,
                &group(&[3], &[0]),
                &mut output,
                &WedgeAssignments::default()
            ),
            Err(AtropisomerError::StereoGroupAtomOutOfRange {
                atom: AtomId::new(99),
                atom_count: 4
            })
        );
        assert_eq!(indices(&output), vec![3]);
    }

    #[test]
    fn missing_actual_csr_row_is_not_silently_zero_degree_and_unused_rows_are_not_read() {
        // The original illegal fixture is now rejected at the source constructor.
        assert_eq!(
            StereoGroup::new(
                StereoGroupKind::Absolute,
                [2, 2].into_iter().map(AtomId::new).collect(),
                Vec::new(),
            ),
            Err(cosmolkit_model::StereoGroupError::DuplicateAtom),
        );
        let isolated = graph(2, &[]);
        assert_eq!(isolated.adjacency.try_neighbors_of(0), Some([].as_slice()));
        assert_eq!(isolated.adjacency.try_neighbors_of(2), None);
        assert_eq!(isolated.adjacency.try_neighbors_of(usize::MAX), None);
        let mut topology = graph(3, &[(0, 1, BondDirection::None)]);
        topology.adjacency = AdjacencyList::default();
        let mut output = vec![AtomId::new(88)];
        assert_eq!(
            collect_stereo_group_atom_ids_source(
                &topology,
                &group(&[2], &[0]),
                &mut output,
                &WedgeAssignments::default()
            ),
            Err(AtropisomerError::SourcePrecondition {
                message: "source stereo-group atom adjacency row is absent"
            })
        );
        assert_eq!(indices(&output), vec![2]);
        collect_stereo_group_atom_ids_source(
            &topology,
            &group(&[2], &[]),
            &mut output,
            &WedgeAssignments::default(),
        )
        .unwrap();
        assert_eq!(indices(&output), vec![2]);
    }

    #[test]
    fn legacy_detached_update_adapter_requires_real_map_write_and_ignores_axial_id_for_type_test() {
        let topology = graph(
            3,
            &[(0, 1, BondDirection::None), (1, 2, BondDirection::None)],
        );
        let update = AtropisomerWedgeUpdate {
            bond: BondId::new(1),
            begin: AtomId::new(1),
            end: AtomId::new(2),
            direction: BondDirection::None,
            atropisomer_bond: BondId::new(99),
        };
        let mut assignment = AtropisomerWedgeAssignment {
            source_map_writes: vec![],
            bond_updates: vec![update],
            diagnostics: vec![],
        };
        assert!(
            stereo_group_atom_ids(&topology, &group(&[], &[0]), &assignment)
                .unwrap()
                .is_empty()
        );
        assignment.source_map_writes.push(update.bond);
        assert_eq!(
            indices(&stereo_group_atom_ids(&topology, &group(&[], &[0]), &assignment).unwrap()),
            vec![1]
        );
    }
}

#[cfg(test)]
mod query_stereo_group_source_tests {
    use super::*;
    use cosmolkit_model::{
        AtomQueryPredicate, BondSpec, QueryAtom, QueryAtomIdentity, QueryBond, QueryGraph,
        QueryNode, StereoGroupKind,
    };
    fn graph(nums: &[u8], edges: &[(usize, usize, BondDirection)]) -> QueryGraph {
        QueryGraph::from_parts(
            nums.iter()
                .enumerate()
                .map(|(i, &n)| {
                    QueryAtom::from_identity_parts(
                        AtomId::new(i),
                        QueryAtomIdentity::from_atomic_number(n),
                        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(n)),
                    )
                })
                .collect(),
            edges
                .iter()
                .enumerate()
                .map(|(i, &(a, b, d))| {
                    QueryBond::new(
                        BondId::new(i),
                        BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single)
                            .with_direction(d),
                    )
                })
                .collect(),
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    fn group(atoms: &[usize], bonds: &[usize]) -> StereoGroup {
        StereoGroup::new(
            StereoGroupKind::Absolute,
            atoms.iter().copied().map(AtomId::new).collect(),
            bonds.iter().copied().map(BondId::new).collect(),
        )
        .expect("valid distinct stereo members")
    }
    fn collect(q: &QueryGraph, g: &StereoGroup, w: &WedgeAssignments) -> Vec<AtomId> {
        let mut out = vec![AtomId::new(999)];
        collect_query_stereo_group_atom_ids_source(q, g, &mut out, w).unwrap();
        out
    }
    #[test]
    fn query_collector_rejects_duplicate_members_and_keeps_reused_output_order() {
        // The original illegal fixture is now rejected at the source constructor.
        assert_eq!(
            StereoGroup::new(
                StereoGroupKind::Absolute,
                [2, 0, 2].into_iter().map(AtomId::new).collect(),
                [0, 0].into_iter().map(BondId::new).collect(),
            ),
            Err(cosmolkit_model::StereoGroupError::DuplicateAtom),
        );
        // The original illegal fixture is now rejected at the source constructor.
        assert_eq!(
            StereoGroup::new(
                StereoGroupKind::Absolute,
                [2, 0].into_iter().map(AtomId::new).collect(),
                [0, 0].into_iter().map(BondId::new).collect(),
            ),
            Err(cosmolkit_model::StereoGroupError::DuplicateBond),
        );
        let q = graph(
            &[6, 6, 6],
            &[
                (0, 1, BondDirection::None),
                (1, 2, BondDirection::BeginWedge),
            ],
        );
        let mut out = Vec::with_capacity(40);
        out.push(AtomId::new(999));
        let pointer = out.as_ptr();
        collect_query_stereo_group_atom_ids_source(
            &q,
            &group(&[2, 0], &[0]),
            &mut out,
            &Default::default(),
        )
        .unwrap();
        assert_eq!(out, vec![AtomId::new(2), AtomId::new(0), AtomId::new(1)]);
        assert_eq!(pointer, out.as_ptr());
    }
    #[test]
    fn actual_begin_and_end_members_are_appended_in_source_order() {
        let q = graph(
            &[6; 4],
            &[
                (0, 1, BondDirection::None),
                (2, 0, BondDirection::BeginWedge),
                (3, 1, BondDirection::BeginDash),
            ],
        );
        assert_eq!(
            collect(&q, &group(&[], &[0]), &Default::default()),
            vec![AtomId::new(0), AtomId::new(1)]
        );
    }
    #[test]
    fn source_group_bond_itself_never_marks_an_endpoint() {
        let q = graph(&[6; 2], &[(0, 1, BondDirection::BeginWedge)]);
        assert!(collect(&q, &group(&[], &[0]), &Default::default()).is_empty());
    }
    #[test]
    fn non_wedge_direction_does_not_mark_endpoint() {
        for direction in [
            BondDirection::None,
            BondDirection::Unknown,
            BondDirection::EndUpRight,
            BondDirection::EndDownRight,
        ] {
            let q = graph(&[6; 3], &[(0, 1, BondDirection::None), (1, 2, direction)]);
            assert!(collect(&q, &group(&[], &[0]), &Default::default()).is_empty());
        }
    }
    #[test]
    fn actual_atrop_info_marks_endpoint_without_direction_guessing() {
        let q = graph(
            &[6; 3],
            &[(0, 1, BondDirection::None), (1, 2, BondDirection::None)],
        );
        let mut wedge = WedgeAssignments::default();
        wedge.insert_atropisomer_source(AtropisomerWedgeUpdate {
            bond: BondId::new(1),
            begin: AtomId::new(2),
            end: AtomId::new(1),
            direction: BondDirection::None,
            atropisomer_bond: BondId::new(999),
        });
        assert_eq!(collect(&q, &group(&[], &[0]), &wedge), vec![AtomId::new(1)]);
    }
    #[test]
    fn source_ids_do_not_require_element_conversion_or_query_rewriting() {
        // The original illegal fixture is now rejected at the source constructor.
        assert_eq!(
            StereoGroup::new(
                StereoGroupKind::Absolute,
                [2, 2].into_iter().map(AtomId::new).collect(),
                [0].into_iter().map(BondId::new).collect(),
            ),
            Err(cosmolkit_model::StereoGroupError::DuplicateAtom),
        );
        let q = graph(
            &[255, 0, 119],
            &[
                (0, 1, BondDirection::None),
                (0, 2, BondDirection::BeginWedge),
            ],
        );
        let before = q.clone();
        assert_eq!(
            collect(&q, &group(&[2], &[0]), &Default::default()),
            vec![AtomId::new(2), AtomId::new(0)]
        );
        assert_eq!(q, before);
    }
    #[test]
    fn later_group_bond_failure_preserves_native_output_prefix() {
        // The original illegal fixture is now rejected at the source constructor.
        assert_eq!(
            StereoGroup::new(
                StereoGroupKind::Absolute,
                [2, 2].into_iter().map(AtomId::new).collect(),
                [0, 99].into_iter().map(BondId::new).collect(),
            ),
            Err(cosmolkit_model::StereoGroupError::DuplicateAtom),
        );
        let q = graph(
            &[6; 3],
            &[
                (0, 1, BondDirection::None),
                (1, 2, BondDirection::BeginWedge),
            ],
        );
        let mut out = vec![AtomId::new(999)];
        let e = collect_query_stereo_group_atom_ids_source(
            &q,
            &group(&[2], &[0, 99]),
            &mut out,
            &Default::default(),
        )
        .unwrap_err();
        assert!(
            matches!(e,AtropisomerError::StereoGroupBondOutOfRange{bond,bond_count:2}if bond==BondId::new(99))
        );
        assert_eq!(out, vec![AtomId::new(2), AtomId::new(1)]);
    }
}

#[cfg(test)]
mod cx_query_bond_config_dependency_source_tests {
    use super::*;

    use cosmolkit_model::{
        AtomQueryPredicate, BondSpec, Conformer2D, QueryAtom, QueryAtomIdentity, QueryBond,
        QueryGraph, QueryNode,
    };
    fn graph(n: usize, edges: &[(usize, usize)]) -> QueryGraph {
        QueryGraph::from_parts(
            (0..n)
                .map(|i| {
                    QueryAtom::from_identity_parts(
                        AtomId::new(i),
                        QueryAtomIdentity::from_atomic_number(if i % 2 == 0 { 0 } else { 119 }),
                        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                    )
                })
                .collect(),
            edges
                .iter()
                .enumerate()
                .map(|(i, &(a, b))| {
                    QueryBond::new(
                        BondId::new(i),
                        BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
                    )
                })
                .collect(),
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }

    #[test]
    fn query_carriers_sort_exactly_two_by_other_atom_without_identity_coercion() {
        let q = graph(6, &[(0, 1), (0, 3), (0, 2), (1, 5), (1, 4)]);
        let before = q.clone();
        let ends = query_atropisomer_carriers_source(&q, BondId::new(0))
            .unwrap()
            .unwrap();
        assert_eq!(ends[0].focus(), AtomId::new(0));
        assert_eq!(ends[0].carrier_bonds(), &[BondId::new(2), BondId::new(1)]);
        assert_eq!(ends[1].carrier_bonds(), &[BondId::new(4), BondId::new(3)]);
        assert_eq!(q, before);
    }
    #[test]
    fn query_three_carriers_keep_source_incident_order_and_empty_side_returns_false() {
        let q = graph(6, &[(0, 1), (0, 4), (0, 3), (0, 2), (1, 5)]);
        let ends = query_atropisomer_carriers_source(&q, BondId::new(0))
            .unwrap()
            .unwrap();
        assert_eq!(
            ends[0].carrier_bonds(),
            &[BondId::new(1), BondId::new(2), BondId::new(3)]
        );
        assert!(
            query_atropisomer_carriers_source(&graph(3, &[(0, 1), (0, 2)]), BondId::new(0))
                .unwrap()
                .is_none()
        );
    }
    #[test]
    fn query_source_output_appends_before_false_and_sets_both_focus_atoms() {
        let q = graph(3, &[(0, 1), (0, 2)]);
        let mut ends = [
            AtropEnd {
                atom: AtomId::new(2),
                bonds: vec![],
            },
            AtropEnd {
                atom: AtomId::new(2),
                bonds: vec![],
            },
        ];
        assert!(!atropisomer_ends(&q, q.bonds()[0].bond(), &mut ends).unwrap());
        assert_eq!(ends[0].atom, AtomId::new(0));
        assert_eq!(ends[1].atom, AtomId::new(1));
        assert_eq!(ends[0].bonds, vec![BondId::new(1)]);
        assert!(ends[1].bonds.is_empty());
    }
    #[test]
    fn query_carrier_range_error_is_not_empty_carrier_fallback() {
        assert!(matches!(
            query_atropisomer_carriers_source(&graph(2, &[(0, 1)]), BondId::new(1)),
            Err(AtropisomerError::AxialBondOutOfRange { bond_count: 1, .. })
        ));
    }
}

/// Query carrier entry to the same source algorithm; actual cache and writes.
#[doc(hidden)]
pub fn wedge_query_bonds_from_atropisomers_source(
    query: &mut cosmolkit_model::QueryGraph,
    rings: &mut RingInfo,
    properties: Option<&mut cosmolkit_model::MoleculeProperties>,
    conformer: Option<AtropisomerConformer<'_>>,
    wedges: &mut WedgeAssignments,
) -> Result<(), AtropisomerError> {
    wedge_bonds_from_atropisomers_graph_source(query, rings, properties, conformer, wedges)
}

#[cfg(test)]
mod source_query_atrop_warning_tests {
    use super::*;
    use cosmolkit_model::{
        AtomQueryPredicate, BondSpec, QueryAtom, QueryAtomIdentity, QueryBond, QueryGraph,
        QueryNode,
    };
    struct State {
        q: QueryGraph,
        wedges: WedgeAssignments,
        warnings: Vec<String>,
    }
    impl AtropWedgeState for State {
        type Graph = QueryGraph;
        fn topology(&self) -> &QueryGraph {
            &self.q
        }
        fn begin(&self, id: BondId) -> AtomId {
            self.q.bonds()[id.index()].begin()
        }
        fn end(&self, id: BondId) -> AtomId {
            self.q.bonds()[id.index()].end()
        }
        fn direction(&self, id: BondId) -> BondDirection {
            self.q.bonds()[id.index()].bond().direction()
        }
        fn occupied(&self, id: BondId) -> bool {
            self.wedges.get(id).is_some()
        }
        fn set_direction(&mut self, id: BondId, dir: BondDirection, axial: BondId) {
            MutableAtropWedgeState {
                topology: &mut self.q,
                wedges: &mut self.wedges,
            }
            .set_direction(id, dir, axial);
        }
        fn orient(&mut self, id: BondId, begin: AtomId, axial: BondId) {
            MutableAtropWedgeState {
                topology: &mut self.q,
                wedges: &mut self.wedges,
            }
            .orient(id, begin, axial);
        }
        fn insert_wedge(&mut self, id: BondId, axial: BondId) {
            MutableAtropWedgeState {
                topology: &mut self.q,
                wedges: &mut self.wedges,
            }
            .insert_wedge(id, axial);
        }
        fn source_warning(&mut self, prefix: &'static str, axial: BondId) {
            self.warnings.push(format!(
                "{prefix} {} {}",
                self.begin(axial).index(),
                self.end(axial).index()
            ));
        }
    }
    fn state() -> State {
        let mut q = QueryGraph::from_parts(
            (0..6)
                .map(|i| {
                    QueryAtom::from_identity_parts(
                        AtomId::new(i),
                        QueryAtomIdentity::from_atomic_number(0),
                        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                    )
                })
                .collect(),
            [(0, 1), (2, 0), (0, 3), (4, 1), (1, 5)]
                .into_iter()
                .enumerate()
                .map(|(i, (a, b))| {
                    QueryBond::new(
                        BondId::new(i),
                        BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
                    )
                })
                .collect(),
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap();
        q.bonds_mut()[0]
            .bond_mut()
            .set_stereo(BondStereo::AtropCw)
            .unwrap();
        State {
            q,
            wedges: WedgeAssignments::default(),
            warnings: vec![],
        }
    }
    fn rings() -> RingInfo {
        RingInfo::new(crate::RingFindType::Sssr, 0, 0)
    }
    #[test]
    fn missing_eligible_carrier_emits_exact_warning_before_false() {
        let mut s = state();
        for b in &mut s.q.bonds_mut()[1..] {
            b.bond_mut().set_order(BondOrder::Double);
        }
        assert!(matches!(
            wedge_one_no_conformer_source(&mut s, &rings(), BondId::new(0)),
            Err(PerceptionError::Rejected(
                AtropisomerRejectionKind::NoUsableWedgeBond
            ))
        ));
        assert_eq!(
            s.warnings,
            ["Failed to find a good bond to set as UP or DOWN for an atropisomer - atoms are: 0 1"]
        );
    }
    #[test]
    fn three_d_direction_conflict_keeps_trial_orientation_before_warning() {
        let mut s = state();
        s.q.bonds_mut()[1]
            .bond_mut()
            .set_direction(BondDirection::EndUpRight);
        let conf = Conformer3D::new(0, vec![[0.0; 3]; 6], true);
        assert!(matches!(
            wedge_one_three_d_source(
                &mut s,
                &rings(),
                BondId::new(0),
                AtropisomerConformer::ThreeD(&conf)
            ),
            Err(PerceptionError::Rejected(
                AtropisomerRejectionKind::DirectionConflict
            ))
        ));
        assert_eq!(s.q.bonds()[1].begin(), AtomId::new(0));
        assert_eq!(
            s.warnings,
            ["Wedge or hash bond found on atropisomer where not expected - atoms are: 0 1"]
        );
    }
    #[test]
    fn zero_axis_uses_pinned_two_d_typo_and_source_order() {
        let mut s = state();
        let conf = Conformer2D::new(0, vec![[0.0; 2]; 6]);
        assert!(matches!(
            wedge_one_two_d_source(
                &mut s,
                &rings(),
                BondId::new(0),
                AtropisomerConformer::TwoD(&conf)
            ),
            Err(PerceptionError::Rejected(
                AtropisomerRejectionKind::ZeroLengthAxis
            ))
        ));
        assert_eq!(
            s.warnings,
            ["Cound not get a frame of reference for an atropisomer bond - atoms are: 0 1"]
        );
    }
    #[test]
    fn native_unknown_carrier_false_is_silent() {
        let mut s = state();
        s.q.bonds_mut()[1]
            .bond_mut()
            .set_direction(BondDirection::Unknown);
        assert!(matches!(
            wedge_one_no_conformer_source(&mut s, &rings(), BondId::new(0)),
            Err(PerceptionError::Rejected(
                AtropisomerRejectionKind::UnknownCarrierDirection
            ))
        ));
        assert!(s.warnings.is_empty());
    }
}
