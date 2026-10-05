//! Source topology bounds over explicit detached graph and chemistry assignments.
use crate::bounds::{BoundsMatrix, BoundsMatrixError};
use cosmolkit_core::{RingInfo, ValenceAssignment};
use cosmolkit_model::{AtomId, Bond, BondId, BondOrder, BondStereo, Hybridization, TopologyBlock};
use std::collections::HashSet;
use std::f64::consts::PI;
const DIST13_TOL: f64 = 0.04;
const MAX_UPPER: f64 = 1000.0;
const DIST12_DELTA: f64 = 0.01;
#[derive(Debug, thiserror::Error)]
pub enum GraphBoundsError {
    #[error("distance bounds storage: {0}")]
    Bounds(#[from] BoundsMatrixError),
    #[error("UFF rest length: {0}")]
    Uff(#[from] cosmolkit_forcefields::UffBoundsError),
    #[error("topological distances: {0}")]
    Matrix(#[from] cosmolkit_core::MatrixError),
    #[error("source valence: {0}")]
    Valence(#[from] cosmolkit_core::ValenceError),
    #[error("invalid distance bounds: {0}")]
    InvalidBounds(String),
    #[error("{0}")]
    Detail(String),
    #[error("{0}")]
    Input(&'static str),
    #[error("atomic number {0} is outside the periodic table")]
    AtomicNumber(u8),
}
// Source-owned bounds accumulator; packed matrices reuse the unique numeric
// storage. Later complete 1-4/1-5 stages add their path outputs here.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum DistType {
    Dist12,
    Dist13,
    Dist14,
}
#[derive(Debug, Clone)]
struct ComputedData {
    paths14: Vec<Path14Configuration>,
    cis_paths: HashSet<u64>,
    trans_paths: HashSet<u64>,
    set15_atoms: Vec<bool>,
    bond_lengths: Vec<f64>,
    bond_angles: crate::numeric::SymmMatrix,
    bond_adj: crate::numeric::SymmMatrix<i32>,
    visited12_bounds: Vec<bool>,
    visited13_bounds: Vec<bool>,
    visited14_bounds: Vec<bool>,
}
impl ComputedData {
    // Source bitsets are packed; these Vec<bool> rows use one byte per
    // flag. This is a known storage regression retained for independent review.
    fn new(n: usize, e: usize) -> Result<Self, GraphBoundsError> {
        // RDKit❗❌:   ComputedData(unsigned int nAtoms, unsigned int nBonds) {
        // RDKit❗❌:     bondLengths.resize(nBonds);
        // RDKit❗❌:     auto *bAdj = new RDNumeric::IntSymmMatrix(nBonds, -1);
        // RDKit❗❌:     bondAdj.reset(bAdj);
        // RDKit❗❌:     auto *bAngles = new RDNumeric::DoubleSymmMatrix(nBonds, -1.0);
        // RDKit❗❌:     bondAngles.reset(bAngles);
        // RDKit❗❌:     set15Atoms.resize(nAtoms * nAtoms);
        // RDKit❗❌:     visited12Bounds.resize(nAtoms * nAtoms);
        // RDKit❗❌:     visited13Bounds.resize(nAtoms * nAtoms);
        // RDKit❗❌:     visited14Bounds.resize(nAtoms * nAtoms);
        // RDKit❗❌:   }

        let size = n
            .checked_mul(n)
            .ok_or(GraphBoundsError::Input("atom pair count overflow"))?;
        e.checked_add(1)
            .and_then(|v| v.checked_mul(e))
            .ok_or(GraphBoundsError::Input("bond pair count overflow"))?;
        Ok(Self {
            paths14: Vec::new(),
            cis_paths: HashSet::new(),
            trans_paths: HashSet::new(),
            set15_atoms: vec![false; size],
            bond_lengths: vec![0.0; e],
            bond_angles: crate::numeric::SymmMatrix::with_value(e, -1.0),
            bond_adj: crate::numeric::SymmMatrix::with_value(e, -1),
            visited12_bounds: vec![false; size],
            visited13_bounds: vec![false; size],
            visited14_bounds: vec![false; size],
        })
    }
    fn get_bond_angle(&self, _n: usize, i: usize, j: usize) -> f64 {
        self.bond_angles.get_val(i, j)
    }
    fn set_bond_angle(&mut self, _n: usize, i: usize, j: usize, v: f64) {
        self.bond_angles.set_val(i, j, v)
    }
    fn get_bond_adj(&self, _n: usize, i: usize, j: usize) -> i32 {
        self.bond_adj.get_val(i, j)
    }
    fn set_bond_adj(&mut self, _n: usize, i: usize, j: usize, v: i32) {
        self.bond_adj.set_val(i, j, v)
    }
    fn visited_bound(&self, pid: usize, max_dist_type: DistType) -> bool {
        // RDKit❗✔️:   bool visitedBound(unsigned int pid, DistType maxDistType) {
        // RDKit❗✔️:     return ((maxDistType >= DistType::DIST12 && visited12Bounds[pid]) ||
        // RDKit❗✔️:             (maxDistType >= DistType::DIST13 && visited13Bounds[pid]) ||
        // RDKit❗✔️:             (maxDistType >= DistType::DIST14 && visited14Bounds[pid]));
        // RDKit❗✔️:   }
        (max_dist_type >= DistType::Dist12 && self.visited12_bounds[pid])
            || (max_dist_type >= DistType::Dist13 && self.visited13_bounds[pid])
            || (max_dist_type >= DistType::Dist14 && self.visited14_bounds[pid])
    }
}
fn set_12_bounds(
    topology: &TopologyBlock,
    rings: &RingInfo,
    valence: &ValenceAssignment,
    hybridizations: &[cosmolkit_model::Hybridization],
    conjugated: &[bool],
    mmat: &mut BoundsMatrix,
    accum: &mut ComputedData,
) -> Result<(), GraphBoundsError> {
    // RDKit❗❌: void set12Bounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat,
    // RDKit❗❌:                  ComputedData &accumData) {
    // RDKit❗❌:   unsigned int npt = mmat->numRows();
    // RDKit❗❌:   CHECK_INVARIANT(npt == mol.getNumAtoms(), "Wrong size metric matrix");
    // RDKit❗❌:   CHECK_INVARIANT(accumData.bondLengths.size() >= mol.getNumBonds(),
    // RDKit❗❌:                   "Wrong size accumData");
    // RDKit❗❌:   auto [atomParams, foundAll] = UFF::getAtomTypes(mol);
    // RDKit❗❌:   CHECK_INVARIANT(atomParams.size() == mol.getNumAtoms(),
    // RDKit❗❌:                   "parameter vector size mismatch");
    // RDKit❗❌:
    // RDKit❗❌:   boost::dynamic_bitset<> squishAtoms(mol.getNumAtoms());
    // RDKit❗❌:   // find larger heteroatoms in conjugated 5 rings, because we need to add a bit
    // RDKit❗❌:   // of extra flex for them
    // RDKit❗❌:   for (const auto bond : mol.bonds()) {
    // RDKit❗❌:     if (bond->getIsConjugated() &&
    // RDKit❗❌:         (bond->getBeginAtom()->getAtomicNum() > 10 ||
    // RDKit❗❌:          bond->getEndAtom()->getAtomicNum() > 10) &&
    // RDKit❗❌:         mol.getRingInfo() && mol.getRingInfo()->isInitialized() &&
    // RDKit❗❌:         mol.getRingInfo()->isBondInRingOfSize(bond->getIdx(), 5)) {
    // RDKit❗❌:       squishAtoms.set(bond->getBeginAtomIdx());
    // RDKit❗❌:       squishAtoms.set(bond->getEndAtomIdx());
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   for (const auto bond : mol.bonds()) {
    // RDKit❗❌:     auto begId = bond->getBeginAtomIdx();
    // RDKit❗❌:     auto endId = bond->getEndAtomIdx();
    // RDKit❗❌:     auto bOrder = bond->getBondTypeAsDouble();
    // RDKit❗❌:     if (atomParams[begId] && atomParams[endId] && bOrder > 0) {
    // RDKit❗❌:       auto bl = ForceFields::UFF::Utils::calcBondRestLength(
    // RDKit❗❌:           bOrder, atomParams[begId], atomParams[endId]);
    // RDKit❗❌:
    // RDKit❗❌:       double extraSquish = 0.0;
    // RDKit❗❌:       if (squishAtoms[begId] || squishAtoms[endId]) {
    // RDKit❗❌:         extraSquish = 0.2;  // empirical
    // RDKit❗❌:       }
    // RDKit❗❌:
    // RDKit❗❌:       accumData.bondLengths[bond->getIdx()] = bl;
    // RDKit❗❌:       mmat->setUpperBound(begId, endId, bl + extraSquish + DIST12_DELTA);
    // RDKit❗❌:       mmat->setLowerBound(begId, endId, bl - extraSquish - DIST12_DELTA);
    // RDKit❗❌:     } else {
    // RDKit❗❌:       // we don't have parameters for one of the atoms... so we're forced to
    // RDKit❗❌:       // use very crude bounds:
    // RDKit❗❌:       auto vw1 = PeriodicTable::getTable()->getRvdw(
    // RDKit❗❌:           mol.getAtomWithIdx(begId)->getAtomicNum());
    // RDKit❗❌:       auto vw2 = PeriodicTable::getTable()->getRvdw(
    // RDKit❗❌:           mol.getAtomWithIdx(endId)->getAtomicNum());
    // RDKit❗❌:       auto bl = (vw1 + vw2) / 2;
    // RDKit❗❌:       accumData.bondLengths[bond->getIdx()] = bl;
    // RDKit❗❌:       mmat->setUpperBound(begId, endId, 1.5 * bl);
    // RDKit❗❌:       mmat->setLowerBound(begId, endId, .5 * bl);
    // RDKit❗❌:     }
    // RDKit❗❌:     unsigned int pid =
    // RDKit❗❌:         std::min(begId, endId) * mol.getNumAtoms() + std::max(begId, endId);
    // RDKit❗❌:
    // RDKit❗❌:     accumData.visited12Bounds.set(pid);
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // Source-shaped ordered two bond loops and existing indexed ring lookup.
    // Extra E rest-length output/chemical assignments are a known additional
    // allocation cost. No topology/table clone, copied chemical logic or formula.
    let n = topology.atoms.len();
    if mmat.dimension() != n {
        return Err(GraphBoundsError::Input("Wrong size metric matrix"));
    }
    if accum.bond_lengths.len() < topology.bonds.len() {
        return Err(GraphBoundsError::Input("Wrong size accumData"));
    }
    let pair_count = n
        .checked_mul(n)
        .ok_or(GraphBoundsError::Input("atom pair count overflow"))?;
    if accum.visited12_bounds.len() < pair_count {
        return Err(GraphBoundsError::Input("Wrong size visited12Bounds"));
    }
    let lengths = cosmolkit_forcefields::uff_bond_rest_lengths(
        topology,
        valence,
        hybridizations,
        conjugated,
    )?;
    let mut squish = vec![false; n];
    for (idx, bond) in topology.bonds.iter().enumerate() {
        let i = bond.begin().index();
        let j = bond.end().index();
        if conjugated[idx]
            && (topology.atoms[i].atomic_number() > 10 || topology.atoms[j].atomic_number() > 10)
            && rings.is_initialized()
            && rings.is_bond_in_ring_of_size(bond.id(), 5)
        {
            squish[i] = true;
            squish[j] = true;
        }
    }
    for (idx, bond) in topology.bonds.iter().enumerate() {
        let i = bond.begin().index();
        let j = bond.end().index();
        match lengths[idx] {
            Some(bl) => {
                let extra = if squish[i] || squish[j] { 0.2 } else { 0.0 };
                accum.bond_lengths[bond.id().index()] = bl;
                mmat.set_upper(i, j, bl + extra + DIST12_DELTA)?;
                mmat.set_lower(i, j, bl - extra - DIST12_DELTA)?;
            }
            None => {
                let z1 = topology.atoms[i].atomic_number();
                let z2 = topology.atoms[j].atomic_number();
                let vw1 = cosmolkit_core::van_der_waals_radius(z1)
                    .ok_or(GraphBoundsError::AtomicNumber(z1))?;
                let vw2 = cosmolkit_core::van_der_waals_radius(z2)
                    .ok_or(GraphBoundsError::AtomicNumber(z2))?;
                let bl = (vw1 + vw2) / 2.0;
                accum.bond_lengths[bond.id().index()] = bl;
                mmat.set_upper(i, j, 1.5 * bl)?;
                mmat.set_lower(i, j, 0.5 * bl)?;
            }
        }
        accum.visited12_bounds[i.min(j) * n + i.max(j)] = true;
    }
    Ok(())
}
fn bond_between_idx_simple(mol: &TopologyBlock, i: usize, j: usize) -> Option<usize> {
    mol.adjacency
        .neighbors_of(i)
        .iter()
        .find(|v| v.atom_index == j)
        .map(|v| v.bond.index())
}
fn compute_13_dist(bl1: f64, bl2: f64, angle: f64) -> f64 {
    // RDKit❗✔️: inline double compute13Dist(double d1, double d2, double angle) {
    // RDKit❗✔️:   double res = d1 * d1 + d2 * d2 - 2 * d1 * d2 * cos(angle);
    // RDKit❗✔️:   return sqrt(res);
    // RDKit❗✔️: }
    (bl1 * bl1 + bl2 * bl2 - 2.0 * bl1 * bl2 * angle.cos()).sqrt()
}
fn is_larger_sp2_atom_idx(mol: &TopologyBlock, rinfo: &RingInfo, idx: usize) -> bool {
    // RDKit❗✔️: bool isLargerSP2Atom(const Atom *atom) {
    // RDKit❗✔️:   return atom->getAtomicNum() > 13 && atom->getHybridization() == Atom::SP2 &&
    // RDKit❗✔️:          atom->getOwningMol().getRingInfo()->numAtomRings(atom->getIdx());
    // RDKit❗✔️: }
    let atom = &mol.atoms[idx];
    atom.atomic_number() > 13
        && atom.hybridization() == Hybridization::Sp2
        && rinfo.num_atom_rings(AtomId::new(idx)) > 0
}
pub(super) fn init_bounds_mat(
    mmat: &mut BoundsMatrix,
    default_min: f64,
    default_max: f64,
) -> Result<(), GraphBoundsError> {
    // RDKit❗✔️: void initBoundsMat(DistGeom::BoundsMatrix *mmat, double defaultMin,
    // RDKit❗✔️:                    double defaultMax) {
    // RDKit❗✔️:   unsigned int npt = mmat->numRows();
    // RDKit❗✔️:
    // RDKit❗✔️:   for (unsigned int i = 1; i < npt; i++) {
    // RDKit❗✔️:     for (unsigned int j = 0; j < i; j++) {
    // RDKit❗✔️:       mmat->setUpperBound(i, j, defaultMax);
    // RDKit❗✔️:       mmat->setLowerBound(i, j, defaultMin);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    for i in 1..mmat.dimension() {
        for j in 0..i {
            mmat.set_upper(i, j, default_max)?;
            mmat.set_lower(i, j, default_min)?;
        }
    }
    Ok(())
}
fn check_and_set_bounds(
    mmat: &mut BoundsMatrix,
    i: usize,
    j: usize,
    lb: f64,
    ub: f64,
    set_if_better: bool,
) -> Result<(), GraphBoundsError> {
    // RDKit❗✔️: void _checkAndSetBounds(unsigned int i, unsigned int j, double lb, double ub,
    // RDKit❗✔️:                         DistGeom::BoundsMatPtr mmat, bool setIfBetter = false) {
    // RDKit❗✔️:   // get the existing bounds
    // RDKit❗✔️:   double clb = mmat->getLowerBound(i, j);
    // RDKit❗✔️:   double cub = mmat->getUpperBound(i, j);
    // RDKit❗✔️:
    // RDKit❗✔️:   CHECK_INVARIANT(ub > lb, "upper bound not greater than lower bound");
    // RDKit❗✔️:   CHECK_INVARIANT(lb > DIST12_DELTA || clb > DIST12_DELTA, "bad lower bound");
    // RDKit❗✔️:
    // RDKit❗✔️:   // Note: setIfBetter should ONLY be set if the distances are consistent;
    // RDKit❗✔️:   // currently this is not the case, therefore, for now, we are pessimistic on
    // RDKit❗✔️:   // the bounds
    // RDKit❗✔️:   if (setIfBetter) {
    // RDKit❗✔️:     double nlb = std::max(clb, lb);
    // RDKit❗✔️:     double nub = std::min(cub, ub);
    // RDKit❗✔️:
    // RDKit❗✔️:     if (nub <= nlb) {
    // RDKit❗✔️:       // if not overlapping ranges -> be conservative
    // RDKit❗✔️:       nlb = std::min(clb, lb);
    // RDKit❗✔️:       nub = std::max(cub, ub);
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     mmat->setLowerBound(i, j, nlb);
    // RDKit❗✔️:     mmat->setUpperBound(i, j, nub);
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     if (clb <= DIST12_DELTA) {
    // RDKit❗✔️:       mmat->setLowerBound(i, j, lb);
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       if ((lb < clb) && (lb > DIST12_DELTA)) {
    // RDKit❗✔️:         mmat->setLowerBound(i, j, lb);  // conservative bound setting
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     if (cub >= MAX_UPPER) {  // FIX this
    // RDKit❗✔️:       mmat->setUpperBound(i, j, ub);
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       if ((ub > cub) && (ub < MAX_UPPER)) {
    // RDKit❗✔️:         mmat->setUpperBound(i, j, ub);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    let clb = mmat.get_lower(i, j)?;
    let cub = mmat.get_upper(i, j)?;
    if !(ub > lb) {
        return Err(invalid_bounds(
            "upper bound not greater than lower bound",
            i,
            j,
            lb,
            ub,
            clb,
            cub,
        ));
    }
    if !(lb > DIST12_DELTA || clb > DIST12_DELTA) {
        return Err(invalid_bounds("bad lower bound", i, j, lb, ub, clb, cub));
    }
    if set_if_better {
        // Comparisons retain std::min/max first-argument NaN ordering.
        let mut nlb = if clb < lb { lb } else { clb };
        let mut nub = if ub < cub { ub } else { cub };
        if nub <= nlb {
            nlb = if lb < clb { lb } else { clb };
            nub = if cub < ub { ub } else { cub };
        }
        mmat.set_lower(i, j, nlb)?;
        mmat.set_upper(i, j, nub)?;
    } else {
        if clb <= DIST12_DELTA || (lb < clb && lb > DIST12_DELTA) {
            mmat.set_lower(i, j, lb)?;
        }
        if cub >= MAX_UPPER || (ub > cub && ub < MAX_UPPER) {
            mmat.set_upper(i, j, ub)?;
        }
    }
    Ok(())
}
fn set_ring_angle(mol: &TopologyBlock, aid2: usize, ring_size: usize) -> f64 {
    // RDKit❗✔️: void _setRingAngle(Atom::HybridizationType aHyb, unsigned int ringSize,
    // RDKit❗✔️:                    double &angle) {
    // RDKit❗✔️:   // NOTE: this assumes that all angles in a ring are equal. This is
    // RDKit❗✔️:   // certainly not always the case, particular in aromatic rings with
    // RDKit❗✔️:   // heteroatoms
    // RDKit❗✔️:   // like s1cncc1. This led to GitHub55, which was fixed elsewhere.
    // RDKit❗✔️:
    // RDKit❗✔️:   if ((aHyb == Atom::SP2 && ringSize <= 8) || (ringSize == 3) ||
    // RDKit❗✔️:       (ringSize == 4)) {
    // RDKit❗✔️:     angle = M_PI * (1 - 2.0 / ringSize);
    // RDKit❗✔️:   } else if (aHyb == Atom::SP3) {
    // RDKit❗✔️:     if (ringSize == 5) {
    // RDKit❗✔️:       angle = 104 * M_PI / 180;
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       angle = 109.5 * M_PI / 180;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   } else if (aHyb == Atom::SP3D) {
    // RDKit❗✔️:     angle = 105.0 * M_PI / 180;
    // RDKit❗✔️:   } else if (aHyb == Atom::SP3D2) {
    // RDKit❗✔️:     angle = 90.0 * M_PI / 180;
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     angle = 120 * M_PI / 180;
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    let hyb = mol.atoms[aid2].hybridization();
    if (hyb == Hybridization::Sp2 && ring_size <= 8) || ring_size == 3 || ring_size == 4 {
        PI * (1.0 - 2.0 / ring_size as f64)
    } else if hyb == Hybridization::Sp3 {
        if ring_size == 5 {
            104.0 * PI / 180.0
        } else {
            109.5 * PI / 180.0
        }
    } else if hyb == Hybridization::Sp3d {
        105.0 * PI / 180.0
    } else if hyb == Hybridization::Sp3d2 {
        90.0 * PI / 180.0
    } else {
        120.0 * PI / 180.0
    }
}
fn set_13_bounds_helper(
    aid1: usize,
    aid: usize,
    aid3: usize,
    angle: f64,
    bond_lengths: &[f64],
    mmat: &mut BoundsMatrix,
    mol: &TopologyBlock,
    rinfo: &RingInfo,
) -> Result<(), GraphBoundsError> {
    // RDKit❗✔️: void _set13BoundsHelper(unsigned int aid1, unsigned int aid, unsigned int aid3,
    // RDKit❗✔️:                         double angle, const ComputedData &accumData,
    // RDKit❗✔️:                         DistGeom::BoundsMatPtr mmat, const ROMol &mol) {
    // RDKit❗✔️:   auto bid1 = mol.getBondBetweenAtoms(aid1, aid)->getIdx();
    // RDKit❗✔️:   auto bid2 = mol.getBondBetweenAtoms(aid, aid3)->getIdx();
    // RDKit❗✔️:   auto dl = RDGeom::compute13Dist(accumData.bondLengths[bid1],
    // RDKit❗✔️:                                   accumData.bondLengths[bid2], angle);
    // RDKit❗✔️:   auto distTol = DIST13_TOL;
    // RDKit❗✔️:   // Now increase the tolerance if we're outside of the first row of the
    // RDKit❗✔️:   // periodic table.
    // RDKit❗✔️:   if (isLargerSP2Atom(mol.getAtomWithIdx(aid1))) {
    // RDKit❗✔️:     distTol *= 2;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (isLargerSP2Atom(mol.getAtomWithIdx(aid))) {
    // RDKit❗✔️:     distTol *= 2;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (isLargerSP2Atom(mol.getAtomWithIdx(aid3))) {
    // RDKit❗✔️:     distTol *= 2;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   auto du = dl + distTol;
    // RDKit❗✔️:   dl -= distTol;
    // RDKit❗✔️:   _checkAndSetBounds(aid1, aid3, dl, du, mmat);
    // RDKit❗✔️: }
    let bid1 = bond_between_idx_simple(mol, aid1, aid)
        .ok_or(GraphBoundsError::Input("missing first 1-3 path bond"))?;
    let bid2 = bond_between_idx_simple(mol, aid, aid3)
        .ok_or(GraphBoundsError::Input("missing second 1-3 path bond"))?;

    let mut dl = compute_13_dist(bond_lengths[bid1], bond_lengths[bid2], angle);
    let mut dist_tol = DIST13_TOL;

    if is_larger_sp2_atom_idx(mol, rinfo, aid1) {
        dist_tol *= 2.0;
    }
    if is_larger_sp2_atom_idx(mol, rinfo, aid) {
        dist_tol *= 2.0;
    }
    if is_larger_sp2_atom_idx(mol, rinfo, aid3) {
        dist_tol *= 2.0;
    }

    let du = dl + dist_tol;
    dl -= dist_tol;
    check_and_set_bounds(mmat, aid1, aid3, dl, du, false)
}
fn set_13_bounds(
    mol: &TopologyBlock,
    mmat: &mut BoundsMatrix,
    accum_data: &mut ComputedData,
    rinfo: &RingInfo,
) -> Result<(), GraphBoundsError> {
    // RDKit❗✔️: void set13Bounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat,
    // RDKit❗✔️:                  ComputedData &accumData) {
    // RDKit❗✔️:   auto npt = mmat->numRows();
    // RDKit❗✔️:   CHECK_INVARIANT(npt == mol.getNumAtoms(), "Wrong size metric matrix");
    // RDKit❗✔️:   CHECK_INVARIANT(accumData.bondAngles->numRows() == mol.getNumBonds(),
    // RDKit❗✔️:                   "Wrong size bond angle matrix");
    // RDKit❗✔️:   CHECK_INVARIANT(accumData.bondAdj->numRows() == mol.getNumBonds(),
    // RDKit❗✔️:                   "Wrong size bond adjacency matrix");
    // RDKit❗✔️:
    // RDKit❗✔️:   // Since most of the special cases arise out of ring system, we will do
    // RDKit❗✔️:   // the following here:
    // RDKit❗✔️:   // - Loop over all the rings and set the 13 distances between atoms in
    // RDKit❗✔️:   // these rings.
    // RDKit❗✔️:   //   While doing this keep track of the ring atoms that have already been
    // RDKit❗✔️:   //   used as the center atom.
    // RDKit❗✔️:   // - Set the 13 distance between atoms that have a ring atom in between;
    // RDKit❗✔️:   // these can be either non-ring atoms,
    // RDKit❗✔️:   //   or a ring atom and a non-ring atom, or ring atoms that belong to
    // RDKit❗✔️:   //   different simple rings
    // RDKit❗✔️:   // - finally set all other 13 distances
    // RDKit❗✔️:   const auto rinfo = mol.getRingInfo();
    // RDKit❗✔️:   CHECK_INVARIANT(rinfo, "");
    // RDKit❗✔️:
    // RDKit❗✔️:   unsigned int aid2, aid1, aid3, bid1, bid2;
    // RDKit❗✔️:   double angle;
    // RDKit❗✔️:
    // RDKit❗✔️:   auto atomRings = rinfo->atomRings();
    // RDKit❗✔️:   std::sort(atomRings.begin(), atomRings.end(), lessVector);
    // RDKit❗✔️:   // sort the rings based on the ring size
    // RDKit❗✔️:   INT_VECT visited(npt, 0);
    // RDKit❗✔️:
    // RDKit❗✔️:   DOUBLE_VECT angleTaken(npt, 0.0);
    // RDKit❗✔️:   auto nb = mol.getNumBonds();
    // RDKit❗✔️:   BIT_SET donePaths(nb * nb);
    // RDKit❗✔️:   // first deal with all rings and atoms in them
    // RDKit❗✔️:   for (const auto &ringi : atomRings) {
    // RDKit❗✔️:     auto rSize = ringi.size();
    // RDKit❗✔️:     aid1 = ringi[rSize - 1];
    // RDKit❗✔️:     for (unsigned int i = 0; i < rSize; i++) {
    // RDKit❗✔️:       aid2 = ringi[i];
    // RDKit❗✔️:       if (i == rSize - 1) {
    // RDKit❗✔️:         aid3 = ringi[0];
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         aid3 = ringi[i + 1];
    // RDKit❗✔️:       }
    // RDKit❗✔️:       const auto b1 = mol.getBondBetweenAtoms(aid1, aid2);
    // RDKit❗✔️:       const auto b2 = mol.getBondBetweenAtoms(aid2, aid3);
    // RDKit❗✔️:       CHECK_INVARIANT(b1, "no bond found");
    // RDKit❗✔️:       CHECK_INVARIANT(b2, "no bond found");
    // RDKit❗✔️:       bid1 = b1->getIdx();
    // RDKit❗✔️:       bid2 = b2->getIdx();
    // RDKit❗✔️:       auto id1 = nb * bid1 + bid2;
    // RDKit❗✔️:       auto id2 = nb * bid2 + bid1;
    // RDKit❗✔️:
    // RDKit❗✔️:       const auto pid =
    // RDKit❗✔️:           std::min(aid1, aid3) * mol.getNumAtoms() + std::max(aid1, aid3);
    // RDKit❗✔️:
    // RDKit❗✔️:       if ((!donePaths[id1]) && (!donePaths[id2])) {
    // RDKit❗✔️:         // this invar stuff is to deal with bridged systems (Issue 215). In
    // RDKit❗✔️:         // bridged
    // RDKit❗✔️:         // systems we may be covering the same 13 (ring) paths multiple
    // RDKit❗✔️:         // times and unnecessarily increasing the angleTaken at the central
    // RDKit❗✔️:         // atom.
    // RDKit❗✔️:         _setRingAngle(mol.getAtomWithIdx(aid2)->getHybridization(), rSize,
    // RDKit❗✔️:                       angle);
    // RDKit❗✔️:
    // RDKit❗✔️:         if (!accumData.visitedBound(pid, DistType::DIST12)) {
    // RDKit❗✔️:           _set13BoundsHelper(aid1, aid2, aid3, angle, accumData, mmat, mol);
    // RDKit❗✔️:           accumData.visited13Bounds.set(pid);
    // RDKit❗✔️:         }
    // RDKit❗✔️:
    // RDKit❗✔️:         accumData.bondAngles->setVal(bid1, bid2, angle);
    // RDKit❗✔️:         accumData.bondAdj->setVal(bid1, bid2, aid2);
    // RDKit❗✔️:         visited[aid2] += 1;
    // RDKit❗✔️:         angleTaken[aid2] += angle;
    // RDKit❗✔️:         donePaths[id1] = 1;
    // RDKit❗✔️:         donePaths[id2] = 1;
    // RDKit❗✔️:         // donePaths.push_back(invar);
    // RDKit❗✔️:       }
    // RDKit❗✔️:       aid1 = aid2;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // now deal with the remaining atoms
    // RDKit❗✔️:   for (aid2 = 0; aid2 < npt; aid2++) {
    // RDKit❗✔️:     const auto atom = mol.getAtomWithIdx(aid2);
    // RDKit❗✔️:     auto deg = atom->getDegree();
    // RDKit❗✔️:     auto n13 = deg * (deg - 1) / 2;
    // RDKit❗✔️:     if (n13 == static_cast<unsigned int>(visited[aid2])) {
    // RDKit❗✔️:       // we are done with this atom
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     auto ahyb = atom->getHybridization();
    // RDKit❗✔️:     auto [beg1, end1] = mol.getAtomBonds(atom);
    // RDKit❗✔️:     if (visited[aid2] >= 1) {
    // RDKit❗✔️:       // deal with atoms that we already visited; i.e. ring atoms. Set 13
    // RDKit❗✔️:       // distances for one of following cases:
    // RDKit❗✔️:       //  1) Non-ring atoms that have a ring atom in-between
    // RDKit❗✔️:       //  2) Non-ring atom and a ring atom that have a ring atom in between
    // RDKit❗✔️:       //  3) Ring atoms that belong to different rings (that are part of a
    // RDKit❗✔️:       //  fused system
    // RDKit❗✔️:
    // RDKit❗✔️:       while (beg1 != end1) {
    // RDKit❗✔️:         const auto bnd1 = mol[*beg1];
    // RDKit❗✔️:         bid1 = bnd1->getIdx();
    // RDKit❗✔️:         aid1 = bnd1->getOtherAtomIdx(aid2);
    // RDKit❗✔️:         auto [beg2, end2] = mol.getAtomBonds(atom);
    // RDKit❗✔️:         while (beg2 != beg1) {
    // RDKit❗✔️:           const auto bnd2 = mol[*beg2];
    // RDKit❗✔️:           bid2 = bnd2->getIdx();
    // RDKit❗✔️:           aid3 = bnd2->getOtherAtomIdx(aid2);
    // RDKit❗✔️:           if (accumData.bondAngles->getVal(bid1, bid2) < 0.0) {
    // RDKit❗✔️:             // if we haven't dealt with these two bonds before
    // RDKit❗✔️:
    // RDKit❗✔️:             // if we have a sp2 atom things are planar - we simply divide
    // RDKit❗✔️:             // the remaining angle among the remaining 13 configurations
    // RDKit❗✔️:             // (and there should only be one)
    // RDKit❗✔️:             if (ahyb == Atom::SP2) {
    // RDKit❗✔️:               angle = (2 * M_PI - angleTaken[aid2]) / (n13 - visited[aid2]);
    // RDKit❗✔️:             } else if (ahyb == Atom::SP3) {
    // RDKit❗✔️:               // in the case of sp3 we will use the tetrahedral angle mostly
    // RDKit❗✔️:               // - but with some special cases
    // RDKit❗✔️:               angle = 109.5 * M_PI / 180;
    // RDKit❗✔️:               // we will special-case a little bit here for 3, 4 members
    // RDKit❗✔️:               // ring atoms that are sp3 hybridized beyond that the angle
    // RDKit❗✔️:               // reasonably close to the tetrahedral angle
    // RDKit❗✔️:               if (rinfo->isAtomInRingOfSize(aid2, 3)) {
    // RDKit❗✔️:                 angle = 116.0 * M_PI / 180;
    // RDKit❗✔️:               } else if (rinfo->isAtomInRingOfSize(aid2, 4)) {
    // RDKit❗✔️:                 angle = 112.0 * M_PI / 180;
    // RDKit❗✔️:               }
    // RDKit❗✔️:             } else if (Chirality::hasNonTetrahedralStereo(atom)) {
    // RDKit❗✔️:               angle = Chirality::getIdealAngleBetweenLigands(
    // RDKit❗✔️:                           atom, mol.getAtomWithIdx(aid1),
    // RDKit❗✔️:                           mol.getAtomWithIdx(aid3)) *
    // RDKit❗✔️:                       M_PI / 180;
    // RDKit❗✔️:             } else {
    // RDKit❗✔️:               // other options we will simply based things on the number of
    // RDKit❗✔️:               // substituent
    // RDKit❗✔️:               if (deg == 5) {
    // RDKit❗✔️:                 angle = 105.0 * M_PI / 180;
    // RDKit❗✔️:               } else if (deg == 6) {
    // RDKit❗✔️:                 angle = 135.0 * M_PI / 180;
    // RDKit❗✔️:               } else {
    // RDKit❗✔️:                 angle = 120.0 * M_PI / 180;  // FIX: this default is probably
    // RDKit❗✔️:                                              // not the best we can do here
    // RDKit❗✔️:               }
    // RDKit❗✔️:             }
    // RDKit❗✔️:
    // RDKit❗✔️:             const unsigned int pid =
    // RDKit❗✔️:                 std::min(aid1, aid3) * mol.getNumAtoms() + std::max(aid1, aid3);
    // RDKit❗✔️:
    // RDKit❗✔️:             if (!accumData.visitedBound(pid, DistType::DIST12)) {
    // RDKit❗✔️:               _set13BoundsHelper(aid1, aid2, aid3, angle, accumData, mmat, mol);
    // RDKit❗✔️:               accumData.visited13Bounds.set(pid);
    // RDKit❗✔️:             }
    // RDKit❗✔️:
    // RDKit❗✔️:             accumData.bondAngles->setVal(bid1, bid2, angle);
    // RDKit❗✔️:             accumData.bondAdj->setVal(bid1, bid2, aid2);
    // RDKit❗✔️:             angleTaken[aid2] += angle;
    // RDKit❗✔️:             visited[aid2] += 1;
    // RDKit❗✔️:           }
    // RDKit❗✔️:           ++beg2;
    // RDKit❗✔️:         }  // while loop over the second bond
    // RDKit❗✔️:         ++beg1;
    // RDKit❗✔️:       }  // while loop over the first bond
    // RDKit❗✔️:     } else if (visited[aid2] == 0) {
    // RDKit❗✔️:       // non-ring atoms - we will simply use angles based on hybridization
    // RDKit❗✔️:       while (beg1 != end1) {
    // RDKit❗✔️:         const auto bnd1 = mol[*beg1];
    // RDKit❗✔️:         bid1 = bnd1->getIdx();
    // RDKit❗✔️:         aid1 = bnd1->getOtherAtomIdx(aid2);
    // RDKit❗✔️:         auto [beg2, end2] = mol.getAtomBonds(atom);
    // RDKit❗✔️:         while (beg2 != beg1) {
    // RDKit❗✔️:           const auto bnd2 = mol[*beg2];
    // RDKit❗✔️:           bid2 = bnd2->getIdx();
    // RDKit❗✔️:           aid3 = bnd2->getOtherAtomIdx(aid2);
    // RDKit❗✔️:           if (Chirality::hasNonTetrahedralStereo(atom)) {
    // RDKit❗✔️:             angle =
    // RDKit❗✔️:                 Chirality::getIdealAngleBetweenLigands(
    // RDKit❗✔️:                     atom, mol.getAtomWithIdx(aid1), mol.getAtomWithIdx(aid3)) *
    // RDKit❗✔️:                 M_PI / 180;
    // RDKit❗✔️:
    // RDKit❗✔️:           } else {
    // RDKit❗✔️:             if (ahyb == Atom::SP) {
    // RDKit❗✔️:               angle = M_PI;
    // RDKit❗✔️:             } else if (ahyb == Atom::SP2) {
    // RDKit❗✔️:               angle = 2 * M_PI / 3;
    // RDKit❗✔️:             } else if (ahyb == Atom::SP3) {
    // RDKit❗✔️:               angle = 109.5 * M_PI / 180;
    // RDKit❗✔️:             } else if (Chirality::hasNonTetrahedralStereo(atom)) {
    // RDKit❗✔️:               angle = Chirality::getIdealAngleBetweenLigands(
    // RDKit❗✔️:                           atom, mol.getAtomWithIdx(aid1),
    // RDKit❗✔️:                           mol.getAtomWithIdx(aid3)) *
    // RDKit❗✔️:                       M_PI / 180;
    // RDKit❗✔️:             } else if (ahyb == Atom::SP3D) {
    // RDKit❗✔️:               // FIX: this and the remaining two hybridization states below
    // RDKit❗✔️:               // should probably be special cased. These defaults below are
    // RDKit❗✔️:               // probably not the best we can do particularly when stereo
    // RDKit❗✔️:               // chemistry is know
    // RDKit❗✔️:               angle = 105.0 * M_PI / 180;
    // RDKit❗✔️:             } else if (ahyb == Atom::SP3D2) {
    // RDKit❗✔️:               angle = 135.0 * M_PI / 180;
    // RDKit❗✔️:             } else {
    // RDKit❗✔️:               angle = 120.0 * M_PI / 180;
    // RDKit❗✔️:             }
    // RDKit❗✔️:           }
    // RDKit❗✔️:           const unsigned int pid =
    // RDKit❗✔️:               std::min(aid1, aid3) * mol.getNumAtoms() + std::max(aid1, aid3);
    // RDKit❗✔️:
    // RDKit❗✔️:           if (!accumData.visitedBound(pid, DistType::DIST12)) {
    // RDKit❗✔️:             if (atom->getDegree() <= 4 ||
    // RDKit❗✔️:                 (Chirality::hasNonTetrahedralStereo(atom) &&
    // RDKit❗✔️:                  atom->hasProp(common_properties::_chiralPermutation))) {
    // RDKit❗✔️:               _set13BoundsHelper(aid1, aid2, aid3, angle, accumData, mmat, mol);
    // RDKit❗✔️:             } else {
    // RDKit❗✔️:               // just use 180 as the max angle and an arbitrary min angle
    // RDKit❗✔️:               auto dmax =
    // RDKit❗✔️:                   accumData.bondLengths[bid1] + accumData.bondLengths[bid2];
    // RDKit❗✔️:               auto dl = 1.0;
    // RDKit❗✔️:               auto du = dmax * 1.2;
    // RDKit❗✔️:               _checkAndSetBounds(aid1, aid3, dl, du, mmat);
    // RDKit❗✔️:             }
    // RDKit❗✔️:             accumData.visited13Bounds.set(pid);
    // RDKit❗✔️:           }
    // RDKit❗✔️:
    // RDKit❗✔️:           accumData.bondAngles->setVal(bid1, bid2, angle);
    // RDKit❗✔️:           accumData.bondAdj->setVal(bid1, bid2, aid2);
    // RDKit❗✔️:           angleTaken[aid2] += angle;
    // RDKit❗✔️:           visited[aid2] += 1;
    // RDKit❗✔️:           ++beg2;
    // RDKit❗✔️:         }  // while loop over second bond
    // RDKit❗✔️:         ++beg1;
    // RDKit❗✔️:       }  // while loop over first bond
    // RDKit❗✔️:     }  // done with non-ring atoms
    // RDKit❗✔️:   }  // done with all atoms
    // RDKit❗✔️: }
    // The same source ring-copy/sort, dense visited path bitset and nested
    // bond iteration. Canonical adjacency borrows preserve bond order and
    // avoid the inherited full-bond scans/temporary neighbor arrays. Packed
    // symmetric integer/f64 storage is the existing numeric implementation.
    let npt = mmat.dimension();
    if npt != mol.atoms.len() {
        return Err(GraphBoundsError::Input("Wrong size metric matrix"));
    }

    let nb = mol.bonds.len();
    if accum_data.bond_angles.num_rows() != nb {
        return Err(GraphBoundsError::Input("Wrong size bond angle matrix"));
    }
    if accum_data.bond_adj.num_rows() != nb {
        return Err(GraphBoundsError::Input("Wrong size bond adjacency matrix"));
    }

    let mut atom_rings: Vec<Vec<usize>> = rinfo
        .atom_rings()
        .iter()
        .map(|ring| ring.iter().map(|aid| aid.index()).collect())
        .collect();
    atom_rings.sort_by_key(Vec::len);

    let mut visited = vec![0usize; npt];
    let mut angle_taken = vec![0.0; npt];
    let done_count = nb
        .checked_mul(nb)
        .ok_or(GraphBoundsError::Input("bond path count overflow"))?;
    let mut done_paths = vec![false; done_count];

    for ring in &atom_rings {
        let r_size = ring.len();
        let mut aid1 = ring[r_size - 1];
        for i in 0..r_size {
            let aid2 = ring[i];
            let aid3 = if i == r_size - 1 {
                ring[0]
            } else {
                ring[i + 1]
            };
            let b1 = bond_between_idx_simple(mol, aid1, aid2)
                .ok_or(GraphBoundsError::Input("no bond found"))?;
            let b2 = bond_between_idx_simple(mol, aid2, aid3)
                .ok_or(GraphBoundsError::Input("no bond found"))?;
            let id1 = nb * b1 + b2;
            let id2 = nb * b2 + b1;
            let pid = aid1.min(aid3) * npt + aid1.max(aid3);

            if !done_paths[id1] && !done_paths[id2] {
                let angle = set_ring_angle(mol, aid2, r_size);
                if !accum_data.visited_bound(pid, DistType::Dist12) {
                    set_13_bounds_helper(
                        aid1,
                        aid2,
                        aid3,
                        angle,
                        &accum_data.bond_lengths,
                        mmat,
                        mol,
                        rinfo,
                    )?;
                    accum_data.visited13_bounds[pid] = true;
                }
                accum_data.set_bond_angle(nb, b1, b2, angle);
                accum_data.set_bond_adj(nb, b1, b2, aid2 as i32);
                visited[aid2] += 1;
                angle_taken[aid2] += angle;
                done_paths[id1] = true;
                done_paths[id2] = true;
            }
            aid1 = aid2;
        }
    }

    for aid2 in 0..npt {
        let atom = &mol.atoms[aid2];
        let nbrs = mol.adjacency.neighbors_of(aid2);
        let deg = nbrs.len();
        let n13 = deg * (deg.saturating_sub(1)) / 2;
        if n13 == visited[aid2] {
            continue;
        }
        let ahyb = atom.hybridization();

        if visited[aid2] >= 1 {
            for left in 0..nbrs.len() {
                let aid1 = nbrs[left].atom_index;
                let bid1 = nbrs[left].bond.index();
                for right in 0..left {
                    let aid3 = nbrs[right].atom_index;
                    let bid2 = nbrs[right].bond.index();
                    if !(accum_data.get_bond_angle(nb, bid1, bid2) < 0.0) {
                        continue;
                    }

                    let angle = if ahyb == Hybridization::Sp2 {
                        (2.0 * std::f64::consts::PI - angle_taken[aid2])
                            / (n13 - visited[aid2]) as f64
                    } else if ahyb == Hybridization::Sp3 {
                        if rinfo.is_atom_in_ring_of_size(AtomId::new(aid2), 3) {
                            116.0 * PI / 180.0
                        } else if rinfo.is_atom_in_ring_of_size(AtomId::new(aid2), 4) {
                            112.0 * PI / 180.0
                        } else {
                            109.5 * PI / 180.0
                        }
                    } else if cosmolkit_core::parser_stereo_order::nontetrahedral_max_neighbors(
                        atom.chiral_tag(),
                    )
                    .is_some()
                    {
                        cosmolkit_core::non_tetrahedral_ideal_angle(mol, aid2, aid1, aid3) * PI
                            / 180.0
                    } else if deg == 5 {
                        105.0 * PI / 180.0
                    } else if deg == 6 {
                        135.0 * PI / 180.0
                    } else {
                        120.0 * PI / 180.0
                    };

                    let pid = aid1.min(aid3) * npt + aid1.max(aid3);
                    if !accum_data.visited_bound(pid, DistType::Dist12) {
                        set_13_bounds_helper(
                            aid1,
                            aid2,
                            aid3,
                            angle,
                            &accum_data.bond_lengths,
                            mmat,
                            mol,
                            rinfo,
                        )?;
                        accum_data.visited13_bounds[pid] = true;
                    }

                    accum_data.set_bond_angle(nb, bid1, bid2, angle);
                    accum_data.set_bond_adj(nb, bid1, bid2, aid2 as i32);
                    angle_taken[aid2] += angle;
                    visited[aid2] += 1;
                }
            }
        } else {
            for left in 0..nbrs.len() {
                let aid1 = nbrs[left].atom_index;
                let bid1 = nbrs[left].bond.index();
                for right in 0..left {
                    let aid3 = nbrs[right].atom_index;
                    let bid2 = nbrs[right].bond.index();

                    let angle =
                        if cosmolkit_core::parser_stereo_order::nontetrahedral_max_neighbors(
                            atom.chiral_tag(),
                        )
                        .is_some()
                        {
                            cosmolkit_core::non_tetrahedral_ideal_angle(mol, aid2, aid1, aid3) * PI
                                / 180.0
                        } else if ahyb == Hybridization::Sp {
                            std::f64::consts::PI
                        } else if ahyb == Hybridization::Sp2 {
                            2.0 * std::f64::consts::PI / 3.0
                        } else if ahyb == Hybridization::Sp3 {
                            109.5 * PI / 180.0
                        } else if ahyb == Hybridization::Sp3d {
                            105.0 * PI / 180.0
                        } else if ahyb == Hybridization::Sp3d2 {
                            135.0 * PI / 180.0
                        } else {
                            120.0 * PI / 180.0
                        };

                    let pid = aid1.min(aid3) * npt + aid1.max(aid3);
                    if !accum_data.visited_bound(pid, DistType::Dist12) {
                        if deg <= 4
                            || (cosmolkit_core::parser_stereo_order::nontetrahedral_max_neighbors(
                                atom.chiral_tag(),
                            )
                            .is_some()
                                && atom.chiral_permutation().is_some())
                        {
                            set_13_bounds_helper(
                                aid1,
                                aid2,
                                aid3,
                                angle,
                                &accum_data.bond_lengths,
                                mmat,
                                mol,
                                rinfo,
                            )?;
                        } else {
                            let dmax =
                                accum_data.bond_lengths[bid1] + accum_data.bond_lengths[bid2];
                            let dl = 1.0;
                            let du = dmax * 1.2;
                            check_and_set_bounds(mmat, aid1, aid3, dl, du, false)?;
                        }
                        accum_data.visited13_bounds[pid] = true;
                    }

                    accum_data.set_bond_angle(nb, bid1, bid2, angle);
                    accum_data.set_bond_adj(nb, bid1, bid2, aid2 as i32);
                    angle_taken[aid2] += angle;
                    visited[aid2] += 1;
                }
            }
        }
    }
    Ok(())
}

const GEN_DIST_TOL: f64 = 0.06;
const MIN_MACROCYCLE_RING_SIZE: usize = 9;
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Path14Kind {
    Cis,
    Trans,
    Other,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Path14Configuration {
    bid1: usize,
    bid2: usize,
    bid3: usize,
    kind: Path14Kind,
}
fn path14_id(nb: usize, bid1: usize, bid2: usize, bid3: usize) -> u64 {
    bid1 as u64 * nb as u64 * nb as u64 + bid2 as u64 * nb as u64 + bid3 as u64
}
fn record_path_flag(paths: &mut HashSet<u64>, id: u64) {
    paths.insert(id);
}
fn has_path_flag(paths: &HashSet<u64>, id: u64) -> bool {
    paths.contains(&id)
}
fn bond_pair_shared_atom(
    mol: &TopologyBlock,
    accum_data: &ComputedData,
    bid1: usize,
    bid2: usize,
) -> Result<usize, GraphBoundsError> {
    let nb = mol.bonds.len();
    let aid = accum_data.get_bond_adj(nb, bid1, bid2);
    if aid < 0 {
        return Err(GraphBoundsError::Detail(format!(
            "missing shared atom for bond pair ({bid1}, {bid2})"
        )));
    }
    let aid = aid as usize;
    if aid >= mol.atoms.len() {
        return Err(GraphBoundsError::Detail(format!(
            "shared atom index {aid} for bond pair ({bid1}, {bid2}) is out of range"
        )));
    }
    Ok(aid)
}
fn validate_bond_angle(
    angle: f64,
    bid1: usize,
    bid2: usize,
    context: &'static str,
) -> Result<(), GraphBoundsError> {
    if angle > 0.0 {
        Ok(())
    } else {
        Err(GraphBoundsError::Detail(format!(
            "{context}: missing or invalid bond angle for bond pair ({bid1}, {bid2}): {angle}"
        )))
    }
}
fn compute_14_dist_3d(d1: f64, d2: f64, d3: f64, ang12: f64, ang23: f64, tor_ang: f64) -> f64 {
    // RDKit❗✔️: inline double compute14Dist3D(double d1, double d2, double d3, double ang12,
    // RDKit❗✔️:                               double ang23, double torAng) {
    // RDKit❗✔️:   // location of atom1
    // RDKit❗✔️:   Point3D p1(d1 * cos(ang12), d1 * sin(ang12), 0.0);
    // RDKit❗✔️:
    // RDKit❗✔️:   // location of atom 4 if the rosion angle was 0
    // RDKit❗✔️:   Point3D p4(d2 - d3 * cos(ang23), d3 * sin(ang23), 0.0);
    // RDKit❗✔️:
    // RDKit❗✔️:   // now we will rotate p4 about the x-axis by the desired torsion angle
    // RDKit❗✔️:   Transform3D trans;
    // RDKit❗✔️:   trans.SetRotation(torAng, X_Axis);
    // RDKit❗✔️:   trans.TransformPoint(p4);
    // RDKit❗✔️:
    // RDKit❗✔️:   // find the distance
    // RDKit❗✔️:   p4 -= p1;
    // RDKit❗✔️:   return p4.length();
    // RDKit❗✔️: }

    let p1x = d1 * ang12.cos();
    let p1y = d1 * ang12.sin();
    let p4x = d2 - d3 * ang23.cos();
    let p4y = d3 * ang23.sin() * tor_ang.cos();
    let p4z = d3 * ang23.sin() * tor_ang.sin();
    let dx = p4x - p1x;
    let dy = p4y - p1y;
    let dz = p4z;
    (dx * dx + dy * dy + dz * dz).sqrt()
}
fn compute_14_dist_cis(d1: f64, d2: f64, d3: f64, ang12: f64, ang23: f64) -> f64 {
    // RDKit❗✔️: inline double compute14DistCis(double d1, double d2, double d3, double ang12,
    // RDKit❗✔️:                                double ang23) {
    // RDKit❗✔️:   double dx = d2 - d3 * cos(ang23) - d1 * cos(ang12);
    // RDKit❗✔️:   double dy = d3 * sin(ang23) - d1 * sin(ang12);
    // RDKit❗✔️:   double res = dx * dx + dy * dy;
    // RDKit❗✔️:   return sqrt(res);
    // RDKit❗✔️: }

    let dx = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy = d3 * ang23.sin() - d1 * ang12.sin();
    (dx * dx + dy * dy).sqrt()
}
fn compute_14_dist_trans(d1: f64, d2: f64, d3: f64, ang12: f64, ang23: f64) -> f64 {
    // RDKit❗✔️: inline double compute14DistTrans(double d1, double d2, double d3, double ang12,
    // RDKit❗✔️:                                  double ang23) {
    // RDKit❗✔️:   double dx = d2 - d3 * cos(ang23) - d1 * cos(ang12);
    // RDKit❗✔️:   double dy = d3 * sin(ang23) + d1 * sin(ang12);
    // RDKit❗✔️:   double res = dx * dx + dy * dy;
    // RDKit❗✔️:   return sqrt(res);
    // RDKit❗✔️: }

    let dx = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy = d3 * ang23.sin() + d1 * ang12.sin();
    (dx * dx + dy * dy).sqrt()
}
fn get_atom_stereo(bond: &Bond, aid1: usize, aid4: usize) -> BondStereo {
    // RDKit❗✔️: Bond::BondStereo _getAtomStereo(const Bond *bnd, unsigned int aid1,
    // RDKit❗✔️:                                 unsigned int aid4) {
    // RDKit❗✔️:   auto stype = bnd->getStereo();
    // RDKit❗✔️:   if (stype > Bond::STEREOANY && bnd->getStereoAtoms().size() >= 2) {
    // RDKit❗✔️:     const auto &stAtoms = bnd->getStereoAtoms();
    // RDKit❗✔️:     if ((static_cast<unsigned int>(stAtoms[0]) != aid1) ^
    // RDKit❗✔️:         (static_cast<unsigned int>(stAtoms[1]) != aid4)) {
    // RDKit❗✔️:       switch (stype) {
    // RDKit❗✔️:         case Bond::STEREOZ:
    // RDKit❗✔️:           stype = Bond::STEREOE;
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         case Bond::STEREOE:
    // RDKit❗✔️:           stype = Bond::STEREOZ;
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         case Bond::STEREOCIS:
    // RDKit❗✔️:           stype = Bond::STEREOTRANS;
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         case Bond::STEREOTRANS:
    // RDKit❗✔️:           stype = Bond::STEREOCIS;
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         default:
    // RDKit❗✔️:           break;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return stype;
    // RDKit❗✔️: }

    let mut stype = bond.stereo();
    if matches!(
        stype,
        BondStereo::Z | BondStereo::E | BondStereo::Cis | BondStereo::Trans
    ) && bond.stereo_atoms().is_some_and(|atoms| atoms.len() >= 2)
    {
        let stereo_atoms = bond.stereo_atoms().expect("checked stereo atoms");
        let needs_flip = (stereo_atoms[0].index() != aid1) ^ (stereo_atoms[1].index() != aid4);
        if needs_flip {
            stype = match stype {
                BondStereo::Z => BondStereo::E,
                BondStereo::E => BondStereo::Z,
                BondStereo::Cis => BondStereo::Trans,
                BondStereo::Trans => BondStereo::Cis,
                other => other,
            };
        }
    }
    stype
}
fn record_14_path(
    mol: &TopologyBlock,
    bid1: usize,
    bid2: usize,
    bid3: usize,
    accum_data: &mut ComputedData,
) -> Result<(), GraphBoundsError> {
    // RDKit❗✔️: void _record14Path(const ROMol &mol, unsigned int bid1, unsigned int bid2,
    // RDKit❗✔️:                    unsigned int bid3, ComputedData &accumData) {
    // RDKit❗✔️:   const Atom *atm2 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid1, bid2));
    // RDKit❗✔️:   PRECONDITION(atm2, "");
    // RDKit❗✔️:   Atom::HybridizationType ahyb2 = atm2->getHybridization();
    // RDKit❗✔️:   const Atom *atm3 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid2, bid3));
    // RDKit❗✔️:   PRECONDITION(atm3, "");
    // RDKit❗✔️:   Atom::HybridizationType ahyb3 = atm3->getHybridization();
    // RDKit❗✔️:   unsigned int nb = mol.getNumBonds();
    // RDKit❗✔️:   Path14Configuration path14;
    // RDKit❗✔️:   path14.bid1 = bid1;
    // RDKit❗✔️:   path14.bid2 = bid2;
    // RDKit❗✔️:   path14.bid3 = bid3;
    // RDKit❗✔️:   if ((ahyb2 == Atom::SP2) && (ahyb3 == Atom::SP2)) {  // FIX: check for trans
    // RDKit❗✔️:     path14.type = Path14Configuration::CIS;
    // RDKit❗✔️:     accumData.cisPaths.insert(static_cast<unsigned long>(bid1) * nb * nb +
    // RDKit❗✔️:                               bid2 * nb + bid3);
    // RDKit❗✔️:     accumData.cisPaths.insert(static_cast<unsigned long>(bid3) * nb * nb +
    // RDKit❗✔️:                               bid2 * nb + bid1);
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     path14.type = Path14Configuration::OTHER;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   accumData.paths14.push_back(path14);
    // RDKit❗✔️: }

    let atm2 = bond_pair_shared_atom(mol, accum_data, bid1, bid2)?;
    let ahyb2 = mol.atoms[atm2].hybridization();
    let atm3 = bond_pair_shared_atom(mol, accum_data, bid2, bid3)?;
    let ahyb3 = mol.atoms[atm3].hybridization();
    let nb = mol.bonds.len();

    let kind = if ahyb2 == Hybridization::Sp2 && ahyb3 == Hybridization::Sp2 {
        record_path_flag(&mut accum_data.cis_paths, path14_id(nb, bid1, bid2, bid3));
        record_path_flag(&mut accum_data.cis_paths, path14_id(nb, bid3, bid2, bid1));
        Path14Kind::Cis
    } else {
        Path14Kind::Other
    };

    accum_data.paths14.push(Path14Configuration {
        bid1,
        bid2,
        bid3,
        kind,
    });
    Ok(())
}
fn set_in_ring_14_bounds(
    mol: &TopologyBlock,
    bid1: usize,
    bid2: usize,
    bid3: usize,
    accum_data: &mut ComputedData,
    mmat: &mut BoundsMatrix,
    dmat: &[f64],
    ring_size: usize,
    ring_info: &RingInfo,
) -> Result<(), GraphBoundsError> {
    // RDKit❗✔️: void _setInRing14Bounds(const ROMol &mol, const Bond *bnd1, const Bond *bnd2,
    // RDKit❗✔️:                         const Bond *bnd3, ComputedData &accumData,
    // RDKit❗✔️:                         DistGeom::BoundsMatPtr mmat, double *dmat,
    // RDKit❗✔️:                         int ringSize) {
    // RDKit❗✔️:   PRECONDITION(bnd1, "");
    // RDKit❗✔️:   PRECONDITION(bnd2, "");
    // RDKit❗✔️:   PRECONDITION(bnd3, "");
    // RDKit❗✔️:   unsigned int bid1, bid2, bid3;
    // RDKit❗✔️:   bid1 = bnd1->getIdx();
    // RDKit❗✔️:   bid2 = bnd2->getIdx();
    // RDKit❗✔️:   bid3 = bnd3->getIdx();
    // RDKit❗✔️:   const Atom *atm2 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid1, bid2));
    // RDKit❗✔️:   PRECONDITION(atm2, "");
    // RDKit❗✔️:   Atom::HybridizationType ahyb2 = atm2->getHybridization();
    // RDKit❗✔️:   const Atom *atm3 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid2, bid3));
    // RDKit❗✔️:   PRECONDITION(atm3, "");
    // RDKit❗✔️:   Atom::HybridizationType ahyb3 = atm3->getHybridization();
    // RDKit❗✔️:
    // RDKit❗✔️:   unsigned int aid1 = bnd1->getOtherAtomIdx(atm2->getIdx());
    // RDKit❗✔️:   unsigned int aid4 = bnd3->getOtherAtomIdx(atm3->getIdx());
    // RDKit❗✔️:
    // RDKit❗✔️:   const unsigned int pid =
    // RDKit❗✔️:       std::min(aid1, aid4) * mol.getNumAtoms() + std::max(aid1, aid4);
    // RDKit❗✔️:
    // RDKit❗✔️:   if (accumData.visitedBound(pid, DistType::DIST13)) {
    // RDKit❗✔️:     // if this is already a 1-3 or 1-2 distance; do not overwrite
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // check that this actually is a 1-4 contact:
    // RDKit❗✔️:   if (dmat[std::max(aid1, aid4) * mmat->numRows() + std::min(aid1, aid4)] <
    // RDKit❗✔️:       2.9) {
    // RDKit❗✔️:     // std::cerr<<"skip: "<<aid1<<"-"<<aid4<<" because
    // RDKit❗✔️:     // d="<<dmat[std::max(aid1,aid4)*mmat->numRows()+std::min(aid1,aid4)]<<std::endl;
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   double bl1 = accumData.bondLengths[bid1];
    // RDKit❗✔️:   double bl2 = accumData.bondLengths[bid2];
    // RDKit❗✔️:   double bl3 = accumData.bondLengths[bid3];
    // RDKit❗✔️:
    // RDKit❗✔️:   double ba12 = accumData.bondAngles->getVal(bid1, bid2);
    // RDKit❗✔️:   double ba23 = accumData.bondAngles->getVal(bid2, bid3);
    // RDKit❗✔️:
    // RDKit❗✔️:   CHECK_INVARIANT(ba12 > 0.0, "");
    // RDKit❗✔️:   CHECK_INVARIANT(ba23 > 0.0, "");
    // RDKit❗✔️:   double dl, du;
    // RDKit❗✔️:   unsigned int nb = mol.getNumBonds();
    // RDKit❗✔️:   // several special cases here
    // RDKit❗✔️:   Path14Configuration path14;
    // RDKit❗✔️:   path14.bid1 = bid1;
    // RDKit❗✔️:   path14.bid2 = bid2;
    // RDKit❗✔️:   path14.bid3 = bid3;
    // RDKit❗✔️:   Bond::BondStereo stype = _getAtomStereo(bnd2, aid1, aid4);
    // RDKit❗✔️:   bool preferCis = false;
    // RDKit❗✔️:   bool preferTrans = false;
    // RDKit❗✔️:
    // RDKit❗✔️:   // we add a check for the ring size here because there's no reason to
    // RDKit❗✔️:   // assume cis bonds in bigger rings. This was part of github #1240:
    // RDKit❗✔️:   // failure to embed larger aromatic rings
    // RDKit❗✔️:   if (ringSize <= 8 && (ahyb2 == Atom::SP2) && (ahyb3 == Atom::SP2) &&
    // RDKit❗✔️:       (stype != Bond::STEREOE && stype != Bond::STEREOTRANS)) {
    // RDKit❗✔️:     // the ring check here was a big part of github #697
    // RDKit❗✔️:     if (mol.getRingInfo()->numBondRings(bid2) > 1) {
    // RDKit❗✔️:       if (mol.getRingInfo()->numBondRings(bid1) == 1 &&
    // RDKit❗✔️:           mol.getRingInfo()->numBondRings(bid3) == 1) {
    // RDKit❗✔️:         for (const auto &br : mol.getRingInfo()->bondRings()) {
    // RDKit❗✔️:           if (std::find(br.begin(), br.end(), bid1) != br.end()) {
    // RDKit❗✔️:             if (std::find(br.begin(), br.end(), bid3) != br.end()) {
    // RDKit❗✔️:               preferCis = true;
    // RDKit❗✔️:             }
    // RDKit❗✔️:             break;
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       preferCis = true;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   } else if (stype == Bond::STEREOZ || stype == Bond::STEREOCIS) {
    // RDKit❗✔️:     preferCis = true;
    // RDKit❗✔️:   } else if (stype == Bond::STEREOE || stype == Bond::STEREOTRANS) {
    // RDKit❗✔️:     preferTrans = true;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   if (preferCis) {
    // RDKit❗✔️:     path14.type = Path14Configuration::CIS;
    // RDKit❗✔️:     accumData.cisPaths.insert(static_cast<unsigned long>(bid1) * nb * nb +
    // RDKit❗✔️:                               bid2 * nb + bid3);
    // RDKit❗✔️:     accumData.cisPaths.insert(static_cast<unsigned long>(bid3) * nb * nb +
    // RDKit❗✔️:                               bid2 * nb + bid1);
    // RDKit❗✔️:   } else if (preferTrans) {
    // RDKit❗✔️:     path14.type = Path14Configuration::TRANS;
    // RDKit❗✔️:     accumData.transPaths.insert(static_cast<unsigned long>(bid1) * nb * nb +
    // RDKit❗✔️:                                 bid2 * nb + bid3);
    // RDKit❗✔️:     accumData.transPaths.insert(static_cast<unsigned long>(bid3) * nb * nb +
    // RDKit❗✔️:                                 bid2 * nb + bid1);
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     path14.type = Path14Configuration::OTHER;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   accumData.paths14.push_back(path14);
    // RDKit❗✔️:
    // RDKit❗✔️:   // std::cerr << "  torsion: " << aid1 << " " << aid4 << ": " << preferCis
    // RDKit❗✔️:   // << "
    // RDKit❗✔️:   // "
    // RDKit❗✔️:   //           << preferTrans << std::endl;
    // RDKit❗✔️:   if (preferCis) {
    // RDKit❗✔️:     dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
    // RDKit❗✔️:     du = dl + 2 * GEN_DIST_TOL;
    // RDKit❗✔️:   } else if (preferTrans) {
    // RDKit❗✔️:     dl = RDGeom::compute14DistTrans(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
    // RDKit❗✔️:     du = dl + 2 * GEN_DIST_TOL;
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     // basically we will assume 0 to 180 allowed
    // RDKit❗✔️:     dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:     du = RDGeom::compute14DistTrans(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:     if (du < dl) {
    // RDKit❗✔️:       std::swap(du, dl);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (fabs(du - dl) < DIST12_DELTA) {
    // RDKit❗✔️:       dl -= GEN_DIST_TOL;
    // RDKit❗✔️:       du += GEN_DIST_TOL;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   // std::cerr << "7: " << aid1 << "-" << aid4 << std::endl;
    // RDKit❗✔️:
    // RDKit❗✔️:   accumData.visited14Bounds.set(pid);
    // RDKit❗✔️:   _checkAndSetBounds(aid1, aid4, dl, du, mmat);
    // RDKit❗✔️: }

    let atm2 = bond_pair_shared_atom(mol, accum_data, bid1, bid2)?;
    let ahyb2 = mol.atoms[atm2].hybridization();
    let atm3 = bond_pair_shared_atom(mol, accum_data, bid2, bid3)?;
    let ahyb3 = mol.atoms[atm3].hybridization();

    let bnd1 = &mol.bonds[bid1];
    let bnd3 = &mol.bonds[bid3];
    let aid1 = if bnd1.begin().index() == atm2 {
        bnd1.end().index()
    } else {
        bnd1.begin().index()
    };
    let aid4 = if bnd3.begin().index() == atm3 {
        bnd3.end().index()
    } else {
        bnd3.begin().index()
    };
    let pid = aid1.min(aid4) * mol.atoms.len() + aid1.max(aid4);
    if accum_data.visited_bound(pid, DistType::Dist13) {
        return Ok(());
    }
    if dmat[aid1.max(aid4) * mmat.dimension() + aid1.min(aid4)] < 2.9 {
        return Ok(());
    }

    let bl1 = accum_data.bond_lengths[bid1];
    let bl2 = accum_data.bond_lengths[bid2];
    let bl3 = accum_data.bond_lengths[bid3];
    let ba12 = accum_data.get_bond_angle(mol.bonds.len(), bid1, bid2);
    let ba23 = accum_data.get_bond_angle(mol.bonds.len(), bid2, bid3);
    validate_bond_angle(ba12, bid1, bid2, "set_in_ring_14_bounds")?;
    validate_bond_angle(ba23, bid2, bid3, "set_in_ring_14_bounds")?;

    let stype = get_atom_stereo(&mol.bonds[bid2], aid1, aid4);
    let mut prefer_cis = false;
    let mut prefer_trans = false;

    if ring_size <= 8
        && ahyb2 == Hybridization::Sp2
        && ahyb3 == Hybridization::Sp2
        && !matches!(stype, BondStereo::E | BondStereo::Trans)
    {
        if ring_info.num_bond_rings(BondId::new(bid2)) > 1 {
            if ring_info.num_bond_rings(BondId::new(bid1)) == 1
                && ring_info.num_bond_rings(BondId::new(bid3)) == 1
            {
                for br in ring_info.bond_rings() {
                    if br.contains(&BondId::new(bid1)) {
                        if br.contains(&BondId::new(bid3)) {
                            prefer_cis = true;
                        }
                        break;
                    }
                }
            }
        } else {
            prefer_cis = true;
        }
    } else if matches!(stype, BondStereo::Z | BondStereo::Cis) {
        prefer_cis = true;
    } else if matches!(stype, BondStereo::E | BondStereo::Trans) {
        prefer_trans = true;
    }

    let nb = mol.bonds.len();
    let kind = if prefer_cis {
        record_path_flag(&mut accum_data.cis_paths, path14_id(nb, bid1, bid2, bid3));
        record_path_flag(&mut accum_data.cis_paths, path14_id(nb, bid3, bid2, bid1));
        Path14Kind::Cis
    } else if prefer_trans {
        record_path_flag(&mut accum_data.trans_paths, path14_id(nb, bid1, bid2, bid3));
        record_path_flag(&mut accum_data.trans_paths, path14_id(nb, bid3, bid2, bid1));
        Path14Kind::Trans
    } else {
        Path14Kind::Other
    };

    accum_data.paths14.push(Path14Configuration {
        bid1,
        bid2,
        bid3,
        kind,
    });

    let (dl, du) = if prefer_cis {
        let dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
        (dl, dl + 2.0 * GEN_DIST_TOL)
    } else if prefer_trans {
        let dl = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
        (dl, dl + 2.0 * GEN_DIST_TOL)
    } else {
        let mut dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23);
        let mut du = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
        if du < dl {
            std::mem::swap(&mut dl, &mut du);
        }
        if (du - dl).abs() < DIST12_DELTA {
            dl -= GEN_DIST_TOL;
            du += GEN_DIST_TOL;
        }
        (dl, du)
    };

    accum_data.visited14_bounds[pid] = true;
    check_and_set_bounds(mmat, aid1, aid4, dl, du, false)
}
fn set_two_in_same_ring_14_bounds(
    mol: &TopologyBlock,
    bid1: usize,
    bid2: usize,
    bid3: usize,
    accum_data: &mut ComputedData,
    mmat: &mut BoundsMatrix,
    dmat: &[f64],
) -> Result<(), GraphBoundsError> {
    // RDKit❗✔️: void _setTwoInSameRing14Bounds(const ROMol &mol, const Bond *bnd1,
    // RDKit❗✔️:                                const Bond *bnd2, const Bond *bnd3,
    // RDKit❗✔️:                                ComputedData &accumData,
    // RDKit❗✔️:                                DistGeom::BoundsMatPtr mmat, double *dmat) {
    // RDKit❗✔️:   PRECONDITION(bnd1, "");
    // RDKit❗✔️:   PRECONDITION(bnd2, "");
    // RDKit❗✔️:   PRECONDITION(bnd3, "");
    // RDKit❗✔️:   unsigned int bid1, bid2, bid3;
    // RDKit❗✔️:   bid1 = bnd1->getIdx();
    // RDKit❗✔️:   bid2 = bnd2->getIdx();
    // RDKit❗✔️:   bid3 = bnd3->getIdx();
    // RDKit❗✔️:   const Atom *atm2 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid1, bid2));
    // RDKit❗✔️:   PRECONDITION(atm2, "");
    // RDKit❗✔️:   const Atom *atm3 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid2, bid3));
    // RDKit❗✔️:   PRECONDITION(atm3, "");
    // RDKit❗✔️:
    // RDKit❗✔️:   unsigned int aid1 = bnd1->getOtherAtomIdx(atm2->getIdx());
    // RDKit❗✔️:   unsigned int aid4 = bnd3->getOtherAtomIdx(atm3->getIdx());
    // RDKit❗✔️:
    // RDKit❗✔️:   const unsigned int pid =
    // RDKit❗✔️:       std::min(aid1, aid4) * mol.getNumAtoms() + std::max(aid1, aid4);
    // RDKit❗✔️:
    // RDKit❗✔️:   if (accumData.visitedBound(pid, DistType::DIST13)) {
    // RDKit❗✔️:     // if this is already a 1-3 or 1-2 distance; do not overwrite
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // check that this actually is a 1-4 contact:
    // RDKit❗✔️:   if (dmat[std::max(aid1, aid4) * mmat->numRows() + std::min(aid1, aid4)] <
    // RDKit❗✔️:       2.9) {
    // RDKit❗✔️:     // std::cerr<<"skip: "<<aid1<<"-"<<aid4<<" because
    // RDKit❗✔️:     // d="<<dmat[std::max(aid1,aid4)*mmat->numRows()+std::min(aid1,aid4)]<<std::endl;
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // when we have fused rings, it can happen that this isn't actually a 1-4
    // RDKit❗✔️:   // contact,
    // RDKit❗✔️:   // (this was the cause of sf.net bug 2835784) check that now:
    // RDKit❗✔️:   if (mol.getBondBetweenAtoms(aid1, atm3->getIdx()) ||
    // RDKit❗✔️:       mol.getBondBetweenAtoms(aid4, atm2->getIdx())) {
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   Atom::HybridizationType ahyb3 = atm3->getHybridization();
    // RDKit❗✔️:   Atom::HybridizationType ahyb2 = atm2->getHybridization();
    // RDKit❗✔️:
    // RDKit❗✔️:   double bl1 = accumData.bondLengths[bid1];
    // RDKit❗✔️:   double bl2 = accumData.bondLengths[bid2];
    // RDKit❗✔️:   double bl3 = accumData.bondLengths[bid3];
    // RDKit❗✔️:
    // RDKit❗✔️:   double ba12 = accumData.bondAngles->getVal(bid1, bid2);
    // RDKit❗✔️:   double ba23 = accumData.bondAngles->getVal(bid2, bid3);
    // RDKit❗✔️:   CHECK_INVARIANT(ba12 > 0.0, "");
    // RDKit❗✔️:   CHECK_INVARIANT(ba23 > 0.0, "");
    // RDKit❗✔️:   double dl, du;
    // RDKit❗✔️:   Path14Configuration path14;
    // RDKit❗✔️:   unsigned int nb = mol.getNumBonds();
    // RDKit❗✔️:
    // RDKit❗✔️:   path14.bid1 = bid1;
    // RDKit❗✔️:   path14.bid2 = bid2;
    // RDKit❗✔️:   path14.bid3 = bid3;
    // RDKit❗✔️:   if ((ahyb2 == Atom::SP2) && (ahyb3 == Atom::SP2)) {  // FIX: check for trans
    // RDKit❗✔️:     // here we will assume 180 degrees: basically flat ring with an external
    // RDKit❗✔️:     // substituent
    // RDKit❗✔️:     dl = RDGeom::compute14DistTrans(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:     du = dl;
    // RDKit❗✔️:     dl -= GEN_DIST_TOL;
    // RDKit❗✔️:     du += GEN_DIST_TOL;
    // RDKit❗✔️:     path14.type = Path14Configuration::TRANS;
    // RDKit❗✔️:     accumData.transPaths.insert(static_cast<unsigned long>(bid1) * nb * nb +
    // RDKit❗✔️:                                 bid2 * nb + bid3);
    // RDKit❗✔️:     accumData.transPaths.insert(static_cast<unsigned long>(bid3) * nb * nb +
    // RDKit❗✔️:                                 bid2 * nb + bid1);
    // RDKit❗✔️:
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     // here we will assume anything is possible
    // RDKit❗✔️:     dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:     du = RDGeom::compute14DistTrans(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:
    // RDKit❗✔️:     // in highly-strained situations these can get mixed up:
    // RDKit❗✔️:     if (du < dl) {
    // RDKit❗✔️:       double tmpD = dl;
    // RDKit❗✔️:       dl = du;
    // RDKit❗✔️:       du = tmpD;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (fabs(du - dl) < DIST12_DELTA) {
    // RDKit❗✔️:       dl -= GEN_DIST_TOL;
    // RDKit❗✔️:       du += GEN_DIST_TOL;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     path14.type = Path14Configuration::OTHER;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // std::cerr << "1: " << aid1 << "-" << aid4 << ": " << dl << " -> " << du
    // RDKit❗✔️:   //           << std::endl;
    // RDKit❗✔️:   _checkAndSetBounds(aid1, aid4, dl, du, mmat);
    // RDKit❗✔️:   accumData.paths14.push_back(path14);
    // RDKit❗✔️:   accumData.visited14Bounds.set(pid);
    // RDKit❗✔️: }

    let atm2 = bond_pair_shared_atom(mol, accum_data, bid1, bid2)?;
    let atm3 = bond_pair_shared_atom(mol, accum_data, bid2, bid3)?;
    let bnd1 = &mol.bonds[bid1];
    let bnd3 = &mol.bonds[bid3];
    let aid1 = if bnd1.begin().index() == atm2 {
        bnd1.end().index()
    } else {
        bnd1.begin().index()
    };
    let aid4 = if bnd3.begin().index() == atm3 {
        bnd3.end().index()
    } else {
        bnd3.begin().index()
    };
    let pid = aid1.min(aid4) * mol.atoms.len() + aid1.max(aid4);

    if accum_data.visited_bound(pid, DistType::Dist13) {
        return Ok(());
    }
    if dmat[aid1.max(aid4) * mmat.dimension() + aid1.min(aid4)] < 2.9 {
        return Ok(());
    }
    if bond_between_idx_simple(mol, aid1, atm3).is_some()
        || bond_between_idx_simple(mol, aid4, atm2).is_some()
    {
        return Ok(());
    }

    let ahyb2 = mol.atoms[atm2].hybridization();
    let ahyb3 = mol.atoms[atm3].hybridization();
    let bl1 = accum_data.bond_lengths[bid1];
    let bl2 = accum_data.bond_lengths[bid2];
    let bl3 = accum_data.bond_lengths[bid3];
    let ba12 = accum_data.get_bond_angle(mol.bonds.len(), bid1, bid2);
    let ba23 = accum_data.get_bond_angle(mol.bonds.len(), bid2, bid3);
    validate_bond_angle(ba12, bid1, bid2, "set_two_in_same_ring_14_bounds")?;
    validate_bond_angle(ba23, bid2, bid3, "set_two_in_same_ring_14_bounds")?;

    let nb = mol.bonds.len();
    let (dl, du, kind) = if ahyb2 == Hybridization::Sp2 && ahyb3 == Hybridization::Sp2 {
        record_path_flag(&mut accum_data.trans_paths, path14_id(nb, bid1, bid2, bid3));
        record_path_flag(&mut accum_data.trans_paths, path14_id(nb, bid3, bid2, bid1));
        let du = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
        (du - GEN_DIST_TOL, du + GEN_DIST_TOL, Path14Kind::Trans)
    } else {
        let mut dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23);
        let mut du = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
        if du < dl {
            std::mem::swap(&mut dl, &mut du);
        }
        if (du - dl).abs() < DIST12_DELTA {
            dl -= GEN_DIST_TOL;
            du += GEN_DIST_TOL;
        }
        (dl, du, Path14Kind::Other)
    };

    check_and_set_bounds(mmat, aid1, aid4, dl, du, false)?;
    accum_data.paths14.push(Path14Configuration {
        bid1,
        bid2,
        bid3,
        kind,
    });
    accum_data.visited14_bounds[pid] = true;
    Ok(())
}
fn set_two_in_diff_ring_14_bounds(
    mol: &TopologyBlock,
    bid1: usize,
    bid2: usize,
    bid3: usize,
    accum_data: &mut ComputedData,
    mmat: &mut BoundsMatrix,
    dmat: &[f64],
    ring_info: &RingInfo,
) -> Result<(), GraphBoundsError> {
    // RDKit❗✔️: void _setTwoInDiffRing14Bounds(const ROMol &mol, const Bond *bnd1,
    // RDKit❗✔️:                                const Bond *bnd2, const Bond *bnd3,
    // RDKit❗✔️:                                ComputedData &accumData,
    // RDKit❗✔️:                                DistGeom::BoundsMatPtr mmat, double *dmat) {
    // RDKit❗✔️:   // this turns out to be very similar to all bonds in the same ring
    // RDKit❗✔️:   // situation.
    // RDKit❗✔️:   // There is probably some fine tuning that can be done when the atoms a2
    // RDKit❗✔️:   // and a3 are not sp2 hybridized, but we will not worry about that now;
    // RDKit❗✔️:   // simple use 0-180 deg for non-sp2 cases.
    // RDKit❗✔️:   _setInRing14Bounds(mol, bnd1, bnd2, bnd3, accumData, mmat, dmat, 0);
    // RDKit❗✔️: }

    set_in_ring_14_bounds(mol, bid1, bid2, bid3, accum_data, mmat, dmat, 0, ring_info)
}
fn set_share_ring_bond_14_bounds(
    mol: &TopologyBlock,
    bid1: usize,
    bid2: usize,
    bid3: usize,
    accum_data: &mut ComputedData,
    mmat: &mut BoundsMatrix,
    dmat: &[f64],
    ring_info: &RingInfo,
) -> Result<(), GraphBoundsError> {
    // RDKit❗✔️: void _setShareRingBond14Bounds(const ROMol &mol, const Bond *bnd1,
    // RDKit❗✔️:                                const Bond *bnd2, const Bond *bnd3,
    // RDKit❗✔️:                                ComputedData &accumData,
    // RDKit❗✔️:                                DistGeom::BoundsMatPtr mmat, double *dmat) {
    // RDKit❗✔️:   // once this turns out to be similar to bonds in the same ring
    // RDKit❗✔️:   _setInRing14Bounds(mol, bnd1, bnd2, bnd3, accumData, mmat, dmat, 0);
    // RDKit❗✔️: }

    set_in_ring_14_bounds(mol, bid1, bid2, bid3, accum_data, mmat, dmat, 0, ring_info)
}
fn check_macrocycle_two_in_same_ring_amide_ester_14(
    mol: &TopologyBlock,
    bnd1_idx: usize,
    bnd3_idx: usize,
    atm1: usize,
    atm2: usize,
    atm3: usize,
    atm4: usize,
) -> bool {
    // RDKit❗✔️: bool _checkMacrocycleTwoInSameRingAmideEster14(
    // RDKit❗✔️:     const Bond *bnd1, const Bond *bnd3, const Atom *atm1, const Atom *atm2,
    // RDKit❗✔️:     const Atom *atm3, const Atom *atm4) {
    // RDKit❗✔️:   unsigned int a1Num = atm1->getAtomicNum();
    // RDKit❗✔️:   unsigned int a2Num = atm2->getAtomicNum();
    // RDKit❗✔️:   unsigned int a3Num = atm3->getAtomicNum();
    // RDKit❗✔️:   unsigned int a4Num = atm4->getAtomicNum();
    // RDKit❗✔️:
    // RDKit❗✔️:   return a1Num != 1 && a3Num == 6 && bnd3->getBondType() == Bond::DOUBLE &&
    // RDKit❗✔️:          (a4Num == 8 || a4Num == 7) && bnd1->getBondType() == Bond::SINGLE &&
    // RDKit❗✔️:          (a2Num == 8 || a2Num == 7);
    // RDKit❗✔️: }

    let a1_num = mol.atoms[atm1].atomic_number();
    let a2_num = mol.atoms[atm2].atomic_number();
    let a3_num = mol.atoms[atm3].atomic_number();
    let a4_num = mol.atoms[atm4].atomic_number();

    a1_num != 1
        && a3_num == 6
        && mol.bonds[bnd3_idx].order() == BondOrder::Double
        && (a4_num == 8 || a4_num == 7)
        && mol.bonds[bnd1_idx].order() == BondOrder::Single
        && (a2_num == 8 || a2_num == 7)
}
fn set_macrocycle_two_in_same_ring_14_bounds(
    mol: &TopologyBlock,
    bid1: usize,
    bid2: usize,
    bid3: usize,
    accum_data: &mut ComputedData,
    mmat: &mut BoundsMatrix,
    dmat: &[f64],
) -> Result<(), GraphBoundsError> {
    // RDKit❗✔️: void _setMacrocycleTwoInSameRing14Bounds(const ROMol &mol, const Bond *bnd1,
    // RDKit❗✔️:                                          const Bond *bnd2, const Bond *bnd3,
    // RDKit❗✔️:                                          ComputedData &accumData,
    // RDKit❗✔️:                                          DistGeom::BoundsMatPtr mmat,
    // RDKit❗✔️:                                          double *dmat) {
    // RDKit❗✔️:   PRECONDITION(bnd1, "");
    // RDKit❗✔️:   PRECONDITION(bnd2, "");
    // RDKit❗✔️:   PRECONDITION(bnd3, "");
    // RDKit❗✔️:   unsigned int bid1, bid2, bid3;
    // RDKit❗✔️:   bid1 = bnd1->getIdx();
    // RDKit❗✔️:   bid2 = bnd2->getIdx();
    // RDKit❗✔️:   bid3 = bnd3->getIdx();
    // RDKit❗✔️:   const Atom *atm2 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid1, bid2));
    // RDKit❗✔️:   PRECONDITION(atm2, "");
    // RDKit❗✔️:   const Atom *atm3 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid2, bid3));
    // RDKit❗✔️:   PRECONDITION(atm3, "");
    // RDKit❗✔️:
    // RDKit❗✔️:   unsigned int aid1 = bnd1->getOtherAtomIdx(atm2->getIdx());
    // RDKit❗✔️:   unsigned int aid4 = bnd3->getOtherAtomIdx(atm3->getIdx());
    // RDKit❗✔️:   const Atom *atm1 = mol.getAtomWithIdx(aid1);
    // RDKit❗✔️:   const Atom *atm4 = mol.getAtomWithIdx(aid4);
    // RDKit❗✔️:
    // RDKit❗✔️:   const unsigned int pid =
    // RDKit❗✔️:       std::min(aid1, aid4) * mol.getNumAtoms() + std::max(aid1, aid4);
    // RDKit❗✔️:
    // RDKit❗✔️:   if (accumData.visitedBound(pid, DistType::DIST13)) {
    // RDKit❗✔️:     // if this is already a 1-3 or 1-2 distance; do not overwrite
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // check that this actually is a 1-4 contact:
    // RDKit❗✔️:   if (dmat[std::max(aid1, aid4) * mmat->numRows() + std::min(aid1, aid4)] <
    // RDKit❗✔️:       2.9) {
    // RDKit❗✔️:     // std::cerr<<"skip: "<<aid1<<"-"<<aid4<<" because
    // RDKit❗✔️:     // d="<<dmat[std::max(aid1,aid4)*mmat->numRows()+std::min(aid1,aid4)]<<std::endl;
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // when we have fused rings, it can happen that this isn't actually a 1-4
    // RDKit❗✔️:   // contact,
    // RDKit❗✔️:   // (this was the cause of sf.net bug 2835784) check that now:
    // RDKit❗✔️:   if (mol.getBondBetweenAtoms(aid1, atm3->getIdx()) ||
    // RDKit❗✔️:       mol.getBondBetweenAtoms(aid4, atm2->getIdx())) {
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   double bl1 = accumData.bondLengths[bid1];
    // RDKit❗✔️:   double bl2 = accumData.bondLengths[bid2];
    // RDKit❗✔️:   double bl3 = accumData.bondLengths[bid3];
    // RDKit❗✔️:
    // RDKit❗✔️:   double ba12 = accumData.bondAngles->getVal(bid1, bid2);
    // RDKit❗✔️:   double ba23 = accumData.bondAngles->getVal(bid2, bid3);
    // RDKit❗✔️:   CHECK_INVARIANT(ba12 > 0.0, "");
    // RDKit❗✔️:   CHECK_INVARIANT(ba23 > 0.0, "");
    // RDKit❗✔️:   double dl, du;
    // RDKit❗✔️:   Path14Configuration path14;
    // RDKit❗✔️:   unsigned int nb = mol.getNumBonds();
    // RDKit❗✔️:
    // RDKit❗✔️:   path14.bid1 = bid1;
    // RDKit❗✔️:   path14.bid2 = bid2;
    // RDKit❗✔️:   path14.bid3 = bid3;
    // RDKit❗✔️:   if ((_checkMacrocycleTwoInSameRingAmideEster14(bnd1, bnd3, atm1, atm2, atm3,
    // RDKit❗✔️:                                                  atm4)) ||
    // RDKit❗✔️:       (_checkMacrocycleTwoInSameRingAmideEster14(bnd3, bnd1, atm4, atm3, atm2,
    // RDKit❗✔️:                                                  atm1))) {
    // RDKit❗✔️:     dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:     path14.type = Path14Configuration::CIS;
    // RDKit❗✔️:     accumData.cisPaths.insert(static_cast<unsigned long>(bid1) * nb * nb +
    // RDKit❗✔️:                               bid2 * nb + bid3);
    // RDKit❗✔️:     accumData.cisPaths.insert(static_cast<unsigned long>(bid3) * nb * nb +
    // RDKit❗✔️:                               bid2 * nb + bid1);
    // RDKit❗✔️:     du = dl;
    // RDKit❗✔️:     dl -= GEN_DIST_TOL;
    // RDKit❗✔️:     du += GEN_DIST_TOL;
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     // here we will assume anything is possible
    // RDKit❗✔️:     dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:     du = RDGeom::compute14DistTrans(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:
    // RDKit❗✔️:     // in highly-strained situations these can get mixed up:
    // RDKit❗✔️:     if (du < dl) {
    // RDKit❗✔️:       double tmpD = dl;
    // RDKit❗✔️:       dl = du;
    // RDKit❗✔️:       du = tmpD;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (fabs(du - dl) < DIST12_DELTA) {
    // RDKit❗✔️:       dl -= GEN_DIST_TOL;
    // RDKit❗✔️:       du += GEN_DIST_TOL;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     path14.type = Path14Configuration::OTHER;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   // std::cerr << "1: " << aid1 << "-" << aid4 << ": " << dl << " -> " << du
    // RDKit❗✔️:   //           << std::endl;
    // RDKit❗✔️:   _checkAndSetBounds(aid1, aid4, dl, du, mmat);
    // RDKit❗✔️:   accumData.paths14.push_back(path14);
    // RDKit❗✔️:   accumData.visited14Bounds.set(pid);
    // RDKit❗✔️: }

    let atm2 = bond_pair_shared_atom(mol, accum_data, bid1, bid2)?;
    let atm3 = bond_pair_shared_atom(mol, accum_data, bid2, bid3)?;
    let bnd1 = &mol.bonds[bid1];
    let bnd3 = &mol.bonds[bid3];
    let aid1 = if bnd1.begin().index() == atm2 {
        bnd1.end().index()
    } else {
        bnd1.begin().index()
    };
    let aid4 = if bnd3.begin().index() == atm3 {
        bnd3.end().index()
    } else {
        bnd3.begin().index()
    };
    let atm1 = aid1;
    let atm4 = aid4;
    let pid = aid1.min(aid4) * mol.atoms.len() + aid1.max(aid4);

    if accum_data.visited_bound(pid, DistType::Dist13) {
        return Ok(());
    }
    if dmat[aid1.max(aid4) * mmat.dimension() + aid1.min(aid4)] < 2.9 {
        return Ok(());
    }
    if bond_between_idx_simple(mol, aid1, atm3).is_some()
        || bond_between_idx_simple(mol, aid4, atm2).is_some()
    {
        return Ok(());
    }

    let bl1 = accum_data.bond_lengths[bid1];
    let bl2 = accum_data.bond_lengths[bid2];
    let bl3 = accum_data.bond_lengths[bid3];
    let ba12 = accum_data.get_bond_angle(mol.bonds.len(), bid1, bid2);
    let ba23 = accum_data.get_bond_angle(mol.bonds.len(), bid2, bid3);
    validate_bond_angle(
        ba12,
        bid1,
        bid2,
        "set_macrocycle_two_in_same_ring_14_bounds",
    )?;
    validate_bond_angle(
        ba23,
        bid2,
        bid3,
        "set_macrocycle_two_in_same_ring_14_bounds",
    )?;

    let nb = mol.bonds.len();
    let (dl, du, kind) = if check_macrocycle_two_in_same_ring_amide_ester_14(
        mol, bid1, bid3, atm1, atm2, atm3, atm4,
    ) || check_macrocycle_two_in_same_ring_amide_ester_14(
        mol, bid3, bid1, atm4, atm3, atm2, atm1,
    ) {
        let dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23);
        record_path_flag(&mut accum_data.cis_paths, path14_id(nb, bid1, bid2, bid3));
        record_path_flag(&mut accum_data.cis_paths, path14_id(nb, bid3, bid2, bid1));
        (dl - GEN_DIST_TOL, dl + GEN_DIST_TOL, Path14Kind::Cis)
    } else {
        let mut dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23);
        let mut du = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
        if du < dl {
            std::mem::swap(&mut dl, &mut du);
        }
        if (du - dl).abs() < DIST12_DELTA {
            dl -= GEN_DIST_TOL;
            du += GEN_DIST_TOL;
        }
        (dl, du, Path14Kind::Other)
    };

    check_and_set_bounds(mmat, aid1, aid4, dl, du, false)?;
    accum_data.paths14.push(Path14Configuration {
        bid1,
        bid2,
        bid3,
        kind,
    });
    accum_data.visited14_bounds[pid] = true;
    Ok(())
}
fn set_macrocycle_all_in_same_ring_14_bounds(
    mol: &TopologyBlock,
    valence: &ValenceAssignment,
    bid1: usize,
    bid2: usize,
    bid3: usize,
    accum_data: &mut ComputedData,
    mmat: &mut BoundsMatrix,
) -> Result<(), GraphBoundsError> {
    // RDKit❗✔️: void _setMacrocycleAllInSameRing14Bounds(const ROMol &mol, const Bond *bnd1,
    // RDKit❗✔️:                                          const Bond *bnd2, const Bond *bnd3,
    // RDKit❗✔️:                                          ComputedData &accumData,
    // RDKit❗✔️:                                          DistGeom::BoundsMatPtr mmat,
    // RDKit❗✔️:                                          double *) {
    // RDKit❗✔️:   // This is adapted from `_setChain14Bounds`, with changes on how trans amide
    // RDKit❗✔️:   // is handled
    // RDKit❗✔️:   PRECONDITION(bnd1, "");
    // RDKit❗✔️:   PRECONDITION(bnd2, "");
    // RDKit❗✔️:   PRECONDITION(bnd3, "");
    // RDKit❗✔️:   unsigned int bid1, bid2, bid3;
    // RDKit❗✔️:   bid1 = bnd1->getIdx();
    // RDKit❗✔️:   bid2 = bnd2->getIdx();
    // RDKit❗✔️:   bid3 = bnd3->getIdx();
    // RDKit❗✔️:   const Atom *atm2 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid1, bid2));
    // RDKit❗✔️:   PRECONDITION(atm2, "");
    // RDKit❗✔️:   const Atom *atm3 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid2, bid3));
    // RDKit❗✔️:   PRECONDITION(atm3, "");
    // RDKit❗✔️:
    // RDKit❗✔️:   unsigned int aid1 = bnd1->getOtherAtomIdx(atm2->getIdx());
    // RDKit❗✔️:   unsigned int aid4 = bnd3->getOtherAtomIdx(atm3->getIdx());
    // RDKit❗✔️:
    // RDKit❗✔️:   const unsigned int pid =
    // RDKit❗✔️:       std::min(aid1, aid4) * mol.getNumAtoms() + std::max(aid1, aid4);
    // RDKit❗✔️:
    // RDKit❗✔️:   if (accumData.visitedBound(pid, DistType::DIST13)) {
    // RDKit❗✔️:     // if this is already a 1-3 or 1-2 distance; do not overwrite
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   const Atom *atm1 = mol.getAtomWithIdx(aid1);
    // RDKit❗✔️:   const Atom *atm4 = mol.getAtomWithIdx(aid4);
    // RDKit❗✔️:
    // RDKit❗✔️:   double bl1 = accumData.bondLengths[bid1];
    // RDKit❗✔️:   double bl2 = accumData.bondLengths[bid2];
    // RDKit❗✔️:   double bl3 = accumData.bondLengths[bid3];
    // RDKit❗✔️:
    // RDKit❗✔️:   double ba12 = accumData.bondAngles->getVal(bid1, bid2);
    // RDKit❗✔️:   double ba23 = accumData.bondAngles->getVal(bid2, bid3);
    // RDKit❗✔️:   CHECK_INVARIANT(ba12 > 0.0, "");
    // RDKit❗✔️:   CHECK_INVARIANT(ba23 > 0.0, "");
    // RDKit❗✔️:   bool setTheBound = true;
    // RDKit❗✔️:   double dl = 0.0, du = 0.0;
    // RDKit❗✔️:
    // RDKit❗✔️:   // if the middle bond is double
    // RDKit❗✔️:   Path14Configuration path14;
    // RDKit❗✔️:   path14.bid1 = bid1;
    // RDKit❗✔️:   path14.bid2 = bid2;
    // RDKit❗✔️:   path14.bid3 = bid3;
    // RDKit❗✔️:   unsigned int nb = mol.getNumBonds();
    // RDKit❗✔️:   switch (bnd2->getBondType()) {
    // RDKit❗✔️:     case Bond::DOUBLE:
    // RDKit❗✔️:       // if any of the other bonds are double - the torsion angle is zero
    // RDKit❗✔️:       // this is CC=C=C situation
    // RDKit❗✔️:       if ((bnd1->getBondType() == Bond::DOUBLE) ||
    // RDKit❗✔️:           (bnd3->getBondType() == Bond::DOUBLE)) {
    // RDKit❗✔️:         dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
    // RDKit❗✔️:         du = dl + 2 * GEN_DIST_TOL;
    // RDKit❗✔️:         path14.type = Path14Configuration::CIS;
    // RDKit❗✔️:         accumData.cisPaths.insert(static_cast<unsigned long>(bid1) * nb * nb +
    // RDKit❗✔️:                                   bid2 * nb + bid3);
    // RDKit❗✔️:         accumData.cisPaths.insert(static_cast<unsigned long>(bid3) * nb * nb +
    // RDKit❗✔️:                                   bid2 * nb + bid1);
    // RDKit❗✔️:         // BOOST_LOG(rdDebugLog) << "Special 5 " << aid1 << " " << aid4 <<
    // RDKit❗✔️:         // "\n";
    // RDKit❗✔️:       } else if (bnd2->getStereo() > Bond::STEREOANY) {
    // RDKit❗✔️:         Bond::BondStereo stype = _getAtomStereo(bnd2, aid1, aid4);
    // RDKit❗✔️:         if (stype == Bond::STEREOZ || stype == Bond::STEREOCIS) {
    // RDKit❗✔️:           dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23) -
    // RDKit❗✔️:                GEN_DIST_TOL;
    // RDKit❗✔️:           du = dl + 2 * GEN_DIST_TOL;
    // RDKit❗✔️:           path14.type = Path14Configuration::CIS;
    // RDKit❗✔️:           // BOOST_LOG(rdDebugLog) << "Special 6 " <<  aid1 << " " << aid4 <<
    // RDKit❗✔️:           // "\n";
    // RDKit❗✔️:           accumData.cisPaths.insert(static_cast<unsigned long>(bid1) * nb * nb +
    // RDKit❗✔️:                                     bid2 * nb + bid3);
    // RDKit❗✔️:           accumData.cisPaths.insert(static_cast<unsigned long>(bid3) * nb * nb +
    // RDKit❗✔️:                                     bid2 * nb + bid1);
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           // BOOST_LOG(rdDebugLog) << "Special 7 " << aid1 << " " << aid4 <<
    // RDKit❗✔️:           // "\n";
    // RDKit❗✔️:           du = RDGeom::compute14DistTrans(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:           dl = du;
    // RDKit❗✔️:           dl -= GEN_DIST_TOL;
    // RDKit❗✔️:           du += GEN_DIST_TOL;
    // RDKit❗✔️:           path14.type = Path14Configuration::TRANS;
    // RDKit❗✔️:           accumData.transPaths.insert(
    // RDKit❗✔️:               static_cast<unsigned long>(bid1) * nb * nb + bid2 * nb + bid3);
    // RDKit❗✔️:           accumData.transPaths.insert(
    // RDKit❗✔️:               static_cast<unsigned long>(bid3) * nb * nb + bid2 * nb + bid1);
    // RDKit❗✔️:         }
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         // double bond with no stereo setting can be 0 or 180
    // RDKit❗✔️:         dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:         du = RDGeom::compute14DistTrans(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:         if (fabs(du - dl) < DIST12_DELTA) {
    // RDKit❗✔️:           dl -= GEN_DIST_TOL;
    // RDKit❗✔️:           du += GEN_DIST_TOL;
    // RDKit❗✔️:         }
    // RDKit❗✔️:         path14.type = Path14Configuration::OTHER;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     case Bond::SINGLE:
    // RDKit❗✔️:       if ((atm2->getAtomicNum() == 16) && (atm3->getAtomicNum() == 16) &&
    // RDKit❗✔️:           (atm2->getDegree() == 2) && (atm3->getDegree() == 2)) {
    // RDKit❗✔️:         // this is *S-S* situation
    // RDKit❗✔️:         dl = RDGeom::compute14Dist3D(bl1, bl2, bl3, ba12, ba23, M_PI / 2) -
    // RDKit❗✔️:              GEN_DIST_TOL;
    // RDKit❗✔️:         du = dl + 2 * GEN_DIST_TOL;
    // RDKit❗✔️:         path14.type = Path14Configuration::OTHER;
    // RDKit❗✔️:         // BOOST_LOG(rdDebugLog) << "Special 9 " << aid1 << " " << aid4 <<
    // RDKit❗✔️:         // "\n";
    // RDKit❗✔️:       } else if ((_checkMacrocycleAllInSameRingAmideEster14(
    // RDKit❗✔️:                      mol, bnd1, bnd3, atm1, atm2, atm3, atm4)) ||
    // RDKit❗✔️:                  (_checkMacrocycleAllInSameRingAmideEster14(
    // RDKit❗✔️:                      mol, bnd3, bnd1, atm4, atm3, atm2, atm1))) {
    // RDKit❗✔️:         dl = RDGeom::compute14DistTrans(bl1, bl2, bl3, ba12, ba23) +
    // RDKit❗✔️:              0.1;  // we saw that the currently defined max distance for trans
    // RDKit❗✔️:                    // is still a bit too short, thus we add an additional 0.1,
    // RDKit❗✔️:                    // which is the max that works without triangular smoothing
    // RDKit❗✔️:                    // error
    // RDKit❗✔️:         path14.type = Path14Configuration::TRANS;
    // RDKit❗✔️:         accumData.transPaths.insert(static_cast<unsigned long>(bid1) * nb * nb +
    // RDKit❗✔️:                                     bid2 * nb + bid3);
    // RDKit❗✔️:         accumData.transPaths.insert(static_cast<unsigned long>(bid3) * nb * nb +
    // RDKit❗✔️:                                     bid2 * nb + bid1);
    // RDKit❗✔️:
    // RDKit❗✔️:         du = dl;
    // RDKit❗✔️:         dl -= GEN_DIST_TOL;
    // RDKit❗✔️:         du += GEN_DIST_TOL;
    // RDKit❗✔️:
    // RDKit❗✔️:         // BOOST_LOG(rdDebugLog) << "  amide: " << aid1 << " " << aid4 << ": "
    // RDKit❗✔️:         // << dl << "->" << du << "\n";
    // RDKit❗✔️:       } else if ((_checkAmideEster15(mol, bnd1, bnd3, atm1, atm2, atm3,
    // RDKit❗✔️:                                      atm4)) ||
    // RDKit❗✔️:                  (_checkAmideEster15(mol, bnd3, bnd1, atm4, atm3, atm2,
    // RDKit❗✔️:                                      atm1))) {
    // RDKit❗✔️: #ifdef FORCE_TRANS_AMIDES
    // RDKit❗✔️:         // amide is trans, we're cis:
    // RDKit❗✔️:         dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:         path14.type = Path14Configuration::CIS;
    // RDKit❗✔️:         accumData.cisPaths.insert(static_cast<unsigned long>(bid1) * nb * nb +
    // RDKit❗✔️:                                   bid2 * nb + bid3);
    // RDKit❗✔️:         accumData.cisPaths.insert(static_cast<unsigned long>(bid3) * nb * nb +
    // RDKit❗✔️:                                   bid2 * nb + bid1);
    // RDKit❗✔️: #else
    // RDKit❗✔️:         // amide is cis, we're trans:
    // RDKit❗✔️:         if (atm2->getAtomicNum() == 7 && atm2->getDegree() == 3 &&
    // RDKit❗✔️:             atm1->getAtomicNum() == 1 && atm2->getTotalNumHs(true) == 1) {
    // RDKit❗✔️:           // secondary amide, this is the H
    // RDKit❗✔️:           setTheBound = false;
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           dl = RDGeom::compute14DistTrans(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:           path14.type = Path14Configuration::TRANS;
    // RDKit❗✔️:           accumData.transPaths.insert(
    // RDKit❗✔️:               static_cast<unsigned long>(bid1) * nb * nb + bid2 * nb + bid3);
    // RDKit❗✔️:           accumData.transPaths.insert(
    // RDKit❗✔️:               static_cast<unsigned long>(bid3) * nb * nb + bid2 * nb + bid1);
    // RDKit❗✔️:         }
    // RDKit❗✔️: #endif
    // RDKit❗✔️:         du = dl;
    // RDKit❗✔️:         dl -= GEN_DIST_TOL;
    // RDKit❗✔️:         du += GEN_DIST_TOL;
    // RDKit❗✔️:         // BOOST_LOG(rdDebugLog) << "    amide neighbor: " << aid1 << " " <<
    // RDKit❗✔️:         // aid4 << ": " << dl << "->" << du << "\n";
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:         du = RDGeom::compute14DistTrans(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:         path14.type = Path14Configuration::OTHER;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     default:
    // RDKit❗✔️:       // BOOST_LOG(rdDebugLog) << "Special 12 " << aid1 << " " << aid4 <<
    // RDKit❗✔️:       // "\n";
    // RDKit❗✔️:       dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:       du = RDGeom::compute14DistTrans(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:
    // RDKit❗✔️:       path14.type = Path14Configuration::OTHER;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   if (setTheBound) {
    // RDKit❗✔️:     if (fabs(du - dl) < DIST12_DELTA) {
    // RDKit❗✔️:       dl -= GEN_DIST_TOL;
    // RDKit❗✔️:       du += GEN_DIST_TOL;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     // std::cerr<<"2: "<<aid1<<"-"<<aid4<<std::endl;
    // RDKit❗✔️:
    // RDKit❗✔️:     // we only overwrite bounds if they are not 1-2 nor 1-3 distances
    // RDKit❗✔️:     _checkAndSetBounds(aid1, aid4, dl, du, mmat);
    // RDKit❗✔️:     accumData.paths14.push_back(path14);
    // RDKit❗✔️:     accumData.visited14Bounds.set(pid);
    // RDKit❗✔️:   }
    // RDKit❗✔️: }

    let atm2 = bond_pair_shared_atom(mol, accum_data, bid1, bid2)?;
    let atm3 = bond_pair_shared_atom(mol, accum_data, bid2, bid3)?;
    let bnd1 = &mol.bonds[bid1];
    let bnd2 = &mol.bonds[bid2];
    let bnd3 = &mol.bonds[bid3];
    let aid1 = if bnd1.begin().index() == atm2 {
        bnd1.end().index()
    } else {
        bnd1.begin().index()
    };
    let aid4 = if bnd3.begin().index() == atm3 {
        bnd3.end().index()
    } else {
        bnd3.begin().index()
    };
    let pid = aid1.min(aid4) * mol.atoms.len() + aid1.max(aid4);
    if accum_data.visited_bound(pid, DistType::Dist13) {
        return Ok(());
    }

    let atm1 = aid1;
    let atm4 = aid4;
    let bl1 = accum_data.bond_lengths[bid1];
    let bl2 = accum_data.bond_lengths[bid2];
    let bl3 = accum_data.bond_lengths[bid3];
    let ba12 = accum_data.get_bond_angle(mol.bonds.len(), bid1, bid2);
    let ba23 = accum_data.get_bond_angle(mol.bonds.len(), bid2, bid3);
    validate_bond_angle(
        ba12,
        bid1,
        bid2,
        "set_macrocycle_all_in_same_ring_14_bounds",
    )?;
    validate_bond_angle(
        ba23,
        bid2,
        bid3,
        "set_macrocycle_all_in_same_ring_14_bounds",
    )?;

    let mut set_the_bound = true;
    let nb = mol.bonds.len();
    let (mut dl, mut du, kind) = match bnd2.order() {
        BondOrder::Double => {
            if bnd1.order() == BondOrder::Double || bnd3.order() == BondOrder::Double {
                let dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
                record_path_flag(&mut accum_data.cis_paths, path14_id(nb, bid1, bid2, bid3));
                record_path_flag(&mut accum_data.cis_paths, path14_id(nb, bid3, bid2, bid1));
                (dl, dl + 2.0 * GEN_DIST_TOL, Path14Kind::Cis)
            } else if matches!(
                bnd2.stereo(),
                BondStereo::Z
                    | BondStereo::E
                    | BondStereo::Cis
                    | BondStereo::Trans
                    | BondStereo::AtropCw
                    | BondStereo::AtropCcw
            ) {
                let stype = get_atom_stereo(bnd2, aid1, aid4);
                if matches!(stype, BondStereo::Z | BondStereo::Cis) {
                    let dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
                    record_path_flag(&mut accum_data.cis_paths, path14_id(nb, bid1, bid2, bid3));
                    record_path_flag(&mut accum_data.cis_paths, path14_id(nb, bid3, bid2, bid1));
                    (dl, dl + 2.0 * GEN_DIST_TOL, Path14Kind::Cis)
                } else {
                    let du = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                    record_path_flag(&mut accum_data.trans_paths, path14_id(nb, bid1, bid2, bid3));
                    record_path_flag(&mut accum_data.trans_paths, path14_id(nb, bid3, bid2, bid1));
                    (du - GEN_DIST_TOL, du + GEN_DIST_TOL, Path14Kind::Trans)
                }
            } else {
                let mut dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                let mut du = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                if (du - dl).abs() < DIST12_DELTA {
                    dl -= GEN_DIST_TOL;
                    du += GEN_DIST_TOL;
                }
                (dl, du, Path14Kind::Other)
            }
        }
        BondOrder::Single => {
            if mol.atoms[atm2].atomic_number() == 16
                && mol.atoms[atm3].atomic_number() == 16
                && mol.adjacency.neighbors_of(atm2).len() == 2
                && mol.adjacency.neighbors_of(atm3).len() == 2
            {
                let dl = compute_14_dist_3d(bl1, bl2, bl3, ba12, ba23, std::f64::consts::PI / 2.0)
                    - GEN_DIST_TOL;
                (dl, dl + 2.0 * GEN_DIST_TOL, Path14Kind::Other)
            } else if check_macrocycle_all_in_same_ring_amide_ester_14(mol, atm1, atm2, atm3, atm4)
                || check_macrocycle_all_in_same_ring_amide_ester_14(mol, atm4, atm3, atm2, atm1)
            {
                let dl = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23) + 0.1;
                record_path_flag(&mut accum_data.trans_paths, path14_id(nb, bid1, bid2, bid3));
                record_path_flag(&mut accum_data.trans_paths, path14_id(nb, bid3, bid2, bid1));
                (dl - GEN_DIST_TOL, dl + GEN_DIST_TOL, Path14Kind::Trans)
            } else if check_amide_ester_15(mol, valence, bid1, bid3, atm2, atm3)?
                || check_amide_ester_15(mol, valence, bid3, bid1, atm3, atm2)?
            {
                let total_hs_atm2 = cosmolkit_core::total_hydrogen_count_from_validated(
                    mol,
                    valence,
                    AtomId::new(atm2),
                    true,
                )?;
                if mol.atoms[atm2].atomic_number() == 7
                    && mol.adjacency.neighbors_of(atm2).len() == 3
                    && mol.atoms[atm1].atomic_number() == 1
                    && total_hs_atm2 == 1
                {
                    set_the_bound = false;
                    (0.0, 0.0, Path14Kind::Other)
                } else {
                    let dl = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                    record_path_flag(&mut accum_data.trans_paths, path14_id(nb, bid1, bid2, bid3));
                    record_path_flag(&mut accum_data.trans_paths, path14_id(nb, bid3, bid2, bid1));
                    (dl - GEN_DIST_TOL, dl + GEN_DIST_TOL, Path14Kind::Trans)
                }
            } else {
                (
                    compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23),
                    compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23),
                    Path14Kind::Other,
                )
            }
        }
        _ => (
            compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23),
            compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23),
            Path14Kind::Other,
        ),
    };

    if set_the_bound {
        if (du - dl).abs() < DIST12_DELTA {
            dl -= GEN_DIST_TOL;
            du += GEN_DIST_TOL;
        }
        check_and_set_bounds(mmat, aid1, aid4, dl, du, false)?;
        accum_data.paths14.push(Path14Configuration {
            bid1,
            bid2,
            bid3,
            kind,
        });
        accum_data.visited14_bounds[pid] = true;
    }
    Ok(())
}
fn set_chain_14_bounds(
    mol: &TopologyBlock,
    valence: &ValenceAssignment,
    bid1: usize,
    bid2: usize,
    bid3: usize,
    accum_data: &mut ComputedData,
    mmat: &mut BoundsMatrix,
    force_trans_amides: bool,
) -> Result<(), GraphBoundsError> {
    // RDKit❗✔️: void _setChain14Bounds(const ROMol &mol, const Bond *bnd1, const Bond *bnd2,
    // RDKit❗✔️:                        const Bond *bnd3, ComputedData &accumData,
    // RDKit❗✔️:                        DistGeom::BoundsMatPtr mmat, double *,
    // RDKit❗✔️:                        bool forceTransAmides) {
    // RDKit❗✔️:   PRECONDITION(bnd1, "");
    // RDKit❗✔️:   PRECONDITION(bnd2, "");
    // RDKit❗✔️:   PRECONDITION(bnd3, "");
    // RDKit❗✔️:   unsigned int bid1, bid2, bid3;
    // RDKit❗✔️:   bid1 = bnd1->getIdx();
    // RDKit❗✔️:   bid2 = bnd2->getIdx();
    // RDKit❗✔️:   bid3 = bnd3->getIdx();
    // RDKit❗✔️:   const Atom *atm2 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid1, bid2));
    // RDKit❗✔️:   PRECONDITION(atm2, "");
    // RDKit❗✔️:   const Atom *atm3 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid2, bid3));
    // RDKit❗✔️:   PRECONDITION(atm3, "");
    // RDKit❗✔️:
    // RDKit❗✔️:   unsigned int aid1 = bnd1->getOtherAtomIdx(atm2->getIdx());
    // RDKit❗✔️:   unsigned int aid4 = bnd3->getOtherAtomIdx(atm3->getIdx());
    // RDKit❗✔️:   const Atom *atm1 = mol.getAtomWithIdx(aid1);
    // RDKit❗✔️:   const Atom *atm4 = mol.getAtomWithIdx(aid4);
    // RDKit❗✔️:
    // RDKit❗✔️:   const unsigned int pid =
    // RDKit❗✔️:       std::min(aid1, aid4) * mol.getNumAtoms() + std::max(aid1, aid4);
    // RDKit❗✔️:
    // RDKit❗✔️:   if (accumData.visitedBound(pid, DistType::DIST13)) {
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   double bl1 = accumData.bondLengths[bid1];
    // RDKit❗✔️:   double bl2 = accumData.bondLengths[bid2];
    // RDKit❗✔️:   double bl3 = accumData.bondLengths[bid3];
    // RDKit❗✔️:
    // RDKit❗✔️:   double ba12 = accumData.bondAngles->getVal(bid1, bid2);
    // RDKit❗✔️:   double ba23 = accumData.bondAngles->getVal(bid2, bid3);
    // RDKit❗✔️:   CHECK_INVARIANT(ba12 > 0.0, "");
    // RDKit❗✔️:   CHECK_INVARIANT(ba23 > 0.0, "");
    // RDKit❗✔️:   bool setTheBound = true;
    // RDKit❗✔️:   double dl = 0.0, du = 0.0;
    // RDKit❗✔️:
    // RDKit❗✔️:   // if the middle bond is double
    // RDKit❗✔️:   Path14Configuration path14;
    // RDKit❗✔️:   path14.bid1 = bid1;
    // RDKit❗✔️:   path14.bid2 = bid2;
    // RDKit❗✔️:   path14.bid3 = bid3;
    // RDKit❗✔️:   unsigned int nb = mol.getNumBonds();
    // RDKit❗✔️:   switch (bnd2->getBondType()) {
    // RDKit❗✔️:     case Bond::DOUBLE:
    // RDKit❗✔️:       // if any of the other bonds are double - the torsion angle is zero
    // RDKit❗✔️:       // this is CC=C=C situation
    // RDKit❗✔️:       if ((bnd1->getBondType() == Bond::DOUBLE) ||
    // RDKit❗✔️:           (bnd3->getBondType() == Bond::DOUBLE)) {
    // RDKit❗✔️:         dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
    // RDKit❗✔️:         du = dl + 2 * GEN_DIST_TOL;
    // RDKit❗✔️:         path14.type = Path14Configuration::CIS;
    // RDKit❗✔️:         accumData.cisPaths.insert(static_cast<unsigned long>(bid1) * nb * nb +
    // RDKit❗✔️:                                   bid2 * nb + bid3);
    // RDKit❗✔️:         accumData.cisPaths.insert(static_cast<unsigned long>(bid3) * nb * nb +
    // RDKit❗✔️:                                   bid2 * nb + bid1);
    // RDKit❗✔️:         // BOOST_LOG(rdDebugLog) << "Special 5 " << aid1 << " " << aid4 <<
    // RDKit❗✔️:         // "\n";
    // RDKit❗✔️:       } else if (bnd2->getStereo() > Bond::STEREOANY) {
    // RDKit❗✔️:         Bond::BondStereo stype = _getAtomStereo(bnd2, aid1, aid4);
    // RDKit❗✔️:         if (stype == Bond::STEREOZ || stype == Bond::STEREOCIS) {
    // RDKit❗✔️:           dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23) -
    // RDKit❗✔️:                GEN_DIST_TOL;
    // RDKit❗✔️:           du = dl + 2 * GEN_DIST_TOL;
    // RDKit❗✔️:           path14.type = Path14Configuration::CIS;
    // RDKit❗✔️:           // BOOST_LOG(rdDebugLog) << "Special 6 " <<  aid1 << " " << aid4
    // RDKit❗✔️:           // <<
    // RDKit❗✔️:           // "\n";
    // RDKit❗✔️:           accumData.cisPaths.insert(static_cast<unsigned long>(bid1) * nb * nb +
    // RDKit❗✔️:                                     bid2 * nb + bid3);
    // RDKit❗✔️:           accumData.cisPaths.insert(static_cast<unsigned long>(bid3) * nb * nb +
    // RDKit❗✔️:                                     bid2 * nb + bid1);
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           // BOOST_LOG(rdDebugLog) << "Special 7 " << aid1 << " " << aid4 <<
    // RDKit❗✔️:           // "\n";
    // RDKit❗✔️:           du = RDGeom::compute14DistTrans(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:           dl = du;
    // RDKit❗✔️:           dl -= GEN_DIST_TOL;
    // RDKit❗✔️:           du += GEN_DIST_TOL;
    // RDKit❗✔️:           path14.type = Path14Configuration::TRANS;
    // RDKit❗✔️:           accumData.transPaths.insert(
    // RDKit❗✔️:               static_cast<unsigned long>(bid1) * nb * nb + bid2 * nb + bid3);
    // RDKit❗✔️:           accumData.transPaths.insert(
    // RDKit❗✔️:               static_cast<unsigned long>(bid3) * nb * nb + bid2 * nb + bid1);
    // RDKit❗✔️:         }
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         // double bond with no stereo setting can be 0 or 180
    // RDKit❗✔️:         dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:         du = RDGeom::compute14DistTrans(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:         if (fabs(du - dl) < DIST12_DELTA) {
    // RDKit❗✔️:           dl -= GEN_DIST_TOL;
    // RDKit❗✔️:           du += GEN_DIST_TOL;
    // RDKit❗✔️:         }
    // RDKit❗✔️:         path14.type = Path14Configuration::OTHER;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     case Bond::SINGLE:
    // RDKit❗✔️:       if ((atm2->getAtomicNum() == 16) && (atm3->getAtomicNum() == 16) &&
    // RDKit❗✔️:           (atm2->getDegree() == 2) && (atm3->getDegree() == 2)) {
    // RDKit❗✔️:         // this is *S-S* situation
    // RDKit❗✔️:         dl = RDGeom::compute14Dist3D(bl1, bl2, bl3, ba12, ba23, M_PI / 2) -
    // RDKit❗✔️:              GEN_DIST_TOL;
    // RDKit❗✔️:         du = dl + 2 * GEN_DIST_TOL;
    // RDKit❗✔️:         path14.type = Path14Configuration::OTHER;
    // RDKit❗✔️:         // BOOST_LOG(rdDebugLog) << "Special 9 " << aid1 << " " << aid4 <<
    // RDKit❗✔️:         // "\n";
    // RDKit❗✔️:       } else if ((_checkAmideEster14(bnd1, bnd3, atm1, atm2, atm3, atm4)) ||
    // RDKit❗✔️:                  (_checkAmideEster14(bnd3, bnd1, atm4, atm3, atm2, atm1))) {
    // RDKit❗✔️:         // It's an amide or ester:
    // RDKit❗✔️:         //
    // RDKit❗✔️:         //        4    <- 4 is the O
    // RDKit❗✔️:         //        |    <- That's the double bond
    // RDKit❗✔️:         //    1   3
    // RDKit❗✔️:         //     \ / \                                         T.S.I.Left Blank
    // RDKit❗✔️:         //      2   5  <- 2 is an oxygen/nitrogen
    // RDKit❗✔️:         //
    // RDKit❗✔️:         // Here we set the distance between atoms 1 and 4,
    // RDKit❗✔️:         //  we'll handle atoms 1 and 5 below.
    // RDKit❗✔️:
    // RDKit❗✔️:         // fix for issue 251 - we were marking this as a cis configuration
    // RDKit❗✔️:         // earlier
    // RDKit❗✔️:         // -------------------------------------------------------
    // RDKit❗✔️:         // Issue284:
    // RDKit❗✔️:         //   As this code originally stood, we forced amide bonds to be trans.
    // RDKit❗✔️:         //   This is convenient a lot of the time for generating nice-looking
    // RDKit❗✔️:         //   structures, but is unfortunately totally bogus.  So here we'll
    // RDKit❗✔️:         //   allow the distance to roam from cis to trans and hope that the
    // RDKit❗✔️:         //   force field planarizes things later.
    // RDKit❗✔️:         //
    // RDKit❗✔️:         //   What we'd really like to be able to do is specify multiple
    // RDKit❗✔️:         //   possible ranges for the distances, but a single bounds matrix
    // RDKit❗✔️:         //   doesn't support this kind of fanciness.
    // RDKit❗✔️:         //
    // RDKit❗✔️:         if (forceTransAmides) {
    // RDKit❗✔️:           if ((atm1->getAtomicNum() == 1 && atm2->getAtomicNum() == 7 &&
    // RDKit❗✔️:                atm2->getDegree() == 3 && atm2->getTotalNumHs(true) == 1) ||
    // RDKit❗✔️:               (atm4->getAtomicNum() == 1 && atm3->getAtomicNum() == 7 &&
    // RDKit❗✔️:                atm3->getDegree() == 3 && atm3->getTotalNumHs(true) == 1)) {
    // RDKit❗✔️:             // secondary amide, this is the H, it should be trans to the O
    // RDKit❗✔️:             dl = RDGeom::compute14DistTrans(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:             path14.type = Path14Configuration::TRANS;
    // RDKit❗✔️:             accumData.transPaths.insert(
    // RDKit❗✔️:                 static_cast<unsigned long>(bid1) * nb * nb + bid2 * nb + bid3);
    // RDKit❗✔️:             accumData.transPaths.insert(
    // RDKit❗✔️:                 static_cast<unsigned long>(bid3) * nb * nb + bid2 * nb + bid1);
    // RDKit❗✔️:           } else {
    // RDKit❗✔️:             dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:             path14.type = Path14Configuration::CIS;
    // RDKit❗✔️:             accumData.cisPaths.insert(
    // RDKit❗✔️:                 static_cast<unsigned long>(bid1) * nb * nb + bid2 * nb + bid3);
    // RDKit❗✔️:             accumData.cisPaths.insert(
    // RDKit❗✔️:                 static_cast<unsigned long>(bid3) * nb * nb + bid2 * nb + bid1);
    // RDKit❗✔️:           }
    // RDKit❗✔️:           du = dl;
    // RDKit❗✔️:           dl -= GEN_DIST_TOL;
    // RDKit❗✔️:           du += GEN_DIST_TOL;
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:           du = RDGeom::compute14DistTrans(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:           path14.type = Path14Configuration::OTHER;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       } else if ((_checkAmideEster15(mol, bnd1, bnd3, atm1, atm2, atm3,
    // RDKit❗✔️:                                      atm4)) ||
    // RDKit❗✔️:                  (_checkAmideEster15(mol, bnd3, bnd1, atm4, atm3, atm2,
    // RDKit❗✔️:                                      atm1))) {
    // RDKit❗✔️:         // it's an amide or ester.
    // RDKit❗✔️:         //
    // RDKit❗✔️:         //        4    <- 4 is the O
    // RDKit❗✔️:         //        |    <- That's the double bond
    // RDKit❗✔️:         //    1   3
    // RDKit❗✔️:         //     \ / \                                          T.S.I.Left Blank
    // RDKit❗✔️:         //      2   5  <- 2 is oxygen or nitrogen
    // RDKit❗✔️:         //
    // RDKit❗✔️:         // we already set the 1-4 contact above, here we are doing 1-5
    // RDKit❗✔️:
    // RDKit❗✔️:         // If we're going to have a hope of getting good geometries
    // RDKit❗✔️:         // out of here we need to set some reasonably smart bounds between 1
    // RDKit❗✔️:         // and 5 (ref Issue355):
    // RDKit❗✔️:
    // RDKit❗✔️:         if (forceTransAmides) {
    // RDKit❗✔️:           if ((atm1->getAtomicNum() == 1 && atm2->getAtomicNum() == 7 &&
    // RDKit❗✔️:                atm2->getDegree() == 3 && atm2->getTotalNumHs(true) == 1) ||
    // RDKit❗✔️:               (atm4->getAtomicNum() == 1 && atm3->getAtomicNum() == 7 &&
    // RDKit❗✔️:                atm3->getDegree() == 3 && atm3->getTotalNumHs(true) == 1)) {
    // RDKit❗✔️:             // secondary amide, this is the H, it's cis to atom 5
    // RDKit❗✔️:             dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:             path14.type = Path14Configuration::CIS;
    // RDKit❗✔️:             accumData.cisPaths.insert(
    // RDKit❗✔️:                 static_cast<unsigned long>(bid1) * nb * nb + bid2 * nb + bid3);
    // RDKit❗✔️:             accumData.cisPaths.insert(
    // RDKit❗✔️:                 static_cast<unsigned long>(bid3) * nb * nb + bid2 * nb + bid1);
    // RDKit❗✔️:           } else {
    // RDKit❗✔️:             dl = RDGeom::compute14DistTrans(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:             path14.type = Path14Configuration::TRANS;
    // RDKit❗✔️:             accumData.transPaths.insert(
    // RDKit❗✔️:                 static_cast<unsigned long>(bid1) * nb * nb + bid2 * nb + bid3);
    // RDKit❗✔️:             accumData.transPaths.insert(
    // RDKit❗✔️:                 static_cast<unsigned long>(bid3) * nb * nb + bid2 * nb + bid1);
    // RDKit❗✔️:           }
    // RDKit❗✔️:           du = dl;
    // RDKit❗✔️:           dl -= GEN_DIST_TOL;
    // RDKit❗✔️:           du += GEN_DIST_TOL;
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:           du = RDGeom::compute14DistTrans(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:           path14.type = Path14Configuration::OTHER;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:         du = RDGeom::compute14DistTrans(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:         path14.type = Path14Configuration::OTHER;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     default:
    // RDKit❗✔️:       dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:       du = RDGeom::compute14DistTrans(bl1, bl2, bl3, ba12, ba23);
    // RDKit❗✔️:
    // RDKit❗✔️:       path14.type = Path14Configuration::OTHER;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (setTheBound) {
    // RDKit❗✔️:     if (fabs(du - dl) < DIST12_DELTA) {
    // RDKit❗✔️:       dl -= GEN_DIST_TOL;
    // RDKit❗✔️:       du += GEN_DIST_TOL;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     _checkAndSetBounds(aid1, aid4, dl, du, mmat);
    // RDKit❗✔️:     accumData.paths14.push_back(path14);
    // RDKit❗✔️:     accumData.visited14Bounds.set(pid);
    // RDKit❗✔️:   }
    // RDKit❗✔️: }

    let atm2 = bond_pair_shared_atom(mol, accum_data, bid1, bid2)?;
    let atm3 = bond_pair_shared_atom(mol, accum_data, bid2, bid3)?;
    let bnd1 = &mol.bonds[bid1];
    let bnd2 = &mol.bonds[bid2];
    let bnd3 = &mol.bonds[bid3];
    let aid1 = if bnd1.begin().index() == atm2 {
        bnd1.end().index()
    } else {
        bnd1.begin().index()
    };
    let aid4 = if bnd3.begin().index() == atm3 {
        bnd3.end().index()
    } else {
        bnd3.begin().index()
    };
    let pid = aid1.min(aid4) * mol.atoms.len() + aid1.max(aid4);
    if accum_data.visited_bound(pid, DistType::Dist13) {
        return Ok(());
    }

    let atm1 = aid1;
    let atm4 = aid4;
    let bl1 = accum_data.bond_lengths[bid1];
    let bl2 = accum_data.bond_lengths[bid2];
    let bl3 = accum_data.bond_lengths[bid3];
    let ba12 = accum_data.get_bond_angle(mol.bonds.len(), bid1, bid2);
    let ba23 = accum_data.get_bond_angle(mol.bonds.len(), bid2, bid3);
    validate_bond_angle(ba12, bid1, bid2, "set_chain_14_bounds")?;
    validate_bond_angle(ba23, bid2, bid3, "set_chain_14_bounds")?;

    let set_the_bound = true;
    let nb = mol.bonds.len();
    let (mut dl, mut du, kind) = match bnd2.order() {
        BondOrder::Double => {
            if bnd1.order() == BondOrder::Double || bnd3.order() == BondOrder::Double {
                let dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
                record_path_flag(&mut accum_data.cis_paths, path14_id(nb, bid1, bid2, bid3));
                record_path_flag(&mut accum_data.cis_paths, path14_id(nb, bid3, bid2, bid1));
                (dl, dl + 2.0 * GEN_DIST_TOL, Path14Kind::Cis)
            } else if matches!(
                bnd2.stereo(),
                BondStereo::Z
                    | BondStereo::E
                    | BondStereo::Cis
                    | BondStereo::Trans
                    | BondStereo::AtropCw
                    | BondStereo::AtropCcw
            ) {
                let stype = get_atom_stereo(bnd2, aid1, aid4);
                if matches!(stype, BondStereo::Z | BondStereo::Cis) {
                    let dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23) - GEN_DIST_TOL;
                    record_path_flag(&mut accum_data.cis_paths, path14_id(nb, bid1, bid2, bid3));
                    record_path_flag(&mut accum_data.cis_paths, path14_id(nb, bid3, bid2, bid1));
                    (dl, dl + 2.0 * GEN_DIST_TOL, Path14Kind::Cis)
                } else {
                    let du = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                    record_path_flag(&mut accum_data.trans_paths, path14_id(nb, bid1, bid2, bid3));
                    record_path_flag(&mut accum_data.trans_paths, path14_id(nb, bid3, bid2, bid1));
                    (du - GEN_DIST_TOL, du + GEN_DIST_TOL, Path14Kind::Trans)
                }
            } else {
                let mut dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                let mut du = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                if (du - dl).abs() < DIST12_DELTA {
                    dl -= GEN_DIST_TOL;
                    du += GEN_DIST_TOL;
                }
                (dl, du, Path14Kind::Other)
            }
        }
        BondOrder::Single => {
            if mol.atoms[atm2].atomic_number() == 16
                && mol.atoms[atm3].atomic_number() == 16
                && mol.adjacency.neighbors_of(atm2).len() == 2
                && mol.adjacency.neighbors_of(atm3).len() == 2
            {
                let dl = compute_14_dist_3d(bl1, bl2, bl3, ba12, ba23, std::f64::consts::PI / 2.0)
                    - GEN_DIST_TOL;
                (dl, dl + 2.0 * GEN_DIST_TOL, Path14Kind::Other)
            } else if check_amide_ester_14(mol, valence, bid1, bid3, atm2, atm3, atm4)?
                || check_amide_ester_14(mol, valence, bid3, bid1, atm3, atm2, atm1)?
            {
                if force_trans_amides {
                    let total_hs_atm2 = cosmolkit_core::total_hydrogen_count_from_validated(
                        mol,
                        valence,
                        AtomId::new(atm2),
                        true,
                    )?;
                    let total_hs_atm3 = cosmolkit_core::total_hydrogen_count_from_validated(
                        mol,
                        valence,
                        AtomId::new(atm3),
                        true,
                    )?;
                    let secondary_left = mol.atoms[atm1].atomic_number() == 1
                        && mol.atoms[atm2].atomic_number() == 7
                        && mol.adjacency.neighbors_of(atm2).len() == 3
                        && total_hs_atm2 == 1;
                    let secondary_right = mol.atoms[atm4].atomic_number() == 1
                        && mol.atoms[atm3].atomic_number() == 7
                        && mol.adjacency.neighbors_of(atm3).len() == 3
                        && total_hs_atm3 == 1;
                    if secondary_left || secondary_right {
                        let dl = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                        record_path_flag(
                            &mut accum_data.trans_paths,
                            path14_id(nb, bid1, bid2, bid3),
                        );
                        record_path_flag(
                            &mut accum_data.trans_paths,
                            path14_id(nb, bid3, bid2, bid1),
                        );
                        (dl - GEN_DIST_TOL, dl + GEN_DIST_TOL, Path14Kind::Trans)
                    } else {
                        let dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                        record_path_flag(
                            &mut accum_data.cis_paths,
                            path14_id(nb, bid1, bid2, bid3),
                        );
                        record_path_flag(
                            &mut accum_data.cis_paths,
                            path14_id(nb, bid3, bid2, bid1),
                        );
                        (dl - GEN_DIST_TOL, dl + GEN_DIST_TOL, Path14Kind::Cis)
                    }
                } else {
                    (
                        compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23),
                        compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23),
                        Path14Kind::Other,
                    )
                }
            } else if check_amide_ester_15(mol, valence, bid1, bid3, atm2, atm3)?
                || check_amide_ester_15(mol, valence, bid3, bid1, atm3, atm2)?
            {
                if force_trans_amides {
                    let total_hs_atm2 = cosmolkit_core::total_hydrogen_count_from_validated(
                        mol,
                        valence,
                        AtomId::new(atm2),
                        true,
                    )?;
                    let total_hs_atm3 = cosmolkit_core::total_hydrogen_count_from_validated(
                        mol,
                        valence,
                        AtomId::new(atm3),
                        true,
                    )?;
                    let secondary_left = mol.atoms[atm1].atomic_number() == 1
                        && mol.atoms[atm2].atomic_number() == 7
                        && mol.adjacency.neighbors_of(atm2).len() == 3
                        && total_hs_atm2 == 1;
                    let secondary_right = mol.atoms[atm4].atomic_number() == 1
                        && mol.atoms[atm3].atomic_number() == 7
                        && mol.adjacency.neighbors_of(atm3).len() == 3
                        && total_hs_atm3 == 1;
                    if secondary_left || secondary_right {
                        let dl = compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23);
                        record_path_flag(
                            &mut accum_data.cis_paths,
                            path14_id(nb, bid1, bid2, bid3),
                        );
                        record_path_flag(
                            &mut accum_data.cis_paths,
                            path14_id(nb, bid3, bid2, bid1),
                        );
                        (dl - GEN_DIST_TOL, dl + GEN_DIST_TOL, Path14Kind::Cis)
                    } else {
                        let dl = compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23);
                        record_path_flag(
                            &mut accum_data.trans_paths,
                            path14_id(nb, bid1, bid2, bid3),
                        );
                        record_path_flag(
                            &mut accum_data.trans_paths,
                            path14_id(nb, bid3, bid2, bid1),
                        );
                        (dl - GEN_DIST_TOL, dl + GEN_DIST_TOL, Path14Kind::Trans)
                    }
                } else {
                    (
                        compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23),
                        compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23),
                        Path14Kind::Other,
                    )
                }
            } else {
                (
                    compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23),
                    compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23),
                    Path14Kind::Other,
                )
            }
        }
        _ => (
            compute_14_dist_cis(bl1, bl2, bl3, ba12, ba23),
            compute_14_dist_trans(bl1, bl2, bl3, ba12, ba23),
            Path14Kind::Other,
        ),
    };

    if set_the_bound {
        if (du - dl).abs() < DIST12_DELTA {
            dl -= GEN_DIST_TOL;
            du += GEN_DIST_TOL;
        }
        check_and_set_bounds(mmat, aid1, aid4, dl, du, false)?;
        accum_data.paths14.push(Path14Configuration {
            bid1,
            bid2,
            bid3,
            kind,
        });
        accum_data.visited14_bounds[pid] = true;
    }
    Ok(())
}
fn check_h2_nx3_h1_ox2(
    mol: &TopologyBlock,
    valence: &ValenceAssignment,
    atom_idx: usize,
) -> Result<bool, GraphBoundsError> {
    // RDKit❗✔️: bool _checkH2NX3H1OX2(const Atom *atm) {
    // RDKit❗✔️:   if ((atm->getAtomicNum() == 6) && (atm->getTotalNumHs(true) == 2)) {
    // RDKit❗✔️:     // CH2
    // RDKit❗✔️:     return true;
    // RDKit❗✔️:   } else if ((atm->getAtomicNum() == 8) && (atm->getTotalNumHs(true) == 0)) {
    // RDKit❗✔️:     // OX2
    // RDKit❗✔️:     return true;
    // RDKit❗✔️:   } else if ((atm->getAtomicNum() == 7) && (atm->getDegree() == 3) &&
    // RDKit❗✔️:              (atm->getTotalNumHs(true) == 1)) {
    // RDKit❗✔️:     // FIX: assuming hydrogen is not in the graph
    // RDKit❗✔️:     // this is NX3H1 situation
    // RDKit❗✔️:     return true;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return false;
    // RDKit❗✔️: }

    let atom = &mol.atoms[atom_idx];
    let total_hs = cosmolkit_core::total_hydrogen_count_from_validated(
        mol,
        valence,
        AtomId::new(atom_idx),
        true,
    )?;
    Ok((atom.atomic_number() == 6 && total_hs == 2)
        || (atom.atomic_number() == 8 && total_hs == 0)
        || (atom.atomic_number() == 7
            && mol.adjacency.neighbors_of(atom_idx).len() == 3
            && total_hs == 1))
}
fn check_nh_ch_ch_nh(
    mol: &TopologyBlock,
    valence: &ValenceAssignment,
    atm1: usize,
    atm2: usize,
    atm3: usize,
    atm4: usize,
) -> Result<bool, GraphBoundsError> {
    // RDKit❗✔️: bool _checkNhChChNh(const Atom *atm1, const Atom *atm2, const Atom *atm3,
    // RDKit❗✔️:                     const Atom *atm4) {
    // RDKit❗✔️:   // checking for [!#1]~$ch!@$ch~[!#1], where ch = [CH2,NX3H1,OX2] situation
    // RDKit❗✔️:   if ((atm1->getAtomicNum() != 1) && (atm4->getAtomicNum() != 1)) {
    // RDKit❗✔️:     // end atom not hydrogens
    // RDKit❗✔️:     if ((_checkH2NX3H1OX2(atm2)) && (_checkH2NX3H1OX2(atm3))) {
    // RDKit❗✔️:       return true;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return false;
    // RDKit❗✔️: }

    Ok(mol.atoms[atm1].atomic_number() != 1
        && mol.atoms[atm4].atomic_number() != 1
        && check_h2_nx3_h1_ox2(mol, valence, atm2)?
        && check_h2_nx3_h1_ox2(mol, valence, atm3)?)
}
fn check_amide_ester_14(
    mol: &TopologyBlock,
    valence: &ValenceAssignment,
    bnd1_idx: usize,
    bnd3_idx: usize,
    atm2: usize,
    atm3: usize,
    atm4: usize,
) -> Result<bool, GraphBoundsError> {
    // RDKit❗✔️: bool _checkAmideEster14(const Bond *bnd1, const Bond *bnd3, const Atom *,
    // RDKit❗✔️:                         const Atom *atm2, const Atom *atm3, const Atom *atm4) {
    // RDKit❗✔️:   unsigned int a2Num = atm2->getAtomicNum();
    // RDKit❗✔️:   unsigned int a3Num = atm3->getAtomicNum();
    // RDKit❗✔️:   unsigned int a4Num = atm4->getAtomicNum();
    // RDKit❗✔️:   // std::cerr << " -> " << atm1->getIdx() << "-" << atm2->getIdx() << "-"
    // RDKit❗✔️:   //           << atm3->getIdx() << "-" << atm4->getIdx()
    // RDKit❗✔️:   //           << " bonds: " << bnd1->getIdx() << "," << bnd3->getIdx()
    // RDKit❗✔️:   //           << std::endl;
    // RDKit❗✔️:   // std::cerr << "   " << a1Num << " " << a3Num << " " <<
    // RDKit❗✔️:   // bnd3->getBondType()
    // RDKit❗✔️:   //           << " " << a4Num << " " << bnd1->getBondType() << " " << a2Num
    // RDKit❗✔️:   //           << " "
    // RDKit❗✔️:   //           << atm2->getTotalNumHs(true) << std::endl;
    // RDKit❗✔️:   if (a3Num == 6 && bnd3->getBondType() == Bond::DOUBLE &&
    // RDKit❗✔️:       (a4Num == 8 || a4Num == 7) && bnd1->getBondType() == Bond::SINGLE &&
    // RDKit❗✔️:       (a2Num == 8 || (a2Num == 7 && atm2->getTotalNumHs(true) == 1))) {
    // RDKit❗✔️:     // std::cerr << " yes!" << std::endl;
    // RDKit❗✔️:     return true;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   // std::cerr << " no!" << std::endl;
    // RDKit❗✔️:   return false;
    // RDKit❗✔️: }

    let bnd1 = &mol.bonds[bnd1_idx];
    let bnd3 = &mol.bonds[bnd3_idx];
    let a2_num = mol.atoms[atm2].atomic_number();
    let a3_num = mol.atoms[atm3].atomic_number();
    let a4_num = mol.atoms[atm4].atomic_number();
    let total_hs_atm2 =
        cosmolkit_core::total_hydrogen_count_from_validated(mol, valence, AtomId::new(atm2), true)?;

    Ok(a3_num == 6
        && bnd3.order() == BondOrder::Double
        && (a4_num == 8 || a4_num == 7)
        && bnd1.order() == BondOrder::Single
        && (a2_num == 8 || (a2_num == 7 && total_hs_atm2 == 1)))
}
fn check_macrocycle_all_in_same_ring_amide_ester_14(
    mol: &TopologyBlock,
    atm1: usize,
    atm2: usize,
    atm3: usize,
    atm4: usize,
) -> bool {
    // RDKit❗✔️: bool _checkMacrocycleAllInSameRingAmideEster14(const ROMol &mol, const Bond *,
    // RDKit❗✔️:                                                const Bond *, const Atom *atm1,
    // RDKit❗✔️:                                                const Atom *atm2,
    // RDKit❗✔️:                                                const Atom *atm3,
    // RDKit❗✔️:                                                const Atom *atm4) {
    // RDKit❗✔️:   //   This is a re-write of `_checkAmideEster14` with more explicit logic
    // RDKit❗✔️:   //   on the checks It is interesting that we find with this function we
    // RDKit❗✔️:   //   get better macrocycle sampling than `_checkAmideEster14`
    // RDKit❗✔️:   unsigned int a2Num = atm2->getAtomicNum();
    // RDKit❗✔️:   unsigned int a3Num = atm3->getAtomicNum();
    // RDKit❗✔️:
    // RDKit❗✔️:   if (a3Num != 6) {
    // RDKit❗✔️:     return false;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   if (a2Num == 7 || a2Num == 8) {
    // RDKit❗✔️:     if (mol.getAtomDegree(atm2) == 3 && mol.getAtomDegree(atm3) == 3) {
    // RDKit❗✔️:       for (auto nbrIdx :
    // RDKit❗✔️:            boost::make_iterator_range(mol.getAtomNeighbors(atm2))) {
    // RDKit❗✔️:         if (nbrIdx != atm1->getIdx() && nbrIdx != atm3->getIdx()) {
    // RDKit❗✔️:           const auto &res = mol.getAtomWithIdx(nbrIdx);
    // RDKit❗✔️:           const auto &resbnd = mol.getBondBetweenAtoms(atm2->getIdx(), nbrIdx);
    // RDKit❗✔️:           if ((res->getAtomicNum() != 6 &&
    // RDKit❗✔️:                res->getAtomicNum() != 1) ||  // check is (methylated)amide
    // RDKit❗✔️:               resbnd->getBondType() != Bond::SINGLE) {
    // RDKit❗✔️:             return false;
    // RDKit❗✔️:           }
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:
    // RDKit❗✔️:       for (auto nbrIdx :
    // RDKit❗✔️:            boost::make_iterator_range(mol.getAtomNeighbors(atm3))) {
    // RDKit❗✔️:         if (nbrIdx != atm2->getIdx() && nbrIdx != atm4->getIdx()) {
    // RDKit❗✔️:           const auto &res = mol.getAtomWithIdx(nbrIdx);
    // RDKit❗✔️:           const auto &resbnd = mol.getBondBetweenAtoms(atm3->getIdx(), nbrIdx);
    // RDKit❗✔️:           if (res->getAtomicNum() != 8 ||  // check for the carbonyl oxygen
    // RDKit❗✔️:               resbnd->getBondType() != Bond::DOUBLE) {
    // RDKit❗✔️:             return false;
    // RDKit❗✔️:           }
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:
    // RDKit❗✔️:       return true;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return false;
    // RDKit❗✔️: }

    let a2_num = mol.atoms[atm2].atomic_number();
    let a3_num = mol.atoms[atm3].atomic_number();
    if a3_num != 6 {
        return false;
    }
    if a2_num != 7 && a2_num != 8 {
        return false;
    }
    if mol.adjacency.neighbors_of(atm2).len() != 3 || mol.adjacency.neighbors_of(atm3).len() != 3 {
        return false;
    }

    for neighbor in mol.adjacency.neighbors_of(atm2) {
        let nbr_idx = neighbor.atom_index;
        if nbr_idx != atm1 && nbr_idx != atm3 {
            let res = &mol.atoms[nbr_idx];
            let res_bnd = &mol.bonds[neighbor.bond.index()];
            if (res.atomic_number() != 6 && res.atomic_number() != 1)
                || res_bnd.order() != BondOrder::Single
            {
                return false;
            }
            break;
        }
    }

    for neighbor in mol.adjacency.neighbors_of(atm3) {
        let nbr_idx = neighbor.atom_index;
        if nbr_idx != atm2 && nbr_idx != atm4 {
            let res = &mol.atoms[nbr_idx];
            let res_bnd = &mol.bonds[neighbor.bond.index()];
            if res.atomic_number() != 8 || res_bnd.order() != BondOrder::Double {
                return false;
            }
            break;
        }
    }

    true
}
fn is_carbonyl(mol: &TopologyBlock, atom_idx: usize) -> bool {
    // RDKit❗✔️: bool _isCarbonyl(const ROMol &mol, const Atom *at) {
    // RDKit❗✔️:   PRECONDITION(at, "bad atom");
    // RDKit❗✔️:   if (at->getAtomicNum() == 6 && at->getDegree() > 2) {
    // RDKit❗✔️:     for (const auto nbr : mol.atomNeighbors(at)) {
    // RDKit❗✔️:       unsigned int atNum = nbr->getAtomicNum();
    // RDKit❗✔️:       if ((atNum == 8 || atNum == 7) &&
    // RDKit❗✔️:           mol.getBondBetweenAtoms(at->getIdx(), nbr->getIdx())->getBondType() ==
    // RDKit❗✔️:               Bond::DOUBLE) {
    // RDKit❗✔️:         return true;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return false;
    // RDKit❗✔️: }

    let atom = &mol.atoms[atom_idx];
    if atom.atomic_number() != 6 || mol.adjacency.neighbors_of(atom_idx).len() <= 2 {
        return false;
    }
    mol.adjacency.neighbors_of(atom_idx).iter().any(|neighbor| {
        let at_num = mol.atoms[neighbor.atom_index].atomic_number();
        (at_num == 8 || at_num == 7)
            && mol.bonds[neighbor.bond.index()].order() == BondOrder::Double
    })
}
fn check_amide_ester_15(
    mol: &TopologyBlock,
    valence: &ValenceAssignment,
    bnd1_idx: usize,
    bnd3_idx: usize,
    atm2: usize,
    atm3: usize,
) -> Result<bool, GraphBoundsError> {
    // RDKit❗✔️: bool _checkAmideEster15(const ROMol &mol, const Bond *bnd1, const Bond *bnd3,
    // RDKit❗✔️:                         const Atom *, const Atom *atm2, const Atom *atm3,
    // RDKit❗✔️:                         const Atom *) {
    // RDKit❗✔️:   unsigned int a2Num = atm2->getAtomicNum();
    // RDKit❗✔️:   if ((a2Num == 8) || ((a2Num == 7) && (atm2->getTotalNumHs(true) == 1))) {
    // RDKit❗✔️:     if ((bnd1->getBondType() == Bond::SINGLE)) {
    // RDKit❗✔️:       if ((atm3->getAtomicNum() == 6) &&
    // RDKit❗✔️:           (bnd3->getBondType() == Bond::SINGLE) && _isCarbonyl(mol, atm3)) {
    // RDKit❗✔️:         return true;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return false;
    // RDKit❗✔️: }

    let a2_num = mol.atoms[atm2].atomic_number();
    let total_hs_atm2 =
        cosmolkit_core::total_hydrogen_count_from_validated(mol, valence, AtomId::new(atm2), true)?;
    Ok(((a2_num == 8) || (a2_num == 7 && total_hs_atm2 == 1))
        && mol.bonds[bnd1_idx].order() == BondOrder::Single
        && mol.atoms[atm3].atomic_number() == 6
        && mol.bonds[bnd3_idx].order() == BondOrder::Single
        && is_carbonyl(mol, atm3))
}
fn set_14_bounds(
    mol: &TopologyBlock,
    valence: &ValenceAssignment,
    mmat: &mut BoundsMatrix,
    accum_data: &mut ComputedData,
    dmat: &[f64],
    use_macrocycle_14config: bool,
    force_trans_amides: bool,
    rinfo: &RingInfo,
) -> Result<(), GraphBoundsError> {
    // RDKit❗✔️: void set14Bounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat,
    // RDKit❗✔️:                  ComputedData &accumData, double *distMatrix,
    // RDKit❗✔️:                  bool useMacrocycle14config, bool forceTransAmides) {
    // RDKit❗✔️:   unsigned int npt = mmat->numRows();
    // RDKit❗✔️:   CHECK_INVARIANT(npt == mol.getNumAtoms(), "Wrong size metric matrix");
    // RDKit❗✔️:   // this is 2.6 million bonds, so it's extremly unlikely to ever occur, but
    // RDKit❗✔️:   // we might as well check:
    // RDKit❗✔️:   const size_t MAX_NUM_BONDS = static_cast<size_t>(
    // RDKit❗✔️:       std::pow(std::numeric_limits<std::uint64_t>::max(), 1. / 3));
    // RDKit❗✔️:   if (mol.getNumBonds() >= MAX_NUM_BONDS) {
    // RDKit❗✔️:     throw ValueErrorException(
    // RDKit❗✔️:         "Too many bonds in the molecule, cannot compute 1-4 bounds");
    // RDKit❗✔️:   }
    // RDKit❗✔️:   const auto rinfo = mol.getRingInfo();  // FIX: make sure we have ring info
    // RDKit❗✔️:   CHECK_INVARIANT(rinfo, "");
    // RDKit❗✔️:   const auto &bondRings = rinfo->bondRings();
    // RDKit❗✔️:
    // RDKit❗✔️:   std::unordered_set<unsigned int> bidIsMacrocycle;
    // RDKit❗✔️:   std::unordered_set<std::uint64_t> ringBondPairs;
    // RDKit❗✔️:   std::unordered_set<std::uint64_t> donePaths;
    // RDKit❗✔️:   std::uint64_t nb = mol.getNumBonds();
    // RDKit❗✔️:   // first we will deal with 1-4 atoms that belong to the same ring
    // RDKit❗✔️:   for (const auto &bring : bondRings) {
    // RDKit❗✔️:     const auto rSize = bring.size();
    // RDKit❗✔️:     if (rSize < 3) {
    // RDKit❗✔️:       continue;  // rings with less than 3 bonds are not useful
    // RDKit❗✔️:     }
    // RDKit❗✔️:     auto bid1 = bring[rSize - 1];
    // RDKit❗✔️:     for (auto i = 0u; i < rSize; i++) {
    // RDKit❗✔️:       auto bid2 = bring[i];
    // RDKit❗✔️:       auto bid3 = bring[(i + 1) % rSize];
    // RDKit❗✔️:       auto pid1 = bid1 * nb + bid2;
    // RDKit❗✔️:       auto pid2 = bid2 * nb + bid1;
    // RDKit❗✔️:       auto id1 = bid1 * nb * nb + bid2 * nb + bid3;
    // RDKit❗✔️:       auto id2 = bid3 * nb * nb + bid2 * nb + bid1;
    // RDKit❗✔️:
    // RDKit❗✔️:       ringBondPairs.insert(pid1);
    // RDKit❗✔️:       ringBondPairs.insert(pid2);
    // RDKit❗✔️:       donePaths.insert(id1);
    // RDKit❗✔️:       donePaths.insert(id2);
    // RDKit❗✔️:
    // RDKit❗✔️:       if (rSize > 5) {
    // RDKit❗✔️:         if (useMacrocycle14config && rSize >= minMacrocycleRingSize) {
    // RDKit❗✔️:           _setMacrocycleAllInSameRing14Bounds(
    // RDKit❗✔️:               mol, mol.getBondWithIdx(bid1), mol.getBondWithIdx(bid2),
    // RDKit❗✔️:               mol.getBondWithIdx(bid3), accumData, mmat, distMatrix);
    // RDKit❗✔️:           bidIsMacrocycle.insert(bid2);
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           _setInRing14Bounds(mol, mol.getBondWithIdx(bid1),
    // RDKit❗✔️:                              mol.getBondWithIdx(bid2), mol.getBondWithIdx(bid3),
    // RDKit❗✔️:                              accumData, mmat, distMatrix, rSize);
    // RDKit❗✔️:         }
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         _record14Path(mol, bid1, bid2, bid3, accumData);
    // RDKit❗✔️:       }
    // RDKit❗✔️:
    // RDKit❗✔️:       bid1 = bid2;
    // RDKit❗✔️:     }  // loop over bonds in the ring
    // RDKit❗✔️:   }  // end of all rings
    // RDKit❗✔️:   for (const auto bond : mol.bonds()) {
    // RDKit❗✔️:     auto bid2 = bond->getIdx();
    // RDKit❗✔️:     auto aid2 = bond->getBeginAtomIdx();
    // RDKit❗✔️:     auto aid3 = bond->getEndAtomIdx();
    // RDKit❗✔️:     for (const auto bnd1 : mol.atomBonds(mol.getAtomWithIdx(aid2))) {
    // RDKit❗✔️:       auto bid1 = bnd1->getIdx();
    // RDKit❗✔️:       if (bid1 != bid2) {
    // RDKit❗✔️:         for (const auto bnd3 : mol.atomBonds(mol.getAtomWithIdx(aid3))) {
    // RDKit❗✔️:           auto bid3 = bnd3->getIdx();
    // RDKit❗✔️:           if (bid3 != bid2) {
    // RDKit❗✔️:             auto id1 = bid1 * nb * nb + bid2 * nb + bid3;
    // RDKit❗✔️:             auto id2 = bid3 * nb * nb + bid2 * nb + bid1;
    // RDKit❗✔️:             if (donePaths.find(id1) == donePaths.end() &&
    // RDKit❗✔️:                 donePaths.find(id2) == donePaths.end()) {
    // RDKit❗✔️:               // we haven't dealt with this path before
    // RDKit❗✔️:               auto pid1 = bid1 * nb + bid2;
    // RDKit❗✔️:               auto pid2 = bid2 * nb + bid1;
    // RDKit❗✔️:               auto pid3 = bid2 * nb + bid3;
    // RDKit❗✔️:               auto pid4 = bid3 * nb + bid2;
    // RDKit❗✔️:
    // RDKit❗✔️:               if (ringBondPairs.find(pid1) != ringBondPairs.end() ||
    // RDKit❗✔️:                   ringBondPairs.find(pid2) != ringBondPairs.end() ||
    // RDKit❗✔️:                   ringBondPairs.find(pid3) != ringBondPairs.end() ||
    // RDKit❗✔️:                   ringBondPairs.find(pid4) != ringBondPairs.end()) {
    // RDKit❗✔️:                 // either (bid1, bid2) or (bid2, bid3) are in the
    // RDKit❗✔️:                 // same ring (note all three cannot be in the same
    // RDKit❗✔️:                 // ring; we dealt with that before)
    // RDKit❗✔️:                 if (useMacrocycle14config &&
    // RDKit❗✔️:                     bidIsMacrocycle.find(bid2) != bidIsMacrocycle.end()) {
    // RDKit❗✔️:                   _setMacrocycleTwoInSameRing14Bounds(
    // RDKit❗✔️:                       mol, bnd1, bond, bnd3, accumData, mmat, distMatrix);
    // RDKit❗✔️:                 } else {
    // RDKit❗✔️:                   _setTwoInSameRing14Bounds(mol, bnd1, bond, bnd3, accumData,
    // RDKit❗✔️:                                             mmat, distMatrix);
    // RDKit❗✔️:                 }
    // RDKit❗✔️:               } else if (((rinfo->numBondRings(bid1) > 0) &&
    // RDKit❗✔️:                           (rinfo->numBondRings(bid2) > 0)) ||
    // RDKit❗✔️:                          ((rinfo->numBondRings(bid2) > 0) &&
    // RDKit❗✔️:                           (rinfo->numBondRings(bid3) > 0))) {
    // RDKit❗✔️:                 // (bid1, bid2) or (bid2, bid3) are ring bonds but
    // RDKit❗✔️:                 // belong to different rings.  Note that the third
    // RDKit❗✔️:                 // bond will not belong to either of these two
    // RDKit❗✔️:                 // rings (if it does, we would have taken care of
    // RDKit❗✔️:                 // it in the previous if block); i.e. if bid1 and
    // RDKit❗✔️:                 // bid2 are ring bonds that belong to ring r1 and
    // RDKit❗✔️:                 // r2, then bid3 is either an external bond or
    // RDKit❗✔️:                 // belongs to a third ring r3.
    // RDKit❗✔️:                 _setTwoInDiffRing14Bounds(mol, bnd1, bond, bnd3, accumData,
    // RDKit❗✔️:                                           mmat, distMatrix);
    // RDKit❗✔️:               } else if (rinfo->numBondRings(bid2) > 0) {
    // RDKit❗✔️:                 // the middle bond is a ring bond and the other
    // RDKit❗✔️:                 // two do not belong to the same ring or are
    // RDKit❗✔️:                 // non-ring bonds
    // RDKit❗✔️:
    // RDKit❗✔️:                 _setShareRingBond14Bounds(mol, bnd1, bond, bnd3, accumData,
    // RDKit❗✔️:                                           mmat, distMatrix);
    // RDKit❗✔️:               } else {
    // RDKit❗✔️:                 // middle bond not a ring
    // RDKit❗✔️:                 _setChain14Bounds(mol, bnd1, bond, bnd3, accumData, mmat,
    // RDKit❗✔️:                                   distMatrix, forceTransAmides);
    // RDKit❗✔️:               }
    // RDKit❗✔️:             }
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }

    if dmat.len()
        != mol
            .atoms
            .len()
            .checked_mul(mol.atoms.len())
            .ok_or(GraphBoundsError::Input("atom pair count overflow"))?
    {
        return Err(GraphBoundsError::Input(
            "Wrong size topological distance matrix",
        ));
    }

    let npt = mmat.dimension();
    if npt != mol.atoms.len() {
        return Err(GraphBoundsError::Input("Wrong size metric matrix"));
    }
    let max_num_bonds = (u64::MAX as f64).powf(1.0 / 3.0) as usize;
    if mol.bonds.len() >= max_num_bonds {
        return Err(GraphBoundsError::Input(
            "Too many bonds in the molecule, cannot compute 1-4 bounds",
        ));
    }
    let bond_rings = rinfo.bond_rings();

    let mut bid_is_macrocycle: HashSet<usize> = HashSet::new();
    let mut ring_bond_pairs: HashSet<u64> = HashSet::new();
    let mut done_paths: HashSet<u64> = HashSet::new();
    let nb = mol.bonds.len() as u64;

    for bring in bond_rings {
        let r_size = bring.len();
        if r_size < 3 {
            continue;
        }
        let mut bid1 = bring[r_size - 1].index();
        for i in 0..r_size {
            let bid2 = bring[i].index();
            let bid3 = bring[(i + 1) % r_size].index();
            let pid1 = bid1 as u64 * nb + bid2 as u64;
            let pid2 = bid2 as u64 * nb + bid1 as u64;
            let id1 = bid1 as u64 * nb * nb + bid2 as u64 * nb + bid3 as u64;
            let id2 = bid3 as u64 * nb * nb + bid2 as u64 * nb + bid1 as u64;

            ring_bond_pairs.insert(pid1);
            ring_bond_pairs.insert(pid2);
            done_paths.insert(id1);
            done_paths.insert(id2);

            if r_size > 5 {
                if use_macrocycle_14config && r_size >= MIN_MACROCYCLE_RING_SIZE {
                    set_macrocycle_all_in_same_ring_14_bounds(
                        mol, valence, bid1, bid2, bid3, accum_data, mmat,
                    )?;
                    bid_is_macrocycle.insert(bid2);
                } else {
                    set_in_ring_14_bounds(
                        mol, bid1, bid2, bid3, accum_data, mmat, dmat, r_size, rinfo,
                    )?;
                }
            } else {
                record_14_path(mol, bid1, bid2, bid3, accum_data)?;
            }
            bid1 = bid2;
        }
    }

    for bond in &mol.bonds {
        let bid2 = bond.id().index();
        let aid2 = bond.begin().index();
        let aid3 = bond.end().index();
        for nbr1 in mol.adjacency.neighbors_of(aid2) {
            let bid1 = nbr1.bond.index();
            if bid1 == bid2 {
                continue;
            }
            for nbr3 in mol.adjacency.neighbors_of(aid3) {
                let bid3 = nbr3.bond.index();
                if bid3 == bid2 {
                    continue;
                }
                let id1 = bid1 as u64 * nb * nb + bid2 as u64 * nb + bid3 as u64;
                let id2 = bid3 as u64 * nb * nb + bid2 as u64 * nb + bid1 as u64;
                if done_paths.contains(&id1) || done_paths.contains(&id2) {
                    continue;
                }

                let pid1 = bid1 as u64 * nb + bid2 as u64;
                let pid2 = bid2 as u64 * nb + bid1 as u64;
                let pid3 = bid2 as u64 * nb + bid3 as u64;
                let pid4 = bid3 as u64 * nb + bid2 as u64;

                if ring_bond_pairs.contains(&pid1)
                    || ring_bond_pairs.contains(&pid2)
                    || ring_bond_pairs.contains(&pid3)
                    || ring_bond_pairs.contains(&pid4)
                {
                    if use_macrocycle_14config && bid_is_macrocycle.contains(&bid2) {
                        set_macrocycle_two_in_same_ring_14_bounds(
                            mol, bid1, bid2, bid3, accum_data, mmat, dmat,
                        )?;
                    } else {
                        set_two_in_same_ring_14_bounds(
                            mol, bid1, bid2, bid3, accum_data, mmat, dmat,
                        )?;
                    }
                } else if (rinfo.num_bond_rings(BondId::new(bid1)) > 0
                    && rinfo.num_bond_rings(BondId::new(bid2)) > 0)
                    || (rinfo.num_bond_rings(BondId::new(bid2)) > 0
                        && rinfo.num_bond_rings(BondId::new(bid3)) > 0)
                {
                    set_two_in_diff_ring_14_bounds(
                        mol, bid1, bid2, bid3, accum_data, mmat, dmat, rinfo,
                    )?;
                } else if rinfo.num_bond_rings(BondId::new(bid2)) > 0 {
                    set_share_ring_bond_14_bounds(
                        mol, bid1, bid2, bid3, accum_data, mmat, dmat, rinfo,
                    )?;
                } else {
                    set_chain_14_bounds(
                        mol,
                        valence,
                        bid1,
                        bid2,
                        bid3,
                        accum_data,
                        mmat,
                        force_trans_amides,
                    )?;
                }
            }
        }
    }
    Ok(())
}

const DIST15_TOL: f64 = 0.08;
const VDW_SCALE_15: f64 = 0.7;
const H_BOND_LENGTH: f64 = 1.8;
fn compute_15_dist_cis_cis(
    d1: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    ang12: f64,
    ang23: f64,
    ang34: f64,
) -> f64 {
    // RDKit❗✔️: double _compute15DistsCisCis(double d1, double d2, double d3, double d4,
    // RDKit❗✔️:                              double ang12, double ang23, double ang34) {
    // RDKit❗✔️:   double dx14 = d2 - d3 * cos(ang23) - d1 * cos(ang12);
    // RDKit❗✔️:   double dy14 = d3 * sin(ang23) - d1 * sin(ang12);
    // RDKit❗✔️:   double d14 = sqrt(dx14 * dx14 + dy14 * dy14);
    // RDKit❗✔️:   double cval = (d3 - d2 * cos(ang23) + d1 * cos(ang12 + ang23)) / d14;
    // RDKit❗✔️:   if (cval > 1.0) {
    // RDKit❗✔️:     cval = 1.0;
    // RDKit❗✔️:   } else if (cval < -1.0) {
    // RDKit❗✔️:     cval = -1.0;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   double ang143 = acos(cval);
    // RDKit❗✔️:   double ang145 = ang34 - ang143;
    // RDKit❗✔️:   double res = RDGeom::compute13Dist(d14, d4, ang145);
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }

    let dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy14 = d3 * ang23.sin() - d1 * ang12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let mut cval = ((d3 - d2 * ang23.cos() + d1 * (ang12 + ang23).cos()) / d14);
    if cval > 1.0 {
        cval = 1.0;
    } else if cval < -1.0 {
        cval = -1.0;
    }
    let ang143 = cval.acos();
    let ang145 = ang34 - ang143;
    compute_13_dist(d14, d4, ang145)
}
fn compute_15_dist_cis_trans(
    d1: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    ang12: f64,
    ang23: f64,
    ang34: f64,
) -> f64 {
    // RDKit❗✔️: double _compute15DistsCisTrans(double d1, double d2, double d3, double d4,
    // RDKit❗✔️:                                double ang12, double ang23, double ang34) {
    // RDKit❗✔️:   double dx14 = d2 - d3 * cos(ang23) - d1 * cos(ang12);
    // RDKit❗✔️:   double dy14 = d3 * sin(ang23) - d1 * sin(ang12);
    // RDKit❗✔️:   double d14 = sqrt(dx14 * dx14 + dy14 * dy14);
    // RDKit❗✔️:   double cval = (d3 - d2 * cos(ang23) + d1 * cos(ang12 + ang23)) / d14;
    // RDKit❗✔️:   if (cval > 1.0) {
    // RDKit❗✔️:     cval = 1.0;
    // RDKit❗✔️:   } else if (cval < -1.0) {
    // RDKit❗✔️:     cval = -1.0;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   double ang143 = acos(cval);
    // RDKit❗✔️:   double ang145 = ang34 + ang143;
    // RDKit❗✔️:   return RDGeom::compute13Dist(d14, d4, ang145);
    // RDKit❗✔️: }

    let dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy14 = d3 * ang23.sin() - d1 * ang12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let mut cval = ((d3 - d2 * ang23.cos() + d1 * (ang12 + ang23).cos()) / d14);
    if cval > 1.0 {
        cval = 1.0;
    } else if cval < -1.0 {
        cval = -1.0;
    }
    let ang143 = cval.acos();
    let ang145 = ang34 + ang143;
    compute_13_dist(d14, d4, ang145)
}
fn compute_15_dist_trans_trans(
    d1: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    ang12: f64,
    ang23: f64,
    ang34: f64,
) -> f64 {
    // RDKit❗✔️: double _compute15DistsTransTrans(double d1, double d2, double d3, double d4,
    // RDKit❗✔️:                                  double ang12, double ang23, double ang34) {
    // RDKit❗✔️:   double dx14 = d2 - d3 * cos(ang23) - d1 * cos(ang12);
    // RDKit❗✔️:   double dy14 = d3 * sin(ang23) + d1 * sin(ang12);
    // RDKit❗✔️:   double d14 = sqrt(dx14 * dx14 + dy14 * dy14);
    // RDKit❗✔️:   double cval = (d3 - d2 * cos(ang23) + d1 * cos(ang12 - ang23)) / d14;
    // RDKit❗✔️:   if (cval > 1.0) {
    // RDKit❗✔️:     cval = 1.0;
    // RDKit❗✔️:   } else if (cval < -1.0) {
    // RDKit❗✔️:     cval = -1.0;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   double ang143 = acos(cval);
    // RDKit❗✔️:   double ang145 = ang34 + ang143;
    // RDKit❗✔️:   return RDGeom::compute13Dist(d14, d4, ang145);
    // RDKit❗✔️: }

    let dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy14 = d3 * ang23.sin() + d1 * ang12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let mut cval = ((d3 - d2 * ang23.cos() + d1 * (ang12 - ang23).cos()) / d14);
    if cval > 1.0 {
        cval = 1.0;
    } else if cval < -1.0 {
        cval = -1.0;
    }
    let ang143 = cval.acos();
    let ang145 = ang34 + ang143;
    compute_13_dist(d14, d4, ang145)
}
fn compute_15_dist_trans_cis(
    d1: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    ang12: f64,
    ang23: f64,
    ang34: f64,
) -> f64 {
    // RDKit❗✔️: double _compute15DistsTransCis(double d1, double d2, double d3, double d4,
    // RDKit❗✔️:                                double ang12, double ang23, double ang34) {
    // RDKit❗✔️:   double dx14 = d2 - d3 * cos(ang23) - d1 * cos(ang12);
    // RDKit❗✔️:   double dy14 = d3 * sin(ang23) + d1 * sin(ang12);
    // RDKit❗✔️:   double d14 = sqrt(dx14 * dx14 + dy14 * dy14);
    // RDKit❗✔️:
    // RDKit❗✔️:   double cval = (d3 - d2 * cos(ang23) + d1 * cos(ang12 - ang23)) / d14;
    // RDKit❗✔️:   if (cval > 1.0) {
    // RDKit❗✔️:     cval = 1.0;
    // RDKit❗✔️:   } else if (cval < -1.0) {
    // RDKit❗✔️:     cval = -1.0;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   double ang143 = acos(cval);
    // RDKit❗✔️:   double ang145 = ang34 - ang143;
    // RDKit❗✔️:   return RDGeom::compute13Dist(d14, d4, ang145);
    // RDKit❗✔️: }

    let dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy14 = d3 * ang23.sin() + d1 * ang12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let mut cval = ((d3 - d2 * ang23.cos() + d1 * (ang12 - ang23).cos()) / d14);
    if cval > 1.0 {
        cval = 1.0;
    } else if cval < -1.0 {
        cval = -1.0;
    }
    let ang143 = cval.acos();
    let ang145 = ang34 - ang143;
    compute_13_dist(d14, d4, ang145)
}
fn set_15_bounds(
    mol: &TopologyBlock,
    mmat: &mut BoundsMatrix,
    accum_data: &mut ComputedData,
    dmat: &[f64],
) -> Result<(), GraphBoundsError> {
    // RDKit❗✔️: void set15Bounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat,
    // RDKit❗✔️:                  ComputedData &accumData, double *distMatrix) {
    // RDKit❗✔️:   PATH14_VECT_CI pti;
    // RDKit❗✔️:   unsigned int bid1, bid2, bid3, type;
    // RDKit❗✔️:   for (pti = accumData.paths14.begin(); pti != accumData.paths14.end(); pti++) {
    // RDKit❗✔️:     bid1 = pti->bid1;
    // RDKit❗✔️:     bid2 = pti->bid2;
    // RDKit❗✔️:     bid3 = pti->bid3;
    // RDKit❗✔️:     type = pti->type;
    // RDKit❗✔️:     // 15 distances going one way with with 14 paths
    // RDKit❗✔️:     _set15BoundsHelper(mol, bid1, bid2, bid3, type, accumData, mmat,
    // RDKit❗✔️:                        distMatrix);
    // RDKit❗✔️:     // going the other way - reverse the 14 path
    // RDKit❗✔️:     _set15BoundsHelper(mol, bid3, bid2, bid1, type, accumData, mmat,
    // RDKit❗✔️:                        distMatrix);
    // RDKit❗✔️:   }
    // RDKit❗✔️: }

    let nb = mol.bonds.len();
    let na = mol.atoms.len();
    for path_idx in 0..accum_data.paths14.len() {
        let path = accum_data.paths14[path_idx];
        set_15_bounds_helper(
            mol, mmat, accum_data, dmat, nb, na, path.bid1, path.bid2, path.bid3, path.kind,
        )?;
        set_15_bounds_helper(
            mol, mmat, accum_data, dmat, nb, na, path.bid3, path.bid2, path.bid1, path.kind,
        )?;
    }
    Ok(())
}
fn set_15_bounds_helper(
    mol: &TopologyBlock,
    mmat: &mut BoundsMatrix,
    accum_data: &mut ComputedData,
    dmat: &[f64],
    nb: usize,
    na: usize,
    bid1: usize,
    bid2: usize,
    bid3: usize,
    kind: Path14Kind,
) -> Result<(), GraphBoundsError> {
    // RDKit❗✔️: void _set15BoundsHelper(const ROMol &mol, unsigned int bid1, unsigned int bid2,
    // RDKit❗✔️:                         unsigned int bid3, unsigned int type,
    // RDKit❗✔️:                         ComputedData &accumData, DistGeom::BoundsMatPtr mmat,
    // RDKit❗✔️:                         double *dmat) {
    // RDKit❗✔️:   unsigned int i, aid1, aid2, aid3, aid4, aid5;
    // RDKit❗✔️:   double d1, d2, d3, d4, ang12, ang23, ang34, du, dl, vw1, vw5;
    // RDKit❗✔️:   unsigned int nb = mol.getNumBonds();
    // RDKit❗✔️:   unsigned int na = mol.getNumAtoms();
    // RDKit❗✔️:
    // RDKit❗✔️:   aid2 = accumData.bondAdj->getVal(bid1, bid2);
    // RDKit❗✔️:   aid1 = mol.getBondWithIdx(bid1)->getOtherAtomIdx(aid2);
    // RDKit❗✔️:   aid3 = accumData.bondAdj->getVal(bid2, bid3);
    // RDKit❗✔️:   aid4 = mol.getBondWithIdx(bid3)->getOtherAtomIdx(aid3);
    // RDKit❗✔️:   d1 = accumData.bondLengths[bid1];
    // RDKit❗✔️:   d2 = accumData.bondLengths[bid2];
    // RDKit❗✔️:   d3 = accumData.bondLengths[bid3];
    // RDKit❗✔️:   ang12 = accumData.bondAngles->getVal(bid1, bid2);
    // RDKit❗✔️:   ang23 = accumData.bondAngles->getVal(bid2, bid3);
    // RDKit❗✔️:   for (i = 0; i < nb; i++) {
    // RDKit❗✔️:     du = -1.0;
    // RDKit❗✔️:     dl = 0.0;
    // RDKit❗✔️:     if (accumData.bondAdj->getVal(bid3, i) == static_cast<int>(aid4)) {
    // RDKit❗✔️:       aid5 = mol.getBondWithIdx(i)->getOtherAtomIdx(aid4);
    // RDKit❗✔️:       // make sure we did not com back to the first atom in the path -
    // RDKit❗✔️:       // possible
    // RDKit❗✔️:       // with 4 membered rings
    // RDKit❗✔️:       // this is a fix for Issue 244
    // RDKit❗✔️:
    // RDKit❗✔️:       const unsigned int pid = std::min(aid1, aid5) * na + std::max(aid1, aid5);
    // RDKit❗✔️:
    // RDKit❗✔️:       if (accumData.visitedBound(pid, DistType::DIST14)) {
    // RDKit❗✔️:         return;
    // RDKit❗✔️:       }
    // RDKit❗✔️:
    // RDKit❗✔️:       // check that this actually is a 1-5 contact:
    // RDKit❗✔️:       if (dmat[std::max(aid1, aid5) * mmat->numRows() + std::min(aid1, aid5)] <
    // RDKit❗✔️:           3.9) {
    // RDKit❗✔️:         // std::cerr<<"skip: "<<aid1<<"-"<<aid5<<" because
    // RDKit❗✔️:         // d="<<dmat[std::max(aid1,aid5)*mmat->numRows()+std::min(aid1,aid5)]<<std::endl;
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:
    // RDKit❗✔️:       if (aid1 != aid5) {  // FIX: do we need this
    // RDKit❗✔️:         unsigned int pid1 = aid1 * na + aid5;
    // RDKit❗✔️:         unsigned int pid2 = aid5 * na + aid1;
    // RDKit❗✔️:         if ((mmat->getLowerBound(aid1, aid5) < DIST12_DELTA) ||
    // RDKit❗✔️:             (accumData.set15Atoms[pid1]) || (accumData.set15Atoms[pid2])) {
    // RDKit❗✔️:           d4 = accumData.bondLengths[i];
    // RDKit❗✔️:           ang34 = accumData.bondAngles->getVal(bid3, i);
    // RDKit❗✔️:           unsigned long pathId =
    // RDKit❗✔️:               static_cast<unsigned long>(bid2) * nb * nb + (bid3)*nb + i;
    // RDKit❗✔️:           if (type == 0) {
    // RDKit❗✔️:             if (accumData.cisPaths.find(pathId) != accumData.cisPaths.end()) {
    // RDKit❗✔️:               dl = _compute15DistsCisCis(d1, d2, d3, d4, ang12, ang23, ang34);
    // RDKit❗✔️:               du = dl + DIST15_TOL;
    // RDKit❗✔️:               dl -= DIST15_TOL;
    // RDKit❗✔️:             } else if (accumData.transPaths.find(pathId) !=
    // RDKit❗✔️:                        accumData.transPaths.end()) {
    // RDKit❗✔️:               dl = _compute15DistsCisTrans(d1, d2, d3, d4, ang12, ang23, ang34);
    // RDKit❗✔️:               du = dl + DIST15_TOL;
    // RDKit❗✔️:               dl -= DIST15_TOL;
    // RDKit❗✔️:             } else {
    // RDKit❗✔️:               dl = _compute15DistsCisCis(d1, d2, d3, d4, ang12, ang23, ang34) -
    // RDKit❗✔️:                    DIST15_TOL;
    // RDKit❗✔️:               du =
    // RDKit❗✔️:                   _compute15DistsCisTrans(d1, d2, d3, d4, ang12, ang23, ang34) +
    // RDKit❗✔️:                   DIST15_TOL;
    // RDKit❗✔️:             }
    // RDKit❗✔️:
    // RDKit❗✔️:           } else if (type == 1) {
    // RDKit❗✔️:             if (accumData.cisPaths.find(pathId) != accumData.cisPaths.end()) {
    // RDKit❗✔️:               dl = _compute15DistsTransCis(d1, d2, d3, d4, ang12, ang23, ang34);
    // RDKit❗✔️:               du = dl + DIST15_TOL;
    // RDKit❗✔️:               dl -= DIST15_TOL;
    // RDKit❗✔️:             } else if (accumData.transPaths.find(pathId) !=
    // RDKit❗✔️:                        accumData.transPaths.end()) {
    // RDKit❗✔️:               dl = _compute15DistsTransTrans(d1, d2, d3, d4, ang12, ang23,
    // RDKit❗✔️:                                              ang34);
    // RDKit❗✔️:               du = dl + DIST15_TOL;
    // RDKit❗✔️:               dl -= DIST15_TOL;
    // RDKit❗✔️:             } else {
    // RDKit❗✔️:               dl =
    // RDKit❗✔️:                   _compute15DistsTransCis(d1, d2, d3, d4, ang12, ang23, ang34) -
    // RDKit❗✔️:                   DIST15_TOL;
    // RDKit❗✔️:               du = _compute15DistsTransTrans(d1, d2, d3, d4, ang12, ang23,
    // RDKit❗✔️:                                              ang34) +
    // RDKit❗✔️:                    DIST15_TOL;
    // RDKit❗✔️:             }
    // RDKit❗✔️:           } else {
    // RDKit❗✔️:             if (accumData.cisPaths.find(pathId) != accumData.cisPaths.end()) {
    // RDKit❗✔️:               dl = _compute15DistsCisCis(d4, d3, d2, d1, ang34, ang23, ang12) -
    // RDKit❗✔️:                    DIST15_TOL;
    // RDKit❗✔️:               du =
    // RDKit❗✔️:                   _compute15DistsCisTrans(d4, d3, d2, d1, ang34, ang23, ang12) +
    // RDKit❗✔️:                   DIST15_TOL;
    // RDKit❗✔️:             } else if (accumData.transPaths.find(pathId) !=
    // RDKit❗✔️:                        accumData.transPaths.end()) {
    // RDKit❗✔️:               dl =
    // RDKit❗✔️:                   _compute15DistsTransCis(d4, d3, d2, d1, ang34, ang23, ang12) -
    // RDKit❗✔️:                   DIST15_TOL;
    // RDKit❗✔️:               du = _compute15DistsTransTrans(d4, d3, d2, d1, ang34, ang23,
    // RDKit❗✔️:                                              ang12) +
    // RDKit❗✔️:                    DIST15_TOL;
    // RDKit❗✔️:             } else {
    // RDKit❗✔️:               vw1 = PeriodicTable::getTable()->getRvdw(
    // RDKit❗✔️:                   mol.getAtomWithIdx(aid1)->getAtomicNum());
    // RDKit❗✔️:               vw5 = PeriodicTable::getTable()->getRvdw(
    // RDKit❗✔️:                   mol.getAtomWithIdx(aid5)->getAtomicNum());
    // RDKit❗✔️:               dl = VDW_SCALE_15 * (vw1 + vw5);
    // RDKit❗✔️:             }
    // RDKit❗✔️:           }
    // RDKit❗✔️:           if (du < 0.0) {
    // RDKit❗✔️:             du = MAX_UPPER;
    // RDKit❗✔️:           }
    // RDKit❗✔️:
    // RDKit❗✔️:           // std::cerr<<"3: "<<aid1<<"-"<<aid5<<std::endl;
    // RDKit❗✔️:           _checkAndSetBounds(aid1, aid5, dl, du, mmat);
    // RDKit❗✔️:           accumData.set15Atoms[aid1 * na + aid5] = 1;
    // RDKit❗✔️:           accumData.set15Atoms[aid5 * na + aid1] = 1;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }

    if na != mol.atoms.len()
        || nb != mol.bonds.len()
        || mmat.dimension() != na
        || dmat.len()
            != na
                .checked_mul(na)
                .ok_or(GraphBoundsError::Input("atom pair count overflow"))?
    {
        return Err(GraphBoundsError::Input(
            "Wrong size metric or distance matrix",
        ));
    }

    let aid2 = bond_pair_shared_atom(mol, accum_data, bid1, bid2)?;
    let aid1 = if mol.bonds[bid1].begin().index() == aid2 {
        mol.bonds[bid1].end().index()
    } else {
        mol.bonds[bid1].begin().index()
    };
    let aid3 = bond_pair_shared_atom(mol, accum_data, bid2, bid3)?;
    let aid4 = if mol.bonds[bid3].begin().index() == aid3 {
        mol.bonds[bid3].end().index()
    } else {
        mol.bonds[bid3].begin().index()
    };

    let d1 = accum_data.bond_lengths[bid1];
    let d2 = accum_data.bond_lengths[bid2];
    let d3 = accum_data.bond_lengths[bid3];
    let ang12 = accum_data.get_bond_angle(nb, bid1, bid2);
    let ang23 = accum_data.get_bond_angle(nb, bid2, bid3);

    for i in 0..nb {
        if accum_data.get_bond_adj(nb, bid3, i) != aid4 as i32 {
            continue;
        }
        let aid5 = if mol.bonds[i].begin().index() == aid4 {
            mol.bonds[i].end().index()
        } else {
            mol.bonds[i].begin().index()
        };

        let pid = aid1.min(aid5) * na + aid1.max(aid5);

        if accum_data.visited_bound(pid, DistType::Dist14) {
            return Ok(());
        }

        if dmat[aid1.max(aid5) * na + aid1.min(aid5)] < 3.9 {
            continue;
        }

        if aid1 == aid5 {
            continue;
        }

        let pid1 = aid1 * na + aid5;
        let pid2 = aid5 * na + aid1;
        if !(mmat.get_lower(aid1, aid5)? < DIST12_DELTA
            || accum_data.set15_atoms[pid1]
            || accum_data.set15_atoms[pid2])
        {
            continue;
        }

        let d4 = accum_data.bond_lengths[i];
        let ang34 = accum_data.get_bond_angle(nb, bid3, i);

        let path_id = bid2 as u64 * nb as u64 * nb as u64 + bid3 as u64 * nb as u64 + i as u64;

        let (dl, mut du) = match kind {
            Path14Kind::Cis => {
                if has_path_flag(&accum_data.cis_paths, path_id) {
                    let base = compute_15_dist_cis_cis(d1, d2, d3, d4, ang12, ang23, ang34);
                    (base - DIST15_TOL, base + DIST15_TOL)
                } else if has_path_flag(&accum_data.trans_paths, path_id) {
                    let base = compute_15_dist_cis_trans(d1, d2, d3, d4, ang12, ang23, ang34);
                    (base - DIST15_TOL, base + DIST15_TOL)
                } else {
                    (
                        compute_15_dist_cis_cis(d1, d2, d3, d4, ang12, ang23, ang34) - DIST15_TOL,
                        compute_15_dist_cis_trans(d1, d2, d3, d4, ang12, ang23, ang34) + DIST15_TOL,
                    )
                }
            }
            Path14Kind::Trans => {
                if has_path_flag(&accum_data.cis_paths, path_id) {
                    let base = compute_15_dist_trans_cis(d1, d2, d3, d4, ang12, ang23, ang34);
                    (base - DIST15_TOL, base + DIST15_TOL)
                } else if has_path_flag(&accum_data.trans_paths, path_id) {
                    let base = compute_15_dist_trans_trans(d1, d2, d3, d4, ang12, ang23, ang34);
                    (base - DIST15_TOL, base + DIST15_TOL)
                } else {
                    (
                        compute_15_dist_trans_cis(d1, d2, d3, d4, ang12, ang23, ang34) - DIST15_TOL,
                        compute_15_dist_trans_trans(d1, d2, d3, d4, ang12, ang23, ang34)
                            + DIST15_TOL,
                    )
                }
            }
            Path14Kind::Other => {
                if has_path_flag(&accum_data.cis_paths, path_id) {
                    (
                        compute_15_dist_cis_cis(d4, d3, d2, d1, ang34, ang23, ang12) - DIST15_TOL,
                        compute_15_dist_cis_trans(d4, d3, d2, d1, ang34, ang23, ang12) + DIST15_TOL,
                    )
                } else if has_path_flag(&accum_data.trans_paths, path_id) {
                    (
                        compute_15_dist_trans_cis(d4, d3, d2, d1, ang34, ang23, ang12) - DIST15_TOL,
                        compute_15_dist_trans_trans(d4, d3, d2, d1, ang34, ang23, ang12)
                            + DIST15_TOL,
                    )
                } else {
                    let vw1 = cosmolkit_core::van_der_waals_radius(mol.atoms[aid1].atomic_number())
                        .ok_or(GraphBoundsError::AtomicNumber(
                            mol.atoms[aid1].atomic_number(),
                        ))?;
                    let vw5 = cosmolkit_core::van_der_waals_radius(mol.atoms[aid5].atomic_number())
                        .ok_or(GraphBoundsError::AtomicNumber(
                            mol.atoms[aid5].atomic_number(),
                        ))?;
                    (VDW_SCALE_15 * (vw1 + vw5), MAX_UPPER)
                }
            }
        };

        if du < 0.0 {
            du = MAX_UPPER;
        }

        check_and_set_bounds(mmat, aid1, aid5, dl, du, false)?;

        accum_data.set15_atoms[pid1] = true;
        accum_data.set15_atoms[pid2] = true;
    }
    Ok(())
}
pub(super) fn collect_bonds_and_angles(
    mol: &TopologyBlock,
    bonds: &mut Vec<(i32, i32)>,
    angles: &mut Vec<Vec<i32>>,
) {
    // RDKit❗✔️: void collectBondsAndAngles(const ROMol &mol,
    // RDKit❗✔️:                            std::vector<std::pair<int, int>> &bonds,
    // RDKit❗✔️:                            std::vector<std::vector<int>> &angles) {
    // RDKit❗✔️:   bonds.resize(0);
    // RDKit❗✔️:   angles.resize(0);
    // RDKit❗✔️:   bonds.reserve(mol.getNumBonds());
    // RDKit❗✔️:   for (const auto bondi : mol.bonds()) {
    // RDKit❗✔️:     bonds.emplace_back(bondi->getBeginAtomIdx(), bondi->getEndAtomIdx());
    // RDKit❗✔️:
    // RDKit❗✔️:     for (unsigned int j = bondi->getIdx() + 1; j < mol.getNumBonds(); ++j) {
    // RDKit❗✔️:       const Bond *bondj = mol.getBondWithIdx(j);
    // RDKit❗✔️:       int aid11 = bondi->getBeginAtomIdx();
    // RDKit❗✔️:       int aid12 = bondi->getEndAtomIdx();
    // RDKit❗✔️:       int aid21 = bondj->getBeginAtomIdx();
    // RDKit❗✔️:       int aid22 = bondj->getEndAtomIdx();
    // RDKit❗✔️:       if (aid11 != aid21 && aid11 != aid22 && aid12 != aid21 &&
    // RDKit❗✔️:           aid12 != aid22) {
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       std::vector<int> tmp(4,
    // RDKit❗✔️:                            0);  // elements: aid1, aid2, flag for triple bonds
    // RDKit❗✔️:
    // RDKit❗✔️:       if (aid12 == aid21) {
    // RDKit❗✔️:         tmp[0] = aid11;
    // RDKit❗✔️:         tmp[1] = aid12;
    // RDKit❗✔️:         tmp[2] = aid22;
    // RDKit❗✔️:       } else if (aid12 == aid22) {
    // RDKit❗✔️:         tmp[0] = aid11;
    // RDKit❗✔️:         tmp[1] = aid12;
    // RDKit❗✔️:         tmp[2] = aid21;
    // RDKit❗✔️:       } else if (aid11 == aid21) {
    // RDKit❗✔️:         tmp[0] = aid12;
    // RDKit❗✔️:         tmp[1] = aid11;
    // RDKit❗✔️:         tmp[2] = aid22;
    // RDKit❗✔️:       } else if (aid11 == aid22) {
    // RDKit❗✔️:         tmp[0] = aid12;
    // RDKit❗✔️:         tmp[1] = aid11;
    // RDKit❗✔️:         tmp[2] = aid21;
    // RDKit❗✔️:       }
    // RDKit❗✔️:
    // RDKit❗✔️:       if (bondi->getBondType() == Bond::TRIPLE ||
    // RDKit❗✔️:           bondj->getBondType() == Bond::TRIPLE) {
    // RDKit❗✔️:         // triple bond
    // RDKit❗✔️:         tmp[3] = 1;
    // RDKit❗✔️:       } else if (bondi->getBondType() == Bond::DOUBLE &&
    // RDKit❗✔️:                  bondj->getBondType() == Bond::DOUBLE &&
    // RDKit❗✔️:                  mol.getAtomWithIdx(tmp[1])->getDegree() == 2) {
    // RDKit❗✔️:         // consecutive double bonds
    // RDKit❗✔️:         tmp[3] = 1;
    // RDKit❗✔️:       }
    // RDKit❗✔️:
    // RDKit❗✔️:       angles.push_back(tmp);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }

    bonds.clear();
    angles.clear();
    bonds.reserve(mol.bonds.len());
    for bondi in &mol.bonds {
        bonds.push((bondi.begin().index() as i32, bondi.end().index() as i32));

        for j in (bondi.id().index() + 1)..mol.bonds.len() {
            let bondj = &mol.bonds[j];
            let aid11 = bondi.begin().index() as i32;
            let aid12 = bondi.end().index() as i32;
            let aid21 = bondj.begin().index() as i32;
            let aid22 = bondj.end().index() as i32;
            if aid11 != aid21 && aid11 != aid22 && aid12 != aid21 && aid12 != aid22 {
                continue;
            }

            let mut tmp = vec![0; 4];
            if aid12 == aid21 {
                tmp[0] = aid11;
                tmp[1] = aid12;
                tmp[2] = aid22;
            } else if aid12 == aid22 {
                tmp[0] = aid11;
                tmp[1] = aid12;
                tmp[2] = aid21;
            } else if aid11 == aid21 {
                tmp[0] = aid12;
                tmp[1] = aid11;
                tmp[2] = aid22;
            } else if aid11 == aid22 {
                tmp[0] = aid12;
                tmp[1] = aid11;
                tmp[2] = aid21;
            }

            if bondi.order() == BondOrder::Triple || bondj.order() == BondOrder::Triple {
                tmp[3] = 1;
            } else if bondi.order() == BondOrder::Double
                && bondj.order() == BondOrder::Double
                && mol.adjacency.neighbors_of(tmp[1] as usize).len() == 2
            {
                tmp[3] = 1;
            }

            angles.push(tmp);
        }
    }
}
fn is_hbond_acceptor(atomic_num: u8) -> bool {
    // RDKit❗✔️: bool isHBondAcceptor(const Atom *atom) {
    // RDKit❗✔️:   return (atom->getAtomicNum() == 7 || atom->getAtomicNum() == 8);
    // RDKit❗✔️: }
    atomic_num == 7 || atomic_num == 8
}
fn is_h_in_hbond_donor(mol: &TopologyBlock, atom_idx: usize) -> bool {
    // RDKit❗✔️: bool isHinHBondDonor(const Atom *atom, const ROMol &mol) {
    // RDKit❗✔️:   if (atom->getAtomicNum() != 1) {
    // RDKit❗✔️:     return false;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   auto nbrs = mol.atomNeighbors(atom);
    // RDKit❗✔️:   return std::any_of(nbrs.begin(), nbrs.end(), [](const Atom *nbr) {
    // RDKit❗✔️:     return nbr->getAtomicNum() == 7 || nbr->getAtomicNum() == 8;
    // RDKit❗✔️:   });
    // RDKit❗✔️: }

    if mol.atoms[atom_idx].atomic_number() != 1 {
        return false;
    }
    mol.adjacency
        .neighbors_of(atom_idx)
        .iter()
        .any(|neighbor| is_hbond_acceptor(mol.atoms[neighbor.atom_index].atomic_number()))
}
fn set_lower_bound_vdw(
    mol: &TopologyBlock,
    mmat: &mut BoundsMatrix,
    _scale_vdw: bool,
    dmat: &[f64],
) -> Result<(), GraphBoundsError> {
    // RDKit❗❌: void setLowerBoundVDW(const ROMol &mol, DistGeom::BoundsMatPtr mmat, bool,
    // RDKit❗❌:                       double *dmat) {
    // RDKit❗❌:   unsigned int npt = mmat->numRows();
    // RDKit❗❌:   PRECONDITION(npt == mol.getNumAtoms(), "Wrong size metric matrix");
    // RDKit❗❌:
    // RDKit❗❌:   boost::dynamic_bitset<> hinHBondDonors(mol.getNumAtoms());
    // RDKit❗❌:   boost::dynamic_bitset<> hBondAcceptors(mol.getNumAtoms());
    // RDKit❗❌:   for (unsigned int i = 1; i < npt; i++) {
    // RDKit❗❌:     const auto atomI = mol.getAtomWithIdx(i);
    // RDKit❗❌:     auto vw1 = PeriodicTable::getTable()->getRvdw(atomI->getAtomicNum());
    // RDKit❗❌:     if (isHinHBondDonor(atomI, mol)) {
    // RDKit❗❌:       hinHBondDonors.set(i);
    // RDKit❗❌:     }
    // RDKit❗❌:     if (isHBondAcceptor(atomI)) {
    // RDKit❗❌:       hBondAcceptors.set(i);
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     for (unsigned int j = 0; j < i; j++) {
    // RDKit❗❌:       const auto atomJ = mol.getAtomWithIdx(j);
    // RDKit❗❌:       auto vw2 = PeriodicTable::getTable()->getRvdw(atomJ->getAtomicNum());
    // RDKit❗❌:       if (mmat->getLowerBound(i, j) < DIST12_DELTA) {
    // RDKit❗❌:         // ok this is what we are going to do
    // RDKit❗❌:         // - for atoms that are 4 or 5 bonds apart (15 or 16 distances), we
    // RDKit❗❌:         // will scale
    // RDKit❗❌:         //   the sum of the VDW radii so that the atoms can get closer
    // RDKit❗❌:         //   For 15 we will use VDW_SCALE_15 and for 16 we will use 1 -
    // RDKit❗❌:         //   0.5*VDW_SCALE_15
    // RDKit❗❌:         // - for all other pairs of atoms more than 5 bonds apart we use the
    // RDKit❗❌:         // sum of the VDW radii
    // RDKit❗❌:         //    as the lower bound
    // RDKit❗❌:         // - if one of the atoms is a H of a H-bond donor and the other is
    // RDKit❗❌:         //    an acceptor we will lower the bound to 1.8A
    // RDKit❗❌:         if ((hinHBondDonors[i] && hBondAcceptors[j]) ||
    // RDKit❗❌:             (hBondAcceptors[i] && hinHBondDonors[j])) {
    // RDKit❗❌:           mmat->setLowerBound(i, j, H_BOND_LENGTH);
    // RDKit❗❌:         } else if (dmat[i * npt + j] == 4.0) {
    // RDKit❗❌:           mmat->setLowerBound(i, j, VDW_SCALE_15 * (vw1 + vw2));
    // RDKit❗❌:         } else if (dmat[i * npt + j] == 5.0) {
    // RDKit❗❌:           mmat->setLowerBound(
    // RDKit❗❌:               i, j, (VDW_SCALE_15 + 0.5 * (1.0 - VDW_SCALE_15)) * (vw1 + vw2));
    // RDKit❗❌:         } else {
    // RDKit❗❌:           mmat->setLowerBound(i, j, (vw1 + vw2));
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌: }

    let n = mol.atoms.len();
    if mmat.dimension() != n
        || dmat.len()
            != n.checked_mul(n)
                .ok_or(GraphBoundsError::Input("atom pair count overflow"))?
    {
        return Err(GraphBoundsError::Input(
            "Wrong size metric or distance matrix",
        ));
    }
    let mut h_in_hbond_donors = vec![false; n];
    let mut hbond_acceptors = vec![false; n];
    for i in 1..n {
        let z1 = mol.atoms[i].atomic_number();
        let vw1 =
            cosmolkit_core::van_der_waals_radius(z1).ok_or(GraphBoundsError::AtomicNumber(z1))?;
        h_in_hbond_donors[i] = is_h_in_hbond_donor(mol, i);
        hbond_acceptors[i] = is_hbond_acceptor(z1);
        for j in 0..i {
            let z2 = mol.atoms[j].atomic_number();
            let vw2 = cosmolkit_core::van_der_waals_radius(z2)
                .ok_or(GraphBoundsError::AtomicNumber(z2))?;
            if mmat.get_lower(i, j)? < DIST12_DELTA {
                let bound = if (h_in_hbond_donors[i] && hbond_acceptors[j])
                    || (hbond_acceptors[i] && h_in_hbond_donors[j])
                {
                    H_BOND_LENGTH
                } else if dmat[i * n + j] == 4.0 {
                    VDW_SCALE_15 * (vw1 + vw2)
                } else if dmat[i * n + j] == 5.0 {
                    (VDW_SCALE_15 + 0.5 * (1.0 - VDW_SCALE_15)) * (vw1 + vw2)
                } else {
                    vw1 + vw2
                };
                mmat.set_lower(i, j, bound)?;
            }
        }
    }
    Ok(())
}
fn set_topol_bounds_stages(
    mol: &TopologyBlock,
    rings: &RingInfo,
    valence: &ValenceAssignment,
    hybridizations: &[Hybridization],
    conjugated: &[bool],
    mmat: &mut BoundsMatrix,
    set15bounds: bool,
    scale_vdw: bool,
    use_macrocycle_14config: bool,
    force_trans_amides: bool,
    set14bounds: bool,
    set13bounds: bool,
) -> Result<(), GraphBoundsError> {
    // RDKit❗❌: void setTopolBounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat,
    // RDKit❗❌:                     bool set15bounds, bool scaleVDW, bool useMacrocycle14config,
    // RDKit❗❌:                     bool forceTransAmides, bool set14bounds, bool set13bounds) {
    // RDKit❗❌:   PRECONDITION(mmat.get(), "bad pointer");
    // RDKit❗❌:   unsigned int nb = mol.getNumBonds();
    // RDKit❗❌:   unsigned int na = mol.getNumAtoms();
    // RDKit❗❌:   if (!na) {
    // RDKit❗❌:     throw ValueErrorException("molecule has no atoms");
    // RDKit❗❌:   }
    // RDKit❗❌:   // this is 2.6 million bonds, so it's extremly unlikely to ever occur, but
    // RDKit❗❌:   // we might as well check:
    // RDKit❗❌:   const size_t MAX_NUM_BONDS = static_cast<size_t>(
    // RDKit❗❌:       std::pow(std::numeric_limits<std::uint64_t>::max(), 1. / 3));
    // RDKit❗❌:   if (mol.getNumBonds() >= MAX_NUM_BONDS) {
    // RDKit❗❌:     throw ValueErrorException(
    // RDKit❗❌:         "Too many bonds in the molecule, cannot compute 1-4 bounds");
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   ComputedData accumData(na, nb);
    // RDKit❗❌:   double *distMatrix = nullptr;
    // RDKit❗❌:   distMatrix = MolOps::getDistanceMat(mol);
    // RDKit❗❌:
    // RDKit❗❌:   set12Bounds(mol, mmat, accumData);
    // RDKit❗❌:   if (set13bounds) {
    // RDKit❗❌:     set13Bounds(mol, mmat, accumData);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (set14bounds) {
    // RDKit❗❌:     set14Bounds(mol, mmat, accumData, distMatrix, useMacrocycle14config,
    // RDKit❗❌:                 forceTransAmides);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (set15bounds) {
    // RDKit❗❌:     set15Bounds(mol, mmat, accumData, distMatrix);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   setLowerBoundVDW(mol, mmat, scaleVDW, distMatrix);
    // RDKit❗❌: }

    let mut accum = ComputedData::new(mol.atoms.len(), mol.bonds.len())?;
    let distances = cosmolkit_core::topological_distance_matrix(
        mol,
        &cosmolkit_core::TopologicalDistanceMatrixParams::default(),
    )?;
    let dmat = distances.values();
    set_12_bounds(
        mol,
        rings,
        valence,
        hybridizations,
        conjugated,
        mmat,
        &mut accum,
    )?;
    if set13bounds {
        set_13_bounds(mol, mmat, &mut accum, rings)?;
    }
    if set14bounds {
        set_14_bounds(
            mol,
            valence,
            mmat,
            &mut accum,
            dmat,
            use_macrocycle_14config,
            force_trans_amides,
            rings,
        )?;
    }
    if set15bounds {
        set_15_bounds(mol, mmat, &mut accum, dmat)?;
    }
    set_lower_bound_vdw(mol, mmat, scale_vdw, dmat)
}
pub(super) fn set_topol_bounds(
    mol: &TopologyBlock,
    rings: &RingInfo,
    valence: &ValenceAssignment,
    hybridizations: &[Hybridization],
    conjugated: &[bool],
    mmat: &mut BoundsMatrix,
    set15bounds: bool,
    scale_vdw: bool,
    use_macrocycle_14config: bool,
    force_trans_amides: bool,
    set14bounds: bool,
    set13bounds: bool,
) -> Result<(), GraphBoundsError> {
    // RDKit❗✔️: void setTopolBounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat,
    // RDKit❗✔️:                     bool set15bounds, bool scaleVDW, bool useMacrocycle14config,
    // RDKit❗✔️:                     bool forceTransAmides, bool set14bounds, bool set13bounds) {
    // RDKit❗✔️:   PRECONDITION(mmat.get(), "bad pointer");
    // RDKit❗✔️:   unsigned int nb = mol.getNumBonds();
    // RDKit❗✔️:   unsigned int na = mol.getNumAtoms();
    // RDKit❗✔️:   if (!na) {
    // RDKit❗✔️:     throw ValueErrorException("molecule has no atoms");
    // RDKit❗✔️:   }
    // RDKit❗✔️:   // this is 2.6 million bonds, so it's extremly unlikely to ever occur, but
    // RDKit❗✔️:   // we might as well check:
    // RDKit❗✔️:   const size_t MAX_NUM_BONDS = static_cast<size_t>(
    // RDKit❗✔️:       std::pow(std::numeric_limits<std::uint64_t>::max(), 1. / 3));
    // RDKit❗✔️:   if (mol.getNumBonds() >= MAX_NUM_BONDS) {
    // RDKit❗✔️:     throw ValueErrorException(
    // RDKit❗✔️:         "Too many bonds in the molecule, cannot compute 1-4 bounds");
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   ComputedData accumData(na, nb);
    // RDKit❗✔️:   double *distMatrix = nullptr;
    // RDKit❗✔️:   distMatrix = MolOps::getDistanceMat(mol);
    // RDKit❗✔️:
    // RDKit❗✔️:   set12Bounds(mol, mmat, accumData);
    // RDKit❗✔️:   if (set13bounds) {
    // RDKit❗✔️:     set13Bounds(mol, mmat, accumData);
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   if (set14bounds) {
    // RDKit❗✔️:     set14Bounds(mol, mmat, accumData, distMatrix, useMacrocycle14config,
    // RDKit❗✔️:                 forceTransAmides);
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   if (set15bounds) {
    // RDKit❗✔️:     set15Bounds(mol, mmat, accumData, distMatrix);
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   setLowerBoundVDW(mol, mmat, scaleVDW, distMatrix);
    // RDKit❗✔️: }

    if mol.atoms.is_empty() {
        return Err(GraphBoundsError::Input("molecule has no atoms"));
    }
    let max_num_bonds = (u64::MAX as f64).powf(1.0 / 3.0) as usize;
    if mol.bonds.len() >= max_num_bonds {
        return Err(GraphBoundsError::Input(
            "Too many bonds in the molecule, cannot compute 1-4 bounds",
        ));
    }
    set_topol_bounds_stages(
        mol,
        rings,
        valence,
        hybridizations,
        conjugated,
        mmat,
        set15bounds,
        scale_vdw,
        use_macrocycle_14config,
        force_trans_amides,
        set14bounds,
        set13bounds,
    )
}
pub(super) fn set_topol_bounds_with_outputs(
    mol: &TopologyBlock,
    rings: &RingInfo,
    valence: &ValenceAssignment,
    hybridizations: &[Hybridization],
    conjugated: &[bool],
    mmat: &mut BoundsMatrix,
    bonds: &mut Vec<(i32, i32)>,
    angles: &mut Vec<Vec<i32>>,
    set15bounds: bool,
    scale_vdw: bool,
    use_macrocycle_14config: bool,
    force_trans_amides: bool,
    set14bounds: bool,
    set13bounds: bool,
) -> Result<(), GraphBoundsError> {
    // RDKit❗✔️: void setTopolBounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat,
    // RDKit❗✔️:                     std::vector<std::pair<int, int>> &bonds,
    // RDKit❗✔️:                     std::vector<std::vector<int>> &angles, bool set15bounds,
    // RDKit❗✔️:                     bool scaleVDW, bool useMacrocycle14config,
    // RDKit❗✔️:                     bool forceTransAmides, bool set14bounds, bool set13bounds) {
    // RDKit❗✔️:   PRECONDITION(mmat.get(), "bad pointer");
    // RDKit❗✔️:   bonds.clear();
    // RDKit❗✔️:   angles.clear();
    // RDKit❗✔️:   unsigned int nb = mol.getNumBonds();
    // RDKit❗✔️:   unsigned int na = mol.getNumAtoms();
    // RDKit❗✔️:   if (!na) {
    // RDKit❗✔️:     throw ValueErrorException("molecule has no atoms");
    // RDKit❗✔️:   }
    // RDKit❗✔️:   ComputedData accumData(na, nb);
    // RDKit❗✔️:   double *distMatrix = nullptr;
    // RDKit❗✔️:   distMatrix = MolOps::getDistanceMat(mol);
    // RDKit❗✔️:
    // RDKit❗✔️:   set12Bounds(mol, mmat, accumData);
    // RDKit❗✔️:
    // RDKit❗✔️:   if (set13bounds) {
    // RDKit❗✔️:     set13Bounds(mol, mmat, accumData);
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   if (set14bounds) {
    // RDKit❗✔️:     set14Bounds(mol, mmat, accumData, distMatrix, useMacrocycle14config,
    // RDKit❗✔️:                 forceTransAmides);
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   if (set15bounds) {
    // RDKit❗✔️:     set15Bounds(mol, mmat, accumData, distMatrix);
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   setLowerBoundVDW(mol, mmat, scaleVDW, distMatrix);
    // RDKit❗✔️:
    // RDKit❗✔️:   collectBondsAndAngles(mol, bonds, angles);
    // RDKit❗✔️: }

    bonds.clear();
    angles.clear();
    if mol.atoms.is_empty() {
        return Err(GraphBoundsError::Input("molecule has no atoms"));
    }
    set_topol_bounds_stages(
        mol,
        rings,
        valence,
        hybridizations,
        conjugated,
        mmat,
        set15bounds,
        scale_vdw,
        use_macrocycle_14config,
        force_trans_amides,
        set14bounds,
        set13bounds,
    )?;
    collect_bonds_and_angles(mol, bonds, angles);
    Ok(())
}

fn invalid_bounds(
    reason: &'static str,
    i: usize,
    j: usize,
    lb: f64,
    ub: f64,
    current_lower: f64,
    current_upper: f64,
) -> GraphBoundsError {
    GraphBoundsError::InvalidBounds(format!(
        "{reason} for atom pair ({i}, {j}); requested lower={lb:.6}, upper={ub:.6}, current lower={current_lower:.6}, current upper={current_upper:.6}"
    ))
}
pub(super) fn build_bounds_matrix(
    molecule: &TopologyBlock,
    rings: &RingInfo,
    valence: &ValenceAssignment,
    hybridizations: &[Hybridization],
    conjugated: &[bool],
    set15bounds: bool,
    scale_vdw: bool,
    do_triangle_smoothing: bool,
    use_macrocycle14config: bool,
) -> Result<BoundsMatrix, GraphBoundsError> {
    // RDKit❗✔️: PyObject *getMolBoundsMatrix(ROMol &mol, bool set15bounds = true,
    // RDKit❗✔️:                              bool scaleVDW = false,
    // RDKit❗✔️:                              bool doTriangleSmoothing = true,
    // RDKit❗✔️:                              bool useMacrocycle14config = false) {
    // RDKit❗✔️:   unsigned int nats = mol.getNumAtoms();
    // RDKit❗✔️:   npy_intp dims[2];
    // RDKit❗✔️:   dims[0] = nats;
    // RDKit❗✔️:   dims[1] = nats;
    // RDKit❗✔️:
    // RDKit❗✔️:   DistGeom::BoundsMatPtr mat(new DistGeom::BoundsMatrix(nats));
    // RDKit❗✔️:   DGeomHelpers::initBoundsMat(mat);
    // RDKit❗✔️:   DGeomHelpers::setTopolBounds(mol, mat, set15bounds, scaleVDW,
    // RDKit❗✔️:                                useMacrocycle14config);
    // RDKit❗✔️:   if (doTriangleSmoothing) {
    // RDKit❗✔️:     DistGeom::triangleSmoothBounds(mat);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   auto *res = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    // RDKit❗✔️:   memcpy(static_cast<void *>(PyArray_DATA(res)),
    // RDKit❗✔️:          static_cast<void *>(mat->getData()), nats * nats * sizeof(double));
    // RDKit❗✔️:
    // RDKit❗✔️:   return PyArray_Return(res);
    // RDKit❗✔️: }
    // This detached owner returns the unique dense buffer; canonical language
    // projections own their result materialization. No second matrix or chemistry.
    let mut bounds = BoundsMatrix::new(molecule.atoms.len())?;
    init_bounds_mat(&mut bounds, 0.0, 1000.0)?;
    set_topol_bounds(
        molecule,
        rings,
        valence,
        hybridizations,
        conjugated,
        &mut bounds,
        set15bounds,
        scale_vdw,
        use_macrocycle14config,
        true,
        true,
        true,
    )?;
    if do_triangle_smoothing {
        let _ = crate::smoothing::triangle_smooth_bounds_shared(&mut bounds, 0.0);
    }
    Ok(bounds)
}
#[cfg(test)]
mod tests {
    use super::*;
    use cosmolkit_core::{
        SanitizeParams, ValenceModel, assign_conjugation_flags,
        assign_hybridization_with_conjugation, assign_valence_for_topology, sanitize_topology,
        symmetrize_sssr_with_options_from_parts,
    };
    use cosmolkit_model::{Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, Element};
    fn fixed_topology(text: &str, sanitize: bool) -> TopologyBlock {
        let record = cosmolkit_smiles::parse_smiles(
            text,
            &cosmolkit_smiles::SmilesParseParams {
                sanitize,
                ..Default::default()
            },
        )
        .expect("original SMILES input");
        if sanitize {
            let topology = sanitize_topology(&record.topology, &SanitizeParams::default())
                .expect("original sanitization preparation")
                .topology;
            // RDKit❗✔️:     // figure out stereochemistry:
            // RDKit❗✔️:     bool cleanIt = true, force = true, flagPossible = true;
            // RDKit❗✔️:     MolOps::assignStereochemistry(*res, cleanIt, force, flagPossible);
            // The fixed original alkene input includes directional stereo.
            // Reuse the existing core stereo owner for detached preparation.
            cosmolkit_core::assign_double_bond_stereo_from_directions(topology)
                .expect("original directional double-bond stereo preparation")
        } else {
            record.topology
        }
    }
    fn run_set12_bounds(topology: &TopologyBlock) -> (BoundsMatrix, ComputedData) {
        let valence = assign_valence_for_topology(topology, ValenceModel::RdkitLike)
            .expect("original valence");
        let conjugated = assign_conjugation_flags(topology, &valence).expect("unique conjugation");
        let hybridization = assign_hybridization_with_conjugation(topology, &valence, &conjugated)
            .expect("unique hybridization");
        let rings = symmetrize_sssr_with_options_from_parts(
            topology.atoms.len(),
            &topology.bonds,
            &topology.adjacency,
            false,
            false,
        )
        .expect("original symmetrized rings");
        let mut bounds = initialized_bounds(topology.atoms.len());
        let mut accum = ComputedData::new(topology.atoms.len(), topology.bonds.len()).unwrap();
        set_12_bounds(
            topology,
            &rings,
            &valence,
            &hybridization.values,
            &conjugated,
            &mut bounds,
            &mut accum,
        )
        .expect("set12Bounds");
        (bounds, accum)
    }
    fn atom_total_valence_for_uff(assignment: &ValenceAssignment, index: usize) -> i32 {
        assignment.explicit_valence[index] + assignment.implicit_hydrogens[index]
    }
    fn vdw_radius(z: u8) -> f64 {
        cosmolkit_core::van_der_waals_radius(z).expect("original supported atomic number")
    }
    #[test]
    fn set12_bounds_uff_total_valence_includes_bracket_explicit_hydrogen() {
        let mol =
            Ok::<_, ()>(fixed_topology("C[SH+][O-]", true)).expect("explicit sulfur hydrogen");
        let assignment =
            assign_valence_for_topology(&mol, ValenceModel::RdkitLike).expect("valence");

        assert_eq!(assignment.explicit_valence[1], 3);
        assert_eq!(assignment.implicit_hydrogens[1], 0);
        assert_eq!(atom_total_valence_for_uff(&assignment, 1), 3);

        let (bounds, accum_data) = run_set12_bounds(&mol);
        assert!((accum_data.bond_lengths[0] - 1.776_838_813_449_356_2).abs() < 1.0e-14);
        assert!((accum_data.bond_lengths[1] - 1.679_472_654_030_108).abs() < 1.0e-14);
        assert!((bounds.get_upper(0, 1).unwrap() - 1.786_838_813_449_356_2).abs() < 1.0e-14);
        assert!((bounds.get_upper(1, 2).unwrap() - 1.689_472_654_030_108).abs() < 1.0e-14);
    }

    #[test]
    fn set12_bounds_uses_uff_rest_length_for_supported_atoms() {
        let mol = Ok::<_, ()>(fixed_topology("CC", false)).expect("ethane skeleton");

        let (mmat, accum_data) = run_set12_bounds(&mol);

        let lower = mmat.get_lower(0, 1).unwrap();
        let upper = mmat.get_upper(0, 1).unwrap();
        let width = upper - lower;
        assert!((width - (2.0 * DIST12_DELTA)).abs() < 1e-9);
        assert!(accum_data.bond_lengths[0] > 1.4);
        assert!(accum_data.bond_lengths[0] < 1.6);
        assert!((lower - (accum_data.bond_lengths[0] - DIST12_DELTA)).abs() < 1e-9);
        assert!((upper - (accum_data.bond_lengths[0] + DIST12_DELTA)).abs() < 1e-9);
        assert!(accum_data.visited12_bounds[1]);
    }

    #[test]
    fn set12_bounds_falls_back_to_vdw_when_uff_params_are_missing() {
        let atoms = [Element::C, Element::DUMMY]
            .into_iter()
            .enumerate()
            .map(|(i, e)| Atom::from_spec(AtomId::new(i), AtomSpec::new(e)))
            .collect();
        let bonds = vec![Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        )];
        let mol = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![])
            .expect("dummy-carbon molecule");

        let (mmat, accum_data) = run_set12_bounds(&mol);
        let expected = (vdw_radius(6) + vdw_radius(0)) / 2.0;

        assert!((accum_data.bond_lengths[0] - expected).abs() < 1e-9);
        assert!((mmat.get_lower(0, 1).unwrap() - (0.5 * expected)).abs() < 1e-9);
        assert!((mmat.get_upper(0, 1).unwrap() - (1.5 * expected)).abs() < 1e-9);
    }

    #[test]
    fn set12_bounds_marks_visited_pid_using_sorted_atom_indices() {
        let mol = Ok::<_, ()>(fixed_topology("CCO", false)).expect("ethanol skeleton");

        let (_mmat, accum_data) = run_set12_bounds(&mol);
        let n = mol.atoms.len();

        assert!(accum_data.visited12_bounds[1]);
        assert!(accum_data.visited12_bounds[n + 2]);
        assert!(!accum_data.visited12_bounds[2]);
    }

    #[test]
    fn set12_bounds_adds_extra_squish_for_conjugated_hetero_five_ring_bonds() {
        let mol = Ok::<_, ()>(fixed_topology("s1cccc1", true)).expect("thiophene");

        let (mmat, _accum_data) = run_set12_bounds(&mol);

        let sulfur_atom = mol
            .atoms
            .iter()
            .position(|atom| atom.atomic_number() == 16)
            .expect("sulfur atom");
        let squished_width = mol
            .bonds
            .iter()
            .find(|bond| bond.begin().index() == sulfur_atom || bond.end().index() == sulfur_atom)
            .map(|bond| {
                mmat.get_upper(bond.begin().index(), bond.end().index())
                    .unwrap()
                    - mmat
                        .get_lower(bond.begin().index(), bond.end().index())
                        .unwrap()
            })
            .expect("sulfur bond");
        let mut sulfur_adjacent = vec![false; mol.atoms.len()];
        sulfur_adjacent[sulfur_atom] = true;
        for bond in &mol.bonds {
            if bond.begin().index() == sulfur_atom {
                sulfur_adjacent[bond.end().index()] = true;
            } else if bond.end().index() == sulfur_atom {
                sulfur_adjacent[bond.begin().index()] = true;
            }
        }
        let unsquished_width = mol
            .bonds
            .iter()
            .find(|bond| {
                !sulfur_adjacent[bond.begin().index()] && !sulfur_adjacent[bond.end().index()]
            })
            .map(|bond| {
                mmat.get_upper(bond.begin().index(), bond.end().index())
                    .unwrap()
                    - mmat
                        .get_lower(bond.begin().index(), bond.end().index())
                        .unwrap()
            })
            .expect("carbon-carbon bond");

        assert!((squished_width - (2.0 * (0.2 + DIST12_DELTA))).abs() < 1e-9);
        assert!((unsquished_width - (2.0 * DIST12_DELTA)).abs() < 1e-9);
    }

    fn initialized_bounds(n: usize) -> BoundsMatrix {
        let mut m = BoundsMatrix::new(n).unwrap();
        init_bounds_mat(&mut m, 0.001, 1000.0).unwrap();
        m
    }
    fn fixed_rings(mol: &TopologyBlock) -> RingInfo {
        symmetrize_sssr_with_options_from_parts(
            mol.atoms.len(),
            &mol.bonds,
            &mol.adjacency,
            false,
            false,
        )
        .expect("original rings")
    }
    fn neighbors_for_atom(mol: &TopologyBlock, i: usize) -> Vec<usize> {
        mol.adjacency
            .neighbors_of(i)
            .iter()
            .map(|n| n.atom_index)
            .collect()
    }
    fn run_set13_bounds(mol: &TopologyBlock) -> (BoundsMatrix, ComputedData) {
        let (mut mat, mut accum) = run_set12_bounds(mol);
        set_13_bounds(mol, &mut mat, &mut accum, &fixed_rings(mol)).expect("original set13Bounds");
        (mat, accum)
    }
    #[test]
    fn set_ring_angle_matches_rdkit_ring_hybridization_special_cases() {
        let make = |hyb| {
            TopologyBlock::try_from_parts(
                vec![Atom::from_spec(
                    AtomId::new(0),
                    AtomSpec::new(Element::C).with_hybridization(hyb),
                )],
                vec![],
                vec![],
                vec![],
            )
            .unwrap()
        };
        let cyclopropane_like = make(Hybridization::Sp2);
        let cyclopentane_like = make(Hybridization::Sp3);

        let tri = set_ring_angle(&cyclopropane_like, 0, 3);
        let five = set_ring_angle(&cyclopentane_like, 0, 5);

        assert!((tri - std::f64::consts::PI / 3.0).abs() < 1e-9);
        assert!((five - (104.0_f64.to_radians())).abs() < 1e-9);
    }

    #[test]
    fn compute13_dist_returns_bond_sum_for_linear_angle() {
        let dist = compute_13_dist(1.4, 1.5, std::f64::consts::PI);

        assert!((dist - 2.9).abs() < 1e-9);
    }

    #[test]
    fn compute13_dist_returns_bond_difference_for_zero_angle() {
        let dist = compute_13_dist(1.5, 1.4, 0.0);

        assert!((dist - 0.1).abs() < 1e-9);
    }

    #[test]
    fn set_13_bounds_helper_doubles_tolerance_for_larger_sp2_ring_atoms() {
        let thiophene = Ok::<_, ()>(fixed_topology("s1cccc1", true)).expect("thiophene");
        let benzene = Ok::<_, ()>(fixed_topology("c1ccccc1", true)).expect("benzene");

        let (_thiophene_mmat_12, thiophene_accum) = run_set12_bounds(&thiophene);
        let (_benzene_mmat_12, benzene_accum) = run_set12_bounds(&benzene);
        let thiophene_rings = &fixed_rings(&thiophene);
        let benzene_rings = &fixed_rings(&benzene);

        let sulfur = thiophene
            .atoms
            .iter()
            .position(|atom| atom.atomic_number() == 16)
            .expect("sulfur");
        let sulfur_neighbors: Vec<usize> = thiophene
            .bonds
            .iter()
            .filter_map(|bond| {
                if bond.begin().index() == sulfur {
                    Some(bond.end().index())
                } else if bond.end().index() == sulfur {
                    Some(bond.begin().index())
                } else {
                    None
                }
            })
            .collect();
        let mut thiophene_mmat = initialized_bounds(thiophene.atoms.len());
        let thiophene_angle = set_ring_angle(&thiophene, sulfur, 5);
        set_13_bounds_helper(
            sulfur_neighbors[0],
            sulfur,
            sulfur_neighbors[1],
            thiophene_angle,
            &thiophene_accum.bond_lengths,
            &mut thiophene_mmat,
            &thiophene,
            thiophene_rings,
        )
        .expect("thiophene 1-3 bounds");

        let mut benzene_mmat = initialized_bounds(benzene.atoms.len());
        let benzene_angle = set_ring_angle(&benzene, 1, 6);
        set_13_bounds_helper(
            0,
            1,
            2,
            benzene_angle,
            &benzene_accum.bond_lengths,
            &mut benzene_mmat,
            &benzene,
            benzene_rings,
        )
        .expect("benzene 1-3 bounds");

        let thiophene_width = thiophene_mmat
            .get_upper(sulfur_neighbors[0], sulfur_neighbors[1])
            .unwrap()
            - thiophene_mmat
                .get_lower(sulfur_neighbors[0], sulfur_neighbors[1])
                .unwrap();
        let benzene_width =
            benzene_mmat.get_upper(0, 2).unwrap() - benzene_mmat.get_lower(0, 2).unwrap();

        assert!((thiophene_width - (4.0 * DIST13_TOL)).abs() < 1e-9);
        assert!((benzene_width - (2.0 * DIST13_TOL)).abs() < 1e-9);
        assert!(is_larger_sp2_atom_idx(&thiophene, thiophene_rings, sulfur));
        assert!(!is_larger_sp2_atom_idx(&benzene, benzene_rings, 1));
    }

    #[test]
    fn visited_bound_obeys_rdkit_dist_type_thresholds() {
        let mut accum = ComputedData::new(4, 2).unwrap();
        accum.visited13_bounds[3] = true;
        accum.visited14_bounds[5] = true;

        assert!(!accum.visited_bound(3, DistType::Dist12));
        assert!(accum.visited_bound(3, DistType::Dist13));
        assert!(accum.visited_bound(3, DistType::Dist14));
        assert!(!accum.visited_bound(5, DistType::Dist13));
        assert!(accum.visited_bound(5, DistType::Dist14));
    }

    #[test]
    fn set_13_bounds_sets_non_ring_sp3_bounds_for_propane_path() {
        let mol = Ok::<_, ()>(fixed_topology("CCC", true)).expect("propane");
        let (mmat, accum_data) = run_set13_bounds(&mol);

        let bid01 = bond_between_idx_simple(&mol, 0, 1).expect("0-1 bond");
        let bid12 = bond_between_idx_simple(&mol, 1, 2).expect("1-2 bond");
        let expected = compute_13_dist(
            accum_data.bond_lengths[bid01],
            accum_data.bond_lengths[bid12],
            109.5_f64.to_radians(),
        );

        assert!((mmat.get_lower(0, 2).unwrap() - (expected - DIST13_TOL)).abs() < 1e-9);
        assert!((mmat.get_upper(0, 2).unwrap() - (expected + DIST13_TOL)).abs() < 1e-9);
        assert!(
            (accum_data.get_bond_angle(mol.bonds.len(), bid01, bid12) - 109.5_f64.to_radians())
                .abs()
                < 1e-9
        );
        assert!(accum_data.visited13_bounds[2]);
    }

    #[test]
    fn set_13_bounds_distributes_remaining_fused_ring_angle_like_rdkit() {
        let mol = Ok::<_, ()>(fixed_topology("c1cccc2ccccc12", true)).expect("naphthalene");
        let (_mmat, accum_data) = run_set13_bounds(&mol);
        let rings = &fixed_rings(&mol);

        let fusion_atom = mol
            .atoms
            .iter()
            .enumerate()
            .find_map(|(idx, atom)| {
                (atom.hybridization() == Hybridization::Sp2
                    && rings.num_atom_rings(atom.id()) > 1
                    && neighbors_for_atom(&mol, idx).len() == 3)
                    .then_some(idx)
            })
            .expect("fusion atom");

        let neighbors = neighbors_for_atom(&mol, fusion_atom);
        let mut pair_angles = Vec::new();
        for left in 0..neighbors.len() {
            let bid1 = bond_between_idx_simple(&mol, fusion_atom, neighbors[left]).expect("bond 1");
            for right in 0..left {
                let bid2 =
                    bond_between_idx_simple(&mol, fusion_atom, neighbors[right]).expect("bond 2");
                pair_angles.push(accum_data.get_bond_angle(mol.bonds.len(), bid1, bid2));
            }
        }

        assert_eq!(pair_angles.len(), 3);
        for angle in pair_angles {
            assert!((angle - 120.0_f64.to_radians()).abs() < 1e-9);
        }
    }

    #[test]
    fn set_13_bounds_uses_wide_bounds_for_non_ring_degree_five_center() {
        let mol = Ok::<_, ()>(fixed_topology("FP(F)(F)(F)F", true)).expect("PF5-like");
        let center = mol
            .atoms
            .iter()
            .position(|atom| atom.atomic_number() == 15)
            .expect("phosphorus center");
        let ligands: Vec<usize> = neighbors_for_atom(&mol, center);
        assert_eq!(ligands.len(), 5);

        let (mmat, accum_data) = run_set13_bounds(&mol);
        let bid1 = bond_between_idx_simple(&mol, center, ligands[0]).expect("P-F1");
        let bid2 = bond_between_idx_simple(&mol, center, ligands[1]).expect("P-F2");
        let dmax = accum_data.bond_lengths[bid1] + accum_data.bond_lengths[bid2];

        assert!((mmat.get_lower(ligands[0], ligands[1]).unwrap() - 1.0).abs() < 1e-9);
        assert!((mmat.get_upper(ligands[0], ligands[1]).unwrap() - (dmax * 1.2)).abs() < 1e-9);
        assert!(accum_data.visited13_bounds[ligands[0] * mol.atoms.len() + ligands[1]]);
    }

    fn fixed_valence(mol: &TopologyBlock) -> ValenceAssignment {
        assign_valence_for_topology(mol, ValenceModel::RdkitLike).expect("original valence")
    }
    fn flatten_topological_distances(mol: &TopologyBlock) -> Vec<f64> {
        cosmolkit_core::topological_distance_matrix(
            mol,
            &cosmolkit_core::TopologicalDistanceMatrixParams::default(),
        )
        .expect("unique source topological distances")
        .values()
        .to_vec()
    }
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    enum Set14DispatchCase {
        TwoSameRing,
        TwoDiffRing,
        ShareRingBond,
        Chain,
    }
    fn run_set14_bounds(
        mol: &TopologyBlock,
        use_macrocycle_14config: bool,
        force_trans_amides: bool,
    ) -> (BoundsMatrix, ComputedData, Vec<f64>) {
        let (mut mmat, mut accum_data) = run_set13_bounds(mol);
        let dmat = flatten_topological_distances(mol);
        set_14_bounds(
            mol,
            &fixed_valence(mol),
            &mut mmat,
            &mut accum_data,
            &dmat,
            use_macrocycle_14config,
            force_trans_amides,
            &fixed_rings(mol),
        )
        .expect("set14Bounds");
        (mmat, accum_data, dmat)
    }
    fn run_set14_same_ring_pass_only(
        mol: &TopologyBlock,
        use_macrocycle_14config: bool,
    ) -> (BoundsMatrix, ComputedData, Vec<f64>) {
        let (mut mmat, mut accum_data) = run_set13_bounds(mol);
        let dmat = flatten_topological_distances(mol);
        let rinfo = fixed_rings(mol);
        for bring in rinfo.bond_rings() {
            let r_size = bring.len();
            if r_size < 3 {
                continue;
            }
            let mut bid1 = bring[r_size - 1].index();
            for i in 0..r_size {
                let bid2 = bring[i].index();
                let bid3 = bring[(i + 1) % r_size].index();
                if r_size > 5 {
                    if use_macrocycle_14config && r_size >= MIN_MACROCYCLE_RING_SIZE {
                        set_macrocycle_all_in_same_ring_14_bounds(
                            mol,
                            &fixed_valence(mol),
                            bid1,
                            bid2,
                            bid3,
                            &mut accum_data,
                            &mut mmat,
                        )
                        .expect("macrocycle same-ring 1-4 bounds");
                    } else {
                        set_in_ring_14_bounds(
                            mol,
                            bid1,
                            bid2,
                            bid3,
                            &mut accum_data,
                            &mut mmat,
                            &dmat,
                            r_size,
                            &rinfo,
                        )
                        .expect("in-ring 1-4 bounds");
                    }
                } else {
                    record_14_path(mol, bid1, bid2, bid3, &mut accum_data).expect("record14Path");
                }
                bid1 = bid2;
            }
        }
        (mmat, accum_data, dmat)
    }
    fn find_same_ring_dispatch_triple(
        mol: &TopologyBlock,
        use_macrocycle_14config: bool,
    ) -> Option<(usize, usize, usize, usize)> {
        let rinfo = fixed_rings(mol);
        for bring in rinfo.bond_rings() {
            let r_size = bring.len();
            if r_size < 3 {
                continue;
            }
            if r_size > 5 && (!use_macrocycle_14config || r_size >= MIN_MACROCYCLE_RING_SIZE) {
                let bid1 = bring[r_size - 1].index();
                let bid2 = bring[0].index();
                let bid3 = bring[1].index();
                return Some((bid1, bid2, bid3, r_size));
            }
        }
        None
    }
    fn find_dispatch_triple(
        mol: &TopologyBlock,
        target: Set14DispatchCase,
        use_macrocycle_14config: bool,
    ) -> Option<(usize, usize, usize)> {
        let rinfo = fixed_rings(mol);
        let mut bid_is_macrocycle: HashSet<usize> = HashSet::new();
        let mut ring_bond_pairs: HashSet<u64> = HashSet::new();
        let mut done_paths: HashSet<u64> = HashSet::new();
        let nb = mol.bonds.len() as u64;

        for bring in rinfo.bond_rings() {
            let r_size = bring.len();
            if r_size < 3 {
                continue;
            }
            let mut bid1 = bring[r_size - 1].index();
            for i in 0..r_size {
                let bid2 = bring[i].index();
                let bid3 = bring[(i + 1) % r_size].index();
                ring_bond_pairs.insert(bid1 as u64 * nb + bid2 as u64);
                ring_bond_pairs.insert(bid2 as u64 * nb + bid1 as u64);
                done_paths.insert(bid1 as u64 * nb * nb + bid2 as u64 * nb + bid3 as u64);
                done_paths.insert(bid3 as u64 * nb * nb + bid2 as u64 * nb + bid1 as u64);
                if use_macrocycle_14config && r_size >= MIN_MACROCYCLE_RING_SIZE {
                    bid_is_macrocycle.insert(bid2);
                }
                bid1 = bid2;
            }
        }

        for bond in &mol.bonds {
            let bid2 = bond.id().index();
            let aid2 = bond.begin().index();
            let aid3 = bond.end().index();
            for nbr1 in neighbors_for_atom(mol, aid2) {
                let Some(bid1) = bond_between_idx_simple(mol, aid2, nbr1) else {
                    continue;
                };
                if bid1 == bid2 {
                    continue;
                }
                for nbr3 in neighbors_for_atom(mol, aid3) {
                    let Some(bid3) = bond_between_idx_simple(mol, aid3, nbr3) else {
                        continue;
                    };
                    if bid3 == bid2 {
                        continue;
                    }
                    let id1 = bid1 as u64 * nb * nb + bid2 as u64 * nb + bid3 as u64;
                    let id2 = bid3 as u64 * nb * nb + bid2 as u64 * nb + bid1 as u64;
                    if done_paths.contains(&id1) || done_paths.contains(&id2) {
                        continue;
                    }
                    let pid1 = bid1 as u64 * nb + bid2 as u64;
                    let pid2 = bid2 as u64 * nb + bid1 as u64;
                    let pid3 = bid2 as u64 * nb + bid3 as u64;
                    let pid4 = bid3 as u64 * nb + bid2 as u64;
                    let case = if ring_bond_pairs.contains(&pid1)
                        || ring_bond_pairs.contains(&pid2)
                        || ring_bond_pairs.contains(&pid3)
                        || ring_bond_pairs.contains(&pid4)
                    {
                        Set14DispatchCase::TwoSameRing
                    } else if (rinfo.num_bond_rings(BondId::new(bid1)) > 0
                        && rinfo.num_bond_rings(BondId::new(bid2)) > 0)
                        || (rinfo.num_bond_rings(BondId::new(bid2)) > 0
                            && rinfo.num_bond_rings(BondId::new(bid3)) > 0)
                    {
                        Set14DispatchCase::TwoDiffRing
                    } else if rinfo.num_bond_rings(BondId::new(bid2)) > 0 {
                        Set14DispatchCase::ShareRingBond
                    } else {
                        Set14DispatchCase::Chain
                    };
                    if case == target {
                        return Some((bid1, bid2, bid3));
                    }
                }
            }
        }
        let _ = bid_is_macrocycle;
        None
    }
    fn ideal_bond_angle(hybridization: &Hybridization, ring_size: Option<usize>) -> f64 {
        const DEG_TO_RAD: f64 = std::f64::consts::PI / 180.0;

        if let Some(rsize) = ring_size {
            match hybridization {
                Hybridization::Sp2 if rsize <= 8 || rsize == 3 || rsize == 4 => {
                    std::f64::consts::PI * (1.0 - 2.0 / rsize as f64)
                }
                Hybridization::Sp3 if rsize == 5 => 104.0 * DEG_TO_RAD,
                Hybridization::Sp3 => 109.5 * DEG_TO_RAD,
                Hybridization::Sp3d => 105.0 * DEG_TO_RAD,
                Hybridization::Sp3d2 => 90.0 * DEG_TO_RAD,
                _ => 120.0 * DEG_TO_RAD,
            }
        } else {
            match hybridization {
                Hybridization::Sp => 180.0 * DEG_TO_RAD,
                Hybridization::Sp2 => 120.0 * DEG_TO_RAD,
                Hybridization::Sp3 => 109.5 * DEG_TO_RAD,
                Hybridization::Sp3d => 90.0 * DEG_TO_RAD,
                Hybridization::Sp3d2 => 90.0 * DEG_TO_RAD,
                Hybridization::Sp2d => 120.0 * DEG_TO_RAD,
                _ => 120.0 * DEG_TO_RAD,
            }
        }
    }

    #[test]
    fn chain_and_carbonyl_classification_helpers_follow_rdkit_patterns() {
        let propane = Ok::<_, ()>(fixed_topology("CCC", true)).expect("propane");
        assert!(
            check_h2_nx3_h1_ox2(&propane, &fixed_valence(&propane), 1)
                .expect("original shared valence")
        );

        let ether = Ok::<_, ()>(fixed_topology("COC", true)).expect("dimethyl ether");
        let oxygen = ether
            .atoms
            .iter()
            .position(|atom| atom.atomic_number() == 8)
            .expect("ether oxygen");
        assert!(
            check_h2_nx3_h1_ox2(&ether, &fixed_valence(&ether), oxygen)
                .expect("original shared valence")
        );

        let butane = Ok::<_, ()>(fixed_topology("CCCC", true)).expect("butane");
        assert!(
            check_nh_ch_ch_nh(&butane, &fixed_valence(&butane), 0, 1, 2, 3)
                .expect("original shared valence")
        );

        let acetate = Ok::<_, ()>(fixed_topology("CC(=O)O", true)).expect("acetate");
        let carbonyl = acetate
            .atoms
            .iter()
            .enumerate()
            .find_map(|(idx, _)| is_carbonyl(&acetate, idx).then_some(idx))
            .expect("carbonyl carbon");
        assert!(is_carbonyl(&acetate, carbonyl));
        assert!(!is_carbonyl(&acetate, 0));
    }
    #[test]
    fn amide_ester_classification_helpers_match_ester_patterns() {
        let ester = Ok::<_, ()>(fixed_topology("COC(=O)C", true)).expect("methyl acetate");
        let carbonyl = ester
            .atoms
            .iter()
            .enumerate()
            .find_map(|(idx, _)| is_carbonyl(&ester, idx).then_some(idx))
            .expect("carbonyl carbon");

        let double_hetero = neighbors_for_atom(&ester, carbonyl)
            .into_iter()
            .find(|&nbr| {
                let bond_idx = bond_between_idx_simple(&ester, carbonyl, nbr).expect("bond");
                ester.bonds[bond_idx].order() == BondOrder::Double
                    && (ester.atoms[nbr].atomic_number() == 8
                        || ester.atoms[nbr].atomic_number() == 7)
            })
            .expect("double bonded hetero");
        let single_hetero = neighbors_for_atom(&ester, carbonyl)
            .into_iter()
            .find(|&nbr| {
                let bond_idx = bond_between_idx_simple(&ester, carbonyl, nbr).expect("bond");
                ester.bonds[bond_idx].order() == BondOrder::Single
                    && (ester.atoms[nbr].atomic_number() == 8
                        || ester.atoms[nbr].atomic_number() == 7)
            })
            .expect("single bonded hetero");
        let atom1 = neighbors_for_atom(&ester, single_hetero)
            .into_iter()
            .find(|&nbr| nbr != carbonyl)
            .expect("atom1");
        let bnd1 = bond_between_idx_simple(&ester, atom1, single_hetero).expect("bond1");
        let bnd3_double =
            bond_between_idx_simple(&ester, carbonyl, double_hetero).expect("bond3 double");
        let carbonyl_substituent = neighbors_for_atom(&ester, carbonyl)
            .into_iter()
            .find(|&nbr| nbr != single_hetero && nbr != double_hetero)
            .expect("carbonyl substituent");
        let bnd3_single =
            bond_between_idx_simple(&ester, carbonyl, carbonyl_substituent).expect("bond3 single");

        assert!(
            check_amide_ester_14(
                &ester,
                &fixed_valence(&ester),
                bnd1,
                bnd3_double,
                single_hetero,
                carbonyl,
                double_hetero,
            )
            .expect("original shared valence")
        );
        assert!(
            check_amide_ester_15(
                &ester,
                &fixed_valence(&ester),
                bnd1,
                bnd3_single,
                single_hetero,
                carbonyl,
            )
            .expect("original shared valence")
        );

        let tertiary_amide =
            Ok::<_, ()>(fixed_topology("CN(C)C(=O)C", true)).expect("tertiary amide");
        let tertiary_carbonyl = tertiary_amide
            .atoms
            .iter()
            .enumerate()
            .find_map(|(idx, _)| is_carbonyl(&tertiary_amide, idx).then_some(idx))
            .expect("tertiary amide carbonyl");
        let tertiary_nitrogen = neighbors_for_atom(&tertiary_amide, tertiary_carbonyl)
            .into_iter()
            .find(|&nbr| {
                let bond_idx =
                    bond_between_idx_simple(&tertiary_amide, tertiary_carbonyl, nbr).expect("bond");
                tertiary_amide.bonds[bond_idx].order() == BondOrder::Single
                    && tertiary_amide.atoms[nbr].atomic_number() == 7
            })
            .expect("amide nitrogen");
        let tertiary_atom1 = neighbors_for_atom(&tertiary_amide, tertiary_nitrogen)
            .into_iter()
            .find(|&nbr| nbr != tertiary_carbonyl)
            .expect("substituent carbon");
        let tertiary_bnd1 =
            bond_between_idx_simple(&tertiary_amide, tertiary_atom1, tertiary_nitrogen)
                .expect("tertiary bond1");
        let tertiary_side = neighbors_for_atom(&tertiary_amide, tertiary_carbonyl)
            .into_iter()
            .find(|&nbr| {
                let bond_idx =
                    bond_between_idx_simple(&tertiary_amide, tertiary_carbonyl, nbr).expect("bond");
                tertiary_amide.bonds[bond_idx].order() == BondOrder::Single
                    && nbr != tertiary_nitrogen
            })
            .expect("carbonyl side");
        let tertiary_bnd3 =
            bond_between_idx_simple(&tertiary_amide, tertiary_carbonyl, tertiary_side)
                .expect("tertiary bond3");

        assert!(
            !check_amide_ester_15(
                &tertiary_amide,
                &fixed_valence(&tertiary_amide),
                tertiary_bnd1,
                tertiary_bnd3,
                tertiary_nitrogen,
                tertiary_carbonyl,
            )
            .expect("original shared valence")
        );
    }
    #[test]
    fn macrocycle_all_in_same_ring_amide_ester_helper_matches_lactone_pattern() {
        let mol = Ok::<_, ()>(fixed_topology("O=C1N(C)CCCCCCCC1", true))
            .expect("macrocyclic tertiary lactam");
        let carbonyl = mol
            .atoms
            .iter()
            .enumerate()
            .find_map(|(idx, _)| is_carbonyl(&mol, idx).then_some(idx))
            .expect("carbonyl carbon");
        let atm2 = neighbors_for_atom(&mol, carbonyl)
            .into_iter()
            .find(|&nbr| {
                let bond_idx = bond_between_idx_simple(&mol, carbonyl, nbr).expect("bond");
                mol.bonds[bond_idx].order() == BondOrder::Single
                    && mol.atoms[nbr].atomic_number() == 7
            })
            .expect("ring amide nitrogen");
        let atm4 = neighbors_for_atom(&mol, carbonyl)
            .into_iter()
            .find(|&nbr| {
                let bond_idx = bond_between_idx_simple(&mol, carbonyl, nbr).expect("bond");
                mol.bonds[bond_idx].order() == BondOrder::Single && nbr != atm2
            })
            .expect("ring carbon neighbor");
        let atm1 = neighbors_for_atom(&mol, atm2)
            .into_iter()
            .find(|&nbr| nbr != carbonyl && mol.atoms[nbr].atomic_number() == 6)
            .expect("preceding ring atom");

        assert!(check_macrocycle_all_in_same_ring_amide_ester_14(
            &mol, atm1, atm2, carbonyl, atm4,
        ));
    }
    #[test]
    fn macrocycle_two_in_same_ring_amide_ester_helper_matches_tertiary_lactam_pattern() {
        let mol = Ok::<_, ()>(fixed_topology("O=C1N(C)CCCCCCCC1", true))
            .expect("macrocyclic tertiary lactam");
        let carbonyl = mol
            .atoms
            .iter()
            .enumerate()
            .find_map(|(idx, _)| is_carbonyl(&mol, idx).then_some(idx))
            .expect("carbonyl carbon");
        let oxygen = neighbors_for_atom(&mol, carbonyl)
            .into_iter()
            .find(|&nbr| {
                let bond_idx = bond_between_idx_simple(&mol, carbonyl, nbr).expect("bond");
                mol.bonds[bond_idx].order() == BondOrder::Double
            })
            .expect("carbonyl oxygen");
        let nitrogen = neighbors_for_atom(&mol, carbonyl)
            .into_iter()
            .find(|&nbr| {
                let bond_idx = bond_between_idx_simple(&mol, carbonyl, nbr).expect("bond");
                mol.bonds[bond_idx].order() == BondOrder::Single
                    && mol.atoms[nbr].atomic_number() == 7
            })
            .expect("amide nitrogen");
        let atom1 = neighbors_for_atom(&mol, nitrogen)
            .into_iter()
            .find(|&nbr| nbr != carbonyl && mol.atoms[nbr].atomic_number() == 6)
            .expect("preceding ring carbon");
        let bnd1 = bond_between_idx_simple(&mol, atom1, nitrogen).expect("bond1");
        let bnd3 = bond_between_idx_simple(&mol, carbonyl, oxygen).expect("bond3");

        assert!(check_macrocycle_two_in_same_ring_amide_ester_14(
            &mol, bnd1, bnd3, atom1, nitrogen, carbonyl, oxygen,
        ));
    }
    #[test]
    fn set_macrocycle_two_in_same_ring_14_bounds_uses_cis_for_macrocycle_amide_path() {
        let mol = Ok::<_, ()>(fixed_topology("O=C1N(C)CCCCCCCC1", true))
            .expect("macrocyclic tertiary lactam");
        let (mut mmat, mut accum_data) = run_set13_bounds(&mol);
        let dmat = flatten_topological_distances(&mol);

        let carbonyl = mol
            .atoms
            .iter()
            .enumerate()
            .find_map(|(idx, _)| is_carbonyl(&mol, idx).then_some(idx))
            .expect("carbonyl carbon");
        let oxygen = neighbors_for_atom(&mol, carbonyl)
            .into_iter()
            .find(|&nbr| {
                let bond_idx = bond_between_idx_simple(&mol, carbonyl, nbr).expect("bond");
                mol.bonds[bond_idx].order() == BondOrder::Double
            })
            .expect("carbonyl oxygen");
        let nitrogen = neighbors_for_atom(&mol, carbonyl)
            .into_iter()
            .find(|&nbr| {
                let bond_idx = bond_between_idx_simple(&mol, carbonyl, nbr).expect("bond");
                mol.bonds[bond_idx].order() == BondOrder::Single
                    && mol.atoms[nbr].atomic_number() == 7
            })
            .expect("amide nitrogen");
        let atom1 = neighbors_for_atom(&mol, nitrogen)
            .into_iter()
            .find(|&nbr| nbr != carbonyl && mol.atoms[nbr].atomic_number() == 6)
            .expect("preceding ring carbon");

        let bid1 = bond_between_idx_simple(&mol, atom1, nitrogen).expect("b1");
        let bid2 = bond_between_idx_simple(&mol, nitrogen, carbonyl).expect("b2");
        let bid3 = bond_between_idx_simple(&mol, carbonyl, oxygen).expect("b3");
        set_macrocycle_two_in_same_ring_14_bounds(
            &mol,
            bid1,
            bid2,
            bid3,
            &mut accum_data,
            &mut mmat,
            &dmat,
        )
        .expect("macrocycle two-in-same-ring bounds");

        let path = accum_data.paths14.last().expect("path");
        assert_eq!(path.kind, Path14Kind::Cis);
        let expected = compute_14_dist_cis(
            accum_data.bond_lengths[bid1],
            accum_data.bond_lengths[bid2],
            accum_data.bond_lengths[bid3],
            accum_data.get_bond_angle(mol.bonds.len(), bid1, bid2),
            accum_data.get_bond_angle(mol.bonds.len(), bid2, bid3),
        );
        assert!((mmat.get_lower(atom1, oxygen).unwrap() - (expected - GEN_DIST_TOL)).abs() < 1e-9);
        assert!((mmat.get_upper(atom1, oxygen).unwrap() - (expected + GEN_DIST_TOL)).abs() < 1e-9);
    }
    #[test]
    fn set_macrocycle_all_in_same_ring_14_bounds_uses_trans_plus_point_one_for_macrocycle_amide() {
        let mol = Ok::<_, ()>(fixed_topology("O=C1N(C)CCCCCCCC1", true))
            .expect("macrocyclic tertiary lactam");
        let (mut mmat, mut accum_data) = run_set13_bounds(&mol);

        let carbonyl = mol
            .atoms
            .iter()
            .enumerate()
            .find_map(|(idx, _)| is_carbonyl(&mol, idx).then_some(idx))
            .expect("carbonyl carbon");
        let nitrogen = neighbors_for_atom(&mol, carbonyl)
            .into_iter()
            .find(|&nbr| {
                let bond_idx = bond_between_idx_simple(&mol, carbonyl, nbr).expect("bond");
                mol.bonds[bond_idx].order() == BondOrder::Single
                    && mol.atoms[nbr].atomic_number() == 7
            })
            .expect("amide nitrogen");
        let atm4 = neighbors_for_atom(&mol, carbonyl)
            .into_iter()
            .find(|&nbr| {
                let bond_idx = bond_between_idx_simple(&mol, carbonyl, nbr).expect("bond");
                mol.bonds[bond_idx].order() == BondOrder::Single && nbr != nitrogen
            })
            .expect("ring carbon neighbor");
        let atom1 = neighbors_for_atom(&mol, nitrogen)
            .into_iter()
            .find(|&nbr| nbr != carbonyl && mol.atoms[nbr].atomic_number() == 6)
            .expect("preceding ring carbon");

        let bid1 = bond_between_idx_simple(&mol, atom1, nitrogen).expect("b1");
        let bid2 = bond_between_idx_simple(&mol, nitrogen, carbonyl).expect("b2");
        let bid3 = bond_between_idx_simple(&mol, carbonyl, atm4).expect("b3");
        set_macrocycle_all_in_same_ring_14_bounds(
            &mol,
            &fixed_valence(&mol),
            bid1,
            bid2,
            bid3,
            &mut accum_data,
            &mut mmat,
        )
        .expect("macrocycle all-in-same-ring bounds");

        let path = accum_data.paths14.last().expect("path");
        assert_eq!(path.kind, Path14Kind::Trans);
        let expected = compute_14_dist_trans(
            accum_data.bond_lengths[bid1],
            accum_data.bond_lengths[bid2],
            accum_data.bond_lengths[bid3],
            accum_data.get_bond_angle(mol.bonds.len(), bid1, bid2),
            accum_data.get_bond_angle(mol.bonds.len(), bid2, bid3),
        ) + 0.1;
        assert!((mmat.get_lower(atom1, atm4).unwrap() - (expected - GEN_DIST_TOL)).abs() < 1e-9);
        assert!((mmat.get_upper(atom1, atm4).unwrap() - (expected + GEN_DIST_TOL)).abs() < 1e-9);
    }
    #[test]
    fn set_macrocycle_all_in_same_ring_14_bounds_uses_other_for_plain_macrocycle_chain() {
        let mol = Ok::<_, ()>(fixed_topology("C1CCCCCCCCC1", true)).expect("cyclodecane");
        let (mut mmat, mut accum_data, _) = run_set14_same_ring_pass_only(&mol, true);
        let (bid1, bid2, bid3, _) =
            find_same_ring_dispatch_triple(&mol, true).expect("macrocycle path");

        let before_paths = accum_data.paths14.len();
        set_macrocycle_all_in_same_ring_14_bounds(
            &mol,
            &fixed_valence(&mol),
            bid1,
            bid2,
            bid3,
            &mut accum_data,
            &mut mmat,
        )
        .expect("macrocycle all-in-same-ring bounds");

        let path = accum_data.paths14.get(before_paths).expect("new path");
        assert_eq!(path.kind, Path14Kind::Other);
    }
    #[test]
    fn set_chain_14_bounds_uses_defined_double_bond_stereo_for_alkenes() {
        let trans = Ok::<_, ()>(fixed_topology("C/C=C/C", true)).expect("trans alkene");
        let (mut trans_mmat, mut trans_accum) = run_set13_bounds(&trans);
        let t_bid1 = bond_between_idx_simple(&trans, 0, 1).expect("0-1");
        let t_bid2 = bond_between_idx_simple(&trans, 1, 2).expect("1-2");
        let t_bid3 = bond_between_idx_simple(&trans, 2, 3).expect("2-3");
        set_chain_14_bounds(
            &trans,
            &fixed_valence(&trans),
            t_bid1,
            t_bid2,
            t_bid3,
            &mut trans_accum,
            &mut trans_mmat,
            false,
        )
        .expect("chain trans bounds");
        let trans_path = trans_accum.paths14.last().expect("trans path");
        assert_eq!(trans_path.kind, Path14Kind::Trans);
        let trans_expected = compute_14_dist_trans(
            trans_accum.bond_lengths[t_bid1],
            trans_accum.bond_lengths[t_bid2],
            trans_accum.bond_lengths[t_bid3],
            trans_accum.get_bond_angle(trans.bonds.len(), t_bid1, t_bid2),
            trans_accum.get_bond_angle(trans.bonds.len(), t_bid2, t_bid3),
        );
        assert!(
            (trans_mmat.get_lower(0, 3).unwrap() - (trans_expected - GEN_DIST_TOL)).abs() < 1e-9
        );
        assert!(
            (trans_mmat.get_upper(0, 3).unwrap() - (trans_expected + GEN_DIST_TOL)).abs() < 1e-9
        );

        let cis = Ok::<_, ()>(fixed_topology("C/C=C\\C", true)).expect("cis alkene");
        let (mut cis_mmat, mut cis_accum) = run_set13_bounds(&cis);
        let c_bid1 = bond_between_idx_simple(&cis, 0, 1).expect("0-1");
        let c_bid2 = bond_between_idx_simple(&cis, 1, 2).expect("1-2");
        let c_bid3 = bond_between_idx_simple(&cis, 2, 3).expect("2-3");
        set_chain_14_bounds(
            &cis,
            &fixed_valence(&cis),
            c_bid1,
            c_bid2,
            c_bid3,
            &mut cis_accum,
            &mut cis_mmat,
            false,
        )
        .expect("chain cis bounds");
        let cis_path = cis_accum.paths14.last().expect("cis path");
        assert_eq!(cis_path.kind, Path14Kind::Cis);
        let cis_expected = compute_14_dist_cis(
            cis_accum.bond_lengths[c_bid1],
            cis_accum.bond_lengths[c_bid2],
            cis_accum.bond_lengths[c_bid3],
            cis_accum.get_bond_angle(cis.bonds.len(), c_bid1, c_bid2),
            cis_accum.get_bond_angle(cis.bonds.len(), c_bid2, c_bid3),
        );
        assert!((cis_mmat.get_lower(0, 3).unwrap() - (cis_expected - GEN_DIST_TOL)).abs() < 1e-9);
        assert!((cis_mmat.get_upper(0, 3).unwrap() - (cis_expected + GEN_DIST_TOL)).abs() < 1e-9);
    }
    #[test]
    fn set_chain_14_bounds_uses_ss_special_case() {
        let mol = Ok::<_, ()>(fixed_topology("CSSC", true)).expect("disulfide");
        let (mut mmat, mut accum_data) = run_set13_bounds(&mol);
        let bid1 = bond_between_idx_simple(&mol, 0, 1).expect("0-1");
        let bid2 = bond_between_idx_simple(&mol, 1, 2).expect("1-2");
        let bid3 = bond_between_idx_simple(&mol, 2, 3).expect("2-3");
        set_chain_14_bounds(
            &mol,
            &fixed_valence(&mol),
            bid1,
            bid2,
            bid3,
            &mut accum_data,
            &mut mmat,
            false,
        )
        .expect("chain S-S bounds");

        let path = accum_data.paths14.last().expect("path");
        assert_eq!(path.kind, Path14Kind::Other);
        let expected = compute_14_dist_3d(
            accum_data.bond_lengths[bid1],
            accum_data.bond_lengths[bid2],
            accum_data.bond_lengths[bid3],
            accum_data.get_bond_angle(mol.bonds.len(), bid1, bid2),
            accum_data.get_bond_angle(mol.bonds.len(), bid2, bid3),
            std::f64::consts::PI / 2.0,
        );
        assert!((mmat.get_lower(0, 3).unwrap() - (expected - GEN_DIST_TOL)).abs() < 1e-9);
        assert!((mmat.get_upper(0, 3).unwrap() - (expected + GEN_DIST_TOL)).abs() < 1e-9);
    }
    #[test]
    fn set_chain_14_bounds_honors_force_trans_amides_for_secondary_amide_h_paths() {
        let hydrogen = 0;
        let nitrogen = 1;
        let carbonyl = 2;
        let oxygen = 3;
        let n_methyl = 4;
        let carbonyl_methyl = 5;
        let atoms = vec![
            AtomSpec::new(Element::H),
            AtomSpec::new(Element::N).with_hybridization(Hybridization::Sp2),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
            AtomSpec::new(Element::O).with_hybridization(Hybridization::Sp2),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
        ]
        .into_iter()
        .enumerate()
        .map(|(i, s)| Atom::from_spec(AtomId::new(i), s))
        .collect();
        let bonds = [
            (hydrogen, nitrogen, BondOrder::Single),
            (nitrogen, carbonyl, BondOrder::Single),
            (nitrogen, n_methyl, BondOrder::Single),
            (carbonyl, oxygen, BondOrder::Double),
            (carbonyl, carbonyl_methyl, BondOrder::Single),
        ]
        .into_iter()
        .enumerate()
        .map(|(i, (a, b, o))| {
            Bond::from_spec(
                BondId::new(i),
                BondSpec::new(AtomId::new(a), AtomId::new(b), o),
            )
        })
        .collect();
        let mol = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![])
            .expect("secondary amide with explicit N-H");
        let (mut mmat, mut accum_data) = run_set12_bounds(&mol);

        let bid_hn = bond_between_idx_simple(&mol, hydrogen, nitrogen).expect("h-n");
        let bid_nc = bond_between_idx_simple(&mol, nitrogen, carbonyl).expect("n-c");
        let bid_co = bond_between_idx_simple(&mol, carbonyl, oxygen).expect("c-o");
        let bid_cm = bond_between_idx_simple(&mol, carbonyl, carbonyl_methyl).expect("c-m");
        let nb = mol.bonds.len();
        accum_data.set_bond_adj(nb, bid_hn, bid_nc, nitrogen as i32);
        accum_data.set_bond_adj(nb, bid_nc, bid_co, carbonyl as i32);
        accum_data.set_bond_adj(nb, bid_nc, bid_cm, carbonyl as i32);
        accum_data.set_bond_angle(
            nb,
            bid_hn,
            bid_nc,
            ideal_bond_angle(&mol.atoms[nitrogen].hybridization(), None),
        );
        accum_data.set_bond_angle(
            nb,
            bid_nc,
            bid_co,
            ideal_bond_angle(&mol.atoms[carbonyl].hybridization(), None),
        );
        accum_data.set_bond_angle(
            nb,
            bid_nc,
            bid_cm,
            ideal_bond_angle(&mol.atoms[carbonyl].hybridization(), None),
        );

        set_chain_14_bounds(
            &mol,
            &fixed_valence(&mol),
            bid_hn,
            bid_nc,
            bid_co,
            &mut accum_data,
            &mut mmat,
            true,
        )
        .expect("chain amide14 bounds");
        let amide14_path = accum_data.paths14.last().expect("amide14 path");
        assert_eq!(amide14_path.kind, Path14Kind::Trans);
        let expected_14 = compute_14_dist_trans(
            accum_data.bond_lengths[bid_hn],
            accum_data.bond_lengths[bid_nc],
            accum_data.bond_lengths[bid_co],
            accum_data.get_bond_angle(mol.bonds.len(), bid_hn, bid_nc),
            accum_data.get_bond_angle(mol.bonds.len(), bid_nc, bid_co),
        );
        assert!(
            (mmat.get_lower(hydrogen, oxygen).unwrap() - (expected_14 - GEN_DIST_TOL)).abs() < 1e-9
        );
        assert!(
            (mmat.get_upper(hydrogen, oxygen).unwrap() - (expected_14 + GEN_DIST_TOL)).abs() < 1e-9
        );

        set_chain_14_bounds(
            &mol,
            &fixed_valence(&mol),
            bid_hn,
            bid_nc,
            bid_cm,
            &mut accum_data,
            &mut mmat,
            true,
        )
        .expect("chain amide15 bounds");
        let amide15_path = accum_data.paths14.last().expect("amide15 path");
        assert_eq!(amide15_path.kind, Path14Kind::Cis);
        let expected_15 = compute_14_dist_cis(
            accum_data.bond_lengths[bid_hn],
            accum_data.bond_lengths[bid_nc],
            accum_data.bond_lengths[bid_cm],
            accum_data.get_bond_angle(mol.bonds.len(), bid_hn, bid_nc),
            accum_data.get_bond_angle(mol.bonds.len(), bid_nc, bid_cm),
        );
        assert!(
            (mmat.get_lower(hydrogen, carbonyl_methyl).unwrap() - (expected_15 - GEN_DIST_TOL))
                .abs()
                < 1e-9
        );
        assert!(
            (mmat.get_upper(hydrogen, carbonyl_methyl).unwrap() - (expected_15 + GEN_DIST_TOL))
                .abs()
                < 1e-9
        );
    }
    #[test]
    fn record_14_path_marks_sp2_sp2_ring_paths_as_cis() {
        let mol = Ok::<_, ()>(fixed_topology("c1ccccc1", true)).expect("benzene");
        let (_mmat, mut accum_data) = run_set13_bounds(&mol);
        let bid1 = bond_between_idx_simple(&mol, 5, 0).expect("5-0");
        let bid2 = bond_between_idx_simple(&mol, 0, 1).expect("0-1");
        let bid3 = bond_between_idx_simple(&mol, 1, 2).expect("1-2");
        record_14_path(&mol, bid1, bid2, bid3, &mut accum_data).expect("record14Path");

        let path = accum_data.paths14.last().expect("path");
        assert_eq!(path.bid1, bid1);
        assert_eq!(path.bid2, bid2);
        assert_eq!(path.bid3, bid3);
        assert_eq!(path.kind, Path14Kind::Cis);
        assert!(has_path_flag(
            &accum_data.cis_paths,
            path14_id(mol.bonds.len(), bid1, bid2, bid3)
        ));
        assert!(has_path_flag(
            &accum_data.cis_paths,
            path14_id(mol.bonds.len(), bid3, bid2, bid1)
        ));
    }
    #[test]
    fn record_14_path_uses_other_for_non_sp2_path_without_cis_flags() {
        let mol = Ok::<_, ()>(fixed_topology("CCCC", true)).expect("butane");
        let (_mmat, mut accum_data) = run_set13_bounds(&mol);
        let bid1 = bond_between_idx_simple(&mol, 0, 1).expect("0-1");
        let bid2 = bond_between_idx_simple(&mol, 1, 2).expect("1-2");
        let bid3 = bond_between_idx_simple(&mol, 2, 3).expect("2-3");

        record_14_path(&mol, bid1, bid2, bid3, &mut accum_data).expect("record14Path");

        let path = accum_data.paths14.last().expect("path");
        assert_eq!(path.kind, Path14Kind::Other);
        assert!(!has_path_flag(
            &accum_data.cis_paths,
            path14_id(mol.bonds.len(), bid1, bid2, bid3)
        ));
        assert!(!has_path_flag(
            &accum_data.cis_paths,
            path14_id(mol.bonds.len(), bid3, bid2, bid1)
        ));
    }
    #[test]
    fn set_in_ring_14_bounds_prefers_cis_for_small_sp2_ring_paths() {
        let mol = Ok::<_, ()>(fixed_topology("c1ccccc1", true)).expect("benzene");
        let (mut mmat, mut accum_data) = run_set13_bounds(&mol);
        let dmat = flatten_topological_distances(&mol);
        let bid1 = bond_between_idx_simple(&mol, 5, 0).expect("5-0");
        let bid2 = bond_between_idx_simple(&mol, 0, 1).expect("0-1");
        let bid3 = bond_between_idx_simple(&mol, 1, 2).expect("1-2");
        let rinfo = fixed_rings(&mol);
        set_in_ring_14_bounds(
            &mol,
            bid1,
            bid2,
            bid3,
            &mut accum_data,
            &mut mmat,
            &dmat,
            6,
            &rinfo,
        )
        .expect("in-ring cis bounds");

        let path = accum_data.paths14.last().expect("path");
        assert_eq!(path.kind, Path14Kind::Cis);
        assert!(has_path_flag(
            &accum_data.cis_paths,
            path14_id(mol.bonds.len(), bid1, bid2, bid3)
        ));

        let aid1 = 5usize;
        let aid4 = 2usize;
        let pid = aid1.min(aid4) * mol.atoms.len() + aid1.max(aid4);
        let expected = compute_14_dist_cis(
            accum_data.bond_lengths[bid1],
            accum_data.bond_lengths[bid2],
            accum_data.bond_lengths[bid3],
            accum_data.get_bond_angle(mol.bonds.len(), bid1, bid2),
            accum_data.get_bond_angle(mol.bonds.len(), bid2, bid3),
        );
        assert!(accum_data.visited14_bounds[pid]);
        assert!((mmat.get_lower(aid1, aid4).unwrap() - (expected - GEN_DIST_TOL)).abs() < 1e-9);
        assert!((mmat.get_upper(aid1, aid4).unwrap() - (expected + GEN_DIST_TOL)).abs() < 1e-9);
    }
    #[test]
    fn set_two_in_same_ring_14_bounds_uses_trans_for_sp2_external_substituent_path() {
        let mol = Ok::<_, ()>(fixed_topology("Cc1ccccc1", true)).expect("toluene");
        let (mut mmat, mut accum_data) = run_set13_bounds(&mol);
        let dmat = flatten_topological_distances(&mol);

        let bid_exo = bond_between_idx_simple(&mol, 0, 1).expect("exo bond");
        let bid_ring_12 = bond_between_idx_simple(&mol, 1, 2).expect("ring bond 1-2");
        let bid_ring_23 = bond_between_idx_simple(&mol, 2, 3).expect("ring bond 2-3");
        set_two_in_same_ring_14_bounds(
            &mol,
            bid_exo,
            bid_ring_12,
            bid_ring_23,
            &mut accum_data,
            &mut mmat,
            &dmat,
        )
        .expect("two-in-same-ring bounds");

        let path = accum_data.paths14.last().expect("path");
        assert_eq!(path.kind, Path14Kind::Trans);
        assert!(has_path_flag(
            &accum_data.trans_paths,
            path14_id(mol.bonds.len(), bid_exo, bid_ring_12, bid_ring_23)
        ));

        let aid1 = 0usize;
        let aid4 = 3usize;
        let expected = compute_14_dist_trans(
            accum_data.bond_lengths[bid_exo],
            accum_data.bond_lengths[bid_ring_12],
            accum_data.bond_lengths[bid_ring_23],
            accum_data.get_bond_angle(mol.bonds.len(), bid_exo, bid_ring_12),
            accum_data.get_bond_angle(mol.bonds.len(), bid_ring_12, bid_ring_23),
        );
        assert!((mmat.get_lower(aid1, aid4).unwrap() - (expected - GEN_DIST_TOL)).abs() < 1e-9);
        assert!((mmat.get_upper(aid1, aid4).unwrap() - (expected + GEN_DIST_TOL)).abs() < 1e-9);
    }
    #[test]
    fn diff_ring14_and_share_ring14_delegate_to_in_ring_helper() {
        let mol = Ok::<_, ()>(fixed_topology("c1ccccc1", true)).expect("benzene");
        let dmat = flatten_topological_distances(&mol);
        let bid1 = bond_between_idx_simple(&mol, 5, 0).expect("5-0");
        let bid2 = bond_between_idx_simple(&mol, 0, 1).expect("0-1");
        let bid3 = bond_between_idx_simple(&mol, 1, 2).expect("1-2");
        let rinfo = fixed_rings(&mol);

        let (mut base_mmat, mut base_accum) = run_set13_bounds(&mol);
        set_in_ring_14_bounds(
            &mol,
            bid1,
            bid2,
            bid3,
            &mut base_accum,
            &mut base_mmat,
            &dmat,
            0,
            &rinfo,
        )
        .expect("base in-ring bounds");

        let (mut diff_mmat, mut diff_accum) = run_set13_bounds(&mol);
        set_two_in_diff_ring_14_bounds(
            &mol,
            bid1,
            bid2,
            bid3,
            &mut diff_accum,
            &mut diff_mmat,
            &dmat,
            &rinfo,
        )
        .expect("diff-ring bounds");

        let (mut share_mmat, mut share_accum) = run_set13_bounds(&mol);
        set_share_ring_bond_14_bounds(
            &mol,
            bid1,
            bid2,
            bid3,
            &mut share_accum,
            &mut share_mmat,
            &dmat,
            &rinfo,
        )
        .expect("share-ring-bond bounds");

        assert_eq!(
            diff_accum.paths14.last().map(|p| p.kind),
            base_accum.paths14.last().map(|p| p.kind)
        );
        assert_eq!(
            share_accum.paths14.last().map(|p| p.kind),
            base_accum.paths14.last().map(|p| p.kind)
        );
        assert!(
            (diff_mmat.get_lower(5, 2).unwrap() - base_mmat.get_lower(5, 2).unwrap()).abs() < 1e-9
        );
        assert!(
            (diff_mmat.get_upper(5, 2).unwrap() - base_mmat.get_upper(5, 2).unwrap()).abs() < 1e-9
        );
        assert!(
            (share_mmat.get_lower(5, 2).unwrap() - base_mmat.get_lower(5, 2).unwrap()).abs() < 1e-9
        );
        assert!(
            (share_mmat.get_upper(5, 2).unwrap() - base_mmat.get_upper(5, 2).unwrap()).abs() < 1e-9
        );
    }
    #[test]
    fn set_14_bounds_entrypoint_same_ring_matches_direct_helper() {
        let mol = Ok::<_, ()>(fixed_topology("c1ccccc1", true)).expect("benzene");
        let (mmat, accum_data, dmat) = run_set14_bounds(&mol, false, false);
        let (bid1, bid2, bid3, ring_size) =
            find_same_ring_dispatch_triple(&mol, false).expect("same-ring triple");
        let rinfo = fixed_rings(&mol);

        let (mut direct_mmat, mut direct_accum) = run_set13_bounds(&mol);
        set_in_ring_14_bounds(
            &mol,
            bid1,
            bid2,
            bid3,
            &mut direct_accum,
            &mut direct_mmat,
            &dmat,
            ring_size,
            &rinfo,
        )
        .expect("direct same-ring bounds");

        let atm2 = bond_pair_shared_atom(&mol, &direct_accum, bid1, bid2).expect("shared atom");
        let atm3 = bond_pair_shared_atom(&mol, &direct_accum, bid2, bid3).expect("shared atom");
        let aid1 = if mol.bonds[bid1].begin().index() == atm2 {
            mol.bonds[bid1].end().index()
        } else {
            mol.bonds[bid1].begin().index()
        };
        let aid4 = if mol.bonds[bid3].begin().index() == atm3 {
            mol.bonds[bid3].end().index()
        } else {
            mol.bonds[bid3].begin().index()
        };
        assert!(
            (mmat.get_lower(aid1, aid4).unwrap() - direct_mmat.get_lower(aid1, aid4).unwrap())
                .abs()
                < 1e-9
        );
        assert!(
            (mmat.get_upper(aid1, aid4).unwrap() - direct_mmat.get_upper(aid1, aid4).unwrap())
                .abs()
                < 1e-9
        );
        assert!(accum_data.visited14_bounds[aid1.min(aid4) * mol.atoms.len() + aid1.max(aid4)]);
    }
    #[test]
    fn set_14_bounds_entrypoint_two_same_ring_matches_direct_helper() {
        let mol = Ok::<_, ()>(fixed_topology("Cc1ccccc1", true)).expect("toluene");
        let (mmat, _accum_data, dmat) = run_set14_bounds(&mol, false, false);
        let (bid1, bid2, bid3) = find_dispatch_triple(&mol, Set14DispatchCase::TwoSameRing, false)
            .expect("two-same-ring triple");

        let (mut direct_mmat, mut direct_accum) = run_set13_bounds(&mol);
        set_two_in_same_ring_14_bounds(
            &mol,
            bid1,
            bid2,
            bid3,
            &mut direct_accum,
            &mut direct_mmat,
            &dmat,
        )
        .expect("direct two-same-ring bounds");

        let atm2 = bond_pair_shared_atom(&mol, &direct_accum, bid1, bid2).expect("shared atom");
        let atm3 = bond_pair_shared_atom(&mol, &direct_accum, bid2, bid3).expect("shared atom");
        let aid1 = if mol.bonds[bid1].begin().index() == atm2 {
            mol.bonds[bid1].end().index()
        } else {
            mol.bonds[bid1].begin().index()
        };
        let aid4 = if mol.bonds[bid3].begin().index() == atm3 {
            mol.bonds[bid3].end().index()
        } else {
            mol.bonds[bid3].begin().index()
        };
        assert!(
            (mmat.get_lower(aid1, aid4).unwrap() - direct_mmat.get_lower(aid1, aid4).unwrap())
                .abs()
                < 1e-9
        );
        assert!(
            (mmat.get_upper(aid1, aid4).unwrap() - direct_mmat.get_upper(aid1, aid4).unwrap())
                .abs()
                < 1e-9
        );
    }
    #[test]
    fn set_14_bounds_entrypoint_two_diff_ring_matches_direct_helper() {
        let mol = Ok::<_, ()>(fixed_topology("C1CCC2(CC1)CCC3CCCCC23", true))
            .expect("two-diff-ring polycycle");
        let (mmat, _accum_data, dmat) = run_set14_bounds(&mol, false, false);
        let (bid1, bid2, bid3) = find_dispatch_triple(&mol, Set14DispatchCase::TwoDiffRing, false)
            .expect("two-diff-ring triple");
        let rinfo = fixed_rings(&mol);

        let (mut direct_mmat, mut direct_accum) = run_set13_bounds(&mol);
        set_two_in_diff_ring_14_bounds(
            &mol,
            bid1,
            bid2,
            bid3,
            &mut direct_accum,
            &mut direct_mmat,
            &dmat,
            &rinfo,
        )
        .expect("direct two-diff-ring bounds");

        let atm2 = bond_pair_shared_atom(&mol, &direct_accum, bid1, bid2).expect("shared atom");
        let atm3 = bond_pair_shared_atom(&mol, &direct_accum, bid2, bid3).expect("shared atom");
        let aid1 = if mol.bonds[bid1].begin().index() == atm2 {
            mol.bonds[bid1].end().index()
        } else {
            mol.bonds[bid1].begin().index()
        };
        let aid4 = if mol.bonds[bid3].begin().index() == atm3 {
            mol.bonds[bid3].end().index()
        } else {
            mol.bonds[bid3].begin().index()
        };
        assert!(
            (mmat.get_lower(aid1, aid4).unwrap() - direct_mmat.get_lower(aid1, aid4).unwrap())
                .abs()
                < 1e-9
        );
        assert!(
            (mmat.get_upper(aid1, aid4).unwrap() - direct_mmat.get_upper(aid1, aid4).unwrap())
                .abs()
                < 1e-9
        );
    }
    #[test]
    fn set_14_bounds_entrypoint_share_ring_bond_matches_direct_helper() {
        let mol = Ok::<_, ()>(fixed_topology("Cc1ccccc1C", true)).expect("xylene");
        let (mmat, _accum_data, dmat) = run_set14_bounds(&mol, false, false);
        let (bid1, bid2, bid3) =
            find_dispatch_triple(&mol, Set14DispatchCase::ShareRingBond, false)
                .expect("share-ring-bond triple");
        let rinfo = fixed_rings(&mol);

        let (mut direct_mmat, mut direct_accum) = run_set13_bounds(&mol);
        set_share_ring_bond_14_bounds(
            &mol,
            bid1,
            bid2,
            bid3,
            &mut direct_accum,
            &mut direct_mmat,
            &dmat,
            &rinfo,
        )
        .expect("direct share-ring-bond bounds");

        let atm2 = bond_pair_shared_atom(&mol, &direct_accum, bid1, bid2).expect("shared atom");
        let atm3 = bond_pair_shared_atom(&mol, &direct_accum, bid2, bid3).expect("shared atom");
        let aid1 = if mol.bonds[bid1].begin().index() == atm2 {
            mol.bonds[bid1].end().index()
        } else {
            mol.bonds[bid1].begin().index()
        };
        let aid4 = if mol.bonds[bid3].begin().index() == atm3 {
            mol.bonds[bid3].end().index()
        } else {
            mol.bonds[bid3].begin().index()
        };
        assert!(
            (mmat.get_lower(aid1, aid4).unwrap() - direct_mmat.get_lower(aid1, aid4).unwrap())
                .abs()
                < 1e-9
        );
        assert!(
            (mmat.get_upper(aid1, aid4).unwrap() - direct_mmat.get_upper(aid1, aid4).unwrap())
                .abs()
                < 1e-9
        );
    }
    #[test]
    fn set_14_bounds_entrypoint_chain_matches_direct_helper() {
        let mol = Ok::<_, ()>(fixed_topology("CSSC", true)).expect("disulfide");
        let (mmat, _accum_data, _dmat) = run_set14_bounds(&mol, false, false);
        let (bid1, bid2, bid3) =
            find_dispatch_triple(&mol, Set14DispatchCase::Chain, false).expect("chain triple");

        let (mut direct_mmat, mut direct_accum) = run_set13_bounds(&mol);
        set_chain_14_bounds(
            &mol,
            &fixed_valence(&mol),
            bid1,
            bid2,
            bid3,
            &mut direct_accum,
            &mut direct_mmat,
            false,
        )
        .expect("direct chain bounds");

        let atm2 = bond_pair_shared_atom(&mol, &direct_accum, bid1, bid2).expect("shared atom");
        let atm3 = bond_pair_shared_atom(&mol, &direct_accum, bid2, bid3).expect("shared atom");
        let aid1 = if mol.bonds[bid1].begin().index() == atm2 {
            mol.bonds[bid1].end().index()
        } else {
            mol.bonds[bid1].begin().index()
        };
        let aid4 = if mol.bonds[bid3].begin().index() == atm3 {
            mol.bonds[bid3].end().index()
        } else {
            mol.bonds[bid3].begin().index()
        };
        assert!(
            (mmat.get_lower(aid1, aid4).unwrap() - direct_mmat.get_lower(aid1, aid4).unwrap())
                .abs()
                < 1e-9
        );
        assert!(
            (mmat.get_upper(aid1, aid4).unwrap() - direct_mmat.get_upper(aid1, aid4).unwrap())
                .abs()
                < 1e-9
        );
    }
    #[test]
    fn set_14_bounds_entrypoint_macrocycle_matches_direct_helper() {
        let mol = Ok::<_, ()>(fixed_topology("O=C1N(C)CCCCCCCC1", true))
            .expect("macrocyclic tertiary lactam");
        let (mmat, _accum_data, _dmat) = run_set14_same_ring_pass_only(&mol, true);
        let (bid1, bid2, bid3, _ring_size) =
            find_same_ring_dispatch_triple(&mol, true).expect("macrocycle same-ring triple");

        let (mut direct_mmat, mut direct_accum) = run_set13_bounds(&mol);
        set_macrocycle_all_in_same_ring_14_bounds(
            &mol,
            &fixed_valence(&mol),
            bid1,
            bid2,
            bid3,
            &mut direct_accum,
            &mut direct_mmat,
        )
        .expect("direct macrocycle bounds");
        let atm2 = bond_pair_shared_atom(&mol, &direct_accum, bid1, bid2).expect("shared atom");
        let atm3 = bond_pair_shared_atom(&mol, &direct_accum, bid2, bid3).expect("shared atom");
        let atom1 = if mol.bonds[bid1].begin().index() == atm2 {
            mol.bonds[bid1].end().index()
        } else {
            mol.bonds[bid1].begin().index()
        };
        let atm4 = if mol.bonds[bid3].begin().index() == atm3 {
            mol.bonds[bid3].end().index()
        } else {
            mol.bonds[bid3].begin().index()
        };
        assert!(
            (mmat.get_lower(atom1, atm4).unwrap() - direct_mmat.get_lower(atom1, atm4).unwrap())
                .abs()
                < 1e-9
        );
        assert!(
            (mmat.get_upper(atom1, atm4).unwrap() - direct_mmat.get_upper(atom1, atm4).unwrap())
                .abs()
                < 1e-9
        );
    }
    #[test]
    fn get_atom_stereo_preserves_stereo_when_stereo_atoms_match_query_order() {
        let a0 = AtomId::new(0);
        let a1 = AtomId::new(1);
        let a2 = AtomId::new(2);
        let a3 = AtomId::new(3);
        let bond = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(a1, a2, BondOrder::Double)
                .with_stereo(BondStereo::Cis)
                .with_stereo_atoms(a0, a3),
        );

        assert_eq!(
            get_atom_stereo(&bond, a0.index(), a3.index()),
            BondStereo::Cis
        );
    }
    #[test]
    fn get_atom_stereo_flips_stereo_when_stereo_atoms_reverse_query_order() {
        let a0 = AtomId::new(0);
        let a1 = AtomId::new(1);
        let a2 = AtomId::new(2);
        let a3 = AtomId::new(3);
        let bond = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(a1, a2, BondOrder::Double)
                .with_stereo(BondStereo::Cis)
                .with_stereo_atoms(a0, a3),
        );

        assert_eq!(
            get_atom_stereo(&bond, a3.index(), a0.index()),
            BondStereo::Cis
        );
    }
    #[test]
    fn get_atom_stereo_flips_stereo_when_single_end_mismatches() {
        let a0 = AtomId::new(0);
        let a1 = AtomId::new(1);
        let a2 = AtomId::new(2);
        let a3 = AtomId::new(3);
        let bond = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(a1, a2, BondOrder::Double)
                .with_stereo(BondStereo::Cis)
                .with_stereo_atoms(a0, a3),
        );

        assert_eq!(
            get_atom_stereo(&bond, a3.index(), a3.index()),
            BondStereo::Trans
        );
    }

    fn matrix_rows(bounds: &BoundsMatrix) -> Vec<Vec<f64>> {
        (0..bounds.dimension())
            .map(|i| {
                (0..bounds.dimension())
                    .map(|j| bounds.get_val(i, j).unwrap())
                    .collect()
            })
            .collect()
    }
    fn fixture_chemistry(
        mol: &TopologyBlock,
    ) -> (RingInfo, ValenceAssignment, Vec<Hybridization>, Vec<bool>) {
        let valence = fixed_valence(mol);
        let conjugated = assign_conjugation_flags(mol, &valence).expect("source conjugation");
        let hyb = assign_hybridization_with_conjugation(mol, &valence, &conjugated)
            .expect("source hybridization")
            .values;
        (fixed_rings(mol), valence, hyb, conjugated)
    }
    fn run_set15_bounds(
        mol: &TopologyBlock,
        use_macrocycle_14config: bool,
        force_trans_amides: bool,
    ) -> (BoundsMatrix, ComputedData, Vec<f64>) {
        let (mut mmat, mut accum_data, dmat) =
            run_set14_bounds(mol, use_macrocycle_14config, force_trans_amides);
        set_15_bounds(mol, &mut mmat, &mut accum_data, &dmat).expect("set15Bounds");
        (mmat, accum_data, dmat)
    }
    fn run_set_topol_bounds(
        mol: &TopologyBlock,
        set15bounds: bool,
        scale_vdw: bool,
        use_macrocycle_14config: bool,
        force_trans_amides: bool,
        set14bounds: bool,
        set13bounds: bool,
    ) -> BoundsMatrix {
        let (rings, valence, hybridizations, conjugated) = fixture_chemistry(mol);

        let mut mmat = initialized_bounds(mol.atoms.len());
        set_topol_bounds(
            mol,
            &rings,
            &valence,
            &hybridizations,
            &conjugated,
            &mut mmat,
            set15bounds,
            scale_vdw,
            use_macrocycle_14config,
            force_trans_amides,
            set14bounds,
            set13bounds,
        )
        .expect("setTopolBounds");
        mmat
    }
    fn run_set_topol_bounds_with_outputs(
        mol: &TopologyBlock,
        set15bounds: bool,
        scale_vdw: bool,
        use_macrocycle_14config: bool,
        force_trans_amides: bool,
        set14bounds: bool,
        set13bounds: bool,
    ) -> (BoundsMatrix, Vec<(i32, i32)>, Vec<Vec<i32>>) {
        let (rings, valence, hybridizations, conjugated) = fixture_chemistry(mol);

        let mut mmat = initialized_bounds(mol.atoms.len());
        let mut bonds = Vec::new();
        let mut angles = Vec::new();
        set_topol_bounds_with_outputs(
            mol,
            &rings,
            &valence,
            &hybridizations,
            &conjugated,
            &mut mmat,
            &mut bonds,
            &mut angles,
            set15bounds,
            scale_vdw,
            use_macrocycle_14config,
            force_trans_amides,
            set14bounds,
            set13bounds,
        )
        .expect("setTopolBounds with outputs");
        (mmat, bonds, angles)
    }

    #[test]
    fn collect_bonds_and_angles_flags_triple_bond_paths_like_rdkit() {
        let mol = Ok::<_, ()>(fixed_topology("CC#N", true)).expect("acetonitrile");
        let mut bonds = Vec::new();
        let mut angles = Vec::new();

        collect_bonds_and_angles(&mol, &mut bonds, &mut angles);

        assert_eq!(bonds, vec![(0, 1), (1, 2)]);
        assert_eq!(angles, vec![vec![0, 1, 2, 1]]);
    }
    #[test]
    fn collect_bonds_and_angles_flags_consecutive_double_bonds_like_rdkit() {
        let mol = Ok::<_, ()>(fixed_topology("C=C=C", true)).expect("allene");
        let mut bonds = Vec::new();
        let mut angles = Vec::new();

        collect_bonds_and_angles(&mol, &mut bonds, &mut angles);

        assert_eq!(bonds, vec![(0, 1), (1, 2)]);
        assert_eq!(angles, vec![vec![0, 1, 2, 1]]);
    }
    #[test]
    fn set_lower_bound_vdw_scales_15_16_and_longer_paths_like_rdkit() {
        let atoms = (0..7)
            .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
            .collect();
        let mol =
            TopologyBlock::try_from_parts(atoms, vec![], vec![], vec![]).expect("seven carbons");
        let mut mmat = initialized_bounds(mol.atoms.len());
        let mut dmat = vec![0.0; mol.atoms.len() * mol.atoms.len()];
        dmat[4 * mol.atoms.len()] = 4.0;
        dmat[5 * mol.atoms.len()] = 5.0;
        dmat[6 * mol.atoms.len()] = 6.0;

        set_lower_bound_vdw(&mol, &mut mmat, true, &dmat).expect("setLowerBoundVDW");

        let vdw_sum = vdw_radius(6) + vdw_radius(6);
        assert!((mmat.get_lower(4, 0).unwrap() - (VDW_SCALE_15 * vdw_sum)).abs() < 1e-9);
        assert!(
            (mmat.get_lower(5, 0).unwrap()
                - ((VDW_SCALE_15 + 0.5 * (1.0 - VDW_SCALE_15)) * vdw_sum))
                .abs()
                < 1e-9
        );
        assert!((mmat.get_lower(6, 0).unwrap() - vdw_sum).abs() < 1e-9);
    }
    #[test]
    fn set_lower_bound_vdw_uses_hbond_floor_before_vdw_scaling() {
        let atoms = [Element::N, Element::H, Element::O]
            .into_iter()
            .enumerate()
            .map(|(i, e)| Atom::from_spec(AtomId::new(i), AtomSpec::new(e)))
            .collect();
        let bonds = vec![Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        )];
        let mol = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).expect("H-N ... O");
        let mut mmat = initialized_bounds(mol.atoms.len());
        let mut dmat = vec![0.0; mol.atoms.len() * mol.atoms.len()];
        dmat[2 * mol.atoms.len() + 1] = 6.0;
        dmat[1 * mol.atoms.len() + 2] = 6.0;

        set_lower_bound_vdw(&mol, &mut mmat, true, &dmat).expect("setLowerBoundVDW");

        assert!((mmat.get_lower(2, 1).unwrap() - H_BOND_LENGTH).abs() < 1e-9);
    }

    #[test]
    fn compute_15_dist_helpers_match_rdkit_cis_and_trans_formulas() {
        let d1: f64 = 1.41;
        let d2: f64 = 1.52;
        let d3: f64 = 1.38;
        let d4: f64 = 1.47;
        let ang12: f64 = 1.91;
        let ang23: f64 = 2.04;
        let ang34: f64 = 1.88;

        let cis_dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
        let cis_dy14 = d3 * ang23.sin() - d1 * ang12.sin();
        let cis_d14 = (cis_dx14 * cis_dx14 + cis_dy14 * cis_dy14).sqrt();
        let cis_cval =
            ((d3 - d2 * ang23.cos() + d1 * (ang12 + ang23).cos()) / cis_d14).clamp(-1.0, 1.0);
        let cis_ang143 = cis_cval.acos();
        let expected_cis_cis = compute_13_dist(cis_d14, d4, ang34 - cis_ang143);
        let expected_cis_trans = compute_13_dist(cis_d14, d4, ang34 + cis_ang143);

        let trans_dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
        let trans_dy14 = d3 * ang23.sin() + d1 * ang12.sin();
        let trans_d14 = (trans_dx14 * trans_dx14 + trans_dy14 * trans_dy14).sqrt();
        let trans_cval =
            ((d3 - d2 * ang23.cos() + d1 * (ang12 - ang23).cos()) / trans_d14).clamp(-1.0, 1.0);
        let trans_ang143 = trans_cval.acos();
        let expected_trans_trans = compute_13_dist(trans_d14, d4, ang34 + trans_ang143);
        let expected_trans_cis = compute_13_dist(trans_d14, d4, ang34 - trans_ang143);

        assert!(
            (compute_15_dist_cis_cis(d1, d2, d3, d4, ang12, ang23, ang34) - expected_cis_cis).abs()
                < 1e-12
        );
        assert!(
            (compute_15_dist_cis_trans(d1, d2, d3, d4, ang12, ang23, ang34) - expected_cis_trans)
                .abs()
                < 1e-12
        );
        assert!(
            (compute_15_dist_trans_trans(d1, d2, d3, d4, ang12, ang23, ang34)
                - expected_trans_trans)
                .abs()
                < 1e-12
        );
        assert!(
            (compute_15_dist_trans_cis(d1, d2, d3, d4, ang12, ang23, ang34) - expected_trans_cis)
                .abs()
                < 1e-12
        );
    }
    #[test]
    fn compute_15_dist_helpers_support_rdkit_reverse_argument_order() {
        let d1: f64 = 1.41;
        let d2: f64 = 1.52;
        let d3: f64 = 1.38;
        let d4: f64 = 1.47;
        let ang12: f64 = 1.91;
        let ang23: f64 = 2.04;
        let ang34: f64 = 1.88;

        let reversed_cis_cis = compute_15_dist_cis_cis(d4, d3, d2, d1, ang34, ang23, ang12);
        let reversed_cis_trans = compute_15_dist_cis_trans(d4, d3, d2, d1, ang34, ang23, ang12);
        let reversed_trans_cis = compute_15_dist_trans_cis(d4, d3, d2, d1, ang34, ang23, ang12);
        let reversed_trans_trans = compute_15_dist_trans_trans(d4, d3, d2, d1, ang34, ang23, ang12);

        assert_ne!(
            compute_15_dist_cis_cis(d1, d2, d3, d4, ang12, ang23, ang34),
            reversed_cis_cis
        );
        assert_ne!(
            compute_15_dist_cis_trans(d1, d2, d3, d4, ang12, ang23, ang34),
            reversed_cis_trans
        );
        assert_ne!(
            compute_15_dist_trans_cis(d1, d2, d3, d4, ang12, ang23, ang34),
            reversed_trans_cis
        );
        assert_ne!(
            compute_15_dist_trans_trans(d1, d2, d3, d4, ang12, ang23, ang34),
            reversed_trans_trans
        );
    }
    #[test]
    fn set_15_bounds_helper_returns_immediately_for_visited_14_pair() {
        let mol = Ok::<_, ()>(fixed_topology("CCCCC", true)).expect("pentane");
        let (mut mmat, mut accum_data) = run_set13_bounds(&mol);
        let dmat = flatten_topological_distances(&mol);
        let nb = mol.bonds.len();
        let na = mol.atoms.len();
        let bid1 = bond_between_idx_simple(&mol, 0, 1).expect("0-1");
        let bid2 = bond_between_idx_simple(&mol, 1, 2).expect("1-2");
        let bid3 = bond_between_idx_simple(&mol, 2, 3).expect("2-3");
        let pid = 0usize * na + 4usize;
        accum_data.visited14_bounds[pid] = true;
        let before_lower = mmat.get_lower(0, 4).unwrap();
        let before_upper = mmat.get_upper(0, 4).unwrap();

        set_15_bounds_helper(
            &mol,
            &mut mmat,
            &mut accum_data,
            &dmat,
            nb,
            na,
            bid1,
            bid2,
            bid3,
            Path14Kind::Other,
        )
        .expect("set15BoundsHelper visited skip");

        assert_eq!(mmat.get_lower(0, 4).unwrap(), before_lower);
        assert_eq!(mmat.get_upper(0, 4).unwrap(), before_upper);
        assert!(!accum_data.set15_atoms[0 * na + 4]);
        assert!(!accum_data.set15_atoms[4 * na + 0]);
    }
    #[test]
    fn set_15_bounds_helper_uses_vdw_fallback_and_marks_set15_atoms_for_unknown_path() {
        let mol = Ok::<_, ()>(fixed_topology("CCCCC", true)).expect("pentane");
        let (mut mmat, mut accum_data) = run_set13_bounds(&mol);
        let dmat = flatten_topological_distances(&mol);
        let nb = mol.bonds.len();
        let na = mol.atoms.len();
        let bid1 = bond_between_idx_simple(&mol, 0, 1).expect("0-1");
        let bid2 = bond_between_idx_simple(&mol, 1, 2).expect("1-2");
        let bid3 = bond_between_idx_simple(&mol, 2, 3).expect("2-3");

        set_15_bounds_helper(
            &mol,
            &mut mmat,
            &mut accum_data,
            &dmat,
            nb,
            na,
            bid1,
            bid2,
            bid3,
            Path14Kind::Other,
        )
        .expect("set15BoundsHelper vdw fallback");

        let expected_lower = VDW_SCALE_15 * (vdw_radius(6) + vdw_radius(6));
        assert!((mmat.get_lower(0, 4).unwrap() - expected_lower).abs() < 1e-12);
        assert_eq!(mmat.get_upper(0, 4).unwrap(), MAX_UPPER);
        assert!(accum_data.set15_atoms[0 * na + 4]);
        assert!(accum_data.set15_atoms[4 * na + 0]);
    }
    #[test]
    fn set_15_bounds_helper_uses_reversed_other_branch_formula_for_cis_path() {
        let mol = Ok::<_, ()>(fixed_topology("CCCCC", true)).expect("pentane");
        let (mut mmat, mut accum_data) = run_set13_bounds(&mol);
        let dmat = flatten_topological_distances(&mol);
        let nb = mol.bonds.len();
        let na = mol.atoms.len();
        let bid1 = bond_between_idx_simple(&mol, 0, 1).expect("0-1");
        let bid2 = bond_between_idx_simple(&mol, 1, 2).expect("1-2");
        let bid3 = bond_between_idx_simple(&mol, 2, 3).expect("2-3");
        let bid4 = bond_between_idx_simple(&mol, 3, 4).expect("3-4");
        let path_id = bid2 as u64 * nb as u64 * nb as u64 + bid3 as u64 * nb as u64 + bid4 as u64;
        record_path_flag(&mut accum_data.cis_paths, path_id);

        set_15_bounds_helper(
            &mol,
            &mut mmat,
            &mut accum_data,
            &dmat,
            nb,
            na,
            bid1,
            bid2,
            bid3,
            Path14Kind::Other,
        )
        .expect("set15BoundsHelper cis path");

        let d1 = accum_data.bond_lengths[bid1];
        let d2 = accum_data.bond_lengths[bid2];
        let d3 = accum_data.bond_lengths[bid3];
        let d4 = accum_data.bond_lengths[bid4];
        let ang12 = accum_data.get_bond_angle(nb, bid1, bid2);
        let ang23 = accum_data.get_bond_angle(nb, bid2, bid3);
        let ang34 = accum_data.get_bond_angle(nb, bid3, bid4);
        let expected_lower =
            compute_15_dist_cis_cis(d4, d3, d2, d1, ang34, ang23, ang12) - DIST15_TOL;
        let expected_upper =
            compute_15_dist_cis_trans(d4, d3, d2, d1, ang34, ang23, ang12) + DIST15_TOL;

        assert!((mmat.get_lower(0, 4).unwrap() - expected_lower).abs() < 1e-12);
        assert!((mmat.get_upper(0, 4).unwrap() - expected_upper).abs() < 1e-12);
        assert!(accum_data.set15_atoms[0 * na + 4]);
        assert!(accum_data.set15_atoms[4 * na + 0]);
    }
    #[test]
    fn set_15_bounds_helper_uses_reversed_other_branch_formula_for_trans_path() {
        let mol = Ok::<_, ()>(fixed_topology("CCCCC", true)).expect("pentane");
        let (mut mmat, mut accum_data) = run_set13_bounds(&mol);
        let dmat = flatten_topological_distances(&mol);
        let nb = mol.bonds.len();
        let na = mol.atoms.len();
        let bid1 = bond_between_idx_simple(&mol, 0, 1).expect("0-1");
        let bid2 = bond_between_idx_simple(&mol, 1, 2).expect("1-2");
        let bid3 = bond_between_idx_simple(&mol, 2, 3).expect("2-3");
        let bid4 = bond_between_idx_simple(&mol, 3, 4).expect("3-4");
        let path_id = bid2 as u64 * nb as u64 * nb as u64 + bid3 as u64 * nb as u64 + bid4 as u64;
        record_path_flag(&mut accum_data.trans_paths, path_id);

        set_15_bounds_helper(
            &mol,
            &mut mmat,
            &mut accum_data,
            &dmat,
            nb,
            na,
            bid1,
            bid2,
            bid3,
            Path14Kind::Other,
        )
        .expect("set15BoundsHelper trans path");

        let d1 = accum_data.bond_lengths[bid1];
        let d2 = accum_data.bond_lengths[bid2];
        let d3 = accum_data.bond_lengths[bid3];
        let d4 = accum_data.bond_lengths[bid4];
        let ang12 = accum_data.get_bond_angle(nb, bid1, bid2);
        let ang23 = accum_data.get_bond_angle(nb, bid2, bid3);
        let ang34 = accum_data.get_bond_angle(nb, bid3, bid4);
        let expected_lower =
            compute_15_dist_trans_cis(d4, d3, d2, d1, ang34, ang23, ang12) - DIST15_TOL;
        let expected_upper =
            compute_15_dist_trans_trans(d4, d3, d2, d1, ang34, ang23, ang12) + DIST15_TOL;

        assert!((mmat.get_lower(0, 4).unwrap() - expected_lower).abs() < 1e-12);
        assert!((mmat.get_upper(0, 4).unwrap() - expected_upper).abs() < 1e-12);
        assert!(!has_path_flag(&accum_data.cis_paths, path_id));
        assert!(accum_data.set15_atoms[0 * na + 4]);
        assert!(accum_data.set15_atoms[4 * na + 0]);
    }
    #[test]
    fn set_15_bounds_entrypoint_matches_two_helper_calls_for_single_path() {
        let mol = Ok::<_, ()>(fixed_topology("CCCCC", true)).expect("pentane");
        let dmat = flatten_topological_distances(&mol);
        let nb = mol.bonds.len();
        let na = mol.atoms.len();
        let bid1 = bond_between_idx_simple(&mol, 0, 1).expect("0-1");
        let bid2 = bond_between_idx_simple(&mol, 1, 2).expect("1-2");
        let bid3 = bond_between_idx_simple(&mol, 2, 3).expect("2-3");
        let bid4 = bond_between_idx_simple(&mol, 3, 4).expect("3-4");
        let path_id = bid2 as u64 * nb as u64 * nb as u64 + bid3 as u64 * nb as u64 + bid4 as u64;

        let (mut entry_mmat, mut entry_accum) = run_set13_bounds(&mol);
        entry_accum.paths14.push(Path14Configuration {
            bid1,
            bid2,
            bid3,
            kind: Path14Kind::Other,
        });
        record_path_flag(&mut entry_accum.cis_paths, path_id);
        set_15_bounds(&mol, &mut entry_mmat, &mut entry_accum, &dmat).expect("set15Bounds");

        let (mut helper_mmat, mut helper_accum) = run_set13_bounds(&mol);
        helper_accum.paths14.push(Path14Configuration {
            bid1,
            bid2,
            bid3,
            kind: Path14Kind::Other,
        });
        record_path_flag(&mut helper_accum.cis_paths, path_id);
        set_15_bounds_helper(
            &mol,
            &mut helper_mmat,
            &mut helper_accum,
            &dmat,
            nb,
            na,
            bid1,
            bid2,
            bid3,
            Path14Kind::Other,
        )
        .expect("set15BoundsHelper forward");
        set_15_bounds_helper(
            &mol,
            &mut helper_mmat,
            &mut helper_accum,
            &dmat,
            nb,
            na,
            bid3,
            bid2,
            bid1,
            Path14Kind::Other,
        )
        .expect("set15BoundsHelper reverse");

        assert_eq!(
            entry_mmat.get_lower(0, 4).unwrap(),
            helper_mmat.get_lower(0, 4).unwrap()
        );
        assert_eq!(
            entry_mmat.get_upper(0, 4).unwrap(),
            helper_mmat.get_upper(0, 4).unwrap()
        );
        assert_eq!(entry_accum.set15_atoms, helper_accum.set15_atoms);
    }
    #[test]
    fn set_15_bounds_uses_paths14_produced_by_set_14_bounds() {
        let mol = Ok::<_, ()>(fixed_topology("CCCCC", true)).expect("pentane");
        let (mmat, accum_data, _dmat) = run_set15_bounds(&mol, false, false);

        assert!(!accum_data.paths14.is_empty());
        assert!(accum_data.set15_atoms[0 * mol.atoms.len() + 4]);
        assert!(accum_data.set15_atoms[4 * mol.atoms.len() + 0]);
        assert!(mmat.get_lower(0, 4).unwrap() > 0.0);
        assert!(mmat.get_upper(0, 4).unwrap() >= mmat.get_lower(0, 4).unwrap());
    }
    #[test]
    fn set_topol_bounds_can_disable_13_and_14_stages_like_rdkit() {
        let mol = Ok::<_, ()>(fixed_topology("CCCC", true)).expect("butane");
        let disabled = run_set_topol_bounds(&mol, false, true, false, false, false, false);
        let enabled = run_set_topol_bounds(&mol, false, true, false, false, true, true);

        assert_eq!(disabled.get_upper(0, 2).unwrap(), MAX_UPPER);
        assert_eq!(disabled.get_upper(0, 3).unwrap(), MAX_UPPER);
        assert!(enabled.get_upper(0, 2).unwrap() < MAX_UPPER);
        assert!(enabled.get_upper(0, 3).unwrap() < MAX_UPPER);
    }
    #[test]
    fn set_topol_bounds_can_enable_15_stage_independently_like_rdkit() {
        let mol = Ok::<_, ()>(fixed_topology("C/C=C/CC", true)).expect("stereo pentene");
        let without_15 = run_set_topol_bounds(&mol, false, true, false, false, true, true);
        let with_15 = run_set_topol_bounds(&mol, true, true, false, false, true, true);

        assert_eq!(without_15.get_upper(0, 4).unwrap(), MAX_UPPER);
        assert!(with_15.get_upper(0, 4).unwrap() < MAX_UPPER);
        assert!(with_15.get_lower(0, 4).unwrap() > without_15.get_lower(0, 4).unwrap());
    }
    #[test]
    fn set_topol_bounds_ignores_scale_vdw_flag_like_current_rdkit_source() {
        let mol = Ok::<_, ()>(fixed_topology("CCCCCC", true)).expect("hexane");
        let scaled = run_set_topol_bounds(&mol, false, true, false, false, false, false);
        let unscaled = run_set_topol_bounds(&mol, false, false, false, false, false, false);

        assert_eq!(
            scaled.get_lower(0, 5).unwrap(),
            unscaled.get_lower(0, 5).unwrap()
        );
        assert_eq!(
            scaled.get_upper(0, 5).unwrap(),
            unscaled.get_upper(0, 5).unwrap()
        );
    }
    #[test]
    fn set_topol_bounds_with_outputs_matches_first_overload_matrix() {
        let mol = Ok::<_, ()>(fixed_topology("C/C=C/CC", true)).expect("stereo pentene");
        let plain = run_set_topol_bounds(&mol, true, true, false, false, true, true);
        let (with_outputs, bonds, angles) =
            run_set_topol_bounds_with_outputs(&mol, true, true, false, false, true, true);

        assert_eq!(matrix_rows(&plain), matrix_rows(&with_outputs));
        assert_eq!(bonds.len(), mol.bonds.len());
        assert!(!angles.is_empty());
    }
    #[test]
    fn set_topol_bounds_with_outputs_emits_exact_rdkit_bonds_and_angles_for_triple_bond_case() {
        let mol = Ok::<_, ()>(fixed_topology("CC#N", true)).expect("acetonitrile");
        let (_mmat, bonds, angles) =
            run_set_topol_bounds_with_outputs(&mol, false, true, false, false, false, false);

        assert_eq!(bonds, vec![(0, 1), (1, 2)]);
        assert_eq!(angles, vec![vec![0, 1, 2, 1]]);
    }

    fn dg_bounds_matrix_with_options(
        mol: &TopologyBlock,
        set15bounds: bool,
        scale_vdw: bool,
        smoothing: bool,
        macrocycle: bool,
    ) -> Result<Vec<Vec<f64>>, GraphBoundsError> {
        if mol.atoms.is_empty() {
            return Err(GraphBoundsError::Input("molecule has no atoms"));
        }
        let (rings, valence, hybridizations, conjugated) = fixture_chemistry(mol);
        build_bounds_matrix(
            mol,
            &rings,
            &valence,
            &hybridizations,
            &conjugated,
            set15bounds,
            scale_vdw,
            smoothing,
            macrocycle,
        )
        .map(|bounds| matrix_rows(&bounds))
    }
    fn dg_bounds_matrix(mol: &TopologyBlock) -> Result<Vec<Vec<f64>>, GraphBoundsError> {
        dg_bounds_matrix_with_options(mol, true, false, true, false)
    }
    fn wrapper_initialized_bounds(n: usize) -> BoundsMatrix {
        let mut m = BoundsMatrix::new(n).unwrap();
        init_bounds_mat(&mut m, 0.0, 1000.0).unwrap();
        m
    }
    #[test]
    fn test_single_atom() {
        let mol = Ok::<_, ()>(fixed_topology("C", true)).expect("methane skeleton");
        let result = dg_bounds_matrix(&mol).expect("dg_bounds");
        assert_eq!(result.len(), 1);
        assert_eq!(result[0][0], 0.0);
    }
    #[test]
    fn test_diatomic() {
        let mol = Ok::<_, ()>(fixed_topology("CC", true)).expect("ethane skeleton");
        let result = dg_bounds_matrix(&mol).expect("dg_bounds");
        assert_eq!(result.len(), 2);
        assert!(result[0][1] > 0.0);
        assert!(result[0][1] < 5.0);
    }
    #[test]
    fn test_ethane() {
        let mol = cosmolkit_core::add_hydrogens_with_params(
            fixed_topology("CC", true),
            cosmolkit_model::CoordinateBlock::default(),
            cosmolkit_model::MoleculeProperties::default(),
            &cosmolkit_core::AddHsParams::default(),
        )
        .expect("explicit hydrogens")
        .topology;
        let result = dg_bounds_matrix(&mol).expect("dg_bounds");
        assert_eq!(result.len(), 8);
        let upper = |i: usize, j: usize| {
            if i < j { result[i][j] } else { result[j][i] }
        };
        let carbon_pair = mol
            .bonds
            .iter()
            .find_map(|bond| {
                let begin = bond.begin().index();
                let end = bond.end().index();
                (mol.atoms[begin].atomic_number() == 6 && mol.atoms[end].atomic_number() == 6)
                    .then_some((begin, end))
            })
            .expect("ethane must contain a carbon-carbon bond");
        let carbon_hydrogen_pair = mol
            .bonds
            .iter()
            .find_map(|bond| {
                let begin = bond.begin().index();
                let end = bond.end().index();
                if mol.atoms[begin].atomic_number() == 6 && mol.atoms[end].atomic_number() == 1 {
                    Some((begin, end))
                } else if mol.atoms[begin].atomic_number() == 1
                    && mol.atoms[end].atomic_number() == 6
                {
                    Some((end, begin))
                } else {
                    None
                }
            })
            .expect("ethane must contain a carbon-hydrogen bond");

        assert!(
            upper(carbon_pair.0, carbon_pair.1) > 1.0 && upper(carbon_pair.0, carbon_pair.1) < 3.0
        );
        assert!(
            upper(carbon_hydrogen_pair.0, carbon_hydrogen_pair.1) > 1.0
                && upper(carbon_hydrogen_pair.0, carbon_hydrogen_pair.1) < 3.0
        );
    }
    #[test]
    fn test_empty() {
        let mol = TopologyBlock::try_from_parts(vec![], vec![], vec![], vec![]).expect("build");
        let err =
            dg_bounds_matrix(&mol).expect_err("empty molecule must now follow setTopolBounds");
        assert!(matches!(
            err,
            GraphBoundsError::Input(message) if message == "molecule has no atoms"
        ));
    }
    #[test]
    fn dg_bounds_matrix_returns_error_instead_of_panicking_for_3rj7_re_complex() {
        let mol2 = include_str!("../../../testdata/mol2/fixtures/3rj7_ligand.mol2");
        let record = cosmolkit_io::read_mol2_detached(mol2)
            .expect("3rj7 mol2 should parse")
            .expect("3rj7 mol2 should contain a molecule");

        // RDKit❗✔️:       constexpr auto sanitizeFlags = MolOps::SanitizeFlags::SANITIZE_ALL ^
        // RDKit❗✔️:                             MolOps::SanitizeFlags::SANITIZE_CLEANUP_ORGANOMETALLICS;
        // RDKit❗✔️:         MolOps::sanitizeMol(*res, failedOp, sanitizeFlags);
        // The detached IO reader leaves the high-level chemistry preparation to
        // its caller. Reuse the source-selected existing sanitizer here.
        let operations = cosmolkit_core::SanitizeOperations::CLEANUP
            | cosmolkit_core::SanitizeOperations::PROPERTIES
            | cosmolkit_core::SanitizeOperations::SYMM_RINGS
            | cosmolkit_core::SanitizeOperations::KEKULIZE
            | cosmolkit_core::SanitizeOperations::FIND_RADICALS
            | cosmolkit_core::SanitizeOperations::SET_AROMATICITY
            | cosmolkit_core::SanitizeOperations::SET_CONJUGATION
            | cosmolkit_core::SanitizeOperations::SET_HYBRIDIZATION
            | cosmolkit_core::SanitizeOperations::CLEANUP_CHIRALITY
            | cosmolkit_core::SanitizeOperations::ADJUST_HS
            | cosmolkit_core::SanitizeOperations::CLEANUP_ATROPISOMERS;
        let prepared = sanitize_topology(&record.topology, &SanitizeParams { operations })
            .expect("original MOL2 topological preparation")
            .topology;
        let err = dg_bounds_matrix(&prepared).expect_err("invalid bounds must be reported");

        assert!(matches!(
            err,
            GraphBoundsError::InvalidBounds(message)
                if message.contains("bad lower bound") && message.contains("atom pair")
        ));
    }
    #[test]
    fn dg_bounds_matrix_matches_source_backed_set_topol_bounds_path() {
        let mol = Ok::<_, ()>(fixed_topology("C/C=C/CC", true)).expect("stereo pentene");
        let (rings, valence, hybridizations, conjugated) = fixture_chemistry(&mol);
        let result = dg_bounds_matrix(&mol).expect("dg_bounds");
        let mut mmat = wrapper_initialized_bounds(mol.atoms.len());
        set_topol_bounds(
            &mol,
            &rings,
            &valence,
            &hybridizations,
            &conjugated,
            &mut mmat,
            true,
            false,
            false,
            false,
            true,
            true,
        )
        .expect("setTopolBounds");
        assert!(crate::smoothing::triangle_smooth_bounds_shared(
            &mut mmat, 0.0
        ));

        assert_eq!(result, matrix_rows(&mmat));
    }
    #[test]
    fn dg_bounds_matrix_uses_rdkit_wrapper_defaults() {
        let mol = Ok::<_, ()>(fixed_topology("CCCCCC", true)).expect("hexane");
        let from_default = dg_bounds_matrix(&mol).expect("default dg_bounds");
        let explicit_default =
            dg_bounds_matrix_with_options(&mol, true, false, true, false).expect("explicit");
        let scaled =
            dg_bounds_matrix_with_options(&mol, true, true, true, false).expect("scaled vdw");

        assert_eq!(from_default, explicit_default);
        assert_eq!(from_default, scaled);
    }
    #[test]
    fn dg_bounds_matrix_with_options_can_skip_triangle_smoothing() {
        let mol = Ok::<_, ()>(fixed_topology("C/C=C/CC", true)).expect("stereo pentene");
        let (rings, valence, hybridizations, conjugated) = fixture_chemistry(&mol);
        let unsmoothed =
            dg_bounds_matrix_with_options(&mol, true, false, false, false).expect("unsmoothed");
        let smoothed =
            dg_bounds_matrix_with_options(&mol, true, false, true, false).expect("smoothed");
        let mut manual_unsmoothed = wrapper_initialized_bounds(mol.atoms.len());
        set_topol_bounds(
            &mol,
            &rings,
            &valence,
            &hybridizations,
            &conjugated,
            &mut manual_unsmoothed,
            true,
            false,
            false,
            false,
            true,
            true,
        )
        .expect("setTopolBounds");
        let mut manual_smoothed = manual_unsmoothed.clone();
        assert!(crate::smoothing::triangle_smooth_bounds_shared(
            &mut manual_smoothed,
            0.0
        ));

        assert_eq!(unsmoothed, matrix_rows(&manual_unsmoothed));
        assert_eq!(smoothed, matrix_rows(&manual_smoothed));
    }
    #[test]
    fn dg_bounds_matrix_with_options_forwards_macrocycle14config_like_rdkit_wrapper() {
        let mol = Ok::<_, ()>(fixed_topology("C1CCCCCCCCC1", true)).expect("cyclodecane");
        let (rings, valence, hybridizations, conjugated) = fixture_chemistry(&mol);
        let wrapper_without_macrocycle =
            dg_bounds_matrix_with_options(&mol, true, false, false, false).expect("wrapper plain");
        let wrapper_with_macrocycle = dg_bounds_matrix_with_options(&mol, true, false, false, true)
            .expect("wrapper macrocycle");

        let mut manual_without_macrocycle = wrapper_initialized_bounds(mol.atoms.len());
        set_topol_bounds(
            &mol,
            &rings,
            &valence,
            &hybridizations,
            &conjugated,
            &mut manual_without_macrocycle,
            true,
            false,
            false,
            false,
            true,
            true,
        )
        .expect("manual plain");

        let mut manual_with_macrocycle = wrapper_initialized_bounds(mol.atoms.len());
        set_topol_bounds(
            &mol,
            &rings,
            &valence,
            &hybridizations,
            &conjugated,
            &mut manual_with_macrocycle,
            true,
            false,
            true,
            false,
            true,
            true,
        )
        .expect("manual macrocycle");

        assert_eq!(
            wrapper_without_macrocycle,
            matrix_rows(&manual_without_macrocycle)
        );
        assert_eq!(
            wrapper_with_macrocycle,
            matrix_rows(&manual_with_macrocycle)
        );
    }

    #[test]
    fn etv2_fixture_bounds_first_boundary_diagnostic() {
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../testdata/conformer/fixtures/rdkit/test_data/torsion.etkdg.v2.mol");
        let (topology, _, _) =
            cosmolkit_io::sdf::read_v2000_detached(&std::fs::read_to_string(path).unwrap())
                .unwrap();
        let topology = cosmolkit_core::sanitize_topology(&topology, &Default::default())
            .unwrap()
            .topology;
        let actual = dg_bounds_matrix(&topology).unwrap();
        println!("etv2_bounds {:?}", actual);
        let expected: Vec<Vec<f64>> = serde_json::from_str(r#"[[0.0, 1.3550727841738393, 2.3994001949364825, 2.8015808378091798, 2.3994001949364825, 3.7778080387314956, 1.3550727841738393, 2.407040464768869, 3.6197111687331383, 2.1455336014957314, 3.397208259243587, 3.892998878135974, 4.56985027850704, 4.56985027850704, 4.56985027850704, 4.432194982292257, 4.432194982292257], [1.3350727841738392, 0.0, 1.389256405400789, 2.4289421708189836, 2.8185128108015784, 4.3145128108015784, 2.3697344020672166, 3.6593534435886563, 4.802764557510261, 1.0914180403267946, 2.1762036346800016, 3.42914944022548, 4.989140369039992, 4.989140369039992, 4.989140369039992, 5.629663703020149, 5.629663703020149], [2.3194001949364824, 1.369256405400789, 0.0, 1.389256405400789, 2.4289421708189836, 3.810195428346906, 2.8185128108015784, 4.216550359588743, 5.251547219002712, 2.1762036346800016, 1.0914180403267946, 2.1762036346800016, 4.59950692012806, 4.59950692012806, 4.59950692012806, 6.07844211175451, 6.07844211175451], [2.6815808378091797, 2.3489421708189835, 1.369256405400789, 0.0, 1.389256405400789, 2.5219587527323406, 2.4289421708189836, 3.7174679728830053, 4.861972272176078, 3.429149440225479, 2.1762036346800016, 1.0914180403267946, 3.3912458481603385, 3.3912458481603385, 3.3912458481603385, 5.6888714717719155, 5.6888714717719155], [2.3194001949364824, 2.6473992346926964, 2.3489421708189835, 1.369256405400789, 0.0, 1.496, 1.389256405400789, 2.4365508859157763, 3.651972415318165, 3.9099308511283732, 3.429149440225479, 2.1762036346800016, 2.1706275582384142, 2.1706275582384142, 2.1706275582384142, 4.461705403439164, 4.461705403439164], [3.6578080387314955, 4.113341878059448, 3.690195428346906, 2.4419587527323405, 1.476, 0.0, 2.5219587527323406, 2.8775527270508237, 4.262743343028898, 5.405930851128373, 4.696825339318566, 2.745916934553264, 1.119400794877744, 1.119400794877744, 1.119400794877744, 4.902707244574211, 4.902707244574211], [1.3350727841738392, 2.2897344020672166, 2.647399234692697, 2.3489421708189835, 1.369256405400789, 2.4419587527323405, 0.0, 1.3980375487871646, 2.433034408201134, 3.370722773835688, 3.9099308511283732, 3.429149440225479, 3.3912458481603385, 3.3912458481603385, 3.3912458481603385, 3.259929300952932, 3.259929300952932], [2.327040464768869, 3.5393534435886562, 4.0153819264891695, 3.5974679728830052, 2.356550885915776, 2.7575527270508235, 1.3780375487871646, 0.0, 1.385190615978074, 4.550974944998207, 5.307968399915538, 4.6101483500685365, 3.9888603940752563, 3.9888603940752563, 3.9888603940752563, 2.025154517523388, 2.025154517523388], [2.748293021721179, 4.013325134150105, 2.8049999999999997, 4.064518960585944, 2.7652633012138907, 2.318908467466266, 2.353034408201134, 1.365190615978074, 0.0, 5.76524477022887, 6.342965259329507, 5.8281760499981665, 5.37405101005333, 5.37405101005333, 5.37405101005333, 1.0544193818795018, 1.0544193818795018], [2.0655336014957313, 1.0714180403267946, 2.0962036346800015, 3.309149440225479, 3.7087710772848337, 3.021923837732654, 3.250722773835688, 4.390974944998207, 3.0057843290201327, 0.0, 2.520674445727583, 4.342013160886561, 6.0805584093667875, 6.0805584093667875, 6.0805584093667875, 6.576129462521594, 6.576129462521594], [3.277208259243587, 2.0962036346800015, 1.0714180403267946, 2.0962036346800015, 3.309149440225479, 4.536825339318565, 3.7087710772848337, 2.923963886162375, 2.8, 2.400674445727583, 0.0, 2.520674445727583, 5.56744948284034, 5.56744948284034, 5.56744948284034, 7.169860152081306, 7.169860152081306], [3.7429536643353405, 3.3091494402254797, 2.0962036346800015, 1.0714180403267946, 2.0962036346800015, 2.625916934553264, 3.309149440225479, 4.4501483500685355, 3.0649577340904615, 4.18201316088656, 2.400674445727583, 0.0, 3.8367137758603302, 3.8367137758603302, 3.8367137758603302, 6.635302867591925, 6.635302867591925], [3.813867062490825, 2.9939410831817046, 3.8480259645306063, 2.550292841931279, 2.090627558238414, 1.099400794877744, 2.5502928419312787, 2.2122192256860433, 2.38, 2.4, 3.417424544440821, 2.200076991737953, 0.0, 1.8519655808531272, 1.8519655808531272, 6.014014911598644, 6.014014911598644], [3.813867062490825, 2.9939410831817046, 3.8480259645306063, 2.550292841931279, 2.090627558238414, 1.099400794877744, 2.5502928419312787, 2.2122192256860433, 2.38, 2.4, 3.417424544440821, 2.200076991737953, 1.7719655808531272, 0.0, 1.8519655808531272, 6.014014911598644, 6.014014911598644], [3.813867062490825, 2.9939410831817046, 3.8480259645306063, 2.550292841931279, 2.090627558238414, 1.099400794877744, 2.5502928419312787, 2.2122192256860433, 2.38, 2.4, 3.417424544440821, 2.200076991737953, 1.7719655808531272, 1.7719655808531272, 0.0, 6.014014911598644, 6.014014911598644], [1.9599999999999997, 2.9589057522706037, 2.9, 3.010099578706442, 2.03, 2.465, 2.427612389561217, 1.9451545175233877, 1.0344193818795018, 2.4, 2.4, 2.4249938325451477, 2.4, 2.4, 2.4, 0.0, 1.7458325365181517], [1.9599999999999997, 2.9589057522706037, 2.9, 3.010099578706442, 2.03, 2.465, 2.427612389561217, 1.9451545175233877, 1.0344193818795018, 2.4, 2.4, 2.4249938325451477, 2.4, 2.4, 2.4, 1.6658325365181517, 0.0]]"#).unwrap();
        let mismatches: Vec<_> = actual
            .iter()
            .enumerate()
            .flat_map(|(i, row)| row.iter().enumerate().map(move |(j, v)| (i, j, *v)))
            .filter(|&(i, j, v)| (v - expected[i][j]).abs() > 1e-10)
            .collect();
        println!(
            "etv2_bounds_mismatches {:?}",
            mismatches
                .iter()
                .map(|&(i, j, v)| (i, j, v, expected[i][j]))
                .collect::<Vec<_>>()
        );
        assert!(
            mismatches.is_empty(),
            "first divergent bounds boundary: {:?}",
            mismatches
        );
    }
}
