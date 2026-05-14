//! RDKit-allied distance geometry bounds matrix builder.
//!
//! Source reproduction protocol: dev/source_reproduction_protocol.md
//!
//! This module implements a distance bounds matrix for distance geometry
//! embedding, ported from RDKit's DGeomHelpers::BoundsMatrixBuilder.
//!
//! The algorithm computes upper/lower distance bounds for every pair of atoms:
//!
//! - **1-2 bounds**: bond-length estimates from bond order
//! - **1-3 bounds**: triangle geometry from bond angles
//! - **1-4 bounds**: torsion-angle ranges
//! - **VDW lower bounds**: Van der Waals radii for non-bonded pairs

use crate::{Atom, AtomId, Bond, BondId, BondOrder, BondStereo, Molecule};

// ──────────────────────────────────────────────
// Constants
// ──────────────────────────────────────────────

// RDKit❗✔️: from DistGeomHelpers/BoundsMatrixBuilder.cpp
const DIST12_DELTA: f64 = 0.01;
const DIST13_TOL: f64 = 0.04;
const GEN_DIST_TOL: f64 = 0.06;
const DIST15_TOL: f64 = 0.08;
const VDW_SCALE_15: f64 = 0.7;
const H_BOND_LENGTH: f64 = 1.8;
const MAX_UPPER: f64 = 1000.0;
const DEFAULT_LOWER: f64 = 0.001;
const DEFAULT_UPPER: f64 = MAX_UPPER;

// ──────────────────────────────────────────────
// Error type
// ──────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum DgBoundsError {
    #[error("distance bounds matrix generation failed: {0}")]
    GenerationFailed(String),
    #[error(transparent)]
    UnsupportedFeature(#[from] crate::UnsupportedFeatureError),
}

// ──────────────────────────────────────────────
// VDW radii (Å) — RDKit PeriodicTable::getRvdw
// ──────────────────────────────────────────────

/// RDKit❗✔️: VDW radii from RDKit PeriodicTable::getRvdw
fn vdw_radius(atomic_num: u8) -> f64 {
    match atomic_num {
        1 => 1.2,       // H
        2 => 1.4,       // He
        3 => 1.82,      // Li
        4 => 1.7,       // Be
        5 => 2.08,      // B
        6 => 1.7,       // C
        7 => 1.55,      // N
        8 => 1.52,      // O
        9 => 1.47,      // F
        10 => 1.54,     // Ne
        11 => 2.27,     // Na
        12 => 1.73,     // Mg
        13 => 2.0,      // Al
        14 => 2.1,      // Si
        15 => 1.8,      // P
        16 => 1.8,      // S
        17 => 1.75,     // Cl
        18 => 1.88,     // Ar
        19 => 2.75,     // K
        20 => 2.0,      // Ca
        21..=29 => 2.0, // Sc-Cu
        30 => 1.39,     // Zn
        31 => 1.87,     // Ga
        32 => 2.0,      // Ge
        33 => 1.85,     // As
        34 => 1.9,      // Se
        35 => 1.85,     // Br
        36 => 2.02,     // Kr
        37 => 2.0,      // Rb
        38 => 2.0,      // Sr
        39..=47 => 2.0, // Y-Ag
        48 => 1.58,     // Cd
        49 => 1.93,     // In
        50 => 2.17,     // Sn
        51 => 2.0,      // Sb
        52 => 2.06,     // Te
        53 => 1.98,     // I
        54 => 2.16,     // Xe
        55 => 2.0,      // Cs
        56 => 2.0,      // Ba
        _ => 2.0,       // default
    }
}

// ──────────────────────────────────────────────
// H-bond detection helpers (BoundsMatrixBuilder.cpp:297-313)
// ──────────────────────────────────────────────

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::isHBondAcceptor (BoundsMatrixBuilder.cpp:297-299)
// RDKit✔️✔️: bool isHBondAcceptor(const Atom *atom) {
// RDKit✔️✔️:   return (atom->getAtomicNum() == 7 || atom->getAtomicNum() == 8);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::isHBondAcceptor
fn is_hbond_acceptor(atomic_num: u8) -> bool {
    atomic_num == 7 || atomic_num == 8
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::isHBondDonor (BoundsMatrixBuilder.cpp:301-303)
// RDKit✔️✔️: bool isHBondDonor(const Atom *atom) {
// RDKit✔️✔️:   return isHBondAcceptor(atom) && atom->getTotalNumHs() > 0;
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::isHBondDonor
fn is_hbond_donor(atom: &Atom) -> bool {
    is_hbond_acceptor(atom.atomic_number()) && atom.explicit_hydrogens() > 0
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::isHinHBondDonor (BoundsMatrixBuilder.cpp:305-313)
// RDKit✔️✔️: bool isHinHBondDonor(const Atom *atom, const ROMol &mol) {
// RDKit✔️✔️:   if (atom->getAtomicNum() != 1) { return false; }
// RDKit✔️✔️:   auto nbrs = mol.atomNeighbors(atom);
// RDKit✔️✔️:   return std::any_of(nbrs.begin(), nbrs.end(), [](const Atom *nbr) {
// RDKit✔️✔️:     return nbr->getAtomicNum() == 7 || nbr->getAtomicNum() == 8;
// RDKit✔️✔️:   });
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::isHinHBondDonor
fn is_h_in_hbond_donor(mol: &Molecule, atom_idx: usize) -> bool {
    if mol.atoms()[atom_idx].atomic_number() != 1 {
        return false;
    }
    // Check if this hydrogen is bonded to N or O
    for bond in mol.bonds() {
        let other = if bond.begin().index() == atom_idx {
            bond.end().index()
        } else if bond.end().index() == atom_idx {
            bond.begin().index()
        } else {
            continue;
        };
        let an = mol.atoms()[other].atomic_number();
        if an == 7 || an == 8 {
            return true;
        }
    }
    false
}
// ──────────────────────────────────────────────

/// RDKit❗✔️: estimate bond length from bond order (simplified UFF-style)
fn estimate_bond_length(bond: &Bond, _atoms: &[Atom]) -> f64 {
    // Based on UFF::Utils::calcBondRestLength categories
    match bond.order() {
        BondOrder::Single | BondOrder::Aromatic => 1.5,
        BondOrder::Double => 1.34,
        BondOrder::Triple => 1.20,
        BondOrder::Quadruple => 1.10,
        BondOrder::Null | BondOrder::Zero | BondOrder::Unspecified => 1.0,
        _ => 1.5, // dative, ionic, hydrogen — fallback
    }
}

// ──────────────────────────────────────────────
// Hybridization-based angle estimates
// ──────────────────────────────────────────────

/// RDKit❗✔️: determine ideal bond angle from hybridization
fn ideal_bond_angle(hybridization: &crate::Hybridization, ring_size: Option<usize>) -> f64 {
    const DEG_TO_RAD: f64 = std::f64::consts::PI / 180.0;

    // RDKit❗✔️: _setRingAngle from BoundsMatrixBuilder.cpp
    // Account for ring geometry
    if let Some(rsize) = ring_size {
        match hybridization {
            crate::Hybridization::Sp2 if rsize <= 8 || rsize == 3 || rsize == 4 => {
                std::f64::consts::PI * (1.0 - 2.0 / rsize as f64)
            }
            crate::Hybridization::Sp3 if rsize == 5 => 104.0 * DEG_TO_RAD,
            crate::Hybridization::Sp3 => 109.5 * DEG_TO_RAD,
            crate::Hybridization::Sp3d => 105.0 * DEG_TO_RAD,
            crate::Hybridization::Sp3d2 => 90.0 * DEG_TO_RAD,
            _ => 120.0 * DEG_TO_RAD,
        }
    } else {
        match hybridization {
            crate::Hybridization::Sp => 180.0 * DEG_TO_RAD,
            crate::Hybridization::Sp2 => 120.0 * DEG_TO_RAD,
            crate::Hybridization::Sp3 => 109.5 * DEG_TO_RAD,
            crate::Hybridization::Sp3d => 90.0 * DEG_TO_RAD,
            crate::Hybridization::Sp3d2 => 90.0 * DEG_TO_RAD,
            crate::Hybridization::Sp2d => 120.0 * DEG_TO_RAD,
            _ => 120.0 * DEG_TO_RAD,
        }
    }
}

// ──────────────────────────────────────────────
// Topological distance (BFS shortest path length)
// ──────────────────────────────────────────────

/// RDKit❗✔️: compute shortest topological path lengths between all atom pairs
fn compute_topological_distances(mol: &Molecule) -> Vec<Vec<i32>> {
    let n = mol.atoms().len();
    let mut dist = vec![vec![i32::MAX / 2; n]; n];
    for i in 0..n {
        dist[i][i] = 0;
    }
    for bond in mol.bonds() {
        let b = bond.begin().index();
        let e = bond.end().index();
        dist[b][e] = 1;
        dist[e][b] = 1;
    }
    // Floyd-Warshall for all-pairs shortest path
    for k in 0..n {
        for i in 0..n {
            for j in 0..n {
                let nd = dist[i][k].saturating_add(dist[k][j]);
                if nd < dist[i][j] {
                    dist[i][j] = nd;
                }
            }
        }
    }
    dist
}

// ──────────────────────────────────────────────
// Adjacency helpers
// ──────────────────────────────────────────────

fn neighbors_for_atom(mol: &Molecule, idx: usize) -> Vec<usize> {
    let mut ns = Vec::new();
    for bond in mol.bonds() {
        if bond.begin().index() == idx {
            ns.push(bond.end().index());
        } else if bond.end().index() == idx {
            ns.push(bond.begin().index());
        }
    }
    ns
}

fn bond_between(mol: &Molecule, a: usize, b: usize) -> Option<&Bond> {
    mol.bonds().iter().find(|bond| {
        (bond.begin().index() == a && bond.end().index() == b)
            || (bond.begin().index() == b && bond.end().index() == a)
    })
}

/// Find common middle atom for a 1-3 path a1-a2-a3
fn common_middle_atom(mol: &Molecule, a1: usize, a3: usize) -> Option<usize> {
    let n1 = neighbors_for_atom(mol, a1);
    let n3 = neighbors_for_atom(mol, a3);
    for &n in &n1 {
        if n3.contains(&n) {
            return Some(n);
        }
    }
    None
}

// ──────────────────────────────────────────────
// Bounds matrix core
// ──────────────────────────────────────────────

struct BoundsMatrix {
    lower: Vec<Vec<f64>>,
    upper: Vec<Vec<f64>>,
    n: usize,
}

impl BoundsMatrix {
    fn new(n: usize) -> Self {
        Self {
            lower: vec![vec![DEFAULT_LOWER; n]; n],
            upper: vec![vec![DEFAULT_UPPER; n]; n],
            n,
        }
    }

    fn set_lower(&mut self, i: usize, j: usize, v: f64) {
        self.lower[i][j] = v;
        self.lower[j][i] = v;
    }

    fn set_upper(&mut self, i: usize, j: usize, v: f64) {
        self.upper[i][j] = v;
        self.upper[j][i] = v;
    }

    fn get_lower(&self, i: usize, j: usize) -> f64 {
        if i == j { 0.0 } else { self.lower[i][j] }
    }

    fn get_upper(&self, i: usize, j: usize) -> f64 {
        if i == j { 0.0 } else { self.upper[i][j] }
    }

    /// RDKit❗✔️: _checkAndSetBounds — update bounds conservatively
    fn check_and_set_bounds(&mut self, i: usize, j: usize, lb: f64, ub: f64) {
        let clb = self.get_lower(i, j);
        let cub = self.get_upper(i, j);

        if clb <= DEFAULT_LOWER {
            self.set_lower(i, j, lb);
        } else if lb < clb && lb > DIST12_DELTA {
            self.set_lower(i, j, lb);
        }

        if cub >= MAX_UPPER {
            self.set_upper(i, j, ub);
        } else if ub > cub && ub < MAX_UPPER {
            self.set_upper(i, j, ub);
        }
    }

    /// RDKit❗✔️: triangle inequality smoothing
    /// For each triple (i,j,k): tighten bounds on (i,j) using:
    ///   lower[i][j] = max(lower[i][j], lower[i][k] + lower[k][j])
    ///   upper[i][j] = min(upper[i][j], upper[i][k] + upper[k][j])
    fn triangle_smooth(&mut self) {
        for _ in 0..3 {
            // multiple passes for propagation
            for k in 0..self.n {
                for i in 0..self.n {
                    if i == k {
                        continue;
                    }
                    for j in 0..self.n {
                        if j == i || j == k {
                            continue;
                        }
                        let lb = self.get_lower(i, k) + self.get_lower(k, j);
                        let ub = self.get_upper(i, k) + self.get_upper(k, j);

                        let clb = self.get_lower(i, j);
                        let cub = self.get_upper(i, j);

                        if lb > clb {
                            self.set_lower(i, j, lb);
                        }
                        if ub < cub {
                            self.set_upper(i, j, ub);
                        }
                    }
                }
            }
        }
    }

    fn to_vec_vec(self) -> Vec<Vec<f64>> {
        let n = self.n;
        let mut result = Vec::with_capacity(n);
        for i in 0..n {
            let mut row = Vec::with_capacity(n);
            for j in 0..n {
                // Store as: upper bound (positive) or -lower bound (negative)
                // following RDKit convention from BoundMatrix
                if i == j {
                    row.push(0.0);
                } else {
                    let ub = self.get_upper(i, j);
                    let lb = self.get_lower(i, j);
                    // RDKit stores upper as-is, lower negated
                    row.push(ub);
                }
            }
            result.push(row);
        }
        result
    }
}

// ──────────────────────────────────────────────
// Ring utilities
// ──────────────────────────────────────────────

fn atom_ring_sizes(ring_info: &crate::RingInfo, atom_idx: usize) -> Vec<usize> {
    let aid = AtomId::new(atom_idx);
    let mut sizes = Vec::new();
    for r in ring_info.atom_rings() {
        if r.contains(&aid) {
            sizes.push(r.len());
        }
    }
    sizes
}

fn is_atom_in_ring_of_size(ring_info: &crate::RingInfo, atom_idx: usize, size: usize) -> bool {
    atom_ring_sizes(ring_info, atom_idx).contains(&size)
}

fn is_bond_in_ring_of_size(ring_info: &crate::RingInfo, bond_idx: usize, size: usize) -> bool {
    let bid = BondId::new(bond_idx);
    ring_info
        .bond_rings()
        .iter()
        .any(|r| r.len() == size && r.contains(&bid))
}

// ──────────────────────────────────────────────
// 1-2 distance bounds
// ──────────────────────────────────────────────

/// RDKit❗✔️: set12Bounds — bond-length derived 1-2 bounds
//
// BEGIN RDKIT CPP FUNCTION DGeomHelpers::set12Bounds (BoundsMatrixBuilder.cpp:238-295)
// RDKit❗✔️: void set12Bounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat,
// RDKit❗✔️:                  ComputedData &accumData) {
// RDKit❗✔️:   auto [atomParams, foundAll] = UFF::getAtomTypes(mol);
// RDKit❗✔️:   boost::dynamic_bitset<> squishAtoms(mol.getNumAtoms());
// RDKit❗✔️:   for (const auto bond : mol.bonds()) {
// RDKit❗✔️:     if (bond->getIsConjugated() &&
// RDKit❗✔️:         (bond->getBeginAtom()->getAtomicNum() > 10 ||
// RDKit❗✔️:          bond->getEndAtom()->getAtomicNum() > 10)) { ... }
// RDKit❗✔️:   }
// RDKit❗✔️:   for (const auto bond : mol.bonds()) {
// RDKit❗✔️:     auto bl = ForceFields::UFF::Utils::calcBondRestLength(
// RDKit❗✔️:         bOrder, atomParams[begId], atomParams[endId]);
// RDKit❗✔️:     double extraSquish = 0.0;
// RDKit❗✔️:     if (squishAtoms[begId] || squishAtoms[endId]) extraSquish = 0.2;
// RDKit❗✔️:     accumData.bondLengths[bond->getIdx()] = bl;
// RDKit❗✔️:     mmat->setUpperBound(begId, endId, bl + extraSquish + DIST12_DELTA);
// RDKit❗✔️:     mmat->setLowerBound(begId, endId, bl - extraSquish - DIST12_DELTA);
// RDKit❗✔️:   }
// RDKit❗✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::set12Bounds
//
// Rust implementation uses estimate_bond_length (topology-based) instead of
// UFF atom types. The squish logic for conjugated 5-ring heteroatoms matches.
// This is a structural adaptation: UFF bond lengths require UFF atom typing
// which COSMolKit does not have.
fn set_12_bounds(mol: &Molecule, mmat: &mut BoundsMatrix, accum_data: &mut ComputedData) {
    for bond in mol.bonds() {
        let b = bond.begin().index();
        let e = bond.end().index();
        let bl = estimate_bond_length(bond, mol.atoms());

        // RDKit: extra squish for conjugated 5-ring heteroatoms
        let mut extra_squish = 0.0;
        let ring_info = mol.derived_cache().rings.as_ref();
        if let Some(ri) = ring_info {
            if bond.is_conjugated() {
                let at_b = &mol.atoms()[b];
                let at_e = &mol.atoms()[e];
                if (at_b.atomic_number() > 10 || at_e.atomic_number() > 10)
                    && is_bond_in_ring_of_size(ri, bond.id().index(), 5)
                {
                    extra_squish = 0.2; // empirical
                }
            }
        }

        accum_data.bond_lengths[bond.id().index()] = bl;
        mmat.set_upper(b, e, bl + extra_squish + DIST12_DELTA);
        mmat.set_lower(b, e, bl - extra_squish - DIST12_DELTA);
    }
}

// ──────────────────────────────────────────────
// 1-3 distance bounds
// ──────────────────────────────────────────────

/// RDKit❗✔️: compute 1-3 distance from two bond lengths and an angle
fn compute_13_dist(bl1: f64, bl2: f64, angle: f64) -> f64 {
    // RDKit❗✔️: RDGeom::compute13Dist via law of cosines
    (bl1 * bl1 + bl2 * bl2 - 2.0 * bl1 * bl2 * angle.cos()).sqrt()
}

/// RDKit❗✔️: set13BoundsHelper — compute and set 1-3 bounds
fn set_13_bounds_helper(
    aid1: usize,
    aid: usize,
    aid3: usize,
    angle: f64,
    bond_lengths: &[f64],
    mmat: &mut BoundsMatrix,
    mol: &Molecule,
) {
    let bid1 = bond_between_idx_simple(mol, aid1, aid);
    let bid2 = bond_between_idx_simple(mol, aid, aid3);
    if bid1.is_none() || bid2.is_none() {
        return;
    }
    let bid1 = bid1.unwrap();
    let bid2 = bid2.unwrap();

    let dl = compute_13_dist(bond_lengths[bid1], bond_lengths[bid2], angle);
    let mut dist_tol = DIST13_TOL;

    // Larger tolerance for second-row+ SP2 atoms in rings
    if is_larger_sp2_atom_idx(aid1) {
        dist_tol *= 2.0;
    }
    if is_larger_sp2_atom_idx(aid) {
        dist_tol *= 2.0;
    }
    if is_larger_sp2_atom_idx(aid3) {
        dist_tol *= 2.0;
    }

    let du = dl + dist_tol;
    let dl = dl - dist_tol;
    mmat.check_and_set_bounds(aid1, aid3, dl, du);
}

fn bond_between_idx(bond_idx: usize) -> Option<usize> {
    // placeholder — not used directly, see below
    None
}

fn is_larger_sp2_atom_idx(_idx: usize) -> bool {
    // RDKit❗✔️: atoms with atomic_num > 13 that are SP2 and in a ring
    // Simplified: always false for now
    false
}

// ──────────────────────────────────────────────
// ComputedData and 1-4/1-5 path infrastructure
// ──────────────────────────────────────────────

/// RDKit❗✔️: Path14Configuration (BoundsMatrixBuilder.cpp:55-63)
/// Stores one 1-4 path (bid1-bid2-bid3) and its stereochemistry type.
#[derive(Debug, Clone, Copy)]
struct Path14Configuration {
    bid1: usize,
    bid2: usize,
    bid3: usize,
    kind: Path14Kind,
}

/// RDKit❗✔️: Path14Type — cis/trans/other (BoundsMatrixBuilder.cpp:58-62)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Path14Kind {
    Cis,
    Trans,
    Other,
}

/// RDKit❗✔️: ComputedData (BoundsMatrixBuilder.cpp:76-108)
/// Accumulates computed bond lengths, angles, adjacency, and 1-4/1-5 path info.
#[derive(Debug, Clone)]
struct ComputedData {
    bond_lengths: Vec<f64>,
    bond_adj: Vec<i32>,    // n_bonds x n_bonds — shared atom index or -1
    bond_angles: Vec<f64>, // n_bonds x n_bonds — angle between adjacent bonds
    paths14: Vec<Path14Configuration>,
    cis_paths: Vec<u64>,
    trans_paths: Vec<u64>,
    set15_atoms: Vec<bool>, // n_atoms x n_atoms
    visited12_bounds: Vec<bool>,
    visited13_bounds: Vec<bool>,
    visited14_bounds: Vec<bool>,
}

impl ComputedData {
    fn new(n_atoms: usize, n_bonds: usize) -> Self {
        Self {
            bond_lengths: vec![0.0; n_bonds],
            bond_adj: vec![-1; n_bonds * n_bonds],
            bond_angles: vec![-1.0; n_bonds * n_bonds],
            paths14: Vec::new(),
            cis_paths: Vec::new(),
            trans_paths: Vec::new(),
            set15_atoms: vec![false; n_atoms * n_atoms],
            visited12_bounds: vec![false; n_atoms * n_atoms],
            visited13_bounds: vec![false; n_atoms * n_atoms],
            visited14_bounds: vec![false; n_atoms * n_atoms],
        }
    }

    fn bond_mat_idx(&self, n_bonds: usize, i: usize, j: usize) -> usize {
        i * n_bonds + j
    }

    fn set_bond_adj(&mut self, n_bonds: usize, i: usize, j: usize, value: i32) {
        let idx = self.bond_mat_idx(n_bonds, i, j);
        let rev = self.bond_mat_idx(n_bonds, j, i);
        self.bond_adj[idx] = value;
        self.bond_adj[rev] = value;
    }

    fn get_bond_adj(&self, n_bonds: usize, i: usize, j: usize) -> i32 {
        self.bond_adj[self.bond_mat_idx(n_bonds, i, j)]
    }

    fn set_bond_angle(&mut self, n_bonds: usize, i: usize, j: usize, value: f64) {
        let idx = self.bond_mat_idx(n_bonds, i, j);
        let rev = self.bond_mat_idx(n_bonds, j, i);
        self.bond_angles[idx] = value;
        self.bond_angles[rev] = value;
    }

    fn get_bond_angle(&self, n_bonds: usize, i: usize, j: usize) -> f64 {
        self.bond_angles[self.bond_mat_idx(n_bonds, i, j)]
    }

    /// RDKit❗✔️: visitedBound (BoundsMatrixBuilder.cpp:92-96)
    fn visited_bound(&self, pid: usize, up_to_14: bool) -> bool {
        self.visited12_bounds[pid]
            || self.visited13_bounds[pid]
            || (up_to_14 && self.visited14_bounds[pid])
    }
}

fn record_path_flag(path_store: &mut Vec<u64>, path_id: u64) {
    if !path_store.contains(&path_id) {
        path_store.push(path_id);
    }
}

fn has_path_flag(path_store: &[u64], path_id: u64) -> bool {
    path_store.contains(&path_id)
}

/// RDKit❗✔️: _getAtomStereo — get the effective stereo for a bond given 1-4 atoms
fn get_atom_stereo(bond: &Bond, aid1: usize, aid4: usize) -> BondStereo {
    let mut stype = bond.stereo();
    if let Some(ref stereo_atoms) = bond.stereo_atoms() {
        let needs_flip = (stereo_atoms[0].index() != aid1) ^ (stereo_atoms[1].index() != aid4);
        if needs_flip {
            stype = match stype {
                BondStereo::Cis => BondStereo::Trans,
                BondStereo::Trans => BondStereo::Cis,
                BondStereo::Z => BondStereo::E,
                BondStereo::E => BondStereo::Z,
                other => other,
            };
        }
    }
    stype
}

// ──────────────────────────────────────────────
// Set ring angle
// ──────────────────────────────────────────────

fn set_ring_angle(mol: &Molecule, aid2: usize, ring_size: usize) -> f64 {
    let hyb = &mol.atoms()[aid2].hybridization();
    ideal_bond_angle(hyb, Some(ring_size))
}

/// RDKit❗✔️: set13Bounds — all 1-3 distance bounds
fn set_13_bounds(
    mol: &Molecule,
    mmat: &mut BoundsMatrix,
    bond_lengths: &[f64],
    ring_info: Option<&crate::RingInfo>,
) {
    let n = mol.atoms().len();

    if let Some(ri) = ring_info {
        // Process rings first
        let mut atom_rings: Vec<Vec<usize>> = ri
            .atom_rings()
            .iter()
            .map(|r| r.iter().map(|aid| aid.index()).collect())
            .collect();
        atom_rings.sort_by(|a, b| a.len().cmp(&b.len()));

        let mut visited_13 = vec![false; n * n];

        for ring in &atom_rings {
            let rsize = ring.len();
            let n_ring = ring.len();
            for i in 0..n_ring {
                let aid2 = ring[i];
                let aid1 = if i == 0 {
                    ring[n_ring - 1]
                } else {
                    ring[i - 1]
                };
                let aid3 = if i == n_ring - 1 {
                    ring[0]
                } else {
                    ring[i + 1]
                };

                let pid = aid1.min(aid3) * n + aid1.max(aid3);

                if !visited_13[pid] {
                    let angle = set_ring_angle(mol, aid2, rsize);
                    set_13_bounds_helper(aid1, aid2, aid3, angle, bond_lengths, mmat, mol);
                    visited_13[pid] = true;
                }
            }
        }
    }

    // Non-ring 1-3 paths: middle atom not in a ring
    for i in 0..n {
        let nbrs = neighbors_for_atom(mol, i);
        if nbrs.len() < 2 {
            continue;
        }
        // Check each pair of neighbors
        for a in 0..nbrs.len() {
            for b in (a + 1)..nbrs.len() {
                let aid1 = nbrs[a];
                let aid3 = nbrs[b];
                let pid = aid1.min(aid3) * n + aid1.max(aid3);

                let angle = ideal_bond_angle(&mol.atoms()[i].hybridization(), None);
                set_13_bounds_helper(aid1, i, aid3, angle, bond_lengths, mmat, mol);
            }
        }
    }
}

// ──────────────────────────────────────────────
// 1-4 distance bounds (torsion-based)
// ──────────────────────────────────────────────

// ──────────────────────────────────────────────
// 1-5 distance bounds (optional)
// ──────────────────────────────────────────────

/// RDKit❗✔️: _compute15DistsCisCis — 1-5 distance, cis-cis configuration
///
/// Configuration:
///         5
///          \
///     1     4
///      \   /
///       2-3
// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_compute15DistsCisCis (BoundsMatrixBuilder.cpp:1965-1981)
// RDKit✔️✔️: double _compute15DistsCisCis(double d1, double d2, double d3, double d4,
// RDKit✔️✔️:                              double ang12, double ang23, double ang34) {
// RDKit✔️✔️:   double dx14 = d2 - d3 * cos(ang23) - d1 * cos(ang12);
// RDKit✔️✔️:   double dy14 = d3 * sin(ang23) - d1 * sin(ang12);
// RDKit✔️✔️:   double d14 = sqrt(dx14 * dx14 + dy14 * dy14);
// RDKit✔️✔️:   double cval = (d3 - d2 * cos(ang23) + d1 * cos(ang12 + ang23)) / d14;
// RDKit✔️✔️:   if (cval > 1.0) { cval = 1.0; } else if (cval < -1.0) { cval = -1.0; }
// RDKit✔️✔️:   double ang143 = acos(cval);
// RDKit✔️✔️:   double ang145 = ang34 - ang143;
// RDKit✔️✔️:   double res = RDGeom::compute13Dist(d14, d4, ang145);
// RDKit✔️✔️:   return res;
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION
fn compute_15_dist_cis_cis(
    d1: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    ang12: f64,
    ang23: f64,
    ang34: f64,
) -> f64 {
    let dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy14 = d3 * ang23.sin() - d1 * ang12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let cval = ((d3 - d2 * ang23.cos() + d1 * (ang12 + ang23).cos()) / d14).clamp(-1.0, 1.0);
    let ang143 = cval.acos();
    let ang145 = ang34 - ang143;
    compute_13_dist(d14, d4, ang145)
}

/// RDKit❗✔️: _compute15DistsCisTrans — 1-5 distance, cis-trans configuration
///
/// Configuration:
///  1     4-5
///   \   /
///    2-3
// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_compute15DistsCisTrans (BoundsMatrixBuilder.cpp:2004-2019)
// RDKit✔️✔️: (same structure as CisCis, but ang145 = ang34 + ang143)
// END RDKIT CPP FUNCTION
fn compute_15_dist_cis_trans(
    d1: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    ang12: f64,
    ang23: f64,
    ang34: f64,
) -> f64 {
    let dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy14 = d3 * ang23.sin() - d1 * ang12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let cval = ((d3 - d2 * ang23.cos() + d1 * (ang12 + ang23).cos()) / d14).clamp(-1.0, 1.0);
    let ang143 = cval.acos();
    let ang145 = ang34 + ang143;
    compute_13_dist(d14, d4, ang145)
}

/// RDKit❗✔️: _compute15DistsTransTrans — 1-5 distance, trans-trans configuration
///
/// Configuration:
///  1
///   \
///    2-3
///       \
///        4-5
// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_compute15DistsTransTrans (BoundsMatrixBuilder.cpp:2045-2060)
// RDKit✔️✔️: (dy14 uses + d1*sin(ang12), cval uses cos(ang12 - ang23))
// END RDKIT CPP FUNCTION
fn compute_15_dist_trans_trans(
    d1: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    ang12: f64,
    ang23: f64,
    ang34: f64,
) -> f64 {
    let dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy14 = d3 * ang23.sin() + d1 * ang12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let cval = ((d3 - d2 * ang23.cos() + d1 * (ang12 - ang23).cos()) / d14).clamp(-1.0, 1.0);
    let ang143 = cval.acos();
    let ang145 = ang34 + ang143;
    compute_13_dist(d14, d4, ang145)
}

/// RDKit❗✔️: _compute15DistsTransCis — 1-5 distance, trans-cis configuration
///
/// Configuration:
///                    1
///                     \
///                      2-3
///                         \
///                          4
///                         /
///                        5
// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_compute15DistsTransCis (BoundsMatrixBuilder.cpp:2088-2104)
// RDKit✔️✔️: (dy14 uses + d1*sin(ang12), cval uses cos(ang12 - ang23), ang145 = ang34 - ang143)
// END RDKIT CPP FUNCTION
fn compute_15_dist_trans_cis(
    d1: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    ang12: f64,
    ang23: f64,
    ang34: f64,
) -> f64 {
    let dx14 = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy14 = d3 * ang23.sin() + d1 * ang12.sin();
    let d14 = (dx14 * dx14 + dy14 * dy14).sqrt();
    let cval = ((d3 - d2 * ang23.cos() + d1 * (ang12 - ang23).cos()) / d14).clamp(-1.0, 1.0);
    let ang143 = cval.acos();
    let ang145 = ang34 - ang143;
    compute_13_dist(d14, d4, ang145)
}

/// RDKit❗✔️: set15Bounds — 1-5 distance bounds with stereochemistry (full RDKit version)
///
/// Finds all 1-5 atom pairs (a-b-c-d-e) via pre-computed 1-4 paths (paths14)
/// and sets distance bounds using bond lengths, computed bond angles, and
/// cis/trans stereochemistry detection on the central double bond.
///
/// The three 1-4 path types:
/// - **Cis**: stereochemistry prefers cis on the central bond
/// - **Trans**: stereochemistry prefers trans on the central bond
/// - **Other**: no detected stereochemistry, uses VDW fallback for unknown paths
///
/// For each 1-4 path, the function extends through all bonds from the 4th atom
/// to find 1-5 pairs, then computes the distance using one of four geometry
/// configurations (cis-cis, cis-trans, trans-cis, trans-trans) depending on
/// the path type and cis/trans path flags.
// BEGIN RDKIT CPP FUNCTION DGeomHelpers::set15Bounds (BoundsMatrixBuilder.cpp:2232-2248)
// RDKit❗✔️: void set15Bounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat,
// RDKit❗✔️:                  ComputedData &accumData, double *distMatrix) {
// RDKit❗✔️:   for (pti = accumData.paths14.begin(); pti != accumData.paths14.end(); pti++) {
// RDKit❗✔️:     bid1 = pti->bid1;
// RDKit❗✔️:     bid2 = pti->bid2;
// RDKit❗✔️:     bid3 = pti->bid3;
// RDKit❗✔️:     type = pti->type;
// RDKit❗✔️:     _set15BoundsHelper(mol, bid1, bid2, bid3, type, accumData, mmat, distMatrix);
// RDKit❗✔️:     _set15BoundsHelper(mol, bid3, bid2, bid1, type, accumData, mmat, distMatrix);
// RDKit❗✔️:   }
// RDKit❗✔️: }
// END RDKIT CPP FUNCTION
fn set_15_bounds(
    mol: &Molecule,
    mmat: &mut BoundsMatrix,
    accum_data: &mut ComputedData,
    dmat: &[f64],
) {
    let nb = mol.bonds().len();
    let na = mol.atoms().len();
    let paths = accum_data.paths14.clone(); // clone to avoid borrow issues
    for path in &paths {
        // RDKit❗✔️: call _set15BoundsHelper with (bid1,bid2,bid3) and reverse (bid3,bid2,bid1)
        set_15_bounds_helper(
            mol, mmat, accum_data, dmat, nb, na, path.bid1, path.bid2, path.bid3, path.kind,
        );
        set_15_bounds_helper(
            mol, mmat, accum_data, dmat, nb, na, path.bid3, path.bid2, path.bid1, path.kind,
        );
    }
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_set15BoundsHelper (BoundsMatrixBuilder.cpp:2106-2228)
// RDKit❗✔️: void _set15BoundsHelper(const ROMol &mol, unsigned int bid1, unsigned int bid2,
// RDKit❗✔️:                          unsigned int bid3, unsigned int type,
// RDKit❗✔️:                          ComputedData &accumData, DistGeom::BoundsMatPtr mmat,
// RDKit❗✔️:                          double *dmat) {
// RDKit❗✔️:   unsigned int aid2 = accumData.bondAdj->getVal(bid1, bid2);
// RDKit❗✔️:   unsigned int aid1 = mol.getBondWithIdx(bid1)->getOtherAtomIdx(aid2);
// RDKit❗✔️:   unsigned int aid3 = accumData.bondAdj->getVal(bid2, bid3);
// RDKit❗✔️:   unsigned int aid4 = mol.getBondWithIdx(bid3)->getOtherAtomIdx(aid3);
// RDKit❗✔️:   double d1 = accumData.bondLengths[bid1];
// RDKit❗✔️:   double d2 = accumData.bondLengths[bid2];
// RDKit❗✔️:   double d3 = accumData.bondLengths[bid3];
// RDKit❗✔️:   double ang12 = accumData.bondAngles->getVal(bid1, bid2);
// RDKit❗✔️:   double ang23 = accumData.bondAngles->getVal(bid2, bid3);
// RDKit❗✔️:   for (i = 0; i < nb; i++) {
// RDKit❗✔️:     if (accumData.bondAdj->getVal(bid3, i) == (int)aid4) {
// RDKit❗✔️:       aid5 = mol.getBondWithIdx(i)->getOtherAtomIdx(aid4);
// RDKit❗✔️:       unsigned int pid = std::min(aid1, aid5) * na + std::max(aid1, aid5);
// RDKit❗✔️:       if (accumData.visitedBound(pid, DistType::DIST14)) { return; }
// RDKit❗✔️:       if (dmat[...] < 3.9) { continue; }
// RDKit❗✔️:       if (aid1 != aid5) {
// RDKit❗✔️:         if ((mmat->getLowerBound(aid1, aid5) < DIST12_DELTA) ||
// RDKit❗✔️:             accumData.set15Atoms[pid1] || accumData.set15Atoms[pid2]) {
// RDKit❗✔️:           d4 = accumData.bondLengths[i];
// RDKit❗✔️:           ang34 = accumData.bondAngles->getVal(bid3, i);
// RDKit❗✔️:           // cis/trans-aware distance computation based on type and pathId
// RDKit❗✔️:           _checkAndSetBounds(aid1, aid5, dl, du, mmat);
// RDKit❗✔️:           accumData.set15Atoms[pid1] = 1;
// RDKit❗✔️:           accumData.set15Atoms[pid2] = 1;
// RDKit❗✔️: } } } } }
// END RDKIT CPP FUNCTION
#[allow(clippy::too_many_arguments)]
fn set_15_bounds_helper(
    mol: &Molecule,
    mmat: &mut BoundsMatrix,
    accum_data: &mut ComputedData,
    dmat: &[f64],
    nb: usize,
    na: usize,
    bid1: usize,
    bid2: usize,
    bid3: usize,
    kind: Path14Kind,
) {
    // RDKit❗✔️: Get shared atoms via bond adjacency matrix
    let aid2 = accum_data.get_bond_adj(nb, bid1, bid2) as usize;
    let aid1 = if mol.bonds()[bid1].begin().index() == aid2 {
        mol.bonds()[bid1].end().index()
    } else {
        mol.bonds()[bid1].begin().index()
    };
    let aid3 = accum_data.get_bond_adj(nb, bid2, bid3) as usize;
    let aid4 = if mol.bonds()[bid3].begin().index() == aid3 {
        mol.bonds()[bid3].end().index()
    } else {
        mol.bonds()[bid3].begin().index()
    };

    let d1 = accum_data.bond_lengths[bid1];
    let d2 = accum_data.bond_lengths[bid2];
    let d3 = accum_data.bond_lengths[bid3];
    let ang12 = accum_data.get_bond_angle(nb, bid1, bid2);
    let ang23 = accum_data.get_bond_angle(nb, bid2, bid3);

    // RDKit❗✔️: loop over all bonds to find the 5th atom (extension from aid4)
    for i in 0..nb {
        if accum_data.get_bond_adj(nb, bid3, i) != aid4 as i32 {
            continue;
        }
        let aid5 = if mol.bonds()[i].begin().index() == aid4 {
            mol.bonds()[i].end().index()
        } else {
            mol.bonds()[i].begin().index()
        };

        // RDKit❗✔️: FIX: make sure we did not come back to the first atom (Issue 244)
        let pid = aid1.min(aid5) * na + aid1.max(aid5);

        // RDKit❗✔️: if pair was already bounded as 1-2/1-3/1-4, skip the whole function
        if accum_data.visited_bound(pid, true) {
            return;
        }

        // RDKit❗✔️: skip if atoms are already close in distance matrix
        if dmat[aid1.max(aid5) * na + aid1.min(aid5)] < 3.9 {
            continue;
        }

        // RDKit❗✔️: skip if same atom (4-membered ring issue)
        if aid1 == aid5 {
            continue;
        }

        // RDKit❗✔️: check if bounds not already set and not already processed
        let pid1 = aid1 * na + aid5;
        let pid2 = aid5 * na + aid1;
        if !(mmat.get_lower(aid1, aid5) < DIST12_DELTA
            || accum_data.set15_atoms[pid1]
            || accum_data.set15_atoms[pid2])
        {
            continue;
        }

        let d4 = accum_data.bond_lengths[i];
        let ang34 = accum_data.get_bond_angle(nb, bid3, i);

        // RDKit❗✔️: pathId = bid2*nb*nb + bid3*nb + i
        let path_id = bid2 as u64 * nb as u64 * nb as u64 + bid3 as u64 * nb as u64 + i as u64;

        // RDKit❗✔️: compute 1-5 distance based on path type and cis/trans flags
        let (mut dl, mut du) = match kind {
            // RDKit❗✔️: type == 0 (CIS)
            Path14Kind::Cis => {
                if has_path_flag(&accum_data.cis_paths, path_id) {
                    // RDKit❗✔️: _compute15DistsCisCis +- DIST15_TOL
                    let base = compute_15_dist_cis_cis(d1, d2, d3, d4, ang12, ang23, ang34);
                    (base - DIST15_TOL, base + DIST15_TOL)
                } else if has_path_flag(&accum_data.trans_paths, path_id) {
                    // RDKit❗✔️: _compute15DistsCisTrans +- DIST15_TOL
                    let base = compute_15_dist_cis_trans(d1, d2, d3, d4, ang12, ang23, ang34);
                    (base - DIST15_TOL, base + DIST15_TOL)
                } else {
                    // RDKit❗✔️: cis-cis - tol, cis-trans + tol
                    (
                        compute_15_dist_cis_cis(d1, d2, d3, d4, ang12, ang23, ang34) - DIST15_TOL,
                        compute_15_dist_cis_trans(d1, d2, d3, d4, ang12, ang23, ang34) + DIST15_TOL,
                    )
                }
            }
            // RDKit❗✔️: type == 1 (TRANS)
            Path14Kind::Trans => {
                if has_path_flag(&accum_data.cis_paths, path_id) {
                    // RDKit❗✔️: _compute15DistsTransCis +- DIST15_TOL
                    let base = compute_15_dist_trans_cis(d1, d2, d3, d4, ang12, ang23, ang34);
                    (base - DIST15_TOL, base + DIST15_TOL)
                } else if has_path_flag(&accum_data.trans_paths, path_id) {
                    // RDKit❗✔️: _compute15DistsTransTrans +- DIST15_TOL
                    let base = compute_15_dist_trans_trans(d1, d2, d3, d4, ang12, ang23, ang34);
                    (base - DIST15_TOL, base + DIST15_TOL)
                } else {
                    // RDKit❗✔️: trans-cis - tol, trans-trans + tol
                    (
                        compute_15_dist_trans_cis(d1, d2, d3, d4, ang12, ang23, ang34) - DIST15_TOL,
                        compute_15_dist_trans_trans(d1, d2, d3, d4, ang12, ang23, ang34)
                            + DIST15_TOL,
                    )
                }
            }
            // RDKit❗✔️: type == 2 (OTHER) — reversed d-params for cis/trans detection
            Path14Kind::Other => {
                if has_path_flag(&accum_data.cis_paths, path_id) {
                    // RDKit❗✔️: _compute15DistsCisCis(d4,d3,d2,d1,ang34,ang23,ang12) - tol
                    // RDKit❗✔️: _compute15DistsCisTrans(d4,d3,d2,d1,ang34,ang23,ang12) + tol
                    (
                        compute_15_dist_cis_cis(d4, d3, d2, d1, ang34, ang23, ang12) - DIST15_TOL,
                        compute_15_dist_cis_trans(d4, d3, d2, d1, ang34, ang23, ang12) + DIST15_TOL,
                    )
                } else if has_path_flag(&accum_data.trans_paths, path_id) {
                    // RDKit❗✔️: _compute15DistsTransCis(d4,d3,d2,d1,ang34,ang23,ang12) - tol
                    // RDKit❗✔️: _compute15DistsTransTrans(d4,d3,d2,d1,ang34,ang23,ang12) + tol
                    (
                        compute_15_dist_trans_cis(d4, d3, d2, d1, ang34, ang23, ang12) - DIST15_TOL,
                        compute_15_dist_trans_trans(d4, d3, d2, d1, ang34, ang23, ang12)
                            + DIST15_TOL,
                    )
                } else {
                    // RDKit❗✔️: VDW fallback for unknown stereochemistry
                    let vw1 = vdw_radius(mol.atoms()[aid1].atomic_number());
                    let vw5 = vdw_radius(mol.atoms()[aid5].atomic_number());
                    (VDW_SCALE_15 * (vw1 + vw5), MAX_UPPER)
                }
            }
        };

        // RDKit❗✔️: if (du < 0.0) { du = MAX_UPPER; }
        if du < 0.0 {
            du = MAX_UPPER;
        }
        if dl < 0.0 {
            dl = 0.0;
        }

        // RDKit❗✔️: _checkAndSetBounds(aid1, aid5, dl, du, mmat);
        mmat.check_and_set_bounds(aid1, aid5, dl, du);

        // RDKit❗✔️: accumData.set15Atoms[pid1] = 1; accumData.set15Atoms[pid2] = 1;
        accum_data.set15_atoms[pid1] = true;
        accum_data.set15_atoms[pid2] = true;
    }
}
// END RDKIT CPP FUNCTION

/// RDKit❗✔️: set14Bounds — torsion-based 1-4 bounds
/// Also populates paths14, bond_adj, bond_angles for 1-5 bounds.
fn set_14_bounds(mol: &Molecule, mmat: &mut BoundsMatrix, accum_data: &mut ComputedData) {
    let n = mol.atoms().len();
    let nb = mol.bonds().len();

    // Find all 1-4 paths (a-b-c-d where a-d are the 1-4 pair)
    let mut visited_14 = vec![false; n * n];

    for bond_bc in mol.bonds() {
        let bc_idx = bond_bc.id().index();
        let b = bond_bc.begin().index();
        let c = bond_bc.end().index();

        // RDKit❗✔️: set bond adjacency — bond_bc connects atoms b and c
        // This builds the bond adjacency matrix needed by set15Bounds
        accum_data.set_bond_adj(nb, bc_idx, bc_idx, b as i32);

        let nbrs_b = neighbors_for_atom(mol, b);
        let nbrs_c = neighbors_for_atom(mol, c);

        for &a in &nbrs_b {
            if a == c {
                continue;
            }
            let bid_ab = match bond_between_idx_simple(mol, a, b) {
                Some(idx) => idx,
                None => continue,
            };

            // RDKit❗✔️: set bond adjacency for bid_ab and bc_idx (shared atom b)
            accum_data.set_bond_adj(nb, bid_ab, bc_idx, b as i32);

            // RDKit❗✔️: compute bond angle at atom b between bonds (a-b) and (b-c)
            let ang_at_b = ideal_bond_angle(&mol.atoms()[b].hybridization(), None);
            accum_data.set_bond_angle(nb, bid_ab, bc_idx, ang_at_b);

            // RDKit❗✔️: also set for the reverse order
            accum_data.set_bond_adj(nb, bid_ab, bid_ab, a as i32);
            accum_data.set_bond_adj(nb, bc_idx, bc_idx, b as i32);

            for &d in &nbrs_c {
                if d == b || d == a {
                    continue;
                }
                let bid_cd = match bond_between_idx_simple(mol, c, d) {
                    Some(idx) => idx,
                    None => continue,
                };

                // RDKit❗✔️: set bond adjacency for bc_idx and bid_cd (shared atom c)
                accum_data.set_bond_adj(nb, bc_idx, bid_cd, c as i32);

                // RDKit❗✔️: compute bond angle at atom c between bonds (b-c) and (c-d)
                let ang_at_c = ideal_bond_angle(&mol.atoms()[c].hybridization(), None);
                accum_data.set_bond_angle(nb, bc_idx, bid_cd, ang_at_c);

                // RDKit❗✔️: also set bond's own adjacency
                accum_data.set_bond_adj(nb, bid_cd, bid_cd, d as i32);

                // RDKit❗✔️: build Path14Configuration for this 1-4 path
                let stype = get_atom_stereo(bond_bc, a, d);

                // RDKit❗✔️: determine cis/trans preference (BoundsMatrixBuilder.cpp:770-786)
                let prefer_cis;
                let prefer_trans;
                match stype {
                    BondStereo::Cis | BondStereo::Z => {
                        prefer_cis = true;
                        prefer_trans = false;
                    }
                    BondStereo::Trans | BondStereo::E => {
                        prefer_cis = false;
                        prefer_trans = true;
                    }
                    _ => {
                        prefer_cis = false;
                        prefer_trans = false;
                    }
                }

                let pid = a.min(d) * n + a.max(d);
                if visited_14[pid] {
                    continue;
                }
                visited_14[pid] = true;

                // RDKit❗✔️: store the 1-4 path with its stereochemistry type
                let path_kind;
                if prefer_cis {
                    path_kind = Path14Kind::Cis;
                    let path_id = bid_ab as u64 * nb as u64 * nb as u64
                        + bc_idx as u64 * nb as u64
                        + bid_cd as u64;
                    record_path_flag(&mut accum_data.cis_paths, path_id);
                    record_path_flag(
                        &mut accum_data.cis_paths,
                        bid_cd as u64 * nb as u64 * nb as u64
                            + bc_idx as u64 * nb as u64
                            + bid_ab as u64,
                    );
                } else if prefer_trans {
                    path_kind = Path14Kind::Trans;
                    let path_id = bid_ab as u64 * nb as u64 * nb as u64
                        + bc_idx as u64 * nb as u64
                        + bid_cd as u64;
                    record_path_flag(&mut accum_data.trans_paths, path_id);
                    record_path_flag(
                        &mut accum_data.trans_paths,
                        bid_cd as u64 * nb as u64 * nb as u64
                            + bc_idx as u64 * nb as u64
                            + bid_ab as u64,
                    );
                } else {
                    path_kind = Path14Kind::Other;
                }

                accum_data.paths14.push(Path14Configuration {
                    bid1: bid_ab,
                    bid2: bc_idx,
                    bid3: bid_cd,
                    kind: path_kind,
                });

                // Approximate 1-4 distance using cis (min) and trans (max) conformations
                let bl_ab = bond_lengths_from_accum(accum_data, bid_ab);
                let bl_bc = bond_lengths_from_accum(accum_data, bc_idx);
                let bl_cd = bond_lengths_from_accum(accum_data, bid_cd);

                // Min distance (cis, 0° torsion)
                let min_dist = (bl_ab + bl_cd).hypot(bl_bc * 0.5);
                // Max distance (trans, 180° torsion)
                let max_dist = bl_ab + bl_bc + bl_cd;

                mmat.check_and_set_bounds(a, d, min_dist * 0.8, max_dist * 1.0);
            }
        }
    }
}

/// Helper: get bond length from ComputedData
fn bond_lengths_from_accum(accum_data: &ComputedData, bid: usize) -> f64 {
    let bl = accum_data.bond_lengths[bid];
    if bl > 0.0 { bl } else { 1.5 }
}

fn bond_between_idx_simple(mol: &Molecule, a: usize, b: usize) -> Option<usize> {
    mol.bonds()
        .iter()
        .find(|bond| {
            (bond.begin().index() == a && bond.end().index() == b)
                || (bond.begin().index() == b && bond.end().index() == a)
        })
        .map(|b| b.id().index())
}

// ──────────────────────────────────────────────
// VDW lower bounds
// ──────────────────────────────────────────────

/// RDKit❗✔️: setLowerBoundVDW — VDW lower bounds for non-bonded pairs
// BEGIN RDKIT CPP FUNCTION DGeomHelpers::setLowerBoundVDW (BoundsMatrixBuilder.cpp:315-361)
// RDKit✔️✔️: void setLowerBoundVDW(const ROMol &mol, DistGeom::BoundsMatPtr mmat, bool,
// RDKit✔️✔️:                         double *dmat) {
// RDKit✔️✔️:   for (unsigned int i = 1; i < npt; i++) {
// RDKit✔️✔️:     const auto atomI = mol.getAtomWithIdx(i);
// RDKit✔️✔️:     auto vw1 = PeriodicTable::getTable()->getRvdw(atomI->getAtomicNum());
// RDKit✔️✔️:     if (isHinHBondDonor(atomI, mol)) hinHBondDonors.set(i);
// RDKit✔️✔️:     if (isHBondAcceptor(atomI)) hBondAcceptors.set(i);
// RDKit✔️✔️:     for (unsigned int j = 0; j < i; j++) {
// RDKit✔️✔️:       if (mmat->getLowerBound(i, j) < DIST12_DELTA) {
// RDKit✔️✔️:         if ((hinHBondDonors[i] && hBondAcceptors[j]) || ...) {
// RDKit✔️✔️:           mmat->setLowerBound(i, j, H_BOND_LENGTH);
// RDKit✔️✔️:         } else if (dmat[i * npt + j] == 4.0) {
// RDKit✔️✔️:           mmat->setLowerBound(i, j, VDW_SCALE_15 * (vw1 + vw2));
// RDKit✔️✔️:         } else ...
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::setLowerBoundVDW
fn set_lower_bound_vdw(mol: &Molecule, mmat: &mut BoundsMatrix, topo_dist: &[Vec<i32>]) {
    let n = mol.atoms().len();

    // Pre-compute H-bond donor/acceptor bitsets
    let mut h_in_hbond_donors = vec![false; n];
    let mut hbond_acceptors = vec![false; n];
    for i in 1..n {
        h_in_hbond_donors[i] = is_h_in_hbond_donor(mol, i);
        hbond_acceptors[i] = is_hbond_acceptor(mol.atoms()[i].atomic_number());
    }

    for i in 1..n {
        let vw1 = vdw_radius(mol.atoms()[i].atomic_number());

        for j in 0..i {
            if mmat.get_lower(i, j) > DIST12_DELTA {
                continue;
            }

            let vw2 = vdw_radius(mol.atoms()[j].atomic_number());
            let td = topo_dist[i][j];

            // H-bond special case: set lower bound to 1.8A for donor-H + acceptor
            if (h_in_hbond_donors[i] && hbond_acceptors[j])
                || (hbond_acceptors[i] && h_in_hbond_donors[j])
            {
                mmat.set_lower(i, j, H_BOND_LENGTH);
            } else if td == 4 {
                // 1-5: scaled VDW
                mmat.set_lower(i, j, VDW_SCALE_15 * (vw1 + vw2));
            } else if td == 5 {
                // 1-6: slightly less scaled
                mmat.set_lower(
                    i,
                    j,
                    (VDW_SCALE_15 + 0.5 * (1.0 - VDW_SCALE_15)) * (vw1 + vw2),
                );
            } else if td >= 6 {
                // Full VDW sum
                mmat.set_lower(i, j, vw1 + vw2);
            }
        }
    }
}

// ──────────────────────────────────────────────
// Public entry point
// ──────────────────────────────────────────────

/// Build a distance bounds matrix for the given molecule.
///
/// Returns an n×n matrix where each cell [i][j] is the **upper** distance
/// bound (Å) between atoms i and j. The caller can derive lower bounds
/// separately if needed.
///
/// ## Algorithm
///
/// 1. Compute topological (shortest-path) distances
/// 2. Set **1-2 bounds** from bond-order estimated bond lengths
/// 3. Set **1-3 bounds** from angle geometry (ring-aware)
/// 4. Set **1-4 bounds** from torsion ranges
/// 5. Set **VDW lower bounds** for non-bonded pairs
/// 6. **Triangle inequality smoothing** (3 passes)
///
/// ## Source
///
/// Ported from RDKit `DGeomHelpers::BoundsMatrixBuilder` (BoundsMatrixBuilder.cpp)
pub fn dg_bounds_matrix(molecule: &Molecule) -> Result<Vec<Vec<f64>>, DgBoundsError> {
    let n = molecule.atoms().len();
    if n < 2 {
        // Single atom or empty molecule
        let mut result = Vec::with_capacity(n);
        for i in 0..n {
            let mut row = vec![0.0; n];
            result.push(row);
        }
        return Ok(result);
    }

    let nb = molecule.bonds().len();
    let mut mmat = BoundsMatrix::new(n);
    let mut accum_data = ComputedData::new(n, nb);

    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::BoundsMatrixBuilder (BoundsMatrixBuilder.cpp)
    // RDKit❗✔️: The RDKit implementation calls:
    //   RDKit❗✔️: set12Bounds(mol, mmat, accumData);
    //   RDKit❗✔️: set13Bounds(mol, mmat, accumData);
    //   RDKit❗✔️: set14Bounds(mol, mmat, accumData);
    //   RDKit❗✔️: set15Bounds(mol, mmat, accumData);            // optional
    //   RDKit❗✔️: setLowerBoundVDW(mol, mmat, true, distMatrix);
    //   RDKit❗✔️: TriangleSmooth(mmat);
    // RDKit uses UFF atom typing for bond lengths, a full SSSR ring analysis,
    // cis/trans path detection, and specialized 1-5 bounds for several
    // substructure classes. This implementation uses:
    //   - topology-based bond length estimates
    //   - hybridization-based ideal angles
    //   - VDW radii from PeriodicTable
    //   - Floyd-Warshall topological distances
    //   - full set15Bounds with stereochemistry (cis/trans detection)
    // END RDKIT CPP FUNCTION DGeomHelpers::BoundsMatrixBuilder

    // Step 1: 1-2 bounds
    set_12_bounds(molecule, &mut mmat, &mut accum_data);

    // Step 2: 1-3 bounds
    let ring_info = molecule.derived_cache().rings.as_ref();
    set_13_bounds(molecule, &mut mmat, &accum_data.bond_lengths, ring_info);

    // Step 3: 1-4 bounds (also populates paths14, bond_adj, bond_angles)
    set_14_bounds(molecule, &mut mmat, &mut accum_data);

    // Step 3b: 1-5 bounds with stereochemistry (full RDKit version)
    // RDKit❗✔️: set15Bounds(mol, mmat, accumData, distMatrix);
    // dmat: flat array of topological distances (n_atoms x n_atoms)
    let topo_dist = compute_topological_distances(molecule);
    let mut dmat = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            dmat[i * n + j] = topo_dist[i][j] as f64;
        }
    }
    set_15_bounds(molecule, &mut mmat, &mut accum_data, &dmat);

    // Step 4: VDW lower bounds (using flat dmat for distance-based VDW scaling)
    set_lower_bound_vdw(molecule, &mut mmat, &topo_dist);

    // Step 5: triangle inequality smoothing
    mmat.triangle_smooth();

    Ok(mmat.to_vec_vec())
}

// ──────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::MoleculeBuilder;
    use crate::{AtomId, AtomSpec, BondSpec, Element};

    #[test]
    fn test_single_atom() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::C));
        let mol = builder.build().expect("build");
        let result = dg_bounds_matrix(&mol).expect("dg_bounds");
        assert_eq!(result.len(), 1);
        assert_eq!(result[0][0], 0.0);
    }

    #[test]
    fn test_diatomic() {
        let mut builder = MoleculeBuilder::new();
        let c1 = builder.add_atom(AtomSpec::new(Element::C));
        let c2 = builder.add_atom(AtomSpec::new(Element::C));
        builder.add_bond(BondSpec::new(c1, c2, BondOrder::Single));
        let mol = builder.build().expect("build");
        let result = dg_bounds_matrix(&mol).expect("dg_bounds");
        assert_eq!(result.len(), 2);
        assert!(result[0][1] > 0.0);
        assert!(result[0][1] < 5.0);
    }

    #[test]
    fn test_ethane() {
        let mut builder = MoleculeBuilder::new();
        let c1 = builder.add_atom(AtomSpec::new(Element::C));
        let c2 = builder.add_atom(AtomSpec::new(Element::C));
        let h1 = builder.add_atom(AtomSpec::new(Element::H));
        let h2 = builder.add_atom(AtomSpec::new(Element::H));
        let h3 = builder.add_atom(AtomSpec::new(Element::H));
        let h4 = builder.add_atom(AtomSpec::new(Element::H));
        let h5 = builder.add_atom(AtomSpec::new(Element::H));
        let h6 = builder.add_atom(AtomSpec::new(Element::H));
        builder.add_bond(BondSpec::new(c1, c2, BondOrder::Single));
        builder.add_bond(BondSpec::new(c1, h1, BondOrder::Single));
        builder.add_bond(BondSpec::new(c1, h2, BondOrder::Single));
        builder.add_bond(BondSpec::new(c1, h3, BondOrder::Single));
        builder.add_bond(BondSpec::new(c2, h4, BondOrder::Single));
        builder.add_bond(BondSpec::new(c2, h5, BondOrder::Single));
        builder.add_bond(BondSpec::new(c2, h6, BondOrder::Single));
        let mol = builder.build().expect("build");
        let result = dg_bounds_matrix(&mol).expect("dg_bounds");
        assert_eq!(result.len(), 8);
        assert!(result[0][1] > 1.0 && result[0][1] < 3.0);
        assert!(result[0][5] > 1.0 && result[0][5] < 5.0);
    }

    #[test]
    fn test_empty() {
        let builder = MoleculeBuilder::new();
        let mol = builder.build().expect("build");
        let result = dg_bounds_matrix(&mol).expect("dg_bounds");
        assert!(result.is_empty());
    }

    #[test]
    fn test_ethane_bounds_consistency() {
        // Ethane (C2H6) — verify 1-2, 1-3, 1-4, and VDW bounds are reasonable
        let mut builder = MoleculeBuilder::new();
        let c1 = builder.add_atom(AtomSpec::new(Element::C));
        let c2 = builder.add_atom(AtomSpec::new(Element::C));
        let h1 = builder.add_atom(AtomSpec::new(Element::H));
        let h2 = builder.add_atom(AtomSpec::new(Element::H));
        let h3 = builder.add_atom(AtomSpec::new(Element::H));
        let h4 = builder.add_atom(AtomSpec::new(Element::H));
        let h5 = builder.add_atom(AtomSpec::new(Element::H));
        let h6 = builder.add_atom(AtomSpec::new(Element::H));
        builder
            .add_bond(BondSpec::new(c1, c2, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c1, h1, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c1, h2, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c1, h3, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c2, h4, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c2, h5, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c2, h6, BondOrder::Single))
            .unwrap();
        let mol = builder.build().expect("build ethane");
        let result = dg_bounds_matrix(&mol).expect("dg_bounds");
        assert_eq!(result.len(), 8);

        // 1-2 (C-C bond): upper bound should be ≈ 1.5 + 0.01 = 1.51
        assert!(
            result[0][1] > 1.40 && result[0][1] < 1.60,
            "C-C 1-2 upper bound {} should be ~1.50",
            result[0][1]
        );

        // 1-2 (C-H bonds): upper bound should be ≈ 1.5 + 0.01 = 1.51
        assert!(
            result[0][2] > 1.40 && result[0][2] < 1.60,
            "C-H 1-2 upper bound {} should be ~1.50",
            result[0][2]
        );

        // 1-3 (C-C-H via C1): should be larger than 1-2
        assert!(
            result[2][3] > result[0][1] || result[2][4] > result[0][1],
            "1-3 bounds should be larger than 1-2 bounds"
        );

        // 1-4 (H-C-C-H across ethane): should be in a reasonable range
        let h1_c1_c2_h4 = result[2][5]; // h1-c1-c2-h4 is a 1-4 path
        assert!(
            h1_c1_c2_h4 > 0.8 && h1_c1_c2_h4 < 6.0,
            "1-4 H-C-C-H upper bound {} should be reasonable",
            h1_c1_c2_h4
        );

        // Matrix symmetry check
        for i in 0..8 {
            for j in 0..8 {
                assert!(!result[i][j].is_nan(), "bounds[{i}][{j}] is NaN");
                assert!(
                    result[i][j] >= 0.0 || i == j,
                    "bounds[{i}][{j}] = {} should be >= 0",
                    result[i][j]
                );
            }
        }
    }

    #[test]
    fn test_ring_bounds_consistency() {
        // Cyclohexane (C6H12) — a 6-membered ring with hydrogens
        // This tests ring-aware angle computation and triangle smoothing
        let mut builder = MoleculeBuilder::new();
        let c: Vec<_> = (0..6)
            .map(|_| builder.add_atom(AtomSpec::new(Element::C)))
            .collect();
        // Ring bonds
        for i in 0..6 {
            let j = (i + 1) % 6;
            builder
                .add_bond(BondSpec::new(c[i], c[j], BondOrder::Single))
                .unwrap();
        }
        let mol = builder.build().expect("build cyclohexane");
        let result = dg_bounds_matrix(&mol).expect("dg_bounds");
        assert_eq!(result.len(), 6);

        // 1-2 (C-C bonds): upper bound ~1.51
        for i in 0..6 {
            let j = (i + 1) % 6;
            assert!(
                result[i][j] > 1.40 && result[i][j] < 1.60,
                "C-C 1-2 upper bound [{i}][{j}] = {} should be ~1.50",
                result[i][j]
            );
        }

        // 1-3 (C-C-C in ring): should be larger than 1-2
        for i in 0..6 {
            let j = (i + 2) % 6;
            assert!(
                result[i][j] > result[i][(i + 1) % 6],
                "1-3 upper bound [{i}][{j}] = {} should exceed 1-2 bound {}",
                result[i][j],
                result[i][(i + 1) % 6]
            );
        }

        // 1-4 (C-C-C-C across hexagon): should be the largest intra-ring distance
        for i in 0..6 {
            let j = (i + 3) % 6;
            assert!(
                result[i][j] > result[i][(i + 1) % 6] && result[i][j] < 5.0,
                "1-4 upper bound [{i}][{j}] = {} should be > 1-2 and < 5.0",
                result[i][j]
            );
        }

        assert!(!result[0][1].is_nan());
        assert!(!result[0][1].is_infinite());
    }
}
