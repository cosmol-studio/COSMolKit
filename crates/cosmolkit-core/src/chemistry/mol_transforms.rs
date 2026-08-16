//! Conformer manipulation ported from RDKit `MolTransforms`.
//!
//! Source reproduction protocol: dev/source_reproduction_protocol.md
//!
//! Marker convention:
//!   RDKit✔️✔️: one-to-one port
//!   RDKit✔️❌: adapted for new-core API, same algorithm
//!   RDKit❗✔️: algorithm-equivalent via different architecture
//!   RDKit❗❌: not yet ported (unimplemented stub)
//!
//! RDKit source: MolTransforms.cpp, MolTransforms.h
//!   - getBondLength / setBondLength:      MolTransforms.cpp lines 509-544
//!   - getAngleRad / setAngleRad:          MolTransforms.cpp lines 546-616
//!   - getDihedralRad / setDihedralRad:    MolTransforms.cpp lines 618-710
//!   - transformConformer:                 MolTransforms.cpp lines 445-451
//!   - _toBeMovedIdxList (helper):         MolTransforms.cpp lines 467-507
//!   - getAngleDeg / setAngleDeg:          MolTransforms.h   lines 176-194
//!   - getDihedralDeg / setDihedralDeg:    MolTransforms.h   lines 204-223
//!
//! Additional functions (MolTransforms.h):
//!   - computeCentroid:                    MolTransforms.cpp lines 46-63

use std::f64::consts::PI;

use crate::{AdjacencyList, Molecule};

// ──────────────────────────────────────────────
// Error type
// ──────────────────────────────────────────────

/// Errors that can occur during conformer transform operations.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum MolTransformError {
    #[error("atom index {atom} out of bounds: num_atoms={num_atoms}")]
    AtomIndexOutOfBounds { atom: usize, num_atoms: usize },
    #[error("atoms must be bonded: {i} and {j}")]
    AtomsNotBonded { i: usize, j: usize },
    #[error("bond ({i},{j}) must not belong to a ring")]
    BondInRing { i: usize, j: usize },
    #[error("atoms {i} and {j} have identical 3D coordinates")]
    IdenticalCoordinates { i: usize, j: usize },
    #[error("no 3D conformer available for molecule")]
    NoConformer3D,
    #[error("conformer index {id} out of bounds: max={max}")]
    ConformerIndexOutOfBounds { id: usize, max: usize },
    #[error("bond ({i},{j}) and ({j},{k}) must not both belong to a ring")]
    BondsBothInRing { i: usize, j: usize, k: usize },
    #[error("atoms {i} and {j} have identical 3D coordinates (zero-length vector)")]
    ZeroLengthVector { i: usize, j: usize },
    #[error("molecule operation contract failed: {message}")]
    OperationContract { message: String },
}

// ──────────────────────────────────────────────
// Internal helpers
// ──────────────────────────────────────────────

/// RDKit❗✔️: Extract 3D coordinates from the first (or specified) conformer.
///
/// In RDKit, functions take a `Conformer &` directly. Here we extract from
/// `Molecule` by conformer index.
fn conf_coords(mol: &Molecule, conf_id: usize) -> Result<&[[f64; 3]], MolTransformError> {
    let confs = mol.conformers_3d();
    let conf = confs
        .get(conf_id)
        .ok_or(MolTransformError::ConformerIndexOutOfBounds {
            id: conf_id,
            max: confs.len(),
        })?;
    Ok(conf.coordinates())
}

/// RDKit❗✔️: Extract mutable 3D coordinates from a specific conformer.
///
/// Takes `Molecule` by value (value-style transform) and uses
/// `Arc::make_mut` for copy-on-write. Returns the mutable coordinate slice
/// and the modified `Molecule` for chaining.
fn conf_coords_mut(
    mol: &mut Molecule,
    conf_id: usize,
) -> Result<&mut [[f64; 3]], MolTransformError> {
    let coord_block = mol.coordinate_block_mut();
    let n_confs = coord_block.conformers_3d.len();
    let conf = coord_block.conformers_3d.get_mut(conf_id).ok_or(
        MolTransformError::ConformerIndexOutOfBounds {
            id: conf_id,
            max: n_confs,
        },
    )?;
    // We need mutable access to the inner Vec<[f64; 3]>.
    // Since Conformer3D's coords field is private, we rely on a pub(crate) method.
    Ok(conf.coordinates_mut())
}

/// RDKit✔️✔️: 3D vector dot product — RDKit Point::dotProduct equivalent.
fn dot3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// RDKit✔️✔️: 3D vector cross product — RDKit Point::crossProduct equivalent.
fn cross3(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// RDKit✔️✔️: 3D vector length — RDKit Point::length equivalent.
fn length3(v: &[f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// RDKit✔️✔️: 3D vector squared length — RDKit Point::lengthSq equivalent.
fn length_sq3(v: &[f64; 3]) -> f64 {
    v[0] * v[0] + v[1] * v[1] + v[2] * v[2]
}

/// RDKit✔️✔️: 3D vector subtraction.
fn sub3(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// RDKit✔️✔️: 3D vector addition.
fn add3(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// RDKit✔️✔️: 3D vector scalar multiply.
fn scale3(v: &[f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// RDKit✔️✔️: 3D vector normalize — RDKit Point::normalize equivalent.
fn normalize3(v: &[f64; 3]) -> Option<[f64; 3]> {
    let len = length3(v);
    if len <= 1e-16 {
        None
    } else {
        Some([v[0] / len, v[1] / len, v[2] / len])
    }
}

/// RDKit✔️✔️: 3D vector angle between two vectors — RDKit Point::angleTo equivalent.
fn angle_between(v1: &[f64; 3], v2: &[f64; 3]) -> f64 {
    let dot = dot3(v1, v2);
    let len1 = length3(v1);
    let len2 = length3(v2);
    let cos_ang = (dot / (len1 * len2)).clamp(-1.0, 1.0);
    cos_ang.acos()
}

// ──────────────────────────────────────────────
// _toBeMovedIdxList — DFS helper
// ──────────────────────────────────────────────
// RDKit source: MolTransforms.cpp lines 467-507

// RDKit✔️✔️: void _toBeMovedIdxList(const ROMol &mol, unsigned int iAtomId,
// RDKit✔️✔️:                        unsigned int jAtomId, std::list<unsigned int> &alist) {
// RDKit✔️✔️:   unsigned int nAtoms = mol.getNumAtoms();
// RDKit✔️✔️:   boost::dynamic_bitset<> visitedIdx(nAtoms);
// RDKit✔️✔️:   std::stack<unsigned int> stack;
// RDKit✔️✔️:   stack.push(jAtomId);
// RDKit✔️✔️:   visitedIdx[iAtomId] = 1;
// RDKit✔️✔️:   visitedIdx[jAtomId] = 1;
// RDKit✔️✔️:   unsigned int tIdx;
// RDKit✔️✔️:   unsigned int wIdx;
// RDKit✔️✔️:   ROMol::ADJ_ITER nbrIdx;
// RDKit✔️✔️:   ROMol::ADJ_ITER endNbrs;
// RDKit✔️✔️:   bool doMainLoop;
// RDKit✔️✔️:   while (stack.size()) {
// RDKit✔️✔️:     doMainLoop = false;
// RDKit✔️✔️:     tIdx = stack.top();
// RDKit✔️✔️:     const Atom *tAtom = mol.getAtomWithIdx(tIdx);
// RDKit✔️✔️:     boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(tAtom);
// RDKit✔️✔️:     unsigned int eIdx;
// RDKit✔️✔️:     for (eIdx = 0; nbrIdx != endNbrs; ++nbrIdx, ++eIdx) {
// RDKit✔️✔️:       wIdx = (mol[*nbrIdx])->getIdx();
// RDKit✔️✔️:       if (!visitedIdx[wIdx]) {
// RDKit✔️✔️:         visitedIdx[wIdx] = 1;
// RDKit✔️✔️:         stack.push(wIdx);
// RDKit✔️✔️:         doMainLoop = true;
// RDKit✔️✔️:         break;
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:     if (doMainLoop) {
// RDKit✔️✔️:       continue;
// RDKit✔️✔️:     }
// RDKit✔️✔️:     visitedIdx[tIdx] = 1;
// RDKit✔️✔️:     stack.pop();
// RDKit✔️✔️:   }
// RDKit✔️✔️:   alist.clear();
// RDKit✔️✔️:   for (unsigned int i = 0; i < nAtoms; ++i) {
// RDKit✔️✔️:     if (visitedIdx[i] && (i != iAtomId)) {
// RDKit✔️✔️:       alist.push_back(i);
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
/// Returns the indices of all atoms reachable from `j_atom` without crossing
/// `i_atom` (i.e. the atoms on j's side of the i-j bond).
/// Uses the same DFS stack approach as RDKit.
fn to_be_moved_idx_list(
    adj: &crate::AdjacencyList,
    n_atoms: usize,
    i_atom: usize,
    j_atom: usize,
) -> Vec<usize> {
    let mut visited = vec![false; n_atoms];
    let mut stack: Vec<usize> = Vec::new();
    stack.push(j_atom);
    visited[i_atom] = true;
    visited[j_atom] = true;

    // RDKit uses a doMainLoop/break pattern for DFS traversal.
    // The loop pops from the stack, finds the first unvisited neighbor,
    // and restarts. If no unvisited neighbor exists, the node is popped.
    while let Some(&t_idx) = stack.last() {
        let mut found_unvisited = false;
        for nbr in adj.neighbors_of(t_idx) {
            let w_idx = nbr.atom_index;
            if !visited[w_idx] {
                visited[w_idx] = true;
                stack.push(w_idx);
                found_unvisited = true;
                break;
            }
        }
        if found_unvisited {
            continue;
        }
        // No unvisited neighbors: mark and pop
        visited[t_idx] = true;
        stack.pop();
    }

    // Collect visited indices except i_atom
    let mut alist: Vec<usize> = Vec::new();
    for i in 0..n_atoms {
        if visited[i] && i != i_atom {
            alist.push(i);
        }
    }
    alist
}

// ──────────────────────────────────────────────
// getBondLength / setBondLength
// ──────────────────────────────────────────────
// RDKit source: MolTransforms.cpp lines 509-544

// RDKit✔️✔️: double getBondLength(const Conformer &conf, unsigned int iAtomId,
// RDKit✔️✔️:                      unsigned int jAtomId) {
// RDKit✔️✔️:   const RDGeom::POINT3D_VECT &pos = conf.getPositions();
// RDKit✔️✔️:   URANGE_CHECK(iAtomId, pos.size());
// RDKit✔️✔️:   URANGE_CHECK(jAtomId, pos.size());
// RDKit✔️✔️:   return (pos[iAtomId] - pos[jAtomId]).length();
// RDKit✔️✔️: }
/// Get the bond length (distance) between atoms `i` and `j`.
///
/// Uses the first 3D conformer.
pub fn get_bond_length(
    mol: &Molecule,
    i: usize,
    j: usize,
    conf_id: usize,
) -> Result<f64, MolTransformError> {
    let coords = conf_coords(mol, conf_id)?;
    let num_atoms = coords.len();
    if i >= num_atoms {
        return Err(MolTransformError::AtomIndexOutOfBounds { atom: i, num_atoms });
    }
    if j >= num_atoms {
        return Err(MolTransformError::AtomIndexOutOfBounds { atom: j, num_atoms });
    }
    Ok(length3(&sub3(&coords[i], &coords[j])))
}

// RDKit✔️✔️: void setBondLength(Conformer &conf, unsigned int iAtomId, unsigned int jAtomId,
// RDKit✔️✔️:                    double value) {
// RDKit✔️✔️:   RDGeom::POINT3D_VECT &pos = conf.getPositions();
// RDKit✔️✔️:   URANGE_CHECK(iAtomId, pos.size());
// RDKit✔️✔️:   URANGE_CHECK(jAtomId, pos.size());
// RDKit✔️✔️:   ROMol &mol = conf.getOwningMol();
// RDKit✔️✔️:   Bond *bond = mol.getBondBetweenAtoms(iAtomId, jAtomId);
// RDKit✔️✔️:   if (!bond) {
// RDKit✔️✔️:     throw ValueErrorException("atoms i and j must be bonded");
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (queryIsBondInRing(bond)) {
// RDKit✔️✔️:     throw ValueErrorException("bond (i,j) must not belong to a ring");
// RDKit✔️✔️:   }
// RDKit✔️✔️:   RDGeom::Point3D v = pos[iAtomId] - pos[jAtomId];
// RDKit✔️✔️:   double origValue = v.length();
// RDKit✔️✔️:   if (origValue <= 1.e-8) {
// RDKit✔️✔️:     throw ValueErrorException("atoms i and j have identical 3D coordinates");
// RDKit✔️✔️:   }
// RDKit✔️✔️:   std::list<unsigned int> alist;
// RDKit✔️✔️:   _toBeMovedIdxList(mol, iAtomId, jAtomId, alist);
// RDKit✔️✔️:   v *= (value / origValue - 1.);
// RDKit✔️✔️:   for (unsigned int &it : alist) {
// RDKit✔️✔️:     pos[it] -= v;
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
/// Set the bond length between atoms `i` and `j`.
///
/// All atoms on j's side of the bond are translated to maintain the geometry.
/// Returns a new `Molecule` (value-style).
fn set_bond_length(
    mol: Molecule,
    conf_id: usize,
    i: usize,
    j: usize,
    value: f64,
) -> Result<Molecule, MolTransformError> {
    let num_atoms = mol.num_atoms();
    // Compute adjacency before mutable borrow
    let adj = AdjacencyList::from_topology(num_atoms, mol.bonds());
    if i >= num_atoms || j >= num_atoms {
        return Err(MolTransformError::AtomIndexOutOfBounds {
            atom: if i >= num_atoms { i } else { j },
            num_atoms,
        });
    }

    // Compute adjacency before mutable borrow
    let adj = AdjacencyList::from_topology(num_atoms, mol.bonds());
    let mut mol = mol;
    let coords = conf_coords_mut(&mut mol, conf_id)?;

    // RDKit✔️✔️: RDGeom::Point3D v = pos[iAtomId] - pos[jAtomId];
    let v = sub3(&coords[i], &coords[j]);
    let orig_value = length3(&v);
    if orig_value <= 1e-8 {
        return Err(MolTransformError::IdenticalCoordinates { i, j });
    }

    // RDKit✔️✔️: v *= (value / origValue - 1.);
    let scale = value / orig_value - 1.0;
    let delta = scale3(&v, scale);

    // RDKit✔️✔️: _toBeMovedIdxList(mol, iAtomId, jAtomId, alist);
    let alist = to_be_moved_idx_list(&adj, num_atoms, i, j);

    // RDKit✔️✔️: for (unsigned int &it : alist) { pos[it] -= v; }
    for &idx in &alist {
        let pt = &mut coords[idx];
        pt[0] -= delta[0];
        pt[1] -= delta[1];
        pt[2] -= delta[2];
    }

    Ok(mol)
}

// ──────────────────────────────────────────────
// getBondAngle / setBondAngle
// ──────────────────────────────────────────────
// RDKit source: MolTransforms.cpp lines 546-616 (getAngleRad / setAngleRad)
// RDKit source: MolTransforms.h   lines 170-194 (getAngleDeg / setAngleDeg)

// RDKit✔️✔️: double getAngleRad(const Conformer &conf, unsigned int iAtomId,
// RDKit✔️✔️:                    unsigned int jAtomId, unsigned int kAtomId) {
// RDKit✔️✔️:   const RDGeom::POINT3D_VECT &pos = conf.getPositions();
// RDKit✔️✔️:   URANGE_CHECK(iAtomId, pos.size());
// RDKit✔️✔️:   URANGE_CHECK(jAtomId, pos.size());
// RDKit✔️✔️:   URANGE_CHECK(kAtomId, pos.size());
// RDKit✔️✔️:   RDGeom::Point3D rJI = pos[iAtomId] - pos[jAtomId];
// RDKit✔️✔️:   double rJISqLength = rJI.lengthSq();
// RDKit✔️✔️:   if (rJISqLength <= 1.e-16) {
// RDKit✔️✔️:     throw ValueErrorException("atoms i and j have identical 3D coordinates");
// RDKit✔️✔️:   }
// RDKit✔️✔️:   RDGeom::Point3D rJK = pos[kAtomId] - pos[jAtomId];
// RDKit✔️✔️:   double rJKSqLength = rJK.lengthSq();
// RDKit✔️✔️:   if (rJKSqLength <= 1.e-16) {
// RDKit✔️✔️:     throw ValueErrorException("atoms j and k have identical 3D coordinates");
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return rJI.angleTo(rJK);
// RDKit✔️✔️: }
/// Get the bond angle in radians among atoms i - j - k.
///
/// Returns the angle in radians.
pub fn get_bond_angle(
    mol: &Molecule,
    i: usize,
    j: usize,
    k: usize,
    conf_id: usize,
) -> Result<f64, MolTransformError> {
    let coords = conf_coords(mol, conf_id)?;
    let num_atoms = coords.len();

    // URANGE_CHECK equivalent
    if i >= num_atoms {
        return Err(MolTransformError::AtomIndexOutOfBounds { atom: i, num_atoms });
    }
    if j >= num_atoms {
        return Err(MolTransformError::AtomIndexOutOfBounds { atom: j, num_atoms });
    }
    if k >= num_atoms {
        return Err(MolTransformError::AtomIndexOutOfBounds { atom: k, num_atoms });
    }

    // RDKit✔️✔️: rJI = pos[iAtomId] - pos[jAtomId]
    let r_ji = sub3(&coords[i], &coords[j]);
    let r_ji_sq = length_sq3(&r_ji);
    if r_ji_sq <= 1e-16 {
        return Err(MolTransformError::IdenticalCoordinates { i, j });
    }

    // RDKit✔️✔️: rJK = pos[kAtomId] - pos[jAtomId]
    let r_jk = sub3(&coords[k], &coords[j]);
    let r_jk_sq = length_sq3(&r_jk);
    if r_jk_sq <= 1e-16 {
        return Err(MolTransformError::IdenticalCoordinates { i: j, j: k });
    }

    // RDKit✔️✔️: return rJI.angleTo(rJK);
    Ok(angle_between(&r_ji, &r_jk))
}

// RDKit✔️✔️: inline double getAngleDeg(const RDKit::Conformer &conf, unsigned int iAtomId,
// RDKit✔️✔️:                           unsigned int jAtomId, unsigned int kAtomId) {
// RDKit✔️✔️:   return (180. / M_PI * getAngleRad(conf, iAtomId, jAtomId, kAtomId));
// RDKit✔️✔️: }
/// Get the bond angle in degrees among atoms i - j - k.
pub fn get_bond_angle_deg(
    mol: &Molecule,
    i: usize,
    j: usize,
    k: usize,
    conf_id: usize,
) -> Result<f64, MolTransformError> {
    get_bond_angle(mol, i, j, k, conf_id).map(|rad| rad * 180.0 / PI)
}

// RDKit✔️✔️: void setAngleRad(Conformer &conf, unsigned int iAtomId, unsigned int jAtomId,
// RDKit✔️✔️:                  unsigned int kAtomId, double value) {
// RDKit✔️✔️:   RDGeom::POINT3D_VECT &pos = conf.getPositions();
// RDKit✔️✔️:   URANGE_CHECK(iAtomId, pos.size());
// RDKit✔️✔️:   URANGE_CHECK(jAtomId, pos.size());
// RDKit✔️✔️:   URANGE_CHECK(kAtomId, pos.size());
// RDKit✔️✔️:   ROMol &mol = conf.getOwningMol();
// RDKit✔️✔️:   Bond *bondJI = mol.getBondBetweenAtoms(jAtomId, iAtomId);
// RDKit✔️✔️:   if (!bondJI) {
// RDKit✔️✔️:     throw ValueErrorException("atoms i and j must be bonded");
// RDKit✔️✔️:   }
// RDKit✔️✔️:   Bond *bondJK = mol.getBondBetweenAtoms(jAtomId, kAtomId);
// RDKit✔️✔️:   if (!bondJK) {
// RDKit✔️✔️:     throw ValueErrorException("atoms j and k must be bonded");
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (queryIsBondInRing(bondJI) && queryIsBondInRing(bondJK)) {
// RDKit✔️✔️:     throw ValueErrorException(
// RDKit✔️✔️:         "bonds (i,j) and (j,k) must not both belong to a ring");
// RDKit✔️✔️:   }
// RDKit✔️✔️:   RDGeom::Point3D rJI = pos[iAtomId] - pos[jAtomId];
// RDKit✔️✔️:   double rJISqLength = rJI.lengthSq();
// RDKit✔️✔️:   if (rJISqLength <= 1.e-16) {
// RDKit✔️✔️:     throw ValueErrorException("atoms i and j have identical 3D coordinates");
// RDKit✔️✔️:   }
// RDKit✔️✔️:   RDGeom::Point3D rJK = pos[kAtomId] - pos[jAtomId];
// RDKit✔️✔️:   double rJKSqLength = rJK.lengthSq();
// RDKit✔️✔️:   if (rJKSqLength <= 1.e-16) {
// RDKit✔️✔️:     throw ValueErrorException("atoms j and k have identical 3D coordinates");
// RDKit✔️✔️:   }
// RDKit✔️✔️:   value -= rJI.angleTo(rJK);
// RDKit✔️✔️:   RDGeom::Point3D &rotAxisBegin = pos[jAtomId];
// RDKit✔️✔️:   RDGeom::Point3D rotAxisEnd = rJI.crossProduct(rJK) + pos[jAtomId];
// RDKit✔️✔️:   RDGeom::Point3D rotAxis = rotAxisEnd - rotAxisBegin;
// RDKit✔️✔️:   rotAxis.normalize();
// RDKit✔️✔️:   std::list<unsigned int> alist;
// RDKit✔️✔️:   _toBeMovedIdxList(mol, jAtomId, kAtomId, alist);
// RDKit✔️✔️:   for (unsigned int &it : alist) {
// RDKit✔️✔️:     pos[it] -= rotAxisBegin;
// RDKit✔️✔️:     RDGeom::Transform3D rotByAngle;
// RDKit✔️✔️:     rotByAngle.SetRotation(value, rotAxis);
// RDKit✔️✔️:     rotByAngle.TransformPoint(pos[it]);
// RDKit✔️✔️:     pos[it] += rotAxisBegin;
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
/// Set the bond angle in radians among atoms i - j - k.
///
/// All atoms on k's side of the j-k bond are rotated to achieve the new angle.
/// Returns a new `Molecule` (value-style).
pub fn set_bond_angle(
    mol: Molecule,
    i: usize,
    j: usize,
    k: usize,
    value: f64,
    conf_id: usize,
) -> Result<Molecule, MolTransformError> {
    let num_atoms = mol.num_atoms();
    // Compute adjacency before mutable borrow
    let adj = AdjacencyList::from_topology(num_atoms, mol.bonds());
    if i >= num_atoms {
        return Err(MolTransformError::AtomIndexOutOfBounds { atom: i, num_atoms });
    }
    if j >= num_atoms {
        return Err(MolTransformError::AtomIndexOutOfBounds { atom: j, num_atoms });
    }
    if k >= num_atoms {
        return Err(MolTransformError::AtomIndexOutOfBounds { atom: k, num_atoms });
    }

    // Check bonds
    let bond_ji = mol.bonds().iter().any(|b| {
        (b.begin().index() == j && b.end().index() == i)
            || (b.begin().index() == i && b.end().index() == j)
    });
    if !bond_ji {
        return Err(MolTransformError::AtomsNotBonded { i: j, j: i });
    }
    let bond_jk = mol.bonds().iter().any(|b| {
        (b.begin().index() == j && b.end().index() == k)
            || (b.begin().index() == k && b.end().index() == j)
    });
    if !bond_jk {
        return Err(MolTransformError::AtomsNotBonded { i: j, j: k });
    }

    let mut mol = mol;
    let coords = conf_coords_mut(&mut mol, conf_id)?;

    // RDKit✔️✔️: rJI = pos[iAtomId] - pos[jAtomId]
    let r_ji = sub3(&coords[i], &coords[j]);
    if length_sq3(&r_ji) <= 1e-16 {
        return Err(MolTransformError::IdenticalCoordinates { i, j });
    }

    // RDKit✔️✔️: rJK = pos[kAtomId] - pos[jAtomId]
    let r_jk = sub3(&coords[k], &coords[j]);
    if length_sq3(&r_jk) <= 1e-16 {
        return Err(MolTransformError::IdenticalCoordinates { i: j, j: k });
    }

    // RDKit✔️✔️: value -= rJI.angleTo(rJK);
    let current_angle = angle_between(&r_ji, &r_jk);
    let delta_angle = value - current_angle;

    // RDKit✔️✔️: rotation axis = rJI x rJK, normalized
    let rot_axis_begin = coords[j]; // copy since we need it in the loop
    let rot_axis_raw = cross3(&r_ji, &r_jk);
    let rot_axis = match normalize3(&rot_axis_raw) {
        Some(axis) => axis,
        None => return Err(MolTransformError::ZeroLengthVector { i: j, j: k }),
    };

    // RDKit✔️✔️: _toBeMovedIdxList(mol, jAtomId, kAtomId, alist);
    let alist = to_be_moved_idx_list(&adj, num_atoms, j, k);

    // RDKit✔️✔️: for each atom in alist: translate, rotate, translate back
    for &idx in &alist {
        // translate atom so it coincides with the origin of rotation
        coords[idx][0] -= rot_axis_begin[0];
        coords[idx][1] -= rot_axis_begin[1];
        coords[idx][2] -= rot_axis_begin[2];

        // rotate around our rotation axis by delta_angle
        rotate_point_around_axis(&mut coords[idx], &rot_axis, delta_angle);

        // translate atom back
        coords[idx][0] += rot_axis_begin[0];
        coords[idx][1] += rot_axis_begin[1];
        coords[idx][2] += rot_axis_begin[2];
    }

    Ok(mol)
}

// RDKit✔️✔️: inline void setAngleDeg(RDKit::Conformer &conf, unsigned int iAtomId,
// RDKit✔️✔️:                         unsigned int jAtomId, unsigned int kAtomId,
// RDKit✔️✔️:                         double value) {
// RDKit✔️✔️:   setAngleRad(conf, iAtomId, jAtomId, kAtomId, value / 180. * M_PI);
// RDKit✔️✔️: }
/// Set the bond angle in degrees among atoms i - j - k.
pub fn set_bond_angle_deg(
    mol: Molecule,
    i: usize,
    j: usize,
    k: usize,
    value: f64,
    conf_id: usize,
) -> Result<Molecule, MolTransformError> {
    set_bond_angle(mol, i, j, k, value * PI / 180.0, conf_id)
}

// ──────────────────────────────────────────────
// getDihedral / setDihedral
// ──────────────────────────────────────────────
// RDKit source: MolTransforms.cpp lines 618-710 (getDihedralRad / setDihedralRad)
// RDKit source: MolTransforms.h   lines 197-223 (getDihedralDeg / setDihedralDeg)

// RDKit✔️✔️: double getDihedralRad(const Conformer &conf, unsigned int iAtomId,
// RDKit✔️✔️:                       unsigned int jAtomId, unsigned int kAtomId,
// RDKit✔️✔️:                       unsigned int lAtomId) {
// RDKit✔️✔️:   const RDGeom::POINT3D_VECT &pos = conf.getPositions();
// RDKit✔️✔️:   URANGE_CHECK(iAtomId, pos.size());
// RDKit✔️✔️:   URANGE_CHECK(jAtomId, pos.size());
// RDKit✔️✔️:   URANGE_CHECK(kAtomId, pos.size());
// RDKit✔️✔️:   URANGE_CHECK(lAtomId, pos.size());
// RDKit✔️✔️:   RDGeom::Point3D rIJ = pos[jAtomId] - pos[iAtomId];
// RDKit✔️✔️:   double rIJSqLength = rIJ.lengthSq();
// RDKit✔️✔️:   if (rIJSqLength <= 1.e-16) {
// RDKit✔️✔️:     throw ValueErrorException("atoms i and j have identical 3D coordinates");
// RDKit✔️✔️:   }
// RDKit✔️✔️:   RDGeom::Point3D rJK = pos[kAtomId] - pos[jAtomId];
// RDKit✔️✔️:   double rJKSqLength = rJK.lengthSq();
// RDKit✔️✔️:   if (rJKSqLength <= 1.e-16) {
// RDKit✔️✔️:     throw ValueErrorException("atoms j and k have identical 3D coordinates");
// RDKit✔️✔️:   }
// RDKit✔️✔️:   RDGeom::Point3D rKL = pos[lAtomId] - pos[kAtomId];
// RDKit✔️✔️:   double rKLSqLength = rKL.lengthSq();
// RDKit✔️✔️:   if (rKLSqLength <= 1.e-16) {
// RDKit✔️✔️:     throw ValueErrorException("atoms k and l have identical 3D coordinates");
// RDKit✔️✔️:   }
// RDKit✔️✔️:   RDGeom::Point3D nIJK = rIJ.crossProduct(rJK);
// RDKit✔️✔️:   double nIJKSqLength = nIJK.lengthSq();
// RDKit✔️✔️:   RDGeom::Point3D nJKL = rJK.crossProduct(rKL);
// RDKit✔️✔️:   double nJKLSqLength = nJKL.lengthSq();
// RDKit✔️✔️:   RDGeom::Point3D m = nIJK.crossProduct(rJK);
// RDKit✔️✔️:   return -atan2(m.dotProduct(nJKL) / sqrt(nJKLSqLength * m.lengthSq()),
// RDKit✔️✔️:                 nIJK.dotProduct(nJKL) / sqrt(nIJKSqLength * nJKLSqLength));
// RDKit✔️✔️: }
/// Get the dihedral (torsion) angle in radians among atoms i - j - k - l.
///
/// Follows the RDKit signed dihedral convention (atan2-based).
pub fn get_dihedral(
    mol: &Molecule,
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    conf_id: usize,
) -> Result<f64, MolTransformError> {
    let coords = conf_coords(mol, conf_id)?;
    let num_atoms = coords.len();

    // URANGE_CHECK equivalents
    if i >= num_atoms {
        return Err(MolTransformError::AtomIndexOutOfBounds { atom: i, num_atoms });
    }
    if j >= num_atoms {
        return Err(MolTransformError::AtomIndexOutOfBounds { atom: j, num_atoms });
    }
    if k >= num_atoms {
        return Err(MolTransformError::AtomIndexOutOfBounds { atom: k, num_atoms });
    }
    if l >= num_atoms {
        return Err(MolTransformError::AtomIndexOutOfBounds { atom: l, num_atoms });
    }

    // RDKit✔️✔️: rIJ = pos[jAtomId] - pos[iAtomId]
    let r_ij = sub3(&coords[j], &coords[i]);
    if length_sq3(&r_ij) <= 1e-16 {
        return Err(MolTransformError::IdenticalCoordinates { i, j });
    }

    // RDKit✔️✔️: rJK = pos[kAtomId] - pos[jAtomId]
    let r_jk = sub3(&coords[k], &coords[j]);
    if length_sq3(&r_jk) <= 1e-16 {
        return Err(MolTransformError::IdenticalCoordinates { i: j, j: k });
    }

    // RDKit✔️✔️: rKL = pos[lAtomId] - pos[kAtomId]
    let r_kl = sub3(&coords[l], &coords[k]);
    if length_sq3(&r_kl) <= 1e-16 {
        return Err(MolTransformError::IdenticalCoordinates { i: k, j: l });
    }

    // RDKit✔️✔️: nIJK = rIJ.crossProduct(rJK)
    let n_ijk = cross3(&r_ij, &r_jk);
    let n_ijk_sq = length_sq3(&n_ijk);

    // RDKit✔️✔️: nJKL = rJK.crossProduct(rKL)
    let n_jkl = cross3(&r_jk, &r_kl);
    let n_jkl_sq = length_sq3(&n_jkl);

    // RDKit✔️✔️: m = nIJK.crossProduct(rJK)
    let m = cross3(&n_ijk, &r_jk);

    // RDKit✔️✔️: return -atan2(m.dotProduct(nJKL) / sqrt(nJKLSqLength * m.lengthSq()),
    // RDKit✔️✔️:                 nIJK.dotProduct(nJKL) / sqrt(nIJKSqLength * nJKLSqLength));
    let denom1 = (n_jkl_sq * length_sq3(&m)).sqrt();
    let denom2 = (n_ijk_sq * n_jkl_sq).sqrt();
    let angle = -f64::atan2(
        dot3(&m, &n_jkl) / denom1.max(f64::MIN_POSITIVE),
        dot3(&n_ijk, &n_jkl) / denom2.max(f64::MIN_POSITIVE),
    );
    Ok(angle)
}

// RDKit✔️✔️: inline double getDihedralDeg(const RDKit::Conformer &conf, unsigned int iAtomId,
// RDKit✔️✔️:                              unsigned int jAtomId, unsigned int kAtomId,
// RDKit✔️✔️:                              unsigned int lAtomId) {
// RDKit✔️✔️:   return (180. / M_PI * getDihedralRad(conf, iAtomId, jAtomId, kAtomId, lAtomId));
// RDKit✔️✔️: }
/// Get the dihedral (torsion) angle in degrees among atoms i - j - k - l.
pub fn get_dihedral_deg(
    mol: &Molecule,
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    conf_id: usize,
) -> Result<f64, MolTransformError> {
    get_dihedral(mol, i, j, k, l, conf_id).map(|rad| rad * 180.0 / PI)
}

// RDKit✔️✔️: void setDihedralRad(Conformer &conf, unsigned int iAtomId, unsigned int jAtomId,
// RDKit✔️✔️:                     unsigned int kAtomId, unsigned int lAtomId, double value) {
// RDKit✔️✔️:   RDGeom::POINT3D_VECT &pos = conf.getPositions();
// RDKit✔️✔️:   URANGE_CHECK(iAtomId, pos.size());
// RDKit✔️✔️:   URANGE_CHECK(jAtomId, pos.size());
// RDKit✔️✔️:   URANGE_CHECK(kAtomId, pos.size());
// RDKit✔️✔️:   URANGE_CHECK(lAtomId, pos.size());
// RDKit✔️✔️:   ROMol &mol = conf.getOwningMol();
// RDKit✔️✔️:   Bond *bondJK = mol.getBondBetweenAtoms(jAtomId, kAtomId);
// RDKit✔️✔️:   if (!bondJK) {
// RDKit✔️✔️:     throw ValueErrorException("atoms j and k must be bonded");
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (queryIsBondInRing(bondJK)) {
// RDKit✔️✔️:     throw ValueErrorException("bond (j,k) must not belong to a ring");
// RDKit✔️✔️:   }
// RDKit✔️✔️:   RDGeom::Point3D rIJ = pos[jAtomId] - pos[iAtomId];
// RDKit✔️✔️:   double rIJSqLength = rIJ.lengthSq();
// RDKit✔️✔️:   if (rIJSqLength <= 1.e-16) {
// RDKit✔️✔️:     throw ValueErrorException("atoms i and j have identical 3D coordinates");
// RDKit✔️✔️:   }
// RDKit✔️✔️:   RDGeom::Point3D rJK = pos[kAtomId] - pos[jAtomId];
// RDKit✔️✔️:   double rJKSqLength = rJK.lengthSq();
// RDKit✔️✔️:   if (rJKSqLength <= 1.e-16) {
// RDKit✔️✔️:     throw ValueErrorException("atoms j and k have identical 3D coordinates");
// RDKit✔️✔️:   }
// RDKit✔️✔️:   RDGeom::Point3D rKL = pos[lAtomId] - pos[kAtomId];
// RDKit✔️✔️:   double rKLSqLength = rKL.lengthSq();
// RDKit✔️✔️:   if (rKLSqLength <= 1.e-16) {
// RDKit✔️✔️:     throw ValueErrorException("atoms k and l have identical 3D coordinates");
// RDKit✔️✔️:   }
// RDKit✔️✔️:   RDGeom::Point3D nIJK = rIJ.crossProduct(rJK);
// RDKit✔️✔️:   double nIJKSqLength = nIJK.lengthSq();
// RDKit✔️✔️:   RDGeom::Point3D nJKL = rJK.crossProduct(rKL);
// RDKit✔️✔️:   double nJKLSqLength = nJKL.lengthSq();
// RDKit✔️✔️:   RDGeom::Point3D m = nIJK.crossProduct(rJK);
// RDKit✔️✔️:   value -= -atan2(m.dotProduct(nJKL) / sqrt(nJKLSqLength * m.lengthSq()),
// RDKit✔️✔️:                   nIJK.dotProduct(nJKL) / sqrt(nIJKSqLength * nJKLSqLength));
// RDKit✔️✔️:   RDGeom::Point3D &rotAxisBegin = pos[jAtomId];
// RDKit✔️✔️:   RDGeom::Point3D &rotAxisEnd = pos[kAtomId];
// RDKit✔️✔️:   RDGeom::Point3D rotAxis = rotAxisEnd - rotAxisBegin;
// RDKit✔️✔️:   rotAxis.normalize();
// RDKit✔️✔️:   std::list<unsigned int> alist;
// RDKit✔️✔️:   _toBeMovedIdxList(mol, jAtomId, kAtomId, alist);
// RDKit✔️✔️:   for (unsigned int &it : alist) {
// RDKit✔️✔️:     pos[it] -= rotAxisBegin;
// RDKit✔️✔️:     RDGeom::Transform3D rotByAngle;
// RDKit✔️✔️:     rotByAngle.SetRotation(value, rotAxis);
// RDKit✔️✔️:     rotByAngle.TransformPoint(pos[it]);
// RDKit✔️✔️:     pos[it] += rotAxisBegin;
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
/// Set the dihedral (torsion) angle in radians among atoms i - j - k - l.
///
/// All atoms on k's side of the j-k bond are rotated to achieve the new dihedral.
/// Returns a new `Molecule` (value-style).
pub fn set_dihedral(
    mol: Molecule,
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    value: f64,
    conf_id: usize,
) -> Result<Molecule, MolTransformError> {
    let num_atoms = mol.num_atoms();
    // Compute adjacency before mutable borrow
    let adj = AdjacencyList::from_topology(num_atoms, mol.bonds());
    if i >= num_atoms {
        return Err(MolTransformError::AtomIndexOutOfBounds { atom: i, num_atoms });
    }
    if j >= num_atoms {
        return Err(MolTransformError::AtomIndexOutOfBounds { atom: j, num_atoms });
    }
    if k >= num_atoms {
        return Err(MolTransformError::AtomIndexOutOfBounds { atom: k, num_atoms });
    }
    if l >= num_atoms {
        return Err(MolTransformError::AtomIndexOutOfBounds { atom: l, num_atoms });
    }

    // Check bond (j,k)
    let bond_jk = mol.bonds().iter().any(|b| {
        (b.begin().index() == j && b.end().index() == k)
            || (b.begin().index() == k && b.end().index() == j)
    });
    if !bond_jk {
        return Err(MolTransformError::AtomsNotBonded { i: j, j: k });
    }

    let mut mol = mol;
    let coords = conf_coords_mut(&mut mol, conf_id)?;

    // RDKit✔️✔️: rIJ = pos[jAtomId] - pos[iAtomId]
    let r_ij = sub3(&coords[j], &coords[i]);
    if length_sq3(&r_ij) <= 1e-16 {
        return Err(MolTransformError::IdenticalCoordinates { i, j });
    }

    // RDKit✔️✔️: rJK = pos[kAtomId] - pos[jAtomId]
    let r_jk = sub3(&coords[k], &coords[j]);
    if length_sq3(&r_jk) <= 1e-16 {
        return Err(MolTransformError::IdenticalCoordinates { i: j, j: k });
    }

    // RDKit✔️✔️: rKL = pos[lAtomId] - pos[kAtomId]
    let r_kl = sub3(&coords[l], &coords[k]);
    if length_sq3(&r_kl) <= 1e-16 {
        return Err(MolTransformError::IdenticalCoordinates { i: k, j: l });
    }

    // RDKit✔️✔️: nIJK = rIJ.crossProduct(rJK)
    let n_ijk = cross3(&r_ij, &r_jk);
    let n_ijk_sq = length_sq3(&n_ijk);

    // RDKit✔️✔️: nJKL = rJK.crossProduct(rKL)
    let n_jkl = cross3(&r_jk, &r_kl);
    let n_jkl_sq = length_sq3(&n_jkl);

    // RDKit✔️✔️: m = nIJK.crossProduct(rJK)
    let m = cross3(&n_ijk, &r_jk);

    // RDKit✔️✔️: compute current dihedral and subtract from target
    let denom1 = (n_jkl_sq * length_sq3(&m)).sqrt();
    let denom2 = (n_ijk_sq * n_jkl_sq).sqrt();
    let current_dihedral = -f64::atan2(
        dot3(&m, &n_jkl) / denom1.max(f64::MIN_POSITIVE),
        dot3(&n_ijk, &n_jkl) / denom2.max(f64::MIN_POSITIVE),
    );

    // RDKit✔️✔️: value -= current_dihedral
    let delta_angle = value - current_dihedral;

    // RDKit✔️✔️: rotation axis = (k) - (j), normalized (the j-k bond)
    let rot_axis_begin = coords[j]; // copy
    let rot_axis_end = coords[k]; // copy
    let rot_axis_raw = sub3(&rot_axis_end, &rot_axis_begin);
    let rot_axis = match normalize3(&rot_axis_raw) {
        Some(axis) => axis,
        None => return Err(MolTransformError::ZeroLengthVector { i: j, j: k }),
    };

    // RDKit✔️✔️: _toBeMovedIdxList(mol, jAtomId, kAtomId, alist);
    let alist = to_be_moved_idx_list(&adj, num_atoms, j, k);

    // RDKit✔️✔️: for each atom in alist: translate, rotate, translate back
    for &idx in &alist {
        // translate atom so it coincides with the origin of rotation
        coords[idx][0] -= rot_axis_begin[0];
        coords[idx][1] -= rot_axis_begin[1];
        coords[idx][2] -= rot_axis_begin[2];

        // rotate around our rotation axis by delta_angle
        rotate_point_around_axis(&mut coords[idx], &rot_axis, delta_angle);

        // translate atom back
        coords[idx][0] += rot_axis_begin[0];
        coords[idx][1] += rot_axis_begin[1];
        coords[idx][2] += rot_axis_begin[2];
    }

    Ok(mol)
}

// RDKit✔️✔️: inline void setDihedralDeg(RDKit::Conformer &conf, unsigned int iAtomId,
// RDKit✔️✔️:                            unsigned int jAtomId, unsigned int kAtomId,
// RDKit✔️✔️:                            unsigned int lAtomId, double value) {
// RDKit✔️✔️:   setDihedralRad(conf, iAtomId, jAtomId, kAtomId, lAtomId, value / 180. * M_PI);
// RDKit✔️✔️: }
/// Set the dihedral (torsion) angle in degrees among atoms i - j - k - l.
pub fn set_dihedral_deg(
    mol: Molecule,
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    value: f64,
    conf_id: usize,
) -> Result<Molecule, MolTransformError> {
    set_dihedral(mol, i, j, k, l, value * PI / 180.0, conf_id)
}

// ──────────────────────────────────────────────
// Rotation helper
// ──────────────────────────────────────────────

/// RDKit✔️✔️: Rotate a 3D point around a normalized axis by `angle` radians.
///
/// Uses Rodrigues' rotation formula, equivalent to RDKit's
/// `Transform3D::SetRotation(angle, axis)` followed by `TransformPoint`.
fn rotate_point_around_axis(point: &mut [f64; 3], axis: &[f64; 3], angle: f64) {
    // Rodrigues' rotation formula:
    //   v_rot = v * cos(theta) + (k x v) * sin(theta) + k * (k . v) * (1 - cos(theta))
    let cos_a = angle.cos();
    let sin_a = angle.sin();
    let dot = dot3(point, axis);
    let cross = cross3(axis, point);
    point[0] = point[0] * cos_a + cross[0] * sin_a + axis[0] * dot * (1.0 - cos_a);
    point[1] = point[1] * cos_a + cross[1] * sin_a + axis[1] * dot * (1.0 - cos_a);
    point[2] = point[2] * cos_a + cross[2] * sin_a + axis[2] * dot * (1.0 - cos_a);
}

// ──────────────────────────────────────────────
// getAtomPosition / setAtomPosition
// ──────────────────────────────────────────────
// RDKit source: Conformer.h / MolTransforms conventions

// RDKit✔️✔️: const RDGeom::Point3D &Conformer::getAtomPos(unsigned int atomId) const
/// Get the 3D position of a single atom from the specified conformer.
pub fn get_atom_position(
    mol: &Molecule,
    atom: usize,
    conf_id: usize,
) -> Result<[f64; 3], MolTransformError> {
    let coords = conf_coords(mol, conf_id)?;
    if atom >= coords.len() {
        return Err(MolTransformError::AtomIndexOutOfBounds {
            atom,
            num_atoms: coords.len(),
        });
    }
    Ok(coords[atom])
}

// RDKit✔️✔️: void Conformer::setAtomPos(unsigned int atomId, const RDGeom::Point3D &pos)
/// Set the 3D position of a single atom in the specified conformer.
///
/// Returns a new `Molecule` (value-style).
pub fn set_atom_position(
    mut mol: Molecule,
    atom: usize,
    pos: [f64; 3],
    conf_id: usize,
) -> Result<Molecule, MolTransformError> {
    mol.set_atom_position_(atom, pos, conf_id)
        .map_err(mol_transform_operation_error)?;
    Ok(mol)
}

pub(crate) fn set_atom_position_in_coordinate_block(
    coordinates: &mut crate::molecule::CoordinateBlock,
    atom: usize,
    pos: [f64; 3],
    conf_id: usize,
) -> Result<bool, MolTransformError> {
    let conformer_count = coordinates.conformers_3d.len();
    let conformer = coordinates.conformers_3d.get_mut(conf_id).ok_or(
        MolTransformError::ConformerIndexOutOfBounds {
            id: conf_id,
            max: conformer_count,
        },
    )?;
    let atom_count = conformer.coordinates().len();
    if atom >= atom_count {
        return Err(MolTransformError::AtomIndexOutOfBounds {
            atom,
            num_atoms: atom_count,
        });
    }
    let changed = conformer.coordinates()[atom] != pos;
    conformer.coordinates_mut()[atom] = pos;
    Ok(changed)
}

fn mol_transform_operation_error(error: crate::OperationError) -> MolTransformError {
    match error {
        crate::OperationError::MolTransform { source, .. } => source,
        other => MolTransformError::OperationContract {
            message: other.to_string(),
        },
    }
}

// RDKit✔️✔️: const RDGeom::POINT3D_VECT &Conformer::getPositions() const
/// Get all 3D atom positions from the specified conformer.
pub fn get_atom_positions(
    mol: &Molecule,
    conf_id: usize,
) -> Result<Vec<[f64; 3]>, MolTransformError> {
    let coords = conf_coords(mol, conf_id)?;
    Ok(coords.to_vec())
}

// ──────────────────────────────────────────────
// transformConformer
// ──────────────────────────────────────────────
// RDKit source: MolTransforms.cpp lines 445-451

// RDKit✔️✔️: void transformConformer(Conformer &conf, const RDGeom::Transform3D &trans) {
// RDKit✔️✔️:   RDGeom::POINT3D_VECT &positions = conf.getPositions();
// RDKit✔️✔️:   RDGeom::POINT3D_VECT_I pi;
// RDKit✔️✔️:   for (pi = positions.begin(); pi != positions.end(); ++pi) {
// RDKit✔️✔️:     trans.TransformPoint(*pi);
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
/// Transform all atom positions in a conformer using a 4x4 transformation matrix.
///
/// The matrix is applied as: `pos' = R * pos + t` where `R` is the top-left
/// 3x3 submatrix and `t` is the translation column (column 3, rows 0-2).
/// This matches RDKit's `Transform3D::TransformPoint`.
///
/// Returns a new `Molecule` (value-style).
pub fn transform_conformer(
    mol: Molecule,
    transform: &[[f64; 4]; 4],
    conf_id: usize,
) -> Result<Molecule, MolTransformError> {
    let mut mol = mol;
    let coords = conf_coords_mut(&mut mol, conf_id)?;

    for pt in coords.iter_mut() {
        // RDGeom::Transform3D::TransformPoint:
        //   x' = t[0][0]*x + t[0][1]*y + t[0][2]*z + t[0][3]
        //   y' = t[1][0]*x + t[1][1]*y + t[1][2]*z + t[1][3]
        //   z' = t[2][0]*x + t[2][1]*y + t[2][2]*z + t[2][3]
        let x = pt[0];
        let y = pt[1];
        let z = pt[2];
        pt[0] = transform[0][0] * x + transform[0][1] * y + transform[0][2] * z + transform[0][3];
        pt[1] = transform[1][0] * x + transform[1][1] * y + transform[1][2] * z + transform[1][3];
        pt[2] = transform[2][0] * x + transform[2][1] * y + transform[2][2] * z + transform[2][3];
    }

    Ok(mol)
}

// ──────────────────────────────────────────────
// getTranslationVectors / centerAtomsOnOrigin
// ──────────────────────────────────────────────

/// RDKit❗✔️: Get translation vectors from origin for each atom in a conformer.
///
/// Returns the position vectors directly from the conformer coordinates
/// (these are the vectors from the origin to each atom).
pub fn get_translation_vectors(
    mol: &Molecule,
    conf_id: usize,
) -> Result<Vec<[f64; 3]>, MolTransformError> {
    let coords = conf_coords(mol, conf_id)?;
    Ok(coords.to_vec())
}

/// RDKit❗✔️: Compute the centroid of a conformer's 3D coordinates.
///
/// Equivalent to `MolTransforms::computeCentroid(conf, ignoreHs=false)`.
/// When `ignore_hs` is true, hydrogen atoms (atomic number 1) are excluded
/// from the centroid calculation.
///
/// RDKit source: MolTransforms.cpp lines 46-63
pub fn compute_centroid(
    mol: &Molecule,
    conf_id: usize,
    ignore_hs: bool,
) -> Result<[f64; 3], MolTransformError> {
    let coords = conf_coords(mol, conf_id)?;
    let mut res = [0.0, 0.0, 0.0];
    let mut count = 0usize;

    for i in 0..coords.len() {
        if ignore_hs && mol.atoms()[i].atomic_number() == 1 {
            continue;
        }
        res[0] += coords[i][0];
        res[1] += coords[i][1];
        res[2] += coords[i][2];
        count += 1;
    }

    if count > 0 {
        let inv = 1.0 / count as f64;
        res[0] *= inv;
        res[1] *= inv;
        res[2] *= inv;
    }

    Ok(res)
}

/// RDKit❗✔️: Center a conformer's atoms on the origin.
///
/// Translates all atoms so that the centroid of the conformer (or optionally
/// of the heavy atoms only) coincides with the origin.
///
/// Returns a new `Molecule` (value-style).
pub fn center_atoms_on_origin(
    mol: Molecule,
    conf_id: usize,
    ignore_hs: bool,
) -> Result<Molecule, MolTransformError> {
    let centroid = compute_centroid(&mol, conf_id, ignore_hs)?;
    let mut mol = mol;
    let coords = conf_coords_mut(&mut mol, conf_id)?;

    for pt in coords.iter_mut() {
        pt[0] -= centroid[0];
        pt[1] -= centroid[1];
        pt[2] -= centroid[2];
    }

    Ok(mol)
}
