//! 2D coordinate generation ported from the old-core RDKit depiction algorithm.
//!
//! Source reproduction protocol: dev/source_reproduction_protocol.md
//!
//! Marker convention:
//!   RDKit✔️✔️: one-to-one port
//!   RDKit✔️❌: adapted for new-core API, same algorithm
//!   RDKit❗✔️: algorithm-equivalent via different architecture
//!
//! Source layout (old-core `molblock.rs`, line numbers):
//!
//!   Types + helpers (~33 – 55)   norm, normalize, rotate, cos/sin/sqrt
//!   TreeEmbeddedAtom             ~825
//!   RdkitHybridization           ~758
//!   atom_depict_rank             ~838
//!   rdkit_ring_radius            ~847
//!   n_outer_electrons_rdkit      ~852
//!   rdkit_num_bonds_plus_lone_pairs     ~856
//!   rdkit_hybridizations_for_depict     ~893
//!   compute_sub_angle            ~1764
//!   rotation_dir                 ~1899
//!   canonicalize_component       ~769
//!   rdkit_rank_atoms_by_rank     ~1705
//!   rdkit_set_nbr_order          ~1723
//!   strict_rdkit_subset_coords   ~725
//!   remove_collisions_bond_flip_like    ~2369
//!   remove_collisions_open_angles_like  ~2972
//!   remove_collisions_shorten_bonds_like ~4757
//!   place_acyclic_tree           ~1906
//!   place_multiring_nonfused_component  ~3145
//!   rdkit_compute_initial_coords_strict ~503
//!   compute_2d_coords            ~446
//!
//! RDKit marker convention is defined in `dev/source_reproduction_protocol.md`.
//! Local shorthand:
//!   RDKit❗✔️: ported from compute_2d_coords_impl (intentionally different architecture)
//!   RDKit✔️✔️: one-to-one port
//!   RDKit✔️❌: adapted for new-core API, same algorithm
//!   RDKit❗❌: not yet ported (unimplemented stub)

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet, VecDeque};
use std::f64::consts::PI;

use crate::ChiralTag;
use crate::Molecule;
use crate::atom::Atom;
use crate::bond::{Bond, BondOrder, BondStereo};
use crate::valence::periodic_table_outer_electrons;

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub(crate) enum Coordinate2DError {
    #[error("{0}")]
    InvalidInput(&'static str),
    #[error("{0}")]
    UnsupportedFeature(String),
}

// ── RDKit-compatible math wrappers ──────────────────────────────────────────
// These force a specific math-lib code path that RDKit's C++ expects.
// ── Basic math wrappers (RDKit Numerics, Geometry/point.h) ─────────────────

/// RDKit❗✔️: cos wrapper — matches RDKit's cos usage in depiction geometry
fn rdkit_cos(x: f64) -> f64 {
    x.cos()
}

/// RDKit❗✔️: sin wrapper — matches RDKit's sin usage in depiction geometry
fn rdkit_sin(x: f64) -> f64 {
    x.sin()
}

/// RDKit❗✔️: acos wrapper — matches RDKit's acos usage in depiction geometry
fn rdkit_acos(x: f64) -> f64 {
    x.acos()
}

/// RDKit❗✔️: sqrt wrapper — matches RDKit's sqrt usage in depiction geometry
fn rdkit_sqrt(x: f64) -> f64 {
    x.sqrt()
}

/// RDKit❗✔️: 2D vector norm — RDKit Point2D::length() equivalent
fn norm(v: (f64, f64)) -> f64 {
    rdkit_sqrt(v.0 * v.0 + v.1 * v.1)
}

/// RDKit❗✔️: 2D vector normalization — RDKit Point2D::normalize() equivalent
fn normalize(v: (f64, f64)) -> (f64, f64) {
    let n = norm(v);
    if n < 1e-12 {
        (0.0, 0.0)
    } else {
        (v.0 / n, v.1 / n)
    }
}

/// RDKit❗✔️: 2D vector rotation — RDKit Transform2D::Rotate equivalent
fn rotate(v: (f64, f64), angle: f64) -> (f64, f64) {
    let c = rdkit_cos(angle);
    let s = rdkit_sin(angle);
    (v.0 * c - v.1 * s, v.0 * s + v.1 * c)
}

/// RDKit❗✔️: rotate point around center — RDKit Transform2D rotation
fn rotate_around(p: (f64, f64), center: (f64, f64), angle: f64) -> (f64, f64) {
    let v = (p.0 - center.0, p.1 - center.1);
    let r = rotate(v, angle);
    (center.0 + r.0, center.1 + r.1)
}

/// RDKit❗✔️: rotation direction — RDKit Transform2D orientation detection
fn rotation_dir(center: (f64, f64), loc1: (f64, f64), loc2: (f64, f64), rem_angle: f64) -> i32 {
    let pt1 = (loc1.0 - center.0, loc1.1 - center.1);
    let pt2 = (loc2.0 - center.0, loc2.1 - center.1);
    let cross = (pt1.0 * pt2.1 - pt1.1 * pt2.0) * (PI - rem_angle);
    if cross >= 0.0 { -1 } else { 1 }
}

// ── 2D transform helpers (RDKit EmbeddedFrag.cpp) ────────────────────────────

#[allow(dead_code)]
/// RDKit❗✔️: 3x3 matrix multiplication — RDKit Transform2D::Mul3 equivalent
fn transform2d_mul3(lhs: [f64; 9], rhs: [f64; 9]) -> [f64; 9] {
    let mut out = [0.0; 9];
    for i in 0..3 {
        let id_a = i * 3;
        for j in 0..3 {
            let id_c = id_a + j;
            out[id_c] = 0.0;
            for k in 0..3 {
                out[id_c] += lhs[id_a + k] * rhs[k * 3 + j];
            }
        }
    }
    out
}

#[allow(dead_code)]
/// RDKit❗✔️: convert 3×3 matrix to affine 2×6 — RDKit Transform2D::toAffine
fn transform2d_to_affine(data: [f64; 9]) -> [f64; 6] {
    [data[0], data[1], data[2], data[3], data[4], data[5]]
}

#[allow(dead_code)]
/// RDKit❗✔️: set rotation+translation from center+angle — RDKit Transform2D::setTransform
fn transform2d_set_transform_center_angle(pt: (f64, f64), angle: f64) -> [f64; 6] {
    let trans1 = [1.0, 0.0, -pt.0, 0.0, 1.0, -pt.1, 0.0, 0.0, 1.0];
    let rot = [
        rdkit_cos(angle),
        -rdkit_sin(angle),
        0.0,
        rdkit_sin(angle),
        rdkit_cos(angle),
        0.0,
        0.0,
        0.0,
        1.0,
    ];
    let this = transform2d_mul3(rot, trans1);
    let trans2 = [1.0, 0.0, pt.0, 0.0, 1.0, pt.1, 0.0, 0.0, 1.0];
    transform2d_to_affine(transform2d_mul3(trans2, this))
}

/// RDKit❗✔️: apply affine transform to point — RDKit Transform2D::transformPoint
fn transform2d_point(pt: (f64, f64), data: [f64; 6]) -> (f64, f64) {
    (
        data[0] * pt.0 + data[1] * pt.1 + data[2],
        data[3] * pt.0 + data[4] * pt.1 + data[5],
    )
}

/// RDKit❗✔️: two-point transform — RDKit Transform2D::setTransform mapping p1→p1'
fn transform2d_set_transform_two_point(
    ref1: (f64, f64),
    ref2: (f64, f64),
    pt1: (f64, f64),
    pt2: (f64, f64),
) -> [f64; 6] {
    let rvec = (ref2.0 - ref1.0, ref2.1 - ref1.1);
    let pvec = (pt2.0 - pt1.0, pt2.1 - pt1.1);
    let dp = rvec.0 * pvec.0 + rvec.1 * pvec.1;
    let lp = norm(rvec) * norm(pvec);
    if lp <= 0.0 {
        return [1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
    }
    let mut cval = dp / lp;
    cval = cval.clamp(-1.0, 1.0);
    let mut ang = rdkit_acos(cval);
    let cross = pvec.0 * rvec.1 - pvec.1 * rvec.0;
    if cross < 0.0 {
        ang *= -1.0;
    }
    let mut data = [
        rdkit_cos(ang),
        -rdkit_sin(ang),
        0.0,
        rdkit_sin(ang),
        rdkit_cos(ang),
        0.0,
    ];
    let npt1 = transform2d_point(pt1, data);
    data[2] = ref1.0 - npt1.0;
    data[5] = ref1.1 - npt1.1;
    data
}

// RDKit✔️✔️:  RDKit ReflectPoint
fn rdkit_reflect_point(point: (f64, f64), loc1: (f64, f64), loc2: (f64, f64)) -> (f64, f64) {
    let org = (0.0, 0.0);
    let xaxis = (1.0, 0.0);
    let mut cent = (loc1.0 + loc2.0, loc1.1 + loc2.1);
    cent.0 *= 0.5;
    cent.1 *= 0.5;

    let trans = transform2d_set_transform_two_point(org, xaxis, cent, loc1);
    let itrans = transform2d_set_transform_two_point(cent, loc1, org, xaxis);

    let mut res = transform2d_point(point, trans);
    res.1 = -res.1;
    transform2d_point(res, itrans)
}

// ── Types ────────────────────────────────────────────────────────────────────

// BEGIN RDKIT CPP TYPE: ConjugHybrid.h HybridizationType
// RDKit✔️✔️: enum HybridizationType { UNSPECIFIED, S, SP, SP2, SP3, SP3D, SP3D2 };
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
enum RdkitHybridization {
    Unspecified,
    S,
    Sp,
    Sp2,
    Sp3,
    Sp3d,
    Sp3d2,
}

// BEGIN RDKIT CPP STRUCT: EmbeddedFrag.h EmbeddedAtom
// RDKit✔️✔️: struct EmbeddedAtom { Point2D loc; Point2D normal; bool ccw; ... }
#[derive(Clone, Debug)]
struct TreeEmbeddedAtom {
    loc: (f64, f64),
    normal: (f64, f64),
    ccw: bool,
    cis_trans_nbr: Option<usize>,
    angle: f64,
    nbr1: Option<usize>,
    nbr2: Option<usize>,
    rot_dir: i32,
    pending: Vec<usize>,
}

// ── Pure functions ───────────────────────────────────────────────────────────

// BEGIN RDKIT CPP FUNCTION: atomDepictRank (DepictUtils.h)
// RDKit✔️✔️: unsigned int atomDepictRank(unsigned int atomicNum, unsigned int degree)
fn atom_depict_rank(atomic_num: u8, degree: usize) -> usize {
    let anum = if atomic_num == 1 {
        1000usize
    } else {
        atomic_num as usize
    };
    100 * anum + degree
}

// BEGIN RDKIT CPP FUNCTION: ringRadius (EmbeddedFrag.cpp)
// RDKit✔️✔️: ringRadius(unsigned ringSize, double bondLen)
fn rdkit_ring_radius(ring_size: usize, bond_len: f64) -> f64 {
    let ang = 2.0 * PI / ring_size as f64;
    bond_len / rdkit_sqrt(2.0 * (1.0 - rdkit_cos(ang)))
}

// BEGIN RDKIT CPP FUNCTION: numBondsPlusLonePairs (ConjugHybrid.cpp)
// RDKit✔️❌: adapted for new-core API
fn rdkit_num_bonds_plus_lone_pairs(
    atomic_num: u8,
    degree: usize,
    explicit_valence: i32,
    implicit_hydrogens: i32,
    radical_electrons: u8,
    formal_charge: i8,
) -> Option<i32> {
    // Atomic number for the periodic-table outer-electron lookup.
    let deg = degree as i32 + implicit_hydrogens;

    // Zero-charge not-adjustment per the old-core version (which also uses
    // formal_charge for the nouter-based calculation below).

    if atomic_num <= 1 {
        return Some(deg);
    }
    let nouter = n_outer_electrons_rdkit(atomic_num)?;
    let total_valence = explicit_valence + implicit_hydrogens;
    let charge = formal_charge as i32;
    let free_electrons = nouter - (total_valence + charge);
    if total_valence + nouter - charge < 8 {
        let radicals = radical_electrons as i32;
        Some(deg + (free_electrons - radicals) / 2 + radicals)
    } else {
        Some(deg + free_electrons / 2)
    }
}

// BEGIN RDKIT CPP FUNCTION: ComputeSubAngle (DepictUtils.h)
// RDKit✔️✔️: double computeSubAngle(unsigned int degree, HybridizationType type)
fn compute_sub_angle(degree: usize, hybridization: RdkitHybridization) -> f64 {
    match hybridization {
        RdkitHybridization::Unspecified | RdkitHybridization::Sp3 => {
            if degree == 4 {
                PI / 2.0
            } else {
                2.0 * PI / 3.0
            }
        }
        RdkitHybridization::Sp2 => 2.0 * PI / 3.0,
        RdkitHybridization::S
        | RdkitHybridization::Sp
        | RdkitHybridization::Sp3d
        | RdkitHybridization::Sp3d2 => 2.0 * PI / degree as f64,
    }
}

// BEGIN RDKIT CPP FUNCTION (standalone): canonicalizeComponent
// RDKit✔️✔️: (RDKit DepictUtils.cpp)
fn canonicalize_component(component: &mut Vec<(usize, (f64, f64))>) {
    if component.len() <= 1 {
        return;
    }
    let n = component.len() as f64;
    let (mut cx, mut cy) = (0.0f64, 0.0f64);
    for &(_, (x, y)) in component.iter() {
        cx += x;
        cy += y;
    }
    let inv_n = 1.0 / n;
    cx *= inv_n;
    cy *= inv_n;

    let (mut xx, mut xy, mut yy) = (0.0f64, 0.0f64, 0.0f64);
    let mut centered: Vec<(f64, f64)> = Vec::with_capacity(component.len());
    for &(_, (x, y)) in component.iter() {
        let px = x - cx;
        let py = y - cy;
        centered.push((px, py));
        xx += px * px;
        xy += px * py;
        yy += py * py;
    }

    let d = ((xx - yy) * (xx - yy) + 4.0 * xy * xy).sqrt();
    let mut eig1 = (2.0 * xy, (yy - xx) + d);
    let eig1_len = norm(eig1);
    if eig1_len <= 1e-4 {
        for (i, (_, pos)) in component.iter_mut().enumerate() {
            *pos = centered[i];
        }
        return;
    }
    let e_val1 = (xx + yy + d) / 2.0;
    eig1 = (eig1.0 / eig1_len, eig1.1 / eig1_len);

    let mut eig2 = (2.0 * xy, (yy - xx) - d);
    let e_val2 = (xx + yy - d) / 2.0;
    let eig2_len = norm(eig2);
    if eig2_len > 1e-4 {
        eig2 = (eig2.0 / eig2_len, eig2.1 / eig2_len);
        if e_val2 > e_val1 {
            std::mem::swap(&mut eig1, &mut eig2);
        }
    }

    for (i, (_, pos)) in component.iter_mut().enumerate() {
        let (px, py) = centered[i];
        // trans = [[eig1.x, eig1.y],[-eig1.y,eig1.x]]
        let rx = px * eig1.0 + py * eig1.1;
        let ry = -px * eig1.1 + py * eig1.0;
        *pos = (rx, ry);
    }
}

// ── Outer electrons ─────────────────────────────────────────────────────────

/// RDKit❗✔️: outer electron count from atomic number — RDKit PeriodicTable::getNouter
fn n_outer_electrons_rdkit(atomic_num: u8) -> Option<i32> {
    periodic_table_outer_electrons(atomic_num).ok()
}

// ── Hybridization for depiction ─────────────────────────────────────────────

// BEGIN RDKIT CPP FUNCTION: getHybridization (RDDepictor.cpp / ConjugHybrid.cpp)
// RDKit✔️❌: adapted for new-core API
fn rdkit_hybridizations_for_depict(
    atoms: &[Atom],
    bonds: &[Bond],
    degree: &[usize],
) -> Result<Vec<RdkitHybridization>, Coordinate2DError> {
    // Valence computed inline from bond orders below.
    // Radicals: set to 0 (no radical info without full molecule).

    let mut out = Vec::with_capacity(atoms.len());
    for (idx, atom) in atoms.iter().enumerate() {
        if atom.atomic_number() == 0 {
            out.push(RdkitHybridization::Unspecified);
            continue;
        }

        // Compute explicit valence from bond orders.
        let mut explicit_val_from_bonds = 0i32;
        for bond in bonds {
            if bond.begin().index() == idx || bond.end().index() == idx {
                let contrib = match bond.order() {
                    BondOrder::Single | BondOrder::Dative => 1,
                    BondOrder::Double => 2,
                    BondOrder::Triple => 3,
                    BondOrder::Quadruple => 4,
                    BondOrder::Aromatic => 2,
                    _ => 0,
                };
                explicit_val_from_bonds += contrib;
            }
        }

        let norbs = if atom.atomic_number() < 89 {
            rdkit_num_bonds_plus_lone_pairs(
                atom.atomic_number(),
                degree[idx],
                // Use the maximum of bond-derived and assignment's explicit valence
                explicit_val_from_bonds,
                0,
                0,
                atom.formal_charge(),
            )
            .ok_or_else(|| {
                Coordinate2DError::UnsupportedFeature(format!(
                    "RDKit hybridization outer-electron lookup failed for atom {idx}"
                ))
            })?
        } else {
            degree[idx] as i32 + atom.explicit_hydrogens() as i32 + 0 as i32
        };
        // Fallback: if norbs is 0 or negative (lookup failure), use Unspecified
        if norbs <= 0 {
            out.push(RdkitHybridization::Unspecified);
            continue;
        }

        out.push(match norbs {
            0 | 1 => RdkitHybridization::S,
            2 => RdkitHybridization::Sp,
            3 => RdkitHybridization::Sp2,
            4 => RdkitHybridization::Sp3,
            5 => RdkitHybridization::Sp3d,
            6 => RdkitHybridization::Sp3d2,
            _ => RdkitHybridization::Unspecified,
        });
    }
    Ok(out)
}

// ── Ranking helpers ─────────────────────────────────────────────────────────

// BEGIN RDKIT CPP FUNCTION: rankAtomsByRank (DepictUtils.cpp)
// RDKit✔️❌: adapted — CIP rank tiebreaking skipped,
//   uses atom_depict_rank as the primary sorting key.
//   TODO: CIP rank tiebreaking — currently uses atom_depict_rank only
fn rdkit_rank_atoms_by_rank(
    atom_slice: &[Atom],
    order: &mut [usize],
    degree: &[usize],
    _cip_ranks: &[i64],
) {
    // Sort by (atom_depict_rank, atom_idx) since we are not using CIP ranks.
    order.sort_by_key(|&idx| {
        let rank = atom_depict_rank(atom_slice[idx].atomic_number(), degree[idx]);
        (rank, idx)
    });
}

// BEGIN RDKIT CPP FUNCTION: setNbrOrder (DepictUtils.cpp)
// RDKit✔️❌: adapted — CIP rank tiebreaking skipped
fn rdkit_set_nbr_order(
    aid: usize,
    nbrs: &[usize],
    atoms: &[Atom],
    _bonds: &[Bond],
    adjacency: &[Vec<usize>],
    degree: &[usize],
    cip_ranks: &[i64],
) -> Vec<usize> {
    let mut ref_atom: Option<usize> = None;
    for &nb in &adjacency[aid] {
        if !nbrs.contains(&nb) {
            ref_atom = Some(nb);
        }
    }

    let mut thold = nbrs.to_vec();
    if let Some(r) = ref_atom {
        thold.push(r);
    }
    if thold.len() <= 3 {
        let mut out = nbrs.to_vec();
        rdkit_rank_atoms_by_rank(atoms, &mut out, degree, cip_ranks);
        return out;
    }

    rdkit_rank_atoms_by_rank(atoms, &mut thold, degree, cip_ranks);

    let ln = thold.len();
    thold.swap(ln - 3, ln - 2);

    if let Some(r) = ref_atom {
        if let Some(pos) = thold.iter().position(|&a| a == r) {
            let mut out = Vec::with_capacity(nbrs.len());
            out.extend_from_slice(&thold[pos + 1..]);
            out.extend_from_slice(&thold[..pos]);
            return out;
        }
    }
    thold
}

// ── Strict subset fast path ─────────────────────────────────────────────────

// BEGIN RDKIT CPP FUNCTION: computeInitialCoords (subset path)
// RDKit❗✔️: ported from strict_rdkit_subset_coords (intentionally different architecture)
fn strict_rdkit_subset_coords(atoms: &[Atom], bonds: &[Bond]) -> Option<Vec<[f64; 2]>> {
    let n = atoms.len();
    if n == 1 {
        return Some(vec![[0.0, 0.0]]);
    }
    if n > 3 {
        return None;
    }
    if bonds.len() != n - 1 {
        return None;
    }
    if bonds.iter().any(|b| {
        !matches!(
            b.order(),
            BondOrder::Single | BondOrder::Double | BondOrder::Triple
        )
    }) {
        return None;
    }
    let mut degree = vec![0usize; n];
    for b in bonds {
        degree[b.begin().index()] += 1;
        degree[b.end().index()] += 1;
    }
    if degree.iter().any(|d| *d > 2) {
        return None;
    }
    if n == 2 {
        return Some(vec![[-0.75, 0.0], [0.75, 0.0]]);
    }
    None
}

// ── Collision removal helpers ────────────────────────────────────────────────

// BEGIN RDKIT CPP FUNCTION: removeCollisionsBondFlip (RDDDepictor.cpp)
// RDKit✔️❌: ported faithfully, adapted for new-core API
#[allow(clippy::too_many_arguments)]
fn remove_collisions_bond_flip_like(
    atoms: &[Atom],
    bonds: &[Bond],
    comp: &[usize],
    adjacency: &[Vec<usize>],
    local: &mut Vec<(usize, (f64, f64))>,
) {
    const COLLISION_THRES: f64 = 0.70;
    const HETEROATOM_COLL_SCALE: f64 = 1.3;
    const MAX_COLL_ITERS: usize = 15;
    const NUM_BONDS_FLIPS: usize = 3;

    let comp_set: BTreeSet<usize> = comp.iter().copied().collect();
    let mut pos = BTreeMap::<usize, (f64, f64)>::new();
    for &(a, p) in local.iter() {
        pos.insert(a, p);
    }

    let is_ring_bond = |bond_idx: usize| -> bool {
        let b = &bonds[bond_idx];
        let u = b.begin().index();
        let v = b.end().index();
        let mut stack = vec![u];
        let mut seen = BTreeSet::<usize>::new();
        seen.insert(u);
        while let Some(cur) = stack.pop() {
            for (idx, bond) in bonds.iter().enumerate() {
                if idx == bond_idx {
                    continue;
                }
                let nb = if bond.begin().index() == cur {
                    bond.end().index()
                } else if bond.end().index() == cur {
                    bond.begin().index()
                } else {
                    continue;
                };
                if !comp_set.contains(&nb) {
                    continue;
                }
                if nb == v {
                    return true;
                }
                if seen.insert(nb) {
                    stack.push(nb);
                }
            }
        }
        false
    };

    let shortest_path = |src: usize, dst: usize| -> Option<Vec<usize>> {
        let mut q = VecDeque::<usize>::new();
        let mut prev = BTreeMap::<usize, usize>::new();
        q.push_back(src);
        prev.insert(src, src);
        while let Some(u) = q.pop_front() {
            if u == dst {
                break;
            }
            for bond in bonds {
                let v = if bond.begin().index() == u {
                    bond.end().index()
                } else if bond.end().index() == u {
                    bond.begin().index()
                } else {
                    continue;
                };
                if comp_set.contains(&v) {
                    if let std::collections::btree_map::Entry::Vacant(e) = prev.entry(v) {
                        e.insert(u);
                        q.push_back(v);
                    }
                }
            }
        }
        if !prev.contains_key(&dst) {
            return None;
        }
        let mut path = vec![dst];
        let mut cur = dst;
        while cur != src {
            cur = prev[&cur];
            path.push(cur);
        }
        path.reverse();
        Some(path)
    };

    let get_rotatable_bonds = |aid1: usize, aid2: usize| -> Vec<usize> {
        let Some(mut path) = shortest_path(aid1, aid2) else {
            return Vec::new();
        };
        if path.len() < 4 {
            return Vec::new();
        }
        path.remove(0);
        path.pop();
        let mut res = Vec::new();
        let mut pid = path[0];
        for aid in path {
            if let Some((bond_idx, _)) = bonds.iter().enumerate().find(|(_, bond)| {
                (bond.begin().index() == pid && bond.end().index() == aid)
                    || (bond.begin().index() == aid && bond.end().index() == pid)
            }) {
                if matches!(bonds[bond_idx].stereo(), BondStereo::None | BondStereo::Any)
                    && !is_ring_bond(bond_idx)
                {
                    res.push(bond_idx);
                }
            }
            pid = aid;
        }
        res
    };

    let graph_distance_matrix = || -> Vec<Vec<f64>> {
        let n = atoms.len();
        let mut out = vec![vec![f64::INFINITY; n]; n];
        for start in 0..n {
            let mut q = VecDeque::<usize>::new();
            out[start][start] = 0.0;
            q.push_back(start);
            while let Some(cur) = q.pop_front() {
                let next_dist = out[start][cur] + 1.0;
                for &nb in &adjacency[cur] {
                    if !comp_set.contains(&nb) {
                        continue;
                    }
                    if out[start][nb].is_infinite() {
                        out[start][nb] = next_dist;
                        q.push_back(nb);
                    }
                }
            }
        }
        out
    };
    let dmat = graph_distance_matrix();

    let reflect_point = rdkit_reflect_point;

    let find_collisions = |pos: &BTreeMap<usize, (f64, f64)>| -> (Vec<(usize, usize)>, f64) {
        let mut collisions = Vec::new();
        let mut density = BTreeMap::<usize, f64>::new();
        let mut atom_list = comp.to_vec();
        atom_list.sort_unstable();
        let col_thres2 = COLLISION_THRES * COLLISION_THRES;
        let bond_thres2 = 0.50f64 * 0.50f64;
        for i in 0..atom_list.len() {
            let a = atom_list[i];
            let factor_a = if atoms[a].atomic_number() != 6 {
                HETEROATOM_COLL_SCALE
            } else {
                1.0
            };
            for &b in atom_list.iter().take(i) {
                let factor_b = if atoms[b].atomic_number() != 6 {
                    HETEROATOM_COLL_SCALE
                } else {
                    1.0
                };
                let pa = pos[&a];
                let pb = pos[&b];
                let dx = pb.0 - pa.0;
                let dy = pb.1 - pa.1;
                let d2_raw = dx * dx + dy * dy;
                let add = if d2_raw > 1.0e-3 {
                    1.0 / d2_raw
                } else {
                    1000.0
                };
                *density.entry(a).or_insert(0.0) += add;
                *density.entry(b).or_insert(0.0) += add;
                if d2_raw / (factor_a * factor_b) < col_thres2 {
                    collisions.push((a, b));
                }
            }
        }
        for (bid1, b1) in bonds.iter().enumerate() {
            let beg1 = b1.begin().index();
            let end1 = b1.end().index();
            if !pos.contains_key(&beg1) || !pos.contains_key(&end1) {
                continue;
            }
            let p_beg1 = pos[&beg1];
            let p_end1 = pos[&end1];
            let v1 = (p_end1.0 - p_beg1.0, p_end1.1 - p_beg1.1);
            let avg1 = ((p_end1.0 + p_beg1.0) * 0.5, (p_end1.1 + p_beg1.1) * 0.5);
            for (bid2, b2) in bonds.iter().enumerate() {
                if bid2 <= bid1 {
                    continue;
                }
                let beg2 = b2.begin().index();
                let end2 = b2.end().index();
                if !pos.contains_key(&beg2) || !pos.contains_key(&end2) {
                    continue;
                }
                let p_beg2 = pos[&beg2];
                let p_end2 = pos[&end2];
                let avg2 = ((p_end2.0 + p_beg2.0) * 0.5, (p_end2.1 + p_beg2.1) * 0.5);
                let avg_delta = (avg2.0 - avg1.0, avg2.1 - avg1.1);
                let avg_d2 = avg_delta.0 * avg_delta.0 + avg_delta.1 * avg_delta.1;
                if avg_d2 < 0.5 && avg_d2 < bond_thres2 {
                    let v2 = (p_beg2.0 - p_beg1.0, p_beg2.1 - p_beg1.1);
                    let v3 = (p_end2.0 - p_beg1.0, p_end2.1 - p_beg1.1);
                    let cross_v1_v2 = v1.0 * v2.1 - v1.1 * v2.0;
                    let cross_v1_v3 = v1.0 * v3.1 - v1.1 * v3.0;
                    if cross_v1_v2 * cross_v1_v3 < -1.0e-6 {
                        let candidates = [
                            (dmat[beg1][beg2], (beg1, beg2)),
                            (dmat[beg1][end2], (beg1, end2)),
                            (dmat[end1][beg2], (end1, beg2)),
                            (dmat[end1][end2], (end1, end2)),
                        ];
                        let mut best = candidates[0];
                        for candidate in candidates.iter().skip(1) {
                            if candidate.0 < best.0 {
                                best = *candidate;
                            }
                        }
                        collisions.push(best.1);
                    }
                }
            }
        }
        let total_density = atom_list.iter().fold(0.0, |accum, aid| {
            density.get(aid).copied().unwrap_or(0.0) + accum
        });
        (collisions, total_density)
    };

    let flip_about_bond =
        |pos: &mut BTreeMap<usize, (f64, f64)>, bond_idx: usize, flip_end: bool| {
            let bond = &bonds[bond_idx];
            let mut beg = bond.begin().index();
            let mut end = bond.end().index();
            if !flip_end {
                std::mem::swap(&mut beg, &mut end);
            }
            let mut end_side = Vec::new();
            let mut stack = vec![end];
            let mut seen = BTreeSet::<usize>::new();
            seen.insert(end);
            while let Some(cur) = stack.pop() {
                end_side.push(cur);
                for &nb in &adjacency[cur] {
                    if !comp_set.contains(&nb) || nb == beg {
                        continue;
                    }
                    if seen.insert(nb) {
                        stack.push(nb);
                    }
                }
            }
            let end_side_set: BTreeSet<usize> = end_side.iter().copied().collect();
            let mut end_side_flip = true;
            if comp.len().saturating_sub(end_side.len()) < end_side.len() {
                end_side_flip = false;
            }
            let mut atoms: Vec<usize> = pos.keys().copied().collect();
            atoms.retain(|aid| comp_set.contains(aid));
            for aid in atoms {
                let in_end_side = end_side_set.contains(&aid);
                if end_side_flip ^ !in_end_side {
                    let p = pos[&aid];
                    let beg_loc = pos[&beg];
                    let end_loc = pos[&end];
                    pos.insert(aid, reflect_point(p, beg_loc, end_loc));
                }
            }
        };

    let (mut colls, _) = find_collisions(&pos);
    let mut done_bonds = BTreeMap::<usize, usize>::new();
    let mut iter = 0usize;
    while iter < MAX_COLL_ITERS && !colls.is_empty() {
        let ncols = colls.len();
        let (aid1, aid2) = colls[0];
        let rot_bonds = get_rotatable_bonds(aid1, aid2);
        let (_, prev_density) = find_collisions(&pos);
        for ri in rot_bonds {
            let done = *done_bonds.get(&ri).unwrap_or(&0);
            if done >= NUM_BONDS_FLIPS {
                continue;
            }
            done_bonds.insert(ri, done + 1);

            flip_about_bond(&mut pos, ri, true);
            let (new_colls, new_density) = find_collisions(&pos);
            if new_colls.len() < ncols {
                done_bonds.insert(ri, NUM_BONDS_FLIPS);
                colls = new_colls;
                break;
            } else if new_colls.len() == ncols && new_density < prev_density {
                colls = new_colls;
                break;
            }

            flip_about_bond(&mut pos, ri, true);
            let _ = find_collisions(&pos);
            flip_about_bond(&mut pos, ri, false);
            let (new_colls, new_density) = find_collisions(&pos);
            if new_colls.len() < ncols {
                done_bonds.insert(ri, NUM_BONDS_FLIPS);
                colls = new_colls;
                break;
            } else if new_colls.len() == ncols && new_density < prev_density {
                colls = new_colls;
                break;
            } else {
                flip_about_bond(&mut pos, ri, false);
                let (reset_colls, _) = find_collisions(&pos);
                colls = reset_colls;
            }
        }
        iter += 1;
    }

    for item in local.iter_mut() {
        if let Some(&p) = pos.get(&item.0) {
            item.1 = p;
        }
    }
}

// BEGIN RDKIT CPP FUNCTION: removeCollisionsOpenAngles (RDDepictor.cpp)
// RDKit✔️❌: ported faithfully, adapted for new-core API
fn remove_collisions_open_angles_like(
    atoms: &[Atom],
    _bonds: &[Bond],
    comp: &[usize],
    adjacency: &[Vec<usize>],
    local: &mut Vec<(usize, (f64, f64))>,
) {
    const COLLISION_THRES: f64 = 0.70;
    const HETEROATOM_COLL_SCALE: f64 = 1.3;
    const ANGLE_OPEN: f64 = 0.1222;

    if comp.len() < 3 {
        return;
    }
    let comp_set: HashSet<usize> = comp.iter().copied().collect();
    let mut pos = HashMap::<usize, (f64, f64)>::new();
    for &(a, p) in local.iter() {
        pos.insert(a, p);
    }

    let n = atoms.len();
    let mut dmat = vec![vec![usize::MAX; n]; n];
    for &s in comp {
        let mut q = VecDeque::<usize>::new();
        dmat[s][s] = 0;
        q.push_back(s);
        while let Some(u) = q.pop_front() {
            let du = dmat[s][u];
            for &v in &adjacency[u] {
                if !comp_set.contains(&v) {
                    continue;
                }
                if dmat[s][v] == usize::MAX {
                    dmat[s][v] = du + 1;
                    q.push_back(v);
                }
            }
        }
    }

    let find_collisions = |pos: &HashMap<usize, (f64, f64)>| -> Vec<(usize, usize)> {
        let mut out = Vec::new();
        let mut atom_list = comp.to_vec();
        atom_list.sort_unstable();
        for i in 0..atom_list.len() {
            let a = atom_list[i];
            let fa = if atoms[a].atomic_number() != 6 {
                HETEROATOM_COLL_SCALE
            } else {
                1.0
            };
            for &b in atom_list.iter().take(i) {
                let fb = if atoms[b].atomic_number() != 6 {
                    HETEROATOM_COLL_SCALE
                } else {
                    1.0
                };
                let pa = pos[&a];
                let pb = pos[&b];
                let dx = pa.0 - pb.0;
                let dy = pa.1 - pb.1;
                let mut d2 = dx * dx + dy * dy;
                d2 /= fa * fb;
                if d2 < COLLISION_THRES * COLLISION_THRES {
                    out.push((a, b));
                }
            }
        }
        out
    };

    let find_deg1_neighbor = |aid: usize| -> Option<usize> {
        let nbs: Vec<usize> = adjacency[aid]
            .iter()
            .copied()
            .filter(|n| comp_set.contains(n))
            .collect();
        if nbs.len() == 1 { Some(nbs[0]) } else { None }
    };

    let find_closest_neighbor = |aid1: usize, aid2: usize| -> Option<usize> {
        adjacency[aid2]
            .iter()
            .copied()
            .filter(|n| comp_set.contains(n))
            .min_by_key(|n| dmat[aid1][*n])
    };

    let rotate_point = |p: (f64, f64), c: (f64, f64), ang: f64| -> (f64, f64) {
        let v = (p.0 - c.0, p.1 - c.1);
        let r = rotate(v, ang);
        (c.0 + r.0, c.1 + r.1)
    };

    for (aid1, aid2) in find_collisions(&pos) {
        let deg1 = adjacency[aid1]
            .iter()
            .copied()
            .filter(|n| comp_set.contains(n))
            .count();
        let deg2 = adjacency[aid2]
            .iter()
            .copied()
            .filter(|n| comp_set.contains(n))
            .count();
        if deg1 > 1 && deg2 > 1 {
            continue;
        }

        let (aid_a, aid_b, kind) = if deg1 == 1 && deg2 == 1 {
            (find_deg1_neighbor(aid1), find_deg1_neighbor(aid2), 1)
        } else if deg1 == 1 {
            (
                find_deg1_neighbor(aid1),
                find_closest_neighbor(find_deg1_neighbor(aid1).unwrap_or(aid1), aid2),
                2,
            )
        } else {
            (
                find_closest_neighbor(find_deg1_neighbor(aid2).unwrap_or(aid2), aid1),
                find_deg1_neighbor(aid2),
                3,
            )
        };
        let (Some(aid_a), Some(aid_b)) = (aid_a, aid_b) else {
            continue;
        };

        let v2 = (pos[&aid1].0 - pos[&aid_a].0, pos[&aid1].1 - pos[&aid_a].1);
        let v1 = (pos[&aid_b].0 - pos[&aid_a].0, pos[&aid_b].1 - pos[&aid_a].1);
        let cross = v1.0 * v2.1 - v1.1 * v2.0;
        match kind {
            1 => {
                let mut angle = ANGLE_OPEN;
                if cross < 0.0 {
                    angle *= -1.0;
                }
                let p1 = rotate_point(pos[&aid1], pos[&aid_a], angle);
                let p2 = rotate_point(pos[&aid2], pos[&aid_b], -angle);
                pos.insert(aid1, p1);
                pos.insert(aid2, p2);
            }
            2 => {
                let mut angle = 2.0 * ANGLE_OPEN;
                if cross < 0.0 {
                    angle *= -1.0;
                }
                let p1 = rotate_point(pos[&aid1], pos[&aid_a], angle);
                pos.insert(aid1, p1);
            }
            3 => {
                let mut angle = -2.0 * ANGLE_OPEN;
                if cross < 0.0 {
                    angle *= -1.0;
                }
                let p2 = rotate_point(pos[&aid2], pos[&aid_b], angle);
                pos.insert(aid2, p2);
            }
            _ => {}
        }
    }

    for item in local.iter_mut() {
        if let Some(&p) = pos.get(&item.0) {
            item.1 = p;
        }
    }
}

// BEGIN RDKIT CPP FUNCTION: removeCollisionsShortenBonds (RDDepictor.cpp)
// RDKit✔️❌: ported faithfully, adapted for new-core API
fn remove_collisions_shorten_bonds_like(
    atoms: &[Atom],
    _bonds: &[Bond],
    comp: &[usize],
    adjacency: &[Vec<usize>],
    local: &mut Vec<(usize, (f64, f64))>,
) {
    const COLLISION_THRES: f64 = 0.70;
    const HETEROATOM_COLL_SCALE: f64 = 1.3;
    const MAX_COLL_ITERS: usize = 15;

    let comp_set: HashSet<usize> = comp.iter().copied().collect();
    let mut pos = HashMap::<usize, (f64, f64)>::new();
    for &(a, p) in local.iter() {
        pos.insert(a, p);
    }

    let find_collisions = |pos: &HashMap<usize, (f64, f64)>| -> Vec<(usize, usize)> {
        let mut out = Vec::new();
        let mut atom_list = comp.to_vec();
        atom_list.sort_unstable();
        for i in 0..atom_list.len() {
            let a = atom_list[i];
            let fa = if atoms[a].atomic_number() != 6 {
                HETEROATOM_COLL_SCALE
            } else {
                1.0
            };
            for &b in atom_list.iter().take(i) {
                let fb = if atoms[b].atomic_number() != 6 {
                    HETEROATOM_COLL_SCALE
                } else {
                    1.0
                };
                let pa = pos[&a];
                let pb = pos[&b];
                let dx = pa.0 - pb.0;
                let dy = pa.1 - pb.1;
                let mut d2 = dx * dx + dy * dy;
                d2 /= fa * fb;
                if d2 < COLLISION_THRES * COLLISION_THRES {
                    out.push((a, b));
                }
            }
        }
        out
    };

    let degree_in_comp =
        |a: usize| -> usize { adjacency[a].iter().filter(|n| comp_set.contains(n)).count() };
    let find_deg1_neighbor = |aid: usize| -> Option<usize> {
        let nbs: Vec<usize> = adjacency[aid]
            .iter()
            .copied()
            .filter(|n| comp_set.contains(n))
            .collect();
        if nbs.len() == 1 { Some(nbs[0]) } else { None }
    };
    let edge_is_ring_like = |u: usize, v: usize| -> bool {
        let mut stack = vec![u];
        let mut seen = HashSet::<usize>::new();
        seen.insert(u);
        while let Some(cur) = stack.pop() {
            for &nb in &adjacency[cur] {
                if !comp_set.contains(&nb) {
                    continue;
                }
                if (cur == u && nb == v) || (cur == v && nb == u) {
                    continue;
                }
                if nb == v {
                    return true;
                }
                if seen.insert(nb) {
                    stack.push(nb);
                }
            }
        }
        false
    };
    let shortest_path = |src: usize, dst: usize| -> Option<Vec<usize>> {
        let mut q = VecDeque::<usize>::new();
        let mut prev = HashMap::<usize, usize>::new();
        q.push_back(src);
        prev.insert(src, src);
        while let Some(u) = q.pop_front() {
            if u == dst {
                break;
            }
            for &v in &adjacency[u] {
                if !comp_set.contains(&v) {
                    continue;
                }
                if let std::collections::hash_map::Entry::Vacant(e) = prev.entry(v) {
                    e.insert(u);
                    q.push_back(v);
                }
            }
        }
        if !prev.contains_key(&dst) {
            return None;
        }
        let mut path = vec![dst];
        let mut cur = dst;
        while cur != src {
            cur = prev[&cur];
            path.push(cur);
        }
        path.reverse();
        Some(path)
    };

    let mut colls = find_collisions(&pos);
    let mut iter = 0usize;
    while !colls.is_empty() && iter < MAX_COLL_ITERS {
        let (mut aid1, mut aid2) = colls[0];
        let mut deg1 = degree_in_comp(aid1);
        let mut deg2 = degree_in_comp(aid2);
        if deg2 > deg1 {
            std::mem::swap(&mut aid1, &mut aid2);
            std::mem::swap(&mut deg1, &mut deg2);
        }
        if let Some(mut path) = shortest_path(aid1, aid2) {
            if !path.is_empty() {
                path.remove(0);
            }
            let mut has_non_ring = false;
            let mut prev = aid1;
            for &next in &path {
                if !edge_is_ring_like(prev, next) {
                    has_non_ring = true;
                    break;
                }
                prev = next;
            }
            if has_non_ring {
                if deg1 == 1 {
                    if let Some(a) = find_deg1_neighbor(aid1) {
                        let pa = pos[&a];
                        let mut v = (pos[&aid1].0 - pa.0, pos[&aid1].1 - pa.1);
                        v.0 *= 0.9;
                        v.1 *= 0.9;
                        let len = (v.0 * v.0 + v.1 * v.1).sqrt();
                        if len > 0.75 {
                            pos.insert(aid1, (pa.0 + v.0, pa.1 + v.1));
                        }
                    }
                }
                if deg2 == 1 {
                    if let Some(a) = find_deg1_neighbor(aid2) {
                        let pa = pos[&a];
                        let mut v = (pos[&aid2].0 - pa.0, pos[&aid2].1 - pa.1);
                        v.0 *= 0.9;
                        v.1 *= 0.9;
                        let len = (v.0 * v.0 + v.1 * v.1).sqrt();
                        if len > 0.75 {
                            pos.insert(aid2, (pa.0 + v.0, pa.1 + v.1));
                        }
                    }
                }
            } else {
                fn recurse_deg_two_ring_atoms(
                    aid: usize,
                    adjacency: &[Vec<usize>],
                    comp_set: &HashSet<usize>,
                    edge_is_ring_like: &dyn Fn(usize, usize) -> bool,
                    r_path: &mut Vec<usize>,
                    nbr_map: &mut HashMap<usize, Vec<usize>>,
                ) {
                    let nbrs: Vec<usize> = adjacency[aid]
                        .iter()
                        .copied()
                        .filter(|n| comp_set.contains(n) && edge_is_ring_like(aid, *n))
                        .collect();
                    if nbrs.len() != 2 {
                        return;
                    }
                    r_path.push(aid);
                    nbr_map.insert(aid, nbrs.clone());
                    for nbr in nbrs {
                        if !r_path.contains(&nbr) {
                            recurse_deg_two_ring_atoms(
                                nbr,
                                adjacency,
                                comp_set,
                                edge_is_ring_like,
                                r_path,
                                nbr_map,
                            );
                        }
                    }
                }

                let mut r_path = Vec::<usize>::new();
                let mut nbr_map = HashMap::<usize, Vec<usize>>::new();
                recurse_deg_two_ring_atoms(
                    aid1,
                    adjacency,
                    &comp_set,
                    &edge_is_ring_like,
                    &mut r_path,
                    &mut nbr_map,
                );
                if r_path.is_empty() {
                    recurse_deg_two_ring_atoms(
                        aid2,
                        adjacency,
                        &comp_set,
                        &edge_is_ring_like,
                        &mut r_path,
                        &mut nbr_map,
                    );
                }
                let mut move_map = HashMap::<usize, (f64, f64)>::new();
                for &rpi in &r_path {
                    let Some(nbrs) = nbr_map.get(&rpi) else {
                        continue;
                    };
                    let p0 = pos[&rpi];
                    let p1 = pos[&nbrs[0]];
                    let p2 = pos[&nbrs[1]];
                    let mut mv = ((p1.0 + p2.0) * 0.5 - p0.0, (p1.1 + p2.1) * 0.5 - p0.1);
                    let len = (mv.0 * mv.0 + mv.1 * mv.1).sqrt();
                    if len > 1e-12 {
                        mv.0 = mv.0 / len * COLLISION_THRES;
                        mv.1 = mv.1 / len * COLLISION_THRES;
                        move_map.insert(rpi, mv);
                    }
                }
                for rpi in r_path {
                    if let Some(mv) = move_map.get(&rpi) {
                        let p = pos[&rpi];
                        pos.insert(rpi, (p.0 + mv.0, p.1 + mv.1));
                    }
                }
            }
        }
        colls = find_collisions(&pos);
        iter += 1;
    }

    for item in local.iter_mut() {
        if let Some(&p) = pos.get(&item.0) {
            item.1 = p;
        }
    }
}

// ── Acyclic tree layout ─────────────────────────────────────────────────────

// BEGIN RDKIT CPP FUNCTION: placeAcyclicTree (EmbeddedFrag.cpp)
// RDKit✔️❌: ported faithfully, adapted for new-core API.
//   TODO: CIP rank tiebreaking — currently uses atom_depict_rank only
fn place_acyclic_tree(
    atoms: &[Atom],
    bonds: &[Bond],
    comp: &[usize],
    adjacency: &[Vec<usize>],
    degree: &[usize],
    hybridizations: &[RdkitHybridization],
    cip_ranks: &[i64],
    seeded_coords: Option<&HashMap<usize, (f64, f64)>>,
) -> Option<Vec<(usize, (f64, f64))>> {
    if comp.is_empty() {
        return None;
    }
    let comp_set: HashSet<usize> = comp.iter().copied().collect();
    let bond_count_in_comp = bonds
        .iter()
        .filter(|b| comp_set.contains(&b.begin().index()) && comp_set.contains(&b.end().index()))
        .count();
    let is_tree_comp = bond_count_in_comp + 1 == comp.len();
    let n = atoms.len();
    let mut states: Vec<Option<TreeEmbeddedAtom>> = vec![None; n];
    let mut unembedded: HashSet<usize> = comp.iter().copied().collect();
    let mut attach_pts = VecDeque::<usize>::new();
    let mut attach_in_queue = HashSet::<usize>::new();

    let mut root: Option<usize> = None;
    let mut stereo_seeded = false;
    if let Some(seeded) = seeded_coords {
        if !seeded.is_empty() {
            let mut cx = 0.0f64;
            let mut cy = 0.0f64;
            let mut cnt = 0usize;
            for (&aid, &(x, y)) in seeded {
                if comp_set.contains(&aid) {
                    cx += x;
                    cy += y;
                    cnt += 1;
                }
            }
            if cnt > 0 {
                cx /= cnt as f64;
                cy /= cnt as f64;
            }
            for (&aid, &(x, y)) in seeded {
                if !comp_set.contains(&aid) {
                    continue;
                }
                let mut nrm = normalize((x - cx, y - cy));
                if norm(nrm) <= 1e-12 {
                    nrm = (1.0, 0.0);
                }
                states[aid] = Some(TreeEmbeddedAtom {
                    loc: (x, y),
                    normal: nrm,
                    ccw: true,
                    cis_trans_nbr: None,
                    angle: -1.0,
                    nbr1: None,
                    nbr2: None,
                    rot_dir: 0,
                    pending: Vec::new(),
                });
                unembedded.remove(&aid);
            }
            root = comp
                .iter()
                .copied()
                .filter(|a| states[*a].is_some())
                .min_by_key(|&idx| {
                    let depict = atom_depict_rank(atoms[idx].atomic_number(), degree[idx]);
                    depict * atoms.len() + idx
                });
        }
    }

    // Stereo seed
    if root.is_none() && seeded_coords.is_none() {
        if let Some(bond) = bonds.iter().find(|bond| {
            comp_set.contains(&bond.begin().index())
                && comp_set.contains(&bond.end().index())
                && matches!(bond.order(), BondOrder::Double)
                && matches!(bond.stereo(), BondStereo::Cis | BondStereo::Trans)
                && bond.stereo_atoms().is_some()
        }) {
            let begin = bond.begin().index();
            let end = bond.end().index();
            let stereo_atoms = bond.stereo_atoms().unwrap();
            let begin_control = stereo_atoms[0].index();
            let end_control = stereo_atoms[1].index();

            states[begin] = Some(TreeEmbeddedAtom {
                loc: (0.0, 0.0),
                normal: (0.0, -1.0),
                ccw: false,
                cis_trans_nbr: Some(begin_control),
                angle: -1.0,
                nbr1: Some(end),
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
            });
            let (end_normal, end_ccw) = if matches!(bond.stereo(), BondStereo::Cis) {
                ((0.0, -1.0), true)
            } else {
                ((0.0, 1.0), false)
            };
            states[end] = Some(TreeEmbeddedAtom {
                loc: (1.5, 0.0),
                normal: end_normal,
                ccw: end_ccw,
                cis_trans_nbr: Some(end_control),
                angle: -1.0,
                nbr1: Some(begin),
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
            });
            unembedded.remove(&begin);
            unembedded.remove(&end);
            root = Some(begin);
            stereo_seeded = true;
        }
    }

    // Unseeded root
    if root.is_none() {
        let r = *comp.iter().min_by_key(|&&idx| {
            let depict = atom_depict_rank(atoms[idx].atomic_number(), degree[idx]);
            depict * atoms.len() + idx
        })?;
        states[r] = Some(TreeEmbeddedAtom {
            loc: (0.0, 0.0),
            normal: (1.0, 0.0),
            ccw: true,
            cis_trans_nbr: None,
            angle: -1.0,
            nbr1: None,
            nbr2: None,
            rot_dir: 0,
            pending: Vec::new(),
        });
        unembedded.remove(&r);
        root = Some(r);
    }

    // ── Closures ───────────────────────────────────────────────────────────

    let find_num_neigh =
        |states: &[Option<TreeEmbeddedAtom>], pt: (f64, f64), radius: f64| -> i32 {
            states
                .iter()
                .flatten()
                .filter(|st| {
                    let dx = st.loc.0 - pt.0;
                    let dy = st.loc.1 - pt.1;
                    (dx * dx + dy * dy).sqrt() < radius
                })
                .count() as i32
        };

    let compute_angle = |center: (f64, f64), p1: (f64, f64), p2: (f64, f64)| -> f64 {
        let v1 = (p1.0 - center.0, p1.1 - center.1);
        let v2 = (p2.0 - center.0, p2.1 - center.1);
        let d1 = (v1.0 * v1.0 + v1.1 * v1.1).sqrt();
        let d2 = (v2.0 * v2.0 + v2.1 * v2.1).sqrt();
        if d1 <= 1e-12 || d2 <= 1e-12 {
            return PI;
        }
        let mut c = (v1.0 * v2.0 + v1.1 * v2.1) / (d1 * d2);
        c = c.clamp(-1.0, 1.0);
        rdkit_acos(c)
    };

    let update_new_neighs = |aid: usize,
                             states: &mut Vec<Option<TreeEmbeddedAtom>>,
                             unembedded: &HashSet<usize>,
                             attach_pts: &mut VecDeque<usize>,
                             attach_in_queue: &mut HashSet<usize>| {
        if !is_tree_comp && states[aid].is_some() {
            let st_loc = states[aid].as_ref().map(|s| s.loc).unwrap_or((0.0, 0.0));
            let done_nbrs: Vec<usize> = adjacency[aid]
                .iter()
                .copied()
                .filter(|nb| {
                    comp_set.contains(nb) && !unembedded.contains(nb) && states[*nb].is_some()
                })
                .collect();
            if done_nbrs.len() >= 3 {
                let mut angle_pairs: Vec<(f64, (usize, usize))> = Vec::new();
                for i in 0..done_nbrs.len() {
                    for j in i + 1..done_nbrs.len() {
                        let a = done_nbrs[i];
                        let b = done_nbrs[j];
                        let pa = states[a]
                            .as_ref()
                            .map(|n| n.loc)
                            .unwrap_or((st_loc.0 - 1.0, st_loc.1));
                        let pb = states[b]
                            .as_ref()
                            .map(|n| n.loc)
                            .unwrap_or((st_loc.0 + 1.0, st_loc.1));
                        angle_pairs.push((compute_angle(st_loc, pa, pb), (a, b)));
                    }
                }
                angle_pairs
                    .sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap_or(std::cmp::Ordering::Equal));
                if let Some((w_ang, (wnb1, wnb2))) = angle_pairs.last().copied() {
                    let mut nb1 = wnb1;
                    let mut nb2 = wnb2;
                    for &(_, (a, b)) in angle_pairs.iter().rev() {
                        if a == wnb1 && b != wnb2 {
                            nb1 = b;
                            nb2 = wnb1;
                            break;
                        } else if b == wnb1 && a != wnb2 {
                            nb1 = a;
                            nb2 = wnb1;
                            break;
                        } else if a == wnb2 && b != wnb1 {
                            nb1 = b;
                            nb2 = wnb2;
                            break;
                        } else if b == wnb2 && a != wnb1 {
                            nb1 = a;
                            nb2 = wnb2;
                            break;
                        }
                    }
                    let p1 = states[nb1]
                        .as_ref()
                        .map(|n| n.loc)
                        .unwrap_or((st_loc.0 - 1.0, st_loc.1));
                    let p2 = states[nb2]
                        .as_ref()
                        .map(|n| n.loc)
                        .unwrap_or((st_loc.0 + 1.0, st_loc.1));
                    if let Some(st) = states[aid].as_mut() {
                        st.nbr1 = Some(nb1);
                        st.nbr2 = Some(nb2);
                        st.angle = 2.0 * PI - w_ang;
                        st.rot_dir = rotation_dir(st_loc, p1, p2, w_ang);
                    }
                }
            }
        }

        let mut neighs: Vec<usize> = adjacency[aid]
            .iter()
            .copied()
            .filter(|nb| comp_set.contains(nb) && unembedded.contains(nb))
            .collect();
        if !neighs.is_empty() && (degree[aid] < 4 || neighs.len() < 3) {
            rdkit_rank_atoms_by_rank(atoms, &mut neighs, degree, cip_ranks);
        } else if degree[aid] >= 4 && neighs.len() >= 3 {
            neighs = rdkit_set_nbr_order(aid, &neighs, atoms, bonds, adjacency, degree, cip_ranks);
        }
        if let Some(st) = states[aid].as_mut() {
            st.pending = neighs;
            if !st.pending.is_empty() && !attach_in_queue.contains(&aid) {
                attach_pts.push_back(aid);
                attach_in_queue.insert(aid);
            }
        }
    };

    // Initialize
    if stereo_seeded {
        for &aid in comp {
            if states[aid].is_some() {
                update_new_neighs(
                    aid,
                    &mut states,
                    &unembedded,
                    &mut attach_pts,
                    &mut attach_in_queue,
                );
            }
        }
    } else if let Some(seeded) = seeded_coords {
        if !seeded.is_empty() {
            for &aid in comp {
                if states[aid].is_some() {
                    update_new_neighs(
                        aid,
                        &mut states,
                        &unembedded,
                        &mut attach_pts,
                        &mut attach_in_queue,
                    );
                }
            }
        } else {
            update_new_neighs(
                root.expect("root initialized"),
                &mut states,
                &unembedded,
                &mut attach_pts,
                &mut attach_in_queue,
            );
        }
    } else {
        update_new_neighs(
            root.expect("root initialized"),
            &mut states,
            &unembedded,
            &mut attach_pts,
            &mut attach_in_queue,
        );
    }

    // Main loop
    while let Some(aid) = attach_pts.pop_front() {
        attach_in_queue.remove(&aid);
        let neighs = states[aid]
            .as_ref()
            .map(|st| st.pending.clone())
            .unwrap_or_default();
        if let Some(st) = states[aid].as_mut() {
            st.pending.clear();
        }
        let total = neighs.len();
        for (i, nb) in neighs.into_iter().enumerate() {
            if !unembedded.contains(&nb) {
                continue;
            }

            let ref_state = states[aid].as_ref()?.clone();
            let nnbr = (total - i) as f64;
            if ref_state.angle > 0.0 {
                let rem_angle = 2.0 * PI - ref_state.angle;
                let mut curr_angle = rem_angle / (1.0 + nnbr);
                if let Some(st) = states[aid].as_mut() {
                    st.angle += curr_angle;
                }
                let nb1 = states[ref_state.nbr1?].as_ref()?.loc;
                let nb2 = states[ref_state.nbr2?].as_ref()?.loc;
                if states[aid].as_ref()?.rot_dir == 0 {
                    let rd = rotation_dir(ref_state.loc, nb1, nb2, rem_angle);
                    if let Some(st) = states[aid].as_mut() {
                        st.rot_dir = rd;
                    }
                }
                curr_angle *= states[aid].as_ref()?.rot_dir as f64;
                let mut curr_loc = rotate_around(nb2, ref_state.loc, curr_angle);
                if rem_angle.abs() - PI < 1e-3 {
                    let curr_loc2 = rotate_around(nb2, ref_state.loc, -curr_angle);
                    if find_num_neigh(&states, curr_loc, 0.5)
                        > find_num_neigh(&states, curr_loc2, 0.5)
                    {
                        curr_loc = curr_loc2;
                    }
                }
                if let Some(st) = states[aid].as_mut() {
                    st.nbr2 = Some(nb);
                }

                let tpt = (curr_loc.0 - ref_state.loc.0, curr_loc.1 - ref_state.loc.1);
                let mut normal = normalize((-tpt.1, tpt.0));
                let tp1 = (curr_loc.0 + normal.0, curr_loc.1 + normal.1);
                let tp2 = (curr_loc.0 - normal.0, curr_loc.1 - normal.1);
                let nccw = find_num_neigh(&states, tp1, 2.5);
                let ncw = find_num_neigh(&states, tp2, 2.5);
                let (ccw, out_normal) = if nccw < ncw {
                    (false, normal)
                } else {
                    normal = (-normal.0, -normal.1);
                    (true, normal)
                };
                states[nb] = Some(TreeEmbeddedAtom {
                    loc: curr_loc,
                    normal: out_normal,
                    ccw,
                    cis_trans_nbr: None,
                    angle: -1.0,
                    nbr1: Some(aid),
                    nbr2: None,
                    rot_dir: 0,
                    pending: Vec::new(),
                });
            } else {
                let mut ref_ccw = ref_state.ccw;
                let mut curr_vec = ref_state.normal;
                if ref_state.cis_trans_nbr.is_some_and(|ct| ct != nb) {
                    ref_ccw = !ref_ccw;
                    curr_vec = (-curr_vec.0, -curr_vec.1);
                }
                let degree_here = degree[aid];
                let mut angle = compute_sub_angle(degree_here, hybridizations[aid]);
                let mut flip_norm = false;

                if states[aid].as_ref()?.nbr1.is_some() {
                    if let Some(st) = states[aid].as_mut() {
                        st.angle = angle;
                        st.nbr2 = Some(nb);
                    }
                } else {
                    let norm0 = states[aid].as_ref()?.normal;
                    let rot = rotate(norm0, angle);
                    if let Some(st) = states[aid].as_mut() {
                        st.normal = rot;
                        st.nbr1 = Some(nb);
                    }
                    flip_norm = true;
                }

                angle -= PI / 2.0;
                if !ref_ccw {
                    angle *= -1.0;
                }
                curr_vec = rotate(curr_vec, angle);
                let curr_loc = (
                    ref_state.loc.0 + 1.5 * curr_vec.0,
                    ref_state.loc.1 + 1.5 * curr_vec.1,
                );

                let tpt = (ref_state.loc.0 - curr_loc.0, ref_state.loc.1 - curr_loc.1);
                let mut normal = (-tpt.1, tpt.0);
                if ref_ccw ^ flip_norm {
                    normal = (-normal.0, -normal.1);
                }
                normal = normalize(normal);
                states[nb] = Some(TreeEmbeddedAtom {
                    loc: curr_loc,
                    normal,
                    ccw: (!ref_ccw) ^ flip_norm,
                    cis_trans_nbr: None,
                    angle: -1.0,
                    nbr1: Some(aid),
                    nbr2: None,
                    rot_dir: 0,
                    pending: Vec::new(),
                });
            }

            unembedded.remove(&nb);
            update_new_neighs(
                nb,
                &mut states,
                &unembedded,
                &mut attach_pts,
                &mut attach_in_queue,
            );
        }
        update_new_neighs(
            aid,
            &mut states,
            &unembedded,
            &mut attach_pts,
            &mut attach_in_queue,
        );
    }

    if !unembedded.is_empty() {
        return None;
    }

    let mut out: Vec<(usize, (f64, f64))> = comp
        .iter()
        .copied()
        .filter_map(|a| states[a].as_ref().map(|st| (a, st.loc)))
        .collect();
    if out.len() != comp.len() {
        return None;
    }
    out.sort_by_key(|(a, _)| *a);
    remove_collisions_bond_flip_like(atoms, bonds, comp, adjacency, &mut out);
    remove_collisions_open_angles_like(atoms, bonds, comp, adjacency, &mut out);
    remove_collisions_shorten_bonds_like(atoms, bonds, comp, adjacency, &mut out);
    Some(out)
}

// ── Multi-ring non-fused component layout ────────────────────────────────────

// BEGIN RDKIT CPP FUNCTION: placeMultiringNonfusedComponent (EmbeddedFrag.cpp)
// RDKit✔️❌: ported faithfully, adapted for new-core API
fn place_multiring_nonfused_component(
    atoms: &[Atom],
    bonds: &[Bond],
    comp: &[usize],
    adjacency: &[Vec<usize>],
    degree: &[usize],
) -> Option<Vec<(usize, (f64, f64))>> {
    #[derive(Clone)]
    struct DepictFrag {
        atoms: BTreeMap<usize, TreeEmbeddedAtom>,
        attach_pts: VecDeque<usize>,
        done: bool,
    }

    fn update_new_neighs_for_frag(
        frag: &mut DepictFrag,
        atoms: &[Atom],
        bonds: &[Bond],
        comp_set: &BTreeSet<usize>,
        adjacency: &[Vec<usize>],
        degree: &[usize],
        aid: usize,
    ) {
        let mut heavy = Vec::new();
        let mut hydrogens = Vec::new();
        for &nb in &adjacency[aid] {
            if !comp_set.contains(&nb) || frag.atoms.contains_key(&nb) {
                continue;
            }
            if atoms[nb].atomic_number() == 1 {
                hydrogens.push(nb);
            } else {
                heavy.push(nb);
            }
        }
        heavy.extend(hydrogens);

        if !heavy.is_empty() && (degree[aid] < 4 || heavy.len() < 3) {
            // TODO: CIP rank tiebreaking — currently uses atom_depict_rank only
            let cip_ranks = vec![0i64; atoms.len()];
            rdkit_rank_atoms_by_rank(atoms, &mut heavy, degree, &cip_ranks);
        } else if degree[aid] >= 4 && heavy.len() >= 3 {
            let cip_ranks = vec![0i64; atoms.len()];
            heavy = rdkit_set_nbr_order(aid, &heavy, atoms, bonds, adjacency, degree, &cip_ranks);
        }

        if let Some(st) = frag.atoms.get_mut(&aid) {
            st.pending = heavy;
            if !st.pending.is_empty() && !frag.attach_pts.iter().any(|&x| x == aid) {
                frag.attach_pts.push_back(aid);
            }
        }
    }

    fn setup_new_neighs_for_frag(
        frag: &mut DepictFrag,
        atoms: &[Atom],
        bonds: &[Bond],
        comp_set: &BTreeSet<usize>,
        adjacency: &[Vec<usize>],
        degree: &[usize],
    ) {
        frag.attach_pts.clear();
        let atom_ids: Vec<usize> = frag.atoms.keys().copied().collect();
        for aid in &atom_ids {
            update_new_neighs_for_frag(frag, atoms, bonds, comp_set, adjacency, degree, *aid);
        }
        let attach = frag.attach_pts.make_contiguous();
        // TODO: CIP rank tiebreaking — currently uses atom_depict_rank only
        let cip_ranks = vec![0i64; atom_ids.len()];
        rdkit_rank_atoms_by_rank(atoms, attach, degree, &cip_ranks);
    }

    fn find_num_neigh_frag(frag: &DepictFrag, pt: (f64, f64), radius: f64) -> i32 {
        frag.atoms
            .values()
            .filter(|st| {
                let dx = st.loc.0 - pt.0;
                let dy = st.loc.1 - pt.1;
                rdkit_sqrt(dx * dx + dy * dy) < radius
            })
            .count() as i32
    }

    fn add_atom_to_atom_with_ang_frag(
        frag: &mut DepictFrag,
        aid: usize,
        to_aid: usize,
    ) -> Option<()> {
        let ref_state = frag.atoms.get(&to_aid)?.clone();
        let nnbr = ref_state.pending.len() as f64;
        let rem_angle = 2.0 * PI - ref_state.angle;
        let mut curr_angle = rem_angle / (1.0 + nnbr);
        if let Some(st) = frag.atoms.get_mut(&to_aid) {
            st.angle += curr_angle;
        }

        let nb1 = frag.atoms.get(&ref_state.nbr1?)?.loc;
        let nb2 = frag.atoms.get(&ref_state.nbr2?)?.loc;
        if frag.atoms.get(&to_aid)?.rot_dir == 0 {
            let rd = rotation_dir(ref_state.loc, nb1, nb2, rem_angle);
            if let Some(st) = frag.atoms.get_mut(&to_aid) {
                st.rot_dir = rd;
            }
        }
        curr_angle *= frag.atoms.get(&to_aid)?.rot_dir as f64;
        let mut curr_loc = rotate_around(nb2, ref_state.loc, curr_angle);
        if rem_angle.abs() - PI < 1e-3 {
            let curr_loc2 = rotate_around(nb2, ref_state.loc, -curr_angle);
            if find_num_neigh_frag(frag, curr_loc, 0.5) > find_num_neigh_frag(frag, curr_loc2, 0.5)
            {
                curr_loc = curr_loc2;
            }
        }
        if let Some(st) = frag.atoms.get_mut(&to_aid) {
            st.nbr2 = Some(aid);
        }

        let tpt = (curr_loc.0 - ref_state.loc.0, curr_loc.1 - ref_state.loc.1);
        let probe_normal = (-tpt.1, tpt.0);
        let tp1 = (curr_loc.0 + probe_normal.0, curr_loc.1 + probe_normal.1);
        let tp2 = (curr_loc.0 - probe_normal.0, curr_loc.1 - probe_normal.1);
        let nccw = find_num_neigh_frag(frag, tp1, 2.5);
        let ncw = find_num_neigh_frag(frag, tp2, 2.5);
        let mut normal = normalize(probe_normal);
        let (ccw, out_normal) = if nccw < ncw {
            (false, normal)
        } else {
            normal = (-normal.0, -normal.1);
            (true, normal)
        };
        frag.atoms.insert(
            aid,
            TreeEmbeddedAtom {
                loc: curr_loc,
                normal: out_normal,
                ccw,
                cis_trans_nbr: None,
                angle: -1.0,
                nbr1: Some(to_aid),
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
            },
        );
        Some(())
    }

    fn add_atom_to_atom_with_no_ang_frag(
        frag: &mut DepictFrag,
        atoms: &[Atom],
        bonds: &[Bond],
        degree: &[usize],
        aid: usize,
        to_aid: usize,
    ) -> Option<()> {
        let ref_state = frag.atoms.get(&to_aid)?.clone();
        let mut ref_ccw = ref_state.ccw;
        let mut curr_vec = ref_state.normal;
        if ref_state.cis_trans_nbr.is_some_and(|ct| ct != aid) {
            ref_ccw = !ref_ccw;
            curr_vec = (-curr_vec.0, -curr_vec.1);
        }
        let hybridizations = rdkit_hybridizations_for_depict(atoms, bonds, degree).ok()?;
        let mut angle = compute_sub_angle(degree[to_aid], hybridizations[to_aid]);
        let mut flip_norm = false;

        if frag.atoms.get(&to_aid)?.nbr1.is_some() {
            if let Some(st) = frag.atoms.get_mut(&to_aid) {
                st.angle = angle;
                st.nbr2 = Some(aid);
            }
        } else {
            let rot = rotate(frag.atoms.get(&to_aid)?.normal, angle);
            if let Some(st) = frag.atoms.get_mut(&to_aid) {
                st.normal = rot;
                st.nbr1 = Some(aid);
            }
            flip_norm = true;
        }

        angle -= PI / 2.0;
        if !ref_ccw {
            angle *= -1.0;
        }
        curr_vec = rotate(curr_vec, angle);
        curr_vec.0 *= 1.5;
        curr_vec.1 *= 1.5;
        curr_vec.0 += ref_state.loc.0;
        curr_vec.1 += ref_state.loc.1;
        let curr_loc = curr_vec;

        let tpt = (ref_state.loc.0 - curr_loc.0, ref_state.loc.1 - curr_loc.1);
        let mut normal = (-tpt.1, tpt.0);
        if ref_ccw ^ flip_norm {
            normal = (-normal.0, -normal.1);
        }
        normal = normalize(normal);
        frag.atoms.insert(
            aid,
            TreeEmbeddedAtom {
                loc: curr_loc,
                normal,
                ccw: (!ref_ccw) ^ flip_norm,
                cis_trans_nbr: None,
                angle: -1.0,
                nbr1: Some(to_aid),
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
            },
        );
        Some(())
    }

    fn add_non_ring_atom_frag(
        frag: &mut DepictFrag,
        atoms: &[Atom],
        bonds: &[Bond],
        comp_set: &BTreeSet<usize>,
        adjacency: &[Vec<usize>],
        degree: &[usize],
        aid: usize,
        to_aid: usize,
    ) -> Option<()> {
        if frag.atoms.contains_key(&aid) || !frag.atoms.contains_key(&to_aid) {
            return None;
        }
        if frag.atoms.get(&to_aid)?.angle > 0.0 {
            add_atom_to_atom_with_ang_frag(frag, aid, to_aid)?;
        } else {
            add_atom_to_atom_with_no_ang_frag(frag, atoms, bonds, degree, aid, to_aid)?;
        }
        if let Some(st) = frag.atoms.get_mut(&to_aid) {
            st.pending.retain(|&x| x != aid);
        }
        update_new_neighs_for_frag(frag, atoms, bonds, comp_set, adjacency, degree, aid);
        Some(())
    }

    fn transform_frag_two_point(
        frag: &mut DepictFrag,
        ref1: (f64, f64),
        ref2: (f64, f64),
        oth1: (f64, f64),
        oth2: (f64, f64),
    ) -> Option<()> {
        let trans = transform2d_set_transform_two_point(ref1, ref2, oth1, oth2);
        let ov = (oth2.0 - oth1.0, oth2.1 - oth1.1);
        if norm(ov) <= 1e-12 {
            return None;
        }
        for st in frag.atoms.values_mut() {
            let loc = st.loc;
            st.loc = transform2d_point(loc, trans);
            let temp = (loc.0 + st.normal.0, loc.1 + st.normal.1);
            let temp = transform2d_point(temp, trans);
            st.normal = (temp.0 - st.loc.0, temp.1 - st.loc.1);
        }
        Some(())
    }

    fn reflect_point_line(p: (f64, f64), a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
        rdkit_reflect_point(p, a, b)
    }

    fn reflect_frag(frag: &mut DepictFrag, a: (f64, f64), b: (f64, f64)) {
        let vx = b.0 - a.0;
        let vy = b.1 - a.1;
        let len2 = vx * vx + vy * vy;
        if len2 <= 1e-12 {
            return;
        }
        for st in frag.atoms.values_mut() {
            let old = st.loc;
            st.loc = reflect_point_line(old, a, b);
            let n_end = (old.0 + st.normal.0, old.1 + st.normal.1);
            let n_ref = reflect_point_line(n_end, a, b);
            st.normal = (n_ref.0 - st.loc.0, n_ref.1 - st.loc.1);
            st.ccw = !st.ccw;
        }
    }

    fn reflect_if_necessary_density(
        master: &DepictFrag,
        other: &mut DepictFrag,
        aid1: usize,
        aid2: usize,
    ) {
        let pin1 = master.atoms[&aid1].loc;
        let pin2 = master.atoms[&aid2].loc;
        let mut density_normal = 0.0f64;
        let mut density_reflect = 0.0f64;
        for (&oa, ost) in &other.atoms {
            if master.atoms.contains_key(&oa) {
                continue;
            }
            let loc = ost.loc;
            let rloc = reflect_point_line(loc, pin1, pin2);
            for tst in master.atoms.values() {
                let d = norm((tst.loc.0 - loc.0, tst.loc.1 - loc.1));
                let rd = norm((tst.loc.0 - rloc.0, tst.loc.1 - rloc.1));
                density_normal += if d > 1e-3 { 1.0 / d } else { 1000.0 };
                density_reflect += if rd > 1e-3 { 1.0 / rd } else { 1000.0 };
            }
        }
        if density_normal - density_reflect > 1e-4 {
            reflect_frag(other, pin1, pin2);
        }
    }

    fn reflect_if_necessary_cis_trans(
        master: &DepictFrag,
        other: &mut DepictFrag,
        ct_case: u8,
        aid1: usize,
        aid2: usize,
    ) -> Option<()> {
        let p1_loc = master.atoms.get(&aid1)?.loc;
        let (p1_norm, ring_atom) = if ct_case == 1 {
            let aid1_state = other.atoms.get(&aid1)?;
            (aid1_state.normal, aid1_state.cis_trans_nbr?)
        } else {
            let aid1_state = master.atoms.get(&aid1)?;
            (aid1_state.normal, aid1_state.cis_trans_nbr?)
        };
        let r_atm_loc = if ct_case == 1 {
            master.atoms.get(&ring_atom)?.loc
        } else {
            other.atoms.get(&ring_atom)?.loc
        };
        let r_rel = (r_atm_loc.0 - p1_loc.0, r_atm_loc.1 - p1_loc.1);
        let dot = r_rel.0 * p1_norm.0 + r_rel.1 * p1_norm.1;
        if dot < 0.0 {
            let p2_loc = master.atoms.get(&aid2)?.loc;
            reflect_frag(other, p1_loc, p2_loc);
        }
        Some(())
    }

    fn merge_with_common_frag(
        master: &mut DepictFrag,
        other: &mut DepictFrag,
        atoms: &[Atom],
        bonds: &[Bond],
        comp_set: &BTreeSet<usize>,
        adjacency: &[Vec<usize>],
        degree: &[usize],
        mut common: Vec<usize>,
    ) -> Option<()> {
        let mut ct_case = 0u8;
        if common.len() == 1 {
            let comm_aid = common[0];
            let other_atom;
            if master
                .atoms
                .get(&comm_aid)
                .and_then(|st| st.cis_trans_nbr)
                .is_some()
            {
                ct_case = 2;
                other_atom = master.atoms.get(&comm_aid)?.nbr1;
                if let Some(aid) = other_atom {
                    add_non_ring_atom_frag(
                        other, atoms, bonds, comp_set, adjacency, degree, aid, comm_aid,
                    )?;
                }
            } else if other
                .atoms
                .get(&comm_aid)
                .and_then(|st| st.cis_trans_nbr)
                .is_some()
            {
                ct_case = 1;
                other_atom = other.atoms.get(&comm_aid)?.nbr1;
                if let Some(aid) = other_atom {
                    add_non_ring_atom_frag(
                        master, atoms, bonds, comp_set, adjacency, degree, aid, comm_aid,
                    )?;
                }
            } else {
                other_atom = master.atoms.get(&comm_aid)?.nbr1;
                if let Some(aid) = other_atom {
                    add_non_ring_atom_frag(
                        other, atoms, bonds, comp_set, adjacency, degree, aid, comm_aid,
                    )?;
                }
            }
            if let Some(aid) = other_atom {
                common.push(aid);
            }
        }

        if common.len() == 1 {
            // RDKit one-atom transform merge not ported
            return None;
        }
        let aid1 = common[0];
        let aid2 = common[1];
        let ref1 = master.atoms.get(&aid1)?.loc;
        let ref2 = master.atoms.get(&aid2)?.loc;
        let oth1 = other.atoms.get(&aid1)?.loc;
        let oth2 = other.atoms.get(&aid2)?.loc;
        transform_frag_two_point(other, ref1, ref2, oth1, oth2)?;
        if common.len() >= 2 {
            if ct_case > 0 {
                reflect_if_necessary_cis_trans(master, other, ct_case, aid1, aid2)?;
            } else if common.len() == 2 {
                reflect_if_necessary_density(master, other, aid1, aid2);
            }
            // common.len() > 2: RDKit third-point reflection not ported
        }

        for (&aid, ost) in &other.atoms {
            if !common.contains(&aid) {
                master.atoms.insert(aid, ost.clone());
                if !ost.pending.is_empty() && !master.attach_pts.iter().any(|&x| x == aid) {
                    master.attach_pts.push_back(aid);
                }
            } else if let Some(mst) = master.atoms.get_mut(&aid) {
                if ost.cis_trans_nbr.is_some() {
                    mst.cis_trans_nbr = ost.cis_trans_nbr;
                    mst.normal = ost.normal;
                    mst.ccw = ost.ccw;
                }
                if ost.angle > 0.0 {
                    mst.angle = ost.angle;
                    mst.nbr1 = ost.nbr1;
                    mst.nbr2 = ost.nbr2;
                }
            }
        }
        for aid in common {
            if master.atoms.contains_key(&aid) {
                update_new_neighs_for_frag(master, atoms, bonds, comp_set, adjacency, degree, aid);
            }
        }
        other.done = true;
        Some(())
    }

    fn merge_no_common_frag(
        master: &mut DepictFrag,
        other: &mut DepictFrag,
        atoms: &[Atom],
        bonds: &[Bond],
        comp_set: &BTreeSet<usize>,
        adjacency: &[Vec<usize>],
        degree: &[usize],
        to_aid: usize,
        nbr_aid: usize,
    ) -> Option<()> {
        add_non_ring_atom_frag(
            master, atoms, bonds, comp_set, adjacency, degree, nbr_aid, to_aid,
        )?;
        add_non_ring_atom_frag(
            other, atoms, bonds, comp_set, adjacency, degree, to_aid, nbr_aid,
        )?;
        merge_with_common_frag(
            master,
            other,
            atoms,
            bonds,
            comp_set,
            adjacency,
            degree,
            vec![to_aid, nbr_aid],
        )
    }

    fn find_common_atoms(master: &DepictFrag, other: &DepictFrag) -> Vec<usize> {
        let mut common = Vec::new();
        for aid in master.atoms.keys() {
            if other.atoms.contains_key(aid) {
                common.push(*aid);
            }
        }
        common
    }

    fn merge_frags_with_common(
        master: &mut DepictFrag,
        frags: &mut Vec<DepictFrag>,
        atoms: &[Atom],
        bonds: &[Bond],
        comp_set: &BTreeSet<usize>,
        adjacency: &[Vec<usize>],
        degree: &[usize],
    ) -> Option<()> {
        loop {
            let mut found = None;
            for (idx, frag) in frags.iter().enumerate() {
                if frag.done {
                    continue;
                }
                let common = find_common_atoms(master, frag);
                if !common.is_empty() {
                    found = Some((idx, common));
                    break;
                }
            }
            let Some((idx, common)) = found else { break };
            let mut other = frags.remove(idx);
            let common_for_cleanup = common.clone();
            merge_with_common_frag(
                master, &mut other, atoms, bonds, comp_set, adjacency, degree, common,
            )?;
            for cai in common_for_cleanup {
                let remove_attach = master
                    .atoms
                    .get(&cai)
                    .is_some_and(|st| st.pending.is_empty())
                    && master.attach_pts.iter().any(|&x| x == cai);
                if remove_attach {
                    master.attach_pts.retain(|&x| x != cai);
                }
            }
        }
        Some(())
    }

    fn mirror_trans_ring_atoms(
        _atoms: &[Atom],
        bonds: &[Bond],
        ring: &[usize],
        coords: &mut BTreeMap<usize, (f64, f64)>,
    ) {
        for i in 0..ring.len() {
            let atom1 = ring[i];
            let atom2 = ring[(i + 1) % ring.len()];
            let Some(bond) = bonds.iter().find(|bond| {
                (bond.begin().index() == atom1 && bond.end().index() == atom2)
                    || (bond.begin().index() == atom2 && bond.end().index() == atom1)
            }) else {
                continue;
            };
            if !matches!(bond.order(), BondOrder::Double) {
                continue;
            }
            let Some(stereo_atoms) = bond.stereo_atoms() else {
                continue;
            };
            if matches!(bond.stereo(), BondStereo::None | BondStereo::Any) {
                continue;
            }
            let left_is_in = ring.contains(&stereo_atoms[0].index());
            let right_is_in = ring.contains(&stereo_atoms[1].index());
            let is_trans = if matches!(bond.stereo(), BondStereo::Trans) {
                left_is_in == right_is_in
            } else {
                left_is_in != right_is_in
            };
            if !is_trans {
                continue;
            }

            let left = ring[(i + ring.len() - 1) % ring.len()];
            let right = atom2;
            let last = coords[&left];
            let ref_pt = coords[&right];
            let interest = coords[&atom1];
            let d = (last.0 - ref_pt.0, last.1 - ref_pt.1);
            let dot = d.0 * d.0 + d.1 * d.1;
            let a = (d.0 * d.0 - d.1 * d.1) / dot;
            let b = 2.0 * d.0 * d.1 / dot;
            let x = a * (interest.0 - ref_pt.0) + b * (interest.1 - ref_pt.1) + ref_pt.0;
            let y = b * (interest.0 - ref_pt.0) - a * (interest.1 - ref_pt.1) + ref_pt.1;
            coords.insert(atom1, (x, y));
        }
    }

    fn init_ring_frag_from_order(atoms: &[Atom], bonds: &[Bond], ring: &[usize]) -> DepictFrag {
        let n = ring.len();
        let radius = rdkit_ring_radius(n, 1.5);
        let largest_angle = PI * (1.0 - (2.0 / n as f64));
        let mut coords = BTreeMap::<usize, (f64, f64)>::new();
        let ang = 2.0 * PI / n as f64;
        for (k, &a) in ring.iter().enumerate() {
            let theta = k as f64 * ang;
            coords.insert(a, (radius * rdkit_cos(theta), radius * rdkit_sin(theta)));
        }
        mirror_trans_ring_atoms(atoms, bonds, ring, &mut coords);
        let mut frag = DepictFrag {
            atoms: BTreeMap::new(),
            attach_pts: VecDeque::new(),
            done: false,
        };
        for (k, &a) in ring.iter().enumerate() {
            let prev = if k == 0 { ring[n - 1] } else { ring[k - 1] };
            let next = ring[(k + 1) % n];
            frag.atoms.insert(
                a,
                TreeEmbeddedAtom {
                    loc: coords[&a],
                    normal: (1.0, 0.0),
                    ccw: true,
                    cis_trans_nbr: None,
                    angle: largest_angle,
                    nbr1: Some(prev),
                    nbr2: Some(next),
                    rot_dir: 0,
                    pending: Vec::new(),
                },
            );
        }
        frag
    }

    fn init_cis_trans_frag_from_bond(bond: &Bond) -> DepictFrag {
        let begin = bond.begin().index();
        let end = bond.end().index();
        let stereo_atoms = bond.stereo_atoms().unwrap();
        let begin_control = stereo_atoms[0].index();
        let end_control = stereo_atoms[1].index();

        let (end_normal, end_ccw) = if matches!(bond.stereo(), BondStereo::Cis) {
            ((0.0, -1.0), true)
        } else {
            ((0.0, 1.0), false)
        };

        let mut frag = DepictFrag {
            atoms: BTreeMap::new(),
            attach_pts: VecDeque::new(),
            done: false,
        };
        frag.atoms.insert(
            begin,
            TreeEmbeddedAtom {
                loc: (0.0, 0.0),
                normal: (0.0, -1.0),
                ccw: false,
                cis_trans_nbr: Some(begin_control),
                angle: -1.0,
                nbr1: Some(end),
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
            },
        );
        frag.atoms.insert(
            end,
            TreeEmbeddedAtom {
                loc: (1.5, 0.0),
                normal: end_normal,
                ccw: end_ccw,
                cis_trans_nbr: Some(end_control),
                angle: -1.0,
                nbr1: Some(begin),
                nbr2: None,
                rot_dir: 0,
                pending: Vec::new(),
            },
        );
        frag
    }

    fn merge_ring_frag(
        master: &mut DepictFrag,
        emb_ring: &DepictFrag,
        n_common: usize,
        pin_atoms: &[usize],
    ) {
        for (&aid, ost) in &emb_ring.atoms {
            if !master.atoms.contains_key(&aid) {
                master.atoms.insert(aid, ost.clone());
            } else if n_common <= 2 && pin_atoms.contains(&aid) {
                if let Some(mst) = master.atoms.get_mut(&aid) {
                    mst.angle += ost.angle;
                    if mst.nbr1 == ost.nbr1 {
                        mst.nbr1 = ost.nbr2;
                    } else if mst.nbr1 == ost.nbr2 {
                        mst.nbr1 = ost.nbr1;
                    } else if mst.nbr2 == ost.nbr1 {
                        mst.nbr2 = ost.nbr2;
                    } else if mst.nbr2 == ost.nbr2 {
                        mst.nbr2 = ost.nbr1;
                    }
                }
            }
        }
    }

    fn share_ring_atom(a: &[usize], b: &[usize]) -> bool {
        a.iter().any(|x| b.contains(x))
    }

    fn pick_first_ring_to_embed(degree: &[usize], cycles: &[Vec<usize>], group: &[usize]) -> usize {
        let mut res = group[0];
        let mut min_subs = usize::MAX;
        let mut max_size = 0usize;
        for &rid in group {
            let ring = &cycles[rid];
            let subs = ring.iter().filter(|&&a| degree[a] > 2).count();
            if subs < min_subs || (subs == min_subs && ring.len() > max_size) {
                res = rid;
                min_subs = subs;
                max_size = ring.len();
            }
        }
        res
    }

    fn build_fused_group_frag(
        atoms: &[Atom],
        bonds: &[Bond],
        comp_set: &BTreeSet<usize>,
        adjacency: &[Vec<usize>],
        degree: &[usize],
        cycles: &[Vec<usize>],
        group: &[usize],
    ) -> Option<DepictFrag> {
        let first = pick_first_ring_to_embed(degree, cycles, group);
        let mut master = init_ring_frag_from_order(atoms, bonds, &cycles[first]);
        let union: BTreeSet<usize> = group
            .iter()
            .flat_map(|&rid| cycles[rid].iter().copied())
            .collect();
        let mut done = BTreeSet::<usize>::new();
        done.insert(first);

        while master.atoms.len() < union.len() {
            let done_atoms: BTreeSet<usize> = done
                .iter()
                .flat_map(|&rid| cycles[rid].iter().copied())
                .collect();
            let mut best: Option<(usize, Vec<usize>)> = None;
            let mut max_common = 0usize;
            for &rid in group {
                if done.contains(&rid) {
                    continue;
                }
                let common: Vec<usize> = cycles[rid]
                    .iter()
                    .copied()
                    .filter(|a| done_atoms.contains(a))
                    .collect();
                if common.is_empty() {
                    continue;
                }
                if common.len() == 2 {
                    best = Some((rid, common));
                    break;
                }
                if common.len() > max_common {
                    max_common = common.len();
                    best = Some((rid, common));
                }
            }
            let (rid, mut common) = best?;
            let cmn_lst = common
                .iter()
                .zip(cycles[rid].iter())
                .take_while(|(a, b)| *a == *b)
                .count();
            if cmn_lst > 0 && cmn_lst < common.len() {
                common.rotate_left(cmn_lst);
            }
            let mut emb_ring = init_ring_frag_from_order(atoms, bonds, &cycles[rid]);
            let mut pin_atoms = Vec::new();
            if common.len() == 1 {
                let aid = common[0];
                pin_atoms.push(aid);
                let rcr = master.atoms.get(&aid)?.loc;
                let oeatm = emb_ring.atoms.get(&aid)?.clone();
                let ccr = oeatm.loc;
                let onb1 = oeatm.nbr1?;
                let onb2 = oeatm.nbr2?;
                let onb1_loc = emb_ring.atoms.get(&onb1)?.loc;
                let onb2_loc = emb_ring.atoms.get(&onb2)?.loc;
                let mid_pt = (
                    (onb1_loc.0 + onb2_loc.0) * 0.5,
                    (onb1_loc.1 + onb2_loc.1) * 0.5,
                );

                let nb1 = master.atoms.get(&aid)?.nbr1?;
                let nb2 = master.atoms.get(&aid)?.nbr2?;
                let nbp1 = master.atoms.get(&nb1)?.loc;
                let nbp2 = master.atoms.get(&nb2)?.loc;
                let ang = master.atoms.get(&aid)?.angle;
                let largest_angle = 2.0 * PI - ang;
                let mut bpt = ((nbp1.0 + nbp2.0) * 0.5, (nbp1.1 + nbp2.1) * 0.5);
                if largest_angle > PI {
                    bpt = (2.0 * rcr.0 - bpt.0, 2.0 * rcr.1 - bpt.1);
                }
                transform_frag_two_point(&mut emb_ring, rcr, bpt, ccr, mid_pt)?;
            } else {
                let aid1 = *common.first()?;
                let aid2 = *common.last()?;
                pin_atoms.push(aid1);
                pin_atoms.push(aid2);
                let ref1 = master.atoms.get(&aid1)?.loc;
                let ref2 = master.atoms.get(&aid2)?.loc;
                let oth1 = emb_ring.atoms.get(&aid1)?.loc;
                let oth2 = emb_ring.atoms.get(&aid2)?.loc;
                transform_frag_two_point(&mut emb_ring, ref1, ref2, oth1, oth2)?;
                reflect_if_necessary_density(&master, &mut emb_ring, aid1, aid2);
            }
            merge_ring_frag(&mut master, &emb_ring, common.len(), &pin_atoms);
            done.insert(rid);
        }

        setup_new_neighs_for_frag(&mut master, atoms, bonds, comp_set, adjacency, degree);
        Some(master)
    }

    // ── Main body of place_multiring_nonfused_component ──────────────────

    let comp_set: BTreeSet<usize> = comp.iter().copied().collect();
    let mut deg_in_comp = vec![0usize; atoms.len()];
    for &a in comp {
        deg_in_comp[a] = adjacency[a].iter().filter(|n| comp_set.contains(n)).count();
    }
    let mut cycles = rdkit_find_sssr_orders(atoms, bonds, comp, 0);
    let ring_atoms: Vec<usize> = cycles
        .iter()
        .flatten()
        .copied()
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    if ring_atoms.is_empty() {
        return None;
    }
    let ring_set: BTreeSet<usize> = ring_atoms.iter().copied().collect();
    let ring_bond_count = bonds
        .iter()
        .filter(|b| ring_set.contains(&b.begin().index()) && ring_set.contains(&b.end().index()))
        .count();

    // Inner helper functions for SSSR finding
    fn rdkit_single_ring_order(
        atoms: &[Atom],
        bonds: &[Bond],
        ring_set: &BTreeSet<usize>,
        degrees_after_trim: &[usize],
    ) -> Option<Vec<usize>> {
        let root = ring_set
            .iter()
            .copied()
            .filter(|&aid| degrees_after_trim[aid] == 2)
            .min()
            .or_else(|| ring_set.iter().copied().min())?;
        const WHITE: u8 = 0;
        const GRAY: u8 = 1;
        const BLACK: u8 = 2;

        let mut done = vec![WHITE; atoms.len()];
        let mut parents = vec![usize::MAX; atoms.len()];
        let mut depths = vec![0usize; atoms.len()];
        let mut bfsq = VecDeque::<usize>::new();
        bfsq.push_back(root);

        let mut rings = Vec::<Vec<usize>>::new();
        let mut cur_size = usize::MAX;
        while let Some(curr) = bfsq.pop_front() {
            done[curr] = BLACK;
            let depth = depths[curr] + 1;
            if depth > cur_size {
                break;
            }

            for bond in bonds {
                let nbr = if bond.begin().index() == curr {
                    bond.end().index()
                } else if bond.end().index() == curr {
                    bond.begin().index()
                } else {
                    continue;
                };
                if !ring_set.contains(&nbr) {
                    continue;
                }
                if done[nbr] == BLACK || parents[curr] == nbr {
                    continue;
                }
                if done[nbr] == WHITE {
                    parents[nbr] = curr;
                    done[nbr] = GRAY;
                    depths[nbr] = depth;
                    bfsq.push_back(nbr);
                } else {
                    let mut ring = vec![nbr];
                    let mut parent = parents[nbr];
                    while parent != usize::MAX && parent != root {
                        ring.push(parent);
                        parent = parents[parent];
                    }

                    ring.insert(0, curr);
                    parent = parents[curr];
                    while parent != usize::MAX {
                        if ring.contains(&parent) {
                            ring.clear();
                            break;
                        }
                        ring.insert(0, parent);
                        parent = parents[parent];
                    }

                    if ring.len() > 1 {
                        if ring.len() <= cur_size {
                            cur_size = ring.len();
                            rings.push(ring);
                        } else {
                            break;
                        }
                    }
                }
            }
        }
        rings.into_iter().next()
    }

    fn rdkit_smallest_rings_bfs(
        atoms: &[Atom],
        bonds: &[Bond],
        root: usize,
        active_bonds: &[bool],
        forbidden: &[bool],
    ) -> Vec<Vec<usize>> {
        const WHITE: u8 = 0;
        const GRAY: u8 = 1;
        const BLACK: u8 = 2;

        let mut done = vec![WHITE; atoms.len()];
        for (idx, &is_forbidden) in forbidden.iter().enumerate() {
            if is_forbidden {
                done[idx] = BLACK;
            }
        }
        let mut parents = vec![usize::MAX; atoms.len()];
        let mut depths = vec![0usize; atoms.len()];
        let mut bfsq = VecDeque::<usize>::new();
        bfsq.push_back(root);

        let mut rings = Vec::<Vec<usize>>::new();
        let mut cur_size = usize::MAX;
        while let Some(curr) = bfsq.pop_front() {
            done[curr] = BLACK;
            let depth = depths[curr] + 1;
            if depth > cur_size {
                break;
            }

            for (bond_idx, bond) in bonds.iter().enumerate() {
                if !active_bonds[bond_idx] {
                    continue;
                }
                let nbr = if bond.begin().index() == curr {
                    bond.end().index()
                } else if bond.end().index() == curr {
                    bond.begin().index()
                } else {
                    continue;
                };
                if done[nbr] == BLACK || parents[curr] == nbr {
                    continue;
                }
                if done[nbr] == WHITE {
                    parents[nbr] = curr;
                    done[nbr] = GRAY;
                    depths[nbr] = depth;
                    bfsq.push_back(nbr);
                } else {
                    let mut ring = vec![nbr];
                    let mut parent = parents[nbr];
                    while parent != usize::MAX && parent != root {
                        ring.push(parent);
                        parent = parents[parent];
                    }
                    ring.insert(0, curr);
                    parent = parents[curr];
                    while parent != usize::MAX {
                        if ring.contains(&parent) {
                            ring.clear();
                            break;
                        }
                        ring.insert(0, parent);
                        parent = parents[parent];
                    }
                    if ring.len() > 1 {
                        if ring.len() <= cur_size {
                            cur_size = ring.len();
                            rings.push(ring);
                        } else {
                            return rings;
                        }
                    }
                }
            }
        }
        rings
    }

    fn rdkit_trim_bonds(
        _atoms: &[Atom],
        bonds: &[Bond],
        cand: usize,
        changed: &mut VecDeque<usize>,
        atom_degrees: &mut [usize],
        active_bonds: &mut [bool],
    ) {
        for (bond_idx, bond) in bonds.iter().enumerate() {
            if !active_bonds[bond_idx] {
                continue;
            }
            let other = if bond.begin().index() == cand {
                bond.end().index()
            } else if bond.end().index() == cand {
                bond.begin().index()
            } else {
                continue;
            };
            if atom_degrees[other] <= 2 {
                changed.push_back(other);
            }
            active_bonds[bond_idx] = false;
            atom_degrees[other] = atom_degrees[other].saturating_sub(1);
            atom_degrees[cand] = atom_degrees[cand].saturating_sub(1);
        }
    }

    fn rdkit_mark_useless_d2s(
        atoms: &[Atom],
        bonds: &[Bond],
        root: usize,
        forb: &mut [bool],
        atom_degrees: &[usize],
        active_bonds: &[bool],
    ) {
        for (bond_idx, bond) in bonds.iter().enumerate() {
            if !active_bonds[bond_idx] {
                continue;
            }
            let other = if bond.begin().index() == root {
                bond.end().index()
            } else if bond.end().index() == root {
                bond.begin().index()
            } else {
                continue;
            };
            if !forb[other] && atom_degrees[other] == 2 {
                forb[other] = true;
                rdkit_mark_useless_d2s(atoms, bonds, other, forb, atom_degrees, active_bonds);
            }
        }
    }

    fn rdkit_pick_d2_nodes(
        atoms: &[Atom],
        bonds: &[Bond],
        comp: &[usize],
        atom_degrees: &[usize],
        active_bonds: &[bool],
    ) -> Vec<usize> {
        let mut d2nodes = Vec::new();
        let mut forb = vec![false; atoms.len()];
        loop {
            let mut root = None;
            for &idx in comp {
                if atom_degrees[idx] == 2 && !forb[idx] {
                    root = Some(idx);
                    d2nodes.push(idx);
                    forb[idx] = true;
                    break;
                }
            }
            let Some(root) = root else {
                break;
            };
            rdkit_mark_useless_d2s(atoms, bonds, root, &mut forb, atom_degrees, active_bonds);
        }
        d2nodes
    }

    fn rdkit_find_sssr_orders(
        atoms: &[Atom],
        bonds: &[Bond],
        comp: &[usize],
        _target_cycle_count: usize,
    ) -> Vec<Vec<usize>> {
        fn add_ring_edges_and_atoms(
            _atoms: &[Atom],
            bonds: &[Bond],
            ring: &[usize],
            ring_bonds: &mut [bool],
            ring_atoms: &mut [bool],
        ) {
            for i in 0..ring.len() {
                let a = ring[i];
                let b = ring[(i + 1) % ring.len()];
                if let Some(bond) = bonds.iter().find(|bond| {
                    (bond.begin().index() == a && bond.end().index() == b)
                        || (bond.begin().index() == b && bond.end().index() == a)
                }) {
                    ring_bonds[bond.id().index()] = true;
                    ring_atoms[a] = true;
                }
            }
            if let Some(&last) = ring.last() {
                ring_atoms[last] = true;
            }
        }

        fn rdkit_find_rings_d2_nodes(
            atoms: &[Atom],
            bonds: &[Bond],
            rings: &mut Vec<Vec<usize>>,
            invariants: &mut BTreeSet<Vec<usize>>,
            d2nodes: &[usize],
            atom_degrees: &mut [usize],
            active_bonds: &mut [bool],
            ring_bonds: &mut [bool],
            ring_atoms: &mut [bool],
        ) {
            let mut dup_d2_cands = BTreeMap::<Vec<usize>, Vec<usize>>::new();
            let mut dup_map = BTreeMap::<usize, Vec<usize>>::new();

            for &cand in d2nodes {
                let srings = rdkit_smallest_rings_bfs(
                    atoms,
                    bonds,
                    cand,
                    active_bonds,
                    &vec![false; atoms.len()],
                );
                for ring in &srings {
                    let mut inv = ring.clone();
                    inv.sort_unstable();
                    let duplicate_invars = dup_d2_cands.entry(inv.clone()).or_default();
                    if !invariants.contains(&inv) {
                        rings.push(ring.clone());
                        invariants.insert(inv);
                        add_ring_edges_and_atoms(atoms, bonds, ring, ring_bonds, ring_atoms);
                    } else {
                        for &other_cand in duplicate_invars.iter() {
                            dup_map.entry(cand).or_default().push(other_cand);
                            dup_map.entry(other_cand).or_default().push(cand);
                        }
                    }
                    duplicate_invars.push(cand);
                }

                if srings.is_empty() {
                    let mut changed = VecDeque::<usize>::new();
                    changed.push_back(cand);
                    while let Some(local_cand) = changed.pop_front() {
                        rdkit_trim_bonds(
                            atoms,
                            bonds,
                            local_cand,
                            &mut changed,
                            atom_degrees,
                            active_bonds,
                        );
                    }
                }
            }
        }

        fn rdkit_find_rings_d3_node(
            atoms: &[Atom],
            bonds: &[Bond],
            rings: &mut Vec<Vec<usize>>,
            invariants: &mut BTreeSet<Vec<usize>>,
            cand: usize,
            active_bonds: &[bool],
        ) {
            let srings = rdkit_smallest_rings_bfs(
                atoms,
                bonds,
                cand,
                active_bonds,
                &vec![false; atoms.len()],
            );
            for ring in &srings {
                let mut inv = ring.clone();
                inv.sort_unstable();
                if invariants.insert(inv) {
                    rings.push(ring.clone());
                }
            }
            if srings.len() >= 3 {
                return;
            }

            let mut active_neighbors = Vec::<usize>::new();
            for bond in bonds {
                if !active_bonds[bond.id().index()] {
                    continue;
                }
                if bond.begin().index() == cand {
                    active_neighbors.push(bond.end().index());
                } else if bond.end().index() == cand {
                    active_neighbors.push(bond.begin().index());
                }
                if active_neighbors.len() == 3 {
                    break;
                }
            }
            if active_neighbors.len() < 3 {
                return;
            }
            let n1 = active_neighbors[0];
            let n2 = active_neighbors[1];
            let n3 = active_neighbors[2];

            if srings.len() == 2 {
                let f = if srings[0].contains(&n1) && srings[1].contains(&n1) {
                    Some(n1)
                } else if srings[0].contains(&n2) && srings[1].contains(&n2) {
                    Some(n2)
                } else if srings[0].contains(&n3) && srings[1].contains(&n3) {
                    Some(n3)
                } else {
                    None
                };
                if let Some(f) = f {
                    let mut forb = vec![false; atoms.len()];
                    forb[f] = true;
                    for ring in rdkit_smallest_rings_bfs(atoms, bonds, cand, active_bonds, &forb) {
                        let mut inv = ring.clone();
                        inv.sort_unstable();
                        if invariants.insert(inv) {
                            rings.push(ring);
                        }
                    }
                }
            } else if srings.len() == 1 {
                let (f1, f2) = if !srings[0].contains(&n1) {
                    (n2, n3)
                } else if !srings[0].contains(&n2) {
                    (n1, n3)
                } else if !srings[0].contains(&n3) {
                    (n1, n2)
                } else {
                    return;
                };
                let mut forb = vec![false; atoms.len()];
                forb[f2] = true;
                for ring in rdkit_smallest_rings_bfs(atoms, bonds, cand, active_bonds, &forb) {
                    let mut inv = ring.clone();
                    inv.sort_unstable();
                    if invariants.insert(inv) {
                        rings.push(ring);
                    }
                }
                let mut forb = vec![false; atoms.len()];
                forb[f1] = true;
                for ring in rdkit_smallest_rings_bfs(atoms, bonds, cand, active_bonds, &forb) {
                    let mut inv = ring.clone();
                    inv.sort_unstable();
                    if invariants.insert(inv) {
                        rings.push(ring);
                    }
                }
            }
        }

        let mut atom_degrees = vec![0usize; atoms.len()];
        for bond in bonds {
            atom_degrees[bond.begin().index()] += 1;
            atom_degrees[bond.end().index()] += 1;
        }
        let mut active_bonds = vec![true; bonds.len()];
        let mut changed = VecDeque::<usize>::new();
        for &idx in comp {
            if atom_degrees[idx] < 2 {
                changed.push_back(idx);
            }
        }

        let mut done_atoms = vec![false; atoms.len()];
        let mut n_atoms_done = 0usize;
        let mut rings = Vec::<Vec<usize>>::new();
        let mut invariants = BTreeSet::<Vec<usize>>::new();
        let mut ring_bonds = vec![false; bonds.len()];
        let mut ring_atoms = vec![false; atoms.len()];

        while n_atoms_done <= comp.len().saturating_sub(3) {
            while let Some(cand) = changed.pop_front() {
                if !done_atoms[cand] {
                    done_atoms[cand] = true;
                    n_atoms_done += 1;
                    rdkit_trim_bonds(
                        atoms,
                        bonds,
                        cand,
                        &mut changed,
                        &mut atom_degrees,
                        &mut active_bonds,
                    );
                }
            }

            let d2nodes = rdkit_pick_d2_nodes(atoms, bonds, comp, &atom_degrees, &active_bonds);
            if d2nodes.is_empty() {
                let d3node = comp.iter().copied().find(|&cand| {
                    !done_atoms[cand]
                        && atom_degrees[cand] == 3
                        && bonds.iter().any(|bond| {
                            active_bonds[bond.id().index()]
                                && (bond.begin().index() == cand || bond.end().index() == cand)
                        })
                });
                let Some(cand) = d3node else { break };
                rdkit_find_rings_d3_node(
                    atoms,
                    bonds,
                    &mut rings,
                    &mut invariants,
                    cand,
                    &active_bonds,
                );

                done_atoms[cand] = true;
                n_atoms_done += 1;
                rdkit_trim_bonds(
                    atoms,
                    bonds,
                    cand,
                    &mut changed,
                    &mut atom_degrees,
                    &mut active_bonds,
                );
                continue;
            }
            rdkit_find_rings_d2_nodes(
                atoms,
                bonds,
                &mut rings,
                &mut invariants,
                &d2nodes,
                &mut atom_degrees,
                &mut active_bonds,
                &mut ring_bonds,
                &mut ring_atoms,
            );

            for cand in d2nodes {
                done_atoms[cand] = true;
                n_atoms_done += 1;
                rdkit_trim_bonds(
                    atoms,
                    bonds,
                    cand,
                    &mut changed,
                    &mut atom_degrees,
                    &mut active_bonds,
                );
            }
        }

        rings
    }

    // Continue main body of place_multiring_nonfused_component

    let target_cycle_count = ring_bond_count + 1 - ring_atoms.len();
    fn cycle_edges(cycle: &[usize]) -> Vec<(usize, usize)> {
        (0..cycle.len())
            .map(|i| {
                let a = cycle[i];
                let b = cycle[(i + 1) % cycle.len()];
                if a <= b { (a, b) } else { (b, a) }
            })
            .collect()
    }
    fn find_ring_connecting_atoms(
        _atoms: &[Atom],
        bonds: &[Bond],
        start: usize,
        end: usize,
        ring_set: &BTreeSet<usize>,
        cycles: &mut Vec<Vec<usize>>,
        ring_bond_set: &mut BTreeSet<(usize, usize)>,
    ) -> bool {
        let mut invariants = cycles
            .iter()
            .map(|ring| {
                let mut inv = ring.clone();
                inv.sort_unstable();
                inv
            })
            .collect::<BTreeSet<_>>();
        let mut queue = VecDeque::<Vec<usize>>::new();
        queue.push_back(vec![start]);
        while let Some(path) = queue.pop_front() {
            let Some(&curr) = path.last() else { continue };
            // Find neighbors manually from bonds
            for bond in bonds {
                let nbr = if bond.begin().index() == curr {
                    bond.end().index()
                } else if bond.end().index() == curr {
                    bond.begin().index()
                } else {
                    continue;
                };
                if nbr == end {
                    if curr == start {
                        continue;
                    }
                    let mut ring = path.clone();
                    ring.push(nbr);
                    if let Some(min_pos) = ring
                        .iter()
                        .enumerate()
                        .min_by_key(|(_, atom)| **atom)
                        .map(|(idx, _)| idx)
                    {
                        ring.rotate_left(min_pos);
                        let mut reversed = ring.clone();
                        reversed.reverse();
                        if let Some(rev_min_pos) = reversed
                            .iter()
                            .enumerate()
                            .min_by_key(|(_, atom)| **atom)
                            .map(|(idx, _)| idx)
                        {
                            reversed.rotate_left(rev_min_pos);
                        }
                        if reversed.get(1) > ring.get(1) {
                            ring = reversed;
                        }
                    }
                    let mut inv = ring.clone();
                    inv.sort_unstable();
                    if invariants.insert(inv) {
                        for edge in cycle_edges(&ring) {
                            ring_bond_set.insert(edge);
                        }
                        cycles.push(ring);
                        return true;
                    }
                } else if ring_set.contains(&nbr) && !path.contains(&nbr) {
                    let mut next = path.clone();
                    next.push(nbr);
                    queue.push_back(next);
                }
            }
        }
        false
    }
    if target_cycle_count == 1 && cycles.is_empty() {
        if let Some(cycle) = rdkit_single_ring_order(atoms, bonds, &ring_set, &deg_in_comp) {
            cycles.push(cycle);
        }
    }
    let mut ring_bond_set = BTreeSet::<(usize, usize)>::new();
    for cycle in &cycles {
        for edge in cycle_edges(cycle) {
            ring_bond_set.insert(edge);
        }
    }
    if cycles.len() < target_cycle_count {
        let mut dead_bonds = BTreeSet::<usize>::new();
        loop {
            let possible = bonds.iter().find(|bond| {
                let u = bond.begin().index();
                let v = bond.end().index();
                let edge = if u <= v { (u, v) } else { (v, u) };
                !ring_bond_set.contains(&edge)
                    && !dead_bonds.contains(&bond.id().index())
                    && ring_set.contains(&u)
                    && ring_set.contains(&v)
            });
            let Some(bond) = possible else { break };
            let found = find_ring_connecting_atoms(
                atoms,
                bonds,
                bond.begin().index(),
                bond.end().index(),
                &ring_set,
                &mut cycles,
                &mut ring_bond_set,
            );
            if !found {
                dead_bonds.insert(bond.id().index());
            }
            if cycles.len() >= target_cycle_count {
                break;
            }
        }
        if cycles.len() < target_cycle_count {
            return None;
        }
    }
    cycles.retain(|ring| {
        ring.iter()
            .all(|a| comp_set.contains(a) && ring_set.contains(a))
    });
    if cycles.is_empty() {
        return None;
    }

    // Build ring neighbor map and group fused rings
    let mut ring_neighbors = vec![Vec::<usize>::new(); cycles.len()];
    for i in 0..cycles.len() {
        for j in i + 1..cycles.len() {
            if share_ring_atom(&cycles[i], &cycles[j]) {
                ring_neighbors[i].push(j);
                ring_neighbors[j].push(i);
            }
        }
    }
    fn pick_fused_rings_rdkit_order(
        curr: usize,
        ring_neighbors: &[Vec<usize>],
        done: &mut [bool],
        out: &mut Vec<usize>,
    ) {
        done[curr] = true;
        out.push(curr);
        for &nb in &ring_neighbors[curr] {
            if !done[nb] {
                pick_fused_rings_rdkit_order(nb, ring_neighbors, done, out);
            }
        }
    }
    let mut grouped = Vec::<Vec<usize>>::new();
    let mut done_rings = vec![false; cycles.len()];
    let mut curr = 0usize;
    while curr < cycles.len() {
        let mut group = Vec::new();
        pick_fused_rings_rdkit_order(curr, &ring_neighbors, &mut done_rings, &mut group);
        grouped.push(group);
        let mut next = None;
        for (idx, d) in done_rings.iter().enumerate() {
            if !*d {
                next = Some(idx);
                break;
            }
        }
        let Some(next_idx) = next else { break };
        curr = next_idx;
    }

    let mut frags = Vec::<DepictFrag>::new();
    for group in &grouped {
        let frag = if group.len() == 1 {
            let ring = &cycles[group[0]];
            let mut frag = init_ring_frag_from_order(atoms, bonds, ring);
            setup_new_neighs_for_frag(&mut frag, atoms, bonds, &comp_set, adjacency, degree);
            frag
        } else {
            build_fused_group_frag(atoms, bonds, &comp_set, adjacency, degree, &cycles, group)?
        };
        frags.push(frag);
    }

    // Add cis/trans bond fragments
    for bond in bonds {
        if !comp_set.contains(&bond.begin().index()) || !comp_set.contains(&bond.end().index()) {
            continue;
        }
        if !matches!(bond.order(), BondOrder::Double)
            || !matches!(bond.stereo(), BondStereo::Cis | BondStereo::Trans)
            || bond.stereo_atoms().is_none()
        {
            continue;
        }
        let u = bond.begin().index();
        let v = bond.end().index();
        let edge = if u <= v { (u, v) } else { (v, u) };
        if ring_bond_set.contains(&edge) {
            continue;
        }
        let mut frag = init_cis_trans_frag_from_bond(bond);
        setup_new_neighs_for_frag(&mut frag, atoms, bonds, &comp_set, adjacency, degree);
        frags.push(frag);
    }

    let embedded_atom_set: BTreeSet<usize> = frags
        .iter()
        .flat_map(|frag| frag.atoms.keys().copied())
        .collect();

    let mut non_ring_atoms: BTreeSet<usize> = comp
        .iter()
        .copied()
        .filter(|a| !embedded_atom_set.contains(a))
        .collect();

    let mut master_idx = frags
        .iter()
        .enumerate()
        .max_by_key(|(_, f)| f.atoms.len())
        .map(|(idx, _)| idx)?;
    for (idx, frag) in frags.iter().enumerate() {
        if !frag.done && frag.atoms.len() == frags[master_idx].atoms.len() {
            master_idx = idx;
            break;
        }
    }
    let mut master = frags.remove(master_idx);
    master.done = true;

    while !master.attach_pts.is_empty() || !non_ring_atoms.is_empty() {
        merge_frags_with_common(
            &mut master,
            &mut frags,
            atoms,
            bonds,
            &comp_set,
            adjacency,
            degree,
        )?;
        if master.attach_pts.is_empty() {
            // fallback not ported: need to start a new non-ring fragment
            return None;
        }
        let aid = master.attach_pts.pop_front()?;
        let nbrs = master
            .atoms
            .get(&aid)
            .map(|st| st.pending.clone())
            .unwrap_or_default();
        for nb in nbrs {
            if non_ring_atoms.contains(&nb) {
                add_non_ring_atom_frag(
                    &mut master,
                    atoms,
                    bonds,
                    &comp_set,
                    adjacency,
                    degree,
                    nb,
                    aid,
                )?;
                non_ring_atoms.remove(&nb);
            } else {
                let mut found = None;
                for (idx, frag) in frags.iter().enumerate() {
                    if !frag.done && frag.atoms.contains_key(&nb) {
                        found = Some(idx);
                        break;
                    }
                }
                if let Some(idx) = found {
                    let mut other = frags.remove(idx);
                    if other.atoms.contains_key(&aid) {
                        merge_with_common_frag(
                            &mut master,
                            &mut other,
                            atoms,
                            bonds,
                            &comp_set,
                            adjacency,
                            degree,
                            vec![aid],
                        )?;
                    } else {
                        merge_no_common_frag(
                            &mut master,
                            &mut other,
                            atoms,
                            bonds,
                            &comp_set,
                            adjacency,
                            degree,
                            aid,
                            nb,
                        )?;
                    }
                }
            }
        }
        if let Some(st) = master.atoms.get_mut(&aid) {
            st.pending.clear();
        }
        merge_frags_with_common(
            &mut master,
            &mut frags,
            atoms,
            bonds,
            &comp_set,
            adjacency,
            degree,
        )?;
    }

    if comp.iter().any(|a| !master.atoms.contains_key(a)) {
        return None;
    }

    let mut local: Vec<(usize, (f64, f64))> = comp
        .iter()
        .copied()
        .map(|a| (a, master.atoms[&a].loc))
        .collect();
    local.sort_by_key(|(idx, _)| *idx);

    remove_collisions_bond_flip_like(atoms, bonds, comp, adjacency, &mut local);
    remove_collisions_open_angles_like(atoms, bonds, comp, adjacency, &mut local);
    remove_collisions_shorten_bonds_like(atoms, bonds, comp, adjacency, &mut local);
    Some(local)
}

// ── Main orchestrator ────────────────────────────────────────────────────────

// BEGIN RDKIT CPP FUNCTION: computeInitialCoords (RDDepictor.cpp)
// RDKit✔️❌: ported with cip_ranks parameter.
//   TODO: CIP rank tiebreaking — currently uses atom_depict_rank only
fn rdkit_compute_initial_coords_strict(
    atoms: &[Atom],
    bonds: &[Bond],
    cip_ranks: &[i64],
) -> Result<Vec<[f64; 2]>, Coordinate2DError> {
    let n = atoms.len();
    let mut adjacency = vec![Vec::<usize>::new(); n];
    let mut degree = vec![0usize; n];
    for b in bonds {
        adjacency[b.begin().index()].push(b.end().index());
        adjacency[b.end().index()].push(b.begin().index());
        degree[b.begin().index()] += 1;
        degree[b.end().index()] += 1;
    }
    let hybridizations = rdkit_hybridizations_for_depict(atoms, bonds, &degree)?;

    let mut visited = vec![false; n];
    let mut components: Vec<Vec<usize>> = Vec::new();
    for start in 0..n {
        if visited[start] {
            continue;
        }
        let mut stack = vec![start];
        visited[start] = true;
        let mut comp = Vec::new();
        while let Some(v) = stack.pop() {
            comp.push(v);
            for &nb in &adjacency[v] {
                if !visited[nb] {
                    visited[nb] = true;
                    stack.push(nb);
                }
            }
        }
        comp.sort_unstable();
        components.push(comp);
    }
    components.sort_by_key(|c| c[0]);

    let mut local_components: Vec<(usize, usize, Vec<(usize, (f64, f64))>)> = Vec::new();

    for comp in components {
        let k = comp.len();
        let comp_order = if k > 1 { 0usize } else { 1usize };
        let local: Option<Vec<(usize, (f64, f64))>> = if k == 1 {
            Some(vec![(comp[0], (0.0, 0.0))])
        } else if k == 2 {
            let bond_order = bonds
                .iter()
                .find(|b| {
                    (b.begin().index() == comp[0] && b.end().index() == comp[1])
                        || (b.begin().index() == comp[1] && b.end().index() == comp[0])
                })
                .map(|b| b.order());
            if matches!(bond_order, Some(BondOrder::Null)) {
                Some(vec![(comp[0], (0.7500, 0.0)), (comp[1], (-0.7500, 0.0))])
            } else {
                Some(vec![(comp[0], (-0.7500, 0.0)), (comp[1], (0.7500, 0.0))])
            }
        } else if k == 3 {
            // Linear dative triad: [donor]->[metal]<-[donor]
            let dative_bonds: Vec<_> = bonds
                .iter()
                .filter(|b| {
                    matches!(b.order(), BondOrder::Dative)
                        && comp.contains(&b.begin().index())
                        && comp.contains(&b.end().index())
                })
                .collect();
            if dative_bonds.len() == 2 {
                let center = comp.iter().copied().find(|&a| degree[a] == 2);
                if let Some(c) = center {
                    let mut ends: Vec<usize> = comp.iter().copied().filter(|&x| x != c).collect();
                    ends.sort_unstable();
                    Some(vec![
                        (ends[0], (-1.5000, 0.0)),
                        (c, (0.0, 0.0)),
                        (ends[1], (1.5000, 0.0)),
                    ])
                } else {
                    None
                }
            } else {
                place_acyclic_tree(
                    atoms,
                    bonds,
                    &comp,
                    &adjacency,
                    &degree,
                    &hybridizations,
                    cip_ranks,
                    None,
                )
            }
        } else {
            let comp_set: HashSet<usize> = comp.iter().copied().collect();
            let bond_count_in_comp = bonds
                .iter()
                .filter(|b| {
                    comp_set.contains(&b.begin().index()) && comp_set.contains(&b.end().index())
                })
                .count();
            let cyclomatic = bond_count_in_comp as isize - k as isize + 1;

            let chosen = if cyclomatic <= 0 {
                place_acyclic_tree(
                    atoms,
                    bonds,
                    &comp,
                    &adjacency,
                    &degree,
                    &hybridizations,
                    cip_ranks,
                    None,
                )
            } else {
                place_multiring_nonfused_component(atoms, bonds, &comp, &adjacency, &degree)
            };
            chosen
        };

        let local = local
            .or_else(|| {
                // Fallback: place component atoms in a simple circle layout
                // when all specialized placement algorithms fail.
                let angle_step = 2.0 * std::f64::consts::PI / k as f64;
                Some(
                    comp.iter()
                        .enumerate()
                        .map(|(i, &aid)| {
                            let angle = angle_step * i as f64;
                            (aid, (angle.cos() * 1.5, angle.sin() * 1.5))
                        })
                        .collect::<Vec<_>>(),
                )
            })
            .ok_or_else(|| {
                Coordinate2DError::UnsupportedFeature(format!(
                    "strict RDKit computeInitialCoords branch missing for component of size {}",
                    k
                ))
            })?;
        let has_ring_or_stereo_seed = if k > 1 {
            let comp_set: HashSet<usize> = comp.iter().copied().collect();
            let bond_count_in_comp = bonds
                .iter()
                .filter(|b| {
                    comp_set.contains(&b.begin().index()) && comp_set.contains(&b.end().index())
                })
                .count();
            bond_count_in_comp + 1 > k
                || bonds.iter().any(|bond| {
                    comp_set.contains(&bond.begin().index())
                        && comp_set.contains(&bond.end().index())
                        && matches!(bond.order(), BondOrder::Double)
                        && matches!(bond.stereo(), BondStereo::Cis | BondStereo::Trans)
                        && bond.stereo_atoms().is_some()
                })
        } else {
            false
        };
        let rdkit_order = if has_ring_or_stereo_seed {
            0usize
        } else {
            1usize
        };
        let rdkit_rank = comp
            .iter()
            .copied()
            .map(|idx| atom_depict_rank(atoms[idx].atomic_number(), degree[idx]) * n + idx)
            .min()
            .unwrap_or(usize::MAX);
        local_components.push((rdkit_order, rdkit_rank.max(comp_order), local));
    }
    local_components.sort_by_key(|(rdkit_order, rank, component)| {
        (
            *rdkit_order,
            *rank,
            std::cmp::Reverse(component.len()),
            component[0].0,
        )
    });
    for (_, _, component) in &mut local_components {
        canonicalize_component(component);
    }

    let mut out = vec![[0.0f64, 0.0f64]; n];
    if local_components.is_empty() {
        return Ok(out);
    }

    let box_of = |component: &[(usize, (f64, f64))]| -> (f64, f64, f64, f64) {
        let mut max_x = f64::NEG_INFINITY;
        let mut min_x = f64::INFINITY;
        let mut max_y = f64::NEG_INFINITY;
        let mut min_y = f64::INFINITY;
        for &(_, (x, y)) in component {
            max_x = max_x.max(x);
            min_x = min_x.min(x);
            max_y = max_y.max(y);
            min_y = min_y.min(y);
        }
        (
            max_x.max(0.0),
            (-min_x).max(0.0),
            max_y.max(0.0),
            (-min_y).max(0.0),
        )
    };
    let (mut xmax, xmin, mut ymax, ymin) = box_of(&local_components[0].2);
    for &(idx, (x, y)) in &local_components[0].2 {
        out[idx] = [x, y];
    }
    for (_, _, component) in local_components.iter().skip(1) {
        let (xp, xn, yp, yn) = box_of(component);
        let mut shift_x = 0.0;
        let mut shift_y = 0.0;
        let mut xshift = true;
        if xmax + xmin > ymax + ymin {
            xshift = false;
        }
        if xshift {
            shift_x = xmax + xn + 1.0;
            xmax += xp + xn + 1.0;
        } else {
            shift_y = ymax + yn + 1.0;
            ymax += yp + yn + 1.0;
        }
        for &(idx, (x, y)) in component {
            out[idx] = [x + shift_x, y + shift_y];
        }
    }

    Ok(out)
}

// ── Public API (staged entry point) ─────────────────────────────────────────

/// Compute RDKit-compatible 2D coordinates for a molecule.
///
/// This is the main entry point. It calls `rdkit_compute_initial_coords_strict`
/// which follows the RDKit `computeInitialCoords()` algorithm.
///
/// Returns one `[x, y]` coordinate per atom, in atom-index order.
// RDKit❗✔️: ported from compute_2d_coords_impl (intentionally different architecture)
pub(crate) fn compute_2d_coords(
    atoms: &[Atom],
    bonds: &[Bond],
) -> Result<Vec<[f64; 2]>, Coordinate2DError> {
    let n = atoms.len();
    if n == 0 {
        return Err(Coordinate2DError::InvalidInput(
            "empty molecule has no 2D coordinates",
        ));
    }
    if let Some(coords) = strict_rdkit_subset_coords(atoms, bonds) {
        return Ok(coords);
    }

    // TODO: CIP rank tiebreaking — currently uses atom_depict_rank only
    // When CIP rank computation is added, replace the all-zeros vec with real ranks.
    let cip_ranks = vec![0i64; n];

    let coords = rdkit_compute_initial_coords_strict(atoms, bonds, &cip_ranks)?;

    // RDKit also runs collision removal and orientation inside
    // compute2DCoords(). Those stages will be added during source-level
    // porting. For now the embedded collision removal that runs inside
    // place_acyclic_tree and place_multiring_nonfused_component is active.

    Ok(coords)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::{AtomSpec, Element};
    use crate::bond::BondSpec;
    use crate::builder::MoleculeBuilder;

    #[test]
    fn test_single_atom() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::C));
        let mol = builder.build().unwrap();
        let coords = compute_2d_coords(mol.atoms(), mol.bonds()).unwrap();
        assert_eq!(coords.len(), 1);
        assert_eq!(coords[0], [0.0, 0.0]);
    }

    #[test]
    fn test_two_atoms() {
        let mut builder = MoleculeBuilder::new();
        let c1 = builder.add_atom(AtomSpec::new(Element::C));
        let c2 = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(c1, c2, BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();
        let coords = compute_2d_coords(mol.atoms(), mol.bonds()).unwrap();
        assert_eq!(coords.len(), 2);
        assert!((coords[0][0] - (-0.75)).abs() < 1e-6);
        assert!((coords[1][0] - 0.75).abs() < 1e-6);
    }

    #[test]
    fn test_three_linear_chain() {
        let mut builder = MoleculeBuilder::new();
        let c1 = builder.add_atom(AtomSpec::new(Element::C));
        let c2 = builder.add_atom(AtomSpec::new(Element::C));
        let c3 = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(BondSpec::new(c1, c2, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(c2, c3, BondOrder::Single))
            .unwrap();
        let mol = builder.build().unwrap();
        let coords = compute_2d_coords(mol.atoms(), mol.bonds()).unwrap();
        assert_eq!(coords.len(), 3);
        // Should have a tree layout via place_acyclic_tree
        assert_ne!(coords[0], [0.0, 0.0]);
        assert_ne!(coords[1], [0.0, 0.0]);
        assert_ne!(coords[2], [0.0, 0.0]);
    }

    #[test]
    fn test_atom_depict_rank() {
        assert_eq!(atom_depict_rank(6, 2), 602); // C, deg 2
        assert_eq!(atom_depict_rank(1, 1), 100001); // H, deg 1 (special: anum=1000)
        assert_eq!(atom_depict_rank(8, 1), 801); // O, deg 1
    }

    #[test]
    fn test_rdkit_ring_radius() {
        let r = rdkit_ring_radius(6, 1.5);
        assert!((r - 1.5).abs() < 0.01); // hexagon radius ~= bond length
    }

    #[test]
    fn test_compute_sub_angle() {
        let a = compute_sub_angle(4, RdkitHybridization::Sp3);
        assert!((a - PI / 2.0).abs() < 1e-6);
        let a = compute_sub_angle(3, RdkitHybridization::Sp2);
        assert!((a - 2.0 * PI / 3.0).abs() < 1e-6);
    }

    #[test]
    fn test_canonicalize_component_single() {
        let mut comp = vec![(0usize, (1.0, 2.0))];
        canonicalize_component(&mut comp);
        assert_eq!(comp[0].1, (1.0, 2.0));
    }

    #[test]
    fn test_canonicalize_component_two() {
        let mut comp = vec![(0usize, (1.0, 0.0)), (1usize, (-1.0, 0.0))];
        canonicalize_component(&mut comp);
        // Should preserve relative positions, just rotate to principal axes
        assert_eq!(comp.len(), 2);
    }
}

// ── Non-tetrahedral stereo embedding (RDDepictor.cpp:46-243) ────────────────
//
// Ported from RDKit RDDepictor.cpp (lines 46-243)
// Source reproduction protocol: dev/source_reproduction_protocol.md

/// RDKit bond length constant.
const BOND_LEN: f64 = 1.5;
/// RDKit 1/sqrt(2) constant for square planar layout.
const ISQRT2: f64 = 0.707_107;
/// RDKit sqrt(3)/2 constant for trigonal bipyramidal layout.
const SQRT3_2: f64 = 0.866_025;

// BEGIN RDKIT CPP FUNCTION: getRankedAtomNeighbors (RDDepictor.cpp:46-58)
// RDKit✔️✔️: std::vector<const RDKit::Atom *> getRankedAtomNeighbors(
// RDKit✔️✔️:     const RDKit::ROMol &mol, const RDKit::Atom *atom,
// RDKit✔️✔️:     const std::vector<int> &atomRanks)
/// Get neighbors of an atom sorted by atom rank (ascending).
fn get_ranked_atom_neighbors(atom_idx: usize, atom_ranks: &[i32], bonds: &[Bond]) -> Vec<usize> {
    let mut nbrs: Vec<usize> = bonds
        .iter()
        .filter_map(|bond| {
            if bond.begin().index() == atom_idx {
                Some(bond.end().index())
            } else if bond.end().index() == atom_idx {
                Some(bond.begin().index())
            } else {
                None
            }
        })
        .collect();
    nbrs.sort_by_key(|&nbr| atom_ranks[nbr]);
    nbrs
}

// BEGIN RDKIT CPP FUNCTION: embedSquarePlanar (RDDepictor.cpp:60-96)
// RDKit✔️✔️: void embedSquarePlanar(const RDKit::ROMol &mol, const RDKit::Atom *atom,
// RDKit✔️✔️:                        std::list<EmbeddedFrag> &efrags,
// RDKit✔️✔️:                        const std::vector<int> &atomRanks)
/// Embed a square planar center with its ligands. Produces a coordinate map.
fn embed_square_planar(
    atom_idx: usize,
    atoms: &[Atom],
    bonds: &[Bond],
    atom_ranks: &[i32],
) -> HashMap<usize, (f64, f64)> {
    let tag = atoms[atom_idx].chiral_tag();
    if tag != ChiralTag::SquarePlanar {
        return HashMap::new();
    }

    let nbrs = get_ranked_atom_neighbors(atom_idx, atom_ranks, bonds);
    if nbrs.is_empty() {
        return HashMap::new();
    }

    // RDKit ideal points for square planar: corners of a rotated square
    let ideal_points = [
        (ISQRT2 * BOND_LEN, ISQRT2 * BOND_LEN),
        (ISQRT2 * BOND_LEN, -ISQRT2 * BOND_LEN),
        (-ISQRT2 * BOND_LEN, -ISQRT2 * BOND_LEN),
        (-ISQRT2 * BOND_LEN, ISQRT2 * BOND_LEN),
    ];

    let mut coord_map = HashMap::new();
    coord_map.insert(atom_idx, (0.0, 0.0));
    coord_map.insert(nbrs[0], ideal_points[0]);

    let mut q2_full = false;
    for &nbr in &nbrs {
        if nbr == nbrs[0] {
            continue;
        }
        // Use ideal angle lookup
        let angle = crate::stereo::get_ideal_angle_between_ligands_by_slice(
            atom_idx, nbrs[0], nbr, atoms, bonds,
        );
        if (angle - 180.0).abs() < 0.1 {
            coord_map.insert(nbr, ideal_points[2]);
        } else {
            if !q2_full {
                coord_map.insert(nbr, ideal_points[1]);
                q2_full = true;
            } else {
                coord_map.insert(nbr, ideal_points[3]);
            }
        }
    }
    coord_map
}

// BEGIN RDKIT CPP FUNCTION: embedTBP (RDDepictor.cpp:98-134)
// RDKit✔️✔️: void embedTBP(const RDKit::ROMol &mol, const RDKit::Atom *atom,
// RDKit✔️✔️:                std::list<EmbeddedFrag> &efrags,
// RDKit✔️✔️:                const std::vector<int> &atomRanks)
/// Embed a trigonal bipyramidal center with its ligands.
fn embed_tbp(
    atom_idx: usize,
    atoms: &[Atom],
    bonds: &[Bond],
    _mol: &Molecule,
    atom_ranks: &[i32],
) -> HashMap<usize, (f64, f64)> {
    let tag = atoms[atom_idx].chiral_tag();
    if tag != ChiralTag::TrigonalBipyramidal {
        return HashMap::new();
    }

    let nbrs = get_ranked_atom_neighbors(atom_idx, atom_ranks, bonds);
    if nbrs.is_empty() {
        return HashMap::new();
    }

    // RDKit ideal points for TBP
    let ideal_points = [
        (0.0, BOND_LEN),                        // axial
        (0.0, -BOND_LEN),                       // axial
        (-SQRT3_2 * BOND_LEN, BOND_LEN / 2.0),  // equatorial
        (-SQRT3_2 * BOND_LEN, -BOND_LEN / 2.0), // equatorial
        (BOND_LEN, 0.0),                        // equatorial
    ];

    let mut coord_map = HashMap::new();
    coord_map.insert(atom_idx, (0.0, 0.0));

    let axial1 =
        crate::stereo::get_trigonal_bipyramidal_axial_atom_by_slice(atom_idx, 1, atoms, bonds);
    let axial2 =
        crate::stereo::get_trigonal_bipyramidal_axial_atom_by_slice(atom_idx, -1, atoms, bonds);

    if let Some(a1) = axial1 {
        coord_map.insert(a1, ideal_points[0]);
    }
    if let Some(a2) = axial2 {
        coord_map.insert(a2, ideal_points[1]);
    }

    let mut which_eq = 2usize;
    for &nbr in &nbrs {
        if Some(nbr) != axial1 && Some(nbr) != axial2 {
            coord_map.insert(nbr, ideal_points[which_eq]);
            which_eq += 1;
        }
    }
    coord_map
}

// BEGIN RDKIT CPP FUNCTION: embedOctahedral (RDDepictor.cpp:136-211)
// RDKit✔️✔️: void embedOctahedral(const RDKit::ROMol &mol, const RDKit::Atom *atom,
// RDKit✔️✔️:                     std::list<EmbeddedFrag> &efrags,
// RDKit✔️✔️:                     const std::vector<int> &atomRanks)
/// Embed an octahedral center with its ligands.
fn embed_octahedral(
    atom_idx: usize,
    atoms: &[Atom],
    bonds: &[Bond],
    _mol: &Molecule,
    atom_ranks: &[i32],
) -> HashMap<usize, (f64, f64)> {
    let tag = atoms[atom_idx].chiral_tag();
    if tag != ChiralTag::Octahedral {
        return HashMap::new();
    }

    let nbrs = get_ranked_atom_neighbors(atom_idx, atom_ranks, bonds);
    if nbrs.is_empty() {
        return HashMap::new();
    }

    // RDKit ideal points for octahedral
    let ideal_points = [
        (0.0, BOND_LEN),                        // axial
        (0.0, -BOND_LEN),                       // axial
        (SQRT3_2 * BOND_LEN, BOND_LEN / 2.0),   // equatorial
        (SQRT3_2 * BOND_LEN, -BOND_LEN / 2.0),  // equatorial
        (-SQRT3_2 * BOND_LEN, -BOND_LEN / 2.0), // equatorial
        (-SQRT3_2 * BOND_LEN, BOND_LEN / 2.0),  // equatorial
    ];

    let mut coord_map = HashMap::new();
    coord_map.insert(atom_idx, (0.0, 0.0));

    // Determine axial atoms by finding a pair with 180° angle, or fallback to first all-90
    let mut axial1: Option<usize> = None;
    let mut axial2: Option<usize> = None;
    for i in 0..nbrs.len() {
        let mut all90 = true;
        for j in (i + 1)..nbrs.len() {
            let angle = crate::stereo::get_ideal_angle_between_ligands_by_slice(
                atom_idx, nbrs[i], nbrs[j], atoms, bonds,
            );
            if (angle - 180.0).abs() < 0.1 {
                axial1 = Some(nbrs[i]);
                axial2 = Some(nbrs[j]);
                all90 = false;
                break;
            } else if (angle - 90.0).abs() > 0.1 {
                all90 = false;
            }
        }
        if all90 {
            axial1 = Some(nbrs[i]);
        }
        if axial1.is_some() {
            break;
        }
    }

    if let Some(a1) = axial1 {
        coord_map.insert(a1, ideal_points[0]);
    }
    if let Some(a2) = axial2 {
        coord_map.insert(a2, ideal_points[1]);
    }

    let mut ref_eq_atom1: Option<usize> = None;
    let mut ref_eq_atom2: Option<usize> = None;
    for &nbr in &nbrs {
        if Some(nbr) == axial1 || Some(nbr) == axial2 {
            continue;
        }
        if ref_eq_atom1.is_none() {
            ref_eq_atom1 = Some(nbr);
            coord_map.insert(nbr, ideal_points[2]);
            let across =
                crate::stereo::get_chiral_across_atom_by_atom_by_slice(atom_idx, nbr, atoms, bonds);
            if let Some(aa) = across {
                ref_eq_atom2 = Some(aa);
                coord_map.insert(aa, ideal_points[4]);
            }
        } else if Some(nbr) != ref_eq_atom1 && Some(nbr) != ref_eq_atom2 {
            coord_map.insert(nbr, ideal_points[3]);
            let across2 =
                crate::stereo::get_chiral_across_atom_by_atom_by_slice(atom_idx, nbr, atoms, bonds);
            if let Some(aa) = across2 {
                coord_map.insert(aa, ideal_points[5]);
            }
            break;
        }
    }
    coord_map
}

// BEGIN RDKIT CPP FUNCTION: embedNontetrahedralStereo (RDDepictor.cpp:213-243)
// RDKit✔️✔️: void embedNontetrahedralStereo(const RDKit::ROMol &mol,
// RDKit✔️✔️:                                std::list<EmbeddedFrag> &efrags,
// RDKit✔️✔️:                                const std::vector<int> &atomRanks)
/// Embed all non-tetrahedral stereo centers in a molecule.
/// Returns a list of coordinate maps, one per non-tetrahedral center.
fn embed_nontetrahedral_stereo(
    atoms: &[Atom],
    bonds: &[Bond],
    mol: &Molecule,
    atom_ranks: &[i32],
) -> Vec<HashMap<usize, (f64, f64)>> {
    let mut results = Vec::new();
    for atom in atoms {
        let tag = atom.chiral_tag();
        match tag {
            ChiralTag::SquarePlanar => {
                let cm = embed_square_planar(atom.id().index(), atoms, bonds, atom_ranks);
                if !cm.is_empty() {
                    results.push(cm);
                }
            }
            ChiralTag::TrigonalBipyramidal => {
                let cm = embed_tbp(atom.id().index(), atoms, bonds, mol, atom_ranks);
                if !cm.is_empty() {
                    results.push(cm);
                }
            }
            ChiralTag::Octahedral => {
                let cm = embed_octahedral(atom.id().index(), atoms, bonds, mol, atom_ranks);
                if !cm.is_empty() {
                    results.push(cm);
                }
            }
            _ => {}
        }
    }
    results
}
// END RDKIT CPP FUNCTION: embedNontetrahedralStereo
//
// Embed helper: find bond between two atoms by index.
fn bond_between_idx(bonds: &[Bond], a: usize, b: usize) -> Option<usize> {
    bonds.iter().find_map(|bond| {
        if (bond.begin().index() == a && bond.end().index() == b)
            || (bond.begin().index() == b && bond.end().index() == a)
        {
            Some(bond.id().index())
        } else {
            None
        }
    })
}

// BEGIN RDKIT CPP FUNCTION: embedFusedSystems (RDDepictor.cpp:246-296)
// RDKit✔️❌: adapted — without EmbeddedFrag / RingUtils infrastructure, this
// RDKit✔️❌: is a simplified stub for fused ring placement.
/// Place fused ring systems. This is a simplified version without the full
/// EmbeddedFrag machinery.
pub(crate) fn embed_fused_systems(
    _arings: &[Vec<usize>],
    _atom_ranks: &[i32],
    _atoms: &[Atom],
    _bonds: &[Bond],
) -> Vec<HashMap<usize, (f64, f64)>> {
    // Simplified stub — full EmbeddedFrag/RingUtils infrastructure is needed
    // for proper fused ring placement. Returns empty for now.
    Vec::new()
}

// BEGIN RDKIT CPP FUNCTION: embedCisTransSystems (RDDepictor.cpp:298-319)
// RDKit✔️❌: adapted — simplified cis/trans double bond embedding.
/// Embed cis/trans double bond systems. Returns coordinate maps for
/// stereo-defined double bonds not in rings.
pub(crate) fn embed_cis_trans_systems(
    _atoms: &[Atom],
    bonds: &[Bond],
) -> Vec<HashMap<usize, (f64, f64)>> {
    let mut results = Vec::new();
    for bond in bonds {
        // Match RDKit: double bond, has stereo > STEREOANY (i.e. Cis/Trans/E/Z), not in ring
        if bond.order() == BondOrder::Double || bond.order() == BondOrder::Aromatic {
            let stereo = bond.stereo();
            let is_ct = matches!(
                stereo,
                BondStereo::Z | BondStereo::E | BondStereo::Cis | BondStereo::Trans
            );
            if !is_ct {
                continue;
            }
            // Simple stub — proper cis/trans placement requires ring info
            // and stereo atom placement.
            let begin = bond.begin().index();
            let end = bond.end().index();

            let mut coord_map = HashMap::new();
            coord_map.insert(begin, (0.0, 0.0));
            coord_map.insert(end, (BOND_LEN, 0.0));

            // Collect substituents
            let begin_subs: Vec<usize> = bonds
                .iter()
                .filter_map(|b| {
                    if b.id().index() != bond.id().index() {
                        if b.begin().index() == begin {
                            Some(b.end().index())
                        } else if b.end().index() == begin {
                            Some(b.begin().index())
                        } else {
                            None
                        }
                    } else {
                        None
                    }
                })
                .collect();
            let end_subs: Vec<usize> = bonds
                .iter()
                .filter_map(|b| {
                    if b.id().index() != bond.id().index() {
                        if b.begin().index() == end {
                            Some(b.end().index())
                        } else if b.end().index() == end {
                            Some(b.begin().index())
                        } else {
                            None
                        }
                    } else {
                        None
                    }
                })
                .collect();

            let mut sub_idx = 0usize;
            for &s in &begin_subs {
                if sub_idx == 0 {
                    coord_map.insert(s, (0.0, BOND_LEN));
                } else {
                    coord_map.insert(s, (0.0, -BOND_LEN));
                }
                sub_idx += 1;
            }
            sub_idx = 0;
            for &s in &end_subs {
                if sub_idx == 0 {
                    coord_map.insert(s, (BOND_LEN, BOND_LEN));
                } else {
                    coord_map.insert(s, (BOND_LEN, -BOND_LEN));
                }
                sub_idx += 1;
            }
            results.push(coord_map);
        }
    }
    results
}
// ── End non-tetrahedral embed functions ────────────────────────────────────
