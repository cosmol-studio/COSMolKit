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

use crate::chemistry::forcefield::crystalff::{
    CrystalFFDetails, TorsionAngleContribs, get_experimental_torsions_without_bonds,
};
use crate::chemistry::forcefield::mmff::nonbonded::NonbondedContrib;
#[cfg(test)]
use crate::chemistry::forcefield::uff::atom_typer::get_atom_label_for_uff as source_get_atom_label_for_uff;
use crate::chemistry::forcefield::uff::atom_typer::get_atom_types_for_uff;
use crate::chemistry::forcefield::uff::inversion::InversionContribs;
use crate::chemistry::forcefield::uff::params::AtomicParams;
use crate::chemistry::forcefield::{
    AngleConstraintContribs, DistanceConstraintContribs, ForceField, ForceFieldContrib,
    ForceFieldVec3,
};
use crate::chemistry::stereo::{get_ideal_angle_between_ligands, has_non_tetrahedral_stereo};
use crate::molecule::CoordinateBlock as MoleculeCoordinateBlock;
use crate::read_parts::MoleculeReadParts;
use crate::{
    Atom, AtomId, Bond, BondId, BondOrder, BondQueryPredicate, BondStereo, ChiralTag, Conformer3D,
    DerivedState, Hybridization, Molecule, MoleculeBuildError, QueryNode, SubstructMatchParams,
    ValenceModel, assign_valence, get_substruct_matches_with_params, rdkit_valence_list,
};
use std::borrow::Cow;
use std::collections::{BTreeMap, HashSet};
use std::env;
use std::sync::{Arc, Mutex, OnceLock};
use std::time::{Duration, Instant};

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
const KNOWN_DIST_TOL: f64 = 0.01;
const KNOWN_DIST_FORCE_CONSTANT: f64 = 100.0;
const MIN_MACROCYCLE_RING_SIZE: usize = 9;
const DEFAULT_LOWER: f64 = 0.001;
const DEFAULT_UPPER: f64 = MAX_UPPER;
const UFF_LAMBDA: f64 = 0.1332;
const EIGVAL_TOL: f64 = 0.001;
const RDKIT_RANDOM_MODULUS: u64 = 2_147_483_647;
const RDKIT_RANDOM_MULTIPLIER: u64 = 48_271;
const EMBEDDER_ERROR_TOL: f64 = 0.00001;
const MAX_MINIMIZED_E_PER_ATOM: f64 = 0.05;
const MIN_TETRAHEDRAL_CHIRAL_VOL: f64 = 0.50;
const TETRAHEDRAL_CENTERINVOLUME_TOL: f64 = 0.30;

thread_local! {
    static RDKIT_DISTGEOM_RNG: std::cell::RefCell<RdkitDistgeomMinStdRand> =
        std::cell::RefCell::new(RdkitDistgeomMinStdRand::default());
    #[cfg(test)]
    static EMBEDDER_TEST_TRACE: std::cell::RefCell<Option<EmbedderTestTrace>> =
        const { std::cell::RefCell::new(None) };
}

fn have_opposite_sign(a: f64, b: f64) -> bool {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::haveOppositeSign (Embedder.cpp:67-69)
    // RDKit✔️✔️: inline bool haveOppositeSign(double a, double b) {
    // RDKit✔️✔️:   return std::signbit(a) ^ std::signbit(b);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::haveOppositeSign
    a.is_sign_negative() ^ b.is_sign_negative()
}

fn row64_distgeom_trace_enabled(num_points: usize) -> bool {
    env::var("RDKIT_ROW64_TRACE").ok().as_deref() == Some("1") && num_points == 106
}

fn row61_distgeom_trace_enabled(num_points: usize) -> bool {
    env::var("RDKIT_ROW61_TRACE").ok().as_deref() == Some("1") && num_points == 26
}

fn distgeom_trace_num_points_enabled(num_points: usize) -> bool {
    env::var("RDKIT_DISTGEOM_TRACE_NUM_POINTS")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        == Some(num_points)
}

fn distgeom_trace_bounds_pairs() -> Vec<(usize, usize)> {
    env::var("COSMOLKIT_DG_TRACE_BOUNDS_PAIRS")
        .or_else(|_| env::var("RDKIT_DG_TRACE_BOUNDS_PAIRS"))
        .ok()
        .map(|value| {
            value
                .split(';')
                .filter_map(|pair| {
                    let (a, b) = pair.split_once(',')?;
                    let a = a.trim().parse::<usize>().ok()?;
                    let b = b.trim().parse::<usize>().ok()?;
                    Some((a, b))
                })
                .collect()
        })
        .unwrap_or_default()
}

fn trace_bounds_stage(stage: &str, mmat: &BoundsMatrix) {
    for (i, j) in distgeom_trace_bounds_pairs() {
        if i < mmat.n && j < mmat.n {
            println!(
                "rust_bounds_stage stage={} pair={},{} lb={:.17} ub={:.17}",
                stage,
                i,
                j,
                mmat.get_lower(i, j),
                mmat.get_upper(i, j)
            );
        }
    }
}

fn row103_distgeom_trace_enabled(num_points: usize) -> bool {
    env::var("RDKIT_ROW103_TRACE").ok().as_deref() == Some("1") && num_points == 45
}

#[allow(dead_code)]
fn failmutex_get() -> &'static Mutex<()> {
    // BEGIN RDKIT CPP FUNCTION failmutex_get (Embedder.cpp:78-82)
    // RDKit✔️✔️: std::mutex &failmutex_get() {
    // RDKit✔️✔️:   // create on demand
    // RDKit✔️✔️:   static std::mutex _mutex;
    // RDKit✔️✔️:   return _mutex;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION failmutex_get
    static MUTEX: OnceLock<Mutex<()>> = OnceLock::new();
    MUTEX.get_or_init(|| Mutex::new(()))
}

#[allow(dead_code)]
fn failmutex_create() {
    // BEGIN RDKIT CPP FUNCTION failmutex_create (Embedder.cpp:84-87)
    // RDKit✔️✔️: void failmutex_create() {
    // RDKit✔️✔️:   std::mutex &mutex = failmutex_get();
    // RDKit✔️✔️:   std::lock_guard<std::mutex> test_lock(mutex);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION failmutex_create
    let _test_lock = failmutex_get()
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
}

#[allow(dead_code)]
fn get_fail_mutex() -> &'static Mutex<()> {
    // BEGIN RDKIT CPP FUNCTION GetFailMutex (Embedder.cpp:89-93)
    // RDKit✔️✔️: std::mutex &GetFailMutex() {
    // RDKit✔️✔️:   static std::once_flag flag;
    // RDKit✔️✔️:   std::call_once(flag, failmutex_create);
    // RDKit✔️✔️:   return failmutex_get();
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION GetFailMutex
    static INIT: OnceLock<()> = OnceLock::new();
    INIT.get_or_init(failmutex_create);
    failmutex_get()
}

fn normalized_point_delta(p0: ForceFieldVec3, p1: ForceFieldVec3) -> ForceFieldVec3 {
    let delta = p0 - p1;
    let length = delta.length();
    delta / length
}

fn embedder_volume_test(chiral_set: &ChiralSet, positions: &[ForceFieldVec3]) -> bool {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::_volumeTest (Embedder.cpp:316-397)
    // RDKit✔️✔️: bool _volumeTest(const DistGeom::ChiralSetPtr &chiralSet,
    // RDKit✔️✔️:                  const RDGeom::PointPtrVect &positions, bool verbose = false) {
    // RDKit✔️✔️:   RDGeom::Point3D p0((*positions[chiralSet->d_idx0])[0],
    // RDKit✔️✔️:                      (*positions[chiralSet->d_idx0])[1],
    // RDKit✔️✔️:                      (*positions[chiralSet->d_idx0])[2]);
    // RDKit✔️✔️:   RDGeom::Point3D p1((*positions[chiralSet->d_idx1])[0],
    // RDKit✔️✔️:                      (*positions[chiralSet->d_idx1])[1],
    // RDKit✔️✔️:                      (*positions[chiralSet->d_idx1])[2]);
    // RDKit✔️✔️:   RDGeom::Point3D p2((*positions[chiralSet->d_idx2])[0],
    // RDKit✔️✔️:                      (*positions[chiralSet->d_idx2])[1],
    // RDKit✔️✔️:                      (*positions[chiralSet->d_idx2])[2]);
    // RDKit✔️✔️:   RDGeom::Point3D p3((*positions[chiralSet->d_idx3])[0],
    // RDKit✔️✔️:                      (*positions[chiralSet->d_idx3])[1],
    // RDKit✔️✔️:                      (*positions[chiralSet->d_idx3])[2]);
    // RDKit✔️✔️:   RDGeom::Point3D p4((*positions[chiralSet->d_idx4])[0],
    // RDKit✔️✔️:                      (*positions[chiralSet->d_idx4])[1],
    // RDKit✔️✔️:                      (*positions[chiralSet->d_idx4])[2]);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // even if we are minimizing in higher dimension the chiral volume is
    // RDKit✔️✔️:   // calculated using only the first 3 dimensions
    // RDKit✔️✔️:   RDGeom::Point3D v1 = p0 - p1;
    // RDKit✔️✔️:   v1.normalize();
    // RDKit✔️✔️:   RDGeom::Point3D v2 = p0 - p2;
    // RDKit✔️✔️:   v2.normalize();
    // RDKit✔️✔️:   RDGeom::Point3D v3 = p0 - p3;
    // RDKit✔️✔️:   v3.normalize();
    // RDKit✔️✔️:   RDGeom::Point3D v4 = p0 - p4;
    // RDKit✔️✔️:   v4.normalize();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // be more tolerant of tethrahedral centers that are involved in multiple
    // RDKit✔️✔️:   // small rings
    // RDKit✔️✔️:   double volScale = 1;
    // RDKit✔️✔️:   if (chiralSet->d_structureFlags &
    // RDKit✔️✔️:       static_cast<std::uint64_t>(
    // RDKit✔️✔️:           DistGeom::ChiralSetStructureFlags::IN_FUSED_SMALL_RINGS)) {
    // RDKit✔️✔️:     volScale = 0.25;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   RDGeom::Point3D crossp = v1.crossProduct(v2);
    // RDKit✔️✔️:   double vol = crossp.dotProduct(v3);
    // RDKit✔️✔️:   if (fabs(vol) < volScale * MIN_TETRAHEDRAL_CHIRAL_VOL) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   crossp = v1.crossProduct(v2);
    // RDKit✔️✔️:   vol = crossp.dotProduct(v4);
    // RDKit✔️✔️:   if (fabs(vol) < volScale * MIN_TETRAHEDRAL_CHIRAL_VOL) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   crossp = v1.crossProduct(v3);
    // RDKit✔️✔️:   vol = crossp.dotProduct(v4);
    // RDKit✔️✔️:   if (fabs(vol) < volScale * MIN_TETRAHEDRAL_CHIRAL_VOL) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   crossp = v2.crossProduct(v3);
    // RDKit✔️✔️:   vol = crossp.dotProduct(v4);
    // RDKit✔️✔️:   return fabs(vol) >= volScale * MIN_TETRAHEDRAL_CHIRAL_VOL;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::_volumeTest
    let p0 = positions[chiral_set.idx0];
    let p1 = positions[chiral_set.idx1];
    let p2 = positions[chiral_set.idx2];
    let p3 = positions[chiral_set.idx3];
    let p4 = positions[chiral_set.idx4];

    let v1 = normalized_point_delta(p0, p1);
    let v2 = normalized_point_delta(p0, p2);
    let v3 = normalized_point_delta(p0, p3);
    let v4 = normalized_point_delta(p0, p4);

    let mut vol_scale = 1.0;
    if chiral_set.structure_flags & ChiralSetStructureFlags::InFusedSmallRings as u64 != 0 {
        vol_scale = 0.25;
    }
    let min_vol = vol_scale * MIN_TETRAHEDRAL_CHIRAL_VOL;

    let mut crossp = v1.cross_product(v2);
    let mut vol = crossp.dot_product(v3);
    if vol.abs() < min_vol {
        return false;
    }
    crossp = v1.cross_product(v2);
    vol = crossp.dot_product(v4);
    if vol.abs() < min_vol {
        return false;
    }
    crossp = v1.cross_product(v3);
    vol = crossp.dot_product(v4);
    if vol.abs() < min_vol {
        return false;
    }
    crossp = v2.cross_product(v3);
    vol = crossp.dot_product(v4);
    vol.abs() >= min_vol
}

fn embedder_same_side(
    v1: ForceFieldVec3,
    v2: ForceFieldVec3,
    v3: ForceFieldVec3,
    v4: ForceFieldVec3,
    p0: ForceFieldVec3,
    tol: f64,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::_sameSide (Embedder.cpp:399-410)
    // RDKit✔️✔️: bool _sameSide(const RDGeom::Point3D &v1, const RDGeom::Point3D &v2,
    // RDKit✔️✔️:                const RDGeom::Point3D &v3, const RDGeom::Point3D &v4,
    // RDKit✔️✔️:                const RDGeom::Point3D &p0, double tol = 0.1) {
    // RDKit✔️✔️:   RDGeom::Point3D normal = (v2 - v1).crossProduct(v3 - v1);
    // RDKit✔️✔️:   double d1 = normal.dotProduct(v4 - v1);
    // RDKit✔️✔️:   double d2 = normal.dotProduct(p0 - v1);
    // RDKit✔️✔️:   // std::cerr << "     " << d1 << " - " << d2 << std::endl;
    // RDKit✔️✔️:   if (fabs(d1) < tol || fabs(d2) < tol) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return !((d1 < 0.) ^ (d2 < 0.));
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::_sameSide
    let normal = (v2 - v1).cross_product(v3 - v1);
    let d1 = normal.dot_product(v4 - v1);
    let d2 = normal.dot_product(p0 - v1);
    if d1.abs() < tol || d2.abs() < tol {
        return false;
    }
    !((d1 < 0.0) ^ (d2 < 0.0))
}

fn embedder_center_in_volume_indices(
    idx0: usize,
    idx1: usize,
    idx2: usize,
    idx3: usize,
    idx4: usize,
    positions: &[ForceFieldVec3],
    tol: f64,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::_centerInVolume index overload (Embedder.cpp:411-437)
    // RDKit✔️✔️: bool _centerInVolume(unsigned int idx0, unsigned int idx1, unsigned int idx2,
    // RDKit✔️✔️:                      unsigned int idx3, unsigned int idx4,
    // RDKit✔️✔️:                      const RDGeom::PointPtrVect &positions, double tol,
    // RDKit✔️✔️:                      bool verbose = false) {
    // RDKit✔️✔️:   RDGeom::Point3D p0((*positions[idx0])[0], (*positions[idx0])[1],
    // RDKit✔️✔️:                      (*positions[idx0])[2]);
    // RDKit✔️✔️:   RDGeom::Point3D p1((*positions[idx1])[0], (*positions[idx1])[1],
    // RDKit✔️✔️:                      (*positions[idx1])[2]);
    // RDKit✔️✔️:   RDGeom::Point3D p2((*positions[idx2])[0], (*positions[idx2])[1],
    // RDKit✔️✔️:                      (*positions[idx2])[2]);
    // RDKit✔️✔️:   RDGeom::Point3D p3((*positions[idx3])[0], (*positions[idx3])[1],
    // RDKit✔️✔️:                      (*positions[idx3])[2]);
    // RDKit✔️✔️:   RDGeom::Point3D p4((*positions[idx4])[0], (*positions[idx4])[1],
    // RDKit✔️✔️:                      (*positions[idx4])[2]);
    // RDKit✔️✔️:   // RDGeom::Point3D centroid = (p1+p2+p3+p4)/4.;
    // RDKit✔️✔️:   bool res = _sameSide(p1, p2, p3, p4, p0, tol) &&
    // RDKit✔️✔️:              _sameSide(p2, p3, p4, p1, p0, tol) &&
    // RDKit✔️✔️:              _sameSide(p3, p4, p1, p2, p0, tol) &&
    // RDKit✔️✔️:              _sameSide(p4, p1, p2, p3, p0, tol);
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::_centerInVolume index overload
    let p0 = positions[idx0];
    let p1 = positions[idx1];
    let p2 = positions[idx2];
    let p3 = positions[idx3];
    let p4 = positions[idx4];
    embedder_same_side(p1, p2, p3, p4, p0, tol)
        && embedder_same_side(p2, p3, p4, p1, p0, tol)
        && embedder_same_side(p3, p4, p1, p2, p0, tol)
        && embedder_same_side(p4, p1, p2, p3, p0, tol)
}

fn embedder_center_in_volume(
    chiral_set: &ChiralSet,
    positions: &[ForceFieldVec3],
    tol: f64,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::_centerInVolume chiral-set overload (Embedder.cpp:439-450)
    // RDKit✔️✔️: bool _centerInVolume(const DistGeom::ChiralSetPtr &chiralSet,
    // RDKit✔️✔️:                      const RDGeom::PointPtrVect &positions, double tol = 0.1,
    // RDKit✔️✔️:                      bool verbose = false) {
    // RDKit✔️✔️:   if (chiralSet->d_idx0 ==
    // RDKit✔️✔️:       chiralSet->d_idx4) {  // this happens for three-coordinate centers
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return _centerInVolume(chiralSet->d_idx0, chiralSet->d_idx1,
    // RDKit✔️✔️:                          chiralSet->d_idx2, chiralSet->d_idx3,
    // RDKit✔️✔️:                          chiralSet->d_idx4, positions, tol, verbose);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::_centerInVolume chiral-set overload
    if chiral_set.idx0 == chiral_set.idx4 {
        return true;
    }
    embedder_center_in_volume_indices(
        chiral_set.idx0,
        chiral_set.idx1,
        chiral_set.idx2,
        chiral_set.idx3,
        chiral_set.idx4,
        positions,
        tol,
    )
}

fn embedder_bounds_fulfilled(
    atoms: &[i32],
    mmat: &BoundsMatrix,
    positions: &[ForceFieldVec3],
) -> bool {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::_boundsFulfilled (Embedder.cpp:451-476)
    // RDKit✔️✔️: bool _boundsFulfilled(const std::vector<int> &atoms,
    // RDKit✔️✔️:                       const DistGeom::BoundsMatrix &mmat,
    // RDKit✔️✔️:                       const RDGeom::PointPtrVect &positions) {
    // RDKit✔️✔️:   // unsigned int N = mmat.numRows();
    // RDKit✔️✔️:   // std::cerr << N << " " << atoms.size() << std::endl;
    // RDKit✔️✔️:   // loop over all pair of atoms
    // RDKit✔️✔️:   for (unsigned int i = 0; i < atoms.size() - 1; ++i) {
    // RDKit✔️✔️:     for (unsigned int j = i + 1; j < atoms.size(); ++j) {
    // RDKit✔️✔️:       int a1 = atoms[i];
    // RDKit✔️✔️:       int a2 = atoms[j];
    // RDKit✔️✔️:       RDGeom::Point3D p0((*positions[a1])[0], (*positions[a1])[1],
    // RDKit✔️✔️:                          (*positions[a1])[2]);
    // RDKit✔️✔️:       RDGeom::Point3D p1((*positions[a2])[0], (*positions[a2])[1],
    // RDKit✔️✔️:                          (*positions[a2])[2]);
    // RDKit✔️✔️:       double d2 = (p0 - p1).length();  // distance
    // RDKit✔️✔️:       double lb = mmat.getLowerBound(a1, a2);
    // RDKit✔️✔️:       double ub = mmat.getUpperBound(a1, a2);  // bounds
    // RDKit✔️✔️:       if (((d2 < lb) && (fabs(d2 - lb) > 0.1 * ub)) ||
    // RDKit✔️✔️:           ((d2 > ub) && (fabs(d2 - ub) > 0.1 * ub))) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::_boundsFulfilled
    if atoms.len() < 2 {
        return true;
    }
    for i in 0..atoms.len() - 1 {
        for j in i + 1..atoms.len() {
            let a1 = atoms[i] as usize;
            let a2 = atoms[j] as usize;
            let d2 = (positions[a1] - positions[a2]).length();
            let lb = mmat.get_lower(a1, a2);
            let ub = mmat.get_upper(a1, a2);
            if (d2 < lb && (d2 - lb).abs() > 0.1 * ub) || (d2 > ub && (d2 - ub).abs() > 0.1 * ub) {
                return false;
            }
        }
    }
    true
}

struct EmbedArgs<'a> {
    mmat: &'a BoundsMatrix,
    chiral_centers: &'a [ChiralSetPtr],
    tetrahedral_carbons: &'a [ChiralSetPtr],
    etkdg_details: Option<&'a CrystalFFDetails>,
    double_bond_ends: Option<&'a [(usize, usize, usize)]>,
    stereo_double_bonds: &'a [(Vec<usize>, i32)],
}

struct EmbedHelperArgs<'a> {
    confs_ok: &'a mut [bool],
    four_d: bool,
    frag_mapping: Option<&'a [usize]>,
    confs: &'a mut [Conformer3D],
    frag_idx: usize,
    mmat: &'a BoundsMatrix,
    chiral_centers: &'a [ChiralSetPtr],
    tetrahedral_carbons: &'a [ChiralSetPtr],
    double_bond_ends: &'a [(usize, usize, usize)],
    stereo_double_bonds: &'a [(Vec<usize>, i32)],
    etkdg_details: &'a CrystalFFDetails,
}

impl<'a> EmbedHelperArgs<'a> {
    fn embed_args(&self) -> EmbedArgs<'_> {
        EmbedArgs {
            mmat: self.mmat,
            chiral_centers: self.chiral_centers,
            tetrahedral_carbons: self.tetrahedral_carbons,
            etkdg_details: Some(self.etkdg_details),
            double_bond_ends: Some(self.double_bond_ends),
            stereo_double_bonds: self.stereo_double_bonds,
        }
    }
}

fn embedder_generate_initial_coords<R: RdkitDoubleRng>(
    positions: &mut [Vec<f64>],
    eargs: &EmbedArgs<'_>,
    embed_params: &EmbedParameters,
    dist_mat: &mut SymmMatrix,
    rng: &mut R,
) -> Result<bool, DgBoundsError> {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::generateInitialCoords (Embedder.cpp:479-520)
    // RDKit✔️✔️: bool generateInitialCoords(RDGeom::PointPtrVect *positions,
    // RDKit✔️✔️:                            const detail::EmbedArgs &eargs,
    // RDKit✔️✔️:                            const EmbedParameters &embedParams,
    // RDKit✔️✔️:                            RDNumeric::DoubleSymmMatrix &distMat,
    // RDKit✔️✔️:                            RDKit::double_source_type *rng) {
    // RDKit✔️✔️:   bool gotCoords = false;
    // RDKit✔️✔️:   if (!embedParams.useRandomCoords) {
    // RDKit✔️✔️:     double largestDistance =
    // RDKit✔️✔️:         DistGeom::pickRandomDistMat(*eargs.mmat, distMat, *rng);
    // RDKit✔️✔️:     RDUNUSED_PARAM(largestDistance);
    // RDKit✔️✔️:     gotCoords = DistGeom::computeInitialCoords(distMat, *positions, *rng,
    // RDKit✔️✔️:                                                embedParams.randNegEig,
    // RDKit✔️✔️:                                                embedParams.numZeroFail);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     double boxSize;
    // RDKit✔️✔️:     if (embedParams.boxSizeMult > 0) {
    // RDKit✔️✔️:       boxSize = 5. * embedParams.boxSizeMult;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       boxSize = -1 * embedParams.boxSizeMult;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     gotCoords = DistGeom::computeRandomCoords(*positions, boxSize, *rng);
    // RDKit✔️✔️:     if (embedParams.useRandomCoords && embedParams.coordMap != nullptr) {
    // RDKit✔️✔️:       for (const auto &v : *embedParams.coordMap) {
    // RDKit✔️✔️:         auto p = positions->at(v.first);
    // RDKit✔️✔️:         for (unsigned int ci = 0; ci < v.second.dimension(); ++ci) {
    // RDKit✔️✔️:           (*p)[ci] = v.second[ci];
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         // zero out any higher dimensional components:
    // RDKit✔️✔️:         for (unsigned int ci = v.second.dimension(); ci < p->dimension();
    // RDKit✔️✔️:              ++ci) {
    // RDKit✔️✔️:           (*p)[ci] = 0.0;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return gotCoords;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::generateInitialCoords
    let got_coords = if !embed_params.use_random_coords {
        let _largest_distance = pick_random_dist_mat_with_rng(eargs.mmat, dist_mat, rng);
        compute_initial_coords_with_rng(
            dist_mat,
            positions,
            rng,
            embed_params.rand_neg_eig,
            embed_params.num_zero_fail as usize,
        )?
    } else {
        let box_size = if embed_params.box_size_mult > 0.0 {
            5.0 * embed_params.box_size_mult
        } else {
            -embed_params.box_size_mult
        };
        let got_coords = compute_random_coords_with_rng(positions, box_size, rng);
        if let Some(coord_map) = &embed_params.coord_map {
            for (idx, mapped_point) in coord_map {
                let point = &mut positions[*idx as usize];
                for ci in 0..3 {
                    point[ci] = mapped_point[ci];
                }
                for coord in point.iter_mut().skip(3) {
                    *coord = 0.0;
                }
            }
        }
        got_coords
    };
    Ok(got_coords)
}

#[cfg(test)]
#[derive(Clone, Copy, Debug)]
pub(crate) enum EmbedderTestStage {
    InitialCoords,
    FirstMinimized,
    FourthDimensionCleaned,
    ExpTorsionMinimized,
    FinalChecked,
}

#[cfg(test)]
#[derive(Clone, Debug, Default, serde::Serialize)]
pub(crate) struct EmbedderTestLowLevelTrace {
    pub(crate) random_dist_preview: Option<Vec<f64>>,
    pub(crate) random_dist_full: Option<Vec<f64>>,
    pub(crate) random_dist_sum: Option<f64>,
    pub(crate) random_dist_sum_sq: Option<f64>,
    pub(crate) random_dist_largest: Option<f64>,
    pub(crate) initial_sum_sq_d2: Option<f64>,
    pub(crate) initial_sq_d0i_preview: Option<Vec<f64>>,
    pub(crate) initial_eig_vals_before_sqrt: Option<Vec<f64>>,
    pub(crate) initial_eig_vals_after_sqrt: Option<Vec<f64>>,
    pub(crate) first_minimization_initial_energy: Option<f64>,
    pub(crate) first_minimization_final_energy: Option<f64>,
    pub(crate) first_minimization_passes: Option<u32>,
    pub(crate) first_minimization_max_contrib: Option<f64>,
}

#[cfg(test)]
#[derive(Clone, Debug, Default)]
pub(crate) struct EmbedderTestTrace {
    pub(crate) initial_coords: Option<Vec<[f64; 3]>>,
    pub(crate) first_minimized: Option<Vec<[f64; 3]>>,
    pub(crate) fourth_dimension_cleaned: Option<Vec<[f64; 3]>>,
    pub(crate) exp_torsion_minimized: Option<Vec<[f64; 3]>>,
    pub(crate) final_checked: Option<Vec<[f64; 3]>>,
    pub(crate) failures: Vec<(String, u32)>,
    pub(crate) low_level: EmbedderTestLowLevelTrace,
}

#[cfg(test)]
fn embedder_test_low_level_trace_reset() {
    EMBEDDER_TEST_TRACE.with(|slot| {
        *slot.borrow_mut() = Some(EmbedderTestTrace::default());
    });
}

#[cfg(test)]
fn embedder_test_low_level_trace_take() -> EmbedderTestTrace {
    EMBEDDER_TEST_TRACE
        .with(|slot| slot.borrow_mut().take())
        .unwrap_or_default()
}

#[cfg(test)]
fn embedder_test_trace_update(update: impl FnOnce(&mut EmbedderTestTrace)) {
    EMBEDDER_TEST_TRACE.with(|slot| {
        if let Some(trace) = slot.borrow_mut().as_mut() {
            update(trace);
        }
    });
}

#[cfg(test)]
fn embedder_test_trace_set_coords(slot: &mut Option<Vec<[f64; 3]>>, positions: &[Vec<f64>]) {
    if slot.is_none() {
        *slot = Some(
            positions
                .iter()
                .map(|point| [point[0], point[1], point[2]])
                .collect(),
        );
    }
}

#[cfg(test)]
fn embedder_test_trace_failures(params: &EmbedParameters) -> Vec<(String, u32)> {
    params
        .failures
        .iter()
        .enumerate()
        .filter_map(|(idx, count)| {
            if *count == 0 {
                return None;
            }
            let cause = EmbedFailureCause::from_rdkit_ordinal(idx as u32)?;
            Some((cause.rdkit_name().to_string(), *count))
        })
        .collect()
}

#[cfg(test)]
pub(crate) fn embedder_stage_coords_for_test(
    mol: &Molecule,
    params: &mut EmbedParameters,
    stage: EmbedderTestStage,
) -> Result<Molecule, DgBoundsError> {
    let (mol_frags, frag_mapping) =
        molecule_fragments_for_embed(mol, params.embed_fragments_separately)?;
    let coord_map_storage = params.coord_map.clone();
    let mut coord_map = coord_map_storage.as_ref();
    if mol_frags.len() > 1 && coord_map.is_some() {
        coord_map = None;
    }

    let mut new_conformers = Vec::new();

    for (frag_idx, piece) in mol_frags.iter().enumerate() {
        let n_atoms = piece.num_atoms();
        let mut etkdg_details = CrystalFFDetails::default();
        embedder_init_etkdg(piece, params, &mut etkdg_details)?;

        let mut mmat;
        if params.bounds_mat.is_none() || mol_frags.len() > 1 {
            mmat = BoundsMatrix::new(n_atoms);
            if !embedder_setup_initial_bounds_matrix(
                piece,
                &mut mmat,
                coord_map,
                params,
                &mut etkdg_details,
            )? {
                continue;
            }
        } else {
            let bounds = params.bounds_mat.as_ref().expect("checked above");
            if bounds.num_rows() != n_atoms {
                return Err(DgBoundsError::GenerationFailed(
                    "size of boundsMat provided does not match the number of atoms in the molecule."
                        .to_string(),
                ));
            }
            collect_bonds_and_angles(piece, &mut etkdg_details.bonds, &mut etkdg_details.angles);
            mmat = (**bounds).clone();
        }

        let piece_storage;
        let piece = if piece
            .derived_cache()
            .rings
            .as_ref()
            .is_some_and(crate::RingInfo::is_initialized)
        {
            piece
        } else {
            let mut ringed_piece = piece.clone();
            let ring_info = ring_info_for_distgeom(piece)?.into_owned();
            ringed_piece.derived_cache_mut().rings = Some(ring_info);
            piece_storage = ringed_piece;
            &piece_storage
        };

        let mut chiral_centers = Vec::new();
        let mut tetrahedral_carbons = Vec::new();
        embedder_find_chiral_sets(
            piece,
            &mut chiral_centers,
            &mut tetrahedral_carbons,
            coord_map,
        );

        let mut double_bond_ends = Vec::new();
        let mut stereo_double_bonds = Vec::new();
        embedder_find_double_bonds(
            piece,
            &mut double_bond_ends,
            &mut stereo_double_bonds,
            coord_map,
        );

        let four_d = params.use_random_coords || !chiral_centers.is_empty();
        let dim = if four_d { 4 } else { 3 };
        let embed_args = EmbedArgs {
            mmat: &mmat,
            chiral_centers: &chiral_centers,
            tetrahedral_carbons: &tetrahedral_carbons,
            etkdg_details: Some(&etkdg_details),
            double_bond_ends: Some(&double_bond_ends),
            stereo_double_bonds: &stereo_double_bonds,
        };
        let mut positions = None;
        let max_iterations = if params.max_iterations == 0 {
            10 * n_atoms as u32
        } else {
            params.max_iterations
        };
        let mut rng = RdkitDistgeomMinStdRand::new_from_embed_points_seed(params.random_seed);
        for _iter in 0..max_iterations {
            let mut attempt_positions = vec![vec![0.0; dim]; n_atoms];
            let mut dist_mat = SymmMatrix::new(n_atoms);
            if !embedder_generate_initial_coords(
                &mut attempt_positions,
                &embed_args,
                params,
                &mut dist_mat,
                &mut rng,
            )? {
                continue;
            }
            if matches!(stage, EmbedderTestStage::InitialCoords) {
                positions = Some(attempt_positions);
                break;
            }

            if !embedder_first_minimization(&mut attempt_positions, &embed_args, params) {
                continue;
            }
            if matches!(stage, EmbedderTestStage::FirstMinimized) {
                positions = Some(attempt_positions);
                break;
            }

            if !embedder_check_tetrahedral_centers(&attempt_positions, &embed_args, params) {
                continue;
            }
            if params.enforce_chirality
                && !embed_args.chiral_centers.is_empty()
                && !embedder_check_chiral_centers(&attempt_positions, &embed_args, params)
            {
                continue;
            }
            if (!embed_args.chiral_centers.is_empty() || params.use_random_coords)
                && !embedder_minimize_fourth_dimension(
                    &mut attempt_positions,
                    &embed_args,
                    params,
                    None,
                )
            {
                continue;
            }
            if matches!(stage, EmbedderTestStage::FourthDimensionCleaned) {
                positions = Some(attempt_positions);
                break;
            }

            if (params.use_exp_torsion_angle_prefs || params.use_basic_knowledge)
                && !embedder_minimize_with_exp_torsions(&mut attempt_positions, &embed_args, params)
            {
                continue;
            }
            if matches!(stage, EmbedderTestStage::ExpTorsionMinimized) {
                positions = Some(attempt_positions);
                break;
            }

            if !embedder_double_bond_geometry_checks(
                &attempt_positions,
                &embed_args,
                params,
                1.0e-3,
            ) {
                continue;
            }
            if params.enforce_chirality {
                if !embed_args.chiral_centers.is_empty()
                    && !embedder_final_chiral_checks(&mut attempt_positions, &embed_args, params)
                {
                    continue;
                }
                if !embed_args.stereo_double_bonds.is_empty()
                    && !embedder_double_bond_stereo_checks(&attempt_positions, &embed_args, params)
                {
                    continue;
                }
            }
            positions = Some(attempt_positions);
            break;
        }
        let Some(positions) = positions else {
            continue;
        };

        let mut conf = Conformer3D::new(frag_idx, vec![[0.0; 3]; n_atoms], true);
        for i in 0..n_atoms {
            conf.coordinates_mut()[i] = [positions[i][0], positions[i][1], positions[i][2]];
        }
        new_conformers.push((frag_idx, conf, frag_mapping.clone()));
    }

    let mut full_conformers = Vec::with_capacity(new_conformers.len());
    for (frag_idx, conformer, mapping) in new_conformers {
        let mut full = Conformer3D::new(frag_idx, vec![[0.0; 3]; mol.num_atoms()], true);
        let mut frag_atom_idx = 0usize;
        for i in 0..mol.num_atoms() {
            if mapping.get(i).copied().unwrap_or_default() == frag_idx {
                full.coordinates_mut()[i] = conformer.coordinates()[frag_atom_idx];
                frag_atom_idx += 1;
            }
        }
        full_conformers.push(full);
    }
    embedder_add_conformers(mol, full_conformers, params.clear_confs).map(|(mol, _)| mol)
}

#[cfg(test)]
pub(crate) fn embedder_trace_for_test(
    mol: &Molecule,
    params: &mut EmbedParameters,
) -> Result<EmbedderTestTrace, DgBoundsError> {
    let (mol_frags, _frag_mapping) =
        molecule_fragments_for_embed(mol, params.embed_fragments_separately)?;
    let coord_map_storage = params.coord_map.clone();
    let mut coord_map = coord_map_storage.as_ref();
    if mol_frags.len() > 1 && coord_map.is_some() {
        coord_map = None;
    }
    let previous_track_failures = params.track_failures;
    params.track_failures = true;
    params
        .failures
        .resize(EmbedFailureCause::EndOfEnum as usize, 0);
    params.failures.fill(0);

    let mut trace = EmbedderTestTrace::default();
    for (frag_idx, piece) in mol_frags.iter().enumerate() {
        if frag_idx > 0 {
            break;
        }
        let n_atoms = piece.num_atoms();
        let mut etkdg_details = CrystalFFDetails::default();
        embedder_init_etkdg(piece, params, &mut etkdg_details)?;

        let mut mmat;
        if params.bounds_mat.is_none() || mol_frags.len() > 1 {
            mmat = BoundsMatrix::new(n_atoms);
            if !embedder_setup_initial_bounds_matrix(
                piece,
                &mut mmat,
                coord_map,
                params,
                &mut etkdg_details,
            )? {
                continue;
            }
        } else {
            let bounds = params.bounds_mat.as_ref().expect("checked above");
            if bounds.num_rows() != n_atoms {
                return Err(DgBoundsError::GenerationFailed(
                    "size of boundsMat provided does not match the number of atoms in the molecule."
                        .to_string(),
                ));
            }
            collect_bonds_and_angles(piece, &mut etkdg_details.bonds, &mut etkdg_details.angles);
            mmat = (**bounds).clone();
        }

        let piece_storage;
        let piece = if piece
            .derived_cache()
            .rings
            .as_ref()
            .is_some_and(crate::RingInfo::is_initialized)
        {
            piece
        } else {
            let mut ringed_piece = piece.clone();
            let ring_info = ring_info_for_distgeom(piece)?.into_owned();
            ringed_piece.derived_cache_mut().rings = Some(ring_info);
            piece_storage = ringed_piece;
            &piece_storage
        };

        let mut chiral_centers = Vec::new();
        let mut tetrahedral_carbons = Vec::new();
        embedder_find_chiral_sets(
            piece,
            &mut chiral_centers,
            &mut tetrahedral_carbons,
            coord_map,
        );

        let mut double_bond_ends = Vec::new();
        let mut stereo_double_bonds = Vec::new();
        embedder_find_double_bonds(
            piece,
            &mut double_bond_ends,
            &mut stereo_double_bonds,
            coord_map,
        );

        let four_d = params.use_random_coords || !chiral_centers.is_empty();
        let dim = if four_d { 4 } else { 3 };
        let embed_args = EmbedArgs {
            mmat: &mmat,
            chiral_centers: &chiral_centers,
            tetrahedral_carbons: &tetrahedral_carbons,
            etkdg_details: Some(&etkdg_details),
            double_bond_ends: Some(&double_bond_ends),
            stereo_double_bonds: &stereo_double_bonds,
        };
        let mut attempt_positions = vec![vec![0.0; dim]; n_atoms];
        embedder_test_low_level_trace_reset();
        let _ = embedder_embed_points(
            &mut attempt_positions,
            embed_args,
            params,
            params.random_seed,
            None,
        )?;
        trace = embedder_test_low_level_trace_take();
    }
    trace.failures = embedder_test_trace_failures(params);
    params.track_failures = previous_track_failures;
    Ok(trace)
}

fn copy_forcefield_positions_to_point_vectors(field: &ForceField, positions: &mut [Vec<f64>]) {
    for (point, field_point) in positions.iter_mut().zip(field.positions()) {
        for ci in 0..point.len() {
            point[ci] = field_point[ci];
        }
    }
}

fn point_vectors_to_forcefield_vec3(positions: &[Vec<f64>]) -> Vec<ForceFieldVec3> {
    positions
        .iter()
        .map(|point| ForceFieldVec3::new(point[0], point[1], point[2]))
        .collect()
}

fn embedder_first_minimization(
    positions: &mut [Vec<f64>],
    eargs: &EmbedArgs<'_>,
    embed_params: &EmbedParameters,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::firstMinimization (Embedder.cpp:521-561)
    // RDKit✔️✔️: bool firstMinimization(RDGeom::PointPtrVect *positions,
    // RDKit✔️✔️:                        const detail::EmbedArgs &eargs,
    // RDKit✔️✔️:                        const EmbedParameters &embedParams) {
    // RDKit✔️✔️:   bool gotCoords = true;
    // RDKit✔️✔️:   boost::dynamic_bitset<> fixedPts(positions->size());
    // RDKit✔️✔️:   if (embedParams.useRandomCoords && embedParams.coordMap != nullptr) {
    // RDKit✔️✔️:     for (const auto &v : *embedParams.coordMap) {
    // RDKit✔️✔️:       fixedPts.set(v.first);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::unique_ptr<ForceFields::ForceField> field(DistGeom::constructForceField(
    // RDKit✔️✔️:       *eargs.mmat, *positions, *eargs.chiralCenters, 1.0, 0.1, nullptr,
    // RDKit✔️✔️:       embedParams.basinThresh, &fixedPts));
    // RDKit✔️✔️:   if (embedParams.useRandomCoords && embedParams.coordMap != nullptr) {
    // RDKit✔️✔️:     for (const auto &v : *embedParams.coordMap) {
    // RDKit✔️✔️:       field->fixedPoints().push_back(v.first);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   field->initialize();
    // RDKit✔️✔️:   if (field->calcEnergy() > ERROR_TOL) {
    // RDKit✔️✔️:     int needMore = 1;
    // RDKit✔️✔️:     while (needMore) {
    // RDKit✔️✔️:       needMore = field->minimize(400, embedParams.optimizerForceTol);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::vector<double> e_contribs;
    // RDKit✔️✔️:   double local_e = field->calcEnergy(&e_contribs);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // check that the energy is not too high (this is part of github #971)
    // RDKit✔️✔️:   if (local_e / positions->size() >= MAX_MINIMIZED_E_PER_ATOM) {
    // RDKit✔️✔️:     gotCoords = false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return gotCoords;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::firstMinimization
    let mut got_coords = true;
    let mut fixed_pts = vec![false; positions.len()];
    if embed_params.use_random_coords
        && let Some(coord_map) = &embed_params.coord_map
    {
        for idx in coord_map.keys() {
            fixed_pts[*idx as usize] = true;
        }
    }
    let mut field = construct_distgeom_forcefield(
        eargs.mmat,
        positions,
        eargs.chiral_centers,
        1.0,
        0.1,
        None,
        embed_params.basin_thresh,
        Some(&fixed_pts),
    );
    if embed_params.use_random_coords
        && let Some(coord_map) = &embed_params.coord_map
    {
        for idx in coord_map.keys() {
            field.fixed_points_mut().push(*idx as usize);
        }
    }
    field.initialize();
    let initial_energy = field.calc_energy_current(None);
    #[cfg(test)]
    embedder_test_trace_update(|trace| {
        trace
            .low_level
            .first_minimization_initial_energy
            .get_or_insert(initial_energy);
    });
    let mut minimize_passes = 0_u32;
    if initial_energy > EMBEDDER_ERROR_TOL {
        let mut need_more = 1;
        while need_more != 0 {
            minimize_passes += 1;
            need_more = field.minimize(400, embed_params.optimizer_force_tol, 1.0e-6);
        }
    }
    copy_forcefield_positions_to_point_vectors(&field, positions);
    let mut e_contribs = Vec::new();
    let local_e = field.calc_energy_current(Some(&mut e_contribs));
    #[cfg(test)]
    embedder_test_trace_update(|trace| {
        trace
            .low_level
            .first_minimization_final_energy
            .get_or_insert(local_e);
        trace
            .low_level
            .first_minimization_passes
            .get_or_insert(minimize_passes);
        trace
            .low_level
            .first_minimization_max_contrib
            .get_or_insert_with(|| {
                e_contribs
                    .iter()
                    .copied()
                    .fold(0.0_f64, |acc, value| acc.max(value))
            });
    });
    if local_e / positions.len() as f64 >= MAX_MINIMIZED_E_PER_ATOM {
        got_coords = false;
    }
    got_coords
}

fn embedder_check_tetrahedral_centers(
    positions: &[Vec<f64>],
    eargs: &EmbedArgs<'_>,
    _embed_params: &EmbedParameters,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::checkTetrahedralCenters (Embedder.cpp:563-584)
    // RDKit✔️✔️: bool checkTetrahedralCenters(const RDGeom::PointPtrVect *positions,
    // RDKit✔️✔️:                              const detail::EmbedArgs &eargs,
    // RDKit✔️✔️:                              const EmbedParameters &) {
    // RDKit✔️✔️:   // for each of the atoms in the "tetrahedralCarbons" list, make sure
    // RDKit✔️✔️:   // that there is a minimum volume around them and that they are inside
    // RDKit✔️✔️:   // that volume. (this is part of github #971)
    // RDKit✔️✔️:   for (const auto &tetSet : *eargs.tetrahedralCarbons) {
    // RDKit✔️✔️:     // it could happen that the centroid is outside the volume defined
    // RDKit✔️✔️:     // by the other
    // RDKit✔️✔️:     // four points. That is also a fail.
    // RDKit✔️✔️:     if (!_volumeTest(tetSet, *positions) ||
    // RDKit✔️✔️:         !_centerInVolume(tetSet, *positions, TETRAHEDRAL_CENTERINVOLUME_TOL)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::checkTetrahedralCenters
    let positions = point_vectors_to_forcefield_vec3(positions);
    let trace_row103 = env::var("RDKIT_ROW103_TRACE").ok().as_deref() == Some("1")
        && positions.len() == 40
        && eargs.tetrahedral_carbons.len() == 3;
    for tet_set in eargs.tetrahedral_carbons {
        let volume_ok = embedder_volume_test(tet_set, &positions);
        let center_ok =
            embedder_center_in_volume(tet_set, &positions, TETRAHEDRAL_CENTERINVOLUME_TOL);
        if trace_row103 {
            println!(
                "row103_tetra_check idx0={} atoms=[{},{},{},{},{}] volume_ok={} center_ok={}",
                tet_set.idx0,
                tet_set.idx0,
                tet_set.idx1,
                tet_set.idx2,
                tet_set.idx3,
                tet_set.idx4,
                volume_ok,
                center_ok
            );
        }
        if !volume_ok || !center_ok {
            return false;
        }
    }
    true
}

fn embedder_check_chiral_centers(
    positions: &[Vec<f64>],
    eargs: &EmbedArgs<'_>,
    _embed_params: &EmbedParameters,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::checkChiralCenters (Embedder.cpp:585-605)
    // RDKit✔️✔️: bool checkChiralCenters(const RDGeom::PointPtrVect *positions,
    // RDKit✔️✔️:                         const detail::EmbedArgs &eargs,
    // RDKit✔️✔️:                         const EmbedParameters &) {
    // RDKit✔️✔️:   // check the chiral volume:
    // RDKit✔️✔️:   for (const auto &chiralSet : *eargs.chiralCenters) {
    // RDKit✔️✔️:     double vol = DistGeom::calcChiralVolume(
    // RDKit✔️✔️:         chiralSet->d_idx1, chiralSet->d_idx2, chiralSet->d_idx3,
    // RDKit✔️✔️:         chiralSet->d_idx4, *positions);
    // RDKit✔️✔️:     double lb = chiralSet->getLowerVolumeBound();
    // RDKit✔️✔️:     double ub = chiralSet->getUpperVolumeBound();
    // RDKit✔️✔️:     if ((lb > 0 && vol < lb && (vol / lb < .8 || haveOppositeSign(vol, lb))) ||
    // RDKit✔️✔️:         (ub < 0 && vol > ub && (vol / ub < .8 || haveOppositeSign(vol, ub)))) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::checkChiralCenters
    let positions = point_vectors_to_forcefield_vec3(positions);
    for chiral_set in eargs.chiral_centers {
        let vol = calc_chiral_volume_points(
            chiral_set.idx1,
            chiral_set.idx2,
            chiral_set.idx3,
            chiral_set.idx4,
            &positions,
        );
        let lb = chiral_set.get_lower_volume_bound();
        let ub = chiral_set.get_upper_volume_bound();
        if (lb > 0.0 && vol < lb && (vol / lb < 0.8 || have_opposite_sign(vol, lb)))
            || (ub < 0.0 && vol > ub && (vol / ub < 0.8 || have_opposite_sign(vol, ub)))
        {
            return false;
        }
    }
    true
}

fn embedder_minimize_fourth_dimension(
    positions: &mut [Vec<f64>],
    eargs: &EmbedArgs<'_>,
    embed_params: &mut EmbedParameters,
    end_time: Option<Instant>,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::minimizeFourthDimension (Embedder.cpp:606-637)
    // RDKit✔️✔️: bool minimizeFourthDimension(RDGeom::PointPtrVect *positions,
    // RDKit✔️✔️:                              const detail::EmbedArgs &eargs,
    // RDKit✔️✔️:                              EmbedParameters &embedParams,
    // RDKit✔️✔️:                              TimePoint *end_time) {
    // RDKit✔️✔️:   // now redo the minimization if we have a chiral center
    // RDKit✔️✔️:   // or have started from random coords. This
    // RDKit✔️✔️:   // time removing the chiral constraints and
    // RDKit✔️✔️:   // increasing the weight on the fourth dimension
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::unique_ptr<ForceFields::ForceField> field2(DistGeom::constructForceField(
    // RDKit✔️✔️:       *eargs.mmat, *positions, *eargs.chiralCenters, 0.2, 1.0, nullptr,
    // RDKit✔️✔️:       embedParams.basinThresh));
    // RDKit✔️✔️:   if (embedParams.useRandomCoords && embedParams.coordMap != nullptr) {
    // RDKit✔️✔️:     for (const auto &v : *embedParams.coordMap) {
    // RDKit✔️✔️:       field2->fixedPoints().push_back(v.first);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   field2->initialize();
    // RDKit✔️✔️:   // std::cerr << "FIELD2 E: " << field2->calcEnergy() << std::endl;
    // RDKit✔️✔️:   if (field2->calcEnergy() > ERROR_TOL) {
    // RDKit✔️✔️:     int needMore = 1;
    // RDKit✔️✔️:     while (needMore) {
    // RDKit✔️✔️:       if (end_time != nullptr && Clock::now() > *end_time) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       needMore = field2->minimize(200, embedParams.optimizerForceTol);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::minimizeFourthDimension
    let mut field = construct_distgeom_forcefield(
        eargs.mmat,
        positions,
        eargs.chiral_centers,
        0.2,
        1.0,
        None,
        embed_params.basin_thresh,
        None,
    );
    if embed_params.use_random_coords
        && let Some(coord_map) = &embed_params.coord_map
    {
        for idx in coord_map.keys() {
            field.fixed_points_mut().push(*idx as usize);
        }
    }

    field.initialize();
    let trace_row64 = row64_distgeom_trace_enabled(positions.len());
    let trace_num_points = distgeom_trace_num_points_enabled(positions.len());
    if trace_row64 || trace_num_points {
        println!(
            "distgeom_fourth_dimension n={} initial_energy={:.15}",
            positions.len(),
            field.calc_energy_current(None)
        );
    }
    if field.calc_energy_current(None) > EMBEDDER_ERROR_TOL {
        let mut need_more = 1;
        let mut pass = 0usize;
        while need_more != 0 {
            if let Some(deadline) = end_time
                && Instant::now() > deadline
            {
                if trace_row64 || trace_num_points {
                    println!(
                        "distgeom_fourth_dimension n={} timeout pass={pass}",
                        positions.len()
                    );
                }
                copy_forcefield_positions_to_point_vectors(&field, positions);
                return false;
            }
            need_more = field.minimize(200, embed_params.optimizer_force_tol, 1.0e-6);
            if trace_row64 || trace_num_points {
                println!(
                    "distgeom_fourth_dimension n={} pass={} need_more={} energy={:.15}",
                    positions.len(),
                    pass,
                    need_more,
                    field.calc_energy_current(None)
                );
            }
            pass += 1;
        }
    }
    copy_forcefield_positions_to_point_vectors(&field, positions);
    if trace_row64 {
        print_debug_like_row64_positions("row64_after_fourth_dimension_internal", positions);
    }
    true
}

fn print_debug_like_row64_positions(stage: &str, positions: &[Vec<f64>]) {
    let coords: Vec<String> = positions
        .iter()
        .map(|point| {
            format!(
                "[{:.15},{:.15},{:.15},{:.15}]",
                point[0],
                point[1],
                point[2],
                *point.get(3).unwrap_or(&0.0)
            )
        })
        .collect();
    println!("{stage} coords={}", coords.join(";"));
}

fn embedder_minimize_with_exp_torsions(
    positions: &mut [Vec<f64>],
    eargs: &EmbedArgs<'_>,
    embed_params: &EmbedParameters,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::minimizeWithExpTorsions (Embedder.cpp:638-701)
    // RDKit✔️✔️: bool minimizeWithExpTorsions(RDGeom::PointPtrVect &positions,
    // RDKit✔️✔️:                              const detail::EmbedArgs &eargs,
    // RDKit✔️✔️:                              const EmbedParameters &embedParams) {
    // RDKit✔️✔️:   PRECONDITION(eargs.etkdgDetails, "bogus etkdgDetails pointer");
    // RDKit✔️✔️:   bool planar = true;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // convert to 3D positions and create coordMap
    // RDKit✔️✔️:   RDGeom::Point3DPtrVect positions3D;
    // RDKit✔️✔️:   for (auto &position : positions) {
    // RDKit✔️✔️:     positions3D.push_back(
    // RDKit✔️✔️:         new RDGeom::Point3D((*position)[0], (*position)[1], (*position)[2]));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // create the force field
    // RDKit✔️✔️:   std::unique_ptr<ForceFields::ForceField> field;
    // RDKit✔️✔️:   if (embedParams.useBasicKnowledge) {  // ETKDG or KDG
    // RDKit✔️✔️:     if (embedParams.CPCI != nullptr) {
    // RDKit✔️✔️:       field.reset(DistGeom::construct3DForceField(
    // RDKit✔️✔️:           *eargs.mmat, positions3D, *eargs.etkdgDetails, *embedParams.CPCI));
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       field.reset(DistGeom::construct3DForceField(*eargs.mmat, positions3D,
    // RDKit✔️✔️:                                                   *eargs.etkdgDetails));
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {  // plain ETDG
    // RDKit✔️✔️:     field.reset(DistGeom::constructPlain3DForceField(*eargs.mmat, positions3D,
    // RDKit✔️✔️:                                                      *eargs.etkdgDetails));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (embedParams.useRandomCoords && embedParams.coordMap != nullptr) {
    // RDKit✔️✔️:     for (const auto &v : *embedParams.coordMap) {
    // RDKit✔️✔️:       field->fixedPoints().push_back(v.first);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // minimize!
    // RDKit✔️✔️:   field->initialize();
    // RDKit✔️✔️:   if (field->calcEnergy() > ERROR_TOL) {
    // RDKit✔️✔️:     // while (needMore) {
    // RDKit✔️✔️:     field->minimize(300, embedParams.optimizerForceTol);
    // RDKit✔️✔️:     //      ++nPasses;
    // RDKit✔️✔️:     //}
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // std::cout << field->calcEnergy() << std::endl;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // check for planarity if ETKDG or KDG
    // RDKit✔️✔️:   if (embedParams.useBasicKnowledge) {
    // RDKit✔️✔️:     // create a force field with only the impropers
    // RDKit✔️✔️:     std::unique_ptr<ForceFields::ForceField> field2(
    // RDKit✔️✔️:         DistGeom::construct3DImproperForceField(*eargs.mmat, positions3D,
    // RDKit✔️✔️:                                                 *eargs.etkdgDetails));
    // RDKit✔️✔️:     if (embedParams.useRandomCoords && embedParams.coordMap != nullptr) {
    // RDKit✔️✔️:       for (const auto &v : *embedParams.coordMap) {
    // RDKit✔️✔️:         field2->fixedPoints().push_back(v.first);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     field2->initialize();
    // RDKit✔️✔️:     // check if the energy is low enough
    // RDKit✔️✔️:     double planarityTolerance = 0.7;
    // RDKit✔️✔️:     if (field2->calcEnergy() >
    // RDKit✔️✔️:         eargs.etkdgDetails->improperAtoms.size() * planarityTolerance) {
    // RDKit✔️✔️:       planar = false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // overwrite positions and delete the 3D ones
    // RDKit✔️✔️:   for (unsigned int i = 0; i < positions3D.size(); ++i) {
    // RDKit✔️✔️:     (*positions[i])[0] = (*positions3D[i])[0];
    // RDKit✔️✔️:     (*positions[i])[1] = (*positions3D[i])[1];
    // RDKit✔️✔️:     (*positions[i])[2] = (*positions3D[i])[2];
    // RDKit✔️✔️:     delete positions3D[i];
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return planar;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::minimizeWithExpTorsions
    let etkdg_details = eargs.etkdg_details.expect("bogus etkdgDetails pointer");
    let mut planar = true;
    let mut positions_3d = point_vectors_to_forcefield_vec3(positions);

    let mut field = if embed_params.use_basic_knowledge {
        if let Some(cpci) = &embed_params.cpci {
            let cpci_usize: BTreeMap<(usize, usize), f64> = cpci
                .iter()
                .map(|(&(idx1, idx2), &charge)| ((idx1 as usize, idx2 as usize), charge))
                .collect();
            construct_3d_forcefield_with_cpci(eargs.mmat, &positions_3d, etkdg_details, &cpci_usize)
        } else {
            construct_3d_forcefield(eargs.mmat, &positions_3d, etkdg_details)
        }
    } else {
        construct_plain_3d_forcefield(eargs.mmat, &positions_3d, etkdg_details)
    };
    if embed_params.use_random_coords
        && let Some(coord_map) = &embed_params.coord_map
    {
        for idx in coord_map.keys() {
            field.fixed_points_mut().push(*idx as usize);
        }
    }

    field.initialize();
    if field.calc_energy_current(None) > EMBEDDER_ERROR_TOL {
        field.minimize(300, embed_params.optimizer_force_tol, 1.0e-6);
    }
    positions_3d.clone_from_slice(field.positions());

    if embed_params.use_basic_knowledge {
        let mut field2 = construct_3d_improper_forcefield(eargs.mmat, &positions_3d, etkdg_details);
        if embed_params.use_random_coords
            && let Some(coord_map) = &embed_params.coord_map
        {
            for idx in coord_map.keys() {
                field2.fixed_points_mut().push(*idx as usize);
            }
        }

        field2.initialize();
        let planarity_tolerance = 0.7;
        if field2.calc_energy_current(None)
            > etkdg_details.improper_atoms.len() as f64 * planarity_tolerance
        {
            planar = false;
        }
    }

    for (point, point3d) in positions.iter_mut().zip(positions_3d.iter()) {
        point[0] = point3d[0];
        point[1] = point3d[1];
        point[2] = point3d[2];
    }
    planar
}

fn embedder_double_bond_geometry_checks(
    positions: &[Vec<f64>],
    eargs: &EmbedArgs<'_>,
    _embed_params: &mut EmbedParameters,
    linear_tol: f64,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::doubleBondGeometryChecks (Embedder.cpp:703-731)
    // RDKit✔️✔️: bool doubleBondGeometryChecks(const RDGeom::PointPtrVect &positions,
    // RDKit✔️✔️:                               const detail::EmbedArgs &eargs, EmbedParameters &,
    // RDKit✔️✔️:                               double linearTol = 1e-3) {
    // RDKit✔️✔️:   if (eargs.doubleBondEnds) {
    // RDKit✔️✔️:     for (const auto &itm : *eargs.doubleBondEnds) {
    // RDKit✔️✔️:       const auto &a0 = *positions[std::get<0>(itm)];
    // RDKit✔️✔️:       const auto &a1 = *positions[std::get<1>(itm)];
    // RDKit✔️✔️:       const auto &a2 = *positions[std::get<2>(itm)];
    // RDKit✔️✔️:       RDGeom::Point3D p0(a0[0], a0[1], a0[2]);
    // RDKit✔️✔️:       RDGeom::Point3D p1(a1[0], a1[1], a1[2]);
    // RDKit✔️✔️:       RDGeom::Point3D p2(a2[0], a2[1], a2[2]);
    // RDKit✔️✔️:
    // RDKit✔️✔️:       // check for a linear arrangement
    // RDKit✔️✔️:
    // RDKit✔️✔️:       auto v1 = p1 - p0;
    // RDKit✔️✔️:       v1.normalize();
    // RDKit✔️✔️:       auto v2 = p1 - p2;
    // RDKit✔️✔️:       v2.normalize();
    // RDKit✔️✔️:       if (v1.dotProduct(v2) + 1.0 < linearTol) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::doubleBondGeometryChecks
    if let Some(double_bond_ends) = eargs.double_bond_ends {
        for &(idx0, idx1, idx2) in double_bond_ends {
            let p0 =
                ForceFieldVec3::new(positions[idx0][0], positions[idx0][1], positions[idx0][2]);
            let p1 =
                ForceFieldVec3::new(positions[idx1][0], positions[idx1][1], positions[idx1][2]);
            let p2 =
                ForceFieldVec3::new(positions[idx2][0], positions[idx2][1], positions[idx2][2]);
            let mut v1 = p1 - p0;
            v1 /= v1.length();
            let mut v2 = p1 - p2;
            v2 /= v2.length();
            if v1.dot_product(v2) + 1.0 < linear_tol {
                return false;
            }
        }
    }
    true
}

fn embedder_double_bond_stereo_checks(
    positions: &[Vec<f64>],
    eargs: &EmbedArgs<'_>,
    _embed_params: &mut EmbedParameters,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::doubleBondStereoChecks (Embedder.cpp:733-768)
    // RDKit✔️✔️: bool doubleBondStereoChecks(const RDGeom::PointPtrVect &positions,
    // RDKit✔️✔️:                             const detail::EmbedArgs &eargs, EmbedParameters &) {
    // RDKit✔️✔️:   for (const auto &itm : *eargs.stereoDoubleBonds) {
    // RDKit✔️✔️:     // itm is a pair with [controlling_atoms], sign
    // RDKit✔️✔️:     // where the sign tells us about cis/trans
    // RDKit✔️✔️:
    // RDKit✔️✔️:     const auto &a0 = *positions[itm.first[0]];
    // RDKit✔️✔️:     const auto &a1 = *positions[itm.first[1]];
    // RDKit✔️✔️:     const auto &a2 = *positions[itm.first[2]];
    // RDKit✔️✔️:     const auto &a3 = *positions[itm.first[3]];
    // RDKit✔️✔️:     RDGeom::Point3D p0(a0[0], a0[1], a0[2]);
    // RDKit✔️✔️:     RDGeom::Point3D p1(a1[0], a1[1], a1[2]);
    // RDKit✔️✔️:     RDGeom::Point3D p2(a2[0], a2[1], a2[2]);
    // RDKit✔️✔️:     RDGeom::Point3D p3(a3[0], a3[1], a3[2]);
    // RDKit✔️✔️:     auto dihedral = RDGeom::computeDihedralAngle(p0, p1, p2, p3);
    // RDKit✔️✔️:     if ((dihedral - M_PI_2) * itm.second < 0) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::doubleBondStereoChecks
    const M_PI_2: f64 = std::f64::consts::FRAC_PI_2;
    for (atoms, sign) in eargs.stereo_double_bonds {
        let p0 = ForceFieldVec3::new(
            positions[atoms[0]][0],
            positions[atoms[0]][1],
            positions[atoms[0]][2],
        );
        let p1 = ForceFieldVec3::new(
            positions[atoms[1]][0],
            positions[atoms[1]][1],
            positions[atoms[1]][2],
        );
        let p2 = ForceFieldVec3::new(
            positions[atoms[2]][0],
            positions[atoms[2]][1],
            positions[atoms[2]][2],
        );
        let p3 = ForceFieldVec3::new(
            positions[atoms[3]][0],
            positions[atoms[3]][1],
            positions[atoms[3]][2],
        );
        let beg_end_vec = p2 - p1;
        let beg_nbr_vec = p0 - p1;
        let crs1 = beg_nbr_vec.cross_product(beg_end_vec);
        let end_nbr_vec = p3 - p2;
        let crs2 = end_nbr_vec.cross_product(beg_end_vec);
        let dihedral = (crs1.dot_product(crs2) / (crs1.length() * crs2.length()))
            .clamp(-1.0, 1.0)
            .acos();
        if (dihedral - M_PI_2) * f64::from(*sign) < 0.0 {
            return false;
        }
    }
    true
}

fn embedder_increment_failure(embed_params: &mut EmbedParameters, cause: EmbedFailureCause) {
    if embed_params.track_failures {
        let _guard = failmutex_get().lock().expect("failure mutex poisoned");
        embed_params.failures[cause as usize] += 1;
    }
}

fn embedder_final_chiral_checks(
    positions: &mut [Vec<f64>],
    eargs: &EmbedArgs<'_>,
    embed_params: &mut EmbedParameters,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::finalChiralChecks (Embedder.cpp:770-836)
    // RDKit✔️✔️: bool finalChiralChecks(RDGeom::PointPtrVect *positions,
    // RDKit✔️✔️:                        const detail::EmbedArgs &eargs,
    // RDKit✔️✔️:                        EmbedParameters &embedParams) {
    // RDKit✔️✔️:   // confirm chiral volumes
    // RDKit✔️✔️:   if (!checkChiralCenters(positions, eargs, embedParams)) {
    // RDKit✔️✔️:     if (embedParams.trackFailures) {
    // RDKit✔️✔️:       embedParams.failures[EmbedFailureCauses::CHECK_CHIRAL_CENTERS2]++;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // "distance matrix" chirality test
    // RDKit✔️✔️:   std::set<int> atoms;
    // RDKit✔️✔️:   for (const auto &chiralSet : *eargs.chiralCenters) {
    // RDKit✔️✔️:     if (chiralSet->d_idx0 != chiralSet->d_idx4) {
    // RDKit✔️✔️:       atoms.insert(chiralSet->d_idx0);
    // RDKit✔️✔️:       atoms.insert(chiralSet->d_idx1);
    // RDKit✔️✔️:       atoms.insert(chiralSet->d_idx2);
    // RDKit✔️✔️:       atoms.insert(chiralSet->d_idx3);
    // RDKit✔️✔️:       atoms.insert(chiralSet->d_idx4);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::vector<int> atomsToCheck(atoms.begin(), atoms.end());
    // RDKit✔️✔️:   if (atomsToCheck.size() > 0) {
    // RDKit✔️✔️:     if (!_boundsFulfilled(atomsToCheck, *eargs.mmat, *positions)) {
    // RDKit✔️✔️:       if (embedParams.trackFailures) {
    // RDKit✔️✔️:         embedParams.failures[EmbedFailureCauses::FINAL_CHIRAL_BOUNDS]++;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // "center in volume" chirality test
    // RDKit✔️✔️:   for (const auto &chiralSet : *eargs.chiralCenters) {
    // RDKit✔️✔️:     // it could happen that the centroid is outside the volume defined
    // RDKit✔️✔️:     // by the other four points. That is also a fail.
    // RDKit✔️✔️:     if (!_centerInVolume(chiralSet, *positions)) {
    // RDKit✔️✔️:       if (embedParams.trackFailures) {
    // RDKit✔️✔️:         embedParams.failures[EmbedFailureCauses::FINAL_CENTER_IN_VOLUME]++;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::finalChiralChecks
    let trace_num_points = distgeom_trace_num_points_enabled(positions.len());
    if !embedder_check_chiral_centers(positions, eargs, embed_params) {
        if trace_num_points {
            println!(
                "distgeom_final_chiral n={} substage=check_chiral_centers2 ok=0",
                positions.len()
            );
        }
        embedder_increment_failure(embed_params, EmbedFailureCause::CheckChiralCenters2);
        return false;
    }

    let mut atoms = std::collections::BTreeSet::new();
    for chiral_set in eargs.chiral_centers {
        if chiral_set.idx0 != chiral_set.idx4 {
            atoms.insert(chiral_set.idx0 as i32);
            atoms.insert(chiral_set.idx1 as i32);
            atoms.insert(chiral_set.idx2 as i32);
            atoms.insert(chiral_set.idx3 as i32);
            atoms.insert(chiral_set.idx4 as i32);
        }
    }
    let atoms_to_check: Vec<i32> = atoms.into_iter().collect();
    if !atoms_to_check.is_empty() {
        let positions_3d = point_vectors_to_forcefield_vec3(positions);
        if !embedder_bounds_fulfilled(&atoms_to_check, eargs.mmat, &positions_3d) {
            if trace_num_points {
                println!(
                    "distgeom_final_chiral n={} substage=final_chiral_bounds ok=0 atoms={atoms_to_check:?}",
                    positions.len()
                );
            }
            embedder_increment_failure(embed_params, EmbedFailureCause::FinalChiralBounds);
            return false;
        }
    }

    let positions_3d = point_vectors_to_forcefield_vec3(positions);
    for chiral_set in eargs.chiral_centers {
        if !embedder_center_in_volume(chiral_set, &positions_3d, 0.1) {
            if trace_num_points {
                println!(
                    "distgeom_final_chiral n={} substage=final_center_in_volume ok=0 center={}",
                    positions.len(),
                    chiral_set.idx0
                );
            }
            embedder_increment_failure(embed_params, EmbedFailureCause::FinalCenterInVolume);
            return false;
        }
    }
    true
}

fn embedder_embed_points_with_rng<R: RdkitDoubleRng>(
    positions: &mut [Vec<f64>],
    eargs: &EmbedArgs<'_>,
    embed_params: &mut EmbedParameters,
    end_time: Option<Instant>,
    rng: &mut R,
) -> Result<bool, DgBoundsError> {
    let mut got_coords = false;
    let mut iter = 0_u32;
    let mut dist_mat = SymmMatrix::new(positions.len());
    let trace_row64 = row64_distgeom_trace_enabled(positions.len());
    let trace_row61 = row61_distgeom_trace_enabled(positions.len());
    let trace_num_points = distgeom_trace_num_points_enabled(positions.len());

    while !got_coords && iter < embed_params.max_iterations {
        if let Some(deadline) = end_time
            && Instant::now() > deadline
        {
            if trace_row64 {
                println!("row64_embed_points timeout_before_iter iter={iter}");
            } else if trace_row61 {
                println!("row61_embed_points timeout_before_iter iter={iter}");
            } else if trace_num_points {
                println!(
                    "distgeom_embed_points n={} timeout_before_iter iter={iter}",
                    positions.len()
                );
            }
            break;
        }

        iter += 1;
        if trace_row64 {
            println!("row64_embed_points iter_start iter={iter}");
        } else if trace_row61 {
            println!("row61_embed_points iter_start iter={iter}");
        } else if trace_num_points {
            println!(
                "distgeom_embed_points n={} iter_start iter={iter}",
                positions.len()
            );
        }
        if let Some(callback) = embed_params.callback {
            callback(iter);
        }

        got_coords =
            embedder_generate_initial_coords(positions, eargs, embed_params, &mut dist_mat, rng)?;
        if !got_coords {
            if trace_row64 {
                println!("row64_embed_points iter={iter} stage=initial_coords ok=0");
            } else if trace_row61 {
                println!("row61_embed_points iter={iter} stage=initial_coords ok=0");
            } else if trace_num_points {
                println!(
                    "distgeom_embed_points n={} iter={iter} stage=initial_coords ok=0",
                    positions.len()
                );
            }
            embedder_increment_failure(embed_params, EmbedFailureCause::InitialCoords);
        } else {
            #[cfg(test)]
            embedder_test_trace_update(|trace| {
                embedder_test_trace_set_coords(&mut trace.initial_coords, positions);
            });
            if trace_row64 {
                println!("row64_embed_points iter={iter} stage=initial_coords ok=1");
            } else if trace_row61 {
                println!("row61_embed_points iter={iter} stage=initial_coords ok=1");
            } else if trace_num_points {
                println!(
                    "distgeom_embed_points n={} iter={iter} stage=initial_coords ok=1",
                    positions.len()
                );
            }
            got_coords = embedder_first_minimization(positions, eargs, embed_params);
            if !got_coords {
                if trace_row64 {
                    println!("row64_embed_points iter={iter} stage=first_minimization ok=0");
                } else if trace_row61 {
                    println!("row61_embed_points iter={iter} stage=first_minimization ok=0");
                } else if trace_num_points {
                    println!(
                        "distgeom_embed_points n={} iter={iter} stage=first_minimization ok=0",
                        positions.len()
                    );
                }
                embedder_increment_failure(embed_params, EmbedFailureCause::FirstMinimization);
            } else {
                #[cfg(test)]
                embedder_test_trace_update(|trace| {
                    embedder_test_trace_set_coords(&mut trace.first_minimized, positions);
                });
                if trace_row64 {
                    println!("row64_embed_points iter={iter} stage=first_minimization ok=1");
                } else if trace_row61 {
                    println!("row61_embed_points iter={iter} stage=first_minimization ok=1");
                } else if trace_num_points {
                    println!(
                        "distgeom_embed_points n={} iter={iter} stage=first_minimization ok=1",
                        positions.len()
                    );
                }
                got_coords = embedder_check_tetrahedral_centers(positions, eargs, embed_params);
                if !got_coords {
                    embedder_increment_failure(
                        embed_params,
                        EmbedFailureCause::CheckTetrahedralCenters,
                    );
                }
            }

            if got_coords && embed_params.enforce_chirality && !eargs.chiral_centers.is_empty() {
                got_coords = embedder_check_chiral_centers(positions, eargs, embed_params);
                if !got_coords {
                    embedder_increment_failure(embed_params, EmbedFailureCause::CheckChiralCenters);
                }
            }

            if got_coords && (!eargs.chiral_centers.is_empty() || embed_params.use_random_coords) {
                got_coords =
                    embedder_minimize_fourth_dimension(positions, eargs, embed_params, end_time);
                if got_coords {
                    #[cfg(test)]
                    embedder_test_trace_update(|trace| {
                        embedder_test_trace_set_coords(
                            &mut trace.fourth_dimension_cleaned,
                            positions,
                        );
                    });
                }
                if trace_row64 {
                    println!(
                        "row64_embed_points iter={iter} stage=fourth_dimension ok={}",
                        i32::from(got_coords)
                    );
                } else if trace_row61 {
                    println!(
                        "row61_embed_points iter={iter} stage=fourth_dimension ok={}",
                        i32::from(got_coords)
                    );
                } else if trace_num_points {
                    println!(
                        "distgeom_embed_points n={} iter={iter} stage=fourth_dimension ok={}",
                        positions.len(),
                        i32::from(got_coords)
                    );
                }
                if !got_coords {
                    if let Some(deadline) = end_time
                        && Instant::now() > deadline
                    {
                        embedder_increment_failure(
                            embed_params,
                            EmbedFailureCause::ExceededTimeout,
                        );
                    }
                    embedder_increment_failure(
                        embed_params,
                        EmbedFailureCause::MinimizeFourthDimension,
                    );
                }
            }

            if got_coords
                && (embed_params.use_exp_torsion_angle_prefs || embed_params.use_basic_knowledge)
            {
                got_coords = embedder_minimize_with_exp_torsions(positions, eargs, embed_params);
                if got_coords {
                    #[cfg(test)]
                    embedder_test_trace_update(|trace| {
                        embedder_test_trace_set_coords(&mut trace.exp_torsion_minimized, positions);
                    });
                }
                if trace_num_points {
                    println!(
                        "distgeom_embed_points n={} iter={iter} stage=etk_minimization ok={}",
                        positions.len(),
                        i32::from(got_coords)
                    );
                }
                if !got_coords {
                    embedder_increment_failure(embed_params, EmbedFailureCause::EtkMinimization);
                }
            }

            if got_coords {
                got_coords =
                    embedder_double_bond_geometry_checks(positions, eargs, embed_params, 1.0e-3);
                if trace_num_points {
                    println!(
                        "distgeom_embed_points n={} iter={iter} stage=double_bond_geometry ok={}",
                        positions.len(),
                        i32::from(got_coords)
                    );
                }
                if !got_coords {
                    embedder_increment_failure(embed_params, EmbedFailureCause::LinearDoubleBond);
                }
            }

            if embed_params.enforce_chirality && got_coords {
                if !eargs.chiral_centers.is_empty() {
                    got_coords = embedder_final_chiral_checks(positions, eargs, embed_params);
                    if trace_num_points {
                        println!(
                            "distgeom_embed_points n={} iter={iter} stage=final_chiral ok={}",
                            positions.len(),
                            i32::from(got_coords)
                        );
                    }
                }
                if got_coords && !eargs.stereo_double_bonds.is_empty() {
                    got_coords = embedder_double_bond_stereo_checks(positions, eargs, embed_params);
                    if trace_num_points {
                        println!(
                            "distgeom_embed_points n={} iter={iter} stage=double_bond_stereo ok={}",
                            positions.len(),
                            i32::from(got_coords)
                        );
                    }
                    if !got_coords {
                        embedder_increment_failure(
                            embed_params,
                            EmbedFailureCause::BadDoubleBondStereo,
                        );
                    }
                }
            }
            if got_coords {
                #[cfg(test)]
                embedder_test_trace_update(|trace| {
                    embedder_test_trace_set_coords(&mut trace.final_checked, positions);
                });
            }
        }
    }

    Ok(got_coords)
}

fn embedder_embed_points(
    positions: &mut [Vec<f64>],
    eargs: EmbedArgs<'_>,
    embed_params: &mut EmbedParameters,
    seed: i32,
    end_time: Option<Instant>,
) -> Result<bool, DgBoundsError> {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::embedPoints (Embedder.cpp:838-1030)
    // RDKit✔️✔️: bool embedPoints(RDGeom::PointPtrVect *positions, detail::EmbedArgs eargs,
    // RDKit✔️✔️:                  EmbedParameters &embedParams, int seed, TimePoint *end_time) {
    // RDKit✔️✔️:   PRECONDITION(positions, "bogus positions");
    // RDKit✔️✔️:   if (embedParams.maxIterations == 0) {
    // RDKit✔️✔️:     embedParams.maxIterations = 10 * positions->size();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   RDNumeric::DoubleSymmMatrix distMat(positions->size(), 0.0);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // The basin threshold just gets us into trouble when we're using
    // RDKit✔️✔️:   // random coordinates since it ends up ignoring 1-4 (and higher)
    // RDKit✔️✔️:   // interactions. This causes us to get folded-up (and self-penetrating)
    // RDKit✔️✔️:   // conformations for large flexible molecules
    // RDKit✔️✔️:   if (embedParams.useRandomCoords) {
    // RDKit✔️✔️:     embedParams.basinThresh = 1e8;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::unique_ptr<RDKit::rng_type> generator;
    // RDKit✔️✔️:   std::unique_ptr<RDKit::uniform_double> distrib;
    // RDKit✔️✔️:   std::unique_ptr<RDKit::double_source_type> rngMgr;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   RDKit::double_source_type *rng = nullptr;
    // RDKit✔️✔️:   CHECK_INVARIANT(seed >= -1,
    // RDKit✔️✔️:                   "random seed must either be positive, zero, or negative one");
    // RDKit✔️✔️:   if (seed > -1) {
    // RDKit✔️✔️:     generator.reset(new RDKit::rng_type(42u));
    // RDKit✔️✔️:     generator->seed(seed);
    // RDKit✔️✔️:     distrib.reset(new RDKit::uniform_double(0.0, 1.0));
    // RDKit✔️✔️:     rngMgr.reset(new RDKit::double_source_type(*generator, *distrib));
    // RDKit✔️✔️:     rng = rngMgr.get();
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     rng = &RDKit::getDoubleRandomSource();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   bool gotCoords = false;
    // RDKit✔️✔️:   unsigned int iter = 0;
    // RDKit✔️✔️:   while (!gotCoords && iter < embedParams.maxIterations) {
    // RDKit✔️✔️:     if (end_time != nullptr && Clock::now() > *end_time) {
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     ++iter;
    // RDKit✔️✔️:     if (embedParams.callback != nullptr) {
    // RDKit✔️✔️:       embedParams.callback(iter);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (ControlCHandler::getGotSignal()) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     gotCoords = EmbeddingOps::generateInitialCoords(positions, eargs,
    // RDKit✔️✔️:                                                     embedParams, distMat, rng);
    // RDKit✔️✔️:     if (!gotCoords) {
    // RDKit✔️✔️:       if (embedParams.trackFailures) {
    // RDKit✔️✔️:         embedParams.failures[EmbedFailureCauses::INITIAL_COORDS]++;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       gotCoords =
    // RDKit✔️✔️:           EmbeddingOps::firstMinimization(positions, eargs, embedParams);
    // RDKit✔️✔️:       if (!gotCoords) {
    // RDKit✔️✔️:         if (embedParams.trackFailures) {
    // RDKit✔️✔️:           embedParams.failures[EmbedFailureCauses::FIRST_MINIMIZATION]++;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         gotCoords = EmbeddingOps::checkTetrahedralCenters(positions, eargs,
    // RDKit✔️✔️:                                                           embedParams);
    // RDKit✔️✔️:         if (!gotCoords) {
    // RDKit✔️✔️:           if (embedParams.trackFailures) {
    // RDKit✔️✔️:             embedParams
    // RDKit✔️✔️:                 .failures[EmbedFailureCauses::CHECK_TETRAHEDRAL_CENTERS]++;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (gotCoords && embedParams.enforceChirality &&
    // RDKit✔️✔️:           eargs.chiralCenters->size() > 0) {
    // RDKit✔️✔️:         gotCoords =
    // RDKit✔️✔️:             EmbeddingOps::checkChiralCenters(positions, eargs, embedParams);
    // RDKit✔️✔️:         if (!gotCoords) {
    // RDKit✔️✔️:           if (embedParams.trackFailures) {
    // RDKit✔️✔️:             embedParams.failures[EmbedFailureCauses::CHECK_CHIRAL_CENTERS]++;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (gotCoords &&
    // RDKit✔️✔️:           (eargs.chiralCenters->size() > 0 || embedParams.useRandomCoords)) {
    // RDKit✔️✔️:         gotCoords = EmbeddingOps::minimizeFourthDimension(
    // RDKit✔️✔️:             positions, eargs, embedParams, end_time);
    // RDKit✔️✔️:         if (!gotCoords) {
    // RDKit✔️✔️:           if (embedParams.trackFailures) {
    // RDKit✔️✔️:             if (end_time != nullptr && Clock::now() > *end_time) {
    // RDKit✔️✔️:               embedParams.failures[EmbedFailureCauses::EXCEEDED_TIMEOUT]++;
    // RDKit✔️✔️:             }
    // RDKit✔️✔️:             embedParams
    // RDKit✔️✔️:                 .failures[EmbedFailureCauses::MINIMIZE_FOURTH_DIMENSION]++;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (gotCoords && (embedParams.useExpTorsionAnglePrefs ||
    // RDKit✔️✔️:                         embedParams.useBasicKnowledge)) {
    // RDKit✔️✔️:         gotCoords = EmbeddingOps::minimizeWithExpTorsions(*positions, eargs,
    // RDKit✔️✔️:                                                           embedParams);
    // RDKit✔️✔️:         if (!gotCoords) {
    // RDKit✔️✔️:           if (embedParams.trackFailures) {
    // RDKit✔️✔️:             embedParams.failures[EmbedFailureCauses::ETK_MINIMIZATION]++;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (gotCoords) {
    // RDKit✔️✔️:         gotCoords = EmbeddingOps::doubleBondGeometryChecks(*positions, eargs,
    // RDKit✔️✔️:                                                            embedParams);
    // RDKit✔️✔️:         if (!gotCoords && embedParams.trackFailures) {
    // RDKit✔️✔️:           embedParams.failures[EmbedFailureCauses::LINEAR_DOUBLE_BOND]++;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (embedParams.enforceChirality && gotCoords) {
    // RDKit✔️✔️:         if (!eargs.chiralCenters->empty()) {
    // RDKit✔️✔️:           gotCoords =
    // RDKit✔️✔️:               EmbeddingOps::finalChiralChecks(positions, eargs, embedParams);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         if (gotCoords && !eargs.stereoDoubleBonds->empty()) {
    // RDKit✔️✔️:           gotCoords = EmbeddingOps::doubleBondStereoChecks(*positions, eargs,
    // RDKit✔️✔️:                                                            embedParams);
    // RDKit✔️✔️:           if (!gotCoords && embedParams.trackFailures) {
    // RDKit✔️✔️:             embedParams.failures[EmbedFailureCauses::BAD_DOUBLE_BOND_STEREO]++;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return gotCoords;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::embedPoints
    assert!(
        seed >= -1,
        "random seed must either be positive, zero, or negative one"
    );
    if embed_params.max_iterations == 0 {
        embed_params.max_iterations = 10 * positions.len() as u32;
    }
    if embed_params.use_random_coords {
        embed_params.basin_thresh = 1.0e8;
    }

    let got_coords = if seed > -1 {
        let mut rng = RdkitDistgeomMinStdRand::new_from_embed_points_seed(seed);
        embedder_embed_points_with_rng(positions, &eargs, embed_params, end_time, &mut rng)
    } else {
        RDKIT_DISTGEOM_RNG.with(|rng| {
            embedder_embed_points_with_rng(
                positions,
                &eargs,
                embed_params,
                end_time,
                &mut *rng.borrow_mut(),
            )
        })
    }?;
    Ok(got_coords)
}

fn embedder_find_double_bonds(
    mol: &Molecule,
    double_bond_ends: &mut Vec<(usize, usize, usize)>,
    stereo_double_bonds: &mut Vec<(Vec<usize>, i32)>,
    coord_map: Option<&BTreeMap<i32, ForceFieldVec3>>,
) {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::findDoubleBonds (Embedder.cpp:1027-1077)
    // RDKit✔️✔️: RDKIT_DISTGEOMHELPERS_EXPORT void findDoubleBonds(
    // RDKit✔️✔️:     const ROMol &mol,
    // RDKit✔️✔️:     std::vector<std::tuple<unsigned int, unsigned int, unsigned int>>
    // RDKit✔️✔️:         &doubleBondEnds,
    // RDKit✔️✔️:     std::vector<std::pair<std::vector<unsigned int>, int>> &stereoDoubleBonds,
    // RDKit✔️✔️:     const std::map<int, RDGeom::Point3D> *coordMap) {
    // RDKit✔️✔️:   doubleBondEnds.clear();
    // RDKit✔️✔️:   stereoDoubleBonds.clear();
    // RDKit✔️✔️:   for (const auto bnd : mol.bonds()) {
    // RDKit✔️✔️:     if (bnd->getBondType() == Bond::BondType::DOUBLE) {
    // RDKit✔️✔️:       for (const auto atm : {bnd->getBeginAtom(), bnd->getEndAtom()}) {
    // RDKit✔️✔️:         if (atm->getDegree() < 2) {
    // RDKit✔️✔️:           continue;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         auto oatm = bnd->getOtherAtom(atm);
    // RDKit✔️✔️:         for (const auto nbr : mol.atomNeighbors(atm)) {
    // RDKit✔️✔️:           if (nbr == oatm) {
    // RDKit✔️✔️:             continue;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           const auto obnd =
    // RDKit✔️✔️:               mol.getBondBetweenAtoms(atm->getIdx(), nbr->getIdx());
    // RDKit✔️✔️:           if (!obnd || (obnd->getBondType() != Bond::BondType::SINGLE &&
    // RDKit✔️✔️:                         atm->getDegree() == 2)) {
    // RDKit✔️✔️:             continue;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           doubleBondEnds.emplace_back(nbr->getIdx(), atm->getIdx(),
    // RDKit✔️✔️:                                       oatm->getIdx());
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // if there's stereo, handle that too:
    // RDKit✔️✔️:       if (bnd->getStereo() > Bond::BondStereo::STEREOANY) {
    // RDKit✔️✔️:         // only do this if the controlling atoms aren't in the coord map
    // RDKit✔️✔️:         if (coordMap &&
    // RDKit✔️✔️:             coordMap->find(bnd->getStereoAtoms()[0]) != coordMap->end() &&
    // RDKit✔️✔️:             coordMap->find(bnd->getStereoAtoms()[1]) != coordMap->end()) {
    // RDKit✔️✔️:           continue;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         int sign = 1;
    // RDKit✔️✔️:         if (bnd->getStereo() == Bond::BondStereo::STEREOCIS ||
    // RDKit✔️✔️:             bnd->getStereo() == Bond::BondStereo::STEREOZ) {
    // RDKit✔️✔️:           sign = -1;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         std::pair<std::vector<unsigned int>, int> elem{
    // RDKit✔️✔️:             {static_cast<unsigned>(bnd->getStereoAtoms()[0]),
    // RDKit✔️✔️:              bnd->getBeginAtomIdx(), bnd->getEndAtomIdx(),
    // RDKit✔️✔️:              static_cast<unsigned>(bnd->getStereoAtoms()[1])},
    // RDKit✔️✔️:             sign};
    // RDKit✔️✔️:         stereoDoubleBonds.push_back(elem);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::findDoubleBonds
    double_bond_ends.clear();
    stereo_double_bonds.clear();

    for bnd in mol.bonds() {
        if bnd.order() != BondOrder::Double {
            continue;
        }
        for atm in [bnd.begin().index(), bnd.end().index()] {
            let atm_degree = neighbors_for_atom(mol, atm).len();
            if atm_degree < 2 {
                continue;
            }
            let oatm = if atm == bnd.begin().index() {
                bnd.end().index()
            } else {
                bnd.begin().index()
            };
            for nbr in neighbors_for_atom(mol, atm) {
                if nbr == oatm {
                    continue;
                }
                let Some(obnd) = bond_between(mol, atm, nbr) else {
                    continue;
                };
                if obnd.order() != BondOrder::Single && atm_degree == 2 {
                    continue;
                }
                double_bond_ends.push((nbr, atm, oatm));
            }
        }

        if !matches!(bnd.stereo(), BondStereo::None | BondStereo::Any) {
            let stereo_atoms = bnd
                .stereo_atoms()
                .expect("stereo double bond requires two controlling atoms");
            if let Some(coord_map) = coord_map
                && coord_map.contains_key(&(stereo_atoms[0].index() as i32))
                && coord_map.contains_key(&(stereo_atoms[1].index() as i32))
            {
                continue;
            }
            let sign = if matches!(bnd.stereo(), BondStereo::Cis | BondStereo::Z) {
                -1
            } else {
                1
            };
            stereo_double_bonds.push((
                vec![
                    stereo_atoms[0].index(),
                    bnd.begin().index(),
                    bnd.end().index(),
                    stereo_atoms[1].index(),
                ],
                sign,
            ));
        }
    }
}

fn embedder_find_chiral_sets(
    mol: &Molecule,
    chiral_centers: &mut VectChiralSet,
    tetrahedral_centers: &mut VectChiralSet,
    coord_map: Option<&BTreeMap<i32, ForceFieldVec3>>,
) {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::findChiralSets (Embedder.cpp:1079-1219)
    // RDKit✔️✔️: void findChiralSets(const ROMol &mol, DistGeom::VECT_CHIRALSET &chiralCenters,
    // RDKit✔️✔️:                     DistGeom::VECT_CHIRALSET &tetrahedralCenters,
    // RDKit✔️✔️:                     const std::map<int, RDGeom::Point3D> *coordMap) {
    // RDKit✔️✔️:   for (const auto &atom : mol.atoms()) {
    // RDKit✔️✔️:     if (atom->getAtomicNum() != 1) {  // skip hydrogens
    // RDKit✔️✔️:       Atom::ChiralType chiralType = atom->getChiralTag();
    // RDKit✔️✔️:       if ((chiralType == Atom::CHI_TETRAHEDRAL_CW ||
    // RDKit✔️✔️:            chiralType == Atom::CHI_TETRAHEDRAL_CCW) ||
    // RDKit✔️✔️:           ((atom->getAtomicNum() == 6 || atom->getAtomicNum() == 7) &&
    // RDKit✔️✔️:            atom->getDegree() == 4)) {
    // RDKit✔️✔️:         // make a chiral set from the neighbors
    // RDKit✔️✔️:         INT_VECT nbrs;
    // RDKit✔️✔️:         nbrs.reserve(4);
    // RDKit✔️✔️:         // find the neighbors of this atom and enter them into the
    // RDKit✔️✔️:         // nbr list
    // RDKit✔️✔️:         ROMol::OEDGE_ITER beg, end;
    // RDKit✔️✔️:         boost::tie(beg, end) = mol.getAtomBonds(atom);
    // RDKit✔️✔️:         while (beg != end) {
    // RDKit✔️✔️:           nbrs.push_back(mol[*beg]->getOtherAtom(atom)->getIdx());
    // RDKit✔️✔️:           ++beg;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         // if we have less than 4 heavy atoms as neighbors,
    // RDKit✔️✔️:         // we need to include the chiral center into the mix
    // RDKit✔️✔️:         // we should at least have 3 though
    // RDKit✔️✔️:         CHECK_INVARIANT(nbrs.size() >= 3, "Cannot be a chiral center");
    // RDKit✔️✔️:
    // RDKit✔️✔️:         double volLowerBound = 5.0;
    // RDKit✔️✔️:         double volUpperBound = 100.0;
    // RDKit✔️✔️:
    // RDKit✔️✔️:         if (nbrs.size() < 4) {
    // RDKit✔️✔️:           // we get lower volumes if there are three neighbors,
    // RDKit✔️✔️:           //  this was github #5883
    // RDKit✔️✔️:           volLowerBound = 2.0;
    // RDKit✔️✔️:           nbrs.insert(nbrs.end(), atom->getIdx());
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:
    // RDKit✔️✔️:         // set a flag for tetrahedral centers that are in multiple small rings
    // RDKit✔️✔️:         auto numSmallRings = 0u;
    // RDKit✔️✔️:         constexpr int smallRingSize = 5;
    // RDKit✔️✔️:         for (const auto sz : mol.getRingInfo()->atomRingSizes(atom->getIdx())) {
    // RDKit✔️✔️:           if (sz < smallRingSize) {
    // RDKit✔️✔️:             ++numSmallRings;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         std::uint64_t structureFlags = 0;
    // RDKit✔️✔️:         if (numSmallRings > 1) {
    // RDKit✔️✔️:           structureFlags = static_cast<std::uint64_t>(
    // RDKit✔️✔️:               DistGeom::ChiralSetStructureFlags::IN_FUSED_SMALL_RINGS);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:
    // RDKit✔️✔️:         // now create a chiral set and set the upper and lower bound on the
    // RDKit✔️✔️:         // volume
    // RDKit✔️✔️:         if (chiralType == Atom::CHI_TETRAHEDRAL_CCW) {
    // RDKit✔️✔️:           // positive chiral volume
    // RDKit✔️✔️:           auto *cset = new DistGeom::ChiralSet(atom->getIdx(), nbrs[0], nbrs[1],
    // RDKit✔️✔️:                                                nbrs[2], nbrs[3], volLowerBound,
    // RDKit✔️✔️:                                                volUpperBound, structureFlags);
    // RDKit✔️✔️:           DistGeom::ChiralSetPtr cptr(cset);
    // RDKit✔️✔️:           chiralCenters.push_back(cptr);
    // RDKit✔️✔️:         } else if (chiralType == Atom::CHI_TETRAHEDRAL_CW) {
    // RDKit✔️✔️:           auto *cset = new DistGeom::ChiralSet(atom->getIdx(), nbrs[0], nbrs[1],
    // RDKit✔️✔️:                                                nbrs[2], nbrs[3], -volUpperBound,
    // RDKit✔️✔️:                                                -volLowerBound, structureFlags);
    // RDKit✔️✔️:           DistGeom::ChiralSetPtr cptr(cset);
    // RDKit✔️✔️:           chiralCenters.push_back(cptr);
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           if ((coordMap && coordMap->find(atom->getIdx()) != coordMap->end()) ||
    // RDKit✔️✔️:               (mol.getRingInfo()->isInitialized() &&
    // RDKit✔️✔️:                (mol.getRingInfo()->numAtomRings(atom->getIdx()) < 2 ||
    // RDKit✔️✔️:                 mol.getRingInfo()->isAtomInRingOfSize(atom->getIdx(), 3)))) {
    // RDKit✔️✔️:             // we only want to these tests for ring atoms that are not part of
    // RDKit✔️✔️:             // the coordMap
    // RDKit✔️✔️:             // there's no sense doing 3-rings because those are a nightmare
    // RDKit✔️✔️:           } else {
    // RDKit✔️✔️:             auto *cset = new DistGeom::ChiralSet(atom->getIdx(), nbrs[0],
    // RDKit✔️✔️:                                                  nbrs[1], nbrs[2], nbrs[3], 0.0,
    // RDKit✔️✔️:                                                  0.0, structureFlags);
    // RDKit✔️✔️:             DistGeom::ChiralSetPtr cptr(cset);
    // RDKit✔️✔️:             tetrahedralCenters.push_back(cptr);
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }  // if block -chirality check
    // RDKit✔️✔️:     }    // if block - heavy atom check
    // RDKit✔️✔️:   }      // for loop over atoms
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // now do atropisomers
    // RDKit✔️✔️:   for (const auto &bond : mol.bonds()) {
    // RDKit✔️✔️:     if (bond->getStereo() != Bond::BondStereo::STEREOATROPCCW &&
    // RDKit✔️✔️:         bond->getStereo() != Bond::BondStereo::STEREOATROPCW) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     Atropisomers::AtropAtomAndBondVec atomsAndBonds[2];
    // RDKit✔️✔️:     Atropisomers::getAtropisomerAtomsAndBonds(bond, atomsAndBonds, mol);
    // RDKit✔️✔️:     // make a chiral set for the atropisomeric bond
    // RDKit✔️✔️:     // we start with only managing cases where there are two exo-substituents on
    // RDKit✔️✔️:     // at least one side
    // RDKit✔️✔️:     if (atomsAndBonds[0].second.size() != 2 &&
    // RDKit✔️✔️:         atomsAndBonds[1].second.size() != 2) {
    // RDKit✔️✔️:       BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:           << "Atropisomer bond stereochemistry not used for bond "
    // RDKit✔️✔️:           << bond->getIdx()
    // RDKit✔️✔️:           << ", which does not have two exo substituents on at least one side."
    // RDKit✔️✔️:           << std::endl;
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     int idx0 = atomsAndBonds[0].first->getIdx();
    // RDKit✔️✔️:     int idx1 = atomsAndBonds[1].first->getIdx();
    // RDKit✔️✔️:
    // RDKit✔️✔️:     int nbr1 = atomsAndBonds[0].second[0]->getOtherAtomIdx(idx0);
    // RDKit✔️✔️:     int nbr2 = 0;
    // RDKit✔️✔️:     int nbr3 = 0;
    // RDKit✔️✔️:     int nbr4 = 0;
    // RDKit✔️✔️:     if (atomsAndBonds[0].second.size() == 2) {
    // RDKit✔️✔️:       nbr2 = atomsAndBonds[0].second[1]->getOtherAtomIdx(idx0);
    // RDKit✔️✔️:       nbr3 = atomsAndBonds[1].second[0]->getOtherAtomIdx(idx1);
    // RDKit✔️✔️:       if (atomsAndBonds[1].second.size() == 2) {
    // RDKit✔️✔️:         nbr4 = atomsAndBonds[1].second[1]->getOtherAtomIdx(idx1);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         nbr4 = idx0;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       nbr2 = atomsAndBonds[1].second[0]->getOtherAtomIdx(idx1);
    // RDKit✔️✔️:       nbr3 = atomsAndBonds[1].second[1]->getOtherAtomIdx(idx1);
    // RDKit✔️✔️:       nbr4 = idx0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     INT_VECT nbrs = {nbr1, nbr2, nbr3, nbr4};
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // FIX: these numbers are empirical and should be revisited
    // RDKit✔️✔️:     double volLowerBound = 1.0;
    // RDKit✔️✔️:     double volUpperBound = 100.0;
    // RDKit✔️✔️:     if (bond->getStereo() == Bond::BondStereo::STEREOATROPCCW) {
    // RDKit✔️✔️:       std::swap(volLowerBound, volUpperBound);
    // RDKit✔️✔️:       volLowerBound *= -1;
    // RDKit✔️✔️:       volUpperBound *= -1;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     auto *cset = new DistGeom::ChiralSet(idx0, nbrs[0], nbrs[1], nbrs[2],
    // RDKit✔️✔️:                                          nbrs[3], volLowerBound, volUpperBound);
    // RDKit✔️✔️:     DistGeom::ChiralSetPtr cptr(cset);
    // RDKit✔️✔️:     chiralCenters.push_back(cptr);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::findChiralSets
    let ring_info = mol.derived_cache().rings.as_ref();
    for atom in mol.atoms() {
        if atom.atomic_number() == 1 {
            continue;
        }
        let chiral_type = atom.chiral_tag();
        let atom_idx = atom.id().index();
        let nbrs0 = neighbors_for_atom(mol, atom_idx);
        if matches!(
            chiral_type,
            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
        ) || ((atom.atomic_number() == 6 || atom.atomic_number() == 7) && nbrs0.len() == 4)
        {
            let mut nbrs = nbrs0;
            assert!(nbrs.len() >= 3, "Cannot be a chiral center");
            let mut vol_lower_bound = 5.0;
            let mut vol_upper_bound = 100.0;
            if nbrs.len() < 4 {
                vol_lower_bound = 2.0;
                nbrs.push(atom_idx);
            }

            let num_small_rings = ring_info
                .map(|rings| {
                    atom_ring_sizes(rings, atom_idx)
                        .into_iter()
                        .filter(|size| *size < 5)
                        .count()
                })
                .unwrap_or(0);
            let structure_flags = if num_small_rings > 1 {
                ChiralSetStructureFlags::InFusedSmallRings as u64
            } else {
                0
            };

            if chiral_type == ChiralTag::TetrahedralCcw {
                chiral_centers.push(Arc::new(ChiralSet::new(
                    atom_idx,
                    nbrs[0],
                    nbrs[1],
                    nbrs[2],
                    nbrs[3],
                    vol_lower_bound,
                    vol_upper_bound,
                    structure_flags,
                )));
            } else if chiral_type == ChiralTag::TetrahedralCw {
                chiral_centers.push(Arc::new(ChiralSet::new(
                    atom_idx,
                    nbrs[0],
                    nbrs[1],
                    nbrs[2],
                    nbrs[3],
                    -vol_upper_bound,
                    -vol_lower_bound,
                    structure_flags,
                )));
            } else {
                let coord_mapped =
                    coord_map.is_some_and(|map| map.contains_key(&(atom_idx as i32)));
                let ring_excluded = ring_info.is_some_and(|rings| {
                    rings.num_atom_rings(AtomId::new(atom_idx)) < 2
                        || rings.is_atom_in_ring_of_size(AtomId::new(atom_idx), 3)
                });
                if !(coord_mapped || ring_excluded) {
                    tetrahedral_centers.push(Arc::new(ChiralSet::new(
                        atom_idx,
                        nbrs[0],
                        nbrs[1],
                        nbrs[2],
                        nbrs[3],
                        0.0,
                        0.0,
                        structure_flags,
                    )));
                }
            }
        }
    }

    for bond in mol.bonds() {
        if !matches!(bond.stereo(), BondStereo::AtropCcw | BondStereo::AtropCw) {
            continue;
        }
        let idx0 = bond.begin().index();
        let idx1 = bond.end().index();
        let atoms_and_bonds0 = embedder_atropisomer_neighbor_bonds(mol, idx0, bond.id().index());
        let atoms_and_bonds1 = embedder_atropisomer_neighbor_bonds(mol, idx1, bond.id().index());
        if atoms_and_bonds0.len() != 2 && atoms_and_bonds1.len() != 2 {
            continue;
        }
        let nbr1 = embedder_bond_other_atom_idx(mol, atoms_and_bonds0[0], idx0);
        let (nbr2, nbr3, nbr4) = if atoms_and_bonds0.len() == 2 {
            let nbr2 = embedder_bond_other_atom_idx(mol, atoms_and_bonds0[1], idx0);
            let nbr3 = embedder_bond_other_atom_idx(mol, atoms_and_bonds1[0], idx1);
            let nbr4 = if atoms_and_bonds1.len() == 2 {
                embedder_bond_other_atom_idx(mol, atoms_and_bonds1[1], idx1)
            } else {
                idx0
            };
            (nbr2, nbr3, nbr4)
        } else {
            let nbr2 = embedder_bond_other_atom_idx(mol, atoms_and_bonds1[0], idx1);
            let nbr3 = embedder_bond_other_atom_idx(mol, atoms_and_bonds1[1], idx1);
            (nbr2, nbr3, idx0)
        };
        let (vol_lower_bound, vol_upper_bound) = if bond.stereo() == BondStereo::AtropCcw {
            (-100.0, -1.0)
        } else {
            (1.0, 100.0)
        };
        chiral_centers.push(Arc::new(ChiralSet::with_default_structure_flags(
            idx0,
            nbr1,
            nbr2,
            nbr3,
            nbr4,
            vol_lower_bound,
            vol_upper_bound,
        )));
    }
}

fn embedder_atropisomer_neighbor_bonds(
    mol: &Molecule,
    focus_atom: usize,
    atrop_bond: usize,
) -> Vec<usize> {
    let mut bonds: Vec<usize> = mol
        .bonds()
        .iter()
        .filter(|bond| {
            bond.id().index() != atrop_bond
                && (bond.begin().index() == focus_atom || bond.end().index() == focus_atom)
        })
        .map(|bond| bond.id().index())
        .collect();
    if bonds.len() == 2 {
        let other0 = embedder_bond_other_atom_idx(mol, bonds[0], focus_atom);
        let other1 = embedder_bond_other_atom_idx(mol, bonds[1], focus_atom);
        if other1 < other0 {
            bonds.swap(0, 1);
        }
    }
    bonds
}

fn embedder_bond_other_atom_idx(mol: &Molecule, bond_idx: usize, atom_idx: usize) -> usize {
    let bond = &mol.bonds()[bond_idx];
    if bond.begin().index() == atom_idx {
        bond.end().index()
    } else {
        assert_eq!(bond.end().index(), atom_idx, "atom is not an endpoint");
        bond.begin().index()
    }
}

fn embedder_adjust_bounds_mat_from_coord_map(
    mmat: &mut BoundsMatrix,
    _num_atoms: usize,
    coord_map: &BTreeMap<i32, ForceFieldVec3>,
) -> Result<(), DgBoundsError> {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::adjustBoundsMatFromCoordMap (Embedder.cpp:1221-1236)
    // RDKit✔️✔️: void adjustBoundsMatFromCoordMap(
    // RDKit✔️✔️:     DistGeom::BoundsMatPtr mmat, unsigned int,
    // RDKit✔️✔️:     const std::map<int, RDGeom::Point3D> *coordMap) {
    // RDKit✔️✔️:   for (auto iIt = coordMap->begin(); iIt != coordMap->end(); ++iIt) {
    // RDKit✔️✔️:     unsigned int iIdx = iIt->first;
    // RDKit✔️✔️:     const RDGeom::Point3D &iPoint = iIt->second;
    // RDKit✔️✔️:     auto jIt = iIt;
    // RDKit✔️✔️:     while (++jIt != coordMap->end()) {
    // RDKit✔️✔️:       unsigned int jIdx = jIt->first;
    // RDKit✔️✔️:       const RDGeom::Point3D &jPoint = jIt->second;
    // RDKit✔️✔️:       double dist = (iPoint - jPoint).length();
    // RDKit✔️✔️:       mmat->setUpperBound(iIdx, jIdx, dist);
    // RDKit✔️✔️:       mmat->setLowerBound(iIdx, jIdx, dist);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::adjustBoundsMatFromCoordMap
    let entries: Vec<(usize, ForceFieldVec3)> = coord_map
        .iter()
        .map(|(idx, point)| {
            let idx = usize::try_from(*idx).expect("coordMap atom index must be non-negative");
            (idx, *point)
        })
        .collect();
    for i in 0..entries.len() {
        let (i_idx, i_point) = entries[i];
        for &(j_idx, j_point) in entries.iter().skip(i + 1) {
            let dist = (i_point - j_point).length();
            mmat.set_upper(i_idx, j_idx, dist)?;
            mmat.set_lower(i_idx, j_idx, dist)?;
        }
    }
    Ok(())
}

fn embedder_init_etkdg(
    mol: &Molecule,
    params: &EmbedParameters,
    etkdg_details: &mut CrystalFFDetails,
) -> Result<(), DgBoundsError> {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::initETKDG (Embedder.cpp:1238-1253)
    // RDKit✔️✔️: void initETKDG(ROMol *mol, const EmbedParameters &params,
    // RDKit✔️✔️:                ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails) {
    // RDKit✔️✔️:   PRECONDITION(mol, "bad molecule");
    // RDKit✔️✔️:   unsigned int nAtoms = mol->getNumAtoms();
    // RDKit✔️✔️:   if (params.useExpTorsionAnglePrefs || params.useBasicKnowledge) {
    // RDKit✔️✔️:     ForceFields::CrystalFF::getExperimentalTorsions(
    // RDKit✔️✔️:         *mol, etkdgDetails, params.useExpTorsionAnglePrefs,
    // RDKit✔️✔️:         params.useSmallRingTorsions, params.useMacrocycleTorsions,
    // RDKit✔️✔️:         params.useBasicKnowledge, params.ETversion, params.verbose);
    // RDKit✔️✔️:     etkdgDetails.atomNums.resize(nAtoms);
    // RDKit✔️✔️:     for (unsigned int i = 0; i < nAtoms; ++i) {
    // RDKit✔️✔️:       etkdgDetails.atomNums[i] = mol->getAtomWithIdx(i)->getAtomicNum();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   etkdgDetails.boundsMatForceScaling = params.boundsMatForceScaling;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::initETKDG
    let n_atoms = mol.num_atoms();
    let trace_row64 = row64_distgeom_trace_enabled(n_atoms);
    if params.use_exp_torsion_angle_prefs || params.use_basic_knowledge {
        let start = trace_row64.then(Instant::now);
        get_experimental_torsions_without_bonds(
            mol,
            etkdg_details,
            params.use_exp_torsion_angle_prefs,
            params.use_small_ring_torsions,
            params.use_macrocycle_torsions,
            params.use_basic_knowledge,
            params.et_version,
            params.verbose,
        )
        .map_err(|err| DgBoundsError::GenerationFailed(err.to_string()))?;
        if let Some(start) = start {
            println!(
                "row64_init_etkdg get_experimental_torsions={:.6}",
                start.elapsed().as_secs_f64()
            );
        }
        etkdg_details.atom_nums.resize(n_atoms, 0);
        for i in 0..n_atoms {
            etkdg_details.atom_nums[i] = i32::from(mol.atoms()[i].atomic_number());
        }
    }
    etkdg_details.bounds_mat_force_scaling = params.bounds_mat_force_scaling;
    Ok(())
}

fn embedder_setup_initial_bounds_matrix(
    mol: &Molecule,
    mmat: &mut BoundsMatrix,
    coord_map: Option<&BTreeMap<i32, ForceFieldVec3>>,
    params: &EmbedParameters,
    etkdg_details: &mut CrystalFFDetails,
) -> Result<bool, DgBoundsError> {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::setupInitialBoundsMatrix (Embedder.cpp:1255-1306)
    // RDKit✔️✔️: bool setupInitialBoundsMatrix(
    // RDKit✔️✔️:     ROMol *mol, DistGeom::BoundsMatPtr mmat,
    // RDKit✔️✔️:     const std::map<int, RDGeom::Point3D> *coordMap,
    // RDKit✔️✔️:     const EmbedParameters &params,
    // RDKit✔️✔️:     ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails) {
    // RDKit✔️✔️:   PRECONDITION(mol, "bad molecule");
    // RDKit✔️✔️:   unsigned int nAtoms = mol->getNumAtoms();
    // RDKit✔️✔️:   if (params.useExpTorsionAnglePrefs || params.useBasicKnowledge) {
    // RDKit✔️✔️:     setTopolBounds(*mol, mmat, etkdgDetails.bonds, etkdgDetails.angles, true,
    // RDKit✔️✔️:                    false, params.useMacrocycle14config,
    // RDKit✔️✔️:                    params.forceTransAmides);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     setTopolBounds(*mol, mmat, true, false, params.useMacrocycle14config,
    // RDKit✔️✔️:                    params.forceTransAmides);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double tol = 0.0;
    // RDKit✔️✔️:   if (coordMap) {
    // RDKit✔️✔️:     adjustBoundsMatFromCoordMap(mmat, nAtoms, coordMap);
    // RDKit✔️✔️:     tol = 0.05;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!DistGeom::triangleSmoothBounds(mmat, tol)) {
    // RDKit✔️✔️:     // ok this bound matrix failed to triangle smooth - re-compute the
    // RDKit✔️✔️:     // bounds matrix without 15 bounds and with VDW scaling
    // RDKit✔️✔️:     initBoundsMat(mmat);
    // RDKit✔️✔️:     setTopolBounds(*mol, mmat, false, true, params.useMacrocycle14config,
    // RDKit✔️✔️:                    params.forceTransAmides);
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (coordMap) {
    // RDKit✔️✔️:       adjustBoundsMatFromCoordMap(mmat, nAtoms, coordMap);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // try triangle smoothing again
    // RDKit✔️✔️:     if (!DistGeom::triangleSmoothBounds(mmat, tol)) {
    // RDKit✔️✔️:       // ok, we're not going to be able to smooth this,
    // RDKit✔️✔️:       if (params.ignoreSmoothingFailures) {
    // RDKit✔️✔️:         // proceed anyway with the more relaxed bounds matrix
    // RDKit✔️✔️:         initBoundsMat(mmat);
    // RDKit✔️✔️:         setTopolBounds(*mol, mmat, false, true, params.useMacrocycle14config,
    // RDKit✔️✔️:                        params.forceTransAmides);
    // RDKit✔️✔️:
    // RDKit✔️✔️:         if (coordMap) {
    // RDKit✔️✔️:           adjustBoundsMatFromCoordMap(mmat, nAtoms, coordMap);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:             << "Could not triangle bounds smooth molecule." << std::endl;
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::EmbeddingOps::setupInitialBoundsMatrix
    let n_atoms = mol.num_atoms();
    let trace_row64 = row64_distgeom_trace_enabled(n_atoms);
    if params.use_exp_torsion_angle_prefs || params.use_basic_knowledge {
        let start = trace_row64.then(Instant::now);
        set_topol_bounds_with_outputs(
            mol,
            mmat,
            &mut etkdg_details.bonds,
            &mut etkdg_details.angles,
            true,
            false,
            params.use_macrocycle14config,
            params.force_trans_amides,
            true,
            true,
        )?;
        if let Some(start) = start {
            println!(
                "row64_setup_bounds set_topol_bounds_with_outputs={:.6}",
                start.elapsed().as_secs_f64()
            );
        }
    } else {
        let start = trace_row64.then(Instant::now);
        set_topol_bounds(
            mol,
            mmat,
            true,
            false,
            params.use_macrocycle14config,
            params.force_trans_amides,
            true,
            true,
        )?;
        if let Some(start) = start {
            println!(
                "row64_setup_bounds set_topol_bounds={:.6}",
                start.elapsed().as_secs_f64()
            );
        }
    }
    trace_bounds_stage("after_initial_set_topol", mmat);
    let mut tol = 0.0;
    if let Some(coord_map) = coord_map {
        embedder_adjust_bounds_mat_from_coord_map(mmat, n_atoms, coord_map)?;
        tol = 0.05;
        trace_bounds_stage("after_coord_map", mmat);
    }
    let smooth_start = trace_row64.then(Instant::now);
    if !triangle_smooth_bounds_shared(mmat, tol) {
        trace_bounds_stage("after_triangle_smooth_first_fail", mmat);
        if let Some(start) = smooth_start {
            println!(
                "row64_setup_bounds triangle_smooth_first={:.6} ok=0",
                start.elapsed().as_secs_f64()
            );
        }
        init_bounds_mat_shared(mmat, DEFAULT_LOWER, DEFAULT_UPPER);
        trace_bounds_stage("after_retry_init", mmat);
        let retry_start = trace_row64.then(Instant::now);
        set_topol_bounds(
            mol,
            mmat,
            false,
            true,
            params.use_macrocycle14config,
            params.force_trans_amides,
            true,
            true,
        )?;
        if let Some(start) = retry_start {
            println!(
                "row64_setup_bounds retry_set_topol_bounds={:.6}",
                start.elapsed().as_secs_f64()
            );
        }
        trace_bounds_stage("after_retry_set_topol", mmat);
        if let Some(coord_map) = coord_map {
            embedder_adjust_bounds_mat_from_coord_map(mmat, n_atoms, coord_map)?;
            trace_bounds_stage("after_retry_coord_map", mmat);
        }
        let second_smooth_start = trace_row64.then(Instant::now);
        if !triangle_smooth_bounds_shared(mmat, tol) {
            trace_bounds_stage("after_triangle_smooth_second_fail", mmat);
            if let Some(start) = second_smooth_start {
                println!(
                    "row64_setup_bounds triangle_smooth_second={:.6} ok=0",
                    start.elapsed().as_secs_f64()
                );
            }
            if params.ignore_smoothing_failures {
                init_bounds_mat_shared(mmat, DEFAULT_LOWER, DEFAULT_UPPER);
                trace_bounds_stage("after_fallback_init", mmat);
                let fallback_start = trace_row64.then(Instant::now);
                set_topol_bounds(
                    mol,
                    mmat,
                    false,
                    true,
                    params.use_macrocycle14config,
                    params.force_trans_amides,
                    true,
                    true,
                )?;
                if let Some(start) = fallback_start {
                    println!(
                        "row64_setup_bounds fallback_set_topol_bounds={:.6}",
                        start.elapsed().as_secs_f64()
                    );
                }
                trace_bounds_stage("after_fallback_set_topol", mmat);
                if let Some(coord_map) = coord_map {
                    embedder_adjust_bounds_mat_from_coord_map(mmat, n_atoms, coord_map)?;
                    trace_bounds_stage("after_fallback_coord_map", mmat);
                }
            } else {
                return Ok(false);
            }
        } else if let Some(start) = second_smooth_start {
            trace_bounds_stage("after_triangle_smooth_second_ok", mmat);
            println!(
                "row64_setup_bounds triangle_smooth_second={:.6} ok=1",
                start.elapsed().as_secs_f64()
            );
        }
    } else if let Some(start) = smooth_start {
        trace_bounds_stage("after_triangle_smooth_first_ok", mmat);
        println!(
            "row64_setup_bounds triangle_smooth_first={:.6} ok=1",
            start.elapsed().as_secs_f64()
        );
    }
    Ok(true)
}

fn embedder_fill_atom_positions(
    pts: &mut [ForceFieldVec3],
    conf: &Conformer3D,
    _mol: &Molecule,
    match_atoms: &[usize],
) {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::_fillAtomPositions (Embedder.cpp:1309-1315)
    // RDKit✔️✔️: void _fillAtomPositions(RDGeom::Point3DConstPtrVect &pts, const Conformer &conf,
    // RDKit✔️✔️:                         const ROMol &, const std::vector<unsigned int> &match) {
    // RDKit✔️✔️:   PRECONDITION(pts.size() == match.size(), "bad pts size");
    // RDKit✔️✔️:   for (unsigned int i = 0; i < match.size(); i++) {
    // RDKit✔️✔️:     pts[i] = &conf.getAtomPos(match[i]);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::_fillAtomPositions
    assert_eq!(pts.len(), match_atoms.len(), "bad pts size");
    for (i, &atom_idx) in match_atoms.iter().enumerate() {
        let coord = conf.coordinates()[atom_idx];
        pts[i] = ForceFieldVec3::new(coord[0], coord[1], coord[2]);
    }
}

fn alignments_weighted_sum_of_points(points: &[ForceFieldVec3]) -> ForceFieldVec3 {
    // BEGIN RDKIT CPP FUNCTION RDNumeric::Alignments::_weightedSumOfPoints unweighted path (AlignPoints.cpp:19-34)
    // RDKit✔️✔️: RDGeom::Point3D _weightedSumOfPoints(const RDGeom::Point3DConstPtrVect &points,
    // RDKit✔️✔️:                                      const DoubleVector *weights) {
    // RDKit✔️✔️:   PRECONDITION(!weights || points.size() == weights->size(), "");
    // RDKit✔️✔️:   RDGeom::Point3DConstPtrVect_CI pti;
    // RDKit✔️✔️:   RDGeom::Point3D tmpPt, res;
    // RDKit✔️✔️:   const double *wData = weights ? weights->getData() : nullptr;
    // RDKit✔️✔️:   unsigned int i = 0;
    // RDKit✔️✔️:   for (pti = points.begin(); pti != points.end(); pti++) {
    // RDKit✔️✔️:     tmpPt = (*(*pti));
    // RDKit✔️✔️:     if (weights) {
    // RDKit✔️✔️:       tmpPt *= wData[i];
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     res += tmpPt;
    // RDKit✔️✔️:     i++;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RDNumeric::Alignments::_weightedSumOfPoints
    points
        .iter()
        .copied()
        .fold(ForceFieldVec3::default(), |acc, point| acc + point)
}

fn alignments_weighted_sum_of_len_sq(points: &[ForceFieldVec3]) -> f64 {
    // BEGIN RDKIT CPP FUNCTION RDNumeric::Alignments::_weightedSumOfLenSq unweighted path (AlignPoints.cpp:36-49)
    // RDKit✔️✔️: double _weightedSumOfLenSq(const RDGeom::Point3DConstPtrVect &points,
    // RDKit✔️✔️:                            const DoubleVector *weights) {
    // RDKit✔️✔️:   PRECONDITION(!weights || (points.size() == weights->size()), "");
    // RDKit✔️✔️:   double res = 0.0;
    // RDKit✔️✔️:   const double *wData = weights ? weights->getData() : nullptr;
    // RDKit✔️✔️:   unsigned int i = 0;
    // RDKit✔️✔️:   for (const auto &pti : points) {
    // RDKit✔️✔️:     auto l = pti->lengthSq();
    // RDKit✔️✔️:     if (weights) {
    // RDKit✔️✔️:       l *= wData[i];
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     res += l;
    // RDKit✔️✔️:     i++;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RDNumeric::Alignments::_weightedSumOfLenSq
    points.iter().map(|point| point.length_sq()).sum()
}

fn alignments_compute_covariance_mat(
    ref_points: &[ForceFieldVec3],
    probe_points: &[ForceFieldVec3],
) -> [[f64; 3]; 3] {
    // BEGIN RDKIT CPP FUNCTION RDNumeric::Alignments::_computeCovarianceMat unweighted path (AlignPoints.cpp:61-88)
    // RDKit✔️✔️: void _computeCovarianceMat(const RDGeom::Point3DConstPtrVect &refPoints,
    // RDKit✔️✔️:                            const RDGeom::Point3DConstPtrVect &probePoints,
    // RDKit✔️✔️:                            const DoubleVector *weights, double covMat[3][3]) {
    // RDKit✔️✔️:   memset(static_cast<void *>(covMat), 0, 9 * sizeof(double));
    // RDKit✔️✔️:   unsigned int npt = refPoints.size();
    // RDKit✔️✔️:   CHECK_INVARIANT(npt == probePoints.size(), "Number of points mismatch");
    // RDKit✔️✔️:   CHECK_INVARIANT(!weights || (npt == weights->size()),
    // RDKit✔️✔️:                   "Number of points and number of weights do not match");
    // RDKit✔️✔️:   const double *wData = weights ? weights->getData() : nullptr;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   const RDGeom::Point3D *rpt, *ppt;
    // RDKit✔️✔️:   for (unsigned int i = 0; i < npt; i++) {
    // RDKit✔️✔️:     rpt = refPoints[i];
    // RDKit✔️✔️:     ppt = probePoints[i];
    // RDKit✔️✔️:     double w = weights ? wData[i] : 1.0;
    assert_eq!(
        ref_points.len(),
        probe_points.len(),
        "Number of points mismatch"
    );
    let mut cov_mat = [[0.0; 3]; 3];
    for (rpt, ppt) in ref_points.iter().zip(probe_points) {
        let w = 1.0;
        // RDKit✔️✔️:     covMat[0][0] += w * (ppt->x) * (rpt->x);
        // RDKit✔️✔️:     covMat[0][1] += w * (ppt->x) * (rpt->y);
        // RDKit✔️✔️:     covMat[0][2] += w * (ppt->x) * (rpt->z);
        cov_mat[0][0] += w * ppt.x * rpt.x;
        cov_mat[0][1] += w * ppt.x * rpt.y;
        cov_mat[0][2] += w * ppt.x * rpt.z;
        // RDKit✔️✔️:     covMat[1][0] += w * (ppt->y) * (rpt->x);
        // RDKit✔️✔️:     covMat[1][1] += w * (ppt->y) * (rpt->y);
        // RDKit✔️✔️:     covMat[1][2] += w * (ppt->y) * (rpt->z);
        cov_mat[1][0] += w * ppt.y * rpt.x;
        cov_mat[1][1] += w * ppt.y * rpt.y;
        cov_mat[1][2] += w * ppt.y * rpt.z;
        // RDKit✔️✔️:     covMat[2][0] += w * (ppt->z) * (rpt->x);
        // RDKit✔️✔️:     covMat[2][1] += w * (ppt->z) * (rpt->y);
        // RDKit✔️✔️:     covMat[2][2] += w * (ppt->z) * (rpt->z);
        cov_mat[2][0] += w * ppt.z * rpt.x;
        cov_mat[2][1] += w * ppt.z * rpt.y;
        cov_mat[2][2] += w * ppt.z * rpt.z;
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RDNumeric::Alignments::_computeCovarianceMat
    cov_mat
}

fn alignments_convert_cov_mat_to_quad(
    cov_mat: [[f64; 3]; 3],
    rpt_sum: ForceFieldVec3,
    ppt_sum: ForceFieldVec3,
    wts_sum: f64,
) -> [[f64; 4]; 4] {
    // BEGIN RDKIT CPP FUNCTION RDNumeric::Alignments::_covertCovMatToQuad (AlignPoints.cpp:90-129)
    // RDKit✔️✔️: void _covertCovMatToQuad(const double covMat[3][3],
    // RDKit✔️✔️:                          const RDGeom::Point3D &rptSum,
    // RDKit✔️✔️:                          const RDGeom::Point3D &pptSum, double wtsSum,
    // RDKit✔️✔️:                          double quad[4][4]) {
    let temp = ppt_sum.x / wts_sum;
    let px_rx = cov_mat[0][0] - temp * rpt_sum.x;
    let px_ry = cov_mat[0][1] - temp * rpt_sum.y;
    let px_rz = cov_mat[0][2] - temp * rpt_sum.z;

    let temp = ppt_sum.y / wts_sum;
    let py_rx = cov_mat[1][0] - temp * rpt_sum.x;
    let py_ry = cov_mat[1][1] - temp * rpt_sum.y;
    let py_rz = cov_mat[1][2] - temp * rpt_sum.z;

    let temp = ppt_sum.z / wts_sum;
    let pz_rx = cov_mat[2][0] - temp * rpt_sum.x;
    let pz_ry = cov_mat[2][1] - temp * rpt_sum.y;
    let pz_rz = cov_mat[2][2] - temp * rpt_sum.z;

    let mut quad = [[0.0; 4]; 4];
    // RDKit✔️✔️:   quad[0][0] = -2.0 * (PxRx + PyRy + PzRz);
    quad[0][0] = -2.0 * (px_rx + py_ry + pz_rz);
    quad[1][1] = -2.0 * (px_rx - py_ry - pz_rz);
    quad[2][2] = -2.0 * (py_ry - pz_rz - px_rx);
    quad[3][3] = -2.0 * (pz_rz - px_rx - py_ry);
    quad[0][1] = 2.0 * (py_rz - pz_ry);
    quad[1][0] = quad[0][1];
    quad[0][2] = 2.0 * (pz_rx - px_rz);
    quad[2][0] = quad[0][2];
    quad[0][3] = 2.0 * (px_ry - py_rx);
    quad[3][0] = quad[0][3];
    quad[1][2] = -2.0 * (px_ry + py_rx);
    quad[2][1] = quad[1][2];
    quad[1][3] = -2.0 * (pz_rx + px_rz);
    quad[3][1] = quad[1][3];
    quad[2][3] = -2.0 * (py_rz + pz_ry);
    quad[3][2] = quad[2][3];
    // END RDKIT CPP FUNCTION RDNumeric::Alignments::_covertCovMatToQuad
    quad
}

fn alignments_jacobi(
    mut quad: [[f64; 4]; 4],
    eigen_vals: &mut [f64; 4],
    eigen_vecs: &mut [[f64; 4]; 4],
    max_iter: u32,
) -> u32 {
    // BEGIN RDKIT CPP FUNCTION RDNumeric::Alignments::jacobi (AlignPoints.cpp:148-270)
    // RDKit✔️✔️: unsigned int jacobi(double quad[4][4], double eigenVals[4],
    // RDKit✔️✔️:                     double eigenVecs[4][4], unsigned int maxIter) {
    const TOLERANCE: f64 = 1.0e-6;
    for j in 0..=3 {
        for row in eigen_vecs.iter_mut().take(4) {
            row[j] = 0.0;
        }
        eigen_vecs[j][j] = 1.0;
        eigen_vals[j] = quad[j][j];
    }
    let mut l = 0;
    'outer: while l < max_iter {
        let mut diag_norm = 0.0;
        let mut off_diag_norm = 0.0;
        for j in 0..=3 {
            diag_norm += eigen_vals[j].abs();
            for row in quad.iter().take(j) {
                off_diag_norm += row[j].abs();
            }
        }
        if diag_norm.abs() > 1.0e-16 && (off_diag_norm / diag_norm) <= TOLERANCE {
            break 'outer;
        }
        for j in 1..=3 {
            for i in 0..=j - 1 {
                let b = quad[i][j];
                if b.abs() > 0.0 {
                    let dma = eigen_vals[j] - eigen_vals[i];
                    let t = if dma.abs() + b.abs() <= dma.abs() {
                        b / dma
                    } else {
                        let q = 0.5 * dma / b;
                        let mut t = 1.0 / (q.abs() + (1.0 + q * q).sqrt());
                        if q < 0.0 {
                            t = -t;
                        }
                        t
                    };
                    let c = 1.0 / (t * t + 1.0).sqrt();
                    let s = t * c;
                    quad[i][j] = 0.0;
                    for k in 0..i {
                        let atemp = c * quad[k][i] - s * quad[k][j];
                        quad[k][j] = s * quad[k][i] + c * quad[k][j];
                        quad[k][i] = atemp;
                    }
                    for k in i + 1..=j - 1 {
                        let atemp = c * quad[i][k] - s * quad[k][j];
                        quad[k][j] = s * quad[i][k] + c * quad[k][j];
                        quad[i][k] = atemp;
                    }
                    for k in j + 1..=3 {
                        let atemp = c * quad[i][k] - s * quad[j][k];
                        quad[j][k] = s * quad[i][k] + c * quad[j][k];
                        quad[i][k] = atemp;
                    }
                    for row in eigen_vecs.iter_mut().take(4) {
                        let vtemp = c * row[i] - s * row[j];
                        row[j] = s * row[i] + c * row[j];
                        row[i] = vtemp;
                    }
                    let dtemp = c * c * eigen_vals[i] + s * s * eigen_vals[j] - 2.0 * c * s * b;
                    eigen_vals[j] = s * s * eigen_vals[i] + c * c * eigen_vals[j] + 2.0 * c * s * b;
                    eigen_vals[i] = dtemp;
                }
            }
        }
        l += 1;
    }
    for j in 0..=2 {
        let mut k = j;
        let mut dtemp = eigen_vals[k];
        for (i, &val) in eigen_vals.iter().enumerate().take(4).skip(j + 1) {
            if val < dtemp {
                k = i;
                dtemp = val;
            }
        }
        if k > j {
            eigen_vals[k] = eigen_vals[j];
            eigen_vals[j] = dtemp;
            for row in eigen_vecs.iter_mut().take(4) {
                let dtemp = row[k];
                row[k] = row[j];
                row[j] = dtemp;
            }
        }
    }
    // RDKit✔️✔️:   return l + 1;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RDNumeric::Alignments::jacobi
    l + 1
}

fn alignments_align_points_ssr(
    ref_points: &[ForceFieldVec3],
    probe_points: &[ForceFieldVec3],
) -> f64 {
    // BEGIN RDKIT CPP FUNCTION RDNumeric::Alignments::AlignPoints unweighted no-reflect SSR path (AlignPoints.cpp:272-329)
    // RDKit✔️✔️: double AlignPoints(const RDGeom::Point3DConstPtrVect &refPoints,
    // RDKit✔️✔️:                    const RDGeom::Point3DConstPtrVect &probePoints,
    // RDKit✔️✔️:                    RDGeom::Transform3D &trans, const DoubleVector *weights,
    // RDKit✔️✔️:                    bool reflect, unsigned int maxIterations) {
    // RDKit✔️✔️:   unsigned int npt = refPoints.size();
    // RDKit✔️✔️:   PRECONDITION(npt == probePoints.size(), "Mismatch in number of points");
    assert_eq!(
        ref_points.len(),
        probe_points.len(),
        "Mismatch in number of points"
    );
    let npt = ref_points.len();
    let wts_sum = npt as f64;
    let rpt_sum = alignments_weighted_sum_of_points(ref_points);
    let ppt_sum = alignments_weighted_sum_of_points(probe_points);
    let rpt_sum_len_sq = alignments_weighted_sum_of_len_sq(ref_points);
    let ppt_sum_len_sq = alignments_weighted_sum_of_len_sq(probe_points);
    let cov_mat = alignments_compute_covariance_mat(ref_points, probe_points);
    let quad = alignments_convert_cov_mat_to_quad(cov_mat, rpt_sum, ppt_sum, wts_sum);
    let mut eigen_vecs = [[0.0; 4]; 4];
    let mut eigen_vals = [0.0; 4];
    alignments_jacobi(quad, &mut eigen_vals, &mut eigen_vecs, 50);
    // RDKit✔️✔️:   double ssr = eigenVals[0] - (pptSum.lengthSq() + rptSum.lengthSq()) / wtsSum +
    // RDKit✔️✔️:                rptSumLenSq + pptSumLenSq;
    let mut ssr = eigen_vals[0] - (ppt_sum.length_sq() + rpt_sum.length_sq()) / wts_sum
        + rpt_sum_len_sq
        + ppt_sum_len_sq;
    if ssr < 0.0 && ssr.abs() < 1.0e-6 {
        ssr = 0.0;
    }
    // RDKit✔️✔️:   return ssr;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RDNumeric::Alignments::AlignPoints
    ssr
}

#[cfg(test)]
pub(crate) fn aligned_rmsd_for_test(ref_points: &[[f64; 3]], probe_points: &[[f64; 3]]) -> f64 {
    assert_eq!(
        ref_points.len(),
        probe_points.len(),
        "Mismatch in number of points"
    );
    assert!(
        !ref_points.is_empty(),
        "alignment requires at least one point"
    );
    let ref_points: Vec<_> = ref_points
        .iter()
        .map(|point| ForceFieldVec3::new(point[0], point[1], point[2]))
        .collect();
    let probe_points: Vec<_> = probe_points
        .iter()
        .map(|point| ForceFieldVec3::new(point[0], point[1], point[2]))
        .collect();
    (alignments_align_points_ssr(&ref_points, &probe_points) / ref_points.len() as f64).sqrt()
}

fn embedder_is_conf_far_from_rest(
    mol: &Molecule,
    conf: &Conformer3D,
    threshold: f64,
    self_matches: &[Vec<usize>],
) -> bool {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::_isConfFarFromRest (Embedder.cpp:1317-1344)
    // RDKit✔️✔️: bool _isConfFarFromRest(
    // RDKit✔️✔️:     const ROMol &mol, const Conformer &conf, double threshold,
    // RDKit✔️✔️:     const std::vector<std::vector<unsigned int>> &selfMatches) {
    // RDKit✔️✔️:   // NOTE: it is tempting to use some triangle inequality to prune
    // RDKit✔️✔️:   // conformations here but some basic testing has shown very
    // RDKit✔️✔️:   // little advantage and given that the time for pruning fades in
    // RDKit✔️✔️:   // comparison to embedding - we will use a simple for loop below
    // RDKit✔️✔️:   // over all conformation until we find a match
    // RDKit✔️✔️:   RDGeom::Point3DConstPtrVect refPoints(selfMatches[0].size());
    // RDKit✔️✔️:   RDGeom::Point3DConstPtrVect prbPoints(selfMatches[0].size());
    // RDKit✔️✔️:   _fillAtomPositions(refPoints, conf, mol, selfMatches[0]);
    assert!(!self_matches.is_empty(), "expected at least one self match");
    let mut ref_points = vec![ForceFieldVec3::default(); self_matches[0].len()];
    let mut prb_points = vec![ForceFieldVec3::default(); self_matches[0].len()];
    embedder_fill_atom_positions(&mut ref_points, conf, mol, &self_matches[0]);
    // RDKit✔️✔️:   double ssrThres = selfMatches[0].size() * threshold * threshold;
    let ssr_thres = self_matches[0].len() as f64 * threshold * threshold;
    // RDKit✔️✔️:   for (const auto &match : selfMatches) {
    // RDKit✔️✔️:     for (auto confi = mol.beginConformers(); confi != mol.endConformers();
    // RDKit✔️✔️:          ++confi) {
    // RDKit✔️✔️:       _fillAtomPositions(prbPoints, *(*confi), mol, match);
    // RDKit✔️✔️:       RDGeom::Transform3D trans;
    // RDKit✔️✔️:       auto ssr =
    // RDKit✔️✔️:           RDNumeric::Alignments::AlignPoints(refPoints, prbPoints, trans);
    // RDKit✔️✔️:       if (ssr < ssrThres) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    for match_atoms in self_matches {
        for existing in mol.conformers_3d() {
            embedder_fill_atom_positions(&mut prb_points, existing, mol, match_atoms);
            let ssr = alignments_align_points_ssr(&ref_points, &prb_points);
            if ssr < ssr_thres {
                return false;
            }
        }
    }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::_isConfFarFromRest
    true
}

fn molecule_without_hs_for_pruning(
    mol: &Molecule,
) -> Result<(Molecule, Vec<usize>), DgBoundsError> {
    let hydrogens: Vec<_> = mol
        .atoms()
        .iter()
        .filter_map(|atom| (atom.atomic_number() == 1).then_some(atom.id()))
        .collect();
    let mut builder = mol.to_builder();
    let old_to_new = builder.remove_atoms_for_construction(&hydrogens);
    let stripped = builder.build()?;
    let new_to_old = old_to_new
        .iter()
        .enumerate()
        .filter_map(|(old_idx, new_idx)| new_idx.map(|idx| (idx.index(), old_idx)))
        .fold(
            vec![0; stripped.num_atoms()],
            |mut acc, (new_idx, old_idx)| {
                acc[new_idx] = old_idx;
                acc
            },
        );
    Ok((stripped, new_to_old))
}

fn embedder_get_mol_self_matches(
    mol: &Molecule,
    params: &EmbedParameters,
) -> Result<Vec<Vec<usize>>, DgBoundsError> {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::detail::getMolSelfMatches (Embedder.cpp:1449-1494)
    // RDKit✔️✔️: std::vector<std::vector<unsigned int>> getMolSelfMatches(
    // RDKit✔️✔️:     const ROMol &mol, const EmbedParameters &params) {
    // RDKit✔️✔️:   std::vector<std::vector<unsigned int>> res;
    let mut res = Vec::new();
    // RDKit✔️✔️:   if (params.pruneRmsThresh && params.useSymmetryForPruning) {
    if params.prune_rms_thresh != 0.0 && params.use_symmetry_for_pruning {
        // RDKit✔️✔️:     RWMol tmol(mol);
        // RDKit✔️✔️:     MolOps::RemoveHsParameters ps;
        // RDKit✔️✔️:     bool sanitize = false;
        // RDKit✔️✔️:     MolOps::removeHs(tmol, ps, sanitize);
        let (tmol, stripped_to_original) = molecule_without_hs_for_pruning(mol)?;

        // RDKit✔️✔️:     std::unique_ptr<RWMol> prbMolSymm;
        // RDKit✔️✔️:     if (params.symmetrizeConjugatedTerminalGroupsForPruning) {
        // RDKit✔️✔️:       prbMolSymm.reset(new RWMol(tmol));
        // RDKit✔️✔️:       MolAlign::details::symmetrizeTerminalAtoms(*prbMolSymm);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     const auto &prbMolForMatch = prbMolSymm ? *prbMolSymm : tmol;
        let prb_mol_symm = params
            .symmetrize_conjugated_terminal_groups_for_pruning
            .then(|| symmetrize_terminal_atoms_for_pruning(&tmol))
            .transpose()?;
        let prb_mol_for_match = prb_mol_symm.as_ref().unwrap_or(&tmol);

        // RDKit✔️✔️:     SubstructMatchParameters sssps;
        // RDKit✔️✔️:     sssps.maxMatches = 1;
        // RDKit✔️✔️:     // provides the atom indices in the molecule corresponding
        // RDKit✔️✔️:     // to the indices in the H-stripped version
        // RDKit✔️✔️:     auto strippedMatch = SubstructMatch(mol, prbMolForMatch, sssps);
        // RDKit✔️✔️:     CHECK_INVARIANT(strippedMatch.size() == 1, "expected match not found");
        let stripped_match = stripped_to_original;

        // RDKit✔️✔️:     sssps.maxMatches = 1000;
        // RDKit✔️✔️:     sssps.uniquify = false;
        // RDKit✔️✔️:     auto heavyAtomMatches = SubstructMatch(tmol, prbMolForMatch, sssps);
        let sssps = SubstructMatchParams {
            max_matches: 1000,
            uniquify: false,
            use_chirality: false,
            specified_stereo_query_matches_unspecified: false,
        };
        let heavy_atom_matches =
            get_substruct_matches_with_params(&tmol, prb_mol_for_match, &sssps);
        // RDKit✔️✔️:     for (const auto &match : heavyAtomMatches) {
        // RDKit✔️✔️:       res.emplace_back(0);
        // RDKit✔️✔️:       res.back().reserve(match.size());
        // RDKit✔️✔️:       for (auto midx : match) {
        // RDKit✔️✔️:         res.back().push_back(strippedMatch[0][midx.second].second);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        for match_result in heavy_atom_matches {
            let mut mapped = Vec::with_capacity(match_result.atom_mapping.len());
            for &midx_second in &match_result.atom_mapping {
                mapped.push(stripped_match[midx_second]);
            }
            res.push(mapped);
        }
    // RDKit✔️✔️:   } else if (params.onlyHeavyAtomsForRMS) {
    } else if params.only_heavy_atoms_for_rms {
        // RDKit✔️✔️:     res.emplace_back(0);
        // RDKit✔️✔️:     for (const auto &at : mol.atoms()) {
        // RDKit✔️✔️:       if (at->getAtomicNum() != 1) {
        // RDKit✔️✔️:         res.back().push_back(at->getIdx());
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        let atoms = mol
            .atoms()
            .iter()
            .filter_map(|atom| (atom.atomic_number() != 1).then_some(atom.id().index()))
            .collect();
        res.push(atoms);
    // RDKit✔️✔️:   } else {
    } else {
        // RDKit✔️✔️:     res.emplace_back(0);
        // RDKit✔️✔️:     res.back().reserve(mol.getNumAtoms());
        // RDKit✔️✔️:     for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
        // RDKit✔️✔️:       res.back().push_back(i);
        // RDKit✔️✔️:     }
        res.push((0..mol.num_atoms()).collect());
    }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::detail::getMolSelfMatches
    Ok(res)
}

fn symmetrize_terminal_atoms_for_pruning(mol: &Molecule) -> Result<Molecule, DgBoundsError> {
    // BEGIN RDKIT CPP FUNCTION MolAlign::details::symmetrizeTerminalAtoms (AlignMolecules.cpp:29-54)
    // RDKit✔️✔️: void symmetrizeTerminalAtoms(RWMol &mol) {
    // RDKit✔️✔️:   // clang-format off
    // RDKit✔️✔️:   static const std::string qsmarts =
    // RDKit✔️✔️:       "[{atomPattern};$([{atomPattern}]-[*]=[{atomPattern}]),$([{atomPattern}]=[*]-[{atomPattern}])]~[*]";
    // RDKit✔️✔️:   static std::map<std::string, std::string> replacements = {
    // RDKit✔️✔️:       {"{atomPattern}", "O,N;D1"}};
    // RDKit✔️✔️:   // clang-format on
    // RDKit✔️✔️:   static SmartsParserParams ps;
    // RDKit✔️✔️:   ps.replacements = &replacements;
    // RDKit✔️✔️:   static const std::unique_ptr<RWMol> qry{SmartsToMol(qsmarts, ps)};
    // RDKit✔️✔️:   CHECK_INVARIANT(qry, "bad query pattern");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto matches = SubstructMatch(mol, *qry);
    // RDKit✔️✔️:   if (matches.empty()) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   QueryBond qb;
    // RDKit✔️✔️:   qb.setQuery(makeSingleOrDoubleBondQuery());
    // RDKit✔️✔️:   for (const auto &match : matches) {
    // RDKit✔️✔️:     mol.getAtomWithIdx(match[0].second)->setFormalCharge(0);
    // RDKit✔️✔️:     auto obond = mol.getBondBetweenAtoms(match[0].second, match[1].second);
    // RDKit✔️✔️:     CHECK_INVARIANT(obond, "could not find expected bond");
    // RDKit✔️✔️:     mol.replaceBond(obond->getIdx(), &qb);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION MolAlign::details::symmetrizeTerminalAtoms
    let mut symmetrized = mol.clone();
    let mut matched_pairs = Vec::new();
    for atom_idx in 0..mol.num_atoms() {
        if !matches_rdkit_terminal_atom_pattern(mol, atom_idx) {
            continue;
        }
        let neighbors = neighbors_for_atom(mol, atom_idx);
        debug_assert_eq!(neighbors.len(), 1);
        let neighbor_idx = neighbors[0];
        if matches_rdkit_terminal_group_symmetry_query(mol, atom_idx, neighbor_idx) {
            matched_pairs.push((atom_idx, neighbor_idx));
        }
    }
    if matched_pairs.is_empty() {
        return Ok(symmetrized);
    }

    {
        let topology = symmetrized.topology_block_mut();
        for (terminal_atom_idx, neighbor_idx) in matched_pairs {
            topology.atoms[terminal_atom_idx].set_formal_charge(0);
            let Some(bond) = topology.bonds.iter_mut().find(|bond| {
                (bond.begin().index() == terminal_atom_idx && bond.end().index() == neighbor_idx)
                    || (bond.begin().index() == neighbor_idx
                        && bond.end().index() == terminal_atom_idx)
            }) else {
                return Err(DgBoundsError::GenerationFailed(
                    "MolAlign::details::symmetrizeTerminalAtoms could not find expected bond"
                        .to_string(),
                ));
            };
            bond.set_query(Some(QueryNode::predicate(BondQueryPredicate::OrderIn(
                vec![BondOrder::Single, BondOrder::Double],
            ))));
        }
    }
    symmetrized
        .derived_cache_mut()
        .invalidate(DerivedState::VALENCE | DerivedState::AROMATICITY | DerivedState::STEREO);
    Ok(symmetrized)
}

fn matches_rdkit_terminal_group_symmetry_query(
    mol: &Molecule,
    terminal_atom_idx: usize,
    neighbor_idx: usize,
) -> bool {
    let Some(terminal_bond) = bond_between(mol, terminal_atom_idx, neighbor_idx) else {
        return false;
    };
    let opposite_order = match terminal_bond.order() {
        BondOrder::Single => BondOrder::Double,
        BondOrder::Double => BondOrder::Single,
        _ => return false,
    };
    for remote_terminal_idx in neighbors_for_atom(mol, neighbor_idx) {
        if remote_terminal_idx == terminal_atom_idx
            || !matches_rdkit_terminal_atom_pattern(mol, remote_terminal_idx)
        {
            continue;
        }
        if let Some(remote_bond) = bond_between(mol, neighbor_idx, remote_terminal_idx)
            && remote_bond.order() == opposite_order
        {
            return true;
        }
    }
    false
}

fn matches_rdkit_terminal_atom_pattern(mol: &Molecule, atom_idx: usize) -> bool {
    let atom = &mol.atoms()[atom_idx];
    matches!(atom.atomic_number(), 7 | 8) && neighbors_for_atom(mol, atom_idx).len() == 1
}

fn embedder_embed_helper(
    thread_id: i32,
    num_threads: i32,
    eargs: &mut EmbedHelperArgs<'_>,
    params: &mut EmbedParameters,
    end_time: Option<Instant>,
) -> Result<(), DgBoundsError> {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::detail::embedHelper_ (Embedder.cpp:1356-1447)
    // RDKit✔️✔️: void embedHelper_(int threadId, int numThreads, EmbedArgs *eargs,
    // RDKit✔️✔️:                   EmbedParameters *params, TimePoint *end_time) {
    // RDKit✔️✔️:   PRECONDITION(eargs, "bogus eargs");
    // RDKit✔️✔️:   PRECONDITION(params, "bogus params");
    // RDKit✔️✔️:   unsigned int nAtoms = eargs->mmat->numRows();
    let n_atoms = eargs.mmat.num_rows();
    // RDKit✔️✔️:   RDGeom::PointPtrVect positions(nAtoms);
    // RDKit✔️✔️:   // we might thrown an exception in a callback
    // RDKit✔️✔️:   // in order to avoid leaking the points we're working with
    // RDKit✔️✔️:   // allocate them with unique_ptrs and then work with the naked
    // RDKit✔️✔️:   // pointers from those
    // RDKit✔️✔️:   std::vector<std::unique_ptr<RDGeom::Point>> positionsStore;
    // RDKit✔️✔️:   positionsStore.reserve(nAtoms);
    let dim = if eargs.four_d { 4 } else { 3 };
    let mut positions = vec![vec![0.0; dim]; n_atoms];
    // RDKit✔️✔️:   for (unsigned int i = 0; i < nAtoms; ++i) {
    // RDKit✔️✔️:     if (eargs->fourD) {
    // RDKit✔️✔️:       positionsStore.emplace_back(new RDGeom::PointND(4));
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       positionsStore.emplace_back(new RDGeom::Point3D());
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     positions[i] = positionsStore[i].get();
    // RDKit✔️✔️:   }

    // RDKit✔️✔️:   for (size_t ci = 0; ci < eargs->confs->size(); ci++) {
    for ci in 0..eargs.confs.len() {
        // RDKit✔️✔️:     if (ControlCHandler::getGotSignal() ||
        // RDKit✔️✔️:         (end_time != nullptr && Clock::now() > *end_time)) {
        // RDKit✔️✔️:       return;
        // RDKit✔️✔️:     }
        if let Some(deadline) = end_time
            && Instant::now() > deadline
        {
            return Ok(());
        }

        // RDKit✔️✔️:     if (rdcast<int>(ci % numThreads) != threadId) {
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        if (ci % num_threads as usize) as i32 != thread_id {
            continue;
        }
        // RDKit✔️✔️:     if (!(*eargs->confsOk)[ci]) {
        // RDKit✔️✔️:       // we call this function for each fragment in a molecule,
        // RDKit✔️✔️:       // if one of the fragments has already failed, there's no
        // RDKit✔️✔️:       // sense in embedding this one
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        if !eargs.confs_ok[ci] {
            continue;
        }

        let new_seed = rdkit_embedder_conformer_seed(
            params.random_seed,
            ci,
            params.enable_sequential_random_seeds,
        );
        let got_coords = embedder_embed_points(
            &mut positions,
            eargs.embed_args(),
            params,
            new_seed,
            end_time,
        )?;

        // RDKit✔️✔️:     // copy the coordinates into the correct conformer
        // RDKit✔️✔️:     if (gotCoords) {
        if got_coords {
            let conf = &mut eargs.confs[ci];
            // RDKit✔️✔️:       unsigned int fragAtomIdx = 0;
            let mut frag_atom_idx = 0;
            // RDKit✔️✔️:       for (unsigned int i = 0; i < conf->getNumAtoms(); ++i) {
            for i in 0..conf.coordinates().len() {
                // RDKit✔️✔️:         if (!eargs->fragMapping ||
                // RDKit✔️✔️:             (*eargs->fragMapping)[i] == static_cast<int>(eargs->fragIdx)) {
                if eargs
                    .frag_mapping
                    .is_none_or(|mapping| mapping[i] == eargs.frag_idx)
                {
                    // RDKit✔️✔️:           conf->setAtomPos(i, RDGeom::Point3D((*positions[fragAtomIdx])[0],
                    // RDKit✔️✔️:                                               (*positions[fragAtomIdx])[1],
                    // RDKit✔️✔️:                                               (*positions[fragAtomIdx])[2]));
                    conf.coordinates_mut()[i] = [
                        positions[frag_atom_idx][0],
                        positions[frag_atom_idx][1],
                        positions[frag_atom_idx][2],
                    ];
                    // RDKit✔️✔️:           ++fragAtomIdx;
                    frag_atom_idx += 1;
                }
            }
        // RDKit✔️✔️:     } else {
        } else {
            // RDKit✔️✔️:       (*eargs->confsOk)[ci] = 0;
            eargs.confs_ok[ci] = false;
        }
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::detail::embedHelper_
    Ok(())
}

// ──────────────────────────────────────────────
// Error type
// ──────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum DgBoundsError {
    #[error("distance bounds matrix generation failed: {0}")]
    GenerationFailed(String),
    #[error("invalid distance bounds: {0}")]
    InvalidBounds(String),
    #[error("invalid embed parameters JSON: {0}")]
    InvalidEmbedParametersJson(String),
    #[error("molecule has no atoms")]
    EmptyMolecule,
    #[error(
        "Only version 1 and 2 of the experimental torsion-angle preferences (ETversion) supported"
    )]
    UnsupportedEtVersion,
    #[error("coordinate block update failed: {0}")]
    CoordinateUpdateFailed(String),
    #[error(
        "RDKit clock-derived implicit conformer seed is unsupported on wasm32; set EmbedParameters.random_seed to a non-negative explicit seed"
    )]
    WasmImplicitClockSeedUnsupported,
    #[error(transparent)]
    UnsupportedFeature(#[from] crate::UnsupportedFeatureError),
}

impl From<MoleculeBuildError> for DgBoundsError {
    fn from(value: MoleculeBuildError) -> Self {
        Self::CoordinateUpdateFailed(value.to_string())
    }
}

impl From<crate::RingFindingError> for DgBoundsError {
    fn from(value: crate::RingFindingError) -> Self {
        Self::GenerationFailed(value.to_string())
    }
}

// ──────────────────────────────────────────────
// Embedder status
// ──────────────────────────────────────────────

// BEGIN RDKIT CPP ENUM DGeomHelpers::EmbedFailureCauses (Embedder.h:25-39)
// RDKit✔️✔️: enum EmbedFailureCauses {
// RDKit✔️✔️:   INITIAL_COORDS = 0,
// RDKit✔️✔️:   FIRST_MINIMIZATION = 1,
// RDKit✔️✔️:   CHECK_TETRAHEDRAL_CENTERS = 2,
// RDKit✔️✔️:   CHECK_CHIRAL_CENTERS = 3,
// RDKit✔️✔️:   MINIMIZE_FOURTH_DIMENSION = 4,
// RDKit✔️✔️:   ETK_MINIMIZATION = 5,
// RDKit✔️✔️:   FINAL_CHIRAL_BOUNDS = 6,
// RDKit✔️✔️:   FINAL_CENTER_IN_VOLUME = 7,
// RDKit✔️✔️:   LINEAR_DOUBLE_BOND = 8,
// RDKit✔️✔️:   BAD_DOUBLE_BOND_STEREO = 9,
// RDKit✔️✔️:   CHECK_CHIRAL_CENTERS2 = 10,
// RDKit✔️✔️:   EXCEEDED_TIMEOUT = 11,
// RDKit✔️✔️:   END_OF_ENUM = 12,
// RDKit✔️✔️: };
// END RDKIT CPP ENUM DGeomHelpers::EmbedFailureCauses
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u32)]
pub enum EmbedFailureCause {
    InitialCoords = 0,
    FirstMinimization = 1,
    CheckTetrahedralCenters = 2,
    CheckChiralCenters = 3,
    MinimizeFourthDimension = 4,
    EtkMinimization = 5,
    FinalChiralBounds = 6,
    FinalCenterInVolume = 7,
    LinearDoubleBond = 8,
    BadDoubleBondStereo = 9,
    CheckChiralCenters2 = 10,
    ExceededTimeout = 11,
    EndOfEnum = 12,
}

impl EmbedFailureCause {
    pub const ALL: [Self; 13] = [
        Self::InitialCoords,
        Self::FirstMinimization,
        Self::CheckTetrahedralCenters,
        Self::CheckChiralCenters,
        Self::MinimizeFourthDimension,
        Self::EtkMinimization,
        Self::FinalChiralBounds,
        Self::FinalCenterInVolume,
        Self::LinearDoubleBond,
        Self::BadDoubleBondStereo,
        Self::CheckChiralCenters2,
        Self::ExceededTimeout,
        Self::EndOfEnum,
    ];

    #[must_use]
    pub const fn rdkit_ordinal(self) -> u32 {
        self as u32
    }

    #[must_use]
    pub const fn from_rdkit_ordinal(value: u32) -> Option<Self> {
        match value {
            0 => Some(Self::InitialCoords),
            1 => Some(Self::FirstMinimization),
            2 => Some(Self::CheckTetrahedralCenters),
            3 => Some(Self::CheckChiralCenters),
            4 => Some(Self::MinimizeFourthDimension),
            5 => Some(Self::EtkMinimization),
            6 => Some(Self::FinalChiralBounds),
            7 => Some(Self::FinalCenterInVolume),
            8 => Some(Self::LinearDoubleBond),
            9 => Some(Self::BadDoubleBondStereo),
            10 => Some(Self::CheckChiralCenters2),
            11 => Some(Self::ExceededTimeout),
            12 => Some(Self::EndOfEnum),
            _ => None,
        }
    }

    #[must_use]
    pub const fn rdkit_name(self) -> &'static str {
        match self {
            Self::InitialCoords => "INITIAL_COORDS",
            Self::FirstMinimization => "FIRST_MINIMIZATION",
            Self::CheckTetrahedralCenters => "CHECK_TETRAHEDRAL_CENTERS",
            Self::CheckChiralCenters => "CHECK_CHIRAL_CENTERS",
            Self::MinimizeFourthDimension => "MINIMIZE_FOURTH_DIMENSION",
            Self::EtkMinimization => "ETK_MINIMIZATION",
            Self::FinalChiralBounds => "FINAL_CHIRAL_BOUNDS",
            Self::FinalCenterInVolume => "FINAL_CENTER_IN_VOLUME",
            Self::LinearDoubleBond => "LINEAR_DOUBLE_BOND",
            Self::BadDoubleBondStereo => "BAD_DOUBLE_BOND_STEREO",
            Self::CheckChiralCenters2 => "CHECK_CHIRAL_CENTERS2",
            Self::ExceededTimeout => "EXCEEDED_TIMEOUT",
            Self::EndOfEnum => "END_OF_ENUM",
        }
    }
}

// BEGIN RDKIT CPP STRUCT DGeomHelpers::EmbedParameters (Embedder.h:122-191)
// RDKit✔️✔️: struct RDKIT_DISTGEOMHELPERS_EXPORT EmbedParameters {
// RDKit✔️✔️:   unsigned int maxIterations{0};
// RDKit✔️✔️:   int numThreads{1};
// RDKit✔️✔️:   int randomSeed{-1};
// RDKit✔️✔️:   bool clearConfs{true};
// RDKit✔️✔️:   bool useRandomCoords{false};
// RDKit✔️✔️:   double boxSizeMult{2.0};
// RDKit✔️✔️:   bool randNegEig{true};
// RDKit✔️✔️:   unsigned int numZeroFail{1};
// RDKit✔️✔️:   const std::map<int, RDGeom::Point3D> *coordMap{nullptr};
// RDKit✔️✔️:   double optimizerForceTol{1e-3};
// RDKit✔️✔️:   bool ignoreSmoothingFailures{false};
// RDKit✔️✔️:   bool enforceChirality{true};
// RDKit✔️✔️:   bool useExpTorsionAnglePrefs{false};
// RDKit✔️✔️:   bool useBasicKnowledge{false};
// RDKit✔️✔️:   bool verbose{false};
// RDKit✔️✔️:   double basinThresh{5.0};
// RDKit✔️✔️:   double pruneRmsThresh{-1.0};
// RDKit✔️✔️:   bool onlyHeavyAtomsForRMS{true};
// RDKit✔️✔️:   unsigned int ETversion{2};
// RDKit✔️✔️:   boost::shared_ptr<const DistGeom::BoundsMatrix> boundsMat;
// RDKit✔️✔️:   bool embedFragmentsSeparately{true};
// RDKit✔️✔️:   bool useSmallRingTorsions{false};
// RDKit✔️✔️:   bool useMacrocycleTorsions{false};
// RDKit✔️✔️:   bool useMacrocycle14config{false};
// RDKit✔️✔️:   unsigned int timeout{0};
// RDKit✔️✔️:   std::shared_ptr<std::map<std::pair<unsigned int, unsigned int>, double>> CPCI;
// RDKit✔️✔️:   void (*callback)(unsigned int);
// RDKit✔️✔️:   bool forceTransAmides{true};
// RDKit✔️✔️:   bool useSymmetryForPruning{true};
// RDKit✔️✔️:   double boundsMatForceScaling{1.0};
// RDKit✔️✔️:   bool trackFailures{false};
// RDKit✔️✔️:   std::vector<unsigned int> failures;
// RDKit✔️✔️:   bool enableSequentialRandomSeeds{false};
// RDKit✔️✔️:   bool symmetrizeConjugatedTerminalGroupsForPruning{true};
#[derive(Clone)]
pub struct EmbedParameters {
    pub max_iterations: u32,
    pub num_threads: i32,
    pub random_seed: i32,
    pub clear_confs: bool,
    pub use_random_coords: bool,
    pub box_size_mult: f64,
    pub rand_neg_eig: bool,
    pub num_zero_fail: u32,
    pub coord_map: Option<BTreeMap<i32, ForceFieldVec3>>,
    pub optimizer_force_tol: f64,
    pub ignore_smoothing_failures: bool,
    pub enforce_chirality: bool,
    pub use_exp_torsion_angle_prefs: bool,
    pub use_basic_knowledge: bool,
    pub verbose: bool,
    pub basin_thresh: f64,
    pub prune_rms_thresh: f64,
    pub only_heavy_atoms_for_rms: bool,
    pub et_version: u32,
    bounds_mat: Option<Arc<BoundsMatrix>>,
    pub embed_fragments_separately: bool,
    pub use_small_ring_torsions: bool,
    pub use_macrocycle_torsions: bool,
    pub use_macrocycle14config: bool,
    pub timeout: u32,
    pub cpci: Option<BTreeMap<(u32, u32), f64>>,
    pub callback: Option<fn(u32)>,
    pub force_trans_amides: bool,
    pub use_symmetry_for_pruning: bool,
    pub bounds_mat_force_scaling: f64,
    pub track_failures: bool,
    pub failures: Vec<u32>,
    pub enable_sequential_random_seeds: bool,
    pub symmetrize_conjugated_terminal_groups_for_pruning: bool,
}

#[derive(Clone)]
pub struct EmbedMoleculeResult {
    pub molecule: Molecule,
    pub conf_id: i32,
    pub params: EmbedParameters,
}

impl EmbedMoleculeResult {
    #[must_use]
    pub fn ok(&self) -> bool {
        self.conf_id >= 0
    }
}

#[derive(Clone)]
pub struct EmbedMultipleConfsResult {
    pub molecule: Molecule,
    pub conf_ids: Vec<i32>,
    pub requested_num_confs: u32,
    pub params: EmbedParameters,
}

impl EmbedMultipleConfsResult {
    #[must_use]
    pub fn generated_count(&self) -> usize {
        self.conf_ids.len()
    }
}

impl Default for EmbedParameters {
    fn default() -> Self {
        // RDKit✔️✔️:   EmbedParameters() : boundsMat(nullptr), CPCI(nullptr), callback(nullptr) {}
        Self {
            max_iterations: 0,
            num_threads: 1,
            random_seed: -1,
            clear_confs: true,
            use_random_coords: false,
            box_size_mult: 2.0,
            rand_neg_eig: true,
            num_zero_fail: 1,
            coord_map: None,
            optimizer_force_tol: 1e-3,
            ignore_smoothing_failures: false,
            enforce_chirality: true,
            use_exp_torsion_angle_prefs: false,
            use_basic_knowledge: false,
            verbose: false,
            basin_thresh: 5.0,
            prune_rms_thresh: -1.0,
            only_heavy_atoms_for_rms: true,
            et_version: 2,
            bounds_mat: None,
            embed_fragments_separately: true,
            use_small_ring_torsions: false,
            use_macrocycle_torsions: false,
            use_macrocycle14config: false,
            timeout: 0,
            cpci: None,
            callback: None,
            force_trans_amides: true,
            use_symmetry_for_pruning: true,
            bounds_mat_force_scaling: 1.0,
            track_failures: false,
            failures: Vec::new(),
            enable_sequential_random_seeds: false,
            symmetrize_conjugated_terminal_groups_for_pruning: true,
        }
    }
}

impl EmbedParameters {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    pub fn set_bounds_matrix(&mut self, bounds: Vec<Vec<f64>>) -> Result<(), DgBoundsError> {
        let n = bounds.len();
        if bounds.iter().any(|row| row.len() != n) {
            return Err(DgBoundsError::GenerationFailed(
                "bounds matrix must be square".to_string(),
            ));
        }
        for row in &bounds {
            for &value in row {
                if !value.is_finite() || value < 0.0 {
                    return Err(DgBoundsError::GenerationFailed(
                        "bounds matrix values must be finite and non-negative".to_string(),
                    ));
                }
            }
        }
        self.bounds_mat = Some(Arc::new(BoundsMatrix { data: bounds, n }));
        Ok(())
    }

    #[must_use]
    pub fn has_bounds_matrix(&self) -> bool {
        self.bounds_mat.is_some()
    }

    #[allow(clippy::too_many_arguments)]
    fn from_rdkit_constructor(
        max_iterations: u32,
        num_threads: i32,
        random_seed: i32,
        clear_confs: bool,
        use_random_coords: bool,
        box_size_mult: f64,
        rand_neg_eig: bool,
        num_zero_fail: u32,
        coord_map: Option<BTreeMap<i32, ForceFieldVec3>>,
        optimizer_force_tol: f64,
        ignore_smoothing_failures: bool,
        enforce_chirality: bool,
        use_exp_torsion_angle_prefs: bool,
        use_basic_knowledge: bool,
        verbose: bool,
        basin_thresh: f64,
        prune_rms_thresh: f64,
        only_heavy_atoms_for_rms: bool,
        et_version: u32,
        bounds_mat: Option<Arc<BoundsMatrix>>,
        embed_fragments_separately: bool,
        use_small_ring_torsions: bool,
        use_macrocycle_torsions: bool,
        use_macrocycle14config: bool,
        timeout: u32,
        cpci: Option<BTreeMap<(u32, u32), f64>>,
        callback: Option<fn(u32)>,
    ) -> Self {
        // RDKit✔️✔️:       : maxIterations(maxIterations),
        // RDKit✔️✔️:         numThreads(numThreads),
        // RDKit✔️✔️:         randomSeed(randomSeed),
        // RDKit✔️✔️:         clearConfs(clearConfs),
        // RDKit✔️✔️:         useRandomCoords(useRandomCoords),
        // RDKit✔️✔️:         boxSizeMult(boxSizeMult),
        // RDKit✔️✔️:         randNegEig(randNegEig),
        // RDKit✔️✔️:         numZeroFail(numZeroFail),
        // RDKit✔️✔️:         coordMap(coordMap),
        // RDKit✔️✔️:         optimizerForceTol(optimizerForceTol),
        // RDKit✔️✔️:         ignoreSmoothingFailures(ignoreSmoothingFailures),
        // RDKit✔️✔️:         enforceChirality(enforceChirality),
        // RDKit✔️✔️:         useExpTorsionAnglePrefs(useExpTorsionAnglePrefs),
        // RDKit✔️✔️:         useBasicKnowledge(useBasicKnowledge),
        // RDKit✔️✔️:         verbose(verbose),
        // RDKit✔️✔️:         basinThresh(basinThresh),
        // RDKit✔️✔️:         pruneRmsThresh(pruneRmsThresh),
        // RDKit✔️✔️:         onlyHeavyAtomsForRMS(onlyHeavyAtomsForRMS),
        // RDKit✔️✔️:         ETversion(ETversion),
        // RDKit✔️✔️:         boundsMat(boundsMat),
        // RDKit✔️✔️:         embedFragmentsSeparately(embedFragmentsSeparately),
        // RDKit✔️✔️:         useSmallRingTorsions(useSmallRingTorsions),
        // RDKit✔️✔️:         useMacrocycleTorsions(useMacrocycleTorsions),
        // RDKit✔️✔️:         useMacrocycle14config(useMacrocycle14config),
        // RDKit✔️✔️:         timeout(timeout),
        // RDKit✔️✔️:         CPCI(std::move(CPCI)),
        // RDKit✔️✔️:         callback(callback) {}
        let mut params = Self {
            max_iterations,
            num_threads,
            random_seed,
            clear_confs,
            use_random_coords,
            box_size_mult,
            rand_neg_eig,
            num_zero_fail,
            coord_map,
            optimizer_force_tol,
            ignore_smoothing_failures,
            enforce_chirality,
            use_exp_torsion_angle_prefs,
            use_basic_knowledge,
            verbose,
            basin_thresh,
            prune_rms_thresh,
            only_heavy_atoms_for_rms,
            et_version,
            bounds_mat,
            embed_fragments_separately,
            use_small_ring_torsions,
            use_macrocycle_torsions,
            use_macrocycle14config,
            timeout,
            cpci,
            callback,
            ..Self::default()
        };
        params.failures = Vec::new();
        params
    }

    #[must_use]
    pub fn kdg() -> Self {
        // RDKit✔️✔️: const EmbedParameters KDG(0,        // maxIterations
        // RDKit✔️✔️:                           1,        // numThreads
        // RDKit✔️✔️:                           -1,       // randomSeed
        // RDKit✔️✔️:                           true,     // clearConfs
        // RDKit✔️✔️:                           false,    // useRandomCoords
        // RDKit✔️✔️:                           2.0,      // boxSizeMult
        // RDKit✔️✔️:                           true,     // randNegEig
        // RDKit✔️✔️:                           1,        // numZeroFail
        // RDKit✔️✔️:                           nullptr,  // coordMap
        // RDKit✔️✔️:                           1e-3,     // optimizerForceTol
        // RDKit✔️✔️:                           false,    // ignoreSmoothingFailures
        // RDKit✔️✔️:                           true,     // enforceChirality
        // RDKit✔️✔️:                           false,    // useExpTorsionAnglePrefs
        // RDKit✔️✔️:                           true,     // useBasicKnowledge
        // RDKit✔️✔️:                           false,    // verbose
        // RDKit✔️✔️:                           5.0,      // basinThresh
        // RDKit✔️✔️:                           -1.0,     // pruneRmsThresh
        // RDKit✔️✔️:                           true,     // onlyHeavyAtomsForRMS
        // RDKit✔️✔️:                           1,        // ETversion
        // RDKit✔️✔️:                           nullptr,  // boundsMat
        // RDKit✔️✔️:                           true,     // embedFragmentsSeparately
        // RDKit✔️✔️:                           false,    // useSmallRingTorsions
        // RDKit✔️✔️:                           false,    // useMacrocycleTorsions
        // RDKit✔️✔️:                           false,    // useMacrocycle14config
        // RDKit✔️✔️:                           0,        // timeout
        // RDKit✔️✔️:                           nullptr,  // CPCI
        // RDKit✔️✔️:                           nullptr   // callback
        Self::from_rdkit_constructor(
            0, 1, -1, true, false, 2.0, true, 1, None, 1e-3, false, true, false, true, false, 5.0,
            -1.0, true, 1, None, true, false, false, false, 0, None, None,
        )
    }

    #[must_use]
    pub fn etdg() -> Self {
        // RDKit✔️✔️: const EmbedParameters ETDG(0,        // maxIterations
        // RDKit✔️✔️:                            1,        // numThreads
        // RDKit✔️✔️:                            -1,       // randomSeed
        // RDKit✔️✔️:                            true,     // clearConfs
        // RDKit✔️✔️:                            false,    // useRandomCoords
        // RDKit✔️✔️:                            2.0,      // boxSizeMult
        // RDKit✔️✔️:                            true,     // randNegEig
        // RDKit✔️✔️:                            1,        // numZeroFail
        // RDKit✔️✔️:                            nullptr,  // coordMap
        // RDKit✔️✔️:                            1e-3,     // optimizerForceTol
        // RDKit✔️✔️:                            false,    // ignoreSmoothingFailures
        // RDKit✔️✔️:                            false,    // enforceChirality
        // RDKit✔️✔️:                            true,     // useExpTorsionAnglePrefs
        // RDKit✔️✔️:                            false,    // useBasicKnowledge
        // RDKit✔️✔️:                            false,    // verbose
        // RDKit✔️✔️:                            5.0,      // basinThresh
        // RDKit✔️✔️:                            -1.0,     // pruneRmsThresh
        // RDKit✔️✔️:                            true,     // onlyHeavyAtomsForRMS
        // RDKit✔️✔️:                            1,        // ETversion
        // RDKit✔️✔️:                            nullptr,  // boundsMat
        // RDKit✔️✔️:                            true,     // embedFragmentsSeparately
        // RDKit✔️✔️:                            false,    // useSmallRingTorsions
        // RDKit✔️✔️:                            false,    // useMacrocycleTorsions
        // RDKit✔️✔️:                            false,    // useMacrocycle14config
        // RDKit✔️✔️:                            0,        // timeout
        // RDKit✔️✔️:                            nullptr,  // CPCI
        // RDKit✔️✔️:                            nullptr   // callback
        Self::from_rdkit_constructor(
            0, 1, -1, true, false, 2.0, true, 1, None, 1e-3, false, false, true, false, false, 5.0,
            -1.0, true, 1, None, true, false, false, false, 0, None, None,
        )
    }

    #[must_use]
    pub fn etdg_v2() -> Self {
        // RDKit✔️✔️: const EmbedParameters ETDGv2(0,        // maxIterations
        // RDKit✔️✔️:                              1,        // numThreads
        // RDKit✔️✔️:                              -1,       // randomSeed
        // RDKit✔️✔️:                              true,     // clearConfs
        // RDKit✔️✔️:                              false,    // useRandomCoords
        // RDKit✔️✔️:                              2.0,      // boxSizeMult
        // RDKit✔️✔️:                              true,     // randNegEig
        // RDKit✔️✔️:                              1,        // numZeroFail
        // RDKit✔️✔️:                              nullptr,  // coordMap
        // RDKit✔️✔️:                              1e-3,     // optimizerForceTol
        // RDKit✔️✔️:                              false,    // ignoreSmoothingFailures
        // RDKit✔️✔️:                              false,    // enforceChirality
        // RDKit✔️✔️:                              true,     // useExpTorsionAnglePrefs
        // RDKit✔️✔️:                              false,    // useBasicKnowledge
        // RDKit✔️✔️:                              false,    // verbose
        // RDKit✔️✔️:                              5.0,      // basinThresh
        // RDKit✔️✔️:                              -1.0,     // pruneRmsThresh
        // RDKit✔️✔️:                              true,     // onlyHeavyAtomsForRMS
        // RDKit✔️✔️:                              2,        // ETversion
        // RDKit✔️✔️:                              nullptr,  // boundsMat
        // RDKit✔️✔️:                              true,     // embedFragmentsSeparately
        // RDKit✔️✔️:                              false,    // useSmallRingTorsions
        // RDKit✔️✔️:                              false,    // useMacrocycleTorsions
        // RDKit✔️✔️:                              false,    // useMacrocycle14config
        // RDKit✔️✔️:                              0,        // timeout
        // RDKit✔️✔️:                              nullptr,  // CPCI
        // RDKit✔️✔️:                              nullptr   // callback
        Self::from_rdkit_constructor(
            0, 1, -1, true, false, 2.0, true, 1, None, 1e-3, false, false, true, false, false, 5.0,
            -1.0, true, 2, None, true, false, false, false, 0, None, None,
        )
    }

    #[must_use]
    pub fn etkdg() -> Self {
        // RDKit✔️✔️: const EmbedParameters ETKDG(0,        // maxIterations
        // RDKit✔️✔️:                             1,        // numThreads
        // RDKit✔️✔️:                             -1,       // randomSeed
        // RDKit✔️✔️:                             true,     // clearConfs
        // RDKit✔️✔️:                             false,    // useRandomCoords
        // RDKit✔️✔️:                             2.0,      // boxSizeMult
        // RDKit✔️✔️:                             true,     // randNegEig
        // RDKit✔️✔️:                             1,        // numZeroFail
        // RDKit✔️✔️:                             nullptr,  // coordMap
        // RDKit✔️✔️:                             1e-3,     // optimizerForceTol
        // RDKit✔️✔️:                             false,    // ignoreSmoothingFailures
        // RDKit✔️✔️:                             true,     // enforceChirality
        // RDKit✔️✔️:                             true,     // useExpTorsionAnglePrefs
        // RDKit✔️✔️:                             true,     // useBasicKnowledge
        // RDKit✔️✔️:                             false,    // verbose
        // RDKit✔️✔️:                             5.0,      // basinThresh
        // RDKit✔️✔️:                             -1.0,     // pruneRmsThresh
        // RDKit✔️✔️:                             true,     // onlyHeavyAtomsForRMS
        // RDKit✔️✔️:                             1,        // ETversion
        // RDKit✔️✔️:                             nullptr,  // boundsMat
        // RDKit✔️✔️:                             true,     // embedFragmentsSeparately
        // RDKit✔️✔️:                             false,    // useSmallRingTorsions
        // RDKit✔️✔️:                             false,    // useMacrocycleTorsions
        // RDKit✔️✔️:                             false,    // useMacrocycle14config
        // RDKit✔️✔️:                             0,        // timeout
        // RDKit✔️✔️:                             nullptr,  // CPCI
        // RDKit✔️✔️:                             nullptr   // callback
        Self::from_rdkit_constructor(
            0, 1, -1, true, false, 2.0, true, 1, None, 1e-3, false, true, true, true, false, 5.0,
            -1.0, true, 1, None, true, false, false, false, 0, None, None,
        )
    }

    #[must_use]
    pub fn etkdg_v2() -> Self {
        // RDKit✔️✔️: const EmbedParameters ETKDGv2(0,        // maxIterations
        // RDKit✔️✔️:                               1,        // numThreads
        // RDKit✔️✔️:                               -1,       // randomSeed
        // RDKit✔️✔️:                               true,     // clearConfs
        // RDKit✔️✔️:                               false,    // useRandomCoords
        // RDKit✔️✔️:                               2.0,      // boxSizeMult
        // RDKit✔️✔️:                               true,     // randNegEig
        // RDKit✔️✔️:                               1,        // numZeroFail
        // RDKit✔️✔️:                               nullptr,  // coordMap
        // RDKit✔️✔️:                               1e-3,     // optimizerForceTol
        // RDKit✔️✔️:                               false,    // ignoreSmoothingFailures
        // RDKit✔️✔️:                               true,     // enforceChirality
        // RDKit✔️✔️:                               true,     // useExpTorsionAnglePrefs
        // RDKit✔️✔️:                               true,     // useBasicKnowledge
        // RDKit✔️✔️:                               false,    // verbose
        // RDKit✔️✔️:                               5.0,      // basinThresh
        // RDKit✔️✔️:                               -1.0,     // pruneRmsThresh
        // RDKit✔️✔️:                               true,     // onlyHeavyAtomsForRMS
        // RDKit✔️✔️:                               2,        // ETversion
        // RDKit✔️✔️:                               nullptr,  // boundsMat
        // RDKit✔️✔️:                               true,     // embedFragmentsSeparately
        // RDKit✔️✔️:                               false,    // useSmallRingTorsions
        // RDKit✔️✔️:                               false,    // useMacrocycleTorsions
        // RDKit✔️✔️:                               false,    // useMacrocycle14config
        // RDKit✔️✔️:                               0,        // timeout
        // RDKit✔️✔️:                               nullptr,  // CPCI
        // RDKit✔️✔️:                               nullptr   // callback
        Self::from_rdkit_constructor(
            0, 1, -1, true, false, 2.0, true, 1, None, 1e-3, false, true, true, true, false, 5.0,
            -1.0, true, 2, None, true, false, false, false, 0, None, None,
        )
    }

    #[must_use]
    pub fn etkdg_v3() -> Self {
        // RDKit✔️✔️: const EmbedParameters ETKDGv3(0,        // maxIterations
        // RDKit✔️✔️:                               1,        // numThreads
        // RDKit✔️✔️:                               -1,       // randomSeed
        // RDKit✔️✔️:                               true,     // clearConfs
        // RDKit✔️✔️:                               false,    // useRandomCoords
        // RDKit✔️✔️:                               2.0,      // boxSizeMult
        // RDKit✔️✔️:                               true,     // randNegEig
        // RDKit✔️✔️:                               1,        // numZeroFail
        // RDKit✔️✔️:                               nullptr,  // coordMap
        // RDKit✔️✔️:                               1e-3,     // optimizerForceTol
        // RDKit✔️✔️:                               false,    // ignoreSmoothingFailures
        // RDKit✔️✔️:                               true,     // enforceChirality
        // RDKit✔️✔️:                               true,     // useExpTorsionAnglePrefs
        // RDKit✔️✔️:                               true,     // useBasicKnowledge
        // RDKit✔️✔️:                               false,    // verbose
        // RDKit✔️✔️:                               5.0,      // basinThresh
        // RDKit✔️✔️:                               -1.0,     // pruneRmsThresh
        // RDKit✔️✔️:                               true,     // onlyHeavyAtomsForRMS
        // RDKit✔️✔️:                               2,        // ETversion
        // RDKit✔️✔️:                               nullptr,  // boundsMat
        // RDKit✔️✔️:                               true,     // embedFragmentsSeparately
        // RDKit✔️✔️:                               false,    // useSmallRingTorsions
        // RDKit✔️✔️:                               true,     // useMacrocycleTorsions
        // RDKit✔️✔️:                               true,     // useMacrocycle14config
        // RDKit✔️✔️:                               0,        // timeout
        // RDKit✔️✔️:                               nullptr,  // CPCI
        // RDKit✔️✔️:                               nullptr   // callback
        Self::from_rdkit_constructor(
            0, 1, -1, true, false, 2.0, true, 1, None, 1e-3, false, true, true, true, false, 5.0,
            -1.0, true, 2, None, true, false, true, true, 0, None, None,
        )
    }

    #[must_use]
    pub fn sr_etkdg_v3() -> Self {
        // RDKit✔️✔️: const EmbedParameters srETKDGv3(0,        // maxIterations
        // RDKit✔️✔️:                                 1,        // numThreads
        // RDKit✔️✔️:                                 -1,       // randomSeed
        // RDKit✔️✔️:                                 true,     // clearConfs
        // RDKit✔️✔️:                                 false,    // useRandomCoords
        // RDKit✔️✔️:                                 2.0,      // boxSizeMult
        // RDKit✔️✔️:                                 true,     // randNegEig
        // RDKit✔️✔️:                                 1,        // numZeroFail
        // RDKit✔️✔️:                                 nullptr,  // coordMap
        // RDKit✔️✔️:                                 1e-3,     // optimizerForceTol
        // RDKit✔️✔️:                                 false,    // ignoreSmoothingFailures
        // RDKit✔️✔️:                                 true,     // enforceChirality
        // RDKit✔️✔️:                                 true,     // useExpTorsionAnglePrefs
        // RDKit✔️✔️:                                 true,     // useBasicKnowledge
        // RDKit✔️✔️:                                 false,    // verbose
        // RDKit✔️✔️:                                 5.0,      // basinThresh
        // RDKit✔️✔️:                                 -1.0,     // pruneRmsThresh
        // RDKit✔️✔️:                                 true,     // onlyHeavyAtomsForRMS
        // RDKit✔️✔️:                                 2,        // ETversion
        // RDKit✔️✔️:                                 nullptr,  // boundsMat
        // RDKit✔️✔️:                                 true,     // embedFragmentsSeparately
        // RDKit✔️✔️:                                 true,     // useSmallRingTorsions
        // RDKit✔️✔️:                                 false,    // useMacrocycleTorsions
        // RDKit✔️✔️:                                 false,    // useMacrocycle14config
        // RDKit✔️✔️:                                 0,        // timeout
        // RDKit✔️✔️:                                 nullptr,  // CPCI
        // RDKit✔️✔️:                                 nullptr   // callback
        Self::from_rdkit_constructor(
            0, 1, -1, true, false, 2.0, true, 1, None, 1e-3, false, true, true, true, false, 5.0,
            -1.0, true, 2, None, true, true, false, false, 0, None, None,
        )
    }

    pub fn update_from_json(&mut self, json: &str) -> Result<(), DgBoundsError> {
        // BEGIN RDKIT CPP FUNCTION DGeomHelpers::updateEmbedParametersFromJSON (EmbedderUtils.cpp:56-87)
        // RDKit✔️✔️: void updateEmbedParametersFromJSON(EmbedParameters &params,
        // RDKit✔️✔️:                                    const std::string &json) {
        // RDKit✔️✔️:   if (json.empty()) {
        // RDKit✔️✔️:     return;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   std::istringstream ss(json);
        // RDKit✔️✔️:   boost::property_tree::ptree pt;
        // RDKit✔️✔️:   boost::property_tree::read_json(ss, pt);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   EMBED_PARAMS_FIELDS(PT_OPT_GET)
        // RDKit✔️✔️:
        // RDKit✔️✔️:   std::map<int, RDGeom::Point3D> *cmap = nullptr;
        // RDKit✔️✔️:   const auto coordMap = pt.get_child_optional("coordMap");
        // RDKit✔️✔️:   if (coordMap) {
        // RDKit✔️✔️:     // NOTE: this leaks since EmbedParameters uses a naked pointer and we don't
        // RDKit✔️✔️:     // have any way to tie the lifetime of the memory we allocate here to the
        // RDKit✔️✔️:     // EmbedParameters object itself.
        // RDKit✔️✔️:     cmap = new std::map<int, RDGeom::Point3D>();
        // RDKit✔️✔️:     for (const auto &entry : *coordMap) {
        // RDKit✔️✔️:       RDGeom::Point3D pt;
        // RDKit✔️✔️:
        // RDKit✔️✔️:       auto itm = entry.second.begin();
        // RDKit✔️✔️:       pt.x = itm->second.get_value<float>();
        // RDKit✔️✔️:       ++itm;
        // RDKit✔️✔️:       pt.y = itm->second.get_value<float>();
        // RDKit✔️✔️:       ++itm;
        // RDKit✔️✔️:       pt.z = itm->second.get_value<float>();
        // RDKit✔️✔️:
        // RDKit✔️✔️:       (*cmap)[boost::lexical_cast<int>(entry.first)] = pt;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     params.coordMap = cmap;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION DGeomHelpers::updateEmbedParametersFromJSON
        if json.is_empty() {
            return Ok(());
        }
        let value: serde_json::Value = serde_json::from_str(json)
            .map_err(|err| DgBoundsError::InvalidEmbedParametersJson(err.to_string()))?;

        update_f64_field(&value, "basinThresh", &mut self.basin_thresh)?;
        update_f64_field(
            &value,
            "boundsMatForceScaling",
            &mut self.bounds_mat_force_scaling,
        )?;
        update_f64_field(&value, "boxSizeMult", &mut self.box_size_mult)?;
        update_bool_field(&value, "clearConfs", &mut self.clear_confs)?;
        update_bool_field(
            &value,
            "embedFragmentsSeparately",
            &mut self.embed_fragments_separately,
        )?;
        update_bool_field(
            &value,
            "enableSequentialRandomSeeds",
            &mut self.enable_sequential_random_seeds,
        )?;
        update_bool_field(&value, "enforceChirality", &mut self.enforce_chirality)?;
        update_u32_field(&value, "ETversion", &mut self.et_version)?;
        update_bool_field(&value, "forceTransAmides", &mut self.force_trans_amides)?;
        update_bool_field(
            &value,
            "ignoreSmoothingFailures",
            &mut self.ignore_smoothing_failures,
        )?;
        update_u32_field(&value, "maxIterations", &mut self.max_iterations)?;
        update_i32_field(&value, "numThreads", &mut self.num_threads)?;
        update_u32_field(&value, "numZeroFail", &mut self.num_zero_fail)?;
        update_bool_field(
            &value,
            "onlyHeavyAtomsForRMS",
            &mut self.only_heavy_atoms_for_rms,
        )?;
        update_f64_field(&value, "optimizerForceTol", &mut self.optimizer_force_tol)?;
        update_f64_field(&value, "pruneRmsThresh", &mut self.prune_rms_thresh)?;
        update_bool_field(&value, "randNegEig", &mut self.rand_neg_eig)?;
        update_i32_field(&value, "randomSeed", &mut self.random_seed)?;
        update_bool_field(
            &value,
            "symmetrizeConjugatedTerminalGroupsForPruning",
            &mut self.symmetrize_conjugated_terminal_groups_for_pruning,
        )?;
        update_u32_field(&value, "timeout", &mut self.timeout)?;
        update_bool_field(&value, "trackFailures", &mut self.track_failures)?;
        update_bool_field(&value, "useBasicKnowledge", &mut self.use_basic_knowledge)?;
        update_bool_field(
            &value,
            "useExpTorsionAnglePrefs",
            &mut self.use_exp_torsion_angle_prefs,
        )?;
        update_bool_field(
            &value,
            "useMacrocycle14config",
            &mut self.use_macrocycle14config,
        )?;
        update_bool_field(
            &value,
            "useMacrocycleTorsions",
            &mut self.use_macrocycle_torsions,
        )?;
        update_bool_field(&value, "useRandomCoords", &mut self.use_random_coords)?;
        update_bool_field(
            &value,
            "useSmallRingTorsions",
            &mut self.use_small_ring_torsions,
        )?;
        update_bool_field(
            &value,
            "useSymmetryForPruning",
            &mut self.use_symmetry_for_pruning,
        )?;
        update_bool_field(&value, "verbose", &mut self.verbose)?;

        if let Some(coord_map) = value.get("coordMap") {
            self.coord_map = Some(parse_embed_parameters_coord_map(coord_map)?);
        }

        Ok(())
    }

    #[must_use]
    pub fn to_json(&self) -> String {
        // BEGIN RDKIT CPP FUNCTION DGeomHelpers::embedParametersToJSON (EmbedderUtils.cpp:90-126)
        // RDKit✔️✔️: std::string embedParametersToJSON(const EmbedParameters &params) {
        // RDKit✔️✔️:   boost::property_tree::ptree pt;
        // RDKit✔️✔️:
        // RDKit✔️✔️:   EMBED_PARAMS_FIELDS(PT_OPT_PUT)
        // RDKit✔️✔️:
        // RDKit✔️✔️:   if (params.coordMap) {
        // RDKit✔️✔️:     boost::property_tree::ptree coordMapPT;
        // RDKit✔️✔️:
        // RDKit✔️✔️:     for (const auto &kv : *params.coordMap) {
        // RDKit✔️✔️:       boost::property_tree::ptree pointPT;
        // RDKit✔️✔️:       pointPT.push_back(
        // RDKit✔️✔️:           {"", boost::property_tree::ptree(std::to_string(kv.second.x))});
        // RDKit✔️✔️:       pointPT.push_back(
        // RDKit✔️✔️:           {"", boost::property_tree::ptree(std::to_string(kv.second.y))});
        // RDKit✔️✔️:       pointPT.push_back(
        // RDKit✔️✔️:           {"", boost::property_tree::ptree(std::to_string(kv.second.z))});
        // RDKit✔️✔️:
        // RDKit✔️✔️:       coordMapPT.add_child(std::to_string(kv.first), pointPT);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:
        // RDKit✔️✔️:     pt.add_child("coordMap", coordMapPT);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   if (params.boundsMat) {
        // RDKit✔️✔️:     boost::property_tree::ptree matrixPT;
        // RDKit✔️✔️:     const unsigned int N = params.boundsMat->numCols();
        // RDKit✔️✔️:     for (unsigned i = 0; i < N; ++i) {
        // RDKit✔️✔️:       boost::property_tree::ptree rowPT;
        // RDKit✔️✔️:
        // RDKit✔️✔️:       for (unsigned j = 0; j < N; ++j) {
        // RDKit✔️✔️:         boost::property_tree::ptree v;
        // RDKit✔️✔️:         v.put("", params.boundsMat->getVal(i, j));
        // RDKit✔️✔️:         rowPT.push_back({"", v});
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:
        // RDKit✔️✔️:       matrixPT.push_back({"", rowPT});
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:
        // RDKit✔️✔️:     pt.add_child("boundsMatrix", matrixPT);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   std::ostringstream ss;
        // RDKit✔️✔️:   boost::property_tree::write_json(ss, pt, false);
        // RDKit✔️✔️:   auto str = ss.str();
        // RDKit✔️✔️:   boost::algorithm::trim(str);
        // RDKit✔️✔️:   return str;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION DGeomHelpers::embedParametersToJSON
        let mut fields = Vec::with_capacity(31);
        push_json_field(&mut fields, "basinThresh", self.basin_thresh);
        push_json_field(
            &mut fields,
            "boundsMatForceScaling",
            self.bounds_mat_force_scaling,
        );
        push_json_field(&mut fields, "boxSizeMult", self.box_size_mult);
        push_json_field(&mut fields, "clearConfs", self.clear_confs);
        push_json_field(
            &mut fields,
            "embedFragmentsSeparately",
            self.embed_fragments_separately,
        );
        push_json_field(
            &mut fields,
            "enableSequentialRandomSeeds",
            self.enable_sequential_random_seeds,
        );
        push_json_field(&mut fields, "enforceChirality", self.enforce_chirality);
        push_json_field(&mut fields, "ETversion", self.et_version);
        push_json_field(&mut fields, "forceTransAmides", self.force_trans_amides);
        push_json_field(
            &mut fields,
            "ignoreSmoothingFailures",
            self.ignore_smoothing_failures,
        );
        push_json_field(&mut fields, "maxIterations", self.max_iterations);
        push_json_field(&mut fields, "numThreads", self.num_threads);
        push_json_field(&mut fields, "numZeroFail", self.num_zero_fail);
        push_json_field(
            &mut fields,
            "onlyHeavyAtomsForRMS",
            self.only_heavy_atoms_for_rms,
        );
        push_json_field(&mut fields, "optimizerForceTol", self.optimizer_force_tol);
        push_json_field(&mut fields, "pruneRmsThresh", self.prune_rms_thresh);
        push_json_field(&mut fields, "randNegEig", self.rand_neg_eig);
        push_json_field(&mut fields, "randomSeed", self.random_seed);
        push_json_field(
            &mut fields,
            "symmetrizeConjugatedTerminalGroupsForPruning",
            self.symmetrize_conjugated_terminal_groups_for_pruning,
        );
        push_json_field(&mut fields, "timeout", self.timeout);
        push_json_field(&mut fields, "trackFailures", self.track_failures);
        push_json_field(&mut fields, "useBasicKnowledge", self.use_basic_knowledge);
        push_json_field(
            &mut fields,
            "useExpTorsionAnglePrefs",
            self.use_exp_torsion_angle_prefs,
        );
        push_json_field(
            &mut fields,
            "useMacrocycle14config",
            self.use_macrocycle14config,
        );
        push_json_field(
            &mut fields,
            "useMacrocycleTorsions",
            self.use_macrocycle_torsions,
        );
        push_json_field(&mut fields, "useRandomCoords", self.use_random_coords);
        push_json_field(
            &mut fields,
            "useSmallRingTorsions",
            self.use_small_ring_torsions,
        );
        push_json_field(
            &mut fields,
            "useSymmetryForPruning",
            self.use_symmetry_for_pruning,
        );
        push_json_field(&mut fields, "verbose", self.verbose);

        if let Some(coord_map) = &self.coord_map {
            let entries = coord_map
                .iter()
                .map(|(atom_idx, point)| {
                    format!(
                        "\"{}\":[\"{:.6}\",\"{:.6}\",\"{:.6}\"]",
                        atom_idx, point.x, point.y, point.z
                    )
                })
                .collect::<Vec<_>>()
                .join(",");
            fields.push(format!("\"coordMap\":{{{entries}}}"));
        }

        if let Some(bounds_mat) = &self.bounds_mat {
            let rows = (0..bounds_mat.n)
                .map(|i| {
                    let row = (0..bounds_mat.n)
                        .map(|j| format!("\"{}\"", bounds_mat.get_val(i, j)))
                        .collect::<Vec<_>>()
                        .join(",");
                    format!("[{row}]")
                })
                .collect::<Vec<_>>()
                .join(",");
            fields.push(format!("\"boundsMatrix\":[{rows}]"));
        }

        format!("{{{}}}", fields.join(","))
    }
}

fn embed_parameters_json_field<'a>(
    value: &'a serde_json::Value,
    name: &str,
) -> Option<&'a serde_json::Value> {
    value.as_object().and_then(|object| object.get(name))
}

fn push_json_field<T: ToString>(fields: &mut Vec<String>, name: &str, value: T) {
    fields.push(format!("\"{}\":\"{}\"", name, value.to_string()));
}

fn embed_parameters_invalid_json(name: &str, expected: &str) -> DgBoundsError {
    DgBoundsError::InvalidEmbedParametersJson(format!("{name} must be {expected}"))
}

fn json_value_as_f64(name: &str, value: &serde_json::Value) -> Result<f64, DgBoundsError> {
    if let Some(number) = value.as_f64() {
        Ok(number)
    } else if let Some(text) = value.as_str() {
        text.parse::<f64>()
            .map_err(|_| embed_parameters_invalid_json(name, "a floating-point number"))
    } else {
        Err(embed_parameters_invalid_json(
            name,
            "a floating-point number",
        ))
    }
}

fn json_value_as_i32(name: &str, value: &serde_json::Value) -> Result<i32, DgBoundsError> {
    if let Some(number) = value.as_i64() {
        i32::try_from(number).map_err(|_| embed_parameters_invalid_json(name, "a 32-bit integer"))
    } else if let Some(text) = value.as_str() {
        text.parse::<i32>()
            .map_err(|_| embed_parameters_invalid_json(name, "a 32-bit integer"))
    } else {
        Err(embed_parameters_invalid_json(name, "a 32-bit integer"))
    }
}

fn json_value_as_u32(name: &str, value: &serde_json::Value) -> Result<u32, DgBoundsError> {
    if let Some(number) = value.as_u64() {
        u32::try_from(number)
            .map_err(|_| embed_parameters_invalid_json(name, "an unsigned 32-bit integer"))
    } else if let Some(text) = value.as_str() {
        text.parse::<u32>()
            .map_err(|_| embed_parameters_invalid_json(name, "an unsigned 32-bit integer"))
    } else {
        Err(embed_parameters_invalid_json(
            name,
            "an unsigned 32-bit integer",
        ))
    }
}

fn json_value_as_bool(name: &str, value: &serde_json::Value) -> Result<bool, DgBoundsError> {
    if let Some(value) = value.as_bool() {
        Ok(value)
    } else if let Some(number) = value.as_u64() {
        match number {
            0 => Ok(false),
            1 => Ok(true),
            _ => Err(embed_parameters_invalid_json(name, "a boolean")),
        }
    } else if let Some(text) = value.as_str() {
        match text {
            "0" => Ok(false),
            "1" => Ok(true),
            "false" => Ok(false),
            "true" => Ok(true),
            _ => Err(embed_parameters_invalid_json(name, "a boolean")),
        }
    } else {
        Err(embed_parameters_invalid_json(name, "a boolean"))
    }
}

fn update_f64_field(
    value: &serde_json::Value,
    name: &str,
    target: &mut f64,
) -> Result<(), DgBoundsError> {
    if let Some(field) = embed_parameters_json_field(value, name) {
        *target = json_value_as_f64(name, field)?;
    }
    Ok(())
}

fn update_i32_field(
    value: &serde_json::Value,
    name: &str,
    target: &mut i32,
) -> Result<(), DgBoundsError> {
    if let Some(field) = embed_parameters_json_field(value, name) {
        *target = json_value_as_i32(name, field)?;
    }
    Ok(())
}

fn update_u32_field(
    value: &serde_json::Value,
    name: &str,
    target: &mut u32,
) -> Result<(), DgBoundsError> {
    if let Some(field) = embed_parameters_json_field(value, name) {
        *target = json_value_as_u32(name, field)?;
    }
    Ok(())
}

fn update_bool_field(
    value: &serde_json::Value,
    name: &str,
    target: &mut bool,
) -> Result<(), DgBoundsError> {
    if let Some(field) = embed_parameters_json_field(value, name) {
        *target = json_value_as_bool(name, field)?;
    }
    Ok(())
}

fn parse_embed_parameters_coord_map(
    value: &serde_json::Value,
) -> Result<BTreeMap<i32, ForceFieldVec3>, DgBoundsError> {
    let object = value
        .as_object()
        .ok_or_else(|| embed_parameters_invalid_json("coordMap", "an object"))?;
    let mut coord_map = BTreeMap::new();
    for (key, point_value) in object {
        let atom_idx = key
            .parse::<i32>()
            .map_err(|_| embed_parameters_invalid_json("coordMap key", "a 32-bit integer"))?;
        let point = point_value
            .as_array()
            .ok_or_else(|| embed_parameters_invalid_json("coordMap value", "an array"))?;
        if point.len() < 3 {
            return Err(embed_parameters_invalid_json(
                "coordMap value",
                "an array with at least three coordinates",
            ));
        }
        coord_map.insert(
            atom_idx,
            ForceFieldVec3::new(
                json_value_as_f64("coordMap x", &point[0])?,
                json_value_as_f64("coordMap y", &point[1])?,
                json_value_as_f64("coordMap z", &point[2])?,
            ),
        );
    }
    Ok(coord_map)
}

// RDKit✔️✔️:   EmbedParameters(
// RDKit✔️✔️:       unsigned int maxIterations, int numThreads, int randomSeed,
// RDKit✔️✔️:       bool clearConfs, bool useRandomCoords, double boxSizeMult,
// RDKit✔️✔️:       bool randNegEig, unsigned int numZeroFail,
// RDKit✔️✔️:       const std::map<int, RDGeom::Point3D> *coordMap, double optimizerForceTol,
// RDKit✔️✔️:       bool ignoreSmoothingFailures, bool enforceChirality,
// RDKit✔️✔️:       bool useExpTorsionAnglePrefs, bool useBasicKnowledge, bool verbose,
// RDKit✔️✔️:       double basinThresh, double pruneRmsThresh, bool onlyHeavyAtomsForRMS,
// RDKit✔️✔️:       unsigned int ETversion = 2,
// RDKit✔️✔️:       const DistGeom::BoundsMatrix *boundsMat = nullptr,
// RDKit✔️✔️:       bool embedFragmentsSeparately = true, bool useSmallRingTorsions = false,
// RDKit✔️✔️:       bool useMacrocycleTorsions = false, bool useMacrocycle14config = false,
// RDKit✔️✔️:       unsigned int timeout = 0,
// RDKit✔️✔️:       std::shared_ptr<std::map<std::pair<unsigned int, unsigned int>, double>>
// RDKit✔️✔️:           CPCI = nullptr,
// RDKit✔️✔️:       void (*callback)(unsigned int) = nullptr)
// RDKit✔️✔️:       : maxIterations(maxIterations),
// RDKit✔️✔️:         numThreads(numThreads),
// RDKit✔️✔️:         randomSeed(randomSeed),
// RDKit✔️✔️:         clearConfs(clearConfs),
// RDKit✔️✔️:         useRandomCoords(useRandomCoords),
// RDKit✔️✔️:         boxSizeMult(boxSizeMult),
// RDKit✔️✔️:         randNegEig(randNegEig),
// RDKit✔️✔️:         numZeroFail(numZeroFail),
// RDKit✔️✔️:         coordMap(coordMap),
// RDKit✔️✔️:         optimizerForceTol(optimizerForceTol),
// RDKit✔️✔️:         ignoreSmoothingFailures(ignoreSmoothingFailures),
// RDKit✔️✔️:         enforceChirality(enforceChirality),
// RDKit✔️✔️:         useExpTorsionAnglePrefs(useExpTorsionAnglePrefs),
// RDKit✔️✔️:         useBasicKnowledge(useBasicKnowledge),
// RDKit✔️✔️:         verbose(verbose),
// RDKit✔️✔️:         basinThresh(basinThresh),
// RDKit✔️✔️:         pruneRmsThresh(pruneRmsThresh),
// RDKit✔️✔️:         onlyHeavyAtomsForRMS(onlyHeavyAtomsForRMS),
// RDKit✔️✔️:         ETversion(ETversion),
// RDKit✔️✔️:         boundsMat(boundsMat),
// RDKit✔️✔️:         embedFragmentsSeparately(embedFragmentsSeparately),
// RDKit✔️✔️:         useSmallRingTorsions(useSmallRingTorsions),
// RDKit✔️✔️:         useMacrocycleTorsions(useMacrocycleTorsions),
// RDKit✔️✔️:         useMacrocycle14config(useMacrocycle14config),
// RDKit✔️✔️:         timeout(timeout),
// RDKit✔️✔️:         CPCI(std::move(CPCI)),
// RDKit✔️✔️:         callback(callback) {}
// RDKit✔️✔️: };
// END RDKIT CPP STRUCT DGeomHelpers::EmbedParameters

// ──────────────────────────────────────────────
// Chiral sets
// ──────────────────────────────────────────────

// BEGIN RDKIT CPP ENUM DistGeom::ChiralSetStructureFlags (ChiralSet.h:18-21)
// RDKit✔️✔️: enum class ChiralSetStructureFlags : std::uint64_t {
// RDKit✔️✔️:   IN_FUSED_SMALL_RINGS =
// RDKit✔️✔️:       1 << 0,  // a chiral center involved in fusing 2 or more small rings
// RDKit✔️✔️: };
// END RDKIT CPP ENUM DistGeom::ChiralSetStructureFlags
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u64)]
pub enum ChiralSetStructureFlags {
    InFusedSmallRings = 1 << 0,
}

/// Class used to store a quartet of points and chiral volume bounds on them.
#[derive(Debug, Clone, PartialEq)]
pub struct ChiralSet {
    pub idx0: usize,
    pub idx1: usize,
    pub idx2: usize,
    pub idx3: usize,
    pub idx4: usize,
    pub volume_lower_bound: f64,
    pub volume_upper_bound: f64,
    pub structure_flags: u64,
}

impl ChiralSet {
    #[must_use]
    pub fn new(
        pid0: usize,
        pid1: usize,
        pid2: usize,
        pid3: usize,
        pid4: usize,
        lower_vol_bound: f64,
        upper_vol_bound: f64,
        structure_flags: u64,
    ) -> Self {
        // BEGIN RDKIT CPP CONSTRUCTOR DistGeom::ChiralSet::ChiralSet (ChiralSet.h:40-57)
        // RDKit✔️✔️: ChiralSet(unsigned int pid0, unsigned int pid1, unsigned int pid2,
        // RDKit✔️✔️:           unsigned int pid3, unsigned int pid4, double lowerVolBound,
        // RDKit✔️✔️:           double upperVolBound, std::uint64_t structureFlags = 0)
        // RDKit✔️✔️:     : d_idx0(pid0),
        // RDKit✔️✔️:       d_idx1(pid1),
        // RDKit✔️✔️:       d_idx2(pid2),
        // RDKit✔️✔️:       d_idx3(pid3),
        // RDKit✔️✔️:       d_idx4(pid4),
        // RDKit✔️✔️:       d_volumeLowerBound(lowerVolBound),
        // RDKit✔️✔️:       d_volumeUpperBound(upperVolBound),
        // RDKit✔️✔️:       d_structureFlags(structureFlags) {
        // RDKit✔️✔️:   CHECK_INVARIANT(lowerVolBound <= upperVolBound, "Inconsistent bounds\n");
        // RDKit✔️✔️:   d_volumeLowerBound = lowerVolBound;
        // RDKit✔️✔️:   d_volumeUpperBound = upperVolBound;
        // RDKit✔️✔️: }
        // END RDKIT CPP CONSTRUCTOR DistGeom::ChiralSet::ChiralSet
        assert!(lower_vol_bound <= upper_vol_bound, "Inconsistent bounds\n");
        Self {
            idx0: pid0,
            idx1: pid1,
            idx2: pid2,
            idx3: pid3,
            idx4: pid4,
            volume_lower_bound: lower_vol_bound,
            volume_upper_bound: upper_vol_bound,
            structure_flags,
        }
    }

    #[must_use]
    pub fn with_default_structure_flags(
        pid0: usize,
        pid1: usize,
        pid2: usize,
        pid3: usize,
        pid4: usize,
        lower_vol_bound: f64,
        upper_vol_bound: f64,
    ) -> Self {
        Self::new(
            pid0,
            pid1,
            pid2,
            pid3,
            pid4,
            lower_vol_bound,
            upper_vol_bound,
            0,
        )
    }

    #[must_use]
    pub fn get_upper_volume_bound(&self) -> f64 {
        // BEGIN RDKIT CPP METHOD DistGeom::ChiralSet::getUpperVolumeBound (ChiralSet.h:59)
        // RDKit✔️✔️: inline double getUpperVolumeBound() const { return d_volumeUpperBound; }
        // END RDKIT CPP METHOD DistGeom::ChiralSet::getUpperVolumeBound
        self.volume_upper_bound
    }

    #[must_use]
    pub fn get_lower_volume_bound(&self) -> f64 {
        // BEGIN RDKIT CPP METHOD DistGeom::ChiralSet::getLowerVolumeBound (ChiralSet.h:61)
        // RDKit✔️✔️: inline double getLowerVolumeBound() const { return d_volumeLowerBound; }
        // END RDKIT CPP METHOD DistGeom::ChiralSet::getLowerVolumeBound
        self.volume_lower_bound
    }
}

// BEGIN RDKIT CPP TYPEDEFS DistGeom::ChiralSetPtr/VECT_CHIRALSET (ChiralSet.h:64-65)
// RDKit✔️✔️: typedef boost::shared_ptr<ChiralSet> ChiralSetPtr;
// RDKit✔️✔️: typedef std::vector<ChiralSetPtr> VECT_CHIRALSET;
// END RDKIT CPP TYPEDEFS DistGeom::ChiralSetPtr/VECT_CHIRALSET
pub type ChiralSetPtr = Arc<ChiralSet>;
pub type VectChiralSet = Vec<ChiralSetPtr>;

#[must_use]
pub fn calc_chiral_volume_flat(
    idx1: usize,
    idx2: usize,
    idx3: usize,
    idx4: usize,
    pos: &[f64],
    dim: usize,
) -> f64 {
    // BEGIN RDKIT CPP FUNCTION DistGeom::calcChiralVolume flat overload (ChiralViolationContribs.cpp:15-35)
    // RDKit✔️✔️: double calcChiralVolume(const unsigned int idx1, const unsigned int idx2,
    // RDKit✔️✔️:                         const unsigned int idx3, const unsigned int idx4,
    // RDKit✔️✔️:                         const double *pos, const unsigned int dim) {
    // RDKit✔️✔️:   // even if we are minimizing in higher dimension the chiral volume is
    // RDKit✔️✔️:   // calculated using only the first 3 dimensions
    // RDKit✔️✔️:   RDGeom::Point3D v1(pos[idx1 * dim] - pos[idx4 * dim],
    // RDKit✔️✔️:                      pos[idx1 * dim + 1] - pos[idx4 * dim + 1],
    // RDKit✔️✔️:                      pos[idx1 * dim + 2] - pos[idx4 * dim + 2]);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   RDGeom::Point3D v2(pos[idx2 * dim] - pos[idx4 * dim],
    // RDKit✔️✔️:                      pos[idx2 * dim + 1] - pos[idx4 * dim + 1],
    // RDKit✔️✔️:                      pos[idx2 * dim + 2] - pos[idx4 * dim + 2]);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   RDGeom::Point3D v3(pos[idx3 * dim] - pos[idx4 * dim],
    // RDKit✔️✔️:                      pos[idx3 * dim + 1] - pos[idx4 * dim + 1],
    // RDKit✔️✔️:                      pos[idx3 * dim + 2] - pos[idx4 * dim + 2]);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   RDGeom::Point3D v2xv3 = v2.crossProduct(v3);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   double vol = v1.dotProduct(v2xv3);
    // RDKit✔️✔️:   return vol;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DistGeom::calcChiralVolume flat overload
    let v1 = ForceFieldVec3::new(
        pos[idx1 * dim] - pos[idx4 * dim],
        pos[idx1 * dim + 1] - pos[idx4 * dim + 1],
        pos[idx1 * dim + 2] - pos[idx4 * dim + 2],
    );
    let v2 = ForceFieldVec3::new(
        pos[idx2 * dim] - pos[idx4 * dim],
        pos[idx2 * dim + 1] - pos[idx4 * dim + 1],
        pos[idx2 * dim + 2] - pos[idx4 * dim + 2],
    );
    let v3 = ForceFieldVec3::new(
        pos[idx3 * dim] - pos[idx4 * dim],
        pos[idx3 * dim + 1] - pos[idx4 * dim + 1],
        pos[idx3 * dim + 2] - pos[idx4 * dim + 2],
    );

    v1.dot_product(v2.cross_product(v3))
}

#[must_use]
pub fn calc_chiral_volume_points(
    idx1: usize,
    idx2: usize,
    idx3: usize,
    idx4: usize,
    pts: &[ForceFieldVec3],
) -> f64 {
    // BEGIN RDKIT CPP FUNCTION DistGeom::calcChiralVolume PointPtrVect overload (ChiralViolationContribs.cpp:36-56)
    // RDKit✔️✔️: double calcChiralVolume(const unsigned int idx1, const unsigned int idx2,
    // RDKit✔️✔️:                         const unsigned int idx3, const unsigned int idx4,
    // RDKit✔️✔️:                         const RDGeom::PointPtrVect &pts) {
    // RDKit✔️✔️:   // even if we are minimizing in higher dimension the chiral volume is
    // RDKit✔️✔️:   // calculated using only the first 3 dimensions
    // RDKit✔️✔️:   RDGeom::Point3D v1((*pts[idx1])[0] - (*pts[idx4])[0],
    // RDKit✔️✔️:                      (*pts[idx1])[1] - (*pts[idx4])[1],
    // RDKit✔️✔️:                      (*pts[idx1])[2] - (*pts[idx4])[2]);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   RDGeom::Point3D v2((*pts[idx2])[0] - (*pts[idx4])[0],
    // RDKit✔️✔️:                      (*pts[idx2])[1] - (*pts[idx4])[1],
    // RDKit✔️✔️:                      (*pts[idx2])[2] - (*pts[idx4])[2]);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   RDGeom::Point3D v3((*pts[idx3])[0] - (*pts[idx4])[0],
    // RDKit✔️✔️:                      (*pts[idx3])[1] - (*pts[idx4])[1],
    // RDKit✔️✔️:                      (*pts[idx3])[2] - (*pts[idx4])[2]);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   RDGeom::Point3D v2xv3 = v2.crossProduct(v3);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   double vol = v1.dotProduct(v2xv3);
    // RDKit✔️✔️:   return vol;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DistGeom::calcChiralVolume PointPtrVect overload
    let v1 = pts[idx1] - pts[idx4];
    let v2 = pts[idx2] - pts[idx4];
    let v3 = pts[idx3] - pts[idx4];

    v1.dot_product(v2.cross_product(v3))
}

#[derive(Debug, Clone, PartialEq)]
pub struct ChiralViolationContribsParams {
    pub idx1: usize,
    pub idx2: usize,
    pub idx3: usize,
    pub idx4: usize,
    pub vol_upper: f64,
    pub vol_lower: f64,
    pub weight: f64,
}

impl ChiralViolationContribsParams {
    #[must_use]
    pub fn new(
        idx1: usize,
        idx2: usize,
        idx3: usize,
        idx4: usize,
        vol_upper: f64,
        vol_lower: f64,
        weight: f64,
    ) -> Self {
        // BEGIN RDKIT CPP CONSTRUCTOR DistGeom::ChiralViolationContribsParams (ChiralViolationContribs.h:25-37)
        // RDKit✔️✔️: ChiralViolationContribsParams(unsigned int i1, unsigned int i2,
        // RDKit✔️✔️:                               unsigned int i3, unsigned int i4, double u,
        // RDKit✔️✔️:                               double l, double w = 1.0)
        // RDKit✔️✔️:     : idx1(i1),
        // RDKit✔️✔️:       idx2(i2),
        // RDKit✔️✔️:       idx3(i3),
        // RDKit✔️✔️:       idx4(i4),
        // RDKit✔️✔️:       volUpper(u),
        // RDKit✔️✔️:       volLower(l),
        // RDKit✔️✔️:       weight(w) {};
        // END RDKIT CPP CONSTRUCTOR DistGeom::ChiralViolationContribsParams
        Self {
            idx1,
            idx2,
            idx3,
            idx4,
            vol_upper,
            vol_lower,
            weight,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ChiralViolationContribs {
    owner: *const ForceField,
    contribs: Vec<ChiralViolationContribsParams>,
}

impl Default for ChiralViolationContribs {
    fn default() -> Self {
        Self {
            owner: std::ptr::null(),
            contribs: Vec::new(),
        }
    }
}

impl ChiralViolationContribs {
    #[must_use]
    pub fn new(owner: &ForceField) -> Self {
        // BEGIN RDKIT CPP CONSTRUCTOR DistGeom::ChiralViolationContribs::ChiralViolationContribs (ChiralViolationContribs.cpp:57-62)
        // RDKit✔️✔️: ChiralViolationContribs::ChiralViolationContribs(
        // RDKit✔️✔️:     ForceFields::ForceField *owner) {
        // RDKit✔️✔️:   PRECONDITION(owner, "bad owner");
        // RDKit✔️✔️:   dp_forceField = owner;
        // RDKit✔️✔️: }
        // END RDKIT CPP CONSTRUCTOR DistGeom::ChiralViolationContribs::ChiralViolationContribs
        Self {
            owner,
            contribs: Vec::new(),
        }
    }

    pub fn add_contrib(&mut self, cset: &ChiralSet, weight: f64) {
        // BEGIN RDKIT CPP METHOD DistGeom::ChiralViolationContribs::addContrib (ChiralViolationContribs.cpp:63-76)
        // RDKit✔️✔️: void ChiralViolationContribs::addContrib(const ChiralSet *cset, double weight) {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(cset, "bad chiral set");
        // RDKit✔️✔️:
        // RDKit✔️✔️:   URANGE_CHECK(cset->d_idx1, dp_forceField->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(cset->d_idx2, dp_forceField->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(cset->d_idx3, dp_forceField->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(cset->d_idx4, dp_forceField->positions().size());
        // RDKit✔️✔️:
        // RDKit✔️✔️:   d_contribs.emplace_back(cset->d_idx1, cset->d_idx2, cset->d_idx3,
        // RDKit✔️✔️:                           cset->d_idx4, cset->getUpperVolumeBound(),
        // RDKit✔️✔️:                           cset->getLowerVolumeBound(), weight);
        // RDKit✔️✔️: }
        // END RDKIT CPP METHOD DistGeom::ChiralViolationContribs::addContrib
        let owner = self.owner();
        assert!(cset.idx1 < owner.positions().len());
        assert!(cset.idx2 < owner.positions().len());
        assert!(cset.idx3 < owner.positions().len());
        assert!(cset.idx4 < owner.positions().len());
        self.contribs.push(ChiralViolationContribsParams::new(
            cset.idx1,
            cset.idx2,
            cset.idx3,
            cset.idx4,
            cset.get_upper_volume_bound(),
            cset.get_lower_volume_bound(),
            weight,
        ));
    }

    #[must_use]
    pub fn empty(&self) -> bool {
        // RDKit✔️✔️: bool empty() const { return d_contribs.empty(); }
        self.contribs.is_empty()
    }

    #[must_use]
    pub fn size(&self) -> usize {
        // RDKit✔️✔️: unsigned int size() const { return d_contribs.size(); }
        self.contribs.len()
    }

    #[must_use]
    pub fn contribs(&self) -> &[ChiralViolationContribsParams] {
        &self.contribs
    }

    fn owner(&self) -> &ForceField {
        assert!(!self.owner.is_null(), "no owner");
        // SAFETY: Matches the existing ForceFieldContrib owner-pointer convention.
        unsafe { &*self.owner }
    }
}

impl ForceFieldContrib for ChiralViolationContribs {
    fn copy(&self) -> Box<dyn ForceFieldContrib> {
        // RDKit✔️✔️: ChiralViolationContribs *copy() const override {
        // RDKit✔️✔️:   return new ChiralViolationContribs(*this);
        // RDKit✔️✔️: }
        Box::new(self.clone())
    }

    fn set_force_field(&mut self, owner: *const ForceField) {
        self.owner = owner;
    }

    fn get_energy(&self, pos: &[f64]) -> f64 {
        // BEGIN RDKIT CPP METHOD DistGeom::ChiralViolationContribs::getEnergy (ChiralViolationContribs.cpp:78-94)
        // RDKit✔️✔️: double ChiralViolationContribs::getEnergy(double *pos) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:
        // RDKit✔️✔️:   const unsigned int dim = dp_forceField->dimension();
        // RDKit✔️✔️:   double res = 0.0;
        // RDKit✔️✔️:   for (const auto &c : d_contribs) {
        // RDKit✔️✔️:     double vol = calcChiralVolume(c.idx1, c.idx2, c.idx3, c.idx4, pos, dim);
        // RDKit✔️✔️:     if (vol < c.volLower) {
        // RDKit✔️✔️:       res += c.weight * (vol - c.volLower) * (vol - c.volLower);
        // RDKit✔️✔️:     } else if (vol > c.volUpper) {
        // RDKit✔️✔️:       res += c.weight * (vol - c.volUpper) * (vol - c.volUpper);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        // END RDKIT CPP METHOD DistGeom::ChiralViolationContribs::getEnergy
        assert!(!pos.is_empty(), "bad vector");
        let dim = self.owner().dimension();
        let mut res = 0.0;
        for c in &self.contribs {
            let vol = calc_chiral_volume_flat(c.idx1, c.idx2, c.idx3, c.idx4, pos, dim);
            if vol < c.vol_lower {
                res += c.weight * (vol - c.vol_lower) * (vol - c.vol_lower);
            } else if vol > c.vol_upper {
                res += c.weight * (vol - c.vol_upper) * (vol - c.vol_upper);
            }
        }
        res
    }

    fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
        // BEGIN RDKIT CPP METHOD DistGeom::ChiralViolationContribs::getGrad (ChiralViolationContribs.cpp:96-174)
        // RDKit✔️✔️: void ChiralViolationContribs::getGrad(double *pos, double *grad) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:
        // RDKit✔️✔️:   const unsigned int dim = dp_forceField->dimension();
        // RDKit✔️✔️:
        // RDKit✔️✔️:   for (const auto &c : d_contribs) {
        // RDKit✔️✔️:     // even if we are minimizing in higher dimension the chiral volume is
        // RDKit✔️✔️:     // calculated using only the first 3 dimensions
        // RDKit✔️✔️:     RDGeom::Point3D v1(pos[c.idx1 * dim] - pos[c.idx4 * dim],
        // RDKit✔️✔️:                        pos[c.idx1 * dim + 1] - pos[c.idx4 * dim + 1],
        // RDKit✔️✔️:                        pos[c.idx1 * dim + 2] - pos[c.idx4 * dim + 2]);
        // RDKit✔️✔️:     RDGeom::Point3D v2(pos[c.idx2 * dim] - pos[c.idx4 * dim],
        // RDKit✔️✔️:                        pos[c.idx2 * dim + 1] - pos[c.idx4 * dim + 1],
        // RDKit✔️✔️:                        pos[c.idx2 * dim + 2] - pos[c.idx4 * dim + 2]);
        // RDKit✔️✔️:     RDGeom::Point3D v3(pos[c.idx3 * dim] - pos[c.idx4 * dim],
        // RDKit✔️✔️:                        pos[c.idx3 * dim + 1] - pos[c.idx4 * dim + 1],
        // RDKit✔️✔️:                        pos[c.idx3 * dim + 2] - pos[c.idx4 * dim + 2]);
        // RDKit✔️✔️:     RDGeom::Point3D v2xv3 = v2.crossProduct(v3);
        // RDKit✔️✔️:     double vol = v1.dotProduct(v2xv3);
        // RDKit✔️✔️:     double preFactor;
        // RDKit✔️✔️:     if (vol < c.volLower) {
        // RDKit✔️✔️:       preFactor = c.weight * (vol - c.volLower);
        // RDKit✔️✔️:     } else if (vol > c.volUpper) {
        // RDKit✔️✔️:       preFactor = c.weight * (vol - c.volUpper);
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        // END RDKIT CPP METHOD DistGeom::ChiralViolationContribs::getGrad
        assert!(!pos.is_empty(), "bad vector");
        assert!(!grad.is_empty(), "bad vector");
        let dim = self.owner().dimension();

        for c in &self.contribs {
            let v1 = ForceFieldVec3::new(
                pos[c.idx1 * dim] - pos[c.idx4 * dim],
                pos[c.idx1 * dim + 1] - pos[c.idx4 * dim + 1],
                pos[c.idx1 * dim + 2] - pos[c.idx4 * dim + 2],
            );
            let v2 = ForceFieldVec3::new(
                pos[c.idx2 * dim] - pos[c.idx4 * dim],
                pos[c.idx2 * dim + 1] - pos[c.idx4 * dim + 1],
                pos[c.idx2 * dim + 2] - pos[c.idx4 * dim + 2],
            );
            let v3 = ForceFieldVec3::new(
                pos[c.idx3 * dim] - pos[c.idx4 * dim],
                pos[c.idx3 * dim + 1] - pos[c.idx4 * dim + 1],
                pos[c.idx3 * dim + 2] - pos[c.idx4 * dim + 2],
            );

            let vol = v1.dot_product(v2.cross_product(v3));
            let pre_factor = if vol < c.vol_lower {
                c.weight * (vol - c.vol_lower)
            } else if vol > c.vol_upper {
                c.weight * (vol - c.vol_upper)
            } else {
                continue;
            };

            // RDKit✔️✔️:     grad[dim * c.idx1] += preFactor * ((v2.y) * (v3.z) - (v3.y) * (v2.z));
            // RDKit✔️✔️:     grad[dim * c.idx1 + 1] += preFactor * ((v3.x) * (v2.z) - (v2.x) * (v3.z));
            // RDKit✔️✔️:     grad[dim * c.idx1 + 2] += preFactor * ((v2.x) * (v3.y) - (v3.x) * (v2.y));
            grad[dim * c.idx1] += pre_factor * (v2.y * v3.z - v3.y * v2.z);
            grad[dim * c.idx1 + 1] += pre_factor * (v3.x * v2.z - v2.x * v3.z);
            grad[dim * c.idx1 + 2] += pre_factor * (v2.x * v3.y - v3.x * v2.y);

            // RDKit✔️✔️:     grad[dim * c.idx2] += preFactor * ((v3.y) * (v1.z) - (v3.z) * (v1.y));
            // RDKit✔️✔️:     grad[dim * c.idx2 + 1] += preFactor * ((v3.z) * (v1.x) - (v3.x) * (v1.z));
            // RDKit✔️✔️:     grad[dim * c.idx2 + 2] += preFactor * ((v3.x) * (v1.y) - (v3.y) * (v1.x));
            grad[dim * c.idx2] += pre_factor * (v3.y * v1.z - v3.z * v1.y);
            grad[dim * c.idx2 + 1] += pre_factor * (v3.z * v1.x - v3.x * v1.z);
            grad[dim * c.idx2 + 2] += pre_factor * (v3.x * v1.y - v3.y * v1.x);

            // RDKit✔️✔️:     grad[dim * c.idx3] += preFactor * ((v2.z) * (v1.y) - (v2.y) * (v1.z));
            // RDKit✔️✔️:     grad[dim * c.idx3 + 1] += preFactor * ((v2.x) * (v1.z) - (v2.z) * (v1.x));
            // RDKit✔️✔️:     grad[dim * c.idx3 + 2] += preFactor * ((v2.y) * (v1.x) - (v2.x) * (v1.y));
            grad[dim * c.idx3] += pre_factor * (v2.z * v1.y - v2.y * v1.z);
            grad[dim * c.idx3 + 1] += pre_factor * (v2.x * v1.z - v2.z * v1.x);
            grad[dim * c.idx3 + 2] += pre_factor * (v2.y * v1.x - v2.x * v1.y);

            // RDKit✔️✔️:     grad[dim * c.idx4] +=
            // RDKit✔️✔️:         preFactor * (pos[c.idx1 * dim + 2] *
            // RDKit✔️✔️:                          (pos[c.idx2 * dim + 1] - pos[c.idx3 * dim + 1]) +
            // RDKit✔️✔️:                      pos[c.idx2 * dim + 2] *
            // RDKit✔️✔️:                          (pos[c.idx3 * dim + 1] - pos[c.idx1 * dim + 1]) +
            // RDKit✔️✔️:                      pos[c.idx3 * dim + 2] *
            // RDKit✔️✔️:                          (pos[c.idx1 * dim + 1] - pos[c.idx2 * dim + 1]));
            grad[dim * c.idx4] += pre_factor
                * (pos[c.idx1 * dim + 2] * (pos[c.idx2 * dim + 1] - pos[c.idx3 * dim + 1])
                    + pos[c.idx2 * dim + 2] * (pos[c.idx3 * dim + 1] - pos[c.idx1 * dim + 1])
                    + pos[c.idx3 * dim + 2] * (pos[c.idx1 * dim + 1] - pos[c.idx2 * dim + 1]));

            // RDKit✔️✔️:     grad[dim * c.idx4 + 1] +=
            // RDKit✔️✔️:         preFactor *
            // RDKit✔️✔️:         (pos[c.idx1 * dim] * (pos[c.idx2 * dim + 2] - pos[c.idx3 * dim + 2]) +
            // RDKit✔️✔️:          pos[c.idx2 * dim] * (pos[c.idx3 * dim + 2] - pos[c.idx1 * dim + 2]) +
            // RDKit✔️✔️:          pos[c.idx3 * dim] * (pos[c.idx1 * dim + 2] - pos[c.idx2 * dim + 2]));
            grad[dim * c.idx4 + 1] += pre_factor
                * (pos[c.idx1 * dim] * (pos[c.idx2 * dim + 2] - pos[c.idx3 * dim + 2])
                    + pos[c.idx2 * dim] * (pos[c.idx3 * dim + 2] - pos[c.idx1 * dim + 2])
                    + pos[c.idx3 * dim] * (pos[c.idx1 * dim + 2] - pos[c.idx2 * dim + 2]));

            // RDKit✔️✔️:     grad[dim * c.idx4 + 2] +=
            // RDKit✔️✔️:         preFactor *
            // RDKit✔️✔️:         (pos[c.idx1 * dim + 1] * (pos[c.idx2 * dim] - pos[c.idx3 * dim]) +
            // RDKit✔️✔️:          pos[c.idx2 * dim + 1] * (pos[c.idx3 * dim] - pos[c.idx1 * dim]) +
            // RDKit✔️✔️:          pos[c.idx3 * dim + 1] * (pos[c.idx1 * dim] - pos[c.idx2 * dim]));
            grad[dim * c.idx4 + 2] += pre_factor
                * (pos[c.idx1 * dim + 1] * (pos[c.idx2 * dim] - pos[c.idx3 * dim])
                    + pos[c.idx2 * dim + 1] * (pos[c.idx3 * dim] - pos[c.idx1 * dim])
                    + pos[c.idx3 * dim + 1] * (pos[c.idx1 * dim] - pos[c.idx2 * dim]));
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct DistViolationContribsParams {
    pub idx1: usize,
    pub idx2: usize,
    pub ub: f64,
    pub lb: f64,
    pub ub2: f64,
    pub lb2: f64,
    pub weight: f64,
}

impl DistViolationContribsParams {
    #[must_use]
    pub fn new(idx1: usize, idx2: usize, ub: f64, lb: f64, weight: f64) -> Self {
        // BEGIN RDKIT CPP CONSTRUCTOR DistGeom::DistViolationContribsParams (DistViolationContribs.h:18-30)
        // RDKit✔️✔️: DistViolationContribsParams(unsigned int i1, unsigned int i2, double u,
        // RDKit✔️✔️:                             double l, double w = 1.0)
        // RDKit✔️✔️:     : idx1(i1), idx2(i2), ub(u), lb(l), ub2(u * u), lb2(l * l), weight(w) {};
        // END RDKIT CPP CONSTRUCTOR DistGeom::DistViolationContribsParams
        Self {
            idx1,
            idx2,
            ub,
            lb,
            ub2: ub * ub,
            lb2: lb * lb,
            weight,
        }
    }
}

// BEGIN RDKIT CPP LOCAL HELPER DistGeom::distance2 (DistViolationContribs.cpp:21-31)
// RDKit✔️✔️: inline double distance2(const unsigned int idx1, const unsigned int idx2,
// RDKit✔️✔️:                         const double *pos, const unsigned int dim) {
// RDKit✔️✔️:   const auto *end1Coords = &(pos[dim * idx1]);
// RDKit✔️✔️:   const auto *end2Coords = &(pos[dim * idx2]);
// RDKit✔️✔️:   double d2 = 0.0;
// RDKit✔️✔️:   for (unsigned int i = 0; i < dim; i++) {
// RDKit✔️✔️:     double d = end1Coords[i] - end2Coords[i];
// RDKit✔️✔️:     d2 += d * d;
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return d2;
// RDKit✔️✔️: }
// END RDKIT CPP LOCAL HELPER DistGeom::distance2
fn dist_violation_distance2(idx1: usize, idx2: usize, pos: &[f64], dim: usize) -> f64 {
    let mut d2 = 0.0;
    for i in 0..dim {
        let d = pos[dim * idx1 + i] - pos[dim * idx2 + i];
        d2 += d * d;
    }
    d2
}

// BEGIN RDKIT CPP LOCAL HELPER DistGeom::distance (DistViolationContribs.cpp:33-36)
// RDKit✔️✔️: inline double distance(const unsigned int idx1, const unsigned int idx2,
// RDKit✔️✔️:                        const double *pos, const unsigned int dim) {
// RDKit✔️✔️:   return sqrt(distance2(idx1, idx2, pos, dim));
// RDKit✔️✔️: }
// END RDKIT CPP LOCAL HELPER DistGeom::distance
fn dist_violation_distance(idx1: usize, idx2: usize, pos: &[f64], dim: usize) -> f64 {
    dist_violation_distance2(idx1, idx2, pos, dim).sqrt()
}

#[derive(Debug, Clone, PartialEq)]
pub struct DistViolationContribs {
    owner: *const ForceField,
    contribs: Vec<DistViolationContribsParams>,
}

impl Default for DistViolationContribs {
    fn default() -> Self {
        Self {
            owner: std::ptr::null(),
            contribs: Vec::new(),
        }
    }
}

impl DistViolationContribs {
    #[must_use]
    pub fn new(owner: &ForceField) -> Self {
        // BEGIN RDKIT CPP CONSTRUCTOR DistGeom::DistViolationContribs::DistViolationContribs (DistViolationContribs.cpp:16-19)
        // RDKit✔️✔️: DistViolationContribs::DistViolationContribs(ForceFields::ForceField *owner) {
        // RDKit✔️✔️:   PRECONDITION(owner, "bad owner");
        // RDKit✔️✔️:   dp_forceField = owner;
        // RDKit✔️✔️: }
        // END RDKIT CPP CONSTRUCTOR DistGeom::DistViolationContribs::DistViolationContribs
        Self {
            owner,
            contribs: Vec::new(),
        }
    }

    pub fn add_contrib(&mut self, idx1: usize, idx2: usize, ub: f64, lb: f64, weight: f64) {
        // BEGIN RDKIT CPP METHOD DistGeom::DistViolationContribs::addContrib (DistViolationContribs.h:49-52)
        // RDKit✔️✔️: void addContrib(unsigned int idx1, unsigned int idx2, double ub, double lb,
        // RDKit✔️✔️:                 double weight = 1.0) {
        // RDKit✔️✔️:   d_contribs.emplace_back(idx1, idx2, ub, lb, weight);
        // RDKit✔️✔️: }
        // END RDKIT CPP METHOD DistGeom::DistViolationContribs::addContrib
        self.contribs
            .push(DistViolationContribsParams::new(idx1, idx2, ub, lb, weight));
    }

    #[must_use]
    pub fn empty(&self) -> bool {
        // RDKit✔️✔️: bool empty() const { return d_contribs.empty(); }
        self.contribs.is_empty()
    }

    #[must_use]
    pub fn size(&self) -> usize {
        // RDKit✔️✔️: unsigned int size() const { return d_contribs.size(); }
        self.contribs.len()
    }

    #[must_use]
    pub fn contribs(&self) -> &[DistViolationContribsParams] {
        &self.contribs
    }

    fn owner(&self) -> &ForceField {
        assert!(!self.owner.is_null(), "no owner");
        // SAFETY: Matches the existing ForceFieldContrib owner-pointer convention.
        unsafe { &*self.owner }
    }
}

impl ForceFieldContrib for DistViolationContribs {
    fn copy(&self) -> Box<dyn ForceFieldContrib> {
        // RDKit✔️✔️: DistViolationContribs *copy() const override {
        // RDKit✔️✔️:   return new DistViolationContribs(*this);
        // RDKit✔️✔️: }
        Box::new(self.clone())
    }

    fn set_force_field(&mut self, owner: *const ForceField) {
        self.owner = owner;
    }

    fn get_energy(&self, pos: &[f64]) -> f64 {
        // BEGIN RDKIT CPP METHOD DistGeom::DistViolationContribs::getEnergy (DistViolationContribs.cpp:38-59)
        // RDKit✔️✔️: double DistViolationContribs::getEnergy(double *pos) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:   double accum = 0.0;
        // RDKit✔️✔️:   auto contrib = [&](const auto &c) {
        // RDKit✔️✔️:     double d2 = distance2(c.idx1, c.idx2, pos, dp_forceField->dimension());
        // RDKit✔️✔️:     double val = 0.0;
        // RDKit✔️✔️:     if (d2 > c.ub2) {
        // RDKit✔️✔️:       val = (d2 / (c.ub2)) - 1.0;
        // RDKit✔️✔️:     } else if (d2 < c.lb2) {
        // RDKit✔️✔️:       val = ((2 * c.lb2) / (c.lb2 + d2)) - 1.0;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     if (val > 0.0) {
        // RDKit✔️✔️:       accum += c.weight * val * val;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   };
        // RDKit✔️✔️:   for (const auto &c : d_contribs) {
        // RDKit✔️✔️:     contrib(c);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return accum;
        // RDKit✔️✔️: }
        // END RDKIT CPP METHOD DistGeom::DistViolationContribs::getEnergy
        assert!(!pos.is_empty(), "bad vector");
        let dim = self.owner().dimension();
        let mut accum = 0.0;
        for c in &self.contribs {
            let d2 = dist_violation_distance2(c.idx1, c.idx2, pos, dim);
            let mut val = 0.0;
            if d2 > c.ub2 {
                val = d2 / c.ub2 - 1.0;
            } else if d2 < c.lb2 {
                val = (2.0 * c.lb2) / (c.lb2 + d2) - 1.0;
            }
            if val > 0.0 {
                accum += c.weight * val * val;
            }
        }
        accum
    }

    fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
        // BEGIN RDKIT CPP METHOD DistGeom::DistViolationContribs::getGrad (DistViolationContribs.cpp:61-96)
        // RDKit✔️✔️: void DistViolationContribs::getGrad(double *pos, double *grad) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:   PRECONDITION(grad, "bad vector");
        // RDKit✔️✔️:   const unsigned int dim = this->dp_forceField->dimension();
        // RDKit✔️✔️:   auto contrib = [&](const auto &c) {
        // RDKit✔️✔️:     double d2 = distance2(c.idx1, c.idx2, pos, dp_forceField->dimension());
        // RDKit✔️✔️:     double d;
        // RDKit✔️✔️:     double preFactor = 0.0;
        // RDKit✔️✔️:     if (d2 > c.ub2) {
        // RDKit✔️✔️:       d = sqrt(d2);
        // RDKit✔️✔️:       preFactor = 4. * (((d * d) / c.ub2) - 1.0) * (d / c.ub2);
        // RDKit✔️✔️:     } else if (d2 < c.lb2) {
        // RDKit✔️✔️:       d = sqrt(d2);
        // RDKit✔️✔️:       double l2d2 = d2 + c.lb2;
        // RDKit✔️✔️:       preFactor = 8. * c.lb2 * d * (1. - 2 * c.lb2 / l2d2) / (l2d2 * l2d2);
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       return;
        // RDKit✔️✔️:     }
        // END RDKIT CPP METHOD DistGeom::DistViolationContribs::getGrad
        assert!(!pos.is_empty(), "bad vector");
        assert!(!grad.is_empty(), "bad vector");
        let dim = self.owner().dimension();
        for c in &self.contribs {
            let d2 = dist_violation_distance2(c.idx1, c.idx2, pos, dim);
            let d;
            let pre_factor;
            if d2 > c.ub2 {
                d = dist_violation_distance(c.idx1, c.idx2, pos, dim);
                pre_factor = 4.0 * ((d * d) / c.ub2 - 1.0) * (d / c.ub2);
            } else if d2 < c.lb2 {
                d = dist_violation_distance(c.idx1, c.idx2, pos, dim);
                let l2d2 = d2 + c.lb2;
                pre_factor = 8.0 * c.lb2 * d * (1.0 - 2.0 * c.lb2 / l2d2) / (l2d2 * l2d2);
            } else {
                continue;
            }
            // RDKit✔️✔️:     for (unsigned int i = 0; i < dim; i++) {
            // RDKit✔️✔️:       const auto p1 = dim * c.idx1 + i;
            // RDKit✔️✔️:       const auto p2 = dim * c.idx2 + i;
            // RDKit✔️✔️:       double dGrad;
            // RDKit✔️✔️:       if (d > 0.0) {
            // RDKit✔️✔️:         dGrad = c.weight * preFactor * (pos[p1] - pos[p2]) / d;
            // RDKit✔️✔️:       } else {
            // RDKit✔️✔️:         // FIX: this likely isn't right
            // RDKit✔️✔️:         dGrad = c.weight * preFactor * (pos[p1] - pos[p2]);
            // RDKit✔️✔️:       }
            // RDKit✔️✔️:       grad[p1] += dGrad;
            // RDKit✔️✔️:       grad[p2] -= dGrad;
            // RDKit✔️✔️:     }
            for i in 0..dim {
                let p1 = dim * c.idx1 + i;
                let p2 = dim * c.idx2 + i;
                let d_grad = if d > 0.0 {
                    c.weight * pre_factor * (pos[p1] - pos[p2]) / d
                } else {
                    c.weight * pre_factor * (pos[p1] - pos[p2])
                };
                grad[p1] += d_grad;
                grad[p2] -= d_grad;
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct FourthDimContribsParams {
    pub idx: usize,
    pub weight: f64,
}

impl FourthDimContribsParams {
    #[must_use]
    pub fn new(idx: usize, weight: f64) -> Self {
        // BEGIN RDKIT CPP CONSTRUCTOR DistGeom::FourthDimContribsParams (FourthDimContribs.h:19-23)
        // RDKit✔️✔️: struct FourthDimContribsParams {
        // RDKit✔️✔️:   unsigned int idx{0};
        // RDKit✔️✔️:   double weight{0.0};
        // RDKit✔️✔️:   FourthDimContribsParams(unsigned int idx, double w) : idx(idx), weight(w) {};
        // RDKit✔️✔️: };
        // END RDKIT CPP CONSTRUCTOR DistGeom::FourthDimContribsParams
        Self { idx, weight }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct FourthDimContribs {
    owner: *const ForceField,
    contribs: Vec<FourthDimContribsParams>,
}

impl Default for FourthDimContribs {
    fn default() -> Self {
        // RDKit✔️✔️: FourthDimContribs() = default;
        Self {
            owner: std::ptr::null(),
            contribs: Vec::new(),
        }
    }
}

impl FourthDimContribs {
    #[must_use]
    pub fn new(owner: &ForceField) -> Self {
        // BEGIN RDKIT CPP CONSTRUCTOR DistGeom::FourthDimContribs::FourthDimContribs (FourthDimContribs.h:36-40)
        // RDKit✔️✔️: FourthDimContribs(ForceFields::ForceField *owner) {
        // RDKit✔️✔️:   PRECONDITION(owner, "bad force field");
        // RDKit✔️✔️:   PRECONDITION(owner->dimension() == 4, "force field has wrong dimension");
        // RDKit✔️✔️:   dp_forceField = owner;
        // RDKit✔️✔️: }
        // END RDKIT CPP CONSTRUCTOR DistGeom::FourthDimContribs::FourthDimContribs
        assert_eq!(owner.dimension(), 4, "force field has wrong dimension");
        Self {
            owner,
            contribs: Vec::new(),
        }
    }

    pub fn add_contrib(&mut self, idx: usize, weight: f64) {
        // BEGIN RDKIT CPP METHOD DistGeom::FourthDimContribs::addContrib (FourthDimContribs.h:42-44)
        // RDKit✔️✔️: void addContrib(unsigned int idx, double weight) {
        // RDKit✔️✔️:   d_contribs.emplace_back(idx, weight);
        // RDKit✔️✔️: }
        // END RDKIT CPP METHOD DistGeom::FourthDimContribs::addContrib
        self.contribs
            .push(FourthDimContribsParams::new(idx, weight));
    }

    #[must_use]
    pub fn empty(&self) -> bool {
        // RDKit✔️✔️: bool empty() const { return d_contribs.empty(); }
        self.contribs.is_empty()
    }

    #[must_use]
    pub fn size(&self) -> usize {
        // RDKit✔️✔️: unsigned int size() const { return d_contribs.size(); }
        self.contribs.len()
    }

    #[must_use]
    pub fn contribs(&self) -> &[FourthDimContribsParams] {
        &self.contribs
    }

    fn owner(&self) -> &ForceField {
        assert!(!self.owner.is_null(), "no owner");
        // SAFETY: Matches the existing ForceFieldContrib owner-pointer convention.
        unsafe { &*self.owner }
    }
}

impl ForceFieldContrib for FourthDimContribs {
    fn copy(&self) -> Box<dyn ForceFieldContrib> {
        // RDKit✔️✔️: FourthDimContribs *copy() const override {
        // RDKit✔️✔️:   return new FourthDimContribs(*this);
        // RDKit✔️✔️: }
        Box::new(self.clone())
    }

    fn set_force_field(&mut self, owner: *const ForceField) {
        self.owner = owner;
    }

    fn get_energy(&self, pos: &[f64]) -> f64 {
        // BEGIN RDKIT CPP METHOD DistGeom::FourthDimContribs::getEnergy (FourthDimContribs.h:47-57)
        // RDKit✔️✔️: double getEnergy(double *pos) const override {
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:   constexpr unsigned int ffdim = 4;
        // RDKit✔️✔️:   double res = 0.0;
        // RDKit✔️✔️:   for (const auto &contrib : d_contribs) {
        // RDKit✔️✔️:     unsigned int pid = contrib.idx * ffdim + 3;
        // RDKit✔️✔️:     res += contrib.weight * pos[pid] * pos[pid];
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        // END RDKIT CPP METHOD DistGeom::FourthDimContribs::getEnergy
        assert!(!pos.is_empty(), "bad vector");
        let _owner = self.owner();
        let ffdim = 4;
        let mut res = 0.0;
        for contrib in &self.contribs {
            let pid = contrib.idx * ffdim + 3;
            res += contrib.weight * pos[pid] * pos[pid];
        }
        res
    }

    fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
        // BEGIN RDKIT CPP METHOD DistGeom::FourthDimContribs::getGrad (FourthDimContribs.h:61-70)
        // RDKit✔️✔️: void getGrad(double *pos, double *grad) const override {
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:   constexpr unsigned int ffdim = 4;
        // RDKit✔️✔️:   for (const auto &contrib : d_contribs) {
        // RDKit✔️✔️:     unsigned int pid = contrib.idx * ffdim + 3;
        // RDKit✔️✔️:     grad[pid] += contrib.weight * pos[pid];
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // END RDKIT CPP METHOD DistGeom::FourthDimContribs::getGrad
        assert!(!pos.is_empty(), "bad vector");
        assert!(!grad.is_empty(), "bad vector");
        let _owner = self.owner();
        let ffdim = 4;
        for contrib in &self.contribs {
            let pid = contrib.idx * ffdim + 3;
            grad[pid] += contrib.weight * pos[pid];
        }
    }
}

// ──────────────────────────────────────────────
// VDW radii (Å) — RDKit PeriodicTable::getRvdw
// ──────────────────────────────────────────────

/// RDKit❗✔️: VDW radii from RDKit PeriodicTable::getRvdw
fn vdw_radius(atomic_num: u8) -> f64 {
    match atomic_num {
        0 => 0.0,   // *
        1 => 1.2,   // H
        2 => 1.4,   // He
        3 => 2.2,   // Li
        4 => 1.9,   // Be
        5 => 1.8,   // B
        6 => 1.7,   // C
        7 => 1.6,   // N
        8 => 1.55,  // O
        9 => 1.5,   // F
        10 => 1.54, // Ne
        11 => 2.4,  // Na
        12 => 2.2,  // Mg
        13 => 2.1,  // Al
        14 => 2.1,  // Si
        15 => 1.95, // P
        16 => 1.8,  // S
        17 => 1.8,  // Cl
        18 => 1.88, // Ar
        19 => 2.8,  // K
        20 => 2.4,  // Ca
        21 => 2.3,  // Sc
        22 => 2.15, // Ti
        23 => 2.05, // V
        24 => 2.05, // Cr
        25 => 2.05, // Mn
        26 => 2.05, // Fe
        27 => 2.0,  // Co
        28 => 2.0,  // Ni
        29 => 2.0,  // Cu
        30 => 2.1,  // Zn
        31 => 2.1,  // Ga
        32 => 2.1,  // Ge
        33 => 2.05, // As
        34 => 1.9,  // Se
        35 => 1.9,  // Br
        36 => 2.02, // Kr
        37 => 2.9,  // Rb
        38 => 2.55, // Sr
        39 => 2.4,  // Y
        40 => 2.3,  // Zr
        41 => 2.15, // Nb
        42 => 2.1,  // Mo
        43 => 2.05, // Tc
        44 => 2.05, // Ru
        45 => 2.0,  // Rh
        46 => 2.05, // Pd
        47 => 2.1,  // Ag
        48 => 2.2,  // Cd
        49 => 2.2,  // In
        50 => 2.25, // Sn
        51 => 2.2,  // Sb
        52 => 2.1,  // Te
        53 => 2.1,  // I
        54 => 2.16, // Xe
        55 => 3.0,  // Cs
        56 => 2.7,  // Ba
        57 => 2.5,  // La
        58 => 2.48, // Ce
        59 => 2.47, // Pr
        60 => 2.45, // Nd
        61 => 2.43, // Pm
        62 => 2.42, // Sm
        63 => 2.4,  // Eu
        64 => 2.38, // Gd
        65 => 2.37, // Tb
        66 => 2.35, // Dy
        67 => 2.33, // Ho
        68 => 2.32, // Er
        69 => 2.3,  // Tm
        70 => 2.28, // Yb
        71 => 2.27, // Lu
        72 => 2.25, // Hf
        73 => 2.2,  // Ta
        74 => 2.1,  // W
        75 => 2.05, // Re
        76 => 2.0,  // Os
        77 => 2.0,  // Ir
        78 => 2.05, // Pt
        79 => 2.1,  // Au
        80 => 2.05, // Hg
        81 => 2.2,  // Tl
        82 => 2.3,  // Pb
        83 => 2.3,  // Bi
        84 => 2.0,  // Po
        85 => 2.0,  // At
        _ => 2.0,   // default
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
fn is_hbond_donor(atom: &Atom, total_num_hydrogens: u32) -> bool {
    is_hbond_acceptor(atom.atomic_number()) && total_num_hydrogens > 0
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

fn rdkit_default_valence(atomic_num: u8) -> Option<i32> {
    let vals = rdkit_valence_list(atomic_num).ok()??;
    vals.iter().copied().find(|v| *v >= 0)
}

fn rdkit_n_outer_electrons(atomic_num: u8) -> Option<i32> {
    crate::chemistry::valence::periodic_table_outer_electrons(atomic_num).ok()
}

fn bond_valence_contrib_for_atom(bond: &Bond, atom_index: usize) -> f64 {
    if bond.begin().index() != atom_index && bond.end().index() != atom_index {
        return 0.0;
    }
    match bond.order() {
        BondOrder::Null | BondOrder::Zero | BondOrder::Ionic => 0.0,
        BondOrder::Single => 1.0,
        BondOrder::Double => 2.0,
        BondOrder::Triple => 3.0,
        BondOrder::Quadruple => 4.0,
        BondOrder::Quintuple => 5.0,
        BondOrder::Hextuple => 6.0,
        BondOrder::OneAndHalf => 1.5,
        BondOrder::TwoAndHalf => 2.5,
        BondOrder::ThreeAndHalf => 3.5,
        BondOrder::FourAndHalf => 4.5,
        BondOrder::FiveAndHalf => 5.5,
        BondOrder::Aromatic => 1.5,
        BondOrder::Hydrogen => 0.0,
        BondOrder::Dative
        | BondOrder::DativeOne
        | BondOrder::DativeLeft
        | BondOrder::DativeRight => {
            if bond.end().index() == atom_index {
                1.0
            } else {
                0.0
            }
        }
        BondOrder::ThreeCenter | BondOrder::Unspecified | BondOrder::Other => 0.0,
    }
}

fn count_atom_electrons_rdkit(
    mol: &Molecule,
    assignment: &crate::ValenceAssignment,
    atom_degree: &[usize],
    atom_index: usize,
) -> i32 {
    let atom = &mol.atoms()[atom_index];
    let Some(dv) = rdkit_default_valence(atom.atomic_number()) else {
        return -1;
    };
    if dv <= 1 {
        return -1;
    }
    let mut degree = atom_degree[atom_index] as i32
        + atom.explicit_hydrogens() as i32
        + assignment.implicit_hydrogens[atom_index] as i32;
    for bond in mol.bonds() {
        if (bond.begin().index() == atom_index || bond.end().index() == atom_index)
            && bond_valence_contrib_for_atom(bond, atom_index) == 0.0
        {
            degree -= 1;
        }
    }
    if degree > 3 {
        return -1;
    }
    let Some(nouter) = rdkit_n_outer_electrons(atom.atomic_number()) else {
        return -1;
    };
    let nlp = (nouter - dv - atom.formal_charge() as i32).max(0);
    let radicals = atom.radical_electrons() as i32;
    let mut res = (dv - degree) + nlp - radicals;
    if res > 1 {
        let n_unsaturations =
            assignment.explicit_valence[atom_index] - atom_degree[atom_index] as i32;
        if n_unsaturations > 1 {
            res = 1;
        }
    }
    res
}

fn is_atom_conjugation_candidate(
    mol: &Molecule,
    assignment: &crate::ValenceAssignment,
    atom_degree: &[usize],
    atom_index: usize,
) -> bool {
    let atom = &mol.atoms()[atom_index];
    if let Ok(Some(vals)) = rdkit_valence_list(atom.atomic_number())
        && atom.formal_charge() == 0
        && !vals.is_empty()
        && vals[0] >= 0
    {
        let total_valence =
            assignment.explicit_valence[atom_index] + assignment.implicit_hydrogens[atom_index];
        if total_valence > vals[0] {
            return false;
        }
    }
    let nouter = rdkit_n_outer_electrons(atom.atomic_number()).unwrap_or(0);
    let total_degree = atom_degree[atom_index]
        + atom.explicit_hydrogens() as usize
        + assignment.implicit_hydrogens[atom_index].max(0) as usize;
    ((atom.atomic_number() <= 10)
        || (nouter != 5 && nouter != 6)
        || (nouter == 6 && total_degree < 2))
        && count_atom_electrons_rdkit(mol, assignment, atom_degree, atom_index) > 0
}

fn compute_conjugated_bonds_for_uff(
    mol: &Molecule,
    assignment: &crate::ValenceAssignment,
    atom_degree: &[usize],
) -> Vec<bool> {
    let mut conjugated = vec![false; mol.bonds().len()];
    for (bi, bond) in mol.bonds().iter().enumerate() {
        conjugated[bi] = bond.is_aromatic();
    }
    for atom_index in 0..mol.atoms().len() {
        if !is_atom_conjugation_candidate(mol, assignment, atom_degree, atom_index) {
            continue;
        }
        let sbo = atom_degree[atom_index]
            + mol.atoms()[atom_index].explicit_hydrogens() as usize
            + assignment.implicit_hydrogens[atom_index].max(0) as usize;
        if !(2..=3).contains(&sbo) {
            continue;
        }
        let incident_bonds: Vec<usize> = mol
            .bonds()
            .iter()
            .enumerate()
            .filter_map(|(bond_index, bond)| {
                if bond.begin().index() == atom_index || bond.end().index() == atom_index {
                    Some(bond_index)
                } else {
                    None
                }
            })
            .collect();
        for &bond1_idx in &incident_bonds {
            let bond1 = &mol.bonds()[bond1_idx];
            if bond_valence_contrib_for_atom(bond1, atom_index) < 1.5 {
                continue;
            }
            let other1 = if bond1.begin().index() == atom_index {
                bond1.end().index()
            } else {
                bond1.begin().index()
            };
            if !is_atom_conjugation_candidate(mol, assignment, atom_degree, other1) {
                continue;
            }
            for &bond2_idx in &incident_bonds {
                if bond1_idx == bond2_idx {
                    continue;
                }
                let bond2 = &mol.bonds()[bond2_idx];
                let other2 = if bond2.begin().index() == atom_index {
                    bond2.end().index()
                } else {
                    bond2.begin().index()
                };
                let sbo2 = atom_degree[other2]
                    + mol.atoms()[other2].explicit_hydrogens() as usize
                    + assignment.implicit_hydrogens[other2].max(0) as usize;
                if sbo2 > 3 {
                    continue;
                }
                if is_atom_conjugation_candidate(mol, assignment, atom_degree, other2) {
                    conjugated[bond1_idx] = true;
                    conjugated[bond2_idx] = true;
                }
            }
        }
    }
    conjugated
}

fn compute_hybridizations_for_uff(
    mol: &Molecule,
    assignment: &crate::ValenceAssignment,
    atom_degree: &[usize],
    atom_has_conjugated_bond: &[bool],
) -> Vec<Hybridization> {
    let mut out = Vec::with_capacity(mol.atoms().len());
    for (atom_index, atom) in mol.atoms().iter().enumerate() {
        if atom.atomic_number() == 0 {
            out.push(Hybridization::Unspecified);
            continue;
        }
        let mut deg = atom_degree[atom_index] as i32
            + atom.explicit_hydrogens() as i32
            + assignment.implicit_hydrogens[atom_index];
        for bond in mol.bonds() {
            if (bond.begin().index() == atom_index || bond.end().index() == atom_index)
                && matches!(bond.order(), BondOrder::Dative | BondOrder::DativeOne)
                && bond.end().index() != atom_index
            {
                deg -= 1;
            }
        }
        // BEGIN RDKIT CPP FUNCTION MolOps::setHybridization chiral-coordination branch
        // RDKit✔️✔️: switch (atom->getChiralTag()) {
        let chiral_hybridization = match atom.chiral_tag() {
            // RDKit✔️✔️:   case Atom::ChiralType::CHI_TETRAHEDRAL:
            // RDKit✔️✔️:   case Atom::ChiralType::CHI_TETRAHEDRAL_CW:
            // RDKit✔️✔️:   case Atom::ChiralType::CHI_TETRAHEDRAL_CCW:
            ChiralTag::Tetrahedral | ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw => {
                // RDKit✔️✔️:     if (atom->getTotalDegree() == 4) {
                // RDKit✔️✔️:       atom->setHybridization(Atom::HybridizationType::SP3);
                // RDKit✔️✔️:       continue;
                // RDKit✔️✔️:     }
                (deg == 4).then_some(Hybridization::Sp3)
                // RDKit✔️✔️:     break;
            }
            // RDKit✔️✔️:   case Atom::ChiralType::CHI_SQUAREPLANAR:
            ChiralTag::SquarePlanar => {
                // RDKit✔️✔️:     if (atom->getTotalDegree() <= 4 && atom->getTotalDegree() >= 2) {
                // RDKit✔️✔️:       atom->setHybridization(Atom::HybridizationType::SP2D);
                // RDKit✔️✔️:       continue;
                // RDKit✔️✔️:     }
                (2..=4).contains(&deg).then_some(Hybridization::Sp2d)
                // RDKit✔️✔️:     break;
            }
            // RDKit✔️✔️:   case Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL:
            ChiralTag::TrigonalBipyramidal => {
                // RDKit✔️✔️:     if (atom->getTotalDegree() <= 5 && atom->getTotalDegree() >= 2) {
                // RDKit✔️✔️:       atom->setHybridization(Atom::HybridizationType::SP3D);
                // RDKit✔️✔️:       continue;
                // RDKit✔️✔️:     }
                (2..=5).contains(&deg).then_some(Hybridization::Sp3d)
                // RDKit✔️✔️:     break;
            }
            // RDKit✔️✔️:   case Atom::ChiralType::CHI_OCTAHEDRAL:
            ChiralTag::Octahedral => {
                // RDKit✔️✔️:     if (atom->getTotalDegree() <= 6 && atom->getTotalDegree() >= 2) {
                // RDKit✔️✔️:       atom->setHybridization(Atom::HybridizationType::SP3D2);
                // RDKit✔️✔️:       continue;
                // RDKit✔️✔️:     }
                (2..=6).contains(&deg).then_some(Hybridization::Sp3d2)
                // RDKit✔️✔️:     break;
            }
            // RDKit✔️✔️:   default:
            // RDKit✔️✔️:     break;
            _ => None,
            // RDKit✔️✔️: }
        };
        // END RDKIT CPP FUNCTION MolOps::setHybridization chiral-coordination branch
        if let Some(hybridization) = chiral_hybridization {
            out.push(hybridization);
            continue;
        }
        let hyb = if atom.atomic_number() <= 1 {
            match deg {
                0 | 1 => Hybridization::S,
                2 => Hybridization::Sp,
                3 => Hybridization::Sp2,
                4 => Hybridization::Sp3,
                5 => Hybridization::Sp3d,
                6 => Hybridization::Sp3d2,
                _ => Hybridization::Unspecified,
            }
        } else {
            let nouter = rdkit_n_outer_electrons(atom.atomic_number()).unwrap_or(0);
            let total_valence =
                assignment.explicit_valence[atom_index] + assignment.implicit_hydrogens[atom_index];
            let num_free = nouter - (total_valence + atom.formal_charge() as i32);
            let norbs = if total_valence + nouter - (atom.formal_charge() as i32) < 8 {
                let radicals = atom.radical_electrons() as i32;
                let lone_pairs = (num_free - radicals) / 2;
                deg + lone_pairs + radicals
            } else {
                let lone_pairs = num_free / 2;
                deg + lone_pairs
            };
            match norbs {
                0 | 1 => Hybridization::S,
                2 => Hybridization::Sp,
                3 => Hybridization::Sp2,
                4 => {
                    let total_degree = atom_degree[atom_index]
                        + atom.explicit_hydrogens() as usize
                        + assignment.implicit_hydrogens[atom_index].max(0) as usize;
                    if total_degree > 3 || !atom_has_conjugated_bond[atom_index] {
                        Hybridization::Sp3
                    } else {
                        Hybridization::Sp2
                    }
                }
                5 => Hybridization::Sp3d,
                6 => Hybridization::Sp3d2,
                _ => Hybridization::Unspecified,
            }
        };
        out.push(hyb);
    }
    out
}

fn atom_total_valence_for_uff(assignment: &crate::ValenceAssignment, atom_index: usize) -> i32 {
    // BEGIN RDKIT CPP FUNCTION Atom::getTotalValence (Atom.cpp:334-336)
    // RDKit✔️✔️: unsigned int Atom::getTotalValence() const {
    // RDKit✔️✔️:   return getValence(ValenceType::EXPLICIT) +
    // RDKit✔️✔️:          getValence(ValenceType::IMPLICIT);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Atom::getTotalValence
    assignment.explicit_valence[atom_index] + assignment.implicit_hydrogens[atom_index]
}

// BEGIN RDKIT CPP FUNCTION Atom::getTotalNumHs (Atom.cpp:283-293)
// RDKit✔️✔️: unsigned int Atom::getTotalNumHs(bool includeNeighbors) const {
// RDKit✔️✔️:   int res = getNumExplicitHs() + getNumImplicitHs();
// RDKit✔️✔️:   if (includeNeighbors && dp_mol) {
// RDKit✔️✔️:     auto nbrs = dp_mol->atomNeighbors(this);
// RDKit✔️✔️:     res += std::count_if(nbrs.begin(), nbrs.end(), [](const auto nbr) {
// RDKit✔️✔️:       return (nbr->getAtomicNum() == 1);
// RDKit✔️✔️:     });
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return res;
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION Atom::getTotalNumHs
fn total_num_hydrogens_for_distgeom(
    mol: &Molecule,
    atom: &Atom,
    assignment: &crate::ValenceAssignment,
    atom_index: usize,
) -> u32 {
    let mut res =
        atom.explicit_hydrogens() as u32 + assignment.implicit_hydrogens[atom_index].max(0) as u32;
    res += neighbors_for_atom(mol, atom_index)
        .into_iter()
        .filter(|&nbr| mol.atoms()[nbr].atomic_number() == 1)
        .count() as u32;
    res
}

#[cfg(test)]
fn get_atom_label_for_uff(
    mol: &Molecule,
    assignment: &crate::ValenceAssignment,
    hybridizations: &[Hybridization],
    atom_has_conjugated_bond: &[bool],
    atom_index: usize,
) -> Result<String, DgBoundsError> {
    let atom = &mol.atoms()[atom_index];
    let total_valence = atom_total_valence_for_uff(assignment, atom_index);
    source_get_atom_label_for_uff(
        atom,
        atom_index,
        total_valence,
        hybridizations[atom_index],
        atom_has_conjugated_bond[atom_index],
        true,
    )
    .map_err(|err| {
        DgBoundsError::GenerationFailed(format!(
            "UFF atom symbol lookup failed for atomic number {}: {err}",
            atom.atomic_number()
        ))
    })
}

fn calc_bond_rest_length(bond_order: f64, end1: &AtomicParams, end2: &AtomicParams) -> f64 {
    let ri = end1.r1;
    let rj = end2.r1;
    let r_bo = -UFF_LAMBDA * (ri + rj) * bond_order.ln();
    let xi = end1.gmp_xi;
    let xj = end2.gmp_xi;
    let root_delta = xi.sqrt() - xj.sqrt();
    let r_en = ri * rj * root_delta * root_delta / (xi * ri + xj * rj);
    ri + rj + r_bo - r_en
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
// Topological distance matrix
// ──────────────────────────────────────────────

use crate::chemistry::matrices::{
    LOCAL_INF_DISTANCE as LOCAL_INF_DIST,
    topological_distance_matrix as compute_topological_distances,
};

fn flatten_topological_distances_matrix(mol: &Molecule) -> Vec<f64> {
    compute_topological_distances(mol).to_vec()
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

#[derive(Clone)]
struct BoundsMatrix {
    data: Vec<Vec<f64>>,
    n: usize,
}

impl BoundsMatrix {
    fn new(n: usize) -> Self {
        let mut matrix = Self {
            data: vec![vec![0.0; n]; n],
            n,
        };
        init_bounds_mat_shared(&mut matrix, DEFAULT_LOWER, DEFAULT_UPPER);
        matrix
    }

    fn get_val(&self, i: usize, j: usize) -> f64 {
        self.data[i][j]
    }

    fn set_val(&mut self, i: usize, j: usize, value: f64) {
        self.data[i][j] = value;
    }

    // BEGIN RDKIT CPP FUNCTION DistGeom::BoundsMatrix::getUpperBound (BoundsMatrix.h:37-43)
    // RDKit✔️✔️: inline double getUpperBound(unsigned int i, unsigned int j) const {
    // RDKit✔️✔️:   if (i < j) {
    // RDKit✔️✔️:     return getVal(i, j);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return getVal(j, i);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DistGeom::BoundsMatrix::getUpperBound
    fn get_upper(&self, i: usize, j: usize) -> f64 {
        if i < j {
            self.get_val(i, j)
        } else {
            self.get_val(j, i)
        }
    }

    // BEGIN RDKIT CPP FUNCTION DistGeom::BoundsMatrix::setUpperBound (BoundsMatrix.h:46-53)
    // RDKit✔️✔️: inline void setUpperBound(unsigned int i, unsigned int j, double val) {
    // RDKit✔️✔️:   CHECK_INVARIANT(val >= 0.0, "Negative upper bound");
    // RDKit✔️✔️:   if (i < j) {
    // RDKit✔️✔️:     setVal(i, j, val);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     setVal(j, i, val);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DistGeom::BoundsMatrix::setUpperBound
    fn invalid_bounds(
        reason: &'static str,
        i: usize,
        j: usize,
        lb: f64,
        ub: f64,
        current_lower: f64,
        current_upper: f64,
    ) -> DgBoundsError {
        DgBoundsError::InvalidBounds(format!(
            "{reason} for atom pair ({i}, {j}); requested lower={lb:.6}, upper={ub:.6}, current lower={current_lower:.6}, current upper={current_upper:.6}"
        ))
    }

    fn set_upper_unchecked(&mut self, i: usize, j: usize, v: f64) {
        if i < j {
            self.set_val(i, j, v);
        } else {
            self.set_val(j, i, v);
        }
    }

    fn set_upper(&mut self, i: usize, j: usize, v: f64) -> Result<(), DgBoundsError> {
        if v < 0.0 || v.is_nan() {
            return Err(Self::invalid_bounds(
                "negative upper bound",
                i,
                j,
                self.get_lower(i, j),
                v,
                self.get_lower(i, j),
                self.get_upper(i, j),
            ));
        }
        self.set_upper_unchecked(i, j, v);
        Ok(())
    }

    // BEGIN RDKIT CPP FUNCTION DistGeom::BoundsMatrix::setUpperBoundIfBetter (BoundsMatrix.h:57-62)
    // RDKit✔️✔️: inline void setUpperBoundIfBetter(unsigned int i, unsigned int j,
    // RDKit✔️✔️:                                     double val) {
    // RDKit✔️✔️:   if ((val < getUpperBound(i, j)) && (val > getLowerBound(i, j))) {
    // RDKit✔️✔️:     setUpperBound(i, j, val);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DistGeom::BoundsMatrix::setUpperBoundIfBetter
    fn set_upper_if_better(&mut self, i: usize, j: usize, val: f64) -> Result<(), DgBoundsError> {
        if val < self.get_upper(i, j) && val > self.get_lower(i, j) {
            self.set_upper(i, j, val)?;
        }
        Ok(())
    }

    // BEGIN RDKIT CPP FUNCTION DistGeom::BoundsMatrix::setLowerBound (BoundsMatrix.h:65-72)
    // RDKit✔️✔️: inline void setLowerBound(unsigned int i, unsigned int j, double val) {
    // RDKit✔️✔️:   CHECK_INVARIANT(val >= 0.0, "Negative lower bound");
    // RDKit✔️✔️:   if (i < j) {
    // RDKit✔️✔️:     setVal(j, i, val);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     setVal(i, j, val);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DistGeom::BoundsMatrix::setLowerBound
    fn set_lower_unchecked(&mut self, i: usize, j: usize, v: f64) {
        if i < j {
            self.set_val(j, i, v);
        } else {
            self.set_val(i, j, v);
        }
    }

    fn set_lower(&mut self, i: usize, j: usize, v: f64) -> Result<(), DgBoundsError> {
        if v < 0.0 || v.is_nan() {
            return Err(Self::invalid_bounds(
                "negative lower bound",
                i,
                j,
                v,
                self.get_upper(i, j),
                self.get_lower(i, j),
                self.get_upper(i, j),
            ));
        }
        self.set_lower_unchecked(i, j, v);
        Ok(())
    }

    // BEGIN RDKIT CPP FUNCTION DistGeom::BoundsMatrix::setLowerBoundIfBetter (BoundsMatrix.h:76-81)
    // RDKit✔️✔️: inline void setLowerBoundIfBetter(unsigned int i, unsigned int j,
    // RDKit✔️✔️:                                     double val) {
    // RDKit✔️✔️:   if ((val > getLowerBound(i, j)) && (val < getUpperBound(i, j))) {
    // RDKit✔️✔️:     setLowerBound(i, j, val);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DistGeom::BoundsMatrix::setLowerBoundIfBetter
    fn set_lower_if_better(&mut self, i: usize, j: usize, val: f64) -> Result<(), DgBoundsError> {
        if val > self.get_lower(i, j) && val < self.get_upper(i, j) {
            self.set_lower(i, j, val)?;
        }
        Ok(())
    }

    // BEGIN RDKIT CPP FUNCTION DistGeom::BoundsMatrix::getLowerBound (BoundsMatrix.h:84-90)
    // RDKit✔️✔️: inline double getLowerBound(unsigned int i, unsigned int j) const {
    // RDKit✔️✔️:   if (i < j) {
    // RDKit✔️✔️:     return getVal(j, i);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return getVal(i, j);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DistGeom::BoundsMatrix::getLowerBound
    fn get_lower(&self, i: usize, j: usize) -> f64 {
        if i < j {
            self.get_val(j, i)
        } else {
            self.get_val(i, j)
        }
    }

    // BEGIN RDKIT CPP FUNCTION DistGeom::BoundsMatrix::checkValid (BoundsMatrix.h:94-103)
    // RDKit✔️✔️: inline bool checkValid() const {
    // RDKit✔️✔️:   unsigned int i, j;
    // RDKit✔️✔️:   for (i = 1; i < d_nRows; i++) {
    // RDKit✔️✔️:     for (j = 0; j < i; j++) {
    // RDKit✔️✔️:       if (getUpperBound(i, j) < getLowerBound(i, j)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DistGeom::BoundsMatrix::checkValid
    fn check_valid(&self) -> bool {
        for i in 1..self.n {
            for j in 0..i {
                if self.get_upper(i, j) < self.get_lower(i, j) {
                    return false;
                }
            }
        }
        true
    }

    fn num_rows(&self) -> usize {
        self.n
    }
}

// ──────────────────────────────────────────────
// Distance-matrix randomization helpers
// ──────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq)]
struct SymmMatrix {
    data: Vec<f64>,
    n: usize,
}

impl SymmMatrix {
    fn new(n: usize) -> Self {
        Self {
            data: vec![0.0; n * (n + 1) / 2],
            n,
        }
    }

    fn with_value(n: usize, value: f64) -> Self {
        // BEGIN RDKIT CPP CONSTRUCTOR RDNumeric::SymmMatrix::SymmMatrix value overload (SymmMatrix.h:40-48)
        // RDKit✔️✔️: SymmMatrix(unsigned int N, TYPE val)
        // RDKit✔️✔️:     : d_size(N), d_dataSize(N * (N + 1) / 2) {
        // RDKit✔️✔️:   TYPE *data = new TYPE[d_dataSize];
        // RDKit✔️✔️:   unsigned int i;
        // RDKit✔️✔️:   for (i = 0; i < d_dataSize; i++) {
        // RDKit✔️✔️:     data[i] = val;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   d_data.reset(data);
        // RDKit✔️✔️: }
        // END RDKIT CPP CONSTRUCTOR RDNumeric::SymmMatrix::SymmMatrix value overload
        Self {
            data: vec![value; n * (n + 1) / 2],
            n,
        }
    }

    fn num_rows(&self) -> usize {
        self.n
    }

    fn get_data_size(&self) -> usize {
        self.data.len()
    }

    fn get_data(&self) -> &[f64] {
        &self.data
    }

    fn get_val(&self, i: usize, j: usize) -> f64 {
        self.data[SymmMatrix::data_index(i, j)]
    }

    fn set_val(&mut self, i: usize, j: usize, value: f64) {
        let idx = SymmMatrix::data_index(i, j);
        self.data[idx] = value;
    }

    fn data_index(i: usize, j: usize) -> usize {
        let (row, col) = if i >= j { (i, j) } else { (j, i) };
        row * (row + 1) / 2 + col
    }
}

#[derive(Debug, Clone, PartialEq)]
struct DoubleMatrix {
    data: Vec<f64>,
    n_rows: usize,
    n_cols: usize,
}

impl DoubleMatrix {
    fn new(n_rows: usize, n_cols: usize) -> Self {
        // BEGIN RDKIT CPP CONSTRUCTOR RDNumeric::Matrix::Matrix size overload (Matrix.h:33-38)
        // RDKit✔️✔️: Matrix(unsigned int nRows, unsigned int nCols)
        // RDKit✔️✔️:     : d_nRows(nRows), d_nCols(nCols), d_dataSize(nRows * nCols) {
        // RDKit✔️✔️:   TYPE *data = new TYPE[d_dataSize];
        // RDKit✔️✔️:   memset(static_cast<void *>(data), 0, d_dataSize * sizeof(TYPE));
        // RDKit✔️✔️:   d_data.reset(data);
        // RDKit✔️✔️: }
        // END RDKIT CPP CONSTRUCTOR RDNumeric::Matrix::Matrix size overload
        Self {
            data: vec![0.0; n_rows * n_cols],
            n_rows,
            n_cols,
        }
    }

    fn num_rows(&self) -> usize {
        self.n_rows
    }

    fn num_cols(&self) -> usize {
        self.n_cols
    }

    fn get_val(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.n_cols + j]
    }
}

trait RdkitDoubleRng {
    fn next_unit_f64(&mut self) -> f64;
}

#[derive(Debug, Clone)]
struct RdkitDistgeomMinStdRand {
    state: u64,
}

impl Default for RdkitDistgeomMinStdRand {
    fn default() -> Self {
        // BEGIN RDKIT CPP STATIC GENERATOR RDKit::generator (RDGeneral/utils.cpp:18)
        // RDKit✔️✔️: static rng_type generator(42u);
        // END RDKIT CPP STATIC GENERATOR RDKit::generator
        Self { state: 42 }
    }
}

impl RdkitDistgeomMinStdRand {
    fn new(seed: i32) -> Self {
        // BEGIN RDKIT CPP FUNCTION RDKit::getRandomGenerator seed behavior (RDGeneral/utils.cpp:24-29)
        // RDKit✔️✔️: rng_type &getRandomGenerator(int seed) {
        // RDKit✔️✔️:   if (seed > 0) {
        // RDKit✔️✔️:     generator.seed(seed);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return generator;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION RDKit::getRandomGenerator seed behavior
        if seed > 0 {
            // BEGIN BOOST CPP METHOD boost::random::linear_congruential_engine::seed (linear_congruential.hpp)
            // RDKit✔️✔️: if(modulus == 0) {
            // RDKit✔️✔️:     _x = x0;
            // RDKit✔️✔️: } else {
            // RDKit✔️✔️:     _x = x0 % modulus;
            // RDKit✔️✔️: }
            // RDKit✔️✔️: if(increment == 0 && _x == 0) {
            // RDKit✔️✔️:     _x = 1;
            // RDKit✔️✔️: }
            // END BOOST CPP METHOD boost::random::linear_congruential_engine::seed
            let mut state = seed as u64 % RDKIT_RANDOM_MODULUS;
            if state == 0 {
                state = 1;
            }
            Self { state }
        } else {
            Self::default()
        }
    }

    fn new_from_embed_points_seed(seed: i32) -> Self {
        assert!(seed >= 0);
        let mut state = seed as u64 % RDKIT_RANDOM_MODULUS;
        if state == 0 {
            state = 1;
        }
        Self { state }
    }

    fn next_raw(&mut self) -> u32 {
        self.state = (self.state * RDKIT_RANDOM_MULTIPLIER) % RDKIT_RANDOM_MODULUS;
        self.state as u32
    }
}

#[cfg(not(target_arch = "wasm32"))]
unsafe extern "C" {
    fn clock() -> libc::clock_t;
}

#[cfg(not(target_arch = "wasm32"))]
fn rdkit_clock_seed() -> Result<i32, DgBoundsError> {
    // BEGIN RDKIT CPP CLOCK SEED SOURCE C library clock() use in conformer generation
    // RDKit✔️✔️: clock()
    // END RDKIT CPP CLOCK SEED SOURCE
    Ok(unsafe { clock() as i32 })
}

#[cfg(target_arch = "wasm32")]
fn rdkit_clock_seed() -> Result<i32, DgBoundsError> {
    // BEGIN RDKIT CPP CLOCK SEED SOURCE C library clock() use in conformer generation
    // RDKit❌✔️: clock()
    // END RDKIT CPP CLOCK SEED SOURCE
    Err(DgBoundsError::WasmImplicitClockSeedUnsupported)
}

impl RdkitDoubleRng for RdkitDistgeomMinStdRand {
    fn next_unit_f64(&mut self) -> f64 {
        // BEGIN BOOST CPP FUNCTION boost::random::detail::generate_uniform_real integral path (uniform_real_distribution.hpp)
        // RDKit✔️✔️: typedef typename Engine::result_type base_result;
        // RDKit✔️✔️: result_type numerator = static_cast<T>(subtract<base_result>()(eng(), (eng.min)()));
        // RDKit✔️✔️: result_type divisor = static_cast<T>(subtract<base_result>()((eng.max)(), (eng.min)())) + 1;
        // RDKit✔️✔️: T result = numerator / divisor * (max_value - min_value) + min_value;
        // RDKit✔️✔️: if(result < max_value) return result;
        // END BOOST CPP FUNCTION boost::random::detail::generate_uniform_real integral path
        let raw = self.next_raw() as f64;
        (raw - 1.0) / (RDKIT_RANDOM_MODULUS as f64 - 1.0)
    }
}

fn rdkit_embedder_multiplication_overflows(a: i32, b: i32) -> bool {
    // BEGIN RDKIT CPP TEMPLATE FUNCTION DGeomHelpers::detail::multiplication_overflows_ (Embedder.cpp:1346-1353)
    // RDKit✔️✔️: template <class T>
    // RDKit✔️✔️: bool multiplication_overflows_(T a, T b) {
    // RDKit✔️✔️:   // a * b > c if and only if a > c / b
    // RDKit✔️✔️:   if (a == 0 || b == 0) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return a > std::numeric_limits<T>::max() / b;
    // RDKit✔️✔️: }
    // END RDKIT CPP TEMPLATE FUNCTION DGeomHelpers::detail::multiplication_overflows_
    if a == 0 || b == 0 {
        return false;
    }
    a > i32::MAX / b
}

fn rdkit_embedder_conformer_seed(
    random_seed: i32,
    conformer_index: usize,
    enable_sequential_random_seeds: bool,
) -> i32 {
    // BEGIN RDKIT CPP BLOCK DGeomHelpers::detail::embedHelper_ per-conformer seed policy (Embedder.cpp:1393-1426)
    // RDKit✔️✔️:     CHECK_INVARIANT(
    // RDKit✔️✔️:         params->randomSeed >= -1,
    // RDKit✔️✔️:         "random seed must either be positive, zero, or negative one");
    // RDKit✔️✔️:     int new_seed = params->randomSeed;
    // RDKit✔️✔️:     if (new_seed > -1) {
    // RDKit✔️✔️:       if (params->enableSequentialRandomSeeds) {
    // RDKit✔️✔️:         new_seed += ci + 1;
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         if (!multiplication_overflows_(rdcast<int>(ci + 1),
    // RDKit✔️✔️:                                        params->randomSeed)) {
    // RDKit✔️✔️:           // old method of computing a new seed
    // RDKit✔️✔️:           new_seed = (ci + 1) * params->randomSeed;
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           // If the above simple multiplication will overflow, use a
    // RDKit✔️✔️:           // cheap and easy way to hash the conformer index and seed
    // RDKit✔️✔️:           // together: for N'ary numerical system, where N is the
    // RDKit✔️✔️:           // maximum possible value of the pair of numbers. The
    // RDKit✔️✔️:           // following will generate unique integers:
    // RDKit✔️✔️:           // hash(a, b) = a + b * N
    // RDKit✔️✔️:           auto big_seed = rdcast<size_t>(params->randomSeed);
    // RDKit✔️✔️:           size_t max_val = std::max(ci + 1, big_seed);
    // RDKit✔️✔️:           size_t big_num = big_seed + max_val * (ci + 1);
    // RDKit✔️✔️:           // only grab the first 31 bits xor'd with the next 31 bits to
    // RDKit✔️✔️:           // make sure its positive, careful, the 'ULL' is important
    // RDKit✔️✔️:           // here, 0x7fffffff is the 'int' type because of C default
    // RDKit✔️✔️:           // number semantics and that we definitely don't want!
    // RDKit✔️✔️:           const size_t positive_int_mask = 0x7fffffffULL;
    // RDKit✔️✔️:           size_t folded_num =
    // RDKit✔️✔️:               (big_num & positive_int_mask) ^ (big_num >> 31ULL);
    // RDKit✔️✔️:           new_seed = rdcast<int>(folded_num & positive_int_mask);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     CHECK_INVARIANT(new_seed >= -1,
    // RDKit✔️✔️:                     "Something went wrong calculating a new seed");
    // END RDKIT CPP BLOCK DGeomHelpers::detail::embedHelper_ per-conformer seed policy
    assert!(
        random_seed >= -1,
        "random seed must either be positive, zero, or negative one"
    );
    let ci_plus_one = i32::try_from(conformer_index + 1).expect("conformer index too large");
    let mut new_seed = random_seed;
    if new_seed > -1 {
        if enable_sequential_random_seeds {
            new_seed = new_seed
                .checked_add(ci_plus_one)
                .expect("sequential conformer seed overflow");
        } else if !rdkit_embedder_multiplication_overflows(ci_plus_one, random_seed) {
            new_seed = ci_plus_one * random_seed;
        } else {
            let big_seed = random_seed as usize;
            let max_val = (conformer_index + 1).max(big_seed);
            let big_num = big_seed + max_val * (conformer_index + 1);
            const POSITIVE_INT_MASK: usize = 0x7fff_ffff;
            let folded_num = (big_num & POSITIVE_INT_MASK) ^ (big_num >> 31);
            new_seed = (folded_num & POSITIVE_INT_MASK) as i32;
        }
    }
    assert!(
        new_seed >= -1,
        "Something went wrong calculating a new seed"
    );
    new_seed
}

fn pick_random_dist_mat(mmat: &BoundsMatrix, dist_mat: &mut SymmMatrix, seed: i32) -> f64 {
    // BEGIN RDKIT CPP FUNCTION DistGeom::pickRandomDistMat seed overload (DistGeomUtils.cpp:35-41)
    // RDKit✔️✔️: double pickRandomDistMat(const BoundsMatrix &mmat,
    // RDKit✔️✔️:                          RDNumeric::SymmMatrix<double> &distMat, int seed) {
    // RDKit✔️✔️:   if (seed > 0) {
    // RDKit✔️✔️:     RDKit::getRandomGenerator(seed);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return pickRandomDistMat(mmat, distMat, RDKit::getDoubleRandomSource());
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DistGeom::pickRandomDistMat seed overload
    if seed > 0 {
        RDKIT_DISTGEOM_RNG.with(|rng| {
            *rng.borrow_mut() = RdkitDistgeomMinStdRand::new(seed);
        });
    }
    RDKIT_DISTGEOM_RNG
        .with(|rng| pick_random_dist_mat_with_rng(mmat, dist_mat, &mut *rng.borrow_mut()))
}

fn pick_random_dist_mat_with_rng<R: RdkitDoubleRng>(
    mmat: &BoundsMatrix,
    dist_mat: &mut SymmMatrix,
    rng: &mut R,
) -> f64 {
    // BEGIN RDKIT CPP FUNCTION DistGeom::pickRandomDistMat RNG overload (DistGeomUtils.cpp:43-68)
    // RDKit✔️✔️: double pickRandomDistMat(const BoundsMatrix &mmat,
    // RDKit✔️✔️:                          RDNumeric::SymmMatrix<double> &distMat,
    // RDKit✔️✔️:                          RDKit::double_source_type &rng) {
    // RDKit✔️✔️:   // make sure the sizes match up
    // RDKit✔️✔️:   unsigned int npt = mmat.numRows();
    // RDKit✔️✔️:   CHECK_INVARIANT(npt == distMat.numRows(), "Size mismatch");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   double largestVal = -1.0;
    // RDKit✔️✔️:   double *ddata = distMat.getData();
    // RDKit✔️✔️:   for (unsigned int i = 1; i < npt; i++) {
    // RDKit✔️✔️:     unsigned int id = i * (i + 1) / 2;
    // RDKit✔️✔️:     for (unsigned int j = 0; j < i; j++) {
    // RDKit✔️✔️:       double ub = mmat.getUpperBound(i, j);
    // RDKit✔️✔️:       double lb = mmat.getLowerBound(i, j);
    // RDKit✔️✔️:       CHECK_INVARIANT(ub >= lb, "");
    // RDKit✔️✔️:       double rval = rng();
    // RDKit✔️✔️:       // std::cerr<<i<<"-"<<j<<": "<<rval<<std::endl;
    // RDKit✔️✔️:       double d = lb + (rval) * (ub - lb);
    // RDKit✔️✔️:       ddata[id + j] = d;
    // RDKit✔️✔️:       if (d > largestVal) {
    // RDKit✔️✔️:         largestVal = d;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return largestVal;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DistGeom::pickRandomDistMat RNG overload
    let npt = mmat.num_rows();
    assert_eq!(npt, dist_mat.num_rows(), "Size mismatch");

    let mut largest_val = -1.0;
    for i in 1..npt {
        let id = i * (i + 1) / 2;
        for j in 0..i {
            let ub = mmat.get_upper(i, j);
            let lb = mmat.get_lower(i, j);
            assert!(ub >= lb);
            let rval = rng.next_unit_f64();
            let d = lb + rval * (ub - lb);
            dist_mat.data[id + j] = d;
            if d > largest_val {
                largest_val = d;
            }
        }
    }
    #[cfg(test)]
    embedder_test_trace_update(|trace| {
        if trace.low_level.random_dist_preview.is_none() {
            trace.low_level.random_dist_preview =
                Some(dist_mat.data.iter().take(16).copied().collect());
            if env::var("COSMOLKIT_DG_TRACE_FULL_RANDOM_DIST")
                .ok()
                .as_deref()
                == Some("1")
            {
                trace.low_level.random_dist_full = Some(dist_mat.data.clone());
            }
            trace.low_level.random_dist_sum = Some(dist_mat.data.iter().sum());
            trace.low_level.random_dist_sum_sq =
                Some(dist_mat.data.iter().map(|value| value * value).sum());
            trace.low_level.random_dist_largest = Some(largest_val);
        }
    });
    largest_val
}

fn rdkit_symm_matrix_vector_multiply(a: &SymmMatrix, x: &[f64], y: &mut [f64]) {
    // BEGIN RDKIT CPP FUNCTION RDNumeric::multiply SymmMatrix-Vector overload (SymmMatrix.h:313-335)
    // RDKit✔️✔️: Vector<TYPE> &multiply(const SymmMatrix<TYPE> &A, const Vector<TYPE> &x,
    // RDKit✔️✔️:                        Vector<TYPE> &y) {
    // RDKit✔️✔️:   unsigned int aSize = A.numRows();
    // RDKit✔️✔️:   CHECK_INVARIANT(aSize == x.size(), "Size mismatch during multiplication");
    // RDKit✔️✔️:   CHECK_INVARIANT(aSize == y.size(), "Size mismatch during multiplication");
    // RDKit✔️✔️:   const TYPE *xData = x.getData();
    // RDKit✔️✔️:   const TYPE *aData = A.getData();
    // RDKit✔️✔️:   TYPE *yData = y.getData();
    // RDKit✔️✔️:   for (unsigned int i = 0; i < aSize; i++) {
    // RDKit✔️✔️:     yData[i] = (TYPE)(0.0);
    // RDKit✔️✔️:     unsigned int idA = i * (i + 1) / 2;
    // RDKit✔️✔️:     for (unsigned int j = 0; j < i + 1; j++) {
    // RDKit✔️✔️:       // idA = i*(i+1)/2 + j;
    // RDKit✔️✔️:       yData[i] += (aData[idA] * xData[j]);
    // RDKit✔️✔️:       idA++;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     idA--;
    // RDKit✔️✔️:     for (unsigned int j = i + 1; j < aSize; j++) {
    // RDKit✔️✔️:       // idA = j*(j+1)/2 + i;
    // RDKit✔️✔️:       idA += j;
    // RDKit✔️✔️:       yData[i] += (aData[idA] * xData[j]);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION RDNumeric::multiply SymmMatrix-Vector overload
    let a_size = a.num_rows();
    assert_eq!(a_size, x.len(), "Size mismatch during multiplication");
    assert_eq!(a_size, y.len(), "Size mismatch during multiplication");
    for i in 0..a_size {
        y[i] = 0.0;
        let mut id_a = i * (i + 1) / 2;
        for xj in x.iter().take(i + 1) {
            y[i] += a.data[id_a] * *xj;
            id_a += 1;
        }
        id_a -= 1;
        for (j, xj) in x.iter().enumerate().take(a_size).skip(i + 1) {
            id_a += j;
            y[i] += a.data[id_a] * *xj;
        }
    }
}

fn rdkit_vector_largest_abs_val_id(data: &[f64]) -> usize {
    // BEGIN RDKIT CPP METHOD RDNumeric::Vector::largestAbsValId (Vector.h:202-213)
    // RDKit✔️✔️: constexpr unsigned int largestAbsValId() const {
    // RDKit✔️✔️:   TYPE res = (TYPE)(-1.0);
    // RDKit✔️✔️:   unsigned int i, id = d_size;
    // RDKit✔️✔️:   TYPE *data = d_data.get();
    // RDKit✔️✔️:   for (i = 0; i < d_size; i++) {
    // RDKit✔️✔️:     if (fabs(data[i]) > res) {
    // RDKit✔️✔️:       res = fabs(data[i]);
    // RDKit✔️✔️:       id = i;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return id;
    // RDKit✔️✔️: }
    // END RDKIT CPP METHOD RDNumeric::Vector::largestAbsValId
    let mut res = -1.0;
    let mut id = data.len();
    for (i, value) in data.iter().enumerate() {
        if value.abs() > res {
            res = value.abs();
            id = i;
        }
    }
    id
}

fn rdkit_vector_normalize(data: &mut [f64]) {
    // BEGIN RDKIT CPP METHOD RDNumeric::Vector::normalize (Vector.h:258-264)
    // RDKit✔️✔️: constexpr void normalize() {
    // RDKit✔️✔️:   TYPE val = this->normL2();
    // RDKit✔️✔️:   if (val < zero_tolerance) {
    // RDKit✔️✔️:     throw std::runtime_error("Cannot normalize a zero length vector");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   (*this) /= val;
    // RDKit✔️✔️: }
    // END RDKIT CPP METHOD RDNumeric::Vector::normalize
    let val = data.iter().map(|v| v * v).sum::<f64>().sqrt();
    assert!(val >= f64::EPSILON, "Cannot normalize a zero length vector");
    for item in data {
        *item /= val;
    }
}

fn rdkit_vector_set_to_random(size: usize, seed: i32) -> Result<Vec<f64>, DgBoundsError> {
    // BEGIN RDKIT CPP METHOD RDNumeric::Vector::setToRandom (Vector.h:267-287)
    // RDKit✔️✔️: void setToRandom(unsigned int seed = 0) {
    // RDKit✔️✔️:   // we want to get our own RNG here instead of using the global
    // RDKit✔️✔️:   // one.  This is related to Issue285.
    // RDKit✔️✔️:   RDKit::rng_type generator(42u);
    // RDKit✔️✔️:   RDKit::uniform_double dist(0, 1.0);
    // RDKit✔️✔️:   RDKit::double_source_type randSource(generator, dist);
    // RDKit✔️✔️:   if (seed > 0) {
    // RDKit✔️✔️:     generator.seed(seed);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     // we can't initialize using only clock(), because it's possible
    // RDKit✔️✔️:     // that we'll get here fast enough that clock() will return 0
    // RDKit✔️✔️:     // and generator.seed(0) is an error:
    // RDKit✔️✔️:     generator.seed(clock() + 1);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   unsigned int i;
    // RDKit✔️✔️:   TYPE *data = d_data.get();
    // RDKit✔️✔️:   for (i = 0; i < d_size; i++) {
    // RDKit✔️✔️:     data[i] = randSource();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   this->normalize();
    // RDKit✔️✔️: }
    // END RDKIT CPP METHOD RDNumeric::Vector::setToRandom
    let effective_seed = if seed > 0 {
        seed
    } else {
        rdkit_clock_seed()?.wrapping_add(1)
    };
    let mut rng = RdkitDistgeomMinStdRand::new(effective_seed);
    let mut data = (0..size).map(|_| rng.next_unit_f64()).collect::<Vec<_>>();
    rdkit_vector_normalize(&mut data);
    Ok(data)
}

fn power_eigen_solver(
    num_eig: usize,
    mat: &mut SymmMatrix,
    eigen_values: &mut [f64],
    mut eigen_vectors: Option<&mut DoubleMatrix>,
    seed: i32,
) -> Result<bool, DgBoundsError> {
    // BEGIN RDKIT CPP FUNCTION RDNumeric::EigenSolvers::powerEigenSolver (PowerEigenSolver.cpp:20-100)
    // RDKit✔️✔️: bool powerEigenSolver(unsigned int numEig, DoubleSymmMatrix &mat,
    // RDKit✔️✔️:                       DoubleVector &eigenValues, DoubleMatrix *eigenVectors,
    // RDKit✔️✔️:                       int seed) {
    // RDKit✔️✔️:   const unsigned int MAX_ITERATIONS = 1000;
    // RDKit✔️✔️:   const double TOLERANCE = 0.001;
    // RDKit✔️✔️:   const double HUGE_EIGVAL = 1.0e10;
    // RDKit✔️✔️:   const double TINY_EIGVAL = 1.0e-10;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // first check all the sizes
    // RDKit✔️✔️:   unsigned int N = mat.numRows();
    // RDKit✔️✔️:   CHECK_INVARIANT(eigenValues.size() >= numEig, "");
    // RDKit✔️✔️:   CHECK_INVARIANT(numEig <= N, "");
    // RDKit✔️✔️:   if (eigenVectors) {
    // RDKit✔️✔️:     unsigned int evRows, evCols;
    // RDKit✔️✔️:     evRows = eigenVectors->numRows();
    // RDKit✔️✔️:     evCols = eigenVectors->numCols();
    // RDKit✔️✔️:     CHECK_INVARIANT(evCols >= N, "");
    // RDKit✔️✔️:     CHECK_INVARIANT(evRows >= numEig, "");
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION RDNumeric::EigenSolvers::powerEigenSolver
    const MAX_ITERATIONS: usize = 1000;
    const TOLERANCE: f64 = 0.001;
    const HUGE_EIGVAL: f64 = 1.0e10;
    const TINY_EIGVAL: f64 = 1.0e-10;

    let n = mat.num_rows();
    assert!(eigen_values.len() >= num_eig);
    assert!(num_eig <= n);
    if let Some(eig_vecs) = eigen_vectors.as_ref() {
        assert!(eig_vecs.num_cols() >= n);
        assert!(eig_vecs.num_rows() >= num_eig);
    }

    // RDKit✔️✔️:   unsigned int ei;
    // RDKit✔️✔️:   double eigVal, prevVal;
    // RDKit✔️✔️:   bool converged = false;
    // RDKit✔️✔️:   unsigned int i, j, id, iter, evalId;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   DoubleVector v(N), z(N);
    // RDKit✔️✔️:   if (seed <= 0) {
    // RDKit✔️✔️:     seed = clock();
    // RDKit✔️✔️:   }
    let mut effective_seed = if seed > 0 { seed } else { rdkit_clock_seed()? };
    let mut converged = false;
    let mut z = vec![0.0; n];
    for (ei, eig_value_slot) in eigen_values.iter_mut().enumerate().take(num_eig) {
        // RDKit✔️✔️:   for (ei = 0; ei < numEig; ei++) {
        // RDKit✔️✔️:     eigVal = -HUGE_EIGVAL;
        // RDKit✔️✔️:     seed += ei;
        // RDKit✔️✔️:     v.setToRandom(seed);
        let mut eig_val = -HUGE_EIGVAL;
        effective_seed += ei as i32;
        let mut v = rdkit_vector_set_to_random(n, effective_seed)?;

        converged = false;
        for _iter in 0..MAX_ITERATIONS {
            // RDKit✔️✔️:     for (iter = 0; iter < MAX_ITERATIONS; iter++) {
            // RDKit✔️✔️:       // z = mat*v
            // RDKit✔️✔️:       multiply(mat, v, z);
            // RDKit✔️✔️:       prevVal = eigVal;
            // RDKit✔️✔️:       evalId = z.largestAbsValId();
            // RDKit✔️✔️:       eigVal = z.getVal(evalId);
            rdkit_symm_matrix_vector_multiply(mat, &v, &mut z);
            let prev_val = eig_val;
            let eval_id = rdkit_vector_largest_abs_val_id(&z);
            eig_val = z[eval_id];

            // RDKit✔️✔️:       if (fabs(eigVal) < TINY_EIGVAL) {
            // RDKit✔️✔️:         break;
            // RDKit✔️✔️:       }
            if eig_val.abs() < TINY_EIGVAL {
                break;
            }

            // RDKit✔️✔️:       // compute the next estimate for the eigen vector
            // RDKit✔️✔️:       v.assign(z);
            // RDKit✔️✔️:       v /= eigVal;
            // RDKit✔️✔️:       if (fabs(eigVal - prevVal) < TOLERANCE) {
            // RDKit✔️✔️:         converged = true;
            // RDKit✔️✔️:         break;
            // RDKit✔️✔️:       }
            v.copy_from_slice(&z);
            for item in &mut v {
                *item /= eig_val;
            }
            if (eig_val - prev_val).abs() < TOLERANCE {
                converged = true;
                break;
            }
        }
        // RDKit✔️✔️:     if (!converged) {
        // RDKit✔️✔️:       break;
        // RDKit✔️✔️:     }
        if !converged {
            break;
        }
        // RDKit✔️✔️:     v.normalize();
        rdkit_vector_normalize(&mut v);

        // RDKit✔️✔️:     // save this is a eigen vector and value
        // RDKit✔️✔️:     // directly access the data instead of setVal so that we save time
        // RDKit✔️✔️:     double *vdata = v.getData();
        // RDKit✔️✔️:     if (eigenVectors) {
        // RDKit✔️✔️:       id = ei * eigenVectors->numCols();
        // RDKit✔️✔️:       double *eigVecData = eigenVectors->getData();
        // RDKit✔️✔️:       for (i = 0; i < N; i++) {
        // RDKit✔️✔️:         eigVecData[id + i] = vdata[i];
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        if let Some(eig_vecs) = eigen_vectors.as_deref_mut() {
            let id = ei * eig_vecs.num_cols();
            eig_vecs.data[id..id + n].copy_from_slice(&v);
        }
        // RDKit✔️✔️:     eigenValues[ei] = eigVal;
        *eig_value_slot = eig_val;

        // RDKit✔️✔️:     // now remove this eigen vector space out of the matrix
        // RDKit✔️✔️:     double *matData = mat.getData();
        // RDKit✔️✔️:     for (i = 0; i < N; i++) {
        // RDKit✔️✔️:       id = i * (i + 1) / 2;
        // RDKit✔️✔️:       for (j = 0; j < i + 1; j++) {
        // RDKit✔️✔️:         matData[id + j] -= (eigVal * vdata[i] * vdata[j]);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        for i in 0..n {
            let id = i * (i + 1) / 2;
            for j in 0..=i {
                mat.data[id + j] -= eig_val * v[i] * v[j];
            }
        }
    }
    // RDKit✔️✔️:   return converged;
    // RDKit✔️✔️: }
    Ok(converged)
}

fn compute_initial_coords(
    dist_mat: &SymmMatrix,
    positions: &mut [Vec<f64>],
    rand_neg_eig: bool,
    num_zero_fail: usize,
    seed: i32,
) -> Result<bool, DgBoundsError> {
    // BEGIN RDKIT CPP FUNCTION DistGeom::computeInitialCoords seed overload (DistGeomUtils.cpp:70-80)
    // RDKit✔️✔️: bool computeInitialCoords(const RDNumeric::SymmMatrix<double> &distMat,
    // RDKit✔️✔️:                           RDGeom::PointPtrVect &positions, bool randNegEig,
    // RDKit✔️✔️:                           unsigned int numZeroFail, int seed) {
    // RDKit✔️✔️:   if (seed > 0) {
    // RDKit✔️✔️:     RDKit::getRandomGenerator(seed);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return computeInitialCoords(distMat, positions,
    // RDKit✔️✔️:                               RDKit::getDoubleRandomSource(), randNegEig,
    // RDKit✔️✔️:                               numZeroFail);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DistGeom::computeInitialCoords seed overload
    if seed > 0 {
        RDKIT_DISTGEOM_RNG.with(|rng| {
            *rng.borrow_mut() = RdkitDistgeomMinStdRand::new(seed);
        });
    }
    RDKIT_DISTGEOM_RNG.with(|rng| {
        compute_initial_coords_with_rng(
            dist_mat,
            positions,
            &mut *rng.borrow_mut(),
            rand_neg_eig,
            num_zero_fail,
        )
    })
}

fn compute_initial_coords_with_rng<R: RdkitDoubleRng>(
    dist_mat: &SymmMatrix,
    positions: &mut [Vec<f64>],
    rng: &mut R,
    rand_neg_eig: bool,
    num_zero_fail: usize,
) -> Result<bool, DgBoundsError> {
    // BEGIN RDKIT CPP FUNCTION DistGeom::computeInitialCoords RNG overload (DistGeomUtils.cpp:81-164)
    // RDKit✔️✔️: bool computeInitialCoords(const RDNumeric::SymmMatrix<double> &distMat,
    // RDKit✔️✔️:                           RDGeom::PointPtrVect &positions,
    // RDKit✔️✔️:                           RDKit::double_source_type &rng, bool randNegEig,
    // RDKit✔️✔️:                           unsigned int numZeroFail) {
    // RDKit✔️✔️:   unsigned int N = distMat.numRows();
    // RDKit✔️✔️:   unsigned int nPt = positions.size();
    // RDKit✔️✔️:   CHECK_INVARIANT(nPt == N, "Size mismatch");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   unsigned int dim = positions.front()->dimension();
    // END RDKIT CPP FUNCTION DistGeom::computeInitialCoords RNG overload
    let n = dist_mat.num_rows();
    assert_eq!(positions.len(), n, "Size mismatch");
    assert!(!positions.is_empty());
    let dim = positions[0].len();
    let trace_row64 = row64_distgeom_trace_enabled(n);
    let trace_row61 = row61_distgeom_trace_enabled(n);

    // RDKit✔️✔️:   const double *data = distMat.getData();
    // RDKit✔️✔️:   RDNumeric::SymmMatrix<double> sqMat(N), T(N, 0.0);
    // RDKit✔️✔️:   RDNumeric::DoubleMatrix eigVecs(dim, N);
    // RDKit✔️✔️:   RDNumeric::DoubleVector eigVals(dim);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   double *sqDat = sqMat.getData();
    let data = dist_mat.get_data();
    let mut sq_mat = SymmMatrix::new(n);
    let mut t = SymmMatrix::with_value(n, 0.0);
    let mut eig_vecs = DoubleMatrix::new(dim, n);
    let mut eig_vals = vec![0.0; dim];

    // RDKit✔️✔️:   unsigned int dSize = distMat.getDataSize();
    // RDKit✔️✔️:   double sumSqD2 = 0.0;
    // RDKit✔️✔️:   for (unsigned int i = 0; i < dSize; i++) {
    // RDKit✔️✔️:     sqDat[i] = data[i] * data[i];
    // RDKit✔️✔️:     sumSqD2 += sqDat[i];
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   sumSqD2 /= (N * N);
    let d_size = dist_mat.get_data_size();
    let mut sum_sq_d2 = 0.0;
    for (i, value) in data.iter().enumerate().take(d_size) {
        sq_mat.data[i] = value * value;
        sum_sq_d2 += sq_mat.data[i];
    }
    sum_sq_d2 /= (n * n) as f64;
    #[cfg(test)]
    embedder_test_trace_update(|trace| {
        trace.low_level.initial_sum_sq_d2.get_or_insert(sum_sq_d2);
    });
    if trace_row64 {
        let preview_len = data.len().min(12);
        let payload = data[..preview_len]
            .iter()
            .map(|v| format!("{v:.15}"))
            .collect::<Vec<_>>()
            .join(",");
        println!(
            "row64_compute_initial_coords dist_data=[{}] sum_sq_d2={:.15}",
            payload, sum_sq_d2
        );
    } else if trace_row61 {
        let preview_len = data.len().min(12);
        let payload = data[..preview_len]
            .iter()
            .map(|v| format!("{v:.15}"))
            .collect::<Vec<_>>()
            .join(",");
        println!(
            "row61_compute_initial_coords dist_data=[{}] sum_sq_d2={:.15}",
            payload, sum_sq_d2
        );
    }

    // RDKit✔️✔️:   RDNumeric::DoubleVector sqD0i(N, 0.0);
    // RDKit✔️✔️:   double *sqD0iData = sqD0i.getData();
    // RDKit✔️✔️:   for (unsigned int i = 0; i < N; i++) {
    // RDKit✔️✔️:     for (unsigned int j = 0; j < N; j++) {
    // RDKit✔️✔️:       sqD0iData[i] += sqMat.getVal(i, j);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     sqD0iData[i] /= N;
    // RDKit✔️✔️:     sqD0iData[i] -= sumSqD2;
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if ((sqD0iData[i] < EIGVAL_TOL) && (N > 3)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    let mut sq_d0i = vec![0.0; n];
    for (i, sq_d0i_value) in sq_d0i.iter_mut().enumerate().take(n) {
        for j in 0..n {
            *sq_d0i_value += sq_mat.get_val(i, j);
        }
        *sq_d0i_value /= n as f64;
        *sq_d0i_value -= sum_sq_d2;
        if *sq_d0i_value < EIGVAL_TOL && n > 3 {
            if trace_row64 {
                let first_fail_val = *sq_d0i_value;
                let preview_len = (i + 1).min(12);
                let payload = sq_d0i[..preview_len]
                    .iter()
                    .map(|v| format!("{v:.15}"))
                    .collect::<Vec<_>>()
                    .join(",");
                println!(
                    "row64_compute_initial_coords sq_d0i=[{}] first_fail_idx={} first_fail_val={:.15}",
                    payload, i, first_fail_val
                );
            } else if trace_row61 {
                let first_fail_val = *sq_d0i_value;
                let preview_len = (i + 1).min(12);
                let payload = sq_d0i[..preview_len]
                    .iter()
                    .map(|v| format!("{v:.15}"))
                    .collect::<Vec<_>>()
                    .join(",");
                println!(
                    "row61_compute_initial_coords sq_d0i=[{}] first_fail_idx={} first_fail_val={:.15}",
                    payload, i, first_fail_val
                );
            }
            return Ok(false);
        }
    }
    if trace_row64 {
        let preview_len = sq_d0i.len().min(12);
        let payload = sq_d0i[..preview_len]
            .iter()
            .map(|v| format!("{v:.15}"))
            .collect::<Vec<_>>()
            .join(",");
        println!("row64_compute_initial_coords sq_d0i=[{}]", payload);
    } else if trace_row61 {
        let preview_len = sq_d0i.len().min(12);
        let payload = sq_d0i[..preview_len]
            .iter()
            .map(|v| format!("{v:.15}"))
            .collect::<Vec<_>>()
            .join(",");
        println!("row61_compute_initial_coords sq_d0i=[{}]", payload);
    }
    #[cfg(test)]
    embedder_test_trace_update(|trace| {
        if trace.low_level.initial_sq_d0i_preview.is_none() {
            trace.low_level.initial_sq_d0i_preview =
                Some(sq_d0i.iter().take(16).copied().collect());
        }
    });

    // RDKit✔️✔️:   for (unsigned int i = 0; i < N; i++) {
    // RDKit✔️✔️:     for (unsigned int j = 0; j <= i; j++) {
    // RDKit✔️✔️:       double val = 0.5 * (sqD0iData[i] + sqD0iData[j] - sqMat.getVal(i, j));
    // RDKit✔️✔️:       T.setVal(i, j, val);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   unsigned int nEigs = (dim < N) ? dim : N;
    // RDKit✔️✔️:   RDNumeric::EigenSolvers::powerEigenSolver(nEigs, T, eigVals, eigVecs,
    // RDKit✔️✔️:                                             (int)(sumSqD2 * N));
    for i in 0..n {
        for j in 0..=i {
            let val = 0.5 * (sq_d0i[i] + sq_d0i[j] - sq_mat.get_val(i, j));
            t.set_val(i, j, val);
        }
    }
    let n_eigs = dim.min(n);
    power_eigen_solver(
        n_eigs,
        &mut t,
        &mut eig_vals,
        Some(&mut eig_vecs),
        (sum_sq_d2 * n as f64) as i32,
    )?;
    #[cfg(test)]
    embedder_test_trace_update(|trace| {
        if trace.low_level.initial_eig_vals_before_sqrt.is_none() {
            trace.low_level.initial_eig_vals_before_sqrt = Some(eig_vals.clone());
        }
    });
    if trace_row64 {
        let payload = eig_vals
            .iter()
            .map(|v| format!("{v:.15}"))
            .collect::<Vec<_>>()
            .join(",");
        println!(
            "row64_compute_initial_coords eig_vals_before_sqrt=[{}]",
            payload
        );
    } else if trace_row61 {
        let payload = eig_vals
            .iter()
            .map(|v| format!("{v:.15}"))
            .collect::<Vec<_>>()
            .join(",");
        println!(
            "row61_compute_initial_coords eig_vals_before_sqrt=[{}]",
            payload
        );
    }
    #[cfg(test)]
    embedder_test_trace_update(|trace| {
        if trace.low_level.initial_eig_vals_before_sqrt.is_none() {
            trace.low_level.initial_eig_vals_before_sqrt = Some(eig_vals.clone());
        }
    });

    // RDKit✔️✔️:   double *eigData = eigVals.getData();
    // RDKit✔️✔️:   bool foundNeg = false;
    // RDKit✔️✔️:   unsigned int zeroEigs = 0;
    // RDKit✔️✔️:   for (unsigned int i = 0; i < dim; i++) {
    // RDKit✔️✔️:     if (eigData[i] > EIGVAL_TOL) {
    // RDKit✔️✔️:       eigData[i] = sqrt(eigData[i]);
    // RDKit✔️✔️:     } else if (fabs(eigData[i]) < EIGVAL_TOL) {
    // RDKit✔️✔️:       eigData[i] = 0.0;
    // RDKit✔️✔️:       zeroEigs++;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       foundNeg = true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    let mut found_neg = false;
    let mut zero_eigs = 0;
    for eig_val in eig_vals.iter_mut().take(dim) {
        if *eig_val > EIGVAL_TOL {
            *eig_val = eig_val.sqrt();
        } else if eig_val.abs() < EIGVAL_TOL {
            *eig_val = 0.0;
            zero_eigs += 1;
        } else {
            found_neg = true;
        }
    }
    if trace_row64 {
        let payload = eig_vals
            .iter()
            .map(|v| format!("{v:.15}"))
            .collect::<Vec<_>>()
            .join(",");
        println!(
            "row64_compute_initial_coords eig_vals_after_sqrt=[{}] found_neg={} zero_eigs={}",
            payload, found_neg, zero_eigs
        );
    } else if trace_row61 {
        let payload = eig_vals
            .iter()
            .map(|v| format!("{v:.15}"))
            .collect::<Vec<_>>()
            .join(",");
        println!(
            "row61_compute_initial_coords eig_vals_after_sqrt=[{}] found_neg={} zero_eigs={}",
            payload, found_neg, zero_eigs
        );
    }
    #[cfg(test)]
    embedder_test_trace_update(|trace| {
        if trace.low_level.initial_eig_vals_after_sqrt.is_none() {
            trace.low_level.initial_eig_vals_after_sqrt = Some(eig_vals.clone());
        }
    });
    // RDKit✔️✔️:   if ((foundNeg) && (!randNegEig)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if ((zeroEigs >= numZeroFail) && (N > 3)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    if found_neg && !rand_neg_eig {
        return Ok(false);
    }
    if zero_eigs >= num_zero_fail && n > 3 {
        return Ok(false);
    }

    // RDKit✔️✔️:   for (unsigned int i = 0; i < N; i++) {
    // RDKit✔️✔️:     RDGeom::Point *pt = positions[i];
    // RDKit✔️✔️:     for (unsigned int j = 0; j < dim; ++j) {
    // RDKit✔️✔️:       if (eigData[j] >= 0.0) {
    // RDKit✔️✔️:         (*pt)[j] = eigData[j] * eigVecs.getVal(j, i);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         // std::cerr<<"!!! "<<i<<"-"<<j<<": "<<eigData[j]<<std::endl;
    // RDKit✔️✔️:         (*pt)[j] = 1.0 - 2.0 * rng();
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    for (i, point) in positions.iter_mut().enumerate().take(n) {
        for j in 0..dim {
            if eig_vals[j] >= 0.0 {
                point[j] = eig_vals[j] * eig_vecs.get_val(j, i);
            } else {
                point[j] = 1.0 - 2.0 * rng.next_unit_f64();
            }
        }
    }
    Ok(true)
}

fn compute_random_coords(positions: &mut [Vec<f64>], box_size: f64, seed: i32) -> bool {
    // BEGIN RDKIT CPP FUNCTION DistGeom::computeRandomCoords seed overload (DistGeomUtils.cpp:166-172)
    // RDKit✔️✔️: bool computeRandomCoords(RDGeom::PointPtrVect &positions, double boxSize,
    // RDKit✔️✔️:                          int seed) {
    // RDKit✔️✔️:   if (seed > 0) {
    // RDKit✔️✔️:     RDKit::getRandomGenerator(seed);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return computeRandomCoords(positions, boxSize,
    // RDKit✔️✔️:                              RDKit::getDoubleRandomSource());
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DistGeom::computeRandomCoords seed overload
    if seed > 0 {
        RDKIT_DISTGEOM_RNG.with(|rng| {
            *rng.borrow_mut() = RdkitDistgeomMinStdRand::new(seed);
        });
    }
    RDKIT_DISTGEOM_RNG
        .with(|rng| compute_random_coords_with_rng(positions, box_size, &mut *rng.borrow_mut()))
}

fn compute_random_coords_with_rng<R: RdkitDoubleRng>(
    positions: &mut [Vec<f64>],
    box_size: f64,
    rng: &mut R,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION DistGeom::computeRandomCoords RNG overload (DistGeomUtils.cpp:173-183)
    // RDKit✔️✔️: bool computeRandomCoords(RDGeom::PointPtrVect &positions, double boxSize,
    // RDKit✔️✔️:                          RDKit::double_source_type &rng) {
    // RDKit✔️✔️:   CHECK_INVARIANT(boxSize > 0.0, "bad boxSize");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (auto pt : positions) {
    // RDKit✔️✔️:     for (unsigned int i = 0; i < pt->dimension(); ++i) {
    // RDKit✔️✔️:       (*pt)[i] = boxSize * (rng() - 0.5);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DistGeom::computeRandomCoords RNG overload
    assert!(box_size > 0.0, "bad boxSize");
    for point in positions {
        for coord in point {
            *coord = box_size * (rng.next_unit_f64() - 0.5);
        }
    }
    true
}

fn forcefield_vec_from_position(point: &[f64]) -> ForceFieldVec3 {
    assert!(
        (3..=4).contains(&point.len()),
        "unsupported point dimension"
    );
    if point.len() == 4 {
        ForceFieldVec3::new4(point[0], point[1], point[2], point[3])
    } else {
        ForceFieldVec3::new(point[0], point[1], point[2])
    }
}

fn construct_distgeom_forcefield(
    mmat: &BoundsMatrix,
    positions: &[Vec<f64>],
    csets: &[ChiralSetPtr],
    weight_chiral: f64,
    weight_fourth_dim: f64,
    extra_weights: Option<&BTreeMap<(usize, usize), f64>>,
    basin_size_tol: f64,
    fixed_pts: Option<&[bool]>,
) -> ForceField {
    // BEGIN RDKIT CPP FUNCTION DistGeom::constructForceField (DistGeomUtils.cpp:184-253)
    // RDKit✔️✔️: ForceFields::ForceField *constructForceField(
    // RDKit✔️✔️:     const BoundsMatrix &mmat, RDGeom::PointPtrVect &positions,
    // RDKit✔️✔️:     const VECT_CHIRALSET &csets, double weightChiral, double weightFourthDim,
    // RDKit✔️✔️:     std::map<std::pair<int, int>, double> *extraWeights, double basinSizeTol,
    // RDKit✔️✔️:     boost::dynamic_bitset<> *fixedPts) {
    // RDKit✔️✔️:   unsigned int N = mmat.numRows();
    // RDKit✔️✔️:   CHECK_INVARIANT(N == positions.size(), "");
    // RDKit✔️✔️:   auto *field = new ForceFields::ForceField(positions[0]->dimension());
    // RDKit✔️✔️:   field->positions().insert(field->positions().begin(), positions.begin(),
    // RDKit✔️✔️:                             positions.end());
    // END RDKIT CPP FUNCTION DistGeom::constructForceField
    let n = mmat.num_rows();
    assert_eq!(n, positions.len());
    assert!(!positions.is_empty());
    let dimension = positions[0].len();
    assert!((3..=4).contains(&dimension), "unsupported point dimension");
    assert!(
        positions.iter().all(|point| point.len() == dimension),
        "inconsistent point dimension"
    );
    if let Some(fixed_pts) = fixed_pts {
        assert!(fixed_pts.len() >= n, "bad fixed point bitset");
    }
    let mut field = ForceField::new(dimension);
    field.positions_mut().extend(
        positions
            .iter()
            .map(|point| forcefield_vec_from_position(point)),
    );

    // RDKit✔️✔️:   auto contrib = new DistViolationContribs(field);
    // RDKit✔️✔️:   for (unsigned int i = 1; i < N; i++) {
    // RDKit✔️✔️:     for (unsigned int j = 0; j < i; j++) {
    // RDKit✔️✔️:       if (fixedPts != nullptr && (*fixedPts)[i] && (*fixedPts)[j]) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       double w = 1.0;
    // RDKit✔️✔️:       double l = mmat.getLowerBound(i, j);
    // RDKit✔️✔️:       double u = mmat.getUpperBound(i, j);
    // RDKit✔️✔️:       bool includeIt = false;
    // RDKit✔️✔️:       if (extraWeights) {
    // RDKit✔️✔️:         std::map<std::pair<int, int>, double>::const_iterator mapIt;
    // RDKit✔️✔️:         mapIt = extraWeights->find(std::make_pair(i, j));
    // RDKit✔️✔️:         if (mapIt != extraWeights->end()) {
    // RDKit✔️✔️:           w = mapIt->second;
    // RDKit✔️✔️:           includeIt = true;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (u - l <= basinSizeTol) {
    // RDKit✔️✔️:         includeIt = true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (includeIt) {
    // RDKit✔️✔️:         contrib->addContrib(i, j, u, l, w);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!contrib->empty()) {
    // RDKit✔️✔️:     field->contribs().push_back(ForceFields::ContribPtr(contrib));
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     delete contrib;
    // RDKit✔️✔️:   }
    let mut dist_contrib = DistViolationContribs::new(&field);
    for i in 1..n {
        for j in 0..i {
            if fixed_pts.is_some_and(|fixed_pts| fixed_pts[i] && fixed_pts[j]) {
                continue;
            }
            let mut weight = 1.0;
            let lower = mmat.get_lower(i, j);
            let upper = mmat.get_upper(i, j);
            let mut include_it = false;
            if let Some(extra_weights) = extra_weights
                && let Some(extra_weight) = extra_weights.get(&(i, j))
            {
                weight = *extra_weight;
                include_it = true;
            }
            if upper - lower <= basin_size_tol {
                include_it = true;
            }
            if include_it {
                dist_contrib.add_contrib(i, j, upper, lower, weight);
            }
        }
    }
    if !dist_contrib.empty() {
        field.add_contrib(Box::new(dist_contrib));
    }

    // RDKit✔️✔️:   // now add chiral constraints
    // RDKit✔️✔️:   if (weightChiral > 1.e-8) {
    // RDKit✔️✔️:     auto contrib = new ChiralViolationContribs(field);
    // RDKit✔️✔️:
    // RDKit✔️✔️:     for (const auto &cset : csets) {
    // RDKit✔️✔️:       contrib->addContrib(cset.get(), weightChiral);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (!contrib->empty()) {
    // RDKit✔️✔️:       field->contribs().push_back(ForceFields::ContribPtr(contrib));
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       delete contrib;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    if weight_chiral > 1.0e-8 {
        let mut chiral_contrib = ChiralViolationContribs::new(&field);
        for cset in csets {
            chiral_contrib.add_contrib(cset, weight_chiral);
        }
        if !chiral_contrib.empty() {
            field.add_contrib(Box::new(chiral_contrib));
        }
    }

    // RDKit✔️✔️:   // finally the contribution from the fourth dimension if we need to
    // RDKit✔️✔️:   if ((field->dimension() == 4) && (weightFourthDim > 1.e-8)) {
    // RDKit✔️✔️:     auto contrib = new FourthDimContribs(field);
    // RDKit✔️✔️:     for (unsigned int i = 0; i < N; i++) {
    // RDKit✔️✔️:       contrib->addContrib(i, weightFourthDim);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (!contrib->empty()) {
    // RDKit✔️✔️:       field->contribs().push_back(ForceFields::ContribPtr(contrib));
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       delete contrib;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return field;
    // RDKit✔️✔️: }  // constructForceField
    if field.dimension() == 4 && weight_fourth_dim > 1.0e-8 {
        let mut fourth_dim_contrib = FourthDimContribs::new(&field);
        for i in 0..n {
            fourth_dim_contrib.add_contrib(i, weight_fourth_dim);
        }
        if !fourth_dim_contrib.empty() {
            field.add_contrib(Box::new(fourth_dim_contrib));
        }
    }
    field
}

fn add_improper_torsion_terms(
    ff: &mut ForceField,
    force_scaling_factor: f64,
    improper_atoms: &[Vec<i32>],
    is_improper_constrained: &mut [bool],
) {
    // BEGIN RDKIT CPP FUNCTION DistGeom::addImproperTorsionTerms (DistGeomUtils.cpp:255-307)
    // RDKit✔️✔️: void addImproperTorsionTerms(ForceFields::ForceField *ff,
    // RDKit✔️✔️:                              double forceScalingFactor,
    // RDKit✔️✔️:                              const std::vector<std::vector<int>> &improperAtoms,
    // RDKit✔️✔️:                              boost::dynamic_bitset<> &isImproperConstrained) {
    // RDKit✔️✔️:   PRECONDITION(ff, "bad force field");
    // RDKit✔️✔️:   auto inversionContribs =
    // RDKit✔️✔️:       std::make_unique<ForceFields::UFF::InversionContribs>(ff);
    // END RDKIT CPP FUNCTION DistGeom::addImproperTorsionTerms
    let mut inversion_contribs = InversionContribs::new(ff);

    // RDKit✔️✔️:   for (const auto &improperAtom : improperAtoms) {
    // RDKit✔️✔️:     std::vector<int> n(4);
    for improper_atom in improper_atoms {
        let mut n = [0_usize; 4];
        // RDKit✔️✔️:     for (unsigned int i = 0; i < 3; ++i) {
        for i in 0..3 {
            // RDKit✔️✔️:       n[1] = 1;
            n[1] = 1;
            // RDKit✔️✔️:       switch (i) {
            match i {
                // RDKit✔️✔️:         case 0:
                // RDKit✔️✔️:           n[0] = 0;
                // RDKit✔️✔️:           n[2] = 2;
                // RDKit✔️✔️:           n[3] = 3;
                // RDKit✔️✔️:           break;
                0 => {
                    n[0] = 0;
                    n[2] = 2;
                    n[3] = 3;
                }
                // RDKit✔️✔️:         case 1:
                // RDKit✔️✔️:           n[0] = 0;
                // RDKit✔️✔️:           n[2] = 3;
                // RDKit✔️✔️:           n[3] = 2;
                // RDKit✔️✔️:           break;
                1 => {
                    n[0] = 0;
                    n[2] = 3;
                    n[3] = 2;
                }
                // RDKit✔️✔️:         case 2:
                // RDKit✔️✔️:           n[0] = 2;
                // RDKit✔️✔️:           n[2] = 3;
                // RDKit✔️✔️:           n[3] = 0;
                // RDKit✔️✔️:           break;
                2 => {
                    n[0] = 2;
                    n[2] = 3;
                    n[3] = 0;
                }
                _ => unreachable!("loop bounds guarantee 0..3"),
            }

            // RDKit✔️✔️:       inversionContribs->addContrib(
            // RDKit✔️✔️:           improperAtom[n[0]], improperAtom[n[1]], improperAtom[n[2]],
            // RDKit✔️✔️:           improperAtom[n[3]], improperAtom[4],
            // RDKit✔️✔️:           static_cast<bool>(improperAtom[5]), forceScalingFactor);
            inversion_contribs.add_contrib(
                improper_atom[n[0]] as usize,
                improper_atom[n[1]] as usize,
                improper_atom[n[2]] as usize,
                improper_atom[n[3]] as usize,
                improper_atom[4],
                improper_atom[5] != 0,
                force_scaling_factor,
            );

            // RDKit✔️✔️:       isImproperConstrained[improperAtom[n[1]]] = 1;
            is_improper_constrained[improper_atom[n[1]] as usize] = true;
        }
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!inversionContribs->empty()) {
    // RDKit✔️✔️:     ff->contribs().push_back(std::move(inversionContribs));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    if !inversion_contribs.empty() {
        ff.add_contrib(Box::new(inversion_contribs));
    }
}

fn add_experimental_torsion_terms(
    ff: &mut ForceField,
    etkdg_details: &CrystalFFDetails,
    atom_pairs: &mut [bool],
    num_atoms: usize,
) {
    // BEGIN RDKIT CPP FUNCTION DistGeom::addExperimentalTorsionTerms (DistGeomUtils.cpp:309-340)
    // RDKit✔️✔️: void addExperimentalTorsionTerms(
    // RDKit✔️✔️:     ForceFields::ForceField *ff,
    // RDKit✔️✔️:     const ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails,
    // RDKit✔️✔️:     boost::dynamic_bitset<> &atomPairs, unsigned int numAtoms) {
    // RDKit✔️✔️:   PRECONDITION(ff, "bad force field");
    // RDKit✔️✔️:   auto torsionContribs =
    // RDKit✔️✔️:       std::make_unique<ForceFields::CrystalFF::TorsionAngleContribs>(ff);
    // END RDKIT CPP FUNCTION DistGeom::addExperimentalTorsionTerms
    let mut torsion_contribs = TorsionAngleContribs::new(ff);

    // RDKit✔️✔️:   for (unsigned int t = 0; t < etkdgDetails.expTorsionAtoms.size(); ++t) {
    for t in 0..etkdg_details.exp_torsion_atoms.len() {
        // RDKit✔️✔️:     int i = etkdgDetails.expTorsionAtoms[t][0];
        // RDKit✔️✔️:     int j = etkdgDetails.expTorsionAtoms[t][1];
        // RDKit✔️✔️:     int k = etkdgDetails.expTorsionAtoms[t][2];
        // RDKit✔️✔️:     int l = etkdgDetails.expTorsionAtoms[t][3];
        let i = etkdg_details.exp_torsion_atoms[t][0] as usize;
        let j = etkdg_details.exp_torsion_atoms[t][1] as usize;
        let k = etkdg_details.exp_torsion_atoms[t][2] as usize;
        let l = etkdg_details.exp_torsion_atoms[t][3] as usize;

        // RDKit✔️✔️:     if (i < l) {
        // RDKit✔️✔️:       atomPairs[i * numAtoms + l] = 1;
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       atomPairs[l * numAtoms + i] = 1;
        // RDKit✔️✔️:     }
        if i < l {
            atom_pairs[i * num_atoms + l] = true;
        } else {
            atom_pairs[l * num_atoms + i] = true;
        }

        // RDKit✔️✔️:     torsionContribs->addContrib(i, j, k, l,
        // RDKit✔️✔️:                                 etkdgDetails.expTorsionAngles[t].second,
        // RDKit✔️✔️:                                 etkdgDetails.expTorsionAngles[t].first);
        torsion_contribs.add_contrib(
            i,
            j,
            k,
            l,
            etkdg_details.exp_torsion_angles[t].1.clone(),
            etkdg_details.exp_torsion_angles[t].0.clone(),
        );
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!torsionContribs->empty()) {
    // RDKit✔️✔️:     ff->contribs().push_back(std::move(torsionContribs));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    if !torsion_contribs.is_empty() {
        ff.add_contrib(Box::new(torsion_contribs));
    }
}

fn add_12_terms(
    ff: &mut ForceField,
    etkdg_details: &CrystalFFDetails,
    atom_pairs: &mut [bool],
    positions: &[ForceFieldVec3],
    force_constant: f64,
    num_atoms: usize,
) {
    // BEGIN RDKIT CPP FUNCTION DistGeom::add12Terms (DistGeomUtils.cpp:342-383)
    // RDKit✔️✔️: void add12Terms(ForceFields::ForceField *ff,
    // RDKit✔️✔️:                 const ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails,
    // RDKit✔️✔️:                 boost::dynamic_bitset<> &atomPairs,
    // RDKit✔️✔️:                 RDGeom::Point3DPtrVect &positions, double forceConstant,
    // RDKit✔️✔️:                 unsigned int numAtoms) {
    // RDKit✔️✔️:   PRECONDITION(ff, "bad force field");
    // RDKit✔️✔️:   auto distContribs =
    // RDKit✔️✔️:       std::make_unique<ForceFields::DistanceConstraintContribs>(ff);
    // END RDKIT CPP FUNCTION DistGeom::add12Terms
    let mut dist_contribs = DistanceConstraintContribs::new(ff);

    // RDKit✔️✔️:   for (const auto &bond : etkdgDetails.bonds) {
    for &(first, second) in &etkdg_details.bonds {
        // RDKit✔️✔️:     unsigned int i = bond.first;
        // RDKit✔️✔️:     unsigned int j = bond.second;
        let i = first as usize;
        let j = second as usize;

        // RDKit✔️✔️:     if (i < j) {
        // RDKit✔️✔️:       atomPairs[i * numAtoms + j] = 1;
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       atomPairs[j * numAtoms + i] = 1;
        // RDKit✔️✔️:     }
        if i < j {
            atom_pairs[i * num_atoms + j] = true;
        } else {
            atom_pairs[j * num_atoms + i] = true;
        }

        // RDKit✔️✔️:     double d = ((*positions[i]) - (*positions[j])).length();
        let d = (positions[i] - positions[j]).length();
        // RDKit✔️✔️:     distContribs->addContrib(i, j, d - KNOWN_DIST_TOL, d + KNOWN_DIST_TOL,
        // RDKit✔️✔️:                              forceConstant);
        dist_contribs.add_contrib(i, j, d - KNOWN_DIST_TOL, d + KNOWN_DIST_TOL, force_constant);
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!distContribs->empty()) {
    // RDKit✔️✔️:     ff->contribs().push_back(std::move(distContribs));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    if !dist_contribs.empty() {
        ff.add_contrib(Box::new(dist_contribs));
    }
}

fn add_13_terms(
    ff: &mut ForceField,
    etkdg_details: &CrystalFFDetails,
    atom_pairs: &mut [bool],
    positions: &[ForceFieldVec3],
    force_constant: f64,
    is_improper_constrained: &[bool],
    use_basic_knowledge: bool,
    mmat: &BoundsMatrix,
    num_atoms: usize,
) {
    // BEGIN RDKIT CPP FUNCTION DistGeom::add13Terms (DistGeomUtils.cpp:385-455)
    // RDKit✔️✔️: void add13Terms(ForceFields::ForceField *ff,
    // RDKit✔️✔️:                 const ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails,
    // RDKit✔️✔️:                 boost::dynamic_bitset<> &atomPairs,
    // RDKit✔️✔️:                 RDGeom::Point3DPtrVect &positions, double forceConstant,
    // RDKit✔️✔️:                 const boost::dynamic_bitset<> &isImproperConstrained,
    // RDKit✔️✔️:                 bool useBasicKnowledge, const BoundsMatrix &mmat,
    // RDKit✔️✔️:                 unsigned int numAtoms) {
    // RDKit✔️✔️:   PRECONDITION(ff, "bad force field");
    // RDKit✔️✔️:   auto distContribs =
    // RDKit✔️✔️:       std::make_unique<ForceFields::DistanceConstraintContribs>(ff);
    // RDKit✔️✔️:   auto angleContribs =
    // RDKit✔️✔️:       std::make_unique<ForceFields::AngleConstraintContribs>(ff);
    // END RDKIT CPP FUNCTION DistGeom::add13Terms
    let mut dist_contribs = DistanceConstraintContribs::new(ff);
    let mut angle_contribs = AngleConstraintContribs::new(ff);

    // RDKit✔️✔️:   for (const auto &angle : etkdgDetails.angles) {
    for angle in &etkdg_details.angles {
        // RDKit✔️✔️:     unsigned int i = angle[0];
        // RDKit✔️✔️:     unsigned int j = angle[1];
        // RDKit✔️✔️:     unsigned int k = angle[2];
        let i = angle[0] as usize;
        let j = angle[1] as usize;
        let k = angle[2] as usize;

        // RDKit✔️✔️:     if (i < k) {
        // RDKit✔️✔️:       atomPairs[i * numAtoms + k] = 1;
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       atomPairs[k * numAtoms + i] = 1;
        // RDKit✔️✔️:     }
        if i < k {
            atom_pairs[i * num_atoms + k] = true;
        } else {
            atom_pairs[k * num_atoms + i] = true;
        }

        // RDKit✔️✔️:     // check for triple bonds
        // RDKit✔️✔️:     if (useBasicKnowledge && angle[3]) {
        // RDKit✔️✔️:       angleContribs->addContrib(i, j, k, 179.0, 180.0, 1);
        if use_basic_knowledge && angle[3] != 0 {
            angle_contribs.add_contrib(i, j, k, 179.0, 180.0, 1.0);
        // RDKit✔️✔️:     } else if (isImproperConstrained[j]) {
        // RDKit✔️✔️:       distContribs->addContrib(i, k, mmat.getLowerBound(i, k),
        // RDKit✔️✔️:                                mmat.getUpperBound(i, k), forceConstant);
        } else if is_improper_constrained[j] {
            dist_contribs.add_contrib(
                i,
                k,
                mmat.get_lower(i, k),
                mmat.get_upper(i, k),
                force_constant,
            );
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       double d = ((*positions[i]) - (*positions[k])).length();
        // RDKit✔️✔️:       distContribs->addContrib(i, k, d - KNOWN_DIST_TOL, d + KNOWN_DIST_TOL,
        // RDKit✔️✔️:                                forceConstant);
        // RDKit✔️✔️:     }
        } else {
            let d = (positions[i] - positions[k]).length();
            dist_contribs.add_contrib(i, k, d - KNOWN_DIST_TOL, d + KNOWN_DIST_TOL, force_constant);
        }
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!angleContribs->empty()) {
    // RDKit✔️✔️:     ff->contribs().push_back(std::move(angleContribs));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!distContribs->empty()) {
    // RDKit✔️✔️:     ff->contribs().push_back(std::move(distContribs));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    if !angle_contribs.empty() {
        ff.add_contrib(Box::new(angle_contribs));
    }
    if !dist_contribs.empty() {
        ff.add_contrib(Box::new(dist_contribs));
    }
}

fn add_long_range_distance_constraints(
    ff: &mut ForceField,
    etkdg_details: &CrystalFFDetails,
    atom_pairs: &[bool],
    positions: &[ForceFieldVec3],
    known_distance_force_constant: f64,
    mmat: &BoundsMatrix,
    num_atoms: usize,
) {
    // BEGIN RDKIT CPP FUNCTION DistGeom::addLongRangeDistanceConstraints (DistGeomUtils.cpp:457-505)
    // RDKit✔️✔️: void addLongRangeDistanceConstraints(
    // RDKit✔️✔️:     ForceFields::ForceField *ff,
    // RDKit✔️✔️:     const ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails,
    // RDKit✔️✔️:     const boost::dynamic_bitset<> &atomPairs, RDGeom::Point3DPtrVect &positions,
    // RDKit✔️✔️:     double knownDistanceForceConstant, const BoundsMatrix &mmat,
    // RDKit✔️✔️:     unsigned int numAtoms) {
    // RDKit✔️✔️:   PRECONDITION(ff, "bad force field");
    // RDKit✔️✔️:   auto distContribs =
    // RDKit✔️✔️:       std::make_unique<ForceFields::DistanceConstraintContribs>(ff);
    // RDKit✔️✔️:   double fdist = knownDistanceForceConstant;
    // END RDKIT CPP FUNCTION DistGeom::addLongRangeDistanceConstraints
    let mut dist_contribs = DistanceConstraintContribs::new(ff);

    // RDKit✔️✔️:   for (unsigned int i = 1; i < numAtoms; ++i) {
    for i in 1..num_atoms {
        // RDKit✔️✔️:     for (unsigned int j = 0; j < i; ++j) {
        for j in 0..i {
            // RDKit✔️✔️:       if (!atomPairs[j * numAtoms + i]) {
            if !atom_pairs[j * num_atoms + i] {
                // RDKit✔️✔️:         fdist = etkdgDetails.boundsMatForceScaling * 10.0;
                let mut fdist = etkdg_details.bounds_mat_force_scaling * 10.0;
                // RDKit✔️✔️:         double l = mmat.getLowerBound(i, j);
                // RDKit✔️✔️:         double u = mmat.getUpperBound(i, j);
                let mut l = mmat.get_lower(i, j);
                let mut u = mmat.get_upper(i, j);

                // RDKit✔️✔️:         if (!etkdgDetails.constrainedAtoms.empty() &&
                // RDKit✔️✔️:             etkdgDetails.constrainedAtoms[i] &&
                // RDKit✔️✔️:             etkdgDetails.constrainedAtoms[j]) {
                if !etkdg_details.constrained_atoms.is_empty()
                    && etkdg_details.constrained_atoms[i]
                    && etkdg_details.constrained_atoms[j]
                {
                    // RDKit✔️✔️:           // we're constrained, so use very tight bounds
                    // RDKit✔️✔️:           l = u = ((*positions[i]) - (*positions[j])).length();
                    let d = (positions[i] - positions[j]).length();
                    l = d;
                    u = d;
                    // RDKit✔️✔️:           l -= KNOWN_DIST_TOL;
                    // RDKit✔️✔️:           u += KNOWN_DIST_TOL;
                    l -= KNOWN_DIST_TOL;
                    u += KNOWN_DIST_TOL;
                    // RDKit✔️✔️:           fdist = knownDistanceForceConstant;
                    fdist = known_distance_force_constant;
                    // RDKit✔️✔️:         }
                }
                // RDKit✔️✔️:         distContribs->addContrib(i, j, l, u, fdist);
                dist_contribs.add_contrib(i, j, l, u, fdist);
            }
            // RDKit✔️✔️:       }
        }
        // RDKit✔️✔️:     }
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!distContribs->empty()) {
    // RDKit✔️✔️:     ff->contribs().push_back(std::move(distContribs));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    if !dist_contribs.empty() {
        ff.add_contrib(Box::new(dist_contribs));
    }
}

fn construct_3d_forcefield(
    mmat: &BoundsMatrix,
    positions: &[ForceFieldVec3],
    etkdg_details: &CrystalFFDetails,
) -> ForceField {
    // BEGIN RDKIT CPP FUNCTION DistGeom::construct3DForceField (DistGeomUtils.cpp:495-524)
    // RDKit✔️✔️: ForceFields::ForceField *construct3DForceField(
    // RDKit✔️✔️:     const BoundsMatrix &mmat, RDGeom::Point3DPtrVect &positions,
    // RDKit✔️✔️:     const ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails) {
    // RDKit✔️✔️:   unsigned int N = mmat.numRows();
    // RDKit✔️✔️:   CHECK_INVARIANT(N == positions.size(), "");
    // RDKit✔️✔️:   CHECK_INVARIANT(etkdgDetails.expTorsionAtoms.size() ==
    // RDKit✔️✔️:                       etkdgDetails.expTorsionAngles.size(),
    // RDKit✔️✔️:                   "");
    let n = mmat.num_rows();
    assert_eq!(n, positions.len());
    assert_eq!(
        etkdg_details.exp_torsion_atoms.len(),
        etkdg_details.exp_torsion_angles.len()
    );

    // RDKit✔️✔️:   auto *field = new ForceFields::ForceField(positions[0]->dimension());
    // RDKit✔️✔️:   field->positions().insert(field->positions().begin(), positions.begin(),
    // RDKit✔️✔️:                             positions.end());
    // Rust stores RDGeom::Point3D-equivalent coordinates by value; Point3D
    // dimension is therefore fixed at three for this overload.
    assert!(!positions.is_empty());
    let mut field = ForceField::new(3);
    field.positions_mut().extend_from_slice(positions);

    // RDKit✔️✔️:   // keep track which atoms are 1,2-, 1,3- or 1,4-restrained
    // RDKit✔️✔️:   boost::dynamic_bitset<> atomPairs(N * N);
    // RDKit✔️✔️:   // don't add 1-3 Distances constraints for angles where the
    // RDKit✔️✔️:   // central atom of the angle is the central atom of an improper torsion.
    // RDKit✔️✔️:   boost::dynamic_bitset<> isImproperConstrained(N);
    let mut atom_pairs = vec![false; n * n];
    let mut is_improper_constrained = vec![false; n];

    // RDKit✔️✔️:   addExperimentalTorsionTerms(field, etkdgDetails, atomPairs, N);
    // RDKit✔️✔️:   addImproperTorsionTerms(field, 10.0, etkdgDetails.improperAtoms,
    // RDKit✔️✔️:                           isImproperConstrained);
    // RDKit✔️✔️:   add12Terms(field, etkdgDetails, atomPairs, positions,
    // RDKit✔️✔️:              KNOWN_DIST_FORCE_CONSTANT, N);
    // RDKit✔️✔️:   add13Terms(field, etkdgDetails, atomPairs, positions,
    // RDKit✔️✔️:              KNOWN_DIST_FORCE_CONSTANT, isImproperConstrained, true, mmat, N);
    // RDKit✔️✔️:   // minimum distance for all other atom pairs that aren't constrained
    // RDKit✔️✔️:   addLongRangeDistanceConstraints(field, etkdgDetails, atomPairs, positions,
    // RDKit✔️✔️:                                   KNOWN_DIST_FORCE_CONSTANT, mmat, N);
    // RDKit✔️✔️:   return field;
    // RDKit✔️✔️: }  // construct3DForceField
    add_experimental_torsion_terms(&mut field, etkdg_details, &mut atom_pairs, n);
    add_improper_torsion_terms(
        &mut field,
        10.0,
        &etkdg_details.improper_atoms,
        &mut is_improper_constrained,
    );
    add_12_terms(
        &mut field,
        etkdg_details,
        &mut atom_pairs,
        positions,
        KNOWN_DIST_FORCE_CONSTANT,
        n,
    );
    add_13_terms(
        &mut field,
        etkdg_details,
        &mut atom_pairs,
        positions,
        KNOWN_DIST_FORCE_CONSTANT,
        &is_improper_constrained,
        true,
        mmat,
        n,
    );
    add_long_range_distance_constraints(
        &mut field,
        etkdg_details,
        &atom_pairs,
        positions,
        KNOWN_DIST_FORCE_CONSTANT,
        mmat,
        n,
    );
    field
}

fn construct_plain_3d_forcefield(
    mmat: &BoundsMatrix,
    positions: &[ForceFieldVec3],
    etkdg_details: &CrystalFFDetails,
) -> ForceField {
    // BEGIN RDKIT CPP FUNCTION DistGeom::constructPlain3DForceField (DistGeomUtils.cpp:545-573)
    // RDKit✔️✔️: ForceFields::ForceField *constructPlain3DForceField(
    // RDKit✔️✔️:     const BoundsMatrix &mmat, RDGeom::Point3DPtrVect &positions,
    // RDKit✔️✔️:     const ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails) {
    // RDKit✔️✔️:   unsigned int N = mmat.numRows();
    // RDKit✔️✔️:   CHECK_INVARIANT(N == positions.size(), "");
    // RDKit✔️✔️:   CHECK_INVARIANT(etkdgDetails.expTorsionAtoms.size() ==
    // RDKit✔️✔️:                       etkdgDetails.expTorsionAngles.size(),
    // RDKit✔️✔️:                   "");
    let n = mmat.num_rows();
    assert_eq!(n, positions.len());
    assert_eq!(
        etkdg_details.exp_torsion_atoms.len(),
        etkdg_details.exp_torsion_angles.len()
    );

    // RDKit✔️✔️:   auto *field = new ForceFields::ForceField(positions[0]->dimension());
    // RDKit✔️✔️:   field->positions().insert(field->positions().begin(), positions.begin(),
    // RDKit✔️✔️:                             positions.end());
    assert!(!positions.is_empty());
    let mut field = ForceField::new(3);
    field.positions_mut().extend_from_slice(positions);

    // RDKit✔️✔️:   // keep track which atoms are 1,2-, 1,3- or 1,4-restrained
    // RDKit✔️✔️:   boost::dynamic_bitset<> atomPairs(N * N);
    // RDKit✔️✔️:   // don't add 1-3 Distances constraints for angles where the
    // RDKit✔️✔️:   // central atom of the angle is the central atom of an improper torsion.
    // RDKit✔️✔️:   boost::dynamic_bitset<> isImproperConstrained(N);
    let mut atom_pairs = vec![false; n * n];
    let is_improper_constrained = vec![false; n];

    // RDKit✔️✔️:   addExperimentalTorsionTerms(field, etkdgDetails, atomPairs, N);
    // RDKit✔️✔️:   add12Terms(field, etkdgDetails, atomPairs, positions,
    // RDKit✔️✔️:              KNOWN_DIST_FORCE_CONSTANT, N);
    // RDKit✔️✔️:   add13Terms(field, etkdgDetails, atomPairs, positions,
    // RDKit✔️✔️:              KNOWN_DIST_FORCE_CONSTANT, isImproperConstrained, false, mmat, N);
    // RDKit✔️✔️:   // minimum distance for all other atom pairs that aren't constrained
    // RDKit✔️✔️:   addLongRangeDistanceConstraints(field, etkdgDetails, atomPairs, positions,
    // RDKit✔️✔️:                                   KNOWN_DIST_FORCE_CONSTANT, mmat, N);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return field;
    // RDKit✔️✔️: }  // constructPlain3DForceField
    add_experimental_torsion_terms(&mut field, etkdg_details, &mut atom_pairs, n);
    add_12_terms(
        &mut field,
        etkdg_details,
        &mut atom_pairs,
        positions,
        KNOWN_DIST_FORCE_CONSTANT,
        n,
    );
    add_13_terms(
        &mut field,
        etkdg_details,
        &mut atom_pairs,
        positions,
        KNOWN_DIST_FORCE_CONSTANT,
        &is_improper_constrained,
        false,
        mmat,
        n,
    );
    add_long_range_distance_constraints(
        &mut field,
        etkdg_details,
        &atom_pairs,
        positions,
        KNOWN_DIST_FORCE_CONSTANT,
        mmat,
        n,
    );
    field
}

fn construct_3d_improper_forcefield_from_parts(
    mmat: &BoundsMatrix,
    positions: &[ForceFieldVec3],
    improper_atoms: &[Vec<i32>],
    angles: &[Vec<i32>],
    atom_nums: &[i32],
) -> ForceField {
    // BEGIN RDKIT CPP FUNCTION DistGeom::construct3DImproperForceField (DistGeomUtils.cpp:575-612)
    // RDKit✔️✔️: ForceFields::ForceField *construct3DImproperForceField(
    // RDKit✔️✔️:     const BoundsMatrix &mmat, RDGeom::Point3DPtrVect &positions,
    // RDKit✔️✔️:     const std::vector<std::vector<int>> &improperAtoms,
    // RDKit✔️✔️:     const std::vector<std::vector<int>> &angles,
    // RDKit✔️✔️:     const std::vector<int> &atomNums) {
    // RDKit✔️✔️:   RDUNUSED_PARAM(atomNums);
    let _ = atom_nums;
    // RDKit✔️✔️:   unsigned int N = mmat.numRows();
    // RDKit✔️✔️:   CHECK_INVARIANT(N == positions.size(), "");
    let n = mmat.num_rows();
    assert_eq!(n, positions.len());

    // RDKit✔️✔️:   auto *field = new ForceFields::ForceField(positions[0]->dimension());
    // RDKit✔️✔️:   field->positions().insert(field->positions().begin(), positions.begin(),
    // RDKit✔️✔️:                             positions.end());
    assert!(!positions.is_empty());
    let mut field = ForceField::new(3);
    field.positions_mut().extend_from_slice(positions);

    // RDKit✔️✔️:   // improper torsions / out-of-plane bend / inversion
    // RDKit✔️✔️:   double oobForceScalingFactor = 10.0;
    // RDKit✔️✔️:   boost::dynamic_bitset<> isImproperConstrained(N);
    // RDKit✔️✔️:   addImproperTorsionTerms(field, oobForceScalingFactor, improperAtoms,
    // RDKit✔️✔️:                           isImproperConstrained);
    let oob_force_scaling_factor = 10.0;
    let mut is_improper_constrained = vec![false; n];
    add_improper_torsion_terms(
        &mut field,
        oob_force_scaling_factor,
        improper_atoms,
        &mut is_improper_constrained,
    );

    // RDKit✔️✔️:   // Check that SP Centers have an angle of 180 degrees.
    // RDKit✔️✔️:   auto angleContribs =
    // RDKit✔️✔️:       std::make_unique<ForceFields::AngleConstraintContribs>(field);
    let mut angle_contribs = AngleConstraintContribs::new(&field);
    // RDKit✔️✔️:   for (const auto &angle : angles) {
    for angle in angles {
        // RDKit✔️✔️:     if (angle[3]) {
        if angle[3] != 0 {
            // RDKit✔️✔️:       angleContribs->addContrib(angle[0], angle[1], angle[2], 179.0, 180.0,
            // RDKit✔️✔️:                                 oobForceScalingFactor);
            angle_contribs.add_contrib(
                angle[0] as usize,
                angle[1] as usize,
                angle[2] as usize,
                179.0,
                180.0,
                oob_force_scaling_factor,
            );
        }
        // RDKit✔️✔️:     }
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!angleContribs->empty()) {
    // RDKit✔️✔️:     field->contribs().push_back(std::move(angleContribs));
    // RDKit✔️✔️:   }
    if !angle_contribs.empty() {
        field.add_contrib(Box::new(angle_contribs));
    }
    // RDKit✔️✔️:   return field;
    // RDKit✔️✔️: }  // construct3DImproperForceField
    field
}

fn construct_3d_improper_forcefield(
    mmat: &BoundsMatrix,
    positions: &[ForceFieldVec3],
    etkdg_details: &CrystalFFDetails,
) -> ForceField {
    // BEGIN RDKIT CPP INLINE OVERLOAD DistGeom::construct3DImproperForceField (DistGeomUtils.h:215-220)
    // RDKit✔️✔️: inline ForceFields::ForceField *construct3DImproperForceField(
    // RDKit✔️✔️:     const BoundsMatrix &mmat, RDGeom::Point3DPtrVect &positions,
    // RDKit✔️✔️:     const ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails) {
    // RDKit✔️✔️:   return construct3DImproperForceField(
    // RDKit✔️✔️:       mmat, positions, etkdgDetails.improperAtoms, etkdgDetails.angles,
    // RDKit✔️✔️:       etkdgDetails.atomNums);
    // RDKit✔️✔️: }
    construct_3d_improper_forcefield_from_parts(
        mmat,
        positions,
        &etkdg_details.improper_atoms,
        &etkdg_details.angles,
        &etkdg_details.atom_nums,
    )
}

fn construct_3d_forcefield_with_cpci(
    mmat: &BoundsMatrix,
    positions: &[ForceFieldVec3],
    etkdg_details: &CrystalFFDetails,
    cpci: &BTreeMap<(usize, usize), f64>,
) -> ForceField {
    // BEGIN RDKIT CPP FUNCTION DistGeom::construct3DForceField CPCI overload (DistGeomUtils.cpp:526-543)
    // RDKit✔️✔️: ForceFields::ForceField *construct3DForceField(
    // RDKit✔️✔️:     const BoundsMatrix &mmat, RDGeom::Point3DPtrVect &positions,
    // RDKit✔️✔️:     const ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails,
    // RDKit✔️✔️:     const std::map<std::pair<unsigned int, unsigned int>, double> &CPCI) {
    // RDKit✔️✔️:   auto *field = construct3DForceField(mmat, positions, etkdgDetails);
    let mut field = construct_3d_forcefield(mmat, positions, etkdg_details);

    // RDKit✔️✔️:   bool is1_4 = false;
    // RDKit✔️✔️:   // double dielConst = 1.0;
    // RDKit✔️✔️:   boost::uint8_t dielModel = 1;
    // RDKit✔️✔️:   auto *contrib = new ForceFields::MMFF::EleContrib(field);
    let is_1_4 = false;
    let diel_model = 1;
    let mut contrib = NonbondedContrib::new(&field);

    // RDKit✔️✔️:   field->contribs().emplace_back(contrib);
    // RDKit✔️✔️:   for (const auto &charge : CPCI) {
    // RDKit✔️✔️:     contrib->addTerm(charge.first.first, charge.first.second, charge.second,
    // RDKit✔️✔️:                      dielModel, is1_4);
    // RDKit✔️✔️:   }
    for (&(idx1, idx2), &charge_term) in cpci {
        contrib.add_term(idx1, idx2, None, true, charge_term, diel_model, is_1_4);
    }

    // RDKit✔️✔️:
    // RDKit✔️✔️:   return field;
    // RDKit✔️✔️: }
    field.add_contrib(Box::new(contrib));
    field
}

impl BoundsMatrix {
    fn check_and_set_bounds(
        &mut self,
        i: usize,
        j: usize,
        lb: f64,
        ub: f64,
    ) -> Result<(), DgBoundsError> {
        self.check_and_set_bounds_with_mode(i, j, lb, ub, false)
    }

    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::_checkAndSetBounds (BoundsMatrixBuilder.cpp:195-227)
    // RDKit✔️✔️: void _checkAndSetBounds(unsigned int i, unsigned int j, double lb, double ub,
    // RDKit✔️✔️:                         DistGeom::BoundsMatPtr mmat, bool setIfBetter = false) {
    // RDKit✔️✔️:   // get the existing bounds
    // RDKit✔️✔️:   double clb = mmat->getLowerBound(i, j);
    // RDKit✔️✔️:   double cub = mmat->getUpperBound(i, j);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   CHECK_INVARIANT(ub > lb, "upper bound not greater than lower bound");
    // RDKit✔️✔️:   CHECK_INVARIANT(lb > DIST12_DELTA || clb > DIST12_DELTA, "bad lower bound");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // Note: setIfBetter should ONLY be set if the distances are consistent;
    // RDKit✔️✔️:   // currently this is not the case, therefore, for now, we are pessimistic on
    // RDKit✔️✔️:   // the bounds
    // RDKit✔️✔️:   if (setIfBetter) {
    // RDKit✔️✔️:     double nlb = std::max(clb, lb);
    // RDKit✔️✔️:     double nub = std::min(cub, ub);
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (nub <= nlb) {
    // RDKit✔️✔️:       // if not overlapping ranges -> be conservative
    // RDKit✔️✔️:       nlb = std::min(clb, lb);
    // RDKit✔️✔️:       nub = std::max(cub, ub);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     mmat->setLowerBound(i, j, nlb);
    // RDKit✔️✔️:     mmat->setUpperBound(i, j, nub);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     if (clb <= DIST12_DELTA) {
    // RDKit✔️✔️:       mmat->setLowerBound(i, j, lb);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       if ((lb < clb) && (lb > DIST12_DELTA)) {
    // RDKit✔️✔️:         mmat->setLowerBound(i, j, lb);  // conservative bound setting
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if (cub >= MAX_UPPER) {  // FIX this
    // RDKit✔️✔️:       mmat->setUpperBound(i, j, ub);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       if ((ub > cub) && (ub < MAX_UPPER)) {
    // RDKit✔️✔️:         mmat->setUpperBound(i, j, ub);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::_checkAndSetBounds
    fn check_and_set_bounds_with_mode(
        &mut self,
        i: usize,
        j: usize,
        lb: f64,
        ub: f64,
        set_if_better: bool,
    ) -> Result<(), DgBoundsError> {
        let clb = self.get_lower(i, j);
        let cub = self.get_upper(i, j);

        if !matches!(ub.partial_cmp(&lb), Some(std::cmp::Ordering::Greater)) {
            return Err(Self::invalid_bounds(
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
            return Err(Self::invalid_bounds(
                "bad lower bound",
                i,
                j,
                lb,
                ub,
                clb,
                cub,
            ));
        }

        if set_if_better {
            let mut nlb = clb.max(lb);
            let mut nub = cub.min(ub);

            if nub <= nlb {
                nlb = clb.min(lb);
                nub = cub.max(ub);
            }

            self.set_lower(i, j, nlb)?;
            self.set_upper(i, j, nub)?;
        } else {
            if clb <= DIST12_DELTA {
                self.set_lower(i, j, lb)?;
            } else if lb < clb && lb > DIST12_DELTA {
                self.set_lower(i, j, lb)?;
            }

            if cub >= MAX_UPPER {
                self.set_upper(i, j, ub)?;
            } else if ub > cub && ub < MAX_UPPER {
                self.set_upper(i, j, ub)?;
            }
        }
        Ok(())
    }

    fn triangle_smooth(&mut self, tol: f64) -> bool {
        triangle_smooth_bounds_ptr(self, tol)
    }

    fn to_vec_vec(self) -> Vec<Vec<f64>> {
        self.data
    }
}

// BEGIN RDKIT CPP FUNCTION DistGeom::triangleSmoothBounds (TriangleSmooth.cpp:11-13)
// RDKit✔️✔️: bool triangleSmoothBounds(BoundsMatPtr boundsMat, double tol) {
// RDKit✔️✔️:   return triangleSmoothBounds(boundsMat.get(), tol);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DistGeom::triangleSmoothBounds
fn triangle_smooth_bounds_shared(bounds_mat: &mut BoundsMatrix, tol: f64) -> bool {
    triangle_smooth_bounds_ptr(bounds_mat, tol)
}

// BEGIN RDKIT CPP FUNCTION DistGeom::triangleSmoothBounds (TriangleSmooth.cpp:14-61)
// RDKit✔️❌: bool triangleSmoothBounds(BoundsMatrix *boundsMat, double tol) {
// RDKit✔️❌:   auto npt = boundsMat->numRows();
// RDKit✔️❌:   for (auto k = 0u; k < npt; k++) {
// RDKit✔️❌:     for (auto i = 0u; i < npt - 1; i++) {
// RDKit✔️❌:       if (i == k) {
// RDKit✔️❌:         continue;
// RDKit✔️❌:       }
// RDKit✔️❌:       auto ii = i;
// RDKit✔️❌:       auto ik = k;
// RDKit✔️❌:       if (ii > ik) {
// RDKit✔️❌:         std::swap(ii, ik);
// RDKit✔️❌:       }
// RDKit✔️❌:
// RDKit✔️❌:       const auto Uik = boundsMat->getValUnchecked(ii, ik);  // upper bound
// RDKit✔️❌:       const auto Lik = boundsMat->getValUnchecked(ik, ii);  // lower bound
// RDKit✔️❌:       for (auto j = i + 1; j < npt; j++) {
// RDKit✔️❌:         if (j == k) {
// RDKit✔️❌:           continue;
// RDKit✔️❌:         }
// RDKit✔️❌:         auto jj = j;
// RDKit✔️❌:         auto jk = k;
// RDKit✔️❌:         if (jj > jk) {
// RDKit✔️❌:           std::swap(jj, jk);
// RDKit✔️❌:         }
// RDKit✔️❌:         const auto Ukj = boundsMat->getValUnchecked(jj, jk);  // upper bound
// RDKit✔️❌:         const auto sumUikUkj = Uik + Ukj;
// RDKit✔️❌:         if (boundsMat->getValUnchecked(i, j) > sumUikUkj) {
// RDKit✔️❌:           // adjust the upper bound
// RDKit✔️❌:           boundsMat->setValUnchecked(i, j, sumUikUkj);
// RDKit✔️❌:         }
// RDKit✔️❌:
// RDKit✔️❌:         const auto diffLikUjk = Lik - Ukj;
// RDKit✔️❌:         const auto diffLjkUik = boundsMat->getValUnchecked(jk, jj) - Uik;
// RDKit✔️❌:         if (boundsMat->getValUnchecked(j, i) < diffLikUjk) {
// RDKit✔️❌:           // adjust the lower bound
// RDKit✔️❌:           boundsMat->setValUnchecked(j, i, diffLikUjk);
// RDKit✔️❌:         } else if (boundsMat->getValUnchecked(j, i) < diffLjkUik) {
// RDKit✔️❌:           // adjust the lower bound
// RDKit✔️❌:           boundsMat->setValUnchecked(j, i, diffLjkUik);
// RDKit✔️❌:         }
// RDKit✔️❌:         const auto lBound = boundsMat->getValUnchecked(j, i);
// RDKit✔️❌:         const auto uBound = boundsMat->getValUnchecked(i, j);
// RDKit✔️❌:         if (tol > 0. && (lBound - uBound) / lBound > 0. &&
// RDKit✔️❌:             (lBound - uBound) / lBound < tol) {
// RDKit✔️❌:           // adjust the upper bound
// RDKit✔️❌:           boundsMat->setValUnchecked(i, j, lBound);
// RDKit✔️❌:         } else if (lBound - uBound > 0.) {
// RDKit✔️❌:           return false;
// RDKit✔️❌:         }
// RDKit✔️❌:       }
// RDKit✔️❌:     }
// RDKit✔️❌:   }
// RDKit✔️❌:   return true;
// RDKit✔️❌: }
// END RDKIT CPP FUNCTION DistGeom::triangleSmoothBounds
//
// Performance note: the Rust port keeps RDKit's O(n^3) control flow and
// raw-triangle access pattern, but `Vec<Vec<f64>>` is less cache-friendly
// than RDKit's contiguous square-matrix backing.
fn triangle_smooth_bounds_ptr(bounds_mat: &mut BoundsMatrix, tol: f64) -> bool {
    let npt = bounds_mat.num_rows();
    if npt < 2 {
        return true;
    }

    for k in 0..npt {
        for i in 0..(npt - 1) {
            if i == k {
                continue;
            }
            let (ii, ik) = if i > k { (k, i) } else { (i, k) };

            let uik = bounds_mat.get_val(ii, ik);
            let lik = bounds_mat.get_val(ik, ii);
            for j in (i + 1)..npt {
                if j == k {
                    continue;
                }
                let (jj, jk) = if j > k { (k, j) } else { (j, k) };
                let ukj = bounds_mat.get_val(jj, jk);
                let sum_uik_ukj = uik + ukj;
                if bounds_mat.get_val(i, j) > sum_uik_ukj {
                    bounds_mat.set_val(i, j, sum_uik_ukj);
                }

                let diff_lik_ujk = lik - ukj;
                let diff_ljk_uik = bounds_mat.get_val(jk, jj) - uik;
                if bounds_mat.get_val(j, i) < diff_lik_ujk {
                    bounds_mat.set_val(j, i, diff_lik_ujk);
                } else if bounds_mat.get_val(j, i) < diff_ljk_uik {
                    bounds_mat.set_val(j, i, diff_ljk_uik);
                }
                let l_bound = bounds_mat.get_val(j, i);
                let u_bound = bounds_mat.get_val(i, j);
                let rel_gap = (l_bound - u_bound) / l_bound;
                if tol > 0.0 && rel_gap > 0.0 && rel_gap < tol {
                    bounds_mat.set_val(i, j, l_bound);
                } else if l_bound - u_bound > 0.0 {
                    return false;
                }
            }
        }
    }
    true
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::initBoundsMat (BoundsMatrixBuilder.cpp:1792-1800)
// RDKit✔️✔️: void initBoundsMat(DistGeom::BoundsMatrix *mmat, double defaultMin,
// RDKit✔️✔️:                    double defaultMax) {
// RDKit✔️✔️:   unsigned int npt = mmat->numRows();
// RDKit✔️✔️:
// RDKit✔️✔️:   for (unsigned int i = 1; i < npt; i++) {
// RDKit✔️✔️:     for (unsigned int j = 0; j < i; j++) {
// RDKit✔️✔️:       mmat->setUpperBound(i, j, defaultMax);
// RDKit✔️✔️:       mmat->setLowerBound(i, j, defaultMin);
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::initBoundsMat
fn init_bounds_mat_ptr(mmat: &mut BoundsMatrix, default_min: f64, default_max: f64) {
    let npt = mmat.num_rows();

    for i in 1..npt {
        for j in 0..i {
            mmat.set_upper_unchecked(i, j, default_max);
            mmat.set_lower_unchecked(i, j, default_min);
        }
    }
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::initBoundsMat (BoundsMatrixBuilder.cpp:1801-1804)
// RDKit✔️✔️: void initBoundsMat(DistGeom::BoundsMatPtr mmat, double defaultMin,
// RDKit✔️✔️:                    double defaultMax) {
// RDKit✔️✔️:   initBoundsMat(mmat.get(), defaultMin, defaultMax);
// RDKit✔️✔️: };
// END RDKIT CPP FUNCTION DGeomHelpers::initBoundsMat
fn init_bounds_mat_shared(mmat: &mut BoundsMatrix, default_min: f64, default_max: f64) {
    init_bounds_mat_ptr(mmat, default_min, default_max);
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

/// RDKit✔️❗: set12Bounds — bond-length derived 1-2 bounds
//
// BEGIN RDKIT CPP FUNCTION DGeomHelpers::set12Bounds (BoundsMatrixBuilder.cpp:238-295)
// RDKit✔️❗: void set12Bounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat,
// RDKit✔️❗:                  ComputedData &accumData) {
// RDKit✔️❗:   unsigned int npt = mmat->numRows();
// RDKit✔️❗:   CHECK_INVARIANT(npt == mol.getNumAtoms(), "Wrong size metric matrix");
// RDKit✔️❗:   CHECK_INVARIANT(accumData.bondLengths.size() >= mol.getNumBonds(),
// RDKit✔️❗:                   "Wrong size accumData");
// RDKit✔️❗:   auto [atomParams, foundAll] = UFF::getAtomTypes(mol);
// RDKit✔️❗:   CHECK_INVARIANT(atomParams.size() == mol.getNumAtoms(),
// RDKit✔️❗:                   "parameter vector size mismatch");
// RDKit✔️❗:   boost::dynamic_bitset<> squishAtoms(mol.getNumAtoms());
// RDKit✔️❗:   for (const auto bond : mol.bonds()) {
// RDKit✔️❗:     if (bond->getIsConjugated() &&
// RDKit✔️❗:         (bond->getBeginAtom()->getAtomicNum() > 10 ||
// RDKit✔️❗:          bond->getEndAtom()->getAtomicNum() > 10) &&
// RDKit✔️❗:         mol.getRingInfo() && mol.getRingInfo()->isInitialized() &&
// RDKit✔️❗:         mol.getRingInfo()->isBondInRingOfSize(bond->getIdx(), 5)) {
// RDKit✔️❗:       squishAtoms.set(bond->getBeginAtomIdx());
// RDKit✔️❗:       squishAtoms.set(bond->getEndAtomIdx());
// RDKit✔️❗:     }
// RDKit✔️❗:   }
// RDKit✔️❗:   for (const auto bond : mol.bonds()) {
// RDKit✔️❗:     auto begId = bond->getBeginAtomIdx();
// RDKit✔️❗:     auto endId = bond->getEndAtomIdx();
// RDKit✔️❗:     auto bOrder = bond->getBondTypeAsDouble();
// RDKit✔️❗:     if (atomParams[begId] && atomParams[endId] && bOrder > 0) {
// RDKit✔️❗:       auto bl = ForceFields::UFF::Utils::calcBondRestLength(
// RDKit✔️❗:           bOrder, atomParams[begId], atomParams[endId]);
// RDKit✔️❗:       double extraSquish = 0.0;
// RDKit✔️❗:       if (squishAtoms[begId] || squishAtoms[endId]) {
// RDKit✔️❗:         extraSquish = 0.2;  // empirical
// RDKit✔️❗:       }
// RDKit✔️❗:       accumData.bondLengths[bond->getIdx()] = bl;
// RDKit✔️❗:       mmat->setUpperBound(begId, endId, bl + extraSquish + DIST12_DELTA);
// RDKit✔️❗:       mmat->setLowerBound(begId, endId, bl - extraSquish - DIST12_DELTA);
// RDKit✔️❗:     } else {
// RDKit✔️❗:       auto vw1 = PeriodicTable::getTable()->getRvdw(
// RDKit✔️❗:           mol.getAtomWithIdx(begId)->getAtomicNum());
// RDKit✔️❗:       auto vw2 = PeriodicTable::getTable()->getRvdw(
// RDKit✔️❗:           mol.getAtomWithIdx(endId)->getAtomicNum());
// RDKit✔️❗:       auto bl = (vw1 + vw2) / 2;
// RDKit✔️❗:       accumData.bondLengths[bond->getIdx()] = bl;
// RDKit✔️❗:       mmat->setUpperBound(begId, endId, 1.5 * bl);
// RDKit✔️❗:       mmat->setLowerBound(begId, endId, .5 * bl);
// RDKit✔️❗:     }
// RDKit✔️❗:     unsigned int pid =
// RDKit✔️❗:         std::min(begId, endId) * mol.getNumAtoms() + std::max(begId, endId);
// RDKit✔️❗:     accumData.visited12Bounds.set(pid);
// RDKit✔️❗:   }
// RDKit✔️❗: }
// END RDKIT CPP FUNCTION DGeomHelpers::set12Bounds
//
// Rust maps RDKit's UFF atom typing through `get_atom_label_for_uff()` and
// `get_atom_types_for_uff()`. Local review found comparable asymptotic work
// but left perf unresolved because the Rust path materializes extra
// intermediate vectors for valence-derived hybridization and conjugation state.
fn set_12_bounds(
    mol: &Molecule,
    mmat: &mut BoundsMatrix,
    accum_data: &mut ComputedData,
) -> Result<(), DgBoundsError> {
    let npt = mmat.num_rows();
    if npt != mol.atoms().len() {
        return Err(DgBoundsError::GenerationFailed(
            "Wrong size metric matrix".to_string(),
        ));
    }
    if accum_data.bond_lengths.len() < mol.bonds().len() {
        return Err(DgBoundsError::GenerationFailed(
            "Wrong size accumData".to_string(),
        ));
    }
    let assignment = assign_valence(mol, ValenceModel::RdkitLike).map_err(|err| {
        DgBoundsError::GenerationFailed(format!(
            "RDKit UFF atom typing valence assignment failed: {err}"
        ))
    })?;
    let mut atom_degree = vec![0usize; mol.atoms().len()];
    for bond in mol.bonds() {
        atom_degree[bond.begin().index()] += 1;
        atom_degree[bond.end().index()] += 1;
    }
    let conjugated = compute_conjugated_bonds_for_uff(mol, &assignment, &atom_degree);
    let mut atom_has_conjugated_bond = vec![false; mol.atoms().len()];
    for (bond_index, bond) in mol.bonds().iter().enumerate() {
        if conjugated[bond_index] {
            atom_has_conjugated_bond[bond.begin().index()] = true;
            atom_has_conjugated_bond[bond.end().index()] = true;
        }
    }
    let hybridizations =
        compute_hybridizations_for_uff(mol, &assignment, &atom_degree, &atom_has_conjugated_bond);
    let total_valences = (0..mol.atoms().len())
        .map(|atom_index| atom_total_valence_for_uff(&assignment, atom_index))
        .collect::<Vec<_>>();
    let (atom_params, _found_all) = get_atom_types_for_uff(
        mol,
        &total_valences,
        &hybridizations,
        &atom_has_conjugated_bond,
    )
    .map_err(|err| {
        DgBoundsError::GenerationFailed(format!("RDKit UFF atom typing failed: {err}"))
    })?;

    let mut squish_atoms = vec![false; mol.atoms().len()];
    let ring_info = ring_info_for_distgeom(mol)?;
    for (bond_index, bond) in mol.bonds().iter().enumerate() {
        if conjugated[bond_index]
            && (mol.atoms()[bond.begin().index()].atomic_number() > 10
                || mol.atoms()[bond.end().index()].atomic_number() > 10)
            && is_bond_in_ring_of_size(ring_info.as_ref(), bond.id().index(), 5)
        {
            squish_atoms[bond.begin().index()] = true;
            squish_atoms[bond.end().index()] = true;
        }
    }

    for bond in mol.bonds() {
        let beg_id = bond.begin().index();
        let end_id = bond.end().index();
        let bond_order =
            crate::chemistry::valence::bond_type_as_double(bond.order()).map_err(|err| {
                DgBoundsError::GenerationFailed(format!(
                    "RDKit set12Bounds bond-order conversion failed: {err}"
                ))
            })?;
        if let (Some(param1), Some(param2)) = (atom_params[beg_id], atom_params[end_id]) {
            if bond_order > 0.0 {
                let bl = calc_bond_rest_length(bond_order, &param1, &param2);
                let extra_squish = if squish_atoms[beg_id] || squish_atoms[end_id] {
                    0.2
                } else {
                    0.0
                };
                accum_data.bond_lengths[bond.id().index()] = bl;
                mmat.set_upper(beg_id, end_id, bl + extra_squish + DIST12_DELTA)?;
                mmat.set_lower(beg_id, end_id, bl - extra_squish - DIST12_DELTA)?;
            } else {
                let bl = (vdw_radius(mol.atoms()[beg_id].atomic_number())
                    + vdw_radius(mol.atoms()[end_id].atomic_number()))
                    / 2.0;
                accum_data.bond_lengths[bond.id().index()] = bl;
                mmat.set_upper(beg_id, end_id, 1.5 * bl)?;
                mmat.set_lower(beg_id, end_id, 0.5 * bl)?;
            }
        } else {
            let bl = (vdw_radius(mol.atoms()[beg_id].atomic_number())
                + vdw_radius(mol.atoms()[end_id].atomic_number()))
                / 2.0;
            accum_data.bond_lengths[bond.id().index()] = bl;
            mmat.set_upper(beg_id, end_id, 1.5 * bl)?;
            mmat.set_lower(beg_id, end_id, 0.5 * bl)?;
        }
        let pid = beg_id.min(end_id) * mol.atoms().len() + beg_id.max(end_id);
        accum_data.visited12_bounds[pid] = true;
    }
    Ok(())
}

// ──────────────────────────────────────────────
// 1-3 distance bounds
// ──────────────────────────────────────────────

// BEGIN RDKIT CPP FUNCTION RDGeom::compute13Dist (Geometry/Utils.h:24-27)
// RDKit✔️✔️: inline double compute13Dist(double d1, double d2, double angle) {
// RDKit✔️✔️:   double res = d1 * d1 + d2 * d2 - 2 * d1 * d2 * cos(angle);
// RDKit✔️✔️:   return sqrt(res);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION RDGeom::compute13Dist
fn compute_13_dist(bl1: f64, bl2: f64, angle: f64) -> f64 {
    let res = bl1 * bl1 + bl2 * bl2 - 2.0 * bl1 * bl2 * angle.cos();
    res.sqrt()
}

// BEGIN RDKIT CPP FUNCTION RDGeom::compute14Dist3D (Geometry/Utils.h:50-65)
// RDKit✔️✔️: inline double compute14Dist3D(double d1, double d2, double d3, double ang12,
// RDKit✔️✔️:                               double ang23, double torAng) {
// RDKit✔️✔️:   Point3D p1(d1 * cos(ang12), d1 * sin(ang12), 0.0);
// RDKit✔️✔️:   Point3D p4(d2 - d3 * cos(ang23), d3 * sin(ang23), 0.0);
// RDKit✔️✔️:   Transform3D trans;
// RDKit✔️✔️:   trans.SetRotation(torAng, X_Axis);
// RDKit✔️✔️:   trans.TransformPoint(p4);
// RDKit✔️✔️:   p4 -= p1;
// RDKit✔️✔️:   return p4.length();
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION RDGeom::compute14Dist3D
fn compute_14_dist_3d(d1: f64, d2: f64, d3: f64, ang12: f64, ang23: f64, tor_ang: f64) -> f64 {
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

// BEGIN RDKIT CPP FUNCTION RDGeom::compute14DistCis (Geometry/Utils.h:74-81)
// RDKit✔️✔️: inline double compute14DistCis(double d1, double d2, double d3, double ang12,
// RDKit✔️✔️:                                double ang23) {
// RDKit✔️✔️:   double dx = d2 - d3 * cos(ang23) - d1 * cos(ang12);
// RDKit✔️✔️:   double dy = d3 * sin(ang23) - d1 * sin(ang12);
// RDKit✔️✔️:   double res = dx * dx + dy * dy;
// RDKit✔️✔️:   return sqrt(res);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION RDGeom::compute14DistCis
fn compute_14_dist_cis(d1: f64, d2: f64, d3: f64, ang12: f64, ang23: f64) -> f64 {
    let dx = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy = d3 * ang23.sin() - d1 * ang12.sin();
    (dx * dx + dy * dy).sqrt()
}

// BEGIN RDKIT CPP FUNCTION RDGeom::compute14DistTrans (Geometry/Utils.h:88-95)
// RDKit✔️✔️: inline double compute14DistTrans(double d1, double d2, double d3, double ang12,
// RDKit✔️✔️:                                  double ang23) {
// RDKit✔️✔️:   double dx = d2 - d3 * cos(ang23) - d1 * cos(ang12);
// RDKit✔️✔️:   double dy = d3 * sin(ang23) + d1 * sin(ang12);
// RDKit✔️✔️:   double res = dx * dx + dy * dy;
// RDKit✔️✔️:   return sqrt(res);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION RDGeom::compute14DistTrans
fn compute_14_dist_trans(d1: f64, d2: f64, d3: f64, ang12: f64, ang23: f64) -> f64 {
    let dx = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy = d3 * ang23.sin() + d1 * ang12.sin();
    (dx * dx + dy * dy).sqrt()
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
    rinfo: &crate::RingInfo,
) -> Result<(), DgBoundsError> {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::_set13BoundsHelper (BoundsMatrixBuilder.cpp:369-391)
    // RDKit✔️✔️: void _set13BoundsHelper(unsigned int aid1, unsigned int aid, unsigned int aid3,
    // RDKit✔️✔️:                         double angle, const ComputedData &accumData,
    // RDKit✔️✔️:                         DistGeom::BoundsMatPtr mmat, const ROMol &mol) {
    // RDKit✔️✔️:   auto bid1 = mol.getBondBetweenAtoms(aid1, aid)->getIdx();
    // RDKit✔️✔️:   auto bid2 = mol.getBondBetweenAtoms(aid, aid3)->getIdx();
    // RDKit✔️✔️:   auto dl = RDGeom::compute13Dist(accumData.bondLengths[bid1],
    // RDKit✔️✔️:                                   accumData.bondLengths[bid2], angle);
    // RDKit✔️✔️:   auto distTol = DIST13_TOL;
    // RDKit✔️✔️:   // Now increase the tolerance if we're outside of the first row of the
    // RDKit✔️✔️:   // periodic table.
    // RDKit✔️✔️:   if (isLargerSP2Atom(mol.getAtomWithIdx(aid1))) {
    // RDKit✔️✔️:     distTol *= 2;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (isLargerSP2Atom(mol.getAtomWithIdx(aid))) {
    // RDKit✔️✔️:     distTol *= 2;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (isLargerSP2Atom(mol.getAtomWithIdx(aid3))) {
    // RDKit✔️✔️:     distTol *= 2;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto du = dl + distTol;
    // RDKit✔️✔️:   dl -= distTol;
    // RDKit✔️✔️:   _checkAndSetBounds(aid1, aid3, dl, du, mmat);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::_set13BoundsHelper
    let Some(bid1) = bond_between_idx_simple(mol, aid1, aid) else {
        return Ok(());
    };
    let Some(bid2) = bond_between_idx_simple(mol, aid, aid3) else {
        return Ok(());
    };

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
    mmat.check_and_set_bounds(aid1, aid3, dl, du)
}

fn bond_between_idx(bond_idx: usize) -> Option<usize> {
    // placeholder — not used directly, see below
    None
}

// BEGIN RDKIT CPP FUNCTION isLargerSP2Atom (BoundsMatrixBuilder.cpp:364-367)
// RDKit✔️✔️: bool isLargerSP2Atom(const Atom *atom) {
// RDKit✔️✔️:   return atom->getAtomicNum() > 13 && atom->getHybridization() == Atom::SP2 &&
// RDKit✔️✔️:          atom->getOwningMol().getRingInfo()->numAtomRings(atom->getIdx());
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION isLargerSP2Atom
fn is_larger_sp2_atom_idx(mol: &Molecule, rinfo: &crate::RingInfo, idx: usize) -> bool {
    let atom = &mol.atoms()[idx];
    atom.atomic_number() > 13
        && atom.hybridization() == Hybridization::Sp2
        && rinfo.num_atom_rings(AtomId::new(idx)) > 0
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum DgPath14Kind {
    Cis,
    Trans,
    Other,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct DgPath14Prior {
    pub atoms: (usize, usize, usize, usize),
    pub kind: DgPath14Kind,
    pub lower_bound: f64,
    pub upper_bound: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DistType {
    Dist12,
    Dist13,
    Dist14,
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

    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::ComputedData::visitedBound (BoundsMatrixBuilder.cpp:92-96)
    // RDKit✔️✔️: bool visitedBound(unsigned int aid, DistType dt) const {
    // RDKit✔️✔️:   return visited12Bounds[aid] ||
    // RDKit✔️✔️:          ((dt >= DIST13) && visited13Bounds[aid]) ||
    // RDKit✔️✔️:          ((dt >= DIST14) && visited14Bounds[aid]);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::ComputedData::visitedBound
    fn visited_bound(&self, pid: usize, dist_type: DistType) -> bool {
        self.visited12_bounds[pid]
            || (matches!(dist_type, DistType::Dist13 | DistType::Dist14)
                && self.visited13_bounds[pid])
            || (matches!(dist_type, DistType::Dist14) && self.visited14_bounds[pid])
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

fn path14_id(nb: usize, bid1: usize, bid2: usize, bid3: usize) -> u64 {
    bid1 as u64 * nb as u64 * nb as u64 + bid2 as u64 * nb as u64 + bid3 as u64
}

fn ring_info_for_distgeom(mol: &Molecule) -> Result<Cow<'_, crate::RingInfo>, DgBoundsError> {
    if let Some(ring_info) = mol.derived_cache().rings.as_ref() {
        Ok(Cow::Borrowed(ring_info))
    } else {
        Ok(Cow::Owned(crate::symmetrize_sssr(mol)?))
    }
}

fn dg_generation_error(message: impl Into<String>) -> DgBoundsError {
    DgBoundsError::GenerationFailed(message.into())
}

fn bond_pair_shared_atom(
    mol: &Molecule,
    accum_data: &ComputedData,
    bid1: usize,
    bid2: usize,
) -> Result<usize, DgBoundsError> {
    let nb = mol.num_bonds();
    let aid = accum_data.get_bond_adj(nb, bid1, bid2);
    if aid < 0 {
        return Err(dg_generation_error(format!(
            "missing shared atom for bond pair ({bid1}, {bid2})"
        )));
    }
    let aid = aid as usize;
    if aid >= mol.num_atoms() {
        return Err(dg_generation_error(format!(
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
) -> Result<(), DgBoundsError> {
    if angle > 0.0 {
        Ok(())
    } else {
        Err(dg_generation_error(format!(
            "{context}: missing or invalid bond angle for bond pair ({bid1}, {bid2}): {angle}"
        )))
    }
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_getAtomStereo (BoundsMatrixBuilder.cpp:660-681)
// RDKit✔️✔️: Bond::BondStereo _getAtomStereo(const Bond *bnd, unsigned int aid1,
// RDKit✔️✔️:                                 unsigned int aid4) {
// RDKit✔️✔️:   auto stype = bnd->getStereo();
// RDKit✔️✔️:   if (stype > Bond::STEREOANY && bnd->getStereoAtoms().size() >= 2) {
// RDKit✔️✔️:     const auto &stAtoms = bnd->getStereoAtoms();
// RDKit✔️✔️:     if ((static_cast<unsigned int>(stAtoms[0]) != aid1) ^
// RDKit✔️✔️:         (static_cast<unsigned int>(stAtoms[1]) != aid4)) {
// RDKit✔️✔️:       switch (stype) {
// RDKit✔️✔️:         case Bond::STEREOZ:
// RDKit✔️✔️:           stype = Bond::STEREOE;
// RDKit✔️✔️:           break;
// RDKit✔️✔️:         case Bond::STEREOE:
// RDKit✔️✔️:           stype = Bond::STEREOZ;
// RDKit✔️✔️:           break;
// RDKit✔️✔️:         case Bond::STEREOCIS:
// RDKit✔️✔️:           stype = Bond::STEREOTRANS;
// RDKit✔️✔️:           break;
// RDKit✔️✔️:         case Bond::STEREOTRANS:
// RDKit✔️✔️:           stype = Bond::STEREOCIS;
// RDKit✔️✔️:           break;
// RDKit✔️✔️:         default:
// RDKit✔️✔️:           break;
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return stype;
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::_getAtomStereo
fn get_atom_stereo(bond: &Bond, aid1: usize, aid4: usize) -> BondStereo {
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

// ──────────────────────────────────────────────
// Set ring angle
// ──────────────────────────────────────────────

fn set_ring_angle(mol: &Molecule, aid2: usize, ring_size: usize) -> f64 {
    crate::source_port_helpers::rdkit_set_ring_angle(mol.atoms()[aid2].hybridization(), ring_size)
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::set13Bounds (BoundsMatrixBuilder.cpp:417-654)
// RDKit✔️❗: void set13Bounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat,
// RDKit✔️❗:                  ComputedData &accumData) {
// RDKit✔️❗:   auto npt = mmat->numRows();
// RDKit✔️❗:   CHECK_INVARIANT(npt == mol.getNumAtoms(), "Wrong size metric matrix");
// RDKit✔️❗:   CHECK_INVARIANT(accumData.bondAngles->numRows() == mol.getNumBonds(),
// RDKit✔️❗:                   "Wrong size bond angle matrix");
// RDKit✔️❗:   CHECK_INVARIANT(accumData.bondAdj->numRows() == mol.getNumBonds(),
// RDKit✔️❗:                   "Wrong size bond adjacency matrix");
// RDKit✔️❗:   // Since most of the special cases arise out of ring system, we will do
// RDKit✔️❗:   // the following here:
// RDKit✔️❗:   // - Loop over all the rings and set the 13 distances between atoms in
// RDKit✔️❗:   // these rings.
// RDKit✔️❗:   //   While doing this keep track of the ring atoms that have already been
// RDKit✔️❗:   //   used as the center atom.
// RDKit✔️❗:   // - Set the 13 distance between atoms that have a ring atom in between;
// RDKit✔️❗:   // these can be either non-ring atoms,
// RDKit✔️❗:   //   or a ring atom and a non-ring atom, or ring atoms that belong to
// RDKit✔️❗:   //   different simple rings
// RDKit✔️❗:   // - finally set all other 13 distances
// RDKit✔️❗:   const auto rinfo = mol.getRingInfo();
// RDKit✔️❗:   CHECK_INVARIANT(rinfo, "");
// RDKit✔️❗:   auto atomRings = rinfo->atomRings();
// RDKit✔️❗:   std::sort(atomRings.begin(), atomRings.end(), lessVector);
// RDKit✔️❗:   INT_VECT visited(npt, 0);
// RDKit✔️❗:   DOUBLE_VECT angleTaken(npt, 0.0);
// RDKit✔️❗:   auto nb = mol.getNumBonds();
// RDKit✔️❗:   BIT_SET donePaths(nb * nb);
// RDKit✔️❗:   // first deal with all rings and atoms in them
// RDKit✔️❗:   for (const auto &ringi : atomRings) { ... }
// RDKit✔️❗:   // now deal with the remaining atoms
// RDKit✔️❗:   for (aid2 = 0; aid2 < npt; aid2++) { ... }
// RDKit✔️❗: }
// END RDKIT CPP FUNCTION DGeomHelpers::set13Bounds
//
// Local review kept the second marker unresolved because this Rust port still
// performs repeated bond lookups through `bond_between_idx_simple()` inside the
// nested loops where RDKit iterates directly over bond iterators and matrix
// cells.
fn set_13_bounds(
    mol: &Molecule,
    mmat: &mut BoundsMatrix,
    accum_data: &mut ComputedData,
    rinfo: &crate::RingInfo,
) -> Result<(), DgBoundsError> {
    let npt = mmat.num_rows();
    if npt != mol.atoms().len() {
        return Err(dg_generation_error("Wrong size metric matrix"));
    }

    let nb = mol.bonds().len();
    if accum_data.bond_angles.len() != nb * nb {
        return Err(dg_generation_error("Wrong size bond angle matrix"));
    }
    if accum_data.bond_adj.len() != nb * nb {
        return Err(dg_generation_error("Wrong size bond adjacency matrix"));
    }

    let mut atom_rings: Vec<Vec<usize>> = rinfo
        .atom_rings()
        .iter()
        .map(|ring| ring.iter().map(|aid| aid.index()).collect())
        .collect();
    atom_rings.sort_by_key(Vec::len);

    let mut visited = vec![0usize; npt];
    let mut angle_taken = vec![0.0; npt];
    let mut done_paths = vec![false; nb * nb];

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
            let b1 = bond_between_idx_simple(mol, aid1, aid2).expect("no bond found");
            let b2 = bond_between_idx_simple(mol, aid2, aid3).expect("no bond found");
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
        let atom = &mol.atoms()[aid2];
        let nbrs = neighbors_for_atom(mol, aid2);
        let deg = nbrs.len();
        let n13 = deg * (deg.saturating_sub(1)) / 2;
        if n13 == visited[aid2] {
            continue;
        }
        let ahyb = atom.hybridization();

        if visited[aid2] >= 1 {
            for left in 0..nbrs.len() {
                let aid1 = nbrs[left];
                let bid1 = bond_between_idx_simple(mol, aid1, aid2).expect("missing bond");
                for right in 0..left {
                    let aid3 = nbrs[right];
                    let bid2 = bond_between_idx_simple(mol, aid3, aid2).expect("missing bond");
                    if accum_data.get_bond_angle(nb, bid1, bid2) >= 0.0 {
                        continue;
                    }

                    let angle = if ahyb == Hybridization::Sp2 {
                        (2.0 * std::f64::consts::PI - angle_taken[aid2])
                            / (n13 - visited[aid2]) as f64
                    } else if ahyb == Hybridization::Sp3 {
                        if rinfo.is_atom_in_ring_of_size(AtomId::new(aid2), 3) {
                            116.0_f64.to_radians()
                        } else if rinfo.is_atom_in_ring_of_size(AtomId::new(aid2), 4) {
                            112.0_f64.to_radians()
                        } else {
                            109.5_f64.to_radians()
                        }
                    } else if has_non_tetrahedral_stereo(atom) {
                        get_ideal_angle_between_ligands(aid2, aid1, aid3, mol).to_radians()
                    } else if deg == 5 {
                        105.0_f64.to_radians()
                    } else if deg == 6 {
                        135.0_f64.to_radians()
                    } else {
                        120.0_f64.to_radians()
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
                let aid1 = nbrs[left];
                let bid1 = bond_between_idx_simple(mol, aid1, aid2).expect("missing bond");
                for right in 0..left {
                    let aid3 = nbrs[right];
                    let bid2 = bond_between_idx_simple(mol, aid3, aid2).expect("missing bond");

                    let angle = if has_non_tetrahedral_stereo(atom) {
                        get_ideal_angle_between_ligands(aid2, aid1, aid3, mol).to_radians()
                    } else if ahyb == Hybridization::Sp {
                        std::f64::consts::PI
                    } else if ahyb == Hybridization::Sp2 {
                        2.0 * std::f64::consts::PI / 3.0
                    } else if ahyb == Hybridization::Sp3 {
                        109.5_f64.to_radians()
                    } else if ahyb == Hybridization::Sp3d {
                        105.0_f64.to_radians()
                    } else if ahyb == Hybridization::Sp3d2 {
                        135.0_f64.to_radians()
                    } else {
                        120.0_f64.to_radians()
                    };

                    let pid = aid1.min(aid3) * npt + aid1.max(aid3);
                    if !accum_data.visited_bound(pid, DistType::Dist12) {
                        if deg <= 4
                            || (has_non_tetrahedral_stereo(atom)
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
                            mmat.check_and_set_bounds(aid1, aid3, dl, du)?;
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

// ──────────────────────────────────────────────
// 1-4 distance bounds (torsion-based)
// ──────────────────────────────────────────────

fn total_num_hydrogens_rdkit_like(mol: &Molecule, atom_idx: usize) -> Option<u32> {
    let assignment = assign_valence(mol, ValenceModel::RdkitLike).ok()?;
    Some(total_num_hydrogens_for_distgeom(
        mol,
        &mol.atoms()[atom_idx],
        &assignment,
        atom_idx,
    ))
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_record14Path (BoundsMatrixBuilder.cpp:1319-1340)
// RDKit✔️❌: void _record14Path(const ROMol &mol, unsigned int bid1, unsigned int bid2,
// RDKit✔️❌:                    unsigned int bid3, ComputedData &accumData) {
// RDKit✔️❌:   const Atom *atm2 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid1, bid2));
// RDKit✔️❌:   PRECONDITION(atm2, "");
// RDKit✔️❌:   Atom::HybridizationType ahyb2 = atm2->getHybridization();
// RDKit✔️❌:   const Atom *atm3 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid2, bid3));
// RDKit✔️❌:   PRECONDITION(atm3, "");
// RDKit✔️❌:   Atom::HybridizationType ahyb3 = atm3->getHybridization();
// RDKit✔️❌:   unsigned int nb = mol.getNumBonds();
// RDKit✔️❌:   Path14Configuration path14;
// RDKit✔️❌:   path14.bid1 = bid1;
// RDKit✔️❌:   path14.bid2 = bid2;
// RDKit✔️❌:   path14.bid3 = bid3;
// RDKit✔️❌:   if ((ahyb2 == Atom::SP2) && (ahyb3 == Atom::SP2)) {  // FIX: check for trans
// RDKit✔️❌:     path14.type = Path14Configuration::CIS;
// RDKit✔️❌:     accumData.cisPaths.insert(static_cast<unsigned long>(bid1) * nb * nb +
// RDKit✔️❌:                               bid2 * nb + bid3);
// RDKit✔️❌:     accumData.cisPaths.insert(static_cast<unsigned long>(bid3) * nb * nb +
// RDKit✔️❌:                               bid2 * nb + bid1);
// RDKit✔️❌:   } else {
// RDKit✔️❌:     path14.type = Path14Configuration::OTHER;
// RDKit✔️❌:   }
// RDKit✔️❌:   accumData.paths14.push_back(path14);
// RDKit✔️❌: }
// END RDKIT CPP FUNCTION DGeomHelpers::_record14Path
//
// Performance note: the Rust path-flag store uses `Vec<u64>` with linear
// duplicate checks, so behavior matches RDKit but insertion cost is worse than
// RDKit's set-based storage.
fn record_14_path(
    mol: &Molecule,
    bid1: usize,
    bid2: usize,
    bid3: usize,
    accum_data: &mut ComputedData,
) -> Result<(), DgBoundsError> {
    let atm2 = bond_pair_shared_atom(mol, accum_data, bid1, bid2)?;
    let ahyb2 = mol.atoms()[atm2].hybridization();
    let atm3 = bond_pair_shared_atom(mol, accum_data, bid2, bid3)?;
    let ahyb3 = mol.atoms()[atm3].hybridization();
    let nb = mol.num_bonds();

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

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_setInRing14Bounds (BoundsMatrixBuilder.cpp:724-826)
// RDKit✔️❌: void _setInRing14Bounds(const ROMol &mol, const Bond *bnd1, const Bond *bnd2,
// RDKit✔️❌:                         const Bond *bnd3, ComputedData &accumData,
// RDKit✔️❌:                         DistGeom::BoundsMatPtr mmat, double *dmat,
// RDKit✔️❌:                         int ringSize) {
// RDKit✔️❌:   PRECONDITION(bnd1, "");
// RDKit✔️❌:   PRECONDITION(bnd2, "");
// RDKit✔️❌:   PRECONDITION(bnd3, "");
// RDKit✔️❌:   unsigned int bid1, bid2, bid3;
// RDKit✔️❌:   bid1 = bnd1->getIdx();
// RDKit✔️❌:   bid2 = bnd2->getIdx();
// RDKit✔️❌:   bid3 = bnd3->getIdx();
// RDKit✔️❌:   const Atom *atm2 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid1, bid2));
// RDKit✔️❌:   PRECONDITION(atm2, "");
// RDKit✔️❌:   Atom::HybridizationType ahyb2 = atm2->getHybridization();
// RDKit✔️❌:   const Atom *atm3 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid2, bid3));
// RDKit✔️❌:   PRECONDITION(atm3, "");
// RDKit✔️❌:   Atom::HybridizationType ahyb3 = atm3->getHybridization();
// RDKit✔️❌:   unsigned int aid1 = bnd1->getOtherAtomIdx(atm2->getIdx());
// RDKit✔️❌:   unsigned int aid4 = bnd3->getOtherAtomIdx(atm3->getIdx());
// RDKit✔️❌:   const unsigned int pid =
// RDKit✔️❌:       std::min(aid1, aid4) * mol.getNumAtoms() + std::max(aid1, aid4);
// RDKit✔️❌:   if (accumData.visitedBound(pid, DistType::DIST13)) {
// RDKit✔️❌:     return;
// RDKit✔️❌:   }
// RDKit✔️❌:   if (dmat[std::max(aid1, aid4) * mmat->numRows() + std::min(aid1, aid4)] < 2.9) {
// RDKit✔️❌:     return;
// RDKit✔️❌:   }
// RDKit✔️❌:   double bl1 = accumData.bondLengths[bid1];
// RDKit✔️❌:   double bl2 = accumData.bondLengths[bid2];
// RDKit✔️❌:   double bl3 = accumData.bondLengths[bid3];
// RDKit✔️❌:   double ba12 = accumData.bondAngles->getVal(bid1, bid2);
// RDKit✔️❌:   double ba23 = accumData.bondAngles->getVal(bid2, bid3);
// RDKit✔️❌:   CHECK_INVARIANT(ba12 > 0.0, "");
// RDKit✔️❌:   CHECK_INVARIANT(ba23 > 0.0, "");
// RDKit✔️❌:   double dl, du;
// RDKit✔️❌:   unsigned int nb = mol.getNumBonds();
// RDKit✔️❌:   Path14Configuration path14;
// RDKit✔️❌:   path14.bid1 = bid1;
// RDKit✔️❌:   path14.bid2 = bid2;
// RDKit✔️❌:   path14.bid3 = bid3;
// RDKit✔️❌:   Bond::BondStereo stype = _getAtomStereo(bnd2, aid1, aid4);
// RDKit✔️❌:   bool preferCis = false;
// RDKit✔️❌:   bool preferTrans = false;
// RDKit✔️❌:   if (ringSize <= 8 && (ahyb2 == Atom::SP2) && (ahyb3 == Atom::SP2) &&
// RDKit✔️❌:       (stype != Bond::STEREOE && stype != Bond::STEREOTRANS)) { ... }
// RDKit✔️❌:   else if (stype == Bond::STEREOZ || stype == Bond::STEREOCIS) { ... }
// RDKit✔️❌:   else if (stype == Bond::STEREOE || stype == Bond::STEREOTRANS) { ... }
// RDKit✔️❌:   if (preferCis) { ... cisPaths insert ... }
// RDKit✔️❌:   else if (preferTrans) { ... transPaths insert ... }
// RDKit✔️❌:   else { path14.type = Path14Configuration::OTHER; }
// RDKit✔️❌:   accumData.paths14.push_back(path14);
// RDKit✔️❌:   if (preferCis) { ... compute14DistCis +- GEN_DIST_TOL ... }
// RDKit✔️❌:   else if (preferTrans) { ... compute14DistTrans +- GEN_DIST_TOL ... }
// RDKit✔️❌:   else { ... cis/trans envelope with DIST12_DELTA widening ... }
// RDKit✔️❌:   accumData.visited14Bounds.set(pid);
// RDKit✔️❌:   _checkAndSetBounds(aid1, aid4, dl, du, mmat);
// RDKit✔️❌: }
// END RDKIT CPP FUNCTION DGeomHelpers::_setInRing14Bounds
//
// Performance note: behavior matches the RDKit source path, but the Rust
// implementation still uses linear path-flag storage and helper-based bond
// access, so this remains slower in hot paths.
fn set_in_ring_14_bounds(
    mol: &Molecule,
    bid1: usize,
    bid2: usize,
    bid3: usize,
    accum_data: &mut ComputedData,
    mmat: &mut BoundsMatrix,
    dmat: &[f64],
    ring_size: usize,
    ring_info: &crate::RingInfo,
) -> Result<(), DgBoundsError> {
    let atm2 = bond_pair_shared_atom(mol, accum_data, bid1, bid2)?;
    let ahyb2 = mol.atoms()[atm2].hybridization();
    let atm3 = bond_pair_shared_atom(mol, accum_data, bid2, bid3)?;
    let ahyb3 = mol.atoms()[atm3].hybridization();

    let bnd1 = &mol.bonds()[bid1];
    let bnd3 = &mol.bonds()[bid3];
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
    let pid = aid1.min(aid4) * mol.num_atoms() + aid1.max(aid4);
    if accum_data.visited_bound(pid, DistType::Dist13) {
        return Ok(());
    }
    if dmat[aid1.max(aid4) * mmat.num_rows() + aid1.min(aid4)] < 2.9 {
        return Ok(());
    }

    let bl1 = accum_data.bond_lengths[bid1];
    let bl2 = accum_data.bond_lengths[bid2];
    let bl3 = accum_data.bond_lengths[bid3];
    let ba12 = accum_data.get_bond_angle(mol.num_bonds(), bid1, bid2);
    let ba23 = accum_data.get_bond_angle(mol.num_bonds(), bid2, bid3);
    validate_bond_angle(ba12, bid1, bid2, "set_in_ring_14_bounds")?;
    validate_bond_angle(ba23, bid2, bid3, "set_in_ring_14_bounds")?;

    let stype = get_atom_stereo(&mol.bonds()[bid2], aid1, aid4);
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

    let nb = mol.num_bonds();
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

    let (mut dl, mut du) = if prefer_cis {
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
    mmat.check_and_set_bounds(aid1, aid4, dl, du)
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_setTwoInSameRing14Bounds (BoundsMatrixBuilder.cpp:828-900)
// RDKit✔️❌: void _setTwoInSameRing14Bounds(const ROMol &mol, const Bond *bnd1,
// RDKit✔️❌:                                const Bond *bnd2, const Bond *bnd3,
// RDKit✔️❌:                                ComputedData &accumData,
// RDKit✔️❌:                                DistGeom::BoundsMatPtr mmat, double *dmat) {
// RDKit✔️❌:   PRECONDITION(bnd1, "");
// RDKit✔️❌:   PRECONDITION(bnd2, "");
// RDKit✔️❌:   PRECONDITION(bnd3, "");
// RDKit✔️❌:   unsigned int bid1, bid2, bid3;
// RDKit✔️❌:   bid1 = bnd1->getIdx();
// RDKit✔️❌:   bid2 = bnd2->getIdx();
// RDKit✔️❌:   bid3 = bnd3->getIdx();
// RDKit✔️❌:   const Atom *atm2 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid1, bid2));
// RDKit✔️❌:   PRECONDITION(atm2, "");
// RDKit✔️❌:   const Atom *atm3 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid2, bid3));
// RDKit✔️❌:   PRECONDITION(atm3, "");
// RDKit✔️❌:   unsigned int aid1 = bnd1->getOtherAtomIdx(atm2->getIdx());
// RDKit✔️❌:   unsigned int aid4 = bnd3->getOtherAtomIdx(atm3->getIdx());
// RDKit✔️❌:   const unsigned int pid =
// RDKit✔️❌:       std::min(aid1, aid4) * mol.getNumAtoms() + std::max(aid1, aid4);
// RDKit✔️❌:   if (accumData.visitedBound(pid, DistType::DIST13)) { return; }
// RDKit✔️❌:   if (dmat[std::max(aid1, aid4) * mmat->numRows() + std::min(aid1, aid4)] < 2.9) {
// RDKit✔️❌:     return;
// RDKit✔️❌:   }
// RDKit✔️❌:   if (mol.getBondBetweenAtoms(aid1, atm3->getIdx()) ||
// RDKit✔️❌:       mol.getBondBetweenAtoms(aid4, atm2->getIdx())) {
// RDKit✔️❌:     return;
// RDKit✔️❌:   }
// RDKit✔️❌:   Atom::HybridizationType ahyb3 = atm3->getHybridization();
// RDKit✔️❌:   Atom::HybridizationType ahyb2 = atm2->getHybridization();
// RDKit✔️❌:   double bl1 = accumData.bondLengths[bid1];
// RDKit✔️❌:   double bl2 = accumData.bondLengths[bid2];
// RDKit✔️❌:   double bl3 = accumData.bondLengths[bid3];
// RDKit✔️❌:   double ba12 = accumData.bondAngles->getVal(bid1, bid2);
// RDKit✔️❌:   double ba23 = accumData.bondAngles->getVal(bid2, bid3);
// RDKit✔️❌:   CHECK_INVARIANT(ba12 > 0.0, "");
// RDKit✔️❌:   CHECK_INVARIANT(ba23 > 0.0, "");
// RDKit✔️❌:   double dl, du;
// RDKit✔️❌:   Path14Configuration path14;
// RDKit✔️❌:   unsigned int nb = mol.getNumBonds();
// RDKit✔️❌:   path14.bid1 = bid1;
// RDKit✔️❌:   path14.bid2 = bid2;
// RDKit✔️❌:   path14.bid3 = bid3;
// RDKit✔️❌:   if ((ahyb2 == Atom::SP2) && (ahyb3 == Atom::SP2)) { ... trans-only branch ... }
// RDKit✔️❌:   else { ... cis/trans envelope with strain swap and widening ... }
// RDKit✔️❌:   _checkAndSetBounds(aid1, aid4, dl, du, mmat);
// RDKit✔️❌:   accumData.paths14.push_back(path14);
// RDKit✔️❌:   accumData.visited14Bounds.set(pid);
// RDKit✔️❌: }
// END RDKIT CPP FUNCTION DGeomHelpers::_setTwoInSameRing14Bounds
fn set_two_in_same_ring_14_bounds(
    mol: &Molecule,
    bid1: usize,
    bid2: usize,
    bid3: usize,
    accum_data: &mut ComputedData,
    mmat: &mut BoundsMatrix,
    dmat: &[f64],
) -> Result<(), DgBoundsError> {
    let atm2 = bond_pair_shared_atom(mol, accum_data, bid1, bid2)?;
    let atm3 = bond_pair_shared_atom(mol, accum_data, bid2, bid3)?;
    let bnd1 = &mol.bonds()[bid1];
    let bnd3 = &mol.bonds()[bid3];
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
    let pid = aid1.min(aid4) * mol.num_atoms() + aid1.max(aid4);

    if accum_data.visited_bound(pid, DistType::Dist13) {
        return Ok(());
    }
    if dmat[aid1.max(aid4) * mmat.num_rows() + aid1.min(aid4)] < 2.9 {
        return Ok(());
    }
    if bond_between(mol, aid1, atm3).is_some() || bond_between(mol, aid4, atm2).is_some() {
        return Ok(());
    }

    let ahyb2 = mol.atoms()[atm2].hybridization();
    let ahyb3 = mol.atoms()[atm3].hybridization();
    let bl1 = accum_data.bond_lengths[bid1];
    let bl2 = accum_data.bond_lengths[bid2];
    let bl3 = accum_data.bond_lengths[bid3];
    let ba12 = accum_data.get_bond_angle(mol.num_bonds(), bid1, bid2);
    let ba23 = accum_data.get_bond_angle(mol.num_bonds(), bid2, bid3);
    validate_bond_angle(ba12, bid1, bid2, "set_two_in_same_ring_14_bounds")?;
    validate_bond_angle(ba23, bid2, bid3, "set_two_in_same_ring_14_bounds")?;

    let nb = mol.num_bonds();
    let (mut dl, mut du, kind) = if ahyb2 == Hybridization::Sp2 && ahyb3 == Hybridization::Sp2 {
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

    mmat.check_and_set_bounds(aid1, aid4, dl, du)?;
    accum_data.paths14.push(Path14Configuration {
        bid1,
        bid2,
        bid3,
        kind,
    });
    accum_data.visited14_bounds[pid] = true;
    Ok(())
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_setTwoInDiffRing14Bounds (BoundsMatrixBuilder.cpp:902-910)
// RDKit✔️✔️: void _setTwoInDiffRing14Bounds(const ROMol &mol, const Bond *bnd1,
// RDKit✔️✔️:                                const Bond *bnd2, const Bond *bnd3,
// RDKit✔️✔️:                                ComputedData &accumData,
// RDKit✔️✔️:                                DistGeom::BoundsMatPtr mmat, double *dmat) {
// RDKit✔️✔️:   // this turns out to be very similar to all bonds in the same ring
// RDKit✔️✔️:   // situation.
// RDKit✔️✔️:   // There is probably some fine tuning that can be done when the atoms a2
// RDKit✔️✔️:   // and a3 are not sp2 hybridized, but we will not worry about that now;
// RDKit✔️✔️:   // simple use 0-180 deg for non-sp2 cases.
// RDKit✔️✔️:   _setInRing14Bounds(mol, bnd1, bnd2, bnd3, accumData, mmat, dmat, 0);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::_setTwoInDiffRing14Bounds
fn set_two_in_diff_ring_14_bounds(
    mol: &Molecule,
    bid1: usize,
    bid2: usize,
    bid3: usize,
    accum_data: &mut ComputedData,
    mmat: &mut BoundsMatrix,
    dmat: &[f64],
    ring_info: &crate::RingInfo,
) -> Result<(), DgBoundsError> {
    set_in_ring_14_bounds(mol, bid1, bid2, bid3, accum_data, mmat, dmat, 0, ring_info)
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_setShareRingBond14Bounds (BoundsMatrixBuilder.cpp:912-920)
// RDKit✔️✔️: void _setShareRingBond14Bounds(const ROMol &mol, const Bond *bnd1,
// RDKit✔️✔️:                                const Bond *bnd2, const Bond *bnd3,
// RDKit✔️✔️:                                ComputedData &accumData,
// RDKit✔️✔️:                                DistGeom::BoundsMatPtr mmat, double *dmat) {
// RDKit✔️✔️:   // once this turns out to be similar to bonds in the same ring
// RDKit✔️✔️:   _setInRing14Bounds(mol, bnd1, bnd2, bnd3, accumData, mmat, dmat, 0);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::_setShareRingBond14Bounds
fn set_share_ring_bond_14_bounds(
    mol: &Molecule,
    bid1: usize,
    bid2: usize,
    bid3: usize,
    accum_data: &mut ComputedData,
    mmat: &mut BoundsMatrix,
    dmat: &[f64],
    ring_info: &crate::RingInfo,
) -> Result<(), DgBoundsError> {
    set_in_ring_14_bounds(mol, bid1, bid2, bid3, accum_data, mmat, dmat, 0, ring_info)
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_checkMacrocycleTwoInSameRingAmideEster14 (BoundsMatrixBuilder.cpp:1342-1361)
// RDKit✔️✔️: bool _checkMacrocycleTwoInSameRingAmideEster14(
// RDKit✔️✔️:     const Bond *bnd1, const Bond *bnd3, const Atom *atm1, const Atom *atm2,
// RDKit✔️✔️:     const Atom *atm3, const Atom *atm4) {
// RDKit✔️✔️:   unsigned int a1Num = atm1->getAtomicNum();
// RDKit✔️✔️:   unsigned int a2Num = atm2->getAtomicNum();
// RDKit✔️✔️:   unsigned int a3Num = atm3->getAtomicNum();
// RDKit✔️✔️:   unsigned int a4Num = atm4->getAtomicNum();
// RDKit✔️✔️:   return a1Num != 1 && a3Num == 6 && bnd3->getBondType() == Bond::DOUBLE &&
// RDKit✔️✔️:          (a4Num == 8 || a4Num == 7) && bnd1->getBondType() == Bond::SINGLE &&
// RDKit✔️✔️:          (a2Num == 8 || a2Num == 7);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::_checkMacrocycleTwoInSameRingAmideEster14
fn check_macrocycle_two_in_same_ring_amide_ester_14(
    mol: &Molecule,
    bnd1_idx: usize,
    bnd3_idx: usize,
    atm1: usize,
    atm2: usize,
    atm3: usize,
    atm4: usize,
) -> bool {
    let a1_num = mol.atoms()[atm1].atomic_number();
    let a2_num = mol.atoms()[atm2].atomic_number();
    let a3_num = mol.atoms()[atm3].atomic_number();
    let a4_num = mol.atoms()[atm4].atomic_number();

    a1_num != 1
        && a3_num == 6
        && mol.bonds()[bnd3_idx].order() == BondOrder::Double
        && (a4_num == 8 || a4_num == 7)
        && mol.bonds()[bnd1_idx].order() == BondOrder::Single
        && (a2_num == 8 || a2_num == 7)
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_setMacrocycleTwoInSameRing14Bounds (BoundsMatrixBuilder.cpp:1363-1454)
// RDKit✔️❌: void _setMacrocycleTwoInSameRing14Bounds(const ROMol &mol, const Bond *bnd1,
// RDKit✔️❌:                                          const Bond *bnd2, const Bond *bnd3,
// RDKit✔️❌:                                          ComputedData &accumData,
// RDKit✔️❌:                                          DistGeom::BoundsMatPtr mmat,
// RDKit✔️❌:                                          double *dmat) {
// RDKit✔️❌:   unsigned int bid1, bid2, bid3;
// RDKit✔️❌:   bid1 = bnd1->getIdx();
// RDKit✔️❌:   bid2 = bnd2->getIdx();
// RDKit✔️❌:   bid3 = bnd3->getIdx();
// RDKit✔️❌:   const Atom *atm2 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid1, bid2));
// RDKit✔️❌:   const Atom *atm3 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid2, bid3));
// RDKit✔️❌:   unsigned int aid1 = bnd1->getOtherAtomIdx(atm2->getIdx());
// RDKit✔️❌:   unsigned int aid4 = bnd3->getOtherAtomIdx(atm3->getIdx());
// RDKit✔️❌:   const Atom *atm1 = mol.getAtomWithIdx(aid1);
// RDKit✔️❌:   const Atom *atm4 = mol.getAtomWithIdx(aid4);
// RDKit✔️❌:   const unsigned int pid =
// RDKit✔️❌:       std::min(aid1, aid4) * mol.getNumAtoms() + std::max(aid1, aid4);
// RDKit✔️❌:   if (accumData.visitedBound(pid, DistType::DIST13)) {
// RDKit✔️❌:     // if this is already a 1-3 or 1-2 distance; do not overwrite
// RDKit✔️❌:     return;
// RDKit✔️❌:   }
// RDKit✔️❌:   if (dmat[std::max(aid1, aid4) * mmat->numRows() + std::min(aid1, aid4)] < 2.9) {
// RDKit✔️❌:     return;
// RDKit✔️❌:   }
// RDKit✔️❌:   if (mol.getBondBetweenAtoms(aid1, atm3->getIdx()) ||
// RDKit✔️❌:       mol.getBondBetweenAtoms(aid4, atm2->getIdx())) {
// RDKit✔️❌:     return;
// RDKit✔️❌:   }
// RDKit✔️❌:   double bl1 = accumData.bondLengths[bid1];
// RDKit✔️❌:   double bl2 = accumData.bondLengths[bid2];
// RDKit✔️❌:   double bl3 = accumData.bondLengths[bid3];
// RDKit✔️❌:   double ba12 = accumData.bondAngles->getVal(bid1, bid2);
// RDKit✔️❌:   double ba23 = accumData.bondAngles->getVal(bid2, bid3);
// RDKit✔️❌:   Path14Configuration path14;
// RDKit✔️❌:   path14.bid1 = bid1;
// RDKit✔️❌:   path14.bid2 = bid2;
// RDKit✔️❌:   path14.bid3 = bid3;
// RDKit✔️❌:   if ((_checkMacrocycleTwoInSameRingAmideEster14(bnd1, bnd3, atm1, atm2, atm3,
// RDKit✔️❌:                                                  atm4)) ||
// RDKit✔️❌:       (_checkMacrocycleTwoInSameRingAmideEster14(bnd3, bnd1, atm4, atm3, atm2,
// RDKit✔️❌:                                                  atm1))) {
// RDKit✔️❌:     dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23);
// RDKit✔️❌:     path14.type = Path14Configuration::CIS;
// RDKit✔️❌:   } else {
// RDKit✔️❌:     // here we will assume anything is possible
// RDKit✔️❌:     dl = RDGeom::compute14DistCis(bl1, bl2, bl3, ba12, ba23);
// RDKit✔️❌:     du = RDGeom::compute14DistTrans(bl1, bl2, bl3, ba12, ba23);
// RDKit✔️❌:     // in highly-strained situations these can get mixed up:
// RDKit✔️❌:     if (du < dl) {
// RDKit✔️❌:       double tmpD = dl;
// RDKit✔️❌:       dl = du;
// RDKit✔️❌:       du = tmpD;
// RDKit✔️❌:     }
// RDKit✔️❌:     path14.type = Path14Configuration::OTHER;
// RDKit✔️❌:   }
// RDKit✔️❌:   _checkAndSetBounds(aid1, aid4, dl, du, mmat);
// RDKit✔️❌:   accumData.paths14.push_back(path14);
// RDKit✔️❌:   accumData.visited14Bounds.set(pid);
// RDKit✔️❌: }
// END RDKIT CPP FUNCTION DGeomHelpers::_setMacrocycleTwoInSameRing14Bounds
fn set_macrocycle_two_in_same_ring_14_bounds(
    mol: &Molecule,
    bid1: usize,
    bid2: usize,
    bid3: usize,
    accum_data: &mut ComputedData,
    mmat: &mut BoundsMatrix,
    dmat: &[f64],
) -> Result<(), DgBoundsError> {
    // Behavior matches RDKit, but duplicate path recording and bond/neighbor lookup
    // still go through Vec-backed helpers instead of RDKit's set and adjacency types.
    let atm2 = bond_pair_shared_atom(mol, accum_data, bid1, bid2)?;
    let atm3 = bond_pair_shared_atom(mol, accum_data, bid2, bid3)?;
    let bnd1 = &mol.bonds()[bid1];
    let bnd3 = &mol.bonds()[bid3];
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
    let pid = aid1.min(aid4) * mol.num_atoms() + aid1.max(aid4);

    if accum_data.visited_bound(pid, DistType::Dist13) {
        return Ok(());
    }
    if dmat[aid1.max(aid4) * mmat.num_rows() + aid1.min(aid4)] < 2.9 {
        return Ok(());
    }
    if bond_between(mol, aid1, atm3).is_some() || bond_between(mol, aid4, atm2).is_some() {
        return Ok(());
    }

    let bl1 = accum_data.bond_lengths[bid1];
    let bl2 = accum_data.bond_lengths[bid2];
    let bl3 = accum_data.bond_lengths[bid3];
    let ba12 = accum_data.get_bond_angle(mol.num_bonds(), bid1, bid2);
    let ba23 = accum_data.get_bond_angle(mol.num_bonds(), bid2, bid3);
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

    let nb = mol.num_bonds();
    let (mut dl, mut du, kind) = if check_macrocycle_two_in_same_ring_amide_ester_14(
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

    mmat.check_and_set_bounds(aid1, aid4, dl, du)?;
    accum_data.paths14.push(Path14Configuration {
        bid1,
        bid2,
        bid3,
        kind,
    });
    accum_data.visited14_bounds[pid] = true;
    Ok(())
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_setMacrocycleAllInSameRing14Bounds (BoundsMatrixBuilder.cpp:1456-1660)
// RDKit✔️❌: void _setMacrocycleAllInSameRing14Bounds(const ROMol &mol, const Bond *bnd1,
// RDKit✔️❌:                                          const Bond *bnd2, const Bond *bnd3,
// RDKit✔️❌:                                          ComputedData &accumData,
// RDKit✔️❌:                                          DistGeom::BoundsMatPtr mmat,
// RDKit✔️❌:                                          double *) {
// RDKit✔️❌:   // This is adapted from `_setChain14Bounds`, with changes on how trans amide
// RDKit✔️❌:   // is handled
// RDKit✔️❌:   unsigned int bid1, bid2, bid3;
// RDKit✔️❌:   bid1 = bnd1->getIdx();
// RDKit✔️❌:   bid2 = bnd2->getIdx();
// RDKit✔️❌:   bid3 = bnd3->getIdx();
// RDKit✔️❌:   const Atom *atm2 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid1, bid2));
// RDKit✔️❌:   const Atom *atm3 = mol.getAtomWithIdx(accumData.bondAdj->getVal(bid2, bid3));
// RDKit✔️❌:   unsigned int aid1 = bnd1->getOtherAtomIdx(atm2->getIdx());
// RDKit✔️❌:   unsigned int aid4 = bnd3->getOtherAtomIdx(atm3->getIdx());
// RDKit✔️❌:   const unsigned int pid =
// RDKit✔️❌:       std::min(aid1, aid4) * mol.getNumAtoms() + std::max(aid1, aid4);
// RDKit✔️❌:   if (accumData.visitedBound(pid, DistType::DIST13)) {
// RDKit✔️❌:     // if this is already a 1-3 or 1-2 distance; do not overwrite
// RDKit✔️❌:     return;
// RDKit✔️❌:   }
// RDKit✔️❌:   bool setTheBound = true;
// RDKit✔️❌:   // if the middle bond is double
// RDKit✔️❌:   switch (bnd2->getBondType()) {
// RDKit✔️❌:     case Bond::DOUBLE:
// RDKit✔️❌:       // if any of the other bonds are double - the torsion angle is zero
// RDKit✔️❌:       // this is CC=C=C situation
// RDKit✔️❌:       if ((bnd1->getBondType() == Bond::DOUBLE) ||
// RDKit✔️❌:           (bnd3->getBondType() == Bond::DOUBLE)) {
// RDKit✔️❌:         path14.type = Path14Configuration::CIS;
// RDKit✔️❌:       } else if (bnd2->getStereo() > Bond::STEREOANY) {
// RDKit✔️❌:         Bond::BondStereo stype = _getAtomStereo(bnd2, aid1, aid4);
// RDKit✔️❌:         if (stype == Bond::STEREOZ || stype == Bond::STEREOCIS) {
// RDKit✔️❌:           path14.type = Path14Configuration::CIS;
// RDKit✔️❌:         } else {
// RDKit✔️❌:           path14.type = Path14Configuration::TRANS;
// RDKit✔️❌:         }
// RDKit✔️❌:       } else {
// RDKit✔️❌:         // double bond with no stereo setting can be 0 or 180
// RDKit✔️❌:         path14.type = Path14Configuration::OTHER;
// RDKit✔️❌:       }
// RDKit✔️❌:       break;
// RDKit✔️❌:     case Bond::SINGLE:
// RDKit✔️❌:       if ((atm2->getAtomicNum() == 16) && (atm3->getAtomicNum() == 16) &&
// RDKit✔️❌:           (atm2->getDegree() == 2) && (atm3->getDegree() == 2)) {
// RDKit✔️❌:         path14.type = Path14Configuration::OTHER;
// RDKit✔️❌:       } else if ((_checkMacrocycleAllInSameRingAmideEster14(
// RDKit✔️❌:                      mol, bnd1, bnd3, atm1, atm2, atm3, atm4)) ||
// RDKit✔️❌:                  (_checkMacrocycleAllInSameRingAmideEster14(
// RDKit✔️❌:                      mol, bnd3, bnd1, atm4, atm3, atm2, atm1))) {
// RDKit✔️❌:         // we saw that the currently defined max distance for trans
// RDKit✔️❌:         // is still a bit too short, thus we add an additional 0.1,
// RDKit✔️❌:         // which is the max that works without triangular smoothing error
// RDKit✔️❌:         path14.type = Path14Configuration::TRANS;
// RDKit✔️❌:       } else if ((_checkAmideEster15(mol, bnd1, bnd3, atm1, atm2, atm3, atm4)) ||
// RDKit✔️❌:                  (_checkAmideEster15(mol, bnd3, bnd1, atm4, atm3, atm2, atm1))) {
// RDKit✔️❌:         // amide is cis, we're trans:
// RDKit✔️❌:         if (atm2->getAtomicNum() == 7 && atm2->getDegree() == 3 &&
// RDKit✔️❌:             atm1->getAtomicNum() == 1 && atm2->getTotalNumHs(true) == 1) {
// RDKit✔️❌:           // secondary amide, this is the H
// RDKit✔️❌:           setTheBound = false;
// RDKit✔️❌:         } else {
// RDKit✔️❌:           path14.type = Path14Configuration::TRANS;
// RDKit✔️❌:         }
// RDKit✔️❌:       } else {
// RDKit✔️❌:         path14.type = Path14Configuration::OTHER;
// RDKit✔️❌:       }
// RDKit✔️❌:       break;
// RDKit✔️❌:     default:
// RDKit✔️❌:       path14.type = Path14Configuration::OTHER;
// RDKit✔️❌:   }
// RDKit✔️❌:   if (setTheBound) {
// RDKit✔️❌:     _checkAndSetBounds(aid1, aid4, dl, du, mmat);
// RDKit✔️❌:     accumData.paths14.push_back(path14);
// RDKit✔️❌:     accumData.visited14Bounds.set(pid);
// RDKit✔️❌:   }
// RDKit✔️❌: }
// END RDKIT CPP FUNCTION DGeomHelpers::_setMacrocycleAllInSameRing14Bounds
fn set_macrocycle_all_in_same_ring_14_bounds(
    mol: &Molecule,
    bid1: usize,
    bid2: usize,
    bid3: usize,
    accum_data: &mut ComputedData,
    mmat: &mut BoundsMatrix,
) -> Result<(), DgBoundsError> {
    // Behavior matches RDKit's branch structure, including the +0.1 trans-amide
    // widening and secondary-amide-H suppression branch, but helper lookups are
    // still repeated over Vec-backed adjacency data.
    let atm2 = bond_pair_shared_atom(mol, accum_data, bid1, bid2)?;
    let atm3 = bond_pair_shared_atom(mol, accum_data, bid2, bid3)?;
    let bnd1 = &mol.bonds()[bid1];
    let bnd2 = &mol.bonds()[bid2];
    let bnd3 = &mol.bonds()[bid3];
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
    let pid = aid1.min(aid4) * mol.num_atoms() + aid1.max(aid4);
    if accum_data.visited_bound(pid, DistType::Dist13) {
        return Ok(());
    }

    let atm1 = aid1;
    let atm4 = aid4;
    let bl1 = accum_data.bond_lengths[bid1];
    let bl2 = accum_data.bond_lengths[bid2];
    let bl3 = accum_data.bond_lengths[bid3];
    let ba12 = accum_data.get_bond_angle(mol.num_bonds(), bid1, bid2);
    let ba23 = accum_data.get_bond_angle(mol.num_bonds(), bid2, bid3);
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
    let nb = mol.num_bonds();
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
            if mol.atoms()[atm2].atomic_number() == 16
                && mol.atoms()[atm3].atomic_number() == 16
                && neighbors_for_atom(mol, atm2).len() == 2
                && neighbors_for_atom(mol, atm3).len() == 2
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
            } else if check_amide_ester_15(mol, bid1, bid3, atm2, atm3)
                || check_amide_ester_15(mol, bid3, bid1, atm3, atm2)
            {
                // COSMolKit models the default RDKit branch where FORCE_TRANS_AMIDES is not defined.
                let total_hs_atm2 = total_num_hydrogens_rdkit_like(mol, atm2).unwrap_or(0);
                if mol.atoms()[atm2].atomic_number() == 7
                    && neighbors_for_atom(mol, atm2).len() == 3
                    && mol.atoms()[atm1].atomic_number() == 1
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
        mmat.check_and_set_bounds(aid1, aid4, dl, du)?;
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

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_setChain14Bounds (BoundsMatrixBuilder.cpp:980-1317)
// RDKit✔️❌: void _setChain14Bounds(const ROMol &mol, const Bond *bnd1, const Bond *bnd2,
// RDKit✔️❌:                        const Bond *bnd3, ComputedData &accumData,
// RDKit✔️❌:                        DistGeom::BoundsMatPtr mmat, double *,
// RDKit✔️❌:                        bool forceTransAmides) {
// RDKit✔️❌:   const unsigned int pid =
// RDKit✔️❌:       std::min(aid1, aid4) * mol.getNumAtoms() + std::max(aid1, aid4);
// RDKit✔️❌:   if (accumData.visitedBound(pid, DistType::DIST13)) {
// RDKit✔️❌:     return;
// RDKit✔️❌:   }
// RDKit✔️❌:   bool setTheBound = true;
// RDKit✔️❌:   switch (bnd2->getBondType()) {
// RDKit✔️❌:     case Bond::DOUBLE:
// RDKit✔️❌:       // if any of the other bonds are double - the torsion angle is zero
// RDKit✔️❌:       // this is CC=C=C situation
// RDKit✔️❌:       if ((bnd1->getBondType() == Bond::DOUBLE) ||
// RDKit✔️❌:           (bnd3->getBondType() == Bond::DOUBLE)) {
// RDKit✔️❌:         path14.type = Path14Configuration::CIS;
// RDKit✔️❌:       } else if (bnd2->getStereo() > Bond::STEREOANY) {
// RDKit✔️❌:         Bond::BondStereo stype = _getAtomStereo(bnd2, aid1, aid4);
// RDKit✔️❌:         if (stype == Bond::STEREOZ || stype == Bond::STEREOCIS) {
// RDKit✔️❌:           path14.type = Path14Configuration::CIS;
// RDKit✔️❌:         } else {
// RDKit✔️❌:           path14.type = Path14Configuration::TRANS;
// RDKit✔️❌:         }
// RDKit✔️❌:       } else {
// RDKit✔️❌:         // double bond with no stereo setting can be 0 or 180
// RDKit✔️❌:         path14.type = Path14Configuration::OTHER;
// RDKit✔️❌:       }
// RDKit✔️❌:       break;
// RDKit✔️❌:     case Bond::SINGLE:
// RDKit✔️❌:       if ((atm2->getAtomicNum() == 16) && (atm3->getAtomicNum() == 16) &&
// RDKit✔️❌:           (atm2->getDegree() == 2) && (atm3->getDegree() == 2)) {
// RDKit✔️❌:         // this is *S-S* situation
// RDKit✔️❌:         path14.type = Path14Configuration::OTHER;
// RDKit✔️❌:       } else if ((_checkAmideEster14(bnd1, bnd3, atm1, atm2, atm3, atm4)) ||
// RDKit✔️❌:                  (_checkAmideEster14(bnd3, bnd1, atm4, atm3, atm2, atm1))) {
// RDKit✔️❌:         if (forceTransAmides) {
// RDKit✔️❌:           // secondary amide H is trans to O; otherwise use cis
// RDKit✔️❌:         } else {
// RDKit✔️❌:           path14.type = Path14Configuration::OTHER;
// RDKit✔️❌:         }
// RDKit✔️❌:       } else if ((_checkAmideEster15(mol, bnd1, bnd3, atm1, atm2, atm3, atm4)) ||
// RDKit✔️❌:                  (_checkAmideEster15(mol, bnd3, bnd1, atm4, atm3, atm2, atm1))) {
// RDKit✔️❌:         if (forceTransAmides) {
// RDKit✔️❌:           // secondary amide H is cis to atom 5; otherwise use trans
// RDKit✔️❌:         } else {
// RDKit✔️❌:           path14.type = Path14Configuration::OTHER;
// RDKit✔️❌:         }
// RDKit✔️❌:       } else {
// RDKit✔️❌:         path14.type = Path14Configuration::OTHER;
// RDKit✔️❌:       }
// RDKit✔️❌:       break;
// RDKit✔️❌:     default:
// RDKit✔️❌:       path14.type = Path14Configuration::OTHER;
// RDKit✔️❌:   }
// RDKit✔️❌:   if (setTheBound) {
// RDKit✔️❌:     _checkAndSetBounds(aid1, aid4, dl, du, mmat);
// RDKit✔️❌:     accumData.paths14.push_back(path14);
// RDKit✔️❌:     accumData.visited14Bounds.set(pid);
// RDKit✔️❌:   }
// RDKit✔️❌: }
// END RDKIT CPP FUNCTION DGeomHelpers::_setChain14Bounds
fn set_chain_14_bounds(
    mol: &Molecule,
    bid1: usize,
    bid2: usize,
    bid3: usize,
    accum_data: &mut ComputedData,
    mmat: &mut BoundsMatrix,
    force_trans_amides: bool,
) -> Result<(), DgBoundsError> {
    // Behavior matches RDKit's branch matrix, but path flags use Vec dedup scans
    // instead of RDKit set containers and several helper predicates rescan local
    // adjacency/state on each branch.
    let atm2 = bond_pair_shared_atom(mol, accum_data, bid1, bid2)?;
    let atm3 = bond_pair_shared_atom(mol, accum_data, bid2, bid3)?;
    let bnd1 = &mol.bonds()[bid1];
    let bnd2 = &mol.bonds()[bid2];
    let bnd3 = &mol.bonds()[bid3];
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
    let pid = aid1.min(aid4) * mol.num_atoms() + aid1.max(aid4);
    if accum_data.visited_bound(pid, DistType::Dist13) {
        return Ok(());
    }

    let atm1 = aid1;
    let atm4 = aid4;
    let bl1 = accum_data.bond_lengths[bid1];
    let bl2 = accum_data.bond_lengths[bid2];
    let bl3 = accum_data.bond_lengths[bid3];
    let ba12 = accum_data.get_bond_angle(mol.num_bonds(), bid1, bid2);
    let ba23 = accum_data.get_bond_angle(mol.num_bonds(), bid2, bid3);
    validate_bond_angle(ba12, bid1, bid2, "set_chain_14_bounds")?;
    validate_bond_angle(ba23, bid2, bid3, "set_chain_14_bounds")?;

    let mut set_the_bound = true;
    let nb = mol.num_bonds();
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
            if mol.atoms()[atm2].atomic_number() == 16
                && mol.atoms()[atm3].atomic_number() == 16
                && neighbors_for_atom(mol, atm2).len() == 2
                && neighbors_for_atom(mol, atm3).len() == 2
            {
                let dl = compute_14_dist_3d(bl1, bl2, bl3, ba12, ba23, std::f64::consts::PI / 2.0)
                    - GEN_DIST_TOL;
                (dl, dl + 2.0 * GEN_DIST_TOL, Path14Kind::Other)
            } else if check_amide_ester_14(mol, bid1, bid3, atm2, atm3, atm4)
                || check_amide_ester_14(mol, bid3, bid1, atm3, atm2, atm1)
            {
                if force_trans_amides {
                    let total_hs_atm2 = total_num_hydrogens_rdkit_like(mol, atm2).unwrap_or(0);
                    let total_hs_atm3 = total_num_hydrogens_rdkit_like(mol, atm3).unwrap_or(0);
                    let secondary_left = mol.atoms()[atm1].atomic_number() == 1
                        && mol.atoms()[atm2].atomic_number() == 7
                        && neighbors_for_atom(mol, atm2).len() == 3
                        && total_hs_atm2 == 1;
                    let secondary_right = mol.atoms()[atm4].atomic_number() == 1
                        && mol.atoms()[atm3].atomic_number() == 7
                        && neighbors_for_atom(mol, atm3).len() == 3
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
            } else if check_amide_ester_15(mol, bid1, bid3, atm2, atm3)
                || check_amide_ester_15(mol, bid3, bid1, atm3, atm2)
            {
                if force_trans_amides {
                    let total_hs_atm2 = total_num_hydrogens_rdkit_like(mol, atm2).unwrap_or(0);
                    let total_hs_atm3 = total_num_hydrogens_rdkit_like(mol, atm3).unwrap_or(0);
                    let secondary_left = mol.atoms()[atm1].atomic_number() == 1
                        && mol.atoms()[atm2].atomic_number() == 7
                        && neighbors_for_atom(mol, atm2).len() == 3
                        && total_hs_atm2 == 1;
                    let secondary_right = mol.atoms()[atm4].atomic_number() == 1
                        && mol.atoms()[atm3].atomic_number() == 7
                        && neighbors_for_atom(mol, atm3).len() == 3
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
        mmat.check_and_set_bounds(aid1, aid4, dl, du)?;
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

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_checkH2NX3H1OX2 (BoundsMatrixBuilder.cpp:845-856)
// RDKit✔️❌: bool _checkH2NX3H1OX2(const Atom *atm) {
// RDKit✔️❌:   if ((atm->getAtomicNum() == 6) && (atm->getTotalNumHs(true) == 2)) {
// RDKit✔️❌:     // CH2
// RDKit✔️❌:     return true;
// RDKit✔️❌:   } else if ((atm->getAtomicNum() == 8) && (atm->getTotalNumHs(true) == 0)) {
// RDKit✔️❌:     // OX2
// RDKit✔️❌:     return true;
// RDKit✔️❌:   } else if ((atm->getAtomicNum() == 7) && (atm->getDegree() == 3) &&
// RDKit✔️❌:              (atm->getTotalNumHs(true) == 1)) {
// RDKit✔️❌:     // FIX: assuming hydrogen is not in the graph
// RDKit✔️❌:     // this is NX3H1 situation
// RDKit✔️❌:     return true;
// RDKit✔️❌:   }
// RDKit✔️❌:   return false;
// RDKit✔️❌: }
// END RDKIT CPP FUNCTION DGeomHelpers::_checkH2NX3H1OX2
fn check_h2_nx3_h1_ox2(mol: &Molecule, atom_idx: usize) -> bool {
    let atom = &mol.atoms()[atom_idx];
    // Behavior matches RDKit, but this helper recomputes total H count and neighbor
    // degree through Rust-side scans instead of RDKit's cached atom state.
    let Some(total_hs) = total_num_hydrogens_rdkit_like(mol, atom_idx) else {
        return false;
    };
    (atom.atomic_number() == 6 && total_hs == 2)
        || (atom.atomic_number() == 8 && total_hs == 0)
        || (atom.atomic_number() == 7
            && neighbors_for_atom(mol, atom_idx).len() == 3
            && total_hs == 1)
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_checkNhChChNh (BoundsMatrixBuilder.cpp:858-866)
// RDKit✔️❌: bool _checkNhChChNh(const Atom *atm1, const Atom *atm2, const Atom *atm3,
// RDKit✔️❌:                     const Atom *atm4) {
// RDKit✔️❌:   // checking for [!#1]~$ch!@$ch~[!#1], where ch = [CH2,NX3H1,OX2] situation
// RDKit✔️❌:   if ((atm1->getAtomicNum() != 1) && (atm4->getAtomicNum() != 1)) {
// RDKit✔️❌:     // end atom not hydrogens
// RDKit✔️❌:     if ((_checkH2NX3H1OX2(atm2)) && (_checkH2NX3H1OX2(atm3))) {
// RDKit✔️❌:       return true;
// RDKit✔️❌:     }
// RDKit✔️❌:   }
// RDKit✔️❌:   return false;
// RDKit✔️❌: }
// END RDKIT CPP FUNCTION DGeomHelpers::_checkNhChChNh
fn check_nh_ch_ch_nh(mol: &Molecule, atm1: usize, atm2: usize, atm3: usize, atm4: usize) -> bool {
    mol.atoms()[atm1].atomic_number() != 1
        && mol.atoms()[atm4].atomic_number() != 1
        && check_h2_nx3_h1_ox2(mol, atm2)
        && check_h2_nx3_h1_ox2(mol, atm3)
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_checkAmideEster14 (BoundsMatrixBuilder.cpp:874-894)
// RDKit✔️✔️: bool _checkAmideEster14(const Bond *bnd1, const Bond *bnd3, const Atom *,
// RDKit✔️✔️:                         const Atom *atm2, const Atom *atm3, const Atom *atm4) {
// RDKit✔️✔️:   unsigned int a2Num = atm2->getAtomicNum();
// RDKit✔️✔️:   unsigned int a3Num = atm3->getAtomicNum();
// RDKit✔️✔️:   unsigned int a4Num = atm4->getAtomicNum();
// RDKit✔️✔️:   if (a3Num == 6 && bnd3->getBondType() == Bond::DOUBLE &&
// RDKit✔️✔️:       (a4Num == 8 || a4Num == 7) && bnd1->getBondType() == Bond::SINGLE &&
// RDKit✔️✔️:       (a2Num == 8 || (a2Num == 7 && atm2->getTotalNumHs(true) == 1))) {
// RDKit✔️✔️:     return true;
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return false;
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::_checkAmideEster14
fn check_amide_ester_14(
    mol: &Molecule,
    bnd1_idx: usize,
    bnd3_idx: usize,
    atm2: usize,
    atm3: usize,
    atm4: usize,
) -> bool {
    let bnd1 = &mol.bonds()[bnd1_idx];
    let bnd3 = &mol.bonds()[bnd3_idx];
    let a2_num = mol.atoms()[atm2].atomic_number();
    let a3_num = mol.atoms()[atm3].atomic_number();
    let a4_num = mol.atoms()[atm4].atomic_number();
    let total_hs_atm2 = total_num_hydrogens_rdkit_like(mol, atm2).unwrap_or(0);

    a3_num == 6
        && bnd3.order() == BondOrder::Double
        && (a4_num == 8 || a4_num == 7)
        && bnd1.order() == BondOrder::Single
        && (a2_num == 8 || (a2_num == 7 && total_hs_atm2 == 1))
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_checkMacrocycleAllInSameRingAmideEster14 (BoundsMatrixBuilder.cpp:903-949)
// RDKit✔️❌: bool _checkMacrocycleAllInSameRingAmideEster14(const ROMol &mol, const Bond *,
// RDKit✔️❌:                                                const Bond *, const Atom *atm1,
// RDKit✔️❌:                                                const Atom *atm2,
// RDKit✔️❌:                                                const Atom *atm3,
// RDKit✔️❌:                                                const Atom *atm4) {
// RDKit✔️❌:   //   This is a re-write of `_checkAmideEster14` with more explicit logic
// RDKit✔️❌:   //   on the checks It is interesting that we find with this function we
// RDKit✔️❌:   //   get better macrocycle sampling than `_checkAmideEster14`
// RDKit✔️❌:   unsigned int a2Num = atm2->getAtomicNum();
// RDKit✔️❌:   unsigned int a3Num = atm3->getAtomicNum();
// RDKit✔️❌:
// RDKit✔️❌:   if (a3Num != 6) {
// RDKit✔️❌:     return false;
// RDKit✔️❌:   }
// RDKit✔️❌:
// RDKit✔️❌:   if (a2Num == 7 || a2Num == 8) {
// RDKit✔️❌:     if (mol.getAtomDegree(atm2) == 3 && mol.getAtomDegree(atm3) == 3) {
// RDKit✔️❌:       for (auto nbrIdx :
// RDKit✔️❌:            boost::make_iterator_range(mol.getAtomNeighbors(atm2))) {
// RDKit✔️❌:         if (nbrIdx != atm1->getIdx() && nbrIdx != atm3->getIdx()) {
// RDKit✔️❌:           const auto &res = mol.getAtomWithIdx(nbrIdx);
// RDKit✔️❌:           const auto &resbnd = mol.getBondBetweenAtoms(atm2->getIdx(), nbrIdx);
// RDKit✔️❌:           if ((res->getAtomicNum() != 6 &&
// RDKit✔️❌:                res->getAtomicNum() != 1) ||  // check is (methylated)amide
// RDKit✔️❌:               resbnd->getBondType() != Bond::SINGLE) {
// RDKit✔️❌:             return false;
// RDKit✔️❌:           }
// RDKit✔️❌:           break;
// RDKit✔️❌:         }
// RDKit✔️❌:       }
// RDKit✔️❌:
// RDKit✔️❌:       for (auto nbrIdx :
// RDKit✔️❌:            boost::make_iterator_range(mol.getAtomNeighbors(atm3))) {
// RDKit✔️❌:         if (nbrIdx != atm2->getIdx() && nbrIdx != atm4->getIdx()) {
// RDKit✔️❌:           const auto &res = mol.getAtomWithIdx(nbrIdx);
// RDKit✔️❌:           const auto &resbnd = mol.getBondBetweenAtoms(atm3->getIdx(), nbrIdx);
// RDKit✔️❌:           if (res->getAtomicNum() != 8 ||  // check for the carbonyl oxygen
// RDKit✔️❌:               resbnd->getBondType() != Bond::DOUBLE) {
// RDKit✔️❌:             return false;
// RDKit✔️❌:           }
// RDKit✔️❌:           break;
// RDKit✔️❌:         }
// RDKit✔️❌:       }
// RDKit✔️❌:
// RDKit✔️❌:       return true;
// RDKit✔️❌:     }
// RDKit✔️❌:   }
// RDKit✔️❌:   return false;
// RDKit✔️❌: }
// END RDKIT CPP FUNCTION DGeomHelpers::_checkMacrocycleAllInSameRingAmideEster14
fn check_macrocycle_all_in_same_ring_amide_ester_14(
    mol: &Molecule,
    atm1: usize,
    atm2: usize,
    atm3: usize,
    atm4: usize,
) -> bool {
    let a2_num = mol.atoms()[atm2].atomic_number();
    let a3_num = mol.atoms()[atm3].atomic_number();
    // Behavior matches RDKit's explicit neighbor scans, but our helper still
    // performs repeated adjacency and bond lookups over Vec-backed structures.
    if a3_num != 6 {
        return false;
    }
    if a2_num != 7 && a2_num != 8 {
        return false;
    }
    if neighbors_for_atom(mol, atm2).len() != 3 || neighbors_for_atom(mol, atm3).len() != 3 {
        return false;
    }

    for nbr_idx in neighbors_for_atom(mol, atm2) {
        if nbr_idx != atm1 && nbr_idx != atm3 {
            let res = &mol.atoms()[nbr_idx];
            let res_bnd = &mol.bonds()[bond_between_idx_simple(mol, atm2, nbr_idx).expect("bond")];
            if (res.atomic_number() != 6 && res.atomic_number() != 1)
                || res_bnd.order() != BondOrder::Single
            {
                return false;
            }
            break;
        }
    }

    for nbr_idx in neighbors_for_atom(mol, atm3) {
        if nbr_idx != atm2 && nbr_idx != atm4 {
            let res = &mol.atoms()[nbr_idx];
            let res_bnd = &mol.bonds()[bond_between_idx_simple(mol, atm3, nbr_idx).expect("bond")];
            if res.atomic_number() != 8 || res_bnd.order() != BondOrder::Double {
                return false;
            }
            break;
        }
    }

    true
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_isCarbonyl (BoundsMatrixBuilder.cpp:951-963)
// RDKit✔️✔️: bool _isCarbonyl(const ROMol &mol, const Atom *at) {
// RDKit✔️✔️:   if (at->getAtomicNum() == 6 && at->getDegree() > 2) {
// RDKit✔️✔️:     for (const auto nbr : mol.atomNeighbors(at)) {
// RDKit✔️✔️:       unsigned int atNum = nbr->getAtomicNum();
// RDKit✔️✔️:       if ((atNum == 8 || atNum == 7) &&
// RDKit✔️✔️:           mol.getBondBetweenAtoms(at->getIdx(), nbr->getIdx())->getBondType() ==
// RDKit✔️✔️:               Bond::DOUBLE) {
// RDKit✔️✔️:         return true;
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return false;
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::_isCarbonyl
fn is_carbonyl(mol: &Molecule, atom_idx: usize) -> bool {
    let atom = &mol.atoms()[atom_idx];
    if atom.atomic_number() != 6 || neighbors_for_atom(mol, atom_idx).len() <= 2 {
        return false;
    }
    neighbors_for_atom(mol, atom_idx)
        .into_iter()
        .any(|nbr_idx| {
            let at_num = mol.atoms()[nbr_idx].atomic_number();
            (at_num == 8 || at_num == 7)
                && mol.bonds()[bond_between_idx_simple(mol, atom_idx, nbr_idx).expect("bond")]
                    .order()
                    == BondOrder::Double
        })
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_checkAmideEster15 (BoundsMatrixBuilder.cpp:965-978)
// RDKit✔️❌: bool _checkAmideEster15(const ROMol &mol, const Bond *bnd1, const Bond *bnd3,
// RDKit✔️❌:                         const Atom *, const Atom *atm2, const Atom *atm3,
// RDKit✔️❌:                         const Atom *) {
// RDKit✔️❌:   unsigned int a2Num = atm2->getAtomicNum();
// RDKit✔️❌:   if ((a2Num == 8) || ((a2Num == 7) && (atm2->getTotalNumHs(true) == 1))) {
// RDKit✔️❌:     if ((bnd1->getBondType() == Bond::SINGLE)) {
// RDKit✔️❌:       if ((atm3->getAtomicNum() == 6) &&
// RDKit✔️❌:           (bnd3->getBondType() == Bond::SINGLE) && _isCarbonyl(mol, atm3)) {
// RDKit✔️❌:         return true;
// RDKit✔️❌:       }
// RDKit✔️❌:     }
// RDKit✔️❌:   }
// RDKit✔️❌:   return false;
// RDKit✔️❌: }
// END RDKIT CPP FUNCTION DGeomHelpers::_checkAmideEster15
fn check_amide_ester_15(
    mol: &Molecule,
    bnd1_idx: usize,
    bnd3_idx: usize,
    atm2: usize,
    atm3: usize,
) -> bool {
    let a2_num = mol.atoms()[atm2].atomic_number();
    // Behavior matches RDKit, but total-H lookup and nested carbonyl detection
    // still require repeated Rust-side scans.
    let total_hs_atm2 = total_num_hydrogens_rdkit_like(mol, atm2).unwrap_or(0);
    ((a2_num == 8) || (a2_num == 7 && total_hs_atm2 == 1))
        && mol.bonds()[bnd1_idx].order() == BondOrder::Single
        && mol.atoms()[atm3].atomic_number() == 6
        && mol.bonds()[bnd3_idx].order() == BondOrder::Single
        && is_carbonyl(mol, atm3)
}

// ──────────────────────────────────────────────
// 1-5 distance bounds (optional)
// ──────────────────────────────────────────────

/// RDKit✔️✔️: _compute15DistsCisCis — 1-5 distance, cis-cis configuration
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

/// RDKit✔️✔️: _compute15DistsCisTrans — 1-5 distance, cis-trans configuration
///
/// Configuration:
///  1     4-5
///   \   /
///    2-3
// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_compute15DistsCisTrans (BoundsMatrixBuilder.cpp:2004-2019)
// RDKit✔️✔️: double _compute15DistsCisTrans(double d1, double d2, double d3, double d4,
// RDKit✔️✔️:                                double ang12, double ang23, double ang34) {
// RDKit✔️✔️:   double dx14 = d2 - d3 * cos(ang23) - d1 * cos(ang12);
// RDKit✔️✔️:   double dy14 = d3 * sin(ang23) - d1 * sin(ang12);
// RDKit✔️✔️:   double d14 = sqrt(dx14 * dx14 + dy14 * dy14);
// RDKit✔️✔️:   double cval = (d3 - d2 * cos(ang23) + d1 * cos(ang12 + ang23)) / d14;
// RDKit✔️✔️:   if (cval > 1.0) { cval = 1.0; } else if (cval < -1.0) { cval = -1.0; }
// RDKit✔️✔️:   double ang143 = acos(cval);
// RDKit✔️✔️:   double ang145 = ang34 + ang143;
// RDKit✔️✔️:   return RDGeom::compute13Dist(d14, d4, ang145);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::_compute15DistsCisTrans
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

/// RDKit✔️✔️: _compute15DistsTransTrans — 1-5 distance, trans-trans configuration
///
/// Configuration:
///  1
///   \
///    2-3
///       \
///        4-5
// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_compute15DistsTransTrans (BoundsMatrixBuilder.cpp:2045-2060)
// RDKit✔️✔️: double _compute15DistsTransTrans(double d1, double d2, double d3, double d4,
// RDKit✔️✔️:                                  double ang12, double ang23, double ang34) {
// RDKit✔️✔️:   double dx14 = d2 - d3 * cos(ang23) - d1 * cos(ang12);
// RDKit✔️✔️:   double dy14 = d3 * sin(ang23) + d1 * sin(ang12);
// RDKit✔️✔️:   double d14 = sqrt(dx14 * dx14 + dy14 * dy14);
// RDKit✔️✔️:   double cval = (d3 - d2 * cos(ang23) + d1 * cos(ang12 - ang23)) / d14;
// RDKit✔️✔️:   if (cval > 1.0) { cval = 1.0; } else if (cval < -1.0) { cval = -1.0; }
// RDKit✔️✔️:   double ang143 = acos(cval);
// RDKit✔️✔️:   double ang145 = ang34 + ang143;
// RDKit✔️✔️:   return RDGeom::compute13Dist(d14, d4, ang145);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::_compute15DistsTransTrans
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

/// RDKit✔️✔️: _compute15DistsTransCis — 1-5 distance, trans-cis configuration
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
// RDKit✔️✔️: double _compute15DistsTransCis(double d1, double d2, double d3, double d4,
// RDKit✔️✔️:                                double ang12, double ang23, double ang34) {
// RDKit✔️✔️:   double dx14 = d2 - d3 * cos(ang23) - d1 * cos(ang12);
// RDKit✔️✔️:   double dy14 = d3 * sin(ang23) + d1 * sin(ang12);
// RDKit✔️✔️:   double d14 = sqrt(dx14 * dx14 + dy14 * dy14);
// RDKit✔️✔️:   double cval = (d3 - d2 * cos(ang23) + d1 * cos(ang12 - ang23)) / d14;
// RDKit✔️✔️:   if (cval > 1.0) { cval = 1.0; } else if (cval < -1.0) { cval = -1.0; }
// RDKit✔️✔️:   double ang143 = acos(cval);
// RDKit✔️✔️:   double ang145 = ang34 - ang143;
// RDKit✔️✔️:   return RDGeom::compute13Dist(d14, d4, ang145);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::_compute15DistsTransCis
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
// RDKit✔️✔️: void set15Bounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat,
// RDKit✔️✔️:                  ComputedData &accumData, double *distMatrix) {
// RDKit✔️✔️:   PATH14_VECT_CI pti;
// RDKit✔️✔️:   unsigned int bid1, bid2, bid3, type;
// RDKit✔️✔️:   for (pti = accumData.paths14.begin(); pti != accumData.paths14.end(); pti++) {
// RDKit✔️✔️:     bid1 = pti->bid1;
// RDKit✔️✔️:     bid2 = pti->bid2;
// RDKit✔️✔️:     bid3 = pti->bid3;
// RDKit✔️✔️:     type = pti->type;
// RDKit✔️✔️:     // 15 distances going one way with with 14 paths
// RDKit✔️✔️:     _set15BoundsHelper(mol, bid1, bid2, bid3, type, accumData, mmat,
// RDKit✔️✔️:                        distMatrix);
// RDKit✔️✔️:     // going the other way - reverse the 14 path
// RDKit✔️✔️:     _set15BoundsHelper(mol, bid3, bid2, bid1, type, accumData, mmat,
// RDKit✔️✔️:                        distMatrix);
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION
fn set_15_bounds(
    mol: &Molecule,
    mmat: &mut BoundsMatrix,
    accum_data: &mut ComputedData,
    dmat: &[f64],
) -> Result<(), DgBoundsError> {
    let nb = mol.bonds().len();
    let na = mol.atoms().len();
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

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::_set15BoundsHelper (BoundsMatrixBuilder.cpp:2106-2253)
// RDKit✔️✔️: void _set15BoundsHelper(const ROMol &mol, unsigned int bid1, unsigned int bid2,
// RDKit✔️✔️:                         unsigned int bid3, unsigned int type,
// RDKit✔️✔️:                         ComputedData &accumData, DistGeom::BoundsMatPtr mmat,
// RDKit✔️✔️:                         double *dmat) {
// RDKit✔️✔️:   unsigned int i, aid1, aid2, aid3, aid4, aid5;
// RDKit✔️✔️:   double d1, d2, d3, d4, ang12, ang23, ang34, du, dl, vw1, vw5;
// RDKit✔️✔️:   unsigned int nb = mol.getNumBonds();
// RDKit✔️✔️:   unsigned int na = mol.getNumAtoms();
// RDKit✔️✔️:   aid2 = accumData.bondAdj->getVal(bid1, bid2);
// RDKit✔️✔️:   aid1 = mol.getBondWithIdx(bid1)->getOtherAtomIdx(aid2);
// RDKit✔️✔️:   aid3 = accumData.bondAdj->getVal(bid2, bid3);
// RDKit✔️✔️:   aid4 = mol.getBondWithIdx(bid3)->getOtherAtomIdx(aid3);
// RDKit✔️✔️:   d1 = accumData.bondLengths[bid1];
// RDKit✔️✔️:   d2 = accumData.bondLengths[bid2];
// RDKit✔️✔️:   d3 = accumData.bondLengths[bid3];
// RDKit✔️✔️:   ang12 = accumData.bondAngles->getVal(bid1, bid2);
// RDKit✔️✔️:   ang23 = accumData.bondAngles->getVal(bid2, bid3);
// RDKit✔️✔️:   for (i = 0; i < nb; i++) {
// RDKit✔️✔️:     du = -1.0;
// RDKit✔️✔️:     dl = 0.0;
// RDKit✔️✔️:     if (accumData.bondAdj->getVal(bid3, i) == static_cast<int>(aid4)) {
// RDKit✔️✔️:       aid5 = mol.getBondWithIdx(i)->getOtherAtomIdx(aid4);
// RDKit✔️✔️:       // make sure we did not com back to the first atom in the path -
// RDKit✔️✔️:       // possible
// RDKit✔️✔️:       // with 4 membered rings
// RDKit✔️✔️:       // this is a fix for Issue 244
// RDKit✔️✔️:       const unsigned int pid = std::min(aid1, aid5) * na + std::max(aid1, aid5);
// RDKit✔️✔️:       if (accumData.visitedBound(pid, DistType::DIST14)) {
// RDKit✔️✔️:         return;
// RDKit✔️✔️:       }
// RDKit✔️✔️:       if (dmat[std::max(aid1, aid5) * mmat->numRows() + std::min(aid1, aid5)] <
// RDKit✔️✔️:           3.9) {
// RDKit✔️✔️:         continue;
// RDKit✔️✔️:       }
// RDKit✔️✔️:       if (aid1 != aid5) {  // FIX: do we need this
// RDKit✔️✔️:         unsigned int pid1 = aid1 * na + aid5;
// RDKit✔️✔️:         unsigned int pid2 = aid5 * na + aid1;
// RDKit✔️✔️:         if ((mmat->getLowerBound(aid1, aid5) < DIST12_DELTA) ||
// RDKit✔️✔️:             (accumData.set15Atoms[pid1]) || (accumData.set15Atoms[pid2])) {
// RDKit✔️✔️:           d4 = accumData.bondLengths[i];
// RDKit✔️✔️:           ang34 = accumData.bondAngles->getVal(bid3, i);
// RDKit✔️✔️:           unsigned long pathId =
// RDKit✔️✔️:               static_cast<unsigned long>(bid2) * nb * nb + (bid3)*nb + i;
// RDKit✔️✔️:           if (type == 0) {
// RDKit✔️✔️:             if (accumData.cisPaths.find(pathId) != accumData.cisPaths.end()) {
// RDKit✔️✔️:               dl = _compute15DistsCisCis(d1, d2, d3, d4, ang12, ang23, ang34);
// RDKit✔️✔️:               du = dl + DIST15_TOL;
// RDKit✔️✔️:               dl -= DIST15_TOL;
// RDKit✔️✔️:             } else if (accumData.transPaths.find(pathId) !=
// RDKit✔️✔️:                        accumData.transPaths.end()) {
// RDKit✔️✔️:               dl = _compute15DistsCisTrans(d1, d2, d3, d4, ang12, ang23, ang34);
// RDKit✔️✔️:               du = dl + DIST15_TOL;
// RDKit✔️✔️:               dl -= DIST15_TOL;
// RDKit✔️✔️:             } else {
// RDKit✔️✔️:               dl = _compute15DistsCisCis(d1, d2, d3, d4, ang12, ang23, ang34) -
// RDKit✔️✔️:                    DIST15_TOL;
// RDKit✔️✔️:               du =
// RDKit✔️✔️:                   _compute15DistsCisTrans(d1, d2, d3, d4, ang12, ang23, ang34) +
// RDKit✔️✔️:                   DIST15_TOL;
// RDKit✔️✔️:             }
// RDKit✔️✔️:           } else if (type == 1) {
// RDKit✔️✔️:             if (accumData.cisPaths.find(pathId) != accumData.cisPaths.end()) {
// RDKit✔️✔️:               dl = _compute15DistsTransCis(d1, d2, d3, d4, ang12, ang23, ang34);
// RDKit✔️✔️:               du = dl + DIST15_TOL;
// RDKit✔️✔️:               dl -= DIST15_TOL;
// RDKit✔️✔️:             } else if (accumData.transPaths.find(pathId) !=
// RDKit✔️✔️:                        accumData.transPaths.end()) {
// RDKit✔️✔️:               dl = _compute15DistsTransTrans(d1, d2, d3, d4, ang12, ang23,
// RDKit✔️✔️:                                              ang34);
// RDKit✔️✔️:               du = dl + DIST15_TOL;
// RDKit✔️✔️:               dl -= DIST15_TOL;
// RDKit✔️✔️:             } else {
// RDKit✔️✔️:               dl =
// RDKit✔️✔️:                   _compute15DistsTransCis(d1, d2, d3, d4, ang12, ang23, ang34) -
// RDKit✔️✔️:                   DIST15_TOL;
// RDKit✔️✔️:               du = _compute15DistsTransTrans(d1, d2, d3, d4, ang12, ang23,
// RDKit✔️✔️:                                              ang34) +
// RDKit✔️✔️:                    DIST15_TOL;
// RDKit✔️✔️:             }
// RDKit✔️✔️:           } else {
// RDKit✔️✔️:             if (accumData.cisPaths.find(pathId) != accumData.cisPaths.end()) {
// RDKit✔️✔️:               dl = _compute15DistsCisCis(d4, d3, d2, d1, ang34, ang23, ang12) -
// RDKit✔️✔️:                    DIST15_TOL;
// RDKit✔️✔️:               du =
// RDKit✔️✔️:                   _compute15DistsCisTrans(d4, d3, d2, d1, ang34, ang23, ang12) +
// RDKit✔️✔️:                   DIST15_TOL;
// RDKit✔️✔️:             } else if (accumData.transPaths.find(pathId) !=
// RDKit✔️✔️:                        accumData.transPaths.end()) {
// RDKit✔️✔️:               dl =
// RDKit✔️✔️:                   _compute15DistsTransCis(d4, d3, d2, d1, ang34, ang23, ang12) -
// RDKit✔️✔️:                   DIST15_TOL;
// RDKit✔️✔️:               du = _compute15DistsTransTrans(d4, d3, d2, d1, ang34, ang23,
// RDKit✔️✔️:                                              ang12) +
// RDKit✔️✔️:                    DIST15_TOL;
// RDKit✔️✔️:             } else {
// RDKit✔️✔️:               vw1 = PeriodicTable::getTable()->getRvdw(
// RDKit✔️✔️:                   mol.getAtomWithIdx(aid1)->getAtomicNum());
// RDKit✔️✔️:               vw5 = PeriodicTable::getTable()->getRvdw(
// RDKit✔️✔️:                   mol.getAtomWithIdx(aid5)->getAtomicNum());
// RDKit✔️✔️:               dl = VDW_SCALE_15 * (vw1 + vw5);
// RDKit✔️✔️:             }
// RDKit✔️✔️:           }
// RDKit✔️✔️:           if (du < 0.0) {
// RDKit✔️✔️:             du = MAX_UPPER;
// RDKit✔️✔️:           }
// RDKit✔️✔️:           _checkAndSetBounds(aid1, aid5, dl, du, mmat);
// RDKit✔️✔️:           accumData.set15Atoms[aid1 * na + aid5] = 1;
// RDKit✔️✔️:           accumData.set15Atoms[aid5 * na + aid1] = 1;
// RDKit✔️✔️:         }
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::_set15BoundsHelper
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
) -> Result<(), DgBoundsError> {
    // Behavior matches RDKit's 1-5 branch matrix, but cis/trans path flags are
    // still searched in Vec-backed stores with linear duplicate checks.
    // RDKit❗✔️: Get shared atoms via bond adjacency matrix
    let aid2 = bond_pair_shared_atom(mol, accum_data, bid1, bid2)?;
    let aid1 = if mol.bonds()[bid1].begin().index() == aid2 {
        mol.bonds()[bid1].end().index()
    } else {
        mol.bonds()[bid1].begin().index()
    };
    let aid3 = bond_pair_shared_atom(mol, accum_data, bid2, bid3)?;
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
        if accum_data.visited_bound(pid, DistType::Dist14) {
            return Ok(());
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

        // RDKit❗✔️: _checkAndSetBounds(aid1, aid5, dl, du, mmat);
        mmat.check_and_set_bounds(aid1, aid5, dl, du)?;

        // RDKit❗✔️: accumData.set15Atoms[pid1] = 1; accumData.set15Atoms[pid2] = 1;
        accum_data.set15_atoms[pid1] = true;
        accum_data.set15_atoms[pid2] = true;
    }
    Ok(())
}

/// RDKit❗✔️: set14Bounds — torsion-based 1-4 bounds
/// Also populates paths14, bond_adj, bond_angles for 1-5 bounds.
fn set_14_bounds(
    mol: &Molecule,
    mmat: &mut BoundsMatrix,
    accum_data: &mut ComputedData,
    dmat: &[f64],
    use_macrocycle_14config: bool,
    force_trans_amides: bool,
    rinfo: &crate::RingInfo,
) -> Result<(), DgBoundsError> {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::set14Bounds (BoundsMatrixBuilder.cpp:1664-1784)
    // RDKit❗✔️: void set14Bounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat,
    // RDKit❗✔️:                  ComputedData &accumData, double *distMatrix,
    // RDKit❗✔️:                  bool useMacrocycle14config, bool forceTransAmides) {
    // RDKit❗✔️:   unsigned int npt = mmat->numRows();
    // RDKit❗✔️:   CHECK_INVARIANT(npt == mol.getNumAtoms(), "Wrong size metric matrix");
    // RDKit❗✔️:   const size_t MAX_NUM_BONDS = static_cast<size_t>(
    // RDKit❗✔️:       std::pow(std::numeric_limits<std::uint64_t>::max(), 1. / 3));
    // RDKit❗✔️:   if (mol.getNumBonds() >= MAX_NUM_BONDS) {
    // RDKit❗✔️:     throw ValueErrorException(
    // RDKit❗✔️:         "Too many bonds in the molecule, cannot compute 1-4 bounds");
    // RDKit❗✔️:   }
    // RDKit❗✔️:   const auto rinfo = mol.getRingInfo();
    // RDKit❗✔️:   const auto &bondRings = rinfo->bondRings();
    // RDKit❗✔️:   std::unordered_set<unsigned int> bidIsMacrocycle;
    // RDKit❗✔️:   std::unordered_set<std::uint64_t> ringBondPairs;
    // RDKit❗✔️:   std::unordered_set<std::uint64_t> donePaths;
    // RDKit❗✔️:   std::uint64_t nb = mol.getNumBonds();
    // RDKit❗✔️:   // first we will deal with 1-4 atoms that belong to the same ring
    // RDKit❗✔️:   for (const auto &bring : bondRings) {
    // RDKit❗✔️:     ... same-ring dispatch to _setMacrocycleAllInSameRing14Bounds,
    // RDKit❗✔️:         _setInRing14Bounds, or _record14Path ...
    // RDKit❗✔️:   }
    // RDKit❗✔️:   for (const auto bond : mol.bonds()) {
    // RDKit❗✔️:     ... donePaths gate, then dispatch to _setMacrocycleTwoInSameRing14Bounds,
    // RDKit❗✔️:         _setTwoInSameRing14Bounds, _setTwoInDiffRing14Bounds,
    // RDKit❗✔️:         _setShareRingBond14Bounds, or _setChain14Bounds ...
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::set14Bounds
    let npt = mmat.num_rows();
    if npt != mol.num_atoms() {
        return Err(dg_generation_error("Wrong size metric matrix"));
    }
    let max_num_bonds = (u64::MAX as f64).cbrt() as usize;
    if mol.num_bonds() >= max_num_bonds {
        return Err(dg_generation_error(
            "Too many bonds in the molecule, cannot compute 1-4 bounds",
        ));
    }
    let bond_rings = rinfo.bond_rings();

    let mut bid_is_macrocycle: HashSet<usize> = HashSet::new();
    let mut ring_bond_pairs: HashSet<u64> = HashSet::new();
    let mut done_paths: HashSet<u64> = HashSet::new();
    let nb = mol.num_bonds() as u64;

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
                        mol, bid1, bid2, bid3, accum_data, mmat,
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

    for bond in mol.bonds() {
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

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::collectBondsAndAngles (BoundsMatrixBuilder.cpp:1847-1898)
// RDKit✔️✔️: void collectBondsAndAngles(const ROMol &mol,
// RDKit✔️✔️:                            std::vector<std::pair<int, int>> &bonds,
// RDKit✔️✔️:                            std::vector<std::vector<int>> &angles) {
// RDKit✔️✔️:   bonds.resize(0);
// RDKit✔️✔️:   angles.resize(0);
// RDKit✔️✔️:   bonds.reserve(mol.getNumBonds());
// RDKit✔️✔️:   for (const auto bondi : mol.bonds()) {
// RDKit✔️✔️:     bonds.emplace_back(bondi->getBeginAtomIdx(), bondi->getEndAtomIdx());
// RDKit✔️✔️:     for (unsigned int j = bondi->getIdx() + 1; j < mol.getNumBonds(); ++j) {
// RDKit✔️✔️:       const Bond *bondj = mol.getBondWithIdx(j);
// RDKit✔️✔️:       int aid11 = bondi->getBeginAtomIdx();
// RDKit✔️✔️:       int aid12 = bondi->getEndAtomIdx();
// RDKit✔️✔️:       int aid21 = bondj->getBeginAtomIdx();
// RDKit✔️✔️:       int aid22 = bondj->getEndAtomIdx();
// RDKit✔️✔️:       if (aid11 != aid21 && aid11 != aid22 && aid12 != aid21 &&
// RDKit✔️✔️:           aid12 != aid22) {
// RDKit✔️✔️:         continue;
// RDKit✔️✔️:       }
// RDKit✔️✔️:       std::vector<int> tmp(4,
// RDKit✔️✔️:                            0);  // elements: aid1, aid2, flag for triple bonds
// RDKit✔️✔️:       if (aid12 == aid21) {
// RDKit✔️✔️:         tmp[0] = aid11;
// RDKit✔️✔️:         tmp[1] = aid12;
// RDKit✔️✔️:         tmp[2] = aid22;
// RDKit✔️✔️:       } else if (aid12 == aid22) {
// RDKit✔️✔️:         tmp[0] = aid11;
// RDKit✔️✔️:         tmp[1] = aid12;
// RDKit✔️✔️:         tmp[2] = aid21;
// RDKit✔️✔️:       } else if (aid11 == aid21) {
// RDKit✔️✔️:         tmp[0] = aid12;
// RDKit✔️✔️:         tmp[1] = aid11;
// RDKit✔️✔️:         tmp[2] = aid22;
// RDKit✔️✔️:       } else if (aid11 == aid22) {
// RDKit✔️✔️:         tmp[0] = aid12;
// RDKit✔️✔️:         tmp[1] = aid11;
// RDKit✔️✔️:         tmp[2] = aid21;
// RDKit✔️✔️:       }
// RDKit✔️✔️:       if (bondi->getBondType() == Bond::TRIPLE ||
// RDKit✔️✔️:           bondj->getBondType() == Bond::TRIPLE) {
// RDKit✔️✔️:         // triple bond
// RDKit✔️✔️:         tmp[3] = 1;
// RDKit✔️✔️:       } else if (bondi->getBondType() == Bond::DOUBLE &&
// RDKit✔️✔️:                  bondj->getBondType() == Bond::DOUBLE &&
// RDKit✔️✔️:                  mol.getAtomWithIdx(tmp[1])->getDegree() == 2) {
// RDKit✔️✔️:         // consecutive double bonds
// RDKit✔️✔️:         tmp[3] = 1;
// RDKit✔️✔️:       }
// RDKit✔️✔️:       angles.push_back(tmp);
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::collectBondsAndAngles
fn collect_bonds_and_angles(
    mol: &Molecule,
    bonds: &mut Vec<(i32, i32)>,
    angles: &mut Vec<Vec<i32>>,
) {
    bonds.clear();
    angles.clear();
    bonds.reserve(mol.num_bonds());
    for bondi in mol.bonds() {
        bonds.push((bondi.begin().index() as i32, bondi.end().index() as i32));

        for j in (bondi.id().index() + 1)..mol.num_bonds() {
            let bondj = &mol.bonds()[j];
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
                && neighbors_for_atom(mol, tmp[1] as usize).len() == 2
            {
                tmp[3] = 1;
            }

            angles.push(tmp);
        }
    }
}

// ──────────────────────────────────────────────
// VDW lower bounds
// ──────────────────────────────────────────────

/// RDKit❗✔️: setLowerBoundVDW — VDW lower bounds for non-bonded pairs
// BEGIN RDKIT CPP FUNCTION DGeomHelpers::setLowerBoundVDW (BoundsMatrixBuilder.cpp:315-361)
// RDKit✔️✔️: void setLowerBoundVDW(const ROMol &mol, DistGeom::BoundsMatPtr mmat, bool,
// RDKit✔️✔️:                       double *dmat) {
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
fn set_lower_bound_vdw(
    mol: &Molecule,
    mmat: &mut BoundsMatrix,
    _scale_vdw: bool,
    dmat: &[f64],
) -> Result<(), DgBoundsError> {
    let n = mol.atoms().len();
    let npt = mmat.num_rows();
    if npt != n {
        return Err(dg_generation_error("Wrong size metric matrix"));
    }

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
            let td = dmat[i * npt + j];

            // H-bond special case: set lower bound to 1.8A for donor-H + acceptor
            if (h_in_hbond_donors[i] && hbond_acceptors[j])
                || (hbond_acceptors[i] && h_in_hbond_donors[j])
            {
                mmat.set_lower(i, j, H_BOND_LENGTH)?;
            } else if td == 4.0 {
                // 1-5: scaled VDW
                mmat.set_lower(i, j, VDW_SCALE_15 * (vw1 + vw2))?;
            } else if td == 5.0 {
                // 1-6: slightly less scaled
                mmat.set_lower(
                    i,
                    j,
                    (VDW_SCALE_15 + 0.5 * (1.0 - VDW_SCALE_15)) * (vw1 + vw2),
                )?;
            } else {
                // Full VDW sum
                mmat.set_lower(i, j, vw1 + vw2)?;
            }
        }
    }
    Ok(())
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::setTopolBounds (BoundsMatrixBuilder.cpp:1808-1845)
// RDKit✔️✔️: void setTopolBounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat,
// RDKit✔️✔️:                     bool set15bounds, bool scaleVDW, bool useMacrocycle14config,
// RDKit✔️✔️:                     bool forceTransAmides, bool set14bounds, bool set13bounds) {
// RDKit✔️✔️:   PRECONDITION(mmat.get(), "bad pointer");
// RDKit✔️✔️:   unsigned int nb = mol.getNumBonds();
// RDKit✔️✔️:   unsigned int na = mol.getNumAtoms();
// RDKit✔️✔️:   if (!na) {
// RDKit✔️✔️:     throw ValueErrorException("molecule has no atoms");
// RDKit✔️✔️:   }
// RDKit✔️✔️:   const size_t MAX_NUM_BONDS = static_cast<size_t>(
// RDKit✔️✔️:       std::pow(std::numeric_limits<std::uint64_t>::max(), 1. / 3));
// RDKit✔️✔️:   if (mol.getNumBonds() >= MAX_NUM_BONDS) {
// RDKit✔️✔️:     throw ValueErrorException(
// RDKit✔️✔️:         "Too many bonds in the molecule, cannot compute 1-4 bounds");
// RDKit✔️✔️:   }
// RDKit✔️✔️:   ComputedData accumData(na, nb);
// RDKit✔️✔️:   double *distMatrix = nullptr;
// RDKit✔️✔️:   distMatrix = MolOps::getDistanceMat(mol);
// RDKit✔️✔️:   set12Bounds(mol, mmat, accumData);
// RDKit✔️✔️:   if (set13bounds) {
// RDKit✔️✔️:     set13Bounds(mol, mmat, accumData);
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (set14bounds) {
// RDKit✔️✔️:     set14Bounds(mol, mmat, accumData, distMatrix, useMacrocycle14config,
// RDKit✔️✔️:                 forceTransAmides);
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (set15bounds) {
// RDKit✔️✔️:     set15Bounds(mol, mmat, accumData, distMatrix);
// RDKit✔️✔️:   }
// RDKit✔️✔️:   setLowerBoundVDW(mol, mmat, scaleVDW, distMatrix);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::setTopolBounds
#[allow(clippy::too_many_arguments)]
fn set_topol_bounds(
    mol: &Molecule,
    mmat: &mut BoundsMatrix,
    set15bounds: bool,
    scale_vdw: bool,
    use_macrocycle_14config: bool,
    force_trans_amides: bool,
    set14bounds: bool,
    set13bounds: bool,
) -> Result<(), DgBoundsError> {
    let nb = mol.num_bonds();
    let na = mol.num_atoms();
    if na == 0 {
        return Err(DgBoundsError::GenerationFailed(
            "molecule has no atoms".to_string(),
        ));
    }
    let max_num_bonds = (u64::MAX as f64).cbrt() as usize;
    if nb >= max_num_bonds {
        return Err(DgBoundsError::GenerationFailed(
            "Too many bonds in the molecule, cannot compute 1-4 bounds".to_string(),
        ));
    }

    let mut accum_data = ComputedData::new(na, nb);
    let dist_matrix = flatten_topological_distances_matrix(mol);
    let ring_info = if set13bounds || set14bounds {
        Some(ring_info_for_distgeom(mol)?)
    } else {
        None
    };

    set_12_bounds(mol, mmat, &mut accum_data)?;
    trace_bounds_stage("after_set12", mmat);
    if set13bounds {
        set_13_bounds(
            mol,
            mmat,
            &mut accum_data,
            ring_info.as_deref().expect("ring info loaded"),
        )?;
        trace_bounds_stage("after_set13", mmat);
    }
    if set14bounds {
        set_14_bounds(
            mol,
            mmat,
            &mut accum_data,
            &dist_matrix,
            use_macrocycle_14config,
            force_trans_amides,
            ring_info.as_deref().expect("ring info loaded"),
        )?;
        trace_bounds_stage("after_set14", mmat);
    }
    if set15bounds {
        set_15_bounds(mol, mmat, &mut accum_data, &dist_matrix)?;
        trace_bounds_stage("after_set15", mmat);
    }
    set_lower_bound_vdw(mol, mmat, scale_vdw, &dist_matrix)?;
    trace_bounds_stage("after_vdw", mmat);
    Ok(())
}

// BEGIN RDKIT CPP FUNCTION DGeomHelpers::setTopolBounds (BoundsMatrixBuilder.cpp:1903-1942)
// RDKit✔️✔️: void setTopolBounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat,
// RDKit✔️✔️:                     std::vector<std::pair<int, int>> &bonds,
// RDKit✔️✔️:                     std::vector<std::vector<int>> &angles, bool set15bounds,
// RDKit✔️✔️:                     bool scaleVDW, bool useMacrocycle14config,
// RDKit✔️✔️:                     bool forceTransAmides, bool set14bounds, bool set13bounds) {
// RDKit✔️✔️:   PRECONDITION(mmat.get(), "bad pointer");
// RDKit✔️✔️:   bonds.clear();
// RDKit✔️✔️:   angles.clear();
// RDKit✔️✔️:   unsigned int nb = mol.getNumBonds();
// RDKit✔️✔️:   unsigned int na = mol.getNumAtoms();
// RDKit✔️✔️:   if (!na) {
// RDKit✔️✔️:     throw ValueErrorException("molecule has no atoms");
// RDKit✔️✔️:   }
// RDKit✔️✔️:   ComputedData accumData(na, nb);
// RDKit✔️✔️:   double *distMatrix = nullptr;
// RDKit✔️✔️:   distMatrix = MolOps::getDistanceMat(mol);
// RDKit✔️✔️:   set12Bounds(mol, mmat, accumData);
// RDKit✔️✔️:   if (set13bounds) {
// RDKit✔️✔️:     set13Bounds(mol, mmat, accumData);
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (set14bounds) {
// RDKit✔️✔️:     set14Bounds(mol, mmat, accumData, distMatrix, useMacrocycle14config,
// RDKit✔️✔️:                 forceTransAmides);
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (set15bounds) {
// RDKit✔️✔️:     set15Bounds(mol, mmat, accumData, distMatrix);
// RDKit✔️✔️:   }
// RDKit✔️✔️:   setLowerBoundVDW(mol, mmat, scaleVDW, distMatrix);
// RDKit✔️✔️:   collectBondsAndAngles(mol, bonds, angles);
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION DGeomHelpers::setTopolBounds
#[allow(clippy::too_many_arguments)]
fn set_topol_bounds_with_outputs(
    mol: &Molecule,
    mmat: &mut BoundsMatrix,
    bonds: &mut Vec<(i32, i32)>,
    angles: &mut Vec<Vec<i32>>,
    set15bounds: bool,
    scale_vdw: bool,
    use_macrocycle_14config: bool,
    force_trans_amides: bool,
    set14bounds: bool,
    set13bounds: bool,
) -> Result<(), DgBoundsError> {
    bonds.clear();
    angles.clear();
    set_topol_bounds(
        mol,
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

pub(crate) fn dg_path14_priors(mol: &Molecule) -> Result<Vec<DgPath14Prior>, DgBoundsError> {
    let nb = mol.num_bonds();
    let na = mol.num_atoms();
    if na == 0 {
        return Err(DgBoundsError::GenerationFailed(
            "molecule has no atoms".to_string(),
        ));
    }
    let max_num_bonds = (u64::MAX as f64).cbrt() as usize;
    if nb >= max_num_bonds {
        return Err(DgBoundsError::GenerationFailed(
            "Too many bonds in the molecule, cannot compute 1-4 bounds".to_string(),
        ));
    }

    let mut mmat = BoundsMatrix::new(na);
    let mut accum_data = ComputedData::new(na, nb);
    let dist_matrix = flatten_topological_distances_matrix(mol);
    let ring_info = ring_info_for_distgeom(mol)?;
    set_12_bounds(mol, &mut mmat, &mut accum_data)?;
    set_13_bounds(mol, &mut mmat, &mut accum_data, ring_info.as_ref())?;
    set_14_bounds(
        mol,
        &mut mmat,
        &mut accum_data,
        &dist_matrix,
        false,
        true,
        ring_info.as_ref(),
    )?;

    let mut priors = Vec::with_capacity(accum_data.paths14.len());
    for path in accum_data.paths14 {
        let Some((i, j)) = bond_atoms_for_path14_outer(mol, path.bid1, path.bid2) else {
            continue;
        };
        let Some((k, l)) = bond_atoms_for_path14_outer(mol, path.bid3, path.bid2) else {
            continue;
        };
        let kind = match path.kind {
            Path14Kind::Cis => DgPath14Kind::Cis,
            Path14Kind::Trans => DgPath14Kind::Trans,
            Path14Kind::Other => DgPath14Kind::Other,
        };
        priors.push(DgPath14Prior {
            atoms: (i, j, k, l),
            kind,
            lower_bound: mmat.get_lower(i, l),
            upper_bound: mmat.get_upper(i, l),
        });
    }
    Ok(priors)
}

fn bond_atoms_for_path14_outer(
    mol: &Molecule,
    outer_bond_idx: usize,
    middle_bond_idx: usize,
) -> Option<(usize, usize)> {
    let outer = mol.bonds().get(outer_bond_idx)?;
    let middle = mol.bonds().get(middle_bond_idx)?;
    let outer_begin = outer.begin().index();
    let outer_end = outer.end().index();
    let middle_begin = middle.begin().index();
    let middle_end = middle.end().index();
    if outer_begin == middle_begin || outer_begin == middle_end {
        Some((outer_end, outer_begin))
    } else if outer_end == middle_begin || outer_end == middle_end {
        Some((outer_begin, outer_end))
    } else {
        None
    }
}

// ──────────────────────────────────────────────
// Public entry point
// ──────────────────────────────────────────────

// BEGIN RDKIT CPP FUNCTION RDKit::getMolBoundsMatrix (rdDistGeom.cpp:223-242)
// RDKit✔️❌: PyObject *getMolBoundsMatrix(ROMol &mol, bool set15bounds = true,
// RDKit✔️❌:                              bool scaleVDW = false,
// RDKit✔️❌:                              bool doTriangleSmoothing = true,
// RDKit✔️❌:                              bool useMacrocycle14config = false) {
// RDKit✔️❌:   unsigned int nats = mol.getNumAtoms();
// RDKit✔️❌:   npy_intp dims[2];
// RDKit✔️❌:   dims[0] = nats;
// RDKit✔️❌:   dims[1] = nats;
// RDKit✔️❌:
// RDKit✔️❌:   DistGeom::BoundsMatPtr mat(new DistGeom::BoundsMatrix(nats));
// RDKit✔️❌:   DGeomHelpers::initBoundsMat(mat);
// RDKit✔️❌:   DGeomHelpers::setTopolBounds(mol, mat, set15bounds, scaleVDW,
// RDKit✔️❌:                                useMacrocycle14config);
// RDKit✔️❌:   if (doTriangleSmoothing) {
// RDKit✔️❌:     DistGeom::triangleSmoothBounds(mat);
// RDKit✔️❌:   }
// RDKit✔️❌:   auto *res = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
// RDKit✔️❌:   memcpy(static_cast<void *>(PyArray_DATA(res)),
// RDKit✔️❌:          static_cast<void *>(mat->getData()), nats * nats * sizeof(double));
// RDKit✔️❌:   return PyArray_Return(res);
// RDKit✔️❌: }
// END RDKIT CPP FUNCTION RDKit::getMolBoundsMatrix
//
// Local performance review keeps the second marker at `❌`: the Rust wrapper
// matches RDKit's control flow and default flag wiring, but it must allocate a
// nested `Vec<Vec<f64>>` and copy row-by-row from `BoundsMatrix` instead of
// handing a single contiguous allocation directly to NumPy with `memcpy()`.
pub fn dg_bounds_matrix_with_options(
    molecule: &Molecule,
    set15bounds: bool,
    scale_vdw: bool,
    do_triangle_smoothing: bool,
    use_macrocycle14config: bool,
) -> Result<Vec<Vec<f64>>, DgBoundsError> {
    let n = molecule.atoms().len();
    let mut mmat = BoundsMatrix::new(n);
    set_topol_bounds(
        molecule,
        &mut mmat,
        set15bounds,
        scale_vdw,
        use_macrocycle14config,
        true,
        true,
        true,
    )?;
    if do_triangle_smoothing {
        // RDKit✔️❌:   if (doTriangleSmoothing) {
        // RDKit✔️❌:     DistGeom::triangleSmoothBounds(mat);
        // RDKit✔️❌:   }
        //
        // The Python wrapper intentionally ignores the boolean return from
        // triangle smoothing and still exposes the mutated matrix.
        let _ = mmat.triangle_smooth(0.0);
    }
    Ok(mmat.to_vec_vec())
}

pub fn dg_bounds_matrix(molecule: &Molecule) -> Result<Vec<Vec<f64>>, DgBoundsError> {
    dg_bounds_matrix_with_options(molecule, true, false, true, false)
}

fn molecule_fragment_mapping(mol: &Molecule) -> Vec<usize> {
    let mut mapping = vec![usize::MAX; mol.num_atoms()];
    let mut frag_idx = 0;
    for atom_idx in 0..mol.num_atoms() {
        if mapping[atom_idx] != usize::MAX {
            continue;
        }
        let mut stack = vec![atom_idx];
        mapping[atom_idx] = frag_idx;
        while let Some(current) = stack.pop() {
            for nbr in neighbors_for_atom(mol, current) {
                if mapping[nbr] == usize::MAX {
                    mapping[nbr] = frag_idx;
                    stack.push(nbr);
                }
            }
        }
        frag_idx += 1;
    }
    mapping
}

fn molecule_fragments_for_embed(
    mol: &Molecule,
    embed_fragments_separately: bool,
) -> Result<(Vec<Molecule>, Vec<usize>), DgBoundsError> {
    if !embed_fragments_separately {
        return Ok((vec![mol.clone()], vec![0; mol.num_atoms()]));
    }
    let frag_mapping = molecule_fragment_mapping(mol);
    let fragment_count = frag_mapping
        .iter()
        .copied()
        .max()
        .map_or(0, |max_frag| max_frag + 1);
    if fragment_count <= 1 {
        return Ok((vec![mol.clone()], vec![0; mol.num_atoms()]));
    }
    let mut fragments = Vec::with_capacity(fragment_count);
    for frag_idx in 0..fragment_count {
        let atoms_to_remove: Vec<_> = mol
            .atoms()
            .iter()
            .filter_map(|atom| (frag_mapping[atom.id().index()] != frag_idx).then_some(atom.id()))
            .collect();
        let mut builder = mol.to_builder();
        builder.remove_atoms_for_construction(&atoms_to_remove);
        fragments.push(builder.build()?);
    }
    Ok((fragments, frag_mapping))
}

fn conformer_from_positions(id: usize, atom_count: usize) -> Conformer3D {
    Conformer3D::new(id, vec![[0.0, 0.0, 0.0]; atom_count], true)
}

fn embedder_add_conformers(
    mol: &Molecule,
    new_conformers: Vec<Conformer3D>,
    clear_confs: bool,
) -> Result<(Molecule, Vec<i32>), DgBoundsError> {
    let mut out = mol.clone();
    let coordinates = out.coordinate_block_mut();
    if clear_confs {
        coordinates.conformers_3d.clear();
    }
    let mut ids = Vec::with_capacity(new_conformers.len());
    for conformer in new_conformers {
        let conf_id = coordinates.conformers_3d.len();
        ids.push(conf_id as i32);
        coordinates.conformers_3d.push(Conformer3D::new(
            conf_id,
            conformer.coordinates().to_vec(),
            conformer.is_3d(),
        ));
    }
    Ok((out, ids))
}

fn molecule_from_read_parts_with_coordinates(
    read: MoleculeReadParts<'_>,
    coordinates: MoleculeCoordinateBlock,
) -> Molecule {
    Molecule::from_operation_blocks(
        read.topology().clone(),
        coordinates,
        read.properties().clone(),
        read.derived_cache().clone(),
        read.capabilities(),
    )
    .expect("an operation-owned coordinate block must preserve molecule invariants")
}

pub(crate) fn embed_multiple_confs_coordinate_update(
    read: MoleculeReadParts<'_>,
    coordinates: &mut MoleculeCoordinateBlock,
    num_confs: u32,
    params: &mut EmbedParameters,
) -> Result<Vec<i32>, DgBoundsError> {
    let owned_coordinates = std::mem::take(coordinates);
    let mut source = molecule_from_read_parts_with_coordinates(read, owned_coordinates);
    match embed_multiple_confs(&source, num_confs, params) {
        Ok((mut embedded, ids)) => {
            drop(source);
            *coordinates = embedded.take_coordinate_block_or_clone();
            Ok(ids)
        }
        Err(error) => {
            *coordinates = source.take_coordinate_block_or_clone();
            Err(error)
        }
    }
}

pub(crate) fn embed_molecule_coordinate_update(
    read: MoleculeReadParts<'_>,
    coordinates: &mut MoleculeCoordinateBlock,
    params: &mut EmbedParameters,
) -> Result<i32, DgBoundsError> {
    let conf_ids = embed_multiple_confs_coordinate_update(read, coordinates, 1, params)?;
    Ok(conf_ids.first().copied().unwrap_or(-1))
}

pub fn embed_multiple_confs(
    mol: &Molecule,
    num_confs: u32,
    params: &mut EmbedParameters,
) -> Result<(Molecule, Vec<i32>), DgBoundsError> {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::EmbedMultipleConfs(ROMol&, INT_VECT&, unsigned int, EmbedParameters&) (Embedder.cpp:1501-1695)
    // RDKit✔️✔️: void EmbedMultipleConfs(ROMol &mol, INT_VECT &res, unsigned int numConfs,
    // RDKit✔️✔️:                         EmbedParameters &params) {
    // RDKit✔️✔️:   TimePoint *end_time = nullptr;
    // RDKit✔️✔️:   TimePoint end_time_storage;
    // RDKit✔️✔️:   if (params.timeout > 0) {
    // RDKit✔️✔️:     end_time_storage = Clock::now() + std::chrono::seconds(params.timeout);
    // RDKit✔️✔️:     end_time = &end_time_storage;
    // RDKit✔️✔️:   }
    let end_time = (params.timeout > 0)
        .then(|| Instant::now() + Duration::from_secs(u64::from(params.timeout)));

    // RDKit✔️✔️:   if (params.trackFailures) {
    // RDKit✔️✔️:     params.failures.resize(EmbedFailureCauses::END_OF_ENUM);
    // RDKit✔️✔️:     std::fill(params.failures.begin(), params.failures.end(), 0);
    // RDKit✔️✔️:   }
    if params.track_failures {
        params
            .failures
            .resize(EmbedFailureCause::EndOfEnum as usize, 0);
        params.failures.fill(0);
    }
    // RDKit✔️✔️:   if (!mol.getNumAtoms()) {
    // RDKit✔️✔️:     throw ValueErrorException("molecule has no atoms");
    // RDKit✔️✔️:   }
    if mol.num_atoms() == 0 {
        return Err(DgBoundsError::EmptyMolecule);
    }
    // RDKit✔️✔️:   if (params.ETversion < 1 || params.ETversion > 2) {
    // RDKit✔️✔️:     throw ValueErrorException(
    // RDKit✔️✔️:         "Only version 1 and 2 of the experimental "
    // RDKit✔️✔️:         "torsion-angle preferences (ETversion) supported");
    // RDKit✔️✔️:   }
    if params.et_version < 1 || params.et_version > 2 {
        return Err(DgBoundsError::UnsupportedEtVersion);
    }

    // RDKit✔️✔️:   if (MolOps::needsHs(mol)) {
    // RDKit✔️✔️:     BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:         << "Molecule does not have explicit Hs. Consider calling AddHs()"
    // RDKit✔️✔️:         << std::endl;
    // RDKit✔️✔️:   }

    // RDKit✔️✔️:   // initialize the conformers we're going to be creating:
    // RDKit✔️✔️:   if (params.clearConfs) {
    // RDKit✔️✔️:     res.clear();
    // RDKit✔️✔️:     mol.clearConformers();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::vector<std::unique_ptr<Conformer>> confs;
    // RDKit✔️✔️:   confs.reserve(numConfs);
    // RDKit✔️✔️:   for (unsigned int i = 0; i < numConfs; ++i) {
    // RDKit✔️✔️:     confs.emplace_back(new Conformer(mol.getNumAtoms()));
    // RDKit✔️✔️:   }
    let mut confs: Vec<Conformer3D> = (0..num_confs as usize)
        .map(|i| conformer_from_positions(i, mol.num_atoms()))
        .collect();
    // RDKit✔️✔️:   boost::dynamic_bitset<> confsOk(numConfs);
    // RDKit✔️✔️:   confsOk.set();
    let mut confs_ok = vec![true; num_confs as usize];

    // RDKit✔️✔️:   INT_VECT fragMapping;
    // RDKit✔️✔️:   std::vector<ROMOL_SPTR> molFrags;
    // RDKit✔️✔️:   if (params.embedFragmentsSeparately) {
    // RDKit✔️✔️:     molFrags = MolOps::getMolFrags(mol, true, &fragMapping);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     molFrags.push_back(ROMOL_SPTR(new ROMol(mol)));
    // RDKit✔️✔️:     fragMapping.resize(mol.getNumAtoms());
    // RDKit✔️✔️:     std::fill(fragMapping.begin(), fragMapping.end(), 0);
    // RDKit✔️✔️:   }
    let (mol_frags, frag_mapping) =
        molecule_fragments_for_embed(mol, params.embed_fragments_separately)?;

    // RDKit✔️✔️:   const std::map<int, RDGeom::Point3D> *coordMap = params.coordMap;
    let coord_map_storage = params.coord_map.clone();
    let mut coord_map = coord_map_storage.as_ref();
    // RDKit✔️✔️:   if (molFrags.size() > 1 && coordMap) {
    // RDKit✔️✔️:     BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:         << "Constrained conformer generation (via the coordMap argument) "
    // RDKit✔️✔️:            "does not work with molecules that have multiple fragments."
    // RDKit✔️✔️:         << std::endl;
    // RDKit✔️✔️:     coordMap = nullptr;
    // RDKit✔️✔️:   }
    if mol_frags.len() > 1 && coord_map.is_some() {
        coord_map = None;
    }
    // RDKit✔️✔️:   boost::dynamic_bitset<> constrainedAtoms(mol.getNumAtoms());
    // RDKit✔️✔️:   if (coordMap) {
    // RDKit✔️✔️:     for (const auto &entry : *coordMap) {
    // RDKit✔️✔️:       constrainedAtoms.set(entry.first);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (molFrags.size() > 1 && params.boundsMat != nullptr) {
    // RDKit✔️✔️:     BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:         << "Conformer generation using a user-provided boundsMat "
    // RDKit✔️✔️:            "does not work with molecules that have multiple fragments. The "
    // RDKit✔️✔️:            "boundsMat will be ignored."
    // RDKit✔️✔️:         << std::endl;
    // RDKit✔️✔️:     coordMap = nullptr;  // FIXME not directly related to ETKDG, but here I
    // RDKit✔️✔️:                          // think it should be params.boundsMat = nullptr
    // RDKit✔️✔️:   }
    if mol_frags.len() > 1 && params.bounds_mat.is_some() {
        coord_map = None;
    }

    // RDKit✔️✔️:   // we will generate conformations for each fragment in the molecule
    // RDKit✔️✔️:   // separately, so loop over them:
    // RDKit✔️✔️:   for (unsigned int fragIdx = 0; fragIdx < molFrags.size(); ++fragIdx) {
    for (frag_idx, piece) in mol_frags.iter().enumerate() {
        // RDKit✔️✔️:     ROMOL_SPTR piece = molFrags[fragIdx];
        // RDKit✔️✔️:     unsigned int nAtoms = piece->getNumAtoms();
        let n_atoms = piece.num_atoms();
        // RDKit✔️✔️:     ForceFields::CrystalFF::CrystalFFDetails etkdgDetails;
        // RDKit✔️✔️:     etkdgDetails.constrainedAtoms = constrainedAtoms;
        // RDKit✔️✔️:     EmbeddingOps::initETKDG(piece.get(), params, etkdgDetails);
        let mut etkdg_details = CrystalFFDetails::default();
        embedder_init_etkdg(piece, params, &mut etkdg_details)?;

        // RDKit✔️✔️:     DistGeom::BoundsMatPtr mmat;
        let mut mmat;
        // RDKit✔️✔️:     if (params.boundsMat == nullptr || molFrags.size() > 1) {
        if params.bounds_mat.is_none() || mol_frags.len() > 1 {
            // RDKit✔️✔️:       mmat.reset(new DistGeom::BoundsMatrix(nAtoms));
            // RDKit✔️✔️:       initBoundsMat(mmat);
            mmat = BoundsMatrix::new(n_atoms);
            // RDKit✔️✔️:       if (!EmbeddingOps::setupInitialBoundsMatrix(piece.get(), mmat, coordMap,
            // RDKit✔️✔️:                                                   params, etkdgDetails)) {
            // RDKit✔️✔️:         // return if we couldn't setup the bounds matrix
            // RDKit✔️✔️:         // possible causes include a triangle smoothing failure
            // RDKit✔️✔️:         return;
            // RDKit✔️✔️:       }
            if !embedder_setup_initial_bounds_matrix(
                piece,
                &mut mmat,
                coord_map,
                params,
                &mut etkdg_details,
            )? {
                return embedder_add_conformers(mol, Vec::new(), params.clear_confs);
            }
        } else {
            // RDKit✔️✔️:     } else {
            // RDKit✔️✔️:       // just use what they gave us
            // RDKit✔️✔️:       // first make sure it's the right size though:
            // RDKit✔️✔️:       if (params.boundsMat->numRows() != nAtoms) {
            // RDKit✔️✔️:         throw ValueErrorException(
            // RDKit✔️✔️:             "size of boundsMat provided does not match the number of atoms in "
            // RDKit✔️✔️:             "the molecule.");
            // RDKit✔️✔️:       }
            let bounds = params.bounds_mat.as_ref().expect("checked above");
            if bounds.num_rows() != n_atoms {
                return Err(DgBoundsError::GenerationFailed(
                    "size of boundsMat provided does not match the number of atoms in the molecule."
                        .to_string(),
                ));
            }
            // RDKit✔️✔️:       collectBondsAndAngles((*piece.get()), etkdgDetails.bonds,
            // RDKit✔️✔️:                             etkdgDetails.angles);
            collect_bonds_and_angles(piece, &mut etkdg_details.bonds, &mut etkdg_details.angles);
            // RDKit✔️✔️:       mmat.reset(new DistGeom::BoundsMatrix(*params.boundsMat));
            mmat = (**bounds).clone();
        }

        // RDKit✔️✔️:       if (!EmbeddingOps::setupInitialBoundsMatrix(piece.get(), mmat, coordMap,
        // RDKit✔️✔️:                                                   params, etkdgDetails)) {
        // RDKit✔️✔️:         // return if we couldn't setup the bounds matrix
        // RDKit✔️✔️:         // possible causes include a triangle smoothing failure
        // RDKit✔️✔️:         return;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:
        // RDKit✔️✔️:     // find all the chiral centers in the molecule
        // RDKit✔️✔️:     MolOps::assignStereochemistry(*piece);
        //
        // RDKit mutates `piece` in place. Its ring info initialized during bounds
        // setup remains attached to the same fragment object and is immediately
        // visible to `findChiralSets()`. COSMolKit fragments are rebuilt value
        // objects with empty derived caches, so persist the same symmetrized ring
        // information on the local fragment before downstream ring-dependent
        // embedding logic.
        let piece_storage;
        let piece = if piece
            .derived_cache()
            .rings
            .as_ref()
            .is_some_and(crate::RingInfo::is_initialized)
        {
            piece
        } else {
            let mut ringed_piece = piece.clone();
            let ring_info = ring_info_for_distgeom(piece)?.into_owned();
            ringed_piece.derived_cache_mut().rings = Some(ring_info);
            piece_storage = ringed_piece;
            &piece_storage
        };

        // RDKit✔️✔️:     // find all the chiral centers in the molecule
        // RDKit✔️✔️:     MolOps::assignStereochemistry(*piece);
        // RDKit✔️✔️:     DistGeom::VECT_CHIRALSET chiralCenters;
        // RDKit✔️✔️:     DistGeom::VECT_CHIRALSET tetrahedralCarbons;
        // RDKit✔️✔️:     EmbeddingOps::findChiralSets(*piece, chiralCenters, tetrahedralCarbons,
        // RDKit✔️✔️:                                  coordMap);
        let mut chiral_centers = Vec::new();
        let mut tetrahedral_carbons = Vec::new();
        embedder_find_chiral_sets(
            piece,
            &mut chiral_centers,
            &mut tetrahedral_carbons,
            coord_map,
        );

        // RDKit✔️✔️:     // find double bonds
        // RDKit✔️✔️:     std::vector<std::tuple<unsigned int, unsigned int, unsigned int>>
        // RDKit✔️✔️:         doubleBondEnds;
        // RDKit✔️✔️:     std::vector<std::pair<std::vector<unsigned int>, int>> stereoDoubleBonds;
        // RDKit✔️✔️:     EmbeddingOps::findDoubleBonds(*piece, doubleBondEnds, stereoDoubleBonds,
        // RDKit✔️✔️:                                   coordMap);
        let mut double_bond_ends = Vec::new();
        let mut stereo_double_bonds = Vec::new();
        embedder_find_double_bonds(
            piece,
            &mut double_bond_ends,
            &mut stereo_double_bonds,
            coord_map,
        );

        // RDKit✔️✔️:     // if we have any chiral centers or are using random coordinates, we
        // RDKit✔️✔️:     // will first embed the molecule in four dimensions, otherwise we will
        // RDKit✔️✔️:     // use 3D
        // RDKit✔️✔️:     bool fourD = false;
        // RDKit✔️✔️:     if (params.useRandomCoords || chiralCenters.size() > 0) {
        // RDKit✔️✔️:       fourD = true;
        // RDKit✔️✔️:     }
        let four_d = params.use_random_coords || !chiral_centers.is_empty();
        // RDKit✔️✔️:     int numThreads = getNumThreadsToUse(params.numThreads);
        let num_threads = if params.num_threads <= 0 {
            std::thread::available_parallelism().map_or(1, std::num::NonZeroUsize::get) as i32
        } else {
            params.num_threads
        };

        // RDKit✔️✔️:     detail::EmbedArgs eargs = {&confsOk,        fourD,
        // RDKit✔️✔️:                                &fragMapping,    &confs,
        // RDKit✔️✔️:                                fragIdx,         mmat,
        // RDKit✔️✔️:                                &chiralCenters,  &tetrahedralCarbons,
        // RDKit✔️✔️:                                &doubleBondEnds, &stereoDoubleBonds,
        // RDKit✔️✔️:                                &etkdgDetails};
        let mut helper_args = EmbedHelperArgs {
            confs_ok: &mut confs_ok,
            four_d,
            frag_mapping: Some(&frag_mapping),
            confs: &mut confs,
            frag_idx,
            mmat: &mmat,
            chiral_centers: &chiral_centers,
            tetrahedral_carbons: &tetrahedral_carbons,
            double_bond_ends: &double_bond_ends,
            stereo_double_bonds: &stereo_double_bonds,
            etkdg_details: &etkdg_details,
        };
        // RDKit✔️✔️:     if (numThreads == 1) {
        // RDKit✔️✔️:       detail::embedHelper_(0, 1, &eargs, &params, end_time);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     else { ... dispatch thread ids ... }
        for thread_id in 0..num_threads {
            embedder_embed_helper(thread_id, num_threads, &mut helper_args, params, end_time)?;
        }

        // RDKit✔️✔️:     if (end_time != nullptr && Clock::now() > *end_time) {
        // RDKit✔️✔️:       if (params.trackFailures) {
        // RDKit✔️✔️:         params.failures[EmbedFailureCauses::EXCEEDED_TIMEOUT]++;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       res.push_back(-1);
        // RDKit✔️✔️:       return;
        // RDKit✔️✔️:     }
        if let Some(deadline) = end_time
            && Instant::now() > deadline
        {
            embedder_increment_failure(params, EmbedFailureCause::ExceededTimeout);
            let (out, mut res) = embedder_add_conformers(mol, Vec::new(), params.clear_confs)?;
            res.push(-1);
            return Ok((out, res));
        }
    }
    // RDKit✔️✔️:   std::vector<std::vector<unsigned int>> selfMatches;
    // RDKit✔️✔️:   if (params.pruneRmsThresh > 0.0) {
    // RDKit✔️✔️:     selfMatches = detail::getMolSelfMatches(mol, params);
    // RDKit✔️✔️:   }
    let self_matches = if params.prune_rms_thresh > 0.0 {
        embedder_get_mol_self_matches(mol, params)?
    } else {
        Vec::new()
    };

    // RDKit✔️✔️:   for (unsigned int ci = 0; ci < confs.size(); ++ci) {
    // RDKit✔️✔️:     auto &conf = confs[ci];
    // RDKit✔️✔️:     if (confsOk[ci]) {
    // RDKit✔️✔️:       // check if we are pruning away conformations and
    // RDKit✔️✔️:       // a close-by conformation has already been chosen :
    // RDKit✔️✔️:       if (params.pruneRmsThresh <= 0.0 ||
    // RDKit✔️✔️:           _isConfFarFromRest(mol, *conf, params.pruneRmsThresh, selfMatches)) {
    // RDKit✔️✔️:         int confId = (int)mol.addConformer(conf.release(), true);
    // RDKit✔️✔️:         res.push_back(confId);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    let mut out = mol.clone();
    if params.clear_confs {
        out.coordinate_block_mut().conformers_3d.clear();
    }
    let mut res = Vec::new();
    for (ci, conf) in confs.into_iter().enumerate() {
        if confs_ok[ci]
            && (params.prune_rms_thresh <= 0.0
                || embedder_is_conf_far_from_rest(
                    &out,
                    &conf,
                    params.prune_rms_thresh,
                    &self_matches,
                ))
        {
            let conf_id = out.coordinate_block_mut().conformers_3d.len();
            out.coordinate_block_mut()
                .conformers_3d
                .push(Conformer3D::new(
                    conf_id,
                    conf.coordinates().to_vec(),
                    conf.is_3d(),
                ));
            res.push(conf_id as i32);
        }
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DGeomHelpers::EmbedMultipleConfs
    Ok((out, res))
}

pub fn embed_multiple_confs_return_vector(
    mol: &Molecule,
    num_confs: u32,
    params: &mut EmbedParameters,
) -> Result<Vec<i32>, DgBoundsError> {
    // BEGIN RDKIT CPP INLINE FUNCTION DGeomHelpers::EmbedMultipleConfs return-vector overload (Embedder.h:216-221)
    // RDKit✔️✔️: inline INT_VECT EmbedMultipleConfs(ROMol &mol, unsigned int numConfs,
    // RDKit✔️✔️:                                    EmbedParameters &params) {
    // RDKit✔️✔️:   INT_VECT res;
    // RDKit✔️✔️:   EmbedMultipleConfs(mol, res, numConfs, params);
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP INLINE FUNCTION DGeomHelpers::EmbedMultipleConfs return-vector overload
    let (_out, res) = embed_multiple_confs(mol, num_confs, params)?;
    Ok(res)
}

pub fn embed_multiple_confs_result(
    mol: &Molecule,
    num_confs: u32,
    params: &mut EmbedParameters,
) -> Result<EmbedMultipleConfsResult, DgBoundsError> {
    let (molecule, conf_ids) = embed_multiple_confs(mol, num_confs, params)?;
    Ok(EmbedMultipleConfsResult {
        molecule,
        conf_ids,
        requested_num_confs: num_confs,
        params: params.clone(),
    })
}

pub fn embed_molecule(
    mol: &Molecule,
    params: &mut EmbedParameters,
) -> Result<(Molecule, i32), DgBoundsError> {
    // BEGIN RDKIT CPP INLINE FUNCTION DGeomHelpers::EmbedMolecule(ROMol&, EmbedParameters&) (Embedder.h:225-236)
    // RDKit✔️✔️: inline int EmbedMolecule(ROMol &mol, EmbedParameters &params) {
    // RDKit✔️✔️:   INT_VECT confIds;
    // RDKit✔️✔️:   EmbedMultipleConfs(mol, confIds, 1, params);
    let (out, conf_ids) = embed_multiple_confs(mol, 1, params)?;
    // RDKit✔️✔️:   int res;
    // RDKit✔️✔️:   if (confIds.size()) {
    // RDKit✔️✔️:     res = confIds[0];
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     res = -1;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP INLINE FUNCTION DGeomHelpers::EmbedMolecule
    Ok((out, conf_ids.first().copied().unwrap_or(-1)))
}

pub fn embed_molecule_result(
    mol: &Molecule,
    params: &mut EmbedParameters,
) -> Result<EmbedMoleculeResult, DgBoundsError> {
    let (molecule, conf_id) = embed_molecule(mol, params)?;
    Ok(EmbedMoleculeResult {
        molecule,
        conf_id,
        params: params.clone(),
    })
}

#[allow(clippy::too_many_arguments)]
pub fn rd_distgeom_embed_molecule_wrapper(
    mol: &Molecule,
    max_attempts: u32,
    seed: i32,
    clear_confs: bool,
    use_random_coords: bool,
    box_size_mult: f64,
    rand_neg_eig: bool,
    num_zero_fail: u32,
    coord_map: BTreeMap<i32, ForceFieldVec3>,
    force_tol: f64,
    ignore_smoothing_failures: bool,
    enforce_chirality: bool,
    use_exp_torsion_angle_prefs: bool,
    use_basic_knowledge: bool,
    print_exp_torsion_angles: bool,
    use_small_ring_torsions: bool,
    use_macrocycle_torsions: bool,
    et_version: u32,
    use_macrocycle14config: bool,
) -> Result<(Molecule, i32), DgBoundsError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::EmbedMolecule wrapper (rdDistGeom.cpp:107-150)
    // RDKit✔️✔️: int EmbedMolecule(ROMol &mol, unsigned int maxAttempts, int seed,
    // RDKit✔️✔️:                   bool clearConfs, bool useRandomCoords, double boxSizeMult,
    // RDKit✔️✔️:                   bool randNegEig, unsigned int numZeroFail,
    // RDKit✔️✔️:                   python::dict &coordMap, double forceTol,
    // RDKit✔️✔️:                   bool ignoreSmoothingFailures, bool enforceChirality,
    // RDKit✔️✔️:                   bool useExpTorsionAnglePrefs, bool useBasicKnowledge,
    // RDKit✔️✔️:                   bool printExpTorsionAngles, bool useSmallRingTorsions,
    // RDKit✔️✔️:                   bool useMacrocycleTorsions, unsigned int ETversion,
    // RDKit✔️✔️:                   bool useMacrocycle14config) {
    // RDKit✔️✔️:   bool verbose = printExpTorsionAngles;
    // RDKit✔️✔️:   int numThreads = 1;
    // RDKit✔️✔️:   double pruneRmsThresh = -1.;
    // RDKit✔️✔️:   const double basinThresh = DGeomHelpers::EmbedParameters().basinThresh;
    // RDKit✔️✔️:   bool onlyHeavyAtomsForRMS = false;
    let mut params = EmbedParameters::from_rdkit_constructor(
        max_attempts,
        1,
        seed,
        clear_confs,
        use_random_coords,
        box_size_mult,
        rand_neg_eig,
        num_zero_fail,
        (!coord_map.is_empty()).then_some(coord_map),
        force_tol,
        ignore_smoothing_failures,
        enforce_chirality,
        use_exp_torsion_angle_prefs,
        use_basic_knowledge,
        print_exp_torsion_angles,
        EmbedParameters::default().basin_thresh,
        -1.0,
        false,
        et_version,
        None,
        true,
        use_small_ring_torsions,
        use_macrocycle_torsions,
        use_macrocycle14config,
        0,
        None,
        None,
    );
    // RDKit✔️✔️:   res = DGeomHelpers::EmbedMolecule(mol, params);
    // END RDKIT CPP FUNCTION RDKit::EmbedMolecule wrapper
    embed_molecule(mol, &mut params)
}

pub fn rd_distgeom_embed_molecule_wrapper_with_params(
    mol: &Molecule,
    params: &mut EmbedParameters,
) -> Result<(Molecule, i32), DgBoundsError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::EmbedMolecule2 wrapper (rdDistGeom.cpp:152-162)
    // RDKit✔️✔️: int EmbedMolecule2(ROMol &mol, DGeomHelpers::EmbedParameters &params) {
    // RDKit✔️✔️:   int res;
    // RDKit✔️✔️:   {
    // RDKit✔️✔️:     NOGIL gil;
    // RDKit✔️✔️:     res = DGeomHelpers::EmbedMolecule(mol, params);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RDKit::EmbedMolecule2 wrapper
    embed_molecule(mol, params)
}

pub fn rd_distgeom_get_kdg() -> EmbedParameters {
    // BEGIN RDKIT CPP FUNCTION RDKit::getKDG/getETDG/getETDGv2/getETKDG/getETKDGv2/getETKDGv3/getsrETKDGv3 (rdDistGeom.cpp:263-269)
    // RDKit✔️✔️: PyEmbedParameters *getKDG() { return new PyEmbedParameters(DGeomHelpers::KDG); }
    // RDKit✔️✔️: PyEmbedParameters *getETDG() { return new PyEmbedParameters(DGeomHelpers::ETDG); }
    // RDKit✔️✔️: PyEmbedParameters *getETDGv2() { return new PyEmbedParameters(DGeomHelpers::ETDGv2); }
    // RDKit✔️✔️: PyEmbedParameters *getETKDG() { return new PyEmbedParameters(DGeomHelpers::ETKDG); }
    // RDKit✔️✔️: PyEmbedParameters *getETKDGv2() { return new PyEmbedParameters(DGeomHelpers::ETKDGv2); }
    // RDKit✔️✔️: PyEmbedParameters *getETKDGv3() { return new PyEmbedParameters(DGeomHelpers::ETKDGv3); }
    // RDKit✔️✔️: PyEmbedParameters *getsrETKDGv3() { return new PyEmbedParameters(DGeomHelpers::srETKDGv3); }
    // END RDKIT CPP FUNCTION RDKit parameter factories
    EmbedParameters::kdg()
}

pub fn rd_distgeom_get_etdg() -> EmbedParameters {
    EmbedParameters::etdg()
}

pub fn rd_distgeom_get_etdg_v2() -> EmbedParameters {
    EmbedParameters::etdg_v2()
}

pub fn rd_distgeom_get_etkdg() -> EmbedParameters {
    EmbedParameters::etkdg()
}

pub fn rd_distgeom_get_etkdg_v2() -> EmbedParameters {
    EmbedParameters::etkdg_v2()
}

pub fn rd_distgeom_get_etkdg_v3() -> EmbedParameters {
    EmbedParameters::etkdg_v3()
}

pub fn rd_distgeom_get_sr_etkdg_v3() -> EmbedParameters {
    EmbedParameters::sr_etkdg_v3()
}

pub fn rd_distgeom_get_exp_tors_helper(
    mol: &Molecule,
    use_exp_torsions: bool,
    use_small_ring_torsions: bool,
    use_macrocycle_torsions: bool,
    use_basic_knowledge: bool,
    et_version: u32,
    verbose: bool,
) -> Result<CrystalFFDetails, DgBoundsError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::getExpTorsHelper (rdDistGeom.cpp:271-295)
    // RDKit✔️✔️: python::tuple getExpTorsHelper(const RDKit::ROMol &mol, bool useExpTorsions,
    // RDKit✔️✔️:                               bool useSmallRingTorsions,
    // RDKit✔️✔️:                               bool useMacrocycleTorsions,
    // RDKit✔️✔️:                               bool useBasicKnowledge, unsigned int ETversion,
    // RDKit✔️✔️:                               bool verbose) {
    // RDKit✔️✔️:   RDKit::ForceFields::CrystalFF::CrystalFFDetails details;
    // RDKit✔️✔️:   RDKit::ForceFields::CrystalFF::getExperimentalTorsions(
    // RDKit✔️✔️:       mol, details, useExpTorsions, useSmallRingTorsions,
    // RDKit✔️✔️:       useMacrocycleTorsions, useBasicKnowledge, ETversion, verbose);
    // END RDKIT CPP FUNCTION RDKit::getExpTorsHelper
    let mut details = CrystalFFDetails::default();
    get_experimental_torsions_without_bonds(
        mol,
        &mut details,
        use_exp_torsions,
        use_small_ring_torsions,
        use_macrocycle_torsions,
        use_basic_knowledge,
        et_version,
        verbose,
    )
    .map_err(|err| DgBoundsError::GenerationFailed(err.to_string()))?;
    Ok(details)
}

pub fn rd_distgeom_get_exp_tors_helper_with_params(
    mol: &Molecule,
    params: &EmbedParameters,
) -> Result<CrystalFFDetails, DgBoundsError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::getExpTorsHelperWithParams (rdDistGeom.cpp:297-302)
    // RDKit✔️✔️: python::tuple getExpTorsHelperWithParams(
    // RDKit✔️✔️:     const RDKit::ROMol &mol, const DGeomHelpers::EmbedParameters &ps) {
    // RDKit✔️✔️:   return getExpTorsHelper(mol, ps.useExpTorsionAnglePrefs,
    // RDKit✔️✔️:                           ps.useSmallRingTorsions, ps.useMacrocycleTorsions,
    // RDKit✔️✔️:                           ps.useBasicKnowledge, ps.ETversion, ps.verbose);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RDKit::getExpTorsHelperWithParams
    rd_distgeom_get_exp_tors_helper(
        mol,
        params.use_exp_torsion_angle_prefs,
        params.use_small_ring_torsions,
        params.use_macrocycle_torsions,
        params.use_basic_knowledge,
        params.et_version,
        params.verbose,
    )
}

pub fn rd_distgeom_embed_parameters_to_json_helper(params: &EmbedParameters) -> String {
    // BEGIN RDKIT CPP FUNCTION RDKit::embedParametersToJSONHelper (rdDistGeom.cpp:304-307)
    // RDKit✔️✔️: python::str embedParametersToJSONHelper(
    // RDKit✔️✔️:     const DGeomHelpers::EmbedParameters &ps) {
    // RDKit✔️✔️:   return python::str(embedParametersToJSON(ps));
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RDKit::embedParametersToJSONHelper
    params.to_json()
}

// ──────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────

#[cfg(test)]
mod tests;
