//! Complete private numeric phases of the pinned detached embedding flow.
use crate::{
    ConformerError, EmbedParams,
    bounds::{BoundsMatrix, BoundsMatrixError},
    numeric::{
        RdkitDoubleRng, SymmMatrix, compute_initial_coords_with_rng,
        compute_random_coords_with_rng, pick_random_dist_mat_with_rng,
    },
};
use cosmolkit_forcefields::{
    ChiralSet, ChiralSetPtr, ChiralSetStructureFlags, ConformerOptimizer, ConformerOptimizerError,
    CrystalFFDetails, DistanceBoundsRead, DistanceGeometryForceFieldParams,
    calc_chiral_volume_rows,
};
use std::collections::BTreeMap;
use web_time::Instant;
const EMBEDDER_ERROR_TOL: f64 = 0.00001;
const MAX_MINIMIZED_E_PER_ATOM: f64 = 0.05;
const MIN_TETRAHEDRAL_CHIRAL_VOL: f64 = 0.50;
const TETRAHEDRAL_CENTERINVOLUME_TOL: f64 = 0.30;
#[derive(Debug, thiserror::Error)]
pub enum GenerationError {
    #[error(transparent)]
    Topology(#[from] cosmolkit_model::TopologyValidationError),
    #[error(transparent)]
    Coordinates(#[from] cosmolkit_model::CoordinateValidationError),
    #[error(transparent)]
    Valence(#[from] cosmolkit_core::ValenceError),
    #[error(transparent)]
    Rings(#[from] cosmolkit_core::RingFindingError),
    #[error(transparent)]
    Conjugation(#[from] cosmolkit_core::ConjugationError),
    #[error(transparent)]
    Hybridization(#[from] cosmolkit_core::HybridizationError),
    #[error(transparent)]
    Stereo(#[from] cosmolkit_core::LegacyStereoError),
    #[error(transparent)]
    Fragments(#[from] cosmolkit_core::MoleculeFragmentsError),
    #[error(transparent)]
    Hydrogens(#[from] cosmolkit_forcefields::MissingExplicitHydrogensError),
    #[error(transparent)]
    Pruning(#[from] crate::PruningError),

    #[error(transparent)]
    Seed(#[from] crate::numeric::EmbedSeedError),
    #[error(transparent)]
    GenerationFailed(#[from] GenerationFailure),
    #[error(transparent)]
    Atropisomer(#[from] cosmolkit_core::AtropisomerError),
    #[error(transparent)]
    Bounds(#[from] BoundsMatrixError),
    #[error(transparent)]
    ForceField(#[from] ConformerOptimizerError),
    #[error(transparent)]
    Numeric(#[from] ConformerError),
    #[error(transparent)]
    ThreadCount(#[from] cosmolkit_core::ThreadCountError),
    #[error("embedding worker launch failed: {0}")]
    WorkerLaunch(#[source] std::io::Error),
    #[error("{0}")]
    Input(&'static str),
}
// The optimizer checks all pair indices against dimension before it reads bounds.
// BoundsMatrix is the sole existing storage; this narrow reader owns no matrix copy.
impl DistanceBoundsRead for BoundsMatrix {
    fn dimension(&self) -> usize {
        self.dimension()
    }
    fn get_lower(&self, i: usize, j: usize) -> f64 {
        BoundsMatrix::get_lower(self, i, j).expect("optimizer validated bounds index")
    }
    fn get_upper(&self, i: usize, j: usize) -> f64 {
        BoundsMatrix::get_upper(self, i, j).expect("optimizer validated bounds index")
    }
}
fn point_at(positions: &[impl AsRef<[f64]>], i: usize) -> Result<[f64; 3], GenerationError> {
    let row = positions
        .get(i)
        .ok_or(GenerationError::Input("point index out of range"))?
        .as_ref();
    if row.len() < 3 {
        return Err(GenerationError::Input("point dimension less than three"));
    }
    Ok([row[0], row[1], row[2]])
}
fn delta(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    // BEGIN RDKIT CPP FUNCTION delta (Geometry/point.cpp:64-70)
    // RDKit❗✔️: Point3D operator-(const Point3D &p1, const Point3D &p2) {
    // RDKit❗✔️:   Point3D res;
    // RDKit❗✔️:   res.x = p1.x - p2.x;
    // RDKit❗✔️:   res.y = p1.y - p2.y;
    // RDKit❗✔️:   res.z = p1.z - p2.z;
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION delta
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    // BEGIN RDKIT CPP FUNCTION dot (Geometry/point.h:169-172)
    // RDKit❗✔️:   constexpr double dotProduct(const Point3D &other) const {
    // RDKit❗✔️:     double res = x * (other.x) + y * (other.y) + z * (other.z);
    // RDKit❗✔️:     return res;
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION dot
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    // BEGIN RDKIT CPP FUNCTION cross (Geometry/point.h:228-235)
    // RDKit❗✔️:   constexpr Point3D crossProduct(const Point3D &other) const {
    // RDKit❗✔️:     Point3D res;
    // RDKit❗✔️:     res.x = y * (other.z) - z * (other.y);
    // RDKit❗✔️:     res.y = -x * (other.z) + z * (other.x);
    // RDKit❗✔️:     res.z = x * (other.y) - y * (other.x);
    // RDKit❗✔️:     return res;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // END RDKIT CPP FUNCTION cross
    [
        a[1] * b[2] - a[2] * b[1],
        -a[0] * b[2] + a[2] * b[0],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn length(a: [f64; 3]) -> f64 {
    // BEGIN RDKIT CPP FUNCTION length (Geometry/point.h:158-161)
    // RDKit❗✔️:   double length() const override {
    // RDKit❗✔️:     double res = x * x + y * y + z * z;
    // RDKit❗✔️:     return sqrt(res);
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION length
    dot(a, a).sqrt()
}
fn normalized_point_delta(p0: [f64; 3], p1: [f64; 3]) -> Result<[f64; 3], GenerationError> {
    // BEGIN RDKIT CPP FUNCTION normalized_point_delta (Geometry/point.h:147-156)
    // RDKit❗✔️:   constexpr void normalize() override {
    // RDKit❗✔️:     double l = this->length();
    // RDKit❗✔️:     if (l < zero_tolerance) {
    // RDKit❗✔️:       throw std::runtime_error("Cannot normalize a zero length vector");
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     x /= l;
    // RDKit❗✔️:     y /= l;
    // RDKit❗✔️:     z /= l;
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION normalized_point_delta

    let v = delta(p0, p1);
    let len = length(v);
    if len < 1e-16 {
        return Err(GenerationError::Input(
            "Cannot normalize a zero length vector",
        ));
    }
    Ok([v[0] / len, v[1] / len, v[2] / len])
}
fn embedder_volume_test(
    chiral_set: &ChiralSet,
    positions: &[impl AsRef<[f64]>],
) -> Result<bool, GenerationError> {
    // BEGIN RDKIT CPP FUNCTION embedder_volume_test (GraphMol/DistGeomHelpers/Embedder.cpp:329-397)
    // RDKit❗✔️: bool _volumeTest(const DistGeom::ChiralSetPtr &chiralSet,
    // RDKit❗✔️:                  const RDGeom::PointPtrVect &positions, bool verbose = false) {
    // RDKit❗✔️:   RDGeom::Point3D p0((*positions[chiralSet->d_idx0])[0],
    // RDKit❗✔️:                      (*positions[chiralSet->d_idx0])[1],
    // RDKit❗✔️:                      (*positions[chiralSet->d_idx0])[2]);
    // RDKit❗✔️:   RDGeom::Point3D p1((*positions[chiralSet->d_idx1])[0],
    // RDKit❗✔️:                      (*positions[chiralSet->d_idx1])[1],
    // RDKit❗✔️:                      (*positions[chiralSet->d_idx1])[2]);
    // RDKit❗✔️:   RDGeom::Point3D p2((*positions[chiralSet->d_idx2])[0],
    // RDKit❗✔️:                      (*positions[chiralSet->d_idx2])[1],
    // RDKit❗✔️:                      (*positions[chiralSet->d_idx2])[2]);
    // RDKit❗✔️:   RDGeom::Point3D p3((*positions[chiralSet->d_idx3])[0],
    // RDKit❗✔️:                      (*positions[chiralSet->d_idx3])[1],
    // RDKit❗✔️:                      (*positions[chiralSet->d_idx3])[2]);
    // RDKit❗✔️:   RDGeom::Point3D p4((*positions[chiralSet->d_idx4])[0],
    // RDKit❗✔️:                      (*positions[chiralSet->d_idx4])[1],
    // RDKit❗✔️:                      (*positions[chiralSet->d_idx4])[2]);
    // RDKit❗✔️:
    // RDKit❗✔️:   // even if we are minimizing in higher dimension the chiral volume is
    // RDKit❗✔️:   // calculated using only the first 3 dimensions
    // RDKit❗✔️:   RDGeom::Point3D v1 = p0 - p1;
    // RDKit❗✔️:   v1.normalize();
    // RDKit❗✔️:   RDGeom::Point3D v2 = p0 - p2;
    // RDKit❗✔️:   v2.normalize();
    // RDKit❗✔️:   RDGeom::Point3D v3 = p0 - p3;
    // RDKit❗✔️:   v3.normalize();
    // RDKit❗✔️:   RDGeom::Point3D v4 = p0 - p4;
    // RDKit❗✔️:   v4.normalize();
    // RDKit❗✔️:
    // RDKit❗✔️:   // be more tolerant of tethrahedral centers that are involved in multiple
    // RDKit❗✔️:   // small rings
    // RDKit❗✔️:   double volScale = 1;
    // RDKit❗✔️:   if (chiralSet->d_structureFlags &
    // RDKit❗✔️:       static_cast<std::uint64_t>(
    // RDKit❗✔️:           DistGeom::ChiralSetStructureFlags::IN_FUSED_SMALL_RINGS)) {
    // RDKit❗✔️:     volScale = 0.25;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   RDGeom::Point3D crossp = v1.crossProduct(v2);
    // RDKit❗✔️:   double vol = crossp.dotProduct(v3);
    // RDKit❗✔️:   if (verbose) {
    // RDKit❗✔️:     std::cerr << "   " << fabs(vol) << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (fabs(vol) < volScale * MIN_TETRAHEDRAL_CHIRAL_VOL) {
    // RDKit❗✔️:     return false;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   crossp = v1.crossProduct(v2);
    // RDKit❗✔️:   vol = crossp.dotProduct(v4);
    // RDKit❗✔️:   if (verbose) {
    // RDKit❗✔️:     std::cerr << "   " << fabs(vol) << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (fabs(vol) < volScale * MIN_TETRAHEDRAL_CHIRAL_VOL) {
    // RDKit❗✔️:     return false;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   crossp = v1.crossProduct(v3);
    // RDKit❗✔️:   vol = crossp.dotProduct(v4);
    // RDKit❗✔️:   if (verbose) {
    // RDKit❗✔️:     std::cerr << "   " << fabs(vol) << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (fabs(vol) < volScale * MIN_TETRAHEDRAL_CHIRAL_VOL) {
    // RDKit❗✔️:     return false;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   crossp = v2.crossProduct(v3);
    // RDKit❗✔️:   vol = crossp.dotProduct(v4);
    // RDKit❗✔️:   if (verbose) {
    // RDKit❗✔️:     std::cerr << "   " << fabs(vol) << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return fabs(vol) >= volScale * MIN_TETRAHEDRAL_CHIRAL_VOL;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION embedder_volume_test

    let p0 = point_at(positions, chiral_set.idx0)?;
    let v1 = normalized_point_delta(p0, point_at(positions, chiral_set.idx1)?)?;
    let v2 = normalized_point_delta(p0, point_at(positions, chiral_set.idx2)?)?;
    let v3 = normalized_point_delta(p0, point_at(positions, chiral_set.idx3)?)?;
    let v4 = normalized_point_delta(p0, point_at(positions, chiral_set.idx4)?)?;
    let vol_scale =
        if chiral_set.structure_flags & ChiralSetStructureFlags::InFusedSmallRings as u64 != 0 {
            0.25
        } else {
            1.0
        };
    let min_vol = vol_scale * MIN_TETRAHEDRAL_CHIRAL_VOL;
    if dot(cross(v1, v2), v3).abs() < min_vol {
        return Ok(false);
    }
    if dot(cross(v1, v2), v4).abs() < min_vol {
        return Ok(false);
    }
    if dot(cross(v1, v3), v4).abs() < min_vol {
        return Ok(false);
    }
    Ok(dot(cross(v2, v3), v4).abs() >= min_vol)
}
fn embedder_same_side(
    v1: [f64; 3],
    v2: [f64; 3],
    v3: [f64; 3],
    v4: [f64; 3],
    p0: [f64; 3],
    tol: f64,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION embedder_same_side (GraphMol/DistGeomHelpers/Embedder.cpp:399-410)
    // RDKit❗✔️: bool _sameSide(const RDGeom::Point3D &v1, const RDGeom::Point3D &v2,
    // RDKit❗✔️:                const RDGeom::Point3D &v3, const RDGeom::Point3D &v4,
    // RDKit❗✔️:                const RDGeom::Point3D &p0, double tol = 0.1) {
    // RDKit❗✔️:   RDGeom::Point3D normal = (v2 - v1).crossProduct(v3 - v1);
    // RDKit❗✔️:   double d1 = normal.dotProduct(v4 - v1);
    // RDKit❗✔️:   double d2 = normal.dotProduct(p0 - v1);
    // RDKit❗✔️:   // std::cerr << "     " << d1 << " - " << d2 << std::endl;
    // RDKit❗✔️:   if (fabs(d1) < tol || fabs(d2) < tol) {
    // RDKit❗✔️:     return false;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return !((d1 < 0.) ^ (d2 < 0.));
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION embedder_same_side

    let normal = cross(delta(v2, v1), delta(v3, v1));
    let d1 = dot(normal, delta(v4, v1));
    let d2 = dot(normal, delta(p0, v1));
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
    positions: &[impl AsRef<[f64]>],
    tol: f64,
) -> Result<bool, GenerationError> {
    // BEGIN RDKIT CPP FUNCTION embedder_center_in_volume_indices (GraphMol/DistGeomHelpers/Embedder.cpp:411-437)
    // RDKit❗✔️: bool _centerInVolume(unsigned int idx0, unsigned int idx1, unsigned int idx2,
    // RDKit❗✔️:                      unsigned int idx3, unsigned int idx4,
    // RDKit❗✔️:                      const RDGeom::PointPtrVect &positions, double tol,
    // RDKit❗✔️:                      bool verbose = false) {
    // RDKit❗✔️:   RDGeom::Point3D p0((*positions[idx0])[0], (*positions[idx0])[1],
    // RDKit❗✔️:                      (*positions[idx0])[2]);
    // RDKit❗✔️:   RDGeom::Point3D p1((*positions[idx1])[0], (*positions[idx1])[1],
    // RDKit❗✔️:                      (*positions[idx1])[2]);
    // RDKit❗✔️:   RDGeom::Point3D p2((*positions[idx2])[0], (*positions[idx2])[1],
    // RDKit❗✔️:                      (*positions[idx2])[2]);
    // RDKit❗✔️:   RDGeom::Point3D p3((*positions[idx3])[0], (*positions[idx3])[1],
    // RDKit❗✔️:                      (*positions[idx3])[2]);
    // RDKit❗✔️:   RDGeom::Point3D p4((*positions[idx4])[0], (*positions[idx4])[1],
    // RDKit❗✔️:                      (*positions[idx4])[2]);
    // RDKit❗✔️:   // RDGeom::Point3D centroid = (p1+p2+p3+p4)/4.;
    // RDKit❗✔️:   if (verbose) {
    // RDKit❗✔️:     std::cerr << _sameSide(p1, p2, p3, p4, p0, tol) << " "
    // RDKit❗✔️:               << _sameSide(p2, p3, p4, p1, p0, tol) << " "
    // RDKit❗✔️:               << _sameSide(p3, p4, p1, p2, p0, tol) << " "
    // RDKit❗✔️:               << _sameSide(p4, p1, p2, p3, p0, tol) << std::endl;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   bool res = _sameSide(p1, p2, p3, p4, p0, tol) &&
    // RDKit❗✔️:              _sameSide(p2, p3, p4, p1, p0, tol) &&
    // RDKit❗✔️:              _sameSide(p3, p4, p1, p2, p0, tol) &&
    // RDKit❗✔️:              _sameSide(p4, p1, p2, p3, p0, tol);
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION embedder_center_in_volume_indices

    let p0 = point_at(positions, idx0)?;
    let p1 = point_at(positions, idx1)?;
    let p2 = point_at(positions, idx2)?;
    let p3 = point_at(positions, idx3)?;
    let p4 = point_at(positions, idx4)?;
    Ok(embedder_same_side(p1, p2, p3, p4, p0, tol)
        && embedder_same_side(p2, p3, p4, p1, p0, tol)
        && embedder_same_side(p3, p4, p1, p2, p0, tol)
        && embedder_same_side(p4, p1, p2, p3, p0, tol))
}
fn embedder_center_in_volume(
    chiral_set: &ChiralSet,
    positions: &[impl AsRef<[f64]>],
    tol: f64,
) -> Result<bool, GenerationError> {
    // BEGIN RDKIT CPP FUNCTION embedder_center_in_volume (GraphMol/DistGeomHelpers/Embedder.cpp:439-449)
    // RDKit❗✔️: bool _centerInVolume(const DistGeom::ChiralSetPtr &chiralSet,
    // RDKit❗✔️:                      const RDGeom::PointPtrVect &positions, double tol = 0.1,
    // RDKit❗✔️:                      bool verbose = false) {
    // RDKit❗✔️:   if (chiralSet->d_idx0 ==
    // RDKit❗✔️:       chiralSet->d_idx4) {  // this happens for three-coordinate centers
    // RDKit❗✔️:     return true;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return _centerInVolume(chiralSet->d_idx0, chiralSet->d_idx1,
    // RDKit❗✔️:                          chiralSet->d_idx2, chiralSet->d_idx3,
    // RDKit❗✔️:                          chiralSet->d_idx4, positions, tol, verbose);
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION embedder_center_in_volume

    if chiral_set.idx0 == chiral_set.idx4 {
        return Ok(true);
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
    positions: &[impl AsRef<[f64]>],
) -> Result<bool, GenerationError> {
    // BEGIN RDKIT CPP FUNCTION embedder_bounds_fulfilled (GraphMol/DistGeomHelpers/Embedder.cpp:450-478)
    // RDKit❗✔️: bool _boundsFulfilled(const std::vector<int> &atoms,
    // RDKit❗✔️:                       const DistGeom::BoundsMatrix &mmat,
    // RDKit❗✔️:                       const RDGeom::PointPtrVect &positions) {
    // RDKit❗✔️:   // unsigned int N = mmat.numRows();
    // RDKit❗✔️:   // std::cerr << N << " " << atoms.size() << std::endl;
    // RDKit❗✔️:   // loop over all pair of atoms
    // RDKit❗✔️:   for (unsigned int i = 0; i < atoms.size() - 1; ++i) {
    // RDKit❗✔️:     for (unsigned int j = i + 1; j < atoms.size(); ++j) {
    // RDKit❗✔️:       int a1 = atoms[i];
    // RDKit❗✔️:       int a2 = atoms[j];
    // RDKit❗✔️:       RDGeom::Point3D p0((*positions[a1])[0], (*positions[a1])[1],
    // RDKit❗✔️:                          (*positions[a1])[2]);
    // RDKit❗✔️:       RDGeom::Point3D p1((*positions[a2])[0], (*positions[a2])[1],
    // RDKit❗✔️:                          (*positions[a2])[2]);
    // RDKit❗✔️:       double d2 = (p0 - p1).length();  // distance
    // RDKit❗✔️:       double lb = mmat.getLowerBound(a1, a2);
    // RDKit❗✔️:       double ub = mmat.getUpperBound(a1, a2);  // bounds
    // RDKit❗✔️:       if (((d2 < lb) && (fabs(d2 - lb) > 0.1 * ub)) ||
    // RDKit❗✔️:           ((d2 > ub) && (fabs(d2 - ub) > 0.1 * ub))) {
    // RDKit❗✔️: #ifdef DEBUG_EMBEDDING
    // RDKit❗✔️:         std::cerr << a1 << " " << a2 << ":" << d2 << " " << lb << " " << ub
    // RDKit❗✔️:                   << " " << fabs(d2 - lb) << " " << fabs(d2 - ub) << std::endl;
    // RDKit❗✔️: #endif
    // RDKit❗✔️:         return false;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return true;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION embedder_bounds_fulfilled

    // ROOT decision conformer-private-empty-bounds-source-precondition.json:
    // retain original empty-list condition as native extension, without RDKit parity
    // credit; finalChiralChecks must preserve its primary nonempty caller guard.
    if atoms.len() < 2 {
        return Ok(true);
    }
    for i in 0..atoms.len() - 1 {
        for j in i + 1..atoms.len() {
            let a1 = usize::try_from(atoms[i])
                .map_err(|_| GenerationError::Input("negative atom index"))?;
            let a2 = usize::try_from(atoms[j])
                .map_err(|_| GenerationError::Input("negative atom index"))?;
            let d2 = length(delta(point_at(positions, a1)?, point_at(positions, a2)?));
            let lb = mmat.get_lower(a1, a2)?;
            let ub = mmat.get_upper(a1, a2)?;
            if (d2 < lb && (d2 - lb).abs() > 0.1 * ub) || (d2 > ub && (d2 - ub).abs() > 0.1 * ub) {
                return Ok(false);
            }
        }
    }
    Ok(true)
}
struct EmbedArgs<'a> {
    mmat: &'a BoundsMatrix,
    chiral_centers: &'a [ChiralSetPtr],
    tetrahedral_carbons: &'a [ChiralSetPtr],
    etkdg_details: Option<&'a CrystalFFDetails>,
    double_bond_ends: Option<&'a [(usize, usize, usize)]>,
    stereo_double_bonds: &'a [(Vec<usize>, i32)],
}
fn coordinate_map_index(i: i32, n: usize) -> Result<usize, GenerationError> {
    let i =
        usize::try_from(i).map_err(|_| GenerationError::Input("negative coordinate map index"))?;
    if i >= n {
        return Err(GenerationError::Input("coordinate map index out of range"));
    }
    Ok(i)
}
fn fixed_points(params: &EmbedParams, n: usize) -> Result<Vec<usize>, GenerationError> {
    if params.use_random_coords {
        if let Some(map) = &params.coord_map {
            return map.keys().map(|&i| coordinate_map_index(i, n)).collect();
        }
    }
    Ok(Vec::new())
}
fn embedder_generate_initial_coords<R: RdkitDoubleRng>(
    positions: &mut [Vec<f64>],
    eargs: &EmbedArgs<'_>,
    embed_params: &EmbedParams,
    dist_mat: &mut SymmMatrix,
    rng: &mut R,
) -> Result<bool, GenerationError> {
    // BEGIN RDKIT CPP FUNCTION embedder_generate_initial_coords (GraphMol/DistGeomHelpers/Embedder.cpp:481-517)
    // RDKit❗✔️: bool generateInitialCoords(RDGeom::PointPtrVect *positions,
    // RDKit❗✔️:                            const detail::EmbedArgs &eargs,
    // RDKit❗✔️:                            const EmbedParameters &embedParams,
    // RDKit❗✔️:                            RDNumeric::DoubleSymmMatrix &distMat,
    // RDKit❗✔️:                            RDKit::double_source_type *rng) {
    // RDKit❗✔️:   bool gotCoords = false;
    // RDKit❗✔️:   if (!embedParams.useRandomCoords) {
    // RDKit❗✔️:     double largestDistance =
    // RDKit❗✔️:         DistGeom::pickRandomDistMat(*eargs.mmat, distMat, *rng);
    // RDKit❗✔️:     RDUNUSED_PARAM(largestDistance);
    // RDKit❗✔️:     gotCoords = DistGeom::computeInitialCoords(distMat, *positions, *rng,
    // RDKit❗✔️:                                                embedParams.randNegEig,
    // RDKit❗✔️:                                                embedParams.numZeroFail);
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     double boxSize;
    // RDKit❗✔️:     if (embedParams.boxSizeMult > 0) {
    // RDKit❗✔️:       boxSize = 5. * embedParams.boxSizeMult;
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       boxSize = -1 * embedParams.boxSizeMult;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     gotCoords = DistGeom::computeRandomCoords(*positions, boxSize, *rng);
    // RDKit❗✔️:     if (embedParams.useRandomCoords && embedParams.coordMap != nullptr) {
    // RDKit❗✔️:       for (const auto &v : *embedParams.coordMap) {
    // RDKit❗✔️:         auto p = positions->at(v.first);
    // RDKit❗✔️:         for (unsigned int ci = 0; ci < v.second.dimension(); ++ci) {
    // RDKit❗✔️:           (*p)[ci] = v.second[ci];
    // RDKit❗✔️:         }
    // RDKit❗✔️:         // zero out any higher dimensional components:
    // RDKit❗✔️:         for (unsigned int ci = v.second.dimension(); ci < p->dimension();
    // RDKit❗✔️:              ++ci) {
    // RDKit❗✔️:           (*p)[ci] = 0.0;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return gotCoords;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION embedder_generate_initial_coords

    if !embed_params.use_random_coords {
        let _largest_distance = pick_random_dist_mat_with_rng(eargs.mmat, dist_mat, rng);
        Ok(compute_initial_coords_with_rng(
            dist_mat,
            positions,
            rng,
            embed_params.rand_neg_eig,
            embed_params.num_zero_fail as usize,
        )?)
    } else {
        let box_size = if embed_params.box_size_mult > 0.0 {
            5.0 * embed_params.box_size_mult
        } else {
            -embed_params.box_size_mult
        };
        let got_coords = compute_random_coords_with_rng(positions, box_size, rng);
        if let Some(map) = &embed_params.coord_map {
            for (&idx, mapped_point) in map {
                let n = positions.len();
                let point = &mut positions[coordinate_map_index(idx, n)?];
                if point.len() < 3 {
                    return Err(GenerationError::Input("point dimension less than three"));
                }
                point[..3].copy_from_slice(mapped_point);
                for coord in point.iter_mut().skip(3) {
                    *coord = 0.0;
                }
            }
        }
        Ok(got_coords)
    }
}
fn embedder_first_minimization(
    positions: &mut [Vec<f64>],
    eargs: &EmbedArgs<'_>,
    embed_params: &EmbedParams,
) -> Result<bool, GenerationError> {
    embedder_first_minimization_with_basin(
        positions,
        eargs,
        embed_params,
        embed_params.basin_thresh,
    )
}
fn embedder_first_minimization_with_basin(
    positions: &mut [Vec<f64>],
    eargs: &EmbedArgs<'_>,
    embed_params: &EmbedParams,
    basin_thresh: f64,
) -> Result<bool, GenerationError> {
    // BEGIN RDKIT CPP FUNCTION embedder_first_minimization (GraphMol/DistGeomHelpers/Embedder.cpp:518-560)
    // RDKit❗❌: bool firstMinimization(RDGeom::PointPtrVect *positions,
    // RDKit❗❌:                        const detail::EmbedArgs &eargs,
    // RDKit❗❌:                        const EmbedParameters &embedParams) {
    // RDKit❗❌:   bool gotCoords = true;
    // RDKit❗❌:   boost::dynamic_bitset<> fixedPts(positions->size());
    // RDKit❗❌:   if (embedParams.useRandomCoords && embedParams.coordMap != nullptr) {
    // RDKit❗❌:     for (const auto &v : *embedParams.coordMap) {
    // RDKit❗❌:       fixedPts.set(v.first);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   std::unique_ptr<ForceFields::ForceField> field(DistGeom::constructForceField(
    // RDKit❗❌:       *eargs.mmat, *positions, *eargs.chiralCenters, 1.0, 0.1, nullptr,
    // RDKit❗❌:       embedParams.basinThresh, &fixedPts));
    // RDKit❗❌:   if (embedParams.useRandomCoords && embedParams.coordMap != nullptr) {
    // RDKit❗❌:     for (const auto &v : *embedParams.coordMap) {
    // RDKit❗❌:       field->fixedPoints().push_back(v.first);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   field->initialize();
    // RDKit❗❌:   if (field->calcEnergy() > ERROR_TOL) {
    // RDKit❗❌:     int needMore = 1;
    // RDKit❗❌:     while (needMore) {
    // RDKit❗❌:       needMore = field->minimize(400, embedParams.optimizerForceTol);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   std::vector<double> e_contribs;
    // RDKit❗❌:   double local_e = field->calcEnergy(&e_contribs);
    // RDKit❗❌:
    // RDKit❗❌: #ifdef DEBUG_EMBEDDING
    // RDKit❗❌:   std::cerr << " Energy : " << local_e / positions->size() << " "
    // RDKit❗❌:             << *(std::max_element(e_contribs.begin(), e_contribs.end()))
    // RDKit❗❌:             << std::endl;
    // RDKit❗❌: #endif
    // RDKit❗❌:
    // RDKit❗❌:   // check that the energy is not too high (this is part of github #971)
    // RDKit❗❌:   if (local_e / positions->size() >= MAX_MINIMIZED_E_PER_ATOM) {
    // RDKit❗❌: #ifdef DEBUG_EMBEDDING
    // RDKit❗❌:     std::cerr << " Energy fail: " << local_e / positions->size() << std::endl;
    // RDKit❗❌: #endif
    // RDKit❗❌:     gotCoords = false;
    // RDKit❗❌:   }
    // RDKit❗❌:   return gotCoords;
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION embedder_first_minimization

    // Source boost dynamic bitset is packed; Vec<bool> is byte storage, a known
    // memory-cost regression retained explicitly in the complexity marker.
    let n = positions.len();
    let fixed = fixed_points(embed_params, n)?;
    let mut fixed_pts = vec![false; n];
    for &i in &fixed {
        fixed_pts[i] = true;
    }
    let mut field = ConformerOptimizer::distance_geometry(
        eargs.mmat,
        positions,
        eargs.chiral_centers,
        DistanceGeometryForceFieldParams {
            weight_chiral: 1.0,
            weight_fourth_dimension: 0.1,
            extra_weights: None,
            basin_size_tolerance: basin_thresh,
            fixed_pair_points: Some(&fixed_pts),
            fixed_points: &fixed,
        },
    )?;
    if field.energy(None)? > EMBEDDER_ERROR_TOL {
        let mut need_more = 1;
        while need_more != 0 {
            need_more = field.minimize(400, embed_params.optimizer_force_tol, 1.0e-6)?;
        }
    }
    let mut e_contribs = Vec::new();
    let local_e = field.energy(Some(&mut e_contribs))?;
    Ok(!(local_e / n as f64 >= MAX_MINIMIZED_E_PER_ATOM))
}
fn embedder_check_tetrahedral_centers(
    positions: &[Vec<f64>],
    eargs: &EmbedArgs<'_>,
    _params: &EmbedParams,
) -> Result<bool, GenerationError> {
    // BEGIN RDKIT CPP FUNCTION embedder_check_tetrahedral_centers (GraphMol/DistGeomHelpers/Embedder.cpp:562-586)
    // RDKit❗✔️: bool checkTetrahedralCenters(const RDGeom::PointPtrVect *positions,
    // RDKit❗✔️:                              const detail::EmbedArgs &eargs,
    // RDKit❗✔️:                              const EmbedParameters &) {
    // RDKit❗✔️:   // for each of the atoms in the "tetrahedralCarbons" list, make sure
    // RDKit❗✔️:   // that there is a minimum volume around them and that they are inside
    // RDKit❗✔️:   // that volume. (this is part of github #971)
    // RDKit❗✔️:   for (const auto &tetSet : *eargs.tetrahedralCarbons) {
    // RDKit❗✔️:     // it could happen that the centroid is outside the volume defined
    // RDKit❗✔️:     // by the other
    // RDKit❗✔️:     // four points. That is also a fail.
    // RDKit❗✔️:     if (!_volumeTest(tetSet, *positions) ||
    // RDKit❗✔️:         !_centerInVolume(tetSet, *positions, TETRAHEDRAL_CENTERINVOLUME_TOL)) {
    // RDKit❗✔️: #ifdef DEBUG_EMBEDDING
    // RDKit❗✔️:       std::cerr << " fail2! (" << tetSet->d_idx0 << ") iter: "  //<< iter
    // RDKit❗✔️:                 << " vol: " << _volumeTest(tetSet, *positions, true)
    // RDKit❗✔️:                 << " center: "
    // RDKit❗✔️:                 << _centerInVolume(tetSet, *positions,
    // RDKit❗✔️:                                    TETRAHEDRAL_CENTERINVOLUME_TOL, true)
    // RDKit❗✔️:                 << std::endl;
    // RDKit❗✔️: #endif
    // RDKit❗✔️:       return false;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return true;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION embedder_check_tetrahedral_centers

    for set in eargs.tetrahedral_carbons {
        if !embedder_volume_test(set, positions)?
            || !embedder_center_in_volume(set, positions, TETRAHEDRAL_CENTERINVOLUME_TOL)?
        {
            return Ok(false);
        }
    }
    Ok(true)
}
fn embedder_check_chiral_centers(
    positions: &[Vec<f64>],
    eargs: &EmbedArgs<'_>,
    _params: &EmbedParams,
) -> Result<bool, GenerationError> {
    // BEGIN RDKIT CPP FUNCTION embedder_check_chiral_centers (GraphMol/DistGeomHelpers/Embedder.cpp:587-607)
    // RDKit❗✔️: bool checkChiralCenters(const RDGeom::PointPtrVect *positions,
    // RDKit❗✔️:                         const detail::EmbedArgs &eargs,
    // RDKit❗✔️:                         const EmbedParameters &) {
    // RDKit❗✔️:   // check the chiral volume:
    // RDKit❗✔️:   for (const auto &chiralSet : *eargs.chiralCenters) {
    // RDKit❗✔️:     double vol = DistGeom::calcChiralVolume(
    // RDKit❗✔️:         chiralSet->d_idx1, chiralSet->d_idx2, chiralSet->d_idx3,
    // RDKit❗✔️:         chiralSet->d_idx4, *positions);
    // RDKit❗✔️:     double lb = chiralSet->getLowerVolumeBound();
    // RDKit❗✔️:     double ub = chiralSet->getUpperVolumeBound();
    // RDKit❗✔️:     if ((lb > 0 && vol < lb && (vol / lb < .8 || haveOppositeSign(vol, lb))) ||
    // RDKit❗✔️:         (ub < 0 && vol > ub && (vol / ub < .8 || haveOppositeSign(vol, ub)))) {
    // RDKit❗✔️: #ifdef DEBUG_EMBEDDING
    // RDKit❗✔️:       std::cerr << " fail! (" << chiralSet->d_idx0 << ") iter: "
    // RDKit❗✔️:                 << " " << vol << " " << lb << "-" << ub << std::endl;
    // RDKit❗✔️: #endif
    // RDKit❗✔️:       return false;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return true;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION embedder_check_chiral_centers

    for set in eargs.chiral_centers {
        let vol = calc_chiral_volume_rows([set.idx1, set.idx2, set.idx3, set.idx4], positions)?;
        let lb = set.get_lower_volume_bound();
        let ub = set.get_upper_volume_bound();
        if (lb > 0.0 && vol < lb && (vol / lb < 0.8 || have_opposite_sign(vol, lb)))
            || (ub < 0.0 && vol > ub && (vol / ub < 0.8 || have_opposite_sign(vol, ub)))
        {
            return Ok(false);
        }
    }
    Ok(true)
}
fn have_opposite_sign(a: f64, b: f64) -> bool {
    // BEGIN RDKIT CPP FUNCTION have_opposite_sign (GraphMol/DistGeomHelpers/Embedder.cpp:67-69)
    // RDKit❗✔️: inline bool haveOppositeSign(double a, double b) {
    // RDKit❗✔️:   return std::signbit(a) ^ std::signbit(b);
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION have_opposite_sign
    a.is_sign_negative() ^ b.is_sign_negative()
}
fn embedder_minimize_fourth_dimension(
    positions: &mut [Vec<f64>],
    eargs: &EmbedArgs<'_>,
    embed_params: &EmbedParams,
    end_time: Option<Instant>,
) -> Result<bool, GenerationError> {
    embedder_minimize_fourth_dimension_with_basin(
        positions,
        eargs,
        embed_params,
        end_time,
        embed_params.basin_thresh,
    )
}
fn embedder_minimize_fourth_dimension_with_basin(
    positions: &mut [Vec<f64>],
    eargs: &EmbedArgs<'_>,
    embed_params: &EmbedParams,
    end_time: Option<Instant>,
    basin_thresh: f64,
) -> Result<bool, GenerationError> {
    // BEGIN RDKIT CPP FUNCTION embedder_minimize_fourth_dimension (GraphMol/DistGeomHelpers/Embedder.cpp:608-638)
    // RDKit❗❌: bool minimizeFourthDimension(RDGeom::PointPtrVect *positions,
    // RDKit❗❌:                              const detail::EmbedArgs &eargs,
    // RDKit❗❌:                              EmbedParameters &embedParams,
    // RDKit❗❌:                              TimePoint *end_time) {
    // RDKit❗❌:   // now redo the minimization if we have a chiral center
    // RDKit❗❌:   // or have started from random coords. This
    // RDKit❗❌:   // time removing the chiral constraints and
    // RDKit❗❌:   // increasing the weight on the fourth dimension
    // RDKit❗❌:
    // RDKit❗❌:   std::unique_ptr<ForceFields::ForceField> field2(DistGeom::constructForceField(
    // RDKit❗❌:       *eargs.mmat, *positions, *eargs.chiralCenters, 0.2, 1.0, nullptr,
    // RDKit❗❌:       embedParams.basinThresh));
    // RDKit❗❌:   if (embedParams.useRandomCoords && embedParams.coordMap != nullptr) {
    // RDKit❗❌:     for (const auto &v : *embedParams.coordMap) {
    // RDKit❗❌:       field2->fixedPoints().push_back(v.first);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   field2->initialize();
    // RDKit❗❌:   // std::cerr << "FIELD2 E: " << field2->calcEnergy() << std::endl;
    // RDKit❗❌:   if (field2->calcEnergy() > ERROR_TOL) {
    // RDKit❗❌:     int needMore = 1;
    // RDKit❗❌:     while (needMore) {
    // RDKit❗❌:       if (end_time != nullptr && Clock::now() > *end_time) {
    // RDKit❗❌:         return false;
    // RDKit❗❌:       }
    // RDKit❗❌:       needMore = field2->minimize(200, embedParams.optimizerForceTol);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return true;
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION embedder_minimize_fourth_dimension

    // A checked detached fixed-index projection allocates a temporary Vec before
    // the existing kernel stores its indices; source uses that kernel list directly.
    let fixed = fixed_points(embed_params, positions.len())?;
    let mut field = ConformerOptimizer::distance_geometry(
        eargs.mmat,
        positions,
        eargs.chiral_centers,
        DistanceGeometryForceFieldParams {
            weight_chiral: 0.2,
            weight_fourth_dimension: 1.0,
            extra_weights: None,
            basin_size_tolerance: basin_thresh,
            fixed_pair_points: None,
            fixed_points: &fixed,
        },
    )?;
    if field.energy(None)? > EMBEDDER_ERROR_TOL {
        let mut need_more = 1;
        while need_more != 0 {
            if let Some(deadline) = end_time
                && Instant::now() > deadline
            {
                return Ok(false);
            }
            need_more = field.minimize(200, embed_params.optimizer_force_tol, 1.0e-6)?;
        }
    }
    Ok(true)
}
fn embedder_minimize_with_exp_torsions(
    positions: &mut [Vec<f64>],
    eargs: &EmbedArgs<'_>,
    embed_params: &EmbedParams,
) -> Result<bool, GenerationError> {
    // BEGIN RDKIT CPP FUNCTION embedder_minimize_with_exp_torsions (GraphMol/DistGeomHelpers/Embedder.cpp:640-718)
    // RDKit❗❌: bool minimizeWithExpTorsions(RDGeom::PointPtrVect &positions,
    // RDKit❗❌:                              const detail::EmbedArgs &eargs,
    // RDKit❗❌:                              const EmbedParameters &embedParams) {
    // RDKit❗❌:   PRECONDITION(eargs.etkdgDetails, "bogus etkdgDetails pointer");
    // RDKit❗❌:   bool planar = true;
    // RDKit❗❌:
    // RDKit❗❌:   // convert to 3D positions and create coordMap
    // RDKit❗❌:   RDGeom::Point3DPtrVect positions3D;
    // RDKit❗❌:   for (auto &position : positions) {
    // RDKit❗❌:     positions3D.push_back(
    // RDKit❗❌:         new RDGeom::Point3D((*position)[0], (*position)[1], (*position)[2]));
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // create the force field
    // RDKit❗❌:   std::unique_ptr<ForceFields::ForceField> field;
    // RDKit❗❌:   if (embedParams.useBasicKnowledge) {  // ETKDG or KDG
    // RDKit❗❌:     if (embedParams.CPCI != nullptr) {
    // RDKit❗❌:       field.reset(DistGeom::construct3DForceField(
    // RDKit❗❌:           *eargs.mmat, positions3D, *eargs.etkdgDetails, *embedParams.CPCI));
    // RDKit❗❌:     } else {
    // RDKit❗❌:       field.reset(DistGeom::construct3DForceField(*eargs.mmat, positions3D,
    // RDKit❗❌:                                                   *eargs.etkdgDetails));
    // RDKit❗❌:     }
    // RDKit❗❌:   } else {  // plain ETDG
    // RDKit❗❌:     field.reset(DistGeom::constructPlain3DForceField(*eargs.mmat, positions3D,
    // RDKit❗❌:                                                      *eargs.etkdgDetails));
    // RDKit❗❌:   }
    // RDKit❗❌:   if (embedParams.useRandomCoords && embedParams.coordMap != nullptr) {
    // RDKit❗❌:     for (const auto &v : *embedParams.coordMap) {
    // RDKit❗❌:       field->fixedPoints().push_back(v.first);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // minimize!
    // RDKit❗❌:   field->initialize();
    // RDKit❗❌:   if (field->calcEnergy() > ERROR_TOL) {
    // RDKit❗❌:     // while (needMore) {
    // RDKit❗❌:     field->minimize(300, embedParams.optimizerForceTol);
    // RDKit❗❌:     //      ++nPasses;
    // RDKit❗❌:     //}
    // RDKit❗❌:   }
    // RDKit❗❌:   // std::cout << field->calcEnergy() << std::endl;
    // RDKit❗❌:
    // RDKit❗❌:   // check for planarity if ETKDG or KDG
    // RDKit❗❌:   if (embedParams.useBasicKnowledge) {
    // RDKit❗❌:     // create a force field with only the impropers
    // RDKit❗❌:     std::unique_ptr<ForceFields::ForceField> field2(
    // RDKit❗❌:         DistGeom::construct3DImproperForceField(*eargs.mmat, positions3D,
    // RDKit❗❌:                                                 *eargs.etkdgDetails));
    // RDKit❗❌:     if (embedParams.useRandomCoords && embedParams.coordMap != nullptr) {
    // RDKit❗❌:       for (const auto &v : *embedParams.coordMap) {
    // RDKit❗❌:         field2->fixedPoints().push_back(v.first);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     field2->initialize();
    // RDKit❗❌:     // check if the energy is low enough
    // RDKit❗❌:     double planarityTolerance = 0.7;
    // RDKit❗❌:     if (field2->calcEnergy() >
    // RDKit❗❌:         eargs.etkdgDetails->improperAtoms.size() * planarityTolerance) {
    // RDKit❗❌: #ifdef DEBUG_EMBEDDING
    // RDKit❗❌:       std::cerr << "   planar fail: " << field2->calcEnergy() << " "
    // RDKit❗❌:                 << eargs.etkdgDetails->improperAtoms.size() * planarityTolerance
    // RDKit❗❌:                 << std::endl;
    // RDKit❗❌: #endif
    // RDKit❗❌:       planar = false;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // overwrite positions and delete the 3D ones
    // RDKit❗❌:   for (unsigned int i = 0; i < positions3D.size(); ++i) {
    // RDKit❗❌:     (*positions[i])[0] = (*positions3D[i])[0];
    // RDKit❗❌:     (*positions[i])[1] = (*positions3D[i])[1];
    // RDKit❗❌:     (*positions[i])[2] = (*positions3D[i])[2];
    // RDKit❗❌:     delete positions3D[i];
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   return planar;
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION embedder_minimize_with_exp_torsions

    let details = eargs
        .etkdg_details
        .ok_or(GenerationError::Input("bogus etkdgDetails pointer"))?;
    let mut positions_3d: Vec<Vec<f64>> = positions
        .iter()
        .enumerate()
        .map(|(i, _)| point_at(positions, i).map(|p| p.to_vec()))
        .collect::<Result<_, _>>()?;
    // A checked detached fixed-index projection allocates a temporary Vec before
    // the existing kernel stores its indices; source uses that kernel list directly.
    let fixed = fixed_points(embed_params, positions.len())?;
    let cpci: Option<BTreeMap<(usize, usize), f64>> = if embed_params.use_basic_knowledge {
        embed_params.cpci.as_ref().map(|m| {
            m.iter()
                .map(|(&(i, j), &q)| ((i as usize, j as usize), q))
                .collect()
        })
    } else {
        None
    };
    {
        let mut field = ConformerOptimizer::torsions(
            eargs.mmat,
            &mut positions_3d,
            details,
            embed_params.use_basic_knowledge,
            cpci.as_ref(),
            &fixed,
        )?;
        if field.energy(None)? > EMBEDDER_ERROR_TOL {
            let _need_more = field.minimize(300, embed_params.optimizer_force_tol, 1.0e-6)?;
        }
    }
    let mut planar = true;
    if embed_params.use_basic_knowledge {
        let mut field =
            ConformerOptimizer::improper(eargs.mmat, &mut positions_3d, details, &fixed)?;
        if field.energy(None)? > details.improper_atoms.len() as f64 * 0.7 {
            planar = false;
        }
    }
    for (point, point3d) in positions.iter_mut().zip(positions_3d) {
        point[..3].copy_from_slice(&point3d);
    }
    Ok(planar)
}

#[cfg(test)]
mod original_numeric_generation_conditions {
    use super::*;
    use std::{sync::Arc, time::Duration};
    #[test]
    fn embedder_volume_test_accepts_well_separated_tetrahedral_center() {
        let positions = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, -1.0, -1.0],
        ];
        let chiral_set = ChiralSet::with_default_structure_flags(0, 1, 2, 3, 4, -1.0, 1.0);

        assert!(embedder_volume_test(&chiral_set, &positions).expect("original checked phase"));
    }

    #[test]
    fn embedder_volume_test_rejects_flat_or_low_volume_center() {
        let flat_positions = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, -1.0, 0.0],
        ];
        let chiral_set = ChiralSet::with_default_structure_flags(0, 1, 2, 3, 4, -1.0, 1.0);

        assert!(
            !embedder_volume_test(&chiral_set, &flat_positions).expect("original checked phase")
        );

        let low_volume_positions = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, -1.0, -0.3],
        ];
        assert!(
            !embedder_volume_test(&chiral_set, &low_volume_positions)
                .expect("original checked phase")
        );
    }

    #[test]
    fn embedder_volume_test_uses_fused_small_ring_relaxed_threshold() {
        let positions = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, -1.0, -0.3],
        ];
        let regular = ChiralSet::with_default_structure_flags(0, 1, 2, 3, 4, -1.0, 1.0);
        let fused = ChiralSet::new(
            0,
            1,
            2,
            3,
            4,
            -1.0,
            1.0,
            ChiralSetStructureFlags::InFusedSmallRings as u64,
        );

        assert!(!embedder_volume_test(&regular, &positions).expect("original checked phase"));
        assert!(embedder_volume_test(&fused, &positions).expect("original checked phase"));
    }

    #[test]
    fn embedder_same_side_matches_plane_side_and_tolerance_rules() {
        let v1 = [0.0, 0.0, 0.0];
        let v2 = [1.0, 0.0, 0.0];
        let v3 = [0.0, 1.0, 0.0];
        let v4 = [0.0, 0.0, 1.0];

        assert!(embedder_same_side(v1, v2, v3, v4, [0.25, 0.25, 0.5], 0.1));
        assert!(!embedder_same_side(v1, v2, v3, v4, [0.25, 0.25, -0.5], 0.1));
        assert!(!embedder_same_side(v1, v2, v3, v4, [0.25, 0.25, 0.05], 0.1));
        assert!(!embedder_same_side(
            v1,
            v2,
            v3,
            [0.0, 0.0, 0.05],
            [0.25, 0.25, 0.5],
            0.1
        ));
    }

    #[test]
    fn embedder_center_in_volume_checks_all_four_faces() {
        let positions = vec![
            [0.1, 0.1, 0.1],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [2.0, 2.0, 2.0],
        ];

        assert!(
            embedder_center_in_volume_indices(0, 1, 2, 3, 4, &positions, 0.01)
                .expect("original checked phase")
        );
        assert!(
            !embedder_center_in_volume_indices(5, 1, 2, 3, 4, &positions, 0.01)
                .expect("original checked phase")
        );
    }

    #[test]
    fn embedder_center_in_volume_chiral_set_overload_handles_three_coordinate_centers() {
        let positions = vec![
            [10.0, 10.0, 10.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let three_coordinate = ChiralSet::with_default_structure_flags(0, 1, 2, 3, 0, -1.0, 1.0);

        assert!(
            embedder_center_in_volume(&three_coordinate, &positions, 0.1)
                .expect("original checked phase")
        );
    }

    #[test]
    fn embedder_center_in_volume_chiral_set_overload_delegates_indices() {
        let positions = vec![
            [0.1, 0.1, 0.1],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let chiral_set = ChiralSet::with_default_structure_flags(0, 1, 2, 3, 4, -1.0, 1.0);

        assert!(
            embedder_center_in_volume(&chiral_set, &positions, 0.01)
                .expect("original checked phase")
        );
    }

    #[test]
    fn embedder_bounds_fulfilled_matches_rdkit_tolerance_rule() {
        let mut mmat = BoundsMatrix::new(3).expect("original bounds dimension");
        mmat.set_lower(0, 1, 0.9).expect("set lower");
        mmat.set_upper(0, 1, 1.1).expect("set upper");
        mmat.set_lower(0, 2, 1.0).expect("set lower");
        mmat.set_upper(0, 2, 2.0).expect("set upper");
        mmat.set_lower(1, 2, 1.0).expect("set lower");
        mmat.set_upper(1, 2, 2.0).expect("set upper");
        let atoms = [0, 1, 2];
        let ok_positions = vec![[0.0, 0.0, 0.0], [1.15, 0.0, 0.0], [0.0, 1.2, 0.0]];
        let bad_positions = vec![[0.0, 0.0, 0.0], [1.25, 0.0, 0.0], [0.0, 1.2, 0.0]];

        assert!(
            embedder_bounds_fulfilled(&atoms, &mmat, &ok_positions)
                .expect("original checked phase")
        );
        assert!(
            !embedder_bounds_fulfilled(&atoms, &mmat, &bad_positions)
                .expect("original checked phase")
        );
        assert!(
            embedder_bounds_fulfilled(&[], &mmat, &bad_positions).expect("original checked phase")
        );
        assert!(
            embedder_bounds_fulfilled(&[0], &mmat, &bad_positions).expect("original checked phase")
        );
    }

    #[test]
    fn embedder_generate_initial_coords_random_branch_applies_coord_map_and_zeroes_higher_dimensions()
     {
        let mut mmat = BoundsMatrix::new(3).expect("original bounds dimension");
        mmat.set_lower(0, 1, 1.0).expect("set lower");
        mmat.set_upper(0, 1, 1.0).expect("set upper");
        let chiral_centers: Vec<ChiralSetPtr> = Vec::new();
        let tetrahedral_carbons: Vec<ChiralSetPtr> = Vec::new();
        let eargs = EmbedArgs {
            mmat: &mmat,
            chiral_centers: &chiral_centers,
            tetrahedral_carbons: &tetrahedral_carbons,
            etkdg_details: None,
            double_bond_ends: None,
            stereo_double_bonds: &[],
        };
        let mut coord_map = BTreeMap::new();
        coord_map.insert(1, [7.0, 8.0, 9.0]);
        let mut params = EmbedParams::default();
        params.use_random_coords = true;
        params.box_size_mult = 2.0;
        params.coord_map = Some(coord_map);
        let mut positions = vec![vec![0.0; 4], vec![0.0; 4], vec![0.0; 4]];
        let mut expected_positions = positions.clone();
        let mut expected_rng = cosmolkit_core::RdkitRandomEngine::from_seed(123);
        assert!(compute_random_coords_with_rng(
            &mut expected_positions,
            10.0,
            &mut expected_rng
        ));
        expected_positions[1][0] = 7.0;
        expected_positions[1][1] = 8.0;
        expected_positions[1][2] = 9.0;
        expected_positions[1][3] = 0.0;

        let mut dist_mat = SymmMatrix::new(3);
        let mut rng = cosmolkit_core::RdkitRandomEngine::from_seed(123);

        assert!(
            embedder_generate_initial_coords(
                &mut positions,
                &eargs,
                &params,
                &mut dist_mat,
                &mut rng
            )
            .expect("generate initial coords")
        );
        assert_eq!(positions, expected_positions);
    }

    #[test]
    fn embedder_generate_initial_coords_distance_matrix_branch_uses_bounds_matrix() {
        let mut mmat = BoundsMatrix::new(2).expect("original bounds dimension");
        mmat.set_lower(0, 1, 2.0).expect("set lower");
        mmat.set_upper(0, 1, 2.0).expect("set upper");
        let chiral_centers: Vec<ChiralSetPtr> = Vec::new();
        let tetrahedral_carbons: Vec<ChiralSetPtr> = Vec::new();
        let eargs = EmbedArgs {
            mmat: &mmat,
            chiral_centers: &chiral_centers,
            tetrahedral_carbons: &tetrahedral_carbons,
            etkdg_details: None,
            double_bond_ends: None,
            stereo_double_bonds: &[],
        };
        let params = EmbedParams::default();
        let mut positions = vec![vec![0.0; 3], vec![0.0; 3]];
        let mut dist_mat = SymmMatrix::new(2);
        let mut rng = cosmolkit_core::RdkitRandomEngine::from_seed(222);

        assert!(
            embedder_generate_initial_coords(
                &mut positions,
                &eargs,
                &params,
                &mut dist_mat,
                &mut rng
            )
            .expect("generate initial coords")
        );
        assert_eq!(dist_mat.get_val(0, 1), 2.0);
        let distance = ((positions[0][0] - positions[1][0]).powi(2)
            + (positions[0][1] - positions[1][1]).powi(2)
            + (positions[0][2] - positions[1][2]).powi(2))
        .sqrt();
        assert!((distance - 2.0).abs() < 1.0e-6);
    }

    #[test]
    fn embedder_first_minimization_keeps_exact_satisfied_two_point_bounds() {
        let mut mmat = BoundsMatrix::new(2).expect("original bounds dimension");
        mmat.set_lower(0, 1, 1.0).expect("set lower");
        mmat.set_upper(0, 1, 1.0).expect("set upper");
        let chiral_centers: Vec<ChiralSetPtr> = Vec::new();
        let tetrahedral_carbons: Vec<ChiralSetPtr> = Vec::new();
        let eargs = EmbedArgs {
            mmat: &mmat,
            chiral_centers: &chiral_centers,
            tetrahedral_carbons: &tetrahedral_carbons,
            etkdg_details: None,
            double_bond_ends: None,
            stereo_double_bonds: &[],
        };
        let params = EmbedParams::default();
        let mut positions = vec![vec![0.0, 0.0, 0.0], vec![1.0, 0.0, 0.0]];

        assert!(
            embedder_first_minimization(&mut positions, &eargs, &params)
                .expect("original checked phase")
        );
        assert_eq!(positions, vec![vec![0.0, 0.0, 0.0], vec![1.0, 0.0, 0.0]]);
    }

    #[test]
    fn embedder_check_tetrahedral_centers_requires_volume_and_center_in_volume() {
        let mmat = BoundsMatrix::new(5).expect("original bounds dimension");
        let tet_set: ChiralSetPtr = Arc::new(ChiralSet::with_default_structure_flags(
            0, 1, 2, 3, 4, -1.0, 1.0,
        ));
        let chiral_centers: Vec<ChiralSetPtr> = Vec::new();
        let tetrahedral_carbons: Vec<ChiralSetPtr> = vec![Arc::clone(&tet_set)];
        let eargs = EmbedArgs {
            mmat: &mmat,
            chiral_centers: &chiral_centers,
            tetrahedral_carbons: &tetrahedral_carbons,
            etkdg_details: None,
            double_bond_ends: None,
            stereo_double_bonds: &[],
        };
        let params = EmbedParams::default();
        let good_positions = vec![
            vec![0.0, 0.0, 0.0],
            vec![1.0, 1.0, 1.0],
            vec![1.0, -1.0, -1.0],
            vec![-1.0, 1.0, -1.0],
            vec![-1.0, -1.0, 1.0],
        ];
        let outside_positions = vec![
            vec![4.0, 4.0, 4.0],
            vec![1.0, 1.0, 1.0],
            vec![1.0, -1.0, -1.0],
            vec![-1.0, 1.0, -1.0],
            vec![-1.0, -1.0, 1.0],
        ];

        assert!(
            embedder_check_tetrahedral_centers(&good_positions, &eargs, &params)
                .expect("original checked phase")
        );
        assert!(
            !embedder_check_tetrahedral_centers(&outside_positions, &eargs, &params)
                .expect("original checked phase")
        );
    }

    #[test]
    fn embedder_check_chiral_centers_matches_rdkit_volume_bound_rule() {
        let mmat = BoundsMatrix::new(5).expect("original bounds dimension");
        let good_set: ChiralSetPtr = Arc::new(ChiralSet::with_default_structure_flags(
            0, 1, 2, 3, 4, 0.5, 2.0,
        ));
        let failing_set: ChiralSetPtr = Arc::new(ChiralSet::with_default_structure_flags(
            0, 1, 2, 3, 4, 2.0, 3.0,
        ));
        let tetrahedral_carbons: Vec<ChiralSetPtr> = Vec::new();
        let positions = vec![
            vec![0.0, 0.0, 0.0],
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
            vec![0.0, 0.0, 0.0],
        ];
        let params = EmbedParams::default();
        let good_chiral_centers: Vec<ChiralSetPtr> = vec![good_set];
        let good_args = EmbedArgs {
            mmat: &mmat,
            chiral_centers: &good_chiral_centers,
            tetrahedral_carbons: &tetrahedral_carbons,
            etkdg_details: None,
            double_bond_ends: None,
            stereo_double_bonds: &[],
        };
        let failing_chiral_centers: Vec<ChiralSetPtr> = vec![failing_set];
        let failing_args = EmbedArgs {
            mmat: &mmat,
            chiral_centers: &failing_chiral_centers,
            tetrahedral_carbons: &tetrahedral_carbons,
            etkdg_details: None,
            double_bond_ends: None,
            stereo_double_bonds: &[],
        };

        assert!(
            embedder_check_chiral_centers(&positions, &good_args, &params)
                .expect("original checked phase")
        );
        assert!(
            !embedder_check_chiral_centers(&positions, &failing_args, &params)
                .expect("original checked phase")
        );
    }

    #[test]
    fn embedder_minimize_fourth_dimension_keeps_exact_satisfied_two_point_bounds() {
        let mut mmat = BoundsMatrix::new(2).expect("original bounds dimension");
        mmat.set_lower(0, 1, 1.0).expect("set lower");
        mmat.set_upper(0, 1, 1.0).expect("set upper");
        let chiral_centers: Vec<ChiralSetPtr> = Vec::new();
        let tetrahedral_carbons: Vec<ChiralSetPtr> = Vec::new();
        let eargs = EmbedArgs {
            mmat: &mmat,
            chiral_centers: &chiral_centers,
            tetrahedral_carbons: &tetrahedral_carbons,
            etkdg_details: None,
            double_bond_ends: None,
            stereo_double_bonds: &[],
        };
        let mut params = EmbedParams::default();
        let mut positions = vec![vec![0.0, 0.0, 0.0, 0.0], vec![1.0, 0.0, 0.0, 0.0]];

        assert!(
            embedder_minimize_fourth_dimension(&mut positions, &eargs, &mut params, None)
                .expect("original checked phase")
        );
        assert_eq!(
            positions,
            vec![vec![0.0, 0.0, 0.0, 0.0], vec![1.0, 0.0, 0.0, 0.0]]
        );
    }

    #[test]
    fn embedder_minimize_fourth_dimension_preserves_random_coord_map_fixed_points() {
        let mut mmat = BoundsMatrix::new(2).expect("original bounds dimension");
        mmat.set_lower(0, 1, 1.0).expect("set lower");
        mmat.set_upper(0, 1, 1.0).expect("set upper");
        let chiral_centers: Vec<ChiralSetPtr> = Vec::new();
        let tetrahedral_carbons: Vec<ChiralSetPtr> = Vec::new();
        let eargs = EmbedArgs {
            mmat: &mmat,
            chiral_centers: &chiral_centers,
            tetrahedral_carbons: &tetrahedral_carbons,
            etkdg_details: None,
            double_bond_ends: None,
            stereo_double_bonds: &[],
        };
        let mut coord_map = BTreeMap::new();
        coord_map.insert(0, [7.0, 8.0, 9.0]);
        let mut params = EmbedParams::default();
        params.use_random_coords = true;
        params.coord_map = Some(coord_map);
        let mut positions = vec![vec![7.0, 8.0, 9.0, 0.0], vec![10.0, 8.0, 9.0, 3.0]];

        assert!(
            embedder_minimize_fourth_dimension(&mut positions, &eargs, &mut params, None)
                .expect("original checked phase")
        );
        assert_eq!(positions[0], vec![7.0, 8.0, 9.0, 0.0]);
    }

    #[test]
    fn embedder_minimize_fourth_dimension_returns_false_after_timeout_deadline() {
        let mut mmat = BoundsMatrix::new(2).expect("original bounds dimension");
        mmat.set_lower(0, 1, 1.0).expect("set lower");
        mmat.set_upper(0, 1, 1.0).expect("set upper");
        let chiral_centers: Vec<ChiralSetPtr> = Vec::new();
        let tetrahedral_carbons: Vec<ChiralSetPtr> = Vec::new();
        let eargs = EmbedArgs {
            mmat: &mmat,
            chiral_centers: &chiral_centers,
            tetrahedral_carbons: &tetrahedral_carbons,
            etkdg_details: None,
            double_bond_ends: None,
            stereo_double_bonds: &[],
        };
        let mut params = EmbedParams::default();
        let mut positions = vec![vec![0.0, 0.0, 0.0, 10.0], vec![3.0, 0.0, 0.0, -10.0]];

        assert!(
            !embedder_minimize_fourth_dimension(
                &mut positions,
                &eargs,
                &mut params,
                Some(Instant::now() - Duration::from_secs(1))
            )
            .expect("original checked phase")
        );
    }

    #[test]
    fn embedder_minimize_with_exp_torsions_plain_etdg_preserves_random_coord_map_fixed_points() {
        let mut mmat = BoundsMatrix::new(2).expect("original bounds dimension");
        mmat.set_lower(0, 1, 1.0).expect("set lower");
        mmat.set_upper(0, 1, 1.0).expect("set upper");
        let chiral_centers: Vec<ChiralSetPtr> = Vec::new();
        let tetrahedral_carbons: Vec<ChiralSetPtr> = Vec::new();
        let details = CrystalFFDetails::default();
        let eargs = EmbedArgs {
            mmat: &mmat,
            chiral_centers: &chiral_centers,
            tetrahedral_carbons: &tetrahedral_carbons,
            etkdg_details: Some(&details),
            double_bond_ends: None,
            stereo_double_bonds: &[],
        };
        let mut coord_map = BTreeMap::new();
        coord_map.insert(0, [3.0, 4.0, 5.0]);
        let mut params = EmbedParams::default();
        params.use_random_coords = true;
        params.coord_map = Some(coord_map);
        let mut positions = vec![vec![3.0, 4.0, 5.0, 9.0], vec![8.0, 4.0, 5.0, -2.0]];

        assert!(
            embedder_minimize_with_exp_torsions(&mut positions, &eargs, &params)
                .expect("original checked bogus etkdgDetails pointer")
        );
        assert_eq!(&positions[0][..3], &[3.0, 4.0, 5.0]);
        assert_eq!(positions[0][3], 9.0);
        assert_eq!(positions[1][3], -2.0);
    }

    #[test]
    fn embedder_minimize_with_exp_torsions_basic_knowledge_accepts_empty_cpci() {
        let mut mmat = BoundsMatrix::new(2).expect("original bounds dimension");
        mmat.set_lower(0, 1, 1.0).expect("set lower");
        mmat.set_upper(0, 1, 1.0).expect("set upper");
        let chiral_centers: Vec<ChiralSetPtr> = Vec::new();
        let tetrahedral_carbons: Vec<ChiralSetPtr> = Vec::new();
        let details = CrystalFFDetails::default();
        let eargs = EmbedArgs {
            mmat: &mmat,
            chiral_centers: &chiral_centers,
            tetrahedral_carbons: &tetrahedral_carbons,
            etkdg_details: Some(&details),
            double_bond_ends: None,
            stereo_double_bonds: &[],
        };
        let mut params = EmbedParams::default();
        params.use_basic_knowledge = true;
        params.cpci = Some(BTreeMap::new());
        let mut positions = vec![vec![0.0, 0.0, 0.0], vec![1.0, 0.0, 0.0]];

        assert!(
            embedder_minimize_with_exp_torsions(&mut positions, &eargs, &params)
                .expect("original checked bogus etkdgDetails pointer")
        );
    }

    #[test]
    #[should_panic(expected = "bogus etkdgDetails pointer")]
    fn embedder_minimize_with_exp_torsions_requires_etkdg_details() {
        let mmat = BoundsMatrix::new(1).expect("original bounds dimension");
        let chiral_centers: Vec<ChiralSetPtr> = Vec::new();
        let tetrahedral_carbons: Vec<ChiralSetPtr> = Vec::new();
        let eargs = EmbedArgs {
            mmat: &mmat,
            chiral_centers: &chiral_centers,
            tetrahedral_carbons: &tetrahedral_carbons,
            etkdg_details: None,
            double_bond_ends: None,
            stereo_double_bonds: &[],
        };
        let params = EmbedParams::default();
        let mut positions = vec![vec![0.0, 0.0, 0.0]];

        let _ = embedder_minimize_with_exp_torsions(&mut positions, &eargs, &params)
            .expect("original checked bogus etkdgDetails pointer");
    }
}

static FAILURE_MUTEX: std::sync::Mutex<()> = std::sync::Mutex::new(());
static SOURCE_INTERRUPTED: std::sync::atomic::AtomicBool =
    std::sync::atomic::AtomicBool::new(false);
fn source_got_signal() -> bool {
    // BEGIN RDKIT CPP FUNCTION source_got_signal (RDGeneral/ControlCHandler.h)
    // RDKit❗✔️:   static bool getGotSignal() { return d_gotSignal; }
    // END RDKIT CPP FUNCTION source_got_signal

    SOURCE_INTERRUPTED.load(std::sync::atomic::Ordering::SeqCst)
}
#[cfg(not(target_arch = "wasm32"))]
extern "C" fn source_signal_handler(signal_number: libc::c_int) {
    // BEGIN RDKIT CPP FUNCTION source_signal_handler (RDGeneral/ControlCHandler.h)
    // RDKit❗✔️:   static void signalHandler(int signalNumber) {
    // RDKit❗✔️:     if (signalNumber == SIGINT) {
    // RDKit❗✔️:       d_gotSignal = true;
    // RDKit❗✔️:       std::signal(SIGINT, d_prev_handler);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION source_signal_handler

    if signal_number == libc::SIGINT {
        SOURCE_INTERRUPTED.store(true, std::sync::atomic::Ordering::SeqCst);
        // Safety: this is the C signal() entrypoint with SIGINT and SIG_DFL, the
        // zero-initialized prior handler in the source reset-only embedding path.
        unsafe {
            libc::signal(libc::SIGINT, libc::SIG_DFL);
        }
    }
}
#[cfg(not(target_arch = "wasm32"))]
fn source_reset_interrupt() -> Result<(), GenerationError> {
    // BEGIN RDKIT CPP FUNCTION source_reset_interrupt (RDGeneral/ControlCHandler.h)
    // RDKit❗✔️:   static void reset() {
    // RDKit❗✔️:     d_gotSignal = false;
    // RDKit❗✔️:     std::signal(SIGINT, signalHandler);
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION source_reset_interrupt

    SOURCE_INTERRUPTED.store(false, std::sync::atomic::Ordering::SeqCst);
    // Safety: the callback has the C signal handler ABI, static lifetime and
    // performs only lock-free atomic storage and the source C signal reset.
    let previous = unsafe {
        libc::signal(
            libc::SIGINT,
            source_signal_handler as *const () as libc::sighandler_t,
        )
    };
    // Windows declares SIG_ERR as c_int, while signal() returns sighandler_t.
    // Convert the -1 sentinel to the return type without narrowing its bits.
    if previous == libc::SIG_ERR as libc::sighandler_t {
        return Err(GenerationError::Input(
            "cannot install embedding SIGINT handler",
        ));
    }
    Ok(())
}
#[cfg(target_arch = "wasm32")]
fn source_reset_interrupt() -> Result<(), GenerationError> {
    // BEGIN RDKIT CPP FUNCTION source_reset_interrupt (RDGeneral/ControlCHandler.h)
    // RDKit❗✔️:   static void reset() {
    // RDKit❗✔️:     d_gotSignal = false;
    // RDKit❗✔️:     std::signal(SIGINT, signalHandler);
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION source_reset_interrupt
    // Approved Web adaptation: browser/Node WASM has no process SIGINT handler.
    // Reset the same interruption state without installing an OS handler;
    // computation and deadline checks remain enabled. Host cancellation (e.g.
    // terminating a Worker) is separate from synchronous chemistry execution.
    SOURCE_INTERRUPTED.store(false, std::sync::atomic::Ordering::SeqCst);
    Ok(())
}
fn embedder_double_bond_geometry_checks(
    positions: &[Vec<f64>],
    eargs: &EmbedArgs<'_>,
    _params: &EmbedParams,
    linear_tol: f64,
) -> Result<bool, GenerationError> {
    // BEGIN RDKIT CPP FUNCTION embedder_double_bond_geometry_checks (GraphMol/DistGeomHelpers/Embedder.cpp)
    // RDKit❗✔️: bool doubleBondGeometryChecks(const RDGeom::PointPtrVect &positions,
    // RDKit❗✔️:                               const detail::EmbedArgs &eargs, EmbedParameters &,
    // RDKit❗✔️:                               double linearTol = 1e-3) {
    // RDKit❗✔️:   if (eargs.doubleBondEnds) {
    // RDKit❗✔️:     for (const auto &itm : *eargs.doubleBondEnds) {
    // RDKit❗✔️:       const auto &a0 = *positions[std::get<0>(itm)];
    // RDKit❗✔️:       const auto &a1 = *positions[std::get<1>(itm)];
    // RDKit❗✔️:       const auto &a2 = *positions[std::get<2>(itm)];
    // RDKit❗✔️:       RDGeom::Point3D p0(a0[0], a0[1], a0[2]);
    // RDKit❗✔️:       RDGeom::Point3D p1(a1[0], a1[1], a1[2]);
    // RDKit❗✔️:       RDGeom::Point3D p2(a2[0], a2[1], a2[2]);
    // RDKit❗✔️:
    // RDKit❗✔️:       // check for a linear arrangement
    // RDKit❗✔️:
    // RDKit❗✔️:       auto v1 = p1 - p0;
    // RDKit❗✔️:       v1.normalize();
    // RDKit❗✔️:       auto v2 = p1 - p2;
    // RDKit❗✔️:       v2.normalize();
    // RDKit❗✔️:       // this is the arrangement:
    // RDKit❗✔️:       //     a0
    // RDKit❗✔️:       //       \       [intentionally left blank]
    // RDKit❗✔️:       //        a1 = a2
    // RDKit❗✔️:       // we want to be sure it's not actually:
    // RDKit❗✔️:       //   ao - a1 = a2
    // RDKit❗✔️:       if (v1.dotProduct(v2) + 1.0 < linearTol) {
    // RDKit❗✔️:         return false;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return true;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION embedder_double_bond_geometry_checks

    if let Some(ends) = eargs.double_bond_ends {
        for &(i, j, k) in ends {
            let p0 = point_at(positions, i)?;
            let p1 = point_at(positions, j)?;
            let p2 = point_at(positions, k)?;
            let v1 = normalized_point_delta(p1, p0)?;
            let v2 = normalized_point_delta(p1, p2)?;
            if dot(v1, v2) + 1.0 < linear_tol {
                return Ok(false);
            }
        }
    }
    Ok(true)
}
fn embedder_double_bond_stereo_checks(
    positions: &[Vec<f64>],
    eargs: &EmbedArgs<'_>,
    _params: &EmbedParams,
) -> Result<bool, GenerationError> {
    // BEGIN RDKIT CPP FUNCTION embedder_double_bond_stereo_checks (GraphMol/DistGeomHelpers/Embedder.cpp)
    // RDKit❗✔️: bool doubleBondStereoChecks(const RDGeom::PointPtrVect &positions,
    // RDKit❗✔️:                             const detail::EmbedArgs &eargs, EmbedParameters &) {
    // RDKit❗✔️:   for (const auto &itm : *eargs.stereoDoubleBonds) {
    // RDKit❗✔️:     // itm is a pair with [controlling_atoms], sign
    // RDKit❗✔️:     // where the sign tells us about cis/trans
    // RDKit❗✔️:
    // RDKit❗✔️:     const auto &a0 = *positions[itm.first[0]];
    // RDKit❗✔️:     const auto &a1 = *positions[itm.first[1]];
    // RDKit❗✔️:     const auto &a2 = *positions[itm.first[2]];
    // RDKit❗✔️:     const auto &a3 = *positions[itm.first[3]];
    // RDKit❗✔️:     RDGeom::Point3D p0(a0[0], a0[1], a0[2]);
    // RDKit❗✔️:     RDGeom::Point3D p1(a1[0], a1[1], a1[2]);
    // RDKit❗✔️:     RDGeom::Point3D p2(a2[0], a2[1], a2[2]);
    // RDKit❗✔️:     RDGeom::Point3D p3(a3[0], a3[1], a3[2]);
    // RDKit❗✔️:
    // RDKit❗✔️:     // check the dihedral and be super permissive. Here's the logic of the
    // RDKit❗✔️:     // check:
    // RDKit❗✔️:     // The second element of the dihedralBond item contains 1 for trans
    // RDKit❗✔️:     //   bonds and -1 for cis bonds.
    // RDKit❗✔️:     // The dihedral is between 0 and 180. subtracting 90 from that gives:
    // RDKit❗✔️:     //   positive values for dihedrals > 90 (closer to trans than cis)
    // RDKit❗✔️:     //   negative values for dihedrals < 90 (closer to cis than trans)
    // RDKit❗✔️:     // So multiplying the result of the subtracion from the second element of
    // RDKit❗✔️:     //   the dihedralBond element will give a positive value if the dihedral is
    // RDKit❗✔️:     //   closer to correct than it is to incorrect and a negative value
    // RDKit❗✔️:     //   otherwise.
    // RDKit❗✔️:     auto dihedral = RDGeom::computeDihedralAngle(p0, p1, p2, p3);
    // RDKit❗✔️:     if ((dihedral - M_PI_2) * itm.second < 0) {
    // RDKit❗✔️:       // closer to incorrect than correct... it's a bad geometry
    // RDKit❗✔️:       return false;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return true;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION embedder_double_bond_stereo_checks

    for (atoms, sign) in eargs.stereo_double_bonds {
        let atoms: [usize; 4] = atoms
            .as_slice()
            .try_into()
            .map_err(|_| GenerationError::Input("wrong controlling double-bond atom count"))?;
        let points = [
            point_at(positions, atoms[0])?,
            point_at(positions, atoms[1])?,
            point_at(positions, atoms[2])?,
            point_at(positions, atoms[3])?,
        ];
        let dihedral = cosmolkit_core::unsigned_dihedral_radians(points);
        if (dihedral - std::f64::consts::FRAC_PI_2) * f64::from(*sign) < 0.0 {
            return Ok(false);
        }
    }
    Ok(true)
}
fn embedder_increment_failure(
    params: &EmbedParams,
    progress: &EmbeddingProgress,
    cause: crate::EmbedFailureCause,
) -> Result<(), GenerationError> {
    // BEGIN RDKIT CPP FUNCTION embedder_increment_failure (GraphMol/DistGeomHelpers/Embedder.cpp source conditional failure increment)
    // RDKit❗✔️:     if (embedParams.trackFailures) {
    // RDKit❗✔️: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit❗✔️:       std::lock_guard<std::mutex> lock(GetFailMutex());
    // RDKit❗✔️: #endif
    // RDKit❗✔️:       embedParams.failures[EmbedFailureCauses::CHECK_CHIRAL_CENTERS2]++;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     return false;
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION embedder_increment_failure

    if params.track_failures {
        let _guard = FAILURE_MUTEX
            .lock()
            .map_err(|_| GenerationError::Input("embedding failure mutex poisoned"))?;
        let mut failures = progress
            .failures
            .lock()
            .map_err(|_| GenerationError::Input("embedding failure progress mutex poisoned"))?;
        let count = failures
            .get_mut(cause as usize)
            .ok_or(GenerationError::Input(
                "embedding failure vector is not initialized",
            ))?;
        *count = count.wrapping_add(1);
    }
    Ok(())
}
fn embedder_final_chiral_checks(
    positions: &mut [Vec<f64>],
    eargs: &EmbedArgs<'_>,
    params: &mut EmbedParams,
) -> Result<bool, GenerationError> {
    with_embedding_progress(params, |settings, progress| {
        embedder_final_chiral_checks_shared(positions, eargs, settings, progress)
    })
}
fn embedder_final_chiral_checks_shared(
    positions: &mut [Vec<f64>],
    eargs: &EmbedArgs<'_>,
    params: &EmbedParams,
    progress: &EmbeddingProgress,
) -> Result<bool, GenerationError> {
    // BEGIN RDKIT CPP FUNCTION embedder_final_chiral_checks (GraphMol/DistGeomHelpers/Embedder.cpp)
    // RDKit❗✔️: bool finalChiralChecks(RDGeom::PointPtrVect *positions,
    // RDKit❗✔️:                        const detail::EmbedArgs &eargs,
    // RDKit❗✔️:                        EmbedParameters &embedParams) {
    // RDKit❗✔️:   // confirm chiral volumes
    // RDKit❗✔️:   if (!checkChiralCenters(positions, eargs, embedParams)) {
    // RDKit❗✔️:     if (embedParams.trackFailures) {
    // RDKit❗✔️: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit❗✔️:       std::lock_guard<std::mutex> lock(GetFailMutex());
    // RDKit❗✔️: #endif
    // RDKit❗✔️:       embedParams.failures[EmbedFailureCauses::CHECK_CHIRAL_CENTERS2]++;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     return false;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // "distance matrix" chirality test
    // RDKit❗✔️:   std::set<int> atoms;
    // RDKit❗✔️:   for (const auto &chiralSet : *eargs.chiralCenters) {
    // RDKit❗✔️:     if (chiralSet->d_idx0 != chiralSet->d_idx4) {
    // RDKit❗✔️:       atoms.insert(chiralSet->d_idx0);
    // RDKit❗✔️:       atoms.insert(chiralSet->d_idx1);
    // RDKit❗✔️:       atoms.insert(chiralSet->d_idx2);
    // RDKit❗✔️:       atoms.insert(chiralSet->d_idx3);
    // RDKit❗✔️:       atoms.insert(chiralSet->d_idx4);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   std::vector<int> atomsToCheck(atoms.begin(), atoms.end());
    // RDKit❗✔️:   if (atomsToCheck.size() > 0) {
    // RDKit❗✔️:     if (!_boundsFulfilled(atomsToCheck, *eargs.mmat, *positions)) {
    // RDKit❗✔️: #ifdef DEBUG_EMBEDDING
    // RDKit❗✔️:       std::cerr << " fail3a! (" << atomsToCheck[0] << ") iter: "  //<< iter
    // RDKit❗✔️:                 << std::endl;
    // RDKit❗✔️: #endif
    // RDKit❗✔️:       if (embedParams.trackFailures) {
    // RDKit❗✔️: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit❗✔️:         std::lock_guard<std::mutex> lock(GetFailMutex());
    // RDKit❗✔️: #endif
    // RDKit❗✔️:         embedParams.failures[EmbedFailureCauses::FINAL_CHIRAL_BOUNDS]++;
    // RDKit❗✔️:       }
    // RDKit❗✔️:
    // RDKit❗✔️:       return false;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // "center in volume" chirality test
    // RDKit❗✔️:   for (const auto &chiralSet : *eargs.chiralCenters) {
    // RDKit❗✔️:     // it could happen that the centroid is outside the volume defined
    // RDKit❗✔️:     // by the other four points. That is also a fail.
    // RDKit❗✔️:     if (!_centerInVolume(chiralSet, *positions)) {
    // RDKit❗✔️: #ifdef DEBUG_EMBEDDING
    // RDKit❗✔️:       std::cerr << " fail3b! (" << chiralSet->d_idx0 << ") iter: "  //<< iter
    // RDKit❗✔️:                 << std::endl;
    // RDKit❗✔️: #endif
    // RDKit❗✔️:       if (embedParams.trackFailures) {
    // RDKit❗✔️: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit❗✔️:         std::lock_guard<std::mutex> lock(GetFailMutex());
    // RDKit❗✔️: #endif
    // RDKit❗✔️:         embedParams.failures[EmbedFailureCauses::FINAL_CENTER_IN_VOLUME]++;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       return false;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   // FIX: do we need some kind of sanity check here for the non-atomic
    // RDKit❗✔️:   // situations (e.g. atropisomers)?
    // RDKit❗✔️:
    // RDKit❗✔️:   return true;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION embedder_final_chiral_checks

    use crate::EmbedFailureCause;
    if !embedder_check_chiral_centers(positions, eargs, params)? {
        embedder_increment_failure(params, progress, EmbedFailureCause::CheckChiralCenters2)?;
        return Ok(false);
    }
    let mut atoms = std::collections::BTreeSet::new();
    for set in eargs.chiral_centers {
        if set.idx0 != set.idx4 {
            for i in [set.idx0, set.idx1, set.idx2, set.idx3, set.idx4] {
                atoms.insert(i32::try_from(i).map_err(|_| {
                    GenerationError::Input("chiral atom index exceeds source signed range")
                })?);
            }
        }
    }
    let atoms_to_check: Vec<i32> = atoms.into_iter().collect();
    // The exact source nonempty caller precondition is retained; private empty-list
    // native extension is excluded from parity by the explicit ROOT decision.
    if !atoms_to_check.is_empty()
        && !embedder_bounds_fulfilled(&atoms_to_check, eargs.mmat, positions)?
    {
        embedder_increment_failure(params, progress, EmbedFailureCause::FinalChiralBounds)?;
        return Ok(false);
    }
    for set in eargs.chiral_centers {
        if !embedder_center_in_volume(set, positions, 0.1)? {
            embedder_increment_failure(params, progress, EmbedFailureCause::FinalCenterInVolume)?;
            return Ok(false);
        }
    }
    Ok(true)
}
fn embedder_embed_points_with_rng<R: RdkitDoubleRng>(
    positions: &mut [Vec<f64>],
    eargs: &EmbedArgs<'_>,
    params: &EmbedParams,
    progress: &EmbeddingProgress,
    max_iterations: u32,
    basin_thresh: f64,
    end_time: Option<Instant>,
    dist_mat: &mut SymmMatrix,
    rng: &mut R,
) -> Result<bool, GenerationError> {
    // BEGIN RDKIT CPP FUNCTION embedder_embed_points_with_rng (GraphMol/DistGeomHelpers/Embedder.cpp complete source retry loop)
    // RDKit❗❌:   bool gotCoords = false;
    // RDKit❗❌:   unsigned int iter = 0;
    // RDKit❗❌:   while (!gotCoords && iter < embedParams.maxIterations) {
    // RDKit❗❌:     if (end_time != nullptr && Clock::now() > *end_time) {
    // RDKit❗❌:       break;
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     ++iter;
    // RDKit❗❌:     if (embedParams.callback != nullptr) {
    // RDKit❗❌:       embedParams.callback(iter);
    // RDKit❗❌:     }
    // RDKit❗❌:     if (ControlCHandler::getGotSignal()) {
    // RDKit❗❌:       return false;
    // RDKit❗❌:     }
    // RDKit❗❌:     gotCoords = EmbeddingOps::generateInitialCoords(positions, eargs,
    // RDKit❗❌:                                                     embedParams, distMat, rng);
    // RDKit❗❌:     if (!gotCoords) {
    // RDKit❗❌:       if (embedParams.trackFailures) {
    // RDKit❗❌: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit❗❌:         std::lock_guard<std::mutex> lock(GetFailMutex());
    // RDKit❗❌: #endif
    // RDKit❗❌:         embedParams.failures[EmbedFailureCauses::INITIAL_COORDS]++;
    // RDKit❗❌:       }
    // RDKit❗❌:     } else {
    // RDKit❗❌:       if (ControlCHandler::getGotSignal()) {
    // RDKit❗❌:         return false;
    // RDKit❗❌:       }
    // RDKit❗❌:       gotCoords =
    // RDKit❗❌:           EmbeddingOps::firstMinimization(positions, eargs, embedParams);
    // RDKit❗❌:       if (!gotCoords) {
    // RDKit❗❌:         if (embedParams.trackFailures) {
    // RDKit❗❌: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit❗❌:           std::lock_guard<std::mutex> lock(GetFailMutex());
    // RDKit❗❌: #endif
    // RDKit❗❌:           embedParams.failures[EmbedFailureCauses::FIRST_MINIMIZATION]++;
    // RDKit❗❌:         }
    // RDKit❗❌:       } else {
    // RDKit❗❌:         gotCoords = EmbeddingOps::checkTetrahedralCenters(positions, eargs,
    // RDKit❗❌:                                                           embedParams);
    // RDKit❗❌:         if (!gotCoords) {
    // RDKit❗❌:           if (embedParams.trackFailures) {
    // RDKit❗❌: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit❗❌:             std::lock_guard<std::mutex> lock(GetFailMutex());
    // RDKit❗❌: #endif
    // RDKit❗❌:             embedParams
    // RDKit❗❌:                 .failures[EmbedFailureCauses::CHECK_TETRAHEDRAL_CENTERS]++;
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:
    // RDKit❗❌:       // Check if any of our chiral centers are badly out of whack.
    // RDKit❗❌:       if (gotCoords && embedParams.enforceChirality &&
    // RDKit❗❌:           eargs.chiralCenters->size() > 0) {
    // RDKit❗❌:         gotCoords =
    // RDKit❗❌:             EmbeddingOps::checkChiralCenters(positions, eargs, embedParams);
    // RDKit❗❌:         if (!gotCoords) {
    // RDKit❗❌:           if (embedParams.trackFailures) {
    // RDKit❗❌: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit❗❌:             std::lock_guard<std::mutex> lock(GetFailMutex());
    // RDKit❗❌: #endif
    // RDKit❗❌:             embedParams.failures[EmbedFailureCauses::CHECK_CHIRAL_CENTERS]++;
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       // redo the minimization if we have a chiral center
    // RDKit❗❌:       // or have started from random coords.
    // RDKit❗❌:       if (gotCoords &&
    // RDKit❗❌:           (eargs.chiralCenters->size() > 0 || embedParams.useRandomCoords)) {
    // RDKit❗❌:         if (ControlCHandler::getGotSignal()) {
    // RDKit❗❌:           return false;
    // RDKit❗❌:         }
    // RDKit❗❌:         gotCoords = EmbeddingOps::minimizeFourthDimension(
    // RDKit❗❌:             positions, eargs, embedParams, end_time);
    // RDKit❗❌:         if (!gotCoords) {
    // RDKit❗❌:           if (embedParams.trackFailures) {
    // RDKit❗❌: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit❗❌:             std::lock_guard<std::mutex> lock(GetFailMutex());
    // RDKit❗❌: #endif
    // RDKit❗❌:             if (end_time != nullptr && Clock::now() > *end_time) {
    // RDKit❗❌:               embedParams.failures[EmbedFailureCauses::EXCEEDED_TIMEOUT]++;
    // RDKit❗❌:             }
    // RDKit❗❌:             embedParams
    // RDKit❗❌:                 .failures[EmbedFailureCauses::MINIMIZE_FOURTH_DIMENSION]++;
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:
    // RDKit❗❌:       // (ET)(K)DG
    // RDKit❗❌:       if (gotCoords && (embedParams.useExpTorsionAnglePrefs ||
    // RDKit❗❌:                         embedParams.useBasicKnowledge)) {
    // RDKit❗❌:         if (ControlCHandler::getGotSignal()) {
    // RDKit❗❌:           return false;
    // RDKit❗❌:         }
    // RDKit❗❌:         gotCoords = EmbeddingOps::minimizeWithExpTorsions(*positions, eargs,
    // RDKit❗❌:                                                           embedParams);
    // RDKit❗❌:         if (!gotCoords) {
    // RDKit❗❌:           if (embedParams.trackFailures) {
    // RDKit❗❌: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit❗❌:             std::lock_guard<std::mutex> lock(GetFailMutex());
    // RDKit❗❌: #endif
    // RDKit❗❌:             embedParams.failures[EmbedFailureCauses::ETK_MINIMIZATION]++;
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       if (gotCoords) {
    // RDKit❗❌:         gotCoords = EmbeddingOps::doubleBondGeometryChecks(*positions, eargs,
    // RDKit❗❌:                                                            embedParams);
    // RDKit❗❌:         if (!gotCoords && embedParams.trackFailures) {
    // RDKit❗❌: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit❗❌:           std::lock_guard<std::mutex> lock(GetFailMutex());
    // RDKit❗❌: #endif
    // RDKit❗❌:           embedParams.failures[EmbedFailureCauses::LINEAR_DOUBLE_BOND]++;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       // test if stereo is correct
    // RDKit❗❌:       if (embedParams.enforceChirality && gotCoords) {
    // RDKit❗❌:         if (!eargs.chiralCenters->empty()) {
    // RDKit❗❌:           // test if chirality is correct. Any additional test failures
    // RDKit❗❌:           // will be tracked there if necessary.
    // RDKit❗❌:           gotCoords =
    // RDKit❗❌:               EmbeddingOps::finalChiralChecks(positions, eargs, embedParams);
    // RDKit❗❌:         }
    // RDKit❗❌:         if (gotCoords && !eargs.stereoDoubleBonds->empty()) {
    // RDKit❗❌:           gotCoords = EmbeddingOps::doubleBondStereoChecks(*positions, eargs,
    // RDKit❗❌:                                                            embedParams);
    // RDKit❗❌:           if (!gotCoords && embedParams.trackFailures) {
    // RDKit❗❌: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit❗❌:             std::lock_guard<std::mutex> lock(GetFailMutex());
    // RDKit❗❌: #endif
    // RDKit❗❌:             embedParams.failures[EmbedFailureCauses::BAD_DOUBLE_BOND_STEREO]++;
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:   }  // while
    // RDKit❗❌:
    // RDKit❗❌:   return gotCoords;
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION embedder_embed_points_with_rng

    use crate::EmbedFailureCause;
    // The global source RNG borrow acquires a Mutex once per seed=-1 operation;
    // this is extra synchronization relative to the source unsynchronized stream.
    let mut got_coords = false;
    let mut iter = 0;
    while !got_coords && iter < max_iterations {
        if let Some(deadline) = end_time
            && Instant::now() > deadline
        {
            break;
        }
        iter += 1;
        if let Some(callback) = params.callback {
            callback(iter);
        }
        if source_got_signal() {
            return Ok(false);
        }
        got_coords = embedder_generate_initial_coords(positions, eargs, params, dist_mat, rng)?;
        if !got_coords {
            embedder_increment_failure(params, progress, EmbedFailureCause::InitialCoords)?;
        } else {
            if source_got_signal() {
                return Ok(false);
            }
            got_coords =
                embedder_first_minimization_with_basin(positions, eargs, params, basin_thresh)?;
            if !got_coords {
                embedder_increment_failure(params, progress, EmbedFailureCause::FirstMinimization)?;
            } else {
                got_coords = embedder_check_tetrahedral_centers(positions, eargs, params)?;
                if !got_coords {
                    embedder_increment_failure(
                        params,
                        progress,
                        EmbedFailureCause::CheckTetrahedralCenters,
                    )?;
                }
            }
            if got_coords && params.enforce_chirality && !eargs.chiral_centers.is_empty() {
                got_coords = embedder_check_chiral_centers(positions, eargs, params)?;
                if !got_coords {
                    embedder_increment_failure(
                        params,
                        progress,
                        EmbedFailureCause::CheckChiralCenters,
                    )?;
                }
            }
            if got_coords && (!eargs.chiral_centers.is_empty() || params.use_random_coords) {
                if source_got_signal() {
                    return Ok(false);
                }
                got_coords = embedder_minimize_fourth_dimension_with_basin(
                    positions,
                    eargs,
                    params,
                    end_time,
                    basin_thresh,
                )?;
                if !got_coords {
                    if let Some(deadline) = end_time
                        && Instant::now() > deadline
                    {
                        embedder_increment_failure(
                            params,
                            progress,
                            EmbedFailureCause::ExceededTimeout,
                        )?;
                    }
                    embedder_increment_failure(
                        params,
                        progress,
                        EmbedFailureCause::MinimizeFourthDimension,
                    )?;
                }
            }
            if got_coords && (params.use_exp_torsion_angle_prefs || params.use_basic_knowledge) {
                if source_got_signal() {
                    return Ok(false);
                }
                got_coords = embedder_minimize_with_exp_torsions(positions, eargs, params)?;
                if !got_coords {
                    embedder_increment_failure(
                        params,
                        progress,
                        EmbedFailureCause::EtkMinimization,
                    )?;
                }
            }
            if got_coords {
                got_coords = embedder_double_bond_geometry_checks(positions, eargs, params, 1e-3)?;
                if !got_coords {
                    embedder_increment_failure(
                        params,
                        progress,
                        EmbedFailureCause::LinearDoubleBond,
                    )?;
                }
            }
            if params.enforce_chirality && got_coords {
                if !eargs.chiral_centers.is_empty() {
                    got_coords =
                        embedder_final_chiral_checks_shared(positions, eargs, params, progress)?;
                }
                if got_coords && !eargs.stereo_double_bonds.is_empty() {
                    got_coords = embedder_double_bond_stereo_checks(positions, eargs, params)?;
                    if !got_coords {
                        embedder_increment_failure(
                            params,
                            progress,
                            EmbedFailureCause::BadDoubleBondStereo,
                        )?;
                    }
                }
            }
        }
    }
    Ok(got_coords)
}
fn embedder_embed_points(
    positions: &mut [Vec<f64>],
    eargs: EmbedArgs<'_>,
    params: &mut EmbedParams,
    seed: i32,
    end_time: Option<Instant>,
) -> Result<bool, GenerationError> {
    with_embedding_progress(params, |settings, progress| {
        embedder_embed_points_shared(positions, eargs, settings, progress, seed, end_time)
    })
}
fn embedder_embed_points_shared(
    positions: &mut [Vec<f64>],
    eargs: EmbedArgs<'_>,
    params: &EmbedParams,
    progress: &EmbeddingProgress,
    seed: i32,
    end_time: Option<Instant>,
) -> Result<bool, GenerationError> {
    // BEGIN RDKIT CPP FUNCTION embedder_embed_points (GraphMol/DistGeomHelpers/Embedder.cpp)
    // RDKit❗❌: bool embedPoints(RDGeom::PointPtrVect *positions, detail::EmbedArgs eargs,
    // RDKit❗❌:                  EmbedParameters &embedParams, int seed, TimePoint *end_time) {
    // RDKit❗❌:   PRECONDITION(positions, "bogus positions");
    // RDKit❗❌:   if (embedParams.maxIterations == 0) {
    // RDKit❗❌:     embedParams.maxIterations = 10 * positions->size();
    // RDKit❗❌:   }
    // RDKit❗❌:   RDNumeric::DoubleSymmMatrix distMat(positions->size(), 0.0);
    // RDKit❗❌:
    // RDKit❗❌:   // The basin threshold just gets us into trouble when we're using
    // RDKit❗❌:   // random coordinates since it ends up ignoring 1-4 (and higher)
    // RDKit❗❌:   // interactions. This causes us to get folded-up (and self-penetrating)
    // RDKit❗❌:   // conformations for large flexible molecules
    // RDKit❗❌:   if (embedParams.useRandomCoords) {
    // RDKit❗❌:     embedParams.basinThresh = 1e8;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   std::unique_ptr<RDKit::rng_type> generator;
    // RDKit❗❌:   std::unique_ptr<RDKit::uniform_double> distrib;
    // RDKit❗❌:   std::unique_ptr<RDKit::double_source_type> rngMgr;
    // RDKit❗❌:
    // RDKit❗❌:   RDKit::double_source_type *rng = nullptr;
    // RDKit❗❌:   CHECK_INVARIANT(seed >= -1,
    // RDKit❗❌:                   "random seed must either be positive, zero, or negative one");
    // RDKit❗❌:   if (seed > -1) {
    // RDKit❗❌:     generator.reset(new RDKit::rng_type(42u));
    // RDKit❗❌:     generator->seed(seed);
    // RDKit❗❌:     distrib.reset(new RDKit::uniform_double(0.0, 1.0));
    // RDKit❗❌:     rngMgr.reset(new RDKit::double_source_type(*generator, *distrib));
    // RDKit❗❌:     rng = rngMgr.get();
    // RDKit❗❌:   } else {
    // RDKit❗❌:     rng = &RDKit::getDoubleRandomSource();
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   bool gotCoords = false;
    // RDKit❗❌:   unsigned int iter = 0;
    // RDKit❗❌:   while (!gotCoords && iter < embedParams.maxIterations) {
    // RDKit❗❌:     if (end_time != nullptr && Clock::now() > *end_time) {
    // RDKit❗❌:       break;
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     ++iter;
    // RDKit❗❌:     if (embedParams.callback != nullptr) {
    // RDKit❗❌:       embedParams.callback(iter);
    // RDKit❗❌:     }
    // RDKit❗❌:     if (ControlCHandler::getGotSignal()) {
    // RDKit❗❌:       return false;
    // RDKit❗❌:     }
    // RDKit❗❌:     gotCoords = EmbeddingOps::generateInitialCoords(positions, eargs,
    // RDKit❗❌:                                                     embedParams, distMat, rng);
    // RDKit❗❌:     if (!gotCoords) {
    // RDKit❗❌:       if (embedParams.trackFailures) {
    // RDKit❗❌: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit❗❌:         std::lock_guard<std::mutex> lock(GetFailMutex());
    // RDKit❗❌: #endif
    // RDKit❗❌:         embedParams.failures[EmbedFailureCauses::INITIAL_COORDS]++;
    // RDKit❗❌:       }
    // RDKit❗❌:     } else {
    // RDKit❗❌:       if (ControlCHandler::getGotSignal()) {
    // RDKit❗❌:         return false;
    // RDKit❗❌:       }
    // RDKit❗❌:       gotCoords =
    // RDKit❗❌:           EmbeddingOps::firstMinimization(positions, eargs, embedParams);
    // RDKit❗❌:       if (!gotCoords) {
    // RDKit❗❌:         if (embedParams.trackFailures) {
    // RDKit❗❌: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit❗❌:           std::lock_guard<std::mutex> lock(GetFailMutex());
    // RDKit❗❌: #endif
    // RDKit❗❌:           embedParams.failures[EmbedFailureCauses::FIRST_MINIMIZATION]++;
    // RDKit❗❌:         }
    // RDKit❗❌:       } else {
    // RDKit❗❌:         gotCoords = EmbeddingOps::checkTetrahedralCenters(positions, eargs,
    // RDKit❗❌:                                                           embedParams);
    // RDKit❗❌:         if (!gotCoords) {
    // RDKit❗❌:           if (embedParams.trackFailures) {
    // RDKit❗❌: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit❗❌:             std::lock_guard<std::mutex> lock(GetFailMutex());
    // RDKit❗❌: #endif
    // RDKit❗❌:             embedParams
    // RDKit❗❌:                 .failures[EmbedFailureCauses::CHECK_TETRAHEDRAL_CENTERS]++;
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:
    // RDKit❗❌:       // Check if any of our chiral centers are badly out of whack.
    // RDKit❗❌:       if (gotCoords && embedParams.enforceChirality &&
    // RDKit❗❌:           eargs.chiralCenters->size() > 0) {
    // RDKit❗❌:         gotCoords =
    // RDKit❗❌:             EmbeddingOps::checkChiralCenters(positions, eargs, embedParams);
    // RDKit❗❌:         if (!gotCoords) {
    // RDKit❗❌:           if (embedParams.trackFailures) {
    // RDKit❗❌: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit❗❌:             std::lock_guard<std::mutex> lock(GetFailMutex());
    // RDKit❗❌: #endif
    // RDKit❗❌:             embedParams.failures[EmbedFailureCauses::CHECK_CHIRAL_CENTERS]++;
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       // redo the minimization if we have a chiral center
    // RDKit❗❌:       // or have started from random coords.
    // RDKit❗❌:       if (gotCoords &&
    // RDKit❗❌:           (eargs.chiralCenters->size() > 0 || embedParams.useRandomCoords)) {
    // RDKit❗❌:         if (ControlCHandler::getGotSignal()) {
    // RDKit❗❌:           return false;
    // RDKit❗❌:         }
    // RDKit❗❌:         gotCoords = EmbeddingOps::minimizeFourthDimension(
    // RDKit❗❌:             positions, eargs, embedParams, end_time);
    // RDKit❗❌:         if (!gotCoords) {
    // RDKit❗❌:           if (embedParams.trackFailures) {
    // RDKit❗❌: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit❗❌:             std::lock_guard<std::mutex> lock(GetFailMutex());
    // RDKit❗❌: #endif
    // RDKit❗❌:             if (end_time != nullptr && Clock::now() > *end_time) {
    // RDKit❗❌:               embedParams.failures[EmbedFailureCauses::EXCEEDED_TIMEOUT]++;
    // RDKit❗❌:             }
    // RDKit❗❌:             embedParams
    // RDKit❗❌:                 .failures[EmbedFailureCauses::MINIMIZE_FOURTH_DIMENSION]++;
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:
    // RDKit❗❌:       // (ET)(K)DG
    // RDKit❗❌:       if (gotCoords && (embedParams.useExpTorsionAnglePrefs ||
    // RDKit❗❌:                         embedParams.useBasicKnowledge)) {
    // RDKit❗❌:         if (ControlCHandler::getGotSignal()) {
    // RDKit❗❌:           return false;
    // RDKit❗❌:         }
    // RDKit❗❌:         gotCoords = EmbeddingOps::minimizeWithExpTorsions(*positions, eargs,
    // RDKit❗❌:                                                           embedParams);
    // RDKit❗❌:         if (!gotCoords) {
    // RDKit❗❌:           if (embedParams.trackFailures) {
    // RDKit❗❌: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit❗❌:             std::lock_guard<std::mutex> lock(GetFailMutex());
    // RDKit❗❌: #endif
    // RDKit❗❌:             embedParams.failures[EmbedFailureCauses::ETK_MINIMIZATION]++;
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       if (gotCoords) {
    // RDKit❗❌:         gotCoords = EmbeddingOps::doubleBondGeometryChecks(*positions, eargs,
    // RDKit❗❌:                                                            embedParams);
    // RDKit❗❌:         if (!gotCoords && embedParams.trackFailures) {
    // RDKit❗❌: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit❗❌:           std::lock_guard<std::mutex> lock(GetFailMutex());
    // RDKit❗❌: #endif
    // RDKit❗❌:           embedParams.failures[EmbedFailureCauses::LINEAR_DOUBLE_BOND]++;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       // test if stereo is correct
    // RDKit❗❌:       if (embedParams.enforceChirality && gotCoords) {
    // RDKit❗❌:         if (!eargs.chiralCenters->empty()) {
    // RDKit❗❌:           // test if chirality is correct. Any additional test failures
    // RDKit❗❌:           // will be tracked there if necessary.
    // RDKit❗❌:           gotCoords =
    // RDKit❗❌:               EmbeddingOps::finalChiralChecks(positions, eargs, embedParams);
    // RDKit❗❌:         }
    // RDKit❗❌:         if (gotCoords && !eargs.stereoDoubleBonds->empty()) {
    // RDKit❗❌:           gotCoords = EmbeddingOps::doubleBondStereoChecks(*positions, eargs,
    // RDKit❗❌:                                                            embedParams);
    // RDKit❗❌:           if (!gotCoords && embedParams.trackFailures) {
    // RDKit❗❌: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit❗❌:             std::lock_guard<std::mutex> lock(GetFailMutex());
    // RDKit❗❌: #endif
    // RDKit❗❌:             embedParams.failures[EmbedFailureCauses::BAD_DOUBLE_BOND_STEREO]++;
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:   }  // while
    // RDKit❗❌:
    // RDKit❗❌:   return gotCoords;
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION embedder_embed_points

    let n = u32::try_from(positions.len())
        .map_err(|_| GenerationError::Input("point count exceeds source unsigned range"))?;
    let &(max_iterations, basin_thresh) = progress.defaults.get_or_init(|| {
        (
            if params.max_iterations == 0 {
                n.wrapping_mul(10)
            } else {
                params.max_iterations
            },
            if params.use_random_coords {
                1e8
            } else {
                params.basin_thresh
            },
        )
    });
    let mut dist_mat = SymmMatrix::new(positions.len());
    if seed < -1 {
        return Err(GenerationError::Input(
            "random seed must either be positive, zero, or negative one",
        ));
    }
    if seed >= 0 {
        let mut rng = cosmolkit_core::RdkitRandomEngine::from_seed(seed as u32);
        embedder_embed_points_with_rng(
            positions,
            &eargs,
            params,
            progress,
            max_iterations,
            basin_thresh,
            end_time,
            &mut dist_mat,
            &mut rng,
        )
    } else {
        cosmolkit_core::with_rdkit_random_generator(-1, |rng| {
            embedder_embed_points_with_rng(
                positions,
                &eargs,
                params,
                progress,
                max_iterations,
                basin_thresh,
                end_time,
                &mut dist_mat,
                rng,
            )
        })
    }
}

#[cfg(test)]
mod original_complete_flow_conditions {
    use super::*;
    use crate::EmbedFailureCause;
    use std::{
        sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        },
        time::Duration,
    };
    static EMBED_POINTS_CALLBACK_COUNT: AtomicUsize = AtomicUsize::new(0);
    fn embed_points_test_callback(iter: u32) {
        EMBED_POINTS_CALLBACK_COUNT.fetch_add(iter as usize, Ordering::SeqCst);
    }
    #[test]
    fn embedder_double_bond_geometry_checks_rejects_linear_arrangement() {
        let mmat = BoundsMatrix::new(3).expect("original bounds dimension");
        let chiral_centers: Vec<ChiralSetPtr> = Vec::new();
        let tetrahedral_carbons: Vec<ChiralSetPtr> = Vec::new();
        let double_bond_ends = vec![(0_usize, 1_usize, 2_usize)];
        let eargs = EmbedArgs {
            mmat: &mmat,
            chiral_centers: &chiral_centers,
            tetrahedral_carbons: &tetrahedral_carbons,
            etkdg_details: None,
            double_bond_ends: Some(&double_bond_ends),
            stereo_double_bonds: &[],
        };
        let mut params = EmbedParams::default();
        let linear = vec![
            vec![0.0, 0.0, 0.0],
            vec![1.0, 0.0, 0.0],
            vec![2.0, 0.0, 0.0],
        ];
        let bent = vec![
            vec![0.0, 0.0, 0.0],
            vec![1.0, 0.0, 0.0],
            vec![1.0, 1.0, 0.0],
        ];

        assert!(
            !embedder_double_bond_geometry_checks(&linear, &eargs, &mut params, 1.0e-3)
                .expect("original checked flow")
        );
        assert!(
            embedder_double_bond_geometry_checks(&bent, &eargs, &mut params, 1.0e-3)
                .expect("original checked flow")
        );
    }

    #[test]
    fn embedder_double_bond_geometry_checks_accepts_missing_double_bond_ends() {
        let mmat = BoundsMatrix::new(3).expect("original bounds dimension");
        let chiral_centers: Vec<ChiralSetPtr> = Vec::new();
        let tetrahedral_carbons: Vec<ChiralSetPtr> = Vec::new();
        let eargs = EmbedArgs {
            mmat: &mmat,
            chiral_centers: &chiral_centers,
            tetrahedral_carbons: &tetrahedral_carbons,
            etkdg_details: None,
            double_bond_ends: None,
            stereo_double_bonds: &[],
        };
        let mut params = EmbedParams::default();
        let positions = vec![
            vec![0.0, 0.0, 0.0],
            vec![1.0, 0.0, 0.0],
            vec![2.0, 0.0, 0.0],
        ];

        assert!(
            embedder_double_bond_geometry_checks(&positions, &eargs, &mut params, 1.0e-3)
                .expect("original checked flow")
        );
    }

    #[test]
    fn embedder_double_bond_stereo_checks_uses_dihedral_sign_rule() {
        let mmat = BoundsMatrix::new(4).expect("original bounds dimension");
        let chiral_centers: Vec<ChiralSetPtr> = Vec::new();
        let tetrahedral_carbons: Vec<ChiralSetPtr> = Vec::new();
        let trans_bonds = vec![(vec![0_usize, 1_usize, 2_usize, 3_usize], 1_i32)];
        let cis_bonds = vec![(vec![0_usize, 1_usize, 2_usize, 3_usize], -1_i32)];
        let trans_args = EmbedArgs {
            mmat: &mmat,
            chiral_centers: &chiral_centers,
            tetrahedral_carbons: &tetrahedral_carbons,
            etkdg_details: None,
            double_bond_ends: None,
            stereo_double_bonds: &trans_bonds,
        };
        let cis_args = EmbedArgs {
            mmat: &mmat,
            chiral_centers: &chiral_centers,
            tetrahedral_carbons: &tetrahedral_carbons,
            etkdg_details: None,
            double_bond_ends: None,
            stereo_double_bonds: &cis_bonds,
        };
        let mut params = EmbedParams::default();
        let cis_positions = vec![
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 0.0],
            vec![1.0, 0.0, 0.0],
            vec![1.0, 1.0, 0.0],
        ];
        let trans_positions = vec![
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 0.0],
            vec![1.0, 0.0, 0.0],
            vec![1.0, -1.0, 0.0],
        ];

        assert!(
            embedder_double_bond_stereo_checks(&cis_positions, &cis_args, &mut params)
                .expect("original checked flow")
        );
        assert!(
            !embedder_double_bond_stereo_checks(&cis_positions, &trans_args, &mut params)
                .expect("original checked flow")
        );
        assert!(
            embedder_double_bond_stereo_checks(&trans_positions, &trans_args, &mut params)
                .expect("original checked flow")
        );
    }

    fn broad_bounds_matrix(size: usize) -> BoundsMatrix {
        let mut mmat = BoundsMatrix::new(size).expect("original bounds dimension");
        for i in 1..size {
            for j in 0..i {
                mmat.set_lower(i, j, 0.0).expect("set lower");
                mmat.set_upper(i, j, 100.0).expect("set upper");
            }
        }
        mmat
    }

    #[test]
    fn embedder_final_chiral_checks_accepts_valid_final_chirality() {
        let mmat = broad_bounds_matrix(5);
        let chiral_set: ChiralSetPtr = Arc::new(ChiralSet::with_default_structure_flags(
            0, 1, 2, 3, 4, -100.0, 100.0,
        ));
        let chiral_centers: Vec<ChiralSetPtr> = vec![chiral_set];
        let tetrahedral_carbons: Vec<ChiralSetPtr> = Vec::new();
        let eargs = EmbedArgs {
            mmat: &mmat,
            chiral_centers: &chiral_centers,
            tetrahedral_carbons: &tetrahedral_carbons,
            etkdg_details: None,
            double_bond_ends: None,
            stereo_double_bonds: &[],
        };
        let mut params = EmbedParams::default();
        let mut positions = vec![
            vec![0.25, 0.25, 0.25],
            vec![0.0, 0.0, 0.0],
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];

        assert!(
            embedder_final_chiral_checks(&mut positions, &eargs, &mut params)
                .expect("original checked flow")
        );
    }

    #[test]
    fn embedder_final_chiral_checks_tracks_volume_failure() {
        let mmat = broad_bounds_matrix(5);
        let chiral_set: ChiralSetPtr = Arc::new(ChiralSet::with_default_structure_flags(
            0, 1, 2, 3, 4, 2.0, 3.0,
        ));
        let chiral_centers: Vec<ChiralSetPtr> = vec![chiral_set];
        let tetrahedral_carbons: Vec<ChiralSetPtr> = Vec::new();
        let eargs = EmbedArgs {
            mmat: &mmat,
            chiral_centers: &chiral_centers,
            tetrahedral_carbons: &tetrahedral_carbons,
            etkdg_details: None,
            double_bond_ends: None,
            stereo_double_bonds: &[],
        };
        let mut params = EmbedParams {
            track_failures: true,
            failures: vec![0; EmbedFailureCause::EndOfEnum as usize],
            ..EmbedParams::default()
        };
        let mut positions = vec![
            vec![0.1, 0.1, 0.1],
            vec![0.0, 0.0, 0.0],
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];

        assert!(
            !embedder_final_chiral_checks(&mut positions, &eargs, &mut params)
                .expect("original checked flow")
        );
        assert_eq!(
            params.failures[EmbedFailureCause::CheckChiralCenters2 as usize],
            1
        );
    }

    #[test]
    fn embedder_final_chiral_checks_tracks_bounds_failure() {
        let mut mmat = broad_bounds_matrix(5);
        mmat.set_upper(0, 1, 0.05).expect("set upper");
        let chiral_set: ChiralSetPtr = Arc::new(ChiralSet::with_default_structure_flags(
            0, 1, 2, 3, 4, -100.0, 100.0,
        ));
        let chiral_centers: Vec<ChiralSetPtr> = vec![chiral_set];
        let tetrahedral_carbons: Vec<ChiralSetPtr> = Vec::new();
        let eargs = EmbedArgs {
            mmat: &mmat,
            chiral_centers: &chiral_centers,
            tetrahedral_carbons: &tetrahedral_carbons,
            etkdg_details: None,
            double_bond_ends: None,
            stereo_double_bonds: &[],
        };
        let mut params = EmbedParams {
            track_failures: true,
            failures: vec![0; EmbedFailureCause::EndOfEnum as usize],
            ..EmbedParams::default()
        };
        let mut positions = vec![
            vec![0.25, 0.25, 0.25],
            vec![0.0, 0.0, 0.0],
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];

        assert!(
            !embedder_final_chiral_checks(&mut positions, &eargs, &mut params)
                .expect("original checked flow")
        );
        assert_eq!(
            params.failures[EmbedFailureCause::FinalChiralBounds as usize],
            1
        );
    }

    #[test]
    fn embedder_final_chiral_checks_tracks_center_in_volume_failure() {
        let mmat = broad_bounds_matrix(5);
        let chiral_set: ChiralSetPtr = Arc::new(ChiralSet::with_default_structure_flags(
            0, 1, 2, 3, 4, -100.0, 100.0,
        ));
        let chiral_centers: Vec<ChiralSetPtr> = vec![chiral_set];
        let tetrahedral_carbons: Vec<ChiralSetPtr> = Vec::new();
        let eargs = EmbedArgs {
            mmat: &mmat,
            chiral_centers: &chiral_centers,
            tetrahedral_carbons: &tetrahedral_carbons,
            etkdg_details: None,
            double_bond_ends: None,
            stereo_double_bonds: &[],
        };
        let mut params = EmbedParams {
            track_failures: true,
            failures: vec![0; EmbedFailureCause::EndOfEnum as usize],
            ..EmbedParams::default()
        };
        let mut positions = vec![
            vec![4.0, 4.0, 4.0],
            vec![0.0, 0.0, 0.0],
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];

        assert!(
            !embedder_final_chiral_checks(&mut positions, &eargs, &mut params)
                .expect("original checked flow")
        );
        assert_eq!(
            params.failures[EmbedFailureCause::FinalCenterInVolume as usize],
            1
        );
    }

    #[test]
    fn embedder_embed_points_sets_default_iterations_and_runs_callback() {
        EMBED_POINTS_CALLBACK_COUNT.store(0, Ordering::SeqCst);
        let mut mmat = BoundsMatrix::new(2).expect("original bounds dimension");
        mmat.set_lower(0, 1, 1.0).expect("set lower");
        mmat.set_upper(0, 1, 1.0).expect("set upper");
        let chiral_centers: Vec<ChiralSetPtr> = Vec::new();
        let tetrahedral_carbons: Vec<ChiralSetPtr> = Vec::new();
        let eargs = EmbedArgs {
            mmat: &mmat,
            chiral_centers: &chiral_centers,
            tetrahedral_carbons: &tetrahedral_carbons,
            etkdg_details: None,
            double_bond_ends: None,
            stereo_double_bonds: &[],
        };
        let mut coord_map = BTreeMap::new();
        coord_map.insert(0, [0.0, 0.0, 0.0]);
        coord_map.insert(1, [1.0, 0.0, 0.0]);
        let mut params = EmbedParams {
            callback: Some(embed_points_test_callback),
            use_random_coords: true,
            coord_map: Some(coord_map),
            ..EmbedParams::default()
        };
        let mut positions = vec![vec![0.0; 3], vec![0.0; 3]];

        assert!(
            embedder_embed_points(&mut positions, eargs, &mut params, 0, None)
                .expect("embed points")
        );
        assert_eq!(params.max_iterations, 20);
        assert_eq!(EMBED_POINTS_CALLBACK_COUNT.load(Ordering::SeqCst), 1);
    }

    #[test]
    fn embedder_embed_points_seed_zero_is_local_and_reproducible() {
        let mut mmat = BoundsMatrix::new(3).expect("original bounds dimension");
        for i in 1..3 {
            for j in 0..i {
                mmat.set_lower(i, j, 1.0).expect("set lower");
                mmat.set_upper(i, j, 2.0).expect("set upper");
            }
        }
        let chiral_centers: Vec<ChiralSetPtr> = Vec::new();
        let tetrahedral_carbons: Vec<ChiralSetPtr> = Vec::new();
        let mut params_a = EmbedParams::default();
        let mut params_b = EmbedParams::default();
        let mut positions_a = vec![vec![0.0; 3], vec![0.0; 3], vec![0.0; 3]];
        let mut positions_b = positions_a.clone();
        let eargs_a = EmbedArgs {
            mmat: &mmat,
            chiral_centers: &chiral_centers,
            tetrahedral_carbons: &tetrahedral_carbons,
            etkdg_details: None,
            double_bond_ends: None,
            stereo_double_bonds: &[],
        };
        let eargs_b = EmbedArgs {
            mmat: &mmat,
            chiral_centers: &chiral_centers,
            tetrahedral_carbons: &tetrahedral_carbons,
            etkdg_details: None,
            double_bond_ends: None,
            stereo_double_bonds: &[],
        };

        assert!(
            embedder_embed_points(&mut positions_a, eargs_a, &mut params_a, 0, None)
                .expect("embed points")
        );
        assert!(
            embedder_embed_points(&mut positions_b, eargs_b, &mut params_b, 0, None)
                .expect("embed points")
        );
        assert_eq!(positions_a, positions_b);
    }

    #[test]
    fn embedder_embed_points_timeout_before_first_iteration_returns_false_without_callback() {
        let mmat = BoundsMatrix::new(2).expect("original bounds dimension");
        let chiral_centers: Vec<ChiralSetPtr> = Vec::new();
        let tetrahedral_carbons: Vec<ChiralSetPtr> = Vec::new();
        let eargs = EmbedArgs {
            mmat: &mmat,
            chiral_centers: &chiral_centers,
            tetrahedral_carbons: &tetrahedral_carbons,
            etkdg_details: None,
            double_bond_ends: None,
            stereo_double_bonds: &[],
        };
        let mut params = EmbedParams::default();
        let mut positions = vec![vec![0.0; 3], vec![0.0; 3]];

        assert!(
            !embedder_embed_points(
                &mut positions,
                eargs,
                &mut params,
                1,
                Some(Instant::now() - Duration::from_secs(1))
            )
            .expect("embed points")
        );
        assert_eq!(params.max_iterations, 20);
    }

    #[test]
    fn embedder_embed_points_tracks_linear_double_bond_failure() {
        let mmat = broad_bounds_matrix(3);
        let chiral_centers: Vec<ChiralSetPtr> = Vec::new();
        let tetrahedral_carbons: Vec<ChiralSetPtr> = Vec::new();
        let double_bond_ends = vec![(0_usize, 1_usize, 2_usize)];
        let eargs = EmbedArgs {
            mmat: &mmat,
            chiral_centers: &chiral_centers,
            tetrahedral_carbons: &tetrahedral_carbons,
            etkdg_details: None,
            double_bond_ends: Some(&double_bond_ends),
            stereo_double_bonds: &[],
        };
        let mut coord_map = BTreeMap::new();
        coord_map.insert(0, [0.0, 0.0, 0.0]);
        coord_map.insert(1, [1.0, 0.0, 0.0]);
        coord_map.insert(2, [2.0, 0.0, 0.0]);
        let mut params = EmbedParams {
            max_iterations: 1,
            use_random_coords: true,
            coord_map: Some(coord_map),
            track_failures: true,
            failures: vec![0; EmbedFailureCause::EndOfEnum as usize],
            ..EmbedParams::default()
        };
        let mut positions = vec![vec![0.0; 3], vec![0.0; 3], vec![0.0; 3]];

        assert!(
            !embedder_embed_points(&mut positions, eargs, &mut params, 7, None)
                .expect("embed points")
        );
        assert_eq!(
            params.failures[EmbedFailureCause::LinearDoubleBond as usize],
            1
        );
    }
}

#[cfg(all(test, not(target_arch = "wasm32")))]
mod source_interrupt_conditions {
    use super::*;
    #[test]
    fn source_interrupt_native_child_reset_and_sigint() {
        const KEY: &str = "COSMOLKIT_CONFORMER_SOURCE_INTERRUPT_CHILD";
        if std::env::var_os(KEY).is_some() {
            source_reset_interrupt().expect("source reset");
            assert!(!source_got_signal());
            // Safety: child process installed the source handler; the synchronous raise
            // delivers SIGINT only to this child and cannot affect the parent test process.
            assert_eq!(unsafe { libc::raise(libc::SIGINT) }, 0);
            assert!(source_got_signal());
            source_reset_interrupt().expect("source second reset");
            assert!(!source_got_signal());
        } else {
            let status=std::process::Command::new(std::env::current_exe().expect("native test image"))
    .args(["--exact","generation::source_interrupt_conditions::source_interrupt_native_child_reset_and_sigint","--nocapture"])
    .env(KEY,"1").status().expect("native child");
            assert!(status.success(), "real native source SIGINT child failed");
        }
    }
}

use cosmolkit_core::{RingInfo, ValenceAssignment};
use cosmolkit_model::{
    AtomId, BondId, BondOrder, BondStereo, ChiralTag, Hybridization, TopologyBlock,
};

/// Prepared chemical rows are borrowed from their unique core owners.
struct PreparedEmbeddingTopology<'a> {
    topology: &'a TopologyBlock,
    rings: &'a RingInfo,
    valence: &'a ValenceAssignment,
    hybridizations: &'a [Hybridization],
    conjugated: &'a [bool],
}
#[derive(Debug, thiserror::Error)]
pub enum GenerationFailure {
    #[error(transparent)]
    TopologyBounds(#[from] crate::graph_bounds::GraphBoundsError),
    #[error(transparent)]
    TorsionPreferences(#[from] cosmolkit_forcefields::CrystalffTorsionPreferencesError),
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PreparationDiagnostic {
    AtropisomerBondSkipped { bond: BondId },
    BoundsSmoothingFailed,
}
impl std::fmt::Display for PreparationDiagnostic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::AtropisomerBondSkipped { bond } => write!(
                f,
                "Atropisomer bond stereochemistry not used for bond {}, which does not have two exo substituents on at least one side.",
                bond.index()
            ),
            Self::BoundsSmoothingFailed => {
                f.write_str("Could not triangle bounds smooth molecule.")
            }
        }
    }
}
fn embedder_find_double_bonds(
    mol: &TopologyBlock,
    double_bond_ends: &mut Vec<(usize, usize, usize)>,
    stereo_double_bonds: &mut Vec<(Vec<usize>, i32)>,
    coord_map: Option<&BTreeMap<i32, [f64; 3]>>,
) -> Result<(), GenerationError> {
    // BEGIN RDKIT CPP FUNCTION embedder_find_double_bonds (GraphMol/DistGeomHelpers/Embedder.cpp)
    // RDKit❗✔️: RDKIT_DISTGEOMHELPERS_EXPORT void findDoubleBonds(
    // RDKit❗✔️:     const ROMol &mol,
    // RDKit❗✔️:     std::vector<std::tuple<unsigned int, unsigned int, unsigned int>>
    // RDKit❗✔️:         &doubleBondEnds,
    // RDKit❗✔️:     std::vector<std::pair<std::vector<unsigned int>, int>> &stereoDoubleBonds,
    // RDKit❗✔️:     const std::map<int, RDGeom::Point3D> *coordMap) {
    // RDKit❗✔️:   doubleBondEnds.clear();
    // RDKit❗✔️:   stereoDoubleBonds.clear();
    // RDKit❗✔️:   for (const auto bnd : mol.bonds()) {
    // RDKit❗✔️:     if (bnd->getBondType() == Bond::BondType::DOUBLE) {
    // RDKit❗✔️:       for (const auto atm : {bnd->getBeginAtom(), bnd->getEndAtom()}) {
    // RDKit❗✔️:         if (atm->getDegree() < 2) {
    // RDKit❗✔️:           continue;
    // RDKit❗✔️:         }
    // RDKit❗✔️:         auto oatm = bnd->getOtherAtom(atm);
    // RDKit❗✔️:         for (const auto nbr : mol.atomNeighbors(atm)) {
    // RDKit❗✔️:           if (nbr == oatm) {
    // RDKit❗✔️:             continue;
    // RDKit❗✔️:           }
    // RDKit❗✔️:           const auto obnd =
    // RDKit❗✔️:               mol.getBondBetweenAtoms(atm->getIdx(), nbr->getIdx());
    // RDKit❗✔️:           if (!obnd || (obnd->getBondType() != Bond::BondType::SINGLE &&
    // RDKit❗✔️:                         atm->getDegree() == 2)) {
    // RDKit❗✔️:             continue;
    // RDKit❗✔️:           }
    // RDKit❗✔️:           doubleBondEnds.emplace_back(nbr->getIdx(), atm->getIdx(),
    // RDKit❗✔️:                                       oatm->getIdx());
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:       // if there's stereo, handle that too:
    // RDKit❗✔️:       if (bnd->getStereo() > Bond::BondStereo::STEREOANY) {
    // RDKit❗✔️:         // only do this if the controlling atoms aren't in the coord map
    // RDKit❗✔️:         if (coordMap &&
    // RDKit❗✔️:             coordMap->find(bnd->getStereoAtoms()[0]) != coordMap->end() &&
    // RDKit❗✔️:             coordMap->find(bnd->getStereoAtoms()[1]) != coordMap->end()) {
    // RDKit❗✔️:           continue;
    // RDKit❗✔️:         }
    // RDKit❗✔️:         int sign = 1;
    // RDKit❗✔️:         if (bnd->getStereo() == Bond::BondStereo::STEREOCIS ||
    // RDKit❗✔️:             bnd->getStereo() == Bond::BondStereo::STEREOZ) {
    // RDKit❗✔️:           sign = -1;
    // RDKit❗✔️:         }
    // RDKit❗✔️:         std::pair<std::vector<unsigned int>, int> elem{
    // RDKit❗✔️:             {static_cast<unsigned>(bnd->getStereoAtoms()[0]),
    // RDKit❗✔️:              bnd->getBeginAtomIdx(), bnd->getEndAtomIdx(),
    // RDKit❗✔️:              static_cast<unsigned>(bnd->getStereoAtoms()[1])},
    // RDKit❗✔️:             sign};
    // RDKit❗✔️:         stereoDoubleBonds.push_back(elem);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION embedder_find_double_bonds

    double_bond_ends.clear();
    stereo_double_bonds.clear();
    for bond in &mol.bonds {
        if bond.order() != BondOrder::Double {
            continue;
        }
        for atom in [bond.begin(), bond.end()] {
            let neighbors = mol.adjacency.neighbors_of(atom.index());
            if neighbors.len() < 2 {
                continue;
            }
            let other = if atom == bond.begin() {
                bond.end()
            } else {
                bond.begin()
            };
            for nbr in neighbors {
                if nbr.atom_index == other.index() {
                    continue;
                }
                let nb = &mol.bonds[nbr.bond.index()];
                if nb.order() != BondOrder::Single && neighbors.len() == 2 {
                    continue;
                }
                double_bond_ends.push((nbr.atom_index, atom.index(), other.index()));
            }
        }
        if !matches!(bond.stereo(), BondStereo::None | BondStereo::Any) {
            let stereo = bond.stereo_atoms().ok_or(GenerationError::Input(
                "stereo double bond requires two controlling atoms",
            ))?;
            if coord_map.is_some_and(|map| {
                map.contains_key(&(stereo[0].index() as i32))
                    && map.contains_key(&(stereo[1].index() as i32))
            }) {
                continue;
            }
            let sign = if matches!(bond.stereo(), BondStereo::Cis | BondStereo::Z) {
                -1
            } else {
                1
            };
            stereo_double_bonds.push((
                vec![
                    stereo[0].index(),
                    bond.begin().index(),
                    bond.end().index(),
                    stereo[1].index(),
                ],
                sign,
            ));
        }
    }
    Ok(())
}
fn embedder_find_chiral_sets(
    mol: &TopologyBlock,
    ring_info: Option<&RingInfo>,
    chiral_centers: &mut Vec<ChiralSetPtr>,
    tetrahedral_centers: &mut Vec<ChiralSetPtr>,
    coord_map: Option<&BTreeMap<i32, [f64; 3]>>,
    diagnostics: &mut Vec<PreparationDiagnostic>,
) -> Result<(), GenerationError> {
    // BEGIN RDKIT CPP FUNCTION embedder_find_chiral_sets (GraphMol/DistGeomHelpers/Embedder.cpp)
    // RDKit❗❌: void findChiralSets(const ROMol &mol, DistGeom::VECT_CHIRALSET &chiralCenters,
    // RDKit❗❌:                     DistGeom::VECT_CHIRALSET &tetrahedralCenters,
    // RDKit❗❌:                     const std::map<int, RDGeom::Point3D> *coordMap) {
    // RDKit❗❌:   for (const auto &atom : mol.atoms()) {
    // RDKit❗❌:     if (atom->getAtomicNum() != 1) {  // skip hydrogens
    // RDKit❗❌:       Atom::ChiralType chiralType = atom->getChiralTag();
    // RDKit❗❌:       if ((chiralType == Atom::CHI_TETRAHEDRAL_CW ||
    // RDKit❗❌:            chiralType == Atom::CHI_TETRAHEDRAL_CCW) ||
    // RDKit❗❌:           ((atom->getAtomicNum() == 6 || atom->getAtomicNum() == 7) &&
    // RDKit❗❌:            atom->getDegree() == 4)) {
    // RDKit❗❌:         // make a chiral set from the neighbors
    // RDKit❗❌:         INT_VECT nbrs;
    // RDKit❗❌:         nbrs.reserve(4);
    // RDKit❗❌:         // find the neighbors of this atom and enter them into the
    // RDKit❗❌:         // nbr list
    // RDKit❗❌:         ROMol::OEDGE_ITER beg, end;
    // RDKit❗❌:         boost::tie(beg, end) = mol.getAtomBonds(atom);
    // RDKit❗❌:         while (beg != end) {
    // RDKit❗❌:           nbrs.push_back(mol[*beg]->getOtherAtom(atom)->getIdx());
    // RDKit❗❌:           ++beg;
    // RDKit❗❌:         }
    // RDKit❗❌:         // if we have less than 4 heavy atoms as neighbors,
    // RDKit❗❌:         // we need to include the chiral center into the mix
    // RDKit❗❌:         // we should at least have 3 though
    // RDKit❗❌:         CHECK_INVARIANT(nbrs.size() >= 3, "Cannot be a chiral center");
    // RDKit❗❌:
    // RDKit❗❌:         double volLowerBound = 5.0;
    // RDKit❗❌:         double volUpperBound = 100.0;
    // RDKit❗❌:
    // RDKit❗❌:         if (nbrs.size() < 4) {
    // RDKit❗❌:           // we get lower volumes if there are three neighbors,
    // RDKit❗❌:           //  this was github #5883
    // RDKit❗❌:           volLowerBound = 2.0;
    // RDKit❗❌:           nbrs.insert(nbrs.end(), atom->getIdx());
    // RDKit❗❌:         }
    // RDKit❗❌:
    // RDKit❗❌:         // set a flag for tetrahedral centers that are in multiple small rings
    // RDKit❗❌:         auto numSmallRings = 0u;
    // RDKit❗❌:         constexpr int smallRingSize = 5;
    // RDKit❗❌:         for (const auto sz : mol.getRingInfo()->atomRingSizes(atom->getIdx())) {
    // RDKit❗❌:           if (sz < smallRingSize) {
    // RDKit❗❌:             ++numSmallRings;
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:         std::uint64_t structureFlags = 0;
    // RDKit❗❌:         if (numSmallRings > 1) {
    // RDKit❗❌:           structureFlags = static_cast<std::uint64_t>(
    // RDKit❗❌:               DistGeom::ChiralSetStructureFlags::IN_FUSED_SMALL_RINGS);
    // RDKit❗❌:         }
    // RDKit❗❌:
    // RDKit❗❌:         // now create a chiral set and set the upper and lower bound on the
    // RDKit❗❌:         // volume
    // RDKit❗❌:         if (chiralType == Atom::CHI_TETRAHEDRAL_CCW) {
    // RDKit❗❌:           // positive chiral volume
    // RDKit❗❌:           auto *cset = new DistGeom::ChiralSet(atom->getIdx(), nbrs[0], nbrs[1],
    // RDKit❗❌:                                                nbrs[2], nbrs[3], volLowerBound,
    // RDKit❗❌:                                                volUpperBound, structureFlags);
    // RDKit❗❌:           DistGeom::ChiralSetPtr cptr(cset);
    // RDKit❗❌:           chiralCenters.push_back(cptr);
    // RDKit❗❌:         } else if (chiralType == Atom::CHI_TETRAHEDRAL_CW) {
    // RDKit❗❌:           auto *cset = new DistGeom::ChiralSet(atom->getIdx(), nbrs[0], nbrs[1],
    // RDKit❗❌:                                                nbrs[2], nbrs[3], -volUpperBound,
    // RDKit❗❌:                                                -volLowerBound, structureFlags);
    // RDKit❗❌:           DistGeom::ChiralSetPtr cptr(cset);
    // RDKit❗❌:           chiralCenters.push_back(cptr);
    // RDKit❗❌:         } else {
    // RDKit❗❌:           if ((coordMap && coordMap->find(atom->getIdx()) != coordMap->end()) ||
    // RDKit❗❌:               (mol.getRingInfo()->isInitialized() &&
    // RDKit❗❌:                (mol.getRingInfo()->numAtomRings(atom->getIdx()) < 2 ||
    // RDKit❗❌:                 mol.getRingInfo()->isAtomInRingOfSize(atom->getIdx(), 3)))) {
    // RDKit❗❌:             // we only want to these tests for ring atoms that are not part of
    // RDKit❗❌:             // the coordMap
    // RDKit❗❌:             // there's no sense doing 3-rings because those are a nightmare
    // RDKit❗❌:           } else {
    // RDKit❗❌:             auto *cset = new DistGeom::ChiralSet(atom->getIdx(), nbrs[0],
    // RDKit❗❌:                                                  nbrs[1], nbrs[2], nbrs[3], 0.0,
    // RDKit❗❌:                                                  0.0, structureFlags);
    // RDKit❗❌:             DistGeom::ChiralSetPtr cptr(cset);
    // RDKit❗❌:             tetrahedralCenters.push_back(cptr);
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }  // if block -chirality check
    // RDKit❗❌:     }    // if block - heavy atom check
    // RDKit❗❌:   }      // for loop over atoms
    // RDKit❗❌:
    // RDKit❗❌:   // now do atropisomers
    // RDKit❗❌:   for (const auto &bond : mol.bonds()) {
    // RDKit❗❌:     if (bond->getStereo() != Bond::BondStereo::STEREOATROPCCW &&
    // RDKit❗❌:         bond->getStereo() != Bond::BondStereo::STEREOATROPCW) {
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:     Atropisomers::AtropAtomAndBondVec atomsAndBonds[2];
    // RDKit❗❌:     Atropisomers::getAtropisomerAtomsAndBonds(bond, atomsAndBonds, mol);
    // RDKit❗❌:     // make a chiral set for the atropisomeric bond
    // RDKit❗❌:     // we start with only managing cases where there are two exo-substituents on
    // RDKit❗❌:     // at least one side
    // RDKit❗❌:     if (atomsAndBonds[0].second.size() != 2 &&
    // RDKit❗❌:         atomsAndBonds[1].second.size() != 2) {
    // RDKit❗❌:       BOOST_LOG(rdWarningLog)
    // RDKit❗❌:           << "Atropisomer bond stereochemistry not used for bond "
    // RDKit❗❌:           << bond->getIdx()
    // RDKit❗❌:           << ", which does not have two exo substituents on at least one side."
    // RDKit❗❌:           << std::endl;
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:     int idx0 = atomsAndBonds[0].first->getIdx();
    // RDKit❗❌:     int idx1 = atomsAndBonds[1].first->getIdx();
    // RDKit❗❌:
    // RDKit❗❌:     int nbr1 = atomsAndBonds[0].second[0]->getOtherAtomIdx(idx0);
    // RDKit❗❌:     int nbr2 = 0;
    // RDKit❗❌:     int nbr3 = 0;
    // RDKit❗❌:     int nbr4 = 0;
    // RDKit❗❌:     if (atomsAndBonds[0].second.size() == 2) {
    // RDKit❗❌:       nbr2 = atomsAndBonds[0].second[1]->getOtherAtomIdx(idx0);
    // RDKit❗❌:       nbr3 = atomsAndBonds[1].second[0]->getOtherAtomIdx(idx1);
    // RDKit❗❌:       if (atomsAndBonds[1].second.size() == 2) {
    // RDKit❗❌:         nbr4 = atomsAndBonds[1].second[1]->getOtherAtomIdx(idx1);
    // RDKit❗❌:       } else {
    // RDKit❗❌:         nbr4 = idx0;
    // RDKit❗❌:       }
    // RDKit❗❌:     } else {
    // RDKit❗❌:       nbr2 = atomsAndBonds[1].second[0]->getOtherAtomIdx(idx1);
    // RDKit❗❌:       nbr3 = atomsAndBonds[1].second[1]->getOtherAtomIdx(idx1);
    // RDKit❗❌:       nbr4 = idx0;
    // RDKit❗❌:     }
    // RDKit❗❌:     INT_VECT nbrs = {nbr1, nbr2, nbr3, nbr4};
    // RDKit❗❌:
    // RDKit❗❌:     // FIX: these numbers are empirical and should be revisited
    // RDKit❗❌:     double volLowerBound = 1.0;
    // RDKit❗❌:     double volUpperBound = 100.0;
    // RDKit❗❌:     if (bond->getStereo() == Bond::BondStereo::STEREOATROPCCW) {
    // RDKit❗❌:       std::swap(volLowerBound, volUpperBound);
    // RDKit❗❌:       volLowerBound *= -1;
    // RDKit❗❌:       volUpperBound *= -1;
    // RDKit❗❌:     }
    // RDKit❗❌:     auto *cset = new DistGeom::ChiralSet(idx0, nbrs[0], nbrs[1], nbrs[2],
    // RDKit❗❌:                                          nbrs[3], volLowerBound, volUpperBound);
    // RDKit❗❌:     DistGeom::ChiralSetPtr cptr(cset);
    // RDKit❗❌:     chiralCenters.push_back(cptr);
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION embedder_find_chiral_sets

    for atom in &mol.atoms {
        if atom.atomic_number() == 1 {
            continue;
        }
        let tag = atom.chiral_tag();
        let id = atom.id();
        let neighbors = mol.adjacency.neighbors_of(id.index());
        if matches!(tag, ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw)
            || (matches!(atom.atomic_number(), 6 | 7) && neighbors.len() == 4)
        {
            if neighbors.len() < 3 {
                return Err(GenerationError::Input("Cannot be a chiral center"));
            }
            let mut nbrs = Vec::with_capacity(4);
            nbrs.extend(neighbors.iter().map(|n| n.atom_index));
            let mut lower = 5.0;
            let upper = 100.0;
            if nbrs.len() < 4 {
                lower = 2.0;
                nbrs.push(id.index());
            }
            let small = ring_info
                .map(|rings| rings.atom_ring_sizes(id).iter().filter(|&&n| n < 5).count())
                .unwrap_or(0);
            let flags = if small > 1 {
                ChiralSetStructureFlags::InFusedSmallRings as u64
            } else {
                0
            };
            if tag == ChiralTag::TetrahedralCcw {
                chiral_centers.push(std::sync::Arc::new(ChiralSet::new(
                    id.index(),
                    nbrs[0],
                    nbrs[1],
                    nbrs[2],
                    nbrs[3],
                    lower,
                    upper,
                    flags,
                )));
            } else if tag == ChiralTag::TetrahedralCw {
                chiral_centers.push(std::sync::Arc::new(ChiralSet::new(
                    id.index(),
                    nbrs[0],
                    nbrs[1],
                    nbrs[2],
                    nbrs[3],
                    -upper,
                    -lower,
                    flags,
                )));
            } else {
                let mapped = coord_map.is_some_and(|map| map.contains_key(&(id.index() as i32)));
                let excluded = ring_info.is_some_and(|rings| {
                    rings.is_initialized()
                        && (rings.num_atom_rings(id) < 2 || rings.is_atom_in_ring_of_size(id, 3))
                });
                if !mapped && !excluded {
                    tetrahedral_centers.push(std::sync::Arc::new(ChiralSet::new(
                        id.index(),
                        nbrs[0],
                        nbrs[1],
                        nbrs[2],
                        nbrs[3],
                        0.0,
                        0.0,
                        flags,
                    )));
                }
            }
        }
    }
    // This narrow shared entrypoint validates topology once, then performs the
    // source adjacency/sorting in core. The selected-id/result vectors and
    // validation are extra work relative to direct C++ bond iteration.
    let axial: Vec<BondId> = mol
        .bonds
        .iter()
        .filter(|b| matches!(b.stereo(), BondStereo::AtropCw | BondStereo::AtropCcw))
        .map(|b| b.id())
        .collect();
    if axial.is_empty() {
        return Ok(());
    }
    let carriers = cosmolkit_core::atropisomer_carriers_for_bonds(mol, &axial)?;
    for (bond_id, ends) in carriers {
        let ends = ends.ok_or(GenerationError::Input(
            "Atropisomer bond has no exo substituents",
        ))?;
        let [left, right] = ends;
        let lb = left.carrier_bonds();
        let rb = right.carrier_bonds();
        if lb.len() != 2 && rb.len() != 2 {
            diagnostics.push(PreparationDiagnostic::AtropisomerBondSkipped { bond: bond_id });
            continue;
        }
        let idx0 = left.focus().index();
        let idx1 = right.focus().index();
        let other = |b: BondId, focus: usize| {
            let bond = &mol.bonds[b.index()];
            if bond.begin().index() == focus {
                bond.end().index()
            } else {
                bond.begin().index()
            }
        };
        let nbr1 = other(lb[0], idx0);
        let (nbr2, nbr3, nbr4) = if lb.len() == 2 {
            (
                other(lb[1], idx0),
                other(rb[0], idx1),
                if rb.len() == 2 {
                    other(rb[1], idx1)
                } else {
                    idx0
                },
            )
        } else {
            (other(rb[0], idx1), other(rb[1], idx1), idx0)
        };
        let (lower, upper) = if mol.bonds[bond_id.index()].stereo() == BondStereo::AtropCcw {
            (-100.0, -1.0)
        } else {
            (1.0, 100.0)
        };
        chiral_centers.push(std::sync::Arc::new(
            ChiralSet::with_default_structure_flags(idx0, nbr1, nbr2, nbr3, nbr4, lower, upper),
        ));
    }
    Ok(())
}
fn embedder_adjust_bounds_mat_from_coord_map(
    mmat: &mut BoundsMatrix,
    _num_atoms: usize,
    coord_map: &BTreeMap<i32, [f64; 3]>,
) -> Result<(), GenerationError> {
    // BEGIN RDKIT CPP FUNCTION embedder_adjust_bounds_mat_from_coord_map (GraphMol/DistGeomHelpers/Embedder.cpp)
    // RDKit❗✔️: void adjustBoundsMatFromCoordMap(
    // RDKit❗✔️:     DistGeom::BoundsMatPtr mmat, unsigned int,
    // RDKit❗✔️:     const std::map<int, RDGeom::Point3D> *coordMap) {
    // RDKit❗✔️:   for (auto iIt = coordMap->begin(); iIt != coordMap->end(); ++iIt) {
    // RDKit❗✔️:     unsigned int iIdx = iIt->first;
    // RDKit❗✔️:     const RDGeom::Point3D &iPoint = iIt->second;
    // RDKit❗✔️:     auto jIt = iIt;
    // RDKit❗✔️:     while (++jIt != coordMap->end()) {
    // RDKit❗✔️:       unsigned int jIdx = jIt->first;
    // RDKit❗✔️:       const RDGeom::Point3D &jPoint = jIt->second;
    // RDKit❗✔️:       double dist = (iPoint - jPoint).length();
    // RDKit❗✔️:       mmat->setUpperBound(iIdx, jIdx, dist);
    // RDKit❗✔️:       mmat->setLowerBound(iIdx, jIdx, dist);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION embedder_adjust_bounds_mat_from_coord_map

    let mut outer = coord_map.iter();
    while let Some((&i, p)) = outer.next() {
        // Iter cloning copies the borrowed tree cursor, not entries or points.
        for (&j, q) in outer.clone() {
            let dist = length(delta(*p, *q));
            // The source converts signed keys to unsigned indices. Invalid
            // indices reach the checked storage error instead of unchecked UB.
            mmat.set_upper(i as u32 as usize, j as u32 as usize, dist)?;
            mmat.set_lower(i as u32 as usize, j as u32 as usize, dist)?;
        }
    }
    Ok(())
}
fn embedder_init_etkdg(
    input: &PreparedEmbeddingTopology<'_>,
    params: &EmbedParams,
    details: &mut CrystalFFDetails,
) -> Result<(), GenerationError> {
    // BEGIN RDKIT CPP FUNCTION embedder_init_etkdg (GraphMol/DistGeomHelpers/Embedder.cpp)
    // RDKit❗✔️: void initETKDG(ROMol *mol, const EmbedParameters &params,
    // RDKit❗✔️:                ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails) {
    // RDKit❗✔️:   PRECONDITION(mol, "bad molecule");
    // RDKit❗✔️:   unsigned int nAtoms = mol->getNumAtoms();
    // RDKit❗✔️:   if (params.useExpTorsionAnglePrefs || params.useBasicKnowledge) {
    // RDKit❗✔️:     ForceFields::CrystalFF::getExperimentalTorsions(
    // RDKit❗✔️:         *mol, etkdgDetails, params.useExpTorsionAnglePrefs,
    // RDKit❗✔️:         params.useSmallRingTorsions, params.useMacrocycleTorsions,
    // RDKit❗✔️:         params.useBasicKnowledge, params.ETversion, params.verbose);
    // RDKit❗✔️:     etkdgDetails.atomNums.resize(nAtoms);
    // RDKit❗✔️:     for (unsigned int i = 0; i < nAtoms; ++i) {
    // RDKit❗✔️:       etkdgDetails.atomNums[i] = mol->getAtomWithIdx(i)->getAtomicNum();
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   etkdgDetails.boundsMatForceScaling = params.boundsMatForceScaling;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION embedder_init_etkdg

    if params.use_exp_torsion_angle_prefs || params.use_basic_knowledge {
        cosmolkit_forcefields::get_experimental_torsions_without_bonds(
            input.topology,
            input.rings,
            input.valence,
            details,
            params.use_exp_torsion_angle_prefs,
            params.use_small_ring_torsions,
            params.use_macrocycle_torsions,
            params.use_basic_knowledge,
            params.et_version,
            params.verbose,
        )
        .map_err(GenerationFailure::from)?;
        details.atom_nums.resize(input.topology.atoms.len(), 0);
        for (i, atom) in input.topology.atoms.iter().enumerate() {
            details.atom_nums[i] = i32::from(atom.atomic_number());
        }
    }
    details.bounds_mat_force_scaling = params.bounds_mat_force_scaling;
    Ok(())
}
fn embedder_setup_initial_bounds_matrix(
    input: &PreparedEmbeddingTopology<'_>,
    mmat: &mut BoundsMatrix,
    coord_map: Option<&BTreeMap<i32, [f64; 3]>>,
    params: &EmbedParams,
    details: &mut CrystalFFDetails,
    diagnostics: &mut Vec<PreparationDiagnostic>,
) -> Result<bool, GenerationError> {
    // BEGIN RDKIT CPP FUNCTION embedder_setup_initial_bounds_matrix (GraphMol/DistGeomHelpers/Embedder.cpp)
    // RDKit❗✔️: bool setupInitialBoundsMatrix(
    // RDKit❗✔️:     ROMol *mol, DistGeom::BoundsMatPtr mmat,
    // RDKit❗✔️:     const std::map<int, RDGeom::Point3D> *coordMap,
    // RDKit❗✔️:     const EmbedParameters &params,
    // RDKit❗✔️:     ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails) {
    // RDKit❗✔️:   PRECONDITION(mol, "bad molecule");
    // RDKit❗✔️:   unsigned int nAtoms = mol->getNumAtoms();
    // RDKit❗✔️:   if (params.useExpTorsionAnglePrefs || params.useBasicKnowledge) {
    // RDKit❗✔️:     setTopolBounds(*mol, mmat, etkdgDetails.bonds, etkdgDetails.angles, true,
    // RDKit❗✔️:                    false, params.useMacrocycle14config,
    // RDKit❗✔️:                    params.forceTransAmides);
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     setTopolBounds(*mol, mmat, true, false, params.useMacrocycle14config,
    // RDKit❗✔️:                    params.forceTransAmides);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   double tol = 0.0;
    // RDKit❗✔️:   if (coordMap) {
    // RDKit❗✔️:     adjustBoundsMatFromCoordMap(mmat, nAtoms, coordMap);
    // RDKit❗✔️:     tol = 0.05;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (!DistGeom::triangleSmoothBounds(mmat, tol)) {
    // RDKit❗✔️:     // ok this bound matrix failed to triangle smooth - re-compute the
    // RDKit❗✔️:     // bounds matrix without 15 bounds and with VDW scaling
    // RDKit❗✔️:     initBoundsMat(mmat);
    // RDKit❗✔️:     setTopolBounds(*mol, mmat, false, true, params.useMacrocycle14config,
    // RDKit❗✔️:                    params.forceTransAmides);
    // RDKit❗✔️:
    // RDKit❗✔️:     if (coordMap) {
    // RDKit❗✔️:       adjustBoundsMatFromCoordMap(mmat, nAtoms, coordMap);
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     // try triangle smoothing again
    // RDKit❗✔️:     if (!DistGeom::triangleSmoothBounds(mmat, tol)) {
    // RDKit❗✔️:       // ok, we're not going to be able to smooth this,
    // RDKit❗✔️:       if (params.ignoreSmoothingFailures) {
    // RDKit❗✔️:         // proceed anyway with the more relaxed bounds matrix
    // RDKit❗✔️:         initBoundsMat(mmat);
    // RDKit❗✔️:         setTopolBounds(*mol, mmat, false, true, params.useMacrocycle14config,
    // RDKit❗✔️:                        params.forceTransAmides);
    // RDKit❗✔️:
    // RDKit❗✔️:         if (coordMap) {
    // RDKit❗✔️:           adjustBoundsMatFromCoordMap(mmat, nAtoms, coordMap);
    // RDKit❗✔️:         }
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:             << "Could not triangle bounds smooth molecule." << std::endl;
    // RDKit❗✔️:         return false;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return true;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION embedder_setup_initial_bounds_matrix

    use crate::graph_bounds::{init_bounds_mat, set_topol_bounds, set_topol_bounds_with_outputs};
    let mol = input.topology;
    let n = mol.atoms.len();
    if params.use_exp_torsion_angle_prefs || params.use_basic_knowledge {
        set_topol_bounds_with_outputs(
            mol,
            input.rings,
            input.valence,
            input.hybridizations,
            input.conjugated,
            mmat,
            &mut details.bonds,
            &mut details.angles,
            true,
            false,
            params.use_macrocycle14config,
            params.force_trans_amides,
            true,
            true,
        )
        .map_err(GenerationFailure::from)?;
    } else {
        set_topol_bounds(
            mol,
            input.rings,
            input.valence,
            input.hybridizations,
            input.conjugated,
            mmat,
            true,
            false,
            params.use_macrocycle14config,
            params.force_trans_amides,
            true,
            true,
        )
        .map_err(GenerationFailure::from)?;
    }
    let mut tol = 0.0;
    if let Some(map) = coord_map {
        embedder_adjust_bounds_mat_from_coord_map(mmat, n, map)?;
        tol = 0.05;
    }
    if !crate::smoothing::triangle_smooth_bounds_shared(mmat, tol) {
        init_bounds_mat(mmat, 0.0, 1000.0).map_err(GenerationFailure::from)?;
        set_topol_bounds(
            mol,
            input.rings,
            input.valence,
            input.hybridizations,
            input.conjugated,
            mmat,
            false,
            true,
            params.use_macrocycle14config,
            params.force_trans_amides,
            true,
            true,
        )
        .map_err(GenerationFailure::from)?;
        if let Some(map) = coord_map {
            embedder_adjust_bounds_mat_from_coord_map(mmat, n, map)?;
        }
        if !crate::smoothing::triangle_smooth_bounds_shared(mmat, tol) {
            if params.ignore_smoothing_failures {
                init_bounds_mat(mmat, 0.0, 1000.0).map_err(GenerationFailure::from)?;
                set_topol_bounds(
                    mol,
                    input.rings,
                    input.valence,
                    input.hybridizations,
                    input.conjugated,
                    mmat,
                    false,
                    true,
                    params.use_macrocycle14config,
                    params.force_trans_amides,
                    true,
                    true,
                )
                .map_err(GenerationFailure::from)?;
                if let Some(map) = coord_map {
                    embedder_adjust_bounds_mat_from_coord_map(mmat, n, map)?;
                }
            } else {
                diagnostics.push(PreparationDiagnostic::BoundsSmoothingFailed);
                return Ok(false);
            }
        }
    }
    Ok(true)
}

#[cfg(test)]
mod original_complete_preparation_conditions {
    use super::*;
    use cosmolkit_model::{Atom, AtomSpec, Bond, BondSpec, Element};
    const DEFAULT_UPPER: f64 = 1000.0;
    fn fixed_atom(atoms: &mut Vec<Atom>, spec: AtomSpec) -> AtomId {
        let id = AtomId::new(atoms.len());
        atoms.push(Atom::from_spec(id, spec));
        id
    }
    fn fixed_bond(
        bonds: &mut Vec<Bond>,
        spec: BondSpec,
    ) -> Result<BondId, cosmolkit_model::BondValueError> {
        spec.validate()?;
        let id = BondId::new(bonds.len());
        bonds.push(Bond::from_spec(id, spec));
        Ok(id)
    }
    fn fixed_smiles(text: &str) -> TopologyBlock {
        let record = cosmolkit_smiles::parse_smiles(text, &Default::default())
            .expect("original SMILES input");
        let topology = cosmolkit_core::sanitize_topology(&record.topology, &Default::default())
            .expect("original sanitization")
            .topology;
        cosmolkit_core::assign_double_bond_stereo_from_directions(topology)
            .expect("original directional stereo")
    }
    fn fixed_chemical_rows(
        mol: &TopologyBlock,
    ) -> (RingInfo, ValenceAssignment, Vec<Hybridization>, Vec<bool>) {
        let rings = cosmolkit_core::symmetrize_sssr_with_options_from_parts(
            mol.atoms.len(),
            &mol.bonds,
            &mol.adjacency,
            false,
            false,
        )
        .expect("original rings");
        let valence = cosmolkit_core::assign_valence_for_topology(
            mol,
            cosmolkit_core::ValenceModel::RdkitLike,
        )
        .expect("original valence");
        let conjugated =
            cosmolkit_core::assign_conjugation_flags(mol, &valence).expect("unique conjugation");
        let hybrid =
            cosmolkit_core::assign_hybridization_with_conjugation(mol, &valence, &conjugated)
                .expect("unique hybridization")
                .values;
        (rings, valence, hybrid, conjugated)
    }
    fn embedder_find_double_bonds(
        mol: &TopologyBlock,
        ends: &mut Vec<(usize, usize, usize)>,
        stereo: &mut Vec<(Vec<usize>, i32)>,
        map: Option<&BTreeMap<i32, [f64; 3]>>,
    ) {
        super::embedder_find_double_bonds(mol, ends, stereo, map)
            .expect("original checked double bonds");
    }
    fn embedder_find_chiral_sets(
        mol: &TopologyBlock,
        chiral: &mut Vec<ChiralSetPtr>,
        tetra: &mut Vec<ChiralSetPtr>,
        map: Option<&BTreeMap<i32, [f64; 3]>>,
    ) {
        // Original builder fixtures have no initialized derived-ring cache. Keep that
        // explicit native fixture state; it is not RDKit uninitialized-ring parity.
        super::embedder_find_chiral_sets(mol, None, chiral, tetra, map, &mut vec![])
            .expect("original checked chiral sets");
    }
    fn embedder_init_etkdg(
        mol: &TopologyBlock,
        params: &EmbedParams,
        details: &mut CrystalFFDetails,
    ) -> Result<(), GenerationError> {
        let (rings, valence, hybrid, conjugated) = fixed_chemical_rows(mol);
        super::embedder_init_etkdg(
            &PreparedEmbeddingTopology {
                topology: mol,
                rings: &rings,
                valence: &valence,
                hybridizations: &hybrid,
                conjugated: &conjugated,
            },
            params,
            details,
        )
    }
    fn embedder_setup_initial_bounds_matrix(
        mol: &TopologyBlock,
        bounds: &mut BoundsMatrix,
        map: Option<&BTreeMap<i32, [f64; 3]>>,
        params: &EmbedParams,
        details: &mut CrystalFFDetails,
    ) -> Result<bool, GenerationError> {
        let (rings, valence, hybrid, conjugated) = fixed_chemical_rows(mol);
        super::embedder_setup_initial_bounds_matrix(
            &PreparedEmbeddingTopology {
                topology: mol,
                rings: &rings,
                valence: &valence,
                hybridizations: &hybrid,
                conjugated: &conjugated,
            },
            bounds,
            map,
            params,
            details,
            &mut vec![],
        )
    }
    #[test]
    fn embedder_find_double_bonds_collects_ends_and_clears_outputs() {
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
        let a0 = fixed_atom(&mut atoms, AtomSpec::new(Element::C));
        let a1 = fixed_atom(&mut atoms, AtomSpec::new(Element::C));
        let a2 = fixed_atom(&mut atoms, AtomSpec::new(Element::C));
        let a3 = fixed_atom(&mut atoms, AtomSpec::new(Element::C));
        fixed_bond(&mut bonds, BondSpec::new(a0, a1, BondOrder::Double)).expect("double");
        fixed_bond(&mut bonds, BondSpec::new(a0, a2, BondOrder::Single)).expect("single left");
        fixed_bond(&mut bonds, BondSpec::new(a1, a3, BondOrder::Single)).expect("single right");
        let mol = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).expect("mol");
        let mut double_bond_ends = vec![(99, 98, 97)];
        let mut stereo_double_bonds = vec![(vec![9, 8, 7, 6], -1)];

        embedder_find_double_bonds(&mol, &mut double_bond_ends, &mut stereo_double_bonds, None);

        assert_eq!(double_bond_ends, vec![(2, 0, 1), (3, 1, 0)]);
        assert!(stereo_double_bonds.is_empty());
    }

    #[test]
    fn embedder_find_double_bonds_skips_non_single_neighbor_only_for_degree_two_atom() {
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
        let a0 = fixed_atom(&mut atoms, AtomSpec::new(Element::C));
        let a1 = fixed_atom(&mut atoms, AtomSpec::new(Element::C));
        let a2 = fixed_atom(&mut atoms, AtomSpec::new(Element::C));
        let a3 = fixed_atom(&mut atoms, AtomSpec::new(Element::C));
        let a4 = fixed_atom(&mut atoms, AtomSpec::new(Element::C));
        fixed_bond(&mut bonds, BondSpec::new(a0, a1, BondOrder::Double)).expect("central double");
        fixed_bond(&mut bonds, BondSpec::new(a0, a2, BondOrder::Double))
            .expect("degree two non-single left");
        fixed_bond(&mut bonds, BondSpec::new(a1, a3, BondOrder::Double))
            .expect("degree three non-single right");
        fixed_bond(&mut bonds, BondSpec::new(a1, a4, BondOrder::Single))
            .expect("degree three single right");
        let mol = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).expect("mol");
        let mut double_bond_ends = Vec::new();
        let mut stereo_double_bonds = Vec::new();

        embedder_find_double_bonds(&mol, &mut double_bond_ends, &mut stereo_double_bonds, None);

        assert!(double_bond_ends.contains(&(3, 1, 0)));
        assert!(double_bond_ends.contains(&(4, 1, 0)));
        assert!(!double_bond_ends.contains(&(2, 0, 1)));
    }

    #[test]
    fn embedder_find_double_bonds_collects_stereo_sign_and_honors_coord_map_skip() {
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
        let a0 = fixed_atom(&mut atoms, AtomSpec::new(Element::C));
        let a1 = fixed_atom(&mut atoms, AtomSpec::new(Element::C));
        let a2 = fixed_atom(&mut atoms, AtomSpec::new(Element::C));
        let a3 = fixed_atom(&mut atoms, AtomSpec::new(Element::C));
        fixed_bond(
            &mut bonds,
            BondSpec::new(a0, a1, BondOrder::Double)
                .with_stereo(BondStereo::Z)
                .with_stereo_atoms(a2, a3),
        )
        .expect("stereo double");
        fixed_bond(&mut bonds, BondSpec::new(a0, a2, BondOrder::Single)).expect("single left");
        fixed_bond(&mut bonds, BondSpec::new(a1, a3, BondOrder::Single)).expect("single right");
        let mol = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).expect("mol");
        let mut double_bond_ends = Vec::new();
        let mut stereo_double_bonds = Vec::new();

        embedder_find_double_bonds(&mol, &mut double_bond_ends, &mut stereo_double_bonds, None);
        assert_eq!(stereo_double_bonds, vec![(vec![2, 0, 1, 3], -1)]);

        let mut coord_map = BTreeMap::new();
        coord_map.insert(2, [0.0, 0.0, 0.0]);
        embedder_find_double_bonds(
            &mol,
            &mut double_bond_ends,
            &mut stereo_double_bonds,
            Some(&coord_map),
        );
        assert_eq!(stereo_double_bonds, vec![(vec![2, 0, 1, 3], -1)]);

        coord_map.insert(3, [1.0, 0.0, 0.0]);
        embedder_find_double_bonds(
            &mol,
            &mut double_bond_ends,
            &mut stereo_double_bonds,
            Some(&coord_map),
        );
        assert!(stereo_double_bonds.is_empty());
    }

    #[test]
    fn embedder_find_double_bonds_uses_positive_sign_for_trans_and_e() {
        for stereo in [BondStereo::Trans, BondStereo::E] {
            let mut atoms = Vec::new();
            let mut bonds = Vec::new();
            let a0 = fixed_atom(&mut atoms, AtomSpec::new(Element::C));
            let a1 = fixed_atom(&mut atoms, AtomSpec::new(Element::C));
            let a2 = fixed_atom(&mut atoms, AtomSpec::new(Element::C));
            let a3 = fixed_atom(&mut atoms, AtomSpec::new(Element::C));
            fixed_bond(
                &mut bonds,
                BondSpec::new(a0, a1, BondOrder::Double)
                    .with_stereo(stereo)
                    .with_stereo_atoms(a2, a3),
            )
            .expect("stereo double");
            fixed_bond(&mut bonds, BondSpec::new(a0, a2, BondOrder::Single)).expect("single left");
            fixed_bond(&mut bonds, BondSpec::new(a1, a3, BondOrder::Single)).expect("single right");
            let mol = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).expect("mol");
            let mut double_bond_ends = Vec::new();
            let mut stereo_double_bonds = Vec::new();

            embedder_find_double_bonds(&mol, &mut double_bond_ends, &mut stereo_double_bonds, None);

            assert_eq!(stereo_double_bonds, vec![(vec![2, 0, 1, 3], 1)]);
        }
    }

    #[test]
    fn embedder_find_chiral_sets_collects_tagged_tetrahedral_center_bounds() {
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
        let center = fixed_atom(
            &mut atoms,
            AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCcw),
        );
        let n1 = fixed_atom(&mut atoms, AtomSpec::new(Element::C));
        let n2 = fixed_atom(&mut atoms, AtomSpec::new(Element::C));
        let n3 = fixed_atom(&mut atoms, AtomSpec::new(Element::C));
        let n4 = fixed_atom(&mut atoms, AtomSpec::new(Element::C));
        for nbr in [n1, n2, n3, n4] {
            fixed_bond(&mut bonds, BondSpec::new(center, nbr, BondOrder::Single)).expect("bond");
        }
        let mol = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).expect("mol");
        let mut chiral_centers = Vec::new();
        let mut tetrahedral_centers = Vec::new();

        embedder_find_chiral_sets(&mol, &mut chiral_centers, &mut tetrahedral_centers, None);

        assert_eq!(chiral_centers.len(), 1);
        assert!(tetrahedral_centers.is_empty());
        let cset = &chiral_centers[0];
        assert_eq!(
            (cset.idx0, cset.idx1, cset.idx2, cset.idx3, cset.idx4),
            (0, 1, 2, 3, 4)
        );
        assert_eq!(cset.volume_lower_bound, 5.0);
        assert_eq!(cset.volume_upper_bound, 100.0);
    }

    #[test]
    fn embedder_find_chiral_sets_uses_center_as_fourth_neighbor_for_three_coordinate_tagged_atom() {
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
        let center = fixed_atom(
            &mut atoms,
            AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCw),
        );
        let n1 = fixed_atom(&mut atoms, AtomSpec::new(Element::C));
        let n2 = fixed_atom(&mut atoms, AtomSpec::new(Element::C));
        let n3 = fixed_atom(&mut atoms, AtomSpec::new(Element::C));
        for nbr in [n1, n2, n3] {
            fixed_bond(&mut bonds, BondSpec::new(center, nbr, BondOrder::Single)).expect("bond");
        }
        let mol = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).expect("mol");
        let mut chiral_centers = Vec::new();
        let mut tetrahedral_centers = Vec::new();

        embedder_find_chiral_sets(&mol, &mut chiral_centers, &mut tetrahedral_centers, None);

        assert_eq!(chiral_centers.len(), 1);
        let cset = &chiral_centers[0];
        assert_eq!(
            (cset.idx0, cset.idx1, cset.idx2, cset.idx3, cset.idx4),
            (0, 1, 2, 3, 0)
        );
        assert_eq!(cset.volume_lower_bound, -100.0);
        assert_eq!(cset.volume_upper_bound, -2.0);
    }

    #[test]
    fn embedder_find_chiral_sets_collects_unmarked_tetrahedral_c_or_n_unless_coord_mapped() {
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
        let center = fixed_atom(&mut atoms, AtomSpec::new(Element::N));
        let ligands = [
            fixed_atom(&mut atoms, AtomSpec::new(Element::C)),
            fixed_atom(&mut atoms, AtomSpec::new(Element::C)),
            fixed_atom(&mut atoms, AtomSpec::new(Element::C)),
            fixed_atom(&mut atoms, AtomSpec::new(Element::C)),
        ];
        for nbr in ligands {
            fixed_bond(&mut bonds, BondSpec::new(center, nbr, BondOrder::Single)).expect("bond");
        }
        let mol = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).expect("mol");
        let mut chiral_centers = Vec::new();
        let mut tetrahedral_centers = Vec::new();

        embedder_find_chiral_sets(&mol, &mut chiral_centers, &mut tetrahedral_centers, None);
        assert!(chiral_centers.is_empty());
        assert_eq!(tetrahedral_centers.len(), 1);
        assert_eq!(tetrahedral_centers[0].volume_lower_bound, 0.0);
        assert_eq!(tetrahedral_centers[0].volume_upper_bound, 0.0);

        let mut coord_map = BTreeMap::new();
        coord_map.insert(center.index() as i32, [0.0, 0.0, 0.0]);
        embedder_find_chiral_sets(
            &mol,
            &mut chiral_centers,
            &mut tetrahedral_centers,
            Some(&coord_map),
        );
        assert_eq!(tetrahedral_centers.len(), 1);
    }

    #[test]
    fn embedder_find_chiral_sets_skips_hydrogen_even_when_tagged() {
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
        let center = fixed_atom(
            &mut atoms,
            AtomSpec::new(Element::H).with_chiral_tag(ChiralTag::TetrahedralCcw),
        );
        let ligands = [
            fixed_atom(&mut atoms, AtomSpec::new(Element::C)),
            fixed_atom(&mut atoms, AtomSpec::new(Element::C)),
            fixed_atom(&mut atoms, AtomSpec::new(Element::C)),
        ];
        for nbr in ligands {
            fixed_bond(&mut bonds, BondSpec::new(center, nbr, BondOrder::Single)).expect("bond");
        }
        let mol = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).expect("mol");
        let mut chiral_centers = Vec::new();
        let mut tetrahedral_centers = Vec::new();

        embedder_find_chiral_sets(&mol, &mut chiral_centers, &mut tetrahedral_centers, None);

        assert!(chiral_centers.is_empty());
        assert!(tetrahedral_centers.is_empty());
    }

    #[test]
    fn embedder_find_chiral_sets_collects_atropisomer_chiral_set_with_sorted_neighbors() {
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
        let a0 = fixed_atom(&mut atoms, AtomSpec::new(Element::C));
        let a1 = fixed_atom(&mut atoms, AtomSpec::new(Element::C));
        let left_hi = fixed_atom(&mut atoms, AtomSpec::new(Element::C));
        let left_lo = fixed_atom(&mut atoms, AtomSpec::new(Element::C));
        let right = fixed_atom(&mut atoms, AtomSpec::new(Element::C));
        fixed_bond(
            &mut bonds,
            BondSpec::new(a0, a1, BondOrder::Single)
                .with_stereo(BondStereo::AtropCcw)
                .with_stereo_atoms(left_lo, right),
        )
        .expect("atrop");
        fixed_bond(&mut bonds, BondSpec::new(a0, left_hi, BondOrder::Single)).expect("left high");
        fixed_bond(&mut bonds, BondSpec::new(a0, left_lo, BondOrder::Single)).expect("left low");
        fixed_bond(&mut bonds, BondSpec::new(a1, right, BondOrder::Single)).expect("right");
        let mol = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).expect("mol");
        let mut chiral_centers = Vec::new();
        let mut tetrahedral_centers = Vec::new();

        embedder_find_chiral_sets(&mol, &mut chiral_centers, &mut tetrahedral_centers, None);

        assert_eq!(chiral_centers.len(), 1);
        let cset = &chiral_centers[0];
        assert_eq!(
            (cset.idx0, cset.idx1, cset.idx2, cset.idx3, cset.idx4),
            (0, 2, 3, 4, 0)
        );
        assert_eq!(cset.volume_lower_bound, -100.0);
        assert_eq!(cset.volume_upper_bound, -1.0);
    }

    #[test]
    fn embedder_adjust_bounds_mat_from_coord_map_sets_exact_pair_distances() {
        let mut mmat = BoundsMatrix::new(4).expect("original checked bounds");
        let original_03_upper = mmat.get_upper(0, 3).expect("original checked bounds");
        let original_03_lower = mmat.get_lower(0, 3).expect("original checked bounds");
        let mut coord_map = BTreeMap::new();
        coord_map.insert(2, [0.0, 0.0, 0.0]);
        coord_map.insert(0, [3.0, 4.0, 0.0]);
        coord_map.insert(1, [3.0, 4.0, 12.0]);

        embedder_adjust_bounds_mat_from_coord_map(&mut mmat, 4, &coord_map)
            .expect("coord map bounds");

        assert_eq!(mmat.get_upper(0, 2).expect("original checked bounds"), 5.0);
        assert_eq!(mmat.get_lower(0, 2).expect("original checked bounds"), 5.0);
        assert_eq!(mmat.get_upper(0, 1).expect("original checked bounds"), 12.0);
        assert_eq!(mmat.get_lower(0, 1).expect("original checked bounds"), 12.0);
        assert_eq!(mmat.get_upper(1, 2).expect("original checked bounds"), 13.0);
        assert_eq!(mmat.get_lower(1, 2).expect("original checked bounds"), 13.0);
        assert_eq!(
            mmat.get_upper(0, 3).expect("original checked bounds"),
            original_03_upper
        );
        assert_eq!(
            mmat.get_lower(0, 3).expect("original checked bounds"),
            original_03_lower
        );
    }

    #[test]
    fn embedder_adjust_bounds_mat_from_coord_map_accepts_empty_and_singleton_maps() {
        let mut mmat = BoundsMatrix::new(2).expect("original checked bounds");
        let before = (
            mmat.get_upper(0, 1).expect("original checked bounds"),
            mmat.get_lower(0, 1).expect("original checked bounds"),
        );
        let empty = BTreeMap::new();
        embedder_adjust_bounds_mat_from_coord_map(&mut mmat, 2, &empty)
            .expect("empty coord map bounds");
        assert_eq!(
            (
                mmat.get_upper(0, 1).expect("original checked bounds"),
                mmat.get_lower(0, 1).expect("original checked bounds")
            ),
            before
        );

        let mut singleton = BTreeMap::new();
        singleton.insert(0, [1.0, 2.0, 3.0]);
        embedder_adjust_bounds_mat_from_coord_map(&mut mmat, 2, &singleton)
            .expect("singleton coord map bounds");
        assert_eq!(
            (
                mmat.get_upper(0, 1).expect("original checked bounds"),
                mmat.get_lower(0, 1).expect("original checked bounds")
            ),
            before
        );
    }

    #[test]
    fn embedder_init_etkdg_without_knowledge_sets_only_force_scaling() {
        let mol = Ok::<_, ()>(fixed_smiles("CCO")).expect("mol");
        let mut details = CrystalFFDetails {
            atom_nums: vec![99],
            bounds_mat_force_scaling: 0.25,
            ..CrystalFFDetails::default()
        };
        let params = EmbedParams {
            bounds_mat_force_scaling: 2.5,
            ..EmbedParams::default()
        };

        embedder_init_etkdg(&mol, &params, &mut details).expect("init");

        assert_eq!(details.atom_nums, vec![99]);
        assert_eq!(details.bounds_mat_force_scaling, 2.5);
    }

    #[test]
    fn embedder_init_etkdg_with_basic_knowledge_populates_atom_numbers() {
        let mol = Ok::<_, ()>(fixed_smiles("CCO")).expect("mol");
        let mut details = CrystalFFDetails::default();
        let params = EmbedParams {
            use_basic_knowledge: true,
            bounds_mat_force_scaling: 3.0,
            ..EmbedParams::default()
        };

        embedder_init_etkdg(&mol, &params, &mut details).expect("init");

        assert_eq!(details.atom_nums, vec![6, 6, 8]);
        assert_eq!(details.bounds_mat_force_scaling, 3.0);
    }

    #[test]
    fn embedder_init_etkdg_propagates_empty_molecule_failure() {
        let mol = TopologyBlock::try_from_parts(vec![], vec![], vec![], vec![]).expect("empty");
        let mut details = CrystalFFDetails::default();
        let params = EmbedParams {
            use_exp_torsion_angle_prefs: true,
            ..EmbedParams::default()
        };

        let err = embedder_init_etkdg(&mol, &params, &mut details).expect_err("empty error");

        assert!(matches!(
            err,
            GenerationError::GenerationFailed(message) if message.to_string() == "RDKit CrystalFF molecule has no atoms"
        ));
    }

    #[test]
    fn embedder_setup_initial_bounds_matrix_sets_topological_bounds_and_smooths() {
        let mol = Ok::<_, ()>(fixed_smiles("CCO")).expect("mol");
        let mut mmat = BoundsMatrix::new(mol.atoms.len()).expect("original checked bounds");
        let mut details = CrystalFFDetails::default();
        let params = EmbedParams::default();

        let ok = embedder_setup_initial_bounds_matrix(&mol, &mut mmat, None, &params, &mut details)
            .expect("setup");

        assert!(ok);
        assert!(mmat.check_valid());
        assert!(details.bonds.is_empty());
        assert!(details.angles.is_empty());
        assert!(mmat.get_upper(0, 1).expect("original checked bounds") < DEFAULT_UPPER);
    }

    #[test]
    fn embedder_setup_initial_bounds_matrix_records_etkdg_bonds_and_angles() {
        let mol = Ok::<_, ()>(fixed_smiles("CCO")).expect("mol");
        let mut mmat = BoundsMatrix::new(mol.atoms.len()).expect("original checked bounds");
        let mut details = CrystalFFDetails::default();
        let params = EmbedParams {
            use_basic_knowledge: true,
            ..EmbedParams::default()
        };

        let ok = embedder_setup_initial_bounds_matrix(&mol, &mut mmat, None, &params, &mut details)
            .expect("setup");

        assert!(ok);
        assert_eq!(details.bonds.len(), mol.bonds.len());
        assert!(!details.angles.is_empty());
    }

    #[test]
    fn embedder_setup_initial_bounds_matrix_applies_coord_map_exact_bounds() {
        let mol = Ok::<_, ()>(fixed_smiles("CCO")).expect("mol");
        let mut mmat = BoundsMatrix::new(mol.atoms.len()).expect("original checked bounds");
        let mut details = CrystalFFDetails::default();
        let mut coord_map = BTreeMap::new();
        coord_map.insert(0, [0.0, 0.0, 0.0]);
        coord_map.insert(1, [1.0, 0.0, 0.0]);
        let params = EmbedParams::default();

        let ok = embedder_setup_initial_bounds_matrix(
            &mol,
            &mut mmat,
            Some(&coord_map),
            &params,
            &mut details,
        )
        .expect("setup");

        assert!(ok);
        assert_eq!(mmat.get_upper(0, 1).expect("original checked bounds"), 1.0);
        assert_eq!(mmat.get_lower(0, 1).expect("original checked bounds"), 1.0);
    }

    #[test]
    fn embedder_setup_initial_bounds_matrix_propagates_topology_errors() {
        let mol = TopologyBlock::try_from_parts(vec![], vec![], vec![], vec![]).expect("empty");
        let mut mmat = BoundsMatrix::new(0).expect("original checked bounds");
        let mut details = CrystalFFDetails::default();
        let params = EmbedParams::default();

        let err =
            embedder_setup_initial_bounds_matrix(&mol, &mut mmat, None, &params, &mut details)
                .expect_err("empty molecule");

        assert!(matches!(
            err,
            GenerationError::GenerationFailed(message) if message.to_string() == "molecule has no atoms"
        ));
    }

    #[test]
    fn source_two_failed_smoothing_attempts_return_false_with_diagnostic() {
        let mol = fixed_smiles("CCO");
        let (rings, valence, hybrid, conjugated) = fixed_chemical_rows(&mol);
        let input = PreparedEmbeddingTopology {
            topology: &mol,
            rings: &rings,
            valence: &valence,
            hybridizations: &hybrid,
            conjugated: &conjugated,
        };
        // Coincident terminal atoms conflict with the different C-C/C-O rest
        // lengths even after the source second relaxed-matrix attempt.
        let map = BTreeMap::from([(0, [0.0, 0.0, 0.0]), (2, [0.0, 0.0, 0.0])]);
        let mut bounds = BoundsMatrix::new(3).unwrap();
        let mut diagnostics = vec![];
        let ok = super::embedder_setup_initial_bounds_matrix(
            &input,
            &mut bounds,
            Some(&map),
            &EmbedParams::default(),
            &mut CrystalFFDetails::default(),
            &mut diagnostics,
        )
        .unwrap();
        assert!(!ok);
        assert_eq!(
            diagnostics,
            vec![PreparationDiagnostic::BoundsSmoothingFailed]
        );
        assert_eq!(
            diagnostics[0].to_string(),
            "Could not triangle bounds smooth molecule."
        );
        assert!(!bounds.check_valid());
    }
    #[test]
    fn source_ignore_smoothing_returns_final_unsmoothed_relaxed_bounds() {
        let mol = fixed_smiles("CCO");
        let (rings, valence, hybrid, conjugated) = fixed_chemical_rows(&mol);
        let input = PreparedEmbeddingTopology {
            topology: &mol,
            rings: &rings,
            valence: &valence,
            hybridizations: &hybrid,
            conjugated: &conjugated,
        };
        let map = BTreeMap::from([(0, [0.0, 0.0, 0.0]), (2, [0.0, 0.0, 0.0])]);
        let mut bounds = BoundsMatrix::new(3).unwrap();
        let mut diagnostics = vec![];
        let params = EmbedParams {
            ignore_smoothing_failures: true,
            ..Default::default()
        };
        let ok = super::embedder_setup_initial_bounds_matrix(
            &input,
            &mut bounds,
            Some(&map),
            &params,
            &mut CrystalFFDetails::default(),
            &mut diagnostics,
        )
        .unwrap();
        assert!(ok);
        assert!(diagnostics.is_empty());
        assert_eq!(bounds.get_lower(0, 2).unwrap(), 0.0);
        assert_eq!(bounds.get_upper(0, 2).unwrap(), 0.0);
        assert!(bounds.check_valid());
        // The final source branch deliberately does not triangle-smooth again.
        // Pair validity and triangle consistency are different guarantees.
        assert!(!crate::smoothing::triangle_smooth_bounds_shared(
            &mut bounds,
            0.05
        ));
    }
}

struct EmbedHelperArgs<'a> {
    confs_ok: &'a mut [bool],
    four_d: bool,
    frag_mapping: Option<&'a [usize]>,
    confs: &'a mut [cosmolkit_model::Conformer3D],
    frag_idx: usize,
    mmat: &'a BoundsMatrix,
    chiral_centers: &'a [ChiralSetPtr],
    tetrahedral_carbons: &'a [ChiralSetPtr],
    double_bond_ends: &'a [(usize, usize, usize)],
    stereo_double_bonds: &'a [(Vec<usize>, i32)],
    etkdg_details: &'a CrystalFFDetails,
}
impl<'a> EmbedHelperArgs<'a> {
    fn worker_geometry(&self) -> WorkerGeometry<'a> {
        WorkerGeometry {
            four_d: self.four_d,
            frag_mapping: self.frag_mapping,
            frag_idx: self.frag_idx,
            mmat: self.mmat,
            chiral_centers: self.chiral_centers,
            tetrahedral_carbons: self.tetrahedral_carbons,
            double_bond_ends: self.double_bond_ends,
            stereo_double_bonds: self.stereo_double_bonds,
            etkdg_details: self.etkdg_details,
        }
    }
}
fn embedder_embed_helper(
    thread_id: i32,
    num_threads: i32,
    eargs: &mut EmbedHelperArgs<'_>,
    params: &mut EmbedParams,
    end_time: Option<Instant>,
) -> Result<(), GenerationError> {
    if num_threads <= 0 || thread_id < 0 || thread_id >= num_threads {
        return Err(GenerationError::Input(
            "invalid embedding worker index or count",
        ));
    }
    if eargs.confs_ok.len() != eargs.confs.len() {
        return Err(GenerationError::Input(
            "embedding conformer status count mismatch",
        ));
    }
    let geometry = eargs.worker_geometry();
    let total = eargs.confs.len();
    with_embedding_progress(params, |settings, progress| {
        let slots = eargs
            .confs_ok
            .iter_mut()
            .zip(eargs.confs.iter_mut())
            .enumerate()
            .filter(|(ci, _)| (*ci % num_threads as usize) as i32 == thread_id)
            .map(|(ci, (ok, conf))| (ci, ok, conf));
        embedder_worker(
            thread_id,
            num_threads,
            total,
            slots,
            geometry,
            settings,
            progress,
            end_time,
        )
    })
}
fn embedder_worker<'a>(
    thread_id: i32,
    num_threads: i32,
    total: usize,
    mut slots: impl Iterator<Item = WorkerSlot<'a>>,
    geometry: WorkerGeometry<'_>,
    params: &EmbedParams,
    progress: &EmbeddingProgress,
    end_time: Option<Instant>,
) -> Result<(), GenerationError> {
    // BEGIN RDKIT CPP FUNCTION embedder_embed_helper (GraphMol/DistGeomHelpers/Embedder.cpp:1356-1447)
    // RDKit❗❌: void embedHelper_(int threadId, int numThreads, EmbedArgs *eargs,
    // RDKit❗❌:                   EmbedParameters *params, TimePoint *end_time) {
    // RDKit❗❌:   PRECONDITION(eargs, "bogus eargs");
    // RDKit❗❌:   PRECONDITION(params, "bogus params");
    // RDKit❗❌:   unsigned int nAtoms = eargs->mmat->numRows();
    // RDKit❗❌:   RDGeom::PointPtrVect positions(nAtoms);
    // RDKit❗❌:   // we might thrown an exception in a callback
    // RDKit❗❌:   // in order to avoid leaking the points we're working with
    // RDKit❗❌:   // allocate them with unique_ptrs and then work with the naked
    // RDKit❗❌:   // pointers from those
    // RDKit❗❌:   std::vector<std::unique_ptr<RDGeom::Point>> positionsStore;
    // RDKit❗❌:   positionsStore.reserve(nAtoms);
    // RDKit❗❌:   for (unsigned int i = 0; i < nAtoms; ++i) {
    // RDKit❗❌:     if (eargs->fourD) {
    // RDKit❗❌:       positionsStore.emplace_back(new RDGeom::PointND(4));
    // RDKit❗❌:     } else {
    // RDKit❗❌:       positionsStore.emplace_back(new RDGeom::Point3D());
    // RDKit❗❌:     }
    // RDKit❗❌:     positions[i] = positionsStore[i].get();
    // RDKit❗❌:   }
    // RDKit❗❌:   for (size_t ci = 0; ci < eargs->confs->size(); ci++) {
    // RDKit❗❌:     if (ControlCHandler::getGotSignal() ||
    // RDKit❗❌:         (end_time != nullptr && Clock::now() > *end_time)) {
    // RDKit❗❌:       return;
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     if (rdcast<int>(ci % numThreads) != threadId) {
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:     if (!(*eargs->confsOk)[ci]) {
    // RDKit❗❌:       // we call this function for each fragment in a molecule,
    // RDKit❗❌:       // if one of the fragments has already failed, there's no
    // RDKit❗❌:       // sense in embedding this one
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     CHECK_INVARIANT(
    // RDKit❗❌:         params->randomSeed >= -1,
    // RDKit❗❌:         "random seed must either be positive, zero, or negative one");
    // RDKit❗❌:     int new_seed = params->randomSeed;
    // RDKit❗❌:     if (new_seed > -1) {
    // RDKit❗❌:       if (params->enableSequentialRandomSeeds) {
    // RDKit❗❌:         new_seed += ci + 1;
    // RDKit❗❌:       } else {
    // RDKit❗❌:         if (!multiplication_overflows_(rdcast<int>(ci + 1),
    // RDKit❗❌:                                        params->randomSeed)) {
    // RDKit❗❌:           // old method of computing a new seed
    // RDKit❗❌:           new_seed = (ci + 1) * params->randomSeed;
    // RDKit❗❌:         } else {
    // RDKit❗❌:           // If the above simple multiplication will overflow, use a
    // RDKit❗❌:           // cheap and easy way to hash the conformer index and seed
    // RDKit❗❌:           // together: for N'ary numerical system, where N is the
    // RDKit❗❌:           // maximum possible value of the pair of numbers. The
    // RDKit❗❌:           // following will generate unique integers:
    // RDKit❗❌:           // hash(a, b) = a + b * N
    // RDKit❗❌:           auto big_seed = rdcast<size_t>(params->randomSeed);
    // RDKit❗❌:           size_t max_val = std::max(ci + 1, big_seed);
    // RDKit❗❌:           size_t big_num = big_seed + max_val * (ci + 1);
    // RDKit❗❌:           // only grab the first 31 bits xor'd with the next 31 bits to
    // RDKit❗❌:           // make sure its positive, careful, the 'ULL' is important
    // RDKit❗❌:           // here, 0x7fffffff is the 'int' type because of C default
    // RDKit❗❌:           // number semantics and that we definitely don't want!
    // RDKit❗❌:           const size_t positive_int_mask = 0x7fffffffULL;
    // RDKit❗❌:           size_t folded_num =
    // RDKit❗❌:               (big_num & positive_int_mask) ^ (big_num >> 31ULL);
    // RDKit❗❌:           new_seed = rdcast<int>(folded_num & positive_int_mask);
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     CHECK_INVARIANT(new_seed >= -1,
    // RDKit❗❌:                     "Something went wrong calculating a new seed");
    // RDKit❗❌:     bool gotCoords = EmbeddingOps::embedPoints(&positions, *eargs, *params,
    // RDKit❗❌:                                                new_seed, end_time);
    // RDKit❗❌:
    // RDKit❗❌:     // copy the coordinates into the correct conformer
    // RDKit❗❌:     if (gotCoords) {
    // RDKit❗❌:       auto &conf = (*eargs->confs)[ci];
    // RDKit❗❌:       unsigned int fragAtomIdx = 0;
    // RDKit❗❌:       for (unsigned int i = 0; i < conf->getNumAtoms(); ++i) {
    // RDKit❗❌:         if (!eargs->fragMapping ||
    // RDKit❗❌:             (*eargs->fragMapping)[i] == static_cast<int>(eargs->fragIdx)) {
    // RDKit❗❌:           conf->setAtomPos(i, RDGeom::Point3D((*positions[fragAtomIdx])[0],
    // RDKit❗❌:                                               (*positions[fragAtomIdx])[1],
    // RDKit❗❌:                                               (*positions[fragAtomIdx])[2]));
    // RDKit❗❌:           ++fragAtomIdx;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     } else {
    // RDKit❗❌:       (*eargs->confsOk)[ci] = 0;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION embedder_embed_helper

    let n = geometry.mmat.dimension();
    let dim = if geometry.four_d { 4 } else { 3 };
    // One reusable owned point vector per worker; no conformer clones and no
    // point allocation inside the source per-conformer loop.
    let mut positions: Vec<Vec<f64>> = (0..n).map(|_| vec![0.0; dim]).collect();
    for ci in 0..total {
        if source_got_signal() || end_time.is_some_and(|t| Instant::now() > t) {
            return Ok(());
        }
        if (ci % num_threads as usize) as i32 != thread_id {
            continue;
        }
        let (slot_ci, ok, conf) = slots
            .next()
            .ok_or(GenerationError::Input("embedding worker slot missing"))?;
        if slot_ci != ci {
            return Err(GenerationError::Input(
                "embedding worker slot index mismatch",
            ));
        }
        if !*ok {
            continue;
        }
        let new_seed = crate::numeric::rdkit_embedder_conformer_seed(
            params.random_seed,
            ci,
            params.enable_sequential_random_seeds,
        )?;
        // Validate only the selected result's destination shape. Source callers
        // already guarantee these lengths; this O(N) pass is an extra safe
        // detached boundary cost compared with direct C++ pointer placement.
        let rows = conf.coordinates().len();
        let selected = if let Some(map) = geometry.frag_mapping {
            if map.len() != rows {
                return Err(GenerationError::Input(
                    "embedding fragment mapping count mismatch",
                ));
            }
            map.iter()
                .filter(|&&frag| frag == geometry.frag_idx)
                .count()
        } else {
            rows
        };
        if selected != n {
            return Err(GenerationError::Input(
                "embedding fragment point count mismatch",
            ));
        }
        let got = embedder_embed_points_shared(
            &mut positions,
            geometry.embed_args(),
            params,
            progress,
            new_seed,
            end_time,
        )?;
        if got {
            let mut frag_atom = 0;
            for (i, p) in conf.coordinates_mut().iter_mut().enumerate() {
                if geometry
                    .frag_mapping
                    .is_none_or(|map| map[i] == geometry.frag_idx)
                {
                    *p = [
                        positions[frag_atom][0],
                        positions[frag_atom][1],
                        positions[frag_atom][2],
                    ];
                    frag_atom += 1;
                }
            }
        } else {
            *ok = false;
        }
    }
    Ok(())
}

#[cfg(test)]
mod original_complete_helper_conditions {
    use super::*;
    use crate::numeric::rdkit_embedder_multiplication_overflows;
    use cosmolkit_model::Conformer3D;
    fn random_embed_params(seed: i32, num_threads: i32) -> EmbedParams {
        EmbedParams {
            max_iterations: 1,
            random_seed: seed,
            use_random_coords: true,
            num_threads,
            box_size_mult: -2.0,
            use_symmetry_for_pruning: false,
            symmetrize_conjugated_terminal_groups_for_pruning: false,
            ..Default::default()
        }
    }
    #[test]
    fn embedder_multiplication_overflows_matches_rdkit_boundaries() {
        assert!(!rdkit_embedder_multiplication_overflows(0, i32::MAX));
        assert!(!rdkit_embedder_multiplication_overflows(1, i32::MAX));
        assert!(!rdkit_embedder_multiplication_overflows(46_340, 46_341));
        assert!(rdkit_embedder_multiplication_overflows(46_341, 46_341));
        assert!(rdkit_embedder_multiplication_overflows(46_342, 46_341));
    }

    #[test]
    fn embedder_helper_uses_thread_scheduling_and_seed_policy() {
        let mut mmat = BoundsMatrix::new(1).expect("original checked bounds");
        let mut confs = vec![
            Conformer3D::new(0, vec![[0.0, 0.0, 0.0]], true),
            Conformer3D::new(1, vec![[0.0, 0.0, 0.0]], true),
        ];
        let mut confs_ok = vec![true, true];
        let details = CrystalFFDetails::default();
        let mut params = random_embed_params(11, 2);
        let mut args = EmbedHelperArgs {
            confs_ok: &mut confs_ok,
            four_d: true,
            frag_mapping: None,
            confs: &mut confs,
            frag_idx: 0,
            mmat: &mmat,
            chiral_centers: &[],
            tetrahedral_carbons: &[],
            double_bond_ends: &[],
            stereo_double_bonds: &[],
            etkdg_details: &details,
        };

        embedder_embed_helper(1, 2, &mut args, &mut params, None).expect("embed helper");

        assert_eq!(args.confs[0].coordinates()[0], [0.0, 0.0, 0.0]);
        assert_ne!(args.confs[1].coordinates()[0], [0.0, 0.0, 0.0]);
    }

    #[test]
    fn source_previously_failed_fragment_is_skipped_without_mutation() {
        let bounds = BoundsMatrix::new(1).unwrap();
        let details = CrystalFFDetails::default();
        let mut confs = [Conformer3D::new(0, vec![[9.0, 8.0, 7.0]], true)];
        let mut ok = [false];
        let mut params = random_embed_params(11, 1);
        let mut args = EmbedHelperArgs {
            confs_ok: &mut ok,
            four_d: true,
            frag_mapping: None,
            confs: &mut confs,
            frag_idx: 0,
            mmat: &bounds,
            chiral_centers: &[],
            tetrahedral_carbons: &[],
            double_bond_ends: &[],
            stereo_double_bonds: &[],
            etkdg_details: &details,
        };
        embedder_embed_helper(0, 1, &mut args, &mut params, None).unwrap();
        assert_eq!(args.confs[0].coordinates(), &[[9.0, 8.0, 7.0]]);
        assert_eq!(args.confs_ok, &[false]);
        assert_eq!(params.basin_thresh, 5.0);
    }
    #[test]
    fn source_helper_places_only_selected_fragment_and_preserves_other_rows() {
        let bounds = BoundsMatrix::new(1).unwrap();
        let details = CrystalFFDetails::default();
        let mut confs = [Conformer3D::new(
            0,
            vec![[9.0, 8.0, 7.0], [0.0; 3], [6.0, 5.0, 4.0]],
            true,
        )];
        let mut ok = [true];
        let mut params = random_embed_params(11, 1);
        params.coord_map = Some(BTreeMap::from([(0, [1.0, 2.0, 3.0])]));
        let mut args = EmbedHelperArgs {
            confs_ok: &mut ok,
            four_d: true,
            frag_mapping: Some(&[0, 1, 0]),
            confs: &mut confs,
            frag_idx: 1,
            mmat: &bounds,
            chiral_centers: &[],
            tetrahedral_carbons: &[],
            double_bond_ends: &[],
            stereo_double_bonds: &[],
            etkdg_details: &details,
        };
        embedder_embed_helper(0, 1, &mut args, &mut params, None).unwrap();
        assert_eq!(
            args.confs[0].coordinates(),
            &[[9.0, 8.0, 7.0], [1.0, 2.0, 3.0], [6.0, 5.0, 4.0]]
        );
        assert_eq!(args.confs_ok, &[true]);
    }
    #[test]
    fn source_helper_failed_embedding_marks_status_and_retains_coordinates_and_cause() {
        let mut bounds = BoundsMatrix::new(3).unwrap();
        for (i, j, d) in [(0, 1, 1.0), (1, 2, 1.0), (0, 2, 2.0)] {
            bounds.set_lower(i, j, d).unwrap();
            bounds.set_upper(i, j, d).unwrap();
        }
        let details = CrystalFFDetails::default();
        let mut confs = [Conformer3D::new(0, vec![[9.0, 8.0, 7.0]; 3], true)];
        let mut ok = [true];
        let mut params = random_embed_params(11, 1);
        params.track_failures = true;
        params.failures = vec![0; crate::EmbedFailureCause::EndOfEnum as usize];
        params.coord_map = Some(BTreeMap::from([
            (0, [0.0, 0.0, 0.0]),
            (1, [1.0, 0.0, 0.0]),
            (2, [2.0, 0.0, 0.0]),
        ]));
        let mut args = EmbedHelperArgs {
            confs_ok: &mut ok,
            four_d: true,
            frag_mapping: None,
            confs: &mut confs,
            frag_idx: 0,
            mmat: &bounds,
            chiral_centers: &[],
            tetrahedral_carbons: &[],
            double_bond_ends: &[(0, 1, 2)],
            stereo_double_bonds: &[],
            etkdg_details: &details,
        };
        embedder_embed_helper(0, 1, &mut args, &mut params, None).unwrap();
        assert_eq!(args.confs[0].coordinates(), &[[9.0, 8.0, 7.0]; 3]);
        assert_eq!(args.confs_ok, &[false]);
        assert_eq!(
            params.failures[crate::EmbedFailureCause::LinearDoubleBond as usize],
            1
        );
        assert_eq!(params.failures.iter().sum::<u32>(), 1);
    }
    #[test]
    fn source_helper_deadline_precedes_scheduling_and_embedding() {
        let bounds = BoundsMatrix::new(1).unwrap();
        let details = CrystalFFDetails::default();
        let mut confs = [Conformer3D::new(0, vec![[9.0, 8.0, 7.0]], true)];
        let mut ok = [true];
        let mut params = random_embed_params(11, 1);
        let mut args = EmbedHelperArgs {
            confs_ok: &mut ok,
            four_d: true,
            frag_mapping: None,
            confs: &mut confs,
            frag_idx: 0,
            mmat: &bounds,
            chiral_centers: &[],
            tetrahedral_carbons: &[],
            double_bond_ends: &[],
            stereo_double_bonds: &[],
            etkdg_details: &details,
        };
        embedder_embed_helper(
            0,
            1,
            &mut args,
            &mut params,
            Some(Instant::now() - std::time::Duration::from_secs(1)),
        )
        .unwrap();
        assert_eq!(args.confs[0].coordinates(), &[[9.0, 8.0, 7.0]]);
        assert_eq!(args.confs_ok, &[true]);
        assert_eq!(params.basin_thresh, 5.0);
    }
}

fn embedder_fill_atom_positions(
    points: &mut [[f64; 3]],
    conf: &cosmolkit_model::Conformer3D,
    matched: &[usize],
) -> Result<(), GenerationError> {
    // BEGIN RDKIT CPP FUNCTION _fillAtomPositions (Embedder.cpp:1309-1315)
    // RDKit❗❌: void _fillAtomPositions(RDGeom::Point3DConstPtrVect &pts, const Conformer &conf,
    // RDKit❗❌:                         const ROMol &, const std::vector<unsigned int> &match) {
    // RDKit❗❌:   PRECONDITION(pts.size() == match.size(), "bad pts size");
    // RDKit❗❌:   for (unsigned int i = 0; i < match.size(); i++) {
    // RDKit❗❌:     pts[i] = &conf.getAtomPos(match[i]);
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION _fillAtomPositions (Embedder.cpp:1309-1315)

    if points.len() != matched.len() {
        return Err(GenerationError::Input("bad pts size"));
    }
    for (point, &atom) in points.iter_mut().zip(matched) {
        *point = *conf
            .coordinates()
            .get(atom)
            .ok_or(GenerationError::Input("atom index out of range"))?;
    }
    Ok(())
}
fn embedder_is_conf_far_from_rest<'a>(
    existing: impl Clone + IntoIterator<Item = &'a cosmolkit_model::Conformer3D>,
    conf: &cosmolkit_model::Conformer3D,
    threshold: f64,
    self_matches: &[Vec<usize>],
) -> Result<bool, GenerationError> {
    // BEGIN RDKIT CPP FUNCTION _isConfFarFromRest (Embedder.cpp:1317-1343)
    // RDKit❗❌: bool _isConfFarFromRest(
    // RDKit❗❌:     const ROMol &mol, const Conformer &conf, double threshold,
    // RDKit❗❌:     const std::vector<std::vector<unsigned int>> &selfMatches) {
    // RDKit❗❌:   // NOTE: it is tempting to use some triangle inequality to prune
    // RDKit❗❌:   // conformations here but some basic testing has shown very
    // RDKit❗❌:   // little advantage and given that the time for pruning fades in
    // RDKit❗❌:   // comparison to embedding - we will use a simple for loop below
    // RDKit❗❌:   // over all conformation until we find a match
    // RDKit❗❌:   RDGeom::Point3DConstPtrVect refPoints(selfMatches[0].size());
    // RDKit❗❌:   RDGeom::Point3DConstPtrVect prbPoints(selfMatches[0].size());
    // RDKit❗❌:   _fillAtomPositions(refPoints, conf, mol, selfMatches[0]);
    // RDKit❗❌:
    // RDKit❗❌:   double ssrThres = selfMatches[0].size() * threshold * threshold;
    // RDKit❗❌:   for (const auto &match : selfMatches) {
    // RDKit❗❌:     for (auto confi = mol.beginConformers(); confi != mol.endConformers();
    // RDKit❗❌:          ++confi) {
    // RDKit❗❌:       _fillAtomPositions(prbPoints, *(*confi), mol, match);
    // RDKit❗❌:       RDGeom::Transform3D trans;
    // RDKit❗❌:       auto ssr =
    // RDKit❗❌:           RDNumeric::Alignments::AlignPoints(refPoints, prbPoints, trans);
    // RDKit❗❌:       if (ssr < ssrThres) {
    // RDKit❗❌:         return false;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return true;
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION _isConfFarFromRest (Embedder.cpp:1317-1343)

    // Source stores two vectors of coordinate pointers. Detached numeric inputs
    // copy three scalars per selected point into two reused contiguous buffers;
    // this adds O(N) copy traffic and triples buffer storage, marked ❗❌.
    let first = self_matches
        .first()
        .ok_or(GenerationError::Input("self matches must not be empty"))?;
    let mut reference = vec![[0.0; 3]; first.len()];
    let mut probe = vec![[0.0; 3]; first.len()];
    embedder_fill_atom_positions(&mut reference, conf, first)?;
    let ssr_threshold = first.len() as f64 * threshold * threshold;
    for matched in self_matches {
        for previous in existing.clone() {
            embedder_fill_atom_positions(&mut probe, previous, matched)?;
            let ssr = cosmolkit_core::alignment_sum_squared_residual(&reference, &probe)
                .map_err(GenerationError::Input)?;
            if ssr < ssr_threshold {
                return Ok(false);
            }
        }
    }
    Ok(true)
}
#[cfg(test)]
mod original_numeric_pruning_conditions {
    use super::*;
    use cosmolkit_model::Conformer3D;
    #[test]
    fn embedder_fill_atom_positions_copies_positions_in_match_order() {
        let conf = Conformer3D::new(
            0,
            vec![[0.0, 1.0, 2.0], [3.0, 4.0, 5.0], [6.0, 7.0, 8.0]],
            true,
        );
        let mut pts = vec![[0.0; 3]; 3];
        embedder_fill_atom_positions(&mut pts, &conf, &[2, 0, 1])
            .expect("original checked source precondition");
        assert_eq!(pts[0], [6.0, 7.0, 8.0]);
        assert_eq!(pts[1], [0.0, 1.0, 2.0]);
        assert_eq!(pts[2], [3.0, 4.0, 5.0]);
    }
    #[test]
    #[should_panic(expected = "bad pts size")]
    fn embedder_fill_atom_positions_rejects_size_mismatch() {
        let conf = Conformer3D::new(0, vec![[0.0, 0.0, 0.0]; 3], true);
        let mut pts = vec![[0.0; 3]; 2];
        embedder_fill_atom_positions(&mut pts, &conf, &[0])
            .expect("original checked source precondition");
    }
    #[test]
    fn embedder_is_conf_far_from_rest_prunes_close_existing_conformer_and_keeps_far_one() {
        let existing = vec![Conformer3D::new(
            0,
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            true,
        )];
        let close = Conformer3D::new(1, vec![[0.01, 0.0, 0.0], [1.01, 0.0, 0.0]], true);
        let far = Conformer3D::new(2, vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]], true);
        let matches = vec![vec![0, 1]];
        assert!(
            !embedder_is_conf_far_from_rest(&existing, &close, 0.5, &matches)
                .expect("original checked source pruning")
        );
        assert!(
            embedder_is_conf_far_from_rest(&existing, &far, 0.5, &matches)
                .expect("original checked source pruning")
        );
    }
    #[test]
    fn source_pruning_strict_threshold_empty_selection_and_late_matches() {
        let conf = Conformer3D::new(0, vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], true);
        let previous = Conformer3D::new(1, vec![[10.0, 0.0, 0.0], [11.0, 0.0, 0.0]], true);
        assert!(
            embedder_is_conf_far_from_rest(
                std::slice::from_ref(&previous),
                &conf,
                0.0,
                &[vec![0, 1]]
            )
            .unwrap()
        );
        assert!(
            !embedder_is_conf_far_from_rest(
                std::slice::from_ref(&previous),
                &conf,
                0.5,
                &[vec![0, 1]]
            )
            .unwrap()
        );
        assert!(
            embedder_is_conf_far_from_rest(std::slice::from_ref(&previous), &conf, 0.5, &[vec![]])
                .unwrap()
        );
        assert!(embedder_is_conf_far_from_rest(&[], &conf, 0.5, &[vec![]]).unwrap());
        let three = Conformer3D::new(
            2,
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            true,
        );
        assert!(
            embedder_is_conf_far_from_rest(std::slice::from_ref(&three), &conf, 0.1, &[vec![0, 1]])
                .unwrap()
        );
        assert!(
            !embedder_is_conf_far_from_rest(
                std::slice::from_ref(&three),
                &conf,
                0.1,
                &[vec![0, 1], vec![0, 2]]
            )
            .unwrap()
        );
    }
}

// Only mutable source progress is shared. Configuration maps, geometry and
// every optimizer remain borrowed or worker-owned; no EmbedParams clone occurs.
struct EmbeddingProgress {
    failures: std::sync::Mutex<Vec<u32>>,
    defaults: std::sync::OnceLock<(u32, f64)>,
}
fn with_embedding_progress<T>(
    params: &mut EmbedParams,
    run: impl FnOnce(&EmbedParams, &EmbeddingProgress) -> T,
) -> T {
    // BEGIN RDKIT CPP BLOCK failure mutex initialization (Embedder.cpp:78-94)
    // RDKit❗❌: std::mutex &failmutex_get() {
    // RDKit❗❌:   // create on demand
    // RDKit❗❌:   static std::mutex _mutex;
    // RDKit❗❌:   return _mutex;
    // RDKit❗❌: }
    // RDKit❗❌:
    // RDKit❗❌: void failmutex_create() {
    // RDKit❗❌:   std::mutex &mutex = failmutex_get();
    // RDKit❗❌:   std::lock_guard<std::mutex> test_lock(mutex);
    // RDKit❗❌: }
    // RDKit❗❌:
    // RDKit❗❌: std::mutex &GetFailMutex() {
    // RDKit❗❌:   static std::once_flag flag;
    // RDKit❗❌:   std::call_once(flag, failmutex_create);
    // RDKit❗❌:   return failmutex_get();
    // RDKit❗❌: }
    // RDKit❗❌: }  // namespace
    // END RDKIT CPP BLOCK failure mutex initialization (Embedder.cpp:78-94)

    // Rust separates source's mutable fields from shared options. Moving the
    // failure vector is O(1). Restore progress on typed errors and callback
    // unwinding; the source does not roll these changes back on exceptions.
    let progress = EmbeddingProgress {
        failures: std::sync::Mutex::new(std::mem::take(&mut params.failures)),
        defaults: std::sync::OnceLock::new(),
    };
    let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| run(params, &progress)));
    if let Some(&(iterations, basin)) = progress.defaults.get() {
        params.max_iterations = iterations;
        params.basin_thresh = basin;
    }
    params.failures = progress
        .failures
        .into_inner()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    match outcome {
        Ok(value) => value,
        Err(payload) => std::panic::resume_unwind(payload),
    }
}

#[derive(Clone, Copy)]
struct WorkerGeometry<'a> {
    four_d: bool,
    frag_mapping: Option<&'a [usize]>,
    frag_idx: usize,
    mmat: &'a BoundsMatrix,
    chiral_centers: &'a [ChiralSetPtr],
    tetrahedral_carbons: &'a [ChiralSetPtr],
    double_bond_ends: &'a [(usize, usize, usize)],
    stereo_double_bonds: &'a [(Vec<usize>, i32)],
    etkdg_details: &'a CrystalFFDetails,
}
impl WorkerGeometry<'_> {
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
type WorkerSlot<'a> = (usize, &'a mut bool, &'a mut cosmolkit_model::Conformer3D);

fn embedder_dispatch_workers(
    args: &mut EmbedHelperArgs<'_>,
    params: &mut EmbedParams,
    end_time: Option<Instant>,
) -> Result<(), GenerationError> {
    // BEGIN RDKIT CPP BLOCK thread/reset/future dispatch (Embedder.cpp:1636-1661)
    // RDKit❗❌:     int numThreads = getNumThreadsToUse(params.numThreads);
    // RDKit❗❌:
    // RDKit❗❌:     ControlCHandler::reset();
    // RDKit❗❌:
    // RDKit❗❌:     // do the embedding, using multiple threads if requested
    // RDKit❗❌:     detail::EmbedArgs eargs = {&confsOk,        fourD,
    // RDKit❗❌:                                &fragMapping,    &confs,
    // RDKit❗❌:                                fragIdx,         mmat,
    // RDKit❗❌:                                &chiralCenters,  &tetrahedralCarbons,
    // RDKit❗❌:                                &doubleBondEnds, &stereoDoubleBonds,
    // RDKit❗❌:                                &etkdgDetails};
    // RDKit❗❌:     if (numThreads == 1) {
    // RDKit❗❌:       detail::embedHelper_(0, 1, &eargs, &params, end_time);
    // RDKit❗❌:     }
    // RDKit❗❌: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit❗❌:     else {
    // RDKit❗❌:       std::vector<std::future<void>> tg;
    // RDKit❗❌:       for (int tid = 0; tid < numThreads; ++tid) {
    // RDKit❗❌:         tg.emplace_back(std::async(std::launch::async, detail::embedHelper_,
    // RDKit❗❌:                                    tid, numThreads, &eargs, &params, end_time));
    // RDKit❗❌:       }
    // RDKit❗❌:       for (auto &fut : tg) {
    // RDKit❗❌:         fut.get();
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌: #endif
    // END RDKIT CPP BLOCK thread/reset/future dispatch (Embedder.cpp:1636-1661)

    let source_count = cosmolkit_core::rdkit_thread_count(params.num_threads)?.get() as i32;
    source_reset_interrupt()?;
    if args.confs_ok.len() != args.confs.len() {
        return Err(GenerationError::Input(
            "embedding conformer status count mismatch",
        ));
    }
    // The unsigned-to-int source assignment may produce a nonpositive count;
    // its async for-loop then launches no workers. Preserve that branch.
    if source_count <= 0 {
        return Ok(());
    }
    let geometry = args.worker_geometry();
    let total = args.confs.len();
    with_embedding_progress(params, |settings, progress| {
        if source_count == 1 {
            let slots = args
                .confs_ok
                .iter_mut()
                .zip(args.confs.iter_mut())
                .enumerate()
                .map(|(ci, (ok, conf))| (ci, ok, conf));
            return embedder_worker(0, 1, total, slots, geometry, settings, progress, end_time);
        }
        // Disjoint mutable borrows are distributed once. This O(C+T) address
        // storage is extra compared with C++ shared pointers (complexity ❌);
        // conformer rows, options, maps and geometry are never cloned.
        let mut partitions: Vec<Vec<WorkerSlot<'_>>> =
            (0..source_count).map(|_| Vec::new()).collect();
        for (ci, (ok, conf)) in args
            .confs_ok
            .iter_mut()
            .zip(args.confs.iter_mut())
            .enumerate()
        {
            partitions[ci % source_count as usize].push((ci, ok, conf));
        }
        std::thread::scope(|scope| {
            let mut handles = Vec::with_capacity(source_count as usize);
            let mut launch_error = None;
            for (tid, slots) in partitions.into_iter().enumerate() {
                match std::thread::Builder::new().spawn_scoped(scope, move || {
                    embedder_worker(
                        tid as i32,
                        source_count,
                        total,
                        slots.into_iter(),
                        geometry,
                        settings,
                        progress,
                        end_time,
                    )
                }) {
                    Ok(handle) => handles.push(handle),
                    Err(cause) => {
                        launch_error = Some(GenerationError::WorkerLaunch(cause));
                        break;
                    }
                }
            }
            // Source futures are observed in launch order. Join all borrowed
            // workers even after an error; scoped ownership cannot outlive them.
            let mut first_error = launch_error;
            let mut first_panic = None;
            for handle in handles {
                match handle.join() {
                    Ok(Ok(())) => (),
                    Ok(Err(cause)) if first_error.is_none() && first_panic.is_none() => {
                        first_error = Some(cause)
                    }
                    Err(payload) if first_error.is_none() && first_panic.is_none() => {
                        first_panic = Some(payload)
                    }
                    _ => (),
                }
            }
            if let Some(payload) = first_panic {
                std::panic::resume_unwind(payload);
            }
            match first_error {
                Some(cause) => Err(cause),
                None => Ok(()),
            }
        })
    })
}

#[cfg(all(test, not(target_arch = "wasm32")))]
mod source_parallel_worker_conditions {
    use super::*;
    use cosmolkit_model::Conformer3D;
    use std::sync::{
        Mutex,
        atomic::{AtomicUsize, Ordering},
    };
    static CALLBACK_THREADS: Mutex<Vec<std::thread::ThreadId>> = Mutex::new(Vec::new());
    static CALLBACK_ARRIVALS: AtomicUsize = AtomicUsize::new(0);
    fn overlap_callback(iter: u32) {
        assert_eq!(iter, 1);
        CALLBACK_THREADS
            .lock()
            .unwrap()
            .push(std::thread::current().id());
        CALLBACK_ARRIVALS.fetch_add(1, Ordering::SeqCst);
        let limit = Instant::now() + std::time::Duration::from_secs(5);
        while CALLBACK_ARRIVALS.load(Ordering::SeqCst) < 2 {
            assert!(Instant::now() < limit, "two source workers did not overlap");
            std::thread::yield_now();
        }
    }
    fn throwing_callback(_: u32) {
        panic!("source callback exception");
    }
    fn params(threads: i32) -> EmbedParams {
        EmbedParams {
            random_seed: 11,
            num_threads: threads,
            use_random_coords: true,
            box_size_mult: -2.0,
            ..Default::default()
        }
    }
    fn run(
        confs: &mut [Conformer3D],
        ok: &mut [bool],
        bounds: &BoundsMatrix,
        params: &mut EmbedParams,
        mapping: Option<&[usize]>,
        frag: usize,
        double_ends: &[(usize, usize, usize)],
        deadline: Option<Instant>,
    ) -> Result<(), GenerationError> {
        let details = CrystalFFDetails::default();
        let mut args = EmbedHelperArgs {
            confs_ok: ok,
            four_d: true,
            frag_mapping: mapping,
            confs,
            frag_idx: frag,
            mmat: bounds,
            chiral_centers: &[],
            tetrahedral_carbons: &[],
            double_bond_ends: double_ends,
            stereo_double_bonds: &[],
            etkdg_details: &details,
        };
        embedder_dispatch_workers(&mut args, params, deadline)
    }
    #[test]
    fn explicit_seed_workers_overlap_and_reproduce_serial_coordinates_in_source_order() {
        CALLBACK_ARRIVALS.store(0, Ordering::SeqCst);
        CALLBACK_THREADS.lock().unwrap().clear();
        let bounds = BoundsMatrix::new(1).unwrap();
        let input = [
            Conformer3D::new(0, vec![[9., 8., 7.]], true),
            Conformer3D::new(1, vec![[9., 8., 7.]], true),
        ];
        for sequential in [false, true] {
            CALLBACK_ARRIVALS.store(0, Ordering::SeqCst);
            CALLBACK_THREADS.lock().unwrap().clear();
            let mut serial = input.clone();
            let mut parallel = input.clone();
            let mut sp = params(1);
            let mut pp = params(2);
            sp.enable_sequential_random_seeds = sequential;
            pp.enable_sequential_random_seeds = sequential;
            pp.callback = Some(overlap_callback);
            let mut sok = [true, true];
            let mut pok = [true, true];
            run(&mut serial, &mut sok, &bounds, &mut sp, None, 0, &[], None).unwrap();
            run(
                &mut parallel,
                &mut pok,
                &bounds,
                &mut pp,
                None,
                0,
                &[],
                None,
            )
            .unwrap();
            assert_eq!(sok, pok);
            assert_eq!(pok, [true, true]);
            for (a, b) in serial.iter().zip(&parallel) {
                assert_eq!(a.coordinates(), b.coordinates());
            }
            assert_eq!(pp.max_iterations, 10);
            assert_eq!(pp.basin_thresh, 1e8);
            let tids = CALLBACK_THREADS.lock().unwrap();
            assert_eq!(tids.len(), 2);
            assert_ne!(tids[0], tids[1]);
        }
    }
    #[test]
    fn parallel_failure_counts_are_shared_and_failed_coordinates_remain_unchanged() {
        let mut bounds = BoundsMatrix::new(3).unwrap();
        for (i, j, d) in [(0, 1, 1.), (1, 2, 1.), (0, 2, 2.)] {
            bounds.set_lower(i, j, d).unwrap();
            bounds.set_upper(i, j, d).unwrap();
        }
        let mut confs: Vec<_> = (0..4)
            .map(|id| Conformer3D::new(id, vec![[9., 8., 7.]; 3], true))
            .collect();
        let mut ok = [true; 4];
        let mut p = params(2);
        p.max_iterations = 1;
        p.track_failures = true;
        p.failures = vec![0; crate::EmbedFailureCause::EndOfEnum as usize];
        p.failures[crate::EmbedFailureCause::LinearDoubleBond as usize] = 7;
        p.coord_map = Some(BTreeMap::from([
            (0, [0., 0., 0.]),
            (1, [1., 0., 0.]),
            (2, [2., 0., 0.]),
        ]));
        run(
            &mut confs,
            &mut ok,
            &bounds,
            &mut p,
            None,
            0,
            &[(0, 1, 2)],
            None,
        )
        .unwrap();
        assert_eq!(ok, [false; 4]);
        for conf in &confs {
            assert_eq!(conf.coordinates(), &[[9., 8., 7.]; 3]);
        }
        assert_eq!(
            p.failures[crate::EmbedFailureCause::LinearDoubleBond as usize],
            11
        );
        assert_eq!(p.failures.iter().sum::<u32>(), 11);
    }
    #[test]
    fn parallel_failed_status_and_prior_deadline_do_not_enter_default_updates() {
        let bounds = BoundsMatrix::new(1).unwrap();
        for deadline in [
            None,
            Some(Instant::now() - std::time::Duration::from_secs(1)),
        ] {
            let mut confs = [Conformer3D::new(0, vec![[9., 8., 7.]], true)];
            let mut ok = [deadline.is_some()];
            let initial = ok;
            let mut p = params(3);
            run(&mut confs, &mut ok, &bounds, &mut p, None, 0, &[], deadline).unwrap();
            assert_eq!(ok, initial);
            assert_eq!(confs[0].coordinates(), &[[9., 8., 7.]]);
            assert_eq!(p.max_iterations, 0);
            assert_eq!(p.basin_thresh, 5.);
        }
    }
    #[test]
    fn parallel_seed_and_progress_shape_errors_propagate_with_prior_state() {
        let bounds = BoundsMatrix::new(1).unwrap();
        let mut p = params(2);
        p.random_seed = -2;
        let mut confs = [Conformer3D::new(0, vec![[9., 8., 7.]], true)];
        let mut ok = [true];
        p.failures = vec![3, 4];
        let e = run(&mut confs, &mut ok, &bounds, &mut p, None, 0, &[], None).unwrap_err();
        assert!(matches!(
            e,
            GenerationError::Seed(crate::numeric::EmbedSeedError::InvalidRandomSeed { seed: -2 })
        ));
        assert_eq!(p.failures, [3, 4]);
        assert_eq!(p.max_iterations, 0);
        assert_eq!(p.basin_thresh, 5.);
        assert_eq!(ok, [true]);
        assert_eq!(confs[0].coordinates(), &[[9., 8., 7.]]);
        p.random_seed = 11;
        let e = run(&mut confs, &mut [], &bounds, &mut p, None, 0, &[], None).unwrap_err();
        assert_eq!(e.to_string(), "embedding conformer status count mismatch");
    }
    #[test]
    fn parallel_fragment_writes_are_disjoint_and_preserve_other_fragment_rows() {
        let bounds = BoundsMatrix::new(1).unwrap();
        let mut p = params(2);
        p.coord_map = Some(BTreeMap::from([(0, [1., 2., 3.])]));
        let mut confs: Vec<_> = (0..2)
            .map(|id| Conformer3D::new(id, vec![[9., 8., 7.], [0.; 3], [6., 5., 4.]], true))
            .collect();
        let mut ok = [true, true];
        run(
            &mut confs,
            &mut ok,
            &bounds,
            &mut p,
            Some(&[0, 1, 0]),
            1,
            &[],
            None,
        )
        .unwrap();
        for conf in &confs {
            assert_eq!(
                conf.coordinates(),
                &[[9., 8., 7.], [1., 2., 3.], [6., 5., 4.]]
            );
        }
        assert_eq!(ok, [true, true]);
    }
    #[test]
    fn parallel_callback_exception_joins_workers_and_preserves_entered_defaults() {
        let bounds = BoundsMatrix::new(1).unwrap();
        let mut p = params(2);
        p.callback = Some(throwing_callback);
        p.failures = vec![4, 5];
        let mut confs: Vec<_> = (0..2)
            .map(|id| Conformer3D::new(id, vec![[9., 8., 7.]], true))
            .collect();
        let mut ok = [true, true];
        let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            run(&mut confs, &mut ok, &bounds, &mut p, None, 0, &[], None)
        }));
        let payload = outcome.unwrap_err();
        assert_eq!(
            payload.downcast_ref::<&str>(),
            Some(&"source callback exception")
        );
        assert_eq!(p.max_iterations, 10);
        assert_eq!(p.basin_thresh, 1e8);
        assert_eq!(p.failures, [4, 5]);
        assert_eq!(ok, [true, true]);
        for conf in &confs {
            assert_eq!(conf.coordinates(), &[[9., 8., 7.]]);
        }
    }
}

/// Source warnings retained as detached information for the public projection.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EmbeddingDiagnostic {
    MissingExplicitHydrogens,
    MultipleFragmentsCoordinateMap,
    MultipleFragmentsBoundsMatrix,
    Preparation(PreparationDiagnostic),
    HydrogenRemoval(cosmolkit_core::HydrogenWarning),
    Interrupted,
}
/// Generated coordinate values and source IDs; the facade alone commits them.
#[derive(Debug)]
pub struct GeneratedConformers {
    pub clear_existing: bool,
    pub conformers: Vec<cosmolkit_model::Conformer3D>,
    pub conf_ids: Vec<i32>,
    pub diagnostics: Vec<EmbeddingDiagnostic>,
}

pub fn generate_conformers(
    topology: &TopologyBlock,
    coordinates: &cosmolkit_model::CoordinateBlock,
    properties: &cosmolkit_model::MoleculeProperties,
    num_confs: u32,
    params: &mut EmbedParams,
    project_query: crate::PruningQueryProjection,
) -> Result<GeneratedConformers, GenerationError> {
    // BEGIN RDKIT CPP FUNCTION DGeomHelpers::EmbedMultipleConfs (Embedder.cpp:1501-1694)
    // RDKit❗❌: void EmbedMultipleConfs(ROMol &mol, INT_VECT &res, unsigned int numConfs,
    // RDKit❗❌:                         EmbedParameters &params) {
    // RDKit❗❌:   TimePoint *end_time = nullptr;
    // RDKit❗❌:   TimePoint end_time_storage;
    // RDKit❗❌:   if (params.timeout > 0) {
    // RDKit❗❌:     end_time_storage = Clock::now() + std::chrono::seconds(params.timeout);
    // RDKit❗❌:     end_time = &end_time_storage;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (params.trackFailures) {
    // RDKit❗❌: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit❗❌:     std::lock_guard<std::mutex> lock(GetFailMutex());
    // RDKit❗❌: #endif
    // RDKit❗❌:     params.failures.resize(EmbedFailureCauses::END_OF_ENUM);
    // RDKit❗❌:     std::fill(params.failures.begin(), params.failures.end(), 0);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (!mol.getNumAtoms()) {
    // RDKit❗❌:     throw ValueErrorException("molecule has no atoms");
    // RDKit❗❌:   }
    // RDKit❗❌:   if (params.ETversion < 1 || params.ETversion > 2) {
    // RDKit❗❌:     throw ValueErrorException(
    // RDKit❗❌:         "Only version 1 and 2 of the experimental "
    // RDKit❗❌:         "torsion-angle preferences (ETversion) supported");
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (MolOps::needsHs(mol)) {
    // RDKit❗❌:     BOOST_LOG(rdWarningLog)
    // RDKit❗❌:         << "Molecule does not have explicit Hs. Consider calling AddHs()"
    // RDKit❗❌:         << std::endl;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // initialize the conformers we're going to be creating:
    // RDKit❗❌:   if (params.clearConfs) {
    // RDKit❗❌:     res.clear();
    // RDKit❗❌:     mol.clearConformers();
    // RDKit❗❌:   }
    // RDKit❗❌:   std::vector<std::unique_ptr<Conformer>> confs;
    // RDKit❗❌:   confs.reserve(numConfs);
    // RDKit❗❌:   for (unsigned int i = 0; i < numConfs; ++i) {
    // RDKit❗❌:     confs.emplace_back(new Conformer(mol.getNumAtoms()));
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   boost::dynamic_bitset<> confsOk(numConfs);
    // RDKit❗❌:   confsOk.set();
    // RDKit❗❌:
    // RDKit❗❌:   INT_VECT fragMapping;
    // RDKit❗❌:   std::vector<ROMOL_SPTR> molFrags;
    // RDKit❗❌:   if (params.embedFragmentsSeparately) {
    // RDKit❗❌:     molFrags = MolOps::getMolFrags(mol, true, &fragMapping);
    // RDKit❗❌:   } else {
    // RDKit❗❌:     molFrags.push_back(ROMOL_SPTR(new ROMol(mol)));
    // RDKit❗❌:     fragMapping.resize(mol.getNumAtoms());
    // RDKit❗❌:     std::fill(fragMapping.begin(), fragMapping.end(), 0);
    // RDKit❗❌:   }
    // RDKit❗❌:   const std::map<int, RDGeom::Point3D> *coordMap = params.coordMap;
    // RDKit❗❌:   if (molFrags.size() > 1 && coordMap) {
    // RDKit❗❌:     BOOST_LOG(rdWarningLog)
    // RDKit❗❌:         << "Constrained conformer generation (via the coordMap argument) "
    // RDKit❗❌:            "does not work with molecules that have multiple fragments."
    // RDKit❗❌:         << std::endl;
    // RDKit❗❌:     coordMap = nullptr;
    // RDKit❗❌:   }
    // RDKit❗❌:   boost::dynamic_bitset<> constrainedAtoms(mol.getNumAtoms());
    // RDKit❗❌:   if (coordMap) {
    // RDKit❗❌:     for (const auto &entry : *coordMap) {
    // RDKit❗❌:       constrainedAtoms.set(entry.first);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (molFrags.size() > 1 && params.boundsMat != nullptr) {
    // RDKit❗❌:     BOOST_LOG(rdWarningLog)
    // RDKit❗❌:         << "Conformer generation using a user-provided boundsMat "
    // RDKit❗❌:            "does not work with molecules that have multiple fragments. The "
    // RDKit❗❌:            "boundsMat will be ignored."
    // RDKit❗❌:         << std::endl;
    // RDKit❗❌:     coordMap = nullptr;  // FIXME not directly related to ETKDG, but here I
    // RDKit❗❌:                          // think it should be params.boundsMat = nullptr
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // we will generate conformations for each fragment in the molecule
    // RDKit❗❌:   // separately, so loop over them:
    // RDKit❗❌:   for (unsigned int fragIdx = 0; fragIdx < molFrags.size(); ++fragIdx) {
    // RDKit❗❌:     ROMOL_SPTR piece = molFrags[fragIdx];
    // RDKit❗❌:     unsigned int nAtoms = piece->getNumAtoms();
    // RDKit❗❌:
    // RDKit❗❌:     ForceFields::CrystalFF::CrystalFFDetails etkdgDetails;
    // RDKit❗❌:     etkdgDetails.constrainedAtoms = constrainedAtoms;
    // RDKit❗❌:     EmbeddingOps::initETKDG(piece.get(), params, etkdgDetails);
    // RDKit❗❌:
    // RDKit❗❌:     DistGeom::BoundsMatPtr mmat;
    // RDKit❗❌:     if (params.boundsMat == nullptr || molFrags.size() > 1) {
    // RDKit❗❌:       // The user didn't provide one, so create and initialize the distance
    // RDKit❗❌:       // bounds matrix
    // RDKit❗❌:       mmat.reset(new DistGeom::BoundsMatrix(nAtoms));
    // RDKit❗❌:       initBoundsMat(mmat);
    // RDKit❗❌:       if (!EmbeddingOps::setupInitialBoundsMatrix(piece.get(), mmat, coordMap,
    // RDKit❗❌:                                                   params, etkdgDetails)) {
    // RDKit❗❌:         // return if we couldn't setup the bounds matrix
    // RDKit❗❌:         // possible causes include a triangle smoothing failure
    // RDKit❗❌:         return;
    // RDKit❗❌:       }
    // RDKit❗❌:     } else {
    // RDKit❗❌:       // just use what they gave us
    // RDKit❗❌:       // first make sure it's the right size though:
    // RDKit❗❌:       if (params.boundsMat->numRows() != nAtoms) {
    // RDKit❗❌:         throw ValueErrorException(
    // RDKit❗❌:             "size of boundsMat provided does not match the number of atoms in "
    // RDKit❗❌:             "the molecule.");
    // RDKit❗❌:       }
    // RDKit❗❌:       collectBondsAndAngles((*piece.get()), etkdgDetails.bonds,
    // RDKit❗❌:                             etkdgDetails.angles);
    // RDKit❗❌:       mmat.reset(new DistGeom::BoundsMatrix(*params.boundsMat));
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     // find all the chiral centers in the molecule
    // RDKit❗❌:     MolOps::assignStereochemistry(*piece);
    // RDKit❗❌:     DistGeom::VECT_CHIRALSET chiralCenters;
    // RDKit❗❌:     DistGeom::VECT_CHIRALSET tetrahedralCarbons;
    // RDKit❗❌:     EmbeddingOps::findChiralSets(*piece, chiralCenters, tetrahedralCarbons,
    // RDKit❗❌:                                  coordMap);
    // RDKit❗❌:
    // RDKit❗❌:     // find double bonds
    // RDKit❗❌:     std::vector<std::tuple<unsigned int, unsigned int, unsigned int>>
    // RDKit❗❌:         doubleBondEnds;
    // RDKit❗❌:     std::vector<std::pair<std::vector<unsigned int>, int>> stereoDoubleBonds;
    // RDKit❗❌:     EmbeddingOps::findDoubleBonds(*piece, doubleBondEnds, stereoDoubleBonds,
    // RDKit❗❌:                                   coordMap);
    // RDKit❗❌:
    // RDKit❗❌:     // if we have any chiral centers or are using random coordinates, we
    // RDKit❗❌:     // will first embed the molecule in four dimensions, otherwise we will
    // RDKit❗❌:     // use 3D
    // RDKit❗❌:     bool fourD = false;
    // RDKit❗❌:     if (params.useRandomCoords || chiralCenters.size() > 0) {
    // RDKit❗❌:       fourD = true;
    // RDKit❗❌:     }
    // RDKit❗❌:     int numThreads = getNumThreadsToUse(params.numThreads);
    // RDKit❗❌:
    // RDKit❗❌:     ControlCHandler::reset();
    // RDKit❗❌:
    // RDKit❗❌:     // do the embedding, using multiple threads if requested
    // RDKit❗❌:     detail::EmbedArgs eargs = {&confsOk,        fourD,
    // RDKit❗❌:                                &fragMapping,    &confs,
    // RDKit❗❌:                                fragIdx,         mmat,
    // RDKit❗❌:                                &chiralCenters,  &tetrahedralCarbons,
    // RDKit❗❌:                                &doubleBondEnds, &stereoDoubleBonds,
    // RDKit❗❌:                                &etkdgDetails};
    // RDKit❗❌:     if (numThreads == 1) {
    // RDKit❗❌:       detail::embedHelper_(0, 1, &eargs, &params, end_time);
    // RDKit❗❌:     }
    // RDKit❗❌: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit❗❌:     else {
    // RDKit❗❌:       std::vector<std::future<void>> tg;
    // RDKit❗❌:       for (int tid = 0; tid < numThreads; ++tid) {
    // RDKit❗❌:         tg.emplace_back(std::async(std::launch::async, detail::embedHelper_,
    // RDKit❗❌:                                    tid, numThreads, &eargs, &params, end_time));
    // RDKit❗❌:       }
    // RDKit❗❌:       for (auto &fut : tg) {
    // RDKit❗❌:         fut.get();
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌: #endif
    // RDKit❗❌:     if (end_time != nullptr && Clock::now() > *end_time) {
    // RDKit❗❌:       if (params.trackFailures) {
    // RDKit❗❌: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit❗❌:         std::lock_guard<std::mutex> lock(GetFailMutex());
    // RDKit❗❌: #endif
    // RDKit❗❌:         params.failures[EmbedFailureCauses::EXCEEDED_TIMEOUT]++;
    // RDKit❗❌:       }
    // RDKit❗❌:       res.push_back(-1);
    // RDKit❗❌:       return;
    // RDKit❗❌:     }
    // RDKit❗❌:     if (ControlCHandler::getGotSignal()) {
    // RDKit❗❌:       BOOST_LOG(rdWarningLog) << INTERRUPT_MESSAGE << std::endl;
    // RDKit❗❌:       return;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   std::vector<std::vector<unsigned int>> selfMatches;
    // RDKit❗❌:   if (params.pruneRmsThresh > 0.0) {
    // RDKit❗❌:     selfMatches = detail::getMolSelfMatches(mol, params);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   for (unsigned int ci = 0; ci < confs.size(); ++ci) {
    // RDKit❗❌:     auto &conf = confs[ci];
    // RDKit❗❌:     if (confsOk[ci]) {
    // RDKit❗❌:       // check if we are pruning away conformations and
    // RDKit❗❌:       // a close-by conformation has already been chosen :
    // RDKit❗❌:       if (params.pruneRmsThresh <= 0.0 ||
    // RDKit❗❌:           _isConfFarFromRest(mol, *conf, params.pruneRmsThresh, selfMatches)) {
    // RDKit❗❌:         int confId = (int)mol.addConformer(conf.release(), true);
    // RDKit❗❌:         res.push_back(confId);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION DGeomHelpers::EmbedMultipleConfs (Embedder.cpp:1501-1694)

    // Chemistry remains in the unique domain/core owners. Only source-required
    // fragment graph copies and per-worker points are owned here. Public runtime
    // authority never crosses this detached boundary.
    let end_time = if params.timeout > 0 {
        Some(
            Instant::now()
                .checked_add(std::time::Duration::from_secs(u64::from(params.timeout)))
                .ok_or(GenerationError::Input(
                    "embedding deadline exceeds monotonic clock range",
                ))?,
        )
    } else {
        None
    };
    if params.track_failures {
        let _guard = FAILURE_MUTEX
            .lock()
            .map_err(|_| GenerationError::Input("embedding failure mutex poisoned"))?;
        params
            .failures
            .resize(crate::EmbedFailureCause::EndOfEnum as usize, 0);
        params.failures.fill(0);
    }
    if topology.atoms.is_empty() {
        return Err(GenerationError::Input("molecule has no atoms"));
    }
    if !(1..=2).contains(&params.et_version) {
        return Err(GenerationError::Input(
            "Only version 1 and 2 of the experimental torsion-angle preferences (ETversion) supported",
        ));
    }
    topology.validate()?;
    coordinates.validate_for_atom_count(topology.atoms.len())?;
    let original_valence = cosmolkit_core::assign_valence_for_topology(
        topology,
        cosmolkit_core::ValenceModel::RdkitLike,
    )?;
    let mut result = GeneratedConformers {
        clear_existing: params.clear_conformers,
        conformers: Vec::new(),
        conf_ids: Vec::new(),
        diagnostics: Vec::new(),
    };
    if cosmolkit_forcefields::needs_explicit_hydrogens(
        topology,
        &original_valence.implicit_hydrogens,
    )? {
        result
            .diagnostics
            .push(EmbeddingDiagnostic::MissingExplicitHydrogens);
    }
    let mut candidates: Vec<_> = (0..num_confs)
        .map(|id| {
            cosmolkit_model::Conformer3D::new(
                id as usize,
                vec![[0.; 3]; topology.atoms.len()],
                true,
            )
        })
        .collect();
    let mut ok = vec![true; num_confs as usize];
    // CK's approved coordinate contract changes only the declared 3D dimension;
    // retained 2D values are borrowed by the fragment owner without clearing.
    let fragment_coordinates =
        cosmolkit_core::FragmentCoordinateView::from_coordinate_block(coordinates);
    let fragment_coordinates = if params.clear_conformers {
        fragment_coordinates.without_3d_conformers()
    } else {
        fragment_coordinates
    };
    let fragments = if params.embed_fragments_separately {
        cosmolkit_core::get_molecule_fragments_with_coordinate_view(
            topology,
            &fragment_coordinates,
            properties,
            true,
            true,
        )?
    } else {
        Vec::new()
    };
    let count = if params.embed_fragments_separately {
        fragments.len()
    } else {
        1
    };
    let mut frag_mapping = vec![0; topology.atoms.len()];
    for (frag, value) in fragments.iter().enumerate() {
        for &atom in value.component_atoms() {
            frag_mapping[atom.index()] = frag;
        }
    }
    let mut coord_map = params.coord_map.as_ref();
    if count > 1 && coord_map.is_some() {
        result
            .diagnostics
            .push(EmbeddingDiagnostic::MultipleFragmentsCoordinateMap);
        coord_map = None;
    }
    let mut constrained_atoms = vec![false; topology.atoms.len()];
    if let Some(map) = coord_map {
        for &key in map.keys() {
            *constrained_atoms
                .get_mut(key as u32 as usize)
                .ok_or(GenerationError::Input(
                    "constrained atom index out of range",
                ))? = true;
        }
    }
    if count > 1 && params.bounds_mat.is_some() {
        result
            .diagnostics
            .push(EmbeddingDiagnostic::MultipleFragmentsBoundsMatrix);
        coord_map = None;
    }
    // Move only the source coordMap handle temporarily out of options so graph
    // preparation can borrow it while worker progress borrows params mutably.
    // Keep the source's separate local coordMap vs params.coordMap distinction.
    let local_map_enabled = coord_map.is_some();
    let mut pieces: Vec<_> = if params.embed_fragments_separately {
        fragments
            .into_iter()
            .map(|f| {
                let (t, c, p) = f.into_parts();
                (t, c, p)
            })
            .collect()
    } else {
        // Source new ROMol(mol) copies the graph. Coordinates/properties needed
        // here are borrowed in the nonfragment path; no full coordinate clone.
        Vec::new()
    };
    for frag_idx in 0..count {
        let (mut piece, piece_properties) = if params.embed_fragments_separately {
            let (t, _, p) = std::mem::take(&mut pieces[frag_idx]);
            (t, p)
        } else {
            (topology.clone(), properties.clone())
        };
        let valence = cosmolkit_core::assign_valence_for_topology(
            &piece,
            cosmolkit_core::ValenceModel::RdkitLike,
        )?;
        let rings = cosmolkit_core::symmetrize_sssr_with_options_from_parts(
            piece.atoms.len(),
            &piece.bonds,
            &piece.adjacency,
            false,
            false,
        )?;
        let conjugated = cosmolkit_core::assign_conjugation_flags(&piece, &valence)?;
        let hybridizations =
            cosmolkit_core::assign_hybridization_with_conjugation(&piece, &valence, &conjugated)?
                .values;
        let prepared = PreparedEmbeddingTopology {
            topology: &piece,
            rings: &rings,
            valence: &valence,
            hybridizations: &hybridizations,
            conjugated: &conjugated,
        };
        let mut details = CrystalFFDetails::default();
        details.constrained_atoms = constrained_atoms.clone();
        embedder_init_etkdg(&prepared, params, &mut details)?;
        let mut preparation = Vec::new();
        let local_map = if local_map_enabled {
            params.coord_map.as_ref()
        } else {
            None
        };
        let bounds = if params.bounds_mat.is_none() || count > 1 {
            let mut bounds = BoundsMatrix::new(piece.atoms.len())?;
            // RDKit❗✔️:       initBoundsMat(mmat);
            // Same existing dense initialization owner and O(N^2) writes as source.
            crate::graph_bounds::init_bounds_mat(&mut bounds, 0.0, 1000.0)
                .map_err(GenerationFailure::from)?;
            if !embedder_setup_initial_bounds_matrix(
                &prepared,
                &mut bounds,
                local_map,
                params,
                &mut details,
                &mut preparation,
            )? {
                result.diagnostics.extend(
                    preparation
                        .into_iter()
                        .map(EmbeddingDiagnostic::Preparation),
                );
                return Ok(result);
            }
            bounds
        } else {
            let bounds = params
                .bounds_mat
                .as_ref()
                .expect("source checked bounds presence");
            if bounds.dimension() != piece.atoms.len() {
                return Err(GenerationError::Input(
                    "size of boundsMat provided does not match the number of atoms in the molecule.",
                ));
            }
            crate::graph_bounds::collect_bonds_and_angles(
                &piece,
                &mut details.bonds,
                &mut details.angles,
            );
            (**bounds).clone()
        };
        if !piece_properties
            .props()
            .contains_key(b"_StereochemDone".as_slice())
        {
            piece = cosmolkit_core::assign_legacy_stereochemistry_with_flags(
                piece, &valence, &rings, false, false,
            )?;
        }
        let mut chiral = Vec::new();
        let mut tetrahedral = Vec::new();
        embedder_find_chiral_sets(
            &piece,
            Some(&rings),
            &mut chiral,
            &mut tetrahedral,
            local_map,
            &mut preparation,
        )?;
        let mut double_ends = Vec::new();
        let mut stereo_double = Vec::new();
        embedder_find_double_bonds(&piece, &mut double_ends, &mut stereo_double, local_map)?;
        result.diagnostics.extend(
            preparation
                .into_iter()
                .map(EmbeddingDiagnostic::Preparation),
        );
        let mut args = EmbedHelperArgs {
            confs_ok: &mut ok,
            four_d: params.use_random_coords || !chiral.is_empty(),
            frag_mapping: Some(&frag_mapping),
            confs: &mut candidates,
            frag_idx,
            mmat: &bounds,
            chiral_centers: &chiral,
            tetrahedral_carbons: &tetrahedral,
            double_bond_ends: &double_ends,
            stereo_double_bonds: &stereo_double,
            etkdg_details: &details,
        };
        embedder_dispatch_workers(&mut args, params, end_time)?;
        if end_time.is_some_and(|time| Instant::now() > time) {
            with_embedding_progress(params, |settings, progress| {
                embedder_increment_failure(
                    settings,
                    progress,
                    crate::EmbedFailureCause::ExceededTimeout,
                )
            })?;
            result.conf_ids.push(-1);
            return Ok(result);
        }
        if source_got_signal() {
            result.diagnostics.push(EmbeddingDiagnostic::Interrupted);
            return Ok(result);
        }
    }
    let matches = if params.prune_rms_thresh > 0. {
        let pruned =
            crate::pruning_self_matches(topology, coordinates, properties, params, project_query)?;
        result.diagnostics.extend(
            pruned
                .hydrogen_warnings
                .into_iter()
                .map(EmbeddingDiagnostic::HydrogenRemoval),
        );
        pruned.self_matches
    } else {
        Vec::new()
    };
    let retained = if params.clear_conformers {
        &[][..]
    } else {
        coordinates.conformers_3d.as_slice()
    };
    // BEGIN RDKIT CPP FUNCTION ROMol::addConformer (ROMol.cpp:665-682)
    // RDKit❗🔝: unsigned int ROMol::addConformer(Conformer *conf, bool assignId) {
    // RDKit❗🔝:   PRECONDITION(conf, "bad conformer");
    // RDKit❗🔝:   PRECONDITION(conf->getNumAtoms() == this->getNumAtoms(),
    // RDKit❗🔝:                "Number of atom mismatch");
    // RDKit❗🔝:   if (assignId) {
    // RDKit❗🔝:     int maxId = -1;
    // RDKit❗🔝:     for (auto cptr : d_confs) {
    // RDKit❗🔝:       maxId = std::max((int)(cptr->getId()), maxId);
    // RDKit❗🔝:     }
    // RDKit❗🔝:     maxId++;
    // RDKit❗🔝:     conf->setId((unsigned int)maxId);
    // RDKit❗🔝:   }
    // RDKit❗🔝:   conf->setOwningMol(this);
    // RDKit❗🔝:   CONFORMER_SPTR nConf(conf);
    // RDKit❗🔝:   d_confs.push_back(nConf);
    // RDKit❗🔝:   return conf->getId();
    // RDKit❗🔝: }
    // RDKit❗🔝:
    // END RDKIT CPP FUNCTION ROMol::addConformer
    // The next source ID is observed only for accepted conformers. For the
    // source valid signed ID range, a running maximum preserves source IDs and
    // reduces repeated O(C+K) scans to one O(C) scan plus O(K) increments.
    let mut next_id = None;
    for (ci, conf) in candidates.into_iter().enumerate() {
        if ok[ci]
            && (params.prune_rms_thresh <= 0.
                || embedder_is_conf_far_from_rest(
                    retained.iter().chain(result.conformers.iter()),
                    &conf,
                    params.prune_rms_thresh,
                    &matches,
                )?)
        {
            let id = match next_id {
                Some(id) => id,
                None => retained.iter().map(|c| c.id()).max().map_or(Ok(0), |id| {
                    id.checked_add(1)
                        .ok_or(GenerationError::Input("conformer id exceeds model range"))
                })?,
            };
            let source_id = i32::try_from(id)
                .map_err(|_| GenerationError::Input("conformer id exceeds source signed range"))?;
            result.conf_ids.push(source_id);
            result.conformers.push(conf.with_id(id));
            next_id = Some(
                id.checked_add(1)
                    .ok_or(GenerationError::Input("conformer id exceeds model range"))?,
            );
        }
    }
    Ok(result)
}

#[cfg(test)]
mod original_complete_generation_conditions {
    use super::*;
    use cosmolkit_model::{
        Atom, AtomSpec, Conformer2D, Conformer3D, CoordinateBlock, Element, MoleculeProperties,
        QueryGraph, QueryGraphError,
    };
    use std::sync::Arc;
    fn unused_projection(
        _: &TopologyBlock,
        _: &CoordinateBlock,
        _: &MoleculeProperties,
    ) -> Result<QueryGraph, QueryGraphError> {
        panic!("source non-symmetry branch must not call query projection")
    }
    fn carbon() -> TopologyBlock {
        let mut topology = TopologyBlock::default();
        topology
            .atoms
            .push(Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)));
        topology.adjacency = cosmolkit_model::AdjacencyList::try_from_topology(1, &[]).unwrap();
        topology
    }
    fn smiles(text: &str) -> TopologyBlock {
        let raw = cosmolkit_smiles::parse_smiles(text, &Default::default()).unwrap();
        cosmolkit_core::sanitize_topology(&raw.topology, &Default::default())
            .unwrap()
            .topology
    }
    fn random_embed_params(seed: i32, num_threads: i32) -> EmbedParams {
        EmbedParams {
            max_iterations: 1,
            random_seed: seed,
            use_random_coords: true,
            num_threads,
            box_size_mult: -2.0,
            use_symmetry_for_pruning: false,
            symmetrize_conjugated_terminal_groups_for_pruning: false,
            ..Default::default()
        }
    }
    // Original distgeom/tests.rs:4999-5010; value transport only. The original
    // public Molecule-returning condition remains required at the facade layer.
    #[test]
    fn embed_multiple_confs_generates_value_style_conformers() {
        let mol = carbon();
        let coordinates = CoordinateBlock::default();
        let mut params = random_embed_params(7, 1);
        let embedded = generate_conformers(
            &mol,
            &coordinates,
            &Default::default(),
            2,
            &mut params,
            unused_projection,
        )
        .expect("embed");
        assert!(coordinates.conformers_3d.is_empty());
        assert_eq!(embedded.conf_ids, vec![0, 1]);
        assert_eq!(embedded.conformers.len(), 2);
    }
    // Original distgeom/tests.rs:5033-5049; fixed source options, diagnostic and
    // expectation retained, replacing only the authorized detached boundary.
    #[test]
    fn conformer_generation_parameter_parity_custom_bounds_matrix_size_check_matches_rdkit() {
        let mol = smiles("CC");
        let mut wrong_size_params = EmbedParams::etkdg();
        wrong_size_params.random_seed = 42;
        wrong_size_params.num_threads = 1;
        wrong_size_params.bounds_mat =
            Some(Arc::new(BoundsMatrix::new(mol.atoms.len() + 1).unwrap()));
        let err = generate_conformers(
            &mol,
            &Default::default(),
            &Default::default(),
            1,
            &mut wrong_size_params,
            unused_projection,
        )
        .expect_err("RDKit custom bounds matrix size mismatch must error");
        assert!(
            err.to_string().contains(
                "size of boundsMat provided does not match the number of atoms in the molecule"
            ),
            "unexpected custom bounds error: {err}"
        );
    }
    #[test]
    fn whole_fragment_mapping_is_serial_parallel_identical_and_does_not_mutate_input() {
        let mol = smiles("C.C");
        let coordinates = CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(7, vec![[1., 2.], [3., 4.]])],
            ..Default::default()
        };
        let before = coordinates.clone();
        let mut serial = random_embed_params(7, 1);
        let mut parallel = random_embed_params(7, 2);
        let a = generate_conformers(
            &mol,
            &coordinates,
            &Default::default(),
            3,
            &mut serial,
            unused_projection,
        )
        .unwrap();
        let b = generate_conformers(
            &mol,
            &coordinates,
            &Default::default(),
            3,
            &mut parallel,
            unused_projection,
        )
        .unwrap();
        assert_eq!(a.conf_ids, vec![0, 1, 2]);
        assert_eq!(a.conf_ids, b.conf_ids);
        assert_eq!(a.conformers, b.conformers);
        for conf in &a.conformers {
            assert_eq!(conf.coordinates()[0], conf.coordinates()[1]);
        }
        assert_eq!(coordinates, before);
        assert!(a.clear_existing);
    }
    #[test]
    fn whole_pruning_retained_ids_and_coordinate_dimensions_follow_source_order() {
        let mol = smiles("CC");
        let coordinates = CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(7, vec![[3., 4.], [5., 6.]])],
            conformers_3d: vec![Conformer3D::new(7, vec![[0., 0., 0.], [1., 0., 0.]], true)],
            ..Default::default()
        };
        let before = coordinates.clone();
        let mut params = random_embed_params(7, 1);
        params.clear_conformers = false;
        params.coord_map = Some(BTreeMap::from([(0, [0., 0., 0.]), (1, [1., 0., 0.])]));
        params.prune_rms_thresh = 0.5;
        let pruned = generate_conformers(
            &mol,
            &coordinates,
            &Default::default(),
            2,
            &mut params,
            unused_projection,
        )
        .unwrap();
        assert!(pruned.conf_ids.is_empty());
        assert!(pruned.conformers.is_empty());
        assert!(!pruned.clear_existing);
        params.prune_rms_thresh = -1.;
        let appended = generate_conformers(
            &mol,
            &coordinates,
            &Default::default(),
            2,
            &mut params,
            unused_projection,
        )
        .unwrap();
        assert_eq!(appended.conf_ids, vec![8, 9]);
        assert_eq!(
            appended.conformers[0].coordinates(),
            coordinates.conformers_3d[0].coordinates()
        );
        params.clear_conformers = true;
        let cleared = generate_conformers(
            &mol,
            &coordinates,
            &Default::default(),
            1,
            &mut params,
            unused_projection,
        )
        .unwrap();
        assert_eq!(cleared.conf_ids, vec![0]);
        assert!(cleared.clear_existing);
        assert_eq!(coordinates, before);
    }
    #[test]
    fn whole_error_order_resets_failure_counts_before_empty_and_et_version_checks() {
        let mut params = random_embed_params(7, 1);
        params.track_failures = true;
        params.failures = vec![9, 8, 7];
        params.et_version = 0;
        let empty = generate_conformers(
            &Default::default(),
            &Default::default(),
            &Default::default(),
            1,
            &mut params,
            unused_projection,
        )
        .unwrap_err();
        assert_eq!(empty.to_string(), "molecule has no atoms");
        assert_eq!(
            params.failures,
            vec![0; crate::EmbedFailureCause::EndOfEnum as usize]
        );
        let invalid = generate_conformers(
            &carbon(),
            &Default::default(),
            &Default::default(),
            1,
            &mut params,
            unused_projection,
        )
        .unwrap_err();
        assert_eq!(
            invalid.to_string(),
            "Only version 1 and 2 of the experimental torsion-angle preferences (ETversion) supported"
        );
    }
    fn slow_callback(_: u32) {
        std::thread::sleep(std::time::Duration::from_millis(1100));
    }
    #[test]
    fn whole_timeout_returns_source_minus_one_and_records_failure() {
        let mut params = random_embed_params(7, 1);
        params.timeout = 1;
        params.callback = Some(slow_callback);
        params.track_failures = true;
        let result = generate_conformers(
            &carbon(),
            &Default::default(),
            &Default::default(),
            1,
            &mut params,
            unused_projection,
        )
        .unwrap();
        assert_eq!(result.conf_ids, vec![-1]);
        assert!(result.conformers.is_empty());
        assert_eq!(
            params.failures[crate::EmbedFailureCause::ExceededTimeout as usize],
            1
        );
    }
    #[test]
    fn whole_etkdgv3_addhs_shared_owners_reproduce_serial_parallel_coordinates() {
        let hydrogens = cosmolkit_core::add_hydrogens_impl(
            smiles("CCO"),
            Default::default(),
            Default::default(),
        )
        .unwrap();
        let mut serial = EmbedParams::etkdg_v3();
        serial.random_seed = 42;
        serial.max_iterations = 3;
        serial.num_threads = 1;
        let mut parallel = serial.clone();
        parallel.num_threads = 2;
        let a = generate_conformers(
            &hydrogens.topology,
            &hydrogens.coordinates,
            &hydrogens.properties,
            3,
            &mut serial,
            unused_projection,
        )
        .unwrap();
        let b = generate_conformers(
            &hydrogens.topology,
            &hydrogens.coordinates,
            &hydrogens.properties,
            3,
            &mut parallel,
            unused_projection,
        )
        .unwrap();
        assert_eq!(a.conf_ids, vec![0, 1, 2]);
        assert_eq!(a.conf_ids, b.conf_ids);
        assert_eq!(a.conformers, b.conformers);
        assert!(hydrogens.coordinates.conformers_3d.is_empty());
    }
    #[test]
    fn no_candidates_do_not_observe_retained_id_or_mutate_worker_defaults() {
        let coordinates = CoordinateBlock {
            conformers_3d: vec![Conformer3D::new(usize::MAX, vec![[0.; 3]], true)],
            ..Default::default()
        };
        let mut params = random_embed_params(7, 1);
        params.max_iterations = 0;
        params.clear_conformers = false;
        let result = generate_conformers(
            &carbon(),
            &coordinates,
            &Default::default(),
            0,
            &mut params,
            unused_projection,
        )
        .unwrap();
        assert!(result.conf_ids.is_empty());
        assert!(result.conformers.is_empty());
        assert_eq!(params.max_iterations, 0);
    }
    #[test]
    #[ignore = "debug helper for row-1 torsion contrib investigation"]
    fn debug_ethene_row1_torsion_collection_vs_single_on_basiconly_coords() {
        let mol = cosmolkit_core::add_hydrogens_impl(
            smiles("C=C"),
            Default::default(),
            Default::default(),
        )
        .expect("add hs");
        let mut params = EmbedParams::etkdg_v3();
        params.use_exp_torsion_angle_prefs = false;
        params.random_seed = 61453;
        params.num_threads = 1;
        params.timeout = 10;
        let embedded = generate_conformers(
            &mol.topology,
            &mol.coordinates,
            &mol.properties,
            1,
            &mut params,
            unused_projection,
        )
        .expect("embed");
        let coords = embedded.conformers[0].coordinates();
        let rings = cosmolkit_core::symmetrized_sssr(&mol.topology, &Default::default())
            .expect("original ring state");
        let valence = cosmolkit_core::assign_valence_for_topology(
            &mol.topology,
            cosmolkit_core::ValenceModel::RdkitLike,
        )
        .expect("original valence state");
        let mut details = CrystalFFDetails::default();
        cosmolkit_forcefields::get_experimental_torsions_without_bonds(
            &mol.topology,
            &rings,
            &valence,
            &mut details,
            true,
            params.use_small_ring_torsions,
            params.use_macrocycle_torsions,
            true,
            params.et_version,
            params.verbose,
        )
        .expect("details");
        let atoms = &details.exp_torsion_atoms[0];
        let (signs, force_constants) = &details.exp_torsion_angles[0];
        let evaluation = cosmolkit_forcefields::evaluate_crystal_torsion_pair(
            coords,
            [
                atoms[0] as usize,
                atoms[1] as usize,
                atoms[2] as usize,
                atoms[3] as usize,
            ],
            force_constants,
            signs,
        )
        .expect("original unique collection/single evaluations");
        let collection_energy = evaluation.collection_energy;
        let single_energy = evaluation.single_energy;
        let collection_grad = evaluation.collection_gradient;
        let single_grad = evaluation.single_gradient;
        println!("collection_energy={collection_energy:.16}");
        println!("single_energy={single_energy:.16}");
        println!("collection_grad={collection_grad:?}");
        println!("single_grad={single_grad:?}");
        assert_eq!(collection_grad, single_grad);
        assert_eq!(collection_energy, single_energy);
    }
}

#[cfg(test)]
mod etv2_first_numeric_boundary_diagnostic {
    use super::*;
    #[test]
    fn etv2_probe3_same_native_numeric_checkpoints() {
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../testdata/conformer/fixtures/rdkit/test_data/torsion.etkdg.v2.mol");
        let (topology, _, _) =
            cosmolkit_io::sdf::read_v2000_detached(&std::fs::read_to_string(path).unwrap())
                .unwrap();
        let topology = cosmolkit_core::sanitize_topology(&topology, &Default::default())
            .unwrap()
            .topology;
        let valence = cosmolkit_core::assign_valence_for_topology(
            &topology,
            cosmolkit_core::ValenceModel::RdkitLike,
        )
        .unwrap();
        let rings = cosmolkit_core::symmetrized_sssr(&topology, &Default::default()).unwrap();
        let conjugated = cosmolkit_core::assign_conjugation_flags(&topology, &valence).unwrap();
        let hybridizations =
            cosmolkit_core::assign_hybridization_with_conjugation(&topology, &valence, &conjugated)
                .unwrap();
        let bounds = crate::graph_bounds::build_bounds_matrix(
            &topology,
            &rings,
            &valence,
            &hybridizations.values,
            &conjugated,
            true,
            false,
            true,
            false,
        )
        .unwrap();
        let params = EmbedParams::etdg_v2();
        let mut positions = vec![vec![0.; 3]; 17];
        let mut dist = SymmMatrix::new(17);
        let mut rng = cosmolkit_core::RdkitRandomEngine::from_seed(42);
        pick_random_dist_mat_with_rng(&bounds, &mut dist, &mut rng);
        assert!(compute_initial_coords_with_rng(&dist, &mut positions, &mut rng, true, 1).unwrap());
        let initial = positions.clone();
        println!("etv2_checkpoint_initial {:?}", initial);
        let fixed = vec![false; 17];
        let (first_energy_before, first_energy_after);
        {
            let mut field = ConformerOptimizer::distance_geometry(
                &bounds,
                &mut positions,
                &[],
                DistanceGeometryForceFieldParams {
                    weight_chiral: 1.,
                    weight_fourth_dimension: 0.1,
                    extra_weights: None,
                    basin_size_tolerance: 5.,
                    fixed_pair_points: Some(&fixed),
                    fixed_points: &[],
                },
            )
            .unwrap();
            first_energy_before = field.energy(None).unwrap();
            if first_energy_before > 1e-5 {
                while field.minimize(400, 1e-3, 1e-6).unwrap() != 0 {}
            }
            first_energy_after = field.energy(None).unwrap();
        }
        let first_minimized = positions.clone();
        println!("etv2_checkpoint_first_minimized {:?}", first_minimized);
        let mut details = CrystalFFDetails::default();
        cosmolkit_forcefields::get_experimental_torsions_without_bonds(
            &topology,
            &rings,
            &valence,
            &mut details,
            true,
            false,
            false,
            false,
            2,
            false,
        )
        .unwrap();
        details.atom_nums = topology
            .atoms
            .iter()
            .map(|a| i32::from(a.atomic_number()))
            .collect();
        crate::graph_bounds::set_topol_bounds_with_outputs(
            &topology,
            &rings,
            &valence,
            &hybridizations.values,
            &conjugated,
            &mut bounds.clone(),
            &mut details.bonds,
            &mut details.angles,
            true,
            false,
            false,
            true,
            true,
            true,
        )
        .unwrap();
        let (torsion_energy_before, torsion_energy_after);
        {
            let mut field =
                ConformerOptimizer::torsions(&bounds, &mut positions, &details, false, None, &[])
                    .unwrap();
            torsion_energy_before = field.energy(None).unwrap();
            if torsion_energy_before > 1e-5 {
                let _ = field
                    .minimize(300, params.optimizer_force_tol, 1e-6)
                    .unwrap();
            }
            torsion_energy_after = field.energy(None).unwrap();
        }
        println!("etv2_checkpoint_final {:?}", positions);
        println!(
            "etv2_checkpoint_energies {:?}",
            (
                first_energy_before,
                first_energy_after,
                torsion_energy_before,
                torsion_energy_after
            )
        );
        let expected: serde_json::Value=serde_json::from_str(r#"{"smoothing":1,"initial_ok":1,"initial":[[0.8398722485412821,-1.2145847914991603,-0.7097566336972684],[1.6407970818902262,-1.0837697982274672,0.03743649066764579],[2.2676597681674866,-0.24070088035107554,-0.13186567594308204],[1.1359658098504126,0.922545014941577,0.8637006022619268],[-0.15252430168387943,0.7631031530416291,0.10787659784358963],[-0.920983182793755,1.5688146158529852,-0.20960920796911744],[-0.4112441961069539,-0.2716961035274239,-0.6237743072596654],[-1.8907659246435478,-0.55927478030847,-1.0831336687694642],[-1.2024147414160316,-2.2670201496534497,0.19063170206276314],[2.5096320653385416,-0.4257961142200684,-1.854873104074399],[2.764570507474988,-0.3417216773199005,1.480437948583766],[0.4776533337342171,0.4794415118069911,1.9245162175948498],[-1.288235042109195,2.399975872738075,0.4506333359778203],[-1.3902138846927112,1.644370968929529,0.06226466209981267],[-0.1962426969498886,2.2913057940338417,-1.355606983541704],[-2.287171651572724,-1.6517150132987442,1.8964626556942887],[-1.896355193028466,-2.013277622938868,-1.0453406315317606]],"first_energy_before":34.21043237438536,"first_grad_before":[1.9537251246547036,-0.6136501080829042,0.5264739363936116,-8.930705954383027,-6.194492110358579,18.96952792883522,0.6306463561454189,-2.6568599680189235,-10.791153742820503,-0.4565837268598285,2.8896658329839755,0.442380036161131,-1.2166233769847035,0.502497930563496,-0.3099359404422125,-2.8848089650983098,-2.483623249349017,3.9276524753226623,-0.12158463360225062,1.2417146841609537,0.4190135960777228,-3.4295315786024094,6.997923551770949,-9.007970917583648,16.470971967684264,-13.08689364355029,-10.42786541815178,8.621304687819078,6.640922256166689,-20.507983766893762,2.3885068080050607,-0.4681657635580031,9.09033116348165,-2.089832779514774,-0.47033008889283445,2.501395013745577,-0.4806612091640632,-0.7576113802629263,0.24108224294792824,1.489068824097914,0.18613203824090102,-0.11095277432014503,2.6164068164967746,2.0693994290830537,-4.067378876127686,-13.416136121044156,5.945910315335591,30.41514862863086,-1.144162239649693,0.2574602737678694,-11.30976358525662],"first_minimized":[[0.7769158159289418,-1.1936041102920902,-1.0887303570001792],[2.0174083098874256,-0.9125572528999197,-0.6183135741224527],[2.0833637403784944,-0.13199285285431517,0.5290155445979855],[1.0743857696038492,0.7487139127195357,0.8663244789324498],[-0.11887322082374396,0.6290312522237073,0.16381184258489245],[-0.9761798801918173,1.8146597741553063,-0.0650828231712111],[-0.2774015425287768,-0.40952317496950424,-0.7484876416686083],[-1.5861952306930736,-0.8997967159249963,-0.8520724224189472],[-1.8586694036885554,-1.9449779632080146,0.018223289833461852],[2.918055463181744,-1.0595145486303588,-1.1848920376726932],[2.9414531737944243,-0.19966533339217885,1.1652260256859006],[1.1344670346634158,1.461343884827805,1.6664130563314952],[-0.7603500065577631,2.5209463963043017,0.7688832890891262],[-2.050727322929517,1.6069498844250805,-0.025947019534324746],[-0.6544180641842868,2.338001861153427,-1.000142352554056],[-2.351954418788671,-1.6332999016006724,0.874717101478015],[-2.31128021705207,-2.7347151120370916,-0.4689464003908526]],"first_energy_after":0.0006775778282256319,"torsion_energy_before":11.959224589345164,"torsion_grad_before":[-0.6298185177312893,3.388243895375178,-8.630050504596063,-0.19134167883691044,0.17431598849791743,0.03536094277162857,-0.31349861348129016,0.043681314304477366,-0.09444276663178007,-0.13746239539280924,-0.18561123031999105,-0.1991419200320264,0,0,0,0.4354518615745667,-0.35589273455110093,-0.09667406086036254,-0.2983874511665009,-2.440020060297632,15.31889990193863,4.042292196698009,-5.9633688029062215,-12.05488101522353,-2.5665525937640057,5.600244365810473,5.922068813338524,-0.06275923793194195,0.06155280451615922,0.06982421301601353,0,0,0,-0.2861851448014905,-0.3027613188150826,-0.2704522789782399,0,0,0,0,0,0,0.008261574833662695,-0.020384221614176492,-0.000511324742790685,0,0,0,0,0,0],"torsion_need_more":0,"final":[[0.6114985012458053,-1.428424782573838,-0.42189385882078884],[1.9252955924757869,-1.0628879893825556,-0.3468123258144795],[2.2105935270737804,0.0621254896097724,0.4198088068363781],[1.2320821775292117,1.0014016328618198,0.7078547390773309],[-0.027145692996548232,0.7787260650721753,0.1882948780665941],[-0.9746767234884467,1.87874821312757,-0.05313698558826245],[-0.3328887198971735,-0.44309520372702144,-0.42607790497117737],[-1.700541020808343,-0.7847075396176533,-0.38098583518836326],[-1.9445004273514868,-2.1410024205420473,-0.3354832557713035],[2.7207946939522785,-1.572957159243348,-0.8417806890841792],[3.1526785172960823,0.12501871483421542,0.94479007626298],[1.4420270130397295,1.9151036600905664,1.2438028367855292],[-0.8367957497105144,2.6192951293982243,0.7804111139311738],[-2.0335690123292745,1.584616832465797,0.006846842668979736],[-0.7124093167120968,2.4356122177489286,-0.9899375528884584],[-2.199383888373151,-2.4916947539758847,0.6092162950350409],[-2.5330594709455956,-2.475878106146647,-1.1049171805370015]],"torsion_energy_after":2.5787096695387443e-07}"#).unwrap();
        let mut failures = Vec::new();
        for (key, rows) in [
            ("initial", &initial),
            ("first_minimized", &first_minimized),
            ("final", &positions),
        ] {
            let expected_rows: Vec<Vec<f64>> =
                serde_json::from_value(expected[key].clone()).unwrap();
            let maximum = rows
                .iter()
                .zip(&expected_rows)
                .flat_map(|(a, b)| a.iter().zip(b))
                .map(|(a, b)| (a - b).abs())
                .fold(0., f64::max);
            println!("etv2_checkpoint_delta {key} {maximum}");
            if maximum > 1e-9 {
                failures.push(key);
            }
        }
        for (key, value) in [
            ("first_energy_before", first_energy_before),
            ("first_energy_after", first_energy_after),
            ("torsion_energy_before", torsion_energy_before),
            ("torsion_energy_after", torsion_energy_after),
        ] {
            let delta = (value - expected[key].as_f64().unwrap()).abs();
            println!("etv2_checkpoint_delta {key} {delta}");
            if delta > 1e-9 {
                failures.push(key);
            }
        }
        assert!(
            failures.is_empty(),
            "first numeric checkpoint divergence {:?}",
            failures
        );
    }
}

#[cfg(test)]
mod complete_generated_bounds_regression {
    use super::*;
    #[test]
    fn whole_etdgv2_fixture_reproduces_original_coordinates() {
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../testdata/conformer/fixtures/rdkit/test_data/torsion.etkdg.v2.mol");
        let (topology, coordinates, properties) =
            cosmolkit_io::sdf::read_v2000_detached(&std::fs::read_to_string(path).unwrap())
                .unwrap();
        let topology = cosmolkit_core::sanitize_topology(&topology, &Default::default())
            .unwrap()
            .topology;
        let mut params = EmbedParams::etdg_v2();
        params.random_seed = 42;
        let embedded = generate_conformers(
            &topology,
            &coordinates,
            &properties,
            1,
            &mut params,
            |_, _, _| unreachable!("pruning disabled in original ETv2 case"),
        )
        .unwrap();
        assert_eq!(embedded.conf_ids, vec![0]);
        let expected: Vec<[f64; 3]> = serde_json::from_str(r#"[[0.6114985012458053, -1.428424782573838, -0.42189385882078884], [1.9252955924757869, -1.0628879893825556, -0.3468123258144795], [2.2105935270737804, 0.0621254896097724, 0.4198088068363781], [1.2320821775292117, 1.0014016328618198, 0.7078547390773309], [-0.027145692996548232, 0.7787260650721753, 0.1882948780665941], [-0.9746767234884467, 1.87874821312757, -0.05313698558826245], [-0.3328887198971735, -0.44309520372702144, -0.42607790497117737], [-1.700541020808343, -0.7847075396176533, -0.38098583518836326], [-1.9445004273514868, -2.1410024205420473, -0.3354832557713035], [2.7207946939522785, -1.572957159243348, -0.8417806890841792], [3.1526785172960823, 0.12501871483421542, 0.94479007626298], [1.4420270130397295, 1.9151036600905664, 1.2438028367855292], [-0.8367957497105144, 2.6192951293982243, 0.7804111139311738], [-2.0335690123292745, 1.584616832465797, 0.006846842668979736], [-0.7124093167120968, 2.4356122177489286, -0.9899375528884584], [-2.199383888373151, -2.4916947539758847, 0.6092162950350409], [-2.5330594709455956, -2.475878106146647, -1.1049171805370015]]"#).unwrap();
        let actual = embedded.conformers[0].coordinates();
        assert_eq!(actual.len(), expected.len());
        for (atom, (actual, expected)) in actual.iter().zip(&expected).enumerate() {
            for axis in 0..3 {
                assert!(
                    (actual[axis] - expected[axis]).abs() <= 1e-6,
                    "unchanged original ETDGv2 atom {atom} axis {axis}: actual={} expected={}",
                    actual[axis],
                    expected[axis]
                );
            }
        }
    }
}
