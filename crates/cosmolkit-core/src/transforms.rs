//! Source-backed detached 3D coordinate transforms.
//!
//! This module owns numerical/domain behavior only. It has no live molecule,
//! operation capability, registry, or cache authority.

use std::f64::consts::PI;

use cosmolkit_model::{
    AtomId, BondId, Conformer3D, CoordinateBlock, CoordinateValidationError, TopologyBlock,
    TopologyValidationError,
};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Transform3D {
    values: [f64; 16],
}

impl Transform3D {
    #[must_use]
    pub const fn identity() -> Self {
        // RDKit✔️✔️: Transform3D() : RDNumeric::SquareMatrix<double>(DIM_3D, 0.0) {
        // RDKit✔️✔️:   for (unsigned int i = 0; i < DIM_3D; i++) {
        // RDKit✔️✔️:     unsigned int id = i * (DIM_3D + 1);
        // RDKit✔️✔️:     d_data[id] = 1.0;
        // RDKit✔️✔️:   }
        Self {
            values: [
                1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            ],
        }
    }

    #[must_use]
    pub const fn values(&self) -> &[f64; 16] {
        &self.values
    }

    #[must_use]
    pub fn transform_point(&self, point: [f64; 3]) -> [f64; 3] {
        // RDKit✔️✔️: double x = data[0] * pt.x + data[1] * pt.y + data[2] * pt.z + data[3];
        // RDKit✔️✔️: double y = data[4] * pt.x + data[5] * pt.y + data[6] * pt.z + data[7];
        // RDKit✔️✔️: double z = data[8] * pt.x + data[9] * pt.y + data[10] * pt.z + data[11];
        [
            self.values[0] * point[0]
                + self.values[1] * point[1]
                + self.values[2] * point[2]
                + self.values[3],
            self.values[4] * point[0]
                + self.values[5] * point[1]
                + self.values[6] * point[2]
                + self.values[7],
            self.values[8] * point[0]
                + self.values[9] * point[1]
                + self.values[10] * point[2]
                + self.values[11],
        ]
    }

    fn rotation(angle: f64, axis: [f64; 3]) -> Self {
        // RDKit✔️✔️: double t = 1 - cosT;
        // RDKit✔️✔️: data[0] = t * X * X + cosT;
        // RDKit✔️✔️: data[1] = t * X * Y - sinT * Z;
        // RDKit✔️✔️: data[2] = t * X * Z + sinT * Y;
        // RDKit✔️✔️: data[4] = t * X * Y + sinT * Z;
        // RDKit✔️✔️: data[5] = t * Y * Y + cosT;
        // RDKit✔️✔️: data[6] = t * Y * Z - sinT * X;
        // RDKit✔️✔️: data[8] = t * X * Z - sinT * Y;
        // RDKit✔️✔️: data[9] = t * Y * Z + sinT * X;
        // RDKit✔️✔️: data[10] = t * Z * Z + cosT;
        let (sin_t, cos_t) = angle.sin_cos();
        let t = 1.0 - cos_t;
        let [x, y, z] = axis;
        let mut result = Self::identity();
        result.values[0] = t * x * x + cos_t;
        result.values[1] = t * x * y - sin_t * z;
        result.values[2] = t * x * z + sin_t * y;
        result.values[4] = t * x * y + sin_t * z;
        result.values[5] = t * y * y + cos_t;
        result.values[6] = t * y * z - sin_t * x;
        result.values[8] = t * x * z - sin_t * y;
        result.values[9] = t * y * z + sin_t * x;
        result.values[10] = t * z * z + cos_t;
        result
    }
}

impl Default for Transform3D {
    fn default() -> Self {
        Self::identity()
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct CentroidParams {
    pub ignore_hydrogens: bool,
    pub weights: Option<Vec<f64>>,
}

impl Default for CentroidParams {
    fn default() -> Self {
        Self {
            ignore_hydrogens: true,
            weights: None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PrincipalAxesKind {
    Inertia,
    Gyration,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PrincipalAxesParams {
    pub kind: PrincipalAxesKind,
    pub ignore_hydrogens: bool,
    pub weights: Option<Vec<f64>>,
}

impl Default for PrincipalAxesParams {
    fn default() -> Self {
        Self {
            kind: PrincipalAxesKind::Inertia,
            ignore_hydrogens: false,
            weights: None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PrincipalAxesAndMoments {
    pub axes: [[f64; 3]; 3],
    pub moments: [f64; 3],
}

#[derive(Debug, Clone, PartialEq)]
pub struct CanonicalTransformParams {
    pub center: Option<[f64; 3]>,
    pub normalize_covariance: bool,
    pub ignore_hydrogens: bool,
}

impl Default for CanonicalTransformParams {
    fn default() -> Self {
        Self {
            center: None,
            normalize_covariance: false,
            ignore_hydrogens: true,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct AtomPositionParams {
    pub conformer_id: Option<usize>,
}

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum TransformError {
    #[error("invalid topology: {0}")]
    InvalidTopology(#[from] TopologyValidationError),
    #[error("invalid coordinates: {0}")]
    InvalidCoordinates(#[from] CoordinateValidationError),
    #[error("atom {atom} in role {role} is out of range for {atom_count} atoms")]
    AtomOutOfRange {
        role: &'static str,
        atom: AtomId,
        atom_count: usize,
    },
    #[error("weight vector has {actual} rows, but at least {required} are required")]
    WeightCount { actual: usize, required: usize },
    #[error("atom {atom} has non-finite weight {value}")]
    NonFiniteWeight { atom: AtomId, value: f64 },
    #[error("selected weights have invalid sum {sum}")]
    InvalidWeightSum { sum: f64 },
    #[error("no atoms remain after selection")]
    NoSelectedAtoms,
    #[error("point {role} has a non-finite {axis} coordinate")]
    NonFinitePoint {
        role: &'static str,
        axis: &'static str,
    },
    #[error("transform element ({row}, {column}) is non-finite")]
    NonFiniteTransform { row: usize, column: usize },
    #[error("symmetric eigenvalue calculation did not converge")]
    EigenDecompositionDidNotConverge,
    #[error("coordinate block has no 3D conformer")]
    No3dConformer,
    #[error("3D conformer id {conformer_id} was not found")]
    ConformerNotFound { conformer_id: usize },
    #[error("atoms {first_role}={first} and {second_role}={second} are not bonded")]
    AtomsNotBonded {
        first_role: &'static str,
        first: AtomId,
        second_role: &'static str,
        second: AtomId,
    },
    #[error("ring bond {bond} cannot be moved")]
    RingBondNotMovable { bond: BondId },
    #[error("angle bonds {first_bond} and {second_bond} both belong to rings")]
    BothAngleBondsInRing {
        first_bond: BondId,
        second_bond: BondId,
    },
    #[error("atoms {first_role}={first} and {second_role}={second} have coincident coordinates")]
    CoincidentCoordinates {
        first_role: &'static str,
        first: AtomId,
        second_role: &'static str,
        second: AtomId,
    },
    #[error("rotation axis is undefined")]
    UndefinedRotationAxis,
    #[error("target {quantity} value {value} is invalid")]
    InvalidTargetValue { quantity: &'static str, value: f64 },
}

pub fn centroid(
    topology: &TopologyBlock,
    conformer: &Conformer3D,
    params: &CentroidParams,
) -> Result<[f64; 3], TransformError> {
    // RDKit✔️✔️: PRECONDITION(!weights || weights->size() >= conf.getNumAtoms(),
    // RDKit✔️✔️:              "bad weights vector");
    // RDKit✔️✔️: for (unsigned int i = 0; i < conf.getNumAtoms(); ++i) {
    // RDKit✔️✔️:   if (ignoreHs && mol.getAtomWithIdx(i)->getAtomicNum() == 1) continue;
    // RDKit✔️✔️:   double w = (weights ? weights->at(i) : 1.0);
    // RDKit✔️✔️:   wSum += w;
    // RDKit✔️✔️:   res += conf.getAtomPos(i) * w;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: res /= wSum;
    validate_pair(topology, conformer)?;
    centroid_unchecked(
        topology,
        conformer,
        params.ignore_hydrogens,
        params.weights.as_deref(),
    )
}

pub fn principal_axes_and_moments(
    topology: &TopologyBlock,
    conformer: &Conformer3D,
    params: &PrincipalAxesParams,
) -> Result<PrincipalAxesAndMoments, TransformError> {
    // RDKit✔️✔️: auto origin = computeCentroid(conf, ignoreHs, weights);
    // RDKit✔️✔️: computeInertiaTerms(conf, origin, sumXX, sumXY, sumXZ, sumYY, sumYZ, sumZZ,
    // RDKit✔️✔️:                     ignoreHs, weights);
    // RDKit✔️✔️: return getEigenValEigenVectHelper(axes, moments, sumXX, sumXY, sumXZ,
    // RDKit✔️✔️:                                   sumYY, sumYZ, sumZZ);
    // RDKit✔️✔️: bool res = getEigenValEigenVectFromCovMat(conf, axes, moments, origin,
    // RDKit✔️✔️:                                           ignoreHs, true, weights);
    validate_pair(topology, conformer)?;
    let center = centroid_unchecked(
        topology,
        conformer,
        params.ignore_hydrogens,
        params.weights.as_deref(),
    )?;
    let tensor = tensor_terms(
        topology,
        conformer,
        center,
        params.ignore_hydrogens,
        params.weights.as_deref(),
        params.kind,
        params.kind == PrincipalAxesKind::Gyration,
    )?;
    let (moments, axes) = symmetric_eigen(tensor)?;
    Ok(PrincipalAxesAndMoments { axes, moments })
}

pub fn canonical_transform(
    topology: &TopologyBlock,
    conformer: &Conformer3D,
    params: &CanonicalTransformParams,
) -> Result<Transform3D, TransformError> {
    // RDKit✔️✔️: if (!center) { origin = computeCentroid(conf, ignoreHs); }
    // RDKit✔️✔️: else { origin = (*center); }
    // RDKit✔️✔️: trans->setToIdentity();
    // RDKit✔️✔️: if (nAtms > 1) {
    // RDKit✔️✔️:   ... std::sort(eigValsSorted.begin(), eigValsSorted.end(), ... a.second > b.second);
    // RDKit✔️✔️:   trans->setVal(col, row, eigVecs(row, colSorted));
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (test < 0.0) { trans->setVal(2, i, trans->getVal(2, i) * -1); }
    // RDKit✔️✔️: origin *= -1.0;
    // RDKit✔️✔️: trans->TransformPoint(origin);
    // RDKit✔️✔️: trans->SetTranslation(origin);
    validate_pair(topology, conformer)?;
    let origin = if let Some(center) = params.center {
        validate_point(center, "center")?;
        center
    } else {
        centroid_unchecked(topology, conformer, params.ignore_hydrogens, None)?
    };
    let mut transform = Transform3D::identity();
    if conformer.coordinates().len() > 1 {
        let tensor = tensor_terms(
            topology,
            conformer,
            origin,
            params.ignore_hydrogens,
            None,
            PrincipalAxesKind::Gyration,
            params.normalize_covariance,
        )?;
        let (moments, axes) = symmetric_eigen(tensor)?;
        let mut order = [0usize, 1, 2];
        order.sort_by(|left, right| moments[*right].total_cmp(&moments[*left]));
        for row in 0..3 {
            for column in 0..3 {
                transform.values[row * 4 + column] = axes[column][order[row]];
            }
        }
    }
    if determinant3(&transform.values) < 0.0 {
        for column in 0..3 {
            transform.values[8 + column] *= -1.0;
        }
    }
    let translated = transform.transform_point([-origin[0], -origin[1], -origin[2]]);
    transform.values[3] = translated[0];
    transform.values[7] = translated[1];
    transform.values[11] = translated[2];
    validate_transform(&transform)?;
    Ok(transform)
}

pub fn transform_conformer(
    conformer: &Conformer3D,
    transform: &Transform3D,
) -> Result<Conformer3D, TransformError> {
    // RDKit✔️✔️: RDGeom::POINT3D_VECT &positions = conf.getPositions();
    // RDKit✔️✔️: for (pi = positions.begin(); pi != positions.end(); ++pi) {
    // RDKit✔️✔️:   trans.TransformPoint(*pi);
    // RDKit✔️✔️: }
    validate_transform(transform)?;
    let mut result = conformer.clone();
    for point in result.coordinates_mut() {
        validate_point(*point, "conformer")?;
        *point = transform.transform_point(*point);
        validate_point(*point, "transformed")?;
    }
    Ok(result)
}

pub fn canonicalize_conformer(
    topology: &TopologyBlock,
    conformer: &Conformer3D,
    params: &CanonicalTransformParams,
) -> Result<Conformer3D, TransformError> {
    // RDKit✔️✔️: RDGeom::Transform3D *trans =
    // RDKit✔️✔️:     computeCanonicalTransform(conf, center, normalizeCovar, ignoreHs);
    // RDKit✔️✔️: transformConformer(conf, *trans);
    let transform = canonical_transform(topology, conformer, params)?;
    transform_conformer(conformer, &transform)
}

pub fn with_atom_position(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    atom: AtomId,
    position: [f64; 3],
    params: &AtomPositionParams,
) -> Result<CoordinateBlock, TransformError> {
    // RDKit✔️✔️: URANGE_CHECK(aid, d_positions.size());
    // RDKit✔️✔️: d_positions[aid] = pos;
    // Detached conformers are fixed-length, so the source resize branch is not
    // reachable and an out-of-range id is reported instead of appending rows.
    topology.validate()?;
    coordinates.validate_for_atom_count(topology.atoms.len())?;
    check_atom(atom, topology.atoms.len(), "atom")?;
    validate_point(position, "position")?;
    let index = match params.conformer_id {
        Some(conformer_id) => coordinates
            .conformers_3d
            .iter()
            .position(|conformer| conformer.id() == conformer_id)
            .ok_or(TransformError::ConformerNotFound { conformer_id })?,
        None => coordinates
            .conformers_3d
            .iter()
            .position(Conformer3D::is_3d)
            .ok_or(TransformError::No3dConformer)?,
    };
    let mut result = coordinates.clone();
    result.conformers_3d[index].coordinates_mut()[atom.index()] = position;
    Ok(result)
}

pub fn bond_length(conformer: &Conformer3D, i: AtomId, j: AtomId) -> Result<f64, TransformError> {
    // RDKit✔️✔️: URANGE_CHECK(iAtomId, pos.size());
    // RDKit✔️✔️: URANGE_CHECK(jAtomId, pos.size());
    // RDKit✔️✔️: return (pos[iAtomId] - pos[jAtomId]).length();
    let [pi, pj] = named_points(conformer, [("i", i), ("j", j)])?;
    Ok(norm(sub(pi, pj)))
}

pub fn with_bond_length(
    topology: &TopologyBlock,
    conformer: &Conformer3D,
    i: AtomId,
    j: AtomId,
    value: f64,
) -> Result<Conformer3D, TransformError> {
    // RDKit✔️✔️: Bond *bond = mol.getBondBetweenAtoms(iAtomId, jAtomId);
    // RDKit✔️✔️: if (!bond) throw ValueErrorException("atoms i and j must be bonded");
    // RDKit✔️✔️: if (queryIsBondInRing(bond)) throw ValueErrorException("bond (i,j) must not belong to a ring");
    // RDKit✔️✔️: if (origValue <= 1.e-8) throw ValueErrorException("atoms i and j have identical 3D coordinates");
    // RDKit✔️✔️: v *= (value / origValue - 1.);
    // RDKit✔️✔️: for (unsigned int &it : alist) { pos[it] -= v; }
    validate_pair(topology, conformer)?;
    validate_target(value, "bond_length", true)?;
    let [pi, pj] = named_points(conformer, [("i", i), ("j", j)])?;
    let bond = require_bond(topology, "i", i, "j", j)?;
    if bond_is_in_ring(topology, bond) {
        return Err(TransformError::RingBondNotMovable { bond });
    }
    let vector = sub(pi, pj);
    let original = norm(vector);
    if original <= 1.0e-8 {
        return Err(coincident("i", i, "j", j));
    }
    let delta = scale(vector, value / original - 1.0);
    let mut result = conformer.clone();
    for atom in moved_atoms(topology, i, j) {
        let current = result.coordinates()[atom];
        result.coordinates_mut()[atom] = sub(current, delta);
    }
    Ok(result)
}

pub fn angle_radians(
    conformer: &Conformer3D,
    i: AtomId,
    j: AtomId,
    k: AtomId,
) -> Result<f64, TransformError> {
    // RDKit✔️✔️: RDGeom::Point3D rJI = pos[iAtomId] - pos[jAtomId];
    // RDKit✔️✔️: if (rJISqLength <= 1.e-16) throw ValueErrorException(...);
    // RDKit✔️✔️: RDGeom::Point3D rJK = pos[kAtomId] - pos[jAtomId];
    // RDKit✔️✔️: if (rJKSqLength <= 1.e-16) throw ValueErrorException(...);
    // RDKit✔️✔️: return rJI.angleTo(rJK);
    let [pi, pj, pk] = named_points(conformer, [("i", i), ("j", j), ("k", k)])?;
    angle_from_points(pi, pj, pk, i, j, k)
}

pub fn angle_degrees(
    conformer: &Conformer3D,
    i: AtomId,
    j: AtomId,
    k: AtomId,
) -> Result<f64, TransformError> {
    // RDKit✔️✔️: return (180. / M_PI * getAngleRad(conf, iAtomId, jAtomId, kAtomId));
    Ok(180.0 / PI * angle_radians(conformer, i, j, k)?)
}

pub fn with_angle_radians(
    topology: &TopologyBlock,
    conformer: &Conformer3D,
    i: AtomId,
    j: AtomId,
    k: AtomId,
    value: f64,
) -> Result<Conformer3D, TransformError> {
    // RDKit✔️✔️: Bond *bondJI = mol.getBondBetweenAtoms(jAtomId, iAtomId);
    // RDKit✔️✔️: Bond *bondJK = mol.getBondBetweenAtoms(jAtomId, kAtomId);
    // RDKit✔️✔️: if (queryIsBondInRing(bondJI) && queryIsBondInRing(bondJK)) throw ValueErrorException(...);
    // RDKit✔️✔️: value -= rJI.angleTo(rJK);
    // RDKit✔️✔️: RDGeom::Point3D rotAxisEnd = rJI.crossProduct(rJK) + pos[jAtomId];
    // RDKit✔️✔️: _toBeMovedIdxList(mol, jAtomId, kAtomId, alist);
    validate_pair(topology, conformer)?;
    validate_target(value, "angle", false)?;
    let [pi, pj, pk] = named_points(conformer, [("i", i), ("j", j), ("k", k)])?;
    let first_bond = require_bond(topology, "i", i, "j", j)?;
    let second_bond = require_bond(topology, "j", j, "k", k)?;
    if bond_is_in_ring(topology, first_bond) && bond_is_in_ring(topology, second_bond) {
        return Err(TransformError::BothAngleBondsInRing {
            first_bond,
            second_bond,
        });
    }
    let current = angle_from_points(pi, pj, pk, i, j, k)?;
    let axis = normalize(cross(sub(pi, pj), sub(pk, pj)))?;
    rotate_moved(topology, conformer, j, k, pj, axis, value - current)
}

pub fn with_angle_degrees(
    topology: &TopologyBlock,
    conformer: &Conformer3D,
    i: AtomId,
    j: AtomId,
    k: AtomId,
    value: f64,
) -> Result<Conformer3D, TransformError> {
    // RDKit✔️✔️: setAngleRad(conf, iAtomId, jAtomId, kAtomId, value / 180. * M_PI);
    with_angle_radians(topology, conformer, i, j, k, value / 180.0 * PI)
}

pub fn dihedral_radians(
    conformer: &Conformer3D,
    i: AtomId,
    j: AtomId,
    k: AtomId,
    l: AtomId,
) -> Result<f64, TransformError> {
    // RDKit✔️✔️: RDGeom::Point3D nIJK = rIJ.crossProduct(rJK);
    // RDKit✔️✔️: RDGeom::Point3D nJKL = rJK.crossProduct(rKL);
    // RDKit✔️✔️: RDGeom::Point3D m = nIJK.crossProduct(rJK);
    // RDKit✔️✔️: return -atan2(m.dotProduct(nJKL) / sqrt(nJKLSqLength * m.lengthSq()),
    // RDKit✔️✔️:               nIJK.dotProduct(nJKL) / sqrt(nIJKSqLength * nJKLSqLength));
    let [pi, pj, pk, pl] = named_points(conformer, [("i", i), ("j", j), ("k", k), ("l", l)])?;
    dihedral_from_points(pi, pj, pk, pl, i, j, k, l)
}

pub fn dihedral_degrees(
    conformer: &Conformer3D,
    i: AtomId,
    j: AtomId,
    k: AtomId,
    l: AtomId,
) -> Result<f64, TransformError> {
    // RDKit✔️✔️: return (180. / M_PI * getDihedralRad(conf, iAtomId, jAtomId, kAtomId, lAtomId));
    Ok(180.0 / PI * dihedral_radians(conformer, i, j, k, l)?)
}

pub fn with_dihedral_radians(
    topology: &TopologyBlock,
    conformer: &Conformer3D,
    i: AtomId,
    j: AtomId,
    k: AtomId,
    l: AtomId,
    value: f64,
) -> Result<Conformer3D, TransformError> {
    // RDKit✔️✔️: Bond *bondJK = mol.getBondBetweenAtoms(jAtomId, kAtomId);
    // RDKit✔️✔️: if (queryIsBondInRing(bondJK)) throw ValueErrorException("bond (j,k) must not belong to a ring");
    // RDKit✔️✔️: value -= -atan2(m.dotProduct(nJKL) / sqrt(nJKLSqLength * m.lengthSq()),
    // RDKit✔️✔️:                 nIJK.dotProduct(nJKL) / sqrt(nIJKSqLength * nJKLSqLength));
    // RDKit✔️✔️: RDGeom::Point3D rotAxis = rotAxisEnd - rotAxisBegin;
    // RDKit✔️✔️: _toBeMovedIdxList(mol, jAtomId, kAtomId, alist);
    validate_pair(topology, conformer)?;
    validate_target(value, "dihedral", false)?;
    let [pi, pj, pk, pl] = named_points(conformer, [("i", i), ("j", j), ("k", k), ("l", l)])?;
    let bond = require_bond(topology, "j", j, "k", k)?;
    if bond_is_in_ring(topology, bond) {
        return Err(TransformError::RingBondNotMovable { bond });
    }
    let current = dihedral_from_points(pi, pj, pk, pl, i, j, k, l)?;
    let axis = normalize(sub(pk, pj))?;
    rotate_moved(topology, conformer, j, k, pj, axis, value - current)
}

pub fn with_dihedral_degrees(
    topology: &TopologyBlock,
    conformer: &Conformer3D,
    i: AtomId,
    j: AtomId,
    k: AtomId,
    l: AtomId,
    value: f64,
) -> Result<Conformer3D, TransformError> {
    // RDKit✔️✔️: setDihedralRad(conf, iAtomId, jAtomId, kAtomId, lAtomId, value / 180. * M_PI);
    with_dihedral_radians(topology, conformer, i, j, k, l, value / 180.0 * PI)
}

fn validate_pair(topology: &TopologyBlock, conformer: &Conformer3D) -> Result<(), TransformError> {
    topology.validate()?;
    conformer.validate_for_atom_count(topology.atoms.len())?;
    Ok(())
}

fn validate_point(point: [f64; 3], role: &'static str) -> Result<(), TransformError> {
    if let Some(axis) = axis_for_non_finite(point) {
        return Err(TransformError::NonFinitePoint { role, axis });
    }
    Ok(())
}

fn axis_for_non_finite(point: [f64; 3]) -> Option<&'static str> {
    [("x", point[0]), ("y", point[1]), ("z", point[2])]
        .into_iter()
        .find_map(|(axis, value)| (!value.is_finite()).then_some(axis))
}

fn validate_transform(transform: &Transform3D) -> Result<(), TransformError> {
    for row in 0..4 {
        for column in 0..4 {
            if !transform.values[row * 4 + column].is_finite() {
                return Err(TransformError::NonFiniteTransform { row, column });
            }
        }
    }
    Ok(())
}

fn check_atom(atom: AtomId, atom_count: usize, role: &'static str) -> Result<(), TransformError> {
    if atom.index() >= atom_count {
        return Err(TransformError::AtomOutOfRange {
            role,
            atom,
            atom_count,
        });
    }
    Ok(())
}

fn named_points<const N: usize>(
    conformer: &Conformer3D,
    atoms: [(&'static str, AtomId); N],
) -> Result<[[f64; 3]; N], TransformError> {
    let count = conformer.coordinates().len();
    let mut points = [[0.0; 3]; N];
    for (slot, (role, atom)) in atoms.into_iter().enumerate() {
        check_atom(atom, count, role)?;
        let point = conformer.coordinates()[atom.index()];
        validate_point(point, role)?;
        points[slot] = point;
    }
    Ok(points)
}

fn selected_weight(
    topology: &TopologyBlock,
    row: usize,
    ignore_hydrogens: bool,
    weights: Option<&[f64]>,
) -> Result<Option<f64>, TransformError> {
    if ignore_hydrogens && topology.atoms[row].atomic_number() == 1 {
        return Ok(None);
    }
    let weight = weights.map_or(1.0, |values| values[row]);
    if !weight.is_finite() {
        return Err(TransformError::NonFiniteWeight {
            atom: AtomId::new(row),
            value: weight,
        });
    }
    Ok(Some(weight))
}

fn validate_weights(weights: Option<&[f64]>, rows: usize) -> Result<(), TransformError> {
    if let Some(weights) = weights
        && weights.len() < rows
    {
        return Err(TransformError::WeightCount {
            actual: weights.len(),
            required: rows,
        });
    }
    Ok(())
}

fn centroid_unchecked(
    topology: &TopologyBlock,
    conformer: &Conformer3D,
    ignore_hydrogens: bool,
    weights: Option<&[f64]>,
) -> Result<[f64; 3], TransformError> {
    validate_weights(weights, conformer.coordinates().len())?;
    let mut sum = [0.0; 3];
    let mut weight_sum = 0.0;
    let mut selected = 0usize;
    for (row, point) in conformer.coordinates().iter().copied().enumerate() {
        let Some(weight) = selected_weight(topology, row, ignore_hydrogens, weights)? else {
            continue;
        };
        selected += 1;
        weight_sum += weight;
        for axis in 0..3 {
            sum[axis] += point[axis] * weight;
        }
    }
    if selected == 0 {
        return Err(TransformError::NoSelectedAtoms);
    }
    if !weight_sum.is_finite() || weight_sum == 0.0 {
        return Err(TransformError::InvalidWeightSum { sum: weight_sum });
    }
    Ok([
        sum[0] / weight_sum,
        sum[1] / weight_sum,
        sum[2] / weight_sum,
    ])
}

fn tensor_terms(
    topology: &TopologyBlock,
    conformer: &Conformer3D,
    center: [f64; 3],
    ignore_hydrogens: bool,
    weights: Option<&[f64]>,
    kind: PrincipalAxesKind,
    normalize_covariance: bool,
) -> Result<[[f64; 3]; 3], TransformError> {
    // RDKit✔️✔️: xx += w * loc.x * loc.x; xy += w * loc.x * loc.y;
    // RDKit✔️✔️: xz += w * loc.x * loc.z; yy += w * loc.y * loc.y;
    // RDKit✔️✔️: yz += w * loc.y * loc.z; zz += w * loc.z * loc.z;
    // RDKit✔️✔️: xx += w * (loc.y * loc.y + loc.z * loc.z);
    // RDKit✔️✔️: yy += w * (loc.x * loc.x + loc.z * loc.z);
    // RDKit✔️✔️: zz += w * (loc.y * loc.y + loc.x * loc.x);
    // RDKit✔️✔️: xy -= w * loc.x * loc.y; xz -= w * loc.x * loc.z; yz -= w * loc.z * loc.y;
    validate_weights(weights, conformer.coordinates().len())?;
    let mut result = [[0.0; 3]; 3];
    let mut weight_sum = 0.0;
    let mut selected = 0usize;
    for (row, point) in conformer.coordinates().iter().copied().enumerate() {
        let Some(weight) = selected_weight(topology, row, ignore_hydrogens, weights)? else {
            continue;
        };
        selected += 1;
        weight_sum += weight;
        let [x, y, z] = sub(point, center);
        match kind {
            PrincipalAxesKind::Gyration => {
                result[0][0] += weight * x * x;
                result[0][1] += weight * x * y;
                result[0][2] += weight * x * z;
                result[1][1] += weight * y * y;
                result[1][2] += weight * y * z;
                result[2][2] += weight * z * z;
            }
            PrincipalAxesKind::Inertia => {
                result[0][0] += weight * (y * y + z * z);
                result[1][1] += weight * (x * x + z * z);
                result[2][2] += weight * (x * x + y * y);
                result[0][1] -= weight * x * y;
                result[0][2] -= weight * x * z;
                result[1][2] -= weight * y * z;
            }
        }
    }
    if selected == 0 {
        return Err(TransformError::NoSelectedAtoms);
    }
    if !weight_sum.is_finite() || weight_sum == 0.0 {
        return Err(TransformError::InvalidWeightSum { sum: weight_sum });
    }
    result[1][0] = result[0][1];
    result[2][0] = result[0][2];
    result[2][1] = result[1][2];
    if normalize_covariance {
        for row in &mut result {
            for value in row {
                *value /= weight_sum;
            }
        }
    }
    Ok(result)
}

fn symmetric_eigen(mut matrix: [[f64; 3]; 3]) -> Result<([f64; 3], [[f64; 3]; 3]), TransformError> {
    // RDKit✔️✔️: Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(mat);
    // RDKit✔️✔️: if (eigensolver.info() != Eigen::Success) return false;
    // RDKit✔️✔️: eigVecs = eigensolver.eigenvectors();
    // RDKit✔️✔️: eigVals = eigensolver.eigenvalues();
    let mut vectors = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let scale = matrix
        .iter()
        .flatten()
        .map(|value| value.abs())
        .fold(0.0, f64::max);
    let tolerance = f64::EPSILON * scale.max(1.0) * 32.0;
    let mut converged = false;
    for _ in 0..64 {
        let mut p = 0usize;
        let mut q = 1usize;
        for (left, right) in [(0, 1), (0, 2), (1, 2)] {
            if matrix[left][right].abs() > matrix[p][q].abs() {
                p = left;
                q = right;
            }
        }
        if matrix[p][q].abs() <= tolerance {
            converged = true;
            break;
        }
        let theta = (matrix[q][q] - matrix[p][p]) / (2.0 * matrix[p][q]);
        let t = if theta >= 0.0 {
            1.0 / (theta + (1.0 + theta * theta).sqrt())
        } else {
            -1.0 / (-theta + (1.0 + theta * theta).sqrt())
        };
        let c = 1.0 / (1.0 + t * t).sqrt();
        let s = t * c;
        for k in 0..3 {
            if k != p && k != q {
                let mkp = matrix[k][p];
                let mkq = matrix[k][q];
                matrix[k][p] = c * mkp - s * mkq;
                matrix[p][k] = matrix[k][p];
                matrix[k][q] = s * mkp + c * mkq;
                matrix[q][k] = matrix[k][q];
            }
        }
        let app = matrix[p][p];
        let aqq = matrix[q][q];
        let apq = matrix[p][q];
        matrix[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
        matrix[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
        matrix[p][q] = 0.0;
        matrix[q][p] = 0.0;
        for row in &mut vectors {
            let vip = row[p];
            let viq = row[q];
            row[p] = c * vip - s * viq;
            row[q] = s * vip + c * viq;
        }
    }
    if !converged {
        return Err(TransformError::EigenDecompositionDidNotConverge);
    }
    let values = [matrix[0][0], matrix[1][1], matrix[2][2]];
    let mut order = [0usize, 1, 2];
    order.sort_by(|left, right| values[*left].total_cmp(&values[*right]));
    let moments = [values[order[0]], values[order[1]], values[order[2]]];
    let mut sorted = [[0.0; 3]; 3];
    for (new_column, old_column) in order.into_iter().enumerate() {
        for row in 0..3 {
            sorted[row][new_column] = vectors[row][old_column];
        }
        let pivot = (0..3)
            .max_by(|left, right| {
                sorted[*left][new_column]
                    .abs()
                    .total_cmp(&sorted[*right][new_column].abs())
            })
            .unwrap_or(0);
        if sorted[pivot][new_column] < 0.0 {
            for row in &mut sorted {
                row[new_column] *= -1.0;
            }
        }
    }
    Ok((moments, sorted))
}

fn require_bond(
    topology: &TopologyBlock,
    first_role: &'static str,
    first: AtomId,
    second_role: &'static str,
    second: AtomId,
) -> Result<BondId, TransformError> {
    check_atom(first, topology.atoms.len(), first_role)?;
    check_atom(second, topology.atoms.len(), second_role)?;
    topology
        .adjacency
        .neighbors_of(first.index())
        .iter()
        .find(|neighbor| neighbor.atom_index == second.index())
        .map(|neighbor| neighbor.bond)
        .ok_or(TransformError::AtomsNotBonded {
            first_role,
            first,
            second_role,
            second,
        })
}

fn bond_is_in_ring(topology: &TopologyBlock, excluded: BondId) -> bool {
    let bond = &topology.bonds[excluded.index()];
    let start = bond.begin().index();
    let target = bond.end().index();
    let mut seen = vec![false; topology.atoms.len()];
    let mut stack = vec![start];
    seen[start] = true;
    while let Some(atom) = stack.pop() {
        for neighbor in topology.adjacency.neighbors_of(atom) {
            if neighbor.bond == excluded {
                continue;
            }
            if neighbor.atom_index == target {
                return true;
            }
            if !seen[neighbor.atom_index] {
                seen[neighbor.atom_index] = true;
                stack.push(neighbor.atom_index);
            }
        }
    }
    false
}

fn moved_atoms(topology: &TopologyBlock, anchor: AtomId, moving: AtomId) -> Vec<usize> {
    // RDKit✔️✔️: stack.push(jAtomId);
    // RDKit✔️✔️: visitedIdx[iAtomId] = 1;
    // RDKit✔️✔️: visitedIdx[jAtomId] = 1;
    // RDKit✔️✔️: for (unsigned int i = 0; i < nAtoms; ++i) {
    // RDKit✔️✔️:   if (visitedIdx[i] && (i != iAtomId)) alist.push_back(i);
    // RDKit✔️✔️: }
    let mut seen = vec![false; topology.atoms.len()];
    seen[anchor.index()] = true;
    seen[moving.index()] = true;
    let mut stack = vec![moving.index()];
    while let Some(atom) = stack.pop() {
        for neighbor in topology.adjacency.neighbors_of(atom) {
            if !seen[neighbor.atom_index] {
                seen[neighbor.atom_index] = true;
                stack.push(neighbor.atom_index);
            }
        }
    }
    seen.into_iter()
        .enumerate()
        .filter_map(|(atom, visited)| (visited && atom != anchor.index()).then_some(atom))
        .collect()
}

fn rotate_moved(
    topology: &TopologyBlock,
    conformer: &Conformer3D,
    anchor: AtomId,
    moving: AtomId,
    origin: [f64; 3],
    axis: [f64; 3],
    angle: f64,
) -> Result<Conformer3D, TransformError> {
    let rotation = Transform3D::rotation(angle, axis);
    let mut result = conformer.clone();
    for atom in moved_atoms(topology, anchor, moving) {
        let relative = sub(result.coordinates()[atom], origin);
        result.coordinates_mut()[atom] = add(rotation.transform_point(relative), origin);
    }
    Ok(result)
}

fn angle_from_points(
    i_point: [f64; 3],
    j_point: [f64; 3],
    k_point: [f64; 3],
    i: AtomId,
    j: AtomId,
    k: AtomId,
) -> Result<f64, TransformError> {
    let ji = sub(i_point, j_point);
    let jk = sub(k_point, j_point);
    let ji_sq = norm_squared(ji);
    let jk_sq = norm_squared(jk);
    if ji_sq <= 1.0e-16 {
        return Err(coincident("i", i, "j", j));
    }
    if jk_sq <= 1.0e-16 {
        return Err(coincident("j", j, "k", k));
    }
    Ok((dot(ji, jk) / (ji_sq * jk_sq).sqrt())
        .clamp(-1.0, 1.0)
        .acos())
}

#[allow(clippy::too_many_arguments)]
fn dihedral_from_points(
    pi: [f64; 3],
    pj: [f64; 3],
    pk: [f64; 3],
    pl: [f64; 3],
    i: AtomId,
    j: AtomId,
    k: AtomId,
    l: AtomId,
) -> Result<f64, TransformError> {
    let ij = sub(pj, pi);
    let jk = sub(pk, pj);
    let kl = sub(pl, pk);
    if norm_squared(ij) <= 1.0e-16 {
        return Err(coincident("i", i, "j", j));
    }
    if norm_squared(jk) <= 1.0e-16 {
        return Err(coincident("j", j, "k", k));
    }
    if norm_squared(kl) <= 1.0e-16 {
        return Err(coincident("k", k, "l", l));
    }
    let n_ijk = cross(ij, jk);
    let n_jkl = cross(jk, kl);
    let m = cross(n_ijk, jk);
    let n_ijk_sq = norm_squared(n_ijk);
    let n_jkl_sq = norm_squared(n_jkl);
    let m_sq = norm_squared(m);
    if n_ijk_sq == 0.0 || n_jkl_sq == 0.0 || m_sq == 0.0 {
        return Err(TransformError::UndefinedRotationAxis);
    }
    Ok(-(dot(m, n_jkl) / (n_jkl_sq * m_sq).sqrt())
        .atan2(dot(n_ijk, n_jkl) / (n_ijk_sq * n_jkl_sq).sqrt()))
}

fn validate_target(
    value: f64,
    quantity: &'static str,
    nonnegative: bool,
) -> Result<(), TransformError> {
    if !value.is_finite() || (nonnegative && value < 0.0) {
        return Err(TransformError::InvalidTargetValue { quantity, value });
    }
    Ok(())
}

fn coincident(
    first_role: &'static str,
    first: AtomId,
    second_role: &'static str,
    second: AtomId,
) -> TransformError {
    TransformError::CoincidentCoordinates {
        first_role,
        first,
        second_role,
        second,
    }
}

fn determinant3(values: &[f64; 16]) -> f64 {
    values[0] * (values[5] * values[10] - values[6] * values[9])
        - values[1] * (values[4] * values[10] - values[6] * values[8])
        + values[2] * (values[4] * values[9] - values[5] * values[8])
}

fn normalize(vector: [f64; 3]) -> Result<[f64; 3], TransformError> {
    let length = norm(vector);
    if !length.is_finite() || length == 0.0 {
        return Err(TransformError::UndefinedRotationAxis);
    }
    Ok(scale(vector, 1.0 / length))
}

fn add(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [left[0] + right[0], left[1] + right[1], left[2] + right[2]]
}

fn sub(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [left[0] - right[0], left[1] - right[1], left[2] - right[2]]
}

fn scale(vector: [f64; 3], factor: f64) -> [f64; 3] {
    [vector[0] * factor, vector[1] * factor, vector[2] * factor]
}

fn dot(left: [f64; 3], right: [f64; 3]) -> f64 {
    left[0] * right[0] + left[1] * right[1] + left[2] * right[2]
}

fn cross(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    ]
}

fn norm_squared(vector: [f64; 3]) -> f64 {
    dot(vector, vector)
}

fn norm(vector: [f64; 3]) -> f64 {
    norm_squared(vector).sqrt()
}
