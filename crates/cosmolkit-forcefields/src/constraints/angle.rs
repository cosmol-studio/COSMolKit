// Copyright (C) 2004-2024 Paolo Tosco and other RDKit contributors
//
// @@ All Rights Reserved @@
// This file is derived from RDKit ForceField/AngleConstraint.cpp.
// The contents are covered by the terms of the BSD license
// which is included in the file license.txt, found at the root
// of the RDKit source tree.

use super::{AngleIndexArgument, AngleRangeBound, ForceField, ForceFieldKernelError};

pub(super) const RAD2DEG: f64 = 180.0 / std::f64::consts::PI;

#[derive(Clone, Debug, PartialEq)]
pub(super) struct AngleConstraintContrib {
    at1_idx: u32,
    at2_idx: u32,
    at3_idx: u32,
    min_angle_deg: f64,
    max_angle_deg: f64,
    force_constant: f64,
}

impl AngleConstraintContrib {
    pub(super) fn new_absolute(
        owner: &ForceField<'_>,
        idx1: u32,
        idx2: u32,
        idx3: u32,
        min_angle_deg: f64,
        max_angle_deg: f64,
        force_constant: f64,
    ) -> Result<Self, ForceFieldKernelError> {
        // BEGIN RDKIT CPP FUNCTION AngleConstraintContrib::AngleConstraintContrib (AngleConstraint.cpp:16-34)
        // RDKit❗✔️: AngleConstraintContrib::AngleConstraintContrib(
        // RDKit❗✔️:     ForceField *owner, unsigned int idx1, unsigned int idx2, unsigned int idx3,
        // RDKit❗✔️:     double minAngleDeg, double maxAngleDeg, double forceConst) {
        // RDKit❗✔️:   PRECONDITION(owner, "bad owner");
        // RDKit❗✔️:   RANGE_CHECK(0.0, minAngleDeg, 180.0);
        // RDKit❗✔️:   RANGE_CHECK(0.0, maxAngleDeg, 180.0);
        // RDKit❗✔️:   URANGE_CHECK(idx1, owner->positions().size());
        // RDKit❗✔️:   URANGE_CHECK(idx2, owner->positions().size());
        // RDKit❗✔️:   URANGE_CHECK(idx3, owner->positions().size());
        // RDKit❗✔️:   PRECONDITION(!(minAngleDeg > maxAngleDeg),
        // RDKit❗✔️:                "minAngleDeg must be <= maxAngleDeg");
        // RDKit❗✔️:   dp_forceField = owner;
        // RDKit❗✔️:   d_at1Idx = idx1;
        // RDKit❗✔️:   d_at2Idx = idx2;
        // RDKit❗✔️:   d_at3Idx = idx3;
        // RDKit❗✔️:   d_minAngleDeg = minAngleDeg;
        // RDKit❗✔️:   d_maxAngleDeg = maxAngleDeg;
        // RDKit❗✔️:   d_forceConstant = forceConst;
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION AngleConstraintContrib::AngleConstraintContrib
        // `&ForceField` represents the source non-null owner precondition.
        // Store no owner reference; the fixed indices/parameters need no heap
        // allocation and remain valid if the forcefield moves.
        validate_angle_range(min_angle_deg, AngleRangeBound::Minimum)?;
        validate_angle_range(max_angle_deg, AngleRangeBound::Maximum)?;
        let point_count = owner.positions().len();
        validate_angle_index(idx1, point_count, AngleIndexArgument::First)?;
        validate_angle_index(idx2, point_count, AngleIndexArgument::Second)?;
        validate_angle_index(idx3, point_count, AngleIndexArgument::Third)?;
        if min_angle_deg > max_angle_deg {
            return Err(ForceFieldKernelError::AngleBoundsOrder);
        }

        Ok(Self {
            at1_idx: idx1,
            at2_idx: idx2,
            at3_idx: idx3,
            min_angle_deg,
            max_angle_deg,
            force_constant,
        })
    }

    pub(super) fn new_relative(
        owner: &ForceField<'_>,
        idx1: u32,
        idx2: u32,
        idx3: u32,
        relative: bool,
        mut min_angle_deg: f64,
        mut max_angle_deg: f64,
        force_constant: f64,
    ) -> Result<Self, ForceFieldKernelError> {
        // BEGIN RDKIT CPP FUNCTION AngleConstraintContrib::AngleConstraintContrib (AngleConstraint.cpp:36-75)
        // RDKit❗✔️: AngleConstraintContrib::AngleConstraintContrib(
        // RDKit❗✔️:     ForceField *owner, unsigned int idx1, unsigned int idx2, unsigned int idx3,
        // RDKit❗✔️:     bool relative, double minAngleDeg, double maxAngleDeg, double forceConst) {
        // RDKit❗✔️:   PRECONDITION(owner, "bad owner");
        // RDKit❗✔️:   const RDGeom::PointPtrVect &pos = owner->positions();
        // RDKit❗✔️:   URANGE_CHECK(idx1, pos.size());
        // RDKit❗✔️:   URANGE_CHECK(idx2, pos.size());
        // RDKit❗✔️:   URANGE_CHECK(idx3, pos.size());
        // RDKit❗✔️:   PRECONDITION(!(minAngleDeg > maxAngleDeg),
        // RDKit❗✔️:                "minAngleDeg must be <= maxAngleDeg");
        // RDKit❗✔️:   if (relative) {
        // RDKit❗✔️:     const RDGeom::Point3D &p1 = *((RDGeom::Point3D *)pos[idx1]);
        // RDKit❗✔️:     const RDGeom::Point3D &p2 = *((RDGeom::Point3D *)pos[idx2]);
        // RDKit❗✔️:     const RDGeom::Point3D &p3 = *((RDGeom::Point3D *)pos[idx3]);
        // RDKit❗✔️:     const RDGeom::Point3D r[2] = {p1 - p2, p3 - p2};
        // RDKit❗✔️:     const double rLengthSq[2] = {std::max(1.0e-5, r[0].lengthSq()),
        // RDKit❗✔️:                                  std::max(1.0e-5, r[1].lengthSq())};
        // RDKit❗✔️:     double cosTheta = r[0].dotProduct(r[1]) / sqrt(rLengthSq[0] * rLengthSq[1]);
        // RDKit❗✔️:     cosTheta = std::clamp(cosTheta, -1.0, 1.0);
        // RDKit❗✔️:     const double angle = RAD2DEG * acos(cosTheta);
        // RDKit❗✔️:     minAngleDeg += angle;
        // RDKit❗✔️:     maxAngleDeg += angle;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   RANGE_CHECK(0.0, minAngleDeg, 180.0);
        // RDKit❗✔️:   RANGE_CHECK(0.0, maxAngleDeg, 180.0);
        // RDKit❗✔️:   dp_forceField = owner;
        // RDKit❗✔️:   d_at1Idx = idx1;
        // RDKit❗✔️:   d_at2Idx = idx2;
        // RDKit❗✔️:   d_at3Idx = idx3;
        // RDKit❗✔️:   d_minAngleDeg = minAngleDeg;
        // RDKit❗✔️:   d_maxAngleDeg = maxAngleDeg;
        // RDKit❗✔️:   d_forceConstant = forceConst;
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION AngleConstraintContrib::AngleConstraintContrib
        let positions = owner.positions();
        let point_count = positions.len();
        validate_angle_index(idx1, point_count, AngleIndexArgument::First)?;
        validate_angle_index(idx2, point_count, AngleIndexArgument::Second)?;
        validate_angle_index(idx3, point_count, AngleIndexArgument::Third)?;
        if min_angle_deg > max_angle_deg {
            return Err(ForceFieldKernelError::AngleBoundsOrder);
        }
        if relative {
            let p1 = &positions[idx1 as usize][..];
            let p2 = &positions[idx2 as usize][..];
            let p3 = &positions[idx3 as usize][..];
            let angle = angle_degrees(p1, p2, p3);
            min_angle_deg += angle;
            max_angle_deg += angle;
        }
        validate_angle_range(min_angle_deg, AngleRangeBound::Minimum)?;
        validate_angle_range(max_angle_deg, AngleRangeBound::Maximum)?;

        Ok(Self {
            at1_idx: idx1,
            at2_idx: idx2,
            at3_idx: idx3,
            min_angle_deg,
            max_angle_deg,
            force_constant,
        })
    }

    pub(super) fn compute_angle_term(&self, angle: f64) -> f64 {
        // BEGIN RDKIT CPP FUNCTION AngleConstraintContrib::computeAngleTerm (AngleConstraint.cpp:77-85)
        // RDKit❗✔️: double AngleConstraintContrib::computeAngleTerm(double angle) const {
        // RDKit❗✔️:   double angleTerm = 0.0;
        // RDKit❗✔️:   if (angle < d_minAngleDeg) {
        // RDKit❗✔️:     angleTerm = angle - d_minAngleDeg;
        // RDKit❗✔️:   } else if (angle > d_maxAngleDeg) {
        // RDKit❗✔️:     angleTerm = angle - d_maxAngleDeg;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return angleTerm;
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION AngleConstraintContrib::computeAngleTerm
        let mut angle_term = 0.0;
        if angle < self.min_angle_deg {
            angle_term = angle - self.min_angle_deg;
        } else if angle > self.max_angle_deg {
            angle_term = angle - self.max_angle_deg;
        }
        angle_term
    }

    pub(super) fn get_energy(&self, pos: &[f64]) -> f64 {
        // BEGIN RDKIT CPP FUNCTION AngleConstraintContrib::getEnergy (AngleConstraint.cpp:87-108)
        // RDKit❗✔️: double AngleConstraintContrib::getEnergy(double *pos) const {
        // RDKit❗✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit❗✔️:   PRECONDITION(pos, "bad vector");
        // RDKit❗✔️:   const RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
        // RDKit❗✔️:                            pos[3 * d_at1Idx + 2]);
        // RDKit❗✔️:   const RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
        // RDKit❗✔️:                            pos[3 * d_at2Idx + 2]);
        // RDKit❗✔️:   const RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
        // RDKit❗✔️:                            pos[3 * d_at3Idx + 2]);
        // RDKit❗✔️:   const RDGeom::Point3D r[2] = {p1 - p2, p3 - p2};
        // RDKit❗✔️:   const double rLengthSq[2] = {std::max(1.0e-5, r[0].lengthSq()),
        // RDKit❗✔️:                                std::max(1.0e-5, r[1].lengthSq())};
        // RDKit❗✔️:   double cosTheta = r[0].dotProduct(r[1]) / sqrt(rLengthSq[0] * rLengthSq[1]);
        // RDKit❗✔️:   cosTheta = std::clamp(cosTheta, -1.0, 1.0);
        // RDKit❗✔️:   const double angle = RAD2DEG * acos(cosTheta);
        // RDKit❗✔️:   const double angleTerm = computeAngleTerm(angle);
        // RDKit❗✔️:   return d_forceConstant * angleTerm * angleTerm;
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION AngleConstraintContrib::getEnergy
        let start1 = 3 * self.at1_idx as usize;
        let start2 = 3 * self.at2_idx as usize;
        let start3 = 3 * self.at3_idx as usize;
        let p1 = &pos[start1..start1 + 3];
        let p2 = &pos[start2..start2 + 3];
        let p3 = &pos[start3..start3 + 3];
        let angle = angle_degrees(p1, p2, p3);
        let angle_term = self.compute_angle_term(angle);
        self.force_constant * angle_term * angle_term
    }

    pub(super) fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
        // BEGIN RDKIT CPP FUNCTION AngleConstraintContrib::getGrad (AngleConstraint.cpp:106-141)
        // RDKit❗✔️: void AngleConstraintContrib::getGrad(double *pos, double *grad) const {
        // RDKit❗✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit❗✔️:   PRECONDITION(pos, "bad vector");
        // RDKit❗✔️:   PRECONDITION(grad, "bad vector");
        // RDKit❗✔️:
        // RDKit❗✔️:   const RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
        // RDKit❗✔️:                            pos[3 * d_at1Idx + 2]);
        // RDKit❗✔️:   const RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
        // RDKit❗✔️:                            pos[3 * d_at2Idx + 2]);
        // RDKit❗✔️:   const RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
        // RDKit❗✔️:                            pos[3 * d_at3Idx + 2]);
        // RDKit❗✔️:   const RDGeom::Point3D r[2] = {p1 - p2, p3 - p2};
        // RDKit❗✔️:   const double rLengthSq[2] = {std::max(1.0e-5, r[0].lengthSq()),
        // RDKit❗✔️:                                std::max(1.0e-5, r[1].lengthSq())};
        // RDKit❗✔️:   double cosTheta = r[0].dotProduct(r[1]) / sqrt(rLengthSq[0] * rLengthSq[1]);
        // RDKit❗✔️:   cosTheta = std::clamp(cosTheta, -1.0, 1.0);
        // RDKit❗✔️:   const double angle = RAD2DEG * acos(cosTheta);
        // RDKit❗✔️:   const double angleTerm = computeAngleTerm(angle);
        // RDKit❗✔️:
        // RDKit❗✔️:   double dE_dTheta = 2.0 * RAD2DEG * d_forceConstant * angleTerm;
        // RDKit❗✔️:
        // RDKit❗✔️:   RDGeom::Point3D rp = r[1].crossProduct(r[0]);
        // RDKit❗✔️:   double prefactor = dE_dTheta / std::max(1.0e-5, rp.length());
        // RDKit❗✔️:   double t[2] = {-prefactor / rLengthSq[0], prefactor / rLengthSq[1]};
        // RDKit❗✔️:   RDGeom::Point3D dedp[3];
        // RDKit❗✔️:   dedp[0] = r[0].crossProduct(rp) * t[0];
        // RDKit❗✔️:   dedp[2] = r[1].crossProduct(rp) * t[1];
        // RDKit❗✔️:   dedp[1] = -dedp[0] - dedp[2];
        // RDKit❗✔️:   double *g[3] = {&(grad[3 * d_at1Idx]), &(grad[3 * d_at2Idx]),
        // RDKit❗✔️:                   &(grad[3 * d_at3Idx])};
        // RDKit❗✔️:   for (unsigned int i = 0; i < 3; ++i) {
        // RDKit❗✔️:     g[i][0] += dedp[i].x;
        // RDKit❗✔️:     g[i][1] += dedp[i].y;
        // RDKit❗✔️:     g[i][2] += dedp[i].z;
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION AngleConstraintContrib::getGrad
        // Construction has already validated the indices against the borrowed
        // owner positions; no owner pointer is retained. Rust slices provide
        // non-null position and gradient references.
        let start1 = 3 * self.at1_idx as usize;
        let start2 = 3 * self.at2_idx as usize;
        let start3 = 3 * self.at3_idx as usize;
        let p1 = [pos[start1], pos[start1 + 1], pos[start1 + 2]];
        let p2 = [pos[start2], pos[start2 + 1], pos[start2 + 2]];
        let p3 = [pos[start3], pos[start3 + 1], pos[start3 + 2]];
        let r = [
            angle_vector_difference(&p1, &p2),
            angle_vector_difference(&p3, &p2),
        ];
        let r0_length_sq = angle_vector_length_sq(&r[0]);
        let r1_length_sq = angle_vector_length_sq(&r[1]);
        let r_length_sq = [
            if 1.0e-5 < r0_length_sq {
                r0_length_sq
            } else {
                1.0e-5
            },
            if 1.0e-5 < r1_length_sq {
                r1_length_sq
            } else {
                1.0e-5
            },
        ];
        let mut cos_theta =
            angle_vector_dot_product(&r[0], &r[1]) / (r_length_sq[0] * r_length_sq[1]).sqrt();
        if cos_theta < -1.0 {
            cos_theta = -1.0;
        } else if cos_theta > 1.0 {
            cos_theta = 1.0;
        }
        let angle = RAD2DEG * cos_theta.acos();
        let angle_term = self.compute_angle_term(angle);

        let d_e_d_theta = 2.0 * RAD2DEG * self.force_constant * angle_term;
        let rp = angle_vector_cross_product(&r[1], &r[0]);
        let rp_length = angle_vector_length(&rp);
        let prefactor = d_e_d_theta
            / if 1.0e-5 < rp_length {
                rp_length
            } else {
                1.0e-5
            };
        let t = [-prefactor / r_length_sq[0], prefactor / r_length_sq[1]];
        let dedp0 = angle_vector_scale(&angle_vector_cross_product(&r[0], &rp), t[0]);
        let dedp2 = angle_vector_scale(&angle_vector_cross_product(&r[1], &rp), t[1]);
        let negative_dedp0 = angle_vector_negated(&dedp0);
        let dedp1 = angle_vector_difference(&negative_dedp0, &dedp2);

        for (atom_idx, delta) in [
            (self.at1_idx, dedp0),
            (self.at2_idx, dedp1),
            (self.at3_idx, dedp2),
        ] {
            let start = 3 * atom_idx as usize;
            grad[start] += delta[0];
            grad[start + 1] += delta[1];
            grad[start + 2] += delta[2];
        }
    }
}

pub(super) fn angle_vector_difference(left: &[f64; 3], right: &[f64; 3]) -> [f64; 3] {
    // BEGIN RDKIT CPP HELPER RDGeom::operator- (Geometry/point.cpp:64-70)
    // RDKit❗✔️: Point3D operator-(const Point3D &p1, const Point3D &p2) {
    // RDKit❗✔️:   Point3D res;
    // RDKit❗✔️:   res.x = p1.x - p2.x;
    // RDKit❗✔️:   res.y = p1.y - p2.y;
    // RDKit❗✔️:   res.z = p1.z - p2.z;
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP HELPER RDGeom::operator-
    [left[0] - right[0], left[1] - right[1], left[2] - right[2]]
}

pub(super) fn angle_vector_length_sq(vector: &[f64; 3]) -> f64 {
    // BEGIN RDKIT CPP HELPER RDGeom::Point3D::lengthSq (Geometry/point.h:163-167)
    // RDKit❗✔️: constexpr double lengthSq() const override {
    // RDKit❗✔️:   // double res = pow(x,2) + pow(y,2) + pow(z,2);
    // RDKit❗✔️:   double res = x * x + y * y + z * z;
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP HELPER RDGeom::Point3D::lengthSq
    vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]
}

pub(super) fn angle_vector_dot_product(left: &[f64; 3], right: &[f64; 3]) -> f64 {
    // BEGIN RDKIT CPP HELPER RDGeom::Point3D::dotProduct (Geometry/point.h:169-172)
    // RDKit❗✔️: constexpr double dotProduct(const Point3D &other) const {
    // RDKit❗✔️:   double res = x * (other.x) + y * (other.y) + z * (other.z);
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP HELPER RDGeom::Point3D::dotProduct
    left[0] * right[0] + left[1] * right[1] + left[2] * right[2]
}

pub(super) fn angle_vector_cross_product(left: &[f64; 3], right: &[f64; 3]) -> [f64; 3] {
    // BEGIN RDKIT CPP HELPER RDGeom::Point3D::crossProduct (Geometry/point.h:228-234)
    // RDKit❗✔️: constexpr Point3D crossProduct(const Point3D &other) const {
    // RDKit❗✔️:   Point3D res;
    // RDKit❗✔️:   res.x = y * (other.z) - z * (other.y);
    // RDKit❗✔️:   res.y = -x * (other.z) + z * (other.x);
    // RDKit❗✔️:   res.z = x * (other.y) - y * (other.x);
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP HELPER RDGeom::Point3D::crossProduct
    [
        left[1] * right[2] - left[2] * right[1],
        -left[0] * right[2] + left[2] * right[0],
        left[0] * right[1] - left[1] * right[0],
    ]
}

pub(super) fn angle_vector_length(vector: &[f64; 3]) -> f64 {
    // BEGIN RDKIT CPP HELPER RDGeom::Point3D::length (Geometry/point.h:158-161)
    // RDKit❗✔️: double length() const override {
    // RDKit❗✔️:   double res = x * x + y * y + z * z;
    // RDKit❗✔️:   return sqrt(res);
    // RDKit❗✔️: }
    // END RDKIT CPP HELPER RDGeom::Point3D::length
    let res = vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2];
    res.sqrt()
}

pub(super) fn angle_vector_scale(vector: &[f64; 3], scale: f64) -> [f64; 3] {
    // BEGIN RDKIT CPP HELPER RDGeom::operator* (Geometry/point.cpp:72-77)
    // RDKit❗✔️: Point3D operator*(const Point3D &p1, double v) {
    // RDKit❗✔️:   Point3D res;
    // RDKit❗✔️:   res.x = p1.x * v;
    // RDKit❗✔️:   res.y = p1.y * v;
    // RDKit❗✔️:   res.z = p1.z * v;
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP HELPER RDGeom::operator*
    [vector[0] * scale, vector[1] * scale, vector[2] * scale]
}

pub(super) fn angle_vector_negated(vector: &[f64; 3]) -> [f64; 3] {
    // BEGIN RDKIT CPP HELPER RDGeom::Point3D::operator- (Geometry/point.h:139-145)
    // RDKit❗✔️: constexpr Point3D operator-() const {
    // RDKit❗✔️:   Point3D res(x, y, z);
    // RDKit❗✔️:   res.x *= -1.0;
    // RDKit❗✔️:   res.y *= -1.0;
    // RDKit❗✔️:   res.z *= -1.0;
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP HELPER RDGeom::Point3D::operator-
    [vector[0] * -1.0, vector[1] * -1.0, vector[2] * -1.0]
}

pub(super) fn validate_angle_index(
    index: u32,
    upper_bound: usize,
    argument: AngleIndexArgument,
) -> Result<(), ForceFieldKernelError> {
    // BEGIN RDKIT CPP HELPER URANGE_CHECK (RDGeneral/Invariant.h:141-151)
    // RDKit❗✔️: #define URANGE_CHECK(x, hi)                                                 \
    // RDKit❗✔️:   if (x >= (hi)) {                                                          \
    // RDKit❗✔️:     std::stringstream errstr;                                               \
    // RDKit❗✔️:     errstr << x << " < " << hi;                                             \
    // RDKit❗✔️:     Invar::Invariant inv("Range Error", #x, errstr.str().c_str(), __FILE__, \
    // RDKit❗✔️:                          __LINE__);                                         \
    // RDKit❗✔️:     BOOST_LOG(rdErrorLog) << "\n\n****\n" << inv << "****\n" << std::endl;  \
    // RDKit❗✔️:     throw inv;                                                              \
    // RDKit❗✔️:   }
    // END RDKIT CPP HELPER URANGE_CHECK
    if index as usize >= upper_bound {
        return Err(ForceFieldKernelError::AngleIndexOutOfRange {
            argument,
            index,
            upper_bound,
        });
    }
    Ok(())
}

pub(super) fn validate_angle_range(
    angle: f64,
    bound: AngleRangeBound,
) -> Result<(), ForceFieldKernelError> {
    // BEGIN RDKIT CPP HELPER RANGE_CHECK (RDGeneral/Invariant.h:131-139)
    // RDKit❗✔️: #define RANGE_CHECK(lo, x, hi)                                              \
    // RDKit❗✔️:   if ((lo) > (hi) || (x) < (lo) || (x) > (hi)) {                          \
    // RDKit❗✔️:     std::stringstream errstr;                                               \
    // RDKit❗✔️:     errstr << lo << " <= " << x << " <= " << hi;                           \
    // RDKit❗✔️:     Invar::Invariant inv("Range Error", #x, errstr.str().c_str(), __FILE__, \
    // RDKit❗✔️:                          __LINE__);                                         \
    // RDKit❗✔️:     BOOST_LOG(rdErrorLog) << "\n\n****\n" << inv << "****\n" << std::endl;  \
    // RDKit❗✔️:     throw inv;                                                              \
    // RDKit❗✔️:   }
    // END RDKIT CPP HELPER RANGE_CHECK
    if 0.0_f64 > 180.0_f64 || angle < 0.0 || angle > 180.0 {
        return Err(ForceFieldKernelError::AngleOutOfRange { bound });
    }
    Ok(())
}

pub(super) fn angle_degrees(p1: &[f64], p2: &[f64], p3: &[f64]) -> f64 {
    // BEGIN RDKIT CPP CONSTANT AngleConstraint.cpp (AngleConstraint.cpp:18-20)
    // RDKit❗✔️: #ifndef M_PI
    // RDKit❗✔️: #define M_PI 3.14159265358979323846
    // RDKit❗✔️: #endif
    // RDKit❗✔️: constexpr double RAD2DEG = 180.0 / M_PI;
    // END RDKIT CPP CONSTANT AngleConstraint.cpp
    // BEGIN RDKIT CPP HELPER RDGeom::operator- (Geometry/point.cpp:64-70)
    // RDKit❗✔️: Point3D operator-(const Point3D &p1, const Point3D &p2) {
    // RDKit❗✔️:   Point3D res;
    // RDKit❗✔️:   res.x = p1.x - p2.x;
    // RDKit❗✔️:   res.y = p1.y - p2.y;
    // RDKit❗✔️:   res.z = p1.z - p2.z;
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP HELPER RDGeom::operator-
    let r0 = [p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]];
    let r1 = [p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2]];

    // BEGIN RDKIT CPP HELPER RDGeom::Point3D::lengthSq (Geometry/point.h:163-167)
    // RDKit❗✔️: constexpr double lengthSq() const override {
    // RDKit❗✔️:   // double res = pow(x,2) + pow(y,2) + pow(z,2);
    // RDKit❗✔️:   double res = x * x + y * y + z * z;
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP HELPER RDGeom::Point3D::lengthSq
    let r0_length_sq = r0[0] * r0[0] + r0[1] * r0[1] + r0[2] * r0[2];
    let r1_length_sq = r1[0] * r1[0] + r1[1] * r1[1] + r1[2] * r1[2];
    let r_length_sq = [
        if 1.0e-5 < r0_length_sq {
            r0_length_sq
        } else {
            1.0e-5
        },
        if 1.0e-5 < r1_length_sq {
            r1_length_sq
        } else {
            1.0e-5
        },
    ];

    // BEGIN RDKIT CPP HELPER RDGeom::Point3D::dotProduct (Geometry/point.h:169-172)
    // RDKit❗✔️: constexpr double dotProduct(const Point3D &other) const {
    // RDKit❗✔️:   double res = x * (other.x) + y * (other.y) + z * (other.z);
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP HELPER RDGeom::Point3D::dotProduct
    let dot_product = r0[0] * r1[0] + r0[1] * r1[1] + r0[2] * r1[2];
    let mut cos_theta = dot_product / (r_length_sq[0] * r_length_sq[1]).sqrt();
    // std::max(1.0e-5, value) and std::clamp retain the source's comparison
    // behavior for NaN; these fixed-size temporaries allocate no heap memory.
    if cos_theta < -1.0 {
        cos_theta = -1.0;
    } else if cos_theta > 1.0 {
        cos_theta = 1.0;
    }
    RAD2DEG * cos_theta.acos()
}

#[cfg(test)]
mod tests {
    use super::{
        AngleConstraintContrib, AngleIndexArgument, AngleRangeBound, ForceField,
        ForceFieldKernelError,
    };

    fn cf3d_f18_force_field<'a>(
        first: &'a mut Vec<f64>,
        middle: &'a mut Vec<f64>,
        third: &'a mut Vec<f64>,
    ) -> ForceField<'a> {
        let mut force_field = ForceField::default();
        force_field.positions_mut().extend([
            first.as_mut_slice(),
            middle.as_mut_slice(),
            third.as_mut_slice(),
        ]);
        force_field
    }

    fn cf3d_f18_flat_third(third: [f64; 3]) -> [f64; 9] {
        [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, third[0], third[1], third[2]]
    }

    fn cf3d_f18_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 1.0e-9,
            "actual={actual}, expected={expected}"
        );
    }

    fn cf3d_f19_assert_gradient(actual: &[f64; 9], expected: [f64; 9]) {
        for index in 0..9 {
            cf3d_f18_close(actual[index], expected[index]);
        }
    }

    #[test]
    fn cf3d_f18_absolute_constructor_preserves_range_index_and_order_checks() {
        // RDKit AngleConstraint.cpp:16-34 and Invariant.h RANGE_CHECK/
        // URANGE_CHECK preserve inclusive bounds and this exact error order.
        let mut first = vec![1.0, 0.0, 0.0];
        let mut middle = vec![0.0, 0.0, 0.0];
        let mut third = vec![0.0, 1.0, 0.0];
        let force_field = cf3d_f18_force_field(&mut first, &mut middle, &mut third);

        for (minimum, maximum) in [(0.0, 0.0), (180.0, 180.0), (0.0, 180.0)] {
            let contribution =
                AngleConstraintContrib::new_absolute(&force_field, 0, 1, 2, minimum, maximum, 2.0)
                    .unwrap();
            assert_eq!(contribution.min_angle_deg, minimum);
            assert_eq!(contribution.max_angle_deg, maximum);
        }

        let minimum_error =
            AngleConstraintContrib::new_absolute(&force_field, 3, 3, 3, -0.1, 180.1, 2.0)
                .unwrap_err();
        assert_eq!(
            minimum_error,
            ForceFieldKernelError::AngleOutOfRange {
                bound: AngleRangeBound::Minimum,
            }
        );
        assert_eq!(minimum_error.source_category(), "Range Error");
        assert_eq!(minimum_error.source_message(), "minAngleDeg");
        assert_eq!(
            minimum_error.source_expression().as_deref(),
            Some("minAngleDeg")
        );

        let maximum_error =
            AngleConstraintContrib::new_absolute(&force_field, 3, 3, 3, 0.0, 180.1, 2.0)
                .unwrap_err();
        assert_eq!(
            maximum_error,
            ForceFieldKernelError::AngleOutOfRange {
                bound: AngleRangeBound::Maximum,
            }
        );
        assert_eq!(maximum_error.source_message(), "maxAngleDeg");
        assert_eq!(
            maximum_error.source_expression().as_deref(),
            Some("maxAngleDeg")
        );

        let first_index_error =
            AngleConstraintContrib::new_absolute(&force_field, 3, 3, 3, 0.0, 180.0, 2.0)
                .unwrap_err();
        assert_eq!(
            first_index_error,
            ForceFieldKernelError::AngleIndexOutOfRange {
                argument: AngleIndexArgument::First,
                index: 3,
                upper_bound: 3,
            }
        );
        assert_eq!(first_index_error.source_category(), "Range Error");
        assert_eq!(first_index_error.source_message(), "idx1");

        let second_index_error =
            AngleConstraintContrib::new_absolute(&force_field, 0, 3, 3, 0.0, 180.0, 2.0)
                .unwrap_err();
        assert_eq!(second_index_error.source_message(), "idx2");

        let third_index_error =
            AngleConstraintContrib::new_absolute(&force_field, 0, 1, 3, 0.0, 180.0, 2.0)
                .unwrap_err();
        assert_eq!(third_index_error.source_message(), "idx3");

        let index_precedes_order =
            AngleConstraintContrib::new_absolute(&force_field, 3, 3, 3, 90.0, 80.0, 2.0)
                .unwrap_err();
        assert_eq!(index_precedes_order, first_index_error);

        let order_error =
            AngleConstraintContrib::new_absolute(&force_field, 0, 1, 2, 90.0, 80.0, 2.0)
                .unwrap_err();
        assert_eq!(order_error, ForceFieldKernelError::AngleBoundsOrder);
        assert_eq!(order_error.source_category(), "Pre-condition Violation");
        assert_eq!(
            order_error.source_message(),
            "minAngleDeg must be <= maxAngleDeg"
        );
        assert_eq!(
            order_error.source_expression().as_deref(),
            Some("!(minAngleDeg > maxAngleDeg)")
        );
    }

    #[test]
    fn cf3d_f18_relative_constructor_offsets_and_checks_adjusted_bounds() {
        // RDKit AngleConstraint.cpp:36-75 offsets each endpoint only for true,
        // after the source interval precondition and before RANGE_CHECK.
        let mut first = vec![1.0, 0.0, 0.0];
        let mut middle = vec![0.0, 0.0, 0.0];
        let mut third = vec![0.0, 1.0, 0.0];
        let force_field = cf3d_f18_force_field(&mut first, &mut middle, &mut third);

        let relative =
            AngleConstraintContrib::new_relative(&force_field, 0, 1, 2, true, -10.0, 10.0, 2.0)
                .unwrap();
        cf3d_f18_close(relative.min_angle_deg, 80.0);
        cf3d_f18_close(relative.max_angle_deg, 100.0);

        let relative_inclusive =
            AngleConstraintContrib::new_relative(&force_field, 0, 1, 2, true, -90.0, -80.0, 2.0)
                .unwrap();
        cf3d_f18_close(relative_inclusive.min_angle_deg, 0.0);
        cf3d_f18_close(relative_inclusive.max_angle_deg, 10.0);

        let absolute =
            AngleConstraintContrib::new_relative(&force_field, 0, 1, 2, false, 10.0, 20.0, 2.0)
                .unwrap();
        assert_eq!(absolute.min_angle_deg, 10.0);
        assert_eq!(absolute.max_angle_deg, 20.0);

        let relative_minimum_error =
            AngleConstraintContrib::new_relative(&force_field, 0, 1, 2, true, 100.0, 110.0, 2.0)
                .unwrap_err();
        assert_eq!(
            relative_minimum_error,
            ForceFieldKernelError::AngleOutOfRange {
                bound: AngleRangeBound::Minimum,
            }
        );

        let relative_maximum_error =
            AngleConstraintContrib::new_relative(&force_field, 0, 1, 2, true, -10.0, 100.0, 2.0)
                .unwrap_err();
        assert_eq!(
            relative_maximum_error,
            ForceFieldKernelError::AngleOutOfRange {
                bound: AngleRangeBound::Maximum,
            }
        );

        let relative_index_precedes_order =
            AngleConstraintContrib::new_relative(&force_field, 3, 3, 3, true, 90.0, 80.0, 2.0)
                .unwrap_err();
        assert_eq!(relative_index_precedes_order.source_message(), "idx1");

        let relative_order_error =
            AngleConstraintContrib::new_relative(&force_field, 0, 1, 2, true, 90.0, 80.0, 2.0)
                .unwrap_err();
        assert_eq!(
            relative_order_error,
            ForceFieldKernelError::AngleBoundsOrder
        );
    }

    #[test]
    fn cf3d_f18_nan_comparisons_and_nan_relative_geometry_match_source() {
        // RANGE_CHECK uses ordered comparisons and the interval check is
        // !(min > max); both therefore accept NaN bounds as pinned.
        let mut first = vec![1.0, 0.0, 0.0];
        let mut middle = vec![0.0, 0.0, 0.0];
        let mut third = vec![0.0, 1.0, 0.0];
        let force_field = cf3d_f18_force_field(&mut first, &mut middle, &mut third);
        let nan_minimum =
            AngleConstraintContrib::new_absolute(&force_field, 0, 1, 2, f64::NAN, 90.0, 2.0)
                .unwrap();
        assert!(nan_minimum.min_angle_deg.is_nan());
        let nan_maximum =
            AngleConstraintContrib::new_relative(&force_field, 0, 1, 2, false, 10.0, f64::NAN, 2.0)
                .unwrap();
        assert!(nan_maximum.max_angle_deg.is_nan());

        let mut nan_first = vec![f64::NAN, 0.0, 0.0];
        let mut nan_middle = vec![0.0, 0.0, 0.0];
        let mut nan_third = vec![0.0, 1.0, 0.0];
        let nan_force_field = cf3d_f18_force_field(&mut nan_first, &mut nan_middle, &mut nan_third);
        let relative_nan =
            AngleConstraintContrib::new_relative(&nan_force_field, 0, 1, 2, true, -10.0, 10.0, 2.0)
                .unwrap();
        assert!(relative_nan.min_angle_deg.is_nan());
        assert!(relative_nan.max_angle_deg.is_nan());

        let order_precedes_nan_offset =
            AngleConstraintContrib::new_relative(&nan_force_field, 0, 1, 2, true, 20.0, 10.0, 2.0)
                .unwrap_err();
        assert_eq!(
            order_precedes_nan_offset,
            ForceFieldKernelError::AngleBoundsOrder
        );
    }

    #[test]
    fn cf3d_f18_energy_uses_strict_flat_bottom_and_trial_coordinates() {
        // RDKit AngleConstraint.cpp:77-108 selects strict violations and
        // evaluates the explicit flat coordinate vector, not owner positions.
        let mut first = vec![1.0, 0.0, 0.0];
        let mut middle = vec![0.0, 0.0, 0.0];
        let mut third = vec![1.0, 0.0, 0.0];
        let force_field = cf3d_f18_force_field(&mut first, &mut middle, &mut third);
        let contribution =
            AngleConstraintContrib::new_absolute(&force_field, 0, 1, 2, 80.0, 100.0, 2.0).unwrap();

        assert_eq!(contribution.compute_angle_term(79.0), -1.0);
        assert_eq!(contribution.compute_angle_term(80.0), 0.0);
        assert_eq!(contribution.compute_angle_term(90.0), 0.0);
        assert_eq!(contribution.compute_angle_term(100.0), 0.0);
        assert_eq!(contribution.compute_angle_term(101.0), 1.0);
        assert_eq!(
            contribution.get_energy(&cf3d_f18_flat_third([0.0, 1.0, 0.0])),
            0.0
        );
        cf3d_f18_close(
            contribution.get_energy(&cf3d_f18_flat_third([1.0, 0.0, 0.0])),
            12_800.0,
        );
        cf3d_f18_close(
            contribution.get_energy(&cf3d_f18_flat_third([-1.0, 0.0, 0.0])),
            12_800.0,
        );

        let point =
            AngleConstraintContrib::new_absolute(&force_field, 0, 1, 2, 90.0, 90.0, 2.0).unwrap();
        cf3d_f18_close(
            point.get_energy(&cf3d_f18_flat_third([1.0, 0.0, 0.0])),
            16_200.0,
        );
        assert_eq!(point.get_energy(&cf3d_f18_flat_third([0.0, 1.0, 0.0])), 0.0);
        cf3d_f18_close(
            point.get_energy(&cf3d_f18_flat_third([-1.0, 0.0, 0.0])),
            16_200.0,
        );

        // AngleConstraint::getEnergy does not require forcefield initialization
        // or read/reset its distance cache.
        assert!(!force_field.initialized);
        assert!(force_field.distance_matrix.is_empty());
    }

    #[test]
    fn cf3d_f18_zero_collinear_and_nan_geometry_keep_source_floors() {
        // Both zero vectors have squared length floored at 1e-5, producing a
        // zero cosine and a 90 degree source angle without allocating scratch.
        let mut first = vec![0.0, 0.0, 0.0];
        let mut middle = vec![0.0, 0.0, 0.0];
        let mut third = vec![0.0, 0.0, 0.0];
        let force_field = cf3d_f18_force_field(&mut first, &mut middle, &mut third);
        let contribution =
            AngleConstraintContrib::new_absolute(&force_field, 0, 1, 2, 100.0, 120.0, 2.0).unwrap();
        cf3d_f18_close(contribution.get_energy(&[0.0; 9]), 200.0);

        let nan_energy = contribution.get_energy(&cf3d_f18_flat_third([f64::NAN, 0.0, 0.0]));
        assert_eq!(nan_energy, 0.0);
    }

    #[test]
    fn cf3d_f19_get_grad_covers_both_active_flat_bottom_sides() {
        // AngleConstraint.cpp:120-139. With a right angle and unit force,
        // 2 * RAD2DEG * 10 degrees gives this source-shaped gradient magnitude.
        const MAGNITUDE: f64 = 1145.9155902616465;
        let mut first = vec![0.0; 3];
        let mut middle = vec![0.0; 3];
        let mut third = vec![0.0; 3];
        let force_field = cf3d_f18_force_field(&mut first, &mut middle, &mut third);
        let lower =
            AngleConstraintContrib::new_absolute(&force_field, 0, 1, 2, 100.0, 120.0, 1.0).unwrap();
        let upper =
            AngleConstraintContrib::new_absolute(&force_field, 0, 1, 2, 60.0, 80.0, 1.0).unwrap();
        let trial = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let mut lower_gradient = [0.0; 9];
        let mut upper_gradient = [0.0; 9];

        lower.get_grad(&trial, &mut lower_gradient);
        upper.get_grad(&trial, &mut upper_gradient);

        cf3d_f19_assert_gradient(
            &lower_gradient,
            [
                0.0, MAGNITUDE, 0.0, -MAGNITUDE, -MAGNITUDE, 0.0, MAGNITUDE, 0.0, 0.0,
            ],
        );
        cf3d_f19_assert_gradient(
            &upper_gradient,
            [
                0.0, -MAGNITUDE, 0.0, MAGNITUDE, MAGNITUDE, 0.0, -MAGNITUDE, 0.0, 0.0,
            ],
        );
    }

    #[test]
    fn cf3d_f19_get_grad_inactive_interval_preserves_existing_values() {
        // computeAngleTerm uses strict comparisons; an interior angle adds
        // zero to each already-populated XYZ component.
        let mut first = vec![0.0; 3];
        let mut middle = vec![0.0; 3];
        let mut third = vec![0.0; 3];
        let force_field = cf3d_f18_force_field(&mut first, &mut middle, &mut third);
        let contribution =
            AngleConstraintContrib::new_absolute(&force_field, 0, 1, 2, 80.0, 100.0, 2.0).unwrap();
        let trial = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let original = [1.0, -2.0, 3.0, -4.0, 5.0, -6.0, 7.0, -8.0, 9.0];
        let mut gradient = original;

        contribution.get_grad(&trial, &mut gradient);

        assert_eq!(gradient, original);
    }

    #[test]
    fn cf3d_f19_get_grad_collinear_and_zero_vectors_use_source_floors() {
        // Both the angle-vector squared lengths and rp.length() have source
        // floors; finite collinear and coincident inputs yield finite zeros.
        let mut first = vec![0.0; 3];
        let mut middle = vec![0.0; 3];
        let mut third = vec![0.0; 3];
        let force_field = cf3d_f18_force_field(&mut first, &mut middle, &mut third);
        let contribution =
            AngleConstraintContrib::new_absolute(&force_field, 0, 1, 2, 10.0, 20.0, 1.0).unwrap();

        for trial in [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0], [0.0; 9]] {
            let mut gradient = [0.0; 9];
            contribution.get_grad(&trial, &mut gradient);
            assert_eq!(gradient, [0.0; 9]);
            assert!(gradient.iter().all(|value| value.is_finite()));
        }
    }

    #[test]
    fn cf3d_f19_get_grad_adds_to_prepopulated_components() {
        // The source uses += for all three point gradients; its input buffer
        // is not cleared before this contribution is accumulated.
        const MAGNITUDE: f64 = 1145.9155902616465;
        let mut first = vec![0.0; 3];
        let mut middle = vec![0.0; 3];
        let mut third = vec![0.0; 3];
        let force_field = cf3d_f18_force_field(&mut first, &mut middle, &mut third);
        let contribution =
            AngleConstraintContrib::new_absolute(&force_field, 0, 1, 2, 100.0, 120.0, 1.0).unwrap();
        let trial = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let mut gradient = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];

        contribution.get_grad(&trial, &mut gradient);

        cf3d_f19_assert_gradient(
            &gradient,
            [
                1.0,
                2.0 + MAGNITUDE,
                3.0,
                4.0 - MAGNITUDE,
                5.0 - MAGNITUDE,
                6.0,
                7.0 + MAGNITUDE,
                8.0,
                9.0,
            ],
        );
    }

    #[test]
    fn cf3d_f19_get_grad_preserves_nan_comparison_and_arithmetic_behavior() {
        // std::max(1e-5, NaN) returns its first argument, std::clamp leaves
        // NaN unchanged, and the source continues into cross products.
        let mut first = vec![0.0; 3];
        let mut middle = vec![0.0; 3];
        let mut third = vec![0.0; 3];
        let force_field = cf3d_f18_force_field(&mut first, &mut middle, &mut third);
        let contribution =
            AngleConstraintContrib::new_absolute(&force_field, 0, 1, 2, 80.0, 100.0, 1.0).unwrap();
        let trial = [f64::NAN, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let mut gradient = [0.0; 9];

        contribution.get_grad(&trial, &mut gradient);

        assert!(gradient.iter().all(|value| value.is_nan()));
    }
}
