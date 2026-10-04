//! Private force-field geometry helpers.

// Copyright (C) 2004-2006 Rational Discovery LLC.
// @@ All Rights Reserved @@
// This file is part of the RDKit.
// The contents are covered by the terms of the BSD license
// which is included in the file license.txt, found at the root
// of the RDKit source tree.

// BEGIN RDKIT CPP FUNCTION RDKit::ForceFieldsHelper::normalizeAngleDeg (ForceField.cpp:19-26)
// RDKit✔️✔️: void normalizeAngleDeg(double &angleDeg) {
// RDKit✔️✔️:   angleDeg = fmod(angleDeg, 360.0);
// RDKit✔️✔️:   if (angleDeg < -180.0) {
// RDKit✔️✔️:     angleDeg += 360.0;
// RDKit✔️✔️:   } else if (angleDeg > 180.0) {
// RDKit✔️✔️:     angleDeg -= 360.0;
// RDKit✔️✔️:   }
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION RDKit::ForceFieldsHelper::normalizeAngleDeg
pub(super) fn normalize_angle_deg(angle_deg: &mut f64) {
    *angle_deg %= 360.0;
    if *angle_deg < -180.0 {
        *angle_deg += 360.0;
    } else if *angle_deg > 180.0 {
        *angle_deg -= 360.0;
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(super) struct Point3 {
    pub(super) x: f64,
    pub(super) y: f64,
    pub(super) z: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum DirectionVectorError {
    LengthBelowSourceTolerance,
}

const RDKIT_ZERO_TOLERANCE: f64 = 1.0e-16;

pub(super) fn direction_vector(
    origin: &Point3,
    other: &Point3,
) -> Result<Point3, DirectionVectorError> {
    // BEGIN RDKIT CPP HELPER RDGeom::Point3D::directionVector (point.h:214-220)
    // RDKit❗✔️: Point3D directionVector(const Point3D &other) const {
    // RDKit❗✔️:   Point3D res;
    // RDKit❗✔️:   res.x = other.x - x;
    // RDKit❗✔️:   res.y = other.y - y;
    // RDKit❗✔️:   res.z = other.z - z;
    // RDKit❗✔️:   res.normalize();
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP HELPER RDGeom::Point3D::directionVector
    // BEGIN RDKIT CPP HELPER RDGeom::Point3D::normalize (point.h:147-155)
    // RDKit❗✔️: constexpr void normalize() override {
    // RDKit❗✔️:   double l = this->length();
    // RDKit❗✔️:   if (l < zero_tolerance) {
    // RDKit❗✔️:     throw std::runtime_error("Cannot normalize a zero length vector");
    // RDKit❗✔️:   }
    // RDKit❗✔️:   x /= l;
    // RDKit❗✔️:   y /= l;
    // RDKit❗✔️:   z /= l;
    // RDKit❗✔️: }
    // END RDKIT CPP HELPER RDGeom::Point3D::normalize
    let mut direction = Point3::difference(other, origin);
    let length = direction.length();
    if length < RDKIT_ZERO_TOLERANCE {
        return Err(DirectionVectorError::LengthBelowSourceTolerance);
    }
    direction.divide_assign(length);
    // Behavior marker — RDKit❗✔️: B16 zero-length and NaN regressions follow.
    // Complexity marker — RDKit✔️✔️: fixed scalar subtraction, length and
    // division; the normalized point is a stack value with no allocation.
    Ok(direction)
}

impl Point3 {
    pub(super) fn difference(p1: &Self, p2: &Self) -> Self {
        // BEGIN RDKIT CPP HELPER RDGeom::operator- (Geometry/point.cpp:64-70)
        // RDKit✔️✔️: Point3D operator-(const Point3D &p1, const Point3D &p2) {
        // RDKit✔️✔️:   Point3D res;
        // RDKit✔️✔️:   res.x = p1.x - p2.x;
        // RDKit✔️✔️:   res.y = p1.y - p2.y;
        // RDKit✔️✔️:   res.z = p1.z - p2.z;
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        // END RDKIT CPP HELPER RDGeom::operator-
        Self {
            x: p1.x - p2.x,
            y: p1.y - p2.y,
            z: p1.z - p2.z,
        }
    }

    fn negated(&self) -> Self {
        // BEGIN RDKIT CPP HELPER RDGeom::Point3D::operator- (Geometry/point.h:139-145)
        // RDKit✔️✔️: constexpr Point3D operator-() const {
        // RDKit✔️✔️:   Point3D res(x, y, z);
        // RDKit✔️✔️:   res.x *= -1.0;
        // RDKit✔️✔️:   res.y *= -1.0;
        // RDKit✔️✔️:   res.z *= -1.0;
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        // END RDKIT CPP HELPER RDGeom::Point3D::operator-
        Self {
            x: self.x * -1.0,
            y: self.y * -1.0,
            z: self.z * -1.0,
        }
    }

    pub(super) fn cross_product(&self, other: &Self) -> Self {
        // BEGIN RDKIT CPP HELPER RDGeom::Point3D::crossProduct (Geometry/point.h:228-234)
        // RDKit✔️✔️: constexpr Point3D crossProduct(const Point3D &other) const {
        // RDKit✔️✔️:   Point3D res;
        // RDKit✔️✔️:   res.x = y * (other.z) - z * (other.y);
        // RDKit✔️✔️:   res.y = -x * (other.z) + z * (other.x);
        // RDKit✔️✔️:   res.z = x * (other.y) - y * (other.x);
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        // END RDKIT CPP HELPER RDGeom::Point3D::crossProduct
        Self {
            x: self.y * other.z - self.z * other.y,
            y: -self.x * other.z + self.z * other.x,
            z: self.x * other.y - self.y * other.x,
        }
    }

    pub(super) fn length(&self) -> f64 {
        // BEGIN RDKIT CPP HELPER RDGeom::Point3D::length (Geometry/point.h:158-161)
        // RDKit✔️✔️: double length() const override {
        // RDKit✔️✔️:   double res = x * x + y * y + z * z;
        // RDKit✔️✔️:   return sqrt(res);
        // RDKit✔️✔️: }
        // END RDKIT CPP HELPER RDGeom::Point3D::length
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    pub(super) fn length_sq(&self) -> f64 {
        // BEGIN RDKIT CPP HELPER RDGeom::Point3D::lengthSq (Geometry/point.h:163-168)
        // RDKit✔️✔️: constexpr double lengthSq() const override {
        // RDKit✔️✔️:   // double res = pow(x,2) + pow(y,2) + pow(z,2);
        // RDKit✔️✔️:   double res = x * x + y * y + z * z;
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        // END RDKIT CPP HELPER RDGeom::Point3D::lengthSq
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    pub(super) fn sum(p1: &Self, p2: &Self) -> Self {
        // BEGIN RDKIT CPP HELPER RDGeom::operator+ (Geometry/point.cpp:56-62)
        // RDKit✔️✔️: Point3D operator+(const Point3D &p1, const Point3D &p2) {
        // RDKit✔️✔️:   Point3D res;
        // RDKit✔️✔️:   res.x = p1.x + p2.x;
        // RDKit✔️✔️:   res.y = p1.y + p2.y;
        // RDKit✔️✔️:   res.z = p1.z + p2.z;
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        // END RDKIT CPP HELPER RDGeom::operator+
        Self {
            x: p1.x + p2.x,
            y: p1.y + p2.y,
            z: p1.z + p2.z,
        }
    }

    pub(super) fn scaled(&self, scale: f64) -> Self {
        // BEGIN RDKIT CPP HELPER RDGeom::operator* (Geometry/point.cpp:72-78)
        // RDKit✔️✔️: Point3D operator*(const Point3D &p1, double v) {
        // RDKit✔️✔️:   Point3D res;
        // RDKit✔️✔️:   res.x = p1.x * v;
        // RDKit✔️✔️:   res.y = p1.y * v;
        // RDKit✔️✔️:   res.z = p1.z * v;
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        // END RDKIT CPP HELPER RDGeom::operator*
        Self {
            x: self.x * scale,
            y: self.y * scale,
            z: self.z * scale,
        }
    }

    pub(super) fn divided(&self, scale: f64) -> Self {
        // BEGIN RDKIT CPP HELPER RDGeom::operator/ (Geometry/point.cpp:80-86)
        // RDKit✔️✔️: Point3D operator/(const Point3D &p1, double v) {
        // RDKit✔️✔️:   Point3D res;
        // RDKit✔️✔️:   res.x = p1.x / v;
        // RDKit✔️✔️:   res.y = p1.y / v;
        // RDKit✔️✔️:   res.z = p1.z / v;
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        // END RDKIT CPP HELPER RDGeom::operator/
        Self {
            x: self.x / scale,
            y: self.y / scale,
            z: self.z / scale,
        }
    }

    pub(super) fn dot_product(&self, other: &Self) -> f64 {
        // BEGIN RDKIT CPP HELPER RDGeom::Point3D::dotProduct (Geometry/point.h:169-172)
        // RDKit✔️✔️: constexpr double dotProduct(const Point3D &other) const {
        // RDKit✔️✔️:   double res = x * (other.x) + y * (other.y) + z * (other.z);
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        // END RDKIT CPP HELPER RDGeom::Point3D::dotProduct
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub(super) fn divide_assign(&mut self, scale: f64) {
        // BEGIN RDKIT CPP HELPER RDGeom::Point3D::operator/= (Geometry/point.h:132-137)
        // RDKit✔️✔️: constexpr Point3D &operator/=(double scale) {
        // RDKit✔️✔️:   x /= scale;
        // RDKit✔️✔️:   y /= scale;
        // RDKit✔️✔️:   z /= scale;
        // RDKit✔️✔️:   return *this;
        // RDKit✔️✔️: }
        // END RDKIT CPP HELPER RDGeom::Point3D::operator/=
        self.x /= scale;
        self.y /= scale;
        self.z /= scale;
    }
}

#[allow(clippy::too_many_arguments)]
fn compute_dihedral_from_points(
    p1: &Point3,
    p2: &Point3,
    p3: &Point3,
    p4: &Point3,
    dihedral: Option<&mut f64>,
    cos_phi: Option<&mut f64>,
    r: Option<&mut [Point3; 4]>,
    t: Option<&mut [Point3; 2]>,
    d: Option<&mut [f64; 2]>,
) {
    // BEGIN RDKIT CPP FUNCTION RDKit::ForceFieldsHelper::computeDihedral (ForceField.cpp:50-92)
    // RDKit✔️✔️: void computeDihedral(const RDGeom::Point3D *p1, const RDGeom::Point3D *p2,
    // RDKit✔️✔️:                      const RDGeom::Point3D *p3, const RDGeom::Point3D *p4,
    // RDKit✔️✔️:                      double *dihedral, double *cosPhi, RDGeom::Point3D r[4],
    // RDKit✔️✔️:                      RDGeom::Point3D t[2], double d[2]) {
    // RDKit✔️✔️:   PRECONDITION(p1, "p1 must not be null");
    // RDKit✔️✔️:   PRECONDITION(p2, "p2 must not be null");
    // RDKit✔️✔️:   PRECONDITION(p3, "p3 must not be null");
    // RDKit✔️✔️:   PRECONDITION(p4, "p4 must not be null");
    // RDKit✔️✔️:   RDGeom::Point3D rLocal[4];
    // RDKit✔️✔️:   RDGeom::Point3D tLocal[2];
    // RDKit✔️✔️:   double dLocal[2];
    // RDKit✔️✔️:   if (!r) {
    // RDKit✔️✔️:     r = rLocal;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!t) {
    // RDKit✔️✔️:     t = tLocal;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!d) {
    // RDKit✔️✔️:     d = dLocal;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   r[0] = *p1 - *p2;
    // RDKit✔️✔️:   r[1] = *p3 - *p2;
    // RDKit✔️✔️:   r[2] = -r[1];
    // RDKit✔️✔️:   r[3] = *p4 - *p3;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   t[0] = r[0].crossProduct(r[1]);
    // RDKit✔️✔️:   d[0] = (std::max)(t[0].length(), 1.0e-5);
    // RDKit✔️✔️:   t[0] /= d[0];
    // RDKit✔️✔️:   t[1] = r[2].crossProduct(r[3]);
    // RDKit✔️✔️:   d[1] = (std::max)(t[1].length(), 1.0e-5);
    // RDKit✔️✔️:   t[1] /= d[1];
    // RDKit✔️✔️:   double cosPhiLocal;
    // RDKit✔️✔️:   if (!cosPhi) {
    // RDKit✔️✔️:     cosPhi = &cosPhiLocal;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   *cosPhi = (std::max)(-1.0, (std::min)(t[0].dotProduct(t[1]), 1.0));
    // RDKit✔️✔️:   // we want a signed dihedral, that's why we use atan2 instead of acos
    // RDKit✔️✔️:   if (dihedral) {
    // RDKit✔️✔️:     RDGeom::Point3D m = t[0].crossProduct(r[1]);
    // RDKit✔️✔️:     double mLength = (std::max)(m.length(), 1.0e-5);
    // RDKit✔️✔️:     *dihedral = -atan2(m.dotProduct(t[1]) / mLength, *cosPhi);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RDKit::ForceFieldsHelper::computeDihedral

    let mut r_local = [Point3::default(); 4];
    let r = r.unwrap_or(&mut r_local);
    let mut t_local = [Point3::default(); 2];
    let t = t.unwrap_or(&mut t_local);
    let mut d_local = [0.0; 2];
    let d = d.unwrap_or(&mut d_local);

    r[0] = Point3::difference(p1, p2);
    r[1] = Point3::difference(p3, p2);
    r[2] = r[1].negated();
    r[3] = Point3::difference(p4, p3);

    t[0] = r[0].cross_product(&r[1]);
    let t0_length = t[0].length();
    d[0] = if t0_length < 1.0e-5 {
        1.0e-5
    } else {
        t0_length
    };
    t[0].divide_assign(d[0]);

    t[1] = r[2].cross_product(&r[3]);
    let t1_length = t[1].length();
    d[1] = if t1_length < 1.0e-5 {
        1.0e-5
    } else {
        t1_length
    };
    t[1].divide_assign(d[1]);

    let mut cos_phi_local = 0.0;
    let cos_phi = cos_phi.unwrap_or(&mut cos_phi_local);
    let normal_dot = t[0].dot_product(&t[1]);
    let clamped_upper = if 1.0 < normal_dot { 1.0 } else { normal_dot };
    *cos_phi = if -1.0 < clamped_upper {
        clamped_upper
    } else {
        -1.0
    };

    if let Some(dihedral) = dihedral {
        let m = t[0].cross_product(&r[1]);
        let m_length = m.length();
        let m_length = if m_length < 1.0e-5 { 1.0e-5 } else { m_length };
        *dihedral = -(m.dot_product(&t[1]) / m_length).atan2(*cos_phi);
    }
}

#[allow(clippy::too_many_arguments)]
fn compute_dihedral_from_position_vec(
    pos: &[Point3],
    idx1: usize,
    idx2: usize,
    idx3: usize,
    idx4: usize,
    dihedral: Option<&mut f64>,
    cos_phi: Option<&mut f64>,
    r: Option<&mut [Point3; 4]>,
    t: Option<&mut [Point3; 2]>,
    d: Option<&mut [f64; 2]>,
) {
    // BEGIN RDKIT CPP FUNCTION RDKit::ForceFieldsHelper::computeDihedral (PointPtrVect overload, ForceField.cpp:28-37)
    // RDKit✔️✔️: void computeDihedral(const RDGeom::PointPtrVect &pos, unsigned int idx1,
    // RDKit✔️✔️:                    unsigned int idx2, unsigned int idx3, unsigned int idx4,
    // RDKit✔️✔️:                    double *dihedral, double *cosPhi, RDGeom::Point3D r[4],
    // RDKit✔️✔️:                    RDGeom::Point3D t[2], double d[2]) {
    // RDKit✔️✔️:   computeDihedral(static_cast<RDGeom::Point3D *>(pos[idx1]),
    // RDKit✔️✔️:                   static_cast<RDGeom::Point3D *>(pos[idx2]),
    // RDKit✔️✔️:                   static_cast<RDGeom::Point3D *>(pos[idx3]),
    // RDKit✔️✔️:                   static_cast<RDGeom::Point3D *>(pos[idx4]), dihedral, cosPhi,
    // RDKit✔️✔️:                   r, t, d);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RDKit::ForceFieldsHelper::computeDihedral (PointPtrVect overload)
    compute_dihedral_from_points(
        &pos[idx1], &pos[idx2], &pos[idx3], &pos[idx4], dihedral, cos_phi, r, t, d,
    );
}

#[allow(clippy::too_many_arguments)]
pub(super) fn compute_dihedral_from_flat(
    pos: &[f64],
    idx1: usize,
    idx2: usize,
    idx3: usize,
    idx4: usize,
    dihedral: Option<&mut f64>,
    cos_phi: Option<&mut f64>,
    r: Option<&mut [Point3; 4]>,
    t: Option<&mut [Point3; 2]>,
    d: Option<&mut [f64; 2]>,
) {
    // BEGIN RDKIT CPP FUNCTION RDKit::ForceFieldsHelper::computeDihedral (flat overload, ForceField.cpp:39-48)
    // RDKit✔️✔️: void computeDihedral(const double *pos, unsigned int idx1, unsigned int idx2,
    // RDKit✔️✔️:                    unsigned int idx3, unsigned int idx4, double *dihedral,
    // RDKit✔️✔️:                    double *cosPhi, RDGeom::Point3D r[4], RDGeom::Point3D t[2],
    // RDKit✔️✔️:                    double d[2]) {
    // RDKit✔️✔️:   RDGeom::Point3D p1(pos[3 * idx1], pos[3 * idx1 + 1], pos[3 * idx1 + 2]);
    // RDKit✔️✔️:   RDGeom::Point3D p2(pos[3 * idx2], pos[3 * idx2 + 1], pos[3 * idx2 + 2]);
    // RDKit✔️✔️:   RDGeom::Point3D p3(pos[3 * idx3], pos[3 * idx3 + 1], pos[3 * idx3 + 2]);
    // RDKit✔️✔️:   RDGeom::Point3D p4(pos[3 * idx4], pos[3 * idx4 + 1], pos[3 * idx4 + 2]);
    // RDKit✔️✔️:   computeDihedral(&p1, &p2, &p3, &p4, dihedral, cosPhi, r, t, d);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RDKit::ForceFieldsHelper::computeDihedral (flat overload)
    let p1 = Point3 {
        x: pos[3 * idx1],
        y: pos[3 * idx1 + 1],
        z: pos[3 * idx1 + 2],
    };
    let p2 = Point3 {
        x: pos[3 * idx2],
        y: pos[3 * idx2 + 1],
        z: pos[3 * idx2 + 2],
    };
    let p3 = Point3 {
        x: pos[3 * idx3],
        y: pos[3 * idx3 + 1],
        z: pos[3 * idx3 + 2],
    };
    let p4 = Point3 {
        x: pos[3 * idx4],
        y: pos[3 * idx4 + 1],
        z: pos[3 * idx4 + 2],
    };
    compute_dihedral_from_points(&p1, &p2, &p3, &p4, dihedral, cos_phi, r, t, d);
}

pub(super) fn compute_dihedral_from_flat_radians(
    pos: &[f64],
    idx1: usize,
    idx2: usize,
    idx3: usize,
    idx4: usize,
) -> f64 {
    // The pinned flat overload and its Point3D calculation remain anchored in
    // `compute_dihedral_from_flat`; this wrapper keeps the private Point3 type
    // out of the sibling kernel interface.
    let mut dihedral = 0.0;
    compute_dihedral_from_flat(
        pos,
        idx1,
        idx2,
        idx3,
        idx4,
        Some(&mut dihedral),
        None,
        None,
        None,
        None,
    );
    dihedral
}

pub(super) fn compute_dihedral_from_flat_with_vectors(
    pos: &[f64],
    idx1: usize,
    idx2: usize,
    idx3: usize,
    idx4: usize,
) -> (f64, [Point3; 4]) {
    let mut dihedral = 0.0;
    let mut r = [Point3::default(); 4];
    compute_dihedral_from_flat(
        pos,
        idx1,
        idx2,
        idx3,
        idx4,
        Some(&mut dihedral),
        None,
        Some(&mut r),
        None,
        None,
    );
    (dihedral, r)
}

pub(super) fn compute_dihedral_from_position_slices_radians(
    p1: &[f64],
    p2: &[f64],
    p3: &[f64],
    p4: &[f64],
) -> f64 {
    // BEGIN RDKIT CPP FUNCTION RDKit::ForceFieldsHelper::computeDihedral (PointPtrVect overload, ForceField.cpp:28-37)
    // RDKit✔️✔️: void computeDihedral(const RDGeom::PointPtrVect &pos, unsigned int idx1,
    // RDKit✔️✔️:                    unsigned int idx2, unsigned int idx3, unsigned int idx4,
    // RDKit✔️✔️:                    double *dihedral, double *cosPhi, RDGeom::Point3D r[4],
    // RDKit✔️✔️:                    RDGeom::Point3D t[2], double d[2]) {
    // RDKit✔️✔️:   computeDihedral(static_cast<RDGeom::Point3D *>(pos[idx1]),
    // RDKit✔️✔️:                   static_cast<RDGeom::Point3D *>(pos[idx2]),
    // RDKit✔️✔️:                   static_cast<RDGeom::Point3D *>(pos[idx3]),
    // RDKit✔️✔️:                   static_cast<RDGeom::Point3D *>(pos[idx4]), dihedral, cosPhi,
    // RDKit✔️✔️:                   r, t, d);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RDKit::ForceFieldsHelper::computeDihedral (PointPtrVect overload)
    // ForceField stores each point as a borrowed scalar row, so materialize
    // only these four source Point3D values on the stack. The full shared
    // computation remains in `compute_dihedral_from_points`.
    let p1 = Point3 {
        x: p1[0],
        y: p1[1],
        z: p1[2],
    };
    let p2 = Point3 {
        x: p2[0],
        y: p2[1],
        z: p2[2],
    };
    let p3 = Point3 {
        x: p3[0],
        y: p3[1],
        z: p3[2],
    };
    let p4 = Point3 {
        x: p4[0],
        y: p4[1],
        z: p4[2],
    };
    let mut dihedral = 0.0;
    compute_dihedral_from_points(
        &p1,
        &p2,
        &p3,
        &p4,
        Some(&mut dihedral),
        None,
        None,
        None,
        None,
    );
    dihedral
}

#[cfg(test)]
mod tests {
    use super::{
        Point3, compute_dihedral_from_flat, compute_dihedral_from_points,
        compute_dihedral_from_position_vec, normalize_angle_deg,
    };

    fn point(x: f64, y: f64, z: f64) -> Point3 {
        Point3 { x, y, z }
    }

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 1.0e-12,
            "expected {expected}, got {actual}"
        );
    }

    fn calculate_all(
        p1: &Point3,
        p2: &Point3,
        p3: &Point3,
        p4: &Point3,
    ) -> (f64, f64, [Point3; 4], [Point3; 2], [f64; 2]) {
        let mut dihedral = 0.0;
        let mut cos_phi = 0.0;
        let mut r = [Point3::default(); 4];
        let mut t = [Point3::default(); 2];
        let mut d = [0.0; 2];
        compute_dihedral_from_points(
            p1,
            p2,
            p3,
            p4,
            Some(&mut dihedral),
            Some(&mut cos_phi),
            Some(&mut r),
            Some(&mut t),
            Some(&mut d),
        );
        (dihedral, cos_phi, r, t, d)
    }

    fn assert_f03_reference(expected: &(f64, f64, [Point3; 4], [Point3; 2], [f64; 2])) {
        assert_eq!(
            expected.2,
            [
                point(1.0, 0.0, 0.0),
                point(0.0, 1.0, 0.0),
                point(0.0, -1.0, 0.0),
                point(0.0, 0.0, 1.0),
            ]
        );
        assert_eq!(expected.3, [point(0.0, 0.0, 1.0), point(-1.0, 0.0, 0.0)]);
        assert_eq!(expected.4, [1.0, 1.0]);
        assert_eq!(expected.1, 0.0);
        assert_close(expected.0, -std::f64::consts::FRAC_PI_2);
    }

    fn assert_f03_outputs(
        output_mask: usize,
        expected: &(f64, f64, [Point3; 4], [Point3; 2], [f64; 2]),
        dihedral: f64,
        cos_phi: f64,
        r: [Point3; 4],
        t: [Point3; 2],
        d: [f64; 2],
    ) {
        let want_dihedral = output_mask & 0b00001 != 0;
        let want_cos_phi = output_mask & 0b00010 != 0;
        let want_r = output_mask & 0b00100 != 0;
        let want_t = output_mask & 0b01000 != 0;
        let want_d = output_mask & 0b10000 != 0;

        if want_dihedral {
            assert_close(dihedral, expected.0);
        } else {
            assert_eq!(dihedral, 17.0, "angle output mask {output_mask}");
        }
        if want_cos_phi {
            assert_close(cos_phi, expected.1);
        } else {
            assert_eq!(cos_phi, 17.0, "cosine output mask {output_mask}");
        }
        assert_eq!(
            r,
            if want_r {
                expected.2
            } else {
                [point(17.0, 17.0, 17.0); 4]
            },
            "r output mask {output_mask}"
        );
        assert_eq!(
            t,
            if want_t {
                expected.3
            } else {
                [point(17.0, 17.0, 17.0); 2]
            },
            "t output mask {output_mask}"
        );
        assert_eq!(
            d,
            if want_d { expected.4 } else { [17.0; 2] },
            "d output mask {output_mask}"
        );
    }

    #[test]
    fn cf3d_f03_position_vector_adapter_uses_noncontiguous_indices_and_all_outputs() {
        let points = [
            point(8.0, 8.0, 8.0),
            point(1.0, 0.0, 0.0),
            point(9.0, 0.0, 0.0),
            point(0.0, 0.0, 0.0),
            point(0.0, 9.0, 0.0),
            point(0.0, 1.0, 0.0),
            point(0.0, 0.0, 9.0),
            point(0.0, 1.0, 1.0),
        ];
        let expected = calculate_all(&points[1], &points[3], &points[5], &points[7]);
        assert_f03_reference(&expected);

        // RDKit ForceField.cpp:28-37 forwards the four indexed points in order.
        for output_mask in 0..32 {
            let mut dihedral = 17.0;
            let mut cos_phi = 17.0;
            let mut r = [point(17.0, 17.0, 17.0); 4];
            let mut t = [point(17.0, 17.0, 17.0); 2];
            let mut d = [17.0; 2];

            compute_dihedral_from_position_vec(
                &points,
                1,
                3,
                5,
                7,
                if output_mask & 0b00001 != 0 {
                    Some(&mut dihedral)
                } else {
                    None
                },
                if output_mask & 0b00010 != 0 {
                    Some(&mut cos_phi)
                } else {
                    None
                },
                if output_mask & 0b00100 != 0 {
                    Some(&mut r)
                } else {
                    None
                },
                if output_mask & 0b01000 != 0 {
                    Some(&mut t)
                } else {
                    None
                },
                if output_mask & 0b10000 != 0 {
                    Some(&mut d)
                } else {
                    None
                },
            );

            assert_f03_outputs(output_mask, &expected, dihedral, cos_phi, r, t, d);
        }
    }

    #[test]
    fn cf3d_f03_flat_adapter_uses_noncontiguous_indices_and_all_outputs() {
        let points = [
            point(8.0, 8.0, 8.0),
            point(1.0, 0.0, 0.0),
            point(9.0, 0.0, 0.0),
            point(0.0, 0.0, 0.0),
            point(0.0, 9.0, 0.0),
            point(0.0, 1.0, 0.0),
            point(0.0, 0.0, 9.0),
            point(0.0, 1.0, 1.0),
        ];
        let flat = points
            .iter()
            .flat_map(|p| [p.x, p.y, p.z])
            .collect::<Vec<_>>();
        let expected = calculate_all(&points[1], &points[3], &points[5], &points[7]);
        assert_f03_reference(&expected);

        // RDKit ForceField.cpp:39-48 reads the selected triples with stride three.
        for output_mask in 0..32 {
            let mut dihedral = 17.0;
            let mut cos_phi = 17.0;
            let mut r = [point(17.0, 17.0, 17.0); 4];
            let mut t = [point(17.0, 17.0, 17.0); 2];
            let mut d = [17.0; 2];

            compute_dihedral_from_flat(
                &flat,
                1,
                3,
                5,
                7,
                if output_mask & 0b00001 != 0 {
                    Some(&mut dihedral)
                } else {
                    None
                },
                if output_mask & 0b00010 != 0 {
                    Some(&mut cos_phi)
                } else {
                    None
                },
                if output_mask & 0b00100 != 0 {
                    Some(&mut r)
                } else {
                    None
                },
                if output_mask & 0b01000 != 0 {
                    Some(&mut t)
                } else {
                    None
                },
                if output_mask & 0b10000 != 0 {
                    Some(&mut d)
                } else {
                    None
                },
            );

            assert_f03_outputs(output_mask, &expected, dihedral, cos_phi, r, t, d);
        }
    }

    #[test]
    fn cf3d_f02_positive_negative_torsions_and_planar_values() {
        let p1 = point(1.0, 0.0, 0.0);
        let p2 = point(0.0, 0.0, 0.0);
        let p3 = point(0.0, 1.0, 0.0);

        let (cis_dihedral, cis_cos, ..) = calculate_all(&p1, &p2, &p3, &point(1.0, 1.0, 0.0));
        assert_close(cis_dihedral, 0.0);
        assert_close(cis_cos, 1.0);

        let (trans_dihedral, trans_cos, ..) = calculate_all(&p1, &p2, &p3, &point(-1.0, 1.0, 0.0));
        assert_close(trans_dihedral.abs(), std::f64::consts::PI);
        assert_close(trans_cos, -1.0);

        let (positive_z_dihedral, positive_z_cos, ..) =
            calculate_all(&p1, &p2, &p3, &point(0.5, 1.0, 0.866_025_403_784_438_6));
        assert_close(positive_z_dihedral, -std::f64::consts::FRAC_PI_3);
        assert_close(positive_z_cos, 0.5);

        let (negative_z_dihedral, negative_z_cos, ..) =
            calculate_all(&p1, &p2, &p3, &point(0.5, 1.0, -0.866_025_403_784_438_6));
        assert_close(negative_z_dihedral, std::f64::consts::FRAC_PI_3);
        assert_close(negative_z_cos, 0.5);
    }

    #[test]
    fn cf3d_f02_collinear_and_zero_vector_floors_match_source() {
        let origin = point(0.0, 0.0, 0.0);
        let (dihedral, cos_phi, r, t, d) = calculate_all(&origin, &origin, &origin, &origin);
        assert_eq!(r, [Point3::default(); 4]);
        assert_eq!(t, [Point3::default(); 2]);
        assert_eq!(d, [1.0e-5; 2]);
        assert_eq!(cos_phi, 0.0);
        assert_eq!(dihedral, 0.0);
        assert!(dihedral.is_sign_negative());

        let (dihedral, cos_phi, _, t, d) = calculate_all(
            &point(1.0, 0.0, 0.0),
            &origin,
            &point(-1.0, 0.0, 0.0),
            &point(-2.0, 0.0, 0.0),
        );
        assert_eq!(t, [Point3::default(); 2]);
        assert_eq!(d, [1.0e-5; 2]);
        assert_eq!(cos_phi, 0.0);
        assert_eq!(dihedral, 0.0);
        assert!(dihedral.is_sign_negative());
    }

    #[test]
    fn cf3d_f02_all_nullable_output_combinations_preserve_source_results() {
        let p1 = point(1.0, 0.0, 0.0);
        let p2 = point(0.0, 0.0, 0.0);
        let p3 = point(0.0, 1.0, 0.0);
        let p4 = point(0.0, 1.0, 1.0);
        let expected = calculate_all(&p1, &p2, &p3, &p4);
        assert_eq!(
            expected.2,
            [
                point(1.0, 0.0, 0.0),
                point(0.0, 1.0, 0.0),
                point(0.0, -1.0, 0.0),
                point(0.0, 0.0, 1.0),
            ]
        );
        assert_eq!(expected.3, [point(0.0, 0.0, 1.0), point(-1.0, 0.0, 0.0)]);
        assert_eq!(expected.4, [1.0, 1.0]);
        assert_eq!(expected.1, 0.0);
        assert_close(expected.0, -std::f64::consts::FRAC_PI_2);

        for output_mask in 0..32 {
            let mut dihedral = 17.0;
            let mut cos_phi = 17.0;
            let mut r = [point(17.0, 17.0, 17.0); 4];
            let mut t = [point(17.0, 17.0, 17.0); 2];
            let mut d = [17.0; 2];
            let want_dihedral = output_mask & 0b00001 != 0;
            let want_cos_phi = output_mask & 0b00010 != 0;
            let want_r = output_mask & 0b00100 != 0;
            let want_t = output_mask & 0b01000 != 0;
            let want_d = output_mask & 0b10000 != 0;

            compute_dihedral_from_points(
                &p1,
                &p2,
                &p3,
                &p4,
                if want_dihedral {
                    Some(&mut dihedral)
                } else {
                    None
                },
                if want_cos_phi {
                    Some(&mut cos_phi)
                } else {
                    None
                },
                if want_r { Some(&mut r) } else { None },
                if want_t { Some(&mut t) } else { None },
                if want_d { Some(&mut d) } else { None },
            );

            if want_dihedral {
                assert!(
                    (dihedral - expected.0).abs() < 1.0e-12,
                    "angle output mask {output_mask}: expected {}, got {dihedral}",
                    expected.0
                );
            } else {
                assert_eq!(dihedral, 17.0, "angle output mask {output_mask}");
            }
            if want_cos_phi {
                assert!(
                    (cos_phi - expected.1).abs() < 1.0e-12,
                    "cosine output mask {output_mask}: expected {}, got {cos_phi}",
                    expected.1
                );
            } else {
                assert_eq!(cos_phi, 17.0, "cosine output mask {output_mask}");
            }
            assert_eq!(
                r,
                if want_r {
                    expected.2
                } else {
                    [point(17.0, 17.0, 17.0); 4]
                },
                "r output mask {output_mask}"
            );
            assert_eq!(
                t,
                if want_t {
                    expected.3
                } else {
                    [point(17.0, 17.0, 17.0); 2]
                },
                "t output mask {output_mask}"
            );
            assert_eq!(
                d,
                if want_d { expected.4 } else { [17.0; 2] },
                "d output mask {output_mask}"
            );
        }
    }

    #[test]
    fn cf3d_f02_nested_source_min_max_preserves_nan_order() {
        let nan_point = point(f64::NAN, 0.0, 0.0);
        let p2 = point(0.0, 0.0, 0.0);
        let p3 = point(0.0, 1.0, 0.0);
        let p4 = point(1.0, 1.0, 0.0);
        let (dihedral, cos_phi, _, _, d) = calculate_all(&nan_point, &p2, &p3, &p4);
        assert_eq!(cos_phi, -1.0);
        assert!(dihedral.is_nan());
        assert!(d[0].is_nan());
    }

    #[test]
    fn cf3d_f01_wraps_positive_and_negative_values_by_one_remainder() {
        for (input, expected) in [
            (-181.0, 179.0),
            (181.0, -179.0),
            (-540.0, -180.0),
            (540.0, 180.0),
            (-725.0, -5.0),
            (725.0, 5.0),
        ] {
            let mut angle = input;
            normalize_angle_deg(&mut angle);
            assert_eq!(angle, expected, "input {input}");
        }
    }

    #[test]
    fn cf3d_f01_preserves_strict_endpoints_and_signed_zero() {
        let mut lower_endpoint = -180.0;
        normalize_angle_deg(&mut lower_endpoint);
        assert_eq!(lower_endpoint, -180.0);

        let mut upper_endpoint = 180.0;
        normalize_angle_deg(&mut upper_endpoint);
        assert_eq!(upper_endpoint, 180.0);

        for input in [-360.0, -720.0, -1080.0, -0.0] {
            let mut angle = input;
            normalize_angle_deg(&mut angle);
            assert_eq!(angle, 0.0);
            assert!(angle.is_sign_negative(), "input {input}");
        }

        for input in [360.0, 720.0, 1080.0, 0.0] {
            let mut angle = input;
            normalize_angle_deg(&mut angle);
            assert_eq!(angle, 0.0);
            assert!(!angle.is_sign_negative(), "input {input}");
        }
    }

    #[test]
    fn cf3d_f01_nonfinite_inputs_follow_fmod_nan_propagation() {
        for input in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut angle = input;
            normalize_angle_deg(&mut angle);
            assert!(angle.is_nan(), "input {input}");
        }
    }
}
