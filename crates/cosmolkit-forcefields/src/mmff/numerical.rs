use super::params::{MmffAngle, MmffBond, MmffOop, MmffStbn, MmffTor};
use crate::geometry::Point3;

pub(super) fn calc_bond_rest_length(params: &MmffBond) -> f64 {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/ForceField/MMFF/BondStretch.cpp:20-24:
    // RDKit❗✔️: double calcBondRestLength(const MMFFBond *mmffBondParams) {
    // RDKit❗✔️:   PRECONDITION(mmffBondParams, "bond parameters not found");
    // RDKit❗✔️:
    // RDKit❗✔️:   return mmffBondParams->r0;
    // RDKit❗✔️: }
    // Behavior: the Rust reference requires a non-null existing value and
    // returns the source r0 field without conversion or fallback.
    // Complexity: one borrowed field read and scalar return, O(1), no allocation.
    params.r0
}

pub(super) fn calc_bond_force_constant(params: &MmffBond) -> f64 {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/ForceField/MMFF/BondStretch.cpp:26-30:
    // RDKit❗✔️: double calcBondForceConstant(const MMFFBond *mmffBondParams) {
    // RDKit❗✔️:   PRECONDITION(mmffBondParams, "bond parameters not found");
    // RDKit❗✔️:
    // RDKit❗✔️:   return mmffBondParams->kb;
    // RDKit❗✔️: }
    // Behavior: the Rust reference requires a non-null existing value and
    // returns the source kb field without conversion or fallback.
    // Complexity: one borrowed field read and scalar return, O(1), no allocation.
    params.kb
}

pub(super) fn calc_angle_rest_value(params: &MmffAngle) -> f64 {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/ForceField/MMFF/AngleBend.cpp:21-25:
    // RDKit❗✔️: double calcAngleRestValue(const MMFFAngle *mmffAngleParams) {
    // RDKit❗✔️:   PRECONDITION(mmffAngleParams, "angle parameters not found");
    // RDKit❗✔️:
    // RDKit❗✔️:   return mmffAngleParams->theta0;
    // RDKit❗✔️: }
    // Behavior: return the source theta0 field from the existing borrowed row.
    // Complexity: one borrowed field read and scalar return, O(1), no allocation.
    params.theta0
}

pub(super) fn calc_angle_force_constant(params: &MmffAngle) -> f64 {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/ForceField/MMFF/AngleBend.cpp:37-41:
    // RDKit❗✔️: double calcAngleForceConstant(const MMFFAngle *mmffAngleParams) {
    // RDKit❗✔️:   PRECONDITION(mmffAngleParams, "angle parameters not found");
    // RDKit❗✔️:
    // RDKit❗✔️:   return mmffAngleParams->ka;
    // RDKit❗✔️: }
    // Behavior: return the source ka field from the existing borrowed row.
    // Complexity: one borrowed field read and scalar return, O(1), no allocation.
    params.ka
}

pub(super) fn calc_stbn_force_constants(params: &MmffStbn) -> (f64, f64) {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/ForceField/MMFF/StretchBend.cpp:22-27:
    // RDKit❗✔️: std::pair<double, double> calcStbnForceConstants(
    // RDKit❗✔️:     const MMFFStbn *mmffStbnParams) {
    // RDKit❗✔️:   PRECONDITION(mmffStbnParams, "stretch-bend parameters not found");
    // RDKit❗✔️:
    // RDKit❗✔️:   return std::make_pair(mmffStbnParams->kbaIJK, mmffStbnParams->kbaKJI);
    // RDKit❗✔️: }
    // Behavior: retain the ordered source pair in a stack tuple from the
    // existing borrowed row. Complexity: two field reads, O(1), no allocation.
    (params.kba_ijk, params.kba_kji)
}

pub(super) fn calc_torsion_force_constant(params: &MmffTor) -> (f64, f64, f64) {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/ForceField/MMFF/TorsionAngle.cpp:40-44:
    // RDKit❗✔️: std::tuple<double, double, double> calcTorsionForceConstant(
    // RDKit❗✔️:     const MMFFTor *mmffTorParams) {
    // RDKit❗✔️:   return std::make_tuple(mmffTorParams->V1, mmffTorParams->V2,
    // RDKit❗✔️:                          mmffTorParams->V3);
    // RDKit❗✔️: }
    // Behavior: retain the source V1/V2/V3 order. Complexity: three borrowed
    // field reads into a stack tuple, O(1), no allocation.
    (params.v1, params.v2, params.v3)
}

pub(super) fn calc_oop_bend_force_constant(params: &MmffOop) -> f64 {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/ForceField/MMFF/OopBend.cpp:36-40:
    // RDKit❗✔️: double calcOopBendForceConstant(const MMFFOop *mmffOopParams) {
    // RDKit❗✔️:   PRECONDITION(mmffOopParams, "no OOP parameters");
    // RDKit❗✔️:
    // RDKit❗✔️:   return mmffOopParams->koop;
    // RDKit❗✔️: }
    // Behavior: return the source koop field from the existing borrowed row.
    // Complexity: one borrowed field read and scalar return, O(1), no allocation.
    params.koop
}

pub(super) fn calc_stretch_bend_energy(
    delta_dist1: f64,
    delta_dist2: f64,
    delta_theta: f64,
    force_constants: (f64, f64),
) -> (f64, f64) {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/ForceField/MMFF/StretchBend.cpp:29-38:
    // RDKit❗✔️: std::pair<double, double> calcStretchBendEnergy(
    // RDKit❗✔️:     const double deltaDist1, const double deltaDist2, const double deltaTheta,
    // RDKit❗✔️:     const std::pair<double, double> forceConstants) {
    // RDKit❗✔️:   double factor = MDYNE_A_TO_KCAL_MOL * DEG2RAD * deltaTheta;
    // RDKit❗✔️:
    // RDKit❗✔️:   return std::make_pair(factor * forceConstants.first * deltaDist1,
    // RDKit❗✔️:                         factor * forceConstants.second * deltaDist2);
    // RDKit❗✔️: }
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/ForceField/MMFF/Params.h:38-40:
    // RDKit❗✔️: constexpr double DEG2RAD = M_PI / 180.0;
    // RDKit❗✔️: constexpr double MDYNE_A_TO_KCAL_MOL = 143.9325;
    // Behavior: retain the source factor and product order for both terms.
    // Complexity: fixed scalar arithmetic, O(1), no allocation.
    let factor = 143.9325 * crate::uff::params::DEG2RAD * delta_theta;

    (
        factor * force_constants.0 * delta_dist1,
        factor * force_constants.1 * delta_dist2,
    )
}

pub(super) fn calc_torsion_energy(v1: f64, v2: f64, v3: f64, cos_phi: f64) -> f64 {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/ForceField/MMFF/TorsionAngle.cpp:46-54:
    // RDKit❗✔️: double calcTorsionEnergy(const double V1, const double V2, const double V3,
    // RDKit❗✔️:                          const double cosPhi) {
    // RDKit❗✔️:   double cos2Phi = 2.0 * cosPhi * cosPhi - 1.0;
    // RDKit❗✔️:   double cos3Phi = cosPhi * (2.0 * cos2Phi - 1.0);
    // RDKit❗✔️:
    // RDKit❗✔️:   return (0.5 *
    // RDKit❗✔️:           (V1 * (1.0 + cosPhi) + V2 * (1.0 - cos2Phi) + V3 * (1.0 + cos3Phi)));
    // RDKit❗✔️: }
    // Behavior: preserve the source recurrence and ordered energy sum, without
    // clipping the supplied cosine. Complexity: fixed scalar arithmetic, O(1).
    let cos2_phi = 2.0 * cos_phi * cos_phi - 1.0;
    let cos3_phi = cos_phi * (2.0 * cos2_phi - 1.0);

    0.5 * (v1 * (1.0 + cos_phi) + v2 * (1.0 - cos2_phi) + v3 * (1.0 + cos3_phi))
}

pub(super) fn calc_oop_bend_energy(chi: f64, koop: f64) -> f64 {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/ForceField/MMFF/OopBend.cpp:42-45:
    // RDKit❗✔️: double calcOopBendEnergy(const double chi, const double koop) {
    // RDKit❗✔️:   double const c2 = MDYNE_A_TO_KCAL_MOL * DEG2RAD * DEG2RAD;
    // RDKit❗✔️:   return (0.5 * c2 * koop * chi * chi);
    // RDKit❗✔️: }
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/ForceField/MMFF/Params.h:38-40:
    // RDKit❗✔️: constexpr double DEG2RAD = M_PI / 180.0;
    // RDKit❗✔️: constexpr double MDYNE_A_TO_KCAL_MOL = 143.9325;
    // Behavior: preserve the source factor and multiplication order; `chi`
    // is the supplied degree value. Complexity: fixed scalar arithmetic, O(1).
    let c2 = 143.9325 * crate::uff::params::DEG2RAD * crate::uff::params::DEG2RAD;
    0.5 * c2 * koop * chi * chi
}

pub(super) fn calc_cos_theta(p1: Point3, p2: Point3, p3: Point3, dist1: f64, dist2: f64) -> f64 {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/ForceField/MMFF/AngleBend.cpp:27-35:
    // RDKit❗✔️: double calcCosTheta(RDGeom::Point3D p1, RDGeom::Point3D p2, RDGeom::Point3D p3,
    // RDKit❗✔️:                     double dist1, double dist2) {
    // RDKit❗✔️:   RDGeom::Point3D p12 = p1 - p2;
    // RDKit❗✔️:   RDGeom::Point3D p32 = p3 - p2;
    // RDKit❗✔️:   double cosTheta = p12.dotProduct(p32) / (dist1 * dist2);
    // RDKit❗✔️:   clipToOne(cosTheta);
    // RDKit❗✔️:
    // RDKit❗✔️:   return cosTheta;
    // RDKit❗✔️: }
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/Geometry/point.cpp:64-70:
    // RDKit❗✔️: Point3D operator-(const Point3D &p1, const Point3D &p2) {
    // RDKit❗✔️:   Point3D res;
    // RDKit❗✔️:   res.x = p1.x - p2.x;
    // RDKit❗✔️:   res.y = p1.y - p2.y;
    // RDKit❗✔️:   res.z = p1.z - p2.z;
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/Geometry/point.h:169-172:
    // RDKit❗✔️: constexpr double dotProduct(const Point3D &other) const {
    // RDKit❗✔️:   double res = x * (other.x) + y * (other.y) + z * (other.z);
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/ForceField/MMFF/Params.h:44:
    // RDKit❗✔️: inline void clipToOne(double &x) { x = std::clamp(x, -1.0, 1.0); }
    // Behavior: use the existing source-shaped point operators, preserve the
    // supplied distances and clip only after division. Complexity: O(1), with
    // stack point values and no heap allocation.
    let p12 = Point3::difference(&p1, &p2);
    let p32 = Point3::difference(&p3, &p2);
    let mut cos_theta = p12.dot_product(&p32) / (dist1 * dist2);
    crate::uff::params::clip_to_one(&mut cos_theta);
    cos_theta
}

pub(crate) fn calc_torsion_cos_phi(i: &Point3, j: &Point3, k: &Point3, l: &Point3) -> f64 {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/ForceField/MMFF/TorsionAngle.cpp:19-38:
    // RDKit❗✔️: double calcTorsionCosPhi(const RDGeom::Point3D &iPoint,
    // RDKit❗✔️:                          const RDGeom::Point3D &jPoint,
    // RDKit❗✔️:                          const RDGeom::Point3D &kPoint,
    // RDKit❗✔️:                          const RDGeom::Point3D &lPoint) {
    // RDKit❗✔️:   RDGeom::Point3D r1 = iPoint - jPoint;
    // RDKit❗✔️:   RDGeom::Point3D r2 = kPoint - jPoint;
    // RDKit❗✔️:   RDGeom::Point3D r3 = jPoint - kPoint;
    // RDKit❗✔️:   RDGeom::Point3D r4 = lPoint - kPoint;
    // RDKit❗✔️:   RDGeom::Point3D t1 = r1.crossProduct(r2);
    // RDKit❗✔️:   RDGeom::Point3D t2 = r3.crossProduct(r4);
    // RDKit❗✔️:   auto t1_len = t1.length();
    // RDKit❗✔️:   auto t2_len = t2.length();
    // RDKit❗✔️:   if (isDoubleZero(t1_len) || isDoubleZero(t2_len)) {
    // RDKit❗✔️:     return 0.0;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   double cosPhi = t1.dotProduct(t2) / (t1_len * t2_len);
    // RDKit❗✔️:   clipToOne(cosPhi);
    // RDKit❗✔️:
    // RDKit❗✔️:   return cosPhi;
    // RDKit❗✔️: }
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/Geometry/point.cpp:64-70:
    // RDKit❗✔️: Point3D operator-(const Point3D &p1, const Point3D &p2) {
    // RDKit❗✔️:   Point3D res;
    // RDKit❗✔️:   res.x = p1.x - p2.x;
    // RDKit❗✔️:   res.y = p1.y - p2.y;
    // RDKit❗✔️:   res.z = p1.z - p2.z;
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/Geometry/point.h:228-234:
    // RDKit❗✔️: constexpr Point3D crossProduct(const Point3D &other) const {
    // RDKit❗✔️:   Point3D res;
    // RDKit❗✔️:   res.x = y * (other.z) - z * (other.y);
    // RDKit❗✔️:   res.y = -x * (other.z) + z * (other.x);
    // RDKit❗✔️:   res.z = x * (other.y) - y * (other.x);
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/Geometry/point.h:158-161:
    // RDKit❗✔️: double length() const override {
    // RDKit❗✔️:   double res = x * x + y * y + z * z;
    // RDKit❗✔️:   return sqrt(res);
    // RDKit❗✔️: }
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/Geometry/point.h:169-172:
    // RDKit❗✔️: constexpr double dotProduct(const Point3D &other) const {
    // RDKit❗✔️:   double res = x * (other.x) + y * (other.y) + z * (other.z);
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/ForceField/UFF/Params.h:29-32:
    // RDKit❗✔️: inline bool isDoubleZero(const double x) {
    // RDKit❗✔️:   return ((x < 1.0e-10) && (x > -1.0e-10));
    // RDKit❗✔️: }
    // RDKit❗✔️: inline void clipToOne(double &x) { x = std::clamp(x, -1.0, 1.0); }
    // Behavior: preserve source vector and normal order, the strict shared
    // zero test and its literal +0.0 return, then clamp only the divided dot.
    // Complexity: constant point arithmetic with four stack vectors and no
    // heap allocation; source read-only point references remain unchanged.
    let r1 = Point3::difference(i, j);
    let r2 = Point3::difference(k, j);
    let r3 = Point3::difference(j, k);
    let r4 = Point3::difference(l, k);
    let t1 = r1.cross_product(&r2);
    let t2 = r3.cross_product(&r4);
    let t1_len = t1.length();
    let t2_len = t2.length();
    if crate::uff::params::is_double_zero(t1_len) || crate::uff::params::is_double_zero(t2_len) {
        return 0.0;
    }
    let mut cos_phi = t1.dot_product(&t2) / (t1_len * t2_len);
    crate::uff::params::clip_to_one(&mut cos_phi);
    cos_phi
}

pub(super) fn calc_oop_chi(i: &Point3, j: &Point3, k: &Point3, l: &Point3) -> f64 {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/ForceField/MMFF/OopBend.cpp:18-34:
    // RDKit❗✔️: double calcOopChi(const RDGeom::Point3D &iPoint, const RDGeom::Point3D &jPoint,
    // RDKit❗✔️:                   const RDGeom::Point3D &kPoint,
    // RDKit❗✔️:                   const RDGeom::Point3D &lPoint) {
    // RDKit❗✔️:   RDGeom::Point3D rJI = iPoint - jPoint;
    // RDKit❗✔️:   RDGeom::Point3D rJK = kPoint - jPoint;
    // RDKit❗✔️:   RDGeom::Point3D rJL = lPoint - jPoint;
    // RDKit❗✔️:   rJI /= rJI.length();
    // RDKit❗✔️:   rJK /= rJK.length();
    // RDKit❗✔️:   rJL /= rJL.length();
    // RDKit❗✔️:
    // RDKit❗✔️:   RDGeom::Point3D n = rJI.crossProduct(rJK);
    // RDKit❗✔️:   n /= n.length();
    // RDKit❗✔️:   double sinChi = n.dotProduct(rJL);
    // RDKit❗✔️:   clipToOne(sinChi);
    // RDKit❗✔️:
    // RDKit❗✔️:   return RAD2DEG * asin(sinChi);
    // RDKit❗✔️: }
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/Geometry/point.cpp:64-70:
    // RDKit❗✔️: Point3D operator-(const Point3D &p1, const Point3D &p2) {
    // RDKit❗✔️:   Point3D res;
    // RDKit❗✔️:   res.x = p1.x - p2.x;
    // RDKit❗✔️:   res.y = p1.y - p2.y;
    // RDKit❗✔️:   res.z = p1.z - p2.z;
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/Geometry/point.h:132-137:
    // RDKit❗✔️: constexpr Point3D &operator/=(double scale) {
    // RDKit❗✔️:   x /= scale;
    // RDKit❗✔️:   y /= scale;
    // RDKit❗✔️:   z /= scale;
    // RDKit❗✔️:   return *this;
    // RDKit❗✔️: }
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/Geometry/point.h:228-234:
    // RDKit❗✔️: constexpr Point3D crossProduct(const Point3D &other) const {
    // RDKit❗✔️:   Point3D res;
    // RDKit❗✔️:   res.x = y * (other.z) - z * (other.y);
    // RDKit❗✔️:   res.y = -x * (other.z) + z * (other.x);
    // RDKit❗✔️:   res.z = x * (other.y) - y * (other.x);
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/Geometry/point.h:158-161:
    // RDKit❗✔️: double length() const override {
    // RDKit❗✔️:   double res = x * x + y * y + z * z;
    // RDKit❗✔️:   return sqrt(res);
    // RDKit❗✔️: }
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/Geometry/point.h:169-172:
    // RDKit❗✔️: constexpr double dotProduct(const Point3D &other) const {
    // RDKit❗✔️:   double res = x * (other.x) + y * (other.y) + z * (other.z);
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/ForceField/MMFF/Params.h:38-44:
    // RDKit❗✔️: constexpr double RAD2DEG = 180.0 / M_PI;
    // RDKit❗✔️: inline void clipToOne(double &x) { x = std::clamp(x, -1.0, 1.0); }
    // Behavior: preserve direct divisions (including zero-length NaNs), the
    // normal/dot/clamp order, source asin and RAD2DEG multiplication. Complexity:
    // fixed point arithmetic, O(1), four stack vectors, no heap allocation.
    let mut r_ji = Point3::difference(i, j);
    let mut r_jk = Point3::difference(k, j);
    let mut r_jl = Point3::difference(l, j);
    let r_ji_length = r_ji.length();
    r_ji.divide_assign(r_ji_length);
    let r_jk_length = r_jk.length();
    r_jk.divide_assign(r_jk_length);
    let r_jl_length = r_jl.length();
    r_jl.divide_assign(r_jl_length);

    let mut normal = r_ji.cross_product(&r_jk);
    let normal_length = normal.length();
    normal.divide_assign(normal_length);
    let mut sin_chi = normal.dot_product(&r_jl);
    crate::uff::params::clip_to_one(&mut sin_chi);
    crate::uff::params::RAD2DEG * sin_chi.asin()
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum MmffGradientError {
    GradientRowOutOfRange {
        slot: usize,
        index: usize,
        rows: usize,
    },
}

pub(super) fn validate_gradient_indices(
    indices: &[usize],
    rows: usize,
) -> Result<(), MmffGradientError> {
    for (slot, index) in indices.iter().copied().enumerate() {
        if index >= rows {
            return Err(MmffGradientError::GradientRowOutOfRange { slot, index, rows });
        }
    }
    Ok(())
}

pub(super) fn calc_angle_bend_grad(
    r: &[Point3; 2],
    dist: &[f64; 2],
    gradient: &mut [[f64; 3]],
    indices: [usize; 3],
    d_e_d_theta: f64,
    cos_theta: f64,
    sin_theta: f64,
) -> Result<(), MmffGradientError> {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/ForceField/MMFF/AngleBend.cpp:59-81:
    // RDKit❗✔️: void calcAngleBendGrad(RDGeom::Point3D *r, double *dist, double **g,
    // RDKit❗✔️:                        double &dE_dTheta, double &cosTheta, double &sinTheta) {
    // RDKit❗✔️:   // -------
    // RDKit❗✔️:   // dTheta/dx is trickier:
    // RDKit❗✔️:   double dCos_dS[6] = {1.0 / dist[0] * (r[1].x - cosTheta * r[0].x),
    // RDKit❗✔️:                        1.0 / dist[0] * (r[1].y - cosTheta * r[0].y),
    // RDKit❗✔️:                        1.0 / dist[0] * (r[1].z - cosTheta * r[0].z),
    // RDKit❗✔️:                        1.0 / dist[1] * (r[0].x - cosTheta * r[1].x),
    // RDKit❗✔️:                        1.0 / dist[1] * (r[0].y - cosTheta * r[1].y),
    // RDKit❗✔️:                        1.0 / dist[1] * (r[0].z - cosTheta * r[1].z)};
    // RDKit❗✔️:
    // RDKit❗✔️:   g[0][0] += dE_dTheta * dCos_dS[0] / (-sinTheta);
    // RDKit❗✔️:   g[0][1] += dE_dTheta * dCos_dS[1] / (-sinTheta);
    // RDKit❗✔️:   g[0][2] += dE_dTheta * dCos_dS[2] / (-sinTheta);
    // RDKit❗✔️:
    // RDKit❗✔️:   g[1][0] += dE_dTheta * (-dCos_dS[0] - dCos_dS[3]) / (-sinTheta);
    // RDKit❗✔️:   g[1][1] += dE_dTheta * (-dCos_dS[1] - dCos_dS[4]) / (-sinTheta);
    // RDKit❗✔️:   g[1][2] += dE_dTheta * (-dCos_dS[2] - dCos_dS[5]) / (-sinTheta);
    // RDKit❗✔️:
    // RDKit❗✔️:   g[2][0] += dE_dTheta * dCos_dS[3] / (-sinTheta);
    // RDKit❗✔️:   g[2][1] += dE_dTheta * dCos_dS[4] / (-sinTheta);
    // RDKit❗✔️:   g[2][2] += dE_dTheta * dCos_dS[5] / (-sinTheta);
    // RDKit❗✔️: }
    // Behavior: preserve the six source derivative expressions and nine
    // sequential accumulator writes. Validate each indexed row in source order
    // before mutation so repeated indices alias one shared accumulator row.
    // Complexity: three fixed bounds checks, six scalar expressions, and nine
    // indexed updates; O(1), no allocation or copied accumulator rows.
    validate_gradient_indices(&indices, gradient.len())?;

    let d_cos_d_s = [
        1.0 / dist[0] * (r[1].x - cos_theta * r[0].x),
        1.0 / dist[0] * (r[1].y - cos_theta * r[0].y),
        1.0 / dist[0] * (r[1].z - cos_theta * r[0].z),
        1.0 / dist[1] * (r[0].x - cos_theta * r[1].x),
        1.0 / dist[1] * (r[0].y - cos_theta * r[1].y),
        1.0 / dist[1] * (r[0].z - cos_theta * r[1].z),
    ];

    gradient[indices[0]][0] += d_e_d_theta * d_cos_d_s[0] / (-sin_theta);
    gradient[indices[0]][1] += d_e_d_theta * d_cos_d_s[1] / (-sin_theta);
    gradient[indices[0]][2] += d_e_d_theta * d_cos_d_s[2] / (-sin_theta);

    gradient[indices[1]][0] += d_e_d_theta * (-d_cos_d_s[0] - d_cos_d_s[3]) / (-sin_theta);
    gradient[indices[1]][1] += d_e_d_theta * (-d_cos_d_s[1] - d_cos_d_s[4]) / (-sin_theta);
    gradient[indices[1]][2] += d_e_d_theta * (-d_cos_d_s[2] - d_cos_d_s[5]) / (-sin_theta);

    gradient[indices[2]][0] += d_e_d_theta * d_cos_d_s[3] / (-sin_theta);
    gradient[indices[2]][1] += d_e_d_theta * d_cos_d_s[4] / (-sin_theta);
    gradient[indices[2]][2] += d_e_d_theta * d_cos_d_s[5] / (-sin_theta);

    Ok(())
}

pub(crate) fn calc_torsion_grad(
    r: &[Point3; 4],
    t: &[Point3; 2],
    d: &[f64; 2],
    gradient: &mut [[f64; 3]],
    indices: [usize; 4],
    sin_term: f64,
    cos_phi: f64,
) -> Result<(), MmffGradientError> {
    // RDKit source pin 351f8f378f8ad6bbd517980c38896e66bf907af8c,
    // Code/ForceField/MMFF/TorsionAngle.cpp:55-93:
    // RDKit❗✔️: void calcTorsionGrad(RDGeom::Point3D *r, RDGeom::Point3D *t, double *d,
    // RDKit❗✔️:                      double **g, double &sinTerm, double &cosPhi) {
    // RDKit❗✔️:   // -------
    // RDKit❗✔️:   // dTheta/dx is trickier:
    // RDKit❗✔️:   double dCos_dT[6] = {1.0 / d[0] * (t[1].x - cosPhi * t[0].x),
    // RDKit❗✔️:                        1.0 / d[0] * (t[1].y - cosPhi * t[0].y),
    // RDKit❗✔️:                        1.0 / d[0] * (t[1].z - cosPhi * t[0].z),
    // RDKit❗✔️:                        1.0 / d[1] * (t[0].x - cosPhi * t[1].x),
    // RDKit❗✔️:                        1.0 / d[1] * (t[0].y - cosPhi * t[1].y),
    // RDKit❗✔️:                        1.0 / d[1] * (t[0].z - cosPhi * t[1].z)};
    // RDKit❗✔️:
    // RDKit❗✔️:   g[0][0] += sinTerm * (dCos_dT[2] * r[1].y - dCos_dT[1] * r[1].z);
    // RDKit❗✔️:   g[0][1] += sinTerm * (dCos_dT[0] * r[1].z - dCos_dT[2] * r[1].x);
    // RDKit❗✔️:   g[0][2] += sinTerm * (dCos_dT[1] * r[1].x - dCos_dT[0] * r[1].y);
    // RDKit❗✔️:
    // RDKit❗✔️:   g[1][0] += sinTerm *
    // RDKit❗✔️:              (dCos_dT[1] * (r[1].z - r[0].z) + dCos_dT[2] * (r[0].y - r[1].y) +
    // RDKit❗✔️:               dCos_dT[4] * (-r[3].z) + dCos_dT[5] * (r[3].y));
    // RDKit❗✔️:   g[1][1] += sinTerm *
    // RDKit❗✔️:              (dCos_dT[0] * (r[0].z - r[1].z) + dCos_dT[2] * (r[1].x - r[0].x) +
    // RDKit❗✔️:               dCos_dT[3] * (r[3].z) + dCos_dT[5] * (-r[3].x));
    // RDKit❗✔️:   g[1][2] += sinTerm *
    // RDKit❗✔️:              (dCos_dT[0] * (r[1].y - r[0].y) + dCos_dT[1] * (r[0].x - r[1].x) +
    // RDKit❗✔️:               dCos_dT[3] * (-r[3].y) + dCos_dT[4] * (r[3].x));
    // RDKit❗✔️:
    // RDKit❗✔️:   g[2][0] += sinTerm *
    // RDKit❗✔️:              (dCos_dT[1] * (r[0].z) + dCos_dT[2] * (-r[0].y) +
    // RDKit❗✔️:               dCos_dT[4] * (r[3].z - r[2].z) + dCos_dT[5] * (r[2].y - r[3].y));
    // RDKit❗✔️:   g[2][1] += sinTerm *
    // RDKit❗✔️:              (dCos_dT[0] * (-r[0].z) + dCos_dT[2] * (r[0].x) +
    // RDKit❗✔️:               dCos_dT[3] * (r[2].z - r[3].z) + dCos_dT[5] * (r[3].x - r[2].x));
    // RDKit❗✔️:   g[2][2] += sinTerm *
    // RDKit❗✔️:              (dCos_dT[0] * (r[0].y) + dCos_dT[1] * (-r[0].x) +
    // RDKit❗✔️:               dCos_dT[3] * (r[3].y - r[2].y) + dCos_dT[4] * (r[2].x - r[3].x));
    // RDKit❗✔️:
    // RDKit❗✔️:   g[3][0] += sinTerm * (dCos_dT[4] * r[2].z - dCos_dT[5] * r[2].y);
    // RDKit❗✔️:   g[3][1] += sinTerm * (dCos_dT[5] * r[2].x - dCos_dT[3] * r[2].z);
    // RDKit❗✔️:   g[3][2] += sinTerm * (dCos_dT[3] * r[2].y - dCos_dT[4] * r[2].x);
    // RDKit❗✔️: }
    // Behavior: preserve the six source derivative expressions and twelve
    // sequential gradient writes in their original component order. The four
    // validated indices address one accumulator, so duplicates alias.
    // Complexity: four fixed bounds checks, fixed scalar arithmetic, and
    // twelve indexed updates; O(1), no allocation or copied rows.
    validate_gradient_indices(&indices, gradient.len())?;

    let d_cos_d_t = [
        1.0 / d[0] * (t[1].x - cos_phi * t[0].x),
        1.0 / d[0] * (t[1].y - cos_phi * t[0].y),
        1.0 / d[0] * (t[1].z - cos_phi * t[0].z),
        1.0 / d[1] * (t[0].x - cos_phi * t[1].x),
        1.0 / d[1] * (t[0].y - cos_phi * t[1].y),
        1.0 / d[1] * (t[0].z - cos_phi * t[1].z),
    ];

    gradient[indices[0]][0] += sin_term * (d_cos_d_t[2] * r[1].y - d_cos_d_t[1] * r[1].z);
    gradient[indices[0]][1] += sin_term * (d_cos_d_t[0] * r[1].z - d_cos_d_t[2] * r[1].x);
    gradient[indices[0]][2] += sin_term * (d_cos_d_t[1] * r[1].x - d_cos_d_t[0] * r[1].y);

    gradient[indices[1]][0] += sin_term
        * (d_cos_d_t[1] * (r[1].z - r[0].z)
            + d_cos_d_t[2] * (r[0].y - r[1].y)
            + d_cos_d_t[4] * (-r[3].z)
            + d_cos_d_t[5] * (r[3].y));
    gradient[indices[1]][1] += sin_term
        * (d_cos_d_t[0] * (r[0].z - r[1].z)
            + d_cos_d_t[2] * (r[1].x - r[0].x)
            + d_cos_d_t[3] * (r[3].z)
            + d_cos_d_t[5] * (-r[3].x));
    gradient[indices[1]][2] += sin_term
        * (d_cos_d_t[0] * (r[1].y - r[0].y)
            + d_cos_d_t[1] * (r[0].x - r[1].x)
            + d_cos_d_t[3] * (-r[3].y)
            + d_cos_d_t[4] * (r[3].x));

    gradient[indices[2]][0] += sin_term
        * (d_cos_d_t[1] * r[0].z
            + d_cos_d_t[2] * (-r[0].y)
            + d_cos_d_t[4] * (r[3].z - r[2].z)
            + d_cos_d_t[5] * (r[2].y - r[3].y));
    gradient[indices[2]][1] += sin_term
        * (d_cos_d_t[0] * (-r[0].z)
            + d_cos_d_t[2] * r[0].x
            + d_cos_d_t[3] * (r[2].z - r[3].z)
            + d_cos_d_t[5] * (r[3].x - r[2].x));
    gradient[indices[2]][2] += sin_term
        * (d_cos_d_t[0] * r[0].y
            + d_cos_d_t[1] * (-r[0].x)
            + d_cos_d_t[3] * (r[3].y - r[2].y)
            + d_cos_d_t[4] * (r[2].x - r[3].x));

    gradient[indices[3]][0] += sin_term * (d_cos_d_t[4] * r[2].z - d_cos_d_t[5] * r[2].y);
    gradient[indices[3]][1] += sin_term * (d_cos_d_t[5] * r[2].x - d_cos_d_t[3] * r[2].z);
    gradient[indices[3]][2] += sin_term * (d_cos_d_t[3] * r[2].y - d_cos_d_t[4] * r[2].x);

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{
        MmffGradientError, calc_angle_bend_grad, calc_angle_force_constant, calc_angle_rest_value,
        calc_bond_force_constant, calc_bond_rest_length, calc_cos_theta, calc_oop_bend_energy,
        calc_oop_bend_force_constant, calc_oop_chi, calc_stbn_force_constants,
        calc_stretch_bend_energy, calc_torsion_cos_phi, calc_torsion_energy,
        calc_torsion_force_constant, calc_torsion_grad,
    };
    use crate::geometry::Point3;
    use crate::mmff::params::{MmffAngle, MmffBond, MmffOop, MmffStbn, MmffTor};

    const SOURCE: &str = include_str!(
        "../../../../testdata/forcefields/expected/rdkit/mmff_numerical/source_literals.tsv"
    );
    const V_BITS: [u64; 3] = [
        0x8000_0000_0000_0000,
        0x3fe0_0000_0000_0000,
        0x4010_0000_0000_0000,
    ];
    const GRADIENT_BASELINES_BITS: [[[u64; 3]; 7]; 2] = [
        [[0; 3]; 7],
        [
            [
                0x3fd0_0000_0000_0000,
                0xbfe0_0000_0000_0000,
                0x8000_0000_0000_0000,
            ],
            [
                0x3fe0_0000_0000_0000,
                0xbff0_0000_0000_0000,
                0x8000_0000_0000_0000,
            ],
            [
                0x3fe8_0000_0000_0000,
                0xbff8_0000_0000_0000,
                0x8000_0000_0000_0000,
            ],
            [
                0x3ff0_0000_0000_0000,
                0xc000_0000_0000_0000,
                0x8000_0000_0000_0000,
            ],
            [
                0x3ff4_0000_0000_0000,
                0xc004_0000_0000_0000,
                0x8000_0000_0000_0000,
            ],
            [
                0x3ff8_0000_0000_0000,
                0xc008_0000_0000_0000,
                0x8000_0000_0000_0000,
            ],
            [
                0x3ffc_0000_0000_0000,
                0xc00c_0000_0000_0000,
                0x8000_0000_0000_0000,
            ],
        ],
    ];

    #[derive(Default)]
    struct Audit {
        discrepancies: Vec<String>,
        calls: usize,
        cells: usize,
        error_controls: usize,
    }

    impl Audit {
        fn record(&mut self, message: impl Into<String>) {
            self.discrepancies.push(message.into());
        }

        fn check_getter<T, Getter, Snapshot>(
            &mut self,
            label: String,
            input: &T,
            getter: Getter,
            snapshot: Snapshot,
        ) where
            Getter: Fn(&T) -> Vec<u64>,
            Snapshot: Fn(&T) -> Vec<u64>,
        {
            let expected = SOURCE.lines().find_map(|line| {
                let mut fields = line.split('\t');
                if fields.next() != Some("U1") || fields.next() != Some(label.as_str()) {
                    return None;
                }
                fields
                    .map(|field| u64::from_str_radix(field, 16).ok())
                    .collect::<Option<Vec<_>>>()
            });
            if expected.is_none() {
                self.record(format!("{label}: missing or invalid frozen U1 row"));
            }

            let frozen_input = snapshot(input);
            let frozen_address = input as *const T as usize;
            let mut previous_output: Option<Vec<u64>> = None;
            for repeat in 0..2 {
                let input_ref = input;
                let address_before = input_ref as *const T as usize;
                let actual_output = getter(input_ref);
                self.calls += 1;

                let address_after = input_ref as *const T as usize;
                let after_input = snapshot(input_ref);
                if after_input != frozen_input {
                    self.record(format!("{label} repeat {repeat}: input fields changed"));
                }
                if address_before != frozen_address || address_after != frozen_address {
                    self.record(format!(
                        "{label} repeat {repeat}: borrowed row address changed"
                    ));
                }

                if let Some(expected_output) = expected.as_ref() {
                    if actual_output.len() != expected_output.len() {
                        self.record(format!(
                            "{label} repeat {repeat}: result width {}, expected {}",
                            actual_output.len(),
                            expected_output.len()
                        ));
                    } else {
                        for (component, (actual, expected)) in
                            actual_output.iter().zip(expected_output).enumerate()
                        {
                            if actual != expected {
                                self.record(format!(
                                    "{label} repeat {repeat} component {component}: {actual:016x}, expected {expected:016x}"
                                ));
                            }
                        }
                    }
                }

                if let Some(previous) = previous_output.as_ref() {
                    if &actual_output != previous {
                        self.record(format!("{label}: repeated result bits differ"));
                    }
                } else {
                    previous_output = Some(actual_output);
                }
            }
            self.cells += 1;
        }

        fn check_stretch_bend_energy(
            &mut self,
            label: String,
            delta_dist1: f64,
            delta_dist2: f64,
            delta_theta: f64,
            force_constants: (f64, f64),
        ) {
            let expected = SOURCE.lines().find_map(|line| {
                let mut fields = line.split('\t');
                if fields.next() != Some("U2") || fields.next() != Some(label.as_str()) {
                    return None;
                }
                fields
                    .map(|field| u64::from_str_radix(field, 16).ok())
                    .collect::<Option<Vec<_>>>()
            });
            if expected.is_none() {
                self.record(format!("{label}: missing or invalid frozen U2 row"));
            }

            let input_before = [
                delta_dist1.to_bits(),
                delta_dist2.to_bits(),
                delta_theta.to_bits(),
                force_constants.0.to_bits(),
                force_constants.1.to_bits(),
            ];
            let mut previous_output: Option<[u64; 2]> = None;
            for repeat in 0..2 {
                let actual = calc_stretch_bend_energy(
                    delta_dist1,
                    delta_dist2,
                    delta_theta,
                    force_constants,
                );
                self.calls += 1;

                let input_after = [
                    delta_dist1.to_bits(),
                    delta_dist2.to_bits(),
                    delta_theta.to_bits(),
                    force_constants.0.to_bits(),
                    force_constants.1.to_bits(),
                ];
                if input_after != input_before {
                    self.record(format!("{label} repeat {repeat}: input bits changed"));
                }

                let actual_bits = [actual.0.to_bits(), actual.1.to_bits()];
                if let Some(expected_output) = expected.as_ref() {
                    if expected_output.len() != actual_bits.len() {
                        self.record(format!(
                            "{label} repeat {repeat}: result width {}, expected {}",
                            actual_bits.len(),
                            expected_output.len()
                        ));
                    } else {
                        for (component, (actual, expected)) in
                            actual_bits.iter().zip(expected_output).enumerate()
                        {
                            if actual != expected {
                                self.record(format!(
                                    "{label} repeat {repeat} component {component}: {actual:016x}, expected {expected:016x}"
                                ));
                            }
                        }
                    }
                }

                if let Some(previous) = previous_output {
                    if actual_bits != previous {
                        self.record(format!("{label}: repeated result bits differ"));
                    }
                } else {
                    previous_output = Some(actual_bits);
                }
            }
            self.cells += 1;
        }

        fn check_torsion_energy(&mut self, label: String, v1: f64, v2: f64, v3: f64, cos_phi: f64) {
            let expected = SOURCE.lines().find_map(|line| {
                let mut fields = line.split('\t');
                if fields.next() != Some("U3") || fields.next() != Some(label.as_str()) {
                    return None;
                }
                fields
                    .map(|field| u64::from_str_radix(field, 16).ok())
                    .collect::<Option<Vec<_>>>()
            });
            if expected.is_none() {
                self.record(format!("{label}: missing or invalid frozen U3 row"));
            }

            let input_before = [v1.to_bits(), v2.to_bits(), v3.to_bits(), cos_phi.to_bits()];
            let mut previous_output: Option<u64> = None;
            for repeat in 0..2 {
                let actual = calc_torsion_energy(v1, v2, v3, cos_phi);
                self.calls += 1;

                let input_after = [v1.to_bits(), v2.to_bits(), v3.to_bits(), cos_phi.to_bits()];
                if input_after != input_before {
                    self.record(format!("{label} repeat {repeat}: input bits changed"));
                }

                let actual_bits = actual.to_bits();
                if let Some(expected_output) = expected.as_ref() {
                    if expected_output.len() != 1 {
                        self.record(format!(
                            "{label} repeat {repeat}: result width 1, expected {}",
                            expected_output.len()
                        ));
                    } else if actual_bits != expected_output[0] {
                        self.record(format!(
                            "{label} repeat {repeat}: {actual_bits:016x}, expected {:016x}",
                            expected_output[0]
                        ));
                    }
                }

                if let Some(previous) = previous_output {
                    if actual_bits != previous {
                        self.record(format!("{label}: repeated result bits differ"));
                    }
                } else {
                    previous_output = Some(actual_bits);
                }
            }
            self.cells += 1;
        }

        fn check_oop_bend_energy(&mut self, label: String, chi: f64, koop: f64) {
            let expected = SOURCE.lines().find_map(|line| {
                let mut fields = line.split('\t');
                if fields.next() != Some("U4") || fields.next() != Some(label.as_str()) {
                    return None;
                }
                fields
                    .map(|field| u64::from_str_radix(field, 16).ok())
                    .collect::<Option<Vec<_>>>()
            });
            if expected.is_none() {
                self.record(format!("{label}: missing or invalid frozen U4 row"));
            }

            let input_before = [chi.to_bits(), koop.to_bits()];
            let mut previous_output: Option<u64> = None;
            for repeat in 0..2 {
                let actual = calc_oop_bend_energy(chi, koop);
                self.calls += 1;

                let input_after = [chi.to_bits(), koop.to_bits()];
                if input_after != input_before {
                    self.record(format!("{label} repeat {repeat}: input bits changed"));
                }

                let actual_bits = actual.to_bits();
                if let Some(expected_output) = expected.as_ref() {
                    if expected_output.len() != 1 {
                        self.record(format!(
                            "{label} repeat {repeat}: result width 1, expected {}",
                            expected_output.len()
                        ));
                    } else if actual_bits != expected_output[0] {
                        self.record(format!(
                            "{label} repeat {repeat}: {actual_bits:016x}, expected {:016x}",
                            expected_output[0]
                        ));
                    }
                }

                if let Some(previous) = previous_output {
                    if actual_bits != previous {
                        self.record(format!("{label}: repeated result bits differ"));
                    }
                } else {
                    previous_output = Some(actual_bits);
                }
            }
            self.cells += 1;
        }

        fn check_cos_theta(
            &mut self,
            label: String,
            p1: Point3,
            p2: Point3,
            p3: Point3,
            dist1: f64,
            dist2: f64,
        ) {
            let expected = SOURCE.lines().find_map(|line| {
                let mut fields = line.split('\t');
                if fields.next() != Some("U5") || fields.next() != Some(label.as_str()) {
                    return None;
                }
                fields
                    .map(|field| u64::from_str_radix(field, 16).ok())
                    .collect::<Option<Vec<_>>>()
            });
            if expected.is_none() {
                self.record(format!("{label}: missing or invalid frozen U5 row"));
            }

            let point_bits =
                |point: Point3| [point.x.to_bits(), point.y.to_bits(), point.z.to_bits()];
            let input_before = [
                point_bits(p1),
                point_bits(p2),
                point_bits(p3),
                [dist1.to_bits(), dist2.to_bits(), 0],
            ];
            let mut previous_output: Option<u64> = None;
            for repeat in 0..2 {
                let actual = calc_cos_theta(p1, p2, p3, dist1, dist2);
                self.calls += 1;

                let input_after = [
                    point_bits(p1),
                    point_bits(p2),
                    point_bits(p3),
                    [dist1.to_bits(), dist2.to_bits(), 0],
                ];
                if input_after != input_before {
                    self.record(format!("{label} repeat {repeat}: input bits changed"));
                }

                let actual_bits = actual.to_bits();
                if let Some(expected_output) = expected.as_ref() {
                    if expected_output.len() != 1 {
                        self.record(format!(
                            "{label} repeat {repeat}: result width 1, expected {}",
                            expected_output.len()
                        ));
                    } else if actual_bits != expected_output[0] {
                        self.record(format!(
                            "{label} repeat {repeat}: {actual_bits:016x}, expected {:016x}",
                            expected_output[0]
                        ));
                    }
                }

                if let Some(previous) = previous_output {
                    if actual_bits != previous {
                        self.record(format!("{label}: repeated result bits differ"));
                    }
                } else {
                    previous_output = Some(actual_bits);
                }
            }
            self.cells += 1;
        }

        fn check_torsion_cos_phi(&mut self, label: String, points: [&Point3; 4]) {
            let expected = SOURCE.lines().find_map(|line| {
                let mut fields = line.split('\t');
                if fields.next() != Some("U6") || fields.next() != Some(label.as_str()) {
                    return None;
                }
                fields
                    .map(|field| u64::from_str_radix(field, 16).ok())
                    .collect::<Option<Vec<_>>>()
            });
            if expected.is_none() {
                self.record(format!("{label}: missing or invalid frozen U6 row"));
            }

            let point_bits =
                |point: &Point3| [point.x.to_bits(), point.y.to_bits(), point.z.to_bits()];
            let input_before = points.map(point_bits);
            let address_before = points.map(|point| point as *const Point3 as usize);
            let mut previous_output: Option<u64> = None;
            for repeat in 0..2 {
                let actual = calc_torsion_cos_phi(points[0], points[1], points[2], points[3]);
                self.calls += 1;

                let input_after = points.map(point_bits);
                let address_after = points.map(|point| point as *const Point3 as usize);
                if input_after != input_before {
                    self.record(format!("{label} repeat {repeat}: point bits changed"));
                }
                if address_after != address_before {
                    self.record(format!(
                        "{label} repeat {repeat}: borrowed point identity changed"
                    ));
                }

                let actual_bits = actual.to_bits();
                if let Some(expected_output) = expected.as_ref() {
                    if expected_output.len() != 1 {
                        self.record(format!(
                            "{label} repeat {repeat}: result width 1, expected {}",
                            expected_output.len()
                        ));
                    } else if actual_bits != expected_output[0] {
                        self.record(format!(
                            "{label} repeat {repeat}: {actual_bits:016x}, expected {:016x}",
                            expected_output[0]
                        ));
                    }
                }

                if let Some(previous) = previous_output {
                    if actual_bits != previous {
                        self.record(format!("{label}: repeated result bits differ"));
                    }
                } else {
                    previous_output = Some(actual_bits);
                }
            }
            self.cells += 1;
        }

        fn check_oop_chi(&mut self, label: String, points: [&Point3; 4]) {
            let expected = SOURCE.lines().find_map(|line| {
                let mut fields = line.split('\t');
                if fields.next() != Some("U7") || fields.next() != Some(label.as_str()) {
                    return None;
                }
                fields
                    .map(|field| u64::from_str_radix(field, 16).ok())
                    .collect::<Option<Vec<_>>>()
            });
            if expected.is_none() {
                self.record(format!("{label}: missing or invalid frozen U7 row"));
            }

            let point_bits =
                |point: &Point3| [point.x.to_bits(), point.y.to_bits(), point.z.to_bits()];
            let input_before = points.map(point_bits);
            let address_before = points.map(|point| point as *const Point3 as usize);
            let mut previous_output: Option<u64> = None;
            for repeat in 0..2 {
                let actual = calc_oop_chi(points[0], points[1], points[2], points[3]);
                self.calls += 1;

                let input_after = points.map(point_bits);
                let address_after = points.map(|point| point as *const Point3 as usize);
                if input_after != input_before {
                    self.record(format!("{label} repeat {repeat}: point bits changed"));
                }
                if address_after != address_before {
                    self.record(format!(
                        "{label} repeat {repeat}: borrowed point identity changed"
                    ));
                }

                let actual_bits = actual.to_bits();
                if let Some(expected_output) = expected.as_ref() {
                    if expected_output.len() != 1 {
                        self.record(format!(
                            "{label} repeat {repeat}: result width 1, expected {}",
                            expected_output.len()
                        ));
                    } else if actual_bits != expected_output[0] {
                        self.record(format!(
                            "{label} repeat {repeat}: {actual_bits:016x}, expected {:016x}",
                            expected_output[0]
                        ));
                    }
                }

                if let Some(previous) = previous_output {
                    if actual_bits != previous {
                        self.record(format!("{label}: repeated result bits differ"));
                    }
                } else {
                    previous_output = Some(actual_bits);
                }
            }
            self.cells += 1;
        }

        fn check_angle_bend_gradient(
            &mut self,
            label: String,
            r: &[Point3; 2],
            dist: &[f64; 2],
            d_e_d_theta: f64,
            cos_theta: f64,
            sin_theta: f64,
            indices: [usize; 3],
            baseline_bits: [[u64; 3]; 7],
        ) {
            let expected = SOURCE.lines().find_map(|line| {
                let mut fields = line.split('\t');
                if fields.next() != Some("U8") || fields.next() != Some(label.as_str()) {
                    return None;
                }
                fields
                    .map(|field| u64::from_str_radix(field, 16).ok())
                    .collect::<Option<Vec<_>>>()
            });
            if expected.is_none() {
                self.record(format!("{label}: missing or invalid frozen U8 row"));
            }

            let point_bits =
                |point: &Point3| [point.x.to_bits(), point.y.to_bits(), point.z.to_bits()];
            let r_bits_before = r.iter().map(point_bits).collect::<Vec<_>>();
            let r_addresses_before = r
                .iter()
                .map(|point| point as *const Point3 as usize)
                .collect::<Vec<_>>();
            let dist_bits_before = dist.map(f64::to_bits);
            let scalar_bits_before = [
                d_e_d_theta.to_bits(),
                cos_theta.to_bits(),
                sin_theta.to_bits(),
            ];
            let indices_before = indices;
            let mut previous_output: Option<Vec<u64>> = None;

            for repeat in 0..2 {
                let mut gradient = baseline_bits.map(|row| row.map(f64::from_bits));
                let baseline_before = gradient.map(|row| row.map(f64::to_bits));
                let gradient_address_before = gradient.as_ptr() as usize;
                let r_address_before = r.as_ptr() as usize;
                let dist_address_before = dist.as_ptr() as usize;
                let result = calc_angle_bend_grad(
                    r,
                    dist,
                    &mut gradient,
                    indices,
                    d_e_d_theta,
                    cos_theta,
                    sin_theta,
                );
                self.calls += 1;

                let r_bits_after = r.iter().map(point_bits).collect::<Vec<_>>();
                let r_addresses_after = r
                    .iter()
                    .map(|point| point as *const Point3 as usize)
                    .collect::<Vec<_>>();
                let dist_bits_after = dist.map(f64::to_bits);
                let scalar_bits_after = [
                    d_e_d_theta.to_bits(),
                    cos_theta.to_bits(),
                    sin_theta.to_bits(),
                ];
                if r_bits_after != r_bits_before
                    || dist_bits_after != dist_bits_before
                    || scalar_bits_after != scalar_bits_before
                    || indices != indices_before
                {
                    self.record(format!("{label} repeat {repeat}: read-only inputs changed"));
                }
                if r_address_before != r.as_ptr() as usize
                    || r_addresses_after != r_addresses_before
                    || dist_address_before != dist.as_ptr() as usize
                    || gradient_address_before != gradient.as_ptr() as usize
                {
                    self.record(format!(
                        "{label} repeat {repeat}: borrowed storage identity changed"
                    ));
                }
                if baseline_before != baseline_bits {
                    self.record(format!(
                        "{label} repeat {repeat}: baseline bits changed before call"
                    ));
                }

                if !matches!(result, Ok(())) {
                    self.record(format!(
                        "{label} repeat {repeat}: valid source input returned {result:?}"
                    ));
                    continue;
                }

                let actual = gradient
                    .iter()
                    .flat_map(|row| row.iter().map(|value| value.to_bits()))
                    .collect::<Vec<_>>();
                if let Some(expected_output) = expected.as_ref() {
                    if expected_output.len() != 21 {
                        self.record(format!(
                            "{label} repeat {repeat}: result width {}, expected 21",
                            expected_output.len()
                        ));
                    } else if &actual != expected_output {
                        for (component, (actual_bits, expected_bits)) in
                            actual.iter().zip(expected_output).enumerate()
                        {
                            if actual_bits != expected_bits {
                                self.record(format!(
                                    "{label} repeat {repeat} component {component}: {actual_bits:016x}, expected {expected_bits:016x}"
                                ));
                            }
                        }
                    }
                }
                if let Some(previous) = previous_output.as_ref() {
                    if &actual != previous {
                        self.record(format!("{label}: repeated full-gradient bits differ"));
                    }
                } else {
                    previous_output = Some(actual);
                }
            }
            self.cells += 1;
        }

        fn check_angle_bend_error(
            &mut self,
            label: &str,
            indices: [usize; 3],
            expected_error: MmffGradientError,
            empty_gradient: bool,
        ) {
            let r = [
                Point3 {
                    x: 1.0,
                    y: 0.0,
                    z: 0.0,
                },
                Point3 {
                    x: 0.0,
                    y: 1.0,
                    z: 0.0,
                },
            ];
            let dist = [1.0, 2.0];
            let d_e_d_theta: f64 = 1.25;
            let cos_theta: f64 = -0.5;
            let sin_theta: f64 = -0.8;
            let mut gradient = if empty_gradient {
                Vec::new()
            } else {
                GRADIENT_BASELINES_BITS[1]
                    .map(|row| row.map(f64::from_bits))
                    .to_vec()
            };
            let r_bits_before =
                r.map(|point| [point.x.to_bits(), point.y.to_bits(), point.z.to_bits()]);
            let dist_bits_before = dist.map(f64::to_bits);
            let scalar_bits_before = [
                d_e_d_theta.to_bits(),
                cos_theta.to_bits(),
                sin_theta.to_bits(),
            ];
            let gradient_bits_before = gradient
                .iter()
                .map(|row| row.map(f64::to_bits))
                .collect::<Vec<_>>();
            let gradient_address_before = gradient.as_ptr() as usize;
            let r_address_before = r.as_ptr() as usize;
            let dist_address_before = dist.as_ptr() as usize;
            let result = calc_angle_bend_grad(
                &r,
                &dist,
                &mut gradient,
                indices,
                d_e_d_theta,
                cos_theta,
                sin_theta,
            );
            self.error_controls += 1;

            let r_bits_after =
                r.map(|point| [point.x.to_bits(), point.y.to_bits(), point.z.to_bits()]);
            let dist_bits_after = dist.map(f64::to_bits);
            let scalar_bits_after = [
                d_e_d_theta.to_bits(),
                cos_theta.to_bits(),
                sin_theta.to_bits(),
            ];
            let gradient_bits_after = gradient
                .iter()
                .map(|row| row.map(f64::to_bits))
                .collect::<Vec<_>>();
            if r_bits_after != r_bits_before
                || dist_bits_after != dist_bits_before
                || scalar_bits_after != scalar_bits_before
            {
                self.record(format!("{label}: read-only inputs changed on bounds error"));
            }
            if gradient_bits_after != gradient_bits_before
                || gradient.as_ptr() as usize != gradient_address_before
                || r.as_ptr() as usize != r_address_before
                || dist.as_ptr() as usize != dist_address_before
            {
                self.record(format!(
                    "{label}: buffer or borrowed storage changed on bounds error"
                ));
            }
            match result {
                Err(actual_error) if actual_error == expected_error => {}
                Err(actual_error) => self.record(format!(
                    "{label}: got {actual_error:?}, expected {expected_error:?}"
                )),
                Ok(()) => self.record(format!("{label}: expected {expected_error:?}, got Ok")),
            }
        }

        fn check_torsion_gradient(
            &mut self,
            label: String,
            r: &[Point3; 4],
            t: &[Point3; 2],
            d: &[f64; 2],
            sin_term: f64,
            cos_phi: f64,
            indices: [usize; 4],
            baseline_bits: [[u64; 3]; 7],
        ) {
            let expected = SOURCE.lines().find_map(|line| {
                let mut fields = line.split('\t');
                if fields.next() != Some("U9") || fields.next() != Some(label.as_str()) {
                    return None;
                }
                fields
                    .map(|field| u64::from_str_radix(field, 16).ok())
                    .collect::<Option<Vec<_>>>()
            });
            if expected.is_none() {
                self.record(format!("{label}: missing or invalid frozen U9 row"));
            }

            let point_bits =
                |point: &Point3| [point.x.to_bits(), point.y.to_bits(), point.z.to_bits()];
            let r_bits_before = r.iter().map(point_bits).collect::<Vec<_>>();
            let t_bits_before = t.iter().map(point_bits).collect::<Vec<_>>();
            let r_addresses_before = r
                .iter()
                .map(|point| point as *const Point3 as usize)
                .collect::<Vec<_>>();
            let t_addresses_before = t
                .iter()
                .map(|point| point as *const Point3 as usize)
                .collect::<Vec<_>>();
            let d_bits_before = d.map(f64::to_bits);
            let scalar_bits_before = [sin_term.to_bits(), cos_phi.to_bits()];
            let indices_before = indices;
            let mut previous_output: Option<Vec<u64>> = None;

            for repeat in 0..2 {
                let mut gradient = baseline_bits.map(|row| row.map(f64::from_bits));
                let baseline_before = gradient.map(|row| row.map(f64::to_bits));
                let gradient_address_before = gradient.as_ptr() as usize;
                let r_address_before = r.as_ptr() as usize;
                let t_address_before = t.as_ptr() as usize;
                let d_address_before = d.as_ptr() as usize;
                let result = calc_torsion_grad(r, t, d, &mut gradient, indices, sin_term, cos_phi);
                self.calls += 1;

                let r_bits_after = r.iter().map(point_bits).collect::<Vec<_>>();
                let t_bits_after = t.iter().map(point_bits).collect::<Vec<_>>();
                let r_addresses_after = r
                    .iter()
                    .map(|point| point as *const Point3 as usize)
                    .collect::<Vec<_>>();
                let t_addresses_after = t
                    .iter()
                    .map(|point| point as *const Point3 as usize)
                    .collect::<Vec<_>>();
                let d_bits_after = d.map(f64::to_bits);
                let scalar_bits_after = [sin_term.to_bits(), cos_phi.to_bits()];
                if r_bits_after != r_bits_before
                    || t_bits_after != t_bits_before
                    || d_bits_after != d_bits_before
                    || scalar_bits_after != scalar_bits_before
                    || indices != indices_before
                {
                    self.record(format!("{label} repeat {repeat}: read-only inputs changed"));
                }
                if r_address_before != r.as_ptr() as usize
                    || t_address_before != t.as_ptr() as usize
                    || d_address_before != d.as_ptr() as usize
                    || r_addresses_after != r_addresses_before
                    || t_addresses_after != t_addresses_before
                    || gradient_address_before != gradient.as_ptr() as usize
                {
                    self.record(format!(
                        "{label} repeat {repeat}: borrowed storage identity changed"
                    ));
                }
                if baseline_before != baseline_bits {
                    self.record(format!(
                        "{label} repeat {repeat}: baseline bits changed before call"
                    ));
                }
                if !matches!(result, Ok(())) {
                    self.record(format!(
                        "{label} repeat {repeat}: valid source input returned {result:?}"
                    ));
                    continue;
                }

                let actual = gradient
                    .iter()
                    .flat_map(|row| row.iter().map(|value| value.to_bits()))
                    .collect::<Vec<_>>();
                if let Some(expected_output) = expected.as_ref() {
                    if expected_output.len() != 21 {
                        self.record(format!(
                            "{label} repeat {repeat}: result width {}, expected 21",
                            expected_output.len()
                        ));
                    } else if &actual != expected_output {
                        for (component, (actual_bits, expected_bits)) in
                            actual.iter().zip(expected_output).enumerate()
                        {
                            if actual_bits != expected_bits {
                                self.record(format!(
                                    "{label} repeat {repeat} component {component}: {actual_bits:016x}, expected {expected_bits:016x}"
                                ));
                            }
                        }
                    }
                }
                if let Some(previous) = previous_output.as_ref() {
                    if &actual != previous {
                        self.record(format!("{label}: repeated full-gradient bits differ"));
                    }
                } else {
                    previous_output = Some(actual);
                }
            }
            self.cells += 1;
        }

        fn check_torsion_gradient_error(
            &mut self,
            label: &str,
            indices: [usize; 4],
            expected_error: MmffGradientError,
            empty_gradient: bool,
        ) {
            let r = [
                Point3 {
                    x: 1.0,
                    y: 0.0,
                    z: 0.0,
                },
                Point3 {
                    x: 0.0,
                    y: 1.0,
                    z: 0.0,
                },
                Point3 {
                    x: 0.25,
                    y: -0.5,
                    z: 1.0,
                },
                Point3 {
                    x: -1.25,
                    y: 0.75,
                    z: 2.0,
                },
            ];
            let t = [r[0], r[1]];
            let d = [1.0, 2.0];
            let sin_term: f64 = 1.25;
            let cos_phi: f64 = -0.5;
            let mut gradient = if empty_gradient {
                Vec::new()
            } else {
                GRADIENT_BASELINES_BITS[1]
                    .map(|row| row.map(f64::from_bits))
                    .to_vec()
            };
            let point_bits =
                |point: &Point3| [point.x.to_bits(), point.y.to_bits(), point.z.to_bits()];
            let r_bits_before = r.map(|point| point_bits(&point));
            let t_bits_before = t.map(|point| point_bits(&point));
            let d_bits_before = d.map(f64::to_bits);
            let scalar_bits_before = [sin_term.to_bits(), cos_phi.to_bits()];
            let indices_before = indices;
            let gradient_bits_before = gradient
                .iter()
                .map(|row| row.map(f64::to_bits))
                .collect::<Vec<_>>();
            let gradient_address_before = gradient.as_ptr() as usize;
            let r_address_before = r.as_ptr() as usize;
            let t_address_before = t.as_ptr() as usize;
            let d_address_before = d.as_ptr() as usize;
            let result = calc_torsion_grad(&r, &t, &d, &mut gradient, indices, sin_term, cos_phi);
            self.error_controls += 1;

            let r_bits_after = r.map(|point| point_bits(&point));
            let t_bits_after = t.map(|point| point_bits(&point));
            let d_bits_after = d.map(f64::to_bits);
            let scalar_bits_after = [sin_term.to_bits(), cos_phi.to_bits()];
            let gradient_bits_after = gradient
                .iter()
                .map(|row| row.map(f64::to_bits))
                .collect::<Vec<_>>();
            if r_bits_after != r_bits_before
                || t_bits_after != t_bits_before
                || d_bits_after != d_bits_before
                || scalar_bits_after != scalar_bits_before
                || indices != indices_before
            {
                self.record(format!("{label}: read-only inputs changed on bounds error"));
            }
            if gradient_bits_after != gradient_bits_before
                || gradient.as_ptr() as usize != gradient_address_before
                || r.as_ptr() as usize != r_address_before
                || t.as_ptr() as usize != t_address_before
                || d.as_ptr() as usize != d_address_before
            {
                self.record(format!(
                    "{label}: buffer or borrowed storage changed on bounds error"
                ));
            }
            match result {
                Err(actual_error) if actual_error == expected_error => {}
                Err(actual_error) => self.record(format!(
                    "{label}: got {actual_error:?}, expected {expected_error:?}"
                )),
                Ok(()) => self.record(format!("{label}: expected {expected_error:?}, got Ok")),
            }
        }
    }

    fn scalar(bits: u64) -> f64 {
        f64::from_bits(bits)
    }

    #[test]
    fn mmff_numerical_u1_getters_match_source_bits_and_preserve_borrows() {
        let mut audit = Audit::default();
        let values = V_BITS.map(scalar);

        let mut ordinal = 0;
        for r0 in values {
            for kb in values {
                let params = MmffBond { kb, r0 };
                audit.check_getter(
                    format!("bond_r0_{ordinal:02}"),
                    &params,
                    |value| vec![calc_bond_rest_length(value).to_bits()],
                    |value| vec![value.kb.to_bits(), value.r0.to_bits()],
                );
                audit.check_getter(
                    format!("bond_kb_{ordinal:02}"),
                    &params,
                    |value| vec![calc_bond_force_constant(value).to_bits()],
                    |value| vec![value.kb.to_bits(), value.r0.to_bits()],
                );
                ordinal += 1;
            }
        }

        ordinal = 0;
        for theta0 in values {
            for ka in values {
                let params = MmffAngle { ka, theta0 };
                audit.check_getter(
                    format!("angle_theta0_{ordinal:02}"),
                    &params,
                    |value| vec![calc_angle_rest_value(value).to_bits()],
                    |value| vec![value.ka.to_bits(), value.theta0.to_bits()],
                );
                audit.check_getter(
                    format!("angle_ka_{ordinal:02}"),
                    &params,
                    |value| vec![calc_angle_force_constant(value).to_bits()],
                    |value| vec![value.ka.to_bits(), value.theta0.to_bits()],
                );
                ordinal += 1;
            }
        }

        ordinal = 0;
        for first in values {
            for second in values {
                let params = MmffStbn {
                    kba_ijk: first,
                    kba_kji: second,
                };
                audit.check_getter(
                    format!("stbn_{ordinal:02}"),
                    &params,
                    |value| {
                        let (first, second) = calc_stbn_force_constants(value);
                        vec![first.to_bits(), second.to_bits()]
                    },
                    |value| vec![value.kba_ijk.to_bits(), value.kba_kji.to_bits()],
                );
                ordinal += 1;
            }
        }

        ordinal = 0;
        for v1 in values {
            for v2 in values {
                for v3 in values {
                    let params = MmffTor { v1, v2, v3 };
                    audit.check_getter(
                        format!("tor_{ordinal:02}"),
                        &params,
                        |value| {
                            let (v1, v2, v3) = calc_torsion_force_constant(value);
                            vec![v1.to_bits(), v2.to_bits(), v3.to_bits()]
                        },
                        |value| vec![value.v1.to_bits(), value.v2.to_bits(), value.v3.to_bits()],
                    );
                    ordinal += 1;
                }
            }
        }

        ordinal = 0;
        for koop in values {
            let params = MmffOop { koop };
            audit.check_getter(
                format!("oop_{ordinal:02}"),
                &params,
                |value| vec![calc_oop_bend_force_constant(value).to_bits()],
                |value| vec![value.koop.to_bits()],
            );
            ordinal += 1;
        }

        if audit.cells != 75 {
            audit.record(format!("expected 75 source cells, got {}", audit.cells));
        }
        if audit.calls != 150 {
            audit.record(format!(
                "expected 150 actual Rust calls, got {}",
                audit.calls
            ));
        }
        assert!(
            audit.discrepancies.is_empty(),
            "MMFF U1 getter discrepancies:\n{}",
            audit.discrepancies.join("\n")
        );
    }

    #[test]
    fn mmff_numerical_u2_stretch_bend_energy_matches_source_bits() {
        const D1_BITS: [u64; 3] = [
            0xbfc0_0000_0000_0000,
            0x0000_0000_0000_0000,
            0x3fe8_0000_0000_0000,
        ];
        const D2_BITS: [u64; 3] = [
            0xbfe0_0000_0000_0000,
            0x0000_0000_0000_0000,
            0x3ff4_0000_0000_0000,
        ];
        const THETA_BITS: [u64; 3] = [
            0xc03e_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x4046_8000_0000_0000,
        ];
        const K1_BITS: [u64; 2] = [0x8000_0000_0000_0000, 0x3fe0_0000_0000_0000];
        const K2_BITS: [u64; 2] = [0xbfe8_0000_0000_0000, 0x3ff4_0000_0000_0000];

        let mut audit = Audit::default();
        let mut ordinal = 0;
        for delta_dist1 in D1_BITS.map(f64::from_bits) {
            for delta_dist2 in D2_BITS.map(f64::from_bits) {
                for delta_theta in THETA_BITS.map(f64::from_bits) {
                    for k1 in K1_BITS.map(f64::from_bits) {
                        for k2 in K2_BITS.map(f64::from_bits) {
                            audit.check_stretch_bend_energy(
                                format!("stretch_bend_{ordinal:03}"),
                                delta_dist1,
                                delta_dist2,
                                delta_theta,
                                (k1, k2),
                            );
                            ordinal += 1;
                        }
                    }
                }
            }
        }

        if audit.cells != 108 {
            audit.record(format!("expected 108 source cells, got {}", audit.cells));
        }
        if audit.calls != 216 {
            audit.record(format!(
                "expected 216 actual Rust calls, got {}",
                audit.calls
            ));
        }
        assert!(
            audit.discrepancies.is_empty(),
            "MMFF U2 stretch-bend energy discrepancies:\n{}",
            audit.discrepancies.join("\n")
        );
    }

    #[test]
    fn mmff_numerical_u3_torsion_energy_matches_source_bits() {
        const V1_BITS: [u64; 2] = [0x8000_0000_0000_0000, 0x3fe0_0000_0000_0000];
        const V2_BITS: [u64; 2] = [0xbff4_0000_0000_0000, 0x3fe8_0000_0000_0000];
        const V3_BITS: [u64; 2] = [0x0000_0000_0000_0000, 0x4000_0000_0000_0000];
        const COS_PHI_BITS: [u64; 6] = [
            0xbff0_0000_0000_0000,
            0xbfe0_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0x3ff0_0000_0000_0000,
            0x3ff4_0000_0000_0000,
        ];

        let mut audit = Audit::default();
        let mut ordinal = 0;
        for v1 in V1_BITS.map(f64::from_bits) {
            for v2 in V2_BITS.map(f64::from_bits) {
                for v3 in V3_BITS.map(f64::from_bits) {
                    for cos_phi in COS_PHI_BITS.map(f64::from_bits) {
                        audit.check_torsion_energy(
                            format!("torsion_energy_{ordinal:03}"),
                            v1,
                            v2,
                            v3,
                            cos_phi,
                        );
                        ordinal += 1;
                    }
                }
            }
        }

        if audit.cells != 48 {
            audit.record(format!("expected 48 source cells, got {}", audit.cells));
        }
        if audit.calls != 96 {
            audit.record(format!(
                "expected 96 actual Rust calls, got {}",
                audit.calls
            ));
        }
        assert!(
            audit.discrepancies.is_empty(),
            "MMFF U3 torsion energy discrepancies:\n{}",
            audit.discrepancies.join("\n")
        );
    }

    #[test]
    fn mmff_numerical_u4_oop_bend_energy_matches_source_bits() {
        const CHI_BITS: [u64; 6] = [
            0xc05e_0000_0000_0000,
            0xc03e_0000_0000_0000,
            0x8000_0000_0000_0000,
            0x0000_0000_0000_0000,
            0x402e_0000_0000_0000,
            0x4056_8000_0000_0000,
        ];
        const KOOP_BITS: [u64; 3] = [
            0x8000_0000_0000_0000,
            0x3fe0_0000_0000_0000,
            0x4000_0000_0000_0000,
        ];

        let mut audit = Audit::default();
        let mut ordinal = 0;
        for chi in CHI_BITS.map(f64::from_bits) {
            for koop in KOOP_BITS.map(f64::from_bits) {
                audit.check_oop_bend_energy(format!("oop_energy_{ordinal:02}"), chi, koop);
                ordinal += 1;
            }
        }

        if audit.cells != 18 {
            audit.record(format!("expected 18 source cells, got {}", audit.cells));
        }
        if audit.calls != 36 {
            audit.record(format!(
                "expected 36 actual Rust calls, got {}",
                audit.calls
            ));
        }
        assert!(
            audit.discrepancies.is_empty(),
            "MMFF U4 OOP bend energy discrepancies:\n{}",
            audit.discrepancies.join("\n")
        );
    }

    #[test]
    fn mmff_numerical_u5_cos_theta_matches_source_bits() {
        const Q_BITS: [[[u64; 3]; 4]; 8] = [
            [
                [0x3ff0_0000_0000_0000, 0, 0],
                [0, 0, 0],
                [0, 0x3ff0_0000_0000_0000, 0],
                [0, 0x3ff0_0000_0000_0000, 0x3ff0_0000_0000_0000],
            ],
            [
                [0x3ff0_0000_0000_0000, 0, 0],
                [0, 0, 0],
                [0xbff0_0000_0000_0000, 0, 0],
                [0xbff0_0000_0000_0000, 0x3ff0_0000_0000_0000, 0],
            ],
            [
                [0x3ff0_0000_0000_0000, 0, 0],
                [0, 0, 0],
                [0x3ff0_0000_0000_0000, 0, 0],
                [0x3ff0_0000_0000_0000, 0x3ff0_0000_0000_0000, 0],
            ],
            [
                [
                    0x3fd0_0000_0000_0000,
                    0xbfe0_0000_0000_0000,
                    0x3ff0_0000_0000_0000,
                ],
                [
                    0xbff4_0000_0000_0000,
                    0x3fe8_0000_0000_0000,
                    0x4000_0000_0000_0000,
                ],
                [
                    0x3fe0_0000_0000_0000,
                    0x3ff8_0000_0000_0000,
                    0xbfe8_0000_0000_0000,
                ],
                [
                    0x4000_0000_0000_0000,
                    0xbfd0_0000_0000_0000,
                    0x3fc0_0000_0000_0000,
                ],
            ],
            [[0, 0, 0]; 4],
            [
                [
                    0x4010_0000_0000_0000,
                    0xc000_0000_0000_0000,
                    0x3fe0_0000_0000_0000,
                ],
                [
                    0x4008_0000_0000_0000,
                    0xc000_0000_0000_0000,
                    0x3fe0_0000_0000_0000,
                ],
                [
                    0x4008_0000_0000_0000,
                    0xbff0_0000_0000_0000,
                    0x3fe0_0000_0000_0000,
                ],
                [
                    0x4008_0000_0000_0000,
                    0xbff0_0000_0000_0000,
                    0x3ff8_0000_0000_0000,
                ],
            ],
            [
                [
                    0x3fa0_0000_0000_0000,
                    0xbfb0_0000_0000_0000,
                    0x3fc0_0000_0000_0000,
                ],
                [
                    0xbfc4_0000_0000_0000,
                    0x3fb8_0000_0000_0000,
                    0x3fd0_0000_0000_0000,
                ],
                [
                    0x3fb0_0000_0000_0000,
                    0x3fc8_0000_0000_0000,
                    0xbfb8_0000_0000_0000,
                ],
                [
                    0x3fd0_0000_0000_0000,
                    0xbfa0_0000_0000_0000,
                    0x3f90_0000_0000_0000,
                ],
            ],
            [
                [
                    0x3ff0_0000_0000_0000,
                    0x8000_0000_0000_0000,
                    0x8000_0000_0000_0000,
                ],
                [0x8000_0000_0000_0000; 3],
                [
                    0x8000_0000_0000_0000,
                    0x3ff0_0000_0000_0000,
                    0x8000_0000_0000_0000,
                ],
                [
                    0x8000_0000_0000_0000,
                    0x3ff0_0000_0000_0000,
                    0x3ff0_0000_0000_0000,
                ],
            ],
        ];
        const DIST_BITS: [[u64; 2]; 3] = [
            [0x3ff0_0000_0000_0000, 0x3ff0_0000_0000_0000],
            [0x3fe0_0000_0000_0000, 0x4000_0000_0000_0000],
            [0x0000_0000_0000_0000, 0x3ff0_0000_0000_0000],
        ];

        let point = |bits: [u64; 3]| Point3 {
            x: f64::from_bits(bits[0]),
            y: f64::from_bits(bits[1]),
            z: f64::from_bits(bits[2]),
        };
        let mut audit = Audit::default();
        for (q, point_bits) in Q_BITS.into_iter().enumerate() {
            let points = point_bits.map(point);
            for (distance_index, [dist1_bits, dist2_bits]) in DIST_BITS.into_iter().enumerate() {
                audit.check_cos_theta(
                    format!("cos_theta_q{q}_d{distance_index}"),
                    points[0],
                    points[1],
                    points[2],
                    f64::from_bits(dist1_bits),
                    f64::from_bits(dist2_bits),
                );
            }
        }

        if audit.cells != 24 {
            audit.record(format!("expected 24 source cells, got {}", audit.cells));
        }
        if audit.calls != 48 {
            audit.record(format!(
                "expected 48 actual Rust calls, got {}",
                audit.calls
            ));
        }
        assert!(
            audit.discrepancies.is_empty(),
            "MMFF U5 cosine geometry discrepancies:\n{}",
            audit.discrepancies.join("\n")
        );
    }

    #[test]
    fn mmff_numerical_u6_torsion_cosine_matches_source_bits_and_preserves_points() {
        const Q_BITS: [[[u64; 3]; 4]; 10] = [
            [
                [0x3ff0_0000_0000_0000, 0, 0],
                [0, 0, 0],
                [0, 0x3ff0_0000_0000_0000, 0],
                [0, 0x3ff0_0000_0000_0000, 0x3ff0_0000_0000_0000],
            ],
            [
                [0x3ff0_0000_0000_0000, 0, 0],
                [0, 0, 0],
                [0xbff0_0000_0000_0000, 0, 0],
                [0xbff0_0000_0000_0000, 0x3ff0_0000_0000_0000, 0],
            ],
            [
                [0x3ff0_0000_0000_0000, 0, 0],
                [0, 0, 0],
                [0x3ff0_0000_0000_0000, 0, 0],
                [0x3ff0_0000_0000_0000, 0x3ff0_0000_0000_0000, 0],
            ],
            [
                [
                    0x3fd0_0000_0000_0000,
                    0xbfe0_0000_0000_0000,
                    0x3ff0_0000_0000_0000,
                ],
                [
                    0xbff4_0000_0000_0000,
                    0x3fe8_0000_0000_0000,
                    0x4000_0000_0000_0000,
                ],
                [
                    0x3fe0_0000_0000_0000,
                    0x3ff8_0000_0000_0000,
                    0xbfe8_0000_0000_0000,
                ],
                [
                    0x4000_0000_0000_0000,
                    0xbfd0_0000_0000_0000,
                    0x3fc0_0000_0000_0000,
                ],
            ],
            [[0, 0, 0]; 4],
            [
                [
                    0x4010_0000_0000_0000,
                    0xc000_0000_0000_0000,
                    0x3fe0_0000_0000_0000,
                ],
                [
                    0x4008_0000_0000_0000,
                    0xc000_0000_0000_0000,
                    0x3fe0_0000_0000_0000,
                ],
                [
                    0x4008_0000_0000_0000,
                    0xbff0_0000_0000_0000,
                    0x3fe0_0000_0000_0000,
                ],
                [
                    0x4008_0000_0000_0000,
                    0xbff0_0000_0000_0000,
                    0x3ff8_0000_0000_0000,
                ],
            ],
            [
                [
                    0x3fa0_0000_0000_0000,
                    0xbfb0_0000_0000_0000,
                    0x3fc0_0000_0000_0000,
                ],
                [
                    0xbfc4_0000_0000_0000,
                    0x3fb8_0000_0000_0000,
                    0x3fd0_0000_0000_0000,
                ],
                [
                    0x3fb0_0000_0000_0000,
                    0x3fc8_0000_0000_0000,
                    0xbfb8_0000_0000_0000,
                ],
                [
                    0x3fd0_0000_0000_0000,
                    0xbfa0_0000_0000_0000,
                    0x3f90_0000_0000_0000,
                ],
            ],
            [
                [
                    0x3ff0_0000_0000_0000,
                    0x8000_0000_0000_0000,
                    0x8000_0000_0000_0000,
                ],
                [0x8000_0000_0000_0000; 3],
                [
                    0x8000_0000_0000_0000,
                    0x3ff0_0000_0000_0000,
                    0x8000_0000_0000_0000,
                ],
                [
                    0x8000_0000_0000_0000,
                    0x3ff0_0000_0000_0000,
                    0x3ff0_0000_0000_0000,
                ],
            ],
            [
                [0x3eb0_c6f7_a0b5_ed8d, 0, 0],
                [0, 0, 0],
                [0, 0x3eb0_c6f7_a0b5_ed8d, 0],
                [0, 0x3eb0_c6f7_a0b5_ed8d, 0x3eb0_c6f7_a0b5_ed8d],
            ],
            [
                [0x3f1a_36e2_eb1c_432d, 0, 0],
                [0, 0, 0],
                [0, 0x3f1a_36e2_eb1c_432d, 0],
                [0, 0x3f1a_36e2_eb1c_432d, 0x3f1a_36e2_eb1c_432d],
            ],
        ];

        let point = |bits: [u64; 3]| Point3 {
            x: f64::from_bits(bits[0]),
            y: f64::from_bits(bits[1]),
            z: f64::from_bits(bits[2]),
        };
        let mut audit = Audit::default();
        for (q, point_bits) in Q_BITS.into_iter().enumerate() {
            let points = point_bits.map(point);
            audit.check_torsion_cos_phi(
                format!("cos_phi_q{q}_order0"),
                [&points[0], &points[1], &points[2], &points[3]],
            );
            let reversed = [&points[3], &points[2], &points[1], &points[0]];
            audit.check_torsion_cos_phi(format!("cos_phi_q{q}_order1"), reversed);
        }

        if audit.cells != 20 {
            audit.record(format!("expected 20 source cells, got {}", audit.cells));
        }
        if audit.calls != 40 {
            audit.record(format!(
                "expected 40 actual Rust calls, got {}",
                audit.calls
            ));
        }
        assert!(
            audit.discrepancies.is_empty(),
            "MMFF U6 torsion cosine discrepancies:\n{}",
            audit.discrepancies.join("\n")
        );
    }

    #[test]
    fn mmff_numerical_u7_oop_chi_matches_source_bits_and_preserves_points() {
        const Q_BITS: [[[u64; 3]; 4]; 8] = [
            [
                [0x3ff0_0000_0000_0000, 0, 0],
                [0, 0, 0],
                [0, 0x3ff0_0000_0000_0000, 0],
                [0, 0x3ff0_0000_0000_0000, 0x3ff0_0000_0000_0000],
            ],
            [
                [0x3ff0_0000_0000_0000, 0, 0],
                [0, 0, 0],
                [0xbff0_0000_0000_0000, 0, 0],
                [0xbff0_0000_0000_0000, 0x3ff0_0000_0000_0000, 0],
            ],
            [
                [0x3ff0_0000_0000_0000, 0, 0],
                [0, 0, 0],
                [0x3ff0_0000_0000_0000, 0, 0],
                [0x3ff0_0000_0000_0000, 0x3ff0_0000_0000_0000, 0],
            ],
            [
                [
                    0x3fd0_0000_0000_0000,
                    0xbfe0_0000_0000_0000,
                    0x3ff0_0000_0000_0000,
                ],
                [
                    0xbff4_0000_0000_0000,
                    0x3fe8_0000_0000_0000,
                    0x4000_0000_0000_0000,
                ],
                [
                    0x3fe0_0000_0000_0000,
                    0x3ff8_0000_0000_0000,
                    0xbfe8_0000_0000_0000,
                ],
                [
                    0x4000_0000_0000_0000,
                    0xbfd0_0000_0000_0000,
                    0x3fc0_0000_0000_0000,
                ],
            ],
            [[0, 0, 0]; 4],
            [
                [
                    0x4010_0000_0000_0000,
                    0xc000_0000_0000_0000,
                    0x3fe0_0000_0000_0000,
                ],
                [
                    0x4008_0000_0000_0000,
                    0xc000_0000_0000_0000,
                    0x3fe0_0000_0000_0000,
                ],
                [
                    0x4008_0000_0000_0000,
                    0xbff0_0000_0000_0000,
                    0x3fe0_0000_0000_0000,
                ],
                [
                    0x4008_0000_0000_0000,
                    0xbff0_0000_0000_0000,
                    0x3ff8_0000_0000_0000,
                ],
            ],
            [
                [
                    0x3fa0_0000_0000_0000,
                    0xbfb0_0000_0000_0000,
                    0x3fc0_0000_0000_0000,
                ],
                [
                    0xbfc4_0000_0000_0000,
                    0x3fb8_0000_0000_0000,
                    0x3fd0_0000_0000_0000,
                ],
                [
                    0x3fb0_0000_0000_0000,
                    0x3fc8_0000_0000_0000,
                    0xbfb8_0000_0000_0000,
                ],
                [
                    0x3fd0_0000_0000_0000,
                    0xbfa0_0000_0000_0000,
                    0x3f90_0000_0000_0000,
                ],
            ],
            [
                [
                    0x3ff0_0000_0000_0000,
                    0x8000_0000_0000_0000,
                    0x8000_0000_0000_0000,
                ],
                [0x8000_0000_0000_0000; 3],
                [
                    0x8000_0000_0000_0000,
                    0x3ff0_0000_0000_0000,
                    0x8000_0000_0000_0000,
                ],
                [
                    0x8000_0000_0000_0000,
                    0x3ff0_0000_0000_0000,
                    0x3ff0_0000_0000_0000,
                ],
            ],
        ];

        let point = |bits: [u64; 3]| Point3 {
            x: f64::from_bits(bits[0]),
            y: f64::from_bits(bits[1]),
            z: f64::from_bits(bits[2]),
        };
        let mut audit = Audit::default();
        for (q, point_bits) in Q_BITS.into_iter().enumerate() {
            let points = point_bits.map(point);
            audit.check_oop_chi(
                format!("oop_chi_q{q}_order0"),
                [&points[0], &points[1], &points[2], &points[3]],
            );
            let reversed = [&points[3], &points[2], &points[1], &points[0]];
            audit.check_oop_chi(format!("oop_chi_q{q}_order1"), reversed);
        }

        if audit.cells != 16 {
            audit.record(format!("expected 16 source cells, got {}", audit.cells));
        }
        if audit.calls != 32 {
            audit.record(format!(
                "expected 32 actual Rust calls, got {}",
                audit.calls
            ));
        }
        assert!(
            audit.discrepancies.is_empty(),
            "MMFF U7 OOP chi discrepancies:\n{}",
            audit.discrepancies.join("\n")
        );
    }

    #[test]
    fn mmff_numerical_u8_angle_gradient_matches_source_bits_and_preserves_aliases() {
        const R_BITS: [[[u64; 3]; 2]; 2] = [
            [[0x3ff0_0000_0000_0000, 0, 0], [0, 0, 0]],
            [
                [
                    0x3fd0_0000_0000_0000,
                    0xbfe0_0000_0000_0000,
                    0x3ff0_0000_0000_0000,
                ],
                [
                    0xbff4_0000_0000_0000,
                    0x3fe8_0000_0000_0000,
                    0x4000_0000_0000_0000,
                ],
            ],
        ];
        const DIST_BITS: [[u64; 2]; 2] = [
            [0x3ff0_0000_0000_0000, 0x4000_0000_0000_0000],
            [0x3fe0_0000_0000_0000, 0x3ff4_0000_0000_0000],
        ];
        const D_E_BITS: [u64; 2] = [0x8000_0000_0000_0000, 0x3ff4_0000_0000_0000];
        const COS_SIN_BITS: [[u64; 2]; 2] = [
            [0xbfe0_0000_0000_0000, 0xbfe9_9999_9999_999a],
            [0x3fe8_0000_0000_0000, 0x3fe3_3333_3333_3333],
        ];
        const INDICES: [[usize; 3]; 3] = [[1, 3, 5], [2, 2, 2], [1, 3, 1]];

        let point = |bits: [u64; 3]| Point3 {
            x: f64::from_bits(bits[0]),
            y: f64::from_bits(bits[1]),
            z: f64::from_bits(bits[2]),
        };
        let r_values = R_BITS.map(|rows| rows.map(point));
        let mut audit = Audit::default();
        let mut ordinal = 0;
        for r in &r_values {
            for dist_bits in DIST_BITS {
                let dist = dist_bits.map(f64::from_bits);
                for d_e_bits in D_E_BITS {
                    for [cos_bits, sin_bits] in COS_SIN_BITS {
                        for indices in INDICES {
                            for baseline_bits in GRADIENT_BASELINES_BITS {
                                audit.check_angle_bend_gradient(
                                    format!("angle_grad_{ordinal:03}"),
                                    r,
                                    &dist,
                                    f64::from_bits(d_e_bits),
                                    f64::from_bits(cos_bits),
                                    f64::from_bits(sin_bits),
                                    indices,
                                    baseline_bits,
                                );
                                ordinal += 1;
                            }
                        }
                    }
                }
            }
        }

        audit.check_angle_bend_error(
            "u8_empty_rows_first_zero",
            [0, 0, 0],
            MmffGradientError::GradientRowOutOfRange {
                slot: 0,
                index: 0,
                rows: 0,
            },
            true,
        );
        audit.check_angle_bend_error(
            "u8_first_index_seven",
            [7, 0, 0],
            MmffGradientError::GradientRowOutOfRange {
                slot: 0,
                index: 7,
                rows: 7,
            },
            false,
        );
        audit.check_angle_bend_error(
            "u8_last_index_seven",
            [0, 0, 7],
            MmffGradientError::GradientRowOutOfRange {
                slot: 2,
                index: 7,
                rows: 7,
            },
            false,
        );
        audit.check_angle_bend_error(
            "u8_first_index_usize_max",
            [usize::MAX, 0, 0],
            MmffGradientError::GradientRowOutOfRange {
                slot: 0,
                index: usize::MAX,
                rows: 7,
            },
            false,
        );

        if audit.cells != 96 {
            audit.record(format!("expected 96 source cells, got {}", audit.cells));
        }
        if audit.calls != 192 {
            audit.record(format!(
                "expected 192 actual valid Rust calls, got {}",
                audit.calls
            ));
        }
        if audit.error_controls != 4 {
            audit.record(format!(
                "expected 4 CK-only bounds controls, got {}",
                audit.error_controls
            ));
        }
        assert!(
            audit.discrepancies.is_empty(),
            "MMFF U8 angle gradient discrepancies:\n{}",
            audit.discrepancies.join("\n")
        );
    }

    #[test]
    fn mmff_numerical_u9_torsion_gradient_matches_source_bits_and_preserves_aliases() {
        const R_BITS: [[[u64; 3]; 4]; 2] = [
            [
                [0x3ff0_0000_0000_0000, 0, 0],
                [0, 0, 0],
                [0, 0x3ff0_0000_0000_0000, 0],
                [0, 0x3ff0_0000_0000_0000, 0x3ff0_0000_0000_0000],
            ],
            [
                [
                    0x3fd0_0000_0000_0000,
                    0xbfe0_0000_0000_0000,
                    0x3ff0_0000_0000_0000,
                ],
                [
                    0xbff4_0000_0000_0000,
                    0x3fe8_0000_0000_0000,
                    0x4000_0000_0000_0000,
                ],
                [
                    0x3fe0_0000_0000_0000,
                    0x3ff8_0000_0000_0000,
                    0xbfe8_0000_0000_0000,
                ],
                [
                    0x4000_0000_0000_0000,
                    0xbfd0_0000_0000_0000,
                    0x3fc0_0000_0000_0000,
                ],
            ],
        ];
        const T_BITS: [[[u64; 3]; 2]; 2] = [
            [[0x3ff0_0000_0000_0000, 0, 0], [0, 0x3ff0_0000_0000_0000, 0]],
            [
                [
                    0x3fd0_0000_0000_0000,
                    0xbfe0_0000_0000_0000,
                    0x3ff0_0000_0000_0000,
                ],
                [
                    0xbff4_0000_0000_0000,
                    0x3fe8_0000_0000_0000,
                    0x4000_0000_0000_0000,
                ],
            ],
        ];
        const D_BITS: [[u64; 2]; 2] = [
            [0x3ff0_0000_0000_0000, 0x4000_0000_0000_0000],
            [0x3fe0_0000_0000_0000, 0x3ff4_0000_0000_0000],
        ];
        const SIN_BITS: [u64; 2] = [0x8000_0000_0000_0000, 0x3ff4_0000_0000_0000];
        const COS_BITS: [u64; 2] = [0xbfe0_0000_0000_0000, 0x3fe8_0000_0000_0000];
        const INDICES: [[usize; 4]; 3] = [[1, 3, 5, 6], [2, 2, 2, 2], [1, 3, 1, 3]];

        let point = |bits: [u64; 3]| Point3 {
            x: f64::from_bits(bits[0]),
            y: f64::from_bits(bits[1]),
            z: f64::from_bits(bits[2]),
        };
        let r_values = R_BITS.map(|rows| rows.map(point));
        let t_values = T_BITS.map(|rows| rows.map(point));
        let mut audit = Audit::default();
        let mut ordinal = 0;
        for r in &r_values {
            for t in &t_values {
                for d_bits in D_BITS {
                    let d = d_bits.map(f64::from_bits);
                    for sin_bits in SIN_BITS {
                        for cos_bits in COS_BITS {
                            for indices in INDICES {
                                for baseline_bits in GRADIENT_BASELINES_BITS {
                                    audit.check_torsion_gradient(
                                        format!("torsion_grad_{ordinal:03}"),
                                        r,
                                        t,
                                        &d,
                                        f64::from_bits(sin_bits),
                                        f64::from_bits(cos_bits),
                                        indices,
                                        baseline_bits,
                                    );
                                    ordinal += 1;
                                }
                            }
                        }
                    }
                }
            }
        }

        audit.check_torsion_gradient_error(
            "u9_empty_rows_first_zero",
            [0, 0, 0, 0],
            MmffGradientError::GradientRowOutOfRange {
                slot: 0,
                index: 0,
                rows: 0,
            },
            true,
        );
        audit.check_torsion_gradient_error(
            "u9_first_index_seven",
            [7, 0, 0, 0],
            MmffGradientError::GradientRowOutOfRange {
                slot: 0,
                index: 7,
                rows: 7,
            },
            false,
        );
        audit.check_torsion_gradient_error(
            "u9_last_index_seven",
            [0, 0, 0, 7],
            MmffGradientError::GradientRowOutOfRange {
                slot: 3,
                index: 7,
                rows: 7,
            },
            false,
        );
        audit.check_torsion_gradient_error(
            "u9_first_index_usize_max",
            [usize::MAX, 0, 0, 0],
            MmffGradientError::GradientRowOutOfRange {
                slot: 0,
                index: usize::MAX,
                rows: 7,
            },
            false,
        );

        if audit.cells != 192 {
            audit.record(format!("expected 192 source cells, got {}", audit.cells));
        }
        if audit.calls != 384 {
            audit.record(format!(
                "expected 384 actual valid Rust calls, got {}",
                audit.calls
            ));
        }
        if audit.error_controls != 4 {
            audit.record(format!(
                "expected 4 CK-only bounds controls, got {}",
                audit.error_controls
            ));
        }
        assert!(
            audit.discrepancies.is_empty(),
            "MMFF U9 torsion gradient discrepancies:\n{}",
            audit.discrepancies.join("\n")
        );
    }
}
