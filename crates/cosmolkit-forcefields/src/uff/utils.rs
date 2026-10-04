// Copyright (C) 2013-2024 Paolo Tosco and other RDKit contributors.
// @@ All Rights Reserved @@
// This file is part of the RDKit.
// The contents are covered by the terms of the BSD license
// which is included in the file license.txt, found at the root
// of the pinned RDKit source tree.

use crate::geometry::Point3;

pub(super) fn calculate_cos_y(
    i_point: &Point3,
    j_point: &Point3,
    k_point: &Point3,
    l_point: &Point3,
) -> f64 {
    // RDKit✔️✔️: double calculateCosY(const RDGeom::Point3D &iPoint,
    // RDKit✔️✔️:                      const RDGeom::Point3D &jPoint,
    // RDKit✔️✔️:                      const RDGeom::Point3D &kPoint,
    // RDKit✔️✔️:                      const RDGeom::Point3D &lPoint) {
    // RDKit✔️✔️:   constexpr double zeroTol = 1.0e-16;
    // RDKit✔️✔️:   RDGeom::Point3D rJI = iPoint - jPoint;
    // RDKit✔️✔️:   RDGeom::Point3D rJK = kPoint - jPoint;
    // RDKit✔️✔️:   RDGeom::Point3D rJL = lPoint - jPoint;
    // RDKit✔️✔️:   auto l2JI = rJI.lengthSq();
    // RDKit✔️✔️:   auto l2JK = rJK.lengthSq();
    // RDKit✔️✔️:   auto l2JL = rJL.lengthSq();
    // RDKit✔️✔️:   if (l2JI < zeroTol || l2JK < zeroTol || l2JL < zeroTol) {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   RDGeom::Point3D n = rJI.crossProduct(rJK);
    // RDKit✔️✔️:   n /= (sqrt(l2JI) * sqrt(l2JK));
    // RDKit✔️✔️:   auto l2n = n.lengthSq();
    // RDKit✔️✔️:   if (l2n < zeroTol) {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return n.dotProduct(rJL) / (sqrt(l2JL) * sqrt(l2n));
    // RDKit✔️✔️: }
    const ZERO_TOL: f64 = 1.0e-16;

    // BEGIN RDKIT CPP HELPER RDGeom::operator- (Geometry/point.cpp:64-70)
    // RDKit✔️✔️: Point3D operator-(const Point3D &p1, const Point3D &p2) {
    // RDKit✔️✔️:   Point3D res;
    // RDKit✔️✔️:   res.x = p1.x - p2.x;
    // RDKit✔️✔️:   res.y = p1.y - p2.y;
    // RDKit✔️✔️:   res.z = p1.z - p2.z;
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP HELPER RDGeom::operator-
    let r_ji = Point3::difference(i_point, j_point);
    let r_jk = Point3::difference(k_point, j_point);
    let r_jl = Point3::difference(l_point, j_point);

    // BEGIN RDKIT CPP HELPER RDGeom::Point3D::lengthSq (Geometry/point.h:163-168)
    // RDKit✔️✔️: constexpr double lengthSq() const override {
    // RDKit✔️✔️:   // double res = pow(x,2) + pow(y,2) + pow(z,2);
    // RDKit✔️✔️:   double res = x * x + y * y + z * z;
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP HELPER RDGeom::Point3D::lengthSq
    let l2_ji = r_ji.length_sq();
    let l2_jk = r_jk.length_sq();
    let l2_jl = r_jl.length_sq();
    if l2_ji < ZERO_TOL || l2_jk < ZERO_TOL || l2_jl < ZERO_TOL {
        return 0.0;
    }

    // BEGIN RDKIT CPP HELPER RDGeom::Point3D::crossProduct (Geometry/point.h:228-234)
    // RDKit✔️✔️: constexpr Point3D crossProduct(const Point3D &other) const {
    // RDKit✔️✔️:   Point3D res;
    // RDKit✔️✔️:   res.x = y * (other.z) - z * (other.y);
    // RDKit✔️✔️:   res.y = -x * (other.z) + z * (other.x);
    // RDKit✔️✔️:   res.z = x * (other.y) - y * (other.x);
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP HELPER RDGeom::Point3D::crossProduct
    let mut n = r_ji.cross_product(&r_jk);

    // BEGIN RDKIT CPP HELPER RDGeom::Point3D::operator/= (Geometry/point.h:132-137)
    // RDKit✔️✔️: constexpr Point3D &operator/=(double scale) {
    // RDKit✔️✔️:   x /= scale;
    // RDKit✔️✔️:   y /= scale;
    // RDKit✔️✔️:   z /= scale;
    // RDKit✔️✔️:   return *this;
    // RDKit✔️✔️: }
    // END RDKIT CPP HELPER RDGeom::Point3D::operator/=
    let normalizer = l2_ji.sqrt() * l2_jk.sqrt();
    n.x /= normalizer;
    n.y /= normalizer;
    n.z /= normalizer;

    let l2_n = n.length_sq();
    if l2_n < ZERO_TOL {
        return 0.0;
    }

    // BEGIN RDKIT CPP HELPER RDGeom::Point3D::dotProduct (Geometry/point.h:169-172)
    // RDKit✔️✔️: constexpr double dotProduct(const Point3D &other) const {
    // RDKit✔️✔️:   double res = x * (other.x) + y * (other.y) + z * (other.z);
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP HELPER RDGeom::Point3D::dotProduct
    let numerator = n.x * r_jl.x + n.y * r_jl.y + n.z * r_jl.z;
    numerator / (l2_jl.sqrt() * l2_n.sqrt())
}

pub(super) fn calc_inversion_coefficients(
    at2_atomic_num: i32,
    is_c_bound_to_o: bool,
) -> (f64, f64, f64, f64) {
    // RDKit✔️✔️: std::tuple<double, double, double, double>
    // RDKit✔️✔️: calcInversionCoefficientsAndForceConstant(int at2AtomicNum, bool isCBoundToO) {
    // RDKit✔️✔️:   double res = 0.0;
    // RDKit✔️✔️:   double C0 = 0.0;
    // RDKit✔️✔️:   double C1 = 0.0;
    // RDKit✔️✔️:   double C2 = 0.0;
    let mut res = 0.0;
    let c0;
    let c1;
    let c2;

    // RDKit✔️✔️:   // if the central atom is sp2 carbon, nitrogen or oxygen
    // RDKit✔️✔️:   if ((at2AtomicNum == 6) || (at2AtomicNum == 7) || (at2AtomicNum == 8)) {
    if at2_atomic_num == 6 || at2_atomic_num == 7 || at2_atomic_num == 8 {
        // RDKit✔️✔️:     C0 = 1.0;
        // RDKit✔️✔️:     C1 = -1.0;
        // RDKit✔️✔️:     C2 = 0.0;
        c0 = 1.0;
        c1 = -1.0;
        c2 = 0.0;
        // RDKit✔️✔️:     res = (isCBoundToO ? 50.0 : 6.0);
        res = if is_c_bound_to_o { 50.0 } else { 6.0 };
        // RDKit✔️✔️:   } else {
    } else {
        // RDKit✔️✔️:     // group 5 elements are not clearly explained in the UFF paper
        // RDKit✔️✔️:     // the following code was inspired by MCCCS Towhee's ffuff.F
        // RDKit✔️✔️:     double w0 = M_PI / 180.0;
        let mut w0 = std::f64::consts::PI / 180.0;
        // RDKit✔️✔️:     switch (at2AtomicNum) {
        match at2_atomic_num {
            // RDKit✔️✔️:       // if the central atom is phosphorous
            // RDKit✔️✔️:       case 15:
            // RDKit✔️✔️:         w0 *= 84.4339;
            // RDKit✔️✔️:         break;
            15 => w0 *= 84.4339,
            // RDKit✔️✔️:
            // RDKit✔️✔️:       // if the central atom is arsenic
            // RDKit✔️✔️:       case 33:
            // RDKit✔️✔️:         w0 *= 86.9735;
            // RDKit✔️✔️:         break;
            33 => w0 *= 86.9735,
            // RDKit✔️✔️:
            // RDKit✔️✔️:       // if the central atom is antimonium
            // RDKit✔️✔️:       case 51:
            // RDKit✔️✔️:         w0 *= 87.7047;
            // RDKit✔️✔️:         break;
            51 => w0 *= 87.7047,
            // RDKit✔️✔️:
            // RDKit✔️✔️:       // if the central atom is bismuth
            // RDKit✔️✔️:       case 83:
            // RDKit✔️✔️:         w0 *= 90.0;
            // RDKit✔️✔️:         break;
            83 => w0 *= 90.0,
            // RDKit✔️✔️:     }
            _ => {}
        }
        // RDKit✔️✔️:     C2 = 1.0;
        c2 = 1.0;
        // RDKit✔️✔️:     C1 = -4.0 * cos(w0);
        c1 = -4.0 * w0.cos();
        // RDKit✔️✔️:     C0 = -(C1 * cos(w0) + C2 * cos(2.0 * w0));
        c0 = -(c1 * w0.cos() + c2 * (2.0 * w0).cos());
        // RDKit✔️✔️:     res = 22.0 / (C0 + C1 + C2);
        res = 22.0 / (c0 + c1 + c2);
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   res /= 3.0;
    res /= 3.0;

    // RDKit✔️✔️:   return std::make_tuple(res, C0, C1, C2);
    // RDKit✔️✔️: }
    (res, c0, c1, c2)
}

#[cfg(test)]
mod tests {
    use super::{calc_inversion_coefficients, calculate_cos_y};
    use crate::geometry::Point3;

    const ZERO_TOL: f64 = 1.0e-16;
    const INVERSION_TOL: f64 = 1.0e-12;

    fn point(x: f64, y: f64, z: f64) -> Point3 {
        Point3 { x, y, z }
    }

    fn assert_same_bits(actual: f64, expected: f64) {
        assert_eq!(actual.to_bits(), expected.to_bits());
    }

    fn assert_tuple_close(actual: (f64, f64, f64, f64), expected: (f64, f64, f64, f64)) {
        for (component, (actual, expected)) in [
            (actual.0, expected.0),
            (actual.1, expected.1),
            (actual.2, expected.2),
            (actual.3, expected.3),
        ]
        .into_iter()
        .enumerate()
        {
            let tolerance = INVERSION_TOL * expected.abs().max(1.0);
            assert!(
                (actual - expected).abs() <= tolerance,
                "inversion component {component}: actual {actual:?}, expected {expected:?}"
            );
        }
    }

    #[test]
    fn cf3d_u02_planar_out_of_plane_and_handedness() {
        // RDKit ForceField/UFF/Utils.cpp::calculateCosY uses the ordered
        // rJI x rJK plane normal and projects rJL onto that signed normal.
        let origin = point(0.0, 0.0, 0.0);
        let i = point(1.0, 0.0, 0.0);
        let k = point(0.0, 1.0, 0.0);
        let planar_l = point(1.0, 1.0, 0.0);
        let out_of_plane_l = point(1.0, 1.0, 1.0);

        assert_same_bits(calculate_cos_y(&i, &origin, &k, &planar_l), 0.0);
        assert_same_bits(
            calculate_cos_y(&i, &origin, &k, &out_of_plane_l),
            1.0 / 3.0_f64.sqrt(),
        );

        let reverse_l = point(1.0, 1.0, -1.0);
        assert_same_bits(
            calculate_cos_y(&i, &origin, &k, &reverse_l),
            -1.0 / 3.0_f64.sqrt(),
        );
        assert_same_bits(
            calculate_cos_y(&k, &origin, &i, &out_of_plane_l),
            -1.0 / 3.0_f64.sqrt(),
        );
    }

    #[test]
    fn cf3d_u02_zero_and_short_vectors_follow_each_source_guard_operand() {
        let origin = point(0.0, 0.0, 0.0);
        let i = point(1.0, 0.0, 0.0);
        let k = point(0.0, 1.0, 0.0);
        let l = point(0.0, 0.0, 1.0);

        // The first source guard tests l2JI, then l2JK, then l2JL.
        assert_same_bits(calculate_cos_y(&origin, &origin, &k, &l), 0.0);
        assert_same_bits(calculate_cos_y(&i, &origin, &origin, &l), 0.0);
        assert_same_bits(calculate_cos_y(&i, &origin, &k, &origin), 0.0);

        let short = point(3.0e-9, 4.0e-9, 0.0);
        assert!(short.x * short.x + short.y * short.y < ZERO_TOL);
        assert_same_bits(calculate_cos_y(&short, &origin, &k, &l), 0.0);
        assert_same_bits(calculate_cos_y(&i, &origin, &short, &l), 0.0);

        let short_l = point(0.0, 3.0e-9, 4.0e-9);
        assert!(short_l.y * short_l.y + short_l.z * short_l.z < ZERO_TOL);
        assert_same_bits(calculate_cos_y(&i, &origin, &k, &short_l), 0.0);

        // These 6:8 vectors have squared length exactly zeroTol in f64. The
        // source uses `<`, so equality proceeds through normalization.
        let cutoff_ji = point(6.0e-9, 8.0e-9, 0.0);
        assert_eq!(
            cutoff_ji.x * cutoff_ji.x + cutoff_ji.y * cutoff_ji.y,
            ZERO_TOL
        );
        let aligned_normal = point(0.8, -0.6, 0.0);
        assert_same_bits(
            calculate_cos_y(&cutoff_ji, &origin, &l, &aligned_normal),
            1.0,
        );

        let cutoff_jk = point(6.0e-9, 8.0e-9, 0.0);
        assert_eq!(
            cutoff_jk.x * cutoff_jk.x + cutoff_jk.y * cutoff_jk.y,
            ZERO_TOL
        );
        let cutoff_i = point(0.0, 0.0, 1.0);
        let other_aligned_normal = point(-0.8, 0.6, 0.0);
        assert_same_bits(
            calculate_cos_y(&cutoff_i, &origin, &cutoff_jk, &other_aligned_normal),
            1.0,
        );

        let cutoff_jl = point(0.0, 6.0e-9, 8.0e-9);
        assert_eq!(
            cutoff_jl.y * cutoff_jl.y + cutoff_jl.z * cutoff_jl.z,
            ZERO_TOL
        );
        assert_same_bits(calculate_cos_y(&i, &origin, &k, &cutoff_jl), 0.8);
    }

    #[test]
    fn cf3d_u02_collinear_guard_uses_strict_normalized_cross_threshold() {
        let origin = point(0.0, 0.0, 0.0);
        let i = point(1.0, 0.0, 0.0);
        let l = point(0.0, 0.0, 1.0);

        let exactly_collinear = point(2.0, 0.0, 0.0);
        assert_same_bits(calculate_cos_y(&i, &origin, &exactly_collinear, &l), 0.0);

        let below_cutoff = point(1.0, 3.0e-9, 4.0e-9);
        assert_same_bits(calculate_cos_y(&i, &origin, &below_cutoff, &l), 0.0);

        // The input length square rounds to one, leaving l2n exactly equal to
        // zeroTol; the source's strict `<` must therefore compute a result.
        let at_cutoff = point(1.0, 6.0e-9, 8.0e-9);
        let unit_normal = point(0.0, -0.8, 0.6);
        let actual = calculate_cos_y(&i, &origin, &at_cutoff, &unit_normal);
        assert_same_bits(actual, 1.0);
    }

    #[test]
    fn cf3d_u02_nan_flows_through_source_comparisons() {
        let nan_i = point(f64::NAN, 0.0, 0.0);
        let origin = point(0.0, 0.0, 0.0);
        let k = point(0.0, 1.0, 0.0);
        let l = point(0.0, 0.0, 1.0);

        assert!(calculate_cos_y(&nan_i, &origin, &k, &l).is_nan());
    }

    #[test]
    fn cf3d_u03_sp2_elements_preserve_the_unfiltered_oxygen_flag_branch() {
        for atomic_num in [6, 7, 8] {
            assert_eq!(
                calc_inversion_coefficients(atomic_num, false),
                (2.0, 1.0, -1.0, 0.0)
            );
            assert_eq!(
                calc_inversion_coefficients(atomic_num, true),
                (50.0 / 3.0, 1.0, -1.0, 0.0)
            );
        }
    }

    #[test]
    fn cf3d_u03_group5_switch_multipliers_match_fixed_source_values() {
        let cases = [
            (
                15,
                (
                    4.496661509237746,
                    1.018815687546345,
                    -0.3879761595391661,
                    1.0,
                ),
            ),
            (
                33,
                (
                    4.08682518341991,
                    1.0055752214991891,
                    -0.21119131609399103,
                    1.0,
                ),
            ),
            (
                51,
                (
                    3.979001005763598,
                    1.0032079774467704,
                    -0.16019931202774357,
                    1.0,
                ),
            ),
            (83, (3.6666666666666674, 1.0, -2.4492935982947064e-16, 1.0)),
        ];

        for (atomic_num, expected) in cases {
            let oxygen_unbound = calc_inversion_coefficients(atomic_num, false);
            let oxygen_bound = calc_inversion_coefficients(atomic_num, true);
            assert_tuple_close(oxygen_unbound, expected);
            assert_tuple_close(oxygen_bound, expected);
            assert_eq!(oxygen_unbound, oxygen_bound);
        }
    }

    #[test]
    fn cf3d_u03_unrecognized_atomic_numbers_fall_through_at_one_degree() {
        let expected = (
            158068015.25333518,
            2.999390827019096,
            -3.999390780625565,
            1.0,
        );
        for atomic_num in [
            i32::MIN,
            -1,
            0,
            1,
            5,
            9,
            14,
            16,
            32,
            34,
            50,
            52,
            82,
            84,
            118,
            i32::MAX,
        ] {
            let oxygen_unbound = calc_inversion_coefficients(atomic_num, false);
            let oxygen_bound = calc_inversion_coefficients(atomic_num, true);
            assert_tuple_close(oxygen_unbound, expected);
            assert_tuple_close(oxygen_bound, expected);
            assert_eq!(oxygen_unbound, oxygen_bound);
        }
    }
}
