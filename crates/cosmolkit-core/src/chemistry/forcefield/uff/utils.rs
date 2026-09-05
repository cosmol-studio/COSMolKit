//! Source-backed RDKit UFF helper calculations.

use std::fmt;

use super::params::{AtomicParams, PARAMS_G, PARAMS_LAMBDA};

#[derive(Debug, Clone, PartialEq)]
pub enum UffUtilsError {
    InvalidBondOrder { bond_order: f64 },
    InvalidRestLength { rest_length: f64 },
}

impl fmt::Display for UffUtilsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidBondOrder { bond_order } => {
                write!(f, "UFF bond order must be positive, got {bond_order}")
            }
            Self::InvalidRestLength { rest_length } => {
                write!(f, "UFF rest length must be positive, got {rest_length}")
            }
        }
    }
}

impl std::error::Error for UffUtilsError {}

pub fn calc_bond_rest_length(
    bond_order: f64,
    end1_params: &AtomicParams,
    end2_params: &AtomicParams,
) -> Result<f64, UffUtilsError> {
    // RDKit✔️✔️: double calcBondRestLength(double bondOrder, const AtomicParams *end1Params,
    // RDKit✔️✔️:                           const AtomicParams *end2Params) {
    // RDKit✔️✔️:   PRECONDITION(bondOrder > 0, "bad bond order");
    if bond_order <= 0.0 {
        return Err(UffUtilsError::InvalidBondOrder { bond_order });
    }

    // RDKit✔️✔️:   double ri = end1Params->r1, rj = end2Params->r1;
    let ri = end1_params.r1;
    let rj = end2_params.r1;

    // RDKit✔️✔️:   // this is the pauling correction:
    // RDKit✔️✔️:   double rBO = -Params::lambda * (ri + rj) * log(bondOrder);
    let r_bo = -PARAMS_LAMBDA * (ri + rj) * bond_order.ln();

    // RDKit✔️✔️:   // O'Keefe and Breese electronegativity correction:
    // RDKit✔️✔️:   double Xi = end1Params->GMP_Xi, Xj = end2Params->GMP_Xi;
    let xi = end1_params.gmp_xi;
    let xj = end2_params.gmp_xi;
    // RDKit✔️✔️:   double rEN = ri * rj * (sqrt(Xi) - sqrt(Xj)) * (sqrt(Xi) - sqrt(Xj)) /
    // RDKit✔️✔️:                (Xi * ri + Xj * rj);
    let sqrt_delta = xi.sqrt() - xj.sqrt();
    let r_en = ri * rj * sqrt_delta * sqrt_delta / (xi * ri + xj * rj);

    // RDKit✔️✔️:   double res = ri + rj + rBO - rEN;
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    Ok(ri + rj + r_bo - r_en)
}

pub fn calc_bond_force_constant(
    rest_length: f64,
    end1_params: &AtomicParams,
    end2_params: &AtomicParams,
) -> Result<f64, UffUtilsError> {
    if rest_length <= 0.0 {
        return Err(UffUtilsError::InvalidRestLength { rest_length });
    }

    // RDKit✔️✔️: double calcBondForceConstant(double restLength, const AtomicParams *end1Params,
    // RDKit✔️✔️:                              const AtomicParams *end2Params) {
    // RDKit✔️✔️:   double res = 2.0 * Params::G * end1Params->Z1 * end2Params->Z1 /
    // RDKit✔️✔️:                (restLength * restLength * restLength);
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    Ok(2.0 * PARAMS_G * end1_params.z1 * end2_params.z1 / (rest_length * rest_length * rest_length))
}

pub fn calc_angle_force_constant(
    theta0: f64,
    bond_order12: f64,
    bond_order23: f64,
    at1_params: &AtomicParams,
    at2_params: &AtomicParams,
    at3_params: &AtomicParams,
) -> Result<f64, UffUtilsError> {
    // RDKit✔️✔️: double calcAngleForceConstant(double theta0, double bondOrder12,
    // RDKit✔️✔️:                               double bondOrder23, const AtomicParams *at1Params,
    // RDKit✔️✔️:                               const AtomicParams *at2Params,
    // RDKit✔️✔️:                               const AtomicParams *at3Params) {
    // RDKit✔️✔️:   double cosTheta0 = cos(theta0);
    let cos_theta0 = theta0.cos();
    // RDKit✔️✔️:   double r12 = calcBondRestLength(bondOrder12, at1Params, at2Params);
    let r12 = calc_bond_rest_length(bond_order12, at1_params, at2_params)?;
    // RDKit✔️✔️:   double r23 = calcBondRestLength(bondOrder23, at2Params, at3Params);
    let r23 = calc_bond_rest_length(bond_order23, at2_params, at3_params)?;
    // RDKit✔️✔️:   double r13 = sqrt(r12 * r12 + r23 * r23 - 2. * r12 * r23 * cosTheta0);
    let r13 = (r12 * r12 + r23 * r23 - 2.0 * r12 * r23 * cos_theta0).sqrt();
    // RDKit✔️✔️:   double beta = 2. * Params::G / (r12 * r23);
    let beta = 2.0 * PARAMS_G / (r12 * r23);

    // RDKit✔️✔️:   double preFactor = beta * at1Params->Z1 * at3Params->Z1 / int_pow<5>(r13);
    let pre_factor = beta * at1_params.z1 * at3_params.z1 / r13.powi(5);
    // RDKit✔️✔️:   double rTerm = r12 * r23;
    let r_term = r12 * r23;
    // RDKit✔️✔️:   double innerBit =
    // RDKit✔️✔️:       3. * rTerm * (1. - cosTheta0 * cosTheta0) - r13 * r13 * cosTheta0;
    let inner_bit = 3.0 * r_term * (1.0 - cos_theta0 * cos_theta0) - r13 * r13 * cos_theta0;
    // RDKit✔️✔️:   double res = preFactor * rTerm * innerBit;
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    Ok(pre_factor * r_term * inner_bit)
}

pub fn calc_nonbonded_minimum(at1_params: &AtomicParams, at2_params: &AtomicParams) -> f64 {
    // RDKit✔️✔️: double calcNonbondedMinimum(const AtomicParams *at1Params,
    // RDKit✔️✔️:                             const AtomicParams *at2Params) {
    // RDKit✔️✔️:   return sqrt(at1Params->x1 * at2Params->x1);
    // RDKit✔️✔️: }
    (at1_params.x1 * at2_params.x1).sqrt()
}

pub fn calc_nonbonded_depth(at1_params: &AtomicParams, at2_params: &AtomicParams) -> f64 {
    // RDKit✔️✔️: double calcNonbondedDepth(const AtomicParams *at1Params,
    // RDKit✔️✔️:                           const AtomicParams *at2Params) {
    // RDKit✔️✔️:   return sqrt(at1Params->D1 * at2Params->D1);
    // RDKit✔️✔️: }
    (at1_params.d1 * at2_params.d1).sqrt()
}

pub fn calc_inversion_coefficients(at2_atomic_num: i32, is_c_bound_to_o: bool) -> (f64, f64, f64, f64) {
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
    use super::*;

    const EPS: f64 = 1.0e-12;

    fn atomic_params(r1: f64, gmp_xi: f64) -> AtomicParams {
        AtomicParams {
            r1,
            theta0: 0.0,
            x1: 0.0,
            d1: 0.0,
            zeta: 0.0,
            z1: 1.0,
            v1: 0.0,
            u1: 0.0,
            gmp_xi,
            gmp_hardness: 0.0,
            gmp_radius: 0.0,
        }
    }

    fn assert_tuple_close(actual: (f64, f64, f64, f64), expected: (f64, f64, f64, f64)) {
        assert!((actual.0 - expected.0).abs() < EPS);
        assert!((actual.1 - expected.1).abs() < EPS);
        assert!((actual.2 - expected.2).abs() < EPS);
        assert!((actual.3 - expected.3).abs() < EPS);
    }

    fn expected_group5_inversion(w0_degrees: f64) -> (f64, f64, f64, f64) {
        let w0 = std::f64::consts::PI / 180.0 * w0_degrees;
        let c2 = 1.0;
        let c1 = -4.0 * w0.cos();
        let c0 = -(c1 * w0.cos() + c2 * (2.0 * w0).cos());
        let res = 22.0 / (c0 + c1 + c2) / 3.0;
        (res, c0, c1, c2)
    }

    #[test]
    fn uff_utils_calc_bond_rest_length_single_bond_same_params_has_no_corrections() {
        let params = atomic_params(0.757, 5.343);

        let rest_length = calc_bond_rest_length(1.0, &params, &params).expect("valid bond order");

        assert!((rest_length - 1.514).abs() < EPS);
    }

    #[test]
    fn uff_utils_calc_bond_rest_length_applies_bond_order_and_electronegativity_terms() {
        let end1 = atomic_params(0.757, 5.343);
        let end2 = atomic_params(0.658, 8.741);

        let rest_length = calc_bond_rest_length(2.0, &end1, &end2).expect("valid bond order");
        let r_bo = -PARAMS_LAMBDA * (end1.r1 + end2.r1) * 2.0_f64.ln();
        let sqrt_delta = end1.gmp_xi.sqrt() - end2.gmp_xi.sqrt();
        let r_en = end1.r1 * end2.r1 * sqrt_delta * sqrt_delta / (end1.gmp_xi * end1.r1 + end2.gmp_xi * end2.r1);
        let expected = end1.r1 + end2.r1 + r_bo - r_en;

        assert!((rest_length - expected).abs() < EPS);
        assert!(rest_length < end1.r1 + end2.r1);
    }

    #[test]
    fn uff_utils_calc_bond_rest_length_rejects_non_positive_bond_order() {
        let params = atomic_params(0.757, 5.343);

        assert_eq!(
            calc_bond_rest_length(0.0, &params, &params),
            Err(UffUtilsError::InvalidBondOrder { bond_order: 0.0 })
        );
        assert_eq!(
            calc_bond_rest_length(-1.0, &params, &params),
            Err(UffUtilsError::InvalidBondOrder { bond_order: -1.0 })
        );
    }

    #[test]
    fn uff_utils_calc_bond_force_constant_uses_source_prefactor_and_effective_charges() {
        let mut end1 = atomic_params(0.757, 5.343);
        let mut end2 = atomic_params(0.658, 8.741);
        end1.z1 = 1.912;
        end2.z1 = 2.3;

        let force_constant = calc_bond_force_constant(1.25, &end1, &end2).expect("valid rest length");
        let expected = 2.0 * PARAMS_G * end1.z1 * end2.z1 / (1.25_f64 * 1.25_f64 * 1.25_f64);

        assert!((force_constant - expected).abs() < EPS);
    }

    #[test]
    fn uff_utils_calc_bond_force_constant_rejects_non_positive_rest_length() {
        let params = atomic_params(0.757, 5.343);

        assert_eq!(
            calc_bond_force_constant(0.0, &params, &params),
            Err(UffUtilsError::InvalidRestLength { rest_length: 0.0 })
        );
        assert_eq!(
            calc_bond_force_constant(-1.0, &params, &params),
            Err(UffUtilsError::InvalidRestLength { rest_length: -1.0 })
        );
    }

    #[test]
    fn uff_utils_calc_angle_force_constant_matches_source_formula() {
        let mut at1 = atomic_params(0.757, 5.343);
        let mut at2 = atomic_params(0.700, 6.899);
        let mut at3 = atomic_params(0.658, 8.741);
        at1.z1 = 1.912;
        at2.z1 = 2.544;
        at3.z1 = 2.3;
        let theta0 = 109.47 * std::f64::consts::PI / 180.0;

        let force_constant = calc_angle_force_constant(theta0, 1.0, 1.0, &at1, &at2, &at3).expect("valid inputs");
        let cos_theta0 = theta0.cos();
        let r12 = calc_bond_rest_length(1.0, &at1, &at2).expect("valid bond order");
        let r23 = calc_bond_rest_length(1.0, &at2, &at3).expect("valid bond order");
        let r13 = (r12 * r12 + r23 * r23 - 2.0 * r12 * r23 * cos_theta0).sqrt();
        let beta = 2.0 * PARAMS_G / (r12 * r23);
        let r_term = r12 * r23;
        let inner_bit = 3.0 * r_term * (1.0 - cos_theta0 * cos_theta0) - r13 * r13 * cos_theta0;
        let expected = beta * at1.z1 * at3.z1 / r13.powi(5) * r_term * inner_bit;

        assert!((force_constant - expected).abs() < EPS);
    }

    #[test]
    fn uff_utils_calc_angle_force_constant_propagates_invalid_bond_order() {
        let params = atomic_params(0.757, 5.343);

        assert_eq!(
            calc_angle_force_constant(1.0, 0.0, 1.0, &params, &params, &params),
            Err(UffUtilsError::InvalidBondOrder { bond_order: 0.0 })
        );
        assert_eq!(
            calc_angle_force_constant(1.0, 1.0, -1.0, &params, &params, &params),
            Err(UffUtilsError::InvalidBondOrder { bond_order: -1.0 })
        );
    }

    #[test]
    fn uff_utils_calc_nonbonded_minimum_uses_geometric_mean() {
        let mut at1 = atomic_params(0.757, 5.343);
        let mut at2 = atomic_params(0.658, 8.741);
        at1.x1 = 3.851;
        at2.x1 = 3.5;

        let minimum = calc_nonbonded_minimum(&at1, &at2);
        let expected = (at1.x1 * at2.x1).sqrt();

        assert!((minimum - expected).abs() < EPS);
    }

    #[test]
    fn uff_utils_calc_nonbonded_minimum_preserves_zero_boundary() {
        let mut at1 = atomic_params(0.757, 5.343);
        let mut at2 = atomic_params(0.658, 8.741);
        at1.x1 = 0.0;
        at2.x1 = 3.5;

        assert_eq!(calc_nonbonded_minimum(&at1, &at2), 0.0);
    }

    #[test]
    fn uff_utils_calc_nonbonded_minimum_preserves_negative_product_nan() {
        let mut at1 = atomic_params(0.757, 5.343);
        let mut at2 = atomic_params(0.658, 8.741);
        at1.x1 = -1.0;
        at2.x1 = 3.5;

        assert!(calc_nonbonded_minimum(&at1, &at2).is_nan());
    }

    #[test]
    fn uff_utils_calc_nonbonded_depth_uses_geometric_mean() {
        let mut at1 = atomic_params(0.757, 5.343);
        let mut at2 = atomic_params(0.658, 8.741);
        at1.d1 = 0.105;
        at2.d1 = 0.06;

        let depth = calc_nonbonded_depth(&at1, &at2);
        let expected = (at1.d1 * at2.d1).sqrt();

        assert!((depth - expected).abs() < EPS);
    }

    #[test]
    fn uff_utils_calc_nonbonded_depth_preserves_zero_boundary() {
        let mut at1 = atomic_params(0.757, 5.343);
        let mut at2 = atomic_params(0.658, 8.741);
        at1.d1 = 0.0;
        at2.d1 = 0.06;

        assert_eq!(calc_nonbonded_depth(&at1, &at2), 0.0);
    }

    #[test]
    fn uff_utils_calc_nonbonded_depth_preserves_negative_product_nan() {
        let mut at1 = atomic_params(0.757, 5.343);
        let mut at2 = atomic_params(0.658, 8.741);
        at1.d1 = -1.0;
        at2.d1 = 0.06;

        assert!(calc_nonbonded_depth(&at1, &at2).is_nan());
    }

    #[test]
    fn uff_utils_calc_inversion_coefficients_handles_sp2_carbon_nitrogen_and_oxygen() {
        let expected = (2.0, 1.0, -1.0, 0.0);

        assert_tuple_close(calc_inversion_coefficients(6, false), expected);
        assert_tuple_close(calc_inversion_coefficients(7, false), expected);
        assert_tuple_close(calc_inversion_coefficients(8, false), expected);
    }

    #[test]
    fn uff_utils_calc_inversion_coefficients_scales_carbon_bound_to_oxygen() {
        assert_tuple_close(calc_inversion_coefficients(6, true), (50.0 / 3.0, 1.0, -1.0, 0.0));
    }

    #[test]
    fn uff_utils_calc_inversion_coefficients_handles_group5_switch_cases() {
        assert_tuple_close(
            calc_inversion_coefficients(15, false),
            expected_group5_inversion(84.4339),
        );
        assert_tuple_close(
            calc_inversion_coefficients(33, false),
            expected_group5_inversion(86.9735),
        );
        assert_tuple_close(
            calc_inversion_coefficients(51, false),
            expected_group5_inversion(87.7047),
        );
        assert_tuple_close(calc_inversion_coefficients(83, false), expected_group5_inversion(90.0));
    }

    #[test]
    fn uff_utils_calc_inversion_coefficients_preserves_default_switch_fallthrough() {
        assert_tuple_close(calc_inversion_coefficients(0, false), expected_group5_inversion(1.0));
    }
}
