// Copyright (C) 2004-2006 Rational Discovery LLC
//
// This file is part of the RDKit-derived force-field implementation and is
// covered by the BSD license in the pinned RDKit source tree.

use super::params::{AtomicParams, PARAMS_G, PARAMS_LAMBDA};
use crate::kernel::{BondIndexArgument, EvaluationContext, ForceFieldKernelError};

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) enum BondMathError {
    InvalidBondOrder { bond_order: f64 },
}

impl std::fmt::Display for BondMathError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Debug::fmt(self, formatter)
    }
}

impl std::error::Error for BondMathError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        None
    }
}

pub(super) fn calc_bond_rest_length(
    bond_order: f64,
    end1_params: &AtomicParams,
    end2_params: &AtomicParams,
) -> Result<f64, BondMathError> {
    // RDKit✔️✔️: double calcBondRestLength(double bondOrder, const AtomicParams *end1Params,
    // RDKit✔️✔️:                           const AtomicParams *end2Params) {
    // RDKit✔️✔️:   PRECONDITION(bondOrder > 0, "bad bond order");
    // RDGeneral/Invariant.h expands PRECONDITION through `if (!(expr))`, so
    // unordered NaN bond orders fail the same source condition.
    if !(bond_order > 0.0) {
        return Err(BondMathError::InvalidBondOrder { bond_order });
    }

    // Rust references enforce the source's two non-null pointer preconditions.
    // RDKit✔️✔️:   PRECONDITION(end1Params, "bad params pointer");
    // RDKit✔️✔️:   PRECONDITION(end2Params, "bad params pointer");
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
    let res = ri + rj + r_bo - r_en;
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    Ok(res)
}

pub(super) fn calc_bond_force_constant(
    rest_length: f64,
    end1_params: &AtomicParams,
    end2_params: &AtomicParams,
) -> f64 {
    // RDKit✔️✔️: double calcBondForceConstant(double restLength, const AtomicParams *end1Params,
    // RDKit✔️✔️:                              const AtomicParams *end2Params) {
    // RDKit✔️✔️:   double res = 2.0 * Params::G * end1Params->Z1 * end2Params->Z1 /
    // RDKit✔️✔️:                (restLength * restLength * restLength);
    let res = 2.0 * PARAMS_G * end1_params.z1 * end2_params.z1
        / (rest_length * rest_length * rest_length);
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    res
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct BondStretchContrib {
    end1_idx: u32,
    end2_idx: u32,
    rest_len: f64,
    force_constant: f64,
}

impl BondStretchContrib {
    #[cfg(test)]
    pub(crate) fn cf3d_frag_accept_stored_identity(
        &self,
    ) -> crate::kernel::Cf3dFragAcceptContributionIdentity {
        crate::kernel::Cf3dFragAcceptContributionIdentity::BondStretch {
            end1_idx: self.end1_idx,
            end2_idx: self.end2_idx,
            rest_len: self.rest_len,
            force_constant: self.force_constant,
        }
    }

    pub(crate) fn new(
        positions: &[&mut [f64]],
        idx1: u32,
        idx2: u32,
        bond_order: f64,
        end1_params: &AtomicParams,
        end2_params: &AtomicParams,
    ) -> Result<Self, ForceFieldKernelError> {
        // RDKit✔️✔️: BondStretchContrib::BondStretchContrib(ForceField *owner, unsigned int idx1,
        // RDKit✔️✔️:                                        unsigned int idx2, double bondOrder,
        // RDKit✔️✔️:                                        const AtomicParams *end1Params,
        // RDKit✔️✔️:                                        const AtomicParams *end2Params) {
        // RDKit✔️✔️:   PRECONDITION(owner, "bad owner");
        // RDKit✔️✔️:   PRECONDITION(end1Params, "bad params pointer");
        // RDKit✔️✔️:   PRECONDITION(end2Params, "bad params pointer");
        // Safe borrows make the owner and parameter-pointer null states
        // unrepresentable; this view reads only the owner's current positions.
        // RDKit✔️✔️:   URANGE_CHECK(idx1, owner->positions().size());
        if idx1 as usize >= positions.len() {
            return Err(ForceFieldKernelError::BondIndexOutOfRange {
                argument: BondIndexArgument::First,
                index: idx1,
                upper_bound: positions.len(),
            });
        }
        // RDKit✔️✔️:   URANGE_CHECK(idx2, owner->positions().size());
        if idx2 as usize >= positions.len() {
            return Err(ForceFieldKernelError::BondIndexOutOfRange {
                argument: BondIndexArgument::Second,
                index: idx2,
                upper_bound: positions.len(),
            });
        }

        // RDKit✔️✔️:   dp_forceField = owner;
        // RDKit✔️✔️:   d_end1Idx = idx1;
        // RDKit✔️✔️:   d_end2Idx = idx2;
        // The Rust contribution owns source indices and computed scalar
        // parameters, and retains no pointer or borrow of its ForceField.
        // RDKit✔️✔️:   d_restLen = Utils::calcBondRestLength(bondOrder, end1Params, end2Params);
        let rest_len = calc_bond_rest_length(bond_order, end1_params, end2_params).map_err(
            |error| match error {
                BondMathError::InvalidBondOrder { .. } => ForceFieldKernelError::BadBondOrder,
            },
        )?;
        // RDKit✔️✔️:   d_forceConstant =
        // RDKit✔️✔️:       Utils::calcBondForceConstant(d_restLen, end1Params, end2Params);
        // RDKit✔️✔️: }
        let force_constant = calc_bond_force_constant(rest_len, end1_params, end2_params);

        Ok(Self {
            end1_idx: idx1,
            end2_idx: idx2,
            rest_len,
            force_constant,
        })
    }

    pub(crate) fn get_energy(
        &self,
        context: &mut EvaluationContext<'_>,
    ) -> Result<f64, ForceFieldKernelError> {
        // RDKit✔️✔️: double BondStretchContrib::getEnergy(double *pos) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // The borrowed context guarantees an owner-backed cache and positions
        // for this evaluation without retaining either on the contribution.
        // RDKit✔️✔️:   double distTerm =
        // RDKit✔️✔️:       dp_forceField->distance(d_end1Idx, d_end2Idx, pos) - d_restLen;
        let dist_term = context.distance(self.end1_idx, self.end2_idx)? - self.rest_len;
        // RDKit✔️✔️:   double res = 0.5 * d_forceConstant * distTerm * distTerm;
        let res = 0.5 * self.force_constant * dist_term * dist_term;
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        Ok(res)
    }

    pub(crate) fn get_grad(
        &self,
        context: &mut EvaluationContext<'_>,
        gradient: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        // BEGIN RDKIT CPP FUNCTION UFF::BondStretchContrib::getGrad (BondStretch.cpp:78-102)
        // RDKit✔️✔️: void BondStretchContrib::getGrad(double *pos, double *grad) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:   PRECONDITION(grad, "bad vector");
        // The borrowed context supplies a valid owner/cache, and Rust slices
        // supply non-null coordinate and gradient pointers.
        // RDKit✔️✔️:   double dist = dp_forceField->distance(d_end1Idx, d_end2Idx, pos);
        let dist = context.distance(self.end1_idx, self.end2_idx)?;
        // RDKit✔️✔️:   double preFactor = d_forceConstant * (dist - d_restLen);
        let pre_factor = self.force_constant * (dist - self.rest_len);

        // RDKit✔️✔️:   // std::cout << "\tDist("<<d_end1Idx<<","<<d_end2Idx<<") " << dist <<
        // RDKit✔️✔️:   // std::endl;
        // RDKit✔️✔️:   double *end1Coords = &(pos[3 * d_end1Idx]);
        // RDKit✔️✔️:   double *end2Coords = &(pos[3 * d_end2Idx]);
        // The source's UFF bond gradient always addresses three coordinates
        // per endpoint, independent of the general force-field dimension.
        let end1_offset = self.end1_idx.wrapping_mul(3) as usize;
        let end2_offset = self.end2_idx.wrapping_mul(3) as usize;
        let coordinates = context.coordinates();
        // RDKit✔️✔️:   for (int i = 0; i < 3; i++) {
        for i in 0..3 {
            // RDKit✔️✔️:     double dGrad;
            // RDKit✔️✔️:     if (dist > 0.0) {
            let d_grad = if dist > 0.0 {
                // RDKit✔️✔️:       dGrad = preFactor * (end1Coords[i] - end2Coords[i]) / dist;
                pre_factor * (coordinates[end1_offset + i] - coordinates[end2_offset + i]) / dist
            // RDKit✔️✔️:     } else {
            } else {
                // RDKit✔️✔️:       // move a small amount in an arbitrary direction
                // RDKit✔️✔️:       dGrad = d_forceConstant * .01;
                self.force_constant * 0.01
            };
            // RDKit✔️✔️:     grad[3 * d_end1Idx + i] += dGrad;
            // RDKit✔️✔️:     grad[3 * d_end2Idx + i] -= dGrad;
            gradient[end1_offset + i] += d_grad;
            gradient[end2_offset + i] -= d_grad;
        }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION UFF::BondStretchContrib::getGrad
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::error::Error as _;

    use super::*;

    #[test]
    fn uff_error_e02_bond_error_preserves_float_payload_bits() {
        // The positive InvalidBondOrder row is synthetic trait-dispatch
        // coverage only; source callers reject zero, negative, and NaN.
        let manual_rows = [
            (
                BondMathError::InvalidBondOrder { bond_order: 1.25 },
                0x3ff4_0000_0000_0000,
            ),
            (
                BondMathError::InvalidBondOrder { bond_order: 0.0 },
                0x0000_0000_0000_0000,
            ),
            (
                BondMathError::InvalidBondOrder { bond_order: -0.0 },
                0x8000_0000_0000_0000,
            ),
            (
                BondMathError::InvalidBondOrder { bond_order: -1.25 },
                0xbff4_0000_0000_0000,
            ),
            (
                BondMathError::InvalidBondOrder {
                    bond_order: f64::from_bits(0x7ff8_0000_0000_0042),
                },
                0x7ff8_0000_0000_0042,
            ),
        ];
        for (error, expected_bits) in manual_rows {
            let erased: &(dyn std::error::Error + 'static) = &error;
            let downcast = erased.downcast_ref::<BondMathError>().unwrap();
            assert!(std::ptr::eq(downcast, &error));
            assert_eq!(
                match downcast {
                    BondMathError::InvalidBondOrder { bond_order } => bond_order.to_bits(),
                },
                expected_bits
            );
            assert!(erased.source().is_none());
            assert_eq!(error.to_string(), format!("{error:?}"));
        }

        let carbon = atomic_params(0.757, 5.343, 1.912);
        for (bond_order, expected_bits) in [
            (0.0, 0x0000_0000_0000_0000),
            (-0.0, 0x8000_0000_0000_0000),
            (-1.25, 0xbff4_0000_0000_0000),
            (f64::from_bits(0x7ff8_0000_0000_0042), 0x7ff8_0000_0000_0042),
        ] {
            let error = calc_bond_rest_length(bond_order, &carbon, &carbon).unwrap_err();
            let erased: &(dyn std::error::Error + 'static) = &error;
            assert!(std::ptr::eq(
                erased.downcast_ref::<BondMathError>().unwrap(),
                &error
            ));
            match error {
                BondMathError::InvalidBondOrder { bond_order: actual } => {
                    assert_eq!(actual.to_bits(), expected_bits);
                }
            }
            assert!(erased.source().is_none());
        }
    }

    fn atomic_params(r1: f64, gmp_xi: f64, z1: f64) -> AtomicParams {
        AtomicParams {
            r1,
            theta0: 0.0,
            x1: 0.0,
            d1: 0.0,
            zeta: 0.0,
            z1,
            v1: 0.0,
            u1: 0.0,
            gmp_xi,
            gmp_hardness: 0.0,
            gmp_radius: 0.0,
        }
    }

    fn assert_close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "expected {expected:.17}, got {actual:.17} (tolerance {tolerance})"
        );
    }

    #[test]
    fn cf3d_u04_rdkit_bond_order_and_electronegativity_fixtures() {
        // Fixed rest-length and force-constant values from pinned RDKit
        // Code/ForceField/UFF/testUFFForceField.cpp::testUFF1.
        let sp3_carbon = atomic_params(0.757, 5.343, 1.912);
        let sp3_c_c = calc_bond_rest_length(1.0, &sp3_carbon, &sp3_carbon)
            .expect("positive source bond order");
        assert_close(sp3_c_c, 1.514, 1.0e-12);
        assert_close(
            calc_bond_force_constant(sp3_c_c, &sp3_carbon, &sp3_carbon),
            699.5918,
            1.0e-3,
        );

        let sp2_carbon = atomic_params(0.732, 5.343, 1.912);
        let sp2_c_c = calc_bond_rest_length(2.0, &sp2_carbon, &sp2_carbon)
            .expect("positive source bond order");
        assert_close(sp2_c_c, 1.32883, 1.0e-5);
        assert_close(
            calc_bond_force_constant(sp2_c_c, &sp2_carbon, &sp2_carbon),
            1034.69,
            1.0e-2,
        );

        let sp3_nitrogen = atomic_params(0.700, 6.899, 2.544);
        let c_n = calc_bond_rest_length(1.0, &sp3_carbon, &sp3_nitrogen)
            .expect("positive source bond order");
        assert_close(c_n, 1.451071, 1.0e-5);
        assert_close(
            calc_bond_force_constant(c_n, &sp3_carbon, &sp3_nitrogen),
            1057.27,
            1.0e-2,
        );

        let resonance_carbon = atomic_params(0.729, 5.343, 1.912);
        let resonance_nitrogen = atomic_params(0.699, 6.899, 2.544);
        let amide = calc_bond_rest_length(1.41, &resonance_carbon, &resonance_nitrogen)
            .expect("positive source amide bond order");
        assert_close(amide, 1.357, 1.0e-3);
        assert_close(
            calc_bond_force_constant(amide, &resonance_carbon, &resonance_nitrogen),
            1293.0,
            1.0,
        );
    }

    #[test]
    fn cf3d_u04_bond_order_source_boundary_and_ieee_cases() {
        let carbon = atomic_params(0.757, 5.343, 1.912);

        // The pinned PRECONDITION is `bondOrder > 0`; its macro rejects
        // unordered NaN as well as zero and negative values.
        for bond_order in [0.0, -0.0, -1.0, f64::NAN] {
            match calc_bond_rest_length(bond_order, &carbon, &carbon) {
                Err(BondMathError::InvalidBondOrder { bond_order: actual }) => {
                    if bond_order.is_nan() {
                        assert!(actual.is_nan());
                    } else {
                        assert_eq!(actual.to_bits(), bond_order.to_bits());
                    }
                }
                Ok(value) => panic!("source precondition accepted {bond_order}: {value}"),
            }
        }

        // Positive values below, at, and above one follow the same source
        // formula; the half-order value fixes the logarithmic correction.
        assert_close(
            calc_bond_rest_length(0.5, &carbon, &carbon).expect("positive bond order"),
            1.6537833875381853,
            1.0e-12,
        );
        assert_close(
            calc_bond_rest_length(1.0, &carbon, &carbon).expect("positive bond order"),
            1.514,
            1.0e-12,
        );
        assert_close(
            calc_bond_rest_length(2.0, &carbon, &carbon).expect("positive bond order"),
            1.3742166124618147,
            1.0e-12,
        );
        assert_close(
            calc_bond_rest_length(f64::MIN_POSITIVE, &carbon, &carbon)
                .expect("smallest positive normal order"),
            144.37262206402536,
            1.0e-12,
        );
        let infinite_order = calc_bond_rest_length(f64::INFINITY, &carbon, &carbon)
            .expect("positive infinity meets the source precondition");
        assert!(infinite_order.is_infinite());
        assert!(infinite_order.is_sign_negative());
    }

    #[test]
    fn cf3d_u04_force_constant_preserves_source_ieee_division() {
        let carbon = atomic_params(0.757, 5.343, 1.912);
        assert_close(
            calc_bond_force_constant(1.0, &carbon, &carbon),
            2427.85270528,
            1.0e-12,
        );
        assert!(calc_bond_force_constant(0.0, &carbon, &carbon).is_infinite());
        assert!(calc_bond_force_constant(0.0, &carbon, &carbon).is_sign_positive());
        assert_close(
            calc_bond_force_constant(-1.0, &carbon, &carbon),
            -2427.85270528,
            1.0e-12,
        );
        let negative_zero = calc_bond_force_constant(-0.0, &carbon, &carbon);
        assert!(negative_zero.is_infinite());
        assert!(negative_zero.is_sign_negative());
        assert_eq!(
            calc_bond_force_constant(f64::INFINITY, &carbon, &carbon),
            0.0
        );
        assert!(calc_bond_force_constant(f64::NAN, &carbon, &carbon).is_nan());
    }

    #[test]
    fn cf3d_u04_atomic_parameter_values_flow_through_source_arithmetic() {
        // RDKit's helpers add no finiteness, sign, or nonzero parameter guard.
        let zero_xi = atomic_params(1.0, 0.0, 0.0);
        assert!(
            calc_bond_rest_length(1.0, &zero_xi, &zero_xi)
                .expect("no source Xi guard")
                .is_nan()
        );

        let negative_xi = atomic_params(1.0, -1.0, 0.0);
        assert!(
            calc_bond_rest_length(1.0, &negative_xi, &negative_xi)
                .expect("no source Xi guard")
                .is_nan()
        );

        let zero_charge = atomic_params(1.0, 1.0, 0.0);
        assert_eq!(
            calc_bond_force_constant(1.0, &zero_charge, &zero_charge),
            0.0
        );
    }
}
