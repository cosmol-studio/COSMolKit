// $Id$
//
// Copyright (C) 2004-2008 Greg Landrum and Rational Discovery LLC
//
// @@ All Rights Reserved @@
// This file is part of the RDKit.
// The contents are covered by the terms of the BSD license
// which is included in the file license.txt, found at the root
// of the RDKit source tree.

use crate::kernel::{
    BondIndexArgument, EvaluationContext, ForceFieldContribution, ForceFieldKernelError,
};

use super::params::AtomicParams;

pub(super) fn calc_nonbonded_minimum(at1_params: &AtomicParams, at2_params: &AtomicParams) -> f64 {
    // RDKit❗✔️: double calcNonbondedMinimum(const AtomicParams *at1Params,
    // RDKit❗✔️:                             const AtomicParams *at2Params) {
    // RDKit❗✔️:   return sqrt(at1Params->x1 * at2Params->x1);
    // RDKit❗✔️: }
    (at1_params.x1 * at2_params.x1).sqrt()
}

pub(super) fn calc_nonbonded_depth(at1_params: &AtomicParams, at2_params: &AtomicParams) -> f64 {
    // RDKit❗✔️: double calcNonbondedDepth(const AtomicParams *at1Params,
    // RDKit❗✔️:                           const AtomicParams *at2Params) {
    // RDKit❗✔️:   return sqrt(at1Params->D1 * at2Params->D1);
    // RDKit❗✔️: }
    (at1_params.d1 * at2_params.d1).sqrt()
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct VdwContrib {
    at1_idx: u32,
    at2_idx: u32,
    x_ij: f64,
    well_depth: f64,
    threshold: f64,
}

impl VdwContrib {
    pub(super) fn new(
        positions: &[&mut [f64]],
        idx1: u32,
        idx2: u32,
        at1_params: &AtomicParams,
        at2_params: &AtomicParams,
    ) -> Result<Self, ForceFieldKernelError> {
        // RDKit❗✔️: vdWContrib(ForceField *owner, unsigned int idx1,
        // RDKit❗✔️:              unsigned int idx2, const AtomicParams *at1Params,
        // RDKit❗✔️:              const AtomicParams *at2Params,
        // RDKit❗✔️:              double threshMultiplier = 10.0);
        // Rust references make the source's non-null owner and parameter
        // preconditions unrepresentable; the source header supplies this
        // default when callers omit the threshold multiplier.
        Self::new_with_threshold_multiplier(positions, idx1, idx2, at1_params, at2_params, 10.0)
    }

    pub(super) fn new_with_threshold_multiplier(
        positions: &[&mut [f64]],
        idx1: u32,
        idx2: u32,
        at1_params: &AtomicParams,
        at2_params: &AtomicParams,
        threshold_multiplier: f64,
    ) -> Result<Self, ForceFieldKernelError> {
        // BEGIN RDKIT CPP FUNCTION UFF::vdWContrib::vdWContrib (Nonbonded.cpp:28-50)
        // RDKit❗✔️: vdWContrib::vdWContrib(ForceField *owner, unsigned int idx1, unsigned int idx2,
        // RDKit❗✔️:                        const AtomicParams *at1Params,
        // RDKit❗✔️:                        const AtomicParams *at2Params, double threshMultiplier) {
        // RDKit❗✔️:   PRECONDITION(owner, "bad owner");
        // RDKit❗✔️:   PRECONDITION(at1Params, "bad params pointer");
        // RDKit❗✔️:   PRECONDITION(at2Params, "bad params pointer");
        // Non-null Rust borrows preserve these source preconditions.

        // RDKit❗✔️:   URANGE_CHECK(idx1, owner->positions().size());
        if idx1 as usize >= positions.len() {
            return Err(ForceFieldKernelError::BondIndexOutOfRange {
                argument: BondIndexArgument::First,
                index: idx1,
                upper_bound: positions.len(),
            });
        }
        // RDKit❗✔️:   URANGE_CHECK(idx2, owner->positions().size());
        if idx2 as usize >= positions.len() {
            return Err(ForceFieldKernelError::BondIndexOutOfRange {
                argument: BondIndexArgument::Second,
                index: idx2,
                upper_bound: positions.len(),
            });
        }

        // RDKit❗✔️:   dp_forceField = owner;
        // RDKit❗✔️:   d_at1Idx = idx1;
        // RDKit❗✔️:   d_at2Idx = idx2;
        // This contribution keeps source indices and computed scalars only;
        // evaluation receives its owner-backed state through a short borrow.
        let at1_idx = idx1;
        let at2_idx = idx2;

        // RDKit❗✔️:   // UFF uses the geometric mean of the vdW parameters:
        // RDKit❗✔️:   d_xij = Utils::calcNonbondedMinimum(at1Params, at2Params);
        let x_ij = calc_nonbonded_minimum(at1_params, at2_params);
        // RDKit❗✔️:   d_wellDepth = Utils::calcNonbondedDepth(at1Params, at2Params);
        let well_depth = calc_nonbonded_depth(at1_params, at2_params);
        // RDKit❗✔️:   d_thresh = threshMultiplier * d_xij;
        let threshold = threshold_multiplier * x_ij;

        // RDKit❗✔️:   // std::cerr << "  non-bonded: " << idx1 << "-" << idx2 << " " << d_xij << " "
        // RDKit❗✔️:   // << d_wellDepth << " " << d_thresh << std::endl;
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION UFF::vdWContrib::vdWContrib
        Ok(Self {
            at1_idx,
            at2_idx,
            x_ij,
            well_depth,
            threshold,
        })
    }

    pub(super) fn get_energy(
        &self,
        context: &mut EvaluationContext<'_>,
    ) -> Result<f64, ForceFieldKernelError> {
        // BEGIN RDKIT CPP FUNCTION UFF::vdWContrib::getEnergy (Nonbonded.cpp:54-69)
        // RDKit❗✔️: double vdWContrib::getEnergy(double *pos) const {
        // RDKit❗✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit❗✔️:   PRECONDITION(pos, "bad vector");
        // The borrowed context supplies the owner-backed positions and cache;
        // the contribution does not keep a raw owner pointer.

        // RDKit❗✔️:   double dist = dp_forceField->distance(d_at1Idx, d_at2Idx, pos);
        let dist = context.distance(self.at1_idx, self.at2_idx)?;
        // RDKit❗✔️:   if (dist > d_thresh || dist <= 0.0) {
        if dist > self.threshold || dist <= 0.0 {
            // RDKit❗✔️:     return 0.0;
            // RDKit❗✔️:   }
            return Ok(0.0);
        }

        // RDKit❗✔️:   double r = d_xij / dist;
        let r = self.x_ij / dist;
        // BEGIN RDKIT CPP HELPER int_pow (RDGeneral/utils.h:86-104)
        // RDKit❗✔️: template <unsigned n>
        // RDKit❗✔️: inline double int_pow(double x) {
        // RDKit❗✔️:   double half = int_pow<n / 2>(x);
        // RDKit❗✔️:   if (n % 2 == 0) {  // even
        // RDKit❗✔️:     return half * half;
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     return half * half * x;
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // RDKit❗✔️: template <>
        // RDKit❗✔️: inline double int_pow<0>(double) {
        // RDKit❗✔️:   return 1;
        // RDKit❗✔️: }
        // RDKit❗✔️: template <>
        // RDKit❗✔️: inline double int_pow<1>(double x) {
        // RDKit❗✔️:   return x;  // this does a series of muls
        // RDKit❗✔️: }
        // END RDKIT CPP HELPER int_pow
        // The source's fixed exponent recursion reduces to these same
        // multiplications for exponent six, without a runtime loop.
        let r3 = r * r * r;
        let r6 = r3 * r3;
        // RDKit❗✔️:   double r12 = r6 * r6;
        let r12 = r6 * r6;
        // RDKit❗✔️:   double res = d_wellDepth * (r12 - 2.0 * r6);
        let res = self.well_depth * (r12 - 2.0 * r6);
        // RDKit❗✔️:   // if(d_at1Idx==12 && d_at2Idx==21 ) std::cerr << "     >: " << d_at1Idx <<
        // RDKit❗✔️:   // "-" << d_at2Idx << " " << r << " = " << res << std::endl;
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION UFF::vdWContrib::getEnergy
        Ok(res)
    }

    pub(super) fn get_grad(
        &self,
        context: &mut EvaluationContext<'_>,
        gradient: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        // BEGIN RDKIT CPP FUNCTION UFF::vdWContrib::getGrad (Nonbonded.cpp:71-103)
        // RDKit❗✔️: void vdWContrib::getGrad(double *pos, double *grad) const {
        // RDKit❗✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit❗✔️:   PRECONDITION(pos, "bad vector");
        // RDKit❗✔️:   PRECONDITION(grad, "bad vector");
        // The borrowed context supplies owner-backed positions/cache, and
        // Rust slices make null coordinate and gradient pointers impossible.

        // RDKit❗✔️:   double dist = dp_forceField->distance(d_at1Idx, d_at2Idx, pos);
        let dist = context.distance(self.at1_idx, self.at2_idx)?;
        // RDKit❗✔️:   if (dist > d_thresh) {
        if dist > self.threshold {
            // RDKit❗✔️:     return;
            // RDKit❗✔️:   }
            return Ok(());
        }

        // RDKit❗✔️:   if (dist <= 0) {
        if dist <= 0.0 {
            // RDKit❗✔️:     for (int i = 0; i < 3; i++) {
            let at1_offset = self.at1_idx.wrapping_mul(3) as usize;
            let at2_offset = self.at2_idx.wrapping_mul(3) as usize;
            for i in 0..3 {
                // RDKit❗✔️:       // move in an arbitrary direction
                // RDKit❗✔️:       double dGrad = 100.0;
                let d_grad = 100.0;
                // RDKit❗✔️:       grad[3 * d_at1Idx + i] += dGrad;
                gradient[at1_offset + i] += d_grad;
                // RDKit❗✔️:       grad[3 * d_at2Idx + i] -= dGrad;
                gradient[at2_offset + i] -= d_grad;
            }
            // RDKit❗✔️:     return;
            // RDKit❗✔️:   }
            return Ok(());
        }

        // RDKit❗✔️:   double r = d_xij / dist;
        let r = self.x_ij / dist;
        // BEGIN RDKIT CPP HELPER int_pow (RDGeneral/utils.h:86-104)
        // RDKit❗✔️: template <unsigned n>
        // RDKit❗✔️: inline double int_pow(double x) {
        // RDKit❗✔️:   double half = int_pow<n / 2>(x);
        // RDKit❗✔️:   if (n % 2 == 0) {  // even
        // RDKit❗✔️:     return half * half;
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     return half * half * x;
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // RDKit❗✔️: template <>
        // RDKit❗✔️: inline double int_pow<0>(double) {
        // RDKit❗✔️:   return 1;
        // RDKit❗✔️: }
        // RDKit❗✔️: template <>
        // RDKit❗✔️: inline double int_pow<1>(double x) {
        // RDKit❗✔️:   return x;  // this does a series of muls
        // RDKit❗✔️: }
        // END RDKIT CPP HELPER int_pow
        // RDKit❗✔️:   double r7 = int_pow<7>(r);
        let r3 = r * r * r;
        let r7 = r3 * r3 * r;
        // RDKit❗✔️:   double r13 = int_pow<13>(r);
        let r6 = r3 * r3;
        let r13 = r6 * r6 * r;
        // RDKit❗✔️:   double preFactor = 12. * d_wellDepth / d_xij * (r7 - r13);
        let pre_factor = 12.0 * self.well_depth / self.x_ij * (r7 - r13);

        // RDKit❗✔️:   double *at1Coords = &(pos[3 * d_at1Idx]);
        let at1_offset = self.at1_idx.wrapping_mul(3) as usize;
        // RDKit❗✔️:   double *at2Coords = &(pos[3 * d_at2Idx]);
        let at2_offset = self.at2_idx.wrapping_mul(3) as usize;
        let coordinates = context.coordinates();
        // RDKit❗✔️:   for (int i = 0; i < 3; i++) {
        for i in 0..3 {
            // RDKit❗✔️:     double dGrad = preFactor * (at1Coords[i] - at2Coords[i]) / dist;
            let d_grad =
                pre_factor * (coordinates[at1_offset + i] - coordinates[at2_offset + i]) / dist;
            // RDKit❗✔️:     grad[3 * d_at1Idx + i] += dGrad;
            gradient[at1_offset + i] += d_grad;
            // RDKit❗✔️:     grad[3 * d_at2Idx + i] -= dGrad;
            gradient[at2_offset + i] -= d_grad;
        }
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION UFF::vdWContrib::getGrad
        Ok(())
    }
}

impl ForceFieldContribution for VdwContrib {
    #[cfg(test)]
    fn cf3d_frag_accept_test_identity(&self) -> crate::kernel::Cf3dFragAcceptContributionIdentity {
        crate::kernel::Cf3dFragAcceptContributionIdentity::Vdw {
            at1_idx: self.at1_idx,
            at2_idx: self.at2_idx,
            x_ij: self.x_ij,
            well_depth: self.well_depth,
            threshold: self.threshold,
        }
    }

    fn get_energy(
        &self,
        context: &mut EvaluationContext<'_>,
    ) -> Result<f64, ForceFieldKernelError> {
        // BEGIN RDKIT CPP CALL ForceFields::ForceField::calcEnergy(pos)
        // (ForceField/ForceField.cpp:320-325)
        // RDKit❗✔️:     double E = (*contrib)->getEnergy(pos);
        // RDKit❗✔️:     res += E;
        // END RDKIT CPP CALL ForceFields::ForceField::calcEnergy(pos)
        // Behavior marker — RDKit❗✔️: delegate to the existing scalar evaluator
        // and retain its typed distance/cache error unchanged.
        // Complexity marker — RDKit✔️✔️: one direct call and Result propagation;
        // this adds no scan, formula work, temporary buffer, or allocation.
        VdwContrib::get_energy(self, context)
    }

    fn get_grad(
        &self,
        context: &mut EvaluationContext<'_>,
        gradient: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        // BEGIN RDKIT CPP CALL ForceFields::ForceField::calcGrad(pos)
        // (ForceField/ForceField.cpp:361-364)
        // RDKit❗✔️:     (*contrib)->getGrad(pos, grad);
        // END RDKIT CPP CALL ForceFields::ForceField::calcGrad(pos)
        // Behavior marker — RDKit❗✔️: the source callback is void; preserve the
        // contribution's additive gradient and return only its typed kernel error.
        // Complexity marker — RDKit✔️✔️: one direct call; no second distance or
        // gradient pass is introduced.
        VdwContrib::get_grad(self, context, gradient)
    }

    fn copy(&self) -> Box<dyn ForceFieldContribution> {
        // BEGIN RDKIT CPP FUNCTION ForceFields::UFF::vdWContrib::copy
        // (ForceField/UFF/Nonbonded.h:49-51)
        // RDKit❗✔️:   vdWContrib *copy() const override { return new vdWContrib(*this); }
        // END RDKIT CPP FUNCTION ForceFields::UFF::vdWContrib::copy
        // Behavior marker — RDKit❗✔️: copy both indices and all three computed
        // scalar parameters; owner-derived state remains borrowed from the kernel.
        // Complexity marker — RDKit✔️✔️: fixed-size copy and one trait allocation,
        // matching the source's heap-allocated virtual copy.
        Box::new(*self)
    }
}

#[cfg(test)]
mod tests {
    use super::{VdwContrib, calc_nonbonded_depth, calc_nonbonded_minimum};
    use crate::kernel::{
        BondIndexArgument, EvaluationContext, ForceField, ForceFieldKernelError,
        cf3d_bld_b05_calc_energy, cf3d_bld_b05_calc_grad, cf3d_bld_b05_copy_force_field,
    };
    use crate::uff::params::AtomicParams;

    fn atomic_params(x1: f64, d1: f64) -> AtomicParams {
        AtomicParams {
            r1: 0.0,
            theta0: 0.0,
            x1,
            d1,
            zeta: 0.0,
            z1: 0.0,
            v1: 0.0,
            u1: 0.0,
            gmp_xi: 0.0,
            gmp_hardness: 0.0,
            gmp_radius: 0.0,
        }
    }

    fn energy_at(contribution: &VdwContrib, distance: f64) -> f64 {
        let coordinates = [0.0, 0.0, 0.0, distance, 0.0, 0.0];
        let mut distance_matrix = [-1.0; 3];
        let mut context = EvaluationContext::for_test(&coordinates, &mut distance_matrix, 2);
        contribution
            .get_energy(&mut context)
            .expect("valid two-point source distance")
    }

    fn gradient_at(
        contribution: &VdwContrib,
        distance: f64,
        initial_gradient: [f64; 6],
    ) -> [f64; 6] {
        let coordinates = [0.0, 0.0, 0.0, distance, 0.0, 0.0];
        let mut distance_matrix = [-1.0; 3];
        let mut context = EvaluationContext::for_test(&coordinates, &mut distance_matrix, 2);
        let mut gradient = initial_gradient;
        contribution
            .get_grad(&mut context, &mut gradient)
            .expect("valid two-point source distance");
        gradient
    }

    fn kernel_vdw_field<'a>(
        first_position: &'a mut [f64],
        second_position: &'a mut [f64],
        threshold_multiplier: f64,
    ) -> ForceField<'a> {
        let first_params = atomic_params(2.0, 3.0);
        let second_params = atomic_params(8.0, 12.0);
        let mut field = ForceField::new(3);
        field
            .positions_mut()
            .extend([first_position, second_position]);
        let contribution = VdwContrib::new_with_threshold_multiplier(
            field.positions(),
            0,
            1,
            &first_params,
            &second_params,
            threshold_multiplier,
        )
        .expect("source-valid pair");
        field.add_contribution(Box::new(contribution));
        field.initialize().expect("two-point kernel field");
        field
    }

    #[test]
    fn cf3d_u18_mixed_parameters_default_threshold_and_minimum_energy() {
        // The pinned calcNonbondedMinimum/Depth helpers use geometric means;
        // the source header defaults threshMultiplier to 10.0.
        let first_params = atomic_params(2.0, 3.0);
        let second_params = atomic_params(8.0, 12.0);
        assert_eq!(calc_nonbonded_minimum(&first_params, &second_params), 4.0);
        assert_eq!(calc_nonbonded_depth(&first_params, &second_params), 6.0);
        assert_eq!(calc_nonbonded_minimum(&second_params, &first_params), 4.0);
        assert_eq!(calc_nonbonded_depth(&second_params, &first_params), 6.0);

        let mut first_position = vec![0.0; 3];
        let mut second_position = vec![0.0; 3];
        let positions = [
            first_position.as_mut_slice(),
            second_position.as_mut_slice(),
        ];
        let contribution = VdwContrib::new(&positions, 0, 1, &first_params, &second_params)
            .expect("source-valid atom pair");

        assert_eq!(contribution.x_ij, 4.0);
        assert_eq!(contribution.well_depth, 6.0);
        assert_eq!(contribution.threshold, 40.0);
        // At the source-preferred distance xij, r=1 and the source well is -D.
        assert_eq!(energy_at(&contribution, 4.0), -6.0);
    }

    #[test]
    fn cf3d_u18_custom_cutoff_includes_equality_and_preserves_zero_branch() {
        let first_params = atomic_params(2.0, 3.0);
        let second_params = atomic_params(8.0, 12.0);
        let mut first_position = vec![0.0; 3];
        let mut second_position = vec![0.0; 3];
        let positions = [
            first_position.as_mut_slice(),
            second_position.as_mut_slice(),
        ];
        let contribution = VdwContrib::new_with_threshold_multiplier(
            &positions,
            0,
            1,
            &first_params,
            &second_params,
            1.0,
        )
        .expect("source-valid atom pair and explicit threshold");

        // Fixed source-formula value: r=2, r6=64, r12=4096, D=6.
        assert_eq!(energy_at(&contribution, 2.0), 23_808.0);
        // The strict `dist > d_thresh` check means equality evaluates.
        assert_eq!(energy_at(&contribution, 4.0), -6.0);
        assert_eq!(energy_at(&contribution, 5.0).to_bits(), 0.0_f64.to_bits());
        assert_eq!(energy_at(&contribution, 0.0).to_bits(), 0.0_f64.to_bits());

        // The source constructor does not validate an explicit multiplier.
        for threshold_multiplier in [0.0, -1.0] {
            let contribution = VdwContrib::new_with_threshold_multiplier(
                &positions,
                0,
                1,
                &first_params,
                &second_params,
                threshold_multiplier,
            )
            .expect("source constructor accepts numeric threshold values");
            assert_eq!(energy_at(&contribution, 2.0).to_bits(), 0.0_f64.to_bits());
        }
    }

    #[test]
    fn cf3d_u18_constructor_preserves_pair_range_check_order() {
        let first_params = atomic_params(2.0, 3.0);
        let second_params = atomic_params(8.0, 12.0);
        let mut first_position = vec![0.0; 3];
        let mut second_position = vec![0.0; 3];
        let positions = [
            first_position.as_mut_slice(),
            second_position.as_mut_slice(),
        ];

        assert_eq!(
            VdwContrib::new(&positions, 2, 2, &first_params, &second_params),
            Err(ForceFieldKernelError::BondIndexOutOfRange {
                argument: BondIndexArgument::First,
                index: 2,
                upper_bound: 2,
            })
        );
        assert_eq!(
            VdwContrib::new(&positions, 0, 2, &first_params, &second_params),
            Err(ForceFieldKernelError::BondIndexOutOfRange {
                argument: BondIndexArgument::Second,
                index: 2,
                upper_bound: 2,
            })
        );
    }

    #[test]
    fn cf3d_u18_preserves_source_nan_flow_without_numeric_guards() {
        let negative_x = atomic_params(-2.0, -3.0);
        let positive = atomic_params(1.0, 1.0);
        assert!(calc_nonbonded_minimum(&negative_x, &positive).is_nan());
        assert!(calc_nonbonded_depth(&negative_x, &positive).is_nan());

        let mut first_position = vec![0.0; 3];
        let mut second_position = vec![0.0; 3];
        let positions = [
            first_position.as_mut_slice(),
            second_position.as_mut_slice(),
        ];
        let contribution = VdwContrib::new(&positions, 0, 1, &negative_x, &positive)
            .expect("source constructor does not reject negative parameter products");
        assert!(energy_at(&contribution, 1.0).is_nan());

        let valid_first = atomic_params(2.0, 3.0);
        let valid_second = atomic_params(8.0, 12.0);
        let valid_contribution = VdwContrib::new(&positions, 0, 1, &valid_first, &valid_second)
            .expect("source-valid atom pair");
        assert!(energy_at(&valid_contribution, f64::NAN).is_nan());
    }

    #[test]
    fn cf3d_u19_gradient_matches_fixed_repulsive_attractive_and_cutoff_values() {
        // Fixed values follow pinned Nonbonded.cpp::getGrad and the r^7/r^13
        // source multiplication tree for xij=4 and well depth=6.
        let first_params = atomic_params(2.0, 3.0);
        let second_params = atomic_params(8.0, 12.0);
        let mut first_position = vec![0.0; 3];
        let mut second_position = vec![0.0; 3];
        let positions = [
            first_position.as_mut_slice(),
            second_position.as_mut_slice(),
        ];
        let cutoff = VdwContrib::new_with_threshold_multiplier(
            &positions,
            0,
            1,
            &first_params,
            &second_params,
            0.5,
        )
        .expect("source-valid pair with cutoff xij/2");

        // At distance 1, r=4, r7=16_384, r13=67_108_864. This is inside
        // the cutoff and exercises the repulsive analytic-gradient regime.
        assert_eq!(
            gradient_at(&cutoff, 1.0, [0.0; 6]),
            [1_207_664_640.0, 0.0, 0.0, -1_207_664_640.0, 0.0, 0.0,]
        );

        // The threshold is exactly 2.0. Strict greater-than leaves equality
        // in the formula; the nonzero result distinguishes it from a return.
        assert_eq!(
            gradient_at(&cutoff, 2.0, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),
            [145_153.0, 2.0, 3.0, -145_148.0, 5.0, 6.0]
        );
        assert_eq!(
            gradient_at(&cutoff, 2.5, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),
            [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        );

        let default_cutoff = VdwContrib::new(&positions, 0, 1, &first_params, &second_params)
            .expect("source default threshold");
        // At distance 8 the source attraction gives a fixed prefactor of
        // 567/4096; initial values prove endpoint writes are additive.
        assert_eq!(
            gradient_at(&default_cutoff, 8.0, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),
            [0.861_572_265_625, 2.0, 3.0, 4.138_427_734_375, 5.0, 6.0]
        );
    }

    #[test]
    fn cf3d_u19_zero_distance_order_additivity_and_nan_flow() {
        let first_params = atomic_params(2.0, 3.0);
        let second_params = atomic_params(8.0, 12.0);
        let mut first_position = vec![0.0; 3];
        let mut second_position = vec![0.0; 3];
        let positions = [
            first_position.as_mut_slice(),
            second_position.as_mut_slice(),
        ];

        let default_cutoff = VdwContrib::new(&positions, 0, 1, &first_params, &second_params)
            .expect("source default threshold");
        // Source's zero-distance branch adds +100/-100 on each of three axes.
        assert_eq!(
            gradient_at(&default_cutoff, 0.0, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),
            [101.0, 102.0, 103.0, -96.0, -95.0, -94.0]
        );

        // Cutoff is tested first: zero threshold enters the zero-distance
        // branch at distance 0, while a negative threshold returns unchanged.
        let zero_threshold = VdwContrib::new_with_threshold_multiplier(
            &positions,
            0,
            1,
            &first_params,
            &second_params,
            0.0,
        )
        .expect("source accepts zero multiplier");
        assert_eq!(
            gradient_at(&zero_threshold, 0.0, [0.0; 6]),
            [100.0, 100.0, 100.0, -100.0, -100.0, -100.0]
        );
        let negative_threshold = VdwContrib::new_with_threshold_multiplier(
            &positions,
            0,
            1,
            &first_params,
            &second_params,
            -1.0,
        )
        .expect("source accepts negative multiplier");
        assert_eq!(
            gradient_at(&negative_threshold, 0.0, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),
            [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        );

        let nan_gradient = gradient_at(&default_cutoff, f64::NAN, [0.0; 6]);
        assert!(nan_gradient.into_iter().all(f64::is_nan));
    }

    #[test]
    fn cf3d_bld_b08_kernel_dispatch_matches_vdw_threshold_and_copy() {
        let mut first_position = [0.0; 3];
        let mut second_position = [0.0; 3];
        let mut field = kernel_vdw_field(&mut first_position, &mut second_position, 2.0);

        let inside = [0.0, 0.0, 0.0, 4.0, 0.0, 0.0];
        assert_eq!(cf3d_bld_b05_calc_energy(&mut field, &inside), Ok(-6.0));
        let mut inside_gradient = [0.0; 6];
        cf3d_bld_b05_calc_grad(&mut field, &inside, &mut inside_gradient)
            .expect("kernel dispatches the VDW gradient");
        assert_eq!(inside_gradient, [0.0; 6]);

        // xij=4, threshold=8. The strict source comparison evaluates equality.
        let at_threshold = [0.0, 0.0, 0.0, 8.0, 0.0, 0.0];
        assert_eq!(
            cf3d_bld_b05_calc_energy(&mut field, &at_threshold),
            Ok(-762.0 / 4096.0)
        );
        let mut threshold_gradient = [0.0; 6];
        cf3d_bld_b05_calc_grad(&mut field, &at_threshold, &mut threshold_gradient)
            .expect("kernel dispatches equality gradient");
        assert_eq!(
            threshold_gradient,
            [-0.138_427_734_375, 0.0, 0.0, 0.138_427_734_375, 0.0, 0.0,]
        );

        let outside = [0.0, 0.0, 0.0, 9.0, 0.0, 0.0];
        assert_eq!(cf3d_bld_b05_calc_energy(&mut field, &outside), Ok(0.0));
        let mut outside_gradient = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        cf3d_bld_b05_calc_grad(&mut field, &outside, &mut outside_gradient)
            .expect("kernel dispatches the outside-cutoff gradient");
        assert_eq!(outside_gradient, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);

        let coincident = [0.0; 6];
        assert_eq!(cf3d_bld_b05_calc_energy(&mut field, &coincident), Ok(0.0));
        let mut coincident_gradient = [0.0; 6];
        cf3d_bld_b05_calc_grad(&mut field, &coincident, &mut coincident_gradient)
            .expect("kernel dispatches the coincident-point gradient");
        assert_eq!(
            coincident_gradient,
            [100.0, 100.0, 100.0, -100.0, -100.0, -100.0]
        );

        let mut copied = cf3d_bld_b05_copy_force_field(&field);
        let mut copied_first = [0.0; 3];
        let mut copied_second = [0.0; 3];
        copied
            .positions_mut()
            .extend([copied_first.as_mut_slice(), copied_second.as_mut_slice()]);
        copied
            .initialize()
            .expect("copied field accepts fresh borrows");
        assert_eq!(
            cf3d_bld_b05_calc_energy(&mut copied, &at_threshold),
            Ok(-762.0 / 4096.0)
        );
        let mut copied_gradient = [0.0; 6];
        cf3d_bld_b05_calc_grad(&mut copied, &at_threshold, &mut copied_gradient)
            .expect("copied contribution remains evaluable");
        assert_eq!(copied_gradient, threshold_gradient);
    }

    #[test]
    fn cf3d_bld_b08_preserves_typed_constructor_and_kernel_failures() {
        let first_params = atomic_params(2.0, 3.0);
        let second_params = atomic_params(8.0, 12.0);
        let mut first_position = [0.0; 3];
        let mut second_position = [0.0; 3];
        let positions = [
            first_position.as_mut_slice(),
            second_position.as_mut_slice(),
        ];

        assert_eq!(
            VdwContrib::new(&positions, 2, 2, &first_params, &second_params),
            Err(ForceFieldKernelError::BondIndexOutOfRange {
                argument: BondIndexArgument::First,
                index: 2,
                upper_bound: 2,
            })
        );
        assert_eq!(
            VdwContrib::new(&positions, 0, 2, &first_params, &second_params),
            Err(ForceFieldKernelError::BondIndexOutOfRange {
                argument: BondIndexArgument::Second,
                index: 2,
                upper_bound: 2,
            })
        );

        let mut uninitialized = ForceField::new(3);
        assert_eq!(
            cf3d_bld_b05_calc_energy(&mut uninitialized, &[0.0; 6]),
            Err(ForceFieldKernelError::NotInitialized)
        );
        assert_eq!(
            cf3d_bld_b05_calc_grad(&mut uninitialized, &[0.0; 6], &mut [0.0; 6]),
            Err(ForceFieldKernelError::NotInitialized)
        );
    }
}
