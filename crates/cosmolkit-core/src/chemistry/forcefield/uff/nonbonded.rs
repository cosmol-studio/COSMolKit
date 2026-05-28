//! Source-backed RDKit UFF van der Waals nonbonded contribution.

use crate::chemistry::forcefield::core::{ForceField, ForceFieldContrib};

use super::params::AtomicParams;
use super::utils::{calc_nonbonded_depth, calc_nonbonded_minimum};

#[derive(Clone, Debug)]
pub struct VdwContrib {
    owner: *const ForceField,
    at1_idx: usize,
    at2_idx: usize,
    x_ij: f64,
    well_depth: f64,
    thresh: f64,
}

impl VdwContrib {
    pub fn new(
        owner: &ForceField,
        idx1: usize,
        idx2: usize,
        at1_params: &AtomicParams,
        at2_params: &AtomicParams,
        thresh_multiplier: f64,
    ) -> Self {
        // RDKit✔️✔️: vdWContrib::vdWContrib(ForceField *owner, unsigned int idx1, unsigned int idx2,
        // RDKit✔️✔️:                        const AtomicParams *at1Params,
        // RDKit✔️✔️:                        const AtomicParams *at2Params, double threshMultiplier) {
        // RDKit✔️✔️:   PRECONDITION(owner, "bad owner");
        // RDKit✔️✔️:   PRECONDITION(at1Params, "bad params pointer");
        // RDKit✔️✔️:   PRECONDITION(at2Params, "bad params pointer");
        // Rust references reproduce RDKit's non-null preconditions.
        // RDKit✔️✔️:   URANGE_CHECK(idx1, owner->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx2, owner->positions().size());
        assert!(idx1 < owner.positions().len());
        assert!(idx2 < owner.positions().len());

        // RDKit✔️✔️:   dp_forceField = owner;
        // RDKit✔️✔️:   d_at1Idx = idx1;
        // RDKit✔️✔️:   d_at2Idx = idx2;
        let owner_ptr = owner as *const ForceField;

        // RDKit✔️✔️:   // UFF uses the geometric mean of the vdW parameters:
        // RDKit✔️✔️:   d_xij = Utils::calcNonbondedMinimum(at1Params, at2Params);
        let x_ij = calc_nonbonded_minimum(at1_params, at2_params);
        // RDKit✔️✔️:   d_wellDepth = Utils::calcNonbondedDepth(at1Params, at2Params);
        let well_depth = calc_nonbonded_depth(at1_params, at2_params);
        // RDKit✔️✔️:   d_thresh = threshMultiplier * d_xij;
        let thresh = thresh_multiplier * x_ij;

        // RDKit✔️✔️:   // std::cerr << "  non-bonded: " << idx1 << "-" << idx2 << " " << d_xij << " "
        // RDKit✔️✔️:   // << d_wellDepth << " " << d_thresh << std::endl;
        // RDKit✔️✔️: }
        Self {
            owner: owner_ptr,
            at1_idx: idx1,
            at2_idx: idx2,
            x_ij,
            well_depth,
            thresh,
        }
    }

    #[must_use]
    pub fn owner(&self) -> *const ForceField {
        self.owner
    }

    #[must_use]
    pub fn at1_idx(&self) -> usize {
        self.at1_idx
    }

    #[must_use]
    pub fn at2_idx(&self) -> usize {
        self.at2_idx
    }

    #[must_use]
    pub fn x_ij(&self) -> f64 {
        self.x_ij
    }

    #[must_use]
    pub fn well_depth(&self) -> f64 {
        self.well_depth
    }

    #[must_use]
    pub fn thresh(&self) -> f64 {
        self.thresh
    }

    fn force_field(&self) -> &ForceField {
        assert!(!self.owner.is_null(), "no owner");
        // SAFETY: UFF contribs follow the same owner-pointer model as
        // ForceFieldContrib objects in core.rs. Constructors and ForceField::add_contrib
        // store a live owning ForceField pointer before energy/gradient evaluation.
        unsafe { &*self.owner }
    }

    #[must_use]
    pub fn get_energy(&self, pos: &[f64]) -> f64 {
        // RDKit✔️✔️: double vdWContrib::getEnergy(double *pos) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        let force_field = self.force_field();

        // RDKit✔️✔️:   double dist = dp_forceField->distance(d_at1Idx, d_at2Idx, pos);
        let dist = force_field.distance_const(self.at1_idx, self.at2_idx, Some(pos));
        // RDKit✔️✔️:   if (dist > d_thresh || dist <= 0.0) {
        if dist > self.thresh || dist <= 0.0 {
            // RDKit✔️✔️:     return 0.0;
            // RDKit✔️✔️:   }
            return 0.0;
        }

        // RDKit✔️✔️:   double r = d_xij / dist;
        let r = self.x_ij / dist;
        // RDKit✔️✔️: template <unsigned n>
        // RDKit✔️✔️: inline double int_pow(double x) {
        // RDKit✔️✔️:   double half = int_pow<n / 2>(x);
        // RDKit✔️✔️:   if (n % 2 == 0) {  // even
        // RDKit✔️✔️:     return half * half;
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     return half * half * x;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // RDKit✔️✔️: template <>
        // RDKit✔️✔️: inline double int_pow<1>(double x) {
        // RDKit✔️✔️:   return x;  // this does a series of muls
        // RDKit✔️✔️: }
        // RDKit✔️✔️:   double r6 = int_pow<6>(r);
        let r3 = r * r * r;
        let r6 = r3 * r3;
        // RDKit✔️✔️:   double r12 = r6 * r6;
        let r12 = r6 * r6;
        // RDKit✔️✔️:   double res = d_wellDepth * (r12 - 2.0 * r6);
        let res = self.well_depth * (r12 - 2.0 * r6);
        // RDKit✔️✔️:   // if(d_at1Idx==12 && d_at2Idx==21 ) std::cerr << "     >: " << d_at1Idx <<
        // RDKit✔️✔️:   // "-" << d_at2Idx << " " << r << " = " << res << std::endl;
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        res
    }

    pub fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
        // RDKit✔️✔️: void vdWContrib::getGrad(double *pos, double *grad) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:   PRECONDITION(grad, "bad vector");
        let force_field = self.force_field();

        // RDKit✔️✔️:   double dist = dp_forceField->distance(d_at1Idx, d_at2Idx, pos);
        let dist = force_field.distance_const(self.at1_idx, self.at2_idx, Some(pos));
        // RDKit✔️✔️:   if (dist > d_thresh) {
        if dist > self.thresh {
            // RDKit✔️✔️:     return;
            // RDKit✔️✔️:   }
            return;
        }

        // RDKit✔️✔️:   if (dist <= 0) {
        if dist <= 0.0 {
            // RDKit✔️✔️:     for (int i = 0; i < 3; i++) {
            for i in 0..3 {
                // RDKit✔️✔️:       // move in an arbitrary direction
                // RDKit✔️✔️:       double dGrad = 100.0;
                let d_grad = 100.0;
                // RDKit✔️✔️:       grad[3 * d_at1Idx + i] += dGrad;
                grad[3 * self.at1_idx + i] += d_grad;
                // RDKit✔️✔️:       grad[3 * d_at2Idx + i] -= dGrad;
                grad[3 * self.at2_idx + i] -= d_grad;
                // RDKit✔️✔️:     }
            }
            // RDKit✔️✔️:     return;
            // RDKit✔️✔️:   }
            return;
        }

        // RDKit✔️✔️:   double r = d_xij / dist;
        let r = self.x_ij / dist;
        // RDKit✔️✔️: template <unsigned n>
        // RDKit✔️✔️: inline double int_pow(double x) {
        // RDKit✔️✔️:   double half = int_pow<n / 2>(x);
        // RDKit✔️✔️:   if (n % 2 == 0) {  // even
        // RDKit✔️✔️:     return half * half;
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     return half * half * x;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // RDKit✔️✔️: template <>
        // RDKit✔️✔️: inline double int_pow<1>(double x) {
        // RDKit✔️✔️:   return x;  // this does a series of muls
        // RDKit✔️✔️: }
        // RDKit✔️✔️:   double r7 = int_pow<7>(r);
        let r3 = r * r * r;
        let r7 = r3 * r3 * r;
        // RDKit✔️✔️:   double r13 = int_pow<13>(r);
        let r6 = r3 * r3;
        let r13 = r6 * r6 * r;
        // RDKit✔️✔️:   double preFactor = 12. * d_wellDepth / d_xij * (r7 - r13);
        let pre_factor = 12.0 * self.well_depth / self.x_ij * (r7 - r13);

        // RDKit✔️✔️:   double *at1Coords = &(pos[3 * d_at1Idx]);
        let at1_base = 3 * self.at1_idx;
        // RDKit✔️✔️:   double *at2Coords = &(pos[3 * d_at2Idx]);
        let at2_base = 3 * self.at2_idx;
        // RDKit✔️✔️:   for (int i = 0; i < 3; i++) {
        for i in 0..3 {
            // RDKit✔️✔️:     double dGrad = preFactor * (at1Coords[i] - at2Coords[i]) / dist;
            let d_grad = pre_factor * (pos[at1_base + i] - pos[at2_base + i]) / dist;
            // RDKit✔️✔️:     grad[3 * d_at1Idx + i] += dGrad;
            grad[at1_base + i] += d_grad;
            // RDKit✔️✔️:     grad[3 * d_at2Idx + i] -= dGrad;
            grad[at2_base + i] -= d_grad;
            // RDKit✔️✔️:   }
        }
        // RDKit✔️✔️: }
    }
}

impl ForceFieldContrib for VdwContrib {
    fn copy(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }

    fn set_force_field(&mut self, ff: *const ForceField) {
        self.owner = ff;
    }

    fn get_energy(&self, pos: &[f64]) -> f64 {
        Self::get_energy(self, pos)
    }

    fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
        Self::get_grad(self, pos, grad);
    }
}

#[cfg(test)]
mod tests {
    use crate::chemistry::forcefield::core::ForceFieldVec3;

    use super::*;

    const EPS: f64 = 1.0e-12;

    fn force_field() -> ForceField {
        let mut ff = ForceField::new(3);
        ff.positions_mut().push(ForceFieldVec3::new(0.0, 0.0, 0.0));
        ff.positions_mut().push(ForceFieldVec3::new(1.0, 0.0, 0.0));
        ff
    }

    fn initialized_force_field() -> ForceField {
        let mut ff = force_field();
        ff.initialize();
        ff
    }

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

    #[test]
    fn uff_vdwcontrib_constructor_stores_indices_owner_and_calculated_terms() {
        let ff = force_field();
        let at1 = atomic_params(3.851, 0.105);
        let at2 = atomic_params(3.5, 0.06);

        let contrib = VdwContrib::new(&ff, 0, 1, &at1, &at2, 10.0);
        let expected_x_ij = (at1.x1 * at2.x1).sqrt();
        let expected_well_depth = (at1.d1 * at2.d1).sqrt();

        assert_eq!(contrib.owner(), &ff as *const ForceField);
        assert_eq!(contrib.at1_idx(), 0);
        assert_eq!(contrib.at2_idx(), 1);
        assert!((contrib.x_ij() - expected_x_ij).abs() < EPS);
        assert!((contrib.well_depth() - expected_well_depth).abs() < EPS);
        assert!((contrib.thresh() - 10.0 * expected_x_ij).abs() < EPS);
    }

    #[test]
    fn uff_vdwcontrib_constructor_uses_supplied_threshold_multiplier() {
        let ff = force_field();
        let at1 = atomic_params(4.0, 0.25);
        let at2 = atomic_params(9.0, 0.36);

        let contrib = VdwContrib::new(&ff, 1, 0, &at1, &at2, 2.5);

        assert!((contrib.x_ij() - 6.0).abs() < EPS);
        assert!((contrib.well_depth() - 0.3).abs() < EPS);
        assert!((contrib.thresh() - 15.0).abs() < EPS);
    }

    #[test]
    fn uff_vdwcontrib_constructor_preserves_nan_from_nonbonded_minimum() {
        let ff = force_field();
        let at1 = atomic_params(-1.0, 0.25);
        let at2 = atomic_params(9.0, 0.36);

        let contrib = VdwContrib::new(&ff, 0, 1, &at1, &at2, 10.0);

        assert!(contrib.x_ij().is_nan());
        assert!(contrib.thresh().is_nan());
        assert!((contrib.well_depth() - 0.3).abs() < EPS);
    }

    #[test]
    fn uff_vdwcontrib_constructor_preserves_nan_from_nonbonded_depth() {
        let ff = force_field();
        let at1 = atomic_params(4.0, -0.25);
        let at2 = atomic_params(9.0, 0.36);

        let contrib = VdwContrib::new(&ff, 0, 1, &at1, &at2, 10.0);

        assert!((contrib.x_ij() - 6.0).abs() < EPS);
        assert!(contrib.well_depth().is_nan());
        assert!((contrib.thresh() - 60.0).abs() < EPS);
    }

    #[test]
    #[should_panic]
    fn uff_vdwcontrib_constructor_rejects_first_index_out_of_range() {
        let ff = force_field();
        let at1 = atomic_params(3.851, 0.105);
        let at2 = atomic_params(3.5, 0.06);

        let _ = VdwContrib::new(&ff, 2, 1, &at1, &at2, 10.0);
    }

    #[test]
    #[should_panic]
    fn uff_vdwcontrib_constructor_rejects_second_index_out_of_range() {
        let ff = force_field();
        let at1 = atomic_params(3.851, 0.105);
        let at2 = atomic_params(3.5, 0.06);

        let _ = VdwContrib::new(&ff, 0, 2, &at1, &at2, 10.0);
    }

    #[test]
    fn uff_vdwcontrib_get_energy_returns_zero_past_threshold() {
        let ff = initialized_force_field();
        let at1 = atomic_params(4.0, 0.25);
        let at2 = atomic_params(9.0, 0.36);
        let contrib = VdwContrib::new(&ff, 0, 1, &at1, &at2, 2.0);
        let pos = [0.0, 0.0, 0.0, 12.1, 0.0, 0.0];

        assert_eq!(contrib.get_energy(&pos), 0.0);
    }

    #[test]
    fn uff_vdwcontrib_get_energy_calculates_at_threshold_boundary() {
        let ff = initialized_force_field();
        let at1 = atomic_params(4.0, 0.25);
        let at2 = atomic_params(9.0, 0.36);
        let contrib = VdwContrib::new(&ff, 0, 1, &at1, &at2, 2.0);
        let pos = [0.0, 0.0, 0.0, 12.0, 0.0, 0.0];
        let r = contrib.x_ij() / 12.0;
        let r6 = r.powi(6);
        let expected = contrib.well_depth() * (r6 * r6 - 2.0 * r6);

        assert!((contrib.get_energy(&pos) - expected).abs() < EPS);
        assert!(contrib.get_energy(&pos) < 0.0);
    }

    #[test]
    fn uff_vdwcontrib_get_energy_returns_zero_at_zero_distance() {
        let ff = initialized_force_field();
        let at1 = atomic_params(4.0, 0.25);
        let at2 = atomic_params(9.0, 0.36);
        let contrib = VdwContrib::new(&ff, 0, 1, &at1, &at2, 2.0);
        let pos = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0];

        assert_eq!(contrib.get_energy(&pos), 0.0);
    }

    #[test]
    fn uff_vdwcontrib_get_energy_uses_lennard_jones_formula() {
        let ff = initialized_force_field();
        let at1 = atomic_params(4.0, 0.25);
        let at2 = atomic_params(9.0, 0.36);
        let contrib = VdwContrib::new(&ff, 0, 1, &at1, &at2, 10.0);
        let pos = [0.0, 0.0, 0.0, 3.0, 0.0, 0.0];
        let r = 2.0;
        let r6 = r * r * r * r * r * r;
        let expected = 0.3 * (r6 * r6 - 2.0 * r6);

        assert!((contrib.get_energy(&pos) - expected).abs() < EPS);
    }

    #[test]
    fn uff_vdwcontrib_get_energy_preserves_nan_well_depth() {
        let ff = initialized_force_field();
        let at1 = atomic_params(4.0, -0.25);
        let at2 = atomic_params(9.0, 0.36);
        let contrib = VdwContrib::new(&ff, 0, 1, &at1, &at2, 10.0);
        let pos = [0.0, 0.0, 0.0, 3.0, 0.0, 0.0];

        assert!(contrib.get_energy(&pos).is_nan());
    }

    #[test]
    fn uff_vdwcontrib_get_grad_leaves_gradient_unchanged_past_threshold() {
        let ff = initialized_force_field();
        let at1 = atomic_params(4.0, 0.25);
        let at2 = atomic_params(9.0, 0.36);
        let contrib = VdwContrib::new(&ff, 0, 1, &at1, &at2, 2.0);
        let pos = [0.0, 0.0, 0.0, 12.1, 0.0, 0.0];
        let mut grad = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let expected = grad;

        contrib.get_grad(&pos, &mut grad);

        assert_eq!(grad, expected);
    }

    #[test]
    fn uff_vdwcontrib_get_grad_uses_fixed_singularity_gradient_at_zero_distance() {
        let ff = initialized_force_field();
        let at1 = atomic_params(4.0, 0.25);
        let at2 = atomic_params(9.0, 0.36);
        let contrib = VdwContrib::new(&ff, 0, 1, &at1, &at2, 2.0);
        let pos = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
        let mut grad = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];

        contrib.get_grad(&pos, &mut grad);

        assert_eq!(grad, [101.0, 102.0, 103.0, -96.0, -95.0, -94.0]);
    }

    #[test]
    fn uff_vdwcontrib_get_grad_calculates_at_threshold_boundary() {
        let ff = initialized_force_field();
        let at1 = atomic_params(4.0, 0.25);
        let at2 = atomic_params(9.0, 0.36);
        let contrib = VdwContrib::new(&ff, 0, 1, &at1, &at2, 2.0);
        let pos = [0.0, 0.0, 0.0, 12.0, 0.0, 0.0];
        let mut grad = [0.0; 6];
        let r = contrib.x_ij() / 12.0;
        let r3 = r * r * r;
        let r7 = r3 * r3 * r;
        let r6 = r3 * r3;
        let r13 = r6 * r6 * r;
        let pre_factor = 12.0 * contrib.well_depth() / contrib.x_ij() * (r7 - r13);
        let expected_x = -pre_factor;

        contrib.get_grad(&pos, &mut grad);

        assert!((grad[0] - expected_x).abs() < EPS);
        assert_eq!(grad[1], 0.0);
        assert_eq!(grad[2], 0.0);
        assert!((grad[3] + expected_x).abs() < EPS);
        assert_eq!(grad[4], 0.0);
        assert_eq!(grad[5], 0.0);
    }

    #[test]
    fn uff_vdwcontrib_get_grad_uses_vdw_gradient_formula() {
        let ff = initialized_force_field();
        let at1 = atomic_params(4.0, 0.25);
        let at2 = atomic_params(9.0, 0.36);
        let contrib = VdwContrib::new(&ff, 0, 1, &at1, &at2, 10.0);
        let pos = [0.0, 0.0, 0.0, 3.0, 0.0, 0.0];
        let mut grad = [0.0; 6];
        let r = 2.0;
        let r3 = r * r * r;
        let r7 = r3 * r3 * r;
        let r6 = r3 * r3;
        let r13 = r6 * r6 * r;
        let pre_factor = 12.0 * 0.3 / 6.0 * (r7 - r13);
        let expected_x = -pre_factor;

        contrib.get_grad(&pos, &mut grad);

        assert!((grad[0] - expected_x).abs() < EPS);
        assert_eq!(grad[1], 0.0);
        assert_eq!(grad[2], 0.0);
        assert!((grad[3] + expected_x).abs() < EPS);
        assert_eq!(grad[4], 0.0);
        assert_eq!(grad[5], 0.0);
    }

    #[test]
    fn uff_vdwcontrib_get_grad_preserves_nan_well_depth() {
        let ff = initialized_force_field();
        let at1 = atomic_params(4.0, -0.25);
        let at2 = atomic_params(9.0, 0.36);
        let contrib = VdwContrib::new(&ff, 0, 1, &at1, &at2, 10.0);
        let pos = [0.0, 0.0, 0.0, 3.0, 0.0, 0.0];
        let mut grad = [0.0; 6];

        contrib.get_grad(&pos, &mut grad);

        assert!(grad[0].is_nan());
        assert!(grad[1].is_nan());
        assert!(grad[2].is_nan());
        assert!(grad[3].is_nan());
        assert!(grad[4].is_nan());
        assert!(grad[5].is_nan());
    }
}
