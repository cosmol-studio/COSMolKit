//! Source-backed RDKit UFF bond stretch contribution.

use crate::chemistry::forcefield::core::{ForceField, ForceFieldContrib};

use super::params::AtomicParams;
use super::utils::{UffUtilsError, calc_bond_force_constant, calc_bond_rest_length};

#[derive(Clone, Debug)]
pub struct BondStretchContrib {
    owner: *const ForceField,
    end1_idx: usize,
    end2_idx: usize,
    rest_len: f64,
    force_constant: f64,
}

impl BondStretchContrib {
    pub fn new(
        owner: &ForceField,
        idx1: usize,
        idx2: usize,
        bond_order: f64,
        end1_params: &AtomicParams,
        end2_params: &AtomicParams,
    ) -> Result<Self, UffUtilsError> {
        // RDKit✔️✔️: BondStretchContrib::BondStretchContrib(ForceField *owner, unsigned int idx1,
        // RDKit✔️✔️:                                        unsigned int idx2, double bondOrder,
        // RDKit✔️✔️:                                        const AtomicParams *end1Params,
        // RDKit✔️✔️:                                        const AtomicParams *end2Params) {
        // RDKit✔️✔️:   PRECONDITION(owner, "bad owner");
        // RDKit✔️✔️:   PRECONDITION(end1Params, "bad params pointer");
        // RDKit✔️✔️:   PRECONDITION(end2Params, "bad params pointer");
        // Rust references reproduce RDKit's non-null preconditions.
        // RDKit✔️✔️:   URANGE_CHECK(idx1, owner->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx2, owner->positions().size());
        assert!(idx1 < owner.positions().len());
        assert!(idx2 < owner.positions().len());

        // RDKit✔️✔️:   dp_forceField = owner;
        // RDKit✔️✔️:   d_end1Idx = idx1;
        // RDKit✔️✔️:   d_end2Idx = idx2;
        let owner_ptr = owner as *const ForceField;

        // RDKit✔️✔️:   d_restLen = Utils::calcBondRestLength(bondOrder, end1Params, end2Params);
        let rest_len = calc_bond_rest_length(bond_order, end1_params, end2_params)?;
        // RDKit✔️✔️:   d_forceConstant =
        // RDKit✔️✔️:       Utils::calcBondForceConstant(d_restLen, end1Params, end2Params);
        // RDKit✔️✔️: }
        let force_constant = calc_bond_force_constant(rest_len, end1_params, end2_params)?;

        Ok(Self {
            owner: owner_ptr,
            end1_idx: idx1,
            end2_idx: idx2,
            rest_len,
            force_constant,
        })
    }

    #[must_use]
    pub fn owner(&self) -> *const ForceField {
        self.owner
    }

    #[must_use]
    pub fn end1_idx(&self) -> usize {
        self.end1_idx
    }

    #[must_use]
    pub fn end2_idx(&self) -> usize {
        self.end2_idx
    }

    #[must_use]
    pub fn rest_len(&self) -> f64 {
        self.rest_len
    }

    #[must_use]
    pub fn force_constant(&self) -> f64 {
        self.force_constant
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
        // RDKit✔️✔️: double BondStretchContrib::getEnergy(double *pos) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        let force_field = self.force_field();

        // RDKit✔️✔️:   double distTerm =
        // RDKit✔️✔️:       dp_forceField->distance(d_end1Idx, d_end2Idx, pos) - d_restLen;
        let dist_term =
            force_field.distance(self.end1_idx, self.end2_idx, Some(pos)) - self.rest_len;
        // RDKit✔️✔️:   double res = 0.5 * d_forceConstant * distTerm * distTerm;
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        0.5 * self.force_constant * dist_term * dist_term
    }

    pub fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
        // RDKit✔️✔️: void BondStretchContrib::getGrad(double *pos, double *grad) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:   PRECONDITION(grad, "bad vector");
        let force_field = self.force_field();

        // RDKit✔️✔️:   double dist = dp_forceField->distance(d_end1Idx, d_end2Idx, pos);
        let dist = force_field.distance(self.end1_idx, self.end2_idx, Some(pos));
        // RDKit✔️✔️:   double preFactor = d_forceConstant * (dist - d_restLen);
        let pre_factor = self.force_constant * (dist - self.rest_len);

        // RDKit✔️✔️:   // std::cout << "\tDist("<<d_end1Idx<<","<<d_end2Idx<<") " << dist <<
        // RDKit✔️✔️:   // std::endl;
        // RDKit✔️✔️:   double *end1Coords = &(pos[3 * d_end1Idx]);
        // RDKit✔️✔️:   double *end2Coords = &(pos[3 * d_end2Idx]);
        let end1_offset = 3 * self.end1_idx;
        let end2_offset = 3 * self.end2_idx;
        // RDKit✔️✔️:   for (int i = 0; i < 3; i++) {
        for i in 0..3 {
            // RDKit✔️✔️:     double dGrad;
            let d_grad;
            // RDKit✔️✔️:     if (dist > 0.0) {
            if dist > 0.0 {
                // RDKit✔️✔️:       dGrad = preFactor * (end1Coords[i] - end2Coords[i]) / dist;
                d_grad = pre_factor * (pos[end1_offset + i] - pos[end2_offset + i]) / dist;
                // RDKit✔️✔️:     } else {
            } else {
                // RDKit✔️✔️:       // move a small amount in an arbitrary direction
                // RDKit✔️✔️:       dGrad = d_forceConstant * .01;
                d_grad = self.force_constant * 0.01;
                // RDKit✔️✔️:     }
            }
            // RDKit✔️✔️:     grad[3 * d_end1Idx + i] += dGrad;
            // RDKit✔️✔️:     grad[3 * d_end2Idx + i] -= dGrad;
            grad[end1_offset + i] += d_grad;
            grad[end2_offset + i] -= d_grad;
        }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
    }
}

impl ForceFieldContrib for BondStretchContrib {
    fn copy(&self) -> Box<dyn ForceFieldContrib> {
        // RDKit✔️✔️: BondStretchContrib *copy() const override {
        // RDKit✔️✔️:   return new BondStretchContrib(*this);
        // RDKit✔️✔️: }
        Box::new(self.clone())
    }

    fn set_force_field(&mut self, owner: *const ForceField) {
        self.owner = owner;
    }

    fn get_energy(&self, pos: &[f64]) -> f64 {
        BondStretchContrib::get_energy(self, pos)
    }

    fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
        BondStretchContrib::get_grad(self, pos, grad);
    }
}

#[cfg(test)]
mod tests {
    use crate::chemistry::forcefield::core::ForceFieldVec3;

    use super::*;
    use crate::chemistry::forcefield::uff::utils::{
        calc_bond_force_constant, calc_bond_rest_length,
    };

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

    #[test]
    fn uff_bondstretchcontrib_constructor_stores_indices_owner_and_calculated_terms() {
        let ff = force_field();
        let end1 = atomic_params(0.757, 5.343, 1.912);
        let end2 = atomic_params(0.658, 8.741, 2.3);

        let contrib =
            BondStretchContrib::new(&ff, 0, 1, 1.5, &end1, &end2).expect("valid constructor");
        let expected_rest_len = calc_bond_rest_length(1.5, &end1, &end2).expect("valid bond order");
        let expected_force_constant =
            calc_bond_force_constant(expected_rest_len, &end1, &end2).expect("valid rest length");

        assert_eq!(contrib.owner(), &ff as *const ForceField);
        assert_eq!(contrib.end1_idx(), 0);
        assert_eq!(contrib.end2_idx(), 1);
        assert!((contrib.rest_len() - expected_rest_len).abs() < EPS);
        assert!((contrib.force_constant() - expected_force_constant).abs() < EPS);
    }

    #[test]
    fn uff_bondstretchcontrib_constructor_returns_invalid_bond_order_error() {
        let ff = force_field();
        let end1 = atomic_params(0.757, 5.343, 1.912);
        let end2 = atomic_params(0.658, 8.741, 2.3);

        assert!(matches!(
            BondStretchContrib::new(&ff, 0, 1, 0.0, &end1, &end2),
            Err(UffUtilsError::InvalidBondOrder { bond_order: 0.0 })
        ));
    }

    #[test]
    #[should_panic]
    fn uff_bondstretchcontrib_constructor_rejects_first_index_out_of_range() {
        let ff = force_field();
        let end1 = atomic_params(0.757, 5.343, 1.912);
        let end2 = atomic_params(0.658, 8.741, 2.3);

        let _ = BondStretchContrib::new(&ff, 2, 1, 1.0, &end1, &end2);
    }

    #[test]
    #[should_panic]
    fn uff_bondstretchcontrib_constructor_rejects_second_index_out_of_range() {
        let ff = force_field();
        let end1 = atomic_params(0.757, 5.343, 1.912);
        let end2 = atomic_params(0.658, 8.741, 2.3);

        let _ = BondStretchContrib::new(&ff, 0, 2, 1.0, &end1, &end2);
    }

    #[test]
    fn uff_bondstretchcontrib_get_energy_returns_zero_at_rest_length() {
        let ff = initialized_force_field();
        let end1 = atomic_params(0.757, 5.343, 1.912);
        let end2 = atomic_params(0.658, 8.741, 2.3);
        let contrib =
            BondStretchContrib::new(&ff, 0, 1, 1.0, &end1, &end2).expect("valid constructor");
        let pos = [0.0, 0.0, 0.0, contrib.rest_len(), 0.0, 0.0];

        assert!(contrib.get_energy(&pos).abs() < EPS);
    }

    #[test]
    fn uff_bondstretchcontrib_get_energy_uses_quadratic_distance_term() {
        let ff = initialized_force_field();
        let end1 = atomic_params(0.757, 5.343, 1.912);
        let end2 = atomic_params(0.658, 8.741, 2.3);
        let contrib =
            BondStretchContrib::new(&ff, 0, 1, 1.0, &end1, &end2).expect("valid constructor");
        let dist = contrib.rest_len() + 0.25;
        let pos = [0.0, 0.0, 0.0, dist, 0.0, 0.0];
        let expected = 0.5 * contrib.force_constant() * 0.25_f64 * 0.25_f64;

        assert!((contrib.get_energy(&pos) - expected).abs() < EPS);
    }

    #[test]
    fn uff_bondstretchcontrib_get_grad_is_zero_at_rest_length() {
        let ff = initialized_force_field();
        let end1 = atomic_params(0.757, 5.343, 1.912);
        let end2 = atomic_params(0.658, 8.741, 2.3);
        let contrib =
            BondStretchContrib::new(&ff, 0, 1, 1.0, &end1, &end2).expect("valid constructor");
        let pos = [0.0, 0.0, 0.0, contrib.rest_len(), 0.0, 0.0];
        let mut grad = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let original = grad;

        contrib.get_grad(&pos, &mut grad);

        assert_eq!(grad, original);
    }

    #[test]
    fn uff_bondstretchcontrib_get_grad_accumulates_normal_distance_branch() {
        let ff = initialized_force_field();
        let end1 = atomic_params(0.757, 5.343, 1.912);
        let end2 = atomic_params(0.658, 8.741, 2.3);
        let contrib =
            BondStretchContrib::new(&ff, 0, 1, 1.0, &end1, &end2).expect("valid constructor");
        let offset = 0.25;
        let dist = contrib.rest_len() + offset;
        let pos = [0.0, 0.0, 0.0, dist, 0.0, 0.0];
        let mut grad = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0];
        let d_grad = -contrib.force_constant() * offset;

        contrib.get_grad(&pos, &mut grad);

        assert!((grad[0] - (10.0 + d_grad)).abs() < EPS);
        assert_eq!(grad[1], 20.0);
        assert_eq!(grad[2], 30.0);
        assert!((grad[3] - (40.0 - d_grad)).abs() < EPS);
        assert_eq!(grad[4], 50.0);
        assert_eq!(grad[5], 60.0);
    }

    #[test]
    fn uff_bondstretchcontrib_get_grad_accumulates_zero_distance_branch() {
        let ff = initialized_force_field();
        let end1 = atomic_params(0.757, 5.343, 1.912);
        let end2 = atomic_params(0.658, 8.741, 2.3);
        let contrib =
            BondStretchContrib::new(&ff, 0, 1, 1.0, &end1, &end2).expect("valid constructor");
        let pos = [0.0; 6];
        let mut grad = [0.0; 6];
        let d_grad = contrib.force_constant() * 0.01;

        contrib.get_grad(&pos, &mut grad);

        for i in 0..3 {
            assert!((grad[i] - d_grad).abs() < EPS);
            assert!((grad[3 + i] + d_grad).abs() < EPS);
        }
    }
}
