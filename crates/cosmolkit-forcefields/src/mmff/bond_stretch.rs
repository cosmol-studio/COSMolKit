//! Packed MMFF bond-stretch contribution; no molecule or builder ownership.
//!
//! Source: RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8,
//! Code/ForceField/MMFF/BondStretch.{cpp,h} (BSD license).

use super::{bonded::calc_bond_stretch_energy, params::MmffBond};
use crate::kernel::{
    BondIndexArgument, EvaluationContext, ForceField, ForceFieldContribution, ForceFieldKernelError,
};

#[derive(Clone, Debug)]
pub(super) struct BondStretchContrib {
    at1_idxs: Vec<u32>,
    at2_idxs: Vec<u32>,
    r0: Vec<f64>,
    kb: Vec<f64>,
}

pub(super) fn source_term_indices(term_count: usize) -> std::ops::Range<i32> {
    // RDKit✔️✔️:   const int numTerms = d_at1Idxs.size();
    // Behavior: BondStretch.cpp uses this conversion in both evaluations.
    // Pinned CMake selects C++20; ROOT's native probe confirms sizeof(int)=4.
    // Narrow modulo 2^32 to signed32, then iterate from zero only while less
    // than that signed end. Negative and zero ends produce no indices;
    // lengths wrapping to a positive end evaluate that prefix in order.
    // Complexity: O(1) conversion/range construction, no allocation or scan.
    0..(term_count as i32)
}

impl BondStretchContrib {
    pub(super) fn new(_owner: &ForceField<'_>) -> Self {
        // RDKit✔️✔️: BondStretchContrib::BondStretchContrib(ForceField *owner) {
        // RDKit✔️✔️:   PRECONDITION(owner, "bad owner");
        // RDKit✔️✔️:
        // RDKit✔️✔️:
        // RDKit✔️✔️:   dp_forceField = owner;
        // RDKit✔️✔️:
        // RDKit✔️✔️: }
        // Behavior: non-null owner is a borrow; the existing contribution
        // trait supplies owner-derived evaluation context instead of storing
        // a pointer. The ownerless C++ default constructor is not evaluable.
        // Complexity: four empty arrays, O(1), no allocation or owner clone.
        Self {
            at1_idxs: Vec::new(),
            at2_idxs: Vec::new(),
            r0: Vec::new(),
            kb: Vec::new(),
        }
    }

    pub(super) fn add_term(
        &mut self,
        positions: &[&mut [f64]],
        idx1: u32,
        idx2: u32,
        params: &MmffBond,
    ) -> Result<(), ForceFieldKernelError> {
        // RDKit✔️✔️: void BondStretchContrib::addTerm(const unsigned int idx1,
        // RDKit✔️✔️:                                  const unsigned int idx2,
        // RDKit✔️✔️:                                  const ForceFields::MMFF::MMFFBond *mmffBondParams) {
        // RDKit✔️✔️:   URANGE_CHECK(idx1, dp_forceField->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx2, dp_forceField->positions().size());
        // RDKit✔️✔️:   PRECONDITION(mmffBondParams, "bond parameters not found");
        // RDKit✔️✔️:
        // RDKit✔️✔️:   d_at1Idxs.push_back(idx1);
        // RDKit✔️✔️:   d_at2Idxs.push_back(idx2);
        // RDKit✔️✔️:   d_r0.push_back(mmffBondParams->r0);
        // RDKit✔️✔️:   d_kb.push_back(mmffBondParams->kb);
        // RDKit✔️✔️: }
        // Behavior: validate against the current owner's position view, in
        // source order, before appending. The parameter borrow cannot be null.
        // Copy parameter values; retain repetitions, orientation and IEEE data.
        // Complexity: two O(1) bounds checks and four amortized O(1) appends,
        // with the same separate-array layout and growth as the source.
        for (argument, index) in [
            (BondIndexArgument::First, idx1),
            (BondIndexArgument::Second, idx2),
        ] {
            if index as usize >= positions.len() {
                return Err(ForceFieldKernelError::BondIndexOutOfRange {
                    argument,
                    index,
                    upper_bound: positions.len(),
                });
            }
        }
        self.at1_idxs.push(idx1);
        self.at2_idxs.push(idx2);
        self.r0.push(params.r0);
        self.kb.push(params.kb);
        Ok(())
    }
}

impl ForceFieldContribution for BondStretchContrib {
    fn get_energy(
        &self,
        context: &mut EvaluationContext<'_>,
    ) -> Result<f64, ForceFieldKernelError> {
        // RDKit✔️✔️: double BondStretchContrib::getEnergy(double *pos) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:
        // RDKit✔️✔️:   const int numTerms = d_at1Idxs.size();
        // RDKit✔️✔️:   double energySum = 0.0;
        // RDKit✔️✔️:   for (int i =0; i < numTerms; i++) {
        // RDKit✔️✔️:     energySum += Utils::calcBondStretchEnergy(
        // RDKit✔️✔️:         d_r0[i], d_kb[i], dp_forceField->distance(d_at1Idxs[i], d_at2Idxs[i], pos));
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return energySum;
        // RDKit✔️✔️: }
        // Behavior: owner/vector borrows are non-null. Preserve source signed32
        // term-count conversion/order, distance cache/error semantics and the
        // single scalar owner. Nonpositive converted counts leave the sum zero.
        // Complexity: O(max(0, signed32(terms))), O(1) temporary space, no allocations,
        // no scans, lookup tables or coordinate/owner clones.
        let mut energy_sum = 0.0;
        for i in source_term_indices(self.at1_idxs.len()) {
            let i = i as usize;
            energy_sum += calc_bond_stretch_energy(
                self.r0[i],
                self.kb[i],
                context.distance(self.at1_idxs[i], self.at2_idxs[i])?,
            );
        }
        Ok(energy_sum)
    }

    fn get_grad(
        &self,
        context: &mut EvaluationContext<'_>,
        gradient: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        // RDKit✔️✔️: void BondStretchContrib::getGrad(double *pos, double *grad) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:   PRECONDITION(grad, "bad vector");
        // RDKit✔️✔️:
        // RDKit✔️✔️:   const int numTerms = d_at1Idxs.size();
        // RDKit✔️✔️:   constexpr double cs = -2.0;
        // RDKit✔️✔️:   constexpr double c1 = MDYNE_A_TO_KCAL_MOL;
        // RDKit✔️✔️:   constexpr double c3 = 7.0 / 12.0;
        // RDKit✔️✔️:   for (int termIdx = 0; termIdx < numTerms; termIdx++) {
        // RDKit✔️✔️:     const int d_at1Idx = d_at1Idxs[termIdx];
        // RDKit✔️✔️:     const int d_at2Idx = d_at2Idxs[termIdx];
        // RDKit✔️✔️:
        // RDKit✔️✔️:     double dist = dp_forceField->distance(d_at1Idx, d_at2Idx, pos);
        // RDKit✔️✔️:
        // RDKit✔️✔️:     double *at1Coords = &(pos[3 * d_at1Idx]);
        // RDKit✔️✔️:     double *at2Coords = &(pos[3 * d_at2Idx]);
        // RDKit✔️✔️:     double *g1 = &(grad[3 * d_at1Idx]);
        // RDKit✔️✔️:     double *g2 = &(grad[3 * d_at2Idx]);
        // RDKit✔️✔️:
        // RDKit✔️✔️:     double distTerm = dist - d_r0[termIdx];
        // RDKit✔️✔️:     double dE_dr =
        // RDKit✔️✔️:         c1 * d_kb[termIdx] * distTerm *
        // RDKit✔️✔️:         (1.0 + 1.5 * cs * distTerm + 2.0 * c3 * cs * cs * distTerm * distTerm);
        // RDKit✔️✔️:     double dGrad;
        // RDKit✔️✔️:     for (unsigned int i = 0; i < 3; ++i) {
        // RDKit✔️✔️:       dGrad = ((dist > 0.0) ? (dE_dr * (at1Coords[i] - at2Coords[i]) / dist)
        // RDKit✔️✔️:                             : d_kb[termIdx] * 0.01);
        // RDKit✔️✔️:       g1[i] += dGrad;
        // RDKit✔️✔️:       g2[i] -= dGrad;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // RDKit✔️✔️: constexpr double MDYNE_A_TO_KCAL_MOL = 143.9325;
        // Behavior: preserve signed32 term-count conversion/order, association
        // of the derivative, three-component
        // addressing and sequential addition/subtraction (even same endpoints).
        // `dist > 0.0` intentionally sends NaN to the source fallback. Valid
        // coordinate/gradient storage is the existing kernel caller contract.
        // Complexity: O(max(0, signed32(terms))), three axis updates per term, O(1)
        // temporary space, indexed cache reuse, no allocations or state clones.
        let cs = -2.0;
        let c1 = 143.9325;
        let c3 = 7.0 / 12.0;
        for term_idx in source_term_indices(self.at1_idxs.len()) {
            let term_idx = term_idx as usize;
            let at1 = self.at1_idxs[term_idx];
            let at2 = self.at2_idxs[term_idx];
            let dist = context.distance(at1, at2)?;
            let at1_offset = 3 * at1 as usize;
            let at2_offset = 3 * at2 as usize;
            let coordinates = context.coordinates();
            let dist_term = dist - self.r0[term_idx];
            let de_dr = c1
                * self.kb[term_idx]
                * dist_term
                * (1.0 + 1.5 * cs * dist_term + 2.0 * c3 * cs * cs * dist_term * dist_term);
            for i in 0..3 {
                let d_grad = if dist > 0.0 {
                    de_dr * (coordinates[at1_offset + i] - coordinates[at2_offset + i]) / dist
                } else {
                    self.kb[term_idx] * 0.01
                };
                gradient[at1_offset + i] += d_grad;
                gradient[at2_offset + i] -= d_grad;
            }
        }
        Ok(())
    }

    fn copy(&self) -> Box<dyn ForceFieldContribution> {
        // RDKit✔️✔️:   BondStretchContrib *copy() const override {
        // RDKit✔️✔️:     return new BondStretchContrib(*this);
        // RDKit✔️✔️:   }
        // Behavior: independent copies of all arrays; evaluation receives the
        // copied field's context through the existing contribution interface.
        // Complexity: O(terms) element copies and four array allocations plus
        // one contribution allocation, matching the source's vector copy.
        Box::new(self.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernel::{
        ForceFieldIndexArgument, cf3d_bld_b05_calc_energy, cf3d_bld_b05_calc_grad,
    };

    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 1.0e-9,
            "{actual:.17} != {expected:.17}"
        );
    }

    fn contribution(terms: &[(u32, u32, f64, f64)]) -> BondStretchContrib {
        let mut rows = [[0.0; 3]; 3];
        let mut owner = ForceField::new(3);
        owner
            .positions_mut()
            .extend(rows.iter_mut().map(|row| row.as_mut_slice()));
        let mut result = BondStretchContrib::new(&owner);
        for &(i, j, r0, kb) in terms {
            result
                .add_term(owner.positions(), i, j, &MmffBond { r0, kb })
                .unwrap();
        }
        result
    }

    fn energy(term: &dyn ForceFieldContribution, coordinates: &[f64]) -> f64 {
        let n = (coordinates.len() / 3) as u32;
        let mut cache = vec![-1.0; (n * (n + 1) / 2) as usize];
        term.get_energy(&mut EvaluationContext::for_test(coordinates, &mut cache, n))
            .unwrap()
    }

    fn gradient(
        term: &dyn ForceFieldContribution,
        coordinates: &[f64],
        initial: &[f64],
    ) -> Vec<f64> {
        let n = (coordinates.len() / 3) as u32;
        let mut cache = vec![-1.0; (n * (n + 1) / 2) as usize];
        let mut result = initial.to_vec();
        term.get_grad(
            &mut EvaluationContext::for_test(coordinates, &mut cache, n),
            &mut result,
        )
        .unwrap();
        result
    }

    #[test]
    fn mmff_bond_stretch_empty_contribution() {
        let term = contribution(&[]);
        assert_eq!(energy(&term, &[0.0; 6]).to_bits(), 0.0_f64.to_bits());
        assert_eq!(gradient(&term, &[0.0; 6], &[1.0; 6]), [1.0; 6]);
        assert_eq!(energy(&*term.copy(), &[0.0; 6]), 0.0);
    }

    #[test]
    fn mmff_bond_stretch_endpoint_errors_precede_append() {
        let mut row = [0.0; 3];
        let mut owner = ForceField::new(3);
        owner.positions_mut().push(&mut row);
        let mut term = BondStretchContrib::new(&owner);
        let params = MmffBond { r0: 1.0, kb: 2.0 };
        term.add_term(owner.positions(), 0, 0, &params).unwrap();
        for (i, j, argument, index) in [
            (1, 1, BondIndexArgument::First, 1),
            (0, 1, BondIndexArgument::Second, 1),
            (u32::MAX, 0, BondIndexArgument::First, u32::MAX),
        ] {
            assert_eq!(
                term.add_term(owner.positions(), i, j, &params),
                Err(ForceFieldKernelError::BondIndexOutOfRange {
                    argument,
                    index,
                    upper_bound: 1
                })
            );
            assert_eq!(term.at1_idxs, [0]);
            assert_eq!(term.at2_idxs, [0]);
            assert_eq!(term.r0, [1.0]);
            assert_eq!(term.kb, [2.0]);
        }
        // The current owner's position count is checked at each append.
        let mut second_row = [0.0; 3];
        owner.positions_mut().push(&mut second_row);
        term.add_term(owner.positions(), 0, 1, &params).unwrap();
        assert_eq!(term.at2_idxs, [0, 1]);
    }

    #[test]
    fn mmff_bond_stretch_parameters_and_copy_are_independent() {
        let mut rows = [[0.0; 3]; 2];
        let mut owner = ForceField::new(3);
        owner
            .positions_mut()
            .extend(rows.iter_mut().map(|row| row.as_mut_slice()));
        let mut term = BondStretchContrib::new(&owner);
        let mut params = MmffBond { r0: 1.0, kb: 2.0 };
        term.add_term(owner.positions(), 0, 1, &params).unwrap();
        params.r0 = 99.0;
        params.kb = -7.0;
        assert_eq!(params, MmffBond { r0: 99.0, kb: -7.0 });
        assert_eq!(term.r0, [1.0]);
        assert_eq!(term.kb, [2.0]);
        let clone = term.clone();
        assert_ne!(term.at1_idxs.as_ptr(), clone.at1_idxs.as_ptr());
        assert_ne!(term.at2_idxs.as_ptr(), clone.at2_idxs.as_ptr());
        assert_ne!(term.r0.as_ptr(), clone.r0.as_ptr());
        assert_ne!(term.kb.as_ptr(), clone.kb.as_ptr());
        let copied = term.copy();
        term.add_term(owner.positions(), 0, 1, &MmffBond { r0: 1.0, kb: 2.0 })
            .unwrap();
        let coordinates = [0.0, 0.0, 0.0, 2.0, 0.0, 0.0];
        close(energy(&term, &coordinates), 383.82);
        close(energy(&*copied, &coordinates), 191.91);
        assert_eq!(
            gradient(&clone, &coordinates, &[0.0; 6]),
            gradient(&*copied, &coordinates, &[0.0; 6])
        );
    }

    #[test]
    fn mmff_bond_stretch_fixed_equilibrium_compression_extension() {
        // Exact arithmetic witnesses from the source polynomial for r0=1,
        // kb=2: delta=0, -1/2, 1. These literals are not refreshed by tests.
        let term = contribution(&[(0, 1, 1.0, 2.0)]);
        for (distance, expected_energy, endpoint_gradient) in [
            (1.0, 0.0, 0.0),
            (0.5, 92.95640625, 527.7525),
            (2.0, 191.91, -767.64),
        ] {
            let coordinates = [0.0, 0.0, 0.0, distance, 0.0, 0.0];
            close(energy(&term, &coordinates), expected_energy);
            let grad = gradient(&term, &coordinates, &[0.0; 6]);
            close(grad[0], endpoint_gradient);
            close(grad[3], -endpoint_gradient);
            assert_eq!([grad[1], grad[2], grad[4], grad[5]], [0.0; 4]);
        }
    }

    #[test]
    fn mmff_bond_stretch_oblique_translation_and_additive_gradient() {
        let term = contribution(&[(0, 1, 4.0, 2.0)]);
        let coordinates = [0.0, 0.0, 0.0, 3.0, 4.0, 0.0];
        close(energy(&term, &coordinates), 191.91);
        let initial = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let expected = [-459.584, -612.112, 3.0, 464.584, 619.112, 6.0];
        for (a, e) in gradient(&term, &coordinates, &initial).iter().zip(expected) {
            close(*a, e);
        }
        let translated = [8.0, -3.0, 7.0, 11.0, 1.0, 7.0];
        close(energy(&term, &translated), 191.91);
        assert_eq!(
            gradient(&term, &coordinates, &initial),
            gradient(&term, &translated, &initial)
        );
        let reversed = contribution(&[(1, 0, 4.0, 2.0)]);
        assert_eq!(
            gradient(&term, &coordinates, &initial),
            gradient(&reversed, &coordinates, &initial)
        );
    }

    #[test]
    fn mmff_bond_stretch_multiple_repeated_reversed_terms() {
        let term = contribution(&[(0, 1, 1.0, 2.0), (1, 2, 1.0, 7.0), (1, 0, 1.0, 1.0)]);
        let coordinates = [0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 3.0, 0.0, 0.0];
        close(energy(&term, &coordinates), 287.865);
        let grad = gradient(&term, &coordinates, &[0.0; 9]);
        close(grad[0], -1151.46);
        close(grad[3], 1151.46);
        assert_eq!(grad[6..], [0.0; 3]);
        assert_eq!(term.at1_idxs, [0, 1, 1]);
        assert_eq!(term.at2_idxs, [1, 2, 0]);
        assert_eq!(term.kb, [2.0, 7.0, 1.0]);
    }

    #[test]
    fn mmff_bond_stretch_zero_and_nan_distance_fallback() {
        let term = contribution(&[(0, 1, 1.0, 2.0)]);
        close(energy(&term, &[0.0; 6]), 767.64);
        assert_eq!(
            gradient(&term, &[0.0; 6], &[0.0; 6]),
            [0.02, 0.02, 0.02, -0.02, -0.02, -0.02]
        );
        let nan_coordinates = [f64::NAN, 0.0, 0.0, 0.0, 0.0, 0.0];
        assert!(energy(&term, &nan_coordinates).is_nan());
        assert_eq!(
            gradient(&term, &nan_coordinates, &[0.0; 6]),
            [0.02, 0.02, 0.02, -0.02, -0.02, -0.02]
        );
    }

    #[test]
    fn mmff_bond_stretch_same_endpoint_updates_cancel_sequentially() {
        let term = contribution(&[(0, 0, 1.0, 2.0)]);
        close(energy(&term, &[0.0; 3]), 767.64);
        let actual = gradient(&term, &[0.0; 3], &[1.0, 2.0, 3.0]);
        for (a, e) in actual.iter().zip([1.0, 2.0, 3.0]) {
            close(*a, e);
        }
        assert_eq!(gradient(&term, &[0.0; 3], &[0.0; 3]), [0.0; 3]);
    }

    #[test]
    fn mmff_bond_stretch_ieee_parameters_follow_source_arithmetic() {
        let coordinates = [0.0, 0.0, 0.0, 2.0, 0.0, 0.0];
        let negative = contribution(&[(0, 1, 1.0, -2.0)]);
        close(energy(&negative, &coordinates), -191.91);
        close(gradient(&negative, &coordinates, &[0.0; 6])[0], 767.64);
        let zero = contribution(&[(0, 1, 1.0, 0.0)]);
        assert_eq!(energy(&zero, &coordinates), 0.0);
        assert_eq!(gradient(&zero, &coordinates, &[0.0; 6]), [0.0; 6]);
        let nan_rest = contribution(&[(0, 1, f64::NAN, 2.0)]);
        assert!(energy(&nan_rest, &coordinates).is_nan());
        assert!(
            gradient(&nan_rest, &coordinates, &[0.0; 6])
                .iter()
                .all(|v| v.is_nan())
        );
        // Degenerate distance uses kb alone even when rest length is NaN.
        assert_eq!(
            gradient(&nan_rest, &[0.0; 6], &[0.0; 6]),
            [0.02, 0.02, 0.02, -0.02, -0.02, -0.02]
        );
        let infinite = contribution(&[(0, 1, 1.0, f64::INFINITY)]);
        assert_eq!(energy(&infinite, &coordinates), f64::INFINITY);
        let grad = gradient(&infinite, &coordinates, &[0.0; 6]);
        assert_eq!(grad[0], f64::NEG_INFINITY);
        assert_eq!(grad[3], f64::INFINITY);
        assert!(grad[1].is_nan());
        let nan_kb = contribution(&[(0, 1, 1.0, f64::NAN)]);
        assert!(
            gradient(&nan_kb, &[0.0; 6], &[0.0; 6])
                .iter()
                .all(|v| v.is_nan())
        );
    }

    #[test]
    fn mmff_bond_stretch_distance_error_propagates() {
        let term = contribution(&[(0, 1, 1.0, 2.0)]);
        let mut cache = [-1.0];
        let mut context = EvaluationContext::for_test(&[0.0; 3], &mut cache, 1);
        let expected = Err(ForceFieldKernelError::IndexOutOfRange {
            argument: ForceFieldIndexArgument::J,
            index: 1,
            upper_bound: 1,
        });
        assert_eq!(term.get_energy(&mut context), expected);
        let mut grad = [3.0; 3];
        assert_eq!(term.get_grad(&mut context, &mut grad), expected.map(|_| ()));
        assert_eq!(grad, [3.0; 3]);
    }

    #[test]
    fn mmff_bond_stretch_real_field_initialization_and_copy() {
        let mut first = [0.0; 3];
        let mut second = [2.0, 0.0, 0.0];
        let mut field = ForceField::new(3);
        field.positions_mut().push(&mut first);
        field.positions_mut().push(&mut second);
        let mut term = BondStretchContrib::new(&field);
        term.add_term(field.positions(), 0, 1, &MmffBond { r0: 1.0, kb: 2.0 })
            .unwrap();
        field.add_contribution(Box::new(term));
        let coordinates = [0.0, 0.0, 0.0, 2.0, 0.0, 0.0];
        assert_eq!(
            cf3d_bld_b05_calc_energy(&mut field, &coordinates),
            Err(ForceFieldKernelError::NotInitialized)
        );
        field.initialize().unwrap();
        close(field.calc_energy_current(None).unwrap(), 191.91);
        let mut copied = field.copy();
        // ForceField.cpp:159-170 resets df_init and leaves positions empty.
        assert!(copied.positions().is_empty());
        assert_eq!(
            cf3d_bld_b05_calc_energy(&mut copied, &coordinates),
            Err(ForceFieldKernelError::NotInitialized)
        );
        let mut copied_first = [0.0; 3];
        let mut copied_second = [2.0, 0.0, 0.0];
        copied.positions_mut().push(&mut copied_first);
        copied.positions_mut().push(&mut copied_second);
        copied.initialize().unwrap();
        close(
            cf3d_bld_b05_calc_energy(&mut copied, &coordinates).unwrap(),
            191.91,
        );
        let mut grad = [0.0; 6];
        cf3d_bld_b05_calc_grad(&mut copied, &coordinates, &mut grad).unwrap();
        close(grad[0], -767.64);
        let equilibrium = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        assert_eq!(
            cf3d_bld_b05_calc_energy(&mut copied, &equilibrium).unwrap(),
            0.0
        );
        close(field.calc_energy_current(None).unwrap(), 191.91);
    }

    #[test]
    fn mmff_bond_stretch_analytical_gradient_matches_energy_derivative() {
        let term = contribution(&[(0, 1, 1.0, 2.0), (1, 2, 1.2, 0.7)]);
        let coordinates = [0.0, 0.1, -0.2, 1.5, 0.3, 0.4, 1.8, 1.4, 0.8];
        let grad = gradient(&term, &coordinates, &[0.0; 9]);
        for i in 0..9 {
            let mut plus = coordinates;
            let mut minus = coordinates;
            plus[i] += 1.0e-6;
            minus[i] -= 1.0e-6;
            let derivative = (energy(&term, &plus) - energy(&term, &minus)) / 2.0e-6;
            assert!(
                (grad[i] - derivative).abs() < 1.0e-6,
                "axis {i}: {} != {derivative}",
                grad[i]
            );
        }
    }

    #[test]
    fn mmff_bond_stretch_signed_term_count_source_boundaries() {
        // Fixed observations from ROOT's retained C++20 scalar native probe:
        // sizeof(int)=4. Test the exact range consumed by both evaluations
        // without allocating contribution arrays or traversing huge counts.
        let cases = [
            (0_usize, 0_i32, 0_usize, None),
            (1, 1, 1, Some(0)),
            (
                2_147_483_647,
                2_147_483_647,
                2_147_483_647,
                Some(2_147_483_646),
            ),
            (2_147_483_648, -2_147_483_648, 0, None),
            (4_294_967_295, -1, 0, None),
            (4_294_967_296, 0, 0, None),
            (4_294_967_297, 1, 1, Some(0)),
        ];
        for (length, signed_end, iterations, last_index) in cases {
            let indices = source_term_indices(length);
            assert_eq!(indices.start, 0, "length {length}");
            assert_eq!(indices.end, signed_end, "length {length}");
            assert_eq!(indices.len(), iterations, "length {length}");
            assert_eq!(
                indices.clone().next(),
                last_index.map(|_| 0),
                "length {length}"
            );
            assert_eq!(indices.last(), last_index, "length {length}");
        }
    }
}
