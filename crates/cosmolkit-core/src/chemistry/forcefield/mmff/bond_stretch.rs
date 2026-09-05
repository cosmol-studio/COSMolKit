//! Source-backed RDKit MMFF bond stretch contribution.

use crate::chemistry::forcefield::core::{ForceField, ForceFieldContrib};

use super::params::MmffBond;

// BEGIN RDKIT CPP CONSTANT ForceFields::MMFF::MDYNE_A_TO_KCAL_MOL (Params.h:37-40)
// RDKit✔️✔️: constexpr double DEG2RAD = M_PI / 180.0;
// RDKit✔️✔️: constexpr double RAD2DEG = 180.0 / M_PI;
// RDKit✔️✔️: constexpr double MDYNE_A_TO_KCAL_MOL = 143.9325;
const MDYNE_A_TO_KCAL_MOL: f64 = 143.9325;

#[derive(Clone, Debug)]
pub struct BondStretchContrib {
    owner: *const ForceField,
    atom1_indices: Vec<usize>,
    atom2_indices: Vec<usize>,
    r0: Vec<f64>,
    kb: Vec<f64>,
}

impl BondStretchContrib {
    #[must_use]
    pub fn new(owner: &ForceField) -> Self {
        // BEGIN RDKIT CPP CONSTRUCTOR ForceFields::MMFF::BondStretchContrib::BondStretchContrib (BondStretch.cpp:45-52)
        // RDKit✔️✔️: BondStretchContrib::BondStretchContrib(ForceField *owner) {
        // RDKit✔️✔️:   PRECONDITION(owner, "bad owner");
        // Rust references reproduce RDKit's non-null owner precondition.
        // RDKit✔️✔️:
        //
        // RDKit✔️✔️:
        // RDKit✔️✔️:   dp_forceField = owner;
        let owner = owner as *const ForceField;
        // RDKit✔️✔️:
        // RDKit✔️✔️: }
        Self {
            owner,
            atom1_indices: Vec::new(),
            atom2_indices: Vec::new(),
            r0: Vec::new(),
            kb: Vec::new(),
        }
    }

    #[must_use]
    pub fn owner(&self) -> *const ForceField {
        self.owner
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.atom1_indices.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.atom1_indices.is_empty() && self.atom2_indices.is_empty() && self.r0.is_empty() && self.kb.is_empty()
    }

    pub fn add_term(&mut self, idx1: usize, idx2: usize, mmff_bond_params: &MmffBond) {
        // BEGIN RDKIT CPP METHOD ForceFields::MMFF::BondStretchContrib::addTerm (BondStretch.cpp:53-64)
        // RDKit✔️✔️: void BondStretchContrib::addTerm(const unsigned int idx1,
        // RDKit✔️✔️:                                  const unsigned int idx2,
        // RDKit✔️✔️:                                  const ForceFields::MMFF::MMFFBond *mmffBondParams) {
        let force_field = self.force_field();
        // RDKit✔️✔️:   URANGE_CHECK(idx1, dp_forceField->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx2, dp_forceField->positions().size());
        assert!(idx1 < force_field.positions().len());
        assert!(idx2 < force_field.positions().len());
        // RDKit✔️✔️:   PRECONDITION(mmffBondParams, "bond parameters not found");
        // Rust references reproduce RDKit's non-null parameter precondition.

        // RDKit✔️✔️:   d_at1Idxs.push_back(idx1);
        // RDKit✔️✔️:   d_at2Idxs.push_back(idx2);
        // RDKit✔️✔️:   d_r0.push_back(mmffBondParams->r0);
        // RDKit✔️✔️:   d_kb.push_back(mmffBondParams->kb);
        self.atom1_indices.push(idx1);
        self.atom2_indices.push(idx2);
        self.r0.push(mmff_bond_params.r0);
        self.kb.push(mmff_bond_params.kb);
        // RDKit✔️✔️: }
    }

    #[must_use]
    pub fn get_energy(&self, pos: &[f64]) -> f64 {
        // BEGIN RDKIT CPP METHOD ForceFields::MMFF::BondStretchContrib::getEnergy (BondStretch.cpp:66-77)
        // RDKit✔️✔️: double BondStretchContrib::getEnergy(double *pos) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        let force_field = self.force_field();
        // Rust slices reproduce RDKit's non-null pos precondition.

        // RDKit✔️✔️:   const int numTerms = d_at1Idxs.size();
        // RDKit✔️✔️:   double energySum = 0.0;
        let num_terms = self.atom1_indices.len();
        let mut energy_sum = 0.0;
        // RDKit✔️✔️:   for (int i =0; i < numTerms; i++) {
        for i in 0..num_terms {
            // RDKit✔️✔️:     energySum += Utils::calcBondStretchEnergy(
            // RDKit✔️✔️:         d_r0[i], d_kb[i], dp_forceField->distance(d_at1Idxs[i], d_at2Idxs[i], pos));
            energy_sum += calc_bond_stretch_energy(
                self.r0[i],
                self.kb[i],
                force_field.distance(self.atom1_indices[i], self.atom2_indices[i], Some(pos)),
            );
            // RDKit✔️✔️:   }
        }
        // RDKit✔️✔️:   return energySum;
        // RDKit✔️✔️: }
        energy_sum
    }

    pub fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
        // BEGIN RDKIT CPP METHOD ForceFields::MMFF::BondStretchContrib::getGrad (BondStretch.cpp:79-117)
        // RDKit✔️✔️: void BondStretchContrib::getGrad(double *pos, double *grad) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:   PRECONDITION(grad, "bad vector");
        let force_field = self.force_field();
        // Rust slices reproduce RDKit's non-null pos and grad preconditions.

        // RDKit✔️✔️:   const int numTerms = d_at1Idxs.size();
        // RDKit✔️✔️:   constexpr double cs = -2.0;
        // RDKit✔️✔️:   constexpr double c1 = MDYNE_A_TO_KCAL_MOL;
        // RDKit✔️✔️:   constexpr double c3 = 7.0 / 12.0;
        let num_terms = self.atom1_indices.len();
        let cs = -2.0;
        let c1 = MDYNE_A_TO_KCAL_MOL;
        let c3 = 7.0 / 12.0;
        // RDKit✔️✔️:   for (int termIdx = 0; termIdx < numTerms; termIdx++) {
        for term_idx in 0..num_terms {
            // RDKit✔️✔️:     const int d_at1Idx = d_at1Idxs[termIdx];
            // RDKit✔️✔️:     const int d_at2Idx = d_at2Idxs[termIdx];
            let atom1_idx = self.atom1_indices[term_idx];
            let atom2_idx = self.atom2_indices[term_idx];

            // RDKit✔️✔️:     double dist = dp_forceField->distance(d_at1Idx, d_at2Idx, pos);
            let dist = force_field.distance(atom1_idx, atom2_idx, Some(pos));

            // RDKit✔️✔️:     double *at1Coords = &(pos[3 * d_at1Idx]);
            // RDKit✔️✔️:     double *at2Coords = &(pos[3 * d_at2Idx]);
            // RDKit✔️✔️:     double *g1 = &(grad[3 * d_at1Idx]);
            // RDKit✔️✔️:     double *g2 = &(grad[3 * d_at2Idx]);
            let atom1_offset = 3 * atom1_idx;
            let atom2_offset = 3 * atom2_idx;

            // RDKit✔️✔️:     double distTerm = dist - d_r0[termIdx];
            let dist_term = dist - self.r0[term_idx];
            // RDKit✔️✔️:     double dE_dr =
            // RDKit✔️✔️:         c1 * d_kb[termIdx] * distTerm *
            // RDKit✔️✔️:         (1.0 + 1.5 * cs * distTerm + 2.0 * c3 * cs * cs * distTerm * distTerm);
            let de_dr = c1
                * self.kb[term_idx]
                * dist_term
                * (1.0 + 1.5 * cs * dist_term + 2.0 * c3 * cs * cs * dist_term * dist_term);
            // RDKit✔️✔️:     double dGrad;
            // RDKit✔️✔️:     for (unsigned int i = 0; i < 3; ++i) {
            for i in 0..3 {
                // RDKit✔️✔️:       dGrad = ((dist > 0.0) ? (dE_dr * (at1Coords[i] - at2Coords[i]) / dist)
                // RDKit✔️✔️:                             : d_kb[termIdx] * 0.01);
                let d_grad = if dist > 0.0 {
                    de_dr * (pos[atom1_offset + i] - pos[atom2_offset + i]) / dist
                } else {
                    self.kb[term_idx] * 0.01
                };
                // RDKit✔️✔️:       g1[i] += dGrad;
                // RDKit✔️✔️:       g2[i] -= dGrad;
                grad[atom1_offset + i] += d_grad;
                grad[atom2_offset + i] -= d_grad;
                // RDKit✔️✔️:     }
            }
            // RDKit✔️✔️:   }
        }
        // RDKit✔️✔️: }
    }

    #[must_use]
    pub fn atom1_indices(&self) -> &[usize] {
        &self.atom1_indices
    }

    #[must_use]
    pub fn atom2_indices(&self) -> &[usize] {
        &self.atom2_indices
    }

    #[must_use]
    pub fn rest_lengths(&self) -> &[f64] {
        &self.r0
    }

    #[must_use]
    pub fn force_constants(&self) -> &[f64] {
        &self.kb
    }

    fn force_field(&self) -> &ForceField {
        assert!(!self.owner.is_null(), "no owner");
        // SAFETY: MMFF contribs follow the same owner-pointer model as
        // ForceFieldContrib objects in core.rs. Constructors store a live
        // ForceField pointer before term insertion or evaluation.
        unsafe { &*self.owner }
    }
}

#[must_use]
fn calc_bond_stretch_energy(r0: f64, kb: f64, distance: f64) -> f64 {
    // BEGIN RDKIT CPP FUNCTION ForceFields::MMFF::Utils::calcBondStretchEnergy (BondStretch.cpp:29-43)
    // RDKit✔️✔️: double calcBondStretchEnergy(const double r0, const double kb,
    // RDKit✔️✔️:                              const double distance) {
    // RDKit✔️✔️:   double distTerm = distance - r0;
    let dist_term = distance - r0;
    // RDKit✔️✔️:   double distTerm2 = distTerm * distTerm;
    let dist_term2 = dist_term * dist_term;
    // RDKit✔️✔️:   double const c1 = MDYNE_A_TO_KCAL_MOL;
    // RDKit✔️✔️:   double const cs = -2.0;
    // RDKit✔️✔️:   double const c3 = 7.0 / 12.0;
    let c1 = MDYNE_A_TO_KCAL_MOL;
    let cs = -2.0;
    let c3 = 7.0 / 12.0;

    // RDKit✔️✔️:   return (0.5 * c1 * kb * distTerm2 *
    // RDKit✔️✔️:           (1.0 + cs * distTerm + c3 * cs * cs * distTerm2));
    // RDKit✔️✔️: }
    0.5 * c1 * kb * dist_term2 * (1.0 + cs * dist_term + c3 * cs * cs * dist_term2)
}

impl ForceFieldContrib for BondStretchContrib {
    fn copy(&self) -> Box<dyn ForceFieldContrib> {
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
    use super::BondStretchContrib;
    use crate::chemistry::forcefield::{
        core::{ForceField, ForceFieldVec3},
        mmff::params::MmffBond,
    };

    fn force_field_with_positions(count: usize) -> ForceField {
        let mut force_field = ForceField::new(3);
        force_field
            .positions_mut()
            .extend((0..count).map(|idx| ForceFieldVec3::new(idx as f64, 0.0, 0.0)));
        force_field
    }

    fn initialized_force_field_with_positions(count: usize) -> ForceField {
        let mut force_field = force_field_with_positions(count);
        force_field.initialize();
        force_field
    }

    fn source_bond_stretch_energy(r0: f64, kb: f64, distance: f64) -> f64 {
        let dist_term = distance - r0;
        let dist_term2 = dist_term * dist_term;
        let c1 = 143.9325;
        let cs = -2.0;
        let c3 = 7.0 / 12.0;
        0.5 * c1 * kb * dist_term2 * (1.0 + cs * dist_term + c3 * cs * cs * dist_term2)
    }

    fn source_bond_stretch_de_dr(r0: f64, kb: f64, distance: f64) -> f64 {
        let dist_term = distance - r0;
        let c1 = 143.9325;
        let cs = -2.0;
        let c3 = 7.0 / 12.0;
        c1 * kb * dist_term * (1.0 + 1.5 * cs * dist_term + 2.0 * c3 * cs * cs * dist_term * dist_term)
    }

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 1.0e-12,
            "actual={actual}, expected={expected}"
        );
    }

    fn assert_slice_close(actual: &[f64], expected: &[f64]) {
        assert_eq!(actual.len(), expected.len());
        for (idx, (actual, expected)) in actual.iter().zip(expected.iter()).enumerate() {
            assert!(
                (actual - expected).abs() < 1.0e-12,
                "idx={idx}, actual={actual}, expected={expected}"
            );
        }
    }

    #[test]
    fn mmff_bondstretchcontrib_constructor_stores_owner_pointer() {
        let force_field = ForceField::new(3);

        let contrib = BondStretchContrib::new(&force_field);

        assert_eq!(contrib.owner(), &force_field as *const ForceField);
    }

    #[test]
    fn mmff_bondstretchcontrib_constructor_initializes_no_terms() {
        let force_field = ForceField::new(3);

        let contrib = BondStretchContrib::new(&force_field);

        assert_eq!(contrib.len(), 0);
        assert!(contrib.is_empty());
    }

    #[test]
    fn mmff_bondstretchcontrib_constructor_accepts_empty_force_field_like_rdkit() {
        let force_field = ForceField::new(3);

        let contrib = BondStretchContrib::new(&force_field);

        assert_eq!(force_field.positions().len(), 0);
        assert!(contrib.is_empty());
    }

    #[test]
    fn mmff_bondstretchcontrib_add_term_pushes_source_fields() {
        let force_field = force_field_with_positions(2);
        let mut contrib = BondStretchContrib::new(&force_field);
        let params = MmffBond { kb: 4.258, r0: 1.508 };

        contrib.add_term(0, 1, &params);

        assert_eq!(contrib.atom1_indices(), &[0]);
        assert_eq!(contrib.atom2_indices(), &[1]);
        assert_eq!(contrib.rest_lengths(), &[1.508]);
        assert_eq!(contrib.force_constants(), &[4.258]);
    }

    #[test]
    fn mmff_bondstretchcontrib_add_term_appends_multiple_terms() {
        let force_field = force_field_with_positions(3);
        let mut contrib = BondStretchContrib::new(&force_field);
        let first = MmffBond { kb: 4.258, r0: 1.508 };
        let second = MmffBond { kb: 5.084, r0: 1.451 };

        contrib.add_term(0, 1, &first);
        contrib.add_term(1, 2, &second);

        assert_eq!(contrib.atom1_indices(), &[0, 1]);
        assert_eq!(contrib.atom2_indices(), &[1, 2]);
        assert_eq!(contrib.rest_lengths(), &[1.508, 1.451]);
        assert_eq!(contrib.force_constants(), &[4.258, 5.084]);
    }

    #[test]
    #[should_panic]
    fn mmff_bondstretchcontrib_add_term_rejects_first_index_out_of_range() {
        let force_field = force_field_with_positions(2);
        let mut contrib = BondStretchContrib::new(&force_field);
        let params = MmffBond { kb: 4.258, r0: 1.508 };

        contrib.add_term(2, 1, &params);
    }

    #[test]
    #[should_panic]
    fn mmff_bondstretchcontrib_add_term_rejects_second_index_out_of_range() {
        let force_field = force_field_with_positions(2);
        let mut contrib = BondStretchContrib::new(&force_field);
        let params = MmffBond { kb: 4.258, r0: 1.508 };

        contrib.add_term(0, 2, &params);
    }

    #[test]
    fn mmff_bondstretchcontrib_get_energy_returns_zero_without_terms() {
        let force_field = ForceField::new(3);
        let contrib = BondStretchContrib::new(&force_field);

        assert_eq!(contrib.get_energy(&[]), 0.0);
    }

    #[test]
    fn mmff_bondstretchcontrib_get_energy_uses_source_formula_for_single_term() {
        let force_field = initialized_force_field_with_positions(2);
        let mut contrib = BondStretchContrib::new(&force_field);
        let params = MmffBond { kb: 4.258, r0: 1.508 };
        contrib.add_term(0, 1, &params);

        let pos = [0.0, 0.0, 0.0, 1.6, 0.0, 0.0];
        let expected = source_bond_stretch_energy(1.508, 4.258, 1.6);

        assert_close(contrib.get_energy(&pos), expected);
    }

    #[test]
    fn mmff_bondstretchcontrib_get_energy_sums_multiple_terms() {
        let force_field = initialized_force_field_with_positions(3);
        let mut contrib = BondStretchContrib::new(&force_field);
        let first = MmffBond { kb: 4.258, r0: 1.508 };
        let second = MmffBond { kb: 5.084, r0: 1.451 };
        contrib.add_term(0, 1, &first);
        contrib.add_term(1, 2, &second);

        let pos = [0.0, 0.0, 0.0, 1.6, 0.0, 0.0, 2.9, 0.0, 0.0];
        let expected = source_bond_stretch_energy(1.508, 4.258, 1.6) + source_bond_stretch_energy(1.451, 5.084, 1.3);

        assert_close(contrib.get_energy(&pos), expected);
    }

    #[test]
    fn mmff_bondstretchcontrib_get_energy_uses_supplied_position_vector() {
        let force_field = initialized_force_field_with_positions(2);
        let mut contrib = BondStretchContrib::new(&force_field);
        let params = MmffBond { kb: 4.258, r0: 1.508 };
        contrib.add_term(0, 1, &params);

        let pos = [0.0, 0.0, 0.0, 2.0, 0.0, 0.0];
        let expected = source_bond_stretch_energy(1.508, 4.258, 2.0);

        assert_close(contrib.get_energy(&pos), expected);
    }

    #[test]
    #[should_panic]
    fn mmff_bondstretchcontrib_get_energy_requires_initialized_force_field_for_terms() {
        let force_field = force_field_with_positions(2);
        let mut contrib = BondStretchContrib::new(&force_field);
        let params = MmffBond { kb: 4.258, r0: 1.508 };
        contrib.add_term(0, 1, &params);

        let _ = contrib.get_energy(&[0.0, 0.0, 0.0, 1.6, 0.0, 0.0]);
    }

    #[test]
    fn mmff_bondstretchcontrib_get_grad_leaves_gradient_without_terms() {
        let force_field = ForceField::new(3);
        let contrib = BondStretchContrib::new(&force_field);
        let mut grad = [1.0, 2.0, 3.0];

        contrib.get_grad(&[], &mut grad);

        assert_eq!(grad, [1.0, 2.0, 3.0]);
    }

    #[test]
    fn mmff_bondstretchcontrib_get_grad_accumulates_nonzero_distance_gradient() {
        let force_field = initialized_force_field_with_positions(2);
        let mut contrib = BondStretchContrib::new(&force_field);
        let params = MmffBond { kb: 4.258, r0: 1.508 };
        contrib.add_term(0, 1, &params);
        let pos = [0.0, 0.0, 0.0, 1.6, 0.0, 0.0];
        let mut grad = [0.0; 6];

        contrib.get_grad(&pos, &mut grad);

        let de_dr = source_bond_stretch_de_dr(1.508, 4.258, 1.6);
        assert_slice_close(&grad, &[-de_dr, 0.0, 0.0, de_dr, 0.0, 0.0]);
    }

    #[test]
    fn mmff_bondstretchcontrib_get_grad_uses_zero_distance_fallback() {
        let force_field = initialized_force_field_with_positions(2);
        let mut contrib = BondStretchContrib::new(&force_field);
        let params = MmffBond { kb: 4.258, r0: 1.508 };
        contrib.add_term(0, 1, &params);
        let pos = [0.0; 6];
        let mut grad = [0.0; 6];

        contrib.get_grad(&pos, &mut grad);

        let fallback = 4.258 * 0.01;
        assert_slice_close(&grad, &[fallback, fallback, fallback, -fallback, -fallback, -fallback]);
    }

    #[test]
    fn mmff_bondstretchcontrib_get_grad_accumulates_multiple_terms() {
        let force_field = initialized_force_field_with_positions(3);
        let mut contrib = BondStretchContrib::new(&force_field);
        let first = MmffBond { kb: 4.258, r0: 1.508 };
        let second = MmffBond { kb: 5.084, r0: 1.451 };
        contrib.add_term(0, 1, &first);
        contrib.add_term(1, 2, &second);
        let pos = [0.0, 0.0, 0.0, 1.6, 0.0, 0.0, 2.9, 0.0, 0.0];
        let mut grad = [0.0; 9];

        contrib.get_grad(&pos, &mut grad);

        let first_de_dr = source_bond_stretch_de_dr(1.508, 4.258, 1.6);
        let second_de_dr = source_bond_stretch_de_dr(1.451, 5.084, 1.3);
        assert_slice_close(
            &grad,
            &[
                -first_de_dr,
                0.0,
                0.0,
                first_de_dr - second_de_dr,
                0.0,
                0.0,
                second_de_dr,
                0.0,
                0.0,
            ],
        );
    }

    #[test]
    #[should_panic]
    fn mmff_bondstretchcontrib_get_grad_requires_initialized_force_field_for_terms() {
        let force_field = force_field_with_positions(2);
        let mut contrib = BondStretchContrib::new(&force_field);
        let params = MmffBond { kb: 4.258, r0: 1.508 };
        contrib.add_term(0, 1, &params);
        let mut grad = [0.0; 6];

        contrib.get_grad(&[0.0, 0.0, 0.0, 1.6, 0.0, 0.0], &mut grad);
    }
}
