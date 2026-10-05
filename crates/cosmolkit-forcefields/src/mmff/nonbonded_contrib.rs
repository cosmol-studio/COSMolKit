//! RDKit MMFF Nonbonded.cpp/Nonbonded.h, BSD license; pinned
//! 351f8f378f8ad6bbd517980c38896e66bf907af8.
//! Migrated packed terms use the sole existing kernel evaluation context.
use super::{
    mol_properties::MmffVdwRijstarEps,
    nonbonded::{calc_ele_energy, calc_vdw_energy},
};
use crate::kernel::{
    BondIndexArgument, EvaluationContext, ForceField, ForceFieldContribution, ForceFieldKernelError,
};
#[derive(Clone, Debug)]
pub(super) struct NonbondedContrib {
    atom1_indices: Vec<i16>,
    atom2_indices: Vec<i16>,
    contrib_types: Vec<u8>,
    r_ij_stars: Vec<f64>,
    well_depths: Vec<f64>,
    charge_terms: Vec<f64>,
    is_1_4s: Vec<u8>,
    diel_models: Vec<u8>,
}

impl NonbondedContrib {
    #[must_use]
    pub(super) fn new(_owner: &ForceField<'_>) -> Self {
        // BEGIN RDKIT CPP CONSTRUCTOR ForceFields::MMFF::NonbondedContrib::NonbondedContrib (Nonbonded.cpp:87-90)
        // RDKit❗✔️: NonbondedContrib::NonbondedContrib(ForceField *owner) {
        // RDKit❗✔️:   PRECONDITION(owner, "bad owner");
        // Rust references reproduce RDKit's non-null owner precondition.
        // RDKit❗✔️:   dp_forceField = owner;

        // RDKit❗✔️: }
        Self {
            atom1_indices: Vec::new(),
            atom2_indices: Vec::new(),
            contrib_types: Vec::new(),
            r_ij_stars: Vec::new(),
            well_depths: Vec::new(),
            charge_terms: Vec::new(),
            is_1_4s: Vec::new(),
            diel_models: Vec::new(),
        }
    }

    #[must_use]
    pub(super) fn len(&self) -> usize {
        self.atom1_indices.len()
    }

    #[must_use]
    pub(super) fn is_empty(&self) -> bool {
        self.atom1_indices.is_empty()
            && self.atom2_indices.is_empty()
            && self.contrib_types.is_empty()
            && self.r_ij_stars.is_empty()
            && self.well_depths.is_empty()
            && self.charge_terms.is_empty()
            && self.is_1_4s.is_empty()
            && self.diel_models.is_empty()
    }

    pub(super) fn add_term(
        &mut self,
        positions: &[&mut [f64]],
        idx1: u32,
        idx2: u32,
        mmff_vdw_constants: Option<&MmffVdwRijstarEps>,
        include_charge: bool,
        charge_term: f64,
        diel_model: u8,
        is_1_4: bool,
    ) -> Result<(), ForceFieldKernelError> {
        // BEGIN RDKIT CPP METHOD ForceFields::MMFF::NonbondedContrib::addTerm (Nonbonded.cpp:268-290)
        // RDKit❗✔️: void NonbondedContrib::addTerm(unsigned int idx1, unsigned int idx2,
        // RDKit❗✔️:                              const MMFFVdWRijstarEps *mmffVdWConstants, bool includeCharge,
        // RDKit❗✔️:                              double chargeTerm, std::uint8_t dielModel, bool is1_4) {
        // RDKit❗✔️:   if (!mmffVdWConstants && !includeCharge) {
        if mmff_vdw_constants.is_none() && !include_charge {
            // RDKit❗✔️:     return;
            return Ok(());
        }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   URANGE_CHECK(idx1, dp_forceField->positions().size());
        // RDKit❗✔️:   URANGE_CHECK(idx2, dp_forceField->positions().size());
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
        // RDKit❗✔️:   d_at1Idxs.push_back(idx1);
        // RDKit❗✔️:   d_at2Idxs.push_back(idx2);
        self.atom1_indices.push(idx1 as i16);
        self.atom2_indices.push(idx2 as i16);
        // RDKit❗✔️:   d_contribTypes.push_back(0);
        self.contrib_types.push(0);
        // RDKit❗✔️:   if (mmffVdWConstants) {
        if let Some(mmff_vdw_constants) = mmff_vdw_constants {
            // RDKit❗✔️:     d_contribTypes.back() |= ContribType::VDW;
            *self.contrib_types.last_mut().expect("term was just pushed") |= 1;
            // RDKit❗✔️:     d_R_ij_stars.push_back(mmffVdWConstants->R_ij_star);
            self.r_ij_stars.push(mmff_vdw_constants.r_ij_star);
            // RDKit❗✔️:     d_wellDepths.push_back(mmffVdWConstants->epsilon);
            self.well_depths.push(mmff_vdw_constants.epsilon);
            // RDKit❗✔️:   } else {
        } else {
            // RDKit❗✔️:     d_R_ij_stars.push_back(0.0);
            self.r_ij_stars.push(0.0);
            // RDKit❗✔️:     d_wellDepths.push_back(0.0);
            self.well_depths.push(0.0);
            // RDKit❗✔️:   }
        }
        // RDKit❗✔️:   if (includeCharge) {
        if include_charge {
            // RDKit❗✔️:     d_contribTypes.back() |= ContribType::ELECTROSTATIC;
            *self.contrib_types.last_mut().expect("term was just pushed") |= 2;
            // RDKit❗✔️:     d_chargeTerms.push_back(chargeTerm);
            self.charge_terms.push(charge_term);
            // RDKit❗✔️:     d_dielModels.push_back(dielModel);
            self.diel_models.push(diel_model);
            // RDKit❗✔️:     d_is_1_4s.push_back(is1_4);
            self.is_1_4s.push(u8::from(is_1_4));
            // RDKit❗✔️:   } else {
        } else {
            // RDKit❗✔️:     d_chargeTerms.push_back(0.0);
            self.charge_terms.push(0.0);
            // RDKit❗✔️:     d_dielModels.push_back(0);
            self.diel_models.push(0);
            // RDKit❗✔️:     d_is_1_4s.push_back(false);
            self.is_1_4s.push(0);
            // RDKit❗✔️:   }
        }
        // RDKit❗✔️: }
        Ok(())
    }

    #[must_use]
    pub(super) fn atom1_indices(&self) -> &[i16] {
        &self.atom1_indices
    }

    #[must_use]
    pub(super) fn atom2_indices(&self) -> &[i16] {
        &self.atom2_indices
    }

    #[must_use]
    pub(super) fn contrib_types(&self) -> &[u8] {
        &self.contrib_types
    }

    #[must_use]
    pub(super) fn r_ij_stars(&self) -> &[f64] {
        &self.r_ij_stars
    }

    #[must_use]
    pub(super) fn well_depths(&self) -> &[f64] {
        &self.well_depths
    }

    #[must_use]
    pub(super) fn charge_terms(&self) -> &[f64] {
        &self.charge_terms
    }

    #[must_use]
    pub(super) fn is_1_4s(&self) -> &[u8] {
        &self.is_1_4s
    }

    #[must_use]
    pub(super) fn diel_models(&self) -> &[u8] {
        &self.diel_models
    }

    #[must_use]
    pub(super) fn energy(
        &self,
        context: &mut EvaluationContext<'_>,
    ) -> Result<f64, ForceFieldKernelError> {
        // BEGIN RDKIT CPP METHOD ForceFields::MMFF::NonbondedContrib::getEnergy (Nonbonded.cpp:295-319)
        // RDKit❗✔️: double NonbondedContrib::getEnergy(double *pos) const {
        // RDKit❗✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit❗✔️:   PRECONDITION(pos, "bad vector");
        // Rust slices reproduce RDKit's non-null pos precondition.
        // RDKit❗✔️:   double energySum = 0.0;
        let mut energy_sum = 0.0;

        // RDKit❗✔️:   const int numPairs = d_at1Idxs.size();
        let num_pairs = self.atom1_indices.len() as i32;
        // RDKit❗✔️:   for (int i = 0; i < numPairs; ++i) {
        for i in 0..num_pairs {
            let i = i as usize;
            // RDKit❗✔️:     unsigned int d_at1Idx = d_at1Idxs[i];
            // RDKit❗✔️:     unsigned int d_at2Idx = d_at2Idxs[i];
            let atom1_idx = self.atom1_indices[i] as u32;
            let atom2_idx = self.atom2_indices[i] as u32;
            // RDKit❗✔️:     double dist = dp_forceField->distance(d_at1Idx, d_at2Idx, pos);
            let dist = context.distance(atom1_idx as u32, atom2_idx as u32)?;

            // RDKit❗✔️:     if (d_contribTypes[i] & ContribType::VDW) {
            if self.contrib_types[i] & 1 != 0 {
                // RDKit❗✔️:       const auto res =
                // RDKit❗✔️:           Utils::calcVdWEnergy(dist, d_R_ij_stars[i], d_wellDepths[i]);
                let res = calc_vdw_energy(dist, self.r_ij_stars[i], self.well_depths[i]);
                // RDKit❗✔️:       energySum += res;
                energy_sum += res;
                // RDKit❗✔️:     }
            }
            // RDKit❗✔️:     if (d_contribTypes[i] & ContribType::ELECTROSTATIC) {
            if self.contrib_types[i] & 2 != 0 {
                // RDKit❗✔️:       const double d_chargeTerm = d_chargeTerms[i];
                let charge_term = self.charge_terms[i];
                // RDKit❗✔️:       const std::uint8_t d_dielModel = d_dielModels[i];
                let diel_model = self.diel_models[i];
                // RDKit❗✔️:       const bool d_is1_4 = d_is_1_4s[i];
                let is_1_4 = self.is_1_4s[i] != 0;
                // RDKit❗✔️:       const auto res = Utils::calcEleEnergy(d_at1Idx, d_at2Idx, dist,
                // RDKit❗✔️:                                             d_chargeTerm, d_dielModel, d_is1_4);
                let res =
                    calc_ele_energy(atom1_idx, atom2_idx, dist, charge_term, diel_model, is_1_4);
                // RDKit❗✔️:       energySum += res;
                energy_sum += res;
                // RDKit❗✔️:     }
            }
            // RDKit❗✔️:   }
        }
        // RDKit❗✔️:   return energySum;
        // RDKit❗✔️: }
        Ok(energy_sum)
    }

    pub(super) fn gradient(
        &self,
        context: &mut EvaluationContext<'_>,
        grad: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        // BEGIN RDKIT CPP METHOD ForceFields::MMFF::NonbondedContrib::getGrad (Nonbonded.cpp:321-380)
        // RDKit❗✔️: void NonbondedContrib::getGrad(double *pos, double *grad) const {
        // RDKit❗✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit❗✔️:   PRECONDITION(pos, "bad vector");
        // RDKit❗✔️:   PRECONDITION(grad, "bad vector");
        // Rust slices reproduce RDKit's non-null pos and grad preconditions.

        // RDKit❗✔️:   constexpr double vdw1 = 1.07;
        let vdw1 = 1.07;
        // RDKit❗✔️:   constexpr double vdw1m1 = vdw1 - 1.0;
        let vdw1m1 = vdw1 - 1.0;
        // RDKit❗✔️:   constexpr double vdw2 = 1.12;
        let vdw2 = 1.12;
        // RDKit❗✔️:   constexpr double vdw2m1 = vdw2 - 1.0;
        let vdw2m1 = vdw2 - 1.0;
        // RDKit❗✔️:   constexpr double vdw2t7 = vdw2 * 7.0;
        let vdw2t7 = vdw2 * 7.0;

        // RDKit❗✔️:   const int numPairs = d_at1Idxs.size();
        let num_pairs = self.atom1_indices.len() as i32;
        // RDKit❗✔️:   for (int pairIdx = 0; pairIdx < numPairs; ++pairIdx) {
        for pair_idx in 0..num_pairs {
            let pair_idx = pair_idx as usize;
            // RDKit❗✔️:     const int d_at1Idx = d_at1Idxs[pairIdx];
            // RDKit❗✔️:     const int d_at2Idx = d_at2Idxs[pairIdx];
            let atom1_idx = i32::from(self.atom1_indices[pair_idx]);
            let atom2_idx = i32::from(self.atom2_indices[pair_idx]);
            // RDKit❗✔️:     const double dist = dp_forceField->distance(d_at1Idx, d_at2Idx, pos);
            let dist = context.distance(atom1_idx as u32, atom2_idx as u32)?;
            // RDKit❗✔️:     const double *at1Coords = &(pos[3 * d_at1Idx]);
            // RDKit❗✔️:     const double *at2Coords = &(pos[3 * d_at2Idx]);
            let atom1_offset = (3 * atom1_idx) as usize;
            let atom2_offset = (3 * atom2_idx) as usize;
            let pos = context.coordinates();
            if atom1_idx < 0
                || atom2_idx < 0
                || atom1_offset
                    .checked_add(2)
                    .is_none_or(|end| end >= pos.len() || end >= grad.len())
                || atom2_offset
                    .checked_add(2)
                    .is_none_or(|end| end >= pos.len() || end >= grad.len())
            {
                return Err(ForceFieldKernelError::BadIndex);
            }

            // RDKit❗✔️:     double vdwGrad = 0.0;
            // RDKit❗✔️:     double eleGrad = 0.0;
            let mut vdw_grad = 0.0;
            let mut ele_grad = 0.0;
            // RDKit❗✔️:     if (dist <= 0.0) {
            if dist <= 0.0 {
                // RDKit❗✔️:       if (d_contribTypes[pairIdx] & ContribType::VDW) {
                if self.contrib_types[pair_idx] & 1 != 0 {
                    // RDKit❗✔️:         const double d_R_ij_star = d_R_ij_stars[pairIdx];
                    let r_ij_star = self.r_ij_stars[pair_idx];
                    // RDKit❗✔️:         for (unsigned int i = 0; i < 3; ++i) {
                    for i in 0..3 {
                        // RDKit❗✔️:           g1[i] += d_R_ij_star * 0.01;
                        // RDKit❗✔️:           g2[i] -= d_R_ij_star * 0.01;
                        grad[atom1_offset + i] += r_ij_star * 0.01;
                        grad[atom2_offset + i] -= r_ij_star * 0.01;
                        // RDKit❗✔️:         }
                    }
                    // RDKit❗✔️:       }
                }
                // RDKit❗✔️:       if (d_contribTypes[pairIdx] & ContribType::ELECTROSTATIC) {
                if self.contrib_types[pair_idx] & 2 != 0 {
                    // RDKit❗✔️:         for (unsigned int i = 0; i < 3; ++i) {
                    for i in 0..3 {
                        // RDKit❗✔️:           g1[i] += 0.02;
                        // RDKit❗✔️:           g2[i] -= 0.02;
                        grad[atom1_offset + i] += 0.02;
                        grad[atom2_offset + i] -= 0.02;
                        // RDKit❗✔️:         }
                    }
                    // RDKit❗✔️:       }
                }
                // RDKit❗✔️:       return;
                return Ok(());
                // RDKit❗✔️:     }
            }
            // RDKit❗✔️:     if (d_contribTypes[pairIdx] & ContribType::VDW) {
            if self.contrib_types[pair_idx] & 1 != 0 {
                // RDKit❗✔️:       const double d_R_ij_star = d_R_ij_stars[pairIdx];
                let r_ij_star = self.r_ij_stars[pair_idx];
                // RDKit❗✔️:       const double d_wellDepth = d_wellDepths[pairIdx];
                let well_depth = self.well_depths[pair_idx];

                // RDKit❗✔️:       const double q = dist / d_R_ij_star;
                let q = dist / r_ij_star;
                // RDKit❗✔️:       const double q2 = q * q;
                let q2 = q * q;
                // RDKit❗✔️:       const double q6 = q2 * q2 * q2;
                let q6 = q2 * q2 * q2;
                // RDKit❗✔️:       const double q7 = q6 * q;
                let q7 = q6 * q;
                // RDKit❗✔️:       const double q7pvdw2m1 = q7 + vdw2m1;
                let q7pvdw2m1 = q7 + vdw2m1;
                // RDKit❗✔️:       const double t = vdw1 / (q + vdw1 - 1.0);
                let t = vdw1 / (q + vdw1 - 1.0);
                // RDKit❗✔️:       const double t2 = t * t;
                let t2 = t * t;
                // RDKit❗✔️:       const double t7 = t2 * t2 * t2 * t;
                let t7 = t2 * t2 * t2 * t;
                // RDKit❗✔️:       const double dE_dr = d_wellDepth / d_R_ij_star * t7 *
                // RDKit❗✔️:                            (-vdw2t7 * q6 / (q7pvdw2m1 * q7pvdw2m1) +
                // RDKit❗✔️:                             ((-vdw2t7 / q7pvdw2m1 + 14.0) / (q + vdw1m1)));
                let de_dr = well_depth / r_ij_star
                    * t7
                    * (-vdw2t7 * q6 / (q7pvdw2m1 * q7pvdw2m1)
                        + ((-vdw2t7 / q7pvdw2m1 + 14.0) / (q + vdw1m1)));
                // RDKit❗✔️:       vdwGrad = dE_dr / dist;
                vdw_grad = de_dr / dist;
                // RDKit❗✔️:     }
            }
            // RDKit❗✔️:     if (d_contribTypes[pairIdx] & ContribType::ELECTROSTATIC) {
            if self.contrib_types[pair_idx] & 2 != 0 {
                // RDKit❗✔️:       const double d_chargeTerm = d_chargeTerms[pairIdx];
                let charge_term = self.charge_terms[pair_idx];
                // RDKit❗✔️:       const std::uint8_t d_dielModel = d_dielModels[pairIdx];
                let diel_model = self.diel_models[pair_idx];
                // RDKit❗✔️:       const bool d_is1_4 = d_is_1_4s[pairIdx];
                let is_1_4 = self.is_1_4s[pair_idx] != 0;

                // RDKit❗✔️:       double corr_dist = dist + 0.05;
                let mut corr_dist = dist + 0.05;
                // RDKit❗✔️:       corr_dist *=
                // RDKit❗✔️:           ((d_dielModel == RDKit::MMFF::DISTANCE) ? corr_dist * corr_dist
                // RDKit❗✔️:                                                   : corr_dist);
                corr_dist *= if diel_model == 2 {
                    corr_dist * corr_dist
                } else {
                    corr_dist
                };
                // RDKit❗✔️:       const double dE_dr = -332.0716 * (double)(d_dielModel)*d_chargeTerm /
                // RDKit❗✔️:                            corr_dist * (d_is1_4 ? 0.75 : 1.0);
                let de_dr = -332.0716 * f64::from(diel_model) * charge_term / corr_dist
                    * if is_1_4 { 0.75 } else { 1.0 };
                // RDKit❗✔️:       eleGrad = dE_dr / dist;
                ele_grad = de_dr / dist;
                // RDKit❗✔️:     }
            }
            // RDKit❗✔️:     const auto dE_dr = vdwGrad + eleGrad;
            let de_dr = vdw_grad + ele_grad;
            // RDKit❗✔️:     for (unsigned int i = 0; i < 3; ++i) {
            for i in 0..3 {
                // RDKit❗✔️:       const double dGrad = dE_dr * (at1Coords[i] - at2Coords[i]);
                let d_grad = de_dr * (pos[atom1_offset + i] - pos[atom2_offset + i]);
                // RDKit❗✔️:       g1[i] += dGrad;
                // RDKit❗✔️:       g2[i] -= dGrad;
                grad[atom1_offset + i] += d_grad;
                grad[atom2_offset + i] -= d_grad;
                // RDKit❗✔️:     }
            }
            // RDKit❗✔️:   }
        }
        // RDKit❗✔️: }
        Ok(())
    }
}

impl ForceFieldContribution for NonbondedContrib {
    fn get_energy(
        &self,
        context: &mut EvaluationContext<'_>,
    ) -> Result<f64, ForceFieldKernelError> {
        self.energy(context)
    }
    fn get_grad(
        &self,
        context: &mut EvaluationContext<'_>,
        gradient: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        self.gradient(context, gradient)
    }
    fn copy(&self) -> Box<dyn ForceFieldContribution> {
        // RDKit❗✔️: return new NonbondedContrib(*this);
        // Behavior: copy all packed parameter arrays; kernel supplies copied owner context.
        // Complexity: eight vector copies, linear in terms, matching native array storage.
        Box::new(self.clone())
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::mmff::mol_properties::MmffVdwRijstarEps;
    use crate::mmff::mol_properties::{MMFF_DIELECTRIC_CONSTANT, MMFF_DIELECTRIC_DISTANCE};

    const TEST_TOLERANCE: f64 = 1.0e-12;

    fn force_field_with_positions(rows: &mut [[f64; 3]]) -> ForceField<'_> {
        let mut force_field = ForceField::new(3);
        for (idx, row) in rows.iter_mut().enumerate() {
            *row = [idx as f64, idx as f64 + 0.25, idx as f64 + 0.5];
            force_field.positions_mut().push(row.as_mut_slice());
        }
        force_field.initialize().unwrap();
        force_field
    }
    fn energy(contrib: &NonbondedContrib, pos: &[f64]) -> f64 {
        let n = (pos.len() / 3) as u32;
        let mut cache = vec![-1.; (n * (n + 1) / 2) as usize];
        contrib
            .energy(&mut EvaluationContext::for_test(pos, &mut cache, n))
            .unwrap()
    }
    fn gradient(contrib: &NonbondedContrib, pos: &[f64], grad: &mut [f64]) {
        let n = (pos.len() / 3) as u32;
        let mut cache = vec![-1.; (n * (n + 1) / 2) as usize];
        contrib
            .gradient(&mut EvaluationContext::for_test(pos, &mut cache, n), grad)
            .unwrap()
    }
    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= TEST_TOLERANCE,
            "expected {expected:.16}, got {actual:.16}"
        );
    }

    fn assert_slice_close(actual: &[f64], expected: &[f64]) {
        assert_eq!(actual.len(), expected.len());
        for (idx, (actual, expected)) in actual.iter().zip(expected.iter()).enumerate() {
            assert!(
                (actual - expected).abs() <= TEST_TOLERANCE,
                "idx={idx}, expected {expected:.16}, got {actual:.16}"
            );
        }
    }

    fn source_vdw_grad_component(r_ij_star: f64, well_depth: f64, dist: f64) -> f64 {
        let vdw1 = 1.07_f64;
        let vdw1m1 = vdw1 - 1.0;
        let vdw2 = 1.12_f64;
        let vdw2m1 = vdw2 - 1.0;
        let vdw2t7 = vdw2 * 7.0;
        let q = dist / r_ij_star;
        let q2 = q * q;
        let q6 = q2 * q2 * q2;
        let q7 = q6 * q;
        let q7pvdw2m1 = q7 + vdw2m1;
        let t = vdw1 / (q + vdw1 - 1.0);
        let t2 = t * t;
        let t7 = t2 * t2 * t2 * t;
        let de_dr = well_depth / r_ij_star
            * t7
            * (-vdw2t7 * q6 / (q7pvdw2m1 * q7pvdw2m1)
                + ((-vdw2t7 / q7pvdw2m1 + 14.0) / (q + vdw1m1)));
        de_dr / dist
    }

    fn source_ele_grad_component(charge_term: f64, diel_model: u8, is_1_4: bool, dist: f64) -> f64 {
        let mut corr_dist = dist + 0.05;
        corr_dist *= if diel_model == MMFF_DIELECTRIC_DISTANCE {
            corr_dist * corr_dist
        } else {
            corr_dist
        };
        let de_dr = -332.0716 * f64::from(diel_model) * charge_term / corr_dist
            * if is_1_4 { 0.75 } else { 1.0 };
        de_dr / dist
    }

    #[test]
    fn mmff_nonbondedcontrib_constructor_evaluates_with_kernel_owner_context() {
        let force_field = ForceField::new(3);

        let contrib = NonbondedContrib::new(&force_field);

        assert_eq!(energy(&contrib, &[0.; 6]), 0.);
        assert_eq!(energy(&contrib, &[0.; 9]), 0.);
    }

    #[test]
    fn mmff_nonbondedcontrib_constructor_initializes_no_terms() {
        let force_field = ForceField::new(3);

        let contrib = NonbondedContrib::new(&force_field);

        assert_eq!(contrib.len(), 0);
        assert!(contrib.is_empty());
    }

    #[test]
    fn mmff_nonbondedcontrib_constructor_accepts_empty_force_field_like_rdkit() {
        let force_field = ForceField::new(3);

        let contrib = NonbondedContrib::new(&force_field);

        assert_eq!(force_field.positions().len(), 0);
        assert!(contrib.is_empty());
    }

    #[test]
    fn mmff_nonbondedcontrib_add_term_pushes_vdw_and_charge_fields() {
        let mut rows = vec![[0.; 3]; 3];
        let force_field = force_field_with_positions(&mut rows);
        let mut contrib = NonbondedContrib::new(&force_field);
        let params = MmffVdwRijstarEps {
            r_ij_star_unscaled: 3.2,
            epsilon_unscaled: 0.15,
            r_ij_star: 3.0,
            epsilon: 0.12,
        };

        contrib
            .add_term(
                force_field.positions(),
                0,
                1,
                Some(&params),
                true,
                -0.25,
                MMFF_DIELECTRIC_DISTANCE,
                true,
            )
            .unwrap();

        assert_eq!(contrib.atom1_indices(), &[0]);
        assert_eq!(contrib.atom2_indices(), &[1]);
        assert_eq!(contrib.contrib_types(), &[3]);
        assert_eq!(contrib.r_ij_stars(), &[3.0]);
        assert_eq!(contrib.well_depths(), &[0.12]);
        assert_eq!(contrib.charge_terms(), &[-0.25]);
        assert_eq!(contrib.diel_models(), &[MMFF_DIELECTRIC_DISTANCE]);
        assert_eq!(contrib.is_1_4s(), &[1]);
    }

    #[test]
    fn mmff_nonbondedcontrib_add_term_stores_vdw_only_term() {
        let mut rows = vec![[0.; 3]; 3];
        let force_field = force_field_with_positions(&mut rows);
        let mut contrib = NonbondedContrib::new(&force_field);
        let params = MmffVdwRijstarEps {
            r_ij_star_unscaled: 3.2,
            epsilon_unscaled: 0.15,
            r_ij_star: 3.0,
            epsilon: 0.12,
        };

        contrib
            .add_term(
                force_field.positions(),
                0,
                1,
                Some(&params),
                false,
                0.75,
                MMFF_DIELECTRIC_DISTANCE,
                true,
            )
            .unwrap();

        assert_eq!(contrib.contrib_types(), &[1]);
        assert_eq!(contrib.r_ij_stars(), &[3.0]);
        assert_eq!(contrib.well_depths(), &[0.12]);
        assert_eq!(contrib.charge_terms(), &[0.0]);
        assert_eq!(contrib.diel_models(), &[0]);
        assert_eq!(contrib.is_1_4s(), &[0]);
    }

    #[test]
    fn mmff_nonbondedcontrib_add_term_stores_charge_only_term_without_vdw_params() {
        let mut rows = vec![[0.; 3]; 3];
        let force_field = force_field_with_positions(&mut rows);
        let mut contrib = NonbondedContrib::new(&force_field);

        contrib
            .add_term(
                force_field.positions(),
                0,
                1,
                None,
                true,
                -0.25,
                MMFF_DIELECTRIC_DISTANCE,
                true,
            )
            .unwrap();

        assert_eq!(contrib.contrib_types(), &[2]);
        assert_eq!(contrib.r_ij_stars(), &[0.0]);
        assert_eq!(contrib.well_depths(), &[0.0]);
        assert_eq!(contrib.charge_terms(), &[-0.25]);
        assert_eq!(contrib.diel_models(), &[MMFF_DIELECTRIC_DISTANCE]);
        assert_eq!(contrib.is_1_4s(), &[1]);
    }

    #[test]
    fn mmff_nonbondedcontrib_add_term_noops_without_vdw_or_charge() {
        let mut rows = vec![[0.; 3]; 3];
        let force_field = force_field_with_positions(&mut rows);
        let mut contrib = NonbondedContrib::new(&force_field);

        contrib
            .add_term(
                force_field.positions(),
                0,
                1,
                None,
                false,
                -0.25,
                MMFF_DIELECTRIC_DISTANCE,
                true,
            )
            .unwrap();

        assert!(contrib.is_empty());
    }

    #[test]
    fn mmff_nonbondedcontrib_add_term_appends_mixed_terms() {
        let mut rows = vec![[0.; 3]; 4];
        let force_field = force_field_with_positions(&mut rows);
        let mut contrib = NonbondedContrib::new(&force_field);
        let params = MmffVdwRijstarEps {
            r_ij_star_unscaled: 3.2,
            epsilon_unscaled: 0.15,
            r_ij_star: 3.0,
            epsilon: 0.12,
        };

        contrib
            .add_term(
                force_field.positions(),
                0,
                1,
                Some(&params),
                true,
                -0.25,
                MMFF_DIELECTRIC_DISTANCE,
                true,
            )
            .unwrap();
        contrib
            .add_term(force_field.positions(), 2, 3, None, true, 0.50, 1, false)
            .unwrap();

        assert_eq!(contrib.atom1_indices(), &[0, 2]);
        assert_eq!(contrib.atom2_indices(), &[1, 3]);
        assert_eq!(contrib.contrib_types(), &[3, 2]);
        assert_eq!(contrib.r_ij_stars(), &[3.0, 0.0]);
        assert_eq!(contrib.well_depths(), &[0.12, 0.0]);
        assert_eq!(contrib.charge_terms(), &[-0.25, 0.50]);
        assert_eq!(contrib.diel_models(), &[MMFF_DIELECTRIC_DISTANCE, 1]);
        assert_eq!(contrib.is_1_4s(), &[1, 0]);
    }

    #[test]
    #[should_panic]
    fn mmff_nonbondedcontrib_add_term_rejects_first_index_out_of_range() {
        let mut rows = vec![[0.; 3]; 2];
        let force_field = force_field_with_positions(&mut rows);
        let mut contrib = NonbondedContrib::new(&force_field);
        let params = MmffVdwRijstarEps {
            r_ij_star_unscaled: 3.2,
            epsilon_unscaled: 0.15,
            r_ij_star: 3.0,
            epsilon: 0.12,
        };

        contrib
            .add_term(
                force_field.positions(),
                2,
                1,
                Some(&params),
                true,
                -0.25,
                MMFF_DIELECTRIC_DISTANCE,
                true,
            )
            .unwrap();
    }

    #[test]
    #[should_panic]
    fn mmff_nonbondedcontrib_add_term_rejects_second_index_out_of_range() {
        let mut rows = vec![[0.; 3]; 2];
        let force_field = force_field_with_positions(&mut rows);
        let mut contrib = NonbondedContrib::new(&force_field);
        let params = MmffVdwRijstarEps {
            r_ij_star_unscaled: 3.2,
            epsilon_unscaled: 0.15,
            r_ij_star: 3.0,
            epsilon: 0.12,
        };

        contrib
            .add_term(
                force_field.positions(),
                0,
                2,
                Some(&params),
                true,
                -0.25,
                MMFF_DIELECTRIC_DISTANCE,
                true,
            )
            .unwrap();
    }

    #[test]
    fn mmff_nonbondedcontrib_get_energy_returns_zero_without_terms() {
        let mut rows = vec![[0.; 3]; 2];
        let force_field = force_field_with_positions(&mut rows);
        let contrib = NonbondedContrib::new(&force_field);
        let pos = [0.0, 0.25, 0.5, 1.0, 1.25, 1.5];

        assert_close(energy(&contrib, &pos), 0.0);
    }

    #[test]
    fn mmff_nonbondedcontrib_get_energy_accumulates_vdw_only_term() {
        let mut rows = vec![[0.; 3]; 2];
        let force_field = force_field_with_positions(&mut rows);
        let mut contrib = NonbondedContrib::new(&force_field);
        let params = MmffVdwRijstarEps {
            r_ij_star_unscaled: 3.2,
            epsilon_unscaled: 0.15,
            r_ij_star: 3.0,
            epsilon: 0.12,
        };
        let pos = [0.0, 0.0, 0.0, 1.5, 0.0, 0.0];

        contrib
            .add_term(
                force_field.positions(),
                0,
                1,
                Some(&params),
                false,
                0.0,
                MMFF_DIELECTRIC_DISTANCE,
                false,
            )
            .unwrap();

        let dist = 1.5_f64;
        let vdw1 = 1.07_f64;
        let vdw1m1 = vdw1 - 1.0;
        let vdw2 = 1.12_f64;
        let vdw2m1 = vdw2 - 1.0;
        let dist2 = dist * dist;
        let dist7 = dist2 * dist2 * dist2 * dist;
        let a_term = vdw1 * 3.0 / (dist + vdw1m1 * 3.0);
        let a_term2 = a_term * a_term;
        let a_term7 = a_term2 * a_term2 * a_term2 * a_term;
        let r_star2 = 3.0 * 3.0;
        let r_star7 = r_star2 * r_star2 * r_star2 * 3.0;
        let b_term = vdw2 * r_star7 / (dist7 + vdw2m1 * r_star7) - 2.0;
        let expected = 0.12 * a_term7 * b_term;

        assert_close(energy(&contrib, &pos), expected);
    }

    #[test]
    fn mmff_nonbondedcontrib_get_energy_accumulates_charge_only_constant_term() {
        let mut rows = vec![[0.; 3]; 2];
        let force_field = force_field_with_positions(&mut rows);
        let mut contrib = NonbondedContrib::new(&force_field);
        let pos = [0.0, 0.0, 0.0, 1.5, 0.0, 0.0];

        contrib
            .add_term(
                force_field.positions(),
                0,
                1,
                None,
                true,
                -0.5,
                MMFF_DIELECTRIC_CONSTANT,
                false,
            )
            .unwrap();

        let expected = 332.0716 * -0.5 / (1.5 + 0.05);

        assert_close(energy(&contrib, &pos), expected);
    }

    #[test]
    fn mmff_nonbondedcontrib_get_energy_accumulates_charge_only_distance_scaled_one_four_term() {
        let mut rows = vec![[0.; 3]; 2];
        let force_field = force_field_with_positions(&mut rows);
        let mut contrib = NonbondedContrib::new(&force_field);
        let pos = [0.0, 0.0, 0.0, 1.5, 0.0, 0.0];

        contrib
            .add_term(
                force_field.positions(),
                0,
                1,
                None,
                true,
                -0.5,
                MMFF_DIELECTRIC_DISTANCE,
                true,
            )
            .unwrap();

        let corr_dist = (1.5_f64 + 0.05).powi(2);
        let expected = 332.0716 * -0.5 / corr_dist * 0.75;

        assert_close(energy(&contrib, &pos), expected);
    }

    #[test]
    fn mmff_nonbondedcontrib_get_energy_sums_vdw_and_electrostatic_contributions() {
        let mut rows = vec![[0.; 3]; 3];
        let force_field = force_field_with_positions(&mut rows);
        let mut contrib = NonbondedContrib::new(&force_field);
        let params = MmffVdwRijstarEps {
            r_ij_star_unscaled: 3.2,
            epsilon_unscaled: 0.15,
            r_ij_star: 3.0,
            epsilon: 0.12,
        };
        let pos = [0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 0.0, 2.0, 0.0];

        contrib
            .add_term(
                force_field.positions(),
                0,
                1,
                Some(&params),
                true,
                -0.25,
                MMFF_DIELECTRIC_DISTANCE,
                true,
            )
            .unwrap();
        contrib
            .add_term(
                force_field.positions(),
                1,
                2,
                None,
                true,
                0.4,
                MMFF_DIELECTRIC_CONSTANT,
                false,
            )
            .unwrap();

        let dist01 = 1.5_f64;
        let vdw1 = 1.07_f64;
        let vdw1m1 = vdw1 - 1.0;
        let vdw2 = 1.12_f64;
        let vdw2m1 = vdw2 - 1.0;
        let dist2 = dist01 * dist01;
        let dist7 = dist2 * dist2 * dist2 * dist01;
        let a_term = vdw1 * 3.0 / (dist01 + vdw1m1 * 3.0);
        let a_term2 = a_term * a_term;
        let a_term7 = a_term2 * a_term2 * a_term2 * a_term;
        let r_star2 = 3.0 * 3.0;
        let r_star7 = r_star2 * r_star2 * r_star2 * 3.0;
        let b_term = vdw2 * r_star7 / (dist7 + vdw2m1 * r_star7) - 2.0;
        let vdw_energy = 0.12 * a_term7 * b_term;
        let ele01 = 332.0716 * -0.25 / (1.5_f64 + 0.05).powi(2) * 0.75;
        let ele12 = 332.0716 * 0.4 / (2.5_f64 + 0.05);
        let expected = vdw_energy + ele01 + ele12;

        assert_close(energy(&contrib, &pos), expected);
    }

    #[test]
    fn mmff_nonbondedcontrib_get_grad_accumulates_vdw_only_term() {
        let mut rows = vec![[0.; 3]; 2];
        let force_field = force_field_with_positions(&mut rows);
        let mut contrib = NonbondedContrib::new(&force_field);
        let params = MmffVdwRijstarEps {
            r_ij_star_unscaled: 3.2,
            epsilon_unscaled: 0.15,
            r_ij_star: 3.0,
            epsilon: 0.12,
        };
        let pos = [0.0, 0.0, 0.0, 1.5, 0.0, 0.0];
        let mut grad = vec![0.0; 6];

        contrib
            .add_term(
                force_field.positions(),
                0,
                1,
                Some(&params),
                false,
                0.0,
                MMFF_DIELECTRIC_DISTANCE,
                false,
            )
            .unwrap();
        gradient(&contrib, &pos, &mut grad);

        let de_dr = source_vdw_grad_component(3.0, 0.12, 1.5);
        let expected = vec![-1.5 * de_dr, 0.0, 0.0, 1.5 * de_dr, 0.0, 0.0];
        assert_slice_close(&grad, &expected);
    }

    #[test]
    fn mmff_nonbondedcontrib_get_grad_accumulates_charge_only_constant_term() {
        let mut rows = vec![[0.; 3]; 2];
        let force_field = force_field_with_positions(&mut rows);
        let mut contrib = NonbondedContrib::new(&force_field);
        let pos = [0.0, 0.0, 0.0, 1.5, 0.0, 0.0];
        let mut grad = vec![0.0; 6];

        contrib
            .add_term(
                force_field.positions(),
                0,
                1,
                None,
                true,
                -0.5,
                MMFF_DIELECTRIC_CONSTANT,
                false,
            )
            .unwrap();
        gradient(&contrib, &pos, &mut grad);

        let de_dr = source_ele_grad_component(-0.5, MMFF_DIELECTRIC_CONSTANT, false, 1.5);
        let expected = vec![-1.5 * de_dr, 0.0, 0.0, 1.5 * de_dr, 0.0, 0.0];
        assert_slice_close(&grad, &expected);
    }

    #[test]
    fn mmff_nonbondedcontrib_get_grad_accumulates_combined_term() {
        let mut rows = vec![[0.; 3]; 2];
        let force_field = force_field_with_positions(&mut rows);
        let mut contrib = NonbondedContrib::new(&force_field);
        let params = MmffVdwRijstarEps {
            r_ij_star_unscaled: 3.2,
            epsilon_unscaled: 0.15,
            r_ij_star: 3.0,
            epsilon: 0.12,
        };
        let pos = [0.0, 0.0, 0.0, 1.5, 0.0, 0.0];
        let mut grad = vec![0.0; 6];

        contrib
            .add_term(
                force_field.positions(),
                0,
                1,
                Some(&params),
                true,
                -0.25,
                MMFF_DIELECTRIC_DISTANCE,
                true,
            )
            .unwrap();
        gradient(&contrib, &pos, &mut grad);

        let de_dr = source_vdw_grad_component(3.0, 0.12, 1.5)
            + source_ele_grad_component(-0.25, MMFF_DIELECTRIC_DISTANCE, true, 1.5);
        let expected = vec![-1.5 * de_dr, 0.0, 0.0, 1.5 * de_dr, 0.0, 0.0];
        assert_slice_close(&grad, &expected);
    }

    #[test]
    fn mmff_nonbondedcontrib_get_grad_uses_zero_distance_fallbacks() {
        let mut rows = vec![[0.; 3]; 2];
        let force_field = force_field_with_positions(&mut rows);
        let mut contrib = NonbondedContrib::new(&force_field);
        let params = MmffVdwRijstarEps {
            r_ij_star_unscaled: 3.2,
            epsilon_unscaled: 0.15,
            r_ij_star: 3.0,
            epsilon: 0.12,
        };
        let pos = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let mut grad = vec![0.0; 6];

        contrib
            .add_term(
                force_field.positions(),
                0,
                1,
                Some(&params),
                true,
                -0.25,
                MMFF_DIELECTRIC_DISTANCE,
                true,
            )
            .unwrap();
        gradient(&contrib, &pos, &mut grad);

        let expected = vec![0.05, 0.05, 0.05, -0.05, -0.05, -0.05];
        assert_slice_close(&grad, &expected);
    }

    #[test]
    fn mmff_nonbondedcontrib_get_grad_zero_distance_returns_before_later_pairs() {
        let mut rows = vec![[0.; 3]; 3];
        let force_field = force_field_with_positions(&mut rows);
        let mut contrib = NonbondedContrib::new(&force_field);
        let params = MmffVdwRijstarEps {
            r_ij_star_unscaled: 3.2,
            epsilon_unscaled: 0.15,
            r_ij_star: 3.0,
            epsilon: 0.12,
        };
        let pos = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.5, 0.0, 0.0];
        let mut grad = vec![0.0; 9];

        contrib
            .add_term(
                force_field.positions(),
                0,
                1,
                Some(&params),
                true,
                -0.25,
                MMFF_DIELECTRIC_DISTANCE,
                true,
            )
            .unwrap();
        contrib
            .add_term(
                force_field.positions(),
                1,
                2,
                None,
                true,
                0.4,
                MMFF_DIELECTRIC_CONSTANT,
                false,
            )
            .unwrap();
        gradient(&contrib, &pos, &mut grad);

        let expected = vec![0.05, 0.05, 0.05, -0.05, -0.05, -0.05, 0.0, 0.0, 0.0];
        assert_slice_close(&grad, &expected);
    }
}
