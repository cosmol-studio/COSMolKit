//! Source-backed RDKit MMFF nonbonded contribution.

use crate::chemistry::forcefield::core::{ForceField, ForceFieldContrib};

use super::params::MmffVdwRijstarEps;

fn calc_vdw_energy(dist: f64, r_star_ij: f64, well_depth: f64) -> f64 {
    // BEGIN RDKIT CPP FUNCTION ForceFields::MMFF::Utils::calcVdWEnergy (Nonbonded.cpp:43-57)
    // RDKit✔️✔️: double calcVdWEnergy(const double dist, const double R_star_ij,
    // RDKit✔️✔️:                      const double wellDepth) {
    // RDKit✔️✔️:   double const vdw1 = 1.07;
    let vdw1 = 1.07;
    // RDKit✔️✔️:   double const vdw1m1 = vdw1 - 1.0;
    let vdw1m1 = vdw1 - 1.0;
    // RDKit✔️✔️:   double const vdw2 = 1.12;
    let vdw2 = 1.12;
    // RDKit✔️✔️:   double const vdw2m1 = vdw2 - 1.0;
    let vdw2m1 = vdw2 - 1.0;
    // RDKit✔️✔️:   double dist2 = dist * dist;
    let dist2 = dist * dist;
    // RDKit✔️✔️:   double dist7 = dist2 * dist2 * dist2 * dist;
    let dist7 = dist2 * dist2 * dist2 * dist;
    // RDKit✔️✔️:   double aTerm = vdw1 * R_star_ij / (dist + vdw1m1 * R_star_ij);
    let a_term = vdw1 * r_star_ij / (dist + vdw1m1 * r_star_ij);
    // RDKit✔️✔️:   double aTerm2 = aTerm * aTerm;
    let a_term2 = a_term * a_term;
    // RDKit✔️✔️:   double aTerm7 = aTerm2 * aTerm2 * aTerm2 * aTerm;
    let a_term7 = a_term2 * a_term2 * a_term2 * a_term;
    // RDKit✔️✔️:   double R_star_ij2 = R_star_ij * R_star_ij;
    let r_star_ij2 = r_star_ij * r_star_ij;
    // RDKit✔️✔️:   double R_star_ij7 = R_star_ij2 * R_star_ij2 * R_star_ij2 * R_star_ij;
    let r_star_ij7 = r_star_ij2 * r_star_ij2 * r_star_ij2 * r_star_ij;
    // RDKit✔️✔️:   double bTerm = vdw2 * R_star_ij7 / (dist7 + vdw2m1 * R_star_ij7) - 2.0;
    let b_term = vdw2 * r_star_ij7 / (dist7 + vdw2m1 * r_star_ij7) - 2.0;
    // RDKit✔️✔️:   double res = wellDepth * aTerm7 * bTerm;
    let res = well_depth * a_term7 * b_term;

    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    res
}

fn calc_ele_energy(
    idx1: usize,
    idx2: usize,
    dist: f64,
    charge_term: f64,
    diel_model: u8,
    is_1_4: bool,
) -> f64 {
    // BEGIN RDKIT CPP FUNCTION ForceFields::MMFF::Utils::calcEleEnergy (Nonbonded.cpp:76-85)
    // RDKit✔️✔️: double calcEleEnergy(unsigned int, unsigned int, double dist, double chargeTerm,
    // RDKit✔️✔️:                      std::uint8_t dielModel, bool is1_4) {
    let _ = (idx1, idx2);
    // RDKit✔️✔️:   double corr_dist = dist + 0.05;
    let mut corr_dist = dist + 0.05;
    // RDKit✔️✔️:   double const diel = 332.0716;
    let diel = 332.0716;
    // RDKit✔️✔️:   double const sc1_4 = 0.75;
    let sc1_4 = 0.75;
    // RDKit✔️✔️:   if (dielModel == RDKit::MMFF::DISTANCE) {
    if diel_model == 2 {
        // RDKit✔️✔️:     corr_dist *= corr_dist;
        corr_dist *= corr_dist;
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return (diel * chargeTerm / corr_dist * (is1_4 ? sc1_4 : 1.0));
    // RDKit✔️✔️: }
    diel * charge_term / corr_dist * if is_1_4 { sc1_4 } else { 1.0 }
}

#[derive(Clone, Debug)]
pub struct NonbondedContrib {
    owner: *const ForceField,
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
    pub fn new(owner: &ForceField) -> Self {
        // BEGIN RDKIT CPP CONSTRUCTOR ForceFields::MMFF::NonbondedContrib::NonbondedContrib (Nonbonded.cpp:87-90)
        // RDKit✔️✔️: NonbondedContrib::NonbondedContrib(ForceField *owner) {
        // RDKit✔️✔️:   PRECONDITION(owner, "bad owner");
        // Rust references reproduce RDKit's non-null owner precondition.
        // RDKit✔️✔️:   dp_forceField = owner;
        let owner = owner as *const ForceField;
        // RDKit✔️✔️: }
        Self {
            owner,
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
    pub fn owner(&self) -> *const ForceField {
        self.owner
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.atom1_indices.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.atom1_indices.is_empty()
            && self.atom2_indices.is_empty()
            && self.contrib_types.is_empty()
            && self.r_ij_stars.is_empty()
            && self.well_depths.is_empty()
            && self.charge_terms.is_empty()
            && self.is_1_4s.is_empty()
            && self.diel_models.is_empty()
    }

    pub fn add_term(
        &mut self,
        idx1: usize,
        idx2: usize,
        mmff_vdw_constants: Option<&MmffVdwRijstarEps>,
        include_charge: bool,
        charge_term: f64,
        diel_model: u8,
        is_1_4: bool,
    ) {
        // BEGIN RDKIT CPP METHOD ForceFields::MMFF::NonbondedContrib::addTerm (Nonbonded.cpp:268-290)
        // RDKit✔️✔️: void NonbondedContrib::addTerm(unsigned int idx1, unsigned int idx2,
        // RDKit✔️✔️:                              const MMFFVdWRijstarEps *mmffVdWConstants, bool includeCharge,
        // RDKit✔️✔️:                              double chargeTerm, std::uint8_t dielModel, bool is1_4) {
        // RDKit✔️✔️:   if (!mmffVdWConstants && !includeCharge) {
        if mmff_vdw_constants.is_none() && !include_charge {
            // RDKit✔️✔️:     return;
            return;
        }
        // RDKit✔️✔️:   }
        let force_field = self.force_field();
        // RDKit✔️✔️:   URANGE_CHECK(idx1, dp_forceField->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx2, dp_forceField->positions().size());
        assert!(idx1 < force_field.positions().len());
        assert!(idx2 < force_field.positions().len());
        // RDKit✔️✔️:   d_at1Idxs.push_back(idx1);
        // RDKit✔️✔️:   d_at2Idxs.push_back(idx2);
        self.atom1_indices.push(idx1 as i16);
        self.atom2_indices.push(idx2 as i16);
        // RDKit✔️✔️:   d_contribTypes.push_back(0);
        self.contrib_types.push(0);
        // RDKit✔️✔️:   if (mmffVdWConstants) {
        if let Some(mmff_vdw_constants) = mmff_vdw_constants {
            // RDKit✔️✔️:     d_contribTypes.back() |= ContribType::VDW;
            *self.contrib_types.last_mut().expect("term was just pushed") |= 1;
            // RDKit✔️✔️:     d_R_ij_stars.push_back(mmffVdWConstants->R_ij_star);
            self.r_ij_stars.push(mmff_vdw_constants.r_ij_star);
            // RDKit✔️✔️:     d_wellDepths.push_back(mmffVdWConstants->epsilon);
            self.well_depths.push(mmff_vdw_constants.epsilon);
            // RDKit✔️✔️:   } else {
        } else {
            // RDKit✔️✔️:     d_R_ij_stars.push_back(0.0);
            self.r_ij_stars.push(0.0);
            // RDKit✔️✔️:     d_wellDepths.push_back(0.0);
            self.well_depths.push(0.0);
            // RDKit✔️✔️:   }
        }
        // RDKit✔️✔️:   if (includeCharge) {
        if include_charge {
            // RDKit✔️✔️:     d_contribTypes.back() |= ContribType::ELECTROSTATIC;
            *self.contrib_types.last_mut().expect("term was just pushed") |= 2;
            // RDKit✔️✔️:     d_chargeTerms.push_back(chargeTerm);
            self.charge_terms.push(charge_term);
            // RDKit✔️✔️:     d_dielModels.push_back(dielModel);
            self.diel_models.push(diel_model);
            // RDKit✔️✔️:     d_is_1_4s.push_back(is1_4);
            self.is_1_4s.push(u8::from(is_1_4));
            // RDKit✔️✔️:   } else {
        } else {
            // RDKit✔️✔️:     d_chargeTerms.push_back(0.0);
            self.charge_terms.push(0.0);
            // RDKit✔️✔️:     d_dielModels.push_back(0);
            self.diel_models.push(0);
            // RDKit✔️✔️:     d_is_1_4s.push_back(false);
            self.is_1_4s.push(0);
            // RDKit✔️✔️:   }
        }
        // RDKit✔️✔️: }
    }

    #[must_use]
    pub fn atom1_indices(&self) -> &[i16] {
        &self.atom1_indices
    }

    #[must_use]
    pub fn atom2_indices(&self) -> &[i16] {
        &self.atom2_indices
    }

    #[must_use]
    pub fn contrib_types(&self) -> &[u8] {
        &self.contrib_types
    }

    #[must_use]
    pub fn r_ij_stars(&self) -> &[f64] {
        &self.r_ij_stars
    }

    #[must_use]
    pub fn well_depths(&self) -> &[f64] {
        &self.well_depths
    }

    #[must_use]
    pub fn charge_terms(&self) -> &[f64] {
        &self.charge_terms
    }

    #[must_use]
    pub fn is_1_4s(&self) -> &[u8] {
        &self.is_1_4s
    }

    #[must_use]
    pub fn diel_models(&self) -> &[u8] {
        &self.diel_models
    }

    #[must_use]
    pub fn get_energy(&self, pos: &[f64]) -> f64 {
        // BEGIN RDKIT CPP METHOD ForceFields::MMFF::NonbondedContrib::getEnergy (Nonbonded.cpp:295-319)
        // RDKit✔️✔️: double NonbondedContrib::getEnergy(double *pos) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        let force_field = self.force_field();
        // Rust slices reproduce RDKit's non-null pos precondition.
        // RDKit✔️✔️:   double energySum = 0.0;
        let mut energy_sum = 0.0;

        // RDKit✔️✔️:   const int numPairs = d_at1Idxs.size();
        let num_pairs = self.atom1_indices.len();
        // RDKit✔️✔️:   for (int i = 0; i < numPairs; ++i) {
        for i in 0..num_pairs {
            // RDKit✔️✔️:     unsigned int d_at1Idx = d_at1Idxs[i];
            // RDKit✔️✔️:     unsigned int d_at2Idx = d_at2Idxs[i];
            let atom1_idx = self.atom1_indices[i] as usize;
            let atom2_idx = self.atom2_indices[i] as usize;
            // RDKit✔️✔️:     double dist = dp_forceField->distance(d_at1Idx, d_at2Idx, pos);
            let dist = force_field.distance(atom1_idx, atom2_idx, Some(pos));

            // RDKit✔️✔️:     if (d_contribTypes[i] & ContribType::VDW) {
            if self.contrib_types[i] & 1 != 0 {
                // RDKit✔️✔️:       const auto res =
                // RDKit✔️✔️:           Utils::calcVdWEnergy(dist, d_R_ij_stars[i], d_wellDepths[i]);
                let res = calc_vdw_energy(dist, self.r_ij_stars[i], self.well_depths[i]);
                // RDKit✔️✔️:       energySum += res;
                energy_sum += res;
                // RDKit✔️✔️:     }
            }
            // RDKit✔️✔️:     if (d_contribTypes[i] & ContribType::ELECTROSTATIC) {
            if self.contrib_types[i] & 2 != 0 {
                // RDKit✔️✔️:       const double d_chargeTerm = d_chargeTerms[i];
                let charge_term = self.charge_terms[i];
                // RDKit✔️✔️:       const std::uint8_t d_dielModel = d_dielModels[i];
                let diel_model = self.diel_models[i];
                // RDKit✔️✔️:       const bool d_is1_4 = d_is_1_4s[i];
                let is_1_4 = self.is_1_4s[i] != 0;
                // RDKit✔️✔️:       const auto res = Utils::calcEleEnergy(d_at1Idx, d_at2Idx, dist,
                // RDKit✔️✔️:                                             d_chargeTerm, d_dielModel, d_is1_4);
                let res =
                    calc_ele_energy(atom1_idx, atom2_idx, dist, charge_term, diel_model, is_1_4);
                // RDKit✔️✔️:       energySum += res;
                energy_sum += res;
                // RDKit✔️✔️:     }
            }
            // RDKit✔️✔️:   }
        }
        // RDKit✔️✔️:   return energySum;
        // RDKit✔️✔️: }
        energy_sum
    }

    pub fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
        // BEGIN RDKIT CPP METHOD ForceFields::MMFF::NonbondedContrib::getGrad (Nonbonded.cpp:321-380)
        // RDKit✔️✔️: void NonbondedContrib::getGrad(double *pos, double *grad) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:   PRECONDITION(grad, "bad vector");
        let force_field = self.force_field();
        // Rust slices reproduce RDKit's non-null pos and grad preconditions.

        // RDKit✔️✔️:   constexpr double vdw1 = 1.07;
        let vdw1 = 1.07;
        // RDKit✔️✔️:   constexpr double vdw1m1 = vdw1 - 1.0;
        let vdw1m1 = vdw1 - 1.0;
        // RDKit✔️✔️:   constexpr double vdw2 = 1.12;
        let vdw2 = 1.12;
        // RDKit✔️✔️:   constexpr double vdw2m1 = vdw2 - 1.0;
        let vdw2m1 = vdw2 - 1.0;
        // RDKit✔️✔️:   constexpr double vdw2t7 = vdw2 * 7.0;
        let vdw2t7 = vdw2 * 7.0;

        // RDKit✔️✔️:   const int numPairs = d_at1Idxs.size();
        let num_pairs = self.atom1_indices.len();
        // RDKit✔️✔️:   for (int pairIdx = 0; pairIdx < numPairs; ++pairIdx) {
        for pair_idx in 0..num_pairs {
            // RDKit✔️✔️:     const int d_at1Idx = d_at1Idxs[pairIdx];
            // RDKit✔️✔️:     const int d_at2Idx = d_at2Idxs[pairIdx];
            let atom1_idx = self.atom1_indices[pair_idx] as usize;
            let atom2_idx = self.atom2_indices[pair_idx] as usize;
            // RDKit✔️✔️:     const double dist = dp_forceField->distance(d_at1Idx, d_at2Idx, pos);
            let dist = force_field.distance(atom1_idx, atom2_idx, Some(pos));
            // RDKit✔️✔️:     const double *at1Coords = &(pos[3 * d_at1Idx]);
            // RDKit✔️✔️:     const double *at2Coords = &(pos[3 * d_at2Idx]);
            let atom1_offset = 3 * atom1_idx;
            let atom2_offset = 3 * atom2_idx;

            // RDKit✔️✔️:     double vdwGrad = 0.0;
            // RDKit✔️✔️:     double eleGrad = 0.0;
            let mut vdw_grad = 0.0;
            let mut ele_grad = 0.0;
            // RDKit✔️✔️:     if (dist <= 0.0) {
            if dist <= 0.0 {
                // RDKit✔️✔️:       if (d_contribTypes[pairIdx] & ContribType::VDW) {
                if self.contrib_types[pair_idx] & 1 != 0 {
                    // RDKit✔️✔️:         const double d_R_ij_star = d_R_ij_stars[pairIdx];
                    let r_ij_star = self.r_ij_stars[pair_idx];
                    // RDKit✔️✔️:         for (unsigned int i = 0; i < 3; ++i) {
                    for i in 0..3 {
                        // RDKit✔️✔️:           g1[i] += d_R_ij_star * 0.01;
                        // RDKit✔️✔️:           g2[i] -= d_R_ij_star * 0.01;
                        grad[atom1_offset + i] += r_ij_star * 0.01;
                        grad[atom2_offset + i] -= r_ij_star * 0.01;
                        // RDKit✔️✔️:         }
                    }
                    // RDKit✔️✔️:       }
                }
                // RDKit✔️✔️:       if (d_contribTypes[pairIdx] & ContribType::ELECTROSTATIC) {
                if self.contrib_types[pair_idx] & 2 != 0 {
                    // RDKit✔️✔️:         for (unsigned int i = 0; i < 3; ++i) {
                    for i in 0..3 {
                        // RDKit✔️✔️:           g1[i] += 0.02;
                        // RDKit✔️✔️:           g2[i] -= 0.02;
                        grad[atom1_offset + i] += 0.02;
                        grad[atom2_offset + i] -= 0.02;
                        // RDKit✔️✔️:         }
                    }
                    // RDKit✔️✔️:       }
                }
                // RDKit✔️✔️:       return;
                return;
                // RDKit✔️✔️:     }
            }
            // RDKit✔️✔️:     if (d_contribTypes[pairIdx] & ContribType::VDW) {
            if self.contrib_types[pair_idx] & 1 != 0 {
                // RDKit✔️✔️:       const double d_R_ij_star = d_R_ij_stars[pairIdx];
                let r_ij_star = self.r_ij_stars[pair_idx];
                // RDKit✔️✔️:       const double d_wellDepth = d_wellDepths[pairIdx];
                let well_depth = self.well_depths[pair_idx];

                // RDKit✔️✔️:       const double q = dist / d_R_ij_star;
                let q = dist / r_ij_star;
                // RDKit✔️✔️:       const double q2 = q * q;
                let q2 = q * q;
                // RDKit✔️✔️:       const double q6 = q2 * q2 * q2;
                let q6 = q2 * q2 * q2;
                // RDKit✔️✔️:       const double q7 = q6 * q;
                let q7 = q6 * q;
                // RDKit✔️✔️:       const double q7pvdw2m1 = q7 + vdw2m1;
                let q7pvdw2m1 = q7 + vdw2m1;
                // RDKit✔️✔️:       const double t = vdw1 / (q + vdw1 - 1.0);
                let t = vdw1 / (q + vdw1 - 1.0);
                // RDKit✔️✔️:       const double t2 = t * t;
                let t2 = t * t;
                // RDKit✔️✔️:       const double t7 = t2 * t2 * t2 * t;
                let t7 = t2 * t2 * t2 * t;
                // RDKit✔️✔️:       const double dE_dr = d_wellDepth / d_R_ij_star * t7 *
                // RDKit✔️✔️:                            (-vdw2t7 * q6 / (q7pvdw2m1 * q7pvdw2m1) +
                // RDKit✔️✔️:                             ((-vdw2t7 / q7pvdw2m1 + 14.0) / (q + vdw1m1)));
                let de_dr = well_depth / r_ij_star
                    * t7
                    * (-vdw2t7 * q6 / (q7pvdw2m1 * q7pvdw2m1)
                        + ((-vdw2t7 / q7pvdw2m1 + 14.0) / (q + vdw1m1)));
                // RDKit✔️✔️:       vdwGrad = dE_dr / dist;
                vdw_grad = de_dr / dist;
                // RDKit✔️✔️:     }
            }
            // RDKit✔️✔️:     if (d_contribTypes[pairIdx] & ContribType::ELECTROSTATIC) {
            if self.contrib_types[pair_idx] & 2 != 0 {
                // RDKit✔️✔️:       const double d_chargeTerm = d_chargeTerms[pairIdx];
                let charge_term = self.charge_terms[pair_idx];
                // RDKit✔️✔️:       const std::uint8_t d_dielModel = d_dielModels[pairIdx];
                let diel_model = self.diel_models[pair_idx];
                // RDKit✔️✔️:       const bool d_is1_4 = d_is_1_4s[pairIdx];
                let is_1_4 = self.is_1_4s[pair_idx] != 0;

                // RDKit✔️✔️:       double corr_dist = dist + 0.05;
                let mut corr_dist = dist + 0.05;
                // RDKit✔️✔️:       corr_dist *=
                // RDKit✔️✔️:           ((d_dielModel == RDKit::MMFF::DISTANCE) ? corr_dist * corr_dist
                // RDKit✔️✔️:                                                   : corr_dist);
                corr_dist *= if diel_model == 2 {
                    corr_dist * corr_dist
                } else {
                    corr_dist
                };
                // RDKit✔️✔️:       const double dE_dr = -332.0716 * (double)(d_dielModel)*d_chargeTerm /
                // RDKit✔️✔️:                            corr_dist * (d_is1_4 ? 0.75 : 1.0);
                let de_dr = -332.0716 * f64::from(diel_model) * charge_term / corr_dist
                    * if is_1_4 { 0.75 } else { 1.0 };
                // RDKit✔️✔️:       eleGrad = dE_dr / dist;
                ele_grad = de_dr / dist;
                // RDKit✔️✔️:     }
            }
            // RDKit✔️✔️:     const auto dE_dr = vdwGrad + eleGrad;
            let de_dr = vdw_grad + ele_grad;
            // RDKit✔️✔️:     for (unsigned int i = 0; i < 3; ++i) {
            for i in 0..3 {
                // RDKit✔️✔️:       const double dGrad = dE_dr * (at1Coords[i] - at2Coords[i]);
                let d_grad = de_dr * (pos[atom1_offset + i] - pos[atom2_offset + i]);
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

    fn force_field(&self) -> &ForceField {
        // RDKit✔️✔️: dp_forceField
        unsafe { &*self.owner }
    }
}

impl ForceFieldContrib for NonbondedContrib {
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
    use super::NonbondedContrib;
    use crate::chemistry::forcefield::core::{ForceField, ForceFieldVec3};
    use crate::chemistry::forcefield::mmff::mol_properties::{
        MMFF_DIELECTRIC_CONSTANT, MMFF_DIELECTRIC_DISTANCE,
    };
    use crate::chemistry::forcefield::mmff::params::MmffVdwRijstarEps;

    const TEST_TOLERANCE: f64 = 1.0e-12;

    fn force_field_with_positions(count: usize) -> ForceField {
        let mut force_field = ForceField::new(count);
        for idx in 0..count {
            force_field.positions_mut().push(ForceFieldVec3::new(
                idx as f64,
                idx as f64 + 0.25,
                idx as f64 + 0.5,
            ));
        }
        force_field.initialize();
        force_field
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
    fn mmff_nonbondedcontrib_constructor_stores_owner_pointer() {
        let force_field = ForceField::new(3);

        let contrib = NonbondedContrib::new(&force_field);

        assert_eq!(contrib.owner(), &force_field as *const ForceField);
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
        let force_field = force_field_with_positions(3);
        let mut contrib = NonbondedContrib::new(&force_field);
        let params = MmffVdwRijstarEps {
            r_ij_star_unscaled: 3.2,
            epsilon_unscaled: 0.15,
            r_ij_star: 3.0,
            epsilon: 0.12,
        };

        contrib.add_term(
            0,
            1,
            Some(&params),
            true,
            -0.25,
            MMFF_DIELECTRIC_DISTANCE,
            true,
        );

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
        let force_field = force_field_with_positions(3);
        let mut contrib = NonbondedContrib::new(&force_field);
        let params = MmffVdwRijstarEps {
            r_ij_star_unscaled: 3.2,
            epsilon_unscaled: 0.15,
            r_ij_star: 3.0,
            epsilon: 0.12,
        };

        contrib.add_term(
            0,
            1,
            Some(&params),
            false,
            0.75,
            MMFF_DIELECTRIC_DISTANCE,
            true,
        );

        assert_eq!(contrib.contrib_types(), &[1]);
        assert_eq!(contrib.r_ij_stars(), &[3.0]);
        assert_eq!(contrib.well_depths(), &[0.12]);
        assert_eq!(contrib.charge_terms(), &[0.0]);
        assert_eq!(contrib.diel_models(), &[0]);
        assert_eq!(contrib.is_1_4s(), &[0]);
    }

    #[test]
    fn mmff_nonbondedcontrib_add_term_stores_charge_only_term_without_vdw_params() {
        let force_field = force_field_with_positions(3);
        let mut contrib = NonbondedContrib::new(&force_field);

        contrib.add_term(0, 1, None, true, -0.25, MMFF_DIELECTRIC_DISTANCE, true);

        assert_eq!(contrib.contrib_types(), &[2]);
        assert_eq!(contrib.r_ij_stars(), &[0.0]);
        assert_eq!(contrib.well_depths(), &[0.0]);
        assert_eq!(contrib.charge_terms(), &[-0.25]);
        assert_eq!(contrib.diel_models(), &[MMFF_DIELECTRIC_DISTANCE]);
        assert_eq!(contrib.is_1_4s(), &[1]);
    }

    #[test]
    fn mmff_nonbondedcontrib_add_term_noops_without_vdw_or_charge() {
        let force_field = force_field_with_positions(3);
        let mut contrib = NonbondedContrib::new(&force_field);

        contrib.add_term(0, 1, None, false, -0.25, MMFF_DIELECTRIC_DISTANCE, true);

        assert!(contrib.is_empty());
    }

    #[test]
    fn mmff_nonbondedcontrib_add_term_appends_mixed_terms() {
        let force_field = force_field_with_positions(4);
        let mut contrib = NonbondedContrib::new(&force_field);
        let params = MmffVdwRijstarEps {
            r_ij_star_unscaled: 3.2,
            epsilon_unscaled: 0.15,
            r_ij_star: 3.0,
            epsilon: 0.12,
        };

        contrib.add_term(
            0,
            1,
            Some(&params),
            true,
            -0.25,
            MMFF_DIELECTRIC_DISTANCE,
            true,
        );
        contrib.add_term(2, 3, None, true, 0.50, 1, false);

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
        let force_field = force_field_with_positions(2);
        let mut contrib = NonbondedContrib::new(&force_field);
        let params = MmffVdwRijstarEps {
            r_ij_star_unscaled: 3.2,
            epsilon_unscaled: 0.15,
            r_ij_star: 3.0,
            epsilon: 0.12,
        };

        contrib.add_term(
            2,
            1,
            Some(&params),
            true,
            -0.25,
            MMFF_DIELECTRIC_DISTANCE,
            true,
        );
    }

    #[test]
    #[should_panic]
    fn mmff_nonbondedcontrib_add_term_rejects_second_index_out_of_range() {
        let force_field = force_field_with_positions(2);
        let mut contrib = NonbondedContrib::new(&force_field);
        let params = MmffVdwRijstarEps {
            r_ij_star_unscaled: 3.2,
            epsilon_unscaled: 0.15,
            r_ij_star: 3.0,
            epsilon: 0.12,
        };

        contrib.add_term(
            0,
            2,
            Some(&params),
            true,
            -0.25,
            MMFF_DIELECTRIC_DISTANCE,
            true,
        );
    }

    #[test]
    fn mmff_nonbondedcontrib_get_energy_returns_zero_without_terms() {
        let force_field = force_field_with_positions(2);
        let contrib = NonbondedContrib::new(&force_field);
        let pos = [0.0, 0.25, 0.5, 1.0, 1.25, 1.5];

        assert_close(contrib.get_energy(&pos), 0.0);
    }

    #[test]
    fn mmff_nonbondedcontrib_get_energy_accumulates_vdw_only_term() {
        let force_field = force_field_with_positions(2);
        let mut contrib = NonbondedContrib::new(&force_field);
        let params = MmffVdwRijstarEps {
            r_ij_star_unscaled: 3.2,
            epsilon_unscaled: 0.15,
            r_ij_star: 3.0,
            epsilon: 0.12,
        };
        let pos = [0.0, 0.0, 0.0, 1.5, 0.0, 0.0];

        contrib.add_term(
            0,
            1,
            Some(&params),
            false,
            0.0,
            MMFF_DIELECTRIC_DISTANCE,
            false,
        );

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

        assert_close(contrib.get_energy(&pos), expected);
    }

    #[test]
    fn mmff_nonbondedcontrib_get_energy_accumulates_charge_only_constant_term() {
        let force_field = force_field_with_positions(2);
        let mut contrib = NonbondedContrib::new(&force_field);
        let pos = [0.0, 0.0, 0.0, 1.5, 0.0, 0.0];

        contrib.add_term(0, 1, None, true, -0.5, MMFF_DIELECTRIC_CONSTANT, false);

        let expected = 332.0716 * -0.5 / (1.5 + 0.05);

        assert_close(contrib.get_energy(&pos), expected);
    }

    #[test]
    fn mmff_nonbondedcontrib_get_energy_accumulates_charge_only_distance_scaled_one_four_term() {
        let force_field = force_field_with_positions(2);
        let mut contrib = NonbondedContrib::new(&force_field);
        let pos = [0.0, 0.0, 0.0, 1.5, 0.0, 0.0];

        contrib.add_term(0, 1, None, true, -0.5, MMFF_DIELECTRIC_DISTANCE, true);

        let corr_dist = (1.5_f64 + 0.05).powi(2);
        let expected = 332.0716 * -0.5 / corr_dist * 0.75;

        assert_close(contrib.get_energy(&pos), expected);
    }

    #[test]
    fn mmff_nonbondedcontrib_get_energy_sums_vdw_and_electrostatic_contributions() {
        let force_field = force_field_with_positions(3);
        let mut contrib = NonbondedContrib::new(&force_field);
        let params = MmffVdwRijstarEps {
            r_ij_star_unscaled: 3.2,
            epsilon_unscaled: 0.15,
            r_ij_star: 3.0,
            epsilon: 0.12,
        };
        let pos = [0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 0.0, 2.0, 0.0];

        contrib.add_term(
            0,
            1,
            Some(&params),
            true,
            -0.25,
            MMFF_DIELECTRIC_DISTANCE,
            true,
        );
        contrib.add_term(1, 2, None, true, 0.4, MMFF_DIELECTRIC_CONSTANT, false);

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

        assert_close(contrib.get_energy(&pos), expected);
    }

    #[test]
    fn mmff_nonbondedcontrib_get_grad_accumulates_vdw_only_term() {
        let force_field = force_field_with_positions(2);
        let mut contrib = NonbondedContrib::new(&force_field);
        let params = MmffVdwRijstarEps {
            r_ij_star_unscaled: 3.2,
            epsilon_unscaled: 0.15,
            r_ij_star: 3.0,
            epsilon: 0.12,
        };
        let pos = [0.0, 0.0, 0.0, 1.5, 0.0, 0.0];
        let mut grad = vec![0.0; 6];

        contrib.add_term(
            0,
            1,
            Some(&params),
            false,
            0.0,
            MMFF_DIELECTRIC_DISTANCE,
            false,
        );
        contrib.get_grad(&pos, &mut grad);

        let de_dr = source_vdw_grad_component(3.0, 0.12, 1.5);
        let expected = vec![-1.5 * de_dr, 0.0, 0.0, 1.5 * de_dr, 0.0, 0.0];
        assert_slice_close(&grad, &expected);
    }

    #[test]
    fn mmff_nonbondedcontrib_get_grad_accumulates_charge_only_constant_term() {
        let force_field = force_field_with_positions(2);
        let mut contrib = NonbondedContrib::new(&force_field);
        let pos = [0.0, 0.0, 0.0, 1.5, 0.0, 0.0];
        let mut grad = vec![0.0; 6];

        contrib.add_term(0, 1, None, true, -0.5, MMFF_DIELECTRIC_CONSTANT, false);
        contrib.get_grad(&pos, &mut grad);

        let de_dr = source_ele_grad_component(-0.5, MMFF_DIELECTRIC_CONSTANT, false, 1.5);
        let expected = vec![-1.5 * de_dr, 0.0, 0.0, 1.5 * de_dr, 0.0, 0.0];
        assert_slice_close(&grad, &expected);
    }

    #[test]
    fn mmff_nonbondedcontrib_get_grad_accumulates_combined_term() {
        let force_field = force_field_with_positions(2);
        let mut contrib = NonbondedContrib::new(&force_field);
        let params = MmffVdwRijstarEps {
            r_ij_star_unscaled: 3.2,
            epsilon_unscaled: 0.15,
            r_ij_star: 3.0,
            epsilon: 0.12,
        };
        let pos = [0.0, 0.0, 0.0, 1.5, 0.0, 0.0];
        let mut grad = vec![0.0; 6];

        contrib.add_term(
            0,
            1,
            Some(&params),
            true,
            -0.25,
            MMFF_DIELECTRIC_DISTANCE,
            true,
        );
        contrib.get_grad(&pos, &mut grad);

        let de_dr = source_vdw_grad_component(3.0, 0.12, 1.5)
            + source_ele_grad_component(-0.25, MMFF_DIELECTRIC_DISTANCE, true, 1.5);
        let expected = vec![-1.5 * de_dr, 0.0, 0.0, 1.5 * de_dr, 0.0, 0.0];
        assert_slice_close(&grad, &expected);
    }

    #[test]
    fn mmff_nonbondedcontrib_get_grad_uses_zero_distance_fallbacks() {
        let force_field = force_field_with_positions(2);
        let mut contrib = NonbondedContrib::new(&force_field);
        let params = MmffVdwRijstarEps {
            r_ij_star_unscaled: 3.2,
            epsilon_unscaled: 0.15,
            r_ij_star: 3.0,
            epsilon: 0.12,
        };
        let pos = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let mut grad = vec![0.0; 6];

        contrib.add_term(
            0,
            1,
            Some(&params),
            true,
            -0.25,
            MMFF_DIELECTRIC_DISTANCE,
            true,
        );
        contrib.get_grad(&pos, &mut grad);

        let expected = vec![0.05, 0.05, 0.05, -0.05, -0.05, -0.05];
        assert_slice_close(&grad, &expected);
    }

    #[test]
    fn mmff_nonbondedcontrib_get_grad_zero_distance_returns_before_later_pairs() {
        let force_field = force_field_with_positions(3);
        let mut contrib = NonbondedContrib::new(&force_field);
        let params = MmffVdwRijstarEps {
            r_ij_star_unscaled: 3.2,
            epsilon_unscaled: 0.15,
            r_ij_star: 3.0,
            epsilon: 0.12,
        };
        let pos = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.5, 0.0, 0.0];
        let mut grad = vec![0.0; 9];

        contrib.add_term(
            0,
            1,
            Some(&params),
            true,
            -0.25,
            MMFF_DIELECTRIC_DISTANCE,
            true,
        );
        contrib.add_term(1, 2, None, true, 0.4, MMFF_DIELECTRIC_CONSTANT, false);
        contrib.get_grad(&pos, &mut grad);

        let expected = vec![0.05, 0.05, 0.05, -0.05, -0.05, -0.05, 0.0, 0.0, 0.0];
        assert_slice_close(&grad, &expected);
    }
}
