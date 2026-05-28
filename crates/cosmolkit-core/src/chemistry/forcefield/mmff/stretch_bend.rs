//! Source-backed RDKit MMFF stretch-bend contribution.

use crate::chemistry::forcefield::core::{ForceField, ForceFieldContrib, ForceFieldVec3};

use super::params::{MmffAngle, MmffBond, MmffStbn};

// BEGIN RDKIT CPP CONSTANT ForceFields::MMFF::MDYNE_A_TO_KCAL_MOL (Params.h:37-40)
// RDKit✔️✔️: constexpr double DEG2RAD = M_PI / 180.0;
// RDKit✔️✔️: constexpr double RAD2DEG = 180.0 / M_PI;
// RDKit✔️✔️: constexpr double MDYNE_A_TO_KCAL_MOL = 143.9325;
const DEG2RAD: f64 = std::f64::consts::PI / 180.0;
const RAD2DEG: f64 = 180.0 / std::f64::consts::PI;
const MDYNE_A_TO_KCAL_MOL: f64 = 143.9325;

#[derive(Clone, Debug)]
pub struct StretchBendContrib {
    owner: *const ForceField,
    atom1_indices: Vec<usize>,
    atom2_indices: Vec<usize>,
    atom3_indices: Vec<usize>,
    rest_lengths1: Vec<f64>,
    rest_lengths2: Vec<f64>,
    theta0: Vec<f64>,
    force_constants1: Vec<f64>,
    force_constants2: Vec<f64>,
}

impl StretchBendContrib {
    #[must_use]
    pub fn new(owner: &ForceField) -> Self {
        // BEGIN RDKIT CPP CONSTRUCTOR ForceFields::MMFF::StretchBendContrib::StretchBendContrib (StretchBend.cpp:39-42)
        // RDKit✔️✔️: StretchBendContrib::StretchBendContrib(ForceField *owner) {
        // RDKit✔️✔️:   PRECONDITION(owner, "bad owner");
        // Rust references reproduce RDKit's non-null owner precondition.
        // RDKit✔️✔️:   dp_forceField = owner;
        let owner = owner as *const ForceField;
        // RDKit✔️✔️: }
        Self {
            owner,
            atom1_indices: Vec::new(),
            atom2_indices: Vec::new(),
            atom3_indices: Vec::new(),
            rest_lengths1: Vec::new(),
            rest_lengths2: Vec::new(),
            theta0: Vec::new(),
            force_constants1: Vec::new(),
            force_constants2: Vec::new(),
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
            && self.atom3_indices.is_empty()
            && self.rest_lengths1.is_empty()
            && self.rest_lengths2.is_empty()
            && self.theta0.is_empty()
            && self.force_constants1.is_empty()
            && self.force_constants2.is_empty()
    }

    pub fn add_term(
        &mut self,
        idx1: usize,
        idx2: usize,
        idx3: usize,
        mmff_stbn_params: &MmffStbn,
        mmff_angle_params: &MmffAngle,
        mmff_bond_params1: &MmffBond,
        mmff_bond_params2: &MmffBond,
    ) {
        // BEGIN RDKIT CPP METHOD ForceFields::MMFF::StretchBendContrib::addTerm (StretchBend.cpp:44-66)
        // RDKit✔️✔️: void StretchBendContrib::addTerm(
        // RDKit✔️✔️:     const unsigned int idx1, const unsigned int idx2, const unsigned int idx3,
        // RDKit✔️✔️:     const ForceFields::MMFF::MMFFStbn *mmffStbnParams,
        // RDKit✔️✔️:     const ForceFields::MMFF::MMFFAngle *mmffAngleParams,
        // RDKit✔️✔️:     const ForceFields::MMFF::MMFFBond *mmffBondParams1,
        // RDKit✔️✔️:     const ForceFields::MMFF::MMFFBond *mmffBondParams2) {
        // RDKit✔️✔️:   PRECONDITION(((idx1 != idx2) && (idx2 != idx3) && (idx1 != idx3)),
        // RDKit✔️✔️:                "degenerate points");
        assert!(idx1 != idx2 && idx2 != idx3 && idx1 != idx3);
        let force_field = self.force_field();
        // RDKit✔️✔️:   URANGE_CHECK(idx1, dp_forceField->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx2, dp_forceField->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx3, dp_forceField->positions().size());
        assert!(idx1 < force_field.positions().len());
        assert!(idx2 < force_field.positions().len());
        assert!(idx3 < force_field.positions().len());

        // RDKit✔️✔️:   d_at1Idxs.push_back(idx1);
        // RDKit✔️✔️:   d_at2Idxs.push_back(idx2);
        // RDKit✔️✔️:   d_at3Idxs.push_back(idx3);
        self.atom1_indices.push(idx1);
        self.atom2_indices.push(idx2);
        self.atom3_indices.push(idx3);
        // RDKit✔️✔️:   d_restLen1s.push_back(Utils::calcBondRestLength(mmffBondParams1));
        // RDKit✔️✔️:   d_restLen2s.push_back(Utils::calcBondRestLength(mmffBondParams2));
        self.rest_lengths1
            .push(calc_bond_rest_length(mmff_bond_params1));
        self.rest_lengths2
            .push(calc_bond_rest_length(mmff_bond_params2));
        // RDKit✔️✔️:   d_theta0s.push_back(Utils::calcAngleRestValue(mmffAngleParams));
        self.theta0.push(calc_angle_rest_value(mmff_angle_params));
        // RDKit✔️✔️:   std::pair<double, double> forceConstants =
        // RDKit✔️✔️:       Utils::calcStbnForceConstants(mmffStbnParams);
        let force_constants = calc_stbn_force_constants(mmff_stbn_params);
        // RDKit✔️✔️:   d_forceConstants1.push_back(forceConstants.first);
        // RDKit✔️✔️:   d_forceConstants2.push_back(forceConstants.second);
        self.force_constants1.push(force_constants.0);
        self.force_constants2.push(force_constants.1);
        // RDKit✔️✔️: }
    }

    #[must_use]
    pub fn get_energy(&self, pos: &[f64]) -> f64 {
        // BEGIN RDKIT CPP METHOD ForceFields::MMFF::StretchBendContrib::getEnergy (StretchBend.cpp:68-99)
        // RDKit✔️✔️: double StretchBendContrib::getEnergy(double *pos) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        let force_field = self.force_field();
        // Rust slices reproduce RDKit's non-null pos precondition.

        // RDKit✔️✔️:   double totalEnergy = 0.0;
        // RDKit✔️✔️:   const int numTerms = d_at1Idxs.size();
        let mut total_energy = 0.0;
        let num_terms = self.atom1_indices.len();
        // RDKit✔️✔️:   for (int i = 0; i < numTerms; i++) {
        for i in 0..num_terms {
            // RDKit✔️✔️:     const int16_t at1Idx = d_at1Idxs[i];
            // RDKit✔️✔️:     const int16_t at2Idx = d_at2Idxs[i];
            // RDKit✔️✔️:     const int16_t at3Idx = d_at3Idxs[i];
            let atom1_idx = self.atom1_indices[i];
            let atom2_idx = self.atom2_indices[i];
            let atom3_idx = self.atom3_indices[i];

            // RDKit✔️✔️:     double dist1 = dp_forceField->distance(at1Idx, at2Idx, pos);
            // RDKit✔️✔️:     double dist2 = dp_forceField->distance(at2Idx, at3Idx, pos);
            let dist1 = force_field.distance_const(atom1_idx, atom2_idx, Some(pos));
            let dist2 = force_field.distance_const(atom2_idx, atom3_idx, Some(pos));

            // RDKit✔️✔️:     RDGeom::Point3D p1(pos[3 * at1Idx], pos[3 * at1Idx + 1],
            // RDKit✔️✔️:                        pos[3 * at1Idx + 2]);
            // RDKit✔️✔️:     RDGeom::Point3D p2(pos[3 * at2Idx], pos[3 * at2Idx + 1],
            // RDKit✔️✔️:                        pos[3 * at2Idx + 2]);
            // RDKit✔️✔️:     RDGeom::Point3D p3(pos[3 * at3Idx], pos[3 * at3Idx + 1],
            // RDKit✔️✔️:                        pos[3 * at3Idx + 2]);
            let p1 = point_from_pos(pos, atom1_idx);
            let p2 = point_from_pos(pos, atom2_idx);
            let p3 = point_from_pos(pos, atom3_idx);
            // RDKit✔️✔️:     std::pair<double, double> forceConstantsPair =
            // RDKit✔️✔️:         std::make_pair(d_forceConstants1[i], d_forceConstants2[i]);
            let force_constants_pair = (self.force_constants1[i], self.force_constants2[i]);
            // RDKit✔️✔️:     std::pair<double, double> stretchBendEnergies =
            // RDKit✔️✔️:         Utils::calcStretchBendEnergy(
            // RDKit✔️✔️:             dist1 - d_restLen1s[i], dist2 - d_restLen2s[i],
            // RDKit✔️✔️:             RAD2DEG * acos(Utils::calcCosTheta(p1, p2, p3, dist1, dist2)) -
            // RDKit✔️✔️:                 d_theta0s[i],
            // RDKit✔️✔️:             forceConstantsPair);
            let stretch_bend_energies = calc_stretch_bend_energy(
                dist1 - self.rest_lengths1[i],
                dist2 - self.rest_lengths2[i],
                RAD2DEG * calc_cos_theta(p1, p2, p3, dist1, dist2).acos() - self.theta0[i],
                force_constants_pair,
            );
            // RDKit✔️✔️:     totalEnergy += (stretchBendEnergies.first + stretchBendEnergies.second);
            total_energy += stretch_bend_energies.0 + stretch_bend_energies.1;
            // RDKit✔️✔️:   }
        }
        // RDKit✔️✔️:   return totalEnergy;
        // RDKit✔️✔️: }
        total_energy
    }

    pub fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
        // BEGIN RDKIT CPP METHOD ForceFields::MMFF::StretchBendContrib::getGrad (StretchBend.cpp:101-175)
        // RDKit✔️✔️: void StretchBendContrib::getGrad(double *pos, double *grad) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:   PRECONDITION(grad, "bad vector");
        let force_field = self.force_field();
        // Rust slices reproduce RDKit's non-null pos and grad preconditions.

        // RDKit✔️✔️:   const int numTerms = d_at1Idxs.size();
        let num_terms = self.atom1_indices.len();
        // RDKit✔️✔️:   for (int i = 0; i < numTerms; ++i) {
        for i in 0..num_terms {
            // RDKit✔️✔️:     const int16_t at1Idx = d_at1Idxs[i];
            // RDKit✔️✔️:     const int16_t at2Idx = d_at2Idxs[i];
            // RDKit✔️✔️:     const int16_t at3Idx = d_at3Idxs[i];
            let atom1_idx = self.atom1_indices[i];
            let atom2_idx = self.atom2_indices[i];
            let atom3_idx = self.atom3_indices[i];
            // RDKit✔️✔️:     const double theta0 = d_theta0s[i];
            // RDKit✔️✔️:     const double forceConstant1 = d_forceConstants1[i];
            // RDKit✔️✔️:     const double forceConstant2 = d_forceConstants2[i];
            // RDKit✔️✔️:     const double restLen1 = d_restLen1s[i];
            // RDKit✔️✔️:     const double restLen2 = d_restLen2s[i];
            let theta0 = self.theta0[i];
            let force_constant1 = self.force_constants1[i];
            let force_constant2 = self.force_constants2[i];
            let rest_len1 = self.rest_lengths1[i];
            let rest_len2 = self.rest_lengths2[i];

            // RDKit✔️✔️:     double dist1 = dp_forceField->distance(at1Idx, at2Idx, pos);
            // RDKit✔️✔️:     double dist2 = dp_forceField->distance(at2Idx, at3Idx, pos);
            let dist1 = force_field.distance_const(atom1_idx, atom2_idx, Some(pos));
            let dist2 = force_field.distance_const(atom2_idx, atom3_idx, Some(pos));

            // RDKit✔️✔️:     RDGeom::Point3D p1(pos[3 * at1Idx], pos[3 * at1Idx + 1],
            // RDKit✔️✔️:                        pos[3 * at1Idx + 2]);
            // RDKit✔️✔️:     RDGeom::Point3D p2(pos[3 * at2Idx], pos[3 * at2Idx + 1],
            // RDKit✔️✔️:                        pos[3 * at2Idx + 2]);
            // RDKit✔️✔️:     RDGeom::Point3D p3(pos[3 * at3Idx], pos[3 * at3Idx + 1],
            // RDKit✔️✔️:                        pos[3 * at3Idx + 2]);
            let p1 = point_from_pos(pos, atom1_idx);
            let p2 = point_from_pos(pos, atom2_idx);
            let p3 = point_from_pos(pos, atom3_idx);
            // RDKit✔️✔️:     double *g1 = &(grad[3 * at1Idx]);
            // RDKit✔️✔️:     double *g2 = &(grad[3 * at2Idx]);
            // RDKit✔️✔️:     double *g3 = &(grad[3 * at3Idx]);
            // Rust writes directly into the flat grad slice at RDKit's offsets.

            // RDKit✔️✔️:     RDGeom::Point3D p12 = (p1 - p2) / dist1;
            // RDKit✔️✔️:     RDGeom::Point3D p32 = (p3 - p2) / dist2;
            let p12 = (p1 - p2) / dist1;
            let p32 = (p3 - p2) / dist2;
            // RDKit✔️✔️:     double const c5 = MDYNE_A_TO_KCAL_MOL * DEG2RAD;
            let c5 = MDYNE_A_TO_KCAL_MOL * DEG2RAD;
            // RDKit✔️✔️:     double cosTheta = p12.dotProduct(p32);
            let mut cos_theta = p12.dot_product(p32);
            // RDKit✔️✔️:     clipToOne(cosTheta);
            cos_theta = cos_theta.clamp(-1.0, 1.0);
            // RDKit✔️✔️:     double sinThetaSq = 1.0 - cosTheta * cosTheta;
            let sin_theta_sq = 1.0 - cos_theta * cos_theta;
            // RDKit✔️✔️:     double sinTheta = std::max(sqrt(sinThetaSq), 1.0e-8);
            let sin_theta = sin_theta_sq.sqrt().max(1.0e-8);
            // RDKit✔️✔️:     double angleTerm = RAD2DEG * acos(cosTheta) - theta0;
            let angle_term = RAD2DEG * cos_theta.acos() - theta0;
            // RDKit✔️✔️:     double distTerm = RAD2DEG * (forceConstant1 * (dist1 - restLen1) +
            // RDKit✔️✔️:                                  forceConstant2 * (dist2 - restLen2));
            let dist_term = RAD2DEG
                * (force_constant1 * (dist1 - rest_len1) + force_constant2 * (dist2 - rest_len2));
            // RDKit✔️✔️:     double dCos_dS1 = 1.0 / dist1 * (p32.x - cosTheta * p12.x);
            // RDKit✔️✔️:     double dCos_dS2 = 1.0 / dist1 * (p32.y - cosTheta * p12.y);
            // RDKit✔️✔️:     double dCos_dS3 = 1.0 / dist1 * (p32.z - cosTheta * p12.z);
            let d_cos_ds1 = (p32.x - cos_theta * p12.x) / dist1;
            let d_cos_ds2 = (p32.y - cos_theta * p12.y) / dist1;
            let d_cos_ds3 = (p32.z - cos_theta * p12.z) / dist1;

            // RDKit✔️✔️:     double dCos_dS4 = 1.0 / dist2 * (p12.x - cosTheta * p32.x);
            // RDKit✔️✔️:     double dCos_dS5 = 1.0 / dist2 * (p12.y - cosTheta * p32.y);
            // RDKit✔️✔️:     double dCos_dS6 = 1.0 / dist2 * (p12.z - cosTheta * p32.z);
            let d_cos_ds4 = (p12.x - cos_theta * p32.x) / dist2;
            let d_cos_ds5 = (p12.y - cos_theta * p32.y) / dist2;
            let d_cos_ds6 = (p12.z - cos_theta * p32.z) / dist2;

            // RDKit✔️✔️:     g1[0] += c5 * (p12.x * forceConstant1 * angleTerm +
            // RDKit✔️✔️:                    dCos_dS1 / (-sinTheta) * distTerm);
            // RDKit✔️✔️:     g1[1] += c5 * (p12.y * forceConstant1 * angleTerm +
            // RDKit✔️✔️:                    dCos_dS2 / (-sinTheta) * distTerm);
            // RDKit✔️✔️:     g1[2] += c5 * (p12.z * forceConstant1 * angleTerm +
            // RDKit✔️✔️:                    dCos_dS3 / (-sinTheta) * distTerm);
            grad[3 * atom1_idx] +=
                c5 * (p12.x * force_constant1 * angle_term + d_cos_ds1 / (-sin_theta) * dist_term);
            grad[3 * atom1_idx + 1] +=
                c5 * (p12.y * force_constant1 * angle_term + d_cos_ds2 / (-sin_theta) * dist_term);
            grad[3 * atom1_idx + 2] +=
                c5 * (p12.z * force_constant1 * angle_term + d_cos_ds3 / (-sin_theta) * dist_term);

            // RDKit✔️✔️:     g2[0] +=
            // RDKit✔️✔️:         c5 * ((-p12.x * forceConstant1 - p32.x * forceConstant2) * angleTerm +
            // RDKit✔️✔️:               (-dCos_dS1 - dCos_dS4) / (-sinTheta) * distTerm);
            // RDKit✔️✔️:     g2[1] +=
            // RDKit✔️✔️:         c5 * ((-p12.y * forceConstant1 - p32.y * forceConstant2) * angleTerm +
            // RDKit✔️✔️:               (-dCos_dS2 - dCos_dS5) / (-sinTheta) * distTerm);
            // RDKit✔️✔️:     g2[2] +=
            // RDKit✔️✔️:         c5 * ((-p12.z * forceConstant1 - p32.z * forceConstant2) * angleTerm +
            // RDKit✔️✔️:               (-dCos_dS3 - dCos_dS6) / (-sinTheta) * distTerm);
            grad[3 * atom2_idx] += c5
                * ((-p12.x * force_constant1 - p32.x * force_constant2) * angle_term
                    + (-d_cos_ds1 - d_cos_ds4) / (-sin_theta) * dist_term);
            grad[3 * atom2_idx + 1] += c5
                * ((-p12.y * force_constant1 - p32.y * force_constant2) * angle_term
                    + (-d_cos_ds2 - d_cos_ds5) / (-sin_theta) * dist_term);
            grad[3 * atom2_idx + 2] += c5
                * ((-p12.z * force_constant1 - p32.z * force_constant2) * angle_term
                    + (-d_cos_ds3 - d_cos_ds6) / (-sin_theta) * dist_term);

            // RDKit✔️✔️:     g3[0] += c5 * (p32.x * forceConstant2 * angleTerm +
            // RDKit✔️✔️:                    dCos_dS4 / (-sinTheta) * distTerm);
            // RDKit✔️✔️:     g3[1] += c5 * (p32.y * forceConstant2 * angleTerm +
            // RDKit✔️✔️:                    dCos_dS5 / (-sinTheta) * distTerm);
            // RDKit✔️✔️:     g3[2] += c5 * (p32.z * forceConstant2 * angleTerm +
            // RDKit✔️✔️:                    dCos_dS6 / (-sinTheta) * distTerm);
            grad[3 * atom3_idx] +=
                c5 * (p32.x * force_constant2 * angle_term + d_cos_ds4 / (-sin_theta) * dist_term);
            grad[3 * atom3_idx + 1] +=
                c5 * (p32.y * force_constant2 * angle_term + d_cos_ds5 / (-sin_theta) * dist_term);
            grad[3 * atom3_idx + 2] +=
                c5 * (p32.z * force_constant2 * angle_term + d_cos_ds6 / (-sin_theta) * dist_term);
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
    pub fn atom3_indices(&self) -> &[usize] {
        &self.atom3_indices
    }

    #[must_use]
    pub fn rest_lengths1(&self) -> &[f64] {
        &self.rest_lengths1
    }

    #[must_use]
    pub fn rest_lengths2(&self) -> &[f64] {
        &self.rest_lengths2
    }

    #[must_use]
    pub fn theta0(&self) -> &[f64] {
        &self.theta0
    }

    #[must_use]
    pub fn force_constants1(&self) -> &[f64] {
        &self.force_constants1
    }

    #[must_use]
    pub fn force_constants2(&self) -> &[f64] {
        &self.force_constants2
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
fn calc_bond_rest_length(mmff_bond_params: &MmffBond) -> f64 {
    // BEGIN RDKIT CPP FUNCTION ForceFields::MMFF::Utils::calcBondRestLength (BondStretch.cpp:18-23)
    // RDKit✔️✔️: double calcBondRestLength(const MMFFBond *mmffBondParams) {
    // RDKit✔️✔️:   PRECONDITION(mmffBondParams, "bond parameters not found");
    // Rust references reproduce RDKit's non-null parameter precondition.
    //
    // RDKit✔️✔️:   return mmffBondParams->r0;
    // RDKit✔️✔️: }
    mmff_bond_params.r0
}

#[must_use]
fn calc_angle_rest_value(mmff_angle_params: &MmffAngle) -> f64 {
    // BEGIN RDKIT CPP FUNCTION ForceFields::MMFF::Utils::calcAngleRestValue (AngleBend.cpp:18-23)
    // RDKit✔️✔️: double calcAngleRestValue(const MMFFAngle *mmffAngleParams) {
    // RDKit✔️✔️:   PRECONDITION(mmffAngleParams, "angle parameters not found");
    // Rust references reproduce RDKit's non-null parameter precondition.
    //
    // RDKit✔️✔️:   return mmffAngleParams->theta0;
    // RDKit✔️✔️: }
    mmff_angle_params.theta0
}

#[must_use]
fn calc_stbn_force_constants(mmff_stbn_params: &MmffStbn) -> (f64, f64) {
    // BEGIN RDKIT CPP FUNCTION ForceFields::MMFF::Utils::calcStbnForceConstants (StretchBend.cpp:18-24)
    // RDKit✔️✔️: std::pair<double, double> calcStbnForceConstants(
    // RDKit✔️✔️:     const MMFFStbn *mmffStbnParams) {
    // RDKit✔️✔️:   PRECONDITION(mmffStbnParams, "stretch-bend parameters not found");
    // Rust references reproduce RDKit's non-null parameter precondition.
    //
    // RDKit✔️✔️:   return std::make_pair(mmffStbnParams->kbaIJK, mmffStbnParams->kbaKJI);
    // RDKit✔️✔️: }
    (mmff_stbn_params.kba_ijk, mmff_stbn_params.kba_kji)
}

#[must_use]
fn point_from_pos(pos: &[f64], atom_idx: usize) -> ForceFieldVec3 {
    let offset = 3 * atom_idx;
    ForceFieldVec3::new(pos[offset], pos[offset + 1], pos[offset + 2])
}

#[must_use]
fn calc_cos_theta(
    p1: ForceFieldVec3,
    p2: ForceFieldVec3,
    p3: ForceFieldVec3,
    dist1: f64,
    dist2: f64,
) -> f64 {
    // BEGIN RDKIT CPP FUNCTION ForceFields::MMFF::Utils::calcCosTheta (AngleBend.cpp:25-33)
    // RDKit✔️✔️: double calcCosTheta(RDGeom::Point3D p1, RDGeom::Point3D p2, RDGeom::Point3D p3,
    // RDKit✔️✔️:                     double dist1, double dist2) {
    // RDKit✔️✔️:   RDGeom::Point3D p12 = p1 - p2;
    // RDKit✔️✔️:   RDGeom::Point3D p32 = p3 - p2;
    let p12 = p1 - p2;
    let p32 = p3 - p2;
    // RDKit✔️✔️:   double cosTheta = p12.dotProduct(p32) / (dist1 * dist2);
    let cos_theta = p12.dot_product(p32) / (dist1 * dist2);
    // RDKit✔️✔️:   clipToOne(cosTheta);
    let cos_theta = cos_theta.clamp(-1.0, 1.0);

    // RDKit✔️✔️:   return cosTheta;
    // RDKit✔️✔️: }
    cos_theta
}

#[must_use]
fn calc_stretch_bend_energy(
    delta_dist1: f64,
    delta_dist2: f64,
    delta_theta: f64,
    force_constants: (f64, f64),
) -> (f64, f64) {
    // BEGIN RDKIT CPP FUNCTION ForceFields::MMFF::Utils::calcStretchBendEnergy (StretchBend.cpp:26-36)
    // RDKit✔️✔️: std::pair<double, double> calcStretchBendEnergy(
    // RDKit✔️✔️:     const double deltaDist1, const double deltaDist2, const double deltaTheta,
    // RDKit✔️✔️:     const std::pair<double, double> forceConstants) {
    // RDKit✔️✔️:   double factor = MDYNE_A_TO_KCAL_MOL * DEG2RAD * deltaTheta;
    let factor = MDYNE_A_TO_KCAL_MOL * DEG2RAD * delta_theta;

    // RDKit✔️✔️:   return std::make_pair(factor * forceConstants.first * deltaDist1,
    // RDKit✔️✔️:                         factor * forceConstants.second * deltaDist2);
    // RDKit✔️✔️: }
    (
        factor * force_constants.0 * delta_dist1,
        factor * force_constants.1 * delta_dist2,
    )
}

impl ForceFieldContrib for StretchBendContrib {
    fn copy(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }

    fn set_force_field(&mut self, owner: *const ForceField) {
        self.owner = owner;
    }

    fn get_energy(&self, pos: &[f64]) -> f64 {
        StretchBendContrib::get_energy(self, pos)
    }

    fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
        StretchBendContrib::get_grad(self, pos, grad);
    }
}

#[cfg(test)]
mod tests {
    use super::{StretchBendContrib, calc_cos_theta};
    use crate::chemistry::forcefield::{
        core::{ForceField, ForceFieldVec3},
        mmff::params::{MmffAngle, MmffBond, MmffStbn},
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

    fn source_stretch_bend_energy(
        dist1: f64,
        rest_len1: f64,
        dist2: f64,
        rest_len2: f64,
        theta: f64,
        theta0: f64,
        force_constants: (f64, f64),
    ) -> f64 {
        let delta_theta = theta - theta0;
        let factor = 143.9325 * (std::f64::consts::PI / 180.0) * delta_theta;
        factor * force_constants.0 * (dist1 - rest_len1)
            + factor * force_constants.1 * (dist2 - rest_len2)
    }

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 1.0e-12,
            "actual={actual}, expected={expected}"
        );
    }

    fn assert_slice_close(actual: &[f64], expected: &[f64]) {
        assert_eq!(actual.len(), expected.len());
        for (index, (actual, expected)) in actual.iter().zip(expected.iter()).enumerate() {
            assert!(
                (actual - expected).abs() < 1.0e-12,
                "index={index}, actual={actual}, expected={expected}"
            );
        }
    }

    fn source_stretch_bend_grad(
        pos: &[f64],
        atom_indices: [usize; 3],
        rest_lengths: (f64, f64),
        theta0: f64,
        force_constants: (f64, f64),
    ) -> Vec<f64> {
        let atom1_idx = atom_indices[0];
        let atom2_idx = atom_indices[1];
        let atom3_idx = atom_indices[2];
        let p1 = ForceFieldVec3::new(
            pos[3 * atom1_idx],
            pos[3 * atom1_idx + 1],
            pos[3 * atom1_idx + 2],
        );
        let p2 = ForceFieldVec3::new(
            pos[3 * atom2_idx],
            pos[3 * atom2_idx + 1],
            pos[3 * atom2_idx + 2],
        );
        let p3 = ForceFieldVec3::new(
            pos[3 * atom3_idx],
            pos[3 * atom3_idx + 1],
            pos[3 * atom3_idx + 2],
        );
        let dist1 = (p1 - p2).length();
        let dist2 = (p3 - p2).length();
        let p12 = (p1 - p2) / dist1;
        let p32 = (p3 - p2) / dist2;
        let c5 = 143.9325 * (std::f64::consts::PI / 180.0);
        let cos_theta = p12.dot_product(p32).clamp(-1.0, 1.0);
        let sin_theta_sq = 1.0 - cos_theta * cos_theta;
        let sin_theta = sin_theta_sq.sqrt().max(1.0e-8);
        let angle_term = (180.0 / std::f64::consts::PI) * cos_theta.acos() - theta0;
        let dist_term = (180.0 / std::f64::consts::PI)
            * (force_constants.0 * (dist1 - rest_lengths.0)
                + force_constants.1 * (dist2 - rest_lengths.1));
        let d_cos_ds1 = (p32.x - cos_theta * p12.x) / dist1;
        let d_cos_ds2 = (p32.y - cos_theta * p12.y) / dist1;
        let d_cos_ds3 = (p32.z - cos_theta * p12.z) / dist1;
        let d_cos_ds4 = (p12.x - cos_theta * p32.x) / dist2;
        let d_cos_ds5 = (p12.y - cos_theta * p32.y) / dist2;
        let d_cos_ds6 = (p12.z - cos_theta * p32.z) / dist2;
        let mut grad = vec![0.0; pos.len()];

        grad[3 * atom1_idx] +=
            c5 * (p12.x * force_constants.0 * angle_term + d_cos_ds1 / (-sin_theta) * dist_term);
        grad[3 * atom1_idx + 1] +=
            c5 * (p12.y * force_constants.0 * angle_term + d_cos_ds2 / (-sin_theta) * dist_term);
        grad[3 * atom1_idx + 2] +=
            c5 * (p12.z * force_constants.0 * angle_term + d_cos_ds3 / (-sin_theta) * dist_term);

        grad[3 * atom2_idx] += c5
            * ((-p12.x * force_constants.0 - p32.x * force_constants.1) * angle_term
                + (-d_cos_ds1 - d_cos_ds4) / (-sin_theta) * dist_term);
        grad[3 * atom2_idx + 1] += c5
            * ((-p12.y * force_constants.0 - p32.y * force_constants.1) * angle_term
                + (-d_cos_ds2 - d_cos_ds5) / (-sin_theta) * dist_term);
        grad[3 * atom2_idx + 2] += c5
            * ((-p12.z * force_constants.0 - p32.z * force_constants.1) * angle_term
                + (-d_cos_ds3 - d_cos_ds6) / (-sin_theta) * dist_term);

        grad[3 * atom3_idx] +=
            c5 * (p32.x * force_constants.1 * angle_term + d_cos_ds4 / (-sin_theta) * dist_term);
        grad[3 * atom3_idx + 1] +=
            c5 * (p32.y * force_constants.1 * angle_term + d_cos_ds5 / (-sin_theta) * dist_term);
        grad[3 * atom3_idx + 2] +=
            c5 * (p32.z * force_constants.1 * angle_term + d_cos_ds6 / (-sin_theta) * dist_term);

        grad
    }

    #[test]
    fn mmff_stretchbendcontrib_constructor_stores_owner_pointer() {
        let force_field = ForceField::new(3);

        let contrib = StretchBendContrib::new(&force_field);

        assert_eq!(contrib.owner(), &force_field as *const ForceField);
    }

    #[test]
    fn mmff_stretchbendcontrib_constructor_initializes_no_terms() {
        let force_field = ForceField::new(3);

        let contrib = StretchBendContrib::new(&force_field);

        assert_eq!(contrib.len(), 0);
        assert!(contrib.is_empty());
    }

    #[test]
    fn mmff_stretchbendcontrib_constructor_accepts_empty_force_field_like_rdkit() {
        let force_field = ForceField::new(3);

        let contrib = StretchBendContrib::new(&force_field);

        assert_eq!(force_field.positions().len(), 0);
        assert!(contrib.is_empty());
    }

    #[test]
    fn mmff_stretchbendcontrib_add_term_pushes_source_fields() {
        let force_field = force_field_with_positions(3);
        let mut contrib = StretchBendContrib::new(&force_field);
        let stbn = MmffStbn {
            kba_ijk: 0.227,
            kba_kji: 0.070,
        };
        let angle = MmffAngle {
            ka: 0.636,
            theta0: 110.549,
        };
        let bond1 = MmffBond {
            kb: 4.258,
            r0: 1.508,
        };
        let bond2 = MmffBond {
            kb: 4.766,
            r0: 1.093,
        };

        contrib.add_term(0, 1, 2, &stbn, &angle, &bond1, &bond2);

        assert_eq!(contrib.atom1_indices(), &[0]);
        assert_eq!(contrib.atom2_indices(), &[1]);
        assert_eq!(contrib.atom3_indices(), &[2]);
        assert_eq!(contrib.rest_lengths1(), &[1.508]);
        assert_eq!(contrib.rest_lengths2(), &[1.093]);
        assert_eq!(contrib.theta0(), &[110.549]);
        assert_eq!(contrib.force_constants1(), &[0.227]);
        assert_eq!(contrib.force_constants2(), &[0.070]);
    }

    #[test]
    fn mmff_stretchbendcontrib_add_term_appends_multiple_terms() {
        let force_field = force_field_with_positions(4);
        let mut contrib = StretchBendContrib::new(&force_field);
        let first_stbn = MmffStbn {
            kba_ijk: 0.227,
            kba_kji: 0.070,
        };
        let second_stbn = MmffStbn {
            kba_ijk: 0.300,
            kba_kji: 0.400,
        };
        let first_angle = MmffAngle {
            ka: 0.636,
            theta0: 110.549,
        };
        let second_angle = MmffAngle {
            ka: 1.006,
            theta0: 110.265,
        };
        let first_bond1 = MmffBond {
            kb: 4.258,
            r0: 1.508,
        };
        let first_bond2 = MmffBond {
            kb: 4.766,
            r0: 1.093,
        };
        let second_bond1 = MmffBond {
            kb: 3.529,
            r0: 1.330,
        };
        let second_bond2 = MmffBond {
            kb: 3.366,
            r0: 1.520,
        };

        contrib.add_term(
            0,
            1,
            2,
            &first_stbn,
            &first_angle,
            &first_bond1,
            &first_bond2,
        );
        contrib.add_term(
            1,
            2,
            3,
            &second_stbn,
            &second_angle,
            &second_bond1,
            &second_bond2,
        );

        assert_eq!(contrib.atom1_indices(), &[0, 1]);
        assert_eq!(contrib.atom2_indices(), &[1, 2]);
        assert_eq!(contrib.atom3_indices(), &[2, 3]);
        assert_eq!(contrib.rest_lengths1(), &[1.508, 1.330]);
        assert_eq!(contrib.rest_lengths2(), &[1.093, 1.520]);
        assert_eq!(contrib.theta0(), &[110.549, 110.265]);
        assert_eq!(contrib.force_constants1(), &[0.227, 0.300]);
        assert_eq!(contrib.force_constants2(), &[0.070, 0.400]);
    }

    #[test]
    #[should_panic]
    fn mmff_stretchbendcontrib_add_term_rejects_first_second_degenerate() {
        let force_field = force_field_with_positions(3);
        let mut contrib = StretchBendContrib::new(&force_field);
        let stbn = MmffStbn {
            kba_ijk: 0.227,
            kba_kji: 0.070,
        };
        let angle = MmffAngle {
            ka: 0.636,
            theta0: 110.549,
        };
        let bond = MmffBond {
            kb: 4.258,
            r0: 1.508,
        };

        contrib.add_term(0, 0, 2, &stbn, &angle, &bond, &bond);
    }

    #[test]
    #[should_panic]
    fn mmff_stretchbendcontrib_add_term_rejects_second_third_degenerate() {
        let force_field = force_field_with_positions(3);
        let mut contrib = StretchBendContrib::new(&force_field);
        let stbn = MmffStbn {
            kba_ijk: 0.227,
            kba_kji: 0.070,
        };
        let angle = MmffAngle {
            ka: 0.636,
            theta0: 110.549,
        };
        let bond = MmffBond {
            kb: 4.258,
            r0: 1.508,
        };

        contrib.add_term(0, 1, 1, &stbn, &angle, &bond, &bond);
    }

    #[test]
    #[should_panic]
    fn mmff_stretchbendcontrib_add_term_rejects_first_third_degenerate() {
        let force_field = force_field_with_positions(3);
        let mut contrib = StretchBendContrib::new(&force_field);
        let stbn = MmffStbn {
            kba_ijk: 0.227,
            kba_kji: 0.070,
        };
        let angle = MmffAngle {
            ka: 0.636,
            theta0: 110.549,
        };
        let bond = MmffBond {
            kb: 4.258,
            r0: 1.508,
        };

        contrib.add_term(0, 1, 0, &stbn, &angle, &bond, &bond);
    }

    #[test]
    #[should_panic]
    fn mmff_stretchbendcontrib_add_term_rejects_first_index_out_of_range() {
        let force_field = force_field_with_positions(3);
        let mut contrib = StretchBendContrib::new(&force_field);
        let stbn = MmffStbn {
            kba_ijk: 0.227,
            kba_kji: 0.070,
        };
        let angle = MmffAngle {
            ka: 0.636,
            theta0: 110.549,
        };
        let bond = MmffBond {
            kb: 4.258,
            r0: 1.508,
        };

        contrib.add_term(3, 1, 2, &stbn, &angle, &bond, &bond);
    }

    #[test]
    #[should_panic]
    fn mmff_stretchbendcontrib_add_term_rejects_second_index_out_of_range() {
        let force_field = force_field_with_positions(3);
        let mut contrib = StretchBendContrib::new(&force_field);
        let stbn = MmffStbn {
            kba_ijk: 0.227,
            kba_kji: 0.070,
        };
        let angle = MmffAngle {
            ka: 0.636,
            theta0: 110.549,
        };
        let bond = MmffBond {
            kb: 4.258,
            r0: 1.508,
        };

        contrib.add_term(0, 3, 2, &stbn, &angle, &bond, &bond);
    }

    #[test]
    #[should_panic]
    fn mmff_stretchbendcontrib_add_term_rejects_third_index_out_of_range() {
        let force_field = force_field_with_positions(3);
        let mut contrib = StretchBendContrib::new(&force_field);
        let stbn = MmffStbn {
            kba_ijk: 0.227,
            kba_kji: 0.070,
        };
        let angle = MmffAngle {
            ka: 0.636,
            theta0: 110.549,
        };
        let bond = MmffBond {
            kb: 4.258,
            r0: 1.508,
        };

        contrib.add_term(0, 1, 3, &stbn, &angle, &bond, &bond);
    }

    #[test]
    fn mmff_stretchbendcontrib_get_energy_returns_zero_without_terms() {
        let force_field = force_field_with_positions(3);
        let contrib = StretchBendContrib::new(&force_field);
        let pos = [0.0, 0.0, 0.0, 1.4, 0.0, 0.0, 1.4, 1.2, 0.0];

        assert_close(contrib.get_energy(&pos), 0.0);
    }

    #[test]
    fn mmff_stretchbendcontrib_get_energy_uses_source_formula() {
        let force_field = initialized_force_field_with_positions(3);
        let mut contrib = StretchBendContrib::new(&force_field);
        let stbn = MmffStbn {
            kba_ijk: 0.227,
            kba_kji: 0.070,
        };
        let angle = MmffAngle {
            ka: 0.636,
            theta0: 110.549,
        };
        let bond1 = MmffBond {
            kb: 4.258,
            r0: 1.508,
        };
        let bond2 = MmffBond {
            kb: 4.766,
            r0: 1.093,
        };
        let pos = [0.0, 0.0, 0.0, 1.6, 0.0, 0.0, 1.6, 1.4, 0.0];

        contrib.add_term(0, 1, 2, &stbn, &angle, &bond1, &bond2);

        let expected =
            source_stretch_bend_energy(1.6, 1.508, 1.4, 1.093, 90.0, 110.549, (0.227, 0.070));
        assert_close(contrib.get_energy(&pos), expected);
    }

    #[test]
    fn mmff_stretchbendcontrib_get_energy_sums_multiple_source_terms() {
        let force_field = initialized_force_field_with_positions(4);
        let mut contrib = StretchBendContrib::new(&force_field);
        let first_stbn = MmffStbn {
            kba_ijk: 0.227,
            kba_kji: 0.070,
        };
        let second_stbn = MmffStbn {
            kba_ijk: 0.300,
            kba_kji: 0.400,
        };
        let first_angle = MmffAngle {
            ka: 0.636,
            theta0: 110.549,
        };
        let second_angle = MmffAngle {
            ka: 1.006,
            theta0: 110.265,
        };
        let first_bond1 = MmffBond {
            kb: 4.258,
            r0: 1.508,
        };
        let first_bond2 = MmffBond {
            kb: 4.766,
            r0: 1.093,
        };
        let second_bond1 = MmffBond {
            kb: 3.529,
            r0: 1.330,
        };
        let second_bond2 = MmffBond {
            kb: 3.366,
            r0: 1.520,
        };
        let pos = [
            0.0, 0.0, 0.0, //
            1.6, 0.0, 0.0, //
            1.6, 1.4, 0.0, //
            2.9, 1.4, 0.0,
        ];

        contrib.add_term(
            0,
            1,
            2,
            &first_stbn,
            &first_angle,
            &first_bond1,
            &first_bond2,
        );
        contrib.add_term(
            1,
            2,
            3,
            &second_stbn,
            &second_angle,
            &second_bond1,
            &second_bond2,
        );

        let first =
            source_stretch_bend_energy(1.6, 1.508, 1.4, 1.093, 90.0, 110.549, (0.227, 0.070));
        let second =
            source_stretch_bend_energy(1.4, 1.330, 1.3, 1.520, 90.0, 110.265, (0.300, 0.400));
        assert_close(contrib.get_energy(&pos), first + second);
    }

    #[test]
    fn mmff_stretchbendcontrib_get_energy_cos_theta_clips_upper_source_bound() {
        let p1 = ForceFieldVec3::new(1.0 + f64::EPSILON, 0.0, 0.0);
        let p2 = ForceFieldVec3::new(0.0, 0.0, 0.0);
        let p3 = ForceFieldVec3::new(1.0, 0.0, 0.0);

        assert_close(calc_cos_theta(p1, p2, p3, 1.0, 1.0), 1.0);
    }

    #[test]
    fn mmff_stretchbendcontrib_get_energy_cos_theta_clips_lower_source_bound() {
        let p1 = ForceFieldVec3::new(-(1.0 + f64::EPSILON), 0.0, 0.0);
        let p2 = ForceFieldVec3::new(0.0, 0.0, 0.0);
        let p3 = ForceFieldVec3::new(1.0, 0.0, 0.0);

        assert_close(calc_cos_theta(p1, p2, p3, 1.0, 1.0), -1.0);
    }

    #[test]
    fn mmff_stretchbendcontrib_get_grad_leaves_gradient_without_terms() {
        let force_field = ForceField::new(3);
        let contrib = StretchBendContrib::new(&force_field);
        let mut grad = vec![1.0, 2.0, 3.0];

        contrib.get_grad(&[], &mut grad);

        assert_eq!(grad, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn mmff_stretchbendcontrib_get_grad_uses_source_formula() {
        let force_field = initialized_force_field_with_positions(3);
        let mut contrib = StretchBendContrib::new(&force_field);
        let stbn = MmffStbn {
            kba_ijk: 0.227,
            kba_kji: 0.070,
        };
        let angle = MmffAngle {
            ka: 0.636,
            theta0: 110.549,
        };
        let bond1 = MmffBond {
            kb: 4.258,
            r0: 1.508,
        };
        let bond2 = MmffBond {
            kb: 4.766,
            r0: 1.093,
        };
        let pos = [0.0, 0.0, 0.0, 1.6, 0.0, 0.0, 1.6, 1.4, 0.0];
        let mut grad = vec![0.0; pos.len()];

        contrib.add_term(0, 1, 2, &stbn, &angle, &bond1, &bond2);
        contrib.get_grad(&pos, &mut grad);

        let expected =
            source_stretch_bend_grad(&pos, [0, 1, 2], (1.508, 1.093), 110.549, (0.227, 0.070));
        assert_slice_close(&grad, &expected);
    }

    #[test]
    fn mmff_stretchbendcontrib_get_grad_accumulates_multiple_terms() {
        let force_field = initialized_force_field_with_positions(4);
        let mut contrib = StretchBendContrib::new(&force_field);
        let first_stbn = MmffStbn {
            kba_ijk: 0.227,
            kba_kji: 0.070,
        };
        let second_stbn = MmffStbn {
            kba_ijk: 0.300,
            kba_kji: 0.400,
        };
        let first_angle = MmffAngle {
            ka: 0.636,
            theta0: 110.549,
        };
        let second_angle = MmffAngle {
            ka: 1.006,
            theta0: 110.265,
        };
        let first_bond1 = MmffBond {
            kb: 4.258,
            r0: 1.508,
        };
        let first_bond2 = MmffBond {
            kb: 4.766,
            r0: 1.093,
        };
        let second_bond1 = MmffBond {
            kb: 3.529,
            r0: 1.330,
        };
        let second_bond2 = MmffBond {
            kb: 3.366,
            r0: 1.520,
        };
        let pos = [
            0.0, 0.0, 0.0, //
            1.6, 0.0, 0.0, //
            1.6, 1.4, 0.0, //
            2.9, 1.4, 0.0,
        ];
        let mut grad = vec![0.0; pos.len()];

        contrib.add_term(
            0,
            1,
            2,
            &first_stbn,
            &first_angle,
            &first_bond1,
            &first_bond2,
        );
        contrib.add_term(
            1,
            2,
            3,
            &second_stbn,
            &second_angle,
            &second_bond1,
            &second_bond2,
        );
        contrib.get_grad(&pos, &mut grad);

        let mut expected =
            source_stretch_bend_grad(&pos, [0, 1, 2], (1.508, 1.093), 110.549, (0.227, 0.070));
        let second =
            source_stretch_bend_grad(&pos, [1, 2, 3], (1.330, 1.520), 110.265, (0.300, 0.400));
        for (expected, second) in expected.iter_mut().zip(second.iter()) {
            *expected += second;
        }
        assert_slice_close(&grad, &expected);
    }

    #[test]
    #[should_panic(expected = "not initialized")]
    fn mmff_stretchbendcontrib_get_grad_requires_initialized_force_field_for_terms() {
        let force_field = force_field_with_positions(3);
        let mut contrib = StretchBendContrib::new(&force_field);
        let stbn = MmffStbn {
            kba_ijk: 0.227,
            kba_kji: 0.070,
        };
        let angle = MmffAngle {
            ka: 0.636,
            theta0: 110.549,
        };
        let bond1 = MmffBond {
            kb: 4.258,
            r0: 1.508,
        };
        let bond2 = MmffBond {
            kb: 4.766,
            r0: 1.093,
        };
        let mut grad = vec![0.0; 9];

        contrib.add_term(0, 1, 2, &stbn, &angle, &bond1, &bond2);
        contrib.get_grad(&[0.0, 0.0, 0.0, 1.6, 0.0, 0.0, 1.6, 1.4, 0.0], &mut grad);
    }

    #[test]
    fn mmff_stretchbendcontrib_get_grad_uses_singularity_sine_clamp() {
        let force_field = initialized_force_field_with_positions(3);
        let mut contrib = StretchBendContrib::new(&force_field);
        let stbn = MmffStbn {
            kba_ijk: 0.227,
            kba_kji: 0.070,
        };
        let angle = MmffAngle {
            ka: 0.636,
            theta0: 110.549,
        };
        let bond1 = MmffBond {
            kb: 4.258,
            r0: 1.508,
        };
        let bond2 = MmffBond {
            kb: 4.766,
            r0: 1.093,
        };
        let pos = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0];
        let mut grad = vec![0.0; pos.len()];

        contrib.add_term(0, 1, 2, &stbn, &angle, &bond1, &bond2);
        contrib.get_grad(&pos, &mut grad);

        assert!(grad.iter().all(|value| value.is_finite()));
    }
}
