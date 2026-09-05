//! Source-backed RDKit MMFF out-of-plane bend contribution.

use crate::chemistry::forcefield::core::{ForceField, ForceFieldContrib, ForceFieldVec3};

use super::params::MmffOop;

// BEGIN RDKIT CPP CONSTANT ForceFields::MMFF oop-bend constants (Params.h:37-39)
// RDKit✔️✔️: constexpr double DEG2RAD = M_PI / 180.0;
const DEG2RAD: f64 = std::f64::consts::PI / 180.0;
// RDKit✔️✔️: constexpr double RAD2DEG = 180.0 / M_PI;
const RAD2DEG: f64 = 180.0 / std::f64::consts::PI;
// RDKit✔️✔️: constexpr double MDYNE_A_TO_KCAL_MOL = 143.9325;
const MDYNE_A_TO_KCAL_MOL: f64 = 143.9325;

fn is_double_zero(value: f64) -> bool {
    // BEGIN RDKIT CPP HELPER ForceFields::isDoubleZero (Params.h via OopBend.cpp)
    // RDKit✔️✔️: inline bool isDoubleZero(const double x) {
    // RDKit✔️✔️:   return ((x < 1.0e-10) && (x > -1.0e-10));
    // RDKit✔️✔️: }
    value < 1.0e-10 && value > -1.0e-10
}

fn clip_to_one(value: &mut f64) {
    // BEGIN RDKIT CPP HELPER ForceFields::clipToOne (Params.h via OopBend.cpp)
    // RDKit✔️✔️: inline void clipToOne(double &x) { x = std::clamp(x, -1.0, 1.0); }
    *value = value.clamp(-1.0, 1.0);
}

fn point_from_pos(pos: &[f64], atom_idx: usize) -> ForceFieldVec3 {
    ForceFieldVec3::new(pos[3 * atom_idx], pos[3 * atom_idx + 1], pos[3 * atom_idx + 2])
}

fn calc_oop_chi(
    i_point: ForceFieldVec3,
    j_point: ForceFieldVec3,
    k_point: ForceFieldVec3,
    l_point: ForceFieldVec3,
) -> f64 {
    // BEGIN RDKIT CPP FUNCTION ForceFields::MMFF::Utils::calcOopChi (OopBend.cpp:16-31)
    // RDKit✔️✔️: double calcOopChi(const RDGeom::Point3D &iPoint, const RDGeom::Point3D &jPoint,
    // RDKit✔️✔️:                   const RDGeom::Point3D &kPoint,
    // RDKit✔️✔️:                   const RDGeom::Point3D &lPoint) {
    // RDKit✔️✔️:   RDGeom::Point3D rJI = iPoint - jPoint;
    let mut r_ji = i_point - j_point;
    // RDKit✔️✔️:   RDGeom::Point3D rJK = kPoint - jPoint;
    let mut r_jk = k_point - j_point;
    // RDKit✔️✔️:   RDGeom::Point3D rJL = lPoint - jPoint;
    let mut r_jl = l_point - j_point;
    // RDKit✔️✔️:   rJI /= rJI.length();
    r_ji /= r_ji.length();
    // RDKit✔️✔️:   rJK /= rJK.length();
    r_jk /= r_jk.length();
    // RDKit✔️✔️:   rJL /= rJL.length();
    r_jl /= r_jl.length();

    // RDKit✔️✔️:   RDGeom::Point3D n = rJI.crossProduct(rJK);
    let mut n = r_ji.cross_product(r_jk);
    // RDKit✔️✔️:   n /= n.length();
    n /= n.length();
    // RDKit✔️✔️:   double sinChi = n.dotProduct(rJL);
    let mut sin_chi = n.dot_product(r_jl);
    // RDKit✔️✔️:   clipToOne(sinChi);
    clip_to_one(&mut sin_chi);

    // RDKit✔️✔️:   return RAD2DEG * asin(sinChi);
    // RDKit✔️✔️: }
    RAD2DEG * sin_chi.asin()
}

fn calc_oop_bend_energy(chi: f64, koop: f64) -> f64 {
    // BEGIN RDKIT CPP FUNCTION ForceFields::MMFF::Utils::calcOopBendEnergy (OopBend.cpp:38-41)
    // RDKit✔️✔️: double calcOopBendEnergy(const double chi, const double koop) {
    // RDKit✔️✔️:   double const c2 = MDYNE_A_TO_KCAL_MOL * DEG2RAD * DEG2RAD;
    let c2 = MDYNE_A_TO_KCAL_MOL * DEG2RAD * DEG2RAD;
    // RDKit✔️✔️:   return (0.5 * c2 * koop * chi * chi);
    // RDKit✔️✔️: }
    0.5 * c2 * koop * chi * chi
}

#[derive(Clone, Debug)]
pub struct OopBendContrib {
    owner: *const ForceField,
    atom1_indices: Vec<usize>,
    atom2_indices: Vec<usize>,
    atom3_indices: Vec<usize>,
    atom4_indices: Vec<usize>,
    koop: Vec<f64>,
}

impl OopBendContrib {
    #[must_use]
    pub fn new(owner: &ForceField) -> Self {
        // BEGIN RDKIT CPP CONSTRUCTOR ForceFields::MMFF::OopBendContrib::OopBendContrib (OopBend.cpp:41-44)
        // RDKit✔️✔️: OopBendContrib::OopBendContrib(ForceField *owner) {
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
            atom4_indices: Vec::new(),
            koop: Vec::new(),
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
            && self.atom4_indices.is_empty()
            && self.koop.is_empty()
    }

    pub fn add_term(&mut self, idx1: usize, idx2: usize, idx3: usize, idx4: usize, mmff_oop_params: &MmffOop) {
        // BEGIN RDKIT CPP METHOD ForceFields::MMFF::OopBendContrib::addTerm (OopBend.cpp:46-64)
        // RDKit✔️✔️: void OopBendContrib::addTerm(unsigned int idx1,
        // RDKit✔️✔️:                              unsigned int idx2,
        // RDKit✔️✔️:                              unsigned int idx3,
        // RDKit✔️✔️:                              unsigned int idx4,
        // RDKit✔️✔️:                              const ForceFields::MMFF::MMFFOop *mmffOopParams) {
        // RDKit✔️✔️:   PRECONDITION(mmffOopParams, "no OOP parameters");
        // Rust references use `&MmffOop`, so the source non-null precondition is enforced by type.
        // RDKit✔️✔️:   PRECONDITION((idx1 != idx2) && (idx1 != idx3) && (idx1 != idx4) &&
        // RDKit✔️✔️:                    (idx2 != idx3) && (idx2 != idx4) && (idx3 != idx4),
        // RDKit✔️✔️:                "degenerate points");
        assert!(
            idx1 != idx2 && idx1 != idx3 && idx1 != idx4 && idx2 != idx3 && idx2 != idx4 && idx3 != idx4,
            "degenerate points"
        );
        let force_field = self.force_field();
        // RDKit✔️✔️:   URANGE_CHECK(idx1, dp_forceField->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx2, dp_forceField->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx3, dp_forceField->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx4, dp_forceField->positions().size());
        assert!(idx1 < force_field.positions().len());
        assert!(idx2 < force_field.positions().len());
        assert!(idx3 < force_field.positions().len());
        assert!(idx4 < force_field.positions().len());

        // RDKit✔️✔️:   d_at1Idxs.push_back(idx1);
        // RDKit✔️✔️:   d_at2Idxs.push_back(idx2);
        // RDKit✔️✔️:   d_at3Idxs.push_back(idx3);
        // RDKit✔️✔️:   d_at4Idxs.push_back(idx4);
        self.atom1_indices.push(idx1);
        self.atom2_indices.push(idx2);
        self.atom3_indices.push(idx3);
        self.atom4_indices.push(idx4);
        // RDKit✔️✔️:   d_koop.push_back(mmffOopParams->koop);
        self.koop.push(mmff_oop_params.koop);
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
    pub fn atom4_indices(&self) -> &[usize] {
        &self.atom4_indices
    }

    #[must_use]
    pub fn koop(&self) -> &[f64] {
        &self.koop
    }

    #[must_use]
    pub fn get_energy(&self, pos: &[f64]) -> f64 {
        // BEGIN RDKIT CPP METHOD ForceFields::MMFF::OopBendContrib::getEnergy (OopBend.cpp:66-88)
        // RDKit✔️✔️: double OopBendContrib::getEnergy(double *pos) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        let _force_field = self.force_field();
        // Rust slices reproduce RDKit's non-null pos precondition.

        // RDKit✔️✔️:   const int numTerms = d_at1Idxs.size();
        // RDKit✔️✔️:   double totalEnergy = 0.0;
        let num_terms = self.atom1_indices.len();
        let mut total_energy = 0.0;
        // RDKit✔️✔️:   for (int i = 0; i < numTerms; ++i) {
        for i in 0..num_terms {
            // RDKit✔️✔️:     const int d_at1Idx = d_at1Idxs[i];
            // RDKit✔️✔️:     const int d_at2Idx = d_at2Idxs[i];
            // RDKit✔️✔️:     const int d_at3Idx = d_at3Idxs[i];
            // RDKit✔️✔️:     const int d_at4Idx = d_at4Idxs[i];
            let atom1_idx = self.atom1_indices[i];
            let atom2_idx = self.atom2_indices[i];
            let atom3_idx = self.atom3_indices[i];
            let atom4_idx = self.atom4_indices[i];

            // RDKit✔️✔️:     RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
            // RDKit✔️✔️:                        pos[3 * d_at1Idx + 2]);
            let p1 = point_from_pos(pos, atom1_idx);
            // RDKit✔️✔️:     RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
            // RDKit✔️✔️:                        pos[3 * d_at2Idx + 2]);
            let p2 = point_from_pos(pos, atom2_idx);
            // RDKit✔️✔️:     RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
            // RDKit✔️✔️:                        pos[3 * d_at3Idx + 2]);
            let p3 = point_from_pos(pos, atom3_idx);
            // RDKit✔️✔️:     RDGeom::Point3D p4(pos[3 * d_at4Idx], pos[3 * d_at4Idx + 1],
            // RDKit✔️✔️:                        pos[3 * d_at4Idx + 2]);
            let p4 = point_from_pos(pos, atom4_idx);

            // RDKit✔️✔️:     totalEnergy += Utils::calcOopBendEnergy(Utils::calcOopChi(p1, p2, p3, p4), d_koop[i]);
            total_energy += calc_oop_bend_energy(calc_oop_chi(p1, p2, p3, p4), self.koop[i]);
            // RDKit✔️✔️:   }
        }
        // RDKit✔️✔️:   return totalEnergy;
        // RDKit✔️✔️: }
        total_energy
    }

    pub fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
        // BEGIN RDKIT CPP METHOD ForceFields::MMFF::OopBendContrib::getGrad (OopBend.cpp:90-100)
        // RDKit✔️✔️: void OopBendContrib::getGrad(double* pos, double* grad) const {
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:   PRECONDITION(grad, "bad vector");
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        let _force_field = self.force_field();
        // Rust slices reproduce RDKit's non-null pos and grad preconditions.

        // RDKit✔️✔️:   const int numTerms = d_at1Idxs.size();
        let num_terms = self.atom1_indices.len();
        // RDKit✔️✔️:   for (int i =0; i < numTerms; i++) {
        for i in 0..num_terms {
            // RDKit✔️✔️:     getSingleGrad(pos, grad, i);
            self.get_single_grad(pos, grad, i);
            // RDKit✔️✔️:   }
        }
        // RDKit✔️✔️: }
    }

    fn force_field(&self) -> &ForceField {
        // RDKit✔️✔️: dp_forceField
        unsafe { &*self.owner }
    }

    fn get_single_grad(&self, pos: &[f64], grad: &mut [f64], term_idx: usize) {
        // BEGIN RDKIT CPP METHOD ForceFields::MMFF::OopBendContrib::getSingleGrad (OopBend.cpp:102-176)
        // RDKit✔️✔️: void OopBendContrib::getSingleGrad(double *pos, double *grad, unsigned int termIdx) const {

        // RDKit✔️✔️:   const int d_at1Idx = d_at1Idxs[termIdx];
        // RDKit✔️✔️:   const int d_at2Idx = d_at2Idxs[termIdx];
        // RDKit✔️✔️:   const int d_at3Idx = d_at3Idxs[termIdx];
        // RDKit✔️✔️:   const int d_at4Idx = d_at4Idxs[termIdx];
        let atom1_idx = self.atom1_indices[term_idx];
        let atom2_idx = self.atom2_indices[term_idx];
        let atom3_idx = self.atom3_indices[term_idx];
        let atom4_idx = self.atom4_indices[term_idx];

        // RDKit✔️✔️:   RDGeom::Point3D iPoint(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
        // RDKit✔️✔️:                          pos[3 * d_at1Idx + 2]);
        let i_point = point_from_pos(pos, atom1_idx);
        // RDKit✔️✔️:   RDGeom::Point3D jPoint(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
        // RDKit✔️✔️:                          pos[3 * d_at2Idx + 2]);
        let j_point = point_from_pos(pos, atom2_idx);
        // RDKit✔️✔️:   RDGeom::Point3D kPoint(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
        // RDKit✔️✔️:                          pos[3 * d_at3Idx + 2]);
        let k_point = point_from_pos(pos, atom3_idx);
        // RDKit✔️✔️:   RDGeom::Point3D lPoint(pos[3 * d_at4Idx], pos[3 * d_at4Idx + 1],
        // RDKit✔️✔️:                          pos[3 * d_at4Idx + 2]);
        let l_point = point_from_pos(pos, atom4_idx);
        let atom1_offset = 3 * atom1_idx;
        let atom2_offset = 3 * atom2_idx;
        let atom3_offset = 3 * atom3_idx;
        let atom4_offset = 3 * atom4_idx;

        // RDKit✔️✔️:   RDGeom::Point3D rJI = iPoint - jPoint;
        let mut r_ji = i_point - j_point;
        // RDKit✔️✔️:   RDGeom::Point3D rJK = kPoint - jPoint;
        let mut r_jk = k_point - j_point;
        // RDKit✔️✔️:   RDGeom::Point3D rJL = lPoint - jPoint;
        let mut r_jl = l_point - j_point;
        // RDKit✔️✔️:   double dJI = rJI.length();
        let d_ji = r_ji.length();
        // RDKit✔️✔️:   double dJK = rJK.length();
        let d_jk = r_jk.length();
        // RDKit✔️✔️:   double dJL = rJL.length();
        let d_jl = r_jl.length();
        // RDKit✔️✔️:   if (isDoubleZero(dJI) || isDoubleZero(dJK) || isDoubleZero(dJL)) {
        // RDKit✔️✔️:     return;
        // RDKit✔️✔️:   }
        if is_double_zero(d_ji) || is_double_zero(d_jk) || is_double_zero(d_jl) {
            return;
        }
        // RDKit✔️✔️:   rJI /= dJI;
        r_ji /= d_ji;
        // RDKit✔️✔️:   rJK /= dJK;
        r_jk /= d_jk;
        // RDKit✔️✔️:   rJL /= dJL;
        r_jl /= d_jl;

        // RDKit✔️✔️:   RDGeom::Point3D n = (-rJI).crossProduct(rJK);
        let mut n = (-r_ji).cross_product(r_jk);
        // RDKit✔️✔️:   n /= n.length();
        n /= n.length();
        // RDKit✔️✔️:   double const c2 = MDYNE_A_TO_KCAL_MOL * DEG2RAD * DEG2RAD;
        let c2 = MDYNE_A_TO_KCAL_MOL * DEG2RAD * DEG2RAD;
        // RDKit✔️✔️:   double sinChi = rJL.dotProduct(n);
        let mut sin_chi = r_jl.dot_product(n);
        // RDKit✔️✔️:   clipToOne(sinChi);
        clip_to_one(&mut sin_chi);
        // RDKit✔️✔️:   double cosChiSq = 1.0 - sinChi * sinChi;
        let cos_chi_sq = 1.0 - sin_chi * sin_chi;
        // RDKit✔️✔️:   double cosChi = std::max(((cosChiSq > 0.0) ? sqrt(cosChiSq) : 0.0), 1.0e-8);
        let cos_chi = (if cos_chi_sq > 0.0 { cos_chi_sq.sqrt() } else { 0.0 }).max(1.0e-8);
        // RDKit✔️✔️:   double chi = RAD2DEG * asin(sinChi);
        let chi = RAD2DEG * sin_chi.asin();
        // RDKit✔️✔️:   double cosTheta = rJI.dotProduct(rJK);
        let mut cos_theta = r_ji.dot_product(r_jk);
        // RDKit✔️✔️:   clipToOne(cosTheta);
        clip_to_one(&mut cos_theta);
        // RDKit✔️✔️:   double sinThetaSq = std::max(1.0 - cosTheta * cosTheta, 1.0e-8);
        let sin_theta_sq = (1.0 - cos_theta * cos_theta).max(1.0e-8);
        // RDKit✔️✔️:   double sinTheta =
        // RDKit✔️✔️:       std::max(((sinThetaSq > 0.0) ? sqrt(sinThetaSq) : 0.0), 1.0e-8);
        let sin_theta = (if sin_theta_sq > 0.0 { sin_theta_sq.sqrt() } else { 0.0 }).max(1.0e-8);

        // RDKit✔️✔️:   double dE_dChi = RAD2DEG * c2 * d_koop[termIdx] * chi;
        let d_e_d_chi = RAD2DEG * c2 * self.koop[term_idx] * chi;
        // RDKit✔️✔️:   RDGeom::Point3D t1 = rJL.crossProduct(rJK);
        let t1 = r_jl.cross_product(r_jk);
        // RDKit✔️✔️:   RDGeom::Point3D t2 = rJI.crossProduct(rJL);
        let t2 = r_ji.cross_product(r_jl);
        // RDKit✔️✔️:   RDGeom::Point3D t3 = rJK.crossProduct(rJI);
        let t3 = r_jk.cross_product(r_ji);
        // RDKit✔️✔️:   double term1 = cosChi * sinTheta;
        let term1 = cos_chi * sin_theta;
        // RDKit✔️✔️:   double term2 = sinChi / (cosChi * sinThetaSq);
        let term2 = sin_chi / (cos_chi * sin_theta_sq);
        // RDKit✔️✔️:   double tg1[3] = {(t1.x / term1 - (rJI.x - rJK.x * cosTheta) * term2) / dJI,
        // RDKit✔️✔️:                    (t1.y / term1 - (rJI.y - rJK.y * cosTheta) * term2) / dJI,
        // RDKit✔️✔️:                    (t1.z / term1 - (rJI.z - rJK.z * cosTheta) * term2) / dJI};
        let tg1 = [
            (t1.x / term1 - (r_ji.x - r_jk.x * cos_theta) * term2) / d_ji,
            (t1.y / term1 - (r_ji.y - r_jk.y * cos_theta) * term2) / d_ji,
            (t1.z / term1 - (r_ji.z - r_jk.z * cos_theta) * term2) / d_ji,
        ];
        // RDKit✔️✔️:   double tg3[3] = {(t2.x / term1 - (rJK.x - rJI.x * cosTheta) * term2) / dJK,
        // RDKit✔️✔️:                    (t2.y / term1 - (rJK.y - rJI.y * cosTheta) * term2) / dJK,
        // RDKit✔️✔️:                    (t2.z / term1 - (rJK.z - rJI.z * cosTheta) * term2) / dJK};
        let tg3 = [
            (t2.x / term1 - (r_jk.x - r_ji.x * cos_theta) * term2) / d_jk,
            (t2.y / term1 - (r_jk.y - r_ji.y * cos_theta) * term2) / d_jk,
            (t2.z / term1 - (r_jk.z - r_ji.z * cos_theta) * term2) / d_jk,
        ];
        // RDKit✔️✔️:   double tg4[3] = {(t3.x / term1 - rJL.x * sinChi / cosChi) / dJL,
        // RDKit✔️✔️:                    (t3.y / term1 - rJL.y * sinChi / cosChi) / dJL,
        // RDKit✔️✔️:                    (t3.z / term1 - rJL.z * sinChi / cosChi) / dJL};
        let tg4 = [
            (t3.x / term1 - r_jl.x * sin_chi / cos_chi) / d_jl,
            (t3.y / term1 - r_jl.y * sin_chi / cos_chi) / d_jl,
            (t3.z / term1 - r_jl.z * sin_chi / cos_chi) / d_jl,
        ];
        // RDKit✔️✔️:   for (unsigned int i = 0; i < 3; ++i) {
        for i in 0..3 {
            // RDKit✔️✔️:     g1[i] += dE_dChi * tg1[i];
            grad[atom1_offset + i] += d_e_d_chi * tg1[i];
            // RDKit✔️✔️:     g2[i] += -dE_dChi * (tg1[i] + tg3[i] + tg4[i]);
            grad[atom2_offset + i] += -d_e_d_chi * (tg1[i] + tg3[i] + tg4[i]);
            // RDKit✔️✔️:     g3[i] += dE_dChi * tg3[i];
            grad[atom3_offset + i] += d_e_d_chi * tg3[i];
            // RDKit✔️✔️:     g4[i] += dE_dChi * tg4[i];
            grad[atom4_offset + i] += d_e_d_chi * tg4[i];
            // RDKit✔️✔️:   }
        }
        // RDKit✔️✔️: }
    }
}

impl ForceFieldContrib for OopBendContrib {
    fn copy(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }

    fn set_force_field(&mut self, owner: *const ForceField) {
        self.owner = owner;
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
    use super::OopBendContrib;
    use crate::chemistry::forcefield::core::{ForceField, ForceFieldVec3};
    use crate::chemistry::forcefield::mmff::params::MmffOop;

    const TEST_EPS: f64 = 1.0e-7;

    fn force_field_with_positions(count: usize) -> ForceField {
        let mut force_field = ForceField::new(count);
        for idx in 0..count {
            force_field
                .positions_mut()
                .push(ForceFieldVec3::new(idx as f64, idx as f64 + 0.25, idx as f64 + 0.5));
        }
        force_field
    }

    fn assert_close(actual: f64, expected: f64) {
        let scale = actual.abs().max(expected.abs()).max(1.0);
        assert!(
            (actual - expected).abs() <= TEST_EPS * scale,
            "actual={actual} expected={expected}"
        );
    }

    fn source_oop_chi(
        i_point: ForceFieldVec3,
        j_point: ForceFieldVec3,
        k_point: ForceFieldVec3,
        l_point: ForceFieldVec3,
    ) -> f64 {
        let mut r_ji = i_point - j_point;
        let mut r_jk = k_point - j_point;
        let mut r_jl = l_point - j_point;
        r_ji /= r_ji.length();
        r_jk /= r_jk.length();
        r_jl /= r_jl.length();

        let mut n = r_ji.cross_product(r_jk);
        n /= n.length();
        let mut sin_chi = n.dot_product(r_jl);
        sin_chi = sin_chi.clamp(-1.0, 1.0);

        (180.0 / std::f64::consts::PI) * sin_chi.asin()
    }

    fn source_oop_bend_energy(chi: f64, koop: f64) -> f64 {
        let c2 = 143.9325 * (std::f64::consts::PI / 180.0) * (std::f64::consts::PI / 180.0);
        0.5 * c2 * koop * chi * chi
    }

    fn assert_slice_close(actual: &[f64], expected: &[f64]) {
        assert_eq!(actual.len(), expected.len());
        for (actual_value, expected_value) in actual.iter().zip(expected.iter()) {
            assert_close(*actual_value, *expected_value);
        }
    }

    fn finite_difference_gradient(contrib: &OopBendContrib, pos: &[f64]) -> Vec<f64> {
        let mut gradient = vec![0.0; pos.len()];
        let step = 1.0e-6;

        for i in 0..pos.len() {
            let mut plus = pos.to_vec();
            plus[i] += step;
            let e_plus = contrib.get_energy(&plus);

            let mut minus = pos.to_vec();
            minus[i] -= step;
            let e_minus = contrib.get_energy(&minus);

            gradient[i] = (e_plus - e_minus) / (2.0 * step);
        }

        gradient
    }

    #[test]
    fn mmff_oopbendcontrib_constructor_stores_owner_pointer() {
        let force_field = ForceField::new(3);

        let contrib = OopBendContrib::new(&force_field);

        assert_eq!(contrib.owner(), &force_field as *const ForceField);
    }

    #[test]
    fn mmff_oopbendcontrib_constructor_initializes_no_terms() {
        let force_field = ForceField::new(3);

        let contrib = OopBendContrib::new(&force_field);

        assert_eq!(contrib.len(), 0);
        assert!(contrib.is_empty());
    }

    #[test]
    fn mmff_oopbendcontrib_constructor_accepts_empty_force_field_like_rdkit() {
        let force_field = ForceField::new(3);

        let contrib = OopBendContrib::new(&force_field);

        assert_eq!(force_field.positions().len(), 0);
        assert!(contrib.is_empty());
    }

    #[test]
    fn mmff_oopbendcontrib_add_term_pushes_source_fields() {
        let force_field = force_field_with_positions(4);
        let mut contrib = OopBendContrib::new(&force_field);
        let oop = MmffOop { koop: 0.045 };

        contrib.add_term(0, 1, 2, 3, &oop);

        assert_eq!(contrib.atom1_indices(), &[0]);
        assert_eq!(contrib.atom2_indices(), &[1]);
        assert_eq!(contrib.atom3_indices(), &[2]);
        assert_eq!(contrib.atom4_indices(), &[3]);
        assert_eq!(contrib.koop(), &[0.045]);
    }

    #[test]
    fn mmff_oopbendcontrib_add_term_appends_multiple_terms() {
        let force_field = force_field_with_positions(5);
        let mut contrib = OopBendContrib::new(&force_field);
        let first = MmffOop { koop: 0.045 };
        let second = MmffOop { koop: -0.020 };

        contrib.add_term(0, 1, 2, 3, &first);
        contrib.add_term(1, 2, 3, 4, &second);

        assert_eq!(contrib.atom1_indices(), &[0, 1]);
        assert_eq!(contrib.atom2_indices(), &[1, 2]);
        assert_eq!(contrib.atom3_indices(), &[2, 3]);
        assert_eq!(contrib.atom4_indices(), &[3, 4]);
        assert_eq!(contrib.koop(), &[0.045, -0.020]);
    }

    #[test]
    #[should_panic(expected = "degenerate points")]
    fn mmff_oopbendcontrib_add_term_rejects_first_second_degenerate() {
        let force_field = force_field_with_positions(4);
        let mut contrib = OopBendContrib::new(&force_field);
        let oop = MmffOop { koop: 0.045 };

        contrib.add_term(0, 0, 2, 3, &oop);
    }

    #[test]
    #[should_panic(expected = "degenerate points")]
    fn mmff_oopbendcontrib_add_term_rejects_first_third_degenerate() {
        let force_field = force_field_with_positions(4);
        let mut contrib = OopBendContrib::new(&force_field);
        let oop = MmffOop { koop: 0.045 };

        contrib.add_term(0, 1, 0, 3, &oop);
    }

    #[test]
    #[should_panic(expected = "degenerate points")]
    fn mmff_oopbendcontrib_add_term_rejects_first_fourth_degenerate() {
        let force_field = force_field_with_positions(4);
        let mut contrib = OopBendContrib::new(&force_field);
        let oop = MmffOop { koop: 0.045 };

        contrib.add_term(0, 1, 2, 0, &oop);
    }

    #[test]
    #[should_panic(expected = "degenerate points")]
    fn mmff_oopbendcontrib_add_term_rejects_second_third_degenerate() {
        let force_field = force_field_with_positions(4);
        let mut contrib = OopBendContrib::new(&force_field);
        let oop = MmffOop { koop: 0.045 };

        contrib.add_term(0, 1, 1, 3, &oop);
    }

    #[test]
    #[should_panic(expected = "degenerate points")]
    fn mmff_oopbendcontrib_add_term_rejects_second_fourth_degenerate() {
        let force_field = force_field_with_positions(4);
        let mut contrib = OopBendContrib::new(&force_field);
        let oop = MmffOop { koop: 0.045 };

        contrib.add_term(0, 1, 2, 1, &oop);
    }

    #[test]
    #[should_panic(expected = "degenerate points")]
    fn mmff_oopbendcontrib_add_term_rejects_third_fourth_degenerate() {
        let force_field = force_field_with_positions(4);
        let mut contrib = OopBendContrib::new(&force_field);
        let oop = MmffOop { koop: 0.045 };

        contrib.add_term(0, 1, 2, 2, &oop);
    }

    #[test]
    #[should_panic]
    fn mmff_oopbendcontrib_add_term_rejects_first_index_out_of_range() {
        let force_field = force_field_with_positions(4);
        let mut contrib = OopBendContrib::new(&force_field);
        let oop = MmffOop { koop: 0.045 };

        contrib.add_term(4, 1, 2, 3, &oop);
    }

    #[test]
    #[should_panic]
    fn mmff_oopbendcontrib_add_term_rejects_second_index_out_of_range() {
        let force_field = force_field_with_positions(4);
        let mut contrib = OopBendContrib::new(&force_field);
        let oop = MmffOop { koop: 0.045 };

        contrib.add_term(0, 4, 2, 3, &oop);
    }

    #[test]
    #[should_panic]
    fn mmff_oopbendcontrib_add_term_rejects_third_index_out_of_range() {
        let force_field = force_field_with_positions(4);
        let mut contrib = OopBendContrib::new(&force_field);
        let oop = MmffOop { koop: 0.045 };

        contrib.add_term(0, 1, 4, 3, &oop);
    }

    #[test]
    #[should_panic]
    fn mmff_oopbendcontrib_add_term_rejects_fourth_index_out_of_range() {
        let force_field = force_field_with_positions(4);
        let mut contrib = OopBendContrib::new(&force_field);
        let oop = MmffOop { koop: 0.045 };

        contrib.add_term(0, 1, 2, 4, &oop);
    }

    #[test]
    fn mmff_oopbendcontrib_get_energy_returns_zero_without_terms() {
        let force_field = ForceField::new(4);
        let contrib = OopBendContrib::new(&force_field);

        assert_eq!(contrib.get_energy(&[]), 0.0);
    }

    #[test]
    fn mmff_oopbendcontrib_get_energy_matches_source_formula_for_single_term() {
        let force_field = force_field_with_positions(4);
        let mut contrib = OopBendContrib::new(&force_field);
        let oop = MmffOop { koop: 0.045 };
        contrib.add_term(0, 1, 2, 3, &oop);
        let pos = [
            1.0, 0.0, 0.0, //
            0.0, 0.0, 0.0, //
            0.0, 1.0, 0.0, //
            0.0, 0.0, 1.0, //
        ];

        let chi = source_oop_chi(
            ForceFieldVec3::new(1.0, 0.0, 0.0),
            ForceFieldVec3::new(0.0, 0.0, 0.0),
            ForceFieldVec3::new(0.0, 1.0, 0.0),
            ForceFieldVec3::new(0.0, 0.0, 1.0),
        );
        let expected = source_oop_bend_energy(chi, 0.045);

        assert_close(contrib.get_energy(&pos), expected);
    }

    #[test]
    fn mmff_oopbendcontrib_get_energy_accumulates_multiple_terms() {
        let force_field = force_field_with_positions(5);
        let mut contrib = OopBendContrib::new(&force_field);
        contrib.add_term(0, 1, 2, 3, &MmffOop { koop: 0.045 });
        contrib.add_term(1, 2, 3, 4, &MmffOop { koop: -0.020 });
        let pos = [
            1.0, 0.0, 0.0, //
            0.0, 0.0, 0.0, //
            0.0, 1.0, 0.0, //
            0.0, 0.0, 1.0, //
            1.0, 1.0, 1.0, //
        ];

        let first = source_oop_bend_energy(
            source_oop_chi(
                ForceFieldVec3::new(1.0, 0.0, 0.0),
                ForceFieldVec3::new(0.0, 0.0, 0.0),
                ForceFieldVec3::new(0.0, 1.0, 0.0),
                ForceFieldVec3::new(0.0, 0.0, 1.0),
            ),
            0.045,
        );
        let second = source_oop_bend_energy(
            source_oop_chi(
                ForceFieldVec3::new(0.0, 0.0, 0.0),
                ForceFieldVec3::new(0.0, 1.0, 0.0),
                ForceFieldVec3::new(0.0, 0.0, 1.0),
                ForceFieldVec3::new(1.0, 1.0, 1.0),
            ),
            -0.020,
        );

        assert_close(contrib.get_energy(&pos), first + second);
    }

    #[test]
    fn mmff_oopbendcontrib_get_grad_leaves_gradient_without_terms() {
        let force_field = ForceField::new(4);
        let contrib = OopBendContrib::new(&force_field);
        let mut grad = vec![1.25; 12];

        contrib.get_grad(&[], &mut grad);

        assert_eq!(grad, vec![1.25; 12]);
    }

    #[test]
    fn mmff_oopbendcontrib_get_grad_matches_source_energy_derivative() {
        let force_field = force_field_with_positions(4);
        let mut contrib = OopBendContrib::new(&force_field);
        let oop = MmffOop { koop: 0.045 };
        let pos = [
            1.0, 0.0, 0.0, //
            0.0, 0.0, 0.0, //
            0.0, 1.0, 0.0, //
            0.0, 0.0, 1.0, //
        ];
        let mut grad = vec![0.0; pos.len()];

        contrib.add_term(0, 1, 2, 3, &oop);
        contrib.get_grad(&pos, &mut grad);

        let expected = finite_difference_gradient(&contrib, &pos);
        assert_slice_close(&grad, &expected);
    }

    #[test]
    fn mmff_oopbendcontrib_get_grad_accumulates_multiple_terms() {
        let force_field = force_field_with_positions(5);
        let mut contrib = OopBendContrib::new(&force_field);
        let pos = [
            1.0, 0.0, 0.0, //
            0.0, 0.0, 0.0, //
            0.0, 1.0, 0.0, //
            0.0, 0.0, 1.0, //
            1.0, 1.0, 1.0, //
        ];
        let mut grad = vec![0.0; pos.len()];

        contrib.add_term(0, 1, 2, 3, &MmffOop { koop: 0.045 });
        contrib.add_term(1, 2, 3, 4, &MmffOop { koop: -0.020 });
        contrib.get_grad(&pos, &mut grad);

        let expected = finite_difference_gradient(&contrib, &pos);
        assert_slice_close(&grad, &expected);
    }

    #[test]
    fn mmff_oopbendcontrib_get_grad_returns_without_changes_for_zero_length_bond_vectors() {
        let force_field = force_field_with_positions(4);
        let mut contrib = OopBendContrib::new(&force_field);
        let pos = [
            0.0, 0.0, 0.0, //
            0.0, 0.0, 0.0, //
            0.0, 1.0, 0.0, //
            0.0, 0.0, 1.0, //
        ];
        let mut grad = vec![0.0; pos.len()];

        contrib.add_term(0, 1, 2, 3, &MmffOop { koop: 0.045 });
        contrib.get_grad(&pos, &mut grad);

        assert_eq!(grad, vec![0.0; pos.len()]);
    }

    #[test]
    fn mmff_oopbendcontrib_get_grad_adds_to_existing_gradient() {
        let force_field = force_field_with_positions(4);
        let mut contrib = OopBendContrib::new(&force_field);
        let oop = MmffOop { koop: 0.045 };
        let pos = [
            1.0, 0.0, 0.0, //
            0.0, 0.0, 0.0, //
            0.0, 1.0, 0.0, //
            0.0, 0.0, 1.0, //
        ];
        let mut grad = vec![0.5; pos.len()];

        contrib.add_term(0, 1, 2, 3, &oop);
        contrib.get_grad(&pos, &mut grad);

        let expected_delta = finite_difference_gradient(&contrib, &pos);
        let expected: Vec<f64> = expected_delta.into_iter().map(|value| value + 0.5).collect();
        assert_slice_close(&grad, &expected);
    }
}
