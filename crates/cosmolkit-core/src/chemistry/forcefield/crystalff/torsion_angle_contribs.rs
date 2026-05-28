//! Source-backed RDKit CrystalFF torsion-angle contribution collection.

use crate::chemistry::forcefield::core::{
    ForceField, ForceFieldContrib, ForceFieldVec3, compute_dihedral_from_flat,
};

fn is_double_zero(value: f64) -> bool {
    // BEGIN RDKIT CPP HELPER ForceFields::MMFF::isDoubleZero (Params.h via TorsionAngleContribs.cpp)
    // RDKit✔️❌: inline bool isDoubleZero(const double x) {
    // RDKit✔️❌:   return ((x < 1.0e-10) && (x > -1.0e-10));
    // RDKit✔️❌: }
    value < 1.0e-10 && value > -1.0e-10
}

fn calc_torsion_grad(
    r: [ForceFieldVec3; 4],
    t: [ForceFieldVec3; 2],
    d: [f64; 2],
    grad: &mut [f64],
    atom_indices: [usize; 4],
    sin_term: f64,
    cos_phi: f64,
) {
    // BEGIN RDKIT CPP FUNCTION ForceFields::MMFF::Utils::calcTorsionGrad (TorsionAngle.cpp:51-95 via TorsionAngleContribs.cpp)
    // RDKit✔️❌: void calcTorsionGrad(RDGeom::Point3D *r, RDGeom::Point3D *t, double *d,
    // RDKit✔️❌:                      double **g, double &sinTerm, double &cosPhi) {
    // RDKit✔️❌:   double dCos_dT[6] = {1.0 / d[0] * (t[1].x - cosPhi * t[0].x),
    // RDKit✔️❌:                        1.0 / d[0] * (t[1].y - cosPhi * t[0].y),
    // RDKit✔️❌:                        1.0 / d[0] * (t[1].z - cosPhi * t[0].z),
    // RDKit✔️❌:                        1.0 / d[1] * (t[0].x - cosPhi * t[1].x),
    // RDKit✔️❌:                        1.0 / d[1] * (t[0].y - cosPhi * t[1].y),
    // RDKit✔️❌:                        1.0 / d[1] * (t[0].z - cosPhi * t[1].z)};
    let d_cos_dt = [
        (t[1].x - cos_phi * t[0].x) / d[0],
        (t[1].y - cos_phi * t[0].y) / d[0],
        (t[1].z - cos_phi * t[0].z) / d[0],
        (t[0].x - cos_phi * t[1].x) / d[1],
        (t[0].y - cos_phi * t[1].y) / d[1],
        (t[0].z - cos_phi * t[1].z) / d[1],
    ];
    let [atom1_idx, atom2_idx, atom3_idx, atom4_idx] = atom_indices;
    let atom1_offset = 3 * atom1_idx;
    let atom2_offset = 3 * atom2_idx;
    let atom3_offset = 3 * atom3_idx;
    let atom4_offset = 3 * atom4_idx;

    // RDKit✔️❌:   g[0][0] += sinTerm * (dCos_dT[2] * r[1].y - dCos_dT[1] * r[1].z);
    grad[atom1_offset] += sin_term * (d_cos_dt[2] * r[1].y - d_cos_dt[1] * r[1].z);
    // RDKit✔️❌:   g[0][1] += sinTerm * (dCos_dT[0] * r[1].z - dCos_dT[2] * r[1].x);
    grad[atom1_offset + 1] += sin_term * (d_cos_dt[0] * r[1].z - d_cos_dt[2] * r[1].x);
    // RDKit✔️❌:   g[0][2] += sinTerm * (dCos_dT[1] * r[1].x - dCos_dT[0] * r[1].y);
    grad[atom1_offset + 2] += sin_term * (d_cos_dt[1] * r[1].x - d_cos_dt[0] * r[1].y);

    // RDKit✔️❌:   g[1][0] += sinTerm *
    grad[atom2_offset] += sin_term
        * ((d_cos_dt[1] * (r[1].z - r[0].z))
            + (d_cos_dt[2] * (r[0].y - r[1].y))
            + (d_cos_dt[4] * (-r[3].z))
            + (d_cos_dt[5] * r[3].y));
    // RDKit✔️❌:   g[1][1] += sinTerm *
    grad[atom2_offset + 1] += sin_term
        * ((d_cos_dt[0] * (r[0].z - r[1].z))
            + (d_cos_dt[2] * (r[1].x - r[0].x))
            + (d_cos_dt[3] * r[3].z)
            + (d_cos_dt[5] * (-r[3].x)));
    // RDKit✔️❌:   g[1][2] += sinTerm *
    grad[atom2_offset + 2] += sin_term
        * ((d_cos_dt[0] * (r[1].y - r[0].y))
            + (d_cos_dt[1] * (r[0].x - r[1].x))
            + (d_cos_dt[3] * (-r[3].y))
            + (d_cos_dt[4] * r[3].x));

    // RDKit✔️❌:   g[2][0] += sinTerm *
    grad[atom3_offset] += sin_term
        * ((d_cos_dt[1] * r[0].z)
            + (d_cos_dt[2] * (-r[0].y))
            + (d_cos_dt[4] * (r[3].z - r[2].z))
            + (d_cos_dt[5] * (r[2].y - r[3].y)));
    // RDKit✔️❌:   g[2][1] += sinTerm *
    grad[atom3_offset + 1] += sin_term
        * ((d_cos_dt[0] * (-r[0].z))
            + (d_cos_dt[2] * r[0].x)
            + (d_cos_dt[3] * (r[2].z - r[3].z))
            + (d_cos_dt[5] * (r[3].x - r[2].x)));
    // RDKit✔️❌:   g[2][2] += sinTerm *
    grad[atom3_offset + 2] += sin_term
        * ((d_cos_dt[0] * r[0].y)
            + (d_cos_dt[1] * (-r[0].x))
            + (d_cos_dt[3] * (r[3].y - r[2].y))
            + (d_cos_dt[4] * (r[2].x - r[3].x)));

    // RDKit✔️❌:   g[3][0] += sinTerm * (dCos_dT[4] * r[2].z - dCos_dT[5] * r[2].y);
    grad[atom4_offset] += sin_term * (d_cos_dt[4] * r[2].z - d_cos_dt[5] * r[2].y);
    // RDKit✔️❌:   g[3][1] += sinTerm * (dCos_dT[5] * r[2].x - dCos_dT[3] * r[2].z);
    grad[atom4_offset + 1] += sin_term * (d_cos_dt[5] * r[2].x - d_cos_dt[3] * r[2].z);
    // RDKit✔️❌:   g[3][2] += sinTerm * (dCos_dT[3] * r[2].y - dCos_dT[4] * r[2].x);
    grad[atom4_offset + 2] += sin_term * (d_cos_dt[3] * r[2].y - d_cos_dt[4] * r[2].x);
    // RDKit✔️❌: }
}

#[derive(Clone, Debug, PartialEq)]
pub struct TorsionAngleContribsParams {
    pub idx1: usize,
    pub idx2: usize,
    pub idx3: usize,
    pub idx4: usize,
    pub force_constants: Vec<f64>,
    pub signs: Vec<i32>,
}

impl TorsionAngleContribsParams {
    #[must_use]
    pub fn new(
        idx1: usize,
        idx2: usize,
        idx3: usize,
        idx4: usize,
        force_constants: Vec<f64>,
        signs: Vec<i32>,
    ) -> Self {
        Self {
            idx1,
            idx2,
            idx3,
            idx4,
            force_constants,
            signs,
        }
    }
}

#[must_use]
pub fn calc_torsion_energy(force_constants: &[f64], signs: &[i32], cos_phi: f64) -> f64 {
    // BEGIN RDKIT CPP FUNCTION ForceFields::CrystalFF::calcTorsionEnergy (TorsionAngleContribs.h:102-104, TorsionAngleContribs.cpp:20-41)
    // RDKit✔️❌: double calcTorsionEnergyM6(const std::vector<double> &forceConstants,
    // RDKit✔️❌:                            const std::vector<int> &signs, const double cosPhi) {
    // RDKit✔️❌:   const double cosPhi2 = cosPhi * cosPhi;
    let cos_phi2 = cos_phi * cos_phi;
    // RDKit✔️❌:   const double cosPhi3 = cosPhi * cosPhi2;
    let cos_phi3 = cos_phi * cos_phi2;
    // RDKit✔️❌:   const double cosPhi4 = cosPhi * cosPhi3;
    let cos_phi4 = cos_phi * cos_phi3;
    // RDKit✔️❌:   const double cosPhi5 = cosPhi * cosPhi4;
    let cos_phi5 = cos_phi * cos_phi4;
    // RDKit✔️❌:   const double cosPhi6 = cosPhi * cosPhi5;
    let cos_phi6 = cos_phi * cos_phi5;

    // RDKit✔️❌:   const double cos2Phi = 2.0 * cosPhi2 - 1.0;
    let cos2_phi = 2.0 * cos_phi2 - 1.0;
    // RDKit✔️❌:   const double cos3Phi = 4.0 * cosPhi3 - 3.0 * cosPhi;
    let cos3_phi = 4.0 * cos_phi3 - 3.0 * cos_phi;
    // RDKit✔️❌:   const double cos4Phi = 8.0 * cosPhi4 - 8.0 * cosPhi2 + 1.0;
    let cos4_phi = 8.0 * cos_phi4 - 8.0 * cos_phi2 + 1.0;
    // RDKit✔️❌:   const double cos5Phi = 16.0 * cosPhi5 - 20.0 * cosPhi3 + 5.0 * cosPhi;
    let cos5_phi = 16.0 * cos_phi5 - 20.0 * cos_phi3 + 5.0 * cos_phi;
    // RDKit✔️❌:   const double cos6Phi = 32.0 * cosPhi6 - 48.0 * cosPhi4 + 18.0 * cosPhi2 - 1.0;
    let cos6_phi = 32.0 * cos_phi6 - 48.0 * cos_phi4 + 18.0 * cos_phi2 - 1.0;

    // RDKit✔️❌:   return (forceConstants[0] * (1.0 + signs[0] * cosPhi) +
    // RDKit✔️❌:           forceConstants[1] * (1.0 + signs[1] * cos2Phi) +
    // RDKit✔️❌:           forceConstants[2] * (1.0 + signs[2] * cos3Phi) +
    // RDKit✔️❌:           forceConstants[3] * (1.0 + signs[3] * cos4Phi) +
    // RDKit✔️❌:           forceConstants[4] * (1.0 + signs[4] * cos5Phi) +
    // RDKit✔️❌:           forceConstants[5] * (1.0 + signs[5] * cos6Phi));
    // RDKit✔️❌: }
    force_constants[0] * (1.0 + f64::from(signs[0]) * cos_phi)
        + force_constants[1] * (1.0 + f64::from(signs[1]) * cos2_phi)
        + force_constants[2] * (1.0 + f64::from(signs[2]) * cos3_phi)
        + force_constants[3] * (1.0 + f64::from(signs[3]) * cos4_phi)
        + force_constants[4] * (1.0 + f64::from(signs[4]) * cos5_phi)
        + force_constants[5] * (1.0 + f64::from(signs[5]) * cos6_phi)
}

#[derive(Clone, Debug)]
pub struct TorsionAngleContribs {
    owner: *const ForceField,
    contribs: Vec<TorsionAngleContribsParams>,
}

impl TorsionAngleContribs {
    #[must_use]
    pub fn new(owner: &ForceField) -> Self {
        // BEGIN RDKIT CPP CONSTRUCTOR ForceFields::CrystalFF::TorsionAngleContribs::TorsionAngleContribs (TorsionAngleContribs.cpp:44-47)
        // RDKit✔️❌: TorsionAngleContribs::TorsionAngleContribs(ForceField *owner) {
        // RDKit✔️❌:   PRECONDITION(owner, "bad owner");
        // Rust references reproduce RDKit's non-null owner precondition.
        // RDKit✔️❌:   dp_forceField = owner;
        // RDKit✔️❌: }
        Self {
            owner: owner as *const ForceField,
            contribs: Vec::new(),
        }
    }

    pub fn add_contrib(
        &mut self,
        idx1: usize,
        idx2: usize,
        idx3: usize,
        idx4: usize,
        force_constants: Vec<f64>,
        signs: Vec<i32>,
    ) {
        // BEGIN RDKIT CPP METHOD ForceFields::CrystalFF::TorsionAngleContribs::addContrib (TorsionAngleContribs.cpp:49-61)
        // RDKit✔️❌: PRECONDITION((idx1 != idx2) && (idx1 != idx3) && (idx1 != idx4) &&
        // RDKit✔️❌:                (idx2 != idx3) && (idx2 != idx4) && (idx3 != idx4),
        // RDKit✔️❌:              "degenerate points");
        assert!(
            idx1 != idx2
                && idx1 != idx3
                && idx1 != idx4
                && idx2 != idx3
                && idx2 != idx4
                && idx3 != idx4,
            "degenerate points"
        );
        // RDKit✔️❌: URANGE_CHECK(idx1, dp_forceField->positions().size());
        // RDKit✔️❌: URANGE_CHECK(idx2, dp_forceField->positions().size());
        // RDKit✔️❌: URANGE_CHECK(idx3, dp_forceField->positions().size());
        // RDKit✔️❌: URANGE_CHECK(idx4, dp_forceField->positions().size());
        let owner = unsafe { self.owner.as_ref() }.expect("bad owner");
        let positions_len = owner.positions().len();
        assert!(idx1 < positions_len);
        assert!(idx2 < positions_len);
        assert!(idx3 < positions_len);
        assert!(idx4 < positions_len);
        // RDKit✔️❌: d_contribs.emplace_back(idx1, idx2, idx3, idx4, std::move(forceConstants),
        // RDKit✔️❌:                         std::move(signs));
        self.contribs.push(TorsionAngleContribsParams::new(
            idx1,
            idx2,
            idx3,
            idx4,
            force_constants,
            signs,
        ));
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.contribs.is_empty()
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.contribs.len()
    }
}

impl ForceFieldContrib for TorsionAngleContribs {
    fn copy(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }

    fn set_force_field(&mut self, owner: *const ForceField) {
        self.owner = owner;
    }

    fn get_energy(&self, pos: &[f64]) -> f64 {
        // BEGIN RDKIT CPP METHOD ForceFields::CrystalFF::TorsionAngleContribs::getEnergy (TorsionAngleContribs.cpp:66-88)
        // RDKit✔️❌: PRECONDITION(dp_forceField, "no owner");
        assert!(!self.owner.is_null(), "no owner");
        // RDKit✔️❌: PRECONDITION(pos, "bad vector");
        assert!(!pos.is_empty(), "bad vector");
        let owner = unsafe { self.owner.as_ref() }.expect("no owner");
        assert!(
            pos.len() >= 3 * owner.positions().len(),
            "bad vector length for force-field positions"
        );
        // RDKit✔️❌: double accum = 0.0;
        let mut accum = 0.0;
        // RDKit✔️❌: for (const auto &contrib : d_contribs) {
        for contrib in &self.contribs {
            // RDKit✔️❌:   const RDGeom::Point3D iPoint(pos[3 * contrib.idx1],
            // RDKit✔️❌:                              pos[3 * contrib.idx1 + 1],
            // RDKit✔️❌:                              pos[3 * contrib.idx1 + 2]);
            // RDKit✔️❌:   const RDGeom::Point3D jPoint(pos[3 * contrib.idx2],
            // RDKit✔️❌:                              pos[3 * contrib.idx2 + 1],
            // RDKit✔️❌:                              pos[3 * contrib.idx2 + 2]);
            // RDKit✔️❌:   const RDGeom::Point3D kPoint(pos[3 * contrib.idx3],
            // RDKit✔️❌:                              pos[3 * contrib.idx3 + 1],
            // RDKit✔️❌:                              pos[3 * contrib.idx3 + 2]);
            // RDKit✔️❌:   const RDGeom::Point3D lPoint(pos[3 * contrib.idx4],
            // RDKit✔️❌:                              pos[3 * contrib.idx4 + 1],
            // RDKit✔️❌:                              pos[3 * contrib.idx4 + 2]);
            // RDKit✔️❌:   accum += calcTorsionEnergyM6(
            // RDKit✔️❌:       contrib.forceConstants, contrib.signs,
            // RDKit✔️❌:       MMFF::Utils::calcTorsionCosPhi(iPoint, jPoint, kPoint, lPoint));
            let cos_phi = compute_dihedral_from_flat(
                pos,
                contrib.idx1,
                contrib.idx2,
                contrib.idx3,
                contrib.idx4,
                false,
            )
            .cos_phi;
            accum += calc_torsion_energy(&contrib.force_constants, &contrib.signs, cos_phi);
        }
        // RDKit✔️❌: return accum;
        // RDKit✔️❌: }
        accum
    }

    fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
        // BEGIN RDKIT CPP METHOD ForceFields::CrystalFF::TorsionAngleContribs::getGrad (TorsionAngleContribs.cpp:90-150)
        // RDKit✔️❌: PRECONDITION(dp_forceField, "no owner");
        assert!(!self.owner.is_null(), "no owner");
        // RDKit✔️❌: PRECONDITION(pos, "bad vector");
        assert!(!pos.is_empty(), "bad vector");
        // RDKit✔️❌: PRECONDITION(grad, "bad vector");
        assert!(!grad.is_empty(), "bad vector");
        let owner = unsafe { self.owner.as_ref() }.expect("no owner");
        assert!(
            pos.len() >= 3 * owner.positions().len(),
            "bad vector length for force-field positions"
        );
        assert!(
            grad.len() >= 3 * owner.positions().len(),
            "bad gradient length for force-field positions"
        );

        // RDKit✔️❌: for (const auto &contrib : d_contribs) {
        for contrib in &self.contribs {
            // RDKit✔️❌:   const RDGeom::Point3D iPoint(pos[3 * contrib.idx1],
            let i_point = ForceFieldVec3::new(
                pos[3 * contrib.idx1],
                pos[3 * contrib.idx1 + 1],
                pos[3 * contrib.idx1 + 2],
            );
            // RDKit✔️❌:   const RDGeom::Point3D jPoint(pos[3 * contrib.idx2],
            let j_point = ForceFieldVec3::new(
                pos[3 * contrib.idx2],
                pos[3 * contrib.idx2 + 1],
                pos[3 * contrib.idx2 + 2],
            );
            // RDKit✔️❌:   const RDGeom::Point3D kPoint(pos[3 * contrib.idx3],
            let k_point = ForceFieldVec3::new(
                pos[3 * contrib.idx3],
                pos[3 * contrib.idx3 + 1],
                pos[3 * contrib.idx3 + 2],
            );
            // RDKit✔️❌:   const RDGeom::Point3D lPoint(pos[3 * contrib.idx4],
            let l_point = ForceFieldVec3::new(
                pos[3 * contrib.idx4],
                pos[3 * contrib.idx4 + 1],
                pos[3 * contrib.idx4 + 2],
            );

            // RDKit✔️❌:   RDGeom::Point3D r[4] = {iPoint - jPoint, kPoint - jPoint, jPoint - kPoint,
            // RDKit✔️❌:                           lPoint - kPoint};
            let r = [
                i_point - j_point,
                k_point - j_point,
                j_point - k_point,
                l_point - k_point,
            ];
            // RDKit✔️❌:   RDGeom::Point3D t[2] = {r[0].crossProduct(r[1]), r[2].crossProduct(r[3])};
            let mut t = [r[0].cross_product(r[1]), r[2].cross_product(r[3])];
            // RDKit✔️❌:   double d[2] = {t[0].length(), t[1].length()};
            let d = [t[0].length(), t[1].length()];
            // RDKit✔️❌:   if (MMFF::isDoubleZero(d[0]) || MMFF::isDoubleZero(d[1])) {
            // RDKit✔️❌:     return;
            // RDKit✔️❌:   }
            if is_double_zero(d[0]) || is_double_zero(d[1]) {
                return;
            }
            // RDKit✔️❌:   t[0] /= d[0];
            // RDKit✔️❌:   t[1] /= d[1];
            t[0] /= d[0];
            t[1] /= d[1];
            // RDKit✔️❌:   double cosPhi = t[0].dotProduct(t[1]);
            let mut cos_phi = t[0].dot_product(t[1]);
            // RDKit✔️❌:   cosPhi = std::clamp(cosPhi, -1.0, 1.0);
            cos_phi = cos_phi.clamp(-1.0, 1.0);
            // RDKit✔️❌:   const double sinPhiSq = 1.0 - cosPhi * cosPhi;
            let sin_phi_sq = 1.0 - cos_phi * cos_phi;
            // RDKit✔️❌:   const double sinPhi = ((sinPhiSq > 0.0) ? sqrt(sinPhiSq) : 0.0);
            let sin_phi = if sin_phi_sq > 0.0 {
                sin_phi_sq.sqrt()
            } else {
                0.0
            };
            // RDKit✔️❌:   const double cosPhi2 = cosPhi * cosPhi;
            let cos_phi2 = cos_phi * cos_phi;
            // RDKit✔️❌:   const double cosPhi3 = cosPhi * cosPhi2;
            let cos_phi3 = cos_phi * cos_phi2;
            // RDKit✔️❌:   const double cosPhi4 = cosPhi * cosPhi3;
            let cos_phi4 = cos_phi * cos_phi3;
            // RDKit✔️❌:   const double cosPhi5 = cosPhi * cosPhi4;
            let cos_phi5 = cos_phi * cos_phi4;
            // RDKit✔️❌:   const double dE_dPhi =
            let d_e_d_phi = -contrib.force_constants[0] * f64::from(contrib.signs[0]) * sin_phi
                - 2.0
                    * contrib.force_constants[1]
                    * f64::from(contrib.signs[1])
                    * (2.0 * cos_phi * sin_phi)
                - 3.0
                    * contrib.force_constants[2]
                    * f64::from(contrib.signs[2])
                    * (4.0 * cos_phi2 * sin_phi - sin_phi)
                - 4.0
                    * contrib.force_constants[3]
                    * f64::from(contrib.signs[3])
                    * (8.0 * cos_phi3 * sin_phi - 4.0 * cos_phi * sin_phi)
                - 5.0
                    * contrib.force_constants[4]
                    * f64::from(contrib.signs[4])
                    * (16.0 * cos_phi4 * sin_phi - 12.0 * cos_phi2 * sin_phi + sin_phi)
                - 6.0
                    * contrib.force_constants[4]
                    * f64::from(contrib.signs[4])
                    * (32.0 * cos_phi5 * sin_phi - 32.0 * cos_phi3 * sin_phi + 6.0 * sin_phi);

            // RDKit✔️❌:   double sinTerm = -dE_dPhi * (MMFF::isDoubleZero(sinPhi) ? (1.0 / cosPhi)
            // RDKit✔️❌:                                                           : (1.0 / sinPhi));
            let sin_term = -d_e_d_phi
                * if is_double_zero(sin_phi) {
                    1.0 / cos_phi
                } else {
                    1.0 / sin_phi
                };

            // RDKit✔️❌:   MMFF::Utils::calcTorsionGrad(r, t, d, g, sinTerm, cosPhi);
            calc_torsion_grad(
                r,
                t,
                d,
                grad,
                [contrib.idx1, contrib.idx2, contrib.idx3, contrib.idx4],
                sin_term,
                cos_phi,
            );
        }
        // RDKit✔️❌: }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Molecule;
    use crate::chemistry::forcefield::core::{ForceField, ForceFieldVec3};
    use crate::chemistry::forcefield::crystalff::torsion_preferences::{
        CrystalFFDetails, get_experimental_torsions_without_bonds,
    };

    const EPS_TEST: f64 = 1.0e-10;

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < EPS_TEST,
            "actual={actual} expected={expected}"
        );
    }

    fn torsion_forcefield() -> ForceField {
        let mut ff = ForceField::new(3);
        ff.positions_mut().push(ForceFieldVec3::new(0.0, 1.5, 0.0));
        ff.positions_mut().push(ForceFieldVec3::new(0.0, 0.0, 0.0));
        ff.positions_mut().push(ForceFieldVec3::new(1.5, 0.0, 0.0));
        ff.positions_mut().push(ForceFieldVec3::new(1.5, 0.0, 1.5));
        ff.initialize();
        ff
    }

    fn flattened_positions(ff: &ForceField) -> Vec<f64> {
        let mut pos = vec![0.0; 3 * ff.positions().len()];
        for (atom_idx, point) in ff.positions().iter().enumerate() {
            pos[3 * atom_idx] = point.x;
            pos[3 * atom_idx + 1] = point.y;
            pos[3 * atom_idx + 2] = point.z;
        }
        pos
    }

    fn degenerate_torsion_forcefield() -> ForceField {
        let mut ff = ForceField::new(3);
        ff.positions_mut().push(ForceFieldVec3::new(0.0, 0.0, 0.0));
        ff.positions_mut().push(ForceFieldVec3::new(1.0, 0.0, 0.0));
        ff.positions_mut().push(ForceFieldVec3::new(2.0, 0.0, 0.0));
        ff.positions_mut().push(ForceFieldVec3::new(3.0, 0.0, 0.0));
        ff.initialize();
        ff
    }

    fn forcefield_for_atom_count(atom_count: usize) -> ForceField {
        let mut ff = ForceField::new(3);
        for idx in 0..atom_count {
            ff.positions_mut().push(ForceFieldVec3::new(
                idx as f64,
                if idx % 2 == 0 { 0.0 } else { 1.0 },
                if idx % 3 == 0 { 0.5 } else { 0.0 },
            ));
        }
        ff.initialize();
        ff
    }

    fn add_details_contribs(contribs: &mut TorsionAngleContribs, details: &CrystalFFDetails) {
        for (atoms, (signs, force_constants)) in details
            .exp_torsion_atoms
            .iter()
            .zip(details.exp_torsion_angles.iter())
        {
            contribs.add_contrib(
                usize::try_from(atoms[0]).expect("non-negative atom index"),
                usize::try_from(atoms[1]).expect("non-negative atom index"),
                usize::try_from(atoms[2]).expect("non-negative atom index"),
                usize::try_from(atoms[3]).expect("non-negative atom index"),
                force_constants.clone(),
                signs.clone(),
            );
        }
    }

    #[test]
    fn crystalff_torsionanglecontribs_calc_torsion_energy_matches_m6_closed_form() {
        let force_constants = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let signs = vec![1, -1, 1, -1, 1, -1];
        let cos_phi = 0.5;

        let energy = calc_torsion_energy(&force_constants, &signs, cos_phi);

        let cos2_phi = 2.0 * cos_phi * cos_phi - 1.0;
        let cos3_phi = 4.0 * cos_phi * cos_phi * cos_phi - 3.0 * cos_phi;
        let cos4_phi = 8.0 * cos_phi.powi(4) - 8.0 * cos_phi.powi(2) + 1.0;
        let cos5_phi = 16.0 * cos_phi.powi(5) - 20.0 * cos_phi.powi(3) + 5.0 * cos_phi;
        let cos6_phi =
            32.0 * cos_phi.powi(6) - 48.0 * cos_phi.powi(4) + 18.0 * cos_phi.powi(2) - 1.0;
        let expected = force_constants[0] * (1.0 + cos_phi)
            + force_constants[1] * (1.0 - cos2_phi)
            + force_constants[2] * (1.0 + cos3_phi)
            + force_constants[3] * (1.0 - cos4_phi)
            + force_constants[4] * (1.0 + cos5_phi)
            + force_constants[5] * (1.0 - cos6_phi);

        assert_close(energy, expected);
    }

    #[test]
    fn crystalff_torsionanglecontribs_constructor_starts_empty_and_supports_m6_gradient_path() {
        let ff = torsion_forcefield();
        let pos = flattened_positions(&ff);
        let mut contribs = TorsionAngleContribs::new(&ff);

        assert!(contribs.is_empty());
        assert_eq!(contribs.len(), 0);

        let force_constants = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let signs = vec![1, -1, 1, -1, 1, -1];
        contribs.add_contrib(0, 1, 2, 3, force_constants.clone(), signs.clone());

        assert_eq!(contribs.len(), 1);

        let energy = contribs.get_energy(&pos);
        let expected_energy = calc_torsion_energy(
            &force_constants,
            &signs,
            compute_dihedral_from_flat(&pos, 0, 1, 2, 3, false).cos_phi,
        );
        let mut grad = vec![0.0; pos.len()];
        contribs.get_grad(&pos, &mut grad);

        assert_close(energy, expected_energy);
        assert!(grad.iter().all(|value| value.is_finite()));
        assert!(grad.iter().any(|value| value.abs() > EPS_TEST));
    }

    #[test]
    fn crystalff_torsionanglecontribs_constructor_accepts_smarts_small_ring_and_macrocycle_terms() {
        let linear = Molecule::from_smiles("CCCCC").expect("pentane");
        let small_ring = Molecule::from_smiles("C1COCC1").expect("tetrahydrofuran-like ring");
        let macrocycle = Molecule::from_smiles("C1COCCCCCCC1").expect("macrocycle");
        let mut linear_details = CrystalFFDetails::default();
        let mut small_ring_details = CrystalFFDetails::default();
        let mut macrocycle_details = CrystalFFDetails::default();

        get_experimental_torsions_without_bonds(
            &linear,
            &mut linear_details,
            true,
            false,
            false,
            false,
            2,
            false,
        )
        .expect("linear SMARTS torsions");
        get_experimental_torsions_without_bonds(
            &small_ring,
            &mut small_ring_details,
            true,
            true,
            false,
            false,
            1,
            false,
        )
        .expect("small-ring torsions");
        get_experimental_torsions_without_bonds(
            &macrocycle,
            &mut macrocycle_details,
            true,
            false,
            true,
            false,
            1,
            false,
        )
        .expect("macrocycle torsions");

        assert_eq!(linear_details.exp_torsion_atoms.len(), 2);
        assert!(!small_ring_details.exp_torsion_atoms.is_empty());
        assert!(!macrocycle_details.exp_torsion_atoms.is_empty());

        let linear_ff = forcefield_for_atom_count(linear.num_atoms());
        let mut linear_contribs = TorsionAngleContribs::new(&linear_ff);
        add_details_contribs(&mut linear_contribs, &linear_details);
        assert_eq!(
            linear_contribs.len(),
            linear_details.exp_torsion_atoms.len()
        );

        let small_ring_ff = forcefield_for_atom_count(small_ring.num_atoms());
        let mut small_ring_contribs = TorsionAngleContribs::new(&small_ring_ff);
        add_details_contribs(&mut small_ring_contribs, &small_ring_details);
        assert_eq!(
            small_ring_contribs.len(),
            small_ring_details.exp_torsion_atoms.len()
        );

        let macrocycle_ff = forcefield_for_atom_count(macrocycle.num_atoms());
        let mut macrocycle_contribs = TorsionAngleContribs::new(&macrocycle_ff);
        add_details_contribs(&mut macrocycle_contribs, &macrocycle_details);
        assert_eq!(
            macrocycle_contribs.len(),
            macrocycle_details.exp_torsion_atoms.len()
        );
    }

    #[test]
    fn crystalff_torsionanglecontribs_addcontrib_accepts_smarts_small_ring_and_macrocycle_terms() {
        let linear = Molecule::from_smiles("CCCCC").expect("pentane");
        let small_ring = Molecule::from_smiles("C1COCC1").expect("tetrahydrofuran-like ring");
        let macrocycle = Molecule::from_smiles("C1COCCCCCCC1").expect("macrocycle");
        let mut linear_details = CrystalFFDetails::default();
        let mut small_ring_details = CrystalFFDetails::default();
        let mut macrocycle_details = CrystalFFDetails::default();

        get_experimental_torsions_without_bonds(
            &linear,
            &mut linear_details,
            true,
            false,
            false,
            false,
            2,
            false,
        )
        .expect("linear SMARTS torsions");
        get_experimental_torsions_without_bonds(
            &small_ring,
            &mut small_ring_details,
            true,
            true,
            false,
            false,
            1,
            false,
        )
        .expect("small-ring torsions");
        get_experimental_torsions_without_bonds(
            &macrocycle,
            &mut macrocycle_details,
            true,
            false,
            true,
            false,
            1,
            false,
        )
        .expect("macrocycle torsions");

        let linear_ff = forcefield_for_atom_count(linear.num_atoms());
        let mut linear_contribs = TorsionAngleContribs::new(&linear_ff);
        add_details_contribs(&mut linear_contribs, &linear_details);
        assert_eq!(
            linear_contribs.len(),
            linear_details.exp_torsion_atoms.len()
        );

        let small_ring_ff = forcefield_for_atom_count(small_ring.num_atoms());
        let mut small_ring_contribs = TorsionAngleContribs::new(&small_ring_ff);
        add_details_contribs(&mut small_ring_contribs, &small_ring_details);
        assert_eq!(
            small_ring_contribs.len(),
            small_ring_details.exp_torsion_atoms.len()
        );

        let macrocycle_ff = forcefield_for_atom_count(macrocycle.num_atoms());
        let mut macrocycle_contribs = TorsionAngleContribs::new(&macrocycle_ff);
        add_details_contribs(&mut macrocycle_contribs, &macrocycle_details);
        assert_eq!(
            macrocycle_contribs.len(),
            macrocycle_details.exp_torsion_atoms.len()
        );
    }

    #[test]
    fn crystalff_torsionanglecontribs_addcontrib_supports_m6_energy_and_gradient_paths() {
        let ff = torsion_forcefield();
        let pos = flattened_positions(&ff);
        let mut contribs = TorsionAngleContribs::new(&ff);
        let force_constants = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let signs = vec![1, -1, 1, -1, 1, -1];

        contribs.add_contrib(0, 1, 2, 3, force_constants.clone(), signs.clone());

        assert_eq!(contribs.len(), 1);
        let energy = contribs.get_energy(&pos);
        let mut grad = vec![0.0; pos.len()];
        contribs.get_grad(&pos, &mut grad);

        assert_close(
            energy,
            calc_torsion_energy(
                &force_constants,
                &signs,
                compute_dihedral_from_flat(&pos, 0, 1, 2, 3, false).cos_phi,
            ),
        );
        assert!(grad.iter().all(|value| value.is_finite()));
        assert!(grad.iter().any(|value| value.abs() > EPS_TEST));
    }

    #[test]
    #[should_panic(expected = "degenerate points")]
    fn crystalff_torsionanglecontribs_addcontrib_panics_for_degenerate_indices() {
        let ff = torsion_forcefield();
        let mut contribs = TorsionAngleContribs::new(&ff);

        contribs.add_contrib(0, 0, 2, 3, vec![0.0; 6], vec![1; 6]);
    }

    #[test]
    #[should_panic(expected = "assertion failed: idx1 < positions_len")]
    fn crystalff_torsionanglecontribs_addcontrib_panics_for_out_of_range_first_index() {
        let ff = torsion_forcefield();
        let mut contribs = TorsionAngleContribs::new(&ff);

        contribs.add_contrib(4, 1, 2, 3, vec![0.0; 6], vec![1; 6]);
    }

    #[test]
    #[should_panic(expected = "assertion failed: idx2 < positions_len")]
    fn crystalff_torsionanglecontribs_addcontrib_panics_for_out_of_range_second_index() {
        let ff = torsion_forcefield();
        let mut contribs = TorsionAngleContribs::new(&ff);

        contribs.add_contrib(0, 4, 2, 3, vec![0.0; 6], vec![1; 6]);
    }

    #[test]
    #[should_panic(expected = "assertion failed: idx3 < positions_len")]
    fn crystalff_torsionanglecontribs_addcontrib_panics_for_out_of_range_third_index() {
        let ff = torsion_forcefield();
        let mut contribs = TorsionAngleContribs::new(&ff);

        contribs.add_contrib(0, 1, 4, 3, vec![0.0; 6], vec![1; 6]);
    }

    #[test]
    #[should_panic(expected = "assertion failed: idx4 < positions_len")]
    fn crystalff_torsionanglecontribs_addcontrib_panics_for_out_of_range_fourth_index() {
        let ff = torsion_forcefield();
        let mut contribs = TorsionAngleContribs::new(&ff);

        contribs.add_contrib(0, 1, 2, 4, vec![0.0; 6], vec![1; 6]);
    }

    #[test]
    fn crystalff_torsionanglecontribs_get_energy_matches_sp3_sp3_source_geometry() {
        let ff = torsion_forcefield();
        let pos = flattened_positions(&ff);
        let mut contrib = TorsionAngleContribs::new(&ff);
        let signs = vec![1, 1, 1, 1, 1, 1];
        let mut force_constants = vec![0.0; 6];
        force_constants[2] = 4.0;
        contrib.add_contrib(0, 1, 2, 3, force_constants, signs);

        let energy = contrib.get_energy(&pos);
        let cos_phi = compute_dihedral_from_flat(&pos, 0, 1, 2, 3, false).cos_phi;

        assert_close(cos_phi, 0.0);
        assert_close(energy, 4.0);
    }

    #[test]
    fn crystalff_torsionanglecontribs_get_energy_accumulates_multiple_contributions() {
        let ff = torsion_forcefield();
        let pos = flattened_positions(&ff);
        let mut contrib = TorsionAngleContribs::new(&ff);
        let mut first = vec![0.0; 6];
        first[2] = 4.0;
        let mut second = vec![0.0; 6];
        second[0] = 1.5;
        second[1] = 0.5;
        contrib.add_contrib(0, 1, 2, 3, first.clone(), vec![1; 6]);
        contrib.add_contrib(0, 1, 2, 3, second.clone(), vec![1, -1, 1, 1, 1, 1]);

        let energy = contrib.get_energy(&pos);
        let expected = calc_torsion_energy(&first, &[1, 1, 1, 1, 1, 1], 0.0)
            + calc_torsion_energy(&second, &[1, -1, 1, 1, 1, 1], 0.0);

        assert_close(energy, expected);
    }

    #[test]
    fn crystalff_torsionanglecontribs_get_energy_returns_zero_when_empty() {
        let ff = torsion_forcefield();
        let pos = flattened_positions(&ff);
        let contrib = TorsionAngleContribs::new(&ff);

        let energy = contrib.get_energy(&pos);

        assert_close(energy, 0.0);
    }

    #[test]
    #[should_panic(expected = "no owner")]
    fn crystalff_torsionanglecontribs_get_energy_panics_without_owner() {
        let ff = torsion_forcefield();
        let pos = flattened_positions(&ff);
        let mut contrib = TorsionAngleContribs::new(&ff);
        contrib.set_force_field(std::ptr::null());

        let _ = contrib.get_energy(&pos);
    }

    #[test]
    #[should_panic(expected = "bad vector")]
    fn crystalff_torsionanglecontribs_get_energy_panics_for_empty_position_vector() {
        let ff = torsion_forcefield();
        let contrib = TorsionAngleContribs::new(&ff);

        let _ = contrib.get_energy(&[]);
    }

    #[test]
    #[should_panic(expected = "bad vector length for force-field positions")]
    fn crystalff_torsionanglecontribs_get_energy_panics_for_short_position_vector() {
        let ff = torsion_forcefield();
        let contrib = TorsionAngleContribs::new(&ff);

        let _ = contrib.get_energy(&vec![0.0; 11]);
    }

    #[test]
    fn crystalff_torsionanglecontribs_get_grad_accumulates_multiple_contributions() {
        let ff = torsion_forcefield();
        let pos = flattened_positions(&ff);
        let mut contrib = TorsionAngleContribs::new(&ff);
        let mut first = vec![0.0; 6];
        first[2] = 4.0;
        let mut second = vec![0.0; 6];
        second[0] = 1.5;
        second[1] = 0.5;
        contrib.add_contrib(0, 1, 2, 3, first, vec![1; 6]);
        contrib.add_contrib(0, 1, 2, 3, second, vec![1, -1, 1, 1, 1, 1]);

        let mut grad = vec![0.0; pos.len()];
        contrib.get_grad(&pos, &mut grad);

        assert!(grad.iter().all(|value| value.is_finite()));
        assert!(
            grad.iter().any(|value| value.abs() > EPS_TEST),
            "expected non-zero accumulated gradient"
        );
    }

    #[test]
    fn crystalff_torsionanglecontribs_get_grad_returns_early_for_degenerate_torsion() {
        let ff = degenerate_torsion_forcefield();
        let pos = flattened_positions(&ff);
        let mut contrib = TorsionAngleContribs::new(&ff);
        let mut force_constants = vec![0.0; 6];
        force_constants[2] = 4.0;
        contrib.add_contrib(0, 1, 2, 3, force_constants, vec![1; 6]);
        let mut grad = vec![3.0; pos.len()];

        contrib.get_grad(&pos, &mut grad);

        assert!(grad.iter().all(|value| *value == 3.0));
    }

    #[test]
    fn crystalff_torsionanglecontribs_get_grad_zero_force_constants_leave_gradient_unchanged() {
        let ff = torsion_forcefield();
        let pos = flattened_positions(&ff);
        let mut contrib = TorsionAngleContribs::new(&ff);
        contrib.add_contrib(0, 1, 2, 3, vec![0.0; 6], vec![1; 6]);
        let mut grad = vec![0.0; pos.len()];

        contrib.get_grad(&pos, &mut grad);

        assert!(grad.iter().all(|value| value.abs() < EPS_TEST));
    }

    #[test]
    #[should_panic(expected = "no owner")]
    fn crystalff_torsionanglecontribs_get_grad_panics_without_owner() {
        let ff = torsion_forcefield();
        let pos = flattened_positions(&ff);
        let mut contrib = TorsionAngleContribs::new(&ff);
        contrib.set_force_field(std::ptr::null());
        let mut grad = vec![0.0; pos.len()];

        contrib.get_grad(&pos, &mut grad);
    }

    #[test]
    #[should_panic(expected = "bad vector")]
    fn crystalff_torsionanglecontribs_get_grad_panics_for_empty_position_vector() {
        let ff = torsion_forcefield();
        let contrib = TorsionAngleContribs::new(&ff);
        let mut grad = vec![0.0; 12];

        contrib.get_grad(&[], &mut grad);
    }

    #[test]
    #[should_panic(expected = "bad vector")]
    fn crystalff_torsionanglecontribs_get_grad_panics_for_empty_gradient_vector() {
        let ff = torsion_forcefield();
        let pos = flattened_positions(&ff);
        let contrib = TorsionAngleContribs::new(&ff);
        let mut grad = vec![];

        contrib.get_grad(&pos, &mut grad);
    }

    #[test]
    #[should_panic(expected = "bad vector length for force-field positions")]
    fn crystalff_torsionanglecontribs_get_grad_panics_for_short_position_vector() {
        let ff = torsion_forcefield();
        let contrib = TorsionAngleContribs::new(&ff);
        let mut grad = vec![0.0; 12];

        contrib.get_grad(&vec![0.0; 11], &mut grad);
    }

    #[test]
    #[should_panic(expected = "bad gradient length for force-field positions")]
    fn crystalff_torsionanglecontribs_get_grad_panics_for_short_gradient_vector() {
        let ff = torsion_forcefield();
        let pos = flattened_positions(&ff);
        let contrib = TorsionAngleContribs::new(&ff);
        let mut grad = vec![0.0; 11];

        contrib.get_grad(&pos, &mut grad);
    }
}
