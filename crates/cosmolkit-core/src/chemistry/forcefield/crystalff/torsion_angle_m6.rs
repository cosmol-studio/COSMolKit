//! Source-backed RDKit CrystalFF single torsion-angle M6 contribution.

use crate::chemistry::forcefield::core::{ForceField, ForceFieldVec3, compute_dihedral_from_flat};
use crate::chemistry::forcefield::crystalff::calc_torsion_energy;

fn is_double_zero(value: f64) -> bool {
    // BEGIN RDKIT CPP HELPER ForceFields::MMFF::isDoubleZero (Params.h via TorsionAngleM6.cpp)
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
    // BEGIN RDKIT CPP FUNCTION ForceFields::MMFF::Utils::calcTorsionGrad (TorsionAngle.cpp:51-95 via TorsionAngleM6.cpp)
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
pub struct TorsionAngleContribM6 {
    owner: *const ForceField,
    at1_idx: usize,
    at2_idx: usize,
    at3_idx: usize,
    at4_idx: usize,
    force_constants: Vec<f64>,
    signs: Vec<i32>,
}

impl TorsionAngleContribM6 {
    #[must_use]
    pub fn new(
        owner: &ForceField,
        idx1: usize,
        idx2: usize,
        idx3: usize,
        idx4: usize,
        force_constants: Vec<f64>,
        signs: Vec<i32>,
    ) -> Self {
        // BEGIN RDKIT CPP CONSTRUCTOR ForceFields::CrystalFF::TorsionAngleContribM6::TorsionAngleContribM6 (TorsionAngleM6.h:49-52, TorsionAngleM6.cpp:46-64)
        // RDKit✔️❌: TorsionAngleContribM6::TorsionAngleContribM6(
        // RDKit✔️❌:     ForceFields::ForceField *owner, unsigned int idx1, unsigned int idx2,
        // RDKit✔️❌:     unsigned int idx3, unsigned int idx4, std::vector<double> V,
        // RDKit✔️❌:     std::vector<int> signs)
        // RDKit✔️❌:     : ForceFieldContrib(owner),
        // RDKit✔️❌:       d_at1Idx(idx1),
        // RDKit✔️❌:       d_at2Idx(idx2),
        // RDKit✔️❌:       d_at3Idx(idx3),
        // RDKit✔️❌:       d_at4Idx(idx4),
        // RDKit✔️❌:       d_V(std::move(V)),
        // RDKit✔️❌:       d_sign(std::move(signs)) {
        // RDKit✔️❌:   PRECONDITION(owner, "bad owner");
        // Rust references reproduce RDKit's non-null owner precondition.
        // RDKit✔️❌:   PRECONDITION((idx1 != idx2) && (idx1 != idx3) && (idx1 != idx4) &&
        // RDKit✔️❌:                    (idx2 != idx3) && (idx2 != idx4) && (idx3 != idx4),
        // RDKit✔️❌:                "degenerate points");
        assert!(
            idx1 != idx2
                && idx1 != idx3
                && idx1 != idx4
                && idx2 != idx3
                && idx2 != idx4
                && idx3 != idx4,
            "degenerate points"
        );
        // RDKit✔️❌:   URANGE_CHECK(idx1, owner->positions().size());
        // RDKit✔️❌:   URANGE_CHECK(idx2, owner->positions().size());
        // RDKit✔️❌:   URANGE_CHECK(idx3, owner->positions().size());
        // RDKit✔️❌:   URANGE_CHECK(idx4, owner->positions().size());
        let positions_len = owner.positions().len();
        assert!(idx1 < positions_len);
        assert!(idx2 < positions_len);
        assert!(idx3 < positions_len);
        assert!(idx4 < positions_len);
        // RDKit✔️❌: };
        Self {
            owner: owner as *const ForceField,
            at1_idx: idx1,
            at2_idx: idx2,
            at3_idx: idx3,
            at4_idx: idx4,
            force_constants,
            signs,
        }
    }

    #[must_use]
    pub const fn owner(&self) -> *const ForceField {
        self.owner
    }

    #[must_use]
    pub const fn atom_indices(&self) -> [usize; 4] {
        [self.at1_idx, self.at2_idx, self.at3_idx, self.at4_idx]
    }

    #[must_use]
    pub fn force_constants(&self) -> &[f64] {
        &self.force_constants
    }

    #[must_use]
    pub fn signs(&self) -> &[i32] {
        &self.signs
    }

    #[must_use]
    pub fn get_energy(&self, pos: &[f64]) -> f64 {
        // BEGIN RDKIT CPP METHOD ForceFields::CrystalFF::TorsionAngleContribM6::getEnergy (TorsionAngleM6.cpp:67-80)
        // RDKit✔️❌: double TorsionAngleContribM6::getEnergy(double *pos) const {
        // RDKit✔️❌:   PRECONDITION(dp_forceField, "no owner");
        assert!(!self.owner.is_null(), "no owner");
        // RDKit✔️❌:   PRECONDITION(pos, "bad vector");
        assert!(!pos.is_empty(), "bad vector");
        let owner = unsafe { self.owner.as_ref() }.expect("no owner");
        assert!(
            pos.len() >= 3 * owner.positions().len(),
            "bad vector length for force-field positions"
        );
        // RDKit✔️❌:
        // RDKit✔️❌:   RDGeom::Point3D iPoint(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
        // RDKit✔️❌:                          pos[3 * d_at1Idx + 2]);
        // RDKit✔️❌:   RDGeom::Point3D jPoint(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
        // RDKit✔️❌:                          pos[3 * d_at2Idx + 2]);
        // RDKit✔️❌:   RDGeom::Point3D kPoint(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
        // RDKit✔️❌:                          pos[3 * d_at3Idx + 2]);
        // RDKit✔️❌:   RDGeom::Point3D lPoint(pos[3 * d_at4Idx], pos[3 * d_at4Idx + 1],
        // RDKit✔️❌:                          pos[3 * d_at4Idx + 2]);
        // RDKit✔️❌:
        // RDKit✔️❌:   return calcTorsionEnergyM6(
        // RDKit✔️❌:       d_V, d_sign, Utils::calcTorsionCosPhi(iPoint, jPoint, kPoint, lPoint));
        let cos_phi = compute_dihedral_from_flat(
            pos,
            self.at1_idx,
            self.at2_idx,
            self.at3_idx,
            self.at4_idx,
            false,
        )
        .cos_phi;
        let energy = calc_torsion_energy(&self.force_constants, &self.signs, cos_phi);
        // RDKit✔️❌: }
        energy
    }

    pub fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
        // BEGIN RDKIT CPP METHOD ForceFields::CrystalFF::TorsionAngleContribM6::getGrad (TorsionAngleM6.cpp:82-142)
        // RDKit✔️❌: void TorsionAngleContribM6::getGrad(double *pos, double *grad) const {
        // RDKit✔️❌:   PRECONDITION(dp_forceField, "no owner");
        assert!(!self.owner.is_null(), "no owner");
        // RDKit✔️❌:   PRECONDITION(pos, "bad vector");
        assert!(!pos.is_empty(), "bad vector");
        // RDKit✔️❌:   PRECONDITION(grad, "bad vector");
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

        // RDKit✔️❌:   RDGeom::Point3D iPoint(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
        let i_point = ForceFieldVec3::new(
            pos[3 * self.at1_idx],
            pos[3 * self.at1_idx + 1],
            pos[3 * self.at1_idx + 2],
        );
        // RDKit✔️❌:   RDGeom::Point3D jPoint(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
        let j_point = ForceFieldVec3::new(
            pos[3 * self.at2_idx],
            pos[3 * self.at2_idx + 1],
            pos[3 * self.at2_idx + 2],
        );
        // RDKit✔️❌:   RDGeom::Point3D kPoint(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
        let k_point = ForceFieldVec3::new(
            pos[3 * self.at3_idx],
            pos[3 * self.at3_idx + 1],
            pos[3 * self.at3_idx + 2],
        );
        // RDKit✔️❌:   RDGeom::Point3D lPoint(pos[3 * d_at4Idx], pos[3 * d_at4Idx + 1],
        let l_point = ForceFieldVec3::new(
            pos[3 * self.at4_idx],
            pos[3 * self.at4_idx + 1],
            pos[3 * self.at4_idx + 2],
        );

        // RDKit✔️❌:   double *g[4] = {&(grad[3 * d_at1Idx]), &(grad[3 * d_at2Idx]),
        // RDKit✔️❌:                   &(grad[3 * d_at3Idx]), &(grad[3 * d_at4Idx])};
        // Gradient slices are addressed directly in `calc_torsion_grad`.

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
        // RDKit✔️❌:   if (isDoubleZero(d[0]) || isDoubleZero(d[1])) {
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
        // RDKit✔️❌:   clipToOne(cosPhi);
        cos_phi = cos_phi.clamp(-1.0, 1.0);
        // RDKit✔️❌:   double sinPhiSq = 1.0 - cosPhi * cosPhi;
        let sin_phi_sq = 1.0 - cos_phi * cos_phi;
        // RDKit✔️❌:   double sinPhi = ((sinPhiSq > 0.0) ? sqrt(sinPhiSq) : 0.0);
        let sin_phi = if sin_phi_sq > 0.0 {
            sin_phi_sq.sqrt()
        } else {
            0.0
        };
        // RDKit✔️❌:   double cosPhi2 = cosPhi * cosPhi;
        let cos_phi2 = cos_phi * cos_phi;
        // RDKit✔️❌:   double cosPhi3 = cosPhi * cosPhi2;
        let cos_phi3 = cos_phi * cos_phi2;
        // RDKit✔️❌:   double cosPhi4 = cosPhi * cosPhi3;
        let cos_phi4 = cos_phi * cos_phi3;
        // RDKit✔️❌:   double cosPhi5 = cosPhi * cosPhi4;
        let cos_phi5 = cos_phi * cos_phi4;
        // RDKit✔️❌:   // dE/dPhi is independent of cartesians:
        // RDKit✔️❌:   double dE_dPhi =
        let d_e_d_phi = -self.force_constants[0] * f64::from(self.signs[0]) * sin_phi
            - 2.0 * self.force_constants[1] * f64::from(self.signs[1]) * (2.0 * cos_phi * sin_phi)
            - 3.0
                * self.force_constants[2]
                * f64::from(self.signs[2])
                * (4.0 * cos_phi2 * sin_phi - sin_phi)
            - 4.0
                * self.force_constants[3]
                * f64::from(self.signs[3])
                * (8.0 * cos_phi3 * sin_phi - 4.0 * cos_phi * sin_phi)
            - 5.0
                * self.force_constants[4]
                * f64::from(self.signs[4])
                * (16.0 * cos_phi4 * sin_phi - 12.0 * cos_phi2 * sin_phi + sin_phi)
            - 6.0
                * self.force_constants[4]
                * f64::from(self.signs[4])
                * (32.0 * cos_phi5 * sin_phi - 32.0 * cos_phi3 * sin_phi + 6.0 * sin_phi);

        // RDKit✔️❌:   // FIX: use a tolerance here
        // RDKit✔️❌:   // this is hacky, but it's per the
        // RDKit✔️❌:   // recommendation from Niketic and Rasmussen:
        // RDKit✔️❌:   double sinTerm =
        // RDKit✔️❌:       -dE_dPhi * (isDoubleZero(sinPhi) ? (1.0 / cosPhi) : (1.0 / sinPhi));
        let sin_term = -d_e_d_phi
            * if is_double_zero(sin_phi) {
                1.0 / cos_phi
            } else {
                1.0 / sin_phi
            };

        // RDKit✔️❌:   Utils::calcTorsionGrad(r, t, d, g, sinTerm, cosPhi);
        calc_torsion_grad(
            r,
            t,
            d,
            grad,
            [self.at1_idx, self.at2_idx, self.at3_idx, self.at4_idx],
            sin_term,
            cos_phi,
        );
        // RDKit✔️❌: }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::forcefield::core::{
        ForceField, ForceFieldContrib, ForceFieldVec3, compute_dihedral_from_flat,
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

    #[test]
    fn crystalff_torsionanglecontribm6_constructor_initializes_fields() {
        let ff = torsion_forcefield();
        let force_constants = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let signs = vec![1, -1, 1, -1, 1, -1];

        let contrib =
            TorsionAngleContribM6::new(&ff, 0, 1, 2, 3, force_constants.clone(), signs.clone());

        assert_eq!(contrib.owner(), &ff as *const ForceField);
        assert_eq!(contrib.atom_indices(), [0, 1, 2, 3]);
        assert_eq!(contrib.force_constants(), force_constants.as_slice());
        assert_eq!(contrib.signs(), signs.as_slice());
    }

    #[test]
    #[should_panic(expected = "degenerate points")]
    fn crystalff_torsionanglecontribm6_constructor_panics_for_degenerate_indices() {
        let ff = torsion_forcefield();

        let _ = TorsionAngleContribM6::new(&ff, 0, 0, 2, 3, vec![0.0; 6], vec![1; 6]);
    }

    #[test]
    #[should_panic]
    fn crystalff_torsionanglecontribm6_constructor_panics_for_out_of_range_first_index() {
        let ff = torsion_forcefield();

        let _ = TorsionAngleContribM6::new(&ff, 4, 1, 2, 3, vec![0.0; 6], vec![1; 6]);
    }

    #[test]
    #[should_panic]
    fn crystalff_torsionanglecontribm6_constructor_panics_for_out_of_range_second_index() {
        let ff = torsion_forcefield();

        let _ = TorsionAngleContribM6::new(&ff, 0, 4, 2, 3, vec![0.0; 6], vec![1; 6]);
    }

    #[test]
    #[should_panic]
    fn crystalff_torsionanglecontribm6_constructor_panics_for_out_of_range_third_index() {
        let ff = torsion_forcefield();

        let _ = TorsionAngleContribM6::new(&ff, 0, 1, 4, 3, vec![0.0; 6], vec![1; 6]);
    }

    #[test]
    #[should_panic]
    fn crystalff_torsionanglecontribm6_constructor_panics_for_out_of_range_fourth_index() {
        let ff = torsion_forcefield();

        let _ = TorsionAngleContribM6::new(&ff, 0, 1, 2, 4, vec![0.0; 6], vec![1; 6]);
    }

    #[test]
    fn crystalff_torsionanglecontribm6_get_energy_matches_sp3_sp3_source_geometry() {
        let ff = torsion_forcefield();
        let pos = flattened_positions(&ff);
        let force_constants = vec![0.0, 0.0, 4.0, 0.0, 0.0, 0.0];
        let signs = vec![1, 1, 1, 1, 1, 1];
        let contrib = TorsionAngleContribM6::new(&ff, 0, 1, 2, 3, force_constants, signs);

        let energy = contrib.get_energy(&pos);
        let cos_phi = compute_dihedral_from_flat(&pos, 0, 1, 2, 3, false).cos_phi;

        assert_eq!(cos_phi, 0.0);
        assert_eq!(energy, 4.0);
    }

    #[test]
    fn crystalff_torsionanglecontribm6_get_energy_matches_m6_closed_form_for_geometry() {
        let ff = torsion_forcefield();
        let pos = flattened_positions(&ff);
        let force_constants = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let signs = vec![1, -1, 1, -1, 1, -1];
        let contrib =
            TorsionAngleContribM6::new(&ff, 0, 1, 2, 3, force_constants.clone(), signs.clone());

        let energy = contrib.get_energy(&pos);
        let cos_phi = compute_dihedral_from_flat(&pos, 0, 1, 2, 3, false).cos_phi;
        let expected = calc_torsion_energy(&force_constants, &signs, cos_phi);

        assert_eq!(energy, expected);
    }

    #[test]
    #[should_panic(expected = "no owner")]
    fn crystalff_torsionanglecontribm6_get_energy_panics_without_owner() {
        let ff = torsion_forcefield();
        let pos = flattened_positions(&ff);
        let mut contrib = TorsionAngleContribM6::new(&ff, 0, 1, 2, 3, vec![0.0; 6], vec![1; 6]);
        contrib.owner = std::ptr::null();

        let _ = contrib.get_energy(&pos);
    }

    #[test]
    #[should_panic(expected = "bad vector")]
    fn crystalff_torsionanglecontribm6_get_energy_panics_for_empty_position_vector() {
        let ff = torsion_forcefield();
        let contrib = TorsionAngleContribM6::new(&ff, 0, 1, 2, 3, vec![0.0; 6], vec![1; 6]);

        let _ = contrib.get_energy(&[]);
    }

    #[test]
    #[should_panic(expected = "bad vector length for force-field positions")]
    fn crystalff_torsionanglecontribm6_get_energy_panics_for_short_position_vector() {
        let ff = torsion_forcefield();
        let contrib = TorsionAngleContribM6::new(&ff, 0, 1, 2, 3, vec![0.0; 6], vec![1; 6]);

        let _ = contrib.get_energy(&vec![0.0; 3 * ff.positions().len() - 1]);
    }

    #[test]
    fn crystalff_torsionanglecontribm6_get_grad_zero_force_constants_leave_gradient_unchanged() {
        let ff = torsion_forcefield();
        let pos = flattened_positions(&ff);
        let mut grad = vec![1.25; 3 * ff.positions().len()];
        let before = grad.clone();
        let contrib = TorsionAngleContribM6::new(&ff, 0, 1, 2, 3, vec![0.0; 6], vec![1; 6]);

        contrib.get_grad(&pos, &mut grad);

        assert_eq!(grad, before);
    }

    #[test]
    fn crystalff_torsionanglecontribm6_get_grad_matches_single_contribs_gradient() {
        let ff = torsion_forcefield();
        let pos = flattened_positions(&ff);
        let force_constants = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let signs = vec![1, -1, 1, -1, 1, -1];
        let contrib =
            TorsionAngleContribM6::new(&ff, 0, 1, 2, 3, force_constants.clone(), signs.clone());
        let mut grad = vec![0.0; pos.len()];
        let mut contribs = crate::chemistry::forcefield::crystalff::TorsionAngleContribs::new(&ff);
        contribs.add_contrib(0, 1, 2, 3, force_constants, signs);
        let mut expected_grad = vec![0.0; pos.len()];

        contrib.get_grad(&pos, &mut grad);
        contribs.get_grad(&pos, &mut expected_grad);

        for (actual, expected) in grad.iter().zip(expected_grad.iter()) {
            assert_close(*actual, *expected);
        }
    }

    #[test]
    fn crystalff_torsionanglecontribm6_get_grad_returns_early_for_degenerate_torsion() {
        let ff = degenerate_torsion_forcefield();
        let pos = flattened_positions(&ff);
        let mut grad = vec![1.0; 3 * ff.positions().len()];
        let before = grad.clone();
        let contrib = TorsionAngleContribM6::new(
            &ff,
            0,
            1,
            2,
            3,
            vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            vec![1, -1, 1, -1, 1, -1],
        );

        contrib.get_grad(&pos, &mut grad);

        assert_eq!(grad, before);
    }

    #[test]
    #[should_panic(expected = "no owner")]
    fn crystalff_torsionanglecontribm6_get_grad_panics_without_owner() {
        let ff = torsion_forcefield();
        let pos = flattened_positions(&ff);
        let mut grad = vec![0.0; 3 * ff.positions().len()];
        let mut contrib = TorsionAngleContribM6::new(&ff, 0, 1, 2, 3, vec![0.0; 6], vec![1; 6]);
        contrib.owner = std::ptr::null();

        contrib.get_grad(&pos, &mut grad);
    }

    #[test]
    #[should_panic(expected = "bad vector")]
    fn crystalff_torsionanglecontribm6_get_grad_panics_for_empty_position_vector() {
        let ff = torsion_forcefield();
        let mut grad = vec![0.0; 3 * ff.positions().len()];
        let contrib = TorsionAngleContribM6::new(&ff, 0, 1, 2, 3, vec![0.0; 6], vec![1; 6]);

        contrib.get_grad(&[], &mut grad);
    }

    #[test]
    #[should_panic(expected = "bad vector")]
    fn crystalff_torsionanglecontribm6_get_grad_panics_for_empty_gradient_vector() {
        let ff = torsion_forcefield();
        let pos = flattened_positions(&ff);
        let mut grad = Vec::new();
        let contrib = TorsionAngleContribM6::new(&ff, 0, 1, 2, 3, vec![0.0; 6], vec![1; 6]);

        contrib.get_grad(&pos, &mut grad);
    }

    #[test]
    #[should_panic(expected = "bad vector length for force-field positions")]
    fn crystalff_torsionanglecontribm6_get_grad_panics_for_short_position_vector() {
        let ff = torsion_forcefield();
        let mut grad = vec![0.0; 3 * ff.positions().len()];
        let contrib = TorsionAngleContribM6::new(&ff, 0, 1, 2, 3, vec![0.0; 6], vec![1; 6]);

        contrib.get_grad(&vec![0.0; 3 * ff.positions().len() - 1], &mut grad);
    }

    #[test]
    #[should_panic(expected = "bad gradient length for force-field positions")]
    fn crystalff_torsionanglecontribm6_get_grad_panics_for_short_gradient_vector() {
        let ff = torsion_forcefield();
        let pos = flattened_positions(&ff);
        let mut grad = vec![0.0; 3 * ff.positions().len() - 1];
        let contrib = TorsionAngleContribM6::new(&ff, 0, 1, 2, 3, vec![0.0; 6], vec![1; 6]);

        contrib.get_grad(&pos, &mut grad);
    }
}
