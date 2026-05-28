//! Source-backed RDKit MMFF torsion-angle contribution.

use crate::chemistry::forcefield::core::{
    ForceField, ForceFieldContrib, ForceFieldVec3, compute_dihedral_from_flat,
};

use super::params::MmffTor;

fn point_from_pos(pos: &[f64], atom_idx: usize) -> ForceFieldVec3 {
    ForceFieldVec3::new(
        pos[3 * atom_idx],
        pos[3 * atom_idx + 1],
        pos[3 * atom_idx + 2],
    )
}

fn is_double_zero(value: f64) -> bool {
    // BEGIN RDKIT CPP HELPER ForceFields::isDoubleZero (Params.h via TorsionAngle.cpp)
    // RDKit✔️✔️: inline bool isDoubleZero(const double x) {
    // RDKit✔️✔️:   return ((x < 1.0e-10) && (x > -1.0e-10));
    // RDKit✔️✔️: }
    value < 1.0e-10 && value > -1.0e-10
}

fn clip_to_one(value: &mut f64) {
    // BEGIN RDKIT CPP HELPER ForceFields::clipToOne (Params.h via TorsionAngle.cpp)
    // RDKit✔️✔️: inline void clipToOne(double &x) { x = std::clamp(x, -1.0, 1.0); }
    *value = value.clamp(-1.0, 1.0);
}

fn calc_torsion_cos_phi(
    i_point: ForceFieldVec3,
    j_point: ForceFieldVec3,
    k_point: ForceFieldVec3,
    l_point: ForceFieldVec3,
) -> f64 {
    // BEGIN RDKIT CPP FUNCTION ForceFields::MMFF::Utils::calcTorsionCosPhi (TorsionAngle.cpp:17-35)
    // RDKit✔️✔️: double calcTorsionCosPhi(const RDGeom::Point3D &iPoint,
    // RDKit✔️✔️:                          const RDGeom::Point3D &jPoint,
    // RDKit✔️✔️:                          const RDGeom::Point3D &kPoint,
    // RDKit✔️✔️:                          const RDGeom::Point3D &lPoint) {
    // RDKit✔️✔️:   RDGeom::Point3D r1 = iPoint - jPoint;
    let r1 = i_point - j_point;
    // RDKit✔️✔️:   RDGeom::Point3D r2 = kPoint - jPoint;
    let r2 = k_point - j_point;
    // RDKit✔️✔️:   RDGeom::Point3D r3 = jPoint - kPoint;
    let r3 = j_point - k_point;
    // RDKit✔️✔️:   RDGeom::Point3D r4 = lPoint - kPoint;
    let r4 = l_point - k_point;
    // RDKit✔️✔️:   RDGeom::Point3D t1 = r1.crossProduct(r2);
    let t1 = r1.cross_product(r2);
    // RDKit✔️✔️:   RDGeom::Point3D t2 = r3.crossProduct(r4);
    let t2 = r3.cross_product(r4);
    // RDKit✔️✔️:   auto t1_len = t1.length();
    let t1_len = t1.length();
    // RDKit✔️✔️:   auto t2_len = t2.length();
    let t2_len = t2.length();
    // RDKit✔️✔️:   if (isDoubleZero(t1_len) || isDoubleZero(t2_len)) {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   }
    if is_double_zero(t1_len) || is_double_zero(t2_len) {
        return 0.0;
    }
    // RDKit✔️✔️:   double cosPhi = t1.dotProduct(t2) / (t1_len * t2_len);
    let mut cos_phi = t1.dot_product(t2) / (t1_len * t2_len);
    // RDKit✔️✔️:   clipToOne(cosPhi);
    clip_to_one(&mut cos_phi);
    // RDKit✔️✔️:   return cosPhi;
    // RDKit✔️✔️: }
    cos_phi
}

fn calc_torsion_energy(v1: f64, v2: f64, v3: f64, cos_phi: f64) -> f64 {
    // BEGIN RDKIT CPP FUNCTION ForceFields::MMFF::Utils::calcTorsionEnergy (TorsionAngle.cpp:43-49)
    // RDKit✔️✔️: double calcTorsionEnergy(const double V1, const double V2, const double V3,
    // RDKit✔️✔️:                          const double cosPhi) {
    // RDKit✔️✔️:   double cos2Phi = 2.0 * cosPhi * cosPhi - 1.0;
    let cos2_phi = 2.0 * cos_phi * cos_phi - 1.0;
    // RDKit✔️✔️:   double cos3Phi = cosPhi * (2.0 * cos2Phi - 1.0);
    let cos3_phi = cos_phi * (2.0 * cos2_phi - 1.0);

    // RDKit✔️✔️:   return (0.5 *
    // RDKit✔️✔️:           (V1 * (1.0 + cosPhi) + V2 * (1.0 - cos2Phi) + V3 * (1.0 + cos3Phi)));
    // RDKit✔️✔️: }
    0.5 * (v1 * (1.0 + cos_phi) + v2 * (1.0 - cos2_phi) + v3 * (1.0 + cos3_phi))
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
    // BEGIN RDKIT CPP FUNCTION ForceFields::MMFF::Utils::calcTorsionGrad (TorsionAngle.cpp:51-95)
    // RDKit✔️✔️: void calcTorsionGrad(RDGeom::Point3D *r, RDGeom::Point3D *t, double *d,
    // RDKit✔️✔️:                      double **g, double &sinTerm, double &cosPhi) {
    // RDKit✔️✔️:   // -------
    // RDKit✔️✔️:   // dTheta/dx is trickier:
    // RDKit✔️✔️:   double dCos_dT[6] = {1.0 / d[0] * (t[1].x - cosPhi * t[0].x),
    // RDKit✔️✔️:                        1.0 / d[0] * (t[1].y - cosPhi * t[0].y),
    // RDKit✔️✔️:                        1.0 / d[0] * (t[1].z - cosPhi * t[0].z),
    // RDKit✔️✔️:                        1.0 / d[1] * (t[0].x - cosPhi * t[1].x),
    // RDKit✔️✔️:                        1.0 / d[1] * (t[0].y - cosPhi * t[1].y),
    // RDKit✔️✔️:                        1.0 / d[1] * (t[0].z - cosPhi * t[1].z)};
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

    // RDKit✔️✔️:   g[0][0] += sinTerm * (dCos_dT[2] * r[1].y - dCos_dT[1] * r[1].z);
    grad[atom1_offset] += sin_term * (d_cos_dt[2] * r[1].y - d_cos_dt[1] * r[1].z);
    // RDKit✔️✔️:   g[0][1] += sinTerm * (dCos_dT[0] * r[1].z - dCos_dT[2] * r[1].x);
    grad[atom1_offset + 1] += sin_term * (d_cos_dt[0] * r[1].z - d_cos_dt[2] * r[1].x);
    // RDKit✔️✔️:   g[0][2] += sinTerm * (dCos_dT[1] * r[1].x - dCos_dT[0] * r[1].y);
    grad[atom1_offset + 2] += sin_term * (d_cos_dt[1] * r[1].x - d_cos_dt[0] * r[1].y);

    // RDKit✔️✔️:   g[1][0] += sinTerm *
    // RDKit✔️✔️:              (dCos_dT[1] * (r[1].z - r[0].z) + dCos_dT[2] * (r[0].y - r[1].y) +
    // RDKit✔️✔️:               dCos_dT[4] * (-r[3].z) + dCos_dT[5] * (r[3].y));
    grad[atom2_offset] += sin_term
        * ((d_cos_dt[1] * (r[1].z - r[0].z))
            + (d_cos_dt[2] * (r[0].y - r[1].y))
            + (d_cos_dt[4] * (-r[3].z))
            + (d_cos_dt[5] * r[3].y));
    // RDKit✔️✔️:   g[1][1] += sinTerm *
    // RDKit✔️✔️:              (dCos_dT[0] * (r[0].z - r[1].z) + dCos_dT[2] * (r[1].x - r[0].x) +
    // RDKit✔️✔️:               dCos_dT[3] * (r[3].z) + dCos_dT[5] * (-r[3].x));
    grad[atom2_offset + 1] += sin_term
        * ((d_cos_dt[0] * (r[0].z - r[1].z))
            + (d_cos_dt[2] * (r[1].x - r[0].x))
            + (d_cos_dt[3] * r[3].z)
            + (d_cos_dt[5] * (-r[3].x)));
    // RDKit✔️✔️:   g[1][2] += sinTerm *
    // RDKit✔️✔️:              (dCos_dT[0] * (r[1].y - r[0].y) + dCos_dT[1] * (r[0].x - r[1].x) +
    // RDKit✔️✔️:               dCos_dT[3] * (-r[3].y) + dCos_dT[4] * (r[3].x));
    grad[atom2_offset + 2] += sin_term
        * ((d_cos_dt[0] * (r[1].y - r[0].y))
            + (d_cos_dt[1] * (r[0].x - r[1].x))
            + (d_cos_dt[3] * (-r[3].y))
            + (d_cos_dt[4] * r[3].x));

    // RDKit✔️✔️:   g[2][0] += sinTerm *
    // RDKit✔️✔️:              (dCos_dT[1] * (r[0].z) + dCos_dT[2] * (-r[0].y) +
    // RDKit✔️✔️:               dCos_dT[4] * (r[3].z - r[2].z) + dCos_dT[5] * (r[2].y - r[3].y));
    grad[atom3_offset] += sin_term
        * ((d_cos_dt[1] * r[0].z)
            + (d_cos_dt[2] * (-r[0].y))
            + (d_cos_dt[4] * (r[3].z - r[2].z))
            + (d_cos_dt[5] * (r[2].y - r[3].y)));
    // RDKit✔️✔️:   g[2][1] += sinTerm *
    // RDKit✔️✔️:              (dCos_dT[0] * (-r[0].z) + dCos_dT[2] * (r[0].x) +
    // RDKit✔️✔️:               dCos_dT[3] * (r[2].z - r[3].z) + dCos_dT[5] * (r[3].x - r[2].x));
    grad[atom3_offset + 1] += sin_term
        * ((d_cos_dt[0] * (-r[0].z))
            + (d_cos_dt[2] * r[0].x)
            + (d_cos_dt[3] * (r[2].z - r[3].z))
            + (d_cos_dt[5] * (r[3].x - r[2].x)));
    // RDKit✔️✔️:   g[2][2] += sinTerm *
    // RDKit✔️✔️:              (dCos_dT[0] * (r[0].y) + dCos_dT[1] * (-r[0].x) +
    // RDKit✔️✔️:               dCos_dT[3] * (r[3].y - r[2].y) + dCos_dT[4] * (r[2].x - r[3].x));
    grad[atom3_offset + 2] += sin_term
        * ((d_cos_dt[0] * r[0].y)
            + (d_cos_dt[1] * (-r[0].x))
            + (d_cos_dt[3] * (r[3].y - r[2].y))
            + (d_cos_dt[4] * (r[2].x - r[3].x)));

    // RDKit✔️✔️:   g[3][0] += sinTerm * (dCos_dT[4] * r[2].z - dCos_dT[5] * r[2].y);
    grad[atom4_offset] += sin_term * (d_cos_dt[4] * r[2].z - d_cos_dt[5] * r[2].y);
    // RDKit✔️✔️:   g[3][1] += sinTerm * (dCos_dT[5] * r[2].x - dCos_dT[3] * r[2].z);
    grad[atom4_offset + 1] += sin_term * (d_cos_dt[5] * r[2].x - d_cos_dt[3] * r[2].z);
    // RDKit✔️✔️:   g[3][2] += sinTerm * (dCos_dT[3] * r[2].y - dCos_dT[4] * r[2].x);
    grad[atom4_offset + 2] += sin_term * (d_cos_dt[3] * r[2].y - d_cos_dt[4] * r[2].x);
    // RDKit✔️✔️: }
}

#[derive(Clone, Debug)]
pub struct TorsionAngleContrib {
    owner: *const ForceField,
    atom1_indices: Vec<usize>,
    atom2_indices: Vec<usize>,
    atom3_indices: Vec<usize>,
    atom4_indices: Vec<usize>,
    v1: Vec<f64>,
    v2: Vec<f64>,
    v3: Vec<f64>,
}

impl TorsionAngleContrib {
    #[must_use]
    pub fn new(owner: &ForceField) -> Self {
        // BEGIN RDKIT CPP CONSTRUCTOR ForceFields::MMFF::TorsionAngleContrib::TorsionAngleContrib (TorsionAngle.cpp:100-103)
        // RDKit✔️✔️: TorsionAngleContrib::TorsionAngleContrib(ForceField *owner) {
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
            v1: Vec::new(),
            v2: Vec::new(),
            v3: Vec::new(),
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
            && self.v1.is_empty()
            && self.v2.is_empty()
            && self.v3.is_empty()
    }

    pub fn add_term(
        &mut self,
        idx1: usize,
        idx2: usize,
        idx3: usize,
        idx4: usize,
        mmff_tor_params: &MmffTor,
    ) {
        // BEGIN RDKIT CPP METHOD ForceFields::MMFF::TorsionAngleContrib::addTerm (TorsionAngle.cpp:105-123)
        // RDKit✔️✔️: void TorsionAngleContrib::addTerm(
        // RDKit✔️✔️:     unsigned int idx1, unsigned int idx2, unsigned int idx3, unsigned int idx4,
        // RDKit✔️✔️:     const ForceFields::MMFF::MMFFTor *mmffTorParams) {
        // RDKit✔️✔️:   PRECONDITION((idx1 != idx2) && (idx1 != idx3) && (idx1 != idx4) &&
        // RDKit✔️✔️:                    (idx2 != idx3) && (idx2 != idx4) && (idx3 != idx4),
        // RDKit✔️✔️:                "degenerate points");
        assert!(
            idx1 != idx2
                && idx1 != idx3
                && idx1 != idx4
                && idx2 != idx3
                && idx2 != idx4
                && idx3 != idx4,
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

        // RDKit✔️✔️:   d_at1Idx.push_back(idx1);
        // RDKit✔️✔️:   d_at2Idx.push_back(idx2);
        // RDKit✔️✔️:   d_at3Idx.push_back(idx3);
        // RDKit✔️✔️:   d_at4Idx.push_back(idx4);
        self.atom1_indices.push(idx1);
        self.atom2_indices.push(idx2);
        self.atom3_indices.push(idx3);
        self.atom4_indices.push(idx4);
        // RDKit✔️✔️:   d_V1.push_back(mmffTorParams->V1);
        // RDKit✔️✔️:   d_V2.push_back(mmffTorParams->V2);
        // RDKit✔️✔️:   d_V3.push_back(mmffTorParams->V3);
        self.v1.push(mmff_tor_params.v1);
        self.v2.push(mmff_tor_params.v2);
        self.v3.push(mmff_tor_params.v3);
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
    pub fn v1(&self) -> &[f64] {
        &self.v1
    }

    #[must_use]
    pub fn v2(&self) -> &[f64] {
        &self.v2
    }

    #[must_use]
    pub fn v3(&self) -> &[f64] {
        &self.v3
    }

    #[must_use]
    pub fn get_energy(&self, pos: &[f64]) -> f64 {
        // BEGIN RDKIT CPP METHOD ForceFields::MMFF::TorsionAngleContrib::getEnergy (TorsionAngle.cpp:125-150)
        // RDKit✔️✔️: double TorsionAngleContrib::getEnergy(double *pos) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        let _ = self.force_field();
        // Rust slices reproduce RDKit's non-null pos precondition.

        // RDKit✔️✔️:   const int numTorsions = d_at1Idx.size();
        // RDKit✔️✔️:   double totalEnergy = 0.0;
        let num_torsions = self.atom1_indices.len();
        let mut total_energy = 0.0;
        // RDKit✔️✔️:   for (int i = 0; i < numTorsions; ++i) {
        for i in 0..num_torsions {
            // RDKit✔️✔️:     const int16_t at1Idx = d_at1Idx[i];
            // RDKit✔️✔️:     const int16_t at2Idx = d_at2Idx[i];
            // RDKit✔️✔️:     const int16_t at3Idx = d_at3Idx[i];
            // RDKit✔️✔️:     const int16_t at4Idx = d_at4Idx[i];
            let atom1_idx = self.atom1_indices[i];
            let atom2_idx = self.atom2_indices[i];
            let atom3_idx = self.atom3_indices[i];
            let atom4_idx = self.atom4_indices[i];

            // RDKit✔️✔️:     RDGeom::Point3D iPoint(pos[3 * at1Idx], pos[3 * at1Idx + 1],
            // RDKit✔️✔️:                            pos[3 * at1Idx + 2]);
            let i_point = point_from_pos(pos, atom1_idx);
            // RDKit✔️✔️:     RDGeom::Point3D jPoint(pos[3 * at2Idx], pos[3 * at2Idx + 1],
            // RDKit✔️✔️:                            pos[3 * at2Idx + 2]);
            let j_point = point_from_pos(pos, atom2_idx);
            // RDKit✔️✔️:     RDGeom::Point3D kPoint(pos[3 * at3Idx], pos[3 * at3Idx + 1],
            // RDKit✔️✔️:                            pos[3 * at3Idx + 2]);
            let k_point = point_from_pos(pos, atom3_idx);
            // RDKit✔️✔️:     RDGeom::Point3D lPoint(pos[3 * at4Idx], pos[3 * at4Idx + 1],
            // RDKit✔️✔️:                            pos[3 * at4Idx + 2]);
            let l_point = point_from_pos(pos, atom4_idx);

            // RDKit✔️✔️:     totalEnergy += Utils::calcTorsionEnergy(
            // RDKit✔️✔️:         d_V1[i], d_V2[i], d_V3[i],
            // RDKit✔️✔️:         Utils::calcTorsionCosPhi(iPoint, jPoint, kPoint, lPoint));
            total_energy += calc_torsion_energy(
                self.v1[i],
                self.v2[i],
                self.v3[i],
                calc_torsion_cos_phi(i_point, j_point, k_point, l_point),
            );
            // RDKit✔️✔️:   }
        }
        // RDKit✔️✔️:   return totalEnergy;
        // RDKit✔️✔️: }
        total_energy
    }

    pub fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
        // BEGIN RDKIT CPP METHOD ForceFields::MMFF::TorsionAngleContrib::getGrad (TorsionAngle.cpp:152-186)
        // RDKit✔️✔️: void TorsionAngleContrib::getGrad(double *pos, double *grad) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:   PRECONDITION(grad, "bad vector");
        let _ = self.force_field();
        // Rust slices reproduce RDKit's non-null pos and grad preconditions.
        // RDKit✔️✔️:   double d[2];
        // RDKit✔️✔️:   double cosPhi;

        // RDKit✔️✔️:   const int numTorsions = d_at1Idx.size();
        let num_torsions = self.atom1_indices.len();
        // RDKit✔️✔️:   for (int i = 0; i < numTorsions; ++i) {
        for i in 0..num_torsions {
            // RDKit✔️✔️:     const int16_t at1Idx = d_at1Idx[i];
            // RDKit✔️✔️:     const int16_t at2Idx = d_at2Idx[i];
            // RDKit✔️✔️:     const int16_t at3Idx = d_at3Idx[i];
            // RDKit✔️✔️:     const int16_t at4Idx = d_at4Idx[i];
            let atom1_idx = self.atom1_indices[i];
            let atom2_idx = self.atom2_indices[i];
            let atom3_idx = self.atom3_indices[i];
            let atom4_idx = self.atom4_indices[i];

            // RDKit✔️✔️:     double *g[4] = {&(grad[3 * at1Idx]), &(grad[3 * at2Idx]),
            // RDKit✔️✔️:                     &(grad[3 * at3Idx]), &(grad[3 * at4Idx])};
            // Rust passes the atom indices directly to calc_torsion_grad and writes to
            // the same flat gradient offsets.

            // RDKit✔️✔️:     RDGeom::Point3D r[4];
            // RDKit✔️✔️:     RDGeom::Point3D t[2];

            // RDKit✔️✔️:     RDKit::ForceFieldsHelper::computeDihedral(
            // RDKit✔️✔️:         pos, at1Idx, at2Idx, at3Idx, at4Idx, nullptr, &cosPhi, r, t, d);
            let dihedral_output =
                compute_dihedral_from_flat(pos, atom1_idx, atom2_idx, atom3_idx, atom4_idx, false);
            let cos_phi = dihedral_output.cos_phi;
            let r = dihedral_output.r;
            let t = dihedral_output.t;
            let d = dihedral_output.d;
            // RDKit✔️✔️:     double sinPhiSq = 1.0 - cosPhi * cosPhi;
            let sin_phi_sq = 1.0 - cos_phi * cos_phi;
            // RDKit✔️✔️:     double sinPhi = ((sinPhiSq > 0.0) ? sqrt(sinPhiSq) : 0.0);
            let sin_phi = if sin_phi_sq > 0.0 {
                sin_phi_sq.sqrt()
            } else {
                0.0
            };
            // RDKit✔️✔️:     double sin2Phi = 2.0 * sinPhi * cosPhi;
            let sin2_phi = 2.0 * sin_phi * cos_phi;
            // RDKit✔️✔️:     double sin3Phi = 3.0 * sinPhi - 4.0 * sinPhi * sinPhiSq;
            let sin3_phi = 3.0 * sin_phi - 4.0 * sin_phi * sin_phi_sq;
            // RDKit✔️✔️:     // dE/dPhi is independent of cartesians:
            // RDKit✔️✔️:     double dE_dPhi = 0.5 * (-(d_V1[i]) * sinPhi + 2.0 * d_V2[i] * sin2Phi -
            // RDKit✔️✔️:                             3.0 * d_V3[i] * sin3Phi);
            let de_dphi = 0.5
                * (-(self.v1[i]) * sin_phi + 2.0 * self.v2[i] * sin2_phi
                    - 3.0 * self.v3[i] * sin3_phi);
            // RDKit✔️✔️:     // FIX: use a tolerance here
            // RDKit✔️✔️:     // this is hacky, but it's per the
            // RDKit✔️✔️:     // recommendation from Niketic and Rasmussen:
            // RDKit✔️✔️:     double sinTerm =
            // RDKit✔️✔️:         -dE_dPhi * (isDoubleZero(sinPhi) ? (1.0 / cosPhi) : (1.0 / sinPhi));
            let sin_term = -de_dphi
                * if is_double_zero(sin_phi) {
                    1.0 / cos_phi
                } else {
                    1.0 / sin_phi
                };

            // RDKit✔️✔️:     Utils::calcTorsionGrad(r, t, d, g, sinTerm, cosPhi);
            calc_torsion_grad(
                r,
                t,
                d,
                grad,
                [atom1_idx, atom2_idx, atom3_idx, atom4_idx],
                sin_term,
                cos_phi,
            );
            // RDKit✔️✔️:   }
        }
        // RDKit✔️✔️: }
    }

    fn force_field(&self) -> &ForceField {
        assert!(!self.owner.is_null(), "no owner");
        // SAFETY: MMFF contribs follow the same owner-pointer model as
        // ForceFieldContrib objects in core.rs. Constructors store a live
        // ForceField pointer before term insertion or evaluation.
        unsafe { &*self.owner }
    }
}

impl ForceFieldContrib for TorsionAngleContrib {
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
    use super::TorsionAngleContrib;
    use crate::chemistry::forcefield::{
        core::{ForceField, ForceFieldVec3},
        mmff::params::MmffTor,
    };

    const EPS: f64 = 1.0e-12;

    fn force_field_with_positions(count: usize) -> ForceField {
        let mut force_field = ForceField::new(3);
        force_field
            .positions_mut()
            .extend((0..count).map(|idx| ForceFieldVec3::new(idx as f64, 0.0, 0.0)));
        force_field
    }

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= EPS,
            "expected {expected:.16}, got {actual:.16}"
        );
    }

    fn assert_slice_close(actual: &[f64], expected: &[f64]) {
        assert_eq!(actual.len(), expected.len());
        for (idx, (actual_value, expected_value)) in actual.iter().zip(expected.iter()).enumerate()
        {
            assert!(
                (actual_value - expected_value).abs() <= 1.0e-6,
                "index {idx}: expected {expected_value:.16}, got {actual_value:.16}"
            );
        }
    }

    fn finite_difference_gradient(contrib: &TorsionAngleContrib, pos: &[f64]) -> Vec<f64> {
        let mut gradient = vec![0.0; pos.len()];
        let step = 1.0e-6;
        for idx in 0..pos.len() {
            let mut plus = pos.to_vec();
            let mut minus = pos.to_vec();
            plus[idx] += step;
            minus[idx] -= step;
            gradient[idx] = (contrib.get_energy(&plus) - contrib.get_energy(&minus)) / (2.0 * step);
        }
        gradient
    }

    #[test]
    fn mmff_torsionanglecontrib_constructor_stores_owner_pointer() {
        let force_field = ForceField::new(3);

        let contrib = TorsionAngleContrib::new(&force_field);

        assert_eq!(contrib.owner(), &force_field as *const ForceField);
    }

    #[test]
    fn mmff_torsionanglecontrib_constructor_initializes_no_terms() {
        let force_field = ForceField::new(3);

        let contrib = TorsionAngleContrib::new(&force_field);

        assert_eq!(contrib.len(), 0);
        assert!(contrib.is_empty());
    }

    #[test]
    fn mmff_torsionanglecontrib_constructor_accepts_empty_force_field_like_rdkit() {
        let force_field = ForceField::new(3);

        let contrib = TorsionAngleContrib::new(&force_field);

        assert_eq!(force_field.positions().len(), 0);
        assert!(contrib.is_empty());
    }

    #[test]
    fn mmff_torsionanglecontrib_add_term_pushes_source_fields() {
        let force_field = force_field_with_positions(4);
        let mut contrib = TorsionAngleContrib::new(&force_field);
        let tor = MmffTor {
            v1: 1.250,
            v2: 1.125,
            v3: 1.000,
        };

        contrib.add_term(0, 1, 2, 3, &tor);

        assert_eq!(contrib.atom1_indices(), &[0]);
        assert_eq!(contrib.atom2_indices(), &[1]);
        assert_eq!(contrib.atom3_indices(), &[2]);
        assert_eq!(contrib.atom4_indices(), &[3]);
        assert_eq!(contrib.v1(), &[1.250]);
        assert_eq!(contrib.v2(), &[1.125]);
        assert_eq!(contrib.v3(), &[1.000]);
    }

    #[test]
    fn mmff_torsionanglecontrib_add_term_appends_multiple_terms() {
        let force_field = force_field_with_positions(5);
        let mut contrib = TorsionAngleContrib::new(&force_field);
        let first = MmffTor {
            v1: 1.250,
            v2: 1.125,
            v3: 1.000,
        };
        let second = MmffTor {
            v1: 0.500,
            v2: 0.250,
            v3: 0.125,
        };

        contrib.add_term(0, 1, 2, 3, &first);
        contrib.add_term(1, 2, 3, 4, &second);

        assert_eq!(contrib.atom1_indices(), &[0, 1]);
        assert_eq!(contrib.atom2_indices(), &[1, 2]);
        assert_eq!(contrib.atom3_indices(), &[2, 3]);
        assert_eq!(contrib.atom4_indices(), &[3, 4]);
        assert_eq!(contrib.v1(), &[1.250, 0.500]);
        assert_eq!(contrib.v2(), &[1.125, 0.250]);
        assert_eq!(contrib.v3(), &[1.000, 0.125]);
    }

    #[test]
    #[should_panic(expected = "degenerate points")]
    fn mmff_torsionanglecontrib_add_term_rejects_first_second_degenerate() {
        let force_field = force_field_with_positions(4);
        let mut contrib = TorsionAngleContrib::new(&force_field);
        let tor = MmffTor {
            v1: 1.250,
            v2: 1.125,
            v3: 1.000,
        };

        contrib.add_term(0, 0, 2, 3, &tor);
    }

    #[test]
    #[should_panic(expected = "degenerate points")]
    fn mmff_torsionanglecontrib_add_term_rejects_first_third_degenerate() {
        let force_field = force_field_with_positions(4);
        let mut contrib = TorsionAngleContrib::new(&force_field);
        let tor = MmffTor {
            v1: 1.250,
            v2: 1.125,
            v3: 1.000,
        };

        contrib.add_term(0, 1, 0, 3, &tor);
    }

    #[test]
    #[should_panic(expected = "degenerate points")]
    fn mmff_torsionanglecontrib_add_term_rejects_first_fourth_degenerate() {
        let force_field = force_field_with_positions(4);
        let mut contrib = TorsionAngleContrib::new(&force_field);
        let tor = MmffTor {
            v1: 1.250,
            v2: 1.125,
            v3: 1.000,
        };

        contrib.add_term(0, 1, 2, 0, &tor);
    }

    #[test]
    #[should_panic(expected = "degenerate points")]
    fn mmff_torsionanglecontrib_add_term_rejects_second_third_degenerate() {
        let force_field = force_field_with_positions(4);
        let mut contrib = TorsionAngleContrib::new(&force_field);
        let tor = MmffTor {
            v1: 1.250,
            v2: 1.125,
            v3: 1.000,
        };

        contrib.add_term(0, 1, 1, 3, &tor);
    }

    #[test]
    #[should_panic(expected = "degenerate points")]
    fn mmff_torsionanglecontrib_add_term_rejects_second_fourth_degenerate() {
        let force_field = force_field_with_positions(4);
        let mut contrib = TorsionAngleContrib::new(&force_field);
        let tor = MmffTor {
            v1: 1.250,
            v2: 1.125,
            v3: 1.000,
        };

        contrib.add_term(0, 1, 2, 1, &tor);
    }

    #[test]
    #[should_panic(expected = "degenerate points")]
    fn mmff_torsionanglecontrib_add_term_rejects_third_fourth_degenerate() {
        let force_field = force_field_with_positions(4);
        let mut contrib = TorsionAngleContrib::new(&force_field);
        let tor = MmffTor {
            v1: 1.250,
            v2: 1.125,
            v3: 1.000,
        };

        contrib.add_term(0, 1, 2, 2, &tor);
    }

    #[test]
    #[should_panic]
    fn mmff_torsionanglecontrib_add_term_rejects_first_index_out_of_range() {
        let force_field = force_field_with_positions(4);
        let mut contrib = TorsionAngleContrib::new(&force_field);
        let tor = MmffTor {
            v1: 1.250,
            v2: 1.125,
            v3: 1.000,
        };

        contrib.add_term(4, 1, 2, 3, &tor);
    }

    #[test]
    #[should_panic]
    fn mmff_torsionanglecontrib_add_term_rejects_second_index_out_of_range() {
        let force_field = force_field_with_positions(4);
        let mut contrib = TorsionAngleContrib::new(&force_field);
        let tor = MmffTor {
            v1: 1.250,
            v2: 1.125,
            v3: 1.000,
        };

        contrib.add_term(0, 4, 2, 3, &tor);
    }

    #[test]
    #[should_panic]
    fn mmff_torsionanglecontrib_add_term_rejects_third_index_out_of_range() {
        let force_field = force_field_with_positions(4);
        let mut contrib = TorsionAngleContrib::new(&force_field);
        let tor = MmffTor {
            v1: 1.250,
            v2: 1.125,
            v3: 1.000,
        };

        contrib.add_term(0, 1, 4, 3, &tor);
    }

    #[test]
    #[should_panic]
    fn mmff_torsionanglecontrib_add_term_rejects_fourth_index_out_of_range() {
        let force_field = force_field_with_positions(4);
        let mut contrib = TorsionAngleContrib::new(&force_field);
        let tor = MmffTor {
            v1: 1.250,
            v2: 1.125,
            v3: 1.000,
        };

        contrib.add_term(0, 1, 2, 4, &tor);
    }

    #[test]
    fn mmff_torsionanglecontrib_get_energy_returns_zero_without_terms() {
        let force_field = force_field_with_positions(4);
        let contrib = TorsionAngleContrib::new(&force_field);
        let pos = [
            0.0, 0.0, 0.0, //
            1.0, 0.0, 0.0, //
            1.0, 1.0, 0.0, //
            2.0, 1.0, 1.0,
        ];

        assert_eq!(contrib.get_energy(&pos), 0.0);
    }

    #[test]
    fn mmff_torsionanglecontrib_get_energy_uses_rdkit_torsion_formula() {
        let force_field = force_field_with_positions(4);
        let mut contrib = TorsionAngleContrib::new(&force_field);
        let tor = MmffTor {
            v1: 1.5,
            v2: 0.5,
            v3: 0.25,
        };
        let pos = [
            0.0, 0.0, 0.0, //
            1.0, 0.0, 0.0, //
            1.0, 1.0, 0.0, //
            2.0, 1.0, 1.0,
        ];

        contrib.add_term(0, 1, 2, 3, &tor);

        // For this geometry RDKit's cross-product definition gives
        // cos(phi) = -1/sqrt(2).
        let cos_phi = -std::f64::consts::FRAC_1_SQRT_2;
        let cos2_phi = 2.0 * cos_phi * cos_phi - 1.0;
        let cos3_phi = cos_phi * (2.0 * cos2_phi - 1.0);
        let expected = 0.5
            * (tor.v1 * (1.0 + cos_phi) + tor.v2 * (1.0 - cos2_phi) + tor.v3 * (1.0 + cos3_phi));
        assert_close(contrib.get_energy(&pos), expected);
    }

    #[test]
    fn mmff_torsionanglecontrib_get_energy_accumulates_multiple_terms() {
        let force_field = force_field_with_positions(5);
        let mut contrib = TorsionAngleContrib::new(&force_field);
        let first = MmffTor {
            v1: 1.5,
            v2: 0.5,
            v3: 0.25,
        };
        let second = MmffTor {
            v1: 0.8,
            v2: 0.4,
            v3: 0.2,
        };
        let pos = [
            0.0, 0.0, 0.0, //
            1.0, 0.0, 0.0, //
            1.0, 1.0, 0.0, //
            2.0, 1.0, 1.0, //
            3.0, 1.0, 1.0,
        ];

        contrib.add_term(0, 1, 2, 3, &first);
        contrib.add_term(1, 2, 3, 4, &second);

        let first_cos_phi = -std::f64::consts::FRAC_1_SQRT_2;
        let first_cos2_phi = 2.0 * first_cos_phi * first_cos_phi - 1.0;
        let first_cos3_phi = first_cos_phi * (2.0 * first_cos2_phi - 1.0);
        let first_expected = 0.5
            * (first.v1 * (1.0 + first_cos_phi)
                + first.v2 * (1.0 - first_cos2_phi)
                + first.v3 * (1.0 + first_cos3_phi));
        let second_expected = 0.5 * (second.v1 + 2.0 * second.v2 + second.v3);
        assert_close(contrib.get_energy(&pos), first_expected + second_expected);
    }

    #[test]
    fn mmff_torsionanglecontrib_get_energy_returns_zero_cosine_for_degenerate_plane() {
        let force_field = force_field_with_positions(4);
        let mut contrib = TorsionAngleContrib::new(&force_field);
        let tor = MmffTor {
            v1: 1.2,
            v2: 0.4,
            v3: 0.6,
        };
        let pos = [
            0.0, 0.0, 0.0, //
            1.0, 0.0, 0.0, //
            2.0, 0.0, 0.0, //
            2.0, 1.0, 0.0,
        ];

        contrib.add_term(0, 1, 2, 3, &tor);

        // RDKit returns cos(phi)=0.0 when either torsion-plane normal length is
        // effectively zero, instead of using the generic dihedral helper clamp.
        let expected = 0.5 * (tor.v1 + 2.0 * tor.v2 + tor.v3);
        assert_close(contrib.get_energy(&pos), expected);
    }

    #[test]
    fn mmff_torsionanglecontrib_get_grad_leaves_gradient_without_terms() {
        let force_field = ForceField::new(3);
        let contrib = TorsionAngleContrib::new(&force_field);
        let mut grad = vec![1.0, 2.0, 3.0];

        contrib.get_grad(&[], &mut grad);

        assert_eq!(grad, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn mmff_torsionanglecontrib_get_grad_matches_source_energy_derivative() {
        let force_field = force_field_with_positions(4);
        let mut contrib = TorsionAngleContrib::new(&force_field);
        let tor = MmffTor {
            v1: 1.5,
            v2: 0.5,
            v3: 0.25,
        };
        let pos = [
            0.0, 0.0, 0.0, //
            1.0, 0.0, 0.0, //
            1.0, 1.0, 0.0, //
            2.0, 1.0, 1.0,
        ];
        let mut grad = vec![0.0; pos.len()];

        contrib.add_term(0, 1, 2, 3, &tor);
        contrib.get_grad(&pos, &mut grad);

        let expected = finite_difference_gradient(&contrib, &pos);
        assert_slice_close(&grad, &expected);
    }

    #[test]
    fn mmff_torsionanglecontrib_get_grad_accumulates_multiple_terms() {
        let force_field = force_field_with_positions(5);
        let mut contrib = TorsionAngleContrib::new(&force_field);
        let first = MmffTor {
            v1: 1.5,
            v2: 0.5,
            v3: 0.25,
        };
        let second = MmffTor {
            v1: 0.8,
            v2: 0.4,
            v3: 0.2,
        };
        let pos = [
            0.0, 0.0, 0.0, //
            1.0, 0.0, 0.0, //
            1.0, 1.0, 0.0, //
            2.0, 1.0, 1.0, //
            3.0, 1.0, 1.0,
        ];
        let mut grad = vec![0.0; pos.len()];

        contrib.add_term(0, 1, 2, 3, &first);
        contrib.add_term(1, 2, 3, 4, &second);
        contrib.get_grad(&pos, &mut grad);

        let expected = finite_difference_gradient(&contrib, &pos);
        assert_slice_close(&grad, &expected);
    }

    #[test]
    fn mmff_torsionanglecontrib_get_grad_uses_zero_sinphi_source_branch() {
        let force_field = force_field_with_positions(4);
        let mut contrib = TorsionAngleContrib::new(&force_field);
        let tor = MmffTor {
            v1: 1.2,
            v2: 0.4,
            v3: 0.6,
        };
        let pos = [
            0.0, 1.0, 0.0, //
            0.0, 0.0, 0.0, //
            1.0, 0.0, 0.0, //
            1.0, 1.0, 0.0,
        ];
        let mut grad = vec![0.0; pos.len()];

        contrib.add_term(0, 1, 2, 3, &tor);
        contrib.get_grad(&pos, &mut grad);

        assert_eq!(grad, vec![0.0; pos.len()]);
    }

    #[test]
    fn mmff_torsionanglecontrib_get_grad_adds_to_existing_gradient() {
        let force_field = force_field_with_positions(4);
        let mut contrib = TorsionAngleContrib::new(&force_field);
        let tor = MmffTor {
            v1: 1.5,
            v2: 0.5,
            v3: 0.25,
        };
        let pos = [
            0.0, 0.0, 0.0, //
            1.0, 0.0, 0.0, //
            1.0, 1.0, 0.0, //
            2.0, 1.0, 1.0,
        ];
        let mut grad = vec![0.5; pos.len()];

        contrib.add_term(0, 1, 2, 3, &tor);
        contrib.get_grad(&pos, &mut grad);

        let expected_delta = finite_difference_gradient(&contrib, &pos);
        let expected: Vec<f64> = expected_delta
            .into_iter()
            .map(|value| value + 0.5)
            .collect();
        assert_slice_close(&grad, &expected);
    }
}
