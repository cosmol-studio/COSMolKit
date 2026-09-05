//! Source-backed RDKit MMFF angle bend contribution.

use crate::chemistry::forcefield::core::{ForceField, ForceFieldContrib, ForceFieldVec3};

use super::params::{MmffAngle, MmffProp};

// BEGIN RDKIT CPP CONSTANT ForceFields::MMFF angle constants (Params.h:37-39)
// RDKit✔️✔️: constexpr double DEG2RAD = M_PI / 180.0;
const DEG2RAD: f64 = std::f64::consts::PI / 180.0;
// RDKit✔️✔️: constexpr double RAD2DEG = 180.0 / M_PI;
const RAD2DEG: f64 = 180.0 / std::f64::consts::PI;
// RDKit✔️✔️: constexpr double MDYNE_A_TO_KCAL_MOL = 143.9325;
const MDYNE_A_TO_KCAL_MOL: f64 = 143.9325;

#[derive(Clone, Debug)]
pub struct AngleBendContrib {
    owner: *const ForceField,
    is_linear: Vec<bool>,
    atom1_indices: Vec<usize>,
    atom2_indices: Vec<usize>,
    atom3_indices: Vec<usize>,
    ka: Vec<f64>,
    theta0: Vec<f64>,
}

impl AngleBendContrib {
    #[must_use]
    pub fn new(owner: &ForceField) -> Self {
        // BEGIN RDKIT CPP CONSTRUCTOR ForceFields::MMFF::AngleBendContrib::AngleBendContrib (AngleBend.cpp:86-89)
        // RDKit✔️✔️: AngleBendContrib::AngleBendContrib(ForceField *owner) {
        // RDKit✔️✔️:   PRECONDITION(owner, "bad owner");
        // Rust references reproduce RDKit's non-null owner precondition.
        // RDKit✔️✔️:   dp_forceField = owner;
        let owner = owner as *const ForceField;
        // RDKit✔️✔️: }
        Self {
            owner,
            is_linear: Vec::new(),
            atom1_indices: Vec::new(),
            atom2_indices: Vec::new(),
            atom3_indices: Vec::new(),
            ka: Vec::new(),
            theta0: Vec::new(),
        }
    }

    pub fn add_term(
        &mut self,
        idx1: usize,
        idx2: usize,
        idx3: usize,
        mmff_angle_params: &MmffAngle,
        mmff_prop_params_central_atom: &MmffProp,
    ) {
        // BEGIN RDKIT CPP METHOD ForceFields::MMFF::AngleBendContrib::addTerm (AngleBend.cpp:91-111)
        // RDKit✔️✔️: void AngleBendContrib::addTerm(unsigned int idx1,
        // RDKit✔️✔️:                                unsigned int idx2,
        // RDKit✔️✔️:                                unsigned int idx3,
        // RDKit✔️✔️:                                const ForceFields::MMFF::MMFFAngle *mmffAngleParams,
        // RDKit✔️✔️:                                const ForceFields::MMFF::MMFFProp *mmffPropParamsCentralAtom) {
        let force_field = self.force_field();
        // RDKit✔️✔️:   URANGE_CHECK(idx1, dp_forceField->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx2, dp_forceField->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx3, dp_forceField->positions().size());
        assert!(idx1 < force_field.positions().len());
        assert!(idx2 < force_field.positions().len());
        assert!(idx3 < force_field.positions().len());
        // RDKit✔️✔️:   PRECONDITION(((idx1 != idx2) && (idx2 != idx3) && (idx1 != idx3)),
        // RDKit✔️✔️:                "degenerate points");
        assert!(
            idx1 != idx2 && idx2 != idx3 && idx1 != idx3,
            "degenerate points"
        );
        // Rust references reproduce RDKit's non-null angle and central-atom
        // parameter preconditions for the dereferences below.
        // RDKit✔️✔️:   d_at1Idxs.push_back(idx1);
        // RDKit✔️✔️:   d_at2Idxs.push_back(idx2);
        // RDKit✔️✔️:   d_at3Idxs.push_back(idx3);
        self.atom1_indices.push(idx1);
        self.atom2_indices.push(idx2);
        self.atom3_indices.push(idx3);
        // RDKit✔️✔️:   d_isLinear.push_back(mmffPropParamsCentralAtom->linh > 0u);
        // RDKit✔️✔️:   d_theta0.push_back(mmffAngleParams->theta0);
        // RDKit✔️✔️:   d_ka.push_back(mmffAngleParams->ka);
        self.is_linear.push(mmff_prop_params_central_atom.linh > 0);
        self.theta0.push(mmff_angle_params.theta0);
        self.ka.push(mmff_angle_params.ka);
        // RDKit✔️✔️: }
    }

    #[must_use]
    pub fn get_energy(&self, pos: &[f64]) -> f64 {
        // BEGIN RDKIT CPP METHOD ForceFields::MMFF::AngleBendContrib::getEnergy (AngleBend.cpp:108-134)
        // RDKit✔️✔️: double AngleBendContrib::getEnergy(double *pos) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        let force_field = self.force_field();
        // Rust slices reproduce RDKit's non-null pos precondition.

        // RDKit✔️✔️:   double res = 0.0;
        // RDKit✔️✔️:   const int numTerms = d_at1Idxs.size();
        let mut res = 0.0;
        let num_terms = self.atom1_indices.len();

        // RDKit✔️✔️:   for (int i = 0; i < numTerms; i++) {
        for i in 0..num_terms {
            // RDKit✔️✔️:     const int d_at1Idx = d_at1Idxs[i];
            // RDKit✔️✔️:     const int d_at2Idx = d_at2Idxs[i];
            // RDKit✔️✔️:     const int d_at3Idx = d_at3Idxs[i];
            let atom1_idx = self.atom1_indices[i];
            let atom2_idx = self.atom2_indices[i];
            let atom3_idx = self.atom3_indices[i];

            // RDKit✔️✔️:     double dist1 = dp_forceField->distance(d_at1Idx, d_at2Idx, pos);
            // RDKit✔️✔️:     double dist2 = dp_forceField->distance(d_at2Idx, d_at3Idx, pos);
            let dist1 = force_field.distance(atom1_idx, atom2_idx, Some(pos));
            let dist2 = force_field.distance(atom2_idx, atom3_idx, Some(pos));

            // RDKit✔️✔️:     RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
            // RDKit✔️✔️:                        pos[3 * d_at1Idx + 2]);
            let p1 = ForceFieldVec3::new(
                pos[3 * atom1_idx],
                pos[3 * atom1_idx + 1],
                pos[3 * atom1_idx + 2],
            );
            // RDKit✔️✔️:     RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
            // RDKit✔️✔️:                        pos[3 * d_at2Idx + 2]);
            let p2 = ForceFieldVec3::new(
                pos[3 * atom2_idx],
                pos[3 * atom2_idx + 1],
                pos[3 * atom2_idx + 2],
            );
            // RDKit✔️✔️:     RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
            // RDKit✔️✔️:                        pos[3 * d_at3Idx + 2]);
            let p3 = ForceFieldVec3::new(
                pos[3 * atom3_idx],
                pos[3 * atom3_idx + 1],
                pos[3 * atom3_idx + 2],
            );

            // RDKit✔️✔️:     res += Utils::calcAngleBendEnergy(
            // RDKit✔️✔️:         d_theta0[i], d_ka[i], d_isLinear[i],
            // RDKit✔️✔️:         Utils::calcCosTheta(p1, p2, p3, dist1, dist2));
            res += calc_angle_bend_energy(
                self.theta0[i],
                self.ka[i],
                self.is_linear[i],
                calc_cos_theta(p1, p2, p3, dist1, dist2),
            );
            // RDKit✔️✔️:   }
        }
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        res
    }

    pub fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
        // BEGIN RDKIT CPP METHOD ForceFields::MMFF::AngleBendContrib::getGrad (AngleBend.cpp:136-179)
        // RDKit✔️✔️: void AngleBendContrib::getGrad(double *pos, double *grad) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:   PRECONDITION(grad, "bad vector");
        let force_field = self.force_field();
        // Rust slices reproduce RDKit's non-null pos and grad preconditions.

        // RDKit✔️✔️:   const int numTerms = d_at1Idxs.size();
        let num_terms = self.atom1_indices.len();
        // RDKit✔️✔️:   for (int i =0; i < numTerms; i++) {
        for i in 0..num_terms {
            // RDKit✔️✔️:     const int d_at1Idx = d_at1Idxs[i];
            // RDKit✔️✔️:     const int d_at2Idx = d_at2Idxs[i];
            // RDKit✔️✔️:     const int d_at3Idx = d_at3Idxs[i];
            let atom1_idx = self.atom1_indices[i];
            let atom2_idx = self.atom2_indices[i];
            let atom3_idx = self.atom3_indices[i];

            // RDKit✔️✔️:     double dist[2] = {dp_forceField->distance(d_at1Idx, d_at2Idx, pos),
            // RDKit✔️✔️:                       dp_forceField->distance(d_at2Idx, d_at3Idx, pos)};
            let dist = [
                force_field.distance(atom1_idx, atom2_idx, Some(pos)),
                force_field.distance(atom2_idx, atom3_idx, Some(pos)),
            ];

            // RDKit✔️✔️:     RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
            // RDKit✔️✔️:                        pos[3 * d_at1Idx + 2]);
            let p1 = ForceFieldVec3::new(
                pos[3 * atom1_idx],
                pos[3 * atom1_idx + 1],
                pos[3 * atom1_idx + 2],
            );
            // RDKit✔️✔️:     RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
            // RDKit✔️✔️:                        pos[3 * d_at2Idx + 2]);
            let p2 = ForceFieldVec3::new(
                pos[3 * atom2_idx],
                pos[3 * atom2_idx + 1],
                pos[3 * atom2_idx + 2],
            );
            // RDKit✔️✔️:     RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
            // RDKit✔️✔️:                        pos[3 * d_at3Idx + 2]);
            let p3 = ForceFieldVec3::new(
                pos[3 * atom3_idx],
                pos[3 * atom3_idx + 1],
                pos[3 * atom3_idx + 2],
            );

            // RDKit✔️✔️:     double *g[3] = {&(grad[3 * d_at1Idx]), &(grad[3 * d_at2Idx]),
            // RDKit✔️✔️:                     &(grad[3 * d_at3Idx])};
            // Rust passes atom indices to calc_angle_bend_grad and indexes the
            // flat gradient slice at the same offsets.
            // RDKit✔️✔️:     RDGeom::Point3D r[2] = {(p1 - p2) / dist[0], (p3 - p2) / dist[1]};
            let r = [(p1 - p2) / dist[0], (p3 - p2) / dist[1]];
            // RDKit✔️✔️:     double cosTheta = r[0].dotProduct(r[1]);
            let mut cos_theta = r[0].dot_product(r[1]);
            // RDKit✔️✔️:     clipToOne(cosTheta);
            cos_theta = cos_theta.clamp(-1.0, 1.0);
            // RDKit✔️✔️:     double sinThetaSq = 1.0 - cosTheta * cosTheta;
            let sin_theta_sq = 1.0 - cos_theta * cos_theta;
            // RDKit✔️✔️:     double sinTheta =
            // RDKit✔️✔️:         std::max(((sinThetaSq > 0.0) ? sqrt(sinThetaSq) : 0.0), 1.0e-8);
            let sin_theta = (if sin_theta_sq > 0.0 {
                sin_theta_sq.sqrt()
            } else {
                0.0
            })
            .max(1.0e-8);

            // RDKit✔️✔️:
            // RDKit✔️✔️:     // use the chain rule:
            // RDKit✔️✔️:     // dE/dx = dE/dTheta * dTheta/dx
            // RDKit✔️✔️:
            // RDKit✔️✔️:     // dE/dTheta is independent of cartesians:
            // RDKit✔️✔️:     double angleTerm = RAD2DEG * acos(cosTheta) - d_theta0[i];
            let angle_term = RAD2DEG * cos_theta.acos() - self.theta0[i];
            // RDKit✔️✔️:     double const cb = -0.006981317;
            // RDKit✔️✔️:     double const c2 = MDYNE_A_TO_KCAL_MOL * DEG2RAD * DEG2RAD;
            let cb = -0.006981317;
            let c2 = MDYNE_A_TO_KCAL_MOL * DEG2RAD * DEG2RAD;

            // RDKit✔️✔️:
            // RDKit✔️✔️:     double dE_dTheta = (d_isLinear[i] ? -MDYNE_A_TO_KCAL_MOL * d_ka[i] * sinTheta
            // RDKit✔️✔️:                                    : RAD2DEG * c2 * d_ka[i] * angleTerm *
            // RDKit✔️✔️:                                          (1.0 + 1.5 * cb * angleTerm));
            let de_dtheta = if self.is_linear[i] {
                -MDYNE_A_TO_KCAL_MOL * self.ka[i] * sin_theta
            } else {
                RAD2DEG * c2 * self.ka[i] * angle_term * (1.0 + 1.5 * cb * angle_term)
            };

            // RDKit✔️✔️:
            // RDKit✔️✔️:     Utils::calcAngleBendGrad(r, dist, g, dE_dTheta, cosTheta, sinTheta);
            calc_angle_bend_grad(
                &r,
                &dist,
                grad,
                [atom1_idx, atom2_idx, atom3_idx],
                de_dtheta,
                cos_theta,
                sin_theta,
            );
            // RDKit✔️✔️:   }
        }
        // RDKit✔️✔️: }
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
        self.is_linear.is_empty()
            && self.atom1_indices.is_empty()
            && self.atom2_indices.is_empty()
            && self.atom3_indices.is_empty()
            && self.ka.is_empty()
            && self.theta0.is_empty()
    }

    #[must_use]
    pub fn is_linear_terms(&self) -> &[bool] {
        &self.is_linear
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
    pub fn force_constants(&self) -> &[f64] {
        &self.ka
    }

    #[must_use]
    pub fn rest_angles(&self) -> &[f64] {
        &self.theta0
    }

    fn force_field(&self) -> &ForceField {
        assert!(!self.owner.is_null(), "no owner");
        // SAFETY: MMFF contribs follow the same owner-pointer model as
        // ForceFieldContrib objects in core.rs. Constructors store a live
        // ForceField pointer before term insertion or evaluation.
        unsafe { &*self.owner }
    }
}

fn calc_cos_theta(
    p1: ForceFieldVec3,
    p2: ForceFieldVec3,
    p3: ForceFieldVec3,
    dist1: f64,
    dist2: f64,
) -> f64 {
    // BEGIN RDKIT CPP HELPER ForceFields::MMFF::Utils::calcCosTheta (AngleBend.cpp:27-35)
    // RDKit✔️✔️: double calcCosTheta(RDGeom::Point3D p1, RDGeom::Point3D p2, RDGeom::Point3D p3,
    // RDKit✔️✔️:                     double dist1, double dist2) {
    // RDKit✔️✔️:   RDGeom::Point3D p12 = p1 - p2;
    // RDKit✔️✔️:   RDGeom::Point3D p32 = p3 - p2;
    let p12 = p1 - p2;
    let p32 = p3 - p2;
    // RDKit✔️✔️:   double cosTheta = p12.dotProduct(p32) / (dist1 * dist2);
    let mut cos_theta = p12.dot_product(p32) / (dist1 * dist2);
    // RDKit✔️✔️:   clipToOne(cosTheta);
    cos_theta = cos_theta.clamp(-1.0, 1.0);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return cosTheta;
    // RDKit✔️✔️: }
    cos_theta
}

fn calc_angle_rest_value(mmff_angle_params: &MmffAngle) -> f64 {
    // BEGIN RDKIT CPP HELPER ForceFields::MMFF::Utils::calcAngleRestValue (AngleBend.cpp:21-25)
    // RDKit✔️✔️: double calcAngleRestValue(const MMFFAngle *mmffAngleParams) {
    // RDKit✔️✔️:   PRECONDITION(mmffAngleParams, "angle parameters not found");
    // Rust references reproduce RDKit's non-null parameter precondition.
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return mmffAngleParams->theta0;
    // RDKit✔️✔️: }
    mmff_angle_params.theta0
}

fn calc_angle_force_constant(mmff_angle_params: &MmffAngle) -> f64 {
    // BEGIN RDKIT CPP HELPER ForceFields::MMFF::Utils::calcAngleForceConstant (AngleBend.cpp:37-41)
    // RDKit✔️✔️: double calcAngleForceConstant(const MMFFAngle *mmffAngleParams) {
    // RDKit✔️✔️:   PRECONDITION(mmffAngleParams, "angle parameters not found");
    // Rust references reproduce RDKit's non-null parameter precondition.
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return mmffAngleParams->ka;
    // RDKit✔️✔️: }
    mmff_angle_params.ka
}

fn calc_angle_bend_energy(theta0: f64, ka: f64, is_linear: bool, cos_theta: f64) -> f64 {
    // BEGIN RDKIT CPP HELPER ForceFields::MMFF::Utils::calcAngleBendEnergy (AngleBend.cpp:43-57)
    // RDKit✔️✔️: double calcAngleBendEnergy(const double theta0, const double ka, bool isLinear,
    // RDKit✔️✔️:                            const double cosTheta) {
    // RDKit✔️✔️:   double angle = RAD2DEG * acos(cosTheta) - theta0;
    let angle = RAD2DEG * cos_theta.acos() - theta0;
    // RDKit✔️✔️:   double const cb = -0.006981317;
    // RDKit✔️✔️:   double const c2 = MDYNE_A_TO_KCAL_MOL * DEG2RAD * DEG2RAD;
    // RDKit✔️✔️:   double res = 0.0;
    let cb = -0.006981317;
    let c2 = MDYNE_A_TO_KCAL_MOL * DEG2RAD * DEG2RAD;
    let res;

    // RDKit✔️✔️:   if (isLinear) {
    if is_linear {
        // RDKit✔️✔️:     res = MDYNE_A_TO_KCAL_MOL * ka * (1.0 + cosTheta);
        res = MDYNE_A_TO_KCAL_MOL * ka * (1.0 + cos_theta);
        // RDKit✔️✔️:   } else {
    } else {
        // RDKit✔️✔️:     res = 0.5 * c2 * ka * angle * angle * (1.0 + cb * angle);
        res = 0.5 * c2 * ka * angle * angle * (1.0 + cb * angle);
        // RDKit✔️✔️:   }
    }

    // RDKit✔️✔️:
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    res
}

fn calc_angle_bend_grad(
    r: &[ForceFieldVec3; 2],
    dist: &[f64; 2],
    grad: &mut [f64],
    idx: [usize; 3],
    de_dtheta: f64,
    cos_theta: f64,
    sin_theta: f64,
) {
    // BEGIN RDKIT CPP HELPER ForceFields::MMFF::Utils::calcAngleBendGrad (AngleBend.cpp:59-81)
    // RDKit✔️✔️: void calcAngleBendGrad(RDGeom::Point3D *r, double *dist, double **g,
    // RDKit✔️✔️:                        double &dE_dTheta, double &cosTheta, double &sinTheta) {
    // RDKit✔️✔️:   // -------
    // RDKit✔️✔️:   // dTheta/dx is trickier:
    // RDKit✔️✔️:   double dCos_dS[6] = {1.0 / dist[0] * (r[1].x - cosTheta * r[0].x),
    // RDKit✔️✔️:                        1.0 / dist[0] * (r[1].y - cosTheta * r[0].y),
    // RDKit✔️✔️:                        1.0 / dist[0] * (r[1].z - cosTheta * r[0].z),
    // RDKit✔️✔️:                        1.0 / dist[1] * (r[0].x - cosTheta * r[1].x),
    // RDKit✔️✔️:                        1.0 / dist[1] * (r[0].y - cosTheta * r[1].y),
    // RDKit✔️✔️:                        1.0 / dist[1] * (r[0].z - cosTheta * r[1].z)};
    let dcos_ds = [
        1.0 / dist[0] * (r[1].x - cos_theta * r[0].x),
        1.0 / dist[0] * (r[1].y - cos_theta * r[0].y),
        1.0 / dist[0] * (r[1].z - cos_theta * r[0].z),
        1.0 / dist[1] * (r[0].x - cos_theta * r[1].x),
        1.0 / dist[1] * (r[0].y - cos_theta * r[1].y),
        1.0 / dist[1] * (r[0].z - cos_theta * r[1].z),
    ];
    let g0 = 3 * idx[0];
    let g1 = 3 * idx[1];
    let g2 = 3 * idx[2];

    // RDKit✔️✔️:   g[0][0] += dE_dTheta * dCos_dS[0] / (-sinTheta);
    // RDKit✔️✔️:   g[0][1] += dE_dTheta * dCos_dS[1] / (-sinTheta);
    // RDKit✔️✔️:   g[0][2] += dE_dTheta * dCos_dS[2] / (-sinTheta);
    grad[g0] += de_dtheta * dcos_ds[0] / (-sin_theta);
    grad[g0 + 1] += de_dtheta * dcos_ds[1] / (-sin_theta);
    grad[g0 + 2] += de_dtheta * dcos_ds[2] / (-sin_theta);

    // RDKit✔️✔️:   g[1][0] += dE_dTheta * (-dCos_dS[0] - dCos_dS[3]) / (-sinTheta);
    // RDKit✔️✔️:   g[1][1] += dE_dTheta * (-dCos_dS[1] - dCos_dS[4]) / (-sinTheta);
    // RDKit✔️✔️:   g[1][2] += dE_dTheta * (-dCos_dS[2] - dCos_dS[5]) / (-sinTheta);
    grad[g1] += de_dtheta * (-dcos_ds[0] - dcos_ds[3]) / (-sin_theta);
    grad[g1 + 1] += de_dtheta * (-dcos_ds[1] - dcos_ds[4]) / (-sin_theta);
    grad[g1 + 2] += de_dtheta * (-dcos_ds[2] - dcos_ds[5]) / (-sin_theta);

    // RDKit✔️✔️:   g[2][0] += dE_dTheta * dCos_dS[3] / (-sinTheta);
    // RDKit✔️✔️:   g[2][1] += dE_dTheta * dCos_dS[4] / (-sinTheta);
    // RDKit✔️✔️:   g[2][2] += dE_dTheta * dCos_dS[5] / (-sinTheta);
    // RDKit✔️✔️: }
    grad[g2] += de_dtheta * dcos_ds[3] / (-sin_theta);
    grad[g2 + 1] += de_dtheta * dcos_ds[4] / (-sin_theta);
    grad[g2 + 2] += de_dtheta * dcos_ds[5] / (-sin_theta);
}

#[cfg(test)]
mod tests {
    use super::{
        AngleBendContrib, calc_angle_bend_energy, calc_angle_bend_grad, calc_angle_force_constant,
        calc_angle_rest_value, calc_cos_theta,
    };
    use crate::chemistry::forcefield::{
        core::{ForceField, ForceFieldVec3},
        mmff::params::{MmffAngle, MmffProp},
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

    fn source_angle_params(ka: f64, theta0: f64) -> MmffAngle {
        MmffAngle { ka, theta0 }
    }

    fn source_central_atom_prop(linh: u8) -> MmffProp {
        MmffProp {
            atno: 6,
            crd: 4,
            val: 4,
            pilp: 0,
            mltb: 0,
            arom: 0,
            linh,
            sbmb: 0,
        }
    }

    fn source_angle_bend_energy(theta0: f64, ka: f64, is_linear: bool, cos_theta: f64) -> f64 {
        let angle = 180.0 / std::f64::consts::PI * cos_theta.acos() - theta0;
        let cb = -0.006981317;
        let c2 = 143.9325 * (std::f64::consts::PI / 180.0) * (std::f64::consts::PI / 180.0);
        if is_linear {
            143.9325 * ka * (1.0 + cos_theta)
        } else {
            0.5 * c2 * ka * angle * angle * (1.0 + cb * angle)
        }
    }

    fn source_angle_bend_grad(
        pos: &[f64],
        idx: [usize; 3],
        theta0: f64,
        ka: f64,
        is_linear: bool,
    ) -> Vec<f64> {
        let p1 = ForceFieldVec3::new(pos[3 * idx[0]], pos[3 * idx[0] + 1], pos[3 * idx[0] + 2]);
        let p2 = ForceFieldVec3::new(pos[3 * idx[1]], pos[3 * idx[1] + 1], pos[3 * idx[1] + 2]);
        let p3 = ForceFieldVec3::new(pos[3 * idx[2]], pos[3 * idx[2] + 1], pos[3 * idx[2] + 2]);
        let dist = [(p1 - p2).length(), (p3 - p2).length()];
        let r = [(p1 - p2) / dist[0], (p3 - p2) / dist[1]];
        let mut cos_theta = r[0].dot_product(r[1]).clamp(-1.0, 1.0);
        if cos_theta.is_nan() {
            cos_theta = f64::NAN;
        }
        let sin_theta_sq = 1.0 - cos_theta * cos_theta;
        let sin_theta = (if sin_theta_sq > 0.0 {
            sin_theta_sq.sqrt()
        } else {
            0.0
        })
        .max(1.0e-8);
        let angle_term = 180.0 / std::f64::consts::PI * cos_theta.acos() - theta0;
        let cb = -0.006981317;
        let c2 = 143.9325 * (std::f64::consts::PI / 180.0) * (std::f64::consts::PI / 180.0);
        let de_dtheta = if is_linear {
            -143.9325 * ka * sin_theta
        } else {
            180.0 / std::f64::consts::PI * c2 * ka * angle_term * (1.0 + 1.5 * cb * angle_term)
        };
        let dcos_ds = [
            1.0 / dist[0] * (r[1].x - cos_theta * r[0].x),
            1.0 / dist[0] * (r[1].y - cos_theta * r[0].y),
            1.0 / dist[0] * (r[1].z - cos_theta * r[0].z),
            1.0 / dist[1] * (r[0].x - cos_theta * r[1].x),
            1.0 / dist[1] * (r[0].y - cos_theta * r[1].y),
            1.0 / dist[1] * (r[0].z - cos_theta * r[1].z),
        ];
        let mut grad = vec![0.0; pos.len()];
        let scale = de_dtheta / (-sin_theta);
        let g0 = 3 * idx[0];
        let g1 = 3 * idx[1];
        let g2 = 3 * idx[2];
        grad[g0] += scale * dcos_ds[0];
        grad[g0 + 1] += scale * dcos_ds[1];
        grad[g0 + 2] += scale * dcos_ds[2];
        grad[g1] += scale * (-dcos_ds[0] - dcos_ds[3]);
        grad[g1 + 1] += scale * (-dcos_ds[1] - dcos_ds[4]);
        grad[g1 + 2] += scale * (-dcos_ds[2] - dcos_ds[5]);
        grad[g2] += scale * dcos_ds[3];
        grad[g2 + 1] += scale * dcos_ds[4];
        grad[g2 + 2] += scale * dcos_ds[5];
        grad
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
                (actual - expected).abs() < 1.0e-10,
                "idx={idx}, actual={actual}, expected={expected}"
            );
        }
    }

    #[test]
    fn mmff_anglebendcontrib_constructor_stores_owner_pointer() {
        let force_field = ForceField::new(3);

        let contrib = AngleBendContrib::new(&force_field);

        assert_eq!(contrib.owner(), &force_field as *const ForceField);
    }

    #[test]
    fn mmff_anglebendcontrib_constructor_initializes_no_terms() {
        let force_field = ForceField::new(3);

        let contrib = AngleBendContrib::new(&force_field);

        assert_eq!(contrib.len(), 0);
        assert!(contrib.is_empty());
        assert!(contrib.is_linear_terms().is_empty());
        assert!(contrib.atom1_indices().is_empty());
        assert!(contrib.atom2_indices().is_empty());
        assert!(contrib.atom3_indices().is_empty());
        assert!(contrib.force_constants().is_empty());
        assert!(contrib.rest_angles().is_empty());
    }

    #[test]
    fn mmff_anglebendcontrib_constructor_accepts_empty_force_field_like_rdkit() {
        let force_field = ForceField::new(3);

        let contrib = AngleBendContrib::new(&force_field);

        assert_eq!(force_field.positions().len(), 0);
        assert!(contrib.is_empty());
    }

    #[test]
    fn mmff_anglebendcontrib_add_term_pushes_source_fields_for_nonlinear_center() {
        let force_field = force_field_with_positions(3);
        let mut contrib = AngleBendContrib::new(&force_field);
        let angle_params = source_angle_params(0.636, 110.549);
        let central_atom_prop = source_central_atom_prop(0);

        contrib.add_term(0, 1, 2, &angle_params, &central_atom_prop);

        assert_eq!(contrib.atom1_indices(), &[0]);
        assert_eq!(contrib.atom2_indices(), &[1]);
        assert_eq!(contrib.atom3_indices(), &[2]);
        assert_eq!(contrib.is_linear_terms(), &[false]);
        assert_eq!(contrib.rest_angles(), &[110.549]);
        assert_eq!(contrib.force_constants(), &[0.636]);
    }

    #[test]
    fn mmff_anglebendcontrib_add_term_marks_linh_positive_as_linear() {
        let force_field = force_field_with_positions(3);
        let mut contrib = AngleBendContrib::new(&force_field);
        let angle_params = source_angle_params(1.006, 180.0);
        let central_atom_prop = source_central_atom_prop(1);

        contrib.add_term(0, 1, 2, &angle_params, &central_atom_prop);

        assert_eq!(contrib.is_linear_terms(), &[true]);
    }

    #[test]
    fn mmff_anglebendcontrib_add_term_appends_multiple_terms_in_source_order() {
        let force_field = force_field_with_positions(4);
        let mut contrib = AngleBendContrib::new(&force_field);
        let first_angle = source_angle_params(0.636, 110.549);
        let second_angle = source_angle_params(1.006, 120.0);
        let nonlinear_prop = source_central_atom_prop(0);
        let linear_prop = source_central_atom_prop(2);

        contrib.add_term(0, 1, 2, &first_angle, &nonlinear_prop);
        contrib.add_term(1, 2, 3, &second_angle, &linear_prop);

        assert_eq!(contrib.atom1_indices(), &[0, 1]);
        assert_eq!(contrib.atom2_indices(), &[1, 2]);
        assert_eq!(contrib.atom3_indices(), &[2, 3]);
        assert_eq!(contrib.is_linear_terms(), &[false, true]);
        assert_eq!(contrib.rest_angles(), &[110.549, 120.0]);
        assert_eq!(contrib.force_constants(), &[0.636, 1.006]);
    }

    #[test]
    #[should_panic]
    fn mmff_anglebendcontrib_add_term_panics_for_out_of_range_idx1() {
        let force_field = force_field_with_positions(3);
        let mut contrib = AngleBendContrib::new(&force_field);
        let angle_params = source_angle_params(0.636, 110.549);
        let central_atom_prop = source_central_atom_prop(0);

        contrib.add_term(3, 1, 2, &angle_params, &central_atom_prop);
    }

    #[test]
    #[should_panic]
    fn mmff_anglebendcontrib_add_term_panics_for_out_of_range_idx2() {
        let force_field = force_field_with_positions(3);
        let mut contrib = AngleBendContrib::new(&force_field);
        let angle_params = source_angle_params(0.636, 110.549);
        let central_atom_prop = source_central_atom_prop(0);

        contrib.add_term(0, 3, 2, &angle_params, &central_atom_prop);
    }

    #[test]
    #[should_panic]
    fn mmff_anglebendcontrib_add_term_panics_for_out_of_range_idx3() {
        let force_field = force_field_with_positions(3);
        let mut contrib = AngleBendContrib::new(&force_field);
        let angle_params = source_angle_params(0.636, 110.549);
        let central_atom_prop = source_central_atom_prop(0);

        contrib.add_term(0, 1, 3, &angle_params, &central_atom_prop);
    }

    #[test]
    #[should_panic(expected = "degenerate points")]
    fn mmff_anglebendcontrib_add_term_panics_for_idx1_equal_idx2() {
        let force_field = force_field_with_positions(3);
        let mut contrib = AngleBendContrib::new(&force_field);
        let angle_params = source_angle_params(0.636, 110.549);
        let central_atom_prop = source_central_atom_prop(0);

        contrib.add_term(1, 1, 2, &angle_params, &central_atom_prop);
    }

    #[test]
    #[should_panic(expected = "degenerate points")]
    fn mmff_anglebendcontrib_add_term_panics_for_idx2_equal_idx3() {
        let force_field = force_field_with_positions(3);
        let mut contrib = AngleBendContrib::new(&force_field);
        let angle_params = source_angle_params(0.636, 110.549);
        let central_atom_prop = source_central_atom_prop(0);

        contrib.add_term(0, 2, 2, &angle_params, &central_atom_prop);
    }

    #[test]
    #[should_panic(expected = "degenerate points")]
    fn mmff_anglebendcontrib_add_term_panics_for_idx1_equal_idx3() {
        let force_field = force_field_with_positions(3);
        let mut contrib = AngleBendContrib::new(&force_field);
        let angle_params = source_angle_params(0.636, 110.549);
        let central_atom_prop = source_central_atom_prop(0);

        contrib.add_term(0, 1, 0, &angle_params, &central_atom_prop);
    }

    #[test]
    fn mmff_anglebendcontrib_get_energy_returns_zero_without_terms() {
        let force_field = ForceField::new(3);
        let contrib = AngleBendContrib::new(&force_field);

        assert_eq!(contrib.get_energy(&[]), 0.0);
    }

    #[test]
    fn mmff_anglebendcontrib_get_energy_uses_source_formula_for_nonlinear_term() {
        let force_field = initialized_force_field_with_positions(3);
        let mut contrib = AngleBendContrib::new(&force_field);
        let angle_params = source_angle_params(0.636, 80.0);
        let central_atom_prop = source_central_atom_prop(0);
        contrib.add_term(0, 1, 2, &angle_params, &central_atom_prop);
        let pos = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0];

        let expected = source_angle_bend_energy(80.0, 0.636, false, 0.0);

        assert_close(contrib.get_energy(&pos), expected);
    }

    #[test]
    fn mmff_anglebendcontrib_get_energy_uses_source_formula_for_linear_term() {
        let force_field = initialized_force_field_with_positions(3);
        let mut contrib = AngleBendContrib::new(&force_field);
        let angle_params = source_angle_params(1.006, 180.0);
        let central_atom_prop = source_central_atom_prop(1);
        contrib.add_term(0, 1, 2, &angle_params, &central_atom_prop);
        let pos = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0];

        let expected = source_angle_bend_energy(180.0, 1.006, true, 0.0);

        assert_close(contrib.get_energy(&pos), expected);
    }

    #[test]
    fn mmff_anglebendcontrib_get_energy_sums_multiple_terms() {
        let force_field = initialized_force_field_with_positions(4);
        let mut contrib = AngleBendContrib::new(&force_field);
        let nonlinear_angle = source_angle_params(0.636, 80.0);
        let linear_angle = source_angle_params(1.006, 180.0);
        let nonlinear_prop = source_central_atom_prop(0);
        let linear_prop = source_central_atom_prop(2);
        contrib.add_term(0, 1, 2, &nonlinear_angle, &nonlinear_prop);
        contrib.add_term(1, 2, 3, &linear_angle, &linear_prop);
        let pos = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 0.0];

        let first = source_angle_bend_energy(80.0, 0.636, false, 0.0);
        let second = source_angle_bend_energy(180.0, 1.006, true, 0.0);

        assert_close(contrib.get_energy(&pos), first + second);
    }

    #[test]
    fn mmff_anglebendcontrib_get_energy_uses_supplied_position_vector() {
        let force_field = initialized_force_field_with_positions(3);
        let mut contrib = AngleBendContrib::new(&force_field);
        let angle_params = source_angle_params(0.636, 100.0);
        let central_atom_prop = source_central_atom_prop(0);
        contrib.add_term(0, 1, 2, &angle_params, &central_atom_prop);
        let pos = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0];

        let expected = source_angle_bend_energy(100.0, 0.636, false, 0.0);

        assert_close(contrib.get_energy(&pos), expected);
    }

    #[test]
    #[should_panic(expected = "not initialized")]
    fn mmff_anglebendcontrib_get_energy_requires_initialized_force_field_for_terms() {
        let force_field = force_field_with_positions(3);
        let mut contrib = AngleBendContrib::new(&force_field);
        let angle_params = source_angle_params(0.636, 100.0);
        let central_atom_prop = source_central_atom_prop(0);
        contrib.add_term(0, 1, 2, &angle_params, &central_atom_prop);

        let _ = contrib.get_energy(&[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]);
    }

    #[test]
    fn mmff_anglebendcontrib_calc_cos_theta_clips_to_source_bounds() {
        let p1 = ForceFieldVec3::new(1.0, 0.0, 0.0);
        let p2 = ForceFieldVec3::new(0.0, 0.0, 0.0);
        let p3 = ForceFieldVec3::new(1.0, 0.0, 0.0);

        assert_eq!(calc_cos_theta(p1, p2, p3, 0.5, 0.5), 1.0);
    }

    #[test]
    fn mmff_anglebendcontrib_calc_cos_theta_uses_source_dot_over_distances() {
        let p1 = ForceFieldVec3::new(1.0, 1.0, 0.0);
        let p2 = ForceFieldVec3::new(0.0, 0.0, 0.0);
        let p3 = ForceFieldVec3::new(1.0, 0.0, 0.0);

        assert_close(
            calc_cos_theta(p1, p2, p3, 2.0_f64.sqrt(), 1.0),
            0.5_f64.sqrt(),
        );
    }

    #[test]
    fn mmff_anglebendcontrib_calc_cos_theta_clips_lower_source_bound() {
        let p1 = ForceFieldVec3::new(-1.0, 0.0, 0.0);
        let p2 = ForceFieldVec3::new(0.0, 0.0, 0.0);
        let p3 = ForceFieldVec3::new(1.0, 0.0, 0.0);

        assert_eq!(calc_cos_theta(p1, p2, p3, 0.5, 0.5), -1.0);
    }

    #[test]
    fn mmff_anglebendcontrib_calc_angle_rest_value_returns_source_theta0() {
        let angle_params = source_angle_params(0.636, 110.549);

        assert_eq!(calc_angle_rest_value(&angle_params), 110.549);
    }

    #[test]
    fn mmff_anglebendcontrib_calc_angle_rest_value_ignores_force_constant_like_source() {
        let low_force_constant = source_angle_params(0.0, 120.0);
        let high_force_constant = source_angle_params(99.0, 120.0);

        assert_eq!(calc_angle_rest_value(&low_force_constant), 120.0);
        assert_eq!(calc_angle_rest_value(&high_force_constant), 120.0);
    }

    #[test]
    fn mmff_anglebendcontrib_calc_angle_force_constant_returns_source_ka() {
        let angle_params = source_angle_params(0.636, 110.549);

        assert_eq!(calc_angle_force_constant(&angle_params), 0.636);
    }

    #[test]
    fn mmff_anglebendcontrib_calc_angle_force_constant_ignores_rest_angle_like_source() {
        let low_rest_angle = source_angle_params(1.006, 0.0);
        let high_rest_angle = source_angle_params(1.006, 180.0);

        assert_eq!(calc_angle_force_constant(&low_rest_angle), 1.006);
        assert_eq!(calc_angle_force_constant(&high_rest_angle), 1.006);
    }

    #[test]
    fn mmff_anglebendcontrib_calc_angle_bend_energy_uses_nonlinear_source_branch() {
        let expected = source_angle_bend_energy(80.0, 0.636, false, 0.25);

        assert_close(calc_angle_bend_energy(80.0, 0.636, false, 0.25), expected);
    }

    #[test]
    fn mmff_anglebendcontrib_calc_angle_bend_energy_uses_linear_source_branch() {
        let expected = source_angle_bend_energy(180.0, 1.006, true, -0.25);

        assert_close(calc_angle_bend_energy(180.0, 1.006, true, -0.25), expected);
    }

    #[test]
    fn mmff_anglebendcontrib_get_grad_leaves_gradient_without_terms() {
        let force_field = ForceField::new(3);
        let contrib = AngleBendContrib::new(&force_field);
        let mut grad = vec![1.0, 2.0, 3.0];

        contrib.get_grad(&[], &mut grad);

        assert_eq!(grad, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn mmff_anglebendcontrib_get_grad_accumulates_nonlinear_source_gradient() {
        let force_field = initialized_force_field_with_positions(3);
        let mut contrib = AngleBendContrib::new(&force_field);
        let angle_params = source_angle_params(0.636, 80.0);
        let central_atom_prop = source_central_atom_prop(0);
        contrib.add_term(0, 1, 2, &angle_params, &central_atom_prop);
        let pos = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let mut grad = vec![0.0; pos.len()];

        contrib.get_grad(&pos, &mut grad);

        let expected = source_angle_bend_grad(&pos, [0, 1, 2], 80.0, 0.636, false);
        assert_slice_close(&grad, &expected);
    }

    #[test]
    fn mmff_anglebendcontrib_get_grad_accumulates_linear_source_gradient() {
        let force_field = initialized_force_field_with_positions(3);
        let mut contrib = AngleBendContrib::new(&force_field);
        let angle_params = source_angle_params(1.006, 180.0);
        let central_atom_prop = source_central_atom_prop(1);
        contrib.add_term(0, 1, 2, &angle_params, &central_atom_prop);
        let pos = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let mut grad = vec![0.0; pos.len()];

        contrib.get_grad(&pos, &mut grad);

        let expected = source_angle_bend_grad(&pos, [0, 1, 2], 180.0, 1.006, true);
        assert_slice_close(&grad, &expected);
    }

    #[test]
    fn mmff_anglebendcontrib_get_grad_accumulates_multiple_terms() {
        let force_field = initialized_force_field_with_positions(4);
        let mut contrib = AngleBendContrib::new(&force_field);
        let nonlinear_angle = source_angle_params(0.636, 80.0);
        let linear_angle = source_angle_params(1.006, 180.0);
        let nonlinear_prop = source_central_atom_prop(0);
        let linear_prop = source_central_atom_prop(2);
        contrib.add_term(0, 1, 2, &nonlinear_angle, &nonlinear_prop);
        contrib.add_term(1, 2, 3, &linear_angle, &linear_prop);
        let pos = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 0.0];
        let mut grad = vec![0.0; pos.len()];

        contrib.get_grad(&pos, &mut grad);

        let mut expected = source_angle_bend_grad(&pos, [0, 1, 2], 80.0, 0.636, false);
        let second = source_angle_bend_grad(&pos, [1, 2, 3], 180.0, 1.006, true);
        for (expected, second) in expected.iter_mut().zip(second.iter()) {
            *expected += second;
        }
        assert_slice_close(&grad, &expected);
    }

    #[test]
    #[should_panic(expected = "not initialized")]
    fn mmff_anglebendcontrib_get_grad_requires_initialized_force_field_for_terms() {
        let force_field = force_field_with_positions(3);
        let mut contrib = AngleBendContrib::new(&force_field);
        let angle_params = source_angle_params(0.636, 100.0);
        let central_atom_prop = source_central_atom_prop(0);
        contrib.add_term(0, 1, 2, &angle_params, &central_atom_prop);
        let mut grad = vec![0.0; 9];

        contrib.get_grad(&[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0], &mut grad);
    }

    #[test]
    fn mmff_anglebendcontrib_calc_angle_bend_grad_accumulates_source_coordinate_paths() {
        let r = [
            ForceFieldVec3::new(0.6, 0.8, -0.1),
            ForceFieldVec3::new(-0.2, 0.5, 0.84),
        ];
        let dist = [1.7, 2.3];
        let mut grad = vec![0.0; 9];

        calc_angle_bend_grad(&r, &dist, &mut grad, [0, 1, 2], 0.75, 0.25, 0.9);

        let dcos_ds = [
            1.0 / dist[0] * (r[1].x - 0.25 * r[0].x),
            1.0 / dist[0] * (r[1].y - 0.25 * r[0].y),
            1.0 / dist[0] * (r[1].z - 0.25 * r[0].z),
            1.0 / dist[1] * (r[0].x - 0.25 * r[1].x),
            1.0 / dist[1] * (r[0].y - 0.25 * r[1].y),
            1.0 / dist[1] * (r[0].z - 0.25 * r[1].z),
        ];
        let expected = [
            0.75 * dcos_ds[0] / -0.9,
            0.75 * dcos_ds[1] / -0.9,
            0.75 * dcos_ds[2] / -0.9,
            0.75 * (-dcos_ds[0] - dcos_ds[3]) / -0.9,
            0.75 * (-dcos_ds[1] - dcos_ds[4]) / -0.9,
            0.75 * (-dcos_ds[2] - dcos_ds[5]) / -0.9,
            0.75 * dcos_ds[3] / -0.9,
            0.75 * dcos_ds[4] / -0.9,
            0.75 * dcos_ds[5] / -0.9,
        ];

        assert_eq!(grad, expected);
    }

    #[test]
    fn mmff_anglebendcontrib_calc_angle_bend_grad_adds_to_existing_gradient() {
        let r = [
            ForceFieldVec3::new(0.6, 0.8, -0.1),
            ForceFieldVec3::new(-0.2, 0.5, 0.84),
        ];
        let dist = [1.7, 2.3];
        let baseline = vec![1.0, -1.0, 0.5, 2.0, -2.0, 1.5, 3.0, -3.0, 2.5];
        let mut grad = baseline.clone();

        calc_angle_bend_grad(&r, &dist, &mut grad, [0, 1, 2], 0.75, 0.25, 0.9);

        let mut delta = vec![0.0; 9];
        calc_angle_bend_grad(&r, &dist, &mut delta, [0, 1, 2], 0.75, 0.25, 0.9);
        let expected = baseline
            .iter()
            .zip(delta.iter())
            .map(|(baseline, delta)| baseline + delta)
            .collect::<Vec<_>>();
        assert_slice_close(&grad, &expected);
    }

    #[test]
    fn mmff_anglebendcontrib_get_grad_uses_singularity_sine_clamp() {
        let force_field = initialized_force_field_with_positions(3);
        let mut contrib = AngleBendContrib::new(&force_field);
        let angle_params = source_angle_params(0.636, 100.0);
        let central_atom_prop = source_central_atom_prop(0);
        contrib.add_term(0, 1, 2, &angle_params, &central_atom_prop);
        let pos = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0];
        let mut grad = vec![0.0; pos.len()];

        contrib.get_grad(&pos, &mut grad);

        assert!(grad.iter().all(|value| value.is_finite()));
    }
}

impl ForceFieldContrib for AngleBendContrib {
    fn copy(&self) -> Box<dyn ForceFieldContrib> {
        Box::new(self.clone())
    }

    fn set_force_field(&mut self, owner: *const ForceField) {
        self.owner = owner;
    }

    fn get_energy(&self, pos: &[f64]) -> f64 {
        AngleBendContrib::get_energy(self, pos)
    }

    fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
        AngleBendContrib::get_grad(self, pos, grad);
    }
}
