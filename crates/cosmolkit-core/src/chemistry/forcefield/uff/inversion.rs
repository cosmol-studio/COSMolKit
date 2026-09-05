//! Source-backed RDKit UFF inversion contribution.

use crate::chemistry::forcefield::core::{ForceField, ForceFieldContrib, ForceFieldVec3};

use super::params::{clip_to_one, is_double_zero};
use super::utils::calc_inversion_coefficients;

#[derive(Clone, Debug)]
pub struct InversionContrib {
    owner: *const ForceField,
    at1_idx: usize,
    at2_idx: usize,
    at3_idx: usize,
    at4_idx: usize,
    force_constant: f64,
    c0: f64,
    c1: f64,
    c2: f64,
}

#[derive(Clone, Debug)]
pub struct InversionContribsParams {
    idx1: usize,
    idx2: usize,
    idx3: usize,
    idx4: usize,
    at2_atomic_num: i32,
    is_c_bound_to_o: bool,
    c0: f64,
    c1: f64,
    c2: f64,
    force_constant: f64,
}

#[derive(Clone, Debug, Default)]
pub struct InversionContribs {
    owner: *const ForceField,
    contribs: Vec<InversionContribsParams>,
}

impl InversionContribsParams {
    #[must_use]
    pub fn idx1(&self) -> usize {
        self.idx1
    }

    #[must_use]
    pub fn idx2(&self) -> usize {
        self.idx2
    }

    #[must_use]
    pub fn idx3(&self) -> usize {
        self.idx3
    }

    #[must_use]
    pub fn idx4(&self) -> usize {
        self.idx4
    }

    #[must_use]
    pub fn at2_atomic_num(&self) -> i32 {
        self.at2_atomic_num
    }

    #[must_use]
    pub fn is_c_bound_to_o(&self) -> bool {
        self.is_c_bound_to_o
    }

    #[must_use]
    pub fn c0(&self) -> f64 {
        self.c0
    }

    #[must_use]
    pub fn c1(&self) -> f64 {
        self.c1
    }

    #[must_use]
    pub fn c2(&self) -> f64 {
        self.c2
    }

    #[must_use]
    pub fn force_constant(&self) -> f64 {
        self.force_constant
    }
}

impl InversionContrib {
    pub fn new(
        owner: &ForceField,
        idx1: usize,
        idx2: usize,
        idx3: usize,
        idx4: usize,
        at2_atomic_num: i32,
        is_c_bound_to_o: bool,
        oob_force_scaling_factor: f64,
    ) -> Self {
        // RDKit✔️✔️: InversionContrib::InversionContrib(ForceField *owner, unsigned int idx1,
        // RDKit✔️✔️:                                    unsigned int idx2, unsigned int idx3,
        // RDKit✔️✔️:                                    unsigned int idx4, int at2AtomicNum,
        // RDKit✔️✔️:                                    bool isCBoundToO,
        // RDKit✔️✔️:                                    double oobForceScalingFactor) {
        // RDKit✔️✔️:   PRECONDITION(owner, "bad owner");
        // Rust references reproduce RDKit's non-null owner precondition.
        // RDKit✔️✔️:   URANGE_CHECK(idx1, owner->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx2, owner->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx3, owner->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx4, owner->positions().size());
        assert!(idx1 < owner.positions().len());
        assert!(idx2 < owner.positions().len());
        assert!(idx3 < owner.positions().len());
        assert!(idx4 < owner.positions().len());

        // RDKit✔️✔️:   dp_forceField = owner;
        // RDKit✔️✔️:   d_at1Idx = idx1;
        // RDKit✔️✔️:   d_at2Idx = idx2;
        // RDKit✔️✔️:   d_at3Idx = idx3;
        // RDKit✔️✔️:   d_at4Idx = idx4;
        let owner_ptr = owner as *const ForceField;

        // RDKit✔️✔️:   auto invCoeffForceCon = Utils::calcInversionCoefficientsAndForceConstant(
        // RDKit✔️✔️:       at2AtomicNum, isCBoundToO);
        let (force_constant, c0, c1, c2) =
            calc_inversion_coefficients(at2_atomic_num, is_c_bound_to_o);
        // RDKit✔️✔️:   d_forceConstant = oobForceScalingFactor * std::get<0>(invCoeffForceCon);
        // RDKit✔️✔️:   d_C0 = std::get<1>(invCoeffForceCon);
        // RDKit✔️✔️:   d_C1 = std::get<2>(invCoeffForceCon);
        // RDKit✔️✔️:   d_C2 = std::get<3>(invCoeffForceCon);
        // RDKit✔️✔️: }
        Self {
            owner: owner_ptr,
            at1_idx: idx1,
            at2_idx: idx2,
            at3_idx: idx3,
            at4_idx: idx4,
            force_constant: oob_force_scaling_factor * force_constant,
            c0,
            c1,
            c2,
        }
    }

    #[must_use]
    pub fn owner(&self) -> *const ForceField {
        self.owner
    }

    #[must_use]
    pub fn at1_idx(&self) -> usize {
        self.at1_idx
    }

    #[must_use]
    pub fn at2_idx(&self) -> usize {
        self.at2_idx
    }

    #[must_use]
    pub fn at3_idx(&self) -> usize {
        self.at3_idx
    }

    #[must_use]
    pub fn at4_idx(&self) -> usize {
        self.at4_idx
    }

    #[must_use]
    pub fn force_constant(&self) -> f64 {
        self.force_constant
    }

    #[must_use]
    pub fn c0(&self) -> f64 {
        self.c0
    }

    #[must_use]
    pub fn c1(&self) -> f64 {
        self.c1
    }

    #[must_use]
    pub fn c2(&self) -> f64 {
        self.c2
    }

    #[must_use]
    pub fn get_energy(&self, pos: &[f64]) -> f64 {
        // RDKit✔️✔️: double InversionContrib::getEnergy(double *pos) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        assert!(!self.owner.is_null(), "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        let required_len = 3 * self
            .at1_idx
            .max(self.at2_idx)
            .max(self.at3_idx)
            .max(self.at4_idx)
            + 3;
        assert!(pos.len() >= required_len, "bad vector");

        // RDKit✔️✔️:   RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
        // RDKit✔️✔️:                      pos[3 * d_at1Idx + 2]);
        let p1 = ForceFieldVec3::new(
            pos[3 * self.at1_idx],
            pos[3 * self.at1_idx + 1],
            pos[3 * self.at1_idx + 2],
        );
        // RDKit✔️✔️:   RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
        // RDKit✔️✔️:                      pos[3 * d_at2Idx + 2]);
        let p2 = ForceFieldVec3::new(
            pos[3 * self.at2_idx],
            pos[3 * self.at2_idx + 1],
            pos[3 * self.at2_idx + 2],
        );
        // RDKit✔️✔️:   RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
        // RDKit✔️✔️:                      pos[3 * d_at3Idx + 2]);
        let p3 = ForceFieldVec3::new(
            pos[3 * self.at3_idx],
            pos[3 * self.at3_idx + 1],
            pos[3 * self.at3_idx + 2],
        );
        // RDKit✔️✔️:   RDGeom::Point3D p4(pos[3 * d_at4Idx], pos[3 * d_at4Idx + 1],
        // RDKit✔️✔️:                      pos[3 * d_at4Idx + 2]);
        let p4 = ForceFieldVec3::new(
            pos[3 * self.at4_idx],
            pos[3 * self.at4_idx + 1],
            pos[3 * self.at4_idx + 2],
        );

        // RDKit✔️✔️:   double cosY = Utils::calculateCosY(p1, p2, p3, p4);
        let cos_y = calculate_cos_y(p1, p2, p3, p4);
        // RDKit✔️✔️:   double sinYSq = 1.0 - cosY * cosY;
        let sin_y_sq = 1.0 - cos_y * cos_y;
        // RDKit✔️✔️:   double sinY = ((sinYSq > 0.0) ? sqrt(sinYSq) : 0.0);
        let sin_y = if sin_y_sq > 0.0 { sin_y_sq.sqrt() } else { 0.0 };
        // RDKit✔️✔️:   // cos(2 * W) = 2 * cos(W) * cos(W) - 1 = 2 * sin(W) * sin(W) - 1
        // RDKit✔️✔️:   double cos2W = 2.0 * sinY * sinY - 1.0;
        let cos2_w = 2.0 * sin_y * sin_y - 1.0;
        // RDKit✔️✔️:   double res = d_forceConstant * (d_C0 + d_C1 * sinY + d_C2 * cos2W);
        let res = self.force_constant * (self.c0 + self.c1 * sin_y + self.c2 * cos2_w);
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        res
    }

    pub fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
        // RDKit✔️✔️: void InversionContrib::getGrad(double *pos, double *grad) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        assert!(!self.owner.is_null(), "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:   PRECONDITION(grad, "bad vector");
        let required_len = 3 * self
            .at1_idx
            .max(self.at2_idx)
            .max(self.at3_idx)
            .max(self.at4_idx)
            + 3;
        assert!(pos.len() >= required_len, "bad vector");
        assert!(grad.len() >= required_len, "bad vector");

        // RDKit✔️✔️:   RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
        // RDKit✔️✔️:                      pos[3 * d_at1Idx + 2]);
        let p1 = ForceFieldVec3::new(
            pos[3 * self.at1_idx],
            pos[3 * self.at1_idx + 1],
            pos[3 * self.at1_idx + 2],
        );
        // RDKit✔️✔️:   RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
        // RDKit✔️✔️:                      pos[3 * d_at2Idx + 2]);
        let p2 = ForceFieldVec3::new(
            pos[3 * self.at2_idx],
            pos[3 * self.at2_idx + 1],
            pos[3 * self.at2_idx + 2],
        );
        // RDKit✔️✔️:   RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
        // RDKit✔️✔️:                      pos[3 * d_at3Idx + 2]);
        let p3 = ForceFieldVec3::new(
            pos[3 * self.at3_idx],
            pos[3 * self.at3_idx + 1],
            pos[3 * self.at3_idx + 2],
        );
        // RDKit✔️✔️:   RDGeom::Point3D p4(pos[3 * d_at4Idx], pos[3 * d_at4Idx + 1],
        // RDKit✔️✔️:                      pos[3 * d_at4Idx + 2]);
        let p4 = ForceFieldVec3::new(
            pos[3 * self.at4_idx],
            pos[3 * self.at4_idx + 1],
            pos[3 * self.at4_idx + 2],
        );
        // RDKit✔️✔️:   double *g1 = &(grad[3 * d_at1Idx]);
        // RDKit✔️✔️:   double *g2 = &(grad[3 * d_at2Idx]);
        // RDKit✔️✔️:   double *g3 = &(grad[3 * d_at3Idx]);
        // RDKit✔️✔️:   double *g4 = &(grad[3 * d_at4Idx]);
        let g1 = 3 * self.at1_idx;
        let g2 = 3 * self.at2_idx;
        let g3 = 3 * self.at3_idx;
        let g4 = 3 * self.at4_idx;

        // RDKit✔️✔️:   RDGeom::Point3D rJI = p1 - p2;
        // RDKit✔️✔️:   RDGeom::Point3D rJK = p3 - p2;
        // RDKit✔️✔️:   RDGeom::Point3D rJL = p4 - p2;
        let mut r_ji = p1 - p2;
        let mut r_jk = p3 - p2;
        let mut r_jl = p4 - p2;
        // RDKit✔️✔️:   double dJI = rJI.length();
        // RDKit✔️✔️:   double dJK = rJK.length();
        // RDKit✔️✔️:   double dJL = rJL.length();
        let d_ji = r_ji.length();
        let d_jk = r_jk.length();
        let d_jl = r_jl.length();
        // RDKit✔️✔️:   if (isDoubleZero(dJI) || isDoubleZero(dJK) || isDoubleZero(dJL)) {
        // RDKit✔️✔️:     return;
        // RDKit✔️✔️:   }
        if is_double_zero(d_ji) || is_double_zero(d_jk) || is_double_zero(d_jl) {
            return;
        }
        // RDKit✔️✔️:   rJI /= dJI;
        // RDKit✔️✔️:   rJK /= dJK;
        // RDKit✔️✔️:   rJL /= dJL;
        r_ji /= d_ji;
        r_jk /= d_jk;
        r_jl /= d_jl;

        // RDKit✔️✔️:   RDGeom::Point3D n = (-rJI).crossProduct(rJK);
        // RDKit✔️✔️:   n /= n.length();
        let mut n = (-r_ji).cross_product(r_jk);
        n /= n.length();
        // RDKit✔️✔️:   double cosY = n.dotProduct(rJL);
        // RDKit✔️✔️:   clipToOne(cosY);
        let mut cos_y = n.dot_product(r_jl);
        clip_to_one(&mut cos_y);
        // RDKit✔️✔️:   double sinYSq = 1.0 - cosY * cosY;
        // RDKit✔️✔️:   double sinY = std::max(sqrt(sinYSq), 1.0e-8);
        let sin_y_sq = 1.0 - cos_y * cos_y;
        let sin_y = sin_y_sq.sqrt().max(1.0e-8);
        // RDKit✔️✔️:   double cosTheta = rJI.dotProduct(rJK);
        // RDKit✔️✔️:   clipToOne(cosTheta);
        let mut cos_theta = r_ji.dot_product(r_jk);
        clip_to_one(&mut cos_theta);
        // RDKit✔️✔️:   double sinThetaSq = 1.0 - cosTheta * cosTheta;
        // RDKit✔️✔️:   double sinTheta = std::max(sqrt(sinThetaSq), 1.0e-8);
        let sin_theta_sq = 1.0 - cos_theta * cos_theta;
        let sin_theta = sin_theta_sq.sqrt().max(1.0e-8);
        // RDKit✔️✔️:   // sin(2 * W) = 2 * sin(W) * cos(W) = 2 * cos(Y) * sin(Y)
        // RDKit✔️✔️:   double dE_dW = -d_forceConstant * (d_C1 * cosY - 4.0 * d_C2 * cosY * sinY);
        let de_dw = -self.force_constant * (self.c1 * cos_y - 4.0 * self.c2 * cos_y * sin_y);
        // RDKit✔️✔️:   RDGeom::Point3D t1 = rJL.crossProduct(rJK);
        // RDKit✔️✔️:   RDGeom::Point3D t2 = rJI.crossProduct(rJL);
        // RDKit✔️✔️:   RDGeom::Point3D t3 = rJK.crossProduct(rJI);
        let t1 = r_jl.cross_product(r_jk);
        let t2 = r_ji.cross_product(r_jl);
        let t3 = r_jk.cross_product(r_ji);
        // RDKit✔️✔️:   double term1 = sinY * sinTheta;
        // RDKit✔️✔️:   double term2 = cosY / (sinY * sinThetaSq);
        let term1 = sin_y * sin_theta;
        let term2 = cos_y / (sin_y * sin_theta_sq);
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
        // RDKit✔️✔️:   double tg4[3] = {(t3.x / term1 - rJL.x * cosY / sinY) / dJL,
        // RDKit✔️✔️:                    (t3.y / term1 - rJL.y * cosY / sinY) / dJL,
        // RDKit✔️✔️:                    (t3.z / term1 - rJL.z * cosY / sinY) / dJL};
        let tg4 = [
            (t3.x / term1 - r_jl.x * cos_y / sin_y) / d_jl,
            (t3.y / term1 - r_jl.y * cos_y / sin_y) / d_jl,
            (t3.z / term1 - r_jl.z * cos_y / sin_y) / d_jl,
        ];
        // RDKit✔️✔️:   for (unsigned int i = 0; i < 3; ++i) {
        // RDKit✔️✔️:     g1[i] += dE_dW * tg1[i];
        // RDKit✔️✔️:     g2[i] += -dE_dW * (tg1[i] + tg3[i] + tg4[i]);
        // RDKit✔️✔️:     g3[i] += dE_dW * tg3[i];
        // RDKit✔️✔️:     g4[i] += dE_dW * tg4[i];
        // RDKit✔️✔️:   }
        for i in 0..3 {
            grad[g1 + i] += de_dw * tg1[i];
            grad[g2 + i] += -de_dw * (tg1[i] + tg3[i] + tg4[i]);
            grad[g3 + i] += de_dw * tg3[i];
            grad[g4 + i] += de_dw * tg4[i];
        }
        // RDKit✔️✔️: }
    }
}

impl ForceFieldContrib for InversionContrib {
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

impl ForceFieldContrib for InversionContribs {
    fn copy(&self) -> Box<dyn ForceFieldContrib> {
        // RDKit✔️✔️: InversionContribs *copy() const override {
        // RDKit✔️✔️:   return new InversionContribs(*this);
        // RDKit✔️✔️: }
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

impl InversionContribs {
    #[must_use]
    pub fn new(owner: &ForceField) -> Self {
        // RDKit✔️✔️: InversionContribs::InversionContribs(ForceField *owner) {
        // RDKit✔️✔️:   PRECONDITION(owner, "bad owner");
        // Rust references reproduce RDKit's non-null owner precondition.
        // RDKit✔️✔️:   dp_forceField = owner;
        // RDKit✔️✔️: }
        Self {
            owner: owner as *const ForceField,
            contribs: Vec::new(),
        }
    }

    #[must_use]
    pub fn owner(&self) -> *const ForceField {
        self.owner
    }

    #[must_use]
    pub fn contribs(&self) -> &[InversionContribsParams] {
        &self.contribs
    }

    #[must_use]
    pub fn empty(&self) -> bool {
        // RDKit✔️✔️: bool empty() const { return d_contribs.empty(); }
        self.contribs.is_empty()
    }

    #[must_use]
    pub fn size(&self) -> usize {
        // RDKit✔️✔️: unsigned int size() const { return d_contribs.size(); }
        self.contribs.len()
    }

    pub fn add_contrib(
        &mut self,
        idx1: usize,
        idx2: usize,
        idx3: usize,
        idx4: usize,
        at2_atomic_num: i32,
        is_c_bound_to_o: bool,
        oob_force_scaling_factor: f64,
    ) {
        // RDKit✔️✔️: void InversionContribs::addContrib(unsigned int idx1, unsigned int idx2,
        // RDKit✔️✔️:                                    unsigned int idx3, unsigned int idx4,
        // RDKit✔️✔️:                                    int at2AtomicNum, bool isCBoundToO,
        // RDKit✔️✔️:                                    double oobForceScalingFactor) {
        assert!(!self.owner.is_null(), "no owner");
        // RDKit✔️✔️:   URANGE_CHECK(idx1, dp_forceField->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx2, dp_forceField->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx3, dp_forceField->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx4, dp_forceField->positions().size());
        let owner = unsafe { &*self.owner };
        assert!(idx1 < owner.positions().len());
        assert!(idx2 < owner.positions().len());
        assert!(idx3 < owner.positions().len());
        assert!(idx4 < owner.positions().len());
        // RDKit✔️✔️:   auto invCoeffForceCon = Utils::calcInversionCoefficientsAndForceConstant(
        // RDKit✔️✔️:       at2AtomicNum, isCBoundToO);
        let (force_constant, c0, c1, c2) =
            calc_inversion_coefficients(at2_atomic_num, is_c_bound_to_o);

        // RDKit✔️✔️:   d_contribs.emplace_back(
        // RDKit✔️✔️:       idx1, idx2, idx3, idx4, at2AtomicNum, isCBoundToO,
        // RDKit✔️✔️:       std::get<1>(invCoeffForceCon), std::get<2>(invCoeffForceCon),
        // RDKit✔️✔️:       std::get<3>(invCoeffForceCon),
        // RDKit✔️✔️:       std::get<0>(invCoeffForceCon) * oobForceScalingFactor);
        // RDKit✔️✔️: }
        self.contribs.push(InversionContribsParams {
            idx1,
            idx2,
            idx3,
            idx4,
            at2_atomic_num,
            is_c_bound_to_o,
            c0,
            c1,
            c2,
            force_constant: force_constant * oob_force_scaling_factor,
        });
    }

    #[must_use]
    pub fn get_energy(&self, pos: &[f64]) -> f64 {
        // RDKit✔️✔️: double InversionContribs::getEnergy(double *pos) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        assert!(!self.owner.is_null(), "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        assert!(!pos.is_empty(), "bad vector");
        // RDKit✔️✔️:   double accum = 0;
        let mut accum = 0.0;
        // RDKit✔️✔️:   for (const auto &contrib : d_contribs) {
        for contrib in &self.contribs {
            let required_len = 3 * contrib
                .idx1
                .max(contrib.idx2)
                .max(contrib.idx3)
                .max(contrib.idx4)
                + 3;
            assert!(pos.len() >= required_len, "bad vector");
            // RDKit✔️✔️:     const RDGeom::Point3D p1(pos[3 * contrib.idx1], pos[3 * contrib.idx1 + 1],
            // RDKit✔️✔️:                              pos[3 * contrib.idx1 + 2]);
            let p1 = ForceFieldVec3::new(
                pos[3 * contrib.idx1],
                pos[3 * contrib.idx1 + 1],
                pos[3 * contrib.idx1 + 2],
            );
            // RDKit✔️✔️:     const RDGeom::Point3D p2(pos[3 * contrib.idx2], pos[3 * contrib.idx2 + 1],
            // RDKit✔️✔️:                              pos[3 * contrib.idx2 + 2]);
            let p2 = ForceFieldVec3::new(
                pos[3 * contrib.idx2],
                pos[3 * contrib.idx2 + 1],
                pos[3 * contrib.idx2 + 2],
            );
            // RDKit✔️✔️:     const RDGeom::Point3D p3(pos[3 * contrib.idx3], pos[3 * contrib.idx3 + 1],
            // RDKit✔️✔️:                              pos[3 * contrib.idx3 + 2]);
            let p3 = ForceFieldVec3::new(
                pos[3 * contrib.idx3],
                pos[3 * contrib.idx3 + 1],
                pos[3 * contrib.idx3 + 2],
            );
            // RDKit✔️✔️:     const RDGeom::Point3D p4(pos[3 * contrib.idx4], pos[3 * contrib.idx4 + 1],
            // RDKit✔️✔️:                              pos[3 * contrib.idx4 + 2]);
            let p4 = ForceFieldVec3::new(
                pos[3 * contrib.idx4],
                pos[3 * contrib.idx4 + 1],
                pos[3 * contrib.idx4 + 2],
            );
            // RDKit✔️✔️:     const double cosY = Utils::calculateCosY(p1, p2, p3, p4);
            let cos_y = calculate_cos_y(p1, p2, p3, p4);
            // RDKit✔️✔️:     const double sinYSq = 1.0 - cosY * cosY;
            let sin_y_sq = 1.0 - cos_y * cos_y;
            // RDKit✔️✔️:     const double sinY = ((sinYSq > 0.0) ? sqrt(sinYSq) : 0.0);
            let sin_y = if sin_y_sq > 0.0 { sin_y_sq.sqrt() } else { 0.0 };
            // RDKit✔️✔️:     // cos(2 * W) = 2 * cos(W) * cos(W) - 1 = 2 * sin(W) * sin(W) - 1
            // RDKit✔️✔️:     const double cos2W = 2.0 * sinY * sinY - 1.0;
            let cos2_w = 2.0 * sin_y * sin_y - 1.0;
            // RDKit✔️✔️:     accum += contrib.forceConstant *
            // RDKit✔️✔️:              (contrib.C0 + contrib.C1 * sinY + contrib.C2 * cos2W);
            accum +=
                contrib.force_constant * (contrib.c0 + contrib.c1 * sin_y + contrib.c2 * cos2_w);
            // RDKit✔️✔️:   }
        }
        // RDKit✔️✔️:   return accum;
        // RDKit✔️✔️: }
        accum
    }

    pub fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
        // RDKit✔️✔️: void InversionContribs::getGrad(double *pos, double *grad) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        assert!(!self.owner.is_null(), "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:   PRECONDITION(grad, "bad vector");
        assert!(!pos.is_empty(), "bad vector");
        assert!(!grad.is_empty(), "bad vector");
        // RDKit✔️✔️:   for (const auto &contrib : d_contribs) {
        for contrib in &self.contribs {
            let required_len = 3 * contrib
                .idx1
                .max(contrib.idx2)
                .max(contrib.idx3)
                .max(contrib.idx4)
                + 3;
            assert!(pos.len() >= required_len, "bad vector");
            assert!(grad.len() >= required_len, "bad vector");
            // RDKit✔️✔️:     const RDGeom::Point3D p1(pos[3 * contrib.idx1], pos[3 * contrib.idx1 + 1],
            // RDKit✔️✔️:                              pos[3 * contrib.idx1 + 2]);
            let p1 = ForceFieldVec3::new(
                pos[3 * contrib.idx1],
                pos[3 * contrib.idx1 + 1],
                pos[3 * contrib.idx1 + 2],
            );
            // RDKit✔️✔️:     const RDGeom::Point3D p2(pos[3 * contrib.idx2], pos[3 * contrib.idx2 + 1],
            // RDKit✔️✔️:                              pos[3 * contrib.idx2 + 2]);
            let p2 = ForceFieldVec3::new(
                pos[3 * contrib.idx2],
                pos[3 * contrib.idx2 + 1],
                pos[3 * contrib.idx2 + 2],
            );
            // RDKit✔️✔️:     const RDGeom::Point3D p3(pos[3 * contrib.idx3], pos[3 * contrib.idx3 + 1],
            // RDKit✔️✔️:                              pos[3 * contrib.idx3 + 2]);
            let p3 = ForceFieldVec3::new(
                pos[3 * contrib.idx3],
                pos[3 * contrib.idx3 + 1],
                pos[3 * contrib.idx3 + 2],
            );
            // RDKit✔️✔️:     const RDGeom::Point3D p4(pos[3 * contrib.idx4], pos[3 * contrib.idx4 + 1],
            // RDKit✔️✔️:                              pos[3 * contrib.idx4 + 2]);
            let p4 = ForceFieldVec3::new(
                pos[3 * contrib.idx4],
                pos[3 * contrib.idx4 + 1],
                pos[3 * contrib.idx4 + 2],
            );
            // RDKit✔️✔️:     double *g1 = &(grad[3 * contrib.idx1]);
            // RDKit✔️✔️:     double *g2 = &(grad[3 * contrib.idx2]);
            // RDKit✔️✔️:     double *g3 = &(grad[3 * contrib.idx3]);
            // RDKit✔️✔️:     double *g4 = &(grad[3 * contrib.idx4]);
            let g1 = 3 * contrib.idx1;
            let g2 = 3 * contrib.idx2;
            let g3 = 3 * contrib.idx3;
            let g4 = 3 * contrib.idx4;
            // RDKit✔️✔️:     RDGeom::Point3D rJI = p1 - p2;
            // RDKit✔️✔️:     RDGeom::Point3D rJK = p3 - p2;
            // RDKit✔️✔️:     RDGeom::Point3D rJL = p4 - p2;
            let mut r_ji = p1 - p2;
            let mut r_jk = p3 - p2;
            let mut r_jl = p4 - p2;
            // RDKit✔️✔️:     const double dJI = rJI.length();
            // RDKit✔️✔️:     const double dJK = rJK.length();
            // RDKit✔️✔️:     const double dJL = rJL.length();
            let d_ji = r_ji.length();
            let d_jk = r_jk.length();
            let d_jl = r_jl.length();
            // RDKit✔️✔️:     if (isDoubleZero(dJI) || isDoubleZero(dJK) || isDoubleZero(dJL)) {
            // RDKit✔️✔️:       return;
            // RDKit✔️✔️:     }
            if is_double_zero(d_ji) || is_double_zero(d_jk) || is_double_zero(d_jl) {
                return;
            }
            // RDKit✔️✔️:     rJI.normalize();
            // RDKit✔️✔️:     rJK.normalize();
            // RDKit✔️✔️:     rJL.normalize();
            r_ji /= d_ji;
            r_jk /= d_jk;
            r_jl /= d_jl;

            // RDKit✔️✔️:     RDGeom::Point3D n = (-rJI).crossProduct(rJK);
            // RDKit✔️✔️:     n.normalize();
            let mut n = (-r_ji).cross_product(r_jk);
            n /= n.length();
            // RDKit✔️✔️:     double cosY = n.dotProduct(rJL);
            // RDKit✔️✔️:     cosY = std::clamp(cosY, -1.0, 1.0);
            let cos_y = n.dot_product(r_jl).clamp(-1.0, 1.0);
            // RDKit✔️✔️:     const double sinYSq = 1.0 - cosY * cosY;
            // RDKit✔️✔️:     const double sinY = std::max(sqrt(sinYSq), 1.0e-8);
            let sin_y_sq = 1.0 - cos_y * cos_y;
            let sin_y = sin_y_sq.sqrt().max(1.0e-8);
            // RDKit✔️✔️:     double cosTheta = rJI.dotProduct(rJK);
            // RDKit✔️✔️:     cosTheta = std::clamp(cosTheta, -1.0, 1.0);
            let cos_theta = r_ji.dot_product(r_jk).clamp(-1.0, 1.0);
            // RDKit✔️✔️:     const double sinThetaSq = 1.0 - cosTheta * cosTheta;
            // RDKit✔️✔️:     const double sinTheta = std::max(sqrt(sinThetaSq), 1.0e-8);
            let sin_theta_sq = 1.0 - cos_theta * cos_theta;
            let sin_theta = sin_theta_sq.sqrt().max(1.0e-8);
            // RDKit✔️✔️:     // sin(2 * W) = 2 * sin(W) * cos(W) = 2 * cos(Y) * sin(Y)
            // RDKit✔️✔️:     const double dE_dW = -contrib.forceConstant *
            // RDKit✔️✔️:                          (contrib.C1 * cosY - 4.0 * contrib.C2 * cosY * sinY);
            let de_dw =
                -contrib.force_constant * (contrib.c1 * cos_y - 4.0 * contrib.c2 * cos_y * sin_y);
            // RDKit✔️✔️:     const RDGeom::Point3D t1 = rJL.crossProduct(rJK);
            // RDKit✔️✔️:     const RDGeom::Point3D t2 = rJI.crossProduct(rJL);
            // RDKit✔️✔️:     const RDGeom::Point3D t3 = rJK.crossProduct(rJI);
            let t1 = r_jl.cross_product(r_jk);
            let t2 = r_ji.cross_product(r_jl);
            let t3 = r_jk.cross_product(r_ji);
            // RDKit✔️✔️:     const double term1 = sinY * sinTheta;
            // RDKit✔️✔️:     const double term2 = cosY / (sinY * sinThetaSq);
            let term1 = sin_y * sin_theta;
            let term2 = cos_y / (sin_y * sin_theta_sq);
            // RDKit✔️✔️:     const double tg1[3] = {
            // RDKit✔️✔️:         (t1.x / term1 - (rJI.x - rJK.x * cosTheta) * term2) / dJI,
            // RDKit✔️✔️:         (t1.y / term1 - (rJI.y - rJK.y * cosTheta) * term2) / dJI,
            // RDKit✔️✔️:         (t1.z / term1 - (rJI.z - rJK.z * cosTheta) * term2) / dJI};
            let tg1 = [
                (t1.x / term1 - (r_ji.x - r_jk.x * cos_theta) * term2) / d_ji,
                (t1.y / term1 - (r_ji.y - r_jk.y * cos_theta) * term2) / d_ji,
                (t1.z / term1 - (r_ji.z - r_jk.z * cos_theta) * term2) / d_ji,
            ];
            // RDKit✔️✔️:     const double tg3[3] = {
            // RDKit✔️✔️:         (t2.x / term1 - (rJK.x - rJI.x * cosTheta) * term2) / dJK,
            // RDKit✔️✔️:         (t2.y / term1 - (rJK.y - rJI.y * cosTheta) * term2) / dJK,
            // RDKit✔️✔️:         (t2.z / term1 - (rJK.z - rJI.z * cosTheta) * term2) / dJK};
            let tg3 = [
                (t2.x / term1 - (r_jk.x - r_ji.x * cos_theta) * term2) / d_jk,
                (t2.y / term1 - (r_jk.y - r_ji.y * cos_theta) * term2) / d_jk,
                (t2.z / term1 - (r_jk.z - r_ji.z * cos_theta) * term2) / d_jk,
            ];
            // RDKit✔️✔️:     const double tg4[3] = {(t3.x / term1 - rJL.x * cosY / sinY) / dJL,
            // RDKit✔️✔️:                            (t3.y / term1 - rJL.y * cosY / sinY) / dJL,
            // RDKit✔️✔️:                            (t3.z / term1 - rJL.z * cosY / sinY) / dJL};
            let tg4 = [
                (t3.x / term1 - r_jl.x * cos_y / sin_y) / d_jl,
                (t3.y / term1 - r_jl.y * cos_y / sin_y) / d_jl,
                (t3.z / term1 - r_jl.z * cos_y / sin_y) / d_jl,
            ];
            // RDKit✔️✔️:     for (unsigned int i = 0; i < 3; ++i) {
            // RDKit✔️✔️:       g1[i] += dE_dW * tg1[i];
            // RDKit✔️✔️:       g2[i] += -dE_dW * (tg1[i] + tg3[i] + tg4[i]);
            // RDKit✔️✔️:       g3[i] += dE_dW * tg3[i];
            // RDKit✔️✔️:       g4[i] += dE_dW * tg4[i];
            // RDKit✔️✔️:     }
            for i in 0..3 {
                grad[g1 + i] += de_dw * tg1[i];
                grad[g2 + i] += -de_dw * (tg1[i] + tg3[i] + tg4[i]);
                grad[g3 + i] += de_dw * tg3[i];
                grad[g4 + i] += de_dw * tg4[i];
            }
            // RDKit✔️✔️:   }
        }
        // RDKit✔️✔️: }
    }
}

fn calculate_cos_y(
    i_point: ForceFieldVec3,
    j_point: ForceFieldVec3,
    k_point: ForceFieldVec3,
    l_point: ForceFieldVec3,
) -> f64 {
    // RDKit✔️✔️: double calculateCosY(const RDGeom::Point3D &iPoint,
    // RDKit✔️✔️:                      const RDGeom::Point3D &jPoint,
    // RDKit✔️✔️:                      const RDGeom::Point3D &kPoint,
    // RDKit✔️✔️:                      const RDGeom::Point3D &lPoint) {
    // RDKit✔️✔️:   constexpr double zeroTol = 1.0e-16;
    const ZERO_TOL: f64 = 1.0e-16;
    // RDKit✔️✔️:   RDGeom::Point3D rJI = iPoint - jPoint;
    // RDKit✔️✔️:   RDGeom::Point3D rJK = kPoint - jPoint;
    // RDKit✔️✔️:   RDGeom::Point3D rJL = lPoint - jPoint;
    let r_ji = i_point - j_point;
    let r_jk = k_point - j_point;
    let r_jl = l_point - j_point;
    // RDKit✔️✔️:   auto l2JI = rJI.lengthSq();
    // RDKit✔️✔️:   auto l2JK = rJK.lengthSq();
    // RDKit✔️✔️:   auto l2JL = rJL.lengthSq();
    let l2_ji = r_ji.length_sq();
    let l2_jk = r_jk.length_sq();
    let l2_jl = r_jl.length_sq();
    // RDKit✔️✔️:   if (l2JI < zeroTol || l2JK < zeroTol || l2JL < zeroTol) {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   }
    if l2_ji < ZERO_TOL || l2_jk < ZERO_TOL || l2_jl < ZERO_TOL {
        return 0.0;
    }
    // RDKit✔️✔️:   RDGeom::Point3D n = rJI.crossProduct(rJK);
    // RDKit✔️✔️:   n /= (sqrt(l2JI) * sqrt(l2JK));
    let mut n = r_ji.cross_product(r_jk);
    n /= l2_ji.sqrt() * l2_jk.sqrt();
    // RDKit✔️✔️:   auto l2n = n.lengthSq();
    let l2n = n.length_sq();
    // RDKit✔️✔️:   if (l2n < zeroTol) {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   }
    if l2n < ZERO_TOL {
        return 0.0;
    }
    // RDKit✔️✔️:   return n.dotProduct(rJL) / (sqrt(l2JL) * sqrt(l2n));
    // RDKit✔️✔️: }
    n.dot_product(r_jl) / (l2_jl.sqrt() * l2n.sqrt())
}

#[cfg(test)]
mod tests {
    use crate::chemistry::forcefield::core::{ForceField, ForceFieldVec3};

    use super::super::utils::calc_inversion_coefficients;
    use super::*;

    const EPS: f64 = 1.0e-12;

    fn force_field() -> ForceField {
        let mut ff = ForceField::new(3);
        ff.positions_mut().push(ForceFieldVec3::new(0.0, 0.0, 0.0));
        ff.positions_mut().push(ForceFieldVec3::new(1.0, 0.0, 0.0));
        ff.positions_mut().push(ForceFieldVec3::new(0.0, 1.0, 0.0));
        ff.positions_mut().push(ForceFieldVec3::new(0.0, 0.0, 1.0));
        ff
    }

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < EPS,
            "actual {actual} expected {expected} diff {}",
            (actual - expected).abs()
        );
    }

    fn assert_close_tol(actual: f64, expected: f64, tol: f64) {
        assert!(
            (actual - expected).abs() < tol,
            "actual {actual} expected {expected} diff {} tol {tol}",
            (actual - expected).abs()
        );
    }

    #[test]
    fn uff_inversioncontrib_constructor_stores_owner_indices_and_coefficients() {
        let ff = force_field();
        let (force_constant, c0, c1, c2) = calc_inversion_coefficients(6, false);

        let contrib = InversionContrib::new(&ff, 0, 1, 2, 3, 6, false, 1.0);

        assert_eq!(contrib.owner(), &ff as *const ForceField);
        assert_eq!(contrib.at1_idx(), 0);
        assert_eq!(contrib.at2_idx(), 1);
        assert_eq!(contrib.at3_idx(), 2);
        assert_eq!(contrib.at4_idx(), 3);
        assert_close(contrib.force_constant(), force_constant);
        assert_close(contrib.c0(), c0);
        assert_close(contrib.c1(), c1);
        assert_close(contrib.c2(), c2);
    }

    #[test]
    fn uff_inversioncontrib_constructor_applies_carbonyl_and_oob_scaling() {
        let ff = force_field();
        let scale = 0.5;
        let (force_constant, c0, c1, c2) = calc_inversion_coefficients(6, true);

        let contrib = InversionContrib::new(&ff, 0, 1, 2, 3, 6, true, scale);

        assert_close(contrib.force_constant(), scale * force_constant);
        assert_close(contrib.c0(), c0);
        assert_close(contrib.c1(), c1);
        assert_close(contrib.c2(), c2);
    }

    #[test]
    fn uff_inversioncontrib_constructor_uses_group5_coefficients() {
        let ff = force_field();
        let (force_constant, c0, c1, c2) = calc_inversion_coefficients(15, false);

        let contrib = InversionContrib::new(&ff, 0, 1, 2, 3, 15, false, 2.0);

        assert_close(contrib.force_constant(), 2.0 * force_constant);
        assert_close(contrib.c0(), c0);
        assert_close(contrib.c1(), c1);
        assert_close(contrib.c2(), c2);
    }

    #[test]
    #[should_panic]
    fn uff_inversioncontrib_constructor_rejects_index_out_of_range() {
        let ff = force_field();

        let _ = InversionContrib::new(&ff, 0, 1, 2, 4, 6, false, 1.0);
    }

    #[test]
    fn uff_inversioncontrib_getenergy_matches_source_formula_for_nonzero_geometry() {
        let ff = force_field();
        let contrib = InversionContrib::new(&ff, 0, 1, 2, 3, 6, false, 1.0);
        let pos = [
            0.0, 0.0, 0.0, //
            1.0, 0.0, 0.0, //
            0.0, 1.0, 0.0, //
            0.0, 0.0, 1.0,
        ];

        let sin_y = 0.5_f64.sqrt();
        let expected = contrib.force_constant() * (contrib.c0() + contrib.c1() * sin_y);

        assert_close(contrib.get_energy(&pos), expected);
    }

    #[test]
    fn uff_inversioncontrib_getenergy_uses_zero_sin_branch_when_cos_y_is_one() {
        let ff = force_field();
        let contrib = InversionContrib::new(&ff, 0, 1, 2, 3, 15, false, 1.0);
        let pos = [
            1.0, 0.0, 0.0, //
            0.0, 0.0, 0.0, //
            0.0, 1.0, 0.0, //
            0.0, 0.0, 1.0,
        ];

        let expected = contrib.force_constant() * (contrib.c0() - contrib.c2());

        assert_close(
            calculate_cos_y(
                ForceFieldVec3::new(pos[0], pos[1], pos[2]),
                ForceFieldVec3::new(pos[3], pos[4], pos[5]),
                ForceFieldVec3::new(pos[6], pos[7], pos[8]),
                ForceFieldVec3::new(pos[9], pos[10], pos[11]),
            ),
            1.0,
        );
        assert_close(contrib.get_energy(&pos), expected);
    }

    #[test]
    fn uff_inversioncontrib_calculate_cos_y_returns_zero_for_degenerate_vectors() {
        let origin = ForceFieldVec3::new(0.0, 0.0, 0.0);
        let x = ForceFieldVec3::new(1.0, 0.0, 0.0);
        let y = ForceFieldVec3::new(0.0, 1.0, 0.0);

        assert_close(calculate_cos_y(origin, origin, y, x), 0.0);
        assert_close(calculate_cos_y(x, origin, y, origin), 0.0);
        assert_close(
            calculate_cos_y(x, origin, ForceFieldVec3::new(2.0, 0.0, 0.0), y),
            0.0,
        );
    }

    #[test]
    #[should_panic(expected = "bad vector")]
    fn uff_inversioncontrib_getenergy_rejects_short_position_vector() {
        let ff = force_field();
        let contrib = InversionContrib::new(&ff, 0, 1, 2, 3, 6, false, 1.0);

        let _ = contrib.get_energy(&[0.0, 0.0, 0.0]);
    }

    #[test]
    fn uff_inversioncontrib_getgrad_matches_getenergy_finite_difference_for_carbon() {
        let ff = force_field();
        let contrib = InversionContrib::new(&ff, 0, 1, 2, 3, 6, false, 1.0);
        let pos = [
            0.1, -0.2, 0.3, //
            1.0, 0.0, 0.1, //
            -0.2, 1.1, 0.2, //
            0.2, -0.1, 1.3,
        ];
        let mut grad = [0.0; 12];

        contrib.get_grad(&pos, &mut grad);

        let step = 1.0e-6;
        for idx in 0..pos.len() {
            let mut plus = pos;
            let mut minus = pos;
            plus[idx] += step;
            minus[idx] -= step;
            let expected = (contrib.get_energy(&plus) - contrib.get_energy(&minus)) / (2.0 * step);
            assert_close_tol(grad[idx], expected, 1.0e-6);
        }
    }

    #[test]
    fn uff_inversioncontrib_getgrad_uses_group5_c2_source_formula() {
        let ff = force_field();
        let contrib = InversionContrib::new(&ff, 0, 1, 2, 3, 15, false, 1.0);
        let pos = [
            0.3, -0.1, 0.2, //
            1.2, 0.1, -0.1, //
            -0.4, 1.0, 0.3, //
            0.1, 0.2, 1.4,
        ];
        let mut grad = [0.0; 12];

        contrib.get_grad(&pos, &mut grad);

        let expected = [
            -2.1944224745433836,
            -0.75237341984344985,
            -7.0848497035257756,
            -1.5219937548842501,
            0.46576402965078911,
            5.5415831496996066,
            -0.46185626241615285,
            -0.15835071854268082,
            -1.4911359329435836,
            4.1782724918437868,
            0.44496010873534159,
            3.034402486769753,
        ];
        for idx in 0..grad.len() {
            assert_close_tol(grad[idx], expected[idx], 1.0e-12);
        }
    }

    #[test]
    fn uff_inversioncontrib_getgrad_returns_without_touching_grad_for_zero_distance() {
        let ff = force_field();
        let contrib = InversionContrib::new(&ff, 0, 1, 2, 3, 6, false, 1.0);
        let pos = [
            0.0, 0.0, 0.0, //
            0.0, 0.0, 0.0, //
            1.0, 0.0, 0.0, //
            0.0, 1.0, 0.0,
        ];
        let mut grad = [2.0; 12];

        contrib.get_grad(&pos, &mut grad);

        assert_eq!(grad, [2.0; 12]);
    }

    #[test]
    fn uff_inversioncontribs_constructor_stores_owner_and_starts_empty() {
        let ff = force_field();

        let contribs = InversionContribs::new(&ff);

        assert_eq!(contribs.owner(), &ff as *const ForceField);
        assert!(contribs.contribs().is_empty());
    }

    #[test]
    fn uff_inversioncontribs_empty_reflects_contrib_vector_state() {
        let ff = force_field();
        let mut contribs = InversionContribs::new(&ff);

        assert!(contribs.empty());

        contribs.add_contrib(0, 1, 2, 3, 6, false, 1.0);

        assert!(!contribs.empty());
    }

    #[test]
    fn uff_inversioncontribs_size_reflects_contrib_vector_len() {
        let ff = force_field();
        let mut contribs = InversionContribs::new(&ff);

        assert_eq!(contribs.size(), 0);

        contribs.add_contrib(0, 1, 2, 3, 6, false, 1.0);
        contribs.add_contrib(3, 1, 2, 0, 15, false, 0.5);

        assert_eq!(contribs.size(), 2);
    }

    #[test]
    fn uff_inversioncontribs_addcontrib_stores_indices_coefficients_and_scaled_force() {
        let ff = force_field();
        let mut contribs = InversionContribs::new(&ff);
        let scale = 0.25;
        let (force_constant, c0, c1, c2) = calc_inversion_coefficients(6, true);

        contribs.add_contrib(0, 1, 2, 3, 6, true, scale);

        assert_eq!(contribs.contribs().len(), 1);
        let contrib = &contribs.contribs()[0];
        assert_eq!(contrib.idx1(), 0);
        assert_eq!(contrib.idx2(), 1);
        assert_eq!(contrib.idx3(), 2);
        assert_eq!(contrib.idx4(), 3);
        assert_eq!(contrib.at2_atomic_num(), 6);
        assert!(contrib.is_c_bound_to_o());
        assert_close(contrib.c0(), c0);
        assert_close(contrib.c1(), c1);
        assert_close(contrib.c2(), c2);
        assert_close(contrib.force_constant(), force_constant * scale);
    }

    #[test]
    fn uff_inversioncontribs_addcontrib_uses_group5_coefficients() {
        let ff = force_field();
        let mut contribs = InversionContribs::new(&ff);
        let (force_constant, c0, c1, c2) = calc_inversion_coefficients(15, false);

        contribs.add_contrib(0, 1, 2, 3, 15, false, 2.0);

        let contrib = &contribs.contribs()[0];
        assert_eq!(contrib.at2_atomic_num(), 15);
        assert!(!contrib.is_c_bound_to_o());
        assert_close(contrib.c0(), c0);
        assert_close(contrib.c1(), c1);
        assert_close(contrib.c2(), c2);
        assert_close(contrib.force_constant(), force_constant * 2.0);
    }

    #[test]
    #[should_panic]
    fn uff_inversioncontribs_addcontrib_rejects_index_out_of_range() {
        let ff = force_field();
        let mut contribs = InversionContribs::new(&ff);

        contribs.add_contrib(0, 1, 2, 4, 6, false, 1.0);
    }

    #[test]
    fn uff_inversioncontribs_getenergy_empty_contribs_returns_zero_accumulator() {
        let ff = force_field();
        let contribs = InversionContribs::new(&ff);

        assert_close(contribs.get_energy(&[0.0, 0.0, 0.0]), 0.0);
    }

    #[test]
    fn uff_inversioncontribs_getenergy_matches_source_formula_for_single_contrib() {
        let ff = force_field();
        let mut contribs = InversionContribs::new(&ff);
        contribs.add_contrib(0, 1, 2, 3, 6, false, 1.0);
        let pos = [
            0.0, 0.0, 0.0, //
            1.0, 0.0, 0.0, //
            0.0, 1.0, 0.0, //
            0.0, 0.0, 1.0,
        ];
        let contrib = &contribs.contribs()[0];
        let sin_y = 0.5_f64.sqrt();
        let cos2_w = 2.0 * sin_y * sin_y - 1.0;
        let expected = contrib.force_constant()
            * (contrib.c0() + contrib.c1() * sin_y + contrib.c2() * cos2_w);

        assert_close(contribs.get_energy(&pos), expected);
    }

    #[test]
    fn uff_inversioncontribs_getenergy_accumulates_multiple_contribs_in_order() {
        let ff = force_field();
        let mut contribs = InversionContribs::new(&ff);
        contribs.add_contrib(0, 1, 2, 3, 6, false, 1.0);
        contribs.add_contrib(3, 1, 2, 0, 15, false, 0.5);
        let pos = [
            0.2, -0.1, 0.4, //
            1.1, 0.0, -0.2, //
            -0.3, 0.9, 0.1, //
            0.4, 0.2, 1.2,
        ];
        let expected = contribs
            .contribs()
            .iter()
            .map(|contrib| {
                let p1 = ForceFieldVec3::new(
                    pos[3 * contrib.idx1()],
                    pos[3 * contrib.idx1() + 1],
                    pos[3 * contrib.idx1() + 2],
                );
                let p2 = ForceFieldVec3::new(
                    pos[3 * contrib.idx2()],
                    pos[3 * contrib.idx2() + 1],
                    pos[3 * contrib.idx2() + 2],
                );
                let p3 = ForceFieldVec3::new(
                    pos[3 * contrib.idx3()],
                    pos[3 * contrib.idx3() + 1],
                    pos[3 * contrib.idx3() + 2],
                );
                let p4 = ForceFieldVec3::new(
                    pos[3 * contrib.idx4()],
                    pos[3 * contrib.idx4() + 1],
                    pos[3 * contrib.idx4() + 2],
                );
                let cos_y = calculate_cos_y(p1, p2, p3, p4);
                let sin_y_sq = 1.0 - cos_y * cos_y;
                let sin_y = if sin_y_sq > 0.0 { sin_y_sq.sqrt() } else { 0.0 };
                let cos2_w = 2.0 * sin_y * sin_y - 1.0;
                contrib.force_constant()
                    * (contrib.c0() + contrib.c1() * sin_y + contrib.c2() * cos2_w)
            })
            .sum::<f64>();

        assert_close(contribs.get_energy(&pos), expected);
    }

    #[test]
    fn uff_inversioncontribs_getenergy_uses_zero_sin_branch_when_cos_y_is_one() {
        let ff = force_field();
        let mut contribs = InversionContribs::new(&ff);
        contribs.add_contrib(0, 1, 2, 3, 15, false, 1.0);
        let pos = [
            1.0, 0.0, 0.0, //
            0.0, 0.0, 0.0, //
            0.0, 1.0, 0.0, //
            0.0, 0.0, 1.0,
        ];
        let contrib = &contribs.contribs()[0];
        let expected = contrib.force_constant() * (contrib.c0() - contrib.c2());

        assert_close(contribs.get_energy(&pos), expected);
    }

    #[test]
    #[should_panic(expected = "bad vector")]
    fn uff_inversioncontribs_getenergy_rejects_empty_position_vector() {
        let ff = force_field();
        let contribs = InversionContribs::new(&ff);

        let _ = contribs.get_energy(&[]);
    }

    #[test]
    #[should_panic(expected = "bad vector")]
    fn uff_inversioncontribs_getenergy_rejects_short_position_vector_for_contrib() {
        let ff = force_field();
        let mut contribs = InversionContribs::new(&ff);
        contribs.add_contrib(0, 1, 2, 3, 6, false, 1.0);

        let _ = contribs.get_energy(&[0.0, 0.0, 0.0]);
    }

    #[test]
    #[should_panic(expected = "no owner")]
    fn uff_inversioncontribs_getenergy_rejects_missing_owner() {
        let contribs = InversionContribs::default();

        let _ = contribs.get_energy(&[0.0, 0.0, 0.0]);
    }

    #[test]
    fn uff_inversioncontribs_getgrad_empty_contribs_leaves_grad_unchanged() {
        let ff = force_field();
        let contribs = InversionContribs::new(&ff);
        let pos = [0.0, 0.0, 0.0];
        let mut grad = [3.0, 4.0, 5.0];

        contribs.get_grad(&pos, &mut grad);

        assert_eq!(grad, [3.0, 4.0, 5.0]);
    }

    #[test]
    fn uff_inversioncontribs_getgrad_matches_single_contrib_source_path() {
        let ff = force_field();
        let mut contribs = InversionContribs::new(&ff);
        contribs.add_contrib(0, 1, 2, 3, 6, false, 1.0);
        let single = InversionContrib::new(&ff, 0, 1, 2, 3, 6, false, 1.0);
        let pos = [
            0.1, -0.2, 0.3, //
            1.0, 0.0, 0.1, //
            -0.2, 1.1, 0.2, //
            0.2, -0.1, 1.3,
        ];
        let mut actual = [0.0; 12];
        let mut expected = [0.0; 12];

        contribs.get_grad(&pos, &mut actual);
        single.get_grad(&pos, &mut expected);

        for idx in 0..actual.len() {
            assert_close_tol(actual[idx], expected[idx], 1.0e-12);
        }
    }

    #[test]
    fn uff_inversioncontribs_getgrad_accumulates_multiple_contribs() {
        let ff = force_field();
        let mut contribs = InversionContribs::new(&ff);
        contribs.add_contrib(0, 1, 2, 3, 6, false, 1.0);
        contribs.add_contrib(3, 1, 2, 0, 15, false, 0.5);
        let single1 = InversionContrib::new(&ff, 0, 1, 2, 3, 6, false, 1.0);
        let single2 = InversionContrib::new(&ff, 3, 1, 2, 0, 15, false, 0.5);
        let pos = [
            0.2, -0.1, 0.4, //
            1.1, 0.0, -0.2, //
            -0.3, 0.9, 0.1, //
            0.4, 0.2, 1.2,
        ];
        let mut actual = [0.0; 12];
        let mut expected = [0.0; 12];

        contribs.get_grad(&pos, &mut actual);
        single1.get_grad(&pos, &mut expected);
        single2.get_grad(&pos, &mut expected);

        for idx in 0..actual.len() {
            assert_close_tol(actual[idx], expected[idx], 1.0e-12);
        }
    }

    #[test]
    fn uff_inversioncontribs_getgrad_returns_from_whole_function_on_zero_distance() {
        let ff = force_field();
        let mut contribs = InversionContribs::new(&ff);
        contribs.add_contrib(0, 1, 2, 3, 6, false, 1.0);
        contribs.add_contrib(3, 1, 2, 0, 15, false, 1.0);
        let pos = [
            0.0, 0.0, 0.0, //
            0.0, 0.0, 0.0, //
            1.0, 0.0, 0.0, //
            0.0, 1.0, 0.0,
        ];
        let mut grad = [2.0; 12];

        contribs.get_grad(&pos, &mut grad);

        assert_eq!(grad, [2.0; 12]);
    }

    #[test]
    #[should_panic(expected = "bad vector")]
    fn uff_inversioncontribs_getgrad_rejects_empty_position_vector() {
        let ff = force_field();
        let contribs = InversionContribs::new(&ff);
        let mut grad = [0.0, 0.0, 0.0];

        contribs.get_grad(&[], &mut grad);
    }

    #[test]
    #[should_panic(expected = "bad vector")]
    fn uff_inversioncontribs_getgrad_rejects_empty_gradient_vector() {
        let ff = force_field();
        let contribs = InversionContribs::new(&ff);
        let mut grad = [];

        contribs.get_grad(&[0.0, 0.0, 0.0], &mut grad);
    }

    #[test]
    #[should_panic(expected = "bad vector")]
    fn uff_inversioncontribs_getgrad_rejects_short_position_vector_for_contrib() {
        let ff = force_field();
        let mut contribs = InversionContribs::new(&ff);
        contribs.add_contrib(0, 1, 2, 3, 6, false, 1.0);
        let mut grad = [0.0; 12];

        contribs.get_grad(&[0.0, 0.0, 0.0], &mut grad);
    }

    #[test]
    #[should_panic(expected = "bad vector")]
    fn uff_inversioncontribs_getgrad_rejects_short_gradient_vector_for_contrib() {
        let ff = force_field();
        let mut contribs = InversionContribs::new(&ff);
        contribs.add_contrib(0, 1, 2, 3, 6, false, 1.0);
        let pos = [
            0.1, -0.2, 0.3, //
            1.0, 0.0, 0.1, //
            -0.2, 1.1, 0.2, //
            0.2, -0.1, 1.3,
        ];
        let mut grad = [0.0; 3];

        contribs.get_grad(&pos, &mut grad);
    }

    #[test]
    #[should_panic(expected = "no owner")]
    fn uff_inversioncontribs_getgrad_rejects_missing_owner() {
        let contribs = InversionContribs::default();
        let pos = [0.0, 0.0, 0.0];
        let mut grad = [0.0, 0.0, 0.0];

        contribs.get_grad(&pos, &mut grad);
    }
}
