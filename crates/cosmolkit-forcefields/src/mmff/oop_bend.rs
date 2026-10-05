//! Complete private MMFF OopBend contribution.
//! Pinned RDKit351f8f378f8ad6bbd517980c38896e66bf907af8 OopBend.cpp/h, BSD.
//! ROOT authorized the private four-path implementation and typed address
//! safety policy; private source conditions passed immediate tests.
//! Independent p1 acceptance and ROOT review remain pending.

use super::{
    angle_bend::{point_at, source_sin_theta},
    bond_stretch::source_term_indices,
    numerical::{calc_oop_bend_energy, calc_oop_bend_force_constant, calc_oop_chi},
    params::MmffOop,
};
use crate::{
    geometry::Point3,
    kernel::{
        EvaluationContext, ForceField, ForceFieldContribution, ForceFieldKernelError,
        TorsionIndexArgument,
    },
    uff::params::{DEG2RAD, RAD2DEG, clip_to_one, is_double_zero},
};

// RDKit✔️✔️:   std::vector<int> d_at1Idxs, d_at2Idxs, d_at3Idxs, d_at4Idxs;
// RDKit✔️✔️:   std::vector<double> d_koop;
// Native source width is i32; u32->i32 is bijective, not a 16-bit narrow.
#[derive(Clone, Debug)]
pub(super) struct OopBendContrib {
    at1_idxs: Vec<i32>,
    at2_idxs: Vec<i32>,
    at3_idxs: Vec<i32>,
    at4_idxs: Vec<i32>,
    koop: Vec<f64>,
}

fn checked_address(index: i32, buffer_len: usize) -> Result<usize, ForceFieldKernelError> {
    // Source implicit signed-address anchors (OopBend.cpp getSingleGrad):
    // RDKit❗✔️:   RDGeom::Point3D iPoint(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
    // RDKit❗✔️:                          pos[3 * d_at1Idx + 2]);
    // RDKit❗✔️:   double *g1 = &(grad[3 * d_at1Idx]);
    // Authorized Native safety, without a C++-defined error or numeric parity oracle:
    // reject negative index, signed 3*index/+2 overflow or absent complete row.
    // Valid defined source inputs keep their exact addresses. ROOT authorized this policy.
    // Complexity: fixed checked integer operations, O(1), no allocation.
    let end = index.checked_mul(3).and_then(|base| base.checked_add(2));
    match end {
        Some(end) if index >= 0 && (end as usize) < buffer_len => Ok(end as usize - 2),
        _ => Err(ForceFieldKernelError::OopAddressOutsideDefinedSource { index, buffer_len }),
    }
}

fn source_theta_factors(cos_theta: f64) -> (f64, f64) {
    // RDKit✔️✔️:   double sinThetaSq = std::max(1.0 - cosTheta * cosTheta, 1.0e-8);
    // RDKit✔️✔️:   double sinTheta =
    // RDKit✔️✔️:       std::max(((sinThetaSq > 0.0) ? sqrt(sinThetaSq) : 0.0), 1.0e-8);
    // Behavior: first stdmax preserves NaN; the conditional sqrt maps it to0.
    // This square-before-root floor gives1e-4 at cos=+/-1; Angle differs here.
    // Complexity: fixed scalar arithmetic/two branches, O(1), no allocation.
    let candidate = 1.0 - cos_theta * cos_theta;
    let squared = if candidate < 1.0e-8 {
        1.0e-8
    } else {
        candidate
    };
    let root = if squared > 0.0 { squared.sqrt() } else { 0.0 };
    let sine = if root < 1.0e-8 { 1.0e-8 } else { root };
    (squared, sine)
}

impl OopBendContrib {
    pub(super) fn new(_owner: &ForceField<'_>) -> Self {
        // RDKit✔️✔️: OopBendContrib::OopBendContrib(ForceField *owner) {
        // RDKit✔️✔️:   PRECONDITION(owner, "bad owner");
        // RDKit✔️✔️:   dp_forceField = owner;
        // RDKit✔️✔️: }
        // Behavior: non-null owner borrow, existing context/rebind lifetime.
        // An ownerless raw-source default is outside the evaluable boundary.
        // Complexity: five empty Vecs, O(1), no allocation or owner clone.
        Self {
            at1_idxs: Vec::new(),
            at2_idxs: Vec::new(),
            at3_idxs: Vec::new(),
            at4_idxs: Vec::new(),
            koop: Vec::new(),
        }
    }

    pub(super) fn add_term(
        &mut self,
        positions: &[&mut [f64]],
        idx1: u32,
        idx2: u32,
        idx3: u32,
        idx4: u32,
        params: Option<&MmffOop>,
    ) -> Result<(), ForceFieldKernelError> {
        // RDKit✔️✔️: void OopBendContrib::addTerm(unsigned int idx1,
        // RDKit✔️✔️:                              unsigned int idx2,
        // RDKit✔️✔️:                              unsigned int idx3,
        // RDKit✔️✔️:                              unsigned int idx4,
        // RDKit✔️✔️:                              const ForceFields::MMFF::MMFFOop *mmffOopParams) {
        // RDKit✔️✔️:   PRECONDITION(mmffOopParams, "no OOP parameters");
        // RDKit✔️✔️:   PRECONDITION((idx1 != idx2) && (idx1 != idx3) && (idx1 != idx4) &&
        // RDKit✔️✔️:                    (idx2 != idx3) && (idx2 != idx4) && (idx3 != idx4),
        // RDKit✔️✔️:                "degenerate points");
        // RDKit✔️✔️:   URANGE_CHECK(idx1, dp_forceField->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx2, dp_forceField->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx3, dp_forceField->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx4, dp_forceField->positions().size());
        // RDKit✔️✔️:
        // RDKit✔️✔️:   d_at1Idxs.push_back(idx1);
        // RDKit✔️✔️:   d_at2Idxs.push_back(idx2);
        // RDKit✔️✔️:   d_at3Idxs.push_back(idx3);
        // RDKit✔️✔️:   d_at4Idxs.push_back(idx4);
        // RDKit✔️✔️:   d_koop.push_back(mmffOopParams->koop);
        // RDKit✔️✔️: }
        // Behavior: missing parameter FIRST, then all6distinct before4rawbounds
        // before any append. Copy the source koop via its existing unique getter.
        // Reuse existing six-pair/four-role error vocabulary with exact same
        // source category, expression and message; new missing-parameter
        // precondition and Native address variants have ROOT authorization.
        // Complexity: fixed checks and five amortized O(1) array appends.
        let params = params.ok_or(ForceFieldKernelError::OopParametersMissing)?;
        if idx1 == idx2
            || idx1 == idx3
            || idx1 == idx4
            || idx2 == idx3
            || idx2 == idx4
            || idx3 == idx4
        {
            return Err(ForceFieldKernelError::TorsionDegeneratePoints);
        }
        for (argument, index) in [
            (TorsionIndexArgument::First, idx1),
            (TorsionIndexArgument::Second, idx2),
            (TorsionIndexArgument::Third, idx3),
            (TorsionIndexArgument::Fourth, idx4),
        ] {
            if index as usize >= positions.len() {
                return Err(ForceFieldKernelError::TorsionIndexOutOfRange {
                    argument,
                    index,
                    upper_bound: positions.len(),
                });
            }
        }
        self.at1_idxs.push(idx1 as i32);
        self.at2_idxs.push(idx2 as i32);
        self.at3_idxs.push(idx3 as i32);
        self.at4_idxs.push(idx4 as i32);
        self.koop.push(calc_oop_bend_force_constant(params));
        Ok(())
    }

    fn get_single_grad(
        &self,
        coordinates: &[f64],
        gradient: &mut [f64],
        term_idx: usize,
    ) -> Result<(), ForceFieldKernelError> {
        // RDKit✔️✔️: void OopBendContrib::getSingleGrad(double *pos, double *grad, unsigned int termIdx) const {
        // RDKit✔️✔️:
        // RDKit✔️✔️:   const int d_at1Idx = d_at1Idxs[termIdx];
        // RDKit✔️✔️:   const int d_at2Idx = d_at2Idxs[termIdx];
        // RDKit✔️✔️:   const int d_at3Idx = d_at3Idxs[termIdx];
        // RDKit✔️✔️:   const int d_at4Idx = d_at4Idxs[termIdx];
        // RDKit✔️✔️:
        // RDKit✔️✔️:   RDGeom::Point3D iPoint(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
        // RDKit✔️✔️:                          pos[3 * d_at1Idx + 2]);
        // RDKit✔️✔️:   RDGeom::Point3D jPoint(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
        // RDKit✔️✔️:                          pos[3 * d_at2Idx + 2]);
        // RDKit✔️✔️:   RDGeom::Point3D kPoint(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
        // RDKit✔️✔️:                          pos[3 * d_at3Idx + 2]);
        // RDKit✔️✔️:   RDGeom::Point3D lPoint(pos[3 * d_at4Idx], pos[3 * d_at4Idx + 1],
        // RDKit✔️✔️:                          pos[3 * d_at4Idx + 2]);
        // RDKit✔️✔️:   double *g1 = &(grad[3 * d_at1Idx]);
        // RDKit✔️✔️:   double *g2 = &(grad[3 * d_at2Idx]);
        // RDKit✔️✔️:   double *g3 = &(grad[3 * d_at3Idx]);
        // RDKit✔️✔️:   double *g4 = &(grad[3 * d_at4Idx]);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   RDGeom::Point3D rJI = iPoint - jPoint;
        // RDKit✔️✔️:   RDGeom::Point3D rJK = kPoint - jPoint;
        // RDKit✔️✔️:   RDGeom::Point3D rJL = lPoint - jPoint;
        // RDKit✔️✔️:   double dJI = rJI.length();
        // RDKit✔️✔️:   double dJK = rJK.length();
        // RDKit✔️✔️:   double dJL = rJL.length();
        // RDKit✔️✔️:   if (isDoubleZero(dJI) || isDoubleZero(dJK) || isDoubleZero(dJL)) {
        // RDKit✔️✔️:     return;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   rJI /= dJI;
        // RDKit✔️✔️:   rJK /= dJK;
        // RDKit✔️✔️:   rJL /= dJL;
        // RDKit✔️✔️:
        // RDKit✔️✔️:   RDGeom::Point3D n = (-rJI).crossProduct(rJK);
        // RDKit✔️✔️:   n /= n.length();
        // RDKit✔️✔️:   double const c2 = MDYNE_A_TO_KCAL_MOL * DEG2RAD * DEG2RAD;
        // RDKit✔️✔️:   double sinChi = rJL.dotProduct(n);
        // RDKit✔️✔️:   clipToOne(sinChi);
        // RDKit✔️✔️:   double cosChiSq = 1.0 - sinChi * sinChi;
        // RDKit✔️✔️:   double cosChi = std::max(((cosChiSq > 0.0) ? sqrt(cosChiSq) : 0.0), 1.0e-8);
        // RDKit✔️✔️:   double chi = RAD2DEG * asin(sinChi);
        // RDKit✔️✔️:   double cosTheta = rJI.dotProduct(rJK);
        // RDKit✔️✔️:   clipToOne(cosTheta);
        // RDKit✔️✔️:   double sinThetaSq = std::max(1.0 - cosTheta * cosTheta, 1.0e-8);
        // RDKit✔️✔️:   double sinTheta =
        // RDKit✔️✔️:       std::max(((sinThetaSq > 0.0) ? sqrt(sinThetaSq) : 0.0), 1.0e-8);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   double dE_dChi = RAD2DEG * c2 * d_koop[termIdx] * chi;
        // RDKit✔️✔️:   RDGeom::Point3D t1 = rJL.crossProduct(rJK);
        // RDKit✔️✔️:   RDGeom::Point3D t2 = rJI.crossProduct(rJL);
        // RDKit✔️✔️:   RDGeom::Point3D t3 = rJK.crossProduct(rJI);
        // RDKit✔️✔️:   double term1 = cosChi * sinTheta;
        // RDKit✔️✔️:   double term2 = sinChi / (cosChi * sinThetaSq);
        // RDKit✔️✔️:   double tg1[3] = {(t1.x / term1 - (rJI.x - rJK.x * cosTheta) * term2) / dJI,
        // RDKit✔️✔️:                    (t1.y / term1 - (rJI.y - rJK.y * cosTheta) * term2) / dJI,
        // RDKit✔️✔️:                    (t1.z / term1 - (rJI.z - rJK.z * cosTheta) * term2) / dJI};
        // RDKit✔️✔️:   double tg3[3] = {(t2.x / term1 - (rJK.x - rJI.x * cosTheta) * term2) / dJK,
        // RDKit✔️✔️:                    (t2.y / term1 - (rJK.y - rJI.y * cosTheta) * term2) / dJK,
        // RDKit✔️✔️:                    (t2.z / term1 - (rJK.z - rJI.z * cosTheta) * term2) / dJK};
        // RDKit✔️✔️:   double tg4[3] = {(t3.x / term1 - rJL.x * sinChi / cosChi) / dJL,
        // RDKit✔️✔️:                    (t3.y / term1 - rJL.y * sinChi / cosChi) / dJL,
        // RDKit✔️✔️:                    (t3.z / term1 - rJL.z * sinChi / cosChi) / dJL};
        // RDKit✔️✔️:   for (unsigned int i = 0; i < 3; ++i) {
        // RDKit✔️✔️:     g1[i] += dE_dChi * tg1[i];
        // RDKit✔️✔️:     g2[i] += -dE_dChi * (tg1[i] + tg3[i] + tg4[i]);
        // RDKit✔️✔️:     g3[i] += dE_dChi * tg3[i];
        // RDKit✔️✔️:     g4[i] += dE_dChi * tg4[i];
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // Behavior: complete independent OOP source body; valid term_idx is
        // guaranteed by the sole signed32 caller loop. Read all four physical
        // points BEFORE four gradient addresses, then calculate all3lengths.
        // Source strict zero-arm checks return without any gradient write.
        // Negation uses existing scaled(-1.0), exactly source component *= -1.
        // cosChi's conditional square/root/floor is the existing Angle helper
        // with sinChi as input; Theta's floor BEFORE sqrt is distinct above.
        // Preserve source three temporary vectors, three tg arrays and twelve
        // ordered writes: axis outer, g1/g2/g3/g4 inner, center left sum.
        // Complexity: O(1) per term; stack values/arrays, no allocation/cache,
        // buffering, scans, gradient copy or whole owner clone. Native address
        // checks are explicit constant work under recorded ROOT authorization.
        let at1 = self.at1_idxs[term_idx];
        let at2 = self.at2_idxs[term_idx];
        let at3 = self.at3_idxs[term_idx];
        let at4 = self.at4_idxs[term_idx];
        let p1 = point_at(coordinates, checked_address(at1, coordinates.len())? / 3);
        let p2 = point_at(coordinates, checked_address(at2, coordinates.len())? / 3);
        let p3 = point_at(coordinates, checked_address(at3, coordinates.len())? / 3);
        let p4 = point_at(coordinates, checked_address(at4, coordinates.len())? / 3);
        let g1 = checked_address(at1, gradient.len())?;
        let g2 = checked_address(at2, gradient.len())?;
        let g3 = checked_address(at3, gradient.len())?;
        let g4 = checked_address(at4, gradient.len())?;
        let mut r_ji = Point3::difference(&p1, &p2);
        let mut r_jk = Point3::difference(&p3, &p2);
        let mut r_jl = Point3::difference(&p4, &p2);
        let d_ji = r_ji.length();
        let d_jk = r_jk.length();
        let d_jl = r_jl.length();
        if is_double_zero(d_ji) || is_double_zero(d_jk) || is_double_zero(d_jl) {
            return Ok(());
        }
        r_ji.divide_assign(d_ji);
        r_jk.divide_assign(d_jk);
        r_jl.divide_assign(d_jl);
        let mut n = r_ji.scaled(-1.0).cross_product(&r_jk);
        let n_length = n.length();
        n.divide_assign(n_length);
        // RDKit✔️✔️: constexpr double MDYNE_A_TO_KCAL_MOL = 143.9325;
        let c2 = 143.9325 * DEG2RAD * DEG2RAD;
        let mut sin_chi = r_jl.dot_product(&n);
        clip_to_one(&mut sin_chi);
        let cos_chi = source_sin_theta(sin_chi);
        let chi = RAD2DEG * sin_chi.asin();
        let mut cos_theta = r_ji.dot_product(&r_jk);
        clip_to_one(&mut cos_theta);
        let (sin_theta_sq, sin_theta) = source_theta_factors(cos_theta);
        let de_dchi = RAD2DEG * c2 * self.koop[term_idx] * chi;
        let t1 = r_jl.cross_product(&r_jk);
        let t2 = r_ji.cross_product(&r_jl);
        let t3 = r_jk.cross_product(&r_ji);
        let term1 = cos_chi * sin_theta;
        let term2 = sin_chi / (cos_chi * sin_theta_sq);
        let tg1 = [
            (t1.x / term1 - (r_ji.x - r_jk.x * cos_theta) * term2) / d_ji,
            (t1.y / term1 - (r_ji.y - r_jk.y * cos_theta) * term2) / d_ji,
            (t1.z / term1 - (r_ji.z - r_jk.z * cos_theta) * term2) / d_ji,
        ];
        let tg3 = [
            (t2.x / term1 - (r_jk.x - r_ji.x * cos_theta) * term2) / d_jk,
            (t2.y / term1 - (r_jk.y - r_ji.y * cos_theta) * term2) / d_jk,
            (t2.z / term1 - (r_jk.z - r_ji.z * cos_theta) * term2) / d_jk,
        ];
        let tg4 = [
            (t3.x / term1 - r_jl.x * sin_chi / cos_chi) / d_jl,
            (t3.y / term1 - r_jl.y * sin_chi / cos_chi) / d_jl,
            (t3.z / term1 - r_jl.z * sin_chi / cos_chi) / d_jl,
        ];
        for axis in 0..3 {
            gradient[g1 + axis] += de_dchi * tg1[axis];
            gradient[g2 + axis] += -de_dchi * (tg1[axis] + tg3[axis] + tg4[axis]);
            gradient[g3 + axis] += de_dchi * tg3[axis];
            gradient[g4 + axis] += de_dchi * tg4[axis];
        }
        Ok(())
    }
}

impl ForceFieldContribution for OopBendContrib {
    fn get_energy(
        &self,
        context: &mut EvaluationContext<'_>,
    ) -> Result<f64, ForceFieldKernelError> {
        // RDKit✔️✔️: double OopBendContrib::getEnergy(double *pos) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:
        // RDKit✔️✔️:   const int numTerms = d_at1Idxs.size();
        // RDKit✔️✔️:   double totalEnergy = 0.0;
        // RDKit✔️✔️:   for (int i = 0; i < numTerms; ++i) {
        // RDKit✔️✔️:     const int d_at1Idx = d_at1Idxs[i];
        // RDKit✔️✔️:     const int d_at2Idx = d_at2Idxs[i];
        // RDKit✔️✔️:     const int d_at3Idx = d_at3Idxs[i];
        // RDKit✔️✔️:     const int d_at4Idx = d_at4Idxs[i];
        // RDKit✔️✔️:
        // RDKit✔️✔️:     RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
        // RDKit✔️✔️:                        pos[3 * d_at1Idx + 2]);
        // RDKit✔️✔️:     RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
        // RDKit✔️✔️:                        pos[3 * d_at2Idx + 2]);
        // RDKit✔️✔️:     RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
        // RDKit✔️✔️:                        pos[3 * d_at3Idx + 2]);
        // RDKit✔️✔️:     RDGeom::Point3D p4(pos[3 * d_at4Idx], pos[3 * d_at4Idx + 1],
        // RDKit✔️✔️:                        pos[3 * d_at4Idx + 2]);
        // RDKit✔️✔️:
        // RDKit✔️✔️:     totalEnergy += Utils::calcOopBendEnergy(Utils::calcOopChi(p1, p2, p3, p4), d_koop[i]);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return totalEnergy;
        // RDKit✔️✔️: }
        // Behavior: source signed32 loop, four i32 endpoints/direct coordinates;
        // source positive normal/asin energy helper unchanged. NO cache call or
        // zero-length fallback. Undefined addresses follow the authorized Native
        // typed boundary above, without a C++ numerical/error oracle.
        // Complexity: O(max(0,signed32terms)), O(1) stack, no allocation/cache.
        let mut total_energy = 0.0;
        for term in source_term_indices(self.at1_idxs.len()) {
            let term = term as usize;
            let coords = context.coordinates();
            let p1 = point_at(
                coords,
                checked_address(self.at1_idxs[term], coords.len())? / 3,
            );
            let p2 = point_at(
                coords,
                checked_address(self.at2_idxs[term], coords.len())? / 3,
            );
            let p3 = point_at(
                coords,
                checked_address(self.at3_idxs[term], coords.len())? / 3,
            );
            let p4 = point_at(
                coords,
                checked_address(self.at4_idxs[term], coords.len())? / 3,
            );
            total_energy += calc_oop_bend_energy(calc_oop_chi(&p1, &p2, &p3, &p4), self.koop[term]);
        }
        Ok(total_energy)
    }
    fn get_grad(
        &self,
        context: &mut EvaluationContext<'_>,
        gradient: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        // RDKit✔️✔️: void OopBendContrib::getGrad(double* pos, double* grad) const {
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:   PRECONDITION(grad, "bad vector");
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:
        // RDKit✔️✔️:   const int numTerms = d_at1Idxs.size();
        // RDKit✔️✔️:   for (int i =0; i < numTerms; i++) {
        // RDKit✔️✔️:     getSingleGrad(pos, grad, i);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // Behavior: pos/grad references and owner-derived context enforce the
        // source ordered non-null valid caller boundary. Shared signed32 loop
        // invokes complete single-term body; no transaction masks partial writes.
        // Complexity: O(max(0,signed32terms))/O(1), no copies or allocation.
        for term in source_term_indices(self.at1_idxs.len()) {
            self.get_single_grad(context.coordinates(), gradient, term as usize)?;
        }
        Ok(())
    }
    fn copy(&self) -> Box<dyn ForceFieldContribution> {
        // RDKit✔️✔️:   OopBendContrib *copy() const override { return new OopBendContrib(*this); }
        // Behavior: five arrays independently cloned; existing ForceField copy
        // owns position reset/initialization/rebind. No owner clone or shared Vec.
        // Complexity: O(terms) elements/storage, five populated Vec allocations
        // plus the contribution box, matching source copy lifetime/shape.
        Box::new(self.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernel::{cf3d_bld_b05_calc_energy, cf3d_bld_b05_calc_grad};

    // Frozen independent pinned-source reductions, not values from Rust.
    const COORDS: [f64; 12] = [2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 3.0, 4.0, 5.0];
    const PLANAR: [f64; 12] = [2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 3.0, 4.0, 0.0];
    const GRAD: [f64; 12] = [
        0.0,
        0.0,
        -67.82659269192193,
        13.565318538384387,
        18.087091384512515,
        105.50803307632302,
        0.0,
        0.0,
        -60.29030461504172,
        -13.565318538384387,
        -18.087091384512515,
        22.60886423064065,
    ];
    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 1e-10,
            "{actual:?} != {expected:?}"
        );
    }
    fn close_grad(actual: &[f64], expected: &[f64]) {
        assert_eq!(actual.len(), expected.len());
        for (&a, &e) in actual.iter().zip(expected) {
            close(a, e);
        }
    }
    fn append(
        term: &mut OopBendContrib,
        field: &ForceField<'_>,
        idx: [u32; 4],
        koop: f64,
    ) -> Result<(), ForceFieldKernelError> {
        term.add_term(
            field.positions(),
            idx[0],
            idx[1],
            idx[2],
            idx[3],
            Some(&MmffOop { koop }),
        )
    }
    fn make(constants: &[f64]) -> OopBendContrib {
        let mut rows = [[0.0; 3]; 4];
        let mut field = ForceField::new(3);
        for row in &mut rows {
            field.positions_mut().push(row);
        }
        let mut term = OopBendContrib::new(&field);
        for &k in constants {
            append(&mut term, &field, [0, 1, 2, 3], k).unwrap();
        }
        term
    }
    fn energy(term: &dyn ForceFieldContribution, coords: &[f64]) -> f64 {
        let mut unused_cache = [];
        term.get_energy(&mut EvaluationContext::for_test(
            coords,
            &mut unused_cache,
            4,
        ))
        .unwrap()
    }
    fn gradient(term: &dyn ForceFieldContribution, coords: &[f64], initial: &[f64]) -> Vec<f64> {
        let mut unused_cache = [];
        let mut result = initial.to_vec();
        term.get_grad(
            &mut EvaluationContext::for_test(coords, &mut unused_cache, 4),
            &mut result,
        )
        .unwrap();
        result
    }
    fn lengths(term: &OopBendContrib) -> [usize; 5] {
        [
            term.at1_idxs.len(),
            term.at2_idxs.len(),
            term.at3_idxs.len(),
            term.at4_idxs.len(),
            term.koop.len(),
        ]
    }

    #[test]
    fn mmff_oop_bend_empty_owner_and_copy() {
        let term = OopBendContrib::new(&ForceField::new(3));
        assert_eq!(lengths(&term), [0; 5]);
        assert_eq!(
            [
                term.at1_idxs.capacity(),
                term.at2_idxs.capacity(),
                term.at3_idxs.capacity(),
                term.at4_idxs.capacity(),
                term.koop.capacity()
            ],
            [0; 5]
        );
        assert_eq!(energy(&term, &[]).to_bits(), 0);
        assert_eq!(
            gradient(&term, &[], &[1.0, -0.0, 3.0])
                .iter()
                .map(|v| v.to_bits())
                .collect::<Vec<_>>(),
            [1.0_f64.to_bits(), (-0.0_f64).to_bits(), 3.0_f64.to_bits()]
        );
        assert_eq!(energy(&term.clone(), &[]).to_bits(), 0);
        assert_eq!(energy(&*term.copy(), &[]).to_bits(), 0);
    }
    #[test]
    fn mmff_oop_bend_null_parameter_before_all_other_checks() {
        let field = ForceField::new(3);
        let mut term = OopBendContrib::new(&field);
        for idx in [[0, 0, 0, 0], [u32::MAX, u32::MAX, 2, 3], [0, 1, 2, 3]] {
            assert_eq!(
                term.add_term(field.positions(), idx[0], idx[1], idx[2], idx[3], None),
                Err(ForceFieldKernelError::OopParametersMissing)
            );
            assert_eq!(lengths(&term), [0; 5]);
        }
        assert_eq!(
            append(&mut term, &field, [0, 0, 0, 0], 2.0),
            Err(ForceFieldKernelError::TorsionDegeneratePoints)
        );
    }
    #[test]
    fn mmff_oop_bend_all_six_distinct_pairs_before_bounds() {
        let field = ForceField::new(3);
        let mut term = OopBendContrib::new(&field);
        for (a, b) in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)] {
            let mut idx = [10, 11, 12, 13];
            idx[b] = idx[a];
            assert_eq!(
                append(&mut term, &field, idx, 2.0),
                Err(ForceFieldKernelError::TorsionDegeneratePoints)
            );
            assert_eq!(lengths(&term), [0; 5]);
        }
    }
    #[test]
    fn mmff_oop_bend_four_raw_bounds_order_no_append() {
        let mut rows = [[0.0; 3]; 4];
        let mut extra = [0.0; 3];
        let mut field = ForceField::new(3);
        for row in &mut rows {
            field.positions_mut().push(row);
        }
        let mut term = OopBendContrib::new(&field);
        for (idx, argument, index) in [
            ([4, 5, 6, 7], TorsionIndexArgument::First, 4),
            ([0, 4, 5, 6], TorsionIndexArgument::Second, 4),
            ([0, 1, 4, 5], TorsionIndexArgument::Third, 4),
            ([0, 1, 2, u32::MAX], TorsionIndexArgument::Fourth, u32::MAX),
        ] {
            assert_eq!(
                append(&mut term, &field, idx, 2.0),
                Err(ForceFieldKernelError::TorsionIndexOutOfRange {
                    argument,
                    index,
                    upper_bound: 4
                })
            );
            assert_eq!(lengths(&term), [0; 5]);
        }
        append(&mut term, &field, [0, 1, 2, 3], 2.0).unwrap();
        assert_eq!(
            append(&mut term, &field, [0, 1, 2, 4], 9.0),
            Err(ForceFieldKernelError::TorsionIndexOutOfRange {
                argument: TorsionIndexArgument::Fourth,
                index: 4,
                upper_bound: 4
            })
        );
        assert_eq!(lengths(&term), [1; 5]);
        assert_eq!(term.koop, [2.0]);
        field.positions_mut().push(&mut extra);
        append(&mut term, &field, [0, 1, 2, 4], 9.0).unwrap();
        assert_eq!(term.at4_idxs, [3, 4]);
    }
    #[test]
    fn mmff_oop_bend_five_arrays_parameter_copy_and_growth() {
        let mut rows = [[0.0; 3]; 4];
        let mut field = ForceField::new(3);
        for row in &mut rows {
            field.positions_mut().push(row);
        }
        let mut term = OopBendContrib::new(&field);
        for i in 0..65 {
            let idx = if i % 2 == 0 {
                [0, 1, 2, 3]
            } else {
                [2, 1, 0, 3]
            };
            append(&mut term, &field, idx, i as f64 - 32.0).unwrap();
        }
        assert_eq!(lengths(&term), [65; 5]);
        for i in 0..65 {
            assert_eq!(
                [
                    term.at1_idxs[i],
                    term.at2_idxs[i],
                    term.at3_idxs[i],
                    term.at4_idxs[i]
                ],
                if i % 2 == 0 {
                    [0, 1, 2, 3]
                } else {
                    [2, 1, 0, 3]
                }
            );
            assert_eq!(term.koop[i], i as f64 - 32.0);
        }
    }
    #[test]
    fn mmff_oop_bend_fixed45_energy_and_all12_gradient() {
        let term = make(&[2.0]);
        close(energy(&term, &COORDS), 88.78480221623713);
        close_grad(&gradient(&term, &COORDS, &[0.0; 12]), &GRAD);
        // Numerical derivatives are supplementary; all12 fixed values above
        // come from the pinned-source closed reduction before implementation.
        for axis in 0..12 {
            let mut plus = COORDS;
            let mut minus = COORDS;
            plus[axis] += 1e-6;
            minus[axis] -= 1e-6;
            let derivative = (energy(&term, &plus) - energy(&term, &minus)) / 2e-6;
            assert!(
                (derivative - GRAD[axis]).abs() < 1e-6,
                "axis{axis}:{derivative}"
            );
        }
    }
    #[test]
    fn mmff_oop_bend_fixed30_energy_and_all12_gradient() {
        let coords = [
            2.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            3.0,
            0.0,
            0.0,
            1.7320508075688772,
            1.0,
        ];
        let term = make(&[2.0]);
        close(energy(&term, &coords), 39.459912096105384);
        close_grad(
            &gradient(&term, &coords, &[0.0; 12]),
            &[
                0.0,
                0.0,
                0.0,
                0.0,
                37.68144038440108,
                -15.024248735625621,
                0.0,
                0.0,
                -50.24192051253477,
                0.0,
                -37.68144038440108,
                65.2661692481604,
            ],
        );
    }
    #[test]
    fn mmff_oop_bend_planar_equilibrium_and_perpendicular_endpoint() {
        let term = make(&[2.0]);
        assert_eq!(energy(&term, &PLANAR).to_bits(), 0);
        assert_eq!(gradient(&term, &PLANAR, &[0.0; 12]), [0.0; 12]);
        let perpendicular = [2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 5.0];
        close(energy(&term, &perpendicular), 355.1392088649485);
        assert_eq!(gradient(&term, &perpendicular, &[0.0; 12]), [0.0; 12]);
    }
    #[test]
    fn mmff_oop_bend_plane_reversal_and_reflection() {
        let term = make(&[2.0]);
        let mut reversed = term.clone();
        reversed.at1_idxs[0] = 2;
        reversed.at3_idxs[0] = 0;
        close(energy(&reversed, &COORDS), 88.78480221623713);
        close_grad(&gradient(&reversed, &COORDS, &[0.0; 12]), &GRAD);
        let mut reflected = COORDS;
        reflected[11] = -5.0;
        let mut expected = GRAD;
        for axis in [2, 5, 8, 11] {
            expected[axis] = -expected[axis];
        }
        close(energy(&term, &reflected), 88.78480221623713);
        close_grad(&gradient(&term, &reflected, &[0.0; 12]), &expected);
        let p = |i| point_at(&COORDS, i);
        close(calc_oop_chi(&p(0), &p(1), &p(2), &p(3)), 45.0);
        close(calc_oop_chi(&p(2), &p(1), &p(0), &p(3)), -45.0);
    }
    #[test]
    fn mmff_oop_bend_translation_and_additive_nonzero_gradient() {
        let term = make(&[2.0]);
        let mut moved = COORDS;
        for row in moved.chunks_exact_mut(3) {
            row[0] += 8.0;
            row[1] -= 3.0;
            row[2] += 7.0;
        }
        let initial = [
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
        ];
        let expected: Vec<_> = initial.iter().zip(GRAD).map(|(a, g)| a + g).collect();
        close_grad(&gradient(&term, &moved, &initial), &expected);
        close(energy(&term, &moved), 88.78480221623713);
    }
    #[test]
    fn mmff_oop_bend_repeated_terms_and_source_total_association() {
        let term = make(&[2.0, 2.0]);
        close(energy(&term, &COORDS), 177.56960443247426);
        let expected: Vec<_> = GRAD.iter().map(|g| 2.0 * g).collect();
        close_grad(&gradient(&term, &COORDS, &[0.0; 12]), &expected);
        // Source per-term total association loses the first term before the
        // opposite huge pair cancels. Do not regroup across source terms.
        assert_eq!(energy(&make(&[2.0, 1e30, -1e30]), &COORDS).to_bits(), 0);
    }
    #[test]
    fn mmff_oop_bend_negative_koop_and_ignored_later_parameter_mutation() {
        let negative = make(&[-2.0]);
        close(energy(&negative, &COORDS), -88.78480221623713);
        close_grad(&gradient(&negative, &COORDS, &[0.0; 12]), &GRAD.map(|g| -g));
        let mut rows = [[0.0; 3]; 4];
        let mut field = ForceField::new(3);
        for row in &mut rows {
            field.positions_mut().push(row);
        }
        let mut term = OopBendContrib::new(&field);
        let mut params = MmffOop { koop: 2.0 };
        term.add_term(field.positions(), 0, 1, 2, 3, Some(&params))
            .unwrap();
        params.koop = f64::NAN;
        assert_eq!(term.koop, [2.0]);
        close(energy(&term, &COORDS), 88.78480221623713);
        term.add_term(field.positions(), 0, 1, 2, 3, Some(&params))
            .unwrap();
        assert!(term.koop[1].is_nan());
    }
    #[test]
    fn mmff_oop_bend_signed_zero_copy_and_gradient_addition() {
        for (k, negative_slots) in [
            (0.0, &[0, 1, 2, 6, 7, 8, 9, 10][..]),
            (-0.0, &[3, 4, 5, 11][..]),
        ] {
            let term = make(&[k]);
            assert_eq!(term.koop[0].to_bits(), k.to_bits());
            assert_eq!(term.clone().koop[0].to_bits(), k.to_bits());
            assert_eq!(energy(&term, &COORDS).to_bits(), 0);
            let g = gradient(&term, &COORDS, &[-0.0; 12]);
            for (i, value) in g.iter().enumerate() {
                assert_eq!(
                    value.to_bits(),
                    if negative_slots.contains(&i) {
                        0x8000000000000000
                    } else {
                        0
                    }
                );
            }
        }
    }
    #[test]
    fn mmff_oop_bend_each_zero_arm_energy_nan_gradient_unchanged() {
        let initial = [
            1.0, -0.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
        ];
        for k in [2.0, 0.0, f64::NAN] {
            let term = make(&[k]);
            for row in [0, 2, 3] {
                let mut coords = COORDS;
                coords[3 * row..3 * row + 3].fill(0.0);
                assert!(energy(&term, &coords).is_nan());
                assert_eq!(
                    gradient(&term, &coords, &initial)
                        .iter()
                        .map(|v| v.to_bits())
                        .collect::<Vec<_>>(),
                    initial.map(f64::to_bits)
                );
            }
            let mut coords = COORDS;
            coords[..3].fill(0.0);
            coords[9] = f64::NAN;
            assert_eq!(
                gradient(&term, &coords, &initial)
                    .iter()
                    .map(|v| v.to_bits())
                    .collect::<Vec<_>>(),
                initial.map(f64::to_bits)
            );
        }
    }
    #[test]
    fn mmff_oop_bend_strict_double_zero_threshold_neighbors() {
        let bits = 1e-10_f64.to_bits();
        let term = make(&[f64::NAN]);
        for (value, is_zero) in [
            (f64::from_bits(bits - 1), true),
            (1e-10, false),
            (f64::from_bits(bits + 1), false),
        ] {
            assert_eq!(is_double_zero(value), is_zero);
            assert_eq!(is_double_zero(-value), is_zero);
            assert_eq!((value * value).sqrt().to_bits(), value.to_bits());
            for row in [0, 2, 3] {
                let mut coords = COORDS;
                coords[3 * row..3 * row + 3].fill(0.0);
                coords[3 * row + if row == 2 { 1 } else { 0 }] = value;
                let g = gradient(&term, &coords, &[7.0; 12]);
                if is_zero {
                    assert_eq!(g, [7.0; 12]);
                } else {
                    assert!(g.iter().all(|v| v.is_nan()));
                }
            }
        }
        assert!(is_double_zero(-0.0));
        assert!(!is_double_zero(f64::NAN));
        assert!(!is_double_zero(f64::INFINITY));
    }
    #[test]
    fn mmff_oop_bend_nan_infinite_geometry() {
        for k in [2.0, 0.0] {
            let term = make(&[k]);
            for (axis, value) in [
                (0, f64::NAN),
                (0, f64::INFINITY),
                (7, f64::NAN),
                (11, f64::INFINITY),
            ] {
                let mut coords = COORDS;
                coords[axis] = value;
                assert!(energy(&term, &coords).is_nan());
                assert!(
                    gradient(&term, &coords, &[1.0; 12])
                        .iter()
                        .all(|v| v.is_nan())
                );
            }
        }
    }
    #[test]
    fn mmff_oop_bend_nan_infinite_koop_and_source_signs() {
        let nan = make(&[f64::NAN]);
        assert!(energy(&nan, &COORDS).is_nan());
        assert!(
            gradient(&nan, &COORDS, &[0.0; 12])
                .iter()
                .all(|v| v.is_nan())
        );
        for k in [f64::INFINITY, f64::NEG_INFINITY] {
            let term = make(&[k]);
            assert_eq!(energy(&term, &COORDS), k);
            let g = gradient(&term, &COORDS, &[0.0; 12]);
            for i in [0, 1, 6, 7] {
                assert!(g[i].is_nan());
            }
            for i in [2, 8, 9, 10] {
                assert_eq!(g[i], -k);
            }
            for i in [3, 4, 5, 11] {
                assert_eq!(g[i], k);
            }
        }
    }
    #[test]
    fn mmff_oop_bend_collinear_nonzero_arms_have_no_early_return() {
        for k in [0.0, 2.0] {
            let term = make(&[k]);
            for sign in [-1.0, 1.0] {
                let coords = [
                    2.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    sign * 3.0,
                    0.0,
                    0.0,
                    3.0,
                    4.0,
                    5.0,
                ];
                assert!(energy(&term, &coords).is_nan());
                assert!(
                    gradient(&term, &coords, &[1.0; 12])
                        .iter()
                        .all(|v| v.is_nan())
                );
            }
        }
    }
    #[test]
    fn mmff_oop_bend_clip_payload_and_conditional_cos_floor() {
        let payload = f64::from_bits(0x7ff8000000000042);
        let mut nan = payload;
        clip_to_one(&mut nan);
        assert_eq!(nan.to_bits(), payload.to_bits());
        assert_eq!(source_sin_theta(payload), 1e-8);
        assert_eq!(source_sin_theta(1.0), 1e-8);
        assert_eq!(source_sin_theta(-1.0), 1e-8);
        assert_eq!(source_sin_theta(0.0), 1.0);
        assert_eq!(
            source_sin_theta(f64::from_bits(1.0_f64.to_bits() + 1)),
            1e-8
        );
        for (mut value, expected) in [
            (f64::INFINITY, 1.0_f64.to_bits()),
            (f64::NEG_INFINITY, (-1.0_f64).to_bits()),
            (-0.0, 0x8000000000000000),
        ] {
            clip_to_one(&mut value);
            assert_eq!(value.to_bits(), expected);
        }
    }
    #[test]
    fn mmff_oop_bend_two_stage_theta_floor_and_nan_operand_order() {
        assert_eq!(source_theta_factors(1.0), (1e-8, 1e-4));
        assert_eq!(source_theta_factors(-1.0), (1e-8, 1e-4));
        assert_eq!(source_theta_factors(0.0), (1.0, 1.0));
        let (sq, sine) = source_theta_factors(f64::NAN);
        assert!(sq.is_nan());
        assert_eq!(sine, 1e-8);
        assert_eq!(
            source_theta_factors(f64::from_bits(1.0_f64.to_bits() + 1)),
            (1e-8, 1e-4)
        );
    }
    #[test]
    fn mmff_oop_bend_no_distance_cache_read_or_write() {
        let term = make(&[2.0]);
        let mut cache = [f64::NAN, -91.0, f64::INFINITY, -0.0];
        let bits = cache.map(f64::to_bits);
        {
            let mut ctx =
                EvaluationContext::for_oop_cache_preservation_test(&COORDS, &mut cache, 4);
            close(term.get_energy(&mut ctx).unwrap(), 88.78480221623713);
            let mut g = [0.0; 12];
            term.get_grad(&mut ctx, &mut g).unwrap();
            close_grad(&g, &GRAD);
        }
        assert_eq!(cache.map(f64::to_bits), bits);
        close(energy(&term, &COORDS), 88.78480221623713); // empty cache valid
    }
    #[test]
    fn mmff_oop_bend_signed32_storage_beyond65535() {
        let mut rows = vec![[0.0; 3]; 65539];
        let mut field = ForceField::new(3);
        for row in &mut rows {
            field.positions_mut().push(row);
        }
        let mut term = OopBendContrib::new(&field);
        append(&mut term, &field, [65536, 1, 65537, 65538], 2.0).unwrap();
        assert_eq!(
            [
                term.at1_idxs[0],
                term.at2_idxs[0],
                term.at3_idxs[0],
                term.at4_idxs[0]
            ],
            [65536, 1, 65537, 65538]
        );
        let mut coords = vec![0.0; 65539 * 3];
        coords[65536 * 3..65536 * 3 + 3].copy_from_slice(&COORDS[..3]);
        coords[65537 * 3..65537 * 3 + 3].copy_from_slice(&COORDS[6..9]);
        coords[65538 * 3..65538 * 3 + 3].copy_from_slice(&COORDS[9..12]);
        close(energy(&term, &coords), 88.78480221623713);
        let g = gradient(&term, &coords, &vec![0.0; coords.len()]);
        close_grad(&g[65536 * 3..65536 * 3 + 3], &GRAD[..3]);
        close_grad(&g[3..6], &GRAD[3..6]);
        close_grad(&g[65537 * 3..65537 * 3 + 3], &GRAD[6..9]);
        close_grad(&g[65538 * 3..65538 * 3 + 3], &GRAD[9..12]);
        assert_eq!(g[..3], [0.0; 3]);
    }
    #[test]
    fn mmff_oop_bend_raw32_narrowing_bijection_and_undefined_limits() {
        for (raw, stored) in [
            (0_u32, 0_i32),
            (65535, 65535),
            (65536, 65536),
            (2147483647, 2147483647),
            (2147483648, -2147483648),
            (u32::MAX, -1),
        ] {
            assert_eq!(raw as i32, stored);
        }
        let raw = [0_u32, 65536, 2147483648, u32::MAX];
        let stored = raw.map(|i| i as i32);
        for i in 0..4 {
            for j in i + 1..4 {
                assert_ne!(stored[i], stored[j]);
            }
        }
        // Raw-distinct u32 indices cannot collapse in i32. Negative offsets,
        // signed-address overflow and pointer aliasing have no source oracle.
        // Billion-row physical arrays are NOT allocated or claimed tested.
    }
    #[test]
    fn mmff_oop_bend_shared_signed32_term_count_scalar_boundaries() {
        for (n, end, count, last) in [
            (0_usize, 0_i32, 0, None),
            (1, 1, 1, Some(0)),
            (2147483647, 2147483647, 2147483647, Some(2147483646)),
            (2147483648, -2147483648, 0, None),
            (4294967295, -1, 0, None),
            (4294967296, 0, 0, None),
            (4294967297, 1, 1, Some(0)),
        ] {
            let range = source_term_indices(n);
            assert_eq!(range.start, 0);
            assert_eq!(range.end, end);
            assert_eq!(range.len(), count);
            assert_eq!(range.last(), last);
        }
    }
    #[test]
    fn mmff_oop_bend_deep_clone_and_trait_copy_all_five_arrays() {
        let term = make(&[2.0]);
        let copy = term.copy();
        let mut clone = term.clone();
        assert_ne!(term.at1_idxs.as_ptr(), clone.at1_idxs.as_ptr());
        assert_ne!(term.at2_idxs.as_ptr(), clone.at2_idxs.as_ptr());
        assert_ne!(term.at3_idxs.as_ptr(), clone.at3_idxs.as_ptr());
        assert_ne!(term.at4_idxs.as_ptr(), clone.at4_idxs.as_ptr());
        assert_ne!(term.koop.as_ptr(), clone.koop.as_ptr());
        clone.at1_idxs[0] = 2;
        clone.at2_idxs[0] = 3;
        clone.at3_idxs[0] = 0;
        clone.at4_idxs[0] = 1;
        clone.koop[0] = 4.0;
        assert_eq!(
            [
                term.at1_idxs[0],
                term.at2_idxs[0],
                term.at3_idxs[0],
                term.at4_idxs[0]
            ],
            [0, 1, 2, 3]
        );
        assert_eq!(term.koop, [2.0]);
        // All five cloned arrays were changed independently above. Restore
        // the central/lateral orientation for the separate fixed-value check.
        clone.at2_idxs[0] = 1;
        clone.at4_idxs[0] = 3;
        close(energy(&clone, &COORDS), 177.56960443247426);
        close(energy(&*copy, &COORDS), 88.78480221623713);
        close_grad(&gradient(&*copy, &COORDS, &[0.0; 12]), &GRAD);
    }
    #[test]
    fn mmff_oop_bend_get_single_gradient_selects_complete_term() {
        let term = make(&[2.0, 4.0]);
        let original_coords = COORDS;
        for (idx, factor) in [(0, 1.0), (1, 2.0)] {
            let mut g = [7.0; 12];
            term.get_single_grad(&COORDS, &mut g, idx).unwrap();
            let expected: Vec<_> = GRAD.iter().map(|v| 7.0 + factor * v).collect();
            close_grad(&g, &expected);
        }
        assert_eq!(COORDS, original_coords);
    }
    #[test]
    fn mmff_oop_bend_real_field_initialization_copy_and_rebind() {
        let mut rows = [[2.0, 0.0, 0.0], [0.0; 3], [0.0, 3.0, 0.0], [3.0, 4.0, 5.0]];
        let mut field = ForceField::new(3);
        for row in &mut rows {
            field.positions_mut().push(row);
        }
        let mut term = OopBendContrib::new(&field);
        append(&mut term, &field, [0, 1, 2, 3], 2.0).unwrap();
        field.add_contribution(Box::new(term));
        assert_eq!(
            cf3d_bld_b05_calc_energy(&mut field, &COORDS),
            Err(ForceFieldKernelError::NotInitialized)
        );
        field.initialize().unwrap();
        close(field.calc_energy_current(None).unwrap(), 88.78480221623713);
        let mut copy = field.copy();
        assert!(copy.positions().is_empty());
        assert_eq!(
            cf3d_bld_b05_calc_energy(&mut copy, &COORDS),
            Err(ForceFieldKernelError::NotInitialized)
        );
        let mut copied_rows = rows_for_copy();
        for row in &mut copied_rows {
            copy.positions_mut().push(row);
        }
        copy.initialize().unwrap();
        close(
            cf3d_bld_b05_calc_energy(&mut copy, &COORDS).unwrap(),
            88.78480221623713,
        );
        let mut g = [0.0; 12];
        cf3d_bld_b05_calc_grad(&mut copy, &COORDS, &mut g).unwrap();
        close_grad(&g, &GRAD);
        assert_eq!(cf3d_bld_b05_calc_energy(&mut copy, &PLANAR), Ok(0.0));
        close(field.calc_energy_current(None).unwrap(), 88.78480221623713);
    }
    fn rows_for_copy() -> [[f64; 3]; 4] {
        [[2.0, 0.0, 0.0], [0.0; 3], [0.0, 3.0, 0.0], [3.0, 4.0, 5.0]]
    }
    #[test]
    fn mmff_oop_bend_physical_coordinates_ignore_stale_distances() {
        let mut rows = rows_for_copy();
        let mut field = ForceField::new(3);
        for row in &mut rows {
            field.positions_mut().push(row);
        }
        let mut term = OopBendContrib::new(&field);
        append(&mut term, &field, [0, 1, 2, 3], 2.0).unwrap();
        field.add_contribution(Box::new(term));
        field.initialize().unwrap();
        // Prime a real stale J-L distance through the accepted private Bond
        // contribution with zero force constant; its energy/gradient are0.
        let mut bond = super::super::bond_stretch::BondStretchContrib::new(&field);
        bond.add_term(
            field.positions(),
            1,
            3,
            &super::super::params::MmffBond { r0: 1.0, kb: 0.0 },
        )
        .unwrap();
        field.add_contribution(Box::new(bond));
        close(field.calc_energy_current(None).unwrap(), 88.78480221623713);
        let mut g = [0.0; 12];
        cf3d_bld_b05_calc_grad(&mut field, &PLANAR, &mut g).unwrap();
        assert_eq!(g, [0.0; 12]);
        cf3d_bld_b05_calc_grad(&mut field, &COORDS, &mut g).unwrap();
        close_grad(&g, &GRAD);
        close(field.calc_energy_current(None).unwrap(), 88.78480221623713);
    }
    #[test]
    fn mmff_oop_bend_native_safe_address_scalar_policy() {
        // PROPOSED NATIVE SAFETY POLICY, pending ROOT; not C++-defined errors.
        assert_eq!(checked_address(0, 3), Ok(0));
        assert_eq!(checked_address(1, 6), Ok(3));
        assert_eq!(checked_address(715827881, usize::MAX), Ok(2147483643));
        for (index, len) in [
            (-1, usize::MAX),
            (i32::MIN, usize::MAX),
            (715827882, usize::MAX),
            (i32::MAX, usize::MAX),
            (0, 2),
            (1, 5),
        ] {
            assert_eq!(
                checked_address(index, len),
                Err(ForceFieldKernelError::OopAddressOutsideDefinedSource {
                    index,
                    buffer_len: len
                })
            );
        }
    }
    #[test]
    fn mmff_oop_bend_native_safe_short_buffers_and_partial_updates() {
        // PROPOSED NATIVE SAFETY POLICY, pending ROOT; no source numerical
        // oracle for short buffers or corrupted negative stored endpoint.
        let term = make(&[2.0]);
        let mut short = [7.0; 11];
        assert_eq!(
            term.get_single_grad(&COORDS, &mut short, 0),
            Err(ForceFieldKernelError::OopAddressOutsideDefinedSource {
                index: 3,
                buffer_len: 11
            })
        );
        assert_eq!(short, [7.0; 11]);
        let mut g = [7.0; 12];
        assert_eq!(
            term.get_single_grad(&COORDS[..11], &mut g, 0),
            Err(ForceFieldKernelError::OopAddressOutsideDefinedSource {
                index: 3,
                buffer_len: 11
            })
        );
        assert_eq!(g, [7.0; 12]);
        let mut partial = make(&[2.0, 4.0]);
        partial.at4_idxs[1] = -1;
        let mut unused_cache = [];
        let mut ctx = EvaluationContext::for_test(&COORDS, &mut unused_cache, 4);
        let mut g = [0.0; 12];
        assert_eq!(
            partial.get_grad(&mut ctx, &mut g),
            Err(ForceFieldKernelError::OopAddressOutsideDefinedSource {
                index: -1,
                buffer_len: 12
            })
        );
        close_grad(&g, &GRAD);
        assert_eq!(
            partial.get_energy(&mut ctx),
            Err(ForceFieldKernelError::OopAddressOutsideDefinedSource {
                index: -1,
                buffer_len: 12
            })
        );
    }
}
