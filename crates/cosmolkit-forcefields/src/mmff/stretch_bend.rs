//! Private MMFF stretch-bend contribution.
//!
//! Source: RDKit351f8f378f8ad6bbd517980c38896e66bf907af8,
//! Code/ForceField/MMFF/StretchBend.{cpp,h}, BSD license.
//! Source cases validate the modeled private 3D boundary; independent acceptance is pending.
//! Native signed-width controls are bounded; huge vectors and public APIs are unvalidated.

use super::{
    angle_bend::point_at,
    bond_stretch::source_term_indices,
    numerical::{
        calc_angle_rest_value, calc_bond_rest_length, calc_cos_theta, calc_stbn_force_constants,
        calc_stretch_bend_energy,
    },
    params::{MmffAngle, MmffBond, MmffStbn},
};
use crate::{
    geometry::Point3,
    kernel::{
        AngleIndexArgument, EvaluationContext, ForceField, ForceFieldContribution,
        ForceFieldKernelError,
    },
    uff::params::{DEG2RAD, RAD2DEG, clip_to_one},
};

// RDKit✔️✔️:   std::vector<int16_t> d_at1Idxs;
// RDKit✔️✔️:   std::vector<int16_t> d_at2Idxs;
// RDKit✔️✔️:   std::vector<int16_t> d_at3Idxs;
// RDKit✔️✔️:   std::vector<double> d_restLen1s;
// RDKit✔️✔️:   std::vector<double> d_restLen2s;
// RDKit✔️✔️:   std::vector<double> d_theta0s;
// RDKit✔️✔️:   std::vector<double> d_forceConstants1;
// RDKit✔️✔️:   std::vector<double> d_forceConstants2;
#[derive(Clone, Debug)]
pub(super) struct StretchBendContrib {
    at1_idxs: Vec<i16>,
    at2_idxs: Vec<i16>,
    at3_idxs: Vec<i16>,
    rest_len1s: Vec<f64>,
    rest_len2s: Vec<f64>,
    theta0s: Vec<f64>,
    force_constants1: Vec<f64>,
    force_constants2: Vec<f64>,
}

fn source_sin_theta(cos_theta: f64) -> f64 {
    // RDKit✔️✔️:     double sinThetaSq = 1.0 - cosTheta * cosTheta;
    // RDKit✔️✔️:     double sinTheta = std::max(sqrt(sinThetaSq), 1.0e-8);
    // Behavior: direct sqrt preserves NaN as std::max's FIRST operand.
    // This differs from AngleBend's conditional square treatment; no f64::max.
    // Complexity: O(1) scalar work, no allocation or copy.
    let sin_theta_sq = 1.0 - cos_theta * cos_theta;
    let candidate = sin_theta_sq.sqrt();
    if candidate < 1e-8 { 1e-8 } else { candidate }
}

impl StretchBendContrib {
    pub(super) fn new(_owner: &ForceField<'_>) -> Self {
        // RDKit✔️✔️: StretchBendContrib::StretchBendContrib(ForceField *owner) {
        // RDKit✔️✔️:   PRECONDITION(owner, "bad owner");
        // RDKit✔️✔️:   dp_forceField = owner;
        // RDKit✔️✔️: }
        // Behavior: owner borrow is non-null; the existing owner-derived
        // context replaces the source pointer. Ownerless default not evaluable.
        // Complexity: eight empty arrays, O(1), no allocation or owner clone.
        Self {
            at1_idxs: Vec::new(),
            at2_idxs: Vec::new(),
            at3_idxs: Vec::new(),
            rest_len1s: Vec::new(),
            rest_len2s: Vec::new(),
            theta0s: Vec::new(),
            force_constants1: Vec::new(),
            force_constants2: Vec::new(),
        }
    }
    pub(super) fn add_term(
        &mut self,
        positions: &[&mut [f64]],
        idx1: u32,
        idx2: u32,
        idx3: u32,
        stbn: &MmffStbn,
        angle: &MmffAngle,
        bond1: &MmffBond,
        bond2: &MmffBond,
    ) -> Result<(), ForceFieldKernelError> {
        // RDKit✔️✔️: void StretchBendContrib::addTerm(
        // RDKit✔️✔️:     const unsigned int idx1, const unsigned int idx2, const unsigned int idx3,
        // RDKit✔️✔️:     const ForceFields::MMFF::MMFFStbn *mmffStbnParams,
        // RDKit✔️✔️:     const ForceFields::MMFF::MMFFAngle *mmffAngleParams,
        // RDKit✔️✔️:     const ForceFields::MMFF::MMFFBond *mmffBondParams1,
        // RDKit✔️✔️:     const ForceFields::MMFF::MMFFBond *mmffBondParams2) {
        // RDKit✔️✔️:   PRECONDITION(((idx1 != idx2) && (idx2 != idx3) && (idx1 != idx3)),
        // RDKit✔️✔️:                "degenerate points");
        // RDKit✔️✔️:   URANGE_CHECK(idx1, dp_forceField->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx2, dp_forceField->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx3, dp_forceField->positions().size());
        // RDKit✔️✔️:
        // RDKit✔️✔️:   d_at1Idxs.push_back(idx1);
        // RDKit✔️✔️:   d_at2Idxs.push_back(idx2);
        // RDKit✔️✔️:   d_at3Idxs.push_back(idx3);
        // RDKit✔️✔️:   d_restLen1s.push_back(Utils::calcBondRestLength(mmffBondParams1));
        // RDKit✔️✔️:   d_restLen2s.push_back(Utils::calcBondRestLength(mmffBondParams2));
        // RDKit✔️✔️:   d_theta0s.push_back(Utils::calcAngleRestValue(mmffAngleParams));
        // RDKit✔️✔️:   std::pair<double, double> forceConstants =
        // RDKit✔️✔️:       Utils::calcStbnForceConstants(mmffStbnParams);
        // RDKit✔️✔️:   d_forceConstants1.push_back(forceConstants.first);
        // RDKit✔️✔️:   d_forceConstants2.push_back(forceConstants.second);
        // RDKit✔️✔️: }
        // Behavior: source DISTINCT BEFORE all raw bounds, before any append.
        // Three source int16_t narrows are modulo2^16 with signed interpretation.
        // Non-null parameter borrows allow reuse of the sole existing helpers;
        // copy five scalar values in source order, no new numeric guards.
        // Complexity: fixed checks, eight amortizedO(1) indexed-array appends;
        // same3i16/5f64 separate storage, no per-term owner/parameter clone.
        // Validation: frozen distinct/bounds, copied IEEE values, width and growth cases pass.
        if idx1 == idx2 || idx2 == idx3 || idx1 == idx3 {
            return Err(ForceFieldKernelError::AngleDegeneratePoints);
        }
        for (argument, index) in [
            (AngleIndexArgument::First, idx1),
            (AngleIndexArgument::Second, idx2),
            (AngleIndexArgument::Third, idx3),
        ] {
            if index as usize >= positions.len() {
                return Err(ForceFieldKernelError::AngleIndexOutOfRange {
                    argument,
                    index,
                    upper_bound: positions.len(),
                });
            }
        }
        self.at1_idxs.push(idx1 as i16);
        self.at2_idxs.push(idx2 as i16);
        self.at3_idxs.push(idx3 as i16);
        self.rest_len1s.push(calc_bond_rest_length(bond1));
        self.rest_len2s.push(calc_bond_rest_length(bond2));
        self.theta0s.push(calc_angle_rest_value(angle));
        let force_constants = calc_stbn_force_constants(stbn);
        self.force_constants1.push(force_constants.0);
        self.force_constants2.push(force_constants.1);
        Ok(())
    }
}

impl ForceFieldContribution for StretchBendContrib {
    fn get_energy(
        &self,
        context: &mut EvaluationContext<'_>,
    ) -> Result<f64, ForceFieldKernelError> {
        // RDKit✔️✔️: double StretchBendContrib::getEnergy(double *pos) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:
        // RDKit✔️✔️:   double totalEnergy = 0.0;
        // RDKit✔️✔️:   const int numTerms = d_at1Idxs.size();
        // RDKit✔️✔️:   for (int i = 0; i < numTerms; i++) {
        // RDKit✔️✔️:     const int16_t at1Idx = d_at1Idxs[i];
        // RDKit✔️✔️:     const int16_t at2Idx = d_at2Idxs[i];
        // RDKit✔️✔️:     const int16_t at3Idx = d_at3Idxs[i];
        // RDKit✔️✔️:
        // RDKit✔️✔️:     double dist1 = dp_forceField->distance(at1Idx, at2Idx, pos);
        // RDKit✔️✔️:     double dist2 = dp_forceField->distance(at2Idx, at3Idx, pos);
        // RDKit✔️✔️:
        // RDKit✔️✔️:     RDGeom::Point3D p1(pos[3 * at1Idx], pos[3 * at1Idx + 1],
        // RDKit✔️✔️:                        pos[3 * at1Idx + 2]);
        // RDKit✔️✔️:     RDGeom::Point3D p2(pos[3 * at2Idx], pos[3 * at2Idx + 1],
        // RDKit✔️✔️:                        pos[3 * at2Idx + 2]);
        // RDKit✔️✔️:     RDGeom::Point3D p3(pos[3 * at3Idx], pos[3 * at3Idx + 1],
        // RDKit✔️✔️:                        pos[3 * at3Idx + 2]);
        // RDKit✔️✔️:     std::pair<double, double> forceConstantsPair =
        // RDKit✔️✔️:         std::make_pair(d_forceConstants1[i], d_forceConstants2[i]);
        // RDKit✔️✔️:     std::pair<double, double> stretchBendEnergies =
        // RDKit✔️✔️:         Utils::calcStretchBendEnergy(
        // RDKit✔️✔️:             dist1 - d_restLen1s[i], dist2 - d_restLen2s[i],
        // RDKit✔️✔️:             RAD2DEG * acos(Utils::calcCosTheta(p1, p2, p3, dist1, dist2)) -
        // RDKit✔️✔️:                 d_theta0s[i],
        // RDKit✔️✔️:             forceConstantsPair);
        // RDKit✔️✔️:     totalEnergy += (stretchBendEnergies.first + stretchBendEnergies.second);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return totalEnergy;
        // RDKit✔️✔️: }
        // Behavior: shared signed32 loop and signed16 endpoints preserve CXX20
        // conversion. BOTH ordered distances precede coordinate addressing;
        // cached distance errors propagate before any invalid position access.
        // Pair energies are summed before the total, preserving association.
        // Complexity: O(max(0,signed32(terms))) / O(1) stack temporary space;
        // two indexed cached distances, no allocation, buffering or state clone.
        // Validation: fixed energies, pair grouping, width/cache/error cases pass.
        let mut total_energy = 0.0;
        for i in source_term_indices(self.at1_idxs.len()) {
            let i = i as usize;
            let at1 = self.at1_idxs[i];
            let at2 = self.at2_idxs[i];
            let at3 = self.at3_idxs[i];
            let dist1 = context.distance(at1 as u32, at2 as u32)?;
            let dist2 = context.distance(at2 as u32, at3 as u32)?;
            let coords = context.coordinates();
            let p1 = point_at(coords, at1 as usize);
            let p2 = point_at(coords, at2 as usize);
            let p3 = point_at(coords, at3 as usize);
            let pair = calc_stretch_bend_energy(
                dist1 - self.rest_len1s[i],
                dist2 - self.rest_len2s[i],
                RAD2DEG * calc_cos_theta(p1, p2, p3, dist1, dist2).acos() - self.theta0s[i],
                (self.force_constants1[i], self.force_constants2[i]),
            );
            total_energy += pair.0 + pair.1;
        }
        Ok(total_energy)
    }
    fn get_grad(
        &self,
        context: &mut EvaluationContext<'_>,
        gradient: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        // RDKit✔️✔️: void StretchBendContrib::getGrad(double *pos, double *grad) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:   PRECONDITION(grad, "bad vector");
        // RDKit✔️✔️:
        // RDKit✔️✔️:   const int numTerms = d_at1Idxs.size();
        // RDKit✔️✔️:   for (int i = 0; i < numTerms; ++i) {
        // RDKit✔️✔️:     const int16_t at1Idx = d_at1Idxs[i];
        // RDKit✔️✔️:     const int16_t at2Idx = d_at2Idxs[i];
        // RDKit✔️✔️:     const int16_t at3Idx = d_at3Idxs[i];
        // RDKit✔️✔️:     const double theta0 = d_theta0s[i];
        // RDKit✔️✔️:     const double forceConstant1 = d_forceConstants1[i];
        // RDKit✔️✔️:     const double forceConstant2 = d_forceConstants2[i];
        // RDKit✔️✔️:     const double restLen1 = d_restLen1s[i];
        // RDKit✔️✔️:     const double restLen2 = d_restLen2s[i];
        // RDKit✔️✔️:
        // RDKit✔️✔️:     double dist1 = dp_forceField->distance(at1Idx, at2Idx, pos);
        // RDKit✔️✔️:     double dist2 = dp_forceField->distance(at2Idx, at3Idx, pos);
        // RDKit✔️✔️:
        // RDKit✔️✔️:     RDGeom::Point3D p1(pos[3 * at1Idx], pos[3 * at1Idx + 1],
        // RDKit✔️✔️:                        pos[3 * at1Idx + 2]);
        // RDKit✔️✔️:     RDGeom::Point3D p2(pos[3 * at2Idx], pos[3 * at2Idx + 1],
        // RDKit✔️✔️:                        pos[3 * at2Idx + 2]);
        // RDKit✔️✔️:     RDGeom::Point3D p3(pos[3 * at3Idx], pos[3 * at3Idx + 1],
        // RDKit✔️✔️:                        pos[3 * at3Idx + 2]);
        // RDKit✔️✔️:     double *g1 = &(grad[3 * at1Idx]);
        // RDKit✔️✔️:     double *g2 = &(grad[3 * at2Idx]);
        // RDKit✔️✔️:     double *g3 = &(grad[3 * at3Idx]);
        // RDKit✔️✔️:
        // RDKit✔️✔️:     RDGeom::Point3D p12 = (p1 - p2) / dist1;
        // RDKit✔️✔️:     RDGeom::Point3D p32 = (p3 - p2) / dist2;
        // RDKit✔️✔️:     double const c5 = MDYNE_A_TO_KCAL_MOL * DEG2RAD;
        // RDKit✔️✔️:     double cosTheta = p12.dotProduct(p32);
        // RDKit✔️✔️:     clipToOne(cosTheta);
        // RDKit✔️✔️:     double sinThetaSq = 1.0 - cosTheta * cosTheta;
        // RDKit✔️✔️:     double sinTheta = std::max(sqrt(sinThetaSq), 1.0e-8);
        // RDKit✔️✔️:     double angleTerm = RAD2DEG * acos(cosTheta) - theta0;
        // RDKit✔️✔️:     double distTerm = RAD2DEG * (forceConstant1 * (dist1 - restLen1) +
        // RDKit✔️✔️:                                  forceConstant2 * (dist2 - restLen2));
        // RDKit✔️✔️:     double dCos_dS1 = 1.0 / dist1 * (p32.x - cosTheta * p12.x);
        // RDKit✔️✔️:     double dCos_dS2 = 1.0 / dist1 * (p32.y - cosTheta * p12.y);
        // RDKit✔️✔️:     double dCos_dS3 = 1.0 / dist1 * (p32.z - cosTheta * p12.z);
        // RDKit✔️✔️:
        // RDKit✔️✔️:     double dCos_dS4 = 1.0 / dist2 * (p12.x - cosTheta * p32.x);
        // RDKit✔️✔️:     double dCos_dS5 = 1.0 / dist2 * (p12.y - cosTheta * p32.y);
        // RDKit✔️✔️:     double dCos_dS6 = 1.0 / dist2 * (p12.z - cosTheta * p32.z);
        // RDKit✔️✔️:
        // RDKit✔️✔️:     g1[0] += c5 * (p12.x * forceConstant1 * angleTerm +
        // RDKit✔️✔️:                    dCos_dS1 / (-sinTheta) * distTerm);
        // RDKit✔️✔️:     g1[1] += c5 * (p12.y * forceConstant1 * angleTerm +
        // RDKit✔️✔️:                    dCos_dS2 / (-sinTheta) * distTerm);
        // RDKit✔️✔️:     g1[2] += c5 * (p12.z * forceConstant1 * angleTerm +
        // RDKit✔️✔️:                    dCos_dS3 / (-sinTheta) * distTerm);
        // RDKit✔️✔️:
        // RDKit✔️✔️:     g2[0] +=
        // RDKit✔️✔️:         c5 * ((-p12.x * forceConstant1 - p32.x * forceConstant2) * angleTerm +
        // RDKit✔️✔️:               (-dCos_dS1 - dCos_dS4) / (-sinTheta) * distTerm);
        // RDKit✔️✔️:     g2[1] +=
        // RDKit✔️✔️:         c5 * ((-p12.y * forceConstant1 - p32.y * forceConstant2) * angleTerm +
        // RDKit✔️✔️:               (-dCos_dS2 - dCos_dS5) / (-sinTheta) * distTerm);
        // RDKit✔️✔️:     g2[2] +=
        // RDKit✔️✔️:         c5 * ((-p12.z * forceConstant1 - p32.z * forceConstant2) * angleTerm +
        // RDKit✔️✔️:               (-dCos_dS3 - dCos_dS6) / (-sinTheta) * distTerm);
        // RDKit✔️✔️:
        // RDKit✔️✔️:     g3[0] += c5 * (p32.x * forceConstant2 * angleTerm +
        // RDKit✔️✔️:                    dCos_dS4 / (-sinTheta) * distTerm);
        // RDKit✔️✔️:     g3[1] += c5 * (p32.y * forceConstant2 * angleTerm +
        // RDKit✔️✔️:                    dCos_dS5 / (-sinTheta) * distTerm);
        // RDKit✔️✔️:     g3[2] += c5 * (p32.z * forceConstant2 * angleTerm +
        // RDKit✔️✔️:                    dCos_dS6 / (-sinTheta) * distTerm);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // RDKit✔️✔️: constexpr double MDYNE_A_TO_KCAL_MOL = 143.9325;
        // Behavior: retain distinct source gradient body, operand association
        // and all nine ordered additions. Reuse Point3 scalar operators and
        // address helper only; AngleBend sin/gradient define different behavior.
        // Both checked distances precede all position/gradient addresses;
        // aliased narrowed endpoints update the same rows in source order.
        // Valid3D position/gradient storage is the unchanged kernel contract.
        // Complexity: O(max(0,signed32(terms))) / O(1) stack temporaries,
        // no allocation, copies, scans, extra registry or gradient buffering.
        // Validation: all9 fixed writes, IEEE/alias/partial-cache lifecycle cases pass.
        for i in source_term_indices(self.at1_idxs.len()) {
            let i = i as usize;
            let at1 = self.at1_idxs[i];
            let at2 = self.at2_idxs[i];
            let at3 = self.at3_idxs[i];
            let theta0 = self.theta0s[i];
            let force_constant1 = self.force_constants1[i];
            let force_constant2 = self.force_constants2[i];
            let rest_len1 = self.rest_len1s[i];
            let rest_len2 = self.rest_len2s[i];
            let dist1 = context.distance(at1 as u32, at2 as u32)?;
            let dist2 = context.distance(at2 as u32, at3 as u32)?;
            let coords = context.coordinates();
            let p1 = point_at(coords, at1 as usize);
            let p2 = point_at(coords, at2 as usize);
            let p3 = point_at(coords, at3 as usize);
            let g1 = 3 * (at1 as usize);
            let g2 = 3 * (at2 as usize);
            let g3 = 3 * (at3 as usize);
            let p12 = Point3::difference(&p1, &p2).divided(dist1);
            let p32 = Point3::difference(&p3, &p2).divided(dist2);
            let c5 = 143.9325 * DEG2RAD;
            let mut cos_theta = p12.dot_product(&p32);
            clip_to_one(&mut cos_theta);
            let sin_theta = source_sin_theta(cos_theta);
            let angle_term = RAD2DEG * cos_theta.acos() - theta0;
            let dist_term = RAD2DEG
                * (force_constant1 * (dist1 - rest_len1) + force_constant2 * (dist2 - rest_len2));
            let dcos_ds1 = 1.0 / dist1 * (p32.x - cos_theta * p12.x);
            let dcos_ds2 = 1.0 / dist1 * (p32.y - cos_theta * p12.y);
            let dcos_ds3 = 1.0 / dist1 * (p32.z - cos_theta * p12.z);
            let dcos_ds4 = 1.0 / dist2 * (p12.x - cos_theta * p32.x);
            let dcos_ds5 = 1.0 / dist2 * (p12.y - cos_theta * p32.y);
            let dcos_ds6 = 1.0 / dist2 * (p12.z - cos_theta * p32.z);
            gradient[g1] +=
                c5 * (p12.x * force_constant1 * angle_term + dcos_ds1 / (-sin_theta) * dist_term);
            gradient[g1 + 1] +=
                c5 * (p12.y * force_constant1 * angle_term + dcos_ds2 / (-sin_theta) * dist_term);
            gradient[g1 + 2] +=
                c5 * (p12.z * force_constant1 * angle_term + dcos_ds3 / (-sin_theta) * dist_term);
            gradient[g2] += c5
                * ((-p12.x * force_constant1 - p32.x * force_constant2) * angle_term
                    + (-dcos_ds1 - dcos_ds4) / (-sin_theta) * dist_term);
            gradient[g2 + 1] += c5
                * ((-p12.y * force_constant1 - p32.y * force_constant2) * angle_term
                    + (-dcos_ds2 - dcos_ds5) / (-sin_theta) * dist_term);
            gradient[g2 + 2] += c5
                * ((-p12.z * force_constant1 - p32.z * force_constant2) * angle_term
                    + (-dcos_ds3 - dcos_ds6) / (-sin_theta) * dist_term);
            gradient[g3] +=
                c5 * (p32.x * force_constant2 * angle_term + dcos_ds4 / (-sin_theta) * dist_term);
            gradient[g3 + 1] +=
                c5 * (p32.y * force_constant2 * angle_term + dcos_ds5 / (-sin_theta) * dist_term);
            gradient[g3 + 2] +=
                c5 * (p32.z * force_constant2 * angle_term + dcos_ds6 / (-sin_theta) * dist_term);
        }
        Ok(())
    }
    fn copy(&self) -> Box<dyn ForceFieldContribution> {
        // RDKit✔️✔️:   StretchBendContrib *copy() const override {
        // RDKit✔️✔️:     return new StretchBendContrib(*this);
        // RDKit✔️✔️:   }
        // Behavior: independently clone all eight arrays; existing ForceField
        // rebind/context lifecycle owns initialization and coordinate borrows.
        // Complexity: O(terms) copied elements, one allocation per populated
        // array; no owner/coordinate clone or shared mutable contribution data.
        // Validation: all8 independent arrays and real field reset/rebind cases pass.
        Box::new(self.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernel::{
        ForceFieldIndexArgument, cf3d_bld_b05_calc_energy, cf3d_bld_b05_calc_grad,
    };

    const ORTHOGONAL: [f64; 9] = [2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0];
    const UNIT: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
    // Independently reduced pinned-source equations: angle delta30, arm
    // deltas1,2, asymmetric constants2,5, c5=143.9325*pi/180.
    const FIXED_GRAD: [f64; 9] = [
        150.7257615376043,
        -863.595,
        0.0,
        425.00423846239573,
        486.78059615598926,
        0.0,
        -575.73,
        376.81440384401077,
        0.0,
    ];

    #[derive(Clone, Copy)]
    struct Spec {
        idx: [u32; 3],
        r0: [f64; 2],
        theta0: f64,
        k: [f64; 2],
    }
    fn spec(r1: f64, r2: f64, theta0: f64, k1: f64, k2: f64) -> Spec {
        Spec {
            idx: [0, 1, 2],
            r0: [r1, r2],
            theta0,
            k: [k1, k2],
        }
    }
    fn base() -> Spec {
        spec(1.0, 1.0, 60.0, 2.0, 5.0)
    }
    fn append(
        term: &mut StretchBendContrib,
        field: &ForceField<'_>,
        s: Spec,
    ) -> Result<(), ForceFieldKernelError> {
        term.add_term(
            field.positions(),
            s.idx[0],
            s.idx[1],
            s.idx[2],
            &MmffStbn {
                kba_ijk: s.k[0],
                kba_kji: s.k[1],
            },
            &MmffAngle {
                theta0: s.theta0,
                ka: f64::NAN,
            },
            &MmffBond {
                r0: s.r0[0],
                kb: f64::NAN,
            },
            &MmffBond {
                r0: s.r0[1],
                kb: f64::INFINITY,
            },
        )
    }
    fn contribution(specs: &[Spec]) -> StretchBendContrib {
        let count = specs
            .iter()
            .flat_map(|s| s.idx)
            .max()
            .map_or(0, |x| x as usize + 1);
        let mut rows = vec![[0.0; 3]; count];
        let mut field = ForceField::new(3);
        for row in &mut rows {
            field.positions_mut().push(row);
        }
        let mut term = StretchBendContrib::new(&field);
        for &s in specs {
            append(&mut term, &field, s).unwrap();
        }
        term
    }
    fn energy(term: &dyn ForceFieldContribution, coords: &[f64]) -> f64 {
        let n = coords.len() / 3;
        let mut cache = vec![-1.0; n * (n + 1) / 2];
        term.get_energy(&mut EvaluationContext::for_test(
            coords, &mut cache, n as u32,
        ))
        .unwrap()
    }
    fn grad(term: &dyn ForceFieldContribution, coords: &[f64], initial: &[f64]) -> Vec<f64> {
        let n = coords.len() / 3;
        let mut cache = vec![-1.0; n * (n + 1) / 2];
        let mut g = initial.to_vec();
        term.get_grad(
            &mut EvaluationContext::for_test(coords, &mut cache, n as u32),
            &mut g,
        )
        .unwrap();
        g
    }
    fn close(a: f64, e: f64) {
        assert!((a - e).abs() < 1e-10, "{a:?} != {e:?}");
    }
    fn close_grad(a: &[f64], e: &[f64]) {
        assert_eq!(a.len(), e.len());
        for (&a, &e) in a.iter().zip(e) {
            close(a, e);
        }
    }

    #[test]
    fn mmff_stretch_bend_empty_owner_and_copy() {
        let term = StretchBendContrib::new(&ForceField::new(3));
        assert_eq!(
            [
                term.at1_idxs.capacity(),
                term.at2_idxs.capacity(),
                term.at3_idxs.capacity(),
                term.rest_len1s.capacity(),
                term.rest_len2s.capacity(),
                term.theta0s.capacity(),
                term.force_constants1.capacity(),
                term.force_constants2.capacity()
            ],
            [0; 8]
        );
        assert_eq!(energy(&term, &[]).to_bits(), 0);
        assert_eq!(grad(&term, &[], &[1.0, -2.0, 3.0]), [1.0, -2.0, 3.0]);
        assert_eq!(energy(&term.clone(), &[]).to_bits(), 0);
        assert_eq!(energy(&*term.copy(), &[]).to_bits(), 0);
    }
    #[test]
    fn mmff_stretch_bend_distinct_before_raw_bounds() {
        let mut rows = [[0.0; 3]; 3];
        let mut extra = [0.0; 3];
        let mut field = ForceField::new(3);
        for row in &mut rows {
            field.positions_mut().push(row);
        }
        let mut term = StretchBendContrib::new(&field);
        for (idx, expected) in [
            (
                [u32::MAX, u32::MAX, 2],
                ForceFieldKernelError::AngleDegeneratePoints,
            ),
            (
                [0, u32::MAX, u32::MAX],
                ForceFieldKernelError::AngleDegeneratePoints,
            ),
            (
                [u32::MAX, 1, u32::MAX],
                ForceFieldKernelError::AngleDegeneratePoints,
            ),
            (
                [0, 0, u32::MAX],
                ForceFieldKernelError::AngleDegeneratePoints,
            ),
            (
                [3, 4, 5],
                ForceFieldKernelError::AngleIndexOutOfRange {
                    argument: AngleIndexArgument::First,
                    index: 3,
                    upper_bound: 3,
                },
            ),
            (
                [0, 3, 4],
                ForceFieldKernelError::AngleIndexOutOfRange {
                    argument: AngleIndexArgument::Second,
                    index: 3,
                    upper_bound: 3,
                },
            ),
            (
                [0, 1, u32::MAX],
                ForceFieldKernelError::AngleIndexOutOfRange {
                    argument: AngleIndexArgument::Third,
                    index: u32::MAX,
                    upper_bound: 3,
                },
            ),
        ] {
            let mut s = base();
            s.idx = idx;
            assert_eq!(append(&mut term, &field, s), Err(expected));
            assert_eq!(
                [
                    term.at1_idxs.len(),
                    term.at2_idxs.len(),
                    term.at3_idxs.len(),
                    term.rest_len1s.len(),
                    term.rest_len2s.len(),
                    term.theta0s.len(),
                    term.force_constants1.len(),
                    term.force_constants2.len()
                ],
                [0; 8]
            );
        }
        append(&mut term, &field, base()).unwrap();
        let old = format!("{term:?}");
        let mut s = base();
        s.idx = [0, 1, 3];
        assert!(append(&mut term, &field, s).is_err());
        assert_eq!(format!("{term:?}"), old);
        field.positions_mut().push(&mut extra);
        append(&mut term, &field, s).unwrap();
        assert_eq!(term.at3_idxs, [2, 3]);
    }
    #[test]
    fn mmff_stretch_bend_parameter_values_copied_and_unused_fields() {
        let mut rows = [[0.0; 3]; 3];
        let mut field = ForceField::new(3);
        for row in &mut rows {
            field.positions_mut().push(row);
        }
        let mut term = StretchBendContrib::new(&field);
        let mut stbn = MmffStbn {
            kba_ijk: 2.0,
            kba_kji: 5.0,
        };
        let mut angle = MmffAngle {
            theta0: 60.0,
            ka: f64::NAN,
        };
        let mut b1 = MmffBond {
            r0: 1.0,
            kb: f64::NAN,
        };
        let mut b2 = MmffBond {
            r0: 1.0,
            kb: f64::INFINITY,
        };
        term.add_term(field.positions(), 0, 1, 2, &stbn, &angle, &b1, &b2)
            .unwrap();
        stbn.kba_ijk = 99.0;
        stbn.kba_kji = -99.0;
        angle.theta0 = f64::NAN;
        b1.r0 = 8.0;
        b2.r0 = 9.0;
        assert_eq!(term.rest_len1s, [1.0]);
        assert_eq!(term.rest_len2s, [1.0]);
        assert_eq!(term.theta0s, [60.0]);
        assert_eq!(term.force_constants1, [2.0]);
        assert_eq!(term.force_constants2, [5.0]);
        close(energy(&term, &ORTHOGONAL), 904.3545692256258);
        close_grad(&grad(&term, &ORTHOGONAL, &[0.0; 9]), &FIXED_GRAD);
        term.add_term(field.positions(), 0, 1, 2, &stbn, &angle, &b1, &b2)
            .unwrap();
        assert!(term.theta0s[1].is_nan());
        assert_eq!(term.force_constants1[1], 99.0);
    }
    #[test]
    fn mmff_stretch_bend_eight_arrays_append_order_and_growth() {
        let mut specs = Vec::new();
        for i in 0..65 {
            let mut s = base();
            s.r0 = [i as f64, 100.0 + i as f64];
            s.theta0 = 60.0 + i as f64;
            s.k = [i as f64 + 1.0, -(i as f64) - 2.0];
            if i % 2 == 1 {
                s.idx = [2, 1, 0];
            }
            specs.push(s);
        }
        let term = contribution(&specs);
        assert_eq!(
            [
                term.at1_idxs.len(),
                term.at2_idxs.len(),
                term.at3_idxs.len(),
                term.rest_len1s.len(),
                term.rest_len2s.len(),
                term.theta0s.len(),
                term.force_constants1.len(),
                term.force_constants2.len()
            ],
            [65; 8]
        );
        for (i, s) in specs.iter().enumerate() {
            assert_eq!(
                [term.at1_idxs[i], term.at2_idxs[i], term.at3_idxs[i]],
                s.idx.map(|x| x as i16)
            );
            assert_eq!([term.rest_len1s[i], term.rest_len2s[i]], s.r0);
            assert_eq!(term.theta0s[i], s.theta0);
            assert_eq!([term.force_constants1[i], term.force_constants2[i]], s.k);
        }
    }
    #[test]
    fn mmff_stretch_bend_fixed_nonzero_energy_and_full_gradient() {
        let term = contribution(&[base()]);
        close(energy(&term, &ORTHOGONAL), 904.3545692256258);
        close_grad(&grad(&term, &ORTHOGONAL, &[0.0; 9]), &FIXED_GRAD);
    }
    #[test]
    fn mmff_stretch_bend_four_rest_and_distortion_combinations() {
        let both = contribution(&[spec(2.0, 3.0, 90.0, 2.0, 5.0)]);
        assert_eq!(energy(&both, &ORTHOGONAL), 0.0);
        assert_eq!(grad(&both, &ORTHOGONAL, &[0.0; 9]), [0.0; 9]);
        let lengths = contribution(&[spec(2.0, 3.0, 60.0, 2.0, 5.0)]);
        assert_eq!(energy(&lengths, &ORTHOGONAL), 0.0);
        close_grad(
            &grad(&lengths, &ORTHOGONAL, &[0.0; 9]),
            &[
                150.7257615376043,
                0.0,
                0.0,
                -150.7257615376043,
                -376.81440384401077,
                0.0,
                0.0,
                376.81440384401077,
                0.0,
            ],
        );
        let angle = contribution(&[spec(1.0, 1.0, 90.0, 2.0, 5.0)]);
        assert_eq!(energy(&angle, &ORTHOGONAL), 0.0);
        close_grad(
            &grad(&angle, &ORTHOGONAL, &[0.0; 9]),
            &[0.0, -863.595, 0.0, 575.73, 863.595, 0.0, -575.73, 0.0, 0.0],
        );
        close(
            energy(&contribution(&[base()]), &ORTHOGONAL),
            904.3545692256258,
        );
    }
    #[test]
    fn mmff_stretch_bend_oblique_derivative_supplements_literals() {
        let coords = [0.0, 0.1, -0.2, 1.5, 0.3, 0.4, 1.8, 1.4, 0.8];
        let term = contribution(&[spec(1.0, 1.2, 60.0, 0.7, 1.3)]);
        let g = grad(&term, &coords, &[0.0; 9]);
        for i in 0..9 {
            let mut plus = coords;
            let mut minus = coords;
            plus[i] += 1e-6;
            minus[i] -= 1e-6;
            let derivative = (energy(&term, &plus) - energy(&term, &minus)) / 2e-6;
            assert!(
                (g[i] - derivative).abs() < 1e-6,
                "axis{i}:{} != {derivative}",
                g[i]
            );
        }
    }
    #[test]
    fn mmff_stretch_bend_reversal_requires_parameter_orientation() {
        let forward = contribution(&[base()]);
        let mut s = base();
        s.idx = [2, 1, 0];
        let unswapped = contribution(&[s]);
        close(energy(&unswapped, &ORTHOGONAL), 678.2659269192194);
        s.r0.reverse();
        s.k.reverse();
        let swapped = contribution(&[s]);
        close(energy(&swapped, &ORTHOGONAL), 904.3545692256258);
        close_grad(&grad(&swapped, &ORTHOGONAL, &[0.0; 9]), &FIXED_GRAD);
        assert_ne!(
            energy(&forward, &ORTHOGONAL),
            energy(&unswapped, &ORTHOGONAL)
        );
    }
    #[test]
    fn mmff_stretch_bend_translation_and_additive_gradient() {
        let translated = [10.0, -3.0, 7.0, 8.0, -3.0, 7.0, 8.0, 0.0, 7.0];
        let term = contribution(&[base()]);
        let initial = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        let expected: Vec<_> = initial.iter().zip(FIXED_GRAD).map(|(a, b)| a + b).collect();
        close_grad(&grad(&term, &ORTHOGONAL, &initial), &expected);
        assert_eq!(
            grad(&term, &ORTHOGONAL, &initial),
            grad(&term, &translated, &initial)
        );
        assert_eq!(energy(&term, &ORTHOGONAL), energy(&term, &translated));
    }
    #[test]
    fn mmff_stretch_bend_repeated_terms_and_energy_pair_grouping() {
        let twice = contribution(&[base(), base()]);
        close(energy(&twice, &ORTHOGONAL), 1808.7091384512516);
        let expected: Vec<_> = FIXED_GRAD.iter().map(|x| 2.0 * x).collect();
        close_grad(&grad(&twice, &ORTHOGONAL, &[0.0; 9]), &expected);
        let coords = [2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0];
        let grouped = contribution(&[
            spec(1.0, 1.0, 60.0, 1.0, 0.0),
            spec(1.0, 1.0, 60.0, 1e16, -1e16),
        ]);
        // Pair energies cancel before being added to the prior total. Flat
        // addition across all arms loses this small surviving first term.
        close(energy(&grouped, &coords), 75.36288076880216);
    }
    #[test]
    fn mmff_stretch_bend_zero_nan_and_infinite_geometry() {
        for k in [[2.0, 5.0], [0.0, 0.0]] {
            let term = contribution(&[spec(1.0, 1.0, 60.0, k[0], k[1])]);
            for coords in [
                [0.0; 9],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0],
                [2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [f64::NAN, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0],
                [f64::INFINITY, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0],
            ] {
                assert!(energy(&term, &coords).is_nan());
                assert!(grad(&term, &coords, &[4.0; 9]).iter().all(|v| v.is_nan()));
            }
        }
    }
    #[test]
    fn mmff_stretch_bend_ieee_parameters_source_association() {
        let negative = contribution(&[spec(1.0, 1.0, 60.0, -2.0, -5.0)]);
        close(energy(&negative, &ORTHOGONAL), -904.3545692256258);
        let negative_g: Vec<_> = FIXED_GRAD.iter().map(|x| -x).collect();
        close_grad(&grad(&negative, &ORTHOGONAL, &[0.0; 9]), &negative_g);
        for s in [
            spec(1.0, 1.0, f64::NAN, 2.0, 5.0),
            spec(1.0, 1.0, 60.0, f64::NAN, 5.0),
            spec(f64::NAN, 1.0, 60.0, 2.0, 5.0),
        ] {
            let term = contribution(&[s]);
            assert!(energy(&term, &ORTHOGONAL).is_nan());
            assert!(
                grad(&term, &ORTHOGONAL, &[0.0; 9])
                    .iter()
                    .all(|v| v.is_nan())
            );
        }
        let theta = contribution(&[spec(1.0, 1.0, f64::INFINITY, 2.0, 5.0)]);
        assert_eq!(energy(&theta, &ORTHOGONAL), f64::NEG_INFINITY);
        let g = grad(&theta, &ORTHOGONAL, &[0.0; 9]);
        for i in [0, 7] {
            assert_eq!(g[i], f64::NEG_INFINITY);
        }
        for i in [3, 4] {
            assert_eq!(g[i], f64::INFINITY);
        }
        for i in [1, 2, 5, 6, 8] {
            assert!(g[i].is_nan());
        }
        let rest = contribution(&[spec(f64::INFINITY, 1.0, 60.0, 2.0, 5.0)]);
        assert_eq!(energy(&rest, &ORTHOGONAL), f64::NEG_INFINITY);
        let g = grad(&rest, &ORTHOGONAL, &[0.0; 9]);
        for i in [1, 6] {
            assert_eq!(g[i], f64::INFINITY);
        }
        for i in [3, 4] {
            assert_eq!(g[i], f64::NEG_INFINITY);
        }
        for i in [0, 2, 5, 7, 8] {
            assert!(g[i].is_nan());
        }
        let force = contribution(&[spec(1.0, 1.0, 60.0, f64::INFINITY, 5.0)]);
        assert_eq!(energy(&force, &ORTHOGONAL), f64::INFINITY);
        let g = grad(&force, &ORTHOGONAL, &[0.0; 9]);
        assert_eq!(g[6], f64::NEG_INFINITY);
        for i in [0, 1, 2, 3, 4, 5, 7, 8] {
            assert!(g[i].is_nan());
        }
    }
    #[test]
    fn mmff_stretch_bend_signed_zero_storage_and_addition() {
        let term = contribution(&[spec(1.0, 1.0, 60.0, -0.0, 0.0)]);
        assert_eq!(term.force_constants1[0].to_bits(), 0x8000000000000000);
        assert_eq!(term.force_constants2[0].to_bits(), 0);
        assert_eq!(energy(&term, &ORTHOGONAL).to_bits(), 0);
        let g = grad(&term, &ORTHOGONAL, &[-0.0; 9]);
        for v in &g[..3] {
            assert_eq!(v.to_bits(), 0x8000000000000000);
        }
        for v in &g[3..] {
            assert_eq!(v.to_bits(), 0);
        }
        let copied = contribution(&[spec(-0.0, 0.0, -0.0, 0.0, -0.0)]);
        assert_eq!(copied.rest_len1s[0].to_bits(), 0x8000000000000000);
        assert_eq!(copied.rest_len2s[0].to_bits(), 0);
        assert_eq!(copied.theta0s[0].to_bits(), 0x8000000000000000);
    }
    #[test]
    fn mmff_stretch_bend_clip_and_source_sine_floor() {
        let term = contribution(&[spec(0.0, 0.0, 0.0, 1.0, 2.0)]);
        for sign in [1.0, -1.0] {
            let coords = [
                0.0001,
                0.0001,
                0.0001,
                0.0,
                0.0,
                0.0,
                sign * 0.0001,
                sign * 0.0001,
                sign * 0.0001,
            ];
            assert!(energy(&term, &coords).is_finite());
            assert!(
                grad(&term, &coords, &[0.0; 9])
                    .iter()
                    .all(|v| v.is_finite())
            );
        }
        assert_eq!(source_sin_theta(1.0), 1e-8);
        assert_eq!(source_sin_theta(-1.0), 1e-8);
        assert!(source_sin_theta(f64::from_bits(1.0_f64.to_bits() - 1)) > 1e-8);
        // Near parallel cos rounds to1, preserving source floor. Fixed scalar
        // source witness: radial0 and distanceTerm=RAD2DEG*3 gives y=-4.317975.
        let near = contribution(&[spec(0.0, 0.0, 0.0, 1.0, 2.0)]);
        let coords = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1e-10, 0.0];
        close_grad(
            &grad(&near, &coords, &[0.0; 9]),
            &[0.0, -4.317975, 0.0, 0.0, 0.0, 0.0, 0.0, 4.317975, 0.0],
        );
    }
    #[test]
    fn mmff_stretch_bend_sqrt_nan_and_ordered_stdmax() {
        let payload = f64::from_bits(0x7ff8000000000042);
        assert!(source_sin_theta(payload).is_nan());
        assert_eq!(source_sin_theta(0.0), 1.0);
        assert_eq!(source_sin_theta(-0.0), 1.0);
        // AngleBend would map this square to zero before max; StretchBend
        // must retain sqrt NaN as the FIRST operand of source std::max.
        assert!(source_sin_theta(f64::from_bits(1.0_f64.to_bits() + 1)).is_nan());
        let mut value = payload;
        clip_to_one(&mut value);
        assert_eq!(value.to_bits(), 0x7ff8000000000042);
        for (mut value, expected) in [
            (-0.0, 0x8000000000000000),
            (0.0, 0),
            (f64::INFINITY, 1.0_f64.to_bits()),
            (f64::NEG_INFINITY, (-1.0_f64).to_bits()),
        ] {
            clip_to_one(&mut value);
            assert_eq!(value.to_bits(), expected);
        }
    }
    #[test]
    fn mmff_stretch_bend_signed_endpoint_storage_and_distance_errors() {
        for (raw, stored, promoted) in [
            (0, 0, 0),
            (1, 1, 1),
            (32767, 32767, 32767),
            (32768, -32768, 4294934528),
            (65535, -1, 4294967295),
            (65536, 0, 0),
            (65537, 1, 1),
            (65538, 2, 2),
            (u32::MAX, -1, u32::MAX),
        ] {
            assert_eq!(raw as i16, stored);
            assert_eq!(stored as u32, promoted);
        }
        for (idx, argument, index) in [
            ([32767, 1, 2], ForceFieldIndexArgument::I, 32767),
            ([32768, 1, 2], ForceFieldIndexArgument::I, 4294934528),
            ([0, 32768, 2], ForceFieldIndexArgument::J, 4294934528),
            ([0, 1, 32768], ForceFieldIndexArgument::J, 4294934528),
            ([65535, 1, 2], ForceFieldIndexArgument::I, u32::MAX),
        ] {
            let mut s = base();
            s.idx = idx;
            let term = contribution(&[s]);
            let expected = ForceFieldKernelError::IndexOutOfRange {
                argument,
                index,
                upper_bound: 3,
            };
            let mut cache = [-1.0; 6];
            let mut context = EvaluationContext::for_test(&ORTHOGONAL[..6], &mut cache, 3);
            assert_eq!(term.get_energy(&mut context), Err(expected));
            let mut g = [7.0; 6];
            assert_eq!(term.get_grad(&mut context, &mut g), Err(expected));
            assert_eq!(g, [7.0; 6]);
        }
    }
    #[test]
    fn mmff_stretch_bend_wrapped_and_collapsed_raw_distinct_endpoints() {
        let mut s = base();
        s.idx = [65536, 65537, 65538];
        let term = contribution(&[s]);
        assert_eq!(
            [term.at1_idxs[0], term.at2_idxs[0], term.at3_idxs[0]],
            [0, 1, 2]
        );
        close(energy(&term, &ORTHOGONAL), 904.3545692256258);
        close_grad(&grad(&term, &ORTHOGONAL, &[0.0; 9]), &FIXED_GRAD);
        s.idx = [0, 65536, 1];
        let collapsed = contribution(&[s]);
        assert_eq!(
            [
                collapsed.at1_idxs[0],
                collapsed.at2_idxs[0],
                collapsed.at3_idxs[0]
            ],
            [0, 0, 1]
        );
        assert!(energy(&collapsed, &ORTHOGONAL).is_nan());
        let g = grad(&collapsed, &ORTHOGONAL, &[42.0; 9]);
        assert!(g[..6].iter().all(|v| v.is_nan()));
        assert_eq!(g[6..], [42.0; 3]);
    }
    #[test]
    fn mmff_stretch_bend_shared_signed_term_count_boundaries() {
        // Exact native64 source-width observations. These are scalar range
        // tests, not claims of physically allocating billion-term vectors.
        for (length, end, iterations, last) in [
            (0_usize, 0_i32, 0, None),
            (1, 1, 1, Some(0)),
            (2147483647, 2147483647, 2147483647, Some(2147483646)),
            (2147483648, -2147483648, 0, None),
            (4294967295, -1, 0, None),
            (4294967296, 0, 0, None),
            (4294967297, 1, 1, Some(0)),
        ] {
            let range = source_term_indices(length);
            assert_eq!(range.start, 0);
            assert_eq!(range.end, end);
            assert_eq!(range.len(), iterations);
            assert_eq!(range.clone().next(), last.map(|_| 0));
            assert_eq!(range.last(), last);
        }
    }
    #[test]
    fn mmff_stretch_bend_clone_and_trait_copy_eight_independent_arrays() {
        let term = contribution(&[base()]);
        let copied = term.copy();
        let mut clone = term.clone();
        assert_ne!(term.at1_idxs.as_ptr(), clone.at1_idxs.as_ptr());
        assert_ne!(term.at2_idxs.as_ptr(), clone.at2_idxs.as_ptr());
        assert_ne!(term.at3_idxs.as_ptr(), clone.at3_idxs.as_ptr());
        assert_ne!(term.rest_len1s.as_ptr(), clone.rest_len1s.as_ptr());
        assert_ne!(term.rest_len2s.as_ptr(), clone.rest_len2s.as_ptr());
        assert_ne!(term.theta0s.as_ptr(), clone.theta0s.as_ptr());
        assert_ne!(
            term.force_constants1.as_ptr(),
            clone.force_constants1.as_ptr()
        );
        assert_ne!(
            term.force_constants2.as_ptr(),
            clone.force_constants2.as_ptr()
        );
        clone.at1_idxs[0] = 2;
        clone.at3_idxs[0] = 0;
        clone.at2_idxs[0] = 1;
        clone.rest_len1s[0] = 3.0;
        clone.rest_len2s[0] = 2.0;
        clone.theta0s[0] = 90.0;
        clone.force_constants1[0] = 9.0;
        clone.force_constants2[0] = 8.0;
        assert_eq!(term.rest_len1s, [1.0]);
        assert_eq!(term.rest_len2s, [1.0]);
        assert_eq!(term.theta0s, [60.0]);
        assert_eq!(term.force_constants1, [2.0]);
        assert_eq!(term.force_constants2, [5.0]);
        close(energy(&*copied, &ORTHOGONAL), 904.3545692256258);
        close_grad(&grad(&*copied, &ORTHOGONAL, &[0.0; 9]), &FIXED_GRAD);
        assert_eq!(energy(&clone, &ORTHOGONAL), 0.0);
    }
    #[test]
    fn mmff_stretch_bend_distance_cache_and_later_error_partial_writes() {
        let good = contribution(&[base()]);
        let mut bad_spec = base();
        bad_spec.idx = [0, 1, 32768];
        let bad = contribution(&[bad_spec]);
        let expected = ForceFieldKernelError::IndexOutOfRange {
            argument: ForceFieldIndexArgument::J,
            index: 4294934528,
            upper_bound: 3,
        };
        for gradient in [false, true] {
            let mut cache = [-1.0; 6];
            let mut g = [9.0; 6];
            {
                let mut ctx = EvaluationContext::for_test(&ORTHOGONAL[..6], &mut cache, 3);
                if gradient {
                    assert_eq!(bad.get_grad(&mut ctx, &mut g), Err(expected));
                } else {
                    assert_eq!(bad.get_energy(&mut ctx), Err(expected));
                }
            }
            assert_eq!(cache, [-1.0, 2.0, -1.0, -1.0, -1.0, -1.0]);
            assert_eq!(g, [9.0; 6]);
        }
        let mut cache = [-1.0; 6];
        {
            let mut ctx = EvaluationContext::for_test(&ORTHOGONAL, &mut cache, 3);
            assert_eq!(ctx.distance(0, 1), Ok(2.0));
            close(good.get_energy(&mut ctx).unwrap(), 904.3545692256258);
            let mut g = [0.0; 9];
            good.get_grad(&mut ctx, &mut g).unwrap();
            close_grad(&g, &FIXED_GRAD);
        }
        assert_eq!(cache, [-1.0, 2.0, -1.0, -1.0, 3.0, -1.0]);
        let partial = contribution(&[base(), bad_spec]);
        let mut cache = [-1.0; 6];
        let mut ctx = EvaluationContext::for_test(&ORTHOGONAL, &mut cache, 3);
        let mut g = [0.0; 9];
        assert_eq!(partial.get_grad(&mut ctx, &mut g), Err(expected));
        close_grad(&g, &FIXED_GRAD);
    }
    #[test]
    fn mmff_stretch_bend_real_field_initialization_copy_rebind_and_cache() {
        let mut p1 = [2.0, 0.0, 0.0];
        let mut p2 = [0.0; 3];
        let mut p3 = [0.0, 3.0, 0.0];
        let mut field = ForceField::new(3);
        field.positions_mut().push(&mut p1);
        field.positions_mut().push(&mut p2);
        field.positions_mut().push(&mut p3);
        let mut term = StretchBendContrib::new(&field);
        append(&mut term, &field, base()).unwrap();
        field.add_contribution(Box::new(term));
        assert_eq!(
            cf3d_bld_b05_calc_energy(&mut field, &ORTHOGONAL),
            Err(ForceFieldKernelError::NotInitialized)
        );
        field.initialize().unwrap();
        close(field.calc_energy_current(None).unwrap(), 904.3545692256258);
        // getGrad reuses cached2,3 even though supplied UNIT coordinates have
        // actual arm lengths1,1; energy(pos) alone invalidates that cache.
        let mut stale = [0.0; 9];
        cf3d_bld_b05_calc_grad(&mut field, &UNIT, &mut stale).unwrap();
        close_grad(
            &stale,
            &[
                75.36288076880216,
                -287.865,
                0.0,
                212.50211923119787,
                162.26019871866308,
                0.0,
                -287.865,
                125.60480128133693,
                0.0,
            ],
        );
        close(
            cf3d_bld_b05_calc_energy(&mut field, &ORTHOGONAL).unwrap(),
            904.3545692256258,
        );
        let mut copied = field.copy();
        assert!(copied.positions().is_empty());
        assert_eq!(
            cf3d_bld_b05_calc_energy(&mut copied, &ORTHOGONAL),
            Err(ForceFieldKernelError::NotInitialized)
        );
        let mut c1 = [2.0, 0.0, 0.0];
        let mut c2 = [0.0; 3];
        let mut c3 = [0.0, 3.0, 0.0];
        copied.positions_mut().push(&mut c1);
        copied.positions_mut().push(&mut c2);
        copied.positions_mut().push(&mut c3);
        copied.initialize().unwrap();
        close(
            cf3d_bld_b05_calc_energy(&mut copied, &ORTHOGONAL).unwrap(),
            904.3545692256258,
        );
        let mut g = [0.0; 9];
        cf3d_bld_b05_calc_grad(&mut copied, &ORTHOGONAL, &mut g).unwrap();
        close_grad(&g, &FIXED_GRAD);
        assert_eq!(cf3d_bld_b05_calc_energy(&mut copied, &UNIT), Ok(0.0));
        close(field.calc_energy_current(None).unwrap(), 904.3545692256258);
    }
    #[test]
    fn mmff_stretch_bend_mixed_accepted_private_contributions() {
        use super::super::{angle_bend::AngleBendContrib, bond_stretch::BondStretchContrib};
        let mut p1 = [2.0, 0.0, 0.0];
        let mut p2 = [0.0; 3];
        let mut p3 = [0.0, 3.0, 0.0];
        let mut field = ForceField::new(3);
        field.positions_mut().push(&mut p1);
        field.positions_mut().push(&mut p2);
        field.positions_mut().push(&mut p3);
        let mut stretch = StretchBendContrib::new(&field);
        append(&mut stretch, &field, base()).unwrap();
        let mut bond = BondStretchContrib::new(&field);
        bond.add_term(field.positions(), 0, 1, &MmffBond { r0: 1.0, kb: 2.0 })
            .unwrap();
        let mut angle = AngleBendContrib::new(&field);
        angle
            .add_term(
                field.positions(),
                0,
                1,
                2,
                &MmffAngle {
                    theta0: 60.0,
                    ka: 1.0,
                },
                &super::super::params::MmffProp {
                    atno: 6,
                    crd: 2,
                    val: 4,
                    pilp: 0,
                    mltb: 0,
                    arom: 0,
                    linh: 0,
                    sbmb: 0,
                },
            )
            .unwrap();
        field.add_contribution(Box::new(stretch));
        field.add_contribution(Box::new(bond));
        field.add_contribution(Box::new(angle));
        field.initialize().unwrap();
        close(
            cf3d_bld_b05_calc_energy(&mut field, &ORTHOGONAL).unwrap(),
            1111.8622929466527,
        );
        // Existing accepted source witnesses: bond191.91,angle15.597723721027004.
        let mut g = [0.0; 9];
        cf3d_bld_b05_calc_grad(&mut field, &ORTHOGONAL, &mut g).unwrap();
        close_grad(
            &g,
            &[
                918.3657615376043,
                -889.4384667690963,
                0.0,
                -325.4067836915401,
                512.6240629250856,
                0.0,
                -592.9589778460642,
                376.81440384401077,
                0.0,
            ],
        );
    }
}
