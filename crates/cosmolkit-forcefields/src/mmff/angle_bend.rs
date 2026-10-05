//! Packed private MMFF angle-bend contribution.
//!
//! Source: RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8,
//! Code/ForceField/MMFF/AngleBend.{cpp,h} (BSD license).

use super::{
    bond_stretch::source_term_indices,
    bonded::calc_angle_bend_energy,
    numerical::{calc_angle_bend_grad, calc_cos_theta},
    params::{MmffAngle, MmffProp},
};
use crate::{
    geometry::Point3,
    kernel::{
        AngleIndexArgument, EvaluationContext, ForceField, ForceFieldContribution,
        ForceFieldKernelError,
    },
    uff::params::{DEG2RAD, RAD2DEG, clip_to_one},
};

// RDKit✔️✔️:   std::vector<bool> d_isLinear;
// RDKit✔️✔️:   std::vector<int16_t> d_at1Idxs, d_at2Idxs, d_at3Idxs;
// RDKit✔️✔️:   std::vector<double> d_ka, d_theta0;
// Packed flags retain the source's 64 flags per native word and word growth.
#[derive(Clone, Debug)]
pub(super) struct AngleBendContrib {
    at1_idxs: Vec<i16>,
    at2_idxs: Vec<i16>,
    at3_idxs: Vec<i16>,
    linear_bits: Vec<u64>,
    theta0: Vec<f64>,
    ka: Vec<f64>,
}

pub(super) fn source_sin_theta(cos_theta: f64) -> f64 {
    // RDKit✔️✔️:     double sinThetaSq = 1.0 - cosTheta * cosTheta;
    // RDKit✔️✔️:     double sinTheta =
    // RDKit✔️✔️:         std::max(((sinThetaSq > 0.0) ? sqrt(sinThetaSq) : 0.0), 1.0e-8);
    // Behavior: NaN square follows the false branch. Preserve std::max's
    // first operand on unordered comparisons; no f64::max NaN rewriting.
    // Complexity: fixed scalar operations and no allocation, O(1).
    let sin_theta_sq = 1.0 - cos_theta * cos_theta;
    let candidate = if sin_theta_sq > 0.0 {
        sin_theta_sq.sqrt()
    } else {
        0.0
    };
    if candidate < 1.0e-8 {
        1.0e-8
    } else {
        candidate
    }
}

pub(super) fn point_at(coordinates: &[f64], index: usize) -> Point3 {
    // RDKit✔️✔️:     RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
    // RDKit✔️✔️:                        pos[3 * d_at1Idx + 2]);
    // Behavior: all three points use this same source addressing only after
    // both checked distances. Storage is the existing valid 3D caller contract.
    // Complexity: three indexed reads and a stack point, O(1), no allocation.
    let offset = 3 * index;
    Point3 {
        x: coordinates[offset],
        y: coordinates[offset + 1],
        z: coordinates[offset + 2],
    }
}

impl AngleBendContrib {
    pub(super) fn new(_owner: &ForceField<'_>) -> Self {
        // RDKit✔️✔️: AngleBendContrib::AngleBendContrib(ForceField *owner) {
        // RDKit✔️✔️:   PRECONDITION(owner, "bad owner");
        // RDKit✔️✔️:   dp_forceField = owner;
        // RDKit✔️✔️: }
        // Behavior: borrowed owner is non-null; evaluation borrows the existing
        // owner-derived context. An ownerless source default is not evaluable.
        // Complexity: six empty arrays; O(1), no allocation or owner clone.
        Self {
            at1_idxs: Vec::new(),
            at2_idxs: Vec::new(),
            at3_idxs: Vec::new(),
            linear_bits: Vec::new(),
            theta0: Vec::new(),
            ka: Vec::new(),
        }
    }

    pub(super) fn add_term(
        &mut self,
        positions: &[&mut [f64]],
        idx1: u32,
        idx2: u32,
        idx3: u32,
        params: &MmffAngle,
        central: &MmffProp,
    ) -> Result<(), ForceFieldKernelError> {
        // RDKit✔️✔️: void AngleBendContrib::addTerm(unsigned int idx1,
        // RDKit✔️✔️:                                unsigned int idx2,
        // RDKit✔️✔️:                                unsigned int idx3,
        // RDKit✔️✔️:                                const ForceFields::MMFF::MMFFAngle *mmffAngleParams,
        // RDKit✔️✔️:                                const ForceFields::MMFF::MMFFProp *mmffPropParamsCentralAtom) {
        // RDKit✔️✔️:   URANGE_CHECK(idx1, dp_forceField->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx2, dp_forceField->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx3, dp_forceField->positions().size());
        // RDKit✔️✔️:   PRECONDITION(((idx1 != idx2) && (idx2 != idx3) && (idx1 != idx3)),
        // RDKit✔️✔️:                "degenerate points");
        // RDKit✔️✔️:   d_at1Idxs.push_back(idx1);
        // RDKit✔️✔️:   d_at2Idxs.push_back(idx2);
        // RDKit✔️✔️:   d_at3Idxs.push_back(idx3);
        // RDKit✔️✔️:   d_isLinear.push_back(mmffPropParamsCentralAtom->linh > 0u);
        // RDKit✔️✔️:   d_theta0.push_back(mmffAngleParams->theta0);
        // RDKit✔️✔️:   d_ka.push_back(mmffAngleParams->ka);
        // RDKit✔️✔️: }
        // Behavior: raw unsigned checks precede the distinct precondition and
        // every append. Source int16_t narrowing is modulo 2^16 with signed
        // interpretation; later evaluation must not recheck distinctness.
        // Values are copied without extra numeric restrictions.
        // Complexity: fixed checks and six amortized O(1) separate-array
        // appends; flags use native packed words rather than byte booleans.
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
        if idx1 == idx2 || idx2 == idx3 || idx1 == idx3 {
            return Err(ForceFieldKernelError::AngleDegeneratePoints);
        }
        let term_index = self.at1_idxs.len();
        self.at1_idxs.push(idx1 as i16);
        self.at2_idxs.push(idx2 as i16);
        self.at3_idxs.push(idx3 as i16);
        self.push_linear(term_index, central.linh > 0);
        self.theta0.push(params.theta0);
        self.ka.push(params.ka);
        Ok(())
    }

    fn push_linear(&mut self, term_index: usize, linear: bool) {
        // RDKit✔️✔️:   d_isLinear.push_back(mmffPropParamsCentralAtom->linh > 0u);
        // Behavior: append the source boolean at its packed bit position.
        // Complexity: 64 flags per u64; observed native capacities grow
        // 0/64/128/256 bits. Word doubling gives amortized O(1) append and
        // equivalent packed payload/allocation growth, no expanded bool array.
        if term_index % 64 == 0 {
            if self.linear_bits.len() == self.linear_bits.capacity() {
                self.linear_bits
                    .reserve_exact(self.linear_bits.len().max(1));
            }
            self.linear_bits.push(0);
        }
        if linear {
            self.linear_bits[term_index / 64] |= 1_u64 << (term_index % 64);
        }
    }

    fn is_linear(&self, term_index: usize) -> bool {
        // RDKit✔️✔️:     double dE_dTheta = (d_isLinear[i] ? -MDYNE_A_TO_KCAL_MOL * d_ka[i] * sinTheta
        // Behavior: read the indexed packed source boolean; also used by energy.
        // Complexity: O(1) word lookup and bit mask, no scan or allocation.
        self.linear_bits[term_index / 64] & (1_u64 << (term_index % 64)) != 0
    }
}

impl ForceFieldContribution for AngleBendContrib {
    fn get_energy(
        &self,
        context: &mut EvaluationContext<'_>,
    ) -> Result<f64, ForceFieldKernelError> {
        // RDKit✔️✔️: double AngleBendContrib::getEnergy(double *pos) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:   double res = 0.0;
        // RDKit✔️✔️:   const int numTerms = d_at1Idxs.size();
        // RDKit✔️✔️:
        // RDKit✔️✔️:   for (int i = 0; i < numTerms; i++) {
        // RDKit✔️✔️:     const int d_at1Idx = d_at1Idxs[i];
        // RDKit✔️✔️:     const int d_at2Idx = d_at2Idxs[i];
        // RDKit✔️✔️:     const int d_at3Idx = d_at3Idxs[i];
        // RDKit✔️✔️:
        // RDKit✔️✔️:     double dist1 = dp_forceField->distance(d_at1Idx, d_at2Idx, pos);
        // RDKit✔️✔️:     double dist2 = dp_forceField->distance(d_at2Idx, d_at3Idx, pos);
        // RDKit✔️✔️:
        // RDKit✔️✔️:     RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
        // RDKit✔️✔️:                        pos[3 * d_at1Idx + 2]);
        // RDKit✔️✔️:     RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
        // RDKit✔️✔️:                        pos[3 * d_at2Idx + 2]);
        // RDKit✔️✔️:     RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
        // RDKit✔️✔️:                        pos[3 * d_at3Idx + 2]);
        // RDKit✔️✔️:
        // RDKit✔️✔️:     res += Utils::calcAngleBendEnergy(
        // RDKit✔️✔️:         d_theta0[i], d_ka[i], d_isLinear[i],
        // RDKit✔️✔️:         Utils::calcCosTheta(p1, p2, p3, dist1, dist2));
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        // Behavior: sign-promote stored int16_t endpoints, then convert to the
        // distance API's unsigned indices. Both calls precede coordinate reads;
        // source errors and cache effects therefore occur in the original order.
        // Reuse the sole signed32 range and source scalar/geometry owners.
        // Complexity: O(max(0, signed32(terms))) time, O(1) stack temporaries;
        // two indexed cached distances per term, no allocation or state clones.
        let mut result = 0.0;
        for i in source_term_indices(self.at1_idxs.len()) {
            let i = i as usize;
            let at1 = self.at1_idxs[i];
            let at2 = self.at2_idxs[i];
            let at3 = self.at3_idxs[i];
            let dist1 = context.distance(at1 as u32, at2 as u32)?;
            let dist2 = context.distance(at2 as u32, at3 as u32)?;
            let coordinates = context.coordinates();
            let p1 = point_at(coordinates, at1 as usize);
            let p2 = point_at(coordinates, at2 as usize);
            let p3 = point_at(coordinates, at3 as usize);
            result += calc_angle_bend_energy(
                self.theta0[i],
                self.ka[i],
                self.is_linear(i),
                calc_cos_theta(p1, p2, p3, dist1, dist2),
            );
        }
        Ok(result)
    }

    fn get_grad(
        &self,
        context: &mut EvaluationContext<'_>,
        gradient: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        // RDKit✔️✔️: void AngleBendContrib::getGrad(double *pos, double *grad) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:   PRECONDITION(grad, "bad vector");
        // RDKit✔️✔️:
        // RDKit✔️✔️:   const int numTerms = d_at1Idxs.size();
        // RDKit✔️✔️:   for (int i =0; i < numTerms; i++) {
        // RDKit✔️✔️:     const int d_at1Idx = d_at1Idxs[i];
        // RDKit✔️✔️:     const int d_at2Idx = d_at2Idxs[i];
        // RDKit✔️✔️:     const int d_at3Idx = d_at3Idxs[i];
        // RDKit✔️✔️:
        // RDKit✔️✔️:     double dist[2] = {dp_forceField->distance(d_at1Idx, d_at2Idx, pos),
        // RDKit✔️✔️:                       dp_forceField->distance(d_at2Idx, d_at3Idx, pos)};
        // RDKit✔️✔️:
        // RDKit✔️✔️:     RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
        // RDKit✔️✔️:                        pos[3 * d_at1Idx + 2]);
        // RDKit✔️✔️:     RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
        // RDKit✔️✔️:                        pos[3 * d_at2Idx + 2]);
        // RDKit✔️✔️:     RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
        // RDKit✔️✔️:                        pos[3 * d_at3Idx + 2]);
        // RDKit✔️✔️:     double *g[3] = {&(grad[3 * d_at1Idx]), &(grad[3 * d_at2Idx]),
        // RDKit✔️✔️:                     &(grad[3 * d_at3Idx])};
        // RDKit✔️✔️:     RDGeom::Point3D r[2] = {(p1 - p2) / dist[0], (p3 - p2) / dist[1]};
        // RDKit✔️✔️:     double cosTheta = r[0].dotProduct(r[1]);
        // RDKit✔️✔️:     clipToOne(cosTheta);
        // RDKit✔️✔️:     double sinThetaSq = 1.0 - cosTheta * cosTheta;
        // RDKit✔️✔️:     double sinTheta =
        // RDKit✔️✔️:         std::max(((sinThetaSq > 0.0) ? sqrt(sinThetaSq) : 0.0), 1.0e-8);
        // RDKit✔️✔️:
        // RDKit✔️✔️:     // use the chain rule:
        // RDKit✔️✔️:     // dE/dx = dE/dTheta * dTheta/dx
        // RDKit✔️✔️:
        // RDKit✔️✔️:     // dE/dTheta is independent of cartesians:
        // RDKit✔️✔️:     double angleTerm = RAD2DEG * acos(cosTheta) - d_theta0[i];
        // RDKit✔️✔️:     double const cb = -0.006981317;
        // RDKit✔️✔️:     double const c2 = MDYNE_A_TO_KCAL_MOL * DEG2RAD * DEG2RAD;
        // RDKit✔️✔️:
        // RDKit✔️✔️:     double dE_dTheta = (d_isLinear[i] ? -MDYNE_A_TO_KCAL_MOL * d_ka[i] * sinTheta
        // RDKit✔️✔️:                                    : RAD2DEG * c2 * d_ka[i] * angleTerm *
        // RDKit✔️✔️:                                          (1.0 + 1.5 * cb * angleTerm));
        // RDKit✔️✔️:
        // RDKit✔️✔️:     Utils::calcAngleBendGrad(r, dist, g, dE_dTheta, cosTheta, sinTheta);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // RDKit✔️✔️: constexpr double MDYNE_A_TO_KCAL_MOL = 143.9325;
        // Behavior: preserve both distances before point/gradient addressing,
        // normalized component divisions, clipping, source sine floor and
        // exact derivative association. Borrow rows without copying; the sole
        // numerical helper retains sequential additive writes even if narrowing
        // aliases endpoints. Valid gradient storage is the kernel caller contract.
        // Complexity: O(max(0, signed32(terms))) time, O(1) stack temporaries,
        // no allocation, cloned gradient rows or duplicate geometry algorithms.
        for i in source_term_indices(self.at1_idxs.len()) {
            let i = i as usize;
            let at1 = self.at1_idxs[i];
            let at2 = self.at2_idxs[i];
            let at3 = self.at3_idxs[i];
            let dist = [
                context.distance(at1 as u32, at2 as u32)?,
                context.distance(at2 as u32, at3 as u32)?,
            ];
            let coordinates = context.coordinates();
            let p1 = point_at(coordinates, at1 as usize);
            let p2 = point_at(coordinates, at2 as usize);
            let p3 = point_at(coordinates, at3 as usize);
            let r = [
                Point3::difference(&p1, &p2).divided(dist[0]),
                Point3::difference(&p3, &p2).divided(dist[1]),
            ];
            let mut cos_theta = r[0].dot_product(&r[1]);
            clip_to_one(&mut cos_theta);
            let sin_theta = source_sin_theta(cos_theta);
            let angle_term = RAD2DEG * cos_theta.acos() - self.theta0[i];
            let cb = -0.006981317;
            let c2 = 143.9325 * DEG2RAD * DEG2RAD;
            let de_dtheta = if self.is_linear(i) {
                -143.9325 * self.ka[i] * sin_theta
            } else {
                RAD2DEG * c2 * self.ka[i] * angle_term * (1.0 + 1.5 * cb * angle_term)
            };
            calc_angle_bend_grad(
                &r,
                &dist,
                gradient.as_chunks_mut::<3>().0,
                [at1 as usize, at2 as usize, at3 as usize],
                de_dtheta,
                cos_theta,
                sin_theta,
            )
            .expect("valid 3D gradient storage is the kernel caller contract");
        }
        Ok(())
    }

    fn copy(&self) -> Box<dyn ForceFieldContribution> {
        // RDKit✔️✔️:   AngleBendContrib *copy() const override {
        // RDKit✔️✔️:     return new AngleBendContrib(*this);
        // RDKit✔️✔️:   }
        // Behavior: six independently owned arrays; existing ForceField copy
        // supplies a new owner-derived context after rebind/initialize.
        // Complexity: O(terms) copying with O(ceil(terms/64)) packed flags;
        // one allocation per populated array, no owner/coordinate clone.
        Box::new(self.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernel::{
        ForceFieldIndexArgument, cf3d_bld_b05_calc_energy, cf3d_bld_b05_calc_grad,
    };

    const RIGHT: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
    const NONLINEAR_GRAD: [f64; 9] = [
        0.0,
        -51.68693353819263,
        0.0,
        51.68693353819263,
        51.68693353819263,
        0.0,
        -51.68693353819263,
        0.0,
        0.0,
    ];
    const LINEAR_GRAD: [f64; 9] = [
        0.0, 143.9325, 0.0, -143.9325, -143.9325, 0.0, 143.9325, 0.0, 0.0,
    ];

    fn central(linh: u8) -> MmffProp {
        MmffProp {
            atno: 6,
            crd: 2,
            val: 4,
            pilp: 0,
            mltb: 0,
            arom: 0,
            linh,
            sbmb: 0,
        }
    }

    fn contribution(specs: &[(u32, u32, u32, f64, f64, u8)]) -> AngleBendContrib {
        let count = specs
            .iter()
            .flat_map(|s| [s.0, s.1, s.2])
            .max()
            .map_or(0, |i| i as usize + 1);
        let mut rows = vec![[0.0; 3]; count];
        let mut field = ForceField::new(3);
        for row in &mut rows {
            field.positions_mut().push(row);
        }
        let mut term = AngleBendContrib::new(&field);
        for &(i, j, k, theta0, ka, linh) in specs {
            term.add_term(
                field.positions(),
                i,
                j,
                k,
                &MmffAngle { theta0, ka },
                &central(linh),
            )
            .unwrap();
        }
        term
    }

    fn energy(term: &dyn ForceFieldContribution, coordinates: &[f64]) -> f64 {
        let points = coordinates.len() / 3;
        let mut cache = vec![-1.0; points * (points + 1) / 2];
        term.get_energy(&mut EvaluationContext::for_test(
            coordinates,
            &mut cache,
            points as u32,
        ))
        .unwrap()
    }

    fn gradient(
        term: &dyn ForceFieldContribution,
        coordinates: &[f64],
        initial: &[f64],
    ) -> Vec<f64> {
        let points = coordinates.len() / 3;
        let mut cache = vec![-1.0; points * (points + 1) / 2];
        let mut result = initial.to_vec();
        term.get_grad(
            &mut EvaluationContext::for_test(coordinates, &mut cache, points as u32),
            &mut result,
        )
        .unwrap();
        result
    }

    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 1.0e-10,
            "{actual:?} != {expected:?}"
        );
    }

    fn close_gradient(actual: &[f64], expected: &[f64]) {
        assert_eq!(actual.len(), expected.len());
        for (&a, &e) in actual.iter().zip(expected) {
            close(a, e);
        }
    }

    #[test]
    fn mmff_angle_bend_empty_owner_and_copy() {
        let term = AngleBendContrib::new(&ForceField::new(3));
        assert_eq!(
            [
                term.at1_idxs.capacity(),
                term.at2_idxs.capacity(),
                term.at3_idxs.capacity(),
                term.linear_bits.capacity(),
                term.theta0.capacity(),
                term.ka.capacity()
            ],
            [0; 6]
        );
        assert_eq!(energy(&term, &[]).to_bits(), 0);
        assert_eq!(gradient(&term, &[], &[1.0, -2.0, 3.0]), [1.0, -2.0, 3.0]);
        assert_eq!(energy(&term.clone(), &[]).to_bits(), 0);
        assert_eq!(energy(&*term.copy(), &[]).to_bits(), 0);
    }

    #[test]
    fn mmff_angle_bend_raw_bounds_then_distinct() {
        let mut rows = [[0.0; 3]; 3];
        let mut extra = [0.0; 3];
        let mut field = ForceField::new(3);
        for row in &mut rows {
            field.positions_mut().push(row);
        }
        let mut term = AngleBendContrib::new(&field);
        let params = MmffAngle {
            theta0: 60.0,
            ka: 1.0,
        };
        for (indices, expected) in [
            (
                [u32::MAX, 3, 3],
                ForceFieldKernelError::AngleIndexOutOfRange {
                    argument: AngleIndexArgument::First,
                    index: u32::MAX,
                    upper_bound: 3,
                },
            ),
            (
                [0, u32::MAX, 3],
                ForceFieldKernelError::AngleIndexOutOfRange {
                    argument: AngleIndexArgument::Second,
                    index: u32::MAX,
                    upper_bound: 3,
                },
            ),
            (
                [0, 0, u32::MAX],
                ForceFieldKernelError::AngleIndexOutOfRange {
                    argument: AngleIndexArgument::Third,
                    index: u32::MAX,
                    upper_bound: 3,
                },
            ),
            (
                [3, 3, 3],
                ForceFieldKernelError::AngleIndexOutOfRange {
                    argument: AngleIndexArgument::First,
                    index: 3,
                    upper_bound: 3,
                },
            ),
            ([0, 0, 2], ForceFieldKernelError::AngleDegeneratePoints),
            ([0, 2, 2], ForceFieldKernelError::AngleDegeneratePoints),
            ([0, 1, 0], ForceFieldKernelError::AngleDegeneratePoints),
        ] {
            assert_eq!(
                term.add_term(
                    field.positions(),
                    indices[0],
                    indices[1],
                    indices[2],
                    &params,
                    &central(0)
                ),
                Err(expected)
            );
            assert_eq!(
                [
                    term.at1_idxs.len(),
                    term.at2_idxs.len(),
                    term.at3_idxs.len(),
                    term.linear_bits.len(),
                    term.theta0.len(),
                    term.ka.len()
                ],
                [0; 6]
            );
        }
        term.add_term(field.positions(), 0, 1, 2, &params, &central(0))
            .unwrap();
        let old = format!("{term:?}");
        assert_eq!(
            term.add_term(field.positions(), 0, 1, 3, &params, &central(1)),
            Err(ForceFieldKernelError::AngleIndexOutOfRange {
                argument: AngleIndexArgument::Third,
                index: 3,
                upper_bound: 3
            })
        );
        assert_eq!(format!("{term:?}"), old);
        field.positions_mut().push(&mut extra);
        term.add_term(field.positions(), 0, 1, 3, &params, &central(1))
            .unwrap();
        assert_eq!(term.at3_idxs, [2, 3]);
    }

    #[test]
    fn mmff_angle_bend_parameter_values_are_copied() {
        let mut rows = [[0.0; 3]; 3];
        let mut field = ForceField::new(3);
        for row in &mut rows {
            field.positions_mut().push(row);
        }
        let mut term = AngleBendContrib::new(&field);
        let mut params = MmffAngle {
            theta0: 60.0,
            ka: 1.0,
        };
        let mut prop = central(0);
        for linh in [0, 1, 255] {
            prop.linh = linh;
            term.add_term(field.positions(), 0, 1, 2, &params, &prop)
                .unwrap();
        }
        params.theta0 = f64::NAN;
        params.ka = f64::INFINITY;
        prop.linh = 0;
        assert_eq!(term.theta0, [60.0; 3]);
        assert_eq!(term.ka, [1.0; 3]);
        assert_eq!(term.linear_bits, [6]);
        term.add_term(field.positions(), 0, 1, 2, &params, &prop)
            .unwrap();
        assert!(term.theta0[3].is_nan());
        assert_eq!(term.ka[3], f64::INFINITY);
    }

    #[test]
    fn mmff_angle_bend_nonlinear_fixed_rest_and_right_angle() {
        let rest = contribution(&[(0, 1, 2, 90.0, 1.0, 0)]);
        assert_eq!(energy(&rest, &RIGHT), 0.0);
        assert_eq!(gradient(&rest, &RIGHT, &[0.0; 9]), [0.0; 9]);
        let term = contribution(&[(0, 1, 2, 60.0, 1.0, 0)]);
        close(energy(&term, &RIGHT), 15.597723721027004);
        close_gradient(&gradient(&term, &RIGHT, &[0.0; 9]), &NONLINEAR_GRAD);
    }

    #[test]
    fn mmff_angle_bend_linear_fixed_angles_and_ignored_rest() {
        for theta0 in [60.0, f64::NAN, f64::INFINITY] {
            let term = contribution(&[(0, 1, 2, theta0, 1.0, 255)]);
            close(energy(&term, &RIGHT), 143.9325);
            close_gradient(&gradient(&term, &RIGHT, &[0.0; 9]), &LINEAR_GRAD);
            let parallel = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0];
            close(energy(&term, &parallel), 287.865);
            assert_eq!(gradient(&term, &parallel, &[0.0; 9]), [0.0; 9]);
            let antiparallel = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0, 0.0, 0.0];
            assert_eq!(energy(&term, &antiparallel), 0.0);
            assert_eq!(gradient(&term, &antiparallel, &[0.0; 9]), [0.0; 9]);
        }
    }

    #[test]
    fn mmff_angle_bend_oblique_energy_derivative() {
        let coordinates = [0.0, 0.1, -0.2, 1.5, 0.3, 0.4, 1.8, 1.4, 0.8];
        for linh in [0, 1] {
            let term = contribution(&[(0, 1, 2, 45.0, 0.7, linh)]);
            let grad = gradient(&term, &coordinates, &[0.0; 9]);
            for i in 0..9 {
                let mut plus = coordinates;
                let mut minus = coordinates;
                plus[i] += 1.0e-6;
                minus[i] -= 1.0e-6;
                let derivative = (energy(&term, &plus) - energy(&term, &minus)) / 2.0e-6;
                assert!(
                    (grad[i] - derivative).abs() < 1.0e-6,
                    "branch {linh}, axis {i}: {} != {derivative}",
                    grad[i]
                );
            }
        }
    }

    #[test]
    fn mmff_angle_bend_translation_reversal_and_addition() {
        let translated = [9.0, -3.0, 7.0, 8.0, -3.0, 7.0, 8.0, -2.0, 7.0];
        let initial = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        for (linh, delta) in [(0, NONLINEAR_GRAD), (1, LINEAR_GRAD)] {
            let term = contribution(&[(0, 1, 2, 60.0, 1.0, linh)]);
            let reversed = contribution(&[(2, 1, 0, 60.0, 1.0, linh)]);
            let expected: Vec<_> = initial.iter().zip(delta).map(|(x, d)| x + d).collect();
            close_gradient(&gradient(&term, &RIGHT, &initial), &expected);
            assert_eq!(
                gradient(&term, &RIGHT, &initial),
                gradient(&term, &translated, &initial)
            );
            assert_eq!(
                gradient(&term, &RIGHT, &initial),
                gradient(&reversed, &RIGHT, &initial)
            );
            assert_eq!(energy(&term, &RIGHT), energy(&term, &translated));
        }
    }

    #[test]
    fn mmff_angle_bend_ordered_repeated_and_reversed_terms() {
        let term = contribution(&[
            (0, 1, 2, 60.0, 1.0, 0),
            (2, 1, 0, 60.0, 2.0, 1),
            (0, 1, 2, 60.0, -1.0, 255),
        ]);
        close(energy(&term, &RIGHT), 15.597723721027004 + 143.9325);
        let expected: Vec<_> = NONLINEAR_GRAD
            .iter()
            .zip(LINEAR_GRAD)
            .map(|(a, b)| a + b)
            .collect();
        close_gradient(&gradient(&term, &RIGHT, &[0.0; 9]), &expected);
        assert_eq!(term.at1_idxs, [0, 2, 0]);
        assert_eq!(term.at2_idxs, [1, 1, 1]);
        assert_eq!(term.at3_idxs, [2, 0, 2]);
        assert_eq!(term.ka, [1.0, 2.0, -1.0]);
        let ordered = contribution(&[
            (0, 1, 2, 60.0, 1.0e16, 0),
            (0, 1, 2, 60.0, -1.0e16, 0),
            (0, 1, 2, 60.0, 1.0, 0),
        ]);
        close(energy(&ordered, &RIGHT), 15.597723721027004);
        close_gradient(&gradient(&ordered, &RIGHT, &[0.0; 9]), &NONLINEAR_GRAD);
    }

    #[test]
    fn mmff_angle_bend_zero_and_nan_geometry() {
        for linh in [0, 1] {
            for ka in [0.0, 1.0] {
                let term = contribution(&[(0, 1, 2, 60.0, ka, linh)]);
                for coordinates in [
                    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
                    [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    [0.0; 9],
                    [f64::NAN, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
                ] {
                    assert!(energy(&term, &coordinates).is_nan());
                    assert!(
                        gradient(&term, &coordinates, &[3.0; 9])
                            .iter()
                            .all(|v| v.is_nan())
                    );
                }
            }
        }
    }

    #[test]
    fn mmff_angle_bend_ieee_parameter_arithmetic() {
        for linh in [0, 1] {
            for ka in [0.0, -0.0] {
                let term = contribution(&[(0, 1, 2, 60.0, ka, linh)]);
                assert_eq!(energy(&term, &RIGHT).to_bits(), 0);
                assert_eq!(gradient(&term, &RIGHT, &[0.0; 9]), [0.0; 9]);
            }
            let nan = contribution(&[(0, 1, 2, 60.0, f64::NAN, linh)]);
            assert!(energy(&nan, &RIGHT).is_nan());
            assert!(gradient(&nan, &RIGHT, &[0.0; 9]).iter().all(|v| v.is_nan()));
            let infinite = contribution(&[(0, 1, 2, 60.0, f64::INFINITY, linh)]);
            assert_eq!(energy(&infinite, &RIGHT), f64::INFINITY);
            let g = gradient(&infinite, &RIGHT, &[0.0; 9]);
            let sign = if linh == 0 { 1.0 } else { -1.0 };
            assert_eq!(g[1], -sign * f64::INFINITY);
            assert_eq!(g[3], sign * f64::INFINITY);
            assert_eq!(g[4], sign * f64::INFINITY);
            assert_eq!(g[6], -sign * f64::INFINITY);
            for i in [0, 2, 5, 7, 8] {
                assert!(g[i].is_nan());
            }
        }
        let negative = contribution(&[(0, 1, 2, 60.0, -1.0, 0)]);
        close(energy(&negative, &RIGHT), -15.597723721027004);
        let expected: Vec<_> = NONLINEAR_GRAD.iter().map(|x| -x).collect();
        close_gradient(&gradient(&negative, &RIGHT, &[0.0; 9]), &expected);
        let nan_rest = contribution(&[(0, 1, 2, f64::NAN, 1.0, 0)]);
        assert!(energy(&nan_rest, &RIGHT).is_nan());
        assert!(
            gradient(&nan_rest, &RIGHT, &[0.0; 9])
                .iter()
                .all(|v| v.is_nan())
        );
        let infinite_rest = contribution(&[(0, 1, 2, f64::INFINITY, 1.0, 0)]);
        // Source left association: finite * -inf * -inf * +inf = +inf.
        // The source derivative is -inf, so only its zero-product axes are NaN.
        assert_eq!(energy(&infinite_rest, &RIGHT), f64::INFINITY);
        let g = gradient(&infinite_rest, &RIGHT, &[0.0; 9]);
        for i in [1, 6] {
            assert_eq!(g[i], f64::INFINITY);
        }
        for i in [3, 4] {
            assert_eq!(g[i], f64::NEG_INFINITY);
        }
        for i in [0, 2, 5, 7, 8] {
            assert!(g[i].is_nan());
        }
    }

    #[test]
    fn mmff_angle_bend_cosine_clip_and_near_parallel_floor() {
        let term = contribution(&[(0, 1, 2, 60.0, 1.0, 0)]);
        for sign in [1.0, -1.0] {
            let coordinates = [
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
            // Fixed native scalar witness: raw and normalized cosine overshoot
            // by one ULP in either orientation. Both contribution paths clip.
            let v = Point3 {
                x: 0.0001,
                y: 0.0001,
                z: 0.0001,
            };
            let dist = (v.x * v.x + v.y * v.y + v.z * v.z).sqrt();
            let raw = v.dot_product(&v) / (dist * dist);
            let normalized = v.divided(dist);
            assert_eq!(raw, 1.0000000000000002);
            assert_eq!(normalized.dot_product(&normalized), 1.0000000000000002);
            assert!(energy(&term, &coordinates).is_finite());
            assert_eq!(gradient(&term, &coordinates, &[0.0; 9]), [0.0; 9]);
        }
        let near = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0e-10, 0.0];
        close_gradient(
            &gradient(&term, &near, &[0.0; 9]),
            &[
                0.0,
                2.454295504600424,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                -2.454295504600424,
                0.0,
            ],
        );
    }

    #[test]
    fn mmff_angle_bend_scalar_control_nan_signed_zero_and_floor() {
        let payload = f64::from_bits(0x7ff8000000000042);
        for (input, expected) in [
            (payload, 0x7ff8000000000042),
            (f64::INFINITY, 1.0_f64.to_bits()),
            (f64::NEG_INFINITY, (-1.0_f64).to_bits()),
            (-0.0, 0x8000000000000000),
            (0.0, 0),
        ] {
            let mut value = input;
            clip_to_one(&mut value);
            assert_eq!(value.to_bits(), expected);
        }
        for input in [payload, f64::INFINITY, f64::NEG_INFINITY, 1.0, -1.0] {
            assert_eq!(source_sin_theta(input).to_bits(), 0x3e45798ee2308c3a);
        }
        assert_eq!(source_sin_theta(0.0), 1.0);
        assert_eq!(source_sin_theta(-0.0), 1.0);
        assert!(source_sin_theta(f64::from_bits(1.0_f64.to_bits() - 1)) > 1.0e-8);
        assert_eq!(
            source_sin_theta(f64::from_bits(1.0_f64.to_bits() + 1)),
            1.0e-8
        );
    }

    #[test]
    fn mmff_angle_bend_signed_endpoint_storage_and_errors() {
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
        for (indices, argument, index) in [
            ([32767, 1, 2], ForceFieldIndexArgument::I, 32767),
            ([32768, 1, 2], ForceFieldIndexArgument::I, 4294934528),
            ([0, 32768, 2], ForceFieldIndexArgument::J, 4294934528),
            ([0, 1, 32768], ForceFieldIndexArgument::J, 4294934528),
            ([65535, 1, 2], ForceFieldIndexArgument::I, u32::MAX),
        ] {
            let term = contribution(&[(indices[0], indices[1], indices[2], 60.0, 1.0, 0)]);
            let expected = ForceFieldKernelError::IndexOutOfRange {
                argument,
                index,
                upper_bound: 3,
            };
            // Only first two coordinate rows exist: invalid third must return
            // its distance error before coordinate or gradient addressing.
            let coordinates = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
            let mut cache = [-1.0; 6];
            let mut context = EvaluationContext::for_test(&coordinates, &mut cache, 3);
            assert_eq!(term.get_energy(&mut context), Err(expected.clone()));
            let mut g = [7.0; 6];
            assert_eq!(term.get_grad(&mut context, &mut g), Err(expected));
            assert_eq!(g, [7.0; 6]);
        }
    }

    #[test]
    fn mmff_angle_bend_wrapped_endpoints_and_collapsed_distinct() {
        let term = contribution(&[(65536, 65537, 65538, 60.0, 1.0, 0)]);
        assert_eq!(term.at1_idxs, [0]);
        assert_eq!(term.at2_idxs, [1]);
        assert_eq!(term.at3_idxs, [2]);
        close(energy(&term, &RIGHT), 15.597723721027004);
        close_gradient(&gradient(&term, &RIGHT, &[0.0; 9]), &NONLINEAR_GRAD);
        let collapsed = contribution(&[(0, 65536, 1, 60.0, 1.0, 0)]);
        assert_eq!(collapsed.at1_idxs, [0]);
        assert_eq!(collapsed.at2_idxs, [0]);
        assert!(energy(&collapsed, &RIGHT).is_nan());
        // Source g points only to rows (0,0,1); row2 is untouched.
        let g = gradient(&collapsed, &RIGHT, &[0.0; 9]);
        assert!(g[..6].iter().all(|v| v.is_nan()));
        assert_eq!(g[6..], [0.0; 3]);
    }

    #[test]
    fn mmff_angle_bend_signed_term_count_shared_range() {
        for (length, end, iterations, last) in [
            (0_usize, 0_i32, 0, None),
            (1, 1, 1, Some(0)),
            (2147483647, 2147483647, 2147483647, Some(2147483646)),
            (2147483648, -2147483648, 0, None),
            (4294967295, -1, 0, None),
            (4294967296, 0, 0, None),
            (4294967297, 1, 1, Some(0)),
        ] {
            let indices = source_term_indices(length);
            assert_eq!(indices.start, 0);
            assert_eq!(indices.end, end);
            assert_eq!(indices.len(), iterations);
            assert_eq!(indices.clone().next(), last.map(|_| 0));
            assert_eq!(indices.last(), last);
        }
    }

    #[test]
    fn mmff_angle_bend_packed_storage_word_boundaries() {
        let mut rows = [[0.0; 3]; 3];
        let mut field = ForceField::new(3);
        for row in &mut rows {
            field.positions_mut().push(row);
        }
        let mut term = AngleBendContrib::new(&field);
        assert_eq!(term.linear_bits.capacity(), 0);
        let mut expected_energy = 0.0;
        for i in 0..129 {
            let linh = if i % 2 == 0 { 1 } else { 0 };
            term.add_term(
                field.positions(),
                0,
                1,
                2,
                &MmffAngle {
                    theta0: 60.0,
                    ka: 1.0,
                },
                &central(linh),
            )
            .unwrap();
            expected_energy += if i % 2 == 0 {
                143.9325
            } else {
                15.597723721027004
            };
            if let Some((_, capacity)) = [
                (1, 64),
                (63, 64),
                (64, 64),
                (65, 128),
                (128, 128),
                (129, 256),
            ]
            .iter()
            .find(|(n, _)| *n == i + 1)
            {
                assert_eq!(term.linear_bits.capacity() * 64, *capacity);
            }
        }
        assert_eq!(
            term.linear_bits,
            [0x5555555555555555, 0x5555555555555555, 1]
        );
        assert_eq!(term.linear_bits.len() * size_of::<u64>(), 24);
        assert!((energy(&term, &RIGHT) - expected_energy).abs() < 1.0e-8);
        for i in 0..129 {
            assert_eq!(term.is_linear(i), i % 2 == 0);
        }
    }

    #[test]
    fn mmff_angle_bend_clone_and_trait_copy_independence() {
        let term = contribution(&[(0, 1, 2, 60.0, 1.0, 0), (2, 1, 0, 60.0, 1.0, 1)]);
        let copied = term.copy();
        let mut clone = term.clone();
        assert_ne!(term.at1_idxs.as_ptr(), clone.at1_idxs.as_ptr());
        assert_ne!(term.at2_idxs.as_ptr(), clone.at2_idxs.as_ptr());
        assert_ne!(term.at3_idxs.as_ptr(), clone.at3_idxs.as_ptr());
        assert_ne!(term.linear_bits.as_ptr(), clone.linear_bits.as_ptr());
        assert_ne!(term.theta0.as_ptr(), clone.theta0.as_ptr());
        assert_ne!(term.ka.as_ptr(), clone.ka.as_ptr());
        clone.ka[0] = 5.0;
        clone.theta0[1] = 90.0;
        clone.linear_bits[0] ^= 1;
        clone.at1_idxs[0] = 2;
        clone.at3_idxs[0] = 0;
        clone.at2_idxs.push(0);
        assert_eq!(term.ka, [1.0; 2]);
        assert_eq!(term.theta0, [60.0; 2]);
        assert_eq!(term.linear_bits, [2]);
        assert_eq!(term.at1_idxs, [0, 2]);
        assert_eq!(term.at2_idxs, [1, 1]);
        assert_eq!(term.at3_idxs, [2, 0]);
        close(energy(&term, &RIGHT), 15.597723721027004 + 143.9325);
        assert_eq!(energy(&term, &RIGHT), energy(&*copied, &RIGHT));
        assert_eq!(
            gradient(&term, &RIGHT, &[0.0; 9]),
            gradient(&*copied, &RIGHT, &[0.0; 9])
        );
    }

    #[test]
    fn mmff_angle_bend_distance_cache_and_partial_update_order() {
        let good = contribution(&[(0, 1, 2, 60.0, 1.0, 0)]);
        let mut cache = [-1.0; 6];
        {
            let mut context = EvaluationContext::for_test(&RIGHT, &mut cache, 3);
            assert_eq!(context.distance(0, 1), Ok(1.0));
            close(good.get_energy(&mut context).unwrap(), 15.597723721027004);
            let mut g = [0.0; 9];
            good.get_grad(&mut context, &mut g).unwrap();
            close_gradient(&g, &NONLINEAR_GRAD);
        }
        assert_eq!(cache, [-1.0, 1.0, -1.0, -1.0, 1.0, -1.0]);
        let bad = contribution(&[(0, 1, 32768, 60.0, 1.0, 0)]);
        let expected = ForceFieldKernelError::IndexOutOfRange {
            argument: ForceFieldIndexArgument::J,
            index: 4294934528,
            upper_bound: 3,
        };
        for evaluate_gradient in [false, true] {
            let mut cache = [-1.0; 6];
            let mut g = [9.0; 6];
            {
                let mut context = EvaluationContext::for_test(&RIGHT[..6], &mut cache, 3);
                if evaluate_gradient {
                    assert_eq!(bad.get_grad(&mut context, &mut g), Err(expected.clone()));
                } else {
                    assert_eq!(bad.get_energy(&mut context), Err(expected.clone()));
                }
            }
            assert_eq!(cache, [-1.0, 1.0, -1.0, -1.0, -1.0, -1.0]);
            assert_eq!(g, [9.0; 6]);
        }
        let partial = contribution(&[(0, 1, 2, 60.0, 1.0, 0), (0, 1, 32768, 60.0, 1.0, 0)]);
        let mut cache = [-1.0; 6];
        let mut g = [0.0; 9];
        let mut context = EvaluationContext::for_test(&RIGHT, &mut cache, 3);
        assert_eq!(partial.get_grad(&mut context, &mut g), Err(expected));
        close_gradient(&g, &NONLINEAR_GRAD);
    }

    #[test]
    fn mmff_angle_bend_real_field_initialize_copy_and_rebind() {
        let mut p1 = [1.0, 0.0, 0.0];
        let mut p2 = [0.0; 3];
        let mut p3 = [0.0, 1.0, 0.0];
        let mut field = ForceField::new(3);
        field.positions_mut().push(&mut p1);
        field.positions_mut().push(&mut p2);
        field.positions_mut().push(&mut p3);
        let mut term = AngleBendContrib::new(&field);
        term.add_term(
            field.positions(),
            0,
            1,
            2,
            &MmffAngle {
                theta0: 60.0,
                ka: 1.0,
            },
            &central(0),
        )
        .unwrap();
        field.add_contribution(Box::new(term));
        assert_eq!(
            cf3d_bld_b05_calc_energy(&mut field, &RIGHT),
            Err(ForceFieldKernelError::NotInitialized)
        );
        field.initialize().unwrap();
        close(field.calc_energy_current(None).unwrap(), 15.597723721027004);
        let mut copied = field.copy();
        assert!(copied.positions().is_empty());
        assert_eq!(
            cf3d_bld_b05_calc_energy(&mut copied, &RIGHT),
            Err(ForceFieldKernelError::NotInitialized)
        );
        let mut c1 = [1.0, 0.0, 0.0];
        let mut c2 = [0.0; 3];
        let mut c3 = [0.0, 1.0, 0.0];
        copied.positions_mut().push(&mut c1);
        copied.positions_mut().push(&mut c2);
        copied.positions_mut().push(&mut c3);
        copied.initialize().unwrap();
        close(
            cf3d_bld_b05_calc_energy(&mut copied, &RIGHT).unwrap(),
            15.597723721027004,
        );
        let mut g = [0.0; 9];
        cf3d_bld_b05_calc_grad(&mut copied, &RIGHT, &mut g).unwrap();
        close_gradient(&g, &NONLINEAR_GRAD);
        let changed = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        assert!(cf3d_bld_b05_calc_energy(&mut copied, &changed).unwrap() != 15.597723721027004);
        close(field.calc_energy_current(None).unwrap(), 15.597723721027004);
    }
}
