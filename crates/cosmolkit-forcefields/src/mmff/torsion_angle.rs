//! Uninstalled complete private MMFF TorsionAngle proposal; source-only.
//! RDKit351f8f378f8ad6bbd517980c38896e66bf907af8, BSD license retained in packet.
//! Requires ROOT approval of all three paths, fixed tests and Native safety.
use super::{
    angle_bend::point_at,
    bond_stretch::source_term_indices,
    numerical::{
        MmffGradientError, calc_torsion_cos_phi, calc_torsion_energy, calc_torsion_force_constant,
        calc_torsion_grad,
    },
    params::MmffTor,
};
use crate::{
    geometry::{Point3, compute_dihedral_from_flat},
    kernel::{
        EvaluationContext, ForceField, ForceFieldContribution, ForceFieldKernelError,
        TorsionIndexArgument,
    },
    uff::params::is_double_zero,
};

// Header storage anchors, verbatim:
// RDKit❗✔️:   std::vector<int16_t> d_at1Idx;
// RDKit❗✔️:   std::vector<int16_t> d_at2Idx;
// RDKit❗✔️:   std::vector<int16_t> d_at3Idx;
// RDKit❗✔️:   std::vector<int16_t> d_at4Idx;
// RDKit❗✔️:   std::vector<double> d_V1;
// RDKit❗✔️:   std::vector<double> d_V2;
// RDKit❗✔️:   std::vector<double> d_V3;
#[derive(Clone, Debug)]
pub(super) struct TorsionAngleContrib {
    at1_idxs: Vec<i16>,
    at2_idxs: Vec<i16>,
    at3_idxs: Vec<i16>,
    at4_idxs: Vec<i16>,
    v1: Vec<f64>,
    v2: Vec<f64>,
    v3: Vec<f64>,
}

impl Default for TorsionAngleContrib {
    fn default() -> Self {
        // RDKit❗✔️: TorsionAngleContrib() {}
        // Behavior: seven empty detached arrays; kernel supplies non-null owner
        // access at evaluation. Raw C++ nullptr evaluation is outside this boundary.
        // Complexity: seven empty Vecs, no allocations or owner cloning.
        Self {
            at1_idxs: Vec::new(),
            at2_idxs: Vec::new(),
            at3_idxs: Vec::new(),
            at4_idxs: Vec::new(),
            v1: Vec::new(),
            v2: Vec::new(),
            v3: Vec::new(),
        }
    }
}

fn checked_row(index: i16, buffer_len: usize) -> Result<usize, ForceFieldKernelError> {
    // RDKit❗✔️:     const int16_t at1Idx = d_at1Idx[i];
    // RDKit❗✔️:     double *g[4] = {&(grad[3 * at1Idx]), &(grad[3 * at2Idx]),
    // RDKit❗✔️:                     &(grad[3 * at3Idx]), &(grad[3 * at4Idx])};
    // Behavior: proposed Native policy outside C++ defined address domain.
    // Reject negative narrowed indices/missing complete rows; ROOT decision required.
    // i16*3+2 fits i32. No typed C++ exception or Unsupported claim.
    // Complexity: constant arithmetic/bounds, O(1), no allocation.
    if index >= 0 && ((3 * i32::from(index) + 2) as usize) < buffer_len {
        Ok(index as usize)
    } else {
        Err(ForceFieldKernelError::TorsionAddressOutsideDefinedSource { index, buffer_len })
    }
}

fn source_sin_term(v1: f64, v2: f64, v3: f64, cos_phi: f64) -> f64 {
    // RDKit❗✔️:     double sinPhiSq = 1.0 - cosPhi * cosPhi;
    // RDKit❗✔️:     double sinPhi = ((sinPhiSq > 0.0) ? sqrt(sinPhiSq) : 0.0);
    // RDKit❗✔️:     double sin2Phi = 2.0 * sinPhi * cosPhi;
    // RDKit❗✔️:     double sin3Phi = 3.0 * sinPhi - 4.0 * sinPhi * sinPhiSq;
    // RDKit❗✔️:     // dE/dPhi is independent of cartesians:
    // RDKit❗✔️:     double dE_dPhi = 0.5 * (-(d_V1[i]) * sinPhi + 2.0 * d_V2[i] * sin2Phi -
    // RDKit❗✔️:                             3.0 * d_V3[i] * sin3Phi);
    // RDKit❗✔️:     // FIX: use a tolerance here
    // RDKit❗✔️:     // this is hacky, but it's per the
    // RDKit❗✔️:     // recommendation from Niketic and Rasmussen:
    // RDKit❗✔️:     double sinTerm =
    // RDKit❗✔️:         -dE_dPhi * (isDoubleZero(sinPhi) ? (1.0 / cosPhi) : (1.0 / sinPhi));
    // Behavior: exact source scalar association and strict 1e-10 reciprocal branch;
    // sin3 uses the original square even after conditional sqrt returns0.
    // Complexity: constant scalar work and one conditional, O(1), no allocation.
    let sin_phi_sq = 1.0 - cos_phi * cos_phi;
    let sin_phi = if sin_phi_sq > 0.0 {
        sin_phi_sq.sqrt()
    } else {
        0.0
    };
    let sin2_phi = 2.0 * sin_phi * cos_phi;
    let sin3_phi = 3.0 * sin_phi - 4.0 * sin_phi * sin_phi_sq;
    let d_e_d_phi = 0.5 * (-v1 * sin_phi + 2.0 * v2 * sin2_phi - 3.0 * v3 * sin3_phi);
    -d_e_d_phi
        * if is_double_zero(sin_phi) {
            1.0 / cos_phi
        } else {
            1.0 / sin_phi
        }
}

impl TorsionAngleContrib {
    pub(super) fn new(_owner: &ForceField<'_>) -> Self {
        // RDKit❗✔️: TorsionAngleContrib::TorsionAngleContrib(ForceField *owner) {
        // RDKit❗✔️:   PRECONDITION(owner, "bad owner");
        // RDKit❗✔️:   dp_forceField = owner;
        // RDKit❗✔️: }
        // Behavior: non-null borrowed owner boundary follows the existing kernel
        // context design; no retained pointer or new lifecycle authority.
        // Complexity: empty sevenVec construction, O(1), no heap/owner clone.
        Self::default()
    }

    pub(super) fn add_term(
        &mut self,
        positions: &[&mut [f64]],
        idx1: u32,
        idx2: u32,
        idx3: u32,
        idx4: u32,
        params: Option<&MmffTor>,
    ) -> Result<(), ForceFieldKernelError> {
        // RDKit❗✔️: void TorsionAngleContrib::addTerm(
        // RDKit❗✔️:     unsigned int idx1, unsigned int idx2, unsigned int idx3, unsigned int idx4,
        // RDKit❗✔️:     const ForceFields::MMFF::MMFFTor *mmffTorParams) {
        // RDKit❗✔️:   PRECONDITION((idx1 != idx2) && (idx1 != idx3) && (idx1 != idx4) &&
        // RDKit❗✔️:                    (idx2 != idx3) && (idx2 != idx4) && (idx3 != idx4),
        // RDKit❗✔️:                "degenerate points");
        // RDKit❗✔️:   URANGE_CHECK(idx1, dp_forceField->positions().size());
        // RDKit❗✔️:   URANGE_CHECK(idx2, dp_forceField->positions().size());
        // RDKit❗✔️:   URANGE_CHECK(idx3, dp_forceField->positions().size());
        // RDKit❗✔️:   URANGE_CHECK(idx4, dp_forceField->positions().size());
        // RDKit❗✔️:
        // RDKit❗✔️:   d_at1Idx.push_back(idx1);
        // RDKit❗✔️:   d_at2Idx.push_back(idx2);
        // RDKit❗✔️:   d_at3Idx.push_back(idx3);
        // RDKit❗✔️:   d_at4Idx.push_back(idx4);
        // RDKit❗✔️:   d_V1.push_back(mmffTorParams->V1);
        // RDKit❗✔️:   d_V2.push_back(mmffTorParams->V2);
        // RDKit❗✔️:   d_V3.push_back(mmffTorParams->V3);
        // RDKit❗✔️: }
        // Behavior: all raw distinct checks, then four ordered raw bounds precede
        // nullable parameter access. There is no source nullable PRECONDITION.
        // Proposed Native null rejection after source checks/before mutation applies
        // only to undefined dereference. Valid rows preserve seven push order and
        // u32 -> i16 low16 interpretation including valid positive narrowed aliases.
        // Complexity: fixed checks and seven amortizedO1 pushes, 32payload bytes/term.
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
        let params = params.ok_or(ForceFieldKernelError::TorsionParametersOutsideDefinedSource)?;
        self.at1_idxs.push(idx1 as i16);
        self.at2_idxs.push(idx2 as i16);
        self.at3_idxs.push(idx3 as i16);
        self.at4_idxs.push(idx4 as i16);
        let (v1, v2, v3) = calc_torsion_force_constant(params);
        self.v1.push(v1);
        self.v2.push(v2);
        self.v3.push(v3);
        Ok(())
    }
}

impl ForceFieldContribution for TorsionAngleContrib {
    fn get_energy(
        &self,
        context: &mut EvaluationContext<'_>,
    ) -> Result<f64, ForceFieldKernelError> {
        // RDKit❗✔️: double TorsionAngleContrib::getEnergy(double *pos) const {
        // RDKit❗✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit❗✔️:   PRECONDITION(pos, "bad vector");
        // RDKit❗✔️:
        // RDKit❗✔️:   const int numTorsions = d_at1Idx.size();
        // RDKit❗✔️:   double totalEnergy = 0.0;
        // RDKit❗✔️:   for (int i = 0; i < numTorsions; ++i) {
        // RDKit❗✔️:     const int16_t at1Idx = d_at1Idx[i];
        // RDKit❗✔️:     const int16_t at2Idx = d_at2Idx[i];
        // RDKit❗✔️:     const int16_t at3Idx = d_at3Idx[i];
        // RDKit❗✔️:     const int16_t at4Idx = d_at4Idx[i];
        // RDKit❗✔️:
        // RDKit❗✔️:     RDGeom::Point3D iPoint(pos[3 * at1Idx], pos[3 * at1Idx + 1],
        // RDKit❗✔️:                            pos[3 * at1Idx + 2]);
        // RDKit❗✔️:     RDGeom::Point3D jPoint(pos[3 * at2Idx], pos[3 * at2Idx + 1],
        // RDKit❗✔️:                            pos[3 * at2Idx + 2]);
        // RDKit❗✔️:     RDGeom::Point3D kPoint(pos[3 * at3Idx], pos[3 * at3Idx + 1],
        // RDKit❗✔️:                            pos[3 * at3Idx + 2]);
        // RDKit❗✔️:     RDGeom::Point3D lPoint(pos[3 * at4Idx], pos[3 * at4Idx + 1],
        // RDKit❗✔️:                            pos[3 * at4Idx + 2]);
        // RDKit❗✔️:
        // RDKit❗✔️:     totalEnergy += Utils::calcTorsionEnergy(
        // RDKit❗✔️:         d_V1[i], d_V2[i], d_V3[i],
        // RDKit❗✔️:         Utils::calcTorsionCosPhi(iPoint, jPoint, kPoint, lPoint));
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return totalEnergy;
        // RDKit❗✔️: }
        // Behavior: kernel references model non-null owner/vector. Signed32 term
        // count loop, four physical points, exact unique cos/E helpers and ordered
        // source total. Does not read context distances or change coordinates/cache.
        // Complexity: O(source count) with O(1) stack points; no allocation/buffering.
        let pos = context.coordinates();
        let mut total_energy = 0.0;
        for term in source_term_indices(self.at1_idxs.len()) {
            let i = term as usize;
            let a = self.at1_idxs[i];
            let b = self.at2_idxs[i];
            let c = self.at3_idxs[i];
            let d = self.at4_idxs[i];
            let p1 = point_at(pos, checked_row(a, pos.len())?);
            let p2 = point_at(pos, checked_row(b, pos.len())?);
            let p3 = point_at(pos, checked_row(c, pos.len())?);
            let p4 = point_at(pos, checked_row(d, pos.len())?);
            total_energy += calc_torsion_energy(
                self.v1[i],
                self.v2[i],
                self.v3[i],
                calc_torsion_cos_phi(&p1, &p2, &p3, &p4),
            );
        }
        Ok(total_energy)
    }

    fn get_grad(
        &self,
        context: &mut EvaluationContext<'_>,
        gradient: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        // RDKit❗✔️: void TorsionAngleContrib::getGrad(double *pos, double *grad) const {
        // RDKit❗✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit❗✔️:   PRECONDITION(pos, "bad vector");
        // RDKit❗✔️:   PRECONDITION(grad, "bad vector");
        // RDKit❗✔️:   double d[2];
        // RDKit❗✔️:   double cosPhi;
        // RDKit❗✔️:
        // RDKit❗✔️:   const int numTorsions = d_at1Idx.size();
        // RDKit❗✔️:   for (int i = 0; i < numTorsions; ++i) {
        // RDKit❗✔️:     const int16_t at1Idx = d_at1Idx[i];
        // RDKit❗✔️:     const int16_t at2Idx = d_at2Idx[i];
        // RDKit❗✔️:     const int16_t at3Idx = d_at3Idx[i];
        // RDKit❗✔️:     const int16_t at4Idx = d_at4Idx[i];
        // RDKit❗✔️:
        // RDKit❗✔️:     double *g[4] = {&(grad[3 * at1Idx]), &(grad[3 * at2Idx]),
        // RDKit❗✔️:                     &(grad[3 * at3Idx]), &(grad[3 * at4Idx])};
        // RDKit❗✔️:
        // RDKit❗✔️:     RDGeom::Point3D r[4];
        // RDKit❗✔️:     RDGeom::Point3D t[2];
        // RDKit❗✔️:
        // RDKit❗✔️:     RDKit::ForceFieldsHelper::computeDihedral(
        // RDKit❗✔️:         pos, at1Idx, at2Idx, at3Idx, at4Idx, nullptr, &cosPhi, r, t, d);
        // RDKit❗✔️:     double sinPhiSq = 1.0 - cosPhi * cosPhi;
        // RDKit❗✔️:     double sinPhi = ((sinPhiSq > 0.0) ? sqrt(sinPhiSq) : 0.0);
        // RDKit❗✔️:     double sin2Phi = 2.0 * sinPhi * cosPhi;
        // RDKit❗✔️:     double sin3Phi = 3.0 * sinPhi - 4.0 * sinPhi * sinPhiSq;
        // RDKit❗✔️:     // dE/dPhi is independent of cartesians:
        // RDKit❗✔️:     double dE_dPhi = 0.5 * (-(d_V1[i]) * sinPhi + 2.0 * d_V2[i] * sin2Phi -
        // RDKit❗✔️:                             3.0 * d_V3[i] * sin3Phi);
        // RDKit❗✔️:     // FIX: use a tolerance here
        // RDKit❗✔️:     // this is hacky, but it's per the
        // RDKit❗✔️:     // recommendation from Niketic and Rasmussen:
        // RDKit❗✔️:     double sinTerm =
        // RDKit❗✔️:         -dE_dPhi * (isDoubleZero(sinPhi) ? (1.0 / cosPhi) : (1.0 / sinPhi));
        // RDKit❗✔️:
        // RDKit❗✔️:     Utils::calcTorsionGrad(r, t, d, g, sinTerm, cosPhi);
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // Behavior: four gradient addresses precede geometry coordinate accesses;
        // unique floored-dihedral helper and scalar source reciprocal rule precede
        // unique numerical helper's twelve sequential writes, including aliased rows.
        // No singular-geometry early return. Native checks do not undo prior terms.
        // Complexity: O(source count), O(1) stack arrays; as_chunks_mut borrows rows
        // in place, no gradient heap buffer or copy. Unused flat tail untouched.
        let pos = context.coordinates();
        let buffer_len = gradient.len();
        let mut d = [0.0; 2];
        let mut cos_phi = 0.0;
        for term in source_term_indices(self.at1_idxs.len()) {
            let i = term as usize;
            let indices = [
                self.at1_idxs[i],
                self.at2_idxs[i],
                self.at3_idxs[i],
                self.at4_idxs[i],
            ];
            let rows = [
                checked_row(indices[0], buffer_len)?,
                checked_row(indices[1], buffer_len)?,
                checked_row(indices[2], buffer_len)?,
                checked_row(indices[3], buffer_len)?,
            ];
            for index in indices {
                checked_row(index, pos.len())?;
            }
            let mut r = [Point3::default(); 4];
            let mut t = [Point3::default(); 2];
            compute_dihedral_from_flat(
                pos,
                rows[0],
                rows[1],
                rows[2],
                rows[3],
                None,
                Some(&mut cos_phi),
                Some(&mut r),
                Some(&mut t),
                Some(&mut d),
            );
            let sin_term = source_sin_term(self.v1[i], self.v2[i], self.v3[i], cos_phi);
            let (gradient_rows, _unused_tail) = gradient.as_chunks_mut::<3>();
            calc_torsion_grad(&r, &t, &d, gradient_rows, rows, sin_term, cos_phi).map_err(
                |error| match error {
                    MmffGradientError::GradientRowOutOfRange { index, .. } => {
                        ForceFieldKernelError::TorsionAddressOutsideDefinedSource {
                            index: index as i16,
                            buffer_len,
                        }
                    }
                },
            )?;
        }
        Ok(())
    }

    fn copy(&self) -> Box<dyn ForceFieldContribution> {
        // RDKit❗✔️: TorsionAngleContrib *copy() const override {
        // RDKit❗✔️:     return new TorsionAngleContrib(*this);
        // RDKit❗✔️:   }
        // Behavior: deep independent copy of all seven Vecs; owner supplied by
        // existing kernel on use. No whole ForceField/position/cache cloning.
        // Complexity: O(stored terms), exactly seven cloned arrays and boxed self.
        Box::new(self.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernel::{cf3d_bld_b05_calc_energy, cf3d_bld_b05_calc_grad};
    // All expected values/tolerances are review proposals. No executed Rust oracle.
    // Pinned source351f8f378f8ad6bbd517980c38896e66bf907af8, full closure in packet.
    const COEFF: MmffTor = MmffTor {
        v1: 1.0,
        v2: 2.0,
        v3: 3.0,
    };
    const POS: [f64; 12] = [
        0.0,
        1.0,
        0.0,
        0.0,
        0.0,
        0.0,
        1.0,
        0.0,
        0.0,
        1.0,
        0.7071067811865476,
        0.7071067811865476,
    ];
    const G: [f64; 12] = [
        0.0,
        0.0,
        1.5355339059327375,
        0.0,
        0.0,
        -1.5355339059327375,
        0.0,
        -1.085786437626906,
        1.0857864376269055,
        0.0,
        1.085786437626906,
        -1.0857864376269055,
    ];
    const E: f64 = 2.292893218813452;
    fn close(a: f64, b: f64) {
        assert!((a - b).abs() < 1e-10, "{a:?} != {b:?}");
    }
    fn close_grad(a: &[f64], b: &[f64]) {
        assert_eq!(a.len(), b.len());
        for (&a, &b) in a.iter().zip(b) {
            close(a, b);
        }
    }
    fn lengths(t: &TorsionAngleContrib) -> [usize; 7] {
        [
            t.at1_idxs.len(),
            t.at2_idxs.len(),
            t.at3_idxs.len(),
            t.at4_idxs.len(),
            t.v1.len(),
            t.v2.len(),
            t.v3.len(),
        ]
    }
    fn append(
        t: &mut TorsionAngleContrib,
        f: &ForceField<'_>,
        i: [u32; 4],
        v: Option<&MmffTor>,
    ) -> Result<(), ForceFieldKernelError> {
        t.add_term(f.positions(), i[0], i[1], i[2], i[3], v)
    }
    fn make(indices: &[[u32; 4]], coefficients: &[MmffTor]) -> TorsionAngleContrib {
        assert_eq!(indices.len(), coefficients.len());
        let n = indices
            .iter()
            .flatten()
            .max()
            .map_or(0, |&i| i as usize + 1);
        let mut rows = vec![[0.0; 3]; n];
        let mut f = ForceField::new(3);
        for row in &mut rows {
            f.positions_mut().push(row);
        }
        let mut t = TorsionAngleContrib::new(&f);
        for (i, c) in indices.iter().zip(coefficients) {
            append(&mut t, &f, *i, Some(c)).unwrap();
        }
        t
    }
    fn one() -> TorsionAngleContrib {
        make(&[[0, 1, 2, 3]], &[COEFF])
    }
    fn energy(t: &dyn ForceFieldContribution, pos: &[f64]) -> f64 {
        let mut cache = [];
        t.get_energy(&mut EvaluationContext::for_test(
            pos,
            &mut cache,
            (pos.len() / 3) as u32,
        ))
        .unwrap()
    }
    fn gradient(t: &dyn ForceFieldContribution, pos: &[f64], initial: &[f64]) -> Vec<f64> {
        let mut cache = [];
        let mut g = initial.to_vec();
        t.get_grad(
            &mut EvaluationContext::for_test(pos, &mut cache, (pos.len() / 3) as u32),
            &mut g,
        )
        .unwrap();
        g
    }

    // Coverage proposal: both constructor bodies; seven empty arrays/no alloc; no-term +0 energy; no-term grad preserves exact bytes; non-null context modeled boundary
    #[test]
    fn mmff_torsion_angle_default_and_owner_empty() {
        for t in [
            TorsionAngleContrib::default(),
            TorsionAngleContrib::new(&ForceField::new(3)),
        ] {
            assert_eq!(lengths(&t), [0; 7]);
            assert_eq!(
                [
                    t.at1_idxs.capacity(),
                    t.at2_idxs.capacity(),
                    t.at3_idxs.capacity(),
                    t.at4_idxs.capacity(),
                    t.v1.capacity(),
                    t.v2.capacity(),
                    t.v3.capacity()
                ],
                [0; 7]
            );
            assert_eq!(energy(&t, &[]).to_bits(), 0);
            assert_eq!(
                gradient(&t, &[], &[-0.0, 1.0])
                    .iter()
                    .map(|v| v.to_bits())
                    .collect::<Vec<_>>(),
                [(-0.0_f64).to_bits(), 1.0_f64.to_bits()]
            );
        }
    }

    // Coverage proposal: all6 distinct pairs, multierror precedence with absent params and invalid raw indices; no mutation
    #[test]
    fn mmff_torsion_angle_six_distinct_before_bounds_and_null() {
        let f = ForceField::new(3);
        let mut t = TorsionAngleContrib::new(&f);
        for (a, b) in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)] {
            let mut i = [10, 11, 12, 13];
            i[b] = i[a];
            for params in [None, Some(&COEFF)] {
                assert_eq!(
                    append(&mut t, &f, i, params),
                    Err(ForceFieldKernelError::TorsionDegeneratePoints)
                );
                assert_eq!(lengths(&t), [0; 7]);
            }
        }
    }

    // Coverage proposal: each ordered raw bounds before nullable dereference, u32max, full role/index/bounds errors
    #[test]
    fn mmff_torsion_angle_four_raw_bounds_order_before_null() {
        let mut rows = [[0.0; 3]; 4];
        let mut f = ForceField::new(3);
        for row in &mut rows {
            f.positions_mut().push(row);
        }
        let mut t = TorsionAngleContrib::new(&f);
        for (i, role, value) in [
            ([4, 5, 6, 7], TorsionIndexArgument::First, 4),
            ([0, 4, 5, 6], TorsionIndexArgument::Second, 4),
            ([0, 1, 4, 5], TorsionIndexArgument::Third, 4),
            ([0, 1, 2, u32::MAX], TorsionIndexArgument::Fourth, u32::MAX),
        ] {
            for v in [None, Some(&COEFF)] {
                assert_eq!(
                    append(&mut t, &f, i, v),
                    Err(ForceFieldKernelError::TorsionIndexOutOfRange {
                        argument: role,
                        index: value,
                        upper_bound: 4
                    })
                );
                assert_eq!(lengths(&t), [0; 7]);
            }
        }
    }

    // Coverage proposal: PROPOSED NATIVE null safety after distinct/bounds; CPP undefined, not precondition parity
    #[test]
    fn mmff_torsion_angle_native_null_after_checks_preserves_arrays() {
        let mut rows = [[0.0; 3]; 4];
        let mut f = ForceField::new(3);
        for row in &mut rows {
            f.positions_mut().push(row);
        }
        let mut t = TorsionAngleContrib::new(&f);
        append(&mut t, &f, [0, 1, 2, 3], Some(&COEFF)).unwrap();
        let before = t.clone();
        assert_eq!(
            append(&mut t, &f, [0, 1, 2, 3], None),
            Err(ForceFieldKernelError::TorsionParametersOutsideDefinedSource)
        );
        assert_eq!(lengths(&t), [1; 7]);
        assert_eq!(t.at1_idxs, before.at1_idxs);
        assert_eq!(t.at2_idxs, before.at2_idxs);
        assert_eq!(t.at3_idxs, before.at3_idxs);
        assert_eq!(t.at4_idxs, before.at4_idxs);
        assert_eq!(t.v1, before.v1);
        assert_eq!(t.v2, before.v2);
        assert_eq!(t.v3, before.v3);
    }

    // Coverage proposal: all seven source arrays, independent copied coefficients, multiple vector growths
    #[test]
    fn mmff_torsion_angle_seven_arrays_growth_and_copied_parameters() {
        let mut rows = [[0.0; 3]; 4];
        let mut f = ForceField::new(3);
        for row in &mut rows {
            f.positions_mut().push(row);
        }
        let mut t = TorsionAngleContrib::new(&f);
        let mut c = COEFF;
        for n in 0..65 {
            c.v1 = n as f64;
            c.v2 = -(n as f64);
            c.v3 = n as f64 + 0.25;
            append(&mut t, &f, [3, 2, 1, 0], Some(&c)).unwrap();
        }
        c.v1 = 999.0;
        c.v2 = 999.0;
        c.v3 = 999.0;
        assert_eq!(lengths(&t), [65; 7]);
        for n in 0..65 {
            assert_eq!(
                [t.at1_idxs[n], t.at2_idxs[n], t.at3_idxs[n], t.at4_idxs[n]],
                [3, 2, 1, 0]
            );
            assert_eq!(
                [t.v1[n], t.v2[n], t.v3[n]],
                [n as f64, -(n as f64), n as f64 + 0.25]
            );
        }
        assert_eq!(c.v1, 999.0);
    }

    // Coverage proposal: E and every12 independent analytic derivative, coordinate preservation
    #[test]
    fn mmff_torsion_angle_independent_phi45_energy_all12_grad() {
        let coords = [
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            1.0,
            0.7071067811865476,
            0.7071067811865476,
        ];
        let t = make(
            &[[0, 1, 2, 3]],
            &[MmffTor {
                v1: 1.0,
                v2: 2.0,
                v3: 3.0,
            }],
        );
        close(energy(&t, &coords), 2.292893218813452);
        close_grad(
            &gradient(&t, &coords, &[0.0; 12]),
            &[
                0.0,
                0.0,
                1.5355339059327375,
                0.0,
                0.0,
                -1.5355339059327375,
                0.0,
                -1.085786437626906,
                1.0857864376269055,
                0.0,
                1.085786437626906,
                -1.0857864376269055,
            ],
        );
        assert_eq!(
            coords,
            [
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                1.0,
                0.7071067811865476,
                0.7071067811865476
            ]
        );
    }

    // Coverage proposal: E and every12 independent analytic derivative, coordinate preservation
    #[test]
    fn mmff_torsion_angle_independent_phi60_energy_all12_grad() {
        let coords = [
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            1.0,
            0.5,
            0.8660254037844386,
        ];
        let t = make(
            &[[0, 1, 2, 3]],
            &[MmffTor {
                v1: 1.0,
                v2: 2.0,
                v3: 3.0,
            }],
        );
        close(energy(&t, &coords), 2.25);
        close_grad(
            &gradient(&t, &coords, &[0.0; 12]),
            &[
                0.0,
                0.0,
                -1.2990381056766565,
                -1.6653345369377343e-16,
                0.0,
                1.2990381056766565,
                1.6653345369377343e-16,
                1.124999999999998,
                -0.6495190528383288,
                0.0,
                -1.124999999999998,
                0.6495190528383288,
            ],
        );
        assert_eq!(
            coords,
            [
                0.0,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
                1.0,
                0.0,
                0.0,
                1.0,
                0.5,
                0.8660254037844386
            ]
        );
    }

    // Coverage proposal: E and every12 independent analytic derivative, coordinate preservation
    #[test]
    fn mmff_torsion_angle_independent_phi90_energy_all12_grad() {
        let coords = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0];
        let t = make(
            &[[0, 1, 2, 3]],
            &[MmffTor {
                v1: 1.0,
                v2: 2.0,
                v3: 3.0,
            }],
        );
        close(energy(&t, &coords), 4.0);
        close_grad(
            &gradient(&t, &coords, &[0.0; 12]),
            &[0.0, 0.0, -4.0, 0.0, 0.0, 4.0, 0.0, 4.0, 0.0, 0.0, -4.0, 0.0],
        );
        assert_eq!(
            coords,
            [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0]
        );
    }

    // Coverage proposal: E and every12 independent analytic derivative, coordinate preservation
    #[test]
    fn mmff_torsion_angle_independent_nonaxis_energy_all12_grad() {
        let coords = [
            2.0, 1.0, -1.0, 0.2, -0.3, 0.4, 1.4, 0.7, 0.8, 2.2, -0.9, 1.7,
        ];
        let t = make(
            &[[0, 1, 2, 3]],
            &[MmffTor {
                v1: 1.0,
                v2: 2.0,
                v3: 3.0,
            }],
        );
        close(energy(&t, &coords), 4.885735553512311);
        close_grad(
            &gradient(&t, &coords, &[0.0; 12]),
            &[
                0.23220239969843842,
                -0.29025299962304796,
                0.029025299962304872,
                0.00837054654529828,
                -0.024399335331696664,
                0.03588669869334693,
                -0.41163472800000833,
                0.3990724350422553,
                0.23722309639438638,
                0.17106178175627146,
                -0.08442010008751066,
                -0.3021350950500383,
            ],
        );
        assert_eq!(
            coords,
            [
                2.0, 1.0, -1.0, 0.2, -0.3, 0.4, 1.4, 0.7, 0.8, 2.2, -0.9, 1.7
            ]
        );
    }

    // Coverage proposal: E and every12 independent analytic derivative, coordinate preservation
    #[test]
    fn mmff_torsion_angle_independent_nonaxis_negative_energy_all12_grad() {
        let coords = [
            2.0, 1.0, -1.0, 0.2, -0.3, 0.4, 1.4, 0.7, 0.8, 2.2, -0.9, 1.7,
        ];
        let t = make(
            &[[0, 1, 2, 3]],
            &[MmffTor {
                v1: -2.0,
                v2: 0.5,
                v3: -3.0,
            }],
        );
        close(energy(&t, &coords), -3.1107895217218613);
        close_grad(
            &gradient(&t, &coords, &[0.0; 12]),
            &[
                0.23615691012458634,
                -0.29519613765573294,
                0.029519613765573334,
                0.008513100686120803,
                -0.024814866894185206,
                0.0364978651771006,
                -0.4186450510016423,
                0.4058688165921982,
                0.24126311152443244,
                0.1739750401909354,
                -0.08585781204227971,
                -0.3072805904671062,
            ],
        );
        assert_eq!(
            coords,
            [
                2.0, 1.0, -1.0, 0.2, -0.3, 0.4, 1.4, 0.7, 0.8, 2.2, -0.9, 1.7
            ]
        );
    }

    // Coverage proposal: supplementary independent finite difference each12 axes h1e-6/tol1e-6; fixed analytic E/G expectations remain primary
    #[test]
    fn mmff_torsion_angle_independent_central_difference_supplement() {
        let t = one();
        let g = gradient(&t, &POS, &[0.0; 12]);
        for axis in 0..12 {
            let mut a = POS;
            let mut b = POS;
            a[axis] += 1e-6;
            b[axis] -= 1e-6;
            assert!(((energy(&t, &a) - energy(&t, &b)) / (2e-6) - g[axis]).abs() < 1e-6);
        }
    }

    // Coverage proposal: cos ±1 endpoints, sin=0 reciprocal1/cos branch, no NaN gradient for finite geometry
    #[test]
    fn mmff_torsion_angle_planar_both_cosine_endpoints() {
        let t = one();
        for y in [-1.0, 1.0] {
            let pos = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, y, 0.0];
            assert_eq!(energy(&t, &pos), if y > 0.0 { 4.0 } else { 0.0 });
            assert_eq!(gradient(&t, &pos, &[0.0; 12]), [0.0; 12]);
        }
    }

    // Coverage proposal: conditional nonpositive/NaN sqrt and sin3 original square, finite sin reciprocal and endpoint branch, infinity/NaN arithmetic
    #[test]
    fn mmff_torsion_angle_source_sine_branch_and_ieee_scalar() {
        for c in [-0.7, 0.0, 0.4, 0.9] {
            close(
                source_sin_term(1.0, 2.0, 3.0, c),
                0.5 * (1.0 - 8.0 * c + 3.0 * (12.0 * c * c - 3.0)),
            );
        }
        // At +/-1 source dE=0 before reciprocal; no derivative-polynomial replacement.
        assert_eq!(source_sin_term(1.0, 2.0, 3.0, 1.0), 0.0);
        assert_eq!(source_sin_term(1.0, 2.0, 3.0, -1.0), 0.0);
        for c in [
            1.00001,
            -1.00001,
            f64::NAN,
            f64::INFINITY,
            f64::NEG_INFINITY,
        ] {
            let value = source_sin_term(1.0, 2.0, 3.0, c);
            if c.is_nan() || c.is_infinite() {
                assert!(value.is_nan());
            } else {
                assert_eq!(value, 0.0);
            }
        }
    }

    // Coverage proposal: original three adjacent cosines near sin1e-6 preserved; all are NONZERO at source strict1e-10; literal reciprocal source witness, no fabricated1e-10 cosine boundary
    #[test]
    fn mmff_torsion_angle_source_sine_strict_threshold_neighbors() {
        // Preserve original cosines near sin1e-6. They use the nonzero1e-10 branch.
        // True1e-10 scalar boundaries and representable cosine endpoints are separate tests.
        let c = (1.0_f64 - 1e-12).sqrt();
        for bits in [c.to_bits() - 1, c.to_bits(), c.to_bits() + 1] {
            let c = f64::from_bits(bits);
            let sq = 1.0 - c * c;
            let s = sq.sqrt();
            let de = 0.5 * (-s + 4.0 * (2.0 * s * c) - 9.0 * (3.0 * s - 4.0 * s * sq));
            let reciprocal = if s < 1e-10 && s > -1e-10 {
                1.0 / c
            } else {
                1.0 / s
            };
            assert_eq!(
                source_sin_term(1.0, 2.0, 3.0, c).to_bits(),
                (-de * reciprocal).to_bits()
            );
        }
        assert!(!is_double_zero(1e-6));
        assert!(!is_double_zero(f64::from_bits(1e-6_f64.to_bits() - 1)));
    }

    // Coverage proposal: separate energy strictepsilon/grad ordered floor; no grad early return; all12source-defined expressions, proposed1e-6abs for up to~1e6 derivatives
    #[test]
    fn mmff_torsion_angle_singular_zero_first() {
        let t = one();
        let coords = [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            1.0,
            0.7071067811865476,
            0.7071067811865476,
        ];
        close(energy(&t, &coords), 4.0);
        let actual = gradient(&t, &coords, &[0.0; 12]);
        let expected = [
            0.0,
            -282842.712474619,
            -282842.712474619,
            -0.0,
            282842.712474619,
            282842.712474619,
            -0.0,
            -0.0,
            -0.0,
            -0.0,
            0.0,
            -0.0,
        ];
        for (&a, &e) in actual.iter().zip(&expected) {
            assert!((a - e).abs() < 1e-6, "{a:?} != {e:?}");
        }
    }

    // Coverage proposal: separate energy strictepsilon/grad ordered floor; no grad early return; all12source-defined expressions, proposed1e-6abs for up to~1e6 derivatives
    #[test]
    fn mmff_torsion_angle_singular_zero_last() {
        let t = one();
        let coords = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        close(energy(&t, &coords), 4.0);
        let actual = gradient(&t, &coords, &[0.0; 12]);
        let expected = [
            -0.0,
            -0.0,
            -0.0,
            0.0,
            -0.0,
            -0.0,
            -0.0,
            399999.99999999994,
            -0.0,
            -0.0,
            -399999.99999999994,
            -0.0,
        ];
        for (&a, &e) in actual.iter().zip(&expected) {
            assert!((a - e).abs() < 1e-6, "{a:?} != {e:?}");
        }
    }

    // Coverage proposal: separate energy strictepsilon/grad ordered floor; no grad early return; all12source-defined expressions, proposed1e-6abs for up to~1e6 derivatives
    #[test]
    fn mmff_torsion_angle_singular_small_norm() {
        let t = one();
        let coords = [
            0.0,
            2e-06,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            1.0,
            1.4142135623730952e-06,
            1.4142135623730952e-06,
        ];
        close(energy(&t, &coords), 2.292893218813452);
        let actual = gradient(&t, &coords, &[0.0; 12]);
        let expected = [
            0.0,
            -55646.299912264374,
            -57964.89574194206,
            0.22258519964905749,
            55646.299912264374,
            57964.89574194206,
            -0.22258519964905749,
            80335.24686580097,
            -1639.4948339959383,
            -0.0,
            -80335.24686580097,
            1639.4948339959383,
        ];
        for (&a, &e) in actual.iter().zip(&expected) {
            assert!((a - e).abs() < 1e-6, "{a:?} != {e:?}");
        }
    }

    // Coverage proposal: separate energy strictepsilon/grad ordered floor; no grad early return; all12source-defined expressions, proposed1e-6abs for up to~1e6 derivatives
    #[test]
    fn mmff_torsion_angle_singular_below_epsilon() {
        let t = one();
        let coords = [
            0.0,
            5e-07,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            1.0,
            1.4142135623730952e-06,
            1.4142135623730952e-06,
        ];
        close(energy(&t, &coords), 2.292893218813452);
        let actual = gradient(&t, &coords, &[0.0; 12]);
        let expected = [
            0.0,
            -56813.425036430286,
            -56955.814572862444,
            0.05574550351318911,
            56813.425036430286,
            56955.814572862444,
            -0.05574550351318911,
            19734.182929112558,
            -402.7384271247461,
            -0.0,
            -19734.182929112558,
            402.7384271247461,
        ];
        for (&a, &e) in actual.iter().zip(&expected) {
            assert!((a - e).abs() < 1e-6, "{a:?} != {e:?}");
        }
    }

    // Coverage proposal: separate energy strictepsilon/grad ordered floor; no grad early return; all12source-defined expressions, proposed1e-6abs for up to~1e6 derivatives
    #[test]
    fn mmff_torsion_angle_singular_at_epsilon() {
        let t = one();
        let coords = [
            0.0,
            1e-06,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            1.0,
            1.4142135623730952e-06,
            1.4142135623730952e-06,
        ];
        close(energy(&t, &coords), 2.292893218813452);
        let actual = gradient(&t, &coords, &[0.0; 12]);
        let expected = [
            0.0,
            -56744.45449861157,
            -57317.63080667836,
            0.11176938007302278,
            56744.45449861157,
            57317.63080667836,
            -0.11176938007302278,
            39719.09171645024,
            -810.5937084989845,
            -0.0,
            -39719.09171645024,
            810.5937084989845,
        ];
        for (&a, &e) in actual.iter().zip(&expected) {
            assert!((a - e).abs() < 1e-6, "{a:?} != {e:?}");
        }
    }

    // Coverage proposal: separate energy strictepsilon/grad ordered floor; no grad early return; all12source-defined expressions, proposed1e-6abs for up to~1e6 derivatives
    #[test]
    fn mmff_torsion_angle_singular_below_floor() {
        let t = one();
        let coords = [
            0.0,
            5e-06,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            1.0,
            1.4142135623730953e-05,
            1.4142135623730953e-05,
        ];
        close(energy(&t, &coords), 2.2928932188134525);
        let actual = gradient(&t, &coords, &[0.0; 12]);
        let expected = [
            0.0,
            -167807.7650307344,
            -223743.68670764583,
            0.8390388251536718,
            167807.7650307344,
            223743.68670764583,
            -0.8390388251536718,
            39552.669529663675,
            -39552.66952966369,
            -0.0,
            -39552.669529663675,
            39552.66952966369,
        ];
        for (&a, &e) in actual.iter().zip(&expected) {
            assert!((a - e).abs() < 1e-6, "{a:?} != {e:?}");
        }
    }

    // Coverage proposal: separate energy strictepsilon/grad ordered floor; no grad early return; all12source-defined expressions, proposed1e-6abs for up to~1e6 derivatives
    #[test]
    fn mmff_torsion_angle_singular_at_floor() {
        let t = one();
        let coords = [
            0.0,
            1e-05,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            1.0,
            1.4142135623730953e-05,
            1.4142135623730953e-05,
        ];
        close(energy(&t, &coords), 2.2928932188134525);
        let actual = gradient(&t, &coords, &[0.0; 12]);
        let expected = [
            -0.0,
            4.821860411516461e-11,
            153553.39059327365,
            -2.4109302057582304e-16,
            -4.821860411516461e-11,
            -153553.39059327365,
            2.4109302057582304e-16,
            -54289.32188134518,
            54289.321881345204,
            0.0,
            54289.32188134518,
            -54289.321881345204,
        ];
        for (&a, &e) in actual.iter().zip(&expected) {
            assert!((a - e).abs() < 1e-6, "{a:?} != {e:?}");
        }
    }

    // Coverage proposal: separate energy strictepsilon/grad ordered floor; no grad early return; all12source-defined expressions, proposed1e-6abs for up to~1e6 derivatives
    #[test]
    fn mmff_torsion_angle_singular_above_floor() {
        let t = one();
        let coords = [
            0.0,
            2e-05,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            1.0,
            1.4142135623730953e-05,
            1.4142135623730953e-05,
        ];
        close(energy(&t, &coords), 2.2928932188134525);
        let actual = gradient(&t, &coords, &[0.0; 12]);
        let expected = [
            -0.0,
            2.4109302057582306e-11,
            76776.69529663683,
            -2.4109302057582304e-16,
            -2.4109302057582306e-11,
            -76776.69529663683,
            2.4109302057582304e-16,
            -54289.32188134518,
            54289.321881345204,
            0.0,
            54289.32188134518,
            -54289.321881345204,
        ];
        for (&a, &e) in actual.iter().zip(&expected) {
            assert!((a - e).abs() < 1e-6, "{a:?} != {e:?}");
        }
    }

    // Coverage proposal: separate energy strictepsilon/grad ordered floor; no grad early return; all12source-defined expressions, proposed1e-6abs for up to~1e6 derivatives
    #[test]
    fn mmff_torsion_angle_singular_true_below_first() {
        let t = one();
        let coords = [
            0.0,
            5e-11,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            1.0,
            0.7071067811865476,
            0.7071067811865476,
        ];
        close(energy(&t, &coords), 4.0);
        let actual = gradient(&t, &coords, &[0.0; 12]);
        let expected = [
            0.0,
            -282843.71245163796,
            -282843.712458709,
            1.4142185622581897e-05,
            282843.71245163796,
            282843.712458709,
            -1.4142185622581897e-05,
            1.0000035354776557e-05,
            -1.0000035354776557e-05,
            -0.0,
            -1.0000035354776557e-05,
            1.0000035354776557e-05,
        ];
        for (&a, &e) in actual.iter().zip(&expected) {
            assert!((a - e).abs() < 1e-6, "{a:?} != {e:?}");
        }
    }

    // Coverage proposal: separate energy strictepsilon/grad ordered floor; no grad early return; all12source-defined expressions, proposed1e-6abs for up to~1e6 derivatives
    #[test]
    fn mmff_torsion_angle_singular_true_at_first() {
        let t = one();
        let coords = [
            0.0,
            1e-10,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            1.0,
            0.7071067811865476,
            0.7071067811865476,
        ];
        close(energy(&t, &coords), 2.292893218813453);
        let actual = gradient(&t, &coords, &[0.0; 12]);
        let expected = [
            0.0,
            -282844.71238269494,
            -282844.7124109794,
            2.8284471238269495e-05,
            282844.71238269494,
            282844.7124109794,
            -2.8284471238269495e-05,
            2.0000141416856237e-05,
            -2.0000141416856237e-05,
            -0.0,
            -2.0000141416856237e-05,
            2.0000141416856237e-05,
        ];
        for (&a, &e) in actual.iter().zip(&expected) {
            assert!((a - e).abs() < 1e-6, "{a:?} != {e:?}");
        }
    }

    // Coverage proposal: separate energy strictepsilon/grad ordered floor; no grad early return; all12source-defined expressions, proposed1e-6abs for up to~1e6 derivatives
    #[test]
    fn mmff_torsion_angle_singular_true_above_first() {
        let t = one();
        let coords = [
            0.0,
            2e-10,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            1.0,
            0.7071067811865476,
            0.7071067811865476,
        ];
        close(energy(&t, &coords), 2.292893218813453);
        let actual = gradient(&t, &coords, &[0.0; 12]);
        let expected = [
            0.0,
            -282846.7121069218,
            -282846.7122200605,
            5.6569342421384365e-05,
            282846.7121069218,
            282846.7122200605,
            -5.6569342421384365e-05,
            4.000056564942494e-05,
            -4.000056564942494e-05,
            -0.0,
            -4.000056564942494e-05,
            4.000056564942494e-05,
        ];
        for (&a, &e) in actual.iter().zip(&expected) {
            assert!((a - e).abs() < 1e-6, "{a:?} != {e:?}");
        }
    }

    // Coverage proposal: separate energy strictepsilon/grad ordered floor; no grad early return; all12source-defined expressions, proposed1e-6abs for up to~1e6 derivatives
    #[test]
    fn mmff_torsion_angle_singular_true_below_last() {
        let t = one();
        let coords = [
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            1.0,
            3.535533905932738e-11,
            3.535533905932738e-11,
        ];
        close(energy(&t, &coords), 4.0);
        let actual = gradient(&t, &coords, &[0.0; 12]);
        let expected = [
            -0.0,
            -0.0,
            -1.4142185622935453e-05,
            1.4142185622581899e-05,
            -0.0,
            1.4142185622935453e-05,
            -1.4142185622581899e-05,
            400001.41418606223,
            -5.000017677388278e-06,
            -0.0,
            -400001.41418606223,
            5.000017677388278e-06,
        ];
        for (&a, &e) in actual.iter().zip(&expected) {
            assert!((a - e).abs() < 1e-6, "{a:?} != {e:?}");
        }
    }

    // Coverage proposal: separate energy strictepsilon/grad ordered floor; no grad early return; all12source-defined expressions, proposed1e-6abs for up to~1e6 derivatives
    #[test]
    fn mmff_torsion_angle_singular_true_at_last() {
        let t = one();
        let coords = [
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            1.0,
            7.071067811865477e-11,
            7.071067811865477e-11,
        ];
        close(energy(&t, &coords), 2.2928932188134525);
        let actual = gradient(&t, &coords, &[0.0; 12]);
        let expected = [
            -0.0,
            -0.0,
            -2.8284471241097944e-05,
            2.8284471238269498e-05,
            -0.0,
            2.8284471241097944e-05,
            -2.8284471238269498e-05,
            400002.8283171246,
            -2.0000141416856237e-05,
            -0.0,
            -400002.8283171246,
            2.0000141416856237e-05,
        ];
        for (&a, &e) in actual.iter().zip(&expected) {
            assert!((a - e).abs() < 1e-6, "{a:?} != {e:?}");
        }
    }

    // Coverage proposal: separate energy strictepsilon/grad ordered floor; no grad early return; all12source-defined expressions, proposed1e-6abs for up to~1e6 derivatives
    #[test]
    fn mmff_torsion_angle_singular_true_above_last() {
        let t = one();
        let coords = [
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            1.0,
            1.4142135623730953e-10,
            1.4142135623730953e-10,
        ];
        close(energy(&t, &coords), 2.2928932188134525);
        let actual = gradient(&t, &coords, &[0.0; 12]);
        let expected = [
            -0.0,
            -0.0,
            -5.656934244401211e-05,
            5.6569342421384365e-05,
            -0.0,
            5.656934244401211e-05,
            -5.6569342421384365e-05,
            400005.6564142482,
            -8.000113129884988e-05,
            -0.0,
            -400005.6564142482,
            8.000113129884988e-05,
        ];
        for (&a, &e) in actual.iter().zip(&expected) {
            assert!((a - e).abs() < 1e-6, "{a:?} != {e:?}");
        }
    }

    // Coverage proposal: remaining zero central arm, both cross-norms0, collinear nonzero arm source paths
    #[test]
    fn mmff_torsion_angle_zero_central_arm_and_both_collinear() {
        let t = one();
        for pos in [
            [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 3.0, 0.0, 0.0],
        ] {
            assert_eq!(energy(&t, &pos), 4.0);
            assert_eq!(gradient(&t, &pos, &[0.0; 12]), [0.0; 12]);
        }
    }

    // Coverage proposal: shared signed32 iteration, source ordered total sum and repeated gradient additions; no parallel reduction
    #[test]
    fn mmff_torsion_angle_repeated_terms_source_total_association() {
        let mut t = one();
        t.v1 = vec![1e16, -1e16, 1.0];
        t.v2 = vec![0.0; 3];
        t.v3 = vec![0.0; 3];
        t.at1_idxs = vec![0; 3];
        t.at2_idxs = vec![1; 3];
        t.at3_idxs = vec![2; 3];
        t.at4_idxs = vec![3; 3];
        let pos = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0];
        assert_eq!(energy(&t, &pos), 1.0);
        let twice = make(&[[0, 1, 2, 3]; 2], &[COEFF; 2]);
        close(energy(&twice, &POS), E + E);
        let g = gradient(&twice, &POS, &[0.0; 12]);
        for (&a, &b) in g.iter().zip(&G) {
            close(a, b + b);
        }
    }

    // Coverage proposal: additive all12 updates, borrowed rows and unused flat tail bits
    #[test]
    fn mmff_torsion_angle_additive_nonzero_gradient_and_unused_tail() {
        let t = one();
        let mut initial = [7.0; 14];
        initial[12] = -0.0;
        initial[13] = f64::from_bits(0x7ff8000000001234);
        let g = gradient(&t, &POS, &initial);
        for i in 0..12 {
            close(g[i], initial[i] + G[i]);
        }
        assert_eq!(g[12].to_bits(), initial[12].to_bits());
        assert_eq!(g[13].to_bits(), initial[13].to_bits());
    }

    // Coverage proposal: physical sign/cosine convention and chain reversal, translated/reflected all12 expected independent symmetry
    #[test]
    fn mmff_torsion_angle_translation_reflection_and_reverse() {
        let t = one();
        let mut translated = POS;
        for row in translated.chunks_mut(3) {
            row[0] += 8.0;
            row[1] -= 4.0;
            row[2] += 2.0;
        }
        close(energy(&t, &translated), E);
        close_grad(&gradient(&t, &translated, &[0.0; 12]), &G);
        let mut mirrored = POS;
        for row in mirrored.chunks_mut(3) {
            row[2] *= -1.0;
        }
        close(energy(&t, &mirrored), E);
        let g = gradient(&t, &mirrored, &[0.0; 12]);
        for i in 0..12 {
            close(g[i], if i % 3 == 2 { -G[i] } else { G[i] });
        }
        let reversed = make(&[[3, 2, 1, 0]], &[COEFF]);
        close(energy(&reversed, &POS), E);
        close_grad(&gradient(&reversed, &POS, &[0.0; 12]), &G);
    }

    // Coverage proposal: negative and copied signed-zero coefficients with no zero-coefficient filtering in private contribution
    #[test]
    fn mmff_torsion_angle_coefficient_negative_zero_and_signed_bits() {
        let t = make(
            &[[0, 1, 2, 3]],
            &[MmffTor {
                v1: -1.0,
                v2: -2.0,
                v3: -3.0,
            }],
        );
        close(energy(&t, &POS), -E);
        let g = gradient(&t, &POS, &[0.0; 12]);
        for i in 0..12 {
            close(g[i], -G[i]);
        }
        let zero = make(
            &[[0, 1, 2, 3]],
            &[MmffTor {
                v1: -0.0,
                v2: 0.0,
                v3: -0.0,
            }],
        );
        assert_eq!(zero.v1[0].to_bits(), (-0.0_f64).to_bits());
        assert_eq!(zero.v3[0].to_bits(), (-0.0_f64).to_bits());
        assert_eq!(energy(&zero, &POS), 0.0);
        assert_eq!(gradient(&zero, &POS, &[0.0; 12]), [0.0; 12]);
    }

    // Coverage proposal: each V field nonfinite, no silent coefficient fallback; exact cos0 energy coefficients multiply strictly positive factors, sign of infinity retained
    #[test]
    fn mmff_torsion_angle_nan_infinite_parameters_no_fallback() {
        for c in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            for slot in 0..3 {
                let mut v = [1.0, 2.0, 3.0];
                v[slot] = c;
                let t = make(
                    &[[0, 1, 2, 3]],
                    &[MmffTor {
                        v1: v[0],
                        v2: v[1],
                        v3: v[2],
                    }],
                );
                let coords = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0];
                let e = energy(&t, &coords);
                if c.is_nan() {
                    assert!(e.is_nan());
                } else {
                    assert_eq!(e, c);
                }
                assert!(
                    gradient(&t, &POS, &[0.0; 12])
                        .iter()
                        .any(|x| x.is_nan() || x.is_infinite())
                );
            }
        }
    }

    // Coverage proposal: every Cartesian input nonfinite; energy NaN vs grad orderedmax/min and subsequent IEEE writes; no unsupported or zero fallback
    #[test]
    fn mmff_torsion_angle_nan_infinite_geometry_distinct_helper_paths() {
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            for axis in 0..12 {
                let mut pos = POS;
                pos[axis] = value;
                assert!(energy(&one(), &pos).is_nan());
                let g = gradient(&one(), &pos, &[0.0; 12]);
                assert!(g.iter().any(|x| x.is_nan()));
            }
        }
    }

    // Coverage proposal: cache no read/write under pathological sentinels; source reads physical coordinate vector; reuse existing cfgtest fixture no new fixture authority
    #[test]
    fn mmff_torsion_angle_cache_bits_and_physical_coordinates() {
        let t = one();
        let mut cache = [
            f64::from_bits(0x7ff800000000beef),
            f64::INFINITY,
            -9.0,
            -0.0,
        ];
        let before = cache.map(f64::to_bits);
        let mut coords = POS;
        for expected in [E, E] {
            {
                let mut ctx =
                    EvaluationContext::for_oop_cache_preservation_test(&coords, &mut cache, 4);
                close(t.get_energy(&mut ctx).unwrap(), expected);
                let mut g = [0.0; 12];
                t.get_grad(&mut ctx, &mut g).unwrap();
                close_grad(&g, &G);
            }
            assert_eq!(cache.map(f64::to_bits), before);
            for row in coords.chunks_mut(3) {
                row[0] += 10.0;
            }
        }
    }

    // Coverage proposal: exact signed16 narrowing including negative/endpoints/aliases; extreme raw scalar only no billion-row validation claim
    #[test]
    fn mmff_torsion_angle_signed16_scalar_width_boundaries() {
        for (raw, expected) in [
            (0u32, 0i16),
            (32767, 32767),
            (32768, -32768),
            (65535, -1),
            (65536, 0),
            (65537, 1),
            (u32::MAX, -1),
        ] {
            assert_eq!(raw as i16, expected);
        }
        assert_eq!((0u32 as i16), (65536u32 as i16));
    }

    // Coverage proposal: actual highraw append bounds and definedpositive narrow address, no Oop32 copy behavior
    #[test]
    fn mmff_torsion_angle_raw_positive_narrowed_high_indices() {
        let t = make(&[[65536, 65537, 65538, 65539]], &[COEFF]);
        assert_eq!(
            [t.at1_idxs[0], t.at2_idxs[0], t.at3_idxs[0], t.at4_idxs[0]],
            [0, 1, 2, 3]
        );
        close(energy(&t, &POS), E);
        close_grad(&gradient(&t, &POS, &[0.0; 12]), &G);
        // append checks real65540 rows, evaluation intentionally reads only narrowed rows.
    }

    // Coverage proposal: six raw-distinct/narrowed-duplicate pairs are permitted; defined alias accumulation, no validation-on-narrowed distinct
    #[test]
    fn mmff_torsion_angle_all_six_narrowed_alias_pairs() {
        for (a, b) in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)] {
            let mut idx = [0, 1, 2, 3];
            idx[b] = idx[a] + 65536;
            let t = make(&[idx], &[COEFF]);
            let mut coords = POS;
            // first/third alias makes central bond geometry collinear; all source
            // aliases have either cross zero (cos0) or equal normals (cos1).
            assert_eq!(energy(&t, &coords), 4.0);
            let g = gradient(&t, &coords, &[0.0; 12]);
            assert_eq!(g, [0.0; 12]);
            coords[11] += 0.0;
            assert_eq!(lengths(&t), [1; 7]);
        }
    }

    // Coverage proposal: collapsed first/second row source sequential += rounding; prohibits aggregated/temporary gradients; hand-reduced sourceJacob
    #[test]
    fn mmff_torsion_angle_alias_sequential_rounding_not_aggregated() {
        let t = make(&[[0, 65536, 1, 2]], &[COEFF]);
        let coords = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0];
        let mut initial = [0.0; 9];
        initial[1] = 0.1;
        let actual = gradient(&t, &coords, &initial);
        let delta = -4.0 * (1.0 / 1e-5);
        let expected: f64 = (0.1 + delta) + (-delta);
        assert_eq!(actual[1].to_bits(), expected.to_bits());
        assert_ne!(expected.to_bits(), 0.1_f64.to_bits());
        for axis in [0, 2, 3, 4, 5, 6, 7, 8] {
            assert_eq!(actual[axis], 0.0);
        }
    }

    // Coverage proposal: seven independent deep arrays and unchanged original; actual clone state
    #[test]
    fn mmff_torsion_angle_deep_clone_all_seven_arrays() {
        let t = one();
        let mut copied = t.clone();
        assert_ne!(t.at1_idxs.as_ptr(), copied.at1_idxs.as_ptr());
        assert_ne!(t.at2_idxs.as_ptr(), copied.at2_idxs.as_ptr());
        assert_ne!(t.at3_idxs.as_ptr(), copied.at3_idxs.as_ptr());
        assert_ne!(t.at4_idxs.as_ptr(), copied.at4_idxs.as_ptr());
        assert_ne!(t.v1.as_ptr(), copied.v1.as_ptr());
        assert_ne!(t.v2.as_ptr(), copied.v2.as_ptr());
        assert_ne!(t.v3.as_ptr(), copied.v3.as_ptr());
        copied.at1_idxs[0] = 9;
        copied.at2_idxs[0] = 8;
        copied.at3_idxs[0] = 7;
        copied.at4_idxs[0] = 6;
        copied.v1[0] = 11.0;
        copied.v2[0] = 12.0;
        copied.v3[0] = 13.0;
        assert_eq!(
            [t.at1_idxs[0], t.at2_idxs[0], t.at3_idxs[0], t.at4_idxs[0]],
            [0, 1, 2, 3]
        );
        assert_eq!([t.v1[0], t.v2[0], t.v3[0]], [1.0, 2.0, 3.0]);
        close(energy(&t, &POS), E);
    }

    // Coverage proposal: headercopy boxed deep clone, data independence and no borrow of original owner/parameter
    #[test]
    fn mmff_torsion_angle_trait_copy_outlives_original() {
        let mut t = one();
        let copied = t.copy();
        t.v1[0] = 100.0;
        t.v2[0] = 200.0;
        t.v3[0] = 300.0;
        drop(t);
        close(energy(&*copied, &POS), E);
        close_grad(&gradient(&*copied, &POS, &[0.0; 12]), &G);
        let empty = TorsionAngleContrib::default().copy();
        assert_eq!(energy(&*empty, &[]), 0.0);
    }

    // Coverage proposal: actual privatekernel initialized callbacks, ForceField copyemptypositions and manually reattached ownrows; no MMFFbuilder
    #[test]
    fn mmff_torsion_angle_real_kernel_initialize_copy_and_rebind() {
        let mut rows = POS
            .chunks_exact(3)
            .map(|r| [r[0], r[1], r[2]])
            .collect::<Vec<_>>();
        let mut f = ForceField::new(3);
        for row in &mut rows {
            f.positions_mut().push(row);
        }
        let mut t = TorsionAngleContrib::new(&f);
        append(&mut t, &f, [0, 1, 2, 3], Some(&COEFF)).unwrap();
        f.add_contribution(Box::new(t));
        f.initialize().unwrap();
        close(f.calc_energy_current(None).unwrap(), E);
        let mut g = [0.0; 12];
        cf3d_bld_b05_calc_grad(&mut f, &POS, &mut g).unwrap();
        close_grad(&g, &G);
        let mut copied = f.copy();
        assert!(copied.positions().is_empty());
        let mut next = POS
            .chunks_exact(3)
            .map(|r| [r[0], r[1], r[2]])
            .collect::<Vec<_>>();
        for row in &mut next {
            copied.positions_mut().push(row);
        }
        copied.initialize().unwrap();
        close(cf3d_bld_b05_calc_energy(&mut copied, &POS).unwrap(), E);
        let mut g2 = [0.0; 12];
        cf3d_bld_b05_calc_grad(&mut copied, &POS, &mut g2).unwrap();
        close_grad(&g2, &G);
        close(f.calc_energy_current(None).unwrap(), E);
    }

    // Coverage proposal: shared count native32 scalar only; actual billions-of-terms allocation UNRUN and not claimed
    #[test]
    fn mmff_torsion_angle_shared_signed32_term_count() {
        for (len, end) in [
            (0usize, 0i32),
            (2, 2),
            (i32::MAX as usize, i32::MAX),
            (i32::MAX as usize + 1, i32::MIN),
            (u32::MAX as usize, -1),
            (u32::MAX as usize + 1, 0),
        ] {
            let range = source_term_indices(len);
            assert_eq!(range.start, 0);
            assert_eq!(range.end, end);
        }
    }

    // Coverage proposal: PROPOSED NATIVE whole-row/negative narrowed bounds, extreme source signed16 multiplication fits32
    #[test]
    fn mmff_torsion_angle_native_signed16_address_bounds() {
        assert_eq!(checked_row(0, 3), Ok(0));
        assert_eq!(checked_row(32767, 98304), Ok(32767));
        for (index, len) in [
            (-1, usize::MAX),
            (i16::MIN, usize::MAX),
            (0, 2),
            (1, 5),
            (32767, 98303),
        ] {
            assert_eq!(
                checked_row(index, len),
                Err(ForceFieldKernelError::TorsionAddressOutsideDefinedSource {
                    index,
                    buffer_len: len
                })
            );
        }
    }

    // Coverage proposal: actual raw accepted then signed16 negative; CPP undefined negative vector domain, explicit Native safety only
    #[test]
    fn mmff_torsion_angle_native_negative_narrowing_actual_append() {
        for (raw, index) in [(32768u32, i16::MIN), (65535, -1)] {
            let t = make(&[[0, 1, 2, raw]], &[COEFF]);
            assert_eq!(t.at4_idxs, [index]);
            let mut cache = [];
            let mut ctx = EvaluationContext::for_test(&POS, &mut cache, 4);
            assert_eq!(
                t.get_energy(&mut ctx),
                Err(ForceFieldKernelError::TorsionAddressOutsideDefinedSource {
                    index,
                    buffer_len: 12
                })
            );
            let mut g = [7.0; 12];
            assert_eq!(
                t.get_grad(&mut ctx, &mut g),
                Err(ForceFieldKernelError::TorsionAddressOutsideDefinedSource {
                    index,
                    buffer_len: 12
                })
            );
            assert_eq!(g, [7.0; 12]);
        }
    }

    // Coverage proposal: PROPOSED NATIVE errors follow original gradientpointer-before-flatcoord ordering; preserves currentterm atomic update
    #[test]
    fn mmff_torsion_angle_native_gradient_addresses_precede_coordinates() {
        let t = one();
        let mut cache = [];
        let mut ctx = EvaluationContext::for_test(&POS[..2], &mut cache, 4);
        let mut grad = [7.0; 11];
        assert_eq!(
            t.get_grad(&mut ctx, &mut grad),
            Err(ForceFieldKernelError::TorsionAddressOutsideDefinedSource {
                index: 3,
                buffer_len: 11
            })
        );
        assert_eq!(grad, [7.0; 11]);
        let mut complete = [7.0; 12];
        assert_eq!(
            t.get_grad(&mut ctx, &mut complete),
            Err(ForceFieldKernelError::TorsionAddressOutsideDefinedSource {
                index: 0,
                buffer_len: 2
            })
        );
        assert_eq!(complete, [7.0; 12]);
    }

    // Coverage proposal: PROPOSED NATIVE each coordinate role and laterterm failure preserves successful prior12writes
    #[test]
    fn mmff_torsion_angle_native_short_coordinates_each_role_and_prior_term() {
        for (len, index) in [(2usize, 0i16), (5, 1), (8, 2), (11, 3)] {
            let t = one();
            let mut cache = [];
            let mut ctx = EvaluationContext::for_test(&POS[..len], &mut cache, 4);
            assert_eq!(
                t.get_energy(&mut ctx),
                Err(ForceFieldKernelError::TorsionAddressOutsideDefinedSource {
                    index,
                    buffer_len: len
                })
            );
        }
        let mut t = make(&[[0, 1, 2, 3]; 2], &[COEFF; 2]);
        t.at4_idxs[1] = -1;
        let mut cache = [];
        let mut ctx = EvaluationContext::for_test(&POS, &mut cache, 4);
        let mut g = [0.0; 12];
        assert_eq!(
            t.get_grad(&mut ctx, &mut g),
            Err(ForceFieldKernelError::TorsionAddressOutsideDefinedSource {
                index: -1,
                buffer_len: 12
            })
        );
        close_grad(&g, &G);
    }

    // Coverage proposal: pinned MMFF/UFF exact strict((-1e-10<x)&&(x<1e-10)), both signs/neighbor bits/zero/nonfinite; original1e-6 inputs retained
    #[test]
    fn mmff_torsion_angle_true_signed_scalar_epsilon_neighbors() {
        let eps = 1e-10_f64;
        let below = f64::from_bits(eps.to_bits() - 1);
        let above = f64::from_bits(eps.to_bits() + 1);
        for x in [0.0, -0.0, below, -below] {
            assert!(is_double_zero(x));
        }
        for x in [
            eps,
            -eps,
            above,
            -above,
            f64::NAN,
            f64::INFINITY,
            f64::NEG_INFINITY,
        ] {
            assert!(!is_double_zero(x));
        }
        // The1e-6 inputs from the original test are retained and explicitly nonzero.
        assert!(!is_double_zero(1e-6));
        assert!(!is_double_zero(f64::from_bits(1e-6_f64.to_bits() - 1)));
    }

    // Coverage proposal: exact next-inner cosine±1 and minimal positive binary64 computed sine2^-26, endpoint-only zero branch; no invented1e-10 cosine neighbor
    #[test]
    fn mmff_torsion_angle_binary64_cosine_representable_sine_domain() {
        let next = f64::from_bits(1.0_f64.to_bits() - 1);
        for c in [-next, next] {
            let sq = 1.0 - c * c;
            let sine = sq.sqrt();
            assert_eq!(sine, 1.4901161193847656e-8);
            assert!(sine > 1e-10);
            assert!(!is_double_zero(sine));
            let de = 0.5 * (-sine + 4.0 * (2.0 * sine * c) - 9.0 * (3.0 * sine - 4.0 * sine * sq));
            assert_eq!(
                source_sin_term(1.0, 2.0, 3.0, c).to_bits(),
                (-de * (1.0 / sine)).to_bits()
            );
        }
        for c in [-1.0_f64, 1.0] {
            assert_eq!(1.0 - c * c, 0.0);
            assert!(is_double_zero(0.0));
            assert_eq!(source_sin_term(1.0, 2.0, 3.0, c), 0.0);
        }
        // No binary64 clamped cosine has positive sine in(0,1e-10).
    }
}
