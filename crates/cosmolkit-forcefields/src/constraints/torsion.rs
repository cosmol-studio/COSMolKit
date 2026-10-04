// Copyright (C) 2013-2024 Paolo Tosco and other RDKit contributors
//
// @@ All Rights Reserved @@
// This file is derived from RDKit ForceField/TorsionConstraint.cpp.
// The contents are covered by the terms of the BSD license
// which is included in the file license.txt, found at the root
// of the RDKit source tree.

use super::angle::RAD2DEG;
use super::{
    EvaluationContext, ForceField, ForceFieldContribution, ForceFieldKernelError,
    TorsionIndexArgument,
};
use crate::geometry::{
    Point3, compute_dihedral_from_flat_radians, compute_dihedral_from_flat_with_vectors,
    compute_dihedral_from_position_slices_radians, normalize_angle_deg,
};

#[derive(Clone, Debug, PartialEq)]
pub(super) struct TorsionConstraintContrib {
    at1_idx: u32,
    at2_idx: u32,
    at3_idx: u32,
    at4_idx: u32,
    min_dihedral_deg: f64,
    max_dihedral_deg: f64,
    force_constant: f64,
}

impl TorsionConstraintContrib {
    pub(super) fn new_absolute(
        owner: &ForceField<'_>,
        idx1: u32,
        idx2: u32,
        idx3: u32,
        idx4: u32,
        min_dihedral_deg: f64,
        max_dihedral_deg: f64,
        force_constant: f64,
    ) -> Result<Self, ForceFieldKernelError> {
        // BEGIN RDKIT CPP FUNCTION TorsionConstraintContrib::TorsionConstraintContrib (TorsionConstraint.cpp:76-84)
        // RDKit❗✔️: TorsionConstraintContrib::TorsionConstraintContrib(
        // RDKit❗✔️:     ForceField *owner, unsigned int idx1, unsigned int idx2, unsigned int idx3,
        // RDKit❗✔️:     unsigned int idx4, double minDihedralDeg, double maxDihedralDeg,
        // RDKit❗✔️:     double forceConst) {
        // RDKit❗✔️:   checkPrecondition(owner, idx1, idx2, idx3, idx4, minDihedralDeg,
        // RDKit❗✔️:                     maxDihedralDeg);
        // RDKit❗✔️:   setParameters(owner, idx1, idx2, idx3, idx4, minDihedralDeg, maxDihedralDeg,
        // RDKit❗✔️:                 forceConst);
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION TorsionConstraintContrib::TorsionConstraintContrib
        validate_torsion_precondition(
            owner,
            idx1,
            idx2,
            idx3,
            idx4,
            min_dihedral_deg,
            max_dihedral_deg,
        )?;
        Ok(Self::set_parameters(
            idx1,
            idx2,
            idx3,
            idx4,
            min_dihedral_deg,
            max_dihedral_deg,
            force_constant,
        ))
    }

    pub(super) fn new_relative(
        owner: &ForceField<'_>,
        idx1: u32,
        idx2: u32,
        idx3: u32,
        idx4: u32,
        relative: bool,
        mut min_dihedral_deg: f64,
        mut max_dihedral_deg: f64,
        force_constant: f64,
    ) -> Result<Self, ForceFieldKernelError> {
        // BEGIN RDKIT CPP FUNCTION TorsionConstraintContrib::TorsionConstraintContrib (TorsionConstraint.cpp:86-102)
        // RDKit❗✔️: TorsionConstraintContrib::TorsionConstraintContrib(
        // RDKit❗✔️:     ForceField *owner, unsigned int idx1, unsigned int idx2, unsigned int idx3,
        // RDKit❗✔️:     unsigned int idx4, bool relative, double minDihedralDeg,
        // RDKit❗✔️:     double maxDihedralDeg, double forceConst) {
        // RDKit❗✔️:   checkPrecondition(owner, idx1, idx2, idx3, idx4, minDihedralDeg,
        // RDKit❗✔️:                     maxDihedralDeg);
        // RDKit❗✔️:   if (relative) {
        // RDKit❗✔️:     double dihedral;
        // RDKit❗✔️:     RDKit::ForceFieldsHelper::computeDihedral(owner->positions(), idx1, idx2,
        // RDKit❗✔️:                                               idx3, idx4, &dihedral);
        // RDKit❗✔️:     dihedral *= RAD2DEG;
        // RDKit❗✔️:     minDihedralDeg += dihedral;
        // RDKit❗✔️:     maxDihedralDeg += dihedral;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   setParameters(owner, idx1, idx2, idx3, idx4, minDihedralDeg, maxDihedralDeg,
        // RDKit❗✔️:                 forceConst);
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION TorsionConstraintContrib::TorsionConstraintContrib
        validate_torsion_precondition(
            owner,
            idx1,
            idx2,
            idx3,
            idx4,
            min_dihedral_deg,
            max_dihedral_deg,
        )?;
        if relative {
            let positions = owner.positions();
            let mut dihedral = compute_dihedral_from_position_slices_radians(
                &positions[idx1 as usize][..],
                &positions[idx2 as usize][..],
                &positions[idx3 as usize][..],
                &positions[idx4 as usize][..],
            );
            dihedral *= RAD2DEG;
            min_dihedral_deg += dihedral;
            max_dihedral_deg += dihedral;
        }
        Ok(Self::set_parameters(
            idx1,
            idx2,
            idx3,
            idx4,
            min_dihedral_deg,
            max_dihedral_deg,
            force_constant,
        ))
    }

    fn set_parameters(
        idx1: u32,
        idx2: u32,
        idx3: u32,
        idx4: u32,
        mut min_dihedral_deg: f64,
        mut max_dihedral_deg: f64,
        force_constant: f64,
    ) -> Self {
        // BEGIN RDKIT CPP FUNCTION TorsionConstraintContrib::setParameters (TorsionConstraint.cpp:60-74)
        // RDKit❗✔️: void TorsionConstraintContrib::setParameters(
        // RDKit❗✔️:     ForceField *owner, unsigned int idx1, unsigned int idx2, unsigned int idx3,
        // RDKit❗✔️:     unsigned int idx4, double minDihedralDeg, double maxDihedralDeg,
        // RDKit❗✔️:     double forceConst) {
        // RDKit❗✔️:   dp_forceField = owner;
        // RDKit❗✔️:   d_at1Idx = idx1;
        // RDKit❗✔️:   d_at2Idx = idx2;
        // RDKit❗✔️:   d_at3Idx = idx3;
        // RDKit❗✔️:   d_at4Idx = idx4;
        // RDKit❗✔️:   RDKit::ForceFieldsHelper::normalizeAngleDeg(minDihedralDeg);
        // RDKit❗✔️:   RDKit::ForceFieldsHelper::normalizeAngleDeg(maxDihedralDeg);
        // RDKit❗✔️:   d_minDihedralDeg = minDihedralDeg;
        // RDKit❗✔️:   d_maxDihedralDeg = maxDihedralDeg;
        // RDKit❗✔️:   d_forceConstant = forceConst;
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION TorsionConstraintContrib::setParameters
        // The owner reference is needed only by source-relative construction;
        // the contribution stores no pointer after validating its indices.
        normalize_angle_deg(&mut min_dihedral_deg);
        normalize_angle_deg(&mut max_dihedral_deg);
        Self {
            at1_idx: idx1,
            at2_idx: idx2,
            at3_idx: idx3,
            at4_idx: idx4,
            min_dihedral_deg,
            max_dihedral_deg,
            force_constant,
        }
    }

    fn compute_dihedral_term(&self, dihedral: f64) -> f64 {
        // BEGIN RDKIT CPP FUNCTION TorsionConstraintContrib::computeDihedralTerm (TorsionConstraint.cpp:40-58)
        // RDKit❗✔️: double dihedralTarget = dihedral;
        // RDKit❗✔️: if (!(dihedral > d_minDihedralDeg && dihedral < d_maxDihedralDeg) &&
        // RDKit❗✔️:     !(dihedral > d_minDihedralDeg && d_minDihedralDeg > d_maxDihedralDeg) &&
        // RDKit❗✔️:     !(dihedral < d_maxDihedralDeg && d_minDihedralDeg > d_maxDihedralDeg)) {
        // RDKit❗✔️:   double dihedralMinTarget = dihedral - d_minDihedralDeg;
        // RDKit❗✔️:   RDKit::ForceFieldsHelper::normalizeAngleDeg(dihedralMinTarget);
        // RDKit❗✔️:   double dihedralMaxTarget = dihedral - d_maxDihedralDeg;
        // RDKit❗✔️:   RDKit::ForceFieldsHelper::normalizeAngleDeg(dihedralMaxTarget);
        // RDKit❗✔️:   if (fabs(dihedralMinTarget) < fabs(dihedralMaxTarget)) {
        // RDKit❗✔️:     dihedralTarget = d_minDihedralDeg;
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     dihedralTarget = d_maxDihedralDeg;
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // RDKit❗✔️: double dihedralTerm = dihedral - dihedralTarget;
        // RDKit❗✔️: RDKit::ForceFieldsHelper::normalizeAngleDeg(dihedralTerm);
        // RDKit❗✔️: return dihedralTerm;
        // END RDKIT CPP FUNCTION TorsionConstraintContrib::computeDihedralTerm
        let mut dihedral_target = dihedral;
        if !(dihedral > self.min_dihedral_deg && dihedral < self.max_dihedral_deg)
            && !(dihedral > self.min_dihedral_deg && self.min_dihedral_deg > self.max_dihedral_deg)
            && !(dihedral < self.max_dihedral_deg && self.min_dihedral_deg > self.max_dihedral_deg)
        {
            let mut dihedral_min_target = dihedral - self.min_dihedral_deg;
            normalize_angle_deg(&mut dihedral_min_target);
            let mut dihedral_max_target = dihedral - self.max_dihedral_deg;
            normalize_angle_deg(&mut dihedral_max_target);
            if dihedral_min_target.abs() < dihedral_max_target.abs() {
                dihedral_target = self.min_dihedral_deg;
            } else {
                dihedral_target = self.max_dihedral_deg;
            }
        }
        let mut dihedral_term = dihedral - dihedral_target;
        normalize_angle_deg(&mut dihedral_term);
        dihedral_term
    }

    pub(super) fn get_energy(&self, context: &EvaluationContext<'_>) -> f64 {
        // BEGIN RDKIT CPP FUNCTION TorsionConstraintContrib::getEnergy (TorsionConstraint.cpp:104-115)
        // RDKit❗✔️: double TorsionConstraintContrib::getEnergy(double *pos) const {
        // RDKit❗✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit❗✔️:   PRECONDITION(pos, "bad vector");
        // RDKit❗✔️:   double dihedral;
        // RDKit❗✔️:   RDKit::ForceFieldsHelper::computeDihedral(pos, d_at1Idx, d_at2Idx, d_at3Idx,
        // RDKit❗✔️:                                             d_at4Idx, &dihedral);
        // RDKit❗✔️:   dihedral *= RAD2DEG;
        // RDKit❗✔️:   double dihedralTerm = computeDihedralTerm(dihedral);
        // RDKit❗✔️:   double res = d_forceConstant * dihedralTerm * dihedralTerm;
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION TorsionConstraintContrib::getEnergy
        // The borrowed context carries the source trial coordinates; neither the
        // owner nor its distance cache participates in this energy calculation.
        let mut dihedral = compute_dihedral_from_flat_radians(
            context.coordinates,
            self.at1_idx as usize,
            self.at2_idx as usize,
            self.at3_idx as usize,
            self.at4_idx as usize,
        );
        dihedral *= RAD2DEG;
        let dihedral_term = self.compute_dihedral_term(dihedral);
        let res = self.force_constant * dihedral_term * dihedral_term;
        res
    }

    fn get_grad(
        &self,
        context: &mut EvaluationContext<'_>,
        gradient: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        // BEGIN RDKIT CPP FUNCTION TorsionConstraintContrib::getGrad (TorsionConstraint.cpp:117-164)
        // RDKit✔️✔️: void TorsionConstraintContrib::getGrad(double *pos, double *grad) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:   PRECONDITION(grad, "bad vector");
        // RDKit✔️✔️:   double *g[4] = {&(grad[3 * d_at1Idx]), &(grad[3 * d_at2Idx]),
        // RDKit✔️✔️:                   &(grad[3 * d_at3Idx]), &(grad[3 * d_at4Idx])};
        // RDKit✔️✔️:   RDGeom::Point3D r[4];
        // RDKit✔️✔️:   RDGeom::Point3D t[2];
        // RDKit✔️✔️:   double d[2];
        // RDKit✔️✔️:   double dihedral;
        // RDKit✔️✔️:   RDKit::ForceFieldsHelper::computeDihedral(
        // RDKit✔️✔️:       pos, d_at1Idx, d_at2Idx, d_at3Idx, d_at4Idx, &dihedral, nullptr, r, t, d);
        // RDKit✔️✔️:   dihedral *= RAD2DEG;
        // RDKit✔️✔️:   double dihedralTerm = computeDihedralTerm(dihedral);
        // RDKit✔️✔️:   double dE_dPhi = 2.0 * RAD2DEG * d_forceConstant * dihedralTerm;
        // RDKit✔️✔️:
        // RDKit✔️✔️:   double d23 = dp_forceField->distance(d_at2Idx, d_at3Idx, pos);
        // RDKit✔️✔️:   RDGeom::Point3D r31(pos[3 * d_at3Idx] - pos[3 * d_at1Idx],
        // RDKit✔️✔️:                       pos[3 * d_at3Idx + 1] - pos[3 * d_at1Idx + 1],
        // RDKit✔️✔️:                       pos[3 * d_at3Idx + 2] - pos[3 * d_at1Idx + 2]);
        // RDKit✔️✔️:   RDGeom::Point3D r42(pos[3 * d_at4Idx] - pos[3 * d_at2Idx],
        // RDKit✔️✔️:                       pos[3 * d_at4Idx + 1] - pos[3 * d_at2Idx + 1],
        // RDKit✔️✔️:                       pos[3 * d_at4Idx + 2] - pos[3 * d_at2Idx + 2]);
        // RDKit✔️✔️:   double prefactor = dE_dPhi / d23;
        // RDKit✔️✔️:   RDGeom::Point3D tt[2] = {r[0].crossProduct(r[1]), r[2].crossProduct(r[3])};
        // RDKit✔️✔️:   RDGeom::Point3D dedt[2] = {
        // RDKit✔️✔️:       tt[0].crossProduct(r[2]) / tt[0].lengthSq() * prefactor,
        // RDKit✔️✔️:       tt[1].crossProduct(r[1]) / tt[1].lengthSq() * prefactor};
        // RDKit✔️✔️:   RDGeom::Point3D dedp[4] = {
        // RDKit✔️✔️:       r[2].crossProduct(dedt[0]),
        // RDKit✔️✔️:       r31.crossProduct(dedt[0]) - r[3].crossProduct(dedt[1]),
        // RDKit✔️✔️:       r[0].crossProduct(dedt[0]) + r42.crossProduct(dedt[1]),
        // RDKit✔️✔️:       r[2].crossProduct(dedt[1])};
        // RDKit✔️✔️:   for (unsigned int i = 0; i < 4; ++i) {
        // RDKit✔️✔️:     g[i][0] += dedp[i].x;
        // RDKit✔️✔️:     g[i][1] += dedp[i].y;
        // RDKit✔️✔️:     g[i][2] += dedp[i].z;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION TorsionConstraintContrib::getGrad
        let at1_idx = self.at1_idx as usize;
        let at2_idx = self.at2_idx as usize;
        let at3_idx = self.at3_idx as usize;
        let at4_idx = self.at4_idx as usize;
        let g = [3 * at1_idx, 3 * at2_idx, 3 * at3_idx, 3 * at4_idx];

        let (mut dihedral, r) = compute_dihedral_from_flat_with_vectors(
            context.coordinates,
            at1_idx,
            at2_idx,
            at3_idx,
            at4_idx,
        );
        dihedral *= RAD2DEG;
        let dihedral_term = self.compute_dihedral_term(dihedral);
        let d_e_d_phi = 2.0 * RAD2DEG * self.force_constant * dihedral_term;

        let d23 = context.distance(self.at2_idx, self.at3_idx)?;
        let r31 = Point3 {
            x: context.coordinates[3 * at3_idx] - context.coordinates[3 * at1_idx],
            y: context.coordinates[3 * at3_idx + 1] - context.coordinates[3 * at1_idx + 1],
            z: context.coordinates[3 * at3_idx + 2] - context.coordinates[3 * at1_idx + 2],
        };
        let r42 = Point3 {
            x: context.coordinates[3 * at4_idx] - context.coordinates[3 * at2_idx],
            y: context.coordinates[3 * at4_idx + 1] - context.coordinates[3 * at2_idx + 1],
            z: context.coordinates[3 * at4_idx + 2] - context.coordinates[3 * at2_idx + 2],
        };
        let prefactor = d_e_d_phi / d23;
        let tt = [r[0].cross_product(&r[1]), r[2].cross_product(&r[3])];
        let dedt = [
            tt[0]
                .cross_product(&r[2])
                .divided(tt[0].length_sq())
                .scaled(prefactor),
            tt[1]
                .cross_product(&r[1])
                .divided(tt[1].length_sq())
                .scaled(prefactor),
        ];
        let dedp = [
            r[2].cross_product(&dedt[0]),
            Point3::difference(&r31.cross_product(&dedt[0]), &r[3].cross_product(&dedt[1])),
            Point3::sum(&r[0].cross_product(&dedt[0]), &r42.cross_product(&dedt[1])),
            r[2].cross_product(&dedt[1]),
        ];
        for i in 0..4 {
            gradient[g[i]] += dedp[i].x;
            gradient[g[i] + 1] += dedp[i].y;
            gradient[g[i] + 2] += dedp[i].z;
        }
        Ok(())
    }
}

impl ForceFieldContribution for TorsionConstraintContrib {
    fn get_energy(
        &self,
        context: &mut EvaluationContext<'_>,
    ) -> Result<f64, ForceFieldKernelError> {
        Ok(TorsionConstraintContrib::get_energy(self, context))
    }

    fn get_grad(
        &self,
        context: &mut EvaluationContext<'_>,
        gradient: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        TorsionConstraintContrib::get_grad(self, context, gradient)
    }

    fn copy(&self) -> Box<dyn ForceFieldContribution> {
        Box::new(self.clone())
    }
}

fn validate_torsion_precondition(
    owner: &ForceField<'_>,
    idx1: u32,
    idx2: u32,
    idx3: u32,
    idx4: u32,
    min_dihedral_deg: f64,
    max_dihedral_deg: f64,
) -> Result<(), ForceFieldKernelError> {
    // BEGIN RDKIT CPP HELPER TorsionConstraintContrib::checkPrecondition (TorsionConstraint.cpp:27-38)
    // RDKit❗✔️: inline void checkPrecondition(const ForceField *owner, unsigned int idx1,
    // RDKit❗✔️:                           unsigned int idx2, unsigned int idx3,
    // RDKit❗✔️:                           unsigned int idx4, double minDihedralDeg,
    // RDKit❗✔️:                           double maxDihedralDeg) {
    // RDKit❗✔️:   PRECONDITION(owner, "bad owner");
    // RDKit❗✔️:   PRECONDITION(!(minDihedralDeg > maxDihedralDeg),
    // RDKit❗✔️:                "minDihedralDeg must be <= maxDihedralDeg");
    // RDKit❗✔️:   URANGE_CHECK(idx1, owner->positions().size());
    // RDKit❗✔️:   URANGE_CHECK(idx2, owner->positions().size());
    // RDKit❗✔️:   URANGE_CHECK(idx3, owner->positions().size());
    // RDKit❗✔️:   URANGE_CHECK(idx4, owner->positions().size());
    // RDKit❗✔️: }
    // END RDKIT CPP HELPER TorsionConstraintContrib::checkPrecondition
    // `owner: &ForceField` makes the non-null precondition structural. Preserve
    // the source's negated strict comparison, which accepts NaN bounds.
    if min_dihedral_deg > max_dihedral_deg {
        return Err(ForceFieldKernelError::TorsionBoundsOrder);
    }
    let point_count = owner.positions().len();
    validate_torsion_index(idx1, point_count, TorsionIndexArgument::First)?;
    validate_torsion_index(idx2, point_count, TorsionIndexArgument::Second)?;
    validate_torsion_index(idx3, point_count, TorsionIndexArgument::Third)?;
    validate_torsion_index(idx4, point_count, TorsionIndexArgument::Fourth)
}

fn validate_torsion_index(
    index: u32,
    upper_bound: usize,
    argument: TorsionIndexArgument,
) -> Result<(), ForceFieldKernelError> {
    // BEGIN RDKIT CPP HELPER URANGE_CHECK (RDGeneral/Invariant.h:141-151)
    // RDKit❗✔️: #define URANGE_CHECK(x, hi)                                                 \
    // RDKit❗✔️:   if (x >= (hi)) {                                                          \
    // RDKit❗✔️:     std::stringstream errstr;                                               \
    // RDKit❗✔️:     errstr << x << " < " << hi;                                             \
    // RDKit❗✔️:     Invar::Invariant inv("Range Error", #x, errstr.str().c_str(), __FILE__, \
    // RDKit❗✔️:                          __LINE__);                                         \
    // RDKit❗✔️:     BOOST_LOG(rdErrorLog) << "\n\n****\n" << inv << "****\n" << std::endl;  \
    // RDKit❗✔️:     throw inv;                                                              \
    // RDKit❗✔️:   }
    // END RDKIT CPP HELPER URANGE_CHECK
    if index as usize >= upper_bound {
        return Err(ForceFieldKernelError::TorsionIndexOutOfRange {
            argument,
            index,
            upper_bound,
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::{ForceField, ForceFieldKernelError, TorsionIndexArgument};
    use super::TorsionConstraintContrib;

    const NEGATIVE_NINETY_DEG: [f64; 12] =
        [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0];
    const POSITIVE_NINETY_DEG: [f64; 12] =
        [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, -1.0];
    const ZERO_DEG: [f64; 12] = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0];
    const NEGATIVE_ONE_EIGHTY_DEG: [f64; 12] =
        [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 0.0];
    const NEGATIVE_179_DEG: [f64; 12] = [
        1.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        1.0,
        0.0,
        -0.9998476951563913,
        1.0,
        0.01745240643728351,
    ];
    const POSITIVE_179_DEG: [f64; 12] = [
        1.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        1.0,
        0.0,
        -0.9998476951563913,
        1.0,
        -0.01745240643728351,
    ];
    const NEGATIVE_160_DEG: [f64; 12] = [
        1.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        1.0,
        0.0,
        -0.9396926207859084,
        1.0,
        0.3420201433256687,
    ];
    const POSITIVE_160_DEG: [f64; 12] = [
        1.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        1.0,
        0.0,
        -0.9396926207859084,
        1.0,
        -0.3420201433256687,
    ];

    fn assert_source_value(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 1.0e-9,
            "expected {expected}, got {actual}"
        );
    }

    fn owner_rows(coordinates: &[f64; 12]) -> [Vec<f64>; 4] {
        std::array::from_fn(|index| coordinates[index * 3..index * 3 + 3].to_vec())
    }

    fn force_field_for_rows(points: &mut [Vec<f64>; 4]) -> ForceField<'_> {
        let mut force_field = ForceField::new(3);
        force_field
            .positions_mut()
            .extend(points.iter_mut().map(Vec::as_mut_slice));
        force_field
    }

    fn energy(
        force_field: &mut ForceField<'_>,
        contribution: &TorsionConstraintContrib,
        coordinates: &[f64; 12],
    ) -> f64 {
        let context = force_field.evaluation_context(coordinates);
        contribution.get_energy(&context)
    }

    fn gradient(
        force_field: &mut ForceField<'_>,
        contribution: &TorsionConstraintContrib,
        coordinates: &[f64; 12],
        gradient: &mut [f64; 12],
    ) -> Result<(), ForceFieldKernelError> {
        let mut context = force_field.evaluation_context(coordinates);
        contribution.get_grad(&mut context, gradient)
    }

    fn assert_gradient(actual: &[f64; 12], expected: &[f64; 12]) {
        for (index, (&actual_component, &expected_component)) in
            actual.iter().zip(expected).enumerate()
        {
            assert!(
                (actual_component - expected_component).abs() <= 1.0e-9,
                "component {index}: expected {expected_component}, got {actual_component}; actual vector {actual:?}; expected vector {expected:?}"
            );
        }
    }

    #[test]
    fn cf3d_f22_absolute_wrap_interval_and_crossing_energy() {
        // The pinned signed-dihedral formula gives -90/+90 for the reflected
        // final point and -180 for the final point at (-1, 1, 0).
        let mut rows = owner_rows(&NEGATIVE_NINETY_DEG);
        let mut force_field = force_field_for_rows(&mut rows);
        let contribution =
            TorsionConstraintContrib::new_absolute(&force_field, 0, 1, 2, 3, -190.0, 190.0, 2.0)
                .unwrap();

        assert_source_value(contribution.min_dihedral_deg, 170.0);
        assert_source_value(contribution.max_dihedral_deg, -170.0);
        assert_source_value(contribution.compute_dihedral_term(179.0), 0.0);
        assert_source_value(contribution.compute_dihedral_term(-179.0), 0.0);
        assert_source_value(contribution.compute_dihedral_term(90.0), -80.0);
        assert_source_value(contribution.compute_dihedral_term(-90.0), 80.0);
        // At the equidistant 0-degree input the source's strict `<` selects
        // the maximum endpoint, yielding +170 degrees.
        assert_source_value(contribution.compute_dihedral_term(0.0), 170.0);
        assert_source_value(energy(&mut force_field, &contribution, &ZERO_DEG), 57_800.0);
        assert_source_value(
            energy(&mut force_field, &contribution, &NEGATIVE_ONE_EIGHTY_DEG),
            0.0,
        );
    }

    #[test]
    fn cf3d_f22_absolute_signed_endpoints_and_inactive_region() {
        let mut rows = owner_rows(&NEGATIVE_NINETY_DEG);
        let mut force_field = force_field_for_rows(&mut rows);
        let contribution =
            TorsionConstraintContrib::new_absolute(&force_field, 0, 1, 2, 3, -10.0, 10.0, 2.0)
                .unwrap();

        assert_source_value(contribution.compute_dihedral_term(-10.0), 0.0);
        assert_source_value(contribution.compute_dihedral_term(10.0), 0.0);
        assert_source_value(contribution.compute_dihedral_term(-30.0), -20.0);
        assert_source_value(contribution.compute_dihedral_term(30.0), 20.0);
        assert_source_value(contribution.compute_dihedral_term(0.0), 0.0);
        assert_source_value(
            energy(&mut force_field, &contribution, &NEGATIVE_NINETY_DEG),
            12_800.0,
        );
        assert_source_value(
            energy(&mut force_field, &contribution, &POSITIVE_NINETY_DEG),
            12_800.0,
        );
        assert_source_value(energy(&mut force_field, &contribution, &ZERO_DEG), 0.0);
    }

    #[test]
    fn cf3d_f22_relative_true_offsets_once_and_false_keeps_absolute_bounds() {
        let mut negative_rows = owner_rows(&NEGATIVE_NINETY_DEG);
        let mut negative_force_field = force_field_for_rows(&mut negative_rows);
        let relative_negative = TorsionConstraintContrib::new_relative(
            &negative_force_field,
            0,
            1,
            2,
            3,
            true,
            -10.0,
            10.0,
            2.0,
        )
        .unwrap();
        let absolute_negative = TorsionConstraintContrib::new_relative(
            &negative_force_field,
            0,
            1,
            2,
            3,
            false,
            -10.0,
            10.0,
            2.0,
        )
        .unwrap();

        assert_source_value(relative_negative.min_dihedral_deg, -100.0);
        assert_source_value(relative_negative.max_dihedral_deg, -80.0);
        assert_source_value(absolute_negative.min_dihedral_deg, -10.0);
        assert_source_value(absolute_negative.max_dihedral_deg, 10.0);
        assert_source_value(
            energy(
                &mut negative_force_field,
                &relative_negative,
                &NEGATIVE_NINETY_DEG,
            ),
            0.0,
        );
        assert_source_value(
            energy(
                &mut negative_force_field,
                &relative_negative,
                &POSITIVE_NINETY_DEG,
            ),
            57_800.0,
        );
        assert_source_value(
            energy(
                &mut negative_force_field,
                &absolute_negative,
                &NEGATIVE_NINETY_DEG,
            ),
            12_800.0,
        );

        let mut positive_rows = owner_rows(&POSITIVE_NINETY_DEG);
        let mut positive_force_field = force_field_for_rows(&mut positive_rows);
        let relative_positive = TorsionConstraintContrib::new_relative(
            &positive_force_field,
            0,
            1,
            2,
            3,
            true,
            -10.0,
            10.0,
            2.0,
        )
        .unwrap();
        assert_source_value(relative_positive.min_dihedral_deg, 80.0);
        assert_source_value(relative_positive.max_dihedral_deg, 100.0);
        assert_source_value(
            energy(
                &mut positive_force_field,
                &relative_positive,
                &POSITIVE_NINETY_DEG,
            ),
            0.0,
        );
    }

    #[test]
    fn cf3d_f22_constructor_keeps_source_error_order_roles_and_nan_bounds() {
        let mut rows = owner_rows(&NEGATIVE_NINETY_DEG);
        let force_field = force_field_for_rows(&mut rows);
        assert_eq!(
            TorsionConstraintContrib::new_absolute(&force_field, 4, 1, 2, 3, 10.0, 0.0, 2.0)
                .unwrap_err(),
            ForceFieldKernelError::TorsionBoundsOrder
        );
        for (indices, argument) in [
            ([4, 1, 2, 3], TorsionIndexArgument::First),
            ([0, 4, 2, 3], TorsionIndexArgument::Second),
            ([0, 1, 4, 3], TorsionIndexArgument::Third),
            ([0, 1, 2, 4], TorsionIndexArgument::Fourth),
        ] {
            assert_eq!(
                TorsionConstraintContrib::new_absolute(
                    &force_field,
                    indices[0],
                    indices[1],
                    indices[2],
                    indices[3],
                    -10.0,
                    10.0,
                    2.0,
                )
                .unwrap_err(),
                ForceFieldKernelError::TorsionIndexOutOfRange {
                    argument,
                    index: 4,
                    upper_bound: 4,
                }
            );
        }
        let nan_bounds =
            TorsionConstraintContrib::new_absolute(&force_field, 0, 1, 2, 3, f64::NAN, 10.0, 2.0)
                .unwrap();
        assert!(nan_bounds.min_dihedral_deg.is_nan());
        assert_source_value(nan_bounds.max_dihedral_deg, 10.0);
        assert_eq!(
            ForceFieldKernelError::TorsionBoundsOrder.source_expression(),
            Some("!(minDihedralDeg > maxDihedralDeg)".to_owned())
        );
    }

    #[test]
    fn cf3d_f23_wrapped_interval_gradients_cover_both_sides_and_violations() {
        // The source wrap predicates in TorsionConstraint.cpp:46-48 retain
        // angles near both +180 and -180 for normalized bounds [170,-170].
        // For fixed r0=(1,0,0), r1=(0,1,0), and unit r3=(a,0,b),
        // getGrad's pinned cross-product formula gives the deltas below.
        let mut rows = owner_rows(&NEGATIVE_NINETY_DEG);
        let mut force_field = force_field_for_rows(&mut rows);
        let contribution =
            TorsionConstraintContrib::new_absolute(&force_field, 0, 1, 2, 3, -190.0, 190.0, 2.0)
                .unwrap();
        force_field.initialize().unwrap();

        for coordinates in [&NEGATIVE_179_DEG, &POSITIVE_179_DEG] {
            let mut actual = [1.5; 12];
            gradient(&mut force_field, &contribution, coordinates, &mut actual).unwrap();
            assert_gradient(&actual, &[1.5; 12]);
        }

        let mut negative_actual = [0.0; 12];
        gradient(
            &mut force_field,
            &contribution,
            &NEGATIVE_160_DEG,
            &mut negative_actual,
        )
        .unwrap();
        assert_gradient(
            &negative_actual,
            &[
                0.0,
                0.0,
                2291.831180523293,
                0.0,
                0.0,
                -2291.831180523293,
                -783.8524288408132,
                0.0,
                -2153.6168484247955,
                783.8524288408132,
                0.0,
                2153.6168484247955,
            ],
        );

        let mut positive_actual = [0.0; 12];
        gradient(
            &mut force_field,
            &contribution,
            &POSITIVE_160_DEG,
            &mut positive_actual,
        )
        .unwrap();
        assert_gradient(
            &positive_actual,
            &[
                0.0,
                0.0,
                -2291.831180523293,
                0.0,
                0.0,
                2291.831180523293,
                -783.8524288408132,
                0.0,
                2153.6168484247955,
                783.8524288408132,
                0.0,
                -2153.6168484247955,
            ],
        );
    }

    #[test]
    fn cf3d_f23_gradient_adds_source_forces_through_trait_dispatch() {
        // TorsionConstraint.cpp:117-164 adds each x/y/z derivative to the
        // existing row values. This goes through the real ForceField trait
        // dispatch, with fixed -90-degree geometry and source bounds [-10,10].
        let mut rows = owner_rows(&NEGATIVE_NINETY_DEG);
        let mut force_field = force_field_for_rows(&mut rows);
        let contribution =
            TorsionConstraintContrib::new_absolute(&force_field, 0, 1, 2, 3, -10.0, 10.0, 2.0)
                .unwrap();
        force_field
            .contributions
            .push(Box::new(contribution.clone()));
        force_field.initialize().unwrap();

        let mut actual = [0.25; 12];
        force_field
            .calc_grad(&NEGATIVE_NINETY_DEG, &mut actual)
            .unwrap();
        let derivative = -18_334.649444186343;
        let mut expected = [0.25; 12];
        expected[2] += derivative;
        expected[5] -= derivative;
        expected[6] -= derivative;
        expected[9] += derivative;
        assert_gradient(&actual, &expected);
    }

    #[test]
    fn cf3d_f23_degenerate_normals_and_zero_central_distance_keep_ieee_writes() {
        // Pinned getGrad has neither the legacy d23 guard nor either legacy
        // zero-normal guard: both cases continue through component division.
        const COLLINEAR: [f64; 12] = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 3.0, 0.0, 0.0];
        const ZERO_CENTRAL_DISTANCE: [f64; 12] =
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0];
        let mut rows = owner_rows(&NEGATIVE_NINETY_DEG);
        let mut force_field = force_field_for_rows(&mut rows);
        let contribution =
            TorsionConstraintContrib::new_absolute(&force_field, 0, 1, 2, 3, 10.0, 20.0, 2.0)
                .unwrap();
        force_field.initialize().unwrap();

        for coordinates in [&COLLINEAR, &ZERO_CENTRAL_DISTANCE] {
            let mut actual = [0.0; 12];
            gradient(&mut force_field, &contribution, coordinates, &mut actual).unwrap();
            for atom_gradient in actual.chunks_exact(3) {
                assert!(
                    atom_gradient.iter().any(|component| !component.is_finite()),
                    "source divisions must reach atom gradient row {atom_gradient:?}"
                );
            }
        }
    }

    #[test]
    fn cf3d_f23_distance_error_precedes_gradient_mutation() {
        let mut rows = owner_rows(&NEGATIVE_NINETY_DEG);
        let mut force_field = force_field_for_rows(&mut rows);
        let contribution =
            TorsionConstraintContrib::new_absolute(&force_field, 0, 1, 2, 3, -10.0, 10.0, 2.0)
                .unwrap();
        let mut actual = [3.5; 12];
        assert_eq!(
            gradient(
                &mut force_field,
                &contribution,
                &NEGATIVE_NINETY_DEG,
                &mut actual,
            ),
            Err(ForceFieldKernelError::NotInitialized)
        );
        assert_eq!(actual, [3.5; 12]);
    }
}
