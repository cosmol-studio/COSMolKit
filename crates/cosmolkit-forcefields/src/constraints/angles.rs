// Copyright (C) 2004-2024 Paolo Tosco and other RDKit contributors
//
// @@ All Rights Reserved @@
// This file is derived from RDKit ForceField/AngleConstraints.cpp.
// The contents are covered by the terms of the BSD license
// which is included in the file license.txt, found at the root
// of the RDKit source tree.

use super::{
    AngleIndexArgument, AngleRangeBound, EvaluationContext, ForceField, ForceFieldContribution,
    ForceFieldKernelError,
    angle::{
        RAD2DEG, angle_degrees, angle_vector_cross_product, angle_vector_difference,
        angle_vector_dot_product, angle_vector_length, angle_vector_length_sq,
        angle_vector_negated, angle_vector_scale, validate_angle_index, validate_angle_range,
    },
};

#[derive(Clone, Copy, Debug, PartialEq)]
struct AngleConstraintContribsParams {
    // BEGIN RDKIT CPP STRUCT ForceFields::AngleConstraintContribsParams (AngleConstraints.h:19-31)
    // RDKit❗✔️: struct AngleConstraintContribsParams {
    // RDKit❗✔️:   unsigned int idx1{0};       //!< index of atom1 of the angle constraint
    // RDKit❗✔️:   unsigned int idx2{0};       //!< index of atom2 of the angle constraint
    // RDKit❗✔️:   unsigned int idx3{0};       //!< index of atom3 of the angle constraint
    // RDKit❗✔️:   double minAngle{0.0};       //!< lower bound of the flat bottom potential
    // RDKit❗✔️:   double maxAngle{0.0};       //!< upper bound of the flat bottom potential
    // RDKit❗✔️:   double forceConstant{1.0};  //!< force constant for angle constraint
    // RDKit❗✔️:   AngleConstraintContribsParams(unsigned int idx1, unsigned int idx2,
    // RDKit❗✔️:                                 unsigned int idx3, double minAngle,
    // RDKit❗✔️:                                 double maxAngle, double forceConstant = 1.0)
    // RDKit❗✔️:       : idx1(idx1),
    // RDKit❗✔️:         idx2(idx2),
    // RDKit❗✔️:         idx3(idx3),
    // RDKit❗✔️:         minAngle(minAngle),
    // RDKit❗✔️:         maxAngle(maxAngle),
    // RDKit❗✔️:         forceConstant(forceConstant) {};
    // RDKit❗✔️: };
    // END RDKIT CPP STRUCT ForceFields::AngleConstraintContribsParams
    idx1: u32,
    idx2: u32,
    idx3: u32,
    min_angle: f64,
    max_angle: f64,
    force_constant: f64,
}

impl Default for AngleConstraintContribsParams {
    fn default() -> Self {
        Self {
            idx1: 0,
            idx2: 0,
            idx3: 0,
            min_angle: 0.0,
            max_angle: 0.0,
            force_constant: 1.0,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct AngleConstraintContribs {
    // RDKit source: ForceField/AngleConstraints.h:91
    // RDKit❗✔️: std::vector<AngleConstraintContribsParams> d_contribs;
    contribs: Vec<AngleConstraintContribsParams>,
}

impl AngleConstraintContribs {
    pub(crate) fn new(_owner: &ForceField<'_>) -> Self {
        // BEGIN RDKIT CPP FUNCTION AngleConstraintContribs::AngleConstraintContribs (AngleConstraints.cpp:23-26)
        // RDKit❗✔️: AngleConstraintContribs::AngleConstraintContribs(ForceField *owner) {
        // RDKit❗✔️:   PRECONDITION(owner, "bad owner");
        // RDKit❗✔️:   dp_forceField = owner;
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION AngleConstraintContribs::AngleConstraintContribs
        // A shared ForceField borrow represents the non-null owner boundary;
        // the group stores no pointer and can safely move with its parameters.
        Self::default()
    }

    pub(crate) fn add_contrib(
        &mut self,
        owner: &ForceField<'_>,
        idx1: u32,
        idx2: u32,
        idx3: u32,
        min_angle_deg: f64,
        max_angle_deg: f64,
        force_constant: f64,
    ) -> Result<(), ForceFieldKernelError> {
        // BEGIN RDKIT CPP FUNCTION AngleConstraintContribs::addContrib (AngleConstraints.cpp:28-39)
        // RDKit❗✔️: void AngleConstraintContribs::addContrib(unsigned int idx1, unsigned int idx2,
        // RDKit❗✔️:                                          unsigned int idx3, double minAngleDeg,
        // RDKit❗✔️:                                          double maxAngleDeg,
        // RDKit❗✔️:                                          double forceConst) {
        // RDKit❗✔️:   RANGE_CHECK(0.0, minAngleDeg, 180.0);
        // RDKit❗✔️:   RANGE_CHECK(0.0, maxAngleDeg, 180.0);
        // RDKit❗✔️:   URANGE_CHECK(idx1, dp_forceField->positions().size());
        // RDKit❗✔️:   URANGE_CHECK(idx2, dp_forceField->positions().size());
        // RDKit❗✔️:   URANGE_CHECK(idx3, dp_forceField->positions().size());
        // RDKit❗✔️:   PRECONDITION(maxAngleDeg >= minAngleDeg,
        // RDKit❗✔️:                "minAngleDeg must be <= maxAngleDeg");
        // RDKit❗✔️:   d_contribs.emplace_back(idx1, idx2, idx3, minAngleDeg, maxAngleDeg,
        // RDKit❗✔️:                           forceConst);
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION AngleConstraintContribs::addContrib
        validate_angle_range(min_angle_deg, AngleRangeBound::Minimum)?;
        validate_angle_range(max_angle_deg, AngleRangeBound::Maximum)?;
        let point_count = owner.positions().len();
        validate_angle_index(idx1, point_count, AngleIndexArgument::First)?;
        validate_angle_index(idx2, point_count, AngleIndexArgument::Second)?;
        validate_angle_index(idx3, point_count, AngleIndexArgument::Third)?;
        if !(max_angle_deg >= min_angle_deg) {
            return Err(ForceFieldKernelError::PackedAngleBoundsOrder);
        }

        self.contribs.push(AngleConstraintContribsParams {
            idx1,
            idx2,
            idx3,
            min_angle: min_angle_deg,
            max_angle: max_angle_deg,
            force_constant,
        });
        Ok(())
    }

    pub(crate) fn add_contrib_relative(
        &mut self,
        owner: &ForceField<'_>,
        idx1: u32,
        idx2: u32,
        idx3: u32,
        relative: bool,
        mut min_angle_deg: f64,
        mut max_angle_deg: f64,
        force_constant: f64,
    ) -> Result<(), ForceFieldKernelError> {
        // BEGIN RDKIT CPP FUNCTION AngleConstraintContribs::addContrib (AngleConstraints.cpp:42-70)
        // RDKit❗✔️: void AngleConstraintContribs::addContrib(unsigned int idx1, unsigned int idx2,
        // RDKit❗✔️:                                          unsigned int idx3, bool relative,
        // RDKit❗✔️:                                          double minAngleDeg, double maxAngleDeg,
        // RDKit❗✔️:                                          double forceConst) {
        // RDKit❗✔️:   const RDGeom::PointPtrVect &pos = dp_forceField->positions();
        // RDKit❗✔️:   URANGE_CHECK(idx1, pos.size());
        // RDKit❗✔️:   URANGE_CHECK(idx2, pos.size());
        // RDKit❗✔️:   URANGE_CHECK(idx3, pos.size());
        // RDKit❗✔️:   PRECONDITION(maxAngleDeg >= minAngleDeg,
        // RDKit❗✔️:                "minAngleDeg must be <= maxAngleDeg");
        // RDKit❗✔️:   if (relative) {
        // RDKit❗✔️:     const RDGeom::Point3D &p1 = *((RDGeom::Point3D *)pos[idx1]);
        // RDKit❗✔️:     const RDGeom::Point3D &p2 = *((RDGeom::Point3D *)pos[idx2]);
        // RDKit❗✔️:     const RDGeom::Point3D &p3 = *((RDGeom::Point3D *)pos[idx3]);
        // RDKit❗✔️:     const RDGeom::Point3D r[2] = {p1 - p2, p3 - p2};
        // RDKit❗✔️:     const double rLengthSq[2] = {std::max(1.0e-5, r[0].lengthSq()),
        // RDKit❗✔️:                                  std::max(1.0e-5, r[1].lengthSq())};
        // RDKit❗✔️:     double cosTheta = r[0].dotProduct(r[1]) / sqrt(rLengthSq[0] * rLengthSq[1]);
        // RDKit❗✔️:     cosTheta = std::clamp(cosTheta, -1.0, 1.0);
        // RDKit❗✔️:     const double angle = RAD2DEG * acos(cosTheta);
        // RDKit❗✔️:     minAngleDeg += angle;
        // RDKit❗✔️:     maxAngleDeg += angle;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   RANGE_CHECK(0.0, minAngleDeg, 180.0);
        // RDKit❗✔️:   RANGE_CHECK(0.0, maxAngleDeg, 180.0);
        // RDKit❗✔️:   d_contribs.emplace_back(idx1, idx2, idx3, minAngleDeg, maxAngleDeg,
        // RDKit❗✔️:                           forceConst);
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION AngleConstraintContribs::addContrib
        let positions = owner.positions();
        let point_count = positions.len();
        validate_angle_index(idx1, point_count, AngleIndexArgument::First)?;
        validate_angle_index(idx2, point_count, AngleIndexArgument::Second)?;
        validate_angle_index(idx3, point_count, AngleIndexArgument::Third)?;
        if !(max_angle_deg >= min_angle_deg) {
            return Err(ForceFieldKernelError::PackedAngleBoundsOrder);
        }
        if relative {
            let angle = angle_degrees(
                &positions[idx1 as usize][..],
                &positions[idx2 as usize][..],
                &positions[idx3 as usize][..],
            );
            min_angle_deg += angle;
            max_angle_deg += angle;
        }
        validate_angle_range(min_angle_deg, AngleRangeBound::Minimum)?;
        validate_angle_range(max_angle_deg, AngleRangeBound::Maximum)?;

        self.contribs.push(AngleConstraintContribsParams {
            idx1,
            idx2,
            idx3,
            min_angle: min_angle_deg,
            max_angle: max_angle_deg,
            force_constant,
        });
        Ok(())
    }

    fn compute_angle_term(&self, angle: f64, contrib: &AngleConstraintContribsParams) -> f64 {
        // BEGIN RDKIT CPP FUNCTION AngleConstraintContribs::computeAngleTerm (AngleConstraints.cpp:71-80)
        // RDKit❗✔️: double AngleConstraintContribs::computeAngleTerm(
        // RDKit❗✔️:     const double &angle, const AngleConstraintContribsParams &contrib) const {
        // RDKit❗✔️:   double angleTerm = 0.0;
        // RDKit❗✔️:   if (angle < contrib.minAngle) {
        // RDKit❗✔️:     angleTerm = angle - contrib.minAngle;
        // RDKit❗✔️:   } else if (angle > contrib.maxAngle) {
        // RDKit❗✔️:     angleTerm = angle - contrib.maxAngle;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return angleTerm;
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION AngleConstraintContribs::computeAngleTerm
        let mut angle_term = 0.0;
        if angle < contrib.min_angle {
            angle_term = angle - contrib.min_angle;
        } else if angle > contrib.max_angle {
            angle_term = angle - contrib.max_angle;
        }
        angle_term
    }

    pub(crate) fn empty(&self) -> bool {
        self.contribs.is_empty()
    }

    pub(crate) fn size(&self) -> usize {
        self.contribs.len()
    }

    pub(crate) fn get_energy(&self, context: &EvaluationContext<'_>) -> f64 {
        // BEGIN RDKIT CPP FUNCTION AngleConstraintContribs::getEnergy (AngleConstraints.cpp:82-103)
        // RDKit❗✔️: double AngleConstraintContribs::getEnergy(double *pos) const {
        // RDKit❗✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit❗✔️:   PRECONDITION(pos, "bad vector");
        // RDKit❗✔️:   double accum = 0.0;
        // RDKit❗✔️:   for (const auto &contrib : d_contribs) {
        // RDKit❗✔️:     const RDGeom::Point3D p1(pos[3 * contrib.idx1], pos[3 * contrib.idx1 + 1],
        // RDKit❗✔️:                              pos[3 * contrib.idx1 + 2]);
        // RDKit❗✔️:     const RDGeom::Point3D p2(pos[3 * contrib.idx2], pos[3 * contrib.idx2 + 1],
        // RDKit❗✔️:                              pos[3 * contrib.idx2 + 2]);
        // RDKit❗✔️:     const RDGeom::Point3D p3(pos[3 * contrib.idx3], pos[3 * contrib.idx3 + 1],
        // RDKit❗✔️:                              pos[3 * contrib.idx3 + 2]);
        // RDKit❗✔️:     const RDGeom::Point3D r[2] = {p1 - p2, p3 - p2};
        // RDKit❗✔️:     const double rLengthSq[2] = {std::max(1.0e-5, r[0].lengthSq()),
        // RDKit❗✔️:                                  std::max(1.0e-5, r[1].lengthSq())};
        // RDKit❗✔️:     double cosTheta = r[0].dotProduct(r[1]) / sqrt(rLengthSq[0] * rLengthSq[1]);
        // RDKit❗✔️:     cosTheta = std::clamp(cosTheta, -1.0, 1.0);
        // RDKit❗✔️:     const double angle = RAD2DEG * acos(cosTheta);
        // RDKit❗✔️:     const double angleTerm = computeAngleTerm(angle, contrib);
        // RDKit❗✔️:     accum += contrib.forceConstant * angleTerm * angleTerm;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return accum;
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION AngleConstraintContribs::getEnergy
        // EvaluationContext is a non-null borrowed owner/coordinate boundary;
        // no owner pointer or cache lookup is retained by the contribution.
        let pos = context.coordinates;
        let mut accum = 0.0;
        for contrib in &self.contribs {
            let start1 = 3 * contrib.idx1 as usize;
            let start2 = 3 * contrib.idx2 as usize;
            let start3 = 3 * contrib.idx3 as usize;
            let p1 = &pos[start1..start1 + 3];
            let p2 = &pos[start2..start2 + 3];
            let p3 = &pos[start3..start3 + 3];
            let angle = angle_degrees(p1, p2, p3);
            let angle_term = self.compute_angle_term(angle, contrib);
            accum += contrib.force_constant * angle_term * angle_term;
        }
        accum
    }

    pub(crate) fn get_grad(
        &self,
        context: &mut EvaluationContext<'_>,
        gradient: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        // BEGIN RDKIT CPP FUNCTION AngleConstraintContribs::getGrad (AngleConstraints.cpp:105-141)
        // RDKit❗✔️: void AngleConstraintContribs::getGrad(double *pos, double *grad) const {
        // RDKit❗✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit❗✔️:   PRECONDITION(pos, "bad vector");
        // RDKit❗✔️:   PRECONDITION(grad, "bad vector");
        // RDKit❗✔️:   for (const auto &contrib : d_contribs) {
        // RDKit❗✔️:     const RDGeom::Point3D p1(pos[3 * contrib.idx1], pos[3 * contrib.idx1 + 1],
        // RDKit❗✔️:                              pos[3 * contrib.idx1 + 2]);
        // RDKit❗✔️:     const RDGeom::Point3D p2(pos[3 * contrib.idx2], pos[3 * contrib.idx2 + 1],
        // RDKit❗✔️:                              pos[3 * contrib.idx2 + 2]);
        // RDKit❗✔️:     const RDGeom::Point3D p3(pos[3 * contrib.idx3], pos[3 * contrib.idx3 + 1],
        // RDKit❗✔️:                              pos[3 * contrib.idx3 + 2]);
        // RDKit❗✔️:     const RDGeom::Point3D r[2] = {p1 - p2, p3 - p2};
        // RDKit❗✔️:     const double rLengthSq[2] = {std::max(1.0e-5, r[0].lengthSq()),
        // RDKit❗✔️:                                  std::max(1.0e-5, r[1].lengthSq())};
        // RDKit❗✔️:     double cosTheta = r[0].dotProduct(r[1]) / sqrt(rLengthSq[0] * rLengthSq[1]);
        // RDKit❗✔️:     cosTheta = std::clamp(cosTheta, -1.0, 1.0);
        // RDKit❗✔️:     const double angle = RAD2DEG * acos(cosTheta);
        // RDKit❗✔️:     const double angleTerm = computeAngleTerm(angle, contrib);
        // RDKit❗✔️:
        // RDKit❗✔️:     const double dE_dTheta = 2.0 * RAD2DEG * contrib.forceConstant * angleTerm;
        // RDKit❗✔️:
        // RDKit❗✔️:     const RDGeom::Point3D rp = r[1].crossProduct(r[0]);
        // RDKit❗✔️:     const double prefactor = dE_dTheta / std::max(1.0e-5, rp.length());
        // RDKit❗✔️:     const double t[2] = {-prefactor / rLengthSq[0], prefactor / rLengthSq[1]};
        // RDKit❗✔️:     RDGeom::Point3D dedp[3];
        // RDKit❗✔️:     dedp[0] = r[0].crossProduct(rp) * t[0];
        // RDKit❗✔️:     dedp[2] = r[1].crossProduct(rp) * t[1];
        // RDKit❗✔️:     dedp[1] = -dedp[0] - dedp[2];
        // RDKit❗✔️:     double *g[3] = {&(grad[3 * contrib.idx1]), &(grad[3 * contrib.idx2]),
        // RDKit❗✔️:                     &(grad[3 * contrib.idx3])};
        // RDKit❗✔️:     for (unsigned int i = 0; i < 3; ++i) {
        // RDKit❗✔️:       g[i][0] += dedp[i].x;
        // RDKit❗✔️:       g[i][1] += dedp[i].y;
        // RDKit❗✔️:       g[i][2] += dedp[i].z;
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION AngleConstraintContribs::getGrad
        // The borrowed context and mutable slice encode the source's non-null
        // owner/position/gradient boundary. Insertion already checked indices;
        // this call uses the supplied trial coordinates without storing owner state.
        let pos = context.coordinates;
        for contrib in &self.contribs {
            let start1 = 3 * contrib.idx1 as usize;
            let start2 = 3 * contrib.idx2 as usize;
            let start3 = 3 * contrib.idx3 as usize;
            let p1 = [pos[start1], pos[start1 + 1], pos[start1 + 2]];
            let p2 = [pos[start2], pos[start2 + 1], pos[start2 + 2]];
            let p3 = [pos[start3], pos[start3 + 1], pos[start3 + 2]];
            let r = [
                angle_vector_difference(&p1, &p2),
                angle_vector_difference(&p3, &p2),
            ];
            let r0_length_sq = angle_vector_length_sq(&r[0]);
            let r1_length_sq = angle_vector_length_sq(&r[1]);
            let r_length_sq = [
                if 1.0e-5 < r0_length_sq {
                    r0_length_sq
                } else {
                    1.0e-5
                },
                if 1.0e-5 < r1_length_sq {
                    r1_length_sq
                } else {
                    1.0e-5
                },
            ];
            let mut cos_theta =
                angle_vector_dot_product(&r[0], &r[1]) / (r_length_sq[0] * r_length_sq[1]).sqrt();
            if cos_theta < -1.0 {
                cos_theta = -1.0;
            } else if cos_theta > 1.0 {
                cos_theta = 1.0;
            }
            let angle = RAD2DEG * cos_theta.acos();
            let angle_term = self.compute_angle_term(angle, contrib);

            let d_e_d_theta = 2.0 * RAD2DEG * contrib.force_constant * angle_term;

            let rp = angle_vector_cross_product(&r[1], &r[0]);
            let rp_length = angle_vector_length(&rp);
            let prefactor = d_e_d_theta
                / if 1.0e-5 < rp_length {
                    rp_length
                } else {
                    1.0e-5
                };
            let t = [-prefactor / r_length_sq[0], prefactor / r_length_sq[1]];
            let dedp0 = angle_vector_scale(&angle_vector_cross_product(&r[0], &rp), t[0]);
            let dedp2 = angle_vector_scale(&angle_vector_cross_product(&r[1], &rp), t[1]);
            let negative_dedp0 = angle_vector_negated(&dedp0);
            let dedp1 = angle_vector_difference(&negative_dedp0, &dedp2);
            let dedp = [dedp0, dedp1, dedp2];
            let atom_indices = [contrib.idx1, contrib.idx2, contrib.idx3];
            for i in 0..3 {
                let gradient_start = 3 * atom_indices[i] as usize;
                gradient[gradient_start] += dedp[i][0];
                gradient[gradient_start + 1] += dedp[i][1];
                gradient[gradient_start + 2] += dedp[i][2];
            }
        }
        Ok(())
    }
}

impl ForceFieldContribution for AngleConstraintContribs {
    fn get_energy(
        &self,
        context: &mut EvaluationContext<'_>,
    ) -> Result<f64, ForceFieldKernelError> {
        Ok(AngleConstraintContribs::get_energy(self, context))
    }

    fn get_grad(
        &self,
        context: &mut EvaluationContext<'_>,
        gradient: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        AngleConstraintContribs::get_grad(self, context, gradient)
    }

    fn copy(&self) -> Box<dyn ForceFieldContribution> {
        // BEGIN RDKIT CPP FUNCTION AngleConstraintContribs::copy (AngleConstraints.h:90-92)
        // RDKit❗✔️: AngleConstraintContribs *copy() const override {
        // RDKit❗✔️:   return new AngleConstraintContribs(*this);
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION AngleConstraintContribs::copy
        // Cloning preserves the ordered parameter vector; the evaluation owner
        // is supplied later through a borrowed context rather than stored.
        Box::new(self.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::{
        AngleConstraintContribs, AngleConstraintContribsParams, AngleIndexArgument,
        AngleRangeBound, ForceField, ForceFieldKernelError,
    };

    fn cf3d_f20_force_field<'a>(positions: &'a mut [Vec<f64>]) -> ForceField<'a> {
        let mut force_field = ForceField::default();
        force_field
            .positions_mut()
            .extend(positions.iter_mut().map(Vec::as_mut_slice));
        force_field.initialize().unwrap();
        force_field
    }

    fn cf3d_f20_energy(
        force_field: &mut ForceField<'_>,
        contribution: &AngleConstraintContribs,
        coordinates: &[f64],
    ) -> f64 {
        let context = force_field.evaluation_context(coordinates);
        contribution.get_energy(&context)
    }

    #[test]
    fn cf3d_f20_empty_collection_has_source_defaults_and_zero_energy() {
        // AngleConstraints.h:19-31 initializes indices and bounds to zero,
        // forceConstant to one; AngleConstraints.cpp:23-26 starts empty.
        let defaults = AngleConstraintContribsParams::default();
        assert_eq!(
            defaults,
            AngleConstraintContribsParams {
                idx1: 0,
                idx2: 0,
                idx3: 0,
                min_angle: 0.0,
                max_angle: 0.0,
                force_constant: 1.0,
            }
        );

        let mut positions = vec![vec![0.0; 3]; 3];
        let mut force_field = cf3d_f20_force_field(&mut positions);
        let contribution = AngleConstraintContribs::new(&force_field);

        assert!(contribution.empty());
        assert_eq!(contribution.size(), 0);
        assert_eq!(
            cf3d_f20_energy(
                &mut force_field,
                &contribution,
                &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
            ),
            0.0
        );
    }

    #[test]
    fn cf3d_f20_absolute_insertion_preserves_check_order_and_typed_errors() {
        // AngleConstraints.cpp:28-39 checks both ranges before idx1/2/3,
        // then the exact maxAngleDeg >= minAngleDeg precondition.
        let mut positions = vec![vec![0.0; 3]; 3];
        let force_field = cf3d_f20_force_field(&mut positions);
        let mut contribution = AngleConstraintContribs::new(&force_field);

        assert_eq!(
            contribution.add_contrib(&force_field, 3, 1, 2, -1.0, 180.0, 1.0),
            Err(ForceFieldKernelError::AngleOutOfRange {
                bound: AngleRangeBound::Minimum
            })
        );
        assert_eq!(
            contribution.add_contrib(&force_field, 3, 1, 2, 0.0, 181.0, 1.0),
            Err(ForceFieldKernelError::AngleOutOfRange {
                bound: AngleRangeBound::Maximum
            })
        );
        assert_eq!(
            contribution.add_contrib(&force_field, 3, 3, 3, 120.0, 100.0, 1.0),
            Err(ForceFieldKernelError::AngleIndexOutOfRange {
                argument: AngleIndexArgument::First,
                index: 3,
                upper_bound: 3,
            })
        );
        assert_eq!(
            contribution.add_contrib(&force_field, 0, 3, 3, 0.0, 180.0, 1.0),
            Err(ForceFieldKernelError::AngleIndexOutOfRange {
                argument: AngleIndexArgument::Second,
                index: 3,
                upper_bound: 3,
            })
        );
        assert_eq!(
            contribution.add_contrib(&force_field, 0, 1, 3, 0.0, 180.0, 1.0),
            Err(ForceFieldKernelError::AngleIndexOutOfRange {
                argument: AngleIndexArgument::Third,
                index: 3,
                upper_bound: 3,
            })
        );
        assert_eq!(
            contribution.add_contrib(&force_field, 3, 1, 2, f64::NAN, 120.0, 1.0),
            Err(ForceFieldKernelError::AngleIndexOutOfRange {
                argument: AngleIndexArgument::First,
                index: 3,
                upper_bound: 3,
            })
        );
        assert_eq!(
            contribution.add_contrib(&force_field, 0, 1, 2, f64::NAN, 120.0, 1.0),
            Err(ForceFieldKernelError::PackedAngleBoundsOrder)
        );
        let order_error = ForceFieldKernelError::PackedAngleBoundsOrder;
        assert_eq!(order_error.source_category(), "Pre-condition Violation");
        assert_eq!(
            order_error.source_expression().as_deref(),
            Some("maxAngleDeg >= minAngleDeg")
        );
        assert!(contribution.empty());

        // RANGE_CHECK includes both endpoints, and forceConst is unchecked.
        contribution
            .add_contrib(&force_field, 0, 1, 2, 0.0, 180.0, f64::NAN)
            .unwrap();
        assert_eq!(contribution.size(), 1);
    }

    #[test]
    fn cf3d_f20_relative_insertion_offsets_only_true_mode_from_owner_geometry() {
        // AngleConstraints.cpp:42-70 measures owner positions only in the
        // relative branch, adds the angle to each bound, and still range-checks
        // both bounds when relative is false.
        let mut positions = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
        ];
        let mut force_field = cf3d_f20_force_field(&mut positions);
        let mut contribution = AngleConstraintContribs::new(&force_field);
        contribution
            .add_contrib_relative(&force_field, 0, 1, 2, true, 10.0, 20.0, 1.0)
            .unwrap();
        contribution
            .add_contrib_relative(&force_field, 0, 1, 2, false, 10.0, 20.0, 2.0)
            .unwrap();

        assert_eq!(contribution.size(), 2);
        assert_eq!(contribution.contribs[0].min_angle, 100.0);
        assert_eq!(contribution.contribs[0].max_angle, 110.0);
        assert_eq!(contribution.contribs[1].min_angle, 10.0);
        assert_eq!(contribution.contribs[1].max_angle, 20.0);

        // Evaluation uses the passed trial coordinates, distinct from the
        // 90-degree owner geometry used to prepare the first bounds.
        assert_eq!(
            cf3d_f20_energy(
                &mut force_field,
                &contribution,
                &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0]
            ),
            10_200.0
        );
    }

    #[test]
    fn cf3d_f20_relative_insertion_preserves_error_precedence_and_nan_bounds() {
        // The relative overload checks indices, then bounds order, then
        // computes the angle and validates final minimum before maximum.
        let mut positions = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
        ];
        let force_field = cf3d_f20_force_field(&mut positions);
        let mut contribution = AngleConstraintContribs::new(&force_field);

        assert_eq!(
            contribution.add_contrib_relative(&force_field, 3, 1, 2, true, 20.0, 10.0, 1.0),
            Err(ForceFieldKernelError::AngleIndexOutOfRange {
                argument: AngleIndexArgument::First,
                index: 3,
                upper_bound: 3,
            })
        );
        assert_eq!(
            contribution.add_contrib_relative(&force_field, 0, 1, 2, true, 20.0, 10.0, 1.0),
            Err(ForceFieldKernelError::PackedAngleBoundsOrder)
        );
        assert_eq!(
            contribution.add_contrib_relative(&force_field, 0, 1, 2, false, -1.0, 20.0, 1.0),
            Err(ForceFieldKernelError::AngleOutOfRange {
                bound: AngleRangeBound::Minimum
            })
        );
        assert_eq!(
            contribution.add_contrib_relative(&force_field, 0, 1, 2, true, 100.0, 120.0, 1.0),
            Err(ForceFieldKernelError::AngleOutOfRange {
                bound: AngleRangeBound::Minimum
            })
        );
        assert_eq!(
            contribution.add_contrib_relative(&force_field, 0, 1, 2, true, 50.0, 100.0, 1.0),
            Err(ForceFieldKernelError::AngleOutOfRange {
                bound: AngleRangeBound::Maximum
            })
        );
        assert!(contribution.empty());

        let mut nan_positions = vec![
            vec![f64::NAN, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
        ];
        let nan_force_field = cf3d_f20_force_field(&mut nan_positions);
        let mut nan_contribution = AngleConstraintContribs::new(&nan_force_field);
        nan_contribution
            .add_contrib_relative(&nan_force_field, 0, 1, 2, true, 10.0, 20.0, 1.0)
            .unwrap();
        assert_eq!(nan_contribution.size(), 1);
        assert!(nan_contribution.contribs[0].min_angle.is_nan());
        assert!(nan_contribution.contribs[0].max_angle.is_nan());
    }

    #[test]
    fn cf3d_f20_energy_keeps_strict_flat_bottom_branches_and_shared_atoms() {
        // AngleConstraints.cpp:82-103 and computeAngleTerm:71-80. Two
        // distinct triples share atom 1; equality and inside angles add zero.
        let mut positions = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![-1.0, 0.0, 0.0],
            vec![0.0, -1.0, 0.0],
        ];
        let mut force_field = cf3d_f20_force_field(&mut positions);
        let mut contribution = AngleConstraintContribs::new(&force_field);
        contribution
            .add_contrib(&force_field, 0, 1, 2, 100.0, 120.0, 1.0)
            .unwrap();
        contribution
            .add_contrib(&force_field, 3, 1, 4, 60.0, 80.0, 2.0)
            .unwrap();
        contribution
            .add_contrib(&force_field, 0, 1, 2, 90.0, 100.0, 3.0)
            .unwrap();
        contribution
            .add_contrib(&force_field, 3, 1, 4, 80.0, 90.0, 4.0)
            .unwrap();
        contribution
            .add_contrib(&force_field, 0, 1, 2, 80.0, 100.0, 5.0)
            .unwrap();
        contribution
            .add_contrib(&force_field, 3, 1, 4, 90.0, 90.0, 6.0)
            .unwrap();
        assert_eq!(contribution.size(), 6);

        let coordinates = [
            1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0, 0.0,
        ];
        // The below/above terms are 1*10^2 and 2*10^2; endpoints,
        // interior values, and equal bounds each produce an exact zero term.
        assert_eq!(
            cf3d_f20_energy(&mut force_field, &contribution, &coordinates),
            300.0
        );
    }

    #[test]
    fn cf3d_f20_energy_accumulates_packed_terms_in_insertion_order() {
        // The source performs a scalar accum += per vector element; it does
        // not reorder or deduplicate repeated triples.
        let mut positions = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
        ];
        let mut force_field = cf3d_f20_force_field(&mut positions);
        let mut contribution = AngleConstraintContribs::new(&force_field);
        for force_constant in [1.0e20, 1.0e3, -1.0e20] {
            contribution
                .add_contrib(&force_field, 0, 1, 2, 91.0, 100.0, force_constant)
                .unwrap();
        }
        assert_eq!(contribution.size(), 3);
        assert_eq!(
            cf3d_f20_energy(
                &mut force_field,
                &contribution,
                &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
            ),
            0.0
        );
    }

    fn cf3d_f21_assert_gradient(actual: &[f64; 9], expected: [f64; 9]) {
        for index in 0..9 {
            assert!(
                (actual[index] - expected[index]).abs() <= 1.0e-9,
                "component {index}: actual={}, expected={}",
                actual[index],
                expected[index]
            );
        }
    }

    #[test]
    fn cf3d_f21_dispatch_accumulates_overlapping_terms_from_fixed_formula() {
        // AngleConstraints.cpp:105-141 visits every packed term in insertion
        // order and adds each point derivative to the caller's existing buffer.
        const MAGNITUDE: f64 = 1145.9155902616465;
        let mut positions = vec![vec![0.0; 3]; 3];
        let mut force_field = cf3d_f20_force_field(&mut positions);
        let mut contribution = AngleConstraintContribs::new(&force_field);

        for (minimum, maximum, force_constant) in [
            (100.0, 120.0, 1.0), // lower violation: +MAGNITUDE at idx1.y
            (110.0, 130.0, 2.0), // overlapping lower violation: +4*MAGNITUDE
            (60.0, 80.0, 2.0),   // upper violation: -2*MAGNITUDE
            (90.0, 90.0, 3.0),   // equal bounds at the angle
            (90.0, 110.0, 3.0),  // equality at the lower boundary
            (70.0, 90.0, 4.0),   // equality at the upper boundary
            (80.0, 100.0, 5.0),  // strict interior
        ] {
            contribution
                .add_contrib(&force_field, 0, 1, 2, minimum, maximum, force_constant)
                .unwrap();
        }
        assert_eq!(contribution.size(), 7);
        force_field.contributions.push(Box::new(contribution));

        let coordinates = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let mut gradient = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        force_field.calc_grad(&coordinates, &mut gradient).unwrap();

        // At a right angle, the three active source terms sum to 3 times
        // 2*RAD2DEG*10. The equality and interior terms add exact zero.
        let total = 3.0 * MAGNITUDE;
        cf3d_f21_assert_gradient(
            &gradient,
            [
                1.0,
                2.0 + total,
                3.0,
                4.0 - total,
                5.0 - total,
                6.0,
                7.0 + total,
                8.0,
                9.0,
            ],
        );
    }

    #[test]
    fn cf3d_f21_degenerate_geometry_uses_source_floors_and_keeps_additive_buffer() {
        // AngleConstraints.cpp:117-128 floors both arm lengths and rp.length()
        // at 1e-5. Coincident, one-zero-arm, and collinear cases remain finite.
        let mut positions = vec![vec![0.0; 3]; 3];
        let mut force_field = cf3d_f20_force_field(&mut positions);
        let mut contribution = AngleConstraintContribs::new(&force_field);
        contribution
            .add_contrib(&force_field, 0, 1, 2, 10.0, 20.0, 1.0)
            .unwrap();
        force_field.contributions.push(Box::new(contribution));

        for coordinates in [
            [0.0; 9],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0],
        ] {
            let original = [1.0, -2.0, 3.0, -4.0, 5.0, -6.0, 7.0, -8.0, 9.0];
            let mut gradient = original;
            force_field.calc_grad(&coordinates, &mut gradient).unwrap();
            assert_eq!(gradient, original);
            assert!(gradient.iter().all(|value| value.is_finite()));
        }
    }
}
