use super::{
    EvaluationContext, ForceField, ForceFieldContribution, ForceFieldIndexArgument,
    ForceFieldKernelError,
};
use crate::geometry::Point3;

#[derive(Clone, Debug, PartialEq)]
pub(super) struct PositionConstraintContrib {
    at_idx: u32,
    max_displ: f64,
    pos0: Point3,
    force_constant: f64,
}

impl PositionConstraintContrib {
    pub(super) fn new(
        owner: &ForceField<'_>,
        idx: u32,
        max_displ: f64,
        force_constant: f64,
    ) -> Result<Self, ForceFieldKernelError> {
        // BEGIN RDKIT CPP FUNCTION ForceFields::PositionConstraintContrib::PositionConstraintContrib (PositionConstraint.cpp:16-31)
        // RDKit✔️✔️: PositionConstraintContrib::PositionConstraintContrib(ForceField *owner,
        // RDKit✔️✔️:                                                      unsigned int idx,
        // RDKit✔️✔️:                                                      double maxDispl,
        // RDKit✔️✔️:                                                      double forceConst) {
        // RDKit✔️✔️:   PRECONDITION(owner, "bad owner");
        // RDKit✔️✔️:   const RDGeom::PointPtrVect &pos = owner->positions();
        // RDKit✔️✔️:   URANGE_CHECK(idx, pos.size());
        // RDKit✔️✔️:
        // RDKit✔️✔️:   dp_forceField = owner;
        // RDKit✔️✔️:   d_atIdx = idx;
        // RDKit✔️✔️:   d_maxDispl = maxDispl;
        // RDKit✔️✔️:   d_pos0 = *((RDGeom::Point3D *)pos[idx]);
        // RDKit✔️✔️:   d_forceConstant = forceConst;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION ForceFields::PositionConstraintContrib::PositionConstraintContrib
        // `owner` is a non-null borrow. Keep the source range check before
        // reading and copying the selected 3D reference point; the source owner
        // pointer is not retained because callbacks receive the live context.
        let positions = owner.positions();
        let point_count = positions.len() as u32;
        if idx >= point_count {
            return Err(ForceFieldKernelError::IndexOutOfRange {
                argument: ForceFieldIndexArgument::I,
                index: idx,
                upper_bound: point_count,
            });
        }
        let point = &positions[idx as usize][..];
        let pos0 = Point3 {
            x: point[0],
            y: point[1],
            z: point[2],
        };

        Ok(Self {
            at_idx: idx,
            max_displ,
            pos0,
            force_constant,
        })
    }

    fn energy(&self, coordinates: &[f64]) -> f64 {
        // BEGIN RDKIT CPP FUNCTION ForceFields::PositionConstraintContrib::getEnergy (PositionConstraint.cpp:34-48)
        // RDKit✔️✔️: double PositionConstraintContrib::getEnergy(double *pos) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:
        // RDKit✔️✔️:   RDGeom::Point3D p(pos[3 * d_atIdx], pos[3 * d_atIdx + 1],
        // RDKit✔️✔️:                     pos[3 * d_atIdx + 2]);
        // RDKit✔️✔️:   double dist = (p - d_pos0).length();
        // RDKit✔️✔️:   double distTerm = std::max(dist - d_maxDispl, 0.0);
        // RDKit✔️✔️:   double res = 0.5 * d_forceConstant * distTerm * distTerm;
        // RDKit✔️✔️:
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION ForceFields::PositionConstraintContrib::getEnergy
        // A borrowed slice is non-null. Spell std::max's first-operand rule
        // directly: when `distTerm < 0.0` is false (including NaN), it returns
        // `distTerm`, unlike Rust's NaN-selecting `f64::max`.
        let offset = 3 * self.at_idx as usize;
        let point = Point3 {
            x: coordinates[offset],
            y: coordinates[offset + 1],
            z: coordinates[offset + 2],
        };
        let dist = Point3::difference(&point, &self.pos0).length();
        let dist_term = dist - self.max_displ;
        let dist_term = if dist_term < 0.0 { 0.0 } else { dist_term };
        let result = 0.5 * self.force_constant * dist_term * dist_term;

        result
    }

    fn gradient(&self, coordinates: &[f64], gradient: &mut [f64]) {
        // BEGIN RDKIT CPP FUNCTION ForceFields::PositionConstraintContrib::getGrad (PositionConstraint.cpp:50-66)
        // RDKit✔️✔️: void PositionConstraintContrib::getGrad(double *pos, double *grad) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:   PRECONDITION(grad, "bad vector");
        // RDKit✔️✔️:
        // RDKit✔️✔️:   RDGeom::Point3D p(pos[3 * d_atIdx], pos[3 * d_atIdx + 1],
        // RDKit✔️✔️:                     pos[3 * d_atIdx + 2]);
        // RDKit✔️✔️:   double dist = (p - d_pos0).length();
        // RDKit✔️✔️:
        // RDKit✔️✔️:   double preFactor = 0.0;
        // RDKit✔️✔️:   if (dist > d_maxDispl) {
        // RDKit✔️✔️:     preFactor = dist - d_maxDispl;
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     return;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   preFactor *= d_forceConstant;
        // RDKit✔️✔️:
        // RDKit✔️✔️:   for (unsigned int i = 0; i < 3; ++i) {
        // RDKit✔️✔️:     double dGrad = preFactor * (p[i] - d_pos0[i]) / std::max(dist, 1.0e-8);
        // RDKit✔️✔️:     grad[3 * d_atIdx + i] += dGrad;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION ForceFields::PositionConstraintContrib::getGrad
        // The borrowed coordinate and gradient slices provide the source's
        // non-null buffers. Preserve the strict flat-bottom early return and
        // keep each source max/divide/add inside its three-component loop.
        let offset = 3 * self.at_idx as usize;
        let point = Point3 {
            x: coordinates[offset],
            y: coordinates[offset + 1],
            z: coordinates[offset + 2],
        };
        let dist = Point3::difference(&point, &self.pos0).length();

        let mut pre_factor = 0.0;
        if dist > self.max_displ {
            pre_factor = dist - self.max_displ;
        } else {
            return;
        }
        pre_factor *= self.force_constant;

        for i in 0..3 {
            let (point_component, reference_component) = match i {
                0 => (point.x, self.pos0.x),
                1 => (point.y, self.pos0.y),
                _ => (point.z, self.pos0.z),
            };
            let denominator = if dist < 1.0e-8 { 1.0e-8 } else { dist };
            let gradient_term = pre_factor * (point_component - reference_component) / denominator;
            gradient[offset + i] += gradient_term;
        }
    }
}

impl ForceFieldContribution for PositionConstraintContrib {
    fn get_energy(
        &self,
        context: &mut EvaluationContext<'_>,
    ) -> Result<f64, ForceFieldKernelError> {
        Ok(self.energy(context.coordinates))
    }

    fn get_grad(
        &self,
        context: &mut EvaluationContext<'_>,
        gradient: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        self.gradient(context.coordinates, gradient);
        Ok(())
    }

    fn copy(&self) -> Box<dyn ForceFieldContribution> {
        Box::new(self.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::super::{ForceField, ForceFieldIndexArgument, ForceFieldKernelError};
    use super::PositionConstraintContrib;

    fn one_point_force_field<'a>(position: &'a mut Vec<f64>) -> ForceField<'a> {
        let mut force_field = ForceField::new(3);
        force_field.positions_mut().push(position);
        force_field.initialize().unwrap();
        force_field
    }

    fn assert_source_close(actual: f64, expected: f64) {
        let tolerance = expected.abs() * 1.0e-14 + f64::MIN_POSITIVE;
        assert!(
            (actual - expected).abs() <= tolerance,
            "expected {expected:?}, got {actual:?}"
        );
    }

    #[test]
    fn cf3d_f24_energy_flat_bottom_matches_zero_inside_threshold_and_exceeded() {
        // RDKit source: ForceField/PositionConstraint.cpp:34-48.
        // The expected values follow distTerm=max(dist-maxDispl,0) and
        // res=0.5*forceConst*distTerm*distTerm for this fixed Point3D input.
        let mut reference = vec![1.0, 2.0, 3.0];
        let mut force_field = one_point_force_field(&mut reference);
        let contribution = PositionConstraintContrib::new(&force_field, 0, 0.5, 10.0).unwrap();
        force_field.contributions.push(Box::new(contribution));

        assert_eq!(force_field.calc_energy(&[1.0, 2.0, 3.0]), Ok(0.0));
        assert_eq!(force_field.calc_energy(&[1.25, 2.0, 3.0]), Ok(0.0));
        assert_eq!(force_field.calc_energy(&[1.5, 2.0, 3.0]), Ok(0.0));
        assert_eq!(force_field.calc_energy(&[2.0, 2.0, 3.0]), Ok(1.25));
    }

    #[test]
    fn cf3d_f24_source_max_keeps_nan_and_unvalidated_scalar_values() {
        // RDKit source: ForceField/PositionConstraint.cpp:34-48,50-66.
        // C++ std::max keeps its first operand when that operand compares
        // false against zero, so NaN remains NaN; source has no scalar guards.
        let mut reference = vec![0.0, 0.0, 0.0];
        let force_field = one_point_force_field(&mut reference);

        let negative_limit = PositionConstraintContrib::new(&force_field, 0, -0.1, 10.0).unwrap();
        assert_source_close(negative_limit.energy(&[0.0, 0.0, 0.0]), 0.05);

        let nan_limit = PositionConstraintContrib::new(&force_field, 0, f64::NAN, 10.0).unwrap();
        assert!(nan_limit.energy(&[0.0, 0.0, 0.0]).is_nan());
        let mut nan_limit_gradient = [1.0, -2.0, 3.0];
        nan_limit.gradient(&[1.0, 0.0, 0.0], &mut nan_limit_gradient);
        assert_eq!(nan_limit_gradient, [1.0, -2.0, 3.0]);

        let negative_force = PositionConstraintContrib::new(&force_field, 0, 0.5, -10.0).unwrap();
        assert_eq!(negative_force.energy(&[1.0, 0.0, 0.0]), -1.25);
        let mut negative_force_gradient = [0.0; 3];
        negative_force.gradient(&[1.0, 0.0, 0.0], &mut negative_force_gradient);
        assert_eq!(negative_force_gradient, [-5.0, 0.0, 0.0]);
    }

    #[test]
    fn cf3d_f24_gradient_is_strict_additive_and_uses_source_denominator_floor() {
        // RDKit source: ForceField/PositionConstraint.cpp:50-66.
        // The active x gradient is (1.0-0.5)*10*(2.0-1.0)/1.0 = 5.0.
        let mut reference = vec![1.0, 2.0, 3.0];
        let mut force_field = one_point_force_field(&mut reference);
        let contribution = PositionConstraintContrib::new(&force_field, 0, 0.5, 10.0).unwrap();
        force_field.contributions.push(Box::new(contribution));

        let mut active = [1.0, -1.0, 2.0];
        force_field
            .calc_grad(&[2.0, 2.0, 3.0], &mut active)
            .unwrap();
        assert_eq!(active, [6.0, -1.0, 2.0]);

        let mut inside = [4.0, 5.0, 6.0];
        force_field
            .calc_grad(&[1.25, 2.0, 3.0], &mut inside)
            .unwrap();
        assert_eq!(inside, [4.0, 5.0, 6.0]);
        let mut threshold = [7.0, 8.0, 9.0];
        force_field
            .calc_grad(&[1.5, 2.0, 3.0], &mut threshold)
            .unwrap();
        assert_eq!(threshold, [7.0, 8.0, 9.0]);

        // The source denominator uses max(dist, 1e-8), even while the
        // source-permitted negative maximum displacement activates at this
        // shorter-than-floor distance.
        let mut zero_reference = vec![0.0, 0.0, 0.0];
        let small_force_field = one_point_force_field(&mut zero_reference);
        let floor = PositionConstraintContrib::new(&small_force_field, 0, -1.0e-8, 2.0).unwrap();
        let mut floor_gradient = [0.0; 3];
        floor.gradient(&[5.0e-9, 0.0, 0.0], &mut floor_gradient);
        assert!((floor_gradient[0] - 1.5e-8).abs() <= 1.0e-22);
        assert_eq!(floor_gradient[1..], [0.0, 0.0]);
    }

    #[test]
    fn cf3d_f24_captures_reference_when_live_positions_move_and_checks_index() {
        // RDKit source: ForceField/PositionConstraint.cpp:16-31,34-48.
        // Constructor copies the current point by value; moving the owner's
        // coordinate row later does not replace d_pos0.
        let mut reference = vec![1.0, 2.0, 3.0];
        let mut force_field = one_point_force_field(&mut reference);
        let contribution = PositionConstraintContrib::new(&force_field, 0, 0.3, 10.0).unwrap();
        force_field.contributions.push(Box::new(contribution));
        force_field.positions_mut()[0][0] = 10.0;

        assert_source_close(force_field.calc_energy(&[1.5, 2.0, 3.0]).unwrap(), 0.2);
        assert_source_close(force_field.calc_energy_current(None).unwrap(), 378.45);

        assert_eq!(
            PositionConstraintContrib::new(&force_field, 1, 0.3, 10.0),
            Err(ForceFieldKernelError::IndexOutOfRange {
                argument: ForceFieldIndexArgument::I,
                index: 1,
                upper_bound: 1,
            })
        );
    }
}
