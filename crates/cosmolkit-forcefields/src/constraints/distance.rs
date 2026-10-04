use super::ForceFieldKernelError;
use super::{EvaluationContext, ForceField, ForceFieldContribution, ForceFieldIndexArgument};

#[derive(Clone, Debug)]
pub(super) struct DistanceConstraintContrib {
    end1_idx: u32,
    end2_idx: u32,
    min_len: f64,
    max_len: f64,
    force_constant: f64,
}

impl DistanceConstraintContrib {
    pub(super) fn new_absolute(
        owner: &ForceField<'_>,
        idx1: u32,
        idx2: u32,
        min_len: f64,
        max_len: f64,
        force_constant: f64,
    ) -> Result<Self, ForceFieldKernelError> {
        // RDKit source: ForceField/DistanceConstraint.cpp:16-30
        // RDKit❗✔️: DistanceConstraintContrib::DistanceConstraintContrib(
        // RDKit❗✔️:     ForceField *owner, unsigned int idx1, unsigned int idx2, double minLen,
        // RDKit❗✔️:     double maxLen, double forceConst) {
        // RDKit❗✔️:   PRECONDITION(owner, "bad owner");
        // RDKit❗✔️:   URANGE_CHECK(idx1, owner->positions().size());
        // RDKit❗✔️:   URANGE_CHECK(idx2, owner->positions().size());
        // RDKit❗✔️:   PRECONDITION(maxLen >= minLen, "bad bounds");
        // RDKit❗✔️:   dp_forceField = owner;
        // RDKit❗✔️:   d_end1Idx = idx1;
        // RDKit❗✔️:   d_end2Idx = idx2;
        // RDKit❗✔️:   d_minLen = minLen;
        // RDKit❗✔️:   d_maxLen = maxLen;
        // RDKit❗✔️:   d_forceConstant = forceConst;
        // RDKit❗✔️: }
        // `owner: &ForceField` represents the non-null owner precondition. Keep
        // the source endpoint and bounds check order; a failed source precondition
        // becomes its corresponding private typed error.
        let point_count = owner.positions().len() as u32;
        if idx1 >= point_count {
            return Err(ForceFieldKernelError::IndexOutOfRange {
                argument: ForceFieldIndexArgument::I,
                index: idx1,
                upper_bound: point_count,
            });
        }
        if idx2 >= point_count {
            return Err(ForceFieldKernelError::IndexOutOfRange {
                argument: ForceFieldIndexArgument::J,
                index: idx2,
                upper_bound: point_count,
            });
        }
        if !(max_len >= min_len) {
            return Err(ForceFieldKernelError::BadBounds);
        }

        Ok(Self {
            end1_idx: idx1,
            end2_idx: idx2,
            min_len,
            max_len,
            force_constant,
        })
    }

    pub(super) fn new_relative(
        owner: &ForceField<'_>,
        idx1: u32,
        idx2: u32,
        relative: bool,
        mut min_len: f64,
        mut max_len: f64,
        force_constant: f64,
    ) -> Result<Self, ForceFieldKernelError> {
        // RDKit source: ForceField/DistanceConstraint.cpp:32-53
        // RDKit❗✔️: DistanceConstraintContrib::DistanceConstraintContrib(
        // RDKit❗✔️:     ForceField *owner, unsigned int idx1, unsigned int idx2, bool relative,
        // RDKit❗✔️:     double minLen, double maxLen, double forceConst) {
        // RDKit❗✔️:   PRECONDITION(owner, "bad owner");
        // RDKit❗✔️:   const RDGeom::PointPtrVect &pos = owner->positions();
        // RDKit❗✔️:   URANGE_CHECK(idx1, pos.size());
        // RDKit❗✔️:   URANGE_CHECK(idx2, pos.size());
        // RDKit❗✔️:   PRECONDITION(maxLen >= minLen, "bad bounds");
        // RDKit❗✔️:   if (relative) {
        // RDKit❗✔️:     RDGeom::Point3D &p1 = *((RDGeom::Point3D *)pos[idx1]);
        // RDKit❗✔️:     RDGeom::Point3D &p2 = *((RDGeom::Point3D *)pos[idx2]);
        // RDKit❗✔️:     const auto dist = (p1 - p2).length();
        // RDKit❗✔️:     minLen = std::max(dist + minLen, 0.0);
        // RDKit❗✔️:     maxLen = std::max(dist + maxLen, 0.0);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   d_minLen = minLen;
        // RDKit❗✔️:   d_maxLen = maxLen;
        // RDKit❗✔️:   dp_forceField = owner;
        // RDKit❗✔️:   d_end1Idx = idx1;
        // RDKit❗✔️:   d_end2Idx = idx2;
        // RDKit❗✔️:   d_forceConstant = forceConst;
        // RDKit❗✔️: }
        // RDKit source: Geometry/point.cpp:64-70, Geometry/point.h:158-161
        // RDKit❗✔️: Point3D operator-(const Point3D &p1, const Point3D &p2) {
        // RDKit❗✔️:   Point3D res;
        // RDKit❗✔️:   res.x = p1.x - p2.x;
        // RDKit❗✔️:   res.y = p1.y - p2.y;
        // RDKit❗✔️:   res.z = p1.z - p2.z;
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // RDKit❗✔️: double length() const override {
        // RDKit❗✔️:   double res = x * x + y * y + z * z;
        // RDKit❗✔️:   return sqrt(res);
        // RDKit❗✔️: }
        let positions = owner.positions();
        let point_count = positions.len() as u32;
        if idx1 >= point_count {
            return Err(ForceFieldKernelError::IndexOutOfRange {
                argument: ForceFieldIndexArgument::I,
                index: idx1,
                upper_bound: point_count,
            });
        }
        if idx2 >= point_count {
            return Err(ForceFieldKernelError::IndexOutOfRange {
                argument: ForceFieldIndexArgument::J,
                index: idx2,
                upper_bound: point_count,
            });
        }
        if !(max_len >= min_len) {
            return Err(ForceFieldKernelError::BadBounds);
        }
        if relative {
            let p1 = &positions[idx1 as usize];
            let p2 = &positions[idx2 as usize];
            let dx = p1[0] - p2[0];
            let dy = p1[1] - p2[1];
            let dz = p1[2] - p2[2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();

            // C++ std::max(a, b) returns a when !(a < b), including NaN and
            // signed-zero cases. Keep the source operand order explicitly.
            let min_with_distance = dist + min_len;
            min_len = if min_with_distance < 0.0 {
                0.0
            } else {
                min_with_distance
            };
            let max_with_distance = dist + max_len;
            max_len = if max_with_distance < 0.0 {
                0.0
            } else {
                max_with_distance
            };
        }

        Ok(Self {
            end1_idx: idx1,
            end2_idx: idx2,
            min_len,
            max_len,
            force_constant,
        })
    }

    pub(super) fn min_len(&self) -> f64 {
        self.min_len
    }

    pub(super) fn max_len(&self) -> f64 {
        self.max_len
    }

    pub(super) fn get_energy(
        &self,
        context: &mut EvaluationContext<'_>,
    ) -> Result<f64, ForceFieldKernelError> {
        // RDKit source: ForceField/DistanceConstraint.cpp:55-69
        // RDKit❗✔️: double DistanceConstraintContrib::getEnergy(double *pos) const {
        // RDKit❗✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit❗✔️:   PRECONDITION(pos, "bad vector");
        // RDKit❗✔️:   double dist = dp_forceField->distance(d_end1Idx, d_end2Idx, pos);
        // RDKit❗✔️:   double distTerm = 0.0;
        // RDKit❗✔️:   if (dist < d_minLen) {
        // RDKit❗✔️:     distTerm = d_minLen - dist;
        // RDKit❗✔️:   } else if (dist > d_maxLen) {
        // RDKit❗✔️:     distTerm = dist - d_maxLen;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   double res = 0.5 * d_forceConstant * distTerm * distTerm;
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // RDKit source: ForceField/ForceField.cpp:172-203
        // RDKit❗✔️: double ForceField::distance(unsigned int i, unsigned int j, double *pos) {
        // RDKit❗✔️:   PRECONDITION(df_init, "not initialized");
        // RDKit❗✔️:   URANGE_CHECK(i, d_numPoints);
        // RDKit❗✔️:   URANGE_CHECK(j, d_numPoints);
        // RDKit❗✔️:   if (j < i) {
        // RDKit❗✔️:     int tmp = j;
        // RDKit❗✔️:     j = i;
        // RDKit❗✔️:     i = tmp;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   unsigned int idx = i + j * (j + 1) / 2;
        // RDKit❗✔️:   CHECK_INVARIANT(idx < d_matSize, "Bad index");
        // RDKit❗✔️:   double &res = dp_distMat[idx];
        // RDKit❗✔️:   if (res < 0.0) {
        // RDKit❗✔️:     // we need to calculate this distance:
        // RDKit❗✔️:     if (!pos) {
        // RDKit❗✔️:       res = 0.0;
        // RDKit❗✔️:       for (unsigned int idx = 0; idx < d_dimension; ++idx) {
        // RDKit❗✔️:         double tmp =
        // RDKit❗✔️:             (*(this->positions()[i]))[idx] - (*(this->positions()[j]))[idx];
        // RDKit❗✔️:         res += tmp * tmp;
        // RDKit❗✔️:       }
        // RDKit❗✔️:     } else {
        // RDKit❗✔️:       res = 0.0;
        // RDKit❗✔️:       double *pi = &(pos[d_dimension * i]), *pj = &(pos[d_dimension * j]);
        // RDKit❗✔️:       for (unsigned int idx = 0; idx < d_dimension; ++idx, ++pi, ++pj) {
        // RDKit❗✔️:         double tmp = *pi - *pj;
        // RDKit❗✔️:         res += tmp * tmp;
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:     res = sqrt(res);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // The borrowed context is the non-null owner and position-vector
        // boundary. `distance` preserves initialized/range/cache checks and
        // uses the current trial coordinates on a cache miss. It allocates no
        // per-term storage; a hit is O(1), a miss is O(dimension).
        let dist = context.distance(self.end1_idx, self.end2_idx)?;
        let mut dist_term = 0.0;
        if dist < self.min_len {
            dist_term = self.min_len - dist;
        } else if dist > self.max_len {
            dist_term = dist - self.max_len;
        }
        Ok(0.5 * self.force_constant * dist_term * dist_term)
    }

    pub(super) fn get_grad(
        &self,
        context: &mut EvaluationContext<'_>,
        gradient: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        // RDKit source: ForceField/DistanceConstraint.cpp:71-95
        // RDKit❗✔️: void DistanceConstraintContrib::getGrad(double *pos, double *grad) const {
        // RDKit❗✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit❗✔️:   PRECONDITION(pos, "bad vector");
        // RDKit❗✔️:   PRECONDITION(grad, "bad vector");
        // RDKit❗✔️:
        // RDKit❗✔️:   double dist = dp_forceField->distance(d_end1Idx, d_end2Idx, pos);
        // RDKit❗✔️:
        // RDKit❗✔️:   double preFactor = 0.0;
        // RDKit❗✔️:   if (dist < d_minLen) {
        // RDKit❗✔️:     preFactor = dist - d_minLen;
        // RDKit❗✔️:   } else if (dist > d_maxLen) {
        // RDKit❗✔️:     preFactor = dist - d_maxLen;
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   preFactor *= d_forceConstant;
        // RDKit❗✔️:
        // RDKit❗✔️:   double *end1Coords = &(pos[3 * d_end1Idx]);
        // RDKit❗✔️:   double *end2Coords = &(pos[3 * d_end2Idx]);
        // RDKit❗✔️:   for (unsigned int i = 0; i < 3; ++i) {
        // RDKit❗✔️:     double dGrad =
        // RDKit❗✔️:         preFactor * (end1Coords[i] - end2Coords[i]) / std::max(dist, 1.0e-8);
        // RDKit❗✔️:     grad[3 * d_end1Idx + i] += dGrad;
        // RDKit❗✔️:     grad[3 * d_end2Idx + i] -= dGrad;
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // RDKit source: ForceField/ForceField.cpp:172-203 (distance cache used above)
        // RDKit❗✔️: double ForceField::distance(unsigned int i, unsigned int j, double *pos) {
        // RDKit❗✔️:   PRECONDITION(df_init, "not initialized");
        // RDKit❗✔️:   URANGE_CHECK(i, d_numPoints);
        // RDKit❗✔️:   URANGE_CHECK(j, d_numPoints);
        // RDKit❗✔️:   if (j < i) {
        // RDKit❗✔️:     int tmp = j;
        // RDKit❗✔️:     j = i;
        // RDKit❗✔️:     i = tmp;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   unsigned int idx = i + j * (j + 1) / 2;
        // RDKit❗✔️:   CHECK_INVARIANT(idx < d_matSize, "Bad index");
        // RDKit❗✔️:   double &res = dp_distMat[idx];
        // RDKit❗✔️:   if (res < 0.0) {
        // RDKit❗✔️:     // we need to calculate this distance:
        // RDKit❗✔️:     if (!pos) {
        // RDKit❗✔️:       res = 0.0;
        // RDKit❗✔️:       for (unsigned int idx = 0; idx < d_dimension; ++idx) {
        // RDKit❗✔️:         double tmp =
        // RDKit❗✔️:             (*(this->positions()[i]))[idx] - (*(this->positions()[j]))[idx];
        // RDKit❗✔️:         res += tmp * tmp;
        // RDKit❗✔️:       }
        // RDKit❗✔️:     } else {
        // RDKit❗✔️:       res = 0.0;
        // RDKit❗✔️:       double *pi = &(pos[d_dimension * i]), *pj = &(pos[d_dimension * j]);
        // RDKit❗✔️:       for (unsigned int idx = 0; idx < d_dimension; ++idx, ++pi, ++pj) {
        // RDKit❗✔️:         double tmp = *pi - *pj;
        // RDKit❗✔️:         res += tmp * tmp;
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:     res = sqrt(res);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // The borrowed owner and slices encode the non-null preconditions.
        // Use the same triangular cache as ForceField::distance; calcGrad
        // deliberately does not reset it. Preserve additive endpoint order
        // without a temporary gradient buffer.
        let dist = context.distance(self.end1_idx, self.end2_idx)?;
        let mut pre_factor = 0.0;
        if dist < self.min_len {
            pre_factor = dist - self.min_len;
        } else if dist > self.max_len {
            pre_factor = dist - self.max_len;
        } else {
            return Ok(());
        }
        pre_factor *= self.force_constant;

        let end1 = 3 * self.end1_idx as usize;
        let end2 = 3 * self.end2_idx as usize;
        for component in 0..3 {
            let denominator = if dist < 1.0e-8 { 1.0e-8 } else { dist };
            let d_grad = pre_factor
                * (context.coordinates[end1 + component] - context.coordinates[end2 + component])
                / denominator;
            gradient[end1 + component] += d_grad;
            gradient[end2 + component] -= d_grad;
        }
        Ok(())
    }
}

impl ForceFieldContribution for DistanceConstraintContrib {
    fn get_energy(
        &self,
        context: &mut EvaluationContext<'_>,
    ) -> Result<f64, ForceFieldKernelError> {
        DistanceConstraintContrib::get_energy(self, context)
    }

    fn get_grad(
        &self,
        context: &mut EvaluationContext<'_>,
        gradient: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        DistanceConstraintContrib::get_grad(self, context, gradient)
    }

    fn copy(&self) -> Box<dyn ForceFieldContribution> {
        Box::new(self.clone())
    }
}
