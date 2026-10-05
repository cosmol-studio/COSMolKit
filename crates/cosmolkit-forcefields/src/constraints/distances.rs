use super::{
    EvaluationContext, ForceField, ForceFieldContribution, ForceFieldIndexArgument,
    ForceFieldKernelError,
};

#[derive(Clone, Copy, Debug, PartialEq)]
struct DistanceConstraintContribsParams {
    // RDKit source: ForceField/DistanceConstraints.h:19-31
    // RDKit❗✔️: unsigned int idx1{0};       //!< index of atom1 of the distance constraint
    // RDKit❗✔️: unsigned int idx2{0};       //!< index of atom2 of the distance constraint
    // RDKit❗✔️: double minLen{0.0};         //!< lower bound of the flat bottom potential
    // RDKit❗✔️: double maxLen{0.0};         //!< upper bound of the flat bottom potential
    // RDKit❗✔️: double forceConstant{1.0};  //!< force constant for distance constraint
    idx1: u32,
    idx2: u32,
    min_len: f64,
    max_len: f64,
    force_constant: f64,
}

impl Default for DistanceConstraintContribsParams {
    fn default() -> Self {
        Self {
            idx1: 0,
            idx2: 0,
            min_len: 0.0,
            max_len: 0.0,
            force_constant: 1.0,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct DistanceConstraintContribs {
    // RDKit source: ForceField/DistanceConstraints.h:91
    // RDKit❗✔️: std::vector<DistanceConstraintContribsParams> d_contribs;
    contribs: Vec<DistanceConstraintContribsParams>,
}

impl DistanceConstraintContribs {
    pub(crate) fn new(_owner: &ForceField<'_>) -> Self {
        // RDKit source: ForceField/DistanceConstraints.cpp:17-20
        // RDKit❗✔️: DistanceConstraintContribs::DistanceConstraintContribs(ForceField *owner) {
        // RDKit❗✔️:   PRECONDITION(owner, "bad owner");
        // RDKit❗✔️:   dp_forceField = owner;
        // RDKit❗✔️: }
        // The borrow encodes the non-null construction boundary. Rust does
        // not retain a self-owner pointer; insertion borrows the current
        // owner and evaluation borrows the shared context.
        Self::default()
    }

    pub(crate) fn add_contrib(
        &mut self,
        owner: &ForceField<'_>,
        idx1: u32,
        idx2: u32,
        min_len: f64,
        max_len: f64,
        force_constant: f64,
    ) -> Result<(), ForceFieldKernelError> {
        // RDKit source: ForceField/DistanceConstraints.cpp:22-31
        // RDKit❗✔️: void DistanceConstraintContribs::addContrib(unsigned int idx1,
        // RDKit❗✔️:                                             unsigned int idx2, double minLen,
        // RDKit❗✔️:                                             double maxLen,
        // RDKit❗✔️:                                             double forceConstant) {
        // RDKit❗✔️:   URANGE_CHECK(idx1, dp_forceField->positions().size());
        // RDKit❗✔️:   URANGE_CHECK(idx2, dp_forceField->positions().size());
        // RDKit❗✔️:   PRECONDITION(maxLen >= minLen, "bad bounds");
        // RDKit❗✔️:   d_contribs.emplace_back(idx1, idx2, minLen, maxLen, forceConstant);
        // RDKit❗✔️: }
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
        self.contribs.push(DistanceConstraintContribsParams {
            idx1,
            idx2,
            min_len,
            max_len,
            force_constant,
        });
        Ok(())
    }

    pub(crate) fn add_contrib_relative(
        &mut self,
        owner: &ForceField<'_>,
        idx1: u32,
        idx2: u32,
        relative: bool,
        mut min_len: f64,
        mut max_len: f64,
        force_constant: f64,
    ) -> Result<(), ForceFieldKernelError> {
        // RDKit source: ForceField/DistanceConstraints.cpp:33-52
        // RDKit❗✔️: void DistanceConstraintContribs::addContrib(unsigned int idx1,
        // RDKit❗✔️:                                             unsigned int idx2, bool relative,
        // RDKit❗✔️:                                             double minLen, double maxLen,
        // RDKit❗✔️:                                             double forceConstant) {
        // RDKit❗✔️:   const RDGeom::PointPtrVect &pos = dp_forceField->positions();
        // RDKit❗✔️:   URANGE_CHECK(idx1, pos.size());
        // RDKit❗✔️:   URANGE_CHECK(idx2, pos.size());
        // RDKit❗✔️:   PRECONDITION(maxLen >= minLen, "bad bounds");
        // RDKit❗✔️:   if (relative) {
        // RDKit❗✔️:     const RDGeom::Point3D p1 = *((RDGeom::Point3D *)pos[idx1]);
        // RDKit❗✔️:     const RDGeom::Point3D p2 = *((RDGeom::Point3D *)pos[idx2]);
        // RDKit❗✔️:     const auto distance = (p1 - p2).length();
        // RDKit❗✔️:     minLen = std::max(minLen + distance, 0.0);
        // RDKit❗✔️:     maxLen = std::max(maxLen + distance, 0.0);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   d_contribs.emplace_back(idx1, idx2, minLen, maxLen, forceConstant);
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
            let distance = (dx * dx + dy * dy + dz * dz).sqrt();

            let min_with_distance = min_len + distance;
            min_len = if min_with_distance < 0.0 {
                0.0
            } else {
                min_with_distance
            };
            let max_with_distance = max_len + distance;
            max_len = if max_with_distance < 0.0 {
                0.0
            } else {
                max_with_distance
            };
        }
        self.contribs.push(DistanceConstraintContribsParams {
            idx1,
            idx2,
            min_len,
            max_len,
            force_constant,
        });
        Ok(())
    }

    pub(crate) fn empty(&self) -> bool {
        // RDKit✔️✔️: bool empty() const { return d_contribs.empty(); }
        self.contribs.is_empty()
    }

    pub(crate) fn size(&self) -> usize {
        // RDKit✔️✔️: unsigned int size() const { return d_contribs.size(); }
        self.contribs.len()
    }

    pub(crate) fn get_energy(
        &self,
        context: &EvaluationContext<'_>,
    ) -> Result<f64, ForceFieldKernelError> {
        // RDKit source: ForceField/DistanceConstraints.cpp:54-72
        // RDKit❗✔️: double DistanceConstraintContribs::getEnergy(double *pos) const {
        // RDKit❗✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit❗✔️:   PRECONDITION(pos, "bad vector");
        // RDKit❗✔️:
        // RDKit❗✔️:   double accum = 0.0;
        // RDKit❗✔️:   for (const auto &contrib : d_contribs) {
        // RDKit❗✔️:     const auto distance2 =
        // RDKit❗✔️:         dp_forceField->distance2(contrib.idx1, contrib.idx2, pos);
        // RDKit❗✔️:     double difference = 0.0;
        // RDKit❗✔️:     if (distance2 < contrib.minLen * contrib.minLen) {
        // RDKit❗✔️:       difference = contrib.minLen - std::sqrt(distance2);
        // RDKit❗✔️:     } else if (distance2 > contrib.maxLen * contrib.maxLen) {
        // RDKit❗✔️:       difference = std::sqrt(distance2) - contrib.maxLen;
        // RDKit❗✔️:     } else {
        // RDKit❗✔️:       continue;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     accum += 0.5 * contrib.forceConstant * difference * difference;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return accum;
        // RDKit❗✔️: }
        // RDKit source: ForceField/ForceField.cpp:206-232
        // RDKit❗✔️: double ForceField::distance2(unsigned int i, unsigned int j,
        // RDKit❗✔️:                              double *pos) const {
        // RDKit❗✔️:   PRECONDITION(df_init, "not initialized");
        // RDKit❗✔️:   URANGE_CHECK(i, d_numPoints);
        // RDKit❗✔️:   URANGE_CHECK(j, d_numPoints);
        // RDKit❗✔️:   if (j < i) {
        // RDKit❗✔️:     int tmp = j;
        // RDKit❗✔️:     j = i;
        // RDKit❗✔️:     i = tmp;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   double res;
        // RDKit❗✔️:   if (!pos) {
        // RDKit❗✔️:     res = 0.0;
        // RDKit❗✔️:     for (unsigned int idx = 0; idx < d_dimension; ++idx) {
        // RDKit❗✔️:       double tmp =
        // RDKit❗✔️:           (*(this->positions()[i]))[idx] - (*(this->positions()[j]))[idx];
        // RDKit❗✔️:       res += tmp * tmp;
        // RDKit❗✔️:     }
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     res = 0.0;
        // RDKit❗✔️:     double *pi = &(pos[d_dimension * i]), *pj = &(pos[d_dimension * j]);
        // RDKit❗✔️:     for (unsigned int idx = 0; idx < d_dimension; ++idx, ++pi, ++pj) {
        // RDKit❗✔️:       double tmp = *pi - *pj;
        // RDKit❗✔️:       res += tmp * tmp;
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // EvaluationContext::distance2 preserves these checks and ordered
        // arithmetic without touching ForceField's cached distance matrix.
        // One ordered pass uses O(1) scratch and performs no per-term allocation.
        let mut accum = 0.0;
        for contrib in &self.contribs {
            let distance2 = context.distance2(contrib.idx1, contrib.idx2)?;
            let mut difference = 0.0;
            if distance2 < contrib.min_len * contrib.min_len {
                difference = contrib.min_len - distance2.sqrt();
            } else if distance2 > contrib.max_len * contrib.max_len {
                difference = distance2.sqrt() - contrib.max_len;
            } else {
                continue;
            }
            accum += 0.5 * contrib.force_constant * difference * difference;
        }
        Ok(accum)
    }

    pub(crate) fn get_grad(
        &self,
        context: &EvaluationContext<'_>,
        gradient: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        // RDKit source: ForceField/DistanceConstraints.cpp:71-96
        // RDKit❗✔️: void DistanceConstraintContribs::getGrad(double *pos, double *grad) const {
        // RDKit❗✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit❗✔️:   PRECONDITION(pos, "bad vector");
        // RDKit❗✔️:   PRECONDITION(grad, "bad vector");
        // RDKit❗✔️:
        // RDKit❗✔️:   for (const auto &contrib : d_contribs) {
        // RDKit❗✔️:     double preFactor = 0.0;
        // RDKit❗✔️:     double distance = 0.0;
        // RDKit❗✔️:     const auto distance2 =
        // RDKit❗✔️:         dp_forceField->distance2(contrib.idx1, contrib.idx2, pos);
        // RDKit❗✔️:     if (distance2 < contrib.minLen * contrib.minLen) {
        // RDKit❗✔️:       distance = std::sqrt(distance2);
        // RDKit❗✔️:       preFactor = distance - contrib.minLen;
        // RDKit❗✔️:     } else if (distance2 > contrib.maxLen * contrib.maxLen) {
        // RDKit❗✔️:       distance = std::sqrt(distance2);
        // RDKit❗✔️:       preFactor = distance - contrib.maxLen;
        // RDKit❗✔️:     } else {
        // RDKit❗✔️:       continue;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     preFactor *= contrib.forceConstant;
        // RDKit❗✔️:     preFactor /= std::max(1.0e-8, distance);
        // RDKit❗✔️:     const double *atom1Coords = &(pos[3 * contrib.idx1]);
        // RDKit❗✔️:     const double *atom2Coords = &(pos[3 * contrib.idx2]);
        // RDKit❗✔️:     for (unsigned int i = 0; i < 3; i++) {
        // RDKit❗✔️:       const double dGrad = preFactor * (atom1Coords[i] - atom2Coords[i]);
        // RDKit❗✔️:       grad[3 * contrib.idx1 + i] += dGrad;
        // RDKit❗✔️:       grad[3 * contrib.idx2 + i] -= dGrad;
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // RDKit source: ForceField/ForceField.cpp:206-232
        // RDKit❗✔️: double ForceField::distance2(unsigned int i, unsigned int j,
        // RDKit❗✔️:                              double *pos) const {
        // RDKit❗✔️:   PRECONDITION(df_init, "not initialized");
        // RDKit❗✔️:   URANGE_CHECK(i, d_numPoints);
        // RDKit❗✔️:   URANGE_CHECK(j, d_numPoints);
        // RDKit❗✔️:   if (j < i) {
        // RDKit❗✔️:     int tmp = j;
        // RDKit❗✔️:     j = i;
        // RDKit❗✔️:     i = tmp;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   double res;
        // RDKit❗✔️:   if (!pos) {
        // RDKit❗✔️:     res = 0.0;
        // RDKit❗✔️:     for (unsigned int idx = 0; idx < d_dimension; ++idx) {
        // RDKit❗✔️:       double tmp =
        // RDKit❗✔️:           (*(this->positions()[i]))[idx] - (*(this->positions()[j]))[idx];
        // RDKit❗✔️:       res += tmp * tmp;
        // RDKit❗✔️:     }
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     res = 0.0;
        // RDKit❗✔️:     double *pi = &(pos[d_dimension * i]), *pj = &(pos[d_dimension * j]);
        // RDKit❗✔️:     for (unsigned int idx = 0; idx < d_dimension; ++idx, ++pi, ++pj) {
        // RDKit❗✔️:       double tmp = *pi - *pj;
        // RDKit❗✔️:       res += tmp * tmp;
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // The short-lived context provides the source owner/position reads;
        // it returns typed initialization/index errors and never fills the
        // triangular cache. No per-term allocation is added to the source
        // ordered loop.
        for contrib in &self.contribs {
            let mut pre_factor = 0.0;
            let mut distance = 0.0;
            let distance2 = context.distance2(contrib.idx1, contrib.idx2)?;
            if distance2 < contrib.min_len * contrib.min_len {
                distance = distance2.sqrt();
                pre_factor = distance - contrib.min_len;
            } else if distance2 > contrib.max_len * contrib.max_len {
                distance = distance2.sqrt();
                pre_factor = distance - contrib.max_len;
            } else {
                continue;
            }
            pre_factor *= contrib.force_constant;
            pre_factor /= 1.0e-8_f64.max(distance);
            let atom1 = contrib.idx1.wrapping_mul(3) as usize;
            let atom2 = contrib.idx2.wrapping_mul(3) as usize;
            for i in 0..3 {
                let d_grad =
                    pre_factor * (context.coordinates[atom1 + i] - context.coordinates[atom2 + i]);
                gradient[atom1 + i] += d_grad;
                gradient[atom2 + i] -= d_grad;
            }
        }
        Ok(())
    }
}

impl ForceFieldContribution for DistanceConstraintContribs {
    fn get_energy(
        &self,
        context: &mut EvaluationContext<'_>,
    ) -> Result<f64, ForceFieldKernelError> {
        DistanceConstraintContribs::get_energy(self, context)
    }

    fn get_grad(
        &self,
        context: &mut EvaluationContext<'_>,
        gradient: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        DistanceConstraintContribs::get_grad(self, context, gradient)
    }

    fn copy(&self) -> Box<dyn ForceFieldContribution> {
        Box::new(self.clone())
    }
}
