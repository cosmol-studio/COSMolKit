// Copyright (C) 2013 Paolo Tosco and other RDKit contributors.
//
// This file contains behavior ported from RDKit and remains subject to its
// BSD license, included in the pinned source tree.

use crate::geometry::Point3;
use crate::kernel::{EvaluationContext, ForceFieldContribution, ForceFieldKernelError};

use super::params::{clip_to_one, is_double_zero};
use super::utils::{calc_inversion_coefficients, calculate_cos_y};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum InversionIndexArgument {
    First,
    Second,
    Third,
    Fourth,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum InversionContributionError {
    IndexOutOfRange {
        argument: InversionIndexArgument,
        index: u32,
        upper_bound: usize,
    },
}

impl InversionContributionError {
    pub(crate) const fn source_category(self) -> &'static str {
        match self {
            Self::IndexOutOfRange { .. } => "Range Error",
        }
    }

    pub(crate) const fn source_message(self) -> &'static str {
        match self {
            Self::IndexOutOfRange { argument, .. } => match argument {
                InversionIndexArgument::First => "idx1",
                InversionIndexArgument::Second => "idx2",
                InversionIndexArgument::Third => "idx3",
                InversionIndexArgument::Fourth => "idx4",
            },
        }
    }

    pub(crate) const fn range_detail(self) -> (u32, usize) {
        match self {
            Self::IndexOutOfRange {
                index, upper_bound, ..
            } => (index, upper_bound),
        }
    }
}

impl std::fmt::Display for InversionContributionError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Debug::fmt(self, formatter)
    }
}

impl std::error::Error for InversionContributionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        None
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct InversionContrib {
    at1_idx: u32,
    at2_idx: u32,
    at3_idx: u32,
    at4_idx: u32,
    force_constant: f64,
    c0: f64,
    c1: f64,
    c2: f64,
}

impl InversionContrib {
    pub(crate) fn new(
        positions: &[&mut [f64]],
        idx1: u32,
        idx2: u32,
        idx3: u32,
        idx4: u32,
        at2_atomic_num: i32,
        is_c_bound_to_o: bool,
    ) -> Result<Self, InversionContributionError> {
        // RDKit's Inversion.h declares oobForceScalingFactor = 1.0. This
        // private Rust entry point preserves that constructor default.
        Self::new_with_scale(
            positions,
            idx1,
            idx2,
            idx3,
            idx4,
            at2_atomic_num,
            is_c_bound_to_o,
            1.0,
        )
    }

    fn new_with_scale(
        positions: &[&mut [f64]],
        idx1: u32,
        idx2: u32,
        idx3: u32,
        idx4: u32,
        at2_atomic_num: i32,
        is_c_bound_to_o: bool,
        oob_force_scaling_factor: f64,
    ) -> Result<Self, InversionContributionError> {
        // BEGIN RDKIT CPP FUNCTION ForceFields::UFF::InversionContrib::InversionContrib (ForceField/UFF/Inversion.cpp:22-45)
        // RDKit❗✔️: InversionContrib::InversionContrib(ForceField *owner, unsigned int idx1,
        // RDKit❗✔️:                                    unsigned int idx2, unsigned int idx3,
        // RDKit❗✔️:                                    unsigned int idx4, int at2AtomicNum,
        // RDKit❗✔️:                                    bool isCBoundToO,
        // RDKit❗✔️:                                    double oobForceScalingFactor) {
        // RDKit❗✔️:   PRECONDITION(owner, "bad owner");
        // RDKit❗✔️:   URANGE_CHECK(idx1, owner->positions().size());
        // RDKit❗✔️:   URANGE_CHECK(idx2, owner->positions().size());
        // RDKit❗✔️:   URANGE_CHECK(idx3, owner->positions().size());
        // RDKit❗✔️:   URANGE_CHECK(idx4, owner->positions().size());
        // RDKit❗✔️:   dp_forceField = owner;
        // RDKit❗✔️:   d_at1Idx = idx1;
        // RDKit❗✔️:   d_at2Idx = idx2;
        // RDKit❗✔️:   d_at3Idx = idx3;
        // RDKit❗✔️:   d_at4Idx = idx4;
        // RDKit❗✔️:   auto invCoeffForceCon = Utils::calcInversionCoefficientsAndForceConstant(
        // RDKit❗✔️:       at2AtomicNum, isCBoundToO);
        // BEGIN RDKIT CPP HELPER ForceFields::UFF::Utils::calcInversionCoefficientsAndForceConstant (ForceField/UFF/Utils.cpp:42-85)
        // RDKit❗✔️: std::tuple<double, double, double, double>
        // RDKit❗✔️: calcInversionCoefficientsAndForceConstant(int at2AtomicNum, bool isCBoundToO) {
        // RDKit❗✔️:   double res = 0.0;
        // RDKit❗✔️:   double C0 = 0.0;
        // RDKit❗✔️:   double C1 = 0.0;
        // RDKit❗✔️:   double C2 = 0.0;
        // RDKit❗✔️:   // if the central atom is sp2 carbon, nitrogen or oxygen
        // RDKit❗✔️:   if ((at2AtomicNum == 6) || (at2AtomicNum == 7) || (at2AtomicNum == 8)) {
        // RDKit❗✔️:     C0 = 1.0;
        // RDKit❗✔️:     C1 = -1.0;
        // RDKit❗✔️:     C2 = 0.0;
        // RDKit❗✔️:     res = (isCBoundToO ? 50.0 : 6.0);
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     // group 5 elements are not clearly explained in the UFF paper
        // RDKit❗✔️:     // the following code was inspired by MCCCS Towhee's ffuff.F
        // RDKit❗✔️:     double w0 = M_PI / 180.0;
        // RDKit❗✔️:     switch (at2AtomicNum) {
        // RDKit❗✔️:       // if the central atom is phosphorous
        // RDKit❗✔️:       case 15:
        // RDKit❗✔️:         w0 *= 84.4339;
        // RDKit❗✔️:         break;
        // RDKit❗✔️:       // if the central atom is arsenic
        // RDKit❗✔️:       case 33:
        // RDKit❗✔️:         w0 *= 86.9735;
        // RDKit❗✔️:         break;
        // RDKit❗✔️:       // if the central atom is antimonium
        // RDKit❗✔️:       case 51:
        // RDKit❗✔️:         w0 *= 87.7047;
        // RDKit❗✔️:         break;
        // RDKit❗✔️:       // if the central atom is bismuth
        // RDKit❗✔️:       case 83:
        // RDKit❗✔️:         w0 *= 90.0;
        // RDKit❗✔️:         break;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     C2 = 1.0;
        // RDKit❗✔️:     C1 = -4.0 * cos(w0);
        // RDKit❗✔️:     C0 = -(C1 * cos(w0) + C2 * cos(2.0 * w0));
        // RDKit❗✔️:     res = 22.0 / (C0 + C1 + C2);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   res /= 3.0;
        // RDKit❗✔️:   return std::make_tuple(res, C0, C1, C2);
        // RDKit❗✔️: }
        // END RDKIT CPP HELPER ForceFields::UFF::Utils::calcInversionCoefficientsAndForceConstant
        // RDKit❗✔️:   d_forceConstant = oobForceScalingFactor * std::get<0>(invCoeffForceCon);
        // RDKit❗✔️:   d_C0 = std::get<1>(invCoeffForceCon);
        // RDKit❗✔️:   d_C1 = std::get<2>(invCoeffForceCon);
        // RDKit❗✔️:   d_C2 = std::get<3>(invCoeffForceCon);
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION ForceFields::UFF::InversionContrib::InversionContrib
        let mut contribution = Self::new_unscaled(
            positions,
            idx1,
            idx2,
            idx3,
            idx4,
            at2_atomic_num,
            is_c_bound_to_o,
        )?;
        contribution.force_constant = oob_force_scaling_factor * contribution.force_constant;
        Ok(contribution)
    }

    pub(crate) fn new_packed_with_scale(
        positions: &[&mut [f64]],
        idx1: u32,
        idx2: u32,
        idx3: u32,
        idx4: u32,
        at2_atomic_num: i32,
        is_c_bound_to_o: bool,
        oob_force_scaling_factor: f64,
    ) -> Result<Self, InversionContributionError> {
        let mut contribution = Self::new_unscaled(
            positions,
            idx1,
            idx2,
            idx3,
            idx4,
            at2_atomic_num,
            is_c_bound_to_o,
        )?;

        // RDKit❗✔️: std::get<0>(invCoeffForceCon) * oobForceScalingFactor
        // Keep this source operation in packed insertion order: force first,
        // scale second, unlike the single-term source constructor.
        contribution.force_constant = contribution.force_constant * oob_force_scaling_factor;
        Ok(contribution)
    }

    fn new_unscaled(
        positions: &[&mut [f64]],
        idx1: u32,
        idx2: u32,
        idx3: u32,
        idx4: u32,
        at2_atomic_num: i32,
        is_c_bound_to_o: bool,
    ) -> Result<Self, InversionContributionError> {
        // BEGIN RDKIT CPP FUNCTION ForceFields::UFF::InversionContrib::InversionContrib (ForceField/UFF/Inversion.cpp:22-45)
        // RDKit❗✔️: InversionContrib::InversionContrib(ForceField *owner, unsigned int idx1,
        // RDKit❗✔️:                                    unsigned int idx2, unsigned int idx3,
        // RDKit❗✔️:                                    unsigned int idx4, int at2AtomicNum,
        // RDKit❗✔️:                                    bool isCBoundToO,
        // RDKit❗✔️:                                    double oobForceScalingFactor) {
        // RDKit❗✔️:   PRECONDITION(owner, "bad owner");
        // A borrowed positions slice replaces the owner and is not retained.
        // RDKit❗✔️:   URANGE_CHECK(idx1, owner->positions().size());
        Self::check_index(idx1, InversionIndexArgument::First, positions.len())?;
        // RDKit❗✔️:   URANGE_CHECK(idx2, owner->positions().size());
        Self::check_index(idx2, InversionIndexArgument::Second, positions.len())?;
        // RDKit❗✔️:   URANGE_CHECK(idx3, owner->positions().size());
        Self::check_index(idx3, InversionIndexArgument::Third, positions.len())?;
        // RDKit❗✔️:   URANGE_CHECK(idx4, owner->positions().size());
        Self::check_index(idx4, InversionIndexArgument::Fourth, positions.len())?;
        // RDKit❗✔️:   d_at1Idx = idx1;
        // RDKit❗✔️:   d_at2Idx = idx2;
        // RDKit❗✔️:   d_at3Idx = idx3;
        // RDKit❗✔️:   d_at4Idx = idx4;

        // BEGIN RDKIT CPP HELPER ForceFields::UFF::Utils::calcInversionCoefficientsAndForceConstant (ForceField/UFF/Utils.cpp:42-85)
        // RDKit❗✔️: std::tuple<double, double, double, double>
        // RDKit❗✔️: calcInversionCoefficientsAndForceConstant(int at2AtomicNum, bool isCBoundToO) {
        // RDKit❗✔️:   double res = 0.0;
        // RDKit❗✔️:   double C0 = 0.0;
        // RDKit❗✔️:   double C1 = 0.0;
        // RDKit❗✔️:   double C2 = 0.0;
        // RDKit❗✔️:   // if the central atom is sp2 carbon, nitrogen or oxygen
        // RDKit❗✔️:   if ((at2AtomicNum == 6) || (at2AtomicNum == 7) || (at2AtomicNum == 8)) {
        // RDKit❗✔️:     C0 = 1.0;
        // RDKit❗✔️:     C1 = -1.0;
        // RDKit❗✔️:     C2 = 0.0;
        // RDKit❗✔️:     res = (isCBoundToO ? 50.0 : 6.0);
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     // group 5 elements are not clearly explained in the UFF paper
        // RDKit❗✔️:     // the following code was inspired by MCCCS Towhee's ffuff.F
        // RDKit❗✔️:     double w0 = M_PI / 180.0;
        // RDKit❗✔️:     switch (at2AtomicNum) {
        // RDKit❗✔️:       // if the central atom is phosphorous
        // RDKit❗✔️:       case 15:
        // RDKit❗✔️:         w0 *= 84.4339;
        // RDKit❗✔️:         break;
        // RDKit❗✔️:       // if the central atom is arsenic
        // RDKit❗✔️:       case 33:
        // RDKit❗✔️:         w0 *= 86.9735;
        // RDKit❗✔️:         break;
        // RDKit❗✔️:       // if the central atom is antimonium
        // RDKit❗✔️:       case 51:
        // RDKit❗✔️:         w0 *= 87.7047;
        // RDKit❗✔️:         break;
        // RDKit❗✔️:       // if the central atom is bismuth
        // RDKit❗✔️:       case 83:
        // RDKit❗✔️:         w0 *= 90.0;
        // RDKit❗✔️:         break;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     C2 = 1.0;
        // RDKit❗✔️:     C1 = -4.0 * cos(w0);
        // RDKit❗✔️:     C0 = -(C1 * cos(w0) + C2 * cos(2.0 * w0));
        // RDKit❗✔️:     res = 22.0 / (C0 + C1 + C2);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   res /= 3.0;
        // RDKit❗✔️:   return std::make_tuple(res, C0, C1, C2);
        // RDKit❗✔️: }
        // END RDKIT CPP HELPER ForceFields::UFF::Utils::calcInversionCoefficientsAndForceConstant
        let (force_constant, c0, c1, c2) =
            calc_inversion_coefficients(at2_atomic_num, is_c_bound_to_o);

        // RDKit❗✔️:   d_C0 = std::get<1>(invCoeffForceCon);
        // RDKit❗✔️:   d_C1 = std::get<2>(invCoeffForceCon);
        // RDKit❗✔️:   d_C2 = std::get<3>(invCoeffForceCon);
        // RDKit❗✔️: }
        // This common private value starts with the helper's unscaled force
        // constant; each source constructor applies its own product order.
        // END RDKIT CPP FUNCTION ForceFields::UFF::InversionContrib::InversionContrib
        Ok(Self {
            at1_idx: idx1,
            at2_idx: idx2,
            at3_idx: idx3,
            at4_idx: idx4,
            force_constant,
            c0,
            c1,
            c2,
        })
    }

    fn check_index(
        index: u32,
        argument: InversionIndexArgument,
        upper_bound: usize,
    ) -> Result<(), InversionContributionError> {
        // BEGIN RDKIT CPP HELPER URANGE_CHECK (RDGeneral/Invariant.h:141-151)
        // RDKit❗✔️: #define URANGE_CHECK(x, hi) \
        // RDKit❗✔️:   if (x >= (hi)) { \
        // RDKit❗✔️:     std::stringstream errstr; \
        // RDKit❗✔️:     errstr << x << " < " << hi; \
        // RDKit❗✔️:     Invar::Invariant inv("Range Error", #x, errstr.str().c_str(), __FILE__, \
        // RDKit❗✔️:                          __LINE__); \
        // RDKit❗✔️:     BOOST_LOG(rdErrorLog) << "\n\n****\n" << inv << "****\n" << std::endl; \
        // RDKit❗✔️:     throw inv; \
        // RDKit❗✔️:   }
        // END RDKIT CPP HELPER URANGE_CHECK
        if index as usize >= upper_bound {
            return Err(InversionContributionError::IndexOutOfRange {
                argument,
                index,
                upper_bound,
            });
        }
        Ok(())
    }

    pub(crate) fn get_energy(&self, context: &mut EvaluationContext<'_>) -> f64 {
        // BEGIN RDKIT CPP FUNCTION ForceFields::UFF::InversionContrib::getEnergy (ForceField/UFF/Inversion.cpp:47-69)
        // RDKit❗✔️: double InversionContrib::getEnergy(double *pos) const {
        // RDKit❗✔️:   PRECONDITION(dp_forceField, "no owner");
        // The borrowed contribution state removes the owner-pointer check;
        // EvaluationContext supplies a nonnull coordinate slice.
        // RDKit❗✔️:   PRECONDITION(pos, "bad vector");
        // RDKit❗✔️:   RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
        // RDKit❗✔️:                      pos[3 * d_at1Idx + 2]);
        // RDKit❗✔️:   RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
        // RDKit❗✔️:                      pos[3 * d_at2Idx + 2]);
        // RDKit❗✔️:   RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
        // RDKit❗✔️:                      pos[3 * d_at3Idx + 2]);
        // RDKit❗✔️:   RDGeom::Point3D p4(pos[3 * d_at4Idx], pos[3 * d_at4Idx + 1],
        // RDKit❗✔️:                      pos[3 * d_at4Idx + 2]);

        // BEGIN RDKIT CPP HELPER ForceFields::UFF::Utils::calculateCosY (ForceField/UFF/Utils.cpp:18-40)
        // RDKit❗✔️: double calculateCosY(const RDGeom::Point3D &iPoint,
        // RDKit❗✔️:                      const RDGeom::Point3D &jPoint,
        // RDKit❗✔️:                      const RDGeom::Point3D &kPoint,
        // RDKit❗✔️:                      const RDGeom::Point3D &lPoint) {
        // RDKit❗✔️:   constexpr double zeroTol = 1.0e-16;
        // RDKit❗✔️:   RDGeom::Point3D rJI = iPoint - jPoint;
        // RDKit❗✔️:   RDGeom::Point3D rJK = kPoint - jPoint;
        // RDKit❗✔️:   RDGeom::Point3D rJL = lPoint - jPoint;
        // RDKit❗✔️:   auto l2JI = rJI.lengthSq();
        // RDKit❗✔️:   auto l2JK = rJK.lengthSq();
        // RDKit❗✔️:   auto l2JL = rJL.lengthSq();
        // RDKit❗✔️:   if (l2JI < zeroTol || l2JK < zeroTol || l2JL < zeroTol) {
        // RDKit❗✔️:     return 0.0;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   RDGeom::Point3D n = rJI.crossProduct(rJK);
        // RDKit❗✔️:   n /= (sqrt(l2JI) * sqrt(l2JK));
        // RDKit❗✔️:   auto l2n = n.lengthSq();
        // RDKit❗✔️:   if (l2n < zeroTol) {
        // RDKit❗✔️:     return 0.0;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return n.dotProduct(rJL) / (sqrt(l2JL) * sqrt(l2n));
        // RDKit❗✔️: }
        // END RDKIT CPP HELPER ForceFields::UFF::Utils::calculateCosY
        let coordinates = context.coordinates();
        let p1_base = self.at1_idx as usize * 3;
        let p2_base = self.at2_idx as usize * 3;
        let p3_base = self.at3_idx as usize * 3;
        let p4_base = self.at4_idx as usize * 3;
        let p1 = Point3 {
            x: coordinates[p1_base],
            y: coordinates[p1_base + 1],
            z: coordinates[p1_base + 2],
        };
        let p2 = Point3 {
            x: coordinates[p2_base],
            y: coordinates[p2_base + 1],
            z: coordinates[p2_base + 2],
        };
        let p3 = Point3 {
            x: coordinates[p3_base],
            y: coordinates[p3_base + 1],
            z: coordinates[p3_base + 2],
        };
        let p4 = Point3 {
            x: coordinates[p4_base],
            y: coordinates[p4_base + 1],
            z: coordinates[p4_base + 2],
        };

        // RDKit❗✔️:   double cosY = Utils::calculateCosY(p1, p2, p3, p4);
        let cos_y = calculate_cos_y(&p1, &p2, &p3, &p4);
        // RDKit❗✔️:   double sinYSq = 1.0 - cosY * cosY;
        let sin_y_sq = 1.0 - cos_y * cos_y;
        // RDKit❗✔️:   double sinY = ((sinYSq > 0.0) ? sqrt(sinYSq) : 0.0);
        let sin_y = if sin_y_sq > 0.0 { sin_y_sq.sqrt() } else { 0.0 };
        // RDKit❗✔️:   // cos(2 * W) = 2 * cos(W) * cos(W) - 1 = 2 * sin(W) * sin(W) - 1
        // RDKit❗✔️:   double cos2W = 2.0 * sinY * sinY - 1.0;
        let cos2_w = 2.0 * sin_y * sin_y - 1.0;
        // RDKit❗✔️:   double res = d_forceConstant * (d_C0 + d_C1 * sinY + d_C2 * cos2W);
        let res = self.force_constant * (self.c0 + self.c1 * sin_y + self.c2 * cos2_w);
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        res
    }

    pub(crate) fn get_grad(
        &self,
        context: &mut EvaluationContext<'_>,
        gradient: &mut [f64],
    ) -> bool {
        // BEGIN RDKIT CPP FUNCTION ForceFields::UFF::InversionContrib::getGrad (ForceField/UFF/Inversion.cpp:71-134)
        // RDKit❗✔️: void InversionContrib::getGrad(double *pos, double *grad) const {
        // RDKit❗✔️:   PRECONDITION(dp_forceField, "no owner");
        // The contribution owns only validated indices and parameter values;
        // no movable ForceField pointer is retained.
        // RDKit❗✔️:   PRECONDITION(pos, "bad vector");
        // RDKit❗✔️:   PRECONDITION(grad, "bad vector");
        // EvaluationContext and the mutable gradient slice are nonnull borrows.
        // Source has no slice-extent checks, so none are added here.
        // RDKit❗✔️:   RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
        // RDKit❗✔️:                      pos[3 * d_at1Idx + 2]);
        // RDKit❗✔️:   RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
        // RDKit❗✔️:                      pos[3 * d_at2Idx + 2]);
        // RDKit❗✔️:   RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
        // RDKit❗✔️:                      pos[3 * d_at3Idx + 2]);
        // RDKit❗✔️:   RDGeom::Point3D p4(pos[3 * d_at4Idx], pos[3 * d_at4Idx + 1],
        // RDKit❗✔️:                      pos[3 * d_at4Idx + 2]);
        let coordinates = context.coordinates();
        let p1_base = self.at1_idx as usize * 3;
        let p2_base = self.at2_idx as usize * 3;
        let p3_base = self.at3_idx as usize * 3;
        let p4_base = self.at4_idx as usize * 3;
        let p1 = Point3 {
            x: coordinates[p1_base],
            y: coordinates[p1_base + 1],
            z: coordinates[p1_base + 2],
        };
        let p2 = Point3 {
            x: coordinates[p2_base],
            y: coordinates[p2_base + 1],
            z: coordinates[p2_base + 2],
        };
        let p3 = Point3 {
            x: coordinates[p3_base],
            y: coordinates[p3_base + 1],
            z: coordinates[p3_base + 2],
        };
        let p4 = Point3 {
            x: coordinates[p4_base],
            y: coordinates[p4_base + 1],
            z: coordinates[p4_base + 2],
        };

        // RDKit❗✔️:   double *g1 = &(grad[3 * d_at1Idx]);
        // RDKit❗✔️:   double *g2 = &(grad[3 * d_at2Idx]);
        // RDKit❗✔️:   double *g3 = &(grad[3 * d_at3Idx]);
        // RDKit❗✔️:   double *g4 = &(grad[3 * d_at4Idx]);
        let g1 = self.at1_idx as usize * 3;
        let g2 = self.at2_idx as usize * 3;
        let g3 = self.at3_idx as usize * 3;
        let g4 = self.at4_idx as usize * 3;

        // BEGIN RDKIT CPP HELPER RDGeom::operator- (Geometry/point.cpp:64-70)
        // RDKit❗✔️: Point3D operator-(const Point3D &p1, const Point3D &p2) {
        // RDKit❗✔️:   Point3D res;
        // RDKit❗✔️:   res.x = p1.x - p2.x;
        // RDKit❗✔️:   res.y = p1.y - p2.y;
        // RDKit❗✔️:   res.z = p1.z - p2.z;
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // END RDKIT CPP HELPER RDGeom::operator-
        // RDKit❗✔️:   RDGeom::Point3D rJI = p1 - p2;
        // RDKit❗✔️:   RDGeom::Point3D rJK = p3 - p2;
        // RDKit❗✔️:   RDGeom::Point3D rJL = p4 - p2;
        let mut r_ji = Point3::difference(&p1, &p2);
        let mut r_jk = Point3::difference(&p3, &p2);
        let mut r_jl = Point3::difference(&p4, &p2);

        // BEGIN RDKIT CPP HELPER RDGeom::Point3D::length (Geometry/point.h:158-161)
        // RDKit❗✔️: double length() const override {
        // RDKit❗✔️:   double res = x * x + y * y + z * z;
        // RDKit❗✔️:   return sqrt(res);
        // RDKit❗✔️: }
        // END RDKIT CPP HELPER RDGeom::Point3D::length
        // RDKit❗✔️:   double dJI = rJI.length();
        // RDKit❗✔️:   double dJK = rJK.length();
        // RDKit❗✔️:   double dJL = rJL.length();
        let d_ji = r_ji.length();
        let d_jk = r_jk.length();
        let d_jl = r_jl.length();

        // BEGIN RDKIT CPP HELPER ForceFields::UFF::isDoubleZero (ForceField/UFF/Params.h:29-31)
        // RDKit❗✔️: inline bool isDoubleZero(const double x) {
        // RDKit❗✔️:   return ((x < 1.0e-10) && (x > -1.0e-10));
        // RDKit❗✔️: }
        // END RDKIT CPP HELPER ForceFields::UFF::isDoubleZero
        // RDKit❗✔️:   if (isDoubleZero(dJI) || isDoubleZero(dJK) || isDoubleZero(dJL)) {
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        if is_double_zero(d_ji) || is_double_zero(d_jk) || is_double_zero(d_jl) {
            // The packed source uses this as a return from its whole loop;
            // carry that control result to its owner without changing the
            // single-term gradient mutations (none have occurred yet).
            return false;
        }

        // BEGIN RDKIT CPP HELPER RDGeom::Point3D::operator/= (Geometry/point.h:132-137)
        // RDKit❗✔️: constexpr Point3D &operator/=(double scale) {
        // RDKit❗✔️:   x /= scale;
        // RDKit❗✔️:   y /= scale;
        // RDKit❗✔️:   z /= scale;
        // RDKit❗✔️:   return *this;
        // RDKit❗✔️: }
        // END RDKIT CPP HELPER RDGeom::Point3D::operator/=
        // RDKit❗✔️:   rJI /= dJI;
        // RDKit❗✔️:   rJK /= dJK;
        // RDKit❗✔️:   rJL /= dJL;
        r_ji.x /= d_ji;
        r_ji.y /= d_ji;
        r_ji.z /= d_ji;
        r_jk.x /= d_jk;
        r_jk.y /= d_jk;
        r_jk.z /= d_jk;
        r_jl.x /= d_jl;
        r_jl.y /= d_jl;
        r_jl.z /= d_jl;

        // BEGIN RDKIT CPP HELPER RDGeom::Point3D::operator- (Geometry/point.h:139-145)
        // RDKit❗✔️: constexpr Point3D operator-() const {
        // RDKit❗✔️:   Point3D res(x, y, z);
        // RDKit❗✔️:   res.x *= -1.0;
        // RDKit❗✔️:   res.y *= -1.0;
        // RDKit❗✔️:   res.z *= -1.0;
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // END RDKIT CPP HELPER RDGeom::Point3D::operator-
        // BEGIN RDKIT CPP HELPER RDGeom::Point3D::crossProduct (Geometry/point.h:228-234)
        // RDKit❗✔️: constexpr Point3D crossProduct(const Point3D &other) const {
        // RDKit❗✔️:   Point3D res;
        // RDKit❗✔️:   res.x = y * (other.z) - z * (other.y);
        // RDKit❗✔️:   res.y = -x * (other.z) + z * (other.x);
        // RDKit❗✔️:   res.z = x * (other.y) - y * (other.x);
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // END RDKIT CPP HELPER RDGeom::Point3D::crossProduct
        // RDKit❗✔️:   RDGeom::Point3D n = (-rJI).crossProduct(rJK);
        // RDKit❗✔️:   n /= n.length();
        let negated_r_ji = Point3 {
            x: r_ji.x * -1.0,
            y: r_ji.y * -1.0,
            z: r_ji.z * -1.0,
        };
        let mut n = negated_r_ji.cross_product(&r_jk);
        let n_length = n.length();
        n.x /= n_length;
        n.y /= n_length;
        n.z /= n_length;

        // BEGIN RDKIT CPP HELPER RDGeom::Point3D::dotProduct (Geometry/point.h:169-172)
        // RDKit❗✔️: constexpr double dotProduct(const Point3D &other) const {
        // RDKit❗✔️:   double res = x * (other.x) + y * (other.y) + z * (other.z);
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // END RDKIT CPP HELPER RDGeom::Point3D::dotProduct
        // RDKit❗✔️:   double cosY = n.dotProduct(rJL);
        // RDKit❗✔️:   clipToOne(cosY);
        let mut cos_y = n.x * r_jl.x + n.y * r_jl.y + n.z * r_jl.z;

        // BEGIN RDKIT CPP HELPER ForceFields::UFF::clipToOne (ForceField/UFF/Params.h:32)
        // RDKit❗✔️: inline void clipToOne(double &x) { x = std::clamp(x, -1.0, 1.0); }
        // END RDKIT CPP HELPER ForceFields::UFF::clipToOne
        clip_to_one(&mut cos_y);
        // RDKit❗✔️:   double sinYSq = 1.0 - cosY * cosY;
        // RDKit❗✔️:   double sinY = std::max(sqrt(sinYSq), 1.0e-8);
        let sin_y_sq = 1.0 - cos_y * cos_y;
        let sin_y_unfloored = sin_y_sq.sqrt();
        // C++ std::max returns its first operand when `a < b` is false. This
        // ordered comparison preserves a first-operand NaN; f64::max would not.
        let sin_y = if sin_y_unfloored < 1.0e-8 {
            1.0e-8
        } else {
            sin_y_unfloored
        };
        // RDKit❗✔️:   double cosTheta = rJI.dotProduct(rJK);
        // RDKit❗✔️:   clipToOne(cosTheta);
        let mut cos_theta = r_ji.x * r_jk.x + r_ji.y * r_jk.y + r_ji.z * r_jk.z;
        clip_to_one(&mut cos_theta);
        // RDKit❗✔️:   double sinThetaSq = 1.0 - cosTheta * cosTheta;
        // RDKit❗✔️:   double sinTheta = std::max(sqrt(sinThetaSq), 1.0e-8);
        let sin_theta_sq = 1.0 - cos_theta * cos_theta;
        let sin_theta_unfloored = sin_theta_sq.sqrt();
        let sin_theta = if sin_theta_unfloored < 1.0e-8 {
            1.0e-8
        } else {
            sin_theta_unfloored
        };

        // RDKit❗✔️:   // sin(2 * W) = 2 * sin(W) * cos(W) = 2 * cos(Y) * sin(Y)
        // RDKit❗✔️:   double dE_dW = -d_forceConstant * (d_C1 * cosY - 4.0 * d_C2 * cosY * sinY);
        let de_dw = -self.force_constant * (self.c1 * cos_y - 4.0 * self.c2 * cos_y * sin_y);
        // RDKit❗✔️:   RDGeom::Point3D t1 = rJL.crossProduct(rJK);
        // RDKit❗✔️:   RDGeom::Point3D t2 = rJI.crossProduct(rJL);
        // RDKit❗✔️:   RDGeom::Point3D t3 = rJK.crossProduct(rJI);
        let t1 = r_jl.cross_product(&r_jk);
        let t2 = r_ji.cross_product(&r_jl);
        let t3 = r_jk.cross_product(&r_ji);
        // RDKit❗✔️:   double term1 = sinY * sinTheta;
        // RDKit❗✔️:   double term2 = cosY / (sinY * sinThetaSq);
        let term1 = sin_y * sin_theta;
        let term2 = cos_y / (sin_y * sin_theta_sq);

        // RDKit❗✔️:   double tg1[3] = {(t1.x / term1 - (rJI.x - rJK.x * cosTheta) * term2) / dJI,
        // RDKit❗✔️:                    (t1.y / term1 - (rJI.y - rJK.y * cosTheta) * term2) / dJI,
        // RDKit❗✔️:                    (t1.z / term1 - (rJI.z - rJK.z * cosTheta) * term2) / dJI};
        let tg1 = [
            (t1.x / term1 - (r_ji.x - r_jk.x * cos_theta) * term2) / d_ji,
            (t1.y / term1 - (r_ji.y - r_jk.y * cos_theta) * term2) / d_ji,
            (t1.z / term1 - (r_ji.z - r_jk.z * cos_theta) * term2) / d_ji,
        ];
        // RDKit❗✔️:   double tg3[3] = {(t2.x / term1 - (rJK.x - rJI.x * cosTheta) * term2) / dJK,
        // RDKit❗✔️:                    (t2.y / term1 - (rJK.y - rJI.y * cosTheta) * term2) / dJK,
        // RDKit❗✔️:                    (t2.z / term1 - (rJK.z - rJI.z * cosTheta) * term2) / dJK};
        let tg3 = [
            (t2.x / term1 - (r_jk.x - r_ji.x * cos_theta) * term2) / d_jk,
            (t2.y / term1 - (r_jk.y - r_ji.y * cos_theta) * term2) / d_jk,
            (t2.z / term1 - (r_jk.z - r_ji.z * cos_theta) * term2) / d_jk,
        ];
        // RDKit❗✔️:   double tg4[3] = {(t3.x / term1 - rJL.x * cosY / sinY) / dJL,
        // RDKit❗✔️:                    (t3.y / term1 - rJL.y * cosY / sinY) / dJL,
        // RDKit❗✔️:                    (t3.z / term1 - rJL.z * cosY / sinY) / dJL};
        let tg4 = [
            (t3.x / term1 - r_jl.x * cos_y / sin_y) / d_jl,
            (t3.y / term1 - r_jl.y * cos_y / sin_y) / d_jl,
            (t3.z / term1 - r_jl.z * cos_y / sin_y) / d_jl,
        ];
        // RDKit❗✔️:   for (unsigned int i = 0; i < 3; ++i) {
        // RDKit❗✔️:     g1[i] += dE_dW * tg1[i];
        // RDKit❗✔️:     g2[i] += -dE_dW * (tg1[i] + tg3[i] + tg4[i]);
        // RDKit❗✔️:     g3[i] += dE_dW * tg3[i];
        // RDKit❗✔️:     g4[i] += dE_dW * tg4[i];
        // RDKit❗✔️:   }
        for i in 0..3 {
            gradient[g1 + i] += de_dw * tg1[i];
            gradient[g2 + i] += -de_dw * (tg1[i] + tg3[i] + tg4[i]);
            gradient[g3 + i] += de_dw * tg3[i];
            gradient[g4 + i] += de_dw * tg4[i];
        }
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION ForceFields::UFF::InversionContrib::getGrad
        true
    }
}

impl ForceFieldContribution for InversionContrib {
    #[cfg(test)]
    fn cf3d_frag_accept_test_identity(&self) -> crate::kernel::Cf3dFragAcceptContributionIdentity {
        crate::kernel::Cf3dFragAcceptContributionIdentity::Inversion {
            at1_idx: self.at1_idx,
            at2_idx: self.at2_idx,
            at3_idx: self.at3_idx,
            at4_idx: self.at4_idx,
            force_constant: self.force_constant,
            c0: self.c0,
            c1: self.c1,
            c2: self.c2,
        }
    }

    fn get_energy(
        &self,
        context: &mut EvaluationContext<'_>,
    ) -> Result<f64, ForceFieldKernelError> {
        // BEGIN RDKIT CPP CALL ForceFields::ForceField::calcEnergy(pos)
        // (ForceField/ForceField.cpp:323-324)
        // RDKit❗✔️:     double E = (*contrib)->getEnergy(pos);
        // RDKit❗✔️:     res += E;
        // END RDKIT CPP CALL ForceFields::ForceField::calcEnergy(pos)
        // Behavior marker — RDKit❗✔️: delegate to the existing scalar inversion
        // evaluator without changing its coordinate reads, arithmetic, or IEEE cases.
        // Complexity marker — RDKit✔️✔️: one direct call and `Ok`, with no additional
        // scan, geometry work, or allocation.
        Ok(InversionContrib::get_energy(self, context))
    }

    fn get_grad(
        &self,
        context: &mut EvaluationContext<'_>,
        gradient: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        // BEGIN RDKIT CPP CALL ForceFields::ForceField::calcGrad(pos)
        // (ForceField/ForceField.cpp:361-364)
        // RDKit❗✔️:     (*contrib)->getGrad(pos, grad);
        // END RDKIT CPP CALL ForceFields::ForceField::calcGrad(pos)
        // Behavior marker — RDKit❗✔️: the scalar source callback returns void. The
        // helper's boolean exists only for the separate packed loop, so discard it
        // while preserving its source early return and additive writes.
        // Complexity marker — RDKit✔️✔️: one direct call and `Ok`; no second geometry
        // pass, buffer, or allocation is introduced.
        InversionContrib::get_grad(self, context, gradient);
        Ok(())
    }

    fn copy(&self) -> Box<dyn ForceFieldContribution> {
        // BEGIN RDKIT CPP FUNCTION ForceFields::UFF::InversionContrib::copy
        // (ForceField/UFF/Inversion.h:44-46)
        // RDKit❗✔️: InversionContrib *copy() const override {
        // RDKit❗✔️:   return new InversionContrib(*this);
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION ForceFields::UFF::InversionContrib::copy
        // Behavior marker — RDKit❗✔️: copy every index and coefficient in this
        // owner-free scalar value, including its stored force constant.
        // Complexity marker — RDKit✔️✔️: fixed-size value copy and one trait-owned
        // allocation, matching the source's `new` copy.
        Box::new(*self)
    }
}

#[cfg(test)]
mod tests {
    use std::error::Error as _;

    use super::{InversionContrib, InversionContributionError, InversionIndexArgument};
    use crate::kernel::{
        EvaluationContext, ForceField, cf3d_bld_b05_calc_energy, cf3d_bld_b05_calc_grad,
        cf3d_bld_b05_copy_force_field,
    };

    #[test]
    fn uff_error_e02_inversion_indices_preserve_exact_payloads() {
        // These manually supplied variants supplement, rather than replace,
        // the source-driven constructor range regressions below.
        let errors = [
            InversionContributionError::IndexOutOfRange {
                argument: InversionIndexArgument::First,
                index: 11,
                upper_bound: 41,
            },
            InversionContributionError::IndexOutOfRange {
                argument: InversionIndexArgument::Second,
                index: 17,
                upper_bound: 43,
            },
            InversionContributionError::IndexOutOfRange {
                argument: InversionIndexArgument::Third,
                index: 23,
                upper_bound: 47,
            },
            InversionContributionError::IndexOutOfRange {
                argument: InversionIndexArgument::Fourth,
                index: 31,
                upper_bound: 53,
            },
        ];
        for (error, (expected_argument, expected_index, expected_upper_bound)) in
            errors.into_iter().zip([
                (InversionIndexArgument::First, 11, 41),
                (InversionIndexArgument::Second, 17, 43),
                (InversionIndexArgument::Third, 23, 47),
                (InversionIndexArgument::Fourth, 31, 53),
            ])
        {
            let erased: &(dyn std::error::Error + 'static) = &error;
            let downcast = erased.downcast_ref::<InversionContributionError>().unwrap();
            assert!(std::ptr::eq(downcast, &error));
            match downcast {
                InversionContributionError::IndexOutOfRange {
                    argument,
                    index,
                    upper_bound,
                } => {
                    assert_eq!(*argument, expected_argument);
                    assert_eq!(*index, expected_index);
                    assert_eq!(*upper_bound, expected_upper_bound);
                }
            }
            assert!(erased.source().is_none());
            assert_eq!(error.to_string(), format!("{error:?}"));
        }
    }

    fn contribution_with_indices(
        indices: [u32; 4],
        atomic_num: i32,
        is_c_bound_to_o: bool,
        scale: f64,
    ) -> Result<InversionContrib, InversionContributionError> {
        let mut position1 = Vec::new();
        let mut position2 = Vec::new();
        let mut position3 = Vec::new();
        let mut position4 = Vec::new();
        let positions = [
            position1.as_mut_slice(),
            position2.as_mut_slice(),
            position3.as_mut_slice(),
            position4.as_mut_slice(),
        ];
        InversionContrib::new_with_scale(
            &positions,
            indices[0],
            indices[1],
            indices[2],
            indices[3],
            atomic_num,
            is_c_bound_to_o,
            scale,
        )
    }

    fn contribution(
        atomic_num: i32,
        is_c_bound_to_o: bool,
    ) -> Result<InversionContrib, InversionContributionError> {
        let mut position1 = Vec::new();
        let mut position2 = Vec::new();
        let mut position3 = Vec::new();
        let mut position4 = Vec::new();
        let positions = [
            position1.as_mut_slice(),
            position2.as_mut_slice(),
            position3.as_mut_slice(),
            position4.as_mut_slice(),
        ];
        InversionContrib::new(&positions, 0, 1, 2, 3, atomic_num, is_c_bound_to_o)
    }

    fn energy_for(contribution: &InversionContrib, coordinates: &[f64]) -> f64 {
        let mut distance_cache = [0.0; 10];
        let mut context = EvaluationContext::for_test(coordinates, &mut distance_cache, 4);
        contribution.get_energy(&mut context)
    }

    fn assert_close(actual: f64, expected: f64) {
        let tolerance = 1.0e-12 * expected.abs().max(1.0);
        assert!(
            (actual - expected).abs() <= tolerance,
            "actual {actual:.17e} differs from expected {expected:.17e} by more than {tolerance:.3e}"
        );
    }

    fn gradient_for(
        contribution: &InversionContrib,
        coordinates: &[f64],
        mut gradient: [f64; 12],
    ) -> [f64; 12] {
        let mut distance_cache = [0.0; 10];
        let mut context = EvaluationContext::for_test(coordinates, &mut distance_cache, 4);
        contribution.get_grad(&mut context, &mut gradient);
        gradient
    }

    fn assert_gradient_close(actual: [f64; 12], expected: [f64; 12]) {
        for (actual, expected) in actual.into_iter().zip(expected) {
            assert_close(actual, expected);
        }
    }

    fn rows_from_coordinates(coordinates: &[f64]) -> [[f64; 3]; 4] {
        [
            [coordinates[0], coordinates[1], coordinates[2]],
            [coordinates[3], coordinates[4], coordinates[5]],
            [coordinates[6], coordinates[7], coordinates[8]],
            [coordinates[9], coordinates[10], coordinates[11]],
        ]
    }

    fn field_with_inversions<'a>(
        rows: &'a mut [[f64; 3]; 4],
        source_options: &[(i32, bool)],
    ) -> ForceField<'a> {
        let mut force_field = ForceField::new(3);
        for row in rows {
            force_field.positions_mut().push(row.as_mut_slice());
        }
        for &(atomic_num, is_c_bound_to_o) in source_options {
            let contribution = InversionContrib::new(
                force_field.positions(),
                0,
                1,
                2,
                3,
                atomic_num,
                is_c_bound_to_o,
            )
            .expect("four rows satisfy the source inversion constructor");
            force_field.add_contribution(Box::new(contribution));
        }
        force_field
            .initialize()
            .expect("four 3D positions initialize the force field");
        force_field
    }

    #[test]
    fn cf3d_u14_constructor_preserves_all_source_coefficient_branches_and_scale() {
        // Fixed values from pinned RDKit UFF/Utils.cpp::
        // calcInversionCoefficientsAndForceConstant (6/7/8, 15/33/51/83,
        // and the default switch arm). The oxygen flag affects only 6/7/8.
        let cases = [
            (6, (2.0, 1.0, -1.0, 0.0)),
            (7, (2.0, 1.0, -1.0, 0.0)),
            (8, (2.0, 1.0, -1.0, 0.0)),
            (
                15,
                (
                    4.496661509237746,
                    1.018815687546345,
                    -0.3879761595391661,
                    1.0,
                ),
            ),
            (
                33,
                (
                    4.08682518341991,
                    1.0055752214991891,
                    -0.21119131609399103,
                    1.0,
                ),
            ),
            (
                51,
                (
                    3.979001005763598,
                    1.0032079774467704,
                    -0.16019931202774357,
                    1.0,
                ),
            ),
            (83, (3.6666666666666674, 1.0, -2.4492935982947064e-16, 1.0)),
            (
                16,
                (
                    158068015.25333518,
                    2.999390827019096,
                    -3.999390780625565,
                    1.0,
                ),
            ),
        ];

        for (atomic_num, (unbound_force, expected_c0, expected_c1, expected_c2)) in cases {
            for is_c_bound_to_o in [false, true] {
                let expected_force = if (6..=8).contains(&atomic_num) && is_c_bound_to_o {
                    50.0 / 3.0
                } else {
                    unbound_force
                };
                let actual = contribution(atomic_num, is_c_bound_to_o).unwrap();
                assert_close(actual.force_constant, expected_force);
                assert_close(actual.c0, expected_c0);
                assert_close(actual.c1, expected_c1);
                assert_close(actual.c2, expected_c2);

                let scaled =
                    contribution_with_indices([0, 1, 2, 3], atomic_num, is_c_bound_to_o, 2.5)
                        .unwrap();
                assert_close(scaled.force_constant, 2.5 * expected_force);
                assert_eq!(scaled.c0, actual.c0);
                assert_eq!(scaled.c1, actual.c1);
                assert_eq!(scaled.c2, actual.c2);
            }
        }
    }

    #[test]
    fn cf3d_u14_energy_matches_fixed_planar_and_nonplanar_source_values() {
        // p1/p2/p3 define the z normal. The planar fourth point has cosY=0;
        // the nonplanar unit fourth point has cosY=0.8 and sinY=0.6.
        let planar = [
            1.0, 0.0, 0.0, // p1
            0.0, 0.0, 0.0, // p2
            0.0, 1.0, 0.0, // p3
            1.0, 0.0, 0.0, // p4
        ];
        let nonplanar = [
            1.0, 0.0, 0.0, // p1
            0.0, 0.0, 0.0, // p2
            0.0, 1.0, 0.0, // p3
            0.6, 0.0, 0.8, // p4
        ];
        // These are fixed source-evaluated values for the pinned coefficient
        // families and the source default arm, not computed by the code under test.
        let cases = [
            (6, 2.0, 50.0 / 3.0, 0.0, 0.8, 20.0 / 3.0),
            (7, 2.0, 50.0 / 3.0, 0.0, 0.8, 20.0 / 3.0),
            (8, 2.0, 50.0 / 3.0, 0.0, 0.8, 20.0 / 3.0),
            (
                15,
                4.496661509237746,
                4.496661509237746,
                7.333333333333334,
                2.27544558674968,
                2.27544558674968,
            ),
            (
                33,
                4.08682518341991,
                4.08682518341991,
                7.333333333333333,
                2.447437894208855,
                2.447437894208855,
            ),
            (
                51,
                3.979001005763598,
                3.979001005763598,
                7.333333333333334,
                2.4951853354283395,
                2.4951853354283395,
            ),
            (
                83,
                3.6666666666666674,
                3.6666666666666674,
                7.333333333333334,
                2.64,
                2.64,
            ),
            (
                16,
                158068015.25333518,
                158068015.25333518,
                7.333333333333333,
                50543252.975452304,
                50543252.975452304,
            ),
        ];

        for (
            atomic_num,
            unbound_force,
            oxygen_bound_force,
            planar_expected,
            nonplanar_unbound_expected,
            nonplanar_oxygen_bound_expected,
        ) in cases
        {
            for is_c_bound_to_o in [false, true] {
                let contribution = contribution(atomic_num, is_c_bound_to_o).unwrap();
                let expected_force = if (6..=8).contains(&atomic_num) && is_c_bound_to_o {
                    oxygen_bound_force
                } else {
                    unbound_force
                };
                assert_close(contribution.force_constant, expected_force);
                assert_close(energy_for(&contribution, &planar), planar_expected);
                assert_close(
                    energy_for(&contribution, &nonplanar),
                    if is_c_bound_to_o && (6..=8).contains(&atomic_num) {
                        nonplanar_oxygen_bound_expected
                    } else {
                        nonplanar_unbound_expected
                    },
                );

                let scaled =
                    contribution_with_indices([0, 1, 2, 3], atomic_num, is_c_bound_to_o, 2.5)
                        .unwrap();
                assert_close(energy_for(&scaled, &planar), 2.5 * planar_expected);
                assert_close(
                    energy_for(&scaled, &nonplanar),
                    2.5 * if is_c_bound_to_o && (6..=8).contains(&atomic_num) {
                        nonplanar_oxygen_bound_expected
                    } else {
                        nonplanar_unbound_expected
                    },
                );
            }
        }
    }

    #[test]
    fn cf3d_u14_constructor_keeps_ordered_range_errors_and_accepts_duplicate_indices() {
        // Source URANGE_CHECK calls fail in constructor argument order.
        let cases = [
            ([4, 4, 4, 4], InversionIndexArgument::First, "idx1"),
            ([0, 4, 4, 4], InversionIndexArgument::Second, "idx2"),
            ([0, 1, 4, 4], InversionIndexArgument::Third, "idx3"),
            ([0, 1, 2, 4], InversionIndexArgument::Fourth, "idx4"),
        ];
        for (indices, expected_argument, expected_source_name) in cases {
            let error = contribution_with_indices(indices, 6, false, 1.0).unwrap_err();
            assert_eq!(
                error,
                InversionContributionError::IndexOutOfRange {
                    argument: expected_argument,
                    index: 4,
                    upper_bound: 4,
                }
            );
            assert_eq!(error.source_category(), "Range Error");
            assert_eq!(error.source_message(), expected_source_name);
            assert_eq!(error.range_detail(), (4, 4));
        }

        // The pinned constructor range-checks but does not reject repeated IDs.
        let duplicate = contribution_with_indices([1, 1, 1, 1], 6, false, 1.0).unwrap();
        assert_eq!(
            [
                duplicate.at1_idx,
                duplicate.at2_idx,
                duplicate.at3_idx,
                duplicate.at4_idx,
            ],
            [1, 1, 1, 1]
        );
    }

    #[test]
    fn cf3d_u15_gradient_matches_fixed_source_values_for_all_parameter_options() {
        // Fixed values evaluated from pinned RDKit 2026.03.1
        // ForceField/UFF/Inversion.cpp::InversionContrib::getGrad for this
        // nonplanar geometry and each U14 coefficient branch. The constructor
        // multiplies only the source force constant by `scale`, so each
        // expected derivative scales linearly. The oxygen flag changes only
        // carbon, nitrogen, and oxygen, as in Utils.cpp.
        let coordinates = [
            0.4, 1.2, -0.3, // p1
            -0.2, 0.1, 0.5, // p2
            1.7, -0.4, 0.2, // p3
            0.6, 1.1, 1.4, // p4
        ];
        let source_first_component = [
            (6, -0.262_966_359_151_671_3),
            (7, -0.262_966_359_151_671_3),
            (8, -0.262_966_359_151_671_3),
            (15, -1.166_607_041_307_889_3),
            (33, -0.965_284_608_498_554_7),
            (51, -0.913_139_584_094_318_1),
            (83, -0.764_229_195_103_490_4),
            (16, -116_065_986.286_788_43),
        ];

        for (atomic_num, unbound_expected) in source_first_component {
            for is_c_bound_to_o in [false, true] {
                let expected = if (6..=8).contains(&atomic_num) && is_c_bound_to_o {
                    -2.191_386_326_263_928
                } else {
                    unbound_expected
                };
                for scale in [1.0, 2.5] {
                    let contribution =
                        contribution_with_indices([0, 1, 2, 3], atomic_num, is_c_bound_to_o, scale)
                            .unwrap();
                    let actual = gradient_for(&contribution, &coordinates, [0.0; 12]);
                    assert_close(actual[0], expected * scale);
                }
            }
        }

        // Full fixed source vectors check all Cartesian components and the
        // four atom contributions for both carbon coefficient cases.
        let unbound = contribution(6, false).unwrap();
        assert_gradient_close(
            gradient_for(&unbound, &coordinates, [0.0; 12]),
            [
                -0.262_966_359_151_671_3,
                -0.482_705_371_593_478_9,
                -0.860_944_655_304_787_5,
                0.979_421_110_870_018_2,
                0.984_363_445_830_566_3,
                0.221_466_459_787_311_05,
                -0.089_443_630_611_437_27,
                -0.164_184_198_656_610_38,
                -0.292_835_996_111_417_3,
                -0.627_011_121_106_909_5,
                -0.337_473_875_580_477_15,
                0.932_314_191_628_893_8,
            ],
        );
        let oxygen_bound = contribution(6, true).unwrap();
        assert_gradient_close(
            gradient_for(&oxygen_bound, &coordinates, [0.0; 12]),
            [
                -2.191_386_326_263_928,
                -4.022_544_763_278_991,
                -7.174_538_794_206_564,
                8.161_842_590_583_486,
                8.203_028_715_254_721,
                1.845_553_831_560_925_5,
                -0.745_363_588_428_644,
                -1.368_201_655_471_753_3,
                -2.440_299_967_595_144_4,
                -5.225_092_675_890_913,
                -2.812_282_296_503_976,
                7.769_284_930_240_783,
            ],
        );
    }

    #[test]
    fn cf3d_u15_planar_and_zero_length_early_return_preserve_gradient() {
        // With p4 in the p1-p2-p3 plane, source cosY is zero and
        // dE_dW is zero; getGrad must leave all accumulated entries intact.
        let planar = [
            1.0, 0.0, 0.0, // p1
            0.0, 0.0, 0.0, // p2
            0.0, 1.0, 0.0, // p3
            1.0, 1.0, 0.0, // p4
        ];
        let initial = [1.25; 12];
        let contribution = contribution(6, false).unwrap();
        assert_eq!(gradient_for(&contribution, &planar, initial), initial);

        // The pinned strict isDoubleZero checks return before any gradient
        // update when one of the three source vector lengths is zero.
        let zero_length = [
            0.0, 0.0, 0.0, // p1
            0.0, 0.0, 0.0, // p2
            0.0, 1.0, 0.0, // p3
            0.0, 0.0, 1.0, // p4
        ];
        assert_eq!(gradient_for(&contribution, &zero_length, initial), initial);
    }

    #[test]
    fn cf3d_u15_degenerate_plane_preserves_source_nan_normalization() {
        // Source divides the zero cross-product normal by its zero length
        // without a guard. This collinear p1/p2/p3 geometry has no zero bond
        // lengths, so it reaches that division and propagates NaNs.
        let collinear = [
            1.0, 0.0, 0.0, // p1
            0.0, 0.0, 0.0, // p2
            2.0, 0.0, 0.0, // p3
            0.0, 1.0, 1.0, // p4
        ];
        let contribution = contribution(6, false).unwrap();
        let actual = gradient_for(&contribution, &collinear, [0.0; 12]);
        assert!(actual.into_iter().all(f64::is_nan));
    }

    #[test]
    fn cf3d_u15_gradient_uses_source_atom_order_and_additive_updates() {
        let p1 = [0.4, 1.2, -0.3];
        let p2 = [-0.2, 0.1, 0.5];
        let p3 = [1.7, -0.4, 0.2];
        let p4 = [0.6, 1.1, 1.4];
        // Source indices [2, 0, 3, 1] address p1,p2,p3,p4 respectively.
        let coordinates = [
            p2[0], p2[1], p2[2], p4[0], p4[1], p4[2], p1[0], p1[1], p1[2], p3[0], p3[1], p3[2],
            8.0, 9.0, 10.0,
        ];
        let contribution = contribution_with_indices([2, 0, 3, 1], 6, false, 1.0).unwrap();
        let mut gradient = [1.25; 15];
        let mut distance_cache = [0.0; 10];
        let mut context = EvaluationContext::for_test(&coordinates, &mut distance_cache, 5);
        contribution.get_grad(&mut context, &mut gradient);

        // The exact source vector is remapped to physical indices; all updates
        // use += and the unreferenced fifth atom remains unchanged.
        let expected = [
            1.25 + 0.979_421_110_870_018_2,
            1.25 + 0.984_363_445_830_566_3,
            1.25 + 0.221_466_459_787_311_05,
            1.25 - 0.627_011_121_106_909_5,
            1.25 - 0.337_473_875_580_477_15,
            1.25 + 0.932_314_191_628_893_8,
            1.25 - 0.262_966_359_151_671_3,
            1.25 - 0.482_705_371_593_478_9,
            1.25 - 0.860_944_655_304_787_5,
            1.25 - 0.089_443_630_611_437_27,
            1.25 - 0.164_184_198_656_610_38,
            1.25 - 0.292_835_996_111_417_3,
            1.25,
            1.25,
            1.25,
        ];
        for (actual, expected) in gradient.into_iter().zip(expected) {
            assert_close(actual, expected);
        }
    }

    #[test]
    fn cf3d_bld_b07_field_energy_sums_source_terms_and_copies_them() {
        // Fixed pinned U14 source energies at cosY=0.8: noncarbonyl carbon
        // 0.8, carbon bound to sp2 oxygen 20/3, and phosphorus 2.27544558674968.
        let coordinates = [
            1.0, 0.0, 0.0, // p1
            0.0, 0.0, 0.0, // p2
            0.0, 1.0, 0.0, // p3
            0.6, 0.0, 0.8, // p4
        ];
        let expected = 0.8 + 20.0 / 3.0 + 2.275_445_586_749_68;
        let mut rows = rows_from_coordinates(&coordinates);
        let mut field = field_with_inversions(&mut rows, &[(6, false), (6, true), (15, false)]);
        assert_close(
            cf3d_bld_b05_calc_energy(&mut field, &coordinates)
                .expect("the scalar source callbacks return energy values"),
            expected,
        );

        let mut copied_field = cf3d_bld_b05_copy_force_field(&field);
        let mut copied_rows = rows_from_coordinates(&coordinates);
        for row in &mut copied_rows {
            copied_field.positions_mut().push(row.as_mut_slice());
        }
        copied_field
            .initialize()
            .expect("the copied field accepts its own borrowed rows");
        assert_close(
            cf3d_bld_b05_calc_energy(&mut copied_field, &coordinates)
                .expect("the copied scalar callbacks retain their source values"),
            expected,
        );
    }

    #[test]
    fn cf3d_bld_b07_field_gradient_uses_the_source_vector_additively() {
        let coordinates = [
            0.4, 1.2, -0.3, // p1
            -0.2, 0.1, 0.5, // p2
            1.7, -0.4, 0.2, // p3
            0.6, 1.1, 1.4, // p4
        ];
        let mut rows = rows_from_coordinates(&coordinates);
        let mut field = field_with_inversions(&mut rows, &[(6, false)]);
        let mut gradient = [1.25; 12];
        cf3d_bld_b05_calc_grad(&mut field, &coordinates, &mut gradient)
            .expect("the scalar inversion gradient callback succeeds");
        assert_gradient_close(
            gradient,
            [
                1.25 - 0.262_966_359_151_671_3,
                1.25 - 0.482_705_371_593_478_9,
                1.25 - 0.860_944_655_304_787_5,
                1.25 + 0.979_421_110_870_018_2,
                1.25 + 0.984_363_445_830_566_3,
                1.25 + 0.221_466_459_787_311_05,
                1.25 - 0.089_443_630_611_437_27,
                1.25 - 0.164_184_198_656_610_38,
                1.25 - 0.292_835_996_111_417_3,
                1.25 - 0.627_011_121_106_909_5,
                1.25 - 0.337_473_875_580_477_15,
                1.25 + 0.932_314_191_628_893_8,
            ],
        );
    }

    #[test]
    fn cf3d_bld_b07_constructor_preserves_source_idx_error_order() {
        let cases = [
            ([4, 4, 4, 4], InversionIndexArgument::First, "idx1"),
            ([0, 4, 4, 4], InversionIndexArgument::Second, "idx2"),
            ([0, 1, 4, 4], InversionIndexArgument::Third, "idx3"),
            ([0, 1, 2, 4], InversionIndexArgument::Fourth, "idx4"),
        ];
        for (indices, expected_argument, expected_name) in cases {
            let error = contribution_with_indices(indices, 6, false, 1.0)
                .expect_err("the first out-of-range source argument must fail");
            assert_eq!(
                error,
                InversionContributionError::IndexOutOfRange {
                    argument: expected_argument,
                    index: 4,
                    upper_bound: 4,
                }
            );
            assert_eq!(error.source_category(), "Range Error");
            assert_eq!(error.source_message(), expected_name);
            assert_eq!(error.range_detail(), (4, 4));
        }
    }

    #[test]
    fn cf3d_bld_b07_field_preserves_zero_nan_and_degenerate_geometry() {
        let zero_coordinates = [0.0; 12];
        let mut zero_rows = rows_from_coordinates(&zero_coordinates);
        let mut zero_field = field_with_inversions(&mut zero_rows, &[(6, false)]);
        assert_eq!(
            cf3d_bld_b05_calc_energy(&mut zero_field, &zero_coordinates),
            Ok(0.0)
        );
        let initial = [1.25; 12];
        let mut zero_gradient = initial;
        cf3d_bld_b05_calc_grad(&mut zero_field, &zero_coordinates, &mut zero_gradient)
            .expect("source zero-length early return is a successful scalar callback");
        assert_eq!(zero_gradient, initial);

        // The source energy ternary maps NaN sinYSq to zero, yielding the
        // fixed unbound-carbon value 2.0 rather than a NaN energy.
        let nan_coordinates = [
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.6,
            f64::NAN,
            0.8,
        ];
        let mut nan_rows = rows_from_coordinates(&nan_coordinates);
        let mut nan_field = field_with_inversions(&mut nan_rows, &[(6, false)]);
        assert_eq!(
            cf3d_bld_b05_calc_energy(&mut nan_field, &nan_coordinates),
            Ok(2.0)
        );

        // A nonzero-bond collinear plane reaches the source's unguarded normal
        // normalization and propagates NaNs through the actual field callback.
        let collinear = [
            1.0, 0.0, 0.0, // p1
            0.0, 0.0, 0.0, // p2
            2.0, 0.0, 0.0, // p3
            0.0, 1.0, 1.0, // p4
        ];
        let mut collinear_rows = rows_from_coordinates(&collinear);
        let mut collinear_field = field_with_inversions(&mut collinear_rows, &[(6, false)]);
        let mut gradient = [0.0; 12];
        cf3d_bld_b05_calc_grad(&mut collinear_field, &collinear, &mut gradient)
            .expect("the source callback propagates degenerate-plane NaNs as values");
        assert!(gradient.into_iter().all(f64::is_nan));
    }
}
