// Copyright (C) 2004-2006 Rational Discovery LLC
//
// This file is part of the RDKit-derived force-field implementation and is
// covered by the BSD license in the pinned RDKit source tree.

use super::bond::{BondMathError, calc_bond_rest_length};
use super::params::{AtomicParams, PARAMS_G, clip_to_one};
use crate::geometry::Point3;
use crate::kernel::{
    AngleIndexArgument, EvaluationContext, ForceFieldContribution, ForceFieldKernelError,
};

const ANGLE_CORRECTION_THRESHOLD: f64 = 0.8660;

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) enum AngleBendError {
    DegeneratePoints,
    IndexOutOfRange {
        argument: AngleIndexArgument,
        index: u32,
        upper_bound: usize,
    },
    BondMath(BondMathError),
    Distance(ForceFieldKernelError),
    BadOrder {
        order: u32,
    },
}

impl AngleBendError {
    const fn source_category(self) -> &'static str {
        match self {
            Self::DegeneratePoints | Self::BondMath(_) | Self::BadOrder { .. } => {
                "Pre-condition Violation"
            }
            Self::IndexOutOfRange { .. } => "Range Error",
            Self::Distance(
                ForceFieldKernelError::IndexOutOfRange { .. }
                | ForceFieldKernelError::AngleIndexOutOfRange { .. },
            ) => "Range Error",
            Self::Distance(ForceFieldKernelError::BadIndex) => "Invariant Violation",
            Self::Distance(ForceFieldKernelError::TransferPostcondition) => {
                "Post-condition Violation"
            }
            Self::Distance(_) => "Pre-condition Violation",
        }
    }
}

impl From<AngleBendError> for ForceFieldKernelError {
    fn from(error: AngleBendError) -> Self {
        // AngleBend.cpp throws the original owner error; the detached bridge
        // preserves kernel failures directly instead of nesting them here.
        match error {
            AngleBendError::Distance(error) => error,
            AngleBendError::DegeneratePoints => Self::AngleDegeneratePoints,
            AngleBendError::IndexOutOfRange {
                argument,
                index,
                upper_bound,
            } => Self::AngleIndexOutOfRange {
                argument,
                index,
                upper_bound,
            },
            AngleBendError::BondMath(BondMathError::InvalidBondOrder { .. }) => Self::BadBondOrder,
            AngleBendError::BadOrder { order } => Self::AngleBadOrder { order },
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct AngleBendContrib {
    at1_idx: u32,
    at2_idx: u32,
    at3_idx: u32,
    order: u32,
    force_constant: f64,
    c0: f64,
    c1: f64,
    c2: f64,
    theta0: f64,
}

impl AngleBendContrib {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn new(
        positions: &[&mut [f64]],
        idx1: u32,
        idx2: u32,
        idx3: u32,
        bond_order12: f64,
        bond_order23: f64,
        at1_params: &AtomicParams,
        at2_params: &AtomicParams,
        at3_params: &AtomicParams,
        order: u32,
    ) -> Result<Self, AngleBendError> {
        // BEGIN RDKIT CPP FUNCTION ForceFields::UFF::AngleBendContrib::AngleBendContrib (AngleBend.cpp:70-122)
        // RDKit❗✔️: AngleBendContrib::AngleBendContrib(ForceField *owner, unsigned int idx1,
        // RDKit❗✔️:                                    unsigned int idx2, unsigned int idx3,
        // RDKit❗✔️:                                    double bondOrder12, double bondOrder23,
        // RDKit❗✔️:                                    const AtomicParams *at1Params,
        // RDKit❗✔️:                                    const AtomicParams *at2Params,
        // RDKit❗✔️:                                    const AtomicParams *at3Params,
        // RDKit❗✔️:                                    unsigned int order) {
        // RDKit❗✔️:   PRECONDITION(owner, "bad owner");
        // RDKit❗✔️:   PRECONDITION(at1Params, "bad params pointer");
        // RDKit❗✔️:   PRECONDITION(at2Params, "bad params pointer");
        // RDKit❗✔️:   PRECONDITION(at3Params, "bad params pointer");
        // RDKit❗✔️:   PRECONDITION((idx1 != idx2 && idx2 != idx3 && idx1 != idx3),
        // RDKit❗✔️:                "degenerate points");
        // RDKit❗✔️:   URANGE_CHECK(idx1, owner->positions().size());
        // RDKit❗✔️:   URANGE_CHECK(idx2, owner->positions().size());
        // RDKit❗✔️:   URANGE_CHECK(idx3, owner->positions().size());
        // RDKit❗✔️:   // the following is a hack to get decent geometries
        // RDKit❗✔️:   // with 3- and 4-membered rings incorporating sp2 atoms
        // RDKit❗✔️:   d_theta0 = at2Params->theta0;
        // RDKit❗✔️:   if (order >= 30) {
        // RDKit❗✔️:     switch (order) {
        // RDKit❗✔️:       case 30:
        // RDKit❗✔️:         d_theta0 = 150.0 / 180.0 * M_PI;
        // RDKit❗✔️:         break;
        // RDKit❗✔️:       case 35:
        // RDKit❗✔️:         d_theta0 = 60.0 / 180.0 * M_PI;
        // RDKit❗✔️:         break;
        // RDKit❗✔️:       case 40:
        // RDKit❗✔️:         d_theta0 = 135.0 / 180.0 * M_PI;
        // RDKit❗✔️:         break;
        // RDKit❗✔️:       case 45:
        // RDKit❗✔️:         d_theta0 = 90.0 / 180.0 * M_PI;
        // RDKit❗✔️:         break;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     order = 0;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   // end of the hack
        // RDKit❗✔️:   dp_forceField = owner;
        // RDKit❗✔️:   d_at1Idx = idx1;
        // RDKit❗✔️:   d_at2Idx = idx2;
        // RDKit❗✔️:   d_at3Idx = idx3;
        // RDKit❗✔️:   d_order = order;
        // RDKit❗✔️:   d_forceConstant = Utils::calcAngleForceConstant(
        // RDKit❗✔️:       d_theta0, bondOrder12, bondOrder23, at1Params, at2Params, at3Params);
        // RDKit❗✔️:   if (order == 0) {
        // RDKit❗✔️:     double sinTheta0 = sin(d_theta0);
        // RDKit❗✔️:     double cosTheta0 = cos(d_theta0);
        // RDKit❗✔️:     d_C2 = 1. / (4. * std::max(sinTheta0 * sinTheta0, 1e-8));
        // RDKit❗✔️:     d_C1 = -4. * d_C2 * cosTheta0;
        // RDKit❗✔️:     d_C0 = d_C2 * (2. * cosTheta0 * cosTheta0 + 1.);
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION ForceFields::UFF::AngleBendContrib::AngleBendContrib
        // Behavior marker — RDKit❗✔️: fixed constructor/error-order regressions follow in U08.
        // Complexity marker — RDKit✔️✔️: three index checks and fixed scalar work; no allocation.
        // Rust references and the borrowed position slice represent non-null owner/parameter
        // preconditions; this value keeps source indices and scalars, not a raw owner pointer.
        if !(idx1 != idx2 && idx2 != idx3 && idx1 != idx3) {
            return Err(AngleBendError::DegeneratePoints);
        }
        if idx1 as usize >= positions.len() {
            return Err(AngleBendError::IndexOutOfRange {
                argument: AngleIndexArgument::First,
                index: idx1,
                upper_bound: positions.len(),
            });
        }
        if idx2 as usize >= positions.len() {
            return Err(AngleBendError::IndexOutOfRange {
                argument: AngleIndexArgument::Second,
                index: idx2,
                upper_bound: positions.len(),
            });
        }
        if idx3 as usize >= positions.len() {
            return Err(AngleBendError::IndexOutOfRange {
                argument: AngleIndexArgument::Third,
                index: idx3,
                upper_bound: positions.len(),
            });
        }

        let mut order = order;
        let mut theta0 = at2_params.theta0;
        if order >= 30 {
            match order {
                30 => theta0 = 150.0 / 180.0 * std::f64::consts::PI,
                35 => theta0 = 60.0 / 180.0 * std::f64::consts::PI,
                40 => theta0 = 135.0 / 180.0 * std::f64::consts::PI,
                45 => theta0 = 90.0 / 180.0 * std::f64::consts::PI,
                _ => {}
            }
            order = 0;
        }

        let force_constant = calc_angle_force_constant(
            theta0,
            bond_order12,
            bond_order23,
            at1_params,
            at2_params,
            at3_params,
        )
        .map_err(AngleBendError::BondMath)?;
        let mut c0 = 0.0;
        let mut c1 = 0.0;
        let mut c2 = 0.0;
        if order == 0 {
            let sin_theta0 = theta0.sin();
            let cos_theta0 = theta0.cos();
            let sin_theta0_squared = sin_theta0 * sin_theta0;
            // std::max(a, b) returns a when `a < b` is false, including NaN.
            let sin_theta0_squared_floor = if sin_theta0_squared < 1.0e-8 {
                1.0e-8
            } else {
                sin_theta0_squared
            };
            c2 = 1.0 / (4.0 * sin_theta0_squared_floor);
            c1 = -4.0 * c2 * cos_theta0;
            c0 = c2 * (2.0 * cos_theta0 * cos_theta0 + 1.0);
        }

        Ok(Self {
            at1_idx: idx1,
            at2_idx: idx2,
            at3_idx: idx3,
            order,
            force_constant,
            c0,
            c1,
            c2,
            theta0,
        })
    }

    fn get_energy_term(&self, cos_theta: f64, sin_theta_sq: f64) -> Result<f64, AngleBendError> {
        // BEGIN RDKIT CPP FUNCTION ForceFields::UFF::AngleBendContrib::getEnergyTerm (AngleBend.cpp:208-241)
        // RDKit❗✔️: double AngleBendContrib::getEnergyTerm(double cosTheta,
        // RDKit❗✔️:                                        double sinThetaSq) const {
        // RDKit❗✔️:   PRECONDITION(d_order == 0 || d_order == 1 || d_order == 2 || d_order == 3 ||
        // RDKit❗✔️:                    d_order == 4,
        // RDKit❗✔️:                "bad order");
        // RDKit❗✔️:   // cos(2x) = cos^2(x) - sin^2(x);
        // RDKit❗✔️:   double cos2Theta = cosTheta * cosTheta - sinThetaSq;
        // RDKit❗✔️:   double res = 0.0;
        // RDKit❗✔️:   if (d_order == 0) {
        // RDKit❗✔️:     res = d_C0 + d_C1 * cosTheta + d_C2 * cos2Theta;
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     switch (d_order) {
        // RDKit❗✔️:       case 1:
        // RDKit❗✔️:         res = -cosTheta;
        // RDKit❗✔️:         break;
        // RDKit❗✔️:       case 2:
        // RDKit❗✔️:         res = cos2Theta;
        // RDKit❗✔️:         break;
        // RDKit❗✔️:       case 3:
        // RDKit❗✔️:         // cos(3x) = cos^3(x) - 3*cos(x)*sin^2(x)
        // RDKit❗✔️:         res = cosTheta * (cosTheta * cosTheta - 3. * sinThetaSq);
        // RDKit❗✔️:         break;
        // RDKit❗✔️:       case 4:
        // RDKit❗✔️:         // cos(4x) = cos^4(x) - 6*cos^2(x)*sin^2(x)+sin^4(x)
        // RDKit❗✔️:         res = int_pow<4>(cosTheta) - 6. * cosTheta * cosTheta * sinThetaSq +
        // RDKit❗✔️:               sinThetaSq * sinThetaSq;
        // RDKit❗✔️:         break;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     res = 1. - res;
        // RDKit❗✔️:     res /= (double)(d_order * d_order);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION ForceFields::UFF::AngleBendContrib::getEnergyTerm
        // BEGIN RDKIT CPP HELPER ::int_pow (RDGeneral/utils.h:89-112)
        // RDKit❗✔️: template <unsigned n>
        // RDKit❗✔️: inline double int_pow(double x) {
        // RDKit❗✔️:   double half = int_pow<n / 2>(x);
        // RDKit❗✔️:   if (n % 2 == 0) {  // even
        // RDKit❗✔️:     return half * half;
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     return half * half * x;
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // RDKit❗✔️: template <>
        // RDKit❗✔️: inline double int_pow<0>(double) {
        // RDKit❗✔️:   return 1;
        // RDKit❗✔️: }
        // RDKit❗✔️: template <>
        // RDKit❗✔️: inline double int_pow<1>(double x) {
        // RDKit❗✔️:   return x;  // this does a series of muls
        // RDKit❗✔️: }
        // END RDKIT CPP HELPER ::int_pow
        // Behavior marker — RDKit❗✔️: fixed all-order and source-order regressions follow in U08.
        // Complexity marker — RDKit✔️✔️: fixed scalar operations and no allocation.
        if !(self.order == 0
            || self.order == 1
            || self.order == 2
            || self.order == 3
            || self.order == 4)
        {
            return Err(AngleBendError::BadOrder { order: self.order });
        }

        let cos2_theta = cos_theta * cos_theta - sin_theta_sq;
        let mut res = 0.0;
        if self.order == 0 {
            res = self.c0 + self.c1 * cos_theta + self.c2 * cos2_theta;
        } else {
            res = match self.order {
                1 => -cos_theta,
                2 => cos2_theta,
                3 => cos_theta * (cos_theta * cos_theta - 3.0 * sin_theta_sq),
                4 => {
                    let cos_theta_squared = cos_theta * cos_theta;
                    let cos_theta_fourth = cos_theta_squared * cos_theta_squared;
                    cos_theta_fourth - 6.0 * cos_theta * cos_theta * sin_theta_sq
                        + sin_theta_sq * sin_theta_sq
                }
                _ => unreachable!("order was validated above"),
            };
            res = 1.0 - res;
            res /= f64::from(self.order * self.order);
        }
        Ok(res)
    }

    fn get_theta_deriv(&self, cos_theta: f64, sin_theta: f64) -> Result<f64, AngleBendError> {
        // BEGIN RDKIT CPP FUNCTION ForceFields::UFF::AngleBendContrib::getThetaDeriv (AngleBend.cpp:243-281)
        // RDKit❗✔️: double AngleBendContrib::getThetaDeriv(double cosTheta, double sinTheta) const {
        // RDKit❗✔️:   PRECONDITION(d_order == 0 || d_order == 1 || d_order == 2 || d_order == 3 ||
        // RDKit❗✔️:                    d_order == 4,
        // RDKit❗✔️:                "bad order");
        // RDKit❗✔️:   double dE_dTheta = 0.0;
        // RDKit❗✔️:   double sin2Theta = 2. * sinTheta * cosTheta;
        // RDKit❗✔️:   if (d_order == 0) {
        // RDKit❗✔️:     dE_dTheta =
        // RDKit❗✔️:         -1. * d_forceConstant * (d_C1 * sinTheta + 2. * d_C2 * sin2Theta);
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     // E = k/n^2 [1-cos(n theta)]
        // RDKit❗✔️:     // dE = - k/n^2 * d cos(n theta)
        // RDKit❗✔️:
        // RDKit❗✔️:     // these all use:
        // RDKit❗✔️:     // d cos(ax) = -a sin(ax)
        // RDKit❗✔️:
        // RDKit❗✔️:     switch (d_order) {
        // RDKit❗✔️:       case 1:
        // RDKit❗✔️:         dE_dTheta = -sinTheta;
        // RDKit❗✔️:         break;
        // RDKit❗✔️:       case 2:
        // RDKit❗✔️:         // sin(2*x) = 2*cos(x)*sin(x)
        // RDKit❗✔️:         dE_dTheta = sin2Theta;
        // RDKit❗✔️:         break;
        // RDKit❗✔️:       case 3:
        // RDKit❗✔️:         // sin(3*x) = 3*sin(x) - 4*sin^3(x)
        // RDKit❗✔️:         dE_dTheta = sinTheta * (3. - 4. * sinTheta * sinTheta);
        // RDKit❗✔️:         break;
        // RDKit❗✔️:       case 4:
        // RDKit❗✔️:         // sin(4*x) = cos(x)*(4*sin(x) - 8*sin^3(x))
        // RDKit❗✔️:         dE_dTheta = cosTheta * sinTheta * (4. - 8. * sinTheta * sinTheta);
        // RDKit❗✔️:         break;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     dE_dTheta *= d_forceConstant / (double)(d_order);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return dE_dTheta;
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION ForceFields::UFF::AngleBendContrib::getThetaDeriv
        // Behavior marker — RDKit❗✔️: fixed all-order derivative regressions follow in U08.
        // Complexity marker — RDKit✔️✔️: fixed scalar operations and no allocation.
        if !(self.order == 0
            || self.order == 1
            || self.order == 2
            || self.order == 3
            || self.order == 4)
        {
            return Err(AngleBendError::BadOrder { order: self.order });
        }

        let mut de_dtheta = 0.0;
        let sin2_theta = 2.0 * sin_theta * cos_theta;
        if self.order == 0 {
            de_dtheta =
                -1.0 * self.force_constant * (self.c1 * sin_theta + 2.0 * self.c2 * sin2_theta);
        } else {
            de_dtheta = match self.order {
                1 => -sin_theta,
                2 => sin2_theta,
                3 => sin_theta * (3.0 - 4.0 * sin_theta * sin_theta),
                4 => cos_theta * sin_theta * (4.0 - 8.0 * sin_theta * sin_theta),
                _ => unreachable!("order was validated above"),
            };
            de_dtheta *= self.force_constant / f64::from(self.order);
        }
        Ok(de_dtheta)
    }

    fn get_energy(&self, context: &mut EvaluationContext<'_>) -> Result<f64, AngleBendError> {
        // BEGIN RDKIT CPP FUNCTION ForceFields::UFF::AngleBendContrib::getEnergy (AngleBend.cpp:123-160)
        // RDKit❗✔️: double AngleBendContrib::getEnergy(double *pos) const {
        // RDKit❗✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit❗✔️:   PRECONDITION(pos, "bad vector");
        // RDKit❗✔️:
        // RDKit❗✔️:   double dist1 = dp_forceField->distance(d_at1Idx, d_at2Idx, pos);
        // RDKit❗✔️:   double dist2 = dp_forceField->distance(d_at2Idx, d_at3Idx, pos);
        // RDKit❗✔️:
        // RDKit❗✔️:   RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
        // RDKit❗✔️:                      pos[3 * d_at1Idx + 2]);
        // RDKit❗✔️:   RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
        // RDKit❗✔️:                      pos[3 * d_at2Idx + 2]);
        // RDKit❗✔️:   RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
        // RDKit❗✔️:                      pos[3 * d_at3Idx + 2]);
        // RDKit❗✔️:   RDGeom::Point3D p12 = p1 - p2;
        // RDKit❗✔️:   RDGeom::Point3D p32 = p3 - p2;
        // RDKit❗✔️:   double cosTheta = p12.dotProduct(p32) / (dist1 * dist2);
        // RDKit❗✔️:   clipToOne(cosTheta);
        // RDKit❗✔️:   // we need sin^2(theta) to get cos(2*theta), so compute that:
        // RDKit❗✔️:   double sinThetaSq = 1. - cosTheta * cosTheta;
        // RDKit❗✔️:
        // RDKit❗✔️:   double angleTerm = getEnergyTerm(cosTheta, sinThetaSq);
        // RDKit❗✔️:   double res = d_forceConstant * angleTerm;
        // RDKit❗✔️:
        // RDKit❗✔️:   // The original UFF does not include any penalty for angles that are zero
        // RDKit❗✔️:   // degrees.
        // RDKit❗✔️:   //   This can lead to overlapping 1-3 atoms (e.g.e Github #7901), which is
        // RDKit❗✔️:   //   obviously bad. We add an empiricial penalty for angles close to zero
        // RDKit❗✔️:   //   borrowed from OpenBabel such that the energy goes up exponentially if the
        // RDKit❗✔️:   //   angle is less than approx theta0,
        // RDKit❗✔️:   // For the sake of efficiency, we only add the penalty if the angle is less
        // RDKit❗✔️:   // than 30 degrees
        // RDKit❗✔️:   if (d_order && d_order < 5 && cosTheta > ANGLE_CORRECTION_THRESHOLD) {
        // RDKit❗✔️:     auto theta = acos(cosTheta);
        // RDKit❗✔️:     res += exp(-20.0 * (theta - d_theta0 + 0.25));
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION ForceFields::UFF::AngleBendContrib::getEnergy
        // BEGIN RDKIT CPP HELPER ForceFields::ForceField::distance (ForceField.cpp:172-203)
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
        // END RDKIT CPP HELPER ForceFields::ForceField::distance
        // BEGIN RDKIT CPP HELPER RDGeom::operator- (Geometry/point.cpp:64-70)
        // RDKit❗✔️: Point3D operator-(const Point3D &p1, const Point3D &p2) {
        // RDKit❗✔️:   Point3D res;
        // RDKit❗✔️:   res.x = p1.x - p2.x;
        // RDKit❗✔️:   res.y = p1.y - p2.y;
        // RDKit❗✔️:   res.z = p1.z - p2.z;
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // END RDKIT CPP HELPER RDGeom::operator-
        // BEGIN RDKIT CPP HELPER RDGeom::Point3D::dotProduct (Geometry/point.h:169-172)
        // RDKit❗✔️: constexpr double dotProduct(const Point3D &other) const {
        // RDKit❗✔️:   double res = x * (other.x) + y * (other.y) + z * (other.z);
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // END RDKIT CPP HELPER RDGeom::Point3D::dotProduct
        // BEGIN RDKIT CPP HELPER ForceFields::UFF::clipToOne (Params.h:32)
        // RDKit❗✔️: inline void clipToOne(double &x) { x = std::clamp(x, -1.0, 1.0); }
        // END RDKIT CPP HELPER ForceFields::UFF::clipToOne
        // BEGIN RDKIT CPP CONSTANT ANGLE_CORRECTION_THRESHOLD (AngleBend.cpp:22)
        // RDKit❗✔️: constexpr double ANGLE_CORRECTION_THRESHOLD = 0.8660;
        // END RDKIT CPP CONSTANT ANGLE_CORRECTION_THRESHOLD
        // Behavior marker — RDKit❗✔️: fixed energy branches and correction cases are added in U09.
        // Complexity marker — RDKit✔️✔️: two cached O(1) distance lookups and fixed stack scalar/vector work; no allocation.
        // Rust references/context replace non-null source owner/position preconditions.
        let dist1 = context
            .distance(self.at1_idx, self.at2_idx)
            .map_err(AngleBendError::Distance)?;
        let dist2 = context
            .distance(self.at2_idx, self.at3_idx)
            .map_err(AngleBendError::Distance)?;

        let coordinates = context.coordinates();
        let p1 = Point3 {
            x: coordinates[3 * self.at1_idx as usize],
            y: coordinates[3 * self.at1_idx as usize + 1],
            z: coordinates[3 * self.at1_idx as usize + 2],
        };
        let p2 = Point3 {
            x: coordinates[3 * self.at2_idx as usize],
            y: coordinates[3 * self.at2_idx as usize + 1],
            z: coordinates[3 * self.at2_idx as usize + 2],
        };
        let p3 = Point3 {
            x: coordinates[3 * self.at3_idx as usize],
            y: coordinates[3 * self.at3_idx as usize + 1],
            z: coordinates[3 * self.at3_idx as usize + 2],
        };
        let p12 = Point3::difference(&p1, &p2);
        let p32 = Point3::difference(&p3, &p2);
        let mut cos_theta = (p12.x * p32.x + p12.y * p32.y + p12.z * p32.z) / (dist1 * dist2);
        clip_to_one(&mut cos_theta);
        let sin_theta_sq = 1.0 - cos_theta * cos_theta;

        let angle_term = self.get_energy_term(cos_theta, sin_theta_sq)?;
        let mut result = self.force_constant * angle_term;
        if self.order != 0 && self.order < 5 && cos_theta > ANGLE_CORRECTION_THRESHOLD {
            let theta = cos_theta.acos();
            result += (-20.0 * (theta - self.theta0 + 0.25)).exp();
        }
        Ok(result)
    }

    fn get_grad(
        &self,
        context: &mut EvaluationContext<'_>,
        gradient: &mut [f64],
    ) -> Result<(), AngleBendError> {
        // BEGIN RDKIT CPP FUNCTION ForceFields::UFF::AngleBendContrib::getGrad (AngleBend.cpp:162-193)
        // RDKit❗✔️: void AngleBendContrib::getGrad(double *pos, double *grad) const {
        // RDKit❗✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit❗✔️:   PRECONDITION(pos, "bad vector");
        // RDKit❗✔️:   PRECONDITION(grad, "bad vector");
        // RDKit❗✔️:
        // RDKit❗✔️:   double dist[2] = {dp_forceField->distance(d_at1Idx, d_at2Idx, pos),
        // RDKit❗✔️:                     dp_forceField->distance(d_at2Idx, d_at3Idx, pos)};
        // RDKit❗✔️:
        // RDKit❗✔️:   RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
        // RDKit❗✔️:                      pos[3 * d_at1Idx + 2]);
        // RDKit❗✔️:   RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
        // RDKit❗✔️:                      pos[3 * d_at2Idx + 2]);
        // RDKit❗✔️:   RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
        // RDKit❗✔️:                      pos[3 * d_at3Idx + 2]);
        // RDKit❗✔️:   double *g[3] = {&(grad[3 * d_at1Idx]), &(grad[3 * d_at2Idx]),
        // RDKit❗✔️:                   &(grad[3 * d_at3Idx])};
        // RDKit❗✔️:   RDGeom::Point3D r[2] = {(p1 - p2) / dist[0], (p3 - p2) / dist[1]};
        // RDKit❗✔️:   double cosTheta = r[0].dotProduct(r[1]);
        // RDKit❗✔️:   clipToOne(cosTheta);
        // RDKit❗✔️:   double sinThetaSq = 1.0 - cosTheta * cosTheta;
        // RDKit❗✔️:   double sinTheta = std::max(sqrt(sinThetaSq), 1.0e-8);
        // RDKit❗✔️:
        // RDKit❗✔️:   // use the chain rule:
        // RDKit❗✔️:   // dE/dx = dE/dTheta * dTheta/dx
        // RDKit❗✔️:
        // RDKit❗✔️:   // dE/dTheta is independent of cartesians:
        // RDKit❗✔️:   double dE_dTheta = getThetaDeriv(cosTheta, sinTheta);
        // RDKit❗✔️:
        // RDKit❗✔️:   // The original UFF does not include any penalty for angles that are zero
        // RDKit❗✔️:   // degrees.
        // RDKit❗✔️:   //   This can lead to overlapping 1-3 atoms (e.g.e Github #7901), which is
        // RDKit❗✔️:   //   obviously bad. We add an empiricial penalty for angles that are close to zero
        // RDKit❗✔️:   //   borrowed from OpenBabel such that the energy goes up exponentially if
        // RDKit❗✔️:   //   the angle is less than approx theta0
        // RDKit❗✔️:   // For the sake of efficiency, we only add the penalty if the angle is less
        // RDKit❗✔️:   // than 30 degrees
        // RDKit❗✔️:   if (d_order && d_order < 5 && cosTheta > ANGLE_CORRECTION_THRESHOLD) {
        // RDKit❗✔️:     auto theta = acos(cosTheta);
        // RDKit❗✔️:     auto corr = -20.0 * exp(-20.0 * (theta - d_theta0 + 0.25));
        // RDKit❗✔️:     dE_dTheta += corr;
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   Utils::calcAngleBendGrad(r, dist, g, dE_dTheta, cosTheta, sinTheta);
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION ForceFields::UFF::AngleBendContrib::getGrad
        // BEGIN RDKIT CPP HELPER ForceFields::ForceField::distance (ForceField.cpp:172-203)
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
        // END RDKIT CPP HELPER ForceFields::ForceField::distance
        // BEGIN RDKIT CPP HELPER RDGeom::operator- (Geometry/point.cpp:64-70)
        // RDKit❗✔️: Point3D operator-(const Point3D &p1, const Point3D &p2) {
        // RDKit❗✔️:   Point3D res;
        // RDKit❗✔️:   res.x = p1.x - p2.x;
        // RDKit❗✔️:   res.y = p1.y - p2.y;
        // RDKit❗✔️:   res.z = p1.z - p2.z;
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // END RDKIT CPP HELPER RDGeom::operator-
        // BEGIN RDKIT CPP HELPER RDGeom::operator/ (Geometry/point.cpp:80-86)
        // RDKit❗✔️: Point3D operator/(const Point3D &p1, double v) {
        // RDKit❗✔️:   Point3D res;
        // RDKit❗✔️:   res.x = p1.x / v;
        // RDKit❗✔️:   res.y = p1.y / v;
        // RDKit❗✔️:   res.z = p1.z / v;
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // END RDKIT CPP HELPER RDGeom::operator/
        // BEGIN RDKIT CPP HELPER RDGeom::Point3D::dotProduct (Geometry/point.h:169-172)
        // RDKit❗✔️: constexpr double dotProduct(const Point3D &other) const {
        // RDKit❗✔️:   double res = x * (other.x) + y * (other.y) + z * (other.z);
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // END RDKIT CPP HELPER RDGeom::Point3D::dotProduct
        // BEGIN RDKIT CPP HELPER ForceFields::UFF::clipToOne (Params.h:32)
        // RDKit❗✔️: inline void clipToOne(double &x) { x = std::clamp(x, -1.0, 1.0); }
        // END RDKIT CPP HELPER ForceFields::UFF::clipToOne
        // BEGIN RDKIT CPP CONSTANT ANGLE_CORRECTION_THRESHOLD (AngleBend.cpp:22)
        // RDKit❗✔️: constexpr double ANGLE_CORRECTION_THRESHOLD = 0.8660;
        // END RDKIT CPP CONSTANT ANGLE_CORRECTION_THRESHOLD
        // Behavior marker — RDKit❗✔️: fixed gradient, correction, clamp, and additive-write regressions follow in U09.
        // Complexity marker — RDKit✔️✔️: two cached O(1) distance lookups and fixed stack vector/scalar work; no allocation.
        // Rust references/context replace non-null source owner/position/gradient preconditions.
        let dist1 = context
            .distance(self.at1_idx, self.at2_idx)
            .map_err(AngleBendError::Distance)?;
        let dist2 = context
            .distance(self.at2_idx, self.at3_idx)
            .map_err(AngleBendError::Distance)?;

        let coordinates = context.coordinates();
        let p1 = Point3 {
            x: coordinates[3 * self.at1_idx as usize],
            y: coordinates[3 * self.at1_idx as usize + 1],
            z: coordinates[3 * self.at1_idx as usize + 2],
        };
        let p2 = Point3 {
            x: coordinates[3 * self.at2_idx as usize],
            y: coordinates[3 * self.at2_idx as usize + 1],
            z: coordinates[3 * self.at2_idx as usize + 2],
        };
        let p3 = Point3 {
            x: coordinates[3 * self.at3_idx as usize],
            y: coordinates[3 * self.at3_idx as usize + 1],
            z: coordinates[3 * self.at3_idx as usize + 2],
        };
        let p12 = Point3::difference(&p1, &p2);
        let p32 = Point3::difference(&p3, &p2);
        let r = [p12.divided(dist1), p32.divided(dist2)];
        let mut cos_theta = r[0].x * r[1].x + r[0].y * r[1].y + r[0].z * r[1].z;
        clip_to_one(&mut cos_theta);
        let sin_theta_sq = 1.0 - cos_theta * cos_theta;
        let raw_sin_theta = sin_theta_sq.sqrt();
        let sin_theta = source_max_angle_sine(raw_sin_theta);

        let mut de_dtheta = self.get_theta_deriv(cos_theta, sin_theta)?;
        if self.order != 0 && self.order < 5 && cos_theta > ANGLE_CORRECTION_THRESHOLD {
            let theta = cos_theta.acos();
            let correction = -20.0 * (-20.0 * (theta - self.theta0 + 0.25)).exp();
            de_dtheta += correction;
        }

        calc_angle_bend_grad(
            &r,
            &[dist1, dist2],
            gradient,
            [
                self.at1_idx as usize,
                self.at2_idx as usize,
                self.at3_idx as usize,
            ],
            de_dtheta,
            cos_theta,
            sin_theta,
        );
        Ok(())
    }
}

impl ForceFieldContribution for AngleBendContrib {
    #[cfg(test)]
    fn cf3d_frag_accept_test_identity(&self) -> crate::kernel::Cf3dFragAcceptContributionIdentity {
        crate::kernel::Cf3dFragAcceptContributionIdentity::AngleBend {
            at1_idx: self.at1_idx,
            at2_idx: self.at2_idx,
            at3_idx: self.at3_idx,
            order: self.order,
            force_constant: self.force_constant,
            c0: self.c0,
            c1: self.c1,
            c2: self.c2,
            theta0: self.theta0,
        }
    }

    fn get_energy(
        &self,
        context: &mut EvaluationContext<'_>,
    ) -> Result<f64, ForceFieldKernelError> {
        // BEGIN RDKIT CPP CALL ForceFields::ForceField::calcEnergy(pos) (ForceField.cpp:323-324)
        // RDKit❗✔️:     double E = (*contrib)->getEnergy(pos);
        // RDKit❗✔️:     res += E;
        // END RDKIT CPP CALL ForceFields::ForceField::calcEnergy(pos)
        // Behavior marker — RDKit❗✔️: delegate to the source-ordered angle evaluator and
        // preserve the first typed distance or order error through the field loop.
        // Complexity marker — RDKit✔️✔️: one direct term call and a fixed-size error match;
        // no extra allocation or scan.
        AngleBendContrib::get_energy(self, context).map_err(ForceFieldKernelError::from)
    }

    fn get_grad(
        &self,
        context: &mut EvaluationContext<'_>,
        gradient: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        // BEGIN RDKIT CPP CALL ForceFields::ForceField::calcGrad(pos) (ForceField.cpp:363)
        // RDKit❗✔️:     (*contrib)->getGrad(pos, grad);
        // END RDKIT CPP CALL ForceFields::ForceField::calcGrad(pos)
        // Behavior marker — RDKit❗✔️: preserve source distance order, additive writes, and
        // immediate typed failure propagation to the existing field loop.
        // Complexity marker — RDKit✔️✔️: one direct term call and a fixed-size error match;
        // no extra allocation or scan.
        AngleBendContrib::get_grad(self, context, gradient).map_err(ForceFieldKernelError::from)
    }

    fn copy(&self) -> Box<dyn ForceFieldContribution> {
        // BEGIN RDKIT CPP FUNCTION ForceFields::UFF::AngleBendContrib::copy (AngleBend.h:53-55)
        // RDKit❗✔️: AngleBendContrib *copy() const override {
        // RDKit❗✔️:   return new AngleBendContrib(*this);
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION ForceFields::UFF::AngleBendContrib::copy
        // The copied Rust value stores no owner pointer; the caller supplies its
        // current EvaluationContext when this term is evaluated.
        // Behavior marker — RDKit❗✔️: copy every stored angle parameter as the source does.
        // Complexity marker — RDKit✔️✔️: one fixed-size value copy and the trait's Box allocation.
        Box::new(*self)
    }
}

fn source_max_angle_sine(raw_sin_theta: f64) -> f64 {
    // BEGIN RDKIT CPP HELPER std::max for AngleBendContrib::getGrad (AngleBend.cpp:180)
    // RDKit❗✔️: double sinTheta = std::max(sqrt(sinThetaSq), 1.0e-8);
    // END RDKIT CPP HELPER std::max for AngleBendContrib::getGrad
    // Behavior marker — RDKit❗✔️: U09 tests cover the zero floor and NaN-first-argument rule.
    // Complexity marker — RDKit✔️✔️: one comparison and one scalar return; no allocation.
    // C++ std::max returns its first argument when `first < second` is false,
    // including when the first argument is NaN.
    if raw_sin_theta < 1.0e-8 {
        1.0e-8
    } else {
        raw_sin_theta
    }
}

pub(super) fn calc_angle_force_constant(
    theta0: f64,
    bond_order12: f64,
    bond_order23: f64,
    at1_params: &AtomicParams,
    at2_params: &AtomicParams,
    at3_params: &AtomicParams,
) -> Result<f64, BondMathError> {
    // BEGIN RDKIT CPP FUNCTION ForceFields::UFF::Utils::calcAngleForceConstant (AngleBend.cpp:27-43)
    // RDKit✔️✔️: double calcAngleForceConstant(double theta0, double bondOrder12,
    // RDKit✔️✔️:                               double bondOrder23, const AtomicParams *at1Params,
    // RDKit✔️✔️:                               const AtomicParams *at2Params,
    // RDKit✔️✔️:                               const AtomicParams *at3Params) {
    // RDKit✔️✔️:   double cosTheta0 = cos(theta0);
    // RDKit✔️✔️:   double r12 = calcBondRestLength(bondOrder12, at1Params, at2Params);
    // RDKit✔️✔️:   double r23 = calcBondRestLength(bondOrder23, at2Params, at3Params);
    // RDKit✔️✔️:   double r13 = sqrt(r12 * r12 + r23 * r23 - 2. * r12 * r23 * cosTheta0);
    // RDKit✔️✔️:   double beta = 2. * Params::G / (r12 * r23);
    // RDKit✔️✔️:   double preFactor = beta * at1Params->Z1 * at3Params->Z1 / int_pow<5>(r13);
    // RDKit✔️✔️:   double rTerm = r12 * r23;
    // RDKit✔️✔️:   double innerBit =
    // RDKit✔️✔️:       3. * rTerm * (1. - cosTheta0 * cosTheta0) - r13 * r13 * cosTheta0;
    // RDKit✔️✔️:   double res = preFactor * rTerm * innerBit;
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION ForceFields::UFF::Utils::calcAngleForceConstant

    // BEGIN RDKIT CPP HELPER ForceFields::UFF::Utils::calcBondRestLength (BondStretch.cpp:20-39)
    // RDKit✔️✔️: double calcBondRestLength(double bondOrder, const AtomicParams *end1Params,
    // RDKit✔️✔️:                           const AtomicParams *end2Params) {
    // RDKit✔️✔️:   PRECONDITION(bondOrder > 0, "bad bond order");
    // RDKit✔️✔️:   PRECONDITION(end1Params, "bad params pointer");
    // RDKit✔️✔️:   PRECONDITION(end2Params, "bad params pointer");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   double ri = end1Params->r1, rj = end2Params->r1;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // this is the pauling correction:
    // RDKit✔️✔️:   double rBO = -Params::lambda * (ri + rj) * log(bondOrder);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // O'Keefe and Breese electronegativity correction:
    // RDKit✔️✔️:   double Xi = end1Params->GMP_Xi, Xj = end2Params->GMP_Xi;
    // RDKit✔️✔️:   double rEN = ri * rj * (sqrt(Xi) - sqrt(Xj)) * (sqrt(Xi) - sqrt(Xj)) /
    // RDKit✔️✔️:                (Xi * ri + Xj * rj);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   double res = ri + rj + rBO - rEN;
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP HELPER ForceFields::UFF::Utils::calcBondRestLength
    // Behavior marker — RDKit✔️✔️: preserve bond 12 then bond 23 validation and typed error propagation.
    // Complexity marker — RDKit✔️✔️: two fixed-cost scalar bond calculations, constant time and no allocation.
    // Rust references rule out null parameter pointers; the existing helper
    // preserves the source positive-bond-order precondition, including NaN.
    let cos_theta0 = theta0.cos();
    let r12 = calc_bond_rest_length(bond_order12, at1_params, at2_params)?;
    let r23 = calc_bond_rest_length(bond_order23, at2_params, at3_params)?;
    let r13 = (r12 * r12 + r23 * r23 - 2.0 * r12 * r23 * cos_theta0).sqrt();
    let beta = 2.0 * PARAMS_G / (r12 * r23);

    // BEGIN RDKIT CPP HELPER ::int_pow (RDGeneral/utils.h:89-112)
    // RDKit✔️✔️: template <unsigned n>
    // RDKit✔️✔️: inline double int_pow(double x) {
    // RDKit✔️✔️:   double half = int_pow<n / 2>(x);
    // RDKit✔️✔️:   if (n % 2 == 0) {  // even
    // RDKit✔️✔️:     return half * half;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return half * half * x;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: template <>
    // RDKit✔️✔️: inline double int_pow<0>(double) {
    // RDKit✔️✔️:   return 1;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: template <>
    // RDKit✔️✔️: inline double int_pow<1>(double x) {
    // RDKit✔️✔️:   return x;  // this does a series of muls
    // RDKit✔️✔️: }
    // END RDKIT CPP HELPER ::int_pow
    // The n=5 instantiation evaluates (r13*r13)*(r13*r13)*r13 in this order.
    let r13_squared = r13 * r13;
    let r13_fourth = r13_squared * r13_squared;
    let r13_fifth = r13_fourth * r13;
    let pre_factor = beta * at1_params.z1 * at3_params.z1 / r13_fifth;
    let r_term = r12 * r23;
    let inner_bit = 3.0 * r_term * (1.0 - cos_theta0 * cos_theta0) - r13 * r13 * cos_theta0;
    let result = pre_factor * r_term * inner_bit;
    Ok(result)
}

fn calc_angle_bend_grad(
    r: &[Point3; 2],
    dist: &[f64; 2],
    grad: &mut [f64],
    idx: [usize; 3],
    de_dtheta: f64,
    cos_theta: f64,
    sin_theta: f64,
) {
    // BEGIN RDKIT CPP FUNCTION ForceFields::UFF::Utils::calcAngleBendGrad (AngleBend.cpp:45-68)
    // RDKit✔️✔️: void calcAngleBendGrad(RDGeom::Point3D *r, double *dist, double **g,
    // RDKit✔️✔️:                        double &dE_dTheta, double &cosTheta, double &sinTheta) {
    // RDKit✔️✔️:   // -------
    // RDKit✔️✔️:   // dTheta/dx is trickier:
    // RDKit✔️✔️:   double dCos_dS[6] = {1.0 / dist[0] * (r[1].x - cosTheta * r[0].x),
    // RDKit✔️✔️:                        1.0 / dist[0] * (r[1].y - cosTheta * r[0].y),
    // RDKit✔️✔️:                        1.0 / dist[0] * (r[1].z - cosTheta * r[0].z),
    // RDKit✔️✔️:                        1.0 / dist[1] * (r[0].x - cosTheta * r[1].x),
    // RDKit✔️✔️:                        1.0 / dist[1] * (r[0].y - cosTheta * r[1].y),
    // RDKit✔️✔️:                        1.0 / dist[1] * (r[0].z - cosTheta * r[1].z)};
    // RDKit✔️✔️:
    // RDKit✔️✔️:   g[0][0] += dE_dTheta * dCos_dS[0] / (-sinTheta);
    // RDKit✔️✔️:   g[0][1] += dE_dTheta * dCos_dS[1] / (-sinTheta);
    // RDKit✔️✔️:   g[0][2] += dE_dTheta * dCos_dS[2] / (-sinTheta);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   g[1][0] += dE_dTheta * (-dCos_dS[0] - dCos_dS[3]) / (-sinTheta);
    // RDKit✔️✔️:   g[1][1] += dE_dTheta * (-dCos_dS[1] - dCos_dS[4]) / (-sinTheta);
    // RDKit✔️✔️:   g[1][2] += dE_dTheta * (-dCos_dS[2] - dCos_dS[5]) / (-sinTheta);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   g[2][0] += dE_dTheta * dCos_dS[3] / (-sinTheta);
    // RDKit✔️✔️:   g[2][1] += dE_dTheta * dCos_dS[4] / (-sinTheta);
    // RDKit✔️✔️:   g[2][2] += dE_dTheta * dCos_dS[5] / (-sinTheta);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION ForceFields::UFF::Utils::calcAngleBendGrad
    // Behavior marker — RDKit✔️✔️: preserve fixed component and row order, direct arithmetic, and additive writes.
    // Complexity marker — RDKit✔️✔️: six stack scalars and nine writes; constant time and no allocation.
    let dcos_ds = [
        1.0 / dist[0] * (r[1].x - cos_theta * r[0].x),
        1.0 / dist[0] * (r[1].y - cos_theta * r[0].y),
        1.0 / dist[0] * (r[1].z - cos_theta * r[0].z),
        1.0 / dist[1] * (r[0].x - cos_theta * r[1].x),
        1.0 / dist[1] * (r[0].y - cos_theta * r[1].y),
        1.0 / dist[1] * (r[0].z - cos_theta * r[1].z),
    ];
    let g0 = 3 * idx[0];
    let g1 = 3 * idx[1];
    let g2 = 3 * idx[2];

    grad[g0] += de_dtheta * dcos_ds[0] / (-sin_theta);
    grad[g0 + 1] += de_dtheta * dcos_ds[1] / (-sin_theta);
    grad[g0 + 2] += de_dtheta * dcos_ds[2] / (-sin_theta);

    grad[g1] += de_dtheta * (-dcos_ds[0] - dcos_ds[3]) / (-sin_theta);
    grad[g1 + 1] += de_dtheta * (-dcos_ds[1] - dcos_ds[4]) / (-sin_theta);
    grad[g1 + 2] += de_dtheta * (-dcos_ds[2] - dcos_ds[5]) / (-sin_theta);

    grad[g2] += de_dtheta * dcos_ds[3] / (-sin_theta);
    grad[g2 + 1] += de_dtheta * dcos_ds[4] / (-sin_theta);
    grad[g2 + 2] += de_dtheta * dcos_ds[5] / (-sin_theta);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernel::ForceField;
    use crate::uff::params::PARAMS_AMIDE_BOND_ORDER;

    fn atomic_params(r1: f64, theta0: f64, gmp_xi: f64, z1: f64) -> AtomicParams {
        AtomicParams {
            r1,
            theta0,
            x1: 0.0,
            d1: 0.0,
            zeta: 0.0,
            z1,
            v1: 0.0,
            u1: 0.0,
            gmp_xi,
            gmp_hardness: 0.0,
            gmp_radius: 0.0,
        }
    }

    fn assert_invalid_bond_order(result: Result<f64, BondMathError>, expected: f64) {
        match result {
            Err(BondMathError::InvalidBondOrder { bond_order }) => {
                if expected.is_nan() {
                    assert!(bond_order.is_nan());
                } else {
                    assert_eq!(bond_order, expected);
                }
            }
            other => panic!("expected source bond-order error, got {other:?}"),
        }
    }

    fn angle_for_test(
        order: u32,
        theta0: f64,
        indices: [u32; 3],
        point_count: usize,
        bond_order12: f64,
        bond_order23: f64,
    ) -> Result<AngleBendContrib, AngleBendError> {
        let mut owned_positions = vec![Vec::new(); point_count];
        let positions = owned_positions
            .iter_mut()
            .map(Vec::as_mut_slice)
            .collect::<Vec<_>>();
        let at1 = atomic_params(0.729, 0.0, 5.343, 1.912);
        let at2 = atomic_params(0.699, theta0, 6.899, 2.544);
        let at3 = atomic_params(0.757, 109.47 * std::f64::consts::PI / 180.0, 5.343, 1.912);

        AngleBendContrib::new(
            &positions,
            indices[0],
            indices[1],
            indices[2],
            bond_order12,
            bond_order23,
            &at1,
            &at2,
            &at3,
            order,
        )
    }

    #[test]
    fn cf3d_u07_angle_force_constant_matches_rdkit_test_uff3_fixtures() {
        // Fixed RDKit testUFF3 parameter inputs and its asserted rest lengths
        // and amide angle constant (the sp3 angle-constant assertion is
        // commented out in the pinned source).
        let p3 = atomic_params(0.757, 109.47 * std::f64::consts::PI / 180.0, 5.343, 1.912);
        let sp3_rest = calc_bond_rest_length(1.0, &p3, &p3).expect("source-valid bond order");
        assert!((sp3_rest - 1.514).abs() <= 1.0e-12);
        let _source_computes_but_does_not_assert_sp3_force =
            calc_angle_force_constant(p3.theta0, 1.0, 1.0, &p3, &p3, &p3)
                .expect("pinned testUFF3 computes the sp3 force constant");

        let p1 = atomic_params(0.729, 0.0, 5.343, 1.912);
        let p2 = atomic_params(0.699, 120.0 * std::f64::consts::PI / 180.0, 6.899, 2.544);
        let amide_rest = calc_bond_rest_length(PARAMS_AMIDE_BOND_ORDER, &p1, &p2)
            .expect("source-valid amide bond order");
        assert!((amide_rest - 1.357).abs() <= 1.0e-3);
        let single_rest = calc_bond_rest_length(1.0, &p2, &p3).expect("source-valid bond order");
        assert!((single_rest - 1.450).abs() <= 1.0e-3);

        let amide_force =
            calc_angle_force_constant(p2.theta0, PARAMS_AMIDE_BOND_ORDER, 1.0, &p1, &p2, &p3)
                .expect("source-valid amide angle");
        assert!((amide_force - 211.0).abs() <= 1.0e-1);
    }

    #[test]
    fn cf3d_u07_angle_force_constant_preserves_int_pow_five_order() {
        // The pinned int_pow<5> instantiation groups (x*x)*(x*x)*x. This
        // fixed source-formula fixture gives r13 ~= 0.655 and distinguishes
        // that product order from the legacy powi(5) expression.
        let p1 = atomic_params(1.1, 0.0, 1.0, 1.2);
        let p2 = atomic_params(0.9, 0.0, 1.0, 1.0);
        let p3 = atomic_params(0.445, 0.0, 1.0, 1.7);
        let actual = calc_angle_force_constant(0.0, 1.0, 1.0, &p1, &p2, &p3)
            .expect("source-valid bond orders");
        assert_eq!(actual.to_bits(), (-4821.174231826022_f64).to_bits());
    }

    #[test]
    fn cf3d_u07_angle_force_constant_preserves_bond_order_error_order() {
        let p1 = atomic_params(0.757, 0.0, 5.343, 1.912);
        let p2 = atomic_params(0.700, 0.0, 6.899, 2.544);
        let p3 = atomic_params(0.658, 0.0, 8.741, 2.300);

        assert_invalid_bond_order(
            calc_angle_force_constant(1.0, 0.0, -1.0, &p1, &p2, &p3),
            0.0,
        );
        assert_invalid_bond_order(
            calc_angle_force_constant(1.0, f64::NAN, 0.0, &p1, &p2, &p3),
            f64::NAN,
        );

        for invalid_second_order in [0.0, -1.0, f64::NAN] {
            assert_invalid_bond_order(
                calc_angle_force_constant(
                    1.0,
                    PARAMS_AMIDE_BOND_ORDER,
                    invalid_second_order,
                    &p1,
                    &p2,
                    &p3,
                ),
                invalid_second_order,
            );
        }
    }

    #[test]
    fn cf3d_u07_angle_gradient_matches_fixed_source_vector_and_adds() {
        // Fixed Point3, distance, index and scalar inputs from the legacy
        // source-order regression, with a nonzero caller gradient to verify
        // the upstream additive writes.
        let r = [
            Point3 {
                x: 0.6,
                y: 0.8,
                z: -0.1,
            },
            Point3 {
                x: 7.516_508_172_266_185,
                y: 0.5,
                z: 0.84,
            },
        ];
        let dist = [1.0, 1.75];
        let mut grad = [
            -0.5, -0.375, -0.25, -0.125, 0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875,
        ];
        let expected: [f64; 12] = [
            -92.755_247_067_864_3,
            -11.610_188_494_806_621,
            -9.439_378_052_468_7,
            -0.125,
            0.0,
            0.125,
            88.480_701_935_396_25,
            6.244_128_318_182_564,
            10.360_135_574_546_707,
            4.649_545_132_468_043,
            6.116_060_176_624_058,
            0.204_242_477_921_992_8,
        ];

        calc_angle_bend_grad(
            &r,
            &dist,
            &mut grad,
            [2, 0, 3],
            -5.714_193_060_994_148,
            0.0,
            0.486_800_828_948_617,
        );

        for (axis, (actual, expected)) in grad.iter().zip(expected).enumerate() {
            assert_eq!(actual.to_bits(), expected.to_bits(), "flat axis {axis}");
        }
    }

    #[test]
    fn cf3d_u07_angle_gradient_zero_length_and_collinear_inputs_flow_ieee() {
        let r = [
            Point3 {
                x: 1.0,
                y: 0.0,
                z: 0.0,
            },
            Point3 {
                x: 2.0,
                y: 0.0,
                z: 0.0,
            },
        ];

        let mut zero_length_grad = [0.0; 9];
        calc_angle_bend_grad(
            &r,
            &[0.0, 2.0],
            &mut zero_length_grad,
            [0, 1, 2],
            1.0,
            1.0,
            0.5,
        );
        assert!(zero_length_grad.iter().any(|value| !value.is_finite()));

        let collinear_r = [
            Point3 {
                x: 1.0,
                y: 0.0,
                z: 0.0,
            },
            Point3 {
                x: 1.0,
                y: 0.0,
                z: 0.0,
            },
        ];
        let mut collinear_grad = [0.0; 9];
        calc_angle_bend_grad(
            &collinear_r,
            &[1.0, 1.0],
            &mut collinear_grad,
            [0, 1, 2],
            1.0,
            1.0,
            0.0,
        );
        assert!(collinear_grad.iter().all(|value| value.is_nan()));
    }

    #[test]
    fn cf3d_u08_constructor_preserves_order_and_coordination_rewrites() {
        // The constructor source is AngleBend.cpp; its documented default is
        // order 0, and pinned testUFF4 covers ordinary periodic orders.
        let center_theta = 1.2345;
        for order in 0..=4 {
            let angle = angle_for_test(order, center_theta, [0, 1, 2], 3, 1.0, 1.0)
                .expect("source-valid order and indices");
            assert_eq!(angle.order, order);
            assert_eq!(angle.theta0.to_bits(), center_theta.to_bits());
            if order == 0 {
                let sin_theta = center_theta.sin();
                let cos_theta = center_theta.cos();
                let sin_squared = sin_theta * sin_theta;
                let floor = if sin_squared < 1.0e-8 {
                    1.0e-8
                } else {
                    sin_squared
                };
                let c2 = 1.0 / (4.0 * floor);
                let c1 = -4.0 * c2 * cos_theta;
                let c0 = c2 * (2.0 * cos_theta * cos_theta + 1.0);
                assert_eq!(angle.c0.to_bits(), c0.to_bits());
                assert_eq!(angle.c1.to_bits(), c1.to_bits());
                assert_eq!(angle.c2.to_bits(), c2.to_bits());
            }
        }

        let pi = std::f64::consts::PI;
        let at1 = atomic_params(0.729, 0.0, 5.343, 1.912);
        let at2 = atomic_params(0.699, center_theta, 6.899, 2.544);
        let at3 = atomic_params(0.757, 109.47 * std::f64::consts::PI / 180.0, 5.343, 1.912);
        for (source_order, degrees) in [(30, 150.0), (35, 60.0), (40, 135.0), (45, 90.0)] {
            let angle = angle_for_test(source_order, center_theta, [0, 1, 2], 3, 1.0, 1.0)
                .expect("source coordination order");
            let expected_theta = degrees / 180.0 * pi;
            let expected_force =
                calc_angle_force_constant(expected_theta, 1.0, 1.0, &at1, &at2, &at3)
                    .expect("source-valid angle force constant");
            assert_eq!(angle.order, 0);
            assert_eq!(angle.theta0.to_bits(), expected_theta.to_bits());
            assert_eq!(angle.force_constant.to_bits(), expected_force.to_bits());
            assert!(angle.c2.is_finite() && angle.c2 > 0.0);
        }

        let ordinary =
            angle_for_test(0, center_theta, [0, 1, 2], 3, 1.0, 1.0).expect("source harmonic order");
        for unknown_order in [31, u32::MAX] {
            let angle = angle_for_test(unknown_order, center_theta, [0, 1, 2], 3, 1.0, 1.0)
                .expect("unknown >=30 order follows source reset");
            assert_eq!(angle.order, 0);
            assert_eq!(angle.theta0.to_bits(), center_theta.to_bits());
            assert_eq!(
                angle.force_constant.to_bits(),
                ordinary.force_constant.to_bits()
            );
            assert_eq!(angle.c0.to_bits(), ordinary.c0.to_bits());
            assert_eq!(angle.c1.to_bits(), ordinary.c1.to_bits());
            assert_eq!(angle.c2.to_bits(), ordinary.c2.to_bits());
        }
    }

    #[test]
    fn cf3d_u08_constructor_preserves_rdkit_uff3_amide_constant() {
        // Fixed parameter values and 211.0 +/- 1e-1 expectation are pinned
        // testUFF3 assertions, now exercised through the actual constructor.
        let p1 = atomic_params(0.729, 0.0, 5.343, 1.912);
        let p2 = atomic_params(0.699, 120.0 * std::f64::consts::PI / 180.0, 6.899, 2.544);
        let p3 = atomic_params(0.757, 109.47 * std::f64::consts::PI / 180.0, 5.343, 1.912);
        let mut owned_positions = vec![Vec::new(); 3];
        let positions = owned_positions
            .iter_mut()
            .map(Vec::as_mut_slice)
            .collect::<Vec<_>>();
        let angle = AngleBendContrib::new(
            &positions,
            0,
            1,
            2,
            PARAMS_AMIDE_BOND_ORDER,
            1.0,
            &p1,
            &p2,
            &p3,
            0,
        )
        .expect("pinned testUFF3 constructor inputs");

        assert_eq!([angle.at1_idx, angle.at2_idx, angle.at3_idx], [0, 1, 2]);
        assert_eq!(angle.theta0.to_bits(), p2.theta0.to_bits());
        assert!((angle.force_constant - 211.0).abs() <= 1.0e-1);
    }

    #[test]
    fn cf3d_u08_constructor_preserves_source_validation_and_bond_error_order() {
        // AngleBend.cpp checks distinct points, then indices 1/2/3, then
        // calcAngleForceConstant (which checks bond orders 12 then 23).
        let duplicate = angle_for_test(0, 1.0, [3, 3, 9], 3, 0.0, -1.0).unwrap_err();
        assert_eq!(duplicate, AngleBendError::DegeneratePoints);
        assert_eq!(duplicate.source_category(), "Pre-condition Violation");

        for (indices, argument, index) in [
            ([3, 1, 2], AngleIndexArgument::First, 3),
            ([0, 3, 2], AngleIndexArgument::Second, 3),
            ([0, 1, 3], AngleIndexArgument::Third, 3),
        ] {
            let error = angle_for_test(0, 1.0, indices, 3, 0.0, -1.0).unwrap_err();
            assert_eq!(
                error,
                AngleBendError::IndexOutOfRange {
                    argument,
                    index,
                    upper_bound: 3,
                }
            );
            assert_eq!(error.source_category(), "Range Error");
        }

        let first_bond_error = angle_for_test(0, 1.0, [0, 1, 2], 3, 0.0, -1.0).unwrap_err();
        assert_eq!(
            first_bond_error,
            AngleBendError::BondMath(BondMathError::InvalidBondOrder { bond_order: 0.0 })
        );
        assert_eq!(
            first_bond_error.source_category(),
            "Pre-condition Violation"
        );

        for invalid_second_order in [0.0, -1.0, f64::NAN] {
            let error = angle_for_test(
                0,
                1.0,
                [0, 1, 2],
                3,
                PARAMS_AMIDE_BOND_ORDER,
                invalid_second_order,
            )
            .unwrap_err();
            match error {
                AngleBendError::BondMath(BondMathError::InvalidBondOrder { bond_order }) => {
                    if invalid_second_order.is_nan() {
                        assert!(bond_order.is_nan());
                    } else {
                        assert_eq!(bond_order, invalid_second_order);
                    }
                }
                other => panic!("expected second-bond source precondition, got {other:?}"),
            }
            assert_eq!(error.source_category(), "Pre-condition Violation");
        }
    }

    #[test]
    fn cf3d_u08_energy_term_covers_all_orders_and_cosine_endpoints() {
        // Fixed scalar inputs exercise every AngleBend.cpp energy polynomial.
        let pi = std::f64::consts::PI;
        for order in 0..=4 {
            let angle = angle_for_test(order, 1.1, [0, 1, 2], 3, 1.0, 1.0)
                .expect("source-valid periodic order");
            for cos_theta in [0.3123456789012345, 1.0, -1.0] {
                let sin_theta_sq = if cos_theta == 1.0 || cos_theta == -1.0 {
                    0.0
                } else {
                    1.0 - cos_theta * cos_theta
                };
                let cos2_theta = cos_theta * cos_theta - sin_theta_sq;
                let expected = if order == 0 {
                    angle.c0 + angle.c1 * cos_theta + angle.c2 * cos2_theta
                } else {
                    let harmonic = match order {
                        1 => -cos_theta,
                        2 => cos2_theta,
                        3 => cos_theta * (cos_theta * cos_theta - 3.0 * sin_theta_sq),
                        4 => {
                            let cos_squared = cos_theta * cos_theta;
                            let cos_fourth = cos_squared * cos_squared;
                            cos_fourth - 6.0 * cos_theta * cos_theta * sin_theta_sq
                                + sin_theta_sq * sin_theta_sq
                        }
                        _ => unreachable!(),
                    };
                    let mut value = 1.0 - harmonic;
                    value /= f64::from(order * order);
                    value
                };
                let actual = angle
                    .get_energy_term(cos_theta, sin_theta_sq)
                    .expect("source order 0 through 4");
                assert_eq!(actual.to_bits(), expected.to_bits(), "order {order}");
            }
        }

        // The fixed non-binary input pins int_pow<4>'s (x*x)*(x*x) grouping.
        let cos_theta = 0.123456789123456;
        let sin_theta_sq = 1.0 - cos_theta * cos_theta;
        let angle =
            angle_for_test(4, pi / 2.0, [0, 1, 2], 3, 1.0, 1.0).expect("source-valid order 4");
        let cos_squared = cos_theta * cos_theta;
        let cos_fourth = cos_squared * cos_squared;
        let harmonic =
            cos_fourth - 6.0 * cos_theta * cos_theta * sin_theta_sq + sin_theta_sq * sin_theta_sq;
        let expected = (1.0 - harmonic) / 16.0;
        assert_eq!(
            angle
                .get_energy_term(cos_theta, sin_theta_sq)
                .expect("source order 4")
                .to_bits(),
            expected.to_bits()
        );
    }

    #[test]
    fn cf3d_u08_theta_derivative_covers_all_orders_and_cosine_endpoints() {
        // Fixed scalar inputs exercise every AngleBend.cpp derivative branch.
        for order in 0..=4 {
            let angle = angle_for_test(order, 1.1, [0, 1, 2], 3, 1.0, 1.0)
                .expect("source-valid periodic order");
            for cos_theta in [0.3123456789012345_f64, 1.0, -1.0] {
                let sin_theta = if cos_theta == 1.0 || cos_theta == -1.0 {
                    0.0
                } else {
                    (1.0 - cos_theta * cos_theta).sqrt()
                };
                let sin2_theta = 2.0 * sin_theta * cos_theta;
                let expected = if order == 0 {
                    -1.0 * angle.force_constant
                        * (angle.c1 * sin_theta + 2.0 * angle.c2 * sin2_theta)
                } else {
                    let mut value = match order {
                        1 => -sin_theta,
                        2 => sin2_theta,
                        3 => sin_theta * (3.0 - 4.0 * sin_theta * sin_theta),
                        4 => cos_theta * sin_theta * (4.0 - 8.0 * sin_theta * sin_theta),
                        _ => unreachable!(),
                    };
                    value *= angle.force_constant / f64::from(order);
                    value
                };
                let actual = angle
                    .get_theta_deriv(cos_theta, sin_theta)
                    .expect("source order 0 through 4");
                assert_eq!(actual.to_bits(), expected.to_bits(), "order {order}");
            }
        }
    }

    #[test]
    fn cf3d_u08_constructor_preserves_coefficient_floor_and_nan_comparison() {
        // The source std::max floor applies below the threshold and preserves
        // NaN in its first argument, unlike Rust f64::max.
        let zero = angle_for_test(0, 0.0, [0, 1, 2], 3, 1.0, 1.0)
            .expect("zero center angle is source-valid");
        let floor_c2: f64 = 1.0 / (4.0 * 1.0e-8);
        assert_eq!(zero.c2.to_bits(), floor_c2.to_bits());

        let pi_half = angle_for_test(0, std::f64::consts::PI / 2.0, [0, 1, 2], 3, 1.0, 1.0)
            .expect("right-angle center is source-valid");
        assert_eq!(pi_half.c2, 0.25);

        for order in [0, 31] {
            let nan = angle_for_test(order, f64::NAN, [0, 1, 2], 3, 1.0, 1.0)
                .expect("source does not reject a NaN theta0");
            assert_eq!(nan.order, 0);
            assert!(nan.theta0.is_nan());
            assert!(nan.force_constant.is_nan());
            assert!(nan.c0.is_nan() && nan.c1.is_nan() && nan.c2.is_nan());
        }
    }

    #[test]
    fn cf3d_u08_scalar_preconditions_defer_orders_five_through_twenty_nine() {
        for order in [5, 29] {
            let angle = angle_for_test(order, 1.0, [0, 1, 2], 3, 1.0, 1.0)
                .expect("constructor accepts orders below 30");
            assert_eq!(angle.order, order);
            let energy_error = angle.get_energy_term(0.25, 0.5).unwrap_err();
            let derivative_error = angle.get_theta_deriv(0.25, 0.5).unwrap_err();
            assert_eq!(energy_error, AngleBendError::BadOrder { order });
            assert_eq!(derivative_error, AngleBendError::BadOrder { order });
            assert_eq!(energy_error.source_category(), "Pre-condition Violation");
            assert_eq!(
                derivative_error.source_category(),
                "Pre-condition Violation"
            );
        }
    }

    fn assert_close_u09(actual: f64, expected: f64) {
        let tolerance = 2.0e-10 * expected.abs().max(1.0);
        assert!(
            (actual - expected).abs() <= tolerance,
            "expected {expected:.17}, got {actual:.17} (tolerance {tolerance})"
        );
    }

    fn assert_vector_close_u09(actual: &[f64], expected: &[f64]) {
        assert_eq!(actual.len(), expected.len());
        for (index, (&actual, &expected)) in actual.iter().zip(expected).enumerate() {
            assert_close_u09(actual, expected);
            assert!(
                actual.is_finite(),
                "gradient component {index} is not finite"
            );
        }
    }

    fn u09_expected_angle_gradient(de_dtheta: f64, cos_theta: f64, second_arm_y: f64) -> [f64; 9] {
        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
        let dcos_second_x = 1.0 - cos_theta * cos_theta;
        [
            0.0,
            de_dtheta * second_arm_y / (-sin_theta),
            0.0,
            de_dtheta * dcos_second_x / sin_theta,
            de_dtheta * (-second_arm_y + cos_theta * second_arm_y) / (-sin_theta),
            0.0,
            de_dtheta * dcos_second_x / (-sin_theta),
            de_dtheta * (-cos_theta * second_arm_y) / (-sin_theta),
            0.0,
        ]
    }

    fn u09_source_derivative(
        angle: &AngleBendContrib,
        cos_theta: f64,
        sin_theta: f64,
        correction_enabled: bool,
    ) -> f64 {
        let mut derivative = angle
            .get_theta_deriv(cos_theta, sin_theta)
            .expect("constructor-normalized source angle order");
        if correction_enabled
            && angle.order != 0
            && angle.order < 5
            && cos_theta > ANGLE_CORRECTION_THRESHOLD
        {
            let theta = cos_theta.acos();
            derivative += -20.0 * (-20.0 * (theta - angle.theta0 + 0.25)).exp();
        }
        derivative
    }

    #[test]
    fn cf3d_u09_energy_and_gradient_cover_equilibrium_and_distorted_orders() {
        // Fixed AngleBend.cpp::getEnergy/getGrad geometries: a 90-degree
        // harmonic equilibrium and a 60-degree distorted angle. Expected
        // periodic terms and derivatives are the pinned source equations.
        let right_angle = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let sixty_degrees = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 3.0_f64.sqrt() / 2.0, 0.0];
        let equilibrium_energy_ratio = [0.0, 1.0, 0.5, 1.0 / 9.0, 0.0];
        let equilibrium_derivative_ratio = [0.0, -1.0, 0.0, -1.0 / 3.0, 0.0];
        let distorted_energy_ratio = [0.125, 1.5, 0.375, 2.0 / 9.0, 3.0 / 32.0];
        let sqrt_three = 3.0_f64.sqrt();
        let distorted_derivative_ratio = [
            -sqrt_three / 4.0,
            -sqrt_three / 2.0,
            sqrt_three / 4.0,
            0.0,
            -sqrt_three / 8.0,
        ];

        for order in 0..=4 {
            let angle = angle_for_test(order, std::f64::consts::FRAC_PI_2, [0, 1, 2], 3, 1.0, 1.0)
                .expect("source-valid angle order");
            let force_constant = angle.force_constant;

            let mut distance_cache = vec![0.0; 6];
            let mut context = EvaluationContext::for_test(&right_angle, &mut distance_cache, 3);
            let energy = angle
                .get_energy(&mut context)
                .expect("right-angle source energy");
            assert_close_u09(
                energy,
                force_constant * equilibrium_energy_ratio[order as usize],
            );

            let mut gradient = [0.0; 9];
            angle
                .get_grad(&mut context, &mut gradient)
                .expect("right-angle source gradient");
            let de_dtheta = force_constant * equilibrium_derivative_ratio[order as usize];
            let expected = [
                0.0, -de_dtheta, 0.0, de_dtheta, de_dtheta, 0.0, -de_dtheta, 0.0, 0.0,
            ];
            assert_vector_close_u09(&gradient, &expected);

            let mut distance_cache = vec![0.0; 6];
            let mut context = EvaluationContext::for_test(&sixty_degrees, &mut distance_cache, 3);
            let energy = angle
                .get_energy(&mut context)
                .expect("distorted source energy");
            assert_close_u09(
                energy,
                force_constant * distorted_energy_ratio[order as usize],
            );

            let mut gradient = [0.0; 9];
            angle
                .get_grad(&mut context, &mut gradient)
                .expect("distorted source gradient");
            let de_dtheta = force_constant * distorted_derivative_ratio[order as usize];
            let expected = u09_expected_angle_gradient(de_dtheta, 0.5, 3.0_f64.sqrt() / 2.0);
            assert_vector_close_u09(&gradient, &expected);
        }
    }

    #[test]
    fn cf3d_u09_correction_uses_all_orders_and_strict_threshold() {
        // AngleBend.cpp enables the borrowed OpenBabel correction only for
        // orders 1..4 and only when cosTheta is strictly greater than 0.8660.
        let threshold = ANGLE_CORRECTION_THRESHOLD;
        for requested_cosine in [threshold - 1.0e-6, threshold, threshold + 1.0e-6] {
            let requested_sine = (1.0 - requested_cosine * requested_cosine).sqrt();
            let coordinates = [
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                requested_cosine,
                requested_sine,
                0.0,
            ];

            for order in 0..=4 {
                let angle = angle_for_test(order, 1.0, [0, 1, 2], 3, 1.0, 1.0)
                    .expect("source-valid angle order");
                let mut distance_cache = vec![0.0; 6];
                let mut context = EvaluationContext::for_test(&coordinates, &mut distance_cache, 3);

                let energy = angle
                    .get_energy(&mut context)
                    .expect("threshold source energy");
                let first_distance = context.distance(0, 1).expect("cached first arm");
                let second_distance = context.distance(1, 2).expect("cached second arm");
                let cos_theta = coordinates[6] / (first_distance * second_distance);
                let sin_theta_sq = 1.0 - cos_theta * cos_theta;
                if requested_cosine == threshold {
                    assert_eq!(cos_theta.to_bits(), threshold.to_bits());
                }
                let source_term = angle
                    .get_energy_term(cos_theta, sin_theta_sq)
                    .expect("source order 0 through 4");
                let base_energy = angle.force_constant * source_term;
                let correction_applies =
                    order != 0 && order < 5 && cos_theta > ANGLE_CORRECTION_THRESHOLD;
                let correction = if correction_applies {
                    let theta = cos_theta.acos();
                    (-20.0 * (theta - angle.theta0 + 0.25)).exp()
                } else {
                    0.0
                };
                assert_close_u09(energy, base_energy + correction);

                let sin_theta = sin_theta_sq.sqrt();
                let de_dtheta =
                    u09_source_derivative(&angle, cos_theta, sin_theta, correction_applies);
                let unit_second_arm_y = coordinates[7] / second_distance;
                let expected = u09_expected_angle_gradient(de_dtheta, cos_theta, unit_second_arm_y);
                let mut gradient = [0.0; 9];
                angle
                    .get_grad(&mut context, &mut gradient)
                    .expect("threshold source gradient");
                assert_vector_close_u09(&gradient, &expected);
            }
        }
    }

    #[test]
    fn cf3d_u09_ring_orders_reach_harmonic_energy_without_empirical_correction() {
        // The constructor's pinned >=30 rewrite selects the harmonic branch,
        // so the coordinate-level correction predicate remains false.
        let ring_orders = [
            (30, 150.0 / 180.0 * std::f64::consts::PI),
            (35, 60.0 / 180.0 * std::f64::consts::PI),
            (40, 135.0 / 180.0 * std::f64::consts::PI),
            (45, 90.0 / 180.0 * std::f64::consts::PI),
            (u32::MAX, 0.7),
        ];
        let cosine: f64 = 0.95;
        let sine = (1.0 - cosine * cosine).sqrt();
        let coordinates = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, cosine, sine, 0.0];

        for (source_order, expected_theta0) in ring_orders {
            let angle = angle_for_test(source_order, 0.7, [0, 1, 2], 3, 1.0, 1.0)
                .expect("source ring-angle order");
            assert_eq!(angle.order, 0);
            assert_eq!(angle.theta0.to_bits(), expected_theta0.to_bits());
            let mut distance_cache = vec![0.0; 6];
            let mut context = EvaluationContext::for_test(&coordinates, &mut distance_cache, 3);
            let energy = angle
                .get_energy(&mut context)
                .expect("ring harmonic energy");
            let first_distance = context.distance(0, 1).expect("cached first arm");
            let second_distance = context.distance(1, 2).expect("cached second arm");
            let cos_theta = coordinates[6] / (first_distance * second_distance);
            assert!(cos_theta > ANGLE_CORRECTION_THRESHOLD);
            let source_term = angle
                .get_energy_term(cos_theta, 1.0 - cos_theta * cos_theta)
                .expect("constructor rewrote to order zero");
            assert_close_u09(energy, angle.force_constant * source_term);
        }
    }

    #[test]
    fn cf3d_u09_collinear_gradient_preserves_source_sine_floor_and_nan_rule() {
        assert_eq!(source_max_angle_sine(0.0), 1.0e-8);
        assert!(source_max_angle_sine(f64::NAN).is_nan());

        let angle =
            angle_for_test(1, 1.0, [0, 1, 2], 3, 1.0, 1.0).expect("source-valid periodic angle");
        let coordinates = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0];
        let mut distance_cache = vec![0.0; 6];
        let mut context = EvaluationContext::for_test(&coordinates, &mut distance_cache, 3);

        let energy = angle
            .get_energy(&mut context)
            .expect("collinear source energy");
        let correction = (-20.0 * (0.0 - angle.theta0 + 0.25)).exp();
        let source_term = angle.get_energy_term(1.0, 0.0).expect("source order one");
        assert_close_u09(energy, angle.force_constant * source_term + correction);

        let mut gradient = [0.0; 9];
        angle
            .get_grad(&mut context, &mut gradient)
            .expect("collinear source gradient");
        assert_vector_close_u09(&gradient, &[0.0; 9]);
    }

    #[test]
    fn cf3d_u09_distance_errors_keep_source_order_and_do_not_write_gradient() {
        let angle = angle_for_test(1, std::f64::consts::FRAC_PI_2, [0, 1, 2], 3, 1.0, 1.0)
            .expect("constructor-valid angle indices");
        let coordinates = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let mut distance_cache = vec![0.0; 3];
        let mut context = EvaluationContext::for_test(&coordinates, &mut distance_cache, 2);

        let energy_error = angle.get_energy(&mut context).unwrap_err();
        assert!(matches!(
            energy_error,
            AngleBendError::Distance(ForceFieldKernelError::IndexOutOfRange {
                index: 2,
                upper_bound: 2,
                ..
            })
        ));
        assert_eq!(energy_error.source_category(), "Range Error");

        let mut distance_cache = vec![0.0; 3];
        let mut context = EvaluationContext::for_test(&coordinates, &mut distance_cache, 2);
        let mut gradient = [7.0; 6];
        let gradient_error = angle.get_grad(&mut context, &mut gradient).unwrap_err();
        assert!(matches!(
            gradient_error,
            AngleBendError::Distance(ForceFieldKernelError::IndexOutOfRange {
                index: 2,
                upper_bound: 2,
                ..
            })
        ));
        assert_eq!(gradient_error.source_category(), "Range Error");
        assert_eq!(gradient, [7.0; 6]);
    }

    #[test]
    fn cf3d_u09_shared_point_gradients_add_to_existing_values() {
        // Two source contributions share their center and third point. The
        // kernel cache is shared, while each call performs the source's nine
        // additive gradient writes into the caller's existing buffer.
        let first = angle_for_test(1, std::f64::consts::FRAC_PI_2, [0, 1, 2], 4, 1.0, 1.0)
            .expect("first source angle");
        let second = angle_for_test(1, std::f64::consts::FRAC_PI_2, [3, 1, 2], 4, 1.0, 1.0)
            .expect("second source angle");
        let coordinates = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0];
        let mut distance_cache = vec![0.0; 10];
        let mut context = EvaluationContext::for_test(&coordinates, &mut distance_cache, 4);
        let mut gradient = [0.25; 12];

        first
            .get_grad(&mut context, &mut gradient)
            .expect("first shared-point gradient");
        second
            .get_grad(&mut context, &mut gradient)
            .expect("second shared-point gradient");

        let first_k = first.force_constant;
        let second_k = second.force_constant;
        let expected = [
            0.25,
            0.25 + first_k,
            0.25,
            0.25 - first_k + second_k,
            0.25 - first_k - second_k,
            0.25,
            0.25 + first_k - second_k,
            0.25,
            0.25,
            0.25,
            0.25 + second_k,
            0.25,
        ];
        assert_vector_close_u09(&gradient, &expected);
    }

    fn cf3d_bld_b05_new_angle(
        positions: &[&mut [f64]],
        indices: [u32; 3],
        order: u32,
        bond_order12: f64,
        bond_order23: f64,
    ) -> Result<AngleBendContrib, AngleBendError> {
        let at1 = atomic_params(0.5, 0.0, 1.0, 1.0);
        let at2 = atomic_params(0.5, std::f64::consts::FRAC_PI_2, 1.0, 1.0);
        let at3 = atomic_params(0.5, 0.0, 1.0, 1.0);
        AngleBendContrib::new(
            positions,
            indices[0],
            indices[1],
            indices[2],
            bond_order12,
            bond_order23,
            &at1,
            &at2,
            &at3,
            order,
        )
    }

    fn cf3d_bld_b05_field<'a>(points: &'a mut [[f64; 3]], order: u32) -> ForceField<'a> {
        let mut field = ForceField::new(3);
        field
            .positions_mut()
            .extend(points.iter_mut().map(|point| &mut point[..]));
        let angle = {
            let positions = field.positions();
            cf3d_bld_b05_new_angle(positions, [0, 1, 2], order, 1.0, 1.0)
                .expect("source-valid B05 field angle")
        };
        field.add_contribution(Box::new(angle));
        field
            .initialize()
            .expect("three-point source field initializes");
        field
    }

    #[test]
    fn cf3d_bld_b05_force_field_energy_gradient_and_copy_cover_each_angle_order() {
        // Source-fixed equal unit bond radii/electronegativities give unit
        // rest lengths; the pinned force-constant expression yields this
        // value at a 90-degree equilibrium angle.
        const FORCE_CONSTANT: f64 = 352.2028166412076;
        let coordinates = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 3.0_f64.sqrt() / 2.0, 0.0];
        let energy_ratios = [0.125, 1.5, 0.375, 2.0 / 9.0, 3.0 / 32.0];
        let derivative_ratios = [
            -3.0_f64.sqrt() / 4.0,
            -3.0_f64.sqrt() / 2.0,
            3.0_f64.sqrt() / 4.0,
            0.0,
            -3.0_f64.sqrt() / 8.0,
        ];

        for order in 0..=4 {
            let mut points = [
                [1.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.5, 3.0_f64.sqrt() / 2.0, 0.0],
            ];
            let mut field = cf3d_bld_b05_field(&mut points, order);
            let expected_energy = FORCE_CONSTANT * energy_ratios[order as usize];
            let energy = crate::kernel::cf3d_bld_b05_calc_energy(&mut field, &coordinates)
                .expect("source field energy");
            assert_close_u09(energy, expected_energy);

            let source_gradient = u09_expected_angle_gradient(
                FORCE_CONSTANT * derivative_ratios[order as usize],
                0.5,
                3.0_f64.sqrt() / 2.0,
            );
            let mut gradient = [0.25; 9];
            crate::kernel::cf3d_bld_b05_calc_grad(&mut field, &coordinates, &mut gradient)
                .expect("source field gradient");
            let expected_gradient = source_gradient.map(|value| value + 0.25);
            assert_vector_close_u09(&gradient, &expected_gradient);

            let mut copied = crate::kernel::cf3d_bld_b05_copy_force_field(&field);
            let mut copied_points = [
                [1.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.5, 3.0_f64.sqrt() / 2.0, 0.0],
            ];
            copied
                .positions_mut()
                .extend(copied_points.iter_mut().map(|point| &mut point[..]));
            copied
                .initialize()
                .expect("copied source field initializes");
            let copied_energy = crate::kernel::cf3d_bld_b05_calc_energy(&mut copied, &coordinates)
                .expect("copied source field energy");
            assert_close_u09(copied_energy, expected_energy);
            let mut copied_gradient = [0.25; 9];
            crate::kernel::cf3d_bld_b05_calc_grad(&mut copied, &coordinates, &mut copied_gradient)
                .expect("copied source field gradient");
            assert_vector_close_u09(&copied_gradient, &expected_gradient);
        }
    }

    #[test]
    fn cf3d_bld_b05_constructor_errors_keep_source_order_and_typed_mapping() {
        let mut three_points = [[0.0; 3]; 3];
        let three_positions = three_points
            .iter_mut()
            .map(|point| &mut point[..])
            .collect::<Vec<_>>();

        let duplicate = cf3d_bld_b05_new_angle(&three_positions, [0, 1, 1], 0, 1.0, 1.0)
            .expect_err("duplicate source indices fail before parameter work");
        assert_eq!(duplicate.source_category(), "Pre-condition Violation");
        let duplicate_kernel = ForceFieldKernelError::from(duplicate);
        assert_eq!(
            duplicate_kernel,
            ForceFieldKernelError::AngleDegeneratePoints
        );

        let mut two_points = [[0.0; 3]; 2];
        let two_positions = two_points
            .iter_mut()
            .map(|point| &mut point[..])
            .collect::<Vec<_>>();
        let third_index = cf3d_bld_b05_new_angle(&two_positions, [0, 1, 2], 0, 1.0, 1.0)
            .expect_err("source checks third index after first two");
        assert_eq!(third_index.source_category(), "Range Error");
        let third_index_kernel = ForceFieldKernelError::from(third_index);
        assert_eq!(
            third_index_kernel,
            ForceFieldKernelError::AngleIndexOutOfRange {
                argument: AngleIndexArgument::Third,
                index: 2,
                upper_bound: 2,
            }
        );

        for (bond_order12, bond_order23, expected_bad_order) in
            [(0.0, -1.0, 0.0), (1.0, 0.0, 0.0), (1.0, f64::NAN, f64::NAN)]
        {
            let bond_error =
                cf3d_bld_b05_new_angle(&three_positions, [0, 1, 2], 0, bond_order12, bond_order23)
                    .expect_err("source bond-order validation fails");
            assert_eq!(bond_error.source_category(), "Pre-condition Violation");
            match bond_error {
                AngleBendError::BondMath(BondMathError::InvalidBondOrder { bond_order }) => {
                    if expected_bad_order.is_nan() {
                        assert!(bond_order.is_nan());
                    } else {
                        assert_eq!(bond_order, expected_bad_order);
                    }
                }
                other => panic!("expected source bond-order error, got {other:?}"),
            }
            let kernel_error = ForceFieldKernelError::from(bond_error);
            assert_eq!(kernel_error, ForceFieldKernelError::BadBondOrder);
        }
    }

    #[test]
    fn cf3d_bld_b05_field_propagates_callback_errors_without_later_gradient_writes() {
        let coordinates = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 3.0_f64.sqrt() / 2.0, 0.0];
        let mut points = [
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.5, 3.0_f64.sqrt() / 2.0, 0.0],
        ];
        let mut invalid_order_field = {
            let mut field = ForceField::new(3);
            field
                .positions_mut()
                .extend(points.iter_mut().map(|point| &mut point[..]));
            let mut angle = {
                let positions = field.positions();
                cf3d_bld_b05_new_angle(positions, [0, 1, 2], 0, 1.0, 1.0)
                    .expect("constructor accepts source order zero")
            };
            angle.order = 5;
            field.add_contribution(Box::new(angle));
            field
                .initialize()
                .expect("three-point source field initializes");
            field
        };
        let order_error =
            crate::kernel::cf3d_bld_b05_calc_energy(&mut invalid_order_field, &coordinates)
                .expect_err("the callback retains the source bad-order failure");
        assert_eq!(
            order_error,
            ForceFieldKernelError::AngleBadOrder { order: 5 }
        );
        let mut untouched_gradient = [7.0; 9];
        let gradient_error = crate::kernel::cf3d_bld_b05_calc_grad(
            &mut invalid_order_field,
            &coordinates,
            &mut untouched_gradient,
        )
        .expect_err("gradient callback retains source bad-order failure");
        assert_eq!(gradient_error, order_error);
        assert_eq!(untouched_gradient, [7.0; 9]);

        let mut distance_points = [
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.5, 3.0_f64.sqrt() / 2.0, 0.0],
        ];
        let mut distance_field = {
            let mut field = ForceField::new(3);
            field
                .positions_mut()
                .extend(distance_points.iter_mut().map(|point| &mut point[..]));
            let mut angle = {
                let positions = field.positions();
                cf3d_bld_b05_new_angle(positions, [0, 1, 2], 0, 1.0, 1.0)
                    .expect("source-valid constructor indices")
            };
            angle.at3_idx = 3;
            field.add_contribution(Box::new(angle));
            field
                .initialize()
                .expect("three-point source field initializes");
            field
        };
        let expected_distance_error = ForceFieldKernelError::IndexOutOfRange {
            argument: crate::kernel::ForceFieldIndexArgument::J,
            index: 3,
            upper_bound: 3,
        };
        let distance_error =
            crate::kernel::cf3d_bld_b05_calc_energy(&mut distance_field, &coordinates)
                .expect_err("the original kernel distance error propagates unchanged");
        assert_eq!(distance_error, expected_distance_error);
        let mut untouched_gradient = [9.0; 9];
        assert_eq!(
            crate::kernel::cf3d_bld_b05_calc_grad(
                &mut distance_field,
                &coordinates,
                &mut untouched_gradient,
            ),
            Err(expected_distance_error)
        );
        assert_eq!(untouched_gradient, [9.0; 9]);
    }
}
