//! Source-backed RDKit UFF angle bend contribution.

use crate::chemistry::forcefield::core::{ForceField, ForceFieldContrib, ForceFieldVec3};

use super::params::{AtomicParams, clip_to_one};
use super::utils::{UffUtilsError, calc_angle_force_constant};

const ANGLE_CORRECTION_THRESHOLD: f64 = 0.8660;

#[derive(Clone, Debug)]
pub struct AngleBendContrib {
    owner: *const ForceField,
    at1_idx: usize,
    at2_idx: usize,
    at3_idx: usize,
    order: u32,
    force_constant: f64,
    c0: f64,
    c1: f64,
    c2: f64,
    theta0: f64,
}

impl AngleBendContrib {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        owner: &ForceField,
        idx1: usize,
        idx2: usize,
        idx3: usize,
        bond_order12: f64,
        bond_order23: f64,
        at1_params: &AtomicParams,
        at2_params: &AtomicParams,
        at3_params: &AtomicParams,
        order: u32,
    ) -> Result<Self, UffUtilsError> {
        // RDKit✔️✔️: AngleBendContrib::AngleBendContrib(ForceField *owner, unsigned int idx1,
        // RDKit✔️✔️:                                    unsigned int idx2, unsigned int idx3,
        // RDKit✔️✔️:                                    double bondOrder12, double bondOrder23,
        // RDKit✔️✔️:                                    const AtomicParams *at1Params,
        // RDKit✔️✔️:                                    const AtomicParams *at2Params,
        // RDKit✔️✔️:                                    const AtomicParams *at3Params,
        // RDKit✔️✔️:                                    unsigned int order) {
        // RDKit✔️✔️:   PRECONDITION(owner, "bad owner");
        // RDKit✔️✔️:   PRECONDITION(at1Params, "bad params pointer");
        // RDKit✔️✔️:   PRECONDITION(at2Params, "bad params pointer");
        // RDKit✔️✔️:   PRECONDITION(at3Params, "bad params pointer");
        // Rust references reproduce RDKit's non-null preconditions.
        // RDKit✔️✔️:   PRECONDITION((idx1 != idx2 && idx2 != idx3 && idx1 != idx3),
        // RDKit✔️✔️:                "degenerate points");
        assert!(idx1 != idx2 && idx2 != idx3 && idx1 != idx3);
        // RDKit✔️✔️:   URANGE_CHECK(idx1, owner->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx2, owner->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx3, owner->positions().size());
        assert!(idx1 < owner.positions().len());
        assert!(idx2 < owner.positions().len());
        assert!(idx3 < owner.positions().len());

        let mut order = order;
        // RDKit✔️✔️:   // the following is a hack to get decent geometries
        // RDKit✔️✔️:   // with 3- and 4-membered rings incorporating sp2 atoms
        // RDKit✔️✔️:   d_theta0 = at2Params->theta0;
        let mut theta0 = at2_params.theta0;
        // RDKit✔️✔️:   if (order >= 30) {
        if order >= 30 {
            // RDKit✔️✔️:     switch (order) {
            match order {
                // RDKit✔️✔️:       case 30:
                // RDKit✔️✔️:         d_theta0 = 150.0 / 180.0 * M_PI;
                // RDKit✔️✔️:         break;
                30 => theta0 = 150.0 / 180.0 * std::f64::consts::PI,
                // RDKit✔️✔️:       case 35:
                // RDKit✔️✔️:         d_theta0 = 60.0 / 180.0 * M_PI;
                // RDKit✔️✔️:         break;
                35 => theta0 = 60.0 / 180.0 * std::f64::consts::PI,
                // RDKit✔️✔️:       case 40:
                // RDKit✔️✔️:         d_theta0 = 135.0 / 180.0 * M_PI;
                // RDKit✔️✔️:         break;
                40 => theta0 = 135.0 / 180.0 * std::f64::consts::PI,
                // RDKit✔️✔️:       case 45:
                // RDKit✔️✔️:         d_theta0 = 90.0 / 180.0 * M_PI;
                // RDKit✔️✔️:         break;
                45 => theta0 = 90.0 / 180.0 * std::f64::consts::PI,
                // RDKit✔️✔️:     }
                _ => {}
            }
            // RDKit✔️✔️:     order = 0;
            order = 0;
            // RDKit✔️✔️:   }
        }
        // RDKit✔️✔️:   // end of the hack
        // RDKit✔️✔️:   dp_forceField = owner;
        // RDKit✔️✔️:   d_at1Idx = idx1;
        // RDKit✔️✔️:   d_at2Idx = idx2;
        // RDKit✔️✔️:   d_at3Idx = idx3;
        // RDKit✔️✔️:   d_order = order;
        let owner_ptr = owner as *const ForceField;
        // RDKit✔️✔️:   d_forceConstant = Utils::calcAngleForceConstant(
        // RDKit✔️✔️:       d_theta0, bondOrder12, bondOrder23, at1Params, at2Params, at3Params);
        let force_constant =
            calc_angle_force_constant(theta0, bond_order12, bond_order23, at1_params, at2_params, at3_params)?;
        let mut c0 = 0.0;
        let mut c1 = 0.0;
        let mut c2 = 0.0;
        // RDKit✔️✔️:   if (order == 0) {
        if order == 0 {
            // RDKit✔️✔️:     double sinTheta0 = sin(d_theta0);
            // RDKit✔️✔️:     double cosTheta0 = cos(d_theta0);
            let sin_theta0 = theta0.sin();
            let cos_theta0 = theta0.cos();
            // RDKit✔️✔️:     d_C2 = 1. / (4. * std::max(sinTheta0 * sinTheta0, 1e-8));
            c2 = 1.0 / (4.0 * (sin_theta0 * sin_theta0).max(1e-8));
            // RDKit✔️✔️:     d_C1 = -4. * d_C2 * cosTheta0;
            c1 = -4.0 * c2 * cos_theta0;
            // RDKit✔️✔️:     d_C0 = d_C2 * (2. * cosTheta0 * cosTheta0 + 1.);
            // RDKit✔️✔️:   }
            // RDKit✔️✔️: }
            c0 = c2 * (2.0 * cos_theta0 * cos_theta0 + 1.0);
        }

        Ok(Self {
            owner: owner_ptr,
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

    #[must_use]
    pub fn owner(&self) -> *const ForceField {
        self.owner
    }

    #[must_use]
    pub fn at1_idx(&self) -> usize {
        self.at1_idx
    }

    #[must_use]
    pub fn at2_idx(&self) -> usize {
        self.at2_idx
    }

    #[must_use]
    pub fn at3_idx(&self) -> usize {
        self.at3_idx
    }

    #[must_use]
    pub fn order(&self) -> u32 {
        self.order
    }

    #[must_use]
    pub fn force_constant(&self) -> f64 {
        self.force_constant
    }

    #[must_use]
    pub fn c0(&self) -> f64 {
        self.c0
    }

    #[must_use]
    pub fn c1(&self) -> f64 {
        self.c1
    }

    #[must_use]
    pub fn c2(&self) -> f64 {
        self.c2
    }

    #[must_use]
    pub fn theta0(&self) -> f64 {
        self.theta0
    }

    fn force_field(&self) -> &ForceField {
        assert!(!self.owner.is_null(), "no owner");
        // SAFETY: UFF contribs follow the same owner-pointer model as
        // ForceFieldContrib objects in core.rs. Constructors and ForceField::add_contrib
        // store a live owning ForceField pointer before energy/gradient evaluation.
        unsafe { &*self.owner }
    }

    #[must_use]
    pub fn get_energy(&self, pos: &[f64]) -> f64 {
        // RDKit✔️✔️: double AngleBendContrib::getEnergy(double *pos) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        let force_field = self.force_field();

        // RDKit✔️✔️:   double dist1 = dp_forceField->distance(d_at1Idx, d_at2Idx, pos);
        // RDKit✔️✔️:   double dist2 = dp_forceField->distance(d_at2Idx, d_at3Idx, pos);
        let dist1 = force_field.distance(self.at1_idx, self.at2_idx, Some(pos));
        let dist2 = force_field.distance(self.at2_idx, self.at3_idx, Some(pos));

        // RDKit✔️✔️:   RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
        // RDKit✔️✔️:                      pos[3 * d_at1Idx + 2]);
        // RDKit✔️✔️:   RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
        // RDKit✔️✔️:                      pos[3 * d_at2Idx + 2]);
        // RDKit✔️✔️:   RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
        // RDKit✔️✔️:                      pos[3 * d_at3Idx + 2]);
        let p1 = ForceFieldVec3::new(
            pos[3 * self.at1_idx],
            pos[3 * self.at1_idx + 1],
            pos[3 * self.at1_idx + 2],
        );
        let p2 = ForceFieldVec3::new(
            pos[3 * self.at2_idx],
            pos[3 * self.at2_idx + 1],
            pos[3 * self.at2_idx + 2],
        );
        let p3 = ForceFieldVec3::new(
            pos[3 * self.at3_idx],
            pos[3 * self.at3_idx + 1],
            pos[3 * self.at3_idx + 2],
        );
        // RDKit✔️✔️:   RDGeom::Point3D p12 = p1 - p2;
        // RDKit✔️✔️:   RDGeom::Point3D p32 = p3 - p2;
        let p12 = p1 - p2;
        let p32 = p3 - p2;
        // RDKit✔️✔️:   double cosTheta = p12.dotProduct(p32) / (dist1 * dist2);
        let mut cos_theta = p12.dot_product(p32) / (dist1 * dist2);
        // RDKit✔️✔️:   clipToOne(cosTheta);
        clip_to_one(&mut cos_theta);
        // RDKit✔️✔️:   // we need sin^2(theta) to get cos(2*theta), so compute that:
        // RDKit✔️✔️:   double sinThetaSq = 1. - cosTheta * cosTheta;
        let sin_theta_sq = 1.0 - cos_theta * cos_theta;

        // RDKit✔️✔️:   double angleTerm = getEnergyTerm(cosTheta, sinThetaSq);
        let angle_term = self.get_energy_term(cos_theta, sin_theta_sq);
        // RDKit✔️✔️:   double res = d_forceConstant * angleTerm;
        let mut res = self.force_constant * angle_term;

        // RDKit✔️✔️:   // The original UFF does not include any penalty for angles that are zero
        // RDKit✔️✔️:   // degrees.
        // RDKit✔️✔️:   //   This can lead to overlapping 1-3 atoms (e.g.e Github #7901), which is
        // RDKit✔️✔️:   //   obviously bad. We add an empiricial penalty for angles close to zero
        // RDKit✔️✔️:   //   borrowed from OpenBabel such that the energy goes up exponentially if the
        // RDKit✔️✔️:   //   angle is less than approx theta0,
        // RDKit✔️✔️:   // For the sake of efficiency, we only add the penalty if the angle is less
        // RDKit✔️✔️:   // than 30 degrees
        // RDKit✔️✔️:   if (d_order && d_order < 5 && cosTheta > ANGLE_CORRECTION_THRESHOLD) {
        if self.order != 0 && self.order < 5 && cos_theta > ANGLE_CORRECTION_THRESHOLD {
            // RDKit✔️✔️:     auto theta = acos(cosTheta);
            let theta = cos_theta.acos();
            // RDKit✔️✔️:     res += exp(-20.0 * (theta - d_theta0 + 0.25));
            res += (-20.0 * (theta - self.theta0 + 0.25)).exp();
            // RDKit✔️✔️:   }
        }

        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        res
    }

    pub fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
        // RDKit✔️✔️: void AngleBendContrib::getGrad(double *pos, double *grad) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:   PRECONDITION(grad, "bad vector");
        let force_field = self.force_field();

        // RDKit✔️✔️:   double dist[2] = {dp_forceField->distance(d_at1Idx, d_at2Idx, pos),
        // RDKit✔️✔️:                     dp_forceField->distance(d_at2Idx, d_at3Idx, pos)};
        let dist = [
            force_field.distance(self.at1_idx, self.at2_idx, Some(pos)),
            force_field.distance(self.at2_idx, self.at3_idx, Some(pos)),
        ];

        // RDKit✔️✔️:   RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
        // RDKit✔️✔️:                      pos[3 * d_at1Idx + 2]);
        // RDKit✔️✔️:   RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
        // RDKit✔️✔️:                      pos[3 * d_at2Idx + 2]);
        // RDKit✔️✔️:   RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
        // RDKit✔️✔️:                      pos[3 * d_at3Idx + 2]);
        let p1 = ForceFieldVec3::new(
            pos[3 * self.at1_idx],
            pos[3 * self.at1_idx + 1],
            pos[3 * self.at1_idx + 2],
        );
        let p2 = ForceFieldVec3::new(
            pos[3 * self.at2_idx],
            pos[3 * self.at2_idx + 1],
            pos[3 * self.at2_idx + 2],
        );
        let p3 = ForceFieldVec3::new(
            pos[3 * self.at3_idx],
            pos[3 * self.at3_idx + 1],
            pos[3 * self.at3_idx + 2],
        );
        // RDKit✔️✔️:   double *g[3] = {&(grad[3 * d_at1Idx]), &(grad[3 * d_at2Idx]),
        // RDKit✔️✔️:                   &(grad[3 * d_at3Idx])};
        // RDKit✔️✔️:   RDGeom::Point3D r[2] = {(p1 - p2) / dist[0], (p3 - p2) / dist[1]};
        let r = [(p1 - p2) / dist[0], (p3 - p2) / dist[1]];
        // RDKit✔️✔️:   double cosTheta = r[0].dotProduct(r[1]);
        let mut cos_theta = r[0].dot_product(r[1]);
        // RDKit✔️✔️:   clipToOne(cosTheta);
        clip_to_one(&mut cos_theta);
        // RDKit✔️✔️:   double sinThetaSq = 1.0 - cosTheta * cosTheta;
        let sin_theta_sq = 1.0 - cos_theta * cos_theta;
        // RDKit✔️✔️:   double sinTheta = std::max(sqrt(sinThetaSq), 1.0e-8);
        let sin_theta = sin_theta_sq.sqrt().max(1.0e-8);

        // RDKit✔️✔️:   // use the chain rule:
        // RDKit✔️✔️:   // dE/dx = dE/dTheta * dTheta/dx
        // RDKit✔️✔️:
        // RDKit✔️✔️:   // dE/dTheta is independent of cartesians:
        // RDKit✔️✔️:   double dE_dTheta = getThetaDeriv(cosTheta, sinTheta);
        let mut de_dtheta = self.get_theta_deriv(cos_theta, sin_theta);

        // RDKit✔️✔️:   // The original UFF does not include any penalty for angles that are zero
        // RDKit✔️✔️:   // degrees.
        // RDKit✔️✔️:   //   This can lead to overlapping 1-3 atoms (e.g.e Github #7901), which is
        // RDKit✔️✔️:   //   obviously bad. We add an empiricial penalty for angles close to zero
        // RDKit✔️✔️:   //   borrowed from OpenBabel such that the energy goes up exponentially if
        // RDKit✔️✔️:   //   the angle is less than approx theta0
        // RDKit✔️✔️:   // For the sake of efficiency, we only add the penalty if the angle is less
        // RDKit✔️✔️:   // than 30 degrees
        // RDKit✔️✔️:   if (d_order && d_order < 5 && cosTheta > ANGLE_CORRECTION_THRESHOLD) {
        if self.order != 0 && self.order < 5 && cos_theta > ANGLE_CORRECTION_THRESHOLD {
            // RDKit✔️✔️:     auto theta = acos(cosTheta);
            let theta = cos_theta.acos();

            // RDKit✔️✔️:     auto corr = -20.0 * exp(-20.0 * (theta - d_theta0 + 0.25));
            // RDKit✔️✔️:     dE_dTheta += corr;
            let corr = -20.0 * (-20.0 * (theta - self.theta0 + 0.25)).exp();
            de_dtheta += corr;
            // RDKit✔️✔️:   }
        }

        // RDKit✔️✔️:   Utils::calcAngleBendGrad(r, dist, g, dE_dTheta, cosTheta, sinTheta);
        calc_angle_bend_grad(
            &r,
            &dist,
            grad,
            [self.at1_idx, self.at2_idx, self.at3_idx],
            de_dtheta,
            cos_theta,
            sin_theta,
        );
        // RDKit✔️✔️: }
    }

    #[must_use]
    fn get_energy_term(&self, cos_theta: f64, sin_theta_sq: f64) -> f64 {
        // RDKit✔️✔️: double AngleBendContrib::getEnergyTerm(double cosTheta,
        // RDKit✔️✔️:                                        double sinThetaSq) const {
        // RDKit✔️✔️:   PRECONDITION(d_order == 0 || d_order == 1 || d_order == 2 || d_order == 3 ||
        // RDKit✔️✔️:                    d_order == 4,
        // RDKit✔️✔️:                "bad order");
        assert!(
            self.order == 0 || self.order == 1 || self.order == 2 || self.order == 3 || self.order == 4,
            "bad order"
        );
        // RDKit✔️✔️:   // cos(2x) = cos^2(x) - sin^2(x);
        // RDKit✔️✔️:   double cos2Theta = cosTheta * cosTheta - sinThetaSq;
        let cos2_theta = cos_theta * cos_theta - sin_theta_sq;

        // RDKit✔️✔️:   double res = 0.0;
        // RDKit✔️✔️:   if (d_order == 0) {
        let mut res = if self.order == 0 {
            // RDKit✔️✔️:     res = d_C0 + d_C1 * cosTheta + d_C2 * cos2Theta;
            self.c0 + self.c1 * cos_theta + self.c2 * cos2_theta
            // RDKit✔️✔️:   } else {
        } else {
            // RDKit✔️✔️:     switch (d_order) {
            match self.order {
                // RDKit✔️✔️:       case 1:
                // RDKit✔️✔️:         res = -cosTheta;
                // RDKit✔️✔️:         break;
                1 => -cos_theta,
                // RDKit✔️✔️:       case 2:
                // RDKit✔️✔️:         res = cos2Theta;
                // RDKit✔️✔️:         break;
                2 => cos2_theta,
                // RDKit✔️✔️:       case 3:
                // RDKit✔️✔️:         // cos(3x) = cos^3(x) - 3*cos(x)*sin^2(x)
                // RDKit✔️✔️:         res = cosTheta * (cosTheta * cosTheta - 3. * sinThetaSq);
                // RDKit✔️✔️:         break;
                3 => cos_theta * (cos_theta * cos_theta - 3.0 * sin_theta_sq),
                // RDKit✔️✔️:       case 4:
                // RDKit✔️✔️:         // cos(4x) = cos^4(x) - 6*cos^2(x)*sin^2(x)+sin^4(x)
                // RDKit✔️✔️:         res = int_pow<4>(cosTheta) - 6. * cosTheta * cosTheta * sinThetaSq +
                // RDKit✔️✔️:               sinThetaSq * sinThetaSq;
                // RDKit✔️✔️:         break;
                4 => cos_theta.powi(4) - 6.0 * cos_theta * cos_theta * sin_theta_sq + sin_theta_sq * sin_theta_sq,
                _ => unreachable!("bad order"),
            }
        };
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     res = 1. - res;
        // RDKit✔️✔️:     res /= (double)(d_order * d_order);
        if self.order != 0 {
            res = 1.0 - res;
            res /= f64::from(self.order * self.order);
        }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        res
    }

    #[must_use]
    fn get_theta_deriv(&self, cos_theta: f64, sin_theta: f64) -> f64 {
        // RDKit✔️✔️: double AngleBendContrib::getThetaDeriv(double cosTheta, double sinTheta) const {
        // RDKit✔️✔️:   PRECONDITION(d_order == 0 || d_order == 1 || d_order == 2 || d_order == 3 ||
        // RDKit✔️✔️:                    d_order == 4,
        // RDKit✔️✔️:                "bad order");
        assert!(
            self.order == 0 || self.order == 1 || self.order == 2 || self.order == 3 || self.order == 4,
            "bad order"
        );

        // RDKit✔️✔️:   double dE_dTheta = 0.0;
        // RDKit✔️✔️:   double sin2Theta = 2. * sinTheta * cosTheta;
        let sin2_theta = 2.0 * sin_theta * cos_theta;

        // RDKit✔️✔️:   if (d_order == 0) {
        let de_dtheta = if self.order == 0 {
            // RDKit✔️✔️:     dE_dTheta =
            // RDKit✔️✔️:         -1. * d_forceConstant * (d_C1 * sinTheta + 2. * d_C2 * sin2Theta);
            -self.force_constant * (self.c1 * sin_theta + 2.0 * self.c2 * sin2_theta)
            // RDKit✔️✔️:   } else {
        } else {
            // RDKit✔️✔️:     // E = k/n^2 [1-cos(n theta)]
            // RDKit✔️✔️:     // dE = - k/n^2 * d cos(n theta)
            // RDKit✔️✔️:
            // RDKit✔️✔️:     // these all use:
            // RDKit✔️✔️:     // d cos(ax) = -a sin(ax)
            // RDKit✔️✔️:
            // RDKit✔️✔️:     switch (d_order) {
            let mut value = match self.order {
                // RDKit✔️✔️:       case 1:
                // RDKit✔️✔️:         dE_dTheta = -sinTheta;
                // RDKit✔️✔️:         break;
                1 => -sin_theta,
                // RDKit✔️✔️:       case 2:
                // RDKit✔️✔️:         // sin(2*x) = 2*cos(x)*sin(x)
                // RDKit✔️✔️:         dE_dTheta = sin2Theta;
                // RDKit✔️✔️:         break;
                2 => sin2_theta,
                // RDKit✔️✔️:       case 3:
                // RDKit✔️✔️:         // sin(3*x) = 3*sin(x) - 4*sin^3(x)
                // RDKit✔️✔️:         dE_dTheta = sinTheta * (3. - 4. * sinTheta * sinTheta);
                // RDKit✔️✔️:         break;
                3 => sin_theta * (3.0 - 4.0 * sin_theta * sin_theta),
                // RDKit✔️✔️:       case 4:
                // RDKit✔️✔️:         // sin(4*x) = cos(x)*(4*sin(x) - 8*sin^3(x))
                // RDKit✔️✔️:         dE_dTheta = cosTheta * sinTheta * (4. - 8. * sinTheta * sinTheta);
                // RDKit✔️✔️:         break;
                4 => cos_theta * sin_theta * (4.0 - 8.0 * sin_theta * sin_theta),
                _ => unreachable!("bad order"),
            };
            // RDKit✔️✔️:     }
            // RDKit✔️✔️:     dE_dTheta *= d_forceConstant / (double)(d_order);
            value *= self.force_constant / f64::from(self.order);
            value
            // RDKit✔️✔️:   }
        };
        // RDKit✔️✔️:   return dE_dTheta;
        // RDKit✔️✔️: }
        de_dtheta
    }
}

fn calc_angle_bend_grad(
    r: &[ForceFieldVec3; 2],
    dist: &[f64; 2],
    grad: &mut [f64],
    idx: [usize; 3],
    de_dtheta: f64,
    cos_theta: f64,
    sin_theta: f64,
) {
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

    // RDKit✔️✔️:   g[0][0] += dE_dTheta * dCos_dS[0] / (-sinTheta);
    // RDKit✔️✔️:   g[0][1] += dE_dTheta * dCos_dS[1] / (-sinTheta);
    // RDKit✔️✔️:   g[0][2] += dE_dTheta * dCos_dS[2] / (-sinTheta);
    grad[g0] += de_dtheta * dcos_ds[0] / (-sin_theta);
    grad[g0 + 1] += de_dtheta * dcos_ds[1] / (-sin_theta);
    grad[g0 + 2] += de_dtheta * dcos_ds[2] / (-sin_theta);

    // RDKit✔️✔️:   g[1][0] += dE_dTheta * (-dCos_dS[0] - dCos_dS[3]) / (-sinTheta);
    // RDKit✔️✔️:   g[1][1] += dE_dTheta * (-dCos_dS[1] - dCos_dS[4]) / (-sinTheta);
    // RDKit✔️✔️:   g[1][2] += dE_dTheta * (-dCos_dS[2] - dCos_dS[5]) / (-sinTheta);
    grad[g1] += de_dtheta * (-dcos_ds[0] - dcos_ds[3]) / (-sin_theta);
    grad[g1 + 1] += de_dtheta * (-dcos_ds[1] - dcos_ds[4]) / (-sin_theta);
    grad[g1 + 2] += de_dtheta * (-dcos_ds[2] - dcos_ds[5]) / (-sin_theta);

    // RDKit✔️✔️:   g[2][0] += dE_dTheta * dCos_dS[3] / (-sinTheta);
    // RDKit✔️✔️:   g[2][1] += dE_dTheta * dCos_dS[4] / (-sinTheta);
    // RDKit✔️✔️:   g[2][2] += dE_dTheta * dCos_dS[5] / (-sinTheta);
    // RDKit✔️✔️: }
    grad[g2] += de_dtheta * dcos_ds[3] / (-sin_theta);
    grad[g2 + 1] += de_dtheta * dcos_ds[4] / (-sin_theta);
    grad[g2 + 2] += de_dtheta * dcos_ds[5] / (-sin_theta);
}

impl ForceFieldContrib for AngleBendContrib {
    fn copy(&self) -> Box<dyn ForceFieldContrib> {
        // RDKit✔️✔️: AngleBendContrib *copy() const override {
        // RDKit✔️✔️:   return new AngleBendContrib(*this);
        // RDKit✔️✔️: }
        Box::new(self.clone())
    }

    fn set_force_field(&mut self, owner: *const ForceField) {
        self.owner = owner;
    }

    fn get_energy(&self, pos: &[f64]) -> f64 {
        AngleBendContrib::get_energy(self, pos)
    }

    fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
        AngleBendContrib::get_grad(self, pos, grad);
    }
}

#[cfg(test)]
mod tests {
    use crate::chemistry::forcefield::core::{ForceField, ForceFieldVec3};
    use crate::chemistry::forcefield::uff::utils::calc_angle_force_constant;

    use super::*;

    const EPS: f64 = 1.0e-12;

    fn force_field() -> ForceField {
        let mut ff = ForceField::new(3);
        ff.positions_mut().push(ForceFieldVec3::new(1.0, 0.0, 0.0));
        ff.positions_mut().push(ForceFieldVec3::new(0.0, 0.0, 0.0));
        ff.positions_mut().push(ForceFieldVec3::new(0.0, 1.0, 0.0));
        ff
    }

    fn initialized_force_field() -> ForceField {
        let mut ff = force_field();
        ff.initialize();
        ff
    }

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

    fn params() -> (AtomicParams, AtomicParams, AtomicParams) {
        (
            atomic_params(0.757, 109.47_f64.to_radians(), 5.343, 1.912),
            atomic_params(0.700, 120.0_f64.to_radians(), 6.899, 2.544),
            atomic_params(0.658, 109.47_f64.to_radians(), 8.741, 2.3),
        )
    }

    fn assert_close(actual: f64, expected: f64) {
        assert!((actual - expected).abs() < EPS);
    }

    fn assert_close_tol(actual: f64, expected: f64, tol: f64) {
        assert!(
            (actual - expected).abs() < tol,
            "actual {actual} expected {expected} diff {} tol {tol}",
            (actual - expected).abs()
        );
    }

    fn angle_energy_with_fresh_distance_cache(order: u32, pos: &[f64; 9]) -> f64 {
        let ff = initialized_force_field();
        let (at1, at2, at3) = params();
        let contrib =
            AngleBendContrib::new(&ff, 0, 1, 2, 1.0, 1.0, &at1, &at2, &at3, order).expect("valid constructor");
        contrib.get_energy(pos)
    }

    fn finite_difference_grad(order: u32, pos: &[f64; 9]) -> [f64; 9] {
        let mut out = [0.0; 9];
        let h = 1.0e-6;
        for i in 0..9 {
            let mut plus = *pos;
            let mut minus = *pos;
            plus[i] += h;
            minus[i] -= h;
            out[i] = (angle_energy_with_fresh_distance_cache(order, &plus)
                - angle_energy_with_fresh_distance_cache(order, &minus))
                / (2.0 * h);
        }
        out
    }

    #[test]
    fn uff_anglebendcontrib_constructor_stores_indices_owner_force_and_order0_coefficients() {
        let ff = force_field();
        let (at1, at2, at3) = params();

        let contrib = AngleBendContrib::new(&ff, 0, 1, 2, 1.0, 1.0, &at1, &at2, &at3, 0).expect("valid constructor");
        let expected_force =
            calc_angle_force_constant(at2.theta0, 1.0, 1.0, &at1, &at2, &at3).expect("valid force constant");
        let sin_theta0 = at2.theta0.sin();
        let cos_theta0 = at2.theta0.cos();
        let c2 = 1.0 / (4.0 * (sin_theta0 * sin_theta0).max(1e-8));
        let c1 = -4.0 * c2 * cos_theta0;
        let c0 = c2 * (2.0 * cos_theta0 * cos_theta0 + 1.0);

        assert_eq!(contrib.owner(), &ff as *const ForceField);
        assert_eq!(contrib.at1_idx(), 0);
        assert_eq!(contrib.at2_idx(), 1);
        assert_eq!(contrib.at3_idx(), 2);
        assert_eq!(contrib.order(), 0);
        assert_close(contrib.theta0(), at2.theta0);
        assert_close(contrib.force_constant(), expected_force);
        assert_close(contrib.c0(), c0);
        assert_close(contrib.c1(), c1);
        assert_close(contrib.c2(), c2);
    }

    #[test]
    fn uff_anglebendcontrib_constructor_preserves_nonzero_order_without_coefficients() {
        let ff = force_field();
        let (at1, at2, at3) = params();

        let contrib = AngleBendContrib::new(&ff, 0, 1, 2, 1.0, 1.0, &at1, &at2, &at3, 3).expect("valid constructor");

        assert_eq!(contrib.order(), 3);
        assert_close(contrib.theta0(), at2.theta0);
        assert_eq!(contrib.c0(), 0.0);
        assert_eq!(contrib.c1(), 0.0);
        assert_eq!(contrib.c2(), 0.0);
    }

    #[test]
    fn uff_anglebendcontrib_constructor_applies_ring_hack_orders() {
        let ff = force_field();
        let (at1, at2, at3) = params();
        let cases = [
            (30, 150.0_f64.to_radians()),
            (35, 60.0_f64.to_radians()),
            (40, 135.0_f64.to_radians()),
            (45, 90.0_f64.to_radians()),
        ];

        for (order, theta0) in cases {
            let contrib =
                AngleBendContrib::new(&ff, 0, 1, 2, 1.0, 1.0, &at1, &at2, &at3, order).expect("valid constructor");
            assert_eq!(contrib.order(), 0);
            assert_close(contrib.theta0(), theta0);
            assert!(contrib.c2() > 0.0);
        }
    }

    #[test]
    fn uff_anglebendcontrib_constructor_resets_unrecognized_high_order_without_theta_override() {
        let ff = force_field();
        let (at1, at2, at3) = params();

        let contrib = AngleBendContrib::new(&ff, 0, 1, 2, 1.0, 1.0, &at1, &at2, &at3, 31).expect("valid constructor");

        assert_eq!(contrib.order(), 0);
        assert_close(contrib.theta0(), at2.theta0);
    }

    #[test]
    fn uff_anglebendcontrib_constructor_propagates_invalid_bond_order() {
        let ff = force_field();
        let (at1, at2, at3) = params();

        assert!(matches!(
            AngleBendContrib::new(&ff, 0, 1, 2, 0.0, 1.0, &at1, &at2, &at3, 0),
            Err(UffUtilsError::InvalidBondOrder { bond_order: 0.0 })
        ));
        assert!(matches!(
            AngleBendContrib::new(&ff, 0, 1, 2, 1.0, -1.0, &at1, &at2, &at3, 0),
            Err(UffUtilsError::InvalidBondOrder { bond_order: -1.0 })
        ));
    }

    #[test]
    #[should_panic]
    fn uff_anglebendcontrib_constructor_rejects_degenerate_indices() {
        let ff = force_field();
        let (at1, at2, at3) = params();

        let _ = AngleBendContrib::new(&ff, 0, 0, 2, 1.0, 1.0, &at1, &at2, &at3, 0);
    }

    #[test]
    #[should_panic]
    fn uff_anglebendcontrib_constructor_rejects_index_out_of_range() {
        let ff = force_field();
        let (at1, at2, at3) = params();

        let _ = AngleBendContrib::new(&ff, 0, 1, 3, 1.0, 1.0, &at1, &at2, &at3, 0);
    }

    #[test]
    fn uff_anglebendcontrib_get_energy_order0_uses_polynomial_term() {
        let ff = initialized_force_field();
        let (at1, at2, at3) = params();
        let contrib = AngleBendContrib::new(&ff, 0, 1, 2, 1.0, 1.0, &at1, &at2, &at3, 0).expect("valid constructor");
        let pos = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let cos_theta = 0.0;
        let sin_theta_sq = 1.0;
        let cos2_theta = cos_theta * cos_theta - sin_theta_sq;
        let expected = contrib.force_constant() * (contrib.c0() + contrib.c1() * cos_theta + contrib.c2() * cos2_theta);

        assert_close(contrib.get_energy(&pos), expected);
    }

    #[test]
    fn uff_anglebendcontrib_get_energy_covers_periodic_orders_without_correction() {
        let ff = initialized_force_field();
        let (at1, at2, at3) = params();
        let pos = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let cases = [(1, 1.0), (2, 0.5), (3, 1.0 / 9.0), (4, 0.0)];

        for (order, angle_term) in cases {
            let contrib =
                AngleBendContrib::new(&ff, 0, 1, 2, 1.0, 1.0, &at1, &at2, &at3, order).expect("valid constructor");

            assert_close(contrib.get_energy(&pos), contrib.force_constant() * angle_term);
        }
    }

    #[test]
    fn uff_anglebendcontrib_get_energy_adds_small_angle_correction_only_for_orders_one_to_four() {
        let ff = initialized_force_field();
        let (at1, at2, at3) = params();
        let pos = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.95, 0.2, 0.0];
        let cos_theta = 0.95 / (0.95_f64 * 0.95 + 0.2_f64 * 0.2).sqrt();
        assert!(cos_theta > super::ANGLE_CORRECTION_THRESHOLD);

        let order1 = AngleBendContrib::new(&ff, 0, 1, 2, 1.0, 1.0, &at1, &at2, &at3, 1).expect("valid constructor");
        let order0 = AngleBendContrib::new(&ff, 0, 1, 2, 1.0, 1.0, &at1, &at2, &at3, 0).expect("valid constructor");
        let sin_theta_sq = 1.0 - cos_theta * cos_theta;
        let order1_base = order1.force_constant() * order1.get_energy_term(cos_theta, sin_theta_sq);
        let correction = (-20.0 * (cos_theta.acos() - order1.theta0() + 0.25)).exp();
        let order0_base = order0.force_constant() * order0.get_energy_term(cos_theta, sin_theta_sq);

        assert_close(order1.get_energy(&pos), order1_base + correction);
        assert_close(order0.get_energy(&pos), order0_base);
    }

    #[test]
    fn uff_anglebendcontrib_get_energy_skips_correction_below_threshold() {
        let ff = initialized_force_field();
        let (at1, at2, at3) = params();
        let contrib = AngleBendContrib::new(&ff, 0, 1, 2, 1.0, 1.0, &at1, &at2, &at3, 2).expect("valid constructor");
        let pos = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 0.0];
        let cos_theta = 0.5 / (0.5_f64 * 0.5 + 1.0_f64).sqrt();
        assert!(cos_theta <= super::ANGLE_CORRECTION_THRESHOLD);
        let sin_theta_sq = 1.0 - cos_theta * cos_theta;
        let expected = contrib.force_constant() * contrib.get_energy_term(cos_theta, sin_theta_sq);

        assert_close(contrib.get_energy(&pos), expected);
    }

    #[test]
    fn uff_anglebendcontrib_get_energy_term_clips_source_boundary_inputs_by_caller_contract() {
        let ff = initialized_force_field();
        let (at1, at2, at3) = params();
        let contrib = AngleBendContrib::new(&ff, 0, 1, 2, 1.0, 1.0, &at1, &at2, &at3, 3).expect("valid constructor");
        let mut above = 1.0 + 1.0e-12;
        let mut below = -1.0 - 1.0e-12;

        clip_to_one(&mut above);
        clip_to_one(&mut below);

        assert_eq!(above, 1.0);
        assert_eq!(below, -1.0);
        assert_close(contrib.get_energy_term(above, 0.0), 0.0);
        assert_close(contrib.get_energy_term(below, 0.0), 2.0 / 9.0);
    }

    #[test]
    fn uff_anglebendcontrib_get_energy_order0_constructor_clamps_singular_theta0_coefficients() {
        let ff = initialized_force_field();
        let at1 = atomic_params(0.757, 109.47_f64.to_radians(), 5.343, 1.912);
        let at2 = atomic_params(0.700, 0.0, 6.899, 2.544);
        let at3 = atomic_params(0.658, 109.47_f64.to_radians(), 8.741, 2.3);

        let contrib = AngleBendContrib::new(&ff, 0, 1, 2, 1.0, 1.0, &at1, &at2, &at3, 0).expect("valid constructor");

        assert_close(contrib.c2(), 1.0 / (4.0 * 1.0e-8));
        assert!(
            contrib
                .get_energy(&[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0])
                .is_finite()
        );
    }

    #[test]
    fn uff_anglebendcontrib_get_grad_matches_energy_finite_difference_for_order0() {
        let ff = initialized_force_field();
        let (at1, at2, at3) = params();
        let contrib = AngleBendContrib::new(&ff, 0, 1, 2, 1.0, 1.0, &at1, &at2, &at3, 0).expect("valid constructor");
        let pos = [1.1, 0.2, 0.0, 0.0, 0.0, 0.1, -0.1, 0.9, 0.2];
        let mut grad = [0.0; 9];

        contrib.get_grad(&pos, &mut grad);
        let expected = finite_difference_grad(0, &pos);

        for i in 0..9 {
            assert_close_tol(grad[i], expected[i], 1.0e-5);
        }
        for axis in 0..3 {
            assert_close_tol(grad[axis] + grad[3 + axis] + grad[6 + axis], 0.0, 1.0e-10);
        }
    }

    #[test]
    fn uff_anglebendcontrib_get_grad_matches_energy_finite_difference_for_periodic_order() {
        let ff = initialized_force_field();
        let (at1, at2, at3) = params();
        let contrib = AngleBendContrib::new(&ff, 0, 1, 2, 1.0, 1.0, &at1, &at2, &at3, 3).expect("valid constructor");
        let pos = [1.2, -0.2, 0.1, 0.0, 0.0, 0.0, -0.2, 1.1, -0.1];
        let mut grad = [0.0; 9];

        contrib.get_grad(&pos, &mut grad);
        let expected = finite_difference_grad(3, &pos);

        for i in 0..9 {
            assert_close_tol(grad[i], expected[i], 1.0e-5);
        }
    }

    #[test]
    fn uff_anglebendcontrib_get_grad_includes_small_angle_correction_derivative() {
        let ff = initialized_force_field();
        let (at1, at2, at3) = params();
        let contrib = AngleBendContrib::new(&ff, 0, 1, 2, 1.0, 1.0, &at1, &at2, &at3, 1).expect("valid constructor");
        let pos = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.95, 0.2, 0.0];
        let mut grad = [0.0; 9];

        contrib.get_grad(&pos, &mut grad);
        let p1 = ForceFieldVec3::new(pos[0], pos[1], pos[2]);
        let p2 = ForceFieldVec3::new(pos[3], pos[4], pos[5]);
        let p3 = ForceFieldVec3::new(pos[6], pos[7], pos[8]);
        let dist = [(p1 - p2).length(), (p3 - p2).length()];
        let r = [(p1 - p2) / dist[0], (p3 - p2) / dist[1]];
        let mut cos_theta = r[0].dot_product(r[1]);
        clip_to_one(&mut cos_theta);
        assert!(cos_theta > super::ANGLE_CORRECTION_THRESHOLD);
        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt().max(1.0e-8);
        let theta = cos_theta.acos();
        let base_de_dtheta = contrib.get_theta_deriv(cos_theta, sin_theta);
        let corrected_de_dtheta = base_de_dtheta + -20.0 * (-20.0 * (theta - contrib.theta0() + 0.25)).exp();
        let mut base_grad = [0.0; 9];
        let mut expected = [0.0; 9];

        calc_angle_bend_grad(
            &r,
            &dist,
            &mut base_grad,
            [0, 1, 2],
            base_de_dtheta,
            cos_theta,
            sin_theta,
        );
        calc_angle_bend_grad(
            &r,
            &dist,
            &mut expected,
            [0, 1, 2],
            corrected_de_dtheta,
            cos_theta,
            sin_theta,
        );

        for i in 0..9 {
            assert_close_tol(grad[i], expected[i], 1.0e-8);
        }
        assert!(
            grad.iter()
                .zip(base_grad.iter())
                .any(|(actual, base)| (actual - base).abs() > 1.0e-8)
        );
    }

    #[test]
    fn uff_anglebendcontrib_calc_angle_bend_grad_preserves_source_evaluation_order() {
        let r = [
            ForceFieldVec3::new(0.6, 0.8, -0.1),
            ForceFieldVec3::new(7.516_508_172_266_185, 0.5, 0.84),
        ];
        let dist = [1.0, 1.75];
        let idx = [2, 0, 3];
        let de_dtheta = -5.714_193_060_994_148;
        let cos_theta = 0.0;
        let sin_theta = 0.486_800_828_948_617;
        let mut grad = [0.0; 12];
        let mut expected = grad;
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
        let mut factorized = grad;
        let scale = de_dtheta / (-sin_theta);

        expected[g0] += de_dtheta * dcos_ds[0] / (-sin_theta);
        expected[g0 + 1] += de_dtheta * dcos_ds[1] / (-sin_theta);
        expected[g0 + 2] += de_dtheta * dcos_ds[2] / (-sin_theta);
        expected[g1] += de_dtheta * (-dcos_ds[0] - dcos_ds[3]) / (-sin_theta);
        expected[g1 + 1] += de_dtheta * (-dcos_ds[1] - dcos_ds[4]) / (-sin_theta);
        expected[g1 + 2] += de_dtheta * (-dcos_ds[2] - dcos_ds[5]) / (-sin_theta);
        expected[g2] += de_dtheta * dcos_ds[3] / (-sin_theta);
        expected[g2 + 1] += de_dtheta * dcos_ds[4] / (-sin_theta);
        expected[g2 + 2] += de_dtheta * dcos_ds[5] / (-sin_theta);

        factorized[g0] += scale * dcos_ds[0];
        factorized[g0 + 1] += scale * dcos_ds[1];
        factorized[g0 + 2] += scale * dcos_ds[2];
        factorized[g1] += scale * (-dcos_ds[0] - dcos_ds[3]);
        factorized[g1 + 1] += scale * (-dcos_ds[1] - dcos_ds[4]);
        factorized[g1 + 2] += scale * (-dcos_ds[2] - dcos_ds[5]);
        factorized[g2] += scale * dcos_ds[3];
        factorized[g2 + 1] += scale * dcos_ds[4];
        factorized[g2 + 2] += scale * dcos_ds[5];

        calc_angle_bend_grad(&r, &dist, &mut grad, idx, de_dtheta, cos_theta, sin_theta);

        for i in 0..grad.len() {
            assert_eq!(grad[i].to_bits(), expected[i].to_bits(), "flat axis {i}");
        }
        assert!(
            expected
                .iter()
                .zip(factorized.iter())
                .any(|(source, optimized)| source.to_bits() != optimized.to_bits()),
            "boundary input must distinguish RDKit's multiply-then-divide order"
        );
    }

    #[test]
    fn uff_anglebendcontrib_calc_angle_bend_grad_uses_caller_singularity_clamp_input() {
        let r = [
            ForceFieldVec3::new(1.0, 0.0, 0.0),
            ForceFieldVec3::new(1.0, 1.0e-8, 0.0),
        ];
        let dist = [1.0, 1.0];
        let mut grad = [0.0; 9];

        calc_angle_bend_grad(&r, &dist, &mut grad, [0, 1, 2], 1.0e-8, 1.0, 1.0e-8);

        assert!(grad.iter().all(|value| value.is_finite()));
        assert_close_tol(grad[1], -1.0e-8, 1.0e-20);
        assert_close_tol(grad[4], 0.0, 1.0e-20);
        assert_close_tol(grad[7], 1.0e-8, 1.0e-20);
    }

    #[test]
    fn uff_anglebendcontrib_get_grad_theta_deriv_covers_all_order_branches() {
        let ff = initialized_force_field();
        let (at1, at2, at3) = params();
        let cos_theta = 0.25_f64;
        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
        let sin2_theta = 2.0 * sin_theta * cos_theta;

        let order0 = AngleBendContrib::new(&ff, 0, 1, 2, 1.0, 1.0, &at1, &at2, &at3, 0).expect("valid constructor");
        assert_close(
            order0.get_theta_deriv(cos_theta, sin_theta),
            -order0.force_constant() * (order0.c1() * sin_theta + 2.0 * order0.c2() * sin2_theta),
        );

        for order in 1..=4 {
            let contrib =
                AngleBendContrib::new(&ff, 0, 1, 2, 1.0, 1.0, &at1, &at2, &at3, order).expect("valid constructor");
            let base = match order {
                1 => -sin_theta,
                2 => sin2_theta,
                3 => sin_theta * (3.0 - 4.0 * sin_theta * sin_theta),
                4 => cos_theta * sin_theta * (4.0 - 8.0 * sin_theta * sin_theta),
                _ => unreachable!(),
            };
            assert_close(
                contrib.get_theta_deriv(cos_theta, sin_theta),
                base * contrib.force_constant() / f64::from(order),
            );
        }
    }

    #[test]
    fn uff_anglebendcontrib_get_grad_singularity_clamp_keeps_collinear_gradient_finite() {
        let ff = initialized_force_field();
        let (at1, at2, at3) = params();
        let contrib = AngleBendContrib::new(&ff, 0, 1, 2, 1.0, 1.0, &at1, &at2, &at3, 2).expect("valid constructor");
        let pos = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0];
        let mut grad = [0.0; 9];

        contrib.get_grad(&pos, &mut grad);

        assert!(grad.iter().all(|value| value.is_finite()));
    }
}
