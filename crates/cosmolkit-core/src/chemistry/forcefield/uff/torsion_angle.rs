//! Source-backed RDKit UFF torsion angle contribution.

use crate::Hybridization;
use crate::chemistry::forcefield::core::{ForceField, ForceFieldContrib, ForceFieldVec3, compute_dihedral_from_flat};

use super::params::{AtomicParams, clip_to_one, is_double_zero};

#[derive(Clone, Debug)]
pub struct TorsionAngleContrib {
    owner: *const ForceField,
    at1_idx: usize,
    at2_idx: usize,
    at3_idx: usize,
    at4_idx: usize,
    order: u32,
    force_constant: f64,
    cos_term: f64,
}

impl TorsionAngleContrib {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        owner: &ForceField,
        idx1: usize,
        idx2: usize,
        idx3: usize,
        idx4: usize,
        bond_order23: f64,
        at_num2: i32,
        at_num3: i32,
        hyb2: Hybridization,
        hyb3: Hybridization,
        at2_params: &AtomicParams,
        at3_params: &AtomicParams,
        end_atom_is_sp2: bool,
    ) -> Self {
        // RDKit✔️✔️: TorsionAngleContrib::TorsionAngleContrib(
        // RDKit✔️✔️:     ForceField *owner, unsigned int idx1, unsigned int idx2, unsigned int idx3,
        // RDKit✔️✔️:     unsigned int idx4, double bondOrder23, int atNum2, int atNum3,
        // RDKit✔️✔️:     RDKit::Atom::HybridizationType hyb2, RDKit::Atom::HybridizationType hyb3,
        // RDKit✔️✔️:     const AtomicParams *at2Params, const AtomicParams *at3Params,
        // RDKit✔️✔️:     bool endAtomIsSP2) {
        // RDKit✔️✔️:   PRECONDITION(owner, "bad owner");
        // RDKit✔️✔️:   PRECONDITION(at2Params, "bad params pointer");
        // RDKit✔️✔️:   PRECONDITION(at3Params, "bad params pointer");
        // Rust references reproduce RDKit's non-null preconditions.
        // RDKit✔️✔️:   PRECONDITION((idx1 != idx2 && idx1 != idx3 && idx1 != idx4 && idx2 != idx3 &&
        // RDKit✔️✔️:                 idx2 != idx4 && idx3 != idx4),
        // RDKit✔️✔️:                "degenerate points");
        assert!(idx1 != idx2 && idx1 != idx3 && idx1 != idx4 && idx2 != idx3 && idx2 != idx4 && idx3 != idx4);
        // RDKit✔️✔️:   URANGE_CHECK(idx1, owner->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx2, owner->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx3, owner->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx4, owner->positions().size());
        assert!(idx1 < owner.positions().len());
        assert!(idx2 < owner.positions().len());
        assert!(idx3 < owner.positions().len());
        assert!(idx4 < owner.positions().len());

        // RDKit✔️✔️:   dp_forceField = owner;
        // RDKit✔️✔️:   d_at1Idx = idx1;
        // RDKit✔️✔️:   d_at2Idx = idx2;
        // RDKit✔️✔️:   d_at3Idx = idx3;
        // RDKit✔️✔️:   d_at4Idx = idx4;
        let owner_ptr = owner as *const ForceField;

        // RDKit✔️✔️:   calcTorsionParams(bondOrder23, atNum2, atNum3, hyb2, hyb3, at2Params,
        // RDKit✔️✔️:                     at3Params, endAtomIsSP2);
        // RDKit✔️✔️: }
        let (force_constant, order, cos_term) = calc_torsion_params(
            bond_order23,
            at_num2,
            at_num3,
            hyb2,
            hyb3,
            at2_params,
            at3_params,
            end_atom_is_sp2,
        );

        Self {
            owner: owner_ptr,
            at1_idx: idx1,
            at2_idx: idx2,
            at3_idx: idx3,
            at4_idx: idx4,
            order,
            force_constant,
            cos_term,
        }
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
    pub fn at4_idx(&self) -> usize {
        self.at4_idx
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
    pub fn cos_term(&self) -> f64 {
        self.cos_term
    }

    pub(crate) fn scale_force_constant(&mut self, count: usize) {
        // RDKit✔️✔️: void scaleForceConstant(unsigned int count) {
        // RDKit✔️✔️:   this->d_forceConstant /= static_cast<double>(count);
        // RDKit✔️✔️: }
        self.force_constant /= count as f64;
    }

    #[must_use]
    pub fn get_energy(&self, pos: &[f64]) -> f64 {
        // RDKit✔️✔️: double TorsionAngleContrib::getEnergy(double *pos) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        assert!(!self.owner.is_null(), "no owner");
        // RDKit✔️✔️:   PRECONDITION(d_order == 2 || d_order == 3 || d_order == 6, "bad order");
        assert!(self.order == 2 || self.order == 3 || self.order == 6);

        // RDKit✔️✔️:   RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
        // RDKit✔️✔️:                      pos[3 * d_at1Idx + 2]);
        // RDKit✔️✔️:   RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
        // RDKit✔️✔️:                      pos[3 * d_at2Idx + 2]);
        // RDKit✔️✔️:   RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
        // RDKit✔️✔️:                      pos[3 * d_at3Idx + 2]);
        // RDKit✔️✔️:   RDGeom::Point3D p4(pos[3 * d_at4Idx], pos[3 * d_at4Idx + 1],
        // RDKit✔️✔️:                      pos[3 * d_at4Idx + 2]);
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
        let p4 = ForceFieldVec3::new(
            pos[3 * self.at4_idx],
            pos[3 * self.at4_idx + 1],
            pos[3 * self.at4_idx + 2],
        );

        // RDKit✔️✔️:   double cosPhi = Utils::calculateCosTorsion(p1, p2, p3, p4);
        let cos_phi = calculate_cos_torsion(p1, p2, p3, p4);
        // RDKit✔️✔️:   double sinPhiSq = 1 - cosPhi * cosPhi;
        let sin_phi_sq = 1.0 - cos_phi * cos_phi;

        // RDKit✔️✔️:   // E(phi) = V/2 * (1 - cos(n*phi_0)*cos(n*phi))
        // RDKit✔️✔️:   double cosNPhi = 0.0;
        // RDKit✔️✔️:   switch (d_order) {
        let cos_n_phi = match self.order {
            // RDKit✔️✔️:     case 2:
            // RDKit✔️✔️:       // cos(2x) = 1 - 2sin^2(x)
            // RDKit✔️✔️:       cosNPhi = 1 - 2 * sinPhiSq;
            // RDKit✔️✔️:       break;
            2 => 1.0 - 2.0 * sin_phi_sq,
            // RDKit✔️✔️:     case 3:
            // RDKit✔️✔️:       // cos(3x) = cos^3(x) - 3*cos(x)*sin^2(x) = 4cos^3(x) -3cos(x)
            // RDKit✔️✔️:       cosNPhi = cosPhi * (cosPhi * cosPhi - 3. * sinPhiSq);
            // RDKit✔️✔️:       break;
            3 => cos_phi * (cos_phi * cos_phi - 3.0 * sin_phi_sq),
            // RDKit✔️✔️:     case 6:
            // RDKit✔️✔️:       // cos(6x) = 1 - 32*sin^6(x) + 48*sin^4(x) - 18*sin^2(x)
            // RDKit✔️✔️:       cosNPhi =
            // RDKit✔️✔️:           1 + sinPhiSq * (-32. * sinPhiSq * sinPhiSq + 48. * sinPhiSq - 18.);
            // RDKit✔️✔️:       break;
            6 => 1.0 + sin_phi_sq * (-32.0 * sin_phi_sq * sin_phi_sq + 48.0 * sin_phi_sq - 18.0),
            _ => unreachable!("bad order"),
        };
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   double res = d_forceConstant / 2.0 * (1. - d_cosTerm * cosNPhi);
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        self.force_constant / 2.0 * (1.0 - self.cos_term * cos_n_phi)
    }

    pub fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
        // RDKit✔️✔️: void TorsionAngleContrib::getGrad(double *pos, double *grad) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // RDKit✔️✔️:   PRECONDITION(grad, "bad vector");
        assert!(!self.owner.is_null(), "no owner");
        assert!(!pos.is_empty(), "bad vector");
        assert!(!grad.is_empty(), "bad vector");

        // RDKit✔️✔️:   double *g[4] = {&(grad[3 * d_at1Idx]), &(grad[3 * d_at2Idx]),
        // RDKit✔️✔️:                   &(grad[3 * d_at3Idx]), &(grad[3 * d_at4Idx])};
        // Rust passes the gradient slice and atom indices to calc_torsion_grad.

        // RDKit✔️✔️:   RDGeom::Point3D r[4];
        // RDKit✔️✔️:   RDGeom::Point3D t[2];
        // RDKit✔️✔️:   double d[2];
        // RDKit✔️✔️:   double cosPhi;
        // RDKit✔️✔️:   RDKit::ForceFieldsHelper::computeDihedral(
        // RDKit✔️✔️:       pos, d_at1Idx, d_at2Idx, d_at3Idx, d_at4Idx, nullptr, &cosPhi, r, t, d);
        let dihedral_output =
            compute_dihedral_from_flat(pos, self.at1_idx, self.at2_idx, self.at3_idx, self.at4_idx, false);
        let r = dihedral_output.r;
        let t = dihedral_output.t;
        let d = dihedral_output.d;
        let cos_phi = dihedral_output.cos_phi;

        // RDKit✔️✔️:   double sinPhiSq = 1.0 - cosPhi * cosPhi;
        let sin_phi_sq = 1.0 - cos_phi * cos_phi;
        // RDKit✔️✔️:   double sinPhi = ((sinPhiSq > 0.0) ? sqrt(sinPhiSq) : 0.0);
        let sin_phi = if sin_phi_sq > 0.0 { sin_phi_sq.sqrt() } else { 0.0 };

        // RDKit✔️✔️:   // dE/dPhi is independent of cartesians:
        // RDKit✔️✔️:   double dE_dPhi = getThetaDeriv(cosPhi, sinPhi);
        let de_dphi = self.get_theta_deriv(cos_phi, sin_phi);
        // RDKit✔️✔️:   double sinTerm =
        // RDKit✔️✔️:       dE_dPhi * (isDoubleZero(sinPhi) ? (1.0 / cosPhi) : (1.0 / sinPhi));
        let sin_term = de_dphi
            * if is_double_zero(sin_phi) {
                1.0 / cos_phi
            } else {
                1.0 / sin_phi
            };

        // RDKit✔️✔️:   Utils::calcTorsionGrad(r, t, d, g, sinTerm, cosPhi);
        calc_torsion_grad(
            &r,
            &t,
            &d,
            grad,
            [self.at1_idx, self.at2_idx, self.at3_idx, self.at4_idx],
            sin_term,
            cos_phi,
        );
        // RDKit✔️✔️: }
    }

    #[must_use]
    fn get_theta_deriv(&self, cos_theta: f64, sin_theta: f64) -> f64 {
        // RDKit✔️✔️: double TorsionAngleContrib::getThetaDeriv(double cosTheta,
        // RDKit✔️✔️:                                           double sinTheta) const {
        // RDKit✔️✔️:   PRECONDITION(d_order == 2 || d_order == 3 || d_order == 6, "bad order");
        assert!(self.order == 2 || self.order == 3 || self.order == 6);
        // RDKit✔️✔️:   double sinThetaSq = sinTheta * sinTheta;
        let sin_theta_sq = sin_theta * sin_theta;
        // RDKit✔️✔️:   // cos(6x) = 1 - 32*sin^6(x) + 48*sin^4(x) - 18*sin^2(x)

        // RDKit✔️✔️:   double res = 0.0;
        // RDKit✔️✔️:   switch (d_order) {
        let mut res = match self.order {
            // RDKit✔️✔️:     case 2:
            // RDKit✔️✔️:       res = 2 * sinTheta * cosTheta;
            // RDKit✔️✔️:       break;
            2 => 2.0 * sin_theta * cos_theta,
            // RDKit✔️✔️:     case 3:
            // RDKit✔️✔️:       // sin(3*x) = 3*sin(x) - 4*sin^3(x)
            // RDKit✔️✔️:       res = sinTheta * (3 - 4 * sinThetaSq);
            // RDKit✔️✔️:       break;
            3 => sin_theta * (3.0 - 4.0 * sin_theta_sq),
            // RDKit✔️✔️:     case 6:
            // RDKit✔️✔️:       // sin(6x) = cos(x) * [ 32*sin^5(x) - 32*sin^3(x) + 6*sin(x) ]
            // RDKit✔️✔️:       res = cosTheta * sinTheta * (32 * sinThetaSq * (sinThetaSq - 1) + 6);
            // RDKit✔️✔️:       break;
            6 => cos_theta * sin_theta * (32.0 * sin_theta_sq * (sin_theta_sq - 1.0) + 6.0),
            _ => unreachable!("bad order"),
        };
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   res *= d_forceConstant / 2.0 * d_cosTerm * -1 * d_order;
        res *= self.force_constant / 2.0 * self.cos_term * -1.0 * f64::from(self.order);

        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        res
    }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn calc_torsion_params(
    bond_order23: f64,
    at_num2: i32,
    at_num3: i32,
    hyb2: Hybridization,
    hyb3: Hybridization,
    at2_params: &AtomicParams,
    at3_params: &AtomicParams,
    end_atom_is_sp2: bool,
) -> (f64, u32, f64) {
    // RDKit✔️✔️: void TorsionAngleContrib::calcTorsionParams(double bondOrder23, int atNum2,
    // RDKit✔️✔️:                                             int atNum3,
    // RDKit✔️✔️:                                             RDKit::Atom::HybridizationType hyb2,
    // RDKit✔️✔️:                                             RDKit::Atom::HybridizationType hyb3,
    // RDKit✔️✔️:                                             const AtomicParams *at2Params,
    // RDKit✔️✔️:                                             const AtomicParams *at3Params,
    // RDKit✔️✔️:                                             bool endAtomIsSP2) {
    // RDKit✔️✔️:   PRECONDITION((hyb2 == RDKit::Atom::SP2 || hyb2 == RDKit::Atom::SP3) &&
    // RDKit✔️✔️:                    (hyb3 == RDKit::Atom::SP2 || hyb3 == RDKit::Atom::SP3),
    // RDKit✔️✔️:                "bad hybridizations");
    assert!(
        (hyb2 == Hybridization::Sp2 || hyb2 == Hybridization::Sp3)
            && (hyb3 == Hybridization::Sp2 || hyb3 == Hybridization::Sp3)
    );

    let force_constant;
    let order;
    let cos_term;

    // RDKit✔️✔️:   if (hyb2 == RDKit::Atom::SP3 && hyb3 == RDKit::Atom::SP3) {
    if hyb2 == Hybridization::Sp3 && hyb3 == Hybridization::Sp3 {
        // RDKit✔️✔️:     // general case:
        // RDKit✔️✔️:     d_forceConstant = sqrt(at2Params->V1 * at3Params->V1);
        force_constant = (at2_params.v1 * at3_params.v1).sqrt();
        // RDKit✔️✔️:     d_order = 3;
        order = 3;
        // RDKit✔️✔️:     d_cosTerm = -1;  // phi0=60
        cos_term = -1.0;

        // RDKit✔️✔️:     // special case for single bonds between group 6 elements:
        // RDKit✔️✔️:     if (bondOrder23 == 1.0 && Utils::isInGroup6(atNum2) &&
        // RDKit✔️✔️:         Utils::isInGroup6(atNum3)) {
        if bond_order23 == 1.0 && is_in_group6(at_num2) && is_in_group6(at_num3) {
            // RDKit✔️✔️:       double V2 = 6.8, V3 = 6.8;
            let mut v2: f64 = 6.8;
            let mut v3: f64 = 6.8;
            // RDKit✔️✔️:       if (atNum2 == 8) {
            // RDKit✔️✔️:         V2 = 2.0;
            // RDKit✔️✔️:       }
            if at_num2 == 8 {
                v2 = 2.0;
            }
            // RDKit✔️✔️:       if (atNum3 == 8) {
            // RDKit✔️✔️:         V3 = 2.0;
            // RDKit✔️✔️:       }
            if at_num3 == 8 {
                v3 = 2.0;
            }
            // RDKit✔️✔️:       d_forceConstant = sqrt(V2 * V3);
            // RDKit✔️✔️:       d_order = 2;
            // RDKit✔️✔️:       d_cosTerm = -1;  // phi0=90
            // RDKit✔️✔️:     }
            return ((v2 * v3).sqrt(), 2, -1.0);
        }
        // RDKit✔️✔️:   } else if (hyb2 == RDKit::Atom::SP2 && hyb3 == RDKit::Atom::SP2) {
    } else if hyb2 == Hybridization::Sp2 && hyb3 == Hybridization::Sp2 {
        // RDKit✔️✔️:     d_forceConstant = Utils::equation17(bondOrder23, at2Params, at3Params);
        force_constant = equation17(bond_order23, at2_params, at3_params);
        // RDKit✔️✔️:     d_order = 2;
        order = 2;
        // RDKit✔️✔️:     // FIX: is this angle term right?
        // RDKit✔️✔️:     d_cosTerm = 1.0;  // phi0= 180
        cos_term = 1.0;
        // RDKit✔️✔️:   } else {
    } else {
        // RDKit✔️✔️:     // SP2 - SP3,  this is, by default, independent of atom type in UFF:
        // RDKit✔️✔️:     d_forceConstant = 1.0;
        force_constant = 1.0;
        // RDKit✔️✔️:     d_order = 6;
        order = 6;
        // RDKit✔️✔️:     d_cosTerm = 1.0;  // phi0 = 0
        cos_term = 1.0;
        // RDKit✔️✔️:     if (bondOrder23 == 1.0) {
        if bond_order23 == 1.0 {
            // RDKit✔️✔️:       // special case between group 6 sp3 and non-group 6 sp2:
            // RDKit✔️✔️:       if ((hyb2 == RDKit::Atom::SP3 && Utils::isInGroup6(atNum2) &&
            // RDKit✔️✔️:            !Utils::isInGroup6(atNum3)) ||
            // RDKit✔️✔️:           (hyb3 == RDKit::Atom::SP3 && Utils::isInGroup6(atNum3) &&
            // RDKit✔️✔️:            !Utils::isInGroup6(atNum2))) {
            if (hyb2 == Hybridization::Sp3 && is_in_group6(at_num2) && !is_in_group6(at_num3))
                || (hyb3 == Hybridization::Sp3 && is_in_group6(at_num3) && !is_in_group6(at_num2))
            {
                // RDKit✔️✔️:         d_forceConstant = Utils::equation17(bondOrder23, at2Params, at3Params);
                // RDKit✔️✔️:         d_order = 2;
                // RDKit✔️✔️:         d_cosTerm = -1;  // phi0 = 90;
                return (equation17(bond_order23, at2_params, at3_params), 2, -1.0);
            }

            // RDKit✔️✔️:       // special case for sp3 - sp2 - sp2
            // RDKit✔️✔️:       // (i.e. the sp2 has another sp2 neighbor, like propene)
            // RDKit✔️✔️:       else if (endAtomIsSP2) {
            if end_atom_is_sp2 {
                // RDKit✔️✔️:         d_forceConstant = 2.0;
                // RDKit✔️✔️:         d_order = 3;
                // RDKit✔️✔️:         d_cosTerm = -1;  // phi0 = 180;
                // RDKit✔️✔️:       }
                return (2.0, 3, -1.0);
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
    }

    (force_constant, order, cos_term)
}

pub(crate) fn is_in_group6(num: i32) -> bool {
    // RDKit✔️✔️: bool isInGroup6(int num) {
    // RDKit✔️✔️:   return (num == 8 || num == 16 || num == 34 || num == 52 || num == 84);
    // RDKit✔️✔️: }
    matches!(num, 8 | 16 | 34 | 52 | 84)
}

fn calculate_cos_torsion(p1: ForceFieldVec3, p2: ForceFieldVec3, p3: ForceFieldVec3, p4: ForceFieldVec3) -> f64 {
    // RDKit✔️✔️: double calculateCosTorsion(const RDGeom::Point3D &p1, const RDGeom::Point3D &p2,
    // RDKit✔️✔️:                            const RDGeom::Point3D &p3,
    // RDKit✔️✔️:                            const RDGeom::Point3D &p4) {
    // RDKit✔️✔️:   RDGeom::Point3D r1 = p1 - p2, r2 = p3 - p2, r3 = p2 - p3, r4 = p4 - p3;
    let r1 = p1 - p2;
    let r2 = p3 - p2;
    let r3 = p2 - p3;
    let r4 = p4 - p3;
    // RDKit✔️✔️:   RDGeom::Point3D t1 = r1.crossProduct(r2);
    // RDKit✔️✔️:   RDGeom::Point3D t2 = r3.crossProduct(r4);
    let t1 = r1.cross_product(r2);
    let t2 = r3.cross_product(r4);
    // RDKit✔️✔️:   double d1 = t1.length();
    // RDKit✔️✔️:   double d2 = t2.length();
    let d1 = t1.length();
    let d2 = t2.length();
    // RDKit✔️✔️:   if (isDoubleZero(d1) || isDoubleZero(d2)) {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   }
    if is_double_zero(d1) || is_double_zero(d2) {
        return 0.0;
    }
    // RDKit✔️✔️:   double cosPhi = t1.dotProduct(t2) / (d1 * d2);
    let mut cos_phi = t1.dot_product(t2) / (d1 * d2);
    // RDKit✔️✔️:   clipToOne(cosPhi);
    clip_to_one(&mut cos_phi);
    // RDKit✔️✔️:   return cosPhi;
    // RDKit✔️✔️: }
    cos_phi
}

pub(crate) fn equation17(bond_order23: f64, at2_params: &AtomicParams, at3_params: &AtomicParams) -> f64 {
    // RDKit✔️✔️: double equation17(double bondOrder23, const AtomicParams *at2Params,
    // RDKit✔️✔️:                   const AtomicParams *at3Params) {
    // RDKit✔️✔️:   return 5. * sqrt(at2Params->U1 * at3Params->U1) *
    // RDKit✔️✔️:          (1. + 4.18 * log(bondOrder23));
    // RDKit✔️✔️: }
    5.0 * (at2_params.u1 * at3_params.u1).sqrt() * (1.0 + 4.18 * bond_order23.ln())
}

fn calc_torsion_grad(
    r: &[ForceFieldVec3; 4],
    t: &[ForceFieldVec3; 2],
    d: &[f64; 2],
    grad: &mut [f64],
    idx: [usize; 4],
    sin_term: f64,
    cos_phi: f64,
) {
    // RDKit✔️✔️: void calcTorsionGrad(RDGeom::Point3D *r, RDGeom::Point3D *t, double *d,
    // RDKit✔️✔️:                      double **g, double &sinTerm, double &cosPhi) {
    // RDKit✔️✔️:   // -------
    // RDKit✔️✔️:   // dTheta/dx is trickier:
    // RDKit✔️✔️:   double dCos_dT[6] = {1.0 / d[0] * (t[1].x - cosPhi * t[0].x),
    // RDKit✔️✔️:                        1.0 / d[0] * (t[1].y - cosPhi * t[0].y),
    // RDKit✔️✔️:                        1.0 / d[0] * (t[1].z - cosPhi * t[0].z),
    // RDKit✔️✔️:                        1.0 / d[1] * (t[0].x - cosPhi * t[1].x),
    // RDKit✔️✔️:                        1.0 / d[1] * (t[0].y - cosPhi * t[1].y),
    // RDKit✔️✔️:                        1.0 / d[1] * (t[0].z - cosPhi * t[1].z)};
    let dcos_dt = [
        1.0 / d[0] * (t[1].x - cos_phi * t[0].x),
        1.0 / d[0] * (t[1].y - cos_phi * t[0].y),
        1.0 / d[0] * (t[1].z - cos_phi * t[0].z),
        1.0 / d[1] * (t[0].x - cos_phi * t[1].x),
        1.0 / d[1] * (t[0].y - cos_phi * t[1].y),
        1.0 / d[1] * (t[0].z - cos_phi * t[1].z),
    ];

    let g0 = 3 * idx[0];
    let g1 = 3 * idx[1];
    let g2 = 3 * idx[2];
    let g3 = 3 * idx[3];

    // RDKit✔️✔️:   g[0][0] += sinTerm * (dCos_dT[2] * r[1].y - dCos_dT[1] * r[1].z);
    // RDKit✔️✔️:   g[0][1] += sinTerm * (dCos_dT[0] * r[1].z - dCos_dT[2] * r[1].x);
    // RDKit✔️✔️:   g[0][2] += sinTerm * (dCos_dT[1] * r[1].x - dCos_dT[0] * r[1].y);
    grad[g0] += sin_term * (dcos_dt[2] * r[1].y - dcos_dt[1] * r[1].z);
    grad[g0 + 1] += sin_term * (dcos_dt[0] * r[1].z - dcos_dt[2] * r[1].x);
    grad[g0 + 2] += sin_term * (dcos_dt[1] * r[1].x - dcos_dt[0] * r[1].y);

    // RDKit✔️✔️:   g[1][0] += sinTerm *
    // RDKit✔️✔️:              (dCos_dT[1] * (r[1].z - r[0].z) + dCos_dT[2] * (r[0].y - r[1].y) +
    // RDKit✔️✔️:               dCos_dT[4] * (-r[3].z) + dCos_dT[5] * (r[3].y));
    // RDKit✔️✔️:   g[1][1] += sinTerm *
    // RDKit✔️✔️:              (dCos_dT[0] * (r[0].z - r[1].z) + dCos_dT[2] * (r[1].x - r[0].x) +
    // RDKit✔️✔️:               dCos_dT[3] * (r[3].z) + dCos_dT[5] * (-r[3].x));
    // RDKit✔️✔️:   g[1][2] += sinTerm *
    // RDKit✔️✔️:              (dCos_dT[0] * (r[1].y - r[0].y) + dCos_dT[1] * (r[0].x - r[1].x) +
    // RDKit✔️✔️:               dCos_dT[3] * (-r[3].y) + dCos_dT[4] * (r[3].x));
    grad[g1] += sin_term
        * (dcos_dt[1] * (r[1].z - r[0].z)
            + dcos_dt[2] * (r[0].y - r[1].y)
            + dcos_dt[4] * (-r[3].z)
            + dcos_dt[5] * r[3].y);
    grad[g1 + 1] += sin_term
        * (dcos_dt[0] * (r[0].z - r[1].z)
            + dcos_dt[2] * (r[1].x - r[0].x)
            + dcos_dt[3] * r[3].z
            + dcos_dt[5] * (-r[3].x));
    grad[g1 + 2] += sin_term
        * (dcos_dt[0] * (r[1].y - r[0].y)
            + dcos_dt[1] * (r[0].x - r[1].x)
            + dcos_dt[3] * (-r[3].y)
            + dcos_dt[4] * r[3].x);

    // RDKit✔️✔️:   g[2][0] += sinTerm *
    // RDKit✔️✔️:              (dCos_dT[1] * (r[0].z) + dCos_dT[2] * (-r[0].y) +
    // RDKit✔️✔️:               dCos_dT[4] * (r[3].z - r[2].z) + dCos_dT[5] * (r[2].y - r[3].y));
    // RDKit✔️✔️:   g[2][1] += sinTerm *
    // RDKit✔️✔️:              (dCos_dT[0] * (-r[0].z) + dCos_dT[2] * (r[0].x) +
    // RDKit✔️✔️:               dCos_dT[3] * (r[2].z - r[3].z) + dCos_dT[5] * (r[3].x - r[2].x));
    // RDKit✔️✔️:   g[2][2] += sinTerm *
    // RDKit✔️✔️:              (dCos_dT[0] * (r[0].y) + dCos_dT[1] * (-r[0].x) +
    // RDKit✔️✔️:               dCos_dT[3] * (r[3].y - r[2].y) + dCos_dT[4] * (r[2].x - r[3].x));
    grad[g2] += sin_term
        * (dcos_dt[1] * r[0].z
            + dcos_dt[2] * (-r[0].y)
            + dcos_dt[4] * (r[3].z - r[2].z)
            + dcos_dt[5] * (r[2].y - r[3].y));
    grad[g2 + 1] += sin_term
        * (dcos_dt[0] * (-r[0].z)
            + dcos_dt[2] * r[0].x
            + dcos_dt[3] * (r[2].z - r[3].z)
            + dcos_dt[5] * (r[3].x - r[2].x));
    grad[g2 + 2] += sin_term
        * (dcos_dt[0] * r[0].y
            + dcos_dt[1] * (-r[0].x)
            + dcos_dt[3] * (r[3].y - r[2].y)
            + dcos_dt[4] * (r[2].x - r[3].x));

    // RDKit✔️✔️:   g[3][0] += sinTerm * (dCos_dT[4] * r[2].z - dCos_dT[5] * r[2].y);
    // RDKit✔️✔️:   g[3][1] += sinTerm * (dCos_dT[5] * r[2].x - dCos_dT[3] * r[2].z);
    // RDKit✔️✔️:   g[3][2] += sinTerm * (dCos_dT[3] * r[2].y - dCos_dT[4] * r[2].x);
    // RDKit✔️✔️: }
    grad[g3] += sin_term * (dcos_dt[4] * r[2].z - dcos_dt[5] * r[2].y);
    grad[g3 + 1] += sin_term * (dcos_dt[5] * r[2].x - dcos_dt[3] * r[2].z);
    grad[g3 + 2] += sin_term * (dcos_dt[3] * r[2].y - dcos_dt[4] * r[2].x);
}

impl ForceFieldContrib for TorsionAngleContrib {
    fn copy(&self) -> Box<dyn ForceFieldContrib> {
        // RDKit✔️✔️: TorsionAngleContrib *copy() const override {
        // RDKit✔️✔️:   return new TorsionAngleContrib(*this);
        // RDKit✔️✔️: }
        Box::new(self.clone())
    }

    fn set_force_field(&mut self, owner: *const ForceField) {
        self.owner = owner;
    }

    fn get_energy(&self, pos: &[f64]) -> f64 {
        TorsionAngleContrib::get_energy(self, pos)
    }

    fn get_grad(&self, pos: &[f64], grad: &mut [f64]) {
        TorsionAngleContrib::get_grad(self, pos, grad);
    }
}

#[cfg(test)]
mod tests {
    use crate::chemistry::forcefield::core::ForceFieldVec3;

    use super::*;

    const EPS: f64 = 1.0e-12;

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < EPS,
            "actual {actual} expected {expected} diff {}",
            (actual - expected).abs()
        );
    }

    fn assert_close_tol(actual: f64, expected: f64, tol: f64) {
        assert!(
            (actual - expected).abs() < tol,
            "actual {actual} expected {expected} diff {} tol {tol}",
            (actual - expected).abs()
        );
    }

    fn force_field() -> ForceField {
        let mut ff = ForceField::new(3);
        ff.positions_mut().push(ForceFieldVec3::new(0.0, 0.0, 0.0));
        ff.positions_mut().push(ForceFieldVec3::new(1.0, 0.0, 0.0));
        ff.positions_mut().push(ForceFieldVec3::new(1.0, 1.0, 0.0));
        ff.positions_mut().push(ForceFieldVec3::new(0.0, 1.0, 1.0));
        ff
    }

    fn atomic_params(v1: f64, u1: f64) -> AtomicParams {
        AtomicParams {
            r1: 0.0,
            theta0: 0.0,
            x1: 0.0,
            d1: 0.0,
            zeta: 0.0,
            z1: 0.0,
            v1,
            u1,
            gmp_xi: 0.0,
            gmp_hardness: 0.0,
            gmp_radius: 0.0,
        }
    }

    fn finite_difference_grad(contrib: &TorsionAngleContrib, pos: &[f64; 12]) -> [f64; 12] {
        let mut out = [0.0; 12];
        let h = 1.0e-6;
        for i in 0..12 {
            let mut plus = *pos;
            let mut minus = *pos;
            plus[i] += h;
            minus[i] -= h;
            out[i] = (contrib.get_energy(&plus) - contrib.get_energy(&minus)) / (2.0 * h);
        }
        out
    }

    #[test]
    fn uff_torsionanglecontrib_constructor_stores_owner_and_indices() {
        let ff = force_field();
        let at2 = atomic_params(2.0, 1.5);
        let at3 = atomic_params(8.0, 2.5);

        let contrib = TorsionAngleContrib::new(
            &ff,
            0,
            1,
            2,
            3,
            1.0,
            6,
            6,
            Hybridization::Sp3,
            Hybridization::Sp3,
            &at2,
            &at3,
            false,
        );

        assert_eq!(contrib.owner(), &ff as *const ForceField);
        assert_eq!(contrib.at1_idx(), 0);
        assert_eq!(contrib.at2_idx(), 1);
        assert_eq!(contrib.at3_idx(), 2);
        assert_eq!(contrib.at4_idx(), 3);
    }

    #[test]
    fn uff_torsionanglecontrib_constructor_sp3_sp3_uses_general_case() {
        let ff = force_field();
        let at2 = atomic_params(2.0, 1.5);
        let at3 = atomic_params(8.0, 2.5);

        let contrib = TorsionAngleContrib::new(
            &ff,
            0,
            1,
            2,
            3,
            1.5,
            6,
            6,
            Hybridization::Sp3,
            Hybridization::Sp3,
            &at2,
            &at3,
            false,
        );

        assert_close(contrib.force_constant(), (at2.v1 * at3.v1).sqrt());
        assert_eq!(contrib.order(), 3);
        assert_close(contrib.cos_term(), -1.0);
    }

    #[test]
    fn uff_torsionanglecontrib_constructor_sp3_sp3_group6_single_bond_special_case() {
        let ff = force_field();
        let at2 = atomic_params(2.0, 1.5);
        let at3 = atomic_params(8.0, 2.5);

        let oxygen_sulfur = TorsionAngleContrib::new(
            &ff,
            0,
            1,
            2,
            3,
            1.0,
            8,
            16,
            Hybridization::Sp3,
            Hybridization::Sp3,
            &at2,
            &at3,
            false,
        );
        let sulfur_sulfur = TorsionAngleContrib::new(
            &ff,
            0,
            1,
            2,
            3,
            1.0,
            16,
            16,
            Hybridization::Sp3,
            Hybridization::Sp3,
            &at2,
            &at3,
            false,
        );

        assert_close(oxygen_sulfur.force_constant(), (2.0_f64 * 6.8).sqrt());
        assert_eq!(oxygen_sulfur.order(), 2);
        assert_close(oxygen_sulfur.cos_term(), -1.0);
        assert_close(sulfur_sulfur.force_constant(), (6.8_f64 * 6.8).sqrt());
        assert_eq!(sulfur_sulfur.order(), 2);
        assert_close(sulfur_sulfur.cos_term(), -1.0);
    }

    #[test]
    fn uff_torsionanglecontrib_constructor_sp2_sp2_uses_equation17() {
        let ff = force_field();
        let at2 = atomic_params(2.0, 1.5);
        let at3 = atomic_params(8.0, 2.5);
        let bond_order = 1.25;

        let contrib = TorsionAngleContrib::new(
            &ff,
            0,
            1,
            2,
            3,
            bond_order,
            6,
            6,
            Hybridization::Sp2,
            Hybridization::Sp2,
            &at2,
            &at3,
            false,
        );

        assert_close(contrib.force_constant(), equation17(bond_order, &at2, &at3));
        assert_eq!(contrib.order(), 2);
        assert_close(contrib.cos_term(), 1.0);
    }

    #[test]
    fn uff_torsionanglecontrib_constructor_sp2_sp3_covers_default_and_special_cases() {
        let ff = force_field();
        let at2 = atomic_params(2.0, 1.5);
        let at3 = atomic_params(8.0, 2.5);

        let default = TorsionAngleContrib::new(
            &ff,
            0,
            1,
            2,
            3,
            1.5,
            6,
            6,
            Hybridization::Sp2,
            Hybridization::Sp3,
            &at2,
            &at3,
            false,
        );
        let group6 = TorsionAngleContrib::new(
            &ff,
            0,
            1,
            2,
            3,
            1.0,
            6,
            8,
            Hybridization::Sp2,
            Hybridization::Sp3,
            &at2,
            &at3,
            false,
        );
        let end_sp2 = TorsionAngleContrib::new(
            &ff,
            0,
            1,
            2,
            3,
            1.0,
            6,
            6,
            Hybridization::Sp2,
            Hybridization::Sp3,
            &at2,
            &at3,
            true,
        );

        assert_close(default.force_constant(), 1.0);
        assert_eq!(default.order(), 6);
        assert_close(default.cos_term(), 1.0);
        assert_close(group6.force_constant(), equation17(1.0, &at2, &at3));
        assert_eq!(group6.order(), 2);
        assert_close(group6.cos_term(), -1.0);
        assert_close(end_sp2.force_constant(), 2.0);
        assert_eq!(end_sp2.order(), 3);
        assert_close(end_sp2.cos_term(), -1.0);
    }

    #[test]
    fn uff_torsionanglecontrib_get_energy_covers_order_branches() {
        let ff = force_field();
        let at2 = atomic_params(2.0, 1.5);
        let at3 = atomic_params(8.0, 2.5);
        let pos = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.25, 1.0, 1.0];
        let p1 = ForceFieldVec3::new(pos[0], pos[1], pos[2]);
        let p2 = ForceFieldVec3::new(pos[3], pos[4], pos[5]);
        let p3 = ForceFieldVec3::new(pos[6], pos[7], pos[8]);
        let p4 = ForceFieldVec3::new(pos[9], pos[10], pos[11]);
        let cos_phi = calculate_cos_torsion(p1, p2, p3, p4);
        let sin_phi_sq = 1.0 - cos_phi * cos_phi;
        let cases = [
            (
                TorsionAngleContrib::new(
                    &ff,
                    0,
                    1,
                    2,
                    3,
                    1.0,
                    6,
                    6,
                    Hybridization::Sp2,
                    Hybridization::Sp2,
                    &at2,
                    &at3,
                    false,
                ),
                1.0 - 2.0 * sin_phi_sq,
            ),
            (
                TorsionAngleContrib::new(
                    &ff,
                    0,
                    1,
                    2,
                    3,
                    1.5,
                    6,
                    6,
                    Hybridization::Sp3,
                    Hybridization::Sp3,
                    &at2,
                    &at3,
                    false,
                ),
                cos_phi * (cos_phi * cos_phi - 3.0 * sin_phi_sq),
            ),
            (
                TorsionAngleContrib::new(
                    &ff,
                    0,
                    1,
                    2,
                    3,
                    1.5,
                    6,
                    6,
                    Hybridization::Sp2,
                    Hybridization::Sp3,
                    &at2,
                    &at3,
                    false,
                ),
                1.0 + sin_phi_sq * (-32.0 * sin_phi_sq * sin_phi_sq + 48.0 * sin_phi_sq - 18.0),
            ),
        ];

        for (contrib, cos_n_phi) in cases {
            let expected = contrib.force_constant() / 2.0 * (1.0 - contrib.cos_term() * cos_n_phi);
            assert_close(contrib.get_energy(&pos), expected);
        }
    }

    #[test]
    fn uff_torsionanglecontrib_get_energy_zero_cross_product_uses_zero_cosine() {
        let ff = force_field();
        let at2 = atomic_params(2.0, 1.5);
        let at3 = atomic_params(8.0, 2.5);
        let contrib = TorsionAngleContrib::new(
            &ff,
            0,
            1,
            2,
            3,
            1.5,
            6,
            6,
            Hybridization::Sp3,
            Hybridization::Sp3,
            &at2,
            &at3,
            false,
        );
        let pos = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 3.0, 1.0, 0.0];
        let cos_phi = calculate_cos_torsion(
            ForceFieldVec3::new(pos[0], pos[1], pos[2]),
            ForceFieldVec3::new(pos[3], pos[4], pos[5]),
            ForceFieldVec3::new(pos[6], pos[7], pos[8]),
            ForceFieldVec3::new(pos[9], pos[10], pos[11]),
        );

        assert_close(cos_phi, 0.0);
        assert_close(contrib.get_energy(&pos), contrib.force_constant() / 2.0);
    }

    #[test]
    fn uff_torsionanglecontrib_get_grad_matches_finite_difference_for_order_branches() {
        let ff = force_field();
        let at2 = atomic_params(2.0, 1.5);
        let at3 = atomic_params(8.0, 2.5);
        let pos = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.25, 1.0, 1.0];
        let cases = [
            TorsionAngleContrib::new(
                &ff,
                0,
                1,
                2,
                3,
                1.0,
                6,
                6,
                Hybridization::Sp2,
                Hybridization::Sp2,
                &at2,
                &at3,
                false,
            ),
            TorsionAngleContrib::new(
                &ff,
                0,
                1,
                2,
                3,
                1.5,
                6,
                6,
                Hybridization::Sp3,
                Hybridization::Sp3,
                &at2,
                &at3,
                false,
            ),
            TorsionAngleContrib::new(
                &ff,
                0,
                1,
                2,
                3,
                1.5,
                6,
                6,
                Hybridization::Sp2,
                Hybridization::Sp3,
                &at2,
                &at3,
                false,
            ),
        ];

        for contrib in cases {
            let expected = finite_difference_grad(&contrib, &pos);
            let mut grad = [0.0; 12];
            contrib.get_grad(&pos, &mut grad);

            for i in 0..12 {
                assert_close_tol(grad[i], expected[i], 1.0e-5);
            }
        }
    }

    #[test]
    fn uff_torsionanglecontrib_get_grad_accumulates_and_handles_zero_sine_branch() {
        let ff = force_field();
        let at2 = atomic_params(2.0, 1.5);
        let at3 = atomic_params(8.0, 2.5);
        let contrib = TorsionAngleContrib::new(
            &ff,
            0,
            1,
            2,
            3,
            1.0,
            6,
            6,
            Hybridization::Sp2,
            Hybridization::Sp2,
            &at2,
            &at3,
            false,
        );
        let pos = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0];
        let mut grad = [1.25; 12];

        contrib.get_grad(&pos, &mut grad);

        for value in grad {
            assert!(value.is_finite());
            assert_close(value, 1.25);
        }
    }

    #[test]
    fn uff_torsionanglecontrib_calc_torsion_grad_accumulates_source_formula() {
        let r = [
            ForceFieldVec3::new(1.0, 2.0, 3.0),
            ForceFieldVec3::new(4.0, 5.0, 6.0),
            ForceFieldVec3::new(7.0, 8.0, 9.0),
            ForceFieldVec3::new(10.0, 11.0, 12.0),
        ];
        let t = [ForceFieldVec3::new(1.0, 0.0, 0.0), ForceFieldVec3::new(0.0, 1.0, 0.0)];
        let d = [2.0, 4.0];
        let mut grad = [1.0; 12];

        calc_torsion_grad(&r, &t, &d, &mut grad, [0, 1, 2, 3], 2.0, 0.25);

        let expected = [
            -5.0, -0.5, 6.25, 5.5, 7.75, -9.5, 3.625, 0.25, 1.375, -0.125, -3.5, 5.875,
        ];
        for i in 0..12 {
            assert_close(grad[i], expected[i]);
        }
    }

    #[test]
    #[should_panic]
    fn uff_torsionanglecontrib_constructor_rejects_degenerate_indices() {
        let ff = force_field();
        let at2 = atomic_params(2.0, 1.5);
        let at3 = atomic_params(8.0, 2.5);

        let _ = TorsionAngleContrib::new(
            &ff,
            0,
            1,
            2,
            2,
            1.0,
            6,
            6,
            Hybridization::Sp3,
            Hybridization::Sp3,
            &at2,
            &at3,
            false,
        );
    }

    #[test]
    #[should_panic]
    fn uff_torsionanglecontrib_constructor_rejects_index_out_of_range() {
        let ff = force_field();
        let at2 = atomic_params(2.0, 1.5);
        let at3 = atomic_params(8.0, 2.5);

        let _ = TorsionAngleContrib::new(
            &ff,
            0,
            1,
            2,
            4,
            1.0,
            6,
            6,
            Hybridization::Sp3,
            Hybridization::Sp3,
            &at2,
            &at3,
            false,
        );
    }

    #[test]
    #[should_panic]
    fn uff_torsionanglecontrib_constructor_rejects_bad_hybridization() {
        let ff = force_field();
        let at2 = atomic_params(2.0, 1.5);
        let at3 = atomic_params(8.0, 2.5);

        let _ = TorsionAngleContrib::new(
            &ff,
            0,
            1,
            2,
            3,
            1.0,
            6,
            6,
            Hybridization::Sp,
            Hybridization::Sp3,
            &at2,
            &at3,
            false,
        );
    }
}
