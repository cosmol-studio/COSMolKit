// Copyright (C) 2004-2006 Rational Discovery LLC.
//
// This file contains behavior ported from RDKit and remains subject to its
// BSD license, included in the pinned source tree.

use crate::geometry::{Point3, compute_dihedral_from_flat};
use crate::kernel::{
    EvaluationContext, ForceFieldContribution, ForceFieldKernelError, TorsionIndexArgument,
};
use cosmolkit_model::Hybridization;

use super::params::{AtomicParams, clip_to_one, is_double_zero};

// BEGIN RDKIT CPP FUNCTION ForceFields::UFF::Utils::isInGroup6 (ForceField/UFF/TorsionAngle.cpp:37-40)
pub(super) fn is_in_group6(num: i32) -> bool {
    // RDKit❗✔️: bool isInGroup6(int num) {
    // RDKit❗✔️:   return (num == 8 || num == 16 || num == 34 || num == 52 || num == 84);
    // RDKit❗✔️: }
    num == 8 || num == 16 || num == 34 || num == 52 || num == 84
}
// END RDKIT CPP FUNCTION ForceFields::UFF::Utils::isInGroup6

// BEGIN RDKIT CPP FUNCTION ForceFields::UFF::Utils::calculateCosTorsion (ForceField/UFF/TorsionAngle.cpp:20-35)
fn calculate_cos_torsion(p1: Point3, p2: Point3, p3: Point3, p4: Point3) -> f64 {
    // RDKit❗✔️: double calculateCosTorsion(const RDGeom::Point3D &p1, const RDGeom::Point3D &p2,
    // RDKit❗✔️:                            const RDGeom::Point3D &p3,
    // RDKit❗✔️:                            const RDGeom::Point3D &p4) {
    // RDKit❗✔️:   RDGeom::Point3D r1 = p1 - p2, r2 = p3 - p2, r3 = p2 - p3, r4 = p4 - p3;

    // BEGIN RDKIT CPP HELPER RDGeom::operator- (Geometry/point.cpp:64-70)
    // RDKit❗✔️: Point3D operator-(const Point3D &p1, const Point3D &p2) {
    // RDKit❗✔️:   Point3D res;
    // RDKit❗✔️:   res.x = p1.x - p2.x;
    // RDKit❗✔️:   res.y = p1.y - p2.y;
    // RDKit❗✔️:   res.z = p1.z - p2.z;
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP HELPER RDGeom::operator-
    let r1 = Point3::difference(&p1, &p2);
    let r2 = Point3::difference(&p3, &p2);
    let r3 = Point3::difference(&p2, &p3);
    let r4 = Point3::difference(&p4, &p3);

    // RDKit❗✔️:   RDGeom::Point3D t1 = r1.crossProduct(r2);
    // RDKit❗✔️:   RDGeom::Point3D t2 = r3.crossProduct(r4);
    // BEGIN RDKIT CPP HELPER RDGeom::Point3D::crossProduct (Geometry/point.h:228-234)
    // RDKit❗✔️: constexpr Point3D crossProduct(const Point3D &other) const {
    // RDKit❗✔️:   Point3D res;
    // RDKit❗✔️:   res.x = y * (other.z) - z * (other.y);
    // RDKit❗✔️:   res.y = -x * (other.z) + z * (other.x);
    // RDKit❗✔️:   res.z = x * (other.y) - y * (other.x);
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP HELPER RDGeom::Point3D::crossProduct
    let t1 = r1.cross_product(&r2);
    let t2 = r3.cross_product(&r4);

    // RDKit❗✔️:   double d1 = t1.length();
    // RDKit❗✔️:   double d2 = t2.length();
    // BEGIN RDKIT CPP HELPER RDGeom::Point3D::length (Geometry/point.h:158-161)
    // RDKit❗✔️: double length() const override {
    // RDKit❗✔️:   double res = x * x + y * y + z * z;
    // RDKit❗✔️:   return sqrt(res);
    // RDKit❗✔️: }
    // END RDKIT CPP HELPER RDGeom::Point3D::length
    let d1 = t1.length();
    let d2 = t2.length();

    // RDKit❗✔️:   if (isDoubleZero(d1) || isDoubleZero(d2)) {
    // RDKit❗✔️:     return 0.0;
    // RDKit❗✔️:   }
    // BEGIN RDKIT CPP HELPER ForceFields::UFF::Utils::isDoubleZero (ForceField/UFF/Params.h:29-31)
    // RDKit❗✔️: inline bool isDoubleZero(const double x) {
    // RDKit❗✔️:   return ((x < 1.0e-10) && (x > -1.0e-10));
    // RDKit❗✔️: }
    // END RDKIT CPP HELPER ForceFields::UFF::Utils::isDoubleZero
    if is_double_zero(d1) || is_double_zero(d2) {
        return 0.0;
    }

    // RDKit❗✔️:   double cosPhi = t1.dotProduct(t2) / (d1 * d2);
    // BEGIN RDKIT CPP HELPER RDGeom::Point3D::dotProduct (Geometry/point.h:169-172)
    // RDKit❗✔️: constexpr double dotProduct(const Point3D &other) const {
    // RDKit❗✔️:   double res = x * (other.x) + y * (other.y) + z * (other.z);
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP HELPER RDGeom::Point3D::dotProduct
    let mut cos_phi = (t1.x * t2.x + t1.y * t2.y + t1.z * t2.z) / (d1 * d2);

    // RDKit❗✔️:   clipToOne(cosPhi);
    // BEGIN RDKIT CPP HELPER ForceFields::UFF::Utils::clipToOne (ForceField/UFF/Params.h:32)
    // RDKit❗✔️: inline void clipToOne(double &x) { x = std::clamp(x, -1.0, 1.0); }
    // END RDKIT CPP HELPER ForceFields::UFF::Utils::clipToOne
    clip_to_one(&mut cos_phi);

    // RDKit❗✔️:   return cosPhi;
    // RDKit❗✔️: }
    cos_phi
}
// END RDKIT CPP FUNCTION ForceFields::UFF::Utils::calculateCosTorsion

// BEGIN RDKIT CPP FUNCTION ForceFields::UFF::Utils::equation17 (ForceField/UFF/TorsionAngle.cpp:42-47)
pub(super) fn equation17(
    bond_order23: f64,
    at2_params: &AtomicParams,
    at3_params: &AtomicParams,
) -> f64 {
    // RDKit❗✔️: double equation17(double bondOrder23, const AtomicParams *at2Params,
    // RDKit❗✔️:                   const AtomicParams *at3Params) {
    // RDKit❗✔️:   return 5. * sqrt(at2Params->U1 * at3Params->U1) *
    // RDKit❗✔️:          (1. + 4.18 * log(bondOrder23));
    // RDKit❗✔️: }
    5.0 * (at2_params.u1 * at3_params.u1).sqrt() * (1.0 + 4.18 * bond_order23.ln())
}
// END RDKIT CPP FUNCTION ForceFields::UFF::Utils::equation17

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum TorsionParamsError {
    BadHybridizations,
}

impl TorsionParamsError {
    const fn source_category(self) -> &'static str {
        match self {
            Self::BadHybridizations => "Pre-condition Violation",
        }
    }

    const fn source_message(self) -> &'static str {
        match self {
            Self::BadHybridizations => "bad hybridizations",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum TorsionContributionError {
    DegeneratePoints,
    IndexOutOfRange {
        argument: TorsionIndexArgument,
        index: u32,
        upper_bound: usize,
    },
    BadHybridizations,
    BadOrder {
        order: u32,
    },
}

impl TorsionContributionError {
    const fn source_category(self) -> &'static str {
        match self {
            Self::DegeneratePoints | Self::BadHybridizations | Self::BadOrder { .. } => {
                "Pre-condition Violation"
            }
            Self::IndexOutOfRange { .. } => "Range Error",
        }
    }

    const fn source_message(self) -> &'static str {
        match self {
            Self::DegeneratePoints => "degenerate points",
            Self::IndexOutOfRange { argument, .. } => match argument {
                TorsionIndexArgument::First => "idx1",
                TorsionIndexArgument::Second => "idx2",
                TorsionIndexArgument::Third => "idx3",
                TorsionIndexArgument::Fourth => "idx4",
            },
            Self::BadHybridizations => "bad hybridizations",
            Self::BadOrder { .. } => "bad order",
        }
    }

    const fn range_detail(self) -> Option<(u32, usize)> {
        match self {
            Self::IndexOutOfRange {
                index, upper_bound, ..
            } => Some((index, upper_bound)),
            _ => None,
        }
    }
}

impl From<TorsionParamsError> for TorsionContributionError {
    fn from(error: TorsionParamsError) -> Self {
        match error {
            TorsionParamsError::BadHybridizations => Self::BadHybridizations,
        }
    }
}

impl From<TorsionContributionError> for ForceFieldKernelError {
    fn from(error: TorsionContributionError) -> Self {
        // BEGIN RDKIT CPP ERROR CONTRACT ForceFields::UFF::TorsionAngleContrib
        // (ForceField/UFF/TorsionAngle.cpp:89-114, 116-123, 178-181, 247-250)
        // RDKit❗✔️: PRECONDITION((idx1 != idx2 && idx1 != idx3 && idx1 != idx4 && idx2 != idx3 &&
        // RDKit❗✔️:               idx2 != idx4 && idx3 != idx4), "degenerate points");
        // RDKit❗✔️: PRECONDITION((hyb2 == RDKit::Atom::SP2 || hyb2 == RDKit::Atom::SP3) &&
        // RDKit❗✔️:               (hyb3 == RDKit::Atom::SP2 || hyb3 == RDKit::Atom::SP3),
        // RDKit❗✔️:               "bad hybridizations");
        // RDKit❗✔️: PRECONDITION(d_order == 2 || d_order == 3 || d_order == 6, "bad order");
        // END RDKIT CPP ERROR CONTRACT ForceFields::UFF::TorsionAngleContrib
        // Behavior marker — RDKit❗✔️: keep each torsion failure in its source category,
        // message, and structured predicate/indices; do not nest or stringify it.
        // Complexity marker — RDKit✔️✔️: constant-time matching with no allocation.
        match error {
            TorsionContributionError::DegeneratePoints => Self::TorsionDegeneratePoints,
            TorsionContributionError::IndexOutOfRange {
                argument,
                index,
                upper_bound,
            } => Self::TorsionIndexOutOfRange {
                argument,
                index,
                upper_bound,
            },
            TorsionContributionError::BadHybridizations => Self::TorsionBadHybridizations,
            TorsionContributionError::BadOrder { order } => Self::TorsionBadOrder { order },
        }
    }
}

// BEGIN RDKIT CPP FUNCTION ForceFields::UFF::TorsionAngleContrib::calcTorsionParams (ForceField/UFF/TorsionAngle.cpp:116-176)
fn calc_torsion_params(
    bond_order23: f64,
    at_num2: i32,
    at_num3: i32,
    hyb2: Hybridization,
    hyb3: Hybridization,
    at2_params: &AtomicParams,
    at3_params: &AtomicParams,
    end_atom_is_sp2: bool,
) -> Result<(f64, u32, f64), TorsionParamsError> {
    // RDKit❗✔️: void TorsionAngleContrib::calcTorsionParams(double bondOrder23, int atNum2,
    // RDKit❗✔️:                                             int atNum3,
    // RDKit❗✔️:                                             RDKit::Atom::HybridizationType hyb2,
    // RDKit❗✔️:                                             RDKit::Atom::HybridizationType hyb3,
    // RDKit❗✔️:                                             const AtomicParams *at2Params,
    // RDKit❗✔️:                                             const AtomicParams *at3Params,
    // RDKit❗✔️:                                             bool endAtomIsSP2) {
    // RDKit❗✔️:   PRECONDITION((hyb2 == RDKit::Atom::SP2 || hyb2 == RDKit::Atom::SP3) &&
    // RDKit❗✔️:                    (hyb3 == RDKit::Atom::SP2 || hyb3 == RDKit::Atom::SP3),
    // RDKit❗✔️:                "bad hybridizations");

    // BEGIN RDKIT CPP HELPER PRECONDITION (RDGeneral/Invariant.h:108-114)
    // RDKit❗✔️: #define PRECONDITION(expr, mess) \
    // RDKit❗✔️:   if (!(expr)) { \
    // RDKit❗✔️:     Invar::Invariant inv("Pre-condition Violation", mess, #expr, __FILE__, \
    // RDKit❗✔️:                          __LINE__); \
    // RDKit❗✔️:     BOOST_LOG(rdErrorLog) << "\n\n****\n" << inv << "****\n" << std::endl; \
    // RDKit❗✔️:     throw inv; \
    // RDKit❗✔️:   }
    // END RDKIT CPP HELPER PRECONDITION
    if !((hyb2 == Hybridization::Sp2 || hyb2 == Hybridization::Sp3)
        && (hyb3 == Hybridization::Sp2 || hyb3 == Hybridization::Sp3))
    {
        return Err(TorsionParamsError::BadHybridizations);
    }

    // RDKit❗✔️:   if (hyb2 == RDKit::Atom::SP3 && hyb3 == RDKit::Atom::SP3) {
    let mut force_constant;
    let mut order;
    let mut cos_term;
    if hyb2 == Hybridization::Sp3 && hyb3 == Hybridization::Sp3 {
        // RDKit❗✔️:     // general case:
        // RDKit❗✔️:     d_forceConstant = sqrt(at2Params->V1 * at3Params->V1);
        force_constant = (at2_params.v1 * at3_params.v1).sqrt();
        // RDKit❗✔️:     d_order = 3;
        order = 3;
        // RDKit❗✔️:     d_cosTerm = -1;  // phi0=60
        cos_term = -1.0;

        // RDKit❗✔️:     // special case for single bonds between group 6 elements:
        // RDKit❗✔️:     if (bondOrder23 == 1.0 && Utils::isInGroup6(atNum2) &&
        // RDKit❗✔️:         Utils::isInGroup6(atNum3)) {
        if bond_order23 == 1.0 && is_in_group6(at_num2) && is_in_group6(at_num3) {
            // RDKit❗✔️:       double V2 = 6.8, V3 = 6.8;
            let mut v2: f64 = 6.8;
            let mut v3: f64 = 6.8;
            // RDKit❗✔️:       if (atNum2 == 8) {
            // RDKit❗✔️:         V2 = 2.0;
            // RDKit❗✔️:       }
            if at_num2 == 8 {
                v2 = 2.0;
            }
            // RDKit❗✔️:       if (atNum3 == 8) {
            // RDKit❗✔️:         V3 = 2.0;
            // RDKit❗✔️:       }
            if at_num3 == 8 {
                v3 = 2.0;
            }
            // RDKit❗✔️:       d_forceConstant = sqrt(V2 * V3);
            force_constant = (v2 * v3).sqrt();
            // RDKit❗✔️:       d_order = 2;
            order = 2;
            // RDKit❗✔️:       d_cosTerm = -1;  // phi0=90
            cos_term = -1.0;
        }
        // RDKit❗✔️:   } else if (hyb2 == RDKit::Atom::SP2 && hyb3 == RDKit::Atom::SP2) {
    } else if hyb2 == Hybridization::Sp2 && hyb3 == Hybridization::Sp2 {
        // RDKit❗✔️:     d_forceConstant = Utils::equation17(bondOrder23, at2Params, at3Params);
        force_constant = equation17(bond_order23, at2_params, at3_params);
        // RDKit❗✔️:     d_order = 2;
        order = 2;
        // RDKit❗✔️:     // FIX: is this angle term right?
        // RDKit❗✔️:     d_cosTerm = 1.0;  // phi0= 180
        cos_term = 1.0;
        // RDKit❗✔️:   } else {
    } else {
        // RDKit❗✔️:     // SP2 - SP3,  this is, by default, independent of atom type in UFF:
        // RDKit❗✔️:     d_forceConstant = 1.0;
        force_constant = 1.0;
        // RDKit❗✔️:     d_order = 6;
        order = 6;
        // RDKit❗✔️:     d_cosTerm = 1.0;  // phi0 = 0
        cos_term = 1.0;
        // RDKit❗✔️:     if (bondOrder23 == 1.0) {
        if bond_order23 == 1.0 {
            // RDKit❗✔️:       // special case between group 6 sp3 and non-group 6 sp2:
            // RDKit❗✔️:       if ((hyb2 == RDKit::Atom::SP3 && Utils::isInGroup6(atNum2) &&
            // RDKit❗✔️:            !Utils::isInGroup6(atNum3)) ||
            // RDKit❗✔️:           (hyb3 == RDKit::Atom::SP3 && Utils::isInGroup6(atNum3) &&
            // RDKit❗✔️:            !Utils::isInGroup6(atNum2))) {
            if (hyb2 == Hybridization::Sp3 && is_in_group6(at_num2) && !is_in_group6(at_num3))
                || (hyb3 == Hybridization::Sp3 && is_in_group6(at_num3) && !is_in_group6(at_num2))
            {
                // RDKit❗✔️:         d_forceConstant = Utils::equation17(bondOrder23, at2Params, at3Params);
                force_constant = equation17(bond_order23, at2_params, at3_params);
                // RDKit❗✔️:         d_order = 2;
                order = 2;
                // RDKit❗✔️:         d_cosTerm = -1;  // phi0 = 90;
                cos_term = -1.0;
            }
            // RDKit❗✔️:       // special case for sp3 - sp2 - sp2
            // RDKit❗✔️:       // (i.e. the sp2 has another sp2 neighbor, like propene)
            // RDKit❗✔️:       else if (endAtomIsSP2) {
            else if end_atom_is_sp2 {
                // RDKit❗✔️:         d_forceConstant = 2.0;
                force_constant = 2.0;
                // RDKit❗✔️:         d_order = 3;
                order = 3;
                // RDKit❗✔️:         d_cosTerm = -1;  // phi0 = 180;
                cos_term = -1.0;
            }
            // RDKit❗✔️:       }
            // RDKit❗✔️:     }
        }
        // RDKit❗✔️:   }
    }

    // RDKit❗✔️: }
    Ok((force_constant, order, cos_term))
}
// END RDKIT CPP FUNCTION ForceFields::UFF::TorsionAngleContrib::calcTorsionParams

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct TorsionAngleContrib {
    at1_idx: u32,
    at2_idx: u32,
    at3_idx: u32,
    at4_idx: u32,
    order: u32,
    force_constant: f64,
    cos_term: f64,
}

impl TorsionAngleContrib {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn new(
        positions: &[&mut [f64]],
        idx1: u32,
        idx2: u32,
        idx3: u32,
        idx4: u32,
        bond_order23: f64,
        at_num2: i32,
        at_num3: i32,
        hyb2: Hybridization,
        hyb3: Hybridization,
        at2_params: &AtomicParams,
        at3_params: &AtomicParams,
        end_atom_is_sp2: bool,
    ) -> Result<Self, TorsionContributionError> {
        // BEGIN RDKIT CPP FUNCTION ForceFields::UFF::TorsionAngleContrib::TorsionAngleContrib (ForceField/UFF/TorsionAngle.cpp:89-114)
        // RDKit❗✔️: TorsionAngleContrib::TorsionAngleContrib(
        // RDKit❗✔️:     ForceField *owner, unsigned int idx1, unsigned int idx2, unsigned int idx3,
        // RDKit❗✔️:     unsigned int idx4, double bondOrder23, int atNum2, int atNum3,
        // RDKit❗✔️:     RDKit::Atom::HybridizationType hyb2, RDKit::Atom::HybridizationType hyb3,
        // RDKit❗✔️:     const AtomicParams *at2Params, const AtomicParams *at3Params,
        // RDKit❗✔️:     bool endAtomIsSP2) {
        // RDKit❗✔️:   PRECONDITION(owner, "bad owner");
        // RDKit❗✔️:   PRECONDITION(at2Params, "bad params pointer");
        // RDKit❗✔️:   PRECONDITION(at3Params, "bad params pointer");
        // References and the borrowed position view make these null pointer
        // states unrepresentable; the contribution retains no owner pointer.
        // RDKit❗✔️:   PRECONDITION((idx1 != idx2 && idx1 != idx3 && idx1 != idx4 && idx2 != idx3 &&
        // RDKit❗✔️:                 idx2 != idx4 && idx3 != idx4),
        // RDKit❗✔️:                "degenerate points");
        // BEGIN RDKIT CPP HELPER PRECONDITION (RDGeneral/Invariant.h:108-114)
        // RDKit❗✔️: #define PRECONDITION(expr, mess) \
        // RDKit❗✔️:   if (!(expr)) { \
        // RDKit❗✔️:     Invar::Invariant inv("Pre-condition Violation", mess, #expr, __FILE__, \
        // RDKit❗✔️:                          __LINE__); \
        // RDKit❗✔️:     BOOST_LOG(rdErrorLog) << "\n\n****\n" << inv << "****\n" << std::endl; \
        // RDKit❗✔️:     throw inv; \
        // RDKit❗✔️:   }
        // END RDKIT CPP HELPER PRECONDITION
        if idx1 == idx2
            || idx1 == idx3
            || idx1 == idx4
            || idx2 == idx3
            || idx2 == idx4
            || idx3 == idx4
        {
            return Err(TorsionContributionError::DegeneratePoints);
        }

        // RDKit❗✔️:   URANGE_CHECK(idx1, owner->positions().size());
        Self::check_index(idx1, TorsionIndexArgument::First, positions.len())?;
        // RDKit❗✔️:   URANGE_CHECK(idx2, owner->positions().size());
        Self::check_index(idx2, TorsionIndexArgument::Second, positions.len())?;
        // RDKit❗✔️:   URANGE_CHECK(idx3, owner->positions().size());
        Self::check_index(idx3, TorsionIndexArgument::Third, positions.len())?;
        // RDKit❗✔️:   URANGE_CHECK(idx4, owner->positions().size());
        Self::check_index(idx4, TorsionIndexArgument::Fourth, positions.len())?;

        // RDKit❗✔️:   dp_forceField = owner;
        // RDKit❗✔️:   d_at1Idx = idx1;
        // RDKit❗✔️:   d_at2Idx = idx2;
        // RDKit❗✔️:   d_at3Idx = idx3;
        // RDKit❗✔️:   d_at4Idx = idx4;
        // RDKit❗✔️:   calcTorsionParams(bondOrder23, atNum2, atNum3, hyb2, hyb3, at2Params,
        // RDKit❗✔️:                     at3Params, endAtomIsSP2);
        // RDKit❗✔️: }
        let (force_constant, order, cos_term) = calc_torsion_params(
            bond_order23,
            at_num2,
            at_num3,
            hyb2,
            hyb3,
            at2_params,
            at3_params,
            end_atom_is_sp2,
        )?;

        Ok(Self {
            at1_idx: idx1,
            at2_idx: idx2,
            at3_idx: idx3,
            at4_idx: idx4,
            order,
            force_constant,
            cos_term,
        })
    }

    fn check_index(
        index: u32,
        argument: TorsionIndexArgument,
        upper_bound: usize,
    ) -> Result<(), TorsionContributionError> {
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
            return Err(TorsionContributionError::IndexOutOfRange {
                argument,
                index,
                upper_bound,
            });
        }
        Ok(())
    }

    pub(super) fn scale_force_constant(&mut self, count: u32) {
        // BEGIN RDKIT CPP HELPER TorsionAngleContrib::scaleForceConstant (ForceField/UFF/TorsionAngle.h:69-71)
        // RDKit❗✔️: void scaleForceConstant(unsigned int count) {
        // RDKit❗✔️:   this->d_forceConstant /= static_cast<double>(count);
        // RDKit❗✔️: }
        self.force_constant /= f64::from(count);
    }

    pub(super) fn get_energy(
        &self,
        context: &mut EvaluationContext<'_>,
    ) -> Result<f64, TorsionContributionError> {
        // BEGIN RDKIT CPP FUNCTION ForceFields::UFF::TorsionAngleContrib::getEnergy (ForceField/UFF/TorsionAngle.cpp:178-218)
        // RDKit❗✔️: double TorsionAngleContrib::getEnergy(double *pos) const {
        // RDKit❗✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit❗✔️:   PRECONDITION(pos, "bad vector");
        // Borrowed `self` and `EvaluationContext` make owner/position null
        // states unrepresentable; construction establishes valid source state.
        // RDKit❗✔️:   PRECONDITION(d_order == 2 || d_order == 3 || d_order == 6, "bad order");
        if !matches!(self.order, 2 | 3 | 6) {
            return Err(TorsionContributionError::BadOrder { order: self.order });
        }

        let order = self.order;
        let coordinates = context.coordinates();
        // RDKit❗✔️:   RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
        // RDKit❗✔️:                      pos[3 * d_at1Idx + 2]);
        // RDKit❗✔️:   RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
        // RDKit❗✔️:                      pos[3 * d_at2Idx + 2]);
        // RDKit❗✔️:   RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
        // RDKit❗✔️:                      pos[3 * d_at3Idx + 2]);
        // RDKit❗✔️:   RDGeom::Point3D p4(pos[3 * d_at4Idx], pos[3 * d_at4Idx + 1],
        // RDKit❗✔️:                      pos[3 * d_at4Idx + 2]);
        // The constructor receives unsigned indices, but the source header
        // stores them as `int` and multiplies in signed `int`; the source has
        // no defined behavior after narrowing or offset overflow. Rust uses
        // `usize` for safe offsets on the source-defined coordinate domain.
        let p1_base = self.at1_idx as usize * 3;
        let p1 = Point3 {
            x: coordinates[p1_base],
            y: coordinates[p1_base + 1],
            z: coordinates[p1_base + 2],
        };
        let p2_base = self.at2_idx as usize * 3;
        let p2 = Point3 {
            x: coordinates[p2_base],
            y: coordinates[p2_base + 1],
            z: coordinates[p2_base + 2],
        };
        let p3_base = self.at3_idx as usize * 3;
        let p3 = Point3 {
            x: coordinates[p3_base],
            y: coordinates[p3_base + 1],
            z: coordinates[p3_base + 2],
        };
        let p4_base = self.at4_idx as usize * 3;
        let p4 = Point3 {
            x: coordinates[p4_base],
            y: coordinates[p4_base + 1],
            z: coordinates[p4_base + 2],
        };

        // RDKit❗✔️:   double cosPhi = Utils::calculateCosTorsion(p1, p2, p3, p4);
        let cos_phi = calculate_cos_torsion(p1, p2, p3, p4);
        // RDKit❗✔️:   double sinPhiSq = 1 - cosPhi * cosPhi;
        let sin_phi_sq = 1.0 - cos_phi * cos_phi;

        // RDKit❗✔️:   // E(phi) = V/2 * (1 - cos(n*phi_0)*cos(n*phi))
        // RDKit❗✔️:   double cosNPhi = 0.0;
        // RDKit❗✔️:   switch (d_order) {
        let cos_n_phi = match order {
            // RDKit❗✔️:     case 2:
            // RDKit❗✔️:       // cos(2x) = 1 - 2sin^2(x)
            // RDKit❗✔️:       cosNPhi = 1 - 2 * sinPhiSq;
            // RDKit❗✔️:       break;
            2 => 1.0 - 2.0 * sin_phi_sq,
            // RDKit❗✔️:     case 3:
            // RDKit❗✔️:       // cos(3x) = cos^3(x) - 3*cos(x)*sin^2(x) = 4cos^3(x) -3cos(x)
            // RDKit❗✔️:       cosNPhi = cosPhi * (cosPhi * cosPhi - 3. * sinPhiSq);
            // RDKit❗✔️:       break;
            3 => cos_phi * (cos_phi * cos_phi - 3.0 * sin_phi_sq),
            // RDKit❗✔️:     case 6:
            // RDKit❗✔️:       // cos(6x) = 1 - 32*sin^6(x) + 48*sin^4(x) - 18*sin^2(x)
            // RDKit❗✔️:       cosNPhi =
            // RDKit❗✔️:           1 + sinPhiSq * (-32. * sinPhiSq * sinPhiSq + 48. * sinPhiSq - 18.);
            // RDKit❗✔️:       break;
            6 => 1.0 + sin_phi_sq * (-32.0 * sin_phi_sq * sin_phi_sq + 48.0 * sin_phi_sq - 18.0),
            _ => unreachable!("order was checked above"),
        };
        // RDKit❗✔️:   }
        // RDKit❗✔️:   double res = d_forceConstant / 2.0 * (1. - d_cosTerm * cosNPhi);
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        Ok(self.force_constant / 2.0 * (1.0 - self.cos_term * cos_n_phi))
    }

    pub(super) fn get_grad(
        &self,
        context: &mut EvaluationContext<'_>,
        gradient: &mut [f64],
    ) -> Result<(), TorsionContributionError> {
        // BEGIN RDKIT CPP FUNCTION ForceFields::UFF::TorsionAngleContrib::getGrad (ForceField/UFF/TorsionAngle.cpp:222-245)
        // RDKit❗✔️: void TorsionAngleContrib::getGrad(double *pos, double *grad) const {
        // RDKit❗✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit❗✔️:   PRECONDITION(pos, "bad vector");
        // RDKit❗✔️:   PRECONDITION(grad, "bad vector");
        // Borrowed self, EvaluationContext, and mutable slice make these
        // non-null source preconditions unrepresentable; `new` validated IDs.
        // RDKit❗✔️:   double *g[4] = {&(grad[3 * d_at1Idx]), &(grad[3 * d_at2Idx]),
        // RDKit❗✔️:                   &(grad[3 * d_at3Idx]), &(grad[3 * d_at4Idx])};
        // Rust keeps the four addresses as source-ordered indices and adds
        // into the caller's flat gradient without allocating pointer arrays.

        // RDKit❗✔️:   RDGeom::Point3D r[4];
        let mut r = [Point3::default(); 4];
        // RDKit❗✔️:   RDGeom::Point3D t[2];
        let mut t = [Point3::default(); 2];
        // RDKit❗✔️:   double d[2];
        let mut d = [0.0; 2];
        // RDKit❗✔️:   double cosPhi;
        let mut cos_phi = 0.0;
        // RDKit❗✔️:   RDKit::ForceFieldsHelper::computeDihedral(
        // RDKit❗✔️:       pos, d_at1Idx, d_at2Idx, d_at3Idx, d_at4Idx, nullptr, &cosPhi, r, t, d);
        // The canonical flat helper retains the pinned computeDihedral
        // normalization, denominator floors, clamps, and stack outputs.
        compute_dihedral_from_flat(
            context.coordinates(),
            self.at1_idx as usize,
            self.at2_idx as usize,
            self.at3_idx as usize,
            self.at4_idx as usize,
            None,
            Some(&mut cos_phi),
            Some(&mut r),
            Some(&mut t),
            Some(&mut d),
        );

        // RDKit❗✔️:   double sinPhiSq = 1.0 - cosPhi * cosPhi;
        let sin_phi_sq = 1.0 - cos_phi * cos_phi;
        // RDKit❗✔️:   double sinPhi = ((sinPhiSq > 0.0) ? sqrt(sinPhiSq) : 0.0);
        let sin_phi = if sin_phi_sq > 0.0 {
            sin_phi_sq.sqrt()
        } else {
            0.0
        };

        // RDKit❗✔️:   // dE/dPhi is independent of cartesians:
        // RDKit❗✔️:   double dE_dPhi = getThetaDeriv(cosPhi, sinPhi);
        let de_dphi = self.get_theta_deriv(cos_phi, sin_phi)?;
        // RDKit❗✔️:   double sinTerm =
        // RDKit❗✔️:       dE_dPhi * (isDoubleZero(sinPhi) ? (1.0 / cosPhi) : (1.0 / sinPhi));
        let sin_term = de_dphi
            * if is_double_zero(sin_phi) {
                1.0 / cos_phi
            } else {
                1.0 / sin_phi
            };

        // RDKit❗✔️:   Utils::calcTorsionGrad(r, t, d, g, sinTerm, cosPhi);
        calc_torsion_grad(
            &r,
            &t,
            &d,
            gradient,
            [self.at1_idx, self.at2_idx, self.at3_idx, self.at4_idx],
            sin_term,
            cos_phi,
        );
        // RDKit❗✔️: }
        Ok(())
    }

    fn get_theta_deriv(
        &self,
        cos_theta: f64,
        sin_theta: f64,
    ) -> Result<f64, TorsionContributionError> {
        // BEGIN RDKIT CPP FUNCTION ForceFields::UFF::TorsionAngleContrib::getThetaDeriv (ForceField/UFF/TorsionAngle.cpp:247-270)
        // RDKit❗✔️: double TorsionAngleContrib::getThetaDeriv(double cosTheta,
        // RDKit❗✔️:                                           double sinTheta) const {
        // RDKit❗✔️:   PRECONDITION(d_order == 2 || d_order == 3 || d_order == 6, "bad order");
        if !matches!(self.order, 2 | 3 | 6) {
            return Err(TorsionContributionError::BadOrder { order: self.order });
        }
        // RDKit❗✔️:   double sinThetaSq = sinTheta * sinTheta;
        let sin_theta_sq = sin_theta * sin_theta;
        // RDKit❗✔️:   // cos(6x) = 1 - 32*sin^6(x) + 48*sin^4(x) - 18*sin^2(x)

        // RDKit❗✔️:   double res = 0.0;
        // RDKit❗✔️:   switch (d_order) {
        let mut res = match self.order {
            // RDKit❗✔️:     case 2:
            // RDKit❗✔️:       res = 2 * sinTheta * cosTheta;
            // RDKit❗✔️:       break;
            2 => 2.0 * sin_theta * cos_theta,
            // RDKit❗✔️:     case 3:
            // RDKit❗✔️:       // sin(3*x) = 3*sin(x) - 4*sin^3(x)
            // RDKit❗✔️:       res = sinTheta * (3 - 4 * sinThetaSq);
            // RDKit❗✔️:       break;
            3 => sin_theta * (3.0 - 4.0 * sin_theta_sq),
            // RDKit❗✔️:     case 6:
            // RDKit❗✔️:       // sin(6x) = cos(x) * [ 32*sin^5(x) - 32*sin^3(x) + 6*sin(x) ]
            // RDKit❗✔️:       res = cosTheta * sinTheta * (32 * sinThetaSq * (sinThetaSq - 1) + 6);
            // RDKit❗✔️:       break;
            6 => cos_theta * sin_theta * (32.0 * sin_theta_sq * (sin_theta_sq - 1.0) + 6.0),
            _ => unreachable!("order was checked above"),
        };
        // RDKit❗✔️:   }
        // RDKit❗✔️:   res *= d_forceConstant / 2.0 * d_cosTerm * -1 * d_order;
        res *= self.force_constant / 2.0 * self.cos_term * -1.0 * f64::from(self.order);

        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        Ok(res)
    }
}

impl ForceFieldContribution for TorsionAngleContrib {
    #[cfg(test)]
    fn cf3d_frag_accept_test_identity(&self) -> crate::kernel::Cf3dFragAcceptContributionIdentity {
        crate::kernel::Cf3dFragAcceptContributionIdentity::TorsionAngle {
            at1_idx: self.at1_idx,
            at2_idx: self.at2_idx,
            at3_idx: self.at3_idx,
            at4_idx: self.at4_idx,
            order: self.order,
            force_constant: self.force_constant,
            cos_term: self.cos_term,
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
        // Behavior marker — RDKit❗✔️: delegate to the existing source-shaped torsion
        // evaluator and preserve its typed failure for the ordered field loop.
        // Complexity marker — RDKit✔️✔️: one direct call and fixed-size error mapping;
        // no extra coordinate read, geometry work, or allocation.
        TorsionAngleContrib::get_energy(self, context).map_err(ForceFieldKernelError::from)
    }

    fn get_grad(
        &self,
        context: &mut EvaluationContext<'_>,
        gradient: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        // BEGIN RDKIT CPP CALL ForceFields::ForceField::calcGrad(pos)
        // (ForceField/ForceField.cpp:363)
        // RDKit❗✔️:     (*contrib)->getGrad(pos, grad);
        // END RDKIT CPP CALL ForceFields::ForceField::calcGrad(pos)
        // Behavior marker — RDKit❗✔️: keep the existing source geometry/derivative order,
        // additive writes, and immediate typed failure propagation.
        // Complexity marker — RDKit✔️✔️: one direct call and fixed-size error mapping;
        // no extra coordinate read, geometry work, or allocation.
        TorsionAngleContrib::get_grad(self, context, gradient).map_err(ForceFieldKernelError::from)
    }

    fn copy(&self) -> Box<dyn ForceFieldContribution> {
        // BEGIN RDKIT CPP FUNCTION ForceFields::UFF::TorsionAngleContrib::copy
        // (ForceField/UFF/TorsionAngle.h:72-74)
        // RDKit❗✔️: TorsionAngleContrib *copy() const override {
        // RDKit❗✔️:   return new TorsionAngleContrib(*this);
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION ForceFields::UFF::TorsionAngleContrib::copy
        // Behavior marker — RDKit❗✔️: copy the complete owner-free value, including
        // its current scaled force constant.
        // Complexity marker — RDKit✔️✔️: fixed-size value copy plus the trait's one Box.
        Box::new(*self)
    }
}

#[allow(clippy::too_many_arguments)]
fn calc_torsion_grad(
    r: &[Point3; 4],
    t: &[Point3; 2],
    d: &[f64; 2],
    gradient: &mut [f64],
    indices: [u32; 4],
    sin_term: f64,
    cos_phi: f64,
) {
    // BEGIN RDKIT CPP FUNCTION ForceFields::UFF::Utils::calcTorsionGrad (ForceField/UFF/TorsionAngle.cpp:48-86)
    // RDKit❗✔️: void calcTorsionGrad(RDGeom::Point3D *r, RDGeom::Point3D *t, double *d,
    // RDKit❗✔️:                      double **g, double &sinTerm, double &cosPhi) {
    // RDKit❗✔️:   // -------
    // RDKit❗✔️:   // dTheta/dx is trickier:
    // RDKit❗✔️:   double dCos_dT[6] = {1.0 / d[0] * (t[1].x - cosPhi * t[0].x),
    // RDKit❗✔️:                        1.0 / d[0] * (t[1].y - cosPhi * t[0].y),
    // RDKit❗✔️:                        1.0 / d[0] * (t[1].z - cosPhi * t[0].z),
    // RDKit❗✔️:                        1.0 / d[1] * (t[0].x - cosPhi * t[1].x),
    // RDKit❗✔️:                        1.0 / d[1] * (t[0].y - cosPhi * t[1].y),
    // RDKit❗✔️:                        1.0 / d[1] * (t[0].z - cosPhi * t[1].z)};
    // RDKit❗✔️:
    // RDKit❗✔️:   g[0][0] += sinTerm * (dCos_dT[2] * r[1].y - dCos_dT[1] * r[1].z);
    // RDKit❗✔️:   g[0][1] += sinTerm * (dCos_dT[0] * r[1].z - dCos_dT[2] * r[1].x);
    // RDKit❗✔️:   g[0][2] += sinTerm * (dCos_dT[1] * r[1].x - dCos_dT[0] * r[1].y);
    // RDKit❗✔️:
    // RDKit❗✔️:   g[1][0] += sinTerm *
    // RDKit❗✔️:              (dCos_dT[1] * (r[1].z - r[0].z) + dCos_dT[2] * (r[0].y - r[1].y) +
    // RDKit❗✔️:               dCos_dT[4] * (-r[3].z) + dCos_dT[5] * (r[3].y));
    // RDKit❗✔️:   g[1][1] += sinTerm *
    // RDKit❗✔️:              (dCos_dT[0] * (r[0].z - r[1].z) + dCos_dT[2] * (r[1].x - r[0].x) +
    // RDKit❗✔️:               dCos_dT[3] * (r[3].z) + dCos_dT[5] * (-r[3].x));
    // RDKit❗✔️:   g[1][2] += sinTerm *
    // RDKit❗✔️:              (dCos_dT[0] * (r[1].y - r[0].y) + dCos_dT[1] * (r[0].x - r[1].x) +
    // RDKit❗✔️:               dCos_dT[3] * (-r[3].y) + dCos_dT[4] * (r[3].x));
    // RDKit❗✔️:
    // RDKit❗✔️:   g[2][0] += sinTerm *
    // RDKit❗✔️:              (dCos_dT[1] * (r[0].z) + dCos_dT[2] * (-r[0].y) +
    // RDKit❗✔️:               dCos_dT[4] * (r[3].z - r[2].z) + dCos_dT[5] * (r[2].y - r[3].y));
    // RDKit❗✔️:   g[2][1] += sinTerm *
    // RDKit❗✔️:              (dCos_dT[0] * (-r[0].z) + dCos_dT[2] * (r[0].x) +
    // RDKit❗✔️:               dCos_dT[3] * (r[2].z - r[3].z) + dCos_dT[5] * (r[3].x - r[2].x));
    // RDKit❗✔️:   g[2][2] += sinTerm *
    // RDKit❗✔️:              (dCos_dT[0] * (r[0].y) + dCos_dT[1] * (-r[0].x) +
    // RDKit❗✔️:               dCos_dT[3] * (r[3].y - r[2].y) + dCos_dT[4] * (r[2].x - r[3].x));
    // RDKit❗✔️:
    // RDKit❗✔️:   g[3][0] += sinTerm * (dCos_dT[4] * r[2].z - dCos_dT[5] * r[2].y);
    // RDKit❗✔️:   g[3][1] += sinTerm * (dCos_dT[5] * r[2].x - dCos_dT[3] * r[2].z);
    // RDKit❗✔️:   g[3][2] += sinTerm * (dCos_dT[3] * r[2].y - dCos_dT[4] * r[2].x);
    // RDKit❗✔️: }
    // Constructor indices are stored as signed `int` in the source; the
    // defined in-range domain is preserved, while source narrowing overflow
    // remains outside the proven domain.
    let d_cos_d_t = [
        1.0 / d[0] * (t[1].x - cos_phi * t[0].x),
        1.0 / d[0] * (t[1].y - cos_phi * t[0].y),
        1.0 / d[0] * (t[1].z - cos_phi * t[0].z),
        1.0 / d[1] * (t[0].x - cos_phi * t[1].x),
        1.0 / d[1] * (t[0].y - cos_phi * t[1].y),
        1.0 / d[1] * (t[0].z - cos_phi * t[1].z),
    ];

    let base0 = indices[0] as usize * 3;
    let base1 = indices[1] as usize * 3;
    let base2 = indices[2] as usize * 3;
    let base3 = indices[3] as usize * 3;
    gradient[base0] += sin_term * (d_cos_d_t[2] * r[1].y - d_cos_d_t[1] * r[1].z);
    gradient[base0 + 1] += sin_term * (d_cos_d_t[0] * r[1].z - d_cos_d_t[2] * r[1].x);
    gradient[base0 + 2] += sin_term * (d_cos_d_t[1] * r[1].x - d_cos_d_t[0] * r[1].y);

    gradient[base1] += sin_term
        * (d_cos_d_t[1] * (r[1].z - r[0].z)
            + d_cos_d_t[2] * (r[0].y - r[1].y)
            + d_cos_d_t[4] * (-r[3].z)
            + d_cos_d_t[5] * r[3].y);
    gradient[base1 + 1] += sin_term
        * (d_cos_d_t[0] * (r[0].z - r[1].z)
            + d_cos_d_t[2] * (r[1].x - r[0].x)
            + d_cos_d_t[3] * r[3].z
            + d_cos_d_t[5] * (-r[3].x));
    gradient[base1 + 2] += sin_term
        * (d_cos_d_t[0] * (r[1].y - r[0].y)
            + d_cos_d_t[1] * (r[0].x - r[1].x)
            + d_cos_d_t[3] * (-r[3].y)
            + d_cos_d_t[4] * r[3].x);

    gradient[base2] += sin_term
        * (d_cos_d_t[1] * r[0].z
            + d_cos_d_t[2] * (-r[0].y)
            + d_cos_d_t[4] * (r[3].z - r[2].z)
            + d_cos_d_t[5] * (r[2].y - r[3].y));
    gradient[base2 + 1] += sin_term
        * (d_cos_d_t[0] * (-r[0].z)
            + d_cos_d_t[2] * r[0].x
            + d_cos_d_t[3] * (r[2].z - r[3].z)
            + d_cos_d_t[5] * (r[3].x - r[2].x));
    gradient[base2 + 2] += sin_term
        * (d_cos_d_t[0] * r[0].y
            + d_cos_d_t[1] * (-r[0].x)
            + d_cos_d_t[3] * (r[3].y - r[2].y)
            + d_cos_d_t[4] * (r[2].x - r[3].x));

    gradient[base3] += sin_term * (d_cos_d_t[4] * r[2].z - d_cos_d_t[5] * r[2].y);
    gradient[base3 + 1] += sin_term * (d_cos_d_t[5] * r[2].x - d_cos_d_t[3] * r[2].z);
    gradient[base3 + 2] += sin_term * (d_cos_d_t[3] * r[2].y - d_cos_d_t[4] * r[2].x);
    // END RDKIT CPP FUNCTION ForceFields::UFF::Utils::calcTorsionGrad
}

#[cfg(test)]
mod tests {
    use super::{
        AtomicParams, TorsionAngleContrib, TorsionContributionError, TorsionParamsError,
        calc_torsion_grad, calc_torsion_params, calculate_cos_torsion, equation17, is_double_zero,
        is_in_group6,
    };
    use crate::geometry::Point3;
    use crate::kernel::{
        EvaluationContext, ForceField, ForceFieldContribution, ForceFieldKernelError,
        TorsionIndexArgument, cf3d_bld_b05_calc_energy, cf3d_bld_b05_calc_grad,
        cf3d_bld_b05_copy_force_field,
    };
    use cosmolkit_model::Hybridization;
    use std::f64::consts::PI;
    use std::sync::Arc;
    use std::sync::atomic::{AtomicUsize, Ordering};

    fn point(x: f64, y: f64, z: f64) -> Point3 {
        Point3 { x, y, z }
    }

    fn atomic_params(u1: f64) -> AtomicParams {
        torsion_params(0.0, u1)
    }

    fn torsion_params(v1: f64, u1: f64) -> AtomicParams {
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

    fn make_contribution(
        indices: [u32; 4],
        bond_order23: f64,
        at_num2: i32,
        at_num3: i32,
        hyb2: Hybridization,
        hyb3: Hybridization,
        end_atom_is_sp2: bool,
    ) -> Result<TorsionAngleContrib, TorsionContributionError> {
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
        let at2_params = torsion_params(4.0, 4.0);
        let at3_params = torsion_params(9.0, 9.0);

        TorsionAngleContrib::new(
            &positions,
            indices[0],
            indices[1],
            indices[2],
            indices[3],
            bond_order23,
            at_num2,
            at_num3,
            hyb2,
            hyb3,
            &at2_params,
            &at3_params,
            end_atom_is_sp2,
        )
    }

    fn torsion_coordinates(phi: f64) -> Vec<f64> {
        let (sin_phi, cos_phi) = phi.sin_cos();
        vec![
            0.0, 1.0, 0.0, // p1
            0.0, 0.0, 0.0, // p2
            1.0, 0.0, 0.0, // p3
            1.0, cos_phi, sin_phi, // p4
        ]
    }

    fn rows_from_coordinates(coordinates: &[f64]) -> [[f64; 3]; 4] {
        [
            [coordinates[0], coordinates[1], coordinates[2]],
            [coordinates[3], coordinates[4], coordinates[5]],
            [coordinates[6], coordinates[7], coordinates[8]],
            [coordinates[9], coordinates[10], coordinates[11]],
        ]
    }

    fn field_with_torsion<'a>(
        rows: &'a mut [[f64; 3]; 4],
        contribution: TorsionAngleContrib,
    ) -> ForceField<'a> {
        let mut field = ForceField::new(3);
        for row in rows {
            field.positions_mut().push(row.as_mut_slice());
        }
        field.add_contribution(Box::new(contribution));
        field.initialize().expect("four 3D positions initialize");
        field
    }

    #[derive(Clone)]
    struct CountContributionCalls(Arc<AtomicUsize>);

    impl ForceFieldContribution for CountContributionCalls {
        fn get_energy(
            &self,
            _context: &mut EvaluationContext<'_>,
        ) -> Result<f64, ForceFieldKernelError> {
            self.0.fetch_add(1, Ordering::Relaxed);
            Ok(0.0)
        }

        fn get_grad(
            &self,
            _context: &mut EvaluationContext<'_>,
            _gradient: &mut [f64],
        ) -> Result<(), ForceFieldKernelError> {
            self.0.fetch_add(1, Ordering::Relaxed);
            Ok(())
        }

        fn copy(&self) -> Box<dyn ForceFieldContribution> {
            Box::new(self.clone())
        }
    }

    fn energy_at(contribution: &TorsionAngleContrib, phi: f64) -> f64 {
        let coordinates = torsion_coordinates(phi);
        let mut distance_cache = vec![0.0; 10];
        let mut context = EvaluationContext::for_test(&coordinates, &mut distance_cache, 4);
        contribution.get_energy(&mut context).unwrap()
    }

    fn assert_close(actual: f64, expected: f64) {
        let tolerance = 1.0e-12 * expected.abs().max(1.0);
        assert!(
            (actual - expected).abs() <= tolerance,
            "actual {actual:.17e} differs from expected {expected:.17e} by more than {tolerance:.3e}"
        );
    }

    fn gradient_for_test(
        contribution: &TorsionAngleContrib,
        coordinates: &[f64],
        mut gradient: [f64; 12],
    ) -> Result<[f64; 12], TorsionContributionError> {
        let mut distance_cache = [0.0; 10];
        let mut context = EvaluationContext::for_test(coordinates, &mut distance_cache, 4);
        contribution.get_grad(&mut context, &mut gradient)?;
        Ok(gradient)
    }

    fn assert_gradient_close(actual: [f64; 12], expected: [f64; 12]) {
        for (actual, expected) in actual.into_iter().zip(expected) {
            assert_close(actual, expected);
        }
    }

    #[test]
    fn cf3d_u10_group_six_membership_is_exact() {
        // RDKit ForceField/UFF/TorsionAngle.cpp::isInGroup6 lists precisely
        // these five atomic numbers and compares against signed `int` values.
        for member in [8, 16, 34, 52, 84] {
            assert!(is_in_group6(member), "{member} must be in group six");
        }
        for nonmember in [-1, 0, 1, 7, 9, 15, 17, 33, 35, 51, 53, 83, 85, 118] {
            assert!(!is_in_group6(nonmember), "{nonmember} is not in group six");
        }
    }

    #[test]
    fn cf3d_u10_cos_torsion_preserves_normal_orientation_signs() {
        // RDKit ForceField/UFF/TorsionAngle.cpp::calculateCosTorsion forms
        // each normal in the stated point order, so reversing p4's side flips
        // the sign while perpendicular normals produce zero.
        let p1 = point(0.0, 1.0, 0.0);
        let p2 = point(0.0, 0.0, 0.0);
        let p3 = point(1.0, 0.0, 0.0);

        assert_eq!(calculate_cos_torsion(p1, p2, p3, point(1.0, 1.0, 0.0)), 1.0);
        assert_eq!(
            calculate_cos_torsion(p1, p2, p3, point(1.0, -1.0, 0.0)),
            -1.0
        );
        assert_eq!(calculate_cos_torsion(p1, p2, p3, point(1.0, 0.0, 1.0)), 0.0);
    }

    #[test]
    fn cf3d_u10_cos_torsion_returns_zero_for_first_degenerate_normal() {
        // RDKit checks isDoubleZero(d1) first and returns the literal 0.0.
        let value = calculate_cos_torsion(
            point(0.0, 0.0, 0.0),
            point(0.0, 0.0, 0.0),
            point(1.0, 0.0, 0.0),
            point(1.0, 1.0, 0.0),
        );

        assert_eq!(value.to_bits(), 0.0_f64.to_bits());
    }

    #[test]
    fn cf3d_u10_cos_torsion_returns_zero_for_second_degenerate_normal() {
        // RDKit reaches isDoubleZero(d2) when the first normal is nonzero.
        let value = calculate_cos_torsion(
            point(0.0, 1.0, 0.0),
            point(0.0, 0.0, 0.0),
            point(1.0, 0.0, 0.0),
            point(2.0, 0.0, 0.0),
        );

        assert_eq!(value.to_bits(), 0.0_f64.to_bits());
    }

    #[test]
    fn cf3d_u10_cos_torsion_keeps_the_strict_zero_threshold() {
        // Params.h::isDoubleZero uses strict comparisons at 1.0e-10.
        let p2 = point(0.0, 0.0, 0.0);
        let p3 = point(1.0, 0.0, 0.0);
        let at_threshold = 1.0e-10;
        let r1 = Point3::difference(&point(0.0, at_threshold, 0.0), &p2);
        let r2 = Point3::difference(&p3, &p2);
        let d1 = r1.cross_product(&r2).length();

        assert_eq!(d1.to_bits(), at_threshold.to_bits());
        assert!(!is_double_zero(d1));
        assert_eq!(
            calculate_cos_torsion(point(0.0, at_threshold, 0.0), p2, p3, point(1.0, 1.0, 0.0),),
            1.0
        );
        assert_eq!(
            calculate_cos_torsion(
                point(0.0, at_threshold * 0.5, 0.0),
                p2,
                p3,
                point(1.0, 1.0, 0.0),
            ),
            0.0
        );
    }

    #[test]
    fn cf3d_u10_cos_torsion_preserves_nan_through_source_clamp() {
        // Params.h::isDoubleZero comparisons and std::clamp leave NaN intact.
        let value = calculate_cos_torsion(
            point(0.0, 1.0, 0.0),
            point(0.0, 0.0, 0.0),
            point(1.0, 0.0, 0.0),
            point(1.0, f64::NAN, 0.0),
        );

        assert!(value.is_nan());
    }

    #[test]
    fn cf3d_u10_equation17_matches_fixed_source_values() {
        // RDKit ForceField/UFF/TorsionAngle.cpp::equation17 with U1=(4,9).
        let at2 = atomic_params(4.0);
        let at3 = atomic_params(9.0);

        assert_eq!(equation17(1.0, &at2, &at3), 30.0);
        assert_eq!(equation17(2.0, &at2, &at3), 116.92065644221714);
    }

    #[test]
    fn cf3d_u10_equation17_preserves_source_ieee_edges() {
        // The source has no guards around sqrt(U1_2*U1_3) or log(bondOrder23).
        let positive = atomic_params(4.0);
        let negative = atomic_params(-4.0);
        let zero = atomic_params(0.0);
        let positive_other = atomic_params(9.0);

        assert_eq!(
            equation17(0.0, &positive, &positive_other),
            f64::NEG_INFINITY
        );
        assert!(equation17(-1.0, &positive, &positive_other).is_nan());
        assert!(equation17(1.0, &negative, &positive_other).is_nan());
        assert!(equation17(0.0, &zero, &positive_other).is_nan());
    }

    #[test]
    fn cf3d_u11_precondition_accepts_only_source_hybridization_pairs() {
        // RDKit Atom.h defines these nine values; calcTorsionParams accepts
        // only SP2/SP3 for each central atom before any parameter assignment.
        let values = [
            Hybridization::Unspecified,
            Hybridization::S,
            Hybridization::Sp,
            Hybridization::Sp2,
            Hybridization::Sp3,
            Hybridization::Sp2d,
            Hybridization::Sp3d,
            Hybridization::Sp3d2,
            Hybridization::Other,
        ];
        let params = atomic_params(4.0);

        for hyb2 in values {
            for hyb3 in values {
                let result = calc_torsion_params(1.0, 6, 6, hyb2, hyb3, &params, &params, false);
                let source_accepts = matches!(hyb2, Hybridization::Sp2 | Hybridization::Sp3)
                    && matches!(hyb3, Hybridization::Sp2 | Hybridization::Sp3);

                assert_eq!(result.is_ok(), source_accepts, "{hyb2:?}/{hyb3:?}");
                if let Err(error) = result {
                    assert_eq!(error, TorsionParamsError::BadHybridizations);
                    assert_eq!(error.source_category(), "Pre-condition Violation");
                    assert_eq!(error.source_message(), "bad hybridizations");
                }
            }
        }
    }

    #[test]
    fn cf3d_u11_sp3_sp3_keeps_general_v1_formula_and_exact_single_guard() {
        // RDKit TorsionAngle.cpp uses sqrt(V1_2*V1_3) unless the exact
        // single-bond/two-group-six condition applies.
        let at2 = torsion_params(4.0, 4.0);
        let at3 = torsion_params(9.0, 9.0);

        for bond_order in [2.0, 0.999_999] {
            assert_eq!(
                calc_torsion_params(
                    bond_order,
                    8,
                    16,
                    Hybridization::Sp3,
                    Hybridization::Sp3,
                    &at2,
                    &at3,
                    true,
                ),
                Ok((6.0, 3, -1.0))
            );
        }
        for (at_num2, at_num3) in [(6, 7), (8, 6), (6, 16)] {
            assert_eq!(
                calc_torsion_params(
                    1.0,
                    at_num2,
                    at_num3,
                    Hybridization::Sp3,
                    Hybridization::Sp3,
                    &at2,
                    &at3,
                    false,
                ),
                Ok((6.0, 3, -1.0)),
                "non-pair {at_num2}/{at_num3}"
            );
        }

        // The source has no guard around sqrt(V1_2*V1_3).
        let negative_v1 = torsion_params(-4.0, 4.0);
        let result = calc_torsion_params(
            2.0,
            6,
            6,
            Hybridization::Sp3,
            Hybridization::Sp3,
            &negative_v1,
            &at3,
            false,
        )
        .unwrap();
        assert!(result.0.is_nan());
        assert_eq!((result.1, result.2), (3, -1.0));
    }

    #[test]
    fn cf3d_u11_sp3_sp3_single_group_six_override_covers_every_pair() {
        // RDKit's per-side oxygen override is independent for all 25 ordered
        // pairs of the five source group-six identities.
        let groups = [8, 16, 34, 52, 84];
        let at2 = torsion_params(4.0, 4.0);
        let at3 = torsion_params(9.0, 9.0);

        for at_num2 in groups {
            for at_num3 in groups {
                let expected_force = match (at_num2 == 8, at_num3 == 8) {
                    (true, true) => 2.0,
                    (true, false) | (false, true) => 3.687817782917155,
                    (false, false) => 6.8,
                };
                assert_eq!(
                    calc_torsion_params(
                        1.0,
                        at_num2,
                        at_num3,
                        Hybridization::Sp3,
                        Hybridization::Sp3,
                        &at2,
                        &at3,
                        false,
                    ),
                    Ok((expected_force, 2, -1.0)),
                    "{at_num2}/{at_num3}"
                );
            }
        }
    }

    #[test]
    fn cf3d_u11_sp2_sp2_uses_equation17_independent_of_bond_order() {
        // RDKit takes the SP2/SP2 branch before the mixed single-bond guard.
        let at2 = atomic_params(4.0);
        let at3 = atomic_params(9.0);

        assert_eq!(
            calc_torsion_params(
                2.0,
                6,
                7,
                Hybridization::Sp2,
                Hybridization::Sp2,
                &at2,
                &at3,
                true,
            ),
            Ok((116.92065644221714, 2, 1.0))
        );
        assert_eq!(
            calc_torsion_params(
                1.0,
                6,
                7,
                Hybridization::Sp2,
                Hybridization::Sp2,
                &at2,
                &at3,
                false,
            ),
            Ok((30.0, 2, 1.0))
        );
    }

    #[test]
    fn cf3d_u11_mixed_group_six_sp3_branch_precedes_end_atom_override() {
        // RDKit's first mixed single-bond special case wins even when
        // endAtomIsSP2 is true, in both argument orientations and for every
        // group-six identity.
        let groups = [8, 16, 34, 52, 84];
        let at2 = atomic_params(4.0);
        let at3 = atomic_params(9.0);

        for group_atom in groups {
            assert_eq!(
                calc_torsion_params(
                    1.0,
                    group_atom,
                    6,
                    Hybridization::Sp3,
                    Hybridization::Sp2,
                    &at2,
                    &at3,
                    true,
                ),
                Ok((30.0, 2, -1.0)),
                "SP3 group-six first: {group_atom}"
            );
            assert_eq!(
                calc_torsion_params(
                    1.0,
                    6,
                    group_atom,
                    Hybridization::Sp2,
                    Hybridization::Sp3,
                    &at2,
                    &at3,
                    true,
                ),
                Ok((30.0, 2, -1.0)),
                "SP3 group-six second: {group_atom}"
            );
        }
    }

    #[test]
    fn cf3d_u11_mixed_defaults_and_end_atom_conjugation_are_source_ordered() {
        // Only an exact single bond reaches either mixed special case. The
        // conjugation branch applies when the prior group-six branch is false.
        let at2 = atomic_params(4.0);
        let at3 = atomic_params(9.0);
        let hyb2 = Hybridization::Sp2;
        let hyb3 = Hybridization::Sp3;

        for bond_order in [2.0, 0.999_999, f64::NAN] {
            assert_eq!(
                calc_torsion_params(bond_order, 6, 7, hyb2, hyb3, &at2, &at3, true),
                Ok((1.0, 6, 1.0))
            );
        }
        assert_eq!(
            calc_torsion_params(1.0, 6, 7, hyb2, hyb3, &at2, &at3, false),
            Ok((1.0, 6, 1.0))
        );
        assert_eq!(
            calc_torsion_params(1.0, 6, 7, hyb2, hyb3, &at2, &at3, true),
            Ok((2.0, 3, -1.0))
        );
        assert_eq!(
            calc_torsion_params(
                2.0,
                16,
                6,
                Hybridization::Sp3,
                Hybridization::Sp2,
                &at2,
                &at3,
                true,
            ),
            Ok((1.0, 6, 1.0))
        );
        assert_eq!(
            calc_torsion_params(
                1.0,
                16,
                8,
                Hybridization::Sp3,
                Hybridization::Sp2,
                &at2,
                &at3,
                true,
            ),
            Ok((2.0, 3, -1.0))
        );
        assert_eq!(
            calc_torsion_params(
                1.0,
                16,
                6,
                Hybridization::Sp2,
                Hybridization::Sp3,
                &at2,
                &at3,
                true,
            ),
            Ok((2.0, 3, -1.0))
        );
    }

    #[test]
    fn cf3d_u12_constructor_phase_and_energy_cover_source_periodicities() {
        // RDKit TorsionAngle.cpp::calcTorsionParams selects all four reachable
        // order/phase pairs; getEnergy evaluates its source cosine polynomial.
        let cases = [
            (
                make_contribution(
                    [0, 1, 2, 3],
                    1.0,
                    6,
                    7,
                    Hybridization::Sp2,
                    Hybridization::Sp2,
                    false,
                )
                .unwrap(),
                2,
                1.0,
                30.0,
                0.0,
                PI / 2.0,
                PI,
                PI / 5.0,
            ),
            (
                make_contribution(
                    [0, 1, 2, 3],
                    1.0,
                    8,
                    16,
                    Hybridization::Sp3,
                    Hybridization::Sp3,
                    false,
                )
                .unwrap(),
                2,
                -1.0,
                3.687_817_782_917_155,
                PI / 2.0,
                0.0,
                PI,
                PI / 4.0,
            ),
            (
                make_contribution(
                    [0, 1, 2, 3],
                    2.0,
                    6,
                    7,
                    Hybridization::Sp3,
                    Hybridization::Sp3,
                    false,
                )
                .unwrap(),
                3,
                -1.0,
                6.0,
                PI / 3.0,
                0.0,
                2.0 * PI / 3.0,
                PI / 5.0,
            ),
            (
                make_contribution(
                    [0, 1, 2, 3],
                    2.0,
                    6,
                    7,
                    Hybridization::Sp2,
                    Hybridization::Sp3,
                    false,
                )
                .unwrap(),
                6,
                1.0,
                1.0,
                0.0,
                PI / 6.0,
                PI / 3.0,
                PI / 12.0,
            ),
        ];

        for (
            contribution,
            expected_order,
            expected_cos_term,
            expected_force_constant,
            equilibrium,
            maximum,
            period,
            sample,
        ) in cases
        {
            assert_eq!(contribution.order, expected_order);
            assert_eq!(contribution.cos_term, expected_cos_term);
            assert_close(contribution.force_constant, expected_force_constant);
            assert_close(energy_at(&contribution, equilibrium), 0.0);
            assert_close(energy_at(&contribution, -equilibrium), 0.0);
            assert_close(energy_at(&contribution, maximum), expected_force_constant);
            assert_close(
                energy_at(&contribution, maximum + period),
                expected_force_constant,
            );
            assert_close(
                energy_at(&contribution, sample),
                energy_at(&contribution, -sample),
            );
            assert_close(
                energy_at(&contribution, equilibrium),
                energy_at(&contribution, equilibrium + period),
            );
        }
    }

    #[test]
    fn cf3d_u12_theta_derivative_uses_each_source_polynomial_and_factor_order() {
        // RDKit TorsionAngle.cpp::getThetaDeriv applies its order polynomial,
        // then forceConstant / 2 * cosTerm * -1 * order in this sequence.
        let cases = [
            (
                make_contribution(
                    [0, 1, 2, 3],
                    1.0,
                    6,
                    7,
                    Hybridization::Sp2,
                    Hybridization::Sp2,
                    false,
                )
                .unwrap(),
                -28.8,
            ),
            (
                make_contribution(
                    [0, 1, 2, 3],
                    1.0,
                    8,
                    16,
                    Hybridization::Sp3,
                    Hybridization::Sp3,
                    false,
                )
                .unwrap(),
                3.540_305_071_600_468_8,
            ),
            (
                make_contribution(
                    [0, 1, 2, 3],
                    2.0,
                    6,
                    7,
                    Hybridization::Sp3,
                    Hybridization::Sp3,
                    false,
                )
                .unwrap(),
                3.168,
            ),
            (
                make_contribution(
                    [0, 1, 2, 3],
                    2.0,
                    6,
                    7,
                    Hybridization::Sp2,
                    Hybridization::Sp3,
                    false,
                )
                .unwrap(),
                1.976_832,
            ),
        ];

        for (contribution, expected) in cases {
            assert_close(contribution.get_theta_deriv(0.6, 0.8).unwrap(), expected);
        }
    }

    #[test]
    fn cf3d_u12_constructor_keeps_duplicate_range_and_parameter_error_order() {
        // RDKit TorsionAngle.cpp checks pairwise-distinct IDs before idx1..idx4
        // range checks, then calcTorsionParams' hybridization precondition.
        let duplicate = make_contribution(
            [4, 4, 2, 3],
            1.0,
            6,
            7,
            Hybridization::Unspecified,
            Hybridization::Sp2,
            false,
        )
        .err()
        .expect("duplicate source indices must be rejected");
        assert_eq!(duplicate, TorsionContributionError::DegeneratePoints);
        assert_eq!(duplicate.source_category(), "Pre-condition Violation");
        assert_eq!(duplicate.source_message(), "degenerate points");
        assert_eq!(duplicate.range_detail(), None);

        for (indices, argument, message) in [
            ([4, 1, 2, 3], TorsionIndexArgument::First, "idx1"),
            ([0, 4, 2, 3], TorsionIndexArgument::Second, "idx2"),
            ([0, 1, 4, 3], TorsionIndexArgument::Third, "idx3"),
            ([0, 1, 2, 4], TorsionIndexArgument::Fourth, "idx4"),
        ] {
            let error = make_contribution(
                indices,
                1.0,
                6,
                7,
                Hybridization::Unspecified,
                Hybridization::Sp2,
                false,
            )
            .err()
            .expect("out-of-range source index must be rejected");
            assert_eq!(
                error,
                TorsionContributionError::IndexOutOfRange {
                    argument,
                    index: 4,
                    upper_bound: 4,
                }
            );
            assert_eq!(error.source_category(), "Range Error");
            assert_eq!(error.source_message(), message);
            assert_eq!(error.range_detail(), Some((4, 4)));
        }

        let invalid_hybridizations = make_contribution(
            [0, 1, 2, 3],
            1.0,
            6,
            7,
            Hybridization::Unspecified,
            Hybridization::Sp2,
            false,
        )
        .err()
        .expect("unsupported source hybridization must be rejected");
        assert_eq!(
            invalid_hybridizations,
            TorsionContributionError::BadHybridizations
        );
        assert_eq!(
            invalid_hybridizations.source_category(),
            "Pre-condition Violation"
        );
        assert_eq!(
            invalid_hybridizations.source_message(),
            "bad hybridizations"
        );
    }

    #[test]
    fn cf3d_u12_scale_force_constant_preserves_unsigned_division_edges() {
        // RDKit TorsionAngle.h::scaleForceConstant divides by the unsigned
        // count after conversion to double and has no zero-count guard.
        let mut scaled = make_contribution(
            [0, 1, 2, 3],
            1.0,
            6,
            7,
            Hybridization::Sp2,
            Hybridization::Sp2,
            false,
        )
        .unwrap();
        scaled.scale_force_constant(3);
        assert_close(scaled.force_constant, 10.0);

        let mut zero_count = make_contribution(
            [0, 1, 2, 3],
            1.0,
            6,
            7,
            Hybridization::Sp2,
            Hybridization::Sp2,
            false,
        )
        .unwrap();
        zero_count.scale_force_constant(0);
        assert_eq!(zero_count.force_constant, f64::INFINITY);
    }

    #[test]
    fn cf3d_u12_invalid_order_keeps_source_precondition_error() {
        // RDKit TorsionAngle.cpp guards both getEnergy and getThetaDeriv even
        // though calcTorsionParams only constructs the three valid orders.
        let contribution = TorsionAngleContrib {
            at1_idx: 0,
            at2_idx: 1,
            at3_idx: 2,
            at4_idx: 3,
            order: 4,
            force_constant: 1.0,
            cos_term: 1.0,
        };
        let expected = Err(TorsionContributionError::BadOrder { order: 4 });
        assert_eq!(contribution.get_theta_deriv(0.6, 0.8), expected);

        let coordinates = torsion_coordinates(0.0);
        let mut distance_cache = vec![0.0; 10];
        let mut context = EvaluationContext::for_test(&coordinates, &mut distance_cache, 4);
        let error = contribution.get_energy(&mut context).unwrap_err();
        assert_eq!(error, TorsionContributionError::BadOrder { order: 4 });
        assert_eq!(error.source_category(), "Pre-condition Violation");
        assert_eq!(error.source_message(), "bad order");
    }

    #[test]
    fn cf3d_u13_get_grad_covers_all_source_orders_at_nonplanar_geometry() {
        // RDKit TorsionAngle.cpp::getGrad reaches getThetaDeriv for each
        // constructor-supported order before applying calcTorsionGrad.
        let root_half = std::f64::consts::FRAC_1_SQRT_2;
        let coordinates = [
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, root_half, root_half,
        ];
        let cases = [
            (
                1.0,
                Hybridization::Sp2,
                Hybridization::Sp2,
                2,
                [
                    0.0,
                    0.0,
                    -30.0,
                    0.0,
                    0.0,
                    30.0,
                    0.0,
                    21.213_203_435_596_423,
                    -21.213_203_435_596_43,
                    0.0,
                    -21.213_203_435_596_423,
                    21.213_203_435_596_43,
                ],
            ),
            (
                2.0,
                Hybridization::Sp3,
                Hybridization::Sp3,
                3,
                [
                    0.0,
                    0.0,
                    6.363_961_030_678_93,
                    0.0,
                    0.0,
                    -6.363_961_030_678_93,
                    0.0,
                    -4.5,
                    4.5,
                    0.0,
                    4.5,
                    -4.5,
                ],
            ),
            (
                2.0,
                Hybridization::Sp2,
                Hybridization::Sp3,
                6,
                [
                    0.0,
                    0.0,
                    3.0,
                    0.0,
                    0.0,
                    -3.0,
                    0.0,
                    -2.121_320_343_559_642_4,
                    2.121_320_343_559_643_3,
                    0.0,
                    2.121_320_343_559_642_4,
                    -2.121_320_343_559_643_3,
                ],
            ),
        ];

        for (bond_order, hyb2, hyb3, expected_order, expected) in cases {
            let contribution = make_contribution([0, 1, 2, 3], bond_order, 6, 6, hyb2, hyb3, false)
                .expect("source-reachable torsion parameters must construct");
            assert_eq!(contribution.order, expected_order);
            let actual = gradient_for_test(&contribution, &coordinates, [0.0; 12])
                .expect("source getGrad accepts the nonplanar geometry");
            assert_gradient_close(actual, expected);
        }
    }

    #[test]
    fn cf3d_u13_get_grad_planar_zero_derivative_preserves_existing_gradient() {
        // With coplanar points RDKit's cosine is 1, sinPhi is 0, and the
        // source order-2 theta derivative leaves the accumulated gradient
        // unchanged.
        let contribution = make_contribution(
            [0, 1, 2, 3],
            1.0,
            6,
            6,
            Hybridization::Sp2,
            Hybridization::Sp2,
            false,
        )
        .expect("source-reachable torsion parameters must construct");
        let coordinates = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0];
        let initial = [1.25; 12];

        let actual = gradient_for_test(&contribution, &coordinates, initial)
            .expect("source getGrad accepts planar points");

        assert_eq!(actual, initial);
    }

    #[test]
    fn cf3d_u13_get_grad_degenerate_normal_uses_source_denominator_floor() {
        // RDKit computeDihedral floors a zero normal length at 1.0e-5 and
        // getGrad continues into the analytic formula without rejecting it.
        let contribution = make_contribution(
            [0, 1, 2, 3],
            2.0,
            6,
            6,
            Hybridization::Sp3,
            Hybridization::Sp3,
            false,
        )
        .expect("source-reachable torsion parameters must construct");
        let coordinates = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0];

        let actual = gradient_for_test(&contribution, &coordinates, [0.0; 12])
            .expect("source getGrad does not reject a degenerate normal");

        assert_gradient_close(
            actual,
            [
                0.0, -900_000.0, 0.0, 0.0, 900_000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            ],
        );
    }

    #[test]
    fn cf3d_u13_get_grad_accumulates_all_four_source_atom_vectors() {
        // Each source atom's three components use +=; no component is reset.
        let contribution = make_contribution(
            [0, 1, 2, 3],
            1.0,
            6,
            6,
            Hybridization::Sp2,
            Hybridization::Sp2,
            false,
        )
        .expect("source-reachable torsion parameters must construct");
        let root_three_over_two = 3.0_f64.sqrt() / 2.0;
        let coordinates = [
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
            root_three_over_two,
        ];

        let actual = gradient_for_test(&contribution, &coordinates, [1.25; 12])
            .expect("source getGrad accepts nonplanar points");

        assert_gradient_close(
            actual,
            [
                1.25,
                1.25,
                -24.730_762_113_533_17,
                1.25,
                1.25,
                27.230_762_113_533_17,
                1.25,
                23.75,
                -11.740_381_056_766_59,
                1.25,
                -21.25,
                14.240_381_056_766_59,
            ],
        );
    }

    #[test]
    fn cf3d_u13_calc_torsion_grad_matches_fixed_legacy_source_vector() {
        // Fixed legacy/source helper case: nontrivial r/t/d, cosPhi=0.25,
        // sinTerm=2, and a nonzero initial gradient expose every += update.
        let r = [
            point(1.0, 2.0, 3.0),
            point(4.0, 5.0, 6.0),
            point(7.0, 8.0, 9.0),
            point(10.0, 11.0, 12.0),
        ];
        let t = [point(1.0, 0.0, 0.0), point(0.0, 1.0, 0.0)];
        let d = [2.0, 4.0];
        let mut gradient = [1.0; 12];

        calc_torsion_grad(&r, &t, &d, &mut gradient, [0, 1, 2, 3], 2.0, 0.25);

        assert_gradient_close(
            gradient,
            [
                -5.0, -0.5, 6.25, 5.5, 7.75, -9.5, 3.625, 0.25, 1.375, -0.125, -3.5, 5.875,
            ],
        );
    }

    #[test]
    fn cf3d_bld_b06_field_energy_scale_and_copy_use_the_real_kernel() {
        // RDKit UFF/TorsionAngle.cpp::getEnergy gives the fixed order-2 maximum
        // at phi=pi/2; the builder scales this stored force constant in place.
        let coordinates = torsion_coordinates(PI / 2.0);
        let mut rows = rows_from_coordinates(&coordinates);
        let unscaled = make_contribution(
            [0, 1, 2, 3],
            1.0,
            6,
            7,
            Hybridization::Sp2,
            Hybridization::Sp2,
            false,
        )
        .expect("source order-2 term constructs");
        let mut unscaled_field = field_with_torsion(&mut rows, unscaled);
        assert_close(
            cf3d_bld_b05_calc_energy(&mut unscaled_field, &coordinates)
                .expect("source field energy succeeds"),
            30.0,
        );

        let mut scaled = make_contribution(
            [0, 1, 2, 3],
            1.0,
            6,
            7,
            Hybridization::Sp2,
            Hybridization::Sp2,
            false,
        )
        .expect("source order-2 term constructs");
        scaled.scale_force_constant(2);
        let mut scaled_rows = rows_from_coordinates(&coordinates);
        let mut scaled_field = field_with_torsion(&mut scaled_rows, scaled);
        assert_close(
            cf3d_bld_b05_calc_energy(&mut scaled_field, &coordinates)
                .expect("scaled source field energy succeeds"),
            15.0,
        );

        // ForceField::copy invokes this contribution's actual trait copy and
        // starts without borrowed positions; attach the copy's own row storage.
        let mut copied_field = cf3d_bld_b05_copy_force_field(&scaled_field);
        let mut copied_rows = rows_from_coordinates(&coordinates);
        for row in &mut copied_rows {
            copied_field.positions_mut().push(row.as_mut_slice());
        }
        copied_field
            .initialize()
            .expect("copied field initializes with its own positions");
        assert_close(
            cf3d_bld_b05_calc_energy(&mut copied_field, &coordinates)
                .expect("copied scaled source field energy succeeds"),
            15.0,
        );
    }

    #[test]
    fn cf3d_bld_b06_field_gradient_keeps_source_values_and_additive_writes() {
        // Fixed RDKit U13 order-2 gradient at phi=pi/4, dispatched through
        // ForceField::calcGrad with a nonzero caller gradient.
        let root_half = std::f64::consts::FRAC_1_SQRT_2;
        let coordinates = [
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, root_half, root_half,
        ];
        let mut rows = rows_from_coordinates(&coordinates);
        let contribution = make_contribution(
            [0, 1, 2, 3],
            1.0,
            6,
            7,
            Hybridization::Sp2,
            Hybridization::Sp2,
            false,
        )
        .expect("source order-2 term constructs");
        let mut field = field_with_torsion(&mut rows, contribution);
        let mut gradient = [1.25; 12];

        cf3d_bld_b05_calc_grad(&mut field, &coordinates, &mut gradient)
            .expect("source field gradient succeeds");

        assert_gradient_close(
            gradient,
            [
                1.25,
                1.25,
                -28.75,
                1.25,
                1.25,
                31.25,
                1.25,
                22.463_203_435_596_423,
                -19.963_203_435_596_43,
                1.25,
                -19.963_203_435_596_423,
                22.463_203_435_596_43,
            ],
        );
    }

    #[test]
    fn cf3d_bld_b06_field_energy_preserves_zero_and_nan_geometry() {
        // calculateCosTorsion returns literal zero for zero normals, and its
        // source clamp leaves NaN unchanged for a NaN normal.
        let zero_coordinates = [0.0; 12];
        let mut zero_rows = rows_from_coordinates(&zero_coordinates);
        let zero_term = make_contribution(
            [0, 1, 2, 3],
            1.0,
            6,
            7,
            Hybridization::Sp2,
            Hybridization::Sp2,
            false,
        )
        .expect("source order-2 term constructs");
        let mut zero_field = field_with_torsion(&mut zero_rows, zero_term);
        assert_eq!(
            cf3d_bld_b05_calc_energy(&mut zero_field, &zero_coordinates),
            Ok(30.0)
        );

        let nan_coordinates = [
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
            f64::NAN,
            0.0,
        ];
        let mut nan_rows = rows_from_coordinates(&nan_coordinates);
        let nan_term = make_contribution(
            [0, 1, 2, 3],
            1.0,
            6,
            7,
            Hybridization::Sp2,
            Hybridization::Sp2,
            false,
        )
        .expect("source order-2 term constructs");
        let mut nan_field = field_with_torsion(&mut nan_rows, nan_term);
        assert!(
            cf3d_bld_b05_calc_energy(&mut nan_field, &nan_coordinates)
                .expect("source field energy propagates NaN as a value")
                .is_nan()
        );
    }

    #[test]
    fn cf3d_bld_b06_field_errors_keep_source_mapping_and_call_order() {
        // Preserve every contribution error's typed source identity at the
        // kernel boundary, including all four range arguments.
        for (source_error, kernel_error) in [
            (
                TorsionContributionError::DegeneratePoints,
                ForceFieldKernelError::TorsionDegeneratePoints,
            ),
            (
                TorsionContributionError::IndexOutOfRange {
                    argument: TorsionIndexArgument::First,
                    index: 8,
                    upper_bound: 4,
                },
                ForceFieldKernelError::TorsionIndexOutOfRange {
                    argument: TorsionIndexArgument::First,
                    index: 8,
                    upper_bound: 4,
                },
            ),
            (
                TorsionContributionError::IndexOutOfRange {
                    argument: TorsionIndexArgument::Second,
                    index: 8,
                    upper_bound: 4,
                },
                ForceFieldKernelError::TorsionIndexOutOfRange {
                    argument: TorsionIndexArgument::Second,
                    index: 8,
                    upper_bound: 4,
                },
            ),
            (
                TorsionContributionError::IndexOutOfRange {
                    argument: TorsionIndexArgument::Third,
                    index: 8,
                    upper_bound: 4,
                },
                ForceFieldKernelError::TorsionIndexOutOfRange {
                    argument: TorsionIndexArgument::Third,
                    index: 8,
                    upper_bound: 4,
                },
            ),
            (
                TorsionContributionError::IndexOutOfRange {
                    argument: TorsionIndexArgument::Fourth,
                    index: 8,
                    upper_bound: 4,
                },
                ForceFieldKernelError::TorsionIndexOutOfRange {
                    argument: TorsionIndexArgument::Fourth,
                    index: 8,
                    upper_bound: 4,
                },
            ),
            (
                TorsionContributionError::BadHybridizations,
                ForceFieldKernelError::TorsionBadHybridizations,
            ),
            (
                TorsionContributionError::BadOrder { order: 4 },
                ForceFieldKernelError::TorsionBadOrder { order: 4 },
            ),
        ] {
            assert_eq!(ForceFieldKernelError::from(source_error), kernel_error);
        }

        let coordinates = torsion_coordinates(PI / 4.0);
        let mut rows = rows_from_coordinates(&coordinates);
        let invalid_order = TorsionAngleContrib {
            at1_idx: 0,
            at2_idx: 1,
            at3_idx: 2,
            at4_idx: 3,
            order: 4,
            force_constant: 30.0,
            cos_term: 1.0,
        };
        let mut field = ForceField::new(3);
        for row in &mut rows {
            field.positions_mut().push(row.as_mut_slice());
        }
        field.add_contribution(Box::new(invalid_order));
        let later_calls = Arc::new(AtomicUsize::new(0));
        field.add_contribution(Box::new(CountContributionCalls(Arc::clone(&later_calls))));
        field.initialize().expect("four 3D positions initialize");

        // getEnergy's order precondition precedes any coordinate read; the
        // source ForceField also stops at this first callback error.
        assert_eq!(
            cf3d_bld_b05_calc_energy(&mut field, &[]),
            Err(ForceFieldKernelError::TorsionBadOrder { order: 4 })
        );
        assert_eq!(later_calls.load(Ordering::Relaxed), 0);

        // getGrad computes source geometry before getThetaDeriv rejects order,
        // and no component is written or later contribution called.
        let mut gradient = [2.5; 12];
        assert_eq!(
            cf3d_bld_b05_calc_grad(&mut field, &coordinates, &mut gradient),
            Err(ForceFieldKernelError::TorsionBadOrder { order: 4 })
        );
        assert_eq!(gradient, [2.5; 12]);
        assert_eq!(later_calls.load(Ordering::Relaxed), 0);
    }
}
