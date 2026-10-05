//! Complete CrystalFF six-term torsion collection and single contribution.
use crate::geometry::Point3;
use crate::kernel::{EvaluationContext, ForceField, ForceFieldContribution, ForceFieldKernelError};
use crate::mmff::{calc_torsion_cos_phi, calc_torsion_grad};
fn point3(x: f64, y: f64, z: f64) -> Point3 {
    Point3 { x, y, z }
}
fn is_double_zero(value: f64) -> bool {
    // BEGIN RDKIT CPP HELPER ForceFields::MMFF::isDoubleZero (Params.h via TorsionAngleContribs.cpp)
    // RDKit✔️✔️: inline bool isDoubleZero(const double x) {
    // RDKit✔️✔️:   return ((x < 1.0e-10) && (x > -1.0e-10));
    // RDKit✔️✔️: }
    value < 1.0e-10 && value > -1.0e-10
}

#[derive(Clone, Debug, PartialEq)]
pub struct TorsionAngleContribsParams {
    pub idx1: usize,
    pub idx2: usize,
    pub idx3: usize,
    pub idx4: usize,
    pub force_constants: Vec<f64>,
    pub signs: Vec<i32>,
}

impl TorsionAngleContribsParams {
    #[must_use]
    pub fn new(
        idx1: usize,
        idx2: usize,
        idx3: usize,
        idx4: usize,
        force_constants: Vec<f64>,
        signs: Vec<i32>,
    ) -> Self {
        Self {
            idx1,
            idx2,
            idx3,
            idx4,
            force_constants,
            signs,
        }
    }
}

#[must_use]
pub fn calc_torsion_energy_m6(force_constants: &[f64], signs: &[i32], cos_phi: f64) -> f64 {
    // BEGIN COMPLETE PINNED CPP TorsionAngleContribs.cpp:23-43
    // RDKit❗✔️: double calcTorsionEnergyM6(const std::vector<double> &forceConstants,
    // RDKit❗✔️:                            const std::vector<int> &signs, const double cosPhi) {
    // RDKit❗✔️:   const double cosPhi2 = cosPhi * cosPhi;
    // RDKit❗✔️:   const double cosPhi3 = cosPhi * cosPhi2;
    // RDKit❗✔️:   const double cosPhi4 = cosPhi * cosPhi3;
    // RDKit❗✔️:   const double cosPhi5 = cosPhi * cosPhi4;
    // RDKit❗✔️:   const double cosPhi6 = cosPhi * cosPhi5;
    // RDKit❗✔️:
    // RDKit❗✔️:   const double cos2Phi = 2.0 * cosPhi2 - 1.0;
    // RDKit❗✔️:   const double cos3Phi = 4.0 * cosPhi3 - 3.0 * cosPhi;
    // RDKit❗✔️:   const double cos4Phi = 8.0 * cosPhi4 - 8.0 * cosPhi2 + 1.0;
    // RDKit❗✔️:   const double cos5Phi = 16.0 * cosPhi5 - 20.0 * cosPhi3 + 5.0 * cosPhi;
    // RDKit❗✔️:   const double cos6Phi = 32.0 * cosPhi6 - 48.0 * cosPhi4 + 18.0 * cosPhi2 - 1.0;
    // RDKit❗✔️:
    // RDKit❗✔️:   return (forceConstants[0] * (1.0 + signs[0] * cosPhi) +
    // RDKit❗✔️:           forceConstants[1] * (1.0 + signs[1] * cos2Phi) +
    // RDKit❗✔️:           forceConstants[2] * (1.0 + signs[2] * cos3Phi) +
    // RDKit❗✔️:           forceConstants[3] * (1.0 + signs[3] * cos4Phi) +
    // RDKit❗✔️:           forceConstants[4] * (1.0 + signs[4] * cos5Phi) +
    // RDKit❗✔️:           forceConstants[5] * (1.0 + signs[5] * cos6Phi));
    // RDKit❗✔️: }
    // END COMPLETE PINNED CPP TorsionAngleContribs.cpp
    // BEGIN COMPLETE PINNED CPP TorsionAngleM6.cpp:26-44
    // RDKit❗✔️: double calcTorsionEnergyM6(const std::vector<double> &V,
    // RDKit❗✔️:                            const std::vector<int> &signs, const double cosPhi) {
    // RDKit❗✔️:   double cosPhi2 = cosPhi * cosPhi;
    // RDKit❗✔️:   double cosPhi3 = cosPhi * cosPhi2;
    // RDKit❗✔️:   double cosPhi4 = cosPhi * cosPhi3;
    // RDKit❗✔️:   double cosPhi5 = cosPhi * cosPhi4;
    // RDKit❗✔️:   double cosPhi6 = cosPhi * cosPhi5;
    // RDKit❗✔️:
    // RDKit❗✔️:   double cos2Phi = 2.0 * cosPhi2 - 1.0;
    // RDKit❗✔️:   double cos3Phi = 4.0 * cosPhi3 - 3.0 * cosPhi;
    // RDKit❗✔️:   double cos4Phi = 8.0 * cosPhi4 - 8.0 * cosPhi2 + 1.0;
    // RDKit❗✔️:   double cos5Phi = 16.0 * cosPhi5 - 20.0 * cosPhi3 + 5.0 * cosPhi;
    // RDKit❗✔️:   double cos6Phi = 32.0 * cosPhi6 - 48.0 * cosPhi4 + 18.0 * cosPhi2 - 1.0;
    // RDKit❗✔️:
    // RDKit❗✔️:   return (
    // RDKit❗✔️:       V[0] * (1.0 + signs[0] * cosPhi) + V[1] * (1.0 + signs[1] * cos2Phi) +
    // RDKit❗✔️:       V[2] * (1.0 + signs[2] * cos3Phi) + V[3] * (1.0 + signs[3] * cos4Phi) +
    // RDKit❗✔️:       V[4] * (1.0 + signs[4] * cos5Phi) + V[5] * (1.0 + signs[5] * cos6Phi));
    // RDKit❗✔️: }
    // END COMPLETE PINNED CPP TorsionAngleM6.cpp

    // BEGIN RDKIT CPP FUNCTION ForceFields::CrystalFF::calcTorsionEnergyM6 (TorsionAngleContribs.cpp:20-41; TorsionAngleM6.cpp:21-43)
    // RDKit✔️✔️: double calcTorsionEnergyM6(const std::vector<double> &forceConstants,
    // RDKit✔️✔️:                            const std::vector<int> &signs, const double cosPhi) {
    // RDKit✔️✔️:   const double cosPhi2 = cosPhi * cosPhi;
    let cos_phi2 = cos_phi * cos_phi;
    // RDKit✔️✔️:   const double cosPhi3 = cosPhi * cosPhi2;
    let cos_phi3 = cos_phi * cos_phi2;
    // RDKit✔️✔️:   const double cosPhi4 = cosPhi * cosPhi3;
    let cos_phi4 = cos_phi * cos_phi3;
    // RDKit✔️✔️:   const double cosPhi5 = cosPhi * cosPhi4;
    let cos_phi5 = cos_phi * cos_phi4;
    // RDKit✔️✔️:   const double cosPhi6 = cosPhi * cosPhi5;
    let cos_phi6 = cos_phi * cos_phi5;

    // RDKit✔️✔️:   const double cos2Phi = 2.0 * cosPhi2 - 1.0;
    let cos2_phi = 2.0 * cos_phi2 - 1.0;
    // RDKit✔️✔️:   const double cos3Phi = 4.0 * cosPhi3 - 3.0 * cosPhi;
    let cos3_phi = 4.0 * cos_phi3 - 3.0 * cos_phi;
    // RDKit✔️✔️:   const double cos4Phi = 8.0 * cosPhi4 - 8.0 * cosPhi2 + 1.0;
    let cos4_phi = 8.0 * cos_phi4 - 8.0 * cos_phi2 + 1.0;
    // RDKit✔️✔️:   const double cos5Phi = 16.0 * cosPhi5 - 20.0 * cosPhi3 + 5.0 * cosPhi;
    let cos5_phi = 16.0 * cos_phi5 - 20.0 * cos_phi3 + 5.0 * cos_phi;
    // RDKit✔️✔️:   const double cos6Phi = 32.0 * cosPhi6 - 48.0 * cosPhi4 + 18.0 * cosPhi2 - 1.0;
    let cos6_phi = 32.0 * cos_phi6 - 48.0 * cos_phi4 + 18.0 * cos_phi2 - 1.0;

    // RDKit✔️✔️:   return (forceConstants[0] * (1.0 + signs[0] * cosPhi) +
    // RDKit✔️✔️:           forceConstants[1] * (1.0 + signs[1] * cos2Phi) +
    // RDKit✔️✔️:           forceConstants[2] * (1.0 + signs[2] * cos3Phi) +
    // RDKit✔️✔️:           forceConstants[3] * (1.0 + signs[3] * cos4Phi) +
    // RDKit✔️✔️:           forceConstants[4] * (1.0 + signs[4] * cos5Phi) +
    // RDKit✔️✔️:           forceConstants[5] * (1.0 + signs[5] * cos6Phi));
    // RDKit✔️✔️: }
    force_constants[0] * (1.0 + f64::from(signs[0]) * cos_phi)
        + force_constants[1] * (1.0 + f64::from(signs[1]) * cos2_phi)
        + force_constants[2] * (1.0 + f64::from(signs[2]) * cos3_phi)
        + force_constants[3] * (1.0 + f64::from(signs[3]) * cos4_phi)
        + force_constants[4] * (1.0 + f64::from(signs[4]) * cos5_phi)
        + force_constants[5] * (1.0 + f64::from(signs[5]) * cos6_phi)
}

#[must_use]
pub fn calc_torsion_energy(force_constants: &[f64], signs: &[i32], cos_phi: f64) -> f64 {
    calc_torsion_energy_m6(force_constants, signs, cos_phi)
}

#[derive(Clone, Debug)]
pub struct TorsionAngleContribs {
    owner_points: Option<usize>,
    contribs: Vec<TorsionAngleContribsParams>,
}

impl TorsionAngleContribs {
    #[must_use]
    pub fn new(owner: &ForceField<'_>) -> Self {
        // BEGIN COMPLETE PINNED CPP TorsionAngleContribs.cpp:46-49
        // RDKit❗✔️: TorsionAngleContribs::TorsionAngleContribs(ForceField *owner) {
        // RDKit❗✔️:   PRECONDITION(owner, "bad owner");
        // RDKit❗✔️:   dp_forceField = owner;
        // RDKit❗✔️: }
        // END COMPLETE PINNED CPP TorsionAngleContribs.cpp

        // BEGIN RDKIT CPP CONSTRUCTOR ForceFields::CrystalFF::TorsionAngleContribs::TorsionAngleContribs (TorsionAngleContribs.cpp:44-47)
        // RDKit✔️✔️: TorsionAngleContribs::TorsionAngleContribs(ForceField *owner) {
        // RDKit✔️✔️:   PRECONDITION(owner, "bad owner");
        // Rust references reproduce RDKit's non-null owner precondition.
        // RDKit✔️✔️:   dp_forceField = owner;
        // RDKit✔️✔️: }
        Self {
            owner_points: Some(owner.positions().len()),
            contribs: Vec::new(),
        }
    }

    pub fn add_contrib(
        &mut self,
        idx1: usize,
        idx2: usize,
        idx3: usize,
        idx4: usize,
        force_constants: Vec<f64>,
        signs: Vec<i32>,
    ) {
        // BEGIN COMPLETE PINNED CPP TorsionAngleContribs::addContrib (TorsionAngleContribs.cpp:51-64)
        // RDKit❗✔️: void TorsionAngleContribs::addContrib(unsigned int idx1, unsigned int idx2,
        // RDKit❗✔️:                                       unsigned int idx3, unsigned int idx4,
        // RDKit❗✔️:                                       std::vector<double> forceConstants,
        // RDKit❗✔️:                                       std::vector<int> signs) {
        // RDKit❗✔️:   PRECONDITION((idx1 != idx2) && (idx1 != idx3) && (idx1 != idx4) &&
        // RDKit❗✔️:                    (idx2 != idx3) && (idx2 != idx4) && (idx3 != idx4),
        // RDKit❗✔️:                "degenerate points");
        // RDKit❗✔️:   URANGE_CHECK(idx1, dp_forceField->positions().size());
        // RDKit❗✔️:   URANGE_CHECK(idx2, dp_forceField->positions().size());
        // RDKit❗✔️:   URANGE_CHECK(idx3, dp_forceField->positions().size());
        // RDKit❗✔️:   URANGE_CHECK(idx4, dp_forceField->positions().size());
        // RDKit❗✔️:   d_contribs.emplace_back(idx1, idx2, idx3, idx4, std::move(forceConstants),
        // RDKit❗✔️:                           std::move(signs));
        // RDKit❗✔️: }
        // END COMPLETE PINNED CPP TorsionAngleContribs::addContrib

        // BEGIN RDKIT CPP METHOD ForceFields::CrystalFF::TorsionAngleContribs::addContrib (TorsionAngleContribs.cpp:49-61)
        // RDKit✔️✔️: PRECONDITION((idx1 != idx2) && (idx1 != idx3) && (idx1 != idx4) &&
        // RDKit✔️✔️:                (idx2 != idx3) && (idx2 != idx4) && (idx3 != idx4),
        // RDKit✔️✔️:              "degenerate points");
        validate_indices(
            self.owner_points.expect("bad owner"),
            idx1,
            idx2,
            idx3,
            idx4,
        );
        // RDKit✔️✔️: d_contribs.emplace_back(idx1, idx2, idx3, idx4, std::move(forceConstants),
        // RDKit✔️✔️:                         std::move(signs));
        self.contribs.push(TorsionAngleContribsParams::new(
            idx1,
            idx2,
            idx3,
            idx4,
            force_constants,
            signs,
        ));
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.contribs.is_empty()
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.contribs.len()
    }
}

impl ForceFieldContribution for TorsionAngleContribs {
    fn copy(&self) -> Box<dyn ForceFieldContribution> {
        Box::new(self.clone())
    }

    fn get_energy(
        &self,
        context: &mut EvaluationContext<'_>,
    ) -> Result<f64, ForceFieldKernelError> {
        energy_terms(self.owner_points, &self.contribs, context)
    }

    fn get_grad(
        &self,
        context: &mut EvaluationContext<'_>,
        grad: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        gradient_terms(self.owner_points, &self.contribs, context, grad)
    }
}

fn energy_terms(
    owner_points: Option<usize>,
    terms: &[TorsionAngleContribsParams],
    context: &mut EvaluationContext<'_>,
) -> Result<f64, ForceFieldKernelError> {
    // BEGIN COMPLETE PINNED CPP TorsionAngleContribs::getEnergy (TorsionAngleContribs.cpp:66-88)
    // RDKit❗✔️: double TorsionAngleContribs::getEnergy(double *pos) const {
    // RDKit❗✔️:   PRECONDITION(dp_forceField, "no owner");
    // RDKit❗✔️:   PRECONDITION(pos, "bad vector");
    // RDKit❗✔️:   double accum = 0.0;
    // RDKit❗✔️:   for (const auto &contrib : d_contribs) {
    // RDKit❗✔️:     const RDGeom::Point3D iPoint(pos[3 * contrib.idx1],
    // RDKit❗✔️:                                  pos[3 * contrib.idx1 + 1],
    // RDKit❗✔️:                                  pos[3 * contrib.idx1 + 2]);
    // RDKit❗✔️:     const RDGeom::Point3D jPoint(pos[3 * contrib.idx2],
    // RDKit❗✔️:                                  pos[3 * contrib.idx2 + 1],
    // RDKit❗✔️:                                  pos[3 * contrib.idx2 + 2]);
    // RDKit❗✔️:     const RDGeom::Point3D kPoint(pos[3 * contrib.idx3],
    // RDKit❗✔️:                                  pos[3 * contrib.idx3 + 1],
    // RDKit❗✔️:                                  pos[3 * contrib.idx3 + 2]);
    // RDKit❗✔️:     const RDGeom::Point3D lPoint(pos[3 * contrib.idx4],
    // RDKit❗✔️:                                  pos[3 * contrib.idx4 + 1],
    // RDKit❗✔️:                                  pos[3 * contrib.idx4 + 2]);
    // RDKit❗✔️:     accum += calcTorsionEnergyM6(
    // RDKit❗✔️:         contrib.forceConstants, contrib.signs,
    // RDKit❗✔️:         MMFF::Utils::calcTorsionCosPhi(iPoint, jPoint, kPoint, lPoint));
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return accum;
    // RDKit❗✔️: }
    // END COMPLETE PINNED CPP TorsionAngleContribs::getEnergy

    let pos = context.coordinates();
    // BEGIN RDKIT CPP METHOD ForceFields::CrystalFF::TorsionAngleContribs::getEnergy (TorsionAngleContribs.cpp:66-88)
    // RDKit✔️✔️: PRECONDITION(dp_forceField, "no owner");
    let owner_points = owner_points.expect("no owner");
    // RDKit✔️✔️: PRECONDITION(pos, "bad vector");
    assert!(!pos.is_empty(), "bad vector");
    assert!(
        pos.len() >= 3 * owner_points,
        "bad vector length for force-field positions"
    );
    // RDKit✔️✔️: double accum = 0.0;
    let mut accum = 0.0;
    // RDKit✔️✔️: for (const auto &contrib : d_contribs) {
    for contrib in terms {
        // RDKit✔️✔️:   const RDGeom::Point3D iPoint(pos[3 * contrib.idx1],
        // RDKit✔️✔️:                              pos[3 * contrib.idx1 + 1],
        // RDKit✔️✔️:                              pos[3 * contrib.idx1 + 2]);
        // RDKit✔️✔️:   const RDGeom::Point3D jPoint(pos[3 * contrib.idx2],
        // RDKit✔️✔️:                              pos[3 * contrib.idx2 + 1],
        // RDKit✔️✔️:                              pos[3 * contrib.idx2 + 2]);
        // RDKit✔️✔️:   const RDGeom::Point3D kPoint(pos[3 * contrib.idx3],
        // RDKit✔️✔️:                              pos[3 * contrib.idx3 + 1],
        // RDKit✔️✔️:                              pos[3 * contrib.idx3 + 2]);
        // RDKit✔️✔️:   const RDGeom::Point3D lPoint(pos[3 * contrib.idx4],
        // RDKit✔️✔️:                              pos[3 * contrib.idx4 + 1],
        // RDKit✔️✔️:                              pos[3 * contrib.idx4 + 2]);
        // RDKit✔️✔️:   accum += calcTorsionEnergyM6(
        // RDKit✔️✔️:       contrib.forceConstants, contrib.signs,
        // RDKit✔️✔️:       MMFF::Utils::calcTorsionCosPhi(iPoint, jPoint, kPoint, lPoint));
        let i_point = point3(
            pos[3 * contrib.idx1],
            pos[3 * contrib.idx1 + 1],
            pos[3 * contrib.idx1 + 2],
        );
        let j_point = point3(
            pos[3 * contrib.idx2],
            pos[3 * contrib.idx2 + 1],
            pos[3 * contrib.idx2 + 2],
        );
        let k_point = point3(
            pos[3 * contrib.idx3],
            pos[3 * contrib.idx3 + 1],
            pos[3 * contrib.idx3 + 2],
        );
        let l_point = point3(
            pos[3 * contrib.idx4],
            pos[3 * contrib.idx4 + 1],
            pos[3 * contrib.idx4 + 2],
        );
        let r = [
            Point3::difference(&i_point, &j_point),
            Point3::difference(&k_point, &j_point),
            Point3::difference(&j_point, &k_point),
            Point3::difference(&l_point, &k_point),
        ];
        let t = [r[0].cross_product(&r[1]), r[2].cross_product(&r[3])];
        let d = [t[0].length(), t[1].length()];
        let cos_phi = calc_torsion_cos_phi(&i_point, &j_point, &k_point, &l_point);
        let energy = calc_torsion_energy_m6(&contrib.force_constants, &contrib.signs, cos_phi);

        accum += energy;
    }
    // RDKit✔️✔️: return accum;
    // RDKit✔️✔️: }
    Ok(accum)
}
fn gradient_terms(
    owner_points: Option<usize>,
    terms: &[TorsionAngleContribsParams],
    context: &mut EvaluationContext<'_>,
    grad: &mut [f64],
) -> Result<(), ForceFieldKernelError> {
    // BEGIN COMPLETE PINNED CPP TorsionAngleContribs::getGrad (TorsionAngleContribs.cpp:90-151)
    // RDKit❗✔️: void TorsionAngleContribs::getGrad(double *pos, double *grad) const {
    // RDKit❗✔️:   PRECONDITION(dp_forceField, "no owner");
    // RDKit❗✔️:   PRECONDITION(pos, "bad vector");
    // RDKit❗✔️:   PRECONDITION(grad, "bad vector");
    // RDKit❗✔️:
    // RDKit❗✔️:   for (const auto &contrib : d_contribs) {
    // RDKit❗✔️:     const RDGeom::Point3D iPoint(pos[3 * contrib.idx1],
    // RDKit❗✔️:                                  pos[3 * contrib.idx1 + 1],
    // RDKit❗✔️:                                  pos[3 * contrib.idx1 + 2]);
    // RDKit❗✔️:     const RDGeom::Point3D jPoint(pos[3 * contrib.idx2],
    // RDKit❗✔️:                                  pos[3 * contrib.idx2 + 1],
    // RDKit❗✔️:                                  pos[3 * contrib.idx2 + 2]);
    // RDKit❗✔️:     const RDGeom::Point3D kPoint(pos[3 * contrib.idx3],
    // RDKit❗✔️:                                  pos[3 * contrib.idx3 + 1],
    // RDKit❗✔️:                                  pos[3 * contrib.idx3 + 2]);
    // RDKit❗✔️:     const RDGeom::Point3D lPoint(pos[3 * contrib.idx4],
    // RDKit❗✔️:                                  pos[3 * contrib.idx4 + 1],
    // RDKit❗✔️:                                  pos[3 * contrib.idx4 + 2]);
    // RDKit❗✔️:     double *g[4] = {&(grad[3 * contrib.idx1]), &(grad[3 * contrib.idx2]),
    // RDKit❗✔️:                     &(grad[3 * contrib.idx3]), &(grad[3 * contrib.idx4])};
    // RDKit❗✔️:
    // RDKit❗✔️:     RDGeom::Point3D r[4] = {iPoint - jPoint, kPoint - jPoint, jPoint - kPoint,
    // RDKit❗✔️:                             lPoint - kPoint};
    // RDKit❗✔️:     RDGeom::Point3D t[2] = {r[0].crossProduct(r[1]), r[2].crossProduct(r[3])};
    // RDKit❗✔️:     double d[2] = {t[0].length(), t[1].length()};
    // RDKit❗✔️:     if (MMFF::isDoubleZero(d[0]) || MMFF::isDoubleZero(d[1])) {
    // RDKit❗✔️:       return;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     t[0] /= d[0];
    // RDKit❗✔️:     t[1] /= d[1];
    // RDKit❗✔️:     double cosPhi = t[0].dotProduct(t[1]);
    // RDKit❗✔️:     cosPhi = std::clamp(cosPhi, -1.0, 1.0);
    // RDKit❗✔️:     const double sinPhiSq = 1.0 - cosPhi * cosPhi;
    // RDKit❗✔️:     const double sinPhi = ((sinPhiSq > 0.0) ? sqrt(sinPhiSq) : 0.0);
    // RDKit❗✔️:     const double cosPhi2 = cosPhi * cosPhi;
    // RDKit❗✔️:     const double cosPhi3 = cosPhi * cosPhi2;
    // RDKit❗✔️:     const double cosPhi4 = cosPhi * cosPhi3;
    // RDKit❗✔️:     const double cosPhi5 = cosPhi * cosPhi4;
    // RDKit❗✔️:     // dE/dPhi is independent of cartesians:
    // RDKit❗✔️:     const double dE_dPhi =
    // RDKit❗✔️:         (-contrib.forceConstants[0] * contrib.signs[0] * sinPhi -
    // RDKit❗✔️:          2.0 * contrib.forceConstants[1] * contrib.signs[1] *
    // RDKit❗✔️:              (2.0 * cosPhi * sinPhi) -
    // RDKit❗✔️:          3.0 * contrib.forceConstants[2] * contrib.signs[2] *
    // RDKit❗✔️:              (4.0 * cosPhi2 * sinPhi - sinPhi) -
    // RDKit❗✔️:          4.0 * contrib.forceConstants[3] * contrib.signs[3] *
    // RDKit❗✔️:              (8.0 * cosPhi3 * sinPhi - 4.0 * cosPhi * sinPhi) -
    // RDKit❗✔️:          5.0 * contrib.forceConstants[4] * contrib.signs[4] *
    // RDKit❗✔️:              (16.0 * cosPhi4 * sinPhi - 12.0 * cosPhi2 * sinPhi + sinPhi) -
    // RDKit❗✔️:          6.0 * contrib.forceConstants[4] * contrib.signs[4] *
    // RDKit❗✔️:              (32.0 * cosPhi5 * sinPhi - 32.0 * cosPhi3 * sinPhi +
    // RDKit❗✔️:               6.0 * sinPhi));
    // RDKit❗✔️:
    // RDKit❗✔️:     // FIX: use a tolerance here
    // RDKit❗✔️:     // this is hacky, but it's per the
    // RDKit❗✔️:     // recommendation from Niketic and Rasmussen:
    // RDKit❗✔️:     double sinTerm = -dE_dPhi * (MMFF::isDoubleZero(sinPhi) ? (1.0 / cosPhi)
    // RDKit❗✔️:                                                             : (1.0 / sinPhi));
    // RDKit❗✔️:
    // RDKit❗✔️:     MMFF::Utils::calcTorsionGrad(r, t, d, g, sinTerm, cosPhi);
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END COMPLETE PINNED CPP TorsionAngleContribs::getGrad

    let pos = context.coordinates();
    // BEGIN RDKIT CPP METHOD ForceFields::CrystalFF::TorsionAngleContribs::getGrad (TorsionAngleContribs.cpp:90-150)
    // RDKit✔️✔️: PRECONDITION(dp_forceField, "no owner");
    let owner_points = owner_points.expect("no owner");
    // RDKit✔️✔️: PRECONDITION(pos, "bad vector");
    assert!(!pos.is_empty(), "bad vector");
    // RDKit✔️✔️: PRECONDITION(grad, "bad vector");
    assert!(!grad.is_empty(), "bad vector");
    assert!(
        pos.len() >= 3 * owner_points,
        "bad vector length for force-field positions"
    );
    assert!(
        grad.len() >= 3 * owner_points,
        "bad gradient length for force-field positions"
    );

    // RDKit✔️✔️: for (const auto &contrib : d_contribs) {
    for contrib in terms {
        // RDKit✔️✔️:   const RDGeom::Point3D iPoint(pos[3 * contrib.idx1],
        let i_point = point3(
            pos[3 * contrib.idx1],
            pos[3 * contrib.idx1 + 1],
            pos[3 * contrib.idx1 + 2],
        );
        // RDKit✔️✔️:   const RDGeom::Point3D jPoint(pos[3 * contrib.idx2],
        let j_point = point3(
            pos[3 * contrib.idx2],
            pos[3 * contrib.idx2 + 1],
            pos[3 * contrib.idx2 + 2],
        );
        // RDKit✔️✔️:   const RDGeom::Point3D kPoint(pos[3 * contrib.idx3],
        let k_point = point3(
            pos[3 * contrib.idx3],
            pos[3 * contrib.idx3 + 1],
            pos[3 * contrib.idx3 + 2],
        );
        // RDKit✔️✔️:   const RDGeom::Point3D lPoint(pos[3 * contrib.idx4],
        let l_point = point3(
            pos[3 * contrib.idx4],
            pos[3 * contrib.idx4 + 1],
            pos[3 * contrib.idx4 + 2],
        );

        // RDKit✔️✔️:   RDGeom::Point3D r[4] = {iPoint - jPoint, kPoint - jPoint, jPoint - kPoint,
        // RDKit✔️✔️:                           lPoint - kPoint};
        let r = [
            Point3::difference(&i_point, &j_point),
            Point3::difference(&k_point, &j_point),
            Point3::difference(&j_point, &k_point),
            Point3::difference(&l_point, &k_point),
        ];
        // RDKit✔️✔️:   RDGeom::Point3D t[2] = {r[0].crossProduct(r[1]), r[2].crossProduct(r[3])};
        let mut t = [r[0].cross_product(&r[1]), r[2].cross_product(&r[3])];
        // RDKit✔️✔️:   double d[2] = {t[0].length(), t[1].length()};
        let d = [t[0].length(), t[1].length()];
        // RDKit✔️✔️:   if (MMFF::isDoubleZero(d[0]) || MMFF::isDoubleZero(d[1])) {
        // RDKit✔️✔️:     return;
        // RDKit✔️✔️:   }
        if is_double_zero(d[0]) || is_double_zero(d[1]) {
            return Ok(());
        }
        // RDKit✔️✔️:   t[0].divide_assign(d[0]);
        // RDKit✔️✔️:   t[1].divide_assign(d[1]);
        t[0].divide_assign(d[0]);
        t[1].divide_assign(d[1]);
        // RDKit✔️✔️:   double cosPhi = t[0].dotProduct(t[1]);
        let mut cos_phi = t[0].dot_product(&t[1]);
        // RDKit✔️✔️:   cosPhi = std::clamp(cosPhi, -1.0, 1.0);
        cos_phi = cos_phi.clamp(-1.0, 1.0);
        // RDKit✔️✔️:   const double sinPhiSq = 1.0 - cosPhi * cosPhi;
        let sin_phi_sq = 1.0 - cos_phi * cos_phi;
        // RDKit✔️✔️:   const double sinPhi = ((sinPhiSq > 0.0) ? sqrt(sinPhiSq) : 0.0);
        let sin_phi = if sin_phi_sq > 0.0 {
            sin_phi_sq.sqrt()
        } else {
            0.0
        };
        // RDKit✔️✔️:   const double cosPhi2 = cosPhi * cosPhi;
        let cos_phi2 = cos_phi * cos_phi;
        // RDKit✔️✔️:   const double cosPhi3 = cosPhi * cosPhi2;
        let cos_phi3 = cos_phi * cos_phi2;
        // RDKit✔️✔️:   const double cosPhi4 = cosPhi * cosPhi3;
        let cos_phi4 = cos_phi * cos_phi3;
        // RDKit✔️✔️:   const double cosPhi5 = cosPhi * cosPhi4;
        let cos_phi5 = cos_phi * cos_phi4;
        // RDKit✔️✔️:   const double dE_dPhi =
        let d_e_d_phi = -contrib.force_constants[0] * f64::from(contrib.signs[0]) * sin_phi
                - 2.0
                    * contrib.force_constants[1]
                    * f64::from(contrib.signs[1])
                    * (2.0 * cos_phi * sin_phi)
                - 3.0
                    * contrib.force_constants[2]
                    * f64::from(contrib.signs[2])
                    * (4.0 * cos_phi2 * sin_phi - sin_phi)
                - 4.0
                    * contrib.force_constants[3]
                    * f64::from(contrib.signs[3])
                    * (8.0 * cos_phi3 * sin_phi - 4.0 * cos_phi * sin_phi)
                - 5.0
                    * contrib.force_constants[4]
                    * f64::from(contrib.signs[4])
                    * (16.0 * cos_phi4 * sin_phi - 12.0 * cos_phi2 * sin_phi + sin_phi)
                // RDKit✔️✔️:   const double dE_dPhi =
                // RDKit✔️✔️:       (-contrib.forceConstants[0] * contrib.signs[0] * sinPhi -
                // RDKit✔️✔️:        2.0 * contrib.forceConstants[1] * contrib.signs[1] *
                // RDKit✔️✔️:            (2.0 * cosPhi * sinPhi) -
                // RDKit✔️✔️:        3.0 * contrib.forceConstants[2] * contrib.signs[2] *
                // RDKit✔️✔️:            (4.0 * cosPhi2 * sinPhi - sinPhi) -
                // RDKit✔️✔️:        4.0 * contrib.forceConstants[3] * contrib.signs[3] *
                // RDKit✔️✔️:            (8.0 * cosPhi3 * sinPhi - 4.0 * cosPhi * sinPhi) -
                // RDKit✔️✔️:        5.0 * contrib.forceConstants[4] * contrib.signs[4] *
                // RDKit✔️✔️:            (16.0 * cosPhi4 * sinPhi - 12.0 * cosPhi2 * sinPhi + sinPhi) -
                // RDKit✔️✔️:        6.0 * contrib.forceConstants[4] * contrib.signs[4] *
                // RDKit✔️✔️:            (32.0 * cosPhi5 * sinPhi - 32.0 * cosPhi3 * sinPhi +
                // RDKit✔️✔️:             6.0 * sinPhi));
                //
                // RDKit source uses index 4 again in the sixth term. Exact parity requires
                // reproducing that source behavior, not the mathematically simplified variant.
                - 6.0
                    * contrib.force_constants[4]
                    * f64::from(contrib.signs[4])
                    * (32.0 * cos_phi5 * sin_phi - 32.0 * cos_phi3 * sin_phi + 6.0 * sin_phi);

        // RDKit✔️✔️:   double sinTerm = -dE_dPhi * (MMFF::isDoubleZero(sinPhi) ? (1.0 / cosPhi)
        // RDKit✔️✔️:                                                           : (1.0 / sinPhi));
        let sin_term = -d_e_d_phi
            * if is_double_zero(sin_phi) {
                1.0 / cos_phi
            } else {
                1.0 / sin_phi
            };

        // RDKit✔️✔️:   MMFF::Utils::calcTorsionGrad(r, t, d, g, sinTerm, cosPhi);
        calc_torsion_grad(
            &r,
            &t,
            &d,
            grad.as_chunks_mut::<3>().0,
            [contrib.idx1, contrib.idx2, contrib.idx3, contrib.idx4],
            sin_term,
            cos_phi,
        )
        .map_err(|_| ForceFieldKernelError::BadIndex)?;
    }
    // RDKit✔️✔️: }
    Ok(())
}
fn validate_indices(positions_len: usize, idx1: usize, idx2: usize, idx3: usize, idx4: usize) {
    assert!(
        idx1 != idx2
            && idx1 != idx3
            && idx1 != idx4
            && idx2 != idx3
            && idx2 != idx4
            && idx3 != idx4,
        "degenerate points"
    );
    // RDKit✔️✔️: URANGE_CHECK(idx1, dp_forceField->positions().size());
    // RDKit✔️✔️: URANGE_CHECK(idx2, dp_forceField->positions().size());
    // RDKit✔️✔️: URANGE_CHECK(idx3, dp_forceField->positions().size());
    // RDKit✔️✔️: URANGE_CHECK(idx4, dp_forceField->positions().size());

    assert!(idx1 < positions_len);
    assert!(idx2 < positions_len);
    assert!(idx3 < positions_len);
    assert!(idx4 < positions_len);
}

#[derive(Clone, Debug)]
pub(crate) struct TorsionAngleContribM6 {
    owner_points: Option<usize>,
    term: TorsionAngleContribsParams,
}
impl TorsionAngleContribM6 {
    pub(crate) fn new(
        owner: &ForceField<'_>,
        idx1: usize,
        idx2: usize,
        idx3: usize,
        idx4: usize,
        force_constants: Vec<f64>,
        signs: Vec<i32>,
    ) -> Self {
        // BEGIN COMPLETE PINNED CPP TorsionAngleContribM6::TorsionAngleContribM6 (TorsionAngleM6.cpp:46-65)
        // RDKit❗✔️: TorsionAngleContribM6::TorsionAngleContribM6(
        // RDKit❗✔️:     ForceFields::ForceField *owner, unsigned int idx1, unsigned int idx2,
        // RDKit❗✔️:     unsigned int idx3, unsigned int idx4, std::vector<double> V,
        // RDKit❗✔️:     std::vector<int> signs)
        // RDKit❗✔️:     : ForceFieldContrib(owner),
        // RDKit❗✔️:       d_at1Idx(idx1),
        // RDKit❗✔️:       d_at2Idx(idx2),
        // RDKit❗✔️:       d_at3Idx(idx3),
        // RDKit❗✔️:       d_at4Idx(idx4),
        // RDKit❗✔️:       d_V(std::move(V)),
        // RDKit❗✔️:       d_sign(std::move(signs)) {
        // RDKit❗✔️:   PRECONDITION(owner, "bad owner");
        // RDKit❗✔️:   PRECONDITION((idx1 != idx2) && (idx1 != idx3) && (idx1 != idx4) &&
        // RDKit❗✔️:                    (idx2 != idx3) && (idx2 != idx4) && (idx3 != idx4),
        // RDKit❗✔️:                "degenerate points");
        // RDKit❗✔️:   URANGE_CHECK(idx1, owner->positions().size());
        // RDKit❗✔️:   URANGE_CHECK(idx2, owner->positions().size());
        // RDKit❗✔️:   URANGE_CHECK(idx3, owner->positions().size());
        // RDKit❗✔️:   URANGE_CHECK(idx4, owner->positions().size());
        // RDKit❗✔️: };
        // END COMPLETE PINNED CPP TorsionAngleContribM6::TorsionAngleContribM6

        // RDKit❗✔️: Behavior pending independent paired review. The exact
        // source checks are shared with the collection. Complexity is O(1)
        // with moved parameter vectors and an inline term, no extra allocation.
        validate_indices(owner.positions().len(), idx1, idx2, idx3, idx4);
        Self {
            owner_points: Some(owner.positions().len()),
            term: TorsionAngleContribsParams::new(idx1, idx2, idx3, idx4, force_constants, signs),
        }
    }
    fn atom_indices(&self) -> [usize; 4] {
        let p = &self.term;
        [p.idx1, p.idx2, p.idx3, p.idx4]
    }
    fn force_constants(&self) -> &[f64] {
        &self.term.force_constants
    }
    fn signs(&self) -> &[i32] {
        &self.term.signs
    }
}
impl ForceFieldContribution for TorsionAngleContribM6 {
    fn get_energy(
        &self,
        context: &mut EvaluationContext<'_>,
    ) -> Result<f64, ForceFieldKernelError> {
        // BEGIN COMPLETE PINNED CPP TorsionAngleContribM6::getEnergy (TorsionAngleM6.cpp:67-82)
        // RDKit❗✔️: double TorsionAngleContribM6::getEnergy(double *pos) const {
        // RDKit❗✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit❗✔️:   PRECONDITION(pos, "bad vector");
        // RDKit❗✔️:
        // RDKit❗✔️:   RDGeom::Point3D iPoint(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
        // RDKit❗✔️:                          pos[3 * d_at1Idx + 2]);
        // RDKit❗✔️:   RDGeom::Point3D jPoint(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
        // RDKit❗✔️:                          pos[3 * d_at2Idx + 2]);
        // RDKit❗✔️:   RDGeom::Point3D kPoint(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
        // RDKit❗✔️:                          pos[3 * d_at3Idx + 2]);
        // RDKit❗✔️:   RDGeom::Point3D lPoint(pos[3 * d_at4Idx], pos[3 * d_at4Idx + 1],
        // RDKit❗✔️:                          pos[3 * d_at4Idx + 2]);
        // RDKit❗✔️:
        // RDKit❗✔️:   return calcTorsionEnergyM6(
        // RDKit❗✔️:       d_V, d_sign, Utils::calcTorsionCosPhi(iPoint, jPoint, kPoint, lPoint));
        // RDKit❗✔️: }
        // END COMPLETE PINNED CPP TorsionAngleContribM6::getEnergy

        energy_terms(self.owner_points, std::slice::from_ref(&self.term), context)
    }
    fn get_grad(
        &self,
        context: &mut EvaluationContext<'_>,
        gradient: &mut [f64],
    ) -> Result<(), ForceFieldKernelError> {
        // BEGIN COMPLETE PINNED CPP TorsionAngleContribM6::getGrad (TorsionAngleM6.cpp:84-136)
        // RDKit❗✔️: void TorsionAngleContribM6::getGrad(double *pos, double *grad) const {
        // RDKit❗✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit❗✔️:   PRECONDITION(pos, "bad vector");
        // RDKit❗✔️:   PRECONDITION(grad, "bad vector");
        // RDKit❗✔️:
        // RDKit❗✔️:   RDGeom::Point3D iPoint(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
        // RDKit❗✔️:                          pos[3 * d_at1Idx + 2]);
        // RDKit❗✔️:   RDGeom::Point3D jPoint(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
        // RDKit❗✔️:                          pos[3 * d_at2Idx + 2]);
        // RDKit❗✔️:   RDGeom::Point3D kPoint(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
        // RDKit❗✔️:                          pos[3 * d_at3Idx + 2]);
        // RDKit❗✔️:   RDGeom::Point3D lPoint(pos[3 * d_at4Idx], pos[3 * d_at4Idx + 1],
        // RDKit❗✔️:                          pos[3 * d_at4Idx + 2]);
        // RDKit❗✔️:   double *g[4] = {&(grad[3 * d_at1Idx]), &(grad[3 * d_at2Idx]),
        // RDKit❗✔️:                   &(grad[3 * d_at3Idx]), &(grad[3 * d_at4Idx])};
        // RDKit❗✔️:
        // RDKit❗✔️:   RDGeom::Point3D r[4] = {iPoint - jPoint, kPoint - jPoint, jPoint - kPoint,
        // RDKit❗✔️:                           lPoint - kPoint};
        // RDKit❗✔️:   RDGeom::Point3D t[2] = {r[0].crossProduct(r[1]), r[2].crossProduct(r[3])};
        // RDKit❗✔️:   double d[2] = {t[0].length(), t[1].length()};
        // RDKit❗✔️:   if (isDoubleZero(d[0]) || isDoubleZero(d[1])) {
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   t[0] /= d[0];
        // RDKit❗✔️:   t[1] /= d[1];
        // RDKit❗✔️:   double cosPhi = t[0].dotProduct(t[1]);
        // RDKit❗✔️:   clipToOne(cosPhi);
        // RDKit❗✔️:   double sinPhiSq = 1.0 - cosPhi * cosPhi;
        // RDKit❗✔️:   double sinPhi = ((sinPhiSq > 0.0) ? sqrt(sinPhiSq) : 0.0);
        // RDKit❗✔️:   double cosPhi2 = cosPhi * cosPhi;
        // RDKit❗✔️:   double cosPhi3 = cosPhi * cosPhi2;
        // RDKit❗✔️:   double cosPhi4 = cosPhi * cosPhi3;
        // RDKit❗✔️:   double cosPhi5 = cosPhi * cosPhi4;
        // RDKit❗✔️:   // dE/dPhi is independent of cartesians:
        // RDKit❗✔️:   double dE_dPhi =
        // RDKit❗✔️:       (-d_V[0] * d_sign[0] * sinPhi -
        // RDKit❗✔️:        2.0 * d_V[1] * d_sign[1] * (2.0 * cosPhi * sinPhi) -
        // RDKit❗✔️:        3.0 * d_V[2] * d_sign[2] * (4.0 * cosPhi2 * sinPhi - sinPhi) -
        // RDKit❗✔️:        4.0 * d_V[3] * d_sign[3] *
        // RDKit❗✔️:            (8.0 * cosPhi3 * sinPhi - 4.0 * cosPhi * sinPhi) -
        // RDKit❗✔️:        5.0 * d_V[4] * d_sign[4] *
        // RDKit❗✔️:            (16.0 * cosPhi4 * sinPhi - 12.0 * cosPhi2 * sinPhi + sinPhi) -
        // RDKit❗✔️:        6.0 * d_V[4] * d_sign[4] *
        // RDKit❗✔️:            (32.0 * cosPhi5 * sinPhi - 32.0 * cosPhi3 * sinPhi + 6.0 * sinPhi));
        // RDKit❗✔️:
        // RDKit❗✔️:   // FIX: use a tolerance here
        // RDKit❗✔️:   // this is hacky, but it's per the
        // RDKit❗✔️:   // recommendation from Niketic and Rasmussen:
        // RDKit❗✔️:   double sinTerm =
        // RDKit❗✔️:       -dE_dPhi * (isDoubleZero(sinPhi) ? (1.0 / cosPhi) : (1.0 / sinPhi));
        // RDKit❗✔️:
        // RDKit❗✔️:   Utils::calcTorsionGrad(r, t, d, g, sinTerm, cosPhi);
        // RDKit❗✔️: }
        // END COMPLETE PINNED CPP TorsionAngleContribM6::getGrad

        gradient_terms(
            self.owner_points,
            std::slice::from_ref(&self.term),
            context,
            gradient,
        )
    }
    fn copy(&self) -> Box<dyn ForceFieldContribution> {
        Box::new(self.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::super::torsion_preferences::CrystalFFDetails;
    use super::*;

    fn fixture_smiles(
        text: &str,
    ) -> Result<cosmolkit_model::TopologyBlock, cosmolkit_smiles::SmilesParseError> {
        let record = cosmolkit_smiles::parse_smiles(text, &Default::default())?;
        Ok(
            cosmolkit_core::sanitize_topology(&record.topology, &Default::default())
                .unwrap()
                .topology,
        )
    }
    fn get_experimental_torsions_without_bonds(
        mol: &cosmolkit_model::TopologyBlock,
        details: &mut CrystalFFDetails,
        exp: bool,
        small: bool,
        macrocycle: bool,
        basic: bool,
        version: u32,
        verbose: bool,
    ) -> Result<(), super::super::torsion_preferences::CrystalffTorsionPreferencesError> {
        let rings = cosmolkit_core::symmetrized_sssr(mol, &Default::default()).unwrap();
        let valence = cosmolkit_core::assign_valence_for_topology(
            mol,
            cosmolkit_core::ValenceModel::RdkitLike,
        )
        .unwrap();
        super::super::torsion_preferences::get_experimental_torsions_without_bonds(
            mol, &rings, &valence, details, exp, small, macrocycle, basic, version, verbose,
        )
    }
    fn torsion_rows() -> Vec<Vec<f64>> {
        vec![
            vec![0., 1.5, 0.],
            vec![0., 0., 0.],
            vec![1.5, 0., 0.],
            vec![1.5, 0., 1.5],
        ]
    }
    fn degenerate_rows() -> Vec<Vec<f64>> {
        vec![
            vec![0., 0., 0.],
            vec![1., 0., 0.],
            vec![2., 0., 0.],
            vec![3., 0., 0.],
        ]
    }
    fn rows_for_atom_count(atom_count: usize) -> Vec<Vec<f64>> {
        (0..atom_count)
            .map(|idx| {
                vec![
                    idx as f64,
                    if idx % 2 == 0 { 0. } else { 1. },
                    if idx % 3 == 0 { 0.5 } else { 0. },
                ]
            })
            .collect()
    }
    fn fixture_field<'a>(rows: &'a mut [Vec<f64>]) -> ForceField<'a> {
        let mut ff = ForceField::new(3);
        ff.positions_mut()
            .extend(rows.iter_mut().map(Vec::as_mut_slice));
        ff.initialize().unwrap();
        ff
    }
    fn flattened_positions(ff: &ForceField<'_>) -> Vec<f64> {
        ff.positions()
            .iter()
            .flat_map(|r| r.iter().copied())
            .collect()
    }
    fn evaluate_energy(contrib: &dyn ForceFieldContribution, pos: &[f64]) -> f64 {
        let n = 4;
        let mut cache = vec![-1.; n * (n + 1) / 2];
        let mut context = EvaluationContext::for_test(pos, &mut cache, n as u32);
        contrib.get_energy(&mut context).unwrap()
    }
    fn evaluate_grad(contrib: &dyn ForceFieldContribution, pos: &[f64], grad: &mut [f64]) {
        let n = 4;
        let mut cache = vec![-1.; n * (n + 1) / 2];
        let mut context = EvaluationContext::for_test(pos, &mut cache, n as u32);
        contrib.get_grad(&mut context, grad).unwrap();
    }
    fn test_cos_phi(pos: &[f64], i: usize, j: usize, k: usize, l: usize) -> f64 {
        let mut cos = 0.;
        crate::geometry::compute_dihedral_from_flat(
            pos,
            i,
            j,
            k,
            l,
            None,
            Some(&mut cos),
            None,
            None,
            None,
        );
        cos
    }
    const EPS_TEST: f64 = 1.0e-10;

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < EPS_TEST,
            "actual={actual} expected={expected}"
        );
    }

    fn add_details_contribs(contribs: &mut TorsionAngleContribs, details: &CrystalFFDetails) {
        for (atoms, (signs, force_constants)) in details
            .exp_torsion_atoms
            .iter()
            .zip(details.exp_torsion_angles.iter())
        {
            contribs.add_contrib(
                usize::try_from(atoms[0]).expect("non-negative atom index"),
                usize::try_from(atoms[1]).expect("non-negative atom index"),
                usize::try_from(atoms[2]).expect("non-negative atom index"),
                usize::try_from(atoms[3]).expect("non-negative atom index"),
                force_constants.clone(),
                signs.clone(),
            );
        }
    }

    #[test]
    fn crystalff_torsionanglecontribs_calc_torsion_energy_matches_m6_closed_form() {
        let force_constants = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let signs = vec![1, -1, 1, -1, 1, -1];
        let cos_phi = 0.5;

        let energy = calc_torsion_energy_m6(&force_constants, &signs, cos_phi);

        let cos2_phi = 2.0 * cos_phi * cos_phi - 1.0;
        let cos3_phi = 4.0 * cos_phi * cos_phi * cos_phi - 3.0 * cos_phi;
        let cos4_phi = 8.0 * cos_phi.powi(4) - 8.0 * cos_phi.powi(2) + 1.0;
        let cos5_phi = 16.0 * cos_phi.powi(5) - 20.0 * cos_phi.powi(3) + 5.0 * cos_phi;
        let cos6_phi =
            32.0 * cos_phi.powi(6) - 48.0 * cos_phi.powi(4) + 18.0 * cos_phi.powi(2) - 1.0;
        let expected = force_constants[0] * (1.0 + cos_phi)
            + force_constants[1] * (1.0 - cos2_phi)
            + force_constants[2] * (1.0 + cos3_phi)
            + force_constants[3] * (1.0 - cos4_phi)
            + force_constants[4] * (1.0 + cos5_phi)
            + force_constants[5] * (1.0 - cos6_phi);

        assert_close(energy, expected);
    }

    #[test]
    fn crystalff_torsion_energy_m6_matches_rdkit_closed_form_for_boundary_and_mixed_signs() {
        let force_constants = [1.25, -2.5, 0.75, 4.0, -5.5, 6.25];
        let signs = [1, -1, -1, 1, 1, -1];

        for cos_phi in [-1.0_f64, -0.25, 0.0, 0.5, 1.0] {
            let cos_phi2 = cos_phi * cos_phi;
            let cos_phi3 = cos_phi * cos_phi2;
            let cos_phi4 = cos_phi * cos_phi3;
            let cos_phi5 = cos_phi * cos_phi4;
            let cos_phi6 = cos_phi * cos_phi5;
            let cos2_phi = 2.0 * cos_phi2 - 1.0;
            let cos3_phi = 4.0 * cos_phi3 - 3.0 * cos_phi;
            let cos4_phi = 8.0 * cos_phi4 - 8.0 * cos_phi2 + 1.0;
            let cos5_phi = 16.0 * cos_phi5 - 20.0 * cos_phi3 + 5.0 * cos_phi;
            let cos6_phi = 32.0 * cos_phi6 - 48.0 * cos_phi4 + 18.0 * cos_phi2 - 1.0;
            let expected = force_constants[0] * (1.0 + f64::from(signs[0]) * cos_phi)
                + force_constants[1] * (1.0 + f64::from(signs[1]) * cos2_phi)
                + force_constants[2] * (1.0 + f64::from(signs[2]) * cos3_phi)
                + force_constants[3] * (1.0 + f64::from(signs[3]) * cos4_phi)
                + force_constants[4] * (1.0 + f64::from(signs[4]) * cos5_phi)
                + force_constants[5] * (1.0 + f64::from(signs[5]) * cos6_phi);

            assert_close(
                calc_torsion_energy_m6(&force_constants, &signs, cos_phi),
                expected,
            );
        }
    }

    #[test]
    fn crystalff_torsion_energy_m6_legacy_wrapper_matches_source_named_function() {
        let force_constants = [1.0, 0.5, 3.0, -2.0, 4.25, 0.125];
        let signs = [-1, -1, 1, 1, -1, 1];
        let cos_phi = 0.375;

        assert_eq!(
            calc_torsion_energy(&force_constants, &signs, cos_phi),
            calc_torsion_energy_m6(&force_constants, &signs, cos_phi)
        );
    }

    #[test]
    #[should_panic]
    fn crystalff_torsion_energy_m6_panics_for_short_force_constants_like_rdkit_indexing() {
        let _ = calc_torsion_energy_m6(&[1.0; 5], &[1; 6], 0.25);
    }

    #[test]
    #[should_panic]
    fn crystalff_torsion_energy_m6_panics_for_short_signs_like_rdkit_indexing() {
        let _ = calc_torsion_energy_m6(&[1.0; 6], &[1; 5], 0.25);
    }

    #[test]
    fn crystalff_torsionanglecontribs_constructor_starts_empty_and_supports_m6_gradient_path() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let pos = flattened_positions(&ff);
        let mut contribs = TorsionAngleContribs::new(&ff);

        assert!(contribs.is_empty());
        assert_eq!(contribs.len(), 0);

        let force_constants = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let signs = vec![1, -1, 1, -1, 1, -1];
        contribs.add_contrib(0, 1, 2, 3, force_constants.clone(), signs.clone());

        assert_eq!(contribs.len(), 1);

        let energy = evaluate_energy(&contribs, &pos);
        let expected_energy =
            calc_torsion_energy(&force_constants, &signs, test_cos_phi(&pos, 0, 1, 2, 3));
        let mut grad = vec![0.0; pos.len()];
        evaluate_grad(&contribs, &pos, &mut grad);

        assert_close(energy, expected_energy);
        assert!(grad.iter().all(|value| value.is_finite()));
        assert!(grad.iter().any(|value| value.abs() > EPS_TEST));
    }

    #[test]
    fn crystalff_torsionanglecontribs_constructor_accepts_smarts_small_ring_and_macrocycle_terms() {
        let linear = fixture_smiles("CCCCC").expect("pentane");
        let small_ring = fixture_smiles("C1COCC1").expect("tetrahydrofuran-like ring");
        let macrocycle = fixture_smiles("C1COCCCCCCC1").expect("macrocycle");
        let mut linear_details = CrystalFFDetails::default();
        let mut small_ring_details = CrystalFFDetails::default();
        let mut macrocycle_details = CrystalFFDetails::default();

        get_experimental_torsions_without_bonds(
            &linear,
            &mut linear_details,
            true,
            false,
            false,
            false,
            2,
            false,
        )
        .expect("linear SMARTS torsions");
        get_experimental_torsions_without_bonds(
            &small_ring,
            &mut small_ring_details,
            true,
            true,
            false,
            false,
            1,
            false,
        )
        .expect("small-ring torsions");
        get_experimental_torsions_without_bonds(
            &macrocycle,
            &mut macrocycle_details,
            true,
            false,
            true,
            false,
            1,
            false,
        )
        .expect("macrocycle torsions");

        assert_eq!(linear_details.exp_torsion_atoms.len(), 2);
        assert!(!small_ring_details.exp_torsion_atoms.is_empty());
        assert!(!macrocycle_details.exp_torsion_atoms.is_empty());

        let mut linear_rows = rows_for_atom_count(linear.atoms.len());
        let linear_ff = fixture_field(&mut linear_rows);
        let mut linear_contribs = TorsionAngleContribs::new(&linear_ff);
        add_details_contribs(&mut linear_contribs, &linear_details);
        assert_eq!(
            linear_contribs.len(),
            linear_details.exp_torsion_atoms.len()
        );

        let mut small_ring_rows = rows_for_atom_count(small_ring.atoms.len());
        let small_ring_ff = fixture_field(&mut small_ring_rows);
        let mut small_ring_contribs = TorsionAngleContribs::new(&small_ring_ff);
        add_details_contribs(&mut small_ring_contribs, &small_ring_details);
        assert_eq!(
            small_ring_contribs.len(),
            small_ring_details.exp_torsion_atoms.len()
        );

        let mut macrocycle_rows = rows_for_atom_count(macrocycle.atoms.len());
        let macrocycle_ff = fixture_field(&mut macrocycle_rows);
        let mut macrocycle_contribs = TorsionAngleContribs::new(&macrocycle_ff);
        add_details_contribs(&mut macrocycle_contribs, &macrocycle_details);
        assert_eq!(
            macrocycle_contribs.len(),
            macrocycle_details.exp_torsion_atoms.len()
        );
    }

    #[test]
    fn crystalff_torsionanglecontribs_addcontrib_accepts_smarts_small_ring_and_macrocycle_terms() {
        let linear = fixture_smiles("CCCCC").expect("pentane");
        let small_ring = fixture_smiles("C1COCC1").expect("tetrahydrofuran-like ring");
        let macrocycle = fixture_smiles("C1COCCCCCCC1").expect("macrocycle");
        let mut linear_details = CrystalFFDetails::default();
        let mut small_ring_details = CrystalFFDetails::default();
        let mut macrocycle_details = CrystalFFDetails::default();

        get_experimental_torsions_without_bonds(
            &linear,
            &mut linear_details,
            true,
            false,
            false,
            false,
            2,
            false,
        )
        .expect("linear SMARTS torsions");
        get_experimental_torsions_without_bonds(
            &small_ring,
            &mut small_ring_details,
            true,
            true,
            false,
            false,
            1,
            false,
        )
        .expect("small-ring torsions");
        get_experimental_torsions_without_bonds(
            &macrocycle,
            &mut macrocycle_details,
            true,
            false,
            true,
            false,
            1,
            false,
        )
        .expect("macrocycle torsions");

        let mut linear_rows = rows_for_atom_count(linear.atoms.len());
        let linear_ff = fixture_field(&mut linear_rows);
        let mut linear_contribs = TorsionAngleContribs::new(&linear_ff);
        add_details_contribs(&mut linear_contribs, &linear_details);
        assert_eq!(
            linear_contribs.len(),
            linear_details.exp_torsion_atoms.len()
        );

        let mut small_ring_rows = rows_for_atom_count(small_ring.atoms.len());
        let small_ring_ff = fixture_field(&mut small_ring_rows);
        let mut small_ring_contribs = TorsionAngleContribs::new(&small_ring_ff);
        add_details_contribs(&mut small_ring_contribs, &small_ring_details);
        assert_eq!(
            small_ring_contribs.len(),
            small_ring_details.exp_torsion_atoms.len()
        );

        let mut macrocycle_rows = rows_for_atom_count(macrocycle.atoms.len());
        let macrocycle_ff = fixture_field(&mut macrocycle_rows);
        let mut macrocycle_contribs = TorsionAngleContribs::new(&macrocycle_ff);
        add_details_contribs(&mut macrocycle_contribs, &macrocycle_details);
        assert_eq!(
            macrocycle_contribs.len(),
            macrocycle_details.exp_torsion_atoms.len()
        );
    }

    #[test]
    fn crystalff_torsionanglecontribs_addcontrib_supports_m6_energy_and_gradient_paths() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let pos = flattened_positions(&ff);
        let mut contribs = TorsionAngleContribs::new(&ff);
        let force_constants = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let signs = vec![1, -1, 1, -1, 1, -1];

        contribs.add_contrib(0, 1, 2, 3, force_constants.clone(), signs.clone());

        assert_eq!(contribs.len(), 1);
        let energy = evaluate_energy(&contribs, &pos);
        let mut grad = vec![0.0; pos.len()];
        evaluate_grad(&contribs, &pos, &mut grad);

        assert_close(
            energy,
            calc_torsion_energy(&force_constants, &signs, test_cos_phi(&pos, 0, 1, 2, 3)),
        );
        assert!(grad.iter().all(|value| value.is_finite()));
        assert!(grad.iter().any(|value| value.abs() > EPS_TEST));
    }

    #[test]
    #[should_panic(expected = "degenerate points")]
    fn crystalff_torsionanglecontribs_addcontrib_panics_for_degenerate_indices() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let mut contribs = TorsionAngleContribs::new(&ff);

        contribs.add_contrib(0, 0, 2, 3, vec![0.0; 6], vec![1; 6]);
    }

    #[test]
    #[should_panic(expected = "assertion failed: idx1 < positions_len")]
    fn crystalff_torsionanglecontribs_addcontrib_panics_for_out_of_range_first_index() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let mut contribs = TorsionAngleContribs::new(&ff);

        contribs.add_contrib(4, 1, 2, 3, vec![0.0; 6], vec![1; 6]);
    }

    #[test]
    #[should_panic(expected = "assertion failed: idx2 < positions_len")]
    fn crystalff_torsionanglecontribs_addcontrib_panics_for_out_of_range_second_index() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let mut contribs = TorsionAngleContribs::new(&ff);

        contribs.add_contrib(0, 4, 2, 3, vec![0.0; 6], vec![1; 6]);
    }

    #[test]
    #[should_panic(expected = "assertion failed: idx3 < positions_len")]
    fn crystalff_torsionanglecontribs_addcontrib_panics_for_out_of_range_third_index() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let mut contribs = TorsionAngleContribs::new(&ff);

        contribs.add_contrib(0, 1, 4, 3, vec![0.0; 6], vec![1; 6]);
    }

    #[test]
    #[should_panic(expected = "assertion failed: idx4 < positions_len")]
    fn crystalff_torsionanglecontribs_addcontrib_panics_for_out_of_range_fourth_index() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let mut contribs = TorsionAngleContribs::new(&ff);

        contribs.add_contrib(0, 1, 2, 4, vec![0.0; 6], vec![1; 6]);
    }

    #[test]
    fn crystalff_torsionanglecontribs_get_energy_matches_sp3_sp3_source_geometry() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let pos = flattened_positions(&ff);
        let mut contrib = TorsionAngleContribs::new(&ff);
        let signs = vec![1, 1, 1, 1, 1, 1];
        let mut force_constants = vec![0.0; 6];
        force_constants[2] = 4.0;
        contrib.add_contrib(0, 1, 2, 3, force_constants, signs);

        let energy = evaluate_energy(&contrib, &pos);
        let cos_phi = test_cos_phi(&pos, 0, 1, 2, 3);

        assert_close(cos_phi, 0.0);
        assert_close(energy, 4.0);
    }

    #[test]
    fn crystalff_torsionanglecontribs_get_energy_accumulates_multiple_contributions() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let pos = flattened_positions(&ff);
        let mut contrib = TorsionAngleContribs::new(&ff);
        let mut first = vec![0.0; 6];
        first[2] = 4.0;
        let mut second = vec![0.0; 6];
        second[0] = 1.5;
        second[1] = 0.5;
        contrib.add_contrib(0, 1, 2, 3, first.clone(), vec![1; 6]);
        contrib.add_contrib(0, 1, 2, 3, second.clone(), vec![1, -1, 1, 1, 1, 1]);

        let energy = evaluate_energy(&contrib, &pos);
        let expected = calc_torsion_energy(&first, &[1, 1, 1, 1, 1, 1], 0.0)
            + calc_torsion_energy(&second, &[1, -1, 1, 1, 1, 1], 0.0);

        assert_close(energy, expected);
    }

    #[test]
    fn crystalff_torsionanglecontribs_get_energy_returns_zero_when_empty() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let pos = flattened_positions(&ff);
        let contrib = TorsionAngleContribs::new(&ff);

        let energy = evaluate_energy(&contrib, &pos);

        assert_close(energy, 0.0);
    }

    #[test]
    #[should_panic(expected = "no owner")]
    fn crystalff_torsionanglecontribs_get_energy_panics_without_owner() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let pos = flattened_positions(&ff);
        let mut contrib = TorsionAngleContribs::new(&ff);
        contrib.owner_points = None;

        let _ = evaluate_energy(&contrib, &pos);
    }

    #[test]
    #[should_panic(expected = "bad vector")]
    fn crystalff_torsionanglecontribs_get_energy_panics_for_empty_position_vector() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let contrib = TorsionAngleContribs::new(&ff);

        let _ = evaluate_energy(&contrib, &[]);
    }

    #[test]
    #[should_panic(expected = "bad vector length for force-field positions")]
    fn crystalff_torsionanglecontribs_get_energy_panics_for_short_position_vector() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let contrib = TorsionAngleContribs::new(&ff);

        let _ = evaluate_energy(&contrib, &vec![0.0; 11]);
    }

    #[test]
    fn crystalff_torsionanglecontribs_get_grad_accumulates_multiple_contributions() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let pos = flattened_positions(&ff);
        let mut contrib = TorsionAngleContribs::new(&ff);
        let mut first = vec![0.0; 6];
        first[2] = 4.0;
        let mut second = vec![0.0; 6];
        second[0] = 1.5;
        second[1] = 0.5;
        contrib.add_contrib(0, 1, 2, 3, first, vec![1; 6]);
        contrib.add_contrib(0, 1, 2, 3, second, vec![1, -1, 1, 1, 1, 1]);

        let mut grad = vec![0.0; pos.len()];
        evaluate_grad(&contrib, &pos, &mut grad);

        assert!(grad.iter().all(|value| value.is_finite()));
        assert!(
            grad.iter().any(|value| value.abs() > EPS_TEST),
            "expected non-zero accumulated gradient"
        );
    }

    #[test]
    fn crystalff_torsion_angle_contribs_get_grad_matches_single_m6_contrib_exactly() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let pos = flattened_positions(&ff);
        let force_constants = vec![1.25, 2.5, 3.75, 4.5, 5.25, 6.0];
        let signs = vec![1, -1, 1, -1, 1, -1];
        let mut contribs = TorsionAngleContribs::new(&ff);
        contribs.add_contrib(0, 1, 2, 3, force_constants.clone(), signs.clone());
        let single = TorsionAngleContribM6::new(&ff, 0, 1, 2, 3, force_constants, signs);
        let mut collection_grad = vec![0.0; pos.len()];
        let mut single_grad = vec![0.0; pos.len()];

        evaluate_grad(&contribs, &pos, &mut collection_grad);
        evaluate_grad(&single, &pos, &mut single_grad);

        assert_eq!(collection_grad, single_grad);
    }

    #[test]
    fn crystalff_torsionanglecontribs_get_grad_returns_early_for_degenerate_torsion() {
        let mut rows = degenerate_rows();
        let ff = fixture_field(&mut rows);
        let pos = flattened_positions(&ff);
        let mut contrib = TorsionAngleContribs::new(&ff);
        let mut force_constants = vec![0.0; 6];
        force_constants[2] = 4.0;
        contrib.add_contrib(0, 1, 2, 3, force_constants, vec![1; 6]);
        let mut grad = vec![3.0; pos.len()];

        evaluate_grad(&contrib, &pos, &mut grad);

        assert!(grad.iter().all(|value| *value == 3.0));
    }

    #[test]
    fn crystalff_torsionanglecontribs_get_grad_zero_force_constants_leave_gradient_unchanged() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let pos = flattened_positions(&ff);
        let mut contrib = TorsionAngleContribs::new(&ff);
        contrib.add_contrib(0, 1, 2, 3, vec![0.0; 6], vec![1; 6]);
        let mut grad = vec![0.0; pos.len()];

        evaluate_grad(&contrib, &pos, &mut grad);

        assert!(grad.iter().all(|value| value.abs() < EPS_TEST));
    }

    #[test]
    #[should_panic(expected = "no owner")]
    fn crystalff_torsionanglecontribs_get_grad_panics_without_owner() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let pos = flattened_positions(&ff);
        let mut contrib = TorsionAngleContribs::new(&ff);
        contrib.owner_points = None;
        let mut grad = vec![0.0; pos.len()];

        evaluate_grad(&contrib, &pos, &mut grad);
    }

    #[test]
    #[should_panic(expected = "bad vector")]
    fn crystalff_torsionanglecontribs_get_grad_panics_for_empty_position_vector() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let contrib = TorsionAngleContribs::new(&ff);
        let mut grad = vec![0.0; 12];

        evaluate_grad(&contrib, &[], &mut grad);
    }

    #[test]
    #[should_panic(expected = "bad vector")]
    fn crystalff_torsionanglecontribs_get_grad_panics_for_empty_gradient_vector() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let pos = flattened_positions(&ff);
        let contrib = TorsionAngleContribs::new(&ff);
        let mut grad = vec![];

        evaluate_grad(&contrib, &pos, &mut grad);
    }

    #[test]
    #[should_panic(expected = "bad vector length for force-field positions")]
    fn crystalff_torsionanglecontribs_get_grad_panics_for_short_position_vector() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let contrib = TorsionAngleContribs::new(&ff);
        let mut grad = vec![0.0; 12];

        evaluate_grad(&contrib, &vec![0.0; 11], &mut grad);
    }

    #[test]
    #[should_panic(expected = "bad gradient length for force-field positions")]
    fn crystalff_torsionanglecontribs_get_grad_panics_for_short_gradient_vector() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let pos = flattened_positions(&ff);
        let contrib = TorsionAngleContribs::new(&ff);
        let mut grad = vec![0.0; 11];

        evaluate_grad(&contrib, &pos, &mut grad);
    }

    #[test]
    fn crystalff_torsionanglecontribm6_constructor_initializes_fields() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let force_constants = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let signs = vec![1, -1, 1, -1, 1, -1];

        let contrib =
            TorsionAngleContribM6::new(&ff, 0, 1, 2, 3, force_constants.clone(), signs.clone());

        assert_eq!(contrib.owner_points, Some(ff.positions().len()));
        assert_eq!(contrib.atom_indices(), [0, 1, 2, 3]);
        assert_eq!(contrib.force_constants(), force_constants.as_slice());
        assert_eq!(contrib.signs(), signs.as_slice());
    }

    #[test]
    #[should_panic(expected = "degenerate points")]
    fn crystalff_torsionanglecontribm6_constructor_panics_for_degenerate_indices() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);

        let _ = TorsionAngleContribM6::new(&ff, 0, 0, 2, 3, vec![0.0; 6], vec![1; 6]);
    }

    #[test]
    #[should_panic]
    fn crystalff_torsionanglecontribm6_constructor_panics_for_out_of_range_first_index() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);

        let _ = TorsionAngleContribM6::new(&ff, 4, 1, 2, 3, vec![0.0; 6], vec![1; 6]);
    }

    #[test]
    #[should_panic]
    fn crystalff_torsionanglecontribm6_constructor_panics_for_out_of_range_second_index() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);

        let _ = TorsionAngleContribM6::new(&ff, 0, 4, 2, 3, vec![0.0; 6], vec![1; 6]);
    }

    #[test]
    #[should_panic]
    fn crystalff_torsionanglecontribm6_constructor_panics_for_out_of_range_third_index() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);

        let _ = TorsionAngleContribM6::new(&ff, 0, 1, 4, 3, vec![0.0; 6], vec![1; 6]);
    }

    #[test]
    #[should_panic]
    fn crystalff_torsionanglecontribm6_constructor_panics_for_out_of_range_fourth_index() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);

        let _ = TorsionAngleContribM6::new(&ff, 0, 1, 2, 4, vec![0.0; 6], vec![1; 6]);
    }

    #[test]
    fn crystalff_torsionanglecontribm6_get_energy_matches_sp3_sp3_source_geometry() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let pos = flattened_positions(&ff);
        let force_constants = vec![0.0, 0.0, 4.0, 0.0, 0.0, 0.0];
        let signs = vec![1, 1, 1, 1, 1, 1];
        let contrib = TorsionAngleContribM6::new(&ff, 0, 1, 2, 3, force_constants, signs);

        let energy = evaluate_energy(&contrib, &pos);
        let cos_phi = test_cos_phi(&pos, 0, 1, 2, 3);

        assert_eq!(cos_phi, 0.0);
        assert_eq!(energy, 4.0);
    }

    #[test]
    fn crystalff_torsionanglecontribm6_get_energy_matches_m6_closed_form_for_geometry() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let pos = flattened_positions(&ff);
        let force_constants = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let signs = vec![1, -1, 1, -1, 1, -1];
        let contrib =
            TorsionAngleContribM6::new(&ff, 0, 1, 2, 3, force_constants.clone(), signs.clone());

        let energy = evaluate_energy(&contrib, &pos);
        let cos_phi = test_cos_phi(&pos, 0, 1, 2, 3);
        let expected = calc_torsion_energy_m6(&force_constants, &signs, cos_phi);

        assert_eq!(energy, expected);
    }

    #[test]
    #[should_panic(expected = "no owner")]
    fn crystalff_torsionanglecontribm6_get_energy_panics_without_owner() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let pos = flattened_positions(&ff);
        let mut contrib = TorsionAngleContribM6::new(&ff, 0, 1, 2, 3, vec![0.0; 6], vec![1; 6]);
        contrib.owner_points = None;

        let _ = evaluate_energy(&contrib, &pos);
    }

    #[test]
    #[should_panic(expected = "bad vector")]
    fn crystalff_torsionanglecontribm6_get_energy_panics_for_empty_position_vector() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let contrib = TorsionAngleContribM6::new(&ff, 0, 1, 2, 3, vec![0.0; 6], vec![1; 6]);

        let _ = evaluate_energy(&contrib, &[]);
    }

    #[test]
    #[should_panic(expected = "bad vector length for force-field positions")]
    fn crystalff_torsionanglecontribm6_get_energy_panics_for_short_position_vector() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let contrib = TorsionAngleContribM6::new(&ff, 0, 1, 2, 3, vec![0.0; 6], vec![1; 6]);

        let _ = evaluate_energy(&contrib, &vec![0.0; 3 * ff.positions().len() - 1]);
    }

    #[test]
    fn crystalff_torsionanglecontribm6_get_grad_zero_force_constants_leave_gradient_unchanged() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let pos = flattened_positions(&ff);
        let mut grad = vec![1.25; 3 * ff.positions().len()];
        let before = grad.clone();
        let contrib = TorsionAngleContribM6::new(&ff, 0, 1, 2, 3, vec![0.0; 6], vec![1; 6]);

        evaluate_grad(&contrib, &pos, &mut grad);

        assert_eq!(grad, before);
    }

    #[test]
    fn crystalff_torsionanglecontribm6_get_grad_matches_single_contribs_gradient() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let pos = flattened_positions(&ff);
        let force_constants = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let signs = vec![1, -1, 1, -1, 1, -1];
        let contrib =
            TorsionAngleContribM6::new(&ff, 0, 1, 2, 3, force_constants.clone(), signs.clone());
        let mut grad = vec![0.0; pos.len()];
        let mut contribs = TorsionAngleContribs::new(&ff);
        contribs.add_contrib(0, 1, 2, 3, force_constants, signs);
        let mut expected_grad = vec![0.0; pos.len()];

        evaluate_grad(&contrib, &pos, &mut grad);
        evaluate_grad(&contribs, &pos, &mut expected_grad);

        for (actual, expected) in grad.iter().zip(expected_grad.iter()) {
            assert_close(*actual, *expected);
        }
    }

    #[test]
    fn crystalff_torsion_angle_m6_get_energy_and_grad_match_collection_contrib() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let pos = flattened_positions(&ff);
        let force_constants = vec![1.25, 2.5, 3.75, 4.5, 5.25, 6.0];
        let signs = vec![1, -1, 1, -1, 1, -1];
        let contrib =
            TorsionAngleContribM6::new(&ff, 0, 1, 2, 3, force_constants.clone(), signs.clone());
        let mut collection = TorsionAngleContribs::new(&ff);
        collection.add_contrib(0, 1, 2, 3, force_constants, signs);
        let mut contrib_grad = vec![0.0; pos.len()];
        let mut collection_grad = vec![0.0; pos.len()];

        let contrib_energy = evaluate_energy(&contrib, &pos);
        let collection_energy = evaluate_energy(&collection, &pos);
        evaluate_grad(&contrib, &pos, &mut contrib_grad);
        evaluate_grad(&collection, &pos, &mut collection_grad);

        assert_eq!(contrib_energy, collection_energy);
        assert_eq!(contrib_grad, collection_grad);
    }

    #[test]
    fn crystalff_torsionanglecontribm6_get_grad_returns_early_for_degenerate_torsion() {
        let mut rows = degenerate_rows();
        let ff = fixture_field(&mut rows);
        let pos = flattened_positions(&ff);
        let mut grad = vec![1.0; 3 * ff.positions().len()];
        let before = grad.clone();
        let contrib = TorsionAngleContribM6::new(
            &ff,
            0,
            1,
            2,
            3,
            vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            vec![1, -1, 1, -1, 1, -1],
        );

        evaluate_grad(&contrib, &pos, &mut grad);

        assert_eq!(grad, before);
    }

    #[test]
    #[should_panic(expected = "no owner")]
    fn crystalff_torsionanglecontribm6_get_grad_panics_without_owner() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let pos = flattened_positions(&ff);
        let mut grad = vec![0.0; 3 * ff.positions().len()];
        let mut contrib = TorsionAngleContribM6::new(&ff, 0, 1, 2, 3, vec![0.0; 6], vec![1; 6]);
        contrib.owner_points = None;

        evaluate_grad(&contrib, &pos, &mut grad);
    }

    #[test]
    #[should_panic(expected = "bad vector")]
    fn crystalff_torsionanglecontribm6_get_grad_panics_for_empty_position_vector() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let mut grad = vec![0.0; 3 * ff.positions().len()];
        let contrib = TorsionAngleContribM6::new(&ff, 0, 1, 2, 3, vec![0.0; 6], vec![1; 6]);

        evaluate_grad(&contrib, &[], &mut grad);
    }

    #[test]
    #[should_panic(expected = "bad vector")]
    fn crystalff_torsionanglecontribm6_get_grad_panics_for_empty_gradient_vector() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let pos = flattened_positions(&ff);
        let mut grad = Vec::new();
        let contrib = TorsionAngleContribM6::new(&ff, 0, 1, 2, 3, vec![0.0; 6], vec![1; 6]);

        evaluate_grad(&contrib, &pos, &mut grad);
    }

    #[test]
    #[should_panic(expected = "bad vector length for force-field positions")]
    fn crystalff_torsionanglecontribm6_get_grad_panics_for_short_position_vector() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let mut grad = vec![0.0; 3 * ff.positions().len()];
        let contrib = TorsionAngleContribM6::new(&ff, 0, 1, 2, 3, vec![0.0; 6], vec![1; 6]);

        evaluate_grad(
            &contrib,
            &vec![0.0; 3 * ff.positions().len() - 1],
            &mut grad,
        );
    }

    #[test]
    #[should_panic(expected = "bad gradient length for force-field positions")]
    fn crystalff_torsionanglecontribm6_get_grad_panics_for_short_gradient_vector() {
        let mut rows = torsion_rows();
        let ff = fixture_field(&mut rows);
        let pos = flattened_positions(&ff);
        let mut grad = vec![0.0; 3 * ff.positions().len() - 1];
        let contrib = TorsionAngleContribM6::new(&ff, 0, 1, 2, 3, vec![0.0; 6], vec![1; 6]);

        evaluate_grad(&contrib, &pos, &mut grad);
    }
}

/// Detached evaluations of the two source torsion contribution forms.
#[doc(hidden)]
#[derive(Debug, Clone, PartialEq)]
pub struct CrystalTorsionPairEvaluation {
    pub collection_energy: f64,
    pub single_energy: f64,
    pub collection_gradient: Vec<f64>,
    pub single_gradient: Vec<f64>,
}
#[derive(Debug)]
enum CrystalTorsionEvaluationFailure {
    Input(&'static str),
    Kernel(ForceFieldKernelError),
}
#[doc(hidden)]
#[derive(Debug)]
pub struct CrystalTorsionEvaluationError(CrystalTorsionEvaluationFailure);
impl std::fmt::Display for CrystalTorsionEvaluationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &self.0 {
            CrystalTorsionEvaluationFailure::Input(message) => f.write_str(message),
            CrystalTorsionEvaluationFailure::Kernel(error) => error.fmt(f),
        }
    }
}
impl std::error::Error for CrystalTorsionEvaluationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match &self.0 {
            CrystalTorsionEvaluationFailure::Input(_) => None,
            CrystalTorsionEvaluationFailure::Kernel(error) => Some(error),
        }
    }
}
impl From<ForceFieldKernelError> for CrystalTorsionEvaluationError {
    fn from(error: ForceFieldKernelError) -> Self {
        Self(CrystalTorsionEvaluationFailure::Kernel(error))
    }
}
/// Checked value evaluation for the conformer-required original diagnostic.
/// No field, positional storage, cache or mutable geometry leaves this owner.
#[doc(hidden)]
pub fn evaluate_crystal_torsion_pair(
    coordinates: &[[f64; 3]],
    indices: [usize; 4],
    force_constants: &[f64],
    signs: &[i32],
) -> Result<CrystalTorsionPairEvaluation, CrystalTorsionEvaluationError> {
    // Mechanical source diagnostic transport: the unique source constructors,
    // energy and full gradient anchors remain in this module's actual owners.
    // RDKit❗✔️: One source field initialization and two original contribution
    // forms, then collection energy/single energy/collection gradient/single
    // gradient in the original order. All chemistry delegates to the unique
    // kernels; typed checks only guard native slice/index representation.
    // Complexity: one O(N) coordinate/flat/position descriptor copy and the
    // same O(N^2) source distance cache, two O(N) zero gradient outputs and two
    // six-term parameter copies. The original diagnostic has these buffers;
    // this boundary never clones topology, molecule state or any conformers.
    let input_error =
        |message| CrystalTorsionEvaluationError(CrystalTorsionEvaluationFailure::Input(message));
    if u32::try_from(coordinates.len()).is_err() {
        return Err(input_error("position count exceeds source unsigned range"));
    }
    let mut positions = coordinates.to_vec();
    let mut field = ForceField::new(3);
    field
        .positions_mut()
        .extend(positions.iter_mut().map(|row| row.as_mut_slice()));
    field.initialize()?;
    let [i, j, k, l] = indices;
    if i == j || i == k || i == l || j == k || j == l || k == l {
        return Err(input_error("degenerate points"));
    }
    if indices.iter().any(|&index| index >= coordinates.len()) {
        return Err(input_error("torsion point index out of range"));
    }
    if force_constants.len() < 6 || signs.len() < 6 {
        return Err(input_error("bad torsion term vector"));
    }
    let mut collection = TorsionAngleContribs::new(&field);
    collection.add_contrib(i, j, k, l, force_constants.to_vec(), signs.to_vec());
    let single =
        TorsionAngleContribM6::new(&field, i, j, k, l, force_constants.to_vec(), signs.to_vec());
    let pos: Vec<f64> = coordinates
        .iter()
        .flat_map(|row| row.iter().copied())
        .collect();
    let mut collection_gradient = vec![0.; pos.len()];
    let mut single_gradient = vec![0.; pos.len()];
    let mut context = field.evaluation_context(&pos);
    let collection_energy = collection.get_energy(&mut context)?;
    let single_energy = single.get_energy(&mut context)?;
    collection.get_grad(&mut context, &mut collection_gradient)?;
    single.get_grad(&mut context, &mut single_gradient)?;
    Ok(CrystalTorsionPairEvaluation {
        collection_energy,
        single_energy,
        collection_gradient,
        single_gradient,
    })
}
