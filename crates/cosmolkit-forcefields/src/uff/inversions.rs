// Copyright (C) 2024 Niels Maeder and other RDKit contributors.
//
// This file contains behavior ported from RDKit and remains subject to its
// BSD license, included in the pinned source tree.

use crate::kernel::EvaluationContext;

use super::inversion::{InversionContrib, InversionContributionError};

#[derive(Clone, Debug, Default)]
pub(super) struct InversionContribs {
    // One ordered owner for each packed source term. The per-term type already
    // owns the indices and coefficients used by both energy and gradient.
    contribs: Vec<InversionContrib>,
}

impl InversionContribs {
    pub(super) fn add_contrib(
        &mut self,
        positions: &[&mut [f64]],
        idx1: u32,
        idx2: u32,
        idx3: u32,
        idx4: u32,
        at2_atomic_num: i32,
        is_c_bound_to_o: bool,
        oob_force_scaling_factor: f64,
    ) -> Result<(), InversionContributionError> {
        // BEGIN RDKIT CPP FUNCTION ForceFields::UFF::InversionContribs::addContrib (ForceField/UFF/Inversions.cpp:27-43)
        // RDKit✔️✔️: void InversionContribs::addContrib(unsigned int idx1, unsigned int idx2,
        // RDKit✔️✔️:                                    unsigned int idx3, unsigned int idx4,
        // RDKit✔️✔️:                                    int at2AtomicNum, bool isCBoundToO,
        // RDKit✔️✔️:                                    double oobForceScalingFactor) {
        // RDKit✔️✔️:   URANGE_CHECK(idx1, dp_forceField->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx2, dp_forceField->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx3, dp_forceField->positions().size());
        // RDKit✔️✔️:   URANGE_CHECK(idx4, dp_forceField->positions().size());
        // The borrowed position list supplies the source count and replaces
        // the stored owner pointer. Each check occurs before term creation.
        // RDKit✔️✔️:   auto invCoeffForceCon = Utils::calcInversionCoefficientsAndForceConstant(
        // RDKit✔️✔️:       at2AtomicNum, isCBoundToO);
        // RDKit✔️✔️:   d_contribs.emplace_back(
        // RDKit✔️✔️:       idx1, idx2, idx3, idx4, at2AtomicNum, isCBoundToO,
        // RDKit✔️✔️:       std::get<1>(invCoeffForceCon), std::get<2>(invCoeffForceCon),
        // RDKit✔️✔️:       std::get<3>(invCoeffForceCon),
        // RDKit✔️✔️:       std::get<0>(invCoeffForceCon) * oobForceScalingFactor);
        // RDKit✔️✔️: }
        let contribution = InversionContrib::new_packed_with_scale(
            positions,
            idx1,
            idx2,
            idx3,
            idx4,
            at2_atomic_num,
            is_c_bound_to_o,
            oob_force_scaling_factor,
        )?;
        self.contribs.push(contribution);
        // END RDKIT CPP FUNCTION ForceFields::UFF::InversionContribs::addContrib
        Ok(())
    }

    pub(super) fn get_energy(&self, context: &mut EvaluationContext<'_>) -> f64 {
        // BEGIN RDKIT CPP FUNCTION ForceFields::UFF::InversionContribs::getEnergy (ForceField/UFF/Inversions.cpp:45-64)
        // RDKit✔️✔️: double InversionContribs::getEnergy(double *pos) const {
        // RDKit✔️✔️:   PRECONDITION(dp_forceField, "no owner");
        // RDKit✔️✔️:   PRECONDITION(pos, "bad vector");
        // The borrowed EvaluationContext and coordinates replace these
        // nonnull raw pointers; no force-field owner is retained.
        // RDKit✔️✔️:   double accum = 0;
        let mut accum = 0.0;
        // RDKit✔️✔️:   for (const auto &contrib : d_contribs) {
        for contribution in &self.contribs {
            // RDKit✔️✔️:     const RDGeom::Point3D p1(pos[3 * contrib.idx1], pos[3 * contrib.idx1 + 1],
            // RDKit✔️✔️:                              pos[3 * contrib.idx1 + 2]);
            // RDKit✔️✔️:     const RDGeom::Point3D p2(pos[3 * contrib.idx2], pos[3 * contrib.idx2 + 1],
            // RDKit✔️✔️:                              pos[3 * contrib.idx2 + 2]);
            // RDKit✔️✔️:     const RDGeom::Point3D p3(pos[3 * contrib.idx3], pos[3 * contrib.idx3 + 1],
            // RDKit✔️✔️:                              pos[3 * contrib.idx3 + 2]);
            // RDKit✔️✔️:     const RDGeom::Point3D p4(pos[3 * contrib.idx4], pos[3 * contrib.idx4 + 1],
            // RDKit✔️✔️:                              pos[3 * contrib.idx4 + 2]);
            // RDKit✔️✔️:     const double cosY = Utils::calculateCosY(p1, p2, p3, p4);
            // RDKit✔️✔️:     const double sinYSq = 1.0 - cosY * cosY;
            // RDKit✔️✔️:     const double sinY = ((sinYSq > 0.0) ? sqrt(sinYSq) : 0.0);
            // RDKit✔️✔️:     // cos(2 * W) = 2 * cos(W) * cos(W) - 1 = 2 * sin(W) * sin(W) - 1
            // RDKit✔️✔️:     const double cos2W = 2.0 * sinY * sinY - 1.0;
            // RDKit✔️✔️:     accum += contrib.forceConstant *
            // RDKit✔️✔️:              (contrib.C0 + contrib.C1 * sinY + contrib.C2 * cos2W);
            // Existing InversionContrib::get_energy owns the equivalent
            // p1–p4/helper/arithmetic path, so the packed loop does not copy or
            // drift a second energy formula.
            accum += contribution.get_energy(context);
            // RDKit✔️✔️:   }
        }
        // RDKit✔️✔️:   return accum;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION ForceFields::UFF::InversionContribs::getEnergy
        accum
    }

    pub(super) fn get_grad(&self, context: &mut EvaluationContext<'_>, gradient: &mut [f64]) {
        // BEGIN RDKIT CPP FUNCTION ForceFields::UFF::InversionContribs::getGrad (ForceField/UFF/Inversions.cpp:69-135)
        // RDKit❗✔️: void InversionContribs::getGrad(double *pos, double *grad) const {
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
        // RDKit❗✔️:     const RDGeom::Point3D p4(pos[3 * contrib.idx4], pos[3 * contrib.idx4 + 1],
        // RDKit❗✔️:                              pos[3 * contrib.idx4 + 2]);
        // RDKit❗✔️:     double *g1 = &(grad[3 * contrib.idx1]);
        // RDKit❗✔️:     double *g2 = &(grad[3 * contrib.idx2]);
        // RDKit❗✔️:     double *g3 = &(grad[3 * contrib.idx3]);
        // RDKit❗✔️:     double *g4 = &(grad[3 * contrib.idx4]);
        // RDKit❗✔️:     RDGeom::Point3D rJI = p1 - p2;
        // RDKit❗✔️:     RDGeom::Point3D rJK = p3 - p2;
        // RDKit❗✔️:     RDGeom::Point3D rJL = p4 - p2;
        // RDKit❗✔️:     const double dJI = rJI.length();
        // RDKit❗✔️:     const double dJK = rJK.length();
        // RDKit❗✔️:     const double dJL = rJL.length();
        // RDKit❗✔️:     if (isDoubleZero(dJI) || isDoubleZero(dJK) || isDoubleZero(dJL)) {
        // RDKit❗✔️:       return;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     rJI.normalize();
        // RDKit❗✔️:     rJK.normalize();
        // RDKit❗✔️:     rJL.normalize();
        // RDKit❗✔️:
        // RDKit❗✔️:     RDGeom::Point3D n = (-rJI).crossProduct(rJK);
        // RDKit❗✔️:     n.normalize();
        // RDKit❗✔️:     double cosY = n.dotProduct(rJL);
        // RDKit❗✔️:     cosY = std::clamp(cosY, -1.0, 1.0);
        // RDKit❗✔️:     const double sinYSq = 1.0 - cosY * cosY;
        // RDKit❗✔️:     const double sinY = std::max(sqrt(sinYSq), 1.0e-8);
        // RDKit❗✔️:     double cosTheta = rJI.dotProduct(rJK);
        // RDKit❗✔️:     cosTheta = std::clamp(cosTheta, -1.0, 1.0);
        // RDKit❗✔️:     const double sinThetaSq = 1.0 - cosTheta * cosTheta;
        // RDKit❗✔️:     const double sinTheta = std::max(sqrt(sinThetaSq), 1.0e-8);
        // RDKit❗✔️:     // sin(2 * W) = 2 * sin(W) * cos(W) = 2 * cos(Y) * sin(Y)
        // RDKit❗✔️:     const double dE_dW = -contrib.forceConstant *
        // RDKit❗✔️:                          (contrib.C1 * cosY - 4.0 * contrib.C2 * cosY * sinY);
        // RDKit❗✔️:     const RDGeom::Point3D t1 = rJL.crossProduct(rJK);
        // RDKit❗✔️:     const RDGeom::Point3D t2 = rJI.crossProduct(rJL);
        // RDKit❗✔️:     const RDGeom::Point3D t3 = rJK.crossProduct(rJI);
        // RDKit❗✔️:     const double term1 = sinY * sinTheta;
        // RDKit❗✔️:     const double term2 = cosY / (sinY * sinThetaSq);
        // RDKit❗✔️:     const double tg1[3] = {
        // RDKit❗✔️:         (t1.x / term1 - (rJI.x - rJK.x * cosTheta) * term2) / dJI,
        // RDKit❗✔️:         (t1.y / term1 - (rJI.y - rJK.y * cosTheta) * term2) / dJI,
        // RDKit❗✔️:         (t1.z / term1 - (rJI.z - rJK.z * cosTheta) * term2) / dJI};
        // RDKit❗✔️:     const double tg3[3] = {
        // RDKit❗✔️:         (t2.x / term1 - (rJK.x - rJI.x * cosTheta) * term2) / dJK,
        // RDKit❗✔️:         (t2.y / term1 - (rJK.y - rJI.y * cosTheta) * term2) / dJK,
        // RDKit❗✔️:         (t2.z / term1 - (rJK.z - rJI.z * cosTheta) * term2) / dJK};
        // RDKit❗✔️:     const double tg4[3] = {(t3.x / term1 - rJL.x * cosY / sinY) / dJL,
        // RDKit❗✔️:                            (t3.y / term1 - rJL.y * cosY / sinY) / dJL,
        // RDKit❗✔️:                            (t3.z / term1 - rJL.z * cosY / sinY) / dJL};
        // RDKit❗✔️:     for (unsigned int i = 0; i < 3; ++i) {
        // RDKit❗✔️:       g1[i] += dE_dW * tg1[i];
        // RDKit❗✔️:       g2[i] += -dE_dW * (tg1[i] + tg3[i] + tg4[i]);
        // RDKit❗✔️:       g3[i] += dE_dW * tg3[i];
        // RDKit❗✔️:       g4[i] += dE_dW * tg4[i];
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION ForceFields::UFF::InversionContribs::getGrad

        // The borrowed context/slice provide nonnull evaluation inputs without
        // retaining an owner pointer or adding source-absent extent checks.
        for contribution in &self.contribs {
            // The canonical term helper performs the anchored p/g loads, zero
            // guard and analytic formula; false carries source's collection-
            // level return to this owner after preserving earlier term writes.
            if !contribution.get_grad(context, gradient) {
                return;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{InversionContribs, InversionContributionError};
    use crate::kernel::EvaluationContext;
    use crate::uff::inversion::InversionIndexArgument;

    fn add_term(
        contribs: &mut InversionContribs,
        positions: &mut [Vec<f64>],
        indices: [u32; 4],
        atomic_num: i32,
        is_c_bound_to_o: bool,
        scale: f64,
    ) -> Result<(), InversionContributionError> {
        let position_refs = positions
            .iter_mut()
            .map(Vec::as_mut_slice)
            .collect::<Vec<_>>();
        contribs.add_contrib(
            &position_refs,
            indices[0],
            indices[1],
            indices[2],
            indices[3],
            atomic_num,
            is_c_bound_to_o,
            scale,
        )
    }

    fn energy(contribs: &InversionContribs, coordinates: &[f64], num_points: u32) -> f64 {
        let matrix_len = (num_points * (num_points + 1) / 2) as usize;
        let mut distance_cache = vec![0.0; matrix_len];
        let mut context = EvaluationContext::for_test(coordinates, &mut distance_cache, num_points);
        contribs.get_energy(&mut context)
    }

    fn packed_gradient(
        contribs: &InversionContribs,
        coordinates: &[f64],
        num_points: u32,
        initial_gradient: &[f64],
    ) -> Vec<f64> {
        let matrix_len = (num_points * (num_points + 1) / 2) as usize;
        let mut distance_cache = vec![0.0; matrix_len];
        let mut context = EvaluationContext::for_test(coordinates, &mut distance_cache, num_points);
        let mut gradient = initial_gradient.to_vec();
        contribs.get_grad(&mut context, &mut gradient);
        gradient
    }

    fn assert_close(actual: f64, expected: f64) {
        let tolerance = 1.0e-12 * expected.abs().max(1.0);
        assert!(
            (actual - expected).abs() <= tolerance,
            "actual {actual:.17e} differs from fixed source value {expected:.17e} by more than {tolerance:.3e}"
        );
    }

    fn assert_gradient_close(actual: &[f64], expected: &[f64]) {
        assert_eq!(actual.len(), expected.len());
        for (&actual, &expected) in actual.iter().zip(expected) {
            assert_close(actual, expected);
        }
    }

    #[test]
    fn cf3d_u16_empty_packed_terms_return_source_positive_zero_without_reads() {
        // RDKit source: Inversions.cpp:45-64 initializes accum to +0 and
        // dereferences coordinates only inside the stored-term loop.
        let contribs = InversionContribs::default();
        let actual = energy(&contribs, &[], 0);
        assert_eq!(actual.to_bits(), 0.0_f64.to_bits());
    }

    #[test]
    fn cf3d_u16_packed_energy_covers_all_source_parameter_options() {
        // Fixed RDKit source values for Inversions.cpp:27-64 and
        // Utils.cpp:42-85, observed through the packed owner. Both geometries,
        // both oxygen-flag values, every source atomic-number branch and the
        // default branch are covered. Scale has no source guard or branch.
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
        let cases = [
            (6, 0.0, 0.8, 20.0 / 3.0),
            (7, 0.0, 0.8, 20.0 / 3.0),
            (8, 0.0, 0.8, 20.0 / 3.0),
            (
                15,
                7.333_333_333_333_334,
                2.275_445_586_749_68,
                2.275_445_586_749_68,
            ),
            (
                33,
                7.333_333_333_333_333,
                2.447_437_894_208_855,
                2.447_437_894_208_855,
            ),
            (
                51,
                7.333_333_333_333_334,
                2.495_185_335_428_339_5,
                2.495_185_335_428_339_5,
            ),
            (83, 7.333_333_333_333_334, 2.64, 2.64),
            (
                16,
                7.333_333_333_333_333,
                50_543_252.975_452_304,
                50_543_252.975_452_304,
            ),
        ];

        for (atomic_num, planar_expected, nonplanar_expected, oxygen_expected) in cases {
            for is_c_bound_to_o in [false, true] {
                for scale in [1.0, 2.5, 0.0, -2.0] {
                    let expected_nonplanar = if (6..=8).contains(&atomic_num) && is_c_bound_to_o {
                        oxygen_expected
                    } else {
                        nonplanar_expected
                    };
                    let mut positions = vec![vec![0.0; 3]; 4];
                    let mut contribs = InversionContribs::default();
                    add_term(
                        &mut contribs,
                        &mut positions,
                        [0, 1, 2, 3],
                        atomic_num,
                        is_c_bound_to_o,
                        scale,
                    )
                    .unwrap();

                    assert_close(energy(&contribs, &planar, 4), scale * planar_expected);
                    assert_close(energy(&contribs, &nonplanar, 4), scale * expected_nonplanar);
                }
            }
        }
    }

    #[test]
    fn cf3d_u16_packed_terms_preserve_source_append_errors_and_duplicate_ids() {
        // RDKit source: Inversions.cpp:27-43 checks idx1 through idx4 before
        // emplacing. URANGE_CHECK does not reject repeated valid indices.
        let mut positions = vec![vec![0.0; 3]; 4];
        let cases = [
            ([4, 4, 4, 4], InversionIndexArgument::First, "idx1"),
            ([0, 4, 4, 4], InversionIndexArgument::Second, "idx2"),
            ([0, 1, 4, 4], InversionIndexArgument::Third, "idx3"),
            ([0, 1, 2, 4], InversionIndexArgument::Fourth, "idx4"),
        ];
        let mut contribs = InversionContribs::default();
        let position_refs = positions
            .iter_mut()
            .map(Vec::as_mut_slice)
            .collect::<Vec<_>>();

        for (indices, expected_argument, expected_source_name) in cases {
            let error = contribs
                .add_contrib(
                    &position_refs,
                    indices[0],
                    indices[1],
                    indices[2],
                    indices[3],
                    6,
                    false,
                    1.0,
                )
                .unwrap_err();
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
            assert!(contribs.contribs.is_empty());
        }

        contribs
            .add_contrib(&position_refs, 1, 1, 1, 1, 6, false, 1.0)
            .unwrap();
        assert_eq!(contribs.contribs.len(), 1);
    }

    #[test]
    fn cf3d_u16_shared_center_energy_keeps_exact_insertion_accumulation_order() {
        // RDKit source: Inversions.cpp:45-64 iterates d_contribs without
        // sorting. Each source term shares center atom 1 and has sinY=0,
        // giving exact fixed energies +1e16, -1e16, and +1. Left-folding in
        // insertion order returns 1; a reordered sum loses the final unit.
        let mut positions = vec![vec![0.0; 3]; 10];
        let mut contribs = InversionContribs::default();
        add_term(
            &mut contribs,
            &mut positions,
            [0, 1, 2, 3],
            6,
            false,
            5.0e15,
        )
        .unwrap();
        add_term(
            &mut contribs,
            &mut positions,
            [4, 1, 5, 6],
            6,
            false,
            -5.0e15,
        )
        .unwrap();
        add_term(&mut contribs, &mut positions, [7, 1, 8, 9], 6, false, 0.5).unwrap();

        let coordinates = [
            1.0, 0.0, 0.0, // atom 0: term 1 p1
            0.0, 0.0, 0.0, // atom 1: shared center
            0.0, 1.0, 0.0, // atom 2: term 1 p3
            0.0, 0.0, 1.0, // atom 3: term 1 p4
            1.0, 0.0, 0.0, // atom 4: term 2 p1
            0.0, 1.0, 0.0, // atom 5: term 2 p3
            0.0, 0.0, 1.0, // atom 6: term 2 p4
            1.0, 0.0, 0.0, // atom 7: term 3 p1
            0.0, 1.0, 0.0, // atom 8: term 3 p3
            0.0, 0.0, 1.0, // atom 9: term 3 p4
        ];
        assert_eq!(energy(&contribs, &coordinates, 10), 1.0);
    }

    #[test]
    fn cf3d_u17_empty_packed_gradient_returns_without_coordinate_reads() {
        // Inversions.cpp::getGrad has no empty special case; the source loop
        // naturally avoids every position and gradient index read.
        let contribs = InversionContribs::default();
        let initial = [1.25, -2.0, 3.5];
        let actual = packed_gradient(&contribs, &[], 0, &initial);
        assert_eq!(actual, initial);
    }

    #[test]
    fn cf3d_u17_packed_gradient_covers_all_source_parameter_options() {
        // Fixed values from pinned Inversion.cpp::getGrad and Utils.cpp. The
        // parameter matrix covers every central-atom branch, both oxygen flag
        // values, and both scales retained by the single-term source tests.
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
                    let mut positions = vec![vec![0.0; 3]; 4];
                    let mut contribs = InversionContribs::default();
                    add_term(
                        &mut contribs,
                        &mut positions,
                        [0, 1, 2, 3],
                        atomic_num,
                        is_c_bound_to_o,
                        scale,
                    )
                    .unwrap();
                    let actual = packed_gradient(&contribs, &coordinates, 4, &[0.0; 12]);
                    assert_close(actual[0], expected * scale);
                }
            }
        }

        // Full independent source vectors check every Cartesian component and
        // all four atom updates for both carbon force-constant families.
        let unbound = [
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
        ];
        let oxygen_bound = [
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
        ];
        for (is_c_bound_to_o, expected) in [(false, unbound), (true, oxygen_bound)] {
            let mut positions = vec![vec![0.0; 3]; 4];
            let mut contribs = InversionContribs::default();
            add_term(
                &mut contribs,
                &mut positions,
                [0, 1, 2, 3],
                6,
                is_c_bound_to_o,
                1.0,
            )
            .unwrap();
            assert_gradient_close(
                &packed_gradient(&contribs, &coordinates, 4, &[0.0; 12]),
                &expected,
            );
        }
    }

    #[test]
    fn cf3d_u17_each_source_zero_length_branch_returns_from_the_packed_loop() {
        // The source return is collection-scoped: retain the first term's
        // gradient, make no current-term writes, and skip the valid third term.
        let source_geometry = [
            0.4, 1.2, -0.3, // p1
            -0.2, 0.1, 0.5, // p2
            1.7, -0.4, 0.2, // p3
            0.6, 1.1, 1.4, // p4
        ];
        let first_term = [
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
        ];
        let mut expected = vec![0.25; 36];
        for (actual, delta) in expected[..12].iter_mut().zip(first_term) {
            *actual += delta;
        }

        let zero_geometries = [
            [
                0.0, 0.0, 0.0, // p1 == p2: dJI is zero
                0.0, 0.0, 0.0, // p2
                1.0, 0.0, 0.0, // p3
                0.0, 1.0, 0.0, // p4
            ],
            [
                1.0, 0.0, 0.0, // p1
                0.0, 0.0, 0.0, // p2 == p3: dJK is zero
                0.0, 0.0, 0.0, // p3
                0.0, 1.0, 0.0, // p4
            ],
            [
                1.0, 0.0, 0.0, // p1
                0.0, 0.0, 0.0, // p2 == p4: dJL is zero
                0.0, 1.0, 0.0, // p3
                0.0, 0.0, 0.0, // p4
            ],
        ];

        for zero_geometry in zero_geometries {
            let mut positions = vec![vec![0.0; 3]; 12];
            let mut contribs = InversionContribs::default();
            for indices in [[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11]] {
                add_term(&mut contribs, &mut positions, indices, 6, false, 1.0).unwrap();
            }
            let mut coordinates = Vec::with_capacity(36);
            coordinates.extend_from_slice(&source_geometry);
            coordinates.extend_from_slice(&zero_geometry);
            coordinates.extend_from_slice(&source_geometry);

            assert_gradient_close(
                &packed_gradient(&contribs, &coordinates, 12, &[0.25; 36]),
                &expected,
            );
        }
    }

    #[test]
    fn cf3d_u17_shared_center_gradient_accumulates_in_source_term_order() {
        // Both terms use the same source geometry and atom 1 as center. Fixed
        // source values verify += across contributions with different scales.
        let p1 = [0.4, 1.2, -0.3];
        let p2 = [-0.2, 0.1, 0.5];
        let p3 = [1.7, -0.4, 0.2];
        let p4 = [0.6, 1.1, 1.4];
        let coordinates = [
            p1[0], p1[1], p1[2], // atom 0: term 1 p1
            p2[0], p2[1], p2[2], // atom 1: shared p2
            p3[0], p3[1], p3[2], // atom 2: term 1 p3
            p4[0], p4[1], p4[2], // atom 3: term 1 p4
            p1[0], p1[1], p1[2], // atom 4: term 2 p1
            p3[0], p3[1], p3[2], // atom 5: term 2 p3
            p4[0], p4[1], p4[2], // atom 6: term 2 p4
        ];
        let mut positions = vec![vec![0.0; 3]; 7];
        let mut contribs = InversionContribs::default();
        add_term(&mut contribs, &mut positions, [0, 1, 2, 3], 6, false, 1.0).unwrap();
        add_term(&mut contribs, &mut positions, [4, 1, 5, 6], 6, false, 2.5).unwrap();

        let g1 = [
            -0.262_966_359_151_671_3,
            -0.482_705_371_593_478_9,
            -0.860_944_655_304_787_5,
        ];
        let g2 = [
            0.979_421_110_870_018_2,
            0.984_363_445_830_566_3,
            0.221_466_459_787_311_05,
        ];
        let g3 = [
            -0.089_443_630_611_437_27,
            -0.164_184_198_656_610_38,
            -0.292_835_996_111_417_3,
        ];
        let g4 = [
            -0.627_011_121_106_909_5,
            -0.337_473_875_580_477_15,
            0.932_314_191_628_893_8,
        ];
        let expected = [
            1.25 + g1[0],
            1.25 + g1[1],
            1.25 + g1[2],
            1.25 + g2[0] + 2.5 * g2[0],
            1.25 + g2[1] + 2.5 * g2[1],
            1.25 + g2[2] + 2.5 * g2[2],
            1.25 + g3[0],
            1.25 + g3[1],
            1.25 + g3[2],
            1.25 + g4[0],
            1.25 + g4[1],
            1.25 + g4[2],
            1.25 + 2.5 * g1[0],
            1.25 + 2.5 * g1[1],
            1.25 + 2.5 * g1[2],
            1.25 + 2.5 * g3[0],
            1.25 + 2.5 * g3[1],
            1.25 + 2.5 * g3[2],
            1.25 + 2.5 * g4[0],
            1.25 + 2.5 * g4[1],
            1.25 + 2.5 * g4[2],
        ];
        assert_gradient_close(
            &packed_gradient(&contribs, &coordinates, 7, &[1.25; 21]),
            &expected,
        );
    }
}
