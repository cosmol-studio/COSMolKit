//! Public UFF parity profiles; no force-field implementation or owner dependency.
use crate::registry::{Input, Operation, Record, SmilesCase, Value};
use cosmolkit::{
    Conformer3D, CoordinateBlock, CoordinateDimension, Molecule, SdfCoordinateMode, SdfReadParams,
    UffConformerOptimizationParams, UffOptimizationParams,
};
use serde::{Deserialize, Serialize};
use std::error::Error;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Profile {
    Coverage {
        add_hydrogens: bool,
    },
    Optimization {
        add_hydrogens: bool,
        max_iterations: i32,
        vdw_threshold: u32,
        ignore_interfragment_interactions: bool,
        conformer_id: Option<usize>,
    },
    ConformerOptimization {
        add_hydrogens: bool,
        max_iterations: i32,
        vdw_threshold: u32,
        ignore_interfragment_interactions: bool,
        conformer_count: usize,
    },
}

pub fn profiles(operation: Operation) -> Vec<Profile> {
    match operation {
        Operation::UffCoverage => [false, true]
            .map(|add_hydrogens| Profile::Coverage { add_hydrogens })
            .into(),
        Operation::UffOptimization => [Profile::Optimization {
            add_hydrogens: true,
            max_iterations: 2,
            vdw_threshold: 100,
            ignore_interfragment_interactions: true,
            conformer_id: None,
        }]
        .into(),
        Operation::UffConformerOptimization => [Profile::ConformerOptimization {
            add_hydrogens: true,
            max_iterations: 2,
            vdw_threshold: 100,
            ignore_interfragment_interactions: true,
            conformer_count: 2,
        }]
        .into(),
        _ => Vec::new(),
    }
}

/// Common input generated before global preflight, NOT an expected CK output.
/// Both libraries read this quantized MolBlock and the same bit-exact rows.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Geometry {
    pub molblock: String,
    pub atom_count: usize,
    #[serde(default)]
    pub coordinate_rows: Vec<CoordinateRow>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CoordinateRow {
    pub conformer_id: usize,
    pub xyz_bits: Vec<[u64; 3]>,
}

/// A recorded source rejection is data; absent preparation is not ready.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum GeometryPreparation {
    Ready(Geometry),
    Rejected {
        stage: crate::molecular::Stage,
        detail: String,
    },
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct UffInput {
    pub case: SmilesCase,
    pub profile: Profile,
    pub preparation: Option<GeometryPreparation>,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum ExpectedErrorReason {
    EmbeddingRejected { case_id: String },
    SourceTbpCenterParamsMissing { center_atom_index: usize },
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Observation {
    Coverage(bool),
    Optimized {
        status: i32,
        energy_bits: u64,
        xyz_bits: Vec<[u64; 3]>,
    },
    OptimizedConformers {
        conformers: Vec<ConformerObservation>,
    },
    Error {
        stage: crate::molecular::Stage,
        detail: String,
        /// Original references retain their raw diagnostics. Actual failures
        /// additionally retain the narrow structured cause, never a success.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        reason: Option<ExpectedErrorReason>,
    },
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ConformerObservation {
    pub conformer_id: usize,
    pub status: i32,
    pub energy_bits: u64,
    pub xyz_bits: Vec<[u64; 3]>,
}

pub fn validate_reference(
    recipe: &UffInput,
    prepared: &UffInput,
    output: &Observation,
) -> Result<(), String> {
    if recipe.case != prepared.case
        || recipe.profile != prepared.profile
        || recipe
            .preparation
            .as_ref()
            .is_some_and(|value| prepared.preparation.as_ref() != Some(value))
    {
        return Err("UFF reference changed recipe identity".into());
    }
    match (recipe.profile, &prepared.preparation, output) {
        (Profile::Coverage { .. }, None, Observation::Coverage(_)) => Ok(()),
        (Profile::Coverage { .. }, None, Observation::Error { detail, .. })
            if !detail.is_empty() =>
        {
            Ok(())
        }
        (
            Profile::Optimization { .. },
            Some(GeometryPreparation::Ready(geometry)),
            Observation::Optimized {
                status,
                energy_bits,
                xyz_bits,
            },
        ) if geometry.atom_count > 0
            && geometry.molblock.contains("M  END")
            && geometry.coordinate_rows.is_empty()
            && matches!(status, 0 | 1)
            && f64::from_bits(*energy_bits).is_finite()
            && xyz_bits.len() == geometry.atom_count
            && xyz_bits
                .iter()
                .flatten()
                .all(|bits| f64::from_bits(*bits).is_finite()) =>
        {
            Ok(())
        }
        (
            Profile::ConformerOptimization {
                conformer_count, ..
            },
            Some(GeometryPreparation::Ready(geometry)),
            Observation::OptimizedConformers { conformers },
        ) if geometry.atom_count > 0
            && geometry.molblock.contains("M  END")
            && matches!(conformer_count, 1 | 2 | 3)
            && valid_coordinate_rows(geometry, conformer_count)
            && conformers.len() == conformer_count
            && conformers.iter().enumerate().all(|(index, row)| {
                row.conformer_id == [7, 3, 11][index]
                    && matches!(row.status, 0 | 1)
                    && f64::from_bits(row.energy_bits).is_finite()
                    && row.xyz_bits.len() == geometry.atom_count
                    && row
                        .xyz_bits
                        .iter()
                        .flatten()
                        .all(|bits| f64::from_bits(*bits).is_finite())
            }) =>
        {
            Ok(())
        }
        (
            Profile::ConformerOptimization {
                conformer_count, ..
            },
            Some(GeometryPreparation::Ready(geometry)),
            Observation::Error { detail, .. },
        ) if geometry.atom_count > 0
            && geometry.molblock.contains("M  END")
            && valid_coordinate_rows(geometry, conformer_count)
            && !detail.is_empty() =>
        {
            Ok(())
        }
        (
            Profile::Optimization { .. },
            Some(GeometryPreparation::Ready(geometry)),
            Observation::Error { detail, .. },
        ) if geometry.atom_count > 0
            && geometry.molblock.contains("M  END")
            && geometry.coordinate_rows.is_empty()
            && !detail.is_empty() =>
        {
            Ok(())
        }
        (
            Profile::Optimization { .. },
            Some(GeometryPreparation::Rejected { stage, detail }),
            Observation::Error {
                stage: output_stage,
                detail: output_detail,
                ..
            },
        ) if matches!(
            stage,
            crate::molecular::Stage::Parse | crate::molecular::Stage::Preparation
        ) && stage == output_stage
            && !detail.is_empty()
            && detail == output_detail =>
        {
            Ok(())
        }
        (
            Profile::ConformerOptimization { .. },
            Some(GeometryPreparation::Rejected { stage, detail }),
            Observation::Error {
                stage: output_stage,
                detail: output_detail,
                ..
            },
        ) if matches!(
            stage,
            crate::molecular::Stage::Parse | crate::molecular::Stage::Preparation
        ) && stage == output_stage
            && !detail.is_empty()
            && detail == output_detail =>
        {
            Ok(())
        }
        _ => Err("UFF reference lacks valid common geometry or typed observation".into()),
    }
}

pub(crate) fn valid_coordinate_rows(geometry: &Geometry, expected_count: usize) -> bool {
    let ids = [7, 3, 11];
    matches!(expected_count, 1 | 2 | 3)
        && geometry.coordinate_rows.len() == expected_count
        && geometry
            .coordinate_rows
            .iter()
            .enumerate()
            .all(|(index, row)| {
                row.conformer_id == ids[index]
                    && row.xyz_bits.len() == geometry.atom_count
                    && row
                        .xyz_bits
                        .iter()
                        .flatten()
                        .all(|bits| f64::from_bits(*bits).is_finite())
            })
}

/// Classify only the existing embedding and source TBP parameter errors.
/// Changing iteration counts does not broaden this set; raw details are retained.
fn reference_error_reason(
    input: &UffInput,
    stage: &crate::molecular::Stage,
    detail: &str,
) -> Option<ExpectedErrorReason> {
    use crate::molecular::Stage;
    let operation = match input.profile {
        Profile::Optimization { .. } => Operation::UffOptimization,
        Profile::ConformerOptimization { .. } => Operation::UffConformerOptimization,
        Profile::Coverage { .. } => return None,
    };
    if !profiles(operation).contains(&input.profile) {
        return None;
    }
    match (&input.preparation, stage) {
        (
            Some(GeometryPreparation::Rejected {
                stage: rejected_stage,
                detail: rejected_detail,
            }),
            Stage::Preparation,
        ) if rejected_stage == stage
            && rejected_detail == detail
            && detail
                == format!(
                    "ValueError: UFF common geometry embedding failed: {}",
                    input.case.id
                ) =>
        {
            Some(ExpectedErrorReason::EmbeddingRejected {
                case_id: input.case.id.clone(),
            })
        }
        (Some(GeometryPreparation::Ready(_)), Stage::Operation) if input.case.id == "line:320" => {
            // Pinned 351f8f378f8ad6bbd517980c38896e66bf907af8:
            // Builder.cpp:324-327 checks both endpoints but passes params[atomIdx].
            // AngleBend.cpp:79: PRECONDITION(at2Params, "bad params pointer");
            // Independent p1 source trace identifies this original center as 1.
            // This is the approved original reference error, not a molecule patch.
            let lines: Vec<_> = detail.lines().map(str::trim).collect();
            if lines.len() == 6
                && lines[0] == "RuntimeError: Pre-condition Violation"
                && lines[1] == "bad params pointer"
                && lines[2]
                    == "Violation occurred on line 79 in file Code/ForceField/UFF/AngleBend.cpp"
                && lines[3] == "Failed Expression: at2Params"
                && lines[4] == "RDKIT: 2026.03.1"
                && lines[5] == "BOOST: 1_85"
            {
                Some(ExpectedErrorReason::SourceTbpCenterParamsMissing {
                    center_atom_index: 1,
                })
            } else {
                None
            }
        }
        _ => None,
    }
}

fn optimization_error_reason(
    error: &cosmolkit::OperationError,
    all_conformers: bool,
) -> Option<ExpectedErrorReason> {
    let cosmolkit::OperationError::UffOptimization(cause) = error else {
        return None;
    };
    let required_kind = if all_conformers {
        cosmolkit::UffOptimizationErrorKind::ConformerOptimization
    } else {
        cosmolkit::UffOptimizationErrorKind::Optimization
    };
    if cause.kind() != required_kind {
        return None;
    }
    // The public boundary preserves Error::source() but keeps builder types private.
    // Decode only the terminal concrete cause's existing Debug field representation;
    // do not inspect the outer display string or accept generic Construction errors.
    let mut leaf: &(dyn Error + 'static) = cause;
    while let Some(source) = leaf.source() {
        leaf = source;
    }
    let diagnostic = format!("{leaf:?}");
    let index = diagnostic
        .strip_prefix("SourceTbpCenterParamsMissing { center_atom_index: ")?
        .strip_suffix(" }")?
        .parse::<usize>()
        .ok()?;
    Some(ExpectedErrorReason::SourceTbpCenterParamsMissing {
        center_atom_index: index,
    })
}

pub fn matches(input: &Input, expected: &Observation, actual: &Observation) -> bool {
    match (expected, actual) {
        (
            Observation::Error { stage, detail, .. },
            Observation::Error {
                stage: actual_stage,
                detail: actual_detail,
                reason: Some(actual_reason),
            },
        ) => {
            let Input::Uff(row) = input else {
                return false;
            };
            if stage != actual_stage || validate_reference(row, row, expected).is_err() {
                return false;
            }
            let Some(expected_reason) = reference_error_reason(row, stage, detail) else {
                return false;
            };
            if &expected_reason != actual_reason {
                return false;
            }
            match expected_reason {
                ExpectedErrorReason::EmbeddingRejected { .. } => detail == actual_detail,
                ExpectedErrorReason::SourceTbpCenterParamsMissing { .. } => {
                    !actual_detail.is_empty()
                }
            }
        }
        (Observation::Error { .. }, _) | (_, Observation::Error { .. }) => false,
        _ => expected == actual,
    }
}

pub fn run(input: &Input) -> Result<Record, String> {
    let Input::Uff(row) = input else {
        return Err("expected UFF input".into());
    };
    let mut stage = crate::molecular::Stage::Parse;
    let mut reason = None;
    let result = (|| -> Result<Observation, String> {
        match row.profile {
            Profile::Coverage { add_hydrogens } => {
                let mol = Molecule::from_smiles(&row.case.smiles).map_err(|e| e.to_string())?;
                stage = crate::molecular::Stage::Preparation;
                let mol = if add_hydrogens {
                    mol.with_hydrogens()
                        .map_err(|e| e.to_string())?
                        .with_assigned_valence()
                        .map_err(|e| e.to_string())?
                } else {
                    mol
                };
                stage = crate::molecular::Stage::Operation;
                Ok(Observation::Coverage(
                    mol.uff_has_all_molecule_params()
                        .map_err(|e| e.to_string())?,
                ))
            }
            Profile::Optimization {
                max_iterations,
                vdw_threshold,
                ignore_interfragment_interactions,
                conformer_id,
                ..
            } => {
                let geometry = match row
                    .preparation
                    .as_ref()
                    .ok_or("UFF common geometry was not prepared")?
                {
                    GeometryPreparation::Ready(geometry) => geometry,
                    GeometryPreparation::Rejected {
                        stage: failed_stage,
                        detail,
                    } => {
                        // No optimization is called for this recorded rejection.
                        // Verify its precise preparation reason separately.
                        stage = failed_stage.clone();
                        reason = reference_error_reason(row, failed_stage, detail);
                        return Err(detail.clone());
                    }
                };
                stage = crate::molecular::Stage::Preparation;
                let mol = Molecule::from_sdf_with_params(
                    &geometry.molblock,
                    &SdfReadParams {
                        sanitize: true,
                        remove_hs: false,
                        coordinate_mode: SdfCoordinateMode::Require3D,
                        ..Default::default()
                    },
                )
                .map_err(|e| e.to_string())?
                .with_assigned_valence()
                .map_err(|e| e.to_string())?;
                if mol.num_atoms() != geometry.atom_count {
                    return Err("common geometry atom count changed".into());
                }
                stage = crate::molecular::Stage::Operation;
                let result = mol
                    .with_uff_optimized_with_params(&UffOptimizationParams {
                        max_iterations,
                        vdw_threshold: f64::from(vdw_threshold),
                        ignore_interfragment_interactions,
                        conformer_id,
                    })
                    .map_err(|e| {
                        reason = optimization_error_reason(&e, false);
                        e.to_string()
                    })?;
                let selected = conformer_id.unwrap_or(0);
                let conf = result
                    .molecule
                    .conformers_3d()
                    .iter()
                    .find(|conf| conf.id() == selected)
                    .ok_or("optimized conformer absent")?;
                Ok(Observation::Optimized {
                    status: result.status,
                    energy_bits: result.energy.to_bits(),
                    xyz_bits: conf
                        .coordinates()
                        .iter()
                        .map(|xyz| xyz.map(f64::to_bits))
                        .collect(),
                })
            }
            Profile::ConformerOptimization {
                max_iterations,
                vdw_threshold,
                ignore_interfragment_interactions,
                conformer_count,
                ..
            } => {
                let geometry = match row
                    .preparation
                    .as_ref()
                    .ok_or("UFF common conformer geometry was not prepared")?
                {
                    GeometryPreparation::Ready(geometry) => geometry,
                    GeometryPreparation::Rejected {
                        stage: failed_stage,
                        detail,
                    } => {
                        stage = failed_stage.clone();
                        reason = reference_error_reason(row, failed_stage, detail);
                        return Err(detail.clone());
                    }
                };
                if !valid_coordinate_rows(geometry, conformer_count) {
                    return Err("common UFF conformer coordinate rows are invalid".into());
                }
                stage = crate::molecular::Stage::Preparation;
                let base = Molecule::from_sdf_with_params(
                    &geometry.molblock,
                    &SdfReadParams {
                        sanitize: true,
                        remove_hs: false,
                        coordinate_mode: SdfCoordinateMode::Require3D,
                        ..Default::default()
                    },
                )
                .map_err(|e| e.to_string())?;
                if base.num_atoms() != geometry.atom_count {
                    return Err("common geometry atom count changed".into());
                }
                let coordinates = CoordinateBlock {
                    conformers_3d: geometry
                        .coordinate_rows
                        .iter()
                        .map(|row| {
                            Conformer3D::new(
                                row.conformer_id,
                                row.xyz_bits
                                    .iter()
                                    .map(|xyz| xyz.map(f64::from_bits))
                                    .collect(),
                                true,
                            )
                        })
                        .collect(),
                    source_coordinate_dim: Some(CoordinateDimension::ThreeD),
                    ..Default::default()
                };
                let mol = Molecule::from_parts(
                    base.topology().clone(),
                    coordinates,
                    base.properties().clone(),
                )
                .map_err(|e| e.to_string())?
                .with_assigned_valence()
                .map_err(|e| e.to_string())?;
                stage = crate::molecular::Stage::Operation;
                let result = mol
                    .with_uff_optimized_conformers_with_params(&UffConformerOptimizationParams {
                        num_threads: 1,
                        max_iterations,
                        vdw_threshold: f64::from(vdw_threshold),
                        ignore_interfragment_interactions,
                    })
                    .map_err(|e| {
                        reason = optimization_error_reason(&e, true);
                        e.to_string()
                    })?;
                let stored_coordinates = result.molecule.conformers_3d();
                if result.conformers.len() != conformer_count
                    || stored_coordinates.len() != conformer_count
                {
                    return Err("all-conformer result count changed".into());
                }
                let conformers = result
                    .conformers
                    .iter()
                    .zip(stored_coordinates)
                    .map(|(outcome, conformer)| {
                        if outcome.conformer_id != conformer.id() {
                            return Err("all-conformer result order changed".into());
                        }
                        Ok(ConformerObservation {
                            conformer_id: outcome.conformer_id,
                            status: outcome.status,
                            energy_bits: outcome.energy.to_bits(),
                            xyz_bits: conformer
                                .coordinates()
                                .iter()
                                .map(|xyz| xyz.map(f64::to_bits))
                                .collect(),
                        })
                    })
                    .collect::<Result<Vec<_>, String>>()?;
                Ok(Observation::OptimizedConformers { conformers })
            }
        }
    })();
    Ok(Record {
        input: input.clone(),
        output: Value::Uff(result.unwrap_or_else(|detail| Observation::Error {
            stage,
            detail,
            reason,
        })),
    })
}
