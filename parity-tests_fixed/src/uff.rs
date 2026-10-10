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

/// Keep the original SMILES topology when MOL's bond type 8 introduces an Any
/// query carrier. Coordinates still come from the same quantized common block;
/// no query predicate is lowered into a concrete molecule.
pub(crate) fn read_geometry(
    case: &SmilesCase,
    geometry: &Geometry,
    add_hydrogens: bool,
) -> Result<Molecule, String> {
    let record = cosmolkit::SdfRecord::from_sdf_with_params(
        &geometry.molblock,
        &SdfReadParams {
            sanitize: true,
            remove_hs: false,
            coordinate_mode: SdfCoordinateMode::Require3D,
            ..Default::default()
        },
    )
    .map_err(|error| error.to_string())?;
    if let Ok(mol) = record.molecule() {
        return Ok(mol.clone());
    }
    let query = record.query_graph().map_err(|error| error.to_string())?;
    let base = Molecule::from_smiles(&case.smiles).map_err(|error| error.to_string())?;
    let base = if add_hydrogens {
        base.with_hydrogens().map_err(|error| error.to_string())?
    } else {
        base
    };
    if base.num_atoms() != query.num_atoms()
        || base.bonds().len() != query.num_bonds()
        || base.atoms().iter().zip(query.atoms()).any(|(a, b)| {
            a.atomic_number() != b.atomic_number() || a.formal_charge() != b.formal_charge()
        })
        || base.bonds().iter().zip(query.bonds()).any(|(a, b)| {
            a.begin() != b.begin() || a.end() != b.end() || a.order() != b.bond().order()
        })
    {
        return Err("common MOL query carrier changed the original SMILES topology".into());
    }
    if query.conformers_3d().len() != 1 {
        return Err("common MOL query carrier has no unique XYZ conformer".into());
    }
    Molecule::from_parts(
        base.topology().clone(),
        CoordinateBlock {
            conformers_3d: query.conformers_3d().to_vec(),
            source_coordinate_dim: Some(CoordinateDimension::ThreeD),
            ..Default::default()
        },
        base.properties().clone(),
    )
    .map_err(|error| error.to_string())
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CoordinateRow {
    pub conformer_id: usize,
    pub xyz_bits: Vec<[u64; 3]>,
}

/// The user-selected reference deadline is distinct from a chemistry error.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PreparationTimeout {
    pub limit_seconds: u64,
    pub mechanism: TimeoutMechanism,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum TimeoutMechanism {
    Native,
    ProcessDeadline,
}

/// A recorded source rejection is data; absent preparation is not ready.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum GeometryPreparation {
    Ready(Geometry),
    TimedOut(PreparationTimeout),
    Rejected {
        stage: crate::molecular::Stage,
        detail: String,
    },
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum PreparationRejection {
    Parse,
    MolWriterAtomicNumberNotFound,
}

pub(crate) fn reproduce_preparation_rejection(
    case: &SmilesCase,
    preparation: &GeometryPreparation,
) -> Result<Option<(PreparationRejection, String)>, String> {
    let GeometryPreparation::Rejected { stage, detail } = preparation else {
        return Ok(None);
    };
    if *stage == crate::molecular::Stage::Parse {
        return match Molecule::from_smiles(&case.smiles) {
            Err(error @ cosmolkit::SmilesError::Construction(_)) => Err(error.to_string()),
            Err(error) => Ok(Some((PreparationRejection::Parse, error.to_string()))),
            Ok(_) => Err("CK accepted the common-input parse rejected by RDKit".into()),
        };
    }
    if *stage == crate::molecular::Stage::Preparation && atomic_number_diagnostic(detail) {
        let mol = Molecule::from_smiles(&case.smiles)
            .map_err(|error| error.to_string())?
            .with_hydrogens()
            .map_err(|error| error.to_string())?;
        return match mol.to_mol_with_params(&cosmolkit::MolBlockWriteParams {
            format: cosmolkit::SdfFormat::V3000,
            ..Default::default()
        }) {
            Err(cosmolkit::MolecularIoError::MolWrite(cosmolkit::MolWriteError::Value(
                message,
            ))) if message == "Atomic number not found" => Ok(Some((
                PreparationRejection::MolWriterAtomicNumberNotFound,
                message,
            ))),
            Err(error) => Err(format!("common-input writer rejection differs: {error}")),
            Ok(_) => Err("CK accepted the common-input writer rejected by RDKit".into()),
        };
    }
    Ok(None)
}

pub(crate) fn atomic_number_diagnostic(detail: &str) -> bool {
    detail.lines().map(str::trim).eq([
        "RuntimeError: Pre-condition Violation",
        "Atomic number not found",
        "Violation occurred on line 159 in file Code/GraphMol/PeriodicTable.h",
        "Failed Expression: atomicNumber < byanum.size()",
        "RDKIT: 2026.03.6",
        "BOOST: 1_85",
    ])
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
    SourceTbpCenterParamsMissing {
        center_atom_index: usize,
    },
    PreparationRejected {
        reason: PreparationRejection,
        detail: String,
    },
    ParseRejected {
        detail: String,
    },
    TimedOut(PreparationTimeout),
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
        (
            Profile::Optimization { .. } | Profile::ConformerOptimization { .. },
            Some(GeometryPreparation::Ready(geometry)),
            Observation::SourceTbpCenterParamsMissing { center_atom_index },
        ) if *center_atom_index < geometry.atom_count => Ok(()),
        (
            _,
            Some(GeometryPreparation::Rejected {
                stage,
                detail: source_detail,
            }),
            Observation::PreparationRejected { reason, detail },
        ) if !detail.is_empty()
            && detail == source_detail
            && match reason {
                PreparationRejection::Parse => *stage == crate::molecular::Stage::Parse,
                PreparationRejection::MolWriterAtomicNumberNotFound => {
                    *stage == crate::molecular::Stage::Preparation
                        && atomic_number_diagnostic(detail)
                }
            } =>
        {
            Ok(())
        }
        (
            Profile::Optimization { .. } | Profile::ConformerOptimization { .. },
            Some(GeometryPreparation::TimedOut(preparation)),
            Observation::TimedOut(observation),
        ) if preparation.limit_seconds == 60 && preparation == observation => Ok(()),
        (Profile::Coverage { .. }, None, Observation::Coverage(_)) => Ok(()),
        (Profile::Coverage { .. }, None, Observation::ParseRejected { detail })
            if !detail.is_empty() =>
        {
            Ok(())
        }
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
            // Pinned 0e0d85f4ca34aeae15dfc0f7cf5503bdb0a8e985:
            // Builder.cpp:323-326 checks both endpoints but passes params[atomIdx].
            // AngleBend.cpp:78: PRECONDITION(at2Params, "bad params pointer");
            // Independent p1 source trace identifies this original center as 1.
            // This is the approved original reference error, not a molecule patch.
            let lines: Vec<_> = detail.lines().map(str::trim).collect();
            if lines.len() == 6
                && lines[0] == "RuntimeError: Pre-condition Violation"
                && lines[1] == "bad params pointer"
                && lines[2]
                    == "Violation occurred on line 78 in file Code/ForceField/UFF/AngleBend.cpp"
                && lines[3] == "Failed Expression: at2Params"
                && lines[4] == "RDKIT: 2026.03.6"
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
            Observation::PreparationRejected { reason: a, .. },
            Observation::PreparationRejected { reason: b, .. },
        ) => a == b,
        (Observation::ParseRejected { .. }, Observation::ParseRejected { .. }) => true,
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
        if let Some(preparation) = &row.preparation {
            if let Some((reason, detail)) = reproduce_preparation_rejection(&row.case, preparation)?
            {
                return Ok(Observation::PreparationRejected { reason, detail });
            }
        }
        match row.profile {
            Profile::Coverage { add_hydrogens } => {
                let mol = match Molecule::from_smiles(&row.case.smiles) {
                    Ok(mol) => mol,
                    Err(error @ cosmolkit::SmilesError::Construction(_)) => {
                        return Err(error.to_string());
                    }
                    Err(error) => {
                        return Ok(Observation::ParseRejected {
                            detail: error.to_string(),
                        });
                    }
                };
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
                    GeometryPreparation::TimedOut(_) => {
                        return Err(
                            "reference preparation timed out; comparison must be skipped".into(),
                        );
                    }
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
                let mol = read_geometry(&row.case, geometry, true)?
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
                    GeometryPreparation::TimedOut(_) => {
                        return Err(
                            "reference preparation timed out; comparison must be skipped".into(),
                        );
                    }
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
                let base = read_geometry(&row.case, geometry, true)?;
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
        output: Value::Uff(result.unwrap_or_else(|detail| match reason {
            Some(ExpectedErrorReason::SourceTbpCenterParamsMissing { center_atom_index }) => {
                Observation::SourceTbpCenterParamsMissing { center_atom_index }
            }
            _ => Observation::Error {
                stage,
                detail,
                reason,
            },
        })),
    })
}

#[cfg(test)]
mod input_tests {
    use super::*;

    #[test]
    fn common_input_rejections_are_reproduced_not_copied() {
        let case = SmilesCase {
            id: "invalid".into(),
            smiles: "C==C".into(),
        };
        let preparation = GeometryPreparation::Rejected {
            stage: crate::molecular::Stage::Parse,
            detail: "native parse rejected".into(),
        };
        assert_eq!(
            reproduce_preparation_rejection(&case, &preparation)
                .unwrap()
                .unwrap()
                .0,
            PreparationRejection::Parse
        );
        assert!(
            reproduce_preparation_rejection(
                &SmilesCase {
                    smiles: "CCO".into(),
                    ..case
                },
                &preparation
            )
            .is_err()
        );
        let detail = [
            "RuntimeError: Pre-condition Violation",
            "Atomic number not found",
            "Violation occurred on line 159 in file Code/GraphMol/PeriodicTable.h",
            "Failed Expression: atomicNumber < byanum.size()",
            "RDKIT: 2026.03.6",
            "BOOST: 1_85",
        ]
        .join("\n");
        let preparation = GeometryPreparation::Rejected {
            stage: crate::molecular::Stage::Preparation,
            detail,
        };
        let case = SmilesCase {
            id: "writer".into(),
            smiles: format!("[C+9]{}F", "(F)".repeat(12)),
        };
        assert_eq!(
            reproduce_preparation_rejection(&case, &preparation)
                .unwrap()
                .unwrap()
                .0,
            PreparationRejection::MolWriterAtomicNumberNotFound
        );
    }

    #[test]
    fn mol_any_bond_transport_preserves_original_topology_and_quantized_xyz() {
        let case = SmilesCase {
            id: "any bond".into(),
            smiles: "C~C".into(),
        };
        let base = Molecule::from_smiles(&case.smiles)
            .unwrap()
            .with_hydrogens()
            .unwrap();
        let coordinates: Vec<_> = (0..base.num_atoms())
            .map(|i| [i as f64 * 0.125, 0.25, 0.5])
            .collect();
        let mol = Molecule::from_parts(
            base.topology().clone(),
            CoordinateBlock {
                conformers_3d: vec![Conformer3D::new(0, coordinates.clone(), true)],
                source_coordinate_dim: Some(CoordinateDimension::ThreeD),
                ..Default::default()
            },
            base.properties().clone(),
        )
        .unwrap();
        let block = mol
            .to_mol_with_params(&cosmolkit::MolBlockWriteParams {
                format: cosmolkit::SdfFormat::V3000,
                ..Default::default()
            })
            .unwrap();
        let geometry = Geometry {
            molblock: block,
            atom_count: mol.num_atoms(),
            coordinate_rows: vec![],
        };
        let read = read_geometry(&case, &geometry, true)
            .unwrap()
            .with_assigned_valence()
            .unwrap();
        assert_eq!(read.coordinates_3d(0).unwrap(), coordinates);
        assert_eq!(read.bonds()[0].order(), base.bonds()[0].order());
        let outcome = run(&Input::Uff(UffInput {
            case,
            profile: profiles(Operation::UffOptimization)[0],
            preparation: Some(GeometryPreparation::Ready(geometry)),
        }))
        .unwrap();
        assert!(
            matches!(
                outcome.output,
                Value::Uff(Observation::SourceTbpCenterParamsMissing {
                    center_atom_index: 0
                })
            ),
            "{outcome:?}"
        );
    }
}
