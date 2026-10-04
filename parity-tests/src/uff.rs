//! Public UFF parity profiles; no force-field implementation or owner dependency.
use crate::registry::{Input, Operation, Record, SmilesCase, Value};
use cosmolkit::{
    Conformer3D, CoordinateBlock, CoordinateDimension, Molecule, SdfCoordinateMode, SdfReadParams,
    UffConformerOptimizationParams, UffOptimizationParams,
};
use serde::{Deserialize, Serialize};

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
            max_iterations: 1,
            vdw_threshold: 100,
            ignore_interfragment_interactions: true,
            conformer_id: None,
        }]
        .into(),
        Operation::UffConformerOptimization => [Profile::ConformerOptimization {
            add_hydrogens: true,
            max_iterations: 1,
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

fn valid_coordinate_rows(geometry: &Geometry, expected_count: usize) -> bool {
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

pub fn matches(expected: &Observation, actual: &Observation) -> bool {
    !matches!(expected, Observation::Error { .. })
        && !matches!(actual, Observation::Error { .. })
        && expected == actual
}

pub fn run(input: &Input) -> Result<Record, String> {
    let Input::Uff(row) = input else {
        return Err("expected UFF input".into());
    };
    let mut stage = crate::molecular::Stage::Parse;
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
                        // There is no common geometry on which either engine
                        // can run. Retain the preparation failure, never pass it.
                        stage = failed_stage.clone();
                        return Err(format!("common input rejected: {detail}"));
                    }
                };
                stage = crate::molecular::Stage::Preparation;
                let mol = Molecule::from_sdf_with_params(
                    &geometry.molblock,
                    &SdfReadParams {
                        sanitize: true,
                        remove_hydrogens: false,
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
                    .with_uff_optimized_coordinates_with_params(&UffOptimizationParams {
                        max_iterations,
                        vdw_threshold: f64::from(vdw_threshold),
                        ignore_interfragment_interactions,
                        conformer_id,
                    })
                    .map_err(|e| e.to_string())?;
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
                        return Err(format!("common input rejected: {detail}"));
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
                        remove_hydrogens: false,
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
                        max_iterations,
                        vdw_threshold: f64::from(vdw_threshold),
                        ignore_interfragment_interactions,
                    })
                    .map_err(|e| e.to_string())?;
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
        output: Value::Uff(result.unwrap_or_else(|detail| Observation::Error { stage, detail })),
    })
}
