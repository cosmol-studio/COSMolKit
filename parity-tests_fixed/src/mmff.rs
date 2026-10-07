//! Public MMFF corpus adapter. Coordinates are prepared once for both libraries.
use crate::uff::{ConformerObservation, Geometry, GeometryPreparation};
use crate::{
    molecular::Stage,
    registry::{Input, Operation, Record, SmilesCase, Value},
};
use cosmolkit::{
    Conformer3D, CoordinateBlock, CoordinateDimension, MmffConformerOptimizationParams,
    MmffEvaluationParams, MmffOptimizationParams, Molecule, SdfCoordinateMode, SdfReadParams,
};
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Variant {
    MMFF94,
    MMFF94s,
}
impl Variant {
    fn name(self) -> &'static str {
        match self {
            Self::MMFF94 => "MMFF94",
            Self::MMFF94s => "MMFF94s",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Profile {
    Coverage {
        add_hydrogens: bool,
    },
    Optimization {
        variant: Variant,
        max_iterations: i32,
    },
    ConformerOptimization {
        variant: Variant,
        max_iterations: i32,
        conformer_count: usize,
    },
}
impl Profile {
    pub fn task_name(self) -> &'static str {
        match self {
            Self::Coverage { .. } => "mmff_has_all_molecule_params",
            Self::Optimization { .. } => "mmff_optimize",
            Self::ConformerOptimization { .. } => "mmff_optimize_conformers",
        }
    }
}
pub fn profiles(operation: Operation) -> Vec<Profile> {
    match operation {
        Operation::MmffCoverage => [false, true]
            .map(|add_hydrogens| Profile::Coverage { add_hydrogens })
            .into(),
        Operation::MmffOptimization => [Variant::MMFF94, Variant::MMFF94s]
            .map(|variant| Profile::Optimization {
                variant,
                max_iterations: 2,
            })
            .into(),
        Operation::MmffConformerOptimization => [Variant::MMFF94, Variant::MMFF94s]
            .map(|variant| Profile::ConformerOptimization {
                variant,
                max_iterations: 2,
                conformer_count: 2,
            })
            .into(),
        _ => vec![],
    }
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MmffInput {
    pub case: SmilesCase,
    pub profile: Profile,
    pub preparation: Option<GeometryPreparation>,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Observation {
    Coverage(bool),
    Optimized {
        status: i32,
        energy_bits: Option<u64>,
        xyz_bits: Vec<[u64; 3]>,
    },
    OptimizedConformers {
        conformers: Vec<ConformerObservation>,
    },
    Error {
        stage: Stage,
        detail: String,
    },
}

fn valid_xyz(rows: &[[u64; 3]], count: usize) -> bool {
    rows.len() == count
        && rows
            .iter()
            .flatten()
            .all(|bits| f64::from_bits(*bits).is_finite())
}
pub fn validate_reference(
    recipe: &MmffInput,
    prepared: &MmffInput,
    output: &Observation,
) -> Result<(), String> {
    if recipe.case != prepared.case
        || recipe.profile != prepared.profile
        || recipe
            .preparation
            .as_ref()
            .is_some_and(|value| prepared.preparation.as_ref() != Some(value))
    {
        return Err("MMFF reference changed recipe identity".into());
    }
    match (&prepared.preparation, prepared.profile, output) {
        (None, Profile::Coverage { .. }, Observation::Coverage(_)) => return Ok(()),
        (None, Profile::Coverage { .. }, Observation::Error { detail, .. })
            if !detail.is_empty() =>
        {
            return Ok(());
        }
        (
            Some(GeometryPreparation::Rejected { stage, detail }),
            Profile::Optimization { .. } | Profile::ConformerOptimization { .. },
            Observation::Error {
                stage: actual_stage,
                detail: actual_detail,
            },
        ) if matches!(stage, Stage::Parse | Stage::Preparation)
            && stage == actual_stage
            && !detail.is_empty()
            && detail == actual_detail =>
        {
            return Ok(());
        }
        _ => (),
    }
    let Some(GeometryPreparation::Ready(geometry)) = &prepared.preparation else {
        return Err("MMFF common geometry absent".into());
    };
    if geometry.atom_count == 0 || !geometry.molblock.contains("M  END") {
        return Err("invalid MMFF common geometry".into());
    }
    let valid_geometry = match prepared.profile {
        Profile::Optimization { .. } => geometry.coordinate_rows.is_empty(),
        Profile::ConformerOptimization {
            conformer_count, ..
        } => crate::uff::valid_coordinate_rows(geometry, conformer_count),
        _ => false,
    };
    let valid_output = match (prepared.profile, output) {
        (
            Profile::Optimization { .. },
            Observation::Optimized {
                status,
                energy_bits,
                xyz_bits,
            },
        ) => {
            matches!(status, -1..=1)
                && valid_xyz(xyz_bits, geometry.atom_count)
                && match energy_bits {
                    Some(bits) => *status != -1 && f64::from_bits(*bits).is_finite(),
                    None => *status == -1,
                }
        }
        (
            Profile::ConformerOptimization {
                conformer_count, ..
            },
            Observation::OptimizedConformers { conformers },
        ) => {
            conformers.len() == conformer_count
                && conformers
                    .iter()
                    .zip(&geometry.coordinate_rows)
                    .all(|(outcome, input)| {
                        outcome.conformer_id == input.conformer_id
                            && matches!(outcome.status, -1..=1)
                            && f64::from_bits(outcome.energy_bits).is_finite()
                            && valid_xyz(&outcome.xyz_bits, geometry.atom_count)
                    })
        }
        (_, Observation::Error { detail, .. }) => !detail.is_empty(),
        _ => false,
    };
    if valid_geometry && valid_output {
        Ok(())
    } else {
        Err("invalid MMFF reference geometry/observation".into())
    }
}

// Same absolute energy/coordinate tolerance as the existing public MMFF
// numerical parity suite, crates/cosmolkit/tests/mmff_numerical_parity.rs.
fn close(a: u64, b: u64) -> bool {
    let (a, b) = (f64::from_bits(a), f64::from_bits(b));
    a.is_finite() && b.is_finite() && (a - b).abs() <= 1e-6
}
fn coordinates_match(a: &[[u64; 3]], b: &[[u64; 3]]) -> bool {
    a.len() == b.len()
        && a.iter()
            .flatten()
            .zip(b.iter().flatten())
            .all(|(&a, &b)| close(a, b))
}
pub fn matches(expected: &Observation, actual: &Observation) -> bool {
    match (expected, actual) {
        (
            Observation::Optimized {
                status: a,
                energy_bits: ae,
                xyz_bits: ax,
            },
            Observation::Optimized {
                status: b,
                energy_bits: be,
                xyz_bits: bx,
            },
        ) => {
            a == b
                && match (ae, be) {
                    (Some(a), Some(b)) => close(*a, *b),
                    (None, None) => true,
                    _ => false,
                }
                && coordinates_match(ax, bx)
        }
        (
            Observation::OptimizedConformers { conformers: a },
            Observation::OptimizedConformers { conformers: b },
        ) => {
            a.len() == b.len()
                && a.iter().zip(b).all(|(a, b)| {
                    a.conformer_id == b.conformer_id
                        && a.status == b.status
                        && close(a.energy_bits, b.energy_bits)
                        && coordinates_match(&a.xyz_bits, &b.xyz_bits)
                })
        }
        _ => expected == actual,
    }
}

fn molecule(geometry: &Geometry) -> Result<Molecule, String> {
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
        return Err("MMFF common geometry atom count changed".into());
    }
    let mol = if geometry.coordinate_rows.is_empty() {
        base
    } else {
        Molecule::from_parts(
            base.topology().clone(),
            CoordinateBlock {
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
            },
            base.properties().clone(),
        )
        .map_err(|e| e.to_string())?
    };
    mol.with_assigned_valence().map_err(|e| e.to_string())
}
pub fn run(row: &MmffInput) -> Result<Record, String> {
    let mut stage = Stage::Parse;
    let result = (|| -> Result<Observation, String> {
        if let Profile::Coverage { add_hydrogens } = row.profile {
            let base = Molecule::from_smiles(&row.case.smiles).map_err(|e| e.to_string())?;
            stage = Stage::Preparation;
            let mol = if add_hydrogens {
                base.with_hydrogens()
                    .map_err(|e| e.to_string())?
                    .with_assigned_valence()
                    .map_err(|e| e.to_string())?
            } else {
                base
            };
            stage = Stage::Operation;
            return mol
                .mmff_has_all_molecule_params()
                .map(Observation::Coverage)
                .map_err(|e| e.to_string());
        }
        stage = Stage::Preparation;
        let geometry = match row
            .preparation
            .as_ref()
            .ok_or("MMFF common geometry was not prepared")?
        {
            GeometryPreparation::Ready(geometry) => geometry,
            GeometryPreparation::Rejected {
                stage: failed_stage,
                detail,
            } => {
                stage = failed_stage.clone();
                return Err(detail.clone());
            }
        };
        let mol = molecule(geometry)?;
        stage = Stage::Operation;
        match row.profile {
            Profile::Optimization {
                variant,
                max_iterations,
            } => {
                let result = mol
                    .with_mmff_optimized_with_params(&MmffOptimizationParams {
                        mmff_variant: variant.name().into(),
                        max_iterations,
                        non_bonded_threshold: 100.,
                        ignore_interfragment_interactions: true,
                        conformer_id: None,
                    })
                    .map_err(|e| e.to_string())?;
                let energy = result
                    .molecule
                    .mmff_energy_gradient_with_params(&MmffEvaluationParams {
                        mmff_variant: variant.name().into(),
                        ..Default::default()
                    })
                    .map_err(|e| e.to_string())?;
                let conf = result
                    .molecule
                    .conformers_3d()
                    .first()
                    .ok_or("MMFF optimized conformer absent")?;
                Ok(Observation::Optimized {
                    status: result.needs_more,
                    energy_bits: energy.map(|value| value.energy.to_bits()),
                    xyz_bits: conf
                        .coordinates()
                        .iter()
                        .map(|xyz| xyz.map(f64::to_bits))
                        .collect(),
                })
            }
            Profile::ConformerOptimization {
                variant,
                max_iterations,
                conformer_count,
            } => {
                let result = mol
                    .with_mmff_optimized_confs_with_params(&MmffConformerOptimizationParams {
                        mmff_variant: variant.name().into(),
                        max_iterations,
                        non_bonded_threshold: 100.,
                        ignore_interfragment_interactions: true,
                        num_threads: 1,
                    })
                    .map_err(|e| e.to_string())?;
                let coordinates = result.molecule.conformers_3d();
                if result.conformer_results.len() != conformer_count
                    || coordinates.len() != conformer_count
                {
                    return Err("MMFF conformer result count changed".into());
                }
                let conformers = result
                    .conformer_results
                    .iter()
                    .zip(coordinates)
                    .map(|(outcome, conf)| ConformerObservation {
                        conformer_id: conf.id(),
                        status: outcome.needs_more,
                        energy_bits: outcome.energy.to_bits(),
                        xyz_bits: conf
                            .coordinates()
                            .iter()
                            .map(|xyz| xyz.map(f64::to_bits))
                            .collect(),
                    })
                    .collect();
                Ok(Observation::OptimizedConformers { conformers })
            }
            Profile::Coverage { .. } => unreachable!(),
        }
    })();
    Ok(Record {
        input: Input::Mmff(row.clone()),
        output: Value::Mmff(result.unwrap_or_else(|detail| Observation::Error { stage, detail })),
    })
}
