//! Owned force-field corpus parity: common arbitrary coordinates, two iterations,
//! and binary64 bits throughout. No embedding, rounded input, or tolerance path.
use crate::registry::{Operation, Record, SmilesCase, Value};
use crate::uff::Geometry;
use cosmolkit::{
    Conformer3D, CoordinateBlock, CoordinateDimension, ForceFieldMinimizeParams,
    MmffForceFieldParams, MolecularForceField, MolecularForceFieldErrorKind, Molecule,
    SdfCoordinateMode, SdfReadParams, UffForceFieldParams,
};
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Kind {
    Mmff,
    Uff,
}
impl Kind {
    pub fn name(self) -> &'static str {
        match self {
            Self::Mmff => "mmff_force_field",
            Self::Uff => "uff_force_field",
        }
    }
    pub fn from_operation(operation: Operation) -> Option<Self> {
        match operation {
            Operation::PersistentMmff => Some(Self::Mmff),
            Operation::PersistentUff => Some(Self::Uff),
            _ => None,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Input {
    pub case: SmilesCase,
    pub kind: Kind,
    pub seed: u64,
    pub max_iterations: u32,
    pub force_tolerance_bits: u64,
    pub energy_tolerance_bits: u64,
    pub preparation: Option<Geometry>,
}
impl Input {
    pub fn new(case: SmilesCase, kind: Kind) -> Self {
        Self {
            case,
            kind,
            seed: 0x434b464620261008,
            max_iterations: 2,
            force_tolerance_bits: 1e-4_f64.to_bits(),
            energy_tolerance_bits: 1e-6_f64.to_bits(),
            preparation: None,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Snapshot {
    pub energy_bits: u64,
    pub gradient_bits: Vec<[u64; 3]>,
    pub positions_bits: Vec<[u64; 3]>,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Observation {
    /// RDKit MMFF's factory returns None for invalid MMFF properties.
    Unavailable,
    SourceTbpCenterParamsMissing {
        center_atom_index: usize,
    },
    Evaluated {
        initial: Snapshot,
        final_state: Snapshot,
        converged: bool,
    },
}

fn valid_geometry(geometry: &Geometry) -> bool {
    geometry.atom_count > 0
        && geometry.coordinate_rows.len() == 1
        && geometry.coordinate_rows[0].conformer_id == 0
        && geometry.coordinate_rows[0].xyz_bits.len() == geometry.atom_count
        && geometry.coordinate_rows[0]
            .xyz_bits
            .iter()
            .flatten()
            .all(|bits| f64::from_bits(*bits).is_finite())
}

pub fn validate_reference(
    recipe: &Input,
    prepared: &Input,
    output: &Observation,
) -> crate::Result<()> {
    if recipe.preparation.is_some()
        || prepared.case != recipe.case
        || prepared.kind != recipe.kind
        || prepared.seed != recipe.seed
        || prepared.max_iterations != recipe.max_iterations
        || prepared.force_tolerance_bits != recipe.force_tolerance_bits
        || prepared.energy_tolerance_bits != recipe.energy_tolerance_bits
    {
        return Err("owned force-field case/parameters changed".into());
    }
    let geometry = prepared
        .preparation
        .as_ref()
        .ok_or("missing force-field geometry")?;
    if !valid_geometry(geometry) {
        return Err("invalid force-field input coordinates".into());
    }
    match output {
        Observation::Unavailable if prepared.kind == Kind::Mmff => Ok(()),
        Observation::SourceTbpCenterParamsMissing { center_atom_index }
            if prepared.kind == Kind::Uff && *center_atom_index < geometry.atom_count =>
        {
            Ok(())
        }
        Observation::Evaluated {
            initial,
            final_state,
            ..
        } => {
            for state in [initial, final_state] {
                if !f64::from_bits(state.energy_bits).is_finite()
                    || state.gradient_bits.len() != geometry.atom_count
                    || state.positions_bits.len() != geometry.atom_count
                    || state
                        .gradient_bits
                        .iter()
                        .chain(&state.positions_bits)
                        .flatten()
                        .any(|bits| !f64::from_bits(*bits).is_finite())
                {
                    return Err("invalid force-field reference snapshot".into());
                }
            }
            if initial.positions_bits != geometry.coordinate_rows[0].xyz_bits {
                return Err("reference changed initial force-field coordinates".into());
            }
            Ok(())
        }
        _ => Err("unexpected force-field reference rejection".into()),
    }
}

fn snapshot(field: &MolecularForceField) -> crate::Result<Snapshot> {
    Ok(Snapshot {
        energy_bits: field.energy().map_err(|e| e.to_string())?.to_bits(),
        gradient_bits: field
            .gradient()
            .map_err(|e| e.to_string())?
            .into_iter()
            .map(|xyz| xyz.map(f64::to_bits))
            .collect(),
        positions_bits: field
            .positions()
            .into_iter()
            .map(|xyz| xyz.map(f64::to_bits))
            .collect(),
    })
}

pub fn run(input: &Input) -> crate::Result<Record> {
    let geometry = input
        .preparation
        .as_ref()
        .ok_or("force-field input not prepared")?;
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
    if base.num_atoms() != geometry.atom_count || !valid_geometry(geometry) {
        return Err("force-field geometry changed".into());
    }
    let xyz = &geometry.coordinate_rows[0].xyz_bits;
    let molecule = Molecule::from_parts(
        base.topology().clone(),
        CoordinateBlock {
            conformers_3d: vec![Conformer3D::new(
                0,
                xyz.iter().map(|v| v.map(f64::from_bits)).collect(),
                true,
            )],
            source_coordinate_dim: Some(CoordinateDimension::ThreeD),
            ..Default::default()
        },
        base.properties().clone(),
    )
    .map_err(|e| e.to_string())?
    .with_assigned_valence()
    .map_err(|e| e.to_string())?;
    let mut rejection = None;
    let field = match input.kind {
        Kind::Mmff => match molecule.mmff_force_field_with_params(&MmffForceFieldParams::new(
            Some(0),
            "MMFF94".into(),
            100.0,
            true,
        )) {
            Ok(field) => Some(field),
            Err(error) if error.kind() == MolecularForceFieldErrorKind::InvalidParameterization => {
                None
            }
            Err(error) => return Err(error.to_string()),
        },
        Kind::Uff => match molecule.uff_force_field_with_params(&UffForceFieldParams::new(
            Some(0),
            10.0,
            true,
        )) {
            Ok(field) => Some(field),
            Err(error) => {
                // As in the existing UFF adapter, only the terminal concrete
                // source cause is admissible; generic construction errors fail.
                use std::error::Error;
                let mut leaf: &(dyn Error + 'static) = &error;
                while let Some(source) = leaf.source() {
                    leaf = source;
                }
                let diagnostic = format!("{leaf:?}");
                let center = diagnostic
                    .strip_prefix("SourceTbpCenterParamsMissing { center_atom_index: ")
                    .and_then(|s| s.strip_suffix(" }"))
                    .and_then(|s| s.parse::<usize>().ok());
                match center {
                    Some(center_atom_index)
                        if error.kind() == MolecularForceFieldErrorKind::Construction
                            && center_atom_index < geometry.atom_count =>
                    {
                        rejection =
                            Some(Observation::SourceTbpCenterParamsMissing { center_atom_index });
                        None
                    }
                    _ => return Err(error.to_string()),
                }
            }
        },
    };
    let output = if let Some(mut field) = field {
        let initial = snapshot(&field)?;
        let result = field
            .minimize_with_params_(&ForceFieldMinimizeParams::new(
                input.max_iterations,
                f64::from_bits(input.force_tolerance_bits),
                f64::from_bits(input.energy_tolerance_bits),
            ))
            .map_err(|e| e.to_string())?;
        if result.iterations() > input.max_iterations {
            return Err("force field exceeded requested iteration limit".into());
        }
        let final_state = snapshot(&field)?;
        if result.energy().to_bits() != final_state.energy_bits {
            return Err("minimize outcome differs from fresh final energy".into());
        }
        Observation::Evaluated {
            initial,
            final_state,
            converged: result.converged(),
        }
    } else {
        rejection.unwrap_or(Observation::Unavailable)
    };
    let source: Vec<_> = molecule.conformers_3d()[0]
        .coordinates()
        .iter()
        .map(|v| v.map(f64::to_bits))
        .collect();
    if &source != xyz {
        return Err("owned force field mutated source molecule".into());
    }
    Ok(Record {
        input: crate::Input::PersistentForceField(input.clone()),
        output: Value::PersistentForceField(output),
    })
}
