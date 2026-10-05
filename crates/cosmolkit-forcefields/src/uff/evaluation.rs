//! Read-only UFF evaluation through the existing source constructor/kernel.
use cosmolkit_core::{RingInfo, ValenceAssignment};
use cosmolkit_model::{CoordinateBlock, MoleculeProperties, TopologyBlock};
use std::{error::Error, fmt};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UffEvaluationParams {
    pub vdw_threshold: f64,
    pub conformer_id: Option<usize>,
    pub ignore_interfragment_interactions: bool,
}
impl Default for UffEvaluationParams {
    fn default() -> Self {
        // RDKit✔️✔️:     ROMol &mol, double vdwThresh = 10.0, int confId = -1,
        // RDKit✔️✔️:     bool ignoreInterfragInteractions = true) {
        // rdForceFields.cpp UFFGetMoleculeForceField source defaults; scalar O(1).
        Self {
            vdw_threshold: 10.0,
            conformer_id: None,
            ignore_interfragment_interactions: true,
        }
    }
}
#[derive(Clone, Debug, PartialEq)]
pub struct UffEnergyGradient {
    pub energy: f64,
    pub gradient: Vec<f64>,
}
impl UffEnergyGradient {
    pub fn energy(&self) -> f64 {
        self.energy
    }
    pub fn gradient(&self) -> &[f64] {
        &self.gradient
    }
}
#[derive(Debug)]
enum Failure {
    MissingConformer { requested: Option<usize> },
    Preparation(super::api::UffParameterError),
    Construction(super::builder::AutomaticForceFieldConstructionError),
    Kernel(crate::kernel::ForceFieldKernelError),
}
#[derive(Debug)]
pub struct UffEvaluationError(Failure);
impl fmt::Display for UffEvaluationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.0 {
            Failure::MissingConformer { requested } => {
                write!(f, "no stored 3D conformer for selector {requested:?}")
            }
            Failure::Preparation(e) => fmt::Display::fmt(e, f),
            Failure::Construction(e) => fmt::Display::fmt(e, f),
            Failure::Kernel(e) => fmt::Display::fmt(e, f),
        }
    }
}
impl Error for UffEvaluationError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match &self.0 {
            Failure::MissingConformer { .. } => None,
            Failure::Preparation(e) => Some(e),
            Failure::Construction(e) => Some(e),
            Failure::Kernel(e) => Some(e),
        }
    }
}
#[allow(clippy::too_many_arguments)]
pub fn evaluate_uff(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    properties: &MoleculeProperties,
    params: &UffEvaluationParams,
) -> Result<UffEnergyGradient, UffEvaluationError> {
    // RDKit❗❌: ForceFields::PyForceField *UFFGetMoleculeForceField(
    // RDKit❗❌:     ROMol &mol, double vdwThresh = 10.0, int confId = -1,
    // RDKit❗❌:     bool ignoreInterfragInteractions = true) {
    // RDKit❗❌:   ForceFields::ForceField *ff = UFF::constructForceField(
    // RDKit❗❌:       mol, vdwThresh, confId, ignoreInterfragInteractions);
    // RDKit❗❌:   auto *res = new ForceFields::PyForceField(ff);
    // RDKit❗❌:   res->initialize();
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    // Pinned rdForceFields.cpp: source construction and initialization reused.
    // Behavior: actual kernel energy and gradient owners preserve formulas,
    // traversal and first errors. No optimization or runtime mutation occurs.
    // Complexity: only selected O(A+P) conformer cloned; all other conformers,
    // topology/cache/properties borrowed. One existing builder and O(3A) result.
    (|| -> Result<UffEnergyGradient, Failure> {
        let prepared =
            super::api::prepare_parameter_query(topology, valence).map_err(Failure::Preparation)?;
        let rows = &coordinates.conformers_3d;
        let index = match params.conformer_id {
            None if !rows.is_empty() => 0,
            Some(id) => {
                rows.iter()
                    .position(|r| r.id() == id)
                    .ok_or(Failure::MissingConformer {
                        requested: Some(id),
                    })?
            }
            None => return Err(Failure::MissingConformer { requested: None }),
        };
        let row = &rows[index];
        let mut selected = row.clone();
        let context = super::builder::UffConformerContext {
            two_d: &coordinates.conformers_2d,
            before: &rows[..index],
            selected_id: row.id(),
            selected_is_3d: row.is_3d(),
            selected_props: row.props(),
            after: &rows[index + 1..],
            source_dimension: coordinates.source_coordinate_dim,
        };
        let mut diagnostics = Vec::new();
        let mut field = super::builder::construct_force_field_with_automatic_typing_from_selected(
            topology,
            &mut selected,
            &context,
            prepared.typing_state,
            rings,
            valence,
            properties,
            &mut diagnostics,
            params.vdw_threshold,
            params.ignore_interfragment_interactions,
        )
        .map_err(Failure::Construction)?;
        field.initialize().map_err(Failure::Kernel)?;
        let energy = field.calc_energy_current(None).map_err(Failure::Kernel)?;
        let mut gradient = vec![0.0; 3 * topology.atoms.len()];
        field
            .calc_grad_current(&mut gradient)
            .map_err(Failure::Kernel)?;
        Ok(UffEnergyGradient { energy, gradient })
    })()
    .map_err(UffEvaluationError)
}
