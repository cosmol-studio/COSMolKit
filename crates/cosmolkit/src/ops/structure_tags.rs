//! Thin 3D structure-tag projection over the detached core owner.

use cosmolkit_macros::mol_op_body;

use super::OperationError;
use crate::{DerivedState, PreservationProof, StructureTagParams};

#[mol_op_body(with_chiral_tags_from_structure, parts)]
pub(crate) fn assign_chiral_tags_from_structure_impl(
    params: &StructureTagParams,
) -> Result<(), OperationError> {
    let topology = parts.checkout_topology()?;
    let coordinates = parts.coordinates()?;
    let valence = cosmolkit_core::assign_valence(
        &topology,
        &cosmolkit_core::ValenceParams {
            model: cosmolkit_core::ValenceModel::RdkitLike,
            strict: false,
        },
    )
    .map_err(OperationError::Valence)?;
    let assignment =
        cosmolkit_core::assign_chiral_tags_from_structure(&topology, coordinates, &valence, params)
            .map_err(OperationError::Stereo)?;

    let expected_conformer = if coordinates.conformers_3d.is_empty() {
        None
    } else if params.conformer_id < 0 {
        Some(&coordinates.conformers_3d[0])
    } else {
        coordinates
            .conformers_3d
            .iter()
            .find(|conformer| conformer.id() == params.conformer_id as usize)
    };
    let expected_conformer_id = expected_conformer.map(|conformer| conformer.id());
    if assignment.selected_conformer_id != expected_conformer_id {
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "with_chiral_tags_from_structure_with_params",
            field: "selected-conformer-id",
            actual: assignment.selected_conformer_id.unwrap_or(usize::MAX),
            expected: expected_conformer_id.unwrap_or(usize::MAX),
        });
    }
    let expected_clear = expected_conformer.is_some_and(|conformer| conformer.is_3d());
    if assignment.clear_stereochem_done != expected_clear {
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "with_chiral_tags_from_structure_with_params",
            field: "clear-stereochem-done",
            actual: usize::from(assignment.clear_stereochem_done),
            expected: usize::from(expected_clear),
        });
    }
    if assignment.topology.atoms.len() != topology.atoms.len() {
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "with_chiral_tags_from_structure_with_params",
            field: "atom",
            actual: assignment.topology.atoms.len(),
            expected: topology.atoms.len(),
        });
    }
    if assignment.topology.bonds.len() != topology.bonds.len() {
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "with_chiral_tags_from_structure_with_params",
            field: "bond",
            actual: assignment.topology.bonds.len(),
            expected: topology.bonds.len(),
        });
    }
    assignment
        .topology
        .validate()
        .map_err(OperationError::InvalidTopology)?;
    parts.install_topology(assignment.topology)?;

    if assignment.clear_stereochem_done {
        let mut properties = parts.checkout_properties()?;
        properties.clear_prop("_StereochemDone");
        parts.install_properties(properties)?;
    }
    parts.clear_cache(
        DerivedState::STEREO
            .union(DerivedState::DRAWING)
            .union(DerivedState::FINGERPRINT),
    )?;
    parts.prove_preserved(
        DerivedState::RINGS
            .union(DerivedState::RING_FAMILIES)
            .union(DerivedState::VALENCE)
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::COORDINATES),
        PreservationProof::StructureTagAssignment {
            clear_stereochem_done: assignment.clear_stereochem_done,
        },
    )?;
    parts.apply_cip_policy()
}
