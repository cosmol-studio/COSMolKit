//! Thin atom-position operation projection over the detached core owner.

use cosmolkit_macros::mol_op_body;

use super::OperationError;
use crate::{AtomId, AtomPositionParams, DerivedState, PreservationProof};

#[mol_op_body(with_atom_position, parts)]
pub(crate) fn with_atom_position_impl(
    atom: AtomId,
    position: [f64; 3],
    params: &AtomPositionParams,
) -> Result<(), OperationError> {
    let topology = parts.checkout_topology()?;
    let atom_count = topology.atoms.len();
    let coordinates = parts.checkout_coordinates()?;
    let updated =
        cosmolkit_core::with_atom_position(&topology, &coordinates, atom, position, params)
            .map_err(OperationError::Transform)?;
    updated
        .validate_for_atom_count(atom_count)
        .map_err(OperationError::InvalidCoordinates)?;
    parts.install_topology(topology)?;
    parts.install_coordinates(updated)?;
    parts.clear_cache(DerivedState::STEREO.union(DerivedState::DRAWING))?;
    parts.prove_preserved(
        DerivedState::RINGS
            .union(DerivedState::RING_FAMILIES)
            .union(DerivedState::VALENCE)
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::FINGERPRINT),
        PreservationProof::CoordinateOnly,
    )?;
    parts.apply_cip_policy()
}
