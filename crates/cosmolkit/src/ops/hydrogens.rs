//! Hydrogen operation bodies; chemistry remains in its final algorithm owner.

use cosmolkit_macros::mol_op_body;

use super::OperationError;
use crate::{AddHsParams, DerivedState, PreservationProof, TopologyEditKind};

#[mol_op_body(with_hydrogens, parts)]
pub(crate) fn add_hydrogens_impl(params: &AddHsParams) -> Result<(), OperationError> {
    let topology = parts.checkout_topology()?;
    let coordinates = parts.checkout_coordinates()?;
    let properties = parts.checkout_properties()?;
    let result =
        cosmolkit_core::add_hydrogens_with_params(topology, coordinates, properties, params)
            .map_err(OperationError::Hydrogen)?;

    parts.install_topology(result.topology)?;
    parts.install_coordinates(result.coordinates)?;
    parts.install_properties(result.properties)?;
    parts.record_topology_edit(TopologyEditKind::Appending)?;
    parts.record_topology_mapping(result.mapping)?;
    parts.apply_runtime_remap()?;
    parts.clear_cache(
        DerivedState::VALENCE
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::STEREO)
            .union(DerivedState::DRAWING)
            .union(DerivedState::FINGERPRINT),
    )?;
    parts.prove_preserved(
        DerivedState::RINGS.union(DerivedState::RING_FAMILIES),
        PreservationProof::LeafAtomAppend,
    )?;
    parts.apply_cip_policy()
}

#[mol_op_body(without_hydrogens, parts)]
pub(crate) fn remove_hydrogens_impl(
    params: &cosmolkit_core::RemoveHsParams,
) -> Result<(), OperationError> {
    let topology = parts.checkout_topology()?;
    let coordinates = parts.checkout_coordinates()?;
    let properties = parts.checkout_properties()?;
    let result =
        cosmolkit_core::remove_hydrogens_with_params(topology, coordinates, properties, params)
            .map_err(OperationError::Hydrogen)?;

    parts.install_topology(result.topology)?;
    parts.install_coordinates(result.coordinates)?;
    parts.install_properties(result.properties)?;
    parts.record_topology_edit(TopologyEditKind::Compacting)?;
    parts.record_topology_mapping(result.mapping)?;
    parts.apply_runtime_remap()?;

    // CK-VALENCE-001: sanitize=false allows an unsanitized molecule, not a stale
    // cache advertised as valid. Clear both the validity bit and stored value.
    // For sanitize=true, only a complete final-topology assignment may be
    // installed. None means unavailable, never a swallowed calculation error.
    // The existing operation_defined effect permits both update and clear;
    // it does not make historical RDKit cache values valid for current state.
    if !params.sanitize {
        parts.clear_cache(DerivedState::VALENCE)?;
    } else if let Some(valence) = result.final_valence {
        let mut cache = parts.checkout_derived_cache()?;
        cache.install_valence_assignment(valence);
        parts.install_derived_cache(cache)?;
        parts.mark_cache_updated(DerivedState::VALENCE)?;
    } else {
        parts.clear_cache(DerivedState::VALENCE)?;
    }
    parts.clear_cache(
        DerivedState::RINGS
            .union(DerivedState::RING_FAMILIES)
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::STEREO)
            .union(DerivedState::DRAWING)
            .union(DerivedState::FINGERPRINT),
    )?;
    parts.apply_cip_policy()
}
