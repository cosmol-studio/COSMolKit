//! Thin generated scaffold transactions; all chemistry is in cosmolkit-core.
use super::OperationError;
use crate::{DerivedState, TopologyEditKind};
use cosmolkit_macros::mol_op_body;

#[mol_op_body(murcko_scaffold, parts)]
pub(crate) fn murcko_scaffold_impl() -> Result<(), OperationError> {
    let topology = parts.checkout_topology()?;
    let coordinates = parts.checkout_coordinates()?;
    let properties = parts.checkout_properties()?;
    let mut cache = parts.checkout_derived_cache()?;
    let result = cosmolkit_core::murcko_scaffold(
        topology,
        coordinates,
        properties,
        cache.valid_ring_info(),
        cache.valence_assignment(),
    )
    .map_err(OperationError::Scaffold)?;
    cache.install_valence_assignment(result.valence.expect("complete MolHash finalization"));
    if let Some(rings) = result.rings {
        cache.install_ring_info(rings);
    }
    parts.install_derived_cache(cache)?;
    parts.mark_cache_updated(DerivedState::VALENCE.union(DerivedState::RINGS))?;
    parts.install_topology(result.topology)?;
    parts.install_coordinates(result.coordinates)?;
    parts.install_properties(result.properties)?;
    parts.record_topology_edit(TopologyEditKind::Compacting)?;
    parts.record_topology_mapping(result.mapping)?;
    parts.apply_runtime_remap()?;
    parts.clear_cache(
        DerivedState::RING_FAMILIES
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::STEREO)
            .union(DerivedState::COORDINATES)
            .union(DerivedState::DRAWING)
            .union(DerivedState::FINGERPRINT),
    )?;
    parts.apply_cip_policy()
}

#[mol_op_body(net_scaffold, parts)]
pub(crate) fn net_scaffold_impl() -> Result<(), OperationError> {
    let topology = parts.checkout_topology()?;
    let coordinates = parts.checkout_coordinates()?;
    let properties = parts.checkout_properties()?;
    let mut cache = parts.checkout_derived_cache()?;
    let result = cosmolkit_core::net_scaffold(
        topology,
        coordinates,
        properties,
        cache.valid_ring_info(),
        cache.valence_assignment(),
    )
    .map_err(OperationError::Scaffold)?;
    cache.install_valence_assignment(result.valence.expect("complete MolHash finalization"));
    if let Some(rings) = result.rings {
        cache.install_ring_info(rings);
    }
    parts.install_derived_cache(cache)?;
    parts.mark_cache_updated(DerivedState::VALENCE.union(DerivedState::RINGS))?;
    parts.install_topology(result.topology)?;
    parts.install_coordinates(result.coordinates)?;
    parts.install_properties(result.properties)?;
    parts.record_topology_edit(TopologyEditKind::Compacting)?;
    parts.record_topology_mapping(result.mapping)?;
    parts.apply_runtime_remap()?;
    parts.clear_cache(
        DerivedState::RING_FAMILIES
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::STEREO)
            .union(DerivedState::COORDINATES)
            .union(DerivedState::DRAWING)
            .union(DerivedState::FINGERPRINT),
    )?;
    parts.apply_cip_policy()
}

#[mol_op_body(murcko_decompose, parts)]
pub(crate) fn murcko_decompose_impl() -> Result<(), OperationError> {
    let topology = parts.checkout_topology()?;
    let coordinates = parts.checkout_coordinates()?;
    let properties = parts.checkout_properties()?;
    let cache = parts.checkout_derived_cache()?;
    let result = cosmolkit_core::murcko_decompose(
        topology,
        coordinates,
        properties,
        cache.valid_ring_info(),
    )
    .map_err(OperationError::Scaffold)?;

    parts.install_derived_cache(cache)?;

    parts.install_topology(result.topology)?;
    parts.install_coordinates(result.coordinates)?;
    parts.install_properties(result.properties)?;
    parts.record_topology_edit(TopologyEditKind::Compacting)?;
    parts.record_topology_mapping(result.mapping)?;
    parts.apply_runtime_remap()?;
    parts.clear_cache(
        DerivedState::VALENCE
            .union(DerivedState::RINGS)
            .union(DerivedState::RING_FAMILIES)
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::STEREO)
            .union(DerivedState::COORDINATES)
            .union(DerivedState::DRAWING)
            .union(DerivedState::FINGERPRINT),
    )?;
    parts.apply_cip_policy()
}
