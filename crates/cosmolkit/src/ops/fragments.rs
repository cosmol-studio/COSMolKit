//! Component chemistry belongs to core; runtime validates every mapped output.
use crate::{DerivedState, OperationError, TopologyEditKind};
use cosmolkit_macros::{mol_multi_op_body, mol_op_body};

#[mol_multi_op_body(fragments, parts)]
pub(crate) fn fragments_impl() -> Result<(), OperationError> {
    let fragments = cosmolkit_core::get_molecule_fragments(
        parts.topology()?,
        parts.coordinates()?,
        parts.properties()?,
        true,
        true,
    )
    .map_err(OperationError::Fragments)?;
    parts.emit_mapped(
        fragments
            .into_iter()
            .map(|fragment| fragment.into_mapped_parts())
            .collect(),
    )
}

#[mol_op_body(largest_fragment, parts)]
pub(crate) fn largest_fragment_impl() -> Result<(), OperationError> {
    let (mapping, prepared) =
        parts.with_mutable_candidate_blocks(|topology, coordinates, properties, _cache| {
            let fragment =
                cosmolkit_core::get_largest_molecule_fragment(topology, coordinates, properties)
                    .map_err(OperationError::Fragments)?
                    .ok_or(OperationError::EmptyFragments)?;
            let (t, c, p, mapping, prepared) = fragment.into_mapped_parts();
            *topology = t;
            *coordinates = c;
            *properties = p;
            let has_prepared = prepared.is_some();
            if let Some((valence, rings)) = prepared {
                _cache.install_valence_assignment(valence);
                _cache.install_ring_info(rings);
            }
            Ok((mapping, has_prepared))
        })?;
    parts.record_topology_edit(TopologyEditKind::Compacting)?;
    parts.record_topology_mapping(mapping)?;
    parts.clear_cache(
        DerivedState::RING_FAMILIES
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::STEREO)
            .union(DerivedState::COORDINATES)
            .union(DerivedState::DRAWING)
            .union(DerivedState::FINGERPRINT),
    )?;
    if prepared {
        parts.mark_cache_updated(DerivedState::VALENCE.union(DerivedState::RINGS))?;
    } else {
        parts.clear_cache(DerivedState::VALENCE.union(DerivedState::RINGS))?;
    }
    parts.apply_cip_policy()
}
