//! Thin modern-CIP-label assignment projection over the detached stereo owner.

use cosmolkit_macros::mol_op_body;

use super::OperationError;
use crate::{CipLabelOptions, DerivedState, PreservationProof, TopologyEditKind};

#[mol_op_body(with_cip_labels, parts)]
pub(crate) fn assign_cip_labels_impl(options: &CipLabelOptions) -> Result<(), OperationError> {
    let topology = parts.checkout_topology()?;
    let properties = parts.checkout_properties()?;
    let assignment = cosmolkit_stereo::assign_cip_labels(topology, properties, options)
        .map_err(OperationError::CipLabeler)?;
    let (topology, properties) = assignment.into_parts();

    topology
        .validate()
        .map_err(OperationError::InvalidTopology)?;
    parts.install_topology(topology)?;
    parts.install_properties(properties)?;
    parts.record_topology_edit(TopologyEditKind::Local)?;
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
        PreservationProof::CipLabelAssignment,
    )?;
    parts.apply_cip_policy()
}
