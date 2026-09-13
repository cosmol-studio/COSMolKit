//! Thin aromaticity-assignment projection over the detached core owner.

use cosmolkit_macros::mol_op_body;

use super::OperationError;
use crate::{AromaticityParams, DerivedState, PreservationProof, TopologyEditKind};

#[mol_op_body(with_assigned_aromaticity, parts)]
pub(crate) fn assign_aromaticity_impl(params: &AromaticityParams) -> Result<(), OperationError> {
    let topology = parts.checkout_topology()?;
    let rings = match cosmolkit_core::symmetrized_sssr(
        &topology,
        &cosmolkit_core::RingSearchParams::default(),
    ) {
        Ok(rings) => rings,
        Err(error) => {
            parts.install_topology(topology)?;
            return Err(OperationError::Rings(error));
        }
    };
    let assignment = match cosmolkit_core::assign_aromaticity(&topology, &rings, params) {
        Ok(assignment) => assignment,
        Err(error) => {
            parts.install_topology(topology)?;
            return Err(OperationError::Aromaticity(error));
        }
    };

    if assignment.topology.atoms.len() != topology.atoms.len() {
        let actual = assignment.topology.atoms.len();
        let expected = topology.atoms.len();
        parts.install_topology(topology)?;
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "with_assigned_aromaticity",
            field: "atom",
            actual,
            expected,
        });
    }
    if assignment.topology.bonds.len() != topology.bonds.len() {
        let actual = assignment.topology.bonds.len();
        let expected = topology.bonds.len();
        parts.install_topology(topology)?;
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "with_assigned_aromaticity",
            field: "bond",
            actual,
            expected,
        });
    }
    if assignment.topology.adjacency != topology.adjacency {
        parts.install_topology(topology)?;
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "with_assigned_aromaticity",
            field: "adjacency",
            actual: 1,
            expected: 0,
        });
    }
    if let Err(error) = assignment.topology.validate() {
        parts.install_topology(topology)?;
        return Err(OperationError::InvalidTopology(error));
    }

    parts.install_topology(assignment.topology)?;
    parts.record_topology_edit(TopologyEditKind::Local)?;
    parts.clear_cache(
        DerivedState::VALENCE
            .union(DerivedState::STEREO)
            .union(DerivedState::DRAWING)
            .union(DerivedState::FINGERPRINT),
    )?;
    parts.mark_cache_updated(DerivedState::AROMATICITY)?;
    parts.prove_preserved(
        DerivedState::RINGS
            .union(DerivedState::RING_FAMILIES)
            .union(DerivedState::COORDINATES),
        PreservationProof::AromaticityAssignment,
    )?;
    parts.apply_cip_policy()
}
