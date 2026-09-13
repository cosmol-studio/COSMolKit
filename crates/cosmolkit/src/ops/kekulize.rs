//! Thin Kekule-bond assignment projection over the detached core owner.

use cosmolkit_macros::mol_op_body;

use super::OperationError;
use crate::{DerivedState, KekulizeParams, PreservationProof, TopologyEditKind};

#[mol_op_body(with_kekulized_bonds, parts)]
pub(crate) fn kekulize_bonds_impl(params: &KekulizeParams) -> Result<(), OperationError> {
    let topology = parts.checkout_topology()?;
    let assignment = match cosmolkit_core::kekulize(&topology, params) {
        Ok(assignment) => assignment,
        Err(error) => {
            parts.install_topology(topology)?;
            return Err(OperationError::Kekulize(error));
        }
    };

    if assignment.topology.atoms.len() != topology.atoms.len() {
        let actual = assignment.topology.atoms.len();
        let expected = topology.atoms.len();
        parts.install_topology(topology)?;
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "with_kekulized_bonds",
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
            operation: "with_kekulized_bonds",
            field: "bond",
            actual,
            expected,
        });
    }
    if assignment.topology.adjacency != topology.adjacency {
        parts.install_topology(topology)?;
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "with_kekulized_bonds",
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
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::STEREO)
            .union(DerivedState::DRAWING)
            .union(DerivedState::FINGERPRINT),
    )?;
    parts.prove_preserved(
        DerivedState::RINGS
            .union(DerivedState::RING_FAMILIES)
            .union(DerivedState::COORDINATES),
        PreservationProof::KekulizeBondAssignment,
    )?;
    parts.apply_cip_policy()
}
