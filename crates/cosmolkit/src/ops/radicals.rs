//! Thin radical-assignment operation projection over the detached core owner.

use cosmolkit_macros::mol_op_body;

use super::OperationError;
use crate::{DerivedState, PreservationProof};

#[mol_op_body(with_assigned_radicals, parts)]
pub(crate) fn assign_radicals_impl() -> Result<(), OperationError> {
    let mut topology = parts.checkout_topology()?;
    let assignment = match cosmolkit_core::assign_radicals(&topology) {
        Ok(assignment) => assignment,
        Err(error) => {
            parts.install_topology(topology)?;
            return Err(OperationError::Radical(error));
        }
    };
    let expected = topology.atoms.len();
    let actual = assignment.radical_electrons.len();
    if actual != expected {
        parts.install_topology(topology)?;
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "with_assigned_radicals",
            field: "radical-electron",
            actual,
            expected,
        });
    }
    for (atom, radical_electrons) in topology.atoms.iter_mut().zip(assignment.radical_electrons) {
        atom.set_radical_electrons(radical_electrons);
    }
    parts.install_topology(topology)?;
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
        PreservationProof::RadicalElectronAssignment,
    )?;
    parts.apply_cip_policy()
}
