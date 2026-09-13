//! Thin valence operation projection over the detached core owner.

use cosmolkit_macros::mol_op_body;

use super::OperationError;
use crate::{DerivedState, Molecule, PreservationProof, ValenceError, ValenceParams};

#[mol_op_body(with_assigned_valence, parts)]
pub(crate) fn assign_valence_impl(params: &ValenceParams) -> Result<(), OperationError> {
    let assignment = cosmolkit_core::assign_valence(parts.topology()?, params)
        .map_err(OperationError::Valence)?;
    let mut cache = parts.checkout_derived_cache()?;
    cache.install_valence_assignment(assignment);
    parts.install_derived_cache(cache)?;
    parts.mark_cache_updated(DerivedState::VALENCE)?;
    parts.prove_preserved(
        DerivedState::RINGS.union(DerivedState::STEREO),
        PreservationProof::UnchangedInput,
    )?;
    parts.apply_cip_policy()
}

impl Molecule {
    /// Reports whether one atom violates the detached owner's valence rules.
    pub fn has_valence_violation(&self, atom_id: crate::AtomId) -> Result<bool, ValenceError> {
        cosmolkit_core::has_valence_violation(self.topology(), atom_id)
    }
}
