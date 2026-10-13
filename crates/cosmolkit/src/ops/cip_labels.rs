//! Thin modern-CIP-label assignment projection over the detached stereo owner.

use cosmolkit_macros::mol_op_body;

use super::OperationError;
use crate::{CipLabelOptions, DerivedState, PreservationProof};

#[mol_op_body(with_cip_labels, parts)]
pub(crate) fn assign_cip_labels_impl(options: &CipLabelOptions) -> Result<(), OperationError> {
    let (_completion, _staged) =
        parts.stage_topology_properties_cow(|topology, properties, _cache| {
            cosmolkit_stereo::assign_cip_labels_cow(topology, properties, options)
                .map(|(topology, properties, completion)| {
                    (completion, Some((topology, properties)))
                })
                .map_err(OperationError::CipLabeler)
        })?;
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
    #[cfg(feature = "op-contracts")]
    return parts.apply_cip_assignment_policy(
        options.atoms().unwrap_or_default(),
        options.bonds().unwrap_or_default(),
        _completion,
    );
    #[cfg(not(feature = "op-contracts"))]
    parts.apply_cip_assignment_policy()
}

impl crate::Molecule {
    /// Reports the source computed marker; its string projection uses BoolTag.
    pub fn cip_computed(&self) -> Result<bool, crate::PropertyValueError> {
        // RDKit❗✔️: mol.setProp(common_properties::_CIPComputed, true, computed);
        // The canonical property owner checks StringVector computed metadata
        // and retains the BoolTag getter error. Neither state is rendered text.
        if !self.properties().is_prop_computed("_CIPComputed")? {
            return Ok(false);
        }
        self.property("_CIPComputed")
            .map(crate::PropertyValue::as_bool)
            .transpose()
            .map(|value| value.unwrap_or(false))
    }
}
