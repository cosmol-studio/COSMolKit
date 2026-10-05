//! Thin guarded atom-code assignment transport, using the sole fingerprint owner.
use crate::{AtomId, DerivedState, Molecule, OperationError, PreservationProof, TopologyEditKind};
use cosmolkit_fingerprints::{AtomCodeInput, AtomCodeOptions, atom_code};
use cosmolkit_macros::{MoleculeResult, mol_op_body};

/// A code and the explicitly returned molecule including source CIP effects.
#[derive(Clone, Debug, PartialEq, MoleculeResult)]
pub struct AtomPairAtomCodeResult<M = Molecule> {
    pub code: u32,
    #[pending_molecule]
    pub molecule: M,
}

#[mol_op_body(with_atom_pair_atom_code, parts)]
pub(crate) fn with_atom_pair_atom_code_impl(
    atom_id: AtomId,
    branch_subtract: u32,
    include_chirality: bool,
    use_legacy_stereo_perception: bool,
) -> Result<
    AtomPairAtomCodeResult<crate::PendingMolecule<super::WithAtomPairAtomCodeAccess>>,
    OperationError,
> {
    // RDKit❗❌:   python::def("GetAtomPairAtomCode", RDKit::AtomPairs::getAtomCode,
    // RDKit❗❌:               (python::arg("atom"), python::arg("branchSubtract") = 0,
    // RDKit❗❌:                python::arg("includeChirality") = false),
    // RDKit❗❌:               docString.c_str());
    // Complexity: each from_cow callback validates the full topology in O(A+B).
    // The source wrapper and no-CIP atom-code path use O(degree) local access.
    // This additional public-call validation remains and is not cost-equivalent.
    // The owner runs numPi, guarded modern CIP and property conversion in their
    // original order. Scoped COW passes existing cached valence, never derives it.
    let (code, changed) = parts.stage_topology_properties(|topology, properties, cache| {
        let valence = cache
            .valence_assignment()
            .and_then(|values| values.explicit_valence.get(atom_id.index()))
            .map(|value| *value as i8);
        let input = AtomCodeInput::from_cow(topology, properties)
            .map_err(OperationError::InvalidTopology)?;
        atom_code(
            input,
            Some(atom_id),
            valence,
            &AtomCodeOptions {
                branch_subtract,
                include_chirality,
                use_legacy_stereo_perception,
            },
        )
        .map(|assignment| assignment.into_optional_owned_parts())
        .map_err(OperationError::AtomCode)
    })?;
    // The declared weak identity discipline requires its trace even on a no-op.
    // Recording it does not detach or replace any shared block.
    parts.record_topology_edit(TopologyEditKind::Local)?;
    if changed {
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
        parts.apply_cip_policy()?;
    }
    Ok(AtomPairAtomCodeResult {
        code,
        molecule: parts.pending_molecule()?,
    })
}
