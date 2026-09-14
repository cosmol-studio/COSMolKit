//! Thin potential-stereochemistry projection over the detached core owner.

use cosmolkit_macros::{MoleculeResult, mol_op_body};

use super::OperationError;
use crate::{
    DerivedState, Molecule, PotentialStereoInfo, PotentialStereoParams, PreservationProof,
    RingStereoRelation,
};

/// Public result of potential-stereochemistry perception.
///
/// Analysis rows never contain a live runtime handle. A cleaned molecule is
/// present only when `PotentialStereoParams::clean` requested the source
/// cleanup branch and the complete detached result passed runtime validation.
#[derive(Clone, Debug, PartialEq, MoleculeResult)]
pub struct PotentialStereoResult<M = Molecule> {
    pub stereo: Vec<PotentialStereoInfo>,
    pub atom_ranks: Vec<u32>,
    pub ring_relations: Vec<RingStereoRelation>,
    #[pending_molecule]
    pub cleaned_molecule: Option<M>,
}

#[mol_op_body(potential_stereo, parts)]
pub(crate) fn potential_stereo_impl(
    params: &PotentialStereoParams,
) -> Result<
    PotentialStereoResult<crate::PendingMolecule<super::PotentialStereoAccess>>,
    OperationError,
> {
    let topology = parts.checkout_topology()?;
    let valence = cosmolkit_core::assign_valence(
        &topology,
        &cosmolkit_core::ValenceParams {
            model: cosmolkit_core::ValenceModel::RdkitLike,
            strict: false,
        },
    )
    .map_err(OperationError::Valence)?;
    let rings =
        cosmolkit_core::symmetrized_sssr(&topology, &cosmolkit_core::RingSearchParams::default())
            .map_err(OperationError::Rings)?;
    let mut assignment = cosmolkit_core::potential_stereo(&topology, &valence, &rings, params)
        .map_err(OperationError::PotentialStereo)?;

    if assignment.atom_ranks.len() != topology.atoms.len() {
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "potential_stereo_with_params",
            field: "atom-rank",
            actual: assignment.atom_ranks.len(),
            expected: topology.atoms.len(),
        });
    }

    let installed_topology = match (params.clean, assignment.cleaned_topology.take()) {
        (true, Some(cleaned)) => cleaned,
        (true, None) => {
            return Err(OperationError::InvalidAlgorithmResult {
                operation: "potential_stereo_with_params",
                field: "cleaned-topology",
                actual: 0,
                expected: 1,
            });
        }
        (false, None) => topology,
        (false, Some(_)) => {
            return Err(OperationError::InvalidAlgorithmResult {
                operation: "potential_stereo_with_params",
                field: "unexpected-cleaned-topology",
                actual: 1,
                expected: 0,
            });
        }
    };
    installed_topology
        .validate()
        .map_err(OperationError::InvalidTopology)?;
    parts.install_topology(installed_topology)?;
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
        PreservationProof::StereoCleanup,
    )?;
    parts.apply_cip_policy()?;

    Ok(PotentialStereoResult {
        stereo: assignment.stereo,
        atom_ranks: assignment.atom_ranks,
        ring_relations: assignment.ring_relations,
        cleaned_molecule: if params.clean {
            Some(parts.pending_molecule()?)
        } else {
            None
        },
    })
}
