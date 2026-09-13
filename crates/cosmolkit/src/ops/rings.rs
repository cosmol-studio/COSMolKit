//! Thin ring operation projections over the detached core owner.

use cosmolkit_macros::mol_op_body;

use super::OperationError;
use crate::{DerivedState, PreservationProof};

#[mol_op_body(with_assigned_rings, parts)]
pub(crate) fn assign_rings_impl() -> Result<(), OperationError> {
    let rings =
        cosmolkit_core::fast_find_rings(parts.topology()?).map_err(OperationError::Rings)?;
    if !rings.is_initialized() {
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "with_assigned_rings",
            field: "initialized-ring-info",
            actual: 0,
            expected: 1,
        });
    }
    if rings.find_type() != cosmolkit_core::RingFindType::Fast {
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "with_assigned_rings",
            field: "fast-ring-info",
            actual: 0,
            expected: 1,
        });
    }
    if rings.are_ring_families_initialized() || rings.num_ring_families() != 0 {
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "with_assigned_rings",
            field: "unexpected-ring-family",
            actual: 1,
            expected: 0,
        });
    }
    if let Some((atoms, bonds)) = rings
        .atom_rings()
        .iter()
        .zip(rings.bond_rings())
        .find(|(atoms, bonds)| atoms.len() != bonds.len())
    {
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "with_assigned_rings",
            field: "ring-bond",
            actual: bonds.len(),
            expected: atoms.len(),
        });
    }
    let mut cache = parts.checkout_derived_cache()?;
    cache.install_ring_info(rings);
    parts.install_derived_cache(cache)?;
    parts.clear_cache(DerivedState::RING_FAMILIES)?;
    parts.mark_cache_updated(DerivedState::RINGS)?;
    parts.prove_preserved(
        DerivedState::VALENCE
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::STEREO)
            .union(DerivedState::COORDINATES)
            .union(DerivedState::DRAWING)
            .union(DerivedState::FINGERPRINT),
        PreservationProof::UnchangedInput,
    )?;
    parts.apply_cip_policy()
}

#[mol_op_body(with_assigned_ring_families, parts)]
pub(crate) fn assign_ring_families_impl(
    params: &cosmolkit_core::RingSearchParams,
) -> Result<(), OperationError> {
    let families = cosmolkit_core::find_ring_families(parts.topology()?, params)
        .map_err(OperationError::Rings)?;
    if !families.is_initialized() || !families.are_ring_families_initialized() {
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "with_assigned_ring_families_with_params",
            field: "initialized-ring-families",
            actual: 0,
            expected: 1,
        });
    }
    if families.find_type() != cosmolkit_core::RingFindType::OtherOrUnknown {
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "with_assigned_ring_families_with_params",
            field: "ring-family-find-type",
            actual: 0,
            expected: 1,
        });
    }
    if !families.atom_rings().is_empty() || !families.bond_rings().is_empty() {
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "with_assigned_ring_families_with_params",
            field: "unexpected-ordinary-rings",
            actual: families.atom_rings().len() + families.bond_rings().len(),
            expected: 0,
        });
    }
    if families.atom_ring_families().len() != families.bond_ring_families().len() {
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "with_assigned_ring_families_with_params",
            field: "ring-family-rows",
            actual: families.bond_ring_families().len(),
            expected: families.atom_ring_families().len(),
        });
    }
    families
        .num_relevant_cycles()
        .map_err(OperationError::Rings)?;
    let topology = parts.topology()?;
    if let Some(atom) = families
        .atom_ring_families()
        .iter()
        .flatten()
        .find(|atom| atom.index() >= topology.atoms.len())
    {
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "with_assigned_ring_families_with_params",
            field: "ring-family-atom-id",
            actual: atom.index(),
            expected: topology.atoms.len(),
        });
    }
    if let Some(bond) = families
        .bond_ring_families()
        .iter()
        .flatten()
        .find(|bond| bond.index() >= topology.bonds.len())
    {
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "with_assigned_ring_families_with_params",
            field: "ring-family-bond-id",
            actual: bond.index(),
            expected: topology.bonds.len(),
        });
    }
    let mut cache = parts.checkout_derived_cache()?;
    cache.install_ring_family_info(families);
    parts.install_derived_cache(cache)?;
    parts.mark_cache_updated(DerivedState::RING_FAMILIES)?;
    parts.prove_preserved(
        DerivedState::RINGS
            .union(DerivedState::VALENCE)
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::STEREO)
            .union(DerivedState::COORDINATES)
            .union(DerivedState::DRAWING)
            .union(DerivedState::FINGERPRINT),
        PreservationProof::UnchangedInput,
    )?;
    parts.apply_cip_policy()
}
