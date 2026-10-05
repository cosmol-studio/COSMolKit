//! Generated manual-coordinate projections; interpretation belongs to core.
use super::OperationError;
use crate::{
    Coordinate2DInputParams, Coordinate3DInputParams, DerivedState, PreservationProof,
    Replace3DCoordinatesParams,
};
use cosmolkit_macros::mol_op_body;

#[mol_op_body(with_2d_coordinate_block, parts)]
pub(crate) fn with_2d_coordinate_block_impl(
    coordinates: Vec<Vec<f64>>,
    params: &Coordinate2DInputParams,
) -> Result<(), OperationError> {
    let atom_count = parts.topology()?.atoms.len();
    let mut block = parts.checkout_coordinates()?;
    let result =
        cosmolkit_core::install_2d_coordinates(&mut block, atom_count, coordinates, params);
    // Return the complete detached block before propagating a domain failure.
    parts.install_coordinates(block)?;
    result.map_err(OperationError::CoordinateInput)?;
    parts.clear_cache(DerivedState::DRAWING)?;
    parts.prove_preserved(
        DerivedState::RINGS
            .union(DerivedState::RING_FAMILIES)
            .union(DerivedState::VALENCE)
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::STEREO)
            .union(DerivedState::FINGERPRINT),
        PreservationProof::CoordinateOnly,
    )?;
    parts.apply_cip_policy()
}
#[mol_op_body(with_3d_coordinates, parts)]
pub(crate) fn with_3d_coordinates_impl(
    coordinates: Vec<Vec<f64>>,
    params: &Replace3DCoordinatesParams,
) -> Result<(), OperationError> {
    let atom_count = parts.topology()?.atoms.len();
    let mut block = parts.checkout_coordinates()?;
    let result =
        cosmolkit_core::replace_3d_coordinates(&mut block, atom_count, coordinates, params);
    parts.install_coordinates(block)?;
    result.map_err(OperationError::CoordinateInput)?;
    parts.clear_cache(DerivedState::DRAWING)?;
    parts.prove_preserved(
        DerivedState::RINGS
            .union(DerivedState::RING_FAMILIES)
            .union(DerivedState::VALENCE)
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::STEREO)
            .union(DerivedState::FINGERPRINT),
        PreservationProof::CoordinateOnly,
    )?;
    parts.apply_cip_policy()
}
#[mol_op_body(with_added_3d_conformer, parts)]
pub(crate) fn with_added_3d_conformer_impl(
    coordinates: Vec<Vec<f64>>,
    params: &Coordinate3DInputParams,
) -> Result<usize, OperationError> {
    let atom_count = parts.topology()?.atoms.len();
    let mut block = parts.checkout_coordinates()?;
    let result = cosmolkit_core::append_3d_conformer(&mut block, atom_count, coordinates, params);
    parts.install_coordinates(block)?;
    let position = result.map_err(OperationError::CoordinateInput)?;
    parts.clear_cache(DerivedState::DRAWING)?;
    parts.prove_preserved(
        DerivedState::RINGS
            .union(DerivedState::RING_FAMILIES)
            .union(DerivedState::VALENCE)
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::STEREO)
            .union(DerivedState::FINGERPRINT),
        PreservationProof::CoordinateOnly,
    )?;
    parts.apply_cip_policy()?;
    Ok(position)
}
#[mol_op_body(with_only_3d_conformer, parts)]
pub(crate) fn with_only_3d_conformer_impl(
    coordinates: Vec<Vec<f64>>,
    params: &Coordinate3DInputParams,
) -> Result<usize, OperationError> {
    let atom_count = parts.topology()?.atoms.len();
    let mut block = parts.checkout_coordinates()?;
    let result =
        cosmolkit_core::install_only_3d_conformer(&mut block, atom_count, coordinates, params);
    parts.install_coordinates(block)?;
    let position = result.map_err(OperationError::CoordinateInput)?;
    parts.clear_cache(DerivedState::DRAWING)?;
    parts.prove_preserved(
        DerivedState::RINGS
            .union(DerivedState::RING_FAMILIES)
            .union(DerivedState::VALENCE)
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::STEREO)
            .union(DerivedState::FINGERPRINT),
        PreservationProof::CoordinateOnly,
    )?;
    parts.apply_cip_policy()?;
    Ok(position)
}
#[mol_op_body(with_cleared_3d_conformers, parts)]
pub(crate) fn with_cleared_3d_conformers_impl() -> Result<(), OperationError> {
    let mut block = parts.checkout_coordinates()?;
    cosmolkit_core::clear_3d_conformers(&mut block);
    parts.install_coordinates(block)?;
    parts.clear_cache(DerivedState::DRAWING)?;
    parts.prove_preserved(
        DerivedState::RINGS
            .union(DerivedState::RING_FAMILIES)
            .union(DerivedState::VALENCE)
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::STEREO)
            .union(DerivedState::FINGERPRINT),
        PreservationProof::CoordinateOnly,
    )?;
    parts.apply_cip_policy()
}
