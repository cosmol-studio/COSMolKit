//! Coordinate-only operation adapters over the sole detached alignment owner.
use crate::{
    AlignmentParameters, AlignmentResult, ConformerAlignmentParameters, ConformerAlignmentReport,
    DerivedState, Molecule, OperationError, PreservationProof,
};
use cosmolkit_macros::mol_op_body;

#[mol_op_body(with_alignment_to, parts)]
pub(crate) fn with_alignment_to_impl(
    reference: &Molecule,
    params: &AlignmentParameters,
) -> Result<AlignmentResult, OperationError> {
    let mut coordinates = parts.checkout_coordinates()?;
    let outcome = (|| {
        let cache = parts.checkout_derived_cache()?;
        let result: Result<AlignmentResult, OperationError> = (|| {
            let result = crate::alignment::alignment_transform_candidate(
                parts.topology()?,
                &coordinates,
                &cache,
                reference,
                params,
            )
            .map_err(OperationError::Alignment)?;
            cosmolkit_alignment::apply_alignment(
                &mut coordinates,
                params.probe_conformer_id,
                &result,
            )
            .map_err(OperationError::Alignment)?;
            Ok(result)
        })();
        parts.install_derived_cache(cache)?;
        result
    })();
    parts.install_coordinates(coordinates)?;
    let result = outcome?;
    parts.clear_cache(DerivedState::DRAWING)?;
    parts.prove_preserved(
        DerivedState::RINGS
            .union(DerivedState::RING_FAMILIES)
            .union(DerivedState::VALENCE)
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::STEREO)
            .union(DerivedState::COORDINATES)
            .union(DerivedState::FINGERPRINT),
        PreservationProof::CoordinateOnly,
    )?;
    parts.apply_cip_policy()?;
    Ok(result)
}

#[mol_op_body(with_aligned_conformers, parts)]
pub(crate) fn with_aligned_conformers_impl(
    params: &ConformerAlignmentParameters,
) -> Result<ConformerAlignmentReport, OperationError> {
    let atom_count = parts.topology()?.atoms.len();
    let mut coordinates = parts.checkout_coordinates()?;
    let result = cosmolkit_alignment::align_conformers(&mut coordinates, atom_count, params);
    parts.install_coordinates(coordinates)?;
    let rmsds = result.map_err(OperationError::Alignment)?;
    parts.clear_cache(DerivedState::DRAWING)?;
    parts.prove_preserved(
        DerivedState::RINGS
            .union(DerivedState::RING_FAMILIES)
            .union(DerivedState::VALENCE)
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::STEREO)
            .union(DerivedState::COORDINATES)
            .union(DerivedState::FINGERPRINT),
        PreservationProof::CoordinateOnly,
    )?;
    parts.apply_cip_policy()?;
    Ok(ConformerAlignmentReport { rmsds })
}
