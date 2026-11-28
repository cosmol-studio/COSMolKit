//! Thin coordinate operations; local value changes and chemistry keep their owners.
use crate::{DerivedState, OperationError, PreservationProof};
use cosmolkit_macros::mol_op_body;

#[mol_op_body(with_3d_conformer, parts)]
pub(crate) fn with_3d_conformer_impl(params: &crate::EmbedParams) -> Result<(), OperationError> {
    let mut run_params = params.clone();
    let mut coordinates = parts.checkout_coordinates()?;
    let generated = cosmolkit_conformer::generate_conformers(
        parts.topology()?,
        &coordinates,
        parts.properties()?,
        1,
        &mut run_params,
        crate::conformer_projection::concrete_pruning_query,
    )
    .map_err(crate::ConformerRunError::generation);
    if let Ok(generated) = &generated {
        emit_diagnostics(&generated.diagnostics);
    }
    let metadata = generated
        .map_err(OperationError::from)
        .and_then(|generated| {
            install_embedding_coordinates(
                &mut coordinates,
                generated.clear_existing,
                generated.conformers,
                params.optimize_conf_with_id,
            )?;
            Ok((generated.conf_ids, generated.termination))
        });
    parts.install_coordinates(coordinates)?;
    let (conf_ids, termination) = metadata?;
    parts.clear_cache(DerivedState::STEREO.union(DerivedState::DRAWING))?;
    parts.prove_preserved(
        DerivedState::RINGS
            .union(DerivedState::RING_FAMILIES)
            .union(DerivedState::VALENCE)
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::FINGERPRINT),
        PreservationProof::CoordinateOnly,
    )?;
    parts.apply_cip_policy()?;
    let _ = (conf_ids, termination);
    Ok(())
}

#[mol_op_body(with_3d_conformer_result, parts)]
pub(crate) fn with_3d_conformer_result_impl(
    params: &crate::EmbedParams,
) -> Result<crate::conformer::EmbedConformerReport, OperationError> {
    let mut run_params = params.clone();
    let mut coordinates = parts.checkout_coordinates()?;
    let generated = cosmolkit_conformer::generate_conformers(
        parts.topology()?,
        &coordinates,
        parts.properties()?,
        1,
        &mut run_params,
        crate::conformer_projection::concrete_pruning_query,
    )
    .map_err(crate::ConformerRunError::generation);
    if let Ok(generated) = &generated {
        emit_diagnostics(&generated.diagnostics);
    }
    let metadata = generated
        .map_err(OperationError::from)
        .and_then(|generated| {
            install_embedding_coordinates(
                &mut coordinates,
                generated.clear_existing,
                generated.conformers,
                params.optimize_conf_with_id,
            )?;
            Ok((generated.conf_ids, generated.termination))
        });
    parts.install_coordinates(coordinates)?;
    let (conf_ids, termination) = metadata?;
    parts.clear_cache(DerivedState::STEREO.union(DerivedState::DRAWING))?;
    parts.prove_preserved(
        DerivedState::RINGS
            .union(DerivedState::RING_FAMILIES)
            .union(DerivedState::VALENCE)
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::FINGERPRINT),
        PreservationProof::CoordinateOnly,
    )?;
    parts.apply_cip_policy()?;
    Ok(crate::conformer::EmbedConformerReport {
        conf_ids,
        termination,
        params: run_params,
        requested_num_confs: 1,
    })
}

#[mol_op_body(with_3d_conformers, parts)]
pub(crate) fn with_3d_conformers_impl(
    num_confs: u32,
    params: &crate::EmbedParams,
) -> Result<(), OperationError> {
    let mut run_params = params.clone();
    let mut coordinates = parts.checkout_coordinates()?;
    let generated = cosmolkit_conformer::generate_conformers(
        parts.topology()?,
        &coordinates,
        parts.properties()?,
        num_confs,
        &mut run_params,
        crate::conformer_projection::concrete_pruning_query,
    )
    .map_err(crate::ConformerRunError::generation);
    if let Ok(generated) = &generated {
        emit_diagnostics(&generated.diagnostics);
    }
    let metadata = generated
        .map_err(OperationError::from)
        .and_then(|generated| {
            install_embedding_coordinates(
                &mut coordinates,
                generated.clear_existing,
                generated.conformers,
                params.optimize_conf_with_id,
            )?;
            Ok((generated.conf_ids, generated.termination))
        });
    parts.install_coordinates(coordinates)?;
    let (conf_ids, termination) = metadata?;
    parts.clear_cache(DerivedState::STEREO.union(DerivedState::DRAWING))?;
    parts.prove_preserved(
        DerivedState::RINGS
            .union(DerivedState::RING_FAMILIES)
            .union(DerivedState::VALENCE)
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::FINGERPRINT),
        PreservationProof::CoordinateOnly,
    )?;
    parts.apply_cip_policy()?;
    let _ = (conf_ids, termination);
    Ok(())
}

#[mol_op_body(with_3d_conformers_result, parts)]
pub(crate) fn with_3d_conformers_result_impl(
    num_confs: u32,
    params: &crate::EmbedParams,
) -> Result<crate::conformer::EmbedConformerReport, OperationError> {
    let mut run_params = params.clone();
    let mut coordinates = parts.checkout_coordinates()?;
    let generated = cosmolkit_conformer::generate_conformers(
        parts.topology()?,
        &coordinates,
        parts.properties()?,
        num_confs,
        &mut run_params,
        crate::conformer_projection::concrete_pruning_query,
    )
    .map_err(crate::ConformerRunError::generation);
    if let Ok(generated) = &generated {
        emit_diagnostics(&generated.diagnostics);
    }
    let metadata = generated
        .map_err(OperationError::from)
        .and_then(|generated| {
            install_embedding_coordinates(
                &mut coordinates,
                generated.clear_existing,
                generated.conformers,
                params.optimize_conf_with_id,
            )?;
            Ok((generated.conf_ids, generated.termination))
        });
    parts.install_coordinates(coordinates)?;
    let (conf_ids, termination) = metadata?;
    parts.clear_cache(DerivedState::STEREO.union(DerivedState::DRAWING))?;
    parts.prove_preserved(
        DerivedState::RINGS
            .union(DerivedState::RING_FAMILIES)
            .union(DerivedState::VALENCE)
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::FINGERPRINT),
        PreservationProof::CoordinateOnly,
    )?;
    parts.apply_cip_policy()?;
    Ok(crate::conformer::EmbedConformerReport {
        conf_ids,
        termination,
        params: run_params,
        requested_num_confs: num_confs,
    })
}

fn install_embedding_coordinates(
    coordinates: &mut cosmolkit_model::CoordinateBlock,
    clear_existing: bool,
    conformers: Vec<cosmolkit_model::Conformer3D>,
    refined_id: Option<u32>,
) -> Result<(), OperationError> {
    if let Some(id) = refined_id {
        // Source copies refined coordinates only after successful completion;
        // retain the selected row's identity/metadata and every other row.
        for refined in conformers {
            let destination = coordinates
                .conformers_3d
                .iter_mut()
                .find(|row| row.id() == id as usize)
                .ok_or_else(|| {
                    OperationError::InvalidCoordinates(
                        cosmolkit_model::CoordinateValidationError::Missing3DConformer {
                            id: id as usize,
                        },
                    )
                })?;
            destination
                .coordinates_mut()
                .copy_from_slice(refined.coordinates());
        }
        Ok(())
    } else {
        coordinates
            .install_generated_3d(clear_existing, conformers)
            .map_err(OperationError::InvalidCoordinates)
    }
}

fn emit_diagnostics(diagnostics: &[cosmolkit_conformer::EmbeddingDiagnostic]) {
    use cosmolkit_conformer::EmbeddingDiagnostic as D;
    for diagnostic in diagnostics {
        match diagnostic {
            D::MissingExplicitHydrogens => {
                eprintln!("Molecule does not have explicit Hs. Consider calling AddHs()")
            }
            D::MultipleFragmentsCoordinateMap => eprintln!(
                "Constrained conformer generation (via the coordMap argument) does not work with molecules that have multiple fragments."
            ),
            D::MultipleFragmentsBoundsMatrix => eprintln!(
                "Conformer generation using a user-provided boundsMat does not work with molecules that have multiple fragments. The boundsMat will be ignored."
            ),
            D::Preparation(diagnostic) => eprintln!("{diagnostic}"),
            D::HydrogenRemoval(diagnostic) => eprintln!("{diagnostic:?}"),
            D::RefinementMultipleFragments => eprintln!(
                "Seeded conformer generation does not work with molecules that have multiple fragments. Aborting."
            ),
            D::RefinementCoordinateMap => eprintln!(
                "optimizeConfWithId and coordMap are provided; initial coordinates use the random-coordinate constraint mode"
            ),
            D::RefinementFailed => {
                eprintln!("Failed to refine the conformer with the provided id, cancelling")
            }
            D::Interrupted => eprintln!("Interrupted, cancelling conformer generation"),
        }
    }
}
