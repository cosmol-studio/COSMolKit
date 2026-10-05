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
    let metadata = generated.map(|generated| {
        coordinates.install_generated_3d(generated.clear_existing, generated.conformers);
        generated.conf_ids
    });
    parts.install_coordinates(coordinates)?;
    let conf_ids = metadata?;
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
    let _ = conf_ids;
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
    let metadata = generated.map(|generated| {
        coordinates.install_generated_3d(generated.clear_existing, generated.conformers);
        generated.conf_ids
    });
    parts.install_coordinates(coordinates)?;
    let conf_ids = metadata?;
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
    let metadata = generated.map(|generated| {
        coordinates.install_generated_3d(generated.clear_existing, generated.conformers);
        generated.conf_ids
    });
    parts.install_coordinates(coordinates)?;
    let conf_ids = metadata?;
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
    let _ = conf_ids;
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
    let metadata = generated.map(|generated| {
        coordinates.install_generated_3d(generated.clear_existing, generated.conformers);
        generated.conf_ids
    });
    parts.install_coordinates(coordinates)?;
    let conf_ids = metadata?;
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
        params: run_params,
        requested_num_confs: num_confs,
    })
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
            D::Interrupted => eprintln!("Interrupted, cancelling conformer generation"),
        }
    }
}
