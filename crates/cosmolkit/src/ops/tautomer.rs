//! Generated capabilities delegate TAU chemistry to its unique detached owner.
use crate::tautomer::{CallbackAdapter, EnumerationMetadata};
use crate::{DerivedState, OperationError, PreservationProof, TautomerParams, TopologyEditKind};
use cosmolkit_macros::{mol_multi_op_body, mol_op_body};
use cosmolkit_tautomer::TautomerRecordView;

#[mol_multi_op_body(enumerate_tautomers_with_params, parts)]
pub(crate) fn enumerate_tautomers_impl(
    params: &TautomerParams,
) -> Result<EnumerationMetadata, OperationError> {
    let cache = parts.derived_cache()?;
    let source = TautomerRecordView {
        topology: parts.topology()?,
        coordinates: parts.coordinates()?,
        properties: parts.properties()?,
        valence: cache.valence_assignment(),
        rings: cache.valid_ring_info(),
    };
    let mut source_rings = cache.valid_ring_info().cloned().unwrap_or_else(|| {
        cosmolkit_core::RingInfo::new(
            cosmolkit_core::RingFindType::OtherOrUnknown,
            source.topology.atoms.len(),
            source.topology.bonds.len(),
        )
    });
    let mut callback = CallbackAdapter { params };
    let result = cosmolkit_tautomer::enumerate_with_catalog(
        cosmolkit_tautomer::TautomerScoreView::new(source, &mut source_rings),
        &params.catalog,
        params.policy,
        Some(&mut callback),
    )
    .map_err(OperationError::Tautomer)?;
    let (keys, candidates): (Vec<_>, Vec<_>) = result
        .entries
        .into_iter()
        .map(|(key, record)| {
            (
                key,
                (
                    record.topology,
                    record.coordinates,
                    record.properties,
                    record.valence,
                    record.rings,
                ),
            )
        })
        .unzip();
    parts.emit_prepared(candidates)?;
    Ok(EnumerationMetadata {
        keys,
        status: result.status,
        modified_atoms: result.modified_atoms,
        modified_bonds: result.modified_bonds,
    })
}

#[mol_op_body(canonical_tautomer_with_params, parts)]
pub(crate) fn canonical_tautomer_impl(params: &TautomerParams) -> Result<(), OperationError> {
    let (atoms, bonds) =
        parts.with_mutable_candidate_blocks(|topology, coordinates, properties, cache| {
            let source = TautomerRecordView {
                topology,
                coordinates,
                properties,
                valence: cache.valence_assignment(),
                rings: cache.valid_ring_info(),
            };
            let mut callback = CallbackAdapter { params };
            let record = if params.finalize_selected {
                cosmolkit_tautomer::finalize_canonical_candidate(source)
            } else {
                let mut score =
                    |view: cosmolkit_tautomer::TautomerScoreView<'_>| params.score_view(view);
                let custom = params.scorer().is_some() || params.score_params.terms.is_some();
                let mut source_rings = cache.valid_ring_info().cloned().unwrap_or_else(|| {
                    cosmolkit_core::RingInfo::new(
                        cosmolkit_core::RingFindType::OtherOrUnknown,
                        source.topology.atoms.len(),
                        source.topology.bonds.len(),
                    )
                });
                cosmolkit_tautomer::canonicalize_with_catalog(
                    cosmolkit_tautomer::TautomerScoreView::new(source, &mut source_rings),
                    &params.catalog,
                    params.policy,
                    Some(&mut callback),
                    if custom { Some(&mut score) } else { None },
                )
            }
            .map_err(OperationError::Tautomer)?;
            let rows = (topology.atoms.len(), topology.bonds.len());
            if let Some(restored) = record.coordinates {
                *coordinates = restored;
            }
            *topology = record.topology;
            *properties = record.properties;
            cache.clear(DerivedState::RINGS.union(DerivedState::VALENCE));
            cache.install_valence_assignment(record.valence);
            cache.install_ring_info(record.rings);
            Ok(rows)
        })?;
    parts.record_topology_edit(TopologyEditKind::Local)?;
    parts.record_topology_mapping(cosmolkit_model::TopologyMapping::identity(atoms, bonds))?;
    parts.clear_cache(
        DerivedState::AROMATICITY
            .union(DerivedState::STEREO)
            .union(DerivedState::RING_FAMILIES)
            .union(DerivedState::DRAWING)
            .union(DerivedState::FINGERPRINT),
    )?;
    parts.mark_cache_updated(DerivedState::RINGS.union(DerivedState::VALENCE))?;
    parts.apply_cip_policy()?;
    parts.mark_cache_updated(DerivedState::COORDINATES)?;
    Ok(())
}

#[mol_op_body(with_assigned_symm_sssr, parts)]
pub(crate) fn assign_symm_sssr_impl() -> Result<(), OperationError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::TautomerScoringFunctions::scoreRings live materialization
    // RDKit❗❌:   auto ringInfo = mol.getRingInfo();
    // RDKit❗❌:   if (!ringInfo->isSymmSssr()) {
    // RDKit❗❌:     MolOps::symmetrizeSSSR(const_cast<ROMol &>(mol));
    // RDKit❗❌:     ringInfo = mol.getRingInfo();
    // RDKit❗❌:   }
    // END RDKIT CPP FUNCTION RDKit::TautomerScoringFunctions::scoreRings live materialization

    let rings = cosmolkit_core::symmetrized_sssr(parts.topology()?, &Default::default())
        .map_err(OperationError::Rings)?;
    if !rings.is_initialized() {
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "with_assigned_symm_sssr",
            field: "initialized-ring-info",
            actual: 0,
            expected: 1,
        });
    }
    if !rings.is_symm_sssr() {
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "with_assigned_symm_sssr",
            field: "symm-sssr-ring-info",
            actual: 0,
            expected: 1,
        });
    }
    if rings.are_ring_families_initialized() || rings.num_ring_families() != 0 {
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "with_assigned_symm_sssr",
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
            operation: "with_assigned_symm_sssr",
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

#[mol_op_body(with_installed_tautomer_score_cache, parts)]
pub(crate) fn install_tautomer_score_cache_impl(
    rings: &cosmolkit_core::RingInfo,
) -> Result<(), OperationError> {
    let rings = rings.clone();
    if !rings.is_initialized() {
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "with_installed_tautomer_score_cache",
            field: "initialized-ring-info",
            actual: 0,
            expected: 1,
        });
    }
    if !rings.is_symm_sssr() {
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "with_installed_tautomer_score_cache",
            field: "symm-sssr-ring-info",
            actual: 0,
            expected: 1,
        });
    }
    if rings.are_ring_families_initialized() || rings.num_ring_families() != 0 {
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "with_installed_tautomer_score_cache",
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
            operation: "with_installed_tautomer_score_cache",
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
