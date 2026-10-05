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
    let callback = CallbackAdapter(params);
    let result = cosmolkit_tautomer::enumerate_with_catalog(
        source,
        &params.catalog,
        params.policy,
        Some(&callback),
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
        parts.with_candidate_blocks(|topology, coordinates, properties, cache| {
            let source = TautomerRecordView {
                topology,
                coordinates,
                properties,
                valence: cache.valence_assignment(),
                rings: cache.valid_ring_info(),
            };
            let callback = CallbackAdapter(params);
            let record = if params.finalize_selected {
                cosmolkit_tautomer::finalize_canonical_candidate(source)
            } else {
                cosmolkit_tautomer::canonicalize_with_catalog(
                    source,
                    &params.catalog,
                    params.policy,
                    Some(&callback),
                    |view| params.score_view(view),
                )
            }
            .map_err(OperationError::Tautomer)?;
            let rows = (topology.atoms.len(), topology.bonds.len());
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
    parts.prove_preserved(
        DerivedState::COORDINATES,
        PreservationProof::StableAtomCoordinates,
    )?;
    Ok(())
}
