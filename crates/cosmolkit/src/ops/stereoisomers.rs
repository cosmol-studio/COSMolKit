//! Thin lazy source-enumeration operation bodies.
use crate::{OperationError, StereoisomerOptions};
use cosmolkit_macros::mol_multi_op_body;
#[mol_multi_op_body(enumerate_stereoisomers, parts)]
pub(crate) fn enumerate_stereoisomers_impl(
    options: &StereoisomerOptions,
) -> Result<(), OperationError> {
    let record = cosmolkit_smiles::SmilesRecord {
        topology: parts.topology()?.clone(),
        coordinates: parts.coordinates()?.clone(),
        properties: parts.properties()?.clone(),
    };
    let changes_coordinates = options.try_embedding;
    let stream = cosmolkit_stereo::enumerate_stereoisomers(
        &record,
        options.clone(),
        crate::conformer_projection::concrete_pruning_query,
    )?;
    parts.emit_lazy_prepared(
        stream.map(move |candidate| transport_candidate(candidate, changes_coordinates)),
    )
}
#[mol_multi_op_body(enumerate_stereoisomers_with_random_bits, parts)]
pub(crate) fn enumerate_stereoisomers_with_random_bits_impl(
    options: &StereoisomerOptions,
    callback: Box<dyn FnMut(usize) -> Result<num_bigint::BigUint, String> + Send + Sync + 'static>,
) -> Result<(), OperationError> {
    let record = cosmolkit_smiles::SmilesRecord {
        topology: parts.topology()?.clone(),
        coordinates: parts.coordinates()?.clone(),
        properties: parts.properties()?.clone(),
    };
    let changes_coordinates = options.try_embedding;
    let stream = cosmolkit_stereo::enumerate_stereoisomers_with_random_bits(
        &record,
        options.clone(),
        callback,
        crate::conformer_projection::concrete_pruning_query,
    )?;
    parts.emit_lazy_prepared(
        stream.map(move |candidate| transport_candidate(candidate, changes_coordinates)),
    )
}
fn transport_candidate(
    candidate: Result<cosmolkit_smiles::SmilesRecord, cosmolkit_stereo::EnumerationError>,
    changes_coordinates: bool,
) -> Result<
    (
        cosmolkit_model::TopologyBlock,
        Option<cosmolkit_model::CoordinateBlock>,
        cosmolkit_model::MoleculeProperties,
        cosmolkit_core::ValenceAssignment,
        cosmolkit_core::RingInfo,
    ),
    OperationError,
> {
    let record = candidate?;
    // Existing canonical owners materialize source-prepared cache facts.
    // This extra O(V+E) preparation is a known adapter cost, not a second algorithm.
    let valence = cosmolkit_core::assign_valence_with_options_for_topology(
        &record.topology,
        cosmolkit_core::ValenceModel::RdkitLike,
        false,
    )
    .map_err(cosmolkit_stereo::EnumerationError::from)?;
    let rings = cosmolkit_core::symmetrized_sssr(&record.topology, &Default::default())
        .map_err(cosmolkit_stereo::EnumerationError::from)?;
    Ok((
        record.topology,
        changes_coordinates.then_some(record.coordinates),
        record.properties,
        valence,
        rings,
    ))
}
