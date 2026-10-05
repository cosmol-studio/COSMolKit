//! Source-defined bounds query preparation over unique core owners.
use cosmolkit_model::TopologyBlock;
#[derive(Debug, thiserror::Error)]
pub enum BoundsQueryError {
    #[error(transparent)]
    Topology(#[from] cosmolkit_model::TopologyValidationError),
    #[error(transparent)]
    Valence(#[from] cosmolkit_core::ValenceError),
    #[error(transparent)]
    Rings(#[from] cosmolkit_core::RingFindingError),
    #[error(transparent)]
    Conjugation(#[from] cosmolkit_core::ConjugationError),
    #[error(transparent)]
    Hybridization(#[from] cosmolkit_core::HybridizationError),
    #[error(transparent)]
    Bounds(#[from] crate::GraphBoundsError),
}
pub fn dg_bounds_matrix(topology: &TopologyBlock) -> Result<Vec<Vec<f64>>, BoundsQueryError> {
    topology.validate()?;
    let valence = cosmolkit_core::assign_valence_for_topology(
        topology,
        cosmolkit_core::ValenceModel::RdkitLike,
    )?;
    let rings = cosmolkit_core::symmetrize_sssr_with_options_from_parts(
        topology.atoms.len(),
        &topology.bonds,
        &topology.adjacency,
        false,
        false,
    )?;
    let conjugated = cosmolkit_core::assign_conjugation_flags(topology, &valence)?;
    let hybridizations =
        cosmolkit_core::assign_hybridization_with_conjugation(topology, &valence, &conjugated)?
            .values;
    let bounds = crate::graph_bounds::build_bounds_matrix(
        topology,
        &rings,
        &valence,
        &hybridizations,
        &conjugated,
        true,
        false,
        true,
        false,
    )?;
    Ok((0..topology.atoms.len())
        .map(|i| {
            (0..topology.atoms.len())
                .map(|j| unsafe { bounds.get_val_unchecked(i, j) })
                .collect()
        })
        .collect())
}
