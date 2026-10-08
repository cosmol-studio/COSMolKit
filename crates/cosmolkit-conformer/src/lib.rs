//! Detached conformer-generation boundaries.

mod bounds;
mod smoothing;

use cosmolkit_model::{Conformer3D, CoordinateBlock, TopologyBlock};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum ConformerError {
    #[error("independent conformer capability is unsupported")]
    Unsupported,
    #[error("Cannot normalize a zero length vector")]
    CannotNormalizeZeroLengthVector,
    #[error("distance bounds matrix generation failed: {0}")]
    GenerationFailed(String),
    #[error("invalid embed parameters JSON: {0}")]
    InvalidEmbedParametersJson(String),
}

mod params;
pub use params::{EmbedFailureCause, EmbedParams};

mod numeric;

mod graph_bounds;

mod generation;

mod pruning;
#[doc(hidden)]
pub use pruning::{PruningError, PruningMatches, PruningQueryProjection, pruning_self_matches};

#[doc(hidden)]
pub use generation::{
    EmbeddingDiagnostic, GeneratedConformers, GenerationError, GenerationFailure,
    PreparationDiagnostic, generate_conformers,
};

#[doc(hidden)]
pub use bounds::{BoundKind, BoundsMatrixError, MatrixAxis};
#[doc(hidden)]
pub use graph_bounds::GraphBoundsError;
#[doc(hidden)]
pub use numeric::EmbedSeedError;

mod public_bounds;
#[doc(hidden)]
pub use public_bounds::{BoundsQueryError, dg_bounds_matrix};
