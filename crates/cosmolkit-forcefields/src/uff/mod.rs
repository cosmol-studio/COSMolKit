mod angle;
pub(crate) mod api;
mod atom_typer;
pub(super) mod bond;
pub(crate) mod builder;
pub(super) mod convenience;
mod inversion;
mod inversions;
mod nonbonded;
mod optimization;
pub(crate) mod params;
mod public;
mod torsion;
mod utils;

pub use api::{UffParameterError, UffParameterErrorKind, uff_has_all_molecule_params};
pub use atom_typer::{UffTypingError, UffTypingInput};
pub use bond::BondMathError as UffBondMathError;
pub use builder::{
    AutomaticForceFieldConstructionError, ForceFieldConstructionError, NonbondedAssemblyError,
    NonbondedFragmentMappingError,
};
pub use builder::{
    DefaultTorsionQueryError as UffDefaultTorsionQueryError,
    PreparedValenceField as UffPreparedValenceField,
    TorsionBondQueryError as UffTorsionBondQueryError, UffBuilderError,
};
#[cfg(not(target_family = "wasm"))]
pub use convenience::DispatchedUffOptimizationError;
pub use convenience::{
    OptimizationStageError, SerialConformerOptimizationError, SerialCoordinateCountError,
    SerialResultCapacityError, SerialUffOptimizationError, SingleConformerOptimizationError,
    UffThreadCountError,
};
pub use inversion::{
    InversionContributionError as UffInversionContributionError,
    InversionIndexArgument as UffInversionIndexArgument,
};
pub use optimization::UffPreparedOptimizationError;
pub use params::UffParamError;
pub use public::{
    UffConformerError, UffConformerOptions, UffConformerOutcome, UffSingleError, UffSingleOptions,
    UffSingleOutcome, optimize_uff_conformers_prepared, optimize_uff_single_prepared,
};

mod evaluation;

pub use evaluation::{
    Failure as UffEvaluationFailure, UffEnergyGradient, UffEvaluationError, UffEvaluationParams,
    evaluate_uff,
};

pub(crate) use inversion::InversionContributionError;
pub(crate) use inversions::InversionContribs;

pub use api::{UffBoundsCause, UffBoundsError, uff_bond_rest_lengths};

pub use api::{MissingExplicitHydrogensError, needs_explicit_hydrogens};
