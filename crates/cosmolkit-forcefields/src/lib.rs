//! Detached UFF/MMFF force-field boundaries.

mod geometry;
mod kernel;
mod mmff;
mod optimizer;
mod uff;

pub use uff::{
    UffConformerError, UffConformerOptions, UffConformerOutcome, UffSingleError, UffSingleOptions,
    UffSingleOutcome, optimize_uff_conformers_prepared, optimize_uff_single_prepared,
};
pub use uff::{UffParameterError, UffParameterErrorKind, uff_has_all_molecule_params};

pub use mmff::mol_properties::{MmffAtomProperties, MmffMolPropertiesError, MmffVariant};
pub use mmff::properties_api::{
    MmffProperties, MmffPropertiesParams, mmff_has_all_molecule_params, mmff_properties,
};

pub use mmff::optimization::{
    MmffConformerOptimizationParams, MmffConformerOutcomes, MmffOptimizationError,
    MmffOptimizationParams, MmffOptimizeMoleculeConfResult, MmffSingleOutcome,
    optimize_mmff_conformers, optimize_mmff_single,
};

pub use mmff::optimization::{MmffEnergyGradient, MmffEvaluationParams, evaluate_mmff};

pub use uff::{UffEnergyGradient, UffEvaluationError, UffEvaluationParams, evaluate_uff};
