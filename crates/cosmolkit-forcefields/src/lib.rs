//! Detached UFF/MMFF force-field boundaries.

mod geometry;
mod kernel;
mod mmff;
mod optimization_outcome;
mod optimizer;
pub use optimization_outcome::{OptimizationOutcome, OptimizationStatus};
mod uff;
pub use kernel::{
    AngleIndexArgument as ForceFieldAngleIndexArgument,
    AngleRangeBound as ForceFieldAngleRangeBound, BondIndexArgument as ForceFieldBondIndexArgument,
    ForceFieldIndexArgument, ForceFieldKernelError,
    TorsionIndexArgument as ForceFieldTorsionIndexArgument,
};
pub use uff::{
    UffBondMathError, UffBuilderError, UffDefaultTorsionQueryError, UffInversionContributionError,
    UffInversionIndexArgument, UffParamError, UffPreparedValenceField, UffTorsionBondQueryError,
    UffTypingError, UffTypingInput,
};

pub use uff::{
    UffConformerError, UffConformerOptions, UffConformerOutcome, UffSingleError, UffSingleOptions,
    UffSingleOutcome, optimize_uff_conformers_prepared, optimize_uff_single_prepared,
};
pub use uff::{UffParameterError, UffParameterErrorKind, uff_has_all_molecule_params};

pub use mmff::mol_properties::{MmffAtomProperties, MmffMolPropertiesError, MmffVariant};
pub use mmff::properties_api::{
    MmffProperties, MmffPropertiesParams, mmff_has_all_molecule_params, mmff_properties,
};
pub use mmff::{
    MmffAngleLookupError, MmffOopLookupError, MmffParamError, MmffParamParseCause,
    MmffParamParseError, MmffParamTable, MmffTorLookupError, MmffVdwParseError,
};

pub use mmff::optimization::{
    MmffConformerOptimizationParams, MmffConformerOutcomes, MmffOptimizationError,
    MmffOptimizationParams, MmffOptimizeMoleculeConfResult, MmffSingleOutcome,
    optimize_mmff_conformers, optimize_mmff_single,
};

pub use mmff::optimization::{MmffEnergyGradient, MmffEvaluationParams, evaluate_mmff};

pub use uff::{UffEnergyGradient, UffEvaluationError, UffEvaluationParams, evaluate_uff};

mod distgeom;

mod crystalff;

pub use crystalff::{
    BoundsTorsionPath, CrystalFFDetails, CrystalffTorsionPreferencesError, ForcedTorsionKind,
    get_experimental_torsions_without_bonds,
};
pub use distgeom::{
    ChiralSet, ChiralSetPtr, ChiralSetStructureFlags, ConformerOptimizer, ConformerOptimizerCause,
    ConformerOptimizerError, DistGeomForceFieldError, DistanceBoundsRead,
    DistanceGeometryForceFieldParams, calc_chiral_volume_rows,
};

pub use uff::{UffBoundsCause, UffBoundsError, uff_bond_rest_lengths};

pub use mmff::{builder::MmffBuilderError, optimization::Failure as MmffOptimizationFailure};
#[cfg(not(target_family = "wasm"))]
pub use uff::DispatchedUffOptimizationError;
pub use uff::{
    AutomaticForceFieldConstructionError, ForceFieldConstructionError, NonbondedAssemblyError,
    NonbondedFragmentMappingError, OptimizationStageError, SerialConformerOptimizationError,
    SerialCoordinateCountError, SerialResultCapacityError, SerialUffOptimizationError,
    SingleConformerOptimizationError, UffEvaluationFailure, UffPreparedOptimizationError,
    UffThreadCountError,
};
pub use uff::{MissingExplicitHydrogensError, needs_explicit_hydrogens};

#[doc(hidden)]
pub use crystalff::{
    CrystalTorsionEvaluationError, CrystalTorsionPairEvaluation, evaluate_crystal_torsion_pair,
};

mod persistent;
pub use persistent::{
    ForceFieldError, MolecularForceFieldErrorKind, PreparedForceField, prepare_mmff_force_field,
    prepare_uff_force_field,
};

#[doc(hidden)]
pub use mmff::optimization::__configuration_MmffConformerOptimizationParams;
#[doc(hidden)]
pub use mmff::optimization::__configuration_MmffEvaluationParams;
#[doc(hidden)]
pub use mmff::optimization::__configuration_MmffOptimizationParams;
#[doc(hidden)]
pub use mmff::properties_api::__configuration_MmffPropertiesParams;
#[doc(hidden)]
pub use uff::__configuration_UffEvaluationParams;
