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

mod distgeom;

mod crystalff;

pub use crystalff::{
    CrystalFFDetails, CrystalffTorsionPreferencesError, get_experimental_torsions_without_bonds,
};
pub use distgeom::{
    ChiralSet, ChiralSetPtr, ChiralSetStructureFlags, ConformerOptimizer, ConformerOptimizerError,
    DistanceBoundsRead, DistanceGeometryForceFieldParams, calc_chiral_volume_rows,
};

pub use uff::{UffBoundsError, uff_bond_rest_lengths};

pub use uff::{MissingExplicitHydrogensError, needs_explicit_hydrogens};

#[doc(hidden)]
pub use crystalff::{
    CrystalTorsionEvaluationError, CrystalTorsionPairEvaluation, evaluate_crystal_torsion_pair,
};
