#![cfg(not(feature = "cap-hydrogens"))]

use cosmolkit_model::{CoordinateBlock, MoleculeProperties, TopologyBlock};

pub use cosmolkit::ops::{
    CipStatePolicy, DerivedEffects, DerivedState, MappingRequirement, ParityPolicy,
    SemanticPreconditionSet,
};
pub use cosmolkit::{
    BlockAccess, BlockSet, FunctionStatus, MoleculeOpKind, MoleculeOpOutput, MoleculeOpSpec,
    OperationDomain, OperationError, TopologyEditKind,
};

mod ops {
    pub use crate::{
        BlockAccess, BlockSet, CipStatePolicy, DerivedEffects, DerivedState, FunctionStatus,
        MappingRequirement, MoleculeOpKind, MoleculeOpOutput, MoleculeOpSpec, OperationDomain,
        OperationError, ParityPolicy, SemanticPreconditionSet, TopologyEditKind,
    };
}

mod strict {
    pub(crate) const RUNTIME_INVARIANTS_ENABLED: bool = cfg!(feature = "runtime-invariants");
    pub(crate) const OPERATION_CONTRACTS_ENABLED: bool = cfg!(feature = "op-contracts");
}

#[path = "support/minimal_runtime_molecule.rs"]
mod molecule;

pub use molecule::Molecule;

#[path = "../src/ops/context.rs"]
mod context;
pub(crate) use context::{PendingMolecule, PendingResult, ResultFinalizer};

#[test]
fn owning_target_uses_the_private_commit_runtime() {
    let source = Molecule::from_parts(
        TopologyBlock::default(),
        CoordinateBlock::default(),
        MoleculeProperties::default(),
    )
    .unwrap();
    assert_eq!(source.runtime_constructions(), 0);
}
