#![cfg(not(feature = "hydrogens"))]

use std::sync::Arc;

use cosmolkit_model::{CoordinateBlock, MoleculeProperties, SdfPropertyListTarget, TopologyBlock};

pub use cosmolkit::ops::{
    CipStatePolicy, DerivedEffects, DerivedState, MappingRequirement, ParityPolicy,
    SemanticPreconditionSet,
};
pub use cosmolkit::{
    BlockAccess, BlockSet, MoleculeOpKind, MoleculeOpOutput, MoleculeOpSpec, OperationDomain,
    OperationError, SupportStatus, TopologyEditKind,
};

mod ops {
    pub use crate::{
        BlockAccess, BlockSet, CipStatePolicy, DerivedEffects, DerivedState, MappingRequirement,
        MoleculeOpKind, MoleculeOpOutput, MoleculeOpSpec, OperationDomain, OperationError,
        ParityPolicy, SemanticPreconditionSet, SupportStatus, TopologyEditKind,
    };
}

mod strict {
    pub(crate) const RUNTIME_INVARIANTS_ENABLED: bool = cfg!(feature = "runtime-invariants");
    pub(crate) const OPERATION_CONTRACTS_ENABLED: bool = cfg!(feature = "op-contracts");
}

mod molecule {
    use super::*;

    #[derive(Clone, Debug, Default, Eq, PartialEq)]
    pub(crate) struct DerivedCacheBlock {
        valid: DerivedState,
    }

    impl DerivedCacheBlock {
        pub(crate) fn valid_states(&self) -> DerivedState {
            self.valid
        }

        pub(crate) fn mark_valid(&mut self, states: DerivedState) {
            self.valid = self.valid.union(states);
        }

        pub(crate) fn clear(&mut self, states: DerivedState) {
            self.valid = self.valid.difference(states);
        }
    }

    #[derive(Clone, Debug, PartialEq)]
    struct MoleculeState {
        topology: TopologyBlock,
        coordinates: CoordinateBlock,
        properties: MoleculeProperties,
        derived_cache: DerivedCacheBlock,
        runtime_constructions: usize,
    }

    #[derive(Clone, Debug, PartialEq)]
    pub struct Molecule {
        state: Arc<MoleculeState>,
    }

    impl Molecule {
        pub fn from_parts(
            topology: TopologyBlock,
            coordinates: CoordinateBlock,
            properties: MoleculeProperties,
        ) -> Result<Self, OperationError> {
            Self::validated(
                topology,
                coordinates,
                properties,
                DerivedCacheBlock::default(),
                0,
            )
        }

        fn validated(
            topology: TopologyBlock,
            coordinates: CoordinateBlock,
            properties: MoleculeProperties,
            derived_cache: DerivedCacheBlock,
            runtime_constructions: usize,
        ) -> Result<Self, OperationError> {
            topology
                .validate()
                .map_err(OperationError::InvalidTopology)?;
            coordinates
                .validate_for_atom_count(topology.atoms.len())
                .map_err(OperationError::InvalidCoordinates)?;
            for property_list in properties.sdf_property_lists() {
                let (target, expected) = match property_list.target() {
                    SdfPropertyListTarget::Atom => ("atom", topology.atoms.len()),
                    SdfPropertyListTarget::Bond => ("bond", topology.bonds.len()),
                };
                if property_list.values().len() != expected {
                    return Err(OperationError::InvalidPropertyList {
                        target,
                        name: property_list.name().to_owned(),
                        values: property_list.values().len(),
                        expected,
                    });
                }
            }
            Ok(Self {
                state: Arc::new(MoleculeState {
                    topology,
                    coordinates,
                    properties,
                    derived_cache,
                    runtime_constructions,
                }),
            })
        }

        pub(crate) fn from_runtime_parts(
            topology: Arc<TopologyBlock>,
            coordinates: Arc<CoordinateBlock>,
            properties: Arc<MoleculeProperties>,
            derived_cache: Arc<DerivedCacheBlock>,
        ) -> Result<Self, OperationError> {
            Self::validated(
                topology.as_ref().clone(),
                coordinates.as_ref().clone(),
                properties.as_ref().clone(),
                derived_cache.as_ref().clone(),
                1,
            )
        }

        pub fn topology(&self) -> &TopologyBlock {
            &self.state.topology
        }

        pub(crate) fn coordinate_block_runtime(&self) -> &CoordinateBlock {
            &self.state.coordinates
        }

        pub fn properties(&self) -> &MoleculeProperties {
            &self.state.properties
        }

        pub fn num_atoms(&self) -> usize {
            self.state.topology.atoms.len()
        }

        pub(crate) fn derived_cache_runtime(&self) -> &DerivedCacheBlock {
            &self.state.derived_cache
        }

        pub(crate) fn topology_arc_runtime(&self) -> Arc<TopologyBlock> {
            Arc::new(self.state.topology.clone())
        }

        pub(crate) fn coordinates_arc_runtime(&self) -> Arc<CoordinateBlock> {
            Arc::new(self.state.coordinates.clone())
        }

        pub(crate) fn properties_arc_runtime(&self) -> Arc<MoleculeProperties> {
            Arc::new(self.state.properties.clone())
        }

        pub(crate) fn derived_cache_arc_runtime(&self) -> Arc<DerivedCacheBlock> {
            Arc::new(self.state.derived_cache.clone())
        }

        pub(crate) fn runtime_constructions(&self) -> usize {
            self.state.runtime_constructions
        }
    }
}

pub use molecule::Molecule;

#[path = "../src/ops/context.rs"]
mod context;

#[test]
fn owning_target_uses_the_private_failure_runtime() {
    let source = Molecule::from_parts(
        TopologyBlock::default(),
        CoordinateBlock::default(),
        MoleculeProperties::default(),
    )
    .unwrap();
    assert_eq!(source.runtime_constructions(), 0);
}
