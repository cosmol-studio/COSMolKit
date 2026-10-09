//! Capability-free fixture for the private runtime's standalone test targets.
//! Block storage mirrors the live value's four Arc-backed blocks. No chemistry
//! cache payload exists in this fixture: derived state is synthetic bit metadata.

use std::sync::Arc;

use cosmolkit_model::{
    Conformer3D, CoordinateBlock, MoleculeProperties, SdfPropertyListTarget, TopologyBlock,
};

use crate::{DerivedState, OperationError};

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

    pub(crate) fn validate_for_topology(
        &self,
        _topology: &TopologyBlock,
    ) -> Result<(), OperationError> {
        // State bits have no atom/bond-indexed payload in this synthetic model.
        // The real cache's payload validation remains covered in its owner.
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq)]
struct MoleculeState {
    topology: Arc<TopologyBlock>,
    coordinates: Arc<CoordinateBlock>,
    properties: Arc<MoleculeProperties>,
    derived_cache: Arc<DerivedCacheBlock>,
    runtime_constructions: usize,
}

#[derive(Clone, Debug)]
pub struct Molecule {
    state: Arc<MoleculeState>,
}

impl PartialEq for Molecule {
    fn eq(&self, other: &Self) -> bool {
        // Construction instrumentation is not molecular state. Keep every
        // modeled block (including synthetic derived-cache bits) in equality;
        // tests assert construction counts independently of value preservation.
        self.state.topology == other.state.topology
            && self.state.coordinates == other.state.coordinates
            && self.state.properties == other.state.properties
            && self.state.derived_cache == other.state.derived_cache
    }
}

impl Molecule {
    pub(crate) fn operation_snapshot_runtime(
        &self,
        _preserve_queries: bool,
        _operation: &'static str,
    ) -> Result<Self, OperationError> {
        Ok(self.clone())
    }

    // This synthetic fixture models no descriptor memo. Actual copy/clear and
    // poisoning are tested on the owning runtime Molecule, never this no-op.
    #[cfg(feature = "cap-descriptors")]
    pub(crate) fn take_descriptor_queries_runtime(&mut self) {}
    #[cfg(feature = "cap-descriptors")]
    pub(crate) fn install_descriptor_queries_runtime(&mut self, _memo: ()) {}

    pub fn from_parts(
        topology: TopologyBlock,
        coordinates: CoordinateBlock,
        properties: MoleculeProperties,
    ) -> Result<Self, OperationError> {
        Self::validated(
            Arc::new(topology),
            Arc::new(coordinates),
            Arc::new(properties),
            Arc::new(DerivedCacheBlock::default()),
            0,
        )
    }

    pub(crate) fn from_runtime_parts(
        topology: Arc<TopologyBlock>,
        coordinates: Arc<CoordinateBlock>,
        properties: Arc<MoleculeProperties>,
        derived_cache: Arc<DerivedCacheBlock>,
    ) -> Result<Self, OperationError> {
        Self::validated(topology, coordinates, properties, derived_cache, 1)
    }

    fn validated(
        topology: Arc<TopologyBlock>,
        coordinates: Arc<CoordinateBlock>,
        properties: Arc<MoleculeProperties>,
        derived_cache: Arc<DerivedCacheBlock>,
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
        derived_cache.validate_for_topology(&topology)?;
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

    pub fn topology(&self) -> &TopologyBlock {
        &self.state.topology
    }

    pub(crate) fn coordinate_block_runtime(&self) -> &CoordinateBlock {
        &self.state.coordinates
    }

    pub fn conformers_3d(&self) -> &[Conformer3D] {
        &self.state.coordinates.conformers_3d
    }

    pub fn properties(&self) -> &MoleculeProperties {
        &self.state.properties
    }

    pub fn num_atoms(&self) -> usize {
        self.state.topology.atoms.len()
    }

    pub fn num_bonds(&self) -> usize {
        self.state.topology.bonds.len()
    }

    pub(crate) fn derived_cache_runtime(&self) -> &DerivedCacheBlock {
        &self.state.derived_cache
    }

    pub(crate) fn topology_arc_runtime(&self) -> Arc<TopologyBlock> {
        Arc::clone(&self.state.topology)
    }

    pub(crate) fn coordinates_arc_runtime(&self) -> Arc<CoordinateBlock> {
        Arc::clone(&self.state.coordinates)
    }

    pub(crate) fn properties_arc_runtime(&self) -> Arc<MoleculeProperties> {
        Arc::clone(&self.state.properties)
    }

    pub(crate) fn derived_cache_arc_runtime(&self) -> Arc<DerivedCacheBlock> {
        Arc::clone(&self.state.derived_cache)
    }

    pub(crate) fn runtime_constructions(&self) -> usize {
        self.state.runtime_constructions
    }
}
