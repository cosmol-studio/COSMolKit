//! Runtime-owned live molecule state.

use std::fmt;
use std::sync::Arc;

use cosmolkit_model::{
    Atom, AtomId, Bond, BondId, Conformer2D, Conformer3D, CoordinateBlock, CoordinateDimension,
    MoleculeProperties, SdfPropertyListTarget, TopologyBlock, TopologyValidationError,
};

use crate::MoleculeBuilder;
use crate::ops::{DerivedState, OperationError};

/// Private runtime authority for derived-state cache storage.
///
/// RUN-effects extends the contents and transition rules. Keeping the cache in
/// the live state now ensures detached model values cannot acquire cache
/// authority and that state replacement also separates cache ownership.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub(crate) struct DerivedCacheBlock {
    valid: DerivedState,
}

impl DerivedCacheBlock {
    fn is_empty(&self) -> bool {
        self.valid == DerivedState::NONE
    }

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

#[derive(Clone)]
struct MoleculeState {
    topology: Arc<TopologyBlock>,
    coordinates: Arc<CoordinateBlock>,
    properties: Arc<MoleculeProperties>,
    derived_cache: Arc<DerivedCacheBlock>,
}

impl MoleculeState {
    fn validate_parts(
        topology: &TopologyBlock,
        coordinates: &CoordinateBlock,
        properties: &MoleculeProperties,
    ) -> Result<(), OperationError> {
        topology
            .validate()
            .map_err(|error: TopologyValidationError| OperationError::InvalidTopology(error))?;
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
        Ok(())
    }

    fn try_new(
        topology: TopologyBlock,
        mut coordinates: CoordinateBlock,
        properties: MoleculeProperties,
    ) -> Result<Self, OperationError> {
        coordinates.source_coordinate_dim = if coordinates.conformers_3d.is_empty() {
            (!coordinates.conformers_2d.is_empty()).then_some(CoordinateDimension::TwoD)
        } else {
            Some(CoordinateDimension::ThreeD)
        };
        Self::validate_parts(&topology, &coordinates, &properties)?;

        Ok(Self {
            topology: Arc::new(topology),
            coordinates: Arc::new(coordinates),
            properties: Arc::new(properties),
            derived_cache: Arc::new(DerivedCacheBlock::default()),
        })
    }

    fn try_new_with_cache(
        topology: Arc<TopologyBlock>,
        coordinates: Arc<CoordinateBlock>,
        properties: Arc<MoleculeProperties>,
        derived_cache: Arc<DerivedCacheBlock>,
    ) -> Result<Self, OperationError> {
        Self::validate_parts(&topology, &coordinates, &properties)?;
        Ok(Self {
            topology,
            coordinates,
            properties,
            derived_cache,
        })
    }
}

/// The authoritative molecule value exposed by the top-level crate.
///
/// The model blocks are deliberately private. Algorithms can only receive
/// detached blocks through an operation context, and only this runtime can
/// install a validated result back into a live molecule.
#[derive(Clone)]
pub struct Molecule {
    state: Arc<MoleculeState>,
}

impl fmt::Debug for Molecule {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("Molecule")
            .field("topology", self.state.topology.as_ref())
            .field("coordinates", self.state.coordinates.as_ref())
            .field("properties", self.state.properties.as_ref())
            .field(
                "derived_cache_is_empty",
                &self.state.derived_cache.is_empty(),
            )
            .finish()
    }
}

impl PartialEq for Molecule {
    fn eq(&self, other: &Self) -> bool {
        self.state.topology == other.state.topology
            && self.state.coordinates == other.state.coordinates
            && self.state.properties == other.state.properties
    }
}

impl Molecule {
    /// Creates an empty, structurally valid molecule.
    #[must_use]
    pub fn new() -> Self {
        Self::from_validated_parts(
            TopologyBlock::default(),
            CoordinateBlock::default(),
            MoleculeProperties::default(),
        )
        .expect("the statically empty molecule blocks are valid")
    }

    /// Constructs a live molecule after validating all local block invariants.
    pub fn from_parts(
        topology: TopologyBlock,
        coordinates: CoordinateBlock,
        properties: MoleculeProperties,
    ) -> Result<Self, OperationError> {
        MoleculeBuilder::from_parts(topology, coordinates, properties).build()
    }

    pub(crate) fn from_validated_parts(
        topology: TopologyBlock,
        coordinates: CoordinateBlock,
        properties: MoleculeProperties,
    ) -> Result<Self, OperationError> {
        Ok(Self {
            state: Arc::new(MoleculeState::try_new(topology, coordinates, properties)?),
        })
    }

    pub(crate) fn from_runtime_parts(
        topology: Arc<TopologyBlock>,
        coordinates: Arc<CoordinateBlock>,
        properties: Arc<MoleculeProperties>,
        derived_cache: Arc<DerivedCacheBlock>,
    ) -> Result<Self, OperationError> {
        Ok(Self {
            state: Arc::new(MoleculeState::try_new_with_cache(
                topology,
                coordinates,
                properties,
                derived_cache,
            )?),
        })
    }

    /// Returns detached semantic blocks for checked construction of a new value.
    #[must_use]
    pub fn to_builder(&self) -> MoleculeBuilder {
        MoleculeBuilder::from_parts(
            self.state.topology.as_ref().clone(),
            self.state.coordinates.as_ref().clone(),
            self.state.properties.as_ref().clone(),
        )
    }

    /// Returns the immutable topology value.
    #[must_use]
    pub fn topology(&self) -> &TopologyBlock {
        self.state.topology.as_ref()
    }

    /// Returns the immutable coordinate value.
    #[must_use]
    pub fn coordinates(&self) -> &CoordinateBlock {
        self.state.coordinates.as_ref()
    }

    /// Returns the immutable molecule properties.
    #[must_use]
    pub fn properties(&self) -> &MoleculeProperties {
        self.state.properties.as_ref()
    }

    /// Returns the number of atoms in the authoritative topology.
    #[must_use]
    pub fn num_atoms(&self) -> usize {
        self.state.topology.atoms.len()
    }

    /// Returns the number of bonds in the authoritative topology.
    #[must_use]
    pub fn num_bonds(&self) -> usize {
        self.state.topology.bonds.len()
    }

    /// Returns all atoms in canonical row order.
    #[must_use]
    pub fn atoms(&self) -> &[Atom] {
        &self.state.topology.atoms
    }

    /// Returns all bonds in canonical row order.
    #[must_use]
    pub fn bonds(&self) -> &[Bond] {
        &self.state.topology.bonds
    }

    /// Returns the atom with the requested stable row identifier.
    #[must_use]
    pub fn atom(&self, atom_id: AtomId) -> Option<&Atom> {
        self.state.topology.atoms.get(atom_id.index())
    }

    /// Returns the bond with the requested stable row identifier.
    #[must_use]
    pub fn bond(&self, bond_id: BondId) -> Option<&Bond> {
        self.state.topology.bonds.get(bond_id.index())
    }

    /// Returns the ordered 2D and 3D conformer slices.
    #[must_use]
    pub fn conformers(&self) -> (&[Conformer2D], &[Conformer3D]) {
        (
            &self.state.coordinates.conformers_2d,
            &self.state.coordinates.conformers_3d,
        )
    }

    /// Returns an ordinary molecule property by key.
    #[must_use]
    pub fn property(&self, key: &str) -> Option<&str> {
        self.state.properties.prop(key)
    }

    pub(crate) fn derived_cache_runtime(&self) -> &DerivedCacheBlock {
        self.state.derived_cache.as_ref()
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

    #[cfg(test)]
    pub(crate) const fn runtime_constructions(&self) -> usize {
        1
    }
}

impl Default for Molecule {
    fn default() -> Self {
        Self::new()
    }
}
