//! Runtime-owned live molecule state.

use cosmolkit_model::{
    CoordinateBlock, MoleculeProperties, TopologyBlock, TopologyValidationError,
};

use crate::ops::OperationError;

/// The authoritative molecule value exposed by the top-level crate.
///
/// The model blocks are deliberately private. Algorithms can only receive
/// detached blocks through an operation context, and only this runtime can
/// install a validated result back into a live molecule.
#[derive(Clone, Debug, PartialEq)]
pub struct Molecule {
    topology: TopologyBlock,
    coordinates: CoordinateBlock,
    properties: MoleculeProperties,
}

impl Molecule {
    /// Creates an empty, structurally valid molecule.
    #[must_use]
    pub fn new() -> Self {
        Self {
            topology: TopologyBlock::default(),
            coordinates: CoordinateBlock::default(),
            properties: MoleculeProperties::default(),
        }
    }

    /// Constructs a live molecule after validating all local block invariants.
    pub fn from_parts(
        topology: TopologyBlock,
        coordinates: CoordinateBlock,
        properties: MoleculeProperties,
    ) -> Result<Self, OperationError> {
        topology
            .validate()
            .map_err(|error: TopologyValidationError| OperationError::InvalidTopology(error))?;
        coordinates
            .validate_for_atom_count(topology.atoms.len())
            .map_err(OperationError::InvalidCoordinates)?;
        Ok(Self {
            topology,
            coordinates,
            properties,
        })
    }

    /// Returns the immutable topology value.
    #[must_use]
    pub const fn topology(&self) -> &TopologyBlock {
        &self.topology
    }

    /// Returns the immutable coordinate value.
    #[must_use]
    pub const fn coordinates(&self) -> &CoordinateBlock {
        &self.coordinates
    }

    /// Returns the immutable molecule properties.
    #[must_use]
    pub const fn properties(&self) -> &MoleculeProperties {
        &self.properties
    }

    /// Returns the number of atoms in the authoritative topology.
    #[must_use]
    pub fn num_atoms(&self) -> usize {
        self.topology.atoms.len()
    }

    /// Returns the number of bonds in the authoritative topology.
    #[must_use]
    pub fn num_bonds(&self) -> usize {
        self.topology.bonds.len()
    }

    #[cfg(feature = "hydrogens")]
    pub(crate) fn install_parts(
        &mut self,
        topology: TopologyBlock,
        coordinates: CoordinateBlock,
        properties: MoleculeProperties,
    ) -> Result<(), OperationError> {
        let replacement = Self::from_parts(topology, coordinates, properties)?;
        *self = replacement;
        Ok(())
    }
}

impl Default for Molecule {
    fn default() -> Self {
        Self::new()
    }
}
