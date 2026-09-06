//! Detached block extraction and installation for operation bodies.

#[cfg(feature = "hydrogens")]
use std::marker::PhantomData;

#[cfg(feature = "hydrogens")]
use cosmolkit_model::{CoordinateBlock, MoleculeProperties, TopologyBlock};

#[cfg(feature = "hydrogens")]
use super::{MoleculeOpSpec, OperationError, OperationOutput};
#[cfg(feature = "hydrogens")]
use crate::Molecule;

#[cfg(feature = "hydrogens")]
pub(crate) struct UnrestrictedAccess;
#[cfg(feature = "hydrogens")]
pub(crate) struct AddHydrogensAccess;
#[cfg(feature = "hydrogens")]
pub(crate) struct RemoveHydrogensAccess;

/// Runtime capability projection for one registered single-output operation.
///
/// The access marker is selected by the operation body macro. The live
/// molecule is copied into detached model blocks and can only be replaced by
/// the runtime after all blocks have been returned and validated.
#[cfg(feature = "hydrogens")]
pub(crate) struct OpParts<'a, Access = UnrestrictedAccess> {
    pub(crate) spec: &'static MoleculeOpSpec,
    pub(crate) topology: Option<TopologyBlock>,
    pub(crate) coordinates: Option<CoordinateBlock>,
    pub(crate) properties: Option<MoleculeProperties>,
    pub(crate) in_place_target: Option<&'a mut Molecule>,
    pub(crate) _access: PhantomData<Access>,
}

#[cfg(feature = "hydrogens")]
impl<'a, Access> OpParts<'a, Access> {
    pub(crate) fn new(
        source: &'a Molecule,
        spec: &'static MoleculeOpSpec,
    ) -> Result<Self, OperationError> {
        if spec.output != OperationOutput::Single {
            return Err(OperationError::Unsupported {
                operation: spec.operation,
            });
        }
        Ok(Self {
            spec,
            topology: Some(source.topology().clone()),
            coordinates: Some(source.coordinates().clone()),
            properties: Some(source.properties().clone()),
            in_place_target: None,
            _access: PhantomData,
        })
    }

    pub(crate) fn new_in_place(
        target: &'a mut Molecule,
        spec: &'static MoleculeOpSpec,
    ) -> Result<Self, OperationError> {
        if spec.output != OperationOutput::Single {
            return Err(OperationError::Unsupported {
                operation: spec.operation,
            });
        }
        Ok(Self {
            spec,
            topology: Some(target.topology().clone()),
            coordinates: Some(target.coordinates().clone()),
            properties: Some(target.properties().clone()),
            in_place_target: Some(target),
            _access: PhantomData,
        })
    }

    pub(crate) fn take_all(
        &mut self,
    ) -> Result<(TopologyBlock, CoordinateBlock, MoleculeProperties), OperationError> {
        let topology = self
            .topology
            .take()
            .ok_or(OperationError::BlockCheckedOut { block: "topology" })?;
        let coordinates = self
            .coordinates
            .take()
            .ok_or(OperationError::BlockCheckedOut {
                block: "coordinates",
            })?;
        let properties = self
            .properties
            .take()
            .ok_or(OperationError::BlockCheckedOut {
                block: "properties",
            })?;
        Ok((topology, coordinates, properties))
    }

    pub(crate) fn put_all(
        &mut self,
        topology: TopologyBlock,
        coordinates: CoordinateBlock,
        properties: MoleculeProperties,
    ) -> Result<(), OperationError> {
        topology
            .validate()
            .map_err(OperationError::InvalidTopology)?;
        coordinates
            .validate_for_atom_count(topology.atoms.len())
            .map_err(OperationError::InvalidCoordinates)?;
        self.topology = Some(topology);
        self.coordinates = Some(coordinates);
        self.properties = Some(properties);
        Ok(())
    }

    pub(crate) fn finish(self) -> Result<Molecule, OperationError> {
        let Self {
            topology,
            coordinates,
            properties,
            ..
        } = self;
        let (topology, coordinates, properties) = topology
            .zip(coordinates)
            .zip(properties)
            .map(|((topology, coordinates), properties)| (topology, coordinates, properties))
            .ok_or(OperationError::IncompleteCommit { block: "molecule" })?;
        Molecule::from_parts(topology, coordinates, properties)
    }

    pub(crate) fn finish_in_place(mut self) -> Result<(), OperationError> {
        let target = self
            .in_place_target
            .take()
            .ok_or(OperationError::IncompleteCommit { block: "target" })?;
        let topology = self
            .topology
            .take()
            .ok_or(OperationError::IncompleteCommit { block: "topology" })?;
        let coordinates = self
            .coordinates
            .take()
            .ok_or(OperationError::IncompleteCommit {
                block: "coordinates",
            })?;
        let properties = self
            .properties
            .take()
            .ok_or(OperationError::IncompleteCommit {
                block: "properties",
            })?;
        target.install_parts(topology, coordinates, properties)
    }
}
