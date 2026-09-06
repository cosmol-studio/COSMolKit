//! Hydrogen operation bodies and their public molecule methods.

use cosmolkit_macros::mol_op_body;

use super::OperationError;
use super::context::{AddHydrogensAccess, OpParts, RemoveHydrogensAccess};
use super::registry::{ADD_HYDROGENS_SPEC, REMOVE_HYDROGENS_SPEC};
use crate::Molecule;

#[mol_op_body(add_hydrogens, context)]
fn add_hydrogens_impl() -> Result<(), OperationError> {
    let _ = context;
    Err(OperationError::Unsupported {
        operation: "add_hydrogens",
    })
}

#[mol_op_body(remove_hydrogens, context)]
fn remove_hydrogens_impl() -> Result<(), OperationError> {
    let _ = context;
    Err(OperationError::Unsupported {
        operation: "remove_hydrogens",
    })
}

impl<'a> OpParts<'a, AddHydrogensAccess> {
    pub(crate) fn extract_all_writable(
        &mut self,
    ) -> Result<
        (
            cosmolkit_model::TopologyBlock,
            cosmolkit_model::CoordinateBlock,
            cosmolkit_model::MoleculeProperties,
        ),
        OperationError,
    > {
        self.take_all()
    }

    pub(crate) fn commit_all_writable(
        &mut self,
        topology: cosmolkit_model::TopologyBlock,
        coordinates: cosmolkit_model::CoordinateBlock,
        properties: cosmolkit_model::MoleculeProperties,
    ) -> Result<(), OperationError> {
        self.put_all(topology, coordinates, properties)
    }
}

impl<'a> OpParts<'a, RemoveHydrogensAccess> {
    pub(crate) fn extract_all_writable(
        &mut self,
    ) -> Result<
        (
            cosmolkit_model::TopologyBlock,
            cosmolkit_model::CoordinateBlock,
            cosmolkit_model::MoleculeProperties,
        ),
        OperationError,
    > {
        self.take_all()
    }

    pub(crate) fn commit_all_writable(
        &mut self,
        topology: cosmolkit_model::TopologyBlock,
        coordinates: cosmolkit_model::CoordinateBlock,
        properties: cosmolkit_model::MoleculeProperties,
    ) -> Result<(), OperationError> {
        self.put_all(topology, coordinates, properties)
    }
}

impl Molecule {
    #[cfg(feature = "hydrogens")]
    pub fn with_hydrogens(&self) -> Result<Self, OperationError> {
        let mut context = OpParts::<AddHydrogensAccess>::new(self, &ADD_HYDROGENS_SPEC)?;
        add_hydrogens_impl(&mut context)?;
        context.finish()
    }

    #[cfg(feature = "hydrogens")]
    pub fn add_hydrogens_(&mut self) -> Result<(), OperationError> {
        let mut context = OpParts::<AddHydrogensAccess>::new_in_place(self, &ADD_HYDROGENS_SPEC)?;
        add_hydrogens_impl(&mut context)?;
        context.finish_in_place()
    }

    #[cfg(feature = "hydrogens")]
    pub fn without_hydrogens(&self) -> Result<Self, OperationError> {
        let mut context = OpParts::<RemoveHydrogensAccess>::new(self, &REMOVE_HYDROGENS_SPEC)?;
        remove_hydrogens_impl(&mut context)?;
        context.finish()
    }

    #[cfg(feature = "hydrogens")]
    pub fn remove_hydrogens_(&mut self) -> Result<(), OperationError> {
        let mut context =
            OpParts::<RemoveHydrogensAccess>::new_in_place(self, &REMOVE_HYDROGENS_SPEC)?;
        remove_hydrogens_impl(&mut context)?;
        context.finish_in_place()
    }
}
