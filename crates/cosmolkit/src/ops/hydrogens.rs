//! Hydrogen operation bodies; chemistry remains in its final algorithm owner.

use cosmolkit_macros::mol_op_body;

use super::OperationError;

#[mol_op_body(with_hydrogens, parts)]
pub(crate) fn add_hydrogens_impl() -> Result<(), OperationError> {
    let _ = parts;
    Err(OperationError::Unsupported {
        operation: "with_hydrogens",
    })
}

#[mol_op_body(without_hydrogens, parts)]
pub(crate) fn remove_hydrogens_impl() -> Result<(), OperationError> {
    let _ = parts;
    Err(OperationError::Unsupported {
        operation: "without_hydrogens",
    })
}
