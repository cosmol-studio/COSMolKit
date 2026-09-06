//! Runtime operation system.

mod context;
mod error;
#[cfg(feature = "hydrogens")]
mod hydrogens;
mod metadata;
mod registry;

pub use error::OperationError;
pub use metadata::{AccessMode, BlockAccess, MoleculeOpSpec, OperationOutput, TopologyEditKind};
pub use registry::MOLECULE_OPS;

#[cfg(feature = "hydrogens")]
pub(crate) use context::{AddHydrogensAccess, OpParts, RemoveHydrogensAccess};
