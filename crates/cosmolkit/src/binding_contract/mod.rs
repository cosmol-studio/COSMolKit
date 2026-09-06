//! Machine-readable cross-language public API contract.

mod registry;
mod types;

pub use registry::BINDING_CONTRACT;
pub use types::{BindingContractEntry, BindingKind, BindingOwner, ReturnKind, StateModel};
