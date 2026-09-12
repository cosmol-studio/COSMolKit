//! Machine-readable cross-language public API contract.

mod registry;
mod types;

pub use registry::BINDING_CONTRACT;
pub use types::{
    BindingCallableContract, BindingContractEntry, BindingDefault, BindingExposure, BindingItem,
    BindingKind, BindingOwner, BindingParameterContract, BindingParity, BindingSupport,
    BindingTypeRole, StateModel,
};
