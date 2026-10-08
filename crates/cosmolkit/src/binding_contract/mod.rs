//! Machine-readable cross-language public API contract.

mod registry;
mod types;

pub use registry::{
    BINDING_CONTRACT, BINDING_CONTRACT_KEYWORDS, BINDING_CONTRACT_PROPERTIES,
    BINDING_CONTRACT_PYTHON_ADAPTERS,
};
pub use types::{
    BindingCallableContract, BindingContractEntry, BindingDefault, BindingItem,
    BindingKeywordContract, BindingKind, BindingOwner, BindingParameterContract,
    BindingPropertyAccess, BindingPropertyContract, BindingPythonAdapterContract, BindingReceiver,
    BindingTypeRole, FunctionStatus, StateModel,
};
