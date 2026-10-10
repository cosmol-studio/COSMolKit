//! Machine-readable cross-language public API contract.

mod registry;
mod types;

pub use registry::{
    BINDING_CONTRACT, BINDING_CONTRACT_KEYWORDS, BINDING_CONTRACT_PROPERTIES,
    BINDING_CONTRACT_PYTHON_ADAPTERS, BINDING_CONTRACT_PYTHON_ALIASES,
    BINDING_CONTRACT_PYTHON_COLLECTIONS,
};
pub use types::{
    BindingCallableContract, BindingConfigurationField, BindingContractEntry, BindingDefault,
    BindingItem, BindingKeywordContract, BindingKind, BindingOwner, BindingParameterContract,
    BindingPropertyAccess, BindingPropertyContract, BindingPythonAdapterContract, BindingReceiver,
    BindingTypeRole, FunctionStatus, StateModel,
};
