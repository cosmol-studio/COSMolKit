//! Public metadata values used by the canonical binding contract registry.

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BindingItem {
    Callable,
    Type,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BindingOwner {
    Molecule,
    Module,
    Type,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BindingKind {
    Instance,
    Static,
    Module,
}

/// Ownership of an instance callable's receiver, independent of its output.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BindingReceiver {
    Shared,
    Mutable,
    Owned,
}

/// Explicit Python descriptor projection of a Rust accessor.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BindingPropertyAccess {
    Getter,
    Setter,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BindingDefault {
    Required,
    Value(&'static str),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BindingParameterContract {
    pub name: &'static str,
    pub type_name: &'static str,
    pub default: BindingDefault,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BindingCallableContract {
    pub kind: BindingKind,
    pub receiver: Option<BindingReceiver>,
    pub parameters: &'static [BindingParameterContract],
    pub output_type: &'static str,
    pub error_type: Option<&'static str>,
    pub state_model: StateModel,
    pub operation_semantic_id: Option<&'static str>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BindingTypeRole {
    Value,
    Parameter,
    Result,
    Error,
}

/// One manually declared behavior commitment, independent of test execution.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum FunctionStatus {
    Parity {
        reference: &'static str,
    },
    ParityWithDifferences {
        reference: &'static str,
        explanation: &'static str,
    },
    Native,
    Experimental,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum StateModel {
    ValueReturning,
    InPlace,
    ReadOnly,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BindingContractEntry {
    pub semantic_id: &'static str,
    pub item: BindingItem,
    pub owner: BindingOwner,
    pub rust_path: &'static str,
    pub python_name: &'static str,
    pub python_property: Option<BindingPropertyAccess>,
    pub javascript_name: &'static str,
    pub feature: &'static str,
    /// Additional capability selectors required together with the owning feature.
    pub required_capabilities: &'static [&'static str],
    pub status: FunctionStatus,
    pub callable: Option<BindingCallableContract>,
    pub type_role: Option<BindingTypeRole>,
}

/// Read-only Python properties generated from the canonical type declaration.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BindingPropertyContract {
    pub type_semantic_id: &'static str,
    pub name: &'static str,
    pub rust_path: &'static str,
    pub output_type: &'static str,
}
/// Keyword convenience calls construct this registered immutable parameter
/// value, then invoke the registered explicit Rust target.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BindingKeywordContract {
    pub semantic_id: &'static str,
    pub parameters_semantic_id: &'static str,
    pub target_semantic_id: &'static str,
}
