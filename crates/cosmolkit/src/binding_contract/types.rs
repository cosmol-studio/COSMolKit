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
pub enum BindingExposure {
    Registered,
    Public,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BindingKind {
    Instance,
    Static,
    Module,
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

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BindingSupport {
    Unsupported,
    PreservedOnly,
    Experimental,
    Supported,
    SupportedWithRdkitParity,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BindingParity {
    NotApplicable,
    RequiredWhenSupported,
    RequiredNow,
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
    pub javascript_name: &'static str,
    pub feature: &'static str,
    pub exposure: BindingExposure,
    pub support: BindingSupport,
    pub parity: BindingParity,
    pub callable: Option<BindingCallableContract>,
    pub type_role: Option<BindingTypeRole>,
}
