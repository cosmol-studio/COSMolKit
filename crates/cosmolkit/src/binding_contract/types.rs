//! Public types used by the binding contract registry.

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BindingOwner {
    Molecule,
    Module,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BindingKind {
    Instance,
    Static,
    Module,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReturnKind {
    ResultMolecule,
    ResultUnit,
    ResultBool,
    ResultF64,
    ResultString,
    String,
    Unit,
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
    pub owner: BindingOwner,
    pub kind: BindingKind,
    pub rust_name: &'static str,
    pub python_name: &'static str,
    pub javascript_name: &'static str,
    pub feature: &'static str,
    pub return_kind: ReturnKind,
    pub state_model: StateModel,
}
