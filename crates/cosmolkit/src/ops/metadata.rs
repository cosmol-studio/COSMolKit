//! Stable operation metadata types.

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct MoleculeOpSpec {
    pub method: &'static str,
    pub operation: &'static str,
    pub output: OperationOutput,
    pub access: BlockAccess,
    pub topology_edit: TopologyEditKind,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum OperationOutput {
    Single,
    Multiple,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AccessMode {
    None,
    Read,
    Write,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BlockAccess {
    pub topology: AccessMode,
    pub coordinates: AccessMode,
    pub properties: AccessMode,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TopologyEditKind {
    None,
    Local,
    Compacting,
    Appending,
    Renumbering,
}
