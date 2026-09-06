//! Operation registry entries.

use super::MoleculeOpSpec;
#[cfg(feature = "hydrogens")]
use super::{AccessMode, BlockAccess, OperationOutput, TopologyEditKind};

#[cfg(feature = "hydrogens")]
pub(crate) const ADD_HYDROGENS_SPEC: MoleculeOpSpec = MoleculeOpSpec {
    method: "with_hydrogens",
    operation: "add_hydrogens",
    output: OperationOutput::Single,
    access: BlockAccess {
        topology: AccessMode::Write,
        coordinates: AccessMode::Write,
        properties: AccessMode::Write,
    },
    topology_edit: TopologyEditKind::Appending,
};

#[cfg(feature = "hydrogens")]
pub(crate) const REMOVE_HYDROGENS_SPEC: MoleculeOpSpec = MoleculeOpSpec {
    method: "without_hydrogens",
    operation: "remove_hydrogens",
    output: OperationOutput::Single,
    access: BlockAccess {
        topology: AccessMode::Write,
        coordinates: AccessMode::Write,
        properties: AccessMode::Write,
    },
    topology_edit: TopologyEditKind::Compacting,
};

/// Registry entries for the runtime-owned operation surface.
pub static MOLECULE_OPS: &[MoleculeOpSpec] = &[
    #[cfg(feature = "hydrogens")]
    ADD_HYDROGENS_SPEC,
    #[cfg(feature = "hydrogens")]
    REMOVE_HYDROGENS_SPEC,
];
