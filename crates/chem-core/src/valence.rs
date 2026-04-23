use crate::Molecule;

/// Valence model selector for future compatibility modes.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum ValenceModel {
    /// Intended RDKit-like valence behavior.
    RdkitLike,
}

/// Per-atom valence assignment result.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValenceAssignment {
    pub explicit_valence: Vec<u8>,
    pub implicit_hydrogens: Vec<u8>,
}

/// Valence handling errors.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ValenceError {
    InvalidValence {
        atom_index: usize,
        atomic_num: u8,
        formal_charge: i8,
    },
    NotImplemented,
}

/// Compute explicit valence and implicit hydrogen assignment.
pub fn assign_valence(
    molecule: &Molecule,
    model: ValenceModel,
) -> Result<ValenceAssignment, ValenceError> {
    let _ = (molecule, model);
    unimplemented!("basic valence handling is not implemented yet")
}
