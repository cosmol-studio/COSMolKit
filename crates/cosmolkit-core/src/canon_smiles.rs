use crate::{Molecule, SmilesWriteError, SmilesWriteParams};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AtomColor {
    Initial,
    Unique,
    Duplicate,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MolStackType {
    Atom,
    Bond,
    Branch,
    RingClosure,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MolStackElem {
    Atom(usize),
    Bond(usize),
    Branch,
    RingClosure(usize),
}

impl MolStackElem {
    #[must_use]
    pub const fn kind(&self) -> MolStackType {
        match self {
            Self::Atom(_) => MolStackType::Atom,
            Self::Bond(_) => MolStackType::Bond,
            Self::Branch => MolStackType::Branch,
            Self::RingClosure(_) => MolStackType::RingClosure,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct FragmentTraversal {
    pub atoms: Vec<usize>,
    pub bonds: Vec<usize>,
}

pub fn rank_mol_atoms(_molecule: &Molecule) -> Result<Vec<usize>, SmilesWriteError> {
    Err(crate::UnsupportedFeatureError::from_spec(&crate::SMILES_WRITE_FEATURE).into())
}

pub fn build_noncanonical_fragment(
    _molecule: &Molecule,
    _params: &SmilesWriteParams,
) -> Result<String, SmilesWriteError> {
    Err(crate::UnsupportedFeatureError::from_spec(&crate::SMILES_WRITE_FEATURE).into())
}

pub fn canonicalize_fragment(
    _molecule: &Molecule,
    _params: &SmilesWriteParams,
) -> Result<String, SmilesWriteError> {
    Err(crate::UnsupportedFeatureError::from_spec(&crate::SMILES_WRITE_FEATURE).into())
}
