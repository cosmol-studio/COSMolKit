//! Project-native ordered tetrahedral ligand records, preserved from d892.
use crate::AtomId;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum LigandRef {
    Atom(AtomId),
    ImplicitHydrogen,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TetrahedralStereo {
    pub center: AtomId,
    pub ligands: [LigandRef; 4],
}
