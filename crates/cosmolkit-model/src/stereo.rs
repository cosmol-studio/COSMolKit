//! Ordered tetrahedral configurations and explicitly undetermined centers.
use crate::AtomId;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum LigandRef {
    /// An atom in the same molecule's index space.
    Atom(AtomId),
    /// A hydrogen without a separate atom row.
    ImplicitHydrogen,
    /// A lone-pair position at a perceived pyramidal center.
    LonePair,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TetrahedralStereo {
    pub center: AtomId,
    /// Ordered configuration, or None for a potential center without a determined configuration.
    pub ligands: Option<[LigandRef; 4]>,
}
