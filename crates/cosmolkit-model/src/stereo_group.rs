//! Stereo-group value metadata shared by parsers and molecule algorithms.

use crate::{AtomId, BondId};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StereoGroupKind {
    Absolute,
    Or,
    And,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StereoGroup {
    id: Option<u32>,
    kind: StereoGroupKind,
    atoms: Vec<AtomId>,
    bonds: Vec<BondId>,
}

impl StereoGroup {
    #[must_use]
    pub fn new(kind: StereoGroupKind, atoms: Vec<AtomId>, bonds: Vec<BondId>) -> Self {
        Self {
            id: None,
            kind,
            atoms,
            bonds,
        }
    }

    #[must_use]
    pub const fn with_id(mut self, id: u32) -> Self {
        self.id = Some(id);
        self
    }

    #[must_use]
    pub const fn id(&self) -> Option<u32> {
        self.id
    }

    #[must_use]
    pub const fn kind(&self) -> StereoGroupKind {
        self.kind
    }

    #[must_use]
    pub fn atoms(&self) -> &[AtomId] {
        &self.atoms
    }

    #[must_use]
    pub fn bonds(&self) -> &[BondId] {
        &self.bonds
    }

    pub fn push_atom(&mut self, atom: AtomId) {
        self.atoms.push(atom);
    }

    pub fn remove_atom(&mut self, atom: AtomId) {
        self.atoms.retain(|candidate| *candidate != atom);
    }

    pub fn remove_bond(&mut self, bond: BondId) {
        self.bonds.retain(|candidate| *candidate != bond);
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty() && self.bonds.is_empty()
    }

    pub fn remapped(&self, atom_map: &[Option<AtomId>], bond_map: &[Option<BondId>]) -> Option<Self> {
        let atoms: Option<Vec<_>> = self
            .atoms
            .iter()
            .map(|atom| atom_map.get(atom.index()).and_then(|x| *x))
            .collect();
        let bonds: Option<Vec<_>> = self
            .bonds
            .iter()
            .map(|bond| bond_map.get(bond.index()).and_then(|x| *x))
            .collect();
        Some(Self {
            id: self.id,
            kind: self.kind,
            atoms: atoms?,
            bonds: bonds?,
        })
    }
}
