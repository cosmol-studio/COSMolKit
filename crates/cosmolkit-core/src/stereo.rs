use crate::{AtomId, BondId, Molecule};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum StereoError {
    #[error("stereochemistry perception is not implemented")]
    NotImplemented,
    #[error(transparent)]
    UnsupportedFeature(#[from] crate::UnsupportedFeatureError),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LigandRef {
    Atom(AtomId),
    ImplicitHydrogen,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TetrahedralStereo {
    pub center: AtomId,
    pub ligands: [LigandRef; 4],
    pub clockwise: bool,
}

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

    pub(crate) fn push_atom(&mut self, atom: AtomId) {
        self.atoms.push(atom);
    }

    pub(crate) fn remapped(
        &self,
        atom_map: &[Option<AtomId>],
        bond_map: &[Option<BondId>],
    ) -> Option<Self> {
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

pub fn tetrahedral_stereo(_molecule: &Molecule) -> Result<Vec<TetrahedralStereo>, StereoError> {
    Err(crate::UnsupportedFeatureError::from_spec(&crate::STEREO_FEATURE).into())
}

pub fn should_detect_double_bond_stereo(
    _molecule: &Molecule,
    _bond: BondId,
) -> Result<bool, StereoError> {
    Err(crate::UnsupportedFeatureError::from_spec(&crate::STEREO_FEATURE).into())
}
