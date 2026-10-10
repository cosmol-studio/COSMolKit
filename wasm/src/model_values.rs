//! Detached snapshots obtained through the canonical public molecule API.
use crate::Molecule;
use cosmolkit as ck;

impl Molecule {
    pub fn from_parts(
        topology: ck::TopologyBlock,
        coordinates: ck::CoordinateBlock,
        properties: ck::MoleculeProperties,
    ) -> Result<Self, ck::OperationError> {
        ck::Molecule::from_parts(topology, coordinates, properties).map(Into::into)
    }
    pub fn to_builder(&self) -> ck::MoleculeBuilder {
        self.inner.borrow().to_builder()
    }
    pub fn topology(&self) -> ck::TopologyBlock {
        self.inner.borrow().topology().clone()
    }
    pub fn properties(&self) -> ck::MoleculeProperties {
        self.inner.borrow().properties().clone()
    }
    pub fn property(&self, key: &str) -> Option<ck::PropertyValue> {
        self.inner.borrow().property(key).cloned()
    }
    pub fn conformers_3d(&self) -> Vec<ck::Conformer3D> {
        self.inner.borrow().conformers_3d().to_vec()
    }
    pub fn has_2d_coordinates(&self) -> bool {
        self.inner.borrow().has_2d_coordinates()
    }
    pub fn has_valence_violation(&self, atom_id: ck::AtomId) -> Result<bool, ck::ValenceError> {
        self.inner.borrow().has_valence_violation(atom_id)
    }
    pub fn atoms(&self) -> Vec<ck::Atom> {
        self.inner.borrow().atoms().to_vec()
    }
    pub fn bonds(&self) -> Vec<ck::Bond> {
        self.inner.borrow().bonds().to_vec()
    }
    pub fn atom(&self, atom: ck::AtomId) -> Option<ck::Atom> {
        self.inner.borrow().atom(atom).cloned()
    }
    pub fn bond(&self, bond: ck::BondId) -> Option<ck::Bond> {
        self.inner.borrow().bond(bond).cloned()
    }
    #[cfg(feature = "cap-valence")]
    pub fn atom_metadata(
        &self,
        recalculate: bool,
    ) -> Result<Vec<ck::AtomMetadata>, ck::ValenceError> {
        self.inner.borrow().atom_metadata(recalculate)
    }
}
