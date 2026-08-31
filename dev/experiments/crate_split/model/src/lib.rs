use std::fmt;

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct AtomId(usize);

impl AtomId {
    #[must_use]
    pub const fn from_index(index: usize) -> Self {
        Self(index)
    }

    #[must_use]
    pub const fn index(self) -> usize {
        self.0
    }
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum Element {
    H,
    C,
    N,
    O,
    Unknown,
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Atom {
    element: Element,
    formal_charge: i8,
}

impl Atom {
    #[must_use]
    pub const fn new(element: Element) -> Self {
        Self {
            element,
            formal_charge: 0,
        }
    }

    #[must_use]
    pub const fn element(&self) -> Element {
        self.element
    }

    #[must_use]
    pub const fn formal_charge(&self) -> i8 {
        self.formal_charge
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Bond {
    source: AtomId,
    target: AtomId,
    order: BondOrder,
}

impl Bond {
    #[must_use]
    pub const fn new(source: AtomId, target: AtomId, order: BondOrder) -> Self {
        Self {
            source,
            target,
            order,
        }
    }

    #[must_use]
    pub const fn source(&self) -> AtomId {
        self.source
    }

    #[must_use]
    pub const fn target(&self) -> AtomId {
        self.target
    }

    #[must_use]
    pub const fn order(&self) -> BondOrder {
        self.order
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TopologyBlock {
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TopologyError {
    InvalidBondEndpoint,
}

impl fmt::Display for TopologyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidBondEndpoint => f.write_str("bond endpoint is outside topology"),
        }
    }
}

impl std::error::Error for TopologyError {}

impl TopologyBlock {
    #[must_use]
    pub fn new() -> Self {
        Self {
            atoms: Vec::new(),
            bonds: Vec::new(),
        }
    }

    #[must_use]
    pub fn atoms(&self) -> &[Atom] {
        &self.atoms
    }

    #[must_use]
    pub fn bonds(&self) -> &[Bond] {
        &self.bonds
    }

    #[must_use]
    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    #[must_use]
    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }

    pub fn add_atom(&mut self, atom: Atom) -> AtomId {
        let id = AtomId(self.atoms.len());
        self.atoms.push(atom);
        id
    }

    pub fn add_bond(
        &mut self,
        source: AtomId,
        target: AtomId,
        order: BondOrder,
    ) -> Result<(), TopologyError> {
        if source.index() >= self.atoms.len() || target.index() >= self.atoms.len() {
            return Err(TopologyError::InvalidBondEndpoint);
        }
        self.bonds.push(Bond::new(source, target, order));
        Ok(())
    }

    pub fn validate(&self) -> Result<(), TopologyError> {
        // This experiment only models the local bond-endpoint check. A
        // production crate split must carry over the complete molecule
        // invariant set, including atom/bond identity, adjacency, coordinate
        // row alignment, and per-entity property alignment.
        for bond in &self.bonds {
            if bond.source().index() >= self.atoms.len()
                || bond.target().index() >= self.atoms.len()
            {
                return Err(TopologyError::InvalidBondEndpoint);
            }
        }
        Ok(())
    }
}

impl Default for TopologyBlock {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Conformer3D {
    coordinates: Vec<[f64; 3]>,
}

impl Conformer3D {
    #[must_use]
    pub fn new(coordinates: Vec<[f64; 3]>) -> Self {
        Self { coordinates }
    }

    #[must_use]
    pub fn coordinates(&self) -> &[[f64; 3]] {
        &self.coordinates
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct CoordinateBlock {
    conformers_3d: Vec<Conformer3D>,
}

impl CoordinateBlock {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    #[must_use]
    pub fn conformers_3d(&self) -> &[Conformer3D] {
        &self.conformers_3d
    }

    pub fn add_conformer_3d(&mut self, conformer: Conformer3D) {
        self.conformers_3d.push(conformer);
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct MoleculeProperties {
    name: Option<String>,
}

impl MoleculeProperties {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    #[must_use]
    pub fn name(&self) -> Option<&str> {
        self.name.as_deref()
    }

    #[must_use]
    pub fn with_name(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct QueryAst {
    source: String,
}

impl QueryAst {
    #[must_use]
    pub fn new(source: impl Into<String>) -> Self {
        Self {
            source: source.into(),
        }
    }

    #[must_use]
    pub fn source(&self) -> &str {
        &self.source
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn topology_builder_keeps_local_endpoint_invariant() {
        let mut topology = TopologyBlock::new();
        let carbon = topology.add_atom(Atom::new(Element::C));
        let oxygen = topology.add_atom(Atom::new(Element::O));
        topology
            .add_bond(carbon, oxygen, BondOrder::Single)
            .expect("valid endpoints");
        assert_eq!(topology.atom_count(), 2);
        assert_eq!(topology.bond_count(), 1);
        assert!(topology.validate().is_ok());
    }
}
