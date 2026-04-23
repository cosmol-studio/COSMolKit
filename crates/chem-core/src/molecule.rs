use crate::{AdjacencyList, Atom, Bond};
use core::fmt;

/// Errors returned by SMILES parsing routines.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SmilesParseError {
    /// Parsing is not implemented yet in bootstrap phase.
    NotImplemented,
    /// Concrete parse error from the RDKit-aligned subset parser.
    ParseError(String),
}

impl fmt::Display for SmilesParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotImplemented => f.write_str("SMILES parser is not implemented yet"),
            Self::ParseError(message) => f.write_str(message),
        }
    }
}

impl std::error::Error for SmilesParseError {}

/// Molecular graph with atom/bond tables and adjacency cache.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub adjacency: Option<AdjacencyList>,
}

impl Molecule {
    /// Create an empty molecule.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Add one atom and return its assigned index.
    pub fn add_atom(&mut self, atom: Atom) -> usize {
        let index = self.atoms.len();
        let mut atom = atom;
        atom.index = index;
        self.atoms.push(atom);
        self.adjacency = None;
        index
    }

    /// Add one bond and return its assigned index.
    pub fn add_bond(&mut self, bond: Bond) -> usize {
        let index = self.bonds.len();
        let mut bond = bond;
        bond.index = index;
        self.bonds.push(bond);
        self.adjacency = None;
        index
    }

    /// Rebuild adjacency representation from current topology.
    pub fn rebuild_adjacency(&mut self) {
        self.adjacency = Some(AdjacencyList::from_topology(self.atoms.len(), &self.bonds));
    }

    /// Construct one molecule from a SMILES string.
    pub fn from_smiles(smiles: &str) -> Result<Self, SmilesParseError> {
        crate::smiles::parse_smiles(smiles)
    }

    /// Return atom atomic numbers in atom-index order.
    #[must_use]
    pub fn atomic_numbers(&self) -> Vec<u8> {
        self.atoms.iter().map(|atom| atom.atomic_num).collect()
    }
}
