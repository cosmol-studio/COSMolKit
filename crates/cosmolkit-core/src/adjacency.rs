/// Adjacency entry: neighboring atom and corresponding bond index.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct NeighborRef {
    pub atom_index: usize,
    pub bond_index: usize,
}

/// Graph adjacency list indexed by atom index.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct AdjacencyList {
    pub neighbors: Vec<Vec<NeighborRef>>,
}

impl AdjacencyList {
    /// Build adjacency representation from atom and bond tables.
    #[must_use]
    pub fn from_topology(atom_count: usize, bonds: &[crate::Bond]) -> Self {
        let mut neighbors = vec![Vec::new(); atom_count];
        for bond in bonds {
            neighbors[bond.begin_atom].push(NeighborRef {
                atom_index: bond.end_atom,
                bond_index: bond.index,
            });
            neighbors[bond.end_atom].push(NeighborRef {
                atom_index: bond.begin_atom,
                bond_index: bond.index,
            });
        }
        Self { neighbors }
    }

    /// Get neighbor list for one atom.
    #[must_use]
    pub fn neighbors_of(&self, atom_index: usize) -> &[NeighborRef] {
        &self.neighbors[atom_index]
    }
}
