use crate::{AtomId, Bond, BondId};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum AdjacencyError {
    #[error("bond {bond} {endpoint} atom index {atom} is out of range for {atom_count} atoms")]
    BondAtomOutOfRange {
        bond: BondId,
        endpoint: &'static str,
        atom: AtomId,
        atom_count: usize,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NeighborRef {
    pub atom_index: usize,
    pub bond: BondId,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct AdjacencyList {
    offsets: Vec<usize>,
    entries: Vec<NeighborRef>,
}

impl AdjacencyList {
    pub fn try_from_topology(atom_count: usize, bonds: &[Bond]) -> Result<Self, AdjacencyError> {
        let mut degrees = vec![0usize; atom_count];
        for bond in bonds {
            let begin = bond.begin();
            let end = bond.end();
            if begin.index() >= atom_count {
                return Err(AdjacencyError::BondAtomOutOfRange {
                    bond: bond.id(),
                    endpoint: "begin",
                    atom: begin,
                    atom_count,
                });
            }
            if end.index() >= atom_count {
                return Err(AdjacencyError::BondAtomOutOfRange {
                    bond: bond.id(),
                    endpoint: "end",
                    atom: end,
                    atom_count,
                });
            }
            degrees[begin.index()] += 1;
            degrees[end.index()] += 1;
        }

        let mut offsets = Vec::with_capacity(atom_count + 1);
        offsets.push(0);
        for degree in degrees {
            offsets.push(offsets.last().copied().expect("offsets is nonempty") + degree);
        }

        let mut entries = vec![
            NeighborRef {
                atom_index: 0,
                bond: BondId::new(0),
            };
            offsets.last().copied().unwrap_or(0)
        ];
        let mut cursor = offsets[..atom_count].to_vec();
        for bond in bonds {
            let begin = bond.begin().index();
            let end = bond.end().index();
            let begin_slot = cursor[begin];
            entries[begin_slot] = NeighborRef {
                atom_index: end,
                bond: bond.id(),
            };
            cursor[begin] += 1;

            let end_slot = cursor[end];
            entries[end_slot] = NeighborRef {
                atom_index: begin,
                bond: bond.id(),
            };
            cursor[end] += 1;
        }

        Ok(Self { offsets, entries })
    }

    #[must_use]
    pub fn from_topology(atom_count: usize, bonds: &[Bond]) -> Self {
        Self::try_from_topology(atom_count, bonds)
            .expect("topology must be valid before building adjacency")
    }

    #[must_use]
    pub fn neighbors_of(&self, atom_index: usize) -> &[NeighborRef] {
        let Some(window) = self.offsets.get(atom_index..atom_index.saturating_add(2)) else {
            return &[];
        };
        if window.len() != 2 {
            return &[];
        }
        &self.entries[window[0]..window[1]]
    }
}
