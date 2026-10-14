use crate::{AtomId, Bond, BondId};
use std::collections::HashMap;

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum AdjacencyError {
    #[error("bond {bond} {endpoint} atom index {atom} is out of range for {atom_count} atoms")]
    BondAtomOutOfRange {
        bond: BondId,
        endpoint: &'static str,
        atom: AtomId,
        atom_count: usize,
    },
    #[error("bond id {bond} is repeated at positions {first_position} and {second_position}")]
    DuplicateBondId {
        bond: BondId,
        first_position: usize,
        second_position: usize,
    },
    #[error("bonds {first_bond} and {second_bond} repeat edge {begin}-{end}")]
    DuplicateEdge {
        first_bond: BondId,
        second_bond: BondId,
        begin: AtomId,
        end: AtomId,
    },
}

#[cfg(all(test, feature = "runtime-invariants"))]
mod validation_tests {
    use super::*;
    use crate::{BondOrder, BondSpec};

    fn bonds(edges: &[(usize, usize)]) -> Vec<Bond> {
        edges
            .iter()
            .enumerate()
            .map(|(id, &(begin, end))| {
                Bond::from_spec(
                    BondId::new(id),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
                )
            })
            .collect()
    }

    fn check_against_rebuild(adjacency: &AdjacencyList, atom_count: usize, bonds: &[Bond]) {
        let expected = AdjacencyList::try_from_topology(atom_count, bonds)
            .is_ok_and(|rebuilt| *adjacency == rebuilt);
        assert_eq!(
            adjacency.matches_validated_topology(atom_count, bonds),
            expected,
            "atoms={atom_count}, bonds={bonds:?}, adjacency={adjacency:?}"
        );
    }

    #[test]
    fn direct_validation_equals_rebuild_for_all_four_atom_graphs_and_row_orders() {
        let edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
        for mask in 0..(1 << edges.len()) {
            let selected: Vec<_> = edges
                .iter()
                .enumerate()
                .filter(|(index, _)| mask & (1 << index) != 0)
                .map(|(_, edge)| *edge)
                .collect();
            for reverse_rows in [false, true] {
                for reverse_endpoints in [false, true] {
                    let mut selected = selected.clone();
                    if reverse_rows {
                        selected.reverse();
                    }
                    if reverse_endpoints {
                        selected
                            .iter_mut()
                            .for_each(|edge| *edge = (edge.1, edge.0));
                    }
                    let bonds = bonds(&selected);
                    let adjacency = AdjacencyList::from_topology(4, &bonds);
                    check_against_rebuild(&adjacency, 4, &bonds);
                }
            }
        }
        for atom_count in [0, 1, 8] {
            check_against_rebuild(
                &AdjacencyList::from_topology(atom_count, &[]),
                atom_count,
                &[],
            );
            check_against_rebuild(&AdjacencyList::default(), atom_count, &[]);
        }
    }

    #[test]
    fn direct_validation_rejects_corrupt_offsets_entries_and_neighbor_order() {
        let bonds = bonds(&[(2, 0), (0, 1), (3, 0), (1, 3)]);
        let original = AdjacencyList::from_topology(5, &bonds);
        for index in 0..original.offsets.len() {
            for value in [
                0,
                1,
                original.entries.len(),
                original.entries.len() + 1,
                usize::MAX,
            ] {
                let mut changed = original.clone();
                changed.offsets[index] = value;
                check_against_rebuild(&changed, 5, &bonds);
            }
        }
        for index in 0..original.entries.len() {
            for value in [0, 1, 2, 3, 4, 5, usize::MAX] {
                let mut changed = original.clone();
                changed.entries[index].atom_index = value;
                check_against_rebuild(&changed, 5, &bonds);
                changed = original.clone();
                changed.entries[index].bond = BondId::new(value);
                check_against_rebuild(&changed, 5, &bonds);
            }
            if index + 1 < original.entries.len() {
                let mut changed = original.clone();
                changed.entries.swap(index, index + 1);
                check_against_rebuild(&changed, 5, &bonds);
            }
        }
        let mut changed = original.clone();
        changed.entries.pop();
        check_against_rebuild(&changed, 5, &bonds);
        changed = original.clone();
        changed.entries.push(original.entries[0]);
        check_against_rebuild(&changed, 5, &bonds);
        changed = original.clone();
        changed.offsets.pop();
        check_against_rebuild(&changed, 5, &bonds);
        changed = original;
        changed.offsets.push(changed.entries.len());
        check_against_rebuild(&changed, 5, &bonds);
    }

    #[test]
    fn direct_validation_rejects_duplicate_edges_with_otherwise_matching_csr() {
        for edges in [[(0, 1), (0, 1)], [(0, 1), (1, 0)]] {
            let bonds = bonds(&edges);
            let adjacency = AdjacencyList {
                offsets: vec![0, 2, 4],
                entries: vec![
                    NeighborRef {
                        atom_index: 1,
                        bond: BondId::new(0),
                    },
                    NeighborRef {
                        atom_index: 1,
                        bond: BondId::new(1),
                    },
                    NeighborRef {
                        atom_index: 0,
                        bond: BondId::new(0),
                    },
                    NeighborRef {
                        atom_index: 0,
                        bond: BondId::new(1),
                    },
                ],
            };
            check_against_rebuild(&adjacency, 2, &bonds);
        }
    }

    #[test]
    fn direct_validation_handles_high_degree_and_missing_endpoint_rows() {
        let edges: Vec<_> = (1..1024).map(|end| (0, end)).collect();
        let bonds = bonds(&edges);
        let adjacency = AdjacencyList::from_topology(1024, &bonds);
        check_against_rebuild(&adjacency, 1024, &bonds);
        let mut missing_endpoint = adjacency.clone();
        missing_endpoint.entries[1023] = missing_endpoint.entries[1024];
        check_against_rebuild(&missing_endpoint, 1024, &bonds);
    }
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

#[cfg(test)]
mod construction_tests {
    use super::*;
    use crate::{BondOrder, BondSpec};

    fn bond(id: usize, begin: usize, end: usize) -> Bond {
        Bond::from_spec(
            BondId::new(id),
            BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
        )
    }

    #[test]
    fn checked_csr_keeps_input_order_for_dense_and_sparse_identifiers() {
        for ids in [[0, 1, 2], [5000, 2, usize::MAX]] {
            let bonds = [bond(ids[0], 2, 0), bond(ids[1], 0, 1), bond(ids[2], 2, 3)];
            let adjacency = AdjacencyList::try_from_topology(5, &bonds).unwrap();
            let row = |atom| {
                adjacency
                    .neighbors_of(atom)
                    .iter()
                    .map(|neighbor| (neighbor.atom_index, neighbor.bond.index()))
                    .collect::<Vec<_>>()
            };
            assert_eq!(row(0), [(2, ids[0]), (1, ids[1])]);
            assert_eq!(row(1), [(0, ids[1])]);
            assert_eq!(row(2), [(0, ids[0]), (3, ids[2])]);
            assert_eq!(row(3), [(2, ids[2])]);
            assert!(row(4).is_empty());
        }
        assert_eq!(AdjacencyList::from_topology(3, &[]).offsets, [0, 0, 0, 0]);
    }

    #[test]
    fn checked_csr_preserves_first_bond_and_within_bond_error_priority() {
        let duplicate = AdjacencyError::DuplicateEdge {
            first_bond: BondId::new(7),
            second_bond: BondId::new(8),
            begin: AtomId::new(0),
            end: AtomId::new(1),
        };
        assert_eq!(
            AdjacencyList::try_from_topology(3, &[bond(7, 0, 1), bond(8, 1, 0), bond(9, 4, 5)]),
            Err(duplicate),
        );
        assert_eq!(
            AdjacencyList::try_from_topology(3, &[bond(7, 0, 1), bond(7, 1, 0)]),
            Err(AdjacencyError::DuplicateBondId {
                bond: BondId::new(7),
                first_position: 0,
                second_position: 1,
            }),
        );
        assert_eq!(
            AdjacencyList::try_from_topology(3, &[bond(7, 0, 1), bond(7, 4, 5)]),
            Err(AdjacencyError::BondAtomOutOfRange {
                bond: BondId::new(7),
                endpoint: "begin",
                atom: AtomId::new(4),
                atom_count: 3,
            }),
        );
        assert_eq!(
            AdjacencyList::try_from_topology(3, &[bond(7, 0, 1), bond(7, 0, 5)]),
            Err(AdjacencyError::BondAtomOutOfRange {
                bond: BondId::new(7),
                endpoint: "end",
                atom: AtomId::new(5),
                atom_count: 3,
            }),
        );
    }
}

impl AdjacencyList {
    /// Compare with the canonical CSR representation without rebuilding it.
    /// The caller has already checked dense bond IDs, endpoints and self-loops.
    #[cfg(feature = "runtime-invariants")]
    pub(crate) fn matches_validated_topology(&self, atom_count: usize, bonds: &[Bond]) -> bool {
        if self.offsets.len().checked_sub(1) != Some(atom_count)
            || self.offsets.first() != Some(&0)
            || self.offsets.last() != Some(&self.entries.len())
            || bonds.len().checked_mul(2) != Some(self.entries.len())
        {
            return false;
        }

        // One reusable O(V) marker array replaces two trees and a temporary
        // CSR. Each row must contain distinct neighbors and ascending bond
        // IDs, exactly as try_from_topology emits them. Every valid entry is
        // one endpoint of its bond; strict ordering permits each endpoint at
        // most once. The 2E entry count therefore also proves none is missing.
        // This is O(V + E), including high-degree graphs, without sorting or
        // changing the observable neighbor order or validation policy.
        let mut seen_neighbors = vec![usize::MAX; atom_count];
        for (atom_index, window) in self.offsets.windows(2).enumerate() {
            let Some(row) = self.entries.get(window[0]..window[1]) else {
                return false;
            };
            let mut previous_bond = None;
            for entry in row {
                let bond_index = entry.bond.index();
                if previous_bond.is_some_and(|previous| bond_index <= previous) {
                    return false;
                }
                let Some(bond) = bonds.get(bond_index) else {
                    return false;
                };
                if !((bond.begin().index() == atom_index && bond.end().index() == entry.atom_index)
                    || (bond.end().index() == atom_index
                        && bond.begin().index() == entry.atom_index))
                {
                    return false;
                }
                let Some(seen) = seen_neighbors.get_mut(entry.atom_index) else {
                    return false;
                };
                if *seen == atom_index {
                    return false;
                }
                *seen = atom_index;
                previous_bond = Some(bond_index);
            }
        }
        true
    }

    pub fn try_from_topology(atom_count: usize, bonds: &[Bond]) -> Result<Self, AdjacencyError> {
        let mut degrees = vec![0usize; atom_count];
        let mut dense_ids = true;
        // Size the final CSR without returning errors yet. Validation below
        // still observes each bond in input order: bounds, ID, then edge.
        // Canonical parser IDs are proven unique by this pass, not assumed.
        for (position, bond) in bonds.iter().enumerate() {
            let begin = bond.begin().index();
            let end = bond.end().index();
            if begin >= atom_count || end >= atom_count {
                // Later rows cannot affect the first failing prefix. Earlier
                // duplicate diagnostics still take precedence in pass two.
                break;
            }
            dense_ids &= bond.id().index() == position;
            degrees[begin] += 1;
            degrees[end] += 1;
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
        // Detached callers may supply arbitrary IDs. Only that case needs a
        // lookup table; its iteration order never affects CSR or diagnostics.
        let mut bond_positions = (!dense_ids).then(HashMap::new);
        for (position, bond) in bonds.iter().enumerate() {
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
            if let Some(first_position) = bond_positions
                .as_mut()
                .and_then(|positions| positions.insert(bond.id(), position))
            {
                return Err(AdjacencyError::DuplicateBondId {
                    bond: bond.id(),
                    first_position,
                    second_position: position,
                });
            }
            // BEGIN RDKIT CPP FUNCTION RWMol::addBond duplicate-edge check
            // RDKit✔️✔️:   PRECONDITION(!(boost::edge(atomIdx1, atomIdx2, d_graph).second),
            // RDKit✔️✔️:                "bond already exists");
            // END RDKIT CPP FUNCTION RWMol::addBond duplicate-edge check
            // Search only the begin atom's already-filled adjacency, like
            // boost::edge on RDKit's vecS graph. No global ordered edge index.
            // For canonical IDs this is O(V+E+sum previous begin degrees),
            // matching the source and linear on bounded-degree graphs.
            if let Some(previous) = entries[offsets[begin.index()]..cursor[begin.index()]]
                .iter()
                .find(|neighbor| neighbor.atom_index == end.index())
            {
                let (edge_begin, edge_end) = if begin <= end {
                    (begin, end)
                } else {
                    (end, begin)
                };
                return Err(AdjacencyError::DuplicateEdge {
                    first_bond: previous.bond,
                    second_bond: bond.id(),
                    begin: edge_begin,
                    end: edge_end,
                });
            }
            let begin = begin.index();
            let end = end.index();
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
    #[doc(hidden)]
    pub fn from_topology(atom_count: usize, bonds: &[Bond]) -> Self {
        Self::try_from_topology(atom_count, bonds)
            .expect("topology must be valid before building adjacency")
    }

    /// Borrow an actual CSR row, preserving absence separately from zero degree.
    #[must_use]
    #[doc(hidden)]
    pub fn try_neighbors_of(&self, atom_index: usize) -> Option<&[NeighborRef]> {
        let end = atom_index.checked_add(2)?;
        let window = self.offsets.get(atom_index..end)?;
        let [start, end] = window else { return None };
        self.entries.get(*start..*end)
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
