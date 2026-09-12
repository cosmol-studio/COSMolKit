use cosmolkit_model::{
    AdjacencyError, AdjacencyList, AtomId, Bond, BondId, BondOrder, BondSpec, NeighborRef,
};

fn bond(id: usize, begin: usize, end: usize) -> Bond {
    Bond::from_spec(
        BondId::new(id),
        BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
    )
}

fn neighbor(atom_index: usize, bond: usize) -> NeighborRef {
    NeighborRef {
        atom_index,
        bond: BondId::new(bond),
    }
}

#[test]
fn empty_isolated_and_missing_rows_are_empty() {
    let empty = AdjacencyList::try_from_topology(0, &[]).unwrap();
    assert!(empty.neighbors_of(0).is_empty());
    assert!(empty.neighbors_of(usize::MAX).is_empty());
    let isolated = AdjacencyList::try_from_topology(2, &[]).unwrap();
    assert!(isolated.neighbors_of(0).is_empty());
    assert!(isolated.neighbors_of(2).is_empty());
}

#[test]
fn chain_preserves_complete_rows_and_input_order() {
    let bonds = [bond(10, 1, 2), bond(11, 0, 1), bond(12, 2, 3)];
    let adjacency = AdjacencyList::try_from_topology(4, &bonds).unwrap();
    assert_eq!(adjacency.neighbors_of(0), &[neighbor(1, 11)]);
    assert_eq!(
        adjacency.neighbors_of(1),
        &[neighbor(2, 10), neighbor(0, 11)]
    );
    assert_eq!(
        adjacency.neighbors_of(2),
        &[neighbor(1, 10), neighbor(3, 12)]
    );
    assert_eq!(adjacency.neighbors_of(3), &[neighbor(2, 12)]);
}

#[test]
fn branch_preserves_complete_rows_and_input_order() {
    let bonds = [bond(20, 0, 3), bond(21, 0, 1), bond(22, 0, 2)];
    let adjacency = AdjacencyList::try_from_topology(4, &bonds).unwrap();
    assert_eq!(
        adjacency.neighbors_of(0),
        &[neighbor(3, 20), neighbor(1, 21), neighbor(2, 22)]
    );
    assert_eq!(adjacency.neighbors_of(1), &[neighbor(0, 21)]);
    assert_eq!(adjacency.neighbors_of(2), &[neighbor(0, 22)]);
    assert_eq!(adjacency.neighbors_of(3), &[neighbor(0, 20)]);
}

#[test]
fn ring_preserves_complete_rows_and_input_order() {
    let bonds = [
        bond(30, 2, 3),
        bond(31, 0, 1),
        bond(32, 3, 0),
        bond(33, 1, 2),
    ];
    let adjacency = AdjacencyList::try_from_topology(4, &bonds).unwrap();
    assert_eq!(
        adjacency.neighbors_of(0),
        &[neighbor(1, 31), neighbor(3, 32)]
    );
    assert_eq!(
        adjacency.neighbors_of(1),
        &[neighbor(0, 31), neighbor(2, 33)]
    );
    assert_eq!(
        adjacency.neighbors_of(2),
        &[neighbor(3, 30), neighbor(1, 33)]
    );
    assert_eq!(
        adjacency.neighbors_of(3),
        &[neighbor(2, 30), neighbor(0, 32)]
    );
}

#[test]
fn both_endpoint_sides_report_exact_errors() {
    let begin = AdjacencyList::try_from_topology(1, &[bond(3, 1, 0)]).unwrap_err();
    assert_eq!(
        begin,
        AdjacencyError::BondAtomOutOfRange {
            bond: BondId::new(3),
            endpoint: "begin",
            atom: AtomId::new(1),
            atom_count: 1,
        }
    );
    let end = AdjacencyList::try_from_topology(1, &[bond(4, 0, 1)]).unwrap_err();
    assert_eq!(
        end,
        AdjacencyError::BondAtomOutOfRange {
            bond: BondId::new(4),
            endpoint: "end",
            atom: AtomId::new(1),
            atom_count: 1,
        }
    );
}

#[test]
fn duplicate_bond_identity_is_structured() {
    let error = AdjacencyList::try_from_topology(3, &[bond(7, 0, 1), bond(7, 1, 2)]).unwrap_err();
    assert_eq!(
        error,
        AdjacencyError::DuplicateBondId {
            bond: BondId::new(7),
            first_position: 0,
            second_position: 1
        }
    );
}

#[test]
fn duplicate_edge_is_normalized_in_both_orientations() {
    for second in [bond(1, 0, 2), bond(1, 2, 0)] {
        let error = AdjacencyList::try_from_topology(3, &[bond(0, 0, 2), second]).unwrap_err();
        assert_eq!(
            error,
            AdjacencyError::DuplicateEdge {
                first_bond: BondId::new(0),
                second_bond: BondId::new(1),
                begin: AtomId::new(0),
                end: AtomId::new(2)
            }
        );
    }
}
