use cosmolkit_core::{
    RingFindType, RingFindingError, RingInfo, fast_find_rings, fast_find_rings_from_parts,
};
use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, StereoGroup, StereoGroupKind,
    SubstanceGroup, SubstanceGroupId, SubstanceGroupKind, TopologyBlock, TopologyValidationError,
};
use cosmolkit_types::{BondOrder, Element};

fn topology(atom_count: usize, edges: &[(usize, usize, BondOrder)]) -> TopologyBlock {
    let atoms = (0..atom_count)
        .map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C)))
        .collect::<Vec<_>>();
    let bonds = edges
        .iter()
        .enumerate()
        .map(|(index, &(begin, end, order))| {
            Bond::from_spec(
                BondId::new(index),
                BondSpec::new(AtomId::new(begin), AtomId::new(end), order),
            )
        })
        .collect::<Vec<_>>();
    TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
}

fn ids(values: &[usize]) -> Vec<AtomId> {
    values.iter().copied().map(AtomId::new).collect()
}

fn bond_ids(values: &[usize]) -> Vec<BondId> {
    values.iter().copied().map(BondId::new).collect()
}

fn assert_fast_info(
    info: &RingInfo,
    atom_count: usize,
    bond_count: usize,
    atom_rows: &[&[usize]],
    bond_rows: &[&[usize]],
) {
    let expected_atom_rows = atom_rows.iter().map(|row| ids(row)).collect::<Vec<_>>();
    let expected_bond_rows = bond_rows
        .iter()
        .map(|row| bond_ids(row))
        .collect::<Vec<_>>();

    assert!(info.is_initialized());
    assert_eq!(info.find_type(), RingFindType::Fast);
    assert!(info.is_find_fast_or_better());
    assert!(!info.is_sssr_or_better());
    assert_eq!(info.num_rings(), atom_rows.len());
    assert_eq!(info.atom_rings(), expected_atom_rows.as_slice());
    assert_eq!(info.bond_rings(), expected_bond_rows.as_slice());
    assert_eq!(info.num_ring_families(), 0);
    assert!(info.atom_ring_families().is_empty());
    assert!(info.bond_ring_families().is_empty());

    for atom in 0..atom_count {
        let expected_members = atom_rows
            .iter()
            .enumerate()
            .filter_map(|(ring, row)| row.contains(&atom).then_some(ring))
            .collect::<Vec<_>>();
        let expected_sizes = expected_members
            .iter()
            .map(|ring| atom_rows[*ring].len())
            .collect::<Vec<_>>();
        assert_eq!(info.atom_members(AtomId::new(atom)), expected_members);
        assert_eq!(
            info.num_atom_rings(AtomId::new(atom)),
            expected_members.len()
        );
        assert_eq!(info.atom_ring_sizes(AtomId::new(atom)), expected_sizes);
    }
    for bond in 0..bond_count {
        let expected_members = bond_rows
            .iter()
            .enumerate()
            .filter_map(|(ring, row)| row.contains(&bond).then_some(ring))
            .collect::<Vec<_>>();
        let expected_sizes = expected_members
            .iter()
            .map(|ring| bond_rows[*ring].len())
            .collect::<Vec<_>>();
        assert_eq!(info.bond_members(BondId::new(bond)), expected_members);
        assert_eq!(
            info.num_bond_rings(BondId::new(bond)),
            expected_members.len()
        );
        assert_eq!(info.bond_ring_sizes(BondId::new(bond)), expected_sizes);
    }
}

#[test]
fn source_empty_isolated_chain_and_tree_cases_are_acyclic() {
    let empty = topology(0, &[]);
    assert_fast_info(&fast_find_rings(&empty).unwrap(), 0, 0, &[], &[]);

    let isolated = topology(3, &[]);
    assert_fast_info(&fast_find_rings(&isolated).unwrap(), 3, 0, &[], &[]);

    let chain = topology(3, &[(0, 1, BondOrder::Single), (1, 2, BondOrder::Single)]);
    assert_fast_info(&fast_find_rings(&chain).unwrap(), 3, 2, &[], &[]);

    let tree = topology(
        5,
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (1, 3, BondOrder::Single),
            (3, 4, BondOrder::Single),
        ],
    );
    assert_fast_info(&fast_find_rings(&tree).unwrap(), 5, 4, &[], &[]);
}

#[test]
fn source_triangle_tail_and_branch_rows_preserve_dfs_and_bond_order() {
    let triangle = topology(
        3,
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (2, 0, BondOrder::Single),
        ],
    );
    assert_fast_info(
        &fast_find_rings(&triangle).unwrap(),
        3,
        3,
        &[&[2, 1, 0]],
        &[&[1, 0, 2]],
    );

    let tail_triangle = topology(
        4,
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Single),
            (3, 1, BondOrder::Single),
        ],
    );
    assert_fast_info(
        &fast_find_rings(&tail_triangle).unwrap(),
        4,
        4,
        &[&[3, 2, 1]],
        &[&[2, 1, 3]],
    );

    let branch_triangle = topology(
        4,
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (1, 3, BondOrder::Single),
            (3, 0, BondOrder::Single),
        ],
    );
    assert_fast_info(
        &fast_find_rings(&branch_triangle).unwrap(),
        4,
        4,
        &[&[3, 1, 0]],
        &[&[2, 0, 3]],
    );
}

#[test]
fn disconnected_cycles_follow_ascending_component_traversal_order() {
    let input = topology(
        6,
        &[
            (3, 4, BondOrder::Single),
            (0, 1, BondOrder::Single),
            (4, 5, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (5, 3, BondOrder::Single),
            (2, 0, BondOrder::Single),
        ],
    );
    assert_fast_info(
        &fast_find_rings(&input).unwrap(),
        6,
        6,
        &[&[2, 1, 0], &[5, 4, 3]],
        &[&[3, 1, 5], &[2, 0, 4]],
    );
}

#[test]
fn source_fused_heteroaromatic_graph_produces_two_complete_fast_rows() {
    // Topology of c1c(=O)nc2[nH]cnn2c1O in source atom/bond creation order.
    let input = topology(
        11,
        &[
            (0, 1, BondOrder::Aromatic),
            (1, 2, BondOrder::Double),
            (1, 3, BondOrder::Aromatic),
            (3, 4, BondOrder::Aromatic),
            (4, 5, BondOrder::Aromatic),
            (5, 6, BondOrder::Aromatic),
            (6, 7, BondOrder::Aromatic),
            (7, 8, BondOrder::Aromatic),
            (8, 4, BondOrder::Aromatic),
            (8, 9, BondOrder::Aromatic),
            (9, 0, BondOrder::Aromatic),
            (9, 10, BondOrder::Single),
        ],
    );
    assert_fast_info(
        &fast_find_rings(&input).unwrap(),
        11,
        12,
        &[&[8, 7, 6, 5, 4], &[9, 8, 7, 6, 5, 4, 3, 1, 0]],
        &[&[7, 6, 5, 4, 8], &[9, 7, 6, 5, 4, 3, 2, 0, 10]],
    );
}

#[test]
fn spiro_fused_and_bridged_graphs_keep_source_fast_rows() {
    let spiro = topology(
        5,
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (2, 0, BondOrder::Single),
            (0, 3, BondOrder::Single),
            (3, 4, BondOrder::Single),
            (4, 0, BondOrder::Single),
        ],
    );
    assert_fast_info(
        &fast_find_rings(&spiro).unwrap(),
        5,
        6,
        &[&[2, 1, 0], &[4, 3, 0]],
        &[&[1, 0, 2], &[4, 3, 5]],
    );

    let fused = topology(
        4,
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (2, 0, BondOrder::Single),
            (1, 3, BondOrder::Single),
            (3, 2, BondOrder::Single),
        ],
    );
    assert_fast_info(
        &fast_find_rings(&fused).unwrap(),
        4,
        5,
        &[&[2, 1, 0], &[3, 2, 1]],
        &[&[1, 0, 2], &[4, 1, 3]],
    );

    let bridged = topology(
        4,
        &[
            (0, 1, BondOrder::Single),
            (0, 2, BondOrder::Single),
            (2, 1, BondOrder::Single),
            (0, 3, BondOrder::Single),
            (3, 1, BondOrder::Single),
        ],
    );
    assert_fast_info(
        &fast_find_rings(&bridged).unwrap(),
        4,
        5,
        &[&[2, 1, 0], &[3, 1, 0]],
        &[&[2, 0, 1], &[4, 0, 3]],
    );
}

#[test]
fn zero_dative_and_hydrogen_bonds_are_not_filtered() {
    for special_order in [BondOrder::Zero, BondOrder::Dative, BondOrder::Hydrogen] {
        let input = topology(
            3,
            &[
                (0, 1, special_order),
                (1, 2, BondOrder::Single),
                (2, 0, BondOrder::Single),
            ],
        );
        assert_fast_info(
            &fast_find_rings(&input).unwrap(),
            3,
            3,
            &[&[2, 1, 0]],
            &[&[1, 0, 2]],
        );
    }
}

#[test]
fn repeated_calls_and_validated_parts_adapter_are_equal_and_non_mutating() {
    let input = topology(
        3,
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (2, 0, BondOrder::Single),
        ],
    );
    let before = input.clone();
    let first = fast_find_rings(&input).unwrap();
    let second = fast_find_rings(&input).unwrap();
    let from_parts =
        fast_find_rings_from_parts(input.atoms.len(), &input.bonds, &input.adjacency).unwrap();

    assert_eq!(first, second);
    assert_eq!(first, from_parts);
    assert_eq!(input, before);
}

#[test]
fn malformed_topology_is_rejected_with_exact_structured_errors() {
    let atom_id_mismatch = TopologyBlock {
        atoms: vec![Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C))],
        adjacency: AdjacencyList::from_topology(1, &[]),
        ..TopologyBlock::default()
    };
    assert_eq!(
        fast_find_rings(&atom_id_mismatch),
        Err(RingFindingError::InvalidTopology(
            TopologyValidationError::AtomIdMismatch {
                position: 0,
                id: AtomId::new(1),
            }
        ))
    );

    let bond_id_mismatch = TopologyBlock {
        atoms: topology(2, &[]).atoms,
        bonds: vec![Bond::from_spec(
            BondId::new(2),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        )],
        ..TopologyBlock::default()
    };
    assert_eq!(
        fast_find_rings(&bond_id_mismatch),
        Err(RingFindingError::InvalidTopology(
            TopologyValidationError::BondIdMismatch {
                position: 0,
                id: BondId::new(2),
            }
        ))
    );

    let endpoint_out_of_range = TopologyBlock {
        atoms: topology(1, &[]).atoms,
        bonds: vec![Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single),
        )],
        ..TopologyBlock::default()
    };
    assert_eq!(
        fast_find_rings(&endpoint_out_of_range),
        Err(RingFindingError::InvalidTopology(
            TopologyValidationError::BondEndpointOutOfRange {
                bond: BondId::new(0),
                endpoint: "end",
                atom: AtomId::new(2),
                atom_count: 1,
            }
        ))
    );

    let self_loop_bond = Bond::from_spec(
        BondId::new(0),
        BondSpec::new(AtomId::new(0), AtomId::new(0), BondOrder::Single),
    );
    let self_loop = TopologyBlock {
        atoms: topology(1, &[]).atoms,
        adjacency: AdjacencyList::from_topology(1, std::slice::from_ref(&self_loop_bond)),
        bonds: vec![self_loop_bond],
        ..TopologyBlock::default()
    };
    assert_eq!(
        fast_find_rings(&self_loop),
        Err(RingFindingError::InvalidTopology(
            TopologyValidationError::SelfLoopBond {
                bond: BondId::new(0),
                atom: AtomId::new(0),
            }
        ))
    );

    let duplicate_edge = TopologyBlock {
        atoms: topology(2, &[]).atoms,
        bonds: vec![
            Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            ),
            Bond::from_spec(
                BondId::new(1),
                BondSpec::new(AtomId::new(1), AtomId::new(0), BondOrder::Double),
            ),
        ],
        ..TopologyBlock::default()
    };
    assert_eq!(
        fast_find_rings(&duplicate_edge),
        Err(RingFindingError::InvalidTopology(
            TopologyValidationError::AdjacencyMismatch
        ))
    );

    let valid_triangle = topology(
        3,
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (2, 0, BondOrder::Single),
        ],
    );
    let stale_adjacency = TopologyBlock {
        adjacency: AdjacencyList::from_topology(3, &[]),
        ..valid_triangle
    };
    assert_eq!(
        fast_find_rings(&stale_adjacency),
        Err(RingFindingError::InvalidTopology(
            TopologyValidationError::AdjacencyMismatch
        ))
    );
}

#[test]
fn invalid_stereo_and_sgroup_references_fail_before_traversal() {
    let stereo_bond = Bond::from_spec(
        BondId::new(0),
        BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Double)
            .with_stereo_atoms(AtomId::new(2), AtomId::new(3)),
    );
    let invalid_bond_stereo = TopologyBlock {
        atoms: topology(2, &[]).atoms,
        adjacency: AdjacencyList::from_topology(2, std::slice::from_ref(&stereo_bond)),
        bonds: vec![stereo_bond],
        ..TopologyBlock::default()
    };
    assert_eq!(
        fast_find_rings(&invalid_bond_stereo),
        Err(RingFindingError::InvalidTopology(
            TopologyValidationError::StereoAtomOutOfRange {
                bond: BondId::new(0),
                begin: AtomId::new(2),
                end: AtomId::new(3),
                atom_count: 2,
            }
        ))
    );

    let invalid_stereo_group = TopologyBlock {
        stereo_groups: vec![StereoGroup::new(
            StereoGroupKind::Absolute,
            vec![AtomId::new(4)],
            vec![],
        )],
        ..topology(1, &[])
    };
    assert_eq!(
        fast_find_rings(&invalid_stereo_group),
        Err(RingFindingError::InvalidTopology(
            TopologyValidationError::StereoGroupAtomOutOfRange {
                atom: AtomId::new(4),
                atom_count: 1,
            }
        ))
    );

    let invalid_sgroup = TopologyBlock {
        substance_groups: vec![
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_atoms(vec![AtomId::new(5)]),
        ],
        ..topology(1, &[])
    };
    assert_eq!(
        fast_find_rings(&invalid_sgroup),
        Err(RingFindingError::InvalidTopology(
            TopologyValidationError::SubstanceGroupAtomOutOfRange {
                sgroup: SubstanceGroupId::new(0),
                atom: AtomId::new(5),
                atom_count: 1,
            }
        ))
    );
}
