use cosmolkit_core::{
    RingFindType, RingFindingError, RingInfo, RingSearchParams, find_sssr, find_sssr_from_parts,
    find_sssr_with_options_from_parts, symmetrized_sssr,
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

fn atom_ids(values: &[usize]) -> Vec<AtomId> {
    values.iter().copied().map(AtomId::new).collect()
}

fn bond_ids(values: &[usize]) -> Vec<BondId> {
    values.iter().copied().map(BondId::new).collect()
}

fn assert_complete_info(
    mut info: RingInfo,
    expected_type: RingFindType,
    atom_count: usize,
    bond_count: usize,
    atom_rows: &[&[usize]],
    bond_rows: &[&[usize]],
) {
    let expected_atom_rows = atom_rows
        .iter()
        .map(|row| atom_ids(row))
        .collect::<Vec<_>>();
    let expected_bond_rows = bond_rows
        .iter()
        .map(|row| bond_ids(row))
        .collect::<Vec<_>>();

    assert!(info.is_initialized());
    assert_eq!(info.find_type(), expected_type);
    assert!(info.is_find_fast_or_better());
    assert_eq!(
        info.is_sssr_or_better(),
        expected_type != RingFindType::Fast
    );
    assert_eq!(info.is_symm_sssr(), expected_type == RingFindType::SymmSssr);
    assert_eq!(info.num_rings(), atom_rows.len());
    assert_eq!(info.atom_rings(), expected_atom_rows.as_slice());
    assert_eq!(info.bond_rings(), expected_bond_rows.as_slice());

    for atom in 0..atom_count {
        let members = atom_rows
            .iter()
            .enumerate()
            .filter_map(|(ring, row)| row.contains(&atom).then_some(ring))
            .collect::<Vec<_>>();
        let sizes = members
            .iter()
            .map(|ring| atom_rows[*ring].len())
            .collect::<Vec<_>>();
        assert_eq!(info.atom_members(AtomId::new(atom)), members);
        assert_eq!(info.num_atom_rings(AtomId::new(atom)), members.len());
        assert_eq!(info.atom_ring_sizes(AtomId::new(atom)), sizes);
        assert_eq!(
            info.min_atom_ring_size(AtomId::new(atom)),
            sizes.iter().copied().min().unwrap_or(0)
        );
    }
    for bond in 0..bond_count {
        let members = bond_rows
            .iter()
            .enumerate()
            .filter_map(|(ring, row)| row.contains(&bond).then_some(ring))
            .collect::<Vec<_>>();
        let sizes = members
            .iter()
            .map(|ring| bond_rows[*ring].len())
            .collect::<Vec<_>>();
        assert_eq!(info.bond_members(BondId::new(bond)), members);
        assert_eq!(info.num_bond_rings(BondId::new(bond)), members.len());
        assert_eq!(info.bond_ring_sizes(BondId::new(bond)), sizes);
        assert_eq!(
            info.min_bond_ring_size(BondId::new(bond)),
            sizes.iter().copied().min().unwrap_or(0)
        );
    }

    for left in 0..atom_count {
        for right in 0..atom_count {
            let shared = atom_rows
                .iter()
                .any(|row| row.contains(&left) && row.contains(&right));
            assert_eq!(
                info.are_atoms_in_same_ring(AtomId::new(left), AtomId::new(right)),
                shared
            );
        }
    }
    for left in 0..bond_count {
        for right in 0..bond_count {
            let shared = bond_rows
                .iter()
                .any(|row| row.contains(&left) && row.contains(&right));
            assert_eq!(
                info.are_bonds_in_same_ring(BondId::new(left), BondId::new(right)),
                shared
            );
        }
    }

    for ring in 0..bond_rows.len() {
        let expected_fused_bonds = bond_rows[ring]
            .iter()
            .filter(|bond| bond_rows.iter().filter(|row| row.contains(bond)).count() > 1)
            .count();
        assert_eq!(info.num_fused_bonds(ring).unwrap(), expected_fused_bonds);
        assert_eq!(
            info.is_ring_fused(ring).unwrap(),
            (0..bond_rows.len()).any(|other| {
                other != ring
                    && bond_rows[ring]
                        .iter()
                        .any(|bond| bond_rows[other].contains(bond))
            })
        );
        for other in 0..bond_rows.len() {
            let shared_bond = ring != other
                && bond_rows[ring]
                    .iter()
                    .any(|bond| bond_rows[other].contains(bond));
            assert_eq!(info.are_rings_fused(ring, other).unwrap(), shared_bond);
        }
    }
}

fn params(dative: bool, hydrogen: bool) -> RingSearchParams {
    RingSearchParams {
        include_dative_bonds: dative,
        include_hydrogen_bonds: hydrogen,
    }
}

#[test]
fn defaults_and_acyclic_inputs_return_initialized_empty_assignments() {
    assert_eq!(RingSearchParams::default(), params(false, false));
    for input in [
        topology(0, &[]),
        topology(3, &[]),
        topology(
            5,
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (1, 3, BondOrder::Single),
                (3, 4, BondOrder::Single),
            ],
        ),
    ] {
        assert_complete_info(
            find_sssr(&input, &RingSearchParams::default()).unwrap(),
            RingFindType::Sssr,
            input.atoms.len(),
            input.bonds.len(),
            &[],
            &[],
        );
        assert_complete_info(
            symmetrized_sssr(&input, &RingSearchParams::default()).unwrap(),
            RingFindType::SymmSssr,
            input.atoms.len(),
            input.bonds.len(),
            &[],
            &[],
        );
    }
}

#[test]
fn source_triangle_square_and_seven_ring_rows_are_exact() {
    for (size, expected_atoms, expected_bonds) in [
        (3, vec![0, 1, 2], vec![0, 1, 2]),
        (4, vec![0, 3, 2, 1], vec![3, 2, 1, 0]),
        (7, vec![0, 1, 2, 3, 4, 5, 6], vec![0, 1, 2, 3, 4, 5, 6]),
    ] {
        let mut edges = (0..size - 1)
            .map(|atom| (atom, atom + 1, BondOrder::Single))
            .collect::<Vec<_>>();
        edges.push((size - 1, 0, BondOrder::Single));
        let input = topology(size, &edges);
        assert_complete_info(
            find_sssr(&input, &RingSearchParams::default()).unwrap(),
            RingFindType::Sssr,
            size,
            size,
            &[expected_atoms.as_slice()],
            &[expected_bonds.as_slice()],
        );
    }
}

#[test]
fn disconnected_components_and_non_ring_rows_keep_source_order() {
    let input = topology(
        8,
        &[
            (4, 5, BondOrder::Single),
            (0, 1, BondOrder::Single),
            (5, 6, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (6, 4, BondOrder::Single),
            (2, 0, BondOrder::Single),
            (2, 3, BondOrder::Single),
            (6, 7, BondOrder::Single),
        ],
    );
    assert_complete_info(
        find_sssr(&input, &RingSearchParams::default()).unwrap(),
        RingFindType::Sssr,
        8,
        8,
        &[&[0, 1, 2], &[4, 5, 6]],
        &[&[1, 3, 5], &[0, 2, 4]],
    );
}

#[test]
fn zero_bonds_are_always_inactive_for_all_parameter_combinations() {
    let input = topology(
        3,
        &[
            (0, 1, BondOrder::Zero),
            (1, 2, BondOrder::Single),
            (2, 0, BondOrder::Single),
        ],
    );
    for dative in [false, true] {
        for hydrogen in [false, true] {
            assert_complete_info(
                find_sssr(&input, &params(dative, hydrogen)).unwrap(),
                RingFindType::Sssr,
                3,
                3,
                &[],
                &[],
            );
        }
    }
}

#[test]
fn every_dative_variant_obeys_only_the_dative_flag() {
    for order in [
        BondOrder::Dative,
        BondOrder::DativeOne,
        BondOrder::DativeLeft,
        BondOrder::DativeRight,
    ] {
        let input = topology(
            3,
            &[
                (0, 1, order),
                (1, 2, BondOrder::Single),
                (2, 0, BondOrder::Single),
            ],
        );
        for hydrogen in [false, true] {
            assert_complete_info(
                find_sssr(&input, &params(false, hydrogen)).unwrap(),
                RingFindType::Sssr,
                3,
                3,
                &[],
                &[],
            );
            assert_complete_info(
                find_sssr(&input, &params(true, hydrogen)).unwrap(),
                RingFindType::Sssr,
                3,
                3,
                &[&[0, 1, 2]],
                &[&[0, 1, 2]],
            );
        }
    }
}

#[test]
fn hydrogen_bonds_obey_only_the_hydrogen_flag_in_both_algorithms() {
    let input = topology(
        3,
        &[
            (0, 1, BondOrder::Hydrogen),
            (1, 2, BondOrder::Single),
            (2, 0, BondOrder::Single),
        ],
    );
    for dative in [false, true] {
        assert_complete_info(
            find_sssr(&input, &params(dative, false)).unwrap(),
            RingFindType::Sssr,
            3,
            3,
            &[],
            &[],
        );
        assert_complete_info(
            symmetrized_sssr(&input, &params(dative, true)).unwrap(),
            RingFindType::SymmSssr,
            3,
            3,
            &[&[0, 1, 2]],
            &[&[0, 1, 2]],
        );
    }
}

#[test]
fn interleaved_dative_hydrogen_and_zero_rows_require_the_exact_active_set() {
    let input = topology(
        5,
        &[
            (0, 1, BondOrder::Dative),
            (3, 4, BondOrder::Zero),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Hydrogen),
            (3, 0, BondOrder::Single),
            (1, 4, BondOrder::Single),
        ],
    );
    for (dative, hydrogen) in [(false, false), (true, false), (false, true)] {
        assert_eq!(
            find_sssr(&input, &params(dative, hydrogen))
                .unwrap()
                .num_rings(),
            0
        );
    }
    assert_complete_info(
        find_sssr(&input, &params(true, true)).unwrap(),
        RingFindType::Sssr,
        5,
        6,
        &[&[0, 3, 2, 1]],
        &[&[4, 3, 2, 0]],
    );
}

#[test]
fn fused_source_case_has_exact_rows_memberships_and_fusion_observations() {
    // Graph of C1CC2C1C2 in source bond creation order.
    let input = topology(
        5,
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Single),
            (3, 0, BondOrder::Single),
            (3, 4, BondOrder::Single),
            (4, 2, BondOrder::Single),
        ],
    );
    assert_complete_info(
        find_sssr(&input, &RingSearchParams::default()).unwrap(),
        RingFindType::Sssr,
        5,
        6,
        &[&[0, 3, 2, 1], &[4, 3, 2]],
        &[&[3, 2, 1, 0], &[4, 2, 5]],
    );
}

#[test]
fn issue_217_distinguishes_sssr_from_three_symmetric_four_rings() {
    // Graph of C=C1C2CC1C2 in source bond creation order.
    let input = topology(
        6,
        &[
            (0, 1, BondOrder::Double),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Single),
            (3, 4, BondOrder::Single),
            (4, 1, BondOrder::Single),
            (4, 5, BondOrder::Single),
            (5, 2, BondOrder::Single),
        ],
    );
    assert_complete_info(
        find_sssr(&input, &RingSearchParams::default()).unwrap(),
        RingFindType::Sssr,
        6,
        7,
        &[&[1, 4, 3, 2], &[1, 4, 5, 2]],
        &[&[4, 3, 2, 1], &[4, 5, 6, 1]],
    );
    assert_complete_info(
        symmetrized_sssr(&input, &RingSearchParams::default()).unwrap(),
        RingFindType::SymmSssr,
        6,
        7,
        &[&[1, 4, 3, 2], &[1, 4, 5, 2], &[3, 4, 5, 2]],
        &[&[4, 3, 2, 1], &[4, 5, 6, 1], &[3, 5, 6, 2]],
    );
}

#[test]
fn cubane_has_five_sssr_rows_and_six_symmetric_rows() {
    let input = topology(
        8,
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Single),
            (3, 0, BondOrder::Single),
            (4, 5, BondOrder::Single),
            (5, 6, BondOrder::Single),
            (6, 7, BondOrder::Single),
            (7, 4, BondOrder::Single),
            (0, 4, BondOrder::Single),
            (1, 5, BondOrder::Single),
            (2, 6, BondOrder::Single),
            (3, 7, BondOrder::Single),
        ],
    );
    let sssr = find_sssr(&input, &RingSearchParams::default()).unwrap();
    assert_eq!(sssr.find_type(), RingFindType::Sssr);
    assert_eq!(sssr.num_rings(), 5);
    assert!(sssr.atom_rings().iter().all(|ring| ring.len() == 4));
    assert!(sssr.bond_rings().iter().all(|ring| ring.len() == 4));

    let symmetric = symmetrized_sssr(&input, &RingSearchParams::default()).unwrap();
    assert_eq!(symmetric.find_type(), RingFindType::SymmSssr);
    assert_eq!(symmetric.num_rings(), 6);
    assert!(symmetric.atom_rings().iter().all(|ring| ring.len() == 4));
    assert!(symmetric.bond_rings().iter().all(|ring| ring.len() == 4));
    assert_eq!(
        &symmetric.atom_rings()[..sssr.num_rings()],
        sssr.atom_rings()
    );
    assert_eq!(
        &symmetric.bond_rings()[..sssr.num_rings()],
        sssr.bond_rings()
    );
}

#[test]
fn figueras_and_highly_fused_source_graphs_keep_their_source_ring_counts() {
    // Graph of C12CC(CC2)CC1 (Figueras figure 4). Its discarded larger
    // candidate must not be appended by symmetrization.
    let figure_four = topology(
        7,
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Single),
            (3, 4, BondOrder::Single),
            (4, 0, BondOrder::Single),
            (2, 5, BondOrder::Single),
            (5, 6, BondOrder::Single),
            (6, 0, BondOrder::Single),
        ],
    );
    let figure_sssr = find_sssr(&figure_four, &RingSearchParams::default()).unwrap();
    assert_eq!(figure_sssr.find_type(), RingFindType::Sssr);
    assert_eq!(figure_sssr.num_rings(), 2);
    assert_eq!(
        symmetrized_sssr(&figure_four, &RingSearchParams::default())
            .unwrap()
            .num_rings(),
        2
    );

    // Graph of C1CC2C1CCC2.
    let mixed_sizes = topology(
        7,
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Single),
            (3, 0, BondOrder::Single),
            (3, 4, BondOrder::Single),
            (4, 5, BondOrder::Single),
            (5, 6, BondOrder::Single),
            (6, 2, BondOrder::Single),
        ],
    );
    let mixed = find_sssr(&mixed_sizes, &RingSearchParams::default()).unwrap();
    assert_eq!(mixed.find_type(), RingFindType::Sssr);
    assert_eq!(
        mixed.atom_rings().iter().map(Vec::len).collect::<Vec<_>>(),
        vec![4, 5]
    );

    // Counterexample from the source ring-perception test:
    // C123C4C5C6(C3)C7C1C8C2C4C5C6C78.
    let counterexample = topology(
        13,
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Single),
            (3, 4, BondOrder::Single),
            (4, 0, BondOrder::Single),
            (3, 5, BondOrder::Single),
            (5, 6, BondOrder::Single),
            (6, 0, BondOrder::Single),
            (6, 7, BondOrder::Single),
            (7, 8, BondOrder::Single),
            (8, 0, BondOrder::Single),
            (8, 9, BondOrder::Single),
            (9, 1, BondOrder::Single),
            (9, 10, BondOrder::Single),
            (10, 2, BondOrder::Single),
            (10, 11, BondOrder::Single),
            (11, 3, BondOrder::Single),
            (11, 12, BondOrder::Single),
            (12, 5, BondOrder::Single),
            (12, 7, BondOrder::Single),
        ],
    );
    let counter_sssr = find_sssr(&counterexample, &RingSearchParams::default()).unwrap();
    assert_eq!(counter_sssr.find_type(), RingFindType::Sssr);
    assert_eq!(counter_sssr.num_rings(), 8 - 1);
    assert!(counter_sssr.atom_rings().iter().all(|ring| ring.len() < 6));
    let counter_symmetric =
        symmetrized_sssr(&counterexample, &RingSearchParams::default()).unwrap();
    assert_eq!(counter_symmetric.num_rings(), 8);
    assert!(
        counter_symmetric
            .atom_rings()
            .iter()
            .all(|ring| ring.len() < 6)
    );
}

#[test]
fn validated_and_parts_entrypoints_are_equal_repeatable_and_non_mutating() {
    let input = topology(
        4,
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Single),
            (3, 0, BondOrder::Single),
        ],
    );
    let before = input.clone();
    let first = find_sssr(&input, &RingSearchParams::default()).unwrap();
    let second = find_sssr(&input, &RingSearchParams::default()).unwrap();
    let default_parts =
        find_sssr_from_parts(input.atoms.len(), &input.bonds, &input.adjacency).unwrap();
    let option_parts = find_sssr_with_options_from_parts(
        input.atoms.len(),
        &input.bonds,
        &input.adjacency,
        false,
        false,
    )
    .unwrap();
    assert_eq!(first, second);
    assert_eq!(first, default_parts);
    assert_eq!(first, option_parts);
    assert_eq!(input, before);
}

#[test]
fn invalid_topology_is_rejected_before_either_algorithm_traverses() {
    let atom_id_mismatch = TopologyBlock {
        atoms: vec![Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C))],
        adjacency: AdjacencyList::from_topology(1, &[]),
        ..TopologyBlock::default()
    };
    let expected_atom = Err(RingFindingError::InvalidTopology(
        TopologyValidationError::AtomIdMismatch {
            position: 0,
            id: AtomId::new(1),
        },
    ));
    assert_eq!(
        find_sssr(&atom_id_mismatch, &RingSearchParams::default()),
        expected_atom
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
        symmetrized_sssr(&stale_adjacency, &RingSearchParams::default()),
        Err(RingFindingError::InvalidTopology(
            TopologyValidationError::AdjacencyMismatch
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
        find_sssr(&invalid_stereo_group, &RingSearchParams::default()),
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
        symmetrized_sssr(&invalid_sgroup, &RingSearchParams::default()),
        Err(RingFindingError::InvalidTopology(
            TopologyValidationError::SubstanceGroupAtomOutOfRange {
                sgroup: SubstanceGroupId::new(0),
                atom: AtomId::new(5),
                atom_count: 1,
            }
        ))
    );
}
