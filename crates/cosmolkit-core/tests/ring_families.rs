use cosmolkit_core::{
    RingFindType, RingFindingError, RingInfo, RingSearchParams, find_ring_families,
    find_ring_families_from_parts,
};
use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, TopologyBlock,
    TopologyValidationError,
};
use cosmolkit_ringdecomposer::RingDecomposerError;
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

fn atom_rows(rows: &[&[usize]]) -> Vec<Vec<AtomId>> {
    rows.iter()
        .map(|row| row.iter().copied().map(AtomId::new).collect())
        .collect()
}

fn bond_rows(rows: &[&[usize]]) -> Vec<Vec<BondId>> {
    rows.iter()
        .map(|row| row.iter().copied().map(BondId::new).collect())
        .collect()
}

fn params(dative: bool, hydrogen: bool) -> RingSearchParams {
    RingSearchParams {
        include_dative_bonds: dative,
        include_hydrogen_bonds: hydrogen,
    }
}

fn assert_family_info(
    info: &RingInfo,
    expected_atoms: &[&[usize]],
    expected_bonds: &[&[usize]],
    relevant_cycles: usize,
) {
    assert!(info.is_initialized());
    assert!(info.are_ring_families_initialized());
    assert_eq!(info.find_type(), RingFindType::OtherOrUnknown);
    assert!(!info.is_find_fast_or_better());
    assert!(!info.is_sssr_or_better());
    assert!(!info.is_symm_sssr());
    assert_eq!(info.num_rings(), 0);
    assert!(info.atom_rings().is_empty());
    assert!(info.bond_rings().is_empty());
    assert_eq!(info.num_ring_families(), expected_atoms.len());
    assert_eq!(info.atom_ring_families(), atom_rows(expected_atoms));
    assert_eq!(info.bond_ring_families(), bond_rows(expected_bonds));
    assert_eq!(info.num_relevant_cycles(), Ok(relevant_cycles));
}

#[test]
fn defaults_empty_and_nonempty_acyclic_states_are_distinct() {
    assert_eq!(RingSearchParams::default(), params(false, false));
    assert_eq!(
        find_ring_families(&topology(0, &[]), &RingSearchParams::default()),
        Err(RingFindingError::RingDecomposer(
            RingDecomposerError::EmptyGraph
        ))
    );

    for input in [
        topology(1, &[]),
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
        assert_family_info(
            &find_ring_families(&input, &RingSearchParams::default()).unwrap(),
            &[],
            &[],
            0,
        );
    }

    let plain = RingInfo::new(RingFindType::OtherOrUnknown, 1, 0);
    assert!(!plain.are_ring_families_initialized());
    assert_eq!(
        plain.num_relevant_cycles(),
        Err(RingFindingError::UnsupportedBranch {
            reason: "URF relevant-cycle data is not modeled"
        })
    );
}

#[test]
fn source_primary_polycycle_has_all_five_exact_family_rows() {
    // c1ccc2c(c1)C3CC3C4CC5CC4CC25 in RDKit atom/bond creation order.
    let input = topology(
        16,
        &[
            (0, 1, BondOrder::Aromatic),
            (1, 2, BondOrder::Aromatic),
            (2, 3, BondOrder::Aromatic),
            (3, 4, BondOrder::Aromatic),
            (4, 5, BondOrder::Aromatic),
            (4, 6, BondOrder::Single),
            (6, 7, BondOrder::Single),
            (7, 8, BondOrder::Single),
            (8, 9, BondOrder::Single),
            (9, 10, BondOrder::Single),
            (10, 11, BondOrder::Single),
            (11, 12, BondOrder::Single),
            (12, 13, BondOrder::Single),
            (13, 14, BondOrder::Single),
            (14, 15, BondOrder::Single),
            (5, 0, BondOrder::Aromatic),
            (15, 3, BondOrder::Single),
            (8, 6, BondOrder::Single),
            (13, 9, BondOrder::Single),
            (15, 11, BondOrder::Single),
        ],
    );
    let info = find_ring_families(&input, &RingSearchParams::default()).unwrap();
    assert_family_info(
        &info,
        &[
            &[6, 7, 8],
            &[9, 10, 11, 12, 13],
            &[11, 12, 13, 14, 15],
            &[0, 1, 2, 3, 4, 5],
            &[3, 4, 6, 8, 9, 10, 11, 13, 14, 15],
        ],
        &[
            &[6, 7, 17],
            &[9, 10, 11, 12, 18],
            &[11, 12, 13, 14, 19],
            &[0, 1, 2, 3, 4, 15],
            &[3, 5, 8, 9, 10, 13, 14, 16, 17, 18, 19],
        ],
        info.num_relevant_cycles().unwrap(),
    );
    assert_eq!(info.atom_ring_families()[4].len(), 10);
    assert_eq!(info.bond_ring_families()[4].len(), 11);
}

#[test]
fn simple_fused_spiro_and_disconnected_cycles_preserve_urf_order() {
    let triangle = topology(
        3,
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (2, 0, BondOrder::Single),
        ],
    );
    assert_family_info(
        &find_ring_families(&triangle, &params(false, false)).unwrap(),
        &[&[0, 1, 2]],
        &[&[0, 1, 2]],
        1,
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
    assert_family_info(
        &find_ring_families(&fused, &params(false, false)).unwrap(),
        &[&[0, 1, 2], &[1, 2, 3]],
        &[&[0, 1, 2], &[1, 3, 4]],
        2,
    );

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
    assert_family_info(
        &find_ring_families(&spiro, &params(false, false)).unwrap(),
        &[&[0, 1, 2], &[0, 3, 4]],
        &[&[0, 1, 2], &[3, 4, 5]],
        2,
    );

    let disconnected = topology(
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
    assert_family_info(
        &find_ring_families(&disconnected, &params(false, false)).unwrap(),
        &[&[0, 1, 2], &[3, 4, 5]],
        &[&[1, 3, 5], &[0, 2, 4]],
        2,
    );
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
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 0, order),
            ],
        );
        assert_family_info(
            &find_ring_families(&input, &params(false, false)).unwrap(),
            &[],
            &[],
            0,
        );
        assert_family_info(
            &find_ring_families(&input, &params(true, false)).unwrap(),
            &[&[0, 1, 2]],
            &[&[0, 1, 2]],
            1,
        );
        assert_family_info(
            &find_ring_families(&input, &params(false, true)).unwrap(),
            &[],
            &[],
            0,
        );
    }
}

#[test]
fn hydrogen_and_zero_bonds_follow_independent_source_policies() {
    let hydrogen = topology(
        3,
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (2, 0, BondOrder::Hydrogen),
        ],
    );
    assert_family_info(
        &find_ring_families(&hydrogen, &params(false, false)).unwrap(),
        &[],
        &[],
        0,
    );
    assert_family_info(
        &find_ring_families(&hydrogen, &params(false, true)).unwrap(),
        &[&[0, 1, 2]],
        &[&[0, 1, 2]],
        1,
    );

    let zero = topology(
        3,
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (2, 0, BondOrder::Zero),
        ],
    );
    for options in [
        params(false, false),
        params(true, false),
        params(false, true),
        params(true, true),
    ] {
        assert_family_info(&find_ring_families(&zero, &options).unwrap(), &[], &[], 0);
    }
}

#[test]
fn four_option_combinations_keep_interleaved_original_bond_ids() {
    let input = topology(
        6,
        &[
            (0, 1, BondOrder::Single),
            (3, 4, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (4, 5, BondOrder::Single),
            (2, 0, BondOrder::Dative),
            (5, 3, BondOrder::Hydrogen),
        ],
    );
    for (options, expected_atoms, expected_bonds, count) in [
        (params(false, false), vec![], vec![], 0),
        (
            params(true, false),
            vec![vec![0, 1, 2]],
            vec![vec![0, 2, 4]],
            1,
        ),
        (
            params(false, true),
            vec![vec![3, 4, 5]],
            vec![vec![1, 3, 5]],
            1,
        ),
        (
            params(true, true),
            vec![vec![0, 1, 2], vec![3, 4, 5]],
            vec![vec![0, 2, 4], vec![1, 3, 5]],
            2,
        ),
    ] {
        let atom_slices = expected_atoms.iter().map(Vec::as_slice).collect::<Vec<_>>();
        let bond_slices = expected_bonds.iter().map(Vec::as_slice).collect::<Vec<_>>();
        assert_family_info(
            &find_ring_families(&input, &options).unwrap(),
            &atom_slices,
            &bond_slices,
            count,
        );
    }
}

#[test]
fn canonical_and_parts_paths_are_equal_repeatable_and_nonmutating() {
    let input = topology(
        3,
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (2, 0, BondOrder::Single),
        ],
    );
    let snapshot = input.clone();
    let canonical = find_ring_families(&input, &RingSearchParams::default()).unwrap();
    let repeated = find_ring_families(&input, &RingSearchParams::default()).unwrap();
    let parts =
        find_ring_families_from_parts(input.atoms.len(), &input.bonds, false, false).unwrap();
    assert_eq!(canonical, repeated);
    assert_eq!(canonical, parts);
    assert_eq!(input, snapshot);
}

#[test]
fn canonical_validation_precedes_graph_construction() {
    let mut input = topology(
        3,
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (2, 0, BondOrder::Single),
        ],
    );
    input.adjacency = AdjacencyList::from_topology(3, &[]);
    assert_eq!(
        find_ring_families(&input, &RingSearchParams::default()),
        Err(RingFindingError::InvalidTopology(
            TopologyValidationError::AdjacencyMismatch
        ))
    );
}

#[test]
fn parts_projection_reports_original_bond_id_out_of_range_with_fields() {
    let bonds = vec![
        Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        ),
        Bond::from_spec(
            BondId::new(1),
            BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Single),
        ),
        Bond::from_spec(
            BondId::new(7),
            BondSpec::new(AtomId::new(2), AtomId::new(0), BondOrder::Single),
        ),
    ];
    assert_eq!(
        find_ring_families_from_parts(3, &bonds, false, false),
        Err(RingFindingError::RingBondOutOfRange {
            bond: 7,
            bond_count: 3
        })
    );
}

#[test]
fn ring_info_family_insertion_accepts_unequal_rows_and_rejects_out_of_range_fields() {
    let mut unequal = RingInfo::new(RingFindType::OtherOrUnknown, 3, 3);
    assert_eq!(unequal.add_ring_family(&[0, 1], &[0]), Ok(1));
    assert_eq!(unequal.atom_ring_families(), atom_rows(&[&[0, 1]]));
    assert_eq!(unequal.bond_ring_families(), bond_rows(&[&[0]]));

    let mut bad_atom = RingInfo::new(RingFindType::OtherOrUnknown, 3, 3);
    assert_eq!(
        bad_atom.add_ring_family(&[0, 1, 4], &[0, 1, 2]),
        Err(RingFindingError::RingAtomOutOfRange {
            atom: 4,
            atom_count: 3
        })
    );

    let mut bad_bond = RingInfo::new(RingFindType::OtherOrUnknown, 3, 3);
    assert_eq!(
        bad_bond.add_ring_family(&[0, 1, 2], &[0, 1, 5]),
        Err(RingFindingError::RingBondOutOfRange {
            bond: 5,
            bond_count: 3
        })
    );
}
