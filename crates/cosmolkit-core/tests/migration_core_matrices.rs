use cosmolkit_core::{
    AdjacencyMatrixParams, DistanceMatrix3dParams, GraphPath, MatrixError, PathRepresentation,
    PathSearchParams, TopologicalDistanceMatrixParams, adjacency_matrix, all_paths_of_length,
    distance_matrix_3d, topological_distance_matrix,
};
use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, Conformer2D, Conformer3D,
    CoordinateBlock, CoordinateValidationError, TopologyBlock, TopologyValidationError,
};
use cosmolkit_types::{BondOrder, Element};

const LOCAL_INF: f64 = 100_000_000.0;

fn atom(index: usize, element: Element) -> Atom {
    Atom::from_spec(AtomId::new(index), AtomSpec::new(element))
}

fn topology_from_specs(elements: &[Element], bond_specs: Vec<BondSpec>) -> TopologyBlock {
    let atoms = elements
        .iter()
        .copied()
        .enumerate()
        .map(|(index, element)| atom(index, element))
        .collect();
    let bonds = bond_specs
        .into_iter()
        .enumerate()
        .map(|(index, spec)| Bond::from_spec(BondId::new(index), spec))
        .collect();
    TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
}

fn topology(elements: &[Element], edges: &[(usize, usize, BondOrder)]) -> TopologyBlock {
    topology_from_specs(
        elements,
        edges
            .iter()
            .map(|&(begin, end, order)| BondSpec::new(AtomId::new(begin), AtomId::new(end), order))
            .collect(),
    )
}

fn carbon_topology(atom_count: usize, edges: &[(usize, usize, BondOrder)]) -> TopologyBlock {
    topology(&vec![Element::C; atom_count], edges)
}

fn atom_ids(values: &[usize]) -> Vec<AtomId> {
    values.iter().copied().map(AtomId::new).collect()
}

fn bond_ids(values: &[usize]) -> Vec<BondId> {
    values.iter().copied().map(BondId::new).collect()
}

#[test]
fn empty_singleton_and_dense_accessors_preserve_row_major_shape() {
    let empty = TopologyBlock::default();
    for matrix in [
        adjacency_matrix(&empty, &AdjacencyMatrixParams::default()).unwrap(),
        topological_distance_matrix(&empty, &TopologicalDistanceMatrixParams::default()).unwrap(),
    ] {
        assert_eq!(matrix.dimension(), 0);
        assert!(matrix.values().is_empty());
        assert_eq!(matrix.get(0, 0), None);
    }

    let singleton = carbon_topology(1, &[]);
    let matrix =
        topological_distance_matrix(&singleton, &TopologicalDistanceMatrixParams::default())
            .unwrap();
    assert_eq!(matrix.dimension(), 1);
    assert_eq!(matrix.values(), &[0.0]);
    assert_eq!(matrix.get(0, 0), Some(0.0));
    assert_eq!(matrix.get(0, 1), None);
    assert_eq!(matrix.get(1, 0), None);
}

#[test]
fn adjacency_locks_exact_rows_filtering_and_source_byte_fill() {
    let input = carbon_topology(
        4,
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Double),
            (2, 3, BondOrder::Triple),
        ],
    );
    let full = adjacency_matrix(&input, &AdjacencyMatrixParams::default()).unwrap();
    assert_eq!(
        full.values(),
        &[
            0.0, 1.0, 0.0, 0.0, // row 0
            1.0, 0.0, 1.0, 0.0, // row 1
            0.0, 1.0, 0.0, 1.0, // row 2
            0.0, 0.0, 1.0, 0.0, // row 3
        ]
    );

    let byte_filled = f64::from_ne_bytes([1_u8; size_of::<f64>()]);
    let filtered = adjacency_matrix(
        &input,
        &AdjacencyMatrixParams {
            empty_value: 1,
            bonds_to_use: Some(bond_ids(&[2, 0])),
            ..AdjacencyMatrixParams::default()
        },
    )
    .unwrap();
    assert_eq!(
        filtered.values(),
        &[
            byte_filled,
            1.0,
            byte_filled,
            byte_filled,
            1.0,
            byte_filled,
            byte_filled,
            byte_filled,
            byte_filled,
            byte_filled,
            byte_filled,
            1.0,
            byte_filled,
            byte_filled,
            1.0,
            byte_filled,
        ]
    );
}

#[test]
fn weighted_adjacency_covers_standard_aromatic_dative_and_hydrogen_orders() {
    let aromatic =
        BondSpec::new(AtomId::new(3), AtomId::new(4), BondOrder::Aromatic).with_aromatic(true);
    let input = topology_from_specs(
        &[Element::C; 8],
        vec![
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Double),
            BondSpec::new(AtomId::new(2), AtomId::new(3), BondOrder::Triple),
            aromatic,
            BondSpec::new(AtomId::new(4), AtomId::new(5), BondOrder::Dative),
            BondSpec::new(AtomId::new(5), AtomId::new(6), BondOrder::DativeOne),
            BondSpec::new(AtomId::new(6), AtomId::new(7), BondOrder::Hydrogen),
        ],
    );
    let matrix = adjacency_matrix(
        &input,
        &AdjacencyMatrixParams {
            use_bond_order: true,
            ..AdjacencyMatrixParams::default()
        },
    )
    .unwrap();
    assert_eq!(matrix.get(0, 1), Some(1.0));
    assert_eq!(matrix.get(1, 0), Some(1.0));
    assert_eq!(matrix.get(1, 2), Some(2.0));
    assert_eq!(matrix.get(2, 1), Some(2.0));
    assert_eq!(matrix.get(2, 3), Some(3.0));
    assert_eq!(matrix.get(3, 2), Some(3.0));
    assert_eq!(matrix.get(3, 4), Some(1.5));
    assert_eq!(matrix.get(4, 3), Some(1.5));
    assert_eq!(matrix.get(4, 5), Some(0.0));
    assert_eq!(matrix.get(5, 4), Some(1.0));
    assert_eq!(matrix.get(5, 6), Some(0.0));
    assert_eq!(matrix.get(6, 5), Some(1.0));
    assert_eq!(matrix.get(6, 7), Some(0.0));
    assert_eq!(matrix.get(7, 6), Some(0.0));
}

#[test]
fn weighted_matrix_rejects_unsupported_orders_and_accepts_source_zero_orders() {
    for order in [
        BondOrder::Unspecified,
        BondOrder::Zero,
        BondOrder::Ionic,
        BondOrder::Hydrogen,
    ] {
        let input = carbon_topology(2, &[(0, 1, order)]);
        let matrix = topological_distance_matrix(
            &input,
            &TopologicalDistanceMatrixParams {
                use_bond_order: true,
                ..TopologicalDistanceMatrixParams::default()
            },
        )
        .unwrap();
        assert!(matrix.get(0, 1).unwrap().is_infinite(), "order={order:?}");
        assert!(matrix.get(1, 0).unwrap().is_infinite(), "order={order:?}");
    }

    for order in [BondOrder::Other, BondOrder::ThreeCenter] {
        let input = carbon_topology(2, &[(0, 1, order)]);
        let expected = Err(MatrixError::UnsupportedBondOrder {
            bond: BondId::new(0),
            order,
        });
        assert_eq!(
            topological_distance_matrix(
                &input,
                &TopologicalDistanceMatrixParams {
                    use_bond_order: true,
                    ..TopologicalDistanceMatrixParams::default()
                }
            ),
            expected
        );
        assert_eq!(
            adjacency_matrix(
                &input,
                &AdjacencyMatrixParams {
                    use_bond_order: true,
                    ..AdjacencyMatrixParams::default()
                }
            ),
            expected
        );
    }
}

#[test]
fn topological_distance_rows_cover_connected_cycle_and_disconnected_graphs() {
    let chain = carbon_topology(3, &[(0, 1, BondOrder::Single), (1, 2, BondOrder::Single)]);
    assert_eq!(
        topological_distance_matrix(&chain, &TopologicalDistanceMatrixParams::default())
            .unwrap()
            .values(),
        &[0.0, 1.0, 2.0, 1.0, 0.0, 1.0, 2.0, 1.0, 0.0]
    );

    let cycle = carbon_topology(
        4,
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Single),
            (3, 0, BondOrder::Single),
        ],
    );
    assert_eq!(
        topological_distance_matrix(&cycle, &TopologicalDistanceMatrixParams::default())
            .unwrap()
            .values(),
        &[
            0.0, 1.0, 2.0, 1.0, 1.0, 0.0, 1.0, 2.0, 2.0, 1.0, 0.0, 1.0, 1.0, 2.0, 1.0, 0.0,
        ]
    );

    let disconnected = carbon_topology(3, &[(0, 1, BondOrder::Single)]);
    assert_eq!(
        topological_distance_matrix(&disconnected, &TopologicalDistanceMatrixParams::default())
            .unwrap()
            .values(),
        &[
            0.0, 1.0, LOCAL_INF, 1.0, 0.0, LOCAL_INF, LOCAL_INF, LOCAL_INF, 0.0,
        ]
    );
}

#[test]
fn topological_bond_order_and_atom_weight_branches_match_source_rows() {
    let input = topology_from_specs(
        &[Element::C, Element::C, Element::O, Element::DUMMY],
        vec![
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Double),
            BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Aromatic).with_aromatic(true),
            BondSpec::new(AtomId::new(2), AtomId::new(3), BondOrder::Single),
        ],
    );
    let matrix = topological_distance_matrix(
        &input,
        &TopologicalDistanceMatrixParams {
            use_bond_order: true,
            use_atom_weights: true,
            ..TopologicalDistanceMatrixParams::default()
        },
    )
    .unwrap();
    let aromatic = 2.0 / 3.0;
    assert_eq!(matrix.get(0, 0), Some(1.0));
    assert_eq!(matrix.get(1, 1), Some(1.0));
    assert_eq!(matrix.get(2, 2), Some(0.75));
    assert!(matrix.get(3, 3).unwrap().is_infinite());
    assert_eq!(matrix.get(0, 1), Some(0.5));
    assert_eq!(matrix.get(1, 2), Some(aromatic));
    assert_eq!(matrix.get(0, 2), Some(0.5 + aromatic));
    assert_eq!(matrix.get(0, 3), Some(1.5 + aromatic));
}

#[test]
fn active_subsets_preserve_order_and_inferred_or_explicit_bond_selection() {
    let input = topology(
        &[Element::C, Element::O, Element::N, Element::C],
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Double),
            (2, 3, BondOrder::Single),
            (0, 3, BondOrder::Single),
        ],
    );
    let inferred = topological_distance_matrix(
        &input,
        &TopologicalDistanceMatrixParams {
            use_atom_weights: true,
            active_atoms: Some(atom_ids(&[2, 1, 0])),
            ..TopologicalDistanceMatrixParams::default()
        },
    )
    .unwrap();
    assert_eq!(inferred.dimension(), 3);
    assert_eq!(
        inferred.values(),
        &[6.0 / 7.0, 1.0, 2.0, 1.0, 0.75, 1.0, 2.0, 1.0, 1.0,]
    );

    let explicit = topological_distance_matrix(
        &input,
        &TopologicalDistanceMatrixParams {
            active_atoms: Some(atom_ids(&[0, 1, 2])),
            active_bonds: Some(bond_ids(&[1])),
            ..TopologicalDistanceMatrixParams::default()
        },
    )
    .unwrap();
    assert_eq!(
        explicit.values(),
        &[
            0.0, LOCAL_INF, LOCAL_INF, LOCAL_INF, 0.0, 1.0, LOCAL_INF, 1.0, 0.0,
        ]
    );

    let none = topological_distance_matrix(
        &input,
        &TopologicalDistanceMatrixParams {
            active_atoms: Some(vec![]),
            ..TopologicalDistanceMatrixParams::default()
        },
    )
    .unwrap();
    assert_eq!(none.dimension(), 0);
    assert!(none.values().is_empty());
}

#[test]
fn active_subset_validation_reports_every_position_id_and_endpoint_field() {
    let input = carbon_topology(3, &[(0, 1, BondOrder::Single), (1, 2, BondOrder::Single)]);
    let run = |active_atoms, active_bonds| {
        topological_distance_matrix(
            &input,
            &TopologicalDistanceMatrixParams {
                active_atoms,
                active_bonds,
                ..TopologicalDistanceMatrixParams::default()
            },
        )
    };
    assert_eq!(
        run(None, Some(bond_ids(&[0]))),
        Err(MatrixError::ActiveBondsWithoutAtoms)
    );
    assert_eq!(
        run(Some(atom_ids(&[0, 9])), None),
        Err(MatrixError::ActiveAtomOutOfRange {
            position: 1,
            atom: AtomId::new(9),
            atom_count: 3,
        })
    );
    assert_eq!(
        run(Some(atom_ids(&[2, 0, 2])), None),
        Err(MatrixError::DuplicateActiveAtom {
            atom: AtomId::new(2),
            first_position: 0,
            second_position: 2,
        })
    );
    assert_eq!(
        run(Some(atom_ids(&[0, 1, 2])), Some(bond_ids(&[7]))),
        Err(MatrixError::ActiveBondOutOfRange {
            position: 0,
            bond: BondId::new(7),
            bond_count: 2,
        })
    );
    assert_eq!(
        run(Some(atom_ids(&[0, 1, 2])), Some(bond_ids(&[1, 0, 1]))),
        Err(MatrixError::DuplicateActiveBond {
            bond: BondId::new(1),
            first_position: 0,
            second_position: 2,
        })
    );
    assert_eq!(
        run(Some(atom_ids(&[1])), Some(bond_ids(&[0]))),
        Err(MatrixError::ActiveBondEndpointMissing {
            bond: BondId::new(0),
            endpoint: "begin",
            atom: AtomId::new(0),
        })
    );
    assert_eq!(
        run(Some(atom_ids(&[0])), Some(bond_ids(&[0]))),
        Err(MatrixError::ActiveBondEndpointMissing {
            bond: BondId::new(0),
            endpoint: "end",
            atom: AtomId::new(1),
        })
    );

    for params in [
        AdjacencyMatrixParams {
            bonds_to_use: Some(bond_ids(&[9])),
            ..AdjacencyMatrixParams::default()
        },
        AdjacencyMatrixParams {
            bonds_to_use: Some(bond_ids(&[1, 0, 1])),
            ..AdjacencyMatrixParams::default()
        },
    ] {
        assert!(matches!(
            adjacency_matrix(&input, &params),
            Err(MatrixError::ActiveBondOutOfRange { .. })
                | Err(MatrixError::DuplicateActiveBond { .. })
        ));
    }
}

#[test]
fn all_matrix_entries_reject_invalid_topology_before_indexing() {
    let valid = carbon_topology(2, &[(0, 1, BondOrder::Single)]);
    let invalid = TopologyBlock {
        adjacency: AdjacencyList::from_topology(2, &[]),
        ..valid
    };
    let expected = Err(MatrixError::InvalidTopology(
        TopologyValidationError::AdjacencyMismatch,
    ));
    assert_eq!(
        adjacency_matrix(&invalid, &AdjacencyMatrixParams::default()),
        expected
    );
    assert_eq!(
        topological_distance_matrix(&invalid, &TopologicalDistanceMatrixParams::default()),
        expected
    );
    assert_eq!(
        distance_matrix_3d(
            &invalid,
            &CoordinateBlock::default(),
            &DistanceMatrix3dParams::default()
        ),
        expected
    );
}

#[test]
fn three_dimensional_selection_preserves_ids_rows_symmetry_and_weights() {
    let topology = topology(&[Element::C, Element::O, Element::DUMMY], &[]);
    let coordinates = CoordinateBlock {
        conformers_3d: vec![
            Conformer3D::new(
                8,
                vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [3.0, 4.0, 0.0]],
                true,
            ),
            Conformer3D::new(
                3,
                vec![[0.0, 0.0, 0.0], [0.0, 0.0, 12.0], [0.0, 5.0, 12.0]],
                true,
            ),
        ],
        ..CoordinateBlock::default()
    };
    let first =
        distance_matrix_3d(&topology, &coordinates, &DistanceMatrix3dParams::default()).unwrap();
    assert_eq!(
        first.values(),
        &[0.0, 3.0, 5.0, 3.0, 0.0, 4.0, 5.0, 4.0, 0.0]
    );

    let selected = distance_matrix_3d(
        &topology,
        &coordinates,
        &DistanceMatrix3dParams {
            conformer_id: Some(3),
            use_atom_weights: true,
        },
    )
    .unwrap();
    assert_eq!(selected.get(0, 0), Some(1.0));
    assert_eq!(selected.get(1, 1), Some(0.75));
    assert!(selected.get(2, 2).unwrap().is_infinite());
    assert_eq!(selected.get(0, 1), Some(12.0));
    assert_eq!(selected.get(1, 0), Some(12.0));
    assert_eq!(selected.get(1, 2), Some(5.0));
    assert_eq!(selected.get(2, 1), Some(5.0));
    assert_eq!(selected.get(0, 2), Some(13.0));
    assert_eq!(selected.get(2, 0), Some(13.0));
}

#[test]
fn three_dimensional_errors_preserve_validation_and_selection_order() {
    let topology = carbon_topology(2, &[]);
    assert_eq!(
        distance_matrix_3d(
            &topology,
            &CoordinateBlock::default(),
            &DistanceMatrix3dParams::default()
        ),
        Err(MatrixError::No3dConformer)
    );
    let two_d_only = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(4, vec![[0.0, 0.0], [1.0, 0.0]])],
        ..CoordinateBlock::default()
    };
    assert_eq!(
        distance_matrix_3d(&topology, &two_d_only, &DistanceMatrix3dParams::default()),
        Err(MatrixError::No3dConformer)
    );

    let valid = CoordinateBlock {
        conformers_3d: vec![Conformer3D::new(
            11,
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            true,
        )],
        ..CoordinateBlock::default()
    };
    assert_eq!(
        distance_matrix_3d(
            &topology,
            &valid,
            &DistanceMatrix3dParams {
                conformer_id: Some(0),
                ..DistanceMatrix3dParams::default()
            }
        ),
        Err(MatrixError::ConformerNotFound { conformer_id: 0 })
    );

    let short = CoordinateBlock {
        conformers_3d: vec![Conformer3D::new(91, vec![[0.0, 0.0, 0.0]], true)],
        ..CoordinateBlock::default()
    };
    assert_eq!(
        distance_matrix_3d(
            &topology,
            &short,
            &DistanceMatrix3dParams {
                conformer_id: Some(404),
                ..DistanceMatrix3dParams::default()
            }
        ),
        Err(MatrixError::InvalidCoordinates(
            CoordinateValidationError::RowCount {
                dimension: "3D",
                conformer: 91,
                rows: 1,
                atom_count: 2,
            }
        ))
    );

    let non_finite = CoordinateBlock {
        conformers_3d: vec![Conformer3D::new(
            6,
            vec![[0.0, 0.0, 0.0], [1.0, f64::NAN, 0.0]],
            true,
        )],
        ..CoordinateBlock::default()
    };
    assert_eq!(
        distance_matrix_3d(&topology, &non_finite, &DistanceMatrix3dParams::default()),
        Err(MatrixError::InvalidCoordinates(
            CoordinateValidationError::NonFiniteCoordinate {
                dimension: "3D",
                conformer: 6,
                atom: 1,
                axis: "y",
            }
        ))
    );
}

#[test]
fn shortest_only_paths_regress_through_the_unique_unweighted_matrix_owner() {
    let input = carbon_topology(
        4,
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Single),
        ],
    );
    let paths = all_paths_of_length(
        &input,
        2,
        &PathSearchParams {
            use_hydrogens: false,
            rooted_at_atom: None,
            only_shortest_paths: true,
            representation: PathRepresentation::Bonds,
        },
    )
    .unwrap();
    assert_eq!(
        paths,
        vec![
            GraphPath::Bonds(bond_ids(&[0, 1])),
            GraphPath::Bonds(bond_ids(&[1, 2])),
        ]
    );
}
