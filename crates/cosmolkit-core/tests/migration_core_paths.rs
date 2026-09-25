use cosmolkit_core::{
    AtomEnvironmentParams, DetachedPathSubgraph, GraphPath, PathError, PathRepresentation,
    PathSearchParams, SubgraphSearchParams, SubtopologyParams, UniqueSubgraphParams,
    all_paths_in_range, all_paths_of_length, all_subgraphs_in_range, all_subgraphs_of_length,
    atom_environment, bond_ids_from_atom_path, connected_components, shortest_path,
    subtopology_from_path, unique_subgraphs_of_length,
};
use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomQueryPredicate, AtomSpec, Bond, BondId, BondQueryPredicate,
    BondSpec, BondStereo, QueryNode, StereoGroup, StereoGroupKind, SubstanceGroup,
    SubstanceGroupId, SubstanceGroupKind, TopologyBlock, TopologyValidationError,
};
use cosmolkit_types::{BondOrder, Element};

fn atom(index: usize, element: Element) -> Atom {
    Atom::from_spec(AtomId::new(index), AtomSpec::new(element))
}

fn topology(elements: &[Element], edges: &[(usize, usize, BondOrder)]) -> TopologyBlock {
    let atoms = elements
        .iter()
        .copied()
        .enumerate()
        .map(|(index, element)| atom(index, element))
        .collect();
    let bonds = edges
        .iter()
        .copied()
        .enumerate()
        .map(|(index, (begin, end, order))| {
            Bond::from_spec(
                BondId::new(index),
                BondSpec::new(AtomId::new(begin), AtomId::new(end), order),
            )
        })
        .collect();
    TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
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
fn connected_components_preserve_source_atom_order_and_isolates() {
    let empty = carbon_topology(0, &[]);
    assert_eq!(
        connected_components(&empty).unwrap(),
        cosmolkit_core::ConnectedComponents {
            atom_to_component: vec![],
            components: vec![],
        }
    );

    let input = carbon_topology(
        7,
        &[
            (4, 5, BondOrder::Double),
            (1, 2, BondOrder::Zero),
            (5, 3, BondOrder::Aromatic),
        ],
    );
    let result = connected_components(&input).unwrap();
    assert_eq!(result.atom_to_component, vec![0, 1, 1, 2, 2, 2, 3]);
    assert_eq!(
        result.components,
        vec![
            atom_ids(&[0]),
            atom_ids(&[1, 2]),
            atom_ids(&[3, 4, 5]),
            atom_ids(&[6]),
        ]
    );
}

#[test]
fn shortest_path_preserves_adjacency_ties_and_reports_all_endpoint_errors() {
    let input = carbon_topology(
        5,
        &[
            (0, 2, BondOrder::Single),
            (0, 1, BondOrder::Single),
            (2, 3, BondOrder::Single),
            (1, 3, BondOrder::Single),
        ],
    );
    assert_eq!(
        shortest_path(&input, AtomId::new(0), AtomId::new(3)).unwrap(),
        atom_ids(&[0, 2, 3])
    );
    assert!(
        shortest_path(&input, AtomId::new(0), AtomId::new(4))
            .unwrap()
            .is_empty()
    );
    assert_eq!(
        shortest_path(&input, AtomId::new(2), AtomId::new(2)),
        Err(PathError::EqualShortestPathEndpoints {
            atom: AtomId::new(2)
        })
    );
    assert!(matches!(
        shortest_path(&input, AtomId::new(9), AtomId::new(0)),
        Err(PathError::AtomOutOfRange { role: "begin", .. })
    ));
    assert!(matches!(
        shortest_path(&input, AtomId::new(0), AtomId::new(9)),
        Err(PathError::AtomOutOfRange { role: "end", .. })
    ));
}

#[test]
fn linear_paths_cover_bond_atom_range_root_ring_and_zero_branches() {
    let chain = carbon_topology(
        4,
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Single),
        ],
    );
    let bonds = all_paths_of_length(&chain, 2, &PathSearchParams::default()).unwrap();
    assert_eq!(
        bonds,
        vec![
            GraphPath::Bonds(bond_ids(&[0, 1])),
            GraphPath::Bonds(bond_ids(&[1, 2])),
        ]
    );
    let atom_params = PathSearchParams {
        representation: PathRepresentation::Atoms,
        ..PathSearchParams::default()
    };
    assert_eq!(
        all_paths_of_length(&chain, 3, &atom_params).unwrap(),
        vec![
            GraphPath::Atoms(atom_ids(&[0, 1, 2])),
            GraphPath::Atoms(atom_ids(&[1, 2, 3])),
        ]
    );
    let rooted = PathSearchParams {
        rooted_at_atom: Some(AtomId::new(0)),
        ..PathSearchParams::default()
    };
    assert_eq!(
        all_paths_of_length(&chain, 2, &rooted).unwrap(),
        vec![GraphPath::Bonds(bond_ids(&[0, 1]))]
    );
    let out_of_range = PathSearchParams {
        rooted_at_atom: Some(AtomId::new(99)),
        ..PathSearchParams::default()
    };
    assert!(
        all_paths_of_length(&chain, 2, &out_of_range)
            .unwrap()
            .is_empty()
    );
    assert!(
        all_paths_of_length(&chain, 0, &PathSearchParams::default())
            .unwrap()
            .is_empty()
    );
    assert_eq!(
        all_paths_in_range(&chain, 3, 2, &PathSearchParams::default()),
        Err(PathError::InvalidLengthRange { lower: 3, upper: 2 })
    );

    let triangle = carbon_topology(
        3,
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (2, 0, BondOrder::Single),
        ],
    );
    assert_eq!(
        all_paths_of_length(&triangle, 3, &PathSearchParams::default()).unwrap(),
        vec![GraphPath::Bonds(bond_ids(&[0, 1, 2]))]
    );
}

#[test]
fn path_hydrogen_filter_and_shortest_mode_preserve_source_interaction() {
    let input = topology(
        &[Element::C, Element::H, Element::O],
        &[(0, 1, BondOrder::Single), (1, 2, BondOrder::Single)],
    );
    assert!(
        all_paths_of_length(&input, 1, &PathSearchParams::default())
            .unwrap()
            .is_empty()
    );
    let with_hydrogen = PathSearchParams {
        use_hydrogens: true,
        ..PathSearchParams::default()
    };
    assert_eq!(
        all_paths_of_length(&input, 1, &with_hydrogen).unwrap(),
        vec![
            GraphPath::Bonds(bond_ids(&[0])),
            GraphPath::Bonds(bond_ids(&[1])),
        ]
    );
    let shortest = PathSearchParams {
        only_shortest_paths: true,
        ..PathSearchParams::default()
    };
    assert_eq!(
        all_paths_of_length(&input, 1, &shortest).unwrap(),
        vec![
            GraphPath::Bonds(bond_ids(&[0])),
            GraphPath::Bonds(bond_ids(&[1])),
        ]
    );
}

#[test]
fn branched_subgraphs_preserve_lifo_rows_range_keys_and_root_scope() {
    let branch = carbon_topology(
        4,
        &[
            (1, 0, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (1, 3, BondOrder::Single),
        ],
    );
    assert_eq!(
        all_subgraphs_of_length(&branch, 2, &SubgraphSearchParams::default()).unwrap(),
        vec![bond_ids(&[0, 2]), bond_ids(&[0, 1]), bond_ids(&[1, 2])]
    );
    let range = all_subgraphs_in_range(&branch, 0, 2, &SubgraphSearchParams::default()).unwrap();
    assert!(range[&0].is_empty());
    assert_eq!(
        range[&1],
        vec![bond_ids(&[0]), bond_ids(&[1]), bond_ids(&[2])]
    );
    assert_eq!(range[&2].len(), 3);
    let rooted = SubgraphSearchParams {
        rooted_at_atom: Some(AtomId::new(0)),
        ..SubgraphSearchParams::default()
    };
    assert_eq!(
        all_subgraphs_of_length(&branch, 2, &rooted).unwrap(),
        vec![bond_ids(&[0, 2]), bond_ids(&[0, 1])]
    );
    let missing_root = SubgraphSearchParams {
        rooted_at_atom: Some(AtomId::new(40)),
        ..SubgraphSearchParams::default()
    };
    assert!(
        all_subgraphs_of_length(&branch, 2, &missing_root)
            .unwrap()
            .is_empty()
    );
    assert_eq!(
        all_subgraphs_in_range(&branch, 5, 4, &SubgraphSearchParams::default()),
        Err(PathError::InvalidLengthRange { lower: 5, upper: 4 })
    );

    let hydrogen = topology(&[Element::C, Element::H], &[(0, 1, BondOrder::Single)]);
    assert!(
        all_subgraphs_of_length(&hydrogen, 1, &SubgraphSearchParams::default())
            .unwrap()
            .is_empty()
    );
    assert_eq!(
        all_subgraphs_of_length(
            &hydrogen,
            1,
            &SubgraphSearchParams {
                use_hydrogens: true,
                rooted_at_atom: None,
            }
        )
        .unwrap(),
        vec![bond_ids(&[0])]
    );
}

#[test]
fn unique_subgraphs_cover_bond_order_extra_invariants_and_length_error() {
    let input = carbon_topology(4, &[(0, 1, BondOrder::Single), (2, 3, BondOrder::Double)]);
    assert_eq!(
        unique_subgraphs_of_length(&input, 1, &UniqueSubgraphParams::default())
            .unwrap()
            .len(),
        2
    );
    let ignore_orders = UniqueSubgraphParams {
        use_bond_orders: false,
        ..UniqueSubgraphParams::default()
    };
    assert_eq!(
        unique_subgraphs_of_length(&input, 1, &ignore_orders)
            .unwrap()
            .len(),
        1
    );

    let equal_bonds = carbon_topology(4, &[(0, 1, BondOrder::Single), (2, 3, BondOrder::Single)]);
    let extra = UniqueSubgraphParams {
        extra_atom_invariants: Some(vec![1, 1, 2, 2]),
        ..UniqueSubgraphParams::default()
    };
    assert_eq!(
        unique_subgraphs_of_length(&equal_bonds, 1, &extra)
            .unwrap()
            .len(),
        2
    );
    let invalid = UniqueSubgraphParams {
        extra_atom_invariants: Some(vec![1]),
        ..UniqueSubgraphParams::default()
    };
    assert_eq!(
        unique_subgraphs_of_length(&equal_bonds, 1, &invalid),
        Err(PathError::ExtraInvariantLength {
            actual: 1,
            expected: 4
        })
    );
}

#[test]
fn unique_discriminator_observes_charge_isotope_and_aromaticity() {
    let base = carbon_topology(4, &[(0, 1, BondOrder::Single), (2, 3, BondOrder::Single)]);
    assert_eq!(
        unique_subgraphs_of_length(&base, 1, &UniqueSubgraphParams::default())
            .unwrap()
            .len(),
        1
    );

    let mut charged = base.clone();
    charged.atoms[2].set_formal_charge(1);
    assert_eq!(
        unique_subgraphs_of_length(&charged, 1, &UniqueSubgraphParams::default())
            .unwrap()
            .len(),
        2
    );
    let mut isotopic = base.clone();
    // The pinned source truncates (exact isotope mass - average atomic mass)
    // to an integer. Carbon-13 differs by less than one and intentionally
    // collides; carbon-14 exercises the non-zero source discriminator branch.
    isotopic.atoms[2].set_isotope(Some(14));
    assert_eq!(
        unique_subgraphs_of_length(&isotopic, 1, &UniqueSubgraphParams::default())
            .unwrap()
            .len(),
        2
    );
    let mut aromatic = base;
    aromatic.atoms[2].set_aromatic(true);
    assert_eq!(
        unique_subgraphs_of_length(&aromatic, 1, &UniqueSubgraphParams::default())
            .unwrap()
            .len(),
        2
    );
}

#[test]
fn atom_environment_covers_radius_hydrogen_and_enforcement_matrix() {
    let chain = carbon_topology(3, &[(0, 1, BondOrder::Single), (1, 2, BondOrder::Single)]);
    assert_eq!(
        atom_environment(&chain, 0, AtomId::new(0), &AtomEnvironmentParams::default()).unwrap(),
        cosmolkit_core::AtomEnvironment {
            bonds: vec![],
            atom_distances: vec![Some(0), None, None],
        }
    );
    let radius_two =
        atom_environment(&chain, 2, AtomId::new(0), &AtomEnvironmentParams::default()).unwrap();
    assert_eq!(radius_two.bonds, bond_ids(&[0, 1]));
    assert_eq!(radius_two.atom_distances, vec![Some(0), Some(1), Some(2)]);
    let cleared =
        atom_environment(&chain, 4, AtomId::new(0), &AtomEnvironmentParams::default()).unwrap();
    assert!(cleared.bonds.is_empty());
    assert_eq!(cleared.atom_distances, vec![None, None, None]);
    let partial = atom_environment(
        &chain,
        4,
        AtomId::new(0),
        &AtomEnvironmentParams {
            enforce_radius: false,
            ..AtomEnvironmentParams::default()
        },
    )
    .unwrap();
    assert_eq!(partial.bonds, bond_ids(&[0, 1]));
    assert_eq!(partial.atom_distances, vec![Some(0), Some(1), Some(2)]);

    let hydrogen = topology(&[Element::C, Element::H], &[(0, 1, BondOrder::Single)]);
    assert!(
        atom_environment(
            &hydrogen,
            1,
            AtomId::new(0),
            &AtomEnvironmentParams::default()
        )
        .unwrap()
        .bonds
        .is_empty()
    );
    assert_eq!(
        atom_environment(
            &hydrogen,
            1,
            AtomId::new(0),
            &AtomEnvironmentParams {
                use_hydrogens: true,
                enforce_radius: true,
            }
        )
        .unwrap()
        .bonds,
        bond_ids(&[0])
    );
    assert!(matches!(
        atom_environment(&chain, 1, AtomId::new(8), &AtomEnvironmentParams::default()),
        Err(PathError::AtomOutOfRange { role: "root", .. })
    ));
}

#[test]
fn atom_path_projection_uses_all_pairs_in_nested_input_order() {
    let triangle = carbon_topology(
        3,
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (2, 0, BondOrder::Single),
        ],
    );
    assert_eq!(
        bond_ids_from_atom_path(&triangle, &atom_ids(&[2, 0, 1])).unwrap(),
        bond_ids(&[2, 1, 0])
    );
    assert!(bond_ids_from_atom_path(&triangle, &[]).unwrap().is_empty());
    assert!(
        bond_ids_from_atom_path(&triangle, &[AtomId::new(1)])
            .unwrap()
            .is_empty()
    );
    assert_eq!(
        bond_ids_from_atom_path(&triangle, &[AtomId::new(0), AtomId::new(7)]),
        Err(PathError::AtomPathOutOfRange {
            position: 1,
            atom: AtomId::new(7),
            atom_count: 3,
        })
    );
}

#[test]
fn detached_subset_coalesces_ids_remaps_state_and_clears_computed_props() {
    let mut input = carbon_topology(
        4,
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Double),
            (2, 3, BondOrder::Single),
        ],
    );
    input.atoms[1].set_prop("kept", "atom").unwrap();
    input.atoms[1]
        .set_computed_prop("computed", "drop")
        .unwrap();
    input.bonds[1].set_prop("kept", "bond").unwrap();
    input.bonds[1]
        .set_computed_prop("computed", "drop")
        .unwrap();
    input.bonds[1].set_stereo_atoms(Some([AtomId::new(0), AtomId::new(3)]));
    input.bonds[1].set_stereo(BondStereo::Cis).unwrap();
    input.substance_groups = vec![
        SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
            .with_atoms(atom_ids(&[1, 2]))
            .with_bonds(bond_ids(&[1])),
        SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Data)
            .with_atoms(atom_ids(&[0, 1])),
    ];
    input.stereo_groups = vec![
        StereoGroup::new(StereoGroupKind::Or, atom_ids(&[0, 1]), bond_ids(&[0, 1])).with_id(7),
    ];
    input.validate().unwrap();
    let snapshot = input.clone();

    let result = subtopology_from_path(
        &input,
        &[BondId::new(1), BondId::new(1), BondId::new(99)],
        &SubtopologyParams::default(),
    )
    .unwrap();
    assert_eq!(input, snapshot);
    assert_eq!(
        result.mapping.atoms.old_to_new,
        vec![None, Some(AtomId::new(0)), Some(AtomId::new(1)), None]
    );
    assert_eq!(
        result.mapping.atoms.new_to_old,
        vec![Some(AtomId::new(1)), Some(AtomId::new(2))]
    );
    assert_eq!(
        result.mapping.bonds.old_to_new,
        vec![None, Some(BondId::new(0)), None]
    );
    assert_eq!(result.mapping.bonds.new_to_old, vec![Some(BondId::new(1))]);
    let DetachedPathSubgraph::Concrete(subset) = result.subgraph else {
        panic!("default subset must be concrete")
    };
    assert_eq!(subset.atoms.len(), 2);
    assert_eq!(subset.bonds.len(), 1);
    assert_eq!(subset.atoms[0].prop("kept"), Some("atom"));
    assert_eq!(subset.atoms[0].prop("computed"), None);
    assert_eq!(subset.bonds[0].prop("kept"), Some("bond"));
    assert_eq!(subset.bonds[0].prop("computed"), None);
    assert_eq!(subset.bonds[0].stereo(), BondStereo::None);
    assert_eq!(subset.bonds[0].stereo_atoms(), None);
    assert_eq!(subset.substance_groups.len(), 1);
    assert_eq!(subset.substance_groups[0].atoms(), atom_ids(&[0, 1]));
    assert_eq!(subset.substance_groups[0].bonds(), bond_ids(&[0]));
    assert_eq!(subset.stereo_groups.len(), 1);
    assert_eq!(subset.stereo_groups[0].id(), Some(7));
    assert_eq!(subset.stereo_groups[0].atoms(), atom_ids(&[0]));
    assert_eq!(subset.stereo_groups[0].bonds(), bond_ids(&[0]));
    subset.validate().unwrap();
}

#[test]
fn query_subset_uses_exact_predicates_and_source_table_order() {
    let input = topology(
        &[Element::N, Element::C, Element::O, Element::S],
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Double),
            (2, 3, BondOrder::Aromatic),
        ],
    );
    let first = subtopology_from_path(
        &input,
        &[BondId::new(2), BondId::new(0)],
        &SubtopologyParams {
            copy_as_query: true,
        },
    )
    .unwrap();
    let second = subtopology_from_path(
        &input,
        &[BondId::new(0), BondId::new(2)],
        &SubtopologyParams {
            copy_as_query: true,
        },
    )
    .unwrap();
    assert_eq!(first, second);
    let DetachedPathSubgraph::Query {
        graph,
        substance_groups,
    } = first.subgraph
    else {
        panic!("copy_as_query must return a query graph")
    };
    assert!(substance_groups.is_empty());
    assert_eq!(graph.num_atoms(), 4);
    assert_eq!(graph.num_bonds(), 2);
    assert_eq!(
        graph
            .atoms()
            .iter()
            .map(|atom| atom
                .element()
                .expect("query carrier has an Element identity"))
            .collect::<Vec<_>>(),
        vec![Element::N, Element::C, Element::O, Element::S]
    );
    assert_eq!(
        graph.atoms()[0].predicate(),
        &QueryNode::predicate(AtomQueryPredicate::AtomicNumber(7))
    );
    assert_eq!(
        graph.bonds()[0].predicate(),
        &QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single))
    );
    assert_eq!(
        graph.bonds()[1].predicate(),
        &QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Aromatic))
    );
    graph.validate().unwrap();
}

#[test]
fn every_entry_rejects_invalid_detached_topology_before_graph_work() {
    let valid = carbon_topology(2, &[(0, 1, BondOrder::Single)]);
    let invalid = TopologyBlock {
        adjacency: AdjacencyList::from_topology(2, &[]),
        ..valid
    };
    assert_eq!(
        connected_components(&invalid),
        Err(PathError::InvalidTopology(
            TopologyValidationError::AdjacencyMismatch
        ))
    );
    assert!(matches!(
        shortest_path(&invalid, AtomId::new(0), AtomId::new(1)),
        Err(PathError::InvalidTopology(_))
    ));
    assert!(matches!(
        all_paths_of_length(&invalid, 1, &PathSearchParams::default()),
        Err(PathError::InvalidTopology(_))
    ));
    assert!(matches!(
        all_subgraphs_of_length(&invalid, 1, &SubgraphSearchParams::default()),
        Err(PathError::InvalidTopology(_))
    ));
    assert!(matches!(
        unique_subgraphs_of_length(&invalid, 1, &UniqueSubgraphParams::default()),
        Err(PathError::InvalidTopology(_))
    ));
    assert!(matches!(
        atom_environment(
            &invalid,
            1,
            AtomId::new(0),
            &AtomEnvironmentParams::default()
        ),
        Err(PathError::InvalidTopology(_))
    ));
    assert!(matches!(
        bond_ids_from_atom_path(&invalid, &[AtomId::new(0)]),
        Err(PathError::InvalidTopology(_))
    ));
    assert!(matches!(
        subtopology_from_path(&invalid, &[BondId::new(0)], &SubtopologyParams::default()),
        Err(PathError::InvalidTopology(_))
    ));
}
