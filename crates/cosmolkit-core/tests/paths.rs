use cosmolkit_core::{
    AtomEnvironmentParams, DetachedPathSubgraph, GraphPath, PathError, PathRepresentation,
    PathSearchParams, SubgraphSearchParams, SubtopologyParams, UniqueSubgraphParams,
    all_paths_in_range, all_paths_of_length, all_subgraphs_in_range, all_subgraphs_of_length,
    atom_environment, bond_ids_from_atom_path, connected_components, shortest_path,
    subtopology_from_path, unique_subgraphs_of_length,
};
use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomQueryPredicate, AtomSpec, Bond, BondId, BondQueryPredicate,
    BondSpec, BondStereo, PropertyValue, QueryNode, StereoGroup, StereoGroupKind, SubstanceGroup,
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
                ignore_atoms: None,
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
        StereoGroup::new(StereoGroupKind::Or, atom_ids(&[0, 1]), bond_ids(&[0, 1]))
            .expect("valid distinct stereo members")
            .with_id(7),
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
    assert_eq!(
        subset.atoms[0].prop("kept"),
        Some(&PropertyValue::String("atom".to_owned().into()))
    );
    assert_eq!(subset.atoms[0].prop("computed"), None);
    assert_eq!(
        subset.bonds[0].prop("kept"),
        Some(&PropertyValue::String("bond".to_owned().into()))
    );
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

#[test]
fn search01_mask_applies_to_single_unique_and_range_before_root() {
    let t = carbon_topology(
        4,
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (1, 3, BondOrder::Single),
        ],
    );
    let mask = [true, false, false, false];
    let p = SubgraphSearchParams {
        ignore_atoms: Some(&mask),
        rooted_at_atom: Some(AtomId::new(1)),
        ..Default::default()
    };
    assert_eq!(
        all_subgraphs_of_length(&t, 1, &p).unwrap(),
        vec![bond_ids(&[1]), bond_ids(&[2])]
    );
    assert_eq!(
        all_subgraphs_of_length(&t, 2, &p).unwrap(),
        vec![bond_ids(&[1, 2])]
    );
    assert_eq!(
        all_subgraphs_in_range(&t, 1, 2, &p).unwrap()[&2],
        vec![bond_ids(&[1, 2])]
    );
    let u = UniqueSubgraphParams {
        ignore_atoms: Some(&mask),
        rooted_at_atom: p.rooted_at_atom,
        ..Default::default()
    };
    assert_eq!(
        unique_subgraphs_of_length(&t, 2, &u).unwrap(),
        vec![bond_ids(&[1, 2])]
    );
    let root_ignored = SubgraphSearchParams {
        rooted_at_atom: Some(AtomId::new(0)),
        ..p
    };
    assert!(
        all_subgraphs_of_length(&t, 1, &root_ignored)
            .unwrap()
            .is_empty()
    );
}
#[test]
fn search01_mask_shape_is_checked_even_when_single_target_is_zero() {
    let t = carbon_topology(2, &[(0, 1, BondOrder::Single)]);
    let p = SubgraphSearchParams {
        ignore_atoms: Some(&[]),
        ..Default::default()
    };
    assert!(matches!(
        all_subgraphs_of_length(&t, 0, &p),
        Err(PathError::IgnoredAtomMaskLength {
            actual: 0,
            expected: 2
        })
    ));
    assert!(matches!(
        all_subgraphs_in_range(&t, 0, 0, &p),
        Err(PathError::IgnoredAtomMaskLength { .. })
    ));
    let path = PathSearchParams {
        ignore_atoms: Some(&[]),
        ..Default::default()
    };
    assert!(matches!(
        all_paths_of_length(&t, 0, &path),
        Err(PathError::IgnoredAtomMaskLength { .. })
    ));
}
#[test]
fn search01_atom_seeds_extensions_and_fullgraph_shortest_distances() {
    let t = carbon_topology(
        4,
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Single),
            (0, 3, BondOrder::Single),
        ],
    );
    let mask = [false, true, false, false];
    let p = PathSearchParams {
        ignore_atoms: Some(&mask),
        representation: PathRepresentation::Atoms,
        ..Default::default()
    };
    assert_eq!(
        all_paths_of_length(&t, 1, &p).unwrap(),
        vec![
            GraphPath::Atoms(atom_ids(&[0])),
            GraphPath::Atoms(atom_ids(&[2])),
            GraphPath::Atoms(atom_ids(&[3]))
        ]
    );
    let rooted = PathSearchParams {
        rooted_at_atom: Some(AtomId::new(0)),
        only_shortest_paths: true,
        ..p
    };
    assert_eq!(
        all_paths_of_length(&t, 3, &rooted).unwrap(),
        vec![GraphPath::Atoms(atom_ids(&[0, 3, 2]))]
    );
    let shortcut = carbon_topology(
        5,
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (0, 4, BondOrder::Single),
            (4, 3, BondOrder::Single),
            (3, 2, BondOrder::Single),
        ],
    );
    let shortcut_mask = [false, true, false, false, false];
    let detour = PathSearchParams {
        ignore_atoms: Some(&shortcut_mask),
        representation: PathRepresentation::Atoms,
        rooted_at_atom: Some(AtomId::new(0)),
        ..Default::default()
    };
    assert_eq!(
        all_paths_of_length(&shortcut, 4, &detour).unwrap(),
        vec![GraphPath::Atoms(atom_ids(&[0, 4, 3, 2]))]
    );
    let shortest_detour = PathSearchParams {
        only_shortest_paths: true,
        ..detour
    };
    assert!(
        all_paths_of_length(&shortcut, 4, &shortest_detour)
            .unwrap()
            .is_empty()
    );
    let ignored_root = PathSearchParams {
        rooted_at_atom: Some(AtomId::new(1)),
        ..p
    };
    assert!(
        all_paths_of_length(&t, 1, &ignored_root)
            .unwrap()
            .is_empty()
    );
}
#[test]
fn search01_some_empty_mask_differs_from_none_only_for_nonempty_graph() {
    let empty = carbon_topology(0, &[]);
    let p = SubgraphSearchParams {
        ignore_atoms: Some(&[]),
        ..Default::default()
    };
    assert!(all_subgraphs_of_length(&empty, 1, &p).unwrap().is_empty());
    let t = carbon_topology(3, &[(0, 1, BondOrder::Single), (1, 2, BondOrder::Single)]);
    let no_ignored = [false; 3];
    for size in 0..=3 {
        assert_eq!(
            all_subgraphs_of_length(&t, size, &Default::default()).unwrap(),
            all_subgraphs_of_length(
                &t,
                size,
                &SubgraphSearchParams {
                    ignore_atoms: Some(&no_ignored),
                    ..Default::default()
                }
            )
            .unwrap()
        );
    }
}

#[test]
fn detached_subset_uses_source_controller_replacement_and_stereo_parity() {
    // RDKit .6 Subset::handleBondStereo: selected original controllers retain
    // parity; each replaced controller toggles it; E/Z become CIS/TRANS.
    // Direct model construction is independent of any SMILES corpus or oracle.
    for (selected, controllers, swaps) in [
        (vec![0, 1, 2], [0, 3], false),
        (vec![1, 2, 3], [3, 2], true),
        (vec![0, 1, 4], [0, 3], true),
        (vec![1, 3, 4], [2, 3], false),
    ] {
        for (source_stereo, unswapped, swapped) in [
            (BondStereo::E, BondStereo::Trans, BondStereo::Cis),
            (BondStereo::Z, BondStereo::Cis, BondStereo::Trans),
            (BondStereo::Cis, BondStereo::Cis, BondStereo::Trans),
            (BondStereo::Trans, BondStereo::Trans, BondStereo::Cis),
        ] {
            for copy_as_query in [false, true] {
                let mut input = carbon_topology(
                    6,
                    &[
                        (0, 1, BondOrder::Single),
                        (1, 2, BondOrder::Double),
                        (2, 3, BondOrder::Single),
                        (1, 4, BondOrder::Single),
                        (2, 5, BondOrder::Single),
                    ],
                );
                input.bonds[1].set_stereo_atoms(Some([AtomId::new(0), AtomId::new(3)]));
                input.bonds[1].set_stereo(source_stereo).unwrap();
                input.validate().unwrap();
                let snapshot = input.clone();
                let result = subtopology_from_path(
                    &input,
                    &bond_ids(&selected),
                    &SubtopologyParams { copy_as_query },
                )
                .unwrap();
                assert_eq!(input, snapshot);
                let row = result.mapping.bonds.old_to_new[1].unwrap().index();
                let expected_stereo = if swaps { swapped } else { unswapped };
                let expected_controllers =
                    Some([AtomId::new(controllers[0]), AtomId::new(controllers[1])]);
                match &result.subgraph {
                    DetachedPathSubgraph::Concrete(subset) => {
                        assert_eq!(subset.bonds[row].stereo(), expected_stereo);
                        assert_eq!(subset.bonds[row].stereo_atoms(), expected_controllers);
                        subset.validate().unwrap();
                    }
                    DetachedPathSubgraph::Query { graph, .. } => {
                        assert_eq!(graph.bonds()[row].bond().stereo(), expected_stereo);
                        assert_eq!(
                            graph.bonds()[row].bond().stereo_atoms(),
                            expected_controllers
                        );
                        assert_eq!(
                            graph.bonds()[row].predicate(),
                            &QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Double))
                        );
                        graph.validate().unwrap();
                    }
                }
            }
        }
    }
}

#[test]
fn detached_subset_clears_all_defined_stereo_when_no_controller_can_be_copied() {
    // Cover both native early returns: source degree < 3, and a degree >= 3
    // alternate that is absent from atomMapping. Source clears E/Z as well.
    for branched in [false, true] {
        for selected in [vec![1], vec![1, 2]] {
            for source_stereo in [
                BondStereo::E,
                BondStereo::Z,
                BondStereo::Cis,
                BondStereo::Trans,
            ] {
                for copy_as_query in [false, true] {
                    let mut edges = vec![
                        (0, 1, BondOrder::Single),
                        (1, 2, BondOrder::Double),
                        (2, 3, BondOrder::Single),
                    ];
                    if branched {
                        edges.extend([(1, 4, BondOrder::Single), (2, 5, BondOrder::Single)]);
                    }
                    let mut input = carbon_topology(if branched { 6 } else { 4 }, &edges);
                    input.bonds[1].set_stereo_atoms(Some([AtomId::new(0), AtomId::new(3)]));
                    input.bonds[1].set_stereo(source_stereo).unwrap();
                    input.validate().unwrap();
                    let snapshot = input.clone();
                    let result = subtopology_from_path(
                        &input,
                        &bond_ids(&selected),
                        &SubtopologyParams { copy_as_query },
                    )
                    .unwrap();
                    assert_eq!(input, snapshot);
                    let row = result.mapping.bonds.old_to_new[1].unwrap().index();
                    match &result.subgraph {
                        DetachedPathSubgraph::Concrete(subset) => {
                            assert_eq!(subset.bonds[row].stereo(), BondStereo::None);
                            assert_eq!(subset.bonds[row].stereo_atoms(), None);
                            subset.validate().unwrap();
                        }
                        DetachedPathSubgraph::Query { graph, .. } => {
                            assert_eq!(graph.bonds()[row].bond().stereo(), BondStereo::None);
                            assert_eq!(graph.bonds()[row].bond().stereo_atoms(), None);
                            graph.validate().unwrap();
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn detached_subset_uses_the_first_source_alternate_before_testing_its_mapping() {
    // getOtherAtomIdx returns the first source neighbor, not the first copied
    // neighbor. Atom 4 precedes atom 5 but is omitted: source must clear stereo.
    for copy_as_query in [false, true] {
        let mut input = carbon_topology(
            7,
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Double),
                (2, 3, BondOrder::Single),
                (1, 4, BondOrder::Single),
                (1, 5, BondOrder::Single),
                (2, 6, BondOrder::Single),
            ],
        );
        input.bonds[1].set_stereo_atoms(Some([AtomId::new(0), AtomId::new(3)]));
        input.bonds[1].set_stereo(BondStereo::E).unwrap();
        input.validate().unwrap();
        let snapshot = input.clone();
        let result = subtopology_from_path(
            &input,
            &bond_ids(&[1, 2, 4]),
            &SubtopologyParams { copy_as_query },
        )
        .unwrap();
        assert_eq!(input, snapshot);
        let row = result.mapping.bonds.old_to_new[1].unwrap().index();
        match &result.subgraph {
            DetachedPathSubgraph::Concrete(subset) => {
                assert_eq!(subset.bonds[row].stereo(), BondStereo::None);
                assert_eq!(subset.bonds[row].stereo_atoms(), None);
                subset.validate().unwrap();
            }
            DetachedPathSubgraph::Query { graph, .. } => {
                assert_eq!(graph.bonds()[row].bond().stereo(), BondStereo::None);
                assert_eq!(graph.bonds()[row].bond().stereo_atoms(), None);
                graph.validate().unwrap();
            }
        }
    }
}
