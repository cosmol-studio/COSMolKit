use cosmolkit_core::{
    KekulizeAttempt, KekulizeError, KekulizeParams, kekulize, kekulize_if_possible,
    kekulize_if_possible_with_query_state, kekulize_with_query_state,
};
use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomQueryPredicate, AtomSpec, Bond, BondId, BondQueryPredicate,
    BondSpec, PropertyValue, QueryAtom, QueryBond, QueryNode, QueryStateRef, TopologyBlock,
    TopologyValidationError,
};
use cosmolkit_types::{BondDirection, BondOrder, Element};

fn atom(id: usize, spec: AtomSpec) -> Atom {
    Atom::from_spec(AtomId::new(id), spec)
}

fn bond(id: usize, spec: BondSpec) -> Bond {
    Bond::from_spec(BondId::new(id), spec)
}

fn aromatic_bond(id: usize, begin: usize, end: usize) -> Bond {
    bond(
        id,
        BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Aromatic)
            .with_aromatic(true),
    )
}

fn topology(atoms: Vec<Atom>, bonds: Vec<Bond>) -> TopologyBlock {
    TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
}

fn aromatic_cycle(specs: Vec<AtomSpec>) -> TopologyBlock {
    let atom_count = specs.len();
    topology(
        specs
            .into_iter()
            .enumerate()
            .map(|(id, spec)| atom(id, spec.with_aromatic(true)))
            .collect(),
        (0..atom_count)
            .map(|id| aromatic_bond(id, id, (id + 1) % atom_count))
            .collect(),
    )
}

fn aromatic_carbon_cycle(size: usize) -> TopologyBlock {
    aromatic_cycle((0..size).map(|_| AtomSpec::new(Element::C)).collect())
}

fn bond_orders(topology: &TopologyBlock) -> Vec<BondOrder> {
    topology.bonds.iter().map(Bond::order).collect()
}

fn assert_kekule_ring(topology: &TopologyBlock, ring_bonds: usize, double_bonds: usize) {
    assert_eq!(
        topology.bonds[..ring_bonds]
            .iter()
            .filter(|bond| bond.order() == BondOrder::Double)
            .count(),
        double_bonds
    );
    assert!(
        topology.bonds[..ring_bonds]
            .iter()
            .all(|bond| matches!(bond.order(), BondOrder::Single | BondOrder::Double))
    );
    topology.validate().unwrap();
}

#[test]
fn kekulize_source_state_independent_flags_parameter_product() {
    // Frozen SOURCE-K A48: Kekulize.cpp selection/markDbondCands/finalization
    // and Atom.cpp isAromaticAtom. False bond flags do not enter makeSingle.
    let mut calls = 0;
    let mut outcomes = Vec::new();
    for graph in 0..3 {
        for mark_atoms_bonds in [false, true] {
            for canonical in [false, true] {
                for max_backtracks in [0, 100] {
                    for if_possible in [false, true] {
                        let size = if graph == 0 { 2 } else { 6 };
                        let input = topology(
                            (0..size)
                                .map(|id| {
                                    let element = if graph == 2 && id == 0 {
                                        Element::N
                                    } else {
                                        Element::C
                                    };
                                    atom(id, AtomSpec::new(element).with_aromatic(graph != 0))
                                })
                                .collect(),
                            (0..if graph == 0 { 1 } else { 6 })
                                .map(|id| {
                                    bond(
                                        id,
                                        BondSpec::new(
                                            AtomId::new(id),
                                            AtomId::new((id + 1) % size),
                                            BondOrder::Aromatic,
                                        )
                                        .with_aromatic(false),
                                    )
                                })
                                .collect(),
                        );
                        for (id, row) in input.atoms.iter().enumerate() {
                            assert_eq!(row.id(), AtomId::new(id));
                            assert_eq!(
                                row.element(),
                                if graph == 2 && id == 0 {
                                    Element::N
                                } else {
                                    Element::C
                                }
                            );
                            assert_eq!(row.is_aromatic(), graph != 0);
                        }
                        for (id, row) in input.bonds.iter().enumerate() {
                            assert_eq!(row.id(), BondId::new(id));
                            assert_eq!(row.begin(), AtomId::new(id));
                            assert_eq!(row.end(), AtomId::new((id + 1) % size));
                            assert_eq!(row.order(), BondOrder::Aromatic);
                            assert!(!row.is_aromatic());
                        }
                        let baseline = input.clone();
                        let mut expected = baseline.clone();
                        if graph != 0 && mark_atoms_bonds {
                            for row in &mut expected.atoms {
                                row.set_aromatic(false);
                            }
                        }
                        let params = KekulizeParams {
                            mark_atoms_bonds,
                            canonical,
                            max_backtracks,
                        };
                        let result = if if_possible {
                            kekulize_if_possible(&input, &params).map(|attempt| match attempt {
                                KekulizeAttempt::Applied(assignment) => Some(assignment.topology),
                                KekulizeAttempt::NotKekulizable { .. } => None,
                            })
                        } else {
                            kekulize(&input, &params).map(|assignment| Some(assignment.topology))
                        };
                        calls += 1;
                        assert_eq!(input, baseline, "retained input changed at call {calls}");
                        outcomes.push((graph, params, if_possible, result, expected));
                    }
                }
            }
        }
    }
    assert_eq!(calls, 48);
    eprintln!("SOURCE-K A48 actual whole-entry calls: {calls}");
    let failures = outcomes.into_iter().filter_map(|(graph, params, if_possible, result, expected)| {
        if result.as_ref().is_ok_and(|output| output.as_ref() == Some(&expected)) {
            None
        } else {
            Some(format!("graph={graph}, params={params:?}, if_possible={if_possible}: {result:?}; expected={expected:?}"))
        }
    }).collect::<Vec<_>>();
    assert!(
        failures.is_empty(),
        "{} of 48 source rows failed:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

#[test]
fn public_defaults_empty_and_nonaromatic_inputs_are_exact_noops() {
    assert_eq!(
        KekulizeParams::default(),
        KekulizeParams {
            mark_atoms_bonds: true,
            canonical: true,
            max_backtracks: 100,
        }
    );

    let empty = TopologyBlock::default();
    assert_eq!(
        kekulize(&empty, &KekulizeParams::default())
            .unwrap()
            .topology,
        empty
    );

    let nonaromatic = topology(
        vec![
            atom(
                0,
                AtomSpec::new(Element::C)
                    .with_prop("atom-note", "left")
                    .unwrap(),
            ),
            atom(1, AtomSpec::new(Element::O)),
        ],
        vec![bond(
            0,
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
                .with_direction(BondDirection::BeginDash)
                .with_prop("bond-note", "unchanged")
                .unwrap(),
        )],
    );
    let snapshot = nonaromatic.clone();
    let assignment = kekulize(&nonaromatic, &KekulizeParams::default()).unwrap();
    assert_eq!(assignment.topology, snapshot);
    assert_eq!(nonaromatic, snapshot);
    assert_eq!(
        kekulize_if_possible(&nonaromatic, &KekulizeParams::default()).unwrap(),
        KekulizeAttempt::Applied(assignment)
    );
}

#[test]
fn mark_true_kekulizes_benzene_and_preserves_unrelated_rows_and_input() {
    let mut input = aromatic_carbon_cycle(6);
    input.atoms[0].set_prop("atom-note", "preserved").unwrap();
    input.bonds[0].set_prop("ring-note", "preserved").unwrap();
    input.atoms.push(atom(6, AtomSpec::new(Element::C)));
    input.atoms.push(atom(7, AtomSpec::new(Element::O)));
    input.bonds.push(bond(
        6,
        BondSpec::new(AtomId::new(6), AtomId::new(7), BondOrder::Single)
            .with_direction(BondDirection::BeginWedge)
            .with_prop("outside", "preserved")
            .unwrap(),
    ));
    input.adjacency = AdjacencyList::try_from_topology(input.atoms.len(), &input.bonds).unwrap();
    input.validate().unwrap();
    let snapshot = input.clone();

    let output = kekulize(&input, &KekulizeParams::default())
        .unwrap()
        .topology;
    assert_kekule_ring(&output, 6, 3);
    assert!(output.atoms[..6].iter().all(|atom| !atom.is_aromatic()));
    assert!(output.bonds[..6].iter().all(|bond| !bond.is_aromatic()));
    assert_eq!(
        output.atoms[0].prop("atom-note"),
        Some(&PropertyValue::String("preserved".to_owned().into()))
    );
    assert_eq!(
        output.bonds[0].prop("ring-note"),
        Some(&PropertyValue::String("preserved".to_owned().into()))
    );
    assert_eq!(output.bonds[6], snapshot.bonds[6]);
    assert_eq!(output.adjacency, snapshot.adjacency);
    assert_eq!(input, snapshot);
}

#[test]
fn mark_false_changes_orders_but_preserves_aromatic_flags_and_is_deterministic() {
    let input = aromatic_carbon_cycle(6);
    let snapshot = input.clone();
    let params = KekulizeParams {
        mark_atoms_bonds: false,
        canonical: false,
        max_backtracks: 100,
    };
    let first = kekulize(&input, &params).unwrap().topology;
    let second = kekulize(&input, &params).unwrap().topology;

    assert_eq!(first, second);
    assert_kekule_ring(&first, 6, 3);
    assert!(first.atoms.iter().all(Atom::is_aromatic));
    assert!(first.bonds.iter().all(Bond::is_aromatic));
    assert_eq!(input, snapshot);

    let canonical = kekulize(
        &input,
        &KekulizeParams {
            canonical: true,
            ..params
        },
    )
    .unwrap()
    .topology;
    assert_kekule_ring(&canonical, 6, 3);
    assert_eq!(
        canonical,
        kekulize(
            &input,
            &KekulizeParams {
                canonical: true,
                ..params
            }
        )
        .unwrap()
        .topology
    );
}

#[test]
fn source_fused_and_disconnected_systems_are_completed_in_stable_order() {
    let graph = topology(
        (0..16)
            .map(|id| atom(id, AtomSpec::new(Element::C).with_aromatic(true)))
            .collect(),
        [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (5, 6),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 4),
            (10, 11),
            (11, 12),
            (12, 13),
            (13, 14),
            (14, 15),
            (15, 10),
        ]
        .into_iter()
        .enumerate()
        .map(|(id, (begin, end))| aromatic_bond(id, begin, end))
        .collect(),
    );
    let snapshot = graph.clone();
    let output = kekulize(&graph, &KekulizeParams::default())
        .unwrap()
        .topology;
    assert_eq!(
        output
            .bonds
            .iter()
            .filter(|bond| bond.order() == BondOrder::Double)
            .count(),
        8
    );
    assert!(output.bonds.iter().all(|bond| !bond.is_aromatic()));
    assert!(output.atoms.iter().all(|atom| !atom.is_aromatic()));
    assert_eq!(output.adjacency, snapshot.adjacency);
    assert_eq!(graph, snapshot);
}

#[test]
fn query_bonds_and_all_dummy_rings_follow_source_exclusion_rules() {
    let query_ring = topology(
        (0..3)
            .map(|id| atom(id, AtomSpec::new(Element::C).with_aromatic(true)))
            .collect(),
        vec![
            aromatic_bond(0, 0, 1),
            aromatic_bond(1, 1, 2),
            aromatic_bond(2, 2, 0),
        ],
    );
    let query_atoms = query_ring
        .atoms
        .iter()
        .map(|carrier| {
            QueryAtom::from_carrier_parts(
                carrier.clone(),
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
            )
        })
        .collect::<Vec<_>>();
    let mut query_bonds = query_ring
        .bonds
        .iter()
        .map(|carrier| {
            QueryBond::from_carrier_parts(
                carrier.clone(),
                QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Aromatic)),
            )
        })
        .collect::<Vec<_>>();
    query_bonds[0] = QueryBond::from_parts(
        query_ring.bonds[0].clone(),
        QueryNode::or(vec![
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Aromatic)),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
        ]),
    );
    let query_state =
        QueryStateRef::try_for_topology(&query_atoms, &query_bonds, &query_ring).unwrap();
    let keep_flags = KekulizeParams {
        mark_atoms_bonds: false,
        ..KekulizeParams::default()
    };
    assert_eq!(
        kekulize_with_query_state(&query_ring, &keep_flags, Some(query_state))
            .unwrap()
            .topology,
        query_ring
    );

    let all_dummy = aromatic_cycle((0..6).map(|_| AtomSpec::new(Element::DUMMY)).collect());
    let output = kekulize(&all_dummy, &KekulizeParams::default())
        .unwrap()
        .topology;
    assert_eq!(bond_orders(&output), vec![BondOrder::Aromatic; 6]);
    assert_eq!(
        all_dummy.bonds.iter().map(Bond::order).collect::<Vec<_>>(),
        vec![BondOrder::Aromatic; 6]
    );
}

#[test]
fn odd_aromatic_failure_has_ordered_atoms_and_if_possible_restores_every_row() {
    let mut input = aromatic_carbon_cycle(5);
    input.atoms[2].set_prop("atom-note", "restore").unwrap();
    input.bonds[3].set_prop("bond-note", "restore").unwrap();
    input.bonds[4].set_direction(BondDirection::BeginDash);
    let snapshot = input.clone();

    assert!(matches!(
        kekulize(&input, &KekulizeParams::default()),
        Err(KekulizeError::NotKekulizable { problem_atoms })
            if problem_atoms == (0..5).map(AtomId::new).collect::<Vec<_>>()
    ));
    assert_eq!(input, snapshot);

    match kekulize_if_possible(&input, &KekulizeParams::default()).unwrap() {
        KekulizeAttempt::NotKekulizable {
            topology,
            problem_atoms,
        } => {
            assert_eq!(topology, snapshot);
            assert_eq!(problem_atoms, (0..5).map(AtomId::new).collect::<Vec<_>>());
        }
        KekulizeAttempt::Applied(_) => panic!("odd carbon ring must not be reported as applied"),
    }
    assert_eq!(input, snapshot);
}

#[test]
fn non_ring_aromatic_failure_reports_atom_and_restores_for_if_possible() {
    let input = topology(
        vec![atom(
            0,
            AtomSpec::new(Element::C)
                .with_aromatic(true)
                .with_prop("source", "unchanged")
                .unwrap(),
        )],
        vec![],
    );
    let snapshot = input.clone();
    assert!(matches!(
        kekulize(&input, &KekulizeParams::default()),
        Err(KekulizeError::AromaticAtomOutsideRing { atom }) if atom == AtomId::new(0)
    ));
    assert_eq!(input, snapshot);
    assert_eq!(
        kekulize_if_possible(&input, &KekulizeParams::default()).unwrap(),
        KekulizeAttempt::NotKekulizable {
            topology: snapshot.clone(),
            problem_atoms: vec![AtomId::new(0)],
        }
    );
    assert_eq!(input, snapshot);
}

#[test]
fn malformed_aromatic_query_and_topology_errors_are_exact_and_not_swallowed() {
    let inconsistent = topology(
        vec![
            atom(0, AtomSpec::new(Element::C)),
            atom(1, AtomSpec::new(Element::C)),
        ],
        vec![bond(
            0,
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Aromatic),
        )],
    );
    // Source accepts independent aromatic order/false flags; preserve the
    // original literal two-atom input and both whole-entry invocations.
    let snapshot = inconsistent.clone();
    let result = kekulize(&inconsistent, &KekulizeParams::default());
    assert_eq!(inconsistent, snapshot);
    assert_eq!(result.unwrap().topology, snapshot);
    let result = kekulize_if_possible(&inconsistent, &KekulizeParams::default());
    assert_eq!(inconsistent, snapshot);
    assert!(
        matches!(result.unwrap(), KekulizeAttempt::Applied(assignment) if assignment.topology == snapshot)
    );

    let compound_query = topology(
        vec![
            atom(0, AtomSpec::new(Element::C)),
            atom(1, AtomSpec::new(Element::C)),
        ],
        vec![bond(
            0,
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        )],
    );
    let compound_atoms = compound_query
        .atoms
        .iter()
        .map(|carrier| {
            QueryAtom::from_carrier_parts(
                carrier.clone(),
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
            )
        })
        .collect::<Vec<_>>();
    let compound_bonds = vec![QueryBond::from_parts(
        compound_query.bonds[0].clone(),
        QueryNode::and(vec![
            QueryNode::predicate(BondQueryPredicate::Any),
            QueryNode::not(QueryNode::predicate(BondQueryPredicate::IsInRing(true))),
        ]),
    )];
    let compound_state =
        QueryStateRef::try_for_topology(&compound_atoms, &compound_bonds, &compound_query).unwrap();
    assert_eq!(
        kekulize_with_query_state(
            &compound_query,
            &KekulizeParams::default(),
            Some(compound_state),
        )
        .unwrap()
        .topology,
        compound_query
    );
    assert!(matches!(
        kekulize_if_possible_with_query_state(
            &compound_query,
            &KekulizeParams::default(),
            Some(compound_state),
        )
        .unwrap(),
        KekulizeAttempt::Applied(assignment) if assignment.topology == compound_query
    ));

    let invalid = TopologyBlock {
        atoms: vec![atom(1, AtomSpec::new(Element::C))],
        bonds: vec![],
        adjacency: AdjacencyList::try_from_topology(1, &[]).unwrap(),
        substance_groups: vec![],
        stereo_groups: vec![],
    };
    assert!(matches!(
        kekulize_if_possible(&invalid, &KekulizeParams::default()),
        Err(KekulizeError::InvalidTopology(TopologyValidationError::AtomIdMismatch {
            position: 0,
            id,
        })) if id == AtomId::new(1)
    ));
}

#[test]
fn source_issue_fixtures_cover_pyrrolic_h_charged_rings_dummy_and_zero_bond() {
    for element in [Element::N, Element::P] {
        let input = aromatic_cycle(
            std::iter::once(
                AtomSpec::new(element)
                    .with_explicit_hydrogens(1)
                    .with_no_implicit(true),
            )
            .chain((0..4).map(|_| AtomSpec::new(Element::C)))
            .collect(),
        );
        let output = kekulize(&input, &KekulizeParams::default())
            .unwrap()
            .topology;
        assert_eq!(output.atoms[0].explicit_hydrogens(), 0);
        assert!(!output.atoms[0].no_implicit());
        assert!(!output.atoms[0].is_aromatic());
        assert_kekule_ring(&output, 5, 2);
    }

    let charged_boron = aromatic_cycle(
        std::iter::once(AtomSpec::new(Element::B).with_formal_charge(-1))
            .chain((0..5).map(|_| AtomSpec::new(Element::C)))
            .collect(),
    );
    assert_kekule_ring(
        &kekulize(&charged_boron, &KekulizeParams::default())
            .unwrap()
            .topology,
        6,
        3,
    );

    let charged_nitrogen = aromatic_cycle(
        std::iter::once(AtomSpec::new(Element::N).with_formal_charge(1))
            .chain((0..5).map(|_| AtomSpec::new(Element::C)))
            .collect(),
    );
    assert_kekule_ring(
        &kekulize(&charged_nitrogen, &KekulizeParams::default())
            .unwrap()
            .topology,
        6,
        3,
    );

    let cyclopropenyl = aromatic_cycle(vec![
        AtomSpec::new(Element::C)
            .with_formal_charge(1)
            .with_explicit_hydrogens(1)
            .with_no_implicit(true),
        AtomSpec::new(Element::C),
        AtomSpec::new(Element::C),
    ]);
    assert_kekule_ring(
        &kekulize(&cyclopropenyl, &KekulizeParams::default())
            .unwrap()
            .topology,
        3,
        1,
    );

    let aromatic_dummy = topology(
        vec![
            atom(0, AtomSpec::new(Element::DUMMY)),
            atom(1, AtomSpec::new(Element::C).with_aromatic(true)),
            atom(2, AtomSpec::new(Element::N).with_aromatic(true)),
            atom(3, AtomSpec::new(Element::C).with_aromatic(true)),
            atom(4, AtomSpec::new(Element::C).with_aromatic(true)),
        ],
        vec![
            bond(
                0,
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            ),
            aromatic_bond(1, 1, 2),
            aromatic_bond(2, 2, 3),
            aromatic_bond(3, 3, 4),
            bond(
                4,
                BondSpec::new(AtomId::new(4), AtomId::new(0), BondOrder::Single),
            ),
        ],
    );
    let dummy_output = kekulize(&aromatic_dummy, &KekulizeParams::default())
        .unwrap()
        .topology;
    assert_eq!(dummy_output.bonds[0].order(), BondOrder::Single);
    assert_eq!(dummy_output.bonds[4].order(), BondOrder::Single);

    let mut zero_bond = aromatic_cycle(vec![
        AtomSpec::new(Element::C),
        AtomSpec::new(Element::C),
        AtomSpec::new(Element::N),
        AtomSpec::new(Element::C),
        AtomSpec::new(Element::N),
        AtomSpec::new(Element::C),
    ]);
    zero_bond.atoms.push(atom(6, AtomSpec::new(Element::FE)));
    zero_bond.bonds.push(bond(
        6,
        BondSpec::new(AtomId::new(5), AtomId::new(6), BondOrder::Zero),
    ));
    zero_bond.adjacency =
        AdjacencyList::try_from_topology(zero_bond.atoms.len(), &zero_bond.bonds).unwrap();
    let zero_output = kekulize(&zero_bond, &KekulizeParams::default())
        .unwrap()
        .topology;
    assert_kekule_ring(&zero_output, 6, 3);
    assert_eq!(zero_output.bonds[6].order(), BondOrder::Zero);
    assert!(!zero_output.bonds[6].is_aromatic());
}

#[test]
fn source_bodies_and_reviewed_two_axis_markers_are_local_to_the_owner() {
    let source = include_str!("../src/kekulize.rs");
    for function in [
        "Kekulize.cpp :: backTrack",
        "Kekulize.cpp :: markDbondCands",
        "Kekulize.cpp :: kekulizeWorker",
        "Kekulize.cpp :: QuestionEnumerator::QuestionEnumerator",
        "Kekulize.cpp :: QuestionEnumerator::next",
        "Kekulize.cpp :: permuteDummiesAndKekulize",
        "Kekulize.cpp :: kekulizeFused",
        "Kekulize.cpp :: KekulizeFragment (2026.03.6 complete)",
        "Kekulize.cpp :: MolOps::Kekulize",
        "Kekulize.cpp :: MolOps::KekulizeIfPossible",
    ] {
        assert!(
            source.contains(&format!(
                "BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/{function}"
            )),
            "missing source body anchor for {function}"
        );
    }
    assert!(source.matches("RDKit✔️✔️:").count() > 400);
    // Source-backed review keeps detached validation/copy/Vec<bool> costs
    // qualified; a blanket ban on qualifications would hide those real costs.
    // Full function anchors and behavioral tests remain separate requirements.
    assert!(source.contains("KekulizeFragment (2026.03.6 complete)"));
    assert!(source.contains("rankFragmentAtoms (2026.03.6 complete)"));
    assert!(source.contains("RDKit❗❌:"));
    assert!(!source.contains("RDKit❌"));
}

#[test]
fn q01_b1_canonical_vector_rank_is_lazy_and_named() {
    use cosmolkit_core::{
        CanonicalRankError, CanonicalRankParams, rank_fragment_atoms_with_params,
        rank_mol_atoms_with_params,
    };
    use cosmolkit_model::PropertyValueKind;
    let mut t = topology(
        vec![
            atom(0, AtomSpec::new(Element::C)),
            atom(1, AtomSpec::new(Element::C)),
        ],
        vec![],
    );
    let control = rank_mol_atoms_with_params(&t, &CanonicalRankParams::default()).unwrap();
    t.atoms[0]
        .set_prop("_CanonicalRankingNumber", vec![1_i32])
        .unwrap();
    assert_eq!(
        rank_mol_atoms_with_params(&t, &CanonicalRankParams::default()).unwrap(),
        control
    );
    let mut params = CanonicalRankParams::default();
    params.use_non_stereo_ranks = true;
    let before = t.clone();
    assert_eq!(
        rank_mol_atoms_with_params(&t, &params),
        Err(CanonicalRankError::InvalidPropertyKind {
            atom_index: 0,
            property: "_CanonicalRankingNumber",
            kind: PropertyValueKind::IntVector
        })
    );
    assert_eq!(t, before);
    t.atoms[0].clear_prop("_CanonicalRankingNumber");
    t.atoms[1]
        .set_prop("_CanonicalRankingNumber", vec![-2_i32])
        .unwrap();
    assert_eq!(
        rank_mol_atoms_with_params(&t, &params),
        Err(CanonicalRankError::InvalidPropertyKind {
            atom_index: 1,
            property: "_CanonicalRankingNumber",
            kind: PropertyValueKind::IntVector
        })
    );
    assert!(rank_fragment_atoms_with_params(&t, &[true, true], &[], None, None, &params).is_ok());
    let singleton = topology(
        vec![atom(
            0,
            AtomSpec::new(Element::C)
                .with_prop("_CanonicalRankingNumber", vec![1_i32])
                .unwrap(),
        )],
        vec![],
    );
    assert!(rank_mol_atoms_with_params(&singleton, &params).is_ok());
}

fn q01_b1_scalar_rank_source_cases() -> Vec<(
    &'static str,
    Option<PropertyValue>,
    Result<i32, cosmolkit_model::PropertyValueKind>,
)> {
    use cosmolkit_model::PropertyValueKind;
    vec![
        (
            "S000",
            Some(PropertyValue::String("0".to_owned().into())),
            Ok(0),
        ),
        (
            "S001",
            Some(PropertyValue::String("-0".to_owned().into())),
            Ok(0),
        ),
        (
            "S002",
            Some(PropertyValue::String("+0".to_owned().into())),
            Ok(0),
        ),
        (
            "S003",
            Some(PropertyValue::String("1".to_owned().into())),
            Ok(1),
        ),
        (
            "S004",
            Some(PropertyValue::String("-1".to_owned().into())),
            Ok(-1),
        ),
        (
            "S005",
            Some(PropertyValue::String("+1".to_owned().into())),
            Ok(1),
        ),
        (
            "S006",
            Some(PropertyValue::String("001".to_owned().into())),
            Ok(1),
        ),
        (
            "S007",
            Some(PropertyValue::String("-001".to_owned().into())),
            Ok(-1),
        ),
        (
            "S008",
            Some(PropertyValue::String("+001".to_owned().into())),
            Ok(1),
        ),
        (
            "S009",
            Some(PropertyValue::String("2147483647".to_owned().into())),
            Ok(2147483647),
        ),
        (
            "S010",
            Some(PropertyValue::String("-2147483648".to_owned().into())),
            Ok(-2147483648),
        ),
        (
            "S011",
            Some(PropertyValue::String("+2147483647".to_owned().into())),
            Ok(2147483647),
        ),
        (
            "S012",
            Some(PropertyValue::String("2147483648".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S013",
            Some(PropertyValue::String("-2147483649".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S014",
            Some(PropertyValue::String("+2147483648".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S015",
            Some(PropertyValue::String(
                ("999999999999999999999999999999999999".to_owned()).into(),
            )),
            Err(PropertyValueKind::String),
        ),
        (
            "S016",
            Some(PropertyValue::String(
                ("000000000000000000000000000000000000000000001".to_owned()).into(),
            )),
            Ok(1),
        ),
        (
            "S017",
            Some(PropertyValue::String(
                ("-000000000000000000000000000000002147483648".to_owned()).into(),
            )),
            Ok(-2147483648),
        ),
        (
            "S018",
            Some(PropertyValue::String("7 ".to_owned().into())),
            Ok(7),
        ),
        (
            "S019",
            Some(PropertyValue::String("7\t".to_owned().into())),
            Ok(7),
        ),
        (
            "S020",
            Some(PropertyValue::String("7\n".to_owned().into())),
            Ok(7),
        ),
        (
            "S021",
            Some(PropertyValue::String("7\r".to_owned().into())),
            Ok(7),
        ),
        (
            "S022",
            Some(PropertyValue::String("7\u{000b}".to_owned().into())),
            Ok(7),
        ),
        (
            "S023",
            Some(PropertyValue::String("7\u{000c}".to_owned().into())),
            Ok(7),
        ),
        (
            "S024",
            Some(PropertyValue::String(
                ("+7 \t\n\r\u{000b}\u{000c}".to_owned()).into(),
            )),
            Ok(7),
        ),
        (
            "S025",
            Some(PropertyValue::String(" 7".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S026",
            Some(PropertyValue::String("\t7".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S027",
            Some(PropertyValue::String("\n7".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S028",
            Some(PropertyValue::String("\r7".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S029",
            Some(PropertyValue::String("\u{000b}7".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S030",
            Some(PropertyValue::String("\u{000c}7".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S031",
            Some(PropertyValue::String(" 7 ".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S032",
            Some(PropertyValue::String("".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S033",
            Some(PropertyValue::String(" ".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S034",
            Some(PropertyValue::String(
                "\t\n\r\u{000b}\u{000c}".to_owned().into(),
            )),
            Err(PropertyValueKind::String),
        ),
        (
            "S035",
            Some(PropertyValue::String("+".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S036",
            Some(PropertyValue::String("-".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S037",
            Some(PropertyValue::String("+-1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S038",
            Some(PropertyValue::String("--1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S039",
            Some(PropertyValue::String("++1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S040",
            Some(PropertyValue::String("1-".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S041",
            Some(PropertyValue::String("1+".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S042",
            Some(PropertyValue::String("1 2".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S043",
            Some(PropertyValue::String("1\t2".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S044",
            Some(PropertyValue::String("1.0".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S045",
            Some(PropertyValue::String("-1.0".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S046",
            Some(PropertyValue::String("1e0".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S047",
            Some(PropertyValue::String("0x10".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S048",
            Some(PropertyValue::String("0b10".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S049",
            Some(PropertyValue::String("1,000".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S050",
            Some(PropertyValue::String("true".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S051",
            Some(PropertyValue::String("false".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S052",
            Some(PropertyValue::String("nan".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S053",
            Some(PropertyValue::String("inf".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S054",
            Some(PropertyValue::String("[1]".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S055",
            Some(PropertyValue::String("[]".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S056",
            Some(PropertyValue::String("1x".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S057",
            Some(PropertyValue::String("1\u{0000}".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S058",
            Some(PropertyValue::String("1\u{0000} ".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S059",
            Some(PropertyValue::String("\u{0000}1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S060",
            Some(PropertyValue::String("1 ".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S061",
            Some(PropertyValue::String(" 1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S062",
            Some(PropertyValue::String("1 ".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S063",
            Some(PropertyValue::String(" 1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S064",
            Some(PropertyValue::String("1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S065",
            Some(PropertyValue::String("1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S066",
            Some(PropertyValue::String("١".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S067",
            Some(PropertyValue::String("１".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S068",
            Some(PropertyValue::String("²".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "S069",
            Some(PropertyValue::String("−1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "I070",
            Some(PropertyValue::Int(-2147483648)),
            Ok(-2147483648),
        ),
        ("I071", Some(PropertyValue::Int(-2)), Ok(-2)),
        ("I072", Some(PropertyValue::Int(-1)), Ok(-1)),
        ("I073", Some(PropertyValue::Int(0)), Ok(0)),
        ("I074", Some(PropertyValue::Int(1)), Ok(1)),
        ("I075", Some(PropertyValue::Int(2)), Ok(2)),
        ("I076", Some(PropertyValue::Int(2147483647)), Ok(2147483647)),
        (
            "B077",
            Some(PropertyValue::Bool(false)),
            Err(PropertyValueKind::Bool),
        ),
        (
            "B078",
            Some(PropertyValue::Bool(true)),
            Err(PropertyValueKind::Bool),
        ),
        (
            "D079",
            Some(PropertyValue::Double(f64::NEG_INFINITY)),
            Err(PropertyValueKind::Double),
        ),
        (
            "D080",
            Some(PropertyValue::Double(-2147483648.0)),
            Err(PropertyValueKind::Double),
        ),
        (
            "D081",
            Some(PropertyValue::Double(-1.5)),
            Err(PropertyValueKind::Double),
        ),
        (
            "D082",
            Some(PropertyValue::Double(-1.0)),
            Err(PropertyValueKind::Double),
        ),
        (
            "D083",
            Some(PropertyValue::Double(-0.0)),
            Err(PropertyValueKind::Double),
        ),
        (
            "D084",
            Some(PropertyValue::Double(0.0)),
            Err(PropertyValueKind::Double),
        ),
        (
            "D085",
            Some(PropertyValue::Double(1.0)),
            Err(PropertyValueKind::Double),
        ),
        (
            "D086",
            Some(PropertyValue::Double(1.5)),
            Err(PropertyValueKind::Double),
        ),
        (
            "D087",
            Some(PropertyValue::Double(2147483647.0)),
            Err(PropertyValueKind::Double),
        ),
        (
            "D088",
            Some(PropertyValue::Double(f64::INFINITY)),
            Err(PropertyValueKind::Double),
        ),
        (
            "D089",
            Some(PropertyValue::Double(f64::NAN)),
            Err(PropertyValueKind::Double),
        ),
        (
            "AL000",
            Some(PropertyValue::String("\u{0000}1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT000",
            Some(PropertyValue::String("1\u{0000}".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL001",
            Some(PropertyValue::String("\u{0001}1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT001",
            Some(PropertyValue::String("1\u{0001}".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL002",
            Some(PropertyValue::String("\u{0002}1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT002",
            Some(PropertyValue::String("1\u{0002}".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL003",
            Some(PropertyValue::String("\u{0003}1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT003",
            Some(PropertyValue::String("1\u{0003}".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL004",
            Some(PropertyValue::String("\u{0004}1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT004",
            Some(PropertyValue::String("1\u{0004}".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL005",
            Some(PropertyValue::String("\u{0005}1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT005",
            Some(PropertyValue::String("1\u{0005}".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL006",
            Some(PropertyValue::String("\u{0006}1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT006",
            Some(PropertyValue::String("1\u{0006}".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL007",
            Some(PropertyValue::String("\u{0007}1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT007",
            Some(PropertyValue::String("1\u{0007}".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL008",
            Some(PropertyValue::String("\u{0008}1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT008",
            Some(PropertyValue::String("1\u{0008}".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL009",
            Some(PropertyValue::String("\t1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT009",
            Some(PropertyValue::String("1\t".to_owned().into())),
            Ok(1),
        ),
        (
            "AL010",
            Some(PropertyValue::String("\n1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT010",
            Some(PropertyValue::String("1\n".to_owned().into())),
            Ok(1),
        ),
        (
            "AL011",
            Some(PropertyValue::String("\u{000b}1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT011",
            Some(PropertyValue::String("1\u{000b}".to_owned().into())),
            Ok(1),
        ),
        (
            "AL012",
            Some(PropertyValue::String("\u{000c}1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT012",
            Some(PropertyValue::String("1\u{000c}".to_owned().into())),
            Ok(1),
        ),
        (
            "AL013",
            Some(PropertyValue::String("\r1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT013",
            Some(PropertyValue::String("1\r".to_owned().into())),
            Ok(1),
        ),
        (
            "AL014",
            Some(PropertyValue::String("\u{000e}1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT014",
            Some(PropertyValue::String("1\u{000e}".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL015",
            Some(PropertyValue::String("\u{000f}1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT015",
            Some(PropertyValue::String("1\u{000f}".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL016",
            Some(PropertyValue::String("\u{0010}1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT016",
            Some(PropertyValue::String("1\u{0010}".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL017",
            Some(PropertyValue::String("\u{0011}1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT017",
            Some(PropertyValue::String("1\u{0011}".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL018",
            Some(PropertyValue::String("\u{0012}1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT018",
            Some(PropertyValue::String("1\u{0012}".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL019",
            Some(PropertyValue::String("\u{0013}1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT019",
            Some(PropertyValue::String("1\u{0013}".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL020",
            Some(PropertyValue::String("\u{0014}1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT020",
            Some(PropertyValue::String("1\u{0014}".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL021",
            Some(PropertyValue::String("\u{0015}1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT021",
            Some(PropertyValue::String("1\u{0015}".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL022",
            Some(PropertyValue::String("\u{0016}1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT022",
            Some(PropertyValue::String("1\u{0016}".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL023",
            Some(PropertyValue::String("\u{0017}1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT023",
            Some(PropertyValue::String("1\u{0017}".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL024",
            Some(PropertyValue::String("\u{0018}1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT024",
            Some(PropertyValue::String("1\u{0018}".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL025",
            Some(PropertyValue::String("\u{0019}1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT025",
            Some(PropertyValue::String("1\u{0019}".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL026",
            Some(PropertyValue::String("\u{001a}1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT026",
            Some(PropertyValue::String("1\u{001a}".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL027",
            Some(PropertyValue::String("\u{001b}1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT027",
            Some(PropertyValue::String("1\u{001b}".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL028",
            Some(PropertyValue::String("\u{001c}1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT028",
            Some(PropertyValue::String("1\u{001c}".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL029",
            Some(PropertyValue::String("\u{001d}1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT029",
            Some(PropertyValue::String("1\u{001d}".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL030",
            Some(PropertyValue::String("\u{001e}1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT030",
            Some(PropertyValue::String("1\u{001e}".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL031",
            Some(PropertyValue::String("\u{001f}1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT031",
            Some(PropertyValue::String("1\u{001f}".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL032",
            Some(PropertyValue::String(" 1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT032",
            Some(PropertyValue::String("1 ".to_owned().into())),
            Ok(1),
        ),
        (
            "AL033",
            Some(PropertyValue::String("!1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT033",
            Some(PropertyValue::String("1!".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL034",
            Some(PropertyValue::String("\"1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT034",
            Some(PropertyValue::String("1\"".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL035",
            Some(PropertyValue::String("#1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT035",
            Some(PropertyValue::String("1#".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL036",
            Some(PropertyValue::String("$1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT036",
            Some(PropertyValue::String("1$".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL037",
            Some(PropertyValue::String("%1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT037",
            Some(PropertyValue::String("1%".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL038",
            Some(PropertyValue::String("&1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT038",
            Some(PropertyValue::String("1&".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL039",
            Some(PropertyValue::String("'1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT039",
            Some(PropertyValue::String("1'".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL040",
            Some(PropertyValue::String("(1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT040",
            Some(PropertyValue::String("1(".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL041",
            Some(PropertyValue::String(")1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT041",
            Some(PropertyValue::String("1)".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL042",
            Some(PropertyValue::String("*1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT042",
            Some(PropertyValue::String("1*".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL043",
            Some(PropertyValue::String("+1".to_owned().into())),
            Ok(1),
        ),
        (
            "AT043",
            Some(PropertyValue::String("1+".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL044",
            Some(PropertyValue::String(",1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT044",
            Some(PropertyValue::String("1,".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL045",
            Some(PropertyValue::String("-1".to_owned().into())),
            Ok(-1),
        ),
        (
            "AT045",
            Some(PropertyValue::String("1-".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL046",
            Some(PropertyValue::String(".1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT046",
            Some(PropertyValue::String("1.".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL047",
            Some(PropertyValue::String("/1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT047",
            Some(PropertyValue::String("1/".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL048",
            Some(PropertyValue::String("01".to_owned().into())),
            Ok(1),
        ),
        (
            "AT048",
            Some(PropertyValue::String("10".to_owned().into())),
            Ok(10),
        ),
        (
            "AL049",
            Some(PropertyValue::String("11".to_owned().into())),
            Ok(11),
        ),
        (
            "AT049",
            Some(PropertyValue::String("11".to_owned().into())),
            Ok(11),
        ),
        (
            "AL050",
            Some(PropertyValue::String("21".to_owned().into())),
            Ok(21),
        ),
        (
            "AT050",
            Some(PropertyValue::String("12".to_owned().into())),
            Ok(12),
        ),
        (
            "AL051",
            Some(PropertyValue::String("31".to_owned().into())),
            Ok(31),
        ),
        (
            "AT051",
            Some(PropertyValue::String("13".to_owned().into())),
            Ok(13),
        ),
        (
            "AL052",
            Some(PropertyValue::String("41".to_owned().into())),
            Ok(41),
        ),
        (
            "AT052",
            Some(PropertyValue::String("14".to_owned().into())),
            Ok(14),
        ),
        (
            "AL053",
            Some(PropertyValue::String("51".to_owned().into())),
            Ok(51),
        ),
        (
            "AT053",
            Some(PropertyValue::String("15".to_owned().into())),
            Ok(15),
        ),
        (
            "AL054",
            Some(PropertyValue::String("61".to_owned().into())),
            Ok(61),
        ),
        (
            "AT054",
            Some(PropertyValue::String("16".to_owned().into())),
            Ok(16),
        ),
        (
            "AL055",
            Some(PropertyValue::String("71".to_owned().into())),
            Ok(71),
        ),
        (
            "AT055",
            Some(PropertyValue::String("17".to_owned().into())),
            Ok(17),
        ),
        (
            "AL056",
            Some(PropertyValue::String("81".to_owned().into())),
            Ok(81),
        ),
        (
            "AT056",
            Some(PropertyValue::String("18".to_owned().into())),
            Ok(18),
        ),
        (
            "AL057",
            Some(PropertyValue::String("91".to_owned().into())),
            Ok(91),
        ),
        (
            "AT057",
            Some(PropertyValue::String("19".to_owned().into())),
            Ok(19),
        ),
        (
            "AL058",
            Some(PropertyValue::String(":1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT058",
            Some(PropertyValue::String("1:".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL059",
            Some(PropertyValue::String(";1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT059",
            Some(PropertyValue::String("1;".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL060",
            Some(PropertyValue::String("<1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT060",
            Some(PropertyValue::String("1<".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL061",
            Some(PropertyValue::String("=1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT061",
            Some(PropertyValue::String("1=".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL062",
            Some(PropertyValue::String(">1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT062",
            Some(PropertyValue::String("1>".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL063",
            Some(PropertyValue::String("?1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT063",
            Some(PropertyValue::String("1?".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL064",
            Some(PropertyValue::String("@1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT064",
            Some(PropertyValue::String("1@".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL065",
            Some(PropertyValue::String("A1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT065",
            Some(PropertyValue::String("1A".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL066",
            Some(PropertyValue::String("B1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT066",
            Some(PropertyValue::String("1B".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL067",
            Some(PropertyValue::String("C1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT067",
            Some(PropertyValue::String("1C".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL068",
            Some(PropertyValue::String("D1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT068",
            Some(PropertyValue::String("1D".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL069",
            Some(PropertyValue::String("E1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT069",
            Some(PropertyValue::String("1E".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL070",
            Some(PropertyValue::String("F1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT070",
            Some(PropertyValue::String("1F".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL071",
            Some(PropertyValue::String("G1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT071",
            Some(PropertyValue::String("1G".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL072",
            Some(PropertyValue::String("H1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT072",
            Some(PropertyValue::String("1H".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL073",
            Some(PropertyValue::String("I1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT073",
            Some(PropertyValue::String("1I".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL074",
            Some(PropertyValue::String("J1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT074",
            Some(PropertyValue::String("1J".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL075",
            Some(PropertyValue::String("K1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT075",
            Some(PropertyValue::String("1K".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL076",
            Some(PropertyValue::String("L1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT076",
            Some(PropertyValue::String("1L".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL077",
            Some(PropertyValue::String("M1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT077",
            Some(PropertyValue::String("1M".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL078",
            Some(PropertyValue::String("N1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT078",
            Some(PropertyValue::String("1N".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL079",
            Some(PropertyValue::String("O1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT079",
            Some(PropertyValue::String("1O".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL080",
            Some(PropertyValue::String("P1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT080",
            Some(PropertyValue::String("1P".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL081",
            Some(PropertyValue::String("Q1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT081",
            Some(PropertyValue::String("1Q".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL082",
            Some(PropertyValue::String("R1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT082",
            Some(PropertyValue::String("1R".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL083",
            Some(PropertyValue::String("S1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT083",
            Some(PropertyValue::String("1S".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL084",
            Some(PropertyValue::String("T1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT084",
            Some(PropertyValue::String("1T".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL085",
            Some(PropertyValue::String("U1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT085",
            Some(PropertyValue::String("1U".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL086",
            Some(PropertyValue::String("V1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT086",
            Some(PropertyValue::String("1V".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL087",
            Some(PropertyValue::String("W1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT087",
            Some(PropertyValue::String("1W".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL088",
            Some(PropertyValue::String("X1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT088",
            Some(PropertyValue::String("1X".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL089",
            Some(PropertyValue::String("Y1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT089",
            Some(PropertyValue::String("1Y".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL090",
            Some(PropertyValue::String("Z1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT090",
            Some(PropertyValue::String("1Z".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL091",
            Some(PropertyValue::String("[1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT091",
            Some(PropertyValue::String("1[".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL092",
            Some(PropertyValue::String("\\1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT092",
            Some(PropertyValue::String("1\\".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL093",
            Some(PropertyValue::String("]1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT093",
            Some(PropertyValue::String("1]".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL094",
            Some(PropertyValue::String("^1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT094",
            Some(PropertyValue::String("1^".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL095",
            Some(PropertyValue::String("_1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT095",
            Some(PropertyValue::String("1_".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL096",
            Some(PropertyValue::String("`1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT096",
            Some(PropertyValue::String("1`".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL097",
            Some(PropertyValue::String("a1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT097",
            Some(PropertyValue::String("1a".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL098",
            Some(PropertyValue::String("b1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT098",
            Some(PropertyValue::String("1b".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL099",
            Some(PropertyValue::String("c1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT099",
            Some(PropertyValue::String("1c".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL100",
            Some(PropertyValue::String("d1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT100",
            Some(PropertyValue::String("1d".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL101",
            Some(PropertyValue::String("e1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT101",
            Some(PropertyValue::String("1e".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL102",
            Some(PropertyValue::String("f1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT102",
            Some(PropertyValue::String("1f".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL103",
            Some(PropertyValue::String("g1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT103",
            Some(PropertyValue::String("1g".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL104",
            Some(PropertyValue::String("h1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT104",
            Some(PropertyValue::String("1h".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL105",
            Some(PropertyValue::String("i1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT105",
            Some(PropertyValue::String("1i".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL106",
            Some(PropertyValue::String("j1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT106",
            Some(PropertyValue::String("1j".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL107",
            Some(PropertyValue::String("k1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT107",
            Some(PropertyValue::String("1k".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL108",
            Some(PropertyValue::String("l1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT108",
            Some(PropertyValue::String("1l".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL109",
            Some(PropertyValue::String("m1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT109",
            Some(PropertyValue::String("1m".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL110",
            Some(PropertyValue::String("n1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT110",
            Some(PropertyValue::String("1n".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL111",
            Some(PropertyValue::String("o1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT111",
            Some(PropertyValue::String("1o".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL112",
            Some(PropertyValue::String("p1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT112",
            Some(PropertyValue::String("1p".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL113",
            Some(PropertyValue::String("q1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT113",
            Some(PropertyValue::String("1q".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL114",
            Some(PropertyValue::String("r1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT114",
            Some(PropertyValue::String("1r".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL115",
            Some(PropertyValue::String("s1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT115",
            Some(PropertyValue::String("1s".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL116",
            Some(PropertyValue::String("t1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT116",
            Some(PropertyValue::String("1t".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL117",
            Some(PropertyValue::String("u1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT117",
            Some(PropertyValue::String("1u".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL118",
            Some(PropertyValue::String("v1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT118",
            Some(PropertyValue::String("1v".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL119",
            Some(PropertyValue::String("w1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT119",
            Some(PropertyValue::String("1w".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL120",
            Some(PropertyValue::String("x1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT120",
            Some(PropertyValue::String("1x".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL121",
            Some(PropertyValue::String("y1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT121",
            Some(PropertyValue::String("1y".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL122",
            Some(PropertyValue::String("z1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT122",
            Some(PropertyValue::String("1z".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL123",
            Some(PropertyValue::String("{1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT123",
            Some(PropertyValue::String("1{".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL124",
            Some(PropertyValue::String("|1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT124",
            Some(PropertyValue::String("1|".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL125",
            Some(PropertyValue::String("}1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT125",
            Some(PropertyValue::String("1}".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL126",
            Some(PropertyValue::String("~1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT126",
            Some(PropertyValue::String("1~".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AL127",
            Some(PropertyValue::String("1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "AT127",
            Some(PropertyValue::String("1".to_owned().into())),
            Err(PropertyValueKind::String),
        ),
        (
            "LZERO",
            Some(PropertyValue::String(
                (format!("{}1", "0".repeat(65_536))).into(),
            )),
            Ok(1),
        ),
        (
            "LOVF",
            Some(PropertyValue::String(("9".repeat(65_536)).into())),
            Err(PropertyValueKind::String),
        ),
        (
            "VEMPTY",
            Some(PropertyValue::IntVector(vec![])),
            Err(PropertyValueKind::IntVector),
        ),
        (
            "VSINGLE",
            Some(PropertyValue::IntVector(vec![1])),
            Err(PropertyValueKind::IntVector),
        ),
        (
            "VORDER",
            Some(PropertyValue::IntVector(vec![1, -2, 1])),
            Err(PropertyValueKind::IntVector),
        ),
        (
            "VBOUNDS",
            Some(PropertyValue::IntVector(vec![-2147483648, 2147483647])),
            Err(PropertyValueKind::IntVector),
        ),
        ("MISSING", None, Ok(0)),
    ]
}

#[test]
fn q01_b1_scalar_rank_casts_cover_full_source_integer_domain() {
    use cosmolkit_core::{CanonicalRankError, CanonicalRankParams, rank_mol_atoms_with_params};
    let mut params = CanonicalRankParams::default();
    params.break_ties = false;
    params.use_non_stereo_ranks = true;
    for (label, value, integer) in q01_b1_scalar_rank_source_cases() {
        for side in 0..2 {
            let mut graph = topology(
                vec![
                    atom(0, AtomSpec::new(Element::C)),
                    atom(1, AtomSpec::new(Element::C)),
                ],
                vec![],
            );
            graph.atoms[0]
                .set_computed_prop("q01_unrelated", vec![1_i32, -2, 1])
                .unwrap();
            if let Some(value) = &value {
                graph.atoms[side]
                    .set_computed_prop("_CanonicalRankingNumber", value.clone())
                    .unwrap();
            }
            match integer {
                Ok(number) => {
                    // The other atom receives the independently frozen exact integer;
                    // equal source classes prove the entire value, not merely its sign.
                    graph.atoms[1 - side]
                        .set_prop("_CanonicalRankingNumber", number)
                        .unwrap();
                    let before = graph.clone();
                    assert_eq!(
                        rank_mol_atoms_with_params(&graph, &params),
                        Ok(vec![0, 0]),
                        "{label}, side{side}"
                    );
                    assert_eq!(graph, before, "{label} input/computed metadata");
                    if let Some(lower) = number.checked_sub(1) {
                        graph.atoms[1 - side]
                            .set_prop("_CanonicalRankingNumber", lower)
                            .unwrap();
                        let before = graph.clone();
                        let expected = if side == 0 { vec![1, 0] } else { vec![0, 1] };
                        assert_eq!(
                            rank_mol_atoms_with_params(&graph, &params),
                            Ok(expected),
                            "{label} lower, side{side}"
                        );
                        assert_eq!(graph, before);
                    }
                    if let Some(upper) = number.checked_add(1) {
                        graph.atoms[1 - side]
                            .set_prop("_CanonicalRankingNumber", upper)
                            .unwrap();
                        let before = graph.clone();
                        let expected = if side == 0 { vec![0, 1] } else { vec![1, 0] };
                        assert_eq!(
                            rank_mol_atoms_with_params(&graph, &params),
                            Ok(expected),
                            "{label} upper, side{side}"
                        );
                        assert_eq!(graph, before);
                    }
                }
                Err(kind) => {
                    let before = graph.clone();
                    assert_eq!(
                        rank_mol_atoms_with_params(&graph, &params),
                        Err(CanonicalRankError::InvalidPropertyKind {
                            atom_index: side,
                            property: "_CanonicalRankingNumber",
                            kind
                        }),
                        "{label}, side{side}"
                    );
                    assert_eq!(graph, before, "{label} failed input/computed metadata");
                }
            }
        }
    }
}

#[test]
fn q01_b1_scalar_rank_error_order_follows_left_then_right_getters() {
    use cosmolkit_core::{CanonicalRankError, CanonicalRankParams, rank_mol_atoms_with_params};
    use cosmolkit_model::PropertyValueKind;
    let mut params = CanonicalRankParams::default();
    params.use_non_stereo_ranks = true;
    for (left, right, index, kind) in [
        (
            Some(PropertyValue::String("not-an-int".to_owned().into())),
            Some(PropertyValue::Double(1.0)),
            0,
            PropertyValueKind::String,
        ),
        (
            Some(PropertyValue::Int(i32::MAX)),
            Some(PropertyValue::Double(1.0)),
            1,
            PropertyValueKind::Double,
        ),
        (
            None,
            Some(PropertyValue::String("2147483648".to_owned().into())),
            1,
            PropertyValueKind::String,
        ),
        (
            Some(PropertyValue::Bool(false)),
            Some(PropertyValue::IntVector(vec![])),
            0,
            PropertyValueKind::Bool,
        ),
        (
            Some(PropertyValue::String("7 	".to_owned().into())),
            Some(PropertyValue::IntVector(vec![1])),
            1,
            PropertyValueKind::IntVector,
        ),
    ] {
        let mut graph = topology(
            vec![
                atom(0, AtomSpec::new(Element::C)),
                atom(1, AtomSpec::new(Element::C)),
            ],
            vec![],
        );
        for (index, value) in [left, right].into_iter().enumerate() {
            if let Some(value) = value {
                graph.atoms[index]
                    .set_prop("_CanonicalRankingNumber", value)
                    .unwrap();
            }
        }
        let before = graph.clone();
        assert_eq!(
            rank_mol_atoms_with_params(&graph, &params),
            Err(CanonicalRankError::InvalidPropertyKind {
                atom_index: index,
                property: "_CanonicalRankingNumber",
                kind
            })
        );
        assert_eq!(graph, before);
    }
}

#[test]
fn q01_b1_scalar_rank_no_read_guards_preserve_all_invalid_kinds() {
    use cosmolkit_core::{
        CanonicalRankError, CanonicalRankParams, ValenceParams, assign_valence,
        rank_fragment_atoms_with_params, rank_fragment_atoms_with_prepared_state,
        rank_mol_atoms_with_params,
    };
    let mut disabled = CanonicalRankParams::default();
    disabled.break_ties = false;
    let mut enabled = disabled;
    enabled.use_non_stereo_ranks = true;
    assert_eq!(
        rank_mol_atoms_with_params(&TopologyBlock::default(), &enabled),
        Ok(vec![])
    );
    for (label, value, integer) in q01_b1_scalar_rank_source_cases() {
        if integer.is_ok() {
            continue;
        }
        let value = value.unwrap();
        let mut graph = topology(
            vec![
                atom(0, AtomSpec::new(Element::C)),
                atom(1, AtomSpec::new(Element::C)),
            ],
            vec![],
        );
        graph.atoms[0]
            .set_computed_prop("_CanonicalRankingNumber", value.clone())
            .unwrap();
        graph.atoms[1]
            .set_prop("_CanonicalRankingNumber", value.clone())
            .unwrap();
        let before = graph.clone();
        assert_eq!(
            rank_mol_atoms_with_params(&graph, &disabled),
            Ok(vec![0, 0]),
            "{label} flagoff"
        );
        assert_eq!(
            rank_fragment_atoms_with_params(&graph, &[true, true], &[], None, None, &enabled),
            Ok(vec![0, 0]),
            "{label} fragment flag"
        );
        let valence = assign_valence(&graph, &ValenceParams::default()).unwrap();
        assert_eq!(
            rank_fragment_atoms_with_prepared_state(
                &graph,
                &valence,
                None,
                &[true, true],
                &[],
                None,
                None,
                &enabled
            ),
            Ok(vec![0, 0]),
            "{label} preparedfragment flag"
        );
        assert_eq!(
            rank_fragment_atoms_with_params(&graph, &[true], &[], None, None, &enabled),
            Err(CanonicalRankError::AtomMaskLength {
                expected: 2,
                actual: 1
            }),
            "{label} masks beforeprops"
        );
        let mut invalid_valence = valence.clone();
        invalid_valence.implicit_hydrogens.pop();
        assert_eq!(
            rank_fragment_atoms_with_prepared_state(
                &graph,
                &invalid_valence,
                None,
                &[true, true],
                &[],
                None,
                None,
                &enabled
            ),
            Err(CanonicalRankError::PreparedValenceLength {
                atom_count: 2,
                explicit_len: 2,
                implicit_len: 1
            }),
            "{label} preparedvalidation beforeprops"
        );
        assert_eq!(graph, before);
        let mut single = topology(vec![atom(0, AtomSpec::new(Element::C))], vec![]);
        single.atoms[0]
            .set_prop("_CanonicalRankingNumber", value)
            .unwrap();
        let before = single.clone();
        assert_eq!(
            rank_mol_atoms_with_params(&single, &enabled),
            Ok(vec![0]),
            "{label} singleton"
        );
        assert_eq!(single, before);
    }
}

#[test]
fn proposed_uint_rank_values_and_positive_overflow_remain_lazy() {
    use cosmolkit_core::{CanonicalRankError, CanonicalRankParams, rank_mol_atoms_with_params};
    let mut params = CanonicalRankParams::default();
    params.break_ties = false;
    params.use_non_stereo_ranks = true;
    for number in [
        0_u32,
        1,
        i32::MAX as u32 - 1,
        i32::MAX as u32,
        i32::MAX as u32 + 1,
        u32::MAX,
    ] {
        for side in 0..2 {
            let mut graph = topology(
                vec![
                    atom(0, AtomSpec::new(Element::C)),
                    atom(1, AtomSpec::new(Element::C)),
                ],
                vec![],
            );
            graph.atoms[side]
                .set_computed_prop("_CanonicalRankingNumber", PropertyValue::UInt(number))
                .unwrap();
            if let Ok(value) = i32::try_from(number) {
                graph.atoms[1 - side]
                    .set_prop("_CanonicalRankingNumber", value)
                    .unwrap();
            }
            let before = graph.clone();
            if number <= i32::MAX as u32 {
                assert_eq!(rank_mol_atoms_with_params(&graph, &params), Ok(vec![0, 0]));
            } else {
                assert_eq!(
                    rank_mol_atoms_with_params(&graph, &params),
                    Err(CanonicalRankError::UnsignedRankOverflow {
                        atom_index: side,
                        property: "_CanonicalRankingNumber",
                        value: number
                    })
                );
                let mut disabled = params;
                disabled.use_non_stereo_ranks = false;
                assert_eq!(
                    rank_mol_atoms_with_params(&graph, &disabled),
                    Ok(vec![0, 0])
                );
            }
            assert_eq!(graph, before);
        }
    }
}
#[test]
fn proposed_uint_public_singleton_fragment_prepared_mask_and_flag_guards() {
    use cosmolkit_core::{
        CanonicalRankError, CanonicalRankParams, ValenceParams, assign_valence,
        rank_fragment_atoms_with_params, rank_fragment_atoms_with_prepared_state,
        rank_mol_atoms_with_params,
    };
    let mut disabled = CanonicalRankParams::default();
    disabled.break_ties = false;
    let mut enabled = disabled;
    enabled.use_non_stereo_ranks = true;
    assert_eq!(
        rank_mol_atoms_with_params(&TopologyBlock::default(), &enabled),
        Ok(vec![])
    );
    for number in [0_u32, 1, 2147483646, 2147483647, 2147483648, 4294967295] {
        let label = number;
        let value = PropertyValue::UInt(number);
        let mut graph = topology(
            vec![
                atom(0, AtomSpec::new(Element::C)),
                atom(1, AtomSpec::new(Element::C)),
            ],
            vec![],
        );
        graph.atoms[0]
            .set_computed_prop("_CanonicalRankingNumber", value.clone())
            .unwrap();
        graph.atoms[1]
            .set_prop("_CanonicalRankingNumber", value.clone())
            .unwrap();
        let before = graph.clone();
        assert_eq!(
            rank_mol_atoms_with_params(&graph, &disabled),
            Ok(vec![0, 0]),
            "{label} flagoff"
        );
        assert_eq!(
            rank_fragment_atoms_with_params(&graph, &[true, true], &[], None, None, &enabled),
            Ok(vec![0, 0]),
            "{label} fragment flag"
        );
        let valence = assign_valence(&graph, &ValenceParams::default()).unwrap();
        assert_eq!(
            rank_fragment_atoms_with_prepared_state(
                &graph,
                &valence,
                None,
                &[true, true],
                &[],
                None,
                None,
                &enabled
            ),
            Ok(vec![0, 0]),
            "{label} preparedfragment flag"
        );
        assert_eq!(
            rank_fragment_atoms_with_params(&graph, &[true], &[], None, None, &enabled),
            Err(CanonicalRankError::AtomMaskLength {
                expected: 2,
                actual: 1
            }),
            "{label} masks beforeprops"
        );
        let mut invalid_valence = valence.clone();
        invalid_valence.implicit_hydrogens.pop();
        assert_eq!(
            rank_fragment_atoms_with_prepared_state(
                &graph,
                &invalid_valence,
                None,
                &[true, true],
                &[],
                None,
                None,
                &enabled
            ),
            Err(CanonicalRankError::PreparedValenceLength {
                atom_count: 2,
                explicit_len: 2,
                implicit_len: 1
            }),
            "{label} preparedvalidation beforeprops"
        );
        assert_eq!(graph, before);
        let mut single = topology(vec![atom(0, AtomSpec::new(Element::C))], vec![]);
        single.atoms[0]
            .set_prop("_CanonicalRankingNumber", value)
            .unwrap();
        let before = single.clone();
        assert_eq!(
            rank_mol_atoms_with_params(&single, &enabled),
            Ok(vec![0]),
            "{label} singleton"
        );
        assert_eq!(single, before);
    }
}
