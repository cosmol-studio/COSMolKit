use cosmolkit_core::{
    KekulizeAttempt, KekulizeError, KekulizeParams, kekulize, kekulize_if_possible,
    kekulize_if_possible_with_query_state, kekulize_with_query_state,
};
use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomQueryPredicate, AtomSpec, Bond, BondId, BondQueryPredicate,
    BondSpec, QueryAtom, QueryBond, QueryNode, QueryStateRef, TopologyBlock,
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
    assert_eq!(output.atoms[0].prop("atom-note"), Some("preserved"));
    assert_eq!(output.bonds[0].prop("ring-note"), Some("preserved"));
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
    for result in [
        kekulize(&inconsistent, &KekulizeParams::default()).map(|_| ()),
        kekulize_if_possible(&inconsistent, &KekulizeParams::default()).map(|_| ()),
    ] {
        assert!(matches!(
            result,
            Err(KekulizeError::AromaticBondStateMismatch {
                bond,
                order: BondOrder::Aromatic,
                is_aromatic: false,
            }) if bond == BondId::new(0)
        ));
    }

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
        "Kekulize.cpp :: KekulizeFragment selection",
        "Kekulize.cpp :: KekulizeFragment ranking and dispatch",
        "Kekulize.cpp :: KekulizeFragment finalization",
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
    assert!(!source.contains("RDKit❗"));
    assert!(!source.contains("RDKit❌"));
}
