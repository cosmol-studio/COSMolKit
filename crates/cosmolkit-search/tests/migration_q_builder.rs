use cosmolkit_model::{Atom, AtomId, AtomSpec, QueryAtom, QueryAtomIdentity, TopologyBlock};
use cosmolkit_search::{
    AtomQueryPredicate, BondQueryPredicate, QueryGraph, QueryNode, SmartsParseParams,
    SmartsWriteParams, match_query, parse_smarts, write_smarts,
};
use cosmolkit_types::{BondOrder, Element};

// Expected outcomes in this file are read from pinned RDKit 2026.03.1
// (351f8f378f8ad6bbd517980c38896e66bf907af8), principally
// SmilesParse/SmartsParseOps.cpp::ClearAtomChemicalProps and the matching
// reductions/carrier constructors in SmilesParse/smarts.yy. They do not come
// from the historical Rust implementation or generated expected-data caches.

fn query_contains<T, F>(node: &QueryNode<T>, predicate: F) -> bool
where
    F: Fn(&T) -> bool + Copy,
{
    match node {
        QueryNode::Predicate(value) => predicate(value),
        QueryNode::And(children) | QueryNode::Or(children) | QueryNode::Xor(children) => children
            .iter()
            .any(|child| query_contains(child, predicate)),
        QueryNode::Not(child) => query_contains(child, predicate),
    }
}

fn recursive_query(graph: &QueryGraph) -> &cosmolkit_model::RecursiveStructureQuery {
    match graph.atom(0).expect("one recursive-query atom").predicate() {
        QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(query)) => query,
        other => panic!("expected recursive SMARTS predicate, got {other:?}"),
    }
}

#[test]
fn smarts_builder_returns_canonical_graph_with_stable_row_identity() {
    let graph: QueryGraph = parse_smarts(
        "[C:7]-[N:8]=[O:9]",
        &SmartsParseParams {
            skip_cleanup: true,
            ..SmartsParseParams::default()
        },
    )
    .expect("source-valid SMARTS graph should finish");

    assert_eq!(graph.num_atoms(), 3);
    assert_eq!(graph.num_bonds(), 2);
    assert_eq!(
        graph
            .atoms()
            .iter()
            .map(|atom| atom.id().index())
            .collect::<Vec<_>>(),
        [0, 1, 2]
    );
    assert_eq!(graph.atom(0).unwrap().atom_map(), Some(7));
    assert_eq!(graph.atom(1).unwrap().atom_map(), Some(8));
    assert_eq!(graph.atom(2).unwrap().atom_map(), Some(9));
    assert_eq!(
        graph
            .bonds()
            .iter()
            .map(|bond| (bond.id().index(), bond.endpoints()))
            .collect::<Vec<_>>(),
        [(0, (0, 1)), (1, (1, 2))]
    );
    assert_eq!(
        graph.bond(0).unwrap().bond().prop("_cxsmilesBondIdx"),
        Some("0")
    );
    assert_eq!(
        graph.bond(1).unwrap().bond().prop("_cxsmilesBondIdx"),
        Some("1")
    );
    assert_eq!(
        graph.adjacency(),
        &[vec![(1, 0)], vec![(0, 0), (2, 1)], vec![(1, 1)]]
    );
    graph
        .validate()
        .expect("finish must produce a valid canonical graph");

    let empty = parse_smarts("", &SmartsParseParams::default())
        .expect("empty source graph must still finish as one QueryGraph");
    assert_eq!((empty.num_atoms(), empty.num_bonds()), (0, 0));
    empty
        .validate()
        .expect("empty finished graph is structurally valid");
}

#[test]
fn smarts_and_reduction_clears_carrier_isotope_without_erasing_query_leaf() {
    // `atom_expr: atom_expr SEMI_TOKEN atom_expr` calls ClearAtomChemicalProps
    // after expanding the query, so isotope 13 remains a query leaf but the
    // QueryAtom carrier isotope is reset to source value zero.
    let graph = parse_smarts("[13C;N]", &SmartsParseParams::default())
        .expect("valid pinned-source isotope conjunction");
    let atom = graph.atom(0).expect("one query atom");
    assert_eq!(atom.isotope(), None);
    assert!(query_contains(atom.predicate(), |value| {
        *value == AtomQueryPredicate::Isotope(13)
    }));
}

#[test]
fn atom_carrier_fields_follow_each_source_reduction_branch() {
    // `atom_expr_and_point_query` keeps the left QueryAtom's chemical
    // properties; AND/OR/semicolon reductions call ClearAtomChemicalProps,
    // while OR and NOT additionally set only the carrier atomic number to 0.
    let cases = [
        ("[13C&N]", 6, Some(13)),
        ("[13CN]", 6, Some(13)),
        ("[13C,N]", 0, None),
        ("[!13C]", 0, None),
    ];
    for (smarts, atomic_number, isotope) in cases {
        let graph = parse_smarts(smarts, &SmartsParseParams::default())
            .unwrap_or_else(|error| panic!("pinned-source SMARTS {smarts:?}: {error}"));
        let atom = graph.atom(0).expect("one query atom");
        assert_eq!(atom.atomic_number(), atomic_number, "{smarts}");
        assert_eq!(atom.isotope(), isotope, "{smarts}");
        assert!(
            query_contains(atom.predicate(), |value| {
                *value == AtomQueryPredicate::Isotope(13)
            }),
            "the source isotope query leaf must survive {smarts}"
        );
    }
}

#[test]
fn clear_chemical_properties_preserves_no_implicit_and_query_predicates() {
    // ClearAtomChemicalProps resets isotope, formal charge, and explicit H,
    // but not the source noImplicit flag or the independent predicate tree.
    let explicit_h = parse_smarts("[CH3;N]", &SmartsParseParams::default()).unwrap();
    let atom = explicit_h.atom(0).expect("one query atom");
    assert_eq!(atom.atomic_number(), 6);
    assert_eq!(atom.explicit_hydrogens(), 0);
    assert!(atom.no_implicit());
    assert!(query_contains(atom.predicate(), |value| {
        *value == AtomQueryPredicate::HydrogenCount(3)
    }));

    let charge = parse_smarts("[N+;C]", &SmartsParseParams::default()).unwrap();
    let atom = charge.atom(0).expect("one query atom");
    assert_eq!(atom.formal_charge(), 0);
    assert!(query_contains(atom.predicate(), |value| {
        *value == AtomQueryPredicate::FormalCharge(1)
    }));
}

#[test]
fn conflicting_hydrogen_and_charge_queries_clear_only_carrier_values() {
    // The point-query helper retains the source masks and query leaves while
    // clearing a conflicting carrier count/charge according to smarts.yy.
    let hydrogen = parse_smarts("[CH3H2]", &SmartsParseParams::default()).unwrap();
    let atom = hydrogen.atom(0).expect("one query atom");
    assert_eq!(atom.explicit_hydrogens(), 0);
    assert!(atom.no_implicit());
    for count in [2, 3] {
        assert!(query_contains(atom.predicate(), |value| {
            *value == AtomQueryPredicate::HydrogenCount(count)
        }));
    }

    let charge = parse_smarts("[N+-]", &SmartsParseParams::default()).unwrap();
    let atom = charge.atom(0).expect("one query atom");
    assert_eq!(atom.formal_charge(), 0);
    for value in [-1, 1] {
        assert!(query_contains(atom.predicate(), |predicate| {
            *predicate == AtomQueryPredicate::FormalCharge(value)
        }));
    }
}

#[test]
fn smarts_charge_primitive_is_carried_and_remains_a_query_leaf() {
    // `charge_spec` sets formal charge on the QueryAtom and creates a separate
    // formal-charge predicate in pinned SmilesParse/smarts.yy.
    let graph = parse_smarts("[N+]", &SmartsParseParams::default())
        .expect("valid pinned-source charged query");
    let atom = graph.atom(0).expect("one query atom");
    assert_eq!(atom.formal_charge(), 1);
    assert!(query_contains(atom.predicate(), |value| {
        *value == AtomQueryPredicate::FormalCharge(1)
    }));
}

#[test]
fn smarts_explicit_hydrogen_count_and_no_implicit_state_reach_carrier() {
    // `number H_TOKEN number` sets explicit H count 3 and noImplicit=true on
    // the carrier, while retaining the independent H-count query leaf.
    let graph = parse_smarts("[CH3]", &SmartsParseParams::default())
        .expect("valid pinned-source explicit-H query");
    let atom = graph.atom(0).expect("one query atom");
    assert_eq!(atom.explicit_hydrogens(), 3);
    assert!(atom.no_implicit());
    assert!(query_contains(atom.predicate(), |value| {
        *value == AtomQueryPredicate::HydrogenCount(3)
    }));
}

#[test]
fn smarts_bond_boolean_expansion_preserves_source_carrier_bond_type() {
    // QueryBond(UNSPECIFIED) is the carrier for `~` and `@`; boolean
    // expansion retains the left QueryBond carrier. The equality predicate
    // on the right of `@=` must not replace that source carrier order.
    let cases = [
        ("C~C", BondOrder::Unspecified),
        ("C@C", BondOrder::Unspecified),
        ("C@=C", BondOrder::Unspecified),
        ("C=@C", BondOrder::Double),
    ];
    let mut actual_orders = Vec::with_capacity(cases.len());
    for (smarts, _) in cases {
        let graph = parse_smarts(smarts, &SmartsParseParams::default())
            .unwrap_or_else(|error| panic!("pinned-source SMARTS {smarts:?}: {error}"));
        assert_eq!(graph.num_bonds(), 1, "{smarts}");
        actual_orders.push(graph.bond(0).unwrap().bond().order());
    }
    let expected_orders = cases.map(|(_, order)| order);
    assert_eq!(actual_orders, expected_orders);

    let ring_bond = parse_smarts("C@=C", &SmartsParseParams::default()).unwrap();
    assert!(query_contains(
        ring_bond.bond(0).unwrap().predicate(),
        |value| { *value == BondQueryPredicate::IsInRing(true) }
    ));
    assert!(query_contains(
        ring_bond.bond(0).unwrap().predicate(),
        |value| { *value == BondQueryPredicate::Order(BondOrder::Double) }
    ));
}

#[test]
fn implicit_bond_carrier_order_uses_both_source_atom_aromatic_flags() {
    // SmilesParseOps::getUnspecifiedQueryBond initially derives the carrier
    // type from the endpoint aromatic flags; SetUnspecifiedBondTypes repeats
    // that rule after ring pairing for bonds marked with _unspecifiedOrder.
    let cases = [
        ("CC", BondOrder::Single),
        ("cC", BondOrder::Single),
        ("cc", BondOrder::Aromatic),
    ];
    for (smarts, expected_order) in cases {
        let graph = parse_smarts(smarts, &SmartsParseParams::default())
            .unwrap_or_else(|error| panic!("pinned-source SMARTS {smarts:?}: {error}"));
        assert_eq!(graph.num_bonds(), 1, "{smarts}");
        assert_eq!(
            graph.bond(0).unwrap().bond().order(),
            expected_order,
            "{smarts}"
        );
    }
}

#[test]
fn bond_query_expansion_preserves_source_predicate_order_and_origin() {
    // smarts.ll gives `~` a null query and `=` a double-order query;
    // smarts.yy gives `@` an in-ring query. QueryBond::expandQuery retains
    // source child order when it builds a non-null composite.
    let cases: [(&str, QueryNode<BondQueryPredicate>); 4] = [
        ("C~C", QueryNode::predicate(BondQueryPredicate::Any)),
        (
            "C@C",
            QueryNode::predicate(BondQueryPredicate::IsInRing(true)),
        ),
        (
            "C@=C",
            QueryNode::and(vec![
                QueryNode::predicate(BondQueryPredicate::IsInRing(true)),
                QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Double)),
            ]),
        ),
        (
            "C=@C",
            QueryNode::and(vec![
                QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Double)),
                QueryNode::predicate(BondQueryPredicate::IsInRing(true)),
            ]),
        ),
    ];
    for (smarts, expected_predicate) in cases {
        let graph = parse_smarts(smarts, &SmartsParseParams::default())
            .unwrap_or_else(|error| panic!("pinned-source SMARTS {smarts:?}: {error}"));
        let bond = graph.bond(0).expect("one query bond");
        assert_eq!(bond.predicate(), &expected_predicate, "{smarts}");
        assert!(!bond.predicate_is_carrier_derived(), "{smarts}");
    }
}

#[test]
fn parser_rows_keep_explicit_query_origin_separate_from_carrier_derived_rows() {
    let graph = parse_smarts("[N+]-C", &SmartsParseParams::default()).unwrap();
    let parsed_atom = graph.atom(0).expect("query atom");
    assert!(!parsed_atom.predicate_is_carrier_derived());
    let carrier_atom = cosmolkit_model::QueryAtom::from_carrier_parts(
        parsed_atom
            .try_to_atom()
            .expect("parsed Element carrier remains concrete"),
        parsed_atom.predicate().clone(),
    );
    assert!(carrier_atom.predicate_is_carrier_derived());
    assert_ne!(parsed_atom, &carrier_atom);

    let parsed_bond = graph.bond(0).expect("query bond");
    assert!(!parsed_bond.predicate_is_carrier_derived());
    let carrier_bond = cosmolkit_model::QueryBond::from_carrier_parts(
        parsed_bond.bond().clone(),
        parsed_bond.predicate().clone(),
    );
    assert!(carrier_bond.predicate_is_carrier_derived());
    assert_ne!(parsed_bond, &carrier_bond);
}

#[test]
fn recursive_query_keeps_its_owned_graph_default_and_explicit_serials() {
    for (smarts, expected_serial) in [("[$(C=O)]", 100), ("[$(C=O)_17]", 17)] {
        let graph = parse_smarts(smarts, &SmartsParseParams::default())
            .unwrap_or_else(|error| panic!("pinned-source SMARTS {smarts:?}: {error}"));
        assert_eq!((graph.num_atoms(), graph.num_bonds()), (1, 0));
        let atom = graph.atom(0).expect("recursive query atom");
        assert_eq!(atom.atomic_number(), 0);
        assert!(!atom.predicate_is_carrier_derived());

        let recursive = recursive_query(&graph);
        assert_eq!(recursive.source_smarts(), Some("$(C=O)"));
        assert_eq!(recursive.serial_number(), expected_serial);
        let inner = recursive.query_graph().expect("owned recursive graph");
        assert_eq!((inner.num_atoms(), inner.num_bonds()), (2, 1));
        assert_eq!(inner.atom(0).unwrap().atomic_number(), 6);
        assert_eq!(inner.atom(1).unwrap().atomic_number(), 8);
        assert_eq!(inner.bond(0).unwrap().bond().order(), BondOrder::Double);
        assert!(!inner.atom(0).unwrap().predicate_is_carrier_derived());
        assert!(!inner.bond(0).unwrap().predicate_is_carrier_derived());
        inner
            .validate()
            .expect("recursive query graph is canonical");

        let copied = graph.clone();
        let copied_recursive = recursive_query(&copied);
        assert_eq!(copied_recursive.serial_number(), expected_serial);
        assert_eq!(copied_recursive.query_graph(), Some(inner));
        assert_eq!(copied_recursive.source_smarts(), Some("$(C=O)"));
    }
}

#[test]
fn source_query_atom_carrier_retains_numeric_atomic_number_119() {
    // Pinned `HASH_TOKEN number` constructs QueryAtom(int), and Atom(unsigned
    // int) stores the value directly in uint8_t d_atomicNum. 119 is accepted
    // even though it is not a periodic-table Element. The typed query carrier
    // preserves that raw identity independently from the predicate.
    let graph = parse_smarts("[#119]", &SmartsParseParams::default())
        .expect("pinned RDKit accepts numeric query atom 119");
    let atom = graph.atom(0).expect("one query atom");
    assert_eq!(
        atom.identity(),
        cosmolkit_model::QueryAtomIdentity::AtomicNumber(119)
    );
    assert_eq!(atom.atomic_number(), 119);
    assert!(query_contains(atom.predicate(), |value| {
        *value == AtomQueryPredicate::AtomicNumber(119)
    }));

    let expected_conversion_error =
        cosmolkit_model::QueryAtomConversionError::NonElementAtomicNumber {
            atom: cosmolkit_model::AtomId::new(0),
            atomic_number: 119,
        };
    assert_eq!(atom.try_to_atom(), Err(expected_conversion_error));

    let cloned = graph.clone();
    let cloned_atom = cloned.atom(0).expect("cloned raw-identity query atom");
    assert_eq!(cloned_atom.identity(), atom.identity());
    assert_eq!(cloned_atom.predicate(), atom.predicate());
    assert_eq!(cloned_atom.try_to_atom(), Err(expected_conversion_error));
}

fn assert_carrier_fields(
    atom: &cosmolkit_model::QueryAtom,
    input: &str,
    atomic_number: u8,
    isotope: Option<u16>,
    aromatic: bool,
    formal_charge: i8,
    explicit_hydrogens: u8,
    no_implicit: bool,
) {
    assert_eq!(atom.atomic_number(), atomic_number, "{input}");
    assert_eq!(atom.isotope(), isotope, "{input}");
    assert_eq!(atom.is_aromatic(), aromatic, "{input}");
    assert_eq!(atom.formal_charge(), formal_charge, "{input}");
    assert_eq!(atom.explicit_hydrogens(), explicit_hydrogens, "{input}");
    assert_eq!(atom.no_implicit(), no_implicit, "{input}");
}

#[test]
fn each_repeated_not_action_clears_carrier_even_when_negation_cancels() {
    // Every recursive point_query NOT reduction clears the QueryAtom carrier;
    // only the query's outer negation toggles back after an even run.
    let isotopic_carbon = QueryNode::and(vec![
        QueryNode::predicate(AtomQueryPredicate::AtomType {
            atomic_number: 6,
            aromatic: false,
        }),
        QueryNode::predicate(AtomQueryPredicate::Isotope(13)),
    ]);
    for (smarts, negated) in [
        ("[!13C]", true),
        ("[!!13C]", false),
        ("[!!!13C]", true),
        ("[!!!!13C]", false),
    ] {
        let graph = parse_smarts(smarts, &SmartsParseParams::default())
            .unwrap_or_else(|error| panic!("pinned-source SMARTS {smarts:?}: {error}"));
        let atom = graph.atom(0).expect("one query atom");
        assert_carrier_fields(atom, smarts, 0, None, false, 0, 0, false);
        let expected = if negated {
            QueryNode::not(isotopic_carbon.clone())
        } else {
            isotopic_carbon.clone()
        };
        assert_eq!(atom.predicate(), &expected, "{smarts}");
        assert!(!atom.predicate_is_carrier_derived(), "{smarts}");
    }

    // The same source helper clears charge and explicit H on a point query.
    // noImplicit is set by H_TOKEN and is not cleared by ClearAtomChemicalProps.
    let charged = parse_smarts("[!!+]", &SmartsParseParams::default()).unwrap();
    let atom = charged.atom(0).expect("one charged query atom");
    assert_carrier_fields(&atom, "[!!+]", 0, None, false, 0, 0, false);
    assert_eq!(
        atom.predicate(),
        &QueryNode::predicate(AtomQueryPredicate::FormalCharge(1))
    );
    assert!(!atom.predicate_is_carrier_derived());

    let explicit_h = parse_smarts("[!!H3]", &SmartsParseParams::default()).unwrap();
    let atom = explicit_h.atom(0).expect("one hydrogen query atom");
    assert_carrier_fields(&atom, "[!!H3]", 0, None, false, 0, 0, true);
    assert_eq!(
        atom.predicate(),
        &QueryNode::predicate(AtomQueryPredicate::HydrogenCount(3))
    );
    assert!(!atom.predicate_is_carrier_derived());
}

#[test]
fn generic_a_and_aromatic_as_tokens_preserve_carrier_and_predicate_identity() {
    let atom_type = |atomic_number, aromatic| {
        QueryNode::predicate(AtomQueryPredicate::AtomType {
            atomic_number,
            aromatic,
        })
    };
    let aromatic_query = QueryNode::predicate(AtomQueryPredicate::IsAromatic(true));
    let nitrogen = atom_type(7, false);
    let arsenic = atom_type(33, true);
    let cases = [
        ("a", 0, true, aromatic_query.clone()),
        ("[a]", 0, true, aromatic_query.clone()),
        ("[as]", 33, true, arsenic.clone()),
        ("[!as]", 0, true, QueryNode::not(arsenic.clone())),
        (
            "[a&N]",
            0,
            true,
            QueryNode::and(vec![aromatic_query.clone(), nitrogen.clone()]),
        ),
        (
            "[a,N]",
            0,
            true,
            QueryNode::or(vec![aromatic_query.clone(), nitrogen.clone()]),
        ),
        (
            "[as&N]",
            33,
            true,
            QueryNode::and(vec![arsenic.clone(), nitrogen.clone()]),
        ),
        ("[as,N]", 0, true, QueryNode::or(vec![arsenic, nitrogen])),
    ];

    for (smarts, atomic_number, aromatic, expected_predicate) in cases {
        let graph = parse_smarts(smarts, &SmartsParseParams::default())
            .unwrap_or_else(|error| panic!("pinned-source SMARTS {smarts:?}: {error}"));
        let atom = graph.atom(0).expect("one query atom");
        assert_carrier_fields(atom, smarts, atomic_number, None, aromatic, 0, 0, false);
        assert_eq!(atom.predicate(), &expected_predicate, "{smarts}");
        assert!(!atom.predicate_is_carrier_derived(), "{smarts}");
    }
}

#[test]
fn generic_a_aromatic_endpoints_create_aromatic_implicit_bonds() {
    for smarts in ["aa", "[a][a]"] {
        let graph = parse_smarts(smarts, &SmartsParseParams::default())
            .unwrap_or_else(|error| panic!("pinned-source SMARTS {smarts:?}: {error}"));
        assert_eq!((graph.num_atoms(), graph.num_bonds()), (2, 1), "{smarts}");
        for atom in graph.atoms() {
            assert_carrier_fields(atom, smarts, 0, None, true, 0, 0, false);
            assert_eq!(
                atom.predicate(),
                &QueryNode::predicate(AtomQueryPredicate::IsAromatic(true)),
                "{smarts}"
            );
            assert!(!atom.predicate_is_carrier_derived(), "{smarts}");
        }
        assert_eq!(
            graph.bond(0).unwrap().bond().order(),
            BondOrder::Aromatic,
            "{smarts}"
        );
    }
}

#[test]
fn widened_query_predicates_match_and_write_without_narrowing() {
    let parse = |smarts: &str| {
        parse_smarts(smarts, &SmartsParseParams::default())
            .unwrap_or_else(|error| panic!("pinned-source SMARTS {smarts:?}: {error}"))
    };
    let target = |formal_charge, isotope| {
        TopologyBlock::try_from_parts(
            vec![Atom::from_spec(
                AtomId::new(0),
                AtomSpec::new(Element::C)
                    .with_formal_charge(formal_charge)
                    .with_isotope(isotope),
            )],
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .expect("one valid detached carbon topology")
    };
    let query_for = |predicate| {
        QueryGraph::from_parts(
            vec![QueryAtom::from_identity_parts(
                AtomId::new(0),
                QueryAtomIdentity::Element(Element::C),
                predicate,
            )],
            Vec::new(),
            Default::default(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .expect("one valid detached carbon query")
    };

    let formal_charge_128 = parse("[C+128]");
    assert!(
        match_query(&formal_charge_128, &target(-128, 0))
            .expect("match widened formal-charge query")
            .is_empty(),
        "FormalCharge(128) must not compare equal to the narrowed carrier -128"
    );
    assert_eq!(
        match_query(&parse("[C+127]"), &target(127, 0))
            .expect("match ordinary-range formal-charge query")
            .len(),
        1
    );

    let negative_formal_charge_128 = query_for(QueryNode::predicate(
        AtomQueryPredicate::NegativeFormalCharge(128),
    ));
    assert_eq!(
        match_query(&negative_formal_charge_128, &target(-128, 0))
            .expect("match widened negative-formal-charge query")
            .len(),
        1,
        "NegativeFormalCharge(128) compares against -(-128) as i32"
    );
    let negative_formal_charge_1 = query_for(QueryNode::predicate(
        AtomQueryPredicate::NegativeFormalCharge(1),
    ));
    assert_eq!(
        match_query(&negative_formal_charge_1, &target(-1, 0))
            .expect("match ordinary-range negative-formal-charge query")
            .len(),
        1
    );

    let isotope_65_536 = parse("[65536C]");
    assert!(
        match_query(&isotope_65_536, &target(0, 0))
            .expect("match widened isotope query")
            .is_empty(),
        "Isotope(65536) must not match its narrowed carrier value 0"
    );
    assert_eq!(
        match_query(&parse("[65535C]"), &target(0, u16::MAX))
            .expect("match ordinary-range isotope query")
            .len(),
        1
    );

    for (smarts, expected) in [
        ("[C+128]", "[C&+128]"),
        ("[C-129]", "[C&-129]"),
        ("[65536C]", "[C&65536*]"),
    ] {
        assert_eq!(
            write_smarts(&parse(smarts), &SmartsWriteParams::default())
                .expect("write full-width predicate value"),
            expected,
            "{smarts}"
        );
    }
}

#[test]
fn q07d_primitive_targets_match_and_write_without_narrowing() {
    let parse = |smarts: &str| {
        parse_smarts(smarts, &SmartsParseParams::default())
            .unwrap_or_else(|error| panic!("pinned-source SMARTS {smarts:?}: {error}"))
    };
    // With no neighbors and no implicit Hs, pinned RDKit's degree, total
    // degree, valence, H-count and heteroatom-neighbor measurements are zero.
    // Each widened equality leaf must reject 256 and accept 0; narrowing the
    // source int target to the carrier's u8 would incorrectly match 256.
    let zero_target = TopologyBlock::try_from_parts(
        vec![Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C).with_no_implicit(true),
        )],
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("one no-implicit detached carbon with zero query counts");
    let target_before = zero_target.clone();

    for (primitive, source_256) in [
        ("D", "[C&D256]"),
        ("X", "[C&X256]"),
        ("v", "[C&v256]"),
        ("h", "[C&h256]"),
        ("H", "[C&H256]"),
        ("z", "[C&z256]"),
        ("Z", "[C&Z256]"),
    ] {
        let widened = parse(source_256);
        let widened_before = widened.clone();
        assert!(
            match_query(&widened, &zero_target)
                .expect("evaluate full-width source target")
                .is_empty(),
            "{source_256} must not wrap to target 0"
        );

        let zero_source = format!("[C&{primitive}0]");
        let zero_query = parse(&zero_source);
        assert_eq!(
            match_query(&zero_query, &zero_target)
                .expect("evaluate ordinary zero target")
                .len(),
            1,
            "{zero_source} is the source-defined zero-measure control"
        );
        assert_eq!(
            write_smarts(&widened, &SmartsWriteParams::default())
                .expect("write full-width source target"),
            source_256,
            "writer must retain the complete target"
        );
        assert_eq!(
            widened, widened_before,
            "matching/writing mutated {source_256}"
        );
    }

    // Range predicates remain their separately modeled narrow variants and
    // retain the pinned inclusive open-end comparisons at degree zero.
    assert_eq!(
        match_query(&parse("[C&D{0-}]"), &zero_target)
            .expect("evaluate explicit-degree upper range")
            .len(),
        1
    );
    assert!(
        match_query(&parse("[C&D{-1}]"), &zero_target)
            .expect("evaluate explicit-degree lower range")
            .is_empty()
    );
    assert_eq!(
        zero_target, target_before,
        "matching mutated the detached target"
    );
}
