use std::collections::BTreeMap;

use cosmolkit_model::{
    Atom, AtomId, AtomQueryPredicate, AtomRangeBounds, AtomRangeDataFunction, AtomRangeQuery,
    AtomSpec, Bond, BondId, BondQueryPredicate, BondSpec, Conformer3D, QueryAtom, QueryBond,
    QueryGraph, QueryNode, RecursiveStructureQuery, StereoGroup, StereoGroupKind,
};
use cosmolkit_search::{
    SmartsWriteError, SmartsWriteParams, query_atom_to_smarts, query_bond_to_smarts,
    query_graph_fragment_to_cx_smarts, query_graph_fragment_to_smarts, query_graph_to_cx_smarts,
    query_graph_to_smarts,
};
use cosmolkit_types::{BondDirection, BondOrder, ChiralTag, Element};

fn write_atom_node(predicate: QueryNode<AtomQueryPredicate>) -> Result<String, SmartsWriteError> {
    write_atom_node_with_params(predicate, &SmartsWriteParams::default())
}

fn write_atom_node_with_params(
    predicate: QueryNode<AtomQueryPredicate>,
    params: &SmartsWriteParams,
) -> Result<String, SmartsWriteError> {
    let atom = QueryAtom::from_parts(
        Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
        predicate,
    );
    query_atom_to_smarts(&atom, params)
}

fn write_atom(predicate: AtomQueryPredicate) -> String {
    write_atom_node(QueryNode::predicate(predicate)).expect("fixed atom predicate is writable")
}

fn write_bond_node(
    predicate: QueryNode<BondQueryPredicate>,
    direction: BondDirection,
    atom_to_left_idx: Option<usize>,
    params: &SmartsWriteParams,
) -> Result<String, SmartsWriteError> {
    let bond = Bond::from_spec(
        BondId::new(0),
        BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single).with_direction(direction),
    );
    query_bond_to_smarts(
        &QueryBond::from_parts(bond, predicate),
        params,
        atom_to_left_idx,
    )
}

fn write_bond(predicate: BondQueryPredicate) -> String {
    write_bond_node(
        QueryNode::predicate(predicate),
        BondDirection::None,
        Some(0),
        &SmartsWriteParams::default(),
    )
    .expect("fixed bond predicate is writable")
}

fn carbon_query_graph(atom_count: usize, endpoints: &[(usize, usize)]) -> QueryGraph {
    let atoms = (0..atom_count)
        .map(|index| QueryAtom::new(AtomId::new(index), AtomSpec::new(Element::C)))
        .collect::<Vec<_>>();
    let bonds = endpoints
        .iter()
        .copied()
        .enumerate()
        .map(|(index, (begin, end))| {
            QueryBond::from_parts(
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
                ),
                QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
            )
        })
        .collect::<Vec<_>>();
    QueryGraph::from_parts(
        atoms,
        bonds,
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("fixed carbon query graph is valid")
}

#[test]
fn q92_simple_atom_predicates_preserve_source_spelling_and_ranges() {
    assert_eq!(write_atom(AtomQueryPredicate::AtomicNumber(6)), "[#6]");
    assert_eq!(write_atom(AtomQueryPredicate::Isotope(13)), "[13*]");

    assert_eq!(write_atom(AtomQueryPredicate::FormalCharge(-1)), "[-]");
    assert_eq!(write_atom(AtomQueryPredicate::FormalCharge(0)), "[+0]");
    assert_eq!(write_atom(AtomQueryPredicate::FormalCharge(2)), "[+2]");
    assert_eq!(
        write_atom(AtomQueryPredicate::NegativeFormalCharge(-1)),
        "[+]"
    );

    assert_eq!(write_atom(AtomQueryPredicate::HydrogenCount(0)), "[H0]");
    assert_eq!(
        write_atom(AtomQueryPredicate::ImplicitHydrogenCount(2)),
        "[h2]"
    );
    assert_eq!(write_atom(AtomQueryPredicate::HasImplicitHydrogen), "[h]");
    assert_eq!(write_atom(AtomQueryPredicate::IsAromatic(true)), "a");
    assert_eq!(write_atom(AtomQueryPredicate::IsAromatic(false)), "A");
    assert_eq!(
        write_atom(AtomQueryPredicate::AtomType {
            atomic_number: 6,
            aromatic: true,
        }),
        "c"
    );

    let degree_range = |bounds| {
        write_atom(AtomQueryPredicate::Range(AtomRangeQuery::new(
            bounds,
            AtomRangeDataFunction::ExplicitDegree,
        )))
    };
    assert_eq!(degree_range(AtomRangeBounds::LessEqual(2)), "[D{2-}]");
    assert_eq!(degree_range(AtomRangeBounds::GreaterEqual(2)), "[D{-2}]");
    assert_eq!(
        degree_range(AtomRangeBounds::Inclusive {
            lower: 2,
            upper: 4,
            lower_open: true,
            upper_open: true,
        }),
        "[D{2-4}]"
    );
}

#[test]
fn q93_atom_boolean_serialization_preserves_precedence_and_negation() {
    let atomic_number = |value| QueryNode::predicate(AtomQueryPredicate::AtomicNumber(value));
    let hydrogen = |value| QueryNode::predicate(AtomQueryPredicate::HydrogenCount(value));

    assert_eq!(
        write_atom_node(QueryNode::and(vec![atomic_number(6), hydrogen(1)])).unwrap(),
        "[#6&H1]"
    );
    assert_eq!(
        write_atom_node(QueryNode::or(vec![atomic_number(6), atomic_number(7)])).unwrap(),
        "[#6,#7]"
    );
    assert_eq!(
        write_atom_node(QueryNode::and(vec![
            QueryNode::or(vec![atomic_number(6), atomic_number(7)]),
            hydrogen(1),
        ]))
        .unwrap(),
        "[#6,#7;H1]"
    );
    assert_eq!(
        write_atom_node(QueryNode::or(vec![
            QueryNode::and(vec![atomic_number(6), hydrogen(1)]),
            atomic_number(7),
        ]))
        .unwrap(),
        "[#6&H1,#7]"
    );

    assert_eq!(
        write_atom_node(QueryNode::not(QueryNode::and(vec![
            atomic_number(6),
            hydrogen(1),
        ])))
        .unwrap(),
        "[!#6,!H1]"
    );
    assert_eq!(
        write_atom_node(QueryNode::not(QueryNode::or(vec![
            atomic_number(6),
            hydrogen(1),
        ])))
        .unwrap(),
        "[!#6&!H1]"
    );
    assert_eq!(
        write_atom_node(QueryNode::and(vec![
            QueryNode::not(atomic_number(6)),
            hydrogen(1),
        ]))
        .unwrap(),
        "[!#6&H1]"
    );

    let non_smartable = QueryNode::or(vec![
        QueryNode::and(vec![
            QueryNode::or(vec![atomic_number(6), atomic_number(7)]),
            hydrogen(1),
        ]),
        atomic_number(8),
    ]);
    assert_eq!(
        write_atom_node(non_smartable).unwrap_err(),
        SmartsWriteError::OrAboveAndBelowAnd
    );
}

fn mapped_carbon_oxygen_query() -> QueryGraph {
    let mut carbon = QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C));
    carbon.set_atom_map(Some(7));
    let mut oxygen = QueryAtom::new(AtomId::new(1), AtomSpec::new(Element::O));
    oxygen.set_atom_map(Some(9));
    let bond = Bond::from_spec(
        BondId::new(0),
        BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
    );
    QueryGraph::from_parts(
        vec![carbon, oxygen],
        vec![QueryBond::from_parts(
            bond,
            QueryNode::predicate(cosmolkit_model::BondQueryPredicate::Order(
                BondOrder::Single,
            )),
        )],
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("fixed recursive query graph is valid")
}

#[test]
fn q94_recursive_serialization_uses_owned_graph_root_and_map_options() {
    let nested_graph = mapped_carbon_oxygen_query();
    let recursive =
        RecursiveStructureQuery::from_query_graph(nested_graph, 77).with_source_smarts("$([N])");
    let node = || QueryNode::predicate(AtomQueryPredicate::RecursiveSmarts(recursive.clone()));

    assert_eq!(write_atom_node(node()).unwrap(), "[$([#6:7]-[#8:9])]");
    assert_eq!(
        write_atom_node(QueryNode::not(node())).unwrap(),
        "[!$([#6:7]-[#8:9])]"
    );

    let no_maps = SmartsWriteParams {
        include_atom_maps: false,
        ..SmartsWriteParams::default()
    };
    assert_eq!(
        write_atom_node_with_params(node(), &no_maps).unwrap(),
        "[$([#6]-[#8])]"
    );
    let rooted_at_oxygen = SmartsWriteParams {
        rooted_at_atom: Some(1),
        ..SmartsWriteParams::default()
    };
    assert_eq!(
        write_atom_node_with_params(node(), &rooted_at_oxygen).unwrap(),
        "[$([#8:9]-[#6:7])]"
    );

    let mut changed_clone = recursive.clone();
    changed_clone
        .query_graph_mut()
        .unwrap()
        .atom_mut(0)
        .unwrap()
        .set_atom_map(Some(88));
    assert_eq!(
        write_atom_node(QueryNode::predicate(AtomQueryPredicate::RecursiveSmarts(
            changed_clone,
        )))
        .unwrap(),
        "[$([#6:88]-[#8:9])]"
    );
    assert_eq!(write_atom_node(node()).unwrap(), "[$([#6:7]-[#8:9])]");
    assert_eq!(recursive.serial_number(), 77);

    assert_eq!(
        write_atom_node(QueryNode::predicate(AtomQueryPredicate::RecursiveSmarts(
            RecursiveStructureQuery::new(),
        )))
        .unwrap_err(),
        SmartsWriteError::MissingRecursiveQueryMolecule
    );
}

#[test]
fn q95_atom_recursion_dispatch_reuses_children_and_rejects_unwritten_leaves() {
    let recursive = QueryNode::predicate(AtomQueryPredicate::RecursiveSmarts(
        RecursiveStructureQuery::from_query_graph(mapped_carbon_oxygen_query(), 5),
    ));
    let composite = QueryNode::and(vec![
        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
        recursive,
        QueryNode::or(vec![
            QueryNode::predicate(AtomQueryPredicate::HydrogenCount(1)),
            QueryNode::predicate(AtomQueryPredicate::IsAromatic(true)),
        ]),
    ]);
    assert_eq!(
        write_atom_node(composite).unwrap(),
        "[#6&$([#6:7]-[#8:9]);H1,a]"
    );

    let unsupported = AtomQueryPredicate::UnsupportedFeature("q95-explicit");
    assert_eq!(
        write_atom_node(QueryNode::predicate(unsupported.clone())).unwrap_err(),
        SmartsWriteError::UnsupportedAtomQuery {
            predicate: unsupported,
        }
    );
    let nested_unsupported = AtomQueryPredicate::HasProperty("probe".to_owned());
    assert_eq!(
        write_atom_node(QueryNode::and(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
            QueryNode::predicate(nested_unsupported.clone()),
        ]))
        .unwrap_err(),
        SmartsWriteError::UnsupportedAtomQuery {
            predicate: nested_unsupported,
        }
    );
    let invalid_type = AtomQueryPredicate::AtomType {
        atomic_number: u8::MAX,
        aromatic: false,
    };
    assert_eq!(
        write_atom_node(QueryNode::predicate(invalid_type.clone())).unwrap_err(),
        SmartsWriteError::UnsupportedAtomQuery {
            predicate: invalid_type,
        }
    );

    let mut symbol_atom = QueryAtom::from_parts(
        Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
    );
    symbol_atom
        .set_prop("smilesSymbol", "C")
        .expect("fixed source writer property is valid");
    assert_eq!(
        query_atom_to_smarts(&symbol_atom, &SmartsWriteParams::default()).unwrap(),
        "[C]"
    );
}

#[test]
fn q96_simple_bond_serialization_preserves_order_direction_ring_and_aromatic_spelling() {
    assert_eq!(write_bond(BondQueryPredicate::Any), "~");
    assert_eq!(write_bond(BondQueryPredicate::IsInRing(true)), "@");
    assert_eq!(
        write_bond(BondQueryPredicate::Order(BondOrder::Single)),
        "-"
    );
    assert_eq!(
        write_bond(BondQueryPredicate::Order(BondOrder::Double)),
        "="
    );
    assert_eq!(
        write_bond(BondQueryPredicate::Order(BondOrder::Triple)),
        "#"
    );
    assert_eq!(
        write_bond(BondQueryPredicate::Order(BondOrder::Quadruple)),
        "$"
    );
    assert_eq!(
        write_bond(BondQueryPredicate::Order(BondOrder::Aromatic)),
        ":"
    );
    assert_eq!(write_bond(BondQueryPredicate::Order(BondOrder::Zero)), "~");

    assert_eq!(
        write_bond(BondQueryPredicate::Direction(BondDirection::EndDownRight)),
        "\\"
    );
    assert_eq!(
        write_bond(BondQueryPredicate::Direction(BondDirection::EndUpRight)),
        "/"
    );

    let single = QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single));
    assert_eq!(
        write_bond_node(
            single.clone(),
            BondDirection::EndDownRight,
            Some(0),
            &SmartsWriteParams::default(),
        )
        .unwrap(),
        "\\"
    );
    let non_isomeric = SmartsWriteParams {
        do_isomeric_smiles: false,
        ..SmartsWriteParams::default()
    };
    assert_eq!(
        write_bond_node(single, BondDirection::EndUpRight, Some(0), &non_isomeric,).unwrap(),
        "-"
    );

    let aromatic = QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Aromatic));
    assert_eq!(
        write_bond_node(
            aromatic.clone(),
            BondDirection::EndUpRight,
            Some(0),
            &SmartsWriteParams::default(),
        )
        .unwrap(),
        "/"
    );
    assert_eq!(
        write_bond_node(
            aromatic,
            BondDirection::EndDownRight,
            Some(0),
            &non_isomeric,
        )
        .unwrap(),
        ":"
    );

    let dative = QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Dative));
    assert_eq!(
        write_bond_node(
            dative.clone(),
            BondDirection::None,
            Some(0),
            &SmartsWriteParams::default(),
        )
        .unwrap(),
        "->"
    );
    assert_eq!(
        write_bond_node(
            dative.clone(),
            BondDirection::None,
            Some(1),
            &SmartsWriteParams::default(),
        )
        .unwrap(),
        "<-"
    );
    let no_dative = SmartsWriteParams {
        include_dative_bonds: false,
        ..SmartsWriteParams::default()
    };
    assert_eq!(
        write_bond_node(dative, BondDirection::None, Some(0), &no_dative).unwrap(),
        "-"
    );

    assert_eq!(
        write_bond(BondQueryPredicate::OrderIn(vec![
            BondOrder::Single,
            BondOrder::Double,
        ])),
        "-,="
    );
    assert_eq!(
        write_bond(BondQueryPredicate::OrderIn(vec![
            BondOrder::Double,
            BondOrder::Aromatic,
        ])),
        "=,:"
    );
    assert_eq!(
        write_bond(BondQueryPredicate::OrderIn(vec![
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Aromatic,
        ])),
        "-,=,:"
    );
    assert_eq!(
        write_bond_node(
            QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                BondOrder::Single,
                BondOrder::Aromatic,
            ])),
            BondDirection::EndDownRight,
            Some(0),
            &non_isomeric,
        )
        .unwrap(),
        "\\"
    );
}

#[test]
fn q97_bond_boolean_serialization_preserves_precedence_negation_and_binary_state() {
    let order = |value| QueryNode::predicate(BondQueryPredicate::Order(value));
    let ring = || QueryNode::predicate(BondQueryPredicate::IsInRing(true));
    let write = |node| {
        write_bond_node(
            node,
            BondDirection::None,
            Some(0),
            &SmartsWriteParams::default(),
        )
    };

    assert_eq!(
        write(QueryNode::and(vec![order(BondOrder::Single), ring()])).unwrap(),
        "-&@"
    );
    assert_eq!(
        write(QueryNode::or(vec![
            order(BondOrder::Single),
            order(BondOrder::Double),
        ]))
        .unwrap(),
        "-,="
    );
    assert_eq!(
        write(QueryNode::not(order(BondOrder::Single))).unwrap(),
        "!-"
    );
    assert_eq!(
        write(QueryNode::not(QueryNode::and(vec![
            order(BondOrder::Single),
            ring(),
        ])))
        .unwrap(),
        "!-,!@"
    );
    assert_eq!(
        write(QueryNode::not(QueryNode::or(vec![
            order(BondOrder::Single),
            order(BondOrder::Double),
        ])))
        .unwrap(),
        "!-&!="
    );

    assert_eq!(
        write(QueryNode::and(vec![
            QueryNode::or(vec![order(BondOrder::Single), order(BondOrder::Double)]),
            ring(),
        ]))
        .unwrap(),
        "-,=;@"
    );
    assert_eq!(
        write(QueryNode::or(vec![
            QueryNode::and(vec![order(BondOrder::Single), ring()]),
            order(BondOrder::Double),
        ]))
        .unwrap(),
        "-&@,="
    );

    // The pinned source assigns a nested second child to csmarts1 and leaves
    // csmarts2 empty. Keep that observable output instead of correcting it.
    assert_eq!(
        write(QueryNode::and(vec![
            order(BondOrder::Single),
            QueryNode::or(vec![order(BondOrder::Double), order(BondOrder::Triple)]),
        ]))
        .unwrap(),
        "=,#"
    );

    let non_smartable = QueryNode::or(vec![
        QueryNode::and(vec![
            QueryNode::or(vec![order(BondOrder::Single), order(BondOrder::Double)]),
            ring(),
        ]),
        order(BondOrder::Triple),
    ]);
    assert_eq!(
        write(non_smartable).unwrap_err(),
        SmartsWriteError::OrAboveAndBelowAnd
    );

    for malformed in [
        QueryNode::and(Vec::new()),
        QueryNode::or(vec![order(BondOrder::Single)]),
        QueryNode::and(vec![
            order(BondOrder::Single),
            order(BondOrder::Double),
            order(BondOrder::Triple),
        ]),
    ] {
        assert_eq!(
            write(malformed).unwrap_err(),
            SmartsWriteError::CompositeChildCount { kind: "bond" }
        );
    }
}

#[test]
fn q98_query_traversal_classification_preserves_roots_components_and_ring_edges() {
    let disconnected = QueryGraph::from_parts(
        vec![
            QueryAtom::new(
                AtomId::new(0),
                AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCw),
            ),
            QueryAtom::new(AtomId::new(1), AtomSpec::new(Element::O)),
        ],
        Vec::new(),
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("fixed disconnected query is valid");
    assert_eq!(
        query_graph_to_smarts(&disconnected, &SmartsWriteParams::default()).unwrap(),
        "[#8].[#6@@]"
    );
    let rooted = SmartsWriteParams {
        rooted_at_atom: Some(0),
        ..SmartsWriteParams::default()
    };
    assert_eq!(
        query_graph_to_smarts(&disconnected, &rooted).unwrap(),
        "[#6@@].[#8]"
    );

    let atoms = [Element::C, Element::O, Element::N, Element::F]
        .into_iter()
        .enumerate()
        .map(|(index, element)| QueryAtom::new(AtomId::new(index), AtomSpec::new(element)))
        .collect::<Vec<_>>();
    let bonds = [(0, 1), (1, 2), (2, 0)]
        .into_iter()
        .enumerate()
        .map(|(index, (begin, end))| {
            QueryBond::from_parts(
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
                ),
                QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
            )
        })
        .collect::<Vec<_>>();
    let ring_and_component = QueryGraph::from_parts(
        atoms,
        bonds,
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("fixed ring query is valid");
    assert_eq!(
        query_graph_to_smarts(&ring_and_component, &SmartsWriteParams::default()).unwrap(),
        "[#6]1-[#8]-[#7]-1.[#9]"
    );
    let rooted_ring = SmartsWriteParams {
        rooted_at_atom: Some(2),
        ..SmartsWriteParams::default()
    };
    assert_eq!(
        query_graph_to_smarts(&ring_and_component, &rooted_ring).unwrap(),
        "[#7]1-[#6]-[#8]-1.[#9]"
    );
}

#[test]
fn q99_ring_numbering_preserves_source_closure_order_reuse_and_label_boundary() {
    let crossing_at_atom = carbon_query_graph(5, &[(0, 1), (1, 2), (2, 3), (3, 4), (2, 0), (4, 2)]);
    assert_eq!(
        query_graph_to_smarts(&crossing_at_atom, &SmartsWriteParams::default()).unwrap(),
        "[#6]1-[#6]-[#6]-12-[#6]-[#6]-2"
    );

    let separate_triangles =
        carbon_query_graph(6, &[(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3)]);
    assert_eq!(
        query_graph_to_smarts(&separate_triangles, &SmartsWriteParams::default()).unwrap(),
        "[#6]1-[#6]-[#6]-1.[#6]1-[#6]-[#6]-1"
    );

    let mut boundary_edges = (0..11).map(|index| (index, index + 1)).collect::<Vec<_>>();
    boundary_edges.extend((0..10).map(|index| (11, index)));
    let ten_open_rings = carbon_query_graph(12, &boundary_edges);
    assert_eq!(
        query_graph_to_smarts(&ten_open_rings, &SmartsWriteParams::default()).unwrap(),
        "[#6]1-[#6]2-[#6]3-[#6]4-[#6]5-[#6]6-[#6]7-[#6]8-[#6]9-[#6]%10-[#6]-[#6]-1-2-3-4-5-6-7-8-9-%10"
    );
}

#[test]
fn q100_graph_emission_preserves_branches_components_and_output_identity() {
    let mut atoms = [Element::C, Element::O, Element::N, Element::F, Element::CL]
        .into_iter()
        .enumerate()
        .map(|(index, element)| QueryAtom::new(AtomId::new(index), AtomSpec::new(element)))
        .collect::<Vec<_>>();
    for (index, atom) in atoms.iter_mut().enumerate() {
        atom.set_atom_map(Some(u32::try_from(index + 10).expect("small fixed map")));
    }
    let bonds = [(0, 1), (0, 2), (0, 3)]
        .into_iter()
        .enumerate()
        .map(|(index, (begin, end))| {
            QueryBond::from_parts(
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
                ),
                QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
            )
        })
        .collect::<Vec<_>>();
    let branched = QueryGraph::from_parts(
        atoms,
        bonds,
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("fixed branched query is valid");

    assert_eq!(
        query_graph_to_smarts(&branched, &SmartsWriteParams::default()).unwrap(),
        "[#6:10](-[#8:11])(-[#7:12])-[#9:13].[#17:14]"
    );
    let rooted = SmartsWriteParams {
        rooted_at_atom: Some(3),
        ..SmartsWriteParams::default()
    };
    assert_eq!(
        query_graph_to_smarts(&branched, &rooted).unwrap(),
        "[#9:13]-[#6:10](-[#8:11])-[#7:12].[#17:14]"
    );
}

#[test]
fn q101_fragment_selection_validates_indices_and_retains_query_state() {
    let atoms = vec![
        QueryAtom::from_parts(
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                QueryNode::predicate(AtomQueryPredicate::FormalCharge(-1)),
            ]),
        ),
        QueryAtom::new(AtomId::new(1), AtomSpec::new(Element::O)),
        QueryAtom::new(AtomId::new(2), AtomSpec::new(Element::N)),
    ];
    let bonds = [(0, 1), (1, 2)]
        .into_iter()
        .enumerate()
        .map(|(index, (begin, end))| {
            QueryBond::from_parts(
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
                ),
                QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
            )
        })
        .collect::<Vec<_>>();
    let query = QueryGraph::from_parts(
        atoms,
        bonds,
        BTreeMap::from([("source".to_owned(), "retained".to_owned())]),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("fixed fragment query is valid");
    let before = query.clone();
    let rooted = SmartsWriteParams {
        rooted_at_atom: Some(2),
        ..SmartsWriteParams::default()
    };

    assert_eq!(
        query_graph_fragment_to_smarts(
            &query,
            &rooted,
            &[AtomId::new(2), AtomId::new(0), AtomId::new(1)],
            Some(&[BondId::new(1)]),
        )
        .unwrap(),
        "[#6&-].[#8]-[#7]"
    );
    assert_eq!(
        query_graph_fragment_to_smarts(
            &query,
            &SmartsWriteParams::default(),
            &[AtomId::new(0), AtomId::new(1)],
            Some(&[BondId::new(1)]),
        )
        .unwrap(),
        "[#6&-].[#8]"
    );

    assert_eq!(
        query_graph_fragment_to_smarts(&query, &rooted, &[], None).unwrap_err(),
        SmartsWriteError::EmptyAtomSelection
    );
    assert_eq!(
        query_graph_fragment_to_smarts(&query, &rooted, &[AtomId::new(0)], Some(&[]),).unwrap_err(),
        SmartsWriteError::EmptyBondSelection
    );
    assert_eq!(
        query_graph_fragment_to_smarts(&query, &rooted, &[AtomId::new(3)], None,).unwrap_err(),
        SmartsWriteError::FragmentAtomOutOfRange { atom: 3 }
    );
    assert_eq!(
        query_graph_fragment_to_smarts(
            &query,
            &rooted,
            &[AtomId::new(0)],
            Some(&[BondId::new(2)]),
        )
        .unwrap_err(),
        SmartsWriteError::FragmentBondOutOfRange { bond: 2 }
    );
    assert_eq!(query, before);
}

#[test]
fn q102_cx_query_output_maps_properties_coordinates_stereo_and_supported_bonds() {
    let mut carbon = QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C));
    carbon
        .set_prop("atomLabel", "left")
        .expect("fixed label property is valid");
    carbon
        .set_prop("a.b", "v.x")
        .expect("fixed ordinary property is valid");
    let mut oxygen = QueryAtom::new(AtomId::new(1), AtomSpec::new(Element::O));
    oxygen
        .set_prop("molFileValue", "mid")
        .expect("fixed molfile value is valid");
    let mut nitrogen = QueryAtom::new(AtomId::new(2), AtomSpec::new(Element::N));
    nitrogen.set_radical_electrons(2);

    let bonds = [BondOrder::Dative, BondOrder::Zero]
        .into_iter()
        .enumerate()
        .map(|(index, order)| {
            QueryBond::from_parts(
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(index), AtomId::new(index + 1), order),
                ),
                QueryNode::predicate(BondQueryPredicate::Order(order)),
            )
        })
        .collect::<Vec<_>>();
    let query = QueryGraph::from_parts(
        vec![carbon, oxygen, nitrogen],
        bonds,
        BTreeMap::new(),
        Vec::new(),
        vec![Conformer3D::new(
            0,
            vec![[0.0, 0.00001, 0.0], [1.25, 2.5, 0.0], [3.0, 4.0, 5.0]],
            true,
        )],
        vec![
            StereoGroup::new(StereoGroupKind::Or, vec![AtomId::new(1)], Vec::new()).with_id(91),
            StereoGroup::new(
                StereoGroupKind::Absolute,
                vec![AtomId::new(0), AtomId::new(2)],
                Vec::new(),
            )
            .with_id(77),
        ],
    )
    .expect("fixed CX query graph is valid");
    let before = query.clone();
    let rooted = SmartsWriteParams {
        rooted_at_atom: Some(2),
        ..SmartsWriteParams::default()
    };

    assert_eq!(
        query_graph_to_cx_smarts(&query, &rooted).unwrap(),
        "[#7]~[#8]-[#6] |(3,4,5;1.25,2.5,;0,0,),$;;left$,$_AV:;mid;$,^2:0,atomProp:2.a&#46;b.v&#46;x,C:2.1,Z:0,a:0,2,o1:1|"
    );
    assert_eq!(
        query_graph_fragment_to_cx_smarts(
            &query,
            &rooted,
            &[AtomId::new(0), AtomId::new(1)],
            Some(&[BondId::new(0)]),
        )
        .unwrap(),
        "[#6]->[#8] |(0,0,;1.25,2.5,),$left;$,$_AV:;mid$,atomProp:0.a&#46;b.v&#46;x,C:0.0,a:0,0,o1:1|"
    );
    assert_eq!(query, before);
}
