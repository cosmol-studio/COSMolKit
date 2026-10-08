// Test-only projection of unchanged UTF-8 spelling fixtures. Raw byte
// boundaries are asserted without conversion in smarts_counted_bytes.rs.
fn fixture_written_text(value: cosmolkit_model::PropertyText) -> String {
    std::str::from_utf8(value.as_bytes())
        .expect("unchanged UTF-8 writer fixture bytes")
        .to_owned()
}

#[test]
fn directional_smarts_writer_clears_unpaired_markers_and_preserves_complete_stereo() {
    // RDKit 2026.03.1: MolToSmarts(FragmentSmartsConstruct -> Canon).
    // Canon.cpp clears input directions, then rebuilds only valid stereo.
    for (input, expected) in [
        ("C/C", "CC"),
        (r"C\C", "CC"),
        ("c/c", "cc"),
        (r"c\c", "cc"),
        ("N/C", "NC"),
        ("O/C", "OC"),
        ("C/C=C", "CC=C"),
        (r"C\C=C", "CC=C"),
        ("C/C=C/C", "C/C=C/C"),
        (r"C/C=C\C", r"C/C=C\C"),
    ] {
        let query = cosmolkit_search::parse_smarts(input, &Default::default()).unwrap();
        let original = query.clone();
        let written = query_graph_to_smarts(&query, &Default::default())
            .unwrap_or_else(|error| panic!("{input}: {error}"));
        assert_eq!(written.as_bytes(), expected.as_bytes(), "{input}");
        assert_eq!(query, original, "writer changed input {input}");
    }
}
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

#[test]
fn stereo_query_writes_accept_the_source_empty_ring_cache() {
    // RDKit 2026.03.1 MolToSmarts/MolToCXSmarts. These queries reach
    // Canon's stereo perception with FragmentSmartsConstruct's empty cache.
    let profiles = [
        SmartsWriteParams::default(),
        SmartsWriteParams {
            do_isomeric_smiles: false,
            ..Default::default()
        },
        SmartsWriteParams {
            include_atom_maps: false,
            ..Default::default()
        },
        SmartsWriteParams {
            include_dative_bonds: false,
            ..Default::default()
        },
        SmartsWriteParams {
            rooted_at_atom: Some(0),
            ..Default::default()
        },
    ];
    for (input, expected) in [
        ("[#7;!@SP1]=[C,N]", ["[!*]=[C,N]"; 5]),
        ("[*;!@TH1]=[C,N]", ["[!*]=[C,N]"; 5]),
        (
            "[*;@@:9]=[C,N]",
            [
                "[C,N]=[*@@:9]",
                "[C,N]=[*:9]",
                "[C,N]=[*@@]",
                "[C,N]=[*@@:9]",
                "[*@@:9]=[C,N]",
            ],
        ),
        ("[A&@SP1]=[C,N]", ["A=[C,N]"; 5]),
        ("[A;@TH1]=[C,N]", ["A=[C,N]"; 5]),
        (
            "[H;@@:9]=[C,N]",
            [
                "[C,N]=[H1@@:9]",
                "[C,N]=[H1:9]",
                "[C,N]=[H1@@]",
                "[C,N]=[H1@@:9]",
                "[H1@@:9]=[C,N]",
            ],
        ),
        ("[O;!@TH2]=[C,N]", ["[!*]=[C,N]"; 5]),
        (
            "[O;@TH2:9]=[C,N]",
            [
                "[O:9]=[C,N]",
                "[O:9]=[C,N]",
                "O=[C,N]",
                "[O:9]=[C,N]",
                "[O:9]=[C,N]",
            ],
        ),
        ("[c&@OH1]=[C,N]", ["c=[C,N]"; 5]),
    ] {
        let query = cosmolkit_search::parse_smarts(input, &Default::default()).unwrap();
        let before = query.clone();
        for (profile, expected) in profiles.iter().zip(expected) {
            assert_eq!(
                query_graph_to_smarts(&query, profile).unwrap().as_bytes(),
                expected.as_bytes(),
                "{input}, {profile:?}"
            );
            assert_eq!(query, before);
        }
        for (profile, expected) in profiles[..2].iter().zip(expected) {
            assert_eq!(
                query_graph_to_cx_smarts(&query, profile)
                    .unwrap()
                    .as_bytes(),
                expected.as_bytes(),
                "CX {input}, {profile:?}"
            );
            assert_eq!(query, before);
        }
    }
}

#[test]
fn typed_query_properties_keep_source_order_and_string_projection() {
    let mut atom = QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C));
    atom.set_prop("z", 7_i32).unwrap();
    atom.set_prop("a", -0.0_f64).unwrap();
    atom.set_prop("flag", true).unwrap();
    atom.set_prop("z", 9_i32).unwrap();
    atom.set_prop("atomLabel", 12_i32).unwrap();
    let query = QueryGraph::from_parts(vec![atom], vec![], BTreeMap::new(), vec![], vec![], vec![])
        .unwrap();
    assert_eq!(
        query_graph_to_cx_smarts(&query, &Default::default())
            .map(fixture_written_text)
            .unwrap(),
        "[#6] |$12$,atomProp:0.z.9:0.a.-0:0.flag.1|"
    );
}

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
    query_atom_to_smarts(&atom, params).map(fixture_written_text)
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
    .map(fixture_written_text)
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
fn q95_atom_recursion_dispatch_reuses_children_and_source_defined_wildcard() {
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
    // SmartsWrite.cpp::getAtomSmarts falls through to res << "*" for
    // HasProp; RDKit 2026.03.1 emits exactly [#6&*], with a warning.
    let property_query = AtomQueryPredicate::HasProperty("probe".to_owned());
    assert_eq!(
        write_atom_node(QueryNode::and(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
            QueryNode::predicate(property_query),
        ]))
        .unwrap(),
        "[#6&*]"
    );
    let invalid_type = AtomQueryPredicate::AtomType {
        atomic_number: u8::MAX,
        aromatic: false,
    };
    assert_eq!(
        write_atom_node(QueryNode::predicate(invalid_type.clone())).unwrap_err(),
        // SmartsWrite.cpp calls getElementSymbol; PeriodicTable.h:54
        // throws "Atomic number not found" for this exact invalid number.
        SmartsWriteError::AtomTypeAtomicNumber { atomic_number: 255 }
    );

    let mut symbol_atom = QueryAtom::from_parts(
        Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
    );
    symbol_atom
        .set_prop("smilesSymbol", "C")
        .expect("fixed source writer property is valid");
    assert_eq!(
        query_atom_to_smarts(&symbol_atom, &SmartsWriteParams::default())
            .map(fixture_written_text)
            .unwrap(),
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
        query_graph_to_smarts(&disconnected, &SmartsWriteParams::default())
            .map(fixture_written_text)
            .unwrap(),
        "[#8].[#6@@]"
    );
    let rooted = SmartsWriteParams {
        rooted_at_atom: Some(0),
        ..SmartsWriteParams::default()
    };
    assert_eq!(
        query_graph_to_smarts(&disconnected, &rooted)
            .map(fixture_written_text)
            .unwrap(),
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
        query_graph_to_smarts(&ring_and_component, &SmartsWriteParams::default())
            .map(fixture_written_text)
            .unwrap(),
        "[#6]1-[#8]-[#7]-1.[#9]"
    );
    let rooted_ring = SmartsWriteParams {
        rooted_at_atom: Some(2),
        ..SmartsWriteParams::default()
    };
    assert_eq!(
        query_graph_to_smarts(&ring_and_component, &rooted_ring)
            .map(fixture_written_text)
            .unwrap(),
        "[#7]1-[#6]-[#8]-1.[#9]"
    );
}

#[test]
fn q99_ring_numbering_preserves_source_closure_order_reuse_and_label_boundary() {
    let crossing_at_atom = carbon_query_graph(5, &[(0, 1), (1, 2), (2, 3), (3, 4), (2, 0), (4, 2)]);
    assert_eq!(
        query_graph_to_smarts(&crossing_at_atom, &SmartsWriteParams::default())
            .map(fixture_written_text)
            .unwrap(),
        "[#6]1-[#6]-[#6]-12-[#6]-[#6]-2"
    );

    let separate_triangles =
        carbon_query_graph(6, &[(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3)]);
    assert_eq!(
        query_graph_to_smarts(&separate_triangles, &SmartsWriteParams::default())
            .map(fixture_written_text)
            .unwrap(),
        "[#6]1-[#6]-[#6]-1.[#6]1-[#6]-[#6]-1"
    );

    let mut boundary_edges = (0..11).map(|index| (index, index + 1)).collect::<Vec<_>>();
    boundary_edges.extend((0..10).map(|index| (11, index)));
    let ten_open_rings = carbon_query_graph(12, &boundary_edges);
    assert_eq!(
        query_graph_to_smarts(&ten_open_rings, &SmartsWriteParams::default())
            .map(fixture_written_text)
            .unwrap(),
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
        query_graph_to_smarts(&branched, &SmartsWriteParams::default())
            .map(fixture_written_text)
            .unwrap(),
        "[#6:10](-[#8:11])(-[#7:12])-[#9:13].[#17:14]"
    );
    let rooted = SmartsWriteParams {
        rooted_at_atom: Some(3),
        ..SmartsWriteParams::default()
    };
    assert_eq!(
        query_graph_to_smarts(&branched, &rooted)
            .map(fixture_written_text)
            .unwrap(),
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
        Vec::from([(
            cosmolkit_model::PropertyText::from("source"),
            cosmolkit_model::PropertyValue::from("retained"),
        )]),
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
        .map(fixture_written_text)
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
        .map(fixture_written_text)
        .unwrap(),
        "[#6&-].[#8]"
    );

    assert_eq!(
        query_graph_fragment_to_smarts(&query, &rooted, &[], None)
            .map(fixture_written_text)
            .unwrap_err(),
        SmartsWriteError::EmptyAtomSelection
    );
    assert_eq!(
        query_graph_fragment_to_smarts(&query, &rooted, &[AtomId::new(0)], Some(&[]),)
            .map(fixture_written_text)
            .unwrap_err(),
        SmartsWriteError::EmptyBondSelection
    );
    assert_eq!(
        query_graph_fragment_to_smarts(&query, &rooted, &[AtomId::new(3)], None,)
            .map(fixture_written_text)
            .unwrap_err(),
        SmartsWriteError::FragmentAtomOutOfRange { atom: 3 }
    );
    assert_eq!(
        query_graph_fragment_to_smarts(
            &query,
            &rooted,
            &[AtomId::new(0)],
            Some(&[BondId::new(2)]),
        ).map(fixture_written_text)
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
        query_graph_to_cx_smarts(&query, &rooted)
            .map(fixture_written_text)
            .unwrap(),
        "[#7]~[#8]-[#6] |(3,4,5;1.25,2.5,;0,0,),$;;left$,$_AV:;mid;$,^2:0,atomProp:2.a&#46;b.v&#46;x,C:2.1,Z:0,a:0,2,o1:1|"
    );
    assert_eq!(
        query_graph_fragment_to_cx_smarts(
            &query,
            &rooted,
            &[AtomId::new(0), AtomId::new(1)],
            Some(&[BondId::new(0)]),
        )
        .map(fixture_written_text)
        .unwrap(),
        "[#6]->[#8] |(0,0,;1.25,2.5,),$left;$,$_AV:;mid$,atomProp:0.a&#46;b.v&#46;x,C:0.0,a:0,0,o1:1|"
    );
    assert_eq!(query, before);
}

#[test]
fn source_cx_coordinates_follow_actual_mixed_conformer_order() {
    use cosmolkit_model::{Conformer2D, CoordinateDimension};
    let mut query = QueryGraph::from_parts(
        vec![QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C))],
        vec![],
        BTreeMap::new(),
        vec![Conformer2D::new(9, vec![[1.0, 2.0]])],
        vec![Conformer3D::new(4, vec![[3.0, 4.0, 9.0]], true)],
        vec![],
    )
    .unwrap();
    query
        .set_source_conformer_order(Some(vec![
            CoordinateDimension::TwoD,
            CoordinateDimension::ThreeD,
        ]))
        .unwrap();
    let before = query.clone();
    assert_eq!(
        query_graph_to_cx_smarts(&query, &Default::default())
            .map(fixture_written_text)
            .unwrap(),
        "[#6] |(1,2,)|"
    );
    assert_eq!(query, before);
    query
        .set_source_conformer_order(Some(vec![
            CoordinateDimension::ThreeD,
            CoordinateDimension::TwoD,
        ]))
        .unwrap();
    assert_eq!(
        query_graph_to_cx_smarts(&query, &Default::default())
            .map(fixture_written_text)
            .unwrap(),
        "[#6] |(3,4,9)|"
    );
}

#[test]
fn source_cx_false_is3d_omits_z_without_discarding_xyz_storage() {
    use cosmolkit_model::{Conformer2D, CoordinateDimension};
    let mut query = QueryGraph::from_parts(
        vec![QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C))],
        vec![],
        BTreeMap::new(),
        vec![Conformer2D::new(9, vec![[1.0, 2.0]])],
        vec![Conformer3D::new(4, vec![[3.0, 4.0, 9.0]], false)],
        vec![],
    )
    .unwrap();
    query
        .set_source_conformer_order(Some(vec![
            CoordinateDimension::ThreeD,
            CoordinateDimension::TwoD,
        ]))
        .unwrap();
    let before = query.clone();
    assert_eq!(
        query_graph_to_cx_smarts(&query, &Default::default())
            .map(fixture_written_text)
            .unwrap(),
        "[#6] |(3,4,)|"
    );
    assert_eq!(query, before);
    assert_eq!(query.conformers_3d()[0].coordinates(), &[[3.0, 4.0, 9.0]]);
}

#[test]
fn source_cx_missing_mixed_conformer_order_is_a_typed_failure() {
    use cosmolkit_model::{Conformer2D, CoordinateValidationError};
    let query = QueryGraph::from_parts(
        vec![QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C))],
        vec![],
        BTreeMap::new(),
        vec![Conformer2D::new(9, vec![[1.0, 2.0]])],
        vec![Conformer3D::new(4, vec![[3.0, 4.0, 9.0]], true)],
        vec![],
    )
    .unwrap();
    let error = query_graph_to_cx_smarts(&query, &Default::default()).unwrap_err();
    assert!(matches!(
        error,
        SmartsWriteError::CxCoordinateSource(
            CoordinateValidationError::MissingSourceConformerOrder
        )
    ));
    assert!(
        std::error::Error::source(&error)
            .unwrap()
            .downcast_ref::<CoordinateValidationError>()
            .is_some()
    );
}

#[test]
fn concrete_smarts_valence_failure_keeps_native_cause_and_input() {
    use cosmolkit_model::{CoordinateBlock, MoleculeProperties, ValenceError};
    use std::error::Error;
    let record = cosmolkit_smiles::parse_smiles("CC", &Default::default()).unwrap();
    let mut topology = record.topology;
    topology.bonds[0].set_order(BondOrder::Other);
    let before = topology.clone();
    let error = cosmolkit_search::topology_to_smarts(
        &topology,
        &CoordinateBlock::default(),
        &MoleculeProperties::default(),
        &Default::default(),
        false,
    )
    .unwrap_err();
    // SmartsWrite.cpp::MolToSmarts calls updatePropertyCache(false) before
    // Canon traversal: retain the concrete cause at that earlier boundary.
    let SmartsWriteError::Valence(ref valence) = error else {
        panic!("wrong error: {error:?}");
    };
    assert!(
        matches!(valence, ValenceError::BadBondType { bond: Some(bond), order: BondOrder::Other } if *bond == BondId::new(0))
    );
    assert!(matches!(
        error.source().unwrap().downcast_ref::<ValenceError>(),
        Some(ValenceError::BadBondType { .. })
    ));
    assert_eq!(topology, before);
}
