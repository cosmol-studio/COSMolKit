use std::collections::BTreeMap;
use std::sync::{
    Arc,
    atomic::{AtomicUsize, Ordering},
};

use cosmolkit_model::{
    Atom, AtomId, AtomQueryPredicate, AtomSpec, Bond, BondId, BondQueryPredicate, BondSpec,
    Conformer3D, CoordinateBlock, QueryAtom, QueryBond, QueryGraph, QueryNode, TopologyBlock,
};
use cosmolkit_search::{
    AtomCoordsMatchFunctor, SearchTarget, SearchTargetAccess, SubstructMatchError,
    SubstructMatchParams, SubstructMatchParamsJsonError, get_substruct_matches_with_params,
    substruct_match_params_to_json, try_get_substruct_matches_with_params,
    update_substruct_match_params_from_json,
};
use cosmolkit_types::{BondOrder, Element};

fn atom(id: usize, element: Element) -> Atom {
    Atom::from_spec(AtomId::new(id), AtomSpec::new(element))
}

fn single_atom_query(element: Element) -> QueryGraph {
    let atomic_number = element.atomic_number();
    QueryGraph::from_parts(
        vec![QueryAtom::from_parts(
            atom(0, element),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(atomic_number)),
        )],
        Vec::new(),
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("fixed single-atom query is valid")
}

fn single_bond_query(order: BondOrder) -> QueryGraph {
    let atoms = vec![
        QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C)),
        QueryAtom::new(AtomId::new(1), AtomSpec::new(Element::O)),
    ];
    let bond = Bond::from_spec(
        BondId::new(0),
        BondSpec::new(AtomId::new(0), AtomId::new(1), order),
    );
    QueryGraph::from_parts(
        atoms,
        vec![QueryBond::from_parts(
            bond,
            QueryNode::predicate(BondQueryPredicate::Order(order)),
        )],
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("fixed single-bond query is valid")
}

fn single_atom_topology(element: Element) -> TopologyBlock {
    TopologyBlock::try_from_parts(vec![atom(0, element)], Vec::new(), Vec::new(), Vec::new())
        .expect("fixed single-atom target is valid")
}

fn single_bond_topology(order: BondOrder) -> TopologyBlock {
    let atoms = vec![atom(0, Element::C), atom(1, Element::O)];
    let bonds = vec![Bond::from_spec(
        BondId::new(0),
        BondSpec::new(AtomId::new(0), AtomId::new(1), order),
    )];
    TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
        .expect("fixed single-bond target is valid")
}

fn two_site_carbon_oxygen_topology() -> TopologyBlock {
    let atoms = vec![
        atom(0, Element::C),
        atom(1, Element::O),
        atom(2, Element::C),
    ];
    let bonds = vec![
        Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        ),
        Bond::from_spec(
            BondId::new(1),
            BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Single),
        ),
    ];
    TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
        .expect("fixed two-site callback target is valid")
}

fn matches(query: &QueryGraph, topology: &TopologyBlock, params: &SubstructMatchParams) -> bool {
    let coordinates = CoordinateBlock::default();
    let target = SearchTarget::new(topology, &coordinates, &topology.stereo_groups, None, None);
    !get_substruct_matches_with_params(&target, query, params).is_empty()
}

fn coordinate_matches(
    query: &QueryGraph,
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    functor: AtomCoordsMatchFunctor,
) -> bool {
    let target = SearchTarget::new(topology, coordinates, &topology.stereo_groups, None, None);
    let params = SubstructMatchParams {
        extra_atom_check: Some(Arc::new(move |query, query_atom, target, target_atom| {
            functor.matches(query, query_atom, target, target_atom)
        })),
        ..SubstructMatchParams::default()
    };
    !get_substruct_matches_with_params(&target, query, &params).is_empty()
}

#[test]
fn q87_atom_pair_callbacks_preserve_source_order_and_rejection() {
    let oxygen_query = single_atom_query(Element::O);
    let carbon_target = single_atom_topology(Element::C);
    let skipped_calls = Arc::new(AtomicUsize::new(0));
    let skipped_calls_in_check = Arc::clone(&skipped_calls);
    let post_params = SubstructMatchParams {
        extra_atom_check: Some(Arc::new(move |_, _, _, _| {
            skipped_calls_in_check.fetch_add(1, Ordering::SeqCst);
            true
        })),
        ..SubstructMatchParams::default()
    };
    assert!(!matches(&oxygen_query, &carbon_target, &post_params));
    assert_eq!(skipped_calls.load(Ordering::SeqCst), 0);

    let carbon_query = single_atom_query(Element::C);
    let ordered_calls = Arc::new(AtomicUsize::new(0));
    let ordered_calls_in_check = Arc::clone(&ordered_calls);
    let rejecting_params = SubstructMatchParams {
        extra_atom_check: Some(Arc::new(move |query, query_atom, target, target_atom| {
            ordered_calls_in_check.fetch_add(1, Ordering::SeqCst);
            assert_eq!(query_atom.id(), AtomId::new(0));
            assert_eq!(query_atom.atomic_number(), Element::C.atomic_number());
            assert_eq!(query.atoms()[0].id(), query_atom.id());
            assert_eq!(target_atom.id(), AtomId::new(0));
            assert_eq!(target_atom.atomic_number(), Element::C.atomic_number());
            assert_eq!(target.atoms()[0].id(), target_atom.id());
            false
        })),
        ..SubstructMatchParams::default()
    };
    assert!(!matches(&carbon_query, &carbon_target, &rejecting_params));
    assert_eq!(ordered_calls.load(Ordering::SeqCst), 1);

    let override_calls = Arc::new(AtomicUsize::new(0));
    let override_calls_in_check = Arc::clone(&override_calls);
    let override_params = SubstructMatchParams {
        extra_atom_check: Some(Arc::new(move |_, query_atom, _, target_atom| {
            override_calls_in_check.fetch_add(1, Ordering::SeqCst);
            assert_eq!(query_atom.atomic_number(), Element::O.atomic_number());
            assert_eq!(target_atom.atomic_number(), Element::C.atomic_number());
            true
        })),
        extra_atom_check_overrides_default_check: true,
        ..SubstructMatchParams::default()
    };
    assert!(matches(&oxygen_query, &carbon_target, &override_params));
    assert_eq!(override_calls.load(Ordering::SeqCst), 1);
}

#[test]
fn q87_bond_pair_callbacks_preserve_source_order_and_rejection() {
    let double_query = single_bond_query(BondOrder::Double);
    let single_target = single_bond_topology(BondOrder::Single);
    let skipped_calls = Arc::new(AtomicUsize::new(0));
    let skipped_calls_in_check = Arc::clone(&skipped_calls);
    let post_params = SubstructMatchParams {
        extra_bond_check: Some(Arc::new(move |_, _| {
            skipped_calls_in_check.fetch_add(1, Ordering::SeqCst);
            true
        })),
        ..SubstructMatchParams::default()
    };
    assert!(!matches(&double_query, &single_target, &post_params));
    assert_eq!(skipped_calls.load(Ordering::SeqCst), 0);

    let single_query = single_bond_query(BondOrder::Single);
    let ordered_calls = Arc::new(AtomicUsize::new(0));
    let ordered_calls_in_check = Arc::clone(&ordered_calls);
    let rejecting_params = SubstructMatchParams {
        extra_bond_check: Some(Arc::new(move |query_bond, target_bond| {
            ordered_calls_in_check.fetch_add(1, Ordering::SeqCst);
            assert_eq!(query_bond.id(), BondId::new(0));
            assert_eq!(query_bond.order(), BondOrder::Single);
            assert_eq!(target_bond.id(), BondId::new(0));
            assert_eq!(target_bond.order(), BondOrder::Single);
            false
        })),
        ..SubstructMatchParams::default()
    };
    assert!(!matches(&single_query, &single_target, &rejecting_params));
    assert_eq!(ordered_calls.load(Ordering::SeqCst), 1);

    let override_calls = Arc::new(AtomicUsize::new(0));
    let override_calls_in_check = Arc::clone(&override_calls);
    let override_params = SubstructMatchParams {
        extra_bond_check: Some(Arc::new(move |query_bond, target_bond| {
            override_calls_in_check.fetch_add(1, Ordering::SeqCst);
            assert_eq!(query_bond.order(), BondOrder::Double);
            assert_eq!(target_bond.order(), BondOrder::Single);
            true
        })),
        extra_bond_check_overrides_default_check: true,
        ..SubstructMatchParams::default()
    };
    assert!(matches(&double_query, &single_target, &override_params));
    assert_eq!(override_calls.load(Ordering::SeqCst), 1);
}

#[test]
fn q87_unsupported_pair_queries_keep_typed_errors_before_callbacks() {
    let atom_query = QueryGraph::from_parts(
        vec![QueryAtom::from_parts(
            atom(0, Element::C),
            QueryNode::predicate(AtomQueryPredicate::UnsupportedFeature(
                "Q87 unsupported atom callback input",
            )),
        )],
        Vec::new(),
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("fixed unsupported atom query is structurally valid");
    let atom_topology = single_atom_topology(Element::C);
    let coordinates = CoordinateBlock::default();
    let atom_target = SearchTarget::new(
        &atom_topology,
        &coordinates,
        &atom_topology.stereo_groups,
        None,
        None,
    );
    let atom_calls = Arc::new(AtomicUsize::new(0));
    let atom_calls_in_check = Arc::clone(&atom_calls);
    let atom_params = SubstructMatchParams {
        extra_atom_check: Some(Arc::new(move |_, _, _, _| {
            atom_calls_in_check.fetch_add(1, Ordering::SeqCst);
            true
        })),
        extra_atom_check_overrides_default_check: true,
        ..SubstructMatchParams::default()
    };
    assert_eq!(
        try_get_substruct_matches_with_params(&atom_target, &atom_query, &atom_params)
            .expect_err("unsupported atom query must fail before its callback"),
        SubstructMatchError::Unsupported {
            branch: "Q87 unsupported atom callback input",
            rdkit_function: "QueryAtom::Match",
        }
    );
    assert_eq!(atom_calls.load(Ordering::SeqCst), 0);

    let bond_query = QueryGraph::from_parts(
        vec![
            QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C)),
            QueryAtom::new(AtomId::new(1), AtomSpec::new(Element::O)),
        ],
        vec![QueryBond::from_parts(
            Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            ),
            QueryNode::predicate(BondQueryPredicate::UnsupportedFeature(
                "Q87 unsupported bond callback input",
            )),
        )],
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("fixed unsupported bond query is structurally valid");
    let bond_topology = single_bond_topology(BondOrder::Single);
    let bond_target = SearchTarget::new(
        &bond_topology,
        &coordinates,
        &bond_topology.stereo_groups,
        None,
        None,
    );
    let bond_calls = Arc::new(AtomicUsize::new(0));
    let bond_calls_in_check = Arc::clone(&bond_calls);
    let bond_params = SubstructMatchParams {
        extra_bond_check: Some(Arc::new(move |_, _| {
            bond_calls_in_check.fetch_add(1, Ordering::SeqCst);
            true
        })),
        extra_bond_check_overrides_default_check: true,
        ..SubstructMatchParams::default()
    };
    assert_eq!(
        try_get_substruct_matches_with_params(&bond_target, &bond_query, &bond_params)
            .expect_err("unsupported bond query must fail before its callback"),
        SubstructMatchError::Unsupported {
            branch: "Q87 unsupported bond callback input",
            rdkit_function: "QueryBond::Match",
        }
    );
    assert_eq!(bond_calls.load(Ordering::SeqCst), 0);
}

#[test]
fn q88_complete_map_callback_accepts_rejects_and_stops_at_source_limit() {
    let query = single_bond_query(BondOrder::Single);
    let topology = two_site_carbon_oxygen_topology();
    let coordinates = CoordinateBlock::default();
    let target = SearchTarget::new(&topology, &coordinates, &topology.stereo_groups, None, None);

    let accept_calls = Arc::new(AtomicUsize::new(0));
    let accept_calls_in_check = Arc::clone(&accept_calls);
    let accept_params = SubstructMatchParams {
        extra_final_check: Some(Arc::new(move |target, mapping| {
            accept_calls_in_check.fetch_add(1, Ordering::SeqCst);
            assert_eq!(target.num_atoms(), 3);
            assert_eq!(mapping.len(), 2);
            assert_eq!(mapping[1], 1);
            assert!(matches!(mapping[0], 0 | 2));
            true
        })),
        ..SubstructMatchParams::default()
    };
    let accepted = get_substruct_matches_with_params(&target, &query, &accept_params);
    assert_eq!(accepted.len(), 2);
    assert_eq!(accept_calls.load(Ordering::SeqCst), 2);

    let selective_calls = Arc::new(AtomicUsize::new(0));
    let selective_calls_in_check = Arc::clone(&selective_calls);
    let selective_params = SubstructMatchParams {
        extra_final_check: Some(Arc::new(move |_, mapping| {
            selective_calls_in_check.fetch_add(1, Ordering::SeqCst);
            mapping[0] == 2
        })),
        ..SubstructMatchParams::default()
    };
    let selected = get_substruct_matches_with_params(&target, &query, &selective_params);
    assert_eq!(selected.len(), 1);
    assert_eq!(selected[0].atom_mapping, [2, 1]);
    assert_eq!(selective_calls.load(Ordering::SeqCst), 2);

    let reject_calls = Arc::new(AtomicUsize::new(0));
    let reject_calls_in_check = Arc::clone(&reject_calls);
    let reject_params = SubstructMatchParams {
        extra_final_check: Some(Arc::new(move |_, _| {
            reject_calls_in_check.fetch_add(1, Ordering::SeqCst);
            false
        })),
        ..SubstructMatchParams::default()
    };
    assert!(get_substruct_matches_with_params(&target, &query, &reject_params).is_empty());
    assert_eq!(reject_calls.load(Ordering::SeqCst), 2);

    let limited_calls = Arc::new(AtomicUsize::new(0));
    let limited_calls_in_check = Arc::clone(&limited_calls);
    let limited_params = SubstructMatchParams {
        max_matches: 1,
        extra_final_check: Some(Arc::new(move |_, _| {
            limited_calls_in_check.fetch_add(1, Ordering::SeqCst);
            true
        })),
        ..SubstructMatchParams::default()
    };
    assert_eq!(
        get_substruct_matches_with_params(&target, &query, &limited_params).len(),
        1
    );
    assert_eq!(limited_calls.load(Ordering::SeqCst), 1);
}

#[test]
fn q89_coordinate_callback_selects_conformers_tolerance_and_three_axes() {
    let mut query = single_atom_query(Element::C);
    query
        .add_conformer_3d(Conformer3D::new(7, vec![[0.0, 0.0, 0.0]], true))
        .expect("first fixed query conformer is valid");
    query
        .add_conformer_3d(Conformer3D::new(9, vec![[0.0, 0.0, 2.0]], false))
        .expect("second fixed query conformer is valid");
    let topology = single_atom_topology(Element::C);
    let coordinates = CoordinateBlock {
        conformers_2d: Vec::new(),
        conformers_3d: vec![
            Conformer3D::new(11, vec![[0.0, 0.0, 0.1]], true),
            Conformer3D::new(13, vec![[0.0, 0.0, 2.0]], false),
        ],
        source_coordinate_dim: None,
    };

    assert!(!coordinate_matches(
        &query,
        &topology,
        &coordinates,
        AtomCoordsMatchFunctor::default(),
    ));
    assert!(coordinate_matches(
        &query,
        &topology,
        &coordinates,
        AtomCoordsMatchFunctor::new(-1, -1, 0.11),
    ));
    assert!(!coordinate_matches(
        &query,
        &topology,
        &coordinates,
        AtomCoordsMatchFunctor::new(-1, -1, 0.09),
    ));

    assert!(coordinate_matches(
        &query,
        &topology,
        &coordinates,
        AtomCoordsMatchFunctor::new(13, 9, 0.0),
    ));
    assert!(!coordinate_matches(
        &query,
        &topology,
        &coordinates,
        AtomCoordsMatchFunctor::new(11, 9, 0.11),
    ));

    let query_without_coordinates = single_atom_query(Element::C);
    assert!(!coordinate_matches(
        &query_without_coordinates,
        &topology,
        &coordinates,
        AtomCoordsMatchFunctor::new(-1, -1, 1.0),
    ));
    assert!(!coordinate_matches(
        &query,
        &topology,
        &CoordinateBlock::default(),
        AtomCoordsMatchFunctor::new(-1, -1, 1.0),
    ));
}

#[test]
fn q90_json_update_preserves_source_types_defaults_and_unknown_keys() {
    let mut params = SubstructMatchParams {
        use_enhanced_stereo: true,
        max_recursive_matches: 77,
        aromatic_matches_single_or_double: true,
        ..SubstructMatchParams::default()
    };
    update_substruct_match_params_from_json(
        &mut params,
        r#"{
            "useChirality": 1,
            "aromaticMatchesConjugated": "true",
            "useQueryQueryMatches": false,
            "recursionPossible": "0",
            "uniquify": 0,
            "maxMatches": "42",
            "numThreads": "-3",
            "unknownSourceKey": {"ignored": true}
        }"#,
    )
    .expect("source-compatible scalar parameter forms update successfully");

    assert!(params.use_chirality);
    assert!(params.use_enhanced_stereo);
    assert!(params.aromatic_matches_conjugated);
    assert!(!params.use_query_query_matches);
    assert!(!params.recursion_possible);
    assert!(!params.uniquify);
    assert_eq!(params.max_matches, 42);
    assert_eq!(params.max_recursive_matches, 77);
    assert_eq!(params.num_threads, -3);
    assert!(params.aromatic_matches_single_or_double);

    let before_empty = params.clone();
    update_substruct_match_params_from_json(&mut params, "")
        .expect("the source empty-string update is a no-op");
    assert_eq!(params.use_chirality, before_empty.use_chirality);
    assert_eq!(params.max_matches, before_empty.max_matches);
    assert_eq!(params.num_threads, before_empty.num_threads);
}

#[test]
fn q90_json_error_keeps_source_ordered_partial_updates_and_integer_width() {
    let mut params = SubstructMatchParams {
        max_matches: 9,
        num_threads: 3,
        ..SubstructMatchParams::default()
    };
    let error = update_substruct_match_params_from_json(
        &mut params,
        r#"{
            "useChirality": true,
            "useEnhancedStereo": "1",
            "maxMatches": "bad",
            "numThreads": -2
        }"#,
    )
    .expect_err("a bad source unsigned field must report its name");
    assert!(matches!(
        error,
        SubstructMatchParamsJsonError::InvalidField {
            field: "maxMatches"
        }
    ));
    assert!(params.use_chirality);
    assert!(params.use_enhanced_stereo);
    assert_eq!(params.max_matches, 9);
    assert_eq!(params.num_threads, 3);

    let overflow_error = update_substruct_match_params_from_json(
        &mut params,
        r#"{"useQueryQueryMatches": 1, "maxMatches": 4294967296}"#,
    )
    .expect_err("source unsigned-int overflow must be rejected");
    assert!(matches!(
        overflow_error,
        SubstructMatchParamsJsonError::InvalidField {
            field: "maxMatches"
        }
    ));
    assert!(params.use_query_query_matches);
    assert_eq!(params.max_matches, 9);

    let malformed_error = update_substruct_match_params_from_json(&mut params, "{")
        .expect_err("malformed JSON fails before any field lookup");
    assert!(matches!(
        malformed_error,
        SubstructMatchParamsJsonError::InvalidJson(_)
    ));
}

#[test]
fn q91_json_serialization_preserves_source_order_types_and_roundtrip() {
    let params = SubstructMatchParams {
        use_chirality: true,
        use_enhanced_stereo: false,
        aromatic_matches_conjugated: true,
        use_query_query_matches: false,
        recursion_possible: false,
        uniquify: true,
        max_matches: 42,
        max_recursive_matches: 77,
        num_threads: -3,
        specified_stereo_query_matches_unspecified: true,
        aromatic_matches_single_or_double: false,
        ..SubstructMatchParams::default()
    };

    let json = substruct_match_params_to_json(&params);
    assert_eq!(
        json,
        concat!(
            "{\n",
            "    \"useChirality\": \"true\",\n",
            "    \"useEnhancedStereo\": \"false\",\n",
            "    \"aromaticMatchesConjugated\": \"true\",\n",
            "    \"useQueryQueryMatches\": \"false\",\n",
            "    \"recursionPossible\": \"false\",\n",
            "    \"uniquify\": \"true\",\n",
            "    \"maxMatches\": \"42\",\n",
            "    \"maxRecursiveMatches\": \"77\",\n",
            "    \"numThreads\": \"-3\",\n",
            "    \"specifiedStereoQueryMatchesUnspecified\": \"true\",\n",
            "    \"aromaticMatchesSingleOrDouble\": \"false\"\n",
            "}\n",
        )
    );

    let parsed: serde_json::Value = serde_json::from_str(&json).expect("source JSON is valid");
    assert!(
        parsed
            .as_object()
            .expect("source output is an object")
            .values()
            .all(serde_json::Value::is_string),
        "Boost property-tree writes every scalar value as a JSON string",
    );

    let mut roundtrip = SubstructMatchParams::default();
    update_substruct_match_params_from_json(&mut roundtrip, &json)
        .expect("source serialization must parse through the source update route");
    assert_eq!(roundtrip.use_chirality, params.use_chirality);
    assert_eq!(roundtrip.use_enhanced_stereo, params.use_enhanced_stereo);
    assert_eq!(
        roundtrip.aromatic_matches_conjugated,
        params.aromatic_matches_conjugated
    );
    assert_eq!(
        roundtrip.use_query_query_matches,
        params.use_query_query_matches
    );
    assert_eq!(roundtrip.recursion_possible, params.recursion_possible);
    assert_eq!(roundtrip.uniquify, params.uniquify);
    assert_eq!(roundtrip.max_matches, params.max_matches);
    assert_eq!(
        roundtrip.max_recursive_matches,
        params.max_recursive_matches
    );
    assert_eq!(roundtrip.num_threads, params.num_threads);
    assert_eq!(
        roundtrip.specified_stereo_query_matches_unspecified,
        params.specified_stereo_query_matches_unspecified
    );
    assert_eq!(
        roundtrip.aromatic_matches_single_or_double,
        params.aromatic_matches_single_or_double
    );
}
