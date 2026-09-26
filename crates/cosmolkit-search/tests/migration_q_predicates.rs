use std::sync::{
    Arc,
    atomic::{AtomicUsize, Ordering},
};

use cosmolkit_core::{RingInfo, ValenceAssignment, find_sssr_from_parts};
use cosmolkit_model::{
    Atom, AtomId, AtomQueryPredicate, AtomRangeBounds, AtomRangeDataFunction, AtomRangeQuery,
    AtomSpec, Bond, BondId, BondQueryPredicate, BondSpec, CoordinateBlock, QueryAtom,
    QueryAtomIdentity, QueryBond, QueryGraph, QueryNode, TopologyBlock,
};
use cosmolkit_search::{
    ExtraAtomCheck, SearchTarget, SubstructMatchError, SubstructMatchParams,
    get_substruct_matches_with_params, try_get_substruct_matches_with_params,
};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo, Element};

fn atom(element: Element, formal_charge: i8, isotope: Option<u16>, radical_electrons: u8) -> Atom {
    let mut spec = AtomSpec::new(element)
        .with_formal_charge(formal_charge)
        .with_radical_electrons(radical_electrons);
    if let Some(isotope) = isotope {
        spec = spec.with_isotope(isotope);
    }
    Atom::from_spec(AtomId::new(0), spec)
}

fn query_graph(atom: QueryAtom) -> QueryGraph {
    QueryGraph::from_parts(
        vec![atom],
        Vec::new(),
        Default::default(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("single-atom query graph is valid")
}

fn topology(atom: Atom) -> TopologyBlock {
    TopologyBlock::try_from_parts(vec![atom], Vec::new(), Vec::new(), Vec::new())
        .expect("single-atom target topology is valid")
}

fn topology_with_edges(atom_count: usize, edges: &[(usize, usize)]) -> TopologyBlock {
    let atoms = (0..atom_count)
        .map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C)))
        .collect();
    let bonds = edges
        .iter()
        .enumerate()
        .map(|(index, &(begin, end))| {
            Bond::from_spec(
                BondId::new(index),
                BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
            )
        })
        .collect();
    TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
        .expect("test graph topology is valid")
}

fn q38_target_with_bond(
    order: BondOrder,
    aromatic: bool,
    direction: BondDirection,
    stereo: BondStereo,
) -> TopologyBlock {
    let mut atoms = (0..4)
        .map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C)))
        .collect::<Vec<_>>();
    atoms[0]
        .set_prop("q37-side", "left")
        .expect("target property is valid");
    atoms[1]
        .set_prop("q37-side", "right")
        .expect("target property is valid");

    let mut spec = BondSpec::new(AtomId::new(0), AtomId::new(1), order)
        .with_aromatic(aromatic)
        .with_direction(direction)
        .with_stereo(stereo);
    if matches!(stereo, BondStereo::Cis | BondStereo::Trans) {
        spec = spec.with_stereo_atoms(AtomId::new(2), AtomId::new(3));
    }
    let bond = Bond::from_spec(BondId::new(0), spec);
    TopologyBlock::try_from_parts(atoms, vec![bond], Vec::new(), Vec::new())
        .expect("Q38 target topology is valid")
}

fn cycle_topology(atom_count: usize) -> TopologyBlock {
    assert!(atom_count >= 3);
    let edges = (0..atom_count)
        .map(|index| (index, (index + 1) % atom_count))
        .collect::<Vec<_>>();
    topology_with_edges(atom_count, &edges)
}

fn query_atom_predicate(predicate: AtomQueryPredicate) -> QueryGraph {
    query_graph(QueryAtom::from_parts(
        atom(Element::C, 0, None, 0),
        QueryNode::predicate(predicate),
    ))
}

fn query_bond_predicate(predicate: BondQueryPredicate) -> QueryGraph {
    let mut left = QueryAtom::from_parts(
        Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
    );
    left.set_prop("q37-side", "left")
        .expect("query property is valid");
    let mut right = QueryAtom::from_parts(
        Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
    );
    right
        .set_prop("q37-side", "right")
        .expect("query property is valid");
    let bond = Bond::from_spec(
        BondId::new(0),
        BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
    );
    QueryGraph::from_parts(
        vec![left, right],
        vec![QueryBond::from_parts(bond, QueryNode::predicate(predicate))],
        Default::default(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("single-bond query graph is valid")
}

fn matches_on_topology(
    query: &QueryGraph,
    topology: &TopologyBlock,
    ring_info: Option<&RingInfo>,
) -> bool {
    matches_on_topology_with_params(query, topology, ring_info, &SubstructMatchParams::default())
}

fn matches_on_topology_with_params(
    query: &QueryGraph,
    topology: &TopologyBlock,
    ring_info: Option<&RingInfo>,
    params: &SubstructMatchParams,
) -> bool {
    let coordinates = CoordinateBlock::default();
    let target = SearchTarget::new(
        topology,
        &coordinates,
        &topology.stereo_groups,
        ring_info,
        None,
    );
    !get_substruct_matches_with_params(&target, query, params).is_empty()
}

fn matches(query: &QueryGraph, target_atom: Atom, params: &SubstructMatchParams) -> bool {
    let topology = topology(target_atom);
    let coordinates = CoordinateBlock::default();
    let target = SearchTarget::new(&topology, &coordinates, &topology.stereo_groups, None, None);
    !get_substruct_matches_with_params(&target, query, params).is_empty()
}

fn matches_with_valence(
    query: &QueryGraph,
    target_atom: Atom,
    valence: &ValenceAssignment,
    params: &SubstructMatchParams,
) -> bool {
    let topology = topology(target_atom);
    let coordinates = CoordinateBlock::default();
    let target = SearchTarget::new(
        &topology,
        &coordinates,
        &topology.stereo_groups,
        None,
        Some(valence),
    );
    !get_substruct_matches_with_params(&target, query, params).is_empty()
}

fn try_matches(
    query: &QueryGraph,
    target_atom: Atom,
    params: &SubstructMatchParams,
) -> Result<bool, SubstructMatchError> {
    let topology = topology(target_atom);
    let coordinates = CoordinateBlock::default();
    let target = SearchTarget::new(&topology, &coordinates, &topology.stereo_groups, None, None);
    Ok(!try_get_substruct_matches_with_params(&target, query, params)?.is_empty())
}

fn counted_atom_check(calls: Arc<AtomicUsize>, result: bool) -> ExtraAtomCheck {
    Arc::new(
        move |_: &QueryGraph, _: &QueryAtom, _: &SearchTarget<'_>, _: &Atom| {
            calls.fetch_add(1, Ordering::SeqCst);
            result
        },
    )
}

#[test]
fn q33_carrier_derived_and_explicit_rows_use_their_source_dispatch() {
    let carrier = QueryAtom::from_carrier_parts(
        atom(Element::C, 0, None, 0),
        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
    )
    .with_identity(QueryAtomIdentity::AtomicNumber(8));
    assert!(carrier.predicate_is_carrier_derived());
    let carrier_query = query_graph(carrier);
    let params = SubstructMatchParams::default();
    assert!(matches(
        &carrier_query,
        atom(Element::O, 0, None, 0),
        &params
    ));
    assert!(!matches(
        &carrier_query,
        atom(Element::C, 0, None, 0),
        &params
    ));

    let explicit = QueryAtom::from_identity_parts(
        AtomId::new(0),
        QueryAtomIdentity::AtomicNumber(8),
        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
    );
    assert!(!explicit.predicate_is_carrier_derived());
    let explicit_query = query_graph(explicit);
    assert!(matches(
        &explicit_query,
        atom(Element::C, 0, None, 0),
        &params
    ));
    assert!(!matches(
        &explicit_query,
        atom(Element::O, 0, None, 0),
        &params
    ));
}

#[test]
fn q33_explicit_negation_and_carrier_default_and_nondefault_fields_match_source() {
    let negated = QueryAtom::from_parts(
        atom(Element::C, 0, None, 0),
        QueryNode::not(QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6))),
    );
    let negated_query = query_graph(negated);
    let params = SubstructMatchParams::default();
    assert!(matches(
        &negated_query,
        atom(Element::O, 0, None, 0),
        &params
    ));
    assert!(!matches(
        &negated_query,
        atom(Element::C, 0, None, 0),
        &params
    ));

    let default_carrier = QueryAtom::from_carrier_parts(
        atom(Element::C, 0, None, 0),
        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
    );
    let default_query = query_graph(default_carrier);
    assert!(matches(
        &default_query,
        atom(Element::C, 1, Some(13), 1),
        &params
    ));

    let constrained_carrier = QueryAtom::from_carrier_parts(
        atom(Element::C, 1, Some(13), 1),
        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
    );
    let constrained_query = query_graph(constrained_carrier);
    assert!(matches(
        &constrained_query,
        atom(Element::C, 1, Some(13), 1),
        &params
    ));
    for target_atom in [
        atom(Element::C, 0, Some(13), 1),
        atom(Element::C, 1, Some(12), 1),
        atom(Element::C, 1, Some(13), 0),
    ] {
        assert!(!matches(&constrained_query, target_atom, &params));
    }
}

#[test]
fn q33_override_and_post_callback_follow_source_property_order() {
    let mut query_atom = QueryAtom::from_parts(
        atom(Element::O, 0, None, 0),
        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(8)),
    );
    query_atom
        .set_prop("gate", "query")
        .expect("query property is valid");
    let query = query_graph(query_atom);
    let mut target_atom = atom(Element::C, 0, None, 0);
    target_atom
        .set_prop("gate", "target")
        .expect("target property is valid");

    let calls = Arc::new(AtomicUsize::new(0));
    let mut params = SubstructMatchParams::default();
    params.atom_properties.push("gate".to_owned());
    params.extra_atom_check = Some(counted_atom_check(Arc::clone(&calls), false));

    // The explicit query fails first, so neither properties nor the post-check run.
    assert!(!matches(&query, target_atom.clone(), &params));
    assert_eq!(calls.load(Ordering::SeqCst), 0);

    // A passing default query reaches property comparison; a property mismatch
    // returns before the post-check callback.
    let mut passing_query_atom = QueryAtom::from_parts(
        atom(Element::C, 0, None, 0),
        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
    );
    passing_query_atom
        .set_prop("gate", "query")
        .expect("query property is valid");
    let passing_query = query_graph(passing_query_atom);
    assert!(!matches(&passing_query, target_atom.clone(), &params));
    assert_eq!(calls.load(Ordering::SeqCst), 0);

    let mut matching_target = atom(Element::C, 0, None, 0);
    matching_target
        .set_prop("gate", "query")
        .expect("target property is valid");
    assert!(!matches(&passing_query, matching_target, &params));
    assert_eq!(calls.load(Ordering::SeqCst), 1);

    // The source override callback returns directly before the default match
    // and requested property checks.
    let override_calls = Arc::new(AtomicUsize::new(0));
    let mut override_params = SubstructMatchParams::default();
    override_params.atom_properties.push("gate".to_owned());
    override_params.extra_atom_check = Some(counted_atom_check(Arc::clone(&override_calls), true));
    override_params.extra_atom_check_overrides_default_check = true;
    assert!(matches(&query, target_atom, &override_params));
    assert_eq!(override_calls.load(Ordering::SeqCst), 1);
}

#[test]
fn q34_mass_predicates_use_only_the_source_unknown_isotope_fallback() {
    let params = SubstructMatchParams::default();
    for (element, source_mass, wrong_mass) in [(Element::C, 999, 1_000), (Element::DUMMY, 0, 999)] {
        let target_atom = atom(element, 0, Some(999), 0);
        let query = query_graph(QueryAtom::from_parts(
            atom(element, 0, None, 0),
            QueryNode::predicate(AtomQueryPredicate::Mass(source_mass)),
        ));
        assert!(
            try_matches(&query, target_atom.clone(), &params)
                .expect("source mass fallback is a successful query evaluation")
        );

        let wrong_query = query_graph(QueryAtom::from_parts(
            atom(element, 0, None, 0),
            QueryNode::predicate(AtomQueryPredicate::Mass(wrong_mass)),
        ));
        assert!(
            !try_matches(&wrong_query, target_atom, &params)
                .expect("a nonmatching mass predicate is not a lookup error")
        );
    }
}

#[test]
fn q35_degree_hydrogen_and_valence_predicates_use_current_topology() {
    let target_atom = atom(Element::C, 0, None, 0);
    let stale_valence = ValenceAssignment {
        explicit_valence: vec![7],
        implicit_hydrogens: vec![0],
    };
    let params = SubstructMatchParams::default();

    // RDKit derives these values while updating the current atom property
    // cache, then QueryOps reads that cache for degree, H, and valence tests.
    // An unbound detached value must not override the isolated carbon topology.
    for predicate in [
        AtomQueryPredicate::HydrogenCount(4),
        AtomQueryPredicate::HasImplicitHydrogen,
        AtomQueryPredicate::ImplicitValence(4),
        AtomQueryPredicate::ExplicitValence(0),
        AtomQueryPredicate::TotalDegree(4),
        AtomQueryPredicate::TotalValence(4),
    ] {
        let query = query_graph(QueryAtom::from_parts(
            atom(Element::C, 0, None, 0),
            QueryNode::predicate(predicate),
        ));
        assert!(matches_with_valence(
            &query,
            target_atom.clone(),
            &stale_valence,
            &params,
        ));
    }

    for predicate in [
        AtomQueryPredicate::HydrogenCount(0),
        AtomQueryPredicate::ImplicitValence(0),
        AtomQueryPredicate::ExplicitValence(7),
        AtomQueryPredicate::TotalDegree(0),
        AtomQueryPredicate::TotalValence(7),
    ] {
        let query = query_graph(QueryAtom::from_parts(
            atom(Element::C, 0, None, 0),
            QueryNode::predicate(predicate),
        ));
        assert!(!matches_with_valence(
            &query,
            target_atom.clone(),
            &stale_valence,
            &params,
        ));
    }
}

#[test]
fn q36_ring_queries_rebuild_current_state_and_preserve_integer_boundaries() {
    let triangle = cycle_topology(3);
    for predicate in [
        AtomQueryPredicate::InRing,
        AtomQueryPredicate::InRingOfSize(3),
        AtomQueryPredicate::NumAtomRings(1),
        AtomQueryPredicate::RingBondCount(2),
        AtomQueryPredicate::SmallestRingSize(3),
    ] {
        assert!(matches_on_topology(
            &query_atom_predicate(predicate),
            &triangle,
            None,
        ));
    }

    let triangle_rings =
        find_sssr_from_parts(triangle.atoms.len(), &triangle.bonds, &triangle.adjacency)
            .expect("triangle SSSR initializes");
    let path = topology_with_edges(3, &[(0, 1), (1, 2)]);
    assert!(!matches_on_topology(
        &query_atom_predicate(AtomQueryPredicate::InRing),
        &path,
        Some(&triangle_rings),
    ));
    assert!(matches_on_topology(
        &query_atom_predicate(AtomQueryPredicate::RingBondCount(0)),
        &path,
        Some(&triangle_rings),
    ));

    // The source range helper returns -1 when no ring size meets a lower-only
    // bound, and INT_MAX when no ring size meets an upper-only bound.
    for predicate in [
        AtomQueryPredicate::InRingOfSizeLessEqual(3),
        AtomQueryPredicate::InRingOfSizeGreaterEqual(3),
    ] {
        assert!(matches_on_topology(
            &query_atom_predicate(predicate),
            &path,
            None,
        ));
    }

    let large_cycle = cycle_topology(256);
    for predicate in [
        AtomQueryPredicate::SmallestRingSize(256),
        AtomQueryPredicate::SmallestRingSizeGreaterEqual(255),
        AtomQueryPredicate::InRingOfSize(256),
    ] {
        assert!(matches_on_topology(
            &query_atom_predicate(predicate),
            &large_cycle,
            None,
        ));
    }
    assert!(!matches_on_topology(
        &query_atom_predicate(AtomQueryPredicate::SmallestRingSizeLessEqual(255)),
        &large_cycle,
        None,
    ));

    let mut flower_edges = Vec::with_capacity(256 * 3);
    for cycle in 0..256 {
        let first = 1 + cycle * 2;
        let second = first + 1;
        flower_edges.extend([(0, first), (first, second), (second, 0)]);
    }
    let flower = topology_with_edges(513, &flower_edges);
    assert!(matches_on_topology(
        &query_atom_predicate(AtomQueryPredicate::NumAtomRings(256)),
        &flower,
        None,
    ));
    assert!(matches_on_topology(
        &query_atom_predicate(AtomQueryPredicate::RingBondCount(512)),
        &flower,
        None,
    ));
    assert!(!matches_on_topology(
        &query_atom_predicate(AtomQueryPredicate::RingBondCountLessEqual(1)),
        &flower,
        None,
    ));
}

#[test]
fn q37_bond_ring_queries_initialize_absent_state_and_keep_fused_ring_counts() {
    let acyclic = topology_with_edges(2, &[(0, 1)]);
    for predicate in [
        BondQueryPredicate::IsInRing(false),
        BondQueryPredicate::NumRingBonds(0),
        BondQueryPredicate::MinRingSize(0),
    ] {
        let query = query_bond_predicate(predicate);
        let mut params = SubstructMatchParams::default();
        params.atom_properties.push("q37-side".to_owned());
        let mut target = acyclic.clone();
        target.atoms[0]
            .set_prop("q37-side", "left")
            .expect("target property is valid");
        target.atoms[1]
            .set_prop("q37-side", "right")
            .expect("target property is valid");
        assert!(matches_on_topology_with_params(
            &query, &target, None, &params,
        ));
    }

    let fused = topology_with_edges(6, &[(0, 1), (0, 2), (2, 3), (3, 1), (0, 4), (4, 5), (5, 1)]);
    let fused_rings = find_sssr_from_parts(fused.atoms.len(), &fused.bonds, &fused.adjacency)
        .expect("fused-ring SSSR initializes");
    assert_eq!(fused_rings.num_bond_rings(BondId::new(0)), 2);
    assert_eq!(fused_rings.min_bond_ring_size(BondId::new(0)), 4);

    let mut fused = fused;
    fused.atoms[0]
        .set_prop("q37-side", "left")
        .expect("target property is valid");
    fused.atoms[1]
        .set_prop("q37-side", "right")
        .expect("target property is valid");
    let mut params = SubstructMatchParams::default();
    params.atom_properties.push("q37-side".to_owned());
    for predicate in [
        BondQueryPredicate::IsInRing(true),
        BondQueryPredicate::NumRingBonds(2),
        BondQueryPredicate::InRingOfSize(4),
        BondQueryPredicate::MinRingSize(4),
        BondQueryPredicate::NumRingBondsGreaterEqual(2),
    ] {
        assert!(matches_on_topology_with_params(
            &query_bond_predicate(predicate),
            &fused,
            None,
            &params,
        ));
    }
    assert!(!matches_on_topology_with_params(
        &query_bond_predicate(BondQueryPredicate::NumRingBonds(1)),
        &fused,
        None,
        &params,
    ));
    assert!(!matches_on_topology_with_params(
        &query_bond_predicate(BondQueryPredicate::NumRingBondsLessEqual(1)),
        &fused,
        None,
        &params,
    ));

    let mut shared_edge_cycles = vec![(0, 1)];
    for cycle in 0..256 {
        let first = 2 + cycle * 2;
        let second = first + 1;
        shared_edge_cycles.extend([(0, first), (first, second), (second, 1)]);
    }
    let mut many_fused = topology_with_edges(514, &shared_edge_cycles);
    many_fused.atoms[0]
        .set_prop("q37-side", "left")
        .expect("target property is valid");
    many_fused.atoms[1]
        .set_prop("q37-side", "right")
        .expect("target property is valid");
    let many_ring_info = find_sssr_from_parts(
        many_fused.atoms.len(),
        &many_fused.bonds,
        &many_fused.adjacency,
    )
    .expect("multi-fused-ring SSSR initializes");
    assert_eq!(many_ring_info.num_bond_rings(BondId::new(0)), 256);
    for predicate in [
        BondQueryPredicate::NumRingBonds(256),
        BondQueryPredicate::NumRingBondsGreaterEqual(255),
        BondQueryPredicate::InRingOfSize(4),
        BondQueryPredicate::MinRingSize(4),
    ] {
        assert!(matches_on_topology_with_params(
            &query_bond_predicate(predicate),
            &many_fused,
            None,
            &params,
        ));
    }
    assert!(!matches_on_topology_with_params(
        &query_bond_predicate(BondQueryPredicate::NumRingBondsLessEqual(255)),
        &many_fused,
        None,
        &params,
    ));
}

#[test]
fn q38_bond_order_sets_aromatic_flag_direction_and_stereo_follow_source() {
    let mut params = SubstructMatchParams::default();
    params.atom_properties.push("q37-side".to_owned());
    let matches_bond = |target: &TopologyBlock, predicate| {
        matches_on_topology_with_params(&query_bond_predicate(predicate), target, None, &params)
    };

    // QueryOps.h compares the bond type for each order set. The aromatic
    // boolean is a separate Bond field and cannot substitute for that type.
    let order_sets = [
        vec![BondOrder::Single, BondOrder::Aromatic],
        vec![BondOrder::Double, BondOrder::Aromatic],
        vec![BondOrder::Single, BondOrder::Double],
        vec![BondOrder::Single, BondOrder::Double, BondOrder::Aromatic],
    ];
    let order_cases = [
        (BondOrder::Single, [true, false, true, true]),
        (BondOrder::Double, [false, true, true, true]),
        (BondOrder::Aromatic, [true, true, false, true]),
        (BondOrder::Triple, [false, false, false, false]),
    ];
    for (order, expected_set_matches) in order_cases {
        let target = q38_target_with_bond(order, false, BondDirection::None, BondStereo::None);
        assert!(matches_bond(&target, BondQueryPredicate::Order(order)));
        for (set, expected) in order_sets.iter().zip(expected_set_matches) {
            assert_eq!(
                matches_bond(&target, BondQueryPredicate::OrderIn(set.clone())),
                expected,
                "bond order {order:?} against set {set:?}"
            );
        }
    }

    let flagged_single = q38_target_with_bond(
        BondOrder::Single,
        true,
        BondDirection::None,
        BondStereo::None,
    );
    assert!(matches_bond(
        &flagged_single,
        BondQueryPredicate::Order(BondOrder::Single)
    ));
    assert!(!matches_bond(
        &flagged_single,
        BondQueryPredicate::Order(BondOrder::Aromatic)
    ));
    assert!(matches_bond(
        &flagged_single,
        BondQueryPredicate::IsAromatic(true)
    ));
    assert!(!matches_bond(
        &flagged_single,
        BondQueryPredicate::IsAromatic(false)
    ));

    let unflagged_aromatic = q38_target_with_bond(
        BondOrder::Aromatic,
        false,
        BondDirection::None,
        BondStereo::None,
    );
    assert!(matches_bond(
        &unflagged_aromatic,
        BondQueryPredicate::Order(BondOrder::Aromatic)
    ));
    assert!(!matches_bond(
        &unflagged_aromatic,
        BondQueryPredicate::IsAromatic(true)
    ));
    assert!(matches_bond(
        &unflagged_aromatic,
        BondQueryPredicate::IsAromatic(false)
    ));

    let directions = [
        BondDirection::None,
        BondDirection::BeginWedge,
        BondDirection::BeginDash,
        BondDirection::EndDownRight,
        BondDirection::EndUpRight,
        BondDirection::EitherDouble,
        BondDirection::Unknown,
    ];
    for (index, direction) in directions.iter().copied().enumerate() {
        let target = q38_target_with_bond(BondOrder::Single, false, direction, BondStereo::None);
        assert!(matches_bond(
            &target,
            BondQueryPredicate::Direction(direction)
        ));
        let different = directions[(index + 1) % directions.len()];
        assert!(!matches_bond(
            &target,
            BondQueryPredicate::Direction(different)
        ));
    }

    let stereo_values = [
        BondStereo::None,
        BondStereo::Any,
        BondStereo::Z,
        BondStereo::E,
        BondStereo::Cis,
        BondStereo::Trans,
        BondStereo::AtropCw,
        BondStereo::AtropCcw,
    ];
    for (index, stereo) in stereo_values.iter().copied().enumerate() {
        let target = q38_target_with_bond(BondOrder::Single, false, BondDirection::None, stereo);
        assert!(matches_bond(&target, BondQueryPredicate::Stereo(stereo)));
        let different = stereo_values[(index + 1) % stereo_values.len()];
        assert!(!matches_bond(
            &target,
            BondQueryPredicate::Stereo(different)
        ));
        assert_eq!(
            matches_bond(&target, BondQueryPredicate::HasStereo),
            stereo != BondStereo::None,
            "RDKit queryBondHasStereo is false only for STEREONONE"
        );
    }
}

#[test]
fn q39_scalar_ranges_tolerance_negation_sets_and_errors_follow_source() {
    let params = SubstructMatchParams::default();
    let range_matches = |bounds: AtomRangeBounds, formal_charge: i8| {
        let query = query_atom_predicate(AtomQueryPredicate::Range(AtomRangeQuery::new(
            bounds,
            AtomRangeDataFunction::FormalCharge,
        )));
        try_matches(&query, atom(Element::C, formal_charge, None, 0), &params)
            .expect("supported scalar range evaluation succeeds")
    };

    // Query's base tolerance defaults to zero, and AtomRangeQuery carries no
    // tolerance field. Pinned LessEqualQuery and GreaterEqualQuery compare
    // their stored threshold with observed data in that order.
    let charges = [-2, -1, 0, 1, 2];
    let range_cases = [
        (
            AtomRangeBounds::LessEqual(0),
            [false, false, true, true, true],
        ),
        (
            AtomRangeBounds::GreaterEqual(0),
            [true, true, true, false, false],
        ),
        (
            AtomRangeBounds::Inclusive {
                lower: -1,
                upper: 1,
                lower_open: true,
                upper_open: true,
            },
            [false, false, true, false, false],
        ),
        (
            AtomRangeBounds::Inclusive {
                lower: -1,
                upper: 1,
                lower_open: true,
                upper_open: false,
            },
            [false, false, true, true, false],
        ),
        (
            AtomRangeBounds::Inclusive {
                lower: -1,
                upper: 1,
                lower_open: false,
                upper_open: true,
            },
            [false, true, true, false, false],
        ),
        (
            AtomRangeBounds::Inclusive {
                lower: -1,
                upper: 1,
                lower_open: false,
                upper_open: false,
            },
            [false, true, true, true, false],
        ),
        (
            AtomRangeBounds::Inclusive {
                lower: 1,
                upper: -1,
                lower_open: false,
                upper_open: false,
            },
            [false; 5],
        ),
    ];
    for (bounds, expected) in range_cases {
        for (formal_charge, expected) in charges.into_iter().zip(expected) {
            assert_eq!(
                range_matches(bounds, formal_charge),
                expected,
                "formal charge {formal_charge} against {bounds:?}"
            );
        }
    }

    let closed_range = AtomRangeQuery::new(
        AtomRangeBounds::Inclusive {
            lower: 0,
            upper: 1,
            lower_open: false,
            upper_open: false,
        },
        AtomRangeDataFunction::FormalCharge,
    );
    let negated = query_graph(QueryAtom::from_parts(
        atom(Element::C, 0, None, 0),
        QueryNode::not(QueryNode::predicate(AtomQueryPredicate::Range(
            closed_range,
        ))),
    ));
    for (formal_charge, expected) in [(-1, true), (0, false), (1, false), (2, true)] {
        assert_eq!(
            try_matches(&negated, atom(Element::C, formal_charge, None, 0), &params)
                .expect("negated supported range evaluation succeeds"),
            expected,
            "negated formal-charge range at {formal_charge}"
        );
    }

    // EqualityQuery's source default tolerance is zero: one atomic-number
    // step away does not match. Set membership retains source set semantics.
    let equal_carbon = query_atom_predicate(AtomQueryPredicate::AtomicNumber(6));
    assert!(matches(
        &equal_carbon,
        atom(Element::C, 0, None, 0),
        &params
    ));
    assert!(!matches(
        &equal_carbon,
        atom(Element::N, 0, None, 0),
        &params
    ));

    let in_set = query_atom_predicate(AtomQueryPredicate::AtomicNumberIn(vec![8, 6, 6]));
    let not_in_set = query_atom_predicate(AtomQueryPredicate::AtomicNumberNotIn(vec![6, 8]));
    assert!(matches(&in_set, atom(Element::C, 0, None, 0), &params));
    assert!(matches(&in_set, atom(Element::O, 0, None, 0), &params));
    assert!(!matches(&in_set, atom(Element::N, 0, None, 0), &params));
    assert!(!matches(&not_in_set, atom(Element::C, 0, None, 0), &params));
    assert!(matches(&not_in_set, atom(Element::N, 0, None, 0), &params));

    // The try matcher preserves the first structured unsupported-leaf error
    // through a Boolean query tree instead of turning it into a scalar miss.
    let range_then_unsupported = query_graph(QueryAtom::from_parts(
        atom(Element::C, 0, None, 0),
        QueryNode::and(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(7)),
            QueryNode::predicate(AtomQueryPredicate::UnsupportedFeature(
                "Q39 error propagation probe",
            )),
        ]),
    ));
    assert_eq!(
        try_matches(
            &range_then_unsupported,
            atom(Element::C, 0, None, 0),
            &params,
        ),
        Err(SubstructMatchError::Unsupported {
            branch: "Q39 error propagation probe",
            rdkit_function: "QueryAtom::Match",
        })
    );
}
