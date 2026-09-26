use std::collections::{BTreeMap, BTreeSet};
use std::sync::{
    Arc,
    atomic::{AtomicUsize, Ordering},
};

use cosmolkit_model::{
    Atom, AtomId, AtomQueryPredicate, AtomSpec, Bond, BondId, BondQueryPredicate, BondSpec,
    CoordinateBlock, QueryAtom, QueryBond, QueryGraph, QueryNode, QueryStateError, QueryStateRef,
    RecursiveStructureQuery, StereoGroup, StereoGroupKind, TopologyBlock, remap_query_rows,
};
use cosmolkit_search::{
    SearchTarget, SmartsParseParams, SubstructMatchError, SubstructMatchParams, compile_query,
    get_substruct_match, get_substruct_matches, get_substruct_matches_with_params, parse_smarts,
    try_get_substruct_matches_with_params,
};
use cosmolkit_types::{BondOrder, BondStereo, ChiralTag, Element};

fn topology(elements: &[Element], edge_endpoints: &[(usize, usize)]) -> TopologyBlock {
    let atoms = elements
        .iter()
        .copied()
        .enumerate()
        .map(|(index, element)| Atom::from_spec(AtomId::new(index), AtomSpec::new(element)))
        .collect();
    let bonds = edge_endpoints
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
        .expect("fixed Q-match target topology is valid")
}

fn topology_with_isotope(
    elements: &[Element],
    edge_endpoints: &[(usize, usize)],
    isotope: u16,
) -> TopologyBlock {
    let atoms = elements
        .iter()
        .copied()
        .enumerate()
        .map(|(index, element)| {
            Atom::from_spec(
                AtomId::new(index),
                AtomSpec::new(element).with_isotope(isotope),
            )
        })
        .collect();
    let bonds = edge_endpoints
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
        .expect("fixed isotope-bearing Q-match topology is structurally valid")
}

fn q80_tetrahedral_target(
    chiral_tag: ChiralTag,
    edge_endpoints: &[(usize, usize)],
    no_implicit: bool,
) -> TopologyBlock {
    let atom_specs = [
        AtomSpec::new(Element::C)
            .with_no_implicit(no_implicit)
            .with_chiral_tag(chiral_tag),
        AtomSpec::new(Element::F),
        AtomSpec::new(Element::CL),
        AtomSpec::new(Element::BR),
    ];
    let atoms = atom_specs
        .into_iter()
        .enumerate()
        .map(|(index, spec)| Atom::from_spec(AtomId::new(index), spec))
        .collect();
    let bonds = edge_endpoints
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
        .expect("fixed Q80 tetrahedral target topology is valid")
}

fn q80_tetrahedral_query(chiral_tag: ChiralTag) -> QueryGraph {
    // smarts.yy's H_TOKEN reduction stores this carrier state and an
    // AtomHCount(1) query leaf. Construct that detached value directly so the
    // regression isolates SubstructMatch's source final-check boundary.
    let center = QueryAtom::from_parts(
        Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_explicit_hydrogens(1)
                .with_no_implicit(true)
                .with_chiral_tag(chiral_tag),
        ),
        QueryNode::and(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomType {
                atomic_number: 6,
                aromatic: false,
            }),
            QueryNode::predicate(AtomQueryPredicate::HydrogenCount(1)),
        ]),
    );
    let atoms = vec![
        center,
        QueryAtom::new(AtomId::new(1), AtomSpec::new(Element::F)),
        QueryAtom::new(AtomId::new(2), AtomSpec::new(Element::CL)),
        QueryAtom::new(AtomId::new(3), AtomSpec::new(Element::BR)),
    ];
    let bonds = [(0, 1), (0, 2), (0, 3)]
        .into_iter()
        .enumerate()
        .map(|(index, (begin, end))| {
            QueryBond::new(
                BondId::new(index),
                BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
            )
        })
        .collect();
    QueryGraph::from_parts(
        atoms,
        bonds,
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("fixed Q80 tetrahedral query graph is valid")
}

fn q81_stereo_bond_spec(
    begin: usize,
    end: usize,
    order: BondOrder,
    stereo: BondStereo,
    stereo_atoms: Option<[usize; 2]>,
) -> BondSpec {
    let spec = BondSpec::new(AtomId::new(begin), AtomId::new(end), order).with_stereo(stereo);
    if let Some([begin_ref, end_ref]) = stereo_atoms {
        spec.with_stereo_atoms(AtomId::new(begin_ref), AtomId::new(end_ref))
    } else {
        spec
    }
}

fn q81_stereo_query(stereo: BondStereo, with_stereo_atoms: bool) -> QueryGraph {
    let atoms = [Element::F, Element::C, Element::C, Element::CL]
        .into_iter()
        .enumerate()
        .map(|(index, element)| QueryAtom::new(AtomId::new(index), AtomSpec::new(element)))
        .collect();
    let stereo_atoms = with_stereo_atoms.then_some([0, 3]);
    let bonds = vec![
        QueryBond::new(
            BondId::new(0),
            q81_stereo_bond_spec(0, 1, BondOrder::Single, BondStereo::None, None),
        ),
        QueryBond::new(
            BondId::new(1),
            q81_stereo_bond_spec(1, 2, BondOrder::Double, stereo, stereo_atoms),
        ),
        QueryBond::new(
            BondId::new(2),
            q81_stereo_bond_spec(2, 3, BondOrder::Single, BondStereo::None, None),
        ),
    ];
    QueryGraph::from_parts(
        atoms,
        bonds,
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("fixed Q81 double-bond stereo query is structurally valid")
}

fn q81_stereo_target(
    stereo: BondStereo,
    stereo_atoms: Option<[usize; 2]>,
    reverse_double_bond: bool,
    extra_right_substituent: bool,
) -> TopologyBlock {
    let mut elements = vec![Element::F, Element::C, Element::C, Element::CL];
    if extra_right_substituent {
        elements.push(Element::BR);
    }
    let atoms = elements
        .into_iter()
        .enumerate()
        .map(|(index, element)| Atom::from_spec(AtomId::new(index), AtomSpec::new(element)))
        .collect();
    let (double_begin, double_end) = if reverse_double_bond { (2, 1) } else { (1, 2) };
    let bonds = vec![
        Bond::from_spec(
            BondId::new(0),
            q81_stereo_bond_spec(0, 1, BondOrder::Single, BondStereo::None, None),
        ),
        Bond::from_spec(
            BondId::new(1),
            q81_stereo_bond_spec(
                double_begin,
                double_end,
                BondOrder::Double,
                stereo,
                stereo_atoms,
            ),
        ),
        Bond::from_spec(
            BondId::new(2),
            q81_stereo_bond_spec(2, 3, BondOrder::Single, BondStereo::None, None),
        ),
    ];
    let bonds = if extra_right_substituent {
        let mut bonds = bonds;
        bonds.push(Bond::from_spec(
            BondId::new(3),
            q81_stereo_bond_spec(2, 4, BondOrder::Single, BondStereo::None, None),
        ));
        bonds
    } else {
        bonds
    };
    TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
        .expect("fixed Q81 double-bond stereo target is structurally valid")
}

fn q82_stereo_group(kind: StereoGroupKind, atom_indices: &[usize]) -> StereoGroup {
    StereoGroup::new(
        kind,
        atom_indices.iter().copied().map(AtomId::new).collect(),
        Vec::new(),
    )
}

fn q82_stereo_spec(index: usize, element: Element, center_tags: [ChiralTag; 2]) -> AtomSpec {
    let spec = AtomSpec::new(element);
    match index {
        0 => spec.with_chiral_tag(center_tags[0]),
        4 => spec.with_chiral_tag(center_tags[1]),
        _ => spec,
    }
}

fn q82_stereo_query(center_tags: [ChiralTag; 2], stereo_groups: Vec<StereoGroup>) -> QueryGraph {
    let elements = [
        Element::C,
        Element::F,
        Element::CL,
        Element::BR,
        Element::C,
        Element::O,
        Element::N,
        Element::S,
    ];
    let atoms = elements
        .into_iter()
        .enumerate()
        .map(|(index, element)| {
            QueryAtom::new(
                AtomId::new(index),
                q82_stereo_spec(index, element, center_tags),
            )
        })
        .collect();
    let bonds = [(0, 1), (0, 2), (0, 3), (4, 5), (4, 6), (4, 7)]
        .into_iter()
        .enumerate()
        .map(|(index, (begin, end))| {
            QueryBond::new(
                BondId::new(index),
                BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
            )
        })
        .collect();
    QueryGraph::from_parts(
        atoms,
        bonds,
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        stereo_groups,
    )
    .expect("fixed Q82 enhanced-stereo query graph is valid")
}

fn q82_stereo_target(
    center_tags: [ChiralTag; 2],
    stereo_groups: Vec<StereoGroup>,
) -> TopologyBlock {
    let elements = [
        Element::C,
        Element::F,
        Element::CL,
        Element::BR,
        Element::C,
        Element::O,
        Element::N,
        Element::S,
    ];
    let atoms = elements
        .into_iter()
        .enumerate()
        .map(|(index, element)| {
            Atom::from_spec(
                AtomId::new(index),
                q82_stereo_spec(index, element, center_tags),
            )
        })
        .collect();
    let bonds = [(0, 1), (0, 2), (0, 3), (4, 5), (4, 6), (4, 7)]
        .into_iter()
        .enumerate()
        .map(|(index, (begin, end))| {
            Bond::from_spec(
                BondId::new(index),
                BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
            )
        })
        .collect();
    TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), stereo_groups)
        .expect("fixed Q82 enhanced-stereo target topology is valid")
}

fn mass_labeled_carbon_pair_query() -> QueryGraph {
    let mass_carbon = QueryAtom::from_parts(
        Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
        QueryNode::and(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
            QueryNode::predicate(AtomQueryPredicate::Mass(u16::MAX)),
        ]),
    );
    let carbon = QueryAtom::new(AtomId::new(1), AtomSpec::new(Element::C));
    let bond = QueryBond::new(
        BondId::new(0),
        BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
    );
    QueryGraph::from_parts(
        vec![mass_carbon, carbon],
        vec![bond],
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("fixed Q72 mass-labeled carbon pair query is valid")
}

fn empty_query() -> QueryGraph {
    QueryGraph::from_parts(
        Vec::new(),
        Vec::new(),
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("empty Q70 query graph is valid")
}

fn recursive_query_atom(
    atom_index: usize,
    inner_query: QueryGraph,
    serial_number: u32,
) -> QueryAtom {
    QueryAtom::from_parts(
        Atom::from_spec(AtomId::new(atom_index), AtomSpec::new(Element::C)),
        QueryNode::predicate(AtomQueryPredicate::RecursiveSmarts(
            RecursiveStructureQuery::from_query_graph(inner_query, serial_number),
        )),
    )
}

fn q86_query_graph(atoms: Vec<QueryAtom>, bonds: Vec<QueryBond>) -> QueryGraph {
    QueryGraph::from_parts(
        atoms,
        bonds,
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("fixed Q86 query graph is valid")
}

#[test]
fn q70_vf2_state_empty_query_and_target_boundaries() {
    let empty_target = topology(&[], &[]);
    let carbon_target = topology(&[Element::C], &[]);
    let compiled_empty = compile_query(&empty_query()).expect("empty query compiles");
    let compiled_carbon = compile_query(
        &parse_smarts("C", &SmartsParseParams::default()).expect("carbon query parses"),
    )
    .expect("carbon query compiles");

    assert!(
        compiled_empty
            .matches(&empty_target)
            .expect("empty query and target are valid")
            .is_empty()
    );
    assert!(
        compiled_empty
            .matches(&carbon_target)
            .expect("empty query against a nonempty target is valid")
            .is_empty()
    );
    assert!(
        compiled_carbon
            .matches(&empty_target)
            .expect("nonempty query against an empty target is valid")
            .is_empty()
    );
}

#[test]
fn q70_vf2_state_initializes_unequal_graphs_and_disconnected_components() {
    let query = parse_smarts("C.C-O", &SmartsParseParams::default())
        .expect("fixed disconnected Q70 query parses");
    let compiled = compile_query(&query).expect("disconnected query compiles");
    let target = topology(&[Element::C, Element::C, Element::O, Element::N], &[(1, 2)]);

    let matches = compiled
        .matches(&target)
        .expect("disconnected query matches a larger detached target");
    assert!(
        matches
            .iter()
            .any(|matched| { matched.atom_mapping == [0, 1, 2] && matched.bond_mapping == [0] })
    );
    assert!(matches.iter().all(|matched| {
        matched.atom_mapping.len() == compiled.num_atoms()
            && matched.bond_mapping.len() == compiled.num_bonds()
    }));
}

#[test]
fn q71_next_pair_enumerates_target_candidates_in_index_order() {
    let query = parse_smarts("C", &SmartsParseParams::default()).expect("carbon query parses");
    let compiled = compile_query(&query).expect("carbon query compiles");
    let target = topology(&[Element::C, Element::C, Element::C], &[]);

    let mappings: Vec<_> = compiled
        .matches(&target)
        .expect("carbon query checks all disconnected target candidates")
        .into_iter()
        .map(|matched| matched.atom_mapping)
        .collect();
    assert_eq!(mappings, [vec![0], vec![1], vec![2]]);
}

#[test]
fn q71_next_pair_walks_frontiers_and_exhausts_neighbor_candidates() {
    let query = parse_smarts("C(-O)-N", &SmartsParseParams::default())
        .expect("fixed Q71 branched query parses");
    let compiled = compile_query(&query).expect("branched query compiles");
    let target = topology(&[Element::C, Element::N, Element::O], &[(0, 1), (0, 2)]);

    let matches = compiled
        .matches(&target)
        .expect("frontier traversal checks target neighbors in adjacency order");
    assert_eq!(matches.len(), 1);
    assert_eq!(matches[0].atom_mapping, [0, 2, 1]);

    let exhausted_target = topology(&[Element::C, Element::C], &[(0, 1)]);
    let carbon_oxygen = compile_query(
        &parse_smarts("C-O", &SmartsParseParams::default()).expect("carbon-oxygen query parses"),
    )
    .expect("carbon-oxygen query compiles");
    assert!(
        carbon_oxygen
            .matches(&exhausted_target)
            .expect("frontier traversal exhausts incompatible target neighbors")
            .is_empty()
    );
}

#[test]
fn q71_next_pair_starts_new_query_components_and_exhausts_target_nodes() {
    let query = parse_smarts("C.O", &SmartsParseParams::default())
        .expect("fixed Q71 disconnected query parses");
    let compiled = compile_query(&query).expect("disconnected query compiles");
    let target = topology(&[Element::C, Element::O, Element::C], &[]);

    let mappings: BTreeSet<_> = compiled
        .matches(&target)
        .expect("all disconnected component candidates are enumerated")
        .into_iter()
        .map(|matched| matched.atom_mapping)
        .collect();
    assert_eq!(mappings, BTreeSet::from([vec![0, 1], vec![2, 1]]));
}

#[test]
fn q72_degree_short_circuit_and_missing_isotope_mass_follow_source() {
    let query = mass_labeled_carbon_pair_query();
    let compiled = compile_query(&query).expect("mass-labeled carbon pair compiles");
    let disconnected = topology_with_isotope(&[Element::C, Element::C], &[], u16::MAX);
    let bonded = topology_with_isotope(&[Element::C, Element::C], &[(0, 1)], u16::MAX);

    assert!(
        compiled
            .matches(&disconnected)
            .expect("degree rejection does not produce a matcher error")
            .is_empty()
    );

    let matches = compiled
        .matches(&bonded)
        .expect("missing isotope mass uses the pinned isotope-number fallback");
    assert_eq!(matches.len(), 1);
    assert_eq!(matches[0].bond_mapping, [0]);
    assert_eq!(matches[0].atom_mapping.len(), 2);
}

#[test]
fn q72_mapped_edge_feasibility_rejects_triangle_in_square() {
    let triangle = compile_query(
        &parse_smarts("C1CC1", &SmartsParseParams::default())
            .expect("fixed Q72 triangle query parses"),
    )
    .expect("triangle query compiles");
    let square = topology(
        &[Element::C, Element::C, Element::C, Element::C],
        &[(0, 1), (1, 2), (2, 3), (3, 0)],
    );

    assert!(
        triangle
            .matches(&square)
            .expect("mapped edge feasibility evaluates a valid triangle-to-square case")
            .is_empty()
    );
}

#[test]
fn q72_atom_and_bond_labels_reject_equal_degree_candidates() {
    let carbon = compile_query(
        &parse_smarts("C", &SmartsParseParams::default()).expect("carbon query parses"),
    )
    .expect("carbon query compiles");
    assert!(
        carbon
            .matches(&topology(&[Element::N], &[]))
            .expect("atom-label mismatch is an ordinary non-match")
            .is_empty()
    );

    let carbon_oxygen_double = compile_query(
        &parse_smarts("C=O", &SmartsParseParams::default())
            .expect("fixed Q72 double-bond query parses"),
    )
    .expect("double-bond query compiles");
    assert!(
        carbon_oxygen_double
            .matches(&topology(&[Element::C, Element::O], &[(0, 1)]))
            .expect("bond-label mismatch is an ordinary non-match")
            .is_empty()
    );
}

#[test]
fn q73_add_pair_tracks_mapping_and_frontier_depth_across_branches() {
    let query = parse_smarts("C(-O)-N", &SmartsParseParams::default())
        .expect("fixed Q73 branched query parses");
    let compiled = compile_query(&query).expect("branched query compiles");
    let target = topology(
        &[Element::C, Element::O, Element::N, Element::O],
        &[(0, 1), (0, 2), (0, 3)],
    );

    let observed: BTreeSet<_> = compiled
        .matches(&target)
        .expect("both query-to-target branch mappings are enumerated")
        .into_iter()
        .map(|matched| (matched.atom_mapping, matched.bond_mapping))
        .collect();
    assert_eq!(
        observed,
        BTreeSet::from([(vec![0, 1, 2], vec![0, 1]), (vec![0, 3, 2], vec![2, 1]),])
    );
}

#[test]
fn q73_add_pair_advances_depth_across_disconnected_components() {
    let query = parse_smarts("C.O", &SmartsParseParams::default())
        .expect("fixed Q73 disconnected query parses");
    let compiled = compile_query(&query).expect("disconnected query compiles");
    let target = topology(&[Element::C, Element::C, Element::O, Element::O], &[]);

    let observed: BTreeSet<_> = compiled
        .matches(&target)
        .expect("component candidates advance the same match state")
        .into_iter()
        .map(|matched| matched.atom_mapping)
        .collect();
    assert_eq!(
        observed,
        BTreeSet::from([vec![0, 2], vec![0, 3], vec![1, 2], vec![1, 3]])
    );
}

#[test]
fn q74_backtrack_restores_state_after_rejected_sibling_branch() {
    let query = parse_smarts("C(-O-C-F)(-N)-S", &SmartsParseParams::default())
        .expect("fixed Q74 branched query parses");
    let compiled = compile_query(&query).expect("branched query compiles");
    let target = topology(
        &[
            Element::C,
            Element::N,
            Element::O,
            Element::C,
            Element::O,
            Element::C,
            Element::F,
            Element::S,
        ],
        &[(0, 1), (0, 2), (2, 3), (0, 7), (0, 4), (4, 5), (5, 6)],
    );

    let matches = compiled
        .matches(&target)
        .expect("rejected sibling branch restores the match state");
    assert_eq!(matches.len(), 1);
    assert_eq!(matches[0].atom_mapping, [0, 4, 5, 6, 1, 7]);
    assert_eq!(matches[0].bond_mapping, [4, 5, 6, 0, 3]);
}

#[test]
fn q75_first_match_uses_first_disconnected_candidate_and_handles_empty_graphs() {
    let carbon_query =
        parse_smarts("C", &SmartsParseParams::default()).expect("carbon query parses");
    let disconnected_target = topology(&[Element::C, Element::C, Element::C], &[]);
    let coordinates = CoordinateBlock::default();
    let target = SearchTarget::new(
        &disconnected_target,
        &coordinates,
        &disconnected_target.stereo_groups,
        None,
        None,
    );

    let first = get_substruct_match(&target, &carbon_query)
        .expect("single-match search returns the first disconnected target candidate");
    assert_eq!(first.atom_mapping, [0]);
    assert!(first.bond_mapping.is_empty());

    let empty_target = topology(&[], &[]);
    let empty_coordinates = CoordinateBlock::default();
    let empty_target_view = SearchTarget::new(
        &empty_target,
        &empty_coordinates,
        &empty_target.stereo_groups,
        None,
        None,
    );
    assert!(get_substruct_match(&empty_target_view, &carbon_query).is_none());

    let one_atom_target = topology(&[Element::C], &[]);
    let one_atom_coordinates = CoordinateBlock::default();
    let one_atom_view = SearchTarget::new(
        &one_atom_target,
        &one_atom_coordinates,
        &one_atom_target.stereo_groups,
        None,
        None,
    );
    assert!(get_substruct_match(&one_atom_view, &empty_query()).is_none());

    let larger_query =
        parse_smarts("C-C", &SmartsParseParams::default()).expect("two-atom query parses");
    assert!(get_substruct_match(&one_atom_view, &larger_query).is_none());
}

#[test]
fn q76_match_all_preserves_source_order_limits_and_exact_completed_mappings() {
    let query =
        parse_smarts("C-C-O", &SmartsParseParams::default()).expect("fixed Q76 chain query parses");
    let target_topology = topology(
        &[
            Element::C,
            Element::C,
            Element::C,
            Element::O,
            Element::C,
            Element::C,
            Element::O,
        ],
        &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
    );
    let coordinates = CoordinateBlock::default();
    let target = SearchTarget::new(
        &target_topology,
        &coordinates,
        &target_topology.stereo_groups,
        None,
        None,
    );
    let expected = vec![
        (vec![1, 2, 3], vec![1, 2]),
        (vec![4, 5, 6], vec![4, 5]),
        (vec![5, 4, 3], vec![4, 3]),
    ];

    let mut params = SubstructMatchParams {
        max_matches: 0,
        ..SubstructMatchParams::default()
    };
    let unlimited = get_substruct_matches_with_params(&target, &query, &params)
        .into_iter()
        .map(|matched| (matched.atom_mapping, matched.bond_mapping))
        .collect::<Vec<_>>();
    assert_eq!(unlimited, expected);

    params.max_matches = 1;
    let limited_to_one = get_substruct_matches_with_params(&target, &query, &params)
        .into_iter()
        .map(|matched| (matched.atom_mapping, matched.bond_mapping))
        .collect::<Vec<_>>();
    assert_eq!(limited_to_one, expected[..1]);

    params.max_matches = 2;
    let limited_to_two = get_substruct_matches_with_params(&target, &query, &params)
        .into_iter()
        .map(|matched| (matched.atom_mapping, matched.bond_mapping))
        .collect::<Vec<_>>();
    assert_eq!(limited_to_two, expected[..2]);

    let compiled = compile_query(&query).expect("query compiles");
    let compiled_matches = compiled
        .matches(&target_topology)
        .expect("compiled query follows the same source result order");
    assert_eq!(
        compiled_matches
            .into_iter()
            .map(|matched| (matched.atom_mapping, matched.bond_mapping))
            .collect::<Vec<_>>(),
        expected
    );
}

#[test]
fn q84_vf2_results_keep_query_to_target_direction_source_order_and_limit() {
    let query = parse_smarts("C-C", &SmartsParseParams::default())
        .expect("fixed Q84 two-carbon query parses");
    let target_topology = topology(&[Element::C, Element::C, Element::C], &[(0, 1), (1, 2)]);
    let coordinates = CoordinateBlock::default();
    let target = SearchTarget::new(
        &target_topology,
        &coordinates,
        &target_topology.stereo_groups,
        None,
        None,
    );
    let expected = vec![
        (vec![0, 1], vec![0]),
        (vec![1, 0], vec![0]),
        (vec![1, 2], vec![1]),
        (vec![2, 1], vec![1]),
    ];
    let mut params = SubstructMatchParams {
        uniquify: false,
        max_matches: 0,
        ..SubstructMatchParams::default()
    };

    let unlimited = get_substruct_matches_with_params(&target, &query, &params)
        .into_iter()
        .map(|matched| (matched.atom_mapping, matched.bond_mapping))
        .collect::<Vec<_>>();
    assert_eq!(unlimited, expected);

    params.max_matches = 2;
    let limited = get_substruct_matches_with_params(&target, &query, &params)
        .into_iter()
        .map(|matched| (matched.atom_mapping, matched.bond_mapping))
        .collect::<Vec<_>>();
    assert_eq!(limited, expected[..2]);
}

#[test]
fn q85_uniquify_keeps_first_mapping_per_target_atom_set_in_vf2_order() {
    let query = parse_smarts("C-C", &SmartsParseParams::default())
        .expect("fixed Q85 two-carbon query parses");
    let target_topology = topology(&[Element::C, Element::C, Element::C], &[(0, 1), (1, 2)]);
    let coordinates = CoordinateBlock::default();
    let target = SearchTarget::new(
        &target_topology,
        &coordinates,
        &target_topology.stereo_groups,
        None,
        None,
    );
    let expected = vec![(vec![0, 1], vec![0]), (vec![1, 2], vec![1])];
    let mut params = SubstructMatchParams {
        uniquify: true,
        max_matches: 0,
        ..SubstructMatchParams::default()
    };

    let unique = get_substruct_matches_with_params(&target, &query, &params)
        .into_iter()
        .map(|matched| (matched.atom_mapping, matched.bond_mapping))
        .collect::<Vec<_>>();
    assert_eq!(unique, expected);

    params.max_matches = 1;
    let limited = get_substruct_matches_with_params(&target, &query, &params)
        .into_iter()
        .map(|matched| (matched.atom_mapping, matched.bond_mapping))
        .collect::<Vec<_>>();
    assert_eq!(limited, expected[..1]);
}

#[test]
fn q77_recursive_cache_uses_root_scope_and_reuses_serial_labels() {
    let mut rooted_inner =
        parse_smarts("C-O", &SmartsParseParams::default()).expect("recursive query parses");
    rooted_inner.set_prop("_queryRootAtom", "1");
    let rooted_query = QueryGraph::from_parts(
        vec![recursive_query_atom(0, rooted_inner, 0)],
        Vec::new(),
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("single rooted recursive atom query is valid");
    let rooted_target_topology = topology(&[Element::C, Element::O, Element::N], &[(0, 1)]);
    let rooted_coordinates = CoordinateBlock::default();
    let rooted_target = SearchTarget::new(
        &rooted_target_topology,
        &rooted_coordinates,
        &rooted_target_topology.stereo_groups,
        None,
        None,
    );
    let rooted_matches = get_substruct_matches_with_params(
        &rooted_target,
        &rooted_query,
        &SubstructMatchParams::default(),
    );
    assert_eq!(rooted_matches.len(), 1);
    assert_eq!(rooted_matches[0].atom_mapping, [1]);

    let carbon_inner =
        parse_smarts("C", &SmartsParseParams::default()).expect("carbon recursive query parses");
    let nitrogen_inner =
        parse_smarts("N", &SmartsParseParams::default()).expect("nitrogen recursive query parses");
    let shared_serial_query = QueryGraph::from_parts(
        vec![
            recursive_query_atom(0, carbon_inner, 77),
            recursive_query_atom(1, nitrogen_inner, 77),
        ],
        vec![QueryBond::new(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        )],
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("shared-serial recursive query graph is valid");
    let serial_target_topology = topology(&[Element::C, Element::C], &[(0, 1)]);
    let serial_coordinates = CoordinateBlock::default();
    let serial_target = SearchTarget::new(
        &serial_target_topology,
        &serial_coordinates,
        &serial_target_topology.stereo_groups,
        None,
        None,
    );
    let serial_matches = get_substruct_matches_with_params(
        &serial_target,
        &shared_serial_query,
        &SubstructMatchParams::default(),
    );
    assert_eq!(serial_matches.len(), 1);
    assert_eq!(serial_matches[0].atom_mapping, [0, 1]);
    assert_eq!(serial_matches[0].bond_mapping, [0]);
}

#[test]
fn q77_recursive_cache_is_rebuilt_for_each_target_call() {
    let inner_query =
        parse_smarts("C", &SmartsParseParams::default()).expect("carbon recursive query parses");
    let query = QueryGraph::from_parts(
        vec![recursive_query_atom(0, inner_query, 0)],
        Vec::new(),
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("recursive atom query is valid");
    let carbon_topology = topology(&[Element::C], &[]);
    let nitrogen_topology = topology(&[Element::N], &[]);
    let coordinates = CoordinateBlock::default();
    let carbon_target = SearchTarget::new(
        &carbon_topology,
        &coordinates,
        &carbon_topology.stereo_groups,
        None,
        None,
    );
    let nitrogen_target = SearchTarget::new(
        &nitrogen_topology,
        &coordinates,
        &nitrogen_topology.stereo_groups,
        None,
        None,
    );
    let params = SubstructMatchParams::default();

    let carbon_matches = get_substruct_matches_with_params(&carbon_target, &query, &params);
    assert_eq!(carbon_matches.len(), 1);
    assert_eq!(carbon_matches[0].atom_mapping, [0]);
    assert!(get_substruct_matches_with_params(&nitrogen_target, &query, &params).is_empty());
}

#[test]
fn q78_recursive_matching_walks_nested_query_children() {
    let deepest_query =
        parse_smarts("O", &SmartsParseParams::default()).expect("oxygen query parses");
    let nested_query = QueryGraph::from_parts(
        vec![recursive_query_atom(0, deepest_query, 0)],
        Vec::new(),
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("nested recursive query graph is valid");
    let query = QueryGraph::from_parts(
        vec![recursive_query_atom(0, nested_query, 0)],
        Vec::new(),
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("outer recursive query graph is valid");
    let target_topology = topology(&[Element::C, Element::O, Element::N], &[(0, 1), (1, 2)]);
    let coordinates = CoordinateBlock::default();
    let target = SearchTarget::new(
        &target_topology,
        &coordinates,
        &target_topology.stereo_groups,
        None,
        None,
    );

    let matches =
        get_substruct_matches_with_params(&target, &query, &SubstructMatchParams::default());

    assert_eq!(
        matches
            .into_iter()
            .map(|matched| matched.atom_mapping)
            .collect::<Vec<_>>(),
        vec![vec![1]],
    );
}

#[test]
fn q78_recursive_match_limit_and_root_projection_follow_source() {
    let mut rooted_inner =
        parse_smarts("C-C", &SmartsParseParams::default()).expect("carbon pair query parses");
    rooted_inner.set_prop("_queryRootAtom", "1");
    let query = QueryGraph::from_parts(
        vec![recursive_query_atom(0, rooted_inner, 0)],
        Vec::new(),
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("rooted recursive query graph is valid");
    let target_topology = topology(&[Element::C, Element::C, Element::C], &[(0, 1), (1, 2)]);
    let coordinates = CoordinateBlock::default();
    let target = SearchTarget::new(
        &target_topology,
        &coordinates,
        &target_topology.stereo_groups,
        None,
        None,
    );
    let one_recursive_match = SubstructMatchParams {
        max_matches: 1,
        max_recursive_matches: 1,
        ..SubstructMatchParams::default()
    };
    let two_recursive_matches = SubstructMatchParams {
        max_matches: 1,
        max_recursive_matches: 2,
        ..SubstructMatchParams::default()
    };

    let one = get_substruct_matches_with_params(&target, &query, &one_recursive_match);
    let two = get_substruct_matches_with_params(&target, &query, &two_recursive_matches);
    assert_eq!(one.len(), 1);
    assert_eq!(one[0].atom_mapping, [1]);
    assert_eq!(two.len(), 1);
    assert_eq!(two[0].atom_mapping, [0]);

    let mut negative_root_inner =
        parse_smarts("C-C", &SmartsParseParams::default()).expect("carbon pair query parses");
    negative_root_inner.set_prop("_queryRootAtom", "-1");
    let negative_root_query = QueryGraph::from_parts(
        vec![recursive_query_atom(0, negative_root_inner, 0)],
        Vec::new(),
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("negative-root recursive query graph is valid");
    assert!(
        get_substruct_matches_with_params(&target, &negative_root_query, &one_recursive_match)
            .is_empty()
    );
}

#[test]
fn q79_atom_and_bond_callback_order_follows_source_short_circuits() {
    let atom_target_topology = topology(&[Element::C], &[]);
    let coordinates = CoordinateBlock::default();
    let atom_target = SearchTarget::new(
        &atom_target_topology,
        &coordinates,
        &atom_target_topology.stereo_groups,
        None,
        None,
    );
    let oxygen_query =
        parse_smarts("O", &SmartsParseParams::default()).expect("oxygen query parses");

    let post_atom_calls = Arc::new(AtomicUsize::new(0));
    let post_atom_calls_in_check = Arc::clone(&post_atom_calls);
    let post_atom_params = SubstructMatchParams {
        extra_atom_check: Some(Arc::new(move |_, _, _, _| {
            post_atom_calls_in_check.fetch_add(1, Ordering::SeqCst);
            true
        })),
        ..SubstructMatchParams::default()
    };
    assert!(
        get_substruct_matches_with_params(&atom_target, &oxygen_query, &post_atom_params)
            .is_empty()
    );
    assert_eq!(post_atom_calls.load(Ordering::SeqCst), 0);

    let override_atom_calls = Arc::new(AtomicUsize::new(0));
    let override_atom_calls_in_check = Arc::clone(&override_atom_calls);
    let override_atom_params = SubstructMatchParams {
        extra_atom_check: Some(Arc::new(move |_, _, _, _| {
            override_atom_calls_in_check.fetch_add(1, Ordering::SeqCst);
            true
        })),
        extra_atom_check_overrides_default_check: true,
        ..SubstructMatchParams::default()
    };
    let override_atom_matches =
        get_substruct_matches_with_params(&atom_target, &oxygen_query, &override_atom_params);
    assert_eq!(override_atom_matches.len(), 1);
    assert_eq!(override_atom_matches[0].atom_mapping, [0]);
    assert_eq!(override_atom_calls.load(Ordering::SeqCst), 1);

    let chiral_query = QueryGraph::from_parts(
        vec![QueryAtom::new(
            AtomId::new(0),
            AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCw),
        )],
        Vec::new(),
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("chiral query graph is valid");
    let chiral_precheck_calls = Arc::new(AtomicUsize::new(0));
    let chiral_precheck_calls_in_check = Arc::clone(&chiral_precheck_calls);
    let chiral_precheck_params = SubstructMatchParams {
        use_chirality: true,
        extra_atom_check: Some(Arc::new(move |_, _, _, _| {
            chiral_precheck_calls_in_check.fetch_add(1, Ordering::SeqCst);
            true
        })),
        extra_atom_check_overrides_default_check: true,
        ..SubstructMatchParams::default()
    };
    assert!(
        get_substruct_matches_with_params(&atom_target, &chiral_query, &chiral_precheck_params)
            .is_empty()
    );
    assert_eq!(chiral_precheck_calls.load(Ordering::SeqCst), 0);

    let bond_target_topology = topology(&[Element::C, Element::O], &[(0, 1)]);
    let bond_target = SearchTarget::new(
        &bond_target_topology,
        &coordinates,
        &bond_target_topology.stereo_groups,
        None,
        None,
    );
    let double_bond_query =
        parse_smarts("C=O", &SmartsParseParams::default()).expect("double-bond query parses");

    let post_bond_calls = Arc::new(AtomicUsize::new(0));
    let post_bond_calls_in_check = Arc::clone(&post_bond_calls);
    let post_bond_params = SubstructMatchParams {
        extra_bond_check: Some(Arc::new(move |_, _| {
            post_bond_calls_in_check.fetch_add(1, Ordering::SeqCst);
            true
        })),
        ..SubstructMatchParams::default()
    };
    assert!(
        get_substruct_matches_with_params(&bond_target, &double_bond_query, &post_bond_params)
            .is_empty()
    );
    assert_eq!(post_bond_calls.load(Ordering::SeqCst), 0);

    let override_bond_calls = Arc::new(AtomicUsize::new(0));
    let override_bond_calls_in_check = Arc::clone(&override_bond_calls);
    let override_bond_params = SubstructMatchParams {
        extra_bond_check: Some(Arc::new(move |_, _| {
            override_bond_calls_in_check.fetch_add(1, Ordering::SeqCst);
            true
        })),
        extra_bond_check_overrides_default_check: true,
        ..SubstructMatchParams::default()
    };
    let override_bond_matches =
        get_substruct_matches_with_params(&bond_target, &double_bond_query, &override_bond_params);
    assert_eq!(override_bond_matches.len(), 1);
    assert_eq!(override_bond_matches[0].atom_mapping, [0, 1]);
    assert_eq!(override_bond_calls.load(Ordering::SeqCst), 1);
}

#[test]
fn q79_unsupported_compatibility_keeps_typed_error_before_callback() {
    let query = QueryGraph::from_parts(
        vec![QueryAtom::from_parts(
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            QueryNode::predicate(AtomQueryPredicate::UnsupportedFeature(
                "Q79 fixed unsupported atom leaf",
            )),
        )],
        Vec::new(),
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("query graph with explicit unsupported leaf is structurally valid");
    let target_topology = topology(&[Element::C], &[]);
    let coordinates = CoordinateBlock::default();
    let target = SearchTarget::new(
        &target_topology,
        &coordinates,
        &target_topology.stereo_groups,
        None,
        None,
    );
    let callback_calls = Arc::new(AtomicUsize::new(0));
    let callback_calls_in_check = Arc::clone(&callback_calls);
    let params = SubstructMatchParams {
        extra_atom_check: Some(Arc::new(move |_, _, _, _| {
            callback_calls_in_check.fetch_add(1, Ordering::SeqCst);
            true
        })),
        extra_atom_check_overrides_default_check: true,
        ..SubstructMatchParams::default()
    };

    let error = try_get_substruct_matches_with_params(&target, &query, &params)
        .expect_err("unsupported query leaves must retain their structured error");
    assert_eq!(
        error,
        SubstructMatchError::Unsupported {
            branch: "Q79 fixed unsupported atom leaf",
            rdkit_function: "QueryAtom::Match",
        }
    );
    assert_eq!(callback_calls.load(Ordering::SeqCst), 0);
}

#[test]
fn q80_tetrahedral_mapping_permutation_and_implicit_hydrogen_follow_source() {
    let query_tag = ChiralTag::TetrahedralCcw;
    let query = q80_tetrahedral_query(query_tag);
    assert_eq!(query.adjacency()[0].len(), 3);

    let query_bond_order = [(0, 1), (0, 2), (0, 3)];
    let coordinates = CoordinateBlock::default();
    let same_order_topology = q80_tetrahedral_target(query_tag, &query_bond_order, false);
    let same_order_target = SearchTarget::new(
        &same_order_topology,
        &coordinates,
        &same_order_topology.stereo_groups,
        None,
        None,
    );
    let use_chirality = SubstructMatchParams {
        use_chirality: true,
        ..SubstructMatchParams::default()
    };
    assert_eq!(
        get_substruct_matches_with_params(&same_order_target, &query, &use_chirality).len(),
        1,
        "the source bond order and one implicit H preserve the specified tetrahedral match"
    );

    let reversed_bond_order = [(0, 3), (0, 2), (0, 1)];
    let same_tag_reversed_topology = q80_tetrahedral_target(query_tag, &reversed_bond_order, false);
    let same_tag_reversed_target = SearchTarget::new(
        &same_tag_reversed_topology,
        &coordinates,
        &same_tag_reversed_topology.stereo_groups,
        None,
        None,
    );
    assert!(
        get_substruct_matches_with_params(&same_tag_reversed_target, &query, &use_chirality)
            .is_empty(),
        "an odd neighbor-order permutation requires the opposite source chiral label"
    );

    let opposite_tag = match query_tag {
        ChiralTag::TetrahedralCw => ChiralTag::TetrahedralCcw,
        ChiralTag::TetrahedralCcw => ChiralTag::TetrahedralCw,
        _ => unreachable!("the parsed Q80 query tag was checked above"),
    };
    let opposite_tag_topology = q80_tetrahedral_target(opposite_tag, &reversed_bond_order, false);
    let opposite_tag_target = SearchTarget::new(
        &opposite_tag_topology,
        &coordinates,
        &opposite_tag_topology.stereo_groups,
        None,
        None,
    );
    assert_eq!(
        get_substruct_matches_with_params(&opposite_tag_target, &query, &use_chirality).len(),
        1,
        "the opposite source label compensates for the odd neighbor-order permutation"
    );

    let no_implicit_h_topology = q80_tetrahedral_target(query_tag, &query_bond_order, true);
    let no_implicit_h_target = SearchTarget::new(
        &no_implicit_h_topology,
        &coordinates,
        &no_implicit_h_topology.stereo_groups,
        None,
        None,
    );
    assert!(
        get_substruct_matches_with_params(&no_implicit_h_target, &query, &use_chirality).is_empty(),
        "the source AtomHCount leaf still requires one target hydrogen"
    );
}

#[test]
fn q80_stereo_options_match_source_unspecified_short_circuit() {
    let query = q80_tetrahedral_query(ChiralTag::TetrahedralCcw);
    let edge_order = [(0, 1), (0, 2), (0, 3)];
    let topology = q80_tetrahedral_target(ChiralTag::Unspecified, &edge_order, false);
    let coordinates = CoordinateBlock::default();
    let target = SearchTarget::new(&topology, &coordinates, &topology.stereo_groups, None, None);

    let strict_stereo = SubstructMatchParams {
        use_chirality: true,
        ..SubstructMatchParams::default()
    };
    assert!(
        get_substruct_matches_with_params(&target, &query, &strict_stereo).is_empty(),
        "specified query stereo rejects an unspecified target by default"
    );

    let allow_unspecified = SubstructMatchParams {
        use_chirality: true,
        specified_stereo_query_matches_unspecified: true,
        ..SubstructMatchParams::default()
    };
    assert_eq!(
        get_substruct_matches_with_params(&target, &query, &allow_unspecified).len(),
        1,
        "the source option skips tetrahedral comparison for an unspecified target"
    );

    let ignore_stereo = SubstructMatchParams {
        use_chirality: false,
        ..SubstructMatchParams::default()
    };
    assert_eq!(
        get_substruct_matches_with_params(&target, &query, &ignore_stereo).len(),
        1,
        "useChirality=false bypasses only the stereo final check"
    );

    let no_h_topology = q80_tetrahedral_target(ChiralTag::Unspecified, &edge_order, true);
    let no_h_target = SearchTarget::new(
        &no_h_topology,
        &coordinates,
        &no_h_topology.stereo_groups,
        None,
        None,
    );
    assert!(
        get_substruct_matches_with_params(&no_h_target, &query, &allow_unspecified).is_empty(),
        "allowing unspecified tetrahedral stereo does not bypass the independent AtomHCount leaf"
    );
}

#[test]
fn q81_double_bond_stereo_uses_mapped_neighbors_and_bond_orientation() {
    let query = q81_stereo_query(BondStereo::E, true);
    let coordinates = CoordinateBlock::default();
    let params = SubstructMatchParams {
        use_chirality: true,
        ..SubstructMatchParams::default()
    };
    let match_count = |topology: &TopologyBlock| {
        let target = SearchTarget::new(topology, &coordinates, &topology.stereo_groups, None, None);
        get_substruct_matches_with_params(&target, &query, &params).len()
    };

    let same_orientation = q81_stereo_target(BondStereo::E, Some([0, 3]), false, false);
    assert_eq!(match_count(&same_orientation), 1);

    let normalized_e_label = q81_stereo_target(BondStereo::Trans, Some([0, 3]), false, false);
    assert_eq!(
        match_count(&normalized_e_label),
        1,
        "source E and TRANS labels translate to the same cis/trans label"
    );

    let different_stereo_two_refs = q81_stereo_target(BondStereo::Z, Some([0, 3]), false, false);
    assert_eq!(
        match_count(&different_stereo_two_refs),
        0,
        "different stereo labels reject when both query references map"
    );

    let same_stereo_one_ref = q81_stereo_target(BondStereo::E, Some([0, 4]), false, true);
    assert_eq!(
        match_count(&same_stereo_one_ref),
        0,
        "equal stereo labels reject when exactly one reference maps"
    );

    let opposite_stereo_one_ref = q81_stereo_target(BondStereo::Z, Some([0, 4]), false, true);
    assert_eq!(
        match_count(&opposite_stereo_one_ref),
        1,
        "opposite stereo labels match when exactly one reference maps"
    );

    let reversed_endpoints_one_ref = q81_stereo_target(BondStereo::Z, Some([4, 0]), true, true);
    assert_eq!(
        match_count(&reversed_endpoints_one_ref),
        1,
        "the source reversed-endpoint branch swaps target stereo-reference slots"
    );
}

#[test]
fn q81_unspecified_bond_stereo_option_matches_source() {
    let query = q81_stereo_query(BondStereo::E, true);
    let coordinates = CoordinateBlock::default();
    for target_stereo in [BondStereo::None, BondStereo::Any] {
        let topology = q81_stereo_target(target_stereo, None, false, false);
        let target =
            SearchTarget::new(&topology, &coordinates, &topology.stereo_groups, None, None);
        let strict = SubstructMatchParams {
            use_chirality: true,
            ..SubstructMatchParams::default()
        };
        assert!(
            get_substruct_matches_with_params(&target, &query, &strict).is_empty(),
            "specified bond stereo rejects target {target_stereo:?} by default"
        );

        let allow_unspecified = SubstructMatchParams {
            use_chirality: true,
            specified_stereo_query_matches_unspecified: true,
            ..SubstructMatchParams::default()
        };
        assert_eq!(
            get_substruct_matches_with_params(&target, &query, &allow_unspecified).len(),
            1,
            "the source option permits target {target_stereo:?} and final check skips absent refs"
        );

        let ignore_chirality = SubstructMatchParams::default();
        assert_eq!(
            get_substruct_matches_with_params(&target, &query, &ignore_chirality).len(),
            1,
            "useChirality=false bypasses the bond stereo precheck and final check"
        );
    }

    let unspecified_query = q81_stereo_query(BondStereo::None, false);
    let specified_target = q81_stereo_target(BondStereo::E, Some([0, 3]), false, false);
    let target = SearchTarget::new(
        &specified_target,
        &coordinates,
        &specified_target.stereo_groups,
        None,
        None,
    );
    let use_chirality = SubstructMatchParams {
        use_chirality: true,
        ..SubstructMatchParams::default()
    };
    assert_eq!(
        get_substruct_matches_with_params(&target, &unspecified_query, &use_chirality).len(),
        1,
        "an unspecified query bond does not enter the source stereo final-check loop"
    );
}

#[test]
fn q82_query_enhanced_stereo_group_kinds_follow_source() {
    let center_tags = [ChiralTag::TetrahedralCw, ChiralTag::TetrahedralCw];
    let coordinates = CoordinateBlock::default();
    let params = SubstructMatchParams {
        use_chirality: true,
        use_enhanced_stereo: true,
        ..SubstructMatchParams::default()
    };
    let match_count = |query_kind, target_kind| {
        let query = q82_stereo_query(center_tags, vec![q82_stereo_group(query_kind, &[0, 4])]);
        let topology = q82_stereo_target(center_tags, vec![q82_stereo_group(target_kind, &[0, 4])]);
        let target =
            SearchTarget::new(&topology, &coordinates, &topology.stereo_groups, None, None);
        get_substruct_matches_with_params(&target, &query, &params).len()
    };

    assert_eq!(
        match_count(StereoGroupKind::Or, StereoGroupKind::Or),
        1,
        "source OR query groups match target OR groups"
    );
    assert_eq!(
        match_count(StereoGroupKind::Or, StereoGroupKind::And),
        1,
        "source OR query groups also match target AND groups"
    );
    assert_eq!(
        match_count(StereoGroupKind::And, StereoGroupKind::And),
        1,
        "source AND query groups match target AND groups"
    );
    assert_eq!(
        match_count(StereoGroupKind::And, StereoGroupKind::Or),
        0,
        "source AND query groups reject target OR groups"
    );
    assert_eq!(
        match_count(StereoGroupKind::Or, StereoGroupKind::Absolute),
        0,
        "source OR query groups reject target absolute stereochemistry"
    );
    assert_eq!(
        match_count(StereoGroupKind::Absolute, StereoGroupKind::Or),
        1,
        "source absolute query groups are skipped by enhanced-group matching"
    );

    let query = q82_stereo_query(
        center_tags,
        vec![q82_stereo_group(StereoGroupKind::And, &[0, 4])],
    );
    let topology = q82_stereo_target(
        center_tags,
        vec![q82_stereo_group(StereoGroupKind::Or, &[0, 4])],
    );
    let target = SearchTarget::new(&topology, &coordinates, &topology.stereo_groups, None, None);
    let enhanced_disabled = SubstructMatchParams {
        use_chirality: true,
        use_enhanced_stereo: false,
        ..SubstructMatchParams::default()
    };
    assert_eq!(
        get_substruct_matches_with_params(&target, &query, &enhanced_disabled).len(),
        1,
        "source useEnhancedStereo=false bypasses group-kind matching"
    );
}

#[test]
fn q82_target_group_members_keep_source_parity_and_query_group_identity() {
    let query_tags = [ChiralTag::TetrahedralCw, ChiralTag::TetrahedralCw];
    let coordinates = CoordinateBlock::default();
    let params = SubstructMatchParams {
        use_chirality: true,
        use_enhanced_stereo: true,
        ..SubstructMatchParams::default()
    };

    let query = q82_stereo_query(
        query_tags,
        vec![q82_stereo_group(StereoGroupKind::Absolute, &[0, 4])],
    );
    let mixed_target = q82_stereo_target(
        [ChiralTag::TetrahedralCcw, ChiralTag::TetrahedralCw],
        vec![q82_stereo_group(StereoGroupKind::Or, &[0, 4])],
    );
    let mixed_target_view = SearchTarget::new(
        &mixed_target,
        &coordinates,
        &mixed_target.stereo_groups,
        None,
        None,
    );
    assert!(
        get_substruct_matches_with_params(&mixed_target_view, &query, &params).is_empty(),
        "source rejects opposite chiral match states within one target stereo group"
    );

    let split_query = q82_stereo_query(
        query_tags,
        vec![
            q82_stereo_group(StereoGroupKind::Or, &[0]),
            q82_stereo_group(StereoGroupKind::Or, &[4]),
        ],
    );
    let same_target = q82_stereo_target(
        query_tags,
        vec![q82_stereo_group(StereoGroupKind::Or, &[0, 4])],
    );
    let same_target_view = SearchTarget::new(
        &same_target,
        &coordinates,
        &same_target.stereo_groups,
        None,
        None,
    );
    assert!(
        get_substruct_matches_with_params(&same_target_view, &split_query, &params).is_empty(),
        "source rejects multiple query groups covering one target stereo group"
    );

    let ignore_chirality = SubstructMatchParams {
        use_chirality: false,
        use_enhanced_stereo: true,
        ..SubstructMatchParams::default()
    };
    assert_eq!(
        get_substruct_matches_with_params(&mixed_target_view, &query, &ignore_chirality).len(),
        1,
        "source useChirality=false returns before enhanced stereo final checks"
    );
}

#[test]
fn q83_extra_final_check_runs_before_stereo_rejection() {
    let edge_order = [(0, 1), (0, 2), (0, 3)];
    let query = q80_tetrahedral_query(ChiralTag::TetrahedralCw);
    let topology = q80_tetrahedral_target(ChiralTag::TetrahedralCcw, &edge_order, false);
    let coordinates = CoordinateBlock::default();
    let target = SearchTarget::new(&topology, &coordinates, &topology.stereo_groups, None, None);
    let final_check_calls = Arc::new(AtomicUsize::new(0));
    let calls_in_final_check = Arc::clone(&final_check_calls);
    let params = SubstructMatchParams {
        use_chirality: true,
        extra_final_check: Some(Arc::new(move |_, mapping| {
            assert_eq!(
                mapping.len(),
                4,
                "source callback receives a completed mapping"
            );
            calls_in_final_check.fetch_add(1, Ordering::SeqCst);
            true
        })),
        ..SubstructMatchParams::default()
    };

    assert!(
        get_substruct_matches_with_params(&target, &query, &params).is_empty(),
        "the later tetrahedral final check rejects the opposite center"
    );
    assert_eq!(
        final_check_calls.load(Ordering::SeqCst),
        1,
        "source extraFinalCheck runs before tetrahedral final-check rejection"
    );
}

#[test]
fn q86_query_target_atom_dispatch_entries_callbacks_and_origins() {
    let mut query_atom = QueryAtom::from_parts(
        Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(8)),
    );
    query_atom
        .set_prop("gate", "query")
        .expect("fixed query property is valid");
    let query = q86_query_graph(vec![query_atom], Vec::new());
    let query_before = query.clone();

    let mut current_atom = Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C));
    current_atom
        .set_prop("gate", "target")
        .expect("fixed target property is valid");
    let topology =
        TopologyBlock::try_from_parts(vec![current_atom], Vec::new(), Vec::new(), Vec::new())
            .expect("fixed Q86 atom target topology is valid");
    let topology_before = topology.clone();
    let coordinates = CoordinateBlock::default();
    let explicit_rows = q86_query_graph(
        vec![QueryAtom::from_parts(
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::O)),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(8)),
        )],
        Vec::new(),
    );
    let explicit_rows_before = explicit_rows.clone();
    let state =
        QueryStateRef::try_for_topology(explicit_rows.atoms(), explicit_rows.bonds(), &topology)
            .expect("query rows align with the current target topology");
    let target = SearchTarget::new(&topology, &coordinates, &topology.stereo_groups, None, None)
        .try_with_query_state(state)
        .expect("validated query rows attach to the target");

    assert!(
        get_substruct_matches_with_params(&target, &query, &SubstructMatchParams::default())
            .is_empty(),
        "with the option disabled, the current carbon carrier rejects the oxygen predicate"
    );
    let query_query = SubstructMatchParams {
        use_query_query_matches: true,
        ..SubstructMatchParams::default()
    };
    assert_eq!(
        try_get_substruct_matches_with_params(&target, &query, &query_query)
            .expect("supported atom query-query comparison succeeds")
            .into_iter()
            .map(|matched| matched.atom_mapping)
            .collect::<Vec<_>>(),
        vec![vec![0]],
        "the enabled both-explicit path compares predicates instead of stale row carriers"
    );

    let carrier_rows = q86_query_graph(
        vec![QueryAtom::from_carrier_parts(
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(8)),
        )],
        Vec::new(),
    );
    let carrier_state =
        QueryStateRef::try_for_topology(carrier_rows.atoms(), carrier_rows.bonds(), &topology)
            .expect("carrier-derived rows align with the target");
    let carrier_target =
        SearchTarget::new(&topology, &coordinates, &topology.stereo_groups, None, None)
            .try_with_query_state(carrier_state)
            .expect("carrier-derived rows attach to the target");
    assert!(
        get_substruct_matches_with_params(&carrier_target, &query, &query_query).is_empty(),
        "a carrier-derived target row retains ordinary current-carrier matching"
    );

    let carrier_query = q86_query_graph(
        vec![QueryAtom::from_carrier_parts(
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(8)),
        )],
        Vec::new(),
    );
    assert_eq!(
        get_substruct_matches_with_params(&target, &carrier_query, &query_query).len(),
        1,
        "a carrier-derived query row retains the ordinary atom compatibility path"
    );

    let post_calls = Arc::new(AtomicUsize::new(0));
    let calls_in_post = Arc::clone(&post_calls);
    let post_params = SubstructMatchParams {
        use_query_query_matches: true,
        atom_properties: vec!["gate".to_owned()],
        extra_atom_check: Some(Arc::new(move |_, _, _, _| {
            calls_in_post.fetch_add(1, Ordering::SeqCst);
            true
        })),
        ..SubstructMatchParams::default()
    };
    assert!(get_substruct_matches_with_params(&target, &query, &post_params).is_empty());
    assert_eq!(post_calls.load(Ordering::SeqCst), 0);

    let override_calls = Arc::new(AtomicUsize::new(0));
    let calls_in_override = Arc::clone(&override_calls);
    let override_params = SubstructMatchParams {
        use_query_query_matches: true,
        atom_properties: vec!["gate".to_owned()],
        extra_atom_check: Some(Arc::new(move |_, _, _, _| {
            calls_in_override.fetch_add(1, Ordering::SeqCst);
            true
        })),
        extra_atom_check_overrides_default_check: true,
        ..SubstructMatchParams::default()
    };
    assert_eq!(
        get_substruct_matches_with_params(&target, &query, &override_params).len(),
        1
    );
    assert_eq!(override_calls.load(Ordering::SeqCst), 1);

    let entry_query = q86_query_graph(
        vec![QueryAtom::from_parts(
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
        )],
        Vec::new(),
    );
    let first = get_substruct_match(&target, &entry_query);
    let all = get_substruct_matches(&target, &entry_query);
    let compiled = compile_query(&entry_query)
        .expect("entry query compiles")
        .matches_target(&target)
        .expect("compiled entry accepts the validated target");
    assert_eq!(first.as_ref(), all.first());
    assert_eq!(compiled, all);

    assert_eq!(query, query_before);
    assert_eq!(topology, topology_before);
    assert_eq!(explicit_rows, explicit_rows_before);
}

#[test]
fn q86_query_target_bond_dispatch_and_mixed_origin_fallback() {
    let query = q86_query_graph(
        vec![
            QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C)),
            QueryAtom::new(AtomId::new(1), AtomSpec::new(Element::C)),
        ],
        vec![QueryBond::from_parts(
            Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            ),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
        )],
    );
    let query_before = query.clone();
    let topology = topology(&[Element::C, Element::C], &[(0, 1)]);
    let topology_before = topology.clone();
    let coordinates = CoordinateBlock::default();
    let explicit_rows = q86_query_graph(
        vec![
            QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C)),
            QueryAtom::new(AtomId::new(1), AtomSpec::new(Element::C)),
        ],
        vec![QueryBond::from_parts(
            Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Double),
            ),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Double)),
        )],
    );
    let state =
        QueryStateRef::try_for_topology(explicit_rows.atoms(), explicit_rows.bonds(), &topology)
            .expect("explicit bond rows align with the current topology");
    let target = SearchTarget::new(&topology, &coordinates, &topology.stereo_groups, None, None)
        .try_with_query_state(state)
        .expect("explicit bond rows attach to the target");
    assert_eq!(get_substruct_matches(&target, &query).len(), 1);
    let query_query = SubstructMatchParams {
        use_query_query_matches: true,
        ..SubstructMatchParams::default()
    };
    assert!(get_substruct_matches_with_params(&target, &query, &query_query).is_empty());

    let carrier_rows = q86_query_graph(
        explicit_rows.atoms().to_vec(),
        vec![QueryBond::from_carrier_parts(
            Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Double),
            ),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Double)),
        )],
    );
    let carrier_state =
        QueryStateRef::try_for_topology(carrier_rows.atoms(), carrier_rows.bonds(), &topology)
            .expect("carrier-derived bond rows align with the current topology");
    let carrier_target =
        SearchTarget::new(&topology, &coordinates, &topology.stereo_groups, None, None)
            .try_with_query_state(carrier_state)
            .expect("carrier-derived bond rows attach to the target");
    assert_eq!(
        get_substruct_matches_with_params(&carrier_target, &query, &query_query).len(),
        1,
        "the mixed-origin bond path reads the current single-bond carrier"
    );
    assert_eq!(query, query_before);
    assert_eq!(topology, topology_before);
}

#[test]
fn q86_query_target_recursive_reuse_keeps_cache_call_local() {
    let inner = q86_query_graph(
        vec![QueryAtom::from_parts(
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
        )],
        Vec::new(),
    );
    let query = q86_query_graph(vec![recursive_query_atom(0, inner, 860)], Vec::new());
    let coordinates = CoordinateBlock::default();
    let carbon_topology = topology(&[Element::C], &[]);
    let nitrogen_topology = topology(&[Element::N], &[]);
    let carbon_rows = q86_query_graph(
        vec![QueryAtom::from_carrier_parts(
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(8)),
        )],
        Vec::new(),
    );
    let nitrogen_rows = q86_query_graph(
        vec![QueryAtom::from_carrier_parts(
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::N)),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
        )],
        Vec::new(),
    );
    let carbon_state =
        QueryStateRef::try_for_topology(carbon_rows.atoms(), carbon_rows.bonds(), &carbon_topology)
            .expect("carbon query rows align");
    let nitrogen_state = QueryStateRef::try_for_topology(
        nitrogen_rows.atoms(),
        nitrogen_rows.bonds(),
        &nitrogen_topology,
    )
    .expect("nitrogen query rows align");
    let carbon_target = SearchTarget::new(
        &carbon_topology,
        &coordinates,
        &carbon_topology.stereo_groups,
        None,
        None,
    )
    .try_with_query_state(carbon_state)
    .expect("carbon query state attaches");
    let nitrogen_target = SearchTarget::new(
        &nitrogen_topology,
        &coordinates,
        &nitrogen_topology.stereo_groups,
        None,
        None,
    )
    .try_with_query_state(nitrogen_state)
    .expect("nitrogen query state attaches");
    let params = SubstructMatchParams {
        use_query_query_matches: true,
        ..SubstructMatchParams::default()
    };

    assert_eq!(
        get_substruct_matches_with_params(&carbon_target, &query, &params).len(),
        1
    );
    assert!(get_substruct_matches_with_params(&nitrogen_target, &query, &params).is_empty());
    assert_eq!(
        get_substruct_matches_with_params(&carbon_target, &query, &params).len(),
        1
    );
}

#[test]
fn q86_query_target_remap_and_attachment_errors_are_typed() {
    let old_topology = topology(&[Element::C, Element::O], &[]);
    let old_rows = q86_query_graph(
        vec![
            QueryAtom::from_parts(
                Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(7)),
            ),
            QueryAtom::from_carrier_parts(
                Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::O)),
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(8)),
            ),
        ],
        Vec::new(),
    );
    let old_state =
        QueryStateRef::try_for_topology(old_rows.atoms(), old_rows.bonds(), &old_topology)
            .expect("old query rows align");
    let (reordered, mapping) = old_topology
        .reordered_atoms(&[AtomId::new(1), AtomId::new(0)])
        .expect("fixed atom permutation is valid");
    let (remapped_atoms, remapped_bonds) =
        remap_query_rows(old_state, &reordered, &mapping).expect("query rows remap");
    assert!(remapped_atoms[0].predicate_is_carrier_derived());
    assert!(!remapped_atoms[1].predicate_is_carrier_derived());
    assert_eq!(
        remapped_atoms[1].predicate(),
        &QueryNode::predicate(AtomQueryPredicate::AtomicNumber(7))
    );
    let remapped_state =
        QueryStateRef::try_for_topology(&remapped_atoms, &remapped_bonds, &reordered)
            .expect("remapped rows revalidate against the reordered topology");
    let coordinates = CoordinateBlock::default();
    let target = SearchTarget::new(
        &reordered,
        &coordinates,
        &reordered.stereo_groups,
        None,
        None,
    )
    .try_with_query_state(remapped_state)
    .expect("revalidated remapped state attaches");
    let nitrogen_query = q86_query_graph(
        vec![QueryAtom::from_parts(
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::N)),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(7)),
        )],
        Vec::new(),
    );
    let params = SubstructMatchParams {
        use_query_query_matches: true,
        ..SubstructMatchParams::default()
    };
    let matches = try_get_substruct_matches_with_params(&target, &nitrogen_query, &params)
        .expect("remapped target query state matches through the real entry point");
    assert_eq!(matches.len(), 1);
    assert_eq!(matches[0].atom_mapping, [1]);

    let short_topology = topology(&[Element::C], &[]);
    let count_error = SearchTarget::new(
        &short_topology,
        &coordinates,
        &short_topology.stereo_groups,
        None,
        None,
    )
    .try_with_query_state(old_state)
    .expect_err("attachment must revalidate row counts before matching");
    assert_eq!(
        count_error,
        QueryStateError::AtomCount {
            actual: 2,
            expected: 1,
        }
    );

    let reverse_topology = topology(&[Element::C, Element::C], &[(1, 0)]);
    let reverse_rows = q86_query_graph(
        vec![
            QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C)),
            QueryAtom::new(AtomId::new(1), AtomSpec::new(Element::C)),
        ],
        vec![QueryBond::new(
            BondId::new(0),
            BondSpec::new(AtomId::new(1), AtomId::new(0), BondOrder::Single),
        )],
    );
    let reverse_state = QueryStateRef::try_for_topology(
        reverse_rows.atoms(),
        reverse_rows.bonds(),
        &reverse_topology,
    )
    .expect("reverse-endpoint rows align with their source topology");
    let forward_topology = topology(&[Element::C, Element::C], &[(0, 1)]);
    let endpoint_error = SearchTarget::new(
        &forward_topology,
        &coordinates,
        &forward_topology.stereo_groups,
        None,
        None,
    )
    .try_with_query_state(reverse_state)
    .expect_err("attachment must revalidate bond endpoints before matching");
    assert_eq!(
        endpoint_error,
        QueryStateError::BondEndpoints {
            position: 0,
            actual: (AtomId::new(1), AtomId::new(0)),
            expected: (AtomId::new(0), AtomId::new(1)),
        }
    );
}
