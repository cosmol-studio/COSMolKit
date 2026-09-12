use std::collections::BTreeMap;

use cosmolkit_model::{
    AtomId, AtomQueryPredicate, AtomRangeBounds, AtomRangeDataFunction, AtomRangeQuery, AtomSpec,
    Bond, BondDirection, BondId, BondOrder, BondQueryPredicate, BondSpec, BondStereo, ChiralTag,
    Conformer2D, Conformer3D, CoordinateValidationError, Element, Hybridization, QueryAtom,
    QueryBond, QueryGraph, QueryGraphError, QueryNode, RecursiveStructureQuery, StereoGroup,
    StereoGroupKind,
};

fn carbon(id: usize) -> QueryAtom {
    QueryAtom::new(AtomId::new(id), AtomSpec::new(Element::C))
}

fn single_bond(id: usize, begin: usize, end: usize) -> QueryBond {
    QueryBond::new(
        BondId::new(id),
        BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
    )
}

fn graph(atom_count: usize, bonds: Vec<QueryBond>) -> QueryGraph {
    QueryGraph::from_parts(
        (0..atom_count).map(carbon).collect(),
        bonds,
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .unwrap()
}

#[test]
fn query_node_covers_all_shapes_children_and_negation_transitions() {
    let carbon = QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6));
    let oxygen = QueryNode::predicate(AtomQueryPredicate::AtomicNumber(8));
    assert!(matches!(carbon, QueryNode::Predicate(_)));

    let mut and = QueryNode::and(vec![carbon.clone()]);
    let mut or = QueryNode::or(Vec::new());
    let mut xor = QueryNode::xor(Vec::new());
    and.add_child(oxygen.clone());
    or.add_child(carbon.clone());
    xor.add_child(oxygen.clone());
    assert_eq!(and, QueryNode::And(vec![carbon.clone(), oxygen.clone()]));
    assert_eq!(or, QueryNode::Or(vec![carbon.clone()]));
    assert_eq!(xor, QueryNode::Xor(vec![oxygen]));

    let original = and.clone();
    and.set_negation(false);
    assert_eq!(and, original);
    and.set_negation(true);
    assert!(and.is_negated());
    assert_eq!(and, QueryNode::not(original.clone()));
    and.set_negation(true);
    assert_eq!(and, QueryNode::not(original.clone()));
    and.set_negation(false);
    assert_eq!(and, original);

    let child_negation = QueryNode::not(carbon);
    let mut nested = QueryNode::not(child_negation.clone());
    nested.set_negation(false);
    assert_eq!(nested, child_negation);
}

#[test]
fn range_query_covers_every_bound_and_data_function_value() {
    let bounds = [
        AtomRangeBounds::LessEqual(2),
        AtomRangeBounds::GreaterEqual(3),
        AtomRangeBounds::Inclusive {
            lower: 1,
            upper: 4,
            lower_open: true,
            upper_open: false,
        },
    ];
    let functions = [
        AtomRangeDataFunction::ExplicitDegree,
        AtomRangeDataFunction::NonHydrogenDegree,
        AtomRangeDataFunction::TotalDegree,
        AtomRangeDataFunction::TotalValence,
        AtomRangeDataFunction::NumAtomRings,
        AtomRangeDataFunction::NumHeteroatomNeighbors,
        AtomRangeDataFunction::NumAliphaticHeteroatomNeighbors,
        AtomRangeDataFunction::MinRingSize,
        AtomRangeDataFunction::RingBondCount,
        AtomRangeDataFunction::ImplicitHydrogenCount,
        AtomRangeDataFunction::FormalCharge,
        AtomRangeDataFunction::NegativeFormalCharge,
        AtomRangeDataFunction::AtomRingSize {
            lower: 3,
            upper: 7,
            lower_open: false,
            upper_open: true,
        },
    ];

    for bound in bounds {
        for function in functions {
            let query = AtomRangeQuery::new(bound, function);
            assert_eq!(query.bounds(), bound);
            assert_eq!(query.data_function(), function);
            assert_eq!(query.writer_parts(), (bound, function));
        }
    }
}

#[test]
fn atom_predicate_family_has_all_audited_variants() {
    let range = AtomRangeQuery::new(
        AtomRangeBounds::LessEqual(2),
        AtomRangeDataFunction::ExplicitDegree,
    );
    let recursive = RecursiveStructureQuery::from_query_graph(graph(0, Vec::new()), 4);
    let predicates = vec![
        AtomQueryPredicate::Any,
        AtomQueryPredicate::AtomicNumber(6),
        AtomQueryPredicate::AtomType {
            atomic_number: 6,
            aromatic: true,
        },
        AtomQueryPredicate::AtomicNumberIn(vec![6, 7]),
        AtomQueryPredicate::AtomicNumberNotIn(vec![8]),
        AtomQueryPredicate::FormalCharge(-1),
        AtomQueryPredicate::NegativeFormalCharge(1),
        AtomQueryPredicate::NumRadicalElectrons(1),
        AtomQueryPredicate::HasChiralTag,
        AtomQueryPredicate::MissingChiralTag,
        AtomQueryPredicate::Isotope(13),
        AtomQueryPredicate::HydrogenCount(2),
        AtomQueryPredicate::HasImplicitHydrogen,
        AtomQueryPredicate::ImplicitHydrogenCount(1),
        AtomQueryPredicate::ImplicitHydrogenCountLessEqual(2),
        AtomQueryPredicate::ImplicitValence(4),
        AtomQueryPredicate::ExplicitValence(3),
        AtomQueryPredicate::ExplicitDegree(2),
        AtomQueryPredicate::ExplicitDegreeLessEqual(3),
        AtomQueryPredicate::NonHydrogenDegree(2),
        AtomQueryPredicate::NonHydrogenDegreeLessEqual(3),
        AtomQueryPredicate::NonHydrogenDegreeGreaterEqual(1),
        AtomQueryPredicate::HeavyAtomDegree(2),
        AtomQueryPredicate::NumHeteroatomNeighbors(1),
        AtomQueryPredicate::HasHeteroatomNeighbors,
        AtomQueryPredicate::NumAliphaticHeteroatomNeighbors(1),
        AtomQueryPredicate::HasAliphaticHeteroatomNeighbors,
        AtomQueryPredicate::RingBondCount(2),
        AtomQueryPredicate::RingBondCountLessEqual(2),
        AtomQueryPredicate::HasRingBond,
        AtomQueryPredicate::IsBridgehead,
        AtomQueryPredicate::IsAromatic(true),
        AtomQueryPredicate::IsUnsaturated,
        AtomQueryPredicate::RecursiveSmarts(recursive),
        AtomQueryPredicate::HasProperty("flag".into()),
        AtomQueryPredicate::PropertyValue {
            name: "key".into(),
            value: "value".into(),
        },
        AtomQueryPredicate::RGroupLabel(2),
        AtomQueryPredicate::MolFileAlias("Q".into()),
        AtomQueryPredicate::HybridizationMatch(Hybridization::Sp3),
        AtomQueryPredicate::TotalDegree(4),
        AtomQueryPredicate::TotalDegreeLessEqual(4),
        AtomQueryPredicate::TotalDegreeGreaterEqual(1),
        AtomQueryPredicate::TotalValence(4),
        AtomQueryPredicate::TotalValenceLessEqual(4),
        AtomQueryPredicate::TotalValenceGreaterEqual(1),
        AtomQueryPredicate::InRing,
        AtomQueryPredicate::NumAtomRings(2),
        AtomQueryPredicate::InRingOfSize(6),
        AtomQueryPredicate::InRingOfSizeLessEqual(6),
        AtomQueryPredicate::InRingOfSizeGreaterEqual(3),
        AtomQueryPredicate::SmallestRingSize(5),
        AtomQueryPredicate::SmallestRingSizeLessEqual(6),
        AtomQueryPredicate::SmallestRingSizeGreaterEqual(3),
        AtomQueryPredicate::Mass(12),
        AtomQueryPredicate::ChiralTagMatch(ChiralTag::TetrahedralCw),
        AtomQueryPredicate::ChiralPermutationMatch(2),
        AtomQueryPredicate::DegreeLessEqual(3),
        AtomQueryPredicate::DegreeGreaterEqual(1),
        AtomQueryPredicate::Range(range),
        AtomQueryPredicate::UnsupportedFeature("future atom leaf"),
    ];
    assert_eq!(predicates.len(), 60);
}

#[test]
fn bond_predicate_family_has_all_audited_variants() {
    let predicates = [
        BondQueryPredicate::Any,
        BondQueryPredicate::Order(BondOrder::Single),
        BondQueryPredicate::OrderIn(vec![BondOrder::Single, BondOrder::Aromatic]),
        BondQueryPredicate::IsAromatic(true),
        BondQueryPredicate::IsInRing(true),
        BondQueryPredicate::Direction(BondDirection::EndUpRight),
        BondQueryPredicate::Stereo(BondStereo::Cis),
        BondQueryPredicate::HasStereo,
        BondQueryPredicate::IsConjugated,
        BondQueryPredicate::NumRingBonds(2),
        BondQueryPredicate::InRingOfSize(6),
        BondQueryPredicate::MinRingSize(5),
        BondQueryPredicate::NumRingBondsGreaterEqual(1),
        BondQueryPredicate::NumRingBondsLessEqual(2),
        BondQueryPredicate::MolFileQueryCode(5),
        BondQueryPredicate::HasProperty("flag".into()),
        BondQueryPredicate::PropertyValue {
            name: "key".into(),
            value: "value".into(),
        },
        BondQueryPredicate::UnsupportedFeature("future bond leaf"),
    ];
    assert_eq!(predicates.len(), 18);
}

#[test]
fn query_atom_and_bond_cover_default_parts_access_and_mutation() {
    let atom_spec = AtomSpec::new(Element::N)
        .with_atom_map(8)
        .with_prop("a", "b")
        .unwrap();
    let mut atom = QueryAtom::new(AtomId::new(0), atom_spec);
    assert_eq!(atom.id(), AtomId::new(0));
    assert_eq!(atom.index(), 0);
    assert_eq!(atom.atom_map(), Some(8));
    assert_eq!(atom.prop("a"), Some("b"));
    assert_eq!(
        atom.predicate(),
        &QueryNode::predicate(AtomQueryPredicate::AtomicNumber(7))
    );
    atom.atom_mut().set_atom_map(Some(9));
    atom.set_predicate(QueryNode::predicate(AtomQueryPredicate::Any));
    *atom.predicate_mut() = QueryNode::not(QueryNode::predicate(AtomQueryPredicate::Any));
    assert_eq!(atom.atom_map(), Some(9));
    assert!(matches!(atom.predicate(), QueryNode::Not(_)));

    let unspecified = QueryBond::new(
        BondId::new(0),
        BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Unspecified),
    );
    assert_eq!(
        unspecified.predicate(),
        &QueryNode::predicate(BondQueryPredicate::Any)
    );

    let concrete = Bond::from_spec(
        BondId::new(1),
        BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Double),
    );
    let mut bond = QueryBond::from_parts(
        concrete,
        QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Double)),
    );
    assert_eq!(bond.id(), BondId::new(1));
    assert_eq!(bond.endpoints(), (1, 2));
    assert_eq!(bond.begin(), AtomId::new(1));
    assert_eq!(bond.end(), AtomId::new(2));
    bond.bond_mut().set_order(BondOrder::Triple);
    bond.set_predicate(QueryNode::predicate(BondQueryPredicate::Any));
    *bond.predicate_mut() = QueryNode::not(QueryNode::predicate(BondQueryPredicate::Any));
    assert_eq!(bond.bond().order(), BondOrder::Triple);
    assert!(matches!(bond.predicate(), QueryNode::Not(_)));
}

#[test]
fn recursive_query_covers_empty_graph_provenance_set_serial_and_deep_clone() {
    let empty = RecursiveStructureQuery::new();
    assert!(empty.query_graph().is_none());
    assert_eq!(empty.source_smarts(), None);
    assert_eq!(empty.serial_number(), 0);

    let mut recursive = RecursiveStructureQuery::from_query_graph(graph(1, Vec::new()), 17)
        .with_source_smarts("[#6]");
    assert_eq!(recursive.query_graph().unwrap().num_atoms(), 1);
    assert_eq!(recursive.source_smarts(), Some("[#6]"));
    assert_eq!(recursive.serial_number(), 17);
    recursive.insert_atom_index(3);
    assert!(recursive.contains_atom_index(3));
    recursive
        .query_graph_mut()
        .unwrap()
        .set_prop("changed", "first");

    let clone = recursive.clone();
    recursive
        .query_graph_mut()
        .unwrap()
        .set_prop("changed", "second");
    assert_eq!(clone.query_graph().unwrap().prop("changed"), Some("first"));
    assert_eq!(clone, clone.clone());

    recursive.set_query_graph(graph(2, vec![single_bond(0, 0, 1)]));
    assert_eq!(recursive.query_graph().unwrap().num_bonds(), 1);
}

#[test]
fn query_graph_covers_accessors_properties_coordinates_and_stereo_groups() {
    let mut graph = graph(2, vec![single_bond(0, 0, 1)])
        .with_name("ethane query")
        .with_prop("origin", "test")
        .with_2d_coordinate_block(vec![[0.0, 0.0], [1.0, 0.0]])
        .unwrap();
    graph.set_prop("mutable", "yes");
    graph
        .add_conformer_3d(Conformer3D::new(
            4,
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            true,
        ))
        .unwrap();
    graph.add_stereo_group(StereoGroup::new(
        StereoGroupKind::Absolute,
        vec![AtomId::new(0)],
        Vec::new(),
    ));

    assert_eq!(graph.num_atoms(), 2);
    assert_eq!(graph.num_bonds(), 1);
    assert_eq!(graph.atoms().len(), 2);
    assert_eq!(graph.atom(0).unwrap().id(), AtomId::new(0));
    assert!(graph.atom(2).is_none());
    assert_eq!(graph.bonds().len(), 1);
    assert!(graph.bond(1).is_none());
    assert_eq!(graph.adjacency(), &[vec![(1, 0)], vec![(0, 0)]]);
    assert_eq!(graph.name(), Some("ethane query"));
    assert_eq!(graph.prop("origin"), Some("test"));
    assert_eq!(
        graph.props().get("mutable").map(String::as_str),
        Some("yes")
    );
    assert_eq!(graph.coordinates_2d().unwrap().len(), 2);
    assert_eq!(graph.conformers_3d()[0].id(), 4);
    assert_eq!(graph.stereo_groups().len(), 1);
    graph.clear_prop("mutable");
    assert_eq!(graph.prop("mutable"), None);
    assert_eq!(graph.validate(), Ok(()));
}

#[test]
fn query_graph_reports_canonical_id_and_endpoint_errors() {
    let atom_error = QueryGraph::from_parts(
        vec![carbon(1)],
        Vec::new(),
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    );
    assert!(matches!(
        atom_error,
        Err(QueryGraphError::AtomIdMismatch { position: 0, id }) if id == AtomId::new(1)
    ));

    let bond_id_error = QueryGraph::from_parts(
        vec![carbon(0), carbon(1)],
        vec![single_bond(1, 0, 1)],
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    );
    assert!(matches!(
        bond_id_error,
        Err(QueryGraphError::BondIdMismatch { position: 0, id }) if id == BondId::new(1)
    ));

    let endpoint_error = QueryGraph::from_parts(
        vec![carbon(0)],
        vec![single_bond(0, 0, 1)],
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    );
    assert_eq!(endpoint_error, Err(QueryGraphError::InvalidBondEndpoint(0)));
}

#[test]
fn query_graph_reports_exact_endpoint_error_after_detached_mutation() {
    let mut mutated = graph(2, vec![single_bond(0, 0, 1)]);
    mutated.bonds_mut()[0]
        .bond_mut()
        .set_endpoints(AtomId::new(0), AtomId::new(2));

    assert_eq!(
        mutated.validate(),
        Err(QueryGraphError::InvalidBondEndpoint(0))
    );
}

#[test]
fn query_graph_reports_stereo_coordinate_and_adjacency_errors() {
    let stereo_atom_error = QueryGraph::from_parts(
        vec![carbon(0), carbon(1)],
        vec![QueryBond::new(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Double)
                .with_stereo_atoms(AtomId::new(0), AtomId::new(2)),
        )],
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    );
    assert!(matches!(
        stereo_atom_error,
        Err(QueryGraphError::StereoAtomOutOfRange { .. })
    ));

    let coordinate_error = QueryGraph::from_parts(
        vec![carbon(0)],
        Vec::new(),
        BTreeMap::new(),
        vec![Conformer2D::new(7, Vec::new())],
        Vec::new(),
        Vec::new(),
    );
    assert!(matches!(
        coordinate_error,
        Err(QueryGraphError::CoordinateValidation(
            CoordinateValidationError::RowCount {
                dimension: "2D",
                conformer: 7,
                ..
            }
        ))
    ));

    let duplicate_coordinate_error = QueryGraph::from_parts(
        vec![carbon(0)],
        Vec::new(),
        BTreeMap::new(),
        Vec::new(),
        vec![
            Conformer3D::new(4, vec![[0.0, 0.0, 0.0]], true),
            Conformer3D::new(4, vec![[1.0, 0.0, 0.0]], true),
        ],
        Vec::new(),
    );
    assert!(matches!(
        duplicate_coordinate_error,
        Err(QueryGraphError::CoordinateValidation(
            CoordinateValidationError::DuplicateConformerId {
                dimension: "3D",
                id: 4
            }
        ))
    ));

    let stereo_group_atom_error = QueryGraph::from_parts(
        vec![carbon(0)],
        Vec::new(),
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        vec![StereoGroup::new(
            StereoGroupKind::Absolute,
            vec![AtomId::new(1)],
            Vec::new(),
        )],
    );
    assert!(matches!(
        stereo_group_atom_error,
        Err(QueryGraphError::StereoGroupAtomOutOfRange { .. })
    ));

    let stereo_group_bond_error = QueryGraph::from_parts(
        vec![carbon(0)],
        Vec::new(),
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        vec![StereoGroup::new(
            StereoGroupKind::Absolute,
            Vec::new(),
            vec![BondId::new(0)],
        )],
    );
    assert!(matches!(
        stereo_group_bond_error,
        Err(QueryGraphError::StereoGroupBondOutOfRange { .. })
    ));

    let mut stale = graph(3, vec![single_bond(0, 0, 1)]);
    stale.bonds_mut()[0]
        .bond_mut()
        .set_endpoints(AtomId::new(1), AtomId::new(2));
    assert_eq!(stale.validate(), Err(QueryGraphError::AdjacencyMismatch));
}

#[test]
fn query_graph_coordinate_builders_return_structured_errors() {
    let mut graph = graph(2, vec![single_bond(0, 0, 1)]);
    let error = graph
        .add_conformer_3d(Conformer3D::new(9, vec![[0.0, 0.0, 0.0]], true))
        .unwrap_err();
    assert!(matches!(
        error,
        QueryGraphError::CoordinateValidation(CoordinateValidationError::RowCount {
            dimension: "3D",
            conformer: 9,
            rows: 1,
            atom_count: 2,
        })
    ));

    let error = graph
        .with_2d_coordinate_block(vec![[0.0, 0.0]])
        .unwrap_err();
    assert!(matches!(
        error,
        QueryGraphError::CoordinateValidation(CoordinateValidationError::RowCount {
            dimension: "2D",
            conformer: 0,
            rows: 1,
            atom_count: 2,
        })
    ));
}
