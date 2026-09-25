use std::collections::BTreeMap;

use cosmolkit_model::{
    Atom, AtomId, AtomPdbResidueInfo, AtomQueryPredicate, AtomRangeBounds, AtomRangeDataFunction,
    AtomRangeQuery, AtomSpec, Bond, BondDirection, BondId, BondOrder, BondQueryPredicate, BondSpec,
    BondStereo, ChiralTag, Conformer2D, Conformer3D, CoordinateValidationError, Element,
    Hybridization, QueryAtom, QueryAtomConversionError, QueryAtomIdentity, QueryBond, QueryGraph,
    QueryGraphError, QueryNode, QueryStateError, QueryStateRef, RecursiveStructureQuery,
    StereoGroup, StereoGroupKind, TemplateAttachment, TemplateAttachmentOrder, TopologyBlock,
    remap_query_rows,
};

fn carbon(id: usize) -> QueryAtom {
    QueryAtom::new(AtomId::new(id), AtomSpec::new(Element::C))
}

#[test]
fn query_origin_equality_compares_representation_not_matching_equivalence() {
    let atom = Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C));
    let predicate = QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6));
    let explicit = QueryAtom::from_parts(atom.clone(), predicate.clone());
    let carrier = QueryAtom::from_carrier_parts(atom, predicate);
    assert_eq!(
        explicit.try_to_atom().unwrap(),
        carrier.try_to_atom().unwrap()
    );
    assert_eq!(explicit.predicate(), carrier.predicate());
    assert_ne!(explicit, carrier);
    assert_eq!(explicit, explicit.clone());
    let make_graph = |row| {
        QueryGraph::from_parts(vec![row], vec![], BTreeMap::new(), vec![], vec![], vec![]).unwrap()
    };
    let explicit_graph = make_graph(explicit);
    let carrier_graph = make_graph(carrier);
    assert_ne!(explicit_graph, carrier_graph);
    assert_eq!(explicit_graph, explicit_graph.clone());

    let bond = Bond::from_spec(
        BondId::new(0),
        BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
    );
    let predicate = QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single));
    let explicit = QueryBond::from_parts(bond.clone(), predicate.clone());
    let carrier = QueryBond::from_carrier_parts(bond, predicate);
    assert_eq!(explicit.bond(), carrier.bond());
    assert_eq!(explicit.predicate(), carrier.predicate());
    assert_ne!(explicit, carrier);
    assert_eq!(explicit, explicit.clone());
}

#[test]
fn query_overlay_remap_uses_current_carriers_without_rewriting_explicit_predicates() {
    let atoms = vec![carbon(0), carbon(1)];
    let bonds = vec![single_bond(0, 0, 1)];
    let mut current = TopologyBlock::try_from_parts(
        atoms.iter().map(|a| a.try_to_atom().unwrap()).collect(),
        bonds.iter().map(|b| b.bond().clone()).collect(),
        vec![],
        vec![],
    )
    .unwrap();
    let state = QueryStateRef::try_for_topology(&atoms, &bonds, &current).unwrap();
    current.atoms[0].set_formal_charge(1);
    current.bonds[0].set_order(BondOrder::Double);
    state.validate_for_topology(&current).unwrap();
    assert_eq!(state.atom_predicate(AtomId::new(0)), atoms[0].predicate());
    assert_eq!(state.bond_predicate(BondId::new(0)), bonds[0].predicate());
    let (updated_atoms, updated_bonds) = remap_query_rows(
        state,
        &current,
        &cosmolkit_model::TopologyMapping::identity(2, 1),
    )
    .unwrap();
    assert_eq!(updated_atoms[0].try_to_atom().unwrap(), current.atoms[0]);
    assert_ne!(
        updated_atoms[0].try_to_atom().unwrap(),
        atoms[0].try_to_atom().unwrap()
    );
    assert_eq!(updated_bonds[0].bond(), &current.bonds[0]);
    assert_ne!(updated_bonds[0].bond(), bonds[0].bond());
    assert_eq!(updated_atoms[0].predicate(), atoms[0].predicate());
    assert_eq!(updated_bonds[0].predicate(), bonds[0].predicate());
    let updated =
        QueryStateRef::try_for_topology(&updated_atoms, &updated_bonds, &current).unwrap();
    assert!(updated.atom_has_query(AtomId::new(0)));
    assert!(updated.bond_has_query(BondId::new(0)));
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
fn q05_query_state_transport_preserves_explicit_and_carrier_origins_and_trees() {
    let explicit_atom_tree = QueryNode::and(vec![
        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(7)),
        QueryNode::not(QueryNode::predicate(AtomQueryPredicate::IsAromatic(true))),
    ]);
    let explicit_atom = QueryAtom::from_parts(
        Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
        explicit_atom_tree.clone(),
    );
    let carrier_atom = QueryAtom::from_carrier_parts(
        Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::H)),
        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(1)),
    );
    let explicit_bond_tree = QueryNode::or(vec![
        QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
        QueryNode::and(vec![
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Double)),
            QueryNode::not(QueryNode::predicate(BondQueryPredicate::IsInRing(true))),
        ]),
    ]);
    let explicit_bond = QueryBond::from_parts(
        Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        ),
        explicit_bond_tree.clone(),
    );

    let cloned_atom = explicit_atom.clone();
    let cloned_carrier = carrier_atom.clone();
    let cloned_bond = explicit_bond.clone();
    let topology = TopologyBlock::try_from_parts(
        vec![
            cloned_atom.try_to_atom().unwrap(),
            cloned_carrier.try_to_atom().unwrap(),
        ],
        vec![cloned_bond.bond().clone()],
        Vec::new(),
        Vec::new(),
    )
    .unwrap();
    let atom_rows = [cloned_atom.clone(), cloned_carrier.clone()];
    let bond_rows = [cloned_bond.clone()];
    let state = QueryStateRef::try_for_topology(&atom_rows, &bond_rows, &topology).unwrap();
    assert!(!cloned_atom.predicate_is_carrier_derived());
    assert!(cloned_carrier.predicate_is_carrier_derived());
    assert!(!cloned_bond.predicate_is_carrier_derived());
    assert_eq!(cloned_atom.predicate(), &explicit_atom_tree);
    assert_eq!(cloned_bond.predicate(), &explicit_bond_tree);
    assert!(state.atom_has_query(AtomId::new(0)));
    assert!(!state.atom_has_query(AtomId::new(1)));
    assert!(state.bond_has_query(BondId::new(0)));
    assert_eq!(state.atom_predicate(AtomId::new(0)), &explicit_atom_tree);
    assert_eq!(state.bond_predicate(BondId::new(0)), &explicit_bond_tree);

    let wrong_atom_rows = [
        QueryAtom::from_parts(
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
            explicit_atom_tree.clone(),
        ),
        cloned_carrier.clone(),
    ];
    assert!(matches!(
        QueryStateRef::try_for_topology(&wrong_atom_rows, &bond_rows, &topology),
        Err(QueryStateError::AtomId { position: 0, .. })
    ));

    let mut promoted = cloned_carrier;
    let prior = promoted.predicate().clone();
    *promoted.predicate_mut() = QueryNode::and(vec![
        prior,
        QueryNode::predicate(AtomQueryPredicate::ExplicitDegree(1)),
    ]);
    assert!(!promoted.predicate_is_carrier_derived());
    assert!(matches!(promoted.predicate(), QueryNode::And(_)));
}

#[test]
fn q05_query_state_transport_distinguishes_explicit_query_h_from_ordinary_carrier_h() {
    let hydrogen = Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::H));
    let explicit = QueryAtom::from_parts(
        hydrogen.clone(),
        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(1)),
    );
    let carrier = QueryAtom::from_carrier_parts(
        hydrogen,
        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(1)),
    );

    assert!(!explicit.predicate_is_carrier_derived());
    assert!(carrier.predicate_is_carrier_derived());
    assert_eq!(explicit.predicate(), carrier.predicate());
    assert!(explicit.prop("_MolFileAtomQuery").is_none());
    assert!(carrier.prop("_MolFileAtomQuery").is_none());
}

#[test]
fn q05_query_state_transport_mapping_baseline_preserves_surviving_typed_rows() {
    let atoms = [Element::C, Element::H, Element::N, Element::O]
        .into_iter()
        .enumerate()
        .map(|(id, element)| Atom::from_spec(AtomId::new(id), AtomSpec::new(element)))
        .collect::<Vec<_>>();
    let bonds = [(0, 1), (0, 2), (2, 3)]
        .into_iter()
        .enumerate()
        .map(|(id, (begin, end))| {
            Bond::from_spec(
                BondId::new(id),
                BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
            )
        })
        .collect::<Vec<_>>();
    let source = TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap();
    let atom_rows = source
        .atoms
        .iter()
        .cloned()
        .enumerate()
        .map(|(index, atom)| {
            if matches!(index, 0 | 3) {
                QueryAtom::from_parts(
                    atom,
                    QueryNode::predicate(AtomQueryPredicate::AtomicNumber(if index == 0 {
                        7
                    } else {
                        8
                    })),
                )
            } else {
                QueryAtom::from_carrier_parts(
                    atom.clone(),
                    QueryNode::predicate(AtomQueryPredicate::AtomicNumber(atom.atomic_number())),
                )
            }
        })
        .collect::<Vec<_>>();
    let compound = QueryNode::and(vec![
        QueryNode::or(vec![
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Double)),
        ]),
        QueryNode::not(QueryNode::predicate(BondQueryPredicate::IsInRing(true))),
    ]);
    let bond_rows = source
        .bonds
        .iter()
        .cloned()
        .enumerate()
        .map(|(index, bond)| {
            if index == 2 {
                QueryBond::from_parts(bond, compound.clone())
            } else {
                QueryBond::from_carrier_parts(
                    bond,
                    QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
                )
            }
        })
        .collect::<Vec<_>>();

    let mut edit = source.begin_batch_edit().unwrap();
    edit.remove_atom(AtomId::new(1)).unwrap();
    let (compacted, mapping) = edit.finish().unwrap();
    let state = QueryStateRef::try_for_topology(&atom_rows, &bond_rows, &source).unwrap();
    let (remapped_atoms, remapped_bonds) = remap_query_rows(state, &compacted, &mapping).unwrap();
    let remapped = QueryGraph::from_parts(
        remapped_atoms,
        remapped_bonds,
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .unwrap();

    assert_eq!(
        mapping.atoms().new_to_old(),
        &[
            Some(AtomId::new(0)),
            Some(AtomId::new(2)),
            Some(AtomId::new(3))
        ]
    );
    assert_eq!(
        remapped.atom(0).unwrap().predicate(),
        &QueryNode::predicate(AtomQueryPredicate::AtomicNumber(7))
    );
    assert!(!remapped.atom(0).unwrap().predicate_is_carrier_derived());
    assert!(remapped.atom(1).unwrap().predicate_is_carrier_derived());
    assert!(!remapped.atom(2).unwrap().predicate_is_carrier_derived());
    assert_eq!(remapped.bond(1).unwrap().predicate(), &compound);
    assert!(!remapped.bond(1).unwrap().predicate_is_carrier_derived());
    assert_eq!(remapped.bond(1).unwrap().endpoints(), (1, 2));
    assert_eq!(remapped.validate(), Ok(()));
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
fn q07e_ring_equality_targets_keep_full_i32_identity_and_carrier_provenance() {
    let atom = Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C));
    for target in [0, 255, 256, 2_147_483_639] {
        let predicates = [
            AtomQueryPredicate::SmallestRingSize(target),
            AtomQueryPredicate::InRingOfSize(target),
            AtomQueryPredicate::RingBondCount(target),
        ];
        assert_ne!(predicates[0], predicates[1], "target {target}");
        assert_ne!(predicates[0], predicates[2], "target {target}");
        assert_ne!(predicates[1], predicates[2], "target {target}");

        for predicate in predicates {
            let query = QueryNode::predicate(predicate.clone());
            assert_eq!(query.clone(), query, "target {target}");
            let explicit = QueryAtom::from_parts(atom.clone(), query.clone());
            let carrier_derived = QueryAtom::from_carrier_parts(atom.clone(), query.clone());
            assert_eq!(explicit.try_to_atom().unwrap(), atom, "target {target}");
            assert_eq!(
                carrier_derived.try_to_atom().unwrap(),
                atom,
                "target {target}"
            );
            assert_eq!(explicit.predicate(), &query, "target {target}");
            assert_eq!(carrier_derived.predicate(), &query, "target {target}");
            assert!(!explicit.predicate_is_carrier_derived(), "target {target}");
            assert!(
                carrier_derived.predicate_is_carrier_derived(),
                "target {target}"
            );
            assert_eq!(explicit.clone(), explicit, "target {target}");
            assert_eq!(carrier_derived.clone(), carrier_derived, "target {target}");
        }
    }

    let bounds = AtomRangeBounds::Inclusive {
        lower: 256,
        upper: 2_147_483_639,
        lower_open: false,
        upper_open: true,
    };
    let data_function = AtomRangeDataFunction::AtomRingSize {
        lower: 256,
        upper: 2_147_483_639,
        lower_open: false,
        upper_open: true,
    };
    let range = AtomRangeQuery::new(bounds, data_function);
    assert_eq!(range.bounds(), bounds);
    assert_eq!(range.data_function(), data_function);
    assert_eq!(range.writer_parts(), (bounds, data_function));
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
    atom.set_atom_map(Some(9));
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

#[test]
fn query_atom_keeps_raw_identity_separate_from_predicate_and_fails_concrete_conversion() {
    let id = AtomId::new(4);
    let predicate = QueryNode::and(vec![
        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(119)),
        QueryNode::predicate(AtomQueryPredicate::Isotope(13)),
    ]);
    let mut atom =
        QueryAtom::from_identity_parts(id, QueryAtomIdentity::AtomicNumber(119), predicate.clone());
    atom.set_isotope(Some(13));
    atom.set_formal_charge(1);

    assert_eq!(atom.id(), id);
    assert_eq!(atom.identity(), QueryAtomIdentity::AtomicNumber(119));
    assert_eq!(atom.atomic_number(), 119);
    assert_eq!(atom.element(), None);
    assert_eq!(atom.isotope(), Some(13));
    assert_eq!(atom.formal_charge(), 1);
    assert_eq!(atom.predicate(), &predicate);
    assert!(!atom.predicate_is_carrier_derived());
    assert_eq!(
        atom.try_to_atom(),
        Err(QueryAtomConversionError::NonElementAtomicNumber {
            atom: id,
            atomic_number: 119,
        })
    );

    assert_eq!(Element::from_atomic_number(119), None);
    let ordinary = Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C));
    assert_eq!(ordinary.element(), Element::C);
}

#[test]
fn widened_query_targets_remain_separate_from_narrow_carriers() {
    let mut positive_charge = QueryAtom::from_identity_parts(
        AtomId::new(10),
        QueryAtomIdentity::AtomicNumber(6),
        QueryNode::predicate(AtomQueryPredicate::FormalCharge(128)),
    );
    positive_charge.set_formal_charge(-128);
    assert_eq!(
        positive_charge.predicate(),
        &QueryNode::predicate(AtomQueryPredicate::FormalCharge(128))
    );
    assert_eq!(positive_charge.formal_charge(), -128);
    assert!(!positive_charge.predicate_is_carrier_derived());
    assert_eq!(positive_charge.clone(), positive_charge);

    let mut negative_charge = QueryAtom::from_identity_parts(
        AtomId::new(11),
        QueryAtomIdentity::AtomicNumber(6),
        QueryNode::predicate(AtomQueryPredicate::FormalCharge(-129)),
    );
    negative_charge.set_formal_charge(127);
    assert_eq!(
        negative_charge.predicate(),
        &QueryNode::predicate(AtomQueryPredicate::FormalCharge(-129))
    );
    assert_eq!(negative_charge.formal_charge(), 127);
    assert!(!negative_charge.predicate_is_carrier_derived());
    assert_eq!(negative_charge.clone(), negative_charge);

    let mut negative_predicate = QueryAtom::from_identity_parts(
        AtomId::new(12),
        QueryAtomIdentity::AtomicNumber(6),
        QueryNode::predicate(AtomQueryPredicate::NegativeFormalCharge(129)),
    );
    negative_predicate.set_formal_charge(-128);
    assert_eq!(
        negative_predicate.predicate(),
        &QueryNode::predicate(AtomQueryPredicate::NegativeFormalCharge(129))
    );
    assert_eq!(negative_predicate.formal_charge(), -128);
    assert!(!negative_predicate.predicate_is_carrier_derived());
    assert_eq!(negative_predicate.clone(), negative_predicate);

    let mut isotope = QueryAtom::from_identity_parts(
        AtomId::new(13),
        QueryAtomIdentity::AtomicNumber(6),
        QueryNode::predicate(AtomQueryPredicate::Isotope(65_536)),
    );
    isotope.set_isotope(Some(0));
    assert_eq!(
        isotope.predicate(),
        &QueryNode::predicate(AtomQueryPredicate::Isotope(65_536))
    );
    assert_eq!(isotope.isotope(), None);
    assert!(!isotope.predicate_is_carrier_derived());
    assert_eq!(isotope.clone(), isotope);

    let ordinary = AtomSpec::new(Element::C)
        .with_formal_charge(-128)
        .with_isotope(u16::MAX);
    assert_eq!(ordinary.formal_charge(), -128);
    assert_eq!(ordinary.isotope(), Some(u16::MAX));
}

#[test]
fn seven_equality_query_targets_store_i32_without_changing_carriers() {
    let targets = [0, 255, 256, 2_147_483_639];
    let predicates: [fn(i32) -> AtomQueryPredicate; 7] = [
        AtomQueryPredicate::HydrogenCount,
        AtomQueryPredicate::ImplicitHydrogenCount,
        AtomQueryPredicate::ExplicitDegree,
        AtomQueryPredicate::NumHeteroatomNeighbors,
        AtomQueryPredicate::NumAliphaticHeteroatomNeighbors,
        AtomQueryPredicate::TotalDegree,
        AtomQueryPredicate::TotalValence,
    ];

    for (target_index, target) in targets.into_iter().enumerate() {
        for (predicate_index, make_predicate) in predicates.into_iter().enumerate() {
            let predicate = make_predicate(target);
            let mut atom = QueryAtom::from_identity_parts(
                AtomId::new(target_index * predicates.len() + predicate_index),
                QueryAtomIdentity::AtomicNumber(6),
                QueryNode::predicate(predicate.clone()),
            );
            atom.set_formal_charge(-7);
            atom.set_isotope(Some(13));
            atom.set_explicit_hydrogens(u8::MAX);
            atom.set_no_implicit(true);

            assert_eq!(atom.predicate(), &QueryNode::predicate(predicate));
            assert!(!atom.predicate_is_carrier_derived());
            assert_eq!(atom.formal_charge(), -7);
            assert_eq!(atom.isotope(), Some(13));
            assert_eq!(atom.explicit_hydrogens(), u8::MAX);
            assert!(atom.no_implicit());

            let cloned = atom.clone();
            assert_eq!(cloned, atom);
            assert_eq!(cloned.predicate(), atom.predicate());
            assert!(!cloned.predicate_is_carrier_derived());
            assert_eq!(cloned.formal_charge(), -7);
            assert_eq!(cloned.isotope(), Some(13));
            assert_eq!(cloned.explicit_hydrogens(), u8::MAX);
            assert!(cloned.no_implicit());
        }
    }
}

#[test]
fn query_graph_retains_raw_identity_and_rejects_concrete_topology_alignment() {
    let query_atom = QueryAtom::from_identity_parts(
        AtomId::new(0),
        QueryAtomIdentity::AtomicNumber(119),
        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(119)),
    );
    let graph = QueryGraph::from_parts(
        vec![query_atom],
        Vec::new(),
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .unwrap();
    let carrier = graph.atom(0).expect("raw identity remains a query row");
    assert_eq!(carrier.identity(), QueryAtomIdentity::AtomicNumber(119));
    assert_eq!(carrier.atomic_number(), 119);
    assert_eq!(
        carrier.predicate(),
        &QueryNode::predicate(AtomQueryPredicate::AtomicNumber(119))
    );

    let concrete_topology = TopologyBlock::try_from_parts(
        vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C))],
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .unwrap();
    assert!(matches!(
        QueryStateRef::try_for_topology(graph.atoms(), &[], &concrete_topology),
        Err(QueryStateError::NonElementAtomIdentity {
            position: 0,
            atom,
            atomic_number: 119,
        }) if atom == AtomId::new(0)
    ));
    assert_eq!(
        carrier.try_to_atom(),
        Err(QueryAtomConversionError::NonElementAtomicNumber {
            atom: AtomId::new(0),
            atomic_number: 119,
        })
    );
}

#[test]
fn query_atom_identity_canonicalizes_element_and_wildcard_numbers() {
    let carbon_predicate = QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6));
    let carbon = QueryAtom::from_identity_parts(
        AtomId::new(0),
        QueryAtomIdentity::AtomicNumber(6),
        carbon_predicate,
    );
    assert_eq!(carbon.identity(), QueryAtomIdentity::Element(Element::C));
    assert_eq!(carbon.atomic_number(), 6);
    assert_eq!(carbon.element(), Some(Element::C));
    assert_eq!(carbon.try_to_atom().unwrap().element(), Element::C);

    let wildcard = QueryAtom::from_identity_parts(
        AtomId::new(0),
        QueryAtomIdentity::from_atomic_number(0),
        QueryNode::predicate(AtomQueryPredicate::Any),
    );
    let dummy = Element::from_atomic_number(0).unwrap();
    assert_eq!(wildcard.identity(), QueryAtomIdentity::Element(dummy));
    assert_eq!(wildcard.atomic_number(), 0);
    assert_eq!(wildcard.element(), Some(dummy));
    assert_eq!(wildcard.try_to_atom().unwrap().element(), dummy);
}

#[test]
fn query_atom_identity_change_clone_mapping_and_common_mutation_preserve_carrier_state() {
    let id = AtomId::new(0);
    let attachment_order =
        TemplateAttachmentOrder::new(vec![TemplateAttachment::new(id, "self")]).unwrap();
    let spec = AtomSpec::new(Element::C)
        .with_formal_charge(1)
        .with_explicit_hydrogens(2)
        .with_chiral_tag(ChiralTag::TetrahedralCw)
        .with_chiral_permutation(3)
        .with_unknown_stereo(true)
        .with_mol_parity(7)
        .with_mol_inversion_flag(1)
        .with_implicit_hydrogen(true)
        .with_tracked_isotopic_hydrogens(vec![2, 3])
        .with_aromatic(true)
        .with_isotope(13)
        .with_atom_map(41)
        .with_no_implicit(true)
        .with_radical_electrons(1)
        .with_hybridization(Hybridization::Sp2)
        .with_prop("user", "kept")
        .unwrap()
        .with_computed_prop("computed", "kept")
        .unwrap()
        .with_pdb_residue_info(
            AtomPdbResidueInfo::new("CA", 12, "GLY", 3, "A", true)
                .with_alt_loc("B")
                .with_insertion_code("C")
                .with_occupancy(0.75)
                .with_temp_factor(12.5)
                .with_secondary_structure(4)
                .with_segment_number(5)
                .with_monomer_class("protein"),
        )
        .with_template_attachment_order(attachment_order.clone());
    let predicate = QueryNode::and(vec![
        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(119)),
        QueryNode::not(QueryNode::predicate(AtomQueryPredicate::IsAromatic(false))),
    ]);
    let element_carrier = Atom::from_spec(id, spec);
    let explicit = QueryAtom::from_parts(element_carrier.clone(), predicate.clone());
    let mut raw = explicit
        .clone()
        .with_identity(QueryAtomIdentity::AtomicNumber(119));

    assert_eq!(raw.identity(), QueryAtomIdentity::AtomicNumber(119));
    assert_eq!(raw.id(), id);
    assert_eq!(raw.formal_charge(), 1);
    assert_eq!(raw.explicit_hydrogens(), 2);
    assert_eq!(raw.chiral_tag(), ChiralTag::TetrahedralCw);
    assert_eq!(raw.chiral_permutation(), Some(3));
    assert!(raw.unknown_stereo());
    assert_eq!(raw.mol_parity(), Some(7));
    assert_eq!(raw.mol_inversion_flag(), Some(1));
    assert!(raw.implicit_hydrogen());
    assert_eq!(raw.tracked_isotopic_hydrogens(), &[2, 3]);
    assert!(raw.is_aromatic());
    assert_eq!(raw.isotope(), Some(13));
    assert_eq!(raw.atom_map(), Some(41));
    assert!(raw.no_implicit());
    assert_eq!(raw.radical_electrons(), 1);
    assert_eq!(raw.hybridization(), Hybridization::Sp2);
    assert_eq!(raw.prop("user"), Some("kept"));
    assert_eq!(raw.prop("computed"), Some("kept"));
    assert!(raw.is_prop_computed("computed"));
    assert_eq!(
        raw.computed_prop_names()
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        ["computed"]
    );
    assert_eq!(raw.pdb_residue_info(), element_carrier.pdb_residue_info());
    assert_eq!(raw.template_attachment_order(), Some(&attachment_order));
    assert_eq!(raw.predicate(), &predicate);
    assert!(!raw.predicate_is_carrier_derived());

    let cloned = raw.clone();
    assert_eq!(cloned, raw);
    let query_graph = QueryGraph::from_parts(
        vec![raw.clone()],
        Vec::new(),
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .unwrap();
    let cloned_graph = query_graph.clone();
    assert_eq!(cloned_graph, query_graph);
    assert_eq!(
        cloned_graph.atom(0).unwrap().identity(),
        QueryAtomIdentity::AtomicNumber(119)
    );

    raw.set_isotope(Some(14));
    raw.set_formal_charge(-2);
    raw.set_prop("added", "value").unwrap();
    assert_eq!(raw.identity(), QueryAtomIdentity::AtomicNumber(119));
    assert_eq!(raw.atomic_number(), 119);
    assert_eq!(raw.isotope(), Some(14));
    assert_eq!(raw.formal_charge(), -2);
    assert_eq!(raw.prop("added"), Some("value"));
    assert_eq!(raw.predicate(), &predicate);
    assert!(!raw.predicate_is_carrier_derived());

    let remapped_id = AtomId::new(1);
    let moved = raw.clone().with_id(remapped_id);
    assert_eq!(moved.id(), remapped_id);
    assert_eq!(moved.identity(), QueryAtomIdentity::AtomicNumber(119));
    assert_eq!(moved.isotope(), Some(14));
    assert_eq!(moved.formal_charge(), -2);
    assert_eq!(moved.predicate(), &predicate);
}
