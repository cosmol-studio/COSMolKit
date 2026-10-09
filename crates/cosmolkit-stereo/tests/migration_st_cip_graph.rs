#![allow(dead_code)]

include!("../src/cip_graph.rs");

#[cfg(test)]
mod migration_tests {
    use std::{cell::Cell, rc::Rc};

    use cosmolkit_model::{AtomId, AtomSpec, BondSpec};

    use super::*;

    fn topology(atom_specs: Vec<AtomSpec>, bond_specs: Vec<BondSpec>) -> TopologyBlock {
        let atoms = atom_specs
            .into_iter()
            .enumerate()
            .map(|(index, spec)| Atom::from_spec(AtomId::new(index), spec))
            .collect();
        let bonds = bond_specs
            .into_iter()
            .enumerate()
            .map(|(index, spec)| Bond::from_spec(BondId::new(index), spec))
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
    }

    fn graph(elements: &[Element], edges: &[(usize, usize, BondOrder)]) -> TopologyBlock {
        topology(
            elements.iter().copied().map(AtomSpec::new).collect(),
            edges
                .iter()
                .map(|(begin, end, order)| {
                    BondSpec::new(AtomId::new(*begin), AtomId::new(*end), *order)
                })
                .collect(),
        )
    }

    fn aromatic_ring(elements: &[Element]) -> TopologyBlock {
        let count = elements.len();
        topology(
            elements
                .iter()
                .copied()
                .map(|element| AtomSpec::new(element).with_aromatic(true))
                .collect(),
            (0..count)
                .map(|begin| {
                    BondSpec::new(
                        AtomId::new(begin),
                        AtomId::new((begin + 1) % count),
                        BondOrder::Aromatic,
                    )
                    .with_aromatic(true)
                })
                .collect(),
        )
    }

    fn edge_for_atom(digraph: &CipDigraph<'_>, edges: &[CipEdgeId], atom: usize) -> CipEdgeId {
        edges
            .iter()
            .copied()
            .find(|edge| digraph.node(digraph.edge(*edge).get_end()).atom_idx() == Some(atom))
            .unwrap()
    }

    fn edge_end_atomic_numbers(digraph: &CipDigraph<'_>, edges: &[CipEdgeId]) -> Vec<u8> {
        edges
            .iter()
            .map(|edge| {
                digraph
                    .node(digraph.edge(*edge).get_end())
                    .get_atomic_num(digraph.mol())
                    .unwrap()
            })
            .collect()
    }

    struct AtomicNumberRule;

    impl CipSequenceRule for AtomicNumberRule {
        fn compare(
            &self,
            digraph: &mut CipDigraph<'_>,
            _context: &mut CipLabelerContext,
            left: CipEdgeId,
            right: CipEdgeId,
        ) -> Result<i32, CipLabelerError> {
            let left = digraph
                .node(digraph.edge(left).get_end())
                .get_atomic_num(digraph.mol())?;
            let right = digraph
                .node(digraph.edge(right).get_end())
                .get_atomic_num(digraph.mol())?;
            Ok(three_way_comparison_i32(i32::from(left), i32::from(right)))
        }
    }

    struct EqualRule;

    impl CipSequenceRule for EqualRule {
        fn compare(
            &self,
            _digraph: &mut CipDigraph<'_>,
            _context: &mut CipLabelerContext,
            _left: CipEdgeId,
            _right: CipEdgeId,
        ) -> Result<i32, CipLabelerError> {
            Ok(0)
        }
    }

    struct PseudoAsymmetricRule;

    impl CipSequenceRule for PseudoAsymmetricRule {
        fn compare(
            &self,
            digraph: &mut CipDigraph<'_>,
            context: &mut CipLabelerContext,
            left: CipEdgeId,
            right: CipEdgeId,
        ) -> Result<i32, CipLabelerError> {
            Ok(AtomicNumberRule.compare(digraph, context, left, right)? * 2)
        }
    }

    struct CountingRule {
        calls: Rc<Cell<usize>>,
    }

    impl CipSequenceRule for CountingRule {
        fn compare(
            &self,
            digraph: &mut CipDigraph<'_>,
            context: &mut CipLabelerContext,
            left: CipEdgeId,
            right: CipEdgeId,
        ) -> Result<i32, CipLabelerError> {
            self.calls.set(self.calls.get() + 1);
            AtomicNumberRule.compare(digraph, context, left, right)
        }
    }

    #[test]
    fn ciplabeler_sort_priority_descriptor_spelling_is_source_ordered() {
        let expected = [
            (Descriptor::None, "NONE"),
            (Descriptor::Unknown, "UNKNOWN"),
            (Descriptor::ns, "ns"),
            (Descriptor::R, "R"),
            (Descriptor::S, "S"),
            (Descriptor::r, "r"),
            (Descriptor::s, "s"),
            (Descriptor::seqTrans, "e"),
            (Descriptor::seqCis, "z"),
            (Descriptor::E, "E"),
            (Descriptor::Z, "Z"),
            (Descriptor::M, "M"),
            (Descriptor::P, "P"),
            (Descriptor::m, "m"),
            (Descriptor::p, "p"),
            (Descriptor::SP_4, "SP_4"),
            (Descriptor::TBPY_5, "TBPY_5"),
            (Descriptor::OC_6, "OC_6"),
        ];
        for (descriptor, spelling) in expected {
            assert_eq!(descriptor_to_string(descriptor), spelling);
        }
        assert!(
            Descriptor::ALL_IN_RDKIT_ORDER
                .windows(2)
                .all(|pair| pair[0] < pair[1])
        );
        assert_eq!(
            public_descriptor_from_source(Descriptor::R),
            Some(CipDescriptor::R)
        );
        assert_eq!(
            public_descriptor_from_source(Descriptor::seqTrans),
            Some(CipDescriptor::LowerE)
        );
        assert_eq!(
            public_descriptor_from_source(Descriptor::seqCis),
            Some(CipDescriptor::LowerZ)
        );
        assert_eq!(
            public_descriptor_from_source(Descriptor::seqTrans)
                .map(|descriptor| descriptor.as_str()),
            Some("e")
        );
        assert_eq!(
            public_descriptor_from_source(Descriptor::seqCis).map(|descriptor| descriptor.as_str()),
            Some("z")
        );
        for descriptor in [
            Descriptor::None,
            Descriptor::Unknown,
            Descriptor::ns,
            Descriptor::SP_4,
            Descriptor::TBPY_5,
            Descriptor::OC_6,
        ] {
            assert_eq!(public_descriptor_from_source(descriptor), None);
        }
    }

    #[test]
    fn ciplabeler_exact_width_rational_mass_and_counter_boundaries_match_source() {
        assert_eq!(CIP_NO_ATOM, u32::MAX);
        assert_eq!(CipNode::NO_ATOM_INDEX, u32::MAX);
        assert!(RationalI32::new(13, 2) < RationalI32::new(7, 1));
        assert_eq!(RationalI32::new(6, 4), RationalI32::new(3, 2));
        assert_eq!(rdkit_atomic_mass(6, Some(999)).unwrap(), 999.0);
        assert_eq!(rdkit_atomic_mass(0, Some(999)).unwrap(), 0.0);
        assert!(rdkit_atomic_mass(6, None).unwrap() > 12.0);
        assert!(matches!(
            rdkit_atomic_mass(u8::MAX, None),
            Err(CipLabelerError::InvalidInternalState { .. })
        ));

        let mut wrapped = CipLabelerContext::with_remaining_call_count(0);
        assert!(wrapped.decrement_remaining_call_count_and_check());
        let mut one = CipLabelerContext::with_remaining_call_count(1);
        assert!(!one.decrement_remaining_call_count_and_check());
        let mut two = CipLabelerContext::with_remaining_call_count(2);
        assert!(two.decrement_remaining_call_count_and_check());
        assert!(!two.decrement_remaining_call_count_and_check());
    }

    #[test]
    fn ciplabeler_cipmol_access_ring_order_mancude_and_error_paths_match_source() {
        let chain = graph(
            &[Element::C, Element::C, Element::O],
            &[(0, 1, BondOrder::Single), (1, 2, BondOrder::Single)],
        );
        let mut view = CipMol::new(&chain);
        assert_eq!(view.get_num_atoms(), 3);
        assert_eq!(view.get_num_bonds(), 2);
        assert_eq!(view.neighbor_indices(1).unwrap(), vec![0, 2]);
        assert_eq!(view.bond_indices_for_atom(1).unwrap(), vec![0, 1]);
        assert!(!view.is_in_ring(0).unwrap());

        let cycle = graph(
            &[Element::C, Element::C, Element::C],
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 0, BondOrder::Single),
            ],
        );
        assert!(CipMol::new(&cycle).is_in_ring(2).unwrap());

        let benzene = aromatic_ring(&[Element::C; 6]);
        let mut benzene_view = CipMol::new(&benzene);
        let mut orders = (0..6)
            .map(|bond| benzene_view.get_bond_order(bond).unwrap())
            .collect::<Vec<_>>();
        orders.sort_unstable();
        assert_eq!(orders, vec![1, 1, 1, 2, 2, 2]);

        let pyridine = aromatic_ring(&[
            Element::C,
            Element::C,
            Element::C,
            Element::N,
            Element::C,
            Element::C,
        ]);
        let mut pyridine_view = CipMol::new(&pyridine);
        let fractions = (0..6)
            .map(|atom| {
                pyridine_view
                    .get_fractional_atomic_num(atom)
                    .unwrap()
                    .tuple()
            })
            .collect::<Vec<_>>();
        assert_eq!(
            fractions,
            vec![(6, 1), (6, 1), (13, 2), (6, 1), (13, 2), (6, 1)]
        );

        let non_integer = graph(&[Element::C, Element::C], &[(0, 1, BondOrder::OneAndHalf)]);
        assert!(matches!(
            CipMol::new(&non_integer).get_bond_order(0),
            Err(CipLabelerError::NonIntegerBondOrder {
                order: BondOrder::OneAndHalf
            })
        ));
        assert!(matches!(
            view.atom(3),
            Err(CipLabelerError::AtomIndexOutOfRange {
                index: 3,
                atom_count: 3
            })
        ));
        assert!(matches!(
            view.bond(2),
            Err(CipLabelerError::BondIndexOutOfRange {
                index: 2,
                bond_count: 2
            })
        ));
        assert!(matches!(
            view.other_atom_idx(0, 2),
            Err(CipLabelerError::BondNotIncident { bond: 0, atom: 2 })
        ));
    }

    #[test]
    fn ciplabeler_node_and_edge_state_duplicate_and_implicit_h_paths_match_source() {
        let isotope = topology(
            vec![
                AtomSpec::new(Element::C)
                    .with_isotope(13)
                    .with_explicit_hydrogens(4)
                    .with_no_implicit(true),
            ],
            Vec::new(),
        );
        let mut view = CipMol::new(&isotope);
        let root = CipNode::new(11, vec![1], Some(0), RationalI32::new(6, 1), 1, 0, &view).unwrap();
        assert_eq!(root.get_digraph(), 11);
        assert_eq!(root.get_atom_idx().unwrap(), 0);
        assert_eq!(root.get_mass_num(&view).unwrap(), 13);
        assert!(root.get_atomic_mass() > 13.0);
        assert!(root.is_visited(0));
        assert!(!root.is_expanded());

        let duplicate = root
            .new_bond_duplicate_child(0, Some(0), &mut view)
            .unwrap();
        assert!(duplicate.is_duplicate());
        assert!(duplicate.is_terminal());
        assert_eq!(duplicate.get_atomic_mass(), 0.0);

        let implicit = root.new_implicit_hydrogen_child(&mut view).unwrap();
        assert_eq!(implicit.atom_idx(), None);
        assert_eq!(implicit.get_atom_idx().unwrap(), u32::MAX);
        assert_eq!(implicit.get_atomic_num(&view).unwrap(), 1);
        assert!(implicit.is_duplicate_or_h());

        let mut edge = CipEdge::new(CipNodeId::new(0), CipNodeId::new(1), Some(7));
        assert_eq!(
            edge.get_other(CipEdgeId::new(4), CipNodeId::new(0))
                .unwrap(),
            CipNodeId::new(1)
        );
        assert!(matches!(
            edge.get_other(CipEdgeId::new(4), CipNodeId::new(2)),
            Err(CipLabelerError::EdgeEndpointMismatch { edge: 4, node: 2 })
        ));
        edge.set_aux(Descriptor::seqTrans);
        edge.flip();
        assert_eq!(edge.get_beg(), CipNodeId::new(1));
        assert_eq!(edge.get_end(), CipNodeId::new(0));
        assert_eq!(edge.get_aux(), Descriptor::seqTrans);
    }

    #[test]
    fn ciplabeler_digraph_expansion_order_duplicates_roots_and_cap_match_source() {
        let ethane = graph(&[Element::C, Element::C], &[(0, 1, BondOrder::Single)]);
        let mut digraph = CipDigraph::new(&ethane, 0, false).unwrap();
        let root = digraph.get_current_root();
        let edges = digraph.node_edges(root).unwrap();
        assert_eq!(edges.len(), 4);
        assert_eq!(digraph.get_num_nodes(), 5);
        assert_eq!(digraph.node_edges(root).unwrap(), edges);
        assert_eq!(
            edges
                .iter()
                .filter(|edge| digraph.edge(**edge).get_bond_idx().is_none())
                .count(),
            3
        );

        let ethene = graph(&[Element::C, Element::C], &[(0, 1, BondOrder::Double)]);
        let mut ordinary = CipDigraph::new(&ethene, 0, false).unwrap();
        let ordinary_root = ordinary.get_current_root();
        assert!(
            !ordinary
                .node_edges(ordinary_root)
                .unwrap()
                .iter()
                .any(|edge| { ordinary.node(ordinary.edge(*edge).get_end()).is_duplicate() })
        );
        let mut atrop = CipDigraph::new(&ethene, 0, true).unwrap();
        let atrop_root = atrop.get_current_root();
        assert_eq!(
            atrop
                .node_edges(atrop_root)
                .unwrap()
                .iter()
                .filter(|edge| atrop.node(atrop.edge(**edge).get_end()).is_duplicate())
                .count(),
            1
        );

        let ring = graph(
            &[Element::C, Element::C, Element::C],
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 0, BondOrder::Single),
            ],
        );
        let mut ring_graph = CipDigraph::new(&ring, 0, false).unwrap();
        let root = ring_graph.get_current_root();
        let root_edges = ring_graph.node_edges(root).unwrap();
        let first = edge_for_atom(&ring_graph, &root_edges, 1);
        let one = ring_graph.edge(first).get_end();
        let one_edges = ring_graph.node_edges(one).unwrap();
        let second = edge_for_atom(&ring_graph, &one_edges, 2);
        let two = ring_graph.edge(second).get_end();
        assert!(ring_graph.node_edges(two).unwrap().iter().any(|edge| {
            ring_graph
                .node(ring_graph.edge(*edge).get_end())
                .is_set(CipNode::RING_DUPLICATE)
        }));
        let repeated = ring_graph.get_nodes(0).unwrap();
        assert!(repeated.len() > 1);
        assert!(repeated.windows(2).all(|pair| {
            ring_graph.node(pair[0]).get_distance() <= ring_graph.node(pair[1]).get_distance()
        }));

        ring_graph.change_root(one).unwrap();
        assert_eq!(ring_graph.get_current_root(), one);
        assert_eq!(ring_graph.get_original_root(), root);

        let mut capped = CipDigraph::new(&ethane, 0, false).unwrap();
        let cap_root = capped.get_current_root();
        let seed = capped.nodes[cap_root.index()].clone();
        capped.nodes.resize(CipDigraph::MAX_NODE_COUNT, seed);
        let error = capped.expand(cap_root).unwrap_err();
        assert_eq!(
            error.to_string(),
            "Digraph generation failed: more than 100000 nodes found."
        );
    }

    #[test]
    fn ciplabeler_sort_priority_unique_tied_pseudo_and_rule_order_match_source() {
        let substituted = graph(
            &[Element::C, Element::F, Element::CL, Element::BR],
            &[
                (0, 1, BondOrder::Single),
                (0, 2, BondOrder::Single),
                (0, 3, BondOrder::Single),
            ],
        );
        let mut digraph = CipDigraph::new(&substituted, 0, false).unwrap();
        let root = digraph.get_current_root();
        let mut edges = digraph.node_edges(root).unwrap();
        let sorter = CipSort::new(&AtomicNumberRule);
        let mut context = CipLabelerContext::new(0);
        let priority = sorter
            .prioritize(&mut digraph, &mut context, root, &mut edges, true)
            .unwrap();
        assert!(priority.is_unique());
        assert_eq!(
            edge_end_atomic_numbers(&digraph, &edges),
            vec![35, 17, 9, 1]
        );

        let tied = graph(
            &[Element::C, Element::C, Element::C, Element::N],
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (1, 3, BondOrder::Single),
            ],
        );
        let mut tied_graph = CipDigraph::new(&tied, 1, false).unwrap();
        let root = tied_graph.get_current_root();
        let mut edges = tied_graph.node_edges(root).unwrap();
        let mut context = CipLabelerContext::new(0);
        let priority = sorter
            .prioritize(&mut tied_graph, &mut context, root, &mut edges, true)
            .unwrap();
        assert!(!priority.is_unique());
        assert_eq!(
            sorter
                .get_groups(&mut tied_graph, &mut context, &edges)
                .unwrap()
                .iter()
                .map(Vec::len)
                .collect::<Vec<_>>(),
            vec![1, 2, 1]
        );

        let calls = Rc::new(Cell::new(0));
        let rules = CipRules::new(vec![
            Box::new(EqualRule),
            Box::new(CountingRule {
                calls: Rc::clone(&calls),
            }),
        ])
        .unwrap();
        let refs = rules.rule_refs();
        let carbon_edges = edges
            .iter()
            .copied()
            .filter(|edge| {
                tied_graph
                    .node(tied_graph.edge(*edge).get_end())
                    .get_atomic_num(tied_graph.mol())
                    .unwrap()
                    == 6
            })
            .collect::<Vec<_>>();
        let mut context = CipLabelerContext::new(0);
        refs[0]
            .recursive_compare_with_sort_rules(
                &refs,
                &mut tied_graph,
                &mut context,
                carbon_edges[0],
                carbon_edges[1],
            )
            .unwrap();
        assert!(calls.get() > 0);
        assert_eq!(cip_constitutional_rules().unwrap().get_num_sub_rules(), 3);
        assert_eq!(cip_all_rules().unwrap().get_num_sub_rules(), 9);

        let f_o_n = graph(
            &[Element::F, Element::O, Element::N],
            &[(0, 1, BondOrder::Single), (1, 2, BondOrder::Single)],
        );
        let mut pseudo_graph = CipDigraph::new(&f_o_n, 1, false).unwrap();
        let pseudo_root = pseudo_graph.get_current_root();
        let mut pseudo_edges = pseudo_graph.node_edges(pseudo_root).unwrap();
        let pseudo_sorter = CipSort::new(&PseudoAsymmetricRule);
        let mut context = CipLabelerContext::new(0);
        let priority = pseudo_sorter
            .prioritize(
                &mut pseudo_graph,
                &mut context,
                pseudo_root,
                &mut pseudo_edges,
                true,
            )
            .unwrap();
        assert!(priority.is_unique());
        assert!(priority.is_pseudo_asymetric());
        assert_eq!(
            edge_end_atomic_numbers(&pseudo_graph, &pseudo_edges),
            vec![9, 7]
        );
    }

    #[test]
    fn ciplabeler_rules_constructor_prefix_and_missing_rule_errors_match_source() {
        let rules = CipRules::new(vec![Box::new(EqualRule), Box::new(AtomicNumberRule)]).unwrap();
        assert_eq!(rules.get_num_sub_rules(), 2);
        assert_eq!(rules.get_sorter().get_rules().len(), 1);

        let mut empty = CipRules::new(Vec::new()).unwrap();
        assert!(matches!(
            empty.add(None),
            Err(CipLabelerError::NoSequenceRuleProvided)
        ));

        let mut empty_pairs = CipPairList::new();
        assert_eq!(empty_pairs.to_rdkit_string(), "");
        assert!(!empty_pairs.add(Descriptor::None));
        let mut tied_left = CipPairList::with_ref(Descriptor::R);
        let mut tied_right = CipPairList::with_ref(Descriptor::S);
        tied_left.add(Descriptor::R);
        tied_right.add(Descriptor::S);
        assert_eq!(tied_left.compare_to(&tied_right).unwrap(), 0);
        tied_right.add(Descriptor::S);
        assert!(matches!(
            tied_left.compare_to(&tied_right),
            Err(CipLabelerError::DescriptorListLengthMismatch)
        ));
    }

    #[test]
    fn ciplabeler_sequence_rule_deep_budget_and_up_edge_errors_match_source() {
        let branches = graph(
            &[Element::C, Element::C, Element::F, Element::C, Element::CL],
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (0, 3, BondOrder::Single),
                (3, 4, BondOrder::Single),
            ],
        );
        let mut digraph = CipDigraph::new(&branches, 0, false).unwrap();
        let root = digraph.get_current_root();
        let edges = digraph.node_edges(root).unwrap();
        let carbons = edges
            .iter()
            .copied()
            .filter(|edge| {
                digraph
                    .node(digraph.edge(*edge).get_end())
                    .get_atomic_num(digraph.mol())
                    .unwrap()
                    == 6
            })
            .collect::<Vec<_>>();
        let mut shallow = CipLabelerContext::new(0);
        assert_eq!(
            AtomicNumberRule
                .get_comparison(&mut digraph, &mut shallow, carbons[0], carbons[1], false)
                .unwrap(),
            0
        );
        let mut deep = CipLabelerContext::new(0);
        assert_eq!(
            AtomicNumberRule
                .get_comparison(&mut digraph, &mut deep, carbons[0], carbons[1], true)
                .unwrap(),
            -1
        );
        let mut exhausted = CipLabelerContext::with_remaining_call_count(1);
        assert!(matches!(
            AtomicNumberRule.recursive_compare(
                &mut digraph,
                &mut exhausted,
                carbons[0],
                carbons[1]
            ),
            Err(CipLabelerError::MaxIterationsExceeded)
        ));

        let chain = graph(
            &[Element::C, Element::C, Element::C],
            &[(0, 1, BondOrder::Single), (1, 2, BondOrder::Single)],
        );
        let mut chain_graph = CipDigraph::new(&chain, 1, false).unwrap();
        let center = chain_graph.get_current_root();
        let center_edges = chain_graph.node_edges(center).unwrap();
        let left_edge = edge_for_atom(&chain_graph, &center_edges, 0);
        let right_edge = edge_for_atom(&chain_graph, &center_edges, 2);
        let left = chain_graph.edge(left_edge).get_end();
        chain_graph.node_edges(left).unwrap();
        assert!(
            AtomicNumberRule
                .are_up_edges(&chain_graph, left, left, left_edge, left_edge)
                .unwrap()
        );
        assert!(matches!(
            AtomicNumberRule.are_up_edges(&chain_graph, left, center, left_edge, right_edge),
            Err(CipLabelerError::UnexpectedUpEdgeOrdering)
        ));
    }

    #[test]
    fn ciplabeler_rules_rule1_rule2_isotope_duplicate_and_dummy_paths_match_source() {
        let isotopes = topology(
            vec![
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::C).with_isotope(12),
                AtomSpec::new(Element::C).with_isotope(13),
            ],
            vec![
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
                BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single),
            ],
        );
        let mut digraph = CipDigraph::new(&isotopes, 0, false).unwrap();
        let root = digraph.get_current_root();
        let edges = digraph.node_edges(root).unwrap();
        let c12 = edge_for_atom(&digraph, &edges, 1);
        let c13 = edge_for_atom(&digraph, &edges, 2);
        let mut context = CipLabelerContext::new(0);
        assert_eq!(
            CipRule2
                .compare(&mut digraph, &mut context, c12, c13)
                .unwrap(),
            -1
        );

        let normal = digraph
            .add_node(Vec::new(), Some(1), RationalI32::new(6, 1), 2, 0)
            .unwrap();
        let fraction = digraph
            .add_node(Vec::new(), Some(1), RationalI32::new(13, 2), 2, 0)
            .unwrap();
        digraph.add_edge(root, Some(0), normal);
        let normal_edge = CipEdgeId::new(digraph.edges.len() - 1);
        digraph.add_edge(root, Some(0), fraction);
        let fraction_edge = CipEdgeId::new(digraph.edges.len() - 1);
        assert_eq!(
            CipRule1a
                .compare(&mut digraph, &mut context, normal_edge, fraction_edge)
                .unwrap(),
            -1
        );

        let ring_duplicate = digraph
            .add_node(
                Vec::new(),
                Some(1),
                RationalI32::new(6, 1),
                3,
                CipNode::RING_DUPLICATE,
            )
            .unwrap();
        digraph.add_edge(root, Some(0), ring_duplicate);
        let duplicate_edge = CipEdgeId::new(digraph.edges.len() - 1);
        assert_eq!(
            CipRule1b
                .compare(&mut digraph, &mut context, duplicate_edge, normal_edge)
                .unwrap(),
            1
        );

        let dummy = topology(
            vec![
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::DUMMY),
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::N),
            ],
            vec![
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
                BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single),
                BondSpec::new(AtomId::new(0), AtomId::new(3), BondOrder::Single),
            ],
        );
        let mut dummy_graph = CipDigraph::new(&dummy, 0, false).unwrap();
        let root = dummy_graph.get_current_root();
        let edges = dummy_graph.node_edges(root).unwrap();
        let dummy_edge = edge_for_atom(&dummy_graph, &edges, 1);
        let carbon_edge = edge_for_atom(&dummy_graph, &edges, 2);
        let nitrogen_edge = edge_for_atom(&dummy_graph, &edges, 3);
        assert_eq!(
            CipRule2
                .compare(&mut dummy_graph, &mut context, dummy_edge, carbon_edge)
                .unwrap(),
            -1
        );
        assert_eq!(
            CipRule2
                .compare(&mut dummy_graph, &mut context, carbon_edge, nitrogen_edge)
                .unwrap(),
            0
        );
    }

    #[test]
    fn ciplabeler_pairlist_reference_pairing_order_and_width_match_source() {
        assert_eq!(CipPairList::ref_descriptor(Descriptor::M), Descriptor::R);
        assert_eq!(CipPairList::ref_descriptor(Descriptor::P), Descriptor::S);
        assert_eq!(CipPairList::ref_descriptor(Descriptor::r), Descriptor::None);

        let mut pairs = CipPairList::with_ref(Descriptor::R);
        for descriptor in [Descriptor::R, Descriptor::P, Descriptor::R] {
            assert!(pairs.add(descriptor));
        }
        assert_eq!(pairs.get_ref_descriptor(), Ok(Descriptor::R));
        assert_eq!(pairs.to_rdkit_string(), "R:lul");

        let mut head = CipPairList::with_ref(Descriptor::R);
        head.add(Descriptor::S);
        let mut tail = CipPairList::new();
        tail.add_all(&[Descriptor::None, Descriptor::seqCis, Descriptor::Unknown]);
        let combined = CipPairList::from_head_tail(&head, &tail);
        assert_eq!(combined.to_rdkit_string(), "R:ul");

        let mut like = CipPairList::with_ref(Descriptor::R);
        like.add(Descriptor::R);
        let mut unlike = CipPairList::with_ref(Descriptor::R);
        unlike.add(Descriptor::S);
        assert_eq!(like.compare_to(&unlike).unwrap(), 1);
        assert!(unlike.less_than(&like).unwrap());
        let mut lists = [unlike.clone(), like.clone()];
        CipPairList::sort_descending(&mut lists).unwrap();
        assert_eq!(lists, [like, unlike]);

        let mut fixed_width = CipPairList::with_ref(Descriptor::R);
        fixed_width.add(Descriptor::R);
        for _ in 0..61 {
            fixed_width.add(Descriptor::S);
        }
        fixed_width.add(Descriptor::R);
        assert_eq!(
            fixed_width.to_rdkit_string(),
            format!("R:l{}l", "u".repeat(61))
        );
        fixed_width.add_all(&[Descriptor::R; 66]);
        assert_eq!(fixed_width.descriptors.len(), 130);
        assert_eq!(
            fixed_width.to_rdkit_string(),
            format!("R:l{}{}", "u".repeat(61), "l".repeat(67))
        );
    }

    #[test]
    fn ciplabeler_rules_rule3_rule4_ordinals_and_invalid_descriptor_match_source() {
        assert_eq!(CipRule3::ord(Descriptor::E), 1);
        assert_eq!(CipRule3::ord(Descriptor::Z), 2);
        assert_eq!(CipRule3::ord(Descriptor::R), 0);
        for descriptor in [Descriptor::Unknown, Descriptor::ns, Descriptor::None] {
            assert_eq!(CipRule4a::ord(descriptor).unwrap(), 0);
        }
        for descriptor in [Descriptor::r, Descriptor::s, Descriptor::m, Descriptor::p] {
            assert_eq!(CipRule4a::ord(descriptor).unwrap(), 1);
        }
        for descriptor in [Descriptor::R, Descriptor::S, Descriptor::M, Descriptor::P] {
            assert_eq!(CipRule4a::ord(descriptor).unwrap(), 2);
        }
        assert!(matches!(
            CipRule4a::ord(Descriptor::SP_4),
            Err(CipLabelerError::InvalidStereoDescriptor)
        ));
        assert_eq!(CipRule4c::ord(Descriptor::m), 2);
        assert_eq!(CipRule4c::ord(Descriptor::p), 1);
        assert_eq!(CipRule4c::ord(Descriptor::R), 0);
    }

    #[test]
    fn ciplabeler_rules_rule4b_rule5new_and_rule6_reference_paths_match_source() {
        let chain = graph(
            &[Element::C, Element::C, Element::C],
            &[(0, 1, BondOrder::Single), (1, 2, BondOrder::Single)],
        );
        let mut digraph = CipDigraph::new(&chain, 1, false).unwrap();
        let root = digraph.get_current_root();
        let branch = digraph
            .add_node(vec![1, 1, 0], Some(0), RationalI32::new(6, 1), 2, 0)
            .unwrap();
        digraph.add_edge(root, Some(0), branch);
        let r_node = digraph
            .add_node(Vec::new(), Some(0), RationalI32::new(6, 1), 3, 0)
            .unwrap();
        let s_node = digraph
            .add_node(Vec::new(), Some(0), RationalI32::new(6, 1), 3, 0)
            .unwrap();
        digraph.nodes[r_node.index()].set_aux(Descriptor::R);
        digraph.nodes[s_node.index()].set_aux(Descriptor::S);
        digraph.add_edge(branch, Some(0), r_node);
        let r_edge = CipEdgeId::new(digraph.edges.len() - 1);
        digraph.add_edge(branch, Some(0), s_node);
        let s_edge = CipEdgeId::new(digraph.edges.len() - 1);
        let mut context = CipLabelerContext::new(0);
        assert_eq!(
            CipRule4b::with_ref(Descriptor::R)
                .compare(&mut digraph, &mut context, r_edge, s_edge)
                .unwrap(),
            1
        );
        assert_eq!(
            CipRule5New::with_ref(Descriptor::S)
                .compare(&mut digraph, &mut context, r_edge, s_edge)
                .unwrap(),
            -1
        );
        assert_eq!(
            CipRule4b::new()
                .compare(&mut digraph, &mut context, r_edge, s_edge)
                .unwrap(),
            0
        );

        let real_edges = digraph.node_edges(root).unwrap();
        let left = edge_for_atom(&digraph, &real_edges, 0);
        let right = edge_for_atom(&digraph, &real_edges, 2);
        assert_eq!(
            CipRule6
                .compare(&mut digraph, &mut context, left, right)
                .unwrap(),
            0
        );
        digraph.set_rule6_ref(Some(0)).unwrap();
        assert_eq!(
            CipRule6
                .compare(&mut digraph, &mut context, left, right)
                .unwrap(),
            2
        );
        digraph.set_rule6_ref(Some(2)).unwrap();
        assert_eq!(
            CipRule6
                .compare(&mut digraph, &mut context, left, right)
                .unwrap(),
            -2
        );
        assert!(matches!(
            digraph.set_rule6_ref(Some(3)),
            Err(CipLabelerError::AtomIndexOutOfRange {
                index: 3,
                atom_count: 3
            })
        ));
    }

    #[test]
    fn ciplabeler_rules_rule4b_rule5new_root_comparisons_restore_root() {
        let chain = graph(
            &[Element::C, Element::C, Element::C, Element::C],
            &[
                (0, 1, BondOrder::Single),
                (1, 2, BondOrder::Single),
                (2, 3, BondOrder::Single),
            ],
        );
        let mut digraph = CipDigraph::new(&chain, 1, false).unwrap();
        let root = digraph.get_current_root();
        let root_edges = digraph.node_edges(root).unwrap();
        let left = edge_for_atom(&digraph, &root_edges, 0);
        let right = edge_for_atom(&digraph, &root_edges, 2);
        let left_node = digraph.edge(left).get_end();
        let right_node = digraph.edge(right).get_end();
        digraph.nodes[left_node.index()].set_aux(Descriptor::R);
        digraph.nodes[right_node.index()].set_aux(Descriptor::S);

        let rule4b = CipRule4b::new();
        let rule4b_refs = [&rule4b as &dyn CipSequenceRule];
        let mut context = CipLabelerContext::new(0);
        assert_eq!(
            rule4b
                .compare_with_sort_rules(
                    Some(&rule4b_refs),
                    &mut digraph,
                    &mut context,
                    left,
                    right,
                )
                .unwrap(),
            0
        );
        assert_eq!(digraph.get_current_root(), root);

        let rule5 = CipRule5New::new();
        let rule5_refs = [&rule5 as &dyn CipSequenceRule];
        assert_eq!(
            rule5
                .compare_with_sort_rules(
                    Some(&rule5_refs),
                    &mut digraph,
                    &mut context,
                    left,
                    right,
                )
                .unwrap(),
            2
        );
        assert_eq!(digraph.get_current_root(), root);
        let unrelated = CipRule1a;
        assert!(matches!(
            rule5.get_ref_sorter(
                Some(&[&unrelated as &dyn CipSequenceRule]),
                &CipRule5New::with_ref(Descriptor::R),
            ),
            Err(CipLabelerError::Rule5NewInstanceNotInRuleSet)
        ));
    }

    #[test]
    fn ciplabeler_exact_width_real_atom_index_conversion_is_fallible() {
        let accepted = CipNode {
            digraph: 0,
            atom_idx: Some(u32::MAX as usize),
            distance: 0,
            atomic_num_fraction: RationalI32::new(6, 1),
            atomic_mass: 12.0,
            aux: Descriptor::None,
            flags: 0,
            edges: Vec::new(),
            visit: Vec::new(),
        };
        assert_eq!(accepted.get_atom_idx().unwrap(), u32::MAX);

        #[cfg(target_pointer_width = "64")]
        {
            let index = u32::MAX as usize + 1;
            let rejected = CipNode {
                atom_idx: Some(index),
                ..accepted
            };
            assert_eq!(
                rejected.get_atom_idx(),
                Err(CipLabelerError::SourceIndexWidthExceeded {
                    kind: "atom",
                    index,
                })
            );
        }
    }

    #[test]
    fn ciplabeler_mancude_seed_relax_partition_and_nonresonant_removal_match_source() {
        let pyridine = aromatic_ring(&[
            Element::C,
            Element::C,
            Element::C,
            Element::N,
            Element::C,
            Element::C,
        ]);
        let mut view = CipMol::new(&pyridine);
        let mut types = vec![MancudeType::Other; 6];
        assert!(seed_types(&mut types, &mut view).unwrap());
        assert_eq!(types[3], MancudeType::Nv3D2);
        relax_types(&mut types, &view).unwrap();
        assert!(types.iter().all(|kind| *kind != MancudeType::Other));
        let mut parts = vec![0; 6];
        assert_eq!(visit_parts(&mut parts, &types, &mut view).unwrap(), 1);
        assert_eq!(parts, vec![1; 6]);

        let propane = graph(
            &[Element::C, Element::C, Element::C],
            &[(0, 1, BondOrder::Single), (1, 2, BondOrder::Single)],
        );
        let view = CipMol::new(&propane);
        let mut chain_types = vec![MancudeType::Cv4D3, MancudeType::Cv4D3, MancudeType::Other];
        relax_types(&mut chain_types, &view).unwrap();
        assert_eq!(chain_types, vec![MancudeType::Other; 3]);

        let charged_ring = topology(
            vec![
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::C).with_formal_charge(-1),
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::C),
            ],
            vec![
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Double),
                BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Single),
                BondSpec::new(AtomId::new(2), AtomId::new(3), BondOrder::Single),
                BondSpec::new(AtomId::new(3), AtomId::new(4), BondOrder::Double),
                BondSpec::new(AtomId::new(4), AtomId::new(0), BondOrder::Single),
            ],
        );
        let mut charged_view = CipMol::new(&charged_ring);
        let fractions = (0..5)
            .map(|atom| {
                charged_view
                    .get_fractional_atomic_num(atom)
                    .unwrap()
                    .tuple()
            })
            .collect::<Vec<_>>();
        assert_eq!(fractions, vec![(0, 1), (6, 1), (6, 1), (4, 1), (9, 2)]);

        let aromatic_triangle = aromatic_ring(&[Element::C; 3]);
        let mut fallback = CipMol::new(&aromatic_triangle);
        assert_eq!(
            (0..3)
                .map(|bond| fallback.get_bond_order(bond).unwrap())
                .collect::<Vec<_>>(),
            vec![1, 1, 1]
        );
    }
}

#[cfg(test)]
mod recovery_chem02_seen {
    use super::*;
    use cosmolkit_model::{AtomId, AtomSpec, BondSpec};
    fn line() -> TopologyBlock {
        let atoms = (0..3)
            .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
            .collect();
        let bonds = (0..2)
            .map(|i| {
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(i), AtomId::new(i + 1), BondOrder::Single),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
    }
    #[test]
    fn pure_created_nodes_only_including_duplicate_and_terminal() {
        let m = line();
        let mut g = CipDigraph::new(&m, 0, false).unwrap();
        let before = (
            g.nodes.clone(),
            g.edges.clone(),
            g.root,
            g.origin,
            g.rule6_ref,
        );
        assert!(g.seen_atom(0));
        assert!(!g.seen_atom(1));
        assert!(!g.seen_atom(2));
        assert!(!g.seen_atom(usize::MAX));
        assert_eq!(
            (
                g.nodes.clone(),
                g.edges.clone(),
                g.root,
                g.origin,
                g.rule6_ref
            ),
            before
        );
        g.node_edges(g.get_current_root()).unwrap();
        assert!(g.seen_atom(1));
        assert!(!g.seen_atom(2));
        g.nodes
            .iter_mut()
            .find(|n| n.atom_idx() == Some(1))
            .unwrap()
            .flags |= CipNode::RING_DUPLICATE;
        let before = (
            g.nodes.clone(),
            g.edges.clone(),
            g.root,
            g.origin,
            g.rule6_ref,
        );
        assert!(g.seen_atom(1));
        assert!(!g.seen_atom(2));
        assert_eq!(
            (
                g.nodes.clone(),
                g.edges.clone(),
                g.root,
                g.origin,
                g.rule6_ref
            ),
            before
        );
    }
}

#[cfg(test)]
mod recovery_chem04 {
    use super::*;
    #[test]
    fn rdkit_2026_03_6_chem04_pairlist_source_boundaries() {
        // Official .6 Pairlist.h UBSan probe: all 18 descriptor kinds and unbounded vectors.
        for d in Descriptor::ALL_IN_RDKIT_ORDER {
            let expected = match d {
                Descriptor::R | Descriptor::M | Descriptor::seqCis => Some(Descriptor::R),
                Descriptor::S | Descriptor::P | Descriptor::seqTrans => Some(Descriptor::S),
                _ => None,
            };
            let mut list = CipPairList::new();
            assert_eq!(list.add(d), expected.is_some());
            match expected {
                Some(d) => {
                    assert_eq!(list.descriptors, vec![d]);
                    assert_eq!(list.get_ref_descriptor(), Ok(d));
                    assert_eq!(
                        list.to_rdkit_string(),
                        format!("{}:", descriptor_to_string(d))
                    );
                }
                None => {
                    let error = list.get_ref_descriptor().unwrap_err();
                    assert_eq!(error, CipLabelerError::EmptyPairListReference);
                    assert_eq!(
                        error.to_string(),
                        "Cannot get a reference from an empty PairList"
                    );
                    assert!(list.descriptors.is_empty());
                    assert_eq!(list.to_rdkit_string(), "");
                }
            }
        }
        for n in [1, 64, 65, 130] {
            let mut same = CipPairList::new();
            let mut alternating = CipPairList::new();
            for i in 0..n {
                same.add(Descriptor::R);
                alternating.add(if i % 2 == 0 {
                    Descriptor::R
                } else {
                    Descriptor::S
                });
            }
            let suffix: String = (1..n).map(|i| if i % 2 == 0 { 'l' } else { 'u' }).collect();
            assert_eq!(same.to_rdkit_string(), format!("R:{}", "l".repeat(n - 1)));
            assert_eq!(alternating.to_rdkit_string(), format!("R:{suffix}"));
            let cmp = if n == 1 { 0 } else { 1 };
            assert_eq!(same.compare_to(&alternating), Ok(cmp));
            assert_eq!(alternating.compare_to(&same), Ok(-cmp));
            assert_eq!(same.compare_to(&same), Ok(0));
            let combined = CipPairList::from_head_tail(&same, &alternating);
            let mut expected = same.descriptors.clone();
            expected.extend_from_slice(&alternating.descriptors);
            assert_eq!(combined.descriptors, expected);
            assert_eq!(combined.descriptors.len(), n * 2);
            assert_eq!(
                combined.to_rdkit_string(),
                format!("R:{}l{suffix}", "l".repeat(n - 1))
            );
            assert_eq!(
                CipPairList::from_head_tail(&CipPairList::new(), &same),
                same
            );
            assert_eq!(
                CipPairList::from_head_tail(&same, &CipPairList::new()),
                same
            );
        }
        let mut empty = vec![CipPairList::new(); 3];
        CipPairList::sort_descending(&mut empty).unwrap();
        assert_eq!(empty, vec![CipPairList::new(); 3]);
        CipPairList::sort_descending(&mut []).unwrap();
    }

    #[test]
    fn fill_pairs_empty_reference_precedes_sorter_and_does_not_mutate_graph() {
        let atoms = (0..2)
            .map(|i| {
                Atom::from_spec(
                    cosmolkit_model::AtomId::new(i),
                    cosmolkit_model::AtomSpec::new(Element::C),
                )
            })
            .collect();
        let b = Bond::from_spec(
            BondId::new(0),
            cosmolkit_model::BondSpec::new(
                cosmolkit_model::AtomId::new(0),
                cosmolkit_model::AtomId::new(1),
                BondOrder::Single,
            ),
        );
        let t = TopologyBlock::try_from_parts(atoms, vec![b], vec![], vec![]).unwrap();
        for missing in [false, true] {
            for five in [false, true] {
                let mut g = CipDigraph::new(&t, 0, false).unwrap();
                let root = g.get_current_root();
                let mut c = CipLabelerContext::new(77);
                let before = (
                    g.nodes.clone(),
                    g.edges.clone(),
                    g.root,
                    g.origin,
                    g.rule6_ref,
                    c.remaining_call_count,
                );
                let mut list = CipPairList::new();
                let unrelated = CipRule1a;
                let refs: [&dyn CipSequenceRule; 1] = [&unrelated];
                let sort = if missing { None } else { Some(refs.as_slice()) };
                let r = if five {
                    CipRule5New::new().fill_pairs(sort, &mut g, &mut c, root, &mut list)
                } else {
                    CipRule4b::new().fill_pairs(sort, &mut g, &mut c, root, &mut list)
                };
                assert_eq!(r, Err(CipLabelerError::EmptyPairListReference));
                assert_eq!(list, CipPairList::new());
                assert_eq!(
                    (
                        g.nodes.clone(),
                        g.edges.clone(),
                        g.root,
                        g.origin,
                        g.rule6_ref,
                        c.remaining_call_count
                    ),
                    before
                );
            }
        }
        assert_eq!(CipPairList::new().compare_to(&CipPairList::new()), Ok(0));
        assert_eq!(
            CipPairList::new().compare_to(&CipPairList::with_ref(Descriptor::R)),
            Err(CipLabelerError::DescriptorListLengthMismatch)
        );
    }
}

#[cfg(test)]
mod recovery_chem06 {
    use super::*;
    #[derive(Debug, PartialEq)]
    struct Chem02PassiveState {
        nodes: Vec<CipNode>,
        edges: Vec<(CipNodeId, CipNodeId, Option<usize>, Descriptor)>,
        root: CipNodeId,
        origin: CipNodeId,
        rule6_ref: Option<usize>,
        atropisomer_mode: bool,
        remaining_calls: u32,
    }

    fn chem02_passive_state(g: &CipDigraph<'_>, c: &CipLabelerContext) -> Chem02PassiveState {
        // Observation must never trigger node_edges/get_nodes/change_root.
        Chem02PassiveState {
            nodes: g.nodes.clone(),
            edges: g
                .edges
                .iter()
                .map(|e| (e.beg, e.end, e.bond_idx, e.aux))
                .collect(),
            root: g.root,
            origin: g.origin,
            rule6_ref: g.rule6_ref,
            atropisomer_mode: g.atropisomer_mode,
            remaining_calls: c.remaining_call_count,
        }
    }

    fn chem06_graph(s: &str) -> TopologyBlock {
        use cosmolkit_model::{AtomId, AtomSpec, BondSpec};
        let (elements, edges) = match s {
            "C(C)(F)Cl" => (
                vec![Element::C, Element::C, Element::F, Element::CL],
                vec![
                    (0, 1, BondOrder::Single),
                    (0, 2, BondOrder::Single),
                    (0, 3, BondOrder::Single),
                ],
            ),
            "C=C" => (
                vec![Element::C, Element::C],
                vec![(0, 1, BondOrder::Double)],
            ),
            "CCC" => (
                vec![Element::C; 3],
                vec![(0, 1, BondOrder::Single), (1, 2, BondOrder::Single)],
            ),
            _ => unreachable!("three fixed typed source graphs"),
        };
        let atoms = elements
            .into_iter()
            .enumerate()
            .map(|(i, e)| Atom::from_spec(AtomId::new(i), AtomSpec::new(e)))
            .collect();
        let bonds = edges
            .into_iter()
            .enumerate()
            .map(|(i, (a, b, o))| {
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(a), AtomId::new(b), o),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
    }
    fn chem06_edge_to(g: &CipDigraph<'_>, edges: &[CipEdgeId], atom: Option<usize>) -> CipEdgeId {
        edges
            .iter()
            .copied()
            .find(|e| {
                let n = g.node(g.edge(*e).get_end());
                n.atom_idx() == atom && !n.is_duplicate()
            })
            .unwrap()
    }

    fn chem06_sort_with_owner(
        owner: usize,
        g: &mut CipDigraph<'_>,
        c: &mut CipLabelerContext,
        root: CipNodeId,
        edges: &mut [CipEdgeId],
        deep: bool,
    ) -> Result<CipPriority, CipLabelerError> {
        match owner {
            0 => CipRule6.sort(g, c, root, edges, deep),
            1 => CipSort::new(&CipRule6).prioritize(g, c, root, edges, deep),
            2 => CipRules::new(vec![Box::new(CipRule6)])
                .unwrap()
                .sort(g, c, root, edges, deep),
            _ => unreachable!("three existing production sorter owners"),
        }
    }

    #[test]
    fn rdkit_2026_03_6_chem06_official_9516_rule_and_sorters() {
        // Official #9516 input and controlled reference state; no full9rule molecular-label claim.
        let molecule = chem06_graph("C(C)(F)Cl");
        let mut g = CipDigraph::new(&molecule, 0, false).unwrap();
        let root = g.get_current_root();
        let edges = g.node_edges(root).unwrap();
        let a = chem06_edge_to(&g, &edges, Some(1));
        let b = chem06_edge_to(&g, &edges, Some(2));
        assert_eq!(
            edges
                .iter()
                .filter(|e| g.node(g.edge(**e).get_end()).atom_idx().is_none())
                .count(),
            1
        );
        g.set_rule6_ref(Some(1)).unwrap();
        let mut c = CipLabelerContext::with_remaining_call_count(100);
        assert_eq!(CipRule6.compare(&mut g, &mut c, a, b), Ok(2));
        assert_eq!(CipRule6.compare(&mut g, &mut c, b, a), Ok(-2));
        for owner in 0..3 {
            let mut order = [b, a];
            let priority =
                chem06_sort_with_owner(owner, &mut g, &mut c, root, &mut order, false).unwrap();
            assert!(priority.is_unique());
            assert!(priority.is_pseudo_asymetric());
            assert_eq!(order, [a, b]);
        }
        // All actual raw forwarding APIs retain ±2; recursive exclusive-ref stops before expansion.
        for (left, right, expected) in [(a, b, 2), (b, a, -2)] {
            let before_nodes = g.nodes.clone();
            let before_edges = g.edges.clone();
            let rules = CipRules::new(vec![Box::new(CipRule6)]).unwrap();
            let refs: [&dyn CipSequenceRule; 1] = [&CipRule6];
            assert_eq!(
                CipRule6.get_comparison(&mut g, &mut c, left, right, true),
                Ok(expected)
            );
            assert_eq!(
                CipRule6.get_comparison_with_sort_rules(
                    Some(&refs),
                    &mut g,
                    &mut c,
                    left,
                    right,
                    true
                ),
                Ok(expected)
            );
            assert_eq!(rules.compare(&mut g, &mut c, left, right), Ok(expected));
            assert_eq!(
                rules.get_comparison(&mut g, &mut c, left, right, false),
                Ok(expected)
            );
            assert_eq!(
                rules.get_comparison_with_sort_rules(None, &mut g, &mut c, left, right, false),
                Ok(expected)
            );
            assert_eq!(g.nodes, before_nodes);
            assert_eq!(g.edges, before_edges);
        }
        let groups = CipSort::new(&CipRule6)
            .get_groups(&mut g, &mut c, &[a, b])
            .unwrap();
        assert_eq!(groups, vec![vec![a], vec![b]]);
        // Direct no-ref0 does not imply deep0: CH3 versus F has different child-edge counts.
        g.set_rule6_ref(None).unwrap();
        assert_eq!(CipRule6.compare(&mut g, &mut c, a, b), Ok(0));
        assert_eq!(CipRule6.get_comparison(&mut g, &mut c, a, b, true), Ok(1));
        assert_eq!(CipRule6.get_comparison(&mut g, &mut c, b, a, true), Ok(-1));
    }

    #[test]
    fn rdkit_2026_03_6_chem06_rule6_reference_identity_matrix() {
        let molecule = chem06_graph("C(C)(F)Cl");
        let mut g = CipDigraph::new(&molecule, 0, false).unwrap();
        let root = g.get_current_root();
        let edges = g.node_edges(root).unwrap();
        let a = chem06_edge_to(&g, &edges, Some(1));
        let b = chem06_edge_to(&g, &edges, Some(2));
        let h = chem06_edge_to(&g, &edges, None);
        let mut c = CipLabelerContext::new(0);
        let rows = [
            (None, a, b, 0),
            (Some(1), a, b, 2),
            (Some(2), a, b, -2),
            (Some(0), a, b, 0),
            (Some(1), a, a, 0),
            (Some(1), a, h, 2),
            (Some(1), h, a, -2),
        ];
        for (reference, left, right, expected) in rows {
            g.set_rule6_ref(reference).unwrap();
            let state = chem02_passive_state(&g, &c);
            assert_eq!(CipRule6.compare(&mut g, &mut c, left, right), Ok(expected));
            assert_eq!(chem02_passive_state(&g, &c), state);
        }
        // Real source graph construction yields distinct nodes sharing an atom plus two None Hs.
        let ethene = chem06_graph("C=C");
        let mut g = CipDigraph::new(&ethene, 0, true).unwrap();
        let root = g.get_current_root();
        let edges = g.node_edges(root).unwrap();
        let normal = chem06_edge_to(&g, &edges, Some(1));
        let duplicate = edges
            .iter()
            .copied()
            .find(|e| {
                let n = g.node(g.edge(*e).get_end());
                n.atom_idx() == Some(1) && n.is_duplicate()
            })
            .unwrap();
        assert_ne!(g.edge(normal).get_end(), g.edge(duplicate).get_end());
        let hs: Vec<_> = edges
            .iter()
            .copied()
            .filter(|e| g.node(g.edge(*e).get_end()).atom_idx().is_none())
            .collect();
        assert_eq!(hs.len(), 2);
        assert_ne!(g.edge(hs[0]).get_end(), g.edge(hs[1]).get_end());
        g.set_rule6_ref(Some(1)).unwrap();
        let state = chem02_passive_state(&g, &c);
        assert_eq!(CipRule6.compare(&mut g, &mut c, normal, duplicate), Ok(0));
        assert_eq!(CipRule6.compare(&mut g, &mut c, hs[0], hs[1]), Ok(0));
        assert_eq!(chem02_passive_state(&g, &c), state);
    }

    #[test]
    fn rdkit_2026_03_6_chem06_pseudo_counts_comparisons_not_any() {
        let molecule = chem06_graph("C(C)(F)Cl");
        // Source insertion comparisons: [2,0] counts1; [0,-2,-2] counts2.
        let cases = [
            (vec![], Some(1), true, false, vec![], 0),
            (vec![1], Some(1), true, false, vec![1], 0),
            (vec![1, 2], Some(1), true, true, vec![1, 2], 1),
            (vec![2, 1], Some(1), true, true, vec![1, 2], 1),
            (vec![1, 2], None, false, false, vec![1, 2], 1),
            (vec![1, 2, 3], Some(1), false, true, vec![1, 2, 3], 2),
            (vec![2, 3, 1], Some(1), false, false, vec![1, 2, 3], 3),
        ];
        for owner in 0..3 {
            for (case, (initial, reference, unique, pseudo, expected, comparisons)) in
                cases.iter().enumerate()
            {
                let mut g = CipDigraph::new(&molecule, 0, false).unwrap();
                let root = g.get_current_root();
                let edges = g.node_edges(root).unwrap();
                g.set_rule6_ref(*reference).unwrap();
                let mut c = CipLabelerContext::with_remaining_call_count(100);
                let mut order: Vec<_> = initial
                    .iter()
                    .map(|atom| chem06_edge_to(&g, &edges, Some(*atom)))
                    .collect();
                // Rules::getComparision explicitly ignores deep; its own sorter is Sort(this).
                // Build exact expected state using only the source-directed child expansions.
                let mut expected_g = CipDigraph::new(&molecule, 0, false).unwrap();
                let expected_root = expected_g.get_current_root();
                let expected_edges = expected_g.node_edges(expected_root).unwrap();
                expected_g.set_rule6_ref(*reference).unwrap();
                let mut expected_c = CipLabelerContext::with_remaining_call_count(100);
                if owner == 2 {
                    expected_c.remaining_call_count -= *comparisons;
                    let expanded_atoms: &[usize] = match case {
                        4 => &[1, 2],     // no-ref recursive CH3/F comparison reaches source sizediff+1
                        5 | 6 => &[2, 3], // F/Cl direct0 recursively expands equal terminal branches
                        _ => &[],         // exclusive reference returns before graph expansion
                    };
                    for atom in expanded_atoms {
                        let edge = chem06_edge_to(&expected_g, &expected_edges, Some(*atom));
                        let end = expected_g.edge(edge).get_end();
                        expected_g.node_edges(end).unwrap();
                    }
                }
                let p =
                    chem06_sort_with_owner(owner, &mut g, &mut c, root, &mut order, false).unwrap();
                let expected_unique = if owner == 2 && case == 4 {
                    true
                } else {
                    *unique
                };
                assert_eq!(p.is_unique(), expected_unique, "owner={owner} case={case}");
                assert_eq!(
                    p.is_pseudo_asymetric(),
                    *pseudo,
                    "owner={owner} case={case}"
                );
                assert_eq!(
                    order
                        .iter()
                        .map(|e| g.node(g.edge(*e).get_end()).atom_idx().unwrap())
                        .collect::<Vec<_>>(),
                    *expected
                );
                assert_eq!(
                    chem02_passive_state(&g, &c),
                    chem02_passive_state(&expected_g, &expected_c),
                    "owner={owner} case={case}"
                );
            }
        }
    }

    #[test]
    fn rdkit_2026_03_6_chem06_incoming_edge_priority_preempts_rule6() {
        let molecule = chem06_graph("CCC");
        let mut g = CipDigraph::new(&molecule, 0, false).unwrap();
        let origin = g.get_current_root();
        let initial = g.node_edges(origin).unwrap();
        let incoming = chem06_edge_to(&g, &initial, Some(1));
        let node = g.edge(incoming).get_end();
        let edges = g.node_edges(node).unwrap();
        let outgoing = chem06_edge_to(&g, &edges, Some(2));
        assert!(!g.edge(incoming).is_beg(node));
        assert!(g.edge(outgoing).is_beg(node));
        g.set_rule6_ref(Some(1)).unwrap();
        let mut c = CipLabelerContext::new(0);
        for (a, b, direct, precheck) in [(incoming, outgoing, 2, 1), (outgoing, incoming, -2, -1)] {
            assert_eq!(CipRule6.compare(&mut g, &mut c, a, b), Ok(direct));
            assert_eq!(
                super::compare_substituents_with_rule(&CipRule6, &mut g, &mut c, node, a, b, false),
                Ok(precheck)
            );
            assert_eq!(
                CipSort::new(&CipRule6).compare_substituents(&mut g, &mut c, node, a, b, false),
                Ok(precheck)
            );
            for owner in 0..3 {
                let mut order = [a, b];
                let p =
                    chem06_sort_with_owner(owner, &mut g, &mut c, node, &mut order, false).unwrap();
                assert_eq!(order, [incoming, outgoing]);
                assert!(p.is_unique());
                assert!(!p.is_pseudo_asymetric());
            }
        }
    }
}
