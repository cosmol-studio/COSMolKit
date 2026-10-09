//! Private representation transport for detached source conformer pruning.
use cosmolkit_model::{
    AtomQueryPredicate, BondQueryPredicate, CoordinateBlock, MoleculeProperties, QueryAtom,
    QueryBond, QueryGraph, QueryGraphError, QueryNode, TopologyBlock,
};
/// Own only carrier representation transport, with no search/chemistry/runtime calls.
pub(crate) fn concrete_pruning_query(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    properties: &MoleculeProperties,
) -> Result<QueryGraph, QueryGraphError> {
    let atoms = topology
        .atoms
        .iter()
        .map(|atom| {
            QueryAtom::from_carrier_parts(
                atom.clone(),
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(atom.atomic_number())),
            )
        })
        .collect();
    let bonds = topology
        .bonds
        .iter()
        .map(|bond| {
            QueryBond::from_carrier_parts(
                bond.clone(),
                QueryNode::predicate(BondQueryPredicate::Order(bond.order())),
            )
        })
        .collect();
    let mut props = properties.props().clone();
    if let Some(name) = properties.name() {
        props.insert(
            "_Name".into(),
            cosmolkit_model::PropertyValue::String(name.clone()),
        );
    }
    let mut query = QueryGraph::from_parts(
        atoms,
        bonds,
        props,
        coordinates.conformers_2d.clone(),
        coordinates.conformers_3d.clone(),
        topology.stereo_groups.clone(),
    )?;
    query.set_source_conformer_order(coordinates.source_conformer_order.clone())?;
    cosmolkit_model::replace_query_substance_groups(&mut query, topology.substance_groups.clone())?;
    Ok(query)
}
#[cfg(all(test, feature = "cap-alignment", feature = "cap-search"))]
mod tests {
    use super::*;
    use cosmolkit_model::{
        AdjacencyList, Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, Conformer3D,
        Element,
    };
    use cosmolkit_search::{SearchTarget, SubstructMatchParams};
    #[test]
    fn pruning_query_preserves_full_source_order_and_first_conformer() {
        let topology = fixed_topology(vec![AtomSpec::new(Element::C)], vec![]);
        let coordinates = CoordinateBlock {
            conformers_2d: vec![cosmolkit_model::Conformer2D::new(9, vec![[1.0, -0.0]])],
            conformers_3d: vec![Conformer3D::new(9, vec![[2.0, 3.0, 4.0]], false)],
            source_conformer_order: Some(vec![
                cosmolkit_model::CoordinateDimension::ThreeD,
                cosmolkit_model::CoordinateDimension::TwoD,
            ]),
            ..Default::default()
        };
        let query = concrete_pruning_query(&topology, &coordinates, &MoleculeProperties::default())
            .unwrap();
        assert_eq!(query.coordinate_block(None), coordinates);
        match query
            .coordinate_block(None)
            .first_source_conformer()
            .unwrap()
            .unwrap()
        {
            cosmolkit_model::CoordinateSourceConformer::ThreeD(first) => {
                assert_eq!(first.id(), 9);
                assert!(!first.is_3d());
            }
            _ => panic!("source first identity must survive representation transport"),
        }
    }

    fn fixed_topology(
        atoms: Vec<AtomSpec>,
        bonds: Vec<(usize, usize, BondOrder)>,
    ) -> TopologyBlock {
        let atoms = atoms
            .into_iter()
            .enumerate()
            .map(|(i, spec)| Atom::from_spec(AtomId::new(i), spec))
            .collect::<Vec<_>>();
        let bonds = bonds
            .into_iter()
            .enumerate()
            .map(|(i, (a, b, order))| {
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(a), AtomId::new(b), order),
                )
            })
            .collect::<Vec<_>>();
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
            atoms,
            bonds,
            ..Default::default()
        };
        topology
            .validate()
            .expect("original checked fixed topology");
        topology
    }
    fn terminal_group_symmetry_topology() -> TopologyBlock {
        fixed_topology(
            vec![
                AtomSpec::new(Element::O).with_formal_charge(-1),
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::O),
            ],
            vec![(0, 1, BondOrder::Single), (1, 2, BondOrder::Double)],
        )
    }
    #[test]
    fn embedder_mol_self_matches_uses_heavy_atoms_or_all_atoms_like_rdkit() {
        let mol = fixed_topology(
            vec![AtomSpec::new(Element::C), AtomSpec::new(Element::H)],
            vec![(0, 1, BondOrder::Single)],
        );
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        let mut heavy = cosmolkit_conformer::EmbedParams::default();
        heavy.only_heavy_atoms_for_rms = true;
        heavy.use_symmetry_for_pruning = false;
        assert_eq!(
            cosmolkit_conformer::pruning_self_matches(
                &mol,
                &coordinates,
                &properties,
                &heavy,
                concrete_pruning_query
            )
            .unwrap()
            .self_matches,
            vec![vec![0]]
        );
        let mut all = cosmolkit_conformer::EmbedParams::default();
        all.only_heavy_atoms_for_rms = false;
        all.use_symmetry_for_pruning = false;
        assert_eq!(
            cosmolkit_conformer::pruning_self_matches(
                &mol,
                &coordinates,
                &properties,
                &all,
                concrete_pruning_query
            )
            .unwrap()
            .self_matches,
            vec![vec![0, 1]]
        );
    }
    #[test]
    fn symmetrize_terminal_atoms_for_pruning_matches_rdkit_terminal_group_query() {
        let mol = terminal_group_symmetry_topology();
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        let target = SearchTarget::new(&mol, &coordinates, &mol.stereo_groups, None, None);
        let terminal_context = cosmolkit_search::build_topology_query_match_context(&mol).unwrap();
        let symmetrized = cosmolkit_alignment::symmetrize_terminal_query_with_context(
            concrete_pruning_query(&mol, &coordinates, &properties).unwrap(),
            &target,
            &terminal_context,
        )
        .expect("symmetrized");
        assert_eq!(symmetrized.atoms()[0].formal_charge(), 0);
        assert_eq!(symmetrized.atoms()[2].formal_charge(), 0);
        assert_eq!(
            symmetrized.bonds()[0].predicate(),
            &QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                BondOrder::Single,
                BondOrder::Double
            ]))
        );
        assert_eq!(
            symmetrized.bonds()[1].predicate(),
            &QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                BondOrder::Single,
                BondOrder::Double
            ]))
        );
        assert!(!symmetrized.bonds()[0].predicate_is_carrier_derived());
        let params = SubstructMatchParams {
            max_matches: 1000,
            uniquify: false,
            use_chirality: false,
            specified_stereo_query_matches_unspecified: false,
            ..Default::default()
        };
        let context = cosmolkit_search::build_topology_query_match_context(&mol).unwrap();
        let matches = cosmolkit_search::try_get_substruct_matches_with_params_and_context(
            &target,
            &symmetrized,
            &params,
            &context,
        )
        .unwrap();
        let atom_mappings: Vec<_> = matches
            .into_iter()
            .map(|matched| matched.atom_mapping)
            .collect();
        assert_eq!(atom_mappings.len(), 2);
        assert!(atom_mappings.contains(&vec![0, 1, 2]));
        assert!(atom_mappings.contains(&vec![2, 1, 0]));
    }
    #[test]
    fn embedder_mol_self_matches_symmetrize_terminal_groups_for_pruning() {
        let mol = terminal_group_symmetry_topology();
        let mut params = cosmolkit_conformer::EmbedParams::default();
        params.prune_rms_thresh = 0.5;
        params.use_symmetry_for_pruning = true;
        params.symmetrize_conjugated_terminal_groups_for_pruning = true;
        let matches = cosmolkit_conformer::pruning_self_matches(
            &mol,
            &CoordinateBlock::default(),
            &MoleculeProperties::default(),
            &params,
            concrete_pruning_query,
        )
        .expect("self matches")
        .self_matches;
        assert_eq!(matches.len(), 2);
        assert!(matches.contains(&vec![0, 1, 2]));
        assert!(matches.contains(&vec![2, 1, 0]));
    }
    #[test]
    fn private_carrier_projection_preserves_origins_and_all_chemistry_rows() {
        let mol = fixed_topology(
            vec![
                AtomSpec::new(Element::C)
                    .with_isotope(13)
                    .with_formal_charge(1),
                AtomSpec::new(Element::O),
            ],
            vec![(0, 1, BondOrder::Double)],
        );
        let coordinates = CoordinateBlock {
            conformers_3d: vec![Conformer3D::new(
                7,
                vec![[0.0, 1.0, 2.0], [3.0, 4.0, 5.0]],
                true,
            )],
            ..Default::default()
        };
        let properties = MoleculeProperties::default().with_name("carrier-original");
        let query = concrete_pruning_query(&mol, &coordinates, &properties).unwrap();
        assert!(
            query
                .atoms()
                .iter()
                .all(|atom| atom.predicate_is_carrier_derived())
        );
        assert!(
            query
                .bonds()
                .iter()
                .all(|bond| bond.predicate_is_carrier_derived())
        );
        assert_eq!(query.atoms()[0].isotope(), Some(13));
        assert_eq!(query.atoms()[0].formal_charge(), 1);
        assert_eq!(query.bonds()[0].bond(), &mol.bonds[0]);
        assert_eq!(query.conformers_3d(), coordinates.conformers_3d.as_slice());
        assert_eq!(
            query.name().unwrap().map(|name| name.as_bytes()),
            Some(b"carrier-original".as_slice())
        );
    }
    #[test]
    fn source_self_matches_keep_degree_zero_hydrogen_and_nonzero_truth_condition() {
        let mol = fixed_topology(
            vec![AtomSpec::new(Element::C), AtomSpec::new(Element::H)],
            vec![],
        );
        let mut params = cosmolkit_conformer::EmbedParams::default();
        params.prune_rms_thresh = -1.0;
        params.only_heavy_atoms_for_rms = true;
        params.use_symmetry_for_pruning = true;
        params.symmetrize_conjugated_terminal_groups_for_pruning = false;
        let matches = cosmolkit_conformer::pruning_self_matches(
            &mol,
            &CoordinateBlock::default(),
            &MoleculeProperties::default(),
            &params,
            concrete_pruning_query,
        )
        .unwrap();
        assert_eq!(matches.self_matches, vec![vec![0, 1]]);
        assert!(!matches.hydrogen_warnings.is_empty());
    }
    #[test]
    fn whole_generation_uses_canonical_projection_after_remove_hs_for_retained_pruning() {
        let record = cosmolkit_smiles::parse_smiles("CCO", &Default::default()).unwrap();
        let topology = cosmolkit_core::sanitize_topology(&record.topology, &Default::default())
            .unwrap()
            .topology;
        let mut molecule =
            cosmolkit_core::add_hydrogens_impl(topology, Default::default(), Default::default())
                .unwrap();
        let mut params = cosmolkit_conformer::EmbedParams::etkdg_v3();
        params.random_seed = 42;
        params.max_iterations = 3;
        params.num_threads = 1;
        params.prune_rms_thresh = 1000.0;
        let first = cosmolkit_conformer::generate_conformers(
            &molecule.topology,
            &molecule.coordinates,
            &molecule.properties,
            2,
            &mut params,
            concrete_pruning_query,
        )
        .unwrap();
        assert_eq!(first.conf_ids, vec![0]);
        assert_eq!(first.conformers.len(), 1);
        molecule.coordinates.conformers_3d = first.conformers;
        let before = molecule.coordinates.clone();
        params.clear_conformers = false;
        let second = cosmolkit_conformer::generate_conformers(
            &molecule.topology,
            &molecule.coordinates,
            &molecule.properties,
            2,
            &mut params,
            concrete_pruning_query,
        )
        .unwrap();
        assert!(second.conf_ids.is_empty());
        assert!(second.conformers.is_empty());
        assert!(!second.clear_existing);
        assert_eq!(molecule.coordinates, before);
    }
}
