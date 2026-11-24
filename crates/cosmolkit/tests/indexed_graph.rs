#![cfg(all(
    feature = "cap-smiles",
    feature = "cap-transforms",
    feature = "cap-matrices",
    feature = "cap-kekulize"
))]
use cosmolkit::{
    AtomId, AtomRankParams, AtomSpec, BondCut, BondId, BondOrder, BondSpec, Conformer3D,
    CoordinateBlock, Element, GraphPath, Molecule, PathEnumerationParams, SubgraphExtractionParams,
    SubgraphSelection,
};

fn mol(s: &str) -> Molecule {
    Molecule::from_smiles(s).unwrap()
}
fn a(i: usize) -> AtomId {
    AtomId::new(i)
}
fn b(i: usize) -> BondId {
    BondId::new(i)
}
fn with_x(s: &str, sets: &[(usize, &[f64])]) -> Molecule {
    let m = mol(s);
    let coords = CoordinateBlock {
        conformers_3d: sets
            .iter()
            .map(|(id, x)| Conformer3D::new(*id, x.iter().map(|&v| [v, 0., 0.]).collect(), true))
            .collect(),
        ..Default::default()
    };
    Molecule::from_parts(m.topology().clone(), coords, m.properties().clone()).unwrap()
}
fn xyz(m: &Molecule, id: usize) -> Vec<f64> {
    m.coordinates_3d(id).unwrap().iter().map(|r| r[0]).collect()
}

#[test]
fn ranks_are_source_indexed_and_fragment_masks_are_explicit() {
    let m = mol("CC(C)O");
    assert_eq!(m.canonical_atom_ranks().unwrap(), [0, 3, 1, 2]);
    let p = AtomRankParams {
        break_ties: false,
        ..Default::default()
    };
    assert_eq!(
        m.canonical_atom_ranks_with_params(&p).unwrap(),
        [0, 3, 0, 2]
    );
    let atoms = [a(1), a(2), a(3)];
    for (bonds, ranks) in [
        (vec![], [-1, 1, 1, 3]),
        (vec![b(1)], [-1, 2, 2, 1]),
        (vec![b(1), b(2)], [-1, 3, 1, 2]),
        (vec![b(0), b(1), b(2)], [-1, 3, 1, 2]),
    ] {
        assert_eq!(
            m.canonical_atom_ranks_in_fragment_with_params(&atoms, &bonds, None, &p)
                .unwrap(),
            ranks
        );
    }
    assert!(m.canonical_atom_ranks_in_fragment(&[], &[]).is_err());
    assert!(
        m.canonical_atom_ranks_in_fragment_with_params(&atoms, &[], Some(&[]), &p)
            .is_err()
    );
    let m = mol("[CH3:1]C");
    let meso = mol("C[C@H](O)[C@@H](C)O");
    assert_eq!(
        meso.canonical_atom_ranks_in_fragment_with_params(&[a(1), a(3)], &[], None, &p)
            .unwrap(),
        [-1, 3, -1, 3, -1, -1]
    );
    for (maps, ranks) in [(true, [1, 0]), (false, [0, 0])] {
        let p = AtomRankParams {
            include_atom_maps: maps,
            ..p.clone()
        };
        assert_eq!(
            m.canonical_atom_ranks_in_fragment_with_params(&[a(0), a(1)], &[b(0)], None, &p)
                .unwrap(),
            ranks
        );
    }
    let cut = mol("C1CC1CO")
        .with_bonds_cut(&[BondCut::at_bond(b(3), None, None).unwrap()])
        .unwrap();
    assert_eq!(cut.canonical_atom_ranks().unwrap(), [4, 5, 6, 3, 2, 0, 1]);
    let m = mol("CC/C=C/C1=C2CC(C)CN=C2CC1");
    let params = AtomRankParams {
        include_chirality: false,
        ..p
    };
    let fragment = m
        .canonical_atom_ranks_in_fragment_with_params(
            &[0, 1, 2, 3, 4, 5, 6, 9, 12].map(a),
            &[0, 2, 3, 4, 6, 8, 11, 12, 14].map(b),
            None,
            &params,
        )
        .unwrap();
    assert_eq!(
        fragment,
        [10, 10, 8, 12, 12, 8, 5, -1, -1, 5, -1, -1, 5, -1]
    );
}

#[test]
fn traversal_keeps_kind_order_and_disconnected_atoms() {
    let m = mol("CC.O");
    assert_eq!(m.bond_between_atoms(a(1), a(0)).unwrap(), Some(b(0)));
    assert_eq!(m.bond_between_atoms(a(0), a(2)).unwrap(), None);
    assert!(m.bond_between_atoms(a(0), a(3)).is_err());
    assert_eq!(m.bond_between_atoms(a(0), a(0)).unwrap(), None);
    assert!(m.shortest_path(a(0), a(2)).unwrap().is_empty());
    assert!(m.shortest_path(a(0), a(0)).is_err());
    assert_eq!(
        mol("C1CCC1").shortest_path(a(0), a(2)).unwrap(),
        [a(0), a(1), a(2)]
    );
    assert_eq!(
        mol("[NH3]->[Cu+2]")
            .adjacency_matrix(true)
            .unwrap()
            .values(),
        [0., 0., 1., 0.]
    );
    let m = mol("CCC");
    assert_eq!(
        m.paths_of_length(2).unwrap(),
        [GraphPath::Bonds(vec![b(0), b(1)])]
    );
    assert_eq!(
        m.paths_of_length_with_params(
            2,
            &PathEnumerationParams {
                use_atom_indices: true,
                ..Default::default()
            }
        )
        .unwrap(),
        [
            GraphPath::Atoms(vec![a(0), a(1)]),
            GraphPath::Atoms(vec![a(1), a(2)])
        ]
    );
    assert_eq!(
        m.subgraphs_in_range(0, 3).unwrap(),
        [
            vec![],
            vec![vec![b(0)], vec![b(1)]],
            vec![vec![b(0), b(1)]],
            vec![]
        ]
    );
    assert!(m.subgraphs_in_range(2, 1).is_err());
    let m = mol("CC(C)C");
    assert!(m.paths_of_length(3).unwrap().is_empty());
    assert_eq!(m.subgraphs_of_length(3).unwrap(), [vec![b(0), b(2), b(1)]]);
    assert_eq!(m.unique_subgraphs_of_length(2).unwrap(), [vec![b(0), b(2)]]);
    let m = mol("CCO");
    assert_eq!(
        m.atom_environment(0, a(1)).unwrap().atom_distances(),
        &[(a(1), 0)].into()
    );
    assert_eq!(m.atom_environment(1, a(1)).unwrap().bonds(), [b(0), b(1)]);
    assert!(
        m.atom_environment(3, a(1))
            .unwrap()
            .atom_distances()
            .is_empty()
    );
}

#[test]
fn permutation_then_extraction_preserves_row_identity() {
    let m = with_x("OCC", &[(4, &[1., 2., 3.]), (9, &[11., 12., 13.])]);
    let mut edit = m.clone();
    for bad in [
        vec![a(0), a(0), a(2)],
        vec![a(1), a(0)],
        vec![a(0), a(1), a(3)],
    ] {
        assert!(edit.renumber_atoms_(&bad).is_err());
        assert_eq!(edit, m);
    }
    let ranks = m.canonical_atom_ranks().unwrap();
    let mut order = [a(0), a(1), a(2)];
    order.sort_by_key(|a| ranks[a.index()]);
    assert_eq!(order, [a(2), a(0), a(1)]);
    edit.renumber_atoms_(&order).unwrap();
    assert_eq!(edit, m.with_atoms_renumbered(&order).unwrap());
    assert_eq!(xyz(&edit, 4), [3., 1., 2.]);
    assert_eq!(xyz(&edit, 9), [13., 11., 12.]);
    assert_eq!(
        (
            edit.topology().bonds[1].begin(),
            edit.topology().bonds[1].end()
        ),
        (a(2), a(0))
    );
    let path = edit.shortest_path(a(0), a(2)).unwrap();
    let selection = SubgraphSelection::bonds(
        path.windows(2)
            .map(|p| edit.bond_between_atoms(p[0], p[1]).unwrap().unwrap())
            .collect(),
    );
    let part = edit.extract_subgraph(&selection).unwrap();
    let map = edit.retain_subgraph_(&selection).unwrap();
    assert_eq!(map.source_atom_indices(), [a(0), a(2)]);
    assert_eq!(map.source_bond_indices(), [b(1)]);
    assert_eq!(&map, part.mapping());
    assert_eq!(part.molecule(), &edit);
    assert_eq!(xyz(&edit, 4), [3., 2.]);
    assert_eq!(xyz(&m, 4), [1., 2., 3.]);
}

#[test]
fn subset_membership_and_options_do_not_depend_on_list_lengths() {
    let m = mol("C1CC1");
    let all = vec![a(0), a(1), a(2)];
    for (selection, n) in [
        (SubgraphSelection::Atoms(all.clone()), 3),
        (
            SubgraphSelection::AtomsAndBonds(all.clone(), vec![b(1), b(0), b(1)]),
            2,
        ),
        (SubgraphSelection::AtomsAndBonds(all, vec![]), 0),
    ] {
        assert_eq!(
            m.extract_subgraph(&selection)
                .unwrap()
                .molecule()
                .num_bonds(),
            n
        );
    }
    let bad = SubgraphSelection::AtomsAndBonds(vec![a(0)], vec![b(0)]);
    let mut target = m.clone();
    assert!(target.retain_subgraph_(&bad).is_err());
    assert_eq!(target, m);
    let m = with_x(
        "CCCO",
        &[(4, &[40., 41., 42., 43.]), (9, &[90., 91., 92., 93.])],
    );
    let sel = SubgraphSelection::Atoms(vec![a(3), a(1)]);
    let part = m.extract_subgraph(&sel).unwrap();
    assert_eq!(xyz(part.molecule(), 4), [41., 43.]);
    assert_eq!(xyz(part.molecule(), 9), [91., 93.]);
    let no_coords = m
        .extract_subgraph_with_params(
            &sel,
            &SubgraphExtractionParams {
                copy_coordinates: false,
                ..Default::default()
            },
        )
        .unwrap();
    assert_eq!(no_coords.mapping(), part.mapping());
    assert!(no_coords.molecule().coordinates_3d(4).is_err());
}

#[test]
fn local_edits_keep_surviving_slots_and_union_matches_conformer_ids() {
    let m = with_x("CCO", &[(1, &[1., 2., 3.])]);
    let deleted = m.without_atom(a(1)).unwrap();
    assert_eq!(deleted.num_bonds(), 0);
    assert_eq!(xyz(&deleted, 1), [1., 3.]);
    let oxygen = AtomSpec::new(Element::O);
    assert_eq!(
        m.with_atom_replaced(a(0), &oxygen, false)
            .unwrap()
            .topology()
            .atoms[0]
            .element(),
        Element::O
    );
    let m = mol("[CH3][CH3]");
    let double = BondSpec::new(a(1), a(0), BondOrder::Double);
    let new = m.with_bond_replaced(b(0), &double, false, true).unwrap();
    assert_eq!(new.topology().atoms[0].explicit_hydrogens(), 2);
    assert_eq!(new.topology().bonds[0].begin(), a(0));
    let mut left = with_x("CC", &[(3, &[31., 32.]), (1, &[11., 12.])]);
    let right = with_x("O", &[(1, &[17.]), (3, &[37.])]);
    let value = left.with_molecule_appended(&right).unwrap();
    let orig = left.clone();
    assert!(left.insert_molecule_(&mol("O")).is_err());
    assert_eq!(left, orig);
    left.insert_molecule_(&right).unwrap();
    assert_eq!(left, value);
    assert_eq!(xyz(&left, 3), [31., 32., 37.]);
    assert_eq!(xyz(&left, 1), [11., 12., 17.]);
    let mut repeated = mol("CC");
    repeated.insert_molecule_(&repeated.clone()).unwrap();
    assert_eq!(repeated.num_atoms(), 4);
    assert_eq!(repeated.topology().bonds[1].begin(), a(2));
}

#[test]
fn subset_group_links_identify_the_same_surviving_members() {
    use cosmolkit::{SGroupAttachPoint, SubstanceGroup, SubstanceGroupId, SubstanceGroupKind};
    let m = mol("CCCC |Sg:n:0::ht:::,Sg:n:1::ht:::,Sg:n:2::ht:::,SgH:2:1|");
    let part = m
        .extract_subgraph(&SubgraphSelection::atoms([1, 2, 3].map(a).to_vec()))
        .unwrap();
    let groups = &part.molecule().topology().substance_groups;
    assert_eq!(groups.len(), 1);
    assert_eq!(groups[0].atoms(), [a(1)]);
    assert_eq!(groups[0].head_crossing_bonds(), [b(0)]);
    let m = mol("CCC |Sg:n:1::ht:::,Sg:n:2::ht:::,SgH:0:1|");
    assert!(
        m.extract_subgraph(&SubgraphSelection::bonds(vec![b(1)]))
            .unwrap()
            .molecule()
            .topology()
            .substance_groups
            .is_empty()
    );

    let m = mol("CCOC");
    let mut t = m.topology().clone();
    t.substance_groups.push(
        SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Superatom)
            .with_atoms(vec![a(2), a(3)])
            .with_attach_points(vec![SGroupAttachPoint {
                atom: a(2),
                leaving_atom: Some(a(1)),
                label: None,
                order: None,
            }]),
    );
    let m = Molecule::from_parts(t, Default::default(), Default::default()).unwrap();
    let part = m
        .extract_subgraph(&SubgraphSelection::atoms([1, 2, 3].map(a).to_vec()))
        .unwrap();
    let attachment = &part.molecule().topology().substance_groups[0].attach_points()[0];
    assert_eq!(
        (attachment.atom, attachment.leaving_atom),
        (a(1), Some(a(0)))
    );
}

#[test]
fn replacement_property_precedence_and_stereo_reference_validation() {
    use cosmolkit::{BondStereo, PropertyValue as V};
    let m = mol("F/C=C/Cl");
    let mut t = m.topology().clone();
    t.atoms[0].set_prop("shared", 1_i32).unwrap();
    t.bonds[1].set_prop("shared", 1_i32).unwrap();
    let mut m = Molecule::from_parts(t, Default::default(), Default::default()).unwrap();
    let atom = AtomSpec::new(Element::F)
        .with_prop("shared", 2_i32)
        .unwrap()
        .with_prop("new_only", 2_i32)
        .unwrap();
    let bond = BondSpec::new(a(0), a(1), BondOrder::Double)
        .with_stereo(BondStereo::Trans)
        .with_stereo_atoms(a(0), a(3))
        .with_prop("shared", 2_i32)
        .unwrap()
        .with_prop("new_only", 2_i32)
        .unwrap();
    for (preserve, expected) in [(false, 2), (true, 1)] {
        let result = m
            .with_atom_replaced(a(0), &atom, preserve)
            .unwrap()
            .with_bond_replaced(b(1), &bond, preserve, true)
            .unwrap();
        assert_eq!(
            result.topology().atoms[0].prop("shared"),
            Some(&V::Int(expected))
        );
        assert_eq!(
            result.topology().bonds[1].prop("shared"),
            Some(&V::Int(expected))
        );
        assert_eq!(
            result.topology().atoms[0].prop("new_only").is_some(),
            !preserve
        );
        assert_eq!(
            result.topology().bonds[1].prop("new_only").is_some(),
            !preserve
        );
        assert_eq!(
            result.topology().bonds[1].stereo_atoms(),
            Some([a(0), a(3)])
        );
    }
    let before = m.clone();
    for invalid in [
        bond.with_stereo_atoms(a(3), a(0)),
        BondSpec::new(a(1), a(2), BondOrder::Double).with_stereo(BondStereo::Cis),
    ] {
        assert!(m.replace_bond_(b(1), &invalid, false, true).is_err());
        assert_eq!(m, before);
    }
    let deleted = m.without_atom(a(0)).unwrap();
    assert_eq!(deleted.topology().bonds[0].stereo(), BondStereo::E);
    assert!(
        deleted.topology().bonds[0]
            .stereo_atom_references()
            .is_empty()
    );
}

#[test]
fn identity_permutations_and_subsets_have_different_property_rules() {
    use cosmolkit::{
        BondStereo, MoleculeProperties, PropertyValue as V, SubstanceGroup, SubstanceGroupId,
        SubstanceGroupKind,
    };
    let mut t = mol("F/C=C/Cl").topology().clone();
    t.substance_groups.push(
        SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
            .with_atoms(vec![a(1)]),
    );
    t.atoms[1].set_prop("tag", 11_i32).unwrap();
    t.bonds[1].set_prop("tag", 22_i32).unwrap();
    t.atoms[1].set_computed_prop("derived", 11_i32).unwrap();
    t.bonds[1].set_computed_prop("derived", 22_i32).unwrap();
    let mut props = MoleculeProperties::default();
    props.set_prop("ordinary", 1_i32).unwrap();
    props.set_computed_prop("derived", 33_i32).unwrap();
    let m = Molecule::from_parts(t, Default::default(), props).unwrap();
    for (order, local, refs) in [([0, 1, 2, 3], 1, [0, 3]), ([3, 2, 1, 0], 2, [3, 0])] {
        let n = m.with_atoms_renumbered(&order.map(a)).unwrap();
        assert!(n.properties().prop("ordinary").is_none());
        assert!(n.topology().substance_groups.is_empty());
        assert_eq!(n.topology().atoms[local].prop("tag"), Some(&V::Int(11)));
        assert_eq!(n.topology().bonds[1].prop("tag"), Some(&V::Int(22)));
        assert_eq!(n.topology().bonds[1].stereo_atoms(), Some(refs.map(a)));
    }
    let selection =
        SubgraphSelection::atoms_and_bonds([0, 1, 2, 3].map(a).to_vec(), [0, 1, 2].map(b).to_vec());
    for clear in [false, true] {
        let result = m
            .extract_subgraph_with_params(
                &selection,
                &SubgraphExtractionParams {
                    clear_computed_props: clear,
                    ..Default::default()
                },
            )
            .unwrap();
        let n = result.molecule();
        assert!(n.properties().prop("ordinary").is_none());
        assert_eq!(n.topology().substance_groups[0].atoms(), [a(1)]);
        assert!(n.topology().atoms[1].prop("derived").is_none());
        assert_eq!(n.topology().bonds[1].prop("derived").is_some(), !clear);
        assert_eq!(n.properties().prop("derived").is_some(), !clear);
        assert_eq!(n.topology().bonds[1].stereo(), BondStereo::Trans);
    }
    let m = mol("C[C@H](O)[C@@H](C)O |&1:1,3|");
    let p = m
        .extract_subgraph(&SubgraphSelection::atoms([0, 1, 2].map(a).to_vec()))
        .unwrap();
    assert_eq!(p.molecule().topology().stereo_groups[0].atoms(), [a(1)]);
}

#[test]
fn union_pairs_each_dimension_independently_and_requires_unique_keys() {
    use cosmolkit::Conformer2D;
    let left = mol("[CH4:7]");
    let right = mol("[OH2:7]");
    let coords = |x| CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(3, vec![[x, 0.]])],
        conformers_3d: vec![Conformer3D::new(3, vec![[x, 1., 2.]], true)],
        ..Default::default()
    };
    let mut left = Molecule::from_parts(
        left.topology().clone(),
        coords(11.),
        left.properties().clone(),
    )
    .unwrap();
    let right = Molecule::from_parts(
        right.topology().clone(),
        coords(22.),
        right.properties().clone(),
    )
    .unwrap();
    let empty = Molecule::default();
    assert_eq!(xyz(&left.with_molecule_appended(&empty).unwrap(), 3), [11.]);
    assert_eq!(
        xyz(&empty.with_molecule_appended(&right).unwrap(), 3),
        [22.]
    );
    let duplicate = with_x("O", &[(3, &[8.]), (3, &[9.])]);
    let mut duplicate_left = duplicate.clone();
    assert!(duplicate_left.insert_molecule_(&duplicate).is_err());
    assert_eq!(duplicate_left, duplicate);
    let peer = left.clone();
    left.insert_molecule_(&right).unwrap();
    assert_eq!(left.coordinates_2d().unwrap(), [[11., 0.], [22., 0.]]);
    assert_eq!(xyz(&left, 3), [11., 22.]);
    assert_eq!(
        left.topology().atoms[0].atom_map(),
        left.topology().atoms[1].atom_map()
    );
    assert_eq!(peer.num_atoms(), 1);
    assert_eq!(right.num_atoms(), 1);
}

#[test]
fn raw_valence_ranking_and_stored_hydrogen_queries_use_the_current_graph() {
    use cosmolkit::{AtomEnvironmentParams, SmilesParseParams, SubgraphEnumerationParams};
    let raw = SmilesParseParams {
        sanitize: false,
        remove_hs: false,
        ..Default::default()
    };
    let m = Molecule::from_smiles_with_params("CCO", &raw).unwrap();
    let peer = m.clone();
    assert_eq!(m.canonical_atom_ranks().unwrap(), [0, 2, 1]);
    assert_eq!(m, peer);
    let m = Molecule::from_smiles_with_params("[H]CC", &raw).unwrap();
    for (include, expected) in [(false, vec![b(1)]), (true, vec![b(0), b(1)])] {
        let p = SubgraphEnumerationParams {
            rooted_at_atom: Some(1),
            use_hydrogens: include,
        };
        assert_eq!(
            m.subgraphs_of_length_with_params(1, &p).unwrap(),
            expected.iter().map(|b| vec![*b]).collect::<Vec<_>>()
        );
        let p = AtomEnvironmentParams {
            use_hydrogens: include,
            enforce_radius: false,
        };
        assert_eq!(
            m.atom_environment_with_params(3, a(1), &p).unwrap().bonds(),
            expected
        );
    }
    assert!(
        m.paths_of_length_with_params(
            0,
            &PathEnumerationParams {
                rooted_at_atom: Some(3),
                ..Default::default()
            }
        )
        .is_err()
    );
}
