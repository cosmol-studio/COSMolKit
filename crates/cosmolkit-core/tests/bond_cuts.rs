use cosmolkit_core::{BondCut, BondCutError, cut_selected_bonds, get_molecule_fragments};
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondDirection, BondId, BondOrder, BondQueryPredicate, BondSpec,
    BondStereo, ChiralTag, Conformer2D, Conformer3D, CoordinateBlock, CoordinateDimension, Element,
    MoleculeProperties, QueryNode, TopologyBlock,
};

fn graph(n: usize, edges: &[(usize, usize)]) -> TopologyBlock {
    let atoms = (0..n)
        .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
        .collect();
    let bonds = edges
        .iter()
        .enumerate()
        .map(|(i, &(a, b))| {
            Bond::from_spec(
                BondId::new(i),
                BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
            )
        })
        .collect();
    TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
}

fn at(index: usize) -> BondCut {
    BondCut::at_bond(BondId::new(index), None, None).unwrap()
}

fn labeled(index: usize, begin: usize, end: usize) -> BondCut {
    BondCut::at_bond(BondId::new(index), Some(begin), Some(end)).unwrap()
}

#[test]
fn attachment_side_labels_and_all_conformer_rows_follow_opposite_source_atoms() {
    let source = graph(3, &[(0, 1), (1, 2)]);
    let coords = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(7, vec![[1.; 2], [2.; 2], [3.; 2]])],
        conformers_3d: vec![Conformer3D::new(9, vec![[8.; 3], [9.; 3], [10.; 3]], true)],
        source_conformer_order: Some(vec![CoordinateDimension::ThreeD, CoordinateDimension::TwoD]),
        ..Default::default()
    };
    let props = MoleculeProperties::default()
        .with_name("sample")
        .with_computed_prop("discard", 8)
        .unwrap();
    let manual = labeled(1, 101, 202);
    let reversed = BondCut::from_matched_ends(
        BondId::new(1),
        [(AtomId::new(2), 202), (AtomId::new(1), 101)],
    );
    for cut in [manual, reversed] {
        let r = cut_selected_bonds(&source, &coords, &props, &[cut]).unwrap();
        assert_eq!(r.topology.atoms[3].isotope(), Some(202));
        assert_eq!(r.topology.atoms[4].isotope(), Some(101));
        assert_eq!(
            r.topology
                .bonds
                .iter()
                .map(|b| (b.begin().index(), b.end().index()))
                .collect::<Vec<_>>(),
            [(0, 1), (3, 2), (1, 4)]
        );
        assert_eq!(
            r.coordinates.conformers_2d[0].coordinates(),
            &[[1.; 2], [2.; 2], [3.; 2], [2.; 2], [3.; 2]]
        );
        assert_eq!(
            r.coordinates.conformers_3d[0].coordinates(),
            &[[8.; 3], [9.; 3], [10.; 3], [9.; 3], [10.; 3]]
        );
        assert_eq!(
            r.coordinates.source_conformer_order,
            coords.source_conformer_order
        );
        assert_eq!(r.properties.name(), props.name());
        assert!(r.properties.prop("discard").is_none());
        assert_eq!(r.mapping.bonds.old_to_new, [Some(BondId::new(0)), None]);
        assert_eq!(
            r.mapping.bonds.new_to_old,
            [Some(BondId::new(0)), None, None]
        );
        assert_eq!(
            r.mapping.atoms.new_to_old,
            [
                Some(AtomId::new(0)),
                Some(AtomId::new(1)),
                Some(AtomId::new(2)),
                None,
                None
            ]
        );
        r.topology.validate().unwrap();
    }
}

#[test]
fn isotope_range_applies_only_to_effective_labels_in_the_selection() {
    let source = graph(65_537, &[(0, 1), (65_535, 65_536), (0, 65_535)]);
    let coords = CoordinateBlock::default();
    let props = MoleculeProperties::default();
    let cut = |cuts: &[BondCut]| cut_selected_bonds(&source, &coords, &props, cuts);
    let error = |cuts: &[BondCut]| cut(cuts).err().expect("selection must be rejected");
    assert_eq!(error(&[]), BondCutError::EmptySelection);
    for (begin, end) in [(Some(0), None), (None, Some(0))] {
        assert_eq!(
            BondCut::at_bond(BondId::new(0), begin, end),
            Err(BondCutError::IncompleteLabels {
                bond: BondId::new(0)
            })
        );
    }
    assert_eq!(
        error(&[at(0), at(0)]),
        BondCutError::Duplicate {
            bond: BondId::new(0)
        }
    );
    assert_eq!(
        error(&[at(0), at(3)]),
        BondCutError::BondOutOfRange {
            bond: BondId::new(3),
            count: 3
        }
    );
    assert!(matches!(
        error(&[BondCut::from_matched_ends(
            BondId::new(0),
            [(AtomId::new(0), 8), (AtomId::new(2), 9)]
        )]),
        BondCutError::EndpointsChanged { .. }
    ));
    assert_eq!(
        error(&[at(0), at(1)]),
        BondCutError::IsotopeOutOfRange {
            bond: BondId::new(1),
            attachment: AtomId::new(65_535),
            label: 65_536
        }
    );
    assert_eq!(
        error(&[at(0), labeled(1, 0, 65_536)]),
        BondCutError::IsotopeOutOfRange {
            bond: BondId::new(1),
            attachment: AtomId::new(65_536),
            label: 65_536
        }
    );

    let r = cut(&[at(0), labeled(1, 65_535, 0), at(2)]).unwrap();
    let labels = r.topology.atoms[65_537..]
        .iter()
        .map(|a| a.isotope().unwrap_or(0))
        .collect::<Vec<_>>();
    assert_eq!(labels, [0, 1, 0, 65_535, 0, 65_535]);
}

#[test]
fn one_or_two_ring_cuts_preserve_component_occurrences() {
    let ring = graph(6, &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]);
    for (cuts, expected) in [(vec![at(0)], vec![8]), (vec![at(0), at(3)], vec![5, 5])] {
        let r = cut_selected_bonds(&ring, &Default::default(), &Default::default(), &cuts).unwrap();
        let sizes = get_molecule_fragments(&r.topology, &r.coordinates, &r.properties, false, true)
            .unwrap()
            .into_iter()
            .map(|f| f.component_atoms().len())
            .collect::<Vec<_>>();
        assert_eq!(sizes, expected);
    }
}

#[test]
fn repeated_cuts_at_a_tetrahedral_center_track_incident_order() {
    let mut source = graph(5, &[(0, 1), (1, 2), (1, 3), (1, 4)]);
    for (i, e) in [Element::N, Element::C, Element::F, Element::CL, Element::BR]
        .into_iter()
        .enumerate()
    {
        source.atoms[i] = Atom::from_spec(AtomId::new(i), AtomSpec::new(e));
    }
    source.atoms[1].set_chiral_tag(ChiralTag::TetrahedralCcw);
    for (cuts, tag) in [
        ([at(0), at(2)], ChiralTag::TetrahedralCw),
        ([at(2), at(0)], ChiralTag::TetrahedralCcw),
    ] {
        let r =
            cut_selected_bonds(&source, &Default::default(), &Default::default(), &cuts).unwrap();
        assert_eq!(r.topology.atoms[1].chiral_tag(), tag);
    }
}

#[test]
fn retain_neighbor_double_bond_stereo_and_source_directions() {
    let mut source = graph(4, &[(0, 1), (1, 2), (2, 3)]);
    let a = AtomId::new;
    source.bonds[1] = Bond::from_spec(
        BondId::new(1),
        BondSpec::new(a(1), a(2), BondOrder::Double)
            .with_stereo_atoms(a(0), a(3))
            .with_stereo(BondStereo::E),
    );
    for i in [0, 2] {
        source.bonds[i].set_direction(BondDirection::EndUpRight);
    }
    let r = cut_selected_bonds(
        &source,
        &Default::default(),
        &Default::default(),
        &[at(0), at(2)],
    )
    .unwrap();
    assert_eq!(r.topology.bonds[0].stereo(), BondStereo::E);
    assert_eq!(r.topology.bonds[0].stereo_atom_references(), &[a(4), a(7)]);
    assert!(
        r.topology.bonds[1..]
            .iter()
            .all(|b| b.direction() == BondDirection::EndUpRight)
    );
}

#[test]
fn copy_query_to_both_bonds_but_wedge_only_to_original_begin_side() {
    let mut source = graph(2, &[(0, 1)]);
    let query = QueryNode::predicate(BondQueryPredicate::Any);
    source.bonds[0] = Bond::from_spec(
        BondId::new(0),
        BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
            .with_query(query.clone())
            .with_direction(BondDirection::BeginWedge),
    );
    let r =
        cut_selected_bonds(&source, &Default::default(), &Default::default(), &[at(0)]).unwrap();
    assert!(r.topology.bonds.iter().all(|b| b.query() == Some(&query)));
    assert_eq!(r.topology.bonds[0].direction(), BondDirection::None);
    assert_eq!(r.topology.bonds[1].direction(), BondDirection::BeginWedge);
}
