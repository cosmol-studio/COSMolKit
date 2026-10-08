use cosmolkit_core::{
    AtropisomerAssignment, AtropisomerBondUpdate, AtropisomerConformer, AtropisomerError,
    AtropisomerRejectionKind, StereoError, detect_atropisomer_chirality,
};
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, Conformer2D, Conformer3D,
    SourceAtomValenceFacts, TopologyBlock,
};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo, Element, Hybridization};

fn facts(explicit: i8, implicit: i8) -> SourceAtomValenceFacts {
    SourceAtomValenceFacts {
        explicit_valence: explicit,
        implicit_valence: implicit,
    }
}
fn topology(
    elements: &[Element],
    edges: &[(usize, usize, BondOrder, BondDirection)],
) -> TopologyBlock {
    let atoms = elements
        .iter()
        .enumerate()
        .map(|(i, &e)| Atom::from_spec(AtomId::new(i), AtomSpec::new(e)))
        .collect();
    let bonds = edges
        .iter()
        .enumerate()
        .map(|(i, &(a, b, o, d))| {
            Bond::from_spec(
                BondId::new(i),
                BondSpec::new(AtomId::new(a), AtomId::new(b), o).with_direction(d),
            )
        })
        .collect();
    TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
}
fn alkane(order: BondOrder, stereo: BondStereo) -> TopologyBlock {
    let mut t = topology(
        &[Element::C; 4],
        &[
            (1, 0, BondOrder::Single, BondDirection::BeginWedge),
            (1, 2, order, BondDirection::None),
            (2, 3, BondOrder::Single, BondDirection::BeginDash),
        ],
    );
    t.bonds[1].set_stereo(stereo).unwrap();
    t
}
// Exact CC(=O)C(=O)C graph, reordered detached IDs only to make the axis 1->2.
fn carbonyl_axis(with_unrelated_nitrogen: bool) -> TopologyBlock {
    let mut elements = vec![
        Element::C,
        Element::C,
        Element::C,
        Element::C,
        Element::O,
        Element::O,
    ];
    if with_unrelated_nitrogen {
        elements.push(Element::N);
    }
    topology(
        &elements,
        &[
            (1, 0, BondOrder::Single, BondDirection::BeginWedge),
            (1, 2, BondOrder::Single, BondDirection::None),
            (2, 3, BondOrder::Single, BondDirection::BeginDash),
            (1, 4, BondOrder::Double, BondDirection::None),
            (2, 5, BondOrder::Double, BondDirection::None),
        ],
    )
}
fn xyz(y: f64, z: f64) -> Conformer3D {
    Conformer3D::new(
        0,
        vec![
            [0.0, y, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.0, z],
            [0.0, -y, 0.0],
            [1.0, 0.0, -z],
        ],
        true,
    )
}
#[test]
fn no_directional_candidates_leave_all_uninitialized_cache_rows_untouched() {
    let mut t = carbonyl_axis(false);
    for b in &mut t.bonds {
        b.set_direction(BondDirection::None);
    }
    let before = t.clone();
    assert_eq!(
        detect_atropisomer_chirality(&t, None).unwrap(),
        AtropisomerAssignment::default()
    );
    assert_eq!(t, before);
}
#[test]
fn failed_total_degree_retains_only_lazy_endpoint_cache_writes_with_exact_implicit_counts() {
    let t = alkane(BondOrder::Single, BondStereo::None);
    let before = t.clone();
    let a = detect_atropisomer_chirality(&t, None).unwrap();
    assert_eq!(
        a.atom_valence_updates,
        vec![(AtomId::new(1), facts(2, 2)), (AtomId::new(2), facts(2, 2))]
    );
    assert!(a.bond_updates.is_empty() && a.diagnostics.is_empty());
    assert_eq!(a.conjugated_bonds, None);
    assert_eq!(a.hybridization, None);
    assert_eq!(t, before);
}
#[test]
fn lazy_endpoint_updates_precede_both_non_single_and_stereo_any_filters() {
    for (order, stereo, expected) in [
        (BondOrder::Double, BondStereo::None, facts(3, 1)),
        (BondOrder::Single, BondStereo::Any, facts(2, 2)),
    ] {
        let t = alkane(order, stereo);
        let before = t.clone();
        let a = detect_atropisomer_chirality(&t, None).unwrap();
        assert_eq!(
            a.atom_valence_updates,
            vec![(AtomId::new(1), expected), (AtomId::new(2), expected)]
        );
        assert!(a.bond_updates.is_empty() && a.diagnostics.is_empty());
        assert_eq!(a.conjugated_bonds, None);
        assert_eq!(a.hybridization, None);
        assert_eq!(t, before);
    }
}
#[test]
fn passing_axis_runs_source_nonlocal_cache_conjugation_and_hybridization_before_detection() {
    let t = carbonyl_axis(true);
    let before = t.clone();
    let a = detect_atropisomer_chirality(&t, None).unwrap();
    let all = [
        facts(1, 3),
        facts(4, 0),
        facts(4, 0),
        facts(1, 3),
        facts(2, 0),
        facts(2, 0),
        facts(0, 3),
    ];
    assert_eq!(a.atom_valence_updates.len(), 9);
    assert_eq!(
        &a.atom_valence_updates[2..],
        &all.iter()
            .copied()
            .enumerate()
            .map(|(i, f)| (AtomId::new(i), f))
            .collect::<Vec<_>>()
    );
    assert_eq!(
        a.conjugated_bonds,
        Some(vec![false, true, false, true, true])
    );
    assert_eq!(
        a.hybridization.unwrap().values,
        vec![
            Hybridization::Sp3,
            Hybridization::Sp2,
            Hybridization::Sp2,
            Hybridization::Sp3,
            Hybridization::Sp2,
            Hybridization::Sp2,
            Hybridization::Sp3
        ]
    );
    assert_eq!(
        a.bond_updates,
        vec![AtropisomerBondUpdate {
            bond: BondId::new(1),
            stereo: BondStereo::AtropCcw
        }]
    );
    assert_eq!(t, before);
}
#[test]
fn healthy_source_cache_and_defined_real_hybridizations_preserve_existing_conjugation_flags() {
    let mut t = carbonyl_axis(false);
    let cached = [
        facts(1, 3),
        facts(4, 0),
        facts(4, 0),
        facts(1, 3),
        facts(2, 0),
        facts(2, 0),
    ];
    for (i, atom) in t.atoms.iter_mut().enumerate() {
        atom.set_source_valence_facts(cached[i]);
        atom.set_hybridization(if i == 0 || i == 3 {
            Hybridization::Sp3
        } else {
            Hybridization::Sp2
        });
    }
    for b in &mut t.bonds {
        b.set_conjugated(true);
    }
    let before = t.clone();
    let a = detect_atropisomer_chirality(&t, None).unwrap();
    assert!(a.atom_valence_updates.is_empty());
    assert_eq!(a.conjugated_bonds, None);
    assert_eq!(a.hybridization, None);
    assert_eq!(
        a.bond_updates,
        vec![AtropisomerBondUpdate {
            bond: BondId::new(1),
            stereo: BondStereo::AtropCcw
        }]
    );
    assert_eq!(t, before);
}
#[test]
fn source_3d_projection_preserves_lengths_and_strict_cross_product_threshold() {
    let t = carbonyl_axis(false);
    let before = t.clone();
    let threshold = 1.0e-7_f64;
    for (y, z, stereo, rejection) in [
        (1.0, 1.0, Some(BondStereo::AtropCw), None),
        (
            1.0e-4,
            1.0e-4,
            None,
            Some(AtropisomerRejectionKind::CoplanarCarriers),
        ),
        (
            1.0,
            threshold,
            None,
            Some(AtropisomerRejectionKind::CoplanarCarriers),
        ),
        (
            1.0,
            f64::from_bits(threshold.to_bits() + 1),
            Some(BondStereo::AtropCw),
            None,
        ),
        (
            1.0,
            f64::from_bits(threshold.to_bits() - 1),
            None,
            Some(AtropisomerRejectionKind::CollinearCarrier),
        ),
    ] {
        let conf = xyz(y, z);
        let original_coords = conf.clone();
        let a =
            detect_atropisomer_chirality(&t, Some(AtropisomerConformer::ThreeD(&conf))).unwrap();
        match stereo {
            Some(stereo) => {
                assert_eq!(
                    a.bond_updates,
                    vec![AtropisomerBondUpdate {
                        bond: BondId::new(1),
                        stereo
                    }]
                );
                assert!(a.diagnostics.is_empty());
            }
            None => {
                assert!(a.bond_updates.is_empty());
                assert_eq!(a.diagnostics.len(), 1);
                assert_eq!(a.diagnostics[0].kind, rejection.unwrap());
            }
        }
        assert_eq!(conf, original_coords);
        assert_eq!(t, before);
    }
}
#[test]
fn source_normalization_exception_is_fatal_instead_of_a_coplanarity_fallback() {
    let t = carbonyl_axis(false);
    let before = t.clone();
    let conf = Conformer3D::new(
        0,
        vec![
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [f64::MAX, 0.0, 0.0],
            [f64::MAX, 0.0, 1.0],
            [0.0, -1.0, 0.0],
            [f64::MAX, 0.0, -1.0],
        ],
        true,
    );
    assert_eq!(
        detect_atropisomer_chirality(&t, Some(AtropisomerConformer::ThreeD(&conf))),
        Err(AtropisomerError::Normalization(
            StereoError::ZeroLengthVector {
                center: AtomId::new(1),
                neighbor: AtomId::new(2)
            }
        ))
    );
    assert_eq!(t, before);
}
#[test]
fn source_2d_direction_failure_precedes_same_side_projection_on_the_same_end() {
    let mut t = carbonyl_axis(false);
    t.bonds[3].set_direction(BondDirection::BeginWedge);
    let before = t.clone();
    let conf = Conformer2D::new(
        0,
        vec![
            [0.0, 1.0],
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, -1.0],
            [0.0, 2.0],
            [1.0, 1.0],
        ],
    );
    let a = detect_atropisomer_chirality(&t, Some(AtropisomerConformer::TwoD(&conf))).unwrap();
    assert!(a.bond_updates.is_empty());
    assert_eq!(a.diagnostics.len(), 1);
    assert_eq!(
        a.diagnostics[0].kind,
        AtropisomerRejectionKind::InconsistentDirections
    );
    assert_eq!(t, before);
}
