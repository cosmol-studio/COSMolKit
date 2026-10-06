//! Preparation regressions at the detached owner boundary.
use super::*;
use cosmolkit_core::{RingSearchParams, symmetrized_sssr};
use cosmolkit_model::{Atom, AtomSpec, Bond, BondSpec, Conformer2D, Conformer3D, Element};

// Literal inputs for the exact legacy draw.rs SMILES identities. Input atom/
// bond rows were transcribed using the existing main RDKit wheel; expectations
// below remain the original legacy/source literals, never generated CK output.
fn literal_topology(
    atoms: &[(Element, bool, u8, bool, ChiralTag)],
    bonds: &[(usize, usize, BondOrder, bool)],
) -> TopologyBlock {
    TopologyBlock::try_from_parts(
        atoms
            .iter()
            .enumerate()
            .map(|(i, &(element, aromatic, hs, no_implicit, tag))| {
                Atom::from_spec(
                    AtomId::new(i),
                    AtomSpec::new(element)
                        .with_aromatic(aromatic)
                        .with_explicit_hydrogens(hs)
                        .with_no_implicit(no_implicit)
                        .with_chiral_tag(tag),
                )
            })
            .collect(),
        bonds
            .iter()
            .enumerate()
            .map(|(i, &(begin, end, order, aromatic))| {
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), order)
                        .with_aromatic(aromatic),
                )
            })
            .collect(),
        vec![],
        vec![],
    )
    .unwrap()
}

fn prepared_preserving(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    rings: Option<&RingInfo>,
) -> PreparedDrawing {
    let properties = MoleculeProperties::default();
    let valence =
        assign_valence_with_options_for_topology(topology, ValenceModel::RdkitLike, false).unwrap();
    let before = (
        topology.clone(),
        coordinates.clone(),
        properties.clone(),
        valence.clone(),
        rings.cloned(),
    );
    let coordinate_bits = bits(coordinates);
    let outcome = prepare(DrawingInput {
        topology,
        coordinates,
        properties: &properties,
        valence: Some(&valence),
        rings,
    });
    assert_eq!(
        before,
        (
            topology.clone(),
            coordinates.clone(),
            properties,
            valence,
            rings.cloned()
        )
    );
    assert_eq!(coordinate_bits, bits(coordinates));
    outcome.unwrap()
}

fn bits(
    coordinates: &CoordinateBlock,
) -> (Vec<(usize, Vec<[u64; 2]>)>, Vec<(usize, Vec<[u64; 3]>)>) {
    (
        coordinates
            .conformers_2d
            .iter()
            .map(|c| {
                (
                    c.id(),
                    c.coordinates()
                        .iter()
                        .map(|p| p.map(f64::to_bits))
                        .collect(),
                )
            })
            .collect(),
        coordinates
            .conformers_3d
            .iter()
            .map(|c| {
                (
                    c.id(),
                    c.coordinates()
                        .iter()
                        .map(|p| p.map(f64::to_bits))
                        .collect(),
                )
            })
            .collect(),
    )
}

// CC1=C(CCNCCN)c2cc3[nH]c(cc4[nH]c(cc5nc(cc1n2)C(C)=C5CCNCCN)c(C)c4CCNCCN)c(CCNCCN)c3C
fn canonical_input() -> TopologyBlock {
    literal_topology(
        &[
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::N, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::N, false, 0, false, ChiralTag::Unspecified),
            (Element::C, true, 0, false, ChiralTag::Unspecified),
            (Element::C, true, 0, false, ChiralTag::Unspecified),
            (Element::C, true, 0, false, ChiralTag::Unspecified),
            (Element::N, true, 1, false, ChiralTag::Unspecified),
            (Element::C, true, 0, false, ChiralTag::Unspecified),
            (Element::C, true, 0, false, ChiralTag::Unspecified),
            (Element::C, true, 0, false, ChiralTag::Unspecified),
            (Element::N, true, 1, false, ChiralTag::Unspecified),
            (Element::C, true, 0, false, ChiralTag::Unspecified),
            (Element::C, true, 0, false, ChiralTag::Unspecified),
            (Element::C, true, 0, false, ChiralTag::Unspecified),
            (Element::N, true, 0, false, ChiralTag::Unspecified),
            (Element::C, true, 0, false, ChiralTag::Unspecified),
            (Element::C, true, 0, false, ChiralTag::Unspecified),
            (Element::C, true, 0, false, ChiralTag::Unspecified),
            (Element::N, true, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::N, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::N, false, 0, false, ChiralTag::Unspecified),
            (Element::C, true, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::C, true, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::N, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::N, false, 0, false, ChiralTag::Unspecified),
            (Element::C, true, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::N, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::N, false, 0, false, ChiralTag::Unspecified),
            (Element::C, true, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
        ],
        &[
            (0, 1, BondOrder::Single, false),
            (1, 2, BondOrder::Double, false),
            (2, 3, BondOrder::Single, false),
            (3, 4, BondOrder::Single, false),
            (4, 5, BondOrder::Single, false),
            (5, 6, BondOrder::Single, false),
            (6, 7, BondOrder::Single, false),
            (7, 8, BondOrder::Single, false),
            (2, 9, BondOrder::Single, false),
            (9, 10, BondOrder::Aromatic, true),
            (10, 11, BondOrder::Aromatic, true),
            (11, 12, BondOrder::Aromatic, true),
            (12, 13, BondOrder::Aromatic, true),
            (13, 14, BondOrder::Aromatic, true),
            (14, 15, BondOrder::Aromatic, true),
            (15, 16, BondOrder::Aromatic, true),
            (16, 17, BondOrder::Aromatic, true),
            (17, 18, BondOrder::Aromatic, true),
            (18, 19, BondOrder::Aromatic, true),
            (19, 20, BondOrder::Aromatic, true),
            (20, 21, BondOrder::Aromatic, true),
            (21, 22, BondOrder::Aromatic, true),
            (22, 23, BondOrder::Aromatic, true),
            (23, 24, BondOrder::Aromatic, true),
            (21, 25, BondOrder::Single, false),
            (25, 26, BondOrder::Single, false),
            (25, 27, BondOrder::Double, false),
            (27, 28, BondOrder::Single, false),
            (28, 29, BondOrder::Single, false),
            (29, 30, BondOrder::Single, false),
            (30, 31, BondOrder::Single, false),
            (31, 32, BondOrder::Single, false),
            (32, 33, BondOrder::Single, false),
            (17, 34, BondOrder::Aromatic, true),
            (34, 35, BondOrder::Single, false),
            (34, 36, BondOrder::Aromatic, true),
            (36, 37, BondOrder::Single, false),
            (37, 38, BondOrder::Single, false),
            (38, 39, BondOrder::Single, false),
            (39, 40, BondOrder::Single, false),
            (40, 41, BondOrder::Single, false),
            (41, 42, BondOrder::Single, false),
            (13, 43, BondOrder::Aromatic, true),
            (43, 44, BondOrder::Single, false),
            (44, 45, BondOrder::Single, false),
            (45, 46, BondOrder::Single, false),
            (46, 47, BondOrder::Single, false),
            (47, 48, BondOrder::Single, false),
            (48, 49, BondOrder::Single, false),
            (43, 50, BondOrder::Aromatic, true),
            (50, 51, BondOrder::Single, false),
            (23, 1, BondOrder::Single, false),
            (24, 9, BondOrder::Aromatic, true),
            (50, 11, BondOrder::Aromatic, true),
            (36, 15, BondOrder::Aromatic, true),
            (27, 19, BondOrder::Single, false),
        ],
    )
}

// O=C1CCC(=O)O[C@@H]2[C@@H](O)[C@H](O)[C@@H](COC(=O)CCC(=O)O[C@@H]3[C@@H](O)[C@H](O)[C@@H](CO1)O[C@@H]3O)O[C@@H]2O
fn polycycle_input() -> TopologyBlock {
    literal_topology(
        &[
            (Element::O, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::O, false, 0, false, ChiralTag::Unspecified),
            (Element::O, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 1, true, ChiralTag::TetrahedralCcw),
            (Element::C, false, 1, true, ChiralTag::TetrahedralCw),
            (Element::O, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 1, true, ChiralTag::TetrahedralCcw),
            (Element::O, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 1, true, ChiralTag::TetrahedralCw),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::O, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::O, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::O, false, 0, false, ChiralTag::Unspecified),
            (Element::O, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 1, true, ChiralTag::TetrahedralCcw),
            (Element::C, false, 1, true, ChiralTag::TetrahedralCw),
            (Element::O, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 1, true, ChiralTag::TetrahedralCcw),
            (Element::O, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 1, true, ChiralTag::TetrahedralCw),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::O, false, 0, false, ChiralTag::Unspecified),
            (Element::O, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 1, true, ChiralTag::TetrahedralCcw),
            (Element::O, false, 0, false, ChiralTag::Unspecified),
            (Element::O, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 1, true, ChiralTag::TetrahedralCcw),
            (Element::O, false, 0, false, ChiralTag::Unspecified),
        ],
        &[
            (0, 1, BondOrder::Double, false),
            (1, 2, BondOrder::Single, false),
            (2, 3, BondOrder::Single, false),
            (3, 4, BondOrder::Single, false),
            (4, 5, BondOrder::Double, false),
            (4, 6, BondOrder::Single, false),
            (6, 7, BondOrder::Single, false),
            (7, 8, BondOrder::Single, false),
            (8, 9, BondOrder::Single, false),
            (8, 10, BondOrder::Single, false),
            (10, 11, BondOrder::Single, false),
            (10, 12, BondOrder::Single, false),
            (12, 13, BondOrder::Single, false),
            (13, 14, BondOrder::Single, false),
            (14, 15, BondOrder::Single, false),
            (15, 16, BondOrder::Double, false),
            (15, 17, BondOrder::Single, false),
            (17, 18, BondOrder::Single, false),
            (18, 19, BondOrder::Single, false),
            (19, 20, BondOrder::Double, false),
            (19, 21, BondOrder::Single, false),
            (21, 22, BondOrder::Single, false),
            (22, 23, BondOrder::Single, false),
            (23, 24, BondOrder::Single, false),
            (23, 25, BondOrder::Single, false),
            (25, 26, BondOrder::Single, false),
            (25, 27, BondOrder::Single, false),
            (27, 28, BondOrder::Single, false),
            (28, 29, BondOrder::Single, false),
            (27, 30, BondOrder::Single, false),
            (30, 31, BondOrder::Single, false),
            (31, 32, BondOrder::Single, false),
            (12, 33, BondOrder::Single, false),
            (33, 34, BondOrder::Single, false),
            (34, 35, BondOrder::Single, false),
            (29, 1, BondOrder::Single, false),
            (34, 7, BondOrder::Single, false),
            (31, 22, BondOrder::Single, false),
        ],
    )
}

#[test]
fn drawing_prepare_canonical_kekule_legacy_regression() {
    let topology = canonical_input();
    let rings = symmetrized_sssr(&topology, &RingSearchParams::default()).unwrap();
    let prepared = prepared_preserving(&topology, &CoordinateBlock::default(), Some(&rings));
    let orders: String = prepared
        .topology
        .bonds
        .iter()
        .map(|b| match b.order() {
            BondOrder::Single => 'S',
            BondOrder::Double => 'D',
            BondOrder::Triple => 'T',
            BondOrder::Aromatic => 'A',
            _ => '?',
        })
        .collect();
    assert_eq!(
        orders,
        "SDSSSSSSSSDSSDSSSSDSDSDSSSDSSSSSSDSSSSSSSSSSSSSSSDSSDSDS"
    );
    let xy = prepared.coordinates.conformers_2d[0].coordinates();
    assert_eq!(xy[0][0].to_bits(), (-0.7271872994113224_f64).to_bits());
    assert_eq!(xy[0][1].to_bits(), (-7.508043887184295_f64).to_bits());
    assert_eq!(xy[49][0].to_bits(), 10.57352991454218_f64.to_bits());
    assert_eq!(xy[49][1].to_bits(), 4.749281830440951_f64.to_bits());
}

#[test]
fn drawing_prepare_polycyclic_chiral_hydrogen_legacy_regression() {
    let topology = polycycle_input();
    assert_eq!(topology.atoms.len(), 36);
    let rings = symmetrized_sssr(&topology, &RingSearchParams::default()).unwrap();
    let original_atom_rings = rings.atom_rings().to_vec();
    let original_bond_rings = rings.bond_rings().to_vec();
    // Both stored dimensions exercise actual AddHs coordinate append as well
    // as the caller's membership resizing. Retained rows/IDs are checked below.
    let coordinates = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(
            17,
            (0..36).map(|i| [i as f64, (i % 5) as f64]).collect(),
        )],
        conformers_3d: vec![Conformer3D::new(
            29,
            (0..36).map(|i| [i as f64, (i % 5) as f64, 2.0]).collect(),
            true,
        )],
        ..Default::default()
    };
    let prepared = prepared_preserving(&topology, &coordinates, Some(&rings));
    assert_eq!(prepared.topology.atoms.len(), 46);
    assert_eq!(
        prepared
            .topology
            .atoms
            .iter()
            .filter(|a| a.atomic_number() == 1)
            .count(),
        10
    );
    assert_eq!(prepared.rings.atom_row_count(), 46);
    assert_eq!(prepared.rings.bond_row_count(), topology.bonds.len() + 10);
    assert_eq!(prepared.rings.atom_rings(), original_atom_rings);
    assert_eq!(prepared.rings.bond_rings(), original_bond_rings);
    assert!(prepared.rings.is_symm_sssr());
    for i in 36..46 {
        assert_eq!(prepared.rings.num_atom_rings(AtomId::new(i)), 0);
    }
    for i in topology.bonds.len()..prepared.topology.bonds.len() {
        assert_eq!(prepared.rings.num_bond_rings(BondId::new(i)), 0);
    }
    assert_eq!(prepared.coordinates.conformers_2d[0].id(), 17);
    assert_eq!(prepared.coordinates.conformers_3d[0].id(), 29);
    assert_eq!(
        &prepared.coordinates.conformers_2d[0].coordinates()[..36],
        coordinates.conformers_2d[0].coordinates()
    );
    assert_eq!(
        &prepared.coordinates.conformers_3d[0].coordinates()[..36],
        coordinates.conformers_3d[0].coordinates()
    );
    prepared.coordinates.validate_for_atom_count(46).unwrap();
}

#[test]
fn drawing_prepare_aromatic_flags_and_unflagged_aromatic_orders() {
    for flag in [false, true] {
        let atoms = vec![(Element::C, true, 0, false, ChiralTag::Unspecified); 6];
        let bonds = (0..6)
            .map(|i| (i, (i + 1) % 6, BondOrder::Aromatic, flag))
            .collect::<Vec<_>>();
        let topology = literal_topology(&atoms, &bonds);
        let coordinates = CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(11, vec![[0.0, 0.0]; 6])],
            ..Default::default()
        };
        let prepared = prepared_preserving(&topology, &coordinates, None);
        assert!(prepared.topology.atoms.iter().all(|a| a.is_aromatic()));
        assert!(
            prepared
                .topology
                .bonds
                .iter()
                .all(|b| b.is_aromatic() == flag)
        );
        let orders = prepared
            .topology
            .bonds
            .iter()
            .map(|b| b.order())
            .collect::<Vec<_>>();
        // Default canonical=true, independently checked by existing-main
        // RDKit KekulizeIfPossible(MolFromSmiles("c1ccccc1"), false).
        // The initial test literal accidentally used noncanonical ordering;
        // its red log is retained. No legacy coordinate/fixture literal changes.
        assert_eq!(
            orders,
            if flag {
                vec![
                    BondOrder::Single,
                    BondOrder::Double,
                    BondOrder::Single,
                    BondOrder::Double,
                    BondOrder::Single,
                    BondOrder::Double,
                ]
            } else {
                vec![BondOrder::Aromatic; 6]
            }
        );
        assert_eq!(prepared.coordinates, coordinates);
    }
}

#[test]
fn drawing_prepare_no_2d_2d_and_3d_only_states() {
    let topology = literal_topology(
        &[
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::O, false, 0, false, ChiralTag::Unspecified),
        ],
        &[
            (0, 1, BondOrder::Single, false),
            (1, 2, BondOrder::Single, false),
        ],
    );
    let supplied = Conformer2D::new(41, vec![[0.0, -0.0], [1.5, 0.2], [3.0, 0.5]])
        .with_prop("layout", "supplied");
    let three_d = Conformer3D::new(
        83,
        vec![[0.0, 1.0, 2.0], [3.0, 4.0, 5.0], [6.0, 7.0, 8.0]],
        true,
    )
    .with_prop("3d", "retained");
    for stored_2d in [false, true] {
        for stored_3d in [false, true] {
            let coordinates = CoordinateBlock {
                conformers_2d: if stored_2d {
                    vec![supplied.clone(), supplied.clone().with_id(42)]
                } else {
                    vec![]
                },
                conformers_3d: if stored_3d {
                    vec![three_d.clone()]
                } else {
                    vec![]
                },
                ..Default::default()
            };
            let prepared = prepared_preserving(&topology, &coordinates, None);
            if stored_2d {
                assert_eq!(
                    prepared.coordinates.conformers_2d,
                    coordinates.conformers_2d
                );
            } else {
                let expected = compute_2d_coordinates(
                    &topology,
                    &MoleculeProperties::default(),
                    &Compute2DCoordinatesParams {
                        canonical_orientation: true,
                        ..Default::default()
                    },
                )
                .unwrap();
                assert_eq!(prepared.coordinates.conformers_2d, vec![expected]);
            }
            assert_eq!(
                prepared.coordinates.conformers_3d,
                coordinates.conformers_3d
            );
        }
    }
}

#[test]
fn drawing_prepare_wedges_use_first_supplied_2d_geometry() {
    let topology = literal_topology(
        &[
            (Element::C, false, 0, false, ChiralTag::TetrahedralCw),
            (Element::F, false, 0, false, ChiralTag::Unspecified),
            (Element::CL, false, 0, false, ChiralTag::Unspecified),
            (Element::BR, false, 0, false, ChiralTag::Unspecified),
        ],
        &[
            (0, 1, BondOrder::Single, false),
            (0, 2, BondOrder::Single, false),
            (0, 3, BondOrder::Single, false),
        ],
    );
    for xy in [
        vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [-1.0, 0.0]],
        vec![[0.0, 0.0], [1.0, 0.0], [0.0, -1.0], [-1.0, 0.0]],
    ] {
        let first = Conformer2D::new(17, xy);
        let second = Conformer2D::new(18, vec![[0.0, 0.0]; 4]);
        let coordinates = CoordinateBlock {
            conformers_2d: vec![first.clone(), second],
            ..Default::default()
        };
        let prepared = prepared_preserving(&topology, &coordinates, None);
        let expected = determine_bond_wedge_state(
            &topology,
            BondId::new(0),
            AtomId::new(0),
            Some(AtropisomerConformer::TwoD(&first)),
        )
        .unwrap();
        assert_eq!(prepared.topology.bonds[0].direction(), expected);
        assert_ne!(expected, BondDirection::None);
        assert_eq!(prepared.coordinates, coordinates);
    }
}

#[test]
fn drawing_prepare_typed_invalid_coordinates_and_memberships() {
    use std::error::Error;
    let topology = literal_topology(
        &[(Element::C, false, 0, false, ChiralTag::Unspecified)],
        &[],
    );
    let properties = MoleculeProperties::default();
    let invalid = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(9, vec![])],
        ..Default::default()
    };
    let error = prepare(DrawingInput {
        topology: &topology,
        coordinates: &invalid,
        properties: &properties,
        valence: None,
        rings: None,
    })
    .err()
    .unwrap();
    assert!(matches!(
        error,
        DrawingError::Coordinates(cosmolkit_model::CoordinateValidationError::RowCount {
            conformer: 9,
            rows: 0,
            atom_count: 1,
            ..
        })
    ));
    assert!(error.source().is_some());
    // A short (even empty) prefix is source-valid; only oversized storage is
    // malformed. Original zero-prefix fixture and its red result are retained
    // in the immutable delivery proposal for independent p1/ROOT review.
    let wrong = RingInfo::new(cosmolkit_core::RingFindType::Sssr, 2, 0);
    let error = prepare(DrawingInput {
        topology: &topology,
        coordinates: &CoordinateBlock::default(),
        properties: &properties,
        valence: None,
        rings: Some(&wrong),
    })
    .err()
    .unwrap();
    assert!(matches!(
        error,
        DrawingError::StateRows {
            field: "ring atom memberships",
            actual: 2,
            expected: 1
        }
    ));
    let wrong = RingInfo::new(cosmolkit_core::RingFindType::Sssr, 1, 1);
    let error = prepare(DrawingInput {
        topology: &topology,
        coordinates: &CoordinateBlock::default(),
        properties: &properties,
        valence: None,
        rings: Some(&wrong),
    })
    .err()
    .unwrap();
    assert!(matches!(
        error,
        DrawingError::StateRows {
            field: "ring bond memberships",
            actual: 1,
            expected: 0
        }
    ));
}

#[test]
fn drawing_prepare_sparse_ring_prefix_preserves_source_membership_reads() {
    // RDKit pin351f8f378f8ad6bbd517980c38896e66bf907af8:
    // RingInfo.cpp atomMembers/bondMembers return empty and num*Rings zero
    // beyond the cached prefix. AddHs preserves that old RingInfo carrier.
    let ethanol = literal_topology(
        &[
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::C, false, 0, false, ChiralTag::Unspecified),
            (Element::O, false, 0, false, ChiralTag::Unspecified),
        ],
        &[
            (0, 1, BondOrder::Single, false),
            (1, 2, BondOrder::Single, false),
        ],
    );
    let benzene = literal_topology(
        &[(Element::C, true, 0, false, ChiralTag::Unspecified); 6],
        &(0..6)
            .map(|i| (i, (i + 1) % 6, BondOrder::Aromatic, true))
            .collect::<Vec<_>>(),
    );
    for topology in [ethanol, benzene] {
        let original = symmetrized_sssr(&topology, &RingSearchParams::default()).unwrap();
        let expanded = add_hydrogens_with_params(
            topology.clone(),
            CoordinateBlock::default(),
            MoleculeProperties::default(),
            &AddHsParams::default(),
        )
        .unwrap();
        assert!(expanded.topology.atoms.len() > topology.atoms.len());
        let supplied = CoordinateBlock {
            conformers_2d: vec![
                compute_2d_coordinates(
                    &expanded.topology,
                    &MoleculeProperties::default(),
                    &Compute2DCoordinatesParams::default(),
                )
                .unwrap(),
            ],
            ..Default::default()
        };
        let prepared = prepared_preserving(&expanded.topology, &supplied, Some(&original));
        assert_eq!(prepared.rings.atom_rings(), original.atom_rings());
        assert_eq!(prepared.rings.bond_rings(), original.bond_rings());
        assert_eq!(prepared.rings.is_symm_sssr(), original.is_symm_sssr());
        for i in 0..expanded.topology.atoms.len() {
            assert_eq!(
                prepared.rings.atom_members(AtomId::new(i)),
                original.atom_members(AtomId::new(i))
            );
            assert_eq!(
                prepared.rings.num_atom_rings(AtomId::new(i)),
                original.num_atom_rings(AtomId::new(i))
            );
        }
        for i in 0..expanded.topology.bonds.len() {
            assert_eq!(
                prepared.rings.bond_members(BondId::new(i)),
                original.bond_members(BondId::new(i))
            );
            assert_eq!(
                prepared.rings.num_bond_rings(BondId::new(i)),
                original.num_bond_rings(BondId::new(i))
            );
        }
        // Borrowed original state and supplied coordinates remain untouched.
        assert_eq!(original.atom_row_count(), topology.atoms.len());
        assert_eq!(original.bond_row_count(), topology.bonds.len());
        assert_eq!(prepared.coordinates, supplied);
    }
    let topology = literal_topology(
        &[(Element::C, false, 0, false, ChiralTag::Unspecified)],
        &[],
    );
    let empty_prefix = RingInfo::new(cosmolkit_core::RingFindType::Sssr, 0, 0);
    let prepared = prepared_preserving(&topology, &CoordinateBlock::default(), Some(&empty_prefix));
    assert!(prepared.rings.atom_rings().is_empty());
    assert!(prepared.rings.bond_rings().is_empty());
    assert_eq!(prepared.rings.num_atom_rings(AtomId::new(0)), 0);
    assert_eq!(empty_prefix.atom_row_count(), 0);
}
