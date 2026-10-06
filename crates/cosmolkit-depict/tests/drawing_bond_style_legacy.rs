//! Frozen source-first bond-style regression. Tests only read fixed references.
use cosmolkit_core::{RingFindType, RingInfo, ValenceAssignment};
use cosmolkit_depict::{DepictOptions, DrawingInput, render_svg};
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, Conformer2D, CoordinateBlock,
    CoordinateDimension, Element, MoleculeProperties, NeighborRef, TopologyBlock,
};
use std::path::PathBuf;

const CASES: [(&str, BondOrder, usize, usize, [i32; 2], [i32; 2]); 12] = [
    ("triple_fwd", BondOrder::Triple, 0, 1, [3, 3], [1, 1]),
    ("triple_rev", BondOrder::Triple, 1, 0, [3, 3], [1, 1]),
    ("quadruple_fwd", BondOrder::Quadruple, 0, 1, [4, 4], [0, 0]),
    ("quadruple_rev", BondOrder::Quadruple, 1, 0, [4, 4], [0, 0]),
    ("hydrogen_fwd", BondOrder::Hydrogen, 0, 1, [0, 0], [4, 4]),
    ("hydrogen_rev", BondOrder::Hydrogen, 1, 0, [0, 0], [4, 4]),
    (
        "unspecified_fwd",
        BondOrder::Unspecified,
        0,
        1,
        [0, 0],
        [4, 4],
    ),
    (
        "unspecified_rev",
        BondOrder::Unspecified,
        1,
        0,
        [0, 0],
        [4, 4],
    ),
    ("dative_fwd", BondOrder::Dative, 0, 1, [0, 1], [4, 3]),
    ("dative_rev", BondOrder::Dative, 1, 0, [1, 0], [3, 4]),
    ("dative_one_fwd", BondOrder::DativeOne, 0, 1, [0, 1], [4, 3]),
    ("dative_one_rev", BondOrder::DativeOne, 1, 0, [1, 0], [3, 4]),
];

fn coordinate_identity(coordinates: &CoordinateBlock) -> Vec<(usize, Vec<[u64; 2]>)> {
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
        .collect()
}

#[test]
fn drawing_bond_style_legacy_product() {
    let reference_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../testdata/depiction/expected/legacy_0_3_0/bond_styles");
    let references: Vec<_> = CASES
        .iter()
        .map(|case| {
            let path = reference_dir.join(format!("{}.svg", case.0));
            std::fs::read(&path)
                .unwrap_or_else(|error| panic!("missing reference {}: {error}", path.display()))
        })
        .collect();
    assert_eq!(references.len(), 12);
    for (a, b) in [(4, 6), (5, 7), (8, 10), (9, 11)] {
        assert_eq!(references[a], references[b], "same-direction source alias");
    }
    let mut calls = 0;
    let mut outcomes = Vec::new();
    let mut preservation_failures = Vec::new();
    for (&(label, order, begin, end, explicit, implicit), expected) in CASES.iter().zip(&references)
    {
        for repeat in 0..2 {
            let topology = TopologyBlock::try_from_parts(
                vec![
                    Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
                    Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
                ],
                vec![Bond::from_spec(
                    BondId::new(0),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), order),
                )],
                vec![],
                vec![],
            )
            .unwrap();
            let coordinates = CoordinateBlock {
                conformers_2d: vec![Conformer2D::new(0, vec![[0.0, 0.0], [1.5, 0.0]])],
                conformers_3d: vec![],
                source_coordinate_dim: Some(CoordinateDimension::TwoD),
                source_conformer_order: None,
            };
            let properties = MoleculeProperties::default();
            let valence = ValenceAssignment {
                explicit_valence: explicit.to_vec(),
                implicit_hydrogens: implicit.to_vec(),
            };
            let rings = RingInfo::new(RingFindType::SymmSssr, 2, 1);
            let options = DepictOptions {
                width: 300,
                height: 300,
            };
            // Complete literal default specs include absence of queries, notes,
            // charges, isotopes, stereo and every ordinary/computed property.
            topology.validate().unwrap();
            assert_eq!(
                topology.atoms,
                vec![
                    Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
                    Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C))
                ]
            );
            assert_eq!(
                topology.bonds,
                vec![Bond::from_spec(
                    BondId::new(0),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), order)
                )]
            );
            assert_eq!(topology.atoms[0].id(), AtomId::new(0));
            assert_eq!(topology.atoms[1].id(), AtomId::new(1));
            assert_eq!(topology.bonds[0].id(), BondId::new(0));
            assert_eq!(topology.bonds[0].begin(), AtomId::new(begin));
            assert_eq!(topology.bonds[0].end(), AtomId::new(end));
            assert_eq!(topology.bonds[0].order(), order);
            assert_eq!(
                topology.adjacency.neighbors_of(0),
                &[NeighborRef {
                    atom_index: 1,
                    bond: BondId::new(0)
                }]
            );
            assert_eq!(
                topology.adjacency.neighbors_of(1),
                &[NeighborRef {
                    atom_index: 0,
                    bond: BondId::new(0)
                }]
            );
            assert!(topology.substance_groups.is_empty());
            assert!(topology.stereo_groups.is_empty());
            for atom in &topology.atoms {
                assert!(atom.props().is_empty());
                assert!(atom.computed_prop_names().is_empty());
            }
            assert!(topology.bonds[0].props().is_empty());
            assert!(topology.bonds[0].computed_prop_names().is_empty());
            assert_eq!(properties, MoleculeProperties::default());
            assert!(properties.props().is_empty());
            assert!(properties.computed_prop_names().is_empty());
            assert!(properties.name().is_none());
            assert!(properties.sdf_data_fields().is_empty());
            assert!(properties.sdf_property_lists().is_empty());
            coordinates.validate_for_atom_count(2).unwrap();
            assert_eq!(coordinates.conformers_2d.len(), 1);
            assert_eq!(
                coordinates.conformers_2d[0].coordinates(),
                &[[0.0, 0.0], [1.5, 0.0]]
            );
            assert!(coordinates.conformers_2d[0].props().is_empty());
            assert!(coordinates.conformers_3d.is_empty());
            assert_eq!(
                coordinates.source_coordinate_dim,
                Some(CoordinateDimension::TwoD)
            );
            assert_eq!(
                coordinate_identity(&coordinates),
                vec![(0, vec![[0, 0], [4609434218613702656, 0]])]
            );
            assert_eq!(valence.explicit_valence, explicit);
            assert_eq!(valence.implicit_hydrogens, implicit);
            assert!(rings.is_initialized());
            assert_eq!(rings.find_type(), RingFindType::SymmSssr);
            assert_eq!(rings.atom_row_count(), 2);
            assert_eq!(rings.bond_row_count(), 1);
            assert!(rings.atom_rings().is_empty());
            assert!(rings.bond_rings().is_empty());
            for id in [AtomId::new(0), AtomId::new(1)] {
                assert_eq!(rings.num_atom_rings(id), 0);
                assert!(rings.atom_members(id).is_empty());
            }
            assert_eq!(rings.num_bond_rings(BondId::new(0)), 0);
            assert!(rings.bond_members(BondId::new(0)).is_empty());
            assert_eq!(rings, RingInfo::new(RingFindType::SymmSssr, 2, 1));

            // Capture fresh whole FIVE inputs and separate coordinate IDs/bits
            // immediately before every actual public call, including repeats.
            let before = (
                topology.clone(),
                coordinates.clone(),
                properties.clone(),
                valence.clone(),
                rings.clone(),
            );
            let before_identity = coordinate_identity(&coordinates);
            let result = render_svg(
                DrawingInput {
                    topology: &topology,
                    coordinates: &coordinates,
                    properties: &properties,
                    valence: Some(&valence),
                    rings: Some(&rings),
                },
                &options,
            );
            calls += 1;
            // Evaluate ALL six checkpoints immediately after Result and BEFORE
            // output/error handling. No short-circuit hides an input comparison.
            let checks = [
                topology == before.0,
                coordinates == before.1,
                properties == before.2,
                valence == before.3,
                rings == before.4,
                coordinate_identity(&coordinates) == before_identity,
            ];
            let preserved = checks.into_iter().all(|equal| equal);
            if !preserved {
                preservation_failures.push((label, repeat, checks));
            }
            eprintln!(
                "bond_style call={calls} label={label} repeat={repeat} preserved={preserved} checks={checks:?}"
            );
            outcomes.push((label, repeat, result, expected));
        }
    }
    assert_eq!(calls, 24, "actual calls counted after invocation");
    assert_eq!(outcomes.len(), 24);
    let mut failures = Vec::new();
    for (label, repeat, result, expected) in outcomes {
        match result {
            Ok(svg) => {
                let actual = svg.as_bytes();
                let equal = actual == expected.as_slice();
                eprintln!(
                    "bond_style outcome label={label} repeat={repeat} expected_bytes={} actual_bytes={} exact={equal}",
                    expected.len(),
                    actual.len()
                );
                if !equal {
                    let first = expected
                        .iter()
                        .zip(actual)
                        .position(|(a, b)| a != b)
                        .unwrap_or(expected.len().min(actual.len()));
                    failures.push(format!("{label}/{repeat}: first={first} expected_len={} actual_len={} expected_byte={:?} actual_byte={:?}",expected.len(),actual.len(),expected.get(first),actual.get(first)));
                }
            }
            Err(error) => {
                eprintln!("bond_style outcome label={label} repeat={repeat} error={error:?}");
                failures.push(format!("{label}/{repeat}: {error:?}"));
            }
        }
    }
    eprintln!(
        "bond_style census calls={calls} preservation_failures={} byte_or_error_failures={}",
        preservation_failures.len(),
        failures.len()
    );
    assert!(
        preservation_failures.is_empty(),
        "five-input/coordinate preservation: {preservation_failures:#?}"
    );
    assert!(
        failures.is_empty(),
        "raw legacy SVG differences/errors: {failures:#?}"
    );
}
