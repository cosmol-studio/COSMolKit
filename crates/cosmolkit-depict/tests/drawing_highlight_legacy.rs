//! Frozen source-first legacy String-highlight regression; no fixture generation.
use std::collections::BTreeMap;
use std::path::PathBuf;

use cosmolkit_core::{RingFindType, RingInfo, ValenceAssignment};
use cosmolkit_depict::{DepictOptions, DrawingInput, render_svg};
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, Conformer2D, CoordinateBlock,
    CoordinateDimension, Element, MoleculeProperties, NeighborRef, PropertyValue, TopologyBlock,
};

#[derive(Clone, Copy)]
struct Case {
    label: &'static str,
    atom: Option<&'static str>,
    bond: Option<&'static str>,
}

const CASES: [Case; 8] = [
    Case {
        label: "absent",
        atom: None,
        bond: None,
    },
    Case {
        label: "atom_one",
        atom: Some("1"),
        bond: None,
    },
    Case {
        label: "atom_true",
        atom: Some("true"),
        bond: None,
    },
    Case {
        label: "bond_empty",
        atom: None,
        bond: Some(""),
    },
    Case {
        label: "bond_zero",
        atom: None,
        bond: Some("0"),
    },
    Case {
        label: "bond_one",
        atom: None,
        bond: Some("1"),
    },
    Case {
        label: "bond_true",
        atom: None,
        bond: Some("true"),
    },
    Case {
        label: "bond_case",
        atom: None,
        bond: Some("True"),
    },
];

fn typed_highlight(
    key: &str,
    note: Option<&str>,
) -> BTreeMap<cosmolkit_model::PropertyText, PropertyValue> {
    note.into_iter()
        .map(|value| (key.into(), PropertyValue::String(value.into())))
        .collect()
}

fn coordinate_identity(coordinates: &CoordinateBlock) -> Vec<(usize, Vec<[u64; 2]>)> {
    coordinates
        .conformers_2d
        .iter()
        .map(|conformer| {
            (
                conformer.id(),
                conformer
                    .coordinates()
                    .iter()
                    .map(|p| [p[0].to_bits(), p[1].to_bits()])
                    .collect(),
            )
        })
        .collect()
}

#[test]
fn drawing_highlight_legacy_product() {
    let reference_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../testdata/depiction/expected/legacy_0_3_0/highlights");
    // Only reads the already frozen source bytes. Missing references are fatal.
    let references: Vec<_> = CASES
        .iter()
        .map(|case| {
            let path = reference_dir.join(format!("{}.svg", case.label));
            std::fs::read(&path)
                .unwrap_or_else(|error| panic!("missing reference {}: {error}", path.display()))
        })
        .collect();
    assert_eq!(references.len(), 8);
    for index in [1, 2, 3, 4, 7] {
        assert_eq!(
            references[index], references[0],
            "source default/case control"
        );
    }
    assert_eq!(references[5], references[6], "source one/true equality");
    assert_ne!(
        references[5], references[0],
        "source bond highlight differs"
    );
    let mut calls = 0;
    let mut outcomes = Vec::new();
    let mut preservation_failures = Vec::new();
    for (case, expected) in CASES.iter().zip(&references) {
        for repeat in 0..2 {
            // Every repetition constructs independent literal detached inputs.
            let mut atom_spec = AtomSpec::new(Element::C);
            if let Some(note) = case.atom {
                atom_spec = atom_spec
                    .with_prop("_highlight", PropertyValue::String(note.into()))
                    .unwrap();
            }
            let mut bond_spec = BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single);
            if let Some(note) = case.bond {
                bond_spec = bond_spec
                    .with_prop("_highlight", PropertyValue::String(note.into()))
                    .unwrap();
            }
            let topology = TopologyBlock::try_from_parts(
                vec![
                    Atom::from_spec(AtomId::new(0), atom_spec.clone()),
                    Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
                ],
                vec![Bond::from_spec(BondId::new(0), bond_spec.clone())],
                vec![],
                vec![],
            )
            .unwrap();
            let coordinates = CoordinateBlock {
                conformers_2d: vec![Conformer2D::new(0, vec![[0.0, 0.0], [1.5, 0.0]])],
                conformers_3d: vec![],
                source_coordinate_dim: Some(CoordinateDimension::TwoD),
                source_conformer_order: Some(vec![CoordinateDimension::TwoD]),
            };
            let properties = MoleculeProperties::default();
            let valence = ValenceAssignment {
                explicit_valence: vec![1, 1],
                implicit_hydrogens: vec![3, 3],
            };
            let rings = RingInfo::new(RingFindType::SymmSssr, 2, 1);
            let options = DepictOptions {
                width: 300,
                height: 300,
            };

            // Validate whole default specs, exact adjacency and ordinary property
            // presence/value variants before taking the immediate call baseline.
            topology.validate().unwrap();
            assert_eq!(
                topology.atoms,
                vec![
                    Atom::from_spec(AtomId::new(0), atom_spec),
                    Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
                ]
            );
            assert_eq!(
                topology.bonds,
                vec![Bond::from_spec(BondId::new(0), bond_spec)]
            );
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
            assert_eq!(
                topology.atoms[0].props(),
                &typed_highlight("_highlight", case.atom)
            );
            assert!(topology.atoms[1].props().is_empty());
            assert_eq!(
                topology.bonds[0].props(),
                &typed_highlight("_highlight", case.bond)
            );
            assert!(topology.atoms.iter().all(|atom| {
                atom.computed_prop_names()
                    .unwrap()
                    .is_none_or(<[_]>::is_empty)
            }));
            assert!(
                topology.bonds[0]
                    .computed_prop_names()
                    .unwrap()
                    .is_none_or(<[_]>::is_empty)
            );
            assert!(properties.props().is_empty());
            assert!(properties.prop("_highlight").is_none());
            assert!(properties.name().is_none());
            assert!(properties.sdf_data_fields().is_empty());
            assert!(properties.sdf_property_lists().is_empty());
            assert!(
                properties
                    .computed_prop_names()
                    .unwrap()
                    .is_none_or(<[_]>::is_empty)
            );
            coordinates.validate_for_atom_count(2).unwrap();
            assert_eq!(coordinates.conformers_2d.len(), 1);
            assert!(coordinates.conformers_3d.is_empty());
            assert_eq!(
                coordinates.source_coordinate_dim,
                Some(CoordinateDimension::TwoD)
            );
            assert!(coordinates.conformers_2d[0].props().is_empty());
            assert_eq!(
                coordinate_identity(&coordinates),
                vec![(0, vec![[0, 0], [4609434218613702656, 0]])]
            );
            assert_eq!(valence.explicit_valence, [1, 1]);
            assert_eq!(valence.implicit_hydrogens, [3, 3]);
            assert!(rings.is_initialized());
            assert_eq!(rings.atom_row_count(), 2);
            assert_eq!(rings.bond_row_count(), 1);
            assert_eq!(rings, RingInfo::new(RingFindType::SymmSssr, 2, 1));

            // Fresh complete five-input and separate bits/IDs checkpoint
            // immediately before the actual public call.
            let before = (
                topology.clone(),
                coordinates.clone(),
                properties.clone(),
                valence.clone(),
                rings.clone(),
            );
            let before_coordinate_identity = coordinate_identity(&coordinates);
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
            // Check every input immediately after Result, even on an error,
            // before inspecting Result or comparing any output bytes.
            let checks = [
                topology == before.0,
                coordinates == before.1,
                properties == before.2,
                valence == before.3,
                rings == before.4,
                coordinate_identity(&coordinates) == before_coordinate_identity,
            ];
            let preserved = checks.into_iter().all(|equal| equal);
            if !preserved {
                preservation_failures.push((case.label, repeat));
            }
            eprintln!(
                "highlight call={calls} label={} repeat={repeat} preserved={preserved}",
                case.label
            );
            outcomes.push((case.label, repeat, result, expected));
        }
    }
    assert_eq!(
        calls, 16,
        "actual render_svg calls counted after invocation"
    );
    assert_eq!(outcomes.len(), 16);
    let mut failures = Vec::new();
    for (label, repeat, result, expected) in outcomes {
        match result {
            Ok(svg) => {
                let actual = svg.as_slice();
                let equal = actual == expected.as_slice();
                eprintln!(
                    "highlight outcome label={label} repeat={repeat} expected_bytes={} actual_bytes={} exact={equal}",
                    expected.len(),
                    actual.len()
                );
                if !equal {
                    let first = expected
                        .iter()
                        .zip(actual)
                        .position(|(a, b)| a != b)
                        .unwrap_or(expected.len().min(actual.len()));
                    failures.push(format!("{label}/{repeat}: bytes differ at {first}; expected={} actual={}; byte expected={:?} actual={:?}",
                        expected.len(), actual.len(), expected.get(first), actual.get(first)));
                }
            }
            Err(error) => {
                eprintln!("highlight outcome label={label} repeat={repeat} error={error:?}");
                failures.push(format!("{label}/{repeat}: {error:?}"));
            }
        }
    }
    eprintln!(
        "highlight census calls={calls} preservation_failures={} byte_or_error_failures={}",
        preservation_failures.len(),
        failures.len()
    );
    assert!(
        preservation_failures.is_empty(),
        "five-input/coordinate preservation failures: {preservation_failures:?}"
    );
    assert!(
        failures.is_empty(),
        "raw legacy SVG mismatches/errors: {failures:#?}"
    );
}
