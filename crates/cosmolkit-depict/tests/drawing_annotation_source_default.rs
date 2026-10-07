//! EXCLUDED additive pinned351f default-condition proposal; original legacy test stays mandatory.
use std::collections::BTreeMap;

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
    molecule: Option<&'static str>,
}

const CASES: [Case; 8] = [
    Case {
        label: "absent",
        atom: None,
        bond: None,
        molecule: None,
    },
    Case {
        label: "empty_atom",
        atom: Some(""),
        bond: None,
        molecule: None,
    },
    Case {
        label: "atom_only",
        atom: Some("atom_A"),
        bond: None,
        molecule: None,
    },
    Case {
        label: "empty_bond",
        atom: None,
        bond: Some(""),
        molecule: None,
    },
    Case {
        label: "bond_only",
        atom: None,
        bond: Some("bond_B"),
        molecule: None,
    },
    Case {
        label: "empty_molecule",
        atom: None,
        bond: None,
        molecule: Some(""),
    },
    Case {
        label: "molecule_only",
        atom: None,
        bond: None,
        molecule: Some("mol_C"),
    },
    Case {
        label: "combined",
        atom: Some("atom_A"),
        bond: Some("bond_B"),
        molecule: Some("mol_C"),
    },
];

fn typed_note(
    key: &str,
    note: Option<&str>,
) -> BTreeMap<cosmolkit_model::PropertyText, PropertyValue> {
    note.into_iter()
        .map(|value| (key.into(), PropertyValue::String(value.to_owned().into())))
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
fn drawing_annotation_pinned_source_default_product() {
    // Source-condition PROPOSAL only: exact native bytes at pinned351f, no oracle calls.
    let references: Vec<Vec<u8>> = vec![
        br###"<?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='300px' height='300px' viewBox='0 0 300 300'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='300.0' height='300.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 48.8,150.0 L 251.2,150.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
</svg>
"###.to_vec(),
        br###"<?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='300px' height='300px' viewBox='0 0 300 300'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='300.0' height='300.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 48.8,150.0 L 251.2,150.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
</svg>
"###.to_vec(),
        br###"<?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='300px' height='300px' viewBox='0 0 300 300'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='300.0' height='300.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 63.2,170.1 L 268.4,170.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<text x='31.6' y='145.9' class='note' style='font-size:20px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#000000' >a</text>
<text x='40.8' y='145.9' class='note' style='font-size:20px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#000000' >t</text>
<text x='45.4' y='145.9' class='note' style='font-size:20px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#000000' >o</text>
<text x='54.6' y='145.9' class='note' style='font-size:20px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#000000' >m</text>
<text x='68.4' y='145.9' class='note' style='font-size:20px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#000000' >_</text>
<text x='77.6' y='145.9' class='note' style='font-size:20px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#000000' >A</text>
</svg>
"###.to_vec(),
        br###"<?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='300px' height='300px' viewBox='0 0 300 300'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='300.0' height='300.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 48.8,150.0 L 251.2,150.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
</svg>
"###.to_vec(),
        br###"<?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='300px' height='300px' viewBox='0 0 300 300'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='300.0' height='300.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 15.0,166.5 L 285.0,166.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<text x='110.6' y='149.5' class='note' style='font-size:20px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#7F7FFF' >b</text>
<text x='122.1' y='149.5' class='note' style='font-size:20px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#7F7FFF' >o</text>
<text x='133.6' y='149.5' class='note' style='font-size:20px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#7F7FFF' >n</text>
<text x='145.1' y='149.5' class='note' style='font-size:20px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#7F7FFF' >d</text>
<text x='156.6' y='149.5' class='note' style='font-size:20px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#7F7FFF' >_</text>
<text x='168.1' y='149.5' class='note' style='font-size:20px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#7F7FFF' >B</text>
</svg>
"###.to_vec(),
        br###"<?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='300px' height='300px' viewBox='0 0 300 300'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='300.0' height='300.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 48.8,150.0 L 251.2,150.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
</svg>
"###.to_vec(),
        br###"<?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='300px' height='300px' viewBox='0 0 300 300'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='300.0' height='300.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 48.8,150.0 L 251.2,150.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<text x='195.4' y='52.0' class='note' style='font-size:40px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#000000' >m</text>
<text x='223.0' y='52.0' class='note' style='font-size:40px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#000000' >o</text>
<text x='241.4' y='52.0' class='note' style='font-size:40px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#000000' >l</text>
<text x='248.8' y='52.0' class='note' style='font-size:40px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#000000' >_</text>
<text x='267.2' y='52.0' class='note' style='font-size:40px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#000000' >C</text>
</svg>
"###.to_vec(),
        br###"<?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='300px' height='300px' viewBox='0 0 300 300'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='300.0' height='300.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 63.2,170.1 L 268.4,170.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<text x='31.6' y='145.9' class='note' style='font-size:20px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#000000' >a</text>
<text x='40.8' y='145.9' class='note' style='font-size:20px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#000000' >t</text>
<text x='45.4' y='145.9' class='note' style='font-size:20px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#000000' >o</text>
<text x='54.6' y='145.9' class='note' style='font-size:20px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#000000' >m</text>
<text x='68.4' y='145.9' class='note' style='font-size:20px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#000000' >_</text>
<text x='77.6' y='145.9' class='note' style='font-size:20px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#000000' >A</text>
<text x='126.4' y='159.6' class='note' style='font-size:20px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#7F7FFF' >b</text>
<text x='137.9' y='159.6' class='note' style='font-size:20px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#7F7FFF' >o</text>
<text x='149.4' y='159.6' class='note' style='font-size:20px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#7F7FFF' >n</text>
<text x='160.9' y='159.6' class='note' style='font-size:20px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#7F7FFF' >d</text>
<text x='172.4' y='159.6' class='note' style='font-size:20px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#7F7FFF' >_</text>
<text x='183.9' y='159.6' class='note' style='font-size:20px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#7F7FFF' >B</text>
<text x='195.4' y='52.0' class='note' style='font-size:40px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#000000' >m</text>
<text x='223.0' y='52.0' class='note' style='font-size:40px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#000000' >o</text>
<text x='241.4' y='52.0' class='note' style='font-size:40px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#000000' >l</text>
<text x='248.8' y='52.0' class='note' style='font-size:40px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#000000' >_</text>
<text x='267.2' y='52.0' class='note' style='font-size:40px;font-style:normal;font-weight:normal;fill-opacity:1;stroke:none;font-family:sans-serif;text-anchor:start;fill:#000000' >C</text>
</svg>
"###.to_vec(),
    ];
    assert_eq!(references.len(), 8);
    for index in [1, 3, 5] {
        assert_eq!(
            references[index], references[0],
            "frozen source empty-note control"
        );
    }
    for index in [2, 4, 6, 7] {
        assert_ne!(
            references[index], references[0],
            "frozen source nonempty control"
        );
    }
    let mut calls = 0;
    let mut outcomes = Vec::new();
    let mut preservation_failures = Vec::new();
    for (case, expected) in CASES.iter().zip(&references) {
        for repeat in 0..2 {
            // Every repetition constructs independent literal detached inputs.
            let mut atom_spec = AtomSpec::new(Element::C);
            if let Some(note) = case.atom {
                atom_spec = atom_spec
                    .with_prop("atomNote", PropertyValue::String(note.to_owned().into()))
                    .unwrap();
            }
            let mut bond_spec = BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single);
            if let Some(note) = case.bond {
                bond_spec = bond_spec
                    .with_prop("bondNote", PropertyValue::String(note.to_owned().into()))
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
                source_conformer_order: Some(vec![cosmolkit_model::CoordinateDimension::TwoD]),
            };
            let mut properties = MoleculeProperties::default();
            if let Some(note) = case.molecule {
                properties = properties.with_prop("molNote", note).unwrap();
            }
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
                &typed_note("atomNote", case.atom)
            );
            assert!(topology.atoms[1].props().is_empty());
            assert_eq!(
                topology.bonds[0].props(),
                &typed_note("bondNote", case.bond)
            );
            assert!(topology.atoms.iter().all(|atom| {
                atom.computed_prop_names()
                    .unwrap()
                    .is_none_or(|names| names.is_empty())
            }));
            assert!(
                topology.bonds[0]
                    .computed_prop_names()
                    .unwrap()
                    .is_none_or(|names| names.is_empty())
            );
            assert_eq!(
                properties.props(),
                &case
                    .molecule
                    .into_iter()
                    .map(|note| (
                        cosmolkit_model::PropertyText::from("molNote"),
                        PropertyValue::String(note.into())
                    ))
                    .collect::<BTreeMap<_, _>>()
            );
            assert_eq!(
                properties.prop("molNote"),
                case.molecule
                    .map(|value| PropertyValue::String(value.into()))
                    .as_ref()
            );
            assert!(properties.prop("atomNote").is_none());
            assert!(properties.name().is_none());
            assert!(properties.sdf_data_fields().is_empty());
            assert!(properties.sdf_property_lists().is_empty());
            assert!(
                properties
                    .computed_prop_names()
                    .unwrap()
                    .is_none_or(|names| names.is_empty())
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
                "annotation call={calls} label={} repeat={repeat} preserved={preserved}",
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
                    "annotation outcome label={label} repeat={repeat} expected_bytes={} actual_bytes={} exact={equal}",
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
                eprintln!("annotation outcome label={label} repeat={repeat} error={error:?}");
                failures.push(format!("{label}/{repeat}: {error:?}"));
            }
        }
    }
    eprintln!(
        "annotation census calls={calls} preservation_failures={} byte_or_error_failures={}",
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
