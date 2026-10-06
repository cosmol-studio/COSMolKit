use std::f64::consts::PI;

use cosmolkit_core::{
    AtomPositionParams, CanonicalTransformParams, CentroidParams, PrincipalAxesKind,
    PrincipalAxesParams, Transform3D, TransformError, angle_degrees, angle_radians, bond_length,
    canonical_transform, canonicalize_conformer, centroid, dihedral_degrees, dihedral_radians,
    principal_axes_and_moments, transform_conformer, with_angle_degrees, with_angle_radians,
    with_atom_position, with_bond_length, with_dihedral_degrees, with_dihedral_radians,
};
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, Conformer3D, CoordinateBlock,
    CoordinateValidationError, TopologyBlock, TopologyValidationError,
};
use cosmolkit_types::{BondOrder, Element};

fn topology(elements: &[Element], edges: &[(usize, usize)]) -> TopologyBlock {
    let atoms = elements
        .iter()
        .copied()
        .enumerate()
        .map(|(index, element)| Atom::from_spec(AtomId::new(index), AtomSpec::new(element)))
        .collect();
    let bonds = edges
        .iter()
        .copied()
        .enumerate()
        .map(|(index, (begin, end))| {
            Bond::from_spec(
                BondId::new(index),
                BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
            )
        })
        .collect();
    TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
}

fn carbons(count: usize, edges: &[(usize, usize)]) -> TopologyBlock {
    topology(&vec![Element::C; count], edges)
}

fn conformer(points: &[[f64; 3]]) -> Conformer3D {
    Conformer3D::new(0, points.to_vec(), true)
}

fn close(left: f64, right: f64) {
    assert!((left - right).abs() <= 1.0e-9, "{left} != {right}");
}

fn point_close(left: [f64; 3], right: [f64; 3]) {
    for axis in 0..3 {
        close(left[axis], right[axis]);
    }
}

#[test]
fn identity_and_single_atom_canonical_translation_preserve_row_major_application() {
    assert_eq!(Transform3D::identity().values()[0], 1.0);
    assert_eq!(Transform3D::identity().values()[5], 1.0);
    assert_eq!(Transform3D::identity().values()[10], 1.0);
    point_close(
        Transform3D::identity().transform_point([1.0, 2.0, 3.0]),
        [1.0, 2.0, 3.0],
    );

    let topology = carbons(1, &[]);
    let input = conformer(&[[2.0, -3.0, 4.0]]);
    let transform = canonical_transform(&topology, &input, &Default::default()).unwrap();
    point_close(transform.transform_point(input.coordinates()[0]), [0.0; 3]);
    assert_eq!(&transform.values()[12..], &[0.0, 0.0, 0.0, 1.0]);

    let centered = canonical_transform(
        &topology,
        &input,
        &CanonicalTransformParams {
            center: Some([1.0, 1.0, 1.0]),
            ..Default::default()
        },
    )
    .unwrap();
    point_close(centered.transform_point([1.0, 1.0, 1.0]), [0.0; 3]);
}

#[test]
fn centroid_covers_hydrogen_filter_weights_extra_rows_and_source_immutability() {
    let topology = topology(&[Element::C, Element::H, Element::O], &[]);
    let input = conformer(&[[0.0, 0.0, 0.0], [100.0, 0.0, 0.0], [4.0, 2.0, 0.0]]);
    let saved = input.clone();
    point_close(
        centroid(&topology, &input, &Default::default()).unwrap(),
        [2.0, 1.0, 0.0],
    );
    point_close(
        centroid(
            &topology,
            &input,
            &CentroidParams {
                ignore_hydrogens: false,
                weights: Some(vec![1.0, 0.0, 3.0, 99.0]),
            },
        )
        .unwrap(),
        [3.0, 1.5, 0.0],
    );
    assert_eq!(input, saved);
}

#[test]
fn centroid_errors_keep_weight_and_selection_fields_distinct() {
    let all_h = topology(&[Element::H], &[]);
    let one = conformer(&[[0.0, 0.0, 0.0]]);
    assert_eq!(
        centroid(&all_h, &one, &Default::default()),
        Err(TransformError::NoSelectedAtoms)
    );

    let carbon = carbons(2, &[]);
    let two = conformer(&[[0.0; 3], [1.0, 0.0, 0.0]]);
    assert_eq!(
        centroid(
            &carbon,
            &two,
            &CentroidParams {
                ignore_hydrogens: false,
                weights: Some(vec![1.0])
            },
        ),
        Err(TransformError::WeightCount {
            actual: 1,
            required: 2
        })
    );
    match centroid(
        &carbon,
        &two,
        &CentroidParams {
            ignore_hydrogens: false,
            weights: Some(vec![1.0, f64::NAN]),
        },
    ) {
        Err(TransformError::NonFiniteWeight { atom, value }) => {
            assert_eq!(atom, AtomId::new(1));
            assert!(value.is_nan());
        }
        other => panic!("expected NonFiniteWeight for atom 1 with a NaN value, got {other:?}"),
    }
    assert_eq!(
        centroid(
            &carbon,
            &two,
            &CentroidParams {
                ignore_hydrogens: false,
                weights: Some(vec![1.0, -1.0])
            },
        ),
        Err(TransformError::InvalidWeightSum { sum: 0.0 })
    );
}

#[test]
fn inertia_and_gyration_principal_moments_are_sorted_and_normalized() {
    let topology = carbons(3, &[]);
    let input = conformer(&[[-1.0, 0.0, 0.0], [0.0; 3], [1.0, 0.0, 0.0]]);
    let inertia = principal_axes_and_moments(&topology, &input, &Default::default()).unwrap();
    assert_eq!(inertia.moments, [0.0, 2.0, 2.0]);
    let gyration = principal_axes_and_moments(
        &topology,
        &input,
        &PrincipalAxesParams {
            kind: PrincipalAxesKind::Gyration,
            ..Default::default()
        },
    )
    .unwrap();
    close(gyration.moments[0], 0.0);
    close(gyration.moments[1], 0.0);
    close(gyration.moments[2], 2.0 / 3.0);
    for column in 0..3 {
        close(
            (0..3).map(|row| gyration.axes[row][column].powi(2)).sum(),
            1.0,
        );
    }
}

#[test]
fn canonicalization_centers_rotates_preserves_distances_and_metadata() {
    let topology = carbons(3, &[]);
    let input = conformer(&[[1.0, 2.0, 0.0], [2.0, 3.0, 0.0], [4.0, 6.0, 0.0]])
        .with_id(7)
        .with_prop("source", "kept");
    let output = canonicalize_conformer(&topology, &input, &Default::default()).unwrap();
    assert_eq!(output.id(), 7);
    assert!(output.is_3d());
    assert_eq!(
        output.props().get("source").map(String::as_str),
        Some("kept")
    );
    point_close(
        centroid(&topology, &output, &Default::default()).unwrap(),
        [0.0; 3],
    );
    for left in 0..3 {
        for right in 0..3 {
            close(
                bond_length(&input, AtomId::new(left), AtomId::new(right)).unwrap(),
                bond_length(&output, AtomId::new(left), AtomId::new(right)).unwrap(),
            );
        }
    }
}

#[test]
fn transform_conformer_rejects_nonfinite_input_without_touching_metadata() {
    let input = Conformer3D::new(4, vec![[f64::NAN, 0.0, 0.0]], false).with_prop("p", "v");
    assert_eq!(
        transform_conformer(&input, &Transform3D::identity()),
        Err(TransformError::NonFinitePoint {
            role: "conformer",
            axis: "x"
        })
    );
}

#[test]
fn atom_position_selects_first_3d_or_exact_stable_id_and_preserves_other_conformers() {
    let topology = carbons(2, &[]);
    let coordinates = CoordinateBlock {
        conformers_3d: vec![
            Conformer3D::new(3, vec![[0.0; 3], [1.0, 0.0, 0.0]], false),
            Conformer3D::new(9, vec![[2.0; 3], [3.0; 3]], true).with_prop("p", "v"),
        ],
        ..Default::default()
    };
    let default_changed = with_atom_position(
        &topology,
        &coordinates,
        AtomId::new(1),
        [8.0, 7.0, 6.0],
        &Default::default(),
    )
    .unwrap();
    assert_eq!(
        default_changed.conformers_3d[0],
        coordinates.conformers_3d[0]
    );
    assert_eq!(
        default_changed.conformers_3d[1].coordinates()[1],
        [8.0, 7.0, 6.0]
    );
    assert_eq!(coordinates.conformers_3d[1].coordinates()[1], [3.0; 3]);

    let exact = with_atom_position(
        &topology,
        &coordinates,
        AtomId::new(0),
        [5.0; 3],
        &AtomPositionParams {
            conformer_id: Some(3),
        },
    )
    .unwrap();
    assert_eq!(exact.conformers_3d[0].coordinates()[0], [5.0; 3]);
    assert_eq!(exact.conformers_3d[1], coordinates.conformers_3d[1]);
}

#[test]
fn atom_position_errors_preserve_id_axis_and_coordinate_validation() {
    let topology = carbons(1, &[]);
    let empty = CoordinateBlock::default();
    assert_eq!(
        with_atom_position(
            &topology,
            &empty,
            AtomId::new(0),
            [0.0; 3],
            &Default::default()
        ),
        Err(TransformError::No3dConformer)
    );
    let valid = CoordinateBlock {
        conformers_3d: vec![conformer(&[[0.0; 3]])],
        ..Default::default()
    };
    assert_eq!(
        with_atom_position(
            &topology,
            &valid,
            AtomId::new(0),
            [0.0, f64::INFINITY, 0.0],
            &Default::default(),
        ),
        Err(TransformError::NonFinitePoint {
            role: "position",
            axis: "y"
        })
    );
    assert_eq!(
        with_atom_position(
            &topology,
            &valid,
            AtomId::new(0),
            [0.0; 3],
            &AtomPositionParams {
                conformer_id: Some(44)
            },
        ),
        Err(TransformError::ConformerNotFound { conformer_id: 44 })
    );
}

#[test]
fn bond_length_setter_moves_only_the_selected_side_and_roundtrips() {
    let topology = carbons(4, &[(0, 1), (1, 2), (2, 3)]);
    let input = conformer(&[
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [3.0, 1.0, 0.0],
    ]);
    let output = with_bond_length(&topology, &input, AtomId::new(1), AtomId::new(2), 2.0).unwrap();
    close(
        bond_length(&output, AtomId::new(1), AtomId::new(2)).unwrap(),
        2.0,
    );
    assert_eq!(output.coordinates()[0], input.coordinates()[0]);
    assert_eq!(output.coordinates()[1], input.coordinates()[1]);
    point_close(output.coordinates()[2], [3.0, 0.0, 0.0]);
    point_close(output.coordinates()[3], [4.0, 1.0, 0.0]);
    assert_eq!(input.coordinates()[2], [2.0, 0.0, 0.0]);
}

#[test]
fn bond_setter_rejects_nonbond_ring_coincident_and_invalid_target() {
    let chain = carbons(3, &[(0, 1), (1, 2)]);
    let input = conformer(&[[0.0; 3], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);
    assert!(matches!(
        with_bond_length(&chain, &input, AtomId::new(0), AtomId::new(2), 1.0),
        Err(TransformError::AtomsNotBonded {
            first_role: "i",
            second_role: "j",
            ..
        })
    ));
    assert!(matches!(
        with_bond_length(&chain, &input, AtomId::new(1), AtomId::new(2), 1.0),
        Err(TransformError::CoincidentCoordinates {
            first_role: "i",
            second_role: "j",
            ..
        })
    ));
    assert!(matches!(
        with_bond_length(&chain, &input, AtomId::new(0), AtomId::new(1), -1.0),
        Err(TransformError::InvalidTargetValue {
            quantity: "bond_length",
            ..
        })
    ));
    let ring = carbons(3, &[(0, 1), (1, 2), (2, 0)]);
    assert_eq!(
        with_bond_length(
            &ring,
            &conformer(&[[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
            AtomId::new(0),
            AtomId::new(1),
            2.0
        ),
        Err(TransformError::RingBondNotMovable {
            bond: BondId::new(0)
        })
    );
}

#[test]
fn angle_radian_degree_setters_roundtrip_and_move_k_side_only() {
    let topology = carbons(4, &[(0, 1), (1, 2), (2, 3)]);
    let input = conformer(&[[1.0, 0.0, 0.0], [0.0; 3], [0.0, 1.0, 0.0], [0.0, 2.0, 0.0]]);
    close(
        angle_radians(&input, AtomId::new(0), AtomId::new(1), AtomId::new(2)).unwrap(),
        PI / 2.0,
    );
    close(
        angle_degrees(&input, AtomId::new(0), AtomId::new(1), AtomId::new(2)).unwrap(),
        90.0,
    );
    let output = with_angle_degrees(
        &topology,
        &input,
        AtomId::new(0),
        AtomId::new(1),
        AtomId::new(2),
        60.0,
    )
    .unwrap();
    close(
        angle_degrees(&output, AtomId::new(0), AtomId::new(1), AtomId::new(2)).unwrap(),
        60.0,
    );
    assert_eq!(output.coordinates()[0], input.coordinates()[0]);
    assert_eq!(output.coordinates()[1], input.coordinates()[1]);
    assert_ne!(output.coordinates()[2], input.coordinates()[2]);
    assert_ne!(output.coordinates()[3], input.coordinates()[3]);
    let radians = with_angle_radians(
        &topology,
        &input,
        AtomId::new(0),
        AtomId::new(1),
        AtomId::new(2),
        PI / 3.0,
    )
    .unwrap();
    point_close(radians.coordinates()[2], output.coordinates()[2]);
}

#[test]
fn angle_setter_rejects_both_ring_bonds_and_undefined_axis() {
    let ring = carbons(3, &[(0, 1), (1, 2), (2, 0)]);
    let bent = conformer(&[[1.0, 0.0, 0.0], [0.0; 3], [0.0, 1.0, 0.0]]);
    assert_eq!(
        with_angle_radians(
            &ring,
            &bent,
            AtomId::new(0),
            AtomId::new(1),
            AtomId::new(2),
            1.0
        ),
        Err(TransformError::BothAngleBondsInRing {
            first_bond: BondId::new(0),
            second_bond: BondId::new(1)
        })
    );
    let chain = carbons(3, &[(0, 1), (1, 2)]);
    let collinear = conformer(&[[-1.0, 0.0, 0.0], [0.0; 3], [1.0, 0.0, 0.0]]);
    assert_eq!(
        with_angle_radians(
            &chain,
            &collinear,
            AtomId::new(0),
            AtomId::new(1),
            AtomId::new(2),
            1.0
        ),
        Err(TransformError::UndefinedRotationAxis)
    );
}

#[test]
fn signed_dihedral_degree_and_radian_setters_roundtrip_with_deterministic_side() {
    let topology = carbons(5, &[(0, 1), (1, 2), (2, 3), (3, 4)]);
    let input = conformer(&[
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 1.0],
        [1.0, 2.0, 1.0],
    ]);
    close(
        dihedral_degrees(
            &input,
            AtomId::new(0),
            AtomId::new(1),
            AtomId::new(2),
            AtomId::new(3),
        )
        .unwrap(),
        dihedral_radians(
            &input,
            AtomId::new(0),
            AtomId::new(1),
            AtomId::new(2),
            AtomId::new(3),
        )
        .unwrap()
            * 180.0
            / PI,
    );
    let output = with_dihedral_degrees(
        &topology,
        &input,
        AtomId::new(0),
        AtomId::new(1),
        AtomId::new(2),
        AtomId::new(3),
        -45.0,
    )
    .unwrap();
    close(
        dihedral_degrees(
            &output,
            AtomId::new(0),
            AtomId::new(1),
            AtomId::new(2),
            AtomId::new(3),
        )
        .unwrap(),
        -45.0,
    );
    assert_eq!(output.coordinates()[0], input.coordinates()[0]);
    assert_eq!(output.coordinates()[1], input.coordinates()[1]);
    assert_ne!(output.coordinates()[3], input.coordinates()[3]);
    assert_ne!(output.coordinates()[4], input.coordinates()[4]);
    let radians = with_dihedral_radians(
        &topology,
        &input,
        AtomId::new(0),
        AtomId::new(1),
        AtomId::new(2),
        AtomId::new(3),
        -PI / 4.0,
    )
    .unwrap();
    point_close(radians.coordinates()[3], output.coordinates()[3]);
}

#[test]
fn dihedral_rejects_ring_central_bond_coincident_and_collinear_planes() {
    let ring = carbons(4, &[(0, 1), (1, 2), (2, 3), (3, 0)]);
    let points = conformer(&[[0.0; 3], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 1.0]]);
    assert!(matches!(
        with_dihedral_radians(
            &ring,
            &points,
            AtomId::new(0),
            AtomId::new(1),
            AtomId::new(2),
            AtomId::new(3),
            0.0
        ),
        Err(TransformError::RingBondNotMovable { .. })
    ));
    let collinear = conformer(&[[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]]);
    assert_eq!(
        dihedral_radians(
            &collinear,
            AtomId::new(0),
            AtomId::new(1),
            AtomId::new(2),
            AtomId::new(3)
        ),
        Err(TransformError::UndefinedRotationAxis)
    );
}

#[test]
fn range_topology_and_coordinate_errors_keep_their_typed_boundaries() {
    let one = conformer(&[[0.0; 3]]);
    assert_eq!(
        angle_radians(&one, AtomId::new(0), AtomId::new(1), AtomId::new(2)),
        Err(TransformError::AtomOutOfRange {
            role: "j",
            atom: AtomId::new(1),
            atom_count: 1
        })
    );

    let invalid_topology = TopologyBlock {
        atoms: vec![Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C))],
        ..Default::default()
    };
    assert!(matches!(
        centroid(&invalid_topology, &one, &Default::default()),
        Err(TransformError::InvalidTopology(
            TopologyValidationError::AtomIdMismatch { .. }
        ))
    ));

    let topology = carbons(2, &[]);
    assert_eq!(
        centroid(&topology, &one, &Default::default()),
        Err(TransformError::InvalidCoordinates(
            CoordinateValidationError::RowCount {
                dimension: "3D",
                conformer: 0,
                rows: 1,
                atom_count: 2,
            }
        ))
    );
}
