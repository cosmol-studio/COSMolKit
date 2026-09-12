use cosmolkit_model::{
    Conformer2D, Conformer3D, CoordinateBlock, CoordinateDimension, CoordinateValidationError,
};

fn bits_2d(rows: &[[f64; 2]]) -> Vec<[u64; 2]> {
    rows.iter()
        .map(|row| [row[0].to_bits(), row[1].to_bits()])
        .collect()
}

fn bits_3d(rows: &[[f64; 3]]) -> Vec<[u64; 3]> {
    rows.iter()
        .map(|row| [row[0].to_bits(), row[1].to_bits(), row[2].to_bits()])
        .collect()
}

#[test]
fn coordinate_dimensions_and_empty_blocks_are_valid() {
    assert_ne!(CoordinateDimension::TwoD, CoordinateDimension::ThreeD);

    let block = CoordinateBlock::default();
    assert!(block.conformers_2d.is_empty());
    assert!(block.conformers_3d.is_empty());
    assert_eq!(block.source_coordinate_dim, None);
    assert_eq!(block.validate_for_atom_count(0), Ok(()));
    assert_eq!(block.validate_for_atom_count(17), Ok(()));

    let zero_atom_block = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(5, Vec::new())],
        conformers_3d: vec![Conformer3D::new(5, Vec::new(), false)],
        source_coordinate_dim: Some(CoordinateDimension::ThreeD),
    };
    assert_eq!(zero_atom_block.validate_for_atom_count(0), Ok(()));
}

#[test]
fn conformer_2d_preserves_bits_properties_ids_and_mutation() {
    let original = [[-0.0, f64::MAX], [f64::MIN, 0.0]];
    let mut conformer = Conformer2D::new(usize::MAX, original.to_vec())
        .with_prop("label", "first")
        .with_prop("stable", "yes")
        .with_prop("label", "replacement")
        .with_id(23);

    assert_eq!(conformer.id(), 23);
    assert_eq!(bits_2d(conformer.coordinates()), bits_2d(&original));
    assert_eq!(
        conformer.props().get("label").map(String::as_str),
        Some("replacement")
    );
    assert_eq!(
        conformer.props().get("stable").map(String::as_str),
        Some("yes")
    );

    conformer.coordinates_mut()[1][0] = -1.25;
    assert_eq!(conformer.coordinates().len(), 2);
    assert_eq!(
        conformer.coordinates()[1][0].to_bits(),
        (-1.25_f64).to_bits()
    );
    assert_eq!(
        conformer.coordinates()[0][0].to_bits(),
        (-0.0_f64).to_bits()
    );
    assert_eq!(conformer.id(), 23);
    assert_eq!(
        conformer.props().get("stable").map(String::as_str),
        Some("yes")
    );
    assert_eq!(conformer.validate_for_atom_count(2), Ok(()));
}

#[test]
fn conformer_3d_preserves_bits_flags_properties_ids_and_mutation() {
    let original = [[-0.0, f64::MAX, f64::MIN], [0.0, 2.0, -3.0]];
    for is_3d in [false, true] {
        let mut conformer = Conformer3D::new(91, original.to_vec(), is_3d)
            .with_prop("name", "before")
            .with_prop("keep", "value")
            .with_prop("name", "after")
            .with_id(7);

        assert_eq!(conformer.id(), 7);
        assert_eq!(conformer.is_3d(), is_3d);
        assert_eq!(bits_3d(conformer.coordinates()), bits_3d(&original));
        assert_eq!(
            conformer.props().get("name").map(String::as_str),
            Some("after")
        );
        assert_eq!(
            conformer.props().get("keep").map(String::as_str),
            Some("value")
        );

        conformer.coordinates_mut()[1][2] = 4.5;
        assert_eq!(conformer.coordinates().len(), 2);
        assert_eq!(conformer.coordinates()[1][2].to_bits(), 4.5_f64.to_bits());
        assert_eq!(
            conformer.coordinates()[0][0].to_bits(),
            (-0.0_f64).to_bits()
        );
        assert_eq!(conformer.id(), 7);
        assert_eq!(conformer.is_3d(), is_3d);
        assert_eq!(conformer.validate_for_atom_count(2), Ok(()));
    }
}

#[test]
fn row_count_errors_are_exact_for_both_dimensions() {
    assert_eq!(
        Conformer2D::new(12, vec![[0.0, 1.0]]).validate_for_atom_count(3),
        Err(CoordinateValidationError::RowCount {
            dimension: "2D",
            conformer: 12,
            rows: 1,
            atom_count: 3,
        })
    );
    assert_eq!(
        Conformer3D::new(19, vec![[0.0, 1.0, 2.0], [3.0, 4.0, 5.0]], true)
            .validate_for_atom_count(1),
        Err(CoordinateValidationError::RowCount {
            dimension: "3D",
            conformer: 19,
            rows: 2,
            atom_count: 1,
        })
    );
}

#[test]
fn duplicate_ids_are_dimension_local_and_exact() {
    let duplicate_2d = CoordinateBlock {
        conformers_2d: vec![
            Conformer2D::new(8, vec![[0.0, 0.0]]),
            Conformer2D::new(8, vec![[1.0, 1.0]]),
        ],
        ..Default::default()
    };
    assert_eq!(
        duplicate_2d.validate_for_atom_count(1),
        Err(CoordinateValidationError::DuplicateConformerId {
            dimension: "2D",
            id: 8,
        })
    );

    let duplicate_3d = CoordinateBlock {
        conformers_3d: vec![
            Conformer3D::new(13, vec![[0.0, 0.0, 0.0]], true),
            Conformer3D::new(13, vec![[1.0, 1.0, 1.0]], false),
        ],
        ..Default::default()
    };
    assert_eq!(
        duplicate_3d.validate_for_atom_count(1),
        Err(CoordinateValidationError::DuplicateConformerId {
            dimension: "3D",
            id: 13,
        })
    );

    let same_id_across_dimensions = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(21, vec![[0.0, 0.0]])],
        conformers_3d: vec![Conformer3D::new(21, vec![[0.0, 0.0, 0.0]], true)],
        source_coordinate_dim: Some(CoordinateDimension::TwoD),
    };
    assert_eq!(same_id_across_dimensions.validate_for_atom_count(1), Ok(()));
}

#[test]
fn nonfinite_2d_axes_report_exact_fields_for_all_classes() {
    for (axis_index, axis) in [(0, "x"), (1, "y")] {
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut rows = vec![[0.0, 0.0], [1.0, 2.0]];
            rows[1][axis_index] = value;
            assert_eq!(
                Conformer2D::new(31, rows).validate_for_atom_count(2),
                Err(CoordinateValidationError::NonFiniteCoordinate {
                    dimension: "2D",
                    conformer: 31,
                    atom: 1,
                    axis,
                })
            );
        }
    }
}

#[test]
fn nonfinite_3d_axes_report_exact_fields_for_all_classes() {
    for (axis_index, axis) in [(0, "x"), (1, "y"), (2, "z")] {
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut rows = vec![[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]];
            rows[1][axis_index] = value;
            assert_eq!(
                Conformer3D::new(37, rows, true).validate_for_atom_count(2),
                Err(CoordinateValidationError::NonFiniteCoordinate {
                    dimension: "3D",
                    conformer: 37,
                    atom: 1,
                    axis,
                })
            );
        }
    }
}

#[test]
fn validation_order_is_rows_then_axes_then_conformers_then_dimensions() {
    assert_eq!(
        Conformer2D::new(40, vec![[f64::NAN, 0.0]]).validate_for_atom_count(2),
        Err(CoordinateValidationError::RowCount {
            dimension: "2D",
            conformer: 40,
            rows: 1,
            atom_count: 2,
        })
    );

    assert_eq!(
        Conformer3D::new(
            41,
            vec![[0.0, f64::INFINITY, 0.0], [f64::NAN, 0.0, 0.0]],
            true,
        )
        .validate_for_atom_count(2),
        Err(CoordinateValidationError::NonFiniteCoordinate {
            dimension: "3D",
            conformer: 41,
            atom: 0,
            axis: "y",
        })
    );

    let earlier_conformer_and_dimension = CoordinateBlock {
        conformers_2d: vec![
            Conformer2D::new(42, vec![[0.0, f64::NEG_INFINITY]]),
            Conformer2D::new(43, vec![[f64::NAN, 0.0]]),
        ],
        conformers_3d: vec![Conformer3D::new(44, vec![[f64::NAN, 0.0, 0.0]], true)],
        source_coordinate_dim: Some(CoordinateDimension::ThreeD),
    };
    assert_eq!(
        earlier_conformer_and_dimension.validate_for_atom_count(1),
        Err(CoordinateValidationError::NonFiniteCoordinate {
            dimension: "2D",
            conformer: 42,
            atom: 0,
            axis: "y",
        })
    );

    let duplicate_precedes_second_validation = CoordinateBlock {
        conformers_2d: vec![
            Conformer2D::new(45, vec![[0.0, 0.0]]),
            Conformer2D::new(45, Vec::new()),
        ],
        ..Default::default()
    };
    assert_eq!(
        duplicate_precedes_second_validation.validate_for_atom_count(1),
        Err(CoordinateValidationError::DuplicateConformerId {
            dimension: "2D",
            id: 45,
        })
    );

    let valid_mixed = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(1, vec![[-0.0, f64::MAX]])],
        conformers_3d: vec![Conformer3D::new(2, vec![[f64::MIN, 0.0, -0.0]], false)],
        source_coordinate_dim: Some(CoordinateDimension::TwoD),
    };
    assert_eq!(valid_mixed.validate_for_atom_count(1), Ok(()));
}

#[test]
fn remap_preserves_order_bits_properties_flags_and_source_dimension() {
    let two_d_first = [[-0.0, 10.0], [11.0, 12.0], [13.0, f64::MAX]];
    let two_d_second = [[20.0, -0.0], [21.0, 22.0], [f64::MIN, 24.0]];
    let three_d_first = [
        [-0.0, 30.0, 31.0],
        [32.0, 33.0, 34.0],
        [35.0, 36.0, f64::MAX],
    ];
    let three_d_second = [
        [40.0, 41.0, -0.0],
        [42.0, 43.0, 44.0],
        [f64::MIN, 46.0, 47.0],
    ];
    let mut block = CoordinateBlock {
        conformers_2d: vec![
            Conformer2D::new(70, two_d_first.to_vec()).with_prop("which", "2d-first"),
            Conformer2D::new(90, two_d_second.to_vec()).with_prop("which", "2d-second"),
        ],
        conformers_3d: vec![
            Conformer3D::new(80, three_d_first.to_vec(), true).with_prop("which", "3d-first"),
            Conformer3D::new(99, three_d_second.to_vec(), false).with_prop("which", "3d-second"),
        ],
        source_coordinate_dim: Some(CoordinateDimension::ThreeD),
    };

    block.remap_topology(&[2, 0]);

    assert_eq!(
        block.source_coordinate_dim,
        Some(CoordinateDimension::ThreeD)
    );
    assert_eq!(
        block
            .conformers_2d
            .iter()
            .map(Conformer2D::id)
            .collect::<Vec<_>>(),
        vec![0, 1]
    );
    assert_eq!(
        block
            .conformers_3d
            .iter()
            .map(Conformer3D::id)
            .collect::<Vec<_>>(),
        vec![0, 1]
    );
    assert_eq!(
        bits_2d(block.conformers_2d[0].coordinates()),
        bits_2d(&[two_d_first[2], two_d_first[0]])
    );
    assert_eq!(
        bits_2d(block.conformers_2d[1].coordinates()),
        bits_2d(&[two_d_second[2], two_d_second[0]])
    );
    assert_eq!(
        bits_3d(block.conformers_3d[0].coordinates()),
        bits_3d(&[three_d_first[2], three_d_first[0]])
    );
    assert_eq!(
        bits_3d(block.conformers_3d[1].coordinates()),
        bits_3d(&[three_d_second[2], three_d_second[0]])
    );
    assert_eq!(
        block.conformers_2d[0]
            .props()
            .get("which")
            .map(String::as_str),
        Some("2d-first")
    );
    assert_eq!(
        block.conformers_2d[1]
            .props()
            .get("which")
            .map(String::as_str),
        Some("2d-second")
    );
    assert_eq!(
        block.conformers_3d[0]
            .props()
            .get("which")
            .map(String::as_str),
        Some("3d-first")
    );
    assert_eq!(
        block.conformers_3d[1]
            .props()
            .get("which")
            .map(String::as_str),
        Some("3d-second")
    );
    assert!(block.conformers_3d[0].is_3d());
    assert!(!block.conformers_3d[1].is_3d());
    assert_eq!(block.validate_for_atom_count(2), Ok(()));
}
