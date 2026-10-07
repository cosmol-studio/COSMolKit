use cosmolkit_core::*;
use cosmolkit_model::{Conformer2D, Conformer3D, CoordinateBlock, CoordinateDimension};

fn xyz() -> Vec<Vec<f64>> {
    vec![vec![0., 0., 0.], vec![1., 2., 3.]]
}
fn mixed() -> CoordinateBlock {
    CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(2, vec![[0., 0.], [1., 0.]])],
        conformers_3d: vec![
            Conformer3D::new(3, vec![[0.; 3]; 2], false).with_prop("kind", "original"),
            Conformer3D::new(7, vec![[1.; 3]; 2], true),
        ],
        source_conformer_order: None,
        source_coordinate_dim: Some(CoordinateDimension::ThreeD),
    }
}

#[test]
fn policy_defaults_casefold_and_original_zero_threshold() {
    assert_eq!(
        Coordinate2DInputParams::default().z_policy,
        CoordinateZPolicy::Ignore
    );
    assert!(Coordinate3DInputParams::default().is_3d);
    assert_eq!(Replace3DCoordinatesParams::default().conformer_id, 0);
    assert_eq!(
        CoordinateZPolicy::from_name("ReQuIrE_ZeRo").unwrap(),
        CoordinateZPolicy::RequireZero
    );
    assert!(matches!(
        CoordinateZPolicy::from_name(" require_zero"),
        Err(CoordinateInputError::UnknownZPolicy { .. })
    ));
    let params = Coordinate2DInputParams {
        z_policy: CoordinateZPolicy::RequireZero,
    };
    assert_eq!(
        coordinates_2d_from_input(2, vec![vec![1., 2., 1e-12], vec![3., 4., -1e-12]], &params)
            .unwrap(),
        vec![[1., 2.], [3., 4.]]
    );
    assert!(matches!(
        coordinates_2d_from_input(1, vec![vec![1., 2., 1.0001e-12]], &params),
        Err(CoordinateInputError::NonZeroZ { row: 0, .. })
    ));
    let error = Coordinate2DInputParams {
        z_policy: CoordinateZPolicy::Error,
    };
    assert!(matches!(
        coordinates_2d_from_input(1, vec![vec![1., 2., 0.]], &error),
        Err(CoordinateInputError::ThreeColumnsForbidden)
    ));
    assert!(
        coordinates_2d_from_input(0, vec![], &error)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn matrix_validation_keeps_row_shape_finite_order_including_ignored_z() {
    assert!(matches!(
        coordinates_3d_from_input(2, vec![vec![f64::NAN]]),
        Err(CoordinateInputError::RowCount {
            expected: 2,
            actual: 1,
            ..
        })
    ));
    assert!(matches!(
        coordinates_3d_from_input(2, vec![vec![f64::NAN; 3], vec![0.; 2]]),
        Err(CoordinateInputError::Shape {
            row: 1,
            columns: 2,
            ..
        })
    ));
    assert!(matches!(
        coordinates_2d_from_input(
            1,
            vec![vec![1., 2., f64::INFINITY]],
            &Coordinate2DInputParams::default()
        ),
        Err(CoordinateInputError::NonFinite {
            row: 0,
            column: 2,
            ..
        })
    ));
}

#[test]
fn two_dimensional_install_preserves_xyz_metadata_and_explicit_provenance() {
    let mut block = mixed();
    let old = block.conformers_3d.clone();
    install_2d_coordinates(&mut block, 2, xyz(), &Coordinate2DInputParams::default()).unwrap();
    assert_eq!(block.conformers_2d.len(), 1);
    assert_eq!(block.conformers_2d[0].id(), 0);
    assert_eq!(block.conformers_2d[0].coordinates(), &[[0., 0.], [1., 2.]]);
    assert_eq!(block.conformers_3d, old);
    assert_eq!(block.source_coordinate_dim, Some(CoordinateDimension::TwoD));
}

#[test]
fn actual_id_replacement_keeps_is3d_and_resets_only_selected_properties() {
    let mut block = mixed();
    let other = block.conformers_3d[1].clone();
    let xy = block.conformers_2d.clone();
    replace_3d_coordinates(
        &mut block,
        2,
        xyz(),
        &Replace3DCoordinatesParams { conformer_id: 3 },
    )
    .unwrap();
    assert_eq!(block.conformers_3d[0].id(), 3);
    assert!(!block.conformers_3d[0].is_3d());
    assert_eq!(
        block.conformers_3d[0],
        Conformer3D::new(3, vec![[0., 0., 0.], [1., 2., 3.]], false)
    );
    assert_eq!(block.conformers_3d[1], other);
    assert_eq!(block.conformers_2d, xy);
    assert_eq!(
        block.source_coordinate_dim,
        Some(CoordinateDimension::ThreeD)
    );
}

#[test]
fn append_report_is_vector_position_ids_use_max_and_only_clear_keep_xy() {
    let mut block = mixed();
    let xy = block.conformers_2d.clone();
    assert_eq!(
        append_3d_conformer(
            &mut block,
            2,
            xyz(),
            &Coordinate3DInputParams { is_3d: false }
        )
        .unwrap(),
        2
    );
    assert_eq!(block.conformers_3d[2].id(), 8);
    assert!(!block.conformers_3d[2].is_3d());
    assert_eq!(
        install_only_3d_conformer(
            &mut block,
            2,
            xyz(),
            &Coordinate3DInputParams { is_3d: false }
        )
        .unwrap(),
        0
    );
    assert_eq!(block.conformers_3d.len(), 1);
    assert_eq!(block.conformers_3d[0].id(), 0);
    assert_eq!(block.source_coordinate_dim, Some(CoordinateDimension::TwoD));
    clear_3d_conformers(&mut block);
    assert_eq!(block.conformers_2d, xy);
    assert_eq!(block.source_coordinate_dim, Some(CoordinateDimension::TwoD));
    block.conformers_2d.clear();
    clear_3d_conformers(&mut block);
    assert_eq!(block.source_coordinate_dim, None);
}

#[test]
fn every_fallible_transition_is_atomic_including_identifier_overflow() {
    let mut block = mixed();
    let before = block.clone();
    assert!(
        install_2d_coordinates(
            &mut block,
            2,
            vec![vec![0.; 2]],
            &Coordinate2DInputParams::default()
        )
        .is_err()
    );
    assert_eq!(block, before);
    assert!(matches!(
        replace_3d_coordinates(
            &mut block,
            2,
            xyz(),
            &Replace3DCoordinatesParams { conformer_id: 99 }
        ),
        Err(CoordinateInputError::ConformerNotFound {
            conformer_id: 99,
            count: 2
        })
    ));
    assert_eq!(block, before);
    assert!(
        append_3d_conformer(
            &mut block,
            2,
            vec![vec![f64::NAN; 3]; 2],
            &Coordinate3DInputParams::default()
        )
        .is_err()
    );
    assert_eq!(block, before);
    assert!(
        install_only_3d_conformer(
            &mut block,
            2,
            vec![vec![0.; 2]; 2],
            &Coordinate3DInputParams::default()
        )
        .is_err()
    );
    assert_eq!(block, before);
    block
        .conformers_3d
        .push(Conformer3D::new(usize::MAX, vec![[0.; 3]; 2], true));
    let before = block.clone();
    assert!(matches!(
        append_3d_conformer(&mut block, 2, xyz(), &Coordinate3DInputParams::default()),
        Err(CoordinateInputError::ConformerIdOverflow { max_id: usize::MAX })
    ));
    assert_eq!(block, before);
}

#[test]
fn readonly_xyz_selects_dimension_scoped_id_without_position_fallback() {
    let block = mixed();
    let before = block.clone();
    assert!(std::ptr::eq(
        coordinates_3d_for_id(&block, 7).unwrap(),
        block.conformers_3d[1].coordinates()
    ));
    assert_eq!(
        coordinates_3d_for_id(&block, 3).unwrap(),
        block.conformers_3d[0].coordinates()
    );
    assert!(matches!(
        coordinates_3d_for_id(&block, 1),
        Err(Coordinate3DReadError::ConformerNotFound {
            conformer_id: 1,
            count: 2
        })
    ));
    assert_eq!(block, before);
    assert!(matches!(
        coordinates_3d_for_id(&CoordinateBlock::default(), 0),
        Err(Coordinate3DReadError::ConformerNotFound {
            conformer_id: 0,
            count: 0
        })
    ));
}

#[test]
fn replacement_id_lookup_rejects_positions_and_missing_default_before_any_mutation() {
    let mut block = mixed();
    block
        .conformers_3d
        .push(Conformer3D::new(11, vec![[2.; 3]; 2], true));
    let before = block.clone();
    for conformer_id in [1, Replace3DCoordinatesParams::default().conformer_id] {
        assert_eq!(
            replace_3d_coordinates(
                &mut block,
                2,
                xyz(),
                &Replace3DCoordinatesParams { conformer_id }
            ),
            Err(CoordinateInputError::ConformerNotFound {
                conformer_id,
                count: 3
            })
        );
        assert_eq!(block, before);
    }
    replace_3d_coordinates(
        &mut block,
        2,
        xyz(),
        &Replace3DCoordinatesParams { conformer_id: 7 },
    )
    .unwrap();
    assert_eq!(
        block.conformers_3d[1],
        Conformer3D::new(7, vec![[0., 0., 0.], [1., 2., 3.]], true)
    );
    assert_eq!(block.conformers_3d[0], before.conformers_3d[0]);
    assert_eq!(block.conformers_3d[2], before.conformers_3d[2]);
    assert_eq!(block.conformers_2d, before.conformers_2d);
    assert_eq!(block.source_coordinate_dim, before.source_coordinate_dim);
}

mod source_order_regression {
    use cosmolkit_core::{Coordinate3DInputParams, append_3d_conformer};
    use cosmolkit_model::{
        Conformer2D, CoordinateBlock, CoordinateDimension, CoordinateSourceConformer,
    };
    #[test]
    fn append_records_actual_source_occurrence() {
        let mut block = CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(9, vec![[1.0, -0.0]])],
            source_conformer_order: Some(vec![CoordinateDimension::TwoD]),
            ..Default::default()
        };
        let position = append_3d_conformer(
            &mut block,
            1,
            vec![vec![2.0, 3.0, 4.0]],
            &Coordinate3DInputParams { is_3d: false },
        )
        .unwrap();
        assert_eq!(position, 0);
        assert_eq!(
            block.conformers_2d[0].coordinates()[0][1].to_bits(),
            (-0.0_f64).to_bits()
        );
        assert_eq!(
            block.source_conformer_order,
            Some(vec![CoordinateDimension::TwoD, CoordinateDimension::ThreeD])
        );
        block.validate_for_atom_count(1).unwrap();
        assert!(matches!(
            block.first_source_conformer().unwrap(),
            Some(CoordinateSourceConformer::TwoD(_))
        ));
    }
}
