use cosmolkit_model::{Conformer2D, Conformer3D, CoordinateValidationError};

#[test]
fn checked_xy_input_rejects_first_nonfinite_without_changing_raw_storage() {
    let nan = f64::from_bits(0xfff8_0000_0000_007b);
    let conformer = Conformer2D::new(7, vec![[1.0, 2.0], [nan, f64::INFINITY]]);
    assert_eq!(conformer.validate_for_atom_count(2), Ok(()));
    assert_eq!(
        conformer.validate_checked_for_atom_count(2),
        Err(CoordinateValidationError::NonFiniteCoordinate {
            dimension: "2D",
            conformer: 7,
            atom: 1,
            axis: "x",
        })
    );
    assert_eq!(conformer.coordinates()[1][0].to_bits(), nan.to_bits());
}

#[test]
fn checked_xyz_input_keeps_row_axis_error_order_and_signed_zero() {
    let conformer = Conformer3D::new(
        19,
        vec![[-0.0, 1.0, f64::NEG_INFINITY], [f64::NAN, 2.0, 3.0]],
        false,
    );
    assert_eq!(conformer.validate_for_atom_count(2), Ok(()));
    assert_eq!(
        conformer.validate_checked_for_atom_count(2),
        Err(CoordinateValidationError::NonFiniteCoordinate {
            dimension: "3D",
            conformer: 19,
            atom: 0,
            axis: "z",
        })
    );
    assert_eq!(
        conformer.coordinates()[0][0].to_bits(),
        (-0.0_f64).to_bits()
    );
    assert!(!conformer.is_3d());
    assert_eq!(
        Conformer3D::new(0, vec![[-0.0, f64::MAX, f64::MIN_POSITIVE]], true)
            .validate_checked_for_atom_count(1),
        Ok(())
    );
}

#[test]
fn checked_input_reports_row_count_before_nonfinite_for_both_dimensions() {
    assert_eq!(
        Conformer2D::new(3, vec![[f64::NAN, 0.0]]).validate_checked_for_atom_count(2),
        Err(CoordinateValidationError::RowCount {
            dimension: "2D",
            conformer: 3,
            rows: 1,
            atom_count: 2,
        })
    );
    assert_eq!(
        Conformer3D::new(4, vec![[f64::NAN, 0.0, 0.0]], true).validate_checked_for_atom_count(2),
        Err(CoordinateValidationError::RowCount {
            dimension: "3D",
            conformer: 4,
            rows: 1,
            atom_count: 2,
        })
    );
}
