use cosmolkit::{
    Coordinate3DInputParams, CoordinateInputError, EmbedParams, Molecule, OperationError,
    Replace3DCoordinatesParams,
};
use cosmolkit_model::{Conformer2D, Conformer3D, CoordinateBlock};

fn input(coords: Vec<[f64; 3]>) -> Vec<Vec<f64>> {
    coords.into_iter().map(|row| row.to_vec()).collect()
}
fn replace(id: usize) -> Replace3DCoordinatesParams {
    Replace3DCoordinatesParams { conformer_id: id }
}
fn flags(is_3d: bool) -> Coordinate3DInputParams {
    Coordinate3DInputParams { is_3d }
}

fn fixture() -> Molecule {
    let water = Molecule::from_smiles("O")
        .unwrap()
        .with_hydrogens()
        .unwrap();
    Molecule::from_parts(
        water.topology().clone(),
        CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(9, vec![[2.0, 3.0]; 3])],
            conformers_3d: vec![
                Conformer3D::new(4, vec![[1.0, 2.0, 3.0]; 3], true),
                Conformer3D::new(9, vec![[4.0, 5.0, 6.0]; 3], false),
            ],
            ..Default::default()
        },
        water.properties().clone(),
    )
    .unwrap()
}
fn preserved_other_blocks(actual: &Molecule, original: &Molecule) {
    assert!(std::ptr::eq(actual.topology(), original.topology()));
    assert!(std::ptr::eq(actual.properties(), original.properties()));
    assert_eq!(
        actual.to_builder().coordinates().conformers_2d,
        original.to_builder().coordinates().conformers_2d
    );
}
fn seeded() -> EmbedParams {
    let mut params = EmbedParams::etkdg_v3();
    params.random_seed = 42;
    params.max_iterations = 3;
    params.num_threads = 1;
    params.track_failures = true;
    params
}

#[test]
fn manual_replacement_selects_actual_id_and_preserves_other_rows_and_2d() {
    let source = fixture();
    let before = source.clone();
    let coords = vec![[7.0, 8.0, 9.0]; 3];
    let value = source
        .with_3d_coordinates_with_params(input(coords.clone()), &replace(9))
        .unwrap();
    assert_eq!(source, before);
    assert_eq!(value.conformers_3d()[0], source.conformers_3d()[0]);
    assert_eq!(value.conformers_3d()[1].coordinates(), coords);
    assert!(!value.conformers_3d()[1].is_3d());
    preserved_other_blocks(&value, &source);
    let mut inplace = source.clone();
    inplace
        .set_3d_coordinates_with_params_(input(coords), &replace(9))
        .unwrap();
    assert_eq!(inplace, value);
    let error = inplace
        .set_3d_coordinates_with_params_(input(vec![[0.0; 3]; 3]), &replace(0))
        .unwrap_err();
    assert_eq!(
        error,
        OperationError::CoordinateInput(CoordinateInputError::ConformerNotFound {
            conformer_id: 0,
            count: 2
        })
    );
    assert_eq!(inplace, value);
}

#[test]
fn manual_append_only_and_clear_preserve_independent_dimension() {
    let source = fixture();
    let coords = vec![[1.0, 0.0, 0.0]; 3];
    let appended = source
        .with_added_3d_conformer_with_params(input(coords.clone()), &flags(true))
        .unwrap();
    assert_eq!(appended.num_3d_conformers(), 3);
    assert_eq!(appended.conformers_3d()[2].id(), 10);
    preserved_other_blocks(&appended, &source);
    let mut inplace = source.clone();
    assert_eq!(
        inplace
            .add_3d_conformer_with_params_(input(coords.clone()), &flags(true))
            .unwrap(),
        2
    );
    assert_eq!(inplace, appended);
    let only = source
        .with_only_3d_conformer_with_params(input(coords.clone()), &flags(false))
        .unwrap();
    assert_eq!(only.num_3d_conformers(), 1);
    assert_eq!(only.conformers_3d()[0].id(), 0);
    assert!(!only.conformers_3d()[0].is_3d());
    preserved_other_blocks(&only, &source);
    assert_eq!(
        inplace
            .set_only_3d_conformer_with_params_(input(coords), &flags(false))
            .unwrap(),
        0
    );
    assert_eq!(inplace, only);
    let cleared = source.with_cleared_3d_conformers().unwrap();
    assert_eq!(cleared.num_3d_conformers(), 0);
    preserved_other_blocks(&cleared, &source);
    inplace.clear_3d_conformers_().unwrap();
    assert_eq!(inplace, cleared);
    assert_eq!(source.num_3d_conformers(), 2);
}

#[test]
fn manual_invalid_rows_nonfinite_and_id_overflow_are_atomic() {
    let source = fixture();
    let mut inplace = source.clone();
    assert!(matches!(
        inplace.add_3d_conformer_with_params_(input(vec![[0.0; 3]; 2]), &flags(true)),
        Err(OperationError::CoordinateInput(
            CoordinateInputError::RowCount {
                expected: 3,
                actual: 2,
                ..
            }
        ))
    ));
    assert_eq!(inplace, source);
    assert!(matches!(
        inplace
            .set_only_3d_conformer_with_params_(input(vec![[f64::NAN, 0.0, 0.0]; 3]), &flags(true)),
        Err(OperationError::CoordinateInput(
            CoordinateInputError::NonFinite {
                row: 0,
                column: 0,
                ..
            }
        ))
    ));
    assert_eq!(inplace, source);
    assert!(matches!(
        inplace.set_3d_coordinates_with_params_(
            input(vec![[0.0, f64::INFINITY, 0.0]; 3]),
            &replace(4)
        ),
        Err(OperationError::CoordinateInput(
            CoordinateInputError::NonFinite {
                row: 0,
                column: 1,
                ..
            }
        ))
    ));
    assert_eq!(inplace, source);
    preserved_other_blocks(&inplace, &source);
    let overflow = Molecule::from_parts(
        source.topology().clone(),
        CoordinateBlock {
            conformers_3d: vec![Conformer3D::new(usize::MAX, vec![[0.0; 3]; 3], true)],
            ..Default::default()
        },
        source.properties().clone(),
    )
    .unwrap();
    let mut retained = overflow.clone();
    assert!(matches!(
        retained.add_3d_conformer_with_params_(input(vec![[0.0; 3]; 3]), &flags(true)),
        Err(OperationError::CoordinateInput(
            CoordinateInputError::ConformerIdOverflow { max_id: usize::MAX }
        ))
    ));
    assert_eq!(retained, overflow);
}

#[test]
fn canonical_reports_append_ids_copy_params_and_commit_inplace_after_validation() {
    let source = fixture();
    let mut params = seeded();
    params.clear_conformers = false;
    let value = source
        .with_3d_conformer_result_with_params(&params)
        .unwrap();
    assert!(value.ok());
    assert_eq!(value.conf_id(), 10);
    assert_eq!(value.molecule().num_3d_conformers(), 3);
    assert_eq!(value.params().random_seed, 42);
    // RDKit .6 EmbedFailureCauses::END_OF_ENUM is 15.
    assert_eq!(value.params().failures.len(), 15);
    assert!(params.failures.is_empty());
    preserved_other_blocks(value.molecule(), &source);
    let mut inplace = source.clone();
    let report = inplace
        .embed_3d_conformer_result_with_params_(&params)
        .unwrap();
    assert_eq!(inplace, value.molecule);
    assert_eq!(report.molecule(), &inplace);
    assert!(std::ptr::eq(
        report.molecule().topology(),
        inplace.topology()
    ));
    let multiple = source
        .with_3d_conformers_result_with_params(3, &params)
        .unwrap();
    assert_eq!(multiple.conf_ids(), &[10, 11, 12]);
    assert_eq!(multiple.requested_num_confs(), 3);
    assert_eq!(multiple.generated_count(), 3);
    assert_eq!(multiple.molecule().num_3d_conformers(), 5);
    let mut inplace = source.clone();
    let report = inplace
        .embed_3d_conformers_result_with_params_(3, &params)
        .unwrap();
    assert_eq!(report.molecule(), &inplace);
    assert_eq!(inplace, multiple.molecule);
    let before = inplace.clone();
    params.et_version = 99;
    let error = inplace
        .embed_3d_conformer_result_with_params_(&params)
        .unwrap_err();
    assert!(matches!(error, OperationError::Conformer(_)));
    assert!(std::error::Error::source(&error).is_some());
    assert_eq!(inplace, before);
}

#[test]
fn all_generation_value_inplace_forms_and_bounds_query_are_real() {
    let source = fixture();
    let params = seeded();
    let value = source.with_3d_conformer_with_params(&params).unwrap();
    let mut inplace = source.clone();
    inplace.embed_3d_conformer_with_params_(&params).unwrap();
    assert_eq!(inplace, value);
    let value = source.with_3d_conformers_with_params(2, &params).unwrap();
    let mut inplace = source.clone();
    inplace
        .embed_3d_conformers_with_params_(2, &params)
        .unwrap();
    assert_eq!(inplace, value);
    assert_eq!(source.with_3d_conformer().unwrap().num_3d_conformers(), 1);
    let mut inplace = source.clone();
    inplace.embed_3d_conformer_().unwrap();
    assert_eq!(inplace.num_3d_conformers(), 1);
    assert_eq!(source.with_3d_conformers(2).unwrap().num_3d_conformers(), 2);
    inplace.embed_3d_conformers_(2).unwrap();
    assert_eq!(inplace.num_3d_conformers(), 2);
    assert!(source.with_3d_conformer_result().unwrap().ok());
    assert!(inplace.embed_3d_conformer_result_().unwrap().ok());
    assert_eq!(
        source
            .with_3d_conformers_result(2)
            .unwrap()
            .generated_count(),
        2
    );
    assert_eq!(
        inplace
            .embed_3d_conformers_result_(2)
            .unwrap()
            .generated_count(),
        2
    );
    let before = source.clone();
    let bounds = source.dg_bounds_matrix().unwrap();
    assert_eq!(bounds.len(), 3);
    assert!(bounds.iter().all(|r| r.len() == 3));
    assert_eq!(source, before);
    assert_eq!(source.num_3d_conformers(), 2);
}
