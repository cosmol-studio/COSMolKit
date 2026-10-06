//! Original d892ec3 batch protocol cases, canonical signature adaptation only.
//! Delivery proposal; original expectations and inputs retained for independent review.
use cosmolkit::{
    BatchErrorMode, BatchParams, BatchRecord, Coordinate2DParams, MoleculeBatch, SmilesParseParams,
};
fn kept(smiles: &[String]) -> MoleculeBatch {
    MoleculeBatch::from_smiles_list_with_params(
        smiles,
        &SmilesParseParams::default(),
        &BatchParams {
            errors: BatchErrorMode::KeepErrors,
            ..Default::default()
        },
    )
    .unwrap()
}

#[test]
fn from_smiles_list_preserves_order_and_keeps_record_errors() {
    let smiles = vec!["CCO".to_string(), "C1".to_string(), "N".to_string()];
    let batch = kept(&smiles);
    assert_eq!(batch.len(), 3);
    assert_eq!(batch.valid_mask(), vec![true, false, true]);
    assert_eq!(batch.valid_count(), 2);
    assert_eq!(batch.invalid_count(), 1);
    let errors = batch.errors();
    assert_eq!(errors.len(), 1);
    assert_eq!(errors[0].index, 1);
    assert_eq!(errors[0].operation, "batch.from_smiles_list");
}

#[test]
fn filter_valid_preserves_surviving_input_order() {
    let smiles = vec!["C".to_string(), "C1".to_string(), "O".to_string()];
    let filtered = kept(&smiles).with_valid_records();
    assert_eq!(filtered.len(), 2);
    assert_eq!(filtered.valid_mask(), vec![true, true]);
}

#[test]
fn batch_configuration_defaults_are_none() {
    let batch = kept(&["CCO".to_string()]);
    assert_eq!(batch.parallel_jobs(), None);
    assert_eq!(batch.progress_bar(), None);
}

#[test]
fn batch_configuration_can_be_set_and_cleared() {
    let batch = kept(&["CCO".to_string()])
        .with_parallel_jobs(Some(4))
        .unwrap()
        .with_progress_bar(Some(true));
    assert_eq!(batch.parallel_jobs(), Some(4));
    assert_eq!(batch.progress_bar(), Some(true));

    let cleared = batch
        .with_parallel_jobs(None)
        .unwrap()
        .with_progress_bar(None);
    assert_eq!(cleared.parallel_jobs(), None);
    assert_eq!(cleared.progress_bar(), None);
}

#[test]
fn batch_configuration_is_preserved_across_transforms() {
    let batch = kept(&["CC".to_string()])
        .with_parallel_jobs(Some(2))
        .unwrap()
        .with_progress_bar(Some(false));

    let transformed = batch
        .with_2d_coordinates()
        .expect("2d coordinates should succeed");

    assert_eq!(transformed.parallel_jobs(), Some(2));
    assert_eq!(transformed.progress_bar(), Some(false));
}

#[test]
fn batch_configuration_is_preserved_across_filter_valid() {
    let batch = kept(&["CCO".to_string(), "C1".to_string()])
        .with_parallel_jobs(Some(8))
        .unwrap()
        .with_progress_bar(Some(true));

    let filtered = batch.with_valid_records();

    assert_eq!(filtered.parallel_jobs(), Some(8));
    assert_eq!(filtered.progress_bar(), Some(true));
}

#[test]
fn transform_options_preserve_batch_configuration() {
    let batch = kept(&["CC".to_string()])
        .with_parallel_jobs(Some(4))
        .unwrap()
        .with_progress_bar(Some(true));

    let transformed = batch
        .with_2d_coordinates_with_params(
            &Coordinate2DParams::default(),
            &BatchParams {
                errors: BatchErrorMode::Strict,
                n_jobs: Some(1),
                progress_bar: Some(false),
            },
        )
        .expect("2d coordinates should succeed");

    assert_eq!(transformed.parallel_jobs(), Some(4));
    assert_eq!(transformed.progress_bar(), Some(true));
}

#[test]
fn strict_aggregates_all_original_parse_errors_and_transform_retains_original_indices() {
    let smiles = vec!["C1".into(), "CCO".into(), "bad".into()];
    let strict = MoleculeBatch::from_smiles_list(&smiles).unwrap_err();
    assert_eq!(strict.errors, 2);
    assert_eq!(
        strict
            .record_errors
            .iter()
            .map(|e| e.index)
            .collect::<Vec<_>>(),
        vec![0, 2]
    );
    assert!(
        strict
            .record_errors
            .iter()
            .all(|e| std::error::Error::source(e).is_some())
    );
    let batch = kept(&smiles)
        .with_parallel_jobs(Some(2))
        .unwrap()
        .with_progress_bar(Some(false));
    assert_eq!(batch.invalid_mask(), vec![true, false, true]);
    let keep = BatchParams {
        errors: BatchErrorMode::KeepErrors,
        n_jobs: Some(1),
        progress_bar: Some(false),
    };
    let transformed = batch
        .with_hydrogens_with_params(&Default::default(), &keep)
        .unwrap();
    assert_eq!(transformed.valid_mask(), batch.valid_mask());
    assert_eq!(
        transformed
            .errors()
            .iter()
            .map(|e| (e.index, e.operation))
            .collect::<Vec<_>>(),
        vec![(0, "batch.from_smiles_list"), (2, "batch.from_smiles_list")]
    );
    assert_eq!(transformed.parallel_jobs(), Some(2));
    assert_eq!(transformed.progress_bar(), Some(false));
    assert!(
        transformed.to_list()[1].as_ref().unwrap().atoms().len()
            > batch.to_list()[1].as_ref().unwrap().atoms().len()
    );
    assert_eq!(batch.to_list()[1].as_ref().unwrap().atoms().len(), 3);
    assert_eq!(
        transformed
            .without_hydrogens_with_params(&Default::default(), &keep)
            .unwrap()
            .to_list()[1]
            .as_ref()
            .unwrap()
            .atoms()
            .len(),
        3
    );
    let strict = batch.with_hydrogens().unwrap_err();
    assert_eq!(
        strict
            .record_errors
            .iter()
            .map(|e| e.index)
            .collect::<Vec<_>>(),
        vec![0, 2]
    );
    assert_eq!(batch.iter().count(), 3);
    assert!(matches!(batch.get(0), Some(BatchRecord::Error(_))));
}
#[test]
fn every_authorized_transform_matches_its_scalar_owner_and_leaves_input_unchanged() {
    let batch = MoleculeBatch::from_smiles_list(&["CCO".into(), "c1ccccc1".into()]).unwrap();
    let input = batch
        .to_list()
        .into_iter()
        .map(Option::unwrap)
        .collect::<Vec<_>>();
    let scalar = input
        .iter()
        .map(|m| m.to_smiles().unwrap())
        .collect::<Vec<_>>();
    for (result, expected) in [
        (
            batch.sanitize().unwrap(),
            input
                .iter()
                .map(|m| m.sanitize().unwrap())
                .collect::<Vec<_>>(),
        ),
        (
            batch.with_hydrogens().unwrap(),
            input.iter().map(|m| m.with_hydrogens().unwrap()).collect(),
        ),
        (
            batch.without_hydrogens().unwrap(),
            input
                .iter()
                .map(|m| m.without_hydrogens().unwrap())
                .collect(),
        ),
        (
            batch.with_kekulized_bonds().unwrap(),
            input
                .iter()
                .map(|m| m.with_kekulized_bonds().unwrap())
                .collect(),
        ),
        (
            batch.with_2d_coordinates().unwrap(),
            input
                .iter()
                .map(|m| m.with_2d_coordinates().unwrap())
                .collect(),
        ),
    ] {
        for (actual, expected) in result.to_list().into_iter().zip(expected) {
            let actual = actual.unwrap();
            assert_eq!(actual.topology(), expected.topology());
            assert_eq!(actual.coordinates_2d(), expected.coordinates_2d());
        }
    }
    assert_eq!(
        batch
            .to_list()
            .iter()
            .map(|m| m.as_ref().unwrap().to_smiles().unwrap())
            .collect::<Vec<_>>(),
        scalar
    );
}
