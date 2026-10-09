use cosmolkit::{BatchErrorMode, BatchParams, MoleculeBatch, SmilesParseParams};
use std::error::Error;
fn duplicate_leaf(error: &(dyn Error + 'static)) -> bool {
    let mut current = Some(error);
    while let Some(error) = current {
        if let Some(error) = error.downcast_ref::<cosmolkit_model::StereoGroupError>() {
            return *error == cosmolkit_model::StereoGroupError::DuplicateAtom;
        }
        current = error.source();
    }
    false
}
#[test]
fn recovery_chem29_natural_batch_strict_preserves_duplicate_leaf_and_row_index() {
    let input = vec!["C".to_owned(), "CC |o1:0,o1:0|".to_owned(), "O".to_owned()];
    let error = MoleculeBatch::from_smiles_list(&input).unwrap_err();
    assert_eq!(error.errors, 1);
    assert_eq!(error.record_errors.len(), 1);
    assert_eq!(error.record_errors[0].index, 1);
    assert!(duplicate_leaf(&error.record_errors[0]));
}
#[test]
fn recovery_chem29_natural_batch_keep_errors_retains_order_and_typed_cause() {
    let input = vec!["C".to_owned(), "CC |o1:0,o1:0|".to_owned(), "O".to_owned()];
    let batch = MoleculeBatch::from_smiles_list_with_params(
        &input,
        &SmilesParseParams::default(),
        &BatchParams {
            errors: Some(BatchErrorMode::KeepErrors),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(batch.valid_mask(), [true, false, true]);
    assert_eq!(batch.errors()[0].index, 1);
    assert!(duplicate_leaf(&batch.errors()[0]));
}
