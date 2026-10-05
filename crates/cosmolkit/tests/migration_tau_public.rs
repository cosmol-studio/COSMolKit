#![cfg(all(feature = "cap-tautomer", feature = "cap-smiles"))]
use cosmolkit::*;
use std::sync::{Arc, Mutex};

#[test]
fn canonical_methods_return_ordered_source_values_and_keep_source_immutable() {
    let source = Molecule::from_smiles("CC(C)=O").unwrap();
    let before = source.clone();
    let result = source.enumerate_tautomers().unwrap();
    assert_eq!(result.canonical_smiles(), ["C=C(C)O", "CC(C)=O"]);
    assert_eq!(result.status(), TautomerEnumerationStatus::Completed);
    assert_eq!(
        result
            .iter()
            .map(|m| m.tautomer_score().unwrap().total())
            .collect::<Vec<_>>(),
        [1, 5]
    );
    assert_eq!(
        source.canonical_tautomer().unwrap().to_smiles().unwrap(),
        "CC(C)=O"
    );
    assert_eq!(source, before);
    assert!(result.get(2).is_none());
    assert_eq!(result[0], *result.get(0).unwrap());
    assert_eq!(result.clone(), result);
}
struct Cancel {
    calls: Mutex<Vec<(usize, Vec<String>)>>,
}
impl TautomerEnumerationCallback for Cancel {
    fn should_continue(
        &self,
        source: TautomerMoleculeView<'_>,
        result: TautomerProgress<'_>,
    ) -> Result<bool, TautomerRunError> {
        assert_eq!(source.num_atoms(), 4);
        assert_eq!(source.to_smiles()?, "CC(C)=O");
        self.calls.lock().unwrap().push((
            result.len(),
            result.entries().map(|(key, _)| key.to_owned()).collect(),
        ));
        Ok(false)
    }
}
#[test]
fn callback_inspects_borrowed_progress_and_cancels_before_apply() {
    let source = Molecule::from_smiles("CC(C)=O").unwrap();
    let before = source.clone();
    let callback = Arc::new(Cancel {
        calls: Mutex::new(Vec::new()),
    });
    let mut params = TautomerParams::default();
    params.set_callback(Some(callback.clone()));
    let result = source.enumerate_tautomers_with_params(&params).unwrap();
    assert_eq!(result.status(), TautomerEnumerationStatus::Canceled);
    assert_eq!(result.canonical_smiles(), ["CC(C)=O"]);
    assert!(result.modified_atoms().is_empty());
    assert!(result.modified_bonds().is_empty());
    assert_eq!(
        callback.calls.lock().unwrap().as_slice(),
        [(1, vec!["CC(C)=O".to_owned()])]
    );
    assert_eq!(source, before);
    params.set_callback(None);
    assert_eq!(
        source
            .enumerate_tautomers_with_params(&params)
            .unwrap()
            .len(),
        2
    );
}
struct Failure;
impl TautomerEnumerationCallback for Failure {
    fn should_continue(
        &self,
        _: TautomerMoleculeView<'_>,
        _: TautomerProgress<'_>,
    ) -> Result<bool, TautomerRunError> {
        Err(TautomerRunError::Callback("callback failure".into()))
    }
}
impl TautomerScorer for Failure {
    fn score(&self, _: TautomerMoleculeView<'_>) -> Result<i32, TautomerRunError> {
        Err(TautomerRunError::Callback("scorer failure".into()))
    }
}
#[test]
fn callback_and_scorer_errors_preserve_source_and_structured_error() {
    let source = Molecule::from_smiles("CC(C)=O").unwrap();
    let before = source.clone();
    let mut params = TautomerParams::default();
    params.set_callback(Some(Arc::new(Failure)));
    for error in [
        source.enumerate_tautomers_with_params(&params).unwrap_err(),
        source.canonical_tautomer_with_params(&params).unwrap_err(),
    ] {
        assert!(
            matches!(error,OperationError::Tautomer(TautomerRunError::Callback(ref value)) if value=="callback failure")
        );
    }
    params.set_callback(None);
    params.set_scorer(Some(Arc::new(Failure)));
    assert!(
        matches!(source.canonical_tautomer_with_params(&params),Err(OperationError::Tautomer(TautomerRunError::Callback(ref value))) if value=="scorer failure")
    );
    assert_eq!(source, before);
}
struct FavorEnol;
impl TautomerScorer for FavorEnol {
    fn score(&self, m: TautomerMoleculeView<'_>) -> Result<i32, TautomerRunError> {
        Ok(if m.to_smiles()? == "C=C(C)O" {
            100
        } else {
            -100
        })
    }
}
#[test]
fn custom_scorer_and_terms_use_canonical_public_configuration() {
    let source = Molecule::from_smiles("CC(C)=O").unwrap();
    let mut params = TautomerParams::default();
    params.set_scorer(Some(Arc::new(FavorEnol)));
    assert_eq!(
        source
            .canonical_tautomer_with_params(&params)
            .unwrap()
            .to_smiles()
            .unwrap(),
        "C=C(C)O"
    );
    assert_eq!(default_tautomer_score_terms().len(), 12);
    let score = source
        .tautomer_score_with_params(&TautomerScoreParams {
            terms: Some(vec![TautomerScoreTerm::new("carbonyl", "C=O", -7)]),
        })
        .unwrap();
    assert_eq!(score.substructure(), -7);
}
#[test]
fn limits_and_current_v1_custom_catalogs_keep_source_defaults() {
    let source = Molecule::from_smiles("CC(C)=O").unwrap();
    let current = TautomerParams::default();
    let v1 = TautomerParams::v1().unwrap();
    assert_eq!(current.transform_count(), 37);
    assert_eq!(v1.transform_count(), 36);
    assert_eq!(current.max_tautomers(), 1000);
    assert_eq!(current.max_transforms(), 1000);
    assert!(
        current.remove_sp3_stereo()
            && current.remove_bond_stereo()
            && current.remove_isotopic_hydrogens()
            && current.reassign_stereo()
    );
    assert_eq!(
        source
            .canonical_tautomer_with_params(&v1)
            .unwrap()
            .to_smiles()
            .unwrap(),
        "CC(C)=O"
    );
    let limited = current.clone().with_max_transforms(0);
    assert_eq!(
        source
            .enumerate_tautomers_with_params(&limited)
            .unwrap()
            .status(),
        TautomerEnumerationStatus::MaxTransformsReached
    );
    assert_eq!(current.max_transforms(), 1000);
    let empty = TautomerParams::from_transform_data(&[]).unwrap();
    assert_eq!(
        source
            .enumerate_tautomers_with_params(&empty)
            .unwrap()
            .len(),
        1
    );
}
