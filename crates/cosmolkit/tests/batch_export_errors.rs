#![cfg(all(feature = "cap-batch", feature = "cap-smiles"))]

use cosmolkit::{BatchErrorMode, BatchExportParams, BatchParams, MoleculeBatch};
use std::error::Error;
use std::path::PathBuf;

struct Directory(PathBuf);
impl Directory {
    fn new() -> Self {
        use std::sync::atomic::{AtomicU64, Ordering};
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "cosmolkit-export-errors-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&path).unwrap();
        Self(path)
    }
}
impl Drop for Directory {
    fn drop(&mut self) {
        std::fs::remove_dir_all(&self.0).unwrap();
    }
}

fn batch(inputs: &[&str]) -> MoleculeBatch {
    MoleculeBatch::from_smiles_list_with_params(
        &inputs
            .iter()
            .map(|input| (*input).into())
            .collect::<Vec<_>>(),
        &Default::default(),
        &BatchParams {
            errors: Some(BatchErrorMode::KeepErrors),
            ..Default::default()
        },
    )
    .unwrap()
}

#[test]
fn sdf_reports_preserve_original_errors_and_count_every_record_once() {
    let root = Directory::new();
    for (case, inputs) in [vec!["CCO", "C1CC", "O"], vec!["C1CC", "N1"], vec![]]
        .iter()
        .enumerate()
    {
        let batch = batch(inputs);
        let original = batch.errors();
        for jobs in [1, 2] {
            let params = BatchExportParams {
                errors: Some(BatchErrorMode::KeepErrors),
                n_jobs: Some(jobs),
                ..Default::default()
            };
            let single = root.0.join(format!("single-{case}-{jobs}.sdf"));
            let reports = [
                batch
                    .write_sdf_with_params(single.to_str().unwrap(), &params, None)
                    .unwrap(),
                batch
                    .write_sdf_files_with_params(
                        root.0
                            .join(format!("files-{case}-{jobs}"))
                            .to_str()
                            .unwrap(),
                        &params,
                        None,
                        None,
                    )
                    .unwrap(),
            ];
            for report in reports {
                assert_eq!(report.total(), inputs.len());
                assert_eq!(report.success(), batch.valid_count());
                assert_eq!(report.failed(), original.len());
                assert_eq!(report.failed(), report.errors().len());
                for (actual, expected) in report.errors().iter().zip(&original) {
                    assert_eq!(actual.index, expected.index);
                    assert_eq!(actual.operation, expected.operation);
                    assert_eq!(actual.message, expected.message);
                    assert!(std::ptr::eq(
                        actual.source().unwrap(),
                        expected.source().unwrap()
                    ));
                }
            }
            assert_eq!(
                std::fs::read_to_string(single)
                    .unwrap()
                    .matches("$$$$")
                    .count(),
                batch.valid_count()
            );
        }
        assert_eq!(batch.errors().len(), original.len());
    }
}

#[test]
fn sdf_file_reports_merge_existing_and_new_errors_in_input_order() {
    let root = Directory::new();
    let batch = batch(&["CCO", "C1CC", "O"]);
    let original = batch.errors();
    let names = [Some("blocked.sdf".into()), None, Some("water.sdf".into())];
    std::fs::create_dir(root.0.join("blocked.sdf")).unwrap();
    for mode in [BatchErrorMode::KeepErrors, BatchErrorMode::Strict] {
        let result = batch.write_sdf_files_with_params(
            root.0.to_str().unwrap(),
            &BatchExportParams {
                errors: Some(mode),
                n_jobs: Some(2),
                ..Default::default()
            },
            Some(&names),
            None,
        );
        let errors = match mode {
            BatchErrorMode::KeepErrors => {
                let report = result.unwrap();
                assert_eq!(
                    (report.total(), report.success(), report.failed()),
                    (3, 1, 2)
                );
                report.errors
            }
            BatchErrorMode::Strict => result.unwrap_err().record_errors,
        };
        assert_eq!(errors.len(), 2);
        assert_eq!(
            (errors[0].index, errors[0].operation),
            (0, "batch.write_sdf_files")
        );
        assert!(errors[0].source().unwrap().is::<std::io::Error>());
        assert_eq!(
            (errors[1].index, errors[1].operation),
            (1, original[0].operation)
        );
        assert_eq!(errors[1].message, original[0].message);
        assert!(std::ptr::eq(
            errors[1].source().unwrap(),
            original[0].source().unwrap()
        ));
    }
    assert!(root.0.join("water.sdf").is_file());
}

#[test]
fn strict_single_sdf_preserves_errors_without_opening_output() {
    let root = Directory::new();
    let batch = batch(&["CCO", "C1CC"]);
    let path = root.0.join("must-not-exist.sdf");
    let error = batch
        .write_sdf_with_params(
            path.to_str().unwrap(),
            &BatchExportParams {
                errors: Some(BatchErrorMode::Strict),
                ..Default::default()
            },
            None,
        )
        .unwrap_err();
    assert_eq!(error.record_errors.len(), 1);
    assert_eq!(error.record_errors[0].index, 1);
    assert_eq!(error.record_errors[0].operation, "batch.from_smiles_list");
    assert!(!path.exists());
}
