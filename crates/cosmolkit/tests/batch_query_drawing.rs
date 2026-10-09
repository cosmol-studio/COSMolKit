#![cfg(all(feature = "cap-batch", feature = "cap-smiles", feature = "cap-depict"))]
//! Source-backed delivery proposals; p1/ROOT alone decide acceptance conditions.
use cosmolkit::{
    BatchErrorMode, BatchImageParams, BatchParams, BatchQueryParams, MoleculeBatch,
    SmilesParseParams, SmilesWriteParams,
};
fn batch() -> MoleculeBatch {
    MoleculeBatch::from_smiles_list_with_params(
        &["CCO".into(), "C1".into(), "c1ccccc1".into()],
        &SmilesParseParams::default(),
        &BatchParams {
            errors: BatchErrorMode::KeepErrors,
            ..Default::default()
        },
    )
    .unwrap()
    .with_parallel_jobs(Some(2))
    .unwrap()
    .with_progress_bar(Some(false))
}
struct Directory(std::path::PathBuf);
impl Directory {
    fn new() -> Self {
        use std::sync::atomic::{AtomicU64, Ordering};
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "cosmolkit-batch-proposal-{}-{}",
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
#[test]
#[cfg(feature = "cap-conformer")]
fn complete_query_values_and_invalid_positions_match_scalar_serial_and_parallel() {
    let batch = batch();
    let originals = batch.to_list();
    for jobs in [1, 2, 4] {
        let query = BatchQueryParams {
            n_jobs: Some(jobs),
            progress_bar: Some(false),
            ..Default::default()
        };
        let strings = batch
            .to_smiles_list_with_params(&SmilesWriteParams::default(), &query)
            .unwrap();
        assert_eq!(
            strings,
            vec![Some("CCO".into()), None, Some("c1ccccc1".into())]
        );
        let bounds = batch.dg_bounds_matrix_list_with_params(&query).unwrap();
        let svg = batch.to_svg_list_with_params(300, 300, &query).unwrap();
        for index in 0..3 {
            match &originals[index] {
                Some(molecule) => {
                    assert_eq!(
                        bounds[index].as_ref().unwrap(),
                        &molecule.dg_bounds_matrix().unwrap()
                    );
                    assert_eq!(
                        svg[index].as_ref().unwrap(),
                        &molecule.to_svg(300, 300).unwrap()
                    );
                }
                None => {
                    assert!(bounds[index].is_none());
                    assert!(svg[index].is_none());
                }
            }
        }
        assert_eq!(batch.parallel_jobs(), Some(2));
        assert_eq!(batch.progress_bar(), Some(false));
    }
    let query = BatchQueryParams {
        n_jobs: Some(0),
        progress_bar: Some(false),
        ..Default::default()
    };
    assert!(
        batch
            .to_smiles_list_with_params(&Default::default(), &query)
            .is_err()
    );
}
#[test]
#[cfg(feature = "cap-depict")]
fn image_exports_preserve_original_filename_rules_complete_bytes_and_reports() {
    let batch = batch();
    let directory = Directory::new();
    for (format, jobs) in [("svg", 1), ("svg", 4), ("png", 2)] {
        let out = directory.0.join(format!("{format}-{jobs}"));
        let report_path = directory.0.join(format!("{format}-{jobs}.json"));
        let params = BatchImageParams {
            format: format.into(),
            execution: BatchParams {
                errors: BatchErrorMode::KeepErrors,
                n_jobs: Some(jobs),
                progress_bar: Some(false),
            },
            filenames: Some(vec![Some(format!(" ethanol.{format} ")), None, None]),
            report_path: Some(report_path.clone()),
            ..Default::default()
        };
        let report = batch.write_images_with_params(&out, &params).unwrap();
        assert_eq!(
            (
                report.total(),
                report.success(),
                report.skipped,
                report.failed()
            ),
            (3, 2, 1, 0)
        );
        assert!(report.errors().is_empty());
        assert_eq!(
            std::fs::read_to_string(report_path).unwrap(),
            "{\n  \"written\": 2,\n  \"skipped\": 1,\n  \"failed\": 0\n}\n"
        );
        let values = batch.to_list();
        for (index, name) in [
            (0, format!("ethanol.{format}")),
            (2, format!("mol_2.{format}")),
        ] {
            let molecule = values[index].as_ref().unwrap();
            let bytes = std::fs::read(out.join(name)).unwrap();
            let expected = if format == "svg" {
                molecule.to_svg(300, 300).unwrap().into_bytes()
            } else {
                molecule.to_png(300, 300).unwrap()
            };
            assert_eq!(bytes, expected);
        }
        assert!(!out.join(format!("mol_1.{format}")).exists());
    }
    for names in [
        vec![Some("../escape".into()), None, None],
        vec![Some("duplicate".into()), None, Some("duplicate.svg".into())],
        vec![Some("wrong.png".into()), None, None],
    ] {
        assert!(
            batch
                .write_images_with_params(
                    &directory.0.join("invalid"),
                    &BatchImageParams {
                        format: "svg".into(),
                        filenames: Some(names),
                        ..Default::default()
                    }
                )
                .is_err()
        );
    }
}
#[test]
#[cfg(feature = "cap-depict")]
fn strict_export_finishes_valid_records_and_orders_new_failures_before_existing_errors() {
    let batch = batch();
    let directory = Directory::new();
    let out = directory.0.join("strict-svg");
    let strict = batch
        .write_images_with_params(
            &out,
            &BatchImageParams {
                format: "svg".into(),
                ..Default::default()
            },
        )
        .unwrap_err();
    assert_eq!(
        strict
            .record_errors
            .iter()
            .map(|e| e.index)
            .collect::<Vec<_>>(),
        vec![1]
    );
    assert!(out.join("mol_0.svg").exists());
    assert!(out.join("mol_2.svg").exists());
    let strict = batch
        .write_images_with_params(
            &directory.0.join("bad-format"),
            &BatchImageParams {
                format: "SVG".into(),
                ..Default::default()
            },
        )
        .unwrap_err();
    assert_eq!(
        strict
            .record_errors
            .iter()
            .map(|e| e.index)
            .collect::<Vec<_>>(),
        vec![0, 2, 1]
    );
    assert_eq!(
        strict
            .record_errors
            .iter()
            .map(|e| e.operation)
            .collect::<Vec<_>>(),
        vec![
            "batch.write_images",
            "batch.write_images",
            "batch.from_smiles_list"
        ]
    );
    let kept = batch
        .write_images_with_params(
            &directory.0.join("keep-format"),
            &BatchImageParams {
                format: "SVG".into(),
                execution: BatchParams {
                    errors: BatchErrorMode::KeepErrors,
                    ..Default::default()
                },
                ..Default::default()
            },
        )
        .unwrap();
    assert_eq!((kept.written, kept.skipped, kept.failed), (0, 1, 2));
    let file = directory.0.join("file");
    std::fs::write(&file, "cannot create directory over file").unwrap();
    let error = batch.write_images(&file).unwrap_err();
    assert!(
        std::error::Error::source(&error.record_errors[0])
            .unwrap()
            .is::<std::io::Error>()
    );
}
