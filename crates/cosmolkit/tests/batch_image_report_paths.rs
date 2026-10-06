//! Independent p1/ROOT acceptance pending.
#[cfg(all(feature = "cap-batch", feature = "cap-depict"))]
use cosmolkit::BatchImageParams;
#[cfg(feature = "cap-batch")]
use cosmolkit::{BatchExportReport, BatchParams, MoleculeBatch};
use std::{cell::Cell, path::PathBuf};
#[path = "../../../python/src/user_path.rs"]
mod user_path;
struct Directory(PathBuf);
impl Directory {
    fn new() -> Self {
        use std::sync::atomic::{AtomicU64, Ordering};
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "cosmolkit-image-report-path-{}-{}",
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
#[cfg(all(feature = "cap-batch", feature = "cap-smiles", feature = "cap-depict"))]
fn batch() -> MoleculeBatch {
    MoleculeBatch::from_smiles_list(&["CCO".into()]).unwrap()
}
#[cfg(all(feature = "cap-batch", feature = "cap-smiles", feature = "cap-depict"))]
fn options(jobs: usize) -> BatchImageParams {
    BatchImageParams {
        format: "svg".into(),
        report_path: Some(PathBuf::from("~/report.json")),
        execution: BatchParams {
            n_jobs: Some(jobs),
            progress_bar: Some(false),
            ..Default::default()
        },
        ..Default::default()
    }
}
#[test]
#[cfg(feature = "cap-batch")]
fn canonical_report_receiver_preserves_exact_source_json_csv_and_literal_rust_paths() {
    let root = Directory::new();
    let report = BatchExportReport {
        written: 2,
        skipped: 3,
        failed: 1,
        errors: vec![],
    };
    for name in [
        "report.json",
        "report.JSON",
        "report",
        "report.txt",
        "report.csv",
        "report.CSV",
    ] {
        let path = root.0.join(name);
        report.write_report(&path).unwrap();
        let expected = if name.to_ascii_lowercase().ends_with(".csv") {
            "written,skipped,failed\n2,3,1\n"
        } else {
            "{\n  \"written\": 2,\n  \"skipped\": 3,\n  \"failed\": 1\n}\n"
        };
        assert_eq!(std::fs::read(path).unwrap(), expected.as_bytes());
    }
    let literal = root.0.join("~");
    std::fs::create_dir(&literal).unwrap();
    report.write_report(&literal.join("report.json")).unwrap();
    assert!(literal.join("report.json").is_file());
    assert_eq!((report.written, report.skipped, report.failed), (2, 3, 1));
}
#[test]
#[cfg(all(feature = "cap-batch", feature = "cap-smiles", feature = "cap-depict"))]
fn actual_image_output_precedes_report_home_lookup_and_uses_the_canonical_writer() {
    let root = Directory::new();
    let batch = batch();
    for jobs in [1, 4] {
        let params = options(jobs);
        let out = format!("~/images-{jobs}");
        let image = root.0.join(format!("images-{jobs}/mol_0.svg"));
        let reads = Cell::new(0);
        let result = user_path::with_image_user_paths(
            &out,
            Some("~/report.json"),
            || {
                reads.set(reads.get() + 1);
                if reads.get() == 2 {
                    assert!(image.is_file());
                }
                Some(root.0.clone().into_os_string())
            },
            |directory| {
                batch.to_images_with_params(
                    directory,
                    &BatchImageParams {
                        report_path: None,
                        ..params.clone()
                    },
                )
            },
            |path, report| report.write_report(path),
        )
        .unwrap();
        assert_eq!(reads.get(), 2);
        assert_eq!((result.written, result.skipped, result.failed), (1, 0, 0));
        assert_eq!(
            std::fs::read(&image).unwrap(),
            batch.to_list()[0]
                .as_ref()
                .unwrap()
                .to_svg(300, 300)
                .unwrap()
                .as_bytes()
        );
        assert_eq!(
            std::fs::read(root.0.join("report.json")).unwrap(),
            b"{\n  \"written\": 1,\n  \"skipped\": 0,\n  \"failed\": 0\n}\n"
        );
        assert_eq!(params.report_path, Some(PathBuf::from("~/report.json")));
        assert_eq!(params.execution.n_jobs, Some(jobs));
    }
}
#[test]
#[cfg(all(feature = "cap-batch", feature = "cap-smiles", feature = "cap-depict"))]
fn missing_report_home_retains_actual_images_and_export_failure_never_reads_report_home() {
    let root = Directory::new();
    let batch = batch();
    let params = options(1);
    let image_dir = root.0.join("images");
    let result = user_path::with_image_user_paths(
        image_dir.to_str().unwrap(),
        Some("~/report.json"),
        || {
            assert!(image_dir.join("mol_0.svg").is_file());
            None
        },
        |directory| {
            batch.to_images_with_params(
                directory,
                &BatchImageParams {
                    report_path: None,
                    ..params.clone()
                },
            )
        },
        |_, _| panic!("missing report HOME must not open report"),
    );
    assert!(matches!(
        result,
        Err(user_path::ImagePathError::Report(user_path::MissingHome))
    ));
    assert!(image_dir.join("mol_0.svg").is_file());
    assert!(!root.0.join("report.json").exists());
    let invalid = options(0);
    let invalid_dir = root.0.join("invalid-jobs");
    let result = user_path::with_image_user_paths(
        invalid_dir.to_str().unwrap(),
        Some("~/report.json"),
        || panic!("export failure must precede report HOME access"),
        |directory| {
            batch.to_images_with_params(
                directory,
                &BatchImageParams {
                    report_path: None,
                    ..invalid.clone()
                },
            )
        },
        |_, _| panic!("export failure must not open report"),
    );
    match result {
        Err(user_path::ImagePathError::Export(error)) => {
            assert_eq!(error.record_errors[0].operation, "n_jobs")
        }
        _ => panic!("expected original export configuration error"),
    }
    assert!(!invalid_dir.exists());
}
#[test]
#[cfg(all(feature = "cap-batch", feature = "cap-smiles", feature = "cap-depict"))]
fn report_io_failure_retains_real_cause_and_completed_images() {
    let root = Directory::new();
    let batch = batch();
    let params = options(1);
    let image_dir = root.0.join("images");
    let report_dir = root.0.join("report.json");
    std::fs::create_dir(&report_dir).unwrap();
    let result = user_path::with_image_user_paths(
        image_dir.to_str().unwrap(),
        Some("~/report.json"),
        || {
            assert!(image_dir.join("mol_0.svg").is_file());
            Some(root.0.clone().into_os_string())
        },
        |directory| {
            batch.to_images_with_params(
                directory,
                &BatchImageParams {
                    report_path: None,
                    ..params.clone()
                },
            )
        },
        |path, report| report.write_report(path),
    );
    match result {
        Err(user_path::ImagePathError::ReportWrite(error)) => {
            assert_eq!(error.errors, 1);
            assert_eq!(error.record_errors.len(), 1);
            let record = &error.record_errors[0];
            assert_eq!(record.index, 0);
            assert_eq!(record.operation, "write error report");
            let cause = std::error::Error::source(record)
                .unwrap()
                .downcast_ref::<std::io::Error>()
                .unwrap();
            assert_eq!(record.message, cause.to_string());
            assert!(cause.raw_os_error().is_some());
        }
        _ => panic!("expected post-image typed report IO error"),
    }
    assert!(image_dir.join("mol_0.svg").is_file());
    assert!(report_dir.is_dir());
}
