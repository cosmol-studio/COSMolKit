//! Literal filesystem paths and post-image report handling through the Rust API.
#[cfg(all(feature = "cap-batch", feature = "cap-depict"))]
use cosmolkit::BatchImageParams;
#[cfg(feature = "cap-batch")]
use cosmolkit::{BatchExportReport, BatchParams, MoleculeBatch};
use std::path::PathBuf;
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
        report_path: None,
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
fn canonical_report_receiver_preserves_two_count_json_csv_and_literal_rust_paths() {
    let root = Directory::new();
    let report = BatchExportReport {
        written: 2,
        failed: 4,
        errors: [1, 3, 4, 5]
            .map(|index| cosmolkit::BatchError::new(index, "fixture", "failed record"))
            .into(),
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
            "written,failed\n2,4\n"
        } else {
            "{\n  \"written\": 2,\n  \"failed\": 4\n}\n"
        };
        assert_eq!(std::fs::read(path).unwrap(), expected.as_bytes());
    }
    let literal = root.0.join("~");
    std::fs::create_dir(&literal).unwrap();
    report.write_report(&literal.join("report.json")).unwrap();
    assert!(literal.join("report.json").is_file());
    assert_eq!((report.written, report.failed, report.total()), (2, 4, 6));
}
#[test]
#[cfg(all(feature = "cap-batch", feature = "cap-smiles", feature = "cap-depict"))]
fn literal_image_paths_use_the_canonical_writer_and_preserve_options() {
    let root = Directory::new();
    let batch = batch();
    for jobs in [1, 4] {
        let out = root.0.join("~").join(format!("images-{jobs}"));
        let report_path = root.0.join("~").join("report.json");
        let params = BatchImageParams {
            report_path: Some(report_path.clone()),
            ..options(jobs)
        };
        let image = out.join("mol_0.svg");
        let result = batch.write_images_with_params(&out, &params).unwrap();
        assert_eq!((result.written, result.failed), (1, 0));
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
            std::fs::read(&report_path).unwrap(),
            b"{\n  \"written\": 1,\n  \"failed\": 0\n}\n"
        );
        assert_eq!(params.report_path, Some(report_path));
        assert_eq!(params.execution.n_jobs, Some(jobs));
    }
}
#[test]
#[cfg(all(feature = "cap-batch", feature = "cap-smiles", feature = "cap-depict"))]
fn missing_report_parent_retains_images_and_invalid_jobs_precede_report_write() {
    let root = Directory::new();
    let batch = batch();
    let image_dir = root.0.join("images");
    let report_path = root.0.join("missing-parent/report.json");
    let params = BatchImageParams {
        report_path: Some(report_path.clone()),
        ..options(1)
    };
    let error = batch
        .write_images_with_params(&image_dir, &params)
        .unwrap_err();
    assert_eq!(error.record_errors[0].operation, "write error report");
    let cause = std::error::Error::source(&error.record_errors[0])
        .unwrap()
        .downcast_ref::<std::io::Error>()
        .unwrap();
    assert_eq!(cause.kind(), std::io::ErrorKind::NotFound);
    assert!(image_dir.join("mol_0.svg").is_file());
    assert!(!report_path.exists());
    let valid_report_path = root.0.join("report.json");
    let invalid = BatchImageParams {
        report_path: Some(valid_report_path.clone()),
        ..options(0)
    };
    let invalid_dir = root.0.join("invalid-jobs");
    let error = batch
        .write_images_with_params(&invalid_dir, &invalid)
        .unwrap_err();
    assert_eq!(error.record_errors[0].operation, "n_jobs");
    assert!(!invalid_dir.exists());
    assert!(!valid_report_path.exists());
}
#[test]
#[cfg(all(feature = "cap-batch", feature = "cap-smiles", feature = "cap-depict"))]
fn report_io_failure_retains_real_cause_and_completed_images() {
    let root = Directory::new();
    let batch = batch();
    let image_dir = root.0.join("images");
    let report_dir = root.0.join("report.json");
    std::fs::create_dir(&report_dir).unwrap();
    let params = BatchImageParams {
        report_path: Some(report_dir.clone()),
        ..options(1)
    };
    let result = batch.write_images_with_params(&image_dir, &params);
    match result {
        Err(error) => {
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
