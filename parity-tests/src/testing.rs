//! Cargo stage: load all references once, compare independent registered tests.
//! No Python path, generator or preparation call is accepted here.
use crate::{Ready, Result, load_suite, registry, root, write_report};
use serde::Serialize;
use std::{
    path::PathBuf,
    sync::{Mutex, OnceLock},
};

struct Loaded {
    keys: Vec<String>,
    ready: Ready,
    reports: Mutex<Vec<TestSummary>>,
    report_dir: PathBuf,
    report_path: PathBuf,
}
static SUITE: OnceLock<Result<Loaded>> = OnceLock::new();

#[derive(Serialize)]
struct TestSummary {
    test: String,
    compared: usize,
    failed: usize,
    details: PathBuf,
}

pub fn run_registered(key: &str, test_keys: &[&str]) -> Result<()> {
    let registered: Vec<_> = registry::TASKS.iter().map(|task| task.key()).collect();
    if registered.iter().map(String::as_str).collect::<Vec<_>>() != test_keys {
        return Err("Cargo tests must match all registered function/corpus pairs".into());
    }
    let suite = SUITE
        .get_or_init(|| {
            let data = std::env::var_os("PARITY_DATA")
                .map(PathBuf::from)
                .unwrap_or_else(|| root().join("target/parity-tests"));
            let data = if data.is_absolute() {
                data
            } else {
                root().join(data)
            };
            // The complete selection is loaded before ANY test gets work.
            let (tasks, ready) = load_suite(&data).map_err(|e| {
                format!(
                    "{e}\nPrepare first: cargo run -p cosmolkit-parity-tests --release -- prepare --data {}",
                    data.display()
                )
            })?;
            let report_dir = tempfile::Builder::new()
                .prefix("run-")
                .tempdir_in(&data)
                .map_err(|e| e.to_string())?
                .keep();
            Ok(Loaded {
                keys: tasks.iter().map(|task| task.key()).collect(),
                ready,
                reports: Mutex::new(Vec::new()),
                report_dir,
                report_path: data.join("rust-report.json"),
            })
        })
        .as_ref()
        .map_err(Clone::clone)?;
    if !suite.keys.iter().any(|selected| selected == key) {
        return Err(format!(
            "{key}: not prepared; prepare this test or filter Cargo to the prepared test"
        ));
    }
    let report = suite.ready.compare_task(key);
    if report.is_empty() {
        return Err(format!("{key}: zero comparisons"));
    }
    let failed = report.iter().filter(|row| !row.matches).count();
    let rows = report.len();
    let details = suite.report_dir.join(format!("{key}.json"));
    write_report(&details, &report)?;
    // Preserve per-test details once; merge only small summaries. Do not
    // reserialize every earlier corpus result whenever another test finishes.
    let mut reports = suite.reports.lock().map_err(|_| "report lock poisoned")?;
    reports.push(TestSummary {
        test: key.into(),
        compared: rows,
        failed,
        details,
    });
    reports.sort_by(|a, b| a.test.cmp(&b.test));
    let summary = serde_json::json!({
        "expected_tests":suite.keys.len(), "finished_tests":reports.len(),
        "expected_cases":suite.ready.len(), "tests":&*reports,
    });
    let mut file = tempfile::NamedTempFile::new_in(&suite.report_dir).map_err(|e| e.to_string())?;
    std::io::Write::write_all(
        &mut file,
        &serde_json::to_vec_pretty(&summary).map_err(|e| e.to_string())?,
    )
    .map_err(|e| e.to_string())?;
    file.persist(&suite.report_path)
        .map_err(|e| e.to_string())?;
    println!("{key}: {}/{rows} matched", rows - failed);
    if failed != 0 {
        return Err(format!(
            "{key}: {failed}/{rows} mismatches; {}",
            suite.report_path.display()
        ));
    }
    Ok(())
}
