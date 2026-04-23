use core::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::process::Command;

use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct GoldenRecord {
    smiles: String,
    rdkit_ok: bool,
    atomic_nums: Option<Vec<u8>>,
    error: Option<String>,
}

#[derive(Debug)]
enum TestDataError {
    Io(std::io::Error),
    Json {
        line_no: usize,
        source: serde_json::Error,
    },
}

impl fmt::Display for TestDataError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(err) => write!(f, "{err}"),
            Self::Json { line_no, source } => {
                write!(f, "invalid jsonl at line {line_no}: {source}")
            }
        }
    }
}

impl std::error::Error for TestDataError {}

impl From<std::io::Error> for TestDataError {
    fn from(value: std::io::Error) -> Self {
        Self::Io(value)
    }
}

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("crates/")
        .parent()
        .expect("repo root")
        .to_path_buf()
}

fn load_smiles() -> Result<Vec<String>, TestDataError> {
    let path = repo_root().join("tests/smiles.txt");
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut rows = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        rows.push(trimmed.to_owned());
    }
    Ok(rows)
}

fn load_golden() -> Result<Vec<GoldenRecord>, TestDataError> {
    let path = repo_root().join("tests/golden/atomic_nums.jsonl");
    ensure_golden_exists(&path);
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut rows = Vec::new();

    for (idx, line) in reader.lines().enumerate() {
        let line_no = idx + 1;
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let record = serde_json::from_str::<GoldenRecord>(&line)
            .map_err(|source| TestDataError::Json { line_no, source })?;
        rows.push(record);
    }
    Ok(rows)
}

fn ensure_golden_exists(golden_path: &PathBuf) {
    if golden_path.exists() {
        return;
    }

    let repo = repo_root();
    let script = repo.join("tests/scripts/gen_rdkit_atomic_nums.py");
    let smiles = repo.join("tests/smiles.txt");

    let candidates = [
        std::env::var("COSMOLKIT_PYTHON").ok(),
        Some(repo.join(".venv/bin/python").display().to_string()),
        Some(String::from("python3")),
    ];

    let mut last_error = String::new();
    for candidate in candidates.iter().flatten() {
        let output = Command::new(candidate)
            .arg(&script)
            .arg("--input")
            .arg(&smiles)
            .arg("--output")
            .arg(golden_path)
            .output();

        match output {
            Ok(out) if out.status.success() => return,
            Ok(out) => {
                last_error = format!(
                    "python={} exit={} stderr={}",
                    candidate,
                    out.status,
                    String::from_utf8_lossy(&out.stderr)
                );
            }
            Err(err) => {
                last_error = format!("python={} spawn error={}", candidate, err);
            }
        }
    }

    panic!(
        "golden file missing and auto-generation failed.\n\
         expected: {}\n\
         tried COSMOLKIT_PYTHON, .venv/bin/python, python3.\n\
         last error: {}\n\
         please run:\n\
         uv sync --group dev && .venv/bin/python tests/scripts/gen_rdkit_atomic_nums.py --input tests/smiles.txt --output tests/golden/atomic_nums.jsonl",
        golden_path.display(),
        last_error
    );
}

#[test]
fn golden_has_entry_for_each_smiles_input() {
    let smiles = load_smiles().expect("should read tests/smiles.txt");
    let golden = load_golden().expect("should read tests/golden/atomic_nums.jsonl");
    assert_eq!(
        golden.len(),
        smiles.len(),
        "golden rows must match input smiles rows"
    );

    for (idx, (record, raw_smiles)) in golden.iter().zip(smiles.iter()).enumerate() {
        assert_eq!(
            record.smiles,
            *raw_smiles,
            "smiles mismatch at row {}",
            idx + 1
        );

        if record.rdkit_ok {
            assert!(
                record.atomic_nums.is_some(),
                "rdkit_ok=true requires atomic_nums at row {}",
                idx + 1
            );
            assert!(
                record.error.is_none(),
                "rdkit_ok=true should not carry error at row {}",
                idx + 1
            );
        } else {
            assert!(
                record.atomic_nums.is_none(),
                "rdkit_ok=false should not carry atomic_nums at row {}",
                idx + 1
            );
            assert!(
                record.error.is_some(),
                "rdkit_ok=false should carry error at row {}",
                idx + 1
            );
        }
    }
}

#[test]
fn atomic_numbers_match_rdkit_golden() {
    let golden = load_golden().expect("should read tests/golden/atomic_nums.jsonl");

    for (idx, record) in golden.iter().enumerate() {
        let ours = cosmolkit_chem_core::Molecule::from_smiles(&record.smiles)
            .map(|mol| mol.atomic_numbers());
        match (record.rdkit_ok, ours) {
            (true, Ok(actual)) => {
                let expected = record.atomic_nums.as_ref().expect("atomic nums missing");
                assert_eq!(
                    &actual,
                    expected,
                    "atomic number mismatch at row {} ({})",
                    idx + 1,
                    record.smiles
                );
            }
            (false, Err(_)) => {}
            (true, Err(err)) => {
                panic!(
                    "unexpected parse error at row {} ({}): {}",
                    idx + 1,
                    record.smiles,
                    err
                );
            }
            (false, Ok(actual)) => {
                panic!(
                    "expected parse failure but got {:?} at row {} ({})",
                    actual,
                    idx + 1,
                    record.smiles
                );
            }
        }
    }
}
