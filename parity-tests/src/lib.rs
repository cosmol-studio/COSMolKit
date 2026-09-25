//! Small Rust-owned parity pilot; no performance or binding claims.
pub mod execute;
pub mod registry;

use registry::{Input, Pair, RDKIT_VERSION, Record, Task};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::{
    fs,
    io::Write,
    path::{Path, PathBuf},
    process::{Command, Stdio},
};

type Result<T> = std::result::Result<T, String>;
const ORACLE: &str = include_str!("../../tools/oracles/rdkit/fingerprint_values_pilot.py");

pub fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .to_path_buf()
}

fn digest(bytes: &[u8]) -> String {
    Sha256::digest(bytes)
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect()
}
fn encode<T: Serialize>(value: &T) -> Result<Vec<u8>> {
    serde_json::to_vec_pretty(value).map_err(|e| e.to_string())
}
fn read(path: &Path) -> Result<Vec<u8>> {
    fs::read(path).map_err(|e| format!("{}: {e}", path.display()))
}

pub fn corpus(path: Option<&Path>) -> Result<Vec<Pair>> {
    match path {
        None => Ok(registry::builtin()),
        Some(p) => serde_json::from_slice(&read(p)?).map_err(|e| format!("{}: {e}", p.display())),
    }
}

#[derive(Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
struct Manifest {
    schema: u32,
    task: String,
    rdkit_version: String,
    registry_sha256: String,
    oracle_sha256: String,
    input_sha256: String,
    reference_sha256: String,
    rows: usize,
}

fn identity(task: &Task, input: &[u8], reference: &[u8], rows: usize) -> Manifest {
    Manifest {
        schema: 1,
        task: task.operation.name().into(),
        rdkit_version: RDKIT_VERSION.into(),
        registry_sha256: digest(include_bytes!("registry.rs")),
        oracle_sha256: digest(ORACLE.as_bytes()),
        input_sha256: digest(input),
        reference_sha256: digest(reference),
        rows,
    }
}

fn check_records(inputs: &[Input], records: &[Record]) -> Result<()> {
    if inputs.len() != records.len() {
        return Err("reference row count mismatch".into());
    }
    for (input, record) in inputs.iter().zip(records) {
        if input != &record.input {
            return Err("reference case/operation/width mismatch".into());
        }
        let entries = &record.output.entries;
        if record.output.length != input.case.length
            || entries.windows(2).any(|w| w[0].0 >= w[1].0)
            || entries.iter().any(|&(key, _)| key >= input.case.length)
        {
            return Err(format!("{}: malformed reference output", input.case.id));
        }
    }
    Ok(())
}

fn oracle(inputs: &[Input], python: &Path) -> Result<Vec<Record>> {
    // Python is solely the pinned RDKit adapter, never the CK executor.
    let mut child = Command::new(python)
        .args(["-c", ORACLE, RDKIT_VERSION])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| format!("RDKit adapter: {e}"))?;
    let payload = encode(&inputs)?;
    let mut stdin = child.stdin.take().ok_or("oracle stdin unavailable")?;
    // Write concurrently so a larger corpus cannot deadlock stdout/stderr pipes.
    let writer = std::thread::spawn(move || stdin.write_all(&payload));
    let output = child.wait_with_output().map_err(|e| e.to_string())?;
    let sent = writer.join().map_err(|_| "oracle input writer panicked")?;
    if !output.status.success() {
        return Err(format!(
            "RDKit adapter failed: {}",
            String::from_utf8_lossy(&output.stderr)
        ));
    }
    sent.map_err(|e| e.to_string())?;
    let records =
        serde_json::from_slice(&output.stdout).map_err(|e| format!("oracle output: {e}"))?;
    Ok(records)
}

#[derive(Debug, Default, PartialEq, Eq)]
pub struct Preparation {
    pub reused_tasks: usize,
    pub generated_tasks: usize,
    pub rows: usize,
}

/// Reuse verified references and prepare missing or invalid generations.
/// Ordinary cargo tests inject a synthetic generator, never the real oracle.
pub fn prepare(tasks: &[&Task], cases: &[Pair], data: &Path, python: &Path) -> Result<Preparation> {
    ensure_ready_with(tasks, cases, data, |inputs| oracle(inputs, python))
        .map(|(_, preparation)| preparation)
}

fn ensure_ready_with(
    tasks: &[&Task],
    cases: &[Pair],
    data: &Path,
    mut generate: impl FnMut(&[Input]) -> Result<Vec<Record>>,
) -> Result<(Ready, Preparation)> {
    registry::validate(cases, tasks)?;
    if tasks.is_empty() || tasks.iter().any(|task| task.widths.is_empty()) {
        return Err("empty task/width selection; 0 Rust operation calls".into());
    }
    fs::create_dir_all(data).map_err(|e| e.to_string())?;
    // Serialize publishers; the OS releases this lock even after interruption.
    let lock = fs::OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .truncate(false)
        .open(data.join(".prepare.lock"))
        .map_err(|e| e.to_string())?;
    lock.lock().map_err(|e| e.to_string())?;
    // Inspect the entire selection before invoking any generator or CK operation.
    let missing: Vec<_> = tasks
        .iter()
        .map(|task| match preflight(&[task], cases, data) {
            Ok(_) => false,
            Err(error) => {
                eprintln!("Preparing {}: {error}", task.operation.name());
                true
            }
        })
        .collect();
    let mut preparation = Preparation::default();
    for (task, needs_data) in tasks.iter().zip(missing) {
        if needs_data {
            let inputs = registry::expand(cases, task);
            let records = generate(&inputs)
                .and_then(|records| {
                    check_records(&inputs, &records)?;
                    Ok(records)
                })
                .map_err(|e| {
                    format!(
                        "{}: preparation failed; 0 Rust operation calls: {e}",
                        task.operation.name()
                    )
                })?;
            publish(data, task, &inputs, &records)?;
            preparation.generated_tasks += 1;
        } else {
            preparation.reused_tasks += 1;
        }
    }
    // Global barrier: never interleave preparation with chemistry execution.
    let ready = preflight(tasks, cases, data)?;
    preparation.rows = ready.len();
    Ok((ready, preparation))
}

fn publish(data: &Path, task: &Task, inputs: &[Input], records: &[Record]) -> Result<()> {
    let input = encode(&inputs)?;
    let reference = encode(&records)?;
    let manifest = encode(&identity(task, &input, &reference, inputs.len()))?;
    let destination = generation(data, task, &input);
    let temporary = tempfile::Builder::new()
        .prefix(".prepare-")
        .tempdir_in(data)
        .map_err(|e| e.to_string())?;
    for (name, bytes) in [
        ("input.json", input),
        ("reference.json", reference),
        ("manifest.json", manifest),
    ] {
        fs::write(temporary.path().join(name), bytes).map_err(|e| e.to_string())?;
    }
    // Preserve corrupt evidence; never delete the previous generation to repair it.
    let backup = destination.exists().then(|| {
        data.join(format!(
            ".invalid-{}-{}",
            task.operation.name(),
            temporary.path().file_name().unwrap().to_string_lossy()
        ))
    });
    if let Some(backup) = &backup {
        fs::rename(&destination, backup).map_err(|e| e.to_string())?;
        eprintln!("Invalid generation preserved at {}", backup.display());
    }
    if let Err(error) = fs::rename(temporary.path(), &destination) {
        if let Some(backup) = &backup {
            fs::rename(backup, &destination).map_err(|rollback| {
                format!(
                    "publish: {error}; rollback: {rollback}; preserved at {}",
                    backup.display()
                )
            })?;
        }
        return Err(format!("publish: {error}"));
    }
    Ok(())
}

fn generation(data: &Path, task: &Task, input: &[u8]) -> PathBuf {
    let identity = format!(
        "{}{}{}",
        digest(input),
        digest(ORACLE.as_bytes()),
        digest(include_bytes!("registry.rs"))
    );
    data.join(format!(
        "{}-{}",
        task.operation.name(),
        digest(identity.as_bytes())
    ))
}

pub struct Ready {
    records: Vec<Record>,
}

impl Ready {
    pub fn len(&self) -> usize {
        self.records.len()
    }
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }
}

/// Load and validate EVERY selected task before returning any executable work.
/// Owned snapshots prevent input changes after preflight affecting this run.
pub fn preflight(tasks: &[&Task], cases: &[Pair], data: &Path) -> Result<Ready> {
    registry::validate(cases, tasks)?;
    let mut all = Vec::new();
    let mut errors = Vec::new();
    for task in tasks {
        let load = || -> Result<Vec<Record>> {
            let inputs = registry::expand(cases, task);
            let expected_input = encode(&inputs)?;
            let dir = generation(data, task, &expected_input);
            let input = read(&dir.join("input.json"))?;
            let reference = read(&dir.join("reference.json"))?;
            let manifest: Manifest = serde_json::from_slice(&read(&dir.join("manifest.json"))?)
                .map_err(|e| e.to_string())?;
            if input != expected_input
                || manifest != identity(task, &input, &reference, inputs.len())
            {
                return Err("stale or corrupted manifest/input/reference".into());
            }
            let records: Vec<Record> =
                serde_json::from_slice(&reference).map_err(|e| e.to_string())?;
            check_records(&inputs, &records)?;
            Ok(records)
        };
        match load() {
            Ok(records) => all.extend(records),
            Err(e) => errors.push(format!("{}: {e}", task.operation.name())),
        }
    }
    if !errors.is_empty() {
        return Err(format!(
            "preflight failed; 0 Rust operation calls\n{}",
            errors.join("\n")
        ));
    }
    Ok(Ready { records: all })
}

#[derive(Debug, Serialize)]
pub struct Comparison {
    pub input: Input,
    pub expected: registry::Value,
    pub actual: std::result::Result<registry::Value, String>,
    pub matches: bool,
}

/// Type equality compares all fields; no per-operation hand-written diff rules.
pub fn compare(ready: Ready, mut run: impl FnMut(&Input) -> Result<Record>) -> Vec<Comparison> {
    ready
        .records
        .into_iter()
        .map(|reference| {
            let actual = run(&reference.input).and_then(|record| {
                if record.input != reference.input {
                    return Err("executor changed case identity".into());
                }
                Ok(record.output)
            });
            let matches = actual.as_ref() == Ok(&reference.output);
            Comparison {
                input: reference.input,
                expected: reference.output,
                actual,
                matches,
            }
        })
        .collect()
}

pub fn run(
    tasks: &[&Task],
    cases: &[Pair],
    data: &Path,
    python: &Path,
    executor: impl FnMut(&Input) -> Result<Record>,
) -> Result<Vec<Comparison>> {
    run_with(
        tasks,
        cases,
        data,
        |inputs| oracle(inputs, python),
        executor,
    )
}

fn run_with(
    tasks: &[&Task],
    cases: &[Pair],
    data: &Path,
    generate: impl FnMut(&[Input]) -> Result<Vec<Record>>,
    executor: impl FnMut(&Input) -> Result<Record>,
) -> Result<Vec<Comparison>> {
    let (ready, preparation) = ensure_ready_with(tasks, cases, data, generate)?;
    eprintln!(
        "Preparation complete: {} reused tasks, {} generated tasks; {} ready cases",
        preparation.reused_tasks, preparation.generated_tasks, preparation.rows
    );
    Ok(compare(ready, executor))
}

pub fn write_report(path: &Path, report: &[Comparison]) -> Result<()> {
    fs::write(path, encode(&report)?).map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests;
