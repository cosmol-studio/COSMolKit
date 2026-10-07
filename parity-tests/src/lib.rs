//! Small Rust-owned parity pilot; no performance or binding claims.
mod bio_reference;
pub mod descriptor_reference;
mod draw_reference;
pub mod execute;
pub mod molecular;
mod native_draw_reference;
pub mod registry;
pub mod search;
pub mod special_regression;
pub mod tautomer;
pub mod tautomer_reference;
pub mod testing;
pub mod uff;

use registry::{Corpus, CorpusType, Input, RDKIT_VERSION, Record, Task};
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
const REFERENCE_PIN: &str = include_str!("../../testdata/reference/rdkit.json");
const GEMMI_PIN: &str = include_str!("../../testdata/reference/gemmi.json");

fn reference_pin(task: &Task) -> &'static str {
    if bio_reference::handles(task) {
        GEMMI_PIN
    } else {
        REFERENCE_PIN
    }
}

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

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CorpusSource {
    pub corpus_type: CorpusType,
    pub path: PathBuf,
}

/// Sources are explicitly typed. A suffix never selects a loader.
pub fn corpus(sources: &[CorpusSource], tasks: &[&Task]) -> Result<Corpus> {
    let mut seen = std::collections::BTreeSet::new();
    for source in sources {
        if !seen.insert(source.corpus_type.name()) {
            return Err("duplicate corpus type".into());
        }
        if !tasks.iter().any(|t| t.corpus_type == source.corpus_type) {
            return Err("selected tasks do not consume the supplied corpus input family".into());
        }
    }
    let mut corpus = Corpus::default();
    let mut loaded_bio = std::collections::BTreeSet::new();
    for task in tasks {
        let source = sources.iter().find(|s| s.corpus_type == task.corpus_type);
        match task.corpus_type {
            CorpusType::Smiles if corpus.molecules.is_empty() => {
                let default = root().join("testdata/smiles/corpus/smiles_small.smi");
                corpus.molecules =
                    molecular::read_corpus(source.map_or(default.as_path(), |s| &s.path))?;
            }
            CorpusType::FingerprintPairs if corpus.fingerprints.is_empty() => {
                corpus.fingerprints = if let Some(source) = source {
                    serde_json::from_slice(&read(&source.path)?).map_err(|e| e.to_string())?
                } else {
                    registry::fingerprint_corpus::generate()
                };
            }
            CorpusType::Smiles | CorpusType::FingerprintPairs => {}
            CorpusType::Pdb | CorpusType::Cif => {
                if loaded_bio.insert(task.corpus_type.name()) {
                    let rows: Vec<registry::BioPdbCase> = if let Some(source) = source {
                        serde_json::from_slice(&read(&source.path)?)
                            .map_err(|e| format!("{}: {e}", source.path.display()))?
                    } else {
                        // Like smiles_small, defaults are existing fixed inputs,
                        // never a claimed 5000-case biological corpus. The task's
                        // declared family selects the format, not the suffix.
                        let (format, path) = if task.corpus_type == CorpusType::Pdb {
                            (
                                registry::BioPdbCorpusFormat::Pdb,
                                "testdata/bio/fixtures/gemmi_full_feature_sample.pdb",
                            )
                        } else {
                            (
                                registry::BioPdbCorpusFormat::Cif,
                                "testdata/bio/fixtures/gemmi_full_feature_sample.cif",
                            )
                        };
                        vec![registry::BioPdbCase {
                            id: path.into(),
                            text: String::from_utf8(read(&root().join(path))?)
                                .map_err(|e| e.to_string())?,
                            format,
                        }]
                    };
                    if rows.iter().any(|row| !row.matches_corpus(task.corpus_type)) {
                        return Err(
                            "BIO corpus row format does not match its explicit input family".into(),
                        );
                    }
                    corpus.bio_cases.extend(rows);
                }
            }
            kind => return Err(format!("unimplemented corpus loader: {}", kind.name())),
        }
    }
    Ok(corpus)
}

fn registry_digest() -> String {
    digest(
        concat!(
            include_str!("registry.rs"),
            include_str!("registry/molecule_plan.rs"),
            include_str!("registry/fingerprint_corpus.rs"),
            include_str!("registry/fingerprint.rs"),
            include_str!("molecular.rs"),
            include_str!("tautomer.rs"),
            include_str!("uff.rs"),
            include_str!("search.rs")
        )
        .as_bytes(),
    )
}

#[derive(Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
struct Manifest {
    schema: u32,
    task: String,
    corpus_type: CorpusType,
    generator: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    rdkit_version: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    gemmi_version: Option<String>,
    reference_pin_sha256: String,
    reference_platform: String,
    registry_sha256: String,
    oracle_sha256: String,
    input_sha256: String,
    reference_sha256: String,
    rows: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    imported_reference: Option<serde_json::Value>,
}

fn adapter_digest(task: &Task) -> String {
    if bio_reference::handles(task) {
        digest(
            concat!(
                include_str!("bio_reference.rs"),
                include_str!("../../tools/oracles/gemmi/bio_pdb_values.py"),
                include_str!("../../tools/testdata/gemmi/pdb_coordinate_oracle.cpp"),
                include_str!("../../tools/testdata/gemmi/to_pdb_probe.cpp")
            )
            .as_bytes(),
        )
    } else if tautomer_reference::handles(task) {
        digest(include_str!("tautomer_reference.rs").as_bytes())
    } else if descriptor_reference::handles(task) {
        digest(include_str!("descriptor_reference.rs").as_bytes())
    } else if draw_reference::handles(task) {
        digest(include_str!("draw_reference.rs").as_bytes())
    } else if native_draw_reference::handles(task) {
        digest(include_str!("native_draw_reference.rs").as_bytes())
    } else {
        digest(ORACLE.as_bytes())
    }
}

fn identity(task: &Task, input: &[u8], reference: &[u8], rows: usize) -> Manifest {
    Manifest {
        schema: 3,
        task: task.key(),
        corpus_type: task.corpus_type,
        generator: task.generator.into(),
        rdkit_version: (!bio_reference::handles(task)).then(|| RDKIT_VERSION.into()),
        gemmi_version: bio_reference::handles(task)
            .then(|| registry::BIO_PDB_REFERENCE.version.into()),
        reference_pin_sha256: digest(reference_pin(task).as_bytes()),
        reference_platform: format!("{}-{}", std::env::consts::ARCH, std::env::consts::OS),
        registry_sha256: registry_digest(),
        oracle_sha256: adapter_digest(task),
        input_sha256: digest(input),
        reference_sha256: digest(reference),
        rows,
        imported_reference: bio_reference::provenance(task)
            .or_else(|| tautomer_reference::provenance(task))
            .or_else(|| descriptor_reference::provenance(task))
            .or_else(|| draw_reference::provenance(task))
            .or_else(|| native_draw_reference::provenance(task)),
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReferenceLabel {
    pub test: String,
    pub corpus_type: CorpusType,
    pub case_id: String,
    pub parameters: serde_json::Value,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct LabeledRecord {
    label: ReferenceLabel,
    record: Record,
}

fn reference_label(task: &Task, input: &Input) -> Result<ReferenceLabel> {
    let (case_id, parameters) = match input {
        Input::Search(row) => (
            &row.case.id,
            serde_json::to_value(&row.profile).map_err(|e| e.to_string())?,
        ),
        Input::Uff(row) => (
            &row.case.id,
            serde_json::to_value(row.profile).map_err(|e| e.to_string())?,
        ),
        Input::Fingerprint(row) => (
            &row.case.id,
            serde_json::json!({ "operation": row.operation, "width": row.width }),
        ),
        Input::Molecular { case, profile } => (
            &case.id,
            serde_json::to_value(profile).map_err(|e| e.to_string())?,
        ),
        Input::BioPdbOutput { case, profile } => (
            &case.id,
            serde_json::to_value(profile).map_err(|e| e.to_string())?,
        ),
    };
    Ok(ReferenceLabel {
        test: task.key(),
        corpus_type: task.corpus_type,
        case_id: case_id.clone(),
        parameters,
    })
}

fn label_records(task: &Task, records: &[Record]) -> Result<Vec<LabeledRecord>> {
    records
        .iter()
        .map(|record| {
            Ok(LabeledRecord {
                label: reference_label(task, &record.input)?,
                record: record.clone(),
            })
        })
        .collect()
}

fn check_records(task: &Task, inputs: &[Input], records: &[Record]) -> Result<()> {
    if inputs.len() != records.len() {
        return Err("reference row count mismatch".into());
    }
    for (input, record) in inputs.iter().zip(records) {
        task.validate_reference(input, &record.input, &record.output)?;
    }
    Ok(())
}

fn oracle(task: &Task, cases: &Corpus, python: &Path, threads: usize) -> Result<Vec<Record>> {
    if threads == 0 {
        return Err("threads must be positive".into());
    }
    if bio_reference::handles(task) {
        return bio_reference::generate(task, cases, python);
    }
    if tautomer_reference::handles(task) {
        return tautomer_reference::generate(task, cases);
    }
    if descriptor_reference::handles(task) {
        return descriptor_reference::generate(task, cases);
    }
    if draw_reference::handles(task) {
        return draw_reference::generate(task, cases);
    }
    if native_draw_reference::handles(task) {
        return native_draw_reference::generate(task, cases);
    }
    let script = root().join("tools/oracles/rdkit/fingerprint_values_pilot.py");
    // Run the exact checksummed source, from a real file so process workers
    // can import their module on spawn-based platforms too.
    if read(&script)? != ORACLE.as_bytes() {
        return Err("oracle source changed; rebuild the preparation binary".into());
    }
    let (corpus, parameters) = match task.operation {
        registry::Operation::SubstructureMatch => (
            serde_json::to_value(&cases.molecules),
            serde_json::to_value(search::profiles()),
        ),
        registry::Operation::Molecular(id) => (
            serde_json::to_value(&cases.molecules),
            serde_json::to_value(id.profiles()),
        ),
        operation @ (registry::Operation::UffCoverage
        | registry::Operation::UffOptimization
        | registry::Operation::UffConformerOptimization) => (
            serde_json::to_value(&cases.molecules),
            serde_json::to_value(uff::profiles(operation)),
        ),
        operation => (
            serde_json::to_value(&cases.fingerprints),
            Ok(serde_json::Value::Array(
                registry::fingerprint::WIDTHS
                    .iter()
                    .map(|width| serde_json::json!({"operation":operation,"width":width}))
                    .collect(),
            )),
        ),
    };
    let request = serde_json::json!({
        "generator":task.generator, "corpus":corpus.map_err(|e|e.to_string())?,
        "parameters":parameters.map_err(|e:serde_json::Error|e.to_string())?, "threads":threads,
    });
    // Python is solely the pinned RDKit adapter, never the CK executor.
    let mut child = Command::new(python)
        .arg(script)
        .arg(RDKIT_VERSION)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| format!("RDKit adapter: {e}"))?;
    let payload = encode(&request)?;
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
/// Framework regressions inject synthetic values; corpus tests never generate.
pub fn prepare(
    tasks: &[&Task],
    cases: &Corpus,
    data: &Path,
    python: &Path,
    threads: usize,
) -> Result<Preparation> {
    if threads == 0 {
        return Err("threads must be positive".into());
    }
    prepare_with(tasks, cases, data, |task, _| {
        oracle(task, cases, python, threads)
    })
    .map(|(_, preparation)| preparation)
}

fn prepare_with(
    tasks: &[&Task],
    cases: &Corpus,
    data: &Path,
    mut generate: impl FnMut(&Task, &[Input]) -> Result<Vec<Record>>,
) -> Result<(Ready, Preparation)> {
    registry::validate(cases, tasks)?;
    if tasks.is_empty() {
        return Err("empty task selection; 0 Rust operation calls".into());
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
            let records = generate(task, &inputs)
                .and_then(|records| {
                    check_records(task, &inputs, &records)?;
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
    let reference = encode(&label_records(task, records)?)?;
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
        "schema3{}{}{}{}{}",
        task.key(),
        digest(input),
        adapter_digest(task),
        registry_digest(),
        digest(reference_pin(task).as_bytes())
    );
    data.join(format!("{}-{}", task.key(), digest(identity.as_bytes())))
}

#[derive(Clone)]
pub struct Ready {
    records: Vec<LabeledRecord>,
}

impl Ready {
    pub fn len(&self) -> usize {
        self.records.len()
    }
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    pub fn compare_task(&self, key: &str) -> Vec<Comparison> {
        compare_rows(
            self.records.iter().filter(|row| row.label.test == key),
            execute::run,
        )
    }
}

/// Load and validate EVERY selected task before returning any executable work.
/// Owned snapshots prevent input changes after preflight affecting this run.
pub fn preflight(tasks: &[&Task], cases: &Corpus, data: &Path) -> Result<Ready> {
    registry::validate(cases, tasks)?;
    let mut all = Vec::new();
    let mut errors = Vec::new();
    for task in tasks {
        let load = || -> Result<Vec<LabeledRecord>> {
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
            let records: Vec<LabeledRecord> =
                serde_json::from_slice(&reference).map_err(|e| e.to_string())?;
            if records.len() != inputs.len() {
                return Err("reference row count mismatch".into());
            }
            for (input, row) in inputs.iter().zip(&records) {
                if row.label != reference_label(task, input)? {
                    return Err("reference label mismatch".into());
                }
                check_records(
                    task,
                    std::slice::from_ref(input),
                    std::slice::from_ref(&row.record),
                )?;
            }
            Ok(records)
        };
        match load() {
            Ok(records) => all.extend(records),
            Err(e) => errors.push(format!("{}: {e}", task.operation.name())),
        }
    }
    if !errors.is_empty() {
        return Err(format!(
            "preflight failed; 0 Rust operation calls\nPrepare references first: cargo run -p cosmolkit-parity-tests --release -- prepare\n{}",
            errors.join("\n")
        ));
    }
    Ok(Ready { records: all })
}

#[derive(Debug, Serialize)]
pub struct Comparison {
    pub label: ReferenceLabel,
    pub input: Input,
    pub expected: registry::Value,
    pub actual: std::result::Result<registry::Value, String>,
    pub matches: bool,
}

/// Exact typed comparison, except the registry's declared 2D numeric tolerance.
pub fn compare(ready: Ready, mut run: impl FnMut(&Input) -> Result<Record>) -> Vec<Comparison> {
    compare_rows(ready.records.iter(), &mut run)
}

fn compare_rows<'a>(
    rows: impl Iterator<Item = &'a LabeledRecord>,
    mut run: impl FnMut(&Input) -> Result<Record>,
) -> Vec<Comparison> {
    rows.map(|row| {
        let reference = &row.record;
        let actual = run(&reference.input).and_then(|record| {
            if record.input != reference.input {
                return Err("executor changed case identity".into());
            }
            Ok(record.output)
        });
        let matches = match (&reference.output, &actual) {
            (registry::Value::Uff(expected), Ok(registry::Value::Uff(actual))) => {
                uff::matches(&reference.input, expected, actual)
            }
            (registry::Value::Molecular(expected), Ok(registry::Value::Molecular(actual))) => {
                if matches!(
                    &reference.input,
                    Input::Molecular {
                        profile: registry::molecule_plan::Profile::SvgDefault,
                        ..
                    }
                ) {
                    molecular::svg_matches(expected, actual)
                } else {
                    molecular::matches(expected, actual)
                }
            }
            _ => actual.as_ref() == Ok(&reference.output),
        };
        Comparison {
            label: row.label.clone(),
            input: reference.input.clone(),
            expected: reference.output.clone(),
            actual,
            matches,
        }
    })
    .collect()
}

/// Read-only test stage. No interpreter or generator is accepted here.
pub fn run(
    tasks: &[&Task],
    cases: &Corpus,
    data: &Path,
    executor: impl FnMut(&Input) -> Result<Record>,
) -> Result<Vec<Comparison>> {
    Ok(compare(preflight(tasks, cases, data)?, executor))
}

pub fn write_report(path: &Path, report: &[Comparison]) -> Result<()> {
    fs::write(path, encode(&report)?).map_err(|e| e.to_string())
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Suite {
    pub tasks: Vec<String>,
    pub sources: Vec<CorpusSource>,
    #[serde(default)]
    pub special_regressions: Vec<String>,
}

/// Save the prepared selection, not a validity certificate: tests revalidate it.
pub fn save_suite(data: &Path, tasks: &[&Task], sources: &[CorpusSource]) -> Result<()> {
    save_suite_with_special(data, tasks, sources, Vec::new())
}

/// Persist both categories for the default complete preparation selection.
pub fn save_suite_with_special(
    data: &Path,
    tasks: &[&Task],
    sources: &[CorpusSource],
    special_regressions: Vec<String>,
) -> Result<()> {
    let sources = sources
        .iter()
        .map(|source| {
            Ok(CorpusSource {
                corpus_type: source.corpus_type,
                path: fs::canonicalize(&source.path).map_err(|e| e.to_string())?,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    let suite = Suite {
        tasks: tasks.iter().map(|task| task.key()).collect(),
        sources,
        special_regressions,
    };
    let mut file = tempfile::NamedTempFile::new_in(data).map_err(|e| e.to_string())?;
    file.write_all(&encode(&suite)?)
        .map_err(|e| e.to_string())?;
    file.persist(data.join("suite.json"))
        .map_err(|e| e.to_string())?;
    Ok(())
}

pub fn load_suite(data: &Path) -> Result<(Vec<&'static Task>, Ready)> {
    let suite: Suite =
        serde_json::from_slice(&read(&data.join("suite.json"))?).map_err(|e| e.to_string())?;
    special_regression::preflight_selection(&suite.special_regressions, data)?;
    let mut keys = std::collections::BTreeSet::new();
    let tasks = suite
        .tasks
        .iter()
        .map(|key| {
            if !keys.insert(key) {
                return Err("duplicate suite test".into());
            }
            let task = registry::select(Some(key))?.remove(0);
            if task.key() != *key {
                return Err("suite requires function/corpus-type keys".into());
            }
            Ok(task)
        })
        .collect::<Result<Vec<_>>>()?;
    let cases = corpus(&suite.sources, &tasks)?;
    let ready = preflight(&tasks, &cases, data)?;
    Ok((tasks, ready))
}

#[cfg(test)]
mod tests;
#[cfg(test)]
mod uff_tests;
