use crate::{
    Corpus, Input, Record, Result, Task, digest, directory, encode, expected, read, reference,
    registry, special_regression,
};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use sha2::Digest;
mod cache_import;
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write},
    path::{Path, PathBuf},
    sync::Mutex,
};

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub enum Selection {
    Corpus(String),
    Special(String),
}
impl Selection {
    pub fn folder(&self) -> Result<PathBuf> {
        let (kind, name) = match self {
            Self::Corpus(n) => ("corpus", n),
            Self::Special(n) => ("special", n),
        };
        if name.is_empty() || !name.chars().all(|c| c.is_ascii_alphanumeric() || c == '_') {
            return Err("selection names use ASCII letters, digits and underscores".into());
        }
        Ok(expected().join(kind).join(name))
    }
    pub fn command(&self) -> String {
        match self {
            Self::Corpus(n) => format!(
                "cargo run -p cosmolkit-parity-tests-fixed --release -- prepare --corpus {n}"
            ),
            Self::Special(_) => {
                "cargo run -p cosmolkit-parity-tests-fixed --release -- prepare --special all"
                    .into()
            }
        }
    }
}

#[derive(Clone, Copy)]
pub(crate) enum Spec {
    Corpus(&'static Task),
    Batch,
    Special(&'static registry::SpecialRegression),
}
impl Spec {
    pub fn key(self) -> &'static str {
        match self {
            Self::Corpus(task) => task.key,
            Self::Batch => "batch_smiles",
            Self::Special(s) => s.key,
        }
    }
}
pub(crate) struct Plan {
    pub selection: Selection,
    pub specs: Vec<Spec>,
    pub cases: Corpus,
}
pub(crate) struct Snapshot {
    pub spec: Spec,
    pub inputs: Value,
    pub rows: Vec<Value>,
}

/// A checked reference owns private bytes. Reopening an expected path after
/// preflight would permit a later preparation to change the checked values.
pub(crate) struct OwnedRows {
    file: Mutex<fs::File>,
    count: usize,
    sha256: String,
}
impl OwnedRows {
    fn copy_from(path: &Path) -> Result<Self> {
        let mut source = fs::File::open(path).map_err(|e| format!("{}: {e}", path.display()))?;
        // An anonymous temporary file is reclaimed by the OS even though
        // OnceLock/static values are not dropped at process exit.
        // References can exceed a host's RAM-backed /tmp quota. Keep private
        // snapshots on the package artifact filesystem, not in the source
        // reference directory (which may be read-only).
        let artifacts = directory().join("reports");
        fs::create_dir_all(&artifacts).map_err(|e| e.to_string())?;
        let mut file = tempfile::tempfile_in(artifacts).map_err(|e| e.to_string())?;
        let mut writer = BufWriter::new(&mut file);
        let mut hash = sha2::Sha256::new();
        std::io::copy(
            &mut source,
            &mut HashingWriter {
                writer: &mut writer,
                hash: &mut hash,
            },
        )
        .map_err(|e| e.to_string())?;
        writer.flush().map_err(|e| e.to_string())?;
        drop(writer);
        let mut rows = Self {
            file: Mutex::new(file),
            count: 0,
            sha256: hex_digest(hash),
        };
        // Preserve the original whole-file UTF-8 check before any JSON error.
        let mut reader = BufReader::new(rows.reader());
        let mut line = Vec::new();
        let mut offset = 0;
        loop {
            line.clear();
            let count = reader
                .read_until(b'\n', &mut line)
                .map_err(|e| e.to_string())?;
            if count == 0 {
                break;
            }
            if let Err(error) = std::str::from_utf8(&line) {
                let index = offset + error.valid_up_to();
                return Err(match error.error_len() {
                    Some(length) => {
                        format!("invalid utf-8 sequence of {length} bytes from index {index}")
                    }
                    None => format!("incomplete utf-8 byte sequence from index {index}"),
                });
            }
            offset += count;
        }
        drop(line);
        drop(reader);
        // Syntax of every line precedes manifest and typed/schema validation.
        let mut count = 0;
        for row in rows.iter()? {
            row?;
            count += 1;
        }
        rows.count = count;
        Ok(rows)
    }

    pub(crate) fn len(&self) -> usize {
        self.count
    }

    fn reader(&self) -> OwnedCursor<'_> {
        OwnedCursor {
            file: &self.file,
            offset: 0,
        }
    }

    pub(crate) fn iter(&self) -> Result<impl Iterator<Item = Result<Value>>> {
        Ok(BufReader::new(self.reader()).lines().map(|line| {
            let line = line.map_err(|e| e.to_string())?;
            serde_json::from_str(&line).map_err(|e| e.to_string())
        }))
    }
}

/// Each iterator has its own offset. A short lock around seek/read makes the
/// anonymous owned file portable without sharing the cursor of File::try_clone.
struct OwnedCursor<'a> {
    file: &'a Mutex<fs::File>,
    offset: u64,
}
impl Read for OwnedCursor<'_> {
    fn read(&mut self, bytes: &mut [u8]) -> std::io::Result<usize> {
        let mut file = self
            .file
            .lock()
            .map_err(|_| std::io::Error::other("owned reference file lock poisoned"))?;
        file.seek(SeekFrom::Start(self.offset))?;
        let count = file.read(bytes)?;
        self.offset += count as u64;
        Ok(count)
    }
}

struct GeneratedRows {
    file: tempfile::NamedTempFile,
    count: usize,
    sha256: String,
}
impl GeneratedRows {
    fn len(&self) -> usize {
        self.count
    }
    fn iter(&self) -> Result<impl Iterator<Item = Result<Value>>> {
        let reader = BufReader::new(self.file.reopen().map_err(|e| e.to_string())?);
        Ok(reader.lines().map(|line| {
            let line = line.map_err(|e| e.to_string())?;
            serde_json::from_str(&line).map_err(|e| e.to_string())
        }))
    }
    fn persist(self, path: &Path) -> Result<()> {
        self.file.persist(path).map_err(|e| e.to_string())?;
        Ok(())
    }
}

struct HashingWriter<'a, W> {
    writer: &'a mut W,
    hash: &'a mut sha2::Sha256,
}
impl<W: Write> Write for HashingWriter<'_, W> {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        let count = self.writer.write(bytes)?;
        self.hash.update(&bytes[..count]);
        Ok(count)
    }
    fn flush(&mut self) -> std::io::Result<()> {
        self.writer.flush()
    }
}
fn hex_digest(hash: sha2::Sha256) -> String {
    hash.finalize().iter().map(|b| format!("{b:02x}")).collect()
}
fn encoded_digest<T: Serialize + ?Sized>(value: &T) -> Result<String> {
    let mut hash = sha2::Sha256::new();
    serde_json::to_writer(
        HashingWriter {
            writer: &mut std::io::sink(),
            hash: &mut hash,
        },
        value,
    )
    .map_err(|e| e.to_string())?;
    Ok(hex_digest(hash))
}

/// Serializes one native array element at a time using the original compact
/// JSONL bytes. No size cutoff or chemistry-specific representation is added.
struct RowWriter {
    file: tempfile::NamedTempFile,
    writer: BufWriter<fs::File>,
    hash: sha2::Sha256,
    count: usize,
}
impl RowWriter {
    fn new_in(folder: &Path) -> Result<Self> {
        let file = tempfile::NamedTempFile::new_in(folder).map_err(|e| e.to_string())?;
        let writer = BufWriter::new(file.reopen().map_err(|e| e.to_string())?);
        Ok(Self {
            file,
            writer,
            hash: sha2::Sha256::new(),
            count: 0,
        })
    }
    fn push(&mut self, row: Value) -> Result<()> {
        let mut output = HashingWriter {
            writer: &mut self.writer,
            hash: &mut self.hash,
        };
        serde_json::to_writer(&mut output, &row).map_err(|e| e.to_string())?;
        output.write_all(b"\n").map_err(|e| e.to_string())?;
        self.count += 1;
        Ok(())
    }
    fn finish(mut self) -> Result<GeneratedRows> {
        self.writer.flush().map_err(|e| e.to_string())?;
        drop(self.writer);
        Ok(GeneratedRows {
            file: self.file,
            count: self.count,
            sha256: hex_digest(self.hash),
        })
    }
}

pub(crate) struct CheckedSnapshot {
    pub spec: Spec,
    // Corpus recipes are fully checked and then released. Only designated
    // special regressions need the complete owned fixture for their consumers.
    inputs: Option<Value>,
    pub rows: OwnedRows,
}

#[derive(Debug, Default)]
pub struct Preparation {
    pub generated_tasks: usize,
    pub reused_tasks: usize,
    pub rows: usize,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct CorpusManifest {
    format: String,
    #[serde(default)]
    input: Option<String>,
    #[serde(default)]
    records: Option<usize>,
    #[serde(default)]
    inputs: Vec<BioSource>,
    #[serde(default)]
    tautomer_profiles: TautomerProfiles,
    #[serde(default)]
    fingerprint_profiles: FingerprintProfiles,
}

#[derive(Default, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
enum FingerprintProfiles {
    #[default]
    Seeded,
    Historical0_3,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct NamedFingerprintProfile {
    name: String,
    params: crate::fingerprints::Params,
}

fn historical_fingerprint_profiles()
-> Result<BTreeMap<crate::fingerprints::Kind, Vec<crate::fingerprints::Params>>> {
    use crate::fingerprints::Kind;
    let fixture: Value = serde_json::from_slice(&read(
        &directory().join("testdata/profiles/fingerprints_historical_0_3.json"),
    )?)
    .map_err(|e| e.to_string())?;
    let recipes: BTreeMap<Kind, Vec<NamedFingerprintProfile>> =
        serde_json::from_value(fixture["profiles"].clone()).map_err(|e| e.to_string())?;
    if recipes.keys().copied().collect::<BTreeSet<_>>()
        != BTreeSet::from([
            Kind::Avalon,
            Kind::Topological,
            Kind::Layered,
            Kind::Pattern,
        ])
    {
        return Err(
            "historical fingerprint selection requires Avalon, Topological, Layered and Pattern recipes".into(),
        );
    }
    let mut selected = BTreeMap::new();
    for (kind, rows) in recipes {
        if rows.is_empty() {
            return Err("empty historical fingerprint profiles".into());
        }
        let mut names = BTreeSet::new();
        let mut profiles = Vec::with_capacity(rows.len());
        for row in rows {
            if row.name.is_empty()
                || !names.insert(row.name)
                || row.params.kind() != kind
                || profiles.contains(&row.params)
            {
                return Err("empty/duplicate/mismatched historical fingerprint profile".into());
            }
            profiles.push(row.params);
        }
        selected.insert(kind, profiles);
    }
    for (kind, file) in [
        (Kind::AtomPair, "atom_pair_historical_0_3.json"),
        (Kind::TopologicalTorsion, "torsion_historical_0_3.json"),
    ] {
        let rows: Vec<NamedFingerprintProfile> = serde_json::from_value(
            serde_json::from_slice::<Value>(&read(
                &directory().join("testdata/profiles").join(file),
            )?)
            .map_err(|e| e.to_string())?["profiles"]
                .clone(),
        )
        .map_err(|e| e.to_string())?;
        let mut names = BTreeSet::new();
        let mut profiles = Vec::new();
        for row in rows {
            if row.name.is_empty()
                || !names.insert(row.name)
                || row.params.kind() != kind
                || profiles.contains(&row.params)
            {
                return Err(format!(
                    "empty/duplicate/mismatched {kind:?} historical profile"
                ));
            }
            profiles.push(row.params);
        }
        if profiles.is_empty() {
            return Err(format!("empty {kind:?} historical profiles"));
        }
        selected.insert(kind, profiles);
    }
    Ok(selected)
}

#[derive(Default, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
enum TautomerProfiles {
    #[default]
    CurrentV1,
    Historical0_3,
}

fn historical_tautomer_profiles() -> Result<Vec<registry::molecule_plan::TautomerProfile>> {
    let fixture: Value = serde_json::from_slice(&read(
        &directory().join("testdata/special/tautomer_focused.json"),
    )?)
    .map_err(|e| e.to_string())?;
    let branches = fixture["branches"]
        .as_array()
        .ok_or("missing tautomer branches")?;
    if branches.is_empty() {
        return Err("empty tautomer profile selection".into());
    }
    let mut names = BTreeSet::new();
    let mut profiles = Vec::with_capacity(branches.len());
    for branch in branches {
        let mut fields = branch.as_object().ok_or("invalid tautomer branch")?.clone();
        let name = fields
            .remove("name")
            .ok_or("missing tautomer branch name")?;
        let name = name.as_str().ok_or("invalid tautomer branch name")?;
        if name.is_empty() || !names.insert(name.to_owned()) {
            return Err("empty/duplicate tautomer branch name".into());
        }
        let profile = serde_json::from_value(Value::Object(fields)).map_err(|e| e.to_string())?;
        if profiles.contains(&profile) {
            return Err("duplicate tautomer parameters".into());
        }
        profiles.push(profile);
    }
    Ok(profiles)
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct BioSource {
    id: String,
    format: registry::BioPdbCorpusFormat,
    input: String,
}

fn input_path(name: &str) -> Result<PathBuf> {
    let path = fs::canonicalize(directory().join("testdata/corpora").join(name))
        .map_err(|e| e.to_string())?;
    let base = fs::canonicalize(directory().join("testdata")).map_err(|e| e.to_string())?;
    if !path.starts_with(base) {
        return Err("corpus inputs must live inside parity-tests_fixed/testdata".into());
    }
    Ok(path)
}

pub(crate) fn plan(selection: &Selection) -> Result<Plan> {
    selection.folder()?;
    match selection {
        Selection::Special(name) => {
            let specs: Vec<_> = registry::SPECIAL_REGRESSIONS
                .iter()
                .filter(|s| name == "all" || s.key == name)
                .map(Spec::Special)
                .collect();
            if specs.is_empty() {
                return Err(format!("unknown special regression: {name}"));
            }
            Ok(Plan {
                selection: selection.clone(),
                specs,
                cases: Corpus::default(),
            })
        }
        Selection::Corpus(name) => {
            let manifest: CorpusManifest =
                serde_json::from_slice(&read(&input_path(&format!("{name}.json"))?)?)
                    .map_err(|e| e.to_string())?;
            let mut cases = Corpus::default();
            if manifest.fingerprint_profiles == FingerprintProfiles::Historical0_3 {
                if manifest.format != "smiles" {
                    return Err("fingerprint profiles require a SMILES corpus".into());
                }
                cases.fingerprint_profiles = historical_fingerprint_profiles()?;
                let fixture: Value = serde_json::from_slice(&read(
                    &directory().join("testdata/profiles/morgan_historical_0_3.json"),
                )?)
                .map_err(|e| e.to_string())?;
                let recipes: Vec<crate::morgan_historical::Recipe> =
                    serde_json::from_value(fixture["recipes"].clone())
                        .map_err(|e| e.to_string())?;
                if recipes != crate::morgan_historical::Recipe::ALL {
                    return Err(
                        "historical Morgan recipes differ from the complete source scope".into(),
                    );
                }
                cases.morgan_historical = Some(recipes);
            }
            if manifest.tautomer_profiles == TautomerProfiles::Historical0_3 {
                if manifest.format != "smiles" {
                    return Err("tautomer profiles require a SMILES corpus".into());
                }
                cases.tautomer_profiles = Some(historical_tautomer_profiles()?);
            }
            let kinds = match manifest.format.as_str() {
                "smiles" => {
                    cases.molecules = crate::molecular::read_corpus(&input_path(
                        manifest
                            .input
                            .as_deref()
                            .ok_or("SMILES corpus needs input")?,
                    )?)?;
                    vec![registry::CorpusType::Smiles]
                }
                "bio" => {
                    for source in manifest.inputs {
                        cases.bio_cases.push(registry::BioPdbCase {
                            id: source.id,
                            format: source.format,
                            text: String::from_utf8(read(&input_path(&source.input)?)?)
                                .map_err(|e| e.to_string())?,
                        });
                    }
                    vec![registry::CorpusType::Pdb, registry::CorpusType::Cif]
                }
                _ => return Err(format!("unsupported corpus format: {}", manifest.format)),
            };
            let count = cases.molecules.len() + cases.bio_cases.len();
            if count == 0 || manifest.records.is_some_and(|n| n != count) {
                return Err(format!("corpus {name}: unexpected record count {count}"));
            }
            let mut specs: Vec<_> = registry::TASKS
                .iter()
                .filter(|t| kinds.contains(&t.corpus_type))
                .map(Spec::Corpus)
                .collect();
            if !cases.molecules.is_empty() {
                specs.push(Spec::Batch);
            }
            let tasks: Vec<_> = specs
                .iter()
                .filter_map(|s| {
                    if let Spec::Corpus(t) = *s {
                        Some(t)
                    } else {
                        None
                    }
                })
                .collect();
            registry::validate(&cases, &tasks)?;
            let keys: BTreeSet<_> = specs.iter().map(|s| s.key()).collect();
            if keys.len() != specs.len() || specs.is_empty() {
                return Err("empty/duplicate registry selection".into());
            }
            Ok(Plan {
                selection: selection.clone(),
                specs,
                cases,
            })
        }
    }
}

pub(crate) fn plan_task(selection: &Selection, task: Option<&str>) -> Result<Plan> {
    let mut selected = plan(selection)?;
    if let Some(key) = task {
        if !matches!(selection, Selection::Corpus(_)) {
            return Err(
                "--task applies to corpora; select special regressions with --special NAME".into(),
            );
        }
        selected.specs.retain(|spec| spec.key() == key);
        if selected.specs.is_empty() {
            return Err(format!(
                "unknown or inapplicable task {key:?} for {selection:?}"
            ));
        }
    }
    Ok(selected)
}

pub(crate) fn inputs(spec: &Spec, cases: &Corpus) -> Result<Value> {
    match *spec {
        Spec::Corpus(task) => {
            serde_json::to_value(registry::expand(cases, task)).map_err(|e| e.to_string())
        }
        Spec::Batch => Ok(
            json!([{"cases": cases.molecules, "workers": 1}, {"cases": cases.molecules, "workers": 4}]),
        ),
        Spec::Special(s) => {
            let mut fixture: Value =
                serde_json::from_slice(&read(&directory().join("testdata").join(s.fixture))?)
                    .map_err(|e| e.to_string())?;
            if matches!(
                s.schema,
                registry::SpecialRegressionSchema::BioMmcifSwitches
            ) {
                fixture["flags"] = json!(crate::bio_mmcif::FLAGS);
                for case in fixture["cases"].as_array_mut().ok_or("missing BIO cases")? {
                    let source = case["input"].as_str().ok_or("missing BIO input")?;
                    case["text"] = json!(
                        String::from_utf8(read(&input_path(&format!("../{source}"))?)?)
                            .map_err(|e| e.to_string())?
                    );
                }
            }
            Ok(fixture)
        }
    }
}

pub(crate) fn jsonl(rows: &[Value]) -> Result<Vec<u8>> {
    let mut output = Vec::new();
    for row in rows {
        output.extend(encode(row)?);
        output.push(b'\n');
    }
    Ok(output)
}

pub(crate) fn validate_rows<I>(spec: &Spec, input: &Value, rows: I, row_count: usize) -> Result<()>
where
    I: IntoIterator<Item = Result<Value>>,
{
    let mut rows = rows.into_iter();
    match *spec {
        Spec::Corpus(task) => {
            let recipes = input
                .as_array()
                .ok_or("reference recipes must be an array")?;
            // The former Vec<Input> conversion checked all recipe types before
            // checking row counts or typed reference rows. Retain that order
            // without retaining a second complete typed recipe collection.
            for recipe in recipes {
                Input::deserialize(recipe).map_err(|e| e.to_string())?;
            }
            if recipes.is_empty() || recipes.len() != row_count {
                return Err("reference row count mismatch".into());
            }
            for recipe in recipes {
                let recipe = Input::deserialize(recipe).map_err(|e| e.to_string())?;
                let row = rows.next().ok_or("reference row count mismatch")??;
                let record: Record = serde_json::from_value(row).map_err(|e| e.to_string())?;
                task.validate_reference(&recipe, &record.input, &record.output)?;
            }
        }
        Spec::Batch => {
            let recipes = input.as_array().ok_or("batch recipes must be an array")?;
            if recipes.is_empty() || recipes.len() != row_count {
                return Err("batch row count mismatch".into());
            }
            for recipe in recipes {
                let row = rows.next().ok_or("batch row count mismatch")??;
                let count = recipe["cases"]
                    .as_array()
                    .ok_or("batch cases missing")?
                    .len();
                let mask: Vec<bool> =
                    Vec::deserialize(&row["output"]["valid_mask"]).map_err(|e| e.to_string())?;
                let smiles: Vec<Option<String>> =
                    Vec::deserialize(&row["output"]["smiles"]).map_err(|e| e.to_string())?;
                let errors: Vec<usize> =
                    Vec::deserialize(&row["output"]["error_indices"]).map_err(|e| e.to_string())?;
                let invalid: Vec<_> = mask
                    .iter()
                    .enumerate()
                    .filter_map(|(i, valid)| (!valid).then_some(i))
                    .collect();
                if row["input"] != *recipe
                    || mask.len() != count
                    || smiles.len() != count
                    || errors != invalid
                    || mask
                        .iter()
                        .zip(smiles)
                        .any(|(valid, text)| *valid != text.is_some())
                {
                    return Err("batch identity/order/error schema mismatch".into());
                }
            }
        }
        Spec::Special(s) => {
            // The existing public special-regression API requires complete
            // owned fixture/rows. Preserve its complete schema validation.
            let rows: Vec<Value> = rows.collect::<Result<_>>()?;
            return special_regression::validate(&encode(input)?, &jsonl(&rows)?, s.rows, s.schema)
                .map(|_| ());
        }
    }
    if rows.next().transpose()?.is_some() {
        return Err("reference row count mismatch".into());
    }
    Ok(())
}

#[derive(Serialize, Deserialize, PartialEq, Eq, Debug)]
#[serde(deny_unknown_fields)]
struct Manifest {
    schema: u32,
    selection: Selection,
    task: String,
    input_sha256: String,
    generator_sha256: String,
    reference_identity: Value,
    platform: String,
    output_sha256: String,
    rows: usize,
}
fn identity(
    selection: &Selection,
    spec: &Spec,
    input: &Value,
    output: &[u8],
    rows: usize,
) -> Result<Manifest> {
    identity_with_output_digest(selection, spec, input, digest(output), rows)
}

fn identity_with_output_digest(
    selection: &Selection,
    spec: &Spec,
    input: &Value,
    output_sha256: String,
    rows: usize,
) -> Result<Manifest> {
    let pin = if reference::uses_gemmi(spec) {
        "gemmi.json"
    } else {
        "rdkit.json"
    };
    Ok(Manifest {
        schema: 1,
        selection: selection.clone(),
        task: spec.key().into(),
        input_sha256: encoded_digest(input)?,
        generator_sha256: reference::source_digest_for_inputs(spec, input)?,
        reference_identity: serde_json::from_slice(&read(
            &directory().join("testdata/reference").join(pin),
        )?)
        .map_err(|e| e.to_string())?,
        platform: format!("{}-{}", std::env::consts::ARCH, std::env::consts::OS),
        output_sha256,
        rows,
    })
}

fn load_reference(plan: &Plan, spec: &Spec, dir: &Path) -> Result<CheckedSnapshot> {
    let input = inputs(spec, &plan.cases)?;
    let stored: Value =
        serde_json::from_slice(&read(&dir.join("input.json"))?).map_err(|e| e.to_string())?;
    if stored != input {
        return Err("prepared input changed".into());
    }
    drop(stored);
    let rows = OwnedRows::copy_from(&dir.join("reference.jsonl"))?;
    let manifest: Manifest =
        serde_json::from_slice(&read(&dir.join("manifest.json"))?).map_err(|e| e.to_string())?;
    if manifest
        != identity_with_output_digest(
            &plan.selection,
            spec,
            &input,
            rows.sha256.clone(),
            rows.len(),
        )?
    {
        return Err("stale/corrupt reference manifest".into());
    }
    validate_rows(spec, &input, rows.iter()?, rows.len())?;
    Ok(CheckedSnapshot {
        spec: *spec,
        inputs: matches!(spec, Spec::Special(_)).then_some(input),
        rows,
    })
}

/// Finish every selected preflight before returning any references to CK.
/// Completion order must never select the reported failure: results are reduced
/// in the original task order after all bounded workers have joined.
fn parallel_preflight<T, F>(count: usize, workers: usize, check: F) -> Result<Vec<T>>
where
    T: Send,
    F: Fn(usize) -> Result<T> + Sync,
{
    if count == 0 {
        return Err("empty test selection".into());
    }
    if workers == 0 {
        return Err("reference preflight workers must be positive".into());
    }
    let workers = workers.min(112).min(count);
    let next = std::sync::atomic::AtomicUsize::new(0);
    let mut ordered: Vec<Option<Result<T>>> = (0..count).map(|_| None).collect();
    std::thread::scope(|scope| -> Result<()> {
        let mut handles = Vec::with_capacity(workers);
        for _ in 0..workers {
            let next = &next;
            let check = &check;
            handles.push(scope.spawn(move || {
                let mut completed = Vec::new();
                loop {
                    let index = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    if index >= count {
                        break;
                    }
                    completed.push((index, check(index)));
                }
                completed
            }));
        }
        let mut panicked = false;
        for handle in handles {
            match handle.join() {
                Ok(completed) => {
                    for (index, result) in completed {
                        ordered[index] = Some(result);
                    }
                }
                Err(_) => panicked = true,
            }
        }
        if panicked {
            return Err("reference preflight worker panicked".into());
        }
        Ok(())
    })?;
    ordered
        .into_iter()
        .map(|result| {
            result.ok_or_else(|| "reference preflight worker omitted a task".to_string())?
        })
        .collect()
}

pub(crate) fn load_references(
    plan: &Plan,
    folder: &Path,
) -> Result<BTreeMap<&'static str, CheckedSnapshot>> {
    let workers = std::thread::available_parallelism()
        .map(usize::from)
        .unwrap_or(1)
        .min(112);
    parallel_preflight(plan.specs.len(), workers, |index| {
        let s = &plan.specs[index];
        load_reference(plan, s, &folder.join(s.key()))
            .map(|snapshot| (s.key(), snapshot))
            .map_err(|e| {
                let mut command = plan.selection.command();
                if matches!(plan.selection, Selection::Corpus(_)) {
                    command.push_str(&format!(" --task {}", s.key()));
                }
                format!("{}: {e}\nPrepare first: {command}", s.key())
            })
    })
    .map(|checked| checked.into_iter().collect())
}

pub(crate) fn load_special_references(
    plan: &Plan,
    folder: &Path,
) -> Result<BTreeMap<&'static str, Snapshot>> {
    // Complete every selected special preflight before materializing consumers.
    load_references(plan, folder)?
        .into_iter()
        .map(|(key, checked)| Ok((key, materialize_special(checked)?)))
        .collect()
}

fn materialize_special(checked: CheckedSnapshot) -> Result<Snapshot> {
    Ok(Snapshot {
        spec: checked.spec,
        inputs: checked.inputs.ok_or("special fixture unavailable")?,
        rows: checked.rows.iter()?.collect::<Result<_>>()?,
    })
}

pub fn prepare(selection: &Selection, threads: usize, task: Option<&str>) -> Result<Preparation> {
    prepare_with_reuse(selection, threads, task, None)
}

pub fn prepare_with_reuse(
    selection: &Selection,
    threads: usize,
    task: Option<&str>,
    reuse_from: Option<&Path>,
) -> Result<Preparation> {
    if threads == 0 {
        return Err("threads must be positive".into());
    }
    let plan = plan_task(selection, task)?;
    let folder = selection.folder()?;
    fs::create_dir_all(&folder).map_err(|e| e.to_string())?;
    let lock = fs::OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .truncate(false)
        .open(folder.join(".prepare.lock"))
        .map_err(|e| e.to_string())?;
    lock.lock().map_err(|e| e.to_string())?;
    let mut result = Preparation::default();
    let mut descriptors = None;
    for (index, spec) in plan.specs.iter().enumerate() {
        eprintln!("Task {}/{}: {}", index + 1, plan.specs.len(), spec.key());
        if let Ok(snapshot) = load_reference(&plan, spec, &folder.join(spec.key())) {
            eprintln!(
                "  [========================] reused {} reference rows",
                snapshot.rows.len()
            );
            result.reused_tasks += 1;
            result.rows += snapshot.rows.len();
            continue;
        }
        let input = inputs(spec, &plan.cases)?;
        if let Some(source) = reuse_from {
            if let Some(rows) = cache_import::import(&plan, spec, &input, source, &folder)? {
                result.reused_tasks += 1;
                result.rows += rows;
                eprintln!("  [========================] imported {rows} verified native rows");
                continue;
            }
        }
        let generator_before = reference::source_digest_for_inputs(spec, &input)?;
        eprintln!("Generating {}", spec.key());
        let temporary = tempfile::Builder::new()
            .prefix(".prepare-")
            .tempdir_in(&folder)
            .map_err(|e| e.to_string())?;
        let mut writer = RowWriter::new_in(temporary.path())?;
        reference::generate(
            spec,
            &plan.cases,
            &input,
            threads,
            &mut descriptors,
            |row| writer.push(row),
        )?;
        let rows = writer.finish()?;
        if reference::source_digest_for_inputs(spec, &input)? != generator_before {
            return Err("reference generator changed during preparation".into());
        }
        validate_rows(spec, &input, rows.iter()?, rows.len())?;
        let manifest =
            identity_with_output_digest(selection, spec, &input, rows.sha256.clone(), rows.len())?;
        let row_count = rows.len();
        let mut input_file = BufWriter::new(
            fs::File::create(temporary.path().join("input.json")).map_err(|e| e.to_string())?,
        );
        serde_json::to_writer(&mut input_file, &input).map_err(|e| e.to_string())?;
        input_file.flush().map_err(|e| e.to_string())?;
        drop(input_file);
        rows.persist(&temporary.path().join("reference.jsonl"))?;
        fs::write(temporary.path().join("manifest.json"), encode(&manifest)?)
            .map_err(|e| e.to_string())?;
        let destination = folder.join(spec.key());
        let backup = folder.join(format!(
            ".invalid-{}-{}",
            spec.key(),
            temporary.path().file_name().unwrap().to_string_lossy()
        ));
        let exists = destination.exists();
        if exists {
            fs::rename(&destination, &backup).map_err(|e| e.to_string())?;
        }
        if let Err(error) = fs::rename(temporary.path(), &destination) {
            if exists {
                fs::rename(&backup, &destination)
                    .map_err(|e| format!("publish {error}; rollback {e}"))?;
            }
            return Err(error.to_string());
        }
        result.generated_tasks += 1;
        result.rows += row_count;
        eprintln!(
            "  [========================] saved {} reference rows",
            row_count
        );
    }
    // Validate all tasks again, dropping each owned temporary snapshot before
    // proceeding. Publishing selection still waits for every final preflight.
    for spec in &plan.specs {
        load_reference(&plan, spec, &folder.join(spec.key())).map_err(|e| {
            let mut command = plan.selection.command();
            if matches!(plan.selection, Selection::Corpus(_)) {
                command.push_str(&format!(" --task {}", spec.key()));
            }
            format!("{}: {e}\nPrepare first: {command}", spec.key())
        })?;
    }
    let lane = if matches!(selection, Selection::Corpus(_)) {
        "corpus"
    } else {
        "special"
    };
    let name = match selection {
        Selection::Corpus(n) | Selection::Special(n) => n,
    };
    let active = expected().join(lane).join("selection.json");
    let bytes = encode(name)?;
    if fs::read(&active).ok().as_deref() == Some(bytes.as_slice()) {
        return Ok(result);
    }
    let mut selected =
        tempfile::NamedTempFile::new_in(active.parent().unwrap()).map_err(|e| e.to_string())?;
    use std::io::Write;
    selected.write_all(&bytes).map_err(|e| e.to_string())?;
    selected.persist(active).map_err(|e| e.to_string())?;
    Ok(result)
}

#[cfg(test)]
mod tests {
    #[test]
    fn historical_morgan_expansion_retains_ordinary_profiles_and_exact_eligibility() {
        use crate::registry::molecule_plan::{Profile, TaskId};
        let mut cases = Corpus::default();
        cases.morgan_historical = Some(crate::morgan_historical::Recipe::ALL.to_vec());
        for (smiles, historical) in [("", 10), ("C", 14), ("CC", 15), ("C.C", 14), ("CC(", 15)] {
            let case = registry::SmilesCase {
                id: smiles.into(),
                smiles: smiles.into(),
            };
            let profiles = cases.molecular_profiles_for_case(TaskId::MorganFingerprint, &case);
            assert_eq!(&profiles[..16], TaskId::MorganFingerprint.profiles());
            assert_eq!(profiles.len(), 16 + historical);
            assert!(
                profiles[16..]
                    .iter()
                    .all(|p| matches!(p, Profile::MorganHistorical { .. }))
            );
            assert_eq!(
                cases.molecular_profiles_for_case(TaskId::MorganSparseFingerprint, &case),
                TaskId::MorganSparseFingerprint.profiles()
            );
            cases.molecules.push(case);
        }
        let task = registry::TASKS
            .iter()
            .find(|t| t.key == "fingerprint_morgan_smiles")
            .unwrap();
        let inputs = registry::expand(&cases, task);
        assert_eq!(task.count(&cases), inputs.len());
        assert_eq!(inputs.len(), 5 * 16 + 10 + 14 + 15 + 14 + 15);
    }
    use super::*;
    #[cfg(unix)]
    #[test]
    fn checked_snapshot_uses_artifact_filesystem_and_survives_source_removal() {
        use std::os::unix::fs::MetadataExt;
        let source = tempfile::NamedTempFile::new().unwrap();
        fs::write(source.path(), b"{\"output\":42}\n").unwrap();
        let snapshot = OwnedRows::copy_from(source.path()).unwrap();
        assert_eq!(
            snapshot.file.lock().unwrap().metadata().unwrap().dev(),
            fs::metadata(directory().join("reports")).unwrap().dev()
        );
        drop(source);
        assert_eq!(
            snapshot.iter().unwrap().next().unwrap().unwrap(),
            json!({"output":42})
        );
    }
    #[test]
    fn registry_keys_and_cargo_declarations_share_one_census() {
        assert_eq!(registry::TASKS.len(), 119);
        let mut unique = BTreeSet::new();
        for task in registry::TASKS {
            assert_eq!(
                task.key,
                format!("{}_{}", task.operation.name(), task.corpus_type.name())
            );
            assert!(unique.insert(task.key));
        }
    }
    #[test]
    fn named_corpora_load_only_owned_inputs_with_exact_counts() {
        let smiles = plan(&Selection::Corpus("smiles_5000".into())).unwrap();
        assert_eq!(smiles.cases.molecules.len(), 5000);
        assert!(smiles.cases.bio_cases.is_empty());
        assert!(smiles.specs.iter().any(|s| s.key() == "batch_smiles"));
        assert!(plan(&Selection::Corpus("fingerprint_5000".into())).is_err());
        let smoke = plan(&Selection::Corpus("smiles_smoke".into())).unwrap();
        assert_eq!(smoke.cases.molecules.len(), 3);
        assert_eq!(smoke.specs.len(), 118);
        let bio = plan(&Selection::Corpus("bio_small".into())).unwrap();
        assert_eq!(bio.cases.bio_cases.len(), 2);
        assert_eq!(bio.specs.len(), 2);
    }

    #[test]
    fn torsion_historical_recipes_preserve_every_form_and_argument() {
        use crate::fingerprints::{Kind, Params};
        use crate::torsion::{Form, Profile};
        let bytes = read(
            &crate::root()
                .join("tools/testdata/rdkit/topological_torsion_fingerprint_profile.json"),
        )
        .unwrap();
        let original: Value = serde_json::from_slice(&bytes).unwrap();
        let fixture: Value = serde_json::from_slice(
            &read(&directory().join("testdata/profiles/torsion_historical_0_3.json")).unwrap(),
        )
        .unwrap();
        assert_eq!(fixture["source_sha256"], digest(&bytes));
        let recipes = historical_fingerprint_profiles().unwrap();
        let profiles = &recipes[&Kind::TopologicalTorsion];
        assert_eq!(profiles.len(), 36);
        for (index, branch) in original["corpus_branches"]
            .as_array()
            .unwrap()
            .iter()
            .enumerate()
        {
            for (offset, form) in [Form::SparseCount, Form::Count, Form::SparseBits, Form::Bits]
                .into_iter()
                .enumerate()
            {
                let expected = Profile {
                    form,
                    torsion_atom_count: branch["torsionAtomCount"].as_u64().unwrap_or(4) as u32,
                    only_shortest_paths: branch["onlyShortestPaths"].as_bool().unwrap_or(false),
                    include_chirality: branch["includeChirality"].as_bool().unwrap_or(false),
                    count_simulation: branch["countSimulation"].as_bool().unwrap_or(true),
                    count_bounds: branch.get("countBounds").map_or(vec![1, 2, 4, 8], |v| {
                        serde_json::from_value(v.clone()).unwrap()
                    }),
                    fp_size: branch["fpSize"].as_u64().unwrap_or(2048) as u32,
                    bits_per_feature: branch["numBitsPerFeature"].as_u64().unwrap_or(1) as u32,
                    from_first: branch["fromAtoms"] == "first",
                    ignore_first: branch["ignoreAtoms"] == "first",
                    index_plus_17: branch["customAtomInvariants"] == "index_plus_17",
                    additional_output: branch["additionalOutput"] == true,
                };
                assert_eq!(
                    profiles[index * 4 + offset],
                    Params::TopologicalTorsion(expected)
                );
                assert_eq!(
                    fixture["profiles"][index * 4 + offset]["name"],
                    format!(
                        "{}.{}",
                        branch["name"].as_str().unwrap(),
                        ["SparseCount", "Count", "SparseBits", "Bits"][offset]
                    )
                );
            }
        }
        let mut corpus = plan(&Selection::Corpus("smiles_smoke".into()))
            .unwrap()
            .cases;
        let task = registry::TASKS
            .iter()
            .find(|t| t.operation == registry::Operation::Fingerprint(Kind::TopologicalTorsion))
            .unwrap();
        corpus
            .fingerprint_profiles
            .insert(Kind::TopologicalTorsion, profiles.clone());
        assert_eq!(task.count(&corpus), corpus.molecules.len() * 36);
        let expanded = registry::expand(&corpus, task);
        assert_eq!(expanded.len(), corpus.molecules.len() * 36);
        for (index, input) in expanded.iter().enumerate() {
            let Input::Fingerprint(row) = input else {
                panic!("wrong input family")
            };
            assert_eq!(row.case, corpus.molecules[index / 36]);
            assert_eq!(row.params, profiles[index % 36]);
        }
    }

    #[test]
    fn atom_pair_historical_recipes_preserve_every_form_and_argument() {
        use crate::atom_pair::{Form, Profile};
        use crate::fingerprints::{Kind, Params};
        let bytes =
            read(&crate::root().join("tools/testdata/rdkit/atom_pair_fingerprint_profile.json"))
                .unwrap();
        let original: Value = serde_json::from_slice(&bytes).unwrap();
        let fixture: Value = serde_json::from_slice(
            &read(&directory().join("testdata/profiles/atom_pair_historical_0_3.json")).unwrap(),
        )
        .unwrap();
        assert_eq!(fixture["source_sha256"], digest(&bytes));
        let recipes = historical_fingerprint_profiles().unwrap();
        let profiles = &recipes[&Kind::AtomPair];
        assert_eq!(profiles.len(), 40);
        for (index, branch) in original["branches"].as_array().unwrap().iter().enumerate() {
            for (offset, form) in [Form::SparseCount, Form::Count, Form::SparseBits, Form::Bits]
                .into_iter()
                .enumerate()
            {
                let expected = Profile {
                    form,
                    min_distance: branch["minDistance"].as_u64().unwrap() as u32,
                    max_distance: branch["maxDistance"].as_u64().unwrap() as u32,
                    include_chirality: branch["includeChirality"].as_bool().unwrap(),
                    use_2d: branch["use2D"].as_bool().unwrap(),
                    count_simulation: branch["countSimulation"].as_bool().unwrap(),
                    count_bounds: serde_json::from_value(branch["countBounds"].clone()).unwrap(),
                    fp_size: branch["fpSize"].as_u64().unwrap() as u32,
                    bits_per_feature: branch["numBitsPerFeature"].as_u64().unwrap() as u32,
                    from_first: branch["fromAtoms"] == "first",
                    ignore_first: branch["ignoreAtoms"] == "first",
                    index_plus_11: branch["customAtomInvariants"] == "index_plus_11",
                    additional_output: branch["additionalOutput"] == true,
                };
                assert_eq!(profiles[index * 4 + offset], Params::AtomPair(expected));
                assert_eq!(
                    fixture["profiles"][index * 4 + offset]["name"],
                    format!(
                        "{}.{}",
                        branch["name"].as_str().unwrap(),
                        ["SparseCount", "Count", "SparseBits", "Bits"][offset]
                    )
                );
            }
        }
        let mut corpus = plan(&Selection::Corpus("smiles_smoke".into()))
            .unwrap()
            .cases;
        let task = registry::TASKS
            .iter()
            .find(|t| t.operation == registry::Operation::Fingerprint(Kind::AtomPair))
            .unwrap();
        corpus
            .fingerprint_profiles
            .insert(Kind::AtomPair, profiles.clone());
        assert_eq!(task.count(&corpus), corpus.molecules.len() * 40);
        let expanded = registry::expand(&corpus, task);
        assert_eq!(expanded.len(), corpus.molecules.len() * 40);
        for (index, input) in expanded.iter().enumerate() {
            let Input::Fingerprint(row) = input else {
                panic!("wrong input family")
            };
            assert_eq!(row.case, corpus.molecules[index / 40]);
            assert_eq!(row.params, profiles[index % 40]);
        }
    }

    #[test]
    fn historical_fingerprint_selection_preserves_recipes_and_censuses() {
        use crate::fingerprints::Kind;
        let selected = plan(&Selection::Corpus("smiles_smoke".into())).unwrap();
        let mut cases = selected.cases;
        let fixture: Value = serde_json::from_slice(
            &read(&directory().join("testdata/profiles/fingerprints_historical_0_3.json")).unwrap(),
        )
        .unwrap();
        let profiles = historical_fingerprint_profiles().unwrap();
        for (kind, file, key, count) in [
            (
                Kind::Avalon,
                "avalon_fingerprint_profile.json",
                "corpus_branches",
                23,
            ),
            (
                Kind::Topological,
                "rdkit_topological_fingerprint_profile.json",
                "branches",
                14,
            ),
            (
                Kind::Layered,
                "layered_fingerprint_profile.json",
                "branches",
                18,
            ),
            (
                Kind::Pattern,
                "pattern_fingerprint_profile.json",
                "corpus_branch_names",
                10,
            ),
        ] {
            let bytes = read(&crate::root().join("tools/testdata/rdkit").join(file)).unwrap();
            let original: Value = serde_json::from_slice(&bytes).unwrap();
            let historical: Vec<&Value> = if kind == Kind::Pattern {
                original[key]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|name| {
                        original["branches"]
                            .as_array()
                            .unwrap()
                            .iter()
                            .find(|branch| branch["name"] == *name)
                            .unwrap()
                    })
                    .collect()
            } else {
                original[key].as_array().unwrap().iter().collect()
            };
            assert_eq!(historical.len(), count);
            assert_eq!(profiles[&kind].len(), count);
            let kind_name = match kind {
                Kind::Avalon => "Avalon",
                Kind::Topological => "Topological",
                Kind::Layered => "Layered",
                Kind::Pattern => "Pattern",
                _ => unreachable!(),
            };
            let provenance = fixture["provenance"]["sources"]
                .as_array()
                .unwrap()
                .iter()
                .find(|source| source["path"].as_str().unwrap().ends_with(file))
                .unwrap();
            assert_eq!(provenance["sha256"], digest(&bytes));
            for (index, branch) in historical.iter().enumerate() {
                assert_eq!(
                    fixture["profiles"][kind_name][index]["name"],
                    branch["name"]
                );
                let parameters = if kind == Kind::Avalon {
                    json!({"Avalon": {
                        "n_bits":branch["n_bits"], "is_query":branch["is_query"],
                        "bit_flags":u32::from_str_radix(branch["bit_flags"].as_str().unwrap().trim_start_matches("0x"),16).unwrap()
                    }})
                } else if kind == Kind::Topological {
                    json!({"Topological": {
                        "min_path":branch["minPath"], "max_path":branch["maxPath"],
                        "fp_size":branch["fpSize"], "bits_per_feature":branch["nBitsPerHash"],
                        "use_hs":branch["useHs"], "density_milli":(branch["tgtDensity"].as_f64().unwrap()*1000.0) as u32,
                        "min_size":branch["minSize"], "branched":branch["branchedPaths"],
                        "bond_order":branch["useBondOrder"],
                        "custom_invariants":branch["atomInvariants"]=="index_plus_one",
                        "roots":if branch["fromAtoms"]=="first" { "First" } else { "All" }
                    }})
                } else if kind == Kind::Layered {
                    json!({"Layered": {
                        "layers":u32::from_str_radix(branch["layerFlags"].as_str().unwrap().trim_start_matches("0x"),16).unwrap(),
                        "min_path":branch["minPath"], "max_path":branch["maxPath"],
                        "fp_size":branch["fpSize"], "branched":branch["branchedPaths"],
                        "roots":match branch["fromAtoms"].as_str() { Some("first") => "First", Some("terminal_pair") => "Terminals", None => "All", _=>panic!("unknown roots") },
                        "counts":match branch["atomCounts"].as_str() { Some("zeros") => "Zero", Some("index_plus_10") => "Seeded", None => "Absent", _=>panic!("unknown counts") },
                        "mask":match branch["setOnlyBits"].as_str() { Some("even") => "Even", Some("mod_three") => "ModThree", None => "Absent", _=>panic!("unknown mask") }
                    }})
                } else {
                    json!({"Pattern": {
                        "fp_size":branch["fpSize"], "tautomeric":branch["tautomericFingerprint"],
                        "source_options":{"atom_counts":branch["atomCounts"], "set_only_bits":branch["setOnlyBits"]}
                    }})
                };
                assert_eq!(
                    serde_json::to_value(&profiles[&kind][index]).unwrap(),
                    parameters
                );
            }
            let task = registry::TASKS
                .iter()
                .find(|task| task.operation == registry::Operation::Fingerprint(kind))
                .unwrap();
            let normal = registry::expand(&cases, task);
            assert_eq!(normal.len(), 3);
            assert_eq!(task.count(&cases), 3);
            cases
                .fingerprint_profiles
                .insert(kind, profiles[&kind].clone());
            let expanded = registry::expand(&cases, task);
            assert_eq!(task.count(&cases), 3 * count);
            assert_eq!(expanded.len(), 3 * count);
            for (index, input) in expanded.iter().enumerate() {
                let Input::Fingerprint(row) = input else {
                    panic!("wrong input family")
                };
                let Input::Fingerprint(sampled) = &normal[index / count] else {
                    panic!("wrong input family")
                };
                assert_eq!(row.case, cases.molecules[index / count]);
                assert_eq!(row.seed, sampled.seed);
                assert_eq!(row.params, profiles[&kind][index % count]);
            }
            cases.fingerprint_profiles.insert(kind, vec![]);
            assert!(registry::validate(&cases, &[task]).is_err());
            cases
                .fingerprint_profiles
                .insert(kind, profiles[&kind].clone());
        }
        let untouched = registry::TASKS
            .iter()
            .find(|task| task.operation == registry::Operation::Fingerprint(Kind::Maccs))
            .unwrap();
        assert_eq!(
            registry::expand(&cases, untouched),
            crate::fingerprints::inputs(&cases.molecules, Kind::Maccs)
        );
        let normal: CorpusManifest = serde_json::from_value(json!({"format":"smiles"})).unwrap();
        assert!(normal.fingerprint_profiles == FingerprintProfiles::Seeded);
        let historical: CorpusManifest = serde_json::from_value(
            json!({"format":"smiles", "fingerprint_profiles":"historical0_3"}),
        )
        .unwrap();
        assert!(historical.fingerprint_profiles == FingerprintProfiles::Historical0_3);
        assert!(
            serde_json::from_value::<CorpusManifest>(
                json!({"format":"smiles", "fingerprint_profiles":"unknown"})
            )
            .is_err()
        );
    }

    #[test]
    fn historical_tautomer_selection_preserves_all_eight_recipes_and_both_censuses() {
        use registry::molecule_plan::{Profile, TaskId};
        let selected = plan(&Selection::Corpus("smiles_smoke".into())).unwrap();
        let mut cases = selected.cases;
        let original: Value = serde_json::from_slice(
            &read(&crate::root().join("tools/testdata/rdkit/tautomer_profile.json")).unwrap(),
        )
        .unwrap();
        let profiles = historical_tautomer_profiles().unwrap();
        assert_eq!(profiles.len(), 8);
        let historical: Vec<_> = original["chembl_branches"]
            .as_array()
            .unwrap()
            .iter()
            .map(|name| {
                let mut branch = original["branches"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .find(|branch| branch["name"] == *name)
                    .unwrap()
                    .clone();
                branch.as_object_mut().unwrap().remove("name");
                serde_json::from_value::<registry::molecule_plan::TautomerProfile>(branch).unwrap()
            })
            .collect();
        assert_eq!(profiles, historical);
        for id in [
            TaskId::TautomerEnumeration,
            TaskId::TautomerCanonicalization,
        ] {
            let task = registry::TASKS
                .iter()
                .find(|task| task.operation == registry::Operation::Molecular(id))
                .unwrap();
            assert_eq!(task.count(&cases), 6);
            assert_eq!(registry::expand(&cases, task).len(), 6);
        }
        cases.tautomer_profiles = Some(profiles.clone());
        for id in [
            TaskId::TautomerEnumeration,
            TaskId::TautomerCanonicalization,
        ] {
            let task = registry::TASKS
                .iter()
                .find(|task| task.operation == registry::Operation::Molecular(id))
                .unwrap();
            assert_eq!(task.count(&cases), 24);
            let expanded = registry::expand(&cases, task);
            assert_eq!(expanded.len(), 24);
            for (row, input) in expanded.iter().enumerate() {
                let Input::Molecular { case, profile } = input else {
                    panic!("wrong input family")
                };
                assert_eq!(case, &cases.molecules[row / 8]);
                let parameters = match profile {
                    Profile::TautomerEnumeration { parameters }
                    | Profile::TautomerCanonicalization { parameters } => parameters,
                    _ => panic!("wrong profile family"),
                };
                assert_eq!(parameters, &profiles[row % 8]);
            }
        }
        assert_eq!(
            cases.molecular_profiles(TaskId::NumHeavyAtoms),
            TaskId::NumHeavyAtoms.profiles()
        );
        let standard: CorpusManifest =
            serde_json::from_value(json!({"format":"smiles", "input":"smiles_smoke.smi"})).unwrap();
        assert!(standard.tautomer_profiles == TautomerProfiles::CurrentV1);
        let historical: CorpusManifest = serde_json::from_value(json!({"format":"smiles", "input":"smiles_smoke.smi", "tautomer_profiles":"historical0_3"})).unwrap();
        assert!(historical.tautomer_profiles == TautomerProfiles::Historical0_3);
        assert!(
            serde_json::from_value::<CorpusManifest>(
                json!({"format":"smiles", "tautomer_profiles":"unknown"})
            )
            .is_err()
        );
        cases.tautomer_profiles = Some(vec![]);
        let task = registry::TASKS
            .iter()
            .find(|task| {
                task.operation == registry::Operation::Molecular(TaskId::TautomerEnumeration)
            })
            .unwrap();
        assert!(registry::validate(&cases, &[task]).is_err());
    }
    #[test]
    fn unknown_and_escaping_selections_fail_instead_of_skipping() {
        assert!(plan(&Selection::Corpus("missing".into())).is_err());
        assert!(plan(&Selection::Corpus("../testdata".into())).is_err());
        assert!(plan(&Selection::Special("unknown".into())).is_err());
        assert!(input_path("../../../Cargo.toml").is_err());
    }
    #[test]
    fn special_regressions_have_separate_fixed_inputs() {
        let special = plan(&Selection::Special("all".into())).unwrap();
        assert_eq!(special.specs.len(), 12);
        assert_eq!(
            special
                .specs
                .iter()
                .map(|s| s.key())
                .collect::<BTreeSet<_>>(),
            BTreeSet::from([
                "conformer_fixed19",
                "conformer_library",
                "forcefield_properties",
                "mcs_upstream",
                "mcs_jnk1",
                "forcefield_optimizers",
                "mmff_builtin",
                "bio_mmcif_switches",
                "molalign_focused",
                "structure_tags",
                "tautomer_long_conjugated",
                "tautomer_focused"
            ])
        );
        let structure = special
            .specs
            .iter()
            .find(|spec| spec.key() == "structure_tags")
            .unwrap();
        let fixture = inputs(structure, &special.cases).unwrap();
        assert_eq!(
            fixture["cases"].as_array().unwrap().len()
                + fixture["octahedral_switch_cases"].as_array().unwrap().len(),
            77
        );
        assert!(special.cases.molecules.is_empty());
        let focused = plan(&Selection::Special("tautomer_focused".into())).unwrap();
        let fixture = inputs(&focused.specs[0], &focused.cases).unwrap();
        assert_eq!(fixture["cases"].as_array().unwrap().len(), 18);
        assert_eq!(fixture["branches"].as_array().unwrap().len(), 8);
    }
    #[test]
    fn expected_data_identity_binds_parameters_output_and_selection() {
        let plan = plan(&Selection::Corpus("smiles_smoke".into())).unwrap();
        let spec = &plan.specs[0];
        let original = identity(&plan.selection, spec, &json!([1]), b"one", 1).unwrap();
        assert_ne!(
            original,
            identity(&plan.selection, spec, &json!([2]), b"one", 1).unwrap()
        );
        assert_ne!(
            original,
            identity(&plan.selection, spec, &json!([1]), b"two", 1).unwrap()
        );
        assert_ne!(
            original,
            identity(
                &Selection::Corpus("other".into()),
                spec,
                &json!([1]),
                b"one",
                1
            )
            .unwrap()
        );
    }
    #[test]
    fn batch_preflight_requires_every_index_error_and_worker_recipe() {
        let spec = Spec::Batch;
        let input = json!([{"cases":[{"id":"a","smiles":"CCO"},{"id":"b","smiles":"invalid"}],"workers":4}]);
        let row = json!({"input":input[0],"output":{"valid_mask":[true,false],"smiles":["CCO",null],"error_indices":[1]}});
        validate_rows(&spec, &input, [Ok(row.clone())], 1).unwrap();
        let mut wrong = row.clone();
        wrong["output"]["error_indices"] = json!([0]);
        assert!(validate_rows(&spec, &input, [Ok(wrong)], 1).is_err());
        let mut wrong = row;
        wrong["input"]["workers"] = json!(1);
        assert!(validate_rows(&spec, &input, [Ok(wrong)], 1).is_err());
        assert!(validate_rows(&spec, &input, std::iter::empty(), 0).is_err());
    }
    #[test]
    fn cargo_does_not_prepare_missing_expectations() {
        let selected = plan(&Selection::Corpus("smiles_smoke".into())).unwrap();
        let isolated = tempfile::tempdir().unwrap();
        assert!(load_reference(&selected, &selected.specs[0], isolated.path()).is_err());
        assert_eq!(fs::read_dir(isolated.path()).unwrap().count(), 0);
        assert!(
            selected
                .selection
                .command()
                .contains("prepare --corpus smiles_smoke")
        );
    }
    #[test]
    fn changed_reference_bytes_and_missing_rows_fail_before_comparison() {
        let mut selected = plan(&Selection::Corpus("smiles_smoke".into())).unwrap();
        selected.cases.molecules.truncate(1);
        let spec = selected
            .specs
            .iter()
            .find(|s| s.key() == "num_heavy_atoms_smiles")
            .unwrap();
        let input = inputs(spec, &selected.cases).unwrap();
        let recipes: Vec<Input> = serde_json::from_value(input.clone()).unwrap();
        let rows: Vec<Value> = recipes
            .iter()
            .map(|recipe| {
                serde_json::to_value(Record {
                    input: recipe.clone(),
                    output: registry::Value::Molecular(crate::molecular::Outcome::Unsigned(3)),
                })
                .unwrap()
            })
            .collect();
        let output = jsonl(&rows).unwrap();
        let isolated = tempfile::tempdir().unwrap();
        fs::write(isolated.path().join("input.json"), encode(&input).unwrap()).unwrap();
        fs::write(isolated.path().join("reference.jsonl"), &output).unwrap();
        fs::write(
            isolated.path().join("manifest.json"),
            encode(&identity(&selected.selection, spec, &input, &output, rows.len()).unwrap())
                .unwrap(),
        )
        .unwrap();
        load_reference(&selected, spec, isolated.path()).unwrap();
        fs::write(isolated.path().join("reference.jsonl"), b"{}\n").unwrap();
        assert!(load_reference(&selected, spec, isolated.path()).is_err());
        fs::write(isolated.path().join("reference.jsonl"), &output).unwrap();
        fs::remove_file(isolated.path().join("manifest.json")).unwrap();
        assert!(load_reference(&selected, spec, isolated.path()).is_err());
    }
    #[test]
    fn task_selection_is_exact_and_must_apply_to_the_corpus() {
        let selection = Selection::Corpus("smiles_smoke".into());
        assert_eq!(plan_task(&selection, None).unwrap().specs.len(), 118);
        for key in [
            "num_heavy_atoms_smiles",
            "batch_smiles",
            "mmff_force_field_smiles",
            "uff_force_field_smiles",
            "murcko_scaffold_smiles",
            "net_scaffold_smiles",
            "murcko_decompose_smiles",
            "remove_hs_smiles",
        ] {
            let selected = plan_task(&selection, Some(key)).unwrap();
            assert_eq!(selected.specs.len(), 1);
            assert_eq!(selected.specs[0].key(), key);
        }
        for key in ["", "unknown", "morgan_", "bio_pdb_output_pdb"] {
            assert!(plan_task(&selection, Some(key)).is_err());
        }
        assert!(plan_task(&Selection::Special("all".into()), Some("structure_tags")).is_err());
    }

    #[test]
    fn focused_preflight_needs_only_selected_references_but_checks_all_selected() {
        let mut selected = plan_task(
            &Selection::Corpus("smiles_smoke".into()),
            Some("num_heavy_atoms_smiles"),
        )
        .unwrap();
        selected.cases.molecules.truncate(1);
        let spec = &selected.specs[0];
        let input = inputs(spec, &selected.cases).unwrap();
        let recipes: Vec<Input> = serde_json::from_value(input.clone()).unwrap();
        let rows: Vec<Value> = recipes
            .iter()
            .map(|recipe| {
                serde_json::to_value(Record {
                    input: recipe.clone(),
                    output: registry::Value::Molecular(crate::molecular::Outcome::Unsigned(3)),
                })
                .unwrap()
            })
            .collect();
        let output = jsonl(&rows).unwrap();
        let isolated = tempfile::tempdir().unwrap();
        let task_dir = isolated.path().join(spec.key());
        fs::create_dir(&task_dir).unwrap();
        fs::write(task_dir.join("input.json"), encode(&input).unwrap()).unwrap();
        fs::write(task_dir.join("reference.jsonl"), &output).unwrap();
        fs::write(
            task_dir.join("manifest.json"),
            encode(&identity(&selected.selection, spec, &input, &output, rows.len()).unwrap())
                .unwrap(),
        )
        .unwrap();
        let snapshots = load_references(&selected, isolated.path()).unwrap();
        assert_eq!(snapshots.len(), 1);
        assert_eq!(fs::read_dir(isolated.path()).unwrap().count(), 1);

        let another = plan(&selected.selection)
            .unwrap()
            .specs
            .into_iter()
            .find(|spec| spec.key() == "smiles_read_smiles")
            .unwrap();
        selected.specs.push(another);
        let error = load_references(&selected, isolated.path()).err().unwrap();
        assert!(error.contains("smiles_read_smiles"));
        assert!(error.contains("--task smiles_read_smiles"));
    }

    #[test]
    fn no_zero_thread_preparation() {
        assert!(
            prepare(&Selection::Corpus("smiles_small".into()), 0, None)
                .unwrap_err()
                .contains("positive")
        );
    }
    fn molecular_reference_fixture() -> (Plan, Spec, tempfile::TempDir, Vec<Value>) {
        let mut selected = plan_task(
            &Selection::Corpus("smiles_smoke".into()),
            Some("num_heavy_atoms_smiles"),
        )
        .unwrap();
        selected.cases.molecules.truncate(1);
        let spec = selected.specs[0];
        let input = inputs(&spec, &selected.cases).unwrap();
        let recipes: Vec<Input> = Vec::deserialize(&input).unwrap();
        let rows: Vec<Value> = recipes
            .into_iter()
            .map(|input| {
                serde_json::to_value(Record {
                    input,
                    output: registry::Value::Molecular(crate::molecular::Outcome::Unsigned(3)),
                })
                .unwrap()
            })
            .collect();
        let output = jsonl(&rows).unwrap();
        let folder = tempfile::tempdir().unwrap();
        fs::write(folder.path().join("input.json"), encode(&input).unwrap()).unwrap();
        fs::write(folder.path().join("reference.jsonl"), &output).unwrap();
        fs::write(
            folder.path().join("manifest.json"),
            encode(&identity(&selected.selection, &spec, &input, &output, rows.len()).unwrap())
                .unwrap(),
        )
        .unwrap();
        (selected, spec, folder, rows)
    }

    #[test]
    fn owned_reference_survives_expected_replacement_and_independent_cursors() {
        let (selected, spec, folder, original) = molecular_reference_fixture();
        let checked = load_reference(&selected, &spec, folder.path()).unwrap();
        assert!(checked.inputs.is_none());
        fs::write(folder.path().join("reference.jsonl"), b"changed\n").unwrap();
        fs::write(folder.path().join("input.json"), b"changed").unwrap();
        fs::remove_file(folder.path().join("manifest.json")).unwrap();
        let mut first = checked.rows.iter().unwrap();
        let mut second = checked.rows.iter().unwrap();
        for expected in &original {
            assert_eq!(first.next().unwrap().unwrap(), *expected);
            assert_eq!(second.next().unwrap().unwrap(), *expected);
        }
        assert!(first.next().is_none());
        assert!(second.next().is_none());
        drop(first);
        drop(second);
        fs::remove_dir_all(folder.path()).unwrap();
        std::thread::scope(|scope| {
            let a = scope.spawn(|| {
                checked
                    .rows
                    .iter()
                    .unwrap()
                    .collect::<Result<Vec<_>>>()
                    .unwrap()
            });
            let b = scope.spawn(|| {
                checked
                    .rows
                    .iter()
                    .unwrap()
                    .collect::<Result<Vec<_>>>()
                    .unwrap()
            });
            assert_eq!(a.join().unwrap(), original);
            assert_eq!(b.join().unwrap(), original);
        });
    }

    #[test]
    fn whole_utf8_and_json_syntax_keep_precedence_over_manifest_and_typed_rows() {
        let (selected, spec, folder, _) = molecular_reference_fixture();
        fs::write(folder.path().join("manifest.json"), b"not-json").unwrap();
        fs::write(folder.path().join("reference.jsonl"), b"{\n\xff\n").unwrap();
        let error = load_reference(&selected, &spec, folder.path())
            .err()
            .unwrap();
        assert!(
            error.contains("invalid utf-8") && error.contains("index 2"),
            "{error}"
        );
        // First {} has no Record fields, but the later syntax error must win.
        fs::write(folder.path().join("reference.jsonl"), b"{}\n{\n").unwrap();
        let error = load_reference(&selected, &spec, folder.path())
            .err()
            .unwrap();
        assert!(error.contains("EOF while parsing"), "{error}");
        assert!(!error.contains("missing field"));
        fs::write(folder.path().join("reference.jsonl"), b"{}\n").unwrap();
        let error = load_reference(&selected, &spec, folder.path())
            .err()
            .unwrap();
        assert!(
            !error.contains("missing field"),
            "manifest must fail before typed row: {error}"
        );
    }

    #[test]
    fn owned_jsonl_preserves_crlf_final_line_empty_file_and_blank_line_rules() {
        let folder = tempfile::tempdir().unwrap();
        let path = folder.path().join("rows.jsonl");
        let bytes = b"{\"a\":1}\r\n{\"b\":2}";
        fs::write(&path, bytes).unwrap();
        let rows = OwnedRows::copy_from(&path).unwrap();
        assert_eq!(rows.sha256, digest(bytes));
        assert_eq!(rows.len(), 2);
        assert_eq!(
            rows.iter().unwrap().collect::<Result<Vec<_>>>().unwrap(),
            vec![json!({"a":1}), json!({"b":2})]
        );
        fs::write(&path, b"").unwrap();
        assert_eq!(OwnedRows::copy_from(&path).unwrap().len(), 0);
        for bytes in [b"{}\n\n".as_slice(), b"\n{}\n".as_slice()] {
            fs::write(&path, bytes).unwrap();
            assert!(OwnedRows::copy_from(&path).is_err());
        }
    }

    #[test]
    fn streamed_generated_jsonl_and_manifest_hashes_equal_original_encoding() {
        let folder = tempfile::tempdir().unwrap();
        let rows = vec![
            json!({"coordinate":0.9389657496748851,"bits":u64::MAX}),
            json!({"coordinate":-0.0,"nested":{"z":[1,2],"a":"value"}}),
        ];
        let expected = jsonl(&rows).unwrap();
        let mut writer = RowWriter::new_in(folder.path()).unwrap();
        for row in &rows {
            writer.push(row.clone()).unwrap();
        }
        let generated = writer.finish().unwrap();
        assert_eq!(generated.len(), rows.len());
        assert_eq!(generated.sha256, digest(&expected));
        assert_eq!(
            encoded_digest(&rows).unwrap(),
            digest(&encode(&rows).unwrap())
        );
        let selected = plan(&Selection::Corpus("smiles_smoke".into())).unwrap();
        assert_eq!(
            identity(
                &selected.selection,
                &selected.specs[0],
                &json!([1]),
                &expected,
                rows.len()
            )
            .unwrap(),
            identity_with_output_digest(
                &selected.selection,
                &selected.specs[0],
                &json!([1]),
                generated.sha256.clone(),
                rows.len()
            )
            .unwrap(),
        );
        let path = folder.path().join("reference.jsonl");
        generated.persist(&path).unwrap();
        assert_eq!(fs::read(path).unwrap(), expected);
    }

    #[test]
    fn special_materialization_retains_complete_owned_fixture_and_rows() {
        let folder = tempfile::tempdir().unwrap();
        let path = folder.path().join("owned.jsonl");
        let rows = vec![
            json!({"case":"first","nested":[1,2]}),
            json!({"case":"second","nested":{"all":"values"}}),
        ];
        fs::write(&path, jsonl(&rows).unwrap()).unwrap();
        let owned = OwnedRows::copy_from(&path).unwrap();
        fs::remove_file(path).unwrap();
        let fixture = json!({"complete":{"cases":[1,2],"branches":["a","b"]}});
        // This boundary consumes an already checked private snapshot; source
        // schema checking remains in load_references before materialization.
        let snapshot = materialize_special(CheckedSnapshot {
            spec: Spec::Special(&registry::SPECIAL_REGRESSIONS[0]),
            inputs: Some(fixture.clone()),
            rows: owned,
        })
        .unwrap();
        assert_eq!(snapshot.inputs, fixture);
        assert_eq!(snapshot.rows, rows);
    }

    #[test]
    fn parallel_preflight_is_bounded_checks_every_selected_task_and_preserves_order() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let active = AtomicUsize::new(0);
        let peak = AtomicUsize::new(0);
        let visited: Vec<_> = (0..9).map(|_| AtomicUsize::new(0)).collect();
        let first_wave = std::sync::Barrier::new(3);
        let observed = parallel_preflight(visited.len(), 3, |index| {
            visited[index].fetch_add(1, Ordering::SeqCst);
            let running = active.fetch_add(1, Ordering::SeqCst) + 1;
            peak.fetch_max(running, Ordering::SeqCst);
            if index < 3 {
                first_wave.wait();
            }
            std::thread::yield_now();
            active.fetch_sub(1, Ordering::SeqCst);
            Ok(index)
        })
        .unwrap();
        assert_eq!(observed, (0..9).collect::<Vec<_>>());
        assert_eq!(peak.load(Ordering::SeqCst), 3);
        assert_eq!(active.load(Ordering::SeqCst), 0);
        assert!(
            visited
                .iter()
                .all(|count| count.load(Ordering::SeqCst) == 1)
        );
    }

    #[test]
    fn parallel_preflight_waits_for_all_tasks_and_returns_first_error_in_selected_order() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let visited: Vec<_> = (0..7).map(|_| AtomicUsize::new(0)).collect();
        let later_completed = (Mutex::new(false), std::sync::Condvar::new());
        let completion_order = Mutex::new(Vec::new());
        let error = parallel_preflight(visited.len(), 3, |index| {
            visited[index].fetch_add(1, Ordering::SeqCst);
            if index == 0 {
                let mut done = later_completed.0.lock().unwrap();
                while !*done {
                    done = later_completed.1.wait(done).unwrap();
                }
            }
            completion_order.lock().unwrap().push(index);
            if index == 4 {
                *later_completed.0.lock().unwrap() = true;
                later_completed.1.notify_all();
            }
            if index == 0 || index == 4 {
                Err(format!("selected task {index} failed"))
            } else {
                Ok(index)
            }
        })
        .unwrap_err();
        assert_eq!(error, "selected task 0 failed");
        assert!(
            visited
                .iter()
                .all(|count| count.load(Ordering::SeqCst) == 1)
        );
        let order = completion_order.lock().unwrap();
        assert_eq!(order.len(), visited.len());
        assert!(
            order.iter().position(|&index| index == 4).unwrap()
                < order.iter().position(|&index| index == 0).unwrap()
        );
    }
}
