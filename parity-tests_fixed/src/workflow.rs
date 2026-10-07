use crate::{
    Corpus, Input, Record, Result, Task, digest, directory, encode, expected, read, reference,
    registry, special_regression,
};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    path::{Path, PathBuf},
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
            serde_json::from_slice(&read(&directory().join("testdata").join(s.fixture))?)
                .map_err(|e| e.to_string())
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

pub(crate) fn validate_rows(spec: &Spec, input: &Value, rows: &[Value]) -> Result<()> {
    match *spec {
        Spec::Corpus(task) => {
            let inputs: Vec<Input> =
                serde_json::from_value(input.clone()).map_err(|e| e.to_string())?;
            if inputs.is_empty() || inputs.len() != rows.len() {
                return Err("reference row count mismatch".into());
            }
            for (recipe, row) in inputs.iter().zip(rows) {
                let record: Record =
                    serde_json::from_value(row.clone()).map_err(|e| e.to_string())?;
                task.validate_reference(recipe, &record.input, &record.output)?;
            }
        }
        Spec::Batch => {
            let recipes = input.as_array().ok_or("batch recipes must be an array")?;
            if recipes.is_empty() || recipes.len() != rows.len() {
                return Err("batch row count mismatch".into());
            }
            for (recipe, row) in recipes.iter().zip(rows) {
                let count = recipe["cases"]
                    .as_array()
                    .ok_or("batch cases missing")?
                    .len();
                let mask: Vec<bool> = serde_json::from_value(row["output"]["valid_mask"].clone())
                    .map_err(|e| e.to_string())?;
                let smiles: Vec<Option<String>> =
                    serde_json::from_value(row["output"]["smiles"].clone())
                        .map_err(|e| e.to_string())?;
                let errors: Vec<usize> =
                    serde_json::from_value(row["output"]["error_indices"].clone())
                        .map_err(|e| e.to_string())?;
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
            special_regression::validate(&encode(input)?, &jsonl(rows)?, s.rows, s.schema)?;
        }
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
    let pin = if matches!(*spec, Spec::Corpus(t) if t.operation == registry::Operation::BioPdbOutput)
    {
        "gemmi.json"
    } else {
        "rdkit.json"
    };
    Ok(Manifest {
        schema: 1,
        selection: selection.clone(),
        task: spec.key().into(),
        input_sha256: digest(&encode(input)?),
        generator_sha256: reference::source_digest(spec)?,
        reference_identity: serde_json::from_slice(&read(
            &directory().join("testdata/reference").join(pin),
        )?)
        .map_err(|e| e.to_string())?,
        platform: format!("{}-{}", std::env::consts::ARCH, std::env::consts::OS),
        output_sha256: digest(output),
        rows,
    })
}

fn load_reference(plan: &Plan, spec: &Spec, dir: &Path) -> Result<Snapshot> {
    let input = inputs(spec, &plan.cases)?;
    let stored: Value =
        serde_json::from_slice(&read(&dir.join("input.json"))?).map_err(|e| e.to_string())?;
    if stored != input {
        return Err("prepared input changed".into());
    }
    let output = read(&dir.join("reference.jsonl"))?;
    let rows: Vec<Value> = String::from_utf8(output.clone())
        .map_err(|e| e.to_string())?
        .lines()
        .map(|line| serde_json::from_str(line).map_err(|e| e.to_string()))
        .collect::<Result<_>>()?;
    let manifest: Manifest =
        serde_json::from_slice(&read(&dir.join("manifest.json"))?).map_err(|e| e.to_string())?;
    if manifest != identity(&plan.selection, spec, &input, &output, rows.len())? {
        return Err("stale/corrupt reference manifest".into());
    }
    validate_rows(spec, &input, &rows)?;
    Ok(Snapshot {
        spec: *spec,
        inputs: input,
        rows,
    })
}

pub(crate) fn load_references(
    plan: &Plan,
    folder: &Path,
) -> Result<BTreeMap<&'static str, Snapshot>> {
    if plan.specs.is_empty() {
        return Err("empty test selection".into());
    }
    plan.specs
        .iter()
        .map(|s| {
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
        .collect()
}

pub fn prepare(selection: &Selection, threads: usize, task: Option<&str>) -> Result<Preparation> {
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
        let generator_before = reference::source_digest(spec)?;
        eprintln!("Generating {}", spec.key());
        let rows = reference::generate(spec, &plan.cases, &input, threads, &mut descriptors)?;
        if reference::source_digest(spec)? != generator_before {
            return Err("reference generator changed during preparation".into());
        }
        validate_rows(spec, &input, &rows)?;
        let bytes = jsonl(&rows)?;
        let temporary = tempfile::Builder::new()
            .prefix(".prepare-")
            .tempdir_in(&folder)
            .map_err(|e| e.to_string())?;
        for (name, bytes) in [
            ("input.json", encode(&input)?),
            ("reference.jsonl", bytes.clone()),
            (
                "manifest.json",
                encode(&identity(selection, spec, &input, &bytes, rows.len())?)?,
            ),
        ] {
            fs::write(temporary.path().join(name), bytes).map_err(|e| e.to_string())?;
        }
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
        result.rows += rows.len();
        eprintln!(
            "  [========================] saved {} reference rows",
            rows.len()
        );
    }
    load_references(&plan, &folder)?;
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
    use super::*;
    #[test]
    fn registry_keys_and_cargo_declarations_share_one_census() {
        assert_eq!(registry::TASKS.len(), 108);
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
        assert_eq!(smoke.specs.len(), 107);
        let bio = plan(&Selection::Corpus("bio_small".into())).unwrap();
        assert_eq!(bio.cases.bio_cases.len(), 2);
        assert_eq!(bio.specs.len(), 2);
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
        assert_eq!(special.specs.len(), 2);
        let fixture = inputs(&special.specs[0], &special.cases).unwrap();
        assert_eq!(
            fixture["cases"].as_array().unwrap().len()
                + fixture["octahedral_switch_cases"].as_array().unwrap().len(),
            77
        );
        assert!(special.cases.molecules.is_empty());
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
        validate_rows(&spec, &input, &[row.clone()]).unwrap();
        let mut wrong = row.clone();
        wrong["output"]["error_indices"] = json!([0]);
        assert!(validate_rows(&spec, &input, &[wrong]).is_err());
        let mut wrong = row;
        wrong["input"]["workers"] = json!(1);
        assert!(validate_rows(&spec, &input, &[wrong]).is_err());
        assert!(validate_rows(&spec, &input, &[]).is_err());
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
        assert_eq!(plan_task(&selection, None).unwrap().specs.len(), 107);
        for key in ["num_heavy_atoms_smiles", "batch_smiles"] {
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

        let another = plan(&selected.selection).unwrap().specs[0];
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
}
