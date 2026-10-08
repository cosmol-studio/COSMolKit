//! Read-only comparison. No oracle invocation or cache repair is reachable here.
use crate::{
    Input, Record, Result, Selection, directory, encode, registry,
    workflow::{self, Snapshot, Spec},
};
use serde_json::{Value, json};
use std::{collections::BTreeMap, path::PathBuf, sync::OnceLock};

struct Loaded {
    snapshots: BTreeMap<&'static str, Snapshot>,
    reports: PathBuf,
}
static CORPUS: OnceLock<Result<Loaded>> = OnceLock::new();
static SPECIAL: OnceLock<Result<BTreeMap<&'static str, Snapshot>>> = OnceLock::new();

/// Mirrors libtest name selection for preflight only. Libtest still owns
/// discovery, scheduling and execution; no tests are skipped inside their bodies.
#[derive(Default)]
struct TestFilter {
    names: Vec<String>,
    skips: Vec<String>,
    exact: bool,
}
impl TestFilter {
    fn parse(args: impl IntoIterator<Item = String>) -> Result<Self> {
        let mut selected = Self::default();
        let mut args = args.into_iter();
        while let Some(arg) = args.next() {
            let (option, inline) = arg
                .split_once('=')
                .map_or((arg.as_str(), None), |(key, value)| (key, Some(value)));
            match option {
                "--exact" => selected.exact = true,
                "--skip" | "--test-threads" | "--logfile" | "--color" | "--format"
                | "--shuffle-seed" | "-Z" => {
                    let value = inline
                        .map(str::to_owned)
                        .or_else(|| args.next())
                        .ok_or_else(|| format!("missing value for libtest option {option}"))?;
                    if option == "--skip" {
                        selected.skips.push(value);
                    }
                }
                _ if arg.starts_with('-') => {}
                _ => selected.names.push(arg),
            }
        }
        Ok(selected)
    }
    fn matches(&self, key: &str) -> bool {
        let matches = |filter: &String| {
            if self.exact {
                key == filter
            } else {
                key.contains(filter.as_str())
            }
        };
        (self.names.is_empty() || self.names.iter().any(&matches))
            && !self.skips.iter().any(matches)
    }
}

fn load_corpus() -> Result<Loaded> {
    let selection = Selection::Corpus(env!("PARITY_DEFAULT_CORPUS").into());
    let mut plan = workflow::plan(&selection)?;
    let filter = TestFilter::parse(std::env::args().skip(1))?;
    plan.specs.retain(|spec| filter.matches(spec.key()));
    let snapshots = workflow::load_references(&plan, &selection.folder()?)?;
    let parent = directory().join("reports");
    std::fs::create_dir_all(&parent).map_err(|e| e.to_string())?;
    let reports = tempfile::Builder::new()
        .prefix("corpus-")
        .tempdir_in(parent)
        .map_err(|e| e.to_string())?
        .keep();
    std::fs::write(
        reports.join("selection.json"),
        encode(&json!({"selection":selection,"tasks":snapshots.keys().collect::<Vec<_>>()}))?,
    )
    .map_err(|e| e.to_string())?;
    Ok(Loaded { snapshots, reports })
}

fn equal(expected: &Record, actual: &Record) -> bool {
    if expected.input != actual.input {
        return false;
    }
    match (&expected.output, &actual.output) {
        (registry::Value::MolAlign(a), registry::Value::MolAlign(b)) => {
            crate::molalign::matches(a, b)
        }
        (registry::Value::Mmff(a), registry::Value::Mmff(b)) => crate::mmff::matches(a, b),
        (registry::Value::Uff(a), registry::Value::Uff(b)) => {
            crate::uff::matches(&expected.input, a, b)
        }
        (registry::Value::Molecular(a), registry::Value::Molecular(b)) => {
            if matches!(
                expected.input,
                Input::Molecular {
                    profile: registry::molecule_plan::Profile::SvgDefault,
                    ..
                }
            ) {
                crate::molecular::svg_matches(a, b)
            } else {
                crate::molecular::matches(a, b)
            }
        }
        _ => expected.output == actual.output,
    }
}

fn batch(input: &Value) -> Result<Value> {
    use cosmolkit::{
        BatchErrorMode, BatchParams, BatchQueryParams, MoleculeBatch, SmilesParseParams,
        SmilesWriteParams,
    };
    let cases: Vec<registry::SmilesCase> =
        serde_json::from_value(input["cases"].clone()).map_err(|e| e.to_string())?;
    let workers = usize::try_from(input["workers"].as_u64().ok_or("missing batch workers")?)
        .map_err(|e| e.to_string())?;
    let strings: Vec<_> = cases.iter().map(|c| c.smiles.clone()).collect();
    let batch = MoleculeBatch::from_smiles_list_with_params(
        &strings,
        &SmilesParseParams::default(),
        &BatchParams {
            errors: BatchErrorMode::KeepErrors,
            n_jobs: Some(workers),
            progress_bar: Some(false),
        },
    )
    .map_err(|e| e.to_string())?;
    let mask = batch.valid_mask();
    let errors: Vec<_> = batch.errors().iter().map(|e| e.index).collect();
    let smiles = batch
        .to_smiles_list_with_params(
            &SmilesWriteParams::default(),
            &BatchQueryParams {
                n_jobs: Some(workers),
                progress_bar: Some(false),
                ..Default::default()
            },
        )
        .map_err(|e| e.to_string())?;
    let smiles = smiles
        .into_iter()
        .map(|text| {
            text.map(|text| String::from_utf8(text.into_bytes()))
                .transpose()
        })
        .collect::<std::result::Result<Vec<_>, _>>()
        .map_err(|e| e.to_string())?;
    if batch.len() != cases.len()
        || batch.valid_mask() != mask
        || batch.errors().iter().map(|e| e.index).collect::<Vec<_>>() != errors
    {
        return Err("batch query changed original records".into());
    }
    Ok(json!({"valid_mask":mask,"smiles":smiles,"error_indices":errors}))
}

pub fn run_corpus(key: &str) -> Result<()> {
    let loaded = CORPUS
        .get_or_init(load_corpus)
        .as_ref()
        .map_err(Clone::clone)?;
    let snapshot = loaded
        .snapshots
        .get(key)
        .ok_or_else(|| format!("unregistered/inapplicable test: {key}"))?;
    let mut results = Vec::new();
    for (index, row) in snapshot.rows.iter().enumerate() {
        let (matches, actual) = match snapshot.spec {
            Spec::Corpus(_) => {
                let expected: Record =
                    serde_json::from_value(row.clone()).map_err(|e| e.to_string())?;
                match crate::execute::run(&expected.input) {
                    Ok(actual) => (equal(&expected, &actual), json!({"record":actual})),
                    Err(error) => (false, json!({"error":error})),
                }
            }
            Spec::Batch => match batch(&row["input"]) {
                Ok(actual) => (actual == row["output"], json!({"output":actual})),
                Err(error) => (false, json!({"error":error})),
            },
            Spec::Special(_) => return Err("special recipe in corpus selection".into()),
        };
        results.push(json!({"index":index,"matches":matches,"expected":row,"actual":actual}));
    }
    if results.is_empty() {
        return Err("zero comparisons".into());
    }
    let failed = results.iter().filter(|r| r["matches"] != true).count();
    let report = loaded.reports.join(format!("{key}.json"));
    std::fs::write(
        &report,
        encode(&json!({"task":key,"compared":results.len(),"failed":failed,"rows":results}))?,
    )
    .map_err(|e| e.to_string())?;
    println!(
        "{key}: {}/{} matched; {}",
        results.len() - failed,
        results.len(),
        report.display()
    );
    if failed != 0 {
        return Err(format!(
            "{key}: {failed}/{} mismatches; {}",
            results.len(),
            report.display()
        ));
    }
    Ok(())
}

pub(crate) fn special_snapshot(key: &str) -> Result<&'static Snapshot> {
    let snapshots = SPECIAL
        .get_or_init(|| {
            let key = env!("PARITY_DEFAULT_SPECIAL").to_owned();
            let selection = Selection::Special(key);
            workflow::load_references(&workflow::plan(&selection)?, &selection.folder()?)
        })
        .as_ref()
        .map_err(Clone::clone)?;
    snapshots
        .get(key)
        .ok_or_else(|| format!("special regression {key} is not selected"))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn persistent_force_field_comparison_rejects_one_ulp_in_every_float_field() {
        use crate::persistent_forcefields::{Input as FieldInput, Kind, Observation, Snapshot};
        let state = Snapshot {
            energy_bits: 1.0_f64.to_bits(),
            gradient_bits: vec![[1.0_f64.to_bits(); 3]],
            positions_bits: vec![[1.0_f64.to_bits(); 3]],
        };
        let expected = Record {
            input: Input::PersistentForceField(FieldInput::new(
                registry::SmilesCase {
                    id: "comparison".into(),
                    smiles: "C".into(),
                },
                Kind::Mmff,
            )),
            output: registry::Value::PersistentForceField(Observation::Evaluated {
                initial: state.clone(),
                final_state: state,
                converged: false,
            }),
        };
        assert!(equal(&expected, &expected));
        for phase in 0..2 {
            for field in 0..3 {
                let mut actual = expected.clone();
                let registry::Value::PersistentForceField(Observation::Evaluated {
                    initial,
                    final_state,
                    ..
                }) = &mut actual.output
                else {
                    unreachable!()
                };
                let state = if phase == 0 { initial } else { final_state };
                match field {
                    0 => state.energy_bits += 1,
                    1 => state.gradient_bits[0][0] += 1,
                    _ => state.positions_bits[0][0] += 1,
                }
                assert!(!equal(&expected, &actual), "phase {phase}, field {field}");
            }
        }
    }
    fn filter(args: &[&str]) -> TestFilter {
        TestFilter::parse(args.iter().map(|s| (*s).to_owned())).unwrap()
    }
    #[test]
    fn libtest_default_substring_exact_and_multiple_filters_select_same_tasks() {
        assert!(filter(&[]).matches("morgan_fingerprint_smiles"));
        let partial = filter(&["morgan_"]);
        assert!(partial.matches("morgan_fingerprint_smiles"));
        assert!(partial.matches("morgan_count_fingerprint_smiles"));
        assert!(!partial.matches("smiles_read_smiles"));
        let exact = filter(&["morgan_fingerprint_smiles", "--exact"]);
        assert!(exact.matches("morgan_fingerprint_smiles"));
        assert!(!exact.matches("morgan_fingerprint_smiles_extra"));
        let multiple = filter(&["morgan_", "smiles_read"]);
        assert!(multiple.matches("smiles_read_smiles"));
        assert!(multiple.matches("morgan_count_fingerprint_smiles"));
        assert!(!multiple.matches("num_heavy_atoms_smiles"));
    }
    #[test]
    fn libtest_skip_respects_exact_and_option_values_are_not_filters() {
        let partial = filter(&["morgan_", "--skip", "count", "--skip=sparse"]);
        assert!(partial.matches("morgan_fingerprint_smiles"));
        assert!(!partial.matches("morgan_count_fingerprint_smiles"));
        assert!(!partial.matches("morgan_sparse_fingerprint_smiles"));
        let exact = filter(&["--exact", "--skip", "morgan_"]);
        assert!(exact.matches("morgan_fingerprint_smiles"));
        let options = filter(&[
            "--test-threads",
            "4",
            "--format=pretty",
            "--color",
            "never",
            "--logfile",
            "morgan_",
            "--shuffle-seed",
            "123",
            "-Z",
            "unstable-options",
            "--nocapture",
            "--show-output",
        ]);
        assert!(options.matches("smiles_read_smiles"));
        assert!(options.matches("morgan_fingerprint_smiles"));
        assert!(TestFilter::parse(["--skip".to_owned()]).is_err());
    }
    #[test]
    fn input_identity_is_not_ignored_by_comparison() {
        let a = Record {
            input: Input::Molecular {
                case: registry::SmilesCase {
                    id: "a".into(),
                    smiles: "CCO".into(),
                },
                profile: registry::molecule_plan::Profile::NumHeavyAtoms {
                    remove_hydrogens: true,
                },
            },
            output: registry::Value::Molecular(crate::molecular::Outcome::Unsigned(3)),
        };
        let mut b = a.clone();
        assert!(equal(&a, &b));
        if let Input::Molecular { case, .. } = &mut b.input {
            case.id = "b".into();
        }
        assert!(!equal(&a, &b));
    }
}
