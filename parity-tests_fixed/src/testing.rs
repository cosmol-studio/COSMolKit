//! Read-only comparison. No oracle invocation or cache repair is reachable here.
use crate::{
    Input, Record, Result, Selection, directory, encode, registry,
    workflow::{self, CheckedSnapshot, Snapshot, Spec},
};
use serde::Deserialize;
use serde_json::{Value, json};
use std::{collections::BTreeMap, io::Write, path::PathBuf, sync::OnceLock};

struct Loaded {
    snapshots: BTreeMap<&'static str, CheckedSnapshot>,
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
        (registry::Value::Fingerprint(a), registry::Value::Fingerprint(b)) => {
            crate::fingerprints::matches(a, b)
        }
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
            errors: Some(BatchErrorMode::KeepErrors),
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

fn preparation_timeout(record: &Record) -> Option<&crate::uff::PreparationTimeout> {
    use crate::uff::GeometryPreparation;
    match (&record.input, &record.output) {
        (Input::Uff(input), registry::Value::Uff(crate::uff::Observation::TimedOut(timeout)))
            if input.preparation.as_ref()
                == Some(&GeometryPreparation::TimedOut(timeout.clone())) =>
        {
            Some(timeout)
        }
        (
            Input::Mmff(input),
            registry::Value::Mmff(crate::mmff::Observation::TimedOut(timeout)),
        ) if input.preparation.as_ref()
            == Some(&GeometryPreparation::TimedOut(timeout.clone())) =>
        {
            Some(timeout)
        }
        _ => None,
    }
}

// User-approved exception: no reference value exists when this pinned native
// Layered branch crashes. Keep the complete observation, never call it equal,
// and do not extend this disposition to ordinary errors or killed processes.
fn known_upstream_crash(record: &Record) -> bool {
    use crate::fingerprints::{Observation, Params, Roots};
    matches!(
        (&record.input, &record.output),
        (
            Input::Fingerprint(crate::fingerprints::FingerprintInput {
                params: Params::Layered {
                    branched: false,
                    roots: Roots::All,
                    ..
                },
                ..
            }),
            registry::Value::Fingerprint(Observation::ReferenceProcessFailure {
                exit_code: -11 | 0xc000_0005,
                process_id: 1..,
                ..
            })
        )
    )
}

#[derive(serde::Serialize)]
#[serde(untagged)]
enum ReportActual<'a> {
    Record { record: &'a Record },
    Error { error: &'a str },
    Output { output: &'a Value },
}

#[derive(serde::Serialize)]
struct ReportRow<'a> {
    // Keep the compact report object's original alphabetical key order.
    actual: Option<ReportActual<'a>>,
    expected: &'a Value,
    index: usize,
    matches: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    skipped: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    stage: Option<&'static str>,
    #[serde(skip_serializing_if = "Option::is_none")]
    timeout: Option<&'a crate::uff::PreparationTimeout>,
}

/// Retain every complete row on private disk, including failures/timeouts.
/// Summary publication copies the array, without constructing a report Value.
struct CorpusReport {
    spool: tempfile::NamedTempFile,
    writer: std::io::BufWriter<std::fs::File>,
    total: usize,
    failed: usize,
    timed_out: usize,
    upstream_crashed: usize,
}
impl CorpusReport {
    fn new() -> Result<Self> {
        let artifacts = directory().join("reports");
        std::fs::create_dir_all(&artifacts).map_err(|e| e.to_string())?;
        let spool = tempfile::NamedTempFile::new_in(artifacts).map_err(|e| e.to_string())?;
        let mut writer = std::io::BufWriter::new(spool.reopen().map_err(|e| e.to_string())?);
        writer.write_all(b"[").map_err(|e| e.to_string())?;
        Ok(Self {
            spool,
            writer,
            total: 0,
            failed: 0,
            timed_out: 0,
            upstream_crashed: 0,
        })
    }
    fn push<T: serde::Serialize + ?Sized>(
        &mut self,
        row: &T,
        matches: Option<bool>,
        skipped: bool,
        upstream_crashed: bool,
    ) -> Result<()> {
        if self.total != 0 {
            self.writer.write_all(b",").map_err(|e| e.to_string())?;
        }
        serde_json::to_writer(&mut self.writer, row).map_err(|e| e.to_string())?;
        self.total += 1;
        self.failed += usize::from(matches == Some(false));
        self.timed_out += usize::from(skipped);
        self.upstream_crashed += usize::from(upstream_crashed);
        Ok(())
    }
    fn finish(mut self, key: &str, report: &std::path::Path) -> Result<()> {
        self.writer.write_all(b"]").map_err(|e| e.to_string())?;
        self.writer.flush().map_err(|e| e.to_string())?;
        drop(self.writer);
        let compared = self.total - self.timed_out - self.upstream_crashed;
        let mut output =
            std::io::BufWriter::new(std::fs::File::create(report).map_err(|e| e.to_string())?);
        write!(
            output,
            "{{\"compared\":{compared},\"failed\":{},\"rows\":",
            self.failed
        )
        .map_err(|e| e.to_string())?;
        std::io::copy(
            &mut self.spool.reopen().map_err(|e| e.to_string())?,
            &mut output,
        )
        .map_err(|e| e.to_string())?;
        output.write_all(b",\"task\":").map_err(|e| e.to_string())?;
        serde_json::to_writer(&mut output, key).map_err(|e| e.to_string())?;
        write!(
            output,
            ",\"timed_out\":{},\"upstream_crashed\":{},\"total\":{}}}",
            self.timed_out, self.upstream_crashed, self.total
        )
        .map_err(|e| e.to_string())?;
        output.flush().map_err(|e| e.to_string())?;
        drop(output);
        println!(
            "{key}: {}/{} matched; {} timed out; {} upstream crashes (not compared); {}",
            compared - self.failed,
            compared,
            self.timed_out,
            self.upstream_crashed,
            report.display()
        );
        // A report I/O error has already returned. Keep complete diagnostics
        // for empty/all-timeout/mismatching cases before rejecting them.
        if compared == 0 {
            return Err(format!(
                "{key}: zero comparisons; {} timed out; {} upstream crashes; {}",
                self.timed_out,
                self.upstream_crashed,
                report.display()
            ));
        }
        if self.failed != 0 {
            return Err(format!(
                "{key}: {}/{compared} mismatches; {} timed out; {}",
                self.failed,
                self.timed_out,
                report.display()
            ));
        }
        Ok(())
    }
}

#[cfg(test)]
fn write_corpus_report(key: &str, results: &[Value], report: &std::path::Path) -> Result<()> {
    let mut writer = CorpusReport::new()?;
    for row in results {
        writer.push(
            row,
            row["matches"].as_bool(),
            row["skipped"] == true && row["stage"] != "UpstreamReferenceCrash",
            row["stage"] == "UpstreamReferenceCrash",
        )?;
    }
    writer.finish(key, report)
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
    let mut report = CorpusReport::new()?;
    for (index, row) in snapshot.rows.iter()?.enumerate() {
        let row = row?;
        match snapshot.spec {
            Spec::Corpus(_) => {
                // Deserialize the borrowed original Value directly: do not
                // create an intermediate cloned JSON tree for each record.
                let expected = Record::deserialize(&row).map_err(|e| e.to_string())?;
                if known_upstream_crash(&expected) {
                    let actual = crate::execute::run(&expected.input);
                    report.push(
                        &ReportRow {
                            actual: Some(match &actual {
                                Ok(record) => ReportActual::Record { record },
                                Err(error) => ReportActual::Error { error },
                            }),
                            expected: &row,
                            index,
                            matches: None,
                            skipped: Some(true),
                            stage: Some("UpstreamReferenceCrash"),
                            timeout: None,
                        },
                        None,
                        false,
                        true,
                    )?;
                    continue;
                }
                if let Some(timeout) = preparation_timeout(&expected) {
                    report.push(
                        &ReportRow {
                            actual: None,
                            expected: &row,
                            index,
                            matches: None,
                            skipped: Some(true),
                            stage: Some("Preparation"),
                            timeout: Some(timeout),
                        },
                        None,
                        true,
                        false,
                    )?;
                    continue;
                }
                match crate::execute::run(&expected.input) {
                    Ok(actual) => {
                        let matches = equal(&expected, &actual);
                        report.push(
                            &ReportRow {
                                actual: Some(ReportActual::Record { record: &actual }),
                                expected: &row,
                                index,
                                matches: Some(matches),
                                skipped: None,
                                stage: None,
                                timeout: None,
                            },
                            Some(matches),
                            false,
                            false,
                        )?;
                    }
                    Err(error) => report.push(
                        &ReportRow {
                            actual: Some(ReportActual::Error { error: &error }),
                            expected: &row,
                            index,
                            matches: Some(false),
                            skipped: None,
                            stage: None,
                            timeout: None,
                        },
                        Some(false),
                        false,
                        false,
                    )?,
                }
            }
            Spec::Batch => match batch(&row["input"]) {
                Ok(actual) => {
                    let matches = actual == row["output"];
                    report.push(
                        &ReportRow {
                            actual: Some(ReportActual::Output { output: &actual }),
                            expected: &row,
                            index,
                            matches: Some(matches),
                            skipped: None,
                            stage: None,
                            timeout: None,
                        },
                        Some(matches),
                        false,
                        false,
                    )?;
                }
                Err(error) => report.push(
                    &ReportRow {
                        actual: Some(ReportActual::Error { error: &error }),
                        expected: &row,
                        index,
                        matches: Some(false),
                        skipped: None,
                        stage: None,
                        timeout: None,
                    },
                    Some(false),
                    false,
                    false,
                )?,
            },
            Spec::Special(_) => return Err("special recipe in corpus selection".into()),
        }
    }
    report.finish(key, &loaded.reports.join(format!("{key}.json")))
}

pub(crate) fn special_snapshot(key: &str) -> Result<&'static Snapshot> {
    let snapshots = SPECIAL
        .get_or_init(|| {
            let key = env!("PARITY_DEFAULT_SPECIAL").to_owned();
            let selection = Selection::Special(key);
            workflow::load_special_references(&workflow::plan(&selection)?, &selection.folder()?)
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
    fn native_reference_crashes_never_compare_equal_and_are_separately_accounted() {
        use crate::fingerprints::{Counts, FingerprintInput, Mask, Observation, Params, Roots};
        let expected = Record {
            input: Input::Fingerprint(FingerprintInput {
                case: registry::SmilesCase {
                    id: "fixed-CCO".into(),
                    smiles: "CCO".into(),
                },
                right: None,
                seed: 0,
                params: Params::Layered {
                    layers: 63,
                    min_path: 1,
                    max_path: 7,
                    fp_size: 2048,
                    branched: false,
                    roots: Roots::All,
                    counts: Counts::Seeded,
                    mask: Mask::Absent,
                },
            }),
            output: registry::Value::Fingerprint(Observation::ReferenceProcessFailure {
                exit_code: -11,
                process_id: 123,
                stderr: "Segmentation fault".into(),
            }),
        };
        assert!(preparation_timeout(&expected).is_none());
        assert!(known_upstream_crash(&expected));
        assert!(!equal(&expected, &expected));
        for code in [-9, 0, 1] {
            let mut other = expected.clone();
            let registry::Value::Fingerprint(Observation::ReferenceProcessFailure {
                exit_code,
                ..
            }) = &mut other.output
            else {
                unreachable!()
            };
            *exit_code = code;
            assert!(!known_upstream_crash(&other));
        }
        let error = crate::execute::run(&expected.input).unwrap_err();
        assert_eq!(error, "enumerated path contains invalid bond index");
        // Also retain the returned-value comparison guard: a native process
        // failure must never match a CK result, whether it is a value or error.
        let mut value_expected = expected.clone();
        let Input::Fingerprint(FingerprintInput {
            params: Params::Layered { branched, .. },
            ..
        }) = &mut value_expected.input
        else {
            panic!("expected the Layered fixture");
        };
        *branched = true;
        assert!(!known_upstream_crash(&value_expected));
        let actual = crate::execute::run(&value_expected.input).unwrap();
        assert!(!equal(&value_expected, &actual));
        let folder = tempfile::tempdir().unwrap();
        let report = folder.path().join("native-failure.json");
        let rows = [
            json!({"index":0,"matches":false,
                   "expected":expected,"actual":{"error":error}}),
            json!({"index":1,"matches":equal(&value_expected,&actual),
                   "expected":value_expected,"actual":{"record":actual}}),
        ];
        assert!(write_corpus_report("layered_fingerprint_smiles", &rows, &report).is_err());
        let saved: Value = serde_json::from_slice(&std::fs::read(&report).unwrap()).unwrap();
        assert_eq!(saved["total"], 2);
        assert_eq!(saved["compared"], 2);
        assert_eq!(saved["failed"], 2);
        assert_eq!(saved["timed_out"], 0);
        assert_eq!(saved["rows"], json!(rows));
        assert_eq!(saved["upstream_crashed"], 0);

        let rows = [
            json!({"index":0,"matches":null,"skipped":true,
                   "stage":"UpstreamReferenceCrash",
                   "expected":expected,"actual":{"error":error}}),
            json!({"index":1,"matches":true,"actual":{"record":actual}}),
        ];
        write_corpus_report("fingerprint_layered_smiles", &rows, &report).unwrap();
        let saved: Value = serde_json::from_slice(&std::fs::read(&report).unwrap()).unwrap();
        assert_eq!(saved["total"], 2);
        assert_eq!(saved["compared"], 1);
        assert_eq!(saved["failed"], 0);
        assert_eq!(saved["timed_out"], 0);
        assert_eq!(saved["upstream_crashed"], 1);
        assert_eq!(saved["rows"], json!(rows));
        assert!(write_corpus_report("fingerprint_layered_smiles", &rows[..1], &report).is_err());
    }

    #[test]
    fn preparation_timeouts_require_matching_typed_metadata_and_optimization_profiles() {
        use crate::uff::{GeometryPreparation, PreparationTimeout, TimeoutMechanism};
        for mechanism in [TimeoutMechanism::Native, TimeoutMechanism::ProcessDeadline] {
            let timeout = PreparationTimeout {
                limit_seconds: 60,
                mechanism,
            };
            let case = registry::SmilesCase {
                id: "timeout-case".into(),
                smiles: "CCO".into(),
            };
            for operation in [
                registry::Operation::UffOptimization,
                registry::Operation::UffConformerOptimization,
            ] {
                let mut recipe = crate::uff::UffInput {
                    case: case.clone(),
                    profile: crate::uff::profiles(operation)[0],
                    preparation: None,
                };
                let mut prepared = recipe.clone();
                prepared.preparation = Some(GeometryPreparation::TimedOut(timeout.clone()));
                let output = crate::uff::Observation::TimedOut(timeout.clone());
                crate::uff::validate_reference(&recipe, &prepared, &output).unwrap();
                let record = Record {
                    input: Input::Uff(prepared.clone()),
                    output: registry::Value::Uff(output.clone()),
                };
                assert_eq!(preparation_timeout(&record), Some(&timeout));
                let mut wrong = timeout.clone();
                wrong.limit_seconds = 59;
                assert!(
                    crate::uff::validate_reference(
                        &recipe,
                        &prepared,
                        &crate::uff::Observation::TimedOut(wrong.clone())
                    )
                    .is_err()
                );
                prepared.preparation = Some(GeometryPreparation::TimedOut(wrong.clone()));
                assert!(
                    crate::uff::validate_reference(
                        &recipe,
                        &prepared,
                        &crate::uff::Observation::TimedOut(wrong)
                    )
                    .is_err()
                );
                prepared.preparation = Some(GeometryPreparation::TimedOut(timeout.clone()));
                recipe.profile = crate::uff::Profile::Coverage {
                    add_hydrogens: true,
                };
                prepared.profile = recipe.profile;
                assert!(crate::uff::validate_reference(&recipe, &prepared, &output).is_err());
            }
            for operation in [
                registry::Operation::MmffOptimization,
                registry::Operation::MmffConformerOptimization,
            ] {
                for profile in crate::mmff::profiles(operation) {
                    let mut recipe = crate::mmff::MmffInput {
                        case: case.clone(),
                        profile,
                        preparation: None,
                    };
                    let mut prepared = recipe.clone();
                    prepared.preparation = Some(GeometryPreparation::TimedOut(timeout.clone()));
                    let output = crate::mmff::Observation::TimedOut(timeout.clone());
                    crate::mmff::validate_reference(&recipe, &prepared, &output).unwrap();
                    let record = Record {
                        input: Input::Mmff(prepared.clone()),
                        output: registry::Value::Mmff(output.clone()),
                    };
                    assert_eq!(preparation_timeout(&record), Some(&timeout));
                    prepared.case.id = "wrong-case".into();
                    assert!(crate::mmff::validate_reference(&recipe, &prepared, &output).is_err());
                    prepared.case = case.clone();
                    let mut wrong = timeout.clone();
                    wrong.limit_seconds = 59;
                    assert!(
                        crate::mmff::validate_reference(
                            &recipe,
                            &prepared,
                            &crate::mmff::Observation::TimedOut(wrong)
                        )
                        .is_err()
                    );
                    recipe.profile = crate::mmff::Profile::Coverage {
                        add_hydrogens: true,
                    };
                    prepared.profile = recipe.profile;
                    assert!(crate::mmff::validate_reference(&recipe, &prepared, &output).is_err());
                }
            }
        }
    }

    #[test]
    fn timeout_reports_preserve_rows_and_exclude_only_explicit_timeouts() {
        let folder = tempfile::tempdir().unwrap();
        let path = folder.path().join("mixed.json");
        let rows = vec![
            json!({"index":0,"matches":null,"skipped":true,"timeout":{"limit_seconds":60,"mechanism":"ProcessDeadline"}}),
            json!({"index":1,"matches":true}),
            json!({"index":2,"matches":false,"actual":{"error":"ordinary preparation rejection"}}),
        ];
        assert!(
            write_corpus_report("mixed", &rows, &path)
                .unwrap_err()
                .contains("1/2 mismatches")
        );
        let report: Value = serde_json::from_slice(&std::fs::read(&path).unwrap()).unwrap();
        assert_eq!(report["total"], 3);
        assert_eq!(report["compared"], 2);
        assert_eq!(report["failed"], 1);
        assert_eq!(report["timed_out"], 1);
        assert_eq!(report["rows"], json!(rows));
        assert!(
            write_corpus_report("all-timeout", &rows[..1], &path)
                .unwrap_err()
                .contains("zero comparisons")
        );
        let report: Value = serde_json::from_slice(&std::fs::read(&path).unwrap()).unwrap();
        assert_eq!(report["compared"], 0);
        assert_eq!(report["failed"], 0);
        assert_eq!(report["timed_out"], 1);
        assert_eq!(report["rows"].as_array().unwrap().len(), 1);
    }

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
        assert!(filter(&[]).matches("fingerprint_morgan_smiles"));
        let partial = filter(&["morgan_"]);
        assert!(partial.matches("fingerprint_morgan_smiles"));
        assert!(partial.matches("fingerprint_morgan_count_smiles"));
        assert!(!partial.matches("smiles_read_smiles"));
        let exact = filter(&["fingerprint_morgan_smiles", "--exact"]);
        assert!(exact.matches("fingerprint_morgan_smiles"));
        assert!(!exact.matches("fingerprint_morgan_smiles_extra"));
        let multiple = filter(&["morgan_", "smiles_read"]);
        assert!(multiple.matches("smiles_read_smiles"));
        assert!(multiple.matches("fingerprint_morgan_count_smiles"));
        assert!(!multiple.matches("num_heavy_atoms_smiles"));
    }
    #[test]
    fn libtest_skip_respects_exact_and_option_values_are_not_filters() {
        let partial = filter(&["morgan_", "--skip", "count", "--skip=sparse"]);
        assert!(partial.matches("fingerprint_morgan_smiles"));
        assert!(!partial.matches("fingerprint_morgan_count_smiles"));
        assert!(!partial.matches("fingerprint_morgan_sparse_smiles"));
        let exact = filter(&["--exact", "--skip", "morgan_"]);
        assert!(exact.matches("fingerprint_morgan_smiles"));
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
        assert!(options.matches("fingerprint_morgan_smiles"));
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
                profile: registry::molecule_plan::Profile::NumHeavyAtoms { remove_hs: true },
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
    #[test]
    fn borrowed_stream_report_retains_complete_expected_actual_and_array_order() {
        let folder = tempfile::tempdir().unwrap();
        let path = folder.path().join("all-rows.json");
        let expected = json!({"nested":{"text":"complete\nvalue","coordinates":[1,2,3]},
                              "native_bits":u64::MAX});
        let actual = json!({"data":[{"z":2,"a":1}],"native_bits":u64::MAX});
        let error = "ordinary source error";
        let timeout = crate::uff::PreparationTimeout {
            limit_seconds: 60,
            mechanism: crate::uff::TimeoutMechanism::ProcessDeadline,
        };
        let mut report = CorpusReport::new().unwrap();
        report
            .push(
                &ReportRow {
                    actual: Some(ReportActual::Output { output: &actual }),
                    expected: &expected,
                    index: 0,
                    matches: Some(true),
                    skipped: None,
                    stage: None,
                    timeout: None,
                },
                Some(true),
                false,
                false,
            )
            .unwrap();
        report
            .push(
                &ReportRow {
                    actual: Some(ReportActual::Error { error }),
                    expected: &expected,
                    index: 1,
                    matches: Some(false),
                    skipped: None,
                    stage: None,
                    timeout: None,
                },
                Some(false),
                false,
                false,
            )
            .unwrap();
        report
            .push(
                &ReportRow {
                    actual: None,
                    expected: &expected,
                    index: 2,
                    matches: None,
                    skipped: Some(true),
                    stage: Some("Preparation"),
                    timeout: Some(&timeout),
                },
                None,
                true,
                false,
            )
            .unwrap();
        assert!(
            report
                .finish("full-task", &path)
                .unwrap_err()
                .contains("1/2 mismatches")
        );
        let saved: Value = serde_json::from_slice(&std::fs::read(path).unwrap()).unwrap();
        assert_eq!(
            saved,
            json!({"task":"full-task","total":3,"compared":2,
            "failed":1,"timed_out":1,"upstream_crashed":0,"rows":[
                {"index":0,"matches":true,"expected":expected,"actual":{"output":actual}},
                {"index":1,"matches":false,"expected":expected,"actual":{"error":error}},
                {"index":2,"matches":null,"expected":expected,"actual":null,
                 "skipped":true,"stage":"Preparation","timeout":timeout}
            ]})
        );
    }

    #[test]
    fn report_io_error_precedes_zero_or_mismatch_rejection_and_empty_is_recorded() {
        let folder = tempfile::tempdir().unwrap();
        let impossible = folder.path().join("missing-parent/report.json");
        let error = CorpusReport::new()
            .unwrap()
            .finish("empty", &impossible)
            .unwrap_err();
        assert!(!error.contains("zero comparisons"));
        let mut failed = CorpusReport::new().unwrap();
        failed
            .push(
                &json!({"index":0,"matches":false}),
                Some(false),
                false,
                false,
            )
            .unwrap();
        let error = failed.finish("mismatch", &impossible).unwrap_err();
        assert!(!error.contains("mismatches"));
        let path = folder.path().join("empty.json");
        assert!(
            CorpusReport::new()
                .unwrap()
                .finish("empty", &path)
                .unwrap_err()
                .contains("zero comparisons")
        );
        let saved: Value = serde_json::from_slice(&std::fs::read(path).unwrap()).unwrap();
        assert_eq!(
            saved,
            json!({"task":"empty","total":0,"compared":0,
                "failed":0,"timed_out":0,"upstream_crashed":0,"rows":[]})
        );
    }

    #[test]
    fn borrowed_typed_actual_report_preserves_complete_record_schema() {
        let actual = Record {
            input: Input::Molecular {
                case: registry::SmilesCase {
                    id: "one".into(),
                    smiles: "CCO".into(),
                },
                profile: registry::molecule_plan::Profile::NumHeavyAtoms { remove_hs: true },
            },
            output: registry::Value::Molecular(crate::molecular::Outcome::Unsigned(3)),
        };
        let expected = serde_json::to_value(&actual).unwrap();
        let folder = tempfile::tempdir().unwrap();
        let path = folder.path().join("typed.json");
        let mut report = CorpusReport::new().unwrap();
        report
            .push(
                &ReportRow {
                    actual: Some(ReportActual::Record { record: &actual }),
                    expected: &expected,
                    index: 0,
                    matches: Some(true),
                    skipped: None,
                    stage: None,
                    timeout: None,
                },
                Some(true),
                false,
                false,
            )
            .unwrap();
        report.finish("typed", &path).unwrap();
        let saved: Value = serde_json::from_slice(&std::fs::read(path).unwrap()).unwrap();
        assert_eq!(
            saved["rows"][0],
            json!({"actual":{"record":actual},
                                          "expected":expected,"index":0,"matches":true})
        );
        assert_eq!(saved["total"], 1);
        assert_eq!(saved["compared"], 1);
        assert_eq!(saved["failed"], 0);
        assert_eq!(saved["timed_out"], 0);
    }
}
