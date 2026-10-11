use crate::{
    Result,
    cli::{Command, Options, select},
    corpus::Corpus,
    storage::{Identity, References},
    tasks, value,
};
use serde_json::json;
use sha2::{Digest, Sha256};
use std::{
    fs::{self, File},
    io::{BufWriter, Write},
    path::Path,
    process::{Command as Process, Stdio},
    time::Instant,
};

pub fn main() -> Result<()> {
    if std::env::args().nth(1).as_deref() == Some("structure-case") {
        return crate::structure_tags::child();
    }
    let options = Options::parse(std::env::args().skip(1))?;
    let available = tasks::available();
    if options.command == Command::List {
        for (name, suite) in &available {
            println!("{name:36} {suite}");
        }
        return Ok(());
    }
    let started = Instant::now();
    let selected = select(&available, &options)?;
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let overrides: serde_json::Value = if let Some(path) = &options.profiles {
        serde_json::from_reader(File::open(profile_path(root, path))?)?
    } else {
        json!({})
    };
    let mut selection = Vec::new();
    let mut inputs = std::collections::HashMap::new();
    for task in selected {
        let mut input_options = options.clone();
        if let Some(path) = tasks::fixed_input(task) {
            if !(options.tasks.len() == 1
                && options.corpus.extension().is_some_and(|e| e == "jsonl"))
            {
                input_options.corpus = root.join(path);
                input_options.offset = 0;
                input_options.limit = None;
            }
        }
        let input_key = (
            input_options.corpus.clone(),
            input_options.offset,
            input_options.limit,
        );
        if !inputs.contains_key(&input_key) {
            inputs.insert(input_key.clone(), Corpus::open(&input_options)?);
        }
        let corpus = inputs[&input_key].clone();
        let profiles = match overrides.get(task) {
            Some(value) => value
                .as_array()
                .filter(|array| !array.is_empty())
                .ok_or("profiles must be a nonempty array")?
                .clone(),
            None => tasks::profiles(task)?,
        };
        let mut profiles = profiles
            .into_iter()
            .map(|profile| tasks::normalize_profile(task, profile))
            .collect::<Result<Vec<_>>>()?;
        profiles.extend(crate::random_profiles::sample(
            task,
            options.random_profiles,
            options.seed,
        )?);
        let identity = Identity {
            rdkit: tasks::RDKIT.into(),
            task: task.into(),
            profiles: json!(profiles),
            input_sha256: corpus.sha256.clone(),
            records: corpus.records,
        };
        let folder = root.join("expected/v1").join(identity.key()?);
        selection.push((identity, folder, profiles, corpus));
    }
    if options.command == Command::Prepare {
        let mut missing = Vec::new();
        for (identity, folder, _, corpus) in &selection {
            match References::open(folder, identity) {
                Ok(_) => println!("{}: cached", identity.task),
                Err(error) => {
                    println!("{}: prepare ({error})", identity.task);
                    missing.push(json!({"identity": identity, "folder": folder,
                        "corpus":corpus.path,"offset":corpus.offset,"limit":corpus.limit}));
                }
            }
        }
        if !missing.is_empty() {
            let mut process = Process::new(root.parent().unwrap().join(".venv/bin/python"))
                .arg(root.join("tools/prepare.py"))
                .stdin(Stdio::piped())
                .spawn()?;
            let request = json!({"rdkit": tasks::RDKIT, "threads": options.threads,
                "tasks": missing});
            serde_json::to_writer(
                process.stdin.take().ok_or("reference stdin unavailable")?,
                &request,
            )?;
            let status = process.wait()?;
            if !status.success() {
                return Err(format!("native preparation failed: {status}").into());
            }
        }
        // Header/identity only; successful payloads are never loaded at prepare exit.
        for (identity, folder, _, _) in &selection {
            References::open(folder, identity)?;
        }
        println!(
            "Prepared {} tasks in {:.3}s",
            selection.len(),
            started.elapsed().as_secs_f64()
        );
        return Ok(());
    }
    let mut references = Vec::new();
    // Missing selection members fail before any CK call. Open handles, not copies.
    for (identity, folder, _, _) in &selection {
        references.push(References::open(folder, identity).map_err(|error| {
            format!(
                "{}: {error}; run prepare with the same corpus/task/profile arguments",
                identity.task
            )
        })?);
    }
    println!(
        "Input identity + {} reference headers: {:.3}s; workers={}",
        selection.len(),
        started.elapsed().as_secs_f64(),
        options.threads
    );
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(options.threads)
        .build()?;
    let reports = root.join("reports");
    fs::create_dir_all(&reports)?;
    let report_dir = tempfile::Builder::new()
        .prefix("run-")
        .tempdir_in(&reports)?
        .keep();
    let mut all_failures = 0u64;
    let mut empty_tasks = 0u64;
    for ((identity, _, profiles, corpus), mut reference) in selection.into_iter().zip(references) {
        let clock = Instant::now();
        let prepared_search = matches!(identity.task.as_str(), "substructure" | "tautomer")
            .then(|| crate::search_tautomer::prepare(&identity.task, &profiles))
            .transpose()?;
        let mut rows = corpus.rows()?;
        let mut completed = 0usize;
        let mut input_digest = Sha256::new();
        input_digest.update((corpus.offset as u64).to_le_bytes());
        let mut matches = 0u64;
        let mut errors = 0u64;
        let mut failures = 0u64;
        let mut undefined = 0u64;
        let mut timeouts = 0u64;
        let mut ck_timeouts = 0u64;
        let mut timeout_cases = 0u64;
        let mut crashes = 0u64;
        let mut report = BufWriter::new(File::create(
            report_dir.join(format!("{}.jsonl", identity.task)),
        )?);
        loop {
            // Storage frames and compute windows are independent: even >64 workers
            // receive enough work, without loading a whole corpus/reference.
            let mut expected_rows = Vec::new();
            let mut expected_bytes = 0usize;
            while expected_rows.len() < options.threads.saturating_mul(4).max(64)
                && expected_bytes < crate::storage::MAX_FRAME_BYTES
            {
                let Some(frame) = reference.next_frame()? else {
                    break;
                };
                expected_bytes += frame.iter().map(Vec::len).sum::<usize>();
                expected_rows.extend(frame);
            }
            if expected_rows.is_empty() {
                break;
            }
            let inputs: Vec<_> = (&mut rows)
                .take(expected_rows.len())
                .collect::<Result<_>>()?;
            if inputs.len() != expected_rows.len() {
                return Err("corpus changed after preflight".into());
            }
            for input in &inputs {
                input_digest.update((input.len() as u64).to_le_bytes());
                input_digest.update(input.as_bytes());
            }
            let actual = if let Some(prepared) = &prepared_search {
                prepared.evaluate(&inputs, &pool)
            } else {
                tasks::evaluate(
                    &identity.task,
                    &profiles,
                    &inputs,
                    completed,
                    options.threads,
                    &pool,
                )?
            };
            if actual.len() != expected_rows.len()
                || actual.iter().any(|row| row.len() != profiles.len())
            {
                return Err(
                    "actual result dimensions do not match selected inputs/profiles".into(),
                );
            }
            for (index, ((input, expected), actual)) in
                inputs.iter().zip(expected_rows).zip(actual).enumerate()
            {
                let expected = tasks::unpack_row(&expected, profiles.len())?;
                for (profile, (expected, actual)) in expected.into_iter().zip(actual).enumerate() {
                    let status = expected[0];
                    let matching_error = tasks::matching_error(&identity.task, expected, &actual);
                    let (reference_timeout, actual_timeout) =
                        timeout_flags(&identity.task, expected, &actual);
                    if reference_timeout || actual_timeout {
                        timeout_cases += 1;
                        timeouts += u64::from(reference_timeout);
                        ck_timeouts += u64::from(actual_timeout);
                    } else {
                        match status {
                            0 | 1 if expected == actual => {
                                if status == 0 {
                                    matches += 1;
                                } else {
                                    errors += 1;
                                }
                            }
                            2 => undefined += 1,
                            4 => crashes += 1,
                            5 if matching_error => errors += 1,
                            0 | 1 | 5 => failures += 1,
                            _ => return Err("unknown reference outcome status".into()),
                        }
                    }
                    if reference_timeout || actual_timeout || !(status <= 1 && expected == actual) {
                        let offset = expected
                            .iter()
                            .zip(&actual)
                            .position(|(a, b)| a != b)
                            .unwrap_or(expected.len().min(actual.len()));
                        let diagnostic = |bytes: &[u8]| -> Result<String> {
                            if bytes.len() <= 8192 {
                                Ok(format!("{:?}", value::decode(&bytes[1..])?))
                            } else {
                                let start = offset.saturating_sub(8).min(bytes.len());
                                let end = start.saturating_add(32).min(bytes.len());
                                Ok(format!(
                                    "{} bytes; sha256={}; bytes[{start}..{end}]={}",
                                    bytes.len(),
                                    crate::storage::hex(&Sha256::digest(bytes)),
                                    crate::storage::hex(&bytes[start..end])
                                ))
                            }
                        };
                        serde_json::to_writer(
                            &mut report,
                            &json!({"index": completed + index,
                            "input": input, "profile": profiles[profile], "reference_status": status,
                            "matching_source_error": matching_error,
                            "reference_timeout":reference_timeout, "actual_timeout":actual_timeout,
                            "actual_status": actual[0], "first_different_byte":offset,
                            "expected":diagnostic(expected)?, "actual":diagnostic(&actual)?}),
                        )?;
                        writeln!(report)?;
                    }
                }
            }
            completed += inputs.len();
            eprint!(
                "\r{} [{completed}/{}] comparisons={} mismatches={failures}",
                identity.task,
                identity.records,
                matches + errors + failures
            );
        }
        eprintln!();
        if rows.next().is_some() || completed as u64 != corpus.records {
            return Err("corpus/reference count mismatch".into());
        }
        if crate::storage::hex(&input_digest.finalize()) != identity.input_sha256 {
            return Err("corpus changed during comparison".into());
        }
        report.flush()?;
        let no_comparable_results = matches + errors + failures == 0;
        // A timeout-only task is incomplete, not a failed or passed comparison.
        empty_tasks += u64::from(no_comparable_results && timeout_cases == 0);
        let summary = json!({"task": identity.task, "records": completed, "profiles": profiles.len(),
            "matched_values": matches, "matched_errors": errors, "mismatches": failures,
            "no_comparable_results": no_comparable_results,
            "upstream_undefined": undefined, "upstream_timeouts": timeouts, "upstream_crashes": crashes,
            "ck_timeouts":ck_timeouts, "timeout_cases":timeout_cases,
            "seconds": clock.elapsed().as_secs_f64(), "threads": options.threads});
        serde_json::to_writer_pretty(
            File::create(report_dir.join(format!("{}.summary.json", identity.task)))?,
            &summary,
        )?;
        println!("{summary}");
        all_failures += failures;
    }
    println!("Reports: {}", report_dir.display());
    if all_failures != 0 || empty_tasks != 0 {
        return Err(format!(
            "{all_failures} mismatches/failures; {empty_tasks} tasks without comparable results"
        )
        .into());
    }
    Ok(())
}

fn timeout_flags(task: &str, expected: &[u8], actual: &[u8]) -> (bool, bool) {
    let timed_out = |bytes: &[u8]| {
        bytes.first() == Some(&3) || (task == "embedding" && embedding_timed_out(bytes))
    };
    (timed_out(expected), timed_out(actual))
}

fn embedding_timed_out(bytes: &[u8]) -> bool {
    // The embedding adapter's exact format is [ids, coordinates, failure_counts].
    // Inspect its fixed-width integer trailer, never decode/copy all coordinates.
    // RDKit Embedder.h: EXCEEDED_TIMEOUT = 11, END_OF_ENUM = 15.
    let count = cosmolkit::EmbedFailureCause::EndOfEnum.rdkit_ordinal() as usize;
    let slot = cosmolkit::EmbedFailureCause::ExceededTimeout.rdkit_ordinal() as usize;
    if !bytes.starts_with(&[0, 7, 3, 0, 0, 0]) || bytes.len() < 5 + count * 9 {
        return false;
    }
    let trailer = &bytes[bytes.len() - (5 + count * 9)..];
    if trailer[0] != 7
        || trailer[1..5] != (count as u32).to_le_bytes()
        || trailer[5..].chunks_exact(9).any(|field| field[0] != 2)
    {
        return false;
    }
    let start = 5 + slot * 9 + 1;
    i64::from_le_bytes(trailer[start..start + 8].try_into().unwrap()) > 0
}

fn profile_path(root: &Path, path: &Path) -> std::path::PathBuf {
    if path.is_absolute() || path.exists() {
        path.to_owned()
    } else {
        root.parent().unwrap().join(path)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn timeouts_are_separate_even_when_both_results_match() {
        use value::Value;
        let result = |slot: usize| {
            let mut counts = vec![Value::Int(0); 15];
            counts[slot] = Value::Int(1);
            tasks::success(Value::List(vec![
                Value::List(vec![Value::Int(-1)]),
                Value::List(vec![]),
                Value::List(counts),
            ]))
        };
        let timeout = result(11);
        let other_failure = result(5);
        assert_eq!(timeout_flags("embedding", &timeout, &timeout), (true, true));
        assert_eq!(
            timeout_flags("embedding", &timeout, &other_failure),
            (true, false)
        );
        assert_eq!(
            timeout_flags("embedding", &other_failure, &timeout),
            (false, true)
        );
        assert_eq!(
            timeout_flags("embedding", &other_failure, &other_failure),
            (false, false)
        );
        assert_eq!(timeout_flags("valence", &timeout, &timeout), (false, false));
        let native_limit = tasks::outcome(3, Value::Text("native worker timed out".into()));
        assert_eq!(
            timeout_flags("smiles_write", &native_limit, &other_failure),
            (true, false)
        );
        assert!(!embedding_timed_out(&timeout[..timeout.len() - 1]));
    }

    #[test]
    fn cargo_profile_path_resolves_from_repository_root() {
        let repository = tempfile::tempdir().unwrap();
        let root = repository.path().join("parity-tests");
        let relative = Path::new("target/parity-test-unique-profiles.json");
        assert_eq!(
            profile_path(&root, relative),
            repository.path().join(relative)
        );
        let absolute = repository.path().join("profiles.json");
        assert_eq!(profile_path(&root, &absolute), absolute);
    }
}
