use cosmolkit_parity_tests::{
    self as parity, CorpusSource,
    registry::{self, CorpusType},
};
use std::path::PathBuf;

fn main() {
    if let Err(error) = entry() {
        eprintln!("{error}");
        std::process::exit(1);
    }
}

fn entry() -> Result<(), String> {
    let mut args = std::env::args().skip(1);
    let command = args.next().unwrap_or_else(|| "help".into());
    if command == "help" || command == "--help" {
        println!(
            "prepare | preflight | list\n  --task FUNCTION_CORPUS (default all corpus tasks)\n  --special-regression KEY (fixed source matrix; separate from corpus tasks)\n  --smiles FILE (default smiles_small.smi)\n  --fingerprint-pairs FILE (default deterministic 5000 pairs)\n  --data DIR (default target/parity-tests)\n  --python PATH (prepare only; default .venv/bin/python)\n  --threads N (prepare only; default 4)\nAfter corpus prepare: cargo test -p cosmolkit-parity-tests --release --test reference_parity\nAfter special-regression prepare: cargo test -p cosmolkit-parity-tests --release --test special_regression_structure_tags\nTests never generate reference values."
        );
        return Ok(());
    }
    if !["prepare", "preflight", "list"].contains(&command.as_str()) {
        return Err(format!(
            "unknown command: {command}; comparison uses cargo test"
        ));
    }
    let mut task = None;
    let mut special = None;
    let mut sources = Vec::new();
    let mut data = parity::root().join("target/parity-tests");
    let mut python = parity::root().join(".venv/bin/python");
    let mut threads = 4;
    let mut seen = std::collections::BTreeSet::new();
    while let Some(flag) = args.next() {
        if !seen.insert(flag.clone()) {
            return Err(format!("duplicate option: {flag}"));
        }
        let value = args
            .next()
            .ok_or_else(|| format!("missing value: {flag}"))?;
        match flag.as_str() {
            "--task" => task = Some(value),
            "--special-regression" => special = Some(value),
            "--smiles" => sources.push(CorpusSource {
                corpus_type: CorpusType::Smiles,
                path: value.into(),
            }),
            "--fingerprint-pairs" => sources.push(CorpusSource {
                corpus_type: CorpusType::FingerprintPairs,
                path: value.into(),
            }),
            "--data" => data = value.into(),
            "--python" if command == "prepare" => python = PathBuf::from(value),
            "--threads" if command == "prepare" => {
                threads = value
                    .parse()
                    .map_err(|_| "threads must be a positive integer")?;
                if threads == 0 {
                    return Err("threads must be positive".into());
                }
            }
            _ => return Err(format!("unsupported option: {flag}")),
        }
    }
    if let Some(key) = special {
        if task.is_some() || !sources.is_empty() || seen.contains("--threads") {
            return Err("special regressions use their fixed inputs; do not combine with corpus/task/thread selection".into());
        }
        if command == "list" {
            let selected = registry::SPECIAL_REGRESSIONS
                .iter()
                .find(|row| row.key == key)
                .ok_or_else(|| format!("unknown special regression: {key}"))?;
            println!(
                "{}: special_regression; {} fixed cases",
                selected.key, selected.rows
            );
        } else if command == "prepare" {
            let result = parity::special_regression::prepare(&key, &data, &python)?;
            println!(
                "Special preparation complete: {} reused, {} generated; {} fixed cases",
                result.reused_tasks, result.generated_tasks, result.rows
            );
        } else {
            let result = parity::special_regression::preflight(&key, &data)?;
            println!(
                "Ready: {} fixed cases; 0 Rust operation calls",
                result.rows.len()
            );
        }
        return Ok(());
    }
    // With no new selection, validate exactly the persisted prepared suite.
    if command == "preflight" && task.is_none() && sources.is_empty() {
        let (_, ready) = parity::load_suite(&data)?;
        println!("Ready: {} cases; 0 Rust operation calls", ready.len());
        return Ok(());
    }
    let tasks = if task.as_deref() == Some("descriptors_existing") {
        registry::TASKS
            .iter()
            .filter(|task| parity::descriptor_reference::handles(task))
            .collect()
    } else if task.as_deref() == Some("tautomers_existing") {
        registry::TASKS
            .iter()
            .filter(|task| parity::tautomer_reference::handles(task))
            .collect()
    } else {
        registry::select(task.as_deref())?
    };
    let cases = parity::corpus(&sources, &tasks)?;
    registry::validate(&cases, &tasks)?;
    for task in &tasks {
        println!(
            "{}: {} cases; {}",
            task.key(),
            task.count(&cases),
            task.generator
        );
    }
    match command.as_str() {
        "list" => {}
        "prepare" => {
            let mut preparation = parity::prepare(&tasks, &cases, &data, &python, threads)?;
            let mut special_keys = Vec::new();
            // Default preparation is complete: newly registered special recipes
            // are included automatically, without a workflow-specific task name.
            if task.is_none() {
                for row in registry::SPECIAL_REGRESSIONS {
                    let prepared = parity::special_regression::prepare(row.key, &data, &python)?;
                    preparation.reused_tasks += prepared.reused_tasks;
                    preparation.generated_tasks += prepared.generated_tasks;
                    preparation.rows += prepared.rows;
                    special_keys.push(row.key.to_owned());
                }
            }
            parity::special_regression::preflight_selection(&special_keys, &data)?;
            parity::save_suite_with_special(&data, &tasks, &sources, special_keys)?;
            println!(
                "Preparation complete: {} reused tasks, {} generated tasks; {} ready cases",
                preparation.reused_tasks, preparation.generated_tasks, preparation.rows
            );
        }
        "preflight" => println!(
            "Ready: {} cases; 0 Rust operation calls",
            parity::preflight(&tasks, &cases, &data)?.len()
        ),
        _ => unreachable!(),
    }
    Ok(())
}
