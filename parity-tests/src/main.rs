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
            "prepare | preflight | list\n  --task FUNCTION_CORPUS (default all)\n  --smiles FILE (default smiles_small.smi)\n  --fingerprint-pairs FILE (default deterministic 5000 pairs)\n  --data DIR (default target/parity-tests)\n  --python PATH (prepare only; default .venv/bin/python)\n  --threads N (prepare only; default 4)\nAfter prepare: cargo test -p cosmolkit-parity-tests --release --test reference_parity\nTests never generate reference values."
        );
        return Ok(());
    }
    if !["prepare", "preflight", "list"].contains(&command.as_str()) {
        return Err(format!(
            "unknown command: {command}; comparison uses cargo test"
        ));
    }
    let mut task = None;
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
    // With no new selection, validate exactly the persisted prepared suite.
    if command == "preflight" && task.is_none() && sources.is_empty() {
        let (_, ready) = parity::load_suite(&data)?;
        println!("Ready: {} cases; 0 Rust operation calls", ready.len());
        return Ok(());
    }
    let tasks = registry::select(task.as_deref())?;
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
            let preparation = parity::prepare(&tasks, &cases, &data, &python, threads)?;
            parity::save_suite(&data, &tasks, &sources)?;
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
