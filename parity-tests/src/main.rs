use cosmolkit_parity_tests::{self as pilot, registry};
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
            "run | prepare | preflight | list\n  run automatically prepares missing/invalid references before global preflight\n  --task fuzzy_and|fuzzy_or (default: both)\n  --corpus FILE (default: builtin branch corpus)\n  --data DIR (default: target/parity-tests)\n  --python PATH (run/prepare; default: .venv/bin/python)\nNo Python/JS binding or performance verification in this pilot."
        );
        return Ok(());
    }
    if !["prepare", "preflight", "run", "list"].contains(&command.as_str()) {
        return Err(format!("unknown command: {command}"));
    }
    let mut task = None;
    let mut corpus = None;
    let mut data = pilot::root().join("target/parity-tests");
    let mut python = pilot::root().join(".venv/bin/python");
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
            "--corpus" => corpus = Some(PathBuf::from(value)),
            "--data" => data = PathBuf::from(value),
            "--python" if command == "prepare" || command == "run" => python = PathBuf::from(value),
            _ => return Err(format!("unsupported option: {flag}")),
        }
    }
    let tasks = registry::select(task.as_deref())?;
    let cases = pilot::corpus(corpus.as_deref())?;
    registry::validate(&cases, &tasks)?;
    for t in &tasks {
        println!(
            "{}: {} inputs x {:?} = {} cases",
            t.operation.name(),
            cases.len(),
            t.widths,
            cases.len() * t.widths.len()
        );
    }
    match command.as_str() {
        "list" => {}
        "prepare" => println!(
            "Preparation: {:?}",
            pilot::prepare(&tasks, &cases, &data, &python)?
        ),
        "preflight" => println!(
            "Ready: {} cases; 0 Rust operation calls",
            pilot::preflight(&tasks, &cases, &data)?.len()
        ),
        "run" => {
            let report = pilot::run(&tasks, &cases, &data, &python, pilot::execute::run)?;
            let failed = report.iter().filter(|row| !row.matches).count();
            let output = data.join("rust-report.json");
            pilot::write_report(&output, &report)?;
            println!(
                "{} compared, {} failed; {}",
                report.len(),
                failed,
                output.display()
            );
            if failed != 0 {
                return Err("parity mismatch; see typed expected/actual results".into());
            }
        }
        _ => unreachable!(),
    }
    Ok(())
}
