use cosmolkit_parity_tests_fixed::{Selection, prepare_with_reuse};

fn main() {
    if let Err(error) = entry() {
        eprintln!("{error}");
        std::process::exit(1);
    }
}
fn entry() -> Result<(), String> {
    let mut args = std::env::args().skip(1);
    let command = args.next().unwrap_or_else(|| "--help".into());
    if matches!(command.as_str(), "--help" | "help") {
        println!(
            "prepare --corpus NAME | --special all|NAME [--threads N] [--task NAME] [--reuse-from CHECKOUT]\nCompare with cargo test --test corpus or --test special_regression."
        );
        return Ok(());
    }
    if command != "prepare" {
        return Err("only prepare is a CLI command; comparison uses cargo test".into());
    }
    let mut selection = None;
    let mut task = None;
    let mut reuse_from = None;
    let mut threads = std::thread::available_parallelism().map_or(4, usize::from);
    let mut seen = std::collections::BTreeSet::new();
    while let Some(flag) = args.next() {
        if !seen.insert(flag.clone()) {
            return Err(format!("duplicate option: {flag}"));
        }
        let value = args
            .next()
            .ok_or_else(|| format!("missing value: {flag}"))?;
        match flag.as_str() {
            "--corpus" | "--special" => {
                if selection.is_some() {
                    return Err("select either --corpus or --special".into());
                }
                selection = Some(if flag == "--corpus" {
                    Selection::Corpus(value)
                } else {
                    Selection::Special(value)
                });
            }
            "--task" => task = Some(value),
            "--reuse-from" => reuse_from = Some(std::path::PathBuf::from(value)),
            "--threads" => {
                threads = value
                    .parse::<usize>()
                    .ok()
                    .filter(|n| *n > 0)
                    .ok_or("threads must be positive")?
            }
            _ => return Err(format!("unknown option: {flag}")),
        }
    }
    let selected = selection.ok_or("select --corpus NAME or --special all|NAME")?;
    let result = prepare_with_reuse(&selected, threads, task.as_deref(), reuse_from.as_deref())?;
    println!(
        "Prepared: {} generated, {} reused, {} rows",
        result.generated_tasks, result.reused_tasks, result.rows
    );
    Ok(())
}
