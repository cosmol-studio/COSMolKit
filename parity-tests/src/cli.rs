use crate::Result;
use std::path::PathBuf;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Command {
    Prepare,
    Compare,
    List,
}

#[derive(Debug, Clone)]
pub struct Options {
    pub command: Command,
    pub corpus: PathBuf,
    pub suites: Vec<String>,
    pub tasks: Vec<String>,
    pub threads: usize,
    pub offset: usize,
    pub limit: Option<usize>,
    pub profiles: Option<PathBuf>,
    pub random_profiles: usize,
    pub seed: u64,
}

impl Options {
    pub fn parse(args: impl IntoIterator<Item = String>) -> Result<Self> {
        let mut options = Self {
            command: Command::Compare,
            corpus: PathBuf::from("smiles_small"),
            suites: vec!["full".into()],
            tasks: Vec::new(),
            threads: std::thread::available_parallelism()?.get(),
            offset: 0,
            limit: None,
            profiles: None,
            random_profiles: 0,
            seed: 20260306,
        };
        let mut args = args.into_iter();
        while let Some(arg) = args.next() {
            match arg.as_str() {
                "prepare" => options.command = Command::Prepare,
                "compare" => options.command = Command::Compare,
                "list" | "--list" => options.command = Command::List,
                "--nocapture" | "--" => (),
                "--corpus" | "--suite" | "--task" | "--threads" | "--offset" | "--limit"
                | "--profiles" | "--random-profiles" | "--seed" => {
                    let value = args.next().ok_or_else(|| format!("{arg} needs a value"))?;
                    match arg.as_str() {
                        "--corpus" => options.corpus = value.into(),
                        "--suite" => options.suites = names(&value)?,
                        "--task" => options.tasks.extend(names(&value)?),
                        "--threads" => options.threads = value.parse()?,
                        "--offset" => options.offset = value.parse()?,
                        "--limit" => options.limit = Some(value.parse()?),
                        "--profiles" => options.profiles = Some(value.into()),
                        "--random-profiles" => options.random_profiles = value.parse()?,
                        "--seed" => options.seed = value.parse()?,
                        _ => unreachable!(),
                    }
                }
                _ if !arg.starts_with('-') => options.tasks.push(arg),
                _ => return Err(format!("unknown option {arg}").into()),
            }
        }
        if options.threads == 0 || options.limit == Some(0) {
            return Err("--threads and --limit must be positive".into());
        }
        Ok(options)
    }
}

fn names(value: &str) -> Result<Vec<String>> {
    let names: Vec<_> = value.split(',').map(str::to_owned).collect();
    if names.iter().any(String::is_empty) {
        return Err("empty task/suite name".into());
    }
    Ok(names)
}

/// Exact names, rather than silent substring matches or unknown-name omissions.
pub fn select<'a>(available: &'a [(&'a str, &'a str)], options: &Options) -> Result<Vec<&'a str>> {
    for name in &options.tasks {
        if !available.iter().any(|(task, _)| task == name) {
            return Err(format!("unknown task {name}; use list").into());
        }
    }
    for name in &options.suites {
        if name != "full" && !available.iter().any(|(_, suite)| suite == name) {
            return Err(format!("unknown suite {name}; use list").into());
        }
    }
    let selected: Vec<_> = available
        .iter()
        .filter(|(task, suite)| {
            (options.suites.iter().any(|s| s == "full" || s == suite))
                && (options.tasks.is_empty() || options.tasks.iter().any(|t| t == task))
        })
        .map(|(task, _)| *task)
        .collect();
    if selected.is_empty() {
        return Err("selection contains no tasks".into());
    }
    Ok(selected)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn parse(args: &[&str]) -> Result<Options> {
        Options::parse(args.iter().map(|arg| (*arg).into()))
    }
    const TASKS: &[(&str, &str)] = &[("cip", "chemistry"), ("morgan", "fingerprints")];

    #[test]
    fn full_single_and_sets_are_runtime_selected() {
        assert_eq!(
            select(TASKS, &parse(&[]).unwrap()).unwrap(),
            ["cip", "morgan"]
        );
        assert_eq!(
            select(TASKS, &parse(&["--task", "cip"]).unwrap()).unwrap(),
            ["cip"]
        );
        assert_eq!(
            select(TASKS, &parse(&["--suite", "fingerprints"]).unwrap()).unwrap(),
            ["morgan"]
        );
        assert_eq!(
            select(
                TASKS,
                &parse(&["--suite", "chemistry,fingerprints"]).unwrap()
            )
            .unwrap()
            .len(),
            2
        );
    }

    #[test]
    fn invalid_or_empty_selections_never_pass() {
        assert!(select(TASKS, &parse(&["--task", "typo"]).unwrap()).is_err());
        assert!(
            select(
                TASKS,
                &parse(&["--suite", "fingerprints", "--task", "cip"]).unwrap()
            )
            .is_err()
        );
        assert!(parse(&["--threads", "0"]).is_err());
        assert!(parse(&["--limit", "0"]).is_err());
        assert!(parse(&["--unknown"]).is_err());
        assert!(parse(&["--task", "cip,"]).is_err());
    }

    #[test]
    fn worker_count_and_corpus_never_require_recompilation() {
        let options = parse(&[
            "prepare",
            "--threads",
            "112",
            "--corpus",
            "other.smi",
            "--offset",
            "20",
            "--limit",
            "100",
        ])
        .unwrap();
        assert_eq!(options.threads, 112);
        assert_eq!(options.command, Command::Prepare);
        assert_eq!(options.corpus, PathBuf::from("other.smi"));
        assert_eq!((options.offset, options.limit), (20, Some(100)));
    }
}
