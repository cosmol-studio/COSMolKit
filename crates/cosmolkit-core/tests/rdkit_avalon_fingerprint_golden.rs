use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use cosmolkit_core::properties::avalon_fingerprint::avalon_fingerprint;
use cosmolkit_core::{AvalonFingerprintFlags, AvalonFingerprintParams, Molecule};
use rayon::prelude::*;
use serde::Deserialize;

mod common;
use common::parity_data;

const GOLDEN_PATH: &str = "tmp/parity-audit/avalon_fingerprint_smiles_5000.jsonl";
const FOCUSED_PATH: &str = "tmp/parity-audit/avalon_fingerprint_focused.json";

#[derive(Debug, Deserialize)]
struct CorpusRecord {
    row: usize,
    smiles: String,
    branches: BTreeMap<String, BranchRecord>,
}

#[derive(Debug, Deserialize)]
struct BranchRecord {
    smiles: String,
    parameters: Option<Parameters>,
    ok: bool,
    n_bits: Option<u32>,
    is_query: Option<bool>,
    bit_flags: Option<String>,
    on_bits: Option<Vec<usize>>,
    error: Option<String>,
}

#[derive(Debug, Clone, Deserialize)]
struct Parameters {
    name: String,
    n_bits: u32,
    is_query: bool,
    bit_flags: String,
    smiles: Option<String>,
}

#[derive(Debug, Deserialize)]
struct FocusedFile {
    cases: Vec<FocusedCase>,
}

#[derive(Debug, Deserialize)]
struct FocusedCase {
    case: Parameters,
    result: BranchRecord,
}

fn read_jsonl() -> Vec<CorpusRecord> {
    let path = parity_data::repo_root().join(GOLDEN_PATH);
    let file = File::open(&path)
        .unwrap_or_else(|error| panic!("failed to open {}: {error}", path.display()));
    BufReader::new(file)
        .lines()
        .enumerate()
        .map(|(line, content)| {
            let content = content.unwrap_or_else(|error| {
                panic!(
                    "failed to read {} line {}: {error}",
                    path.display(),
                    line + 1
                )
            });
            serde_json::from_str(&content).unwrap_or_else(|error| {
                panic!(
                    "failed to parse {} line {}: {error}",
                    path.display(),
                    line + 1
                )
            })
        })
        .collect()
}

fn flags(value: &str) -> AvalonFingerprintFlags {
    let bits = u32::from_str_radix(value.trim_start_matches("0x"), 16)
        .unwrap_or_else(|error| panic!("invalid Avalon flags {value}: {error}"));
    AvalonFingerprintFlags::from_bits(bits)
        .unwrap_or_else(|| panic!("undefined Avalon flags {value}"))
}

fn assert_branch(smiles: &str, branch_name: &str, branch: &BranchRecord, molecule: &Molecule) {
    if !branch.ok {
        assert!(
            avalon_fingerprint(molecule, &AvalonFingerprintParams::default()).is_err(),
            "row {smiles} branch {branch_name} unexpectedly succeeded; source error: {:?}",
            branch.error
        );
        return;
    }
    let parameters = branch
        .parameters
        .as_ref()
        .unwrap_or_else(|| panic!("{smiles} {branch_name} is missing source parameters"));
    let params = AvalonFingerprintParams {
        n_bits: branch.n_bits.unwrap_or(parameters.n_bits),
        is_query: branch.is_query.unwrap_or(parameters.is_query),
        bit_flags: flags(branch.bit_flags.as_deref().unwrap_or(&parameters.bit_flags)),
    };
    assert_eq!(parameters.name, branch_name, "branch metadata changed");
    let actual = avalon_fingerprint(molecule, &params);
    let actual =
        actual.unwrap_or_else(|error| panic!("row {smiles} branch {branch_name} failed: {error}"));
    assert_eq!(
        actual.n_bits(),
        branch.n_bits.unwrap() as usize,
        "{smiles} {branch_name} size mismatch"
    );
    assert_eq!(
        actual.on_bits(),
        branch.on_bits.clone().unwrap_or_default(),
        "{smiles} {branch_name} exact bits mismatch"
    );
}

#[test]
fn avalon_fingerprint_matches_every_active_5000_row_profile_exactly() {
    let records = read_jsonl();
    assert_eq!(
        records.len(),
        5000,
        "Avalon corpus golden must contain every row"
    );
    let mut failures: Vec<(usize, String)> = records
        .par_iter()
        .enumerate()
        .filter_map(|(row, record)| {
            let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                assert_eq!(record.row, row, "Avalon row index changed");
                let molecule = Molecule::from_smiles(&record.smiles);
                for (branch_name, branch) in &record.branches {
                    if branch.ok {
                        let molecule = molecule.as_ref().unwrap_or_else(|error| {
                            panic!("row {row} ({}) failed to parse: {error}", record.smiles)
                        });
                        assert_branch(&record.smiles, branch_name, branch, molecule);
                    } else if let Ok(molecule) = &molecule {
                        assert_branch(&record.smiles, branch_name, branch, molecule);
                    }
                }
            }));
            result.err().map(|payload| {
                let message = payload
                    .downcast_ref::<String>()
                    .cloned()
                    .or_else(|| {
                        payload
                            .downcast_ref::<&str>()
                            .map(|value| (*value).to_owned())
                    })
                    .unwrap_or_else(|| "non-string panic".to_owned());
                (row, message)
            })
        })
        .collect();
    failures.sort_by_key(|(row, _)| *row);
    assert!(
        failures.is_empty(),
        "Avalon exact parity failures: {failures:?}"
    );
}

#[test]
fn avalon_fingerprint_matches_focused_profiles_and_size_boundaries() {
    let path = parity_data::repo_root().join(FOCUSED_PATH);
    let file = File::open(&path)
        .unwrap_or_else(|error| panic!("failed to open {}: {error}", path.display()));
    let focused: FocusedFile = serde_json::from_reader(file).expect("focused Avalon JSON is valid");
    assert_eq!(focused.cases.len(), 8);
    for case in focused.cases {
        let mut result = case.result;
        result.parameters = Some(case.case.clone());
        assert_branch(
            &result.smiles,
            &case.case.name,
            &result,
            &Molecule::from_smiles(&result.smiles).unwrap_or_else(|error| {
                panic!("focused {} failed to parse: {error}", case.case.name)
            }),
        );
    }
}
