use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use cosmolkit_core::{Molecule, SmilesWriteParams};
use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct BranchResult {
    ok: bool,
    smiles: Option<String>,
    error: Option<String>,
}

#[derive(Debug, Deserialize)]
struct SmilesWriterRecord {
    smiles: String,
    rdkit_ok: bool,
    branches: BTreeMap<String, BranchResult>,
    error: Option<String>,
}

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn load_golden() -> Vec<SmilesWriterRecord> {
    let path = repo_root().join("tests/golden/smiles_writer.jsonl");
    let file = File::open(&path).unwrap_or_else(|err| {
        panic!(
            "failed to open {}; regenerate with tests/scripts/gen_rdkit_smiles_writer_golden.py: {err}",
            path.display()
        )
    });
    BufReader::new(file)
        .lines()
        .enumerate()
        .map(|(idx, line)| {
            let line = line.unwrap_or_else(|err| {
                panic!("failed to read {} line {}: {err}", path.display(), idx + 1)
            });
            serde_json::from_str(&line).unwrap_or_else(|err| {
                panic!("failed to parse {} line {}: {err}", path.display(), idx + 1)
            })
        })
        .collect()
}

fn branch_params(name: &str) -> SmilesWriteParams {
    match name {
        "iso_true_kek_false_can_true" => SmilesWriteParams {
            do_isomeric_smiles: true,
            do_kekule: false,
            canonical: true,
            ..Default::default()
        },
        "iso_false_kek_false_can_true" => SmilesWriteParams {
            do_isomeric_smiles: false,
            do_kekule: false,
            canonical: true,
            ..Default::default()
        },
        "iso_true_kek_false_can_false" => SmilesWriteParams {
            do_isomeric_smiles: true,
            do_kekule: false,
            canonical: false,
            ..Default::default()
        },
        "iso_false_kek_false_can_false" => SmilesWriteParams {
            do_isomeric_smiles: false,
            do_kekule: false,
            canonical: false,
            ..Default::default()
        },
        "iso_true_kek_true_can_true" => SmilesWriteParams {
            do_isomeric_smiles: true,
            do_kekule: true,
            canonical: true,
            ..Default::default()
        },
        "iso_false_kek_true_can_true" => SmilesWriteParams {
            do_isomeric_smiles: false,
            do_kekule: true,
            canonical: true,
            ..Default::default()
        },
        "iso_true_kek_true_can_false" => SmilesWriteParams {
            do_isomeric_smiles: true,
            do_kekule: true,
            canonical: false,
            ..Default::default()
        },
        "iso_false_kek_true_can_false" => SmilesWriteParams {
            do_isomeric_smiles: false,
            do_kekule: true,
            canonical: false,
            ..Default::default()
        },
        other => panic!("unknown smiles writer branch '{other}'"),
    }
}

#[test]
fn smiles_writer_golden_has_one_record_per_smiles() {
    let smiles_path = repo_root().join("tests/smiles.smi");
    let expected = std::fs::read_to_string(&smiles_path)
        .unwrap_or_else(|err| panic!("failed to read {}: {err}", smiles_path.display()))
        .lines()
        .filter(|line| {
            let line = line.trim();
            !line.is_empty() && !line.starts_with('#')
        })
        .count();
    let records = load_golden();
    assert_eq!(
        records.len(),
        expected,
        "smiles writer golden row count must match tests/smiles.smi"
    );
}

#[test]
fn smiles_writer_matches_rdkit_golden_across_param_branches() {
    let records = load_golden();
    for (row_idx, record) in records.iter().enumerate() {
        if !record.rdkit_ok {
            assert!(
                record.error.is_some(),
                "row {} ({}) is rdkit not ok but has no error",
                row_idx + 1,
                record.smiles
            );
            continue;
        }

        let mol = Molecule::from_smiles(&record.smiles).unwrap_or_else(|err| {
            panic!(
                "cosmolkit failed to parse row {} ({}): {err}",
                row_idx + 1,
                record.smiles
            )
        });

        for (branch_name, expected_branch) in &record.branches {
            let params = branch_params(branch_name);
            match mol.to_smiles_with_params(&params) {
                Ok(actual) => {
                    assert!(
                        expected_branch.ok,
                        "row {} ({}) branch {}: cosmolkit succeeded with '{}' but RDKit golden recorded error {:?}",
                        row_idx + 1,
                        record.smiles,
                        branch_name,
                        actual,
                        expected_branch.error
                    );
                    let expected = expected_branch.smiles.as_ref().unwrap_or_else(|| {
                        panic!(
                            "row {} ({}) branch {}: RDKit branch ok but missing smiles",
                            row_idx + 1,
                            record.smiles,
                            branch_name
                        )
                    });
                    assert_eq!(
                        actual,
                        *expected,
                        "smiles writer mismatch at row {} ({}) branch {}",
                        row_idx + 1,
                        record.smiles,
                        branch_name
                    );
                }
                Err(err) => {
                    panic!(
                        "smiles writer unsupported for row {} ({}) branch {} with params {:?}; RDKit golden = {:?}; error = {}",
                        row_idx + 1,
                        record.smiles,
                        branch_name,
                        params,
                        expected_branch.smiles,
                        err
                    );
                }
            }
        }
    }
}
