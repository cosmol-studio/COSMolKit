use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::panic::AssertUnwindSafe;
use std::path::PathBuf;

use cosmolkit_core::{
    BatchErrorMode, Molecule, MoleculeBatch, SmilesWriteError, SmilesWriteParams,
};
use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct BranchResult {
    params: SmilesBranchParams,
    ok: bool,
    smiles: Option<String>,
    error: Option<String>,
}

#[derive(Debug, Deserialize)]
struct SmilesBranchParams {
    do_isomeric_smiles: bool,
    do_kekule: bool,
    canonical: bool,
    clean_stereo: bool,
    all_bonds_explicit: bool,
    all_hs_explicit: bool,
    include_dative_bonds: bool,
    ignore_atom_map_numbers: bool,
    rooted_at_atom: Option<String>,
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
            "failed to open {}; regenerate all RDKit goldens with `.venv/bin/python tests/scripts/gen_all_rdkit_goldens.py --python .venv/bin/python --clean --jobs 4`: {err}",
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

fn branch_params(branch: &BranchResult, mol: &Molecule) -> SmilesWriteParams {
    let rooted_at_atom = match branch.params.rooted_at_atom.as_deref() {
        None => None,
        Some("first") => (!mol.atoms().is_empty()).then_some(0),
        Some("last") => mol.atoms().len().checked_sub(1),
        Some(other) => panic!("unknown rooted_at_atom branch value '{other}'"),
    };
    SmilesWriteParams {
        do_isomeric_smiles: branch.params.do_isomeric_smiles,
        do_kekule: branch.params.do_kekule,
        canonical: branch.params.canonical,
        clean_stereo: branch.params.clean_stereo,
        all_bonds_explicit: branch.params.all_bonds_explicit,
        all_hydrogens_explicit: branch.params.all_hs_explicit,
        include_dative_bonds: branch.params.include_dative_bonds,
        ignore_atom_map_numbers: branch.params.ignore_atom_map_numbers,
        rooted_at_atom,
        ..Default::default()
    }
}

fn common_branch_names() -> &'static [&'static str] {
    &[
        "iso1_kek0_can1_clean1_bond0_hs0_dat1_map0_root_none",
        "iso1_kek0_can0_clean1_bond0_hs0_dat1_map0_root_none",
        "iso1_kek0_can1_clean1_bond0_hs0_dat1_map0_root_last",
        "iso1_kek1_can1_clean1_bond0_hs0_dat1_map0_root_none",
        "iso1_kek1_can0_clean1_bond0_hs0_dat1_map0_root_none",
        "iso1_kek0_can1_clean0_bond0_hs0_dat1_map0_root_none",
        "iso1_kek0_can1_clean1_bond1_hs0_dat1_map0_root_none",
        "iso1_kek0_can1_clean1_bond0_hs1_dat1_map0_root_none",
        "iso1_kek0_can1_clean1_bond0_hs0_dat0_map0_root_none",
        "iso1_kek0_can1_clean1_bond0_hs0_dat1_map1_root_none",
    ]
}

fn assert_branch_names_present(record: &SmilesWriterRecord, branch_names: &[&str]) {
    for branch_name in branch_names {
        assert!(
            record.branches.contains_key(*branch_name),
            "smiles writer golden for '{}' is missing expected branch '{}'",
            record.smiles,
            branch_name
        );
    }
}

fn assert_supported_or_explicitly_unimplemented(
    row_idx: usize,
    smiles: &str,
    branch_name: &str,
    expected_branch: &BranchResult,
    params: &SmilesWriteParams,
    err: &SmilesWriteError,
) {
    assert!(
        !expected_branch.ok || matches!(err, SmilesWriteError::UnsupportedFeature(_)),
        "smiles writer failed without an explicit unsupported error at row {} ({}) branch {} params {:?}; RDKit golden = {:?}; error = {}",
        row_idx + 1,
        smiles,
        branch_name,
        params,
        expected_branch.smiles,
        err
    );
}

fn run_smiles_writer_parity(branch_names: Option<&[&str]>) {
    let records = load_golden();
    let mut total_branches = 0_usize;
    let mut passed_branches = 0_usize;
    let mut unsupported_branches = 0_usize;
    let mut mismatches = Vec::new();

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

        let selected_branches = if let Some(branch_names) = branch_names {
            assert_branch_names_present(record, branch_names);
            branch_names
                .iter()
                .map(|branch_name| {
                    (
                        *branch_name,
                        record.branches.get(*branch_name).unwrap_or_else(|| {
                            panic!(
                                "row {} ({}) missing branch {} after presence check",
                                row_idx + 1,
                                record.smiles,
                                branch_name
                            )
                        }),
                    )
                })
                .collect::<Vec<_>>()
        } else {
            record
                .branches
                .iter()
                .map(|(branch_name, branch)| (branch_name.as_str(), branch))
                .collect::<Vec<_>>()
        };

        let mol = match Molecule::from_smiles(&record.smiles) {
            Ok(mol) => mol,
            Err(err) => {
                total_branches += selected_branches.len();
                mismatches.push(format!(
                    "row {} ({}) parse failed before writer: {}",
                    row_idx + 1,
                    record.smiles,
                    err
                ));
                continue;
            }
        };

        for (branch_name, expected_branch) in selected_branches {
            total_branches += 1;
            let params = branch_params(expected_branch, &mol);
            match mol.to_smiles_with_params(&params) {
                Ok(actual) => {
                    if !expected_branch.ok {
                        mismatches.push(format!(
                            "row {} ({}) branch {}: cosmolkit succeeded with '{}' but RDKit golden recorded error {:?}",
                            row_idx + 1,
                            record.smiles,
                            branch_name,
                            actual,
                            expected_branch.error
                        ));
                        continue;
                    }
                    let expected = expected_branch.smiles.as_ref().unwrap_or_else(|| {
                        panic!(
                            "row {} ({}) branch {}: RDKit branch ok but missing smiles",
                            row_idx + 1,
                            record.smiles,
                            branch_name
                        )
                    });
                    if actual == *expected {
                        passed_branches += 1;
                    } else {
                        mismatches.push(format!(
                            "row {} ({}) branch {} mismatch: expected '{}' got '{}'",
                            row_idx + 1,
                            record.smiles,
                            branch_name,
                            expected,
                            actual
                        ));
                    }
                }
                Err(err) => {
                    if expected_branch.ok {
                        if matches!(err, SmilesWriteError::UnsupportedFeature(_)) {
                            unsupported_branches += 1;
                        }
                        mismatches.push(format!(
                            "row {} ({}) branch {} error with params {:?}: {}",
                            row_idx + 1,
                            record.smiles,
                            branch_name,
                            params,
                            err
                        ));
                    } else {
                        let failed = std::panic::catch_unwind(AssertUnwindSafe(|| {
                            assert_supported_or_explicitly_unimplemented(
                                row_idx,
                                &record.smiles,
                                branch_name,
                                expected_branch,
                                &params,
                                &err,
                            );
                        }));
                        if failed.is_ok() {
                            passed_branches += 1;
                        } else {
                            mismatches.push(format!(
                                "row {} ({}) branch {} error with params {:?}: {}",
                                row_idx + 1,
                                record.smiles,
                                branch_name,
                                params,
                                err
                            ));
                        }
                    }
                }
            }
        }
    }

    let failed_branches = total_branches - passed_branches;
    let mismatch_preview = mismatches
        .iter()
        .take(10)
        .map(String::as_str)
        .collect::<Vec<_>>()
        .join("\n");
    assert!(
        mismatches.is_empty(),
        "smiles writer parity matched {passed_branches}/{total_branches} branches ({:.2}%), failed {failed_branches}, unsupported-on-rdkit-ok {unsupported_branches}\nfirst mismatches:\n{}",
        (passed_branches as f64 / total_branches as f64) * 100.0,
        mismatch_preview
    );
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
fn smiles_writer_matches_rdkit_golden_for_common_param_branches() {
    run_smiles_writer_parity(Some(common_branch_names()));
}

#[test]
fn smiles_writer_matches_rdkit_golden_for_common_root_none_param_branches_in_parallel_batch() {
    let records = load_golden();
    let smiles = records
        .iter()
        .map(|record| record.smiles.clone())
        .collect::<Vec<_>>();
    let batch = MoleculeBatch::from_smiles_list(&smiles).with_parallel_jobs(Some(4));
    let template_mol = Molecule::from_smiles("CC").expect("template molecule should parse");

    for branch_name in common_branch_names() {
        let first_branch = records
            .iter()
            .find_map(|record| record.branches.get(*branch_name))
            .unwrap_or_else(|| panic!("missing common branch {branch_name}"));
        if first_branch.params.rooted_at_atom.is_some() {
            continue;
        }
        let params = branch_params(first_branch, &template_mol);
        let actual = batch
            .to_smiles_list_with_params_and_options(
                &params,
                BatchErrorMode::KeepErrors,
                Some(4),
                Some(false),
            )
            .unwrap_or_else(|err| {
                panic!("parallel batch SMILES writer failed for branch {branch_name}: {err}")
            });

        assert_eq!(actual.len(), records.len());
        for (row_idx, (record, actual_smiles)) in records.iter().zip(actual.iter()).enumerate() {
            if !record.rdkit_ok {
                assert!(
                    record.error.is_some(),
                    "row {} ({}) is rdkit not ok but has no error",
                    row_idx + 1,
                    record.smiles
                );
                continue;
            }
            let expected_branch = record.branches.get(*branch_name).unwrap_or_else(|| {
                panic!(
                    "row {} ({}) missing common branch {}",
                    row_idx + 1,
                    record.smiles,
                    branch_name
                )
            });
            if !expected_branch.ok {
                assert_eq!(
                    actual_smiles,
                    "?",
                    "parallel batch branch {} should keep error marker at row {} ({})",
                    branch_name,
                    row_idx + 1,
                    record.smiles
                );
                continue;
            }
            let expected = expected_branch.smiles.as_ref().unwrap_or_else(|| {
                panic!(
                    "row {} ({}) branch {}: RDKit branch ok but missing smiles",
                    row_idx + 1,
                    record.smiles,
                    branch_name
                )
            });
            assert_eq!(
                actual_smiles,
                expected,
                "parallel batch SMILES mismatch at row {} ({}) branch {}",
                row_idx + 1,
                record.smiles,
                branch_name
            );
        }
    }
}

#[test]
#[ignore = "full RDKit smiles writer branch matrix is expensive; run explicitly for exhaustive parity"]
fn smiles_writer_matches_rdkit_golden_across_param_branches() {
    run_smiles_writer_parity(None);
}
