use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use cosmolkit_core::{BatchErrorMode, Molecule, MoleculeBatch, SmilesWriteParams};
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
        all_hs_explicit: branch.params.all_hs_explicit,
        include_dative_bonds: branch.params.include_dative_bonds,
        ignore_atom_map_numbers: branch.params.ignore_atom_map_numbers,
        rooted_at_atom,
        ..Default::default()
    }
}

fn assert_supported_or_explicitly_unimplemented(
    row_idx: usize,
    smiles: &str,
    branch_name: &str,
    expected_branch: &BranchResult,
    params: &SmilesWriteParams,
    err: &cosmolkit_core::SmilesWriteError,
) {
    assert!(
        !expected_branch.ok || matches!(err, cosmolkit_core::SmilesWriteError::UnsupportedPath(_)),
        "smiles writer failed without an explicit unsupported error at row {} ({}) branch {} params {:?}; RDKit golden = {:?}; error = {}",
        row_idx + 1,
        smiles,
        branch_name,
        params,
        expected_branch.smiles,
        err
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
            let params = branch_params(expected_branch, &mol);
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
                    assert_supported_or_explicitly_unimplemented(
                        row_idx,
                        &record.smiles,
                        branch_name,
                        expected_branch,
                        &params,
                        &err,
                    );
                }
            }
        }
    }

    let smiles = records
        .iter()
        .map(|record| record.smiles.clone())
        .collect::<Vec<_>>();
    let batch = MoleculeBatch::from_smiles_list(&smiles, BatchErrorMode::Keep)
        .expect("batch SMILES parse should not raise in keep mode");
    for branch_name in records
        .iter()
        .flat_map(|record| record.branches.keys())
        .map(String::as_str)
        .collect::<std::collections::BTreeSet<_>>()
    {
        let first_ok_record = records
            .iter()
            .find(|record| record.rdkit_ok)
            .expect("SMILES writer golden has no RDKit-ok rows");
        if first_ok_record.branches[branch_name]
            .params
            .rooted_at_atom
            .is_some()
        {
            continue;
        }
        let first_ok_mol = Molecule::from_smiles(&first_ok_record.smiles)
            .expect("first RDKit-ok SMILES should parse");
        let params = branch_params(&first_ok_record.branches[branch_name], &first_ok_mol);
        let Ok(actual) = batch.to_smiles_list_with_params(&params) else {
            continue;
        };
        for (row_idx, (record, actual)) in records.iter().zip(actual).enumerate() {
            if !record.rdkit_ok {
                continue;
            }
            let expected = record.branches[branch_name]
                .smiles
                .as_ref()
                .expect("RDKit ok branch should have SMILES");
            assert_eq!(
                actual.as_ref(),
                Some(expected),
                "batch smiles writer mismatch at row {} ({}) branch {}",
                row_idx + 1,
                record.smiles,
                branch_name
            );
        }
    }
}
