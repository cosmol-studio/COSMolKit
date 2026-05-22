use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use cosmolkit_core::{BatchRecord, Molecule, MoleculeBatch};
use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct DgBoundsRecord {
    smiles: String,
    rdkit_ok: bool,
    bounds: Option<Vec<Vec<f64>>>,
    error: Option<String>,
}

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn load_golden() -> Vec<DgBoundsRecord> {
    let path = repo_root().join("tests/golden/dg_bounds_matrix.jsonl");
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

#[test]
fn dg_bounds_golden_has_one_record_per_smiles() {
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
        "DG bounds golden row count must match tests/smiles.smi"
    );
}

#[test]
fn dg_bounds_matrix_matches_rdkit_golden() {
    let records = load_golden();
    let row_filter = std::env::var("COSMOLKIT_DG_ROW_FILTER")
        .ok()
        .and_then(|s| s.parse::<usize>().ok());
    for (row_idx, record) in records.iter().enumerate() {
        if let Some(filter) = row_filter {
            if row_idx + 1 != filter {
                continue;
            }
        }
        if !record.rdkit_ok {
            assert!(
                record.error.is_some(),
                "row {} ({}) is rdkit not ok but has no error",
                row_idx + 1,
                record.smiles
            );
            continue;
        }
        let expected = record.bounds.as_ref().unwrap_or_else(|| {
            panic!(
                "row {} ({}) is rdkit_ok but has no bounds",
                row_idx + 1,
                record.smiles
            )
        });
        let mol = Molecule::from_smiles(&record.smiles).unwrap_or_else(|err| {
            panic!(
                "cosmolkit failed to parse row {} ({}): {err}",
                row_idx + 1,
                record.smiles
            )
        });
        let actual = mol.dg_bounds_matrix().unwrap_or_else(|err| {
            panic!(
                "DG bounds matrix is not implemented for row {} ({}), expected RDKit golden comparison: {err}",
                row_idx + 1,
                record.smiles
            )
        });

        assert_eq!(
            actual.len(),
            expected.len(),
            "DG bounds row count mismatch at row {} ({})",
            row_idx + 1,
            record.smiles
        );
        for (i, (actual_row, expected_row)) in actual.iter().zip(expected).enumerate() {
            assert_eq!(
                actual_row.len(),
                expected_row.len(),
                "DG bounds column count mismatch at row {} ({}) matrix row {}",
                row_idx + 1,
                record.smiles,
                i
            );
            for (j, (&a, &e)) in actual_row.iter().zip(expected_row).enumerate() {
                assert!(
                    (a - e).abs() <= 1e-8,
                    "DG bounds mismatch at row {} ({}) matrix[{}][{}]: ours={} expected={}",
                    row_idx + 1,
                    record.smiles,
                    i,
                    j,
                    a,
                    e
                );
            }
        }
    }

    if row_filter.is_none() {
        let smiles = records
            .iter()
            .map(|record| record.smiles.clone())
            .collect::<Vec<_>>();
        let batch = MoleculeBatch::from_smiles_list(&smiles);
        for (row_idx, (record, batch_record)) in records.iter().zip(batch.iter()).enumerate() {
            if !record.rdkit_ok {
                continue;
            }
            let molecule = match batch_record {
                BatchRecord::Molecule(molecule) => molecule,
                BatchRecord::Error(error) => {
                    panic!(
                        "batch SMILES parse failed at row {} ({}): {}",
                        row_idx + 1,
                        record.smiles,
                        error.message
                    )
                }
            };
            let actual = molecule.dg_bounds_matrix().unwrap_or_else(|err| {
                panic!(
                    "batch DG bounds missing at row {} ({}): {err}",
                    row_idx + 1,
                    record.smiles
                )
            });
            let expected = record.bounds.as_ref().expect("RDKit ok row has bounds");
            for (i, (actual_row, expected_row)) in actual.iter().zip(expected).enumerate() {
                for (j, (&a, &e)) in actual_row.iter().zip(expected_row).enumerate() {
                    assert!(
                        (a - e).abs() <= 1e-8,
                        "batch DG bounds mismatch at row {} ({}) matrix[{}][{}]: ours={} expected={}",
                        row_idx + 1,
                        record.smiles,
                        i,
                        j,
                        a,
                        e
                    );
                }
            }
        }
    }
}

#[test]
fn dg_bounds_matrix_matches_rdkit_golden_in_parallel_batch() {
    let records = load_golden();
    let smiles = records
        .iter()
        .map(|record| record.smiles.clone())
        .collect::<Vec<_>>();
    let batch = MoleculeBatch::from_smiles_list(&smiles).with_parallel_jobs(Some(4));
    let actual = batch
        .dg_bounds_matrix_list_with_progress(None)
        .expect("parallel batch DG bounds should collect without molecule errors");

    assert_eq!(actual.len(), records.len());
    for (row_idx, (record, actual)) in records.iter().zip(actual.iter()).enumerate() {
        if !record.rdkit_ok {
            assert!(
                record.error.is_some(),
                "row {} ({}) is rdkit not ok but has no error",
                row_idx + 1,
                record.smiles
            );
            assert!(
                actual.is_none(),
                "parallel batch DG bounds should not produce row {} ({})",
                row_idx + 1,
                record.smiles
            );
            continue;
        }
        let actual = actual.as_ref().unwrap_or_else(|| {
            panic!(
                "parallel batch DG bounds missing at row {} ({})",
                row_idx + 1,
                record.smiles
            )
        });
        let expected = record.bounds.as_ref().expect("RDKit ok row has bounds");
        assert_eq!(
            actual.len(),
            expected.len(),
            "parallel batch DG bounds row count mismatch at row {} ({})",
            row_idx + 1,
            record.smiles
        );
        for (i, (actual_row, expected_row)) in actual.iter().zip(expected).enumerate() {
            assert_eq!(
                actual_row.len(),
                expected_row.len(),
                "parallel batch DG bounds column count mismatch at row {} ({}) matrix row {}",
                row_idx + 1,
                record.smiles,
                i
            );
            for (j, (&a, &e)) in actual_row.iter().zip(expected_row).enumerate() {
                assert!(
                    (a - e).abs() <= 1e-8,
                    "parallel batch DG bounds mismatch at row {} ({}) matrix[{}][{}]: ours={} expected={}",
                    row_idx + 1,
                    record.smiles,
                    i,
                    j,
                    a,
                    e
                );
            }
        }
    }
}
