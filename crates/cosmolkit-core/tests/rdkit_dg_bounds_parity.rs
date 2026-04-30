use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use cosmolkit_core::Molecule;
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
            "failed to open {}; regenerate with tests/scripts/gen_rdkit_dg_bounds_golden.py: {err}",
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
    for (row_idx, record) in records.iter().enumerate() {
        if let Some(filter) = std::env::var("COSMOLKIT_DG_ROW_FILTER")
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
        {
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
}
