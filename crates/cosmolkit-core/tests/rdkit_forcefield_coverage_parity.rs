use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use cosmolkit_core::{
    MmffMolProperties, Molecule, mmff_has_all_molecule_params, uff_has_all_molecule_params,
};
use serde::Deserialize;

mod common;
use common::parity_data;

const CHARGE_TOLERANCE: f64 = 1.0e-12;

#[derive(Debug, Deserialize)]
struct ForcefieldResult {
    ok: bool,
    has_all: Option<bool>,
    #[serde(default)]
    atom_types: Option<Vec<u8>>,
    #[serde(default)]
    formal_charges: Option<Vec<f64>>,
    #[serde(default)]
    partial_charges: Option<Vec<f64>>,
    error: Option<String>,
}

#[derive(Debug, Deserialize)]
struct ForcefieldCoverageRecord {
    smiles: String,
    rdkit_ok: bool,
    uff: ForcefieldResult,
    mmff: ForcefieldResult,
    uff_explicit_h: ForcefieldResult,
    mmff_explicit_h: ForcefieldResult,
    error: Option<String>,
}

fn load_golden() -> Vec<ForcefieldCoverageRecord> {
    let path = std::env::var_os("COSMOLKIT_FORCEFIELD_COVERAGE_GOLDEN")
        .map(PathBuf::from)
        .unwrap_or_else(|| parity_data::golden_path("forcefield_coverage.jsonl"));
    let file = File::open(&path).unwrap_or_else(|err| {
        panic!(
            "failed to open {}; prepare it with `{}`: {err}",
            path.display(),
            parity_data::rdkit_prepare_command("forcefield_coverage")
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

fn assert_uff_coverage(
    row: usize,
    smiles: &str,
    mol: &Molecule,
    expected: &ForcefieldResult,
    surface: &str,
) {
    assert!(
        expected.ok,
        "row {row} ({smiles}) has RDKit {surface} UFF error: {:?}",
        expected.error
    );
    let expected_has_all = expected.has_all.unwrap_or_else(|| {
        panic!("row {row} ({smiles}) has RDKit {surface} UFF result without has_all")
    });
    let actual_has_all = uff_has_all_molecule_params(mol).unwrap_or_else(|err| {
        panic!("row {row} ({smiles}) COSMolKit {surface} UFF coverage errored: {err}")
    });
    assert_eq!(
        actual_has_all, expected_has_all,
        "row {row} ({smiles}) {surface} UFF parameter coverage mismatch"
    );
}

fn assert_mmff_coverage(
    row: usize,
    smiles: &str,
    mol: &Molecule,
    expected: &ForcefieldResult,
    surface: &str,
) {
    assert!(
        expected.ok,
        "row {row} ({smiles}) has RDKit {surface} MMFF error: {:?}",
        expected.error
    );
    let expected_has_all = expected.has_all.unwrap_or_else(|| {
        panic!("row {row} ({smiles}) has RDKit {surface} MMFF result without has_all")
    });
    let actual_has_all = mmff_has_all_molecule_params(mol).unwrap_or_else(|err| {
        panic!("row {row} ({smiles}) COSMolKit {surface} MMFF coverage errored: {err}")
    });
    assert_eq!(
        actual_has_all, expected_has_all,
        "row {row} ({smiles}) {surface} MMFF parameter coverage mismatch"
    );

    let Some(expected_atom_types) = expected.atom_types.as_ref() else {
        assert!(
            !expected_has_all,
            "row {row} ({smiles}) has RDKit {surface} MMFF parameters without atom types"
        );
        assert!(expected.formal_charges.is_none());
        assert!(expected.partial_charges.is_none());
        return;
    };

    let props = MmffMolProperties::new(mol, "MMFF94", 0).unwrap_or_else(|err| {
        panic!("row {row} ({smiles}) COSMolKit {surface} MMFF properties errored: {err}")
    });
    let actual_atom_types = (0..mol.num_atoms())
        .map(|idx| {
            props.get_mmff_atom_type(idx).unwrap_or_else(|err| {
                panic!(
                    "row {row} ({smiles}) COSMolKit {surface} MMFF atom {idx} type errored: {err}"
                )
            })
        })
        .collect::<Vec<_>>();
    assert_eq!(
        actual_atom_types, *expected_atom_types,
        "row {row} ({smiles}) {surface} MMFF atom type mismatch"
    );

    let expected_formal = expected.formal_charges.as_ref().unwrap_or_else(|| {
        panic!("row {row} ({smiles}) RDKit {surface} MMFF result has no formal charges")
    });
    let expected_partial = expected.partial_charges.as_ref().unwrap_or_else(|| {
        panic!("row {row} ({smiles}) RDKit {surface} MMFF result has no partial charges")
    });
    assert_eq!(expected_formal.len(), mol.num_atoms());
    assert_eq!(expected_partial.len(), mol.num_atoms());
    for idx in 0..mol.num_atoms() {
        let actual_formal = props.get_mmff_formal_charge(idx).unwrap();
        let actual_partial = props.get_mmff_partial_charge(idx).unwrap();
        assert!(
            (actual_formal - expected_formal[idx]).abs() <= CHARGE_TOLERANCE,
            "row {row} ({smiles}) {surface} MMFF formal charge mismatch at atom {idx}: actual={actual_formal} expected={}",
            expected_formal[idx]
        );
        assert!(
            (actual_partial - expected_partial[idx]).abs() <= CHARGE_TOLERANCE,
            "row {row} ({smiles}) {surface} MMFF partial charge mismatch at atom {idx}: actual={actual_partial} expected={}",
            expected_partial[idx]
        );
    }
}

#[test]
fn forcefield_coverage_matches_rdkit_for_every_active_profile_row() {
    let records = load_golden();
    let start_row = std::env::var("COSMOLKIT_FORCEFIELD_COVERAGE_START_ROW")
        .ok()
        .map(|value| {
            value.parse::<usize>().unwrap_or_else(|err| {
                panic!("COSMOLKIT_FORCEFIELD_COVERAGE_START_ROW must be a positive integer: {err}")
            })
        })
        .unwrap_or(1);
    let shard_count = std::env::var("COSMOLKIT_FORCEFIELD_COVERAGE_SHARD_COUNT")
        .ok()
        .map(|value| {
            value.parse::<usize>().unwrap_or_else(|err| {
                panic!(
                    "COSMOLKIT_FORCEFIELD_COVERAGE_SHARD_COUNT must be a positive integer: {err}"
                )
            })
        })
        .unwrap_or(1);
    let shard_index = std::env::var("COSMOLKIT_FORCEFIELD_COVERAGE_SHARD_INDEX")
        .ok()
        .map(|value| {
            value.parse::<usize>().unwrap_or_else(|err| {
                panic!("COSMOLKIT_FORCEFIELD_COVERAGE_SHARD_INDEX must be a non-negative integer: {err}")
            })
        })
        .unwrap_or(0);
    assert!(start_row > 0, "coverage start row must be one-based");
    assert!(shard_count > 0, "coverage shard count must be positive");
    assert!(
        shard_index < shard_count,
        "coverage shard index must be less than shard count"
    );
    assert_eq!(
        records.len(),
        parity_data::count_smiles_rows(),
        "force-field coverage golden row count must match the active corpus"
    );

    for (row_idx, record) in records.iter().enumerate() {
        let row = row_idx + 1;
        if row < start_row || row_idx % shard_count != shard_index {
            continue;
        }
        if !record.rdkit_ok {
            assert!(
                record.error.is_some(),
                "row {row} ({}) is RDKit-not-ok without an error",
                record.smiles
            );
            continue;
        }

        let mol = Molecule::from_smiles(&record.smiles).unwrap_or_else(|err| {
            panic!(
                "row {row} ({}) COSMolKit parse errored: {err}",
                record.smiles
            )
        });
        assert_uff_coverage(row, &record.smiles, &mol, &record.uff, "implicit-H");
        assert_mmff_coverage(row, &record.smiles, &mol, &record.mmff, "implicit-H");

        let explicit_h_mol = mol.with_hydrogens().unwrap_or_else(|err| {
            panic!(
                "row {row} ({}) COSMolKit explicit-H construction errored: {err}",
                record.smiles
            )
        });
        assert_uff_coverage(
            row,
            &record.smiles,
            &explicit_h_mol,
            &record.uff_explicit_h,
            "explicit-H",
        );
        assert_mmff_coverage(
            row,
            &record.smiles,
            &explicit_h_mol,
            &record.mmff_explicit_h,
            "explicit-H",
        );
    }
}
