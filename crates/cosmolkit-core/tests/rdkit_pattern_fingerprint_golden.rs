use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use cosmolkit_core::{Molecule, PatternFingerprintParams, SmartsParseParams, mol_from_smarts};
use rayon::prelude::*;
use serde::Deserialize;

mod common;
use common::parity_data;

const OUTPUT_NAME: &str = "pattern_fingerprint.jsonl";

#[derive(Debug, Deserialize)]
struct GoldenRecord {
    row: usize,
    smiles: String,
    input_kind: String,
    rdkit_ok: bool,
    branches: BTreeMap<String, GoldenBranch>,
    error: Option<String>,
}

#[derive(Debug, Deserialize)]
struct GoldenBranch {
    parameters: GoldenParameters,
    ok: bool,
    n_bits: Option<usize>,
    on_bits: Option<Vec<usize>>,
    atom_counts_before: Option<Vec<u32>>,
    atom_counts_after: Option<Vec<u32>>,
    set_only_bits: Option<Vec<usize>>,
    error: Option<String>,
}

#[derive(Debug, Deserialize)]
struct GoldenParameters {
    name: String,
    #[serde(rename = "fpSize")]
    fp_size: usize,
    #[serde(rename = "tautomericFingerprint")]
    tautomeric_fingerprint: bool,
    #[serde(rename = "atomCounts")]
    atom_counts: String,
    #[serde(rename = "setOnlyBits")]
    set_only_bits: String,
}

fn read_records(profile: &str) -> Vec<GoldenRecord> {
    let path = parity_data::expected_path_for_profile("fingerprint", "rdkit", profile, OUTPUT_NAME);
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

fn parse_record(record: &GoldenRecord) -> Result<Molecule, String> {
    match record.input_kind.as_str() {
        "smiles" => Molecule::from_smiles(&record.smiles).map_err(|error| error.to_string()),
        "smarts" => mol_from_smarts(&record.smiles, &SmartsParseParams::default())
            .map_err(|error| error.to_string()),
        other => Err(format!("unknown Pattern input kind {other:?}")),
    }
}

fn assert_profile(profile: &str, expected_records: usize, permits_query_inputs: bool) {
    let records = read_records(profile);
    assert_eq!(
        records.len(),
        expected_records,
        "Pattern profile {profile} row count changed"
    );
    if permits_query_inputs {
        assert!(
            records.iter().any(|record| record.input_kind == "smarts"),
            "focused Pattern profile must exercise query-bearing graphs"
        );
    } else {
        assert!(
            records.iter().all(|record| record.input_kind == "smiles"),
            "corpus Pattern profiles must retain every SMILES row"
        );
    }

    records.par_iter().enumerate().for_each(|(row, record)| {
        assert_eq!(record.row, row, "Pattern row identity changed in {profile}");
        let molecule = parse_record(record);
        if !record.rdkit_ok {
            assert!(
                molecule.is_err(),
                "{profile} row {row} ({}) parsed only in COSMolKit; RDKit error: {:?}",
                record.smiles,
                record.error
            );
            assert!(record.branches.is_empty());
            return;
        }
        assert!(record.error.is_none(), "RDKit-success row has an error");
        let molecule = molecule.unwrap_or_else(|error| {
            panic!(
                "{profile} row {row} ({}) failed to parse in COSMolKit: {error}",
                record.smiles
            )
        });
        assert_eq!(
            record.branches.len(),
            11,
            "{profile} row {row} Pattern branch matrix changed"
        );

        for (branch_name, branch) in &record.branches {
            let context = format!(
                "{profile} row {row} ({}) branch {branch_name}",
                record.smiles
            );
            assert_eq!(branch.parameters.name, *branch_name, "{context}");
            if branch.parameters.atom_counts == "none" {
                assert!(branch.atom_counts_before.is_none(), "{context}");
                assert!(branch.atom_counts_after.is_none(), "{context}");
            } else {
                assert_eq!(
                    branch.atom_counts_after, branch.atom_counts_before,
                    "{context}: pinned ordinary overload must leave atomCounts inert"
                );
            }

            if branch.parameters.set_only_bits == "wrong_width" {
                assert!(!branch.ok, "{context}: invalid mask unexpectedly succeeded");
                assert!(
                    branch
                        .error
                        .as_deref()
                        .is_some_and(|error| error.contains("bad setOnlyBits size")),
                    "{context}: source mask validation error changed: {:?}",
                    branch.error
                );
                continue;
            }

            assert!(branch.ok, "{context}: RDKit error: {:?}", branch.error);
            if branch.parameters.set_only_bits == "none" {
                assert!(branch.set_only_bits.is_none(), "{context}");
            } else {
                assert!(branch.set_only_bits.is_some(), "{context}");
            }
            let params = PatternFingerprintParams {
                n_bits: branch.parameters.fp_size,
                tautomeric: branch.parameters.tautomeric_fingerprint,
            };
            let actual = molecule
                .pattern_fingerprint(&params)
                .unwrap_or_else(|error| panic!("{context}: COSMolKit fingerprint failed: {error}"));
            assert_eq!(actual.n_bits(), branch.n_bits.unwrap(), "{context}");
            assert_eq!(
                actual.on_bits(),
                branch.on_bits.as_deref().unwrap(),
                "{context}: exact ordered on-bit set differs"
            );
        }
    });
}

#[test]
fn pattern_fingerprint_matches_every_focused_profile_row_exactly() {
    assert_profile("pattern_focused", 18, true);
}

#[test]
fn pattern_fingerprint_matches_every_smiles_small_profile_row_exactly() {
    assert_profile("smiles_small", 152, false);
}

#[test]
fn pattern_fingerprint_matches_every_smiles_5000_profile_row_exactly() {
    assert_profile("smiles_5000", 5_000, false);
}
