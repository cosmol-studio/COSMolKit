use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use cosmolkit::{
    Fingerprint, Molecule, PatternFingerprintParams, QueryGraph, SmartsParseParams,
    parse_smarts_with_params, pattern_query_fingerprint_with_params,
};
use serde::Deserialize;
use sha2::{Digest, Sha256};

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
    let (variable, expected_sha) = match profile {
        "pattern_focused" => (
            "COSMOLKIT_PATTERN_ORIGINAL_FOCUSED_GOLDEN",
            "0e5a58bc9cf7ce07836f71fc32dca80b25ae8b3cc038f4e729cf976b8528cd74",
        ),
        "smiles_small" => (
            "COSMOLKIT_PATTERN_ORIGINAL_SMALL_GOLDEN",
            "05908b2765e0d7fa8665f5ac421dd82c8c3ad1f8e8799acde9d5e916c4199eab",
        ),
        "smiles_5000" => (
            "COSMOLKIT_PATTERN_ORIGINAL5000_GOLDEN",
            "b4df21695a5b0caf71b5b4505f43c25f67782c4c72a248c9a7abce7f27767a85",
        ),
        _ => panic!("unknown original Pattern profile {profile}"),
    };
    let path = std::path::PathBuf::from(std::env::var(variable).unwrap_or_else(|_| {
        panic!("prepare original pinned Pattern references and set {variable}")
    }));
    assert_eq!(
        Sha256::digest(std::fs::read(&path).unwrap())
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect::<String>(),
        expected_sha,
        "immutable original reference bytes"
    );
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

enum ParsedInput {
    Concrete(Molecule),
    Query(QueryGraph),
}
impl ParsedInput {
    fn pattern_fingerprint(
        &self,
        p: &PatternFingerprintParams,
    ) -> Result<Fingerprint, cosmolkit::PatternFingerprintError> {
        match self {
            Self::Concrete(m) => m.pattern_fingerprint_with_params(p),
            Self::Query(q) => pattern_query_fingerprint_with_params(q, p),
        }
    }
}
fn parse_record(record: &GoldenRecord) -> Result<ParsedInput, String> {
    match record.input_kind.as_str() {
        "smiles" => Molecule::from_smiles(&record.smiles)
            .map(ParsedInput::Concrete)
            .map_err(|error| error.to_string()),
        "smarts" => parse_smarts_with_params(&record.smiles, &SmartsParseParams::default())
            .map(ParsedInput::Query)
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

    records.iter().enumerate().for_each(|(row, record)| {
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
            assert_eq!(
                actual.n_bits() as usize,
                branch.n_bits.unwrap(),
                "{context}"
            );
            assert_eq!(
                actual
                    .on_bits()
                    .iter()
                    .map(|&b| b as usize)
                    .collect::<Vec<_>>(),
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
