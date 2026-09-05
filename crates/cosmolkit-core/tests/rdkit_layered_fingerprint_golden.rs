use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use cosmolkit_core::{Fingerprint, LayeredFingerprintLayers, LayeredFingerprintParams, Molecule};
use rayon::prelude::*;
use serde::Deserialize;

mod common;
use common::parity_data;

const OUTPUT_NAME: &str = "layered_fingerprint.jsonl";
const PROFILE_BRANCH_COUNT: usize = 18;

#[derive(Debug, Deserialize)]
struct GoldenRecord {
    row: usize,
    smiles: String,
    rdkit_ok: bool,
    branches: BTreeMap<String, GoldenBranch>,
    error: Option<String>,
}

#[derive(Debug, Deserialize)]
struct GoldenBranch {
    parameters: GoldenParameters,
    resolved_arguments: Option<GoldenArguments>,
    ok: bool,
    num_bits: Option<usize>,
    on_bits: Option<Vec<usize>>,
    atom_counts: Option<Vec<u32>>,
    error: Option<String>,
}

#[derive(Debug, Deserialize)]
struct GoldenParameters {
    name: String,
    #[serde(rename = "layerFlags")]
    layer_flags: String,
    #[serde(rename = "minPath")]
    min_path: u32,
    #[serde(rename = "maxPath")]
    max_path: u32,
    #[serde(rename = "fpSize")]
    fp_size: u32,
    #[serde(rename = "branchedPaths")]
    branched_paths: bool,
    #[serde(default, rename = "fromAtoms")]
    from_atoms: Option<String>,
    #[serde(default, rename = "atomCounts")]
    atom_counts: Option<String>,
    #[serde(default, rename = "setOnlyBits")]
    set_only_bits: Option<String>,
}

#[derive(Debug, Deserialize)]
struct GoldenArguments {
    #[serde(rename = "layerFlags")]
    layer_flags: u32,
    #[serde(rename = "fromAtoms")]
    from_atoms: Option<Vec<u32>>,
    #[serde(rename = "atomCounts")]
    atom_counts: Option<Vec<u32>>,
    #[serde(rename = "setOnlyBits")]
    set_only_bits: Option<Vec<usize>>,
}

fn read_records(profile: &str) -> Vec<GoldenRecord> {
    let path = parity_data::expected_path_for_profile("fingerprint", "rdkit", profile, OUTPUT_NAME);
    let file = File::open(&path).unwrap_or_else(|error| panic!("failed to open {}: {error}", path.display()));
    BufReader::new(file)
        .lines()
        .enumerate()
        .map(|(line, content)| {
            let content =
                content.unwrap_or_else(|error| panic!("failed to read {} line {}: {error}", path.display(), line + 1));
            serde_json::from_str(&content)
                .unwrap_or_else(|error| panic!("failed to parse {} line {}: {error}", path.display(), line + 1))
        })
        .collect()
}

fn read_smiles(path: &Path) -> Vec<String> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|error| panic!("failed to read {}: {error}", path.display()))
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(str::to_owned)
        .collect()
}

fn parse_source_u32(value: &str) -> u32 {
    value
        .strip_prefix("0x")
        .and_then(|digits| u32::from_str_radix(digits, 16).ok())
        .unwrap_or_else(|| panic!("invalid source unsigned value {value:?}"))
}

fn params_from_golden(parameters: &GoldenParameters, arguments: &GoldenArguments) -> LayeredFingerprintParams {
    let layer_flags = parse_source_u32(&parameters.layer_flags);
    assert_eq!(arguments.layer_flags, layer_flags, "resolved layer flags");
    assert_eq!(parameters.from_atoms.is_some(), arguments.from_atoms.is_some());
    assert_eq!(parameters.atom_counts.is_some(), arguments.atom_counts.is_some());
    assert_eq!(parameters.set_only_bits.is_some(), arguments.set_only_bits.is_some());
    LayeredFingerprintParams {
        layers: LayeredFingerprintLayers::from_bits_retain(layer_flags),
        min_path: parameters.min_path,
        max_path: parameters.max_path,
        fp_size: parameters.fp_size,
        atom_counts: arguments.atom_counts.clone(),
        set_only_bits: arguments
            .set_only_bits
            .as_ref()
            .map(|bits| Fingerprint::from_on_bits(parameters.fp_size as usize, bits.iter().copied())),
        branched_paths: parameters.branched_paths,
        from_atoms: arguments.from_atoms.clone(),
    }
}

fn assert_branch(
    profile: &str,
    row: usize,
    smiles: &str,
    branch_name: &str,
    branch: &GoldenBranch,
    molecule: &Molecule,
) {
    let context = format!("{profile} row {row} ({smiles}) branch {branch_name}");
    assert_eq!(branch.parameters.name, branch_name, "{context}: identity");
    let arguments = branch
        .resolved_arguments
        .as_ref()
        .unwrap_or_else(|| panic!("{context}: missing resolved arguments"));
    let actual = molecule.layered_fingerprint_with_output(&params_from_golden(&branch.parameters, arguments));
    if !branch.ok {
        assert!(
            actual.is_err(),
            "{context}: unexpectedly succeeded; RDKit error: {:?}",
            branch.error
        );
        return;
    }
    assert!(branch.error.is_none(), "{context}: success has an error");
    let actual = actual.unwrap_or_else(|error| panic!("{context}: {error}"));
    assert_eq!(
        actual.fingerprint.n_bits(),
        branch.num_bits.unwrap_or_else(|| panic!("{context}: missing num_bits")),
        "{context}: vector size"
    );
    let actual_bits = actual.fingerprint.on_bits();
    let expected_bits = branch
        .on_bits
        .as_ref()
        .unwrap_or_else(|| panic!("{context}: missing on_bits"));
    assert_eq!(
        actual_bits.as_slice(),
        expected_bits.as_slice(),
        "{context}: exact bits"
    );
    assert_eq!(
        actual.atom_counts.as_ref(),
        branch.atom_counts.as_ref(),
        "{context}: exact atom counts"
    );
}

fn assert_profile(profile: &str, corpus_path: &Path) {
    let records = read_records(profile);
    let corpus = read_smiles(corpus_path);
    assert_eq!(records.len(), corpus.len(), "{profile}: record count");
    let failures: Vec<(usize, String)> = records
        .par_iter()
        .zip(corpus.par_iter())
        .enumerate()
        .filter_map(|(row, (record, expected_smiles))| {
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                assert_eq!(record.row, row, "{profile}: row identity");
                assert_eq!(record.smiles, *expected_smiles, "{profile}: input identity");
                if !record.rdkit_ok {
                    assert!(Molecule::from_smiles(expected_smiles).is_err());
                    assert!(record.branches.is_empty());
                    assert!(record.error.is_some());
                    return;
                }
                assert!(record.error.is_none(), "{profile} row {row}: parse error");
                assert_eq!(
                    record.branches.len(),
                    PROFILE_BRANCH_COUNT,
                    "{profile} row {row}: incomplete branch matrix"
                );
                let molecule = Molecule::from_smiles(expected_smiles)
                    .unwrap_or_else(|error| panic!("{profile} row {row} ({expected_smiles}) parse failed: {error}"));
                for (branch_name, branch) in &record.branches {
                    assert_branch(profile, row, expected_smiles, branch_name, branch, &molecule);
                }
            }))
            .err()
            .map(|payload| {
                let message = payload
                    .downcast_ref::<String>()
                    .cloned()
                    .or_else(|| payload.downcast_ref::<&str>().map(|value| (*value).to_owned()))
                    .unwrap_or_else(|| "non-string panic".to_owned());
                (row, message)
            })
        })
        .collect();
    assert!(failures.is_empty(), "{profile} Layered parity failures: {failures:?}");
}

#[test]
fn layered_fingerprint_matches_complete_focused_and_small_profiles_exactly() {
    assert_profile(
        "layered_focused",
        &parity_data::repo_root().join("testdata/fingerprint/fixtures/rdkit/layered_fingerprint_focused.smi"),
    );
    assert_profile(
        "smiles_small",
        &parity_data::repo_root().join("testdata/smiles/corpus/smiles_small.smi"),
    );
}

#[test]
fn layered_fingerprint_matches_complete_5000_profile_exactly() {
    assert_profile(
        "smiles_5000",
        &parity_data::repo_root().join("testdata/smiles/corpus/smiles_5000.smi"),
    );
}
