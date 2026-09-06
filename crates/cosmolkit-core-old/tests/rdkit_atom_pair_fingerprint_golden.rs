use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use cosmolkit_core::{
    AdditionalOutput, AtomPairFingerprintGenerator, AtomPairFingerprintParams,
    FingerprintFuncArguments, Molecule,
};
use rayon::prelude::*;
use serde::Deserialize;

mod common;
use common::parity_data;

const OUTPUT_NAME: &str = "atom_pair_fingerprint.jsonl";

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
    resolved_arguments: GoldenArguments,
    sparse_count: GoldenOutput,
    count: GoldenOutput,
    sparse_bit: GoldenOutput,
    explicit_bit: GoldenOutput,
    ok: bool,
    error: Option<String>,
}

#[derive(Debug, Deserialize)]
struct GoldenParameters {
    name: String,
    #[serde(rename = "minDistance")]
    min_distance: u32,
    #[serde(rename = "maxDistance")]
    max_distance: u32,
    #[serde(rename = "includeChirality")]
    include_chirality: bool,
    #[serde(rename = "use2D")]
    use_2d: bool,
    #[serde(rename = "countSimulation")]
    count_simulation: bool,
    #[serde(rename = "countBounds")]
    count_bounds: Vec<u32>,
    #[serde(rename = "fpSize")]
    fp_size: usize,
    #[serde(rename = "numBitsPerFeature")]
    num_bits_per_feature: u32,
}

#[derive(Debug, Default, Deserialize)]
struct GoldenArguments {
    #[serde(default, rename = "fromAtoms")]
    from_atoms: Option<Vec<usize>>,
    #[serde(default, rename = "ignoreAtoms")]
    ignore_atoms: Option<Vec<usize>>,
    #[serde(default, rename = "customAtomInvariants")]
    custom_atom_invariants: Option<Vec<u32>>,
}

#[derive(Debug, Deserialize)]
struct GoldenOutput {
    ok: bool,
    length: Option<u64>,
    #[serde(default)]
    nonzero_elements: Option<BTreeMap<u64, i32>>,
    #[serde(default)]
    on_bits: Option<Vec<u64>>,
    #[serde(default)]
    additional_output: Option<GoldenAdditionalOutput>,
    error: Option<String>,
}

#[derive(Debug, Deserialize)]
struct GoldenAdditionalOutput {
    atom_counts: Vec<u32>,
    atom_to_bits: Vec<Vec<u64>>,
    bit_info_map: BTreeMap<u64, Vec<(u32, u32)>>,
    atoms_per_bit: BTreeMap<u64, Vec<Vec<usize>>>,
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

fn read_smiles(path: &Path) -> Vec<String> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|error| panic!("failed to read {}: {error}", path.display()))
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(str::to_owned)
        .collect()
}

fn params(parameters: &GoldenParameters) -> AtomPairFingerprintParams {
    AtomPairFingerprintParams {
        n_bits: parameters.fp_size,
        min_distance: parameters.min_distance,
        max_distance: parameters.max_distance,
        use_2d: parameters.use_2d,
        use_chirality: parameters.include_chirality,
        count_simulation: parameters.count_simulation,
        count_bounds: parameters.count_bounds.clone(),
        num_bits_per_feature: parameters.num_bits_per_feature,
        ..Default::default()
    }
}

fn arguments(expected: &GoldenOutput, arguments: &GoldenArguments) -> FingerprintFuncArguments {
    let mut additional_output = expected.additional_output.as_ref().map(|_| {
        let mut output = AdditionalOutput::new();
        output.allocate_atom_counts();
        output.allocate_atom_to_bits();
        output.allocate_bit_info_map();
        output.allocate_atoms_per_bit();
        output
    });
    FingerprintFuncArguments {
        from_atoms: arguments.from_atoms.clone(),
        ignore_atoms: arguments.ignore_atoms.clone(),
        custom_atom_invariants: arguments.custom_atom_invariants.clone(),
        additional_output: additional_output.take(),
        ..Default::default()
    }
}

fn assert_additional_output(
    context: &str,
    expected: Option<&GoldenAdditionalOutput>,
    actual: Option<&AdditionalOutput>,
) {
    match (expected, actual) {
        (None, None) => {}
        (Some(expected), Some(actual)) => {
            assert_eq!(
                actual.atom_counts.as_ref(),
                Some(&expected.atom_counts),
                "{context}"
            );
            assert_eq!(
                actual.atom_to_bits.as_ref(),
                Some(&expected.atom_to_bits),
                "{context}"
            );
            assert_eq!(
                actual.bit_info_map.as_ref(),
                Some(&expected.bit_info_map),
                "{context}"
            );
            assert_eq!(
                actual.atoms_per_bit.as_ref(),
                Some(&expected.atoms_per_bit),
                "{context}"
            );
        }
        _ => panic!("{context}: AdditionalOutput presence mismatch"),
    }
}

fn assert_branch(
    row: usize,
    smiles: &str,
    branch_name: &str,
    branch: &GoldenBranch,
    molecule: &Molecule,
) {
    let context = format!("row {row} ({smiles}) branch {branch_name}");
    assert_eq!(
        branch.parameters.name, branch_name,
        "{context}: option identity"
    );
    let generator = AtomPairFingerprintGenerator::new(&params(&branch.parameters))
        .unwrap_or_else(|error| panic!("{context}: generator creation failed: {error}"));

    let mut args = arguments(&branch.sparse_count, &branch.resolved_arguments);
    let actual = generator.sparse_count_fingerprint(molecule, &mut args);
    if branch.sparse_count.ok {
        let actual = actual.unwrap_or_else(|error| panic!("{context} sparse count: {error}"));
        assert_eq!(
            actual.size(),
            branch.sparse_count.length.unwrap(),
            "{context}"
        );
        assert_eq!(
            actual.nonzero_elements(),
            branch.sparse_count.nonzero_elements.as_ref().unwrap(),
            "{context}: sparse count elements"
        );
        assert_additional_output(
            &format!("{context}: sparse count provenance"),
            branch.sparse_count.additional_output.as_ref(),
            args.additional_output.as_ref(),
        );
    } else {
        assert!(
            actual.is_err(),
            "{context}: sparse count unexpectedly succeeded; RDKit error: {:?}",
            branch.sparse_count.error
        );
    }

    let mut args = arguments(&branch.count, &branch.resolved_arguments);
    let actual = generator.count_fingerprint(molecule, &mut args);
    if branch.count.ok {
        let actual = actual.unwrap_or_else(|error| panic!("{context} count: {error}"));
        assert_eq!(actual.size(), branch.count.length.unwrap(), "{context}");
        assert_eq!(
            actual.nonzero_elements(),
            branch.count.nonzero_elements.as_ref().unwrap(),
            "{context}: count elements"
        );
        assert_additional_output(
            &format!("{context}: count provenance"),
            branch.count.additional_output.as_ref(),
            args.additional_output.as_ref(),
        );
    } else {
        assert!(
            actual.is_err(),
            "{context}: count unexpectedly succeeded; RDKit error: {:?}",
            branch.count.error
        );
    }

    let mut args = arguments(&branch.sparse_bit, &branch.resolved_arguments);
    let actual = generator.sparse_bit_fingerprint(molecule, &mut args);
    if branch.sparse_bit.ok {
        let actual = actual.unwrap_or_else(|error| panic!("{context} sparse bit: {error}"));
        assert_eq!(
            actual.size(),
            branch.sparse_bit.length.unwrap(),
            "{context}"
        );
        let actual_bits: Vec<u64> = actual.on_bits().iter().copied().collect();
        assert_eq!(
            actual_bits,
            branch.sparse_bit.on_bits.as_ref().unwrap().as_slice(),
            "{context}: sparse bits"
        );
        assert_additional_output(
            &format!("{context}: sparse bit provenance"),
            branch.sparse_bit.additional_output.as_ref(),
            args.additional_output.as_ref(),
        );
    } else {
        assert!(
            actual.is_err(),
            "{context}: sparse bit unexpectedly succeeded; RDKit error: {:?}",
            branch.sparse_bit.error
        );
    }

    let mut args = arguments(&branch.explicit_bit, &branch.resolved_arguments);
    let actual = generator.fingerprint(molecule, &mut args);
    if branch.explicit_bit.ok {
        let actual = actual.unwrap_or_else(|error| panic!("{context} explicit bit: {error}"));
        assert_eq!(
            actual.n_bits() as u64,
            branch.explicit_bit.length.unwrap(),
            "{context}"
        );
        let actual_bits: Vec<u64> = actual.on_bits().iter().map(|&bit| bit as u64).collect();
        assert_eq!(
            actual_bits,
            branch.explicit_bit.on_bits.as_ref().unwrap().as_slice(),
            "{context}: explicit bits"
        );
        assert_additional_output(
            &format!("{context}: explicit provenance"),
            branch.explicit_bit.additional_output.as_ref(),
            args.additional_output.as_ref(),
        );
    } else {
        assert!(
            actual.is_err(),
            "{context}: explicit bit unexpectedly succeeded; RDKit error: {:?}",
            branch.explicit_bit.error
        );
    }

    let actual_ok =
        branch.sparse_count.ok && branch.count.ok && branch.sparse_bit.ok && branch.explicit_bit.ok;
    assert_eq!(branch.ok, actual_ok, "{context}: aggregate status");
    assert_eq!(
        branch.error.is_none(),
        branch.ok,
        "{context}: aggregate error"
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
                    assert!(record.error.is_some());
                    return;
                }
                assert!(record.error.is_none());
                let molecule = Molecule::from_smiles(expected_smiles).unwrap_or_else(|error| {
                    panic!("{profile} row {row} ({expected_smiles}) parse failed: {error}")
                });
                for (branch_name, branch) in &record.branches {
                    assert_branch(row, expected_smiles, branch_name, branch, &molecule);
                }
            }))
            .err()
            .map(|payload| {
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
    assert!(
        failures.is_empty(),
        "{profile} AtomPair parity failures: {failures:?}"
    );
}

#[test]
fn focused_and_small() {
    assert_profile(
        "atom_pair_focused",
        &parity_data::repo_root()
            .join("testdata/fingerprint/fixtures/rdkit/atom_pair_fingerprint_focused.smi"),
    );
    assert_profile(
        "smiles_small",
        &parity_data::repo_root().join("testdata/smiles/corpus/smiles_small.smi"),
    );
}

#[test]
fn smiles_5000() {
    assert_profile(
        "smiles_5000",
        &parity_data::repo_root().join("testdata/smiles/corpus/smiles_5000.smi"),
    );
}
