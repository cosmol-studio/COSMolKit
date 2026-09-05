use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use cosmolkit_core::fingerprint::topological_fingerprint;
use cosmolkit_core::{Molecule, TopologicalFingerprintParams};
use rayon::prelude::*;
use serde::Deserialize;

mod common;
use common::parity_data;

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
    ok: bool,
    on_bits: Option<Vec<usize>>,
    num_bits: Option<usize>,
    num_on_bits: Option<usize>,
    error: Option<String>,
}

#[derive(Debug, Deserialize)]
struct GoldenParameters {
    name: String,
    #[serde(rename = "minPath")]
    min_path: u32,
    #[serde(rename = "maxPath")]
    max_path: u32,
    #[serde(rename = "fpSize")]
    fp_size: u32,
    #[serde(rename = "nBitsPerHash")]
    num_bits_per_feature: u32,
    #[serde(rename = "useHs")]
    use_hs: bool,
    #[serde(rename = "tgtDensity")]
    target_density: f64,
    #[serde(rename = "minSize")]
    min_size: u32,
    #[serde(rename = "branchedPaths")]
    branched_paths: bool,
    #[serde(rename = "useBondOrder")]
    use_bond_order: bool,
    #[serde(rename = "fromAtoms")]
    from_atoms: Option<String>,
    #[serde(rename = "atomInvariants")]
    atom_invariants: Option<String>,
}

fn load_golden() -> Vec<GoldenRecord> {
    let path = parity_data::golden_path("rdkit_topological_fingerprint.jsonl");
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

fn corpus_smiles() -> Vec<String> {
    let path = parity_data::smiles_path();
    std::fs::read_to_string(&path)
        .unwrap_or_else(|error| panic!("failed to read {}: {error}", path.display()))
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(str::to_owned)
        .collect()
}

fn params_from_golden(
    parameters: &GoldenParameters,
    atom_count: usize,
) -> TopologicalFingerprintParams {
    let mut params = TopologicalFingerprintParams {
        min_path: parameters.min_path,
        max_path: parameters.max_path,
        fp_size: parameters.fp_size,
        num_bits_per_feature: parameters.num_bits_per_feature,
        use_hs: parameters.use_hs,
        target_density: parameters.target_density,
        min_size: parameters.min_size,
        branched_paths: parameters.branched_paths,
        use_bond_order: parameters.use_bond_order,
        ..Default::default()
    };
    if parameters.from_atoms.as_deref() == Some("first") && atom_count > 0 {
        params.from_atoms = Some(vec![0]);
    }
    if parameters.atom_invariants.as_deref() == Some("index_plus_one") {
        params.atom_invariants = Some((1..=atom_count as u32).collect());
    }
    params
}

#[test]
fn rdkit_topological_fingerprint_golden_has_one_record_per_active_corpus_row() {
    let corpus = corpus_smiles();
    let golden = load_golden();
    assert_eq!(
        golden.len(),
        corpus.len(),
        "RDKFingerprint golden must have one row per active corpus input"
    );
    for (row, (record, smiles)) in golden.iter().zip(&corpus).enumerate() {
        assert_eq!(record.row, row, "golden row index changed at row {row}");
        assert_eq!(record.smiles, *smiles, "golden SMILES changed at row {row}");
        if record.rdkit_ok {
            assert!(
                record.error.is_none(),
                "RDKit-success row {row} has an error"
            );
            assert!(
                !record.branches.is_empty(),
                "RDKit-success row {row} has no branches"
            );
        } else {
            assert!(
                record.branches.is_empty(),
                "RDKit-failed row {row} has branches"
            );
            assert!(
                record.error.is_some(),
                "RDKit-failed row {row} has no error"
            );
        }
    }
}

#[test]
fn rdkit_topological_fingerprint_matches_every_active_golden_profile_exactly() {
    let corpus = corpus_smiles();
    let golden = load_golden();
    assert_eq!(golden.len(), corpus.len());

    let expected_branch_names = golden
        .first()
        .map(|record| record.branches.keys().cloned().collect::<Vec<_>>())
        .unwrap_or_default();
    assert!(
        !expected_branch_names.is_empty(),
        "RDKFingerprint profile is empty"
    );

    golden
        .par_iter()
        .zip(corpus.par_iter())
        .enumerate()
        .for_each(|(row, (record, smiles))| {
            if !record.rdkit_ok {
                let parse_result = Molecule::from_smiles(smiles);
                assert!(parse_result.is_err(), "row {row} RDKit rejected but COSMolKit parsed");
                return;
            }

            let molecule = Molecule::from_smiles(smiles)
                .unwrap_or_else(|error| panic!("row {row} ({smiles}) failed to parse: {error}"));
            assert_eq!(
                record.branches.keys().cloned().collect::<Vec<_>>(),
                expected_branch_names,
                "RDKFingerprint branch set changed at row {row}"
            );

            for branch_name in &expected_branch_names {
                let branch = record
                    .branches
                    .get(branch_name)
                    .unwrap_or_else(|| panic!("row {row} missing branch {branch_name}"));
                assert_eq!(
                    branch.parameters.name, *branch_name,
                    "row {row} branch metadata mismatch"
                );
                let params = params_from_golden(&branch.parameters, molecule.num_atoms());
                let actual = topological_fingerprint(&molecule, &params);

                if !branch.ok {
                    assert!(
                        actual.is_err(),
                        "row {row} ({smiles}) branch {branch_name} unexpectedly succeeded; RDKit error: {:?}",
                        branch.error
                    );
                    continue;
                }

                let actual =
                    actual.unwrap_or_else(|error| panic!("row {row} ({smiles}) branch {branch_name} failed: {error}"));
                assert_eq!(
                    actual.n_bits(),
                    branch.num_bits.unwrap_or_else(|| panic!("row {row} missing num_bits")),
                    "row {row} ({smiles}) branch {branch_name} fingerprint size mismatch"
                );
                assert_eq!(
                    actual.on_bits(),
                    branch
                        .on_bits
                        .clone()
                        .unwrap_or_else(|| panic!("row {row} missing on_bits")),
                    "row {row} ({smiles}) branch {branch_name} exact bits mismatch"
                );
                assert_eq!(
                    actual.on_bits().len(),
                    branch
                        .num_on_bits
                        .unwrap_or_else(|| panic!("row {row} missing num_on_bits")),
                    "row {row} ({smiles}) branch {branch_name} on-bit count mismatch"
                );
            }
        });
}
