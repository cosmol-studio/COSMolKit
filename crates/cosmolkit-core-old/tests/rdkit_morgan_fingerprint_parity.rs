use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use cosmolkit_core::{
    BatchErrorMode, Molecule, MoleculeBatch, MorganAtomInvariantsGenerator,
    MorganBondInvariantsGenerator, MorganFingerprintParams,
};
use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct MorganBranch {
    ok: bool,
    on_bits: Option<Vec<usize>>,
    num_on_bits: Option<usize>,
    tanimoto_to_previous: Option<f64>,
    additional_output: Option<MorganAdditionalOutputGolden>,
    error: Option<String>,
}

#[derive(Debug, Deserialize)]
struct MorganAdditionalOutputGolden {
    atom_counts: Vec<u32>,
    atom_to_bits: Vec<Vec<usize>>,
    bit_info_map: BTreeMap<String, Vec<[usize; 2]>>,
    atoms_per_bit: BTreeMap<String, Vec<Vec<usize>>>,
}

#[derive(Debug, Deserialize)]
struct MorganRecord {
    smiles: String,
    rdkit_ok: bool,
    branches: BTreeMap<String, MorganBranch>,
    error: Option<String>,
}

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn load_golden() -> Vec<MorganRecord> {
    let path = repo_root().join("tests/golden/morgan_fingerprint.jsonl");
    ensure_golden_exists(&path);
    let file = File::open(&path).unwrap_or_else(|err| {
        panic!(
            "failed to open {}; regenerate with tests/scripts/gen_rdkit_morgan_fingerprint_golden.py: {err}",
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

fn ensure_golden_exists(golden_path: &PathBuf) {
    assert!(
        golden_path.exists(),
        "missing RDKit Morgan fingerprint golden: {}. Generate it before running tests:\n\
         uv sync --group dev && .venv/bin/python tests/scripts/gen_rdkit_morgan_fingerprint_golden.py --input tests/smiles.smi --output tests/golden/morgan_fingerprint.jsonl",
        golden_path.display()
    );
}

fn branch_params(name: &str, mol: &Molecule) -> MorganFingerprintParams {
    match name {
        "r2_n2048_chiral_false_bonds_true" => MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            use_chirality: false,
            use_bond_types: true,
            ..Default::default()
        },
        "r2_n2048_chiral_true_bonds_true" => MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            use_chirality: true,
            use_bond_types: true,
            ..Default::default()
        },
        "r3_n2048_chiral_false_bonds_true" => MorganFingerprintParams {
            radius: 3,
            n_bits: 2048,
            use_chirality: false,
            use_bond_types: true,
            ..Default::default()
        },
        "r2_n1024_chiral_false_bonds_false" => MorganFingerprintParams {
            radius: 2,
            n_bits: 1024,
            use_chirality: false,
            use_bond_types: false,
            ..Default::default()
        },
        "r2_n2048_count_sim_bounds_1_2_4_8" => MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            count_simulation: true,
            count_bounds: vec![1, 2, 4, 8],
            ..Default::default()
        },
        "r2_n2048_redundant_true" => MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            include_redundant_environments: true,
            ..Default::default()
        },
        "r2_n2048_ring_membership_false" => MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            // RDKit Fingerprints/Wrap/MorganWrapper.cpp accepts
            // includeRingMembership but intentionally does not forward it to
            // getMorganGenerator(); parity here follows that wrapper behavior.
            include_ring_membership: true,
            ..Default::default()
        },
        "r2_n2048_from_atom_0" => MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            from_atoms: (!mol.atoms().is_empty()).then(|| vec![0]),
            ..Default::default()
        },
        "r2_n2048_ignore_atom_0" => MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            ignore_atoms: (!mol.atoms().is_empty()).then(|| vec![0]),
            ..Default::default()
        },
        "r2_n2048_custom_atom_invariants" => MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            custom_atom_invariants: Some(
                (0..mol.atoms().len()).map(|idx| idx as u32 + 1).collect(),
            ),
            ..Default::default()
        },
        "r2_n2048_custom_bond_invariants" => MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            custom_bond_invariants: Some(
                (0..mol.bonds().len()).map(|idx| idx as u32 + 7).collect(),
            ),
            ..Default::default()
        },
        "r2_n2048_only_nonzero_custom_atom" => MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            only_nonzero_invariants: true,
            custom_atom_invariants: Some(
                (0..mol.atoms().len())
                    .map(|idx| if idx % 2 == 0 { 0 } else { idx as u32 + 1 })
                    .collect(),
            ),
            ..Default::default()
        },
        "r2_n2048_atom_inv_gen_ring_false" => MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            atom_invariants_generator: MorganAtomInvariantsGenerator::Connectivity {
                include_ring_membership: false,
            },
            ..Default::default()
        },
        "r2_n2048_feature_atom_inv_gen" => MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            atom_invariants_generator: MorganAtomInvariantsGenerator::Feature,
            ..Default::default()
        },
        "r2_n2048_bond_inv_gen_no_bond_types" => MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            bond_invariants_generator: Some(MorganBondInvariantsGenerator {
                use_bond_types: false,
                use_chirality: false,
            }),
            ..Default::default()
        },
        "r2_n2048_num_bits_per_feature_2" => MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            num_bits_per_feature: 2,
            ..Default::default()
        },
        "r2_n2048_additional_output" => MorganFingerprintParams {
            radius: 2,
            n_bits: 2048,
            collect_additional_output: true,
            ..Default::default()
        },
        other => panic!("unknown Morgan fingerprint branch '{other}'"),
    }
}

#[test]
fn morgan_fingerprint_golden_has_one_record_per_smiles() {
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
    let branch_names = records
        .iter()
        .flat_map(|record| record.branches.keys().map(String::as_str))
        .collect::<std::collections::BTreeSet<_>>();
    for required in [
        "r2_n2048_count_sim_bounds_1_2_4_8",
        "r2_n2048_redundant_true",
        "r2_n2048_ring_membership_false",
        "r2_n2048_from_atom_0",
        "r2_n2048_ignore_atom_0",
        "r2_n2048_custom_atom_invariants",
        "r2_n2048_custom_bond_invariants",
        "r2_n2048_only_nonzero_custom_atom",
        "r2_n2048_atom_inv_gen_ring_false",
        "r2_n2048_feature_atom_inv_gen",
        "r2_n2048_bond_inv_gen_no_bond_types",
        "r2_n2048_num_bits_per_feature_2",
        "r2_n2048_additional_output",
    ] {
        assert!(
            branch_names.contains(required),
            "Morgan fingerprint golden is missing required branch {required}"
        );
    }
    assert_eq!(
        records.len(),
        expected,
        "Morgan fingerprint golden row count must match tests/smiles.smi"
    );
    assert!(
        records.iter().any(|record| {
            record.rdkit_ok
                && record
                    .branches
                    .values()
                    .any(|branch| branch.tanimoto_to_previous.is_some())
        }),
        "Morgan fingerprint golden should include adjacent-row Tanimoto baselines"
    );
}

#[test]
fn morgan_fingerprint_matches_rdkit_golden_across_param_branches() {
    let records = load_golden();
    let mut previous_by_branch: BTreeMap<String, cosmolkit_core::Fingerprint> = BTreeMap::new();

    for (row_idx, record) in records.iter().enumerate() {
        if !record.rdkit_ok {
            assert!(
                record.error.is_some(),
                "row {} ({}) is rdkit not ok but has no error",
                row_idx + 1,
                record.smiles
            );
            previous_by_branch.clear();
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
            let params = branch_params(branch_name, &mol);
            let actual_output = mol.morgan_fingerprint_with_output(&params).unwrap_or_else(|err| {
                panic!(
                    "Morgan fingerprint unsupported for row {} ({}) branch {} with params {:?}; RDKit golden has {:?} on bits; error = {}",
                    row_idx + 1,
                    record.smiles,
                    branch_name,
                    params,
                    expected_branch.num_on_bits,
                    err
                )
            });
            let actual = actual_output.fingerprint.clone();

            assert!(
                expected_branch.ok,
                "row {} ({}) branch {}: cosmolkit succeeded but RDKit golden recorded error {:?}",
                row_idx + 1,
                record.smiles,
                branch_name,
                expected_branch.error
            );
            let expected_on_bits = expected_branch.on_bits.as_ref().unwrap_or_else(|| {
                panic!(
                    "row {} ({}) branch {}: RDKit branch ok but missing on_bits",
                    row_idx + 1,
                    record.smiles,
                    branch_name
                )
            });
            assert_eq!(
                actual.on_bits(),
                *expected_on_bits,
                "Morgan fingerprint bit mismatch at row {} ({}) branch {}",
                row_idx + 1,
                record.smiles,
                branch_name
            );
            assert_eq!(
                actual.on_bits().len(),
                expected_branch
                    .num_on_bits
                    .expect("RDKit ok branch should include num_on_bits"),
                "Morgan fingerprint on-bit count mismatch at row {} ({}) branch {}",
                row_idx + 1,
                record.smiles,
                branch_name
            );

            if let Some(expected_tanimoto) = expected_branch.tanimoto_to_previous {
                let previous = previous_by_branch.get(branch_name).unwrap_or_else(|| {
                    panic!(
                        "row {} ({}) branch {} has Tanimoto golden but no previous fingerprint",
                        row_idx + 1,
                        record.smiles,
                        branch_name
                    )
                });
                let actual_tanimoto = actual
                    .tanimoto(previous)
                    .expect("same branch fingerprints should have same bit length");
                assert!(
                    (actual_tanimoto - expected_tanimoto).abs() <= 1e-12,
                    "Morgan Tanimoto mismatch at row {} ({}) branch {}: actual {}, expected {}",
                    row_idx + 1,
                    record.smiles,
                    branch_name,
                    actual_tanimoto,
                    expected_tanimoto
                );
            }

            if let Some(expected_output) = expected_branch.additional_output.as_ref() {
                let actual_additional = actual_output
                    .additional_output
                    .as_ref()
                    .expect("branch with additional_output golden should collect output");
                assert_eq!(
                    actual_additional.atom_counts,
                    expected_output.atom_counts,
                    "Morgan atomCounts mismatch at row {} ({}) branch {}",
                    row_idx + 1,
                    record.smiles,
                    branch_name
                );
                assert_eq!(
                    actual_additional.atom_to_bits,
                    expected_output.atom_to_bits,
                    "Morgan atomToBits mismatch at row {} ({}) branch {}",
                    row_idx + 1,
                    record.smiles,
                    branch_name
                );
                let expected_bit_info = expected_output
                    .bit_info_map
                    .iter()
                    .map(|(bit, pairs)| {
                        (
                            bit.parse::<usize>().expect("bit id should parse"),
                            pairs
                                .iter()
                                .map(|pair| (pair[0], pair[1] as u32))
                                .collect::<Vec<_>>(),
                        )
                    })
                    .collect::<BTreeMap<_, _>>();
                assert_eq!(
                    actual_additional.bit_info_map,
                    expected_bit_info,
                    "Morgan bitInfoMap mismatch at row {} ({}) branch {}",
                    row_idx + 1,
                    record.smiles,
                    branch_name
                );
                let expected_atoms_per_bit = expected_output
                    .atoms_per_bit
                    .iter()
                    .map(|(bit, atoms_per_bit)| {
                        (
                            bit.parse::<usize>().expect("bit id should parse"),
                            atoms_per_bit.clone(),
                        )
                    })
                    .collect::<BTreeMap<_, _>>();
                assert_eq!(
                    actual_additional.atoms_per_bit,
                    expected_atoms_per_bit,
                    "Morgan atomsPerBit mismatch at row {} ({}) branch {}",
                    row_idx + 1,
                    record.smiles,
                    branch_name
                );
            }

            let batch_smiles = vec![record.smiles.clone(), record.smiles.clone()];
            let batch = MoleculeBatch::from_smiles_list(&batch_smiles, BatchErrorMode::Keep)
                .expect("batch Morgan SMILES parse should not raise in keep mode");
            let batch_outputs = batch
                .morgan_fingerprint_with_output_list(&params)
                .expect("batch Morgan fingerprint should succeed after scalar branch passed");
            for (batch_idx, batch_output) in batch_outputs.into_iter().enumerate() {
                let batch_output = batch_output.unwrap_or_else(|| {
                    panic!(
                        "batch Morgan fingerprint missing at row {} ({}) branch {} duplicate {}",
                        row_idx + 1,
                        record.smiles,
                        branch_name,
                        batch_idx
                    )
                });
                assert_eq!(
                    batch_output.fingerprint.on_bits(),
                    *expected_on_bits,
                    "batch Morgan fingerprint bit mismatch at row {} ({}) branch {} duplicate {}",
                    row_idx + 1,
                    record.smiles,
                    branch_name,
                    batch_idx
                );
                if let Some(expected_output) = expected_branch.additional_output.as_ref() {
                    let batch_additional = batch_output
                        .additional_output
                        .as_ref()
                        .expect("batch branch with additional_output golden should collect output");
                    assert_eq!(
                        batch_additional.atom_counts,
                        expected_output.atom_counts,
                        "batch Morgan atomCounts mismatch at row {} ({}) branch {} duplicate {}",
                        row_idx + 1,
                        record.smiles,
                        branch_name,
                        batch_idx
                    );
                    assert_eq!(
                        batch_additional.atom_to_bits,
                        expected_output.atom_to_bits,
                        "batch Morgan atomToBits mismatch at row {} ({}) branch {} duplicate {}",
                        row_idx + 1,
                        record.smiles,
                        branch_name,
                        batch_idx
                    );
                    let expected_bit_info = expected_output
                        .bit_info_map
                        .iter()
                        .map(|(bit, pairs)| {
                            (
                                bit.parse::<usize>().expect("bit id should parse"),
                                pairs
                                    .iter()
                                    .map(|pair| (pair[0], pair[1] as u32))
                                    .collect::<Vec<_>>(),
                            )
                        })
                        .collect::<BTreeMap<_, _>>();
                    assert_eq!(
                        batch_additional.bit_info_map,
                        expected_bit_info,
                        "batch Morgan bitInfoMap mismatch at row {} ({}) branch {} duplicate {}",
                        row_idx + 1,
                        record.smiles,
                        branch_name,
                        batch_idx
                    );
                    let expected_atoms_per_bit = expected_output
                        .atoms_per_bit
                        .iter()
                        .map(|(bit, atoms_per_bit)| {
                            (
                                bit.parse::<usize>().expect("bit id should parse"),
                                atoms_per_bit.clone(),
                            )
                        })
                        .collect::<BTreeMap<_, _>>();
                    assert_eq!(
                        batch_additional.atoms_per_bit,
                        expected_atoms_per_bit,
                        "batch Morgan atomsPerBit mismatch at row {} ({}) branch {} duplicate {}",
                        row_idx + 1,
                        record.smiles,
                        branch_name,
                        batch_idx
                    );
                }
            }

            previous_by_branch.insert(branch_name.clone(), actual);
        }
    }
}
