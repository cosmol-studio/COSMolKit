use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use cosmolkit_core::io::{
    molblock::{
        MolBlockWriteParams, SdfFormat, mol_to_2d_sdf_record_with_params,
        mol_to_3d_sdf_record_with_params,
    },
    molfile::read_mol_record_from_str,
};
use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct SdfWriteRecord {
    smiles: String,
    rdkit_ok: bool,
    source_2d_molblock: Option<String>,
    source_3d_molblock: Option<String>,
    branches: BTreeMap<String, SdfWriteBranch>,
    error: Option<String>,
}

#[derive(Debug, Deserialize)]
struct SdfWriteBranch {
    params: SdfWriteBranchParams,
    ok: bool,
    body: Option<String>,
    error: Option<String>,
}

#[derive(Debug, Deserialize)]
struct SdfWriteBranchParams {
    dimension: String,
    include_stereo: bool,
    kekulize: bool,
    force_v3000: bool,
}

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn load_golden() -> Vec<SdfWriteRecord> {
    let path = repo_root().join("tests/golden/sdf_write.jsonl");
    let file = File::open(&path).unwrap_or_else(|err| {
        panic!(
            "failed to open {}; regenerate with tests/scripts/gen_rdkit_sdf_write_golden.py: {err}",
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

fn molblock_body(record: &str) -> String {
    record
        .lines()
        .skip(3)
        .filter(|line| *line != "$$$$")
        .collect::<Vec<_>>()
        .join("\n")
}

fn branch_params(branch: &SdfWriteBranch) -> MolBlockWriteParams {
    MolBlockWriteParams {
        format: if branch.params.force_v3000 {
            SdfFormat::V3000
        } else {
            SdfFormat::Auto
        },
        include_stereo: branch.params.include_stereo,
        kekulize: branch.params.kekulize,
    }
}

#[test]
fn sdf_write_golden_has_one_record_per_smiles_and_expected_branches() {
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
        "SDF write golden row count must match tests/smiles.smi"
    );
    let first_ok = records
        .iter()
        .find(|record| record.rdkit_ok)
        .expect("SDF write golden has no RDKit-ok records");
    assert_eq!(
        first_ok.branches.len(),
        16,
        "SDF write golden should cover 2 dimensions x includeStereo x kekulize x V2000/V3000"
    );
}

#[test]
fn sdf_write_matches_rdkit_for_supported_branches_and_reports_unimplemented_ones() {
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

        let mol_2d = record
            .source_2d_molblock
            .as_ref()
            .and_then(|block| read_mol_record_from_str(block).ok())
            .map(|record| record.molecule);
        let mol_3d = record
            .source_3d_molblock
            .as_ref()
            .and_then(|block| read_mol_record_from_str(block).ok())
            .map(|record| record.molecule);

        for (branch_name, branch) in &record.branches {
            let params = branch_params(branch);
            let actual = match branch.params.dimension.as_str() {
                "2d" => mol_2d
                    .as_ref()
                    .map(|mol| mol_to_2d_sdf_record_with_params(mol, &params)),
                "3d" => mol_3d
                    .as_ref()
                    .map(|mol| mol_to_3d_sdf_record_with_params(mol, &params)),
                other => panic!("unknown SDF write dimension branch '{other}'"),
            };
            let Some(actual) = actual else {
                assert!(
                    !branch.ok,
                    "row {} ({}) branch {} has RDKit output but no COSMolKit source molecule",
                    row_idx + 1,
                    record.smiles,
                    branch_name
                );
                continue;
            };

            match actual {
                Ok(actual) => {
                    assert!(
                        branch.ok,
                        "row {} ({}) branch {}: COSMolKit succeeded but RDKit recorded {:?}",
                        row_idx + 1,
                        record.smiles,
                        branch_name,
                        branch.error
                    );
                    let expected = branch.body.as_ref().unwrap_or_else(|| {
                        panic!(
                            "row {} ({}) branch {}: RDKit ok but missing body",
                            row_idx + 1,
                            record.smiles,
                            branch_name
                        )
                    });
                    assert_eq!(
                        molblock_body(&actual),
                        *expected,
                        "SDF write molblock body mismatch at row {} ({}) branch {}",
                        row_idx + 1,
                        record.smiles,
                        branch_name
                    );
                }
                Err(err) => {
                    assert!(
                        err.to_string().contains("not ported yet"),
                        "row {} ({}) branch {} should expose explicit unimplemented error, got: {err}",
                        row_idx + 1,
                        record.smiles,
                        branch_name
                    );
                }
            }
        }
    }
}
