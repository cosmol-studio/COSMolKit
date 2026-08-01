use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::Command;

use cosmolkit_core::{
    BondOrder, Molecule,
    io::{
        mol2::{Mol2ReadParams, Mol2Type, read_mol2_from_str_with_params},
        molfile::read_mol_record_from_str_with_params,
        sdf::SdfReadParams,
    },
};
use serde::Deserialize;

mod common;
use common::parity_data;

#[derive(Debug, Deserialize)]
struct BuiltinFixtureRecord {
    kind: String,
    fixture: String,
    record_index: Option<usize>,
    smiles: Option<String>,
    rdkit_ok: Option<bool>,
    atom_count: Option<usize>,
    bond_count: Option<usize>,
    atomic_numbers: Option<Vec<u8>>,
    formal_charges: Option<Vec<i8>>,
    bond_types: Option<Vec<String>>,
    error: Option<String>,
    byte_len: Option<usize>,
    nonempty: Option<bool>,
}

fn repo_root() -> PathBuf {
    parity_data::repo_root()
}

fn fixture_root() -> PathBuf {
    repo_root().join("testdata/rdkit_builtin/fixtures")
}

fn load_golden() -> Vec<BuiltinFixtureRecord> {
    let path = parity_data::golden_path("rdkit_builtin_fixture_migration.jsonl");
    let file = File::open(&path).unwrap_or_else(|err| {
        panic!(
            "failed to open {}; regenerate RDKit goldens with `{}`: {err}",
            path.display(),
            parity_data::regenerate_command()
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

fn fixture_path(record: &BuiltinFixtureRecord) -> PathBuf {
    fixture_root().join(&record.fixture)
}

fn read_text(path: &Path) -> String {
    std::fs::read_to_string(path)
        .unwrap_or_else(|err| panic!("failed to read fixture {}: {err}", path.display()))
}

fn read_gzip_text(path: &Path) -> String {
    let output = Command::new("gzip")
        .arg("-dc")
        .arg(path)
        .output()
        .unwrap_or_else(|err| panic!("failed to run gzip for {}: {err}", path.display()));
    assert!(
        output.status.success(),
        "gzip -dc failed for {} with status {:?}: {}",
        path.display(),
        output.status.code(),
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout)
        .unwrap_or_else(|err| panic!("gzip output for {} is not UTF-8: {err}", path.display()))
}

fn split_sdf_records(text: &str) -> Vec<String> {
    text.split("$$$$")
        .filter(|chunk| !chunk.trim().is_empty())
        .map(|chunk| {
            let mut record = chunk.trim_end_matches('\n').to_owned();
            record.push_str("\n$$$$\n");
            record
        })
        .collect()
}

fn sdf_record_text(path: &Path, index: usize) -> Option<String> {
    let text = if path
        .file_name()
        .and_then(|name| name.to_str())
        .is_some_and(|name| name.ends_with(".gz"))
    {
        read_gzip_text(path)
    } else {
        read_text(path)
    };
    split_sdf_records(&text).into_iter().nth(index)
}

fn bond_order_name(order: BondOrder) -> &'static str {
    order.rdkit_name()
}

fn assert_molecule_summary(record: &BuiltinFixtureRecord, molecule: &Molecule, context: &str) {
    assert_eq!(
        molecule.num_atoms(),
        record.atom_count.expect("RDKit ok record has atom_count"),
        "atom count mismatch for {context}"
    );
    assert_eq!(
        molecule.num_bonds(),
        record.bond_count.expect("RDKit ok record has bond_count"),
        "bond count mismatch for {context}"
    );
    assert_eq!(
        molecule.atomic_numbers(),
        *record
            .atomic_numbers
            .as_ref()
            .expect("RDKit ok record has atomic_numbers"),
        "atomic number mismatch for {context}"
    );
    let actual_formal_charges = molecule
        .atoms()
        .iter()
        .map(|atom| atom.formal_charge())
        .collect::<Vec<_>>();
    assert_eq!(
        actual_formal_charges,
        *record
            .formal_charges
            .as_ref()
            .expect("RDKit ok record has formal_charges"),
        "formal charge mismatch for {context}"
    );
    let actual_bond_types = molecule
        .bonds()
        .iter()
        .map(|bond| bond_order_name(bond.order()).to_owned())
        .collect::<Vec<_>>();
    assert_eq!(
        actual_bond_types,
        *record
            .bond_types
            .as_ref()
            .expect("RDKit ok record has bond_types"),
        "bond type mismatch for {context}"
    );
}

fn assert_rdkit_error_has_context(record: &BuiltinFixtureRecord, row_idx: usize) {
    assert!(
        record.error.is_some(),
        "row {} fixture {} is RDKit-failing without error detail",
        row_idx + 1,
        record.fixture
    );
}

#[test]
fn rdkit_builtin_fixture_golden_covers_migrated_corpora() {
    let records = load_golden();
    assert!(
        records.len() >= 1000,
        "RDKit built-in fixture migration golden should cover the migrated corpora"
    );
    for prefix in [
        "Code/GraphMol/FileParsers/",
        "Code/GraphMol/ForceFieldHelpers_UFF/",
        "Code/GraphMol/ForceFieldHelpers_MMFF/",
        "Code/GraphMol/DistGeomHelpers/",
        "Code/GraphMol/Depictor/",
        "Code/GraphMol/GraphMol/",
        "Data/NCI/",
        "Regress/Data/",
    ] {
        assert!(
            records
                .iter()
                .any(|record| record.fixture.starts_with(prefix)),
            "RDKit built-in fixture migration golden is missing corpus prefix {prefix}"
        );
    }
}

#[test]
fn rdkit_builtin_fixture_migration_matches_rdkit_parse_status_and_summary() {
    let records = load_golden();
    for (row_idx, record) in records.iter().enumerate() {
        let path = fixture_path(record);
        assert!(
            path.exists(),
            "row {} fixture was not migrated locally: {}",
            row_idx + 1,
            path.display()
        );
        let context = format!(
            "row {} kind {} fixture {} record {:?}",
            row_idx + 1,
            record.kind,
            record.fixture,
            record.record_index
        );
        match record.kind.as_str() {
            "mol" => {
                let molblock = read_text(&path);
                let actual = read_mol_record_from_str_with_params(
                    &molblock,
                    SdfReadParams {
                        sanitize: true,
                        remove_hs: false,
                        process_property_lists: false,
                        ..Default::default()
                    },
                );
                match record.rdkit_ok {
                    Some(true) => assert_molecule_summary(
                        record,
                        &actual
                            .unwrap_or_else(|err| {
                                panic!("COSMolKit failed to parse {context}: {err}")
                            })
                            .molecule,
                        &context,
                    ),
                    Some(false) => {
                        assert_rdkit_error_has_context(record, row_idx);
                        assert!(
                            actual.is_err(),
                            "COSMolKit parsed {context}, but RDKit rejected it"
                        );
                    }
                    None => panic!("{context} has no rdkit_ok field"),
                }
            }
            "sdf" => {
                let record_text =
                    sdf_record_text(&path, record.record_index.expect("SDF record has index"));
                let actual = record_text.as_ref().map(|text| {
                    // The golden generator uses RDKit MolFromMolBlock() on each
                    // split SDF record, not SDMolSupplier. Match that API
                    // boundary here so SD data-field parsing is not part of
                    // this migrated fixture summary.
                    read_mol_record_from_str_with_params(
                        text,
                        SdfReadParams {
                            sanitize: true,
                            remove_hs: false,
                            process_property_lists: false,
                            ..Default::default()
                        },
                    )
                });
                match record.rdkit_ok {
                    Some(true) => assert_molecule_summary(
                        record,
                        &actual
                            .unwrap_or_else(|| panic!("missing local SDF record for {context}"))
                            .unwrap_or_else(|err| {
                                panic!("COSMolKit failed to parse {context}: {err}")
                            })
                            .molecule,
                        &context,
                    ),
                    Some(false) => {
                        assert_rdkit_error_has_context(record, row_idx);
                        assert!(
                            actual.is_none_or(|result| result.is_err()),
                            "COSMolKit parsed {context}, but RDKit rejected it"
                        );
                    }
                    None => panic!("{context} has no rdkit_ok field"),
                }
            }
            "mol2" => {
                let mol2 = read_text(&path);
                let actual = read_mol2_from_str_with_params(
                    &mol2,
                    Mol2ReadParams {
                        sanitize: true,
                        remove_hs: false,
                        variant: Mol2Type::Corina,
                        cleanup_substructures: true,
                    },
                );
                match record.rdkit_ok {
                    Some(true) => {
                        let parsed = actual.unwrap_or_else(|err| {
                            panic!("COSMolKit failed to parse {context}: {err}")
                        });
                        assert_molecule_summary(
                            record,
                            &parsed
                                .unwrap_or_else(|| {
                                    panic!("COSMolKit returned no molecule for {context}")
                                })
                                .molecule,
                            &context,
                        );
                    }
                    Some(false) => {
                        assert_rdkit_error_has_context(record, row_idx);
                        assert!(
                            actual.is_err() || actual.as_ref().is_ok_and(Option::is_none),
                            "COSMolKit parsed {context}, but RDKit rejected it"
                        );
                    }
                    None => panic!("{context} has no rdkit_ok field"),
                }
            }
            "xyz" => {
                let xyz = read_text(&path);
                let actual = Molecule::from_xyz_block(&xyz);
                match record.rdkit_ok {
                    Some(true) => assert_molecule_summary(
                        record,
                        &actual.unwrap_or_else(|err| {
                            panic!("COSMolKit failed to parse {context}: {err}")
                        }),
                        &context,
                    ),
                    Some(false) => {
                        assert_rdkit_error_has_context(record, row_idx);
                        assert!(
                            actual.is_err(),
                            "COSMolKit parsed {context}, but RDKit rejected it"
                        );
                    }
                    None => panic!("{context} has no rdkit_ok field"),
                }
            }
            "smi" => {
                let smiles = record
                    .smiles
                    .as_deref()
                    .unwrap_or_else(|| panic!("{context} has no smiles field"));
                let actual = Molecule::from_smiles(smiles);
                match record.rdkit_ok {
                    Some(true) => assert_molecule_summary(
                        record,
                        &actual.unwrap_or_else(|err| {
                            panic!("COSMolKit failed to parse {context}: {err}")
                        }),
                        &context,
                    ),
                    Some(false) => {
                        assert_rdkit_error_has_context(record, row_idx);
                        assert!(
                            actual.is_err(),
                            "COSMolKit parsed {context}, but RDKit rejected it"
                        );
                    }
                    None => panic!("{context} has no rdkit_ok field"),
                }
            }
            "inventory" => {
                let bytes = std::fs::read(&path)
                    .unwrap_or_else(|err| panic!("failed to read {context}: {err}"));
                assert_eq!(
                    bytes.len(),
                    record.byte_len.expect("inventory record has byte_len"),
                    "byte length mismatch for {context}"
                );
                assert_eq!(
                    !bytes.is_empty(),
                    record.nonempty.expect("inventory record has nonempty"),
                    "nonempty mismatch for {context}"
                );
                assert!(
                    record.nonempty == Some(true),
                    "empty migrated fixture: {context}"
                );
            }
            other => panic!("unknown RDKit built-in fixture migration kind {other}: {context}"),
        }
    }
}
