use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use cosmolkit_core::{BondOrder, Molecule};
use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct KekulizeFlagRecord {
    smiles: String,
    parse_ok: bool,
    parse_error: Option<String>,
    kekulize_ok: bool,
    kekulize_error: Option<String>,
    atom_is_aromatic: Option<Vec<bool>>,
    bond_types: Option<Vec<String>>,
    bond_is_aromatic: Option<Vec<bool>>,
}

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn load_golden() -> Vec<KekulizeFlagRecord> {
    let path = repo_root().join("tests/golden/kekulize_clear_flags_false.jsonl");
    let file = File::open(&path).unwrap_or_else(|err| {
        panic!(
            "failed to open {}; regenerate with tests/scripts/gen_rdkit_kekulize_flags_golden.py: {err}",
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

fn bond_type_name(order: BondOrder) -> &'static str {
    match order {
        BondOrder::Null => "UNSPECIFIED",
        BondOrder::Single => "SINGLE",
        BondOrder::Double => "DOUBLE",
        BondOrder::Triple => "TRIPLE",
        BondOrder::Quadruple => "QUADRUPLE",
        BondOrder::Aromatic => "AROMATIC",
        BondOrder::Dative => "DATIVE",
    }
}

#[test]
fn kekulize_clear_flags_false_golden_has_one_record_per_smiles() {
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
    assert_eq!(records.len(), expected);
}

#[test]
fn kekulize_clear_flags_false_matches_rdkit_golden() {
    let records = load_golden();
    for (row_idx, record) in records.iter().enumerate() {
        let parsed = Molecule::from_smiles(&record.smiles);
        if record.parse_ok {
            assert!(
                parsed.is_ok(),
                "parse should succeed at row {} ({})",
                row_idx + 1,
                record.smiles
            );
        } else {
            assert!(
                parsed.is_err(),
                "parse should fail at row {} ({}) with RDKit parse_error {:?}",
                row_idx + 1,
                record.smiles,
                record.parse_error
            );
            continue;
        }
        let mut mol = parsed.expect("checked above");
        let ours = cosmolkit_core::kekulize::kekulize_in_place(&mut mol, false);
        if record.kekulize_ok {
            ours.unwrap_or_else(|err| {
                panic!(
                    "clearAromaticFlags=false kekulize failed at row {} ({}) vs RDKit success: {err}",
                    row_idx + 1,
                    record.smiles
                )
            });
        } else {
            assert!(
                ours.is_err(),
                "clearAromaticFlags=false kekulize should fail at row {} ({}) with RDKit error {:?}",
                row_idx + 1,
                record.smiles,
                record.kekulize_error
            );
            continue;
        }

        let expected_atom_is_aromatic = record.atom_is_aromatic.as_ref().unwrap_or_else(|| {
            panic!(
                "row {} ({}) missing atom_is_aromatic after successful RDKit kekulize",
                row_idx + 1,
                record.smiles
            )
        });
        let expected_bond_types = record.bond_types.as_ref().unwrap_or_else(|| {
            panic!(
                "row {} ({}) missing bond_types after successful RDKit kekulize",
                row_idx + 1,
                record.smiles
            )
        });
        let expected_bond_is_aromatic = record.bond_is_aromatic.as_ref().unwrap_or_else(|| {
            panic!(
                "row {} ({}) missing bond_is_aromatic after successful RDKit kekulize",
                row_idx + 1,
                record.smiles
            )
        });

        let actual_atom_is_aromatic: Vec<bool> =
            mol.atoms.iter().map(|atom| atom.is_aromatic).collect();
        let actual_bond_types: Vec<&'static str> = mol
            .bonds
            .iter()
            .map(|bond| bond_type_name(bond.order))
            .collect();
        let actual_bond_is_aromatic: Vec<bool> =
            mol.bonds.iter().map(|bond| bond.is_aromatic).collect();

        assert_eq!(
            actual_atom_is_aromatic,
            *expected_atom_is_aromatic,
            "atom aromatic flags mismatch at row {} ({})",
            row_idx + 1,
            record.smiles
        );
        assert_eq!(
            actual_bond_types,
            *expected_bond_types,
            "bond types mismatch at row {} ({})",
            row_idx + 1,
            record.smiles
        );
        assert_eq!(
            actual_bond_is_aromatic,
            *expected_bond_is_aromatic,
            "bond aromatic flags mismatch at row {} ({})",
            row_idx + 1,
            record.smiles
        );
    }
}
