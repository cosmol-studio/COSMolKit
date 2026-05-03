use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use cosmolkit_core::{BatchErrorMode, BondOrder, Molecule, MoleculeBatch, PreparedDrawMolecule};
use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct PreparedDrawAtomRecord {
    idx: usize,
    atomic_num: u8,
    x: f64,
    y: f64,
}

#[derive(Debug, Deserialize)]
struct PreparedDrawBondRecord {
    idx: usize,
    begin: usize,
    end: usize,
    bond_type: String,
    is_aromatic: bool,
    dir: String,
}

#[derive(Debug, Deserialize)]
struct PreparedDrawRecord {
    smiles: String,
    rdkit_ok: bool,
    atoms: Option<Vec<PreparedDrawAtomRecord>>,
    bonds: Option<Vec<PreparedDrawBondRecord>>,
    error: Option<String>,
}

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn load_golden() -> Vec<PreparedDrawRecord> {
    let path = repo_root().join("tests/golden/prepared_draw_molecule.jsonl");
    let file = File::open(&path).unwrap_or_else(|err| {
        panic!(
            "failed to open {}; regenerate with tests/scripts/gen_rdkit_prepared_draw_golden.py: {err}",
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

fn bond_order_name(order: BondOrder) -> &'static str {
    match order {
        BondOrder::Single => "SINGLE",
        BondOrder::Double => "DOUBLE",
        BondOrder::Triple => "TRIPLE",
        BondOrder::Quadruple => "QUADRUPLE",
        BondOrder::Aromatic => "AROMATIC",
        BondOrder::Dative => "DATIVE",
        BondOrder::Null => "UNSPECIFIED",
    }
}

#[test]
fn prepared_draw_golden_has_one_record_per_smiles() {
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
        "prepared draw golden row count must match tests/smiles.smi"
    );
}

#[test]
fn prepared_draw_molecule_matches_rdkit_golden() {
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

        let mol = Molecule::from_smiles(&record.smiles).unwrap_or_else(|err| {
            panic!(
                "cosmolkit failed to parse row {} ({}): {err}",
                row_idx + 1,
                record.smiles
            )
        });
        let actual = mol.prepare_for_drawing_parity().unwrap_or_else(|err| {
            panic!(
                "cosmolkit failed to prepare row {} ({}): {err}",
                row_idx + 1,
                record.smiles
            )
        });
        let expected_atoms = record.atoms.as_ref().expect("rdkit ok row has atoms");
        let expected_bonds = record.bonds.as_ref().expect("rdkit ok row has bonds");

        assert_eq!(
            actual.atoms.len(),
            expected_atoms.len(),
            "row {} ({}) atom count mismatch",
            row_idx + 1,
            record.smiles
        );
        assert_eq!(
            actual.bonds.len(),
            expected_bonds.len(),
            "row {} ({}) bond count mismatch",
            row_idx + 1,
            record.smiles
        );

        for (atom_idx, (actual_atom, expected_atom)) in
            actual.atoms.iter().zip(expected_atoms).enumerate()
        {
            assert_eq!(
                actual_atom.index,
                expected_atom.idx,
                "row {} atom {atom_idx} index",
                row_idx + 1
            );
            assert_eq!(
                actual_atom.atomic_num,
                expected_atom.atomic_num,
                "row {} atom {atom_idx} atomic number",
                row_idx + 1
            );
            assert!(
                (actual_atom.x - expected_atom.x).abs() <= 1e-8,
                "row {} ({}) atom {atom_idx} x mismatch: expected {}, got {}",
                row_idx + 1,
                record.smiles,
                expected_atom.x,
                actual_atom.x
            );
            assert!(
                (actual_atom.y - expected_atom.y).abs() <= 1e-8,
                "row {} ({}) atom {atom_idx} y mismatch: expected {}, got {}",
                row_idx + 1,
                record.smiles,
                expected_atom.y,
                actual_atom.y
            );
        }

        for (bond_idx, (actual_bond, expected_bond)) in
            actual.bonds.iter().zip(expected_bonds).enumerate()
        {
            assert_eq!(
                actual_bond.index,
                expected_bond.idx,
                "row {} bond {bond_idx} index",
                row_idx + 1
            );
            assert_eq!(
                actual_bond.begin_atom,
                expected_bond.begin,
                "row {} bond {bond_idx} begin",
                row_idx + 1
            );
            assert_eq!(
                actual_bond.end_atom,
                expected_bond.end,
                "row {} bond {bond_idx} end",
                row_idx + 1
            );
            assert_eq!(
                bond_order_name(actual_bond.bond_type),
                expected_bond.bond_type,
                "row {} ({}) bond {bond_idx} type",
                row_idx + 1,
                record.smiles
            );
            assert_eq!(
                actual_bond.is_aromatic,
                expected_bond.is_aromatic,
                "row {} ({}) bond {bond_idx} aromatic flag",
                row_idx + 1,
                record.smiles
            );
            assert_eq!(
                actual_bond.rdkit_direction_name,
                expected_bond.dir,
                "row {} ({}) bond {bond_idx} direction",
                row_idx + 1,
                record.smiles
            );
        }
    }

    let smiles = records
        .iter()
        .map(|record| record.smiles.clone())
        .collect::<Vec<_>>();
    let batch = MoleculeBatch::from_smiles_list(&smiles, BatchErrorMode::Keep)
        .expect("batch SMILES parse should not raise in keep mode");
    let prepared = batch
        .prepare_for_drawing_parity_list()
        .expect("batch prepared drawing should succeed");
    for (row_idx, (record, actual)) in records.iter().zip(prepared).enumerate() {
        if !record.rdkit_ok {
            continue;
        }
        let actual =
            actual.unwrap_or_else(|| panic!("batch prepared draw missing at row {}", row_idx + 1));
        assert_prepared_draw_matches(row_idx, record, &actual);
    }
}

fn assert_prepared_draw_matches(
    row_idx: usize,
    record: &PreparedDrawRecord,
    actual: &PreparedDrawMolecule,
) {
    let expected_atoms = record.atoms.as_ref().expect("rdkit ok row has atoms");
    let expected_bonds = record.bonds.as_ref().expect("rdkit ok row has bonds");

    assert_eq!(
        actual.atoms.len(),
        expected_atoms.len(),
        "row {} ({}) atom count mismatch",
        row_idx + 1,
        record.smiles
    );
    assert_eq!(
        actual.bonds.len(),
        expected_bonds.len(),
        "row {} ({}) bond count mismatch",
        row_idx + 1,
        record.smiles
    );

    for (atom_idx, (actual_atom, expected_atom)) in
        actual.atoms.iter().zip(expected_atoms).enumerate()
    {
        assert_eq!(
            actual_atom.index,
            expected_atom.idx,
            "row {} atom {atom_idx} index",
            row_idx + 1
        );
        assert_eq!(
            actual_atom.atomic_num,
            expected_atom.atomic_num,
            "row {} atom {atom_idx} atomic number",
            row_idx + 1
        );
        assert!(
            (actual_atom.x - expected_atom.x).abs() <= 1e-8,
            "row {} ({}) atom {atom_idx} x mismatch: expected {}, got {}",
            row_idx + 1,
            record.smiles,
            expected_atom.x,
            actual_atom.x
        );
        assert!(
            (actual_atom.y - expected_atom.y).abs() <= 1e-8,
            "row {} ({}) atom {atom_idx} y mismatch: expected {}, got {}",
            row_idx + 1,
            record.smiles,
            expected_atom.y,
            actual_atom.y
        );
    }

    for (bond_idx, (actual_bond, expected_bond)) in
        actual.bonds.iter().zip(expected_bonds).enumerate()
    {
        assert_eq!(
            actual_bond.index,
            expected_bond.idx,
            "row {} bond {bond_idx} index",
            row_idx + 1
        );
        assert_eq!(
            actual_bond.begin_atom,
            expected_bond.begin,
            "row {} bond {bond_idx} begin",
            row_idx + 1
        );
        assert_eq!(
            actual_bond.end_atom,
            expected_bond.end,
            "row {} bond {bond_idx} end",
            row_idx + 1
        );
        assert_eq!(
            bond_order_name(actual_bond.bond_type),
            expected_bond.bond_type,
            "row {} ({}) bond {bond_idx} type",
            row_idx + 1,
            record.smiles
        );
        assert_eq!(
            actual_bond.is_aromatic,
            expected_bond.is_aromatic,
            "row {} ({}) bond {bond_idx} aromatic flag",
            row_idx + 1,
            record.smiles
        );
        assert_eq!(
            actual_bond.rdkit_direction_name,
            expected_bond.dir,
            "row {} ({}) bond {bond_idx} direction",
            row_idx + 1,
            record.smiles
        );
    }
}
