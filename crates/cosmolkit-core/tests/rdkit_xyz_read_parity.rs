use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use cosmolkit_core::Molecule;
use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct XyzReadRecord {
    smiles: String,
    rdkit_ok: bool,
    xyz_block: Option<String>,
    atomic_numbers: Option<Vec<u8>>,
    coordinates: Option<Vec<[f64; 3]>>,
    comment: Option<String>,
    num_bonds: Option<usize>,
    error: Option<String>,
}

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn load_golden() -> Vec<XyzReadRecord> {
    let path = repo_root().join("tests/golden/xyz_read.jsonl");
    let file = File::open(&path).unwrap_or_else(|err| {
        panic!(
            "failed to open {}; regenerate all RDKit goldens with `.venv/bin/python tests/scripts/gen_all_rdkit_goldens.py --python .venv/bin/python --clean --jobs 4`: {err}",
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

#[test]
fn xyz_read_golden_has_one_record_per_smiles() {
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
        "XYZ read golden row count must match tests/smiles.smi"
    );
}

#[test]
fn xyz_reader_matches_rdkit_golden_for_smiles_generated_xyz_blocks() {
    let records = load_golden();
    let row_filter = std::env::var("COSMOLKIT_XYZ_ROW_FILTER")
        .ok()
        .and_then(|s| s.parse::<usize>().ok());

    for (row_idx, record) in records.iter().enumerate() {
        if let Some(filter) = row_filter
            && row_idx + 1 != filter
        {
            continue;
        }
        if !record.rdkit_ok {
            assert!(
                record.error.is_some(),
                "row {} ({}) is rdkit not ok but has no error",
                row_idx + 1,
                record.smiles
            );
            continue;
        }

        let xyz_block = record.xyz_block.as_ref().unwrap_or_else(|| {
            panic!(
                "row {} ({}) is rdkit_ok but has no xyz block",
                row_idx + 1,
                record.smiles
            )
        });
        let expected_atomic_numbers = record.atomic_numbers.as_ref().unwrap_or_else(|| {
            panic!(
                "row {} ({}) is rdkit_ok but has no atomic numbers",
                row_idx + 1,
                record.smiles
            )
        });
        let expected_coordinates = record.coordinates.as_ref().unwrap_or_else(|| {
            panic!(
                "row {} ({}) is rdkit_ok but has no coordinates",
                row_idx + 1,
                record.smiles
            )
        });

        let molecule = Molecule::from_xyz_block(xyz_block).unwrap_or_else(|err| {
            panic!(
                "COSMolKit failed to parse RDKit XYZ block at row {} ({}): {err}",
                row_idx + 1,
                record.smiles
            )
        });

        assert_eq!(
            molecule.num_atoms(),
            expected_atomic_numbers.len(),
            "atom count mismatch at row {} ({})",
            row_idx + 1,
            record.smiles
        );
        assert_eq!(
            molecule.num_bonds(),
            record.num_bonds.expect("RDKit ok row has num_bonds"),
            "XYZ parser must not infer bonds at row {} ({})",
            row_idx + 1,
            record.smiles
        );
        assert_eq!(
            molecule.atomic_numbers(),
            *expected_atomic_numbers,
            "atomic number mismatch at row {} ({})",
            row_idx + 1,
            record.smiles
        );
        assert_eq!(
            molecule.properties().prop("_FileComments"),
            record.comment.as_deref(),
            "_FileComments mismatch at row {} ({})",
            row_idx + 1,
            record.smiles
        );
        assert_eq!(
            molecule.conformers_3d().len(),
            usize::from(!expected_coordinates.is_empty()),
            "3D conformer count mismatch at row {} ({})",
            row_idx + 1,
            record.smiles
        );
        if !expected_coordinates.is_empty() {
            let actual_coordinates = molecule.conformers_3d()[0].coords();
            assert_eq!(
                actual_coordinates.len(),
                expected_coordinates.len(),
                "coordinate row count mismatch at row {} ({})",
                row_idx + 1,
                record.smiles
            );
            for (atom_idx, (actual, expected)) in actual_coordinates
                .iter()
                .zip(expected_coordinates.iter())
                .enumerate()
            {
                for axis in 0..3 {
                    assert!(
                        (actual[axis] - expected[axis]).abs() <= 1e-10,
                        "coordinate mismatch at row {} ({}) atom {} axis {}: ours={} expected={}",
                        row_idx + 1,
                        record.smiles,
                        atom_idx,
                        axis,
                        actual[axis],
                        expected[axis]
                    );
                }
            }
        }
    }
}
