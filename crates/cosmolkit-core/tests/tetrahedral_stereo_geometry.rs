use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use cosmolkit_core::{BatchRecord, LigandRef, Molecule, MoleculeBatch};
use serde::Deserialize;

mod common;
use common::parity_data;

#[derive(Debug, Deserialize)]
struct GeometryRecord {
    smiles: String,
    rdkit_ok: bool,
    centers: Vec<usize>,
    positions: Option<Vec<[f64; 3]>>,
    error: Option<String>,
}

fn golden_path() -> PathBuf {
    parity_data::golden_path("tetrahedral_stereo_geometry.jsonl")
}

fn ensure_golden_exists() {
    let path = golden_path();
    assert!(
        path.exists(),
        "missing RDKit tetrahedral stereo geometry golden: {}. Generate it before running tests:\n\
         uv sync --group dev && .venv/bin/python tools/testdata/rdkit/generate_all.py --python .venv/bin/python --profile {} --suite all --clean --jobs 4",
        path.display(),
        parity_data::profile_name()
    );
}

fn load_golden() -> Vec<GeometryRecord> {
    ensure_golden_exists();
    let file = File::open(golden_path()).expect("should read tetrahedral stereo geometry golden");
    BufReader::new(file)
        .lines()
        .map(|line| line.expect("golden line should be readable"))
        .filter(|line| !line.trim().is_empty())
        .map(|line| serde_json::from_str(&line).expect("golden line should deserialize"))
        .collect()
}

fn sub(lhs: [f64; 3], rhs: [f64; 3]) -> [f64; 3] {
    [lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2]]
}

fn cross(lhs: [f64; 3], rhs: [f64; 3]) -> [f64; 3] {
    [
        lhs[1] * rhs[2] - lhs[2] * rhs[1],
        lhs[2] * rhs[0] - lhs[0] * rhs[2],
        lhs[0] * rhs[1] - lhs[1] * rhs[0],
    ]
}

fn dot(lhs: [f64; 3], rhs: [f64; 3]) -> f64 {
    lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2]
}

fn ligand_position(ligand: LigandRef, positions: &[[f64; 3]]) -> Option<[f64; 3]> {
    match ligand {
        LigandRef::Atom(index) => positions.get(index.index()).copied(),
        LigandRef::ImplicitHydrogen => None,
    }
}

fn oriented_volume(center: [f64; 3], ligands: [[f64; 3]; 3]) -> f64 {
    let v0 = sub(ligands[0], center);
    let v1 = sub(ligands[1], center);
    let v2 = sub(ligands[2], center);
    dot(v0, cross(v1, v2))
}

#[test]
fn tetrahedral_stereo_ordered_ligands_match_rdkit_etkdg_positive_volume() {
    // The ligand-order contract is defined in dev/tetrahedral_stereo.md.
    let golden = load_golden();
    assert!(
        golden.len() >= 7,
        "golden should cover tetrahedral stereocenters in testdata/smiles/corpus/smiles_small.smi"
    );

    for (row_idx, record) in golden.iter().enumerate() {
        if !record.rdkit_ok {
            let _ = &record.error;
            continue;
        }
        let positions = record.positions.as_ref().expect("positions missing");
        let mol = Molecule::from_smiles(&record.smiles).expect("COSMolKit should parse golden SMILES");
        let stereos = mol
            .tetrahedral_stereo()
            .expect("COSMolKit should report tetrahedral stereo");
        let centers: Vec<usize> = stereos.iter().map(|stereo| stereo.center.index()).collect();
        assert_eq!(
            centers,
            record.centers,
            "tetrahedral center mismatch at row {} ({})",
            row_idx + 1,
            record.smiles
        );

        for stereo in stereos {
            let center = positions
                .get(stereo.center.index())
                .copied()
                .expect("center coordinate missing");
            let first_three = [
                ligand_position(stereo.ligands[0], positions).expect("ligand 0 coordinate missing"),
                ligand_position(stereo.ligands[1], positions).expect("ligand 1 coordinate missing"),
                ligand_position(stereo.ligands[2], positions).expect("ligand 2 coordinate missing"),
            ];
            let volume = oriented_volume(center, first_three);
            assert!(
                volume > 1.0e-8,
                "ordered tetrahedral ligands should produce positive RDKit ETKDG volume at row {} center {} ({}), got {} from {:?}",
                row_idx + 1,
                stereo.center.index(),
                record.smiles,
                volume,
                stereo.ligands
            );
        }

        let batch_smiles = vec![record.smiles.clone(), record.smiles.clone()];
        let batch = MoleculeBatch::from_smiles_list(&batch_smiles);
        for (batch_idx, batch_record) in batch.iter().enumerate() {
            let BatchRecord::Molecule(batch_mol) = batch_record else {
                panic!(
                    "batch tetrahedral stereo molecule missing at row {} ({}) duplicate {}",
                    row_idx + 1,
                    record.smiles,
                    batch_idx
                );
            };
            let batch_stereos = batch_mol
                .tetrahedral_stereo()
                .expect("batch molecule should report tetrahedral stereo");
            let batch_centers: Vec<usize> = batch_stereos.iter().map(|stereo| stereo.center.index()).collect();
            assert_eq!(
                batch_centers,
                record.centers,
                "batch tetrahedral center mismatch at row {} ({}) duplicate {}",
                row_idx + 1,
                record.smiles,
                batch_idx
            );
            for stereo in batch_stereos {
                let center = positions
                    .get(stereo.center.index())
                    .copied()
                    .expect("center coordinate missing");
                let first_three = [
                    ligand_position(stereo.ligands[0], positions).expect("ligand 0 coordinate missing"),
                    ligand_position(stereo.ligands[1], positions).expect("ligand 1 coordinate missing"),
                    ligand_position(stereo.ligands[2], positions).expect("ligand 2 coordinate missing"),
                ];
                let volume = oriented_volume(center, first_three);
                assert!(
                    volume > 1.0e-8,
                    "batch ordered tetrahedral ligands should produce positive RDKit ETKDG volume at row {} center {} ({}) duplicate {}, got {} from {:?}",
                    row_idx + 1,
                    stereo.center.index(),
                    record.smiles,
                    batch_idx,
                    volume,
                    stereo.ligands
                );
            }
        }
    }
}
