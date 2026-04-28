use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::process::Command;

use cosmolkit_core::{LigandRef, Molecule};
use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct GeometryRecord {
    smiles: String,
    rdkit_ok: bool,
    centers: Vec<usize>,
    positions: Option<Vec<[f64; 3]>>,
    error: Option<String>,
}

fn golden_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/golden/tetrahedral_stereo_geometry.jsonl")
}

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn ensure_golden_exists() {
    let path = golden_path();
    if path.exists() {
        return;
    }

    let root = repo_root();
    let status = Command::new(root.join(".venv/bin/python"))
        .arg("tests/scripts/gen_rdkit_tetrahedral_stereo_geometry.py")
        .current_dir(&root)
        .status()
        .expect("failed to run RDKit tetrahedral stereo geometry generator");
    assert!(
        status.success(),
        "RDKit tetrahedral stereo geometry generator failed with status {status}"
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
        LigandRef::Atom(index) => positions.get(index).copied(),
        LigandRef::ImplicitH => None,
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
    // The ligand-order contract is defined in tetrahedral_stereo_representation.md.
    let golden = load_golden();
    assert!(
        golden.len() >= 7,
        "golden should cover all tetrahedral stereocenters in tests/smiles.smi"
    );

    for (row_idx, record) in golden.iter().enumerate() {
        assert!(
            record.rdkit_ok,
            "RDKit ETKDG golden failed at row {} ({}) with error {:?}",
            row_idx + 1,
            record.smiles,
            record.error
        );
        let positions = record.positions.as_ref().expect("positions missing");
        let mol =
            Molecule::from_smiles(&record.smiles).expect("COSMolKit should parse golden SMILES");
        let stereos = mol.tetrahedral_stereo();
        let centers: Vec<usize> = stereos.iter().map(|stereo| stereo.center).collect();
        assert_eq!(
            centers,
            record.centers,
            "tetrahedral center mismatch at row {} ({})",
            row_idx + 1,
            record.smiles
        );

        for stereo in stereos {
            let center = positions
                .get(stereo.center)
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
                stereo.center,
                record.smiles,
                volume,
                stereo.ligands
            );
        }
    }
}
