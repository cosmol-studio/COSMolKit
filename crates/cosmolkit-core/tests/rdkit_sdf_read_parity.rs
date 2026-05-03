use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::process::Command;

use cosmolkit_core::{
    BatchErrorMode, BondDirection, BondOrder, BondStereo, ChiralTag, MoleculeBatch,
    SmilesWriteParams,
    io::sdf::{SdfCoordinateMode, read_sdf_record_from_str},
};
use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct SdfReadRecord {
    smiles: String,
    case_id: String,
    dimension: String,
    format: String,
    stereo_markers: String,
    rdkit_ok: bool,
    sdf: Option<String>,
    atoms: Option<Vec<AtomRecord>>,
    bonds: Option<Vec<BondRecord>>,
    positions: Option<Vec<[f64; 3]>>,
    chiral_tags: Option<Vec<String>>,
    smiles_out: Option<SmilesOut>,
    error: Option<String>,
}

#[derive(Debug, Deserialize, PartialEq, Eq)]
struct AtomRecord {
    atomic_num: u8,
    isotope: Option<u16>,
    formal_charge: i8,
    is_aromatic: bool,
    atom_map_num: Option<u32>,
}

#[derive(Debug, Deserialize, PartialEq, Eq)]
struct BondRecord {
    begin: usize,
    end: usize,
    bond_type: String,
    is_aromatic: bool,
    direction: String,
    stereo: String,
    stereo_atoms: Vec<usize>,
}

#[derive(Debug, Deserialize)]
struct SmilesOut {
    canonical: String,
    noncanonical: String,
}

fn golden_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/sdf_read.jsonl")
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
        .arg("tests/scripts/gen_rdkit_sdf_read_golden.py")
        .current_dir(&root)
        .status()
        .expect("failed to run RDKit SDF read golden generator");
    assert!(
        status.success(),
        "RDKit SDF read golden generator failed with status {status}"
    );
}

fn load_golden() -> Vec<SdfReadRecord> {
    ensure_golden_exists();
    let path = golden_path();
    let file = File::open(&path).expect("should read SDF read golden");
    BufReader::new(file)
        .lines()
        .enumerate()
        .map(|(idx, line)| {
            let line = line.unwrap_or_else(|error| {
                panic!(
                    "failed to read {} line {}: {error}",
                    path.display(),
                    idx + 1
                )
            });
            serde_json::from_str(&line).unwrap_or_else(|error| {
                panic!(
                    "failed to parse {} line {}: {error}",
                    path.display(),
                    idx + 1
                )
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

fn bond_direction_name(direction: BondDirection) -> &'static str {
    match direction {
        BondDirection::None => "NONE",
        BondDirection::EndUpRight => "ENDUPRIGHT",
        BondDirection::EndDownRight => "ENDDOWNRIGHT",
    }
}

fn bond_stereo_name(stereo: BondStereo) -> &'static str {
    match stereo {
        BondStereo::None => "STEREONONE",
        BondStereo::Any => "STEREOANY",
        BondStereo::Cis => "STEREOZ",
        BondStereo::Trans => "STEREOE",
    }
}

fn chiral_tag_name(tag: ChiralTag) -> &'static str {
    match tag {
        ChiralTag::Unspecified => "CHI_UNSPECIFIED",
        ChiralTag::TetrahedralCw => "CHI_TETRAHEDRAL_CW",
        ChiralTag::TetrahedralCcw => "CHI_TETRAHEDRAL_CCW",
        ChiralTag::TrigonalBipyramidal => "CHI_TRIGONALBIPYRAMIDAL",
    }
}

fn parsed_record(record: &SdfReadRecord, row_idx: usize) -> cosmolkit_core::io::sdf::SdfRecord {
    let sdf = record
        .sdf
        .as_ref()
        .unwrap_or_else(|| panic!("rdkit_ok row {} has no SDF", row_idx + 1));
    read_sdf_record_from_str(sdf).unwrap_or_else(|error| {
        panic!(
            "COSMolKit should parse SDF row {} case {} {}: {error}",
            row_idx + 1,
            record.case_id,
            record.smiles
        )
    })
}

fn assert_case_matrix(records: &[SdfReadRecord]) {
    let expected_cases = [
        "2d_v2000_with_markers",
        "2d_v2000_coords_only",
        "2d_v3000_with_markers",
        "2d_v3000_coords_only",
        "3d_v2000_with_markers",
        "3d_v2000_coords_only",
        "3d_v3000_with_markers",
        "3d_v3000_coords_only",
    ];
    assert!(
        records.iter().any(|record| record.rdkit_ok),
        "SDF read golden should include at least one successful RDKit case"
    );
    for case_id in expected_cases {
        assert!(
            records.iter().any(|record| record.case_id == case_id),
            "SDF read golden is missing case {case_id}"
        );
    }
    for format in ["V2000", "V3000"] {
        assert!(
            records.iter().any(|record| record.format == format),
            "SDF read golden is missing format {format}"
        );
    }
    assert!(
        records.iter().any(|record| record.dimension == "3D"
            && record.stereo_markers == "coords_only"
            && record
                .chiral_tags
                .as_ref()
                .is_some_and(|tags| tags.iter().any(|tag| tag != "CHI_UNSPECIFIED"))),
        "SDF read golden must cover coordinate-inferred 3D chirality"
    );
    assert!(
        records
            .iter()
            .any(|record| record.dimension == "2D" && record.stereo_markers == "with_markers"),
        "SDF read golden must cover wedge/CFG stereo-marker parsing"
    );
}

#[test]
fn sdf_read_golden_covers_expected_case_matrix() {
    let records = load_golden();
    assert_case_matrix(&records);
}

#[test]
fn sdf_read_topology_and_atom_fields_match_rdkit() {
    let records = load_golden();
    for (row_idx, record) in records.iter().enumerate() {
        if !record.rdkit_ok {
            assert!(
                record.error.is_some(),
                "row {} {} ({}) is rdkit not ok but has no error",
                row_idx + 1,
                record.case_id,
                record.smiles
            );
            continue;
        }
        let parsed = parsed_record(record, row_idx);
        let expected_atoms = record.atoms.as_ref().expect("rdkit_ok row has atoms");
        let expected_bonds = record.bonds.as_ref().expect("rdkit_ok row has bonds");

        let actual_atoms = parsed
            .molecule
            .atoms()
            .iter()
            .map(|atom| AtomRecord {
                atomic_num: atom.atomic_num,
                isotope: atom.isotope,
                formal_charge: atom.formal_charge,
                is_aromatic: atom.is_aromatic,
                atom_map_num: atom.atom_map_num,
            })
            .collect::<Vec<_>>();
        assert_eq!(
            actual_atoms,
            *expected_atoms,
            "atom field mismatch at row {} case {} ({})",
            row_idx + 1,
            record.case_id,
            record.smiles
        );

        let actual_bonds = parsed
            .molecule
            .bonds()
            .iter()
            .map(|bond| BondRecord {
                begin: bond.begin_atom,
                end: bond.end_atom,
                bond_type: bond_order_name(bond.order).to_owned(),
                is_aromatic: bond.is_aromatic,
                direction: bond_direction_name(bond.direction).to_owned(),
                stereo: bond_stereo_name(bond.stereo).to_owned(),
                stereo_atoms: bond.stereo_atoms.clone(),
            })
            .collect::<Vec<_>>();
        assert_eq!(
            actual_bonds,
            *expected_bonds,
            "bond field mismatch at row {} case {} ({})",
            row_idx + 1,
            record.case_id,
            record.smiles
        );
    }

    let batch = batch_from_sdf_records(&records);
    for (row_idx, (record, batch_record)) in records.iter().zip(batch.records.iter()).enumerate() {
        if !record.rdkit_ok {
            continue;
        }
        let cosmolkit_core::BatchRecord::Valid(molecule) = batch_record else {
            panic!(
                "batch SDF read missing valid molecule at row {}",
                row_idx + 1
            );
        };
        let expected_atoms = record.atoms.as_ref().expect("rdkit_ok row has atoms");
        let actual_atoms = molecule
            .atoms()
            .iter()
            .map(|atom| AtomRecord {
                atomic_num: atom.atomic_num,
                isotope: atom.isotope,
                formal_charge: atom.formal_charge,
                is_aromatic: atom.is_aromatic,
                atom_map_num: atom.atom_map_num,
            })
            .collect::<Vec<_>>();
        assert_eq!(
            actual_atoms,
            *expected_atoms,
            "batch atom field mismatch at row {} case {} ({})",
            row_idx + 1,
            record.case_id,
            record.smiles
        );
        let expected_bonds = record.bonds.as_ref().expect("rdkit_ok row has bonds");
        let actual_bonds = molecule
            .bonds()
            .iter()
            .map(|bond| BondRecord {
                begin: bond.begin_atom,
                end: bond.end_atom,
                bond_type: bond_order_name(bond.order).to_owned(),
                is_aromatic: bond.is_aromatic,
                direction: bond_direction_name(bond.direction).to_owned(),
                stereo: bond_stereo_name(bond.stereo).to_owned(),
                stereo_atoms: bond.stereo_atoms.clone(),
            })
            .collect::<Vec<_>>();
        assert_eq!(
            actual_bonds,
            *expected_bonds,
            "batch bond field mismatch at row {} case {} ({})",
            row_idx + 1,
            record.case_id,
            record.smiles
        );
    }
}

#[test]
fn sdf_read_coordinates_match_rdkit_for_2d_and_3d_records() {
    let records = load_golden();
    for (row_idx, record) in records.iter().enumerate() {
        if !record.rdkit_ok {
            continue;
        }
        let parsed = parsed_record(record, row_idx);
        let expected_positions = record
            .positions
            .as_ref()
            .expect("rdkit_ok row should have positions");

        if record.dimension == "2D" {
            let coords = parsed.molecule.coords_2d().unwrap_or_else(|| {
                panic!(
                    "row {} {} should preserve 2D coords",
                    row_idx + 1,
                    record.case_id
                )
            });
            assert_eq!(coords.len(), expected_positions.len());
            for (atom_idx, (actual, expected)) in
                coords.iter().zip(expected_positions.iter()).enumerate()
            {
                assert!(
                    (actual.x - expected[0]).abs() <= 1e-12
                        && (actual.y - expected[1]).abs() <= 1e-12
                        && expected[2].abs() <= 1e-12,
                    "2D coordinate mismatch at row {} case {} atom {} ({})",
                    row_idx + 1,
                    record.case_id,
                    atom_idx,
                    record.smiles
                );
            }
        } else {
            let coords = parsed.molecule.coords_3d().unwrap_or_else(|| {
                panic!(
                    "row {} {} should preserve 3D coords",
                    row_idx + 1,
                    record.case_id
                )
            });
            assert_eq!(coords.len(), expected_positions.len());
            for (atom_idx, (actual, expected)) in
                coords.iter().zip(expected_positions.iter()).enumerate()
            {
                assert!(
                    (actual.x - expected[0]).abs() <= 1e-12
                        && (actual.y - expected[1]).abs() <= 1e-12
                        && (actual.z - expected[2]).abs() <= 1e-12,
                    "3D coordinate mismatch at row {} case {} atom {} ({})",
                    row_idx + 1,
                    record.case_id,
                    atom_idx,
                    record.smiles
                );
            }
        }
    }

    let batch = batch_from_sdf_records(&records);
    for (row_idx, (record, batch_record)) in records.iter().zip(batch.records.iter()).enumerate() {
        if !record.rdkit_ok {
            continue;
        }
        let cosmolkit_core::BatchRecord::Valid(molecule) = batch_record else {
            panic!(
                "batch SDF read missing valid molecule at row {}",
                row_idx + 1
            );
        };
        let expected_positions = record
            .positions
            .as_ref()
            .expect("rdkit_ok row should have positions");
        if record.dimension == "2D" {
            let coords = molecule.coords_2d().unwrap_or_else(|| {
                panic!(
                    "batch row {} {} should preserve 2D coords",
                    row_idx + 1,
                    record.case_id
                )
            });
            assert_eq!(coords.len(), expected_positions.len());
            for (atom_idx, (actual, expected)) in
                coords.iter().zip(expected_positions.iter()).enumerate()
            {
                assert!(
                    (actual.x - expected[0]).abs() <= 1e-12
                        && (actual.y - expected[1]).abs() <= 1e-12
                        && expected[2].abs() <= 1e-12,
                    "batch 2D coordinate mismatch at row {} case {} atom {} ({})",
                    row_idx + 1,
                    record.case_id,
                    atom_idx,
                    record.smiles
                );
            }
        } else {
            let coords = molecule.coords_3d().unwrap_or_else(|| {
                panic!(
                    "batch row {} {} should preserve 3D coords",
                    row_idx + 1,
                    record.case_id
                )
            });
            assert_eq!(coords.len(), expected_positions.len());
            for (atom_idx, (actual, expected)) in
                coords.iter().zip(expected_positions.iter()).enumerate()
            {
                assert!(
                    (actual.x - expected[0]).abs() <= 1e-12
                        && (actual.y - expected[1]).abs() <= 1e-12
                        && (actual.z - expected[2]).abs() <= 1e-12,
                    "batch 3D coordinate mismatch at row {} case {} atom {} ({})",
                    row_idx + 1,
                    record.case_id,
                    atom_idx,
                    record.smiles
                );
            }
        }
    }
}

#[test]
fn sdf_read_chirality_matches_rdkit_for_markers_and_coordinates() {
    let records = load_golden();
    for (row_idx, record) in records.iter().enumerate() {
        if !record.rdkit_ok {
            continue;
        }
        let parsed = parsed_record(record, row_idx);
        let expected_tags = record
            .chiral_tags
            .as_ref()
            .expect("rdkit_ok row should have chiral tags");
        let actual_tags = parsed
            .molecule
            .atoms()
            .iter()
            .map(|atom| chiral_tag_name(atom.chiral_tag).to_owned())
            .collect::<Vec<_>>();
        assert_eq!(
            actual_tags,
            *expected_tags,
            "chirality mismatch at row {} case {} ({})",
            row_idx + 1,
            record.case_id,
            record.smiles
        );
    }

    let batch = batch_from_sdf_records(&records);
    for (row_idx, (record, batch_record)) in records.iter().zip(batch.records.iter()).enumerate() {
        if !record.rdkit_ok {
            continue;
        }
        let cosmolkit_core::BatchRecord::Valid(molecule) = batch_record else {
            panic!(
                "batch SDF read missing valid molecule at row {}",
                row_idx + 1
            );
        };
        let expected_tags = record
            .chiral_tags
            .as_ref()
            .expect("rdkit_ok row should have chiral tags");
        let actual_tags = molecule
            .atoms()
            .iter()
            .map(|atom| chiral_tag_name(atom.chiral_tag).to_owned())
            .collect::<Vec<_>>();
        assert_eq!(
            actual_tags,
            *expected_tags,
            "batch chirality mismatch at row {} case {} ({})",
            row_idx + 1,
            record.case_id,
            record.smiles
        );
    }
}

#[test]
fn sdf_read_to_smiles_matches_rdkit_canonical_and_noncanonical() {
    let records = load_golden();
    for (row_idx, record) in records.iter().enumerate() {
        if !record.rdkit_ok {
            continue;
        }
        let parsed = parsed_record(record, row_idx);
        let expected = record
            .smiles_out
            .as_ref()
            .expect("rdkit_ok row should have SMILES output");

        let canonical = parsed
            .molecule
            .to_smiles_with_params(&SmilesWriteParams {
                canonical: true,
                do_isomeric_smiles: true,
                ..Default::default()
            })
            .unwrap_or_else(|error| {
                panic!(
                    "canonical SMILES write failed at row {} case {} ({}): {error}",
                    row_idx + 1,
                    record.case_id,
                    record.smiles
                )
            });
        assert_eq!(
            canonical,
            expected.canonical,
            "canonical SMILES mismatch at row {} case {} ({})",
            row_idx + 1,
            record.case_id,
            record.smiles
        );

        let noncanonical = parsed
            .molecule
            .to_smiles_with_params(&SmilesWriteParams {
                canonical: false,
                do_isomeric_smiles: true,
                ..Default::default()
            })
            .unwrap_or_else(|error| {
                panic!(
                    "noncanonical SMILES write failed at row {} case {} ({}): {error}",
                    row_idx + 1,
                    record.case_id,
                    record.smiles
                )
            });
        assert_eq!(
            noncanonical,
            expected.noncanonical,
            "noncanonical SMILES mismatch at row {} case {} ({})",
            row_idx + 1,
            record.case_id,
            record.smiles
        );
    }

    let batch = batch_from_sdf_records(&records);
    let canonical = batch
        .to_smiles_list_with_params(&SmilesWriteParams {
            canonical: true,
            do_isomeric_smiles: true,
            ..Default::default()
        })
        .expect("batch canonical SMILES should serialize");
    let noncanonical = batch
        .to_smiles_list_with_params(&SmilesWriteParams {
            canonical: false,
            do_isomeric_smiles: true,
            ..Default::default()
        })
        .expect("batch noncanonical SMILES should serialize");
    for (row_idx, ((record, canonical), noncanonical)) in
        records.iter().zip(canonical).zip(noncanonical).enumerate()
    {
        if !record.rdkit_ok {
            continue;
        }
        let expected = record
            .smiles_out
            .as_ref()
            .expect("rdkit_ok row should have SMILES output");
        assert_eq!(
            canonical.as_ref(),
            Some(&expected.canonical),
            "batch canonical SMILES mismatch at row {} case {} ({})",
            row_idx + 1,
            record.case_id,
            record.smiles
        );
        assert_eq!(
            noncanonical.as_ref(),
            Some(&expected.noncanonical),
            "batch noncanonical SMILES mismatch at row {} case {} ({})",
            row_idx + 1,
            record.case_id,
            record.smiles
        );
    }
}

fn batch_from_sdf_records(records: &[SdfReadRecord]) -> MoleculeBatch {
    let sdfs = records
        .iter()
        .map(|record| record.sdf.clone().unwrap_or_default())
        .collect::<Vec<_>>();
    MoleculeBatch::from_sdf_record_strings(&sdfs, SdfCoordinateMode::Auto, BatchErrorMode::Keep)
        .expect("batch SDF read should not raise in keep mode")
}
