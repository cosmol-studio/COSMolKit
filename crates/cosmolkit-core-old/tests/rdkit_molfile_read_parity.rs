use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use cosmolkit_core::{
    BondDirection, BondOrder, BondStereo, ChiralTag, SmilesWriteParams,
    io::molfile::read_mol_record_from_str,
};
use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct MolFileReadRecord {
    smiles: String,
    case_id: String,
    dimension: String,
    format: String,
    stereo_markers: String,
    rdkit_ok: bool,
    molblock: Option<String>,
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
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/molfile_read.jsonl")
}

fn ensure_golden_exists() {
    let path = golden_path();
    assert!(
        path.exists(),
        "missing RDKit molfile read golden: {}. Generate it before running tests:\n\
         uv sync --group dev && .venv/bin/python tests/scripts/gen_rdkit_molfile_read_golden.py --input tests/smiles.smi --output tests/golden/molfile_read.jsonl",
        path.display()
    );
}

fn load_golden() -> Vec<MolFileReadRecord> {
    ensure_golden_exists();
    let path = golden_path();
    let file = File::open(&path).expect("should read molfile read golden");
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
        BondOrder::Hydrogen => "HYDROGEN",
        BondOrder::Null => "UNSPECIFIED",
    }
}

fn bond_direction_name(direction: BondDirection) -> &'static str {
    match direction {
        BondDirection::None => "NONE",
        BondDirection::EndUpRight => "ENDUPRIGHT",
        BondDirection::EndDownRight => "ENDDOWNRIGHT",
        BondDirection::Unknown => "UNKNOWN",
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

fn parsed_record(
    record: &MolFileReadRecord,
    row_idx: usize,
) -> cosmolkit_core::io::molfile::MolFileRecord {
    let molblock = record
        .molblock
        .as_ref()
        .unwrap_or_else(|| panic!("rdkit_ok row {} has no molblock", row_idx + 1));
    read_mol_record_from_str(molblock).unwrap_or_else(|error| {
        panic!(
            "COSMolKit should parse molfile row {} case {} {}: {error}",
            row_idx + 1,
            record.case_id,
            record.smiles
        )
    })
}

fn assert_case_matrix(records: &[MolFileReadRecord]) {
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
        "molfile read golden should include at least one successful RDKit case"
    );
    for case_id in expected_cases {
        assert!(
            records.iter().any(|record| record.case_id == case_id),
            "molfile read golden is missing case {case_id}"
        );
    }
    for format in ["V2000", "V3000"] {
        assert!(
            records.iter().any(|record| record.format == format),
            "molfile read golden is missing format {format}"
        );
    }
    assert!(
        records.iter().any(|record| record.dimension == "3D"
            && record.stereo_markers == "coords_only"
            && record
                .chiral_tags
                .as_ref()
                .is_some_and(|tags| tags.iter().any(|tag| tag != "CHI_UNSPECIFIED"))),
        "molfile read golden must cover coordinate-inferred 3D chirality"
    );
    assert!(
        records
            .iter()
            .any(|record| record.dimension == "2D" && record.stereo_markers == "with_markers"),
        "molfile read golden must cover wedge/CFG stereo-marker parsing"
    );
}

#[test]
fn molfile_read_golden_covers_expected_case_matrix() {
    let records = load_golden();
    assert_case_matrix(&records);
}

#[test]
fn molfile_read_topology_and_atom_fields_match_rdkit() {
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
}

#[test]
fn molfile_read_coordinates_match_rdkit_for_2d_and_3d_records() {
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
}

#[test]
fn molfile_read_chirality_matches_rdkit_for_markers_and_coordinates() {
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
}

#[test]
fn molfile_read_to_smiles_matches_rdkit_canonical_and_noncanonical() {
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

        match parsed.molecule.to_smiles_with_params(&SmilesWriteParams {
            canonical: true,
            do_isomeric_smiles: true,
            ..Default::default()
        }) {
            Ok(canonical) => assert_eq!(
                canonical,
                expected.canonical,
                "canonical SMILES mismatch at row {} case {} ({})",
                row_idx + 1,
                record.case_id,
                record.smiles
            ),
            Err(error) => assert!(
                error.to_string().contains("unsupported path"),
                "canonical SMILES write failed without explicit unsupported error at row {} case {} ({}): {error}",
                row_idx + 1,
                record.case_id,
                record.smiles
            ),
        }

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
}
