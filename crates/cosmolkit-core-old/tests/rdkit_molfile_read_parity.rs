use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use cosmolkit_core::{
    BondDirection, BondOrder, BondStereo, ChiralTag, Molecule, SmilesWriteParams,
    io::{
        molfile::{MolFileRecord, read_mol_file_with_params, read_mol_record_from_str_with_params},
        sdf::SdfReadParams,
    },
};
use serde::Deserialize;

mod common;
use common::parity_data;

#[derive(Debug, Deserialize)]
struct MolFileReadRecord {
    smiles: Option<String>,
    case_id: String,
    dimension: String,
    format: String,
    stereo_markers: String,
    #[serde(default)]
    api: Option<String>,
    #[serde(default)]
    operation: Option<String>,
    #[serde(default)]
    sanitize: Option<bool>,
    #[serde(default)]
    remove_hs: Option<bool>,
    #[serde(default)]
    strict_parsing: Option<bool>,
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
    parity_data::golden_path("molfile_read.jsonl")
}

fn ensure_golden_exists() {
    let path = golden_path();
    assert!(
        path.exists(),
        "missing RDKit molfile read golden: {}. Generate it before running tests:\n\
         uv sync --group dev && {}",
        path.display(),
        parity_data::regenerate_command()
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
    order.rdkit_name()
}

fn record_label(record: &MolFileReadRecord) -> &str {
    record.smiles.as_deref().unwrap_or("<no-smiles>")
}

fn bond_direction_name(direction: BondDirection) -> &'static str {
    match direction {
        BondDirection::None => "NONE",
        BondDirection::BeginWedge => "BEGINWEDGE",
        BondDirection::BeginDash => "BEGINDASH",
        BondDirection::EndUpRight => "ENDUPRIGHT",
        BondDirection::EndDownRight => "ENDDOWNRIGHT",
        BondDirection::EitherDouble => "EITHERDOUBLE",
        BondDirection::Unknown => "UNKNOWN",
    }
}

fn bond_stereo_name(stereo: BondStereo) -> &'static str {
    match stereo {
        BondStereo::None => "STEREONONE",
        BondStereo::Any => "STEREOANY",
        BondStereo::Z => "STEREOZ",
        BondStereo::E => "STEREOE",
        BondStereo::Cis => "STEREOZ",
        BondStereo::Trans => "STEREOE",
        BondStereo::AtropCw => "STEREOATROPCW",
        BondStereo::AtropCcw => "STEREOATROPCCW",
    }
}

fn chiral_tag_name(tag: ChiralTag) -> &'static str {
    tag.rdkit_name()
}

fn parsed_record(record: &MolFileReadRecord, row_idx: usize) -> MolFileRecord {
    let molblock = record
        .molblock
        .as_ref()
        .unwrap_or_else(|| panic!("rdkit_ok row {} has no molblock", row_idx + 1));
    let params = SdfReadParams {
        sanitize: record.sanitize.unwrap_or(true),
        remove_hs: record.remove_hs.unwrap_or(false),
        strict_parsing: record.strict_parsing.unwrap_or(true),
        process_property_lists: false,
        ..Default::default()
    };
    let read_result = match record.api.as_deref().unwrap_or("MolFromMolBlock") {
        "MolFromMolBlock" => read_mol_record_from_str_with_params(molblock, params),
        "MolFromMolFile" => {
            let mut temp = tempfile::NamedTempFile::new().expect("should create temp molfile");
            std::io::Write::write_all(&mut temp, molblock.as_bytes())
                .expect("should write temp molfile");
            read_mol_file_with_params(temp.path(), params)
        }
        other => panic!(
            "unsupported molfile golden API {other} at row {}",
            row_idx + 1
        ),
    };
    let mut parsed = read_result.unwrap_or_else(|error| {
        panic!(
            "COSMolKit should parse molfile row {} case {} {}: {error}",
            row_idx + 1,
            record.case_id,
            record_label(record)
        )
    });
    parsed.molecule = apply_delayed_operation(parsed.molecule, record, row_idx);
    parsed
}

fn apply_delayed_operation(
    molecule: Molecule,
    record: &MolFileReadRecord,
    row_idx: usize,
) -> Molecule {
    match record.operation.as_deref().unwrap_or("read") {
        "read" => molecule,
        "delayed_sanitize" => molecule.sanitize().unwrap_or_else(|error| {
            panic!(
                "COSMolKit delayed sanitize failed at row {} case {} ({}): {error}",
                row_idx + 1,
                record.case_id,
                record_label(record)
            )
        }),
        "delayed_remove_hs" => molecule.without_hydrogens().unwrap_or_else(|error| {
            panic!(
                "COSMolKit delayed removeHs failed at row {} case {} ({}): {error}",
                row_idx + 1,
                record.case_id,
                record_label(record)
            )
        }),
        "failure" => molecule,
        other => panic!(
            "unsupported molfile golden operation {other} at row {}",
            row_idx + 1
        ),
    }
}

fn assert_error_matches_rdkit(record: &MolFileReadRecord, row_idx: usize) {
    let Some(molblock) = record.molblock.as_ref() else {
        assert!(
            record.error.is_some(),
            "row {} {} ({}) is rdkit not ok but has no error",
            row_idx + 1,
            record.case_id,
            record_label(record)
        );
        return;
    };
    let params = SdfReadParams {
        sanitize: record.sanitize.unwrap_or(true),
        remove_hs: record.remove_hs.unwrap_or(false),
        strict_parsing: record.strict_parsing.unwrap_or(true),
        process_property_lists: false,
        ..Default::default()
    };
    let result = match record.api.as_deref().unwrap_or("MolFromMolBlock") {
        "MolFromMolBlock" => read_mol_record_from_str_with_params(molblock, params),
        "MolFromMolFile" => {
            let mut temp = tempfile::NamedTempFile::new().expect("should create temp molfile");
            std::io::Write::write_all(&mut temp, molblock.as_bytes())
                .expect("should write temp molfile");
            read_mol_file_with_params(temp.path(), params)
        }
        other => panic!(
            "unsupported molfile golden API {other} at row {}",
            row_idx + 1
        ),
    };
    assert!(
        result.is_err(),
        "COSMolKit unexpectedly parsed RDKit failure row {} case {} ({})",
        row_idx + 1,
        record.case_id,
        record_label(record)
    );
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
            assert_error_matches_rdkit(record, row_idx);
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
                atomic_num: atom.atomic_number(),
                isotope: atom.isotope(),
                formal_charge: atom.formal_charge(),
                is_aromatic: atom.is_aromatic(),
                atom_map_num: atom.atom_map(),
            })
            .collect::<Vec<_>>();
        assert_eq!(
            actual_atoms,
            *expected_atoms,
            "atom field mismatch at row {} case {} ({})",
            row_idx + 1,
            record.case_id,
            record_label(record)
        );

        let actual_bonds = parsed
            .molecule
            .bonds()
            .iter()
            .map(|bond| BondRecord {
                begin: bond.begin().index(),
                end: bond.end().index(),
                bond_type: bond_order_name(bond.order()).to_owned(),
                is_aromatic: bond.is_aromatic(),
                direction: bond_direction_name(bond.direction()).to_owned(),
                stereo: bond_stereo_name(bond.stereo()).to_owned(),
                stereo_atoms: bond
                    .stereo_atoms()
                    .map(|atoms| atoms.into_iter().map(|atom| atom.index()).collect())
                    .unwrap_or_default(),
            })
            .collect::<Vec<_>>();
        assert_eq!(
            actual_bonds,
            *expected_bonds,
            "bond field mismatch at row {} case {} ({})",
            row_idx + 1,
            record.case_id,
            record_label(record)
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
        let Some(expected_positions) = record.positions.as_ref() else {
            continue;
        };

        if record.dimension == "2D" {
            let coords = parsed.molecule.coordinates_2d().unwrap_or_else(|| {
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
                    (actual[0] - expected[0]).abs() <= 1e-12
                        && (actual[1] - expected[1]).abs() <= 1e-12
                        && expected[2].abs() <= 1e-12,
                    "2D coordinate mismatch at row {} case {} atom {} ({})",
                    row_idx + 1,
                    record.case_id,
                    atom_idx,
                    record_label(record)
                );
            }
        } else {
            let coords = parsed
                .molecule
                .conformers_3d()
                .first()
                .map(|c| c.coordinates())
                .unwrap_or_else(|| {
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
                    (actual[0] - expected[0]).abs() <= 1e-12
                        && (actual[1] - expected[1]).abs() <= 1e-12
                        && (actual[2] - expected[2]).abs() <= 1e-12,
                    "3D coordinate mismatch at row {} case {} atom {} ({})",
                    row_idx + 1,
                    record.case_id,
                    atom_idx,
                    record_label(record)
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
            .map(|atom| chiral_tag_name(atom.chiral_tag()).to_owned())
            .collect::<Vec<_>>();
        assert_eq!(
            actual_tags,
            *expected_tags,
            "chirality mismatch at row {} case {} ({})",
            row_idx + 1,
            record.case_id,
            record_label(record)
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
                record_label(record)
            ),
            Err(error) => assert!(
                error.to_string().contains("unsupported path"),
                "canonical SMILES write failed without explicit unsupported error at row {} case {} ({}): {error}",
                row_idx + 1,
                record.case_id,
                record_label(record)
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
                    record_label(record)
                )
            });
        assert_eq!(
            noncanonical,
            expected.noncanonical,
            "noncanonical SMILES mismatch at row {} case {} ({})",
            row_idx + 1,
            record.case_id,
            record_label(record)
        );
    }
}
