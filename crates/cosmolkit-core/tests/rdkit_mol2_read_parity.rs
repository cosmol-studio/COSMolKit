use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use cosmolkit_core::{
    BondDirection, BondOrder, BondStereo, ChiralTag, SmilesWriteParams,
    io::mol2::{Mol2ReadParams, Mol2Type, read_mol2_from_str_with_params},
};
use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct Mol2ReadRecord {
    fixture: String,
    case_id: String,
    sanitize: bool,
    remove_hs: bool,
    cleanup_substructures: bool,
    variant: String,
    mol2: String,
    rdkit_ok: bool,
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
    chiral_tag: String,
    tripos_atom_name: Option<String>,
    tripos_atom_type: Option<String>,
    tripos_partial_charge: Option<String>,
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
    tripos_bond_type: Option<String>,
}

#[derive(Debug, Deserialize)]
struct SmilesOut {
    canonical: String,
    noncanonical: String,
}

fn golden_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/golden/mol2_read.jsonl")
}

fn ensure_golden_exists() {
    let path = golden_path();
    assert!(
        path.exists(),
        "missing RDKit MOL2 read golden: {}. Generate it before running tests:\n\
         uv sync --group dev && .venv/bin/python tests/scripts/gen_all_rdkit_goldens.py --python .venv/bin/python --clean --jobs 4",
        path.display()
    );
}

fn load_golden() -> Vec<Mol2ReadRecord> {
    ensure_golden_exists();
    let path = golden_path();
    let file = File::open(&path).expect("should read MOL2 read golden");
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

fn read_params(record: &Mol2ReadRecord) -> Mol2ReadParams {
    assert_eq!(
        record.variant, "CORINA",
        "MOL2 parity currently covers RDKit CORINA variant only"
    );
    Mol2ReadParams {
        sanitize: record.sanitize,
        remove_hs: record.remove_hs,
        variant: Mol2Type::Corina,
        cleanup_substructures: record.cleanup_substructures,
    }
}

fn parsed_record(
    record: &Mol2ReadRecord,
    row_idx: usize,
) -> Option<cosmolkit_core::io::mol2::Mol2Record> {
    read_mol2_from_str_with_params(&record.mol2, read_params(record)).unwrap_or_else(|error| {
        panic!(
            "COSMolKit should parse MOL2 row {} fixture {} case {}: {error}",
            row_idx + 1,
            record.fixture,
            record.case_id
        )
    })
}

fn assert_case_matrix(records: &[Mol2ReadRecord]) {
    assert!(
        records.iter().any(|record| record.rdkit_ok),
        "MOL2 read golden should include at least one successful RDKit case"
    );
    for fixture in [
        "pyrazole_pyridine.mol2",
        "benzene.mol2",
        "mol_noatoms.mol2",
        "mol_nomol.mol2",
        "lonePairMol.mol2",
        "symmetricGuanidine.mol2",
        "Noxide.mol2",
        "fusedRing.mol2",
        "EZ_mol2_issue114.mol2",
    ] {
        assert!(
            records.iter().any(|record| record.fixture == fixture),
            "MOL2 read golden is missing fixture {fixture}"
        );
    }
    for case_id in [
        "sanitize_true_remove_hs_true_cleanup_true",
        "sanitize_true_remove_hs_false_cleanup_true",
        "sanitize_false_remove_hs_false_cleanup_true",
        "sanitize_false_remove_hs_false_cleanup_false",
    ] {
        assert!(
            records.iter().any(|record| record.case_id == case_id),
            "MOL2 read golden is missing case {case_id}"
        );
    }
    assert!(
        records
            .iter()
            .any(|record| record.rdkit_ok && record.positions.is_some()),
        "MOL2 read golden must cover 3D coordinate preservation"
    );
}

#[test]
fn mol2_read_golden_covers_expected_case_matrix() {
    let records = load_golden();
    assert_case_matrix(&records);
}

#[test]
fn mol2_read_topology_atom_and_bond_fields_match_rdkit() {
    let records = load_golden();
    for (row_idx, record) in records.iter().enumerate() {
        if !record.rdkit_ok {
            assert!(
                record.error.is_some(),
                "row {} {} ({}) is rdkit not ok but has no error",
                row_idx + 1,
                record.fixture,
                record.case_id
            );
            continue;
        }
        let parsed = parsed_record(record, row_idx).unwrap_or_else(|| {
            panic!(
                "COSMolKit returned null for RDKit-ok row {} fixture {} case {}",
                row_idx + 1,
                record.fixture,
                record.case_id
            )
        });
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
                chiral_tag: chiral_tag_name(atom.chiral_tag()).to_owned(),
                tripos_atom_name: atom.prop("_TriposAtomName").map(ToOwned::to_owned),
                tripos_atom_type: atom.prop("_TriposAtomType").map(ToOwned::to_owned),
                tripos_partial_charge: atom.prop("_TriposPartialCharge").map(ToOwned::to_owned),
            })
            .collect::<Vec<_>>();
        assert_eq!(
            actual_atoms,
            *expected_atoms,
            "atom field mismatch at row {} fixture {} case {}",
            row_idx + 1,
            record.fixture,
            record.case_id
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
                tripos_bond_type: bond.prop("_TriposBondType").map(ToOwned::to_owned),
            })
            .collect::<Vec<_>>();
        assert_eq!(
            actual_bonds,
            *expected_bonds,
            "bond field mismatch at row {} fixture {} case {}",
            row_idx + 1,
            record.fixture,
            record.case_id
        );
    }
}

#[test]
fn mol2_read_coordinates_and_chirality_match_rdkit() {
    let records = load_golden();
    for (row_idx, record) in records.iter().enumerate() {
        if !record.rdkit_ok {
            continue;
        }
        let parsed = parsed_record(record, row_idx).unwrap_or_else(|| {
            panic!(
                "COSMolKit returned null for RDKit-ok row {} fixture {} case {}",
                row_idx + 1,
                record.fixture,
                record.case_id
            )
        });
        let expected_positions = record
            .positions
            .as_ref()
            .expect("rdkit_ok row should have positions");
        let coords = parsed
            .molecule
            .conformers_3d()
            .first()
            .map(|c| c.coordinates())
            .unwrap_or_else(|| {
                panic!(
                    "row {} fixture {} case {} should preserve 3D coords",
                    row_idx + 1,
                    record.fixture,
                    record.case_id
                )
            });
        assert_eq!(coords.len(), expected_positions.len());
        for (atom_idx, (actual, expected)) in coords.iter().zip(expected_positions).enumerate() {
            assert!(
                (actual[0] - expected[0]).abs() <= 1e-12
                    && (actual[1] - expected[1]).abs() <= 1e-12
                    && (actual[2] - expected[2]).abs() <= 1e-12,
                "3D coordinate mismatch at row {} fixture {} case {} atom {}",
                row_idx + 1,
                record.fixture,
                record.case_id,
                atom_idx
            );
        }

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
            "chirality mismatch at row {} fixture {} case {}",
            row_idx + 1,
            record.fixture,
            record.case_id
        );
    }
}

#[test]
fn mol2_read_to_smiles_matches_rdkit_canonical_and_noncanonical() {
    let records = load_golden();
    for (row_idx, record) in records.iter().enumerate() {
        if !record.rdkit_ok {
            continue;
        }
        let parsed = parsed_record(record, row_idx).unwrap_or_else(|| {
            panic!(
                "COSMolKit returned null for RDKit-ok row {} fixture {} case {}",
                row_idx + 1,
                record.fixture,
                record.case_id
            )
        });
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
                "canonical SMILES mismatch at row {} fixture {} case {}",
                row_idx + 1,
                record.fixture,
                record.case_id
            ),
            Err(error) => assert!(
                error.to_string().contains("unsupported path"),
                "canonical SMILES write failed without explicit unsupported error at row {} fixture {} case {}: {error}",
                row_idx + 1,
                record.fixture,
                record.case_id
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
                    "noncanonical SMILES write failed at row {} fixture {} case {}: {error}",
                    row_idx + 1,
                    record.fixture,
                    record.case_id
                )
            });
        assert_eq!(
            noncanonical,
            expected.noncanonical,
            "noncanonical SMILES mismatch at row {} fixture {} case {}",
            row_idx + 1,
            record.fixture,
            record.case_id
        );
    }
}
