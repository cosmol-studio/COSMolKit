//! Public facade tests for BioStructure PDB output (BIO-PDB-WRITE Steps 51-54).

use cosmolkit::{BioPdbWriteParams, BioStructure};

fn pdb_text() -> &'static str {
    "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \nTER       2      ALA A   1                                                      \nEND"
}

fn mmcif_text() -> &'static str {
    r#"data_test
_entry.id test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
ATOM 1 C CA ALA A 1 1.000 2.000 3.000 1.00 20.00
"#
}

fn all_pdb_controls() -> BioPdbWriteParams {
    BioPdbWriteParams {
        ter_records: true,
        numbered_ter: true,
        ter_ignores_type: true,
        preserve_serial: true,
        end_record: true,
    }
}

#[test]
fn bio_pdb_output_pdb_all_controls_smoke() {
    let structure = BioStructure::from_pdb(pdb_text()).expect("parses PDB");
    let text = structure
        .to_pdb_with_params(&all_pdb_controls())
        .expect("writes PDB with all controls enabled");
    assert!(!text.is_empty(), "PDB output is nonempty");
    assert!(text.contains("ATOM"), "contains ATOM");
    assert!(text.contains("TER"), "contains TER");
}

#[test]
fn bio_pdb_output_mmcif_all_controls_smoke() {
    let structure = BioStructure::from_mmcif(mmcif_text()).expect("parses mmCIF");
    let text = structure
        .to_pdb_with_params(&all_pdb_controls())
        .expect("writes PDB from mmCIF with all controls enabled");
    assert!(!text.is_empty(), "mmCIF-to-PDB output is nonempty");
}

#[test]
fn bio_pdb_output_pdb_explicit_ter_controls_smoke() {
    let structure = BioStructure::from_pdb(pdb_text()).expect("parses PDB");
    let text = structure
        .to_pdb_with_params(&BioPdbWriteParams {
            ter_records: true,
            numbered_ter: true,
            ter_ignores_type: false,
            preserve_serial: false,
            end_record: true,
        })
        .expect("writes PDB with explicit TER controls");
    assert!(!text.is_empty(), "PDB output is nonempty");
    assert!(text.contains("ATOM"), "contains ATOM");
    assert!(text.contains("TER"), "contains TER");
}

#[test]
fn bio_pdb_output_to_pdb_default() {
    let structure = BioStructure::from_pdb(pdb_text()).expect("parses");
    let text = structure.to_pdb().expect("writes");
    assert!(!text.is_empty());
    assert!(text.contains("ATOM"));
    assert!(text.contains("TER"));
}

#[test]
fn bio_pdb_output_to_pdb_explicit_equality() {
    let structure = BioStructure::from_pdb(pdb_text()).expect("parses");
    let default = structure.to_pdb().unwrap();
    let explicit = structure
        .to_pdb_with_params(&BioPdbWriteParams::default())
        .unwrap();
    assert_eq!(default, explicit, "default == explicit default params");
}

#[test]
fn bio_pdb_output_to_pdb_end_false() {
    let structure = BioStructure::from_pdb(pdb_text()).expect("parses");
    let params = BioPdbWriteParams {
        end_record: false,
        ..Default::default()
    };
    let text = structure.to_pdb_with_params(&params).expect("writes");
    let default = structure.to_pdb().unwrap();
    assert!(text.len() < default.len(), "end=false shorter than default");
}

#[test]
fn bio_pdb_output_write_pdb_file() {
    let structure = BioStructure::from_pdb(pdb_text()).expect("parses");
    let dir = std::env::temp_dir().join("ck_pdb_facade_test");
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("test.pdb");
    structure.write_pdb(&path).expect("writes file");
    let bytes = std::fs::read(&path).expect("reads back");
    let text = structure.to_pdb().unwrap();
    assert_eq!(bytes, text.as_bytes(), "file == text");
    let _ = std::fs::remove_file(&path);
}

#[test]
fn bio_pdb_output_write_pdb_file_with_params() {
    let structure = BioStructure::from_pdb(pdb_text()).expect("parses");
    let dir = std::env::temp_dir().join("ck_pdb_facade_test");
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("test_params.pdb");
    structure
        .write_pdb_with_params(&path, &BioPdbWriteParams::default())
        .expect("writes file with params");
    let _ = std::fs::remove_file(&path);
}

#[test]
fn bio_pdb_output_write_pdb_directory_fails() {
    let structure = BioStructure::from_pdb(pdb_text()).expect("parses");
    let dir = std::env::temp_dir().join("ck_pdb_facade_test_dir_target");
    std::fs::create_dir_all(&dir).unwrap();
    let result = structure.write_pdb(&dir);
    assert!(result.is_err(), "writing to a directory fails");
    let _ = std::fs::remove_dir(&dir);
}

#[test]
fn bio_pdb_output_write_pdb_missing_dir_fails() {
    let structure = BioStructure::from_pdb(pdb_text()).expect("parses");
    let path = std::env::temp_dir().join("ck_pdb_facade_test/nonexistent/test.pdb");
    let result = structure.write_pdb(&path);
    assert!(result.is_err(), "writing to missing directory fails");
}

#[test]
fn bio_pdb_output_error_preserves_existing_file() {
    let dir = std::env::temp_dir().join("ck_pdb_facade_test");
    std::fs::create_dir_all(&dir).unwrap();
    let existing = dir.join("existing.pdb");
    std::fs::write(&existing, b"EXISTING").unwrap();

    // Attempt to write to a path that is a directory.
    let dir_target = dir.join("target_dir.pdb");
    std::fs::create_dir_all(&dir_target).unwrap();
    let structure = BioStructure::from_pdb(pdb_text()).expect("parses");
    let _ = structure.write_pdb(&dir_target); // Fails.

    assert_eq!(
        std::fs::read(&existing).unwrap(),
        b"EXISTING",
        "existing file unchanged after failed write"
    );
    let _ = std::fs::remove_file(&existing);
    let _ = std::fs::remove_dir(&dir_target);
}

/// The FROZEN 24 public-call census (§13.4).
#[test]
fn bio_pdb_output_frozen24_public_calls() {
    let dir = std::env::temp_dir().join("ck_pdb_facade_test");
    std::fs::create_dir_all(&dir).unwrap();
    let mut calls = 0usize;
    let structure = BioStructure::from_pdb(pdb_text()).expect("parses");
    let profiles = [
        BioPdbWriteParams::default(),
        BioPdbWriteParams {
            end_record: false,
            ..Default::default()
        },
        BioPdbWriteParams {
            ter_records: false,
            ..Default::default()
        },
        BioPdbWriteParams {
            preserve_serial: true,
            ..Default::default()
        },
    ];
    for (i, params) in profiles.iter().enumerate() {
        let _ = structure.to_pdb();
        calls += 1;
        let _ = structure.to_pdb_with_params(params);
        calls += 1;
        let p = dir.join(format!("f24_{i}.pdb"));
        let _ = structure.write_pdb(&p);
        calls += 1;
        let _ = std::fs::remove_file(&p);
        let p2 = dir.join(format!("f24p_{i}.pdb"));
        let _ = structure.write_pdb_with_params(&p2, params);
        calls += 1;
        let _ = std::fs::remove_file(&p2);
    }
    let s2 = BioStructure::from_pdb(pdb_text()).expect("re-parses");
    for (i, s) in [&structure, &s2].iter().enumerate() {
        let _ = s.to_pdb();
        calls += 1;
        let _ = s.to_pdb_with_params(&BioPdbWriteParams::default());
        calls += 1;
        let p = dir.join(format!("prov_{i}.pdb"));
        let _ = s.write_pdb(&p);
        calls += 1;
        let _ = std::fs::remove_file(&p);
        let p2 = dir.join(format!("provp_{i}.pdb"));
        let _ = s.write_pdb_with_params(&p2, &BioPdbWriteParams::default());
        calls += 1;
        let _ = std::fs::remove_file(&p2);
    }
    assert_eq!(calls, 24, "exact 24 frozen public calls");
}

/// The FROZEN 8 file-control census (§13.4).
#[test]
fn bio_pdb_output_frozen8_file_controls() {
    let dir = std::env::temp_dir().join("ck_pdb_facade_test");
    std::fs::create_dir_all(&dir).unwrap();
    let structure = BioStructure::from_pdb(pdb_text()).expect("parses");
    let mut controls = 0usize;

    let p1 = dir.join("c1.pdb");
    structure.write_pdb(&p1).expect("c1");
    controls += 1;
    let _ = std::fs::remove_file(&p1);

    let p2 = dir.join("c2.pdb");
    std::fs::write(&p2, b"EXISTING").unwrap();
    let d2 = dir.join("c2_dir");
    std::fs::create_dir_all(&d2).unwrap();
    let _ = structure.write_pdb(&d2);
    controls += 1;
    assert_eq!(std::fs::read(&p2).unwrap(), b"EXISTING");
    let _ = std::fs::remove_file(&p2);
    let _ = std::fs::remove_dir(&d2);

    let p3 = dir.join("nonexistent/c3.pdb");
    assert!(structure.write_pdb(&p3).is_err(), "c3");
    controls += 1;

    let p4 = dir.join("c4_dir");
    std::fs::create_dir_all(&p4).unwrap();
    assert!(structure.write_pdb(&p4).is_err(), "c4");
    controls += 1;
    let _ = std::fs::remove_dir(&p4);

    let p5 = dir.join("c5.pdb");
    structure
        .write_pdb_with_params(&p5, &BioPdbWriteParams::default())
        .expect("c5");
    controls += 1;
    let _ = std::fs::remove_file(&p5);

    let p6 = dir.join("nonexistent/c6.pdb");
    assert!(
        structure
            .write_pdb_with_params(&p6, &BioPdbWriteParams::default())
            .is_err(),
        "c6"
    );
    controls += 1;

    let p7 = dir.join("c7.pdb");
    structure.write_pdb(&p7).expect("c7");
    let file = std::fs::read(&p7).unwrap();
    let text = structure.to_pdb().unwrap();
    controls += 1;
    assert_eq!(file, text.as_bytes(), "c7 file==text");
    let _ = std::fs::remove_file(&p7);

    let p8 = dir.join("c8.pdb");
    structure
        .write_pdb_with_params(
            &p8,
            &BioPdbWriteParams {
                end_record: false,
                ..Default::default()
            },
        )
        .expect("c8");
    controls += 1;
    let _ = std::fs::remove_file(&p8);

    assert_eq!(controls, 8, "exact 8 frozen file controls");
}
