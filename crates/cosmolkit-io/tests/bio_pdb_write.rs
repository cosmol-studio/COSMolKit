//! IO public text/file regressions for the coordinate PDB writer
//! (BIO-PDB-WRITE Steps 45-48).

use cosmolkit_io::{BioPdbWriteParams, bio_structure_to_pdb_text, write_bio_structure_pdb_file};

// Re-use the C02 fixture construction approach from the unit tests but
// through the PUBLIC detached API surface only.

fn c02_structure() -> cosmolkit_bio::BioStructureData {
    use cosmolkit_bio::*;
    use cosmolkit_types::Element;

    BioStructureData::from_parts(BioStructureParts {
        input_format: BioCoordinateFormat::Pdb,
        models: vec![BioModelRow::new(BioRowSpan::new(0, 1).unwrap(), Some(1))],
        chains: vec![BioChainRow::new(
            BioModelId::new(0),
            None,
            BioRowSpan::new(0, 1).unwrap(),
            ChainKind::default(),
            ChainSourceIds::new(Some(PdbChainId::from_ascii(b"A").unwrap()), None),
        )],
        residues: vec![BioResidueRow::new(
            BioChainId::new(0),
            BioRowSpan::new(0, 1).unwrap(),
            ResidueName::from_ascii(b"ALA").unwrap(),
            ResidueInfoKind::Unknown,
            EntityKind::Polymer,
            None,
            Some(b'A'),
            ResidueSourceIds::new(Some(PdbSeqId::new(1, None)), None, None, None, None).unwrap(),
            BioSiftsUnpResidue::default(),
        )],
        atoms: vec![BioAtomRow::new(
            BioResidueId::new(0),
            AtomName::from_ascii(b"CA").unwrap(),
            Element::C,
            None,
            None,
            0,
            BioCalcFlag::NotSet,
            1.0,
            20.0,
            [0.0; 6],
            -1,
            0.0,
            AtomSourceIds::new(None),
        )],
        entities: vec![],
        connections: vec![],
        cispeps: vec![],
        mod_residues: vec![],
        helices: vec![],
        sheets: vec![],
        metadata: BioMetadata::default(),
        source_state: BioStructureSourceState::default(),
        coordinates: BioCoordinateBlock::new(vec![[1.0, 2.0, 3.0]]),
        crystal: None,
        ncs_operators: vec![],
        assemblies: vec![],
    })
    .unwrap()
}

#[test]
fn bio_pdb_write_public_text_default() {
    let data = c02_structure();
    let text = bio_structure_to_pdb_text(&data, &BioPdbWriteParams::default())
        .expect("default params write C02");
    // C02 default: 1 ATOM + 1 TER + 1 END = 3×81 = 243 bytes
    assert_eq!(text.len(), 243, "C02 default output length");
    assert!(text.starts_with("ATOM      1  CA  ALA A   1"));
    assert!(text.contains("TER       2"));
    assert!(text.trim_end().ends_with("END"), "ends with padded END");
}

#[test]
fn bio_pdb_write_public_text_end_false() {
    let data = c02_structure();
    let params = BioPdbWriteParams {
        end_record: false,
        ..Default::default()
    };
    let text = bio_structure_to_pdb_text(&data, &params).expect("end=false writes C02");
    // Without END: 2×81 = 162 bytes
    assert_eq!(text.len(), 162, "end=false output length");
    assert!(!text.trim_end().ends_with("END"));
}

#[test]
fn bio_pdb_write_public_file_roundtrip() {
    let data = c02_structure();
    let dir = std::env::temp_dir().join("ck_pdb_write_test");
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("c02.pdb");
    write_bio_structure_pdb_file(&data, &path, &BioPdbWriteParams::default())
        .expect("write C02 to file");
    let read_back = std::fs::read(&path).expect("read back");
    let expected = bio_structure_to_pdb_text(&data, &BioPdbWriteParams::default()).unwrap();
    assert_eq!(read_back, expected.as_bytes(), "file bytes == text bytes");
    let _ = std::fs::remove_file(&path);
}

#[test]
fn bio_pdb_write_public_file_error_preserves_existing() {
    let dir = std::env::temp_dir().join("ck_pdb_write_test");
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("existing.pdb");
    std::fs::write(&path, b"EXISTING CONTENT").unwrap();

    // Attempt to write to a path that's a directory → IO error, existing
    // file must remain unchanged (the writer produces text first, then
    // fails on the IO, not before).
    let dir_path = dir.join("subdir.pdb");
    std::fs::create_dir_all(&dir_path).unwrap();
    let data = c02_structure();
    let result = write_bio_structure_pdb_file(&data, &dir_path, &BioPdbWriteParams::default());
    assert!(result.is_err(), "writing to a directory fails");

    // The EXISTING file in the same directory is untouched.
    assert_eq!(
        std::fs::read(&path).unwrap(),
        b"EXISTING CONTENT",
        "existing file unchanged"
    );
    let _ = std::fs::remove_file(&path);
    let _ = std::fs::remove_dir(&dir_path);
}
