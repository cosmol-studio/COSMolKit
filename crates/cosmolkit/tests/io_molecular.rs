#![cfg(feature = "cap-io")]
//! Proposed fixed source-boundary regressions; independent acceptance pending.
use cosmolkit::{
    Mol2ReadParams, MolecularIoError, Molecule, PropertyText, PropertyValue, SdfDataset, SdfError,
    SdfReadParams, SdfReader, SdfRecord, SdfRecordStream,
};
use std::io::Cursor;

const MOL: &str = "carbon\n  RDKit          2D\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    1.2500   -2.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";
const MOL2: &str = "@<TRIPOS>MOLECULE\ncarbon\n1 0 0 0 0\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n1 C1 1.25 -2.5 0.0 C.3 1 MOL 0.0\n@<TRIPOS>BOND\n";

struct FileFixture(std::path::PathBuf);
impl FileFixture {
    fn new(text: &str) -> Self {
        let id = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path =
            std::env::temp_dir().join(format!("cosmolkit-io44-{}-{id}.sdf", std::process::id()));
        std::fs::write(&path, text).unwrap();
        Self(path)
    }
    fn path(&self) -> &str {
        self.0.to_str().unwrap()
    }
}
impl Drop for FileFixture {
    fn drop(&mut self) {
        std::fs::remove_file(&self.0).unwrap();
    }
}

#[test]
fn xyz_keeps_disconnected_graph_coordinates_and_comment() {
    let input = "3\noriginal comment\nO 0 0 0\nH 0 0 1\nH 1 0 0\n";
    let molecule = Molecule::from_xyz_block(input).unwrap();
    assert_eq!(molecule.num_atoms(), 3);
    assert_eq!(molecule.num_bonds(), 0);
    assert_eq!(molecule.conformers_3d().len(), 1);
    let written = molecule.to_xyz().unwrap();
    let reread = Molecule::from_xyz_block(&written).unwrap();
    assert_eq!(reread.num_atoms(), 3);
    assert_eq!(reread.num_bonds(), 0);
    assert_eq!(reread.conformers_3d(), molecule.conformers_3d());
    assert_eq!(
        molecule.properties().prop("_FileComments"),
        Some(&PropertyValue::from("original comment"))
    );
    assert_eq!(molecule.properties().name(), None);
    assert_eq!(written.lines().nth(1), Some(""));
}

#[test]
fn mol_ignores_unread_sdf_fields_but_sdf_preserves_duplicate_order() {
    let text = format!("{MOL}>  <ID>\nfirst\n\n>  <ID>\nsecond\n\n$$$$\n");
    let mol = Molecule::from_mol(&text).unwrap();
    assert_eq!(mol.num_atoms(), 1);
    assert!(mol.properties().sdf_data_fields().is_empty());
    let record = SdfRecord::from_sdf(&text).unwrap();
    assert_eq!(
        record.data_fields(),
        &[
            ("ID".into(), "first".into()),
            ("ID".into(), "second".into())
        ]
    );
    assert_eq!(record.data_field("ID"), Some(&PropertyText::from("first")));
    assert_eq!(record.data_field("absent"), None);
    assert_eq!(record.title(), Some(&PropertyText::from("carbon")));
    assert_eq!(record.index(), 0);
}

#[test]
fn mol_does_not_interpret_a_mismatched_trailing_atom_property_list() {
    let text = format!("{MOL}>  <atom.iprop.count>\n1 2\n\n$$$$\n");
    assert_eq!(Molecule::from_mol(&text).unwrap().num_atoms(), 1);
    assert!(matches!(SdfRecord::from_sdf(&text), Err(SdfError::Read(
        cosmolkit::SdfReadError::PropertyListCount { target: "atom", name, actual: 2, expected: 1 }
    )) if name == "atom.iprop.count"));
    let params = SdfReadParams {
        strict_parsing: false,
        ..Default::default()
    };
    let record = SdfRecord::from_sdf_with_params(&text, &params).unwrap();
    assert_eq!(
        record.data_field("atom.iprop.count"),
        Some(&PropertyText::from("1 2"))
    );
}

#[test]
fn forward_reader_consumes_bad_record_then_recovers_with_original_index() {
    let bad = "broken\n$$$$\n";
    let good = format!("{MOL}$$$$\n");
    let original_text = format!("{bad}{good}");
    let mut original_reader = SdfRecordStream::new(Cursor::new(original_text.as_bytes()));
    assert!(original_reader.next_record().is_err());
    assert_eq!(original_reader.bytes_consumed(), original_text.len() as u64);
    assert_eq!(original_reader.records_consumed(), 1);
    assert!(original_reader.next_record().unwrap().is_none());
    assert!(original_reader.is_end());
    // The parser consumes header/count lines before recovery; the delimiter
    // in the original short malformed header is not a recovery boundary.
    // A third physical record remains available after the next delimiter.
    let text = format!("{bad}{good}{good}");
    let mut reader = SdfRecordStream::new(Cursor::new(text.as_bytes()));
    assert_eq!(reader.records_consumed(), 0);
    assert_eq!(reader.bytes_consumed(), 0);
    assert!(reader.next_record().is_err());
    assert_eq!(reader.records_consumed(), 1);
    assert_eq!(reader.bytes_consumed(), (bad.len() + good.len()) as u64);
    let record = reader.next_record().unwrap().unwrap();
    assert_eq!(record.index(), 1);
    assert_eq!(record.title(), Some(&PropertyText::from("carbon")));
    assert_eq!(reader.bytes_consumed(), text.len() as u64);
    assert_eq!(reader.records_consumed(), 2);
    assert!(reader.next_record().unwrap().is_none());
    assert!(reader.is_end());
    assert!(reader.next().is_none());
}

#[test]
fn dataset_indexes_without_parsing_and_retains_half_open_byte_line_ranges() {
    let bad = "broken\n$$$$\n";
    let good = format!("{MOL}>  <NAME>\nα\n\n$$$$\n");
    let text = format!("{bad}{good}");
    let file = FileFixture::new(&text);
    let dataset = SdfDataset::open(file.path()).unwrap();
    assert_eq!(dataset.len(), 2);
    assert!(!dataset.is_empty());
    assert_eq!(dataset.path(), file.0.as_path());
    assert!(dataset.record(0).is_err());
    let metadata = dataset.metadata(1).unwrap();
    assert_eq!(metadata.index(), 1);
    assert_eq!(metadata.byte_offset(), bad.len() as u64);
    assert_eq!(metadata.byte_len(), good.len() as u64);
    assert_eq!(metadata.byte_range(), (bad.len() as u64, text.len() as u64));
    assert_eq!(metadata.line_range(), (2, text.lines().count()));
    assert_eq!(metadata.title(), Some("carbon"));
    assert_eq!(dataset.record_text(1).unwrap(), good);
    assert!(dataset.metadata(2).is_none());
    assert!(dataset.record(2).is_err());
    let mut iter = dataset.iter();
    assert_eq!(iter.len(), 2);
    assert!(iter.next().unwrap().is_err());
    assert_eq!(iter.len(), 1);
    assert_eq!(iter.next().unwrap().unwrap().index(), 1);
    assert_eq!(iter.len(), 0);
    assert!(iter.next().is_none());
    assert!(iter.next().is_none());
}

#[test]
fn molecule_file_readers_preserve_first_record_and_real_io_source() {
    let file = FileFixture::new(&format!("{MOL}$$$$\nbroken\n$$$$\n"));
    assert_eq!(Molecule::read_mol(file.path()).unwrap().num_atoms(), 1);
    assert_eq!(Molecule::read_sdf(file.path()).unwrap().num_atoms(), 1);
    let missing = format!("{}.missing", file.path());
    let error = Molecule::read_mol(&missing).unwrap_err();
    assert!(
        matches!(&error,MolecularIoError::Io{path,source} if path.to_str()==Some(&missing) && source.kind()==std::io::ErrorKind::NotFound)
    );
    assert!(std::error::Error::source(&error).is_some());
    let empty = FileFixture::new("");
    assert!(matches!(
        Molecule::read_sdf(empty.path()),
        Err(MolecularIoError::NoRecord { format: "SDF" })
    ));
}

#[test]
fn mol2_finalization_runs_all_public_sanitize_remove_hydrogens_branches() {
    for sanitize in [false, true] {
        for remove_hs in [false, true] {
            for cleanup_substructures in [false, true] {
                let params = Mol2ReadParams {
                    sanitize,
                    remove_hs,
                    cleanup_substructures,
                    ..Default::default()
                };
                let molecule = Molecule::from_mol2_with_params(MOL2, &params).unwrap();
                assert_eq!(molecule.num_atoms(), 1);
                assert_eq!(molecule.conformers_3d().len(), 1);
                assert_eq!(
                    molecule.properties().name(),
                    Some(&PropertyText::from("carbon"))
                );
            }
        }
    }
    assert!(
        matches!(Molecule::from_mol2(""),Err(MolecularIoError::Mol2Read(cosmolkit::Mol2ReadError::Parse(message))) if message=="No MOLECULE block found in Mol2 data")
    );
}

#[test]
fn file_reader_source_does_not_open_the_file_before_batches() {
    let path = std::env::temp_dir().join(format!(
        "cosmolkit-io44-missing-reader-source-{}.sdf",
        std::process::id()
    ));
    assert!(!path.exists());
    let source = SdfReader::open(path.to_str().unwrap()).unwrap();
    assert_eq!(source.path(), path.as_path());
    assert_eq!(source.params(), &SdfReadParams::default());
}

#[test]
fn io44_empty_mol_title_survives_record_live_and_batch_boundaries() {
    let text = MOL.replacen("carbon", "", 1);
    for sanitize in [false, true] {
        let params = SdfReadParams {
            sanitize,
            ..Default::default()
        };
        let molecule = Molecule::from_mol_with_params(&text, &params).unwrap();
        assert_eq!(molecule.properties().name(), Some(&PropertyText::from("")));
        let sdf = format!("{text}$$$$\n");
        let record = SdfRecord::from_sdf_with_params(&sdf, &params).unwrap();
        assert_eq!(record.title(), Some(&PropertyText::from("")));
        assert_eq!(
            record.molecule().unwrap().properties().name(),
            Some(&PropertyText::from(""))
        );
        let reread = SdfRecord::from_sdf_with_params(&record.to_sdf().unwrap(), &params).unwrap();
        assert_eq!(reread.title(), Some(&PropertyText::from("")));
        #[cfg(feature = "cap-batch")]
        {
            let batch = SdfRecordStream::with_params(Cursor::new(sdf.as_bytes()), params)
                .batches(1, cosmolkit::BatchErrorMode::KeepErrors, Some(1))
                .unwrap()
                .next_batch()
                .unwrap()
                .unwrap();
            let cosmolkit::BatchRecord::Molecule(molecule) = batch.get(0).unwrap() else {
                panic!("valid carbon");
            };
            assert_eq!(molecule.properties().name(), Some(&PropertyText::from("")));
        }
    }
}
