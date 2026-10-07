#![cfg(feature = "cap-io")]
//! Proposed non-strict boundary coverage. Original strict source differences remain in receipts.
use cosmolkit::*;

const INPUT: &str =
    include_str!("../../../testdata/molblock/fixtures/io44-sdf-row34-property-list-policy.sdf");

fn assert_non_strict_properties(molecule: &Molecule) {
    assert_eq!(molecule.num_atoms(), 2);
    assert_eq!(molecule.num_bonds(), 1);
    assert_eq!(
        molecule
            .atom(AtomId::new(0))
            .unwrap()
            .prop("note")
            .map(PropertyValue::as_string)
            .transpose()
            .unwrap(),
        Some(&PropertyText::from("alpha"))
    );
    assert_eq!(
        molecule
            .atom(AtomId::new(1))
            .unwrap()
            .prop("note")
            .map(PropertyValue::as_string)
            .transpose()
            .unwrap(),
        Some(&PropertyText::from("beta"))
    );
    for atom in [AtomId::new(0), AtomId::new(1)] {
        assert_eq!(molecule.atom(atom).unwrap().prop("rank"), None);
    }
    assert_eq!(molecule.bond(BondId::new(0)).unwrap().prop("flag"), None);
    let fields = molecule.properties().sdf_data_fields();
    assert!(
        fields
            .iter()
            .any(|(name, value)| name.as_bytes() == b"atom.iprop.rank"
                && value.as_bytes() == b"10 20 30 40 50 60 70 80 90 100")
    );
    assert!(
        fields
            .iter()
            .any(|(name, value)| name.as_bytes() == b"bond.prop.flag"
                && value.as_bytes() == b"one two three four five six seven eight nine")
    );
}

#[test]
fn non_strict_row34_preserves_raw_lists_and_expands_only_matching_counts() {
    let params = SdfReadParams {
        strict_parsing: false,
        ..Default::default()
    };
    let record = SdfRecord::from_sdf_with_params(INPUT, &params).unwrap();
    assert_non_strict_properties(record.molecule().unwrap());
    let molecule = Molecule::from_sdf_with_params(INPUT, &params).unwrap();
    assert_non_strict_properties(&molecule);
    assert_eq!(
        record.data_fields(),
        molecule.properties().sdf_data_fields()
    );
    #[cfg(feature = "cap-batch")]
    {
        let batch = SdfRecordStream::with_params(std::io::Cursor::new(INPUT.as_bytes()), params)
            .batches(1, BatchErrorMode::KeepErrors, Some(1))
            .unwrap()
            .next_batch()
            .unwrap()
            .unwrap();
        let BatchRecord::Molecule(molecule) = batch.get(0).unwrap() else {
            panic!("source valid non-strict record");
        };
        assert_non_strict_properties(molecule);
        assert_eq!(
            record.data_fields(),
            molecule.properties().sdf_data_fields()
        );
    }
}

fn strict_property_cause_context(
    error: &(dyn std::error::Error + 'static),
) -> Vec<(usize, u64, usize)> {
    let mut current = Some(error);
    let mut records = Vec::new();
    while let Some(cause) = current {
        if let Some(cause) = cause.downcast_ref::<SdfReadError>() {
            match cause {
                SdfReadError::Record {
                    index,
                    byte_offset,
                    line_offset,
                    ..
                } => {
                    records.push((*index, *byte_offset, *line_offset));
                }
                SdfReadError::PropertyListCount {
                    target,
                    name,
                    actual,
                    expected,
                } => {
                    assert_eq!(
                        (*target, name.as_str(), *actual, *expected),
                        ("atom", "atom.iprop.rank", 10, 2)
                    );
                    return records;
                }
                other => panic!("unexpected structured reader cause: {other:?}"),
            }
        }
        current = cause.source();
    }
    panic!("missing concrete property-list count cause: {error}");
}

#[test]
fn strict_row34_preserves_concrete_cause_and_record_context_across_readers() {
    use std::io::{Cursor, Write};
    let params = SdfReadParams::default();
    let error = SdfRecord::from_sdf_with_params(INPUT, &params).unwrap_err();
    assert!(strict_property_cause_context(&error).is_empty());
    let error = Molecule::from_sdf_with_params(INPUT, &params).unwrap_err();
    assert!(strict_property_cause_context(&error).is_empty());
    let input = INPUT.repeat(2);
    let context = |index: usize| {
        (
            index,
            (index * INPUT.len()) as u64,
            index * INPUT.lines().count(),
        )
    };
    let mut stream = SdfRecordStream::with_params(Cursor::new(input.as_bytes()), params);
    for index in 0..2 {
        let error = stream.next_record().unwrap_err();
        assert_eq!(strict_property_cause_context(&error), vec![context(index)]);
        assert_eq!(stream.records_consumed(), index + 1);
    }
    assert!(stream.next_record().unwrap().is_none());
    assert_eq!(stream.bytes_consumed(), input.len() as u64);
    assert_eq!(stream.lines_consumed(), input.lines().count());
    let path = std::env::temp_dir().join(format!(
        "cosmolkit-sdf-row34-cause-{}.sdf",
        std::process::id()
    ));
    let mut file = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&path)
        .unwrap();
    file.write_all(input.as_bytes()).unwrap();
    drop(file);
    let dataset = SdfDataset::open_with_params(path.to_str().unwrap(), &params).unwrap();
    assert_eq!(dataset.len(), 2);
    for index in 0..2 {
        let metadata = dataset.metadata(index).unwrap();
        assert_eq!(
            (metadata.index, metadata.byte_offset, metadata.line_offset),
            context(index)
        );
        assert_eq!(dataset.record_text(index).unwrap(), INPUT);
        let error = dataset.record(index).unwrap_err();
        assert_eq!(strict_property_cause_context(&error), vec![context(index)]);
    }
    #[cfg(feature = "cap-batch")]
    {
        let batch = SdfRecordStream::with_params(Cursor::new(input.as_bytes()), params)
            .batches(2, BatchErrorMode::KeepErrors, Some(1))
            .unwrap()
            .next_batch()
            .unwrap()
            .unwrap();
        assert_eq!(batch.len(), 2);
        for index in 0..2 {
            let BatchRecord::Error(error) = batch.get(index).unwrap() else {
                panic!("strict count mismatch must remain a failed record");
            };
            assert_eq!((error.index, error.operation), (index, "read_sdf"));
            assert_eq!(strict_property_cause_context(error), vec![context(index)]);
        }
    }
    drop(dataset);
    std::fs::remove_file(path).unwrap();
}

#[test]
fn nested_sdf_record_sources_preserve_each_concrete_context_and_display() {
    use std::error::Error;
    let mut error = SdfReadError::PropertyListCount {
        target: "atom",
        name: "atom.iprop.rank".into(),
        actual: 10,
        expected: 2,
    };
    for index in 0..3 {
        error = SdfReadError::Record {
            index,
            byte_offset: index as u64 * 100,
            line_offset: index * 10,
            source: Box::new(error),
        };
    }
    let before = error.clone();
    assert_eq!(
        strict_property_cause_context(&error),
        vec![(2, 200, 20), (1, 100, 10), (0, 0, 0)]
    );
    assert!(error.source().unwrap().is::<SdfReadError>());
    assert_eq!(
        error.to_string(),
        "SDF record 2 at byte 200, line 20 failed: SDF record 1 at byte 100, line 10 failed: SDF record 0 at byte 0, line 0 failed: SDF atom property list 'atom.iprop.rank' has 10 values, expected 2"
    );
    assert_eq!(error, before);
}
