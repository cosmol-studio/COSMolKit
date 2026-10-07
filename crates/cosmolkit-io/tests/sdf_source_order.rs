#![cfg(feature = "molecule")]
use cosmolkit_io::*;
use cosmolkit_model::{AtomId, PropertyValue};
use std::error::Error;
use std::io::Cursor;

fn mol(atoms: &[&str], bonds: &[(usize, usize)], properties: &str) -> String {
    let mut text = format!(
        "source-order\n  COSMolKit         2D\n\n{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000\n",
        atoms.len(),
        bonds.len()
    );
    for (index, symbol) in atoms.iter().enumerate() {
        text.push_str(&format!(
            "{:>10.4}{:>10.4}{:>10.4} {symbol:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n",
            index as f64, 0.0, 0.0
        ));
    }
    for &(first, second) in bonds {
        text.push_str(&format!("{first:>3}{second:>3}  1  0  0  0  0\n"));
    }
    text.push_str(properties);
    text.push_str("M  END\n");
    text
}
fn finalized() -> SdfDataReadParams {
    SdfDataReadParams {
        mol_post: Some(MolPostParams::default()),
        ..Default::default()
    }
}
fn root(error: &SdfReadError) -> &SdfReadError {
    match error {
        SdfReadError::Record { source, .. } => root(source),
        other => other,
    }
}

#[test]
fn mol_reader_finalizes_once_without_interpreting_trailing_sdf_fields() {
    let block = mol(&["C", "H"], &[(1, 2)], "");
    let text = format!("{block}>  <atom.iprop.rank>\n1 2 3\n\nmalformed\n$$$$\n");
    let parsed = read_mol_graph_record_detached_with_params(&text, finalized())
        .unwrap()
        .into_concrete()
        .unwrap();
    assert_eq!(parsed.topology.atoms.len(), 1);
    assert!(parsed.data_fields.is_empty());
    assert_eq!(parsed.topology.atoms[0].prop("rank"), None);
    let syntax = read_mol_graph_record_detached_with_params(&text, Default::default())
        .unwrap()
        .into_concrete()
        .unwrap();
    assert_eq!(syntax.topology.atoms.len(), 2);
    assert!(syntax.data_fields.is_empty());
    let bad = mol(
        &["C", "C", "C", "C", "C", "C"],
        &[(1, 2), (1, 3), (1, 4), (1, 5), (1, 6)],
        "",
    );
    assert!(matches!(
        read_mol_graph_record_detached_with_params(&bad, finalized()),
        Err(SdfReadError::MolPost(MolPostError::Processing(_)))
    ));
    assert!(read_mol_graph_record_detached_with_params(&bad, Default::default()).is_ok());
}

#[test]
fn lists_count_the_source_finalized_graph_and_syntax_only_policy_is_explicit() {
    // MolFromMolDataStream finishes removeHs before readMolProps expands lists.
    let block = mol(&["C", "H"], &[(1, 2)], "");
    let one = format!("{block}>  <atom.iprop.rank>\n7\n\n$$$$\n");
    let parsed = read_sdf_graph_record_detached_with_params(&one, finalized())
        .unwrap()
        .into_concrete()
        .unwrap();
    assert_eq!(parsed.topology.atoms.len(), 1);
    assert_eq!(
        parsed.topology.atoms[0].prop("rank"),
        Some(&PropertyValue::Int(7))
    );
    assert_eq!(parsed.data_fields, [("atom.iprop.rank".into(), "7".into())]);
    assert!(matches!(
        read_sdf_graph_record_detached_with_params(&one, Default::default()),
        Err(SdfReadError::PropertyListCount {
            actual: 1,
            expected: 2,
            ..
        })
    ));
    let two = format!("{block}>  <atom.iprop.rank>\n7 8\n\n$$$$\n");
    let params = SdfDataReadParams {
        mol_post: Some(MolPostParams {
            remove_hs: false,
            ..Default::default()
        }),
        ..Default::default()
    };
    let parsed = read_sdf_graph_record_detached_with_params(&two, params)
        .unwrap()
        .into_concrete()
        .unwrap();
    assert_eq!(parsed.topology.atoms.len(), 2);
    assert_eq!(
        parsed.topology.atoms[1].prop("rank"),
        Some(&PropertyValue::Int(8))
    );
}

#[test]
fn attachment_promotion_precedes_query_property_list_target_count() {
    let block = mol(&["C"], &[], "M  APO  1   1   1\n");
    let text = format!("{block}>  <atom.iprop.rank>\n7 9\n\n$$$$\n");
    let params = SdfDataReadParams {
        mol_post: Some(MolPostParams {
            sanitize: false,
            remove_hs: false,
            expand_attachment_points: true,
        }),
        ..Default::default()
    };
    let parsed = read_sdf_graph_record_detached_with_params(&text, params).unwrap();
    let MolBlockRecord::Query(query) = parsed.mol_block else {
        panic!("source attachment promotes to a query graph");
    };
    assert_eq!(query.query.num_atoms(), 2);
    assert_eq!(
        query.query.atoms()[0].prop("rank"),
        Some(&PropertyValue::Int(7))
    );
    assert_eq!(
        query.query.atoms()[1].prop("rank"),
        Some(&PropertyValue::Int(9))
    );
}

#[test]
fn property_list_failure_precedes_later_malformed_field_and_raw_mode_retains_source() {
    let block = mol(&["C"], &[], "");
    let text = format!("{block}>  <atom.iprop.rank>\n1 2\n\nmalformed\n$$$$\n");
    let error = read_sdf_graph_record_detached_with_params(&text, finalized()).unwrap_err();
    assert!(
        matches!(error, SdfReadError::PropertyListCount { target: "atom", ref name, actual: 2, expected: 1 } if name == "atom.iprop.rank")
    );
    let raw = SdfDataReadParams {
        process_property_lists: false,
        ..finalized()
    };
    assert!(
        matches!(read_sdf_graph_record_detached_with_params(&text, raw), Err(SdfReadError::Parse(message)) if message == "Problems encountered parsing data fields")
    );
    let permissive = SdfDataReadParams {
        strict_parsing: false,
        ..finalized()
    };
    let record = read_sdf_graph_record_detached_with_params(&text, permissive)
        .unwrap()
        .into_concrete()
        .unwrap();
    assert_eq!(
        record.properties.prop("atom.iprop.rank"),
        Some(&PropertyValue::from("1 2"))
    );
    assert_eq!(record.topology.atoms[0].prop("rank"), None);
}

#[test]
fn duplicate_fields_apply_in_encounter_order_and_keep_every_raw_field() {
    let block = mol(&["C", "C"], &[(1, 2)], "");
    let text = format!("{block}>  <atom.iprop.rank>\n1 2\n\n>  <atom.iprop.rank>\nn/a 9\n\n$$$$\n");
    let record = read_sdf_graph_record_detached_with_params(&text, finalized())
        .unwrap()
        .into_concrete()
        .unwrap();
    assert_eq!(
        record.topology.atoms[0].prop("rank"),
        Some(&PropertyValue::Int(1))
    );
    assert_eq!(
        record.topology.atoms[1].prop("rank"),
        Some(&PropertyValue::Int(9))
    );
    assert_eq!(
        record.properties.prop("atom.iprop.rank"),
        Some(&PropertyValue::from("n/a 9"))
    );
    assert_eq!(
        record.data_fields,
        [
            ("atom.iprop.rank".into(), "1 2".into()),
            ("atom.iprop.rank".into(), "n/a 9".into())
        ]
    );
}

#[test]
fn failed_field_recovery_starts_before_unread_later_field_values() {
    let block = mol(&["C"], &[], "");
    let text = format!("{block}>  <atom.iprop.rank>\n1 2\n\n>  <later>\n$$$$\n{block}$$$$\n");
    let first_boundary = text.find("$$$$\n").unwrap() + 5;
    let mut reader = SdfGraphReader::with_params(Cursor::new(text.as_bytes()), finalized());
    let error = reader.next_record().unwrap_err();
    assert!(matches!(
        root(&error),
        SdfReadError::PropertyListCount {
            actual: 2,
            expected: 1,
            ..
        }
    ));
    assert_eq!(reader.bytes_consumed(), first_boundary as u64);
    assert_eq!(reader.records_consumed(), 1);
    let next = reader
        .next_record()
        .unwrap()
        .unwrap()
        .into_concrete()
        .unwrap();
    assert_eq!(next.topology.atoms.len(), 1);
    assert!(next.data_fields.is_empty());
    assert!(reader.next_record().unwrap().is_none());
}

#[test]
fn finalization_error_precedes_data_syntax_and_retains_typed_cause_during_recovery() {
    let bad = mol(
        &["C", "C", "C", "C", "C", "C"],
        &[(1, 2), (1, 3), (1, 4), (1, 5), (1, 6)],
        "",
    );
    let good = mol(&["C"], &[], "");
    let text = format!("{bad}malformed\n$$$$\n{good}$$$$\n");
    let error = read_sdf_graph_record_detached_with_params(&text, finalized()).unwrap_err();
    assert!(matches!(
        &error,
        SdfReadError::MolPost(MolPostError::Processing(_))
    ));
    assert!(
        error
            .source()
            .unwrap()
            .downcast_ref::<MolPostError>()
            .is_some()
    );
    let mut reader = SdfGraphReader::with_params(Cursor::new(text.as_bytes()), finalized());
    let error = reader.next_record().unwrap_err();
    assert!(matches!(
        root(&error),
        SdfReadError::MolPost(MolPostError::Processing(_))
    ));
    let next = reader
        .next_record()
        .unwrap()
        .unwrap()
        .into_concrete()
        .unwrap();
    assert_eq!(next.topology.atoms[0].id(), AtomId::new(0));
    assert!(reader.next_record().unwrap().is_none());
}
