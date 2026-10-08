//! Counted-byte boundary controls from the pinned SMARTS preprocessing/scanner.

use cosmolkit_model::{PropertyText, RecursiveStructureQuery};
use cosmolkit_search::{SmartsParseError, SmartsParseParams, parse_smarts};

#[test]
fn counted_name_keeps_invalid_utf8_nul_and_vertical_tab_trim_order() {
    let params = SmartsParseParams::default();
    let graph = parse_smarts(b"CO name\xff\0\x80\x0b", &params).unwrap();
    assert_eq!(graph.num_atoms(), 2);
    assert_eq!(graph.name().unwrap().unwrap().as_bytes(), b"name\xff\0\x80");
    let no_cx = SmartsParseParams {
        allow_cxsmiles: false,
        ..params
    };
    let graph = parse_smarts(b"CO\tname\xff\0\x80\x0b", &no_cx).unwrap();
    assert_eq!(graph.name().unwrap().unwrap().as_bytes(), b"name\xff\0\x80");
}

#[test]
fn counted_cx_label_data_and_suffix_are_exact_source_bytes() {
    let graph = parse_smarts(b"C |$\xff$| tail\x80", &SmartsParseParams::default()).unwrap();
    assert_eq!(graph.num_atoms(), 1);
    let label = graph
        .atom(0)
        .unwrap()
        .props()
        .get(b"atomLabel".as_slice())
        .unwrap();
    assert_eq!(label.as_string().unwrap().as_bytes(), b"\xff");
    let prefix = graph.props().get(b"_CXSMILES_Data".as_slice()).unwrap();
    assert_eq!(prefix.as_string().unwrap().as_bytes(), b"|$\xff$|");
    assert_eq!(graph.name().unwrap().unwrap().as_bytes(), b"tail\x80");
}

#[test]
fn counted_internal_bad_byte_has_scanner_offset_not_unicode_rejection() {
    let error = parse_smarts(b"C[\xff]", &SmartsParseParams::default()).unwrap_err();
    assert!(matches!(error, SmartsParseError::UnexpectedCharacter {
        position: 3, character, ..
    } if character == 0xff));
}

#[test]
fn counted_embedded_nul_is_a_lexer_byte_not_a_c_string_terminator() {
    let error = parse_smarts(b"C\0C", &SmartsParseParams::default()).unwrap_err();
    assert!(matches!(
        error,
        SmartsParseError::UnexpectedCharacter {
            position: 2,
            character: b'\0',
            ..
        }
    ));
}

#[test]
fn counted_source_signed_char_window_trims_high_edge_bytes() {
    let graph = parse_smarts(b"\xffC\xfe", &SmartsParseParams::default()).unwrap();
    assert_eq!(graph.num_atoms(), 1);
    assert_eq!(graph.atom(0).unwrap().atomic_number(), 6);
}

#[test]
fn counted_chemical_error_precedes_valid_raw_cx_metadata() {
    let error = parse_smarts(b"C1C |$\xff$|", &SmartsParseParams::default()).unwrap_err();
    assert_eq!(error, SmartsParseError::Parse("unclosed ring".to_owned()));
}

#[test]
fn counted_replacements_preserve_raw_name_and_nonoverlapping_source_ranges() {
    let mut params = SmartsParseParams::default();
    params.replacements.insert("{X}".into(), "N".into());
    let graph = parse_smarts(b"{X}O name\xff", &params).unwrap();
    assert_eq!(graph.num_atoms(), 2);
    assert_eq!(graph.atom(0).unwrap().atomic_number(), 7);
    assert_eq!(graph.atom(1).unwrap().atomic_number(), 8);
    assert_eq!(graph.name().unwrap().unwrap().as_bytes(), b"name\xff");
}

#[test]
fn recursive_source_provenance_uses_one_canonical_byte_value_and_deep_clone() {
    let graph = parse_smarts("CC", &SmartsParseParams::default()).unwrap();
    let bytes = b"CC name\xff\0";
    let recursive = RecursiveStructureQuery::from_query_graph(graph, 17)
        .with_source_smarts(PropertyText::from_bytes(bytes));
    assert_eq!(recursive.source_smarts().unwrap().as_bytes(), bytes);
    let mut copied = recursive.clone();
    // QueryOps.h::copy uses quickCopy=true; ROMol.cpp initializes the empty
    // __computedProps entry even when the original has no molecule props.
    let mut expected = recursive.query_graph().unwrap().clone();
    expected
        .set_prop(
            "__computedProps",
            cosmolkit_model::PropertyValue::StringVector(vec![]),
        )
        .unwrap();
    assert_eq!(copied.query_graph(), Some(&expected));
    copied.set_query_graph(parse_smarts("N", &SmartsParseParams::default()).unwrap());
    assert_eq!(recursive.serial_number(), 17);
    assert_eq!(copied.serial_number(), 17);
    assert_eq!(copied.source_smarts().unwrap().as_bytes(), bytes);
    assert_eq!(copied.query_graph().unwrap().num_atoms(), 1);
    assert_eq!(recursive.query_graph().unwrap().num_atoms(), 2);
}

#[test]
fn all_trimmed_nonempty_input_copies_the_source_c_string_nul_into_scanner_data() {
    for input in [b" \t\r\n".as_slice(), b"\xff\x80".as_slice()] {
        let error = parse_smarts(input, &SmartsParseParams::default()).unwrap_err();
        assert!(matches!(
            error,
            SmartsParseError::UnexpectedCharacter {
                position: 1,
                character: b'\0',
                ..
            }
        ));
    }
    assert_eq!(
        parse_smarts(b"", &SmartsParseParams::default())
            .unwrap()
            .num_atoms(),
        0
    );
}
