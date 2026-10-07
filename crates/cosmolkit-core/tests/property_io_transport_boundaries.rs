//! Proposed fixed boundaries for the existing CORE numeric transport owner.
//! Source anchors reside in the implementing functions; independent acceptance pending.
use cosmolkit_core::{
    DoubleLexicalReadError, DoubleLexicalReadErrorKind, PropertyDoubleReadError,
    PropertyUIntReadError, SourceNumericStreamState, property_value_to_double,
    property_value_to_string, property_value_to_uint, source_field_double, source_lexical_double,
    source_unsigned_stream_array, source_unsigned_stream_read,
};
use cosmolkit_model::{PropertyText, PropertyValue};

#[test]
fn counted_string_getter_keeps_nul_and_non_utf8_bytes() {
    let bytes = b"name\0\xff<sub>x</sub>";
    let value = PropertyValue::String(PropertyText::from(bytes.as_slice()));
    assert_eq!(property_value_to_string(&value).unwrap().as_bytes(), bytes);
    assert_eq!(
        value,
        PropertyValue::String(PropertyText::from(bytes.as_slice()))
    );
}

#[test]
fn source_double_tags_keep_bits_and_complete_failed_property_bytes() {
    for value in [-0.0, 1.25, f64::NEG_INFINITY, f64::INFINITY] {
        assert_eq!(
            property_value_to_double(&PropertyValue::Double(value))
                .unwrap()
                .to_bits(),
            value.to_bits()
        );
    }
    for value in [
        PropertyValue::Int(1),
        PropertyValue::UInt(1),
        PropertyValue::Bool(true),
        PropertyValue::IntVector(vec![1]),
        PropertyValue::StringVector(vec![]),
    ] {
        assert_eq!(
            property_value_to_double(&value),
            Err(PropertyDoubleReadError::InvalidKind { kind: value.kind() })
        );
    }
    let original = PropertyText::from(b" 1.25\0\xff \t".as_slice());
    match property_value_to_double(&PropertyValue::String(original.clone())).unwrap_err() {
        PropertyDoubleReadError::Lexical { value, source } => {
            assert_eq!(value, original);
            assert_eq!(
                source,
                DoubleLexicalReadError {
                    kind: DoubleLexicalReadErrorKind::Syntax,
                    consumed: 0,
                    state: SourceNumericStreamState {
                        eof: false,
                        fail: true,
                        bad: false
                    },
                    value_bits: Some(0),
                }
            );
        }
        other => panic!("wrong source cause: {other:?}"),
    }
    assert_eq!(
        property_value_to_double(&PropertyValue::String("1.25 \t\r\n".into())).unwrap(),
        1.25
    );
}

#[test]
fn counted_special_double_payloads_and_finite_grammar_remain_distinct() {
    assert!(source_lexical_double(b"nan(\0\xff)").unwrap().is_nan());
    assert!(
        source_lexical_double(b"-NaN(\xff)")
            .unwrap()
            .is_sign_negative()
    );
    assert_eq!(source_lexical_double(b"+INFINITY").unwrap(), f64::INFINITY);
    assert_eq!(
        source_lexical_double(b"-1e-9999").unwrap().to_bits(),
        (-0.0_f64).to_bits()
    );
    assert_eq!(
        source_lexical_double(b"1e309"),
        Err(DoubleLexicalReadError {
            kind: DoubleLexicalReadErrorKind::Overflow,
            consumed: 5,
            state: SourceNumericStreamState {
                eof: true,
                fail: true,
                bad: false
            },
            value_bits: Some(f64::MAX.to_bits()),
        })
    );
    for (input, kind, consumed, eof, fail, value_bits) in [
        (
            b" 1".as_slice(),
            DoubleLexicalReadErrorKind::Syntax,
            0,
            false,
            true,
            0,
        ),
        (
            b"1 ".as_slice(),
            DoubleLexicalReadErrorKind::TrailingByte,
            2,
            false,
            false,
            1.0_f64.to_bits(),
        ),
        (
            b"1e+".as_slice(),
            DoubleLexicalReadErrorKind::Syntax,
            3,
            true,
            true,
            0,
        ),
        (
            b"0x1p0".as_slice(),
            DoubleLexicalReadErrorKind::TrailingByte,
            2,
            false,
            false,
            0,
        ),
        (
            b"1\0".as_slice(),
            DoubleLexicalReadErrorKind::TrailingByte,
            2,
            false,
            false,
            1.0_f64.to_bits(),
        ),
        (
            b"1\xff".as_slice(),
            DoubleLexicalReadErrorKind::TrailingByte,
            2,
            false,
            false,
            1.0_f64.to_bits(),
        ),
    ] {
        assert_eq!(
            source_lexical_double(input),
            Err(DoubleLexicalReadError {
                kind,
                consumed,
                state: SourceNumericStreamState {
                    eof,
                    fail,
                    bad: false
                },
                value_bits: Some(value_bits),
            }),
            "{input:?}",
        );
    }
}

#[test]
fn field_double_uses_its_source_space_and_c_string_boundary() {
    assert_eq!(source_field_double(b"  -1.25  \0\xffjunk").unwrap(), -1.25);
    assert!(source_field_double(b"\t1.25").is_err());
    assert!(source_field_double(b"1.25\t").is_err());
    assert_eq!(
        source_field_double(b"     "),
        Err(DoubleLexicalReadError {
            kind: DoubleLexicalReadErrorKind::Syntax,
            consumed: 0,
            state: SourceNumericStreamState {
                eof: true,
                fail: true,
                bad: false
            },
            value_bits: Some(0),
        })
    );
    assert!(
        property_value_to_double(&PropertyValue::String(PropertyText::from(
            b"1.25\0junk".as_slice()
        )))
        .is_err()
    );
}

#[test]
fn unsigned_property_getter_uses_the_existing_full_width_owner() {
    assert_eq!(
        property_value_to_uint(&PropertyValue::UInt(u32::MAX)).unwrap(),
        u32::MAX
    );
    assert_eq!(
        property_value_to_uint(&PropertyValue::String("4294967295 \t".into())).unwrap(),
        u32::MAX
    );
    assert_eq!(
        property_value_to_uint(&PropertyValue::String("-1".into())).unwrap(),
        u32::MAX
    );
    assert!(property_value_to_uint(&PropertyValue::Int(-1)).is_err());
    let wrong = PropertyValue::StringVector(vec![]);
    assert_eq!(
        property_value_to_uint(&wrong),
        Err(PropertyUIntReadError::InvalidKind { kind: wrong.kind() })
    );
}

#[test]
fn unsigned_stream_retains_conversion_assignments_and_sticky_failure() {
    let mut position = 0;
    let mut failed = false;
    let input = b" -1 0\xff 2";
    assert_eq!(
        source_unsigned_stream_read(input, &mut position, &mut failed),
        Some(u32::MAX)
    );
    assert_eq!(
        source_unsigned_stream_read(input, &mut position, &mut failed),
        Some(0)
    );
    assert!(!failed);
    assert_eq!(
        source_unsigned_stream_read(input, &mut position, &mut failed),
        Some(0)
    );
    assert!(failed);
    let stopped = position;
    assert_eq!(
        source_unsigned_stream_read(input, &mut position, &mut failed),
        None
    );
    assert_eq!(position, stopped);
    for (input, assignment) in [(b"-x".as_slice(), 0), (b"4294967296", u32::MAX)] {
        position = 0;
        failed = false;
        assert_eq!(
            source_unsigned_stream_read(input, &mut position, &mut failed),
            Some(assignment)
        );
        assert!(failed);
        assert_eq!(
            source_unsigned_stream_read(input, &mut position, &mut failed),
            None
        );
    }
}

#[test]
fn unsigned_counted_array_retains_order_duplicates_and_both_warning_facts() {
    assert_eq!(
        source_unsigned_stream_array(b"(3 7 11 7)").unwrap(),
        (vec![7, 11, 7], true, true)
    );
    assert_eq!(
        source_unsigned_stream_array(b"X2 9 4)").unwrap(),
        (vec![9, 4], false, true)
    );
    assert_eq!(
        source_unsigned_stream_array(b"(2 9 4X").unwrap(),
        (vec![9, 4], true, false)
    );
    assert_eq!(
        source_unsigned_stream_array(b"(3 7)").unwrap(),
        (vec![7, 0, 0], true, false)
    );
    assert_eq!(
        source_unsigned_stream_array(b"(0)").unwrap(),
        (vec![], true, true)
    );
}
