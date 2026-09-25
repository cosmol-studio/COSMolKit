// Behavior expectations in this file are frozen from the vendored Gemmi
// implementation at 5cc1c23c6007e0e6cbd69289c6f7c0bff50e943e, principally
// cif.hpp::{rules,read_input,try_parse}, cifdoc.hpp::{char_table,quote},
// and sprintf.hpp::{to_str,to_str_prec,to_chars_z} through bundled stb.

use cosmolkit_io::cif::{
    CifCheckLevel, CifItem, CifReadErrorKind, check_cif_syntax, cif_as_char, cif_as_f64,
    cif_as_i32, cif_as_string, cif_is_null, format_cif_f32, format_cif_f64,
    format_cif_f64_precision, format_cif_i32, format_cif_usize, quote_cif_value, read_cif_document,
};

#[test]
fn gemmi_cif_control_words_are_case_insensitive_with_source_boundaries() {
    let input = "# leading comment\n\
        DATA_parent\n\
        _loop loop_thing\n\
        _stop stop_thing\n\
        _global global_thing\n\
        LOOP_\n\
        _loop.tag value\n\
        STOP_\n\
        SAVE_note\n\
        _inside yes\n\
        SAVE_\n\
        DATA_child\n\
        _value done\n";
    let document = read_cif_document(input, "controls.cif", CifCheckLevel::Syntax).unwrap();

    assert_eq!(document.source(), "controls.cif");
    assert_eq!(document.blocks().len(), 2);
    let parent = &document.blocks()[0];
    assert_eq!(parent.name(), "parent");
    assert_eq!(parent.find_value("_loop").unwrap().raw(), "loop_thing");
    assert_eq!(parent.find_value("_stop").unwrap().raw(), "stop_thing");
    assert_eq!(parent.find_value("_global").unwrap().raw(), "global_thing");
    assert!(matches!(parent.items()[3], CifItem::Loop(_)));
    let CifItem::Frame(frame) = &parent.items()[4] else {
        panic!("SAVE_note must open a save frame");
    };
    assert_eq!(frame.name(), "note");
    assert_eq!(frame.find_value("_inside").unwrap().decoded(), "yes");
    assert_eq!(document.blocks()[1].name(), "child");
    assert_eq!(
        document.blocks()[1].find_value("_value").unwrap().decoded(),
        "done"
    );
}

#[test]
fn gemmi_cif_comments_quotes_and_crlf_text_fields_preserve_raw_tokens() {
    let input = concat!(
        "# initial comment\r\n",
        "DaTa_text\r\n",
        "_single 'can't'\r\n",
        "_double \"a\"b\"\r\n",
        "_hash alpha#beta # trailing comment\r\n",
        "_text\r\n",
        ";first line\r\n",
        "second line\r\n",
        ";\r\n",
    );
    let document = read_cif_document(input, "tokens.cif", CifCheckLevel::Syntax).unwrap();
    let block = document.sole_block().unwrap();

    assert_eq!(block.name(), "text");
    let single = block.find_value("_single").unwrap();
    assert_eq!(single.raw(), "'can't'");
    assert_eq!(single.decoded(), "can't");
    let double = block.find_value("_double").unwrap();
    assert_eq!(double.raw(), "\"a\"b\"");
    assert_eq!(double.decoded(), "a\"b");
    let hash = block.find_value("_hash").unwrap();
    assert_eq!(hash.raw(), "alpha#beta");
    assert_eq!(
        block.find_value("_text").unwrap().raw(),
        ";first line\r\nsecond line\r\n;"
    );
    // Gemmi stores the item tag's line; it does not retain a separate value-token line.
    assert_eq!(block.find_pair("_text").unwrap().line(), 6);
}

#[test]
fn gemmi_cif_read_check_levels_keep_missing_value_scope() {
    let input = "data_a\n_x\n_next yes\n";
    assert!(read_cif_document(input, "missing.cif", CifCheckLevel::Syntax).is_ok());
    for level in [CifCheckLevel::Default, CifCheckLevel::Strict] {
        let error = read_cif_document(input, "missing.cif", level).unwrap_err();
        assert_eq!(error.kind(), CifReadErrorKind::MissingValue);
        assert_eq!(error.line(), 2);
        assert!(error.message().contains("_x has no value"));
    }
}

#[test]
fn gemmi_cif_read_check_levels_reject_case_folded_duplicate_tags() {
    let input = "data_a\n_x 1\n_X 2\n";
    assert!(read_cif_document(input, "duplicate.cif", CifCheckLevel::Syntax).is_ok());
    for level in [CifCheckLevel::Default, CifCheckLevel::Strict] {
        let error = read_cif_document(input, "duplicate.cif", level).unwrap_err();
        assert_eq!(error.kind(), CifReadErrorKind::DuplicateName);
        assert_eq!(error.line(), 3);
        assert!(error.message().contains("duplicate tag _X"));
    }
}

#[test]
fn gemmi_cif_strict_level_alone_rejects_bare_blocks_and_empty_loops() {
    let bare = "data_\n_x yes\n";
    for level in [CifCheckLevel::Syntax, CifCheckLevel::Default] {
        assert_eq!(
            read_cif_document(bare, "bare.cif", level).unwrap().blocks()[0].name(),
            " "
        );
    }
    assert_eq!(
        read_cif_document(bare, "bare.cif", CifCheckLevel::Strict)
            .unwrap_err()
            .kind(),
        CifReadErrorKind::InvalidValue
    );

    let empty_loop = "data_a\nloop_\n_x\n";
    for level in [CifCheckLevel::Syntax, CifCheckLevel::Default] {
        assert!(read_cif_document(empty_loop, "empty-loop.cif", level).is_ok());
    }
    let error = read_cif_document(empty_loop, "empty-loop.cif", CifCheckLevel::Strict).unwrap_err();
    assert_eq!(error.kind(), CifReadErrorKind::InvalidLoop);
    assert!(error.message().contains("empty loop with _x"));
}

#[test]
fn gemmi_cif_read_width_action_differs_from_syntax_check_action() {
    let input = "data_a\nloop_\n_x _y\nvalue\n";
    for level in [
        CifCheckLevel::Syntax,
        CifCheckLevel::Default,
        CifCheckLevel::Strict,
    ] {
        let error = read_cif_document(input, "width.cif", level).unwrap_err();
        // The pinned loop action throws `pegtl::parse_error`, which maps to
        // the syntax category rather than a document-validation InvalidLoop.
        assert_eq!(error.kind(), CifReadErrorKind::Syntax);
        assert_eq!(error.line(), 2);
        // Pinned PEGTL positions report a zero-based column at the loop action.
        assert_eq!(error.column(), 0);
        assert!(
            error
                .message()
                .contains("Wrong number of values in loop _*")
        );
    }
    assert!(check_cif_syntax(input, "width.cif").is_ok());
}

#[test]
fn gemmi_cif_check_action_reports_missing_pair_value_at_next_line() {
    let input = "data_a\n_x\n_next yes\n";
    let error = check_cif_syntax(input, "check-action.cif").unwrap_err();
    assert_eq!(error.kind(), CifReadErrorKind::Syntax);
    assert_eq!(error.source(), "check-action.cif");
    assert_eq!(error.line(), 3);
    assert_eq!(error.column(), 0);
    assert_eq!(error.message(), "tag without value");
}

#[test]
fn gemmi_cif_parse_error_reports_zero_based_location_after_crlf() {
    let input = "data_a\r\n_x yes\r\nbad\r\n";
    let error = read_cif_document(input, "bad-location.cif", CifCheckLevel::Syntax).unwrap_err();
    assert_eq!(error.kind(), CifReadErrorKind::Syntax);
    assert_eq!(error.source(), "bad-location.cif");
    assert_eq!(error.line(), 3);
    assert_eq!(error.column(), 0);
}

#[test]
fn gemmi_cif_unterminated_quote_reports_reference_eof_location() {
    let input = "data_a\n_x 'unterminated";
    let error = read_cif_document(input, "unterminated.cif", CifCheckLevel::Syntax).unwrap_err();
    assert_eq!(error.kind(), CifReadErrorKind::Syntax);
    assert_eq!(error.source(), "unterminated.cif");
    assert_eq!(error.line(), 2);
    assert_eq!(error.column(), 16);
    assert_eq!(error.message(), "unterminated 'string'");
}

#[test]
fn gemmi_cif_document_retains_block_item_and_row_major_value_order() {
    let input = concat!(
        "data_first\n",
        "_pair.before 'one value'\n",
        "loop_\n",
        "_order.second\n",
        "_order.first\n",
        "second0 first0 second1 first1\n",
        "save_note\n",
        "_frame.value inside\n",
        "save_\n",
        "_pair.after two\n",
        "data_second\n",
        "_next end\n",
    );
    let document = read_cif_document(input, "ordered.cif", CifCheckLevel::Syntax).unwrap();

    assert_eq!(
        document
            .blocks()
            .iter()
            .map(|block| block.name())
            .collect::<Vec<_>>(),
        ["first", "second"]
    );
    let first = &document.blocks()[0];
    assert_eq!(first.items().len(), 4);
    assert!(matches!(first.items()[0], CifItem::Pair(_)));
    assert!(matches!(first.items()[1], CifItem::Loop(_)));
    assert!(matches!(first.items()[2], CifItem::Frame(_)));
    assert!(matches!(first.items()[3], CifItem::Pair(_)));
    assert_eq!(
        first
            .find_pair("_pair.before")
            .unwrap()
            .value()
            .unwrap()
            .raw(),
        "'one value'"
    );
    assert_eq!(
        first
            .find_pair("_pair.after")
            .unwrap()
            .value()
            .unwrap()
            .raw(),
        "two"
    );

    let CifItem::Loop(loop_) = &first.items()[1] else {
        panic!("second item must remain the source loop");
    };
    assert_eq!(loop_.tags(), ["_order.second", "_order.first"]);
    assert_eq!(
        loop_
            .values()
            .iter()
            .map(|value| value.raw())
            .collect::<Vec<_>>(),
        ["second0", "first0", "second1", "first1"]
    );
    assert_eq!(loop_.value(1, 1).unwrap().raw(), "first1");

    let CifItem::Frame(frame) = &first.items()[2] else {
        panic!("third item must remain the source save frame");
    };
    assert_eq!(frame.name(), "note");
    assert_eq!(frame.find_value("_frame.value").unwrap().raw(), "inside");
}

#[test]
fn gemmi_cif_duplicate_scopes_cover_pair_loop_frame_and_block_names() {
    let duplicate_pair_loop = concat!(
        "data_first\n",
        "_CaseTag one\n",
        "loop_\n",
        "_casetag\n",
        "_other\n",
        "two three\n",
    );
    assert!(
        read_cif_document(
            duplicate_pair_loop,
            "duplicate-pair-loop.cif",
            CifCheckLevel::Syntax,
        )
        .is_ok()
    );
    let error = read_cif_document(
        duplicate_pair_loop,
        "duplicate-pair-loop.cif",
        CifCheckLevel::Default,
    )
    .unwrap_err();
    assert_eq!(error.kind(), CifReadErrorKind::DuplicateName);
    assert!(error.message().contains("duplicate tag _casetag"));

    let duplicate_frames = concat!(
        "data_first\n",
        "save_note\n",
        "_frame.value one\n",
        "save_\n",
        "save_NOTE\n",
        "_frame.value two\n",
        "save_\n",
    );
    assert_eq!(
        read_cif_document(
            duplicate_frames,
            "duplicate-frames.cif",
            CifCheckLevel::Default,
        )
        .unwrap_err()
        .kind(),
        CifReadErrorKind::DuplicateName
    );

    let duplicate_blocks = "data_First\n_x one\ndata_first\n_x two\n";
    assert_eq!(
        read_cif_document(
            duplicate_blocks,
            "duplicate-blocks.cif",
            CifCheckLevel::Default,
        )
        .unwrap_err()
        .kind(),
        CifReadErrorKind::DuplicateName
    );
    let repeated_globals = "global_\n_x one\nglobal_\n_x two\n";
    assert_eq!(
        read_cif_document(repeated_globals, "globals.cif", CifCheckLevel::Default,)
            .unwrap()
            .blocks()
            .len(),
        2
    );
}

#[test]
fn gemmi_cif_optional_tables_retain_loop_positions_and_pair_rows() {
    let input = concat!(
        "data_tables\n",
        "_pair.name 'pair value'\n",
        "_pair.maybe ?\n",
        "loop_\n",
        "_sample.second\n",
        "_sample.first\n",
        "_sample.present\n",
        "second0 first0 ? second1 first1 present1\n",
    );
    let block = read_cif_document(input, "tables.cif", CifCheckLevel::Syntax)
        .unwrap()
        .sole_block()
        .unwrap()
        .clone();

    let table = block
        .find("_sample.", &["first", "?absent", "second", "present"])
        .unwrap();
    assert!(table.is_present());
    assert_eq!((table.len(), table.width()), (2, 4));
    assert_eq!(table.tag(0), Some("_sample.first"));
    assert_eq!(table.tag(1), None);
    assert_eq!(table.tag(2), Some("_sample.second"));
    assert!(!table.has_column(1));
    let first = table.row(0).unwrap();
    assert_eq!(first.get(0).unwrap().raw(), "first0");
    assert!(!first.has(1));
    assert!(!first.has_value(1));
    assert!(first.has(3));
    assert!(!first.has_value(3));
    assert_eq!(first.one_of(1, 2).unwrap(), "second0");
    let second = table.row(1).unwrap();
    assert_eq!(second.decoded(0).as_deref(), Some("first1"));
    assert!(second.has_value(3));

    let pair_table = block.find("_pair.", &["name", "?absent"]).unwrap();
    assert!(pair_table.is_present());
    assert_eq!((pair_table.len(), pair_table.width()), (1, 2));
    let pair_row = pair_table.one().unwrap();
    assert_eq!(pair_row.decoded(0).as_deref(), Some("pair value"));
    assert!(pair_row.has(1) == false);
    assert_eq!(pair_row.get(1), None);
}

#[test]
fn gemmi_cif_row_one_of_preserves_source_presence_and_null_branches() {
    let input = concat!(
        "data_choices\n",
        "loop_\n",
        "_choice.anchor\n",
        "_choice.primary\n",
        "_choice.fallback\n",
        "a0 primary0 fallback0\n",
        "a1 ? fallback1\n",
        "a2 ? ?\n",
        "a3 ? .\n",
        "a4 . fallback4\n",
    );
    let document = read_cif_document(input, "one-of.cif", CifCheckLevel::Syntax).unwrap();
    let block = document.sole_block().unwrap();
    let table = block
        .find("_choice.", &["anchor", "primary", "fallback"])
        .unwrap();
    let selected = |row: isize| table.row(row).unwrap().one_of(1, 2).unwrap();
    assert_eq!(selected(0), "primary0");
    assert_eq!(selected(1), "fallback1");
    // Gemmi checks column presence for the fallback, not its nullness.
    assert_eq!(selected(2), "?");
    assert_eq!(selected(3), ".");
    assert_eq!(selected(4), "fallback4");

    let pair_input = concat!(
        "data_pair_choices\n",
        "_choice.anchor a\n",
        "_choice.primary\n",
        "_choice.fallback value\n",
        "_choice.null ?\n",
    );
    let pair_document =
        read_cif_document(pair_input, "one-of-pair.cif", CifCheckLevel::Syntax).unwrap();
    let pair_block = pair_document.sole_block().unwrap();
    let pair_table = pair_block
        .find("_choice.", &["anchor", "primary", "fallback", "null"])
        .unwrap();
    let pair = pair_table.one().unwrap();
    assert_eq!(pair.one_of(1, 2).unwrap(), "");
    assert_eq!(pair.one_of(3, 1).unwrap(), "");
    assert_eq!(pair.one_of(2, 9).unwrap(), "value");
    assert!(matches!(
        pair.one_of(9, 2),
        Err(error) if error.kind() == CifReadErrorKind::OutOfRange
    ));
    assert!(matches!(
        pair.one_of(3, 9),
        Err(error) if error.kind() == CifReadErrorKind::OutOfRange
    ));

    let absent_positions = block
        .find(
            "_choice.",
            &["anchor", "?missing_primary", "?missing_fallback"],
        )
        .unwrap();
    assert_eq!((absent_positions.width(), absent_positions.len()), (3, 5));
    for row_index in 0..5 {
        let row = absent_positions.row(row_index).unwrap();
        assert!(!row.has(1));
        assert!(!row.has(2));
        assert_eq!(
            row.one_of(1, 2).unwrap(),
            ".",
            "Gemmi Row::one_of returns its static raw '.' when both optional positions are absent"
        );
    }
}

#[test]
fn gemmi_cif_category_and_table_lookup_preserve_order_and_errors() {
    let input = concat!(
        "data_categories\n",
        "_PAIR.first pair0\n",
        "_PAIR.second pair1\n",
        "loop_\n",
        "_LOOP.second\n",
        "_LOOP.first\n",
        "second0 first0 second1 first1\n",
        "save_note\n",
        "_FRAME.value ignored\n",
        "save_\n",
        "_MISC.value m0\n",
    );
    let document = read_cif_document(input, "categories.cif", CifCheckLevel::Syntax).unwrap();
    let block = document.sole_block().unwrap();

    assert_eq!(
        block
            .mmcif_category_names()
            .iter()
            .map(String::as_str)
            .collect::<Vec<_>>(),
        ["_PAIR.", "_LOOP.", "_MISC."]
    );

    // Gemmi appends a missing dot and compares category/tag names
    // case-insensitively, while retaining loop column and row order.
    let loop_table = block.find_mmcif_category("_loop").unwrap();
    assert_eq!((loop_table.width(), loop_table.len()), (2, 2));
    assert_eq!(loop_table.tag(0), Some("_LOOP.second"));
    assert_eq!(loop_table.tag(1), Some("_LOOP.first"));
    assert_eq!(loop_table.row(0).unwrap().get(0).unwrap().raw(), "second0");
    assert_eq!(loop_table.row(0).unwrap().get(1).unwrap().raw(), "first0");
    assert_eq!(loop_table.row(1).unwrap().get(0).unwrap().raw(), "second1");
    assert_eq!(loop_table.row(1).unwrap().get(1).unwrap().raw(), "first1");

    let pair_table = block.find_mmcif_category("_pair.").unwrap();
    assert_eq!((pair_table.width(), pair_table.len()), (2, 1));
    let pair = pair_table.one().unwrap();
    assert_eq!(pair.get(0).unwrap().raw(), "pair0");
    assert_eq!(pair.get(1).unwrap().raw(), "pair1");

    let absent = block.find_mmcif_category("_missing").unwrap();
    assert!(!absent.is_present());
    assert_eq!((absent.width(), absent.len()), (0, 0));
    assert_eq!(block.find("_LOOP.", &["absent"]).unwrap().len(), 0);
    assert_eq!(
        block.find("_LOOP.", &["?first"]).unwrap_err().kind(),
        CifReadErrorKind::InvalidValue
    );
    assert_eq!(
        block
            .find_mmcif_category("category-without-underscore")
            .unwrap_err()
            .kind(),
        CifReadErrorKind::InvalidValue
    );

    let mixed = "data_mixed\nloop_\n_mixed.one\n_other.two\na b\n";
    let mixed_block = read_cif_document(mixed, "mixed-category.cif", CifCheckLevel::Syntax)
        .unwrap()
        .sole_block()
        .unwrap()
        .clone();
    let error = mixed_block.find_mmcif_category("_mixed").unwrap_err();
    assert_eq!(error.kind(), CifReadErrorKind::InvalidLoop);
    assert_eq!(error.message(), "Tag _other.two in loop with _mixed.");
}

#[test]
fn gemmi_cif_null_string_integer_and_number_conversions_follow_source() {
    assert!(cif_is_null("."));
    assert!(cif_is_null("?"));
    assert!(!cif_is_null("??"));
    assert_eq!(cif_as_string("'one two'"), "one two");
    assert_eq!(cif_as_string(";body\n;"), "body");
    assert_eq!(cif_as_char("?", '?').unwrap(), '?');
    assert_eq!(cif_as_char("", '?').unwrap(), '\0');
    assert_eq!(cif_as_char("A", '?').unwrap(), 'A');
    assert_eq!(cif_as_char("''", '?').unwrap(), '\0');
    assert_eq!(cif_as_char("'A'", '?').unwrap(), 'A');
    assert!(cif_as_char("é", '?').is_err());

    assert_eq!(cif_as_i32(" \t\n\u{000b}\u{000c}+42\r ").unwrap(), 42);
    assert_eq!(cif_as_i32("-42").unwrap(), -42);
    assert_eq!(cif_as_i32("-2147483648").unwrap(), i32::MIN);
    assert_eq!(cif_as_i32("2147483647").unwrap(), i32::MAX);
    assert!(cif_as_i32("42tail").is_err());
    assert!(cif_as_i32("+").is_err());

    let null = -77.0;
    assert_eq!(cif_as_f64("+2.5", null).unwrap(), 2.5);
    assert_eq!(cif_as_f64("1.25(3)", null).unwrap(), 1.25);
    assert_eq!(cif_as_f64("1.25()", null).unwrap(), 1.25);
    let negative_zero = cif_as_f64("-0", null).unwrap();
    assert_eq!(negative_zero, 0.0);
    assert!(negative_zero.is_sign_negative());
    for invalid in ["", "1.25(3)x", "1e999", "inf", "nan", "not-a-number", "?"] {
        assert_eq!(
            cif_as_f64(invalid, null).unwrap(),
            null,
            "Gemmi returns the caller's null result for {invalid:?}"
        );
    }
}

#[test]
fn gemmi_cif_numeric_scanner_matches_source_sign_whitespace_and_uncertainty() {
    let null = -77.25_f64;
    let null_bits = 0xc053_5000_0000_0000;
    assert_eq!(null.to_bits(), null_bits);

    for (input, expected_bits) in [
        ("++1", null_bits),
        ("+1", 0x3ff0_0000_0000_0000),
        ("-1", 0xbff0_0000_0000_0000),
        ("+-1", 0xbff0_0000_0000_0000),
        ("-+1", null_bits),
        (" +1", null_bits),
        ("1 ", null_bits),
        (".5", 0x3fe0_0000_0000_0000),
        ("+.5", 0x3fe0_0000_0000_0000),
        ("-.5", 0xbfe0_0000_0000_0000),
        ("1.", 0x3ff0_0000_0000_0000),
        ("+.5(2)", 0x3fe0_0000_0000_0000),
        ("1.25()", 0x3ff4_0000_0000_0000),
        ("1.25(3)", 0x3ff4_0000_0000_0000),
        ("1.25(000)", 0x3ff4_0000_0000_0000),
        ("1.25(3)x", null_bits),
        ("1.25(3))", null_bits),
        ("1.25(+3)", null_bits),
    ] {
        let actual = cif_as_f64(input, null).unwrap();
        assert_eq!(actual.to_bits(), expected_bits, "{input:?}");
    }
}

#[test]
fn gemmi_cif_underflow_range_errors_fall_back_but_exact_zero_and_subnormals_survive() {
    // Oracle: pinned Gemmi 0.7.5 / fast_float at
    // 5cc1c23c6007e0e6cbd69289c6f7c0bff50e943e, GNU C++ 15.2.0,
    // Ubuntu 15.2.0-16ubuntu1, -std=c++17, default general format.
    // Expectations are from the recorded fixed-source probe, not Rust parsing.
    let null = -77.25_f64;
    let null_bits = 0xc053_5000_0000_0000;
    assert_eq!(null.to_bits(), null_bits);

    for (input, expected_bits) in [
        ("1e-999", null_bits),
        ("-1e-999", null_bits),
        ("+1e-999", null_bits),
        ("1e-324", null_bits),
        ("-1e-324", null_bits),
        ("2e-324", null_bits),
        ("2.4703282292062326e-324", null_bits),
        ("2.4703282292062327e-324", null_bits),
        ("-2.4703282292062327e-324", null_bits),
        ("2.4703282292062328e-324", 0x0000_0000_0000_0001),
        ("2.4703282292062329e-324", 0x0000_0000_0000_0001),
        ("-2.4703282292062328e-324", 0x8000_0000_0000_0001),
        ("3e-324", 0x0000_0000_0000_0001),
        ("-3e-324", 0x8000_0000_0000_0001),
        ("1e-323", 0x0000_0000_0000_0002),
        ("-1e-323", 0x8000_0000_0000_0002),
        ("0e-999", 0x0000_0000_0000_0000),
        ("-0e-999", 0x8000_0000_0000_0000),
        ("0.000e+999", 0x0000_0000_0000_0000),
        ("-0.000e+999", 0x8000_0000_0000_0000),
        ("1e-323(4)", 0x0000_0000_0000_0002),
        ("1e-999(4)", null_bits),
        ("0e-999(4)", 0x0000_0000_0000_0000),
        ("-0e-999(4)", 0x8000_0000_0000_0000),
    ] {
        let actual = cif_as_f64(input, null).unwrap();
        assert_eq!(actual.to_bits(), expected_bits, "{input:?}");
    }
}

#[test]
fn gemmi_cif_quote_matches_source_branch_priority() {
    for (input, expected) in [
        ("plain", "plain"),
        ("", "''"),
        (".", "'.'"),
        ("?", "'?'"),
        ("with space", "'with space'"),
        ("can't", "\"can't\""),
        ("a\"b", "'a\"b'"),
        ("both '\"", ";both '\"\n;"),
        ("line1\nline2", ";line1\nline2\n;"),
        ("é", "'é'"),
        ("_tag", "'_tag'"),
        ("loop_", "'loop_'"),
        ("#comment", "'#comment'"),
    ] {
        assert_eq!(quote_cif_value(input.to_owned()), expected, "{input:?}");
    }
}

#[test]
fn gemmi_cif_typed_numeric_formatting_matches_bundled_source_profile() {
    for (value, expected) in [
        (0.0, "0"),
        (-0.0, "-0"),
        (1.234_567_895, "1.2345679"),
        (9.999_999_995, "9.99999999"),
        (1e-4, "0.0001"),
        (1e-5, "1e-05"),
        (9.999_999_995e-5, "0.0001"),
        (99_999_999.0, "99999999"),
        (1e8, "100000000"),
        (100_000_001.0, "100000001"),
        (-1e8, "-100000000"),
        (f64::MAX, "1.79769313e+308"),
        (f64::from_bits(1), "4.94065646e-324"),
        (f64::INFINITY, "Inf"),
        (f64::NEG_INFINITY, "-Inf"),
        (f64::NAN, "NaN"),
    ] {
        assert_eq!(format_cif_f64(value), expected, "f64 {value:?}");
    }

    for (value, expected) in [
        (0.0_f32, "0"),
        (-0.0_f32, "-0"),
        (1.234_567_89_f32, "1.23457"),
        (9.999_995_f32, "10"),
        (f32::MAX, "3.40282e+38"),
    ] {
        assert_eq!(format_cif_f32(value), expected, "f32 {value:?}");
    }

    assert_eq!(format_cif_f64_precision::<0>(2.5), "3");
    assert_eq!(format_cif_f64_precision::<0>(-2.5), "-3");
    assert_eq!(format_cif_f64_precision::<0>(-0.0), "-0");
    assert_eq!(format_cif_f64_precision::<3>(1.23456), "1.235");
    assert_eq!(format_cif_f64_precision::<3>(-0.0), "-0.000");
    assert_eq!(format_cif_f64_precision::<3>(9.9995), "9.999");
    assert_eq!(format_cif_f64_precision::<3>(1e-4), "0.000");
    assert_eq!(format_cif_f64_precision::<3>(1e-5), "0.000");
    assert_eq!(format_cif_f64_precision::<3>(99_999_999.0), "99999999.000");
    assert_eq!(format_cif_f64_precision::<3>(1e8), "1e+08");
    assert_eq!(format_cif_f64_precision::<3>(-1e8), "-1e+08");
    assert_eq!(format_cif_f64_precision::<6>(1e-5), "0.000010");
    assert_eq!(format_cif_i32(i32::MIN), "-2147483648");
    assert_eq!(format_cif_i32(i32::MAX), "2147483647");
    #[cfg(target_pointer_width = "64")]
    assert_eq!(format_cif_usize(usize::MAX), "18446744073709551615");
    #[cfg(target_pointer_width = "32")]
    assert_eq!(format_cif_usize(usize::MAX), "4294967295");
}

#[test]
fn gemmi_cif_general_format_matches_bundled_stb_rounding_and_signs() {
    let f32_tie = 123_456.5_f32;
    assert_eq!(
        format_cif_f32(f32::from_bits(f32_tie.to_bits() - 1)),
        "123456"
    );
    assert_eq!(format_cif_f32(f32_tie), "123457");
    assert_eq!(
        format_cif_f32(f32::from_bits(f32_tie.to_bits() + 1)),
        "123457"
    );

    let f64_tie = 123_456_788.5_f64;
    assert_eq!(
        format_cif_f64(f64::from_bits(f64_tie.to_bits() - 1)),
        "123456788"
    );
    assert_eq!(format_cif_f64(f64_tie), "123456789");
    assert_eq!(
        format_cif_f64(f64::from_bits(f64_tie.to_bits() + 1)),
        "123456789"
    );

    let negative_tie = -f64_tie;
    assert_eq!(
        format_cif_f64(f64::from_bits(negative_tie.to_bits() + 1)),
        "-123456789"
    );
    assert_eq!(format_cif_f64(negative_tie), "-123456789");
    assert_eq!(
        format_cif_f64(f64::from_bits(negative_tie.to_bits() - 1)),
        "-123456788"
    );

    let f32_carry = 999_999.5_f32;
    assert_eq!(
        format_cif_f32(f32::from_bits(f32_carry.to_bits() - 1)),
        "999999"
    );
    assert_eq!(format_cif_f32(f32_carry), "1e+06");
    assert_eq!(
        format_cif_f32(f32::from_bits(f32_carry.to_bits() + 1)),
        "1e+06"
    );

    let f64_carry = 999_999_999.5_f64;
    assert_eq!(
        format_cif_f64(f64::from_bits(f64_carry.to_bits() - 1)),
        "999999999"
    );
    assert_eq!(format_cif_f64(f64_carry), "1e+09");
    assert_eq!(
        format_cif_f64(f64::from_bits(f64_carry.to_bits() + 1)),
        "1e+09"
    );
    assert_eq!(
        format_cif_f64(f64::from_bits((-f64_carry).to_bits() + 1)),
        "-1e+09"
    );
    assert_eq!(format_cif_f64(-f64_carry), "-1e+09");
    assert_eq!(
        format_cif_f64(f64::from_bits((-f64_carry).to_bits() - 1)),
        "-999999999"
    );

    assert_eq!(format_cif_f64(-0.0), "-0");
    assert_eq!(format_cif_f32(-0.0), "-0");
    assert_eq!(
        format_cif_f64(f64::from_bits(0xfff8_0000_0000_0000)),
        "-NaN"
    );
    assert_eq!(format_cif_f32(f32::from_bits(0xffc0_0000)), "-NaN");
}
