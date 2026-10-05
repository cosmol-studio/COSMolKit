//! Retained catalog regressions migrated through the detached public boundary.
//! Original assertions remain in immutable TAU-source-audit/ready. Corrected
//! constructor/default/511 boundaries cite full source and ROOT8426 decision.
use cosmolkit_model::QueryGraph;
use cosmolkit_search::{CompiledQuery, SmartsParseParams, parse_smarts};
use cosmolkit_tautomer::{
    TautomerCatalog, TautomerCatalogError, TautomerTransform, TautomerTransformError,
};
use cosmolkit_types::BondOrder;
use std::io::{self, BufRead, Cursor, Read};
use std::path::Path;
// Source-derived fixed73-row data; no runtime source parsing or oracle call.
mod source_rows {
    include!("../src/transforms.rs");
}
use source_rows::{
    CURRENT_TAUTOMER_TRANSFORM_DEFINITIONS, TautomerTransformDefinition,
    V1_TAUTOMER_TRANSFORM_DEFINITIONS,
};
fn transform_error<T>(
    result: Result<T, TautomerCatalogError>,
) -> Result<T, TautomerTransformError> {
    match result {
        Ok(value) => Ok(value),
        Err(TautomerCatalogError::Transform(error)) => Err(error),
        Err(error) => panic!("unexpected catalog boundary: {error}"),
    }
}
fn transform_from_fields(
    name: &str,
    smarts: &str,
    bonds: &str,
    charges: &str,
) -> Result<TautomerTransform, TautomerTransformError> {
    transform_error(
        TautomerCatalog::from_data(&[(name, smarts, bonds, charges)]).and_then(|c| c.transform(0)),
    )
}
fn string_to_bond_types(text: &str) -> Vec<BondOrder> {
    transform_from_fields("tokens", "", text, "")
        .expect("bond tokens")
        .bond_types()
        .to_vec()
}
fn string_to_charges(text: &str) -> Result<Vec<i32>, TautomerTransformError> {
    transform_from_fields("tokens", "", "", text).map(|t| t.charges().to_vec())
}
fn transform_from_line(line: &str) -> Result<Option<TautomerTransform>, TautomerTransformError> {
    transform_error(
        TautomerCatalog::from_reader(&mut Cursor::new(line.as_bytes()), -1)
            .map(|c| c.transforms().first().cloned()),
    )
}
fn read_transforms(
    reader: &mut impl BufRead,
    n: i32,
) -> Result<Vec<TautomerTransform>, TautomerCatalogError> {
    TautomerCatalog::from_reader(reader, n).map(|c| c.transforms().to_vec())
}
fn read_transforms_from_file(
    path: impl AsRef<Path>,
) -> Result<Vec<TautomerTransform>, TautomerCatalogError> {
    TautomerCatalog::from_file(path).map(|c| c.transforms().to_vec())
}
fn read_transforms_from_definitions(
    rows: &[TautomerTransformDefinition<'_>],
) -> Result<Vec<TautomerTransform>, TautomerCatalogError> {
    TautomerCatalog::from_data(rows).map(|c| c.transforms().to_vec())
}
fn current_builtin_transforms() -> Result<Vec<TautomerTransform>, TautomerCatalogError> {
    TautomerCatalog::current().map(|c| c.transforms().to_vec())
}
fn v1_builtin_transforms() -> Result<Vec<TautomerTransform>, TautomerCatalogError> {
    TautomerCatalog::v1().map(|c| c.transforms().to_vec())
}

fn query(smarts: &str) -> QueryGraph {
    parse_smarts(smarts, &SmartsParseParams::default()).expect("compile transform SMARTS")
}

#[test]
fn catalog_tokens_map_every_source_bond_symbol_in_order() {
    assert_eq!(
        string_to_bond_types("-=#:"),
        vec![
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Triple,
            BondOrder::Aromatic,
        ]
    );
}

#[test]
fn catalog_tokens_map_every_source_charge_symbol_in_order() {
    assert_eq!(
        string_to_charges("+0-").expect("valid charges"),
        &[1, 0, -1]
    );
}

#[test]
fn catalog_tokens_accept_empty_bond_and_charge_fields() {
    assert!(string_to_bond_types("").is_empty());
    assert!(string_to_charges("").expect("empty charges").is_empty());
}

#[test]
fn catalog_tokens_bond_parser_ignores_whitespace_and_unknown_bytes() {
    assert_eq!(
        string_to_bond_types(" - x\t=\n?#:µ"),
        vec![
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Triple,
            BondOrder::Aromatic,
        ]
    );
}

#[test]
fn catalog_tokens_charge_parser_rejects_whitespace() {
    assert_eq!(
        string_to_charges("+ -"),
        Err(TautomerTransformError::ChargeSymbolNotRecognised)
    );
}

#[test]
fn catalog_tokens_charge_parser_rejects_invalid_ascii_and_non_ascii_bytes() {
    for invalid in ["1", "9", "x", "µ"] {
        assert_eq!(
            string_to_charges(invalid),
            Err(TautomerTransformError::ChargeSymbolNotRecognised),
            "invalid token {invalid:?}"
        );
    }
}

#[test]
fn catalog_tokens_charge_parser_does_not_accept_signed_integer_spellings() {
    assert_eq!(
        string_to_charges("+1"),
        Err(TautomerTransformError::ChargeSymbolNotRecognised)
    );
    assert_eq!(
        string_to_charges("-1"),
        Err(TautomerTransformError::ChargeSymbolNotRecognised)
    );
}

#[test]
fn catalog_tokens_charge_parser_returns_at_first_invalid_byte() {
    assert_eq!(
        string_to_charges("+0x-"),
        Err(TautomerTransformError::ChargeSymbolNotRecognised)
    );
}

#[test]
fn transform_lines_skip_empty_comment_and_too_short_lines() {
    for line in ["", "// comment", " ", "only-a-name"] {
        assert!(
            transform_from_line(line)
                .expect("skipped line is not an error")
                .is_none(),
            "line {line:?}"
        );
    }
}

#[test]
fn transform_lines_parse_two_columns_with_source_defaults() {
    let transform = transform_from_line("Simple name\t[O]-[C]")
        .expect("valid line")
        .expect("transform");

    assert_eq!(transform.name(), "Simplename");
    assert_eq!(transform.query().num_atoms(), 2);
    assert!(transform.bond_types().is_empty());
    assert!(transform.charges().is_empty());
}

#[test]
fn transform_lines_parse_three_columns_and_remove_ascii_spaces() {
    let transform = transform_from_line("Bond edit\t[O] - [C] = [C]\t= -")
        .expect("valid line")
        .expect("transform");

    assert_eq!(transform.name(), "Bondedit");
    assert_eq!(
        transform.bond_types(),
        &[BondOrder::Double, BondOrder::Single]
    );
    assert!(transform.charges().is_empty());
}

#[test]
fn transform_lines_parse_four_columns_in_source_order() {
    let transform = transform_from_line("Charged\t[N]-[C]=[O]\t= -\t+ 0 -")
        .expect("valid line")
        .expect("transform");

    assert_eq!(
        transform.bond_types(),
        &[BondOrder::Double, BondOrder::Single]
    );
    assert_eq!(transform.charges(), &[1, 0, -1]);
}

#[test]
fn transform_lines_collapse_omitted_optional_tab_columns_like_boost_tokenizer() {
    let transform = transform_from_line("Collapsed\t[O]-[C]\t\t+")
        .expect("source treats the fourth token as the third nonempty token")
        .expect("transform");

    assert!(transform.bond_types().is_empty());
    assert!(transform.charges().is_empty());
}

#[test]
fn transform_lines_do_not_treat_indented_slashes_as_a_comment() {
    let error = transform_from_line(" // name\tinvalid smarts")
        .expect_err("only columns beginning with two slashes are comments");

    assert!(matches!(
        error,
        TautomerTransformError::CannotParseSmarts { .. }
    ));
}

#[test]
fn transform_lines_report_malformed_smarts_with_the_source_text() {
    let error =
        transform_from_line("broken\t[C").expect_err("malformed SMARTS must remain visible");

    assert!(matches!(
        error,
        TautomerTransformError::CannotParseSmarts { ref smarts, .. } if smarts == "[C"
    ));
}

#[test]
fn transform_lines_report_invalid_charge_symbols_before_smarts_parsing() {
    let error = transform_from_fields("broken", "[", "", "x")
        .expect_err("source converts charge tokens before parsing SMARTS");

    assert_eq!(error, TautomerTransformError::ChargeSymbolNotRecognised);
}

#[test]
fn transform_lines_over_four_columns_store_the_source_empty_transform() {
    // Original3483 expected MissingDonorOrAcceptor{actual:0}; kept in archived source.
    // TautomerCatalogUtils.cpp40–86 leaves fields empty. SmilesParse.cpp247–249
    // returns an empty query; ROOT8426 approves the source constructor result.
    let transform = transform_from_line("one\ttwo\tthree\tfour\tfive")
        .expect("source empty query")
        .expect("transform");
    assert_eq!(transform.name(), "");
    assert_eq!(transform.query().num_atoms(), 0);
    assert_eq!(transform.query().num_bonds(), 0);
    assert!(transform.bond_types().is_empty());
    assert!(transform.charges().is_empty());
}

#[test]
fn catalog_reader_reads_stream_comments_blank_lines_and_early_eof() {
    let data = "// header\n\nfirst\t[O]-[C]\n // not a comment\nsecond\t[N]-[C]";
    let mut reader = Cursor::new(data.as_bytes());
    let transforms = read_transforms(&mut reader, -1).expect("read stream");

    assert_eq!(
        transforms
            .iter()
            .map(TautomerTransform::name)
            .collect::<Vec<_>>(),
        ["first", "second"]
    );
}

#[test]
fn catalog_reader_bounded_stream_counts_only_emitted_transforms() {
    let data = "// header\nfirst\t[O]-[C]\n\n// middle\nsecond\t[N]-[C]\nthird\t[S]-[C]\n";
    let mut reader = Cursor::new(data.as_bytes());
    let transforms = read_transforms(&mut reader, 2).expect("read bounded stream");

    assert_eq!(
        transforms
            .iter()
            .map(TautomerTransform::name)
            .collect::<Vec<_>>(),
        ["first", "second"]
    );
}

#[test]
fn catalog_reader_zero_limit_reads_nothing() {
    let mut reader = Cursor::new(b"first\t[O]-[C]\n".as_slice());
    assert!(
        read_transforms(&mut reader, 0)
            .expect("zero limit")
            .is_empty()
    );
    assert_eq!(reader.position(), 0);
}

#[test]
fn catalog_reader_file_entry_preserves_source_order() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("tautomers.in");
    std::fs::write(&path, "first\t[O]-[C]\nsecond\t[N]-[C]\nthird\t[S]-[C]\n")
        .expect("write fixture");

    let transforms = read_transforms_from_file(&path).expect("read file");
    assert_eq!(
        transforms
            .iter()
            .map(TautomerTransform::name)
            .collect::<Vec<_>>(),
        ["first", "second", "third"]
    );
}

#[test]
fn catalog_reader_embedded_definitions_use_the_same_constructor() {
    let definitions = [
        ("first", "[O]-[C]", "=", "+0"),
        ("second", "[N]-[C]", "-", "0-"),
    ];
    let transforms =
        read_transforms_from_definitions(&definitions).expect("read embedded definitions");

    assert_eq!(transforms[0].name(), "first");
    assert_eq!(transforms[0].bond_types(), &[BondOrder::Double]);
    assert_eq!(transforms[0].charges(), &[1, 0]);
    assert_eq!(transforms[1].name(), "second");
    assert_eq!(transforms[1].bond_types(), &[BondOrder::Single]);
    assert_eq!(transforms[1].charges(), &[0, -1]);
}

#[test]
fn catalog_reader_missing_file_returns_structured_io_error() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("missing.in");
    let error = read_transforms_from_file(&path).expect_err("missing file must fail");

    assert!(matches!(
        error,
        TautomerCatalogError::BadInputFile { path: actual, .. } if actual == path
    ));
}

struct BrokenReader;

impl Read for BrokenReader {
    fn read(&mut self, _buffer: &mut [u8]) -> io::Result<usize> {
        Err(io::Error::other("broken catalog stream"))
    }
}

impl BufRead for BrokenReader {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        Err(io::Error::other("broken catalog stream"))
    }

    fn consume(&mut self, _amount: usize) {}
}

#[test]
fn catalog_reader_stream_failure_returns_structured_io_error() {
    let error = read_transforms(&mut BrokenReader, -1).expect_err("bad stream must fail");
    assert!(matches!(error, TautomerCatalogError::BadStreamContents(_)));
}

#[test]
fn catalog_reader_rejects_non_utf8_catalog_data() {
    let mut reader = Cursor::new([0xff, b'\n']);
    let error = read_transforms(&mut reader, -1).expect_err("catalog text must be UTF-8");
    assert!(matches!(error, TautomerCatalogError::InvalidUtf8(_)));
}

#[test]
fn catalog_reader_source_line_limit_sets_fail_state_after_truncated_line() {
    // Archived3610 used511 but getline accepts511 followed by newline.
    // Standard-library probe retained; 512 is the source failbit case.
    let mut bytes = vec![b'/'; 512];
    bytes.extend_from_slice(b"\nvalid\t[O]-[C]\n");
    let mut reader = Cursor::new(bytes);

    let transforms = read_transforms(&mut reader, -1).expect("read truncated source line");
    assert!(transforms.is_empty());
}

fn assert_builtin_catalog(
    actual: &[TautomerTransform],
    expected: &[TautomerTransformDefinition<'_>],
) {
    assert_eq!(actual.len(), expected.len());
    for (index, (transform, &(name, smarts, bonds, charges))) in
        actual.iter().zip(expected).enumerate()
    {
        assert_eq!(transform.name(), name, "name at source row {index}");
        let expected_query = query(smarts);
        let plan = CompiledQuery::compile(transform.query().clone())
            .expect("each of73 real compiled plans");
        assert_eq!(plan.num_atoms(), expected_query.num_atoms());
        assert_eq!(plan.num_bonds(), expected_query.num_bonds());
        assert_eq!(plan.atom_order().len(), expected_query.num_atoms());
        assert_eq!(
            transform.query().atoms(),
            expected_query.atoms(),
            "compiled query atoms at source row {index}: {name}"
        );
        assert_eq!(
            transform.query().bonds(),
            expected_query.bonds(),
            "compiled query bonds at source row {index}: {name}"
        );
        let expected_bonds: Vec<_> = bonds
            .bytes()
            .filter_map(|byte| match byte {
                b'-' => Some(BondOrder::Single),
                b'=' => Some(BondOrder::Double),
                b'#' => Some(BondOrder::Triple),
                b':' => Some(BondOrder::Aromatic),
                _ => None,
            })
            .collect();
        assert_eq!(
            transform.bond_types(),
            expected_bonds,
            "bond edits at source row {index}: {name}"
        );
        let expected_charges: Vec<_> = charges
            .bytes()
            .map(|byte| match byte {
                b'+' => 1,
                b'0' => 0,
                b'-' => -1,
                _ => panic!("invalid generated charge byte at source row {index}"),
            })
            .collect();
        assert_eq!(
            transform.charges(),
            expected_charges,
            "charge edits at source row {index}: {name}"
        );
    }
}

#[test]
fn builtin_catalogs_current_contains_every_compiled_source_row_in_order() {
    assert_eq!(CURRENT_TAUTOMER_TRANSFORM_DEFINITIONS.len(), 37);
    let transforms = current_builtin_transforms().expect("compile current catalog");
    assert_builtin_catalog(&transforms, CURRENT_TAUTOMER_TRANSFORM_DEFINITIONS);
}

#[test]
fn builtin_catalogs_v1_contains_every_compiled_source_row_in_order() {
    assert_eq!(V1_TAUTOMER_TRANSFORM_DEFINITIONS.len(), 36);
    let transforms = v1_builtin_transforms().expect("compile V1 catalog");
    assert_builtin_catalog(&transforms, V1_TAUTOMER_TRANSFORM_DEFINITIONS);
}

#[test]
fn catalog_object_empty_default_and_empty_file_current_are_distinct() {
    // Original3686 Default/current expectation retained in immutable archive.
    // TautomerCatalogParams.h no-arg ctor clears transforms; .cpp file="" selects37.
    let empty = TautomerCatalog::empty();
    let default = TautomerCatalog::default();
    let current = TautomerCatalog::from_file("").expect("current");
    assert!(empty.transforms().is_empty());
    assert_eq!(default, empty);
    assert_eq!(current.transforms().len(), 37);
    assert_eq!(
        current,
        TautomerCatalog::current().expect("explicit current")
    );
}

#[test]
fn catalog_object_constructs_from_file_data_and_v1_sources() {
    let directory = tempfile::tempdir().expect("temporary directory");
    let path = directory.path().join("custom.in");
    std::fs::write(&path, "file-first\t[O]-[C]\nfile-second\t[N]-[C]\n").expect("write catalog");
    let file_catalog = TautomerCatalog::from_file(&path).expect("file catalog");
    let data_catalog =
        TautomerCatalog::from_data(&[("data", "[S]-[C]", "=", "+0")]).expect("data catalog");
    let v1_catalog = TautomerCatalog::v1().expect("V1 catalog");

    assert_eq!(
        file_catalog
            .transforms()
            .iter()
            .map(TautomerTransform::name)
            .collect::<Vec<_>>(),
        ["file-first", "file-second"]
    );
    assert_eq!(data_catalog.transforms()[0].name(), "data");
    assert_eq!(
        data_catalog.transforms()[0].bond_types(),
        &[BondOrder::Double]
    );
    assert_eq!(data_catalog.transforms()[0].charges(), &[1, 0]);
    assert_eq!(v1_catalog.transforms().len(), 36);
}

#[test]
fn catalog_object_indexing_returns_an_independent_value_and_checks_bounds() {
    let catalog = TautomerCatalog::from_data(&[("one", "[O]-[C]", "-", "")]).expect("catalog");
    let returned = catalog.transform(0).expect("returned");
    let mut edited = returned.query().clone();
    edited.atom_mut(0).expect("atom").set_formal_charge(-1);
    let edited = TautomerTransform::new(
        returned.name(),
        edited,
        vec![BondOrder::Double],
        returned.charges().to_vec(),
    )
    .expect("edited independent value");
    assert_eq!(edited.query().atoms()[0].formal_charge(), -1);
    assert_eq!(
        catalog.transforms()[0].query().atoms()[0].formal_charge(),
        0
    );
    assert_eq!(catalog.transforms()[0].bond_types(), &[BondOrder::Single]);
    assert!(matches!(
        catalog.transform(1),
        Err(TautomerCatalogError::TransformIndexOutOfRange { index: 1, len: 1 })
    ));
}

#[test]
fn catalog_object_clone_and_clone_from_have_independent_value_semantics() {
    let source = TautomerCatalog::from_data(&[("source", "[O]-[C]", "-", "")]).expect("source");
    let cloned = source.clone();
    let returned = cloned.transform(0).expect("clone row");
    let changed = TautomerTransform::new(
        returned.name(),
        returned.query().clone(),
        vec![BondOrder::Double],
        returned.charges().to_vec(),
    )
    .expect("changed row");
    assert_eq!(changed.bond_types(), &[BondOrder::Double]);
    assert_eq!(source.transforms()[0].bond_types(), &[BondOrder::Single]);
    let mut target = TautomerCatalog::from_data(&[("target", "[N]-[C]", "=", "")]).expect("target");
    target.clone_from(&source);
    assert_eq!(target, source);
}

#[test]
fn catalog_object_serialization_is_exactly_count_and_newline() {
    let catalog =
        TautomerCatalog::from_data(&[("first", "[O]-[C]", "", ""), ("second", "[N]-[C]", "", "")])
            .expect("catalog");
    let mut bytes = Vec::new();

    catalog.write_to(&mut bytes).expect("write catalog count");
    assert_eq!(bytes, b"2\n");
    assert_eq!(catalog.serialize(), "2\n");
}

#[test]
fn catalog_object_deserialization_is_a_structured_upstream_boundary() {
    assert!(matches!(
        TautomerCatalog::deserialize("37\n"),
        Err(TautomerCatalogError::DeserializationUnderConstruction)
    ));
}

#[test]
fn transform_construction_preserves_name_query_and_source_ordered_edits() {
    let transform = TautomerTransform::new(
        "ordered",
        query("[O]-[C]=[C]"),
        vec![BondOrder::Double, BondOrder::Single],
        vec![-1, 0, 1],
    )
    .expect("valid transform");

    assert_eq!(transform.name(), "ordered");
    assert_eq!(transform.query().num_atoms(), 3);
    assert_eq!(transform.query().num_bonds(), 2);
    assert_eq!(
        transform.bond_types(),
        &[BondOrder::Double, BondOrder::Single]
    );
    assert_eq!(transform.charges(), &[-1, 0, 1]);
}

#[test]
fn transform_clone_has_independent_value_semantics() {
    let original = TautomerTransform::new(
        "original",
        query("[O]-[C]=[C]"),
        vec![BondOrder::Double, BondOrder::Single],
        Vec::new(),
    )
    .expect("valid");
    let cloned = original.clone();
    let mut edited = cloned.query().clone();
    edited.atom_mut(0).expect("atom").set_formal_charge(-1);
    let cloned = TautomerTransform::new(
        "clone",
        edited,
        vec![BondOrder::Triple, BondOrder::Single],
        cloned.charges().to_vec(),
    )
    .expect("independent change");
    assert_eq!(original.name(), "original");
    assert_eq!(original.query().atoms()[0].formal_charge(), 0);
    assert_eq!(original.bond_types()[0], BondOrder::Double);
    assert_eq!(cloned.name(), "clone");
    assert_eq!(cloned.query().atoms()[0].formal_charge(), -1);
    assert_eq!(cloned.bond_types()[0], BondOrder::Triple);
}

#[test]
fn transform_clone_from_replaces_every_source_owned_field() {
    let source = TautomerTransform::new(
        "source",
        query("[N]-[C]=[O]"),
        vec![BondOrder::Double, BondOrder::Single],
        vec![1, 0, -1],
    )
    .expect("valid source transform");
    let mut target = TautomerTransform::new("target", query("[O]-[C]=[C]"), Vec::new(), Vec::new())
        .expect("valid target transform");

    target.clone_from(&source);

    assert_eq!(target, source);
}

#[test]
fn transform_keeps_and_reuses_the_compiled_query_value() {
    let compiled = query("[O]-[C]=[C]").with_prop("compiled-sentinel", "present");
    let shared_atoms = compiled.atoms().as_ptr();

    let transform = TautomerTransform::new("compiled", compiled, Vec::new(), Vec::new())
        .expect("valid transform");

    assert_eq!(transform.query().prop("compiled-sentinel"), Some("present"));
    assert_eq!(transform.query().atoms().as_ptr(), shared_atoms);
}

#[test]
fn transform_accepts_empty_bond_and_charge_edit_vectors() {
    let transform =
        TautomerTransform::new("alternating", query("[O]-[C]=[C]"), Vec::new(), Vec::new())
            .expect("empty edits select source defaults");

    assert!(transform.bond_types().is_empty());
    assert!(transform.charges().is_empty());
}

#[test]
fn transform_accepts_source_empty_and_single_atom_queries() {
    // Original6641 one-atom rejection retained. Params.h constructor has no >=2 guard.
    for (smarts, n) in [("", 0), ("[O]", 1)] {
        let transform = TautomerTransform::new("allowed", query(smarts), Vec::new(), Vec::new())
            .expect("source stores query");
        assert_eq!(transform.query().num_atoms(), n);
        assert!(transform.bond_types().is_empty());
        assert!(transform.charges().is_empty());
    }
}

#[test]
fn transform_stores_non_aligned_source_bond_edits() {
    // Original6652 BondEditCount rejection retained; Params.h simply moves the vector.
    let transform = TautomerTransform::new(
        "unaligned",
        query("[O]-[C]=[C]"),
        vec![BondOrder::Double],
        Vec::new(),
    )
    .expect("source constructor");
    assert_eq!(transform.query().num_bonds(), 2);
    assert_eq!(transform.bond_types(), &[BondOrder::Double]);
}

#[test]
fn transform_stores_non_aligned_source_charge_edits() {
    // Original6671 ChargeEditCount rejection retained; Params.h simply moves the vector.
    let transform =
        TautomerTransform::new("unaligned", query("[O]-[C]=[C]"), Vec::new(), vec![1, -1])
            .expect("source constructor");
    assert_eq!(transform.query().num_atoms(), 3);
    assert_eq!(transform.charges(), &[1, -1]);
}
#[test]
fn exactly_511_payload_with_newline_continues_to_next_line() {
    let mut bytes = vec![b'/'; 511];
    bytes.extend_from_slice(b"\nvalid\t[O]-[C]\n");
    let mut r = Cursor::new(bytes);
    let c = TautomerCatalog::from_reader(&mut r, -1).expect("getline source");
    assert_eq!(c.transforms().len(), 1);
    assert_eq!(c.transforms()[0].name(), "valid");
}
#[test]
fn exactly_511_payload_at_eof_terminates_without_extra_transform() {
    let mut r = Cursor::new(vec![b'/'; 511]);
    let c = TautomerCatalog::from_reader(&mut r, -1).expect("EOF");
    assert!(c.transforms().is_empty());
    assert_eq!(r.position(), 511);
}
#[test]
fn overlong_line_parses_truncated_prefix_and_leaves_unread_suffix() {
    let mut bytes = b"first\t[O]-[C]".to_vec();
    bytes.resize(512, b' ');
    bytes.extend_from_slice(b"\nsecond\t[N]-[C]\n");
    let mut r = Cursor::new(bytes);
    let c = TautomerCatalog::from_reader(&mut r, -1).expect("truncated source prefix");
    assert_eq!(c.transforms().len(), 1);
    assert_eq!(c.transforms()[0].name(), "first");
    assert_eq!(r.position(), 511);
}
#[test]
fn source_c_string_nul_truncates_parse_but_consumes_full_physical_line() {
    let mut bytes = b"first\t[O]-[C]\0ignored".to_vec();
    bytes.push(0xff);
    bytes.extend_from_slice(b"\nsecond\t[N]-[C]\n");
    let mut r = Cursor::new(bytes);
    let c = TautomerCatalog::from_reader(&mut r, -1).expect("source invisible suffix");
    assert_eq!(
        c.transforms()
            .iter()
            .map(TautomerTransform::name)
            .collect::<Vec<_>>(),
        ["first", "second"]
    );
}
#[test]
fn direct_tuple_fields_keep_spaces_and_source_name_metadata() {
    let c = TautomerCatalog::from_data(&[("name with spaces", "[O]-[C]", "-", "")]).expect("tuple");
    assert_eq!(c.transforms()[0].name(), "name with spaces");
    assert_eq!(
        c.transforms()[0].query().prop("_Name"),
        Some("name with spaces")
    );
}
#[test]
fn empty_definition_array_is_an_empty_catalog() {
    assert_eq!(
        TautomerCatalog::from_data(&[]).expect("empty data"),
        TautomerCatalog::empty()
    );
}
#[test]
fn catalog_stream_deserialization_is_source_unimplemented_and_does_not_consume() {
    let mut r = Cursor::new(b"37\n");
    assert!(matches!(
        TautomerCatalog::deserialize_from_reader(&mut r),
        Err(TautomerCatalogError::DeserializationUnderConstruction)
    ));
    assert_eq!(r.position(), 0);
}
#[test]
fn zero_stream_limit_does_not_observe_unread_io_failure() {
    let c = TautomerCatalog::from_reader(&mut BrokenReader, 0)
        .expect("BufRead has no preexisting badbit; no read");
    assert!(c.transforms().is_empty());
}
#[test]
fn count_stream_write_propagates_io_error() {
    struct BrokenWriter;
    impl io::Write for BrokenWriter {
        fn write(&mut self, _: &[u8]) -> io::Result<usize> {
            Err(io::Error::new(
                io::ErrorKind::BrokenPipe,
                "broken count stream",
            ))
        }
        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }
    let error = TautomerCatalog::empty()
        .write_to(&mut BrokenWriter)
        .expect_err("write error");
    assert_eq!(error.kind(), io::ErrorKind::BrokenPipe);
}
#[test]
fn stream_error_retains_the_original_io_cause() {
    use std::error::Error;
    let error = TautomerCatalog::from_reader(&mut BrokenReader, -1).expect_err("IO");
    assert_eq!(
        error.source().expect("cause").to_string(),
        "broken catalog stream"
    );
}
