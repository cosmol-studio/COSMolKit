use cosmolkit_core::property_value_to_string;
use cosmolkit_io::{
    MolBlockRecord, SdfDataReadParams, SdfReadError, read_sdf_graph_record_detached,
    read_sdf_record_detached, read_sdf_record_detached_with_params,
};
use cosmolkit_model::{MoleculeProperties, PropertyValue, SdfPropertyListTarget};

fn v2000_atom(symbol: &str) -> String {
    format!("    0.0000    0.0000    0.0000 {symbol:<3} 0  0  0  0  0  0  0  0  0  0  0  0")
}

fn concrete_record(fields: &str) -> String {
    format!(
        concat!(
            "property lists\n  COSMolKit         2D\n\n",
            "  3  2  0  0  0  0  0  0  0  0999 V2000\n",
            "{}\n{}\n{}\n",
            "  1  2  1  0  0  0  0\n",
            "  2  3  1  0  0  0  0\n",
            "M  END\n{}$$$$\n",
        ),
        v2000_atom("C"),
        v2000_atom("N"),
        v2000_atom("O"),
        fields,
    )
}

fn zero_record(fields: &str) -> String {
    format!(
        concat!(
            "zero property lists\n  COSMolKit         2D\n\n",
            "  0  0  0  0  0  0  0  0  0  0999 V2000\n",
            "M  END\n{}$$$$\n",
        ),
        fields,
    )
}

fn assert_list(
    properties: &MoleculeProperties,
    index: usize,
    target: SdfPropertyListTarget,
    name: &str,
    values: &[Option<PropertyValue>],
) {
    let list = &properties.sdf_property_lists()[index];
    assert_eq!(list.target(), target);
    assert_eq!(list.name(), name);
    assert_eq!(list.values(), values);
}

fn assert_typed_property(actual: Option<&PropertyValue>, expected: PropertyValue, projected: &str) {
    assert_eq!(actual, Some(&expected));
    assert_eq!(
        property_value_to_string(actual.expect("typed property must exist")),
        Ok(projected.to_owned())
    );
}

#[test]
fn sdf_props_all_eight_prefixes_apply_typed_values_and_continue_after_bad_items() {
    // Pinned applyMolListProp dispatches the eight case-sensitive prefixes,
    // leaves missing and bad lexical items absent, and continues later rows.
    let input = concrete_record(concat!(
        ">  <atom.prop.Text>\nalpha n/a omega\n\n",
        ">  <atom.iprop.Integer>\n+7 invalid -3\n\n",
        ">  <atom.dprop.Real>\n1e2 broken -5e-1\n\n",
        ">  <atom.bprop.Flag>\n1 0 1\n\n",
        ">  <bond.prop.Text>\nleft right\n\n",
        ">  <bond.iprop.Integer>\n-9 +11\n\n",
        ">  <bond.dprop.Real>\n6.25e-1 invalid\n\n",
        ">  <bond.bprop.Flag>\n0 1\n\n",
    ));
    let record = read_sdf_record_detached(&input).expect("all property-list kinds");

    assert_typed_property(
        record.topology.atoms[0].prop("Text"),
        PropertyValue::String("alpha".to_owned()),
        "alpha",
    );
    assert_eq!(record.topology.atoms[1].prop("Text"), None);
    assert_typed_property(
        record.topology.atoms[2].prop("Text"),
        PropertyValue::String("omega".to_owned()),
        "omega",
    );
    assert_typed_property(
        record.topology.atoms[0].prop("Integer"),
        PropertyValue::Int(7),
        "7",
    );
    assert_eq!(record.topology.atoms[1].prop("Integer"), None);
    assert_typed_property(
        record.topology.atoms[2].prop("Integer"),
        PropertyValue::Int(-3),
        "-3",
    );
    assert_typed_property(
        record.topology.atoms[0].prop("Real"),
        PropertyValue::Double(100.0),
        "100",
    );
    assert_eq!(record.topology.atoms[1].prop("Real"), None);
    assert_typed_property(
        record.topology.atoms[2].prop("Real"),
        PropertyValue::Double(-0.5),
        "-0.5",
    );
    assert_typed_property(
        record.topology.atoms[0].prop("Flag"),
        PropertyValue::Bool(true),
        "1",
    );
    assert_typed_property(
        record.topology.atoms[1].prop("Flag"),
        PropertyValue::Bool(false),
        "0",
    );
    assert_typed_property(
        record.topology.atoms[2].prop("Flag"),
        PropertyValue::Bool(true),
        "1",
    );

    assert_typed_property(
        record.topology.bonds[0].prop("Text"),
        PropertyValue::String("left".to_owned()),
        "left",
    );
    assert_typed_property(
        record.topology.bonds[1].prop("Text"),
        PropertyValue::String("right".to_owned()),
        "right",
    );
    assert_typed_property(
        record.topology.bonds[0].prop("Integer"),
        PropertyValue::Int(-9),
        "-9",
    );
    assert_typed_property(
        record.topology.bonds[1].prop("Integer"),
        PropertyValue::Int(11),
        "11",
    );
    assert_typed_property(
        record.topology.bonds[0].prop("Real"),
        PropertyValue::Double(0.625),
        "0.625",
    );
    assert_eq!(record.topology.bonds[1].prop("Real"), None);
    assert_typed_property(
        record.topology.bonds[0].prop("Flag"),
        PropertyValue::Bool(false),
        "0",
    );
    assert_typed_property(
        record.topology.bonds[1].prop("Flag"),
        PropertyValue::Bool(true),
        "1",
    );

    let lists = record.properties.sdf_property_lists();
    assert_eq!(lists.len(), 8);
    assert_list(
        &record.properties,
        0,
        SdfPropertyListTarget::Atom,
        "Text",
        &[
            Some(PropertyValue::String("alpha".to_owned())),
            None,
            Some(PropertyValue::String("omega".to_owned())),
        ],
    );
    assert_list(
        &record.properties,
        1,
        SdfPropertyListTarget::Atom,
        "Integer",
        &[
            Some(PropertyValue::Int(7)),
            None,
            Some(PropertyValue::Int(-3)),
        ],
    );
    assert_list(
        &record.properties,
        2,
        SdfPropertyListTarget::Atom,
        "Real",
        &[
            Some(PropertyValue::Double(100.0)),
            None,
            Some(PropertyValue::Double(-0.5)),
        ],
    );
    assert_list(
        &record.properties,
        3,
        SdfPropertyListTarget::Atom,
        "Flag",
        &[
            Some(PropertyValue::Bool(true)),
            Some(PropertyValue::Bool(false)),
            Some(PropertyValue::Bool(true)),
        ],
    );
    assert_list(
        &record.properties,
        4,
        SdfPropertyListTarget::Bond,
        "Text",
        &[
            Some(PropertyValue::String("left".to_owned())),
            Some(PropertyValue::String("right".to_owned())),
        ],
    );
    assert_list(
        &record.properties,
        5,
        SdfPropertyListTarget::Bond,
        "Integer",
        &[Some(PropertyValue::Int(-9)), Some(PropertyValue::Int(11))],
    );
    assert_list(
        &record.properties,
        6,
        SdfPropertyListTarget::Bond,
        "Real",
        &[Some(PropertyValue::Double(0.625)), None],
    );
    assert_list(
        &record.properties,
        7,
        SdfPropertyListTarget::Bond,
        "Flag",
        &[
            Some(PropertyValue::Bool(false)),
            Some(PropertyValue::Bool(true)),
        ],
    );
    assert_eq!(
        record.properties.prop("atom.iprop.Integer"),
        Some("+7 invalid -3")
    );
}

#[test]
fn sdf_props_numeric_lexemes_retain_values_and_source_projection_separately() {
    let input = concrete_record(concat!(
        ">  <atom.iprop.Number>\n+7 007 -3\n\n",
        ">  <atom.dprop.Real>\n1.00 1e2 -5e-1\n\n",
    ));
    let record = read_sdf_record_detached(&input).expect("typed numeric spellings");

    for (row, expected, projected) in [(0, 7, "7"), (1, 7, "7"), (2, -3, "-3")] {
        assert_typed_property(
            record.topology.atoms[row].prop("Number"),
            PropertyValue::Int(expected),
            projected,
        );
    }
    for (row, expected, projected) in [(0, 1.0, "1"), (1, 100.0, "100"), (2, -0.5, "-0.5")] {
        assert_typed_property(
            record.topology.atoms[row].prop("Real"),
            PropertyValue::Double(expected),
            projected,
        );
    }
    assert_eq!(
        record.properties.prop("atom.iprop.Number"),
        Some("+7 007 -3")
    );
    assert_eq!(
        record.properties.prop("atom.dprop.Real"),
        Some("1.00 1e2 -5e-1")
    );
    assert_eq!(
        record.properties.sdf_data_fields(),
        &[
            ("atom.iprop.Number".to_owned(), "+7 007 -3".to_owned()),
            ("atom.dprop.Real".to_owned(), "1.00 1e2 -5e-1".to_owned(),),
        ]
    );
}

#[test]
fn sdf_props_custom_empty_markers_and_compressed_multiline_delimiters_preserve_rows() {
    let input = concrete_record(concat!(
        ">  <atom.prop.Custom>\n[?]\tfirst\t?\n third\n\n",
        ">  <atom.prop.EmptyMarker>\n[]  one\ttwo three\n\n",
    ));
    let record = read_sdf_record_detached(&input).expect("custom missing markers");

    assert_list(
        &record.properties,
        0,
        SdfPropertyListTarget::Atom,
        "Custom",
        &[
            Some(PropertyValue::String("first".to_owned())),
            None,
            Some(PropertyValue::String("third".to_owned())),
        ],
    );
    assert_list(
        &record.properties,
        1,
        SdfPropertyListTarget::Atom,
        "EmptyMarker",
        &[
            Some(PropertyValue::String("one".to_owned())),
            Some(PropertyValue::String("two".to_owned())),
            Some(PropertyValue::String("three".to_owned())),
        ],
    );
    assert_typed_property(
        record.topology.atoms[2].prop("Custom"),
        PropertyValue::String("third".to_owned()),
        "third",
    );
    assert_eq!(
        record.properties.prop("atom.prop.Custom"),
        Some("[?]\tfirst\t?\n third")
    );
}

#[test]
fn sdf_props_query_records_keep_query_carrier_and_apply_atom_and_bond_lists() {
    let input = concat!(
        "query properties\n  COSMolKit         2D\n\n",
        "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
        "M  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\n",
        "M  V30 BEGIN ATOM\n",
        "M  V30 1 * 0 0 0 0\nM  V30 2 O 1 0 0 0\n",
        "M  V30 END ATOM\nM  V30 BEGIN BOND\n",
        "M  V30 1 5 1 2\n",
        "M  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ">  <atom.prop.Text>\nquery n/a\n\n",
        ">  <atom.iprop.Count>\n+7 invalid\n\n",
        ">  <atom.dprop.Weight>\n2.5 invalid\n\n",
        ">  <atom.bprop.Active>\n0 1\n\n",
        ">  <bond.prop.Text>\nedge\n\n",
        ">  <bond.iprop.Count>\n007\n\n",
        ">  <bond.dprop.Weight>\n-5e-1\n\n",
        ">  <bond.bprop.Selected>\n1\n\n$$$$\n",
    );
    let record = read_sdf_graph_record_detached(input).expect("query property lists");
    let MolBlockRecord::Query(query) = record.mol_block else {
        panic!("property expansion must not lower the query graph");
    };

    assert_typed_property(
        query.query.atoms()[0].prop("Text"),
        PropertyValue::String("query".to_owned()),
        "query",
    );
    assert_eq!(query.query.atoms()[1].prop("Text"), None);
    assert_typed_property(
        query.query.atoms()[0].prop("Count"),
        PropertyValue::Int(7),
        "7",
    );
    assert_eq!(query.query.atoms()[1].prop("Count"), None);
    assert_typed_property(
        query.query.atoms()[0].prop("Weight"),
        PropertyValue::Double(2.5),
        "2.5",
    );
    assert_eq!(query.query.atoms()[1].prop("Weight"), None);
    assert_typed_property(
        query.query.atoms()[0].prop("Active"),
        PropertyValue::Bool(false),
        "0",
    );
    assert_typed_property(
        query.query.atoms()[1].prop("Active"),
        PropertyValue::Bool(true),
        "1",
    );
    assert_typed_property(
        query.query.bonds()[0].bond().prop("Text"),
        PropertyValue::String("edge".to_owned()),
        "edge",
    );
    assert_typed_property(
        query.query.bonds()[0].bond().prop("Count"),
        PropertyValue::Int(7),
        "7",
    );
    assert_typed_property(
        query.query.bonds()[0].bond().prop("Weight"),
        PropertyValue::Double(-0.5),
        "-0.5",
    );
    assert_typed_property(
        query.query.bonds()[0].bond().prop("Selected"),
        PropertyValue::Bool(true),
        "1",
    );
    assert_eq!(query.properties.sdf_property_lists().len(), 8);
    assert_list(
        &query.properties,
        2,
        SdfPropertyListTarget::Atom,
        "Weight",
        &[Some(PropertyValue::Double(2.5)), None],
    );
    assert_list(
        &query.properties,
        7,
        SdfPropertyListTarget::Bond,
        "Selected",
        &[Some(PropertyValue::Bool(true))],
    );
    assert_eq!(
        query.properties.prop("atom.iprop.Count"),
        Some("+7 invalid")
    );
    assert_eq!(query.properties.prop("bond.dprop.Weight"), Some("-5e-1"));
}

#[test]
fn sdf_props_strict_count_mismatches_are_structured_for_atom_and_bond_targets() {
    let atom_short =
        concrete_record(">  <atom.prop.ValidBefore>\na b c\n\n>  <atom.iprop.Short>\n1 2\n\n");
    // Source applyMolListProp warns and returns before assigning any item.
    let record =
        read_sdf_record_detached(&atom_short).expect("source count mismatch preserves raw field");
    assert_eq!(record.properties.prop("atom.iprop.Short"), Some("1 2"));
    assert_eq!(
        record.data_fields.last(),
        Some(&("atom.iprop.Short".to_owned(), "1 2".to_owned()))
    );
    assert_eq!(record.topology.atoms.len(), 3);
    assert_eq!(record.topology.bonds.len(), 2);
    assert_eq!(record.properties.sdf_property_lists().len(), 1);
    assert_list(
        &record.properties,
        0,
        SdfPropertyListTarget::Atom,
        "ValidBefore",
        &[
            Some(PropertyValue::String("a".into())),
            Some(PropertyValue::String("b".into())),
            Some(PropertyValue::String("c".into())),
        ],
    );
    assert!(
        record
            .topology
            .atoms
            .iter()
            .all(|item| item.prop("Short").is_none())
    );

    let bond_long = concrete_record(">  <bond.prop.Long>\na b c\n\n");
    // Source applyMolListProp warns and returns before assigning any item.
    let record =
        read_sdf_record_detached(&bond_long).expect("source count mismatch preserves raw field");
    assert_eq!(record.properties.prop("bond.prop.Long"), Some("a b c"));
    assert_eq!(
        record.data_fields.last(),
        Some(&("bond.prop.Long".to_owned(), "a b c".to_owned()))
    );
    assert_eq!(record.topology.atoms.len(), 3);
    assert_eq!(record.topology.bonds.len(), 2);
    assert!(record.properties.sdf_property_lists().is_empty());
    assert!(
        record
            .topology
            .bonds
            .iter()
            .all(|item| item.prop("Long").is_none())
    );
}

#[test]
fn sdf_props_nonstrict_count_mismatches_preserve_raw_fields_without_partial_expansion() {
    let input = concrete_record(concat!(
        ">  <atom.prop.Short>\nfirst second\n\n",
        ">  <bond.iprop.Long>\n1 2 3\n\n",
    ));
    let record = read_sdf_record_detached_with_params(
        &input,
        SdfDataReadParams {
            strict_parsing: false,
            ..SdfDataReadParams::default()
        },
    )
    .expect("nonstrict mismatches retain only raw fields");

    assert_eq!(
        record.properties.prop("atom.prop.Short"),
        Some("first second")
    );
    assert_eq!(record.properties.prop("bond.iprop.Long"), Some("1 2 3"));
    assert!(record.properties.sdf_property_lists().is_empty());
    assert!(
        record
            .topology
            .atoms
            .iter()
            .all(|atom| atom.prop("Short").is_none())
    );
    assert!(
        record
            .topology
            .bonds
            .iter()
            .all(|bond| bond.prop("Long").is_none())
    );
}

#[test]
fn sdf_props_zero_tables_and_boundary_empty_tokens_follow_count_rules() {
    let empty_value = zero_record(">  <atom.prop.Empty>\n\n");
    // Source applyMolListProp warns and returns before assigning any item.
    let record =
        read_sdf_record_detached(&empty_value).expect("source count mismatch preserves raw field");
    assert_eq!(record.properties.prop("atom.prop.Empty"), Some(""));
    assert_eq!(
        record.data_fields.last(),
        Some(&("atom.prop.Empty".to_owned(), "".to_owned()))
    );
    assert_eq!(record.topology.atoms.len(), 0);
    assert_eq!(record.topology.bonds.len(), 0);
    assert!(record.properties.sdf_property_lists().is_empty());
    assert!(
        record
            .topology
            .atoms
            .iter()
            .all(|item| item.prop("Empty").is_none())
    );

    // boost::split(token_compress_on) preserves boundary empty tokens. They
    // therefore participate in the source count check instead of being trim.
    let leading = concrete_record(">  <atom.prop.Leading>\n one two three\n\n");
    // Source applyMolListProp warns and returns before assigning any item.
    let record =
        read_sdf_record_detached(&leading).expect("source count mismatch preserves raw field");
    assert_eq!(
        record.properties.prop("atom.prop.Leading"),
        Some(" one two three")
    );
    assert_eq!(
        record.data_fields.last(),
        Some(&("atom.prop.Leading".to_owned(), " one two three".to_owned()))
    );
    assert_eq!(record.topology.atoms.len(), 3);
    assert_eq!(record.topology.bonds.len(), 2);
    assert!(record.properties.sdf_property_lists().is_empty());
    assert!(
        record
            .topology
            .atoms
            .iter()
            .all(|item| item.prop("Leading").is_none())
    );
    let trailing = concrete_record(">  <bond.prop.Trailing>\none two \n\n");
    // Source applyMolListProp warns and returns before assigning any item.
    let record =
        read_sdf_record_detached(&trailing).expect("source count mismatch preserves raw field");
    assert_eq!(
        record.properties.prop("bond.prop.Trailing"),
        Some("one two ")
    );
    assert_eq!(
        record.data_fields.last(),
        Some(&("bond.prop.Trailing".to_owned(), "one two ".to_owned()))
    );
    assert_eq!(record.topology.atoms.len(), 3);
    assert_eq!(record.topology.bonds.len(), 2);
    assert!(record.properties.sdf_property_lists().is_empty());
    assert!(
        record
            .topology
            .bonds
            .iter()
            .all(|item| item.prop("Trailing").is_none())
    );
}

#[test]
fn sdf_props_unrecognized_names_remain_raw_and_processing_can_be_disabled() {
    let fields = concat!(
        ">  <atom.prop.>\na b c\n\n",
        ">  <Atom.prop.Case>\na b c\n\n",
        ">  <unrelated>\nvalue\n\n",
        ">  <atom.prop.Label>\nfirst second third\n\n",
    );
    let input = concrete_record(fields);
    let record = read_sdf_record_detached(&input).expect("recognized and raw fields");
    assert_eq!(record.properties.sdf_property_lists().len(), 1);
    assert_eq!(record.properties.prop("atom.prop."), Some("a b c"));
    assert_eq!(record.properties.prop("Atom.prop.Case"), Some("a b c"));
    assert_eq!(record.properties.prop("unrelated"), Some("value"));

    let raw = read_sdf_record_detached_with_params(
        &input,
        SdfDataReadParams {
            process_property_lists: false,
            ..SdfDataReadParams::default()
        },
    )
    .expect("property-list processing disabled");
    assert!(raw.properties.sdf_property_lists().is_empty());
    assert!(
        raw.topology
            .atoms
            .iter()
            .all(|atom| atom.prop("Label").is_none())
    );
    assert_eq!(raw.properties.sdf_data_fields().len(), 4);
    assert_eq!(
        raw.properties.prop("atom.prop.Label"),
        Some("first second third")
    );
}

#[test]
fn sdf_props_repeated_lists_keep_encounter_order_and_replace_item_properties() {
    let input = concrete_record(concat!(
        ">  <atom.prop.Label>\none two three\n\n",
        ">  <atom.prop.Label>\nun deux trois\n\n",
    ));
    let record = read_sdf_record_detached(&input).expect("repeated property lists");

    assert_eq!(record.properties.sdf_data_fields().len(), 2);
    assert_eq!(record.properties.sdf_property_lists().len(), 2);
    assert_list(
        &record.properties,
        0,
        SdfPropertyListTarget::Atom,
        "Label",
        &[
            Some(PropertyValue::String("one".to_owned())),
            Some(PropertyValue::String("two".to_owned())),
            Some(PropertyValue::String("three".to_owned())),
        ],
    );
    assert_list(
        &record.properties,
        1,
        SdfPropertyListTarget::Atom,
        "Label",
        &[
            Some(PropertyValue::String("un".to_owned())),
            Some(PropertyValue::String("deux".to_owned())),
            Some(PropertyValue::String("trois".to_owned())),
        ],
    );
    assert_typed_property(
        record.topology.atoms[0].prop("Label"),
        PropertyValue::String("un".to_owned()),
        "un",
    );
    assert_typed_property(
        record.topology.atoms[2].prop("Label"),
        PropertyValue::String("trois".to_owned()),
        "trois",
    );
    assert_eq!(
        record.properties.prop("atom.prop.Label"),
        Some("un deux trois")
    );
}
