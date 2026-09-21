use cosmolkit_io::{
    MolBlockRecord, SdfDataReadParams, SdfReadError, read_sdf_graph_record_detached,
    read_sdf_record_detached, read_sdf_record_detached_with_params,
};
use cosmolkit_model::{MoleculeProperties, SdfPropertyListTarget};

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
    values: &[Option<&str>],
) {
    let list = &properties.sdf_property_lists()[index];
    assert_eq!(list.target(), target);
    assert_eq!(list.name(), name);
    assert_eq!(
        list.values(),
        values
            .iter()
            .map(|value| value.map(str::to_owned))
            .collect::<Vec<_>>()
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

    assert_eq!(record.topology.atoms[0].prop("Text"), Some("alpha"));
    assert_eq!(record.topology.atoms[1].prop("Text"), None);
    assert_eq!(record.topology.atoms[2].prop("Text"), Some("omega"));
    assert_eq!(record.topology.atoms[0].prop("Integer"), Some("+7"));
    assert_eq!(record.topology.atoms[1].prop("Integer"), None);
    assert_eq!(record.topology.atoms[2].prop("Integer"), Some("-3"));
    assert_eq!(record.topology.atoms[0].prop("Real"), Some("1e2"));
    assert_eq!(record.topology.atoms[1].prop("Real"), None);
    assert_eq!(record.topology.atoms[2].prop("Real"), Some("-5e-1"));
    assert_eq!(record.topology.atoms[0].prop("Flag"), Some("true"));
    assert_eq!(record.topology.atoms[1].prop("Flag"), Some("false"));
    assert_eq!(record.topology.atoms[2].prop("Flag"), Some("true"));

    assert_eq!(record.topology.bonds[0].prop("Text"), Some("left"));
    assert_eq!(record.topology.bonds[1].prop("Text"), Some("right"));
    assert_eq!(record.topology.bonds[0].prop("Integer"), Some("-9"));
    assert_eq!(record.topology.bonds[1].prop("Integer"), Some("+11"));
    assert_eq!(record.topology.bonds[0].prop("Real"), Some("6.25e-1"));
    assert_eq!(record.topology.bonds[1].prop("Real"), None);
    assert_eq!(record.topology.bonds[0].prop("Flag"), Some("false"));
    assert_eq!(record.topology.bonds[1].prop("Flag"), Some("true"));

    let lists = record.properties.sdf_property_lists();
    assert_eq!(lists.len(), 8);
    assert_list(
        &record.properties,
        0,
        SdfPropertyListTarget::Atom,
        "Text",
        &[Some("alpha"), None, Some("omega")],
    );
    assert_list(
        &record.properties,
        1,
        SdfPropertyListTarget::Atom,
        "Integer",
        &[Some("+7"), None, Some("-3")],
    );
    assert_list(
        &record.properties,
        6,
        SdfPropertyListTarget::Bond,
        "Real",
        &[Some("6.25e-1"), None],
    );
    assert_eq!(
        record.properties.prop("atom.iprop.Integer"),
        Some("+7 invalid -3")
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
        &[Some("first"), None, Some("third")],
    );
    assert_list(
        &record.properties,
        1,
        SdfPropertyListTarget::Atom,
        "EmptyMarker",
        &[Some("one"), Some("two"), Some("three")],
    );
    assert_eq!(record.topology.atoms[2].prop("Custom"), Some("third"));
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
        ">  <atom.dprop.Weight>\n2.5 invalid\n\n",
        ">  <bond.bprop.Selected>\n1\n\n$$$$\n",
    );
    let record = read_sdf_graph_record_detached(input).expect("query property lists");
    let MolBlockRecord::Query(query) = record.mol_block else {
        panic!("property expansion must not lower the query graph");
    };

    assert_eq!(query.query.atoms()[0].prop("Weight"), Some("2.5"));
    assert_eq!(query.query.atoms()[1].prop("Weight"), None);
    assert_eq!(query.query.bonds()[0].bond().prop("Selected"), Some("true"));
    assert_list(
        &query.properties,
        0,
        SdfPropertyListTarget::Atom,
        "Weight",
        &[Some("2.5"), None],
    );
    assert_list(
        &query.properties,
        1,
        SdfPropertyListTarget::Bond,
        "Selected",
        &[Some("true")],
    );
}

#[test]
fn sdf_props_strict_count_mismatches_are_structured_for_atom_and_bond_targets() {
    let atom_short =
        concrete_record(">  <atom.prop.ValidBefore>\na b c\n\n>  <atom.iprop.Short>\n1 2\n\n");
    assert_eq!(
        read_sdf_record_detached(&atom_short),
        Err(SdfReadError::PropertyListCount {
            target: "atom",
            name: "atom.iprop.Short".to_owned(),
            actual: 2,
            expected: 3,
        })
    );

    let bond_long = concrete_record(">  <bond.prop.Long>\na b c\n\n");
    assert_eq!(
        read_sdf_record_detached(&bond_long),
        Err(SdfReadError::PropertyListCount {
            target: "bond",
            name: "bond.prop.Long".to_owned(),
            actual: 3,
            expected: 2,
        })
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
    assert_eq!(
        read_sdf_record_detached(&empty_value),
        Err(SdfReadError::PropertyListCount {
            target: "atom",
            name: "atom.prop.Empty".to_owned(),
            actual: 1,
            expected: 0,
        })
    );

    // boost::split(token_compress_on) preserves boundary empty tokens. They
    // therefore participate in the source count check instead of being trim.
    let leading = concrete_record(">  <atom.prop.Leading>\n one two three\n\n");
    assert_eq!(
        read_sdf_record_detached(&leading),
        Err(SdfReadError::PropertyListCount {
            target: "atom",
            name: "atom.prop.Leading".to_owned(),
            actual: 4,
            expected: 3,
        })
    );
    let trailing = concrete_record(">  <bond.prop.Trailing>\none two \n\n");
    assert_eq!(
        read_sdf_record_detached(&trailing),
        Err(SdfReadError::PropertyListCount {
            target: "bond",
            name: "bond.prop.Trailing".to_owned(),
            actual: 3,
            expected: 2,
        })
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
        &[Some("one"), Some("two"), Some("three")],
    );
    assert_list(
        &record.properties,
        1,
        SdfPropertyListTarget::Atom,
        "Label",
        &[Some("un"), Some("deux"), Some("trois")],
    );
    assert_eq!(record.topology.atoms[0].prop("Label"), Some("un"));
    assert_eq!(record.topology.atoms[2].prop("Label"), Some("trois"));
    assert_eq!(
        record.properties.prop("atom.prop.Label"),
        Some("un deux trois")
    );
}
