//! V3000 read behavior regressions.
//!
//! Named behavior slices map to pinned RDKit 2026.03.1
//! `MolFileParser.cpp` / `MolSGroupParsing.cpp` and are exercised through the
//! public detached reader; private helpers stay in their owning module.

use cosmolkit_io::{MolBlockRecord, SdfReadError, read_mol_block_detached};
use cosmolkit_model::{
    AtomId, AtomQueryPredicate, BondId, BondQueryPredicate, CoordinateDimension, QueryNode,
    SGroupBondRole, SGroupBracketStyle, SGroupConnection, StereoGroupKind, query_substance_groups,
};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo};

#[test]
fn review_utf8_attchord_byte_truncation_is_a_parse_error() {
    // RDKit substr(1, size - 2) cuts a byte. The Rust label must remain UTF-8;
    // invalid byte truncation is a structured error, not a panic or lossy label.
    let block = v3000_block(&["M  V30 1 C 0 0 0 0 ATTCHORD=(2 1 é"], 1);
    assert!(matches!(
        read_mol_block_detached(&block),
        Err(SdfReadError::Parse(_))
    ));
    let block = v3000_block(&["M  V30 1 C 0 0 0 0 ATTCHORD=(2 1 é)"], 1);
    let MolBlockRecord::Concrete { topology, .. } = read_mol_block_detached(&block).unwrap() else {
        panic!("concrete record expected");
    };
    assert_eq!(
        topology.atoms[0]
            .template_attachment_order()
            .unwrap()
            .entries()[0]
            .label(),
        "é"
    );
}

#[test]
fn review_utf8_unknown_blocks_and_collections_are_skipped() {
    // Pinned ParseV3000CTAB skips unknown blocks; parseEnhancedStereo skips
    // nonmatching collection lines. Neither operation slices Unicode strings.
    for extra in [
        "M  V30 BEGIN😀ABCDEFG\nM  V30 END UNKNOWN\n",
        "M  V30 BEGIN COLLECTION\nM  V30 MDLV30/STE😀 ATOMS=(1 1)\nM  V30 END COLLECTION\n",
    ] {
        let block = v3000_block(&["M  V30 1 C 0 0 0 0"], 1)
            .replace("M  V30 END CTAB", &format!("{extra}M  V30 END CTAB"));
        let MolBlockRecord::Concrete { topology, .. } = read_mol_block_detached(&block).unwrap()
        else {
            panic!("concrete record expected");
        };
        assert_eq!(topology.atoms.len(), 1);
        assert!(topology.stereo_groups.is_empty());
    }
}

fn v3000_block(atom_lines: &[&str], atom_count: usize) -> String {
    let mut block = String::from(concat!(
        "\n",
        "  COSMolKit\n",
        "\n",
        "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
        "M  V30 BEGIN CTAB\n",
    ));
    block.push_str(&format!("M  V30 COUNTS {atom_count} 0 0 0 0\n"));
    block.push_str("M  V30 BEGIN ATOM\n");
    for line in atom_lines {
        block.push_str(line);
        block.push('\n');
    }
    block.push_str("M  V30 END ATOM\nM  V30 END CTAB\nM  END\n");
    block
}

fn v3k_sg_array_record(labels: &str, strict_parsing: bool) -> Result<MolBlockRecord, SdfReadError> {
    let block = v3000_block(&["M  V30 1 C 0 0 0 0", "M  V30 2 N 1 0 0 0"], 2)
        .replace("M  V30 COUNTS 2 0 0 0 0", "M  V30 COUNTS 2 0 1 0 0")
        .replace(
            "M  V30 END CTAB",
            &format!(
                "M  V30 BEGIN SGROUP\nM  V30 1 SUP 0 {labels}\nM  V30 END SGROUP\nM  V30 END CTAB"
            ),
        );
    cosmolkit_io::read_mol_block_detached_with_params(
        &block,
        cosmolkit_io::MolBlockReadParams {
            strict_parsing,
            ..cosmolkit_io::MolBlockReadParams::default()
        },
    )
}

fn v3k_sg_string_record(
    labels: &str,
    strict_parsing: bool,
) -> Result<MolBlockRecord, SdfReadError> {
    v3k_sg_array_record(labels, strict_parsing)
}

fn v3k_sg_default_record(
    defaults: Option<&str>,
    labels: &str,
    strict_parsing: bool,
) -> Result<MolBlockRecord, SdfReadError> {
    let default_line = defaults
        .map(|value| format!("M  V30 DEFAULT {value}\n"))
        .unwrap_or_default();
    let block = format!(
        concat!(
            "\n",
            "  COSMolKit\n",
            "\n",
            "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 1 0 1 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 1 C 0 0 0 0\n",
            "M  V30 END ATOM\n",
            "M  V30 BEGIN BOND\n",
            "M  V30 END BOND\n",
            "M  V30 BEGIN SGROUP\n",
            "{default_line}",
            "M  V30 1 DAT 0 {labels}\n",
            "M  V30 END SGROUP\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        ),
        default_line = default_line,
        labels = labels,
    );
    cosmolkit_io::read_mol_block_detached_with_params(
        &block,
        cosmolkit_io::MolBlockReadParams {
            strict_parsing,
            ..cosmolkit_io::MolBlockReadParams::default()
        },
    )
}

fn v3k_sg_identity_record(
    rows: &[&str],
    declared_count: usize,
    strict_parsing: bool,
) -> Result<MolBlockRecord, SdfReadError> {
    let rows = rows
        .iter()
        .map(|row| format!("M  V30 {row}\n"))
        .collect::<String>();
    let block = format!(
        concat!(
            "\n",
            "  COSMolKit\n",
            "\n",
            "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 1 0 {declared_count} 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 1 C 0 0 0 0\n",
            "M  V30 END ATOM\n",
            "M  V30 BEGIN BOND\n",
            "M  V30 END BOND\n",
            "M  V30 BEGIN SGROUP\n",
            "{rows}",
            "M  V30 END SGROUP\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        ),
        declared_count = declared_count,
        rows = rows,
    );
    cosmolkit_io::read_mol_block_detached_with_params(
        &block,
        cosmolkit_io::MolBlockReadParams {
            strict_parsing,
            ..cosmolkit_io::MolBlockReadParams::default()
        },
    )
}

fn v3k_sg_member_record(
    labels: &str,
    strict_parsing: bool,
) -> Result<MolBlockRecord, SdfReadError> {
    let block = format!(
        concat!(
            "\n",
            "  COSMolKit\n",
            "\n",
            "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 3 2 1 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 10 C 0 0 0 0\n",
            "M  V30 20 N 1 0 0 0\n",
            "M  V30 30 O 2 0 0 0\n",
            "M  V30 END ATOM\n",
            "M  V30 BEGIN BOND\n",
            "M  V30 50 1 10 20\n",
            "M  V30 60 1 20 30\n",
            "M  V30 END BOND\n",
            "M  V30 BEGIN SGROUP\n",
            "M  V30 1 SUP 0 {labels}\n",
            "M  V30 END SGROUP\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        ),
        labels = labels,
    );
    cosmolkit_io::read_mol_block_detached_with_params(
        &block,
        cosmolkit_io::MolBlockReadParams {
            strict_parsing,
            ..cosmolkit_io::MolBlockReadParams::default()
        },
    )
}

fn v3k_sg_cstate_record(
    kind: &str,
    labels: &str,
    strict_parsing: bool,
) -> Result<MolBlockRecord, SdfReadError> {
    let block = format!(
        concat!(
            "\n",
            "  COSMolKit\n",
            "\n",
            "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 3 2 1 0 0\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 10 C 0 0 0 0\n",
            "M  V30 20 N 1 0 0 0\n",
            "M  V30 30 O 2 0 0 0\n",
            "M  V30 END ATOM\n",
            "M  V30 BEGIN BOND\n",
            "M  V30 50 1 10 20\n",
            "M  V30 60 1 20 30\n",
            "M  V30 END BOND\n",
            "M  V30 BEGIN SGROUP\n",
            "M  V30 1 {kind} 0 {labels}\n",
            "M  V30 END SGROUP\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        ),
        kind = kind,
        labels = labels,
    );
    cosmolkit_io::read_mol_block_detached_with_params(
        &block,
        cosmolkit_io::MolBlockReadParams {
            strict_parsing,
            ..cosmolkit_io::MolBlockReadParams::default()
        },
    )
}

fn v3k_sg_data_record(labels: &str, strict_parsing: bool) -> Result<MolBlockRecord, SdfReadError> {
    v3k_sg_cstate_record("DAT", labels, strict_parsing)
}

fn v3k_sg_sap_record(labels: &str, strict_parsing: bool) -> Result<MolBlockRecord, SdfReadError> {
    v3k_sg_member_record(labels, strict_parsing)
}

#[test]
fn v3k_sg_members_resolve_nonsequential_bookmarks_and_preserve_multiplicity() {
    // Pinned ParseV3000ParseLabel resolves molecule bookmarks through the
    // owning molecule and each SubstanceGroup add helper appends without
    // sorting or deduplicating.
    for strict_parsing in [true, false] {
        let MolBlockRecord::Concrete { topology, .. } = v3k_sg_member_record(
            "ATOMS=(3 20 10 20) PATOMS=(2 20 20) CBONDS=(2 50 50) XBONDS=(1 60)",
            strict_parsing,
        )
        .expect("ordered duplicate SGroup membership") else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        let group = &topology.substance_groups[0];
        assert_eq!(
            group
                .atoms()
                .iter()
                .map(|id| id.index())
                .collect::<Vec<_>>(),
            vec![1, 0, 1]
        );
        assert_eq!(
            group
                .parent_atoms()
                .iter()
                .map(|id| id.index())
                .collect::<Vec<_>>(),
            vec![1, 1]
        );
        assert_eq!(
            group
                .bonds()
                .iter()
                .map(|id| id.index())
                .collect::<Vec<_>>(),
            vec![0, 0, 1]
        );
        assert_eq!(group.bond_role(group.bonds()[0]), SGroupBondRole::Contained);
        assert_eq!(group.bond_role(group.bonds()[2]), SGroupBondRole::Crossing);
    }
}

#[test]
fn v3k_sg_members_empty_arrays_preserve_an_empty_typed_group() {
    // Source arrays with a zero declared count call no add helper and leave a
    // valid group with empty typed membership vectors.
    for strict_parsing in [true, false] {
        let MolBlockRecord::Concrete { topology, .. } =
            v3k_sg_member_record("ATOMS=(0) PATOMS=(0) CBONDS=(0) XBONDS=(0)", strict_parsing)
                .expect("empty SGroup memberships")
        else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        let group = &topology.substance_groups[0];
        assert!(group.atoms().is_empty());
        assert!(group.parent_atoms().is_empty());
        assert!(group.bonds().is_empty());
    }
}

#[test]
fn v3k_sg_members_missing_bookmarks_fail_strict_and_drop_the_whole_group_non_strict() {
    // Each source add-with-bookmark helper requires a unique existing molecule
    // bookmark. The enclosing parser invalidates the complete SGroup on a
    // non-strict failure instead of retaining earlier partial membership.
    for labels in [
        "ATOMS=(1 99)",
        "ATOMS=(1 10) PATOMS=(1 99)",
        "ATOMS=(2 10 20) CBONDS=(1 99)",
        "ATOMS=(2 10 20) XBONDS=(1 99)",
    ] {
        assert!(matches!(
            v3k_sg_member_record(labels, true),
            Err(SdfReadError::Parse(_))
        ));
        let MolBlockRecord::Concrete { topology, .. } =
            v3k_sg_member_record(labels, false).expect("invalid SGroup is discarded")
        else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        assert!(topology.substance_groups.is_empty());
    }
}

#[test]
fn v3k_sg_members_parent_atoms_must_already_be_atom_members() {
    // SubstanceGroup::addParentAtomWithBookmark performs a linear membership
    // check at the point of the PATOMS label, so both an unrelated atom and a
    // PATOMS-before-ATOMS ordering are invalid.
    for labels in ["ATOMS=(1 10) PATOMS=(1 20)", "PATOMS=(1 10) ATOMS=(1 10)"] {
        assert!(matches!(
            v3k_sg_member_record(labels, true),
            Err(SdfReadError::Parse(_))
        ));
        let MolBlockRecord::Concrete { topology, .. } =
            v3k_sg_member_record(labels, false).expect("invalid SGroup is discarded")
        else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        assert!(topology.substance_groups.is_empty());
    }
}

#[test]
fn v3k_sg_bond_refs_use_one_based_rows_not_bookmarks_and_preserve_order() {
    // Pinned ParseV3000ParseLabel parses XBHEAD/XBCORR as unsigned arrays and
    // subtracts one. The values are bond-table positions, so bookmark 50 does
    // not refer to the first bond here; source values 1 and 2 do.
    for strict_parsing in [true, false] {
        let MolBlockRecord::Concrete { topology, .. } = v3k_sg_member_record(
            "ATOMS=(3 10 20 30) XBONDS=(2 50 60) XBHEAD=(2 2 1) XBCORR=(2 1 1)",
            strict_parsing,
        )
        .expect("typed crossing-bond row references") else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        let group = &topology.substance_groups[0];
        assert_eq!(
            group
                .head_crossing_bonds()
                .iter()
                .map(|bond| bond.index())
                .collect::<Vec<_>>(),
            vec![1, 0]
        );
        assert_eq!(
            group
                .crossing_bond_correspondence()
                .iter()
                .map(|bond| bond.index())
                .collect::<Vec<_>>(),
            vec![0, 0]
        );
        assert!(!group.props().contains_key("XBHEAD"));
        assert!(!group.props().contains_key("XBCORR"));
    }
}

#[test]
fn v3k_sg_bond_refs_zero_and_out_of_range_fail_the_canonical_reference_boundary() {
    // RDKit's source-level unsigned subtraction retains 0 as UINT_MAX and can
    // retain other out-of-range property values. COSMolKit's unique canonical
    // SGroup representation uses BondId and therefore rejects those invalid
    // persistent references instead of flattening them into raw properties.
    for labels in [
        "ATOMS=(3 10 20 30) XBHEAD=(1 0)",
        "ATOMS=(3 10 20 30) XBHEAD=(1 3)",
        "ATOMS=(3 10 20 30) XBCORR=(1 0)",
        "ATOMS=(3 10 20 30) XBCORR=(1 3)",
    ] {
        assert!(matches!(
            v3k_sg_member_record(labels, true),
            Err(SdfReadError::Parse(_))
        ));
        let MolBlockRecord::Concrete { topology, .. } =
            v3k_sg_member_record(labels, false).expect("invalid SGroup is discarded")
        else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        assert!(topology.substance_groups.is_empty());
    }
}

#[test]
fn v3k_sg_bond_refs_reject_declared_xbcorr_count_above_the_bond_table() {
    // Pinned ParseV3000ParseLabel passes getNumBonds() as ParseV3000Array's
    // maximum. The declared count is checked before elements are read, so
    // duplicate references do not exempt a count above the two-bond table.
    assert!(matches!(
        v3k_sg_member_record("ATOMS=(3 10 20 30) XBCORR=(3 1 1 2)", true),
        Err(SdfReadError::Parse(message)) if message == "invalid count value"
    ));
    let MolBlockRecord::Concrete { topology, .. } =
        v3k_sg_member_record("ATOMS=(3 10 20 30) XBCORR=(3 1 1 2)", false)
            .expect("non-strict over-limit array warns and returns empty")
    else {
        panic!("ordinary SGroup carrier must remain concrete");
    };
    assert_eq!(topology.substance_groups.len(), 1);
    assert!(
        topology.substance_groups[0]
            .crossing_bond_correspondence()
            .is_empty()
    );
}

#[test]
fn v3k_sg_bond_refs_empty_arrays_remain_empty_typed_vectors() {
    for strict_parsing in [true, false] {
        let MolBlockRecord::Concrete { topology, .. } =
            v3k_sg_member_record("ATOMS=(3 10 20 30) XBHEAD=(0) XBCORR=(0)", strict_parsing)
                .expect("empty typed bond-reference arrays")
        else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        let group = &topology.substance_groups[0];
        assert!(group.head_crossing_bonds().is_empty());
        assert!(group.crossing_bond_correspondence().is_empty());
    }
}

#[test]
fn v3k_sg_data_strict_fielddata_uses_the_source_200_byte_boundary() {
    for length in [199_usize, 200, 201] {
        let value = "a".repeat(length);
        let labels = format!("ATOMS=(1 20) FIELDDATA=\"{value}\"");
        let MolBlockRecord::Concrete { topology, .. } =
            v3k_sg_data_record(&labels, true).expect("strict ASCII DAT field")
        else {
            panic!("ordinary DAT carrier must remain concrete");
        };
        assert_eq!(
            topology.substance_groups[0].data().unwrap().values[0],
            "a".repeat(length.min(200)),
            "source strict boundary at {length} bytes"
        );

        let MolBlockRecord::Concrete { topology, .. } =
            v3k_sg_data_record(&labels, false).expect("non-strict ASCII DAT field")
        else {
            panic!("ordinary DAT carrier must remain concrete");
        };
        assert_eq!(
            topology.substance_groups[0].data().unwrap().values[0],
            value,
            "non-strict source path preserves {length} bytes"
        );
    }
}

#[test]
fn v3k_sg_data_utf8_cut_is_an_explicit_rust_text_boundary() {
    let exact = format!("{}é", "a".repeat(198));
    let labels = format!("ATOMS=(1 20) FIELDDATA=\"{exact}\"");
    let MolBlockRecord::Concrete { topology, .. } =
        v3k_sg_data_record(&labels, true).expect("200-byte UTF-8 value")
    else {
        panic!("ordinary DAT carrier must remain concrete");
    };
    assert_eq!(topology.substance_groups[0].data().unwrap().values, [exact]);

    let split = format!("{}é", "a".repeat(199));
    let labels = format!("ATOMS=(1 20) FIELDDATA=\"{split}\"");
    assert!(matches!(
        v3k_sg_data_record(&labels, true),
        Err(SdfReadError::Parse(message))
            if message.contains("FIELDDATA 200-byte truncation splits UTF-8")
    ));
    let MolBlockRecord::Concrete { topology, .. } =
        v3k_sg_data_record(&labels, false).expect("non-strict UTF-8 value")
    else {
        panic!("ordinary DAT carrier must remain concrete");
    };
    assert_eq!(topology.substance_groups[0].data().unwrap().values, [split]);
}

#[test]
fn v3k_sg_data_preserves_typed_metadata_empty_values_escapes_and_order() {
    let labels = concat!(
        "ATOMS=(1 20) FIELDNAME=\"FIELD NAME\" FIELDTYPE=TYPE ",
        "FIELDINFO=\"UNIT INFO\" FIELDDISP=\"DISPLAY SPEC\" ",
        "QUERYTYPE=PQ QUERYOP== FIELDDATA=\"\" FIELDDATA=\"a\"\"b\""
    );
    for strict_parsing in [true, false] {
        let MolBlockRecord::Concrete { topology, .. } =
            v3k_sg_data_record(labels, strict_parsing).expect("typed DAT metadata")
        else {
            panic!("ordinary DAT carrier must remain concrete");
        };
        let data = topology.substance_groups[0].data().unwrap();
        assert_eq!(data.field_name.as_deref(), Some("FIELD NAME"));
        assert_eq!(data.field_type.as_deref(), Some("TYPE"));
        assert_eq!(data.field_info.as_deref(), Some("UNIT INFO"));
        assert_eq!(data.field_display.as_deref(), Some("DISPLAY SPEC"));
        assert_eq!(data.query_type.as_deref(), Some("PQ"));
        assert_eq!(data.query_op.as_deref(), Some("="));
        assert_eq!(data.values, ["", "a\"b"]);
    }
}

#[test]
fn v3k_sg_data_special_fields_remain_typed_until_mol_postprocessing() {
    // Pinned ParseV3000ParseLabel stores this DAT row. The later
    // processSGroups pass interprets MRV_IMPLICIT_H and removes the group;
    // the detached reader must not perform that later mol-post behavior.
    let MolBlockRecord::Concrete { topology, .. } =
        v3k_sg_data_record("ATOMS=(1 20) FIELDNAME=MRV_IMPLICIT_H FIELDDATA=5", true)
            .expect("detached special DAT state")
    else {
        panic!("ordinary DAT carrier must remain concrete");
    };
    assert_eq!(topology.substance_groups.len(), 1);
    let data = topology.substance_groups[0].data().unwrap();
    assert_eq!(data.field_name.as_deref(), Some("MRV_IMPLICIT_H"));
    assert_eq!(data.values, ["5"]);
    assert_eq!(topology.atoms[1].explicit_hydrogens(), 0);
    assert_eq!(topology.atoms[1].prop("_ZBO_H"), None);
}

#[test]
fn v3k_sg_text_labels_preserve_canonical_known_values_without_raw_duplicates() {
    // Pinned ParseV3000ParseLabel sends every remaining label through
    // ParseV3000StringPropLabel and setProp. COSMolKit projects the known
    // canonical values into their unique typed fields rather than keeping a
    // second raw representation.
    for strict_parsing in [true, false] {
        let MolBlockRecord::Concrete { topology, .. } = v3k_sg_string_record(
            "LABEL=first LABEL=\"final label\" ESTATE=E BRKTYP=PAREN",
            strict_parsing,
        )
        .expect("canonical SGroup string labels") else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        let group = &topology.substance_groups[0];
        assert_eq!(group.label(), Some("final label"));
        assert_eq!(group.expansion_state(), Some("E"));
        assert_eq!(
            group.bracket_style(),
            Some(&SGroupBracketStyle::Parenthesis)
        );
        for key in ["LABEL", "ESTATE", "BRKTYP"] {
            assert!(!group.props().contains_key(key), "duplicate raw {key}");
        }
    }
}

#[test]
fn v3k_sg_text_labels_preserve_unvalidated_and_unknown_raw_metadata() {
    // The pinned source explicitly does not validate or interpret
    // NATREPLACE. Unknown string properties follow the same source tail.
    let MolBlockRecord::Concrete { topology, .. } = v3k_sg_string_record(
        "NATREPLACE=\"not/a/known/template\" VENDOR=\"a\"\"b\" EMPTY=",
        true,
    )
    .expect("raw SGroup metadata") else {
        panic!("ordinary SGroup carrier must remain concrete");
    };
    let props = topology.substance_groups[0].props();
    assert_eq!(
        props.get("NATREPLACE").map(String::as_str),
        Some("not/a/known/template")
    );
    assert_eq!(props.get("VENDOR").map(String::as_str), Some("a\"b"));
    assert_eq!(props.get("EMPTY").map(String::as_str), Some(""));
}

#[test]
fn v3k_sg_text_labels_preserve_source_bracket_style_values() {
    for (source, expected) in [
        ("BRACKET", SGroupBracketStyle::Bracket),
        ("PAREN", SGroupBracketStyle::Parenthesis),
        ("", SGroupBracketStyle::None),
        (
            "vendor-style",
            SGroupBracketStyle::Unknown("vendor-style".to_owned()),
        ),
    ] {
        let labels = format!("BRKTYP={source}");
        let MolBlockRecord::Concrete { topology, .. } =
            v3k_sg_string_record(&labels, true).expect("source bracket style")
        else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        let group = &topology.substance_groups[0];
        assert_eq!(group.bracket_style(), Some(&expected), "{source:?}");
        assert!(!group.props().contains_key("BRKTYP"));
    }
}

#[test]
fn v3k_sg_text_labels_never_flatten_typed_bond_references() {
    let MolBlockRecord::Concrete { topology, .. } = v3k_sg_member_record(
        "ATOMS=(3 10 20 30) XBHEAD=(2 2 1) XBCORR=(2 1 1) NATREPLACE=AA/X",
        true,
    )
    .expect("typed crossing-bond references") else {
        panic!("ordinary SGroup carrier must remain concrete");
    };
    let group = &topology.substance_groups[0];
    assert_eq!(
        group
            .head_crossing_bonds()
            .iter()
            .map(|bond| bond.index())
            .collect::<Vec<_>>(),
        [1, 0]
    );
    assert_eq!(
        group
            .crossing_bond_correspondence()
            .iter()
            .map(|bond| bond.index())
            .collect::<Vec<_>>(),
        [0, 0]
    );
    assert!(!group.props().contains_key("XBHEAD"));
    assert!(!group.props().contains_key("XBCORR"));
    assert_eq!(
        group.props().get("NATREPLACE").map(String::as_str),
        Some("AA/X")
    );

    assert!(matches!(
        v3k_sg_member_record("XBHEAD=(1 0)", true),
        Err(SdfReadError::Parse(message)) if message.contains("bond-row position 0")
    ));
}

#[test]
fn v3k_sg_defaults_apply_only_when_explicit_labels_are_absent() {
    // Pinned ParseV3000SGroupsBlock parses the saved DEFAULT string after the
    // row and calls ParseV3000ParseLabel only for labels absent from the row's
    // `parsedLabels` set.
    for strict_parsing in [true, false] {
        let MolBlockRecord::Concrete { topology, .. } = v3k_sg_default_record(
            Some("CLASS=AA LABEL=default ESTATE=E BRKTYP=PAREN"),
            "ATOMS=(1 1)",
            strict_parsing,
        )
        .expect("default-only canonical values") else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        let group = &topology.substance_groups[0];
        assert_eq!(group.class(), Some("AA"));
        assert_eq!(group.label(), Some("default"));
        assert_eq!(group.expansion_state(), Some("E"));
        assert_eq!(
            group.bracket_style(),
            Some(&SGroupBracketStyle::Parenthesis)
        );

        let MolBlockRecord::Concrete { topology, .. } =
            v3k_sg_default_record(None, "LABEL=explicit", strict_parsing)
                .expect("absent DEFAULT line")
        else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        assert_eq!(topology.substance_groups[0].label(), Some("explicit"));
    }
}

#[test]
fn v3k_sg_defaults_preserve_explicit_precedence_and_default_encounter_order() {
    // `parsedLabels` contains explicit row labels only. Repeated defaults that
    // are not explicit are therefore all parsed, and the later source property
    // assignment wins; every explicit occurrence suppresses every default.
    for strict_parsing in [true, false] {
        let MolBlockRecord::Concrete { topology, .. } = v3k_sg_default_record(
            Some("LABEL=default VENDOR=first VENDOR=second ESTATE=E"),
            "LABEL=early LABEL=explicit",
            strict_parsing,
        )
        .expect("explicit/default precedence") else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        let group = &topology.substance_groups[0];
        assert_eq!(group.label(), Some("explicit"));
        assert_eq!(group.expansion_state(), Some("E"));
        assert_eq!(
            group.props().get("VENDOR").map(String::as_str),
            Some("second")
        );
    }
}

#[test]
fn v3k_sg_defaults_skip_overridden_quote_and_parenthesis_with_source_cursor_rules() {
    // The source skip branch is not ParseV3000StringPropLabel: a quoted value
    // stops at its first closing quote, while a parenthesized value consumes
    // one additional byte after `)`. With a following label that byte is its
    // separator, so the next outer read sees `E` and invalidates the group.
    for strict_parsing in [true, false] {
        let MolBlockRecord::Concrete { topology, .. } = v3k_sg_default_record(
            Some("LABEL=\"ignored value\" ESTATE=E"),
            "LABEL=explicit",
            strict_parsing,
        )
        .expect("quoted overridden default") else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        let group = &topology.substance_groups[0];
        assert_eq!(group.label(), Some("explicit"));
        assert_eq!(group.expansion_state(), Some("E"));
    }

    assert!(matches!(
        v3k_sg_default_record(Some("LABEL=(2 a b) ESTATE=E"), "LABEL=explicit", true,),
        Err(SdfReadError::Parse(_))
    ));
    let MolBlockRecord::Concrete { topology, .. } =
        v3k_sg_default_record(Some("LABEL=(2 a b) ESTATE=E"), "LABEL=explicit", false)
            .expect("non-strict invalid SGroup is discarded")
    else {
        panic!("ordinary SGroup carrier must remain concrete");
    };
    assert!(topology.substance_groups.is_empty());
}

#[test]
fn v3k_sg_defaults_empty_overridden_value_invalidates_the_complete_group() {
    // In the pinned source, a literal space immediately after an overridden
    // DEFAULT `=` is an explicit parse error rather than an empty string.
    assert!(matches!(
        v3k_sg_default_record(Some("LABEL= ESTATE=E"), "LABEL=explicit", true),
        Err(SdfReadError::Parse(message))
            if message.contains("unexpected whitespace at DEFAULT label LABEL")
    ));
    let MolBlockRecord::Concrete { topology, .. } =
        v3k_sg_default_record(Some("LABEL= ESTATE=E"), "LABEL=explicit", false)
            .expect("non-strict invalid SGroup is discarded")
    else {
        panic!("ordinary SGroup carrier must remain concrete");
    };
    assert!(topology.substance_groups.is_empty());
}

#[test]
fn v3k_sg_identity_sorts_sequence_keys_and_assigns_compact_canonical_ids() {
    // Pinned ParseV3000SGroupsBlock stores rows in std::map keyed by source
    // sequence ID, then installs them in sorted key order. Source sequence and
    // external IDs remain metadata; canonical IDs are compact output rows.
    for strict_parsing in [true, false] {
        let MolBlockRecord::Concrete { topology, .. } = v3k_sg_identity_record(
            &["20 DAT 200 LABEL=twenty", "3 DAT 30 LABEL=three"],
            2,
            strict_parsing,
        )
        .expect("unsorted source sequence IDs") else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        assert_eq!(topology.substance_groups.len(), 2);
        assert_eq!(topology.substance_groups[0].id().index(), 0);
        assert_eq!(topology.substance_groups[1].id().index(), 1);
        assert_eq!(topology.substance_groups[0].rdkit_sequence_id(), Some(3));
        assert_eq!(topology.substance_groups[1].rdkit_sequence_id(), Some(20));
        assert_eq!(topology.substance_groups[0].external_id(), Some(30));
        assert_eq!(topology.substance_groups[1].external_id(), Some(200));
        assert_eq!(topology.substance_groups[0].label(), Some("three"));
        assert_eq!(topology.substance_groups[1].label(), Some("twenty"));
    }
}

#[test]
fn v3k_sg_identity_duplicate_sequence_is_first_wins_and_count_checked() {
    // std::map::emplace retains the first row. The source then compares map
    // size with the declared count: strict parsing fails, non-strict parsing
    // warns and installs the first row only.
    let rows = ["1 DAT 10 LABEL=first", "1 DAT 20 LABEL=second"];
    assert!(matches!(
        v3k_sg_identity_record(&rows, 2, true),
        Err(SdfReadError::Parse(message)) if message.contains("Found 1 SGroups when 2 were expected")
    ));
    let MolBlockRecord::Concrete { topology, .. } =
        v3k_sg_identity_record(&rows, 2, false).expect("non-strict first-wins map")
    else {
        panic!("ordinary SGroup carrier must remain concrete");
    };
    assert_eq!(topology.substance_groups.len(), 1);
    assert_eq!(topology.substance_groups[0].external_id(), Some(10));
    assert_eq!(topology.substance_groups[0].label(), Some("first"));
}

#[test]
fn v3k_sg_identity_external_lookup_does_not_see_temporary_same_block_rows() {
    // isSubstanceGroupIdFree searches mol.d_sgroups. V3000 rows remain in the
    // temporary map until the whole block is parsed, so equal positive IDs in
    // two distinct source rows are both retained by pinned RDKit.
    for strict_parsing in [true, false] {
        let MolBlockRecord::Concrete { topology, .. } = v3k_sg_identity_record(
            &["1 DAT 7 LABEL=one", "2 DAT 7 LABEL=two"],
            2,
            strict_parsing,
        )
        .expect("same-block duplicate external IDs") else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        assert_eq!(topology.substance_groups.len(), 2);
        assert_eq!(topology.substance_groups[0].external_id(), Some(7));
        assert_eq!(topology.substance_groups[1].external_id(), Some(7));
    }
}

#[test]
fn v3k_sg_identity_formatted_signs_and_malformed_boundaries_are_explicit() {
    // Formatted unsigned extraction accepts both signs; a representable
    // negative magnitude is assigned modulo the u32 destination width.
    let MolBlockRecord::Concrete { topology, .. } =
        v3k_sg_identity_record(&["-1 DAT +7 LABEL=signed"], 1, true)
            .expect("source formatted signs")
    else {
        panic!("ordinary SGroup carrier must remain concrete");
    };
    assert_eq!(
        topology.substance_groups[0].rdkit_sequence_id(),
        Some(u32::MAX)
    );
    assert_eq!(topology.substance_groups[0].external_id(), Some(7));

    // The pinned C++ destinations are uninitialized on no conversion and its
    // overflow failbit prevents later header extraction. COSMolKit rejects
    // those undefined/uninitialized cases instead of fabricating IDs/types.
    for row in [
        "not-an-id DAT 0 LABEL=x",
        "4294967296 DAT 0 LABEL=x",
        "1 DAT not-an-id LABEL=x",
        "1 DAT 4294967296 LABEL=x",
    ] {
        for strict_parsing in [true, false] {
            assert!(
                matches!(
                    v3k_sg_identity_record(&[row], 1, strict_parsing),
                    Err(SdfReadError::Parse(_))
                ),
                "row={row:?}, strict={strict_parsing}"
            );
        }
    }
}

#[test]
fn v3k_sg_identity_unknown_count_and_early_end_have_safe_source_shaped_termination() {
    // A zero declared count is rejected strictly by the CTAB caller. In
    // non-strict mode the source switches to UINT_MAX and terminates when the
    // END row cannot start an SGroup; all valid rows before it survive.
    assert!(matches!(
        v3k_sg_identity_record(&["9 DAT 0 LABEL=nine"], 0, true),
        Err(SdfReadError::Parse(_))
    ));
    let MolBlockRecord::Concrete { topology, .. } =
        v3k_sg_identity_record(&["9 DAT 0 LABEL=nine"], 0, false)
            .expect("unknown-count non-strict block")
    else {
        panic!("ordinary SGroup carrier must remain concrete");
    };
    assert_eq!(topology.substance_groups.len(), 1);
    assert_eq!(topology.substance_groups[0].rdkit_sequence_id(), Some(9));

    // For a known count, an early END is a strict count failure. The source's
    // non-strict path reads uninitialized header destinations and may install
    // an empty artifact; the canonical model deliberately retains only the
    // valid prefix and never persists that undefined source state.
    assert!(matches!(
        v3k_sg_identity_record(&["1 DAT 0 LABEL=one"], 2, true),
        Err(SdfReadError::Parse(message)) if message.contains("Found 1 SGroups when 2 were expected")
    ));
    let MolBlockRecord::Concrete { topology, .. } =
        v3k_sg_identity_record(&["1 DAT 0 LABEL=one"], 2, false)
            .expect("safe non-strict early END")
    else {
        panic!("ordinary SGroup carrier must remain concrete");
    };
    assert_eq!(topology.substance_groups.len(), 1);
    assert_eq!(topology.substance_groups[0].label(), Some("one"));
}

#[test]
fn v3k_sg_identity_type_validation_preserves_strictness_policy() {
    assert!(matches!(
        v3k_sg_identity_record(&["1 VENDOR 0 LABEL=x"], 1, true),
        Err(SdfReadError::Parse(message)) if message.contains("Unsupported SGroup type 'VENDOR'")
    ));
    let MolBlockRecord::Concrete { topology, .. } =
        v3k_sg_identity_record(&["1 VENDOR 0 LABEL=x"], 1, false)
            .expect("non-strict unknown SGroup type")
    else {
        panic!("ordinary SGroup carrier must remain concrete");
    };
    assert_eq!(topology.substance_groups.len(), 1);
    assert_eq!(topology.substance_groups[0].label(), Some("x"));
}

#[test]
fn v3k_sg_parent_resolves_forward_unsorted_sequence_to_compact_typed_id() {
    // PARENT stores a source sequence ID. Rows are sorted by sequence before
    // canonical IDs are assigned, so a forward relation from 20 to 3 resolves
    // to the compact ID of sequence 3 rather than either source integer.
    for strict_parsing in [true, false] {
        let MolBlockRecord::Concrete { topology, .. } = v3k_sg_identity_record(
            &["20 DAT 0 PARENT=3 LABEL=child", "3 DAT 0 LABEL=parent"],
            2,
            strict_parsing,
        )
        .expect("forward unsorted parent") else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        assert_eq!(topology.substance_groups.len(), 2);
        assert_eq!(topology.substance_groups[0].rdkit_sequence_id(), Some(3));
        assert_eq!(topology.substance_groups[0].parent(), None);
        assert_eq!(topology.substance_groups[1].rdkit_sequence_id(), Some(20));
        assert_eq!(
            topology.substance_groups[1].parent(),
            Some(topology.substance_groups[0].id())
        );
        assert!(!topology.substance_groups[1].props().contains_key("PARENT"));
    }
}

#[test]
fn v3k_sg_parent_missing_target_errors_strictly_and_drops_complete_child_non_strictly() {
    let rows = ["1 DAT 0 LABEL=sibling", "2 DAT 0 PARENT=99 LABEL=child"];
    assert!(matches!(
        v3k_sg_identity_record(&rows, 2, true),
        Err(SdfReadError::Parse(message))
            if message.contains("SGroup 2 references missing parent SGroup 99")
    ));
    let MolBlockRecord::Concrete { topology, .. } =
        v3k_sg_identity_record(&rows, 2, false).expect("non-strict missing parent")
    else {
        panic!("ordinary SGroup carrier must remain concrete");
    };
    assert_eq!(topology.substance_groups.len(), 1);
    assert_eq!(topology.substance_groups[0].label(), Some("sibling"));
    assert_eq!(topology.substance_groups[0].id().index(), 0);
}

#[test]
fn v3k_sg_parent_dropped_parent_cascades_without_removing_independent_siblings() {
    // PATOMS invalidates sequence 1 because atom bookmark 2 is absent. The
    // canonical no-dangling boundary then drops its child and grandchild as
    // complete groups in non-strict mode, retaining the independent sibling.
    let rows = [
        "1 DAT 0 ATOMS=(1 1) PATOMS=(1 2) LABEL=invalid-parent",
        "2 DAT 0 PARENT=1 LABEL=child",
        "3 DAT 0 PARENT=2 LABEL=grandchild",
        "4 DAT 0 LABEL=sibling",
    ];
    assert!(matches!(
        v3k_sg_identity_record(&rows, 4, true),
        Err(SdfReadError::Parse(_))
    ));
    let MolBlockRecord::Concrete { topology, .. } =
        v3k_sg_identity_record(&rows, 4, false).expect("cascading non-strict removal")
    else {
        panic!("ordinary SGroup carrier must remain concrete");
    };
    assert_eq!(topology.substance_groups.len(), 1);
    assert_eq!(topology.substance_groups[0].rdkit_sequence_id(), Some(4));
    assert_eq!(topology.substance_groups[0].label(), Some("sibling"));
}

#[test]
fn v3k_sg_parent_malformed_child_drops_only_that_group_non_strictly() {
    let rows = [
        "1 DAT 0 LABEL=parent",
        "2 DAT 0 PARENT=bad LABEL=invalid-child",
        "3 DAT 0 LABEL=sibling",
    ];
    assert!(matches!(
        v3k_sg_identity_record(&rows, 3, true),
        Err(SdfReadError::Parse(message)) if message.contains("Invalid PARENT label")
    ));
    let MolBlockRecord::Concrete { topology, .. } =
        v3k_sg_identity_record(&rows, 3, false).expect("non-strict malformed child")
    else {
        panic!("ordinary SGroup carrier must remain concrete");
    };
    assert_eq!(topology.substance_groups.len(), 2);
    assert_eq!(
        topology
            .substance_groups
            .iter()
            .map(|group| group.label())
            .collect::<Vec<_>>(),
        vec![Some("parent"), Some("sibling")]
    );
}

#[test]
fn v3k_sg_parent_formatted_signs_overflow_and_duplicate_sequence_follow_owner_rules() {
    let MolBlockRecord::Concrete { topology, .. } = v3k_sg_identity_record(
        &["1 DAT 0 LABEL=parent", "2 DAT 0 PARENT=+1 LABEL=child"],
        2,
        true,
    )
    .expect("formatted plus parent") else {
        panic!("ordinary SGroup carrier must remain concrete");
    };
    assert_eq!(
        topology.substance_groups[1].parent(),
        Some(topology.substance_groups[0].id())
    );

    let overflow_rows = ["1 DAT 0 LABEL=parent", "2 DAT 0 PARENT=4294967296"];
    assert!(matches!(
        v3k_sg_identity_record(&overflow_rows, 2, true),
        Err(SdfReadError::Parse(message)) if message.contains("Invalid PARENT label")
    ));
    let MolBlockRecord::Concrete { topology, .. } =
        v3k_sg_identity_record(&overflow_rows, 2, false)
            .expect("non-strict overflow invalidates only the child")
    else {
        panic!("ordinary SGroup carrier must remain concrete");
    };
    assert_eq!(topology.substance_groups.len(), 1);
    assert_eq!(topology.substance_groups[0].label(), Some("parent"));

    // std::map::emplace retains the first sequence-2 row and therefore its
    // parent relation. The duplicate row cannot overwrite temporary metadata.
    let rows = [
        "1 DAT 0 LABEL=first-parent",
        "2 DAT 0 PARENT=1 LABEL=first-child",
        "2 DAT 0 PARENT=3 LABEL=duplicate-child",
        "3 DAT 0 LABEL=second-parent",
    ];
    assert!(matches!(
        v3k_sg_identity_record(&rows, 4, true),
        Err(SdfReadError::Parse(message)) if message.contains("Found 3 SGroups when 4 were expected")
    ));
    let MolBlockRecord::Concrete { topology, .. } =
        v3k_sg_identity_record(&rows, 4, false).expect("non-strict duplicate sequence")
    else {
        panic!("ordinary SGroup carrier must remain concrete");
    };
    assert_eq!(topology.substance_groups[1].label(), Some("first-child"));
    assert_eq!(
        topology.substance_groups[1].parent(),
        Some(topology.substance_groups[0].id())
    );
}

#[test]
fn v3k_sg_parent_cycles_are_retained_as_in_range_typed_source_state() {
    // RDKit stores PARENT without hierarchy-cycle validation and the canonical
    // model validates range rather than acyclicity. Preserve self and two-row
    // cycles as typed relations instead of inventing a new parser policy.
    let MolBlockRecord::Concrete { topology, .. } = v3k_sg_identity_record(
        &[
            "1 DAT 0 PARENT=2 LABEL=one",
            "2 DAT 0 PARENT=1 LABEL=two",
            "3 DAT 0 PARENT=3 LABEL=self",
        ],
        3,
        true,
    )
    .expect("source parent cycles") else {
        panic!("ordinary SGroup carrier must remain concrete");
    };
    assert_eq!(
        topology.substance_groups[0].parent().map(|id| id.index()),
        Some(1)
    );
    assert_eq!(
        topology.substance_groups[1].parent().map(|id| id.index()),
        Some(0)
    );
    assert_eq!(
        topology.substance_groups[2].parent().map(|id| id.index()),
        Some(2)
    );
}

#[test]
fn v3k_sg_bond_refs_removed_target_drops_the_complete_group() {
    // Canonical typed references participate in the model's existing atomic
    // remapping policy: removal of a referenced bond cannot leave a stale row.
    let MolBlockRecord::Concrete { topology, .. } = v3k_sg_member_record(
        "ATOMS=(3 10 20 30) XBONDS=(2 50 60) XBHEAD=(1 2) XBCORR=(1 2)",
        true,
    )
    .expect("valid typed references before compaction") else {
        panic!("ordinary SGroup carrier must remain concrete");
    };
    let mut edit = topology.begin_batch_edit().expect("valid parsed topology");
    edit.remove_bond(BondId::new(1))
        .expect("referenced bond exists");
    let (compacted, _) = edit.finish().expect("atomic SGroup removal");
    assert!(compacted.substance_groups.is_empty());
    compacted
        .validate()
        .expect("no stale bond reference remains");
}

#[test]
fn v3k_sg_brackets_preserve_all_xyz_components_and_append_in_source_order() {
    // Pinned ParseV3000ParseLabel fills three Point3D values from each nine
    // element BRKXYZ array and addBracket appends repeated labels in order.
    for strict_parsing in [true, false] {
        let MolBlockRecord::Concrete { topology, .. } = v3k_sg_array_record(
            concat!(
                "ATOMS=(1 1) BRKXYZ=(9 1 2 3 4 5 6 7 8 9) ",
                "BRKXYZ=(9 -1 -2 -3 -4 -5 -6 -7 -8 -9)"
            ),
            strict_parsing,
        )
        .expect("two complete source brackets") else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        let brackets = topology.substance_groups[0]
            .display()
            .expect("typed display metadata")
            .brackets();
        assert_eq!(brackets.len(), 2);
        assert_eq!(
            brackets[0].points(),
            &[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
        );
        assert_eq!(
            brackets[1].points(),
            &[[-1.0, -2.0, -3.0], [-4.0, -5.0, -6.0], [-7.0, -8.0, -9.0],]
        );
    }
}

#[test]
fn v3k_sg_formatted_double_brkxyz_preserves_stream_failure_destination_reuse() {
    // Pinned RDKit 2026.03.1 calls ParseV3000Array<double>, whose single
    // `double value` destination is reused for every formatted extraction.
    // In the fixed libstdc++ environment an incomplete exponent assigns zero
    // and failbit, while overflow assigns DBL_MAX and failbit; later sentry
    // failures leave that assigned destination unchanged.
    for (value, expected) in [
        ("1e 2 3 4 5 6 7 8 9", [0.0; 9]),
        ("1e+ 2 3 4 5 6 7 8 9", [0.0; 9]),
        (
            "1 2 1e 4 5 6 7 8 9",
            [1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ),
        ("1e309 2 3 4 5 6 7 8 9", [f64::MAX; 9]),
        (
            "1 2 1e309 4 5 6 7 8 9",
            [
                1.0,
                2.0,
                f64::MAX,
                f64::MAX,
                f64::MAX,
                f64::MAX,
                f64::MAX,
                f64::MAX,
                f64::MAX,
            ],
        ),
    ] {
        for strict_parsing in [true, false] {
            let labels = format!("ATOMS=(1 1) BRKXYZ=(9 {value}) LABEL=unreached");
            let MolBlockRecord::Concrete { topology, .. } =
                v3k_sg_array_record(&labels, strict_parsing)
                    .expect("source-formatted BRKXYZ extraction")
            else {
                panic!("ordinary SGroup carrier must remain concrete");
            };
            let group = &topology.substance_groups[0];
            let actual = group.display().unwrap().brackets()[0]
                .points()
                .iter()
                .flatten()
                .copied()
                .collect::<Vec<_>>();
            assert_eq!(actual, expected, "value={value:?}, strict={strict_parsing}");
            assert_eq!(group.label(), None, "failbit stops following labels");
        }
    }
}

#[test]
fn v3k_sg_formatted_double_accepts_source_finite_spellings_and_c_locale_whitespace() {
    let labels = "ATOMS=(1 1) BRKXYZ=(9\t+1\x0b-2\x0c.5\t1.\x0b1e2\x0c-1E-2 0 -0 3) LABEL=after";
    let MolBlockRecord::Concrete { topology, .. } =
        v3k_sg_array_record(labels, true).expect("ordinary formatted doubles")
    else {
        panic!("ordinary SGroup carrier must remain concrete");
    };
    let group = &topology.substance_groups[0];
    let actual = group.display().unwrap().brackets()[0]
        .points()
        .iter()
        .flatten()
        .map(|value| value.to_bits())
        .collect::<Vec<_>>();
    let expected = [1.0_f64, -2.0, 0.5, 1.0, 100.0, -0.01, 0.0, -0.0, 3.0].map(f64::to_bits);
    assert_eq!(actual, expected);
    assert_eq!(group.label(), Some("after"));
}

#[test]
fn v3k_sg_brackets_require_exactly_nine_values_and_invalidate_the_whole_group() {
    // The source branch rejects both a parsed vector shorter than nine and an
    // over-limit declared array; non-strict parsing drops the invalid SGroup.
    for labels in [
        "ATOMS=(1 1) BRKXYZ=(8 1 2 3 4 5 6 7 8)",
        "ATOMS=(1 1) BRKXYZ=(10 1 2 3 4 5 6 7 8 9 10)",
    ] {
        assert!(matches!(
            v3k_sg_array_record(labels, true),
            Err(SdfReadError::Parse(_))
        ));
        let MolBlockRecord::Concrete { topology, .. } =
            v3k_sg_array_record(labels, false).expect("invalid SGroup is discarded")
        else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        assert!(topology.substance_groups.is_empty());
    }
}

#[test]
fn v3k_sg_brackets_do_not_reclassify_effective_2d_molecular_coordinates() {
    // BRKXYZ z components are SGroup display metadata, not atom coordinates.
    let MolBlockRecord::Concrete {
        topology,
        coordinates,
        ..
    } = v3k_sg_array_record("ATOMS=(1 1) BRKXYZ=(9 1 2 3 4 5 6 7 8 9)", true)
        .expect("3D bracket metadata on an effective-2D molecule")
    else {
        panic!("ordinary SGroup carrier must remain concrete");
    };
    assert_eq!(
        coordinates.source_coordinate_dim,
        Some(CoordinateDimension::TwoD)
    );
    assert_eq!(coordinates.conformers_2d.len(), 1);
    assert!(coordinates.conformers_3d.is_empty());
    assert_eq!(
        topology.substance_groups[0].display().unwrap().brackets()[0].points()[2],
        [7.0, 8.0, 9.0]
    );
}

#[test]
fn v3k_sg_cstate_preserves_sup_xyz_and_non_sup_uses_the_source_zero_vector() {
    // ParseV3000CStateLabel reads all three Point3D components only for SUP;
    // every other type passes the default-constructed zero Point3D.
    for strict_parsing in [true, false] {
        let MolBlockRecord::Concrete {
            topology,
            coordinates,
            ..
        } = v3k_sg_cstate_record(
            "SUP",
            "ATOMS=(1 20) XBONDS=(1 50) CSTATE=(4 50 1.25 -2.5 3.75)",
            strict_parsing,
        )
        .expect("complete SUP CSTATE")
        else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        let cstate = &topology.substance_groups[0].cstates()[0];
        assert_eq!(cstate.bond(), BondId::new(0));
        assert_eq!(cstate.vector(), &[1.25, -2.5, 3.75]);
        assert_eq!(
            coordinates.source_coordinate_dim,
            Some(CoordinateDimension::TwoD)
        );
        assert_eq!(coordinates.conformers_2d.len(), 1);
        assert!(coordinates.conformers_3d.is_empty());

        let MolBlockRecord::Concrete { topology, .. } = v3k_sg_cstate_record(
            "DAT",
            "ATOMS=(1 20) XBONDS=(1 50) CSTATE=(1 50)",
            strict_parsing,
        )
        .expect("non-SUP default CSTATE vector") else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        assert_eq!(
            topology.substance_groups[0].cstates()[0].vector(),
            &[0.0, 0.0, 0.0]
        );
    }
}

#[test]
fn v3k_sg_cstate_formatted_extraction_zero_fills_from_the_first_missing_or_bad_component() {
    // Point3D begins at zero and C++ formatted extraction leaves the stream in
    // fail state, so later components cannot resume after the first failure.
    for (value, expected) in [
        ("CSTATE=(4 50 1 2)", [1.0, 2.0, 0.0]),
        ("CSTATE=(4 50 1)", [1.0, 0.0, 0.0]),
        ("CSTATE=(4 50)", [0.0, 0.0, 0.0]),
        ("CSTATE=(4 50 x 2 3)", [0.0, 0.0, 0.0]),
        ("CSTATE=(4 50 1 x 3)", [1.0, 0.0, 0.0]),
        ("CSTATE=(4 50 1 2 x)", [1.0, 2.0, 0.0]),
    ] {
        for strict_parsing in [true, false] {
            let labels = format!("ATOMS=(1 20) XBONDS=(1 50) {value}");
            let MolBlockRecord::Concrete { topology, .. } =
                v3k_sg_cstate_record("SUP", &labels, strict_parsing)
                    .expect("source-formatted CSTATE extraction")
            else {
                panic!("ordinary SGroup carrier must remain concrete");
            };
            assert_eq!(
                topology.substance_groups[0].cstates()[0].vector(),
                &expected,
                "{value}, strict={strict_parsing}"
            );
        }
    }
}

#[test]
fn v3k_sg_formatted_double_cstate_preserves_malformed_exponent_and_overflow_fail_state() {
    for (components, expected) in [
        ("1e 2 3", [0.0; 3]),
        ("1e+ 2 3", [0.0; 3]),
        ("1 1e 3", [1.0, 0.0, 0.0]),
        ("1e309 2 3", [f64::MAX, 0.0, 0.0]),
        ("1 1e309 3", [1.0, f64::MAX, 0.0]),
    ] {
        for strict_parsing in [true, false] {
            let labels =
                format!("ATOMS=(1 20) XBONDS=(1 50) CSTATE=(4 50 {components}) LABEL=unreached");
            let MolBlockRecord::Concrete { topology, .. } =
                v3k_sg_cstate_record("SUP", &labels, strict_parsing)
                    .expect("source-formatted CSTATE extraction")
            else {
                panic!("ordinary SGroup carrier must remain concrete");
            };
            let group = &topology.substance_groups[0];
            assert_eq!(
                group.cstates()[0].vector(),
                &expected,
                "components={components:?}, strict={strict_parsing}"
            );
            assert_eq!(group.label(), None, "failbit stops following labels");
        }
    }
}

#[test]
fn v3k_sg_cstate_wrong_counts_fail_strict_and_drop_the_whole_group_non_strict() {
    for (kind, labels) in [
        ("SUP", "ATOMS=(1 20) XBONDS=(1 50) CSTATE=(3 50 1 2 3)"),
        ("DAT", "ATOMS=(1 20) XBONDS=(1 50) CSTATE=(4 50 1 2 3)"),
    ] {
        assert!(matches!(
            v3k_sg_cstate_record(kind, labels, true),
            Err(SdfReadError::Parse(_))
        ));
        let MolBlockRecord::Concrete { topology, .. } =
            v3k_sg_cstate_record(kind, labels, false).expect("invalid SGroup is discarded")
        else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        assert!(topology.substance_groups.is_empty());
    }
}

#[test]
fn v3k_sg_cstate_requires_a_unique_listed_actual_crossing_bond() {
    // addCState/getBondType checks resolved bond membership and endpoints; an
    // XBONDS spelling cannot make a contained bond into an actual crossing.
    for labels in [
        "ATOMS=(1 20) XBONDS=(1 50) CSTATE=(4 99 1 2 3)",
        "ATOMS=(1 20) CSTATE=(4 50 1 2 3)",
        "ATOMS=(2 10 20) CBONDS=(1 50) CSTATE=(4 50 1 2 3)",
        "ATOMS=(2 10 20) XBONDS=(1 50) CSTATE=(4 50 1 2 3)",
        "CSTATE=(4 50 1 2 3) ATOMS=(1 20) XBONDS=(1 50)",
    ] {
        assert!(matches!(
            v3k_sg_cstate_record("SUP", labels, true),
            Err(SdfReadError::Parse(_))
        ));
        let MolBlockRecord::Concrete { topology, .. } =
            v3k_sg_cstate_record("SUP", labels, false).expect("invalid SGroup is discarded")
        else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        assert!(topology.substance_groups.is_empty());
    }
}

#[test]
fn v3k_sg_sap_resolves_aidx_zero_and_explicit_nonsequential_bookmarks() {
    // Pinned ParseV3000SAPLabel uppercases AIDX, treats toInt zero as no
    // leaving atom, and resolves all nonzero source bookmarks before the
    // append-only addAttachPoint call. The attachment need not be in ATOMS.
    let MolBlockRecord::Concrete { topology, .. } =
        v3k_sg_sap_record("SAP=(2 30 aidx lo) SAP=(3 10 0 ZZ) SAP=(4 20 10 XY)", true)
            .expect("source-valid SAP rows")
    else {
        panic!("concrete record expected");
    };
    let group = &topology.substance_groups[0];
    assert!(group.atoms().is_empty());
    assert_eq!(group.attach_points().len(), 3);
    assert_eq!(group.attach_points()[0].atom, AtomId::new(2));
    assert_eq!(group.attach_points()[0].leaving_atom, Some(AtomId::new(2)));
    assert_eq!(group.attach_points()[0].label.as_deref(), Some("lo"));
    assert_eq!(group.attach_points()[1].atom, AtomId::new(0));
    assert_eq!(group.attach_points()[1].leaving_atom, None);
    assert_eq!(group.attach_points()[1].label.as_deref(), Some("ZZ"));
    assert_eq!(group.attach_points()[2].atom, AtomId::new(1));
    assert_eq!(group.attach_points()[2].leaving_atom, Some(AtomId::new(0)));
    assert_eq!(group.attach_points()[2].label.as_deref(), Some("XY"));
}

#[test]
fn v3k_sg_sap_preserves_formatted_conversion_and_label_last_byte_behavior() {
    // The count is deliberately ignored; stream.get() discards any opening
    // byte; toInt does not accept leading '+'. Missing ')' still causes the
    // source pop_back to remove the last label byte, after which parsing
    // continues at the next outer label.
    let MolBlockRecord::Concrete { topology, .. } =
        v3k_sg_sap_record("SAP=[7 +10 +20 PP) SAP=(3 10 20 AP SAP=(3 20 0 AB))", true)
            .expect("source-valid formatted SAP edge behavior")
    else {
        panic!("concrete record expected");
    };
    let points = topology.substance_groups[0].attach_points();
    assert_eq!(points.len(), 3);
    assert_eq!(points[0].atom, AtomId::new(0));
    assert_eq!(points[0].leaving_atom, None);
    assert_eq!(points[0].label.as_deref(), Some("PP"));
    assert_eq!(points[1].atom, AtomId::new(0));
    assert_eq!(points[1].leaving_atom, Some(AtomId::new(1)));
    assert_eq!(points[1].label.as_deref(), Some("A"));
    assert_eq!(points[2].atom, AtomId::new(1));
    assert_eq!(points[2].leaving_atom, None);
    assert_eq!(points[2].label.as_deref(), Some("AB)"));
}

#[test]
fn v3k_sg_sap_empty_and_utf8_labels_use_the_canonical_text_boundary() {
    let MolBlockRecord::Concrete { topology, .. } =
        v3k_sg_sap_record("SAP=(3 10 20 ) SAP=(3 20 0 é)", true)
            .expect("empty and UTF-8 labels with closing parentheses")
    else {
        panic!("concrete record expected");
    };
    let points = topology.substance_groups[0].attach_points();
    assert_eq!(points[0].label.as_deref(), Some(""));
    assert_eq!(points[1].label.as_deref(), Some("é"));

    assert!(matches!(
        v3k_sg_sap_record("SAP=(3 10 20 é", true),
        Err(SdfReadError::Parse(_))
    ));
    let MolBlockRecord::Concrete { topology, .. } =
        v3k_sg_sap_record("SAP=(3 10 20 é", false).expect("non-strict invalid group is discarded")
    else {
        panic!("concrete record expected");
    };
    assert!(topology.substance_groups.is_empty());
}

#[test]
fn v3k_sg_sap_malformed_fields_and_missing_bookmarks_invalidate_the_whole_group() {
    for labels in [
        "SAP=(x 10 20 AP)",
        "SAP=(3 x 20 AP)",
        "SAP=(3 10 x20 AP)",
        "SAP=(3 10 20",
        "SAP=(3 99 20 AP)",
        "SAP=(3 10 99 AP)",
    ] {
        assert!(
            matches!(v3k_sg_sap_record(labels, true), Err(SdfReadError::Parse(_))),
            "strict case {labels:?}"
        );
        let MolBlockRecord::Concrete { topology, .. } =
            v3k_sg_sap_record(labels, false).expect("non-strict invalid group is discarded")
        else {
            panic!("concrete record expected");
        };
        assert!(
            topology.substance_groups.is_empty(),
            "non-strict case {labels:?}"
        );
    }
}

#[test]
fn v3k_sg_scalar_labels_component_number_accepts_zero_and_256_but_rejects_257() {
    // Pinned ParseV3000ParseLabel performs formatted unsigned extraction and
    // rejects only values greater than 256 after conversion.
    for (source, expected) in [("0", 0), ("256", 256), ("+256", 256)] {
        for strict_parsing in [true, false] {
            let MolBlockRecord::Concrete { topology, .. } =
                v3k_sg_string_record(&format!("COMPNO={source}"), strict_parsing)
                    .expect("source-valid component number")
            else {
                panic!("concrete record expected");
            };
            assert_eq!(
                topology.substance_groups[0].component_number(),
                Some(expected),
                "source value {source:?}"
            );
        }
    }

    assert!(matches!(
        v3k_sg_string_record("COMPNO=257", true),
        Err(SdfReadError::Parse(message)) if message.contains("over 256")
    ));
    let MolBlockRecord::Concrete { topology, .. } =
        v3k_sg_string_record("COMPNO=257", false).expect("invalid group is discarded")
    else {
        panic!("concrete record expected");
    };
    assert!(topology.substance_groups.is_empty());
}

#[test]
fn v3k_sg_scalar_labels_accept_every_source_vocabulary_value() {
    // These are the complete case-sensitive vectors in pinned
    // SubstanceGroupChecks. The canonical model retains subtype/class text
    // and represents each accepted connection as its typed enum state.
    for subtype in ["ALT", "RAN", "BLO"] {
        let MolBlockRecord::Concrete { topology, .. } =
            v3k_sg_string_record(&format!("SUBTYPE={subtype}"), true)
                .expect("source-valid subtype")
        else {
            panic!("concrete record expected");
        };
        assert_eq!(topology.substance_groups[0].subtype(), Some(subtype));
    }

    for (source, expected) in [
        ("HH", SGroupConnection::HeadToHead),
        ("HT", SGroupConnection::HeadToTail),
        ("EU", SGroupConnection::Either),
    ] {
        let MolBlockRecord::Concrete { topology, .. } =
            v3k_sg_string_record(&format!("CONNECT={source}"), true)
                .expect("source-valid connection")
        else {
            panic!("concrete record expected");
        };
        assert_eq!(topology.substance_groups[0].connection(), Some(&expected));
    }

    for class in [
        "AA",
        "dAA",
        "DNA",
        "RNA",
        "SUGAR",
        "BASE",
        "PHOSPHATE",
        "LINKER",
        "CHEM",
        "LGRP",
        "MODAA",
        "MODdAA",
        "MODDNA",
        "MODRNA",
        "XLINKAA",
        "XLINKdAA",
        "XLINKDNA",
        "XLINKRNA",
    ] {
        let MolBlockRecord::Concrete { topology, .. } =
            v3k_sg_string_record(&format!("CLASS={class}"), true)
                .expect("source-valid template class")
        else {
            panic!("concrete record expected");
        };
        assert_eq!(topology.substance_groups[0].class(), Some(class));
    }
}

#[test]
fn v3k_sg_scalar_labels_invalid_case_or_spelling_invalidates_the_whole_group() {
    // Validation uses exact std::string equality; it does not uppercase or
    // otherwise normalize these source values.
    for labels in [
        "SUBTYPE=alt",
        "SUBTYPE=RAND",
        "CONNECT=hh",
        "CONNECT=HEADTAIL",
        "CLASS=DAA",
        "CLASS=XLINKrna",
    ] {
        assert!(
            matches!(
                v3k_sg_string_record(labels, true),
                Err(SdfReadError::Parse(_))
            ),
            "strict case {labels:?}"
        );
        let MolBlockRecord::Concrete { topology, .. } =
            v3k_sg_string_record(labels, false).expect("invalid group is discarded")
        else {
            panic!("concrete record expected");
        };
        assert!(
            topology.substance_groups.is_empty(),
            "non-strict case {labels:?}"
        );
    }
}

#[test]
fn v3k_sg_strings_empty_value_retains_the_following_label() {
    // Pinned ParseV3000StringPropLabel returns without consuming a literal
    // space, leaving that separator and the following label to the outer loop.
    for strict_parsing in [true, false] {
        let MolBlockRecord::Concrete { topology, .. } =
            v3k_sg_string_record("LABEL= ESTATE=E", strict_parsing)
                .expect("source empty-value separator")
        else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        let group = &topology.substance_groups[0];
        assert_eq!(group.label(), Some(""));
        assert_eq!(group.expansion_state(), Some("E"));
    }
}

#[test]
fn v3k_sg_strings_doubled_quotes_equals_and_trailing_trim_match_source() {
    // The pinned string helper collapses doubled double quotes, retains equals
    // signs as value bytes, trims the quoted value's trailing C whitespace,
    // and leaves the next label visible to the enclosing parser.
    for strict_parsing in [true, false] {
        let MolBlockRecord::Concrete { topology, .. } =
            v3k_sg_string_record("LABEL=\"a\"\"b=c   \" ESTATE=E", strict_parsing)
                .expect("quoted source string")
        else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        let group = &topology.substance_groups[0];
        assert_eq!(group.label(), Some("a\"b=c"));
        assert_eq!(group.expansion_state(), Some("E"));

        let MolBlockRecord::Concrete { topology, .. } =
            v3k_sg_string_record("LABEL=a=b ESTATE=E", strict_parsing)
                .expect("unquoted equals sign")
        else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        let group = &topology.substance_groups[0];
        assert_eq!(group.label(), Some("a=b"));
        assert_eq!(group.expansion_state(), Some("E"));
    }
}

#[test]
fn v3k_sg_strings_single_quote_consumption_matches_source_error_boundary() {
    // Source getline starts on the opening single quote, returns an empty
    // value, and leaves the next byte for the outer literal-space check.
    assert!(v3k_sg_string_record("LABEL='abc' ESTATE=E", true).is_err());
    let MolBlockRecord::Concrete { topology, .. } =
        v3k_sg_string_record("LABEL='abc' ESTATE=E", false)
            .expect("non-strict invalid SGroup is discarded")
    else {
        panic!("ordinary SGroup carrier must remain concrete");
    };
    assert!(topology.substance_groups.is_empty());
}

#[test]
fn v3k_sg_strings_tabs_are_value_whitespace_but_not_label_separators() {
    // Formatted string extraction skips C-locale tab whitespace after '=',
    // while the enclosing label loop accepts only a literal ASCII space.
    for strict_parsing in [true, false] {
        let MolBlockRecord::Concrete { topology, .. } =
            v3k_sg_string_record("LABEL=\tvalue ESTATE=E", strict_parsing)
                .expect("tab before an unquoted string value")
        else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        let group = &topology.substance_groups[0];
        assert_eq!(group.label(), Some("value"));
        assert_eq!(group.expansion_state(), Some("E"));
    }

    assert!(v3k_sg_string_record("LABEL=value\tESTATE=E", true).is_err());
    let MolBlockRecord::Concrete { topology, .. } =
        v3k_sg_string_record("LABEL=value\tESTATE=E", false)
            .expect("non-strict invalid SGroup is discarded")
    else {
        panic!("ordinary SGroup carrier must remain concrete");
    };
    assert!(topology.substance_groups.is_empty());
}

#[test]
fn v3k_sg_arrays_zero_and_max_counts_preserve_typed_state_in_both_modes() {
    // Pinned ParseV3000Array reads the count from the shared stream and accepts
    // zero as well as a count equal to the caller's maximum.
    for strict_parsing in [true, false] {
        let MolBlockRecord::Concrete { topology, .. } =
            v3k_sg_array_record("ATOMS=(0) LABEL=empty", strict_parsing).expect("zero-count array")
        else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        assert_eq!(topology.substance_groups.len(), 1);
        assert!(topology.substance_groups[0].atoms().is_empty());
        assert_eq!(topology.substance_groups[0].label(), Some("empty"));

        let MolBlockRecord::Concrete { topology, .. } =
            v3k_sg_array_record("ATOMS=(2 1 2) LABEL=full", strict_parsing)
                .expect("maximum-count array")
        else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        let group = &topology.substance_groups[0];
        assert_eq!(
            group
                .atoms()
                .iter()
                .map(|atom| atom.index())
                .collect::<Vec<_>>(),
            vec![0, 1]
        );
        assert_eq!(group.label(), Some("full"));

        let MolBlockRecord::Concrete { topology, .. } =
            v3k_sg_array_record("BRKXYZ=(9 1 2 3 4 5 6 7 8 9)", strict_parsing)
                .expect("maximum double array")
        else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        assert_eq!(
            topology.substance_groups[0].display().unwrap().brackets[0].points,
            [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
        );
    }
}

#[test]
fn v3k_sg_arrays_over_limit_counts_follow_strict_warning_boundary() {
    // SGroupWarnOrThrow throws in strict mode and returns an empty vector in
    // non-strict mode without consuming the declared elements.
    assert!(matches!(
        v3k_sg_array_record("ATOMS=(3 1 2 1)", true),
        Err(SdfReadError::Parse(message)) if message == "invalid count value"
    ));
    let MolBlockRecord::Concrete { topology, .. } =
        v3k_sg_array_record("ATOMS=(3 1 2 1)", false).expect("warned count")
    else {
        panic!("ordinary SGroup carrier must remain concrete");
    };
    assert_eq!(topology.substance_groups.len(), 1);
    assert!(topology.substance_groups[0].atoms().is_empty());

    assert!(v3k_sg_array_record("BRKXYZ=(10 1 2 3 4 5 6 7 8 9 10)", true).is_err());
    let MolBlockRecord::Concrete { topology, .. } =
        v3k_sg_array_record("BRKXYZ=(10 1 2 3 4 5 6 7 8 9 10)", false)
            .expect("non-strict invalid bracket group is discarded")
    else {
        panic!("ordinary SGroup carrier must remain concrete");
    };
    assert!(topology.substance_groups.is_empty());
}

#[test]
fn v3k_sg_arrays_parentheses_and_shared_cursor_match_source() {
    // Parentheses are consumed-and-warned delimiters, not validation gates.
    // The final get() and the outer label loop share one cursor.
    for strict_parsing in [true, false] {
        for labels in [
            "ATOMS=[2 1 2) LABEL=kept",
            "ATOMS=(2\t1\u{000b}2) LABEL=kept",
        ] {
            let MolBlockRecord::Concrete { topology, .. } =
                v3k_sg_array_record(labels, strict_parsing).expect("source-shaped delimiters")
            else {
                panic!("ordinary SGroup carrier must remain concrete");
            };
            let group = &topology.substance_groups[0];
            assert_eq!(
                group
                    .atoms()
                    .iter()
                    .map(|atom| atom.index())
                    .collect::<Vec<_>>(),
                vec![0, 1]
            );
            assert_eq!(group.label(), Some("kept"));
        }

        let MolBlockRecord::Concrete { topology, .. } =
            v3k_sg_array_record("ATOMS=2 1 2", strict_parsing)
                .expect("missing opening parenthesis only warns")
        else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        assert_eq!(topology.substance_groups[0].atoms()[0].index(), 1);
    }

    assert!(v3k_sg_array_record("ATOMS=(2 1 2)\tLABEL=lost", true).is_err());
    let MolBlockRecord::Concrete { topology, .. } =
        v3k_sg_array_record("ATOMS=(2 1 2)\tLABEL=lost", false)
            .expect("non-strict separator failure discards group")
    else {
        panic!("ordinary SGroup carrier must remain concrete");
    };
    assert!(topology.substance_groups.is_empty());

    assert!(v3k_sg_array_record("ATOMS=(2 1 2 LABEL=lost", true).is_err());
    let MolBlockRecord::Concrete { topology, .. } =
        v3k_sg_array_record("ATOMS=(2 1 2 LABEL=lost", false)
            .expect("consumed non-parenthesis exposes separator failure")
    else {
        panic!("ordinary SGroup carrier must remain concrete");
    };
    assert!(topology.substance_groups.is_empty());
}

#[test]
fn v3k_sg_arrays_too_few_and_extra_values_do_not_fabricate_success() {
    // A failed arithmetic extraction assigns zero and sets failbit, so
    // bookmark resolution rejects the result. An extra scalar remains for the
    // outer separator loop and invalidates the group.
    for labels in ["ATOMS=(2 1)", "ATOMS=(1 x)", "ATOMS=(1 1 2)"] {
        assert!(v3k_sg_array_record(labels, true).is_err());
        let MolBlockRecord::Concrete { topology, .. } =
            v3k_sg_array_record(labels, false).expect("invalid group is discarded")
        else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        assert!(topology.substance_groups.is_empty());
    }

    for strict_parsing in [true, false] {
        let MolBlockRecord::Concrete { topology, .. } =
            v3k_sg_array_record("ATOMS=(2 1 2", strict_parsing)
                .expect("missing trailing parenthesis only warns at EOF")
        else {
            panic!("ordinary SGroup carrier must remain concrete");
        };
        assert_eq!(topology.substance_groups[0].atoms().len(), 2);
    }
}

#[test]
fn v3k_attachment_props_scalar_values_remain_distinct_and_share_one_atom() {
    // ParseV3000AtomProps stores ATTCHPT and scalar ATTCHORD in distinct atom
    // properties. The literal ATTCHPT value `0` alone is skipped; other
    // spellings still pass through FileParserUtils::toInt.
    let block = v3000_block(&["M  V30 42 C 0 0 0 0 ATTCHPT=3 ATTCHORD=4"], 1);
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("scalar attachment properties")
    else {
        panic!("concrete attachment carrier expected");
    };
    let atom = &topology.atoms[0];
    assert_eq!(atom.prop("molAttachPoint"), Some("3"));
    assert_eq!(atom.prop("molAttachOrder"), Some("4"));
    assert!(atom.template_attachment_order().is_none());

    let literal_zero = v3000_block(&["M  V30 1 C 0 0 0 0 ATTCHPT=0"], 1);
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&literal_zero).expect("literal zero is skipped")
    else {
        panic!("concrete attachment carrier expected");
    };
    assert_eq!(topology.atoms[0].prop("molAttachPoint"), None);

    let converted_zero = v3000_block(&["M  V30 1 C 0 0 0 0 ATTCHPT=00"], 1);
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&converted_zero).expect("nonliteral zero is stored")
    else {
        panic!("concrete attachment carrier expected");
    };
    assert_eq!(topology.atoms[0].prop("molAttachPoint"), Some("0"));
}

#[test]
fn v3k_attachment_props_template_order_uses_rows_not_nonsequential_bookmarks() {
    // The frozen SCSRMolFileParser consumer treats the stored one-based values
    // as atom-table rows. They are not V3000 atom bookmarks, and forward rows
    // remain valid until the complete atom table is available.
    let block = v3000_block(
        &[
            "M  V30 42 C 0 0 0 0 ATTCHPT=2 ATTCHORD=(4 3 Br 2 Al)",
            "M  V30 7 N 1 0 0 0",
            "M  V30 99 O 2 0 0 0",
        ],
        3,
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("ordered forward template attachments")
    else {
        panic!("concrete attachment carrier expected");
    };
    let carrier = &topology.atoms[0];
    assert_eq!(carrier.prop("molAttachPoint"), Some("2"));
    let entries = carrier
        .template_attachment_order()
        .expect("typed template attachment order")
        .entries();
    assert_eq!(entries.len(), 2);
    assert_eq!((entries[0].target().index(), entries[0].label()), (2, "Br"));
    assert_eq!((entries[1].target().index(), entries[1].label()), (1, "Al"));
}

#[test]
fn v3k_attachment_props_template_empty_labels_and_repeated_separators_follow_source() {
    // boost::split without token compression preserves the final empty label,
    // while an extra separator creates an extra field and fails the count.
    let block = v3000_block(
        &["M  V30 1 C 0 0 0 0 ATTCHORD=(2 2 )", "M  V30 2 N 1 0 0 0"],
        2,
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("empty template label")
    else {
        panic!("concrete attachment carrier expected");
    };
    let entries = topology.atoms[0]
        .template_attachment_order()
        .expect("typed template attachment order")
        .entries();
    assert_eq!((entries[0].target().index(), entries[0].label()), (1, ""));

    for value in ["(2  2 Al)", "(2\t\t2 Al)", "(2 \t 2 Al)"] {
        let block = v3000_block(
            &[
                &format!("M  V30 1 C 0 0 0 0 ATTCHORD={value}"),
                "M  V30 2 N 1 0 0 0",
            ],
            2,
        );
        assert!(
            matches!(
                read_mol_block_detached(&block),
                Err(SdfReadError::Parse(message)) if message.contains("Invalid ATTCHORD value")
            ),
            "repeated separators in {value:?} must be rejected"
        );
    }
}

#[test]
fn v3k_attachment_props_duplicate_policies_are_source_shaped() {
    let duplicate_point = v3000_block(&["M  V30 1 C 0 0 0 0 ATTCHPT=1 ATTCHPT=2"], 1);
    assert!(matches!(
        read_mol_block_detached(&duplicate_point),
        Err(SdfReadError::Parse(message)) if message.contains("Multiple ATTCHPT values")
    ));
    let params = cosmolkit_io::MolBlockReadParams {
        strict_parsing: false,
        ..cosmolkit_io::MolBlockReadParams::default()
    };
    let MolBlockRecord::Concrete { topology, .. } =
        cosmolkit_io::read_mol_block_detached_with_params(&duplicate_point, params)
            .expect("non-strict duplicate ATTCHPT keeps the first value")
    else {
        panic!("concrete attachment carrier expected");
    };
    assert_eq!(topology.atoms[0].prop("molAttachPoint"), Some("1"));

    for value in ["(4 2 Al 2 Br)", "(4 2 Al 3 Al)"] {
        let block = v3000_block(
            &[
                &format!("M  V30 1 C 0 0 0 0 ATTCHORD={value}"),
                "M  V30 2 N 1 0 0 0",
                "M  V30 3 O 2 0 0 0",
            ],
            3,
        );
        assert!(matches!(
            read_mol_block_detached(&block),
            Err(SdfReadError::Parse(message)) if message.contains("Invalid ATTCHORD value")
        ));
    }
}

#[test]
fn v3k_attachment_props_targets_and_query_carriers_are_typed_or_fail_structurally() {
    let query = v3000_block(
        &[
            "M  V30 10 * 0 0 0 0 ATTCHPT=5 ATTCHORD=(2 2 Cx)",
            "M  V30 80 N 1 0 0 0",
        ],
        2,
    );
    let MolBlockRecord::Query(record) =
        read_mol_block_detached(&query).expect("query attachment carrier")
    else {
        panic!("wildcard carrier must produce query topology");
    };
    let carrier = &record.query.atoms()[0];
    assert_eq!(carrier.prop("molAttachPoint"), Some("5"));
    let entries = carrier
        .template_attachment_order()
        .expect("query carrier preserves typed attachment order")
        .entries();
    assert_eq!((entries[0].target().index(), entries[0].label()), (1, "Cx"));

    for value in ["(2 0 Al)", "(2 3 Al)"] {
        let concrete = v3000_block(
            &[
                &format!("M  V30 5 C 0 0 0 0 ATTCHORD={value}"),
                "M  V30 90 N 1 0 0 0",
            ],
            2,
        );
        let result = read_mol_block_detached(&concrete);
        if value.contains(" 0 ") {
            assert!(matches!(result, Err(SdfReadError::Parse(_))));
        } else {
            assert!(matches!(result, Err(SdfReadError::Topology(_))));
        }
    }

    let query_out_of_range = v3000_block(
        &[
            "M  V30 5 * 0 0 0 0 ATTCHORD=(2 3 Al)",
            "M  V30 90 N 1 0 0 0",
        ],
        2,
    );
    assert!(matches!(
        read_mol_block_detached(&query_out_of_range),
        Err(SdfReadError::QueryGraph(_))
    ));
}

#[test]
fn v3k_atom_text_props_raw_empty_and_lowercase_names_follow_source() {
    // splitAssignToken uppercases property names only. Values are not trimmed
    // or normalized: CLASS stores empty text, while SEQNAME skips it.
    let lower = v3000_block(&["M  V30 1 C 0 0 0 0 class=alpha seqname=gly"], 1);
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&lower).expect("lowercase property names")
    else {
        panic!("concrete text-property carrier expected");
    };
    assert_eq!(topology.atoms[0].prop("molAtomClass"), Some("alpha"));
    assert_eq!(topology.atoms[0].prop("molAtomSeqName"), Some("gly"));

    let empty = v3000_block(&["M  V30 1 C 0 0 0 0 CLASS= SEQNAME="], 1);
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&empty).expect("empty text-property values")
    else {
        panic!("concrete text-property carrier expected");
    };
    assert_eq!(topology.atoms[0].prop("molAtomClass"), Some(""));
    assert_eq!(topology.atoms[0].prop("molAtomSeqName"), None);
}

#[test]
fn v3k_atom_text_props_repeated_labels_replace_in_source_order() {
    let block = v3000_block(
        &["M  V30 1 C 0 0 0 0 CLASS=A SEQNAME=B CLASS=C SEQNAME=D"],
        1,
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("repeated text properties")
    else {
        panic!("concrete text-property carrier expected");
    };
    assert_eq!(topology.atoms[0].prop("molAtomClass"), Some("C"));
    assert_eq!(topology.atoms[0].prop("molAtomSeqName"), Some("D"));
}

#[test]
fn v3k_atom_text_props_quoted_and_parenthesized_tokens_are_not_normalized() {
    // tokenizeV3000Line does not turn quoted assignment values into ordinary
    // CLASS/SEQNAME assignments. Parenthesized values, by contrast, remain
    // one raw token and retain their delimiters verbatim.
    let quoted = v3000_block(&["M  V30 1 C 0 0 0 0 CLASS=\"AA\" SEQNAME=\"GLY\""], 1);
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&quoted).expect("source-shaped quoted tokens")
    else {
        panic!("concrete text-property carrier expected");
    };
    assert_eq!(topology.atoms[0].prop("molAtomClass"), None);
    assert_eq!(topology.atoms[0].prop("molAtomSeqName"), None);

    let parenthesized = v3000_block(&["M  V30 1 C 0 0 0 0 CLASS=(A B)"], 1);
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&parenthesized).expect("raw parenthesized class")
    else {
        panic!("concrete text-property carrier expected");
    };
    assert_eq!(topology.atoms[0].prop("molAtomClass"), Some("(A B)"));
}

#[test]
fn v3k_atom_text_props_unknown_labels_are_ignored_but_invalid_assignments_fail() {
    let unknown = v3000_block(&["M  V30 1 C 0 0 0 0 UNKNOWN=value CLASS=A OTHER="], 1);
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&unknown).expect("valid unknown assignments are ignored")
    else {
        panic!("concrete text-property carrier expected");
    };
    assert_eq!(topology.atoms[0].prop("molAtomClass"), Some("A"));
    assert_eq!(topology.atoms[0].prop("UNKNOWN"), None);
    assert_eq!(topology.atoms[0].prop("OTHER"), None);

    for token in ["NOEQUALS", "CLASS=A=B"] {
        let block = v3000_block(&[&format!("M  V30 1 C 0 0 0 0 {token}")], 1);
        assert!(matches!(
            read_mol_block_detached(&block),
            Err(SdfReadError::Parse(message))
                if message.contains("Invalid atom property") && message.contains(token)
        ));
    }
}

fn review_collection_record(blocks: &[&str]) -> Result<MolBlockRecord, SdfReadError> {
    let mut extra = String::new();
    for lines in blocks {
        extra.push_str("M  V30 BEGIN COLLECTION\n");
        for line in lines.split('\n').filter(|line| !line.is_empty()) {
            extra.push_str(&format!("M  V30 {line}\n"));
        }
        extra.push_str("M  V30 END COLLECTION\n");
    }
    let block = v3000_block(&["M  V30 1 C 0 0 0 0", "M  V30 2 O 1 0 0 0"], 2)
        .replace("M  V30 END CTAB", &format!("{extra}M  V30 END CTAB"));
    read_mol_block_detached(&block)
}

fn v3k_collections_record(
    atom_lines: &[&str],
    blocks: &[&str],
    strict_parsing: bool,
) -> Result<MolBlockRecord, SdfReadError> {
    let mut extra = String::new();
    for lines in blocks {
        extra.push_str("M  V30 BEGIN COLLECTION\n");
        for line in lines.split('\n').filter(|line| !line.is_empty()) {
            extra.push_str("M  V30 ");
            extra.push_str(line);
            extra.push('\n');
        }
        extra.push_str("M  V30 END COLLECTION\n");
    }
    let block = v3000_block(atom_lines, atom_lines.len())
        .replace("M  V30 END CTAB", &format!("{extra}M  V30 END CTAB"));
    cosmolkit_io::read_mol_block_detached_with_params(
        &block,
        cosmolkit_io::MolBlockReadParams {
            strict_parsing,
            ..cosmolkit_io::MolBlockReadParams::default()
        },
    )
}

fn v3k_collections_raw_record(
    atom_lines: &[&str],
    collection_text: &str,
) -> Result<MolBlockRecord, SdfReadError> {
    let block = v3000_block(atom_lines, atom_lines.len()).replace(
        "M  V30 END CTAB",
        &format!("{collection_text}M  V30 END CTAB"),
    );
    read_mol_block_detached(&block)
}

fn v3k_collection_groups(record: MolBlockRecord) -> Vec<(StereoGroupKind, u32, Vec<usize>)> {
    let MolBlockRecord::Concrete { topology, .. } = record else {
        panic!("ordinary collection record must remain concrete");
    };
    topology
        .stereo_groups
        .iter()
        .map(|group| {
            (
                group.kind(),
                group.id().expect("V3000 collection groups carry an id"),
                group.atoms().iter().map(|atom| atom.index()).collect(),
            )
        })
        .collect()
}

#[test]
fn v3k_collections_uppercase_only_the_first_logical_payload_line() {
    // Pinned parseEnhancedStereo uppercases the first getV3000Line result
    // before entering the loop. Its bottom-of-loop read is not normalized,
    // so the case-sensitive regex recognizes a lowercase/mixed-case first
    // payload but skips the same spellings on subsequent logical lines.
    for first in ["mdlv30/sterel1 atoms=(1 2)", "MdLv30/StErEl1 AtOmS=(1 2)"] {
        let record = v3k_collections_record(
            &["M  V30 1 F 0 0 0 0", "M  V30 2 C 1 0 0 0"],
            &[first],
            true,
        )
        .expect("first logical collection payload is normalized");
        assert_eq!(
            v3k_collection_groups(record),
            vec![(StereoGroupKind::Or, 1, vec![1])]
        );
    }

    for later in ["mdlv30/sterel1 atoms=(1 2)", "MdLv30/StErEl1 AtOmS=(1 2)"] {
        let rows = format!("MDLV30/HILITE ATOMS=(1 2)\n{later}");
        let record = v3k_collections_record(
            &["M  V30 1 F 0 0 0 0", "M  V30 2 C 1 0 0 0"],
            &[&rows],
            true,
        )
        .expect("unrecognized later-case collection row is skipped");
        assert!(v3k_collection_groups(record).is_empty());
    }

    let record = v3k_collections_record(
        &["M  V30 1 F 0 0 0 0", "M  V30 2 C 1 0 0 0"],
        &[concat!(
            "MDLV30/HILITE ATOMS=(1 2)\n",
            "MDLV30/STEREL1 ATOMS=(1 2)"
        )],
        true,
    )
    .expect("uppercase subsequent stereo row is recognized");
    assert_eq!(
        v3k_collection_groups(record),
        vec![(StereoGroupKind::Or, 1, vec![1])]
    );
}

#[test]
fn v3k_collections_apply_first_line_case_rule_after_continuation_and_per_block() {
    // getV3000Line assembles continuation fragments without changing case;
    // parseEnhancedStereo then uppercases only its first assembled payload.
    let first_continued = v3k_collections_raw_record(
        &["M  V30 1 F 0 0 0 0", "M  V30 2 C 1 0 0 0"],
        concat!(
            "M  V30 BEGIN COLLECTION\n",
            "M  V30 mdlv30/sterel1 atoms=(1 -\n",
            "M  V30 2)\n",
            "M  V30 END COLLECTION\n",
        ),
    )
    .expect("continued first payload is normalized after assembly");
    assert_eq!(
        v3k_collection_groups(first_continued),
        vec![(StereoGroupKind::Or, 1, vec![1])]
    );

    let later_continued = v3k_collections_raw_record(
        &["M  V30 1 F 0 0 0 0", "M  V30 2 C 1 0 0 0"],
        concat!(
            "M  V30 BEGIN COLLECTION\n",
            "M  V30 MDLV30/HILITE ATOMS=(1 2)\n",
            "M  V30 mdlv30/sterel1 atoms=(1 -\n",
            "M  V30 2)\n",
            "M  V30 END COLLECTION\n",
        ),
    )
    .expect("continued subsequent lowercase payload is skipped");
    assert!(v3k_collection_groups(later_continued).is_empty());

    let repeated = v3k_collections_raw_record(
        &["M  V30 1 F 0 0 0 0", "M  V30 2 C 1 0 0 0"],
        concat!(
            "M  V30 BEGIN COLLECTION\n",
            "M  V30 MDLV30/STEREL1 ATOMS=(1 1)\n",
            "M  V30 END COLLECTION\n",
            "M  V30 BEGIN COLLECTION\n",
            "M  V30 mdlv30/sterac9 atoms=(1 2)\n",
            "M  V30 END COLLECTION\n",
        ),
    )
    .expect("each collection block normalizes its own first payload");
    assert_eq!(
        v3k_collection_groups(repeated),
        vec![(StereoGroupKind::And, 9, vec![1])]
    );

    let preserves_prior = v3k_collections_raw_record(
        &["M  V30 1 F 0 0 0 0", "M  V30 2 C 1 0 0 0"],
        concat!(
            "M  V30 BEGIN COLLECTION\n",
            "M  V30 MDLV30/STEREL1 ATOMS=(1 1)\n",
            "M  V30 END COLLECTION\n",
            "M  V30 BEGIN COLLECTION\n",
            "M  V30 mdlv30/hilite atoms=(1 2)\n",
            "M  V30 END COLLECTION\n",
        ),
    )
    .expect("empty later collection does not replace existing groups");
    assert_eq!(
        v3k_collection_groups(preserves_prior),
        vec![(StereoGroupKind::Or, 1, vec![0])]
    );
}

#[test]
fn v3k_collections_end_prefix_case_depends_on_logical_position() {
    // The first payload is normalized before the source's three-byte END
    // prefix test, so lowercase END terminates and leaves the CTAB cursor at
    // the following row.
    let first_lower_end = v3k_collections_raw_record(
        &["M  V30 1 C 0 0 0 0"],
        "M  V30 BEGIN COLLECTION\nM  V30 end collection\n",
    )
    .expect("lowercase END is recognized in first payload position");
    assert!(v3k_collection_groups(first_lower_end).is_empty());

    // A subsequent lowercase END is merely an unrecognized collection row.
    // With no later uppercase END, parseEnhancedStereo consumes END CTAB as
    // its terminator and then getV3000Line rejects the following M  END row.
    assert!(matches!(
        v3k_collections_raw_record(
            &["M  V30 1 C 0 0 0 0"],
            concat!(
                "M  V30 BEGIN COLLECTION\n",
                "M  V30 MDLV30/HILITE ATOMS=(1 1)\n",
                "M  V30 end collection\n",
            ),
        ),
        Err(SdfReadError::Parse(message))
            if message.contains("does not start with 'M  V30 '")
    ));

    let later_upper_end = v3k_collections_raw_record(
        &["M  V30 1 C 0 0 0 0"],
        concat!(
            "M  V30 BEGIN COLLECTION\n",
            "M  V30 MDLV30/HILITE ATOMS=(1 1)\n",
            "M  V30 end collection\n",
            "M  V30 END COLLECTION\n",
        ),
    )
    .expect("later uppercase END terminates after skipped lowercase END");
    assert!(v3k_collection_groups(later_upper_end).is_empty());
}

#[test]
fn v3k_collections_use_row_positions_and_source_unsigned_group_ids() {
    // Pinned parseEnhancedStereo calls getAtomWithIdx(index - 1), so these
    // memberships follow atom-table rows and never the nonsequential V3000
    // bookmarks. REL/RAC ids use toUnsigned's initialized-zero overflow.
    let record = v3k_collections_record(
        &[
            "M  V30 42 C 0 0 0 0",
            "M  V30 7 N 1 0 0 0",
            "M  V30 99 O 2 0 0 0",
        ],
        &[concat!(
            "MDLV30/STEABS ATOMS=(2 3 1)\n",
            "MDLV30/STEREL4294967295 ATOMS=(1 2)\n",
            "MDLV30/STERAC4294967296 ATOMS=(1 1)\n",
            "MDLV30/STEREL ATOMS=(1 3)"
        )],
        true,
    )
    .expect("source-shaped collection rows");
    let MolBlockRecord::Concrete { topology, .. } = record else {
        panic!("ordinary collection record must remain concrete");
    };
    let groups = &topology.stereo_groups;
    assert_eq!(groups.len(), 4);
    assert_eq!(groups[0].kind(), StereoGroupKind::Absolute);
    assert_eq!(groups[0].id(), Some(0));
    assert_eq!(
        groups[0]
            .atoms()
            .iter()
            .map(|atom| atom.index())
            .collect::<Vec<_>>(),
        vec![2, 0]
    );
    assert_eq!(
        (groups[1].kind(), groups[1].id()),
        (StereoGroupKind::Or, Some(u32::MAX))
    );
    assert_eq!(groups[1].atoms()[0].index(), 1);
    assert_eq!(
        (groups[2].kind(), groups[2].id()),
        (StereoGroupKind::And, Some(0))
    );
    assert_eq!(groups[2].atoms()[0].index(), 0);
    assert_eq!(
        (groups[3].kind(), groups[3].id()),
        (StereoGroupKind::Or, Some(0))
    );
    assert_eq!(groups[3].atoms()[0].index(), 2);
}

#[test]
fn v3k_collections_distinguish_recognized_errors_from_nonmatching_rows() {
    // regex_match skips malformed/nonmatching collection rows, but an exact
    // stereo-label match with an unknown three-byte tag reaches the source
    // error branch. A plus-prefixed group id is outside the regex and skipped.
    let skipped = v3k_collections_record(
        &["M  V30 1 C 0 0 0 0"],
        &[concat!(
            "MDLV30/HILITE ATOMS=(1 1)\n",
            "MDLV30/STEREL+1 ATOMS=(1 1)\n",
            "MDLV30/STEREL1\tATOMS=(1 1)\n",
            "MDLV30/STEREL1 ATOMS=(1 1) trailing"
        )],
        true,
    )
    .expect("nonmatching collection rows are skipped");
    let MolBlockRecord::Concrete { topology, .. } = skipped else {
        panic!("ordinary collection record must remain concrete");
    };
    assert!(topology.stereo_groups.is_empty());

    assert!(matches!(
        v3k_collections_record(
            &["M  V30 1 C 0 0 0 0"],
            &["MDLV30/STEXYZ1 ATOMS=(1 1)"],
            true,
        ),
        Err(SdfReadError::Parse(message))
            if message.contains("Unrecognized stereogroup type")
    ));
}

#[test]
fn v3k_collections_preserve_source_array_extraction_and_bounds() {
    // The source permits a zero declared count and ignores extra values. At
    // EOF after one initialized value, formatted extraction retains that value
    // for the remaining declared rows. COSMolKit rejects only the source's
    // undefined first-uninitialized extraction and invalid typed row bounds.
    let record = v3k_collections_record(
        &["M  V30 10 C 0 0 0 0", "M  V30 20 N 1 0 0 0"],
        &[concat!(
            "MDLV30/STEABS ATOMS=(0 )\n",
            "MDLV30/STEREL1 ATOMS=(3 2)\n",
            "MDLV30/STERAC2 ATOMS=(1 1 2)"
        )],
        true,
    )
    .expect("source-shaped count handling");
    let MolBlockRecord::Concrete { topology, .. } = record else {
        panic!("ordinary collection record must remain concrete");
    };
    assert!(topology.stereo_groups[0].atoms().is_empty());
    assert_eq!(
        topology.stereo_groups[1]
            .atoms()
            .iter()
            .map(|atom| atom.index())
            .collect::<Vec<_>>(),
        vec![1, 1, 1]
    );
    assert_eq!(
        topology.stereo_groups[2]
            .atoms()
            .iter()
            .map(|atom| atom.index())
            .collect::<Vec<_>>(),
        vec![0]
    );

    for atoms in ["1 ", "1 0", "1 3", "1 4294967295"] {
        let line = format!("MDLV30/STEREL1 ATOMS=({atoms})");
        assert!(
            matches!(
                v3k_collections_record(
                    &["M  V30 10 C 0 0 0 0", "M  V30 20 N 1 0 0 0"],
                    &[&line],
                    true,
                ),
                Err(SdfReadError::Parse(_))
            ),
            "{atoms:?} must fail its typed collection boundary"
        );
    }
}

#[test]
fn v3k_collections_duplicate_absolute_strictness_and_replacement_match_source() {
    let duplicate = concat!("MDLV30/STEABS ATOMS=(1 1)\n", "MDLV30/STEABS ATOMS=(1 2)");
    assert!(matches!(
        v3k_collections_record(
            &["M  V30 1 C 0 0 0 0", "M  V30 2 N 1 0 0 0"],
            &[duplicate],
            true,
        ),
        Err(SdfReadError::Parse(message)) if message.contains("second ABS stereo group")
    ));
    let non_strict = v3k_collections_record(
        &["M  V30 1 C 0 0 0 0", "M  V30 2 N 1 0 0 0"],
        &[duplicate],
        false,
    )
    .expect("non-strict duplicate ABS groups are retained");
    let MolBlockRecord::Concrete { topology, .. } = non_strict else {
        panic!("ordinary collection record must remain concrete");
    };
    assert_eq!(topology.stereo_groups.len(), 2);
    assert_eq!(topology.stereo_groups[0].atoms()[0].index(), 0);
    assert_eq!(topology.stereo_groups[1].atoms()[0].index(), 1);

    let repeated = v3k_collections_record(
        &["M  V30 1 C 0 0 0 0", "M  V30 2 N 1 0 0 0"],
        &[
            "MDLV30/STEREL1 ATOMS=(1 1)",
            "MDLV30/HILITE ATOMS=(1 2)",
            "MDLV30/STERAC9 ATOMS=(1 2)",
        ],
        true,
    )
    .expect("later nonempty collection replaces prior stereo groups");
    let MolBlockRecord::Concrete { topology, .. } = repeated else {
        panic!("ordinary collection record must remain concrete");
    };
    assert_eq!(topology.stereo_groups.len(), 1);
    assert_eq!(topology.stereo_groups[0].kind(), StereoGroupKind::And);
    assert_eq!(topology.stereo_groups[0].id(), Some(9));
    assert_eq!(topology.stereo_groups[0].atoms()[0].index(), 1);
}

#[test]
fn v3k_collections_missing_end_is_a_structured_error() {
    let block = v3000_block(&["M  V30 1 C 0 0 0 0"], 1).replace(
        "M  V30 END CTAB\nM  END\n",
        concat!(
            "M  V30 BEGIN COLLECTION\n",
            "M  V30 MDLV30/STEREL1 ATOMS=(1 1)\n"
        ),
    );
    assert!(matches!(
        read_mol_block_detached(&block),
        Err(SdfReadError::Parse(_))
    ));
}

fn v3k_optional_blocks_record(
    inner_counts: &str,
    trailing: &str,
    strict_parsing: bool,
) -> Result<MolBlockRecord, SdfReadError> {
    let block =
        v3000_with_outer_and_blocks(ZERO_OUTER, inner_counts, &["M  V30 10 C 1 2 0 0"], &[])
            .replace("M  V30 END CTAB", &format!("{trailing}M  V30 END CTAB"));
    cosmolkit_io::read_mol_block_detached_with_params(
        &block,
        cosmolkit_io::MolBlockReadParams {
            strict_parsing,
            ..cosmolkit_io::MolBlockReadParams::default()
        },
    )
}

#[test]
fn v3k_optional_blocks_sgroup_occurrence_matrix_matches_source() {
    // ParseV3000CTAB rejects duplicate SGROUP blocks in both modes, while an
    // undeclared or declared-but-missing block is fatal only in strict mode.
    let group = concat!(
        "M  V30 BEGIN SGROUP\n",
        "M  V30 1 SUP 0 ATOMS=(1 10)\n",
        "M  V30 END SGROUP\n"
    );
    let MolBlockRecord::Concrete { topology, .. } =
        v3k_optional_blocks_record("1 0 1 0 0", group, true).expect("declared SGROUP block")
    else {
        panic!("ordinary SGROUP carrier must remain concrete");
    };
    assert_eq!(topology.substance_groups.len(), 1);
    assert_eq!(topology.substance_groups[0].atoms()[0].index(), 0);

    assert!(matches!(
        v3k_optional_blocks_record("1 0 0 0 0", group, true),
        Err(SdfReadError::Parse(message)) if message.contains("Sgroups NOT expected")
    ));
    let MolBlockRecord::Concrete { topology, .. } =
        v3k_optional_blocks_record("1 0 0 0 0", group, false)
            .expect("non-strict undeclared SGROUP block")
    else {
        panic!("ordinary SGROUP carrier must remain concrete");
    };
    assert_eq!(topology.substance_groups.len(), 1);

    assert!(matches!(
        v3k_optional_blocks_record("1 0 1 0 0", "", true),
        Err(SdfReadError::Parse(message)) if message.contains("BEGIN SGROUP line not found")
    ));
    assert!(v3k_optional_blocks_record("1 0 1 0 0", "", false).is_ok());

    for strict in [true, false] {
        let duplicate = format!("{group}{group}");
        assert!(
            matches!(
                v3k_optional_blocks_record("1 0 1 0 0", &duplicate, strict),
                Err(SdfReadError::Parse(message))
                    if message.contains("BEGIN SGROUP found more than once")
            ),
            "duplicate SGROUP must fail with strict={strict}"
        );
    }
}

#[test]
fn v3k_optional_blocks_obj3d_consumes_declared_rows_without_coordinates() {
    // The pinned parser consumes exactly n3DConstraints logical rows and
    // deliberately ignores their contents instead of constructing conformers.
    let block = concat!(
        "M  V30 BEGIN OBJ3D\n",
        "M  V30 arbitrary constraint one\n",
        "M  V30 99 98 97\n",
        "M  V30 END OBJ3D\n"
    );
    let MolBlockRecord::Concrete { coordinates, .. } =
        v3k_optional_blocks_record("1 0 0 2 0", block, true).expect("declared OBJ3D block")
    else {
        panic!("ordinary OBJ3D carrier must remain concrete");
    };
    assert_eq!(coordinates.conformers_2d.len(), 1);
    assert_eq!(coordinates.conformers_2d[0].coordinates(), &[[1.0, 2.0]]);
    assert!(coordinates.conformers_3d.is_empty());

    assert!(matches!(
        v3k_optional_blocks_record(
            "1 0 0 0 0",
            "M  V30 BEGIN OBJ3D\nM  V30 END OBJ3D\n",
            true,
        ),
        Err(SdfReadError::Parse(message)) if message.contains("OBJ3D")
    ));
    assert!(
        v3k_optional_blocks_record("1 0 0 0 0", "M  V30 BEGIN OBJ3D\nM  V30 END OBJ3D\n", false,)
            .is_ok()
    );

    assert!(matches!(
        v3k_optional_blocks_record("1 0 0 1 0", "", true),
        Err(SdfReadError::Parse(message)) if message.contains("BEGIN OBJ3D line not found")
    ));
    assert!(v3k_optional_blocks_record("1 0 0 1 0", "", false).is_ok());
}

#[test]
fn v3k_optional_blocks_obj3d_duplicate_and_end_policy_match_source() {
    // Give each block its one declared payload row. With a zero declaration,
    // strict source behavior rejects the first block as undeclared before the
    // duplicate-occurrence branch can be observed.
    let declared = concat!(
        "M  V30 BEGIN OBJ3D\n",
        "M  V30 ignored constraint\n",
        "M  V30 END OBJ3D\n"
    );
    for strict in [true, false] {
        let duplicate = format!("{declared}{declared}");
        assert!(
            matches!(
                v3k_optional_blocks_record("1 0 0 1 0", &duplicate, strict),
                Err(SdfReadError::Parse(message))
                    if message.contains("BEGIN OBJ3D found more than once")
            ),
            "duplicate OBJ3D must fail with strict={strict}"
        );
    }

    let malformed = concat!(
        "M  V30 BEGIN OBJ3D\n",
        "M  V30 payload\n",
        "M  V30 NOT THE END\n"
    );
    assert!(matches!(
        v3k_optional_blocks_record("1 0 0 1 0", malformed, true),
        Err(SdfReadError::Parse(message)) if message.contains("END OBJ3D line not found")
    ));
    // Non-strict mode consumes the malformed end row and resumes at END CTAB.
    assert!(v3k_optional_blocks_record("1 0 0 1 0", malformed, false).is_ok());

    let truncated =
        v3000_with_outer_and_blocks(ZERO_OUTER, "1 0 0 1 0", &["M  V30 10 C 1 2 0 0"], &[])
            .replace(
                "M  V30 END CTAB\nM  END\n",
                "M  V30 BEGIN OBJ3D\nM  V30 payload\n",
            );
    assert!(matches!(
        read_mol_block_detached(&truncated),
        Err(SdfReadError::Parse(_))
    ));
}

#[test]
fn v3k_optional_blocks_unknown_skip_is_case_sensitive_after_begin() {
    // ParseV3000CTAB uppercases the BEGIN dispatch row, but the unknown-block
    // loop tests freshly read rows without uppercasing and stops at the first
    // case-sensitive END prefix. The following COLLECTION proves cursor
    // placement after the consumed terminator.
    let trailing = concat!(
        "M  V30 begin vendor\n",
        "M  V30 payload\n",
        "M  V30 end vendor\n",
        "M  V30 still payload\n",
        "M  V30 END ACTUAL\n",
        "M  V30 BEGIN COLLECTION\n",
        "M  V30 MDLV30/STEREL4 ATOMS=(1 1)\n",
        "M  V30 END COLLECTION\n"
    );
    let MolBlockRecord::Concrete { topology, .. } =
        v3k_optional_blocks_record("1 0 0 0 0", trailing, true)
            .expect("unknown block followed by typed collection")
    else {
        panic!("ordinary unknown-block carrier must remain concrete");
    };
    assert_eq!(topology.stereo_groups.len(), 1);
    assert_eq!(topology.stereo_groups[0].id(), Some(4));
    assert_eq!(topology.stereo_groups[0].atoms()[0].index(), 0);

    let truncated =
        v3000_with_outer_and_blocks(ZERO_OUTER, "1 0 0 0 0", &["M  V30 10 C 1 2 0 0"], &[])
            .replace(
                "M  V30 END CTAB\nM  END\n",
                "M  V30 BEGIN VENDOR\nM  V30 payload\n",
            );
    assert!(matches!(
        read_mol_block_detached(&truncated),
        Err(SdfReadError::Parse(_))
    ));
}

#[test]
fn v3k_optional_blocks_malformed_sgroup_terminator_is_not_empty_success() {
    let malformed = concat!(
        "M  V30 BEGIN SGROUP\n",
        "M  V30 1 SUP 0 ATOMS=(1 10)\n",
        "M  V30 NOT END SGROUP\n"
    );
    for strict in [true, false] {
        assert!(
            v3k_optional_blocks_record("1 0 1 0 0", malformed, strict).is_err(),
            "malformed SGROUP termination must fail with strict={strict}"
        );
    }
}

#[test]
fn v3k_ctab_end_strictness_and_cursor_match_source() {
    let valid = v3000_block(&["M  V30 1 C 0 0 0 0"], 1);
    for strict_parsing in [true, false] {
        assert!(
            cosmolkit_io::read_mol_block_detached_with_params(
                &valid,
                cosmolkit_io::MolBlockReadParams {
                    strict_parsing,
                    ..cosmolkit_io::MolBlockReadParams::default()
                },
            )
            .is_ok()
        );
    }

    // ParseV3000CTAB warns about this candidate in non-strict mode, then
    // reads the following raw M END line once instead of reusing it.
    let malformed_ctab_end = valid.replace("M  V30 END CTAB", "M  V30 NOT END CTAB");
    assert!(matches!(
        read_mol_block_detached(&malformed_ctab_end),
        Err(SdfReadError::Parse(message)) if message.contains("END CTAB line not found")
    ));
    assert!(
        cosmolkit_io::read_mol_block_detached_with_params(
            &malformed_ctab_end,
            cosmolkit_io::MolBlockReadParams {
                strict_parsing: false,
                ..cosmolkit_io::MolBlockReadParams::default()
            },
        )
        .is_ok()
    );
}

#[test]
fn v3k_ctab_end_requires_raw_m_end_in_both_modes() {
    let valid = v3000_block(&["M  V30 1 C 0 0 0 0"], 1);
    for malformed in [
        valid.replace("M  END\n", ""),
        valid.replace("M  END", "m  END"),
        valid.replace("M  END", "M END"),
    ] {
        for strict_parsing in [true, false] {
            assert!(
                matches!(
                    cosmolkit_io::read_mol_block_detached_with_params(
                        &malformed,
                        cosmolkit_io::MolBlockReadParams {
                            strict_parsing,
                            ..cosmolkit_io::MolBlockReadParams::default()
                        },
                    ),
                    Err(SdfReadError::Parse(message))
                        if message.contains("M  END line not found")
                ),
                "raw M END must fail with strict={strict_parsing}: {malformed:?}"
            );
        }
    }
}

#[test]
fn v3k_ctab_end_uses_source_prefix_rules_and_ignores_later_lines() {
    let valid = v3000_block(&["M  V30 1 C 0 0 0 0"], 1)
        .replace("M  V30 END CTAB", "M  V30 end ctab suffix")
        .replace("M  END\n", "M  END suffix\ntrailing data\n");
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&valid).expect("source prefix checks accept suffixes")
    else {
        panic!("ordinary terminated record must remain concrete");
    };
    assert_eq!(topology.atoms.len(), 1);
}

#[test]
fn v3k_ctab_end_truncated_or_malformed_envelopes_never_succeed_empty() {
    let valid = v3000_block(&["M  V30 1 C 0 0 0 0"], 1);
    for malformed in [
        valid.replace("M  V30 END CTAB\nM  END\n", ""),
        valid.replace("M  V30 END CTAB\nM  END\n", "M  V30 END CTAB\n"),
        valid.replace("M  V30 END CTAB\nM  END\n", "M  V30 END OTHER\nM  END\n"),
    ] {
        assert!(matches!(
            read_mol_block_detached(&malformed),
            Err(SdfReadError::Parse(_))
        ));
    }
}

#[test]
fn review_collection_formatted_unsigned_extraction() {
    // RDKit 2026.03.1 parseEnhancedStereo: stringstream >> unsigned int.
    // Expected rows verified against the pinned wheel with sanitize=False.
    for (text, expected) in [
        ("1 1x", vec![0]),
        ("2 1+2", vec![0, 1]),
        ("2 1", vec![0, 0]),
        ("2 1 ", vec![0, 0]),
        ("1 -4294967295", vec![0]),
        ("1 +1", vec![0]),
        ("1 01", vec![0]),
        ("2 1\t2", vec![0, 1]),
        ("2 1\u{b}2", vec![0, 1]),
        ("1 1 2", vec![0]),
    ] {
        let line = format!("MDLV30/STEREL1 ATOMS=({text})");
        let MolBlockRecord::Concrete { topology, .. } = review_collection_record(&[&line]).unwrap()
        else {
            panic!("concrete record expected");
        };
        let actual: Vec<_> = topology.stereo_groups[0]
            .atoms()
            .iter()
            .map(|id| id.index())
            .collect();
        assert_eq!(actual, expected, "{text:?}");
    }
    for text in [
        "2 1x",
        "1 0x1",
        "1 4294967296",
        "1 -1",
        "1 ",
        "1 +",
        "1 3",
        "2 1\u{a0}2",
    ] {
        let line = format!("MDLV30/STEREL1 ATOMS=({text})");
        assert!(
            matches!(
                review_collection_record(&[&line]),
                Err(SdfReadError::Parse(_))
            ),
            "{text:?}"
        );
    }
}

#[test]
fn review_collection_only_nonempty_stereo_replaces_previous_groups() {
    // parseEnhancedStereo installs only when !groups.empty(). A recognized
    // group with zero members is still one group and must replace the old one.
    for (second, expected_id, expected_atoms) in [
        ("", 1, vec![0]),
        ("MDLV30/HILITE ATOMS=(1 2)", 1, vec![0]),
        ("MDLV30/STEREL2 ATOMS=(1 2)", 2, vec![1]),
        ("MDLV30/STEREL2 ATOMS=(0 )", 2, vec![]),
    ] {
        let MolBlockRecord::Concrete { topology, .. } =
            review_collection_record(&["MDLV30/STEREL1 ATOMS=(1 1)", second]).unwrap()
        else {
            panic!("concrete record expected");
        };
        assert_eq!(topology.stereo_groups.len(), 1);
        assert_eq!(topology.stereo_groups[0].id(), Some(expected_id));
        let actual: Vec<_> = topology.stereo_groups[0]
            .atoms()
            .iter()
            .map(|id| id.index())
            .collect();
        assert_eq!(actual, expected_atoms, "{second}");
    }
}

#[test]
fn v3k_lines_joins_single_continuation_without_separator() {
    // Physical line 7 is `M  V30 1 C 0 0 -`; its continuation marker drops the
    // trailing dash and concatenates the next payload with no separator.
    let block = v3000_block(&["M  V30 1 C 0 0 -", "M  V30 0 0"], 1);
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("continued atom row")
    else {
        panic!("ordinary atom must produce concrete topology");
    };
    assert_eq!(topology.atoms.len(), 1);
    assert_eq!(topology.atoms[0].atomic_number(), 6);
}

#[test]
fn v3k_lines_joins_multiple_continuations_in_order() {
    // Payload fragments `1 C 0 ` + `0 ` + `0 0` concatenate to a valid row.
    let block = v3000_block(&["M  V30 1 C 0 -", "M  V30 0 -", "M  V30 0 0"], 1);
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("multi-continued atom row")
    else {
        panic!("ordinary atom must produce concrete topology");
    };
    assert_eq!(topology.atoms.len(), 1);
    assert_eq!(topology.atoms[0].atomic_number(), 6);
}

#[test]
fn v3k_lines_rejects_unprefixed_continuation_with_line_number() {
    // The third physical ATOM line lacks the `M  V30 ` prefix.
    let block = v3000_block(&["M  V30 1 C 0 -", "not a v30 line", "M  V30 0 0"], 1);
    match read_mol_block_detached(&block) {
        Err(SdfReadError::Parse(message)) => {
            assert!(
                message.contains("Line 9 does not start with 'M  V30 '"),
                "unexpected diagnostic: {message}"
            );
        }
        other => panic!("expected a structured prefix failure, got {other:?}"),
    }
}

#[test]
fn v3k_lines_reports_eof_on_truncated_continuation() {
    // The final physical line ends with a continuation marker with no next
    // line; the failure must carry the missing physical line number.
    let mut block = v3000_block(&["M  V30 1 C 0 -"], 1);
    // Drop the generated `M  V30 END ATOM\n...` tail so the record ends on the
    // continuation marker.
    if let Some(index) = block.find("M  V30 END ATOM") {
        block.truncate(index);
    }
    match read_mol_block_detached(&block) {
        Err(SdfReadError::Parse(message)) => {
            assert!(
                message.contains("Line 9 does not start with 'M  V30 '"),
                "unexpected diagnostic: {message}"
            );
        }
        other => panic!("expected a structured truncation failure, got {other:?}"),
    }
}

/// Build a V3000 record whose ATOM/BOND blocks are emitted only when the
/// corresponding slice is non-empty, so zero-count block omission is exercised
/// with genuinely absent blocks rather than empty BEGIN/END pairs.
fn v3000_with_outer_and_blocks(
    outer: &str,
    inner_counts: &str,
    atom_lines: &[&str],
    bond_lines: &[&str],
) -> String {
    let mut block = String::new();
    block.push_str("\n  COSMolKit\n\n");
    block.push_str(&format!("{outer} V3000\n"));
    block.push_str("M  V30 BEGIN CTAB\n");
    block.push_str(&format!("M  V30 COUNTS {inner_counts}\n"));
    if !atom_lines.is_empty() {
        block.push_str("M  V30 BEGIN ATOM\n");
        for line in atom_lines {
            block.push_str(line);
            block.push('\n');
        }
        block.push_str("M  V30 END ATOM\n");
    }
    if !bond_lines.is_empty() {
        block.push_str("M  V30 BEGIN BOND\n");
        for line in bond_lines {
            block.push_str(line);
            block.push('\n');
        }
        block.push_str("M  V30 END BOND\n");
    }
    block.push_str("M  V30 END CTAB\nM  END\n");
    block
}

fn v3000_with_outer_and_inner(outer: &str, inner_counts: &str, atom_lines: &[&str]) -> String {
    v3000_with_outer_and_blocks(outer, inner_counts, atom_lines, &[])
}

#[test]
fn v3k_bond_rows_nonsequential_bookmarks_resolve_to_contiguous_typed_rows() {
    // ParseV3000BondBlock uses raw unsigned `from_chars` for all four leading
    // fields, resolves atom bookmarks through the molecule bookmark table, and
    // lets addBond assign canonical row-order ids. Numeric suffixes are not
    // screened because conversion status and end pointers are ignored.
    let block = v3000_with_outer_and_blocks(
        ZERO_OUTER,
        "3 2 0 0 0",
        &[
            "M  V30 10 C 0 0 0 0",
            "M  V30 30 O 1 0 0 0",
            "M  V30 20 N 2 0 0 0",
        ],
        &["M  V30 90tail 1tail 10x 20x", "M  V30 7 2 20 30"],
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("source-shaped bond rows")
    else {
        panic!("ordinary bond rows must remain concrete");
    };
    assert_eq!(topology.bonds.len(), 2);
    assert_eq!(topology.bonds[0].id().index(), 0);
    assert_eq!(topology.bonds[0].begin().index(), 0);
    assert_eq!(topology.bonds[0].end().index(), 2);
    assert_eq!(topology.bonds[0].order(), BondOrder::Single);
    assert_eq!(topology.bonds[0].prop("_MolFileBondType"), Some("1"));
    assert_eq!(topology.bonds[1].id().index(), 1);
    assert_eq!(topology.bonds[1].begin().index(), 2);
    assert_eq!(topology.bonds[1].end().index(), 1);
    assert_eq!(topology.bonds[1].order(), BondOrder::Double);
}

#[test]
fn v3k_bond_rows_missing_endpoints_and_no_conversion_are_structured() {
    let atoms = ["M  V30 10 C 0 0 0 0", "M  V30 20 O 1 0 0 0"];
    for (row, expected_kind, expected_value) in [
        ("M  V30 1 1 99 20", "bond begin atom bookmark", "99"),
        ("M  V30 1 1 10 99", "bond end atom bookmark", "99"),
        // Unsigned std::from_chars does not accept a leading plus; the
        // initialized endpoint remains zero and lookup fails without changing
        // raw atoi semantics elsewhere.
        ("M  V30 1 1 +10 20", "bond begin atom bookmark", "+10"),
        (
            "M  V30 1 1 4294967296 20",
            "bond begin atom bookmark",
            "4294967296",
        ),
    ] {
        let block = v3000_with_outer_and_blocks(ZERO_OUTER, "2 1 0 0 0", &atoms, &[row]);
        assert!(matches!(
            read_mol_block_detached(&block),
            Err(SdfReadError::Field { kind, value, .. })
                if kind == expected_kind && value == expected_value
        ));
    }
}

#[test]
fn v3k_bond_rows_duplicate_external_bookmarks_follow_frozen_unique_policy() {
    // RDKit can retain more than one bond under a bookmark. The detached model
    // intentionally requires an unambiguous external-bookmark map for later
    // SGroup and collection references, so the second row fails structurally.
    let block = v3000_with_outer_and_blocks(
        ZERO_OUTER,
        "3 2 0 0 0",
        &[
            "M  V30 10 C 0 0 0 0",
            "M  V30 20 O 1 0 0 0",
            "M  V30 30 N 2 0 0 0",
        ],
        &["M  V30 44 1 10 20", "M  V30 44 1 20 30"],
    );
    assert!(matches!(
        read_mol_block_detached(&block),
        Err(SdfReadError::Parse(message))
            if message.contains("duplicate V3000 bond index 44")
    ));
}

#[test]
fn v3k_bond_rows_declared_count_rejects_short_extra_and_missing_rows() {
    let atoms = ["M  V30 10 C 0 0 0 0", "M  V30 20 O 1 0 0 0"];

    let short = v3000_with_outer_and_blocks(ZERO_OUTER, "2 1 0 0 0", &atoms, &["M  V30 1 1 10"]);
    assert!(matches!(
        read_mol_block_detached(&short),
        Err(SdfReadError::Parse(message)) if message.contains("is too short")
    ));

    let extra = v3000_with_outer_and_blocks(
        ZERO_OUTER,
        "2 1 0 0 0",
        &atoms,
        &["M  V30 1 1 10 20", "M  V30 2 1 10 20"],
    );
    assert!(matches!(
        read_mol_block_detached(&extra),
        Err(SdfReadError::Parse(message)) if message.contains("END BOND line not found")
    ));

    let missing =
        v3000_with_outer_and_blocks(ZERO_OUTER, "2 2 0 0 0", &atoms, &["M  V30 1 1 10 20"]);
    assert!(matches!(
        read_mol_block_detached(&missing),
        Err(SdfReadError::Parse(message)) if message.contains("is too short")
    ));
}

#[test]
fn v3k_bond_rows_self_duplicate_edges_and_bad_end_marker_fail_structurally() {
    // ROMol::addBond rejects a self bond and an already present undirected
    // edge. The detached reader must report the same invalid graph boundary,
    // never panic while constructing adjacency.
    let atoms = ["M  V30 10 C 0 0 0 0", "M  V30 20 O 1 0 0 0"];
    let self_bond =
        v3000_with_outer_and_blocks(ZERO_OUTER, "2 1 0 0 0", &atoms, &["M  V30 1 1 10 10"]);
    assert!(read_mol_block_detached(&self_bond).is_err());

    let duplicate_edge = v3000_with_outer_and_blocks(
        ZERO_OUTER,
        "2 2 0 0 0",
        &atoms,
        &["M  V30 1 1 10 20", "M  V30 2 1 20 10"],
    );
    assert!(read_mol_block_detached(&duplicate_edge).is_err());

    let bad_end =
        v3000_with_outer_and_blocks(ZERO_OUTER, "2 1 0 0 0", &atoms, &["M  V30 1 1 10 20"])
            .replace("M  V30 END BOND", "M  V30 END BONX");
    assert!(matches!(
        read_mol_block_detached(&bad_end),
        Err(SdfReadError::Parse(message)) if message.contains("END BOND line not found")
    ));
}

#[test]
fn v3k_bond_orders_complete_source_table_preserves_typed_state_and_classification() {
    // Pinned ParseV3000BondBlock has one concrete switch table plus a query
    // default branch. `_MolFileBondType` is unconditional; only 5/6/7 set
    // `_MolFileBondQuery`, while 8 and unknown unsigned types use BondNull.
    let cases = [
        (0, BondOrder::Unspecified, false, None, false),
        (1, BondOrder::Single, false, None, false),
        (2, BondOrder::Double, false, None, false),
        (3, BondOrder::Triple, false, None, false),
        (4, BondOrder::Aromatic, true, None, false),
        (9, BondOrder::Dative, false, None, false),
        (10, BondOrder::Hydrogen, false, None, false),
        (
            5,
            BondOrder::Unspecified,
            false,
            Some(BondQueryPredicate::OrderIn(vec![
                BondOrder::Single,
                BondOrder::Double,
            ])),
            true,
        ),
        (
            6,
            BondOrder::Unspecified,
            false,
            Some(BondQueryPredicate::OrderIn(vec![
                BondOrder::Single,
                BondOrder::Aromatic,
            ])),
            true,
        ),
        (
            7,
            BondOrder::Unspecified,
            false,
            Some(BondQueryPredicate::OrderIn(vec![
                BondOrder::Double,
                BondOrder::Aromatic,
            ])),
            true,
        ),
        (
            8,
            BondOrder::Unspecified,
            false,
            Some(BondQueryPredicate::Any),
            false,
        ),
        (
            11,
            BondOrder::Unspecified,
            false,
            Some(BondQueryPredicate::Any),
            false,
        ),
        (
            u32::MAX,
            BondOrder::Unspecified,
            false,
            Some(BondQueryPredicate::Any),
            false,
        ),
    ];

    for (bond_type, order, aromatic, query_predicate, has_query_prop) in cases {
        let row = format!("M  V30 1 {bond_type} 10 20");
        let block = v3000_with_outer_and_blocks(
            ZERO_OUTER,
            "2 1 0 0 0",
            &["M  V30 10 C 0 0 0 0", "M  V30 20 O 1 0 0 0"],
            &[&row],
        );
        match (read_mol_block_detached(&block), query_predicate) {
            (Ok(MolBlockRecord::Concrete { topology, .. }), None) => {
                let bond = &topology.bonds[0];
                let type_text = bond_type.to_string();
                assert_eq!(bond.order(), order, "type {bond_type}");
                assert_eq!(bond.is_aromatic(), aromatic, "type {bond_type}");
                assert_eq!(
                    bond.prop("_MolFileBondType"),
                    Some(type_text.as_str()),
                    "type {bond_type}"
                );
                assert_eq!(bond.prop("_MolFileBondQuery"), None, "type {bond_type}");
            }
            (Ok(MolBlockRecord::Query(record)), Some(predicate)) => {
                let bond = record.query.bond(0).expect("typed query bond");
                assert_eq!(bond.bond().order(), order, "type {bond_type}");
                assert_eq!(bond.bond().is_aromatic(), aromatic, "type {bond_type}");
                let type_text = bond_type.to_string();
                assert_eq!(
                    bond.bond().prop("_MolFileBondType"),
                    Some(type_text.as_str()),
                    "type {bond_type}"
                );
                assert_eq!(
                    bond.bond().prop("_MolFileBondQuery"),
                    has_query_prop.then_some("1"),
                    "type {bond_type}"
                );
                assert_eq!(
                    bond.predicate(),
                    &QueryNode::predicate(predicate),
                    "type {bond_type}"
                );
            }
            (result, expected) => panic!(
                "type {bond_type} classification mismatch: result={result:?}, expected={expected:?}"
            ),
        }
    }
}

#[test]
fn v3k_bond_cfg_all_values_preserve_direction_stereo_and_stored_value() {
    // Pinned ParseV3000BondBlock applies CFG to the raw molfile bond type:
    // 1/3 wedge or dash every carrier, while CFG=2 distinguishes only raw
    // single and double types. Every accepted value is stored afterward.
    let concrete_cases = [
        (1, 0, BondDirection::None, BondStereo::None),
        (1, 1, BondDirection::BeginWedge, BondStereo::None),
        (1, 2, BondDirection::Unknown, BondStereo::None),
        (1, 3, BondDirection::BeginDash, BondStereo::None),
        (2, 0, BondDirection::None, BondStereo::None),
        (2, 1, BondDirection::BeginWedge, BondStereo::None),
        (2, 2, BondDirection::EitherDouble, BondStereo::Any),
        (2, 3, BondDirection::BeginDash, BondStereo::None),
        (3, 2, BondDirection::None, BondStereo::None),
    ];
    for (bond_type, cfg, direction, stereo) in concrete_cases {
        let row = format!("M  V30 1 {bond_type} 10 20 CFG={cfg}");
        let block = v3000_with_outer_and_blocks(
            ZERO_OUTER,
            "2 1 0 0 0",
            &["M  V30 10 C 0 0 0 0", "M  V30 20 O 1 0 0 0"],
            &[&row],
        );
        let MolBlockRecord::Concrete { topology, .. } =
            read_mol_block_detached(&block).expect("concrete CFG row")
        else {
            panic!("bond type {bond_type} must remain concrete");
        };
        let bond = &topology.bonds[0];
        assert_eq!(bond.direction(), direction, "type {bond_type}, CFG={cfg}");
        assert_eq!(bond.stereo(), stereo, "type {bond_type}, CFG={cfg}");
        let cfg_text = cfg.to_string();
        assert_eq!(
            bond.prop("_MolFileBondCfg"),
            Some(cfg_text.as_str()),
            "type {bond_type}, CFG={cfg}"
        );
    }

    for (cfg, direction) in [
        (1, BondDirection::BeginWedge),
        (2, BondDirection::None),
        (3, BondDirection::BeginDash),
    ] {
        let row = format!("M  V30 1 5 10 20 CFG={cfg}");
        let block = v3000_with_outer_and_blocks(
            ZERO_OUTER,
            "2 1 0 0 0",
            &["M  V30 10 C 0 0 0 0", "M  V30 20 O 1 0 0 0"],
            &[&row],
        );
        let MolBlockRecord::Query(record) = read_mol_block_detached(&block).expect("query CFG row")
        else {
            panic!("bond type 5 must remain a query carrier");
        };
        let bond = record.query.bond(0).expect("typed query bond").bond();
        assert_eq!(bond.direction(), direction, "query CFG={cfg}");
        assert_eq!(bond.stereo(), BondStereo::None, "query CFG={cfg}");
        let cfg_text = cfg.to_string();
        assert_eq!(bond.prop("_MolFileBondCfg"), Some(cfg_text.as_str()));
    }
}

#[test]
fn v3k_bond_cfg_raw_unsigned_prefix_and_boundaries_follow_source() {
    // The source destination is initialized unsigned int and conversion status
    // is ignored: signs/no-conversion/overflow leave zero, while a leading
    // decimal run is retained even when suffix text remains.
    for (value, expected, direction) in [
        ("+1", "0", BondDirection::None),
        ("-1", "0", BondDirection::None),
        ("x", "0", BondDirection::None),
        ("", "0", BondDirection::None),
        ("4294967296", "0", BondDirection::None),
        ("1x", "1", BondDirection::BeginWedge),
    ] {
        let row = format!("M  V30 1 1 10 20 CFG={value}");
        let block = v3000_with_outer_and_blocks(
            ZERO_OUTER,
            "2 1 0 0 0",
            &["M  V30 10 C 0 0 0 0", "M  V30 20 O 1 0 0 0"],
            &[&row],
        );
        let MolBlockRecord::Concrete { topology, .. } =
            read_mol_block_detached(&block).expect("raw CFG conversion")
        else {
            panic!("ordinary CFG carrier must remain concrete");
        };
        let bond = &topology.bonds[0];
        assert_eq!(
            bond.prop("_MolFileBondCfg"),
            Some(expected),
            "CFG={value:?}"
        );
        assert_eq!(bond.direction(), direction, "CFG={value:?}");
    }

    for value in ["4", "4294967295"] {
        let row = format!("M  V30 1 1 10 20 CFG={value}");
        let block = v3000_with_outer_and_blocks(
            ZERO_OUTER,
            "2 1 0 0 0",
            &["M  V30 10 C 0 0 0 0", "M  V30 20 O 1 0 0 0"],
            &[&row],
        );
        assert!(matches!(
            read_mol_block_detached(&block),
            Err(SdfReadError::Parse(message)) if message.contains("bad bond CFG")
        ));
    }

    let malformed = v3000_with_outer_and_blocks(
        ZERO_OUTER,
        "2 1 0 0 0",
        &["M  V30 10 C 0 0 0 0", "M  V30 20 O 1 0 0 0"],
        &["M  V30 1 1 10 20 CFG"],
    );
    assert!(matches!(
        read_mol_block_detached(&malformed),
        Err(SdfReadError::Parse(message)) if message.contains("bad bond property")
    ));
}

#[test]
fn v3k_bond_cfg_wedge_and_dash_propagate_chirality_possible() {
    // calculate3dFlag treats a zero-Z record explicitly marked 3D as 2D only
    // when the parser propagated chiralityPossible. CFG=1/3 do so; CFG=0/2
    // do not, even though single-bond CFG=2 stores UNKNOWN direction.
    for (cfg, expect_2d) in [(0, false), (1, true), (2, false), (3, true)] {
        let row = format!("M  V30 1 1 10 20 CFG={cfg}");
        let mut block = v3000_with_outer_and_blocks(
            ZERO_OUTER,
            "2 1 0 0 0",
            &["M  V30 10 C 0 0 0 0", "M  V30 20 O 1 0 0 0"],
            &[&row],
        );
        block = block.replacen("  COSMolKit", "  COSMolKit         3D", 1);
        let MolBlockRecord::Concrete { coordinates, .. } =
            read_mol_block_detached(&block).expect("marked-3D CFG record")
        else {
            panic!("ordinary CFG carrier must remain concrete");
        };
        assert_eq!(coordinates.conformers_2d.len(), usize::from(expect_2d));
        assert_eq!(coordinates.conformers_3d.len(), usize::from(!expect_2d));
    }
}

#[test]
fn v3k_dimension_header_z_and_cfg_cross_product_matches_source() {
    // Pinned `calculate3dFlag()` gives nonzero Z (strictly greater than the
    // conformer tolerance) precedence, otherwise an exact 3D info-line label
    // remains 3D unless wedge/dash CFG made chirality possible.
    let cases = [
        ("", 0.0_f64, None, CoordinateDimension::TwoD),
        ("2D", 0.0, None, CoordinateDimension::TwoD),
        ("3D", 0.0, None, CoordinateDimension::ThreeD),
        ("", 0.001, None, CoordinateDimension::TwoD),
        ("2D", -0.0005, None, CoordinateDimension::TwoD),
        ("", -0.0, None, CoordinateDimension::TwoD),
        ("3D", 0.001, Some(1), CoordinateDimension::TwoD),
        ("2D", 0.0011, None, CoordinateDimension::ThreeD),
        ("3D", -0.0011, None, CoordinateDimension::ThreeD),
        ("3D", 0.0, Some(1), CoordinateDimension::TwoD),
        ("3D", 0.0, Some(3), CoordinateDimension::TwoD),
        ("3D", 0.0011, Some(1), CoordinateDimension::ThreeD),
    ];

    for (header, z, cfg, expected) in cases {
        let atom1 = format!("M  V30 10 C 1 2 {z} 0");
        let bond = cfg.map(|value| format!("M  V30 1 1 10 20 CFG={value}"));
        let bond_rows = bond.as_deref().into_iter().collect::<Vec<_>>();
        let mut block = v3000_with_outer_and_blocks(
            ZERO_OUTER,
            if cfg.is_some() {
                "2 1 0 0 0"
            } else {
                "2 0 0 0 0"
            },
            &[&atom1, "M  V30 20 O 4 5 0 0"],
            &bond_rows,
        );
        if !header.is_empty() {
            block = block.replacen("  COSMolKit", &format!("  COSMolKit         {header}"), 1);
        }

        let MolBlockRecord::Concrete { coordinates, .. } =
            read_mol_block_detached(&block).expect("dimension cross-product record")
        else {
            panic!("ordinary dimension record must remain concrete");
        };
        assert_eq!(
            coordinates.source_coordinate_dim,
            Some(expected),
            "header={header:?}, z={z}, cfg={cfg:?}"
        );
        match expected {
            CoordinateDimension::TwoD if z.to_bits() == 0 => {
                assert_eq!(coordinates.conformers_2d.len(), 1);
                assert!(coordinates.conformers_3d.is_empty());
                assert_eq!(coordinates.conformers_2d[0].coordinates()[0], [1.0, 2.0]);
            }
            _ => {
                assert!(coordinates.conformers_2d.is_empty());
                assert_eq!(coordinates.conformers_3d.len(), 1);
                let conformer = &coordinates.conformers_3d[0];
                assert_eq!(conformer.is_3d(), expected == CoordinateDimension::ThreeD);
                assert_eq!(
                    conformer.coordinates()[0].map(f64::to_bits),
                    [1.0, 2.0, z].map(f64::to_bits)
                );
            }
        }
    }
}

#[test]
fn v3k_dimension_preserves_empty_conformer_and_atom_row_alignment() {
    // ParseV3000CTAB attaches one conformer even for an empty atom block; for
    // nonempty input its rows retain atom-block order in the chosen storage.
    for (header, expected) in [
        ("", CoordinateDimension::TwoD),
        ("2D", CoordinateDimension::TwoD),
        ("3D", CoordinateDimension::ThreeD),
    ] {
        let mut block = v3000_with_outer_and_blocks(ZERO_OUTER, "0 0 0 0 0", &[], &[]);
        if !header.is_empty() {
            block = block.replacen("  COSMolKit", &format!("  COSMolKit         {header}"), 1);
        }
        let MolBlockRecord::Concrete { coordinates, .. } =
            read_mol_block_detached(&block).expect("empty dimension record")
        else {
            panic!("empty atom table must remain concrete");
        };
        assert_eq!(coordinates.source_coordinate_dim, Some(expected));
        match expected {
            CoordinateDimension::TwoD => {
                assert_eq!(coordinates.conformers_2d.len(), 1);
                assert!(coordinates.conformers_2d[0].coordinates().is_empty());
                assert!(coordinates.conformers_3d.is_empty());
            }
            CoordinateDimension::ThreeD => {
                assert!(coordinates.conformers_2d.is_empty());
                assert_eq!(coordinates.conformers_3d.len(), 1);
                assert!(coordinates.conformers_3d[0].coordinates().is_empty());
            }
        }
    }

    let block = v3000_with_outer_and_blocks(
        ZERO_OUTER,
        "3 0 0 0 0",
        &[
            "M  V30 90 C 1 2 3 0",
            "M  V30 10 N 4 5 -0.0011 0",
            "M  V30 70 O 7 8 9 0",
        ],
        &[],
    );
    let MolBlockRecord::Concrete { coordinates, .. } =
        read_mol_block_detached(&block).expect("3D row-alignment record")
    else {
        panic!("ordinary atom rows must remain concrete");
    };
    assert_eq!(
        coordinates.source_coordinate_dim,
        Some(CoordinateDimension::ThreeD)
    );
    assert_eq!(
        coordinates.conformers_3d[0].coordinates(),
        &[[1.0, 2.0, 3.0], [4.0, 5.0, -0.0011], [7.0, 8.0, 9.0]]
    );
}

#[test]
fn v3k_dimension_query_records_preserve_classified_storage() {
    // ParseV3000CTAB changes the conformer flag, not its XYZ rows. Query
    // conversion must retain both, even when the effective flag is 2D.
    for (header, z, expected) in [
        ("2D", 0.0_f64, CoordinateDimension::TwoD),
        ("2D", 0.001, CoordinateDimension::TwoD),
        ("", -0.0005, CoordinateDimension::TwoD),
        ("", -0.0, CoordinateDimension::TwoD),
        ("", 0.0011, CoordinateDimension::ThreeD),
    ] {
        let atom = format!("M  V30 10 * 1 2 {z} 0");
        let mut block = v3000_with_outer_and_blocks(ZERO_OUTER, "1 0 0 0 0", &[&atom], &[]);
        if !header.is_empty() {
            block = block.replacen("  COSMolKit", &format!("  COSMolKit         {header}"), 1);
        }
        let MolBlockRecord::Query(record) =
            read_mol_block_detached(&block).expect("query dimension record")
        else {
            panic!("wildcard atom must produce query state");
        };
        assert_eq!(record.source_coordinate_dim, Some(expected));
        match expected {
            CoordinateDimension::TwoD if z.to_bits() == 0 => {
                assert_eq!(record.query.coordinates_2d(), Some(&[[1.0, 2.0]][..]));
                assert!(record.query.conformers_3d().is_empty());
            }
            _ => {
                assert!(record.query.coordinates_2d().is_none());
                assert_eq!(record.query.conformers_3d().len(), 1);
                assert_eq!(
                    record.query.conformers_3d()[0].is_3d(),
                    expected == CoordinateDimension::ThreeD
                );
                assert_eq!(
                    record.query.conformers_3d()[0].coordinates()[0].map(f64::to_bits),
                    [1.0, 2.0, z].map(f64::to_bits)
                );
            }
        }
    }
}

#[test]
fn v3k_linknodes_zero_and_empty_rows_do_not_install_a_property() {
    // Pinned ParseV3000CTAB ignores LINKNODE text with no payload beyond the
    // separating byte and leaves the molecule property absent.
    for rows in ["", "M  V30 LINKNODE\n", "M  V30 LINKNODE \n"] {
        let block =
            v3000_with_outer_and_blocks(ZERO_OUTER, "1 0 0 0 0", &["M  V30 10 C 0 0 0 0"], &[])
                .replace("M  V30 END CTAB", &format!("{rows}M  V30 END CTAB"));
        let MolBlockRecord::Concrete { properties, .. } =
            read_mol_block_detached(&block).expect("empty LINKNODE record")
        else {
            panic!("ordinary atom must remain concrete");
        };
        assert_eq!(properties.prop("_MolFileLinkNodes"), None);
    }
}

#[test]
fn v3k_linknodes_uppercase_and_accumulate_in_source_order() {
    // Pinned ParseV3000CTAB uppercases the complete logical row before taking
    // the payload and appends successive payloads with a single `|`.
    let block = v3000_with_outer_and_blocks(ZERO_OUTER, "1 0 0 0 0", &["M  V30 10 C 0 0 0 0"], &[])
        .replace(
            "M  V30 END CTAB",
            concat!(
                "M  V30 linknode first mixedCase\n",
                "M  V30 LINKNODE second lower\n",
                "M  V30 END CTAB",
            ),
        );
    let MolBlockRecord::Concrete { properties, .. } =
        read_mol_block_detached(&block).expect("ordered LINKNODE rows")
    else {
        panic!("ordinary atom must remain concrete");
    };
    assert_eq!(
        properties.prop("_MolFileLinkNodes"),
        Some("FIRST MIXEDCASE|SECOND LOWER")
    );
}

#[test]
fn v3k_linknodes_continuations_preserve_other_detached_properties() {
    // getV3000Line joins continued physical payloads without inserting text;
    // the LINKNODE loop uppercases that assembled row and does not disturb
    // title, info-line, comment, or CTAB chiral-flag properties.
    let mut block =
        v3000_with_outer_and_blocks(ZERO_OUTER, "1 0 0 0 1", &["M  V30 10 C 0 0 0 0"], &[]);
    block = block.replacen("\n  COSMolKit\n\n", "named\n  Generator 2D\ncomment\n", 1);
    block = block.replace(
        "M  V30 END CTAB",
        concat!(
            "M  V30 linknode alpha -\n",
            "M  V30 beta\n",
            "M  V30 END CTAB",
        ),
    );
    let MolBlockRecord::Concrete { properties, .. } =
        read_mol_block_detached(&block).expect("continued LINKNODE row")
    else {
        panic!("ordinary atom must remain concrete");
    };
    assert_eq!(properties.prop("_MolFileLinkNodes"), Some("ALPHA BETA"));
    assert_eq!(properties.name(), Some("named"));
    assert_eq!(properties.prop("_MolFileInfo"), Some("  Generator 2D"));
    assert_eq!(properties.prop("_MolFileComments"), Some("comment"));
    assert_eq!(properties.prop("_MolFileChiralFlag"), Some("1"));
}

#[test]
fn v3k_bond_topo_concrete_conversion_preserves_order_and_negation() {
    for (bond_type, order, value, in_ring) in [
        (0, BondOrder::Unspecified, "1", true),
        (1, BondOrder::Single, "1", true),
        (2, BondOrder::Double, "2", false),
    ] {
        let row = format!("M  V30 1 {bond_type} 10 20 TOPO={value}");
        let block = v3000_with_outer_and_blocks(
            ZERO_OUTER,
            "2 1 0 0 0",
            &["M  V30 10 C 0 0 0 0", "M  V30 20 O 1 0 0 0"],
            &[&row],
        );
        let MolBlockRecord::Query(record) = read_mol_block_detached(&block).unwrap() else {
            panic!("nonzero TOPO must convert a concrete bond to query state");
        };
        assert_eq!(
            record.query.bond(0).unwrap().predicate(),
            &QueryNode::and(vec![
                QueryNode::predicate(BondQueryPredicate::Order(order)),
                QueryNode::predicate(BondQueryPredicate::IsInRing(in_ring)),
            ])
        );
    }
}

#[test]
fn v3k_bond_topo_existing_queries_compose_with_source_null_algebra() {
    for (bond_type, value, expected) in [
        (
            5,
            "1",
            QueryNode::and(vec![
                QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                    BondOrder::Single,
                    BondOrder::Double,
                ])),
                QueryNode::predicate(BondQueryPredicate::IsInRing(true)),
            ]),
        ),
        (
            8,
            "2",
            QueryNode::predicate(BondQueryPredicate::IsInRing(false)),
        ),
    ] {
        let row = format!("M  V30 1 {bond_type} 10 20 TOPO={value}");
        let block = v3000_with_outer_and_blocks(
            ZERO_OUTER,
            "2 1 0 0 0",
            &["M  V30 10 C 0 0 0 0", "M  V30 20 O 1 0 0 0"],
            &[&row],
        );
        let MolBlockRecord::Query(record) = read_mol_block_detached(&block).unwrap() else {
            panic!("query bond must remain query state");
        };
        assert_eq!(record.query.bond(0).unwrap().predicate(), &expected);
    }
}

#[test]
fn v3k_bond_topo_uses_exact_literal_zero_one_two_only() {
    let literal_zero = v3000_with_outer_and_blocks(
        ZERO_OUTER,
        "2 1 0 0 0",
        &["M  V30 10 C 0 0 0 0", "M  V30 20 O 1 0 0 0"],
        &["M  V30 1 1 10 20 TOPO=0"],
    );
    assert!(matches!(
        read_mol_block_detached(&literal_zero),
        Ok(MolBlockRecord::Concrete { .. })
    ));

    let query_zero = v3000_with_outer_and_blocks(
        ZERO_OUTER,
        "2 1 0 0 0",
        &["M  V30 10 C 0 0 0 0", "M  V30 20 O 1 0 0 0"],
        &["M  V30 1 5 10 20 TOPO=0"],
    );
    let MolBlockRecord::Query(record) = read_mol_block_detached(&query_zero).unwrap() else {
        panic!("literal zero must preserve the existing query");
    };
    assert_eq!(
        record.query.bond(0).unwrap().predicate(),
        &QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
            BondOrder::Single,
            BondOrder::Double
        ]))
    );

    for value in ["00", "01", "+1", "1x", "", "3"] {
        let row = format!("M  V30 1 1 10 20 TOPO={value}");
        let block = v3000_with_outer_and_blocks(
            ZERO_OUTER,
            "2 1 0 0 0",
            &["M  V30 10 C 0 0 0 0", "M  V30 20 O 1 0 0 0"],
            &[&row],
        );
        assert!(matches!(
            read_mol_block_detached(&block),
            Err(SdfReadError::Parse(message)) if message.contains("bad bond TOPO")
        ));
    }
}

#[test]
fn v3k_bond_props_rxctr_and_raw_values_follow_source_contracts() {
    // Pinned ParseV3000BondBlock sends RXCTR through FileParserUtils::toInt,
    // but installs the tokenized STBOX, ENDPTS, and ATTACH strings verbatim.
    let block = v3000_with_outer_and_blocks(
        ZERO_OUTER,
        "2 1 0 0 0",
        &["M  V30 10 C 0 0 0 0", "M  V30 20 O 1 0 0 0"],
        &["M  V30 1 1 10 20 RXCTR=-7 STBOX=raw ENDPTS=(2 20 10) ATTACH=ANY"],
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("remaining bond properties")
    else {
        panic!("ordinary bond properties must remain concrete");
    };
    let bond = &topology.bonds[0];
    assert_eq!(bond.prop("molReactStatus"), Some("-7"));
    assert_eq!(bond.prop("molStereoCare"), Some("raw"));
    assert_eq!(bond.prop("_MolFileBondEndPts"), Some("(2 20 10)"));
    assert_eq!(bond.prop("_MolFileBondAttach"), Some("ANY"));

    for (value, expected) in [("7", "7"), ("+7", "0"), ("", "0"), ("2147483648", "0")] {
        let row = format!("M  V30 1 1 10 20 RXCTR={value}");
        let block = v3000_with_outer_and_blocks(
            ZERO_OUTER,
            "2 1 0 0 0",
            &["M  V30 10 C 0 0 0 0", "M  V30 20 O 1 0 0 0"],
            &[&row],
        );
        let MolBlockRecord::Concrete { topology, .. } =
            read_mol_block_detached(&block).expect("source-shaped RXCTR")
        else {
            panic!("RXCTR must not create query state");
        };
        assert_eq!(topology.bonds[0].prop("molReactStatus"), Some(expected));
    }

    let invalid = v3000_with_outer_and_blocks(
        ZERO_OUTER,
        "2 1 0 0 0",
        &["M  V30 10 C 0 0 0 0", "M  V30 20 O 1 0 0 0"],
        &["M  V30 1 1 10 20 RXCTR=7x"],
    );
    assert!(matches!(
        read_mol_block_detached(&invalid),
        Err(SdfReadError::Field {
            kind: "V3000 bond reaction status",
            value,
            ..
        }) if value == "7x"
    ));
}

#[test]
fn v3k_bond_props_endpoint_stereo_care_requires_equal_pair_without_override() {
    // Source propagation happens only after property parsing and only when
    // both endpoint atom properties exist, compare equal, and no bond STBOX
    // (including an explicit empty value) was installed.
    for (atom_one, atom_two, bond_props, expected) in [
        ("STBOX=1", "STBOX=1", "", Some("1")),
        ("STBOX=1", "STBOX=2", "", None),
        ("STBOX=1", "", "", None),
        ("STBOX=1", "STBOX=1", "STBOX=explicit", Some("explicit")),
        ("STBOX=1", "STBOX=1", "STBOX=", Some("")),
        ("STBOX=0", "STBOX=0", "", None),
    ] {
        let atom_one = format!("M  V30 10 C 0 0 0 0 {atom_one}");
        let atom_two = format!("M  V30 20 O 1 0 0 0 {atom_two}");
        let bond = format!("M  V30 1 1 10 20 {bond_props}");
        let block =
            v3000_with_outer_and_blocks(ZERO_OUTER, "2 1 0 0 0", &[&atom_one, &atom_two], &[&bond]);
        let MolBlockRecord::Concrete { topology, .. } =
            read_mol_block_detached(&block).expect("endpoint stereo care")
        else {
            panic!("stereo-care properties must not create query state");
        };
        assert_eq!(
            topology.bonds[0].prop("molStereoCare"),
            expected,
            "atoms=({atom_one:?}, {atom_two:?}), bond={bond:?}"
        );
    }
}

#[test]
fn v3k_bond_props_malformed_assignments_fail_structurally() {
    for token in ["RXCTR", "STBOX=a=b", "ENDPTS=(2 10 20)=x", "ATTACH"] {
        let row = format!("M  V30 1 1 10 20 {token}");
        let block = v3000_with_outer_and_blocks(
            ZERO_OUTER,
            "2 1 0 0 0",
            &["M  V30 10 C 0 0 0 0", "M  V30 20 O 1 0 0 0"],
            &[&row],
        );
        assert!(matches!(
            read_mol_block_detached(&block),
            Err(SdfReadError::Parse(message))
                if message.contains("bad bond property") && message.contains(token)
        ));
    }
}

#[test]
fn v3k_counts_rejects_nonzero_outer_counts_in_strict_mode() {
    let block = v3000_with_outer_and_inner("  1  1  0  0  0  0  0  0  0  0999", "0 0 0 0 0", &[]);
    match read_mol_block_detached(&block) {
        Err(SdfReadError::Parse(message)) => assert!(
            message.contains("should have 0s"),
            "unexpected diagnostic: {message}"
        ),
        other => panic!("expected structured outer-count failure, got {other:?}"),
    }
}

#[test]
fn v3k_counts_allows_nonzero_outer_counts_in_non_strict_mode() {
    let block = v3000_with_outer_and_inner("  1  1  0  0  0  0  0  0  0  0999", "0 0 0 0 0", &[]);
    let params = cosmolkit_io::MolBlockReadParams {
        strict_parsing: false,
        ..cosmolkit_io::MolBlockReadParams::default()
    };
    let result = cosmolkit_io::read_mol_block_detached_with_params(&block, params);
    assert!(
        result.is_ok(),
        "non-strict outer counts must parse: {result:?}"
    );
}

#[test]
fn v3k_counts_requires_two_inner_count_fields() {
    let missing = v3000_with_outer_and_inner("  0  0  0  0  0  0  0  0  0  0999", "0", &[]);
    assert!(
        matches!(read_mol_block_detached(&missing), Err(SdfReadError::Counts)),
        "a single COUNTS field must be a structured Counts error"
    );
    let two_fields = v3000_with_outer_and_inner("  0  0  0  0  0  0  0  0  0  0999", "0 0", &[]);
    assert!(read_mol_block_detached(&two_fields).is_ok());
}

#[test]
fn v3k_counts_omits_zero_count_blocks_without_panic() {
    let block = v3000_with_outer_and_inner("  0  0  0  0  0  0  0  0  0  0999", "0 0 0 0 0", &[]);
    let MolBlockRecord::Concrete {
        topology,
        coordinates,
        properties,
    } = read_mol_block_detached(&block).expect("zero-count record")
    else {
        panic!("zero-count record must be concrete");
    };
    assert!(topology.atoms.is_empty());
    assert!(topology.bonds.is_empty());
    // RDKit still constructs the declared `Conformer(0)`; every coordinate row
    // (across both dimensions) must be absent, matching the atom table.
    assert!(
        coordinates
            .conformers_2d
            .iter()
            .all(|conformer| conformer.coordinates().is_empty())
    );
    assert!(coordinates.conformers_3d.is_empty());
    assert_eq!(properties.prop("_MolFileChiralFlag"), Some("0"));
}

#[test]
fn v3k_counts_rejects_malformed_counts_without_panic() {
    let block = v3000_with_outer_and_inner("  0  0  0  0  0  0  0  0  0  0999", "abc 0 0 0 0", &[]);
    assert!(
        matches!(
            read_mol_block_detached(&block),
            Err(SdfReadError::Field { .. })
        ),
        "malformed count text must be a structured field error"
    );
}

#[test]
fn v3k_counts_leading_plus_is_zero_atoms() {
    // Pinned `FileParserUtils::toUnsigned` passes the text to
    // `std::from_chars`, which does not recognize a leading `+`; RDKit therefore
    // reads `+1` as 0 atoms, so no ATOM block is required.
    let block = v3000_with_outer_and_inner("  0  0  0  0  0  0  0  0  0  0999", "+1 0 0 0 0", &[]);
    let MolBlockRecord::Concrete {
        topology,
        properties,
        ..
    } = read_mol_block_detached(&block).expect("leading plus must parse as zero atoms")
    else {
        panic!("zero-atom record must be concrete");
    };
    assert!(topology.atoms.is_empty());
    assert_eq!(properties.prop("_MolFileChiralFlag"), Some("0"));
}

#[test]
fn v3k_counts_unsigned_overflow_is_zero_atoms() {
    // `std::from_chars` reports out-of-range without modifying its initialized
    // target, so `2^32` reads as 0 atoms and must not allocate.
    let block = v3000_with_outer_and_inner(
        "  0  0  0  0  0  0  0  0  0  0999",
        "4294967296 0 0 0 0",
        &[],
    );
    let MolBlockRecord::Concrete {
        topology,
        properties,
        ..
    } = read_mol_block_detached(&block).expect("overflow must parse as zero atoms")
    else {
        panic!("zero-atom record must be concrete");
    };
    assert!(topology.atoms.is_empty());
    assert_eq!(properties.prop("_MolFileChiralFlag"), Some("0"));
}

#[test]
fn v3k_counts_chiral_flag_unsigned_boundaries() {
    // The maximum `unsigned int` is accepted and preserved. It is applied to
    // the chiral flag only: using it as an atom count would be an allocation
    // probe, which this suite intentionally avoids.
    let max = v3000_with_outer_and_inner(
        "  0  0  0  0  0  0  0  0  0  0999",
        "0 0 0 0 4294967295",
        &[],
    );
    let MolBlockRecord::Concrete { properties, .. } =
        read_mol_block_detached(&max).expect("max chiral flag")
    else {
        panic!("zero-atom record must be concrete");
    };
    assert_eq!(properties.prop("_MolFileChiralFlag"), Some("4294967295"));

    let overflow = v3000_with_outer_and_inner(
        "  0  0  0  0  0  0  0  0  0  0999",
        "0 0 0 0 4294967296",
        &[],
    );
    let MolBlockRecord::Concrete { properties, .. } =
        read_mol_block_detached(&overflow).expect("overflow chiral flag")
    else {
        panic!("zero-atom record must be concrete");
    };
    assert_eq!(properties.prop("_MolFileChiralFlag"), Some("0"));
}

#[test]
fn v3k_counts_optional_fields_absent_and_present() {
    let outer = "  0  0  0  0  0  0  0  0  0  0999";
    let atom = "M  V30 1 C 0.0 0.0 0.0 0";

    let absent = v3000_with_outer_and_inner(outer, "1 0", &[atom]);
    let MolBlockRecord::Concrete {
        topology,
        coordinates,
        properties,
    } = read_mol_block_detached(&absent).expect("absent optional COUNTS fields")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(topology.atoms.len(), 1);
    assert_eq!(coordinates.conformers_2d.len(), 1);
    assert_eq!(coordinates.conformers_2d[0].coordinates().len(), 1);
    assert_eq!(properties.prop("_MolFileChiralFlag"), Some("0"));

    let present = v3000_with_outer_and_inner(outer, "1 0 0 0 1", &[atom]);
    let MolBlockRecord::Concrete { properties, .. } =
        read_mol_block_detached(&present).expect("present optional COUNTS fields")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(properties.prop("_MolFileChiralFlag"), Some("1"));
}

#[test]
fn v3k_counts_zero_bond_record_omits_bond_block() {
    // Declared bond count 0 with no BEGIN BOND block must remain valid and must
    // not fabricate bonds.
    let block = v3000_with_outer_and_blocks(
        "  0  0  0  0  0  0  0  0  0  0999",
        "1 0 0 0 0",
        &["M  V30 1 C 0.0 0.0 0.0 0"],
        &[],
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("zero-bond record")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(topology.atoms.len(), 1);
    assert!(topology.bonds.is_empty());
}

#[test]
fn v3k_counts_lowercase_keyword_is_accepted() {
    // The counts line is uppercased before `COUNTS ` is matched, so a lowercase
    // keyword in the source block is valid.
    let block = concat!(
        "\n  COSMolKit\n\n",
        "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
        "M  V30 BEGIN CTAB\n",
        "M  V30 counts 1 0 0 0 0\n",
        "M  V30 BEGIN ATOM\n",
        "M  V30 1 C 0.0 0.0 0.0 0\n",
        "M  V30 END ATOM\n",
        "M  V30 END CTAB\n",
        "M  END\n",
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(block).expect("lowercase counts keyword")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(topology.atoms.len(), 1);
}

#[test]
fn v3k_counts_rejects_invalid_optional_fields() {
    let outer = "  0  0  0  0  0  0  0  0  0  0999";
    for invalid in ["0 0 x 0 0", "0 0 0 x 0", "0 0 0 0 x"] {
        let block = v3000_with_outer_and_inner(outer, invalid, &[]);
        assert!(
            matches!(
                read_mol_block_detached(&block),
                Err(SdfReadError::Field { .. })
            ),
            "invalid optional field {invalid:?} must be a structured field error"
        );
    }
}

fn v3000_block_with_line_ending(line_ending: &str) -> String {
    let lines = [
        "",
        "  COSMolKit",
        "",
        "  0  0  0  0  0  0  0  0  0  0999 V3000",
        "M  V30 BEGIN CTAB",
        "M  V30 COUNTS 1 0 0 0 0",
        "M  V30 BEGIN ATOM",
        "M  V30 1 C 1-",
        "M  V30 0 0 0 0",
        "M  V30 END ATOM",
        "M  V30 END CTAB",
        "M  END",
    ];
    let mut block = lines.join(line_ending);
    block.push_str(line_ending);
    block
}

#[test]
fn v3k_lines_joins_without_inserting_whitespace() {
    // `1-` + `0 0 0 0` concatenates to `10 0 0 0`, so the x coordinate is 10.
    // If a separator were inserted the row would read `1 0 0 0 0` (x = 1).
    let block = v3000_block_with_line_ending("\n");
    let MolBlockRecord::Concrete { coordinates, .. } =
        read_mol_block_detached(&block).expect("continued split token")
    else {
        panic!("ordinary atom must produce concrete topology");
    };
    assert_eq!(coordinates.conformers_2d.len(), 1);
    assert_eq!(coordinates.conformers_2d[0].coordinates()[0][0], 10.0);
}

#[test]
fn v3k_lines_handles_crlf_line_endings() {
    let block = v3000_block_with_line_ending("\r\n");
    let MolBlockRecord::Concrete {
        topology,
        coordinates,
        ..
    } = read_mol_block_detached(&block).expect("CRLF record")
    else {
        panic!("ordinary atom must produce concrete topology");
    };
    assert_eq!(topology.atoms.len(), 1);
    assert_eq!(coordinates.conformers_2d.len(), 1);
    assert_eq!(coordinates.conformers_2d[0].coordinates()[0][0], 10.0);
}

const ZERO_OUTER: &str = "  0  0  0  0  0  0  0  0  0  0999";

#[test]
fn v3k_atom_rows_nonsequential_bookmarks_resolve_by_bookmark() {
    // Bookmarks 10/20/30 are not row positions; the bond must resolve them
    // through the bookmark map, not by using the parsed value as a row index.
    let block = v3000_with_outer_and_blocks(
        ZERO_OUTER,
        "3 1 0 0 0",
        &[
            "M  V30 10 C 0.0 0.0 0.0 0",
            "M  V30 20 O 1.0 0.0 0.0 0",
            "M  V30 30 N 2.0 0.0 0.0 0",
        ],
        &["M  V30 1 1 10 20"],
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("nonsequential bookmarks")
    else {
        panic!("atom/bond record must be concrete");
    };
    assert_eq!(topology.atoms.len(), 3);
    assert_eq!(topology.bonds.len(), 1);
    assert_eq!(topology.bonds[0].begin().index(), 0);
    assert_eq!(topology.bonds[0].end().index(), 1);
}

#[test]
fn v3k_atom_rows_duplicate_bookmark_keeps_first_atom() {
    // `setAtomBookmark` appends and `getAtomWithBookmark` returns the first
    // inserted atom, so the duplicate bookmark 1 resolves to row 0.
    let block = v3000_with_outer_and_blocks(
        ZERO_OUTER,
        "3 1 0 0 0",
        &[
            "M  V30 1 C 0.0 0.0 0.0 0",
            "M  V30 1 O 1.0 0.0 0.0 0",
            "M  V30 2 N 2.0 0.0 0.0 0",
        ],
        &["M  V30 1 1 1 2"],
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("duplicate bookmark keeps first")
    else {
        panic!("atom/bond record must be concrete");
    };
    assert_eq!(topology.atoms.len(), 3);
    assert_eq!(topology.bonds[0].begin().index(), 0);
    assert_eq!(topology.bonds[0].end().index(), 2);
    assert_eq!(topology.atoms[0].atomic_number(), 6);
}

#[test]
fn v3k_atom_rows_zero_and_malformed_bookmarks_are_zero() {
    // `std::from_chars` leaves the initialized target at 0 for `abc`; the first
    // bookmark-0 association (row 0) is the one a bond resolves to.
    let block = v3000_with_outer_and_blocks(
        ZERO_OUTER,
        "3 1 0 0 0",
        &[
            "M  V30 0 C 0.0 0.0 0.0 0",
            "M  V30 abc O 1.0 0.0 0.0 0",
            "M  V30 2 N 2.0 0.0 0.0 0",
        ],
        &["M  V30 1 1 0 2"],
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("malformed bookmark reads as zero")
    else {
        panic!("atom/bond record must be concrete");
    };
    assert_eq!(topology.bonds[0].begin().index(), 0);
    assert_eq!(topology.bonds[0].end().index(), 2);
}

#[test]
fn v3k_atom_rows_prefix_and_signed_bookmarks_follow_from_chars() {
    // `1x` parses as 1 (maximal decimal-digit prefix); `+1` is not a recognized
    // sign and parses as 0, so bookmark 1 resolves to row 0 and bookmark 2 to
    // row 2, leaving the `+1` row unreferenced at index 1.
    let block = v3000_with_outer_and_blocks(
        ZERO_OUTER,
        "3 1 0 0 0",
        &[
            "M  V30 1x C 0.0 0.0 0.0 0",
            "M  V30 +1 O 1.0 0.0 0.0 0",
            "M  V30 2 N 2.0 0.0 0.0 0",
        ],
        &["M  V30 1 1 1 2"],
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("prefix/sign bookmark")
    else {
        panic!("atom/bond record must be concrete");
    };
    assert_eq!(topology.bonds[0].begin().index(), 0);
    assert_eq!(topology.bonds[0].end().index(), 2);
}

#[test]
fn v3k_atom_rows_short_row_is_bad_atom_line() {
    let block =
        v3000_with_outer_and_blocks(ZERO_OUTER, "1 0 0 0 0", &["M  V30 1 C 0.0 0.0 0.0"], &[]);
    match read_mol_block_detached(&block) {
        Err(SdfReadError::Parse(message)) => assert!(
            message.contains("Bad atom line"),
            "unexpected diagnostic: {message}"
        ),
        other => panic!("expected a short-row failure, got {other:?}"),
    }
}

#[test]
fn v3k_atom_rows_early_end_atom_is_bad_atom_line() {
    // Two atoms declared but only one row before `END ATOM`; the marker line is
    // consumed as the second atom row and fails the required-field check.
    let block =
        v3000_with_outer_and_blocks(ZERO_OUTER, "2 0 0 0 0", &["M  V30 1 C 0.0 0.0 0.0 0"], &[]);
    match read_mol_block_detached(&block) {
        Err(SdfReadError::Parse(message)) => assert!(
            message.contains("Bad atom line"),
            "unexpected diagnostic: {message}"
        ),
        other => panic!("expected an early-END-ATOM failure, got {other:?}"),
    }
}

#[test]
fn v3k_atom_rows_late_end_atom_is_reported() {
    // One atom declared but two rows present; the extra row is read by the
    // `END ATOM` check and reported.
    let block = v3000_with_outer_and_blocks(
        ZERO_OUTER,
        "1 0 0 0 0",
        &["M  V30 1 C 0.0 0.0 0.0 0", "M  V30 2 O 1.0 0.0 0.0 0"],
        &[],
    );
    match read_mol_block_detached(&block) {
        Err(SdfReadError::Parse(message)) => assert!(
            message.contains("END ATOM"),
            "unexpected diagnostic: {message}"
        ),
        other => panic!("expected a late-END-ATOM failure, got {other:?}"),
    }
}

#[test]
fn v3k_atom_rows_missing_bookmark_reference_is_structured_error() {
    // RDKit's `getAtomWithBookmark` is a debug-only precondition; release
    // behavior is undefined, so the detached reader uses a structured error.
    let block = v3000_with_outer_and_blocks(
        ZERO_OUTER,
        "2 1 0 0 0",
        &["M  V30 1 C 0.0 0.0 0.0 0", "M  V30 2 O 1.0 0.0 0.0 0"],
        &["M  V30 1 1 1 99"],
    );
    match read_mol_block_detached(&block) {
        Err(SdfReadError::Field { kind, .. }) => {
            assert_eq!(kind, "bond end atom bookmark");
        }
        other => panic!("expected a missing-reference failure, got {other:?}"),
    }
}

#[test]
fn v3k_atom_numbers_prefix_suffix_and_exponent() {
    // A non-zero z selects 3D storage, so all three axes are observable.
    let block = v3000_with_outer_and_inner(
        ZERO_OUTER,
        "2 0 0 0 0",
        &["M  V30 1 C 1.5xyz 2.5 0 0", "M  V30 2 O 1e2 2E-1 3.0e+1 0"],
    );
    let MolBlockRecord::Concrete { coordinates, .. } =
        read_mol_block_detached(&block).expect("suffix/exponent coordinates")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(coordinates.conformers_3d.len(), 1);
    let rows = coordinates.conformers_3d[0].coordinates();
    assert_eq!(rows[0], [1.5, 2.5, 0.0]);
    assert_eq!(rows[1], [100.0, 0.2, 30.0]);
}

#[test]
fn v3k_atom_numbers_leading_signs_and_invalid_text() {
    let block = v3000_with_outer_and_inner(
        ZERO_OUTER,
        "3 0 0 0 0",
        &[
            "M  V30 1 C +1.5 -2.5 0 0",
            "M  V30 2 O .5 0 0 0",
            "M  V30 3 N abc def 0 0",
        ],
    );
    let MolBlockRecord::Concrete { coordinates, .. } =
        read_mol_block_detached(&block).expect("signs and invalid text")
    else {
        panic!("atom record must be concrete");
    };
    let rows = coordinates.conformers_2d[0].coordinates();
    assert_eq!(rows[0], [1.5, -2.5]);
    assert_eq!(rows[1], [0.5, 0.0]);
    assert_eq!(rows[2], [0.0, 0.0]);
}

#[test]
fn v3k_atom_numbers_nonfinite_coordinates_are_rejected() {
    for value in ["inf", "nan", "1e999"] {
        let block = v3000_with_outer_and_inner(
            ZERO_OUTER,
            "1 0 0 0 0",
            &[&format!("M  V30 1 C {value} 0 0 0")],
        );
        assert!(
            matches!(
                read_mol_block_detached(&block),
                Err(SdfReadError::Coordinates(_))
            ),
            "non-finite coordinate {value:?} must fail the finite-coordinate boundary"
        );
    }
}

#[test]
fn v3k_atom_numbers_atom_map_positive_only_with_prefix_text() {
    let block = v3000_with_outer_and_inner(
        ZERO_OUTER,
        "5 0 0 0 0",
        &[
            "M  V30 1 C 0 0 0 5x",
            "M  V30 2 O 0 0 0 5.9",
            "M  V30 3 N 0 0 0 -3",
            "M  V30 4 F 0 0 0 0",
            "M  V30 5 Cl 0 0 0 +7",
        ],
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("atom map conversions")
    else {
        panic!("atom record must be concrete");
    };
    let maps = topology
        .atoms
        .iter()
        .map(|atom| atom.atom_map())
        .collect::<Vec<_>>();
    assert_eq!(maps, vec![Some(5), Some(5), None, None, Some(7)]);
}

#[test]
fn v3k_atom_numbers_atom_map_c_locale_whitespace_and_no_conversion() {
    // Quoting keeps whitespace inside the V3000 token, after which RDKit's
    // C-locale atoi accepts SP/HT/VT/FF/CR. LF is covered at the private
    // helper boundary because it terminates a physical molfile line.
    let block = v3000_with_outer_and_inner(
        ZERO_OUTER,
        "9 0 0 0 0",
        &[
            "M  V30 1 C 0 0 0 \" 7\"",
            "M  V30 2 C 0 0 0 \"\t7\"",
            "M  V30 3 C 0 0 0 \"\u{000b}7\"",
            "M  V30 4 C 0 0 0 \"\u{000c}7\"",
            "M  V30 5 C 0 0 0 \"\r7\"",
            "M  V30 6 C 0 0 0 \"\"",
            "M  V30 7 C 0 0 0 \"+\"",
            "M  V30 8 C 0 0 0 \"-\"",
            "M  V30 9 C 0 0 0 abc",
        ],
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("C-locale map conversion")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(
        topology
            .atoms
            .iter()
            .map(|atom| atom.atom_map())
            .collect::<Vec<_>>(),
        vec![
            Some(7),
            Some(7),
            Some(7),
            Some(7),
            Some(7),
            None,
            None,
            None,
            None,
        ]
    );
}

#[test]
fn v3k_atom_numbers_atom_map_fixed_glibc_int_width_boundaries() {
    // Pinned RDKit 2026.03.1 on x86_64 glibc 2.43 was observed on every row.
    // Values outside C `int` range are fixed-environment regressions, not a
    // claim that C defines portable atoi behavior for those inputs.
    let block = v3000_with_outer_and_inner(
        ZERO_OUTER,
        "9 0 0 0 0",
        &[
            "M  V30 1 C 0 0 0 2147483647",
            "M  V30 2 C 0 0 0 2147483648",
            "M  V30 3 C 0 0 0 4294967296",
            "M  V30 4 C 0 0 0 4294967297",
            "M  V30 5 C 0 0 0 9223372036854775807",
            "M  V30 6 C 0 0 0 99999999999999999999",
            "M  V30 7 C 0 0 0 -2147483649",
            "M  V30 8 C 0 0 0 -4294967295",
            "M  V30 9 C 0 0 0 -99999999999999999999",
        ],
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("fixed glibc map boundaries")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(
        topology
            .atoms
            .iter()
            .map(|atom| atom.atom_map())
            .collect::<Vec<_>>(),
        vec![
            Some(i32::MAX as u32),
            None,
            None,
            Some(1),
            None,
            None,
            Some(i32::MAX as u32),
            Some(1),
            None,
        ]
    );
}

#[test]
fn v3k_atom_numbers_coordinate_row_alignment() {
    let block = v3000_with_outer_and_inner(
        ZERO_OUTER,
        "3 0 0 0 0",
        &[
            "M  V30 1 C 0.0 1.0 0 0",
            "M  V30 2 O 2.0 3.0 0 0",
            "M  V30 3 N 4.0 5.0 0 0",
        ],
    );
    let MolBlockRecord::Concrete {
        topology,
        coordinates,
        ..
    } = read_mol_block_detached(&block).expect("coordinate row alignment")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(topology.atoms.len(), 3);
    assert_eq!(coordinates.conformers_2d[0].coordinates().len(), 3);
    assert_eq!(coordinates.conformers_2d[0].coordinates()[2], [4.0, 5.0]);
}

#[test]
fn v3k_symbols_elements_mixed_case_second_char_is_normalized() {
    let block = v3000_with_outer_and_inner(
        ZERO_OUTER,
        "3 0 0 0 0",
        &[
            "M  V30 1 CL 0 0 0 0",
            "M  V30 2 BR 0 0 0 0",
            "M  V30 3 NA 0 0 0 0",
        ],
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("mixed-case elements")
    else {
        panic!("atom record must be concrete");
    };
    let numbers = topology
        .atoms
        .iter()
        .map(|atom| atom.element().atomic_number())
        .collect::<Vec<_>>();
    assert_eq!(numbers, vec![17, 35, 11]);
}

#[test]
fn v3k_symbols_elements_deuterium_and_tritium_isotopes() {
    let block = v3000_with_outer_and_inner(
        ZERO_OUTER,
        "2 0 0 0 0",
        &["M  V30 1 D 0 0 0 0", "M  V30 2 T 0 0 0 0"],
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("D/T shorthand")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(topology.atoms[0].element().atomic_number(), 1);
    assert_eq!(topology.atoms[0].isotope(), Some(2));
    assert_eq!(topology.atoms[1].element().atomic_number(), 1);
    assert_eq!(topology.atoms[1].isotope(), Some(3));
}

#[test]
fn v3k_symbols_elements_dummy_labels_for_pol_and_mod() {
    let block = v3000_with_outer_and_inner(
        ZERO_OUTER,
        "2 0 0 0 0",
        &["M  V30 1 Pol 0 0 0 0", "M  V30 2 Mod 0 0 0 0"],
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("Pol/Mod dummy labels")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(topology.atoms[0].element().atomic_number(), 0);
    assert_eq!(topology.atoms[0].prop("dummyLabel"), Some("Pol"));
    assert_eq!(topology.atoms[1].element().atomic_number(), 0);
    assert_eq!(topology.atoms[1].prop("dummyLabel"), Some("Mod"));
}

#[test]
fn v3k_symbols_elements_unknown_symbol_strictness() {
    let unknown = v3000_with_outer_and_inner(
        ZERO_OUTER,
        "2 0 0 0 0",
        &["M  V30 1 Zz 0 0 0 0", "M  V30 2 cl 0 0 0 0"],
    );
    match read_mol_block_detached(&unknown) {
        Err(SdfReadError::Parse(message)) => {
            assert!(
                message.contains("Element 'Zz' not found"),
                "unexpected diagnostic: {message}"
            );
        }
        other => panic!("expected a strict unknown-symbol failure, got {other:?}"),
    }

    let params = cosmolkit_io::MolBlockReadParams {
        strict_parsing: false,
        ..cosmolkit_io::MolBlockReadParams::default()
    };
    let MolBlockRecord::Concrete { topology, .. } =
        cosmolkit_io::read_mol_block_detached_with_params(&unknown, params)
            .expect("non-strict unknown symbols")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(topology.atoms[0].element().atomic_number(), 0);
    assert_eq!(topology.atoms[0].prop("dummyLabel"), Some("Zz"));
    assert_eq!(topology.atoms[1].element().atomic_number(), 0);
    assert_eq!(topology.atoms[1].prop("dummyLabel"), Some("cl"));
}

#[test]
fn v3k_symbols_elements_empty_symbol_row_is_malformed() {
    // The doubled-quote tokenizer rule preserves `""` as the symbol token; it is
    // not a valid element, so strict parsing reports a structured parse error
    // rather than fabricating an atom.
    let block = v3000_with_outer_and_inner(ZERO_OUTER, "1 0 0 0 0", &["M  V30 1 \"\" 0 0 0 0"]);
    assert!(
        matches!(read_mol_block_detached(&block), Err(SdfReadError::Parse(_))),
        "an empty quoted symbol must fail as a structured parse error"
    );
}

fn v3000_single_atom(token: &str) -> String {
    let atom = format!("M  V30 1 {token} 0 0 0 0");
    v3000_with_outer_and_inner(ZERO_OUTER, "1 0 0 0 0", &[&atom])
}

#[test]
fn v3k_symbols_lists_positive_with_spaces_and_empty_members() {
    // Quoting keeps `[C,  N,,]` as a single token; the inner spaces and empty
    // members are stripped/skipped.
    let block = v3000_single_atom("\"[C,  N,,]\"");
    let MolBlockRecord::Query(record) = read_mol_block_detached(&block).expect("atom list query")
    else {
        panic!("atom list must produce a query record");
    };
    let query_atom = record.query.atom(0).expect("query atom");
    let debug = format!("{:?}", query_atom.predicate());
    assert!(debug.contains("AtomicNumberIn([6, 7])"), "{debug}");
    assert_eq!(query_atom.atomic_number(), 0);
    assert!(!query_atom.no_implicit());
}

#[test]
fn v3k_symbols_lists_not_case_variants() {
    for token in ["NOT[C,N]", "not[C,N]"] {
        let block = v3000_single_atom(token);
        let MolBlockRecord::Query(record) = read_mol_block_detached(&block).expect("NOT atom list")
        else {
            panic!("NOT atom list must produce a query record");
        };
        let debug = format!(
            "{:?}",
            record.query.atom(0).expect("query atom").predicate()
        );
        assert!(
            debug.contains("AtomicNumberNotIn([6, 7])"),
            "{token}: {debug}"
        );
    }
}

#[test]
fn v3k_symbols_lists_not_on_non_list_is_rejected() {
    let block = v3000_single_atom("NOTC");
    match read_mol_block_detached(&block) {
        Err(SdfReadError::Parse(message)) => assert!(
            message.contains("NOT tokens only supported for atom lists"),
            "unexpected diagnostic: {message}"
        ),
        other => panic!("expected a NOT-on-non-list failure, got {other:?}"),
    }
}

#[test]
fn v3k_symbols_lists_all_empty_is_structurally_rejected() {
    for token in ["[]", "[,,]"] {
        let block = v3000_single_atom(token);
        assert!(
            matches!(read_mol_block_detached(&block), Err(SdfReadError::Parse(_))),
            "an all-empty atom list {token:?} must be a structured error"
        );
    }
}

#[test]
fn v3k_symbols_lists_missing_bracket_is_rejected() {
    let block = v3000_single_atom("[C,N");
    match read_mol_block_detached(&block) {
        Err(SdfReadError::Parse(message)) => {
            assert!(
                message.contains("Bad atom token"),
                "unexpected diagnostic: {message}"
            );
        }
        other => panic!("expected a malformed atom-list failure, got {other:?}"),
    }
}

#[test]
fn v3k_symbols_queries_complex_names_and_wildcard_are_queries_without_implicit_hydrogens() {
    // RDKit 2026.03.1 ParseV3000AtomSymbol routes all eight complex query
    // names through convertComplexNameToQuery and `*` through AtomNull; only
    // these symbol-query branches explicitly set noImplicit.
    let symbols = ["A", "AH", "Q", "QH", "X", "XH", "M", "MH", "*"];
    let atom_lines = symbols
        .iter()
        .enumerate()
        .map(|(index, symbol)| format!("M  V30 {} {symbol} 0 0 0 0", index + 1))
        .collect::<Vec<_>>();
    let atom_refs = atom_lines.iter().map(String::as_str).collect::<Vec<_>>();
    let block = v3000_with_outer_and_inner(
        ZERO_OUTER,
        &format!("{} 0 0 0 0", symbols.len()),
        &atom_refs,
    );
    let MolBlockRecord::Query(record) =
        read_mol_block_detached(&block).expect("complex V3000 query symbols")
    else {
        panic!("complex symbols must produce a query record");
    };
    for (index, symbol) in symbols.iter().enumerate() {
        let atom = record.query.atom(index).expect("query atom");
        assert!(atom.no_implicit(), "{symbol}");
        assert_eq!(atom.atomic_number(), 0, "{symbol}");
    }
    let wildcard = format!("{:?}", record.query.atom(8).unwrap().predicate());
    assert!(wildcard.contains("Any"), "{wildcard}");
}

#[test]
fn v3k_symbols_queries_r_groups_remain_concrete_with_source_metadata() {
    // ParseV3000AtomSymbol constructs ordinary atomic-number-zero atoms for
    // R/R#/R-labels, stores the full dummy label, and stores positive parsed
    // suffixes as isotopes. The source's string comparison also admits R100.
    let symbols = ["R", "R#", "R0", "R7", "R99", "R100"];
    let atom_lines = symbols
        .iter()
        .enumerate()
        .map(|(index, symbol)| format!("M  V30 {} {symbol} 0 0 0 0", index + 1))
        .collect::<Vec<_>>();
    let atom_refs = atom_lines.iter().map(String::as_str).collect::<Vec<_>>();
    let block = v3000_with_outer_and_inner(
        ZERO_OUTER,
        &format!("{} 0 0 0 0", symbols.len()),
        &atom_refs,
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("R-group symbols")
    else {
        panic!("R-group symbols must remain concrete atoms");
    };
    let expected_isotopes = [None, None, None, Some(7), Some(99), Some(100)];
    for ((atom, symbol), expected_isotope) in
        topology.atoms.iter().zip(symbols).zip(expected_isotopes)
    {
        assert_eq!(atom.atomic_number(), 0, "{symbol}");
        assert_eq!(atom.prop("dummyLabel"), Some(symbol), "{symbol}");
        assert_eq!(atom.isotope(), expected_isotope, "{symbol}");
        assert!(!atom.no_implicit(), "{symbol}");
    }
}

#[test]
fn v3k_symbols_queries_generic_labels_preserve_the_complete_source_key_set() {
    // This is the complete key set of RDKit 2026.03.1
    // GenericGroups::genericMatchers. Parsing only preserves the query atom
    // and atomLabel; it does not run or duplicate generic matching.
    let labels = [
        "Group",
        "G",
        "GroupH",
        "GH",
        "Group*",
        "G*",
        "GroupH*",
        "GH*",
        "Alkyl",
        "ALK",
        "AlkylH",
        "ALH",
        "Alkenyl",
        "AEL",
        "AlkenylH",
        "AEH",
        "Alkynyl",
        "AYL",
        "AlkynylH",
        "AYH",
        "Carbocyclic",
        "CBC",
        "CarbocyclicH",
        "CBH",
        "Carbocycloalkyl",
        "CAL",
        "CarbocycloalkylH",
        "CAH",
        "Carbocycloalkenyl",
        "CEL",
        "CarbocycloalkenylH",
        "CEH",
        "Carboaryl",
        "ARY",
        "CarboarylH",
        "ARH",
        "Cyclic",
        "CYC",
        "CyclicH",
        "CYH",
        "Acyclic",
        "ACY",
        "AcyclicH",
        "ACH",
        "Carboacyclic",
        "ABC",
        "CarboacyclicH",
        "ABH",
        "Heteroacyclic",
        "AHC",
        "HeteroacyclicH",
        "AHH",
        "Alkoxy",
        "AOX",
        "AlkoxyH",
        "AOH",
        "Heterocyclic",
        "CHC",
        "HeterocyclicH",
        "CHH",
        "Heteroaryl",
        "HAR",
        "HeteroarylH",
        "HAH",
        "NoCarbonRing",
        "CXX",
        "NoCarbonRingH",
        "CXH",
    ];
    let atom_lines = labels
        .iter()
        .enumerate()
        .map(|(index, label)| format!("M  V30 {} {label} 0 0 0 0", index + 1))
        .collect::<Vec<_>>();
    let atom_refs = atom_lines.iter().map(String::as_str).collect::<Vec<_>>();
    let block =
        v3000_with_outer_and_inner(ZERO_OUTER, &format!("{} 0 0 0 0", labels.len()), &atom_refs);
    let MolBlockRecord::Query(record) =
        read_mol_block_detached(&block).expect("generic group symbol table")
    else {
        panic!("generic labels must produce a query record");
    };
    for (index, label) in labels.iter().enumerate() {
        let atom = record.query.atom(index).expect("generic query atom");
        assert_eq!(atom.prop("atomLabel"), Some(*label), "{label}");
        assert_eq!(atom.atomic_number(), 0, "{label}");
        assert!(!atom.no_implicit(), "{label}");
        let debug = format!("{:?}", atom.predicate());
        assert!(debug.contains("AtomicNumber(0)"), "{label}: {debug}");
    }
}

#[test]
fn v3k_symbols_queries_negative_lookalikes_are_not_special_symbols() {
    // Source matching is case-sensitive and exact. R-label recognition uses
    // the pinned lexicographic condition; these tokens fall outside it.
    for symbol in ["a", "qh", "r7", "R9x", "R999", "group", "Alkylx", "ALKY"] {
        let block = v3000_single_atom(symbol);
        assert!(
            matches!(read_mol_block_detached(&block), Err(SdfReadError::Parse(_))),
            "lookalike {symbol:?} must use the normal strict element path"
        );
    }
}

#[test]
fn v3k_charge_concrete_atoms_preserve_signed_model_boundaries() {
    // RDKit 2026.03.1 ParseV3000AtomProps parses CHG with toInt and calls
    // setFormalCharge for atoms without a query. These values also exercise
    // both inclusive boundaries of the detached i8 charge representation.
    let block = v3000_block(
        &[
            "M  V30 1 C 0 0 0 0 CHG=7",
            "M  V30 2 N 0 0 0 0 CHG=-2",
            "M  V30 3 O 0 0 0 0 CHG=127",
            "M  V30 4 F 0 0 0 0 CHG=-128",
        ],
        4,
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("signed concrete charges")
    else {
        panic!("ordinary element rows must remain concrete");
    };
    let charges = topology
        .atoms
        .iter()
        .map(|atom| atom.formal_charge())
        .collect::<Vec<_>>();
    assert_eq!(charges, [7, -2, 127, -128]);
}

#[test]
fn v3k_charge_leading_plus_is_screened_but_not_converted_for_concrete_and_query_atoms() {
    // RDKit 2026.03.1 FileParserUtils::toInt permits '+' during character
    // screening, but signed std::from_chars does not accept a leading plus.
    // Its initially-zero destination is retained because the error is ignored.
    let concrete = v3000_block(&["M  V30 1 C 0 0 0 0 CHG=+7"], 1);
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&concrete).expect("screened concrete leading plus")
    else {
        panic!("ordinary element row must remain concrete");
    };
    assert_eq!(topology.atoms[0].formal_charge(), 0);

    let query = v3000_block(&["M  V30 1 * 0 0 0 0 CHG=+7"], 1);
    let MolBlockRecord::Query(record) =
        read_mol_block_detached(&query).expect("screened query leading plus")
    else {
        panic!("wildcard row must remain a query");
    };
    let atom = record.query.atom(0).expect("query atom");
    assert_eq!(atom.formal_charge(), 0);
    let QueryNode::And(children) = atom.predicate() else {
        panic!("expandQuery must conjoin the formal-charge predicate");
    };
    assert_eq!(children.len(), 2);
    assert!(matches!(
        &children[0],
        QueryNode::Predicate(AtomQueryPredicate::Any)
    ));
    assert!(matches!(
        &children[1],
        QueryNode::Predicate(AtomQueryPredicate::FormalCharge(0))
    ));
}

#[test]
fn v3k_charge_query_atoms_expand_the_query_without_setting_concrete_charge() {
    // ParseV3000AtomProps uses expandQuery(makeAtomFormalChargeQuery(charge))
    // when the atom already has a query, leaving its concrete charge field at
    // the source default.
    let block = v3000_block(&["M  V30 1 * 0 0 0 0 CHG=-3"], 1);
    let MolBlockRecord::Query(record) = read_mol_block_detached(&block).expect("query charge")
    else {
        panic!("wildcard with CHG must remain a query record");
    };
    let atom = record.query.atom(0).expect("query atom");
    assert_eq!(atom.formal_charge(), 0);
    let QueryNode::And(children) = atom.predicate() else {
        panic!("expandQuery must conjoin the formal-charge predicate");
    };
    assert_eq!(children.len(), 2);
    assert!(matches!(
        &children[0],
        QueryNode::Predicate(AtomQueryPredicate::Any)
    ));
    assert!(matches!(
        &children[1],
        QueryNode::Predicate(AtomQueryPredicate::FormalCharge(-3))
    ));
}

#[test]
fn v3k_charge_to_int_empty_invalid_and_source_overflow_behavior() {
    // FileParserUtils::toInt initializes its destination to zero and ignores
    // from_chars' result. Empty text and an out-of-int-range decimal therefore
    // produce zero, while the preceding character screen rejects other bytes.
    for (property, expected) in [
        ("CHG=", 0),
        ("CHG=2147483648", 0),
        ("CHG=--1", 0),
        ("\"CHG=   -7\"", -7),
        ("\"CHG=   \"", 0),
        ("\"CHG=7-\"", 7),
        ("\"CHG=7 8\"", 7),
        ("\"CHG=-2+9\"", -2),
    ] {
        let atom = format!("M  V30 1 C 0 0 0 0 {property}");
        let block = v3000_block(&[&atom], 1);
        let MolBlockRecord::Concrete { topology, .. } =
            read_mol_block_detached(&block).unwrap_or_else(|error| panic!("{property}: {error:?}"))
        else {
            panic!("{property}: concrete record expected");
        };
        assert_eq!(topology.atoms[0].formal_charge(), expected, "{property}");
    }

    for symbol in ["C", "*"] {
        for property in ["CHG=1x", "\"CHG=\t7\""] {
            let atom = format!("M  V30 1 {symbol} 0 0 0 0 {property}");
            let block = v3000_block(&[&atom], 1);
            assert!(matches!(
                read_mol_block_detached(&block),
                Err(SdfReadError::Field {
                    kind: "V3000 atom charge",
                    ..
                })
            ));
        }
    }
}

#[test]
fn v3k_charge_rejects_unrepresentable_concrete_and_query_values_without_narrowing() {
    // Pinned MolFileParser.cpp::ParseV3000AtomProps parses CHG into int,
    // narrows only through concrete Atom::setFormalCharge, and passes the
    // original int to makeAtomFormalChargeQuery for existing query atoms.
    for (charge, expected) in [("128", 128), ("-129", -129)] {
        let concrete_atom = format!("M  V30 1 C 0 0 0 0 CHG={charge}");
        let concrete_block = v3000_block(&[&concrete_atom], 1);
        assert!(
            matches!(
                read_mol_block_detached(&concrete_block),
                Err(SdfReadError::Unsupported(_))
            ),
            "concrete CHG={charge} keeps the existing i8 model boundary"
        );

        let query_atom = format!("M  V30 1 * 0 0 0 0 CHG={charge}");
        let query_block = v3000_block(&[&query_atom], 1);
        let MolBlockRecord::Query(record) =
            read_mol_block_detached(&query_block).expect("wide query charge remains representable")
        else {
            panic!("wildcard CHG={charge} must remain a query record");
        };
        let atom = record.query.atom(0).expect("query atom");
        assert_eq!(atom.formal_charge(), 0, "query carrier for CHG={charge}");
        assert_eq!(
            atom.predicate(),
            &QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::Any),
                QueryNode::predicate(AtomQueryPredicate::FormalCharge(expected)),
            ]),
            "full-width query target for CHG={charge}"
        );
    }
}

#[test]
fn v3k_radical_all_four_source_values_preserve_typed_atom_state() {
    // RDKit 2026.03.1 ParseV3000AtomProps treats RAD=0 as a no-op and maps
    // CTAB radical codes 1/2/3 to 2/1/2 radical electrons respectively.
    let block = v3000_block(
        &[
            "M  V30 1 C 0 0 0 0 RAD=0",
            "M  V30 2 N 0 0 0 0 RAD=1",
            "M  V30 3 O 0 0 0 0 RAD=2",
            "M  V30 4 F 0 0 0 0 RAD=3",
        ],
        4,
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("all V3000 RAD codes")
    else {
        panic!("ordinary element rows must remain concrete");
    };
    let radicals = topology
        .atoms
        .iter()
        .map(|atom| atom.radical_electrons())
        .collect::<Vec<_>>();
    assert_eq!(radicals, [0, 2, 1, 2]);
}

#[test]
fn v3k_radical_repeats_zero_noop_and_query_atoms_do_not_gain_query_chemistry() {
    // The source processes properties in order. RAD=0 executes `break`
    // without calling the setter, while later nonzero values overwrite the
    // atom field. Its explicit query FIXME does not add a query predicate.
    let concrete = v3000_block(
        &[
            "M  V30 1 C 0 0 0 0 RAD=1 RAD=0",
            "M  V30 2 N 0 0 0 0 RAD=2 RAD=3",
            "M  V30 3 O 0 0 0 0 RAD=+7",
            "M  V30 4 F 0 0 0 0 RAD=",
        ],
        4,
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&concrete).expect("repeated and no-conversion RAD values")
    else {
        panic!("ordinary element rows must remain concrete");
    };
    let radicals = topology
        .atoms
        .iter()
        .map(|atom| atom.radical_electrons())
        .collect::<Vec<_>>();
    assert_eq!(radicals, [2, 2, 0, 0]);

    let query = v3000_block(&["M  V30 1 * 0 0 0 0 RAD=2"], 1);
    let MolBlockRecord::Query(record) =
        read_mol_block_detached(&query).expect("query atom RAD state")
    else {
        panic!("wildcard row must remain a query");
    };
    let atom = record.query.atom(0).expect("query atom");
    assert_eq!(atom.radical_electrons(), 1);
    assert!(matches!(
        atom.predicate(),
        QueryNode::Predicate(AtomQueryPredicate::Any)
    ));
}

#[test]
fn v3k_radical_invalid_positive_negative_and_screened_text_are_structured_errors() {
    // ParseV3000AtomProps rejects parsed enum values outside 0..=3, while
    // FileParserUtils::toInt rejects bytes outside its screened character set.
    for value in ["-1", "4"] {
        let atom = format!("M  V30 1 C 0 0 0 0 RAD={value}");
        let block = v3000_block(&[&atom], 1);
        assert!(
            matches!(read_mol_block_detached(&block), Err(SdfReadError::Parse(_))),
            "RAD={value} must be an invalid-enum parse error"
        );
    }

    for symbol in ["C", "*"] {
        let atom = format!("M  V30 1 {symbol} 0 0 0 0 RAD=1x");
        let block = v3000_block(&[&atom], 1);
        assert!(matches!(
            read_mol_block_detached(&block),
            Err(SdfReadError::Field {
                kind: "V3000 radical",
                ..
            })
        ));
    }
}

#[test]
fn v3k_mass_integer_fractional_and_integer_conversion_edges_are_source_shaped() {
    // RDKit 2026.03.1 ParseV3000AtomProps tries FileParserUtils::toInt first
    // and only falls back to toDouble/floor after the integer character screen
    // throws. Thus leading plus and signed-int overflow retain toInt's zero,
    // while a decimal point selects the floating fallback.
    let block = v3000_block(
        &[
            "M  V30 1 C 0 0 0 0 MASS=13",
            "M  V30 2 N 0 0 0 0 MASS=13.9",
            "M  V30 3 O 0 0 0 0 MASS=+13",
            "M  V30 4 F 0 0 0 0 MASS=2147483648",
            "M  V30 5 Cl 0 0 0 0 \"MASS=   14.9\"",
            "M  V30 6 Br 0 0 0 0 MASS=1-2",
        ],
        6,
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("source-shaped V3000 MASS conversion")
    else {
        panic!("ordinary element rows must remain concrete");
    };
    let isotopes = topology
        .atoms
        .iter()
        .map(|atom| atom.isotope())
        .collect::<Vec<_>>();
    assert_eq!(
        isotopes,
        [Some(13), Some(13), Some(0), Some(0), Some(14), Some(1)]
    );
}

#[test]
fn v3k_mass_query_atoms_expand_typed_isotope_predicates() {
    // The source calls expandQuery(makeAtomIsotopeQuery(v)) for an existing
    // query atom. Assert the typed conjunction and leave concrete isotope
    // state unchanged instead of relying on a Debug rendering.
    for (property, expected) in [
        ("MASS=13", 13),
        ("MASS=13.9", 13),
        ("MASS=+13", 0),
        ("MASS=2147483648", 0),
        ("MASS=1-2", 1),
    ] {
        let atom_line = format!("M  V30 1 * 0 0 0 0 {property}");
        let block = v3000_block(&[&atom_line], 1);
        let MolBlockRecord::Query(record) =
            read_mol_block_detached(&block).unwrap_or_else(|error| panic!("{property}: {error:?}"))
        else {
            panic!("{property}: wildcard row must remain a query");
        };
        let atom = record.query.atom(0).expect("query atom");
        assert_eq!(atom.isotope(), None, "{property}");
        let QueryNode::And(children) = atom.predicate() else {
            panic!("{property}: isotope expansion must conjoin the predicate");
        };
        assert_eq!(children.len(), 2, "{property}");
        assert!(matches!(
            &children[0],
            QueryNode::Predicate(AtomQueryPredicate::Any)
        ));
        assert!(matches!(
            &children[1],
            QueryNode::Predicate(AtomQueryPredicate::Isotope(value)) if *value == expected
        ));
    }
}

#[test]
fn v3k_mass_negative_and_malformed_values_are_parse_errors() {
    // Negative integers and negative floored fractions reach the source's
    // explicit v < 0 failure. Text rejected by both toInt and toDouble is
    // converted to the same source MASS parse failure.
    for symbol in ["C", "*"] {
        for value in ["-1", "-0.1", "nope", "1x"] {
            let atom_line = format!("M  V30 1 {symbol} 0 0 0 0 MASS={value}");
            let block = v3000_block(&[&atom_line], 1);
            assert!(
                matches!(read_mol_block_detached(&block), Err(SdfReadError::Parse(_))),
                "{symbol} MASS={value} must be a source parse failure"
            );
        }
    }
}

#[test]
fn v3k_mass_unrepresentable_model_values_are_checked_without_narrowing() {
    // Pinned MolFileParser.cpp::ParseV3000AtomProps parses MASS to int,
    // narrows only through concrete Atom::setIsotope, and passes the original
    // int to makeAtomIsotopeQuery for existing query atoms.
    for (value, expected) in [("65536", 65_536), ("2147483647", i32::MAX)] {
        let concrete_atom = format!("M  V30 1 C 0 0 0 0 MASS={value}");
        let concrete_block = v3000_block(&[&concrete_atom], 1);
        assert!(
            matches!(
                read_mol_block_detached(&concrete_block),
                Err(SdfReadError::Unsupported(_))
            ),
            "concrete MASS={value} keeps the existing u16 model boundary"
        );

        let query_atom = format!("M  V30 1 * 0 0 0 0 MASS={value}");
        let query_block = v3000_block(&[&query_atom], 1);
        let MolBlockRecord::Query(record) = read_mol_block_detached(&query_block)
            .expect("wide query isotope remains representable")
        else {
            panic!("wildcard MASS={value} must remain a query record");
        };
        let atom = record.query.atom(0).expect("query atom");
        assert_eq!(atom.isotope(), None, "query carrier for MASS={value}");
        assert_eq!(
            atom.predicate(),
            &QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::Any),
                QueryNode::predicate(AtomQueryPredicate::Isotope(expected)),
            ]),
            "full-width query target for MASS={value}"
        );
    }
}

#[test]
fn v3k_atom_cfg_all_source_values_store_only_typed_parity_metadata() {
    // RDKit 2026.03.1 ParseV3000AtomProps treats CFG=0 as a no-op and stores
    // molParity for CFG=1/2/3; stereo perception belongs to later processing.
    let block = v3000_block(
        &[
            "M  V30 1 C 0 0 0 0 CFG=0",
            "M  V30 2 N 0 0 0 0 CFG=1",
            "M  V30 3 O 0 0 0 0 CFG=2",
            "M  V30 4 F 0 0 0 0 CFG=3",
        ],
        4,
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("all V3000 atom CFG values")
    else {
        panic!("ordinary element rows must remain concrete");
    };
    let parities = topology
        .atoms
        .iter()
        .map(|atom| atom.mol_parity())
        .collect::<Vec<_>>();
    assert_eq!(parities, [None, Some(1), Some(2), Some(3)]);

    let query = v3000_block(&["M  V30 1 * 0 0 0 0 CFG=2"], 1);
    let MolBlockRecord::Query(record) =
        read_mol_block_detached(&query).expect("query atom CFG metadata")
    else {
        panic!("wildcard row must remain a query");
    };
    let atom = record.query.atom(0).expect("query atom");
    assert_eq!(atom.mol_parity(), Some(2));
    assert_eq!(atom.prop("molParity"), Some("2"));
    assert!(matches!(
        atom.predicate(),
        QueryNode::Predicate(AtomQueryPredicate::Any)
    ));
}

#[test]
fn v3k_atom_cfg_repeats_and_source_integer_conversion_preserve_order() {
    let block = v3000_block(
        &[
            "M  V30 1 C 0 0 0 0 CFG=1 CFG=0",
            "M  V30 2 N 0 0 0 0 CFG=1 CFG=2",
            "M  V30 3 O 0 0 0 0 CFG=+2",
            "M  V30 4 F 0 0 0 0 CFG=",
            "M  V30 5 Cl 0 0 0 0 CFG=1-2",
        ],
        5,
    );
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("ordered CFG properties")
    else {
        panic!("ordinary element rows must remain concrete");
    };
    let parities = topology
        .atoms
        .iter()
        .map(|atom| atom.mol_parity())
        .collect::<Vec<_>>();
    assert_eq!(parities, [Some(1), Some(2), None, None, Some(1)]);
}

#[test]
fn v3k_atom_cfg_invalid_values_are_structured_errors() {
    for symbol in ["C", "*"] {
        for value in ["-1", "4"] {
            let atom_line = format!("M  V30 1 {symbol} 0 0 0 0 CFG={value}");
            let block = v3000_block(&[&atom_line], 1);
            assert!(matches!(
                read_mol_block_detached(&block),
                Err(SdfReadError::Parse(_))
            ));
        }
        let atom_line = format!("M  V30 1 {symbol} 0 0 0 0 CFG=1x");
        let block = v3000_block(&[&atom_line], 1);
        assert!(matches!(
            read_mol_block_detached(&block),
            Err(SdfReadError::Field {
                kind: "V3000 atom CFG",
                ..
            })
        ));
    }
}

#[test]
fn v3k_hcount_literal_zero_is_the_only_source_noop() {
    // RDKit 2026.03.1 ParseV3000AtomProps guards HCOUNT with the exact
    // textual comparison `val != "0"`; this spelling never converts a
    // concrete atom to a query and never expands an existing query atom.
    let concrete = v3000_block(&["M  V30 1 C 0 0 0 0 HCOUNT=0"], 1);
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&concrete).expect("literal-zero concrete HCOUNT")
    else {
        panic!("literal HCOUNT=0 must leave a concrete atom concrete");
    };
    assert_eq!(topology.atoms[0].atomic_number(), 6);

    let wildcard = v3000_block(&["M  V30 1 * 0 0 0 0 HCOUNT=0"], 1);
    let MolBlockRecord::Query(record) =
        read_mol_block_detached(&wildcard).expect("literal-zero wildcard HCOUNT")
    else {
        panic!("wildcard atom must remain a query");
    };
    assert_eq!(
        record.query.atoms()[0].predicate(),
        &QueryNode::predicate(AtomQueryPredicate::Any)
    );
}

#[test]
fn v3k_hcount_alternate_zero_and_nonpositive_values_expand_equality_zero() {
    // Every non-literal-zero spelling enters FileParserUtils::toInt. The
    // source normalizes -1 explicitly and sends every other nonpositive
    // result, including -2 and numeric zero text, to the same equality query.
    // QueryAtom::expandQuery replaces a non-negated AtomNull wildcard under
    // AND, so wildcard rows contain only the added equality predicate.
    for value in ["00", "-1", "-2"] {
        for symbol in ["C", "*"] {
            let atom_line = format!("M  V30 1 {symbol} 0 0 0 0 HCOUNT={value}");
            let block = v3000_block(&[&atom_line], 1);
            let MolBlockRecord::Query(record) = read_mol_block_detached(&block)
                .unwrap_or_else(|error| panic!("{symbol} HCOUNT={value}: {error:?}"))
            else {
                panic!("{symbol} HCOUNT={value} must produce query topology");
            };
            let hydrogen = QueryNode::predicate(AtomQueryPredicate::ImplicitHydrogenCount(0));
            let expected = if symbol == "C" {
                QueryNode::and(vec![
                    QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                    hydrogen,
                ])
            } else {
                hydrogen
            };
            assert_eq!(
                record.query.atoms()[0].predicate(),
                &expected,
                "{symbol} HCOUNT={value}"
            );
        }
    }
}

#[test]
fn v3k_hcount_positive_values_expand_less_equal_with_checked_model_width() {
    // The source constructs ATOM_LESSEQUAL_QUERY with the parsed C++ int.
    // COSMolKit preserves representable u8 targets and fails closed at 256
    // instead of silently narrowing the independent source query target. A
    // wildcard AtomNull is replaced by the added less-equal predicate.
    for value in [1_u8, 255] {
        for symbol in ["C", "*"] {
            let atom_line = format!("M  V30 1 {symbol} 0 0 0 0 HCOUNT={value}");
            let block = v3000_block(&[&atom_line], 1);
            let MolBlockRecord::Query(record) = read_mol_block_detached(&block)
                .unwrap_or_else(|error| panic!("{symbol} HCOUNT={value}: {error:?}"))
            else {
                panic!("{symbol} HCOUNT={value} must produce query topology");
            };
            let hydrogen =
                QueryNode::predicate(AtomQueryPredicate::ImplicitHydrogenCountLessEqual(value));
            let expected = if symbol == "C" {
                QueryNode::and(vec![
                    QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                    hydrogen,
                ])
            } else {
                hydrogen
            };
            assert_eq!(
                record.query.atoms()[0].predicate(),
                &expected,
                "{symbol} HCOUNT={value}"
            );
        }
    }

    for symbol in ["C", "*"] {
        let atom_line = format!("M  V30 1 {symbol} 0 0 0 0 HCOUNT=256");
        let block = v3000_block(&[&atom_line], 1);
        assert!(matches!(
            read_mol_block_detached(&block),
            Err(SdfReadError::Unsupported(message))
                if message.contains("outside the detached implicit-hydrogen-count query model")
        ));
    }
}

#[test]
fn v3k_hcount_concrete_conversion_skips_zero_isotope_state() {
    // RDKit 2026.03.1 QueryAtom(const Atom&) guards isotope expansion with
    // `if (other.getIsotope())`. MASS=0 stores zero on the concrete Atom, but
    // conversion for HCOUNT must therefore retain only atomic number before
    // appending the H-count predicate. MASS=+13 follows source toInt
    // no-conversion semantics and stores the same zero value.
    let expected = QueryNode::and(vec![
        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
        QueryNode::predicate(AtomQueryPredicate::ImplicitHydrogenCountLessEqual(1)),
    ]);
    for properties in ["MASS=0 HCOUNT=1", "MASS=+13 HCOUNT=1"] {
        let atom_line = format!("M  V30 1 C 0 0 0 0 {properties}");
        let block = v3000_block(&[&atom_line], 1);
        let MolBlockRecord::Query(record) = read_mol_block_detached(&block)
            .unwrap_or_else(|error| panic!("{properties}: {error:?}"))
        else {
            panic!("{properties} must produce query topology");
        };
        assert_eq!(
            record.query.atoms()[0].predicate(),
            &expected,
            "{properties}"
        );
    }
}

#[test]
fn v3k_hcount_concrete_conversion_retains_nonzero_constructor_predicates() {
    // QueryAtom(const Atom&) expands nonzero isotope, formal charge, and
    // radical predicates in that order before ParseV3000AtomProps appends the
    // HCOUNT predicate. Assert the complete typed tree so the zero-isotope
    // correction cannot discard meaningful constructor state.
    let cases = [
        (
            "MASS=13 HCOUNT=1",
            QueryNode::and(vec![
                QueryNode::and(vec![
                    QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                    QueryNode::predicate(AtomQueryPredicate::Isotope(13)),
                ]),
                QueryNode::predicate(AtomQueryPredicate::ImplicitHydrogenCountLessEqual(1)),
            ]),
        ),
        (
            "CHG=1 MASS=13 RAD=2 HCOUNT=1",
            QueryNode::and(vec![
                QueryNode::and(vec![
                    QueryNode::and(vec![
                        QueryNode::and(vec![
                            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                            QueryNode::predicate(AtomQueryPredicate::Isotope(13)),
                        ]),
                        QueryNode::predicate(AtomQueryPredicate::FormalCharge(1)),
                    ]),
                    QueryNode::predicate(AtomQueryPredicate::NumRadicalElectrons(1)),
                ]),
                QueryNode::predicate(AtomQueryPredicate::ImplicitHydrogenCountLessEqual(1)),
            ]),
        ),
    ];
    for (properties, expected) in cases {
        let atom_line = format!("M  V30 1 C 0 0 0 0 {properties}");
        let block = v3000_block(&[&atom_line], 1);
        let MolBlockRecord::Query(record) = read_mol_block_detached(&block)
            .unwrap_or_else(|error| panic!("{properties}: {error:?}"))
        else {
            panic!("{properties} must produce query topology");
        };
        assert_eq!(
            record.query.atoms()[0].predicate(),
            &expected,
            "{properties}"
        );
    }
}

#[test]
fn v3k_hcount_then_mass_zero_retains_source_query_expansion() {
    // Property order is observable: after HCOUNT has replaced the concrete
    // atom with a QueryAtom, the MASS branch directly expands Isotope(0).
    // This is distinct from constructor conversion of pre-existing zero state.
    let block = v3000_block(&["M  V30 1 C 0 0 0 0 HCOUNT=1 MASS=0"], 1);
    let MolBlockRecord::Query(record) =
        read_mol_block_detached(&block).expect("HCOUNT before MASS=0")
    else {
        panic!("HCOUNT=1 MASS=0 must produce query topology");
    };
    assert_eq!(
        record.query.atoms()[0].predicate(),
        &QueryNode::and(vec![
            QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                QueryNode::predicate(AtomQueryPredicate::ImplicitHydrogenCountLessEqual(1)),
            ]),
            QueryNode::predicate(AtomQueryPredicate::Isotope(0)),
        ])
    );
}

#[test]
fn v3k_unsat_activates_only_for_exact_raw_one() {
    // Pinned ParseV3000AtomProps compares the raw string_view to "1" and
    // never numerically normalizes UNSAT. QueryAtom's null algebra reduces
    // wildcard AND UNSAT to the UNSAT predicate itself.
    let concrete = v3000_block(&["M  V30 1 C 0 0 0 0 UNSAT=1"], 1);
    let MolBlockRecord::Query(record) =
        read_mol_block_detached(&concrete).expect("exact UNSAT=1 on carbon")
    else {
        panic!("exact UNSAT=1 must convert concrete carbon to a query");
    };
    assert_eq!(
        record.query.atoms()[0].predicate(),
        &QueryNode::and(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
            QueryNode::predicate(AtomQueryPredicate::IsUnsaturated),
        ])
    );

    let wildcard = v3000_block(&["M  V30 1 * 0 0 0 0 UNSAT=1"], 1);
    let MolBlockRecord::Query(record) =
        read_mol_block_detached(&wildcard).expect("exact UNSAT=1 on wildcard")
    else {
        panic!("wildcard must remain a query");
    };
    assert_eq!(
        record.query.atoms()[0].predicate(),
        &QueryNode::predicate(AtomQueryPredicate::IsUnsaturated)
    );

    let list = v3000_block(&["M  V30 1 [C,N] 0 0 0 0 UNSAT=1"], 1);
    let MolBlockRecord::Query(record) =
        read_mol_block_detached(&list).expect("exact UNSAT=1 on atom list")
    else {
        panic!("atom list must remain a query");
    };
    assert_eq!(
        record.query.atoms()[0].predicate(),
        &QueryNode::and(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomicNumberIn(vec![6, 7])),
            QueryNode::predicate(AtomQueryPredicate::IsUnsaturated),
        ])
    );
}

#[test]
fn v3k_unsat_nonactivating_spellings_preserve_concrete_and_query_state() {
    // `01`, numeric-looking alternatives, empty text, and arbitrary text are
    // all source no-ops because the branch performs no conversion.
    for value in ["01", "0", "+1", "-1", "", "other"] {
        let concrete_line = format!("M  V30 1 C 0 0 0 0 UNSAT={value}");
        let concrete = v3000_block(&[&concrete_line], 1);
        let MolBlockRecord::Concrete { topology, .. } = read_mol_block_detached(&concrete)
            .unwrap_or_else(|error| panic!("concrete UNSAT={value}: {error:?}"))
        else {
            panic!("concrete UNSAT={value} must remain concrete");
        };
        assert_eq!(topology.atoms[0].atomic_number(), 6, "UNSAT={value}");

        let wildcard_line = format!("M  V30 1 * 0 0 0 0 UNSAT={value}");
        let wildcard = v3000_block(&[&wildcard_line], 1);
        let MolBlockRecord::Query(record) = read_mol_block_detached(&wildcard)
            .unwrap_or_else(|error| panic!("wildcard UNSAT={value}: {error:?}"))
        else {
            panic!("wildcard UNSAT={value} must remain a query");
        };
        assert_eq!(
            record.query.atoms()[0].predicate(),
            &QueryNode::predicate(AtomQueryPredicate::Any),
            "UNSAT={value}"
        );
    }
}

#[test]
fn v3k_unsat_concrete_conversion_preserves_source_constructor_state() {
    // QueryAtom(const Atom&) skips numeric isotope zero but retains nonzero
    // isotope, formal charge, and radical predicates in constructor order.
    let without_isotope = QueryNode::and(vec![
        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
        QueryNode::predicate(AtomQueryPredicate::IsUnsaturated),
    ]);
    for properties in ["MASS=0 UNSAT=1", "MASS=+13 UNSAT=1"] {
        let atom_line = format!("M  V30 1 C 0 0 0 0 {properties}");
        let block = v3000_block(&[&atom_line], 1);
        let MolBlockRecord::Query(record) = read_mol_block_detached(&block)
            .unwrap_or_else(|error| panic!("{properties}: {error:?}"))
        else {
            panic!("{properties} must produce query topology");
        };
        assert_eq!(
            record.query.atoms()[0].predicate(),
            &without_isotope,
            "{properties}"
        );
    }

    let block = v3000_block(&["M  V30 1 C 0 0 0 0 MASS=13 CHG=1 RAD=2 UNSAT=1"], 1);
    let MolBlockRecord::Query(record) =
        read_mol_block_detached(&block).expect("nonzero constructor state before UNSAT")
    else {
        panic!("nonzero constructor state with UNSAT must produce query topology");
    };
    assert_eq!(
        record.query.atoms()[0].predicate(),
        &QueryNode::and(vec![
            QueryNode::and(vec![
                QueryNode::and(vec![
                    QueryNode::and(vec![
                        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                        QueryNode::predicate(AtomQueryPredicate::Isotope(13)),
                    ]),
                    QueryNode::predicate(AtomQueryPredicate::FormalCharge(1)),
                ]),
                QueryNode::predicate(AtomQueryPredicate::NumRadicalElectrons(1)),
            ]),
            QueryNode::predicate(AtomQueryPredicate::IsUnsaturated),
        ])
    );
}

#[test]
fn v3k_unsat_then_mass_zero_retains_explicit_query_expansion() {
    // Once UNSAT has converted the atom, the later MASS branch directly
    // expands Isotope(0); constructor zero suppression must not erase it.
    let block = v3000_block(&["M  V30 1 C 0 0 0 0 UNSAT=1 MASS=0"], 1);
    let MolBlockRecord::Query(record) =
        read_mol_block_detached(&block).expect("UNSAT before MASS=0")
    else {
        panic!("UNSAT=1 MASS=0 must produce query topology");
    };
    assert_eq!(
        record.query.atoms()[0].predicate(),
        &QueryNode::and(vec![
            QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                QueryNode::predicate(AtomQueryPredicate::IsUnsaturated),
            ]),
            QueryNode::predicate(AtomQueryPredicate::Isotope(0)),
        ])
    );
}

#[test]
fn v3k_rbcnt_literal_zero_is_the_only_source_noop() {
    // Pinned RDKit 2026.03.1 ParseV3000AtomProps guards RBCNT with the exact
    // textual comparison `val != "0"`. The literal spelling neither stores
    // molRingBondCount nor converts a concrete Atom to QueryAtom.
    let block = v3000_block(&["M  V30 1 C 0 0 0 0 RBCNT=0"], 1);
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("literal-zero RBCNT")
    else {
        panic!("literal RBCNT=0 must leave a concrete atom concrete");
    };
    assert_eq!(topology.atoms[0].atomic_number(), 6);
    assert_eq!(topology.atoms[0].prop("molRingBondCount"), None);
}

#[test]
fn v3k_rbcnt_source_values_build_equality_queries_and_retain_metadata() {
    // ParseV3000AtomProps maps -1 to equality zero, -2 to the unsigned
    // 0xDEADBEEF equality sentinel, and values above four to equality four.
    // Unlike the V2000 M  RBC path, V3000 never builds LESS-EQUAL here.
    for (value, target) in [("00", 0), ("-1", 0), ("1", 1), ("4", 4), ("5", 4)] {
        let atom_line = format!("M  V30 1 C 0 0 0 0 RBCNT={value}");
        let block = v3000_block(&[&atom_line], 1);
        let MolBlockRecord::Query(record) = read_mol_block_detached(&block)
            .unwrap_or_else(|error| panic!("RBCNT={value}: {error:?}"))
        else {
            panic!("RBCNT={value} must produce query topology");
        };
        assert_eq!(
            record.query.atoms()[0].predicate(),
            &QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                QueryNode::predicate(AtomQueryPredicate::RingBondCount(target)),
            ]),
            "RBCNT={value}"
        );
        assert_eq!(
            record.query.atoms()[0].prop("molRingBondCount"),
            Some(match value {
                "00" => "0",
                "5" => "5",
                other => other,
            }),
            "RBCNT={value}"
        );
        assert_eq!(record.properties.prop("_NeedsQueryScan"), None);
        assert_eq!(record.query.prop("_NeedsQueryScan"), None);
    }

    // FileParserUtils::toInt uses from_chars: leading '+' performs no
    // conversion and leaves the initially-zero destination unchanged.
    let block = v3000_block(&["M  V30 1 C 0 0 0 0 RBCNT=+1"], 1);
    let MolBlockRecord::Query(record) =
        read_mol_block_detached(&block).expect("leading-plus RBCNT")
    else {
        panic!("non-literal RBCNT=+1 must produce query topology");
    };
    assert_eq!(
        record.query.atoms()[0].predicate(),
        &QueryNode::and(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
            QueryNode::predicate(AtomQueryPredicate::RingBondCount(0)),
        ])
    );
    assert_eq!(record.query.atoms()[0].prop("molRingBondCount"), Some("0"));
}

#[test]
fn v3k_rbcnt_minus_two_sets_deferred_scan_state_and_unsigned_sentinel() {
    // The source stores the original signed property, assigns the sentinel in
    // an unsigned destination, and marks both molecule/query record state for
    // the later ring-info scan.
    let block = v3000_block(&["M  V30 1 C 0 0 0 0 RBCNT=-2"], 1);
    let MolBlockRecord::Query(record) = read_mol_block_detached(&block).expect("RBCNT=-2") else {
        panic!("RBCNT=-2 must produce query topology");
    };
    assert_eq!(
        record.query.atoms()[0].predicate(),
        &QueryNode::and(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
            QueryNode::predicate(AtomQueryPredicate::RingBondCount(0xDEAD_BEEF_u32 as i32)),
        ])
    );
    assert_eq!(record.query.atoms()[0].prop("molRingBondCount"), Some("-2"));
    assert_eq!(record.properties.prop("_NeedsQueryScan"), Some("1"));
    assert_eq!(record.query.prop("_NeedsQueryScan"), Some("1"));
}

#[test]
fn v3k_rbcnt_preserves_constructor_and_property_order_semantics() {
    // Concrete-to-query conversion uses QueryAtom(const Atom&): numeric
    // isotope zero is omitted, while nonzero isotope/charge/radical state is
    // retained in constructor order before the RBCNT equality predicate.
    let block = v3000_block(&["M  V30 1 C 0 0 0 0 MASS=13 CHG=1 RAD=2 RBCNT=1"], 1);
    let MolBlockRecord::Query(record) =
        read_mol_block_detached(&block).expect("nonzero state before RBCNT")
    else {
        panic!("nonzero state with RBCNT must produce query topology");
    };
    assert_eq!(
        record.query.atoms()[0].predicate(),
        &QueryNode::and(vec![
            QueryNode::and(vec![
                QueryNode::and(vec![
                    QueryNode::and(vec![
                        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                        QueryNode::predicate(AtomQueryPredicate::Isotope(13)),
                    ]),
                    QueryNode::predicate(AtomQueryPredicate::FormalCharge(1)),
                ]),
                QueryNode::predicate(AtomQueryPredicate::NumRadicalElectrons(1)),
            ]),
            QueryNode::predicate(AtomQueryPredicate::RingBondCount(1)),
        ])
    );

    let before = v3000_block(&["M  V30 1 C 0 0 0 0 MASS=0 RBCNT=1"], 1);
    let MolBlockRecord::Query(record) =
        read_mol_block_detached(&before).expect("MASS=0 before RBCNT")
    else {
        panic!("MASS=0 RBCNT=1 must produce query topology");
    };
    assert_eq!(
        record.query.atoms()[0].predicate(),
        &QueryNode::and(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
            QueryNode::predicate(AtomQueryPredicate::RingBondCount(1)),
        ])
    );

    // Once RBCNT has converted the atom, the later MASS branch explicitly
    // expands Isotope(0); constructor zero suppression must not erase it.
    let after = v3000_block(&["M  V30 1 C 0 0 0 0 RBCNT=1 MASS=0"], 1);
    let MolBlockRecord::Query(record) =
        read_mol_block_detached(&after).expect("RBCNT before MASS=0")
    else {
        panic!("RBCNT=1 MASS=0 must produce query topology");
    };
    assert_eq!(
        record.query.atoms()[0].predicate(),
        &QueryNode::and(vec![
            QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                QueryNode::predicate(AtomQueryPredicate::RingBondCount(1)),
            ]),
            QueryNode::predicate(AtomQueryPredicate::Isotope(0)),
        ])
    );
}

#[test]
fn v3k_rbcnt_wildcard_null_query_and_signed_negative_targets_match_source() {
    // QueryAtom null-query algebra replaces wildcard AtomNull under AND.
    let block = v3000_block(&["M  V30 1 * 0 0 0 0 RBCNT=1"], 1);
    let MolBlockRecord::Query(record) = read_mol_block_detached(&block).expect("wildcard RBCNT")
    else {
        panic!("wildcard RBCNT must remain query topology");
    };
    assert_eq!(
        record.query.atoms()[0].predicate(),
        &QueryNode::predicate(AtomQueryPredicate::RingBondCount(1))
    );

    // Negative targets use AtomRingQuery's source-defined nonzero-count
    // branch and remain signed all the way through the detached query target.
    for (value, target) in [("-3", -3), ("-2147483648", i32::MIN)] {
        let atom_line = format!("M  V30 1 C 0 0 0 0 RBCNT={value}");
        let block = v3000_block(&[&atom_line], 1);
        let MolBlockRecord::Query(record) =
            read_mol_block_detached(&block).expect("signed RBCNT query")
        else {
            panic!("RBCNT={value} must produce query topology");
        };
        assert_eq!(
            record.query.atoms()[0].predicate(),
            &QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                QueryNode::predicate(AtomQueryPredicate::RingBondCount(target)),
            ])
        );
    }
}

#[test]
fn v3k_atom_int_props_literal_zero_skips_each_distinct_property() {
    // Pinned RDKit 2026.03.1 ParseV3000AtomProps compares the raw value to
    // exactly "0" before calling toInt for every scalar branch in this table.
    let properties = [
        ("VAL", Some("molTotValence")),
        ("STBOX", Some("molStereoCare")),
        ("SUBST", Some("molSubstCount")),
        ("EXACHG", Some("molRxnExactChange")),
        ("INVRET", None),
        ("SEQID", Some("molAtomSeqId")),
    ];
    for (label, key) in properties {
        let atom_line = format!("M  V30 1 C 0 0 0 0 {label}=0");
        let block = v3000_block(&[&atom_line], 1);
        let MolBlockRecord::Concrete { topology, .. } =
            read_mol_block_detached(&block).unwrap_or_else(|error| panic!("{label}=0: {error:?}"))
        else {
            panic!("{label}=0 must remain concrete");
        };
        let atom = &topology.atoms[0];
        if let Some(key) = key {
            assert_eq!(atom.prop(key), None, "{label}=0");
        } else {
            assert_eq!(atom.mol_inversion_flag(), None, "{label}=0");
        }
    }
}

#[test]
fn v3k_atom_int_props_share_source_conversion_without_aliasing_keys() {
    // The non-literal spellings all enter FileParserUtils::toInt. These rows
    // cover numeric zero, ordinary signed values, leading-plus no-conversion,
    // and out-of-range from_chars no-conversion while asserting each branch's
    // independent detached storage.
    let properties = [
        ("VAL", Some("molTotValence")),
        ("STBOX", Some("molStereoCare")),
        ("SUBST", Some("molSubstCount")),
        ("EXACHG", Some("molRxnExactChange")),
        ("INVRET", None),
        ("SEQID", Some("molAtomSeqId")),
    ];
    for (text, expected) in [
        ("00", 0),
        ("7", 7),
        ("-2", -2),
        ("+7", 0),
        ("2147483648", 0),
    ] {
        for (label, key) in properties {
            let atom_line = format!("M  V30 1 C 0 0 0 0 {label}={text}");
            let block = v3000_block(&[&atom_line], 1);
            let MolBlockRecord::Concrete { topology, .. } = read_mol_block_detached(&block)
                .unwrap_or_else(|error| panic!("{label}={text}: {error:?}"))
            else {
                panic!("{label}={text} must remain concrete");
            };
            let atom = &topology.atoms[0];
            if let Some(key) = key {
                let expected_text = expected.to_string();
                assert_eq!(
                    atom.prop(key),
                    Some(expected_text.as_str()),
                    "{label}={text}"
                );
                for other_key in [
                    "molTotValence",
                    "molStereoCare",
                    "molSubstCount",
                    "molRxnExactChange",
                    "molAtomSeqId",
                ] {
                    if other_key != key {
                        assert_eq!(
                            atom.prop(other_key),
                            None,
                            "{label} must not set {other_key}"
                        );
                    }
                }
                assert_eq!(atom.mol_inversion_flag(), None, "{label}={text}");
            } else {
                assert_eq!(atom.mol_inversion_flag(), Some(expected), "{label}={text}");
                for other_key in [
                    "molTotValence",
                    "molStereoCare",
                    "molSubstCount",
                    "molRxnExactChange",
                    "molAtomSeqId",
                ] {
                    assert_eq!(
                        atom.prop(other_key),
                        None,
                        "INVRET must not set {other_key}"
                    );
                }
            }
        }
    }
}

#[test]
fn v3k_atom_int_props_repeated_values_replace_and_invalid_text_is_structured() {
    // RDProps::setProp replaces an existing value with the same key. The
    // typed inversion field has the same last-write-wins behavior.
    let properties = [
        ("VAL", Some("molTotValence"), "V3000 total valence"),
        ("STBOX", Some("molStereoCare"), "V3000 stereo care"),
        ("SUBST", Some("molSubstCount"), "V3000 substitution count"),
        (
            "EXACHG",
            Some("molRxnExactChange"),
            "V3000 exact-change flag",
        ),
        ("INVRET", None, "V3000 inversion flag"),
        ("SEQID", Some("molAtomSeqId"), "V3000 sequence id"),
    ];
    for (label, key, kind) in properties {
        let atom_line = format!("M  V30 1 C 0 0 0 0 {label}=1 {label}=-2");
        let block = v3000_block(&[&atom_line], 1);
        let MolBlockRecord::Concrete { topology, .. } = read_mol_block_detached(&block)
            .unwrap_or_else(|error| panic!("repeated {label}: {error:?}"))
        else {
            panic!("repeated {label} must remain concrete");
        };
        if let Some(key) = key {
            assert_eq!(topology.atoms[0].prop(key), Some("-2"), "{label}");
        } else {
            assert_eq!(topology.atoms[0].mol_inversion_flag(), Some(-2));
        }

        let invalid_line = format!("M  V30 1 C 0 0 0 0 {label}=12x");
        let invalid = v3000_block(&[&invalid_line], 1);
        assert!(matches!(
            read_mol_block_detached(&invalid),
            Err(SdfReadError::Field {
                kind: actual_kind,
                value,
                ..
            }) if actual_kind == kind && value == "12x"
        ));
    }
}

#[test]
fn v3k_rgroups_zero_count_is_an_exact_noop_for_concrete_and_query_atoms() {
    // Pinned ParseV3000RGroups loops exactly nRs times. A zero count ignores
    // following tokens and performs no concrete-to-query replacement.
    let concrete = v3000_block(&["M  V30 1 C 0 0 0 0 RGROUPS=(0 7)"], 1);
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&concrete).expect("zero concrete RGROUPS")
    else {
        panic!("zero RGROUPS must leave concrete carbon concrete");
    };
    let atom = &topology.atoms[0];
    assert_eq!(atom.atomic_number(), 6);
    assert_eq!(atom.isotope(), None);
    assert_eq!(atom.prop("_MolFileRLabel"), None);
    assert_eq!(atom.prop("dummyLabel"), None);

    let wildcard = v3000_block(&["M  V30 1 * 0 0 0 0 RGROUPS=(0)"], 1);
    let MolBlockRecord::Query(record) =
        read_mol_block_detached(&wildcard).expect("zero query RGROUPS")
    else {
        panic!("wildcard atom remains a query");
    };
    assert_eq!(
        record.query.atoms()[0].predicate(),
        &QueryNode::predicate(AtomQueryPredicate::Any)
    );
    assert_eq!(record.query.atoms()[0].isotope(), None);
}

#[test]
fn v3k_rgroups_one_and_multiple_labels_apply_source_order_with_last_label_winning() {
    // Every source-loop iteration overwrites both label properties, isotope,
    // and the query with AtomNull. The selected labels remain ordered even
    // though the last iteration supplies the observable scalar state.
    for (value, isotope, label) in [("(1 7)", 7, "7"), ("(2 7 12)", 12, "12")] {
        let atom_line = format!("M  V30 1 C 0 0 0 0 RGROUPS={value}");
        let block = v3000_block(&[&atom_line], 1);
        let MolBlockRecord::Query(record) =
            read_mol_block_detached(&block).expect("counted R-group list")
        else {
            panic!("nonempty RGROUPS must construct query state");
        };
        let atom = &record.query.atoms()[0];
        assert_eq!(atom.atomic_number(), 6, "{value}");
        assert_eq!(atom.isotope(), Some(isotope), "{value}");
        assert_eq!(atom.prop("_MolFileRLabel"), Some(label), "{value}");
        let dummy = format!("R{label}");
        assert_eq!(atom.prop("dummyLabel"), Some(dummy.as_str()), "{value}");
        assert_eq!(
            atom.predicate(),
            &QueryNode::predicate(AtomQueryPredicate::Any),
            "{value}"
        );
    }
}

#[test]
fn v3k_rgroups_uses_unsigned_lexical_cast_sign_semantics_and_checked_isotope_width() {
    // boost::lexical_cast<unsigned int> accepts signs and converts a negative
    // magnitude modulo 2^32. -4294967295 therefore yields label/isotope 1.
    let block = v3000_block(&["M  V30 1 C 0 0 0 0 RGROUPS=(1 -4294967295)"], 1);
    let MolBlockRecord::Query(record) =
        read_mol_block_detached(&block).expect("source-shaped signed label")
    else {
        panic!("nonempty RGROUPS must construct query state");
    };
    let atom = &record.query.atoms()[0];
    assert_eq!(atom.isotope(), Some(1));
    assert_eq!(atom.prop("_MolFileRLabel"), Some("1"));
    assert_eq!(atom.prop("dummyLabel"), Some("R1"));

    // RDKit narrows larger unsigned labels into Atom::d_isotope (uint16_t).
    // The detached model intentionally keeps a checked boundary instead of
    // silently reproducing that narrowing.
    for value in ["65536", "-1"] {
        let atom_line = format!("M  V30 1 C 0 0 0 0 RGROUPS=(1 {value})");
        let block = v3000_block(&[&atom_line], 1);
        assert!(matches!(
            read_mol_block_detached(&block),
            Err(SdfReadError::Unsupported(message))
                if message.contains("R-group labels outside the detached u16 isotope model")
        ));
    }
}

#[test]
fn v3k_rgroups_count_syntax_and_integer_failures_are_structured() {
    // ParseV3000RGroups requires parens, preserves empty tokens from repeated
    // literal spaces, consumes selected integers in full, and checks the
    // declared count before entering its update loop.
    for value in ["1", "(2 7)", "(one)", "(1 nope)", "(1  7)", "(1 7x)"] {
        let atom_line = format!("M  V30 1 C 0 0 0 0 RGROUPS={value}");
        let block = v3000_block(&[&atom_line], 1);
        assert!(
            matches!(read_mol_block_detached(&block), Err(SdfReadError::Parse(_))),
            "{value}"
        );
    }
}

#[test]
fn v3k_rgroups_replaces_existing_query_state_and_later_properties_expand_from_null() {
    // replaceAtomWithQueryAtom is a no-op for an existing QueryAtom, but each
    // RGROUPS iteration then setQuery(AtomNull), discarding the old tree.
    let before = v3000_block(&["M  V30 1 C 0 0 0 0 HCOUNT=1 RGROUPS=(1 +7)"], 1);
    let MolBlockRecord::Query(record) = read_mol_block_detached(&before).unwrap() else {
        panic!("HCOUNT/RGROUPS record must be a query");
    };
    let atom = &record.query.atoms()[0];
    assert_eq!(atom.isotope(), Some(7));
    assert_eq!(atom.prop("dummyLabel"), Some("R7"));
    assert_eq!(
        atom.predicate(),
        &QueryNode::predicate(AtomQueryPredicate::Any)
    );

    // A later HCOUNT expands against AtomNull. Pinned QueryAtom null algebra
    // replaces the null tree instead of retaining Any as an AND child.
    let after = v3000_block(&["M  V30 1 C 0 0 0 0 RGROUPS=(1 7) HCOUNT=1"], 1);
    let MolBlockRecord::Query(record) = read_mol_block_detached(&after).unwrap() else {
        panic!("RGROUPS/HCOUNT record must be a query");
    };
    let atom = &record.query.atoms()[0];
    assert_eq!(atom.isotope(), Some(7));
    assert_eq!(
        atom.predicate(),
        &QueryNode::predicate(AtomQueryPredicate::ImplicitHydrogenCountLessEqual(1))
    );
}

fn v3000_single_atom_x(x: &str) -> String {
    let atom = format!("M  V30 1 C {x} 0 0 0");
    v3000_with_outer_and_inner(ZERO_OUTER, "1 0 0 0 0", &[&atom])
}

#[test]
fn v3k_atom_numbers_hexadecimal_forms_exact() {
    // A non-zero z selects 3D storage so all three axes are observable.
    let block = v3000_with_outer_and_inner(
        ZERO_OUTER,
        "2 0 0 0 0",
        &[
            "M  V30 1 C 0x1p+1 0x1.8p+1 0x.8p0 0",
            "M  V30 2 O 0X1P+1 0x10 -0x1p+1 0",
        ],
    );
    let MolBlockRecord::Concrete { coordinates, .. } =
        read_mol_block_detached(&block).expect("hexadecimal coordinates")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(coordinates.conformers_3d.len(), 1);
    let rows = coordinates.conformers_3d[0].coordinates();
    assert_eq!(rows[0], [2.0, 3.0, 0.5]);
    assert_eq!(rows[1], [2.0, 16.0, -2.0]);
}

#[test]
fn v3k_atom_numbers_hexadecimal_partial_and_prefix() {
    let cases = [
        ("0x10", 16.0),
        ("0x1", 1.0),
        ("0x1p", 1.0),
        ("0x1p-1", 0.5),
        ("0x.8p0", 0.5),
        ("0x", 0.0),
        ("0xg", 0.0),
    ];
    for (token, expected) in cases {
        let block = v3000_single_atom_x(token);
        let MolBlockRecord::Concrete { coordinates, .. } =
            read_mol_block_detached(&block).unwrap_or_else(|error| panic!("{token}: {error:?}"))
        else {
            panic!("{token}: atom record must be concrete");
        };
        assert_eq!(
            coordinates.conformers_2d[0].coordinates()[0][0],
            expected,
            "{token}"
        );
    }
}

#[test]
fn v3k_atom_numbers_signed_zero_and_subnormal_boundaries() {
    let block = v3000_with_outer_and_inner(
        ZERO_OUTER,
        "4 0 0 0 0",
        &[
            "M  V30 1 C -0 0 0 0",
            "M  V30 2 O 0x1p-1074 0 0 0",
            "M  V30 3 N 0x1.8p-1075 0 0 0",
            "M  V30 4 F 0x1p-1075 0 0 0",
        ],
    );
    let MolBlockRecord::Concrete { coordinates, .. } =
        read_mol_block_detached(&block).expect("subnormal boundaries")
    else {
        panic!("atom record must be concrete");
    };
    let rows = coordinates.conformers_2d[0].coordinates();
    assert_eq!(rows[0][0], 0.0);
    assert!(rows[0][0].is_sign_negative());
    // 2^-1074 is the smallest positive subnormal.
    assert_eq!(rows[1][0], f64::from_bits(1));
    // 1.5 * 2^-1075 rounds up to the smallest subnormal.
    assert_eq!(rows[2][0], f64::from_bits(1));
    // Exactly 2^-1075 is a tie and rounds to even (zero).
    assert_eq!(rows[3][0], 0.0);
}

#[test]
fn v3k_atom_numbers_overflow_and_nonfinite_rejected() {
    for token in ["0x1p+2000", "0x1p+1024", "1e400"] {
        let block = v3000_single_atom_x(token);
        assert!(
            matches!(
                read_mol_block_detached(&block),
                Err(SdfReadError::Coordinates(_))
            ),
            "overflowing coordinate {token:?} must fail the finite-coordinate boundary"
        );
    }
    // Gradual underflow to zero is a finite coordinate and is accepted.
    let block = v3000_single_atom_x("0x1p-2000");
    let MolBlockRecord::Concrete { coordinates, .. } =
        read_mol_block_detached(&block).expect("underflow to zero")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(coordinates.conformers_2d[0].coordinates()[0][0], 0.0);
}

#[test]
fn v3k_symbols_elements_multibyte_symbol_does_not_panic() {
    // U+1F600 is four bytes long, so byte index 3 is not a character boundary.
    let block = v3000_single_atom("\u{1F600}");
    assert!(
        matches!(read_mol_block_detached(&block), Err(SdfReadError::Parse(_))),
        "strict unknown multibyte symbol must be a structured parse error"
    );
    let params = cosmolkit_io::MolBlockReadParams {
        strict_parsing: false,
        ..cosmolkit_io::MolBlockReadParams::default()
    };
    let MolBlockRecord::Concrete { topology, .. } =
        cosmolkit_io::read_mol_block_detached_with_params(&block, params)
            .expect("non-strict multibyte dummy")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(topology.atoms[0].element().atomic_number(), 0);
    assert_eq!(topology.atoms[0].prop("dummyLabel"), Some("\u{1F600}"));
}

#[test]
fn v3k_symbols_elements_multibyte_prefixes_do_not_panic() {
    // Each token has text before the four-byte code point, so the old
    // `symbol[..3]` prefix check split a code point; all must be structured.
    for token in ["N\u{1F600}", "NO\u{1F600}", "NOT\u{1F600}"] {
        let block = v3000_single_atom(token);
        assert!(
            read_mol_block_detached(&block).is_err(),
            "{token:?} must produce a structured error, not a panic"
        );
    }
}

#[test]
fn v3k_symbols_lists_multibyte_entry_does_not_panic() {
    let block = v3000_single_atom("\"[\u{1F600},C]\"");
    assert!(
        matches!(read_mol_block_detached(&block), Err(SdfReadError::Parse(_))),
        "an unknown multibyte list entry must be a structured parse error"
    );
}

#[test]
fn v3k_atom_numbers_atof_port_counterexamples() {
    // 0x10p followed by fifty 9s overflows to +inf and must fail the existing
    // finite-coordinate boundary (previously a silent 0.0).
    let overflow = format!("0x10p{}", "9".repeat(50));
    let block = v3000_single_atom_x(&overflow);
    assert!(
        matches!(
            read_mol_block_detached(&block),
            Err(SdfReadError::Coordinates(_))
        ),
        "huge positive hex exponent must be rejected as non-finite"
    );

    // 0x0.1p- followed by fifty 9s underflows to +0.0 and is accepted as a
    // finite coordinate (previously a spurious non-finite error).
    let underflow = format!("0x0.1p-{}", "9".repeat(50));
    let block = v3000_single_atom_x(&underflow);
    let MolBlockRecord::Concrete { coordinates, .. } =
        read_mol_block_detached(&block).expect("huge negative hex exponent underflows to zero")
    else {
        panic!("atom record must be concrete");
    };
    let value = coordinates.conformers_2d[0].coordinates()[0][0];
    assert_eq!(value, 0.0);
    assert!(!value.is_sign_negative());
}

#[test]
fn v3k_atom_numbers_atof_port_finite_bits_and_signed_zero() {
    let cases = [
        ("0x1.0000000000001p0", f64::from_bits(0x3ff0_0000_0000_0001)),
        ("0x1p-1074", f64::from_bits(0x0000_0000_0000_0001)),
        ("0x1.8p+1", 3.0),
    ];
    for (token, expected) in cases {
        let block = v3000_single_atom_x(token);
        let MolBlockRecord::Concrete { coordinates, .. } =
            read_mol_block_detached(&block).unwrap_or_else(|error| panic!("{token}: {error:?}"))
        else {
            panic!("{token}: atom record must be concrete");
        };
        assert_eq!(
            coordinates.conformers_2d[0].coordinates()[0][0].to_bits(),
            expected.to_bits(),
            "{token}"
        );
    }

    // Negative zero stays negative through the detached reader.
    let block = v3000_single_atom_x("-0");
    let MolBlockRecord::Concrete { coordinates, .. } =
        read_mol_block_detached(&block).expect("negative zero coordinate")
    else {
        panic!("atom record must be concrete");
    };
    let value = coordinates.conformers_2d[0].coordinates()[0][0];
    assert_eq!(value.to_bits(), 0x8000_0000_0000_0000);
}

#[test]
fn v3k_atom_numbers_atof_port_suffix_and_malformed_prefixes() {
    let cases = [
        ("12xyz", 12.0),
        ("0x", 0.0),
        ("abc", 0.0),
        ("1e", 1.0),
        ("  3.5", 3.5),
        ("1,5", 1.0),
    ];
    for (token, expected) in cases {
        let block = v3000_single_atom_x(token);
        let MolBlockRecord::Concrete { coordinates, .. } =
            read_mol_block_detached(&block).unwrap_or_else(|error| panic!("{token}: {error:?}"))
        else {
            panic!("{token}: atom record must be concrete");
        };
        assert_eq!(
            coordinates.conformers_2d[0].coordinates()[0][0],
            expected,
            "{token}"
        );
    }
}

#[test]
fn v3k_atom_numbers_atof_port_nan_rejected_and_alignment_unchanged() {
    for token in ["nan", "-nan", "0x1p+1024", "1e400"] {
        let block = v3000_single_atom_x(token);
        assert!(
            matches!(
                read_mol_block_detached(&block),
                Err(SdfReadError::Coordinates(_))
            ),
            "non-finite coordinate {token:?} must fail the finite-coordinate boundary"
        );
    }

    // Row alignment and atom ordering are unchanged by the conversion swap.
    let block = v3000_with_outer_and_inner(
        ZERO_OUTER,
        "3 0 0 0 0",
        &[
            "M  V30 1 C 0x1p+1 1.0 0 0",
            "M  V30 2 O 2.0 0X1P+1 0 0",
            "M  V30 3 N 4.0 5.0 0 0",
        ],
    );
    let MolBlockRecord::Concrete {
        topology,
        coordinates,
        ..
    } = read_mol_block_detached(&block).expect("conversion row alignment")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(topology.atoms.len(), 3);
    let rows = coordinates.conformers_2d[0].coordinates();
    assert_eq!(rows.len(), 3);
    assert_eq!(rows[0], [2.0, 1.0]);
    assert_eq!(rows[1], [2.0, 2.0]);
    assert_eq!(rows[2], [4.0, 5.0]);
}

#[test]
fn v3k_atom_numbers_atof_port_musl_long_hex_source_value() {
    // Retained excluded-scope discrepancy, NOT RDKit acceptance evidence:
    // this musl-derived implementation rounds a long hexadecimal significand
    // one ULP below glibc/RDKit; see IO-atof-port.md. Narrowing the coordinate
    // contract must not erase this regression or claim all-input equivalence.
    let block = v3000_single_atom_x("0x1814D89994224121");
    let MolBlockRecord::Concrete { coordinates, .. } =
        read_mol_block_detached(&block).expect("long hex significand")
    else {
        panic!("atom record must be concrete");
    };
    assert_eq!(
        coordinates.conformers_2d[0].coordinates()[0][0].to_bits(),
        0x43b8_14d8_9994_2242
    );
}

#[test]
fn v3k_atom_numbers_coordinate_contract_decimal_and_missing_fields() {
    // Pinned ParseV3000AtomBlock uses raw atof (unlike V2000 toDouble), but
    // requires each coordinate/map token to exist before conversion.
    for (token, expected) in [
        ("+1.25e2", 125.0_f64),
        ("-1.25E-2", -0.0125),
        ("-0.000e+12", -0.0),
        ("1.25suffix", 1.25),
        ("abc", 0.0),
    ] {
        let MolBlockRecord::Concrete { coordinates, .. } =
            read_mol_block_detached(&v3000_single_atom_x(token)).expect("V3000 prefix")
        else {
            panic!("expected concrete molecule");
        };
        assert_eq!(
            coordinates.conformers_2d[0].coordinates()[0][0].to_bits(),
            expected.to_bits(),
            "{token}"
        );
    }
    for row in [
        "M  V30 1 C",
        "M  V30 1 C 1",
        "M  V30 1 C 1 2",
        "M  V30 1 C 1 2 3",
    ] {
        let block = v3000_with_outer_and_inner(ZERO_OUTER, "1 0 0 0 0", &[row]);
        let error = read_mol_block_detached(&block).expect_err("missing required token");
        assert!(matches!(error, SdfReadError::Parse(_)), "{row}: {error:?}");
        assert!(error.to_string().contains("Bad atom line"), "{error:?}");
    }
}

fn v3k_record_finish_block(query: bool) -> String {
    if query {
        return concat!(
            "query finish\n",
            "  COSMolKit 2D\n",
            "query comments\n",
            "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS 2 1 1 0 1\n",
            "M  V30 BEGIN ATOM\n",
            "M  V30 10 * 0 0 0 0 RBCNT=-2 ATTCHORD=(2 2 L)\n",
            "M  V30 20 N 1 0 0 0\n",
            "M  V30 END ATOM\n",
            "M  V30 BEGIN BOND\n",
            "M  V30 50 1 10 20\n",
            "M  V30 END BOND\n",
            "M  V30 BEGIN SGROUP\n",
            "M  V30 1 DAT 0 ATOMS=(1 10) FIELDNAME=note FIELDDATA=value\n",
            "M  V30 END SGROUP\n",
            "M  V30 BEGIN COLLECTION\n",
            "M  V30 MDLV30/STERAC1 ATOMS=(1 1)\n",
            "M  V30 END COLLECTION\n",
            "M  V30 END CTAB\n",
            "M  END\n",
        )
        .to_owned();
    }

    concat!(
        "concrete finish\n",
        "  COSMolKit 3D\n",
        "concrete comments\n",
        "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
        "M  V30 BEGIN CTAB\n",
        "M  V30 COUNTS 3 2 1 0 1\n",
        "M  V30 BEGIN ATOM\n",
        "M  V30 10 C 1 2 3 0 ATTCHORD=(2 2 A)\n",
        "M  V30 20 C 4 5 6 0 CFG=1\n",
        "M  V30 30 O 7 8 9 0\n",
        "M  V30 END ATOM\n",
        "M  V30 BEGIN BOND\n",
        "M  V30 50 1 10 20\n",
        "M  V30 60 1 20 30\n",
        "M  V30 END BOND\n",
        "M  V30 BEGIN SGROUP\n",
        "M  V30 1 SUP 0 ATOMS=(1 20) PATOMS=(1 20) XBONDS=(2 50 60) -\n",
        "M  V30 XBHEAD=(2 2 1) XBCORR=(2 1 1) -\n",
        "M  V30 BRKXYZ=(9 1 2 3 4 5 6 7 8 9) -\n",
        "M  V30 CSTATE=(4 50 1.25 -2.5 3.75) SAP=(4 20 10 AP) LABEL=unit\n",
        "M  V30 END SGROUP\n",
        "M  V30 BEGIN COLLECTION\n",
        "M  V30 MDLV30/STEREL1 ATOMS=(1 2)\n",
        "M  V30 END COLLECTION\n",
        "M  V30 END CTAB\n",
        "M  END\n",
    )
    .to_owned()
}

#[test]
fn v3k_record_finish_concrete_retains_complete_typed_state_and_properties() {
    // Pinned ParseV3000CTAB installs the parsed atom/bond/SGroup/collection
    // state, sets the molfile properties, classifies the conformer, and only
    // then returns a complete record. External bookmarks never become row IDs.
    let MolBlockRecord::Concrete {
        topology,
        coordinates,
        properties,
    } = read_mol_block_detached(&v3k_record_finish_block(false))
        .expect("cross-feature concrete V3000 record")
    else {
        panic!("ordinary atoms must finish as a concrete record");
    };

    assert_eq!(
        topology
            .atoms
            .iter()
            .map(|atom| atom.id())
            .collect::<Vec<_>>(),
        [AtomId::new(0), AtomId::new(1), AtomId::new(2),]
    );
    assert_eq!(
        topology
            .bonds
            .iter()
            .map(|bond| bond.id())
            .collect::<Vec<_>>(),
        [BondId::new(0), BondId::new(1),]
    );
    assert_eq!(topology.bonds[0].begin(), AtomId::new(0));
    assert_eq!(topology.bonds[0].end(), AtomId::new(1));
    let attachment = topology.atoms[0]
        .template_attachment_order()
        .expect("typed template attachment");
    assert_eq!(attachment.entries()[0].target(), AtomId::new(1));
    assert_eq!(attachment.entries()[0].label(), "A");

    let group = &topology.substance_groups[0];
    assert_eq!(group.atoms(), [AtomId::new(1)]);
    assert_eq!(group.parent_atoms(), [AtomId::new(1)]);
    assert_eq!(group.bonds(), [BondId::new(0), BondId::new(1)]);
    assert_eq!(
        group.head_crossing_bonds(),
        [BondId::new(1), BondId::new(0)]
    );
    assert_eq!(
        group.crossing_bond_correspondence(),
        [BondId::new(0), BondId::new(0)]
    );
    assert_eq!(
        group.display().unwrap().brackets()[0].points()[2],
        [7.0, 8.0, 9.0]
    );
    assert_eq!(group.cstates()[0].bond(), BondId::new(0));
    assert_eq!(group.cstates()[0].vector(), &[1.25, -2.5, 3.75]);
    assert_eq!(group.attach_points()[0].atom, AtomId::new(1));
    assert_eq!(group.attach_points()[0].leaving_atom, Some(AtomId::new(0)));
    assert_eq!(group.label(), Some("unit"));
    assert_eq!(topology.stereo_groups.len(), 1);
    assert_eq!(topology.stereo_groups[0].kind(), StereoGroupKind::Or);
    assert_eq!(topology.stereo_groups[0].id(), Some(1));
    assert_eq!(topology.stereo_groups[0].atoms(), [AtomId::new(1)]);

    assert!(coordinates.conformers_2d.is_empty());
    assert_eq!(coordinates.conformers_3d.len(), 1);
    assert!(coordinates.conformers_3d[0].is_3d());
    assert_eq!(
        coordinates.conformers_3d[0].coordinates(),
        &[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
    );
    assert_eq!(
        coordinates.source_coordinate_dim,
        Some(CoordinateDimension::ThreeD)
    );
    assert_eq!(properties.name(), Some("concrete finish"));
    assert_eq!(properties.prop("_MolFileInfo"), Some("  COSMolKit 3D"));
    assert_eq!(
        properties.prop("_MolFileComments"),
        Some("concrete comments")
    );
    assert_eq!(properties.prop("_MolFileChiralFlag"), Some("1"));
}

#[test]
fn v3k_record_finish_query_retains_carriers_metadata_and_deferred_scan_state() {
    let MolBlockRecord::Query(record) = read_mol_block_detached(&v3k_record_finish_block(true))
        .expect("cross-feature query V3000 record")
    else {
        panic!("wildcard/RBCNT state must finish as a query record");
    };

    assert_eq!(record.query.atoms().len(), 2);
    assert_eq!(record.query.atoms()[0].id(), AtomId::new(0));
    assert_eq!(record.query.atoms()[1].id(), AtomId::new(1));
    assert_eq!(record.query.bonds()[0].id(), BondId::new(0));
    assert_eq!(record.query.bonds()[0].begin(), AtomId::new(0));
    assert_eq!(record.query.bonds()[0].end(), AtomId::new(1));
    let attachment = record.query.atoms()[0]
        .template_attachment_order()
        .expect("query carrier keeps typed attachment");
    assert_eq!(attachment.entries()[0].target(), AtomId::new(1));
    assert_eq!(attachment.entries()[0].label(), "L");
    assert_eq!(
        record.query.coordinates_2d(),
        Some(&[[0.0, 0.0], [1.0, 0.0]][..])
    );
    assert!(record.query.conformers_3d().is_empty());
    assert_eq!(
        record.source_coordinate_dim,
        Some(CoordinateDimension::TwoD)
    );
    assert_eq!(record.query.stereo_groups()[0].kind(), StereoGroupKind::And);
    assert_eq!(record.query.stereo_groups()[0].atoms(), [AtomId::new(0)]);
    let groups = query_substance_groups(&record.query);
    assert_eq!(groups[0].atoms(), [AtomId::new(0)]);
    assert_eq!(
        groups[0].data().unwrap().field_name.as_deref(),
        Some("note")
    );
    assert_eq!(groups[0].data().unwrap().values, ["value"]);
    assert_eq!(record.properties.name(), Some("query finish"));
    assert_eq!(
        record.properties.prop("_MolFileInfo"),
        Some("  COSMolKit 2D")
    );
    assert_eq!(
        record.properties.prop("_MolFileComments"),
        Some("query comments")
    );
    assert_eq!(record.properties.prop("_MolFileChiralFlag"), Some("1"));
    assert_eq!(record.properties.prop("_NeedsQueryScan"), Some("1"));
    assert_eq!(record.query.name(), Some("query finish"));
    assert_eq!(record.query.prop("_NeedsQueryScan"), Some("1"));
}

#[test]
fn v3k_record_finish_validation_errors_return_no_partial_record() {
    // The canonical constructors validate attachment and stereo references
    // before any MolBlockRecord variant can be returned.
    let invalid_attachment =
        v3k_record_finish_block(true).replace("ATTCHORD=(2 2 L)", "ATTCHORD=(2 3 L)");
    assert!(matches!(
        read_mol_block_detached(&invalid_attachment),
        Err(SdfReadError::QueryGraph(_))
    ));

    let invalid_stereo = v3k_record_finish_block(false)
        .replace("MDLV30/STEREL1 ATOMS=(1 2)", "MDLV30/STEREL1 ATOMS=(1 99)");
    assert!(matches!(
        read_mol_block_detached(&invalid_stereo),
        Err(SdfReadError::Parse(_))
    ));
}
