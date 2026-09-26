use std::error::Error;

use cosmolkit_cx::{
    CxAtomConstraint, CxAtomProperty, CxBondReference, CxCoordinateBondKind, CxCoordinateBonds,
    CxCoordinates, CxCountConstraint, CxDataSGroup, CxDoubleBondStereo, CxDoubleBondStereoKind,
    CxEnhancedStereo, CxLinkNode, CxParseError, CxPolymerSGroup, CxRadical, CxRecord, CxRingBond,
    CxSGroupHierarchy, CxStereoGroupKind, CxVariableAttachment, CxWedgeBond, CxWedgeDirection,
    ParsedCxExtensions, parse_cx_extensions,
};

#[test]
fn public_record_vocabulary_is_complete_and_value_semantic() {
    let records = vec![
        CxRecord::Coordinates(CxCoordinates {
            conformer: 2,
            values: vec![Some([1.0, 2.0, 3.0]), None],
            is_3d: true,
        }),
        CxRecord::AtomLabels(vec![Some("label".into()), None]),
        CxRecord::AtomValues(vec![None, Some("value".into())]),
        CxRecord::AtomProperties(vec![CxAtomProperty {
            atom: 1,
            name: "key".into(),
            value: "value".into(),
        }]),
        CxRecord::CoordinateBonds(CxCoordinateBonds {
            kind: CxCoordinateBondKind::Dative,
            bonds: vec![CxBondReference { atom: 2, bond: 3 }],
        }),
        CxRecord::ZeroBonds(vec![4]),
        CxRecord::EnhancedStereo(CxEnhancedStereo {
            kind: CxStereoGroupKind::And,
            group_id: 5,
            atoms: vec![1, 3],
        }),
        CxRecord::Unsaturation(vec![6]),
        CxRecord::RingBonds(vec![CxRingBond {
            atom: 7,
            constraint: CxCountConstraint::LessEqual(4),
        }]),
        CxRecord::LinkNodes(vec![CxLinkNode {
            atom: 8,
            start_repetitions: 1,
            end_repetitions: 4,
            outer_atoms: Some([2, 9]),
        }]),
        CxRecord::DataSGroup(CxDataSGroup {
            atoms: vec![1, 0],
            field_name: "field".into(),
            data: "data".into(),
            query_op: "=".into(),
            field_info: "unit".into(),
            field_tag: "tag".into(),
            coordinates: Some("(1,2)".into()),
        }),
        CxRecord::SGroupHierarchy(vec![CxSGroupHierarchy {
            parent: 1,
            children: vec![2, 3],
        }]),
        CxRecord::PolymerSGroup(CxPolymerSGroup {
            type_code: "n".into(),
            atoms: vec![2, 1],
            label: "repeat".into(),
            connect: "ht".into(),
            head_crossings: vec![4],
            tail_crossings: vec![5],
        }),
        CxRecord::Substitution(vec![CxAtomConstraint {
            atom: 9,
            constraint: CxCountConstraint::QueryScan,
        }]),
        CxRecord::VariableAttachments(vec![CxVariableAttachment {
            atom: 2,
            endpoints: vec![3, 5, 4],
        }]),
        CxRecord::WedgedBonds(vec![CxWedgeBond {
            atom: 1,
            bond: 0,
            direction: CxWedgeDirection::BeginDash,
            configuration: 3,
        }]),
        CxRecord::DoubleBondStereo(CxDoubleBondStereo {
            stereo: CxDoubleBondStereoKind::Trans,
            bonds: vec![2],
        }),
        CxRecord::Radicals(vec![CxRadical {
            atom: 4,
            electrons: 1,
        }]),
        CxRecord::Unknown("future:data".into()),
    ];

    let parsed = ParsedCxExtensions::new(records.clone(), 37);
    assert_eq!(parsed.records(), records.as_slice());
    assert_eq!(parsed.consumed(), 37);
    assert_eq!(parsed.clone().into_records(), records);
}

#[test]
fn positioned_error_has_exact_public_fields_display_and_error_contract() {
    fn accepts_error(_: &dyn Error) {}

    let error = CxParseError::new(12, "bad CX token");
    accepts_error(&error);
    assert_eq!(error.offset, 12);
    assert_eq!(error.message, "bad CX token");
    assert_eq!(error.to_string(), "bad CX token at byte 12");
    assert_eq!(error.clone(), error);
}

#[test]
fn scanner_decodes_source_decimal_entities_in_every_text_family() {
    let labels = parse_cx_extensions("|$alpha&#59;beta;naïve&#38;x;a&#;b$|").unwrap();
    assert_eq!(
        labels.records(),
        &[CxRecord::AtomLabels(vec![
            Some("alpha;beta".into()),
            Some("naïve&x".into()),
            Some("ab".into()),
        ])]
    );

    let values = parse_cx_extensions("|$_AV:a&#44;b$|").unwrap();
    assert_eq!(
        values.records(),
        &[CxRecord::AtomValues(vec![Some("a,b".into())])]
    );

    let properties = parse_cx_extensions("|atomProp:0.na&#109;e.v&#38;x|").unwrap();
    assert_eq!(
        properties.records(),
        &[CxRecord::AtomProperties(vec![CxAtomProperty {
            atom: 0,
            name: "name".into(),
            value: "v&x".into(),
        }])]
    );

    let sgroup = parse_cx_extensions("|SgD:0:F&#73;ELD:d&#59;x::::|").unwrap();
    let CxRecord::DataSGroup(sgroup) = &sgroup.records()[0] else {
        panic!("expected data SGroup")
    };
    assert_eq!(sgroup.field_name, "FIELD");
    assert_eq!(sgroup.data, "d;x");
}

#[test]
fn scanner_reproduces_pinned_int_to_char_narrowing_and_string_lift() {
    let input = "|$plain;x&#0;y;x&#256;y;x&#321;y;&#128;;&#255;;&#2147483647;;&#;$| suffix";
    let parsed = parse_cx_extensions(input).unwrap();
    assert_eq!(
        parsed.records(),
        &[CxRecord::AtomLabels(vec![
            Some("plain".into()),
            Some("x\0y".into()),
            Some("x\0y".into()),
            Some("xAy".into()),
            Some("\u{80}".into()),
            Some("\u{ff}".into()),
            Some("\u{ff}".into()),
            None,
        ])]
    );
    assert_eq!(parsed.consumed(), input.find(" suffix").unwrap());
}

#[test]
fn scanner_applies_entity_conversion_in_property_and_sgroup_text_consumers() {
    let properties = parse_cx_extensions("|atomProp:0.n&#321;me.v&#256;x|").unwrap();
    assert_eq!(
        properties.records(),
        &[CxRecord::AtomProperties(vec![CxAtomProperty {
            atom: 0,
            name: "nAme".into(),
            value: "v\0x".into(),
        }])]
    );

    let sgroup = parse_cx_extensions("|SgD:0:F&#255;:d&#321;ta::::|").unwrap();
    let CxRecord::DataSGroup(sgroup) = &sgroup.records()[0] else {
        panic!("expected data SGroup")
    };
    assert_eq!(sgroup.field_name, "F\u{ff}");
    assert_eq!(sgroup.data, "dAta");
}

#[test]
fn scanner_reports_exact_entity_integer_and_pair_failure_offsets() {
    let unterminated = parse_cx_extensions("|$x&#59$|").unwrap_err();
    assert_eq!(unterminated.offset, 3);
    assert_eq!(
        unterminated.message,
        "failure parsing CXSMILES extensions: quoted block not terminated with ';'"
    );

    let character_overflow = parse_cx_extensions("|$x&#2147483648;$|").unwrap_err();
    assert_eq!(character_overflow.offset, 3);
    assert_eq!(character_overflow.message, "invalid CX character code");

    let integer_overflow = parse_cx_extensions("|C:4294967296.0|").unwrap_err();
    assert_eq!(integer_overflow.offset, 3);
    assert_eq!(integer_overflow.message, "invalid CX integer");

    let missing_pair_member = parse_cx_extensions("|C:1.|").unwrap_err();
    assert_eq!(missing_pair_member.offset, 5);
    assert_eq!(missing_pair_member.message, "invalid CX integer");
}

#[test]
fn list_scans_preserve_order_and_separate_adjacent_records() {
    let parsed = parse_cx_extensions("|u:3,1,,rb:2:0,s:4:*|").unwrap();
    assert_eq!(
        parsed.records(),
        &[
            CxRecord::Unsaturation(vec![3, 1]),
            CxRecord::Unknown(",".into()),
            CxRecord::RingBonds(vec![CxRingBond {
                atom: 2,
                constraint: CxCountConstraint::Exact(0),
            }]),
            CxRecord::Substitution(vec![CxAtomConstraint {
                atom: 4,
                constraint: CxCountConstraint::QueryScan,
            }]),
        ]
    );
}

#[test]
fn unknown_records_are_lossless_stop_before_known_records_and_make_progress() {
    let parsed = parse_cx_extensions("|###:x,y,rb:2:3|").unwrap();
    assert_eq!(
        parsed.records(),
        &[
            CxRecord::Unknown("###:x,y,".into()),
            CxRecord::RingBonds(vec![CxRingBond {
                atom: 2,
                constraint: CxCountConstraint::Exact(3),
            }]),
        ]
    );

    let consecutive = parse_cx_extensions("|,,,###,,,|").unwrap();
    assert_eq!(
        consecutive.records(),
        &[CxRecord::Unknown(",,,###,,,".into())]
    );

    let embedded_dispatch = parse_cx_extensions("|future:a,b,rb:2:3|").unwrap_err();
    assert_eq!(embedded_dispatch.offset, 3);
    assert_eq!(
        embedded_dispatch.message,
        "expected ':', found CX syntax mismatch"
    );
}

#[test]
fn byte_offsets_and_consumption_remain_valid_around_utf8_payloads() {
    let input = "|$α;β$| molecule-name";
    let parsed = parse_cx_extensions(input).unwrap();
    assert_eq!(
        parsed.records(),
        &[CxRecord::AtomLabels(vec![
            Some("α".into()),
            Some("β".into()),
        ])]
    );
    assert_eq!(parsed.consumed(), input.find(" molecule-name").unwrap());

    let malformed = "|$α&#59$|";
    let error = parse_cx_extensions(malformed).unwrap_err();
    assert_eq!(error.offset, malformed.find('&').unwrap());
}

#[test]
fn dependency_direction_keeps_cx_below_both_real_consumers() {
    let cx_manifest = include_str!("../Cargo.toml");
    let dependencies = cx_manifest.split_once("[dependencies]").unwrap().1.trim();
    assert!(
        dependencies.is_empty(),
        "unexpected CX dependency: {dependencies}"
    );

    let smiles_manifest = include_str!("../../cosmolkit-smiles/Cargo.toml");
    let search_manifest = include_str!("../../cosmolkit-search/Cargo.toml");
    assert!(smiles_manifest.contains("cosmolkit-cx ="));
    assert!(search_manifest.contains("cosmolkit-cx ="));
    for forbidden in [
        "cosmolkit-smiles",
        "cosmolkit-search",
        "cosmolkit-core",
        "cosmolkit-model",
    ] {
        assert!(!dependencies.contains(forbidden));
    }
}
