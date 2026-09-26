use cosmolkit_cx::{
    CxCoordinateBondKind, CxCountConstraint, CxDoubleBondStereoKind, CxParseError, CxRecord,
    CxStereoGroupKind, CxWedgeDirection, parse_cx_extensions,
};

fn assert_error(input: &str, offset: usize, message: &str) {
    assert_eq!(
        parse_cx_extensions(input).unwrap_err(),
        CxParseError {
            offset,
            message: message.to_owned(),
        }
    );
}

#[test]
fn empty_blocks_pipe_contract_and_suffix_consumption_are_exact() {
    let empty = parse_cx_extensions("").unwrap();
    assert!(empty.records().is_empty());
    assert_eq!(empty.consumed(), 0);

    let block = parse_cx_extensions("||suffix").unwrap();
    assert!(block.records().is_empty());
    assert_eq!(block.consumed(), 2);

    let utf8 = parse_cx_extensions("|$é$|尾").unwrap();
    assert_eq!(utf8.consumed(), "|$é$|".len());
    assert_eq!(
        utf8.records(),
        &[CxRecord::AtomLabels(vec![Some("é".to_owned())])]
    );

    assert_error("u:0|", 0, "CXSMILES extension does not start with |");
    assert_error("|u:0", 4, "failure parsing CXSMILES extensions");
}

#[test]
fn coordinates_preserve_omissions_ordinals_threshold_and_extra_components() {
    let input = "|(;,,2;1,,;3,4,0.001;5,6,0.0011,99)(7,8)|tail";
    let parsed = parse_cx_extensions(input).unwrap();
    assert_eq!(parsed.consumed(), input.find("tail").unwrap());
    assert_eq!(parsed.records().len(), 2);

    let CxRecord::Coordinates(first) = &parsed.records()[0] else {
        panic!("expected first coordinate record");
    };
    assert_eq!(first.conformer, 0);
    assert_eq!(
        first.values,
        vec![
            None,
            Some([0.0, 0.0, 2.0]),
            Some([1.0, 0.0, 0.0]),
            Some([3.0, 4.0, 0.001]),
            Some([5.0, 6.0, 0.0011]),
        ]
    );
    assert!(first.is_3d);

    let CxRecord::Coordinates(second) = &parsed.records()[1] else {
        panic!("expected second coordinate record");
    };
    assert_eq!(second.conformer, 1);
    assert_eq!(second.values, vec![Some([7.0, 8.0, 0.0])]);
    assert!(!second.is_3d);
}

#[test]
fn coordinate_errors_report_exact_source_byte_positions() {
    assert_error("|(x,0)|", 2, "invalid CX coordinate");
    assert_error("|(0,0", 1, "unterminated CX coordinate record");
}

#[test]
fn labels_and_values_preserve_slots_whitespace_and_decoded_entities() {
    let parsed = parse_cx_extensions("|$a;;b&#44;c; space $,$_AV:x&#58;y;;z$|suffix").unwrap();
    assert_eq!(
        parsed.records(),
        &[
            CxRecord::AtomLabels(vec![
                Some("a".to_owned()),
                None,
                Some("b,c".to_owned()),
                Some(" space ".to_owned()),
            ]),
            CxRecord::Unknown(",".to_owned()),
            CxRecord::AtomValues(vec![Some("x:y".to_owned()), None, Some("z".to_owned()),]),
        ]
    );
}

#[test]
fn atom_properties_keep_two_structural_dots_and_value_dot_content() {
    let parsed =
        parse_cx_extensions("|atomProp:0.na&#46;me.value.with.dots&#44;x:1.empty.:2.k.v,u:3|")
            .unwrap();
    let CxRecord::AtomProperties(properties) = &parsed.records()[0] else {
        panic!("expected atom properties");
    };
    assert_eq!(properties.len(), 2);
    assert_eq!(properties[0].atom, 0);
    assert_eq!(properties[0].name, "na.me");
    assert_eq!(properties[0].value, "value.with.dots,x");
    assert_eq!(properties[1].atom, 2);
    assert_eq!(properties[1].name, "k");
    assert_eq!(properties[1].value, "v");
    assert!(matches!(parsed.records()[1], CxRecord::Unsaturation(ref rows) if rows == &[3]));

    let empty_name = parse_cx_extensions("|atomProp:0..ignored|").unwrap();
    assert!(matches!(
        empty_name.records(),
        [CxRecord::AtomProperties(rows)] if rows.is_empty()
    ));
    assert_error(
        "|atomProp:0name.value|",
        11,
        "expected '.', found CX syntax mismatch",
    );
}

#[test]
fn coordinate_hydrogen_and_zero_bond_records_preserve_order() {
    let parsed = parse_cx_extensions("|C:1.2,3.4,H:5.6,Z:7,8|").unwrap();
    let CxRecord::CoordinateBonds(dative) = &parsed.records()[0] else {
        panic!("expected dative bonds");
    };
    assert_eq!(dative.kind, CxCoordinateBondKind::Dative);
    assert_eq!(
        dative
            .bonds
            .iter()
            .map(|row| (row.atom, row.bond))
            .collect::<Vec<_>>(),
        vec![(1, 2), (3, 4)]
    );
    let CxRecord::CoordinateBonds(hydrogen) = &parsed.records()[1] else {
        panic!("expected hydrogen bonds");
    };
    assert_eq!(hydrogen.kind, CxCoordinateBondKind::Hydrogen);
    assert_eq!(
        hydrogen
            .bonds
            .iter()
            .map(|row| (row.atom, row.bond))
            .collect::<Vec<_>>(),
        vec![(5, 6)]
    );
    assert!(matches!(parsed.records()[2], CxRecord::ZeroBonds(ref rows) if rows == &[7, 8]));
    assert_error("|C:1.|", 5, "invalid CX integer");
}

#[test]
fn radicals_cover_all_source_marker_classes_and_repeated_sections() {
    let parsed = parse_cx_extensions("|^1:0,^2:1^3:2^4:3^5:4^6:5^7:6|").unwrap();
    let CxRecord::Radicals(rows) = &parsed.records()[0] else {
        panic!("expected radicals");
    };
    assert_eq!(
        rows.iter()
            .map(|row| (row.atom, row.electrons))
            .collect::<Vec<_>>(),
        vec![(0, 1), (1, 2), (2, 2), (3, 2), (4, 3), (5, 3), (6, 3)]
    );
    assert_error("|^0:1|", 2, "invalid CX radical marker");
}

#[test]
fn enhanced_stereo_covers_absolute_or_and_ids_and_repetition() {
    let parsed = parse_cx_extensions("|a:0,1,o:2,o7:3,&:4,&9:5|").unwrap();
    let expected = [
        (CxStereoGroupKind::Absolute, 0, vec![0, 1]),
        (CxStereoGroupKind::Or, 0, vec![2]),
        (CxStereoGroupKind::Or, 7, vec![3]),
        (CxStereoGroupKind::And, 0, vec![4]),
        (CxStereoGroupKind::And, 9, vec![5]),
    ];
    assert_eq!(parsed.records().len(), expected.len());
    for (record, (kind, group_id, atoms)) in parsed.records().iter().zip(expected) {
        let CxRecord::EnhancedStereo(stereo) = record else {
            panic!("expected enhanced stereo");
        };
        assert_eq!(
            (stereo.kind, stereo.group_id, &stereo.atoms),
            (kind, group_id, &atoms)
        );
    }
}

#[test]
fn ring_unsaturation_and_substitution_constraints_cover_all_classifications() {
    let parsed = parse_cx_extensions("|rb:0:0,1:2,2:3,3:4,4:*,u:5,6,s:7:0,8:*|").unwrap();
    let CxRecord::RingBonds(rings) = &parsed.records()[0] else {
        panic!("expected ring constraints");
    };
    assert_eq!(
        rings
            .iter()
            .map(|row| (row.atom, row.constraint))
            .collect::<Vec<_>>(),
        vec![
            (0, CxCountConstraint::Exact(0)),
            (1, CxCountConstraint::Exact(2)),
            (2, CxCountConstraint::Exact(3)),
            (3, CxCountConstraint::LessEqual(4)),
            (4, CxCountConstraint::QueryScan),
        ]
    );
    assert!(matches!(parsed.records()[1], CxRecord::Unsaturation(ref rows) if rows == &[5, 6]));
    let CxRecord::Substitution(substitution) = &parsed.records()[2] else {
        panic!("expected substitution constraints");
    };
    assert_eq!(substitution[0].constraint, CxCountConstraint::Exact(0));
    assert_eq!(substitution[1].constraint, CxCountConstraint::QueryScan);
    assert_error("|rb:0:1|", 7, "unrecognized CX ring-bond count");
    assert_error("|rb:0:|", 6, "invalid CX integer");
}

#[test]
fn link_nodes_preserve_explicit_and_deferred_outer_atoms() {
    let parsed = parse_cx_extensions("|LN:1:1.3,4:1.4.3.6|").unwrap();
    let CxRecord::LinkNodes(rows) = &parsed.records()[0] else {
        panic!("expected link nodes");
    };
    assert_eq!(rows.len(), 2);
    assert_eq!(
        (
            rows[0].atom,
            rows[0].start_repetitions,
            rows[0].end_repetitions
        ),
        (1, 1, 3)
    );
    assert_eq!(rows[0].outer_atoms, None);
    assert_eq!(rows[1].outer_atoms, Some([3, 6]));
    assert_error("|LN:1:1.3.2|", 12, "invalid CX integer");
}

#[test]
fn data_and_hierarchy_sgroups_preserve_all_fields_and_order() {
    let parsed =
        parse_cx_extensions("|SgD:3,2,1,0:name:data&#58;x:like:unit:t:(1.,1.),SgH:1:0.2,3:4|")
            .unwrap();
    let CxRecord::DataSGroup(data) = &parsed.records()[0] else {
        panic!("expected data SGroup");
    };
    assert_eq!(data.atoms, vec![3, 2, 1, 0]);
    assert_eq!(data.field_name, "name");
    assert_eq!(data.data, "data:x");
    assert_eq!(data.query_op, "like");
    assert_eq!(data.field_info, "unit");
    assert_eq!(data.field_tag, "t");
    assert_eq!(data.coordinates.as_deref(), Some("(1.,1."));

    let CxRecord::Unknown(separator) = &parsed.records()[1] else {
        panic!("expected source-advanced separator byte");
    };
    assert_eq!(separator, ",");
    let CxRecord::SGroupHierarchy(hierarchy) = &parsed.records()[2] else {
        panic!("expected SGroup hierarchy");
    };
    assert_eq!(hierarchy.len(), 2);
    assert_eq!(
        (hierarchy[0].parent, &hierarchy[0].children),
        (1, &vec![0, 2])
    );
    assert_eq!((hierarchy[1].parent, &hierarchy[1].children), (3, &vec![4]));
    assert_error("|SgH:1:0,Sg:n:8|", 9, "invalid CX integer");

    let omitted = parse_cx_extensions("|SgD:2,1:FIELD:info::::|").unwrap();
    let CxRecord::DataSGroup(omitted) = &omitted.records()[0] else {
        panic!("expected data SGroup");
    };
    assert_eq!(omitted.coordinates, None);
    assert_error(
        "|SgD:0:a:b:c:d:e:(x|",
        20,
        "failure parsing CXSMILES extensions",
    );
}

#[test]
fn polymer_sgroups_cover_absent_partial_and_complete_optional_fields() {
    let absent = parse_cx_extensions("|Sg:n:4,1,2,3|").unwrap();
    let CxRecord::PolymerSGroup(absent) = &absent.records()[0] else {
        panic!("expected polymer SGroup");
    };
    assert_eq!(absent.atoms, vec![4, 1, 2, 3]);
    assert_eq!((absent.label.as_str(), absent.connect.as_str()), ("", ""));
    assert!(absent.head_crossings.is_empty());
    assert!(absent.tail_crossings.is_empty());

    let partial = parse_cx_extensions("|Sg:alt:0,1:n|").unwrap();
    let CxRecord::PolymerSGroup(partial) = &partial.records()[0] else {
        panic!("expected polymer SGroup");
    };
    assert_eq!(partial.type_code, "alt");
    assert_eq!(partial.label, "n");
    assert_eq!(partial.connect, "");

    let complete = parse_cx_extensions("|Sg:ran:0,1:n:ht:2,3:4,5|").unwrap();
    let CxRecord::PolymerSGroup(complete) = &complete.records()[0] else {
        panic!("expected polymer SGroup");
    };
    assert_eq!(complete.connect, "ht");
    assert_eq!(complete.head_crossings, vec![2, 3]);
    assert_eq!(complete.tail_crossings, vec![4, 5]);
}

#[test]
fn attachments_wedges_and_double_bonds_cover_every_source_form() {
    let parsed =
        parse_cx_extensions("|m:2:3.5.4,6:7,w:1.0,wU:2.1,wD:3.2,ctu:4,5,c:6,t:7,8|").unwrap();
    let CxRecord::VariableAttachments(attachments) = &parsed.records()[0] else {
        panic!("expected attachments");
    };
    assert_eq!(attachments[0].atom, 2);
    assert_eq!(attachments[0].endpoints, vec![3, 5, 4]);
    assert_eq!(attachments[1].endpoints, vec![7]);

    for (index, (direction, configuration)) in [
        (CxWedgeDirection::Unknown, 2),
        (CxWedgeDirection::BeginWedge, 1),
        (CxWedgeDirection::BeginDash, 3),
    ]
    .into_iter()
    .enumerate()
    {
        let CxRecord::WedgedBonds(rows) = &parsed.records()[index + 1] else {
            panic!("expected wedge record");
        };
        assert_eq!(
            (rows[0].direction, rows[0].configuration),
            (direction, configuration)
        );
    }

    for (index, (kind, bonds)) in [
        (CxDoubleBondStereoKind::Any, vec![4, 5]),
        (CxDoubleBondStereoKind::Cis, vec![6]),
        (CxDoubleBondStereoKind::Trans, vec![7, 8]),
    ]
    .into_iter()
    .enumerate()
    {
        let CxRecord::DoubleBondStereo(row) = &parsed.records()[index + 4] else {
            panic!("expected double-bond stereo record");
        };
        assert_eq!((row.stereo, &row.bonds), (kind, &bonds));
    }
}

#[test]
fn mixed_every_family_block_preserves_complete_dispatch_order() {
    let input = concat!(
        "|(0,0),$L$,$_AV:V$,atomProp:0.p.v,C:0.0,H:1.1,Z:2,^1:3,a:4,",
        "rb:5:2,LN:6:1.2,SgD:7:F:D:Q:I:T:,Sg:n:8,u:9,s:10:*,",
        "m:11:12,w:13.14,ctu:15,c:16,t:17,SgH:1:0|suffix"
    );
    let parsed = parse_cx_extensions(input).unwrap();
    assert_eq!(parsed.consumed(), input.find("suffix").unwrap());
    let typed_records = parsed
        .records()
        .iter()
        .filter(|record| !matches!(record, CxRecord::Unknown(_)))
        .collect::<Vec<_>>();
    assert_eq!(typed_records.len(), 21);
    assert!(
        parsed
            .records()
            .iter()
            .filter_map(|record| match record {
                CxRecord::Unknown(raw) => Some(raw.as_str()),
                _ => None,
            })
            .all(|raw| raw == ",")
    );
    assert!(matches!(typed_records[0], CxRecord::Coordinates(_)));
    assert!(matches!(typed_records[1], CxRecord::AtomLabels(_)));
    assert!(matches!(typed_records[2], CxRecord::AtomValues(_)));
    assert!(matches!(typed_records[3], CxRecord::AtomProperties(_)));
    assert!(matches!(typed_records[4], CxRecord::CoordinateBonds(_)));
    assert!(matches!(typed_records[5], CxRecord::CoordinateBonds(_)));
    assert!(matches!(typed_records[6], CxRecord::ZeroBonds(_)));
    assert!(matches!(typed_records[7], CxRecord::Radicals(_)));
    assert!(matches!(typed_records[8], CxRecord::EnhancedStereo(_)));
    assert!(matches!(typed_records[9], CxRecord::RingBonds(_)));
    assert!(matches!(typed_records[10], CxRecord::LinkNodes(_)));
    assert!(matches!(typed_records[11], CxRecord::DataSGroup(_)));
    assert!(matches!(typed_records[12], CxRecord::PolymerSGroup(_)));
    assert!(matches!(typed_records[13], CxRecord::Unsaturation(_)));
    assert!(matches!(typed_records[14], CxRecord::Substitution(_)));
    assert!(matches!(
        typed_records[15],
        CxRecord::VariableAttachments(_)
    ));
    assert!(matches!(typed_records[16], CxRecord::WedgedBonds(_)));
    assert!(matches!(typed_records[17], CxRecord::DoubleBondStereo(_)));
    assert!(matches!(typed_records[18], CxRecord::DoubleBondStereo(_)));
    assert!(matches!(typed_records[19], CxRecord::DoubleBondStereo(_)));
    assert!(matches!(typed_records[20], CxRecord::SGroupHierarchy(_)));
}

#[test]
fn unknown_records_are_lossless_interleaved_and_always_make_progress() {
    let input = "|###:α,rb:0:2,???,z,s:1:*|tail";
    let parsed = parse_cx_extensions(input).unwrap();
    assert_eq!(parsed.consumed(), input.find("tail").unwrap());
    assert_eq!(parsed.records().len(), 4);
    assert!(matches!(parsed.records()[0], CxRecord::Unknown(ref raw) if raw == "###:α,"));
    assert!(matches!(parsed.records()[1], CxRecord::RingBonds(_)));
    assert!(matches!(parsed.records()[2], CxRecord::Unknown(ref raw) if raw == "???,z,"));
    assert!(matches!(parsed.records()[3], CxRecord::Substitution(_)));
    assert_error(
        "|vendor:α,rb:0:2|",
        6,
        "expected ':', found CX syntax mismatch",
    );
}

#[test]
fn malformed_supported_syntax_reports_stable_utf8_byte_offsets_and_messages() {
    assert_error(
        "|$é&#12$|",
        4,
        "failure parsing CXSMILES extensions: quoted block not terminated with ';'",
    );
    assert_error("|^8:0|", 2, "invalid CX radical marker");
    assert_error("|wX:1.0|", 2, "invalid CX wedge marker");
}
