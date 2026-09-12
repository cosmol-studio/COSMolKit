use cosmolkit_model::{
    AtomId, BondId, SGroupAttachPoint, SGroupBondRole, SGroupBracket, SGroupBracketStyle,
    SGroupCState, SGroupConnection, SGroupData, SGroupDisplay, StereoGroup, StereoGroupKind,
    SubstanceGroup, SubstanceGroupId, SubstanceGroupKind,
};

fn atom(index: usize) -> AtomId {
    AtomId::new(index)
}

fn bond(index: usize) -> BondId {
    BondId::new(index)
}

fn sgroup(index: usize) -> SubstanceGroupId {
    SubstanceGroupId::new(index)
}

#[test]
fn value_vocabularies_cover_every_frozen_branch() {
    let kinds = [
        SubstanceGroupKind::Data,
        SubstanceGroupKind::Superatom,
        SubstanceGroupKind::MultipleGroup,
        SubstanceGroupKind::StructuralRepeatUnit,
        SubstanceGroupKind::Monomer,
        SubstanceGroupKind::Copolymer,
        SubstanceGroupKind::Crosslink,
        SubstanceGroupKind::Graft,
        SubstanceGroupKind::Modification,
        SubstanceGroupKind::Mer,
        SubstanceGroupKind::AnyPolymer,
        SubstanceGroupKind::MixtureComponent,
        SubstanceGroupKind::Mixture,
        SubstanceGroupKind::Formulation,
        SubstanceGroupKind::Generic("VENDOR".into()),
    ];
    assert_eq!(kinds.len(), 15);
    assert_eq!(kinds[14], SubstanceGroupKind::Generic("VENDOR".into()));

    assert_eq!(
        [SGroupBondRole::Crossing, SGroupBondRole::Contained],
        [SGroupBondRole::Crossing, SGroupBondRole::Contained]
    );
    assert_eq!(
        [
            SGroupConnection::HeadToHead,
            SGroupConnection::HeadToTail,
            SGroupConnection::Either,
            SGroupConnection::Unknown("vendor".into()),
        ]
        .len(),
        4
    );
    assert_eq!(
        [
            SGroupBracketStyle::Bracket,
            SGroupBracketStyle::Parenthesis,
            SGroupBracketStyle::None,
            SGroupBracketStyle::Unknown("vendor".into()),
        ]
        .len(),
        4
    );
    assert_eq!(
        [
            StereoGroupKind::Absolute,
            StereoGroupKind::Or,
            StereoGroupKind::And,
        ],
        [
            StereoGroupKind::Absolute,
            StereoGroupKind::Or,
            StereoGroupKind::And,
        ]
    );
}

#[test]
fn substance_group_defaults_and_builders_preserve_all_typed_state() {
    let empty = SubstanceGroup::new(sgroup(0), SubstanceGroupKind::Data);
    assert_eq!(empty.id(), sgroup(0));
    assert_eq!(empty.kind(), &SubstanceGroupKind::Data);
    assert_eq!(empty.rdkit_sequence_id(), None);
    assert_eq!(empty.external_id(), None);
    assert!(empty.atoms().is_empty());
    assert!(empty.bonds().is_empty());
    assert!(empty.parent_atoms().is_empty());
    assert_eq!(empty.parent(), None);
    assert_eq!(empty.label(), None);
    assert_eq!(empty.connection(), None);
    assert_eq!(empty.subtype(), None);
    assert_eq!(empty.bracket_style(), None);
    assert_eq!(empty.expansion_state(), None);
    assert_eq!(empty.class(), None);
    assert_eq!(empty.component_number(), None);
    assert_eq!(empty.display(), None);
    assert_eq!(empty.data(), None);
    assert!(empty.attach_points().is_empty());
    assert!(empty.cstates().is_empty());
    assert!(empty.props().is_empty());
    assert!(empty.data_fields().is_empty());

    let display = SGroupDisplay {
        brackets: vec![SGroupBracket {
            p1: [1.0, 2.0],
            p2: [3.0, 4.0],
        }],
        field_position: Some([5.0, 6.0]),
        display_tag: Some("DA".into()),
    };
    let data = SGroupData {
        field_name: Some("FIELD".into()),
        field_type: Some("TYPE".into()),
        field_info: Some("INFO".into()),
        field_display: Some("DISPLAY".into()),
        units: Some("UNIT".into()),
        query_type: Some("QT".into()),
        query_op: Some("QO".into()),
        values: vec!["first".into(), "second".into()],
    };
    let attach_point = SGroupAttachPoint {
        atom: atom(2),
        leaving_atom: Some(atom(3)),
        label: Some("AP".into()),
        order: Some(7),
    };
    let cstate = SGroupCState {
        bond: bond(4),
        vector: [8.0, 9.0],
    };
    let group = SubstanceGroup::new(sgroup(1), SubstanceGroupKind::Generic("X".into()))
        .with_rdkit_sequence_id(101)
        .with_external_id(202)
        .with_atoms(vec![atom(2), atom(1)])
        .with_bonds(vec![bond(4), bond(3)])
        .with_bond_role(bond(4), SGroupBondRole::Contained)
        .with_parent_atoms(vec![atom(1)])
        .with_parent(sgroup(0))
        .with_label("label")
        .with_connection(SGroupConnection::Unknown("CX".into()))
        .with_subtype("ALT")
        .with_bracket_style(SGroupBracketStyle::Unknown("BX".into()))
        .with_expansion_state("expanded")
        .with_class("CHEM")
        .with_component_number(3)
        .with_display(display.clone())
        .with_data(data.clone())
        .with_attach_points(vec![attach_point.clone()])
        .with_cstates(vec![cstate])
        .with_prop("vendor", "kept")
        .with_data_field("raw-1")
        .with_data_field("raw-2");

    assert_eq!(group.id(), sgroup(1));
    assert_eq!(group.rdkit_sequence_id(), Some(101));
    assert_eq!(group.external_id(), Some(202));
    assert_eq!(group.kind(), &SubstanceGroupKind::Generic("X".into()));
    assert_eq!(group.atoms(), &[atom(2), atom(1)]);
    assert_eq!(group.bonds(), &[bond(4), bond(3)]);
    assert_eq!(group.bond_role(bond(4)), SGroupBondRole::Contained);
    assert_eq!(group.bond_role(bond(3)), SGroupBondRole::Crossing);
    assert_eq!(group.parent_atoms(), &[atom(1)]);
    assert_eq!(group.parent(), Some(sgroup(0)));
    assert_eq!(group.label(), Some("label"));
    assert_eq!(
        group.connection(),
        Some(&SGroupConnection::Unknown("CX".into()))
    );
    assert_eq!(group.subtype(), Some("ALT"));
    assert_eq!(
        group.bracket_style(),
        Some(&SGroupBracketStyle::Unknown("BX".into()))
    );
    assert_eq!(group.expansion_state(), Some("expanded"));
    assert_eq!(group.class(), Some("CHEM"));
    assert_eq!(group.component_number(), Some(3));
    assert_eq!(group.display(), Some(&display));
    assert_eq!(group.data(), Some(&data));
    assert_eq!(group.attach_points(), &[attach_point]);
    assert_eq!(group.cstates(), &[cstate]);
    assert_eq!(
        group.props().get("vendor").map(String::as_str),
        Some("kept")
    );
    assert_eq!(group.data_fields(), &["raw-1", "raw-2"]);
}

#[test]
fn substance_group_mutators_preserve_order_and_remove_only_first_match() {
    let mut group = SubstanceGroup::new(sgroup(0), SubstanceGroupKind::Superatom);
    group.set_id(sgroup(4));
    group.set_rdkit_sequence_id(11);
    group.set_external_id(12);
    group.set_parent(sgroup(3));
    group.set_label("SUP");
    group.set_connection(SGroupConnection::HeadToTail);
    group.set_subtype("RAN");
    group.set_bracket_style(SGroupBracketStyle::Parenthesis);
    group.set_expansion_state("contracted");
    group.set_class("AA");
    group.set_component_number(8);
    group.push_atom(atom(1));
    group.push_atom(atom(2));
    group.push_atom(atom(1));
    group.push_parent_atom(atom(1));
    group.push_parent_atom(atom(2));
    group.push_parent_atom(atom(1));
    group.push_bond_with_role(bond(1), SGroupBondRole::Contained);
    group.push_bond(bond(2));
    group.push_bond(bond(1));
    group.push_data_field("a");
    group.push_data_field("b");
    group.push_attach_point(SGroupAttachPoint {
        atom: atom(1),
        leaving_atom: Some(atom(9)),
        label: None,
        order: None,
    });
    group.push_attach_point(SGroupAttachPoint {
        atom: atom(2),
        leaving_atom: Some(atom(9)),
        label: Some("two".into()),
        order: Some(2),
    });
    group.push_cstate(SGroupCState {
        bond: bond(2),
        vector: [2.0, 3.0],
    });
    group.display_mut().display_tag = Some("tag".into());
    group.data_mut().values.push("value".into());
    group.set_prop("keep", "yes");
    group.set_prop("clear", "me");
    group.clear_prop("clear");

    group.remove_atom(atom(1));
    group.remove_parent_atom(atom(1));
    group.remove_bond(bond(1));
    assert_eq!(group.atoms(), &[atom(2), atom(1)]);
    assert_eq!(group.parent_atoms(), &[atom(2), atom(1)]);
    assert_eq!(group.bonds(), &[bond(2), bond(1)]);
    assert_eq!(group.bond_role(bond(1)), SGroupBondRole::Contained);

    group.remove_bond(bond(1));
    assert_eq!(group.bonds(), &[bond(2)]);
    assert_eq!(group.bond_role(bond(1)), SGroupBondRole::Crossing);
    let unchanged = group
        .clone()
        .with_bond_role(bond(99), SGroupBondRole::Contained);
    assert_eq!(unchanged.bond_role(bond(99)), SGroupBondRole::Crossing);

    group.clear_attach_point_leaving_atom(atom(9));
    assert!(
        group
            .attach_points()
            .iter()
            .all(|point| point.leaving_atom.is_none())
    );
    assert_eq!(group.id(), sgroup(4));
    assert_eq!(group.rdkit_sequence_id(), Some(11));
    assert_eq!(group.external_id(), Some(12));
    assert_eq!(group.parent(), Some(sgroup(3)));
    assert_eq!(group.label(), Some("SUP"));
    assert_eq!(group.connection(), Some(&SGroupConnection::HeadToTail));
    assert_eq!(group.subtype(), Some("RAN"));
    assert_eq!(
        group.bracket_style(),
        Some(&SGroupBracketStyle::Parenthesis)
    );
    assert_eq!(group.expansion_state(), Some("contracted"));
    assert_eq!(group.class(), Some("AA"));
    assert_eq!(group.component_number(), Some(8));
    assert_eq!(group.data_fields(), &["a", "b"]);
    assert_eq!(
        group
            .display()
            .and_then(|value| value.display_tag.as_deref()),
        Some("tag")
    );
    assert_eq!(
        group.data().expect("data").values.as_slice(),
        &[String::from("value")]
    );
    assert_eq!(group.props().get("keep").map(String::as_str), Some("yes"));
    assert!(!group.props().contains_key("clear"));
}

#[test]
fn membership_predicates_include_every_reference_category() {
    let group = SubstanceGroup::new(sgroup(0), SubstanceGroupKind::Data)
        .with_atoms(vec![atom(1)])
        .with_parent_atoms(vec![atom(2)])
        .with_bonds(vec![bond(3)])
        .with_attach_points(vec![SGroupAttachPoint {
            atom: atom(4),
            leaving_atom: Some(atom(5)),
            label: None,
            order: None,
        }])
        .with_cstates(vec![SGroupCState {
            bond: bond(6),
            vector: [0.0, 0.0],
        }]);

    for referenced in [atom(1), atom(2), atom(4), atom(5)] {
        assert!(group.includes_atom(referenced));
    }
    assert!(!group.includes_atom(atom(0)));
    assert!(group.includes_bond(bond(3)));
    assert!(group.includes_bond(bond(6)));
    assert!(!group.includes_bond(bond(0)));
}

fn fully_populated_group() -> SubstanceGroup {
    SubstanceGroup::new(sgroup(0), SubstanceGroupKind::MultipleGroup)
        .with_rdkit_sequence_id(31)
        .with_external_id(32)
        .with_atoms(vec![atom(0), atom(2)])
        .with_bonds(vec![bond(0), bond(2)])
        .with_bond_role(bond(2), SGroupBondRole::Contained)
        .with_parent_atoms(vec![atom(1)])
        .with_parent(sgroup(1))
        .with_label("MUL")
        .with_connection(SGroupConnection::Either)
        .with_subtype("ALT")
        .with_bracket_style(SGroupBracketStyle::Bracket)
        .with_expansion_state("expanded")
        .with_class("CHEM")
        .with_component_number(2)
        .with_display(SGroupDisplay {
            brackets: vec![SGroupBracket {
                p1: [1.0, 1.5],
                p2: [2.0, 2.5],
            }],
            field_position: Some([3.0, 3.5]),
            display_tag: Some("tag".into()),
        })
        .with_data(SGroupData {
            field_name: Some("name".into()),
            field_type: Some("type".into()),
            field_info: Some("info".into()),
            field_display: Some("display".into()),
            units: Some("units".into()),
            query_type: Some("query-type".into()),
            query_op: Some("query-op".into()),
            values: vec!["one".into(), "two".into()],
        })
        .with_attach_points(vec![SGroupAttachPoint {
            atom: atom(2),
            leaving_atom: Some(atom(3)),
            label: Some("AP".into()),
            order: Some(4),
        }])
        .with_cstates(vec![SGroupCState {
            bond: bond(0),
            vector: [4.0, 5.0],
        }])
        .with_prop("vendor", "property")
        .with_data_field("raw")
}

#[test]
fn substance_group_remap_maps_every_reference_and_preserves_other_state() {
    let group = fully_populated_group();
    let atom_map = [
        Some(atom(10)),
        Some(atom(11)),
        Some(atom(12)),
        Some(atom(13)),
    ];
    let bond_map = [Some(bond(20)), Some(bond(21)), Some(bond(22))];
    let sgroup_map = [Some(sgroup(4)), Some(sgroup(5))];
    assert!(group.can_remap_without_parent(&atom_map, &bond_map));

    let remapped = group
        .remapped(sgroup(9), &atom_map, &bond_map, &sgroup_map)
        .expect("all references are mapped");
    let expected = SubstanceGroup::new(sgroup(9), SubstanceGroupKind::MultipleGroup)
        .with_rdkit_sequence_id(31)
        .with_external_id(32)
        .with_atoms(vec![atom(10), atom(12)])
        .with_bonds(vec![bond(20), bond(22)])
        .with_bond_role(bond(22), SGroupBondRole::Contained)
        .with_parent_atoms(vec![atom(11)])
        .with_parent(sgroup(5))
        .with_label("MUL")
        .with_connection(SGroupConnection::Either)
        .with_subtype("ALT")
        .with_bracket_style(SGroupBracketStyle::Bracket)
        .with_expansion_state("expanded")
        .with_class("CHEM")
        .with_component_number(2)
        .with_display(group.display().expect("display").clone())
        .with_data(group.data().expect("data").clone())
        .with_attach_points(vec![SGroupAttachPoint {
            atom: atom(12),
            leaving_atom: Some(atom(13)),
            label: Some("AP".into()),
            order: Some(4),
        }])
        .with_cstates(vec![SGroupCState {
            bond: bond(20),
            vector: [4.0, 5.0],
        }])
        .with_prop("vendor", "property")
        .with_data_field("raw");
    assert_eq!(remapped, expected);
}

#[test]
fn substance_group_remap_rejects_each_missing_reference_category() {
    let atom_map = [
        Some(atom(10)),
        Some(atom(11)),
        Some(atom(12)),
        Some(atom(13)),
    ];
    let bond_map = [Some(bond(20)), Some(bond(21)), Some(bond(22))];
    let sgroup_map = [Some(sgroup(4)), Some(sgroup(5))];

    let cases = [
        SubstanceGroup::new(sgroup(0), SubstanceGroupKind::Data).with_atoms(vec![atom(4)]),
        SubstanceGroup::new(sgroup(0), SubstanceGroupKind::Data).with_parent_atoms(vec![atom(4)]),
        SubstanceGroup::new(sgroup(0), SubstanceGroupKind::Data).with_attach_points(vec![
            SGroupAttachPoint {
                atom: atom(4),
                leaving_atom: None,
                label: None,
                order: None,
            },
        ]),
        SubstanceGroup::new(sgroup(0), SubstanceGroupKind::Data).with_attach_points(vec![
            SGroupAttachPoint {
                atom: atom(0),
                leaving_atom: Some(atom(4)),
                label: None,
                order: None,
            },
        ]),
        SubstanceGroup::new(sgroup(0), SubstanceGroupKind::Data).with_bonds(vec![bond(3)]),
        SubstanceGroup::new(sgroup(0), SubstanceGroupKind::Data)
            .with_bonds(vec![bond(3)])
            .with_bond_role(bond(3), SGroupBondRole::Contained),
        SubstanceGroup::new(sgroup(0), SubstanceGroupKind::Data).with_cstates(vec![SGroupCState {
            bond: bond(3),
            vector: [0.0, 0.0],
        }]),
    ];
    for (case_index, group) in cases.into_iter().enumerate() {
        assert!(
            !group.can_remap_without_parent(&atom_map, &bond_map),
            "case {case_index} must fail the local precheck"
        );
        assert_eq!(
            group.remapped(sgroup(9), &atom_map, &bond_map, &sgroup_map),
            None,
            "case {case_index} must fail all-or-none remapping"
        );
    }

    let missing_parent =
        SubstanceGroup::new(sgroup(0), SubstanceGroupKind::Data).with_parent(sgroup(2));
    assert!(missing_parent.can_remap_without_parent(&atom_map, &bond_map));
    assert_eq!(
        missing_parent.remapped(sgroup(9), &atom_map, &bond_map, &sgroup_map),
        None
    );
}

#[test]
fn stereo_group_mutation_and_remap_are_ordered_and_all_or_none() {
    let mut group = StereoGroup::new(
        StereoGroupKind::Or,
        vec![atom(0), atom(1), atom(0)],
        vec![bond(0), bond(1), bond(0)],
    )
    .with_id(17);
    group.push_atom(atom(2));
    group.push_bond(bond(2));
    group.remove_atom(atom(0));
    group.remove_bond(bond(0));
    assert_eq!(group.id(), Some(17));
    assert_eq!(group.kind(), StereoGroupKind::Or);
    assert_eq!(group.atoms(), &[atom(1), atom(0), atom(2)]);
    assert_eq!(group.bonds(), &[bond(1), bond(0), bond(2)]);
    assert!(!group.is_empty());
    assert!(StereoGroup::new(StereoGroupKind::Absolute, vec![], vec![]).is_empty());

    let atom_map = [Some(atom(10)), Some(atom(11)), Some(atom(12))];
    let bond_map = [Some(bond(20)), Some(bond(21)), Some(bond(22))];
    assert_eq!(
        group.remapped(&atom_map, &bond_map),
        Some(
            StereoGroup::new(
                StereoGroupKind::Or,
                vec![atom(11), atom(10), atom(12)],
                vec![bond(21), bond(20), bond(22)],
            )
            .with_id(17)
        )
    );
    assert_eq!(group.remapped(&atom_map[..2], &bond_map), None);
    assert_eq!(group.remapped(&atom_map, &bond_map[..2]), None);
}
