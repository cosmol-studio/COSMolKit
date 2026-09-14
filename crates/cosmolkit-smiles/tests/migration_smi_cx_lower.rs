use cosmolkit_model::{
    AtomId, BondDirection, BondId, BondOrder, BondStereo, CoordinateDimension, SGroupBondRole,
    SGroupConnection, StereoGroupKind, SubstanceGroupKind,
};
use cosmolkit_smiles::{SmilesParseError, SmilesParseParams, parse_smiles};

fn parse(input: &str) -> cosmolkit_smiles::SmilesRecord {
    parse_smiles(input, &SmilesParseParams::default()).unwrap_or_else(|error| {
        panic!("failed to parse {input:?}: {error}");
    })
}

#[test]
fn coordinates_lower_into_dimension_specific_typed_conformers() {
    let record = parse("CCC |(0,0;1,2,0.001;3,4)(5,6,0.0011)|");
    assert_eq!(record.coordinates.conformers_2d.len(), 1);
    assert_eq!(record.coordinates.conformers_3d.len(), 1);
    assert_eq!(
        record.coordinates.source_coordinate_dim,
        Some(CoordinateDimension::ThreeD)
    );
    assert_eq!(record.coordinates.conformers_2d[0].id(), 0);
    assert_eq!(
        record.coordinates.conformers_2d[0].coordinates(),
        &[[0.0, 0.0], [1.0, 2.0], [3.0, 4.0]]
    );
    assert_eq!(record.coordinates.conformers_3d[0].id(), 1);
    assert_eq!(
        record.coordinates.conformers_3d[0].coordinates(),
        &[[5.0, 6.0, 0.0011], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    );

    let truncated = parse("CC |(1,2;3,4;5,6)|");
    assert_eq!(
        truncated.coordinates.conformers_2d[0].coordinates(),
        &[[1.0, 2.0], [3.0, 4.0]]
    );
    assert_eq!(
        truncated.coordinates.source_coordinate_dim,
        Some(CoordinateDimension::TwoD)
    );

    let error = parse_smiles("C |(1e309,0)|", &Default::default()).unwrap_err();
    assert!(matches!(error, SmilesParseError::Model(message) if message.contains("non-finite")));
}

#[test]
fn labels_values_and_properties_preserve_source_order_and_special_label_policy() {
    let record = parse(concat!(
        "[*][*]C |$Pol_p;Mod_p;ordinary;ignored$,$_AV:first;;third;ignored$,",
        "atomProp:0.key.old:0.key.new:2.kind.value:9.skip.value|"
    ));
    assert_eq!(record.topology.atoms[0].prop("atomLabel"), None);
    assert_eq!(record.topology.atoms[0].prop("dummyLabel"), Some("Pol"));
    assert_eq!(record.topology.atoms[1].prop("atomLabel"), None);
    assert_eq!(record.topology.atoms[1].prop("dummyLabel"), Some("Mod"));
    assert_eq!(record.topology.atoms[2].prop("atomLabel"), Some("ordinary"));
    assert_eq!(record.topology.atoms[0].prop("molFileValue"), Some("first"));
    assert_eq!(record.topology.atoms[1].prop("molFileValue"), None);
    assert_eq!(record.topology.atoms[2].prop("molFileValue"), Some("third"));
    assert_eq!(record.topology.atoms[0].prop("key"), Some("new"));
    assert_eq!(record.topology.atoms[2].prop("kind"), Some("value"));

    for label in [
        "star_e", "Q_e", "QH_p", "AH_p", "X_p", "XH_p", "M_p", "MH_p",
    ] {
        let input = format!("* |${label}$|");
        let error = parse_smiles(&input, &Default::default()).unwrap_err();
        assert!(
            matches!(error, SmilesParseError::UnsupportedCx(_)),
            "{label}: {error:?}"
        );
    }
}

#[test]
fn coordinate_hydrogen_and_zero_bonds_validate_ranges_and_orientation() {
    let record = parse("NOC |C:1.0,H:2.1,Z:9,0|");
    assert_eq!(record.topology.bonds[0].order(), BondOrder::Zero);
    assert_eq!(record.topology.bonds[1].order(), BondOrder::Hydrogen);
    assert_eq!(record.topology.bonds[1].begin(), AtomId::new(2));
    assert_eq!(record.topology.bonds[1].end(), AtomId::new(1));

    let ignored = parse("CC |C:8.0,C:0.8,Z:8|");
    assert_eq!(ignored.topology.bonds[0].order(), BondOrder::Single);

    let error = parse_smiles("CCC |C:2.0|", &Default::default()).unwrap_err();
    assert!(matches!(error, SmilesParseError::Cx(_)));
}

#[test]
fn enhanced_stereo_merges_repeated_groups_and_filters_invalid_atoms() {
    let record = parse("CCCC |a:0,9,o:1,o:2,9,&7:3,&8:9|");
    assert_eq!(record.topology.stereo_groups.len(), 3);
    assert_eq!(
        record.topology.stereo_groups[0].kind(),
        StereoGroupKind::Absolute
    );
    assert_eq!(record.topology.stereo_groups[0].id(), Some(0));
    assert_eq!(record.topology.stereo_groups[0].atoms(), &[AtomId::new(0)]);
    assert_eq!(record.topology.stereo_groups[1].kind(), StereoGroupKind::Or);
    assert_eq!(record.topology.stereo_groups[1].id(), Some(0));
    assert_eq!(
        record.topology.stereo_groups[1].atoms(),
        &[AtomId::new(1), AtomId::new(2)]
    );
    assert_eq!(
        record.topology.stereo_groups[2].kind(),
        StereoGroupKind::And
    );
    assert_eq!(record.topology.stereo_groups[2].id(), Some(7));
    assert_eq!(record.topology.stereo_groups[2].atoms(), &[AtomId::new(3)]);
}

#[test]
fn wedges_cover_all_directions_duplicate_rejection_and_cleanup() {
    for (record_text, direction, cfg) in [
        ("w:1.0", BondDirection::Unknown, "2"),
        ("wU:1.0", BondDirection::BeginWedge, "1"),
        ("wD:1.0", BondDirection::BeginDash, "3"),
    ] {
        let input = format!("CC |{record_text}|");
        let record = parse(&input);
        assert_eq!(record.topology.bonds[0].begin(), AtomId::new(1));
        assert_eq!(record.topology.bonds[0].direction(), direction);
        assert_eq!(record.topology.bonds[0].prop("_MolFileBondCfg"), Some(cfg));
        assert_eq!(record.properties.prop("_needsDetectAtomStereo"), None);
    }

    let error = parse_smiles("CC |wU:0.0,wD:0.0|", &Default::default()).unwrap_err();
    assert!(matches!(error, SmilesParseError::Cx(_)));

    let terminal = parse("C=C |wU:0.0|");
    assert_eq!(terminal.topology.bonds[0].order(), BondOrder::Double);
    assert_eq!(
        terminal.topology.bonds[0].prop("_MolFileBondCfg"),
        Some("1")
    );
    assert_eq!(terminal.properties.prop("_needsDetectAtomStereo"), None);
}

#[test]
fn double_bond_stereo_uses_source_control_order_and_terminal_noop() {
    for (marker, expected) in [
        ("ctu", BondStereo::Any),
        ("c", BondStereo::Cis),
        ("t", BondStereo::Trans),
    ] {
        let input = format!("CC=CC |{marker}:1|");
        let record = parse(&input);
        assert_eq!(record.topology.bonds[1].stereo(), expected);
        assert_eq!(
            record.topology.bonds[1].stereo_atoms(),
            Some([AtomId::new(0), AtomId::new(3)])
        );
    }

    let terminal = parse("C=C |c:0|");
    assert_eq!(terminal.topology.bonds[0].stereo(), BondStereo::None);
    assert_eq!(terminal.topology.bonds[0].stereo_atoms(), None);
    let ignored = parse("CC=CC |t:9|");
    assert_eq!(ignored.topology.bonds[1].stereo(), BondStereo::None);
}

#[test]
fn radicals_link_nodes_and_variable_attachments_preserve_order_and_errors() {
    let radicals = parse("CCCCCCC |^1:0,^2:1,^3:2,^4:3,^5:4,^6:5,^7:6,^7:9|");
    assert_eq!(
        radicals
            .topology
            .atoms
            .iter()
            .map(|atom| atom.radical_electrons())
            .collect::<Vec<_>>(),
        vec![1, 2, 2, 2, 3, 3, 3]
    );

    let link = parse("C1CC1.CCC |LN:1:1.3,4:2.5.3.5|");
    assert_eq!(
        link.properties.prop("_MolFileLinkNodes"),
        Some("1 3 2 2 1 2 3|2 5 2 5 4 5 6")
    );
    assert!(matches!(
        parse_smiles("CC |LN:0:1.2|", &Default::default()),
        Err(SmilesParseError::Cx(_))
    ));

    let attachment = parse("CO*.C1=CC=NC=C1 |m:2:3.99.5.4,m:99:0|");
    assert_eq!(
        attachment.topology.bonds[1].prop("_MolFileBondEndPts"),
        Some("(3 4 6 5)")
    );
    assert_eq!(
        attachment.topology.bonds[1].prop("_MolFileBondAttach"),
        Some("ANY")
    );
    assert!(matches!(
        parse_smiles("CCC |m:1:0.2|", &Default::default()),
        Err(SmilesParseError::Cx(_))
    ));
}

#[test]
fn data_polymer_and_hierarchy_records_install_typed_sgroups() {
    let data = parse("CCC |SgD:2,9,1:FIELD:value:like:unit:tag:(1.,2.)|");
    let group = &data.topology.substance_groups[0];
    assert_eq!(group.kind(), &SubstanceGroupKind::Data);
    assert_eq!(group.atoms(), &[AtomId::new(2), AtomId::new(1)]);
    assert_eq!(group.data().unwrap().field_name.as_deref(), Some("FIELD"));
    assert_eq!(group.data().unwrap().field_info.as_deref(), Some("unit"));
    assert_eq!(group.data().unwrap().query_op.as_deref(), Some("like"));
    assert_eq!(group.data().unwrap().values, ["value"]);
    assert_eq!(group.data_fields(), &["value"]);

    let expected_kinds = [
        ("n", SubstanceGroupKind::StructuralRepeatUnit),
        ("mon", SubstanceGroupKind::Monomer),
        ("mer", SubstanceGroupKind::Mer),
        ("co", SubstanceGroupKind::Copolymer),
        ("xl", SubstanceGroupKind::Crosslink),
        ("mod", SubstanceGroupKind::Modification),
        ("mix", SubstanceGroupKind::MixtureComponent),
        ("f", SubstanceGroupKind::Formulation),
        ("any", SubstanceGroupKind::AnyPolymer),
        ("gen", SubstanceGroupKind::Generic("GEN".to_owned())),
        ("c", SubstanceGroupKind::Generic("COM".to_owned())),
        ("grf", SubstanceGroupKind::Graft),
        ("alt", SubstanceGroupKind::Copolymer),
        ("ran", SubstanceGroupKind::Copolymer),
        ("blk", SubstanceGroupKind::Copolymer),
    ];
    for (code, expected) in expected_kinds {
        let input = format!("CC |Sg:{code}:0,1:label:ht|");
        let record = parse(&input);
        let group = &record.topology.substance_groups[0];
        assert_eq!(group.kind(), &expected, "{code}");
        assert_eq!(group.label(), Some("label"), "{code}");
        assert_eq!(
            group.connection(),
            Some(&SGroupConnection::HeadToTail),
            "{code}"
        );
        match code {
            "alt" => assert_eq!(group.subtype(), Some("ALT")),
            "ran" => assert_eq!(group.subtype(), Some("RAN")),
            "blk" => assert_eq!(group.subtype(), Some("BLO")),
            _ => {}
        }
    }

    let hierarchy = parse(concat!(
        "CC |SgD:0:PARENT:p::::,SgD:1:CHILD:c::::,SgH:0:1,",
        "SgH:9:1|"
    ));
    assert_eq!(
        hierarchy.topology.substance_groups[1]
            .parent()
            .unwrap()
            .index(),
        0
    );
    assert_eq!(
        hierarchy.topology.substance_groups[1]
            .props()
            .get("PARENT")
            .map(String::as_str),
        Some("1")
    );
    assert!(
        hierarchy
            .topology
            .substance_groups
            .iter()
            .all(|group| group.props().get("_cxsmilesindex").is_none())
    );

    let missing_child = parse("CC |SgD:9:SKIP:x::::,SgD:0:PARENT:p::::,SgH:1:0|");
    assert_eq!(missing_child.topology.substance_groups.len(), 1);
    assert_eq!(missing_child.topology.substance_groups[0].parent(), None);

    assert!(matches!(
        parse_smiles(
            "CC |SgD:0:PARENT:p::::,SgD:1:CHILD:c::::,SgH:0:9|",
            &Default::default()
        ),
        Err(SmilesParseError::Cx(message))
            if message == "child id references non-existent SGroup"
    ));
}

#[test]
fn polymer_crossings_are_typed_and_invalid_explicit_crossings_drop_the_group() {
    let inferred = parse("CCCC |Sg:n:1,2::ht|");
    let group = &inferred.topology.substance_groups[0];
    assert_eq!(group.bonds(), &[BondId::new(0), BondId::new(2)]);
    assert_eq!(group.bond_role(BondId::new(0)), SGroupBondRole::Crossing);
    assert_eq!(group.bond_role(BondId::new(2)), SGroupBondRole::Crossing);

    let single_atom = parse("CCC |Sg:n:1::eu|");
    assert_eq!(
        single_atom.topology.substance_groups[0].bonds(),
        &[BondId::new(0), BondId::new(1)]
    );
    assert_eq!(
        single_atom.topology.substance_groups[0].connection(),
        Some(&SGroupConnection::Either)
    );

    let dropped = parse("CCC |Sg:n:1::ht:9:|");
    assert!(dropped.topology.substance_groups.is_empty());
}

#[test]
fn query_only_and_late_lowering_failures_are_atomic_in_both_parser_modes() {
    for query in ["u:0", "rb:0:0", "s:0:*"] {
        let input = format!("CC |$first;second$,{query}| named");
        assert!(matches!(
            parse_smiles(&input, &Default::default()),
            Err(SmilesParseError::UnsupportedCx(_))
        ));
        let recovered = parse_smiles(
            &input,
            &SmilesParseParams {
                strict_cxsmiles: false,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(recovered.topology.atoms[0].prop("atomLabel"), None);
        assert_eq!(recovered.topology.atoms[1].prop("atomLabel"), None);
        assert_eq!(recovered.properties.prop("_CXSMILES_Data"), Some(""));
        assert_eq!(recovered.properties.name(), None);
    }

    let late = "CCC |$first;second;third$,C:2.0| named";
    assert!(matches!(
        parse_smiles(late, &Default::default()),
        Err(SmilesParseError::Cx(_))
    ));
    let recovered = parse_smiles(
        late,
        &SmilesParseParams {
            strict_cxsmiles: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(
        recovered
            .topology
            .atoms
            .iter()
            .all(|atom| atom.prop("atomLabel").is_none())
    );
    assert_eq!(recovered.topology.bonds[0].order(), BondOrder::Single);
    assert_eq!(recovered.properties.prop("_CXSMILES_Data"), Some(""));
}

#[test]
fn unknown_records_are_ignored_and_parser_only_indices_are_removed() {
    let record = parse("CC |vendor:opaque,Sg:n:0::ht|");
    assert_eq!(record.topology.substance_groups.len(), 1);
    assert!(
        record
            .topology
            .bonds
            .iter()
            .all(|bond| bond.prop("_cxsmilesBondIdx").is_none())
    );
    assert!(
        record
            .topology
            .substance_groups
            .iter()
            .all(|group| group.props().get("_cxsmilesindex").is_none())
    );
}
