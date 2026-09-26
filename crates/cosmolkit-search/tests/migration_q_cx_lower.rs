use cosmolkit_cx::{
    CxAtomConstraint, CxCountConstraint, CxRecord, CxRingBond, ParsedCxExtensions,
    parse_cx_extensions,
};
use cosmolkit_model::{
    Atom, AtomId, AtomQueryPredicate, AtomSpec, BondId, Element, QueryAtom, QueryAtomIdentity,
    QueryNode, SGroupConnection, StereoGroup, StereoGroupKind, SubstanceGroupId,
    SubstanceGroupKind, query_substance_groups,
};
use cosmolkit_search::{QueryGraph, SmartsParseParams, apply_cx_to_query_graph, parse_smarts};

fn parse_query(smarts: &str) -> QueryGraph {
    parse_smarts(smarts, &SmartsParseParams::default())
        .unwrap_or_else(|error| panic!("pinned CXSMARTS {smarts:?}: {error}"))
}

#[test]
fn q26_atom_label_and_value_slots_keep_decoded_text_and_source_indices() {
    let labels = parse_query("CCC |$left&#59;semi;right\\raw;;ignored$|");
    assert_eq!(labels.num_atoms(), 3);
    assert_eq!(labels.atom(0).unwrap().prop("atomLabel"), Some("left;semi"));
    assert_eq!(
        labels.atom(1).unwrap().prop("atomLabel"),
        Some(r"right\raw")
    );
    assert_eq!(labels.atom(2).unwrap().prop("atomLabel"), None);
    assert!(
        labels
            .atoms()
            .iter()
            .all(|atom| { atom.identity() == QueryAtomIdentity::Element(Element::C) })
    );

    let values = parse_query("CCC |$_AV:first&#59;middle;second\\raw;;ignored$|");
    assert_eq!(
        values.atom(0).unwrap().prop("molFileValue"),
        Some("first;middle")
    );
    assert_eq!(
        values.atom(1).unwrap().prop("molFileValue"),
        Some(r"second\raw")
    );
    assert_eq!(values.atom(2).unwrap().prop("molFileValue"), None);
}

#[test]
fn q26_atom_properties_keep_decoded_values_overwrite_order_and_skip_invalid_indices() {
    let graph = parse_query(
        "CC |atomProp:0.escaped.first&#58;part:0.note.before:0.note.after:1.note.right\\raw:9.note.ignored|",
    );
    let first = graph.atom(0).unwrap();
    let second = graph.atom(1).unwrap();
    assert_eq!(first.prop("escaped"), Some("first:part"));
    assert_eq!(first.prop("note"), Some("after"));
    assert_eq!(second.prop("note"), Some(r"right\raw"));
}

#[test]
fn q26_labels_and_atom_properties_follow_record_order_before_label_processing() {
    let property_after_label = parse_query("C |$Q_e$atomProp:0.atomLabel.final|");
    let atom = property_after_label.atom(0).unwrap();
    assert_eq!(atom.identity(), QueryAtomIdentity::Element(Element::C));
    assert_eq!(atom.prop("atomLabel"), Some("final"));

    let label_after_property = parse_query("C |atomProp:0.atomLabel.before,$Q_e$|");
    let atom = label_after_property.atom(0).unwrap();
    assert_eq!(atom.identity(), QueryAtomIdentity::Element(Element::DUMMY));
    assert_eq!(atom.prop("atomLabel"), Some("Q_e"));
}

#[test]
fn q27_cx_query_predicates_preserve_origin_and_raw_identity() {
    let carrier_atom = QueryAtom::from_carrier_parts(
        Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
    );
    let explicit_atom = QueryAtom::from_identity_parts(
        AtomId::new(1),
        QueryAtomIdentity::AtomicNumber(246),
        QueryNode::and(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(246)),
            QueryNode::predicate(AtomQueryPredicate::FormalCharge(-1)),
        ]),
    );
    let null_query_atom = QueryAtom::from_identity_parts(
        AtomId::new(2),
        QueryAtomIdentity::AtomicNumber(247),
        QueryNode::predicate(AtomQueryPredicate::Any),
    );
    let mut graph = QueryGraph::from_parts(
        vec![carrier_atom, explicit_atom, null_query_atom],
        Vec::new(),
        Default::default(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("valid detached query graph");
    let carrier_before = graph.atom(0).unwrap().try_to_atom().unwrap();
    assert!(graph.atom(0).unwrap().predicate_is_carrier_derived());
    assert!(!graph.atom(1).unwrap().predicate_is_carrier_derived());

    let parsed = ParsedCxExtensions::new(
        vec![
            CxRecord::Unsaturation(vec![0]),
            CxRecord::RingBonds(vec![CxRingBond {
                atom: 1,
                constraint: CxCountConstraint::Exact(3),
            }]),
            CxRecord::Substitution(vec![CxAtomConstraint {
                atom: 2,
                constraint: CxCountConstraint::QueryScan,
            }]),
        ],
        0,
    );
    apply_cx_to_query_graph(&mut graph, &parsed).expect("source query predicate expansion");

    let carrier = graph.atom(0).unwrap();
    assert_eq!(carrier.identity(), QueryAtomIdentity::Element(Element::C));
    assert_eq!(carrier.try_to_atom().unwrap(), carrier_before);
    assert!(!carrier.predicate_is_carrier_derived());
    assert_eq!(
        carrier.predicate(),
        &QueryNode::and(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
            QueryNode::predicate(AtomQueryPredicate::IsUnsaturated),
        ])
    );

    let explicit = graph.atom(1).unwrap();
    assert_eq!(explicit.identity(), QueryAtomIdentity::AtomicNumber(246));
    assert!(!explicit.predicate_is_carrier_derived());
    assert_eq!(
        explicit.predicate(),
        &QueryNode::and(vec![
            QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(246)),
                QueryNode::predicate(AtomQueryPredicate::FormalCharge(-1)),
            ]),
            QueryNode::predicate(AtomQueryPredicate::RingBondCount(3)),
        ])
    );

    let null_query = graph.atom(2).unwrap();
    assert_eq!(null_query.identity(), QueryAtomIdentity::AtomicNumber(247));
    assert!(!null_query.predicate_is_carrier_derived());
    assert_eq!(
        null_query.predicate(),
        &QueryNode::predicate(AtomQueryPredicate::NonHydrogenDegree(0xDEAD_BEEF_u32))
    );
}

#[test]
fn q28_coordinate_projection_maps_slots_and_zero_fills_missing_rows() {
    let mut graph = parse_query("CCC");
    let parsed =
        parse_cx_extensions("|(1,2;;5,6)|").expect("pinned coordinate rows with an empty slot");
    apply_cx_to_query_graph(&mut graph, &parsed).expect("coordinate projection");

    assert_eq!(graph.conformers_3d().len(), 1);
    assert_eq!(
        graph.conformers_3d()[0].coordinates(),
        &[[1.0, 2.0, 0.0], [0.0, 0.0, 0.0], [5.0, 6.0, 0.0]]
    );
    assert!(!graph.conformers_3d()[0].is_3d());

    let mut missing_trailing_rows = parse_query("CCC");
    let parsed =
        parse_cx_extensions("|(4,5)|").expect("pinned coordinate rows with missing trailing slots");
    apply_cx_to_query_graph(&mut missing_trailing_rows, &parsed)
        .expect("missing coordinate rows are zero-filled");
    assert_eq!(
        missing_trailing_rows.conformers_3d()[0].coordinates(),
        &[[4.0, 5.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    );

    let mut excess_row = parse_query("C");
    let parsed =
        parse_cx_extensions("|(1,2;8,9,4)|").expect("pinned coordinate rows with an excess slot");
    apply_cx_to_query_graph(&mut excess_row, &parsed)
        .expect("invalid source atom rows are skipped");
    assert_eq!(
        excess_row.conformers_3d()[0].coordinates(),
        &[[1.0, 2.0, 0.0]]
    );
    assert!(!excess_row.conformers_3d()[0].is_3d());
}

#[test]
fn q28_coordinate_dimension_requires_projected_z_above_source_tolerance() {
    for (cx, expected_is_3d) in [
        ("|(0,0,0)|", false),
        ("|(0,0,0.001)|", false),
        ("|(0,0,0.0011)|", true),
        ("|(0,0,-0.0011)|", true),
    ] {
        let mut graph = parse_query("C");
        let parsed = parse_cx_extensions(cx).expect("pinned coordinate dimension record");
        apply_cx_to_query_graph(&mut graph, &parsed).expect("coordinate projection");
        assert_eq!(graph.conformers_3d()[0].is_3d(), expected_is_3d, "{cx}");
    }
}

#[test]
fn q29_enhanced_stereo_projects_group_kinds_ids_and_valid_members() {
    let mut graph = parse_query("CCO");
    let parsed = parse_cx_extensions("|a:0,99,a:0,o7:1,2,&3:2,o7:0,&3:1,o8:99|")
        .expect("pinned enhanced stereo records");
    apply_cx_to_query_graph(&mut graph, &parsed).expect("source enhanced stereo projection");

    assert_eq!(
        graph.stereo_groups(),
        &[
            StereoGroup::new(
                StereoGroupKind::Absolute,
                vec![AtomId::new(0), AtomId::new(0)],
                Vec::new(),
            )
            .with_id(0),
            StereoGroup::new(
                StereoGroupKind::Or,
                vec![AtomId::new(1), AtomId::new(2), AtomId::new(0)],
                Vec::new(),
            )
            .with_id(7),
            StereoGroup::new(
                StereoGroupKind::And,
                vec![AtomId::new(2), AtomId::new(1)],
                Vec::new(),
            )
            .with_id(3),
        ]
    );
}

#[test]
fn q30_data_sgroups_preserve_members_payload_and_source_indexes() {
    let mut graph = parse_query("CC");
    let parsed = parse_cx_extensions(
        "|SgD:9:IGNORED:::::SgD:1,0,1:FIELD:value,with,comma:=:unit:tag:(1,2)|",
    )
    .expect("pinned data SGroup fields and members");
    apply_cx_to_query_graph(&mut graph, &parsed).expect("source data SGroup projection");

    let groups = query_substance_groups(&graph);
    assert_eq!(groups.len(), 1);
    let group = &groups[0];
    assert_eq!(group.id(), SubstanceGroupId::new(0));
    assert_eq!(group.rdkit_sequence_id(), Some(1));
    assert_eq!(
        group.atoms(),
        &[AtomId::new(1), AtomId::new(0), AtomId::new(1)]
    );
    assert_eq!(
        group.props().get("_cxsmilesindex").map(String::as_str),
        Some("1")
    );
    assert_eq!(group.props().get("index").map(String::as_str), Some("1"));
    assert_eq!(
        group.props().get("FIELDNAME").map(String::as_str),
        Some("FIELD")
    );
    assert_eq!(
        group.props().get("DATAFIELDS").map(String::as_str),
        Some("value,with,comma")
    );
    assert_eq!(group.props().get("QUERYOP").map(String::as_str), Some("="));
    assert_eq!(
        group.props().get("FIELDINFO").map(String::as_str),
        Some("unit")
    );
    assert_eq!(
        group.props().get("FIELDTAG").map(String::as_str),
        Some("tag")
    );
    assert_eq!(
        group.props().get("COORDS").map(String::as_str),
        Some("(1,2")
    );
    assert_eq!(group.data_fields(), &["value,with,comma"]);

    let data = group.data().expect("typed data SGroup payload");
    assert_eq!(data.field_name.as_deref(), Some("FIELD"));
    assert_eq!(data.query_op.as_deref(), Some("="));
    assert_eq!(data.field_info.as_deref(), Some("unit"));
    assert_eq!(
        data.field_display.as_deref(),
        Some("    0.0000    0.0000    DR    ALL  0       0")
    );
    assert_eq!(data.values, ["value,with,comma"]);
}

#[test]
fn q31_polymer_sgroups_preserve_crossings_hierarchy_and_source_order() {
    let mut graph = parse_query("CCCC");
    let parsed =
        parse_cx_extensions("|Sg:n:1,0,1:parent:hh:2,0:1,2,Sg:n:3,2:child:ht:1,0:2,SgH:0:1|")
            .expect("pinned polymer and hierarchy records");
    apply_cx_to_query_graph(&mut graph, &parsed).expect("source polymer SGroup projection");

    let groups = query_substance_groups(&graph);
    assert_eq!(groups.len(), 2);

    let parent = &groups[0];
    assert_eq!(parent.id(), SubstanceGroupId::new(0));
    assert_eq!(parent.rdkit_sequence_id(), Some(0));
    assert_eq!(parent.kind(), &SubstanceGroupKind::StructuralRepeatUnit);
    assert_eq!(parent.label(), Some("parent"));
    assert_eq!(
        parent.atoms(),
        &[AtomId::new(1), AtomId::new(0), AtomId::new(1)]
    );
    assert_eq!(
        parent.bonds(),
        &[
            BondId::new(2),
            BondId::new(0),
            BondId::new(1),
            BondId::new(2),
        ]
    );
    assert_eq!(
        parent.head_crossing_bonds(),
        &[BondId::new(2), BondId::new(0)]
    );
    assert_eq!(
        parent.crossing_bond_correspondence(),
        &[
            BondId::new(2),
            BondId::new(1),
            BondId::new(0),
            BondId::new(2),
        ]
    );
    assert_eq!(parent.connection(), Some(&SGroupConnection::HeadToHead));
    assert_eq!(parent.parent(), None);
    assert_eq!(
        parent.props().get("_cxsmilesindex").map(String::as_str),
        Some("0")
    );
    assert_eq!(parent.props().get("index").map(String::as_str), Some("1"));

    let child = &groups[1];
    assert_eq!(child.id(), SubstanceGroupId::new(1));
    assert_eq!(child.rdkit_sequence_id(), Some(1));
    assert_eq!(child.kind(), &SubstanceGroupKind::StructuralRepeatUnit);
    assert_eq!(child.label(), Some("child"));
    assert_eq!(child.atoms(), &[AtomId::new(3), AtomId::new(2)]);
    assert_eq!(
        child.bonds(),
        &[BondId::new(1), BondId::new(0), BondId::new(2),]
    );
    assert_eq!(
        child.head_crossing_bonds(),
        &[BondId::new(1), BondId::new(0)]
    );
    assert_eq!(
        child.crossing_bond_correspondence(),
        &[BondId::new(1), BondId::new(2)]
    );
    assert_eq!(child.connection(), Some(&SGroupConnection::HeadToTail));
    assert_eq!(child.parent(), Some(SubstanceGroupId::new(0)));
    assert_eq!(
        child.props().get("_cxsmilesindex").map(String::as_str),
        Some("1")
    );
    assert_eq!(child.props().get("index").map(String::as_str), Some("2"));
    assert_eq!(child.props().get("PARENT").map(String::as_str), Some("1"));
}

#[test]
fn q32_link_node_projection_keeps_source_atom_references_without_topology_rewrite() {
    let graph = parse_query("CCC |LN:1:1.3,0:2.4.8.9,99:1.3|");

    assert_eq!(graph.num_atoms(), 3);
    assert_eq!(graph.num_bonds(), 2);
    assert_eq!(
        graph.prop("molFileLinkNodes"),
        Some("1 3 2 2 1 2 3|2 4 2 1 9 1 10")
    );
}

#[test]
fn q32_variable_attachment_projection_keeps_source_bond_and_endpoint_references() {
    let graph = parse_query("CCC |m:0:1.1.99|");

    assert_eq!(graph.num_atoms(), 3);
    assert_eq!(graph.num_bonds(), 2);
    let attached_bond = graph.bond(0).expect("anchor atom's single incident bond");
    assert_eq!(attached_bond.bond().begin(), AtomId::new(0));
    assert_eq!(attached_bond.bond().end(), AtomId::new(1));
    assert_eq!(
        attached_bond.bond().prop("_MolFileBondEndPts"),
        Some("(2 2 2)")
    );
    assert_eq!(attached_bond.bond().prop("_MolFileBondAttach"), Some("ANY"));
    assert_eq!(
        graph
            .bond(1)
            .expect("unrelated bond remains")
            .bond()
            .prop("_MolFileBondEndPts"),
        None
    );
}

#[test]
fn q32_variable_attachment_rejects_a_valid_anchor_with_nonunit_degree() {
    let error = parse_smarts("CCC |m:1:0|", &SmartsParseParams::default())
        .expect_err("pinned m: lowering rejects an anchor with degree two");

    assert!(
        error
            .to_string()
            .contains("position variation bond to atom with more than one bond")
    );
}
