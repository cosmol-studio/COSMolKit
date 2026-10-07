#![cfg(feature = "cap-io")]
//! Same-input source oracle and canonical detached writer proposals, p1 review pending.
use cosmolkit::*;
use std::collections::BTreeMap;
fn record(predicate: QueryNode<BondQueryPredicate>) -> SdfRecord {
    let atoms = (0..2)
        .map(|i| {
            QueryAtom::from_carrier_parts(
                Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)),
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
            )
        })
        .collect();
    let bond = QueryBond::from_parts(
        Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Unspecified),
        ),
        predicate,
    );
    let graph = QueryGraph::from_parts(
        atoms,
        vec![bond],
        BTreeMap::new(),
        vec![
            Conformer2D::new(0, vec![[1.25, -2.5], [-3.75, 4.0]]),
            Conformer2D::new(1, vec![[7.0, 8.0], [9.0, 10.0]]),
        ],
        vec![],
        vec![],
    )
    .unwrap();
    SdfRecord::from_query_graph(
        graph,
        MoleculeProperties::default()
            .with_name("query-boundary")
            .with_sdf_data_field("ID", "original"),
    )
    .unwrap()
}
fn params(format: SdfFormat) -> MolBlockWriteParams {
    MolBlockWriteParams {
        format,
        kekulize: false,
        include_stereo: false,
        coordinate_selection: MolCoordinateSelection::TwoD { id: 0 },
        ..Default::default()
    }
}
#[test]
fn io44_query_negated_order_and_ring_topology_matches_source_in_both_child_orders() {
    for in_ring in [true, false] {
        for reverse in [false, true] {
            let not_order = QueryNode::Not(Box::new(QueryNode::predicate(
                BondQueryPredicate::Order(BondOrder::Single),
            )));
            let ring = QueryNode::predicate(BondQueryPredicate::IsInRing(in_ring));
            let children = if reverse {
                vec![ring, not_order]
            } else {
                vec![not_order, ring]
            };
            let rec = record(QueryNode::and(children));
            let before = rec.query_graph().unwrap().clone();
            for format in [SdfFormat::V2000, SdfFormat::V3000] {
                let text = rec.to_mol_with_params(&params(format)).unwrap();
                let topo = if in_ring { 1 } else { 2 };
                let expected = if format == SdfFormat::V2000 {
                    format!("  1  2  0  0  0  {topo}\n")
                } else {
                    format!("M  V30 1 0 1 2 TOPO={topo}\n")
                };
                assert!(text.contains(&expected), "{text}");
            }
            assert_eq!(rec.query_graph().unwrap(), &before);
        }
    }
}
#[test]
fn io44_query_writer_preserves_detached_payload_and_data_fields_on_selection_success_and_failure() {
    let rec = record(QueryNode::predicate(BondQueryPredicate::Any));
    let before = rec.query_graph().unwrap().clone();
    assert!(rec.to_mol().is_err());
    for format in [SdfFormat::V2000, SdfFormat::V3000] {
        let text = rec.to_sdf_with_params(&params(format)).unwrap();
        assert!(text.starts_with("query-boundary\n"));
        assert!(text.ends_with(">  <ID>  \noriginal\n\n$$$$\n"));
        let reread = SdfRecord::from_sdf_with_params(
            &text,
            &SdfReadParams {
                sanitize: false,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(matches!(reread.graph(), SdfGraph::Query(_)));
        assert_eq!(
            reread.data_fields(),
            &[(PropertyText::from("ID"), PropertyText::from("original"))]
        );
        assert_eq!(
            reread.query_graph().unwrap().coordinates_2d().unwrap(),
            &[[1.25, -2.5], [-3.75, 4.0]]
        );
    }
    assert_eq!(rec.query_graph().unwrap(), &before);
}
