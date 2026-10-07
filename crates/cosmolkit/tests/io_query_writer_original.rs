#![cfg(feature = "cap-io")]
//! Source-corrected proposals: original d892 bytes and actual failures are immutable.
//! Only two expectations change under ROOT io44-pinned-query-writer-conflict; p1 review pending.
use cosmolkit::*;
use std::collections::BTreeMap;
fn mol_to_v2000_block(record: &SdfRecord) -> Result<String, MolecularIoError> {
    record.to_mol()
}
fn atom_list_query_molecule() -> SdfRecord {
    let mut first = QueryAtom::from_parts(
        Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::DUMMY)),
        QueryNode::predicate(AtomQueryPredicate::AtomicNumberIn(vec![6, 7])),
    );
    first.set_prop("molFileValue", "[#6,#7]").unwrap();
    let second = QueryAtom::from_parts(
        Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
        QueryNode::predicate(AtomQueryPredicate::AtomicNumberNotIn(vec![8, 16])),
    );
    let query = QueryGraph::from_parts(
        vec![first, second],
        vec![],
        BTreeMap::new(),
        vec![Conformer2D::new(0, vec![[0.0, 0.0], [1.0, 0.0]])],
        vec![],
        vec![],
    )
    .unwrap();
    SdfRecord::from_query_graph(query, MoleculeProperties::default().with_name("atom-list"))
        .unwrap()
}
fn query_bond_molecule() -> SdfRecord {
    let atoms = (0..10)
        .map(|i| {
            QueryAtom::from_carrier_parts(
                Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)),
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
            )
        })
        .collect();
    let definitions = vec![
        (
            0,
            1,
            BondOrder::Unspecified,
            QueryNode::predicate(BondQueryPredicate::Any),
        ),
        (
            2,
            3,
            BondOrder::Unspecified,
            QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                BondOrder::Single,
                BondOrder::Double,
            ])),
        ),
        (
            4,
            5,
            BondOrder::Unspecified,
            QueryNode::predicate(BondQueryPredicate::MolFileQueryCode(42)),
        ),
        (
            6,
            7,
            BondOrder::Single,
            QueryNode::predicate(BondQueryPredicate::IsInRing(true)),
        ),
        (
            8,
            9,
            BondOrder::Unspecified,
            QueryNode::and(vec![
                QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                    BondOrder::Single,
                    BondOrder::Aromatic,
                ])),
                QueryNode::predicate(BondQueryPredicate::IsInRing(false)),
            ]),
        ),
    ];
    let bonds = definitions
        .into_iter()
        .enumerate()
        .map(|(i, (a, z, order, predicate))| {
            QueryBond::from_parts(
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(a), AtomId::new(z), order),
                ),
                predicate,
            )
        })
        .collect();
    let query = QueryGraph::from_parts(
        atoms,
        bonds,
        BTreeMap::new(),
        vec![Conformer2D::new(
            0,
            (0..10).map(|i| [i as f64, 0.0]).collect(),
        )],
        vec![],
        vec![],
    )
    .unwrap();
    SdfRecord::from_query_graph(query, MoleculeProperties::default().with_name("query-bond"))
        .unwrap()
}

#[test]
fn mol_to_v2000_block_writes_atom_list_query_lines() {
    let molecule = atom_list_query_molecule();

    let block = mol_to_v2000_block(&molecule).unwrap();

    assert!(block.contains("    0.0000    0.0000    0.0000 L   0"));
    assert!(block.contains("V    1 [#6,#7]\n"));
    assert!(block.contains("M  ALS   1  2 F C   N   \n"));
    assert!(block.contains("M  ALS   2  2 T O   S   \n"));
}

#[test]
fn mol_to_v2000_block_writes_supported_query_bond_type_codes() {
    let molecule = query_bond_molecule();

    let block = mol_to_v2000_block(&molecule).unwrap();

    assert!(block.contains("  1  2  8  0\n"));
    assert!(block.contains("  3  4  5  0\n"));
    assert!(block.contains("  5  6  8  0\n"));
    assert!(block.contains("  7  8  1  0  0  1\n"));
    assert!(block.contains("  9 10  0  0  0  2\n"));
}
