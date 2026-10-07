//! Actual SMARTSQ counted-byte dependency through pinned processSMARTSQ.
use cosmolkit_io::{MolBlockRecord, MolPostParams, finish_mol_block_record};
use cosmolkit_model::{
    Atom, AtomId, AtomQueryPredicate, AtomSpec, CoordinateBlock, MoleculeProperties, PropertyText,
    PropertyValue, QueryNode, SGroupData, SubstanceGroup, SubstanceGroupId, SubstanceGroupKind,
    TopologyBlock, query_substance_groups,
};
use cosmolkit_types::Element;

fn source_record(query_type: &str, smarts: &[u8]) -> MolBlockRecord {
    let group = SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
        .with_atoms(vec![AtomId::new(0)])
        .with_data(SGroupData {
            query_type: Some(query_type.into()),
            values: vec![PropertyText::from_bytes(smarts), "[#8]".into()],
            ..SGroupData::default()
        });
    MolBlockRecord::Concrete {
        topology: TopologyBlock::try_from_parts(
            vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C))],
            vec![],
            vec![group],
            vec![],
        )
        .expect("fixed SMARTSQ producer topology is valid"),
        coordinates: CoordinateBlock::default(),
        properties: MoleculeProperties::default(),
    }
}

fn finish(original: &MolBlockRecord) -> cosmolkit_io::QueryMolBlockRecord {
    let MolBlockRecord::Query(result) = finish_mol_block_record(
        original.clone(),
        false,
        MolPostParams {
            sanitize: false,
            remove_hs: false,
            expand_attachment_points: false,
        },
        None,
    )
    .expect("source SMARTSQ consumes counted bytes through existing SEARCH") else {
        panic!("source SMARTSQ must promote its concrete input");
    };
    result
}

#[test]
fn smartsq_preserves_counted_recursive_source_and_original_producer() {
    let raw = b"CC source\xff\0\x80";
    let original = source_record("SMARTSQ", raw);
    let result = finish(&original);
    let atom = result.query.atom(0).expect("original atom retained");
    let QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(recursive)) = atom.predicate()
    else {
        panic!("source multiple-atom SMARTS creates the existing recursive query");
    };
    assert_eq!(
        recursive.source_smarts().map(PropertyText::as_bytes),
        Some(raw.as_slice())
    );
    assert_eq!(recursive.serial_number(), 0);
    let inner = recursive.query_graph().expect("recursive graph is owned");
    assert_eq!((inner.num_atoms(), inner.num_bonds()), (2, 1));
    assert_eq!(
        inner.name().unwrap().map(PropertyText::as_bytes),
        Some(b"source\xff\0\x80".as_slice())
    );
    assert_eq!(
        atom.prop("MRV SMA")
            .unwrap()
            .as_string()
            .unwrap()
            .as_bytes(),
        raw
    );
    assert_eq!(atom.prop("_MolFileAtomQuery"), Some(&PropertyValue::Int(1)));
    assert!(query_substance_groups(&result.query).is_empty());
    let MolBlockRecord::Concrete { topology, .. } = &original else {
        panic!("original stays concrete");
    };
    assert_eq!(topology.substance_groups.len(), 1);
    assert_eq!(
        topology.substance_groups[0].data().unwrap().values[0].as_bytes(),
        raw
    );
    assert_eq!(
        topology.substance_groups[0].data().unwrap().values[1].as_bytes(),
        b"[#8]"
    );
}

#[test]
fn sq_preserves_single_atom_predicate_and_exact_counted_mrv_sma() {
    let raw = b"[#7] source\xff\0\x80";
    let original = source_record("SQ", raw);
    let result = finish(&original);
    let atom = result.query.atom(0).expect("original atom retained");
    assert_eq!(
        atom.predicate(),
        &QueryNode::predicate(AtomQueryPredicate::AtomicNumber(7))
    );
    assert_eq!(
        atom.prop("MRV SMA")
            .unwrap()
            .as_string()
            .unwrap()
            .as_bytes(),
        raw
    );
    assert_eq!(atom.prop("_MolFileAtomQuery"), Some(&PropertyValue::Int(1)));
    assert!(query_substance_groups(&result.query).is_empty());
}
