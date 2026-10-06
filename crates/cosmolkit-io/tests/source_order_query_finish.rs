use cosmolkit_io::{MolBlockRecord, MolPostParams, QueryMolBlockRecord, finish_mol_block_record};
use cosmolkit_model::{
    AtomId, AtomSpec, Conformer2D, Conformer3D, CoordinateDimension, CoordinateSourceConformer,
    Element, MoleculeProperties, QueryAtom, QueryGraph,
};

#[test]
fn query_finish_retains_actual_first_and_complete_source_order() {
    // Pinned ROMol::getConformer(-1) uses front(); finishMolProcessing
    // edits the same res and never clears or reorders its conformers.
    let mut query = QueryGraph::from_parts(
        vec![QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C))],
        vec![],
        Default::default(),
        vec![Conformer2D::new(9, vec![[3.0, 4.0]])],
        vec![
            Conformer3D::new(9, vec![[1.0, 2.0, -0.0]], false),
            Conformer3D::new(0, vec![[5.0, 6.0, 7.0]], true),
        ],
        vec![],
    )
    .unwrap();
    let expected = vec![
        CoordinateDimension::ThreeD,
        CoordinateDimension::TwoD,
        CoordinateDimension::ThreeD,
    ];
    query
        .set_source_conformer_order(Some(expected.clone()))
        .unwrap();
    let input = query.coordinate_block(None);
    input.validate_for_atom_count(1).unwrap();
    let record = MolBlockRecord::Query(QueryMolBlockRecord {
        query,
        properties: MoleculeProperties::default(),
        source_coordinate_dim: Some(CoordinateDimension::TwoD),
    });
    let output = finish_mol_block_record(
        record,
        false,
        MolPostParams {
            sanitize: false,
            remove_hs: false,
            expand_attachment_points: false,
        },
        None,
    )
    .unwrap();
    let MolBlockRecord::Query(output) = output else {
        panic!("query identity must survive");
    };
    let coordinates = output.query.coordinate_block(output.source_coordinate_dim);
    assert_eq!(coordinates.conformers_2d, input.conformers_2d);
    assert_eq!(coordinates.conformers_3d, input.conformers_3d);
    assert_eq!(coordinates.source_conformer_order, Some(expected));
    match coordinates.first_source_conformer().unwrap().unwrap() {
        CoordinateSourceConformer::ThreeD(first) => {
            assert_eq!(first.id(), 9);
            assert!(!first.is_3d());
            assert_eq!(first.coordinates()[0][2].to_bits(), (-0.0_f64).to_bits());
        }
        _ => panic!("actual first source conformer remains ThreeD even with false is3D"),
    }
}
