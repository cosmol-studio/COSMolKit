use cosmolkit_core::get_molecule_fragments;
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, Conformer2D, Conformer3D,
    CoordinateBlock, CoordinateDimension, Element, MoleculeProperties, TopologyBlock,
};
fn check(three_d: bool) {
    // Pinned MolOps.cpp704-831 slow route clones ROMol conformers then only
    // deletes atom rows; ROMol.cpp130-134 and RWMol.cpp901-913 have no finite guard.
    let topology = TopologyBlock::try_from_parts(
        (0..4)
            .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
            .collect(),
        vec![
            Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            ),
            Bond::from_spec(
                BondId::new(1),
                BondSpec::new(AtomId::new(2), AtomId::new(3), BondOrder::Single),
            ),
        ],
        vec![],
        vec![],
    )
    .unwrap();
    let nan = f64::from_bits(0xfff8_0000_0000_007b);
    let coordinates = if three_d {
        CoordinateBlock {
            conformers_3d: vec![Conformer3D::new(
                9,
                vec![
                    [nan, 1.0, -0.0],
                    [2.0, 3.0, 4.0],
                    [5.0, 6.0, 7.0],
                    [8.0, 9.0, 10.0],
                ],
                false,
            )],
            source_conformer_order: Some(vec![CoordinateDimension::ThreeD]),
            ..Default::default()
        }
    } else {
        CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(
                9,
                vec![[nan, -0.0], [2.0, 3.0], [4.0, 5.0], [6.0, 7.0]],
            )],
            source_conformer_order: Some(vec![CoordinateDimension::TwoD]),
            ..Default::default()
        }
    };
    coordinates.validate_for_atom_count(4).unwrap();
    let fragments = get_molecule_fragments(
        &topology,
        &coordinates,
        &MoleculeProperties::default(),
        false,
        true,
    )
    .unwrap();
    assert_eq!(fragments.len(), 2);
    let first = fragments[0].coordinates();
    assert_eq!(
        first.source_conformer_order,
        coordinates.source_conformer_order
    );
    let value = if three_d {
        first.conformers_3d[0].coordinates()[0][0]
    } else {
        first.conformers_2d[0].coordinates()[0][0]
    };
    assert_eq!(value.to_bits(), nan.to_bits());
    assert_eq!(fragments[0].topology().atoms.len(), 2);
    assert_eq!(fragments[1].topology().atoms.len(), 2);
}
#[test]
fn fragment_copy_transports_source_nan_2d_bits() {
    check(false);
}
#[test]
fn fragment_copy_transports_source_nan_3d_bits() {
    check(true);
}
