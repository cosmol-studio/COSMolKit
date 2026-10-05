use cosmolkit_core::atom_metadata;
use cosmolkit_model::{Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, TopologyBlock};
use cosmolkit_types::{BondOrder, Element};
#[test]
fn metadata_retains_degree_explicit_and_implicit_hydrogen_distinction() {
    let atoms = vec![
        Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
        Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::O)),
        Atom::from_spec(AtomId::new(2), AtomSpec::new(Element::H)),
    ];
    let bonds = vec![
        Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        ),
        Bond::from_spec(
            BondId::new(1),
            BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Single),
        ),
    ];
    let topology = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap();
    let before = topology.clone();
    let rows = atom_metadata(&topology).unwrap();
    assert_eq!(rows.len(), 3);
    assert_eq!(
        (
            rows[0].degree,
            rows[0].explicit_valence,
            rows[0].implicit_hydrogens,
            rows[0].total_hydrogens,
            rows[0].total_valence
        ),
        (1, 1, 3, 3, 4)
    );
    assert_eq!(
        (
            rows[1].degree,
            rows[1].explicit_valence,
            rows[1].implicit_hydrogens,
            rows[1].total_hydrogens,
            rows[1].total_valence
        ),
        (2, 2, 0, 0, 2)
    );
    assert_eq!(topology, before);
}
#[test]
fn metadata_propagates_invalid_valence_and_empty_is_present_empty() {
    assert!(atom_metadata(&TopologyBlock::default()).unwrap().is_empty());
    let atom = Atom::from_spec(
        AtomId::new(0),
        AtomSpec::new(Element::C).with_explicit_hydrogens(5),
    );
    let topology = TopologyBlock::try_from_parts(vec![atom], vec![], vec![], vec![]).unwrap();
    assert!(atom_metadata(&topology).is_err());
}
