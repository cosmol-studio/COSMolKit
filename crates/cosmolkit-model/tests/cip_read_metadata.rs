use cosmolkit_model::{Atom, AtomId, AtomSpec, PropertyValue};
use cosmolkit_types::Element;
#[test]
fn cip_metadata_preserves_order_unsigned_rank_absence_and_invalid_values() {
    let make = |v| {
        Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_prop("_CIPNeighborOrder", v)
                .unwrap(),
        )
    };
    assert_eq!(
        make(PropertyValue::String("[3,2,0]".into()))
            .cip_neighbor_order()
            .unwrap(),
        Some(vec![3, 2, 0])
    );
    assert_eq!(
        make(PropertyValue::IntVector(vec![3, 2, 0]))
            .cip_neighbor_order()
            .unwrap(),
        Some(vec![3, 2, 0])
    );
    assert!(
        make(PropertyValue::String("[3,no]".into()))
            .cip_neighbor_order()
            .is_err()
    );
    assert!(
        make(PropertyValue::IntVector(vec![-1]))
            .cip_neighbor_order()
            .is_err()
    );
    assert!(
        make(PropertyValue::Bool(true))
            .cip_neighbor_order()
            .is_err()
    );
    let atom = Atom::from_spec(
        AtomId::new(0),
        AtomSpec::new(Element::C)
            .with_prop("_CIPRank", PropertyValue::UInt(u32::MAX))
            .unwrap(),
    );
    assert_eq!(atom.cip_rank().unwrap(), Some(u32::MAX));
    assert_eq!(atom.cip_neighbor_order().unwrap(), None);
    let wrong = Atom::from_spec(
        AtomId::new(0),
        AtomSpec::new(Element::C)
            .with_prop("_CIPRank", PropertyValue::Int(1))
            .unwrap(),
    );
    assert!(wrong.cip_rank().is_err());
}
