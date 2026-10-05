use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, CipDescriptorError, PropertyValue,
};
use cosmolkit_types::{BondOrder, Element};
#[test]
fn present_nonstring_cip_atom_and_bond_properties_are_typed_errors() {
    let atom = Atom::from_spec(
        AtomId::new(0),
        AtomSpec::new(Element::C)
            .with_prop("_CIPCode", PropertyValue::UInt(7))
            .unwrap(),
    );
    assert!(matches!(
        atom.cip_descriptor(),
        Err(CipDescriptorError::Property(_))
    ));
    let bond = Bond::from_spec(
        BondId::new(0),
        BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
            .with_prop("_CIPCode", PropertyValue::Int(1))
            .unwrap(),
    );
    assert!(matches!(
        bond.cip_descriptor(),
        Err(CipDescriptorError::Property(_))
    ));
}
