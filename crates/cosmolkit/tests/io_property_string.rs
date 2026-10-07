#![cfg(feature = "cap-io")]
//! Public source-property projection and missing-lookup regressions.
use cosmolkit::{
    AtomId, AtomSpec, BondId, BondOrder, BondSpec, Element, MoleculeBuilder, PropertyText,
    PropertyValue, property_value_to_text,
};
#[test]
fn atom_and_bond_source_strings_preserve_kinds_and_missing_lookup() {
    let cases = [
        (PropertyValue::String(" x\n".into()), " x\n"),
        (PropertyValue::Int(i32::MIN), "-2147483648"),
        (PropertyValue::UInt(u32::MAX), "4294967295"),
        (PropertyValue::IntVector(vec![-2, 0, 3]), "[-2,0,3]"),
        (PropertyValue::Double(1.25), "1.25"),
        (PropertyValue::Bool(true), "1"),
    ];
    for (value, expected) in cases {
        let mut builder = MoleculeBuilder::new();
        let a = builder.add_atom(
            AtomSpec::new(Element::C)
                .with_prop("source", value.clone())
                .unwrap(),
        );
        let b = builder.add_atom(AtomSpec::new(Element::C));
        let bond = builder
            .add_bond(
                BondSpec::new(a, b, BondOrder::Single)
                    .with_prop("source", value.clone())
                    .unwrap(),
            )
            .unwrap();
        let molecule = builder.build().unwrap();
        assert_eq!(
            molecule
                .atom(a)
                .and_then(|atom| atom.prop("source"))
                .map(property_value_to_text)
                .transpose()
                .unwrap()
                .as_ref()
                .map(PropertyText::as_bytes),
            Some(expected.as_bytes())
        );
        assert_eq!(
            molecule
                .bond(bond)
                .and_then(|bond| bond.prop("source"))
                .map(property_value_to_text)
                .transpose()
                .unwrap()
                .as_ref()
                .map(PropertyText::as_bytes),
            Some(expected.as_bytes())
        );
        assert_eq!(
            molecule
                .atom(a)
                .and_then(|atom| atom.prop("missing"))
                .map(property_value_to_text)
                .transpose()
                .unwrap(),
            None
        );
        assert_eq!(
            molecule
                .bond(bond)
                .and_then(|bond| bond.prop("missing"))
                .map(property_value_to_text)
                .transpose()
                .unwrap(),
            None
        );
        assert_eq!(
            molecule
                .atom(AtomId::new(19))
                .and_then(|atom| atom.prop("source"))
                .map(property_value_to_text)
                .transpose()
                .unwrap(),
            None
        );
        assert_eq!(
            molecule
                .bond(BondId::new(19))
                .and_then(|bond| bond.prop("source"))
                .map(property_value_to_text)
                .transpose()
                .unwrap(),
            None
        );
        assert_eq!(molecule.atom(a).unwrap().prop("source"), Some(&value));
        assert_eq!(molecule.bond(bond).unwrap().prop("source"), Some(&value));
    }
}

#[test]
fn canonical_property_text_projection_preserves_counted_string_and_vector_bytes() {
    let values = [
        (
            PropertyValue::String(PropertyText::from_bytes(b"a\0\xffb")),
            b"a\0\xffb".as_slice(),
        ),
        (
            PropertyValue::StringVector(vec![
                PropertyText::from_bytes(b"\0\xff"),
                PropertyText::from_bytes(b"A\x80"),
            ]),
            b"[\0\xff,A\x80]".as_slice(),
        ),
    ];
    for (value, bytes) in values {
        let before = value.clone();
        assert_eq!(property_value_to_text(&value).unwrap().as_bytes(), bytes);
        assert_eq!(value, before);
    }
    let _: fn(&PropertyValue) -> Result<PropertyText, cosmolkit::PropertyStringError> =
        property_value_to_text;
    let entry = cosmolkit::BINDING_CONTRACT
        .iter()
        .find(|entry| entry.semantic_id == "module.property_value_to_text")
        .unwrap();
    assert_eq!(entry.feature, "cap-io");
    assert_eq!(
        entry.callable.unwrap().state_model,
        cosmolkit::StateModel::ReadOnly
    );
    assert_eq!(entry.callable.unwrap().operation_semantic_id, None);
}
