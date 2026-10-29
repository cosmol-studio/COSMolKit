use cosmolkit_wasm::{Element, element_info};

#[test]
fn element_identity_matches_public_facade_without_consuming_copy_values() {
    for number in 0..=u8::MAX {
        let actual = Element::from_atomic_number(number);
        let expected = cosmolkit::Element::from_atomic_number(number);
        assert_eq!(
            actual.is_some(),
            expected.is_some(),
            "atomic number {number}"
        );
        if let (Some(actual), Some(expected)) = (actual, expected) {
            assert_eq!(actual.atomic_number(), expected.atomic_number());
            assert_eq!(actual.symbol(), expected.symbol());
            assert_eq!(actual.atomic_number(), number);
        }
    }
    for symbol in ["*", "C", "Og", "Uut", "Uup", "c", "", "C\0", " C", "é"] {
        let actual = Element::from_symbol(symbol);
        let expected = cosmolkit::Element::from_symbol(symbol);
        assert_eq!(
            actual.map(|element| element.atomic_number()),
            expected.map(|element| element.atomic_number())
        );
    }
}

#[test]
fn element_info_preserves_every_field_and_owns_returned_array_storage() {
    for number in 0..=118 {
        let element = Element::from_atomic_number(number).unwrap();
        let actual = element_info(&element);
        let expected =
            cosmolkit::element_info(cosmolkit::Element::from_atomic_number(number).unwrap());
        assert_eq!(
            actual.element().atomic_number(),
            expected.element.atomic_number()
        );
        assert_eq!(actual.symbol(), expected.symbol);
        assert_eq!(actual.atomic_number(), expected.atomic_number);
        assert_eq!(actual.period(), expected.period);
        assert_eq!(actual.outer_electrons(), expected.outer_electrons);
        assert_eq!(actual.valences(), expected.valences);
        assert_eq!(actual.rb0().to_bits(), expected.rb0.to_bits());
        assert_eq!(
            actual.atomic_weight().to_bits(),
            expected.atomic_weight.to_bits()
        );
        let mut detached = actual.valences();
        detached.push(i32::MAX);
        assert_eq!(actual.valences(), expected.valences);
        assert_eq!(element.atomic_number(), number);
    }
}
