#![cfg(all(feature = "cap-smiles", feature = "cap-transforms"))]
use cosmolkit::{
    AtomId, AtomPropertyReadError, Molecule, OperationError, PropertyValue, PropertyValueKind,
};

#[test]
fn exact_types_missing_values_errors_and_cow() {
    let source = Molecule::from_smiles("CCO").unwrap();
    let atom = AtomId::new(0);
    macro_rules! roundtrip {
        ($write:ident, $set:ident, $read:ident, $value:expr, $expected:expr) => {{
            let mut tagged = source.$write(atom, "tag", $value).unwrap();
            assert_eq!(tagged.$read(atom, "tag").unwrap(), Some($expected));
            assert_eq!(source.$read(atom, "tag").unwrap(), None);
            let peer = tagged.clone();
            tagged.$set(atom, "other", $value).unwrap();
            assert_eq!(tagged.$read(atom, "other").unwrap(), Some($expected));
            assert_eq!(peer.$read(atom, "other").unwrap(), None);
        }};
    }
    roundtrip!(
        with_atom_property_int,
        set_atom_property_int_,
        atom_property_int,
        i32::MIN,
        i32::MIN
    );
    roundtrip!(
        with_atom_property_uint,
        set_atom_property_uint_,
        atom_property_uint,
        u32::MAX,
        u32::MAX
    );
    roundtrip!(
        with_atom_property_bool,
        set_atom_property_bool_,
        atom_property_bool,
        true,
        true
    );
    roundtrip!(
        with_atom_property_double,
        set_atom_property_double_,
        atom_property_double,
        -0.0,
        -0.0
    );
    roundtrip!(
        with_atom_property_string,
        set_atom_property_string_,
        atom_property_string,
        "α\0β",
        "α\0β".to_owned()
    );
    roundtrip!(
        with_atom_property_int_vector,
        set_atom_property_int_vector_,
        atom_property_int_vector,
        &[],
        Vec::<i32>::new()
    );
    roundtrip!(
        with_atom_property_string_vector,
        set_atom_property_string_vector_,
        atom_property_string_vector,
        &["α".into(), "".into()],
        vec!["α".to_owned(), "".to_owned()]
    );
    let number = source.with_atom_property_int(atom, "tag", 42).unwrap();
    match number.atom_property_string(atom, "tag").unwrap_err() {
        AtomPropertyReadError::Type(e) => {
            assert_eq!(e.expected(), PropertyValueKind::String);
            assert_eq!(e.actual(), PropertyValueKind::Int);
        }
        other => panic!("{other:?}"),
    }
    assert!(matches!(
        number.atom_property_uint(atom, "tag"),
        Err(AtomPropertyReadError::Type(_))
    ));
    assert!(matches!(
        source.atom_property_int(AtomId::new(9), "missing"),
        Err(AtomPropertyReadError::Operation(
            OperationError::AtomPropertyIndex { .. }
        ))
    ));
    let bytes = source
        .with_atom_property(atom, "tag", &PropertyValue::String(vec![0xff].into()))
        .unwrap();
    assert!(matches!(
        bytes.atom_property_string(atom, "tag"),
        Err(AtomPropertyReadError::Utf8(_))
    ));
    let vectors = source
        .with_atom_property(
            atom,
            "tag",
            &PropertyValue::StringVector(vec![vec![0xff].into()]),
        )
        .unwrap();
    assert!(matches!(
        vectors.atom_property_string_vector(atom, "tag"),
        Err(AtomPropertyReadError::Utf8(_))
    ));
    let minus_zero = source.with_atom_property_double(atom, "tag", -0.0).unwrap();
    assert_eq!(
        minus_zero
            .atom_property_double(atom, "tag")
            .unwrap()
            .unwrap()
            .to_bits(),
        (-0.0_f64).to_bits()
    );
    let mut unchanged = source.clone();
    assert!(matches!(
        unchanged.set_atom_property_int_(atom, "_CIPCode", 1),
        Err(OperationError::ReservedAtomPropertyKey { .. })
    ));
    assert_eq!(unchanged, source);
}
