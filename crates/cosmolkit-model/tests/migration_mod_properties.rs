use std::str::FromStr;

use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, CipDescriptor, CipDescriptorError,
    Element, MoleculeProperties, MoleculePropertyError, SdfPropertyList, SdfPropertyListTarget,
};

const CIP_CASES: [(CipDescriptor, &str); 10] = [
    (CipDescriptor::R, "R"),
    (CipDescriptor::S, "S"),
    (CipDescriptor::LowerR, "r"),
    (CipDescriptor::LowerS, "s"),
    (CipDescriptor::E, "E"),
    (CipDescriptor::Z, "Z"),
    (CipDescriptor::M, "M"),
    (CipDescriptor::P, "P"),
    (CipDescriptor::LowerM, "m"),
    (CipDescriptor::LowerP, "p"),
];

#[test]
fn metadata_defaults_names_and_ordered_sdf_records_are_exact() {
    let default = MoleculeProperties::default();
    assert_eq!(default.name(), None);
    assert!(default.sdf_data_fields().is_empty());
    assert!(default.sdf_property_lists().is_empty());
    assert!(default.props().is_empty());
    assert!(default.computed_prop_names().is_empty());

    let atom_list = SdfPropertyList::new(
        SdfPropertyListTarget::Atom,
        "shared",
        vec![Some(" a ".to_owned()), None],
    );
    let bond_list = SdfPropertyList::new(
        SdfPropertyListTarget::Bond,
        "shared",
        vec![Some(String::new())],
    );
    let properties = default
        .with_name("first")
        .with_sdf_data_field("duplicate", "one")
        .with_sdf_data_field("duplicate", "two")
        .with_sdf_property_list(atom_list.clone())
        .with_sdf_property_list(bond_list.clone())
        .with_name("");

    assert_eq!(properties.name(), Some(""));
    assert_eq!(
        properties.sdf_data_fields(),
        &[
            ("duplicate".to_owned(), "one".to_owned()),
            ("duplicate".to_owned(), "two".to_owned()),
        ]
    );
    assert_eq!(properties.sdf_property_lists(), &[atom_list, bond_list]);
    assert_eq!(
        properties.sdf_property_lists()[0].target(),
        SdfPropertyListTarget::Atom
    );
    assert_eq!(properties.sdf_property_lists()[0].name(), "shared");
    assert_eq!(
        properties.sdf_property_lists()[0].values(),
        &[Some(" a ".to_owned()), None]
    );
}

#[test]
fn ordinary_and_computed_properties_cover_errors_overwrites_and_clears() {
    assert_eq!(
        MoleculeProperties::default().with_prop("", "value"),
        Err(MoleculePropertyError::EmptyKey)
    );
    assert_eq!(
        MoleculeProperties::default().with_computed_prop("", "value"),
        Err(MoleculePropertyError::EmptyKey)
    );

    let mut properties = MoleculeProperties::default()
        .with_name("unchanged")
        .with_sdf_data_field("raw", "record")
        .with_prop("ordinary", "first")
        .unwrap()
        .with_prop("ordinary", "second")
        .unwrap()
        .with_computed_prop("computed", "first")
        .unwrap()
        .with_computed_prop("computed", "second")
        .unwrap();
    let before_error = properties.clone();
    assert_eq!(
        properties.set_prop("", "not stored"),
        Err(MoleculePropertyError::EmptyKey)
    );
    assert_eq!(properties, before_error);
    assert_eq!(
        properties.set_computed_prop("", "not stored"),
        Err(MoleculePropertyError::EmptyKey)
    );
    assert_eq!(properties, before_error);

    assert_eq!(properties.prop("ordinary"), Some("second"));
    assert_eq!(properties.prop("missing"), None);
    assert_eq!(properties.prop("computed"), Some("second"));
    assert!(properties.is_prop_computed("computed"));
    assert_eq!(properties.computed_prop_names().len(), 1);

    properties
        .set_prop("computed", "ordinary overwrite")
        .unwrap();
    assert_eq!(properties.prop("computed"), Some("ordinary overwrite"));
    assert!(properties.is_prop_computed("computed"));
    properties.clear_prop("absent");
    properties.clear_prop("computed");
    assert_eq!(properties.prop("computed"), None);
    assert!(!properties.is_prop_computed("computed"));

    properties.set_computed_prop("temporary", "gone").unwrap();
    properties.set_prop("cache_like_name", "kept").unwrap();
    properties.clear_computed_props();
    properties.clear_computed_props();
    assert_eq!(properties.prop("temporary"), None);
    assert_eq!(properties.prop("ordinary"), Some("second"));
    assert_eq!(properties.prop("cache_like_name"), Some("kept"));
    assert!(properties.computed_prop_names().is_empty());
    assert_eq!(properties.name(), Some("unchanged"));
    assert_eq!(
        properties.sdf_data_fields(),
        &[("raw".to_owned(), "record".to_owned())]
    );
}

#[test]
fn sdf_property_lists_remap_each_target_in_order_and_preserve_metadata() {
    let mut properties = MoleculeProperties::default()
        .with_name("molecule")
        .with_sdf_data_field("raw", "field")
        .with_prop("ordinary", "value")
        .unwrap()
        .with_computed_prop("computed", "cached")
        .unwrap()
        .with_sdf_property_list(SdfPropertyList::new(
            SdfPropertyListTarget::Bond,
            "bond-list",
            vec![None, Some("b1".to_owned())],
        ))
        .with_sdf_property_list(SdfPropertyList::new(
            SdfPropertyListTarget::Atom,
            "atom-list",
            vec![Some("a0".to_owned()), None, Some("a2".to_owned())],
        ))
        .with_sdf_property_list(SdfPropertyList::new(
            SdfPropertyListTarget::Atom,
            "atom-subset",
            vec![Some("x0".to_owned()), Some("x1".to_owned())],
        ));

    properties.remap_topology(
        &[
            Some(AtomId::new(2)),
            Some(AtomId::new(1)),
            None,
            Some(AtomId::new(99)),
            Some(AtomId::new(0)),
        ],
        &[
            Some(BondId::new(1)),
            Some(BondId::new(0)),
            None,
            Some(BondId::new(8)),
        ],
    );

    let lists = properties.sdf_property_lists();
    assert_eq!(lists.len(), 3);
    assert_eq!(lists[0].target(), SdfPropertyListTarget::Bond);
    assert_eq!(lists[0].name(), "bond-list");
    assert_eq!(
        lists[0].values(),
        &[Some("b1".to_owned()), None, None, None]
    );
    assert_eq!(lists[1].target(), SdfPropertyListTarget::Atom);
    assert_eq!(lists[1].name(), "atom-list");
    assert_eq!(
        lists[1].values(),
        &[
            Some("a2".to_owned()),
            None,
            None,
            None,
            Some("a0".to_owned()),
        ]
    );
    assert_eq!(lists[2].target(), SdfPropertyListTarget::Atom);
    assert_eq!(lists[2].name(), "atom-subset");
    assert_eq!(
        lists[2].values(),
        &[
            None,
            Some("x1".to_owned()),
            None,
            None,
            Some("x0".to_owned())
        ]
    );

    assert_eq!(properties.name(), Some("molecule"));
    assert_eq!(
        properties.sdf_data_fields(),
        &[("raw".to_owned(), "field".to_owned())]
    );
    assert_eq!(properties.prop("ordinary"), Some("value"));
    assert_eq!(properties.prop("computed"), Some("cached"));
    assert!(properties.is_prop_computed("computed"));
}

#[test]
fn cip_descriptor_ten_spelling_matrix_is_exact_and_case_sensitive() {
    for (descriptor, spelling) in CIP_CASES {
        assert_eq!(descriptor.as_str(), spelling);
        assert_eq!(descriptor.to_string(), spelling);
        assert_eq!(CipDescriptor::from_str(spelling), Ok(descriptor));
    }
    assert_ne!(CipDescriptor::R, CipDescriptor::LowerR);
    assert_ne!(CipDescriptor::S, CipDescriptor::LowerS);
    assert_ne!(CipDescriptor::M, CipDescriptor::LowerM);
    assert_ne!(CipDescriptor::P, CipDescriptor::LowerP);
}

#[test]
fn cip_descriptor_invalid_spellings_preserve_the_exact_input() {
    for value in [
        "", " R", "R ", "rR", "NONE", "UNKNOWN", "e", "z", "SP_4", "TBPY_5", "OC_6", "ns",
        "seqTrans", "seqCis",
    ] {
        assert_eq!(
            CipDescriptor::from_str(value),
            Err(CipDescriptorError::InvalidStoredDescriptor {
                value: value.to_owned(),
            })
        );
    }
}

#[test]
fn atom_cip_projection_covers_absent_all_supported_and_invalid_values() {
    let absent = Atom::from_spec(AtomId::new(7), AtomSpec::new(Element::C));
    assert_eq!(absent.cip_descriptor(), Ok(None));

    for (descriptor, spelling) in CIP_CASES {
        let atom = Atom::from_spec(
            AtomId::new(7),
            AtomSpec::new(Element::C)
                .with_prop("_CIPCode", spelling)
                .unwrap(),
        );
        assert_eq!(atom.cip_descriptor(), Ok(Some(descriptor)));
        assert_eq!(atom.prop("_CIPCode"), Some(spelling));
    }

    let invalid = Atom::from_spec(
        AtomId::new(7),
        AtomSpec::new(Element::C)
            .with_prop("_CIPCode", "UNKNOWN")
            .unwrap(),
    );
    assert_eq!(
        invalid.cip_descriptor(),
        Err(CipDescriptorError::InvalidStoredDescriptor {
            value: "UNKNOWN".to_owned(),
        })
    );
    assert_eq!(invalid.prop("_CIPCode"), Some("UNKNOWN"));
}

#[test]
fn bond_cip_projection_covers_absent_all_supported_and_invalid_values() {
    let bond_spec = || BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Double);
    let absent = Bond::from_spec(BondId::new(4), bond_spec());
    assert_eq!(absent.cip_descriptor(), Ok(None));

    for (descriptor, spelling) in CIP_CASES {
        let bond = Bond::from_spec(
            BondId::new(4),
            bond_spec().with_prop("_CIPCode", spelling).unwrap(),
        );
        assert_eq!(bond.cip_descriptor(), Ok(Some(descriptor)));
        assert_eq!(bond.prop("_CIPCode"), Some(spelling));
    }

    let invalid = Bond::from_spec(
        BondId::new(4),
        bond_spec().with_prop("_CIPCode", "seqCis").unwrap(),
    );
    assert_eq!(
        invalid.cip_descriptor(),
        Err(CipDescriptorError::InvalidStoredDescriptor {
            value: "seqCis".to_owned(),
        })
    );
    assert_eq!(invalid.prop("_CIPCode"), Some("seqCis"));
}
