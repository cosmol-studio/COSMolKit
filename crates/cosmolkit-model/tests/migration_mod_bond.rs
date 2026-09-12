use cosmolkit_model::{
    AtomId, Bond, BondDirection, BondId, BondOrder, BondSpec, BondStereo, BondValueError,
};

fn full_spec() -> BondSpec {
    BondSpec::new(AtomId::new(1), AtomId::new(3), BondOrder::Double)
        .with_order(BondOrder::Triple)
        .with_aromatic(true)
        .with_conjugated(true)
        .with_direction(BondDirection::BeginWedge)
        .with_stereo_atoms(AtomId::new(0), AtomId::new(4))
        .with_stereo(BondStereo::Cis)
        .with_unknown_stereo(true)
        .with_prop("ordinary", "kept")
        .unwrap()
        .with_computed_prop("computed", "cached")
        .unwrap()
}

#[test]
fn bond_id_is_a_stable_ordered_display_value() {
    let first = BondId::new(2);
    let second = BondId::new(9);
    assert_eq!(first.index(), 2);
    assert!(first < second);
    assert_eq!(first.to_string(), "2");
}

#[test]
fn bond_spec_covers_defaults_builders_validation_and_remapping() {
    let default = BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single);
    assert_eq!(default.begin(), AtomId::new(0));
    assert_eq!(default.end(), AtomId::new(1));
    assert_eq!(default.order(), BondOrder::Single);
    assert!(!default.is_aromatic());
    assert!(!default.is_conjugated());
    assert_eq!(default.direction(), BondDirection::None);
    assert_eq!(default.stereo(), BondStereo::None);
    assert_eq!(default.stereo_atoms(), None);
    assert!(!default.unknown_stereo());
    assert!(default.props().is_empty());
    assert!(default.computed_prop_names().is_empty());
    assert_eq!(default.validate(), Ok(()));

    let full = full_spec();
    assert_eq!(full.begin(), AtomId::new(1));
    assert_eq!(full.end(), AtomId::new(3));
    assert_eq!(full.order(), BondOrder::Triple);
    assert!(full.is_aromatic());
    assert!(full.is_conjugated());
    assert_eq!(full.direction(), BondDirection::BeginWedge);
    assert_eq!(full.stereo(), BondStereo::Cis);
    assert_eq!(full.stereo_atoms(), Some([AtomId::new(0), AtomId::new(4)]));
    assert!(full.unknown_stereo());
    assert_eq!(full.prop("ordinary"), Some("kept"));
    assert!(full.is_prop_computed("computed"));
    assert_eq!(full.validate(), Ok(()));

    let remapped = full.remapped_endpoints(
        AtomId::new(4),
        AtomId::new(5),
        Some([AtomId::new(2), AtomId::new(7)]),
    );
    assert_eq!(remapped.begin(), AtomId::new(4));
    assert_eq!(remapped.end(), AtomId::new(5));
    assert_eq!(
        remapped.stereo_atoms(),
        Some([AtomId::new(2), AtomId::new(7)])
    );
    assert_eq!(remapped.without_stereo_atoms().stereo_atoms(), None);
}

#[test]
fn bond_from_spec_preserves_facts_and_all_detached_mutators() {
    let mut bond = Bond::from_spec(BondId::new(6), full_spec());
    assert_eq!(bond.id(), BondId::new(6));
    assert_eq!(bond.begin(), AtomId::new(1));
    assert_eq!(bond.end(), AtomId::new(3));
    assert_eq!(bond.order(), BondOrder::Triple);
    assert!(bond.is_aromatic());
    assert!(bond.is_conjugated());
    assert_eq!(bond.direction(), BondDirection::BeginWedge);
    assert_eq!(bond.stereo(), BondStereo::Cis);
    assert_eq!(bond.stereo_atoms(), Some([AtomId::new(0), AtomId::new(4)]));
    assert!(bond.unknown_stereo());
    assert_eq!(bond.prop("ordinary"), Some("kept"));
    assert!(bond.is_prop_computed("computed"));
    assert_eq!(bond.validate(), Ok(()));

    bond = bond.remapped(BondId::new(2), AtomId::new(5), AtomId::new(8), None);
    assert_eq!(bond.id(), BondId::new(2));
    assert_eq!(bond.begin(), AtomId::new(5));
    assert_eq!(bond.end(), AtomId::new(8));
    assert_eq!(bond.stereo_atoms(), None);
    assert_eq!(bond.validate(), Err(BondValueError::StereoAtomsRequired));

    bond.set_id_for_construction(BondId::new(3));
    bond.set_endpoints(AtomId::new(0), AtomId::new(1));
    bond.set_order(BondOrder::Single);
    bond.set_aromatic(false);
    bond.set_conjugated(false);
    bond.set_direction(BondDirection::None);
    bond.set_unknown_stereo(false);
    assert_eq!(bond.set_stereo(BondStereo::None), Ok(()));
    bond.set_stereo_atoms(Some([AtomId::new(2), AtomId::new(3)]));
    assert_eq!(bond.set_stereo(BondStereo::Trans), Ok(()));

    assert_eq!(bond.id(), BondId::new(3));
    assert_eq!((bond.begin(), bond.end()), (AtomId::new(0), AtomId::new(1)));
    assert_eq!(bond.order(), BondOrder::Single);
    assert!(!bond.is_aromatic());
    assert!(!bond.is_conjugated());
    assert_eq!(bond.direction(), BondDirection::None);
    assert_eq!(bond.stereo(), BondStereo::Trans);
    assert!(!bond.unknown_stereo());
    assert_eq!(bond.validate(), Ok(()));
}

#[test]
fn stereo_validation_covers_every_value_and_failure_is_atomic() {
    let unrestricted = [
        BondStereo::None,
        BondStereo::Any,
        BondStereo::Z,
        BondStereo::E,
        BondStereo::AtropCw,
        BondStereo::AtropCcw,
    ];
    for stereo in unrestricted {
        let spec =
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Double).with_stereo(stereo);
        assert_eq!(spec.validate(), Ok(()));
        let mut bond = Bond::from_spec(BondId::new(0), spec);
        assert_eq!(bond.set_stereo(stereo), Ok(()));
    }

    for stereo in [BondStereo::Cis, BondStereo::Trans] {
        let spec =
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Double).with_stereo(stereo);
        assert_eq!(spec.validate(), Err(BondValueError::StereoAtomsRequired));

        let mut bond = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Double),
        );
        assert_eq!(
            bond.set_stereo(stereo),
            Err(BondValueError::StereoAtomsRequired)
        );
        assert_eq!(bond.stereo(), BondStereo::None);
        bond.set_stereo_atoms(Some([AtomId::new(2), AtomId::new(3)]));
        assert_eq!(bond.set_stereo(stereo), Ok(()));
    }
}

#[test]
fn checked_bond_properties_cover_empty_overwrite_membership_and_clear() {
    let base = || BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single);
    assert_eq!(
        base().with_prop("", "value"),
        Err(BondValueError::EmptyPropertyKey)
    );
    assert_eq!(
        base().with_computed_prop("", "value"),
        Err(BondValueError::EmptyPropertyKey)
    );

    let spec = base()
        .with_computed_prop("cache", "first")
        .unwrap()
        .with_computed_prop("cache", "second")
        .unwrap()
        .with_prop("ordinary", "first")
        .unwrap()
        .with_prop("ordinary", "second")
        .unwrap();
    assert_eq!(spec.prop("cache"), Some("second"));
    assert_eq!(spec.prop("ordinary"), Some("second"));
    assert_eq!(spec.computed_prop_names().len(), 1);

    let mut bond = Bond::from_spec(BondId::new(0), spec);
    assert_eq!(
        bond.set_prop("", "value"),
        Err(BondValueError::EmptyPropertyKey)
    );
    assert_eq!(
        bond.set_computed_prop("", "value"),
        Err(BondValueError::EmptyPropertyKey)
    );
    bond.set_prop("cache", "ordinary overwrite").unwrap();
    assert!(bond.is_prop_computed("cache"));
    bond.set_computed_prop("cache", "computed overwrite")
        .unwrap();
    assert_eq!(bond.computed_prop_names().len(), 1);

    bond.clear_prop("missing");
    bond.clear_prop("cache");
    assert_eq!(bond.prop("cache"), None);
    assert!(!bond.is_prop_computed("cache"));
    bond.set_computed_prop("temporary", "gone").unwrap();
    bond.clear_computed_props();
    assert_eq!(bond.prop("temporary"), None);
    assert_eq!(bond.prop("ordinary"), Some("second"));
    assert!(bond.computed_prop_names().is_empty());
}
