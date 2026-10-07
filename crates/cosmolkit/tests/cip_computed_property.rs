use cosmolkit::{CipLabelOptions, Molecule};
#[test]
fn modern_cip_computed_bool_property_matches_source_string_projection() {
    let original = Molecule::from_smiles("C[C@H](F)Cl").unwrap();
    assert!(!original.cip_computed().unwrap());
    let labeled = original.with_cip_labels().unwrap();
    assert_eq!(
        labeled.property("_CIPComputed"),
        Some(&cosmolkit::PropertyValue::Bool(true))
    );
    assert!(labeled.cip_computed().unwrap());
    assert!(!original.cip_computed().unwrap());
    let empty = original
        .with_cip_labels_with_options(&CipLabelOptions::default().with_atoms([]).with_bonds([]))
        .unwrap();
    assert!(empty.cip_computed().unwrap());
}

#[test]
fn canonical_cip_query_reports_wrong_source_tags_and_assignment_is_atomic() {
    use cosmolkit::{
        CoordinateBlock, MoleculeProperties, PropertyValue, PropertyValueKind, TopologyBlock,
    };
    let construct = |properties| {
        Molecule::from_parts(
            TopologyBlock::default(),
            CoordinateBlock::default(),
            properties,
        )
        .unwrap()
    };
    let wrong_bool = construct(
        MoleculeProperties::default()
            .with_computed_prop("_CIPComputed", "1")
            .unwrap(),
    );
    let error = wrong_bool.cip_computed().unwrap_err();
    assert_eq!(error.expected(), PropertyValueKind::Bool);
    assert_eq!(error.actual(), PropertyValueKind::String);
    let wrong_computed = construct(
        MoleculeProperties::default()
            .with_prop("__computedProps", PropertyValue::Int(7))
            .unwrap(),
    );
    let observer = wrong_computed.clone();
    let error = wrong_computed.cip_computed().unwrap_err();
    assert_eq!(error.expected(), PropertyValueKind::StringVector);
    assert_eq!(error.actual(), PropertyValueKind::Int);
    assert!(matches!(
        wrong_computed.with_cip_labels(),
        Err(cosmolkit::OperationError::CipLabeler(_))
    ));
    assert_eq!(wrong_computed, observer);
    let missing = construct(MoleculeProperties::default());
    assert!(!missing.cip_computed().unwrap());
    let ordinary = construct(
        MoleculeProperties::default()
            .with_prop("_CIPComputed", true)
            .unwrap(),
    );
    assert!(!ordinary.cip_computed().unwrap());
}
