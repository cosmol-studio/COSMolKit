use cosmolkit::{CipLabelOptions, Molecule};
#[test]
fn modern_cip_computed_bool_property_matches_source_string_projection() {
    let original = Molecule::from_smiles("C[C@H](F)Cl").unwrap();
    assert!(!original.cip_computed());
    let labeled = original.with_cip_labels().unwrap();
    assert_eq!(labeled.property("_CIPComputed"), Some("1"));
    assert!(labeled.cip_computed());
    assert!(!original.cip_computed());
    let empty = original
        .with_cip_labels_with_options(&CipLabelOptions::default().with_atoms([]).with_bonds([]))
        .unwrap();
    assert!(empty.cip_computed());
}
