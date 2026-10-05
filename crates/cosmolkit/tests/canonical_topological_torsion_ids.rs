#![cfg(feature = "cap-fingerprints")]
use cosmolkit::Molecule;

#[test]
fn original_native_row35_ids_preserve_order_and_duplicate_frequencies() {
    // Immutable original source-native TT corpus row35, pinned RDKit351f8f.
    let molecule = Molecule::from_smiles("B1BCCCC1").unwrap();
    let expected = vec![
        4429185057, 4437573633, 4437573633, 4437590017, 4437590017, 4437590049,
    ];
    assert_eq!(molecule.topological_torsion_ids().unwrap(), expected);
    assert_eq!(
        molecule.topological_torsion_ids_with_params(4).unwrap(),
        expected
    );
}

#[test]
fn original_native_no_torsion_and_unfolded_ids_remain_distinct_values() {
    assert!(
        Molecule::from_smiles("[He]")
            .unwrap()
            .topological_torsion_ids()
            .unwrap()
            .is_empty()
    );
    let molecule = Molecule::from_smiles("CCCCO").unwrap();
    assert_eq!(
        molecule.topological_torsion_ids().unwrap(),
        [4437590048, 12893306913]
    );
    let mut copied = molecule.topological_torsion_ids().unwrap();
    copied.clear();
    assert_eq!(
        molecule.topological_torsion_ids().unwrap(),
        [4437590048, 12893306913]
    );
}
