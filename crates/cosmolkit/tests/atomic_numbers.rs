#![cfg(feature = "cap-smiles")]

use cosmolkit::Molecule;

#[test]
fn atomic_numbers_are_ordered_independent_snapshots() {
    assert!(Molecule::new().atomic_numbers().is_empty());
    let molecule = Molecule::from_smiles("ClC[13CH2]N.*.[Na+]").unwrap();
    let expected = vec![17, 6, 6, 7, 0, 11];
    let mut numbers = molecule.atomic_numbers();
    assert_eq!(numbers, expected);
    assert_eq!(
        numbers,
        molecule
            .atoms()
            .iter()
            .map(|a| a.atomic_number())
            .collect::<Vec<_>>()
    );
    numbers[0] = 0;
    assert_eq!(molecule.atomic_numbers(), expected);
    let retained = molecule.atomic_numbers();
    drop(molecule);
    assert_eq!(retained, expected);
}
