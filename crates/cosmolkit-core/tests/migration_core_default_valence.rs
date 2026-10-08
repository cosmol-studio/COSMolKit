use cosmolkit_core::{ValenceError, rdkit_default_valence};

#[test]
fn default_valence_is_first_signed_entry_of_all_119_pinned_rows() {
    // RDKit 2026.03.1 atomic_data.cpp periodicTableAtomData, first numerical
    // row for each atomic number. Literal source expectations include dummy,
    // noble gas, transition metal, and multivalence rows (not inferred values).
    let expected: [i32; 119] = [
        -1, 1, 0, 1, 2, 3, 4, 3, 2, 1, 0, 1, 2, 3, 4, 3, 2, 1, 0, 1, 2, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, 3, 4, 3, 2, 1, 0, 1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 3, 2, 3, 2, 1,
        0, 1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, 2, 3, 2, 1, 0, 1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    ];
    for (atomic_number, expected) in expected.into_iter().enumerate() {
        assert_eq!(
            rdkit_default_valence(atomic_number as u8),
            Ok(expected),
            "Z={atomic_number}"
        );
    }
}

#[test]
fn default_valence_outside_pinned_rows_preserves_checked_lookup_error() {
    for atomic_number in 119..=u8::MAX {
        assert_eq!(
            rdkit_default_valence(atomic_number),
            Err(ValenceError::PeriodicTableLookup {
                atomic_number,
                field: "valences",
            })
        );
    }
}
