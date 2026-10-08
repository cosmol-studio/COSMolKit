use cosmolkit_core::{ValenceError, rdkit_valence_list, required_valence_list};

#[test]
fn complete_signed_valence_lists_preserve_pinned_rows_order_and_shared_borrow() {
    // Literal complete valence suffixes of the119first numeric rows from pinned
    // RDKit2026.03.1 atomic_data.cpp, including every terminal -1 sentinel.
    let expected: [&[i32]; 119] = [
        &[-1],
        &[1],
        &[0],
        &[1, -1],
        &[2],
        &[3],
        &[4],
        &[3],
        &[2],
        &[1],
        &[0],
        &[1, -1],
        &[2, -1],
        &[3],
        &[4],
        &[3, 5],
        &[2, 4, 6],
        &[1],
        &[0],
        &[1, -1],
        &[2, -1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[3],
        &[4],
        &[3, 5],
        &[2, 4, 6],
        &[1],
        &[0],
        &[1, -1],
        &[2, -1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[3],
        &[2, 4],
        &[3, 5],
        &[2, 4, 6],
        &[1, 3, 5],
        &[0, 2, 4, 6],
        &[1],
        &[2, -1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[2, 4],
        &[3, 5],
        &[2, 4, 6],
        &[1, 3, 5],
        &[0],
        &[1],
        &[2, -1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
        &[-1],
    ];
    for (atomic_number, expected) in expected.into_iter().enumerate() {
        let atomic_number = atomic_number as u8;
        let actual = required_valence_list(atomic_number).unwrap();
        assert_eq!(actual, expected, "Z={atomic_number}");
        assert_eq!(rdkit_valence_list(atomic_number), Ok(Some(expected)));
        assert!(std::ptr::eq(
            actual,
            required_valence_list(atomic_number).unwrap()
        ));
        assert!(std::ptr::eq(
            actual,
            rdkit_valence_list(atomic_number).unwrap().unwrap()
        ));
    }
}

#[test]
fn valence_list_out_of_range_returns_lookup_error_from_both_existing_projections() {
    for atomic_number in 119..=u8::MAX {
        let error = ValenceError::PeriodicTableLookup {
            atomic_number,
            field: "valences",
        };
        assert_eq!(required_valence_list(atomic_number), Err(error.clone()));
        assert_eq!(rdkit_valence_list(atomic_number), Err(error));
    }
}
