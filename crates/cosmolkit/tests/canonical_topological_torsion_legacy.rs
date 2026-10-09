#![cfg(feature = "cap-fingerprints")]
//! Legacy source-native fixed conditions through the canonical public boundary.
use cosmolkit::{LegacyTopologicalTorsionParams, Molecule};

#[test]
fn legacy_default_unfolded_and_hashed_source_values_are_distinct() {
    let molecule = Molecule::from_smiles("CCCCO").unwrap();
    let unfolded = molecule
        .fingerprint_topological_torsion_sparse_count_legacy()
        .unwrap();
    assert_eq!(unfolded.length(), (1_u64 << 36) - 1);
    assert_eq!(
        unfolded
            .nonzero_elements()
            .iter()
            .map(|(&k, &v)| (k, v))
            .collect::<Vec<_>>(),
        [(4437590048, 1), (12893306913, 1)]
    );
    let params = LegacyTopologicalTorsionParams {
        fp_size: 1000,
        ..Default::default()
    };
    let hashed = molecule
        .fingerprint_topological_torsion_count_legacy_with_params(&params)
        .unwrap();
    assert_eq!(hashed.length(), 1000);
    assert_eq!(
        hashed
            .nonzero_elements()
            .iter()
            .map(|(&k, &v)| (k, v))
            .collect::<Vec<_>>(),
        [(24, 1), (288, 1)]
    );
    let defaults = LegacyTopologicalTorsionParams::default();
    assert_eq!(
        molecule
            .fingerprint_topological_torsion_count_legacy()
            .unwrap(),
        molecule
            .fingerprint_topological_torsion_count_legacy_with_params(&defaults)
            .unwrap()
    );
    assert_eq!(
        molecule.fingerprint_topological_torsion_legacy().unwrap(),
        molecule
            .fingerprint_topological_torsion_legacy_with_params(&defaults)
            .unwrap()
    );
    assert_eq!(
        unfolded,
        molecule
            .fingerprint_topological_torsion_sparse_count_legacy_with_params(&defaults)
            .unwrap()
    );
}

#[test]
fn legacy_optional_roots_are_distinct_and_failure_keeps_input_usable() {
    let molecule = Molecule::from_smiles("CCCCO").unwrap();
    let empty = LegacyTopologicalTorsionParams {
        from_atoms: Some(vec![]),
        ..Default::default()
    };
    assert!(
        molecule
            .fingerprint_topological_torsion_sparse_count_legacy_with_params(&empty)
            .unwrap()
            .nonzero_elements()
            .is_empty()
    );
    assert!(
        molecule
            .fingerprint_topological_torsion_count_legacy_with_params(&empty)
            .unwrap()
            .nonzero_elements()
            .is_empty()
    );
    assert!(
        molecule
            .fingerprint_topological_torsion_legacy_with_params(&empty)
            .unwrap()
            .on_bits()
            .is_empty()
    );
    let invalid = LegacyTopologicalTorsionParams {
        custom_atom_invariants: Some(vec![1]),
        ..Default::default()
    };
    for error in [
        molecule
            .fingerprint_topological_torsion_sparse_count_legacy_with_params(&invalid)
            .unwrap_err(),
        molecule
            .fingerprint_topological_torsion_count_legacy_with_params(&invalid)
            .unwrap_err(),
        molecule
            .fingerprint_topological_torsion_legacy_with_params(&invalid)
            .unwrap_err(),
    ] {
        assert!(error.to_string().contains("bad atomInvariants size"));
        assert!(std::error::Error::source(&error).is_some());
    }
    assert_eq!(
        molecule
            .fingerprint_topological_torsion_sparse_count_legacy()
            .unwrap()
            .nonzero_elements()
            .values()
            .sum::<i32>(),
        2
    );
}

#[test]
fn source_four_and_nonfour_entry_thresholds_survive_public_transport() {
    let molecule = Molecule::from_smiles("CCCCCCCCCCCC").unwrap();
    for entry in [1, 2, 4, 6] {
        let params = LegacyTopologicalTorsionParams::new(
            4,
            false,
            16 * entry,
            entry,
            None,
            None,
            Some(vec![7; 12]),
        );
        let block_params = LegacyTopologicalTorsionParams {
            fp_size: 16,
            ..params.clone()
        };
        let counts = molecule
            .fingerprint_topological_torsion_count_legacy_with_params(&block_params)
            .unwrap();
        let (&block, &count) = counts.nonzero_elements().first_key_value().unwrap();
        assert_eq!(counts.nonzero_elements().len(), 1);
        assert_eq!(count, 9);
        let bit_vector = molecule
            .fingerprint_topological_torsion_legacy_with_params(&params)
            .unwrap();
        assert_eq!(
            bit_vector.on_bits(),
            (0..entry)
                .map(|offset| (block as u32) * entry + offset)
                .collect::<Vec<_>>()
        );
        assert_eq!(bit_vector.n_bits(), 16 * entry);
    }
}
