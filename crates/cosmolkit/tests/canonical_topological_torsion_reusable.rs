#![cfg(feature = "cap-fingerprints")]
//! Public-chain proposal using pinned native source values and state contracts.
use cosmolkit::{
    FingerprintAdditionalOutput, Molecule, TopologicalTorsionCallParams,
    TopologicalTorsionFingerprintGenerator,
};

#[test]
fn source_four_bulk_forms_preserve_none_order_and_independent_workers() {
    let generator = TopologicalTorsionFingerprintGenerator::new(None, None).unwrap();
    let molecules = ["CCCCO", "CCCC", "C1CC1"].map(|s| Molecule::from_smiles(s).unwrap());
    let rows = [
        Some(&molecules[0]),
        None,
        Some(&molecules[1]),
        Some(&molecules[2]),
        None,
    ];
    for workers in [1, 2, 7] {
        let bits = generator.fingerprints(&rows, workers).unwrap();
        assert!(bits[1].is_none() && bits[4].is_none());
        assert_eq!(bits[0].as_ref().unwrap().on_bits(), [0, 1]);
        assert_eq!(bits[2].as_ref().unwrap().on_bits(), [60]);
        assert_eq!(bits[3].as_ref().unwrap().on_bits(), [284]);
        let sparse = generator.sparse_fingerprints(&rows, workers).unwrap();
        let counts = generator.counts(&rows, workers).unwrap();
        let unfolded = generator.sparse_counts(&rows, workers).unwrap();
        for (index, molecule) in rows.iter().enumerate() {
            assert_eq!(
                sparse[index],
                molecule.map(|m| m
                    .topological_torsion_sparse_fingerprint_with_generator(&generator, None, None)
                    .unwrap())
            );
            assert_eq!(
                counts[index],
                molecule.map(|m| m
                    .topological_torsion_count_fingerprint_with_generator(&generator, None, None)
                    .unwrap())
            );
            assert_eq!(
                unfolded[index],
                molecule.map(|m| m
                    .topological_torsion_sparse_count_fingerprint_with_generator(
                        &generator, None, None
                    )
                    .unwrap())
            );
        }
        assert!(generator.fingerprints(&[], workers).unwrap().is_empty());
        assert_eq!(
            generator.counts(&[None, None], workers).unwrap(),
            [None, None]
        );
    }
}

#[test]
fn bound_settings_alias_snapshot_and_lifetime_are_canonical() {
    let generator = TopologicalTorsionFingerprintGenerator::new(None, None).unwrap();
    let clone = generator.clone();
    let mut settings = generator.settings();
    let alias = clone.settings();
    let before = settings.params().unwrap();
    settings.set_torsion_atom_count(3).unwrap();
    settings.set_fp_size(777).unwrap();
    settings.set_count_bounds(vec![1, 3, 5]).unwrap();
    assert_eq!(alias.torsion_atom_count().unwrap(), 3);
    assert_eq!(alias.fp_size().unwrap(), 777);
    assert_eq!(alias.count_bounds().unwrap(), [1, 3, 5]);
    assert_eq!(before.torsion_atom_count, 4);
    assert_eq!(before.fp_size, 2048);
    let molecule = Molecule::from_smiles("CCCCO").unwrap();
    assert_eq!(
        molecule
            .topological_torsion_count_fingerprint_with_generator(&generator, None, None)
            .unwrap(),
        molecule
            .topological_torsion_count_fingerprint_with_generator(&clone, None, None)
            .unwrap()
    );
    settings.set_include_chirality(true).unwrap();
    let json = generator.to_json().unwrap();
    assert!(json.contains("\"includeChirality\":\"true\""));
    assert!(json.contains("\"includeChirality\":\"false\""));
    drop(generator);
    drop(clone);
    settings.set_fp_size(1000).unwrap();
    assert_eq!(alias.fp_size().unwrap(), 1000);
}

#[test]
fn source_json_restoration_and_empty_call_roots_reset_owned_output() {
    let generator = TopologicalTorsionFingerprintGenerator::new(None, None).unwrap();
    let restored =
        TopologicalTorsionFingerprintGenerator::from_json(&generator.to_json().unwrap()).unwrap();
    assert_eq!(
        generator.info_string().unwrap(),
        restored.info_string().unwrap()
    );
    assert_eq!(generator.to_json().unwrap(), restored.to_json().unwrap());
    let molecule = Molecule::from_smiles("CCCCO").unwrap();
    let call = TopologicalTorsionCallParams {
        custom_atom_invariants: Some(vec![17, 18, 19, 20, 21]),
        ..Default::default()
    };
    assert_eq!(
        molecule
            .topological_torsion_count_fingerprint_with_generator(&generator, Some(&call), None)
            .unwrap(),
        molecule
            .topological_torsion_count_fingerprint_with_generator(&restored, Some(&call), None)
            .unwrap()
    );
    let mut output = FingerprintAdditionalOutput::new();
    output.allocate_atom_counts();
    output.allocate_atom_to_bits();
    output.allocate_bit_paths();
    molecule
        .topological_torsion_sparse_count_fingerprint_with_generator(
            &generator,
            None,
            Some(&mut output),
        )
        .unwrap();
    assert_eq!(
        output
            .bit_paths()
            .unwrap()
            .keys()
            .copied()
            .collect::<Vec<_>>(),
        [4437590048, 12893306913]
    );
    let empty = TopologicalTorsionCallParams {
        from_atoms: Some(vec![]),
        ..Default::default()
    };
    molecule
        .topological_torsion_sparse_count_fingerprint_with_generator(
            &generator,
            Some(&empty),
            Some(&mut output),
        )
        .unwrap();
    assert!(output.bit_paths().unwrap().is_empty());
    assert_eq!(output.atom_counts().unwrap(), [0; 5]);
    assert!(output.atom_to_bits().unwrap().iter().all(Vec::is_empty));
    assert!(TopologicalTorsionCallParams::default().from_atoms.is_none());
}

#[test]
fn live_dense_source_bounds_error_keeps_state_usable_and_intmin_structured() {
    let generator = TopologicalTorsionFingerprintGenerator::new(None, None).unwrap();
    let molecule = Molecule::from_smiles("CCCC").unwrap();
    let mut settings = generator.settings();
    settings.set_count_bounds(vec![]).unwrap();
    let error = generator
        .fingerprints(&[None, Some(&molecule)], 2)
        .unwrap_err();
    assert!(error.to_string().contains("Count bounds are empty"));
    assert!(std::error::Error::source(&error).is_some());
    assert!(generator.sparse_counts(&[Some(&molecule)], 2).is_ok());
    settings.set_count_bounds(vec![1, 2, 4, 8]).unwrap();
    assert_eq!(
        generator.fingerprints(&[Some(&molecule)], 2).unwrap()[0]
            .as_ref()
            .unwrap()
            .on_bits(),
        [60]
    );
    assert!(
        generator
            .counts(&[Some(&molecule)], i32::MIN)
            .unwrap_err()
            .to_string()
            .contains("INT_MIN")
    );
}
