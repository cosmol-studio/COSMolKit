//! Source-backed public reusable Morgan regression tests.
#![cfg(feature = "cap-fingerprints")]
use cosmolkit::{
    FingerprintAdditionalOutput, FingerprintPreparationError, Molecule,
    MorganAtomInvariantsGenerator, MorganBondInvariantsGenerator, MorganCallParams,
    MorganFingerprintGenerator, MorganParams, MorganReadError,
};
use serde_json::{Value, json};
use std::collections::BTreeMap;
fn default_generator() -> MorganFingerprintGenerator {
    MorganFingerprintGenerator::new(None, None, None).unwrap()
}
fn complete_output() -> FingerprintAdditionalOutput {
    let mut output = FingerprintAdditionalOutput::new();
    output.allocate_atom_counts();
    output.allocate_atom_to_bits();
    output.allocate_bit_info_map();
    output.allocate_bit_paths();
    output
}
#[test]
fn fixed_source_cco_full_four_outputs_and_complete_metadata() {
    let molecule = Molecule::from_smiles("CCO").unwrap();
    let generator = default_generator();
    assert_eq!(
        molecule
            .morgan_sparse_count_fingerprint_with_generator(&generator, None, None)
            .unwrap()
            .nonzero_elements(),
        &BTreeMap::from([
            (864662311, 1),
            (1535166686, 1),
            (2245384272, 1),
            (2246728737, 1),
            (3542456614, 1),
            (4018048386, 1)
        ])
    );
    assert_eq!(
        molecule
            .morgan_sparse_fingerprint_with_generator(&generator, None, None)
            .unwrap()
            .on_bits(),
        [
            -2049583024,
            -2048238559,
            -752510682,
            -276918910,
            864662311,
            1535166686
        ]
    );
    assert_eq!(
        molecule
            .morgan_count_fingerprint_with_generator(&generator, None, None)
            .unwrap()
            .nonzero_elements(),
        &BTreeMap::from([(80, 1), (222, 1), (294, 1), (807, 1), (1057, 1), (1410, 1)])
    );
    assert_eq!(
        molecule
            .morgan_fingerprint_with_generator(&generator, None, None)
            .unwrap()
            .on_bits(),
        [80, 222, 294, 807, 1057, 1410]
    );
    assert_eq!(
        generator.info_string().unwrap(),
        "Common arguments : countSimulation=0 fpSize=2048 bitsPerFeature=1 includeChirality=0 --- MorganArguments onlyNonzeroInvariants=0 radius=3 --- MorganEnvironmentGenerator --- MorganInvariantGenerator includeRingMembership=1 --- MorganInvariantGenerator useBondTypes=1 useChirality=0"
    );
    let value: Value = serde_json::from_str(&generator.to_json().unwrap()).unwrap();
    assert_eq!(
        value,
        json!({"name":"FingerprintGenerator","fingerprintArguments":{"type":"MorganArguments","onlyNonzeroInvariants":"false","radius":"3","countSimulation":"false","fpSize":"2048","numBitsPerFeature":"1","includeChirality":"false","countBounds":["1","2","4","8"]},"atomEnvironmentGenerator":{"type":"MorganEnvGenerator"},"atomInvariantsGenerator":{"type":"MorganAtomInvGenerator","includeRingMembership":"true"},"bondInvariantsGenerator":{"type":"MorganBondInvGenerator","useBondTypes":"true","useChirality":"false"}})
    );
}
#[test]
fn settings_alias_lifetime_and_immutable_snapshot_cover_all_live_fields() {
    let generator = default_generator();
    let alias = generator.clone();
    let mut settings = generator.settings();
    let another = alias.settings();
    let before = settings.params().unwrap();
    settings.set_radius(0).unwrap();
    settings.set_only_nonzero_invariants(true).unwrap();
    settings.set_include_redundant_environments(true).unwrap();
    settings.set_include_chirality(true).unwrap();
    settings.set_count_simulation(true).unwrap();
    settings.set_fp_size(1000).unwrap();
    settings.set_bits_per_feature(2).unwrap();
    settings.set_count_bounds(vec![1, 3, 5]).unwrap();
    assert_eq!(another.radius().unwrap(), 0);
    assert!(another.only_nonzero_invariants().unwrap());
    assert!(another.include_redundant_environments().unwrap());
    assert!(another.include_chirality().unwrap());
    assert!(another.count_simulation().unwrap());
    assert_eq!(another.fp_size().unwrap(), 1000);
    assert_eq!(another.bits_per_feature().unwrap(), 2);
    assert_eq!(another.count_bounds().unwrap(), [1, 3, 5]);
    assert_eq!(before, MorganParams::default());
    let mut owned = another.count_bounds().unwrap();
    owned.push(99);
    assert_eq!(another.count_bounds().unwrap(), [1, 3, 5]);
    let value: Value = serde_json::from_str(&generator.to_json().unwrap()).unwrap();
    assert_eq!(value["fingerprintArguments"]["includeChirality"], "true");
    assert_eq!(value["bondInvariantsGenerator"]["useChirality"], "false");
    drop(generator);
    drop(alias);
    settings.set_radius(1).unwrap();
    assert_eq!(another.radius().unwrap(), 1);
}
#[test]
fn explicit_source_provider_precedence_and_independent_query_copy_lifetime() {
    let ring = Molecule::from_smiles("C1CCCCC1").unwrap();
    let params = MorganParams {
        include_ring_membership: false,
        ..Default::default()
    };
    let default_false = MorganFingerprintGenerator::new(Some(&params), None, None).unwrap();
    let atom = MorganAtomInvariantsGenerator::connectivity(true);
    let bond = MorganBondInvariantsGenerator::new(false, true);
    let explicit =
        MorganFingerprintGenerator::new(Some(&params), Some(&atom), Some(&bond)).unwrap();
    assert!(!bond.use_bond_types());
    assert!(bond.include_chirality());
    let value: Value = serde_json::from_str(&explicit.to_json().unwrap()).unwrap();
    assert_eq!(
        value["atomInvariantsGenerator"]["includeRingMembership"],
        "true"
    );
    assert_eq!(value["bondInvariantsGenerator"]["useBondTypes"], "false");
    assert_eq!(value["bondInvariantsGenerator"]["useChirality"], "true");
    assert_ne!(
        ring.morgan_sparse_count_fingerprint_with_generator(&default_false, None, None)
            .unwrap(),
        ring.morgan_sparse_count_fingerprint_with_generator(&explicit, None, None)
            .unwrap()
    );
    let molecule = Molecule::from_smiles("CCO").unwrap();
    let patterns = vec![
        cosmolkit::search::from_smarts("[C]").unwrap(),
        cosmolkit::search::from_smarts("[O]").unwrap(),
    ];
    let provider = MorganAtomInvariantsGenerator::features(Some(patterns));
    let params = MorganParams {
        radius: 0,
        ..Default::default()
    };
    let generator = MorganFingerprintGenerator::new(Some(&params), Some(&provider), None).unwrap();
    drop(provider);
    assert_eq!(
        molecule
            .morgan_sparse_count_fingerprint_with_generator(&generator, None, None)
            .unwrap()
            .nonzero_elements(),
        &BTreeMap::from([(1, 2), (2, 1)])
    );
    let empty = MorganAtomInvariantsGenerator::features(Some(vec![]));
    let generator = MorganFingerprintGenerator::new(Some(&params), Some(&empty), None).unwrap();
    assert_eq!(
        molecule
            .morgan_sparse_count_fingerprint_with_generator(&generator, None, None)
            .unwrap()
            .nonzero_elements(),
        &BTreeMap::from([(0, 3)])
    );
    let call = MorganCallParams {
        custom_atom_invariants: Some(vec![11, 12, 13]),
        ..Default::default()
    };
    assert_eq!(
        molecule
            .morgan_sparse_count_fingerprint_with_generator(&generator, Some(&call), None)
            .unwrap()
            .nonzero_elements(),
        &BTreeMap::from([(11, 1), (12, 1), (13, 1)])
    );
}
#[test]
fn full_raw_additional_output_and_present_empty_roots_reinitialize_owner_values() {
    let molecule = Molecule::from_smiles("CCO").unwrap();
    let generator = default_generator();
    let mut output = complete_output();
    molecule
        .morgan_sparse_count_fingerprint_with_generator(&generator, None, Some(&mut output))
        .unwrap();
    assert_eq!(output.atom_counts().unwrap(), [2, 2, 2]);
    assert_eq!(
        output.atom_to_bits().unwrap(),
        [
            vec![2246728737, 3542456614],
            vec![2245384272, 4018048386],
            vec![864662311, 1535166686]
        ]
    );
    assert_eq!(
        output.bit_info_map().unwrap(),
        &BTreeMap::from([
            (864662311, vec![(2, 0)]),
            (1535166686, vec![(2, 1)]),
            (2245384272, vec![(1, 0)]),
            (2246728737, vec![(0, 0)]),
            (3542456614, vec![(0, 1)]),
            (4018048386, vec![(1, 1)])
        ])
    );
    assert!(output.bit_paths().unwrap().is_empty());
    let empty = MorganCallParams {
        from_atoms: Some(vec![]),
        ..Default::default()
    };
    molecule
        .morgan_sparse_count_fingerprint_with_generator(&generator, Some(&empty), Some(&mut output))
        .unwrap();
    assert_eq!(output.atom_counts().unwrap(), [0; 3]);
    assert!(output.atom_to_bits().unwrap().iter().all(Vec::is_empty));
    assert!(output.bit_info_map().unwrap().is_empty());
    assert!(output.bit_paths().unwrap().is_empty());
    assert!(MorganCallParams::default().from_atoms.is_none());
}
#[test]
fn source_json_missing_providers_keep_null_preconditions_and_custom_precedence() {
    let generator = default_generator();
    let mut value: Value = serde_json::from_str(&generator.to_json().unwrap()).unwrap();
    value
        .as_object_mut()
        .unwrap()
        .remove("atomInvariantsGenerator");
    value
        .as_object_mut()
        .unwrap()
        .remove("bondInvariantsGenerator");
    let restored = MorganFingerprintGenerator::from_json(&value.to_string()).unwrap();
    let molecule = Molecule::from_smiles("CCO").unwrap();
    let error = molecule
        .morgan_sparse_count_fingerprint_with_generator(&restored, None, None)
        .unwrap_err();
    assert!(error.to_string().contains("atom invariants"));
    assert!(std::error::Error::source(&error).is_some());
    let atoms = MorganCallParams {
        custom_atom_invariants: Some(vec![11, 12, 13]),
        ..Default::default()
    };
    let error = molecule
        .morgan_sparse_count_fingerprint_with_generator(&restored, Some(&atoms), None)
        .unwrap_err();
    assert!(error.to_string().contains("bond invariants"));
    let full = MorganCallParams {
        custom_bond_invariants: Some(vec![1, 1]),
        ..atoms
    };
    assert!(
        molecule
            .morgan_sparse_count_fingerprint_with_generator(&restored, Some(&full), None)
            .is_ok()
    );
    assert!(restored.sparse_counts(&[None, Some(&molecule)], 3).is_err());
}
#[test]
fn source_json_arguments_and_provider_flags_survive_restore_without_factory_validation() {
    let generator = default_generator();
    let mut settings = generator.settings();
    settings.set_include_chirality(true).unwrap();
    settings.set_include_redundant_environments(true).unwrap();
    settings.set_count_bounds(vec![]).unwrap();
    settings.set_bits_per_feature(0).unwrap();
    let restored = MorganFingerprintGenerator::from_json(&generator.to_json().unwrap()).unwrap();
    assert!(restored.settings().include_chirality().unwrap());
    assert!(
        !restored
            .settings()
            .include_redundant_environments()
            .unwrap()
    );
    assert_eq!(restored.settings().bits_per_feature().unwrap(), 0);
    assert!(restored.settings().count_bounds().unwrap().is_empty());
    assert!(restored.info_string().unwrap().contains("useChirality=0"));
    let molecule = Molecule::from_smiles("CCO").unwrap();
    assert_eq!(
        molecule
            .morgan_sparse_count_fingerprint_with_generator(&restored, None, None)
            .unwrap()
            .nonzero_elements()
            .len(),
        6
    );
    let json = generator
        .to_json()
        .unwrap()
        .replace("\"fpSize\":\"2048\"", "\"fpSize\":1,\"fpSize\":2")
        .replace(
            "\"includeChirality\":\"true\"",
            "\"includeChirality\":\"invalid\"",
        );
    let restored = MorganFingerprintGenerator::from_json(&json).unwrap();
    assert_eq!(restored.settings().fp_size().unwrap(), 1);
    assert!(!restored.settings().include_chirality().unwrap());
}
#[test]
fn bulk_four_forms_preserve_source_order_none_slots_and_empty_input() {
    let generator = default_generator();
    let molecules = ["CCO", "C", "C1CCCCC1"].map(|s| Molecule::from_smiles(s).unwrap());
    let rows = [
        Some(&molecules[0]),
        None,
        Some(&molecules[1]),
        Some(&molecules[2]),
        None,
    ];
    for workers in [1, 3, 7] {
        let bits = generator.fingerprints(&rows, workers).unwrap();
        let counts = generator.counts(&rows, workers).unwrap();
        let sparse = generator.sparse_fingerprints(&rows, workers).unwrap();
        let raw = generator.sparse_counts(&rows, workers).unwrap();
        for (i, m) in rows.iter().enumerate() {
            assert_eq!(
                bits[i],
                m.map(|m| m
                    .morgan_fingerprint_with_generator(&generator, None, None)
                    .unwrap())
            );
            assert_eq!(
                counts[i],
                m.map(|m| m
                    .morgan_count_fingerprint_with_generator(&generator, None, None)
                    .unwrap())
            );
            assert_eq!(
                sparse[i],
                m.map(|m| m
                    .morgan_sparse_fingerprint_with_generator(&generator, None, None)
                    .unwrap())
            );
            assert_eq!(
                raw[i],
                m.map(|m| m
                    .morgan_sparse_count_fingerprint_with_generator(&generator, None, None)
                    .unwrap())
            );
        }
        assert!(generator.fingerprints(&[], workers).unwrap().is_empty());
        assert_eq!(
            generator.counts(&[None, None], workers).unwrap(),
            [None, None]
        );
    }
    assert!(
        generator
            .counts(&[Some(&molecules[0])], i32::MIN)
            .unwrap_err()
            .to_string()
            .contains("INT_MIN")
    );
}
#[test]
fn shared_preparation_category_and_registry_cover_the_canonical_boundary() {
    let generator = default_generator();
    let error = Molecule::new()
        .morgan_fingerprint_with_generator(&generator, None, None)
        .unwrap_err();
    assert!(matches!(
        error,
        MorganReadError::Preparation(FingerprintPreparationError::MissingPreparedValence)
    ));
    let cause = std::error::Error::source(&error).unwrap();
    assert!(
        cause
            .downcast_ref::<FingerprintPreparationError>()
            .is_some()
    );
    for name in [
        "types.MorganAtomInvariantsGenerator",
        "types.MorganBondInvariantsGenerator",
        "types.MorganFingerprintGenerator",
        "types.MorganSettings",
        "types.MorganCallParams",
        "Molecule.morgan_fingerprint_with_generator",
        "Molecule.morgan_sparse_fingerprint_with_generator",
        "Molecule.morgan_count_fingerprint_with_generator",
        "Molecule.morgan_sparse_count_fingerprint_with_generator",
    ] {
        assert!(
            cosmolkit::binding_contract::BINDING_CONTRACT
                .iter()
                .any(|entry| entry.semantic_id == name),
            "{name}"
        );
    }
}
