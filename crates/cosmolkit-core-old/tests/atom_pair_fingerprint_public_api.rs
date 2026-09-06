use cosmolkit_core::{
    AdditionalOutput, AtomPairFingerprintGenerator, AtomPairFingerprintParams, FingerprintError,
    FingerprintFuncArguments, Molecule, MorganFingerprintParams, atom_pair_fingerprint,
    atom_pair_fingerprint_with_output,
};

#[test]
fn defaults_and_all_result_forms_are_available_through_one_generator() {
    let params = AtomPairFingerprintParams::default();
    assert_eq!(params.n_bits, 2048);
    assert_eq!(params.min_distance, 1);
    assert_eq!(params.max_distance, 30);
    assert!(params.use_2d);
    assert!(!params.use_chirality);
    assert!(params.count_simulation);
    assert_eq!(params.count_bounds, [1, 2, 4, 8]);
    assert_eq!(params.num_bits_per_feature, 1);
    assert_eq!(params.conformer_id, -1);

    let molecule = Molecule::from_smiles("CCCO").unwrap();
    let generator = AtomPairFingerprintGenerator::new(&params).unwrap();
    let sparse_count = generator
        .sparse_count_fingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .unwrap();
    let count = generator
        .count_fingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .unwrap();
    let sparse_bits = generator
        .sparse_bit_fingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .unwrap();
    let bits = generator
        .fingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .unwrap();

    assert_eq!(sparse_count.size(), 1 << 23);
    assert_eq!(count.size(), 2048);
    assert_eq!(sparse_bits.size(), 1 << 23);
    assert_eq!(bits.n_bits(), 2048);
    assert!(!sparse_count.nonzero_elements().is_empty());
    assert!(!count.nonzero_elements().is_empty());
    assert!(!sparse_bits.on_bits().is_empty());
    assert!(!bits.on_bits().is_empty());
}

#[test]
fn molecule_method_free_function_and_output_wrapper_have_value_semantics() {
    let molecule = Molecule::from_smiles("CCCO").unwrap();
    let source_snapshot = molecule.clone();
    let params = AtomPairFingerprintParams {
        collect_additional_output: true,
        ..Default::default()
    };
    let params_snapshot = params.clone();

    let from_method = molecule.atom_pair_fingerprint(&params).unwrap();
    let from_function = atom_pair_fingerprint(&molecule, &params).unwrap();
    let with_output = molecule.atom_pair_fingerprint_with_output(&params).unwrap();
    let free_output = atom_pair_fingerprint_with_output(&molecule, &params).unwrap();

    assert_eq!(from_method, from_function);
    assert_eq!(with_output, free_output);
    assert_eq!(with_output.fingerprint, from_method);
    let output = with_output.additional_output.unwrap();
    assert_eq!(output.atom_counts.as_ref().map(Vec::len), Some(4));
    assert_eq!(output.atom_to_bits.as_ref().map(Vec::len), Some(4));
    assert!(output.bit_info_map.is_some());
    assert!(output.atoms_per_bit.is_some());
    assert!(output.bit_paths.is_none());
    assert_eq!(molecule, source_snapshot);
    assert_eq!(params, params_snapshot);
}

#[test]
fn generator_clone_json_and_repeated_calls_are_deterministic() {
    let molecule = Molecule::from_smiles("c1ccccc1O").unwrap();
    let params = AtomPairFingerprintParams {
        count_simulation: false,
        num_bits_per_feature: 2,
        ..Default::default()
    };
    let generator = AtomPairFingerprintGenerator::new(&params).unwrap();
    let cloned = generator.clone();
    assert_eq!(generator, cloned);

    let restored = AtomPairFingerprintGenerator::from_json(&generator.to_json()).unwrap();
    assert_eq!(restored.to_json(), generator.to_json());
    assert_eq!(restored.info_string(), generator.info_string());

    let first = generator
        .fingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .unwrap();
    let second = generator
        .fingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .unwrap();
    let restored_result = restored
        .fingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .unwrap();
    assert_eq!(first, second);
    assert_eq!(first, restored_result);
}

#[test]
fn invalid_public_arguments_return_typed_errors() {
    let molecule = Molecule::from_smiles("CCO").unwrap();
    let zero_size = AtomPairFingerprintParams {
        n_bits: 0,
        ..Default::default()
    };
    assert!(matches!(
        molecule.atom_pair_fingerprint(&zero_size),
        Err(FingerprintError::EmptyFingerprint)
    ));

    let bad_distances = AtomPairFingerprintParams {
        min_distance: 4,
        max_distance: 3,
        ..Default::default()
    };
    assert!(matches!(
        AtomPairFingerprintGenerator::new(&bad_distances),
        Err(FingerprintError::InvalidArguments { .. })
    ));

    let three_dimensional = AtomPairFingerprintParams {
        use_2d: false,
        ..Default::default()
    };
    let generator =
        AtomPairFingerprintGenerator::new(&three_dimensional).expect("3D generator parameters");
    assert!(matches!(
        generator.fingerprint(&molecule, &mut FingerprintFuncArguments::default()),
        Err(FingerprintError::AtomPair { .. })
    ));

    let generator = AtomPairFingerprintGenerator::new(&AtomPairFingerprintParams::default())
        .expect("default generator");
    let mut bad_indices = FingerprintFuncArguments {
        from_atoms: Some(vec![molecule.num_atoms()]),
        ..Default::default()
    };
    let filtered = generator
        .fingerprint(&molecule, &mut bad_indices)
        .expect("source-compatible unmatched fromAtoms filter");
    assert!(filtered.on_bits().is_empty());
}

#[test]
fn mixed_family_call_order_does_not_change_atom_pair_results() {
    let molecule = Molecule::from_smiles("CC(O)C1=CC=CC=C1").unwrap();
    let atom_pair_params = AtomPairFingerprintParams::default();
    let before = molecule.atom_pair_fingerprint(&atom_pair_params).unwrap();
    let _morgan = molecule
        .morgan_fingerprint(&MorganFingerprintParams::default())
        .unwrap();
    let after = molecule.atom_pair_fingerprint(&atom_pair_params).unwrap();
    assert_eq!(before, after);
}

#[test]
fn public_surface_does_not_add_rdkit_style_atom_pair_duplicates() {
    let atom_pair_source = include_str!("../src/properties/fingerprint/atom_pair.rs");
    let crate_exports = include_str!("../src/lib.rs");
    for forbidden in [
        "pub fn getAtomPairGenerator",
        "pub fn get_atom_pair_generator",
        "pub fn getAtomPairFingerprint",
        "pub fn getHashedAtomPairFingerprint",
    ] {
        assert!(!atom_pair_source.contains(forbidden), "{forbidden}");
        assert!(!crate_exports.contains(forbidden), "{forbidden}");
    }
}

#[test]
fn caller_selected_additional_output_allocations_are_preserved() {
    let molecule = Molecule::from_smiles("CCO").unwrap();
    let generator = AtomPairFingerprintGenerator::new(&AtomPairFingerprintParams::default())
        .expect("default generator");
    let mut output = AdditionalOutput::new();
    output.allocate_atom_counts();
    output.allocate_bit_info_map();
    let mut arguments = FingerprintFuncArguments {
        additional_output: Some(output),
        ..Default::default()
    };
    generator
        .fingerprint(&molecule, &mut arguments)
        .expect("AtomPair fingerprint with selected provenance");
    let output = arguments.additional_output.unwrap();
    assert!(output.atom_counts.is_some());
    assert!(output.bit_info_map.is_some());
    assert!(output.atom_to_bits.is_none());
    assert!(output.atoms_per_bit.is_none());
    assert!(output.bit_paths.is_none());
}
