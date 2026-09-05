use cosmolkit_core::{
    FingerprintError, Molecule, PATTERN_FINGERPRINT_VERSION, PatternFingerprintParams, pattern_fingerprint,
};

#[test]
fn pattern_fingerprint_defaults_version_and_exact_ethane_bits_are_public() {
    let params = PatternFingerprintParams::default();
    assert_eq!(params.n_bits, 2048);
    assert!(!params.tautomeric);
    assert_eq!(PATTERN_FINGERPRINT_VERSION, "1.0.0");

    let molecule = Molecule::from_smiles("CC").expect("ethane");
    let fingerprint = molecule.pattern_fingerprint(&params).expect("Pattern fingerprint");
    assert_eq!(fingerprint.n_bits(), 2048);
    assert_eq!(fingerprint.on_bits(), vec![429, 778, 1022, 1061, 1236, 1295]);
}

#[test]
fn pattern_fingerprint_free_function_method_and_repeated_calls_are_value_style() {
    let molecule = Molecule::from_smiles("c1ccccc1O").expect("phenol");
    let source_snapshot = molecule.clone();
    let params = PatternFingerprintParams {
        n_bits: 997,
        tautomeric: true,
    };
    let params_snapshot = params;

    let from_method = molecule
        .pattern_fingerprint(&params)
        .expect("method Pattern fingerprint");
    let from_function = pattern_fingerprint(&molecule, &params).expect("free Pattern fingerprint");
    let repeated = molecule
        .pattern_fingerprint(&params)
        .expect("repeated Pattern fingerprint");

    assert_eq!(from_method, from_function);
    assert_eq!(from_method, repeated);
    assert_eq!(from_method.n_bits(), 997);
    assert_eq!(molecule, source_snapshot);
    assert_eq!(params, params_snapshot);
}

#[test]
fn pattern_fingerprint_supports_boundary_widths_and_rejects_zero() {
    let molecule = Molecule::from_smiles("CCC").expect("propane");
    for n_bits in [1, 63, 64, 65, 127, 2048, 4093] {
        let fingerprint = molecule
            .pattern_fingerprint(&PatternFingerprintParams {
                n_bits,
                tautomeric: false,
            })
            .expect("boundary Pattern width");
        assert_eq!(fingerprint.n_bits(), n_bits);
        assert!(fingerprint.on_bits().iter().all(|&bit| bit < n_bits));
    }

    assert_eq!(
        molecule.pattern_fingerprint(&PatternFingerprintParams {
            n_bits: 0,
            tautomeric: false,
        }),
        Err(FingerprintError::EmptyFingerprint)
    );
}

#[test]
fn pattern_fingerprint_empty_and_query_molecules_are_deterministic() {
    let empty = Molecule::new();
    let empty_fingerprint = empty
        .pattern_fingerprint(&PatternFingerprintParams::default())
        .expect("empty Pattern fingerprint");
    assert!(empty_fingerprint.on_bits().is_empty());

    let query =
        cosmolkit_core::mol_from_smarts("C~C", &cosmolkit_core::SmartsParseParams::default()).expect("query molecule");
    let first = query
        .pattern_fingerprint(&PatternFingerprintParams::default())
        .expect("query Pattern fingerprint");
    let second = query
        .pattern_fingerprint(&PatternFingerprintParams::default())
        .expect("repeated query Pattern fingerprint");
    assert_eq!(first, second);
}
