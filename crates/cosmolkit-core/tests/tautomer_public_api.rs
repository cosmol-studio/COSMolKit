use std::sync::{
    Arc,
    atomic::{AtomicUsize, Ordering},
};

use cosmolkit_core::{
    ENUMERATE_TAUTOMERS_WITH_OPTIONS_SPEC, MaccsFingerprintParams, Molecule, OperationError,
    SupportStatus, TAUTOMER_ENUMERATION_FEATURE, TautomerCanonicalizationError, TautomerCatalog,
    TautomerCatalogError, TautomerEnumeration, TautomerEnumerationCallback,
    TautomerEnumerationError, TautomerEnumerationStatus, TautomerEnumerator, TautomerOptions,
    TautomerRunError, TautomerScore, TautomerScoreError, TautomerScoreTerm, TautomerTransform,
    TautomerTransformError, calc_exact_mol_wt, default_tautomer_score_terms, mol_from_binary,
    mol_to_binary, score_tautomer,
};

fn acetone() -> Molecule {
    Molecule::from_smiles("CC(C)=O").expect("parse acetone")
}

fn assert_send_sync<T: Send + Sync>() {}

#[test]
fn root_exports_defaults_and_registered_molecule_entry_points_are_discoverable() {
    assert_send_sync::<TautomerCatalog>();
    assert_send_sync::<TautomerOptions>();
    assert_send_sync::<TautomerEnumeration>();
    assert_send_sync::<TautomerEnumerator<'static>>();

    let options = TautomerOptions::default();
    assert_eq!(options.max_tautomers(), 1000);
    assert_eq!(options.max_transforms(), 1000);
    assert!(options.remove_sp3_stereo());
    assert!(options.remove_bond_stereo());
    assert!(options.remove_isotopic_hydrogens());
    assert!(options.reassign_stereo());

    assert_eq!(TAUTOMER_ENUMERATION_FEATURE.name, "tautomer.enumeration");
    assert_eq!(
        TAUTOMER_ENUMERATION_FEATURE.status,
        SupportStatus::Experimental
    );
    assert_eq!(
        ENUMERATE_TAUTOMERS_WITH_OPTIONS_SPEC.method,
        "enumerate_tautomers_with_options"
    );

    let source = acetone();
    let default_result = source.enumerate_tautomers().expect("default entry point");
    let enumerator = TautomerEnumerator::new();
    let configured_result = source
        .enumerate_tautomers_with_options(&enumerator)
        .expect("configured molecule entry point");
    let enumerator_result = enumerator
        .enumerate(&source)
        .expect("enumerator entry point");

    assert_eq!(default_result, configured_result);
    assert_eq!(configured_result, enumerator_result);
}

#[test]
fn options_result_sequence_and_canonicalization_form_one_value_style_surface() {
    let options = TautomerOptions::default()
        .with_max_tautomers(17)
        .with_max_transforms(23)
        .with_remove_sp3_stereo(false)
        .with_remove_bond_stereo(false)
        .with_remove_isotopic_hydrogens(false)
        .with_reassign_stereo(false);
    let mut enumerator = TautomerEnumerator::from_options(options);
    assert_eq!(enumerator.options(), options);
    assert_eq!(enumerator.max_tautomers(), 17);
    assert_eq!(enumerator.max_transforms(), 23);
    assert!(!enumerator.remove_sp3_stereo());
    assert!(!enumerator.remove_bond_stereo());
    assert!(!enumerator.remove_isotopic_hydrogens());
    assert!(!enumerator.reassign_stereo());

    enumerator.set_max_tautomers(1000);
    enumerator.set_max_transforms(1000);
    enumerator.set_remove_sp3_stereo(true);
    enumerator.set_remove_bond_stereo(true);
    enumerator.set_remove_isotopic_hydrogens(true);
    enumerator.set_reassign_stereo(true);

    let source = acetone().with_prop("record", "source");
    let before = source.clone();
    let result = enumerator.enumerate(&source).expect("enumerate");
    assert_eq!(result.status(), TautomerEnumerationStatus::Completed);
    assert_eq!(result.len(), 2);
    assert_eq!(result.canonical_smiles(), ["C=C(C)O", "CC(C)=O"]);
    assert_eq!(result.get(2), None);
    assert!(matches!(
        result.at(2),
        Err(TautomerEnumerationError::IndexOutOfRange { index: 2, len: 2 })
    ));
    assert_eq!(result.iter().count(), 2);
    assert_eq!(result.iter().rev().count(), 2);
    assert_eq!((&result).into_iter().count(), 2);
    assert_eq!(result[0].to_smiles(true).unwrap(), "C=C(C)O");

    let selected = enumerator
        .pick_canonical(&result)
        .expect("pick canonical tautomer");
    let canonical = enumerator
        .canonicalize(&source)
        .expect("canonicalize source");
    assert_eq!(
        selected.to_smiles(true).unwrap(),
        canonical.to_smiles(true).unwrap()
    );
    assert_eq!(source, before);
    assert_eq!(canonical.prop("record"), Some("source"));
}

#[test]
fn custom_catalog_and_custom_scoring_reuse_the_public_core_types() {
    let definitions = [(
        "1,3 (thio)keto/enol f",
        "[CX4!H0R{0-2}]-[C;z{1-2}]=[O,S,Se,Te;X1]",
        "",
        "",
    )];
    let catalog = TautomerCatalog::from_data(&definitions).expect("custom catalog");
    assert_eq!(catalog.transforms().len(), 1);
    assert_eq!(catalog.transform(0).unwrap().name(), definitions[0].0);
    let enumerator = TautomerEnumerator::from_catalog_and_options(
        catalog,
        TautomerOptions::default().with_max_transforms(100),
    );
    assert_eq!(enumerator.catalog().transforms().len(), 1);

    let source = acetone();
    let result = source.enumerate_tautomers().expect("enumerate for scoring");
    let custom = enumerator
        .pick_canonical_with(&result, |candidate| {
            Ok(if candidate.to_smiles(true)? == "C=C(C)O" {
                10
            } else {
                -10
            })
        })
        .expect("custom canonical scorer");
    assert_eq!(custom.to_smiles(true).unwrap(), "C=C(C)O");

    let term = TautomerScoreTerm::new("carbonyl", "[C]=[O]", 7);
    assert_eq!(term.name(), "carbonyl");
    assert_eq!(term.smarts(), "[C]=[O]");
    assert_eq!(term.score(), 7);
    assert_eq!(default_tautomer_score_terms().len(), 12);
    let score: TautomerScore = score_tautomer(&custom).expect("score custom selection");
    assert_eq!(
        score.total(),
        score
            .ring()
            .wrapping_add(score.substructure())
            .wrapping_add(score.hetero_hydrogen())
    );
}

#[derive(Default)]
struct CancelImmediately {
    calls: AtomicUsize,
}

impl TautomerEnumerationCallback for CancelImmediately {
    fn should_continue(&self, source: &Molecule, result: &TautomerEnumeration) -> bool {
        assert_eq!(source.to_smiles(true).unwrap(), "CC(C)=O");
        assert_eq!(result.status(), TautomerEnumerationStatus::Completed);
        self.calls.fetch_add(1, Ordering::SeqCst);
        false
    }
}

#[test]
fn borrowed_callback_cancels_at_the_source_callback_boundary() {
    let callback = CancelImmediately::default();
    let mut enumerator = TautomerEnumerator::new();
    enumerator.set_callback(Some(&callback));

    let result = enumerator.enumerate(&acetone()).expect("canceled run");

    assert_eq!(callback.calls.load(Ordering::SeqCst), 1);
    assert_eq!(result.status(), TautomerEnumerationStatus::Canceled);
    assert_eq!(result.len(), 1);
    assert!(enumerator.callback().is_some());
}

#[test]
fn public_failures_retain_their_structured_error_categories() {
    let invalid = Molecule::from_smiles_with_sanitize("c1cccc1", false)
        .expect("parse deliberately invalid aromatic ring");
    let run_error = invalid
        .enumerate_tautomers()
        .expect_err("invalid aromatic input must fail");
    assert!(matches!(
        run_error,
        OperationError::Tautomer { source, .. }
            if matches!(*source, TautomerRunError::Kekulize(_))
    ));

    let canonical_error = TautomerEnumerator::new()
        .canonicalize(&invalid)
        .expect_err("canonicalization must retain the operation error");
    assert!(matches!(
        canonical_error,
        TautomerCanonicalizationError::Enumeration { source }
            if matches!(*source, OperationError::Tautomer { .. })
    ));

    let empty = TautomerEnumeration::default();
    assert!(matches!(
        TautomerEnumerator::new().pick_canonical(&empty),
        Err(TautomerRunError::NoCanonicalTautomer)
    ));
    assert!(matches!(
        TautomerCatalog::from_data(&[("bad", "[O]-[C]", "-", "x")]),
        Err(TautomerCatalogError::Transform(
            TautomerTransformError::ChargeSymbolNotRecognised
        ))
    ));

    fn assert_score_error(_: &TautomerScoreError) {}
    let unsanitized_phosphorus =
        Molecule::from_smiles_with_sanitize("P", false).expect("parse unsanitized phosphorus");
    let score_error =
        score_tautomer(&unsanitized_phosphorus).expect_err("phosphorus requires valence state");
    assert_score_error(&score_error);
}

#[test]
fn emitted_tautomers_compose_with_descriptors_fingerprints_and_binary_io() {
    let result = acetone().enumerate_tautomers().expect("enumerate");

    for tautomer in &result {
        assert!(calc_exact_mol_wt(tautomer, false).unwrap().is_finite());
        assert_eq!(
            tautomer
                .maccs_fingerprint(&MaccsFingerprintParams::default())
                .unwrap()
                .n_bits(),
            166
        );
        let binary = mol_to_binary(tautomer).expect("serialize tautomer");
        let restored = mol_from_binary(&binary).expect("restore tautomer");
        assert_eq!(
            restored.to_smiles(true).unwrap(),
            tautomer.to_smiles(true).unwrap()
        );
    }
}

#[test]
fn independent_parallel_runs_are_deterministic_and_share_no_mutable_state() {
    let enumerator = Arc::new(TautomerEnumerator::new());
    let source = Arc::new(acetone());
    let expected = source.enumerate_tautomers().unwrap().canonical_smiles();
    let handles = (0..8)
        .map(|_| {
            let enumerator = Arc::clone(&enumerator);
            let source = Arc::clone(&source);
            std::thread::spawn(move || enumerator.enumerate(&source).unwrap().canonical_smiles())
        })
        .collect::<Vec<_>>();

    for handle in handles {
        assert_eq!(handle.join().expect("tautomer worker panicked"), expected);
    }
}

#[test]
fn public_surface_has_no_rdkit_style_or_duplicate_in_place_names() {
    let tautomer_source = include_str!("../src/chemistry/tautomer.rs");
    let operation_source = include_str!("../src/operations/ops.rs");
    for forbidden in [
        "pub fn pickCanonical",
        "pub fn canonicalizeInPlace",
        "pub fn enumerateTautomers",
        "inplace_method: canonicalize",
        "default_inplace_method: canonicalize",
    ] {
        assert!(!tautomer_source.contains(forbidden));
        assert!(!operation_source.contains(forbidden));
    }
}

#[test]
fn root_transform_export_constructs_checked_custom_values() {
    let query =
        cosmolkit_core::mol_from_smarts("[O]-[C]", &cosmolkit_core::SmartsParseParams::default())
            .expect("compile transform query");
    let transform = TautomerTransform::new(
        "custom",
        query,
        vec![cosmolkit_core::BondOrder::Single],
        Vec::new(),
    )
    .expect("construct transform");
    assert_eq!(transform.name(), "custom");
    assert_eq!(transform.bond_types(), [cosmolkit_core::BondOrder::Single]);
}
