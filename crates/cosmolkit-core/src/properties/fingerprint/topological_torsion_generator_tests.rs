use super::{
    AtomPairAtomInvGenerator, FPType, Fingerprint, FingerprintFuncArguments, SparseBitFingerprint,
    SparseCountFingerprint, TopologicalTorsionArguments, TypedFingerprintGenerator,
    generatorFromJSON, generatorToJSON, getCountFP, getCountFPBulk, getFP, getFPBulk,
    getSparseCountFP, getSparseCountFPBulk, getSparseFP, getSparseFPBulk,
    getTopologicalTorsionFingerprint, getTopologicalTorsionGenerator,
};
use crate::{
    AtomQueryPredicate, AtomSpec, BondOrder, BondQueryPredicate, BondSpec, Element, Molecule,
    QueryNode,
};

fn sparse_total(fingerprint: &SparseCountFingerprint) -> i32 {
    fingerprint.nonzero_elements().values().sum()
}

fn scalar_sparse_counts(
    generator: &super::TopologicalTorsionFingerprintGenerator,
    molecules: &[Option<&Molecule>],
) -> Vec<Option<SparseCountFingerprint>> {
    molecules
        .iter()
        .map(|molecule| {
            molecule.map(|molecule| {
                generator
                    .getSparseCountFingerprint(molecule, &mut FingerprintFuncArguments::default())
                    .expect("sparse count")
            })
        })
        .collect()
}

fn scalar_sparse_bits(
    generator: &super::TopologicalTorsionFingerprintGenerator,
    molecules: &[Option<&Molecule>],
) -> Vec<Option<SparseBitFingerprint>> {
    molecules
        .iter()
        .map(|molecule| {
            molecule.map(|molecule| {
                generator
                    .getSparseFingerprint(molecule, &mut FingerprintFuncArguments::default())
                    .expect("sparse bit")
            })
        })
        .collect()
}

fn scalar_counts(
    generator: &super::TopologicalTorsionFingerprintGenerator,
    molecules: &[Option<&Molecule>],
) -> Vec<Option<SparseCountFingerprint>> {
    molecules
        .iter()
        .map(|molecule| {
            molecule.map(|molecule| {
                generator
                    .getCountFingerprint(molecule, &mut FingerprintFuncArguments::default())
                    .expect("count")
            })
        })
        .collect()
}

fn scalar_bits(
    generator: &super::TopologicalTorsionFingerprintGenerator,
    molecules: &[Option<&Molecule>],
) -> Vec<Option<Fingerprint>> {
    molecules
        .iter()
        .map(|molecule| {
            molecule.map(|molecule| {
                generator
                    .getFingerprint(molecule, &mut FingerprintFuncArguments::default())
                    .expect("bit")
            })
        })
        .collect()
}

#[test]
fn modern_generator_produces_all_four_scalar_vector_forms() {
    // Pinned RDKit testFingerprintGenerators.cpp::testGitHubIssue25 fixes the
    // unfolded ids and the 1000-bit folded count ids for CCCCO.
    let molecule = Molecule::from_smiles("CCCCO").expect("molecule");
    let default_generator =
        getTopologicalTorsionGenerator(&TopologicalTorsionArguments::default(), None, true)
            .expect("default generator");

    let sparse_count = default_generator
        .getSparseCountFingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .expect("sparse count");
    assert_eq!(sparse_total(&sparse_count), 2);
    assert_eq!(
        sparse_count
            .nonzero_elements()
            .keys()
            .copied()
            .collect::<Vec<_>>(),
        vec![4_437_590_048, 12_893_306_913]
    );

    let sparse_bit = default_generator
        .getSparseFingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .expect("sparse bit");
    assert_eq!(sparse_bit.size(), u64::from(u32::MAX));
    assert_eq!(sparse_bit.on_bits().len(), 2);

    let explicit_bit = default_generator
        .getFingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .expect("explicit bit");
    assert_eq!(explicit_bit.n_bits(), 2048);
    assert_eq!(explicit_bit.on_bits().len(), 2);

    let folded_arguments = TopologicalTorsionArguments::new(false, 4, true, vec![1, 2, 4, 8], 1000)
        .expect("folded arguments");
    let folded_generator =
        getTopologicalTorsionGenerator(&folded_arguments, None, true).expect("folded generator");
    let count = folded_generator
        .getCountFingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .expect("folded count");
    assert_eq!(count.size(), 1000);
    assert_eq!(sparse_total(&count), 2);
    assert_eq!(
        count.nonzero_elements().keys().copied().collect::<Vec<_>>(),
        vec![24, 288]
    );
}

#[test]
fn all_four_bulk_forms_are_ordered_parallel_scalar_equivalents() {
    let first = Molecule::from_smiles("CCCC").expect("first");
    let second = Molecule::from_smiles("CCCCO").expect("second");
    let third = Molecule::from_smiles("c1ccccc1").expect("third");
    let molecules = [Some(&first), None, Some(&second), Some(&third)];
    let generator =
        getTopologicalTorsionGenerator(&TopologicalTorsionArguments::default(), None, true)
            .expect("generator");

    let sparse_counts = scalar_sparse_counts(&generator, &molecules);
    let sparse_bits = scalar_sparse_bits(&generator, &molecules);
    let counts = scalar_counts(&generator, &molecules);
    let bits = scalar_bits(&generator, &molecules);

    assert_eq!(
        generator
            .getSparseCountFingerprints(&molecules, 2)
            .expect("bulk sparse counts"),
        sparse_counts
    );
    assert_eq!(
        generator
            .getSparseFingerprints(&molecules, 2)
            .expect("bulk sparse bits"),
        sparse_bits
    );
    assert_eq!(
        generator
            .getCountFingerprints(&molecules, 2)
            .expect("bulk counts"),
        counts
    );
    assert_eq!(
        generator.getFingerprints(&molecules, 2).expect("bulk bits"),
        bits
    );

    assert_eq!(
        getSparseCountFPBulk(&molecules, FPType::TopologicalTorsionFP)
            .expect("typed sparse count bulk"),
        sparse_counts
    );
    assert_eq!(
        getSparseFPBulk(&molecules, FPType::TopologicalTorsionFP).expect("typed sparse bit bulk"),
        sparse_bits
    );
    assert_eq!(
        getCountFPBulk(&molecules, FPType::TopologicalTorsionFP).expect("typed count bulk"),
        counts
    );
    assert_eq!(
        getFPBulk(&molecules, FPType::TopologicalTorsionFP).expect("typed bit bulk"),
        bits
    );

    assert_eq!(
        getSparseCountFP(&first, FPType::TopologicalTorsionFP).unwrap(),
        sparse_counts[0].clone().unwrap()
    );
    assert_eq!(
        getSparseFP(&first, FPType::TopologicalTorsionFP).unwrap(),
        sparse_bits[0].clone().unwrap()
    );
    assert_eq!(
        getCountFP(&first, FPType::TopologicalTorsionFP).unwrap(),
        counts[0].clone().unwrap()
    );
    assert_eq!(
        getFP(&first, FPType::TopologicalTorsionFP).unwrap(),
        bits[0].clone().unwrap()
    );
}

#[test]
fn mutable_options_change_only_the_shared_generator_configuration() {
    let molecule = Molecule::from_smiles("CCCCO.Cl").expect("molecule");
    let mut generator =
        getTopologicalTorsionGenerator(&TopologicalTorsionArguments::default(), None, true)
            .expect("generator");

    let four_atom = generator
        .getSparseCountFingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .expect("four atom");
    assert_eq!(sparse_total(&four_atom), 2);
    assert_eq!(four_atom.size(), 1_u64 << 36);

    generator.fingerprint_arguments.d_torsion_atom_count = 3;
    generator.fingerprint_arguments.df_only_shortest_paths = true;
    generator
        .fingerprint_arguments
        .fingerprint_arguments
        .d_fp_size = 1024;
    let three_atom = generator
        .getSparseCountFingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .expect("three atom");
    let folded = generator
        .getCountFingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .expect("folded");

    assert_eq!(sparse_total(&three_atom), 3);
    assert_eq!(three_atom.size(), 1_u64 << 27);
    assert_ne!(three_atom, four_atom);
    assert_eq!(folded.size(), 1024);
    assert_eq!(sparse_total(&folded), 3);

    generator
        .fingerprint_arguments
        .fingerprint_arguments
        .df_include_chirality = true;
    let chiral_three_atom = generator
        .getSparseCountFingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .expect("chiral three atom");
    assert_eq!(chiral_three_atom.size(), 1_u64 << 33);
}

#[test]
fn live_invalid_options_return_structured_errors_without_breaking_unfolded_counts() {
    let molecule = Molecule::from_smiles("CCCC").expect("molecule");
    let mut generator =
        getTopologicalTorsionGenerator(&TopologicalTorsionArguments::default(), None, true)
            .expect("generator");
    generator
        .fingerprint_arguments
        .fingerprint_arguments
        .d_count_bounds
        .clear();

    assert!(
        generator
            .getSparseCountFingerprint(&molecule, &mut FingerprintFuncArguments::default())
            .is_ok(),
        "count bounds are not used by the unfolded sparse-count form"
    );
    assert!(matches!(
        generator.getSparseFingerprint(&molecule, &mut FingerprintFuncArguments::default()),
        Err(super::FingerprintError::InvalidArguments {
            reason: "Count bounds are empty"
        })
    ));
    assert!(matches!(
        generator.getFingerprint(&molecule, &mut FingerprintFuncArguments::default()),
        Err(super::FingerprintError::InvalidArguments {
            reason: "Count bounds are empty"
        })
    ));

    generator
        .fingerprint_arguments
        .fingerprint_arguments
        .df_count_simulation = false;
    generator
        .fingerprint_arguments
        .fingerprint_arguments
        .d_fp_size = 0;
    assert!(
        generator
            .getSparseCountFingerprint(&molecule, &mut FingerprintFuncArguments::default())
            .is_ok()
    );
    assert!(
        generator
            .getSparseFingerprint(&molecule, &mut FingerprintFuncArguments::default())
            .is_ok()
    );
    assert!(matches!(
        generator.getCountFingerprint(&molecule, &mut FingerprintFuncArguments::default()),
        Err(super::FingerprintError::InvalidArguments { .. })
    ));
    assert!(matches!(
        generator.getFingerprint(&molecule, &mut FingerprintFuncArguments::default()),
        Err(super::FingerprintError::InvalidArguments { .. })
    ));

    generator.fingerprint_arguments.d_torsion_atom_count = 8;
    assert!(matches!(
        generator.getSparseCountFingerprint(&molecule, &mut FingerprintFuncArguments::default()),
        Err(super::FingerprintError::InvalidArguments {
            reason: "topological torsion result-size shift must be less than 64 bits"
        })
    ));
}

#[test]
fn json_restored_empty_count_bounds_return_structured_generation_errors() {
    let molecule = Molecule::from_smiles("CCCC").expect("molecule");
    let mut arguments = TopologicalTorsionArguments::default();
    arguments
        .fromJSON(r#"{"countSimulation":true,"countBounds":[]}"#)
        .expect("RDKit accepts the restored live option state");
    let generator = getTopologicalTorsionGenerator(&arguments, None, true).expect("generator");

    assert!(matches!(
        generator.getSparseFingerprint(&molecule, &mut FingerprintFuncArguments::default()),
        Err(super::FingerprintError::InvalidArguments {
            reason: "Count bounds are empty"
        })
    ));
    assert!(matches!(
        generator.getFingerprint(&molecule, &mut FingerprintFuncArguments::default()),
        Err(super::FingerprintError::InvalidArguments {
            reason: "Count bounds are empty"
        })
    ));
}

#[test]
#[allow(deprecated)]
fn legacy_unfolded_api_uses_source_atom_codes_for_low_code_endpoints() {
    let bccc = Molecule::from_smiles("BCCC").expect("BCCC");
    let bcccc = Molecule::from_smiles("BCCCC").expect("BCCCC");

    let bccc_fingerprint =
        getTopologicalTorsionFingerprint(&bccc, 4, None, None, None, false).expect("BCCC fp");
    assert_eq!(
        bccc_fingerprint.nonzero_elements(),
        &[(4_303_372_288, 1)].into_iter().collect()
    );

    let bcccc_fingerprint =
        getTopologicalTorsionFingerprint(&bcccc, 4, None, None, None, false).expect("BCCCC fp");
    assert_eq!(
        bcccc_fingerprint.nonzero_elements(),
        &[(4_437_590_016, 1), (4_437_590_048, 1)]
            .into_iter()
            .collect()
    );
}

#[test]
#[allow(deprecated)]
fn unsanitized_inputs_require_the_source_explicit_valence_cache() {
    let molecule = Molecule::from_smiles_with_sanitize("CCCC", false).expect("raw molecule");
    let generator =
        getTopologicalTorsionGenerator(&TopologicalTorsionArguments::default(), None, true)
            .expect("generator");

    assert!(matches!(
        generator.getSparseCountFingerprint(&molecule, &mut FingerprintFuncArguments::default()),
        Err(super::FingerprintError::Valence(
            crate::ValenceError::ExplicitValenceCacheNotInitialized { .. }
        ))
    ));
    assert!(matches!(
        getTopologicalTorsionFingerprint(&molecule, 4, None, None, None, false),
        Err(super::FingerprintError::Valence(
            crate::ValenceError::ExplicitValenceCacheNotInitialized { .. }
        ))
    ));

    let cached = molecule
        .with_assigned_valence_strict(false)
        .expect("property cache");
    assert!(
        generator
            .getSparseCountFingerprint(&cached, &mut FingerprintFuncArguments::default())
            .is_ok()
    );
}

fn carbon_chain_query_molecule() -> Molecule {
    let mut builder = Molecule::builder();
    let atoms = (0..4)
        .map(|_| {
            builder.add_atom(
                AtomSpec::new(Element::C)
                    .with_query(QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6))),
            )
        })
        .collect::<Vec<_>>();
    for pair in atoms.windows(2) {
        builder
            .add_bond(
                BondSpec::new(pair[0], pair[1], BondOrder::Single).with_query(
                    QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
                ),
            )
            .expect("query bond");
    }
    builder
        .build()
        .expect("query molecule")
        .with_assigned_valence_strict(false)
        .expect("query property cache")
}

#[test]
#[allow(deprecated)]
fn query_dative_hypervalent_and_explicit_h_graphs_match_full_rdkit_apis() {
    let query = carbon_chain_query_molecule();
    let dative = Molecule::from_smiles("N->[Cu]<-N-C").expect("dative molecule");
    let hypervalent = Molecule::from_smiles("FS(F)(F)(F)(F)F").expect("hypervalent molecule");
    let explicit_h = Molecule::from_smiles("CC")
        .expect("ethane")
        .with_hydrogens()
        .expect("explicit-H graph with source-defined valence state");

    for (name, molecule, expected) in [
        ("query", &query, &[(4_303_372_320, 1)][..]),
        ("dative", &dative, &[(8_715_796_512, 1)][..]),
        ("hypervalent", &hypervalent, &[][..]),
        ("explicit_h", &explicit_h, &[][..]),
    ] {
        let arguments = TopologicalTorsionArguments::new(false, 4, false, vec![1, 2, 4, 8], 2048)
            .expect("arguments");
        let generator = getTopologicalTorsionGenerator(&arguments, None, true).expect("generator");
        let modern = generator
            .getSparseCountFingerprint(molecule, &mut FingerprintFuncArguments::default())
            .unwrap_or_else(|error| panic!("{name} modern fingerprint failed: {error}"));
        let legacy = getTopologicalTorsionFingerprint(molecule, 4, None, None, None, false)
            .unwrap_or_else(|error| panic!("{name} legacy fingerprint failed: {error}"));
        let expected = expected.iter().copied().collect();

        assert_eq!(modern.size(), 1_u64 << 36, "{name} modern size");
        assert_eq!(modern.nonzero_elements(), &expected, "{name} modern ids");
        assert_eq!(legacy.nonzero_elements(), &expected, "{name} legacy ids");
    }
}

#[test]
fn generator_json_roundtrip_preserves_torsion_configuration_and_output() {
    let molecule = Molecule::from_smiles("C1CC1").expect("triangle");
    let mut arguments =
        TopologicalTorsionArguments::new(true, 3, true, vec![1, 3, 5], 1536).expect("arguments");
    arguments.df_only_shortest_paths = true;
    arguments.fingerprint_arguments.d_num_bits_per_feature = 2;
    let generator = getTopologicalTorsionGenerator(&arguments, None, true).expect("generator");
    let typed = TypedFingerprintGenerator::TopologicalTorsion(generator);

    let json = generatorToJSON(&typed);
    let restored = generatorFromJSON(&json).expect("round trip");
    let TypedFingerprintGenerator::TopologicalTorsion(restored_generator) = restored else {
        panic!("expected topological torsion generator")
    };
    assert_eq!(restored_generator.fingerprint_arguments, arguments);
    assert_eq!(
        restored_generator.atom_invariants_generator,
        AtomPairAtomInvGenerator::new(true, true)
    );
    assert!(restored_generator.owns_atom_invariants_generator);

    let expected = typed
        .getFingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .expect("expected");
    let actual = restored_generator
        .getFingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .expect("actual");
    assert_eq!(actual, expected);
}

#[test]
fn supplied_atom_invariant_generator_and_ownership_are_retained() {
    let supplied = AtomPairAtomInvGenerator::new(false, false);
    let generator = getTopologicalTorsionGenerator(
        &TopologicalTorsionArguments::default(),
        Some(supplied.clone()),
        false,
    )
    .expect("generator");
    assert_eq!(generator.atom_invariants_generator, supplied);
    assert!(!generator.owns_atom_invariants_generator);

    let molecule = Molecule::from_smiles("CCCCO").expect("molecule");
    let supplied_output = generator
        .getSparseCountFingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .expect("supplied output");
    let default_output =
        getTopologicalTorsionGenerator(&TopologicalTorsionArguments::default(), None, false)
            .expect("default generator")
            .getSparseCountFingerprint(&molecule, &mut FingerprintFuncArguments::default())
            .expect("default output");
    assert_ne!(supplied_output, default_output);
}

#[test]
fn count_bounds_expand_colliding_counts_at_exact_thresholds() {
    let molecule = Molecule::from_smiles("CC(C)C").expect("branched molecule");
    let arguments =
        TopologicalTorsionArguments::new(false, 3, true, vec![1, 2, 4, 8], 64).expect("arguments");
    let generator = getTopologicalTorsionGenerator(&arguments, None, true).expect("generator");
    let mut call = FingerprintFuncArguments {
        custom_atom_invariants: Some(vec![7; molecule.num_atoms()]),
        ..FingerprintFuncArguments::default()
    };
    let count = generator
        .getCountFingerprint(&molecule, &mut call)
        .expect("count");
    assert_eq!(count.nonzero_elements().len(), 1);
    assert_eq!(sparse_total(&count), 3);

    let mut call = FingerprintFuncArguments {
        custom_atom_invariants: Some(vec![7; molecule.num_atoms()]),
        ..FingerprintFuncArguments::default()
    };
    let bits = generator
        .getFingerprint(&molecule, &mut call)
        .expect("bits");
    assert_eq!(bits.on_bits().len(), 2);
    assert_eq!(bits.on_bits()[1], bits.on_bits()[0] + 1);
}

#[test]
fn extra_bits_are_deterministic_and_use_the_shared_rng_driver() {
    let molecule = Molecule::from_smiles("CCCCO").expect("molecule");
    let mut arguments = TopologicalTorsionArguments::new(false, 4, false, vec![1, 2, 4, 8], 1000)
        .expect("arguments");
    arguments.fingerprint_arguments.d_num_bits_per_feature = 3;
    let generator = getTopologicalTorsionGenerator(&arguments, None, true).expect("generator");

    let first = generator
        .getCountFingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .expect("first");
    let second = generator
        .getCountFingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .expect("second");
    assert_eq!(first, second);
    assert_eq!(sparse_total(&first), 6);
    assert!(first.nonzero_elements().len() >= 2);
}

#[test]
fn chirality_distinguishes_enantiomers_only_when_enabled() {
    let clockwise = Molecule::from_smiles("CC[C@H](F)Cl").expect("clockwise");
    let anticlockwise = Molecule::from_smiles("CC[C@@H](F)Cl").expect("anticlockwise");

    let achiral =
        getTopologicalTorsionGenerator(&TopologicalTorsionArguments::default(), None, true)
            .expect("achiral");
    let clockwise_achiral = achiral
        .getSparseCountFingerprint(&clockwise, &mut FingerprintFuncArguments::default())
        .expect("clockwise achiral");
    let anticlockwise_achiral = achiral
        .getSparseCountFingerprint(&anticlockwise, &mut FingerprintFuncArguments::default())
        .expect("anticlockwise achiral");
    assert_eq!(clockwise_achiral, anticlockwise_achiral);

    let chiral_arguments = TopologicalTorsionArguments::new(true, 4, true, vec![1, 2, 4, 8], 4096)
        .expect("chiral arguments");
    let chiral = getTopologicalTorsionGenerator(&chiral_arguments, None, true).expect("chiral");
    let clockwise_chiral = chiral
        .getSparseCountFingerprint(&clockwise, &mut FingerprintFuncArguments::default())
        .expect("clockwise chiral");
    let anticlockwise_chiral = chiral
        .getSparseCountFingerprint(&anticlockwise, &mut FingerprintFuncArguments::default())
        .expect("anticlockwise chiral");
    assert_ne!(clockwise_chiral, anticlockwise_chiral);
}

#[test]
fn shortest_path_option_prunes_non_shortest_ring_torsions() {
    let molecule = Molecule::from_smiles("C1CC1").expect("triangle");
    let mut arguments = TopologicalTorsionArguments::new(false, 3, false, vec![1, 2, 4, 8], 2048)
        .expect("arguments");
    let all_generator =
        getTopologicalTorsionGenerator(&arguments, None, true).expect("all generator");
    let all = all_generator
        .getSparseCountFingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .expect("all");
    assert_eq!(sparse_total(&all), 3);

    arguments.df_only_shortest_paths = true;
    let shortest_generator =
        getTopologicalTorsionGenerator(&arguments, None, true).expect("shortest generator");
    let shortest = shortest_generator
        .getSparseCountFingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .expect("shortest");
    assert!(shortest.nonzero_elements().is_empty());
}

#[test]
fn every_scalar_form_is_deterministic_across_repeated_calls() {
    let molecule = Molecule::from_smiles("CC(C)CCO").expect("molecule");
    let generator =
        getTopologicalTorsionGenerator(&TopologicalTorsionArguments::default(), None, true)
            .expect("generator");

    assert_eq!(
        generator
            .getSparseCountFingerprint(&molecule, &mut FingerprintFuncArguments::default())
            .unwrap(),
        generator
            .getSparseCountFingerprint(&molecule, &mut FingerprintFuncArguments::default())
            .unwrap()
    );
    assert_eq!(
        generator
            .getSparseFingerprint(&molecule, &mut FingerprintFuncArguments::default())
            .unwrap(),
        generator
            .getSparseFingerprint(&molecule, &mut FingerprintFuncArguments::default())
            .unwrap()
    );
    assert_eq!(
        generator
            .getCountFingerprint(&molecule, &mut FingerprintFuncArguments::default())
            .unwrap(),
        generator
            .getCountFingerprint(&molecule, &mut FingerprintFuncArguments::default())
            .unwrap()
    );
    assert_eq!(
        generator
            .getFingerprint(&molecule, &mut FingerprintFuncArguments::default())
            .unwrap(),
        generator
            .getFingerprint(&molecule, &mut FingerprintFuncArguments::default())
            .unwrap()
    );
}
