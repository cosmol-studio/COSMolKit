//! Original d892 public cases, canonical rc.10 facade transport only.
#![cfg(all(
    feature = "cap-stereoisomers",
    feature = "cap-smiles",
    feature = "cap-io",
    feature = "cap-fingerprints"
))]
use cosmolkit::{
    CoordinateBlock, Molecule, MoleculeBuilder, MorganFingerprintParams, OperationError,
    StereoisomerOptions,
};
use num_bigint::{BigInt, BigUint};
fn mol_to_binary(molecule: &Molecule) -> Result<Vec<u8>, cosmolkit::PickleError> {
    molecule.to_binary()
}
fn stereoisomer_count(
    molecule: &Molecule,
    options: &StereoisomerOptions,
) -> Result<BigUint, cosmolkit::EnumerationError> {
    molecule.stereoisomer_count_with_options(options)
}
fn enumerate_stereoisomers(
    molecule: &Molecule,
    options: StereoisomerOptions,
) -> Result<cosmolkit::StereoisomerIterator, OperationError> {
    molecule.enumerate_stereoisomers_with_options(&options)
}

fn canonical_isomeric_smiles(molecule: &Molecule) -> String {
    std::str::from_utf8(
        molecule
            .to_smiles()
            .expect("canonical isomeric SMILES")
            .as_bytes(),
    )
    .expect("ASCII fixture SMILES")
    .to_owned()
}

fn collect_smiles(
    molecule: &Molecule,
    options: StereoisomerOptions,
) -> Result<Vec<String>, OperationError> {
    enumerate_stereoisomers(molecule, options)?
        .map(|result| result.map(|isomer| canonical_isomeric_smiles(&isomer)))
        .collect()
}

#[test]
fn stereoisomer_options_match_the_pinned_python_defaults() {
    let options = StereoisomerOptions::default();
    assert!(!options.try_embedding);
    assert!(options.only_unassigned);
    assert!(!options.only_stereo_groups);
    assert_eq!(options.max_isomers, 1024);
    assert_eq!(options.random_source, None);
    assert!(options.unique);
}

#[test]
fn no_center_public_methods_yield_one_clone_and_preserve_the_source() {
    let base = Molecule::from_smiles("CC").expect("ethane");
    let source = MoleculeBuilder::from_parts(
        base.topology().clone(),
        CoordinateBlock::default(),
        base.properties().clone(),
    )
    .with_name("source".into())
    .with_property("record".into(), "kept".into())
    .unwrap()
    .build()
    .unwrap();
    let before = mol_to_binary(&source).expect("source binary snapshot");
    let options = StereoisomerOptions::default();

    assert_eq!(
        source.stereoisomer_count_with_options(&options).unwrap(),
        BigUint::from(1_u8)
    );
    assert_eq!(
        stereoisomer_count(&source, &options).unwrap(),
        BigUint::from(1_u8)
    );

    let mut iterator = source
        .enumerate_stereoisomers_with_options(&options)
        .unwrap();
    assert_eq!(iterator.yielded_count(), 0);
    let only = iterator.next().expect("one result").expect("valid result");
    assert_eq!(canonical_isomeric_smiles(&only), "CC");
    assert_eq!(iterator.yielded_count(), 1);
    assert!(iterator.next().is_none());
    assert!(
        iterator.next().is_none(),
        "the public iterator must be fused"
    );
    assert_eq!(
        mol_to_binary(&source).expect("source after enumeration"),
        before,
        "enumeration changed the source molecule"
    );
}

#[test]
fn atom_and_double_bond_outputs_match_pinned_rdkit_order_exactly() {
    let source = Molecule::from_smiles("CC(F)=CC(Cl)C").expect("mixed stereo fixture");
    let before = mol_to_binary(&source).expect("source binary snapshot");
    let expected = vec![
        r"C/C(F)=C\[C@@H](C)Cl".to_string(),
        r"C/C(F)=C\[C@H](C)Cl".to_string(),
        r"C/C(F)=C/[C@@H](C)Cl".to_string(),
        r"C/C(F)=C/[C@H](C)Cl".to_string(),
    ];

    assert_eq!(
        source
            .stereoisomer_count_with_options(&StereoisomerOptions::default())
            .unwrap(),
        BigUint::from(4_u8)
    );
    assert_eq!(
        collect_smiles(&source, StereoisomerOptions::default()).unwrap(),
        expected
    );
    assert_eq!(mol_to_binary(&source).unwrap(), before);
}

#[test]
fn option_matrix_matches_assignment_and_uniqueness_branches() {
    let assigned = Molecule::from_smiles("C/C(F)=C/[C@@H](C)Cl").expect("assigned fixture");
    assert_eq!(
        collect_smiles(&assigned, StereoisomerOptions::default()).unwrap(),
        vec!["C/C(F)=C/[C@@H](C)Cl"]
    );

    let all_assigned = StereoisomerOptions {
        only_unassigned: false,
        ..Default::default()
    };
    assert_eq!(
        collect_smiles(&assigned, all_assigned).unwrap(),
        vec![
            r"C/C(F)=C\[C@H](C)Cl",
            r"C/C(F)=C\[C@@H](C)Cl",
            r"C/C(F)=C/[C@H](C)Cl",
            r"C/C(F)=C/[C@@H](C)Cl",
        ]
    );

    let meso = Molecule::from_smiles("FC(Cl)C(Cl)F").expect("meso fixture");
    assert_eq!(
        collect_smiles(&meso, StereoisomerOptions::default()).unwrap(),
        vec![
            "F[C@H](Cl)[C@@H](F)Cl",
            "F[C@@H](Cl)[C@@H](F)Cl",
            "F[C@H](Cl)[C@H](F)Cl",
        ]
    );
    let non_unique = StereoisomerOptions {
        unique: false,
        ..Default::default()
    };
    assert_eq!(
        collect_smiles(&meso, non_unique).unwrap(),
        vec![
            "F[C@H](Cl)[C@@H](F)Cl",
            "F[C@@H](Cl)[C@@H](F)Cl",
            "F[C@H](Cl)[C@H](F)Cl",
            "F[C@H](Cl)[C@@H](F)Cl",
        ]
    );
}

#[test]
fn random_iteration_is_lazy_and_supports_more_than_machine_word_centers() {
    let smiles = format!("Br{}F", "[CH](Cl)".repeat(101));
    let source = Molecule::from_smiles(&smiles).expect("101-center fixture");
    let options = StereoisomerOptions {
        max_isomers: 3,
        random_source: Some(cosmolkit::StereoisomerRandomSource::from_integer_seed(
            BigInt::from(61_453),
        )),
        ..Default::default()
    };

    assert_eq!(
        source.stereoisomer_count_with_options(&options).unwrap(),
        BigUint::from(1_u8) << 101_usize
    );
    let mut iterator = source
        .enumerate_stereoisomers_with_options(&options)
        .unwrap();
    assert_eq!(iterator.yielded_count(), 0);
    for expected_count in 1..=3 {
        let isomer = iterator
            .next()
            .expect("bounded random result")
            .expect("valid random result");
        assert_eq!(isomer.num_atoms(), source.num_atoms());
        assert_eq!(iterator.yielded_count(), expected_count);
    }
    assert!(iterator.next().is_none());
}

#[test]
fn repeated_composed_and_parallel_calls_match_scalar_results() {
    let inputs = ["CC(F)=CC(Cl)C", "FC(Cl)C(Cl)F", "CC1CC(C)C1", "CCC(C)Br"];
    let molecules = inputs
        .iter()
        .map(|smiles| Molecule::from_smiles(smiles).expect("fixture molecule"))
        .collect::<Vec<_>>();
    let scalar = molecules
        .iter()
        .map(|molecule| collect_smiles(molecule, StereoisomerOptions::default()).unwrap())
        .collect::<Vec<_>>();

    for _ in 0..3 {
        let repeated = molecules
            .iter()
            .map(|molecule| collect_smiles(molecule, StereoisomerOptions::default()).unwrap())
            .collect::<Vec<_>>();
        assert_eq!(repeated, scalar);
    }

    let parallel = std::thread::scope(|scope| {
        molecules
            .iter()
            .map(|molecule| {
                scope.spawn(move || {
                    collect_smiles(molecule, StereoisomerOptions::default()).unwrap()
                })
            })
            .collect::<Vec<_>>()
            .into_iter()
            .map(|task| task.join().unwrap())
            .collect::<Vec<_>>()
    });
    assert_eq!(parallel, scalar);

    for smiles in scalar.iter().flatten() {
        let roundtrip = Molecule::from_smiles(smiles).expect("enumerated SMILES roundtrip");
        assert_eq!(canonical_isomeric_smiles(&roundtrip), *smiles);
        assert!(
            roundtrip
                .fingerprint_morgan_with_params(&MorganFingerprintParams::default(), None)
                .expect("fingerprint composition")
                .on_bits()
                .len()
                <= MorganFingerprintParams::default().generator.fp_size as usize
        );
    }
}

#[test]
fn default_count_preserves_the_upper_bound_and_never_consumes_random_sources() {
    let source = Molecule::from_smiles("FC(Cl)C(Cl)F").expect("meso fixture");
    let before = source.to_binary().unwrap();
    assert_eq!(source.stereoisomer_count().unwrap(), BigUint::from(4_u8));
    let mut options = StereoisomerOptions {
        max_isomers: 1,
        random_source: Some(cosmolkit::StereoisomerRandomSource::from_random_bits(
            Box::new(|_| panic!("count must not consume the random source")),
        )),
        ..Default::default()
    };
    assert_eq!(
        source.stereoisomer_count_with_options(&options).unwrap(),
        BigUint::from(4_u8)
    );
    options.only_stereo_groups = true;
    assert_eq!(
        source.stereoisomer_count_with_options(&options).unwrap(),
        BigUint::from(1_u8)
    );
    assert_eq!(source.to_binary().unwrap(), before);
}
