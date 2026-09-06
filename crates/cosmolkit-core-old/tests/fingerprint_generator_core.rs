use std::collections::{BTreeMap, BTreeSet};

use cosmolkit_core::fingerprint::{
    AdditionalOutput, FingerprintFuncArguments, getMorganGeneratorWithParams,
};
use cosmolkit_core::{AtomSpec, BondOrder, BondSpec, Element, Molecule};

fn morgan_generator(
    fp_size: u32,
    count_simulation: bool,
) -> cosmolkit_core::fingerprint::MorganFingerprintGenerator {
    getMorganGeneratorWithParams(
        1,
        count_simulation,
        false,
        true,
        false,
        false,
        None,
        None,
        fp_size,
        vec![1, 2, 4, 8],
        false,
        false,
    )
    .expect("Morgan generator")
}

fn all_additional_output() -> AdditionalOutput {
    let mut output = AdditionalOutput::new();
    output.allocate_atom_to_bits();
    output.allocate_bit_info_map();
    output.allocate_bit_paths();
    output.allocate_atom_counts();
    output.allocate_atoms_per_bit();
    output
}

fn molecule_without_stereochem_done() -> Molecule {
    let mut builder = Molecule::builder();
    let left = builder.add_atom(AtomSpec::new(Element::C));
    let right = builder.add_atom(AtomSpec::new(Element::C));
    builder
        .add_bond(BondSpec::new(left, right, BondOrder::Single))
        .expect("C-C bond");
    builder
        .build()
        .expect("molecule")
        .with_assigned_valence()
        .expect("valence-initialized molecule")
}

#[test]
fn shared_generator_core_preserves_all_four_morgan_scalar_outputs() {
    let molecule = Molecule::from_smiles("CC").expect("ethane");
    let generator = morgan_generator(2048, false);

    let sparse_count = generator
        .getSparseCountFingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .expect("sparse count");
    let sparse_bit = generator
        .getSparseFingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .expect("sparse bit");
    let folded_count = generator
        .getCountFingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .expect("folded count");
    let explicit_bit = generator
        .getFingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .expect("explicit bit");

    assert_eq!(sparse_count.size(), u64::from(u32::MAX));
    assert_eq!(
        sparse_count.nonzero_elements(),
        &BTreeMap::from([(2_246_728_737, 2), (3_545_175_291, 1)])
    );
    assert_eq!(sparse_bit.size(), u64::from(u32::MAX));
    assert_eq!(
        sparse_bit.on_bits(),
        &BTreeSet::from([2_246_728_737, 3_545_175_291])
    );
    assert_eq!(
        folded_count.nonzero_elements(),
        &BTreeMap::from([(1057, 2), (1275, 1)])
    );
    assert_eq!(explicit_bit.on_bits(), vec![1057, 1275]);
}

#[test]
fn shared_generator_core_handles_collisions_count_simulation_and_extra_bits() {
    let molecule = Molecule::from_smiles("CC").expect("ethane");

    let collision_generator = morgan_generator(1, false);
    let collision_count = collision_generator
        .getCountFingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .expect("collision count");
    let collision_bits = collision_generator
        .getFingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .expect("collision bits");
    assert_eq!(
        collision_count.nonzero_elements(),
        &BTreeMap::from([(0, 3)])
    );
    assert_eq!(collision_bits.on_bits(), vec![0]);

    let count_generator = morgan_generator(2048, true);
    let mut count_args = FingerprintFuncArguments {
        additional_output: Some(all_additional_output()),
        ..Default::default()
    };
    let simulated = count_generator
        .getFingerprint(&molecule, &mut count_args)
        .expect("count-simulated bit vector");
    assert_eq!(simulated.on_bits(), vec![132, 133, 1004]);

    let mut random_generator = morgan_generator(2048, false);
    random_generator
        .fingerprint_arguments
        .fingerprint_arguments
        .d_num_bits_per_feature = 2;
    let first = random_generator
        .getCountFingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .expect("extra-bit count");
    let second = random_generator
        .getCountFingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .expect("repeated extra-bit count");
    assert_eq!(first, second, "extra-bit RNG must be deterministic");
    assert_eq!(first.nonzero_elements().values().sum::<i32>(), 6);
    assert!(first.nonzero_elements().len() > 2);
}

#[test]
fn shared_generator_core_reinitializes_all_outputs_but_retains_atoms_per_bit() {
    let molecule = Molecule::from_smiles("CC").expect("ethane");
    let generator = morgan_generator(2048, false);
    let mut args = FingerprintFuncArguments {
        additional_output: Some(all_additional_output()),
        ..Default::default()
    };

    generator
        .getSparseCountFingerprint(&molecule, &mut args)
        .expect("first fingerprint");
    let first = args
        .additional_output
        .as_ref()
        .expect("first additional output")
        .clone();
    generator
        .getSparseCountFingerprint(&molecule, &mut args)
        .expect("second fingerprint");
    let second = args
        .additional_output
        .as_ref()
        .expect("second additional output");

    assert_eq!(second.atom_counts, first.atom_counts);
    assert_eq!(second.atom_to_bits, first.atom_to_bits);
    assert_eq!(second.bit_info_map, first.bit_info_map);
    assert_eq!(second.bit_paths, Some(BTreeMap::new()));
    for (bit_id, first_paths) in first.atoms_per_bit.expect("first atomsPerBit") {
        let second_paths = second
            .atoms_per_bit
            .as_ref()
            .expect("second atomsPerBit")
            .get(&bit_id)
            .expect("retained atomsPerBit id");
        assert_eq!(second_paths.len(), first_paths.len() * 2);
        assert_eq!(&second_paths[..first_paths.len()], first_paths.as_slice());
        assert_eq!(&second_paths[first_paths.len()..], first_paths.as_slice());
    }
}

#[test]
fn shared_generator_core_prepares_a_private_stereo_copy_without_mutating_source() {
    let molecule = molecule_without_stereochem_done();
    assert_eq!(molecule.prop("_StereochemDone"), None);
    let source_atoms = molecule.atoms().to_vec();
    let source_bonds = molecule.bonds().to_vec();
    let source_smiles = molecule.to_smiles(true).expect("source SMILES");

    let mut generator = morgan_generator(2048, false);
    generator
        .fingerprint_arguments
        .fingerprint_arguments
        .df_include_chirality = true;
    generator.bond_invariants_generator.use_chirality = true;

    let prepared_result = generator
        .getSparseCountFingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .expect("private-copy stereo preparation");
    let premarked = molecule.clone().with_prop("_StereochemDone", "1");
    let premarked_result = generator
        .getSparseCountFingerprint(&premarked, &mut FingerprintFuncArguments::default())
        .expect("pre-marked stereo molecule");

    assert_eq!(prepared_result, premarked_result);
    assert_eq!(molecule.prop("_StereochemDone"), None);
    assert_eq!(molecule.atoms(), source_atoms.as_slice());
    assert_eq!(molecule.bonds(), source_bonds.as_slice());
    assert_eq!(
        molecule.to_smiles(true).expect("final SMILES"),
        source_smiles
    );
}
