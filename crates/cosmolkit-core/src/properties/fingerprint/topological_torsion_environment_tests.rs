use super::{
    AdditionalOutput, AtomPairAtomInvGenerator, FingerprintError, TopologicalTorsionArguments,
    TopologicalTorsionAtomEnv, TopologicalTorsionEnvGenerator, get_topological_torsion_code,
    get_topological_torsion_hash,
};
use crate::Molecule;

fn default_invariants(molecule: &Molecule) -> Vec<u32> {
    AtomPairAtomInvGenerator::new(false, true)
        .getAtomInvariants(molecule)
        .expect("topological torsion atom invariants")
}

fn environments(
    molecule: &Molecule,
    arguments: &TopologicalTorsionArguments,
    from_atoms: Option<&[usize]>,
    ignore_atoms: Option<&[usize]>,
    atom_invariants: &[u32],
    hash_results: bool,
) -> Vec<TopologicalTorsionAtomEnv> {
    TopologicalTorsionEnvGenerator::new()
        .getEnvironments(
            molecule,
            arguments,
            from_atoms,
            ignore_atoms,
            atom_invariants,
            hash_results,
        )
        .expect("valid environments")
}

#[test]
fn null_and_empty_from_atom_selections_are_distinct() {
    let molecule = Molecule::from_smiles("CCCCC").expect("chain");
    let arguments = TopologicalTorsionArguments::default();
    let invariants = default_invariants(&molecule);

    let all = environments(&molecule, &arguments, None, None, &invariants, false);
    let empty = environments(&molecule, &arguments, Some(&[]), None, &invariants, false);

    assert_eq!(all.len(), 2);
    assert!(empty.is_empty());
    assert_eq!(all[0].atom_path, vec![0, 1, 2, 3]);
    assert_eq!(all[1].atom_path, vec![1, 2, 3, 4]);
}

#[test]
fn from_atoms_match_only_path_endpoints() {
    let molecule = Molecule::from_smiles("CCCCC").expect("chain");
    let arguments = TopologicalTorsionArguments::default();
    let invariants = default_invariants(&molecule);

    let internal_only = environments(&molecule, &arguments, Some(&[2]), None, &invariants, false);
    let endpoint_and_internal =
        environments(&molecule, &arguments, Some(&[3]), None, &invariants, false);

    assert!(internal_only.is_empty());
    assert_eq!(endpoint_and_internal.len(), 1);
    assert_eq!(endpoint_and_internal[0].atom_path, vec![0, 1, 2, 3]);
}

#[test]
fn ignored_atoms_remove_paths_when_the_atom_is_internal_or_terminal() {
    let molecule = Molecule::from_smiles("CCCCC").expect("chain");
    let arguments = TopologicalTorsionArguments::default();
    let invariants = default_invariants(&molecule);

    assert!(environments(&molecule, &arguments, None, Some(&[2]), &invariants, false,).is_empty());
    let terminal = environments(&molecule, &arguments, None, Some(&[0]), &invariants, false);
    assert_eq!(terminal.len(), 1);
    assert_eq!(terminal[0].atom_path, vec![1, 2, 3, 4]);
}

#[test]
fn custom_invariants_receive_endpoint_and_internal_source_corrections() {
    let molecule = Molecule::from_smiles("CCCCC").expect("chain");
    let arguments = TopologicalTorsionArguments::default();
    let invariants = [10, 20, 30, 40, 50];
    let generated = environments(&molecule, &arguments, None, None, &invariants, false);

    assert_eq!(
        generated
            .iter()
            .map(|env| env.getBitId())
            .collect::<Vec<_>>(),
        vec![
            get_topological_torsion_code(&[11, 20, 30, 41], false).unwrap(),
            get_topological_torsion_code(&[21, 30, 40, 51], false).unwrap(),
        ]
    );
}

#[test]
fn legacy_unfolded_atom_codes_bypass_the_modern_base_code_modulus() {
    let molecule = Molecule::from_smiles("CCCC").expect("chain");
    let mut arguments = TopologicalTorsionArguments::default();
    arguments.fingerprint_arguments.df_include_chirality = true;
    let invariants = [600, 600, 600, 600];

    let modern = TopologicalTorsionEnvGenerator::new()
        .getEnvironments(&molecule, &arguments, None, None, &invariants, false)
        .expect("modern environments");
    let mut legacy_generator = TopologicalTorsionEnvGenerator::new();
    legacy_generator.use_legacy_unfolded_atom_codes();
    let legacy = legacy_generator
        .getEnvironments(&molecule, &arguments, None, None, &invariants, false)
        .expect("legacy environments");

    assert_eq!(modern.len(), 1);
    assert_eq!(legacy.len(), 1);
    assert_eq!(
        modern[0].getBitId(),
        get_topological_torsion_code(&[90, 89, 89, 90], true).unwrap()
    );
    assert_eq!(
        legacy[0].getBitId(),
        get_topological_torsion_code(&[601, 600, 600, 601], true).unwrap()
    );
    assert_ne!(modern[0].getBitId(), legacy[0].getBitId());
}

#[test]
fn full_ring_closure_is_kept_but_a_later_illegal_repeat_is_not_enumerated() {
    let molecule = Molecule::from_smiles("C1CCC1").expect("four-membered ring");
    let invariants = default_invariants(&molecule);
    let mut arguments = TopologicalTorsionArguments::default();
    arguments.d_torsion_atom_count = 5;

    let closed = environments(&molecule, &arguments, None, None, &invariants, false);
    assert_eq!(closed.len(), 1);
    assert_eq!(closed[0].atom_path, vec![0, 1, 2, 3, 0]);

    arguments.d_torsion_atom_count = 6;
    assert!(environments(&molecule, &arguments, None, None, &invariants, false).is_empty());
}

#[test]
fn hashed_and_unhashed_ids_preserve_source_collisions() {
    let molecule = Molecule::from_smiles("CC(C)C").expect("branched molecule");
    let mut arguments = TopologicalTorsionArguments::default();
    arguments.d_torsion_atom_count = 3;
    let invariants = vec![7; molecule.num_atoms()];

    let packed = environments(&molecule, &arguments, None, None, &invariants, false);
    let hashed = environments(&molecule, &arguments, None, None, &invariants, true);
    assert_eq!(packed.len(), 3);
    assert_eq!(hashed.len(), 3);
    assert!(
        packed
            .iter()
            .all(|env| env.getBitId() == packed[0].getBitId())
    );
    assert!(
        hashed
            .iter()
            .all(|env| env.getBitId() == hashed[0].getBitId())
    );
    assert_eq!(
        packed[0].getBitId(),
        get_topological_torsion_code(&[8, 7, 8], false).unwrap()
    );
    assert_eq!(
        hashed[0].getBitId(),
        u64::from(get_topological_torsion_hash(&[8, 7, 8]).unwrap())
    );
}

#[test]
fn atom_environment_updates_every_supported_provenance_allocation() {
    let environment = TopologicalTorsionAtomEnv::new(99, vec![0, 2, 4]);
    let mut output = AdditionalOutput::new();
    output.allocate_atom_to_bits();
    output.allocate_atom_counts();
    output.allocate_bit_paths();
    output.allocate_atoms_per_bit();
    output.reset_for_atom_count(5);

    environment.updateAdditionalOutput(&mut output, 7);
    assert_eq!(environment.getBitId(), 99);
    assert_eq!(
        output.atom_to_bits,
        Some(vec![vec![7], vec![], vec![7], vec![], vec![7]])
    );
    assert_eq!(output.atom_counts, Some(vec![1, 0, 1, 0, 1]));
    assert_eq!(
        output.bit_paths.as_ref().unwrap().get(&7),
        Some(&vec![vec![0, 2, 4]])
    );
    assert_eq!(
        output.atoms_per_bit.as_ref().unwrap().get(&7),
        Some(&vec![vec![0, 2, 4]])
    );
}

#[test]
fn provenance_keeps_duplicate_colliding_paths_and_repeated_atoms() {
    let mut output = AdditionalOutput::new();
    output.allocate_atom_to_bits();
    output.allocate_atom_counts();
    output.allocate_bit_paths();
    output.allocate_atoms_per_bit();
    output.reset_for_atom_count(3);

    TopologicalTorsionAtomEnv::new(4, vec![0, 1, 0]).updateAdditionalOutput(&mut output, 5);
    TopologicalTorsionAtomEnv::new(4, vec![1, 2]).updateAdditionalOutput(&mut output, 5);

    assert_eq!(
        output.atom_to_bits,
        Some(vec![vec![5, 5], vec![5, 5], vec![5]])
    );
    assert_eq!(output.atom_counts, Some(vec![2, 2, 1]));
    assert_eq!(
        output.bit_paths.as_ref().unwrap().get(&5),
        Some(&vec![vec![0, 1, 0], vec![1, 2]])
    );
    assert_eq!(output.atoms_per_bit, output.bit_paths);
}

#[test]
fn environment_metadata_and_invalid_inputs_are_structured() {
    let mut generator = TopologicalTorsionEnvGenerator::new();
    assert_eq!(generator.infoString(), "TopologicalTorsionEnvGenerator");
    assert_eq!(
        serde_json::from_str::<serde_json::Value>(&generator.toJSON()).unwrap(),
        serde_json::json!({"type": "TopologicalTorsionEnvGenerator"})
    );
    generator.fromJSON("").expect("empty property tree");
    generator
        .fromJSON(r#"{"type":"DifferentGenerator","extra":1}"#)
        .expect("the stateless source unit ignores fields");
    assert!(matches!(
        generator.fromJSON("[]"),
        Err(FingerprintError::InvalidArgumentsJson(_))
    ));

    let molecule = Molecule::from_smiles("CCCC").expect("chain");
    let arguments = TopologicalTorsionArguments::default();
    assert!(matches!(
        generator.getEnvironments(&molecule, &arguments, None, None, &[1, 2], false),
        Err(FingerprintError::InvalidArguments { .. })
    ));
    assert!(matches!(
        generator.getEnvironments(
            &molecule,
            &arguments,
            Some(&[molecule.num_atoms()]),
            None,
            &[1, 2, 3, 4],
            false,
        ),
        Err(FingerprintError::InvalidArguments { .. })
    ));
}
