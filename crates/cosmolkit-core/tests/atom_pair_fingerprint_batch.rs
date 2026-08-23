use std::sync::{
    Arc,
    atomic::{AtomicUsize, Ordering},
};

use cosmolkit_core::{
    AtomPairFingerprintGenerator, AtomPairFingerprintParams, BatchRecord, FingerprintFuncArguments,
    Molecule, MoleculeBatch, MorganFingerprintParams,
};

fn smiles(values: &[&str]) -> Vec<String> {
    values.iter().map(|value| (*value).to_string()).collect()
}

fn arguments(params: &AtomPairFingerprintParams) -> FingerprintFuncArguments {
    FingerprintFuncArguments {
        from_atoms: params.from_atoms.clone(),
        ignore_atoms: params.ignore_atoms.clone(),
        conf_id: params.conformer_id,
        custom_atom_invariants: params.custom_atom_invariants.clone(),
        ..Default::default()
    }
}

#[test]
fn every_batch_result_form_matches_ordered_scalar_calls() {
    let source = smiles(&["C", "CCO", "c1ccccc1", "C[C@H](O)F"]);
    let molecules: Vec<_> = source
        .iter()
        .map(|value| Molecule::from_smiles(value).unwrap())
        .collect();
    let batch = MoleculeBatch::from_smiles_list(&source);
    let params = AtomPairFingerprintParams {
        n_bits: 256,
        count_simulation: true,
        count_bounds: vec![1, 3, 5],
        num_bits_per_feature: 2,
        use_chirality: true,
        ..Default::default()
    };
    let generator = AtomPairFingerprintGenerator::new(&params).unwrap();

    let sparse_count = batch
        .atom_pair_sparse_count_fingerprint_list_with_options(&params, Some(4), Some(false))
        .unwrap();
    let expected_sparse_count: Vec<_> = molecules
        .iter()
        .map(|molecule| {
            Some(
                generator
                    .sparse_count_fingerprint(molecule, &mut arguments(&params))
                    .unwrap(),
            )
        })
        .collect();
    assert_eq!(sparse_count, expected_sparse_count);

    let count = batch
        .atom_pair_count_fingerprint_list_with_options(&params, Some(3), Some(false))
        .unwrap();
    let expected_count: Vec<_> = molecules
        .iter()
        .map(|molecule| {
            Some(
                generator
                    .count_fingerprint(molecule, &mut arguments(&params))
                    .unwrap(),
            )
        })
        .collect();
    assert_eq!(count, expected_count);

    let sparse_bits = batch
        .atom_pair_sparse_bit_fingerprint_list_with_options(&params, Some(2), Some(false))
        .unwrap();
    let expected_sparse_bits: Vec<_> = molecules
        .iter()
        .map(|molecule| {
            Some(
                generator
                    .sparse_bit_fingerprint(molecule, &mut arguments(&params))
                    .unwrap(),
            )
        })
        .collect();
    assert_eq!(sparse_bits, expected_sparse_bits);

    let fingerprints = batch
        .atom_pair_fingerprint_list_with_options(&params, Some(4), Some(false))
        .unwrap();
    let expected_fingerprints: Vec<_> = molecules
        .iter()
        .map(|molecule| Some(molecule.atom_pair_fingerprint(&params).unwrap()))
        .collect();
    assert_eq!(fingerprints, expected_fingerprints);

    let output_params = AtomPairFingerprintParams {
        collect_additional_output: true,
        ..params
    };
    let outputs = batch
        .atom_pair_fingerprint_with_output_list_with_options(&output_params, Some(4), Some(false))
        .unwrap();
    let expected_outputs: Vec<_> = molecules
        .iter()
        .map(|molecule| {
            Some(
                molecule
                    .atom_pair_fingerprint_with_output(&output_params)
                    .unwrap(),
            )
        })
        .collect();
    assert_eq!(outputs, expected_outputs);
}

#[test]
fn batch_order_thread_count_progress_and_repeated_calls_are_deterministic() {
    let source = smiles(&["C", "CC", "CCC", "CCCC", "CCCCC", "CCCCCC"]);
    let batch = MoleculeBatch::from_smiles_list(&source).with_parallel_jobs(Some(4));
    let params = AtomPairFingerprintParams::default();
    let baseline = batch
        .atom_pair_fingerprint_list_with_options(&params, Some(1), Some(false))
        .unwrap();
    for n_jobs in [1, 2, 4] {
        for _ in 0..3 {
            assert_eq!(
                batch
                    .atom_pair_fingerprint_list_with_options(&params, Some(n_jobs), Some(false))
                    .unwrap(),
                baseline
            );
        }
    }

    let ticks = AtomicUsize::new(0);
    let progress = || {
        ticks.fetch_add(1, Ordering::Relaxed);
    };
    assert_eq!(
        batch
            .atom_pair_fingerprint_list_with_progress(&params, Some(&progress))
            .unwrap(),
        baseline
    );
    assert_eq!(ticks.load(Ordering::Relaxed), source.len());
}

#[test]
fn invalid_input_records_keep_their_indices_and_operation_errors_are_indexed() {
    let source = smiles(&["CC", "C1", "CCC"]);
    let batch = MoleculeBatch::from_smiles_list(&source);
    assert_eq!(batch.errors().len(), 1);
    assert_eq!(batch.errors()[0].index, 1);

    let values = batch
        .atom_pair_fingerprint_list_with_options(
            &AtomPairFingerprintParams::default(),
            Some(2),
            Some(false),
        )
        .unwrap();
    assert!(values[0].is_some());
    assert!(values[1].is_none());
    assert!(values[2].is_some());

    let no_conformer = AtomPairFingerprintParams {
        use_2d: false,
        ..Default::default()
    };
    let error = batch
        .atom_pair_fingerprint_list_with_options(&no_conformer, Some(2), Some(false))
        .unwrap_err();
    assert_eq!(
        error
            .record_errors
            .iter()
            .map(|error| error.index)
            .collect::<Vec<_>>(),
        vec![0, 2]
    );
    assert!(
        error
            .record_errors
            .iter()
            .all(|error| error.operation == "batch.atom_pair_fingerprint")
    );

    let bad_config = AtomPairFingerprintParams {
        n_bits: 0,
        ..Default::default()
    };
    let error = batch
        .atom_pair_fingerprint_list_with_options(&bad_config, Some(2), Some(false))
        .unwrap_err();
    assert_eq!(error.record_errors[0].index, 0);
    assert_eq!(
        error.record_errors[0].operation,
        "batch.atom_pair_fingerprint"
    );
}

#[test]
fn shared_batch_is_immutable_and_safe_across_mixed_family_threads() {
    let source = smiles(&["CCO", "c1ccccc1O", "C[C@H](O)F", "CC(C)C(=O)O"]);
    let batch = Arc::new(MoleculeBatch::from_smiles_list(&source));
    let snapshot = (*batch).clone();
    let atom_pair_params = AtomPairFingerprintParams::default();
    let expected = batch
        .atom_pair_fingerprint_list_with_options(&atom_pair_params, Some(2), Some(false))
        .unwrap();

    let atom_pair_batch = Arc::clone(&batch);
    let atom_pair = std::thread::spawn(move || {
        atom_pair_batch
            .atom_pair_fingerprint_list_with_options(
                &AtomPairFingerprintParams::default(),
                Some(4),
                Some(false),
            )
            .unwrap()
    });
    let morgan_batch = Arc::clone(&batch);
    let morgan = std::thread::spawn(move || {
        morgan_batch
            .morgan_fingerprint_list_with_options(
                &MorganFingerprintParams::default(),
                Some(3),
                Some(false),
            )
            .unwrap()
    });

    assert_eq!(atom_pair.join().unwrap(), expected);
    assert_eq!(morgan.join().unwrap().len(), source.len());
    assert_eq!(*batch, snapshot);
    assert!(
        batch
            .iter()
            .all(|record| matches!(record, BatchRecord::Molecule(_)))
    );
}
