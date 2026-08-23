use cosmolkit_core::{
    BatchRecord, BatchRecordError, Molecule, MoleculeBatch,
    TopologicalTorsionFingerprintOutputRequest, TopologicalTorsionFingerprintParams,
    TopologicalTorsionFingerprintVector, topological_torsion_count_fingerprint,
    topological_torsion_fingerprint, topological_torsion_fingerprint_with_output,
    topological_torsion_sparse_count_fingerprint, topological_torsion_sparse_fingerprint,
};

fn molecule(smiles: &str) -> Molecule {
    Molecule::from_smiles(smiles).unwrap_or_else(|error| panic!("parse {smiles}: {error}"))
}

fn valid_molecules() -> Vec<Molecule> {
    ["CCCC", "CCCCC", "CC(C)CC", "c1ccccc1", "C[C@H](F)Cl"]
        .into_iter()
        .map(molecule)
        .collect()
}

fn batch_from_molecules(molecules: &[Molecule]) -> MoleculeBatch {
    MoleculeBatch::new(
        molecules
            .iter()
            .cloned()
            .map(BatchRecord::Molecule)
            .collect(),
    )
}

#[test]
fn ordered_batch_vectors_match_the_single_scalar_core_for_every_thread_count() {
    let molecules = valid_molecules();
    let batch = batch_from_molecules(&molecules);
    let batch_before = batch.clone();
    let molecules_before = molecules.clone();
    let params = TopologicalTorsionFingerprintParams {
        include_chirality: true,
        num_bits_per_feature: 3,
        ..TopologicalTorsionFingerprintParams::default()
    };
    let output_request = TopologicalTorsionFingerprintOutputRequest {
        vector: TopologicalTorsionFingerprintVector::Count,
        atom_to_bits: true,
        atom_counts: true,
        bit_paths: true,
        atoms_per_bit: true,
    };

    let expected_sparse_counts = molecules
        .iter()
        .map(|molecule| {
            Some(
                topological_torsion_sparse_count_fingerprint(molecule, &params)
                    .expect("scalar sparse-count fingerprint"),
            )
        })
        .collect::<Vec<_>>();
    let expected_sparse_bits = molecules
        .iter()
        .map(|molecule| {
            Some(
                topological_torsion_sparse_fingerprint(molecule, &params)
                    .expect("scalar sparse-bit fingerprint"),
            )
        })
        .collect::<Vec<_>>();
    let expected_counts = molecules
        .iter()
        .map(|molecule| {
            Some(
                topological_torsion_count_fingerprint(molecule, &params)
                    .expect("scalar count fingerprint"),
            )
        })
        .collect::<Vec<_>>();
    let expected_bits = molecules
        .iter()
        .map(|molecule| {
            Some(
                topological_torsion_fingerprint(molecule, &params).expect("scalar bit fingerprint"),
            )
        })
        .collect::<Vec<_>>();
    let expected_with_output = molecules
        .iter()
        .map(|molecule| {
            Some(
                topological_torsion_fingerprint_with_output(molecule, &params, output_request)
                    .expect("scalar fingerprint with output"),
            )
        })
        .collect::<Vec<_>>();

    for n_jobs in [None, Some(1), Some(2), Some(4)] {
        assert_eq!(
            batch
                .topological_torsion_sparse_count_fingerprint_list_with_options(
                    &params,
                    n_jobs,
                    Some(false),
                )
                .expect("batch sparse-count fingerprints"),
            expected_sparse_counts,
            "sparse-count mismatch for n_jobs={n_jobs:?}"
        );
        assert_eq!(
            batch
                .topological_torsion_sparse_fingerprint_list_with_options(
                    &params,
                    n_jobs,
                    Some(false),
                )
                .expect("batch sparse-bit fingerprints"),
            expected_sparse_bits,
            "sparse-bit mismatch for n_jobs={n_jobs:?}"
        );
        assert_eq!(
            batch
                .topological_torsion_count_fingerprint_list_with_options(
                    &params,
                    n_jobs,
                    Some(false),
                )
                .expect("batch count fingerprints"),
            expected_counts,
            "count mismatch for n_jobs={n_jobs:?}"
        );
        assert_eq!(
            batch
                .topological_torsion_fingerprint_list_with_options(&params, n_jobs, Some(false),)
                .expect("batch bit fingerprints"),
            expected_bits,
            "bit mismatch for n_jobs={n_jobs:?}"
        );
        assert_eq!(
            batch
                .topological_torsion_fingerprint_with_output_list_with_options(
                    &params,
                    output_request,
                    n_jobs,
                    Some(false),
                )
                .expect("batch fingerprints with output"),
            expected_with_output,
            "additional-output mismatch for n_jobs={n_jobs:?}"
        );
    }

    let configured = batch.clone().with_parallel_jobs(Some(3));
    assert_eq!(
        configured
            .topological_torsion_fingerprint_list_with_options(&params, None, Some(false))
            .expect("configured parallel batch"),
        expected_bits
    );

    for _ in 0..8 {
        assert_eq!(
            batch
                .topological_torsion_fingerprint_with_output_list_with_options(
                    &params,
                    output_request,
                    Some(4),
                    Some(false),
                )
                .expect("repeated parallel batch"),
            expected_with_output
        );
    }

    assert_eq!(
        batch, batch_before,
        "batch generation mutated the source batch"
    );
    assert_eq!(
        molecules, molecules_before,
        "batch generation mutated a source molecule"
    );
}

#[test]
fn batch_invalid_records_stay_aligned_and_new_errors_keep_their_input_index() {
    let first = molecule("CCCC");
    let last = molecule("CC");
    let original_error = BatchRecordError::new(1, "fixture.parse", "invalid source record");
    let batch = MoleculeBatch::new(vec![
        BatchRecord::Molecule(first.clone()),
        BatchRecord::Error(original_error.clone()),
        BatchRecord::Molecule(last.clone()),
    ]);
    let batch_before = batch.clone();

    let values = batch
        .topological_torsion_fingerprint_list_with_options(
            &TopologicalTorsionFingerprintParams::default(),
            Some(4),
            Some(false),
        )
        .expect("invalid input records remain optional outputs");
    assert_eq!(values.len(), 3);
    assert_eq!(
        values[0],
        Some(
            topological_torsion_fingerprint(
                &first,
                &TopologicalTorsionFingerprintParams::default(),
            )
            .expect("first scalar fingerprint"),
        )
    );
    assert!(values[1].is_none());
    assert_eq!(
        values[2],
        Some(
            topological_torsion_fingerprint(
                &last,
                &TopologicalTorsionFingerprintParams::default(),
            )
            .expect("last scalar fingerprint"),
        )
    );
    assert_eq!(batch.errors(), vec![original_error]);

    let invalid_params = TopologicalTorsionFingerprintParams {
        from_atoms: Some(vec![3]),
        ..TopologicalTorsionFingerprintParams::default()
    };
    for n_jobs in [Some(1), Some(2), Some(4)] {
        let error = batch
            .topological_torsion_fingerprint_list_with_options(&invalid_params, n_jobs, Some(false))
            .expect_err("two-atom record must reject atom index three");
        assert_eq!(error.errors, 1);
        assert_eq!(error.record_errors.len(), 1);
        assert_eq!(error.record_errors[0].index, 2);
        assert_eq!(
            error.record_errors[0].operation,
            "batch.topological_torsion_fingerprint"
        );
        assert!(
            error.record_errors[0]
                .message
                .contains("atom index outside the molecule"),
            "unexpected indexed error: {}",
            error.record_errors[0].message
        );
    }

    assert_eq!(batch, batch_before);
    assert_eq!(first, molecule("CCCC"));
    assert_eq!(last, molecule("CC"));
}
