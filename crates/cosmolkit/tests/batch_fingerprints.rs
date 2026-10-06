#![cfg(all(
    feature = "cap-batch",
    feature = "cap-smiles",
    feature = "cap-fingerprints"
))]
use cosmolkit::{
    AtomPairFingerprintParams, AtomPairParams, BatchErrorMode, BatchParams, BatchQueryParams,
    BatchRecord, FingerprintAdditionalOutput, LayeredFingerprintParams, Molecule, MoleculeBatch,
    MorganFingerprintParams, MorganParams, PatternFingerprintParams,
};
use std::error::Error;
use std::sync::{
    Arc,
    atomic::{AtomicUsize, Ordering},
};
fn smiles(values: &[&str]) -> Vec<String> {
    values.iter().map(|value| (*value).into()).collect()
}
fn batch(values: &[&str]) -> MoleculeBatch {
    MoleculeBatch::from_smiles_list_with_params(
        &smiles(values),
        &Default::default(),
        &BatchParams {
            errors: BatchErrorMode::KeepErrors,
            ..Default::default()
        },
    )
    .unwrap()
}
fn execution(jobs: usize) -> BatchQueryParams {
    BatchQueryParams {
        n_jobs: Some(jobs),
        progress_bar: Some(false),
        ..Default::default()
    }
}
fn additional_output() -> FingerprintAdditionalOutput {
    let mut output = FingerprintAdditionalOutput::new();
    output.allocate_atom_counts();
    output.allocate_atom_to_bits();
    output.allocate_bit_info_map();
    output.allocate_atoms_per_bit();
    output
}
#[test]
fn every_batch_result_form_matches_ordered_scalar_calls() {
    let source = smiles(&["C", "CCO", "c1ccccc1", "C[C@H](O)F"]);
    let molecules: Vec<_> = source
        .iter()
        .map(|s| Molecule::from_smiles(s).unwrap())
        .collect();
    let batch = MoleculeBatch::from_smiles_list(&source).unwrap();
    let params = AtomPairFingerprintParams {
        generator: AtomPairParams {
            fp_size: 256,
            count_simulation: true,
            count_bounds: vec![1, 3, 5],
            bits_per_feature: 1,
            include_chirality: true,
            ..Default::default()
        },
        ..Default::default()
    };
    let sparse_count = batch
        .fingerprint_atom_pair_sparse_count_list_with_params(&params, &execution(4))
        .unwrap();
    let expected_sparse_count: Vec<_> = molecules
        .iter()
        .map(|m| {
            Some(
                m.atom_pair_sparse_count_fingerprint_with_params(&params, None)
                    .unwrap(),
            )
        })
        .collect();
    assert_eq!(sparse_count, expected_sparse_count);
    let count = batch
        .fingerprint_atom_pair_count_list_with_params(&params, &execution(3))
        .unwrap();
    let expected_count: Vec<_> = molecules
        .iter()
        .map(|m| {
            Some(
                m.atom_pair_count_fingerprint_with_params(&params, None)
                    .unwrap(),
            )
        })
        .collect();
    assert_eq!(count, expected_count);
    let sparse_bits = batch
        .fingerprint_atom_pair_sparse_bits_list_with_params(&params, &execution(2))
        .unwrap();
    let expected_sparse_bits: Vec<_> = molecules
        .iter()
        .map(|m| {
            Some(
                m.atom_pair_sparse_fingerprint_with_params(&params, None)
                    .unwrap(),
            )
        })
        .collect();
    assert_eq!(sparse_bits, expected_sparse_bits);
    let fingerprints = batch
        .fingerprint_atom_pair_list_with_params(&params, &execution(4))
        .unwrap();
    let expected_fingerprints: Vec<_> = molecules
        .iter()
        .map(|m| Some(m.atom_pair_fingerprint_with_params(&params, None).unwrap()))
        .collect();
    assert_eq!(fingerprints, expected_fingerprints);
    let outputs = batch
        .fingerprint_atom_pair_with_output_list_with_params(&params, true, &execution(4))
        .unwrap();
    let expected_outputs: Vec<_> = molecules
        .iter()
        .map(|m| {
            let mut output = additional_output();
            let fingerprint = m
                .atom_pair_fingerprint_with_params(&params, Some(&mut output))
                .unwrap();
            Some(cosmolkit_fingerprints::batch_fingerprint_output(
                fingerprint,
                Some(&output),
                "AtomPair",
            ))
        })
        .collect();
    assert_eq!(outputs, expected_outputs);
}
#[test]
fn unfolded_extra_bits_preserve_source_index_errors_and_batch_positions() {
    let source = smiles(&["C", "CCO", "c1ccccc1", "C[C@H](O)F"]);
    let molecules: Vec<_> = source
        .iter()
        .map(|s| Molecule::from_smiles(s).unwrap())
        .collect();
    let batch = MoleculeBatch::from_smiles_list(&source).unwrap();
    let params = AtomPairFingerprintParams {
        generator: AtomPairParams {
            bits_per_feature: 2,
            include_chirality: true,
            ..Default::default()
        },
        ..Default::default()
    };
    assert!(
        molecules[0]
            .atom_pair_sparse_count_fingerprint_with_params(&params, None)
            .is_ok(),
        "a graph with no atom-pair environment does not consume an extra bit"
    );
    let expected_indices = [1_795_012_513_u64, 601_524_142, 525_510_353];
    for (molecule, expected_index) in molecules[1..].iter().zip(expected_indices) {
        let error = molecule
            .atom_pair_sparse_count_fingerprint_with_params(&params, None)
            .expect_err("RDKit's unfolded extra bit must exceed the chiral result width");
        assert!(
            error.to_string().contains(&expected_index.to_string()),
            "unexpected scalar source-parity error: {error}"
        );
    }
    let error = batch
        .fingerprint_atom_pair_sparse_count_list_with_params(&params, &execution(4))
        .expect_err("batch collection must retain every source-defined index error");
    assert_eq!(error.errors, 3);
    assert_eq!(
        error
            .record_errors
            .iter()
            .map(|e| e.index)
            .collect::<Vec<_>>(),
        vec![1, 2, 3]
    );
    for (error, expected_index) in error.record_errors.iter().zip(expected_indices) {
        assert_eq!(error.operation, "batch.atom_pair_sparse_count_fingerprint");
        assert!(
            error.message.contains(&expected_index.to_string()),
            "unexpected indexed source-parity error: {}",
            error.message
        );
        assert!(error.source().is_some());
    }
}
#[test]
fn batch_order_thread_count_progress_and_repeated_calls_are_deterministic() {
    let source = smiles(&["C", "CC", "CCC", "CCCC", "CCCCC", "CCCCCC"]);
    let batch = MoleculeBatch::from_smiles_list(&source)
        .unwrap()
        .with_parallel_jobs(Some(4))
        .unwrap();
    let params = AtomPairFingerprintParams::default();
    let baseline = batch
        .fingerprint_atom_pair_list_with_params(&params, &execution(1))
        .unwrap();
    for n_jobs in [1, 2, 4] {
        for _ in 0..3 {
            assert_eq!(
                batch
                    .fingerprint_atom_pair_list_with_params(&params, &execution(n_jobs))
                    .unwrap(),
                baseline
            );
        }
    }
    let ticks = Arc::new(AtomicUsize::new(0));
    let counter = Arc::clone(&ticks);
    let progress = BatchQueryParams {
        progress_callback: Some(Arc::new(move || {
            counter.fetch_add(1, Ordering::Relaxed);
        })),
        ..Default::default()
    };
    assert_eq!(
        batch
            .fingerprint_atom_pair_list_with_params(&params, &progress)
            .unwrap(),
        baseline
    );
    assert_eq!(ticks.load(Ordering::Relaxed), source.len());
    assert_eq!(batch.parallel_jobs(), Some(4));
}
#[test]
fn invalid_input_records_keep_their_indices_and_operation_errors_are_indexed() {
    let batch = batch(&["CC", "C1", "CCC"]);
    assert_eq!(batch.errors().len(), 1);
    assert_eq!(batch.errors()[0].index, 1);
    let values = batch
        .fingerprint_atom_pair_list_with_params(&Default::default(), &execution(2))
        .unwrap();
    assert!(values[0].is_some());
    assert!(values[1].is_none());
    assert!(values[2].is_some());
    let no_conformer = AtomPairFingerprintParams {
        generator: AtomPairParams {
            use_2d: false,
            ..Default::default()
        },
        ..Default::default()
    };
    let error = batch
        .fingerprint_atom_pair_list_with_params(&no_conformer, &execution(2))
        .unwrap_err();
    assert_eq!(
        error
            .record_errors
            .iter()
            .map(|e| e.index)
            .collect::<Vec<_>>(),
        vec![0, 2]
    );
    assert!(
        error
            .record_errors
            .iter()
            .all(|e| e.operation == "batch.atom_pair_fingerprint")
    );
    let bad_config = AtomPairFingerprintParams {
        generator: AtomPairParams {
            fp_size: 0,
            ..Default::default()
        },
        ..Default::default()
    };
    let error = batch
        .fingerprint_atom_pair_list_with_params(&bad_config, &execution(2))
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
    let batch = Arc::new(MoleculeBatch::from_smiles_list(&source).unwrap());
    let snapshot = batch.to_list();
    let expected = batch
        .fingerprint_atom_pair_list_with_params(&Default::default(), &execution(2))
        .unwrap();
    let atom_pair_batch = Arc::clone(&batch);
    let atom_pair = std::thread::spawn(move || {
        atom_pair_batch
            .fingerprint_atom_pair_list_with_params(&Default::default(), &execution(4))
            .unwrap()
    });
    let morgan_batch = Arc::clone(&batch);
    let morgan = std::thread::spawn(move || {
        morgan_batch
            .fingerprint_morgan_list_with_params(&Default::default(), &execution(3))
            .unwrap()
    });
    assert_eq!(atom_pair.join().unwrap(), expected);
    assert_eq!(morgan.join().unwrap().len(), source.len());
    for (before, after) in snapshot.iter().zip(batch.to_list()) {
        let before = before.as_ref().unwrap();
        let after = after.unwrap();
        assert_eq!(after.topology(), before.topology());
        assert_eq!(after.properties(), before.properties());
        assert_eq!(after.coordinates_2d(), before.coordinates_2d());
        assert_eq!(after.conformers_3d(), before.conformers_3d());
    }
    assert!(batch.iter().all(|r| matches!(r, BatchRecord::Molecule(_))));
}
#[test]
fn all_families_return_full_ordered_bits_counts_and_metadata_masks() {
    let batch = batch(&["CCO", "C1", "c1ccccc1O", "C[C@H](O)F"]);
    let originals = batch.to_list();
    let ap = AtomPairFingerprintParams::default();
    let morgan = MorganFingerprintParams {
        generator: MorganParams {
            radius: 2,
            fp_size: 128,
            ..Default::default()
        },
        ..Default::default()
    };
    let layered = LayeredFingerprintParams {
        fp_size: 128,
        atom_counts: Some(vec![0; 12]),
        ..Default::default()
    };
    let pattern = PatternFingerprintParams {
        n_bits: 128,
        ..Default::default()
    };
    for jobs in [1, 2, 4] {
        let query = execution(jobs);
        for collect in [false, true] {
            let ap_results = batch
                .fingerprint_atom_pair_with_output_list_with_params(&ap, collect, &query)
                .unwrap();
            let mo_results = batch
                .fingerprint_morgan_with_output_list_with_params(&morgan, collect, &query)
                .unwrap();
            for (index, molecule) in originals.iter().enumerate() {
                if let Some(m) = molecule {
                    let mut ao = collect.then(additional_output);
                    let fp = m
                        .atom_pair_fingerprint_with_params(&ap, ao.as_mut())
                        .unwrap();
                    assert_eq!(
                        ap_results[index],
                        Some(cosmolkit_fingerprints::batch_fingerprint_output(
                            fp,
                            ao.as_ref(),
                            "AtomPair"
                        ))
                    );
                    let mut ao = collect.then(additional_output);
                    let fp = m
                        .morgan_fingerprint_with_params(&morgan, ao.as_mut())
                        .unwrap();
                    assert_eq!(
                        mo_results[index],
                        Some(cosmolkit_fingerprints::batch_fingerprint_output(
                            fp,
                            ao.as_ref(),
                            "Morgan"
                        ))
                    );
                    if collect {
                        for output in [&ap_results[index], &mo_results[index]] {
                            let output =
                                output.as_ref().unwrap().additional_output.as_ref().unwrap();
                            assert!(output.atom_counts().is_some());
                            assert!(output.atom_to_bits().is_some());
                            assert!(output.bit_info_map().is_some());
                            assert!(output.atoms_per_bit().is_some());
                            assert!(output.bit_paths().is_none());
                        }
                    }
                } else {
                    assert!(ap_results[index].is_none());
                    assert!(mo_results[index].is_none());
                }
            }
        }
        let bits = batch
            .fingerprint_layered_list_with_params(&layered, &query)
            .unwrap();
        let output = batch
            .fingerprint_layered_with_output_list_with_params(&layered, &query)
            .unwrap();
        let pattern_bits = batch
            .pattern_fingerprint_list_with_params(&pattern, &query)
            .unwrap();
        for (index, molecule) in originals.iter().enumerate() {
            if let Some(m) = molecule {
                assert_eq!(
                    bits[index],
                    Some(m.layered_fingerprint_with_params(&layered).unwrap())
                );
                assert_eq!(
                    output[index],
                    Some(
                        m.layered_fingerprint_with_output_with_params(&layered)
                            .unwrap()
                    )
                );
                assert_eq!(
                    pattern_bits[index],
                    Some(m.pattern_fingerprint_with_params(&pattern).unwrap())
                );
            } else {
                assert!(bits[index].is_none());
                assert!(output[index].is_none());
                assert!(pattern_bits[index].is_none());
            }
        }
    }
}
#[test]
fn ap_preflight_runs_on_empty_invalid_rows_before_scheduling_and_morgan_preserves_per_record_validation()
 {
    for values in [vec![], vec!["C1", "C2"]] {
        let batch = batch(&values);
        let params = AtomPairFingerprintParams {
            generator: AtomPairParams {
                fp_size: 0,
                min_distance: 31,
                max_distance: 30,
                ..Default::default()
            },
            ..Default::default()
        };
        let error = batch
            .fingerprint_atom_pair_list_with_params(&params, &execution(0))
            .unwrap_err();
        let owner =
            cosmolkit_fingerprints::validate_atom_pair_params(&params.generator).unwrap_err();
        assert_eq!(error.record_errors[0].index, 0);
        assert!(error.record_errors[0].message.contains(&owner.to_string()));
        assert!(
            error.record_errors[0]
                .source()
                .unwrap()
                .downcast_ref::<cosmolkit_fingerprints::FingerprintError>()
                .is_some()
        );
        assert_eq!(
            error.record_errors[0].operation,
            "batch.atom_pair_fingerprint"
        );
        let bad = MorganFingerprintParams {
            generator: MorganParams {
                fp_size: 0,
                ..Default::default()
            },
            ..Default::default()
        };
        assert_eq!(
            batch
                .fingerprint_morgan_list_with_params(&bad, &execution(1))
                .unwrap()
                .len(),
            values.len()
        );
    }
    let valid = batch(&["CCO"]);
    assert!(
        valid
            .fingerprint_morgan_list_with_params(
                &MorganFingerprintParams {
                    generator: MorganParams {
                        fp_size: 0,
                        ..Default::default()
                    },
                    ..Default::default()
                },
                &execution(1)
            )
            .is_err()
    );
}
#[test]
fn original_morgan_explicit_providers_defaults_and_empty_options_match_scalar_owner() {
    use cosmolkit::{
        MorganAtomInvariantsGenerator, MorganBondInvariantsGenerator, MorganCallParams,
        MorganFingerprintGenerator,
    };
    let batch = batch(&["CCO", "C1", "c1ccccc1O", "C[C@H](O)F"]);
    let originals = batch.to_list();
    for atom in [
        MorganAtomInvariantsGenerator::connectivity(false),
        MorganAtomInvariantsGenerator::features(None),
    ] {
        for bond in [None, Some(MorganBondInvariantsGenerator::new(false, true))] {
            let params = MorganParams {
                radius: 2,
                fp_size: 128,
                ..Default::default()
            };
            let call = MorganCallParams::default();
            let generator =
                MorganFingerprintGenerator::new(Some(&params), Some(&atom), bond.as_ref()).unwrap();
            let outputs = batch
                .fingerprint_morgan_with_output_list_with_generator_params(
                    &params,
                    Some(&atom),
                    bond.as_ref(),
                    &call,
                    true,
                    &execution(4),
                )
                .unwrap();
            let bits = batch
                .fingerprint_morgan_list_with_generator_params(
                    &params,
                    Some(&atom),
                    bond.as_ref(),
                    &call,
                    &execution(1),
                )
                .unwrap();
            for (index, molecule) in originals.iter().enumerate() {
                if let Some(m) = molecule {
                    let mut ao = additional_output();
                    let fingerprint = m
                        .morgan_fingerprint_with_generator(&generator, Some(&call), Some(&mut ao))
                        .unwrap();
                    assert_eq!(bits[index], Some(fingerprint.clone()));
                    assert_eq!(
                        outputs[index],
                        Some(cosmolkit_fingerprints::batch_fingerprint_output(
                            fingerprint,
                            Some(&ao),
                            "Morgan"
                        ))
                    );
                } else {
                    assert!(outputs[index].is_none());
                    assert!(bits[index].is_none());
                }
            }
        }
    }
    let mut defaults = MorganFingerprintParams::default();
    defaults.generator.radius = 2;
    assert_eq!(
        batch.fingerprint_morgan_list().unwrap(),
        batch
            .fingerprint_morgan_list_with_params(&defaults, &execution(1))
            .unwrap()
    );
}

#[test]
fn legacy_morgan_batch_zero_size_and_zero_extra_bits_follow_wrapper_order() {
    let batch = batch(&["CCO", "C1CC", "C"]);
    let zero = MorganFingerprintParams {
        generator: MorganParams {
            fp_size: 0,
            count_simulation: false,
            ..Default::default()
        },
        ..Default::default()
    };
    let error = batch
        .fingerprint_morgan_list_with_params(&zero, &execution(2))
        .unwrap_err();
    assert_eq!(
        error
            .record_errors
            .iter()
            .map(|e| e.index)
            .collect::<Vec<_>>(),
        [0, 2]
    );
    assert!(
        error
            .record_errors
            .iter()
            .all(|e| e.message.contains("fingerprint requires n_bits > 0"))
    );
    let zero_bits = MorganFingerprintParams {
        generator: MorganParams {
            bits_per_feature: 0,
            radius: 2,
            ..Default::default()
        },
        ..Default::default()
    };
    let actual = batch
        .fingerprint_morgan_list_with_params(&zero_bits, &execution(2))
        .unwrap();
    let generator = cosmolkit::MorganFingerprintGenerator::new(
        Some(&MorganParams {
            radius: 2,
            ..Default::default()
        }),
        None,
        None,
    )
    .unwrap();
    generator.settings().set_bits_per_feature(0).unwrap();
    for (i, item) in actual.iter().enumerate() {
        if let Some(BatchRecord::Molecule(m)) = batch.get(i) {
            assert_eq!(
                *item,
                Some(
                    m.morgan_fingerprint_with_generator(&generator, None, None)
                        .unwrap()
                )
            );
        } else {
            assert!(item.is_none());
        }
    }
}
