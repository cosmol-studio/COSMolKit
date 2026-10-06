//! Author supplement required by ROOT CKebb; independent acceptance pending.
//! Canonical public width errors, never a public invalid-mask call.
#![cfg(all(feature = "cap-smiles", feature = "cap-fingerprints"))]

#[cfg(feature = "cap-batch")]
use cosmolkit::{BatchErrorMode, BatchParams, BatchQueryParams, MoleculeBatch};
use cosmolkit::{Molecule, PatternFingerprintError, PatternFingerprintParams};
use std::error::Error;

fn widths() -> Vec<(usize, PatternFingerprintError)> {
    let mut values = vec![(0, PatternFingerprintError::EmptyFingerprint)];
    if let Ok(width) = usize::try_from(u64::from(u32::MAX) + 1) {
        values.push((
            width,
            PatternFingerprintError::InvalidArguments {
                reason: "Pattern n_bits exceeds source unsigned int width",
            },
        ));
    }
    values
}

#[test]
fn scalar_pattern_zero_and_source_unsigned_width_errors_remain_structured() {
    // Pinned RDKit351f8f378f8ad6bbd517980c38896e66bf907af8 source
    // PatternFingerprints.cpp uses unsigned int fpSize and rejects fpSize==0.
    // Existing canonical detached owner pattern.rs:245-250 rejects zero and
    // checked narrowing outside that source width before allocating. This
    // test exercises public projections of those errors, not new chemistry.
    let molecule = Molecule::from_smiles("CCC").unwrap();
    for (n_bits, expected) in widths() {
        assert_eq!(
            molecule.pattern_fingerprint_with_params(&PatternFingerprintParams {
                n_bits,
                tautomeric: false,
            }),
            Err(expected)
        );
    }
}

#[test]
#[cfg(feature = "cap-batch")]
fn batch_pattern_zero_and_source_unsigned_width_errors_retain_all_indices_and_causes() {
    let input = vec!["CCC".to_owned(), "CC".to_owned(), "C".to_owned()];
    let batch = MoleculeBatch::from_smiles_list_with_params(
        &input,
        &Default::default(),
        &BatchParams {
            errors: BatchErrorMode::KeepErrors,
            n_jobs: Some(1),
            progress_bar: Some(false),
        },
    )
    .unwrap();
    for (n_bits, expected) in widths() {
        for n_jobs in [1, 4] {
            let error = batch
                .pattern_fingerprint_list_with_params(
                    &PatternFingerprintParams {
                        n_bits,
                        tautomeric: false,
                    },
                    &BatchQueryParams {
                        n_jobs: Some(n_jobs),
                        progress_bar: Some(false),
                        ..Default::default()
                    },
                )
                .expect_err("all three valid records must retain the same structured width error");
            assert_eq!(error.errors, 3);
            assert_eq!(
                error
                    .record_errors
                    .iter()
                    .map(|r| r.index)
                    .collect::<Vec<_>>(),
                vec![0, 1, 2]
            );
            for record in &error.record_errors {
                let mut source: Option<&(dyn Error + 'static)> = Some(record);
                let mut actual = None;
                while let Some(value) = source {
                    if let Some(cause) = value.downcast_ref::<PatternFingerprintError>() {
                        actual = Some(cause);
                    }
                    source = value.source();
                }
                assert_eq!(actual, Some(&expected));
            }
        }
    }
}
