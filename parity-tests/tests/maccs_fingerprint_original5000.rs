use std::fs::File;
use std::io::{BufRead, BufReader};

use cosmolkit::{MaccsFingerprintError, MaccsFingerprintParams, Molecule};
use serde::Deserialize;
use sha2::{Digest, Sha256};

fn source_bits(value: &cosmolkit::Fingerprint) -> Vec<usize> {
    value
        .on_bits()
        .into_iter()
        .map(|bit| bit as usize)
        .collect()
}

#[derive(Debug, Deserialize)]
struct MaccsRecord {
    record_type: String,
    label: Option<String>,
    smiles: String,
    rdkit_ok: bool,
    raw_n_bits: usize,
    public_n_bits: usize,
    raw_on_bits: Option<Vec<usize>>,
    public_on_bits: Option<Vec<usize>>,
    error: Option<String>,
}

fn load_golden() -> Vec<MaccsRecord> {
    let path = std::path::PathBuf::from(
        std::env::var("COSMOLKIT_MACCS_ORIGINAL5000_GOLDEN")
            .expect("prepare unchanged original MACCS source reference"),
    );
    let bytes = std::fs::read(&path).unwrap();
    assert_eq!(
        Sha256::digest(&bytes)
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect::<String>(),
        "5e94d3b504cbdfab3a0c6b65f8d37c696f444b58625d1b894a69e0358f675687",
        "immutable unchanged original source reference"
    );
    let file = File::open(&path).unwrap_or_else(|err| {
        panic!("failed to open {}: {err}", path.display());
    });
    BufReader::new(file)
        .lines()
        .enumerate()
        .map(|(idx, line)| {
            let line = line.unwrap_or_else(|err| {
                panic!("failed to read {} line {}: {err}", path.display(), idx + 1)
            });
            serde_json::from_str(&line).unwrap_or_else(|err| {
                panic!("failed to parse {} line {}: {err}", path.display(), idx + 1)
            })
        })
        .collect()
}

fn record_context(record: &MaccsRecord, row_idx: usize) -> String {
    match (&record.record_type[..], &record.label) {
        ("fixture", Some(label)) => {
            format!("fixture {label} row {} ({})", row_idx + 1, record.smiles)
        }
        _ => format!(
            "{} row {} ({})",
            record.record_type,
            row_idx + 1,
            record.smiles
        ),
    }
}

#[test]
fn maccs_fingerprint_golden_has_profile_corpus_and_targeted_fixtures() {
    let records = load_golden();
    let corpus_count = records
        .iter()
        .filter(|record| record.record_type == "corpus")
        .count();
    let fixture_count = records
        .iter()
        .filter(|record| record.record_type == "fixture")
        .count();

    assert_eq!(
        corpus_count, 5000,
        "MACCS corpus golden row count must match the selected SMILES profile"
    );
    assert!(
        fixture_count >= 80,
        "MACCS golden should include targeted key fixtures, found {fixture_count}"
    );
}

#[test]
fn maccs_fingerprint_matches_rdkit_raw_and_public_golden() {
    let records = load_golden();

    for (row_idx, record) in records.iter().enumerate() {
        let context = record_context(record, row_idx);
        assert_eq!(
            record.raw_n_bits, 167,
            "RDKit MACCS raw vector size mismatch in {context}"
        );
        assert_eq!(
            record.public_n_bits, 166,
            "COSMolKit MACCS public projection size mismatch in {context}"
        );

        if !record.rdkit_ok {
            assert!(
                record.error.is_some(),
                "RDKit-not-ok MACCS record has no error in {context}"
            );
            continue;
        }

        let mol = Molecule::from_smiles(&record.smiles)
            .unwrap_or_else(|err| panic!("COSMolKit failed to parse {context}: {err}"));

        let expected_raw = record
            .raw_on_bits
            .as_ref()
            .unwrap_or_else(|| panic!("RDKit-ok MACCS record missing raw_on_bits in {context}"));
        let raw = mol.maccs_fingerprint_raw().unwrap_or_else(|err| {
            panic!("COSMolKit failed raw MACCS generation in {context}: {err}")
        });
        assert_eq!(raw.n_bits(), 167, "raw MACCS n_bits mismatch in {context}");
        assert!(
            !raw.on_bits().contains(&0),
            "RDKit MACCS raw bit 0 must remain unused in {context}"
        );
        assert_eq!(
            source_bits(&raw),
            *expected_raw,
            "MACCS raw 167-bit vector must be RDKit bit-identical in {context}"
        );

        let expected_public = record
            .public_on_bits
            .as_ref()
            .unwrap_or_else(|| panic!("RDKit-ok MACCS record missing public_on_bits in {context}"));
        let public = mol
            .maccs_fingerprint_with_params(&MaccsFingerprintParams::default())
            .expect("MACCS public fingerprint");
        assert_eq!(
            public.n_bits(),
            166,
            "public MACCS n_bits mismatch in {context}"
        );
        assert_eq!(
            source_bits(&public),
            *expected_public,
            "MACCS public 166-bit projection must be RDKit bit-identical in {context}"
        );

        let err = mol
            .maccs_fingerprint_with_params(&MaccsFingerprintParams { n_bits: 64 })
            .unwrap_err();
        assert!(
            matches!(
                err,
                MaccsFingerprintError::UnsupportedOption {
                    option: "MaccsFingerprintParams.n_bits",
                    ..
                }
            ),
            "MACCS non-166-bit request must be rejected at its option boundary in {context}: {err}"
        );
    }
}
