use std::{
    any::Any,
    fs::File,
    io::{BufRead, BufReader},
    panic::{AssertUnwindSafe, catch_unwind},
    sync::Mutex,
    sync::atomic::{AtomicUsize, Ordering},
};

use rayon::iter::{ParallelBridge, ParallelIterator};

mod common;
use common::{
    parity_data,
    tautomer_parity::{GoldenRecord, assert_branch, assert_endpoint_input, parse_record},
};

fn assert_record(record: GoldenRecord) {
    let expected_parse = record.parse["ok"].as_bool().unwrap();
    match parse_record(&record) {
        Ok(molecule) => {
            assert!(expected_parse, "row {} unexpectedly parsed", record.row);
            assert_eq!(record.branches.len(), 2, "row {} profile coverage", record.row);
            for (name, branch) in &record.branches {
                assert_branch(&record, name, branch, &molecule);
            }

            let input_tautomer = record
                .input_tautomer
                .as_ref()
                .unwrap_or_else(|| panic!("row {} missing paired-input evidence", record.row));
            assert_eq!(
                Some(input_tautomer.smiles.as_str()),
                record.expected_canonical_smiles.as_deref(),
                "row {} PCS paired-input identity",
                record.row
            );
            assert_endpoint_input(record.row, "PCS paired input", input_tautomer);

            let permutation = record
                .atom_order_permutation
                .as_ref()
                .unwrap_or_else(|| panic!("row {} missing atom-order permutation evidence", record.row));
            assert_endpoint_input(record.row, "atom-order permutation", permutation);
        }
        Err(error) => {
            assert!(!expected_parse, "row {} parse failed: {error}", record.row);
            assert!(
                record.branches.is_empty(),
                "row {} failed parse must have no branches",
                record.row
            );
        }
    }
}

fn assert_profile(profile: &str, expected_records: usize) {
    let path = parity_data::expected_path_for_profile("tautomer", "rdkit", profile, "tautomer.jsonl");
    let file = File::open(&path).unwrap_or_else(|error| panic!("failed to open {}: {error}", path.display()));
    let observed = AtomicUsize::new(0);
    let first_failure = Mutex::new(None::<(usize, String)>);
    BufReader::new(file)
        .lines()
        .enumerate()
        .par_bridge()
        .for_each(|(index, line)| {
            if first_failure.lock().unwrap().is_some() {
                return;
            }
            let line =
                line.unwrap_or_else(|error| panic!("failed to read {} line {}: {error}", path.display(), index + 1));
            let record = serde_json::from_str::<GoldenRecord>(&line)
                .unwrap_or_else(|error| panic!("failed to parse {} line {}: {error}", path.display(), index + 1));
            let row = record.row;
            match catch_unwind(AssertUnwindSafe(|| assert_record(record))) {
                Ok(()) => {
                    observed.fetch_add(1, Ordering::Relaxed);
                }
                Err(payload) => {
                    let mut failure = first_failure.lock().unwrap();
                    if failure.is_none() {
                        *failure = Some((row, panic_message(payload.as_ref())));
                    }
                }
            }
        });
    if let Some((row, message)) = first_failure.into_inner().unwrap() {
        panic!("{profile} first observed failure at row {row}: {message}");
    }
    assert_eq!(
        observed.load(Ordering::Relaxed),
        expected_records,
        "{profile} record count"
    );
}

fn panic_message(payload: &(dyn Any + Send)) -> String {
    if let Some(message) = payload.downcast_ref::<String>() {
        message.clone()
    } else if let Some(message) = payload.downcast_ref::<&str>() {
        (*message).to_string()
    } else {
        "non-string panic payload".to_string()
    }
}

#[test]
#[ignore = "full 1,000-row PCS tautomer matrix; run explicitly in the exhaustive parity tier"]
fn pcs_1k_matches_pinned_rdkit_for_every_full_enumeration_and_endpoint() {
    assert_profile("tautomer_pcs_1k", 1_000);
}

#[test]
#[ignore = "full 100,000-row PCS tautomer matrix; run explicitly in the exhaustive parity tier"]
fn pcs_100k_matches_pinned_rdkit_for_every_full_enumeration_and_endpoint() {
    assert_profile("tautomer_pcs_100k", 99_991);
}
