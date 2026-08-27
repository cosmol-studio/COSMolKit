use std::{
    fs::File,
    io::{BufRead, BufReader},
};

use cosmolkit_core::{MaccsFingerprintParams, calc_exact_mol_wt, mol_from_binary, mol_to_binary};
use rayon::prelude::*;

mod common;
use common::{
    parity_data,
    tautomer_parity::{
        GoldenBranch, GoldenRecord, assert_branch, configured_enumerator, molecule_state,
        parse_record,
    },
};

fn profile_records(profile: &str) -> impl Iterator<Item = GoldenRecord> {
    let path =
        parity_data::expected_path_for_profile("tautomer", "rdkit", profile, "tautomer.jsonl");
    let file = File::open(&path)
        .unwrap_or_else(|error| panic!("failed to open {}: {error}", path.display()));
    BufReader::new(file)
        .lines()
        .enumerate()
        .map(move |(index, line)| {
            let line = line.unwrap_or_else(|error| {
                panic!(
                    "failed to read {} line {}: {error}",
                    path.display(),
                    index + 1
                )
            });
            serde_json::from_str(&line).unwrap_or_else(|error| {
                panic!(
                    "failed to parse {} line {}: {error}",
                    path.display(),
                    index + 1
                )
            })
        })
}

fn assert_record(record: &GoldenRecord, expected_branch_count: usize) {
    let expected_parse = record.parse["ok"].as_bool().unwrap();
    match parse_record(record) {
        Ok(molecule) => {
            assert!(expected_parse, "row {} unexpectedly parsed", record.row);
            assert_eq!(
                record.branches.len(),
                expected_branch_count,
                "row {} profile coverage",
                record.row
            );
            let input_state = molecule_state(&molecule);
            for (name, branch) in &record.branches {
                assert_branch(record, name, branch, &molecule);
            }
            assert_eq!(
                molecule_state(&molecule),
                input_state,
                "row {} source molecule changed",
                record.row
            );
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

fn assert_composition(record: &GoldenRecord, branch_name: &str, branch: &GoldenBranch) {
    if !branch.ok {
        return;
    }
    let molecule = parse_record(record)
        .unwrap_or_else(|error| panic!("row {} composition parse failed: {error}", record.row));
    let source_state = molecule_state(&molecule);
    let enumerator = configured_enumerator(&branch.parameters);
    let result = enumerator.enumerate(&molecule).unwrap_or_else(|error| {
        panic!(
            "row {} branch {branch_name} composition enumeration failed: {error}",
            record.row
        )
    });
    assert_eq!(result.len(), branch.molecule_states.len());

    for (index, (tautomer, expected_state)) in
        result.iter().zip(branch.molecule_states.iter()).enumerate()
    {
        assert_eq!(
            molecule_state(tautomer),
            *expected_state,
            "row {} branch {branch_name} tautomer {index} pre-composition state",
            record.row
        );

        let first_mass = calc_exact_mol_wt(tautomer, false).unwrap_or_else(|error| {
            panic!(
                "row {} branch {branch_name} tautomer {index} exact mass failed: {error}",
                record.row
            )
        });
        let second_mass = calc_exact_mol_wt(tautomer, false).unwrap();
        assert_eq!(first_mass.to_bits(), second_mass.to_bits());

        let first_fingerprint = tautomer
            .maccs_fingerprint(&MaccsFingerprintParams::default())
            .unwrap_or_else(|error| {
                panic!(
                    "row {} branch {branch_name} tautomer {index} fingerprint failed: {error}",
                    record.row
                )
            });
        let second_fingerprint = tautomer
            .maccs_fingerprint(&MaccsFingerprintParams::default())
            .unwrap();
        assert_eq!(first_fingerprint, second_fingerprint);

        let binary = mol_to_binary(tautomer).unwrap_or_else(|error| {
            panic!(
                "row {} branch {branch_name} tautomer {index} binary write failed: {error}",
                record.row
            )
        });
        let restored = mol_from_binary(&binary).unwrap_or_else(|error| {
            panic!(
                "row {} branch {branch_name} tautomer {index} binary read failed: {error}",
                record.row
            )
        });
        assert_eq!(
            molecule_state(&restored),
            *expected_state,
            "row {} branch {branch_name} tautomer {index} binary composition state",
            record.row
        );
        assert_eq!(
            molecule_state(tautomer),
            *expected_state,
            "row {} branch {branch_name} tautomer {index} post-composition state",
            record.row
        );
    }

    assert_eq!(
        molecule_state(&molecule),
        source_state,
        "row {} branch {branch_name} composition changed the source molecule",
        record.row
    );
}

#[test]
fn focused_profiles_are_repeatable_and_compose_without_state_drift() {
    let records = profile_records("tautomer_focused").collect::<Vec<_>>();
    assert_eq!(records.len(), 18, "focused fixture row count changed");
    for record in &records {
        assert_record(record, 8);
        assert_record(record, 8);
        for (name, branch) in &record.branches {
            assert_composition(record, name, branch);
        }
    }
}

#[test]
#[ignore = "full 5,000-row tautomer matrix; run explicitly in the exhaustive parity tier"]
fn smiles_5000_profiles_match_every_pinned_rdkit_observation() {
    let mut observed = 0usize;
    for record in profile_records("smiles_5000") {
        assert_record(&record, 2);
        observed += 1;
    }
    assert_eq!(observed, 5_000, "5,000-row corpus size changed");
}

#[test]
#[ignore = "full 5,000-row tautomer batch matrix; run explicitly in the exhaustive parity tier"]
fn smiles_5000_parallel_batches_preserve_order_and_exact_results() {
    let mut records = profile_records("smiles_5000");
    let mut observed_rows = Vec::with_capacity(5_000);
    loop {
        let batch = records.by_ref().take(128).collect::<Vec<_>>();
        if batch.is_empty() {
            break;
        }
        observed_rows.extend(
            batch
                .par_iter()
                .map(|record| {
                    assert_record(record, 2);
                    record.row
                })
                .collect::<Vec<_>>(),
        );
    }
    assert_eq!(
        observed_rows,
        (0..5_000).collect::<Vec<_>>(),
        "parallel batch collection changed corpus order"
    );
}
