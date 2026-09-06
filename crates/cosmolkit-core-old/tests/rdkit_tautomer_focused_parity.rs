use std::{
    fs::File,
    io::{BufRead, BufReader},
};

use cosmolkit_core::{BondOrder, TautomerCatalog};
use serde::Deserialize;

mod common;
use common::{
    parity_data,
    tautomer_parity::{GoldenRecord, assert_branch, parse_record},
};

const PROFILE: &str = "tautomer_focused";

#[derive(Debug, Deserialize)]
struct CatalogRecord {
    catalog: String,
    index: usize,
    source_path: String,
    source_line: usize,
    name: String,
    smarts: String,
    bonds: String,
    charges: String,
}

fn read_jsonl<T: for<'de> Deserialize<'de>>(output: &str) -> Vec<T> {
    let path = parity_data::expected_path_for_profile("tautomer", "rdkit", PROFILE, output);
    let file = File::open(&path)
        .unwrap_or_else(|error| panic!("failed to open {}: {error}", path.display()));
    BufReader::new(file)
        .lines()
        .enumerate()
        .map(|(index, line)| {
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
        .collect()
}

#[test]
fn focused_tautomer_profiles_match_every_pinned_rdkit_observation() {
    let records = read_jsonl::<GoldenRecord>("tautomer.jsonl");
    assert_eq!(records.len(), 18, "focused fixture row count changed");
    for record in records {
        let expected_parse = record.parse["ok"].as_bool().unwrap();
        match parse_record(&record) {
            Ok(molecule) => {
                assert!(expected_parse, "row {} unexpectedly parsed", record.row);
                assert_eq!(
                    record.branches.len(),
                    8,
                    "row {} profile coverage",
                    record.row
                );
                for (name, branch) in &record.branches {
                    assert_branch(&record, name, branch, &molecule);
                }
            }
            Err(error) => {
                assert!(!expected_parse, "row {} parse failed: {error}", record.row);
                assert!(
                    record.branches.is_empty(),
                    "failed parse must have no branches"
                );
            }
        }
    }
}

fn expected_bond_types(source: &str) -> Vec<BondOrder> {
    source
        .bytes()
        .filter_map(|byte| match byte {
            b'-' => Some(BondOrder::Single),
            b'=' => Some(BondOrder::Double),
            b'#' => Some(BondOrder::Triple),
            b':' => Some(BondOrder::Aromatic),
            _ => None,
        })
        .collect()
}

fn expected_charges(source: &str) -> Vec<i32> {
    source
        .bytes()
        .map(|byte| match byte {
            b'+' => 1,
            b'0' => 0,
            b'-' => -1,
            _ => panic!("golden catalog contains an invalid charge byte"),
        })
        .collect()
}

#[test]
fn current_and_v1_catalogs_match_every_pinned_rdkit_transform() {
    let records = read_jsonl::<CatalogRecord>("tautomer_catalog.jsonl");
    let catalogs = [
        (
            "current",
            TautomerCatalog::current().expect("compile current catalog"),
            37,
        ),
        ("v1", TautomerCatalog::v1().expect("compile V1 catalog"), 36),
    ];
    for (name, catalog, expected_len) in catalogs {
        let expected = records
            .iter()
            .filter(|record| record.catalog == name)
            .collect::<Vec<_>>();
        assert_eq!(expected.len(), expected_len, "{name} golden count");
        assert_eq!(
            catalog.transforms().len(),
            expected_len,
            "{name} compiled count"
        );
        for (index, (actual, expected)) in catalog.transforms().iter().zip(expected).enumerate() {
            let context = format!(
                "{name} transform {index} at {}:{} ({})",
                expected.source_path, expected.source_line, expected.smarts
            );
            assert_eq!(expected.index, index, "{context}: source order");
            assert_eq!(actual.name(), expected.name, "{context}: name");
            assert_eq!(
                actual.bond_types(),
                expected_bond_types(&expected.bonds),
                "{context}: bonds"
            );
            assert_eq!(
                actual.charges(),
                expected_charges(&expected.charges),
                "{context}: charges"
            );
        }
    }
}
