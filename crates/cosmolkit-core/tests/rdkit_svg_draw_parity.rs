use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use cosmolkit_core::{Molecule, MoleculeBatch};
use serde::Deserialize;

mod common;
use common::parity_data;

#[derive(Debug, Deserialize)]
struct SvgDrawRecord {
    smiles: String,
    rdkit_ok: bool,
    width: u32,
    height: u32,
    svg: Option<String>,
    error: Option<String>,
}

fn load_golden() -> Vec<SvgDrawRecord> {
    let path = parity_data::golden_path("svg_drawer.jsonl");
    let file = File::open(&path).unwrap_or_else(|err| {
        panic!(
            "failed to open {}; regenerate RDKit goldens with `{}`: {err}",
            path.display(),
            parity_data::regenerate_command()
        )
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

fn normalize_svg_tool_identifiers(svg: &str) -> String {
    svg.replace(
        "xmlns:rdkit='http://www.rdkit.org/xml'",
        "xmlns:tool='__tool_namespace__'",
    )
    .replace(
        "xmlns:cosmolkit='https://www.cosmol.org'",
        "xmlns:tool='__tool_namespace__'",
    )
    .replace("rdkit:", "tool:")
    .replace("cosmolkit:", "tool:")
}

fn maybe_dump_svg_debug(row_idx: usize, expected_svg: &str, actual_svg: &str) {
    let Some(target_row) = std::env::var("COSMOLKIT_SVG_DEBUG_DUMP_ROW")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
    else {
        return;
    };
    if row_idx + 1 != target_row {
        return;
    }
    let dump_dir =
        std::env::var("COSMOLKIT_SVG_DEBUG_DUMP_DIR").unwrap_or_else(|_| "/tmp".to_string());
    let dump_dir = PathBuf::from(dump_dir);
    std::fs::create_dir_all(&dump_dir).unwrap_or_else(|err| {
        panic!(
            "failed to create svg debug dump dir {}: {err}",
            dump_dir.display()
        )
    });
    let expected_path = dump_dir.join(format!("row{}_expected.svg", target_row));
    let actual_path = dump_dir.join(format!("row{}_actual.svg", target_row));
    std::fs::write(&expected_path, expected_svg)
        .unwrap_or_else(|err| panic!("failed to write {}: {err}", expected_path.display()));
    std::fs::write(&actual_path, actual_svg)
        .unwrap_or_else(|err| panic!("failed to write {}: {err}", actual_path.display()));
}

#[test]
fn svg_draw_golden_has_one_record_per_smiles() {
    let expected = parity_data::count_smiles_rows();
    let records = load_golden();
    assert_eq!(
        records.len(),
        expected,
        "SVG golden row count must match the active parity corpus"
    );
}

#[test]
fn svg_drawer_matches_rdkit_golden_except_tool_identifiers() {
    let records = load_golden();
    let row_filter = std::env::var("COSMOLKIT_SVG_ROW_FILTER")
        .ok()
        .and_then(|s| s.parse::<usize>().ok());

    for (row_idx, record) in records.iter().enumerate() {
        if let Some(filter) = row_filter {
            if row_idx + 1 != filter {
                continue;
            }
        }
        if !record.rdkit_ok {
            assert!(
                record.error.is_some(),
                "row {} ({}) is rdkit not ok but has no error",
                row_idx + 1,
                record.smiles
            );
            continue;
        }

        let mol = Molecule::from_smiles(&record.smiles).unwrap_or_else(|err| {
            panic!(
                "cosmolkit failed to parse row {} ({}): {err}",
                row_idx + 1,
                record.smiles
            )
        });
        let actual_svg = mol
            .to_svg(record.width, record.height)
            .unwrap_or_else(|err| {
                panic!(
                    "cosmolkit failed to draw row {} ({}): {err}",
                    row_idx + 1,
                    record.smiles
                )
            });
        let expected_svg = record.svg.as_ref().unwrap_or_else(|| {
            panic!(
                "row {} ({}) is rdkit ok but has no svg payload",
                row_idx + 1,
                record.smiles
            )
        });
        let normalized_actual = normalize_svg_tool_identifiers(&actual_svg);
        let normalized_expected = normalize_svg_tool_identifiers(expected_svg);
        maybe_dump_svg_debug(row_idx, &normalized_expected, &normalized_actual);

        assert_eq!(
            normalized_actual,
            normalized_expected,
            "SVG mismatch at row {} ({}) after normalizing tool identifiers",
            row_idx + 1,
            record.smiles
        );
    }
}

#[test]
fn svg_drawer_matches_rdkit_golden_except_tool_identifiers_in_parallel_batch() {
    let records = load_golden();
    let smiles = records
        .iter()
        .map(|record| record.smiles.clone())
        .collect::<Vec<_>>();
    let batch = MoleculeBatch::from_smiles_list(&smiles).with_parallel_jobs(Some(4));
    let actual_svgs = batch
        .to_svg_list_with_options(300, 300, Some(4), Some(false))
        .expect("parallel batch SVG drawing should collect without draw errors");

    assert_eq!(actual_svgs.len(), records.len());
    for (row_idx, (record, actual_svg)) in records.iter().zip(actual_svgs.iter()).enumerate() {
        if !record.rdkit_ok {
            assert!(
                record.error.is_some(),
                "row {} ({}) is rdkit not ok but has no error",
                row_idx + 1,
                record.smiles
            );
            assert!(
                actual_svg.is_none(),
                "parallel batch SVG should not produce row {} ({})",
                row_idx + 1,
                record.smiles
            );
            continue;
        }

        assert_eq!(
            (record.width, record.height),
            (300, 300),
            "parallel batch SVG golden has unexpected dimensions at row {} ({})",
            row_idx + 1,
            record.smiles
        );
        let actual_svg = actual_svg.as_ref().unwrap_or_else(|| {
            panic!(
                "parallel batch SVG missing row {} ({})",
                row_idx + 1,
                record.smiles
            )
        });
        let expected_svg = record.svg.as_ref().unwrap_or_else(|| {
            panic!(
                "row {} ({}) is rdkit ok but has no svg payload",
                row_idx + 1,
                record.smiles
            )
        });
        let normalized_actual = normalize_svg_tool_identifiers(actual_svg);
        let normalized_expected = normalize_svg_tool_identifiers(expected_svg);
        maybe_dump_svg_debug(row_idx, &normalized_expected, &normalized_actual);

        assert_eq!(
            normalized_actual,
            normalized_expected,
            "parallel batch SVG mismatch at row {} ({}) after normalizing tool identifiers",
            row_idx + 1,
            record.smiles
        );
    }
}
