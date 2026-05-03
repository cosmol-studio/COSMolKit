use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use cosmolkit_core::{BatchErrorMode, Molecule, MoleculeBatch};
use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct SvgDrawerRecord {
    smiles: String,
    rdkit_ok: bool,
    width: u32,
    height: u32,
    svg: Option<String>,
    error: Option<String>,
}

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn load_golden() -> Vec<SvgDrawerRecord> {
    let path = repo_root().join("tests/golden/svg_drawer.jsonl");
    let file = File::open(&path).unwrap_or_else(|err| {
        panic!(
            "failed to open {}; regenerate with tests/scripts/gen_rdkit_svg_golden.py: {err}",
            path.display()
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

fn normalize_svg_identity(svg: &str) -> String {
    svg.replace(
        "xmlns:rdkit='http://www.rdkit.org/xml'",
        "xmlns:cosmolkit='https://www.cosmol.org'",
    )
    .replace(
        "font-family:sans-serif",
        "font-family:\"Noto Sans\",sans-serif",
    )
}

fn svg_diff(actual: &str, expected: &str) -> String {
    const MAX_DIFF_LINES: usize = 80;
    const CONTEXT: usize = 2;

    let actual_lines: Vec<&str> = actual.lines().collect();
    let expected_lines: Vec<&str> = expected.lines().collect();
    let mut first_diff = 0usize;
    while first_diff < actual_lines.len().min(expected_lines.len())
        && actual_lines[first_diff] == expected_lines[first_diff]
    {
        first_diff += 1;
    }

    let start = first_diff.saturating_sub(CONTEXT);
    let end = (first_diff + CONTEXT + 1)
        .max(actual_lines.len().min(expected_lines.len()))
        .min(start + MAX_DIFF_LINES)
        .min(actual_lines.len().max(expected_lines.len()));

    let mut out = String::new();
    out.push_str(&format!(
        "first differing line: {}\nactual length: {}, expected length: {}\n--- actual\n+++ expected\n",
        first_diff + 1,
        actual.len(),
        expected.len()
    ));
    for idx in start..end {
        let actual_line = actual_lines.get(idx).copied();
        let expected_line = expected_lines.get(idx).copied();
        match (actual_line, expected_line) {
            (Some(a), Some(e)) if a == e => {
                out.push_str(&format!(" {}\n", a));
            }
            (Some(a), Some(e)) => {
                out.push_str(&format!("-{}\n", a));
                out.push_str(&format!("+{}\n", e));
            }
            (Some(a), None) => out.push_str(&format!("-{}\n", a)),
            (None, Some(e)) => out.push_str(&format!("+{}\n", e)),
            (None, None) => {}
        }
    }
    if end < actual_lines.len().max(expected_lines.len()) {
        out.push_str("... diff truncated ...\n");
    }
    out
}

#[test]
fn svg_drawer_golden_has_one_record_per_smiles() {
    let smiles_path = repo_root().join("tests/smiles.smi");
    let expected = std::fs::read_to_string(&smiles_path)
        .unwrap_or_else(|err| panic!("failed to read {}: {err}", smiles_path.display()))
        .lines()
        .filter(|line| {
            let line = line.trim();
            !line.is_empty() && !line.starts_with('#')
        })
        .count();
    let records = load_golden();
    assert_eq!(
        records.len(),
        expected,
        "SVG drawer golden row count must match tests/smiles.smi"
    );
}

#[test]
fn svg_drawer_matches_rdkit_golden() {
    let records = load_golden();
    for (row_idx, record) in records.iter().enumerate() {
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
        let actual = mol
            .to_svg(record.width, record.height)
            .unwrap_or_else(|err| {
                panic!(
                    "cosmolkit failed to draw SVG for row {} ({}): {err}",
                    row_idx + 1,
                    record.smiles
                )
            });
        let expected = record.svg.as_ref().unwrap_or_else(|| {
            panic!(
                "row {} ({}) is rdkit ok but missing SVG",
                row_idx + 1,
                record.smiles
            )
        });
        let actual = normalize_svg_identity(&actual);
        let expected = normalize_svg_identity(expected);
        assert!(
            actual == expected,
            "SVG drawer mismatch at row {} ({})\n{}",
            row_idx + 1,
            record.smiles,
            svg_diff(&actual, &expected)
        );
    }

    let smiles = records
        .iter()
        .map(|record| record.smiles.clone())
        .collect::<Vec<_>>();
    let batch = MoleculeBatch::from_smiles_list(&smiles, BatchErrorMode::Keep)
        .expect("batch SMILES parse should not raise in keep mode");
    let batch_svgs = batch
        .to_svg_list(300, 300)
        .expect("batch SVG rendering should succeed");
    for (row_idx, (record, actual)) in records.iter().zip(batch_svgs).enumerate() {
        if !record.rdkit_ok {
            continue;
        }
        let actual = normalize_svg_identity(
            actual
                .as_deref()
                .unwrap_or_else(|| panic!("batch SVG missing at row {}", row_idx + 1)),
        );
        let expected =
            normalize_svg_identity(record.svg.as_ref().expect("RDKit ok row should have SVG"));
        assert!(
            actual == expected,
            "batch SVG drawer mismatch at row {} ({})\n{}",
            row_idx + 1,
            record.smiles,
            svg_diff(&actual, &expected)
        );
    }
}
