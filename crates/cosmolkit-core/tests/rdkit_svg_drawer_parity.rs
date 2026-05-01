use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use cosmolkit_core::Molecule;
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
        assert_eq!(
            normalize_svg_identity(&actual),
            normalize_svg_identity(expected),
            "SVG drawer mismatch at row {} ({})",
            row_idx + 1,
            record.smiles
        );
    }
}
