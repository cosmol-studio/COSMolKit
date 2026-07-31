use std::fs::File;
use std::io::{BufRead, BufReader};

use cosmolkit_core::{InchiErrorKind, Molecule, mol_to_inchi};
use serde::Deserialize;

mod common;
use common::parity_data;

#[derive(Debug, Deserialize)]
struct InchiGoldenRecord {
    row: usize,
    smiles: String,
    rdkit_ok: bool,
    inchi: Option<String>,
    error: Option<String>,
}

#[derive(Debug)]
enum ActualInchi {
    Output(Vec<u8>),
    UnsupportedState,
    ParseFailed(String),
    GenerationFailed(String),
}

impl ActualInchi {
    fn description(&self) -> String {
        match self {
            Self::Output(bytes) => format!("output {:?}", String::from_utf8_lossy(bytes)),
            Self::UnsupportedState => "structured UnsupportedState".to_string(),
            Self::ParseFailed(error) => format!("parse error {error}"),
            Self::GenerationFailed(error) => format!("generation error {error}"),
        }
    }
}

fn load_corpus() -> Vec<String> {
    let path = parity_data::smiles_path();
    std::fs::read_to_string(&path)
        .unwrap_or_else(|error| panic!("failed to read {}: {error}", path.display()))
        .lines()
        .filter_map(|raw| {
            let line = raw.trim();
            (!line.is_empty() && !line.starts_with('#')).then(|| line.to_string())
        })
        .collect()
}

fn load_golden() -> Vec<InchiGoldenRecord> {
    let path = parity_data::golden_path("inchi.jsonl");
    let file = File::open(&path).unwrap_or_else(|error| {
        panic!(
            "failed to open {}; regenerate RDKit goldens with `{}`: {error}",
            path.display(),
            parity_data::regenerate_command()
        )
    });
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
fn inchi_matches_pinned_rdkit_for_every_active_profile_row() {
    let corpus = load_corpus();
    let golden = load_golden();
    let mut mismatches = Vec::new();
    assert_eq!(
        golden.len(),
        corpus.len(),
        "InChI golden row count must match the active parity corpus"
    );
    assert_eq!(
        golden.len(),
        parity_data::count_smiles_rows(),
        "InChI test must execute every active parity corpus row"
    );

    for (index, (smiles, expected)) in corpus.iter().zip(&golden).enumerate() {
        let row = index + 1;
        assert_eq!(expected.row, row, "golden row index changed at row {row}");
        assert_eq!(
            expected.smiles, *smiles,
            "golden SMILES differs from corpus row {row}"
        );

        let actual = match Molecule::from_smiles(smiles) {
            Ok(molecule) => match mol_to_inchi(&molecule, None) {
                Ok(output) => ActualInchi::Output(output.inchi),
                Err(error) if error.kind == InchiErrorKind::UnsupportedState => {
                    ActualInchi::UnsupportedState
                }
                Err(error) => ActualInchi::GenerationFailed(error.to_string()),
            },
            Err(error) => ActualInchi::ParseFailed(error.to_string()),
        };
        match (expected.rdkit_ok, expected.inchi.as_deref(), actual) {
            (true, Some(expected_inchi), ActualInchi::Output(actual)) => {
                if actual.as_slice() != expected_inchi.as_bytes() {
                    mismatches.push(format!(
                        "row {row} ({smiles}): expected {expected_inchi}, got {}",
                        String::from_utf8_lossy(&actual)
                    ));
                }
            }
            (true, Some(expected_inchi), actual) => mismatches.push(format!(
                "row {row} ({smiles}): RDKit produced {expected_inchi}, but COSMolKit returned {}",
                actual.description()
            )),
            (false, None, ActualInchi::UnsupportedState) => {
                if expected.error.is_none() {
                    mismatches.push(format!(
                        "row {row} ({smiles}): failed RDKit record has no error"
                    ));
                }
            }
            (false, None, ActualInchi::Output(actual)) if actual.is_empty() => {}
            (false, None, ActualInchi::ParseFailed(_))
                if expected.error.as_deref() == Some("MolFromSmiles returned None") => {}
            (false, None, actual) => mismatches.push(format!(
                "row {row} ({smiles}): RDKit produced no InChI ({:?}), but COSMolKit returned {}",
                expected.error,
                actual.description()
            )),
            _ => panic!("row {row} ({smiles}): inconsistent RDKit golden record: {expected:?}"),
        }
    }

    assert!(
        mismatches.is_empty(),
        "{} InChI parity rows failed:\n{}",
        mismatches.len(),
        mismatches.join("\n")
    );
}
