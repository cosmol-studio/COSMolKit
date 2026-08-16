use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use cosmolkit_core::{
    InchiErrorKind, Molecule, cached_valence_assignment, mol_from_inchi, mol_to_inchi,
};
use serde::Deserialize;
use serde_json::json;
use sha2::{Digest, Sha256};

mod common;
use common::parity_data;

#[derive(Debug, Deserialize)]
struct InchiGoldenRecord {
    row: usize,
    smiles: String,
    rdkit_ok: bool,
    inchi: Option<String>,
    error: Option<String>,
    mol_from_inchi_branches: BTreeMap<String, MolFromInchiBranchGolden>,
}

#[derive(Debug, Deserialize)]
struct MolFromInchiBranchGolden {
    status: String,
    digest: Option<String>,
    #[allow(dead_code)]
    error_kind: Option<String>,
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

fn mol_from_inchi_digest(molecule: &Molecule) -> Result<String, String> {
    let valence = cached_valence_assignment(molecule)
        .ok_or_else(|| "MolFromInchi result has no source property-cache valence".to_owned())?;
    let mut degrees = vec![0_usize; molecule.num_atoms()];
    for bond in molecule.bonds() {
        degrees[bond.begin().index()] += 1;
        degrees[bond.end().index()] += 1;
    }
    let atoms = molecule
        .atoms()
        .iter()
        .map(|atom| {
            let index = atom.id().index();
            let explicit = valence.explicit_valence[index];
            let implicit = valence.implicit_hydrogens[index];
            json!([
                atom.atomic_number(),
                atom.formal_charge(),
                atom.chiral_tag().rdkit_name(),
                atom.isotope(),
                atom.is_aromatic(),
                atom.explicit_hydrogens(),
                degrees[index],
                explicit,
                implicit,
                i32::from(atom.explicit_hydrogens()) + implicit,
                explicit + implicit,
                atom.radical_electrons(),
                atom.no_implicit(),
            ])
        })
        .collect::<Vec<_>>();
    let bonds = molecule
        .bonds()
        .iter()
        .map(|bond| {
            json!([
                bond.begin().index(),
                bond.end().index(),
                bond.order().rdkit_name(),
                bond.direction().rdkit_name(),
                bond.stereo().rdkit_name(),
                bond.stereo_atoms()
                    .map(|atoms| atoms.map(|atom| atom.index()).to_vec())
                    .unwrap_or_default(),
                bond.is_aromatic(),
            ])
        })
        .collect::<Vec<_>>();
    let smiles = molecule
        .to_smiles(true)
        .map_err(|error| format!("MolFromInchi result cannot be serialized: {error}"))?;
    let encoded = serde_json::to_vec(&json!([smiles, atoms, bonds]))
        .map_err(|error| format!("MolFromInchi state cannot be encoded: {error}"))?;
    Ok(format!("{:x}", Sha256::digest(encoded)))
}

fn compare_mol_from_inchi_branches(
    row: usize,
    inchi: &str,
    expected: &BTreeMap<String, MolFromInchiBranchGolden>,
    mismatches: &mut Vec<String>,
) {
    for sanitize in [false, true] {
        for remove_hs in [false, true] {
            let branch = format!(
                "sanitize{}_remove_hs{}",
                u8::from(sanitize),
                u8::from(remove_hs)
            );
            let golden = expected
                .get(&branch)
                .unwrap_or_else(|| panic!("row {row}: InChI golden is missing branch {branch}"));
            let actual = match mol_from_inchi(inchi.as_bytes(), sanitize, remove_hs) {
                Ok(output) => match output.molecule {
                    Some(molecule) => match mol_from_inchi_digest(&molecule) {
                        Ok(digest) => ("ok", Some(digest)),
                        Err(error) => {
                            mismatches.push(format!("row {row} branch {branch}: {error}"));
                            continue;
                        }
                    },
                    None => ("none", None),
                },
                Err(error) if error.kind == InchiErrorKind::SanitizeFailed => ("none", None),
                Err(_) => ("error", None),
            };
            if actual.0 != golden.status || actual.1.as_deref() != golden.digest.as_deref() {
                mismatches.push(format!(
                    "row {row} branch {branch}: expected status {} digest {:?}, got status {} digest {:?}",
                    golden.status, golden.digest, actual.0, actual.1
                ));
            }
        }
    }
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
                } else {
                    compare_mol_from_inchi_branches(
                        row,
                        expected_inchi,
                        &expected.mol_from_inchi_branches,
                        &mut mismatches,
                    );
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
