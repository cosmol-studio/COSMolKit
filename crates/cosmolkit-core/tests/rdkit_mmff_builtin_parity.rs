use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use cosmolkit_core::{MmffMolProperties, Molecule, mmff_has_all_molecule_params};
use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct MmffBuiltinRecord {
    fixture: String,
    line_number: usize,
    row_name: String,
    smiles: String,
    variant: String,
    rdkit_ok: bool,
    num_atoms: Option<usize>,
    has_all: Option<bool>,
    props_ok: bool,
    atom_types: Option<Vec<u8>>,
    error: Option<String>,
}

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn load_golden() -> Vec<MmffBuiltinRecord> {
    let path = repo_root().join("tests/golden/mmff_builtin.jsonl");
    let file = File::open(&path).unwrap_or_else(|err| {
        panic!(
            "failed to open {}; regenerate the MMFF built-in RDKit golden with `.venv/bin/python tests/scripts/gen_rdkit_mmff_builtin_golden.py`: {err}",
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

fn expected_fixture_rows(fixture: &str) -> usize {
    let path = repo_root()
        .join("tests/fixtures/forcefield/mmff/rdkit")
        .join(fixture);
    std::fs::read_to_string(&path)
        .unwrap_or_else(|err| panic!("failed to read {}: {err}", path.display()))
        .lines()
        .filter(|line| {
            let line = line.trim();
            !line.is_empty() && !line.starts_with('#') && line != "SMILES Name"
        })
        .count()
}

#[test]
fn mmff_builtin_golden_covers_rdkit_fixture_rows() {
    let records = load_golden();
    let fixtures = [
        "MMFF94_dative.smi",
        "MMFF94_hypervalent.smi",
        "MMFF94s_dative.smi",
        "MMFF94s_hypervalent.smi",
    ];
    let expected_total = fixtures
        .iter()
        .map(|fixture| expected_fixture_rows(fixture))
        .sum::<usize>();
    assert_eq!(
        records.len(),
        expected_total,
        "MMFF built-in golden row count must match copied RDKit fixture rows"
    );
    for fixture in fixtures {
        let actual = records
            .iter()
            .filter(|record| record.fixture == fixture)
            .count();
        assert_eq!(
            actual,
            expected_fixture_rows(fixture),
            "MMFF built-in golden row count mismatch for {fixture}"
        );
    }
}

#[test]
fn mmff_builtin_has_all_matches_rdkit_golden() {
    let records = load_golden();
    let row_filter = std::env::var("COSMOLKIT_MMFF_BUILTIN_ROW_FILTER")
        .ok()
        .and_then(|s| s.parse::<usize>().ok());

    for (row_idx, record) in records.iter().enumerate() {
        if let Some(filter) = row_filter
            && row_idx + 1 != filter
        {
            continue;
        }

        if !record.rdkit_ok {
            assert!(
                record.error.is_some(),
                "row {} ({}:{} {}) is RDKit not ok but has no error",
                row_idx + 1,
                record.fixture,
                record.line_number,
                record.row_name
            );
            continue;
        }

        let expected = record.has_all.unwrap_or_else(|| {
            panic!(
                "row {} ({}:{} {}) has RDKit ok result without has_all",
                row_idx + 1,
                record.fixture,
                record.line_number,
                record.row_name
            )
        });
        let mol = Molecule::from_smiles(&record.smiles).unwrap_or_else(|err| {
            panic!(
                "COSMolKit failed to parse row {} ({}:{} {} {}): {err}",
                row_idx + 1,
                record.fixture,
                record.line_number,
                record.row_name,
                record.smiles
            )
        });
        let actual = mmff_has_all_molecule_params(&mol).unwrap_or_else(|err| {
            panic!(
                "COSMolKit MMFF parameter coverage errored at row {} ({}:{} {} {}): {err}",
                row_idx + 1,
                record.fixture,
                record.line_number,
                record.row_name,
                record.smiles
            )
        });
        assert_eq!(
            actual,
            expected,
            "MMFF parameter coverage mismatch at row {} ({}:{} {} {})",
            row_idx + 1,
            record.fixture,
            record.line_number,
            record.row_name,
            record.smiles
        );
    }
}

#[test]
fn mmff_builtin_atom_types_match_rdkit_golden() {
    let records = load_golden();
    let row_filter = std::env::var("COSMOLKIT_MMFF_BUILTIN_ROW_FILTER")
        .ok()
        .and_then(|s| s.parse::<usize>().ok());

    for (row_idx, record) in records.iter().enumerate() {
        if let Some(filter) = row_filter
            && row_idx + 1 != filter
        {
            continue;
        }

        if !record.rdkit_ok || !record.props_ok {
            assert!(
                record.error.is_some() || record.atom_types.is_none(),
                "row {} ({}:{} {}) lacks RDKit MMFF properties without error detail",
                row_idx + 1,
                record.fixture,
                record.line_number,
                record.row_name
            );
            continue;
        }

        let expected_atom_types = record.atom_types.as_ref().unwrap_or_else(|| {
            panic!(
                "row {} ({}:{} {}) has RDKit MMFF properties but no atom_types",
                row_idx + 1,
                record.fixture,
                record.line_number,
                record.row_name
            )
        });
        let expected_num_atoms = record.num_atoms.unwrap_or_else(|| {
            panic!(
                "row {} ({}:{} {}) has RDKit MMFF properties but no num_atoms",
                row_idx + 1,
                record.fixture,
                record.line_number,
                record.row_name
            )
        });
        assert_eq!(
            expected_atom_types.len(),
            expected_num_atoms,
            "RDKit atom type count must match RDKit atom count at row {} ({}:{} {})",
            row_idx + 1,
            record.fixture,
            record.line_number,
            record.row_name
        );

        let mol = Molecule::from_smiles(&record.smiles).unwrap_or_else(|err| {
            panic!(
                "COSMolKit failed to parse row {} ({}:{} {} {}): {err}",
                row_idx + 1,
                record.fixture,
                record.line_number,
                record.row_name,
                record.smiles
            )
        });
        let props = MmffMolProperties::new(&mol, &record.variant, 0).unwrap_or_else(|err| {
            panic!(
                "COSMolKit MMFF atom typing errored at row {} ({}:{} {} {} variant {}), expected RDKit atom types {:?}: {err}",
                row_idx + 1,
                record.fixture,
                record.line_number,
                record.row_name,
                record.smiles,
                record.variant,
                expected_atom_types
            )
        });
        let actual_atom_types = (0..mol.num_atoms())
            .map(|idx| {
                props.get_mmff_atom_type(idx).unwrap_or_else(|err| {
                    panic!(
                        "COSMolKit MMFF atom type lookup failed at row {} ({}:{} {} {}) atom {}: {err}",
                        row_idx + 1,
                        record.fixture,
                        record.line_number,
                        record.row_name,
                        record.smiles,
                        idx
                    )
                })
            })
            .collect::<Vec<_>>();

        assert_eq!(
            actual_atom_types,
            *expected_atom_types,
            "MMFF atom type mismatch at row {} ({}:{} {} {} variant {})",
            row_idx + 1,
            record.fixture,
            record.line_number,
            record.row_name,
            record.smiles,
            record.variant
        );
    }
}
