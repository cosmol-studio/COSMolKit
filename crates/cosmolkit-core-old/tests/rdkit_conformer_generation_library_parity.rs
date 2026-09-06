use std::fs::File;
use std::io::{BufRead, BufReader};

use cosmolkit_core::{EmbedParameters, Molecule, chemistry::distgeom::embed_molecule};
use serde::Deserialize;

mod common;
use common::parity_data;

const COORD_TOLERANCE: f64 = 1.0e-6;
const CONFORMER_LIBRARY_SEED: i32 = 61453;
// Keep library parity on a deterministic embed-attempt budget. RDKit's
// timeout is wall-clock based, so it is unsuitable for CI-stable parity rows.
const CONFORMER_LIBRARY_MAX_ITERATIONS: u32 = 3;
const CONFORMER_LIBRARY_TIMEOUT: u32 = 0;

#[derive(Debug, Deserialize)]
struct ConformerGenerationLibraryRecord {
    smiles: String,
    seed: i32,
    preset: String,
    max_iterations: u32,
    timeout: u32,
    rdkit_parse_ok: bool,
    rdkit_add_hs_ok: bool,
    rdkit_embed_ok: bool,
    status: Option<i32>,
    coords: Option<Vec<[f64; 3]>>,
    error_stage: Option<String>,
    error: Option<String>,
}

fn load_golden() -> Vec<ConformerGenerationLibraryRecord> {
    let path = parity_data::golden_path("conformer_generation_library.jsonl");
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

#[test]
fn conformer_generation_library_golden_has_one_record_per_smiles() {
    let expected = parity_data::count_smiles_rows();
    let records = load_golden();
    assert_eq!(
        records.len(),
        expected,
        "conformer-generation library golden row count must match the active parity corpus"
    );
}

#[test]
fn conformer_generation_library_matches_rdkit_golden_under_fixed_iteration_budget() {
    let records = load_golden();

    for (row_idx, record) in records.iter().enumerate() {
        assert_eq!(
            record.seed,
            CONFORMER_LIBRARY_SEED,
            "row {} ({}) unexpected library seed",
            row_idx + 1,
            record.smiles
        );
        assert_eq!(
            record.preset,
            "ETKDGv3",
            "row {} ({}) unexpected library preset",
            row_idx + 1,
            record.smiles
        );
        assert_eq!(
            record.max_iterations,
            CONFORMER_LIBRARY_MAX_ITERATIONS,
            "row {} ({}) unexpected library max_iterations",
            row_idx + 1,
            record.smiles
        );
        assert_eq!(
            record.timeout,
            CONFORMER_LIBRARY_TIMEOUT,
            "row {} ({}) unexpected library timeout",
            row_idx + 1,
            record.smiles
        );

        let parsed = Molecule::from_smiles(&record.smiles);
        if !record.rdkit_parse_ok {
            assert!(
                parsed.is_err(),
                "row {} ({}) RDKit parse failed but COSMolKit parsed successfully",
                row_idx + 1,
                record.smiles
            );
            assert_eq!(
                record.error_stage.as_deref(),
                Some("parse"),
                "row {} ({}) RDKit parse-failure row must be tagged with parse stage",
                row_idx + 1,
                record.smiles
            );
            assert!(
                record.error.is_some(),
                "row {} ({}) RDKit parse-failure row has no error detail",
                row_idx + 1,
                record.smiles
            );
            continue;
        }

        let parsed = parsed.unwrap_or_else(|err| {
            panic!(
                "row {} ({}) COSMolKit parse failed where RDKit parse succeeded: {err}",
                row_idx + 1,
                record.smiles
            )
        });
        let with_h = parsed.with_hydrogens().unwrap_or_else(|err| {
            panic!(
                "row {} ({}) COSMolKit AddHs failed where RDKit AddHs succeeded: {err}",
                row_idx + 1,
                record.smiles
            )
        });

        assert!(
            record.rdkit_add_hs_ok,
            "row {} ({}) RDKit parse succeeded but AddHs did not",
            row_idx + 1,
            record.smiles
        );

        let mut params = EmbedParameters::etkdg_v3();
        params.max_iterations = CONFORMER_LIBRARY_MAX_ITERATIONS;
        params.random_seed = CONFORMER_LIBRARY_SEED;
        params.num_threads = 1;
        params.timeout = CONFORMER_LIBRARY_TIMEOUT;

        let embedded = embed_molecule(&with_h, &mut params);
        match (record.rdkit_embed_ok, embedded) {
            (true, Ok((embedded, status))) => {
                let expected_status = record.status.unwrap_or_else(|| {
                    panic!(
                        "row {} ({}) RDKit success row missing status",
                        row_idx + 1,
                        record.smiles
                    )
                });
                assert_eq!(
                    status,
                    expected_status,
                    "row {} ({}) embed status mismatch",
                    row_idx + 1,
                    record.smiles
                );
                let expected_coords = record.coords.as_ref().unwrap_or_else(|| {
                    panic!(
                        "row {} ({}) RDKit success row missing coords",
                        row_idx + 1,
                        record.smiles
                    )
                });
                let actual_conf = embedded.conformers_3d().first().unwrap_or_else(|| {
                    panic!(
                        "row {} ({}) COSMolKit succeeded but produced no conformer",
                        row_idx + 1,
                        record.smiles
                    )
                });
                assert_eq!(
                    actual_conf.coordinates().len(),
                    expected_coords.len(),
                    "row {} ({}) atom count mismatch after AddHs/embed",
                    row_idx + 1,
                    record.smiles
                );
                for (atom_idx, (actual_xyz, expected_xyz)) in actual_conf
                    .coordinates()
                    .iter()
                    .zip(expected_coords)
                    .enumerate()
                {
                    for axis in 0..3 {
                        let a = actual_xyz[axis];
                        let e = expected_xyz[axis];
                        assert!(
                            (a - e).abs() <= COORD_TOLERANCE,
                            "row {} ({}) atom {} axis {} mismatch: ours={} expected={}",
                            row_idx + 1,
                            record.smiles,
                            atom_idx,
                            axis,
                            a,
                            e
                        );
                    }
                }
            }
            (true, Err(err)) => {
                panic!(
                    "row {} ({}) COSMolKit embed failed where RDKit succeeded: {err}",
                    row_idx + 1,
                    record.smiles
                );
            }
            (false, Ok((_, status))) => {
                let expected_status = record.status.unwrap_or_else(|| {
                    panic!(
                        "row {} ({}) RDKit failed embed row missing status/error detail",
                        row_idx + 1,
                        record.smiles
                    )
                });
                assert_eq!(
                    status,
                    expected_status,
                    "row {} ({}) expected non-success embed status to match RDKit",
                    row_idx + 1,
                    record.smiles
                );
            }
            (false, Err(_)) => {
                assert_eq!(
                    record.error_stage.as_deref(),
                    Some("embed"),
                    "row {} ({}) RDKit failed embed row must be tagged with embed stage",
                    row_idx + 1,
                    record.smiles
                );
                assert!(
                    record.error.is_some(),
                    "row {} ({}) RDKit failed embed row has no error detail",
                    row_idx + 1,
                    record.smiles
                );
            }
        }
    }
}
