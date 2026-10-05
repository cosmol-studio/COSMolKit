use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use std::collections::BTreeMap;

use cosmolkit::{EmbedParams as EmbedParameters, Molecule, SdfReadParams};
use serde::Deserialize;

mod parity_data {
    use std::path::PathBuf;
    pub fn repo_root() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
    }
    pub fn golden_path(name: &str) -> PathBuf {
        let directory = std::env::var_os("COSMOLKIT_CONFORMER_REFERENCE_ROOT")
            .map(PathBuf::from)
            .unwrap_or_else(|| repo_root().join("testdata/conformer/expected/rdkit/smiles_small"));
        let path = directory.join(name);
        use sha2::{Digest, Sha256};
        let bytes = std::fs::read(&path).unwrap();
        assert_eq!(
            Sha256::digest(&bytes)
                .iter()
                .map(|byte| format!("{byte:02x}"))
                .collect::<String>(),
            "710255e3ace6d263cbbe47e7e68dd71db2cb5ca8986456c78c042bbd9dc6c352",
            "original fixed19 golden bytes changed"
        );
        path
    }
    pub fn regenerate_command() -> &'static str {
        "use the unchanged original fixed19 input and generator; do not regenerate expectations to fit implementation"
    }
}
fn embed_molecule(
    molecule: &Molecule,
    params: &mut EmbedParameters,
) -> Result<(Molecule, i32), cosmolkit::OperationError> {
    let result = molecule.with_3d_conformer_result_with_params(params)?;
    *params = result.params;
    Ok((result.molecule, result.conf_id))
}
fn embed_multiple_confs(
    molecule: &Molecule,
    count: u32,
    params: &mut EmbedParameters,
) -> Result<(Molecule, Vec<i32>), cosmolkit::OperationError> {
    let result = molecule.with_3d_conformers_result_with_params(count, params)?;
    *params = result.params;
    Ok((result.molecule, result.conf_ids))
}

const COORD_TOLERANCE: f64 = 1.0e-6;

#[derive(Debug, Deserialize)]
struct ConformerGenerationGoldenRecord {
    case_id: String,
    mode: String,
    source_kind: String,
    source: String,
    preset: String,
    rdkit_ok: bool,
    status: Option<i32>,
    ids: Option<Vec<i32>>,
    conformers: Option<Vec<Vec<[f64; 3]>>>,
    failure_counts: Option<Vec<u32>>,
    error: Option<String>,
}

#[derive(Debug, Deserialize)]
struct ConformerGenerationFixtureRecord {
    fixture: Option<String>,
    source: String,
}

fn repo_root() -> PathBuf {
    parity_data::repo_root()
}

fn inventory_path() -> PathBuf {
    repo_root().join("testdata/conformer/fixtures/rdkit_inventory.jsonl")
}

fn fixture_root() -> PathBuf {
    repo_root().join("testdata/conformer/fixtures")
}

fn load_inventory() -> Vec<ConformerGenerationFixtureRecord> {
    let path = inventory_path();
    let file =
        File::open(&path).unwrap_or_else(|err| panic!("failed to open {}: {err}", path.display()));
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

fn vendored_fixture_path(source: &str) -> PathBuf {
    let record = load_inventory()
        .into_iter()
        .find(|record| record.source == source)
        .unwrap_or_else(|| panic!("fixture inventory missing source {source}"));
    let fixture = record.fixture.unwrap_or_else(|| {
        panic!("fixture inventory source {source} does not have a vendored fixture path")
    });
    fixture_root().join(fixture)
}

fn load_golden() -> Vec<ConformerGenerationGoldenRecord> {
    let path = parity_data::golden_path("conformer_generation.jsonl");
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

fn read_rdkit_mol_fixture(source: &str) -> Molecule {
    let path = vendored_fixture_path(source);
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|err| panic!("failed to read RDKit mol fixture {}: {err}", path.display()));
    Molecule::from_sdf_with_params(
        &text,
        &SdfReadParams {
            sanitize: true,
            remove_hydrogens: false,
            process_property_lists: false,
            ..Default::default()
        },
    )
    .unwrap_or_else(|err| {
        panic!(
            "failed to parse RDKit mol fixture {}: {err}",
            path.display()
        )
    })
}

fn load_case_molecule(source_kind: &str, source: &str) -> Molecule {
    match source_kind {
        "fixture_mol" => read_rdkit_mol_fixture(source),
        "smiles" => Molecule::from_smiles(source)
            .unwrap_or_else(|err| panic!("failed to parse SMILES {source}: {err}")),
        "smiles_with_hydrogens" => Molecule::from_smiles(source)
            .unwrap_or_else(|err| panic!("failed to parse SMILES {source}: {err}"))
            .with_hydrogens()
            .unwrap_or_else(|err| panic!("failed to add hydrogens to SMILES {source}: {err}")),
        other => panic!("unsupported source_kind {other}"),
    }
}

fn params_for_case(case_id: &str, preset: &str) -> (EmbedParameters, Option<usize>) {
    let mut params = match preset {
        "DG" => EmbedParameters::default(),
        "KDG" => EmbedParameters::kdg(),
        "ETDG" => EmbedParameters::etdg(),
        "ETDGv2" => EmbedParameters::etdg_v2(),
        "ETKDG" => EmbedParameters::etkdg(),
        "ETKDGv2" => EmbedParameters::etkdg_v2(),
        "ETKDGv3" => EmbedParameters::etkdg_v3(),
        "srETKDGv3" => EmbedParameters::sr_etkdg_v3(),
        other => panic!("unsupported preset {other}"),
    };
    let num_confs = match case_id {
        "multi_embed_fragments_separately" => Some(3),
        "multi_etkdg_seeded" => Some(10),
        "multi_etkdg_pruned_symmetry" => Some(10),
        "multi_etkdg_sequential_seeds" => Some(4),
        "multi_force_trans_amides_true" => Some(10),
        "multi_force_trans_amides_false" => Some(10),
        _ => None,
    };

    match case_id {
        "single_dg_simple_torsion"
        | "single_kdg_simple_torsion"
        | "single_etdg_simple_torsion"
        | "single_etdgv2_torsion"
        | "single_etkdg_simple_torsion"
        | "single_etkdgv2_torsion"
        | "single_etkdgv3_macrocycle"
        | "single_sretkdgv3_smallring" => {
            params.random_seed = 42;
            params.num_threads = 1;
        }
        "single_chirality_failure_fixture" => {
            params.random_seed = 0xf00d;
            params.num_threads = 1;
            params.track_failures = true;
            params.max_iterations = 50;
        }
        "single_coordmap_randomcoords" => {
            params.random_seed = 0xC0FFEE;
            params.num_threads = 1;
            params.use_random_coords = true;
            params.coord_map = Some(BTreeMap::from([
                (0, [0.0, 0.0, 0.0]),
                (1, [0.0, 0.0, 1.5]),
                (2, [0.0, 1.5, 1.5]),
            ]));
        }
        "single_cpci_etkdgv3" => {
            params.random_seed = 0xC0FFEE;
            params.num_threads = 1;
            params.cpci = Some(BTreeMap::from([((0, 3), 0.5), ((1, 4), -0.25)]));
        }
        "single_etkdgv3_x0_ring_connectivity_first"
        | "single_etkdgv3_x0_ring_connectivity_second" => {
            params.max_iterations = 3;
            params.random_seed = 61453;
            params.num_threads = 1;
            params.timeout = 0;
        }
        "multi_embed_fragments_separately" => {
            params.random_seed = 0xf00d;
            params.num_threads = 1;
            params.embed_fragments_separately = true;
        }
        "multi_etkdg_seeded" => {
            params.random_seed = 0xf00d;
            params.num_threads = 1;
            params.track_failures = true;
            params.timeout = 1;
        }
        "multi_etkdg_pruned_symmetry" => {
            params.random_seed = 0xf00d;
            params.num_threads = 1;
            params.prune_rms_thresh = 0.5;
            params.use_symmetry_for_pruning = true;
        }
        "multi_etkdg_sequential_seeds" => {
            params.random_seed = 0xf00d;
            params.num_threads = 1;
            params.enable_sequential_random_seeds = true;
        }
        "multi_force_trans_amides_true" => {
            params.random_seed = 0xf00d;
            params.num_threads = 1;
            params.force_trans_amides = true;
            params.use_exp_torsion_angle_prefs = false;
            params.use_basic_knowledge = true;
        }
        "multi_force_trans_amides_false" => {
            params.random_seed = 0xf00d;
            params.num_threads = 1;
            params.force_trans_amides = false;
            params.use_exp_torsion_angle_prefs = false;
            params.use_basic_knowledge = true;
        }
        other => panic!("unsupported case_id {other}"),
    }

    (params, num_confs)
}

fn assert_conformer_coords_match(
    case_id: &str,
    actual: &[Vec<[f64; 3]>],
    expected: &[Vec<[f64; 3]>],
) {
    assert_eq!(
        actual.len(),
        expected.len(),
        "{case_id}: conformer count mismatch"
    );
    for (conf_idx, (actual_conf, expected_conf)) in actual.iter().zip(expected).enumerate() {
        assert_eq!(
            actual_conf.len(),
            expected_conf.len(),
            "{case_id}: atom count mismatch in conformer {conf_idx}"
        );
        for (atom_idx, (actual_xyz, expected_xyz)) in
            actual_conf.iter().zip(expected_conf).enumerate()
        {
            for axis in 0..3 {
                let a = actual_xyz[axis];
                let e = expected_xyz[axis];
                assert!(
                    (a - e).abs() <= COORD_TOLERANCE,
                    "{case_id}: conformer {conf_idx} atom {atom_idx} axis {axis} mismatch: ours={a} expected={e}"
                );
            }
        }
    }
}

fn conformer_coords(mol: &Molecule) -> Vec<Vec<[f64; 3]>> {
    mol.conformers_3d()
        .iter()
        .map(|conformer| conformer.coordinates().to_vec())
        .collect()
}

#[test]
fn conformer_generation_golden_has_expected_case_count() {
    let records = load_golden();
    assert_eq!(
        records.len(),
        19,
        "unexpected conformer-generation golden case count"
    );
}

#[test]
fn conformer_generation_exhaustive_parity_matches_exact_coordinates() {
    let records = load_golden();
    for record in records {
        assert!(
            record.rdkit_ok,
            "{}: RDKit golden generation failed: {:?}",
            record.case_id, record.error
        );

        let expected_conformers = record
            .conformers
            .as_ref()
            .unwrap_or_else(|| panic!("{}: RDKit golden has no conformer payload", record.case_id));
        let expected_failures = record.failure_counts.as_deref().unwrap_or(&[]);

        let mol = load_case_molecule(&record.source_kind, &record.source);
        let (mut params, num_confs) = params_for_case(&record.case_id, &record.preset);

        match record.mode.as_str() {
            "single" => {
                let (embedded, status) = embed_molecule(&mol, &mut params).unwrap_or_else(|err| {
                    panic!(
                        "{}: COSMolKit single embedding failed: {err}",
                        record.case_id
                    )
                });
                assert_eq!(
                    Some(status),
                    record.status,
                    "{}: status mismatch",
                    record.case_id
                );
                assert_eq!(
                    params.failures.as_slice(),
                    expected_failures,
                    "{}: failure counts mismatch",
                    record.case_id
                );
                let actual_conformers = conformer_coords(&embedded);
                assert_conformer_coords_match(
                    &record.case_id,
                    &actual_conformers,
                    expected_conformers,
                );
            }
            "multi" => {
                let (embedded, ids) = embed_multiple_confs(
                    &mol,
                    u32::try_from(num_confs.expect("multi case has count"))
                        .expect("multi conformer count fits u32"),
                    &mut params,
                )
                .unwrap_or_else(|err| {
                    panic!(
                        "{}: COSMolKit multi embedding failed: {err}",
                        record.case_id
                    )
                });
                assert_eq!(
                    Some(ids),
                    record.ids,
                    "{}: conformer id mismatch",
                    record.case_id
                );
                assert_eq!(
                    params.failures.as_slice(),
                    expected_failures,
                    "{}: failure counts mismatch",
                    record.case_id
                );
                let actual_conformers = conformer_coords(&embedded);
                assert_conformer_coords_match(
                    &record.case_id,
                    &actual_conformers,
                    expected_conformers,
                );
            }
            other => panic!("{}: unsupported mode {other}", record.case_id),
        }
    }
}

fn diagnose_one_original_record(record: ConformerGenerationGoldenRecord) {
    assert!(
        record.rdkit_ok,
        "{}: RDKit golden generation failed: {:?}",
        record.case_id, record.error
    );

    let expected_conformers = record
        .conformers
        .as_ref()
        .unwrap_or_else(|| panic!("{}: RDKit golden has no conformer payload", record.case_id));
    let expected_failures = record.failure_counts.as_deref().unwrap_or(&[]);

    let mol = load_case_molecule(&record.source_kind, &record.source);
    let (mut params, num_confs) = params_for_case(&record.case_id, &record.preset);

    match record.mode.as_str() {
        "single" => {
            let (embedded, status) = embed_molecule(&mol, &mut params).unwrap_or_else(|err| {
                panic!(
                    "{}: COSMolKit single embedding failed: {err}",
                    record.case_id
                )
            });
            assert_eq!(
                Some(status),
                record.status,
                "{}: status mismatch",
                record.case_id
            );
            assert_eq!(
                params.failures.as_slice(),
                expected_failures,
                "{}: failure counts mismatch",
                record.case_id
            );
            let actual_conformers = conformer_coords(&embedded);
            assert_conformer_coords_match(&record.case_id, &actual_conformers, expected_conformers);
        }
        "multi" => {
            let (embedded, ids) = embed_multiple_confs(
                &mol,
                u32::try_from(num_confs.expect("multi case has count"))
                    .expect("multi conformer count fits u32"),
                &mut params,
            )
            .unwrap_or_else(|err| {
                panic!(
                    "{}: COSMolKit multi embedding failed: {err}",
                    record.case_id
                )
            });
            assert_eq!(
                Some(ids),
                record.ids,
                "{}: conformer id mismatch",
                record.case_id
            );
            assert_eq!(
                params.failures.as_slice(),
                expected_failures,
                "{}: failure counts mismatch",
                record.case_id
            );
            let actual_conformers = conformer_coords(&embedded);
            assert_conformer_coords_match(&record.case_id, &actual_conformers, expected_conformers);
        }
        other => panic!("{}: unsupported mode {other}", record.case_id),
    }
}

#[test]
fn diagnostic_all_original19_retains_each_original_assertion() {
    let records = load_golden();
    assert_eq!(records.len(), 19);
    let mut failures = Vec::new();
    let mut matched = 0;
    for record in records {
        let case_id = record.case_id.clone();
        match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            diagnose_one_original_record(record)
        })) {
            Ok(()) => {
                eprintln!("ORIGINAL19_MATCH {case_id}");
                matched += 1;
            }
            Err(cause) => {
                let text = cause
                    .downcast_ref::<String>()
                    .cloned()
                    .or_else(|| cause.downcast_ref::<&str>().map(|s| s.to_string()))
                    .unwrap_or_else(|| "non-string source failure".to_string());
                eprintln!("ORIGINAL19_FAILURE {case_id}: {text}");
                failures.push((case_id, text));
            }
        }
    }
    eprintln!(
        "ORIGINAL19_ALL actual_rows=19 matched={matched} failed={}",
        failures.len()
    );
    assert!(failures.is_empty(), "all original19 failures: {failures:?}");
}
