use std::path::Path;
use std::process::Command;

use cosmolkit_test_support::{
    KnownFailure, KnownFailureOutcome, TestDataKind, classify_known_failure, count_smiles_rows,
    expected_family_dir, expected_path, expected_path_for_profile, golden_path,
    load_known_failures, profile_name, rdkit_expected_domain, rdkit_prepare_command,
    regenerate_command, repo_root, smiles_path, testdata_path, validate_expected_family,
};

fn run_profile_probe(profile: Option<&str>, smiles: Option<&str>, expectation: &str) -> bool {
    let mut command = Command::new(std::env::current_exe().expect("test executable should exist"));
    command
        .arg("--exact")
        .arg("profile_probe")
        .arg("--nocapture")
        .env("COSMOLKIT_SUP_PROFILE_PROBE", expectation)
        .env_remove("COSMOLKIT_PARITY_PROFILE")
        .env_remove("COSMOLKIT_PARITY_SMILES");
    if let Some(profile) = profile {
        command.env("COSMOLKIT_PARITY_PROFILE", profile);
    }
    if let Some(smiles) = smiles {
        command.env("COSMOLKIT_PARITY_SMILES", smiles);
    }
    command
        .status()
        .expect("profile probe should launch")
        .success()
}

#[test]
fn profile_probe() {
    let Ok(expectation) = std::env::var("COSMOLKIT_SUP_PROFILE_PROBE") else {
        return;
    };
    match expectation.as_str() {
        "small" => {
            assert_eq!(profile_name(), "smiles_small");
            assert!(smiles_path().ends_with("testdata/smiles/corpus/smiles_small.smi"));
        }
        "strict" => {
            assert_eq!(profile_name(), "smiles_5000");
            assert!(smiles_path().ends_with("testdata/smiles/corpus/smiles_5000.smi"));
        }
        "override" => assert_eq!(smiles_path(), Path::new("isolated-corpus.smi")),
        "invalid" => {
            let _ = profile_name();
            panic!("an invalid profile unexpectedly returned");
        }
        other => panic!("unknown probe expectation {other}"),
    }
}

#[test]
fn profile_aliases_override_and_invalid_value_are_isolated() {
    for alias in [None, Some("small"), Some("smiles_small")] {
        assert!(run_profile_probe(alias, None, "small"));
    }
    for alias in [Some("strict"), Some("5000"), Some("smiles_5000")] {
        assert!(run_profile_probe(alias, None, "strict"));
    }
    assert!(run_profile_probe(
        Some("small"),
        Some("isolated-corpus.smi"),
        "override"
    ));
    assert!(!run_profile_probe(Some("not-a-profile"), None, "invalid"));
}

#[test]
fn repository_paths_counts_and_commands_are_canonical() {
    let root = repo_root();
    assert!(root.join("Cargo.toml").is_file());
    assert!(count_smiles_rows() > 0);
    assert!(expected_family_dir("smiles", "rdkit").starts_with(root.join("testdata/smiles")));
    assert!(rdkit_prepare_command("smiles").contains("--suite smiles --jobs 4"));
    assert_eq!(regenerate_command(), rdkit_prepare_command("all"));
}

#[test]
fn every_declared_rdkit_output_route_has_one_domain() {
    for (domain, files) in [
        ("alignment", &["molalign.jsonl"][..]),
        ("inchi", &["inchi.jsonl"]),
        (
            "molblock",
            &[
                "molblock_v2000_kekulized.jsonl",
                "molblock_v2000_minimal.jsonl",
                "molfile_read.jsonl",
            ],
        ),
        ("sdf", &["sdf_read.jsonl", "sdf_write.jsonl"]),
        (
            "fingerprint",
            &[
                "morgan_fingerprint.jsonl",
                "atom_pair_fingerprint.jsonl",
                "pattern_fingerprint.jsonl",
                "maccs_fingerprint.jsonl",
                "rdkit_topological_fingerprint.jsonl",
                "layered_fingerprint.jsonl",
                "topological_torsion_fingerprint.jsonl",
                "avalon_fingerprint.jsonl",
            ],
        ),
        (
            "conformer",
            &[
                "conformer_generation.jsonl",
                "conformer_generation_library.jsonl",
                "confseq_embed_template.jsonl",
            ],
        ),
        (
            "forcefield",
            &["forcefield_params.jsonl", "mmff_builtin.jsonl"],
        ),
        ("forcefield_coverage", &["forcefield_coverage.jsonl"]),
        ("smiles", &["smiles_writer.jsonl", "isomeric_smiles.jsonl"]),
        (
            "depiction",
            &["svg_drawer.jsonl", "prepared_draw_molecule.jsonl"],
        ),
        ("graph", &["graph_features.jsonl"]),
        ("descriptors", &["molecular_descriptors.jsonl"]),
        (
            "stereo",
            &[
                "tetrahedral_stereo_geometry.jsonl",
                "assign_atom_chiral_tags_from_structure.jsonl",
                "python_stereoisomer_corpus.jsonl",
            ],
        ),
        ("distgeom", &["dg_bounds_matrix.jsonl"]),
        ("mol2", &["mol2_read.jsonl"]),
        ("xyz", &["xyz_read.jsonl"]),
        ("kekulize", &["kekulize_clear_flags_false.jsonl"]),
        (
            "substructure",
            &[
                "delete_substructs.jsonl",
                "delete_substructs_onlyfrags_chirality.jsonl",
            ],
        ),
        ("smarts", &["smarts.jsonl"]),
        ("rdkit_builtin", &["rdkit_builtin_fixture_migration.jsonl"]),
    ] {
        for file in files {
            assert_eq!(rdkit_expected_domain(file), domain, "route for {file}");
        }
    }
    assert!(std::panic::catch_unwind(|| rdkit_expected_domain("unknown.jsonl")).is_err());
}

#[test]
fn fixture_and_corpus_paths_are_confined_and_regular_files() {
    let fixture = testdata_path(
        "molblock",
        TestDataKind::Fixture,
        "molblock_stereo_direction_components.jsonl",
    )
    .expect("real fixture should resolve");
    let corpus = testdata_path("topology", TestDataKind::Corpus, "core.csv")
        .expect("real corpus should resolve");
    assert!(fixture.is_file());
    assert!(corpus.is_file());

    for error in [
        testdata_path("", TestDataKind::Fixture, "x").unwrap_err(),
        testdata_path("../molblock", TestDataKind::Fixture, "x").unwrap_err(),
        testdata_path("molblock", TestDataKind::Fixture, "../README.md").unwrap_err(),
        testdata_path("molblock", TestDataKind::Fixture, "/etc/passwd").unwrap_err(),
        testdata_path("molblock", TestDataKind::Fixture, ".").unwrap_err(),
        testdata_path("molblock", TestDataKind::Fixture, "missing.file").unwrap_err(),
    ] {
        assert!(!error.is_empty());
    }
}

#[test]
fn real_known_failure_is_loaded_losslessly_and_never_skipped() {
    let records = load_known_failures("topology_invariants.jsonl")
        .expect("committed known failures should load");
    assert_eq!(records.len(), 1);
    let record = &records[0];
    assert_eq!(record.case_id, "stereo_tetra_h_001");
    assert_eq!(record.operation.as_deref(), Some("without_hydrogens"));
    assert_eq!(record.invariant.as_deref(), Some("stereo_remapping"));
    assert_eq!(record.expected_failure_kind, "InvalidStereoReference");
    assert_eq!(
        record.expected_error_kind.as_deref(),
        Some("InvalidStereoReference")
    );
    assert!(matches!(
        classify_known_failure(record, Some("InvalidStereoReference")),
        KnownFailureOutcome::ExpectedFailure(found) if found == record
    ));
    assert!(matches!(
        classify_known_failure(record, None),
        KnownFailureOutcome::UnexpectedPass(found) if found == record
    ));
    assert!(matches!(
        classify_known_failure(record, Some("DifferentFailure")),
        KnownFailureOutcome::UnexpectedFailureKind {
            record: found,
            actual_failure_kind: "DifferentFailure"
        } if found == record
    ));
}

#[test]
fn known_failure_schema_normalizes_common_and_domain_kinds_strictly() {
    let common: KnownFailure = serde_json::from_str(
        r#"{"case_id":"a","feature":"smiles","expected_failure_kind":"Mismatch","reason":"reason","created_at":"2026-01-01","rdkit_version":"2026.03.1","branch":"writer"}"#,
    )
    .expect("common schema should parse");
    assert_eq!(common.expected_failure_kind, "Mismatch");
    assert_eq!(common.expected_error_kind, None);

    let domain: KnownFailure = serde_json::from_str(
        r#"{"case_id":"b","expected_error_kind":"InvalidState","reason":"reason","created_at":"2026-01-02","operation":"sanitize","invariant":"indices"}"#,
    )
    .expect("domain schema should parse");
    assert_eq!(domain.expected_failure_kind, "InvalidState");
    assert_eq!(domain.expected_error_kind.as_deref(), Some("InvalidState"));

    for invalid in [
        r#"{"case_id":"","expected_failure_kind":"X","reason":"r","created_at":"d"}"#,
        r#"{"case_id":"a","reason":"r","created_at":"d"}"#,
        r#"{"case_id":"a","expected_failure_kind":"X","expected_error_kind":"Y","reason":"r","created_at":"d"}"#,
        r#"{"case_id":"a","expected_failure_kind":" ","reason":"r","created_at":"d"}"#,
        r#"{"case_id":"a","expected_failure_kind":"X","reason":"","created_at":"d"}"#,
    ] {
        assert!(
            serde_json::from_str::<KnownFailure>(invalid).is_err(),
            "{invalid}"
        );
    }
}

#[test]
fn known_failure_file_names_reject_traversal_and_missing_files() {
    for name in [
        "",
        "../topology_invariants.jsonl",
        "/tmp/failure.jsonl",
        "missing.jsonl",
    ] {
        assert!(load_known_failures(name).is_err(), "{name}");
    }
}

#[test]
fn expected_data_entrypoints_fail_closed_when_data_is_absent_or_invalid() {
    let unique = format!(
        "cosmolkit-sup-invalid-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock should follow the epoch")
            .as_nanos()
    );
    let family = std::env::temp_dir().join(unique).join("smiles_small");
    std::fs::create_dir_all(&family).expect("temporary family should be created");
    std::fs::write(family.join("manifest.json"), b"{}")
        .expect("invalid manifest should be written");
    let error = validate_expected_family(&family, "smiles", "rdkit")
        .expect_err("schema-invalid manifest must be rejected");
    assert!(error.contains("failed to parse"), "{error}");
    std::fs::remove_dir_all(family.parent().expect("family should have a parent"))
        .expect("temporary family should be removed");

    let prepared_golden = golden_path("smiles_writer.jsonl");
    assert!(prepared_golden.is_file());
    assert!(prepared_golden.starts_with(repo_root().join("testdata/smiles/expected/rdkit")));
    assert!(
        std::panic::catch_unwind(|| {
            expected_path("smiles", "rdkit", "missing-sup-output.jsonl")
        })
        .is_err()
    );
    assert!(
        std::panic::catch_unwind(|| {
            expected_path_for_profile(
                "smiles",
                "rdkit",
                "missing-sup-profile",
                "missing-sup-output.jsonl",
            )
        })
        .is_err()
    );
}
