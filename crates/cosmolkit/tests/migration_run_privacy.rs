//! Real-module privacy probes share compilations, not just dependencies.
use std::path::PathBuf;
use std::process::{Command, Output};

type Case = (&'static str, &'static [&'static str]);

const ALLOWED: &[&str] = &[
    "allowed",
    "cip_cow_allowed",
    "atom_code_cow_allowed",
    "pending_allowed",
    "preserve_allowed",
    "sanitize_allowed",
];

const ACCESS: &[Case] = &[
    (
        "cip_cow_missing_writes",
        &["E0599", "stage_topology_properties_cow"],
    ),
    (
        "cip_cow_runtime_private",
        &["private", "stage_topology_properties_cow_runtime"],
    ),
    (
        "atom_code_cow_missing_writes",
        &["E0599", "stage_topology_properties"],
    ),
    (
        "atom_code_cow_runtime_private",
        &["private", "stage_topology_properties"],
    ),
    (
        "preserve_forbidden",
        &[
            "checkout_derived_cache",
            "install_derived_cache",
            "clear_cache",
            "mark_cache_updated",
            "read_derived_cache_runtime",
            "checkout_derived_cache_runtime",
            "prove_preserved_runtime",
            "private",
        ],
    ),
    ("preserve_no_read", &["E0599", "derived_cache"]),
    (
        "pending_forbidden",
        &[
            "pending_molecule",
            "pending_molecule_runtime",
            "ensure_unsealed_runtime",
            "finish_result",
            "identity",
            "topology",
            "coordinates",
            "properties",
            "derived_cache",
            "clone",
            "num_atoms",
            "Molecule",
            "private",
        ],
    ),
    (
        "pending_wrong_marker",
        &["mismatched types", "WithHydrogensAccess"],
    ),
    (
        "sanitize_forbidden",
        &[
            "coordinates",
            "read_topology_runtime",
            "read_properties_runtime",
            "checkout_topology_runtime",
            "install_topology_runtime",
            "spec",
            "source",
            "derived_cache",
            "in_place_target",
            "new_in_place",
            "finish",
            "abort_in_place",
            "finish_in_place",
        ],
    ),
    (
        "forbidden",
        &[
            "properties",
            "read_properties_runtime",
            "checkout_topology_runtime",
            "checkout_coordinates_runtime",
            "checkout_properties_runtime",
            "checkout_derived_cache_runtime",
            "install_topology_runtime",
            "install_coordinates_runtime",
            "install_properties_runtime",
            "install_derived_cache_runtime",
            "source_topology_runtime",
            "source_coordinates_runtime",
            "source_properties_runtime",
            "emit_all_runtime",
            "spec",
            "source",
            "topology",
            "coordinates",
            "derived_cache",
            "in_place_target",
            "new",
            "new_in_place",
            "finish",
            "abort_in_place",
            "finish_in_place",
        ],
    ),
];

// Type/access errors stop rustc before borrow checking. Keep these in a
// separate compilation so the lifetime, mutability and move probes run too.
const BORROWS: &[Case] = &[
    (
        "cip_cow_borrow_escape",
        &["lifetime may not live long enough"],
    ),
    ("atom_code_cow_borrow_escape", &["lifetime"]),
    ("preserve_immutable", &["E0596", "mutable"]),
    ("pending_reuse", &["use of moved value", "pending"]),
];

// Struct-field privacy is checked after type/borrow checking succeeds.
const FINALIZER: &[Case] = &[("pending_finalizer", &["private", "parts", "operation"])];

#[test]
fn real_operation_module_capability_boundaries_default_and_strict() {
    for strict in [false, true] {
        let allowed = compile_probe(ALLOWED, strict);
        assert!(allowed.status.success(), "{}", output_text(&allowed));
        for cases in [ACCESS, BORROWS, FINALIZER] {
            let selected: Vec<_> = cases.iter().map(|(case, _)| *case).collect();
            let output = compile_probe(&selected, strict);
            assert!(
                !output.status.success(),
                "{selected:?} unexpectedly compiled"
            );
            let messages: Vec<serde_json::Value> = String::from_utf8_lossy(&output.stdout)
                .lines()
                .map(|line| serde_json::from_str(line).expect("Cargo JSON diagnostic"))
                .filter(|event: &serde_json::Value| {
                    event["reason"] == "compiler-message" && event["message"]["level"] == "error"
                })
                .map(|event| event["message"].clone())
                .collect();
            for &(case, expected) in cases {
                // Attribute every error to its own probe's primary source span.
                // An unrelated failure must not satisfy another probe.
                let errors = messages
                    .iter()
                    .filter(|message| error_belongs_to(message, case))
                    .map(|message| message["rendered"].as_str().expect("rendered error"))
                    .collect::<Vec<_>>()
                    .join("\n");
                assert!(
                    !errors.is_empty(),
                    "strict={strict} case={case}: no probe error\n{}",
                    output_text(&output)
                );
                for diagnostic in expected {
                    assert!(
                        errors.contains(diagnostic),
                        "strict={strict} case={case}: missing {diagnostic:?}\n{errors}"
                    );
                }
            }
        }
    }
    println!(
        "8 compiler invocations; 15 rejected probe groups verified in both default and strict modes"
    );
}

fn error_belongs_to(message: &serde_json::Value, case: &str) -> bool {
    let source = include_str!("../src/ops/runtime_privacy_probe.rs");
    message["spans"]
        .as_array()
        .expect("diagnostic spans")
        .iter()
        .any(|span| {
            span["is_primary"] == true
                && span["file_name"]
                    .as_str()
                    .is_some_and(|path| path.ends_with("/runtime_privacy_probe.rs"))
                && span["line_start"].as_u64().is_some_and(|line| {
                    source
                        .lines()
                        .take(line as usize)
                        .collect::<Vec<_>>()
                        .into_iter()
                        .rev()
                        .find_map(|line| {
                            line.strip_prefix("#[cfg(cosmolkit_runtime_privacy_case = \"")
                                .and_then(|line| line.strip_suffix("\")]"))
                        })
                        == Some(case)
                })
        })
}

fn output_text(output: &Output) -> String {
    let errors = String::from_utf8_lossy(&output.stdout)
        .lines()
        .filter_map(|line| serde_json::from_str::<serde_json::Value>(line).ok())
        .filter(|event| {
            event["reason"] == "compiler-message" && event["message"]["level"] == "error"
        })
        .filter_map(|event| event["message"]["rendered"].as_str().map(str::to_owned))
        .collect::<Vec<_>>()
        .join("\n");
    format!("{errors}\n{}", String::from_utf8_lossy(&output.stderr))
}

fn compile_probe(cases: &[&str], strict: bool) -> Output {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let workspace = manifest_dir
        .parent()
        .and_then(|path| path.parent())
        .unwrap();
    let features = if strict {
        "cap-hydrogens,cap-stereo,cap-stereoisomers,cap-fingerprints,cap-sanitize,op-contracts-strict"
    } else {
        "cap-hydrogens,cap-stereo,cap-stereoisomers,cap-fingerprints,cap-sanitize"
    };
    let mut command = Command::new(std::env::var("CARGO").unwrap_or_else(|_| "cargo".into()));
    command
        .current_dir(workspace)
        .env(
            "CARGO_TARGET_DIR",
            workspace.join("target/runtime-privacy-compile"),
        )
        .env("CARGO_INCREMENTAL", "0")
        .args([
            "rustc",
            "--quiet",
            "--locked",
            "--message-format=json",
            "-p",
            "cosmolkit",
            "--lib",
            "--no-default-features",
            "--features",
            features,
            "--",
            "--emit=metadata",
            "--cfg",
            "cosmolkit_runtime_privacy_probe",
        ]);
    for case in cases {
        command.args([
            "--cfg",
            &format!("cosmolkit_runtime_privacy_case=\"{case}\""),
        ]);
    }
    command.output().expect("run real-module privacy probes")
}

// This integration target is an external crate. Unlike the sibling runtime
// probes above, it cannot use pub(crate) generated Molecule methods.
#[cfg(cosmolkit_external_tautomer_privacy = "allowed")]
fn public_tautomer_wrappers_compile() {
    let _: fn(
        &cosmolkit::Molecule,
    ) -> Result<cosmolkit::TautomerEnumeration, cosmolkit::OperationError> =
        cosmolkit::Molecule::enumerate_tautomers;
    let _ = cosmolkit::Molecule::enumerate_tautomers_with_params;
    let _: fn(&cosmolkit::Molecule) -> Result<cosmolkit::Molecule, cosmolkit::OperationError> =
        cosmolkit::Molecule::canonical_tautomer;
    let _ = cosmolkit::Molecule::canonical_tautomer_with_params;
    let _: fn(&cosmolkit::Molecule) -> Result<cosmolkit::TautomerScore, cosmolkit::OperationError> =
        cosmolkit::Molecule::tautomer_score;
    let _ = cosmolkit::Molecule::tautomer_score_with_params;
    let _ = cosmolkit::TautomerEnumeration::canonical_tautomer;
    let _ = cosmolkit::TautomerEnumeration::canonical_tautomer_with_params;
}

#[cfg(cosmolkit_external_tautomer_privacy = "forbidden")]
fn private_tautomer_generated_methods_do_not_compile() {
    let _ = cosmolkit::Molecule::assign_symm_sssr_;
    let _ = cosmolkit::Molecule::install_tautomer_score_cache_;
    let _ = cosmolkit::Molecule::with_assigned_symm_sssr;
    let _ = cosmolkit::Molecule::with_installed_tautomer_score_cache;
}

#[cfg(feature = "cap-tautomer")]
#[test]
fn private_tautomer_generated_methods_remain_private_default_and_strict() {
    const METHODS: [&str; 4] = [
        "assign_symm_sssr_",
        "install_tautomer_score_cache_",
        "with_assigned_symm_sssr",
        "with_installed_tautomer_score_cache",
    ];
    for strict in [false, true] {
        let allowed = compile_external_tautomer_probe(strict, "allowed");
        assert!(allowed.status.success(), "{}", output_text(&allowed));
        let forbidden = compile_external_tautomer_probe(strict, "forbidden");
        assert!(
            !forbidden.status.success(),
            "private methods unexpectedly compiled"
        );
        let messages: Vec<serde_json::Value> = String::from_utf8_lossy(&forbidden.stdout)
            .lines()
            .map(|line| serde_json::from_str(line).expect("Cargo JSON diagnostic"))
            .filter(|event: &serde_json::Value| {
                event["reason"] == "compiler-message" && event["message"]["level"] == "error"
            })
            .map(|event| event["message"].clone())
            .collect();
        for method in METHODS {
            let expected_line = format!("let _ = cosmolkit::Molecule::{method};");
            let errors = messages
                .iter()
                .filter(|message| {
                    message["spans"]
                        .as_array()
                        .expect("diagnostic spans")
                        .iter()
                        .any(|span| {
                            span["is_primary"] == true
                                && span["file_name"]
                                    .as_str()
                                    .is_some_and(|path| path.ends_with("/migration_run_privacy.rs"))
                                && span["line_start"].as_u64().is_some_and(|line| {
                                    include_str!("migration_run_privacy.rs")
                                        .lines()
                                        .nth(line as usize - 1)
                                        .is_some_and(|source| source.trim() == expected_line)
                                })
                        })
                })
                .collect::<Vec<_>>();
            assert!(
                !errors.is_empty(),
                "strict={strict} method={method}: no attributed error\n{}",
                output_text(&forbidden)
            );
            assert!(
                errors.iter().any(|message| {
                    message["code"]["code"] == "E0624"
                        && message["rendered"]
                            .as_str()
                            .is_some_and(|text| text.contains(method) && text.contains("private"))
                }),
                "strict={strict} method={method}: expected private-method E0624\n{}",
                output_text(&forbidden)
            );
        }
    }
    println!(
        "4 external compiler invocations; all 4 private tautomer methods rejected in default and strict modes; all 8 public wrappers compile"
    );
}

#[cfg(feature = "cap-tautomer")]
fn compile_external_tautomer_probe(strict: bool, case: &str) -> Output {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let workspace = manifest_dir
        .parent()
        .and_then(|path| path.parent())
        .unwrap();
    let features = if strict {
        "cap-tautomer,cap-hydrogens,cap-stereo,cap-stereoisomers,cap-fingerprints,cap-sanitize,op-contracts-strict"
    } else {
        "cap-tautomer,cap-hydrogens,cap-stereo,cap-stereoisomers,cap-fingerprints,cap-sanitize"
    };
    Command::new(std::env::var("CARGO").unwrap_or_else(|_| "cargo".into()))
        .current_dir(workspace)
        .env(
            "CARGO_TARGET_DIR",
            workspace.join("target/runtime-privacy-compile"),
        )
        .env("CARGO_INCREMENTAL", "0")
        .args([
            "rustc",
            "--quiet",
            "--locked",
            "--message-format=json",
            "-p",
            "cosmolkit",
            "--test",
            "migration_run_privacy",
            "--no-default-features",
            "--features",
            features,
            "--",
            "--emit=metadata",
            "--check-cfg",
            "cfg(cosmolkit_external_tautomer_privacy, values(\"allowed\", \"forbidden\"))",
            "--cfg",
            &format!("cosmolkit_external_tautomer_privacy=\"{case}\""),
        ])
        .output()
        .expect("run external tautomer privacy probes")
}
