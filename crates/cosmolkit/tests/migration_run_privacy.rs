use std::path::PathBuf;
use std::process::{Command, Output};

fn cargo_check(case: &str, strict: bool) -> Output {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let workspace = manifest_dir
        .parent()
        .and_then(|path| path.parent())
        .expect("cosmolkit crate must be nested under the workspace crates directory");
    let target = workspace.join("target/runtime-privacy-compile");
    let inherited = std::env::var("RUSTFLAGS").unwrap_or_default();
    let rustflags = format!(
        "{inherited} --cfg cosmolkit_runtime_privacy_probe \
         --cfg=cosmolkit_runtime_privacy_case=\"{case}\""
    );
    let features = if strict {
        "hydrogens,stereo,op-contracts-strict"
    } else {
        "hydrogens,stereo"
    };

    Command::new(std::env::var("CARGO").unwrap_or_else(|_| "cargo".into()))
        .current_dir(workspace)
        .env("CARGO_TARGET_DIR", target)
        .env("CARGO_INCREMENTAL", "0")
        .env("RUSTFLAGS", rustflags)
        .args([
            "check",
            "--quiet",
            "-p",
            "cosmolkit",
            "--lib",
            "--features",
            features,
        ])
        .output()
        .expect("run the real cosmolkit compile-privacy probe")
}

#[test]
fn pending_results_are_registry_scoped_and_finalization_is_wrapper_private() {
    for strict in [false, true] {
        let allowed = cargo_check("pending_allowed", strict);
        assert!(
            allowed.status.success(),
            "{}",
            String::from_utf8_lossy(&allowed.stderr)
        );
        let forbidden = cargo_check("pending_forbidden", strict);
        assert!(!forbidden.status.success());
        let errors = String::from_utf8_lossy(&forbidden.stderr);
        for surface in [
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
        ] {
            assert!(
                errors.contains(surface),
                "missing rejection of {surface}: {errors}"
            );
        }
        let wrong_marker = cargo_check("pending_wrong_marker", strict);
        assert!(!wrong_marker.status.success());
        let errors = String::from_utf8_lossy(&wrong_marker.stderr);
        assert!(
            errors.contains("mismatched types") && errors.contains("WithHydrogensAccess"),
            "{errors}"
        );
        let finalizer = cargo_check("pending_finalizer", strict);
        assert!(!finalizer.status.success());
        let errors = String::from_utf8_lossy(&finalizer.stderr);
        assert!(
            errors.contains("private") && errors.contains("parts") && errors.contains("operation"),
            "{errors}"
        );
        let reuse = cargo_check("pending_reuse", strict);
        assert!(!reuse.status.success());
        let errors = String::from_utf8_lossy(&reuse.stderr);
        assert!(
            errors.contains("use of moved value") && errors.contains("pending"),
            "{errors}"
        );
    }
}

#[test]
fn real_operation_body_module_only_sees_generated_capabilities() {
    for strict in [false, true] {
        let mode = if strict { "strict" } else { "default" };

        let allowed = cargo_check("allowed", strict);
        assert!(
            allowed.status.success(),
            "authorized generated methods did not compile in {mode} mode:\n{}",
            String::from_utf8_lossy(&allowed.stderr)
        );

        let forbidden = cargo_check("forbidden", strict);
        assert!(
            !forbidden.status.success(),
            "runtime internals unexpectedly compiled in {mode} mode"
        );
        let stderr = String::from_utf8_lossy(&forbidden.stderr);
        for rejected_surface in [
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
        ] {
            assert!(
                stderr.contains(rejected_surface),
                "{mode} compile failure did not prove `{rejected_surface}` is hidden:\n{stderr}"
            );
        }
    }
}
