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
        "hydrogens,op-contracts-strict"
    } else {
        "hydrogens"
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
