//! Proposed extra compile checks on the real body/runtime layout; old cases retained.
use std::path::PathBuf;
use std::process::{Command, Output};

fn cargo_check(case: &str, strict: bool) -> Output {
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let workspace = manifest.parent().unwrap().parent().unwrap();
    let inherited = std::env::var("RUSTFLAGS").unwrap_or_default();
    let flags = format!(
        "{inherited} --cfg cosmolkit_runtime_privacy_probe --cfg=cosmolkit_runtime_privacy_case=\"{case}\""
    );
    let features = if strict {
        "cap-fingerprints,cap-hydrogens,cap-stereo,op-contracts-strict"
    } else {
        "cap-fingerprints,cap-hydrogens,cap-stereo"
    };
    Command::new(std::env::var("CARGO").unwrap_or_else(|_| "cargo".into()))
        .current_dir(workspace)
        .env(
            "CARGO_TARGET_DIR",
            workspace.join("target/runtime-atom-code-cow-privacy-compile"),
        )
        .env("CARGO_INCREMENTAL", "0")
        .env("RUSTFLAGS", flags)
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
        .expect("run actual real-module scoped COW privacy compile")
}

#[test]
fn generated_scoped_cow_capability_and_borrow_lifetimes_are_compile_time_boundaries() {
    for strict in [false, true] {
        let pass = cargo_check("atom_code_cow_allowed", strict);
        assert!(
            pass.status.success(),
            "{}",
            String::from_utf8_lossy(&pass.stderr)
        );
        for (case, diagnostic) in [
            ("atom_code_cow_missing_writes", "E0599"),
            ("atom_code_cow_runtime_private", "private"),
            ("atom_code_cow_borrow_escape", "lifetime"),
        ] {
            let fail = cargo_check(case, strict);
            assert!(!fail.status.success(), "case={case} strict={strict}");
            let errors = String::from_utf8_lossy(&fail.stderr);
            assert!(
                errors.contains(diagnostic),
                "case={case} strict={strict}: {errors}"
            );
            if case != "atom_code_cow_borrow_escape" {
                assert!(errors.contains("stage_topology_properties"), "{errors}");
            }
        }
    }
}
