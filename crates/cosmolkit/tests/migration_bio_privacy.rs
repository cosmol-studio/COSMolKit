use std::process::Command;

#[test]
fn bio_real_operation_modules_expose_only_declared_fields_in_default_and_strict() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap();
    for strict in [false, true] {
        for case in [
            "allowed",
            "undeclared",
            "storage",
            "replacement_readonly",
            "metadata_readonly",
            "replacement_undeclared",
            "replacement_required",
        ] {
            // Probe cfgs belong to this facade only. Dependencies do not need
            // to be rebuilt for each deliberately invalid operation body.
            let case_cfg = format!("cosmolkit_bio_privacy_case=\"{case}\"");
            let output = Command::new(std::env::var("CARGO").unwrap_or_else(|_| "cargo".into()))
                .current_dir(root)
                .env(
                    "CARGO_TARGET_DIR",
                    root.join("target/runtime-privacy-compile"),
                )
                .env("CARGO_BUILD_JOBS", "4")
                .env("CARGO_INCREMENTAL", "0")
                .args([
                    "rustc",
                    "--quiet",
                    "--locked",
                    "-p",
                    "cosmolkit",
                    "--lib",
                    "--no-default-features",
                    "--features",
                    if strict {
                        "cap-bio,op-contracts-strict"
                    } else {
                        "cap-bio"
                    },
                    "--",
                    "--emit=metadata",
                    "--cfg",
                    "cosmolkit_bio_privacy_probe",
                    "--cfg",
                    &case_cfg,
                ])
                .output()
                .unwrap();
            let errors = String::from_utf8_lossy(&output.stderr);
            if case == "allowed" {
                assert!(output.status.success(), "{errors}");
            } else {
                assert!(!output.status.success(), "{case} unexpectedly compiled");
                if case == "undeclared" {
                    assert!(errors.contains("no field `atoms`"), "{errors}");
                } else if case == "storage" {
                    assert!(
                        errors.contains("private")
                            && errors.contains("BioStructure")
                            && errors.contains("Protein"),
                        "{errors}"
                    );
                } else {
                    let expected = match case {
                        "replacement_readonly" => "E0594",
                        "metadata_readonly" => "E0308",
                        "replacement_undeclared" => "no field `assemblies`",
                        "replacement_required" => "missing field `source_state`",
                        _ => unreachable!(),
                    };
                    assert!(errors.contains(expected), "{case}: {errors}");
                    assert!(errors.contains("bio.rs"), "{errors}");
                }
            }
        }
    }
}
