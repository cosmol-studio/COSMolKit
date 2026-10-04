use std::collections::BTreeSet;
use std::path::PathBuf;
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};

const CORE: &[&str] = &[
    "cap-smiles",
    "cap-hydrogens",
    "cap-valence",
    "cap-radicals",
    "cap-rings",
    "cap-matrices",
    "cap-transforms",
    "cap-stereo",
    "cap-kekulize",
    "cap-aromaticity",
    "cap-sanitize",
];
const BUNDLES: &[(&str, &[&str])] = &[
    ("core", CORE),
    ("bio", &["cap-bio"]),
    ("descriptors", &["cap-descriptors"]),
    ("tautomer", &["cap-tautomer"]),
    ("io", &["cap-io", "cap-serialization"]),
    (
        "conformer",
        &["cap-conformer", "cap-confseq", "cap-alignment"],
    ),
    ("forcefields", &["cap-forcefields"]),
    ("fingerprints", &["cap-fingerprints", "cap-hashing"]),
    ("search", &["cap-search"]),
    ("depict", &["cap-depict"]),
    ("inchi", &["cap-inchi"]),
    ("batch", &["cap-batch"]),
];
const FOUR_CAPS: &[&str] = &["cap-io", "cap-kekulize", "cap-sanitize", "cap-hydrogens"];

#[test]
fn drawing_queries_and_error_follow_cap_depict_default_and_strict() {
    let probe = Probe::new();
    for strict in [false, true] {
        let mut features = vec!["cap-depict"];
        if strict {
            features.push("op-contracts-strict");
        }
        probe.configure(false, &features);
        let output = probe.check_source("pub fn probe() {\nlet _: fn(&cosmolkit::Molecule,u32,u32)->Result<String,cosmolkit::DrawingError> = cosmolkit::Molecule::to_svg;\nlet _: fn(&cosmolkit::Molecule,u32,u32)->Result<Vec<u8>,cosmolkit::DrawingError> = cosmolkit::Molecule::to_png;\n}\n");
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let output = probe.check_source("pub fn probe() {\nlet _ = cosmolkit::Molecule::from_smiles;\nlet _ = cosmolkit::Molecule::with_hydrogens;\nlet _ = cosmolkit::Molecule::with_kekulized_bonds;\nlet _ = cosmolkit::Molecule::with_assigned_rings;\n}\n");
        assert!(!output.status.success());
        assert_eq!(
            String::from_utf8_lossy(&output.stderr)
                .matches("error[E0599]")
                .count(),
            4
        );
        let absent = if strict {
            vec!["op-contracts-strict"]
        } else {
            Vec::new()
        };
        probe.configure(false, &absent);
        let output = probe.check_source("pub fn probe() {\nlet _ = cosmolkit::Molecule::to_svg;\nlet _ = cosmolkit::Molecule::to_png;\nlet _: Option<cosmolkit::DrawingError> = None;\n}\n");
        assert!(!output.status.success());
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(stderr.matches("error[E0599]").count(), 2, "{stderr}");
        assert_eq!(stderr.matches("error[E0425]").count(), 1, "{stderr}");
        assert!(stderr.contains("DrawingError"), "{stderr}");
    }
}

// Independent expected membership, not a second production feature registry.
fn all_caps() -> BTreeSet<&'static str> {
    BUNDLES
        .iter()
        .flat_map(|(_, caps)| caps.iter().copied())
        // Enumeration is an optional advanced capability, not basic stereo.
        // `full` retains it even though no coarse bundle currently groups it.
        .chain(["cap-stereoisomers"])
        .collect()
}

struct Probe(PathBuf, AtomicUsize);

impl Probe {
    fn new() -> Self {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let dir = std::env::temp_dir().join(format!(
            "ck-feature-probe-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&dir).unwrap();
        std::fs::create_dir(dir.join("src")).unwrap();
        Self(dir, AtomicUsize::new(0))
    }

    fn configure(&self, defaults: bool, features: &[&str]) {
        let dependency = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        std::fs::write(self.0.join("Cargo.toml"), format!(
            "[package]\nname = \"ck-feature-probe\"\nversion = \"0.0.0\"\nedition = \"2024\"\n\
             [dependencies]\ncosmolkit = {{ path = {:?}, default-features = {defaults}, features = {features:?} }}\n\
             [workspace]\n", dependency
        )).unwrap();
        std::fs::write(
            self.0.join("src/lib.rs"),
            "pub fn baseline(m: &cosmolkit::Molecule) -> usize { m.num_atoms() }\n",
        )
        .unwrap();
    }

    fn cargo(&self, args: &[&str]) -> Output {
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .parent()
            .unwrap()
            .to_owned();
        let program = std::env::var("CARGO").unwrap_or_else(|_| "cargo".into());
        let is_check = args.first() == Some(&"check");
        let mut command = Command::new(&program);
        command
            .current_dir(&self.0)
            .env(
                "CARGO_TARGET_DIR",
                // Each independently mutable probe needs its own artifacts.
                // Parallel tests use the same package/crate name and rewrite
                // src/lib.rs between positive and negative compilation checks.
                // Sharing their target directory cannot establish source
                // freshness for an individual check.
                root.join("target/feature-selection-compile")
                    .join(self.0.file_name().expect("named probe directory")),
            )
            .env("CARGO_BUILD_JOBS", "4")
            .args(args);
        if is_check {
            command.env("CARGO_LOG", "cargo::core::compiler::fingerprint=info");
        }
        let output = command.output().unwrap();
        if is_check {
            let sequence = self.1.fetch_add(1, Ordering::Relaxed);
            let evidence_dir = self.0.join(format!("cargo-check-{sequence:04}"));
            std::fs::create_dir(&evidence_dir).unwrap();
            let evidence = serde_json::json!({
                "program": program,
                "args": args,
                "cwd": self.0.display().to_string(),
                "target_dir": root.join("target/feature-selection-compile")
                    .join(self.0.file_name().expect("named probe directory"))
                    .display().to_string(),
                "profile": "release",
                "jobs": 4,
                "status_debug": format!("{:?}", output.status),
                "status_success": output.status.success(),
                "status_code": output.status.code(),
                "CARGO_LOG": "cargo::core::compiler::fingerprint=info",
            });
            std::fs::write(
                evidence_dir.join("command.json"),
                serde_json::to_vec_pretty(&evidence).unwrap(),
            )
            .unwrap();
            std::fs::write(evidence_dir.join("stdout.bin"), &output.stdout).unwrap();
            std::fs::write(evidence_dir.join("stderr.bin"), &output.stderr).unwrap();
        }
        output
    }

    fn metadata(&self) -> serde_json::Value {
        let output = self.cargo(&["metadata", "--offline", "--format-version", "1"]);
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        serde_json::from_slice(&output.stdout).unwrap()
    }

    fn resolved_caps(&self) -> (BTreeSet<String>, Output) {
        let output = self.cargo(&["metadata", "--offline", "--format-version", "1"]);
        assert!(
            output.status.success(),
            "Cargo metadata failed: {}\n{}",
            String::from_utf8_lossy(&output.stderr),
            self.failure_context(&output)
        );
        let metadata: serde_json::Value =
            serde_json::from_slice(&output.stdout).unwrap_or_else(|error| {
                panic!(
                    "Cargo metadata JSON parse failed: {error}\n{}",
                    self.failure_context(&output)
                )
            });
        let package = metadata["packages"]
            .as_array()
            .unwrap_or_else(|| {
                panic!(
                    "Cargo metadata packages is not an array\n{}",
                    self.failure_context(&output)
                )
            })
            .iter()
            .find(|p| p["name"] == "cosmolkit")
            .unwrap_or_else(|| {
                panic!(
                    "Cargo metadata has no cosmolkit package\n{}",
                    self.failure_context(&output)
                )
            });
        let node = metadata["resolve"]["nodes"]
            .as_array()
            .unwrap_or_else(|| {
                panic!(
                    "Cargo metadata resolve.nodes is not an array\n{}",
                    self.failure_context(&output)
                )
            })
            .iter()
            .find(|n| n["id"] == package["id"])
            .unwrap_or_else(|| {
                panic!(
                    "Cargo metadata resolve has no cosmolkit node\n{}",
                    self.failure_context(&output)
                )
            });
        let caps: BTreeSet<String> = node["features"]
            .as_array()
            .unwrap_or_else(|| {
                panic!(
                    "Cargo metadata cosmolkit features is not an array\n{}",
                    self.failure_context(&output)
                )
            })
            .iter()
            .map(|feature| {
                feature.as_str().unwrap_or_else(|| {
                    panic!(
                        "Cargo metadata contains a non-string cosmolkit feature\n{}",
                        self.failure_context(&output)
                    )
                })
            })
            .filter(|f| f.starts_with("cap-"))
            .map(str::to_owned)
            .collect();
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .parent()
            .unwrap()
            .to_owned();
        let evidence = serde_json::json!({
            "features": caps,
            "program": std::env::var("CARGO").unwrap_or_else(|_| "cargo".into()),
            "args": ["metadata", "--offline", "--format-version", "1"],
            "cwd": self.0.display().to_string(),
            "target_dir": root.join("target/feature-selection-compile")
                .join(self.0.file_name().expect("named probe directory"))
                .display().to_string(),
            "profile": "metadata",
            "jobs": 4,
            "status_debug": format!("{:?}", output.status),
            "status_success": output.status.success(),
            "status_code": output.status.code(),
        });
        std::fs::write(
            self.0.join("resolved-caps.json"),
            serde_json::to_vec_pretty(&evidence).unwrap(),
        )
        .unwrap();
        (caps, output)
    }

    fn failure_context(&self, output: &Output) -> String {
        let read_text = |path: &std::path::Path| match std::fs::read(path) {
            Ok(bytes) => String::from_utf8_lossy(&bytes).into_owned(),
            Err(error) => format!("<missing/unreadable {}: {error}>", path.display()),
        };
        let check_count = self.1.load(Ordering::Relaxed);
        let command_evidence = match check_count.checked_sub(1) {
            Some(sequence) => read_text(
                &self
                    .0
                    .join(format!("cargo-check-{sequence:04}"))
                    .join("command.json"),
            ),
            None => "<no Cargo check invocation captured for this Probe>".into(),
        };
        let output_text = |bytes: &[u8]| match std::str::from_utf8(bytes) {
            Ok(text) => text.to_owned(),
            Err(error) => format!("<invalid UTF-8 at byte {}: {bytes:?}>", error.valid_up_to()),
        };
        format!(
            "Probe diagnostics at {}\nmanifest:\n{}\ncurrent source:\n{}\nlast resolved feature/metadata evidence:\n{}\nlast Cargo check command evidence:\n{}\nOutput status: {:?} (success={}, code={:?})\nOutput stdout:\n{}\nOutput stderr:\n{}",
            self.0.display(),
            read_text(&self.0.join("Cargo.toml")),
            read_text(&self.0.join("src/lib.rs")),
            read_text(&self.0.join("resolved-caps.json")),
            command_evidence,
            output.status,
            output.status.success(),
            output.status.code(),
            output_text(&output.stdout),
            output_text(&output.stderr),
        )
    }

    fn resolved_dependencies(&self) -> BTreeSet<String> {
        let metadata = self.metadata();
        let nodes = metadata["resolve"]["nodes"].as_array().unwrap();
        let packages = metadata["packages"].as_array().unwrap();
        let root = metadata["resolve"]["root"].as_str().unwrap();
        let mut pending = vec![root];
        let mut visited = BTreeSet::new();
        let mut names = BTreeSet::new();
        while let Some(id) = pending.pop() {
            if !visited.insert(id) {
                continue;
            }
            let node = nodes.iter().find(|node| node["id"] == id).unwrap();
            let package = packages.iter().find(|package| package["id"] == id).unwrap();
            names.insert(package["name"].as_str().unwrap().to_owned());
            pending.extend(
                node["deps"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|dep| dep["pkg"].as_str().unwrap()),
            );
        }
        names
    }

    fn check_source(&self, source: &str) -> Output {
        std::fs::write(self.0.join("src/lib.rs"), source).unwrap();
        let _ = self.resolved_caps();
        self.cargo(&["check", "--offline", "--release", "--quiet", "--lib"])
    }
}

impl Drop for Probe {
    fn drop(&mut self) {
        if std::thread::panicking() {
            eprintln!(
                "retained failing feature Probe evidence at {}",
                self.0.display()
            );
        } else {
            std::fs::remove_dir_all(&self.0).unwrap();
        }
    }
}

#[test]
fn cargo_resolves_exact_bundle_and_individual_capability_sets() {
    let probe = Probe::new();
    let mut cases = vec![
        (true, vec![], all_caps()),
        (false, vec!["full"], all_caps()),
        (false, vec![], BTreeSet::new()),
        (
            false,
            FOUR_CAPS.to_vec(),
            FOUR_CAPS.iter().copied().collect(),
        ),
    ];
    for &(bundle, caps) in BUNDLES {
        cases.push((false, vec![bundle], caps.iter().copied().collect()));
    }
    for cap in all_caps() {
        cases.push((false, vec![cap], [cap].into_iter().collect()));
    }
    // Empty/all-four/single-capability cases above already cover six subsets.
    // Add all ten two/three-capability subsets, completing the 16-way product.
    for mask in 0_u8..16 {
        if !matches!(mask.count_ones(), 2 | 3) {
            continue;
        }
        let selected: Vec<_> = FOUR_CAPS
            .iter()
            .enumerate()
            .filter(|(i, _)| mask & (1 << i) != 0)
            .map(|(_, cap)| *cap)
            .collect();
        cases.push((false, selected.clone(), selected.into_iter().collect()));
    }
    cases.push((
        false,
        vec!["core", "bio"],
        CORE.iter().copied().chain(["cap-bio"]).collect(),
    ));
    assert_eq!(cases.len(), 54); // 43 base + 10 subsets + the README core/bio example.
    for strict in [false, true] {
        for (defaults, features, expected) in &cases {
            let mut selected = features.clone();
            if strict {
                selected.push("op-contracts-strict");
            }
            probe.configure(*defaults, &selected);
            let (resolved, metadata_output) = probe.resolved_caps();
            assert_eq!(
                resolved,
                expected.iter().map(|s| (*s).to_owned()).collect(),
                "default-features={defaults}, features={selected:?}\n{}",
                probe.failure_context(&metadata_output)
            );
            let compiled = probe.cargo(&["check", "--offline", "--release", "--quiet", "--lib"]);
            assert!(
                compiled.status.success(),
                "default-features={defaults}, features={selected:?}: {}\n{}",
                String::from_utf8_lossy(&compiled.stderr),
                probe.failure_context(&compiled)
            );
        }
    }
}

#[test]
fn lightweight_bundles_have_isolated_build_dependencies() {
    // Walk resolved edges in a separate consuming project: neither lockfile
    // entries nor a facade workspace's feature-unified graph prove isolation.
    let cases: &[(&[&str], &[&str], &[&str])] = &[
        (
            &[],
            &["cosmolkit-model", "cosmolkit-macros"],
            &[
                "cosmolkit-core",
                "cosmolkit-io",
                "cosmolkit-bio",
                "cosmolkit-smiles",
                "cosmolkit-stereo",
                "cosmolkit-descriptors",
                "cosmolkit-search",
                "cosmolkit-tautomer",
            ],
        ),
        (
            &["core"],
            &["cosmolkit-core", "cosmolkit-smiles", "cosmolkit-stereo"],
            &[
                "cosmolkit-io",
                "cosmolkit-bio",
                "cosmolkit-descriptors",
                "cosmolkit-search",
                "cosmolkit-tautomer",
            ],
        ),
        (
            &["bio"],
            &["cosmolkit-bio", "cosmolkit-io"],
            &[
                "cosmolkit-core",
                "cosmolkit-smiles",
                "cosmolkit-stereo",
                "cosmolkit-descriptors",
                "cosmolkit-search",
                "cosmolkit-tautomer",
            ],
        ),
        (
            &["bio", "core"],
            &[
                "cosmolkit-bio",
                "cosmolkit-io",
                "cosmolkit-core",
                "cosmolkit-smiles",
                "cosmolkit-stereo",
            ],
            &[
                "cosmolkit-descriptors",
                "cosmolkit-search",
                "cosmolkit-tautomer",
            ],
        ),
        (
            &["descriptors"],
            &[
                "cosmolkit-descriptors",
                "cosmolkit-core",
                "cosmolkit-search",
            ],
            &[
                "cosmolkit-bio",
                "cosmolkit-io",
                "cosmolkit-smiles",
                "cosmolkit-tautomer",
            ],
        ),
        (
            &["io"],
            &["cosmolkit-io", "cosmolkit-core", "cosmolkit-search"],
            &[
                "cosmolkit-bio",
                "cosmolkit-descriptors",
                "cosmolkit-smiles",
                "cosmolkit-tautomer",
            ],
        ),
        (
            &["tautomer"],
            &["cosmolkit-tautomer"],
            &[
                "cosmolkit-core",
                "cosmolkit-bio",
                "cosmolkit-io",
                "cosmolkit-smiles",
                "cosmolkit-stereo",
                "cosmolkit-descriptors",
                "cosmolkit-search",
            ],
        ),
        (
            &["cap-serialization"],
            &["cosmolkit-io", "cosmolkit-core", "cosmolkit-search"],
            &[
                "cosmolkit-bio",
                "cosmolkit-descriptors",
                "cosmolkit-tautomer",
            ],
        ),
        (
            &["full"],
            &[
                "cosmolkit-bio",
                "cosmolkit-io",
                "cosmolkit-core",
                "cosmolkit-smiles",
                "cosmolkit-stereo",
                "cosmolkit-descriptors",
                "cosmolkit-search",
                "cosmolkit-tautomer",
            ],
            &[],
        ),
    ];
    let probe = Probe::new();
    let mut calls = 0;
    for strict in [false, true] {
        for &(features, present, absent) in cases {
            let mut selected = features.to_vec();
            if strict {
                selected.push("op-contracts-strict");
            }
            probe.configure(false, &selected);
            let packages = probe.resolved_dependencies();
            calls += 1;
            for name in present {
                assert!(
                    packages.contains(*name),
                    "{selected:?}: missing {name}: {packages:?}"
                );
            }
            for name in absent {
                assert!(
                    !packages.contains(*name),
                    "{selected:?}: unexpected {name}: {packages:?}"
                );
            }
            if !features.contains(&"full") {
                assert!(
                    !probe.resolved_caps().0.contains("cap-search"),
                    "{selected:?}"
                );
            }
        }
    }
    assert_eq!(calls, 18);

    // A real BIO-only consuming build executes constructors/CID and checks
    // method availability: a lightweight graph must still deliver BIO behavior.
    for strict in [false, true] {
        let mut features = vec!["bio"];
        if strict {
            features.push("op-contracts-strict");
        }
        probe.configure(false, &features);
        let source = r#"
#[test]
fn bio_only_runtime_control() {
    let structure = cosmolkit::BioStructure::from_pdb("HEADER    TEST\n").unwrap();
    let selection = cosmolkit::BioSelection::from_cid("//*").unwrap();
    assert!(structure.selected_atom_ids(&selection).unwrap().is_empty());
    let _ = cosmolkit::BioStructure::from_mmcif;
}
"#;
        std::fs::write(probe.0.join("src/lib.rs"), source).unwrap();
        let output = probe.cargo(&["test", "--offline", "--release", "--quiet", "--lib"]);
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(String::from_utf8_lossy(&output.stdout).contains("1 passed"));
    }
}

#[test]
fn enabled_registry_rows_use_capability_names_not_bundle_names() {
    let enabled = [
        ("cap-alignment", cfg!(feature = "cap-alignment")),
        ("cap-batch", cfg!(feature = "cap-batch")),
        ("cap-bio", cfg!(feature = "cap-bio")),
        ("cap-conformer", cfg!(feature = "cap-conformer")),
        ("cap-confseq", cfg!(feature = "cap-confseq")),
        ("cap-depict", cfg!(feature = "cap-depict")),
        ("cap-fingerprints", cfg!(feature = "cap-fingerprints")),
        ("cap-forcefields", cfg!(feature = "cap-forcefields")),
        ("cap-hashing", cfg!(feature = "cap-hashing")),
        ("cap-inchi", cfg!(feature = "cap-inchi")),
        ("cap-io", cfg!(feature = "cap-io")),
        ("cap-search", cfg!(feature = "cap-search")),
        ("cap-serialization", cfg!(feature = "cap-serialization")),
        ("cap-smiles", cfg!(feature = "cap-smiles")),
        ("cap-stereoisomers", cfg!(feature = "cap-stereoisomers")),
        ("cap-tautomer", cfg!(feature = "cap-tautomer")),
        ("cap-hydrogens", cfg!(feature = "cap-hydrogens")),
        ("cap-descriptors", cfg!(feature = "cap-descriptors")),
        ("cap-valence", cfg!(feature = "cap-valence")),
        ("cap-radicals", cfg!(feature = "cap-radicals")),
        ("cap-rings", cfg!(feature = "cap-rings")),
        ("cap-matrices", cfg!(feature = "cap-matrices")),
        ("cap-transforms", cfg!(feature = "cap-transforms")),
        ("cap-stereo", cfg!(feature = "cap-stereo")),
        ("cap-kekulize", cfg!(feature = "cap-kekulize")),
        ("cap-aromaticity", cfg!(feature = "cap-aromaticity")),
        ("cap-sanitize", cfg!(feature = "cap-sanitize")),
    ];
    for row in cosmolkit::BINDING_CONTRACT {
        if row.feature.starts_with("cap-") {
            assert!(
                enabled.contains(&(row.feature, true)),
                "{}: {}",
                row.semantic_id,
                row.feature
            );
        } else {
            assert!(
                matches!(row.feature, "runtime" | "metadata"),
                "{}",
                row.feature
            );
        }
    }
    for row in cosmolkit::SUPPORT_MATRIX {
        assert!(enabled.contains(&(row.feature.name, true)));
    }
}

#[test]
fn public_methods_compile_only_with_their_own_capability_default_and_strict() {
    let cases: &[(&[&str], &[&str], &[&str])] = &[
        (
            &[],
            &[],
            &[
                "from_smiles",
                "from_sdf",
                "with_hydrogens",
                "with_kekulized_bonds",
                "sanitize",
                "molecular_weight",
            ],
        ),
        (
            FOUR_CAPS,
            &[
                "from_sdf",
                "with_hydrogens",
                "with_kekulized_bonds",
                "sanitize",
            ],
            &[
                "from_smiles",
                "molecular_weight",
                "with_assigned_aromaticity",
                "potential_stereo",
            ],
        ),
        (
            &["cap-smiles"],
            &["from_smiles"],
            &[
                "from_sdf",
                "with_hydrogens",
                "with_kekulized_bonds",
                "sanitize",
            ],
        ),
        (&["cap-bio"], &[], &["from_sdf", "from_smiles"]),
        (
            &["cap-forcefields"],
            &[
                "with_uff_optimized_coordinates",
                "with_uff_optimized_coordinates_with_params",
                "with_uff_optimized_conformers",
                "with_uff_optimized_conformers_with_params",
            ],
            &["from_smiles"],
        ),
        (
            &["core"],
            &[
                "from_smiles",
                "with_hydrogens",
                "with_kekulized_bonds",
                "sanitize",
                "potential_stereo",
            ],
            &["from_sdf", "molecular_weight", "with_2d_coordinates"],
        ),
        (
            &["full"],
            &[
                "from_smiles",
                "from_sdf",
                "with_hydrogens",
                "with_kekulized_bonds",
                "sanitize",
                "molecular_weight",
                "with_2d_coordinates",
            ],
            &[],
        ),
    ];
    let mut cases: Vec<_> = cases
        .iter()
        .map(|(features, allowed, forbidden)| {
            (features.to_vec(), allowed.to_vec(), forbidden.to_vec())
        })
        .collect();
    let four_methods = [
        "from_sdf",
        "with_kekulized_bonds",
        "sanitize",
        "with_hydrogens",
    ];
    for mask in 1_u8..15 {
        let mut features = Vec::new();
        let mut allowed = Vec::new();
        let mut forbidden = vec![
            "from_smiles",
            "molecular_weight",
            "with_assigned_aromaticity",
            "potential_stereo",
        ];
        for (i, cap) in FOUR_CAPS.iter().enumerate() {
            if mask & (1 << i) != 0 {
                features.push(*cap);
                allowed.push(four_methods[i]);
            } else {
                forbidden.push(four_methods[i]);
            }
        }
        cases.push((features, allowed, forbidden));
    }
    assert_eq!(cases.len(), 21); // Full 16-way four-cap product + five other API boundaries.
    let probe = Probe::new();
    for strict in [false, true] {
        for (features, allowed, forbidden) in &cases {
            let mut selected = features.to_vec();
            if strict {
                selected.push("op-contracts-strict");
            }
            probe.configure(false, &selected);
            let mut source =
                String::from("pub fn probe() { let _ = cosmolkit::Molecule::num_atoms;\n");
            for method in allowed {
                source.push_str(&format!("let _ = cosmolkit::Molecule::{method};\n"));
            }
            if features.contains(&"cap-bio") {
                source.push_str("let _ = cosmolkit::BioStructure::from_pdb;\n");
            }
            source.push_str("}\n");
            let output = probe.check_source(&source);
            assert!(
                output.status.success(),
                "{selected:?}: {}\n{}",
                String::from_utf8_lossy(&output.stderr),
                probe.failure_context(&output)
            );
            if forbidden.is_empty() {
                continue;
            }
            source = "pub fn probe() {\n".into();
            for method in forbidden {
                source.push_str(&format!("let _ = cosmolkit::Molecule::{method};\n"));
            }
            source.push_str("}\n");
            let output = probe.check_source(&source);
            assert!(
                !output.status.success(),
                "unselected methods compiled: {selected:?}\n{}",
                probe.failure_context(&output)
            );
            let stderr = String::from_utf8_lossy(&output.stderr);
            assert_eq!(
                stderr.matches("error[E0599]").count(),
                forbidden.len(),
                "{stderr}\n{}",
                probe.failure_context(&output)
            );
            for method in forbidden {
                assert!(
                    stderr.contains(&format!("`{method}`")),
                    "{stderr}\n{}",
                    probe.failure_context(&output)
                );
            }
        }
    }
}

#[test]
fn descriptor_query_feature_gates_are_exact_for_five_queries_and_errors() {
    // Real external Rust method-item lookups with positive controls: the
    // five descriptor count queries and both public error types compile
    // under cap-descriptors ALONE (defaults disabled), fail as missing
    // API with NO capabilities (E0599 for methods, E0412/E0432 as
    // applicable for the error types), and cap-descriptors does not
    // accidentally expose SMILES, hydrogen, stereo, ring or search
    // public methods.
    const QUERIES: &[&str] = &[
        "num_heavy_atoms",
        "total_atom_count",
        "lipinski_hba",
        "lipinski_hbd",
        "fraction_csp3",
    ];
    const UNRELATED: &[&str] = &[
        "from_smiles",
        "with_hydrogens",
        "with_kekulized_bonds",
        "sanitize",
        "potential_stereo",
        "symmorph",
    ];
    let probe = Probe::new();

    // Compile-pass under cap-descriptors alone, defaults disabled: the
    // five method items plus both error types as real type paths.
    probe.configure(false, &["cap-descriptors"]);
    let mut source = String::from("pub fn probe() {\n");
    for method in QUERIES {
        source.push_str(&format!("let _ = cosmolkit::Molecule::{method};\n"));
    }
    source.push_str("let _: Option<cosmolkit::DescriptorError> = None;\n");
    source.push_str("let _: Option<cosmolkit::DescriptorReadError> = None;\n");
    source.push_str("}\n");
    let output = probe.check_source(&source);
    assert!(
        output.status.success(),
        "cap-descriptors alone must compile the five queries + errors: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Unrelated capability methods stay absent under cap-descriptors.
    let mut source = String::from("pub fn probe() {\n");
    for method in UNRELATED {
        source.push_str(&format!("let _ = cosmolkit::Molecule::{method};\n"));
    }
    source.push_str("}\n");
    let output = probe.check_source(&source);
    assert!(
        !output.status.success(),
        "cap-descriptors must not expose unrelated capability methods"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(
        stderr.matches("error[E0599]").count(),
        UNRELATED.len(),
        "{stderr}"
    );

    // No capabilities at all: the five methods fail E0599 and both error
    // types fail as unresolved paths (E0412/E0432 as applicable).
    probe.configure(false, &[]);
    let mut source = String::from("pub fn probe() {\n");
    for method in QUERIES {
        source.push_str(&format!("let _ = cosmolkit::Molecule::{method};\n"));
    }
    source.push_str("let _: Option<cosmolkit::DescriptorError> = None;\n");
    source.push_str("let _: Option<cosmolkit::DescriptorReadError> = None;\n");
    source.push_str("}\n");
    let output = probe.check_source(&source);
    assert!(
        !output.status.success(),
        "no capabilities must not compile the queries or error types"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(
        stderr.matches("error[E0599]").count(),
        QUERIES.len(),
        "{stderr}"
    );
    for method in QUERIES {
        assert!(stderr.contains(&format!("`{method}`")), "{stderr}");
    }
    assert!(
        stderr.contains("DescriptorError") && stderr.contains("DescriptorReadError"),
        "both error types must fail as unresolved: {stderr}"
    );
}

/// RING-LIVE-PUBLIC F1: the eleven ring queries follow the exact same
/// capability discipline with real external compiler evidence, and the
/// neighboring capabilities stay independently gated.
#[test]
fn ring_query_feature_gates_are_exact_for_eleven_queries_and_errors() {
    const RING_QUERIES: &[&str] = &[
        "num_rings",
        "num_heterocycles",
        "num_aromatic_rings",
        "num_saturated_rings",
        "num_aliphatic_rings",
        "num_aromatic_heterocycles",
        "num_aromatic_carbocycles",
        "num_aliphatic_heterocycles",
        "num_aliphatic_carbocycles",
        "num_saturated_heterocycles",
        "num_saturated_carbocycles",
    ];
    let probe = Probe::new();

    // cap-descriptors ALONE (defaults disabled): all eleven compile with
    // the typed error; constructor/assignment/mutator APIs stay absent.
    probe.configure(false, &["cap-descriptors"]);
    let mut source = String::from("pub fn probe() {\n");
    for method in RING_QUERIES {
        source.push_str(&format!("let _ = cosmolkit::Molecule::{method};\n"));
    }
    source.push_str("let _: Option<cosmolkit::DescriptorReadError> = None;\n");
    source.push_str("}\n");
    let output = probe.check_source(&source);
    assert!(
        output.status.success(),
        "cap-descriptors alone must compile the eleven ring queries: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    const NOT_DESCRIPTORS: &[&str] = &["from_smiles", "with_assigned_rings", "add_hydrogens_"];
    let mut source = String::from("pub fn probe() {\n");
    for method in NOT_DESCRIPTORS {
        source.push_str(&format!("let _ = cosmolkit::Molecule::{method};\n"));
    }
    source.push_str("}\n");
    let output = probe.check_source(&source);
    assert!(
        !output.status.success(),
        "cap-descriptors must not expose constructor/assignment/mutators"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(
        stderr.matches("error[E0599]").count(),
        NOT_DESCRIPTORS.len(),
        "{stderr}"
    );

    // Empty selection: the eleven methods fail E0599.
    probe.configure(false, &[]);
    let mut source = String::from("pub fn probe() {\n");
    for method in RING_QUERIES {
        source.push_str(&format!("let _ = cosmolkit::Molecule::{method};\n"));
    }
    source.push_str("}\n");
    let output = probe.check_source(&source);
    assert!(
        !output.status.success(),
        "empty selection must not compile the ring queries"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(
        stderr.matches("error[E0599]").count(),
        RING_QUERIES.len(),
        "{stderr}"
    );

    // cap-smiles ALONE: the constructor compiles while the ring queries
    // and the hydrogen API stay absent.
    probe.configure(false, &["cap-smiles"]);
    let mut source = String::from("pub fn probe() {\n");
    source.push_str("let _ = cosmolkit::Molecule::from_smiles;\n");
    source.push_str("}\n");
    let output = probe.check_source(&source);
    assert!(
        output.status.success(),
        "cap-smiles alone must compile the constructor: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let mut source = String::from("pub fn probe() {\n");
    for method in RING_QUERIES {
        source.push_str(&format!("let _ = cosmolkit::Molecule::{method};\n"));
    }
    source.push_str("let _ = cosmolkit::Molecule::with_hydrogens;\n");
    source.push_str("}\n");
    let output = probe.check_source(&source);
    assert!(
        !output.status.success(),
        "cap-smiles alone must not expose the queries or hydrogen API: stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(
        stderr.matches("error[E0599]").count(),
        RING_QUERIES.len() + 1,
        "{stderr}"
    );

    // cap-rings ALONE: the assignment API still compiles.
    probe.configure(false, &["cap-rings"]);
    let source = "pub fn probe() {\nlet _ = cosmolkit::Molecule::with_assigned_rings;\n}\n";
    let output = probe.check_source(source);
    assert!(
        output.status.success(),
        "cap-rings alone must compile the assignment API: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // cap-smiles + cap-descriptors WITHOUT cap-rings: constructor and all
    // eleven queries compile together.
    probe.configure(false, &["cap-smiles", "cap-descriptors"]);
    let mut source = String::from("pub fn probe() {\n");
    source.push_str("let _ = cosmolkit::Molecule::from_smiles;\n");
    for method in RING_QUERIES {
        source.push_str(&format!("let _ = cosmolkit::Molecule::{method};\n"));
    }
    source.push_str("}\n");
    let output = probe.check_source(&source);
    assert!(
        output.status.success(),
        "cap-smiles+cap-descriptors must compile constructor + queries: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}
#[test]
fn feature_probe_evidence_normal_drop_removes_its_exact_root() {
    let probe = Probe::new();
    let root = probe.0.clone();
    probe.configure(false, &[]);
    assert!(root.is_dir());
    drop(probe);
    assert!(
        !root.exists(),
        "normal Drop left its Probe root at {}",
        root.display()
    );
}

#[test]
fn feature_probe_evidence_panic_drop_retains_exact_payloads() {
    let probe = Probe::new();
    let root = probe.0.clone();
    probe.configure(false, &[]);

    let dependency = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let expected_manifest = format!(
        "[package]\nname = \"ck-feature-probe\"\nversion = \"0.0.0\"\nedition = \"2024\"\n[dependencies]\ncosmolkit = {{ path = {:?}, default-features = false, features = [] }}\n[workspace]\n",
        dependency
    );
    let expected_source = b"pub fn baseline(m: &cosmolkit::Molecule) -> usize { m.num_atoms() }\n";
    assert_eq!(
        std::fs::read(root.join("Cargo.toml")).unwrap(),
        expected_manifest.as_bytes()
    );
    assert_eq!(
        std::fs::read(root.join("src/lib.rs")).unwrap(),
        expected_source
    );

    let (resolved, metadata_output) = probe.resolved_caps();
    let check_count_before_format = probe.1.load(Ordering::Relaxed);
    assert_eq!(check_count_before_format, 0);
    let metadata_context = probe.failure_context(&metadata_output);
    assert_eq!(probe.1.load(Ordering::Relaxed), check_count_before_format);
    assert!(metadata_context.contains(&root.display().to_string()));
    assert!(metadata_context.contains(&expected_manifest));
    assert!(metadata_context.contains(std::str::from_utf8(expected_source).unwrap()));
    assert!(metadata_context.contains("<no Cargo check invocation captured for this Probe>"));
    assert!(metadata_context.contains(String::from_utf8_lossy(&metadata_output.stdout).as_ref()));
    let caps_path = root.join("resolved-caps.json");
    let caps_bytes = std::fs::read(&caps_path).unwrap();
    let caps_json: serde_json::Value = serde_json::from_slice(&caps_bytes).unwrap();
    assert_eq!(
        caps_json["features"],
        serde_json::to_value(&resolved).unwrap()
    );
    assert_eq!(caps_bytes, serde_json::to_vec_pretty(&caps_json).unwrap());

    let check_output = probe.cargo(&["check", "--offline", "--release", "--quiet", "--lib"]);
    assert!(
        check_output.status.success(),
        "{}",
        probe.failure_context(&check_output)
    );
    let workspace_root = dependency.parent().unwrap().parent().unwrap();
    let evidence_dir = root.join("cargo-check-0000");
    let payload_paths = [
        evidence_dir.join("command.json"),
        evidence_dir.join("stdout.bin"),
        evidence_dir.join("stderr.bin"),
        caps_path,
    ];
    let mut unique_paths = BTreeSet::new();
    for path in &payload_paths {
        unique_paths.insert(path.clone());
    }
    assert_eq!(unique_paths.len(), 4);
    let payload_bytes: Vec<Vec<u8>> = payload_paths
        .iter()
        .map(|path| std::fs::read(path).unwrap())
        .collect();
    assert_eq!(payload_bytes[1], check_output.stdout);
    assert_eq!(payload_bytes[2], check_output.stderr);

    let command_json: serde_json::Value = serde_json::from_slice(&payload_bytes[0]).unwrap();
    let expected_command = serde_json::json!({
        "program": std::env::var("CARGO").unwrap_or_else(|_| "cargo".into()),
        "args": ["check", "--offline", "--release", "--quiet", "--lib"],
        "cwd": root.display().to_string(),
        "target_dir": workspace_root.join("target/feature-selection-compile").join(root.file_name().unwrap()).display().to_string(),
        "profile": "release",
        "jobs": 4,
        "status_debug": format!("{:?}", check_output.status),
        "status_success": check_output.status.success(),
        "status_code": check_output.status.code(),
        "CARGO_LOG": "cargo::core::compiler::fingerprint=info",
    });
    assert_eq!(command_json, expected_command);
    assert_eq!(
        payload_bytes[0],
        serde_json::to_vec_pretty(&expected_command).unwrap()
    );

    let expected_payloads = payload_bytes.clone();
    let caught = std::panic::catch_unwind(std::panic::AssertUnwindSafe(move || {
        let _probe = probe;
        panic!("intentional feature Probe evidence-retention fixture");
    }));
    assert!(caught.is_err());
    assert!(root.is_dir(), "panic Drop deleted its Probe root");
    assert_eq!(
        std::fs::read(root.join("Cargo.toml")).unwrap(),
        expected_manifest.as_bytes()
    );
    assert_eq!(
        std::fs::read(root.join("src/lib.rs")).unwrap(),
        expected_source
    );
    for (path, expected) in payload_paths.iter().zip(expected_payloads) {
        assert_eq!(std::fs::read(path).unwrap(), expected, "{}", path.display());
    }
    std::fs::remove_dir_all(&root).unwrap();
}

#[test]
fn morgan_public_feature_external_values_and_capability_isolation() {
    const MORGAN_VALUES_SOURCE: &str = r#"
use cosmolkit::{
    FingerprintError, MorganFingerprintParams, MorganInvariants, MorganParams, Molecule,
    SparseCountFingerprint, SparseCountFingerprint32,
};

pub fn probe(molecule: &Molecule) -> usize {
    let _: MorganParams = MorganParams::default();
    let _: MorganFingerprintParams = MorganFingerprintParams::default();
    let _: MorganInvariants = MorganInvariants::Connectivity;
    let _ = MorganInvariants::Features;
    let _ = MorganInvariants::FeaturePatterns(Vec::new());
    let _: Option<&FingerprintError> = None;
    let _: Option<&SparseCountFingerprint> = None;
    let _: Option<&SparseCountFingerprint32> = None;
    molecule.num_atoms()
}

pub fn inspect_ring_error(error: &cosmolkit::OperationError) {
    match error {
        cosmolkit::OperationError::Rings(cause) => {
            let _: &(dyn std::error::Error + 'static) = cause;
        }
        _ => {}
    }
    let _ = std::fmt::format(format_args!("{error}"));
    let _ = std::error::Error::source(error);
}
"#;

    const METHOD_CAPS: &[(&str, &str)] = &[
        ("with_assigned_valence", "cap-valence"),
        ("with_assigned_rings", "cap-rings"),
        ("with_hydrogens", "cap-hydrogens"),
    ];

    for strict in [false, true] {
        let mut fingerprint_features = vec!["cap-fingerprints"];
        if strict {
            fingerprint_features.push("op-contracts-strict");
        }

        let probe = Probe::new();
        probe.configure(false, &fingerprint_features);
        let (resolved, metadata_output) = probe.resolved_caps();
        let metadata_context = probe.failure_context(&metadata_output);
        assert!(
            metadata_output.status.success(),
            "Morgan-only metadata command failed (strict={strict})\n{metadata_context}"
        );
        assert_eq!(
            resolved,
            BTreeSet::from(["cap-fingerprints".to_owned()]),
            "Morgan-only resolved capabilities (strict={strict})\n{metadata_context}"
        );

        let output = probe.check_source(MORGAN_VALUES_SOURCE);
        assert!(
            output.status.success(),
            "installed Morgan facade values did not compile (strict={strict}): {}\n{}",
            String::from_utf8_lossy(&output.stderr),
            probe.failure_context(&output)
        );

        for &(method, capability) in METHOD_CAPS {
            let method_source =
                format!("pub fn probe() {{ let _ = cosmolkit::Molecule::{method}; }}\n");

            probe.configure(false, &fingerprint_features);
            let output = probe.check_source(&method_source);
            let stderr = String::from_utf8_lossy(&output.stderr);
            assert!(
                !output.status.success(),
                "{method} unexpectedly compiled without {capability} (strict={strict})\n{}",
                probe.failure_context(&output)
            );
            assert_eq!(
                stderr.matches("error[E0599]").count(),
                1,
                "expected one missing method-item error for {method} (strict={strict})\n{stderr}\n{}",
                probe.failure_context(&output)
            );
            assert!(
                stderr.contains(&format!("`{method}`")),
                "missing method diagnostic did not name {method} (strict={strict})\n{stderr}\n{}",
                probe.failure_context(&output)
            );

            let mut enabled_features = vec!["cap-fingerprints", capability];
            if strict {
                enabled_features.push("op-contracts-strict");
            }
            probe.configure(false, &enabled_features);
            let output = probe.check_source(&method_source);
            assert!(
                output.status.success(),
                "{method} failed with only its owning capability enabled (strict={strict}): {}\n{}",
                String::from_utf8_lossy(&output.stderr),
                probe.failure_context(&output)
            );
        }
    }

    let mut search_metadata_observations = 0;
    for strict in [false, true] {
        for search_enabled in [false, true] {
            let mut features = vec!["cap-fingerprints"];
            if search_enabled {
                features.push("cap-search");
            }
            if strict {
                features.push("op-contracts-strict");
            }

            let probe = Probe::new();
            probe.configure(false, &features);
            let (resolved, metadata_output) = probe.resolved_caps();
            let metadata_context = probe.failure_context(&metadata_output);
            assert!(
                metadata_output.status.success(),
                "resolved-cap metadata command failed for {features:?}\n{metadata_context}"
            );
            let expected = if search_enabled {
                BTreeSet::from(["cap-fingerprints".to_owned(), "cap-search".to_owned()])
            } else {
                BTreeSet::from(["cap-fingerprints".to_owned()])
            };
            assert_eq!(
                resolved, expected,
                "unexpected resolved capabilities for {features:?}\n{metadata_context}"
            );
            search_metadata_observations += 1;
        }
    }
    assert_eq!(search_metadata_observations, 4);
}
