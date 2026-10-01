use std::collections::BTreeSet;
use std::path::PathBuf;
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};

const CORE: &[&str] = &[
    "cap-smiles",
    "cap-io",
    "cap-serialization",
    "cap-descriptors",
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
    "cap-stereoisomers",
    "cap-tautomer",
];
const BUNDLES: &[(&str, &[&str])] = &[
    ("core", CORE),
    ("bio", &["cap-bio"]),
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

// Independent expected membership, not a second production feature registry.
fn all_caps() -> BTreeSet<&'static str> {
    BUNDLES
        .iter()
        .flat_map(|(_, caps)| caps.iter().copied())
        .collect()
}

struct Probe(PathBuf);

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
        Self(dir)
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
        Command::new(std::env::var("CARGO").unwrap_or_else(|_| "cargo".into()))
            .current_dir(&self.0)
            .env(
                "CARGO_TARGET_DIR",
                root.join("target/feature-selection-compile"),
            )
            .env("CARGO_BUILD_JOBS", "4")
            .args(args)
            .output()
            .unwrap()
    }

    fn resolved_caps(&self) -> BTreeSet<String> {
        let output = self.cargo(&["metadata", "--offline", "--format-version", "1"]);
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let metadata: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
        let package = metadata["packages"]
            .as_array()
            .unwrap()
            .iter()
            .find(|p| p["name"] == "cosmolkit")
            .unwrap();
        let node = metadata["resolve"]["nodes"]
            .as_array()
            .unwrap()
            .iter()
            .find(|n| n["id"] == package["id"])
            .unwrap();
        node["features"]
            .as_array()
            .unwrap()
            .iter()
            .map(|f| f.as_str().unwrap())
            .filter(|f| f.starts_with("cap-"))
            .map(str::to_owned)
            .collect()
    }

    fn check_source(&self, source: &str) -> Output {
        std::fs::write(self.0.join("src/lib.rs"), source).unwrap();
        self.cargo(&["check", "--offline", "--release", "--quiet", "--lib"])
    }
}

impl Drop for Probe {
    fn drop(&mut self) {
        std::fs::remove_dir_all(&self.0).unwrap();
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
    assert_eq!(cases.len(), 51); // 40 base + 10 subsets + the README core/bio example.
    for strict in [false, true] {
        for (defaults, features, expected) in &cases {
            let mut selected = features.clone();
            if strict {
                selected.push("op-contracts-strict");
            }
            probe.configure(*defaults, &selected);
            assert_eq!(
                probe.resolved_caps(),
                expected.iter().map(|s| (*s).to_owned()).collect(),
                "default-features={defaults}, features={selected:?}"
            );
            let compiled = probe.cargo(&["check", "--offline", "--release", "--quiet", "--lib"]);
            assert!(
                compiled.status.success(),
                "default-features={defaults}, features={selected:?}: {}",
                String::from_utf8_lossy(&compiled.stderr)
            );
        }
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
            &["core"],
            &[
                "from_smiles",
                "from_sdf",
                "with_hydrogens",
                "with_kekulized_bonds",
                "sanitize",
                "molecular_weight",
                "potential_stereo",
            ],
            &["with_2d_coordinates"],
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
    assert_eq!(cases.len(), 20); // Full 16-way four-cap product + four other API boundaries.
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
                "{selected:?}: {}",
                String::from_utf8_lossy(&output.stderr)
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
                "unselected methods compiled: {selected:?}"
            );
            let stderr = String::from_utf8_lossy(&output.stderr);
            assert_eq!(
                stderr.matches("error[E0599]").count(),
                forbidden.len(),
                "{stderr}"
            );
            for method in forbidden {
                assert!(stderr.contains(&format!("`{method}`")), "{stderr}");
            }
        }
    }
}
