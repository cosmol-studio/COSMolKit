//! Feature regressions use the current Cargo build, never nested builds.
//! CI selects build profiles. Missing public APIs use native rustdoc tests.
use std::collections::BTreeSet;

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

fn manifest() -> toml::Value {
    toml::from_str(include_str!("../Cargo.toml")).unwrap()
}

// Check local feature declarations, not a simulation of Cargo's cross-package
// resolver. Actual builds remain in the CI feature matrix.
fn declared_caps(manifest: &toml::Value, defaults: bool, selected: &[&str]) -> BTreeSet<String> {
    let table = manifest["features"].as_table().unwrap();
    let mut pending: Vec<&str> = selected.to_vec();
    if defaults {
        pending.push("default");
    }
    let mut visited = BTreeSet::new();
    let mut caps = BTreeSet::new();
    while let Some(feature) = pending.pop() {
        if !visited.insert(feature) {
            continue;
        }
        if feature.starts_with("cap-") {
            caps.insert(feature.to_owned());
        }
        for value in table
            .get(feature)
            .unwrap_or_else(|| panic!("unknown feature {feature}"))
            .as_array()
            .unwrap()
        {
            let value = value.as_str().unwrap();
            if !value.starts_with("dep:") && !value.contains('/') {
                pending.push(value);
            }
        }
    }
    caps
}

#[test]
fn bundle_and_capability_declarations_have_exact_membership() {
    let manifest = manifest();
    let all: BTreeSet<_> = BUNDLES
        .iter()
        .flat_map(|(_, caps)| caps.iter().copied())
        .chain(["cap-stereoisomers"])
        .collect();
    let declared: BTreeSet<_> = manifest["features"]
        .as_table()
        .unwrap()
        .keys()
        .filter(|name| name.starts_with("cap-"))
        .map(String::as_str)
        .collect();
    assert_eq!(declared, all);
    let mut cases = vec![
        (true, vec![], all.clone()),
        (false, vec!["full"], all.clone()),
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
    for &cap in &all {
        cases.push((false, vec![cap], [cap].into_iter().collect()));
    }
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
    assert_eq!(cases.len(), 54);
    for strict in [false, true] {
        for (defaults, selected, expected) in &cases {
            let mut selected = selected.clone();
            if strict {
                selected.push("op-contracts-strict");
            }
            assert_eq!(
                declared_caps(&manifest, *defaults, &selected),
                expected.iter().map(|s| (*s).to_owned()).collect(),
                "defaults={defaults}, selected={selected:?}",
            );
        }
    }
}

#[test]
fn optional_dependencies_and_io_branches_stay_independent() {
    let manifest = manifest();
    let dependencies = manifest["dependencies"].as_table().unwrap();
    for (name, dependency) in dependencies {
        assert_eq!(
            dependency
                .get("optional")
                .and_then(toml::Value::as_bool)
                .unwrap_or(false),
            !matches!(name.as_str(), "cosmolkit-model" | "cosmolkit-macros"),
            "{name}",
        );
    }
    assert_eq!(
        dependencies["cosmolkit-io"]["default-features"].as_bool(),
        Some(false)
    );
    for (feature, expected) in [
        (
            "cap-bio",
            &["dep:cosmolkit-bio", "dep:cosmolkit-io", "cosmolkit-io/bio"][..],
        ),
        (
            "cap-io",
            &[
                "dep:cosmolkit-io",
                "cosmolkit-io/molecule",
                "dep:cosmolkit-core",
            ][..],
        ),
        (
            "cap-serialization",
            &[
                "dep:cosmolkit-io",
                "cosmolkit-io/molecule",
                "dep:cosmolkit-core",
            ][..],
        ),
        (
            "cap-fingerprints",
            &[
                "dep:cosmolkit-fingerprints",
                "dep:cosmolkit-core",
                "dep:cosmolkit-io",
                "cosmolkit-io/molecule",
            ][..],
        ),
    ] {
        let actual: BTreeSet<_> = manifest["features"][feature]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        assert_eq!(actual, expected.iter().copied().collect(), "{feature}");
    }
    let io: toml::Value = toml::from_str(include_str!("../../cosmolkit-io/Cargo.toml")).unwrap();
    assert_eq!(
        io["features"]["default"].as_array().unwrap()[0].as_str(),
        Some("molecule")
    );
    assert_eq!(io["features"]["default"].as_array().unwrap().len(), 1);
    assert_eq!(
        io["features"]["bio"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect::<Vec<_>>(),
        ["dep:cosmolkit-bio"],
    );
    let molecule: BTreeSet<_> = io["features"]["molecule"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap())
        .collect();
    for dependency in ["dep:cosmolkit-core", "dep:cosmolkit-search"] {
        assert!(molecule.contains(dependency));
    }
    assert!(!molecule.contains("bio"));
    assert!(!molecule.contains("dep:cosmolkit-bio"));
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
                "{}",
                row.semantic_id
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

// Method references compile once with this target, without executing chemistry.
#[test]
fn selected_public_methods_are_available() {
    macro_rules! methods {
        ($feature:literal; $($method:ident),+ $(,)?) => {
            #[cfg(feature = $feature)]
            { $(let _ = cosmolkit::Molecule::$method;)+ }
        };
    }
    let _ = cosmolkit::Molecule::num_atoms;
    let _ = cosmolkit::Molecule::new;
    methods!("cap-smiles"; from_smiles);
    methods!("cap-io"; from_sdf);
    methods!("cap-hydrogens"; with_hydrogens, add_hydrogens_);
    methods!("cap-kekulize"; with_kekulized_bonds);
    methods!("cap-sanitize"; sanitize);
    methods!("cap-aromaticity"; with_assigned_aromaticity);
    methods!("cap-stereo"; potential_stereo);
    methods!("cap-rings"; with_assigned_rings);
    methods!("cap-valence"; with_assigned_valence);
    methods!("cap-depict"; with_2d_coordinates, to_svg, to_png);
    methods!("cap-forcefields";
        with_uff_optimized, with_uff_optimized_with_params,
        with_uff_optimized_confs, with_uff_optimized_confs_with_params,
        uff_has_all_molecule_params,
    );
    methods!("cap-descriptors";
        molecular_weight, num_heavy_atoms, total_atom_count, lipinski_hba,
        lipinski_hbd, fraction_csp3, num_rings, num_heterocycles,
        num_aromatic_rings, num_saturated_rings, num_aliphatic_rings,
        num_aromatic_heterocycles, num_aromatic_carbocycles,
        num_aliphatic_heterocycles, num_aliphatic_carbocycles,
        num_saturated_heterocycles, num_saturated_carbocycles,
    );
    #[cfg(feature = "cap-bio")]
    {
        let _ = cosmolkit::BioStructure::from_pdb;
        let _ = cosmolkit::BioStructure::from_mmcif;
    }
    #[cfg(feature = "cap-depict")]
    {
        let _: fn(&cosmolkit::Molecule, u32, u32) -> Result<String, cosmolkit::DrawingError> =
            cosmolkit::Molecule::to_svg;
        let _: fn(&cosmolkit::Molecule, u32, u32) -> Result<Vec<u8>, cosmolkit::DrawingError> =
            cosmolkit::Molecule::to_png;
    }
    #[cfg(feature = "cap-descriptors")]
    {
        let _: Option<cosmolkit::DescriptorError> = None;
        let _: Option<cosmolkit::DescriptorReadError> = None;
    }
    #[cfg(feature = "cap-fingerprints")]
    {
        let _: cosmolkit::MorganParams = cosmolkit::MorganParams::default();
        let _: cosmolkit::MorganFingerprintParams = cosmolkit::MorganFingerprintParams::default();
        let _ = [
            cosmolkit::MorganInvariants::Connectivity,
            cosmolkit::MorganInvariants::Features,
            cosmolkit::MorganInvariants::FeaturePatterns(Vec::new()),
        ];
        let _: Option<cosmolkit::FingerprintError> = None;
        let _: Option<cosmolkit::SparseCountFingerprint> = None;
        let _: Option<cosmolkit::SparseCountFingerprint32> = None;
        fn inspect_ring_error(error: &cosmolkit::OperationError) {
            if let cosmolkit::OperationError::Rings(cause) = error {
                let _: &(dyn std::error::Error + 'static) = cause;
            }
            let _ = std::fmt::format(format_args!("{error}"));
            let _ = std::error::Error::source(error);
        }
        let _ = inspect_ring_error;
    }
}

#[cfg(feature = "cap-bio")]
#[test]
fn bio_constructors_and_cid_selection_work() {
    let structure = cosmolkit::BioStructure::from_pdb("HEADER    TEST\n").unwrap();
    let selection = cosmolkit::BioSelection::from_cid("//*").unwrap();
    assert!(structure.selected_atom_ids(&selection).unwrap().is_empty());
}
