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
    "cap-io",
    "cap-serialization",
    "cap-batch",
];
const BUNDLES: &[(&str, &[&str])] = &[
    ("core", CORE),
    ("bio", &["cap-bio"]),
    ("descriptors", &["cap-descriptors"]),
    ("tautomer", &["cap-tautomer"]),
    (
        "conformer",
        &[
            "cap-conformer",
            "cap-confseq",
            "cap-alignment",
            "cap-forcefields",
        ],
    ),
    ("fingerprints", &["cap-fingerprints", "cap-hashing"]),
    ("search", &["cap-search"]),
    ("reaction", &["cap-reaction"]),
    ("stereoisomers", &["cap-stereoisomers"]),
    ("depict", &["cap-depict"]),
    ("inchi", &["cap-inchi"]),
];
const FOUR_CAPS: &[&str] = &["cap-io", "cap-kekulize", "cap-sanitize", "cap-hydrogens"];

// Independent expectations for public functionality included by prerequisites.
const PREREQUISITES: &[(&str, &[&str])] = &[
    ("cap-batch", &["cap-io"]),
    (
        "cap-conformer",
        &["cap-forcefields", "cap-alignment", "cap-io"],
    ),
    ("cap-confseq", &["cap-conformer"]),
    // Avalon converts through the molecular IO owner, without enabling depiction.
    ("cap-fingerprints", &["cap-io"]),
    ("cap-search", &["cap-smiles"]),
    ("cap-serialization", &["cap-io"]),
    ("cap-stereoisomers", &["cap-stereo"]),
    ("cap-tautomer", &["cap-search"]),
    ("cap-reaction", &["cap-search"]),
];

fn expected_caps<'a>(selected: impl IntoIterator<Item = &'a str>) -> BTreeSet<&'a str> {
    let mut pending: Vec<_> = selected.into_iter().collect();
    let mut caps = BTreeSet::new();
    while let Some(cap) = pending.pop() {
        if caps.insert(cap) {
            if let Some((_, dependencies)) = PREREQUISITES.iter().find(|(name, _)| *name == cap) {
                pending.extend(dependencies.iter().copied());
            }
        }
    }
    caps
}

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
        .chain(["cap-stereoisomers", "cap-reaction", "cap-serialization"])
        .collect();
    let declared: BTreeSet<_> = manifest["features"]
        .as_table()
        .unwrap()
        .keys()
        .filter(|name| name.starts_with("cap-"))
        .map(String::as_str)
        .collect();
    // Public reaction bindings enable the same capability in the full bundle.
    assert_eq!(declared, all);
    let mut cases = vec![
        (true, vec![], all.clone()),
        (false, vec!["full"], all.clone()),
        (false, vec![], BTreeSet::new()),
        (
            false,
            FOUR_CAPS.to_vec(),
            expected_caps(FOUR_CAPS.iter().copied()),
        ),
    ];
    for &(bundle, caps) in BUNDLES {
        let core = if bundle == "bio" { &[][..] } else { CORE };
        cases.push((
            false,
            vec![bundle],
            expected_caps(caps.iter().chain(core).copied()),
        ));
    }
    for &cap in &all {
        cases.push((false, vec![cap], expected_caps([cap])));
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
        cases.push((false, selected.clone(), expected_caps(selected)));
    }
    cases.push((
        false,
        vec!["core", "bio"],
        expected_caps(CORE.iter().copied().chain(["cap-bio"])),
    ));
    assert_eq!(cases.len(), 54);
    assert!(!manifest["features"].as_table().unwrap().contains_key("io"));
    assert!(declared_caps(&manifest, false, &["core"]).contains("cap-serialization"));
    assert!(declared_caps(&manifest, false, &["core"]).contains("cap-batch"));
    for bundle in [
        "core",
        "depict",
        "descriptors",
        "fingerprints",
        "conformer",
        "stereoisomers",
        "inchi",
    ] {
        let caps = declared_caps(&manifest, false, &[bundle]);
        assert!(!caps.contains("cap-search"), "{bundle}");
        if bundle != "depict" {
            assert!(!caps.contains("cap-depict"), "{bundle}");
        }
    }
    for removed in ["forcefields", "serialization", "batch", "io"] {
        assert!(
            !manifest["features"]
                .as_table()
                .unwrap()
                .contains_key(removed)
        );
    }
    assert!(declared_caps(&manifest, false, &["reaction"]).contains("cap-search"));
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
        ("cap-serialization", &["cap-io", "cosmolkit-io/binary"][..]),
        (
            "cap-fingerprints",
            &["dep:cosmolkit-fingerprints", "dep:cosmolkit-core", "cap-io"][..],
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
    assert_eq!(io["features"]["default"].as_array().unwrap().len(), 4);
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
    for dependency in ["dep:cosmolkit-core"] {
        assert!(molecule.contains(dependency));
    }
    assert!(!molecule.contains("bio"));
    assert!(!molecule.contains("dep:cosmolkit-bio"));
    for dependency in [
        "dep:cosmolkit-search",
        "dep:cosmolkit-depict",
        "dep:musli",
        "dep:postcard",
    ] {
        assert!(!molecule.contains(dependency));
    }
}

#[test]
fn enabled_registry_rows_use_capability_names_not_bundle_names() {
    let enabled = [
        ("cap-alignment", cfg!(feature = "cap-alignment")),
        ("cap-batch", cfg!(feature = "cap-batch")),
        ("cap-reaction", cfg!(feature = "cap-reaction")),
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

#[test]
fn reaction_object_and_execution_are_gated_together() {
    for id in ["types.Reaction", "Reaction.run"] {
        let entry = cosmolkit::BINDING_CONTRACT
            .iter()
            .find(|entry| entry.semantic_id == id);
        assert_eq!(entry.is_some(), cfg!(feature = "cap-reaction"), "{id}");
    }
    assert_eq!(
        cosmolkit::MOLECULE_OPS.iter().any(|op| op.method == "run"),
        cfg!(feature = "cap-reaction")
    );
}

#[cfg(feature = "cap-reaction")]
#[test]
fn reaction_selection_exposes_search_and_smiles_without_extra_features() {
    let _reaction = cosmolkit::Reaction::from_smirks("[C:1]>>[C:1]").unwrap();
    let query = cosmolkit::parse_smarts("[#6]").unwrap();
    let compiled = cosmolkit::compile_query(&query).unwrap();
    let molecule = cosmolkit::Molecule::from_smiles("CCO").unwrap();
    assert_eq!(molecule.num_atoms(), 3);
    assert_eq!(molecule.substruct_matches(&query).unwrap().len(), 2);
    let _ = compiled;
    for id in ["search.parse_smarts", "Molecule.from_smiles"] {
        assert!(
            cosmolkit::BINDING_CONTRACT
                .iter()
                .any(|entry| entry.semantic_id == id),
            "{id}"
        );
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
        with_uff_optimized_conformers, with_uff_optimized_conformers_with_params,
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

#[cfg(feature = "core")]
#[test]
fn core_parses_sdf_text_with_explicit_coordinate_policy() {
    let sdf = "carbon\n     RDKit          3D\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    1.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n$$$$\n";
    let record = cosmolkit::SdfRecord::from_sdf_with_params(
        sdf,
        &cosmolkit::SdfReadParams {
            coordinate_mode: cosmolkit::SdfCoordinateMode::Require3D,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(record.molecule().unwrap().num_atoms(), 1);
    assert_eq!(record.molecule().unwrap().conformers_3d().len(), 1);
}

#[cfg(all(feature = "core", not(target_arch = "wasm32")))]
#[test]
fn core_includes_native_archives_and_ordered_batch_without_optional_domains() {
    let molecule = cosmolkit::Molecule::from_smiles("CCO").unwrap();
    let bytes = molecule.to_binary().unwrap();
    let restored = cosmolkit::Molecule::from_binary(&bytes).unwrap();
    assert_eq!(restored.to_binary().unwrap(), bytes);
    assert_eq!(restored.to_smiles().unwrap(), molecule.to_smiles().unwrap());
    let batch = cosmolkit::MoleculeBatch::from_smiles_list(&["CCO".into(), "N".into()]).unwrap();
    assert_eq!(
        batch.to_smiles_list().unwrap(),
        vec![Some("CCO".into()), Some("N".into())]
    );
}

#[cfg(all(feature = "core", not(target_arch = "wasm32")))]
#[test]
fn io_coordinate_generation_requires_depict_without_silent_stereo_loss() {
    let molecule = cosmolkit::Molecule::from_smiles("N[C@@H](C)C(=O)O").unwrap();
    let before = molecule.to_binary().unwrap();
    let result = molecule.to_mol();
    if cfg!(feature = "cap-depict") {
        let restored = cosmolkit::Molecule::from_mol(&result.unwrap()).unwrap();
        assert_eq!(restored.to_smiles().unwrap(), molecule.to_smiles().unwrap());
    } else {
        assert!(matches!(
            result,
            Err(cosmolkit::MolecularIoError::MolWrite(
                cosmolkit::MolWriteError::MissingCapability("depict")
            ))
        ));
        let params = cosmolkit::MolBlockWriteParams {
            include_stereo: false,
            ..Default::default()
        };
        assert!(
            molecule
                .to_mol_with_params(&params)
                .unwrap()
                .contains("M  END")
        );
    }
    assert_eq!(molecule.to_binary().unwrap(), before);
}

#[cfg(feature = "core")]
#[test]
fn io_query_records_require_search_instead_of_lowering_to_concrete_atoms() {
    let mol = "query\n  COSMolKit\n\n  0  0  0  0  0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 [C,N] 0 0 0 0\nM  V30 END ATOM\nM  V30 END CTAB\nM  END\n$$$$\n";
    let result = cosmolkit::SdfRecord::from_sdf_with_params(
        mol,
        &cosmolkit::SdfReadParams {
            sanitize: false,
            remove_hs: false,
            ..Default::default()
        },
    );
    if cfg!(feature = "cap-search") {
        assert!(matches!(
            result.unwrap().graph(),
            cosmolkit::SdfGraph::Query(_)
        ));
    } else {
        assert!(matches!(
            result,
            Err(cosmolkit::SdfError::Read(
                cosmolkit::SdfReadError::MissingCapability("search")
            ))
        ));
    }
}

#[cfg(feature = "cap-bio")]
#[test]
fn bio_constructors_and_cid_selection_work() {
    let structure = cosmolkit::BioStructure::from_pdb("HEADER    TEST\n").unwrap();
    let selection = cosmolkit::BioSelection::from_cid("//*").unwrap();
    assert!(structure.selected_atom_ids(&selection).unwrap().is_empty());
}
