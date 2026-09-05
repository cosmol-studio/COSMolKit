use std::{
    collections::BTreeMap,
    fs,
    path::{Path, PathBuf},
};

fn crate_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn collect_rust_files(directory: &Path, files: &mut Vec<PathBuf>) {
    let mut entries = fs::read_dir(directory)
        .unwrap_or_else(|error| panic!("failed to read {}: {error}", directory.display()))
        .collect::<Result<Vec<_>, _>>()
        .unwrap_or_else(|error| panic!("failed to enumerate {}: {error}", directory.display()));
    entries.sort_by_key(|entry| entry.path());

    for entry in entries {
        let path = entry.path();
        if path.is_dir() {
            collect_rust_files(&path, files);
        } else if path.extension().is_some_and(|extension| extension == "rs") {
            files.push(path);
        }
    }
}

fn production_source(path: &Path) -> String {
    let source = fs::read_to_string(path).unwrap_or_else(|error| panic!("failed to read {}: {error}", path.display()));
    source
        .find("#[cfg(test)]\nmod ")
        .map_or(source.as_str(), |test_module| &source[..test_module])
        .to_owned()
}

fn source_inventory(needle: &str) -> BTreeMap<String, usize> {
    let root = crate_root();
    let mut files = Vec::new();
    collect_rust_files(&root.join("src"), &mut files);

    files
        .into_iter()
        .filter_map(|path| {
            let count = production_source(&path).matches(needle).count();
            (count != 0).then(|| {
                (
                    path.strip_prefix(&root)
                        .expect("source path is below the crate root")
                        .to_string_lossy()
                        .replace('\\', "/"),
                    count,
                )
            })
        })
        .collect()
}

fn assert_owned_once(needle: &str, owner: &str) {
    assert_eq!(
        source_inventory(needle),
        BTreeMap::from([(owner.to_owned(), 1)]),
        "stereoisomer-enumeration ownership changed for `{needle}`; extend the canonical core instead of adding a parallel implementation",
    );
}

#[test]
fn potential_stereo_flippers_and_configurations_have_one_owner_each() {
    assert_owned_once(
        "pub(crate) fn find_potential_stereo_in_workspace(",
        "src/chemistry/potential_stereo.rs",
    );
    assert_owned_once("pub(crate) enum StereoFlipper {", "src/chemistry/stereo_enumerate.rs");
    assert_owned_once(
        "pub(crate) fn select_stereo_flippers(",
        "src/chemistry/stereo_enumerate.rs",
    );
    assert_owned_once("enum ConfigurationSource {", "src/chemistry/stereo_enumerate.rs");
    assert_owned_once("pub struct StereoisomerIterator {", "src/chemistry/stereo_enumerate.rs");
}

#[test]
fn public_surfaces_are_value_style_delegates_to_the_lazy_core() {
    let root = crate_root();
    let enumeration = production_source(&root.join("src/chemistry/stereo_enumerate.rs"));
    let molecule = production_source(&root.join("src/model/molecule.rs"));
    let facade = production_source(&root.join("src/lib.rs"));

    assert_eq!(enumeration.matches("pub fn enumerate_stereoisomers(").count(), 1,);
    assert_eq!(
        enumeration
            .matches("StereoisomerIterator::new(molecule, options)")
            .count(),
        1,
    );
    assert_eq!(
        molecule
            .matches("crate::stereo_enumerate::enumerate_stereoisomers(self, options)")
            .count(),
        1,
    );
    assert_eq!(
        molecule
            .matches("crate::stereo_enumerate::stereoisomer_count(self, options)")
            .count(),
        1,
    );
    assert!(facade.contains("EnumerationError, StereoisomerIterator, StereoisomerOptions, enumerate_stereoisomers,"));

    assert!(!enumeration.contains("pub fn enumerate_stereoisomers(\n    molecule: &mut Molecule"));
    assert!(!enumeration.contains("pub fn stereoisomer_count(\n    molecule: &mut Molecule"));
}

#[test]
fn retired_enumeration_strategies_cannot_reenter_production() {
    let retired = [
        "find_tetrahedral_centers",
        "find_stereo_bonds",
        "is_valid_tetrahedral_config",
        "is_valid_double_bond_config",
        "generate_combinations",
        "build_tetrahedral_isomer",
        "build_double_bond_isomer",
        "EnumerationStrategy",
        "EnumerationParams",
        "enum_stereoisomers",
        "enum_double_bond_stereoisomers",
        "enum_all_stereoisomers",
        "count_double_bond_stereoisomers",
        "count_all_stereoisomers",
        "max_tries",
        "xorshift",
        "maximum supported is 20",
    ];

    for needle in retired {
        assert_eq!(
            source_inventory(needle),
            BTreeMap::new(),
            "retired stereoisomer-enumeration strategy `{needle}` reappeared",
        );
    }
}

#[test]
fn enumeration_source_functions_are_not_duplicated() {
    for (needle, owner) in [
        (
            "fn theoretical_configuration_count(",
            "src/chemistry/stereo_enumerate.rs",
        ),
        ("fn default_python_random_seed(", "src/chemistry/stereo_enumerate.rs"),
        (
            "fn clear_computed_props_preserving_rings(",
            "src/chemistry/stereo_enumerate.rs",
        ),
        ("fn canonical_isomeric_smiles(", "src/model/read_parts.rs"),
    ] {
        assert_owned_once(needle, owner);
    }
}
