use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};

#[derive(Debug)]
struct ExactOccurrence {
    category: &'static str,
    path: &'static str,
    needle: &'static str,
    count: usize,
}

const EXACT_OCCURRENCES: &[ExactOccurrence] = &[
    ExactOccurrence {
        category: "parser",
        path: "src/search/smarts_parse.rs",
        needle: "fn parse_smarts(",
        count: 1,
    },
    ExactOccurrence {
        category: "parser",
        path: "src/search/smarts_parse.rs",
        needle: "struct SmartsParser<'a>",
        count: 1,
    },
    ExactOccurrence {
        category: "recursive string representation",
        path: "src/search/query.rs",
        needle: "RecursiveSmarts(RecursiveStructureQuery)",
        count: 1,
    },
];

const EXPECTED_SMARTS_FUNCTIONS: &[&str] = &[
    "src/chemistry/coordinates.rs::build_rdkit_template_query_molecule",
    "src/chemistry/coordinates.rs::expand_template_smarts_bonds",
    "src/chemistry/coordinates.rs::parse_rdkit_template_graph_model",
    "src/chemistry/coordinates.rs::rdkit_default_template_smarts",
    "src/chemistry/forcefield/crystalff/torsion_preferences.rs::expand_crystalff_smarts_bonds",
    "src/chemistry/forcefield/crystalff/torsion_preferences.rs::map_pattern_atom_indices",
    "src/chemistry/forcefield/crystalff/torsion_preferences.rs::smarts",
    "src/chemistry/tautomer.rs::smarts",
    "src/io/sdf.rs::parse_marvin_smarts_line",
    "src/io/sdf/postprocess.rs::process_smartsq",
    "src/notation/smiles_write.rs::canonicalize_fragment_stack_for_smarts",
    "src/properties/descriptors.rs::rdkit_cached_smarts_matcher",
    "src/properties/descriptors.rs::rdkit_count_smarts_matches",
    "src/properties/fingerprint.rs::default_feature_smarts",
    "src/properties/fingerprint.rs::from_smarts_patterns",
    "src/search/query.rs::from_smarts",
    "src/search/query.rs::source_smarts",
    "src/search/query_graph.rs::to_smarts",
    "src/search/query_graph.rs::to_cx_smarts",
    "src/search/query_graph.rs::fragment_to_smarts",
    "src/search/query_graph.rs::from_smarts",
    "src/search/query_graph.rs::fragment_to_cx_smarts",
    "src/search/smarts_parse.rs::atom_from_smarts",
    "src/search/smarts_parse.rs::bond_from_smarts",
    "src/search/smarts_parse.rs::materialize_smarts_atom_state",
    "src/search/smarts_parse.rs::mol_from_smarts",
    "src/search/smarts_parse.rs::parse_atom_token",
    "src/search/smarts_parse.rs::parse_smarts",
    "src/search/smarts_parse.rs::parse_smarts_chain",
    "src/search/smarts_parse.rs::parse_smarts_molecule",
    "src/search/smarts_parse.rs::parse_smarts_with_params",
    "src/search/smarts_parse.rs::preprocess_smarts",
    "src/search/smarts_parse.rs::smarts_atom_parse",
    "src/search/smarts_parse.rs::smarts_bond_parse",
    "src/search/smarts_parse.rs::smarts_parse_entry",
    "src/search/smarts_parse.rs::smarts_parse_helper",
    "src/search/smarts_parse.rs::unspecified_smarts_bond_query",
    "src/search/smarts_write.rs::combine_child_smarts",
    "src/search/smarts_write.rs::fragment_smarts_construct",
    "src/search/smarts_write.rs::get_atom_smarts",
    "src/search/smarts_write.rs::get_atom_smarts_simple",
    "src/search/smarts_write.rs::get_bond_smarts",
    "src/search/smarts_write.rs::get_bond_smarts_simple",
    "src/search/smarts_write.rs::get_non_query_atom_smarts",
    "src/search/smarts_write.rs::get_non_query_bond_smarts",
    "src/search/smarts_write.rs::get_recursive_structure_query_smarts",
    "src/search/smarts_write.rs::mol_fragment_to_cx_smarts",
    "src/search/smarts_write.rs::mol_fragment_to_smarts",
    "src/search/smarts_write.rs::mol_to_cx_smarts",
    "src/search/smarts_write.rs::mol_to_smarts",
    "src/search/smarts_write.rs::mol_to_smarts_internal",
    "src/search/smarts_write.rs::recurse_bond_smarts",
    "src/search/smarts_write.rs::recurse_get_smarts",
    "src/search/substruct.rs::recursive_smarts_root_matches",
];

fn core_root() -> PathBuf {
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

fn without_cfg_test_modules(source: &str) -> String {
    let lines = source.lines().collect::<Vec<_>>();
    let mut production = String::with_capacity(source.len());
    let mut line_idx = 0;

    while line_idx < lines.len() {
        let line = lines[line_idx];
        let next_is_test_module = line.trim() == "#[cfg(test)]"
            && lines
                .get(line_idx + 1)
                .is_some_and(|next| next.trim_start().starts_with("mod "));
        if !next_is_test_module {
            production.push_str(line);
            production.push('\n');
            line_idx += 1;
            continue;
        }

        line_idx += 1;
        let mut depth = 0_isize;
        let mut entered_module = false;
        while line_idx < lines.len() {
            let module_line = lines[line_idx];
            depth += module_line.matches('{').count() as isize;
            depth -= module_line.matches('}').count() as isize;
            entered_module |= module_line.contains('{');
            line_idx += 1;
            if entered_module && depth == 0 {
                break;
            }
        }
    }

    production
}

fn production_source(path: &Path) -> String {
    let source = fs::read_to_string(path)
        .unwrap_or_else(|error| panic!("failed to read {}: {error}", path.display()));
    without_cfg_test_modules(&source)
}

fn function_name(line: &str) -> Option<&str> {
    let line = line.trim_start();
    if line.starts_with("//") {
        return None;
    }
    let function = line.find("fn ")?;
    let prefix = line[..function].trim();
    if !prefix.is_empty()
        && !prefix.split_whitespace().all(|part| {
            matches!(
                part,
                "pub" | "pub(crate)" | "pub(super)" | "async" | "unsafe"
            )
        })
    {
        return None;
    }
    let name = &line[function + 3..];
    let end = name
        .find(|character: char| character == '(' || character == '<' || character.is_whitespace())
        .unwrap_or(name.len());
    (end != 0).then_some(&name[..end])
}

fn smarts_function_inventory(root: &Path) -> BTreeSet<String> {
    let source_root = root.join("src");
    let mut files = Vec::new();
    collect_rust_files(&source_root, &mut files);

    files
        .into_iter()
        .filter(|path| path.file_name().is_none_or(|name| name != "tests.rs"))
        .flat_map(|path| {
            let relative = path
                .strip_prefix(root)
                .expect("source path is below crate root")
                .to_string_lossy()
                .replace('\\', "/");
            production_source(&path)
                .lines()
                .filter_map(function_name)
                .filter(|name| {
                    name.contains("smarts")
                        || matches!(
                            *name,
                            "parse_rdkit_template_graph_model"
                                | "build_rdkit_template_query_molecule"
                                | "map_pattern_atom_indices"
                                | "parse_atom_token"
                        )
                })
                .map(|name| format!("{relative}::{name}"))
                .collect::<Vec<_>>()
        })
        .collect()
}

fn structural_hit_inventory(root: &Path, needle: &str) -> BTreeMap<String, usize> {
    let source_root = root.join("src");
    let mut files = Vec::new();
    collect_rust_files(&source_root, &mut files);

    files
        .into_iter()
        .filter_map(|path| {
            let count = production_source(&path).matches(needle).count();
            (count != 0).then(|| {
                (
                    path.strip_prefix(root)
                        .expect("source path is below crate root")
                        .to_string_lossy()
                        .replace('\\', "/"),
                    count,
                )
            })
        })
        .collect()
}

#[test]
fn smarts_duplicate_inventory_baseline() {
    let root = core_root();

    for expected in EXACT_OCCURRENCES {
        let source = production_source(&root.join(expected.path));
        assert_eq!(
            source.matches(expected.needle).count(),
            expected.count,
            "{} inventory changed for {} in {}; migrate an existing branch or update the explicit architecture baseline instead of adding a parallel SMARTS implementation",
            expected.category,
            expected.needle,
            expected.path,
        );
    }

    let expected_functions = EXPECTED_SMARTS_FUNCTIONS
        .iter()
        .map(|entry| (*entry).to_owned())
        .collect::<BTreeSet<_>>();
    assert_eq!(
        smarts_function_inventory(&root),
        expected_functions,
        "SMARTS production-function inventory changed; parser or decoder branches must remain within the canonical ownership baseline",
    );

    assert_eq!(
        structural_hit_inventory(&root, "struct SmartsParser<'a>"),
        BTreeMap::from([("src/search/smarts_parse.rs".to_owned(), 1)]),
        "SMARTS parser-engine inventory changed",
    );
    assert_eq!(
        structural_hit_inventory(&root, "build_query_molecule(inner)"),
        BTreeMap::new(),
        "alternate recursive SMARTS conversion appeared",
    );
    assert_eq!(
        structural_hit_inventory(&root, "struct ParserState<'a>"),
        BTreeMap::new(),
        "consumer-local SMARTS decoder inventory changed",
    );
}

#[test]
fn production_inventory_keeps_code_after_interleaved_test_modules() {
    let source = r#"
fn first_smarts_owner() {}

#[cfg(test)]
mod tests {
    fn test_smarts_helper() {
        assert_eq!("{query}", "{query}");
    }
}

fn second_smarts_owner() {}
"#;
    let production = without_cfg_test_modules(source);
    assert!(production.contains("fn first_smarts_owner()"));
    assert!(production.contains("fn second_smarts_owner()"));
    assert!(!production.contains("fn test_smarts_helper()"));
}

#[test]
fn smarts_single_core_architecture() {
    let root = core_root();
    let mut files = Vec::new();
    collect_rust_files(&root.join("src"), &mut files);
    let production = files
        .iter()
        .map(|path| production_source(path))
        .collect::<String>();

    assert!(!production.contains("SmartsMolecule"));
    assert!(!production.contains("build_query_molecule"));
    assert_eq!(production.matches("struct SmartsParser<'a>").count(), 1);
    assert!(!production.contains("struct ParserState<'a>"));
}

#[test]
fn smarts_existing_consumer_regression_matrix() {
    let root = core_root();
    let required = [
        (
            "src/properties/descriptors.rs",
            "smarts_consumer_descriptor_patterns",
        ),
        (
            "src/properties/fingerprint.rs",
            "smarts_consumer_fingerprint_patterns",
        ),
        (
            "src/chemistry/forcefield/torsion_query.rs",
            "smarts_consumer_torsion_query",
        ),
        (
            "src/chemistry/coordinates.rs",
            "smarts_consumer_coordinate_template_query",
        ),
        ("src/io/sdf.rs", "smarts_consumer_sdf_smartsq"),
    ];
    for (path, test_name) in required {
        let source = fs::read_to_string(root.join(path)).unwrap();
        assert!(
            source.contains(test_name),
            "missing SMARTS consumer regression {test_name}"
        );
    }
}

#[test]
fn smarts_canonical_module_baseline() {
    let root = core_root();
    let canonical = production_source(&root.join("src/search/smarts_parse.rs"));
    let query_model = production_source(&root.join("src/search/query.rs"));
    let facade = production_source(&root.join("src/lib.rs"));

    assert!(
        canonical.contains("This module is the sole canonical SMARTS parser/compiler owner."),
        "the canonical SMARTS parser/compiler owner must remain explicit",
    );
    assert!(
        query_model.contains("define the predicate vocabulary."),
        "the shared typed predicate vocabulary owner must remain explicit",
    );

    for (needle, expected) in [
        (
            "pub enum QueryNode<",
            BTreeMap::from([("src/search/query.rs".to_owned(), 1)]),
        ),
        (
            "pub enum AtomQueryPredicate",
            BTreeMap::from([("src/search/query.rs".to_owned(), 1)]),
        ),
        (
            "pub enum BondQueryPredicate",
            BTreeMap::from([("src/search/query.rs".to_owned(), 1)]),
        ),
        (
            "pub enum SmartsParseError",
            BTreeMap::from([("src/search/query.rs".to_owned(), 1)]),
        ),
        (
            "pub fn mol_from_smarts",
            BTreeMap::from([("src/search/smarts_parse.rs".to_owned(), 1)]),
        ),
    ] {
        assert_eq!(
            structural_hit_inventory(&root, needle),
            expected,
            "canonical SMARTS ownership changed for {needle}",
        );
    }

    assert_eq!(
        structural_hit_inventory(&root, "crate::search::query::parse_smarts("),
        BTreeMap::new(),
        "production callers must target the canonical SMARTS parser",
    );
    assert_eq!(
        structural_hit_inventory(&root, "use crate::search::query::parse_smarts"),
        BTreeMap::new(),
        "production callers must import the canonical SMARTS parser",
    );
    assert!(
        facade.contains("pub use search::{query, smarts_parse, substruct};")
            && facade.contains(
                "pub use search::query::{AtomQueryPredicate, BondQueryPredicate, QueryNode, SmartsParseError};",
            )
            && facade.contains("CompiledQuery")
            && facade.contains("QueryGraph"),
        "the public facade must expose the canonical parser module and first-class query graph model",
    );
}
