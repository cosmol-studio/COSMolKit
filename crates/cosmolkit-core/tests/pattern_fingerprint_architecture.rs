use std::fs;
use std::path::{Path, PathBuf};

const PATTERN_RELATIVE: &str = "src/properties/fingerprint/pattern.rs";
const PATTERN_SIGNATURES: [&str; 13] = [
    "[*]~[*]",
    "[*]~[*]~[*]",
    "[R]~1~[R]~[R]~1",
    "[*]~[*](~[*])~[*]",
    "[R]~1[R]~[R]~[R]~1",
    "[*]~[*]~[*](~[*])~[*]",
    "[R]~1~[R]~[R]~[R]~[R]~1",
    "[R]~1~[R]~[R]~[R]~[R]~[R]~1",
    "[R](@[R])(@[R])~[R]~[R](@[R])(@[R])",
    "[R](@[R])(@[R])~[R]@[R]~[R](@[R])(@[R])",
    "[*]~[R](@[R])@[R](@[R])~[*]",
    "[*]~[R](@[R])@[R]@[R](@[R])~[*]",
    "[*]",
];

fn rust_sources(root: &Path) -> Vec<PathBuf> {
    fn visit(directory: &Path, paths: &mut Vec<PathBuf>) {
        for entry in fs::read_dir(directory).expect("read Rust source directory") {
            let path = entry.expect("source directory entry").path();
            if path.is_dir() {
                visit(&path, paths);
            } else if path.extension().is_some_and(|extension| extension == "rs") {
                paths.push(path);
            }
        }
    }

    let mut paths = Vec::new();
    visit(root, &mut paths);
    paths.sort();
    paths
}

fn source(path: &Path) -> String {
    fs::read_to_string(path).unwrap_or_else(|error| panic!("read {}: {error}", path.display()))
}

fn occurrences(text: &str, needle: &str) -> usize {
    text.match_indices(needle).count()
}

#[test]
fn pattern_table_and_compilation_have_one_owner() {
    let crate_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let source_root = crate_root.join("src");
    let pattern_path = crate_root.join(PATTERN_RELATIVE);
    let pattern_source = source(&pattern_path);

    assert_eq!(
        occurrences(&pattern_source, "const PATTERN_FINGERPRINT_SMARTS:"),
        1,
        "Pattern must retain one fixed source table"
    );
    assert_eq!(
        occurrences(&pattern_source, "static CACHE: OnceLock"),
        1,
        "Pattern must retain one process-lifetime compiled table"
    );
    assert_eq!(
        occurrences(&pattern_source, "SsMatcher::try_new(pattern)"),
        1,
        "the fixed Pattern table must have one compilation site"
    );

    for path in rust_sources(&source_root) {
        if path == pattern_path {
            continue;
        }
        let text = source(&path);
        let signature_count = PATTERN_SIGNATURES
            .iter()
            .filter(|signature| text.contains(**signature))
            .count();
        assert!(
            signature_count < 3,
            "{} contains a second Pattern-like SMARTS table ({signature_count} signatures)",
            path.display()
        );
        assert!(
            !text.contains("PATTERN_FINGERPRINT_SMARTS"),
            "{} reaches into the private Pattern table",
            path.display()
        );
        assert!(
            !text.contains("compiled_pattern_fingerprint_queries"),
            "{} adds another Pattern compilation/cache caller",
            path.display()
        );
    }
}

#[test]
fn pattern_uses_the_shared_query_match_hash_and_vector_cores() {
    let crate_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let pattern_source = source(&crate_root.join(PATTERN_RELATIVE));
    let fingerprint_source = source(&crate_root.join("src/properties/fingerprint.rs"));

    for (definition, expected) in [
        ("fn is_pattern_complex_query(", 1),
        ("fn is_tautomer_bond_query(", 1),
        ("fn update_pattern_fingerprint_impl(", 1),
    ] {
        assert_eq!(
            rust_sources(&crate_root.join("src"))
                .iter()
                .map(|path| occurrences(&source(path), definition))
                .sum::<usize>(),
            expected,
            "Pattern architecture requires exactly {expected} `{definition}` definition"
        );
    }

    assert!(pattern_source.contains("is_complex_atom_query(atom)"));
    assert_eq!(
        occurrences(&pattern_source, "try_get_substruct_matches_with_params_and_context("),
        1,
        "Pattern must have one match call through the shared VF2 core"
    );
    assert!(!pattern_source.contains("try_get_substruct_matches_with_params("));
    assert!(!pattern_source.contains("get_substruct_matches_with_params("));
    assert!(!pattern_source.contains("fn hash_combine("));
    assert!(!pattern_source.contains("struct PatternFingerprint {"));
    assert_eq!(
        occurrences(&fingerprint_source, "pub struct Fingerprint {"),
        1,
        "the crate must retain one explicit bit-vector representation"
    );
    assert!(pattern_source.contains("hash_combine(&mut match_index"));
    assert!(pattern_source.contains("Fingerprint::zeroed(params.n_bits)"));
    assert!(pattern_source.contains("fingerprint.set_bit("));
}

#[test]
fn pattern_public_scalar_and_batch_surfaces_are_thin_delegates() {
    let crate_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let molecule_source = source(&crate_root.join("src/model/molecule.rs"));
    let batch_source = source(&crate_root.join("src/properties/batch.rs"));
    let python_source = source(&crate_root.join("../../python/src/lib.rs"));

    assert_eq!(
        occurrences(
            &molecule_source,
            "crate::fingerprint::pattern_fingerprint(self, params)"
        ),
        1,
        "Molecule must delegate once to the sole Pattern facade"
    );
    assert_eq!(
        occurrences(&batch_source, ".pattern_fingerprint(params)"),
        1,
        "the Rust batch path must map the scalar Molecule method"
    );
    assert_eq!(
        occurrences(&python_source, ".pattern_fingerprint(&params)"),
        1,
        "the Python scalar path must delegate to the Rust Molecule method"
    );
    assert_eq!(
        occurrences(&python_source, ".pattern_fingerprint_list_with_options("),
        1,
        "the Python batch path must delegate to the Rust ordered batch method"
    );
    for source in [&molecule_source, &batch_source, &python_source] {
        assert!(!source.contains("PATTERN_FINGERPRINT_SMARTS"));
        assert!(!source.contains("SsMatcher::try_new"));
        assert!(!source.contains("try_get_substruct_matches"));
        assert!(!source.contains("fn hash_combine"));
    }
}
