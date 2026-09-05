use std::{
    collections::BTreeMap,
    fs,
    path::{Path, PathBuf},
};

fn core_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn repository_root() -> PathBuf {
    core_root()
        .parent()
        .and_then(Path::parent)
        .expect("core crate is under repository/crates")
        .to_path_buf()
}

fn collect_rust_files(directory: &Path, files: &mut Vec<PathBuf>) {
    for entry in fs::read_dir(directory).expect("read source directory") {
        let path = entry.expect("source entry").path();
        if path.is_dir() {
            collect_rust_files(&path, files);
        } else if path.extension().is_some_and(|extension| extension == "rs") {
            files.push(path);
        }
    }
}

fn core_inventory(needle: &str) -> BTreeMap<String, usize> {
    let root = core_root();
    let mut files = Vec::new();
    collect_rust_files(&root.join("src"), &mut files);
    files.sort();
    files
        .into_iter()
        .filter_map(|path| {
            let source = fs::read_to_string(&path).expect("read Rust source");
            let count = source.matches(needle).count();
            (count != 0).then(|| {
                (
                    path.strip_prefix(&root)
                        .expect("source below crate root")
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
        core_inventory(needle),
        BTreeMap::from([(owner.to_owned(), 1)]),
        "Layered architecture inventory changed for {needle}; extend the canonical core instead of adding a parallel implementation",
    );
}

#[test]
fn layered_algorithm_helpers_have_one_production_owner_each() {
    for needle in [
        "fn enumerate_fingerprint_paths_for_root(",
        "pub(crate) fn enumerate_fingerprint_paths(",
        "pub(crate) fn hash_range(",
        "fn prepare_layered_fingerprint<'a>(",
        "fn layered_topology_hash(",
        "fn layered_bond_order_hash(",
        "fn layered_atom_type_hash(",
        "fn layered_aromaticity_hash(",
        "fn layered_ring_presence_hash(",
        "fn layered_min_ring_size_hash(",
        "fn project_layered_path(",
    ] {
        assert_owned_once(needle, "src/properties/fingerprint.rs");
    }
}

#[test]
fn layered_query_helpers_remain_owned_by_the_shared_query_core() {
    for needle in [
        "pub(crate) fn complex_atom_query_helper(",
        "pub(crate) fn is_complex_atom_query(",
        "pub(crate) fn is_complex_bond_query(",
        "pub(crate) fn is_atom_aromatic(",
        "pub(crate) fn query_is_bond_in_ring(",
        "pub(crate) fn query_bond_min_ring_size(",
    ] {
        assert_owned_once(needle, "src/search/query.rs");
    }
}

#[test]
fn layered_uses_the_shared_vector_and_one_scalar_core() {
    assert_owned_once("pub struct Fingerprint {", "src/properties/fingerprint.rs");
    assert_eq!(
        core_inventory("pub struct LayeredFingerprint {"),
        BTreeMap::new(),
        "Layered must use the shared Fingerprint vector",
    );
    assert_eq!(
        core_inventory("pub fn layered_fingerprint("),
        BTreeMap::from([
            ("src/model/molecule.rs".to_owned(), 1),
            ("src/properties/fingerprint.rs".to_owned(), 1),
        ]),
        "Layered must retain one algorithm core and one Molecule delegate",
    );
    assert_eq!(
        core_inventory("pub fn layered_fingerprint_with_output("),
        BTreeMap::from([
            ("src/model/molecule.rs".to_owned(), 1),
            ("src/properties/fingerprint.rs".to_owned(), 1),
        ]),
        "Layered output must retain one algorithm core and one Molecule delegate",
    );
}

#[test]
fn layered_batch_and_python_surfaces_delegate_without_algorithm_copies() {
    let batch =
        fs::read_to_string(core_root().join("src/properties/batch.rs")).expect("read batch source");
    assert_eq!(
        batch
            .matches("fn layered_fingerprint_list_with_runtime(")
            .count(),
        1
    );
    assert_eq!(
        batch
            .matches("fn layered_fingerprint_with_output_list_with_runtime(")
            .count(),
        1
    );
    assert_eq!(batch.matches(".layered_fingerprint(params)").count(), 1);
    assert_eq!(
        batch
            .matches(".layered_fingerprint_with_output(params)")
            .count(),
        1
    );
    assert!(!batch.contains("prepare_layered_fingerprint"));
    assert!(!batch.contains("project_layered_path"));

    let python = fs::read_to_string(repository_root().join("python/src/lib.rs"))
        .expect("read Python binding source");
    assert_eq!(python.matches("fn fingerprint_layered(").count(), 1);
    assert_eq!(
        python
            .matches("fn fingerprint_layered_with_output(")
            .count(),
        1
    );
    assert_eq!(python.matches(".layered_fingerprint(&params)").count(), 1);
    assert_eq!(
        python
            .matches(".layered_fingerprint_with_output(&params)")
            .count(),
        1
    );
    assert!(!python.contains("prepare_layered_fingerprint"));
    assert!(!python.contains("project_layered_path"));
}
