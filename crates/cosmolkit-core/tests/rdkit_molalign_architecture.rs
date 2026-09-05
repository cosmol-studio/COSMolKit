use std::{
    collections::BTreeMap,
    fs,
    path::{Path, PathBuf},
};

fn core_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
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

fn inventory(needle: &str) -> BTreeMap<String, usize> {
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

#[test]
fn molalign_has_one_private_pure_rust_alignment_kernel() {
    assert_eq!(
        inventory("pub(crate) fn align_points("),
        BTreeMap::from([("src/chemistry/numerics/alignment.rs".to_owned(), 1)])
    );
    let alignment =
        fs::read_to_string(core_root().join("src/chemistry/numerics/alignment.rs")).expect("read alignment kernel");
    let mol_align = fs::read_to_string(core_root().join("src/chemistry/mol_align.rs")).expect("read MolAlign");
    for source in [&alignment, &mol_align] {
        assert!(!source.contains("extern \"C\""));
        assert!(!source.contains("AlignmentBackend"));
        assert!(!source.contains("alignment_ffi"));
    }
}

#[test]
fn depiction_distgeom_molalign_and_moltransforms_delegate_to_the_shared_kernel() {
    let coordinates =
        fs::read_to_string(core_root().join("src/chemistry/coordinates.rs")).expect("read coordinates source");
    let distgeom = fs::read_to_string(core_root().join("src/chemistry/distgeom.rs")).expect("read distgeom source");
    let mol_align = fs::read_to_string(core_root().join("src/chemistry/mol_align.rs")).expect("read MolAlign");
    let mol_transforms = fs::read_to_string(core_root().join("src/chemistry/mol_transforms.rs"))
        .expect("read molecule transforms source");
    for (owner, source) in [
        ("coordinates", coordinates),
        ("distgeom", distgeom),
        ("MolAlign", mol_align),
        ("MolTransforms", mol_transforms),
    ] {
        assert!(
            source.contains("chemistry::numerics::alignment") || source.contains("chemistry::numerics::alignment::{"),
            "{owner} stopped delegating to the sole alignment kernel"
        );
    }
    for forbidden in [
        "fn rdkit_transform3d_identity",
        "fn rdkit_transform3d_mul",
        "fn rdkit_transform3d_transform_point",
        "fn rdkit_transform3d_set_translation",
        "fn rdkit_transform3d_set_rotation_from_quaternion",
        "fn rdkit_transform3d_reflect",
    ] {
        assert!(
            !inventory(forbidden)
                .keys()
                .any(|path| { path != "src/chemistry/numerics/alignment.rs" }),
            "{forbidden} reintroduced a duplicate 3D alignment transform helper"
        );
    }
}

#[test]
fn ordinary_molalign_surface_does_not_absorb_o3a() {
    let mol_align = fs::read_to_string(core_root().join("src/chemistry/mol_align.rs")).expect("read MolAlign");
    for excluded in ["O3A", "Open3DAlign", "CrippenO3A", "MMFFO3A"] {
        assert!(
            !mol_align.contains(excluded),
            "{excluded} belongs to a separate capability boundary"
        );
    }
}
