use std::{
    collections::BTreeMap,
    fs,
    path::{Path, PathBuf},
};
fn workspace() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}
fn collect(directory: &Path, files: &mut Vec<PathBuf>) {
    for e in fs::read_dir(directory).unwrap() {
        let p = e.unwrap().path();
        if p.is_dir() {
            collect(&p, files)
        } else if p.extension().is_some_and(|x| x == "rs") {
            files.push(p)
        }
    }
}
fn inventory(needle: &str) -> BTreeMap<String, usize> {
    let root = workspace();
    let mut files = vec![];
    for name in ["cosmolkit-core", "cosmolkit-alignment", "cosmolkit"] {
        collect(&root.join("crates").join(name).join("src"), &mut files)
    }
    files
        .into_iter()
        .filter_map(|p| {
            let s = fs::read_to_string(&p).unwrap();
            let n = s.matches(needle).count();
            (n > 0).then(|| {
                (
                    p.strip_prefix(&root).unwrap().to_string_lossy().to_string(),
                    n,
                )
            })
        })
        .collect()
}
#[test]
fn molalign_has_one_private_pure_rust_alignment_kernel() {
    // Split ownership requires a narrow detached CORE export for domain reuse.
    assert_eq!(
        inventory("pub fn align_points("),
        BTreeMap::from([("crates/cosmolkit-core/src/alignment.rs".into(), 1)])
    );
    for name in [
        "cosmolkit-core/src/alignment.rs",
        "cosmolkit-alignment/src/lib.rs",
    ] {
        let s = fs::read_to_string(workspace().join("crates").join(name)).unwrap();
        for excluded in ["extern \"C\"", "AlignmentBackend", "alignment_ffi"] {
            assert!(!s.contains(excluded))
        }
    }
}
#[test]
fn depiction_distgeom_molalign_and_moltransforms_delegate_to_the_shared_kernel() {
    // The original co-located coordinates/distgeom files belong to other split
    // owners. This lane retains the sole-kernel/no-duplicate condition and tests
    // actual ALIGNMENT -> CORE, facade -> ALIGNMENT, and CORE transform reuse.
    let root = workspace();
    let owner = fs::read_to_string(root.join("crates/cosmolkit-alignment/src/lib.rs")).unwrap();
    assert!(owner.contains("align_points("));
    assert!(owner.contains("cosmolkit_core::{"));
    let facade = fs::read_to_string(root.join("crates/cosmolkit/src/alignment.rs")).unwrap();
    assert!(facade.contains("cosmolkit_alignment::AlignmentInput"));
    assert!(!facade.contains("align_points("));
    let kernel = fs::read_to_string(root.join("crates/cosmolkit-core/src/alignment.rs")).unwrap();
    assert!(
        kernel.contains("crate::transforms::Transform3D") || kernel.contains("crate::Transform3D")
    );
    for forbidden in [
        "fn rdkit_transform3d_identity",
        "fn rdkit_transform3d_mul",
        "fn rdkit_transform3d_transform_point",
        "fn rdkit_transform3d_set_translation",
        "fn rdkit_transform3d_set_rotation_from_quaternion",
        "fn rdkit_transform3d_reflect",
    ] {
        assert!(
            inventory(forbidden).is_empty(),
            "duplicate helper {forbidden}"
        );
    }
}
#[test]
fn ordinary_molalign_surface_does_not_absorb_o3a() {
    let s = fs::read_to_string(workspace().join("crates/cosmolkit-alignment/src/lib.rs")).unwrap();
    for excluded in ["O3A", "Open3DAlign", "CrippenO3A", "MMFFO3A"] {
        assert!(!s.contains(excluded))
    }
}
