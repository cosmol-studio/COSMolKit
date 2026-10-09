#[cfg(any(parity_special_mcs_upstream, parity_special_mcs_jnk1))]
#[path = "support/mcs.rs"]
mod mcs;

#[cfg(parity_special_structure_tags)]
#[path = "support/structure_tags.rs"]
mod structure_tags;

#[cfg(parity_special_molalign_focused)]
#[test]
fn molalign_focused() {
    let snapshot = cosmolkit_parity_tests_fixed::special_regression::preflight(
        "molalign_focused",
        &cosmolkit_parity_tests_fixed::expected(),
    )
    .unwrap();
    cosmolkit_parity_tests_fixed::molalign::compare_rows(&snapshot.rows);
}

#[cfg(any(
    parity_special_tautomer_long_conjugated,
    parity_special_tautomer_focused
))]
#[allow(dead_code)]
#[path = "support/tautomer.rs"]
mod tautomer_oracle;

#[cfg(parity_special_tautomer_long_conjugated)]
#[test]
fn tautomer_long_conjugated() {
    let snapshot = cosmolkit_parity_tests_fixed::special_regression::preflight(
        "tautomer_long_conjugated",
        &cosmolkit_parity_tests_fixed::expected(),
    )
    .unwrap();
    let directory = cosmolkit_parity_tests_fixed::directory().join("reports");
    std::fs::create_dir_all(&directory).unwrap();
    let report = tempfile::Builder::new()
        .prefix("tautomer-")
        .tempdir_in(directory)
        .unwrap()
        .keep();
    tautomer_oracle::compare_rows(&snapshot.rows, 2, &report, "long-conjugated");
}

#[cfg(parity_special_tautomer_focused)]
#[test]
fn tautomer_focused() {
    let snapshot = cosmolkit_parity_tests_fixed::special_regression::preflight(
        "tautomer_focused",
        &cosmolkit_parity_tests_fixed::expected(),
    )
    .unwrap();
    assert_eq!(snapshot.rows.len(), 18);
    let directory = cosmolkit_parity_tests_fixed::directory().join("reports");
    std::fs::create_dir_all(&directory).unwrap();
    let report = tempfile::Builder::new()
        .prefix("tautomer-focused-")
        .tempdir_in(directory)
        .unwrap()
        .keep();
    tautomer_oracle::compare_rows(&snapshot.rows, 136, &report, "focused");
}
#[cfg(parity_special_bio_mmcif_switches)]
#[test]
fn bio_mmcif_switches() {
    let snapshot = cosmolkit_parity_tests_fixed::special_regression::preflight(
        "bio_mmcif_switches",
        &cosmolkit_parity_tests_fixed::expected(),
    )
    .unwrap();
    cosmolkit_parity_tests_fixed::bio_mmcif::compare(&snapshot).unwrap();
}
#[cfg(parity_special_forcefield_optimizers)]
#[test]
fn forcefield_optimizers() {
    let snapshot = cosmolkit_parity_tests_fixed::special_regression::preflight(
        "forcefield_optimizers",
        &cosmolkit_parity_tests_fixed::expected(),
    )
    .unwrap();
    cosmolkit_parity_tests_fixed::forcefield_regression::compare_optimizers(&snapshot);
}

#[cfg(parity_special_mmff_builtin)]
#[test]
fn mmff_builtin() {
    let snapshot = cosmolkit_parity_tests_fixed::special_regression::preflight(
        "mmff_builtin",
        &cosmolkit_parity_tests_fixed::expected(),
    )
    .unwrap();
    cosmolkit_parity_tests_fixed::forcefield_regression::compare_builtin(&snapshot);
}
