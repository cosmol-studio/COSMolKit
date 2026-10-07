#[cfg(parity_special_structure_tags)]
#[path = "support/structure_tags.rs"]
mod structure_tags;

#[cfg(parity_special_tautomer_long_conjugated)]
#[allow(dead_code)]
#[path = "../../crates/cosmolkit-tautomer/tests/support/oracle.rs"]
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
