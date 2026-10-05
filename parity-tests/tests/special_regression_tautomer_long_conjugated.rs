//! Fixed source regression, not expansion over the SMILES corpus.
//! Reuse the original detached comparison boundary and every compared field.
#[allow(dead_code)]
#[path = "../../crates/cosmolkit-tautomer/tests/support/oracle.rs"]
mod oracle;

#[test]
fn original_case1399_default_and_v1_keep_every_ordered_candidate_state() {
    let data = std::env::var_os("PARITY_DATA")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| cosmolkit_parity_tests::root().join("target/parity-tests"));
    let data = cosmolkit_parity_tests::root().join(data);
    // All identities, input/parameter bindings and complete branch tables are
    // checked before catalog construction or any CK operation.
    let snapshot =
        cosmolkit_parity_tests::special_regression::preflight("tautomer_long_conjugated", &data)
            .expect("complete special regression preflight before CK calls");
    assert_eq!(snapshot.rows.len(), 1);
    for branch in ["default", "v1"] {
        let expected = &snapshot.rows[0]["branches"][branch];
        assert_eq!(expected["ordered_smiles"].as_array().unwrap().len(), 42);
        assert_eq!(expected["molecule_states"].as_array().unwrap().len(), 42);
        assert_eq!(expected["scores"].as_array().unwrap().len(), 42);
    }
    let report = tempfile::Builder::new()
        .prefix("special-run-")
        .tempdir_in(&data)
        .expect("actual comparison report directory")
        .keep();
    oracle::compare_rows(&snapshot.rows, 2, &report, "long-conjugated");
    println!(
        "tautomer_long_conjugated: 1 fixed case / 2 branches matched; {}",
        report.display()
    );
}
