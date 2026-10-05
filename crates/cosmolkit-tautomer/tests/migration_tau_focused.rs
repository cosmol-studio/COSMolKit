#[path = "support/oracle.rs"]
mod oracle;
#[test]
fn focused_source_oracle_all_branches() {
    oracle::validate_oracle(
        "../../testdata/tautomer/expected/rdkit/tautomer_focused/tautomer.jsonl",
        18,
        136,
        "focused",
    );
}
