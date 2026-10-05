#[path = "support/oracle.rs"]
mod oracle;
#[test]
fn original_case1399_default_and_v1_keep_every_ordered_candidate_state() {
    oracle::validate_oracle(
        "../../testdata/tautomer/expected/rdkit/regressions/long_conjugated.jsonl",
        1,
        2,
        "long-conjugated",
    );
}
