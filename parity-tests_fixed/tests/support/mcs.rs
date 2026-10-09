//! Fixed upstream and complete JNK1 pair regressions, scheduled by Cargo.
use cosmolkit_parity_tests_fixed::{
    Result,
    special_regression::{self, Snapshot},
};
use std::sync::OnceLock;

fn compare_case(key: &str, id: &str) {
    static UPSTREAM: OnceLock<Result<Snapshot>> = OnceLock::new();
    static JNK1: OnceLock<Result<Snapshot>> = OnceLock::new();
    let cache = match key {
        "mcs_upstream" => &UPSTREAM,
        "mcs_jnk1" => &JNK1,
        _ => panic!("unknown MCS regression {key}"),
    };
    let snapshot = cache
        .get_or_init(|| {
            special_regression::preflight(key, &cosmolkit_parity_tests_fixed::expected())
        })
        .as_ref()
        .unwrap();
    cosmolkit_parity_tests_fixed::mcs::compare_case(key, snapshot, id).unwrap();
}

include!(concat!(env!("OUT_DIR"), "/mcs_cases.rs"));
