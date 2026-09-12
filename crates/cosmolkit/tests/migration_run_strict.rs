const RUNTIME_INVARIANTS: bool = cfg!(feature = "runtime-invariants");
const OPERATION_CONTRACTS: bool = cfg!(feature = "op-contracts");
const STRICT_SELECTOR: bool = cfg!(feature = "op-contracts-strict");

#[test]
fn active_parent_feature_tuple_is_an_expected_configuration() {
    match (RUNTIME_INVARIANTS, OPERATION_CONTRACTS, STRICT_SELECTOR) {
        (false, false, false) => assert_eq!(active_configuration(), "default-non-strict"),
        (true, false, false) => assert_eq!(active_configuration(), "runtime-invariants-only"),
        (false, true, false) => assert_eq!(active_configuration(), "op-contracts-only"),
        (true, true, true) => assert_eq!(active_configuration(), "combined-strict"),
        unexpected => panic!("unexpected parent feature configuration: {unexpected:?}"),
    }
}

#[test]
fn combined_strict_selector_enables_both_parent_checks() {
    if STRICT_SELECTOR {
        assert!(RUNTIME_INVARIANTS);
        assert!(OPERATION_CONTRACTS);
    }

    assert_eq!(
        STRICT_SELECTOR,
        RUNTIME_INVARIANTS && OPERATION_CONTRACTS,
        "the four required invocations must not confuse a leaf selector with combined strict",
    );
}

#[test]
fn parent_manifest_declares_the_exact_strict_feature_composition_once() {
    let manifest = include_str!("../Cargo.toml");

    assert_eq!(manifest.matches("runtime-invariants = []").count(), 1);
    assert_eq!(manifest.matches("op-contracts = []").count(), 1);
    assert_eq!(
        manifest
            .matches("op-contracts-strict = [\"runtime-invariants\", \"op-contracts\"]")
            .count(),
        1,
    );
    assert_eq!(env!("CARGO_PKG_NAME"), "cosmolkit");
}

fn active_configuration() -> &'static str {
    match (RUNTIME_INVARIANTS, OPERATION_CONTRACTS, STRICT_SELECTOR) {
        (false, false, false) => "default-non-strict",
        (true, false, false) => "runtime-invariants-only",
        (false, true, false) => "op-contracts-only",
        (true, true, true) => "combined-strict",
        _ => "unexpected",
    }
}
