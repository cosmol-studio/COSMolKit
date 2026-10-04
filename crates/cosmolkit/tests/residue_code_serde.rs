//! BIO-SERDE1-28: facade-level Serde re-export proof for ResidueCode.
//!
//! Facade crates verify exports without repeating core suites (repository
//! organization policy §7): this file covers only the exact wire behavior
//! of a representative subset plus one invalid input, importing ONLY the
//! public `cosmolkit::ResidueCode` — never BIO internals. The complete
//! 368-tuple suite lives in the owning `cosmolkit-bio` crate.

use cosmolkit::ResidueCode;

/// Facade subset: exact serialize/deserialize for ALA, the digit-start
/// R0TD, UNK and the distinct UNKNOWN token, plus a representative
/// invalid string.
#[test]
fn residue_serde_facade() {
    for (code, wire) in [
        (ResidueCode::ALA, "ALA"),
        (ResidueCode::R0TD, "0TD"),
        (ResidueCode::UNK, "UNK"),
        (ResidueCode::UNKNOWN, "UNKNOWN"),
    ] {
        let json = serde_json::to_string(&code).unwrap();
        assert_eq!(json, format!("\"{wire}\""), "{wire}: exact serialize");
        let parsed: ResidueCode = serde_json::from_str(&json).unwrap();
        assert_eq!(parsed, code, "{wire}: roundtrip");
    }
    // Facade-visible strictness sample: the Rust-safe digit-prefixed
    // VARIANT name is not a wire token.
    assert!(serde_json::from_str::<ResidueCode>("\"R0TD\"").is_err());
}
