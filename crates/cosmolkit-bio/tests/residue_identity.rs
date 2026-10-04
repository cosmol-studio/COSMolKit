//! Frozen raw-name ResidueIdentity products (BIO-IDENTITY1-36).
//!
//! Identity16 literal annex is source-derived (BIO-CLOSURE report annex),
//! not CK output. Construction census 4 routes x 16 = 64 real calls; wire
//! census 16 Serialize + 16 Deserialize round trips against independently
//! encoded raw JSON literals, six malformed-wire errors, and six
//! strict-Code boundary controls. Per-call mismatches are COLLECTED until
//! all product calls complete; censuses assert only after every call ran.

use cosmolkit_bio::ResidueCode;
use cosmolkit_bio::ResidueIdentity;
use serde::Deserialize;
use serde::Serialize;
use std::convert::Infallible;
use std::str::FromStr;

/// Frozen literal annex: (raw name, expected ordinal, tabulated).
const LITERAL16: [(&str, u16, bool); 16] = [
    ("ALA", 0, true),
    ("ala", 0, true),
    ("MSE", 17, true),
    ("mse", 17, true),
    ("HOH", 154, true),
    ("wat", 154, true),
    ("H2O", 154, true),
    ("UNK", 25, true),
    ("UNKNOWN", 367, false),
    ("unknown", 367, false),
    ("CUSTOM", 367, false),
    ("", 367, false),
    ("A", 327, true),
    ("AL", 343, true),
    ("al", 367, false),
    ("DA", 335, true),
];

#[test]
fn identity_construction_four_routes_exact_64() {
    let mut calls = 0usize;
    let mut discrepancies: Vec<String> = Vec::new();
    for (raw, ordinal, tabulated) in LITERAL16 {
        // Per-constructor caller-input baselines: fresh input bytes are
        // captured immediately BEFORE each actual constructor call and
        // compared immediately AFTER that SAME call, so each route is
        // bracketed independently instead of sharing one prefix check.
        let input_bytes = raw.as_bytes().to_vec();
        let via_new = ResidueIdentity::new(raw);
        calls += 1;
        if input_bytes != raw.as_bytes() {
            discrepancies.push(format!("new:{raw}: caller input mutated"));
        }
        let input_bytes = raw.as_bytes().to_vec();
        let via_ref: ResidueIdentity = ResidueIdentity::from(raw);
        calls += 1;
        if input_bytes != raw.as_bytes() {
            discrepancies.push(format!("from&:{raw}: caller input mutated"));
        }
        // Owned-route move proof: allocate the input ONCE, capture its
        // original buffer pointer BEFORE moving it into the constructor,
        // and compare the stored-name pointer AFTER. Nonempty only, since
        // an empty String has no stable distinct allocation. The owned
        // input is never cloned and no extra construction calls are added.
        let owned_input = raw.to_string();
        let owned_bytes = owned_input.as_bytes().to_vec();
        let owned_ptr = owned_input.as_ptr();
        let via_owned: ResidueIdentity = ResidueIdentity::from(owned_input);
        calls += 1;
        let owned_checkpoint = via_owned.clone();
        let stored_name = via_owned.name();
        if via_owned != owned_checkpoint {
            discrepancies.push(format!("fromString:{raw}: name() mutated identity"));
        }
        if stored_name.as_bytes() != owned_bytes {
            discrepancies.push(format!("fromString:{raw}: owned input bytes changed"));
        }
        if !raw.is_empty() && stored_name.as_ptr() != owned_ptr {
            discrepancies.push(format!("fromString:{raw}: owned buffer not retained"));
        }
        let input_bytes = raw.as_bytes().to_vec();
        let via_parse: Result<ResidueIdentity, Infallible> = raw.parse();
        calls += 1;
        if input_bytes != raw.as_bytes() {
            discrepancies.push(format!("parse:{raw}: caller input mutated"));
        }
        let via_parse = via_parse.unwrap_or_else(|never| match never {});

        // ALL FOUR results get their own frozen checks; the parsed route is
        // checked directly rather than only compared to another route.
        let parse_reference = via_parse.clone();
        for (route, id) in [
            ("new", via_new),
            ("from&", via_ref),
            ("fromString", via_owned),
            ("parse", via_parse.clone()),
        ] {
            let label = format!("{route}:{raw}");
            let checkpoint = id.clone();
            let name = id.name();
            if id != checkpoint {
                discrepancies.push(format!("{label}: name() mutated identity"));
            }
            if name != raw {
                discrepancies.push(format!("{label}: raw name {name} != {raw}"));
            }
            let checkpoint = id.clone();
            let actual_ordinal = u16::from(id.code().as_u16());
            if id != checkpoint {
                discrepancies.push(format!("{label}: code() mutated identity"));
            }
            if actual_ordinal != ordinal {
                discrepancies.push(format!("{label}: ordinal {actual_ordinal} != {ordinal}"));
            }
            let checkpoint = id.clone();
            let actual_tabulated = id.is_tabulated();
            if id != checkpoint {
                discrepancies.push(format!("{label}: is_tabulated() mutated identity"));
            }
            if actual_tabulated != tabulated {
                discrepancies.push(format!(
                    "{label}: tabulated {actual_tabulated} != {tabulated}"
                ));
            }
            let checkpoint = id.clone();
            let actual_info_ordinal = u16::from(id.info().code.as_u16());
            if id != checkpoint {
                discrepancies.push(format!("{label}: info() mutated identity"));
            }
            if actual_info_ordinal != ordinal {
                discrepancies.push(format!("{label}: info.code disagrees"));
            }
            let checkpoint = id.clone();
            let displayed = id.to_string();
            if id != checkpoint {
                discrepancies.push(format!("{label}: Display mutated identity"));
            }
            if displayed != raw {
                discrepancies.push(format!("{label}: Display {displayed} != {raw}"));
            }
            // Cross-route agreement remains an ADDITIONAL check; the
            // per-route frozen checks above no longer depend on it.
            if id != parse_reference {
                discrepancies.push(format!("{label}: route disagrees with FromStr"));
            }
            let before = id.clone();
            let serialized = serde_json::to_string(&id);
            if id != before {
                discrepancies.push(format!("{label}: serialization mutated identity"));
            }
            let serialized = serialized.expect("serialize identity");
            let expected_json = serde_json::to_string(raw).expect("encode raw literal");
            if serialized != expected_json {
                discrepancies.push(format!(
                    "{label}: wire {serialized} != raw literal {expected_json}"
                ));
            }
        }
    }
    assert_eq!(calls, 64, "exact 64-call construction census");
    assert!(
        discrepancies.is_empty(),
        "collected discrepancies after all 64 calls: {discrepancies:?}"
    );
}

#[test]
fn identity_from_string_retains_original_buffer_pointer() {
    let mut buffer = String::from("MSE");
    let buffer_ptr = buffer.as_ptr();
    let identity = ResidueIdentity::from(buffer);
    assert_eq!(identity.name().as_ptr(), buffer_ptr);
    assert_eq!(identity.name(), "MSE");
    assert_eq!(u16::from(identity.code().as_u16()), 17);
    let borrowed_input = "wat";
    let bytes_before = borrowed_input.as_bytes().to_vec();
    let from_borrowed = ResidueIdentity::from(borrowed_input);
    assert_eq!(bytes_before, borrowed_input.as_bytes());
    assert_eq!(from_borrowed.name(), borrowed_input);
    assert_eq!(u16::from(from_borrowed.code().as_u16()), 154);
}

#[test]
fn identity_raw_wire_roundtrips_exact_16() {
    let mut serialize_calls = 0usize;
    let mut deserialize_calls = 0usize;
    let mut discrepancies: Vec<String> = Vec::new();
    for (raw, ordinal, _tabulated) in LITERAL16 {
        let identity = ResidueIdentity::new(raw);
        // Fresh raw-input and full original-identity snapshots immediately
        // BEFORE the serializer; both are compared immediately AFTER that
        // same invocation, before any output checks.
        let raw_before = raw.as_bytes().to_vec();
        let identity_before = identity.clone();
        let wire = serde_json::to_string(&identity);
        serialize_calls += 1;
        if raw.as_bytes() != raw_before || identity != identity_before {
            discrepancies.push(format!("{raw}: serialize mutated raw identity"));
        }
        let wire = wire.expect("serialize identity");
        let expected_wire = serde_json::to_string(raw).expect("encode raw literal");
        if wire != expected_wire {
            discrepancies.push(format!("{raw}: wire {wire} != {expected_wire}"));
        }
        // Wire bytes and the original identity are captured immediately
        // BEFORE the deserializer; both are compared immediately AFTER the
        // actual result exists, BEFORE any output checks.
        let wire_bytes = wire.as_bytes().to_vec();
        let original = identity.clone();
        let round = ResidueIdentity::deserialize(&mut serde_json::Deserializer::from_str(&wire));
        deserialize_calls += 1;
        if wire_bytes != wire.as_bytes() {
            discrepancies.push(format!("{raw}: deserialize consumed the wire literal"));
        }
        if identity != original {
            discrepancies.push(format!("{raw}: deserialize mutated the source identity"));
        }
        let round = round.expect("deserialize identity");
        if round != identity {
            discrepancies.push(format!("{raw}: roundtrip identity differs"));
        }
        if round.name() != raw {
            discrepancies.push(format!(
                "{raw}: roundtrip raw name {} differs",
                round.name()
            ));
        }
        if u16::from(round.code().as_u16()) != ordinal {
            discrepancies.push(format!("{raw}: roundtrip ordinal differs"));
        }
        let mut replaced = round.clone();
        replaced = ResidueIdentity::new("HOH");
        if round.name() != raw || u16::from(round.code().as_u16()) != ordinal {
            discrepancies.push(format!("{raw}: replacement mutated the sibling clone"));
        }
        if replaced.name() != "HOH" || u16::from(replaced.code().as_u16()) != 154 {
            discrepancies.push(format!("{raw}: replacement value incorrect"));
        }
    }
    assert_eq!(serialize_calls, 16, "exact 16 Serialize calls");
    assert_eq!(deserialize_calls, 16, "exact 16 actual Deserialize calls");
    assert!(
        discrepancies.is_empty(),
        "collected discrepancies after all wire calls: {discrepancies:?}"
    );
}

#[test]
fn identity_malformed_wire_six_errors() {
    // Frozen six-product malformed literals (annex order): the first is
    // the authorized literal `null` (a JSON null that is not a string),
    // restored per ROOT; the other five are unchanged.
    let malformed = ["null", "true", "1", "[]", "{}", "{\"name\":\"ALA\"}"];
    let mut errors = 0usize;
    let mut discrepancies: Vec<String> = Vec::new();
    for payload in malformed {
        let result: Result<ResidueIdentity, _> = serde_json::from_str(payload);
        match result {
            Ok(value) => {
                discrepancies.push(format!(
                    "{payload}: unexpectedly accepted as {}",
                    value.name()
                ));
            }
            Err(_) => errors += 1,
        }
    }
    assert_eq!(errors, 6, "exactly six malformed-wire deserialize errors");
    assert!(
        discrepancies.is_empty(),
        "malformed wire must never produce an identity: {discrepancies:?}"
    );
}

/// Supplementary control (separately named and counted, NOT part of the
/// frozen six-product): the previous truncated-JSON literal retained as
/// history so the correction does not erase what ran before.
#[test]
fn identity_malformed_wire_supplementary_truncated_json_rejected() {
    let result: Result<ResidueIdentity, _> = serde_json::from_str("[null,true,1");
    let error = result.expect_err("truncated JSON must be rejected");
    let message = error.to_string();
    assert!(
        message.contains("EOF") || message.contains("expected") || message.contains("trailing"),
        "unexpected truncated-JSON message: {message}"
    );
}

#[test]
fn identity_accepts_what_strict_code_rejects_six() {
    let controls = ["wat", "ala", "CUSTOM", "", "al", "unknown"];
    let mut identity_calls = 0usize;
    let mut code_rejections = 0usize;
    let mut discrepancies: Vec<String> = Vec::new();
    for raw in controls {
        let identity = ResidueIdentity::new(raw);
        identity_calls += 1;
        if identity.name() != raw {
            discrepancies.push(format!("{raw}: identity lost raw name"));
        }
        let wire = serde_json::to_string(raw).expect("raw literal");
        let strict_code: Result<ResidueCode, _> = serde_json::from_str(&wire);
        match strict_code {
            Ok(code) => {
                discrepancies.push(format!("{raw}: strict Code unexpectedly accepted {code:?}"));
            }
            Err(_) => code_rejections += 1,
        }
        let round: ResidueIdentity = serde_json::from_str(&wire).expect("identity wire");
        if round.name() != raw {
            discrepancies.push(format!("{raw}: identity roundtrip lost raw name"));
        }
    }
    assert_eq!(identity_calls, 6);
    assert_eq!(
        code_rejections, 6,
        "all six strict Code rejections retained"
    );
    assert!(
        discrepancies.is_empty(),
        "raw-identity/strict-Code separation violated: {discrepancies:?}"
    );
}
