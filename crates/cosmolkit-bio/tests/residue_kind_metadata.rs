//! BIO-KIND1-28: frozen regressions for the restored legacy residue-kind
//! names and value traits. Expected values are pasted from the ROOT-frozen
//! 12-row table — never computed by invoking the methods under test.

use cosmolkit_bio::ResidueInfoKind;

/// Frozen 12-row table (variant, ordinal, stable name); immutable
/// expected data from the packet.
const KINDS: [(ResidueInfoKind, u8, &str); 12] = [
    (ResidueInfoKind::Unknown, 0, "UNKNOWN"),
    (ResidueInfoKind::Aa, 1, "AA"),
    (ResidueInfoKind::Aad, 2, "AAD"),
    (ResidueInfoKind::Paa, 3, "PAA"),
    (ResidueInfoKind::Maa, 4, "MAA"),
    (ResidueInfoKind::Rna, 5, "RNA"),
    (ResidueInfoKind::Dna, 6, "DNA"),
    (ResidueInfoKind::Buf, 7, "BUF"),
    (ResidueInfoKind::Hoh, 8, "HOH"),
    (ResidueInfoKind::Pyr, 9, "PYR"),
    (ResidueInfoKind::Ket, 10, "KET"),
    (ResidueInfoKind::Els, 11, "ELS"),
];

/// Frozen 24-call owner name product: twelve literal enum/name/ordinal rows
/// x2 repeats with stable copied input checked before/after each call; the
/// census increments only after invocation. Const/fn-pointer controls are
/// compile-time, not extra runtime product. No transmute or invalid
/// discriminants.
#[test]
fn residue_kind_name_all_rows() {
    // Signature + const-evaluation controls (compile-time only).
    let name_fn: fn(ResidueInfoKind) -> &'static str = ResidueInfoKind::name;
    const _CONST_PROOF: &str = ResidueInfoKind::Aa.name();
    assert_eq!(_CONST_PROOF, "AA");

    let mut calls = 0usize;
    for (kind, ordinal, name) in KINDS {
        for repeat in 0..2 {
            let label = format!("{name}#{repeat}");
            let before = kind;
            assert_eq!(kind as u8, ordinal, "{label}: ordinal");
            assert_eq!(name_fn(kind), name, "{label}: literal name");
            calls += 1;
            assert_eq!(kind, before, "{label}: input unchanged");
        }
    }
    assert_eq!(calls, 24, "exact 24-call census");
}

/// Frozen 48-call Display product: all 12 kinds x2 repeats x2 formats
/// (default `format!` and the exact frozen width path
/// `format!("{:>12}", kind)`). BOTH formats expect the SAME unpadded
/// frozen name, because the legacy `write_str` mapping owns the
/// formatting — a set width adds no padding. Each invocation has its OWN
/// fresh copied-value/ordinal checkpoints immediately before/after that
/// call; census after invocation.
#[test]
fn residue_kind_display_all_rows_both_formats_unpadded() {
    let mut calls = 0usize;
    for (kind, ordinal, name) in KINDS {
        for repeat in 0..2 {
            let label = format!("{name}#{repeat}");

            // Invocation 1 (default format): fresh checkpoints around THIS call.
            let before_default = kind;
            assert_eq!(kind as u8, ordinal, "{label}: ordinal (default)");
            let default = format!("{}", kind);
            calls += 1;
            assert_eq!(default, name, "{label}: default format");
            assert_eq!(kind, before_default, "{label}: input unchanged (default)");

            // Invocation 2 (exact frozen width path): fresh checkpoints
            // around THIS call. write_str ignores the width: the frozen
            // name stays unpadded.
            let before_width = kind;
            assert_eq!(kind as u8, ordinal, "{label}: ordinal (width)");
            let wide = format!("{:>12}", kind);
            calls += 1;
            assert_eq!(wide, name, "{label}: width format stays unpadded");
            assert_eq!(kind, before_width, "{label}: input unchanged (width)");
        }
    }
    assert_eq!(calls, 48, "exact 48-call Display census");
}

/// Frozen 24-call Serialize product + 12 derived-consumer codec calls.
/// Each direct `serde_json::to_string` yields the literal expected JSON
/// string and a `Value::String`; each derived wrapper (test-only
/// `#[derive(serde::Serialize)]`, one field) yields the exact
/// `{"kind":"..."}` literal. Copied enum/ordinal before/after every call;
/// the two censuses are counted separately.
#[test]
fn residue_kind_serialize_all_rows_and_derived_consumer() {
    #[derive(serde::Serialize)]
    struct KindWrapper {
        kind: ResidueInfoKind,
    }

    let mut direct_calls = 0usize;
    let mut derived_calls = 0usize;
    for (kind, ordinal, name) in KINDS {
        for repeat in 0..2 {
            let label = format!("{name}#{repeat}");
            let before = kind;
            assert_eq!(kind as u8, ordinal, "{label}: ordinal");

            let json =
                serde_json::to_string(&kind).unwrap_or_else(|error| panic!("{label}: {error}"));
            direct_calls += 1;
            let expected_json = format!("\"{name}\"");
            assert_eq!(json, expected_json, "{label}: literal JSON");
            let value: serde_json::Value = serde_json::from_str(&json).unwrap();
            assert!(
                value == serde_json::Value::String(name.to_string()),
                "{label}: Value::String proof"
            );
            assert_eq!(kind, before, "{label}: input unchanged");
        }
    }
    assert_eq!(direct_calls, 24, "exact 24-call direct census");

    for (kind, _, name) in KINDS {
        let label = format!("{name}: derived");
        let before = kind;
        let wrapper = KindWrapper { kind };
        let json =
            serde_json::to_string(&wrapper).unwrap_or_else(|error| panic!("{label}: {error}"));
        derived_calls += 1;
        let expected = format!("{{\"kind\":\"{name}\"}}");
        assert_eq!(json, expected, "{label}: exact literal document");
        assert_eq!(kind, before, "{label}: input unchanged");
    }
    assert_eq!(derived_calls, 12, "exact 12-call derived census");
}
