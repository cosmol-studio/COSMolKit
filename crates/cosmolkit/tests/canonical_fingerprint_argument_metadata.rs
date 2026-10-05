//! Proposal: pinned RDKit351f argument metadata and immutable public projection.
#![cfg(feature = "cap-fingerprints")]
use cosmolkit::{AtomPairParams, FingerprintJsonError, MorganParams, TopologicalTorsionParams};
use serde_json::{Value, json};

#[test]
fn default_source_info_and_complete_quoted_leaf_metadata() {
    let m = MorganParams::default();
    let a = AtomPairParams::default();
    let t = TopologicalTorsionParams::default();
    assert_eq!(
        m.info_string(),
        "MorganArguments onlyNonzeroInvariants=0 radius=3"
    );
    assert_eq!(
        a.info_string(),
        "AtomPairArguments use2D=1 minDistance=1 maxDistance=30"
    );
    assert_eq!(
        t.info_string(),
        "TopologicalTorsionArguments torsionAtomCount=4 onlyShortestPaths=0"
    );
    for (actual, expected) in [
        (
            m.to_json(),
            json!({"type":"MorganArguments","onlyNonzeroInvariants":"false","radius":"3","countSimulation":"false","fpSize":"2048","numBitsPerFeature":"1","includeChirality":"false","countBounds":["1","2","4","8"]}),
        ),
        (
            a.to_json(),
            json!({"type":"AtomPairArguments","use2D":"true","minDistance":"1","maxDistance":"30","countSimulation":"true","fpSize":"2048","numBitsPerFeature":"1","includeChirality":"false","countBounds":["1","2","4","8"]}),
        ),
        (
            t.to_json(),
            json!({"type":"TopologicalTorsionArguments","torsionAtomCount":"4","onlyShortestPaths":"false","countSimulation":"true","fpSize":"2048","numBitsPerFeature":"1","includeChirality":"false","countBounds":["1","2","4","8"]}),
        ),
    ] {
        assert_eq!(serde_json::from_str::<Value>(&actual).unwrap(), expected);
    }
    assert_eq!(m.with_json(&m.to_json()).unwrap(), m);
    assert_eq!(a.with_json(&a.to_json()).unwrap(), a);
    assert_eq!(t.with_json(&t.to_json()).unwrap(), t);
}
#[test]
fn omitted_morgan_policies_preserve_receiver_and_do_not_appear_in_source_json() {
    let original = MorganParams {
        use_bond_types: false,
        include_ring_membership: false,
        include_redundant_environments: true,
        ..Default::default()
    };
    let updated=original.with_json(r#"{"radius":"4294967295","onlyNonzeroInvariants":"1","useBondTypes":true,"includeRedundantEnvironments":false,"includeRingMembership":true}"#).unwrap();
    assert_eq!(updated.radius, u32::MAX);
    assert!(updated.only_nonzero_invariants);
    assert!(!updated.use_bond_types);
    assert!(!updated.include_ring_membership);
    assert!(updated.include_redundant_environments);
    assert_eq!(original.radius, 3);
    assert_eq!(original.count_bounds, [1, 2, 4, 8]);
    assert!(updated.count_bounds.is_empty());
    let v: Value = serde_json::from_str(&updated.to_json()).unwrap();
    assert!(v.get("useBondTypes").is_none());
    assert!(v.get("includeRedundantEnvironments").is_none());
    assert!(v.get("includeRingMembership").is_none());
    assert_eq!(v["countBounds"], "");
}
#[test]
fn common_update_and_absent_empty_bounds_have_source_semantics() {
    let full = r#"{"includeChirality":"true","countSimulation":false,"fpSize":"1000","numBitsPerFeature":3,"countBounds":[1,"3",7]}"#;
    let m = MorganParams::default().with_json(full).unwrap();
    let a = AtomPairParams::default().with_json(full).unwrap();
    let t = TopologicalTorsionParams::default().with_json(full).unwrap();
    for (chi, sim, size, bits, bounds) in [
        (
            m.include_chirality,
            m.count_simulation,
            m.fp_size,
            m.bits_per_feature,
            m.count_bounds,
        ),
        (
            a.include_chirality,
            a.count_simulation,
            a.fp_size,
            a.bits_per_feature,
            a.count_bounds,
        ),
        (
            t.include_chirality,
            t.count_simulation,
            t.fp_size,
            t.bits_per_feature,
            t.count_bounds,
        ),
    ] {
        assert!(chi);
        assert!(!sim);
        assert_eq!(size, 1000);
        assert_eq!(bits, 3);
        assert_eq!(bounds, [1, 3, 7]);
    }
    for text in ["{}", r#"{"countBounds":[]}"#, r#"{"countBounds":""}"#] {
        assert!(
            MorganParams::default()
                .with_json(text)
                .unwrap()
                .count_bounds
                .is_empty()
        );
        assert!(
            AtomPairParams::default()
                .with_json(text)
                .unwrap()
                .count_bounds
                .is_empty()
        );
        assert!(
            TopologicalTorsionParams::default()
                .with_json(text)
                .unwrap()
                .count_bounds
                .is_empty()
        );
    }
    for text in ["", "  \n\t"] {
        assert_eq!(
            MorganParams::default().with_json(text).unwrap(),
            MorganParams::default()
        );
        assert_eq!(
            AtomPairParams::default().with_json(text).unwrap(),
            AtomPairParams::default()
        );
        assert_eq!(
            TopologicalTorsionParams::default().with_json(text).unwrap(),
            TopologicalTorsionParams::default()
        );
    }
}
#[test]
fn derived_updates_do_not_repeat_factory_preconditions() {
    let a=AtomPairParams::default().with_json(r#"{"use2D":0,"minDistance":30,"maxDistance":1,"numBitsPerFeature":0,"countBounds":[]}"#).unwrap();
    assert!(!a.use_2d);
    assert_eq!(a.min_distance, 30);
    assert_eq!(a.max_distance, 1);
    assert_eq!(a.bits_per_feature, 0);
    assert!(a.count_bounds.is_empty());
    let t = TopologicalTorsionParams::default()
        .with_json(r#"{"torsionAtomCount":"8","onlyShortestPaths":"true","numBitsPerFeature":0}"#)
        .unwrap();
    assert_eq!(t.torsion_atom_count, 8);
    assert!(t.only_shortest_paths);
    assert_eq!(t.bits_per_feature, 0);
    let m = MorganParams::default()
        .with_json(r#"{"countSimulation":true,"countBounds":[],"numBitsPerFeature":0}"#)
        .unwrap();
    assert!(m.count_simulation);
    assert!(m.count_bounds.is_empty());
    assert_eq!(m.bits_per_feature, 0);
}
#[test]
fn failure_keeps_original_input_and_has_structured_source_error() {
    let m = MorganParams::default();
    let a = AtomPairParams::default();
    let t = TopologicalTorsionParams::default();
    // Original v9 error inputs retained, corrected by actual pinned Boost helper.
    for text in [
        "[1]",
        r#"{"fpSize":-1}"#,
        r#"{"includeChirality":"invalid"}"#,
        r#"{"countBounds":[-1]}"#,
    ] {
        let values = [
            m.with_json(text).unwrap().to_json(),
            a.with_json(text).unwrap().to_json(),
            t.with_json(text).unwrap().to_json(),
        ];
        for encoded in values {
            let value: Value = serde_json::from_str(&encoded).unwrap();
            assert_eq!(value["includeChirality"], "false");
            assert_eq!(
                value["fpSize"],
                if text.contains("fpSize") {
                    "4294967295"
                } else {
                    "2048"
                }
            );
            assert_eq!(
                value["countBounds"],
                if text.contains("countBounds") {
                    json!(["4294967295"])
                } else {
                    json!("")
                }
            );
        }
    }
    // Strict get_value<T>() on bounds children still fails, unlike get<T>(default).
    for text in [
        r#"{"countBounds":[1,"invalid"]}"#,
        r#"{"countBounds":[{}]}"#,
        r#"{"countBounds":[4294967296]}"#,
    ] {
        assert!(matches!(
            m.with_json(text),
            Err(FingerprintJsonError::Invalid(_))
        ));
        assert!(matches!(
            a.with_json(text),
            Err(FingerprintJsonError::Invalid(_))
        ));
        assert!(matches!(
            t.with_json(text),
            Err(FingerprintJsonError::Invalid(_))
        ));
    }
    for error in [
        m.with_json("{").unwrap_err(),
        a.with_json("{").unwrap_err(),
        t.with_json("{").unwrap_err(),
    ] {
        assert!(matches!(error, FingerprintJsonError::Parse(_)));
        assert!(std::error::Error::source(&error).is_some());
    }
    let updated = m
        .with_json(r#"{"radius":99,"countBounds":[1,-1]}"#)
        .unwrap();
    assert_eq!(updated.radius, 99);
    assert_eq!(updated.count_bounds, [1, u32::MAX]);
    let updated = a.with_json(r#"{"use2D":false,"minDistance":-1}"#).unwrap();
    assert!(!updated.use_2d);
    assert_eq!(updated.min_distance, u32::MAX);
    let updated = t
        .with_json(r#"{"torsionAtomCount":7,"onlyShortestPaths":"invalid"}"#)
        .unwrap();
    assert_eq!(updated.torsion_atom_count, 7);
    assert!(!updated.only_shortest_paths);
    assert_eq!(m, MorganParams::default());
    assert_eq!(a, AtomPairParams::default());
    assert_eq!(t, TopologicalTorsionParams::default());
}

#[test]
fn independent_p1_seven_source_conditions_preserve_first_keys_and_arbitrary_children() {
    for (text, size, bounds) in [
        (r#"{"includeChirality":"invalid"}"#, 2048, vec![]),
        (r#"{"fpSize":-1}"#, u32::MAX, vec![]),
        (r#"{"fpSize":"invalid"}"#, 2048, vec![]),
        (r#"{"fpSize":4294967296}"#, 2048, vec![]),
        (r#"{"countBounds":{"x":1,"y":3}}"#, 2048, vec![1, 3]),
        (r#"{"countBounds":"bad"}"#, 2048, vec![]),
        (r#"{"fpSize":1,"fpSize":2}"#, 1, vec![]),
    ] {
        let values = [
            MorganParams::default().with_json(text).unwrap().to_json(),
            AtomPairParams::default().with_json(text).unwrap().to_json(),
            TopologicalTorsionParams::default()
                .with_json(text)
                .unwrap()
                .to_json(),
        ];
        for encoded in values {
            let value: Value = serde_json::from_str(&encoded).unwrap();
            assert_eq!(value["fpSize"], size.to_string());
            assert_eq!(value["includeChirality"], "false");
            let expected = if bounds.is_empty() {
                json!("")
            } else {
                json!(bounds.iter().map(u32::to_string).collect::<Vec<_>>())
            };
            assert_eq!(value["countBounds"], expected);
        }
    }
    let m = MorganParams {
        fp_size: 19,
        include_chirality: true,
        ..Default::default()
    }
    .with_json(r#"{"fpSize":"bad","includeChirality":"bad","countBounds":{"z":3,"a":1,"z":7}}"#)
    .unwrap();
    assert_eq!(m.fp_size, 19);
    assert!(m.include_chirality);
    assert_eq!(m.count_bounds, [3, 1, 7]);
}
#[test]
fn source_lexical_numbers_and_boolean_stream_retry_are_not_json_normalized() {
    for text in [
        r#"{"fpSize":1e0,"includeChirality":1e0}"#,
        r#"{"fpSize":1.0,"includeChirality":1.0}"#,
    ] {
        let m = MorganParams::default().with_json(text).unwrap();
        assert_eq!(m.fp_size, 2048);
        assert!(!m.include_chirality);
    }
    for (text, flag) in [
        ("2true", true),
        ("2 true", true),
        ("-1true", true),
        ("+true", true),
        ("+false", false),
        ("1true", false),
        ("+1", true),
        ("01", true),
        ("-0", false),
    ] {
        let encoded = json!({"includeChirality":text}).to_string();
        let m = MorganParams::default().with_json(&encoded).unwrap();
        assert_eq!(m.include_chirality, flag, "{text}");
    }
    for text in ["[1]", "true", "null", "1"] {
        assert!(
            MorganParams::default()
                .with_json(text)
                .unwrap()
                .count_bounds
                .is_empty()
        );
    }
    let m = MorganParams::default()
        .with_json(r#"{"fpSize":"-4294967295","countBounds":["+1","-4294967295",-1]}"#)
        .unwrap();
    assert_eq!(m.fp_size, 1);
    assert_eq!(m.count_bounds, [1, 1, u32::MAX]);
}
