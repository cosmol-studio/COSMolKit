use super::{FingerprintError, TopologicalTorsionArguments, TopologicalTorsionEnvGenerator};

#[test]
fn defaults_and_constructor_match_the_pinned_source_contract() {
    let arguments = TopologicalTorsionArguments::default();
    let common = &arguments.fingerprint_arguments;

    assert_eq!(arguments.d_torsion_atom_count, 4);
    assert!(!arguments.df_only_shortest_paths);
    assert!(common.df_count_simulation);
    assert!(!common.df_include_chirality);
    assert_eq!(common.d_count_bounds, vec![1, 2, 4, 8]);
    assert_eq!(common.d_fp_size, 2048);
    assert_eq!(common.d_num_bits_per_feature, 1);

    let constructed = TopologicalTorsionArguments::new(true, 5, false, vec![2, 3, 7], 4096)
        .expect("valid custom arguments");
    assert_eq!(constructed.d_torsion_atom_count, 5);
    assert!(!constructed.df_only_shortest_paths);
    assert!(constructed.fingerprint_arguments.df_include_chirality);
    assert!(!constructed.fingerprint_arguments.df_count_simulation);
    assert_eq!(
        constructed.fingerprint_arguments.d_count_bounds,
        vec![2, 3, 7]
    );
    assert_eq!(constructed.fingerprint_arguments.d_fp_size, 4096);
    assert_eq!(constructed.fingerprint_arguments.d_num_bits_per_feature, 1);

    assert!(matches!(
        TopologicalTorsionArguments::new(false, 4, true, Vec::new(), 2048),
        Err(FingerprintError::InvalidArguments { .. })
    ));
}

#[test]
fn mutable_options_and_information_strings_use_the_shared_common_arguments() {
    let mut arguments = TopologicalTorsionArguments::default();
    arguments.d_torsion_atom_count = 3;
    arguments.df_only_shortest_paths = true;
    arguments.fingerprint_arguments.df_count_simulation = false;
    arguments.fingerprint_arguments.df_include_chirality = true;
    arguments.fingerprint_arguments.d_count_bounds = vec![3, 9];
    arguments.fingerprint_arguments.d_fp_size = 1024;
    arguments.fingerprint_arguments.d_num_bits_per_feature = 2;

    assert_eq!(
        arguments.info_string(),
        "TopologicalTorsionArguments torsionAtomCount=3 onlyShortestPaths=1"
    );
    assert_eq!(arguments.infoString(), arguments.info_string());
    assert_eq!(
        arguments.fingerprint_arguments.common_arguments_string(),
        "Common arguments : countSimulation=0 fpSize=1024 bitsPerFeature=2 includeChirality=1"
    );
}

#[test]
fn json_roundtrip_updates_derived_and_common_fields() {
    let mut arguments = TopologicalTorsionArguments::default();
    arguments.d_torsion_atom_count = 5;
    arguments.df_only_shortest_paths = true;
    arguments.fingerprint_arguments.df_count_simulation = false;
    arguments.fingerprint_arguments.df_include_chirality = true;
    arguments.fingerprint_arguments.d_count_bounds = vec![2, 6, 10];
    arguments.fingerprint_arguments.d_fp_size = 4096;
    arguments.fingerprint_arguments.d_num_bits_per_feature = 3;

    let json = arguments.to_json();
    let value: serde_json::Value = serde_json::from_str(&json).expect("valid JSON");
    assert_eq!(value["type"], "TopologicalTorsionArguments");
    assert_eq!(value["torsionAtomCount"], "5");
    assert_eq!(value["onlyShortestPaths"], "true");
    assert_eq!(value["countSimulation"], "false");
    assert_eq!(value["includeChirality"], "true");
    assert_eq!(value["countBounds"], serde_json::json!(["2", "6", "10"]));
    assert_eq!(value["fpSize"], "4096");
    assert_eq!(value["numBitsPerFeature"], "3");
    assert_eq!(arguments.toJSON(), json);

    let mut restored = TopologicalTorsionArguments::default();
    restored.from_json(&json).expect("roundtrip");
    assert_eq!(restored, arguments);
}

#[test]
fn partial_json_preserves_scalar_defaults_and_clears_missing_count_bounds_like_source() {
    let mut arguments = TopologicalTorsionArguments::default();
    arguments
        .fromJSON(r#"{"torsionAtomCount":5,"onlyShortestPaths":true,"fpSize":1024}"#)
        .expect("partial update");

    assert_eq!(arguments.d_torsion_atom_count, 5);
    assert!(arguments.df_only_shortest_paths);
    assert_eq!(arguments.fingerprint_arguments.d_fp_size, 1024);
    assert!(arguments.fingerprint_arguments.df_count_simulation);
    assert!(!arguments.fingerprint_arguments.df_include_chirality);
    assert_eq!(arguments.fingerprint_arguments.d_num_bits_per_feature, 1);
    assert!(arguments.fingerprint_arguments.d_count_bounds.is_empty());

    arguments
        .from_json(r#"{"countBounds":[1,"4",8],"includeChirality":"true"}"#)
        .expect("property-tree-compatible string values");
    assert_eq!(
        arguments.fingerprint_arguments.d_count_bounds,
        vec![1, 4, 8]
    );
    assert!(arguments.fingerprint_arguments.df_include_chirality);
}

#[test]
fn malformed_json_and_invalid_field_types_return_structured_errors() {
    let mut arguments = TopologicalTorsionArguments::default();
    assert!(matches!(
        arguments.from_json("{"),
        Err(FingerprintError::InvalidArgumentsJson(_))
    ));
    assert!(matches!(
        arguments.from_json("[]"),
        Err(FingerprintError::InvalidArgumentsJson(_))
    ));
    assert!(matches!(
        arguments.from_json(r#"{"torsionAtomCount":-1}"#),
        Err(FingerprintError::InvalidArgumentsJson(_))
    ));
    assert!(matches!(
        arguments.from_json(r#"{"onlyShortestPaths":"unknown"}"#),
        Err(FingerprintError::InvalidArgumentsJson(_))
    ));
    assert!(matches!(
        arguments.from_json(r#"{"countBounds":{}}"#),
        Err(FingerprintError::InvalidArgumentsJson(_))
    ));
}

#[test]
fn argument_unit_ignores_unknown_type_tags_and_fields_like_the_source_method() {
    let mut arguments = TopologicalTorsionArguments::default();
    arguments
        .from_json(
            r#"{"type":"UnknownArguments","unknownOption":17,"torsionAtomCount":2,"countBounds":[1]}"#,
        )
        .expect("the argument unit does not dispatch on type");

    assert_eq!(arguments.d_torsion_atom_count, 2);
    assert_eq!(arguments.fingerprint_arguments.d_count_bounds, vec![1]);
}

#[test]
fn result_sizes_cover_zero_default_chiral_and_maximum_defined_shifts() {
    let generator = TopologicalTorsionEnvGenerator::new();
    let mut arguments = TopologicalTorsionArguments::default();

    assert_eq!(generator.get_result_size(&arguments).unwrap(), 1_u64 << 36);
    assert_eq!(generator.getResultSize(&arguments).unwrap(), 1_u64 << 36);

    arguments.d_torsion_atom_count = 0;
    assert_eq!(generator.get_result_size(&arguments).unwrap(), 1);

    arguments.d_torsion_atom_count = 7;
    assert_eq!(generator.get_result_size(&arguments).unwrap(), 1_u64 << 63);

    arguments.fingerprint_arguments.df_include_chirality = true;
    arguments.d_torsion_atom_count = 5;
    assert_eq!(generator.get_result_size(&arguments).unwrap(), 1_u64 << 55);
}

#[test]
fn result_size_rejects_every_source_undefined_width_mapping() {
    let generator = TopologicalTorsionEnvGenerator::new();
    let mut arguments = TopologicalTorsionArguments::default();

    arguments.d_torsion_atom_count = 8;
    assert!(matches!(
        generator.get_result_size(&arguments),
        Err(FingerprintError::InvalidArguments { .. })
    ));

    arguments.fingerprint_arguments.df_include_chirality = true;
    arguments.d_torsion_atom_count = 6;
    assert!(matches!(
        generator.get_result_size(&arguments),
        Err(FingerprintError::InvalidArguments { .. })
    ));

    arguments.d_torsion_atom_count = u32::MAX;
    assert!(matches!(
        generator.get_result_size(&arguments),
        Err(FingerprintError::InvalidArguments { .. })
    ));
}
