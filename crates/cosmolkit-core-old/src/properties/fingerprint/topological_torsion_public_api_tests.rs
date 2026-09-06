use super::{
    FingerprintError, TopologicalFingerprintParams, TopologicalTorsionFingerprintOutputRequest,
    TopologicalTorsionFingerprintParams, TopologicalTorsionFingerprintValue,
    TopologicalTorsionFingerprintVector, TopologicalTorsionLegacyKind,
    TopologicalTorsionLegacyParams, TopologicalTorsionLegacyResult, topological_fingerprint,
    topological_torsion_count_fingerprint, topological_torsion_fingerprint,
    topological_torsion_fingerprint_with_output, topological_torsion_generator,
    topological_torsion_legacy_fingerprint, topological_torsion_sparse_count_fingerprint,
    topological_torsion_sparse_fingerprint,
};
use crate::Molecule;

#[test]
fn public_defaults_match_generator_defaults_and_remain_named_separately_from_rdkfingerprint() {
    let params = TopologicalTorsionFingerprintParams::default();
    assert!(!params.include_chirality);
    assert_eq!(params.torsion_atom_count, 4);
    assert!(params.count_simulation);
    assert_eq!(params.count_bounds, vec![1, 2, 4, 8]);
    assert_eq!(params.fp_size, 2048);
    assert_eq!(params.num_bits_per_feature, 1);
    assert!(!params.only_shortest_paths);
    assert!(params.from_atoms.is_none());
    assert!(params.ignore_atoms.is_none());
    assert!(params.custom_atom_invariants.is_none());

    fn accepts_rdkfingerprint(_: &TopologicalFingerprintParams) {}
    fn accepts_torsion(_: &TopologicalTorsionFingerprintParams) {}
    accepts_rdkfingerprint(&TopologicalFingerprintParams::default());
    accepts_torsion(&params);

    let molecule = Molecule::from_smiles("CCCCO").expect("molecule");
    let rdk = topological_fingerprint(&molecule, &TopologicalFingerprintParams::default())
        .expect("RDKFingerprint");
    let torsion = topological_torsion_fingerprint(&molecule, &params).expect("torsion");
    assert_ne!(rdk, torsion);
}

#[test]
fn public_scalar_conveniences_expose_all_four_shared_vector_forms() {
    let molecule = Molecule::from_smiles("CCCCO").expect("molecule");
    let params = TopologicalTorsionFingerprintParams::default();

    let sparse_count =
        topological_torsion_sparse_count_fingerprint(&molecule, &params).expect("sparse count");
    let sparse_bit =
        topological_torsion_sparse_fingerprint(&molecule, &params).expect("sparse bit");
    let count = topological_torsion_count_fingerprint(&molecule, &params).expect("folded count");
    let bit = topological_torsion_fingerprint(&molecule, &params).expect("explicit bit");

    assert_eq!(sparse_count.nonzero_elements().values().sum::<i32>(), 2);
    assert_eq!(sparse_bit.on_bits().len(), 2);
    assert_eq!(count.size(), 2048);
    assert_eq!(count.nonzero_elements().values().sum::<i32>(), 2);
    assert_eq!(bit.n_bits(), 2048);
    assert_eq!(bit.on_bits().len(), 2);

    for vector in [
        TopologicalTorsionFingerprintVector::SparseCount,
        TopologicalTorsionFingerprintVector::SparseBit,
        TopologicalTorsionFingerprintVector::Count,
        TopologicalTorsionFingerprintVector::Bit,
    ] {
        let result = topological_torsion_fingerprint_with_output(
            &molecule,
            &params,
            TopologicalTorsionFingerprintOutputRequest {
                vector,
                ..Default::default()
            },
        )
        .expect("requested vector");
        assert!(result.additional_output.is_none());
        assert!(matches!(
            (vector, result.fingerprint),
            (
                TopologicalTorsionFingerprintVector::SparseCount,
                TopologicalTorsionFingerprintValue::SparseCount(_)
            ) | (
                TopologicalTorsionFingerprintVector::SparseBit,
                TopologicalTorsionFingerprintValue::SparseBit(_)
            ) | (
                TopologicalTorsionFingerprintVector::Count,
                TopologicalTorsionFingerprintValue::Count(_)
            ) | (
                TopologicalTorsionFingerprintVector::Bit,
                TopologicalTorsionFingerprintValue::Bit(_)
            )
        ));
    }
}

#[test]
fn public_output_request_returns_every_supported_provenance_allocation() {
    let molecule = Molecule::from_smiles("CCCCO").expect("molecule");
    let result = topological_torsion_fingerprint_with_output(
        &molecule,
        &TopologicalTorsionFingerprintParams::default(),
        TopologicalTorsionFingerprintOutputRequest {
            vector: TopologicalTorsionFingerprintVector::Bit,
            atom_to_bits: true,
            atom_counts: true,
            bit_paths: true,
            atoms_per_bit: true,
        },
    )
    .expect("result");
    let output = result.additional_output.expect("provenance");
    assert_eq!(
        output.atom_to_bits.as_ref().unwrap().len(),
        molecule.num_atoms()
    );
    assert_eq!(
        output.atom_counts.as_ref().unwrap().len(),
        molecule.num_atoms()
    );
    assert!(!output.bit_paths.as_ref().unwrap().is_empty());
    assert!(!output.atoms_per_bit.as_ref().unwrap().is_empty());
    assert!(output.bit_info_map.is_none());
}

#[test]
fn public_api_returns_structured_errors_for_invalid_parameters() {
    let molecule = Molecule::from_smiles("CCCCO").expect("molecule");
    let mut params = TopologicalTorsionFingerprintParams {
        num_bits_per_feature: 0,
        ..Default::default()
    };
    assert!(matches!(
        topological_torsion_generator(&params),
        Err(FingerprintError::InvalidArguments { .. })
    ));

    params.num_bits_per_feature = 1;
    params.fp_size = 0;
    assert!(matches!(
        topological_torsion_fingerprint(&molecule, &params),
        Err(FingerprintError::InvalidArguments { .. })
    ));

    params.fp_size = 2048;
    params.custom_atom_invariants = Some(vec![1]);
    assert!(matches!(
        topological_torsion_count_fingerprint(&molecule, &params),
        Err(FingerprintError::InvalidArguments { .. })
    ));

    params.custom_atom_invariants = None;
    params.torsion_atom_count = 8;
    assert!(matches!(
        topological_torsion_generator(&params),
        Err(FingerprintError::InvalidArguments { .. })
    ));
}

#[test]
fn public_calls_preserve_source_molecule_and_are_repeatably_deterministic() {
    let molecule = Molecule::from_smiles("CC[C@H](F)Cl")
        .expect("molecule")
        .with_prop("source", "unchanged");
    let before = molecule.clone();
    let params = TopologicalTorsionFingerprintParams {
        include_chirality: true,
        num_bits_per_feature: 3,
        ..Default::default()
    };

    let first = topological_torsion_fingerprint_with_output(
        &molecule,
        &params,
        TopologicalTorsionFingerprintOutputRequest {
            vector: TopologicalTorsionFingerprintVector::Count,
            bit_paths: true,
            ..Default::default()
        },
    )
    .expect("first");
    let second = topological_torsion_fingerprint_with_output(
        &molecule,
        &params,
        TopologicalTorsionFingerprintOutputRequest {
            vector: TopologicalTorsionFingerprintVector::Count,
            bit_paths: true,
            ..Default::default()
        },
    )
    .expect("second");

    assert_eq!(first, second);
    assert_eq!(molecule, before);
    assert_eq!(molecule.prop("source"), Some("unchanged"));
}

#[test]
fn public_legacy_compatibility_type_routes_all_three_legacy_forms() {
    let molecule = Molecule::from_smiles("CCCCO").expect("molecule");
    let unfolded = topological_torsion_legacy_fingerprint(
        &molecule,
        &TopologicalTorsionLegacyParams::default(),
    )
    .expect("unfolded");
    assert!(matches!(
        unfolded,
        TopologicalTorsionLegacyResult::SparseCount(_)
    ));

    let hashed = topological_torsion_legacy_fingerprint(
        &molecule,
        &TopologicalTorsionLegacyParams {
            kind: TopologicalTorsionLegacyKind::HashedCount,
            n_bits: 1000,
            ..Default::default()
        },
    )
    .expect("hashed");
    let TopologicalTorsionLegacyResult::SparseCount(hashed) = hashed else {
        panic!("expected hashed count")
    };
    assert_eq!(
        hashed
            .nonzero_elements()
            .keys()
            .copied()
            .collect::<Vec<_>>(),
        vec![24, 288]
    );

    let bits = topological_torsion_legacy_fingerprint(
        &molecule,
        &TopologicalTorsionLegacyParams {
            kind: TopologicalTorsionLegacyKind::HashedBit,
            ..Default::default()
        },
    )
    .expect("bits");
    assert!(matches!(bits, TopologicalTorsionLegacyResult::Bit(_)));
}
