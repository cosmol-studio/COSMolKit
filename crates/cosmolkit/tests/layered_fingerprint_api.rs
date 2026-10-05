//! Original source-backed Layered observations through the canonical public API.
#![cfg(feature = "cap-fingerprints")]
use cosmolkit::{
    BINDING_CONTRACT, Fingerprint, LayeredFingerprintError, LayeredFingerprintLayers as L,
    LayeredFingerprintParams as P, Molecule, layered_query_fingerprint_with_params,
};
#[test]
fn original_default_mask_seed_and_root_observations_use_public_facade() {
    let mol = Molecule::from_smiles("CCO").unwrap();
    let before = mol.to_smiles().unwrap();
    assert_eq!(
        mol.layered_fingerprint().unwrap().on_bits(),
        [92, 360, 596, 610, 611, 674, 867, 1044, 1111, 1783, 1784]
    );
    assert!(
        mol.layered_fingerprint_with_output()
            .unwrap()
            .atom_counts()
            .is_none()
    );
    let p = P {
        atom_counts: Some(vec![10, 20, 30]),
        set_only_bits: Some(Fingerprint::from_on_bits(2048, [674]).unwrap()),
        ..Default::default()
    };
    let result = mol.layered_fingerprint_with_output_with_params(&p).unwrap();
    assert_eq!(result.fingerprint().on_bits(), [674]);
    assert_eq!(result.atom_counts(), Some([11, 22, 31].as_slice()));
    assert_eq!(p.atom_counts, Some(vec![10, 20, 30]));
    assert_eq!(
        mol.layered_fingerprint_with_params(&P {
            layers: L::TOPOLOGY,
            ..Default::default()
        })
        .unwrap()
        .on_bits(),
        [674, 867]
    );
    assert_eq!(mol.to_smiles().unwrap(), before);
}
#[test]
fn source_distinguishes_empty_roots_high_flags_and_count_absence() {
    let mol = Molecule::from_smiles("CCO").unwrap();
    for p in [
        P {
            from_atoms: Some(vec![]),
            ..Default::default()
        },
        P {
            layers: L::from_bits_retain(0xffff_ffc0),
            ..Default::default()
        },
    ] {
        assert!(
            mol.layered_fingerprint_with_params(&p)
                .unwrap()
                .on_bits()
                .is_empty()
        );
    }
    let p = P {
        branched_paths: false,
        from_atoms: Some(vec![0]),
        ..Default::default()
    };
    assert_eq!(
        mol.layered_fingerprint_with_params(&p).unwrap().on_bits(),
        [360, 596, 610, 611, 674, 867, 1044, 1111, 1783, 1784]
    );
    assert!(matches!(
        mol.layered_fingerprint_with_params(&P {
            from_atoms: Some(vec![3]),
            ..Default::default()
        }),
        Err(LayeredFingerprintError::InvalidArguments {
            reason: "fromAtoms contains atom index out of range"
        })
    ));
}
#[test]
fn canonical_query_graph_retains_non_element_identity_and_original_masks() {
    for (text, bits) in [
        ("C-C", vec![20, 98, 99, 162, 360]),
        ("C~C", vec![20, 98, 99, 162]),
        ("[C,N]-C", vec![162, 360]),
        ("[C,N]~C", vec![162]),
        ("C-[C,N]", vec![162, 360]),
        ("C~[C,N]", vec![162]),
        ("[C,N]-[C,N]", vec![162, 360]),
        ("[C,N]~[C,N]", vec![162]),
    ] {
        let graph = cosmolkit::search::parse_smarts(text).unwrap();
        let before = graph.clone();
        let result = layered_query_fingerprint_with_params(
            &graph,
            &P {
                fp_size: 512,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(result.on_bits(), bits, "{text}");
        assert_eq!(graph, before);
    }
}
#[test]
fn declared_callable_contracts_resolve_to_real_layered_implementations() {
    for id in [
        "Molecule.layered_fingerprint",
        "Molecule.layered_fingerprint_with_params",
        "Molecule.layered_fingerprint_with_output",
        "Molecule.layered_fingerprint_with_output_with_params",
        "LayeredFingerprintResult.fingerprint",
        "LayeredFingerprintResult.atom_counts",
        "layered_query_fingerprint_with_params",
        "layered_query_fingerprint_with_output_with_params",
    ] {
        let row = BINDING_CONTRACT
            .iter()
            .find(|r| r.semantic_id == id)
            .unwrap();
        assert!(row.callable.is_some());
        assert_eq!(row.feature, "cap-fingerprints");
    }
}
