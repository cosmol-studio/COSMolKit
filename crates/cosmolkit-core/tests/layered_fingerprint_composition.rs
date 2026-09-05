use std::sync::Arc;

use cosmolkit_core::{
    AtomPairFingerprintParams, AvalonFingerprintParams, Fingerprint, LayeredFingerprintLayers,
    LayeredFingerprintParams, MaccsFingerprintParams, Molecule, MorganFingerprintParams, TopologicalFingerprintParams,
    TopologicalTorsionFingerprintParams, topological_torsion_fingerprint,
};

fn fingerprint_key(fingerprint: Fingerprint) -> (usize, Vec<usize>) {
    (fingerprint.n_bits(), fingerprint.on_bits())
}

fn every_family_key(molecule: &Molecule) -> Vec<(usize, Vec<usize>)> {
    vec![
        fingerprint_key(
            molecule
                .morgan_fingerprint(&MorganFingerprintParams::default())
                .expect("Morgan fingerprint"),
        ),
        fingerprint_key(
            molecule
                .maccs_fingerprint(&MaccsFingerprintParams::default())
                .expect("MACCS fingerprint"),
        ),
        fingerprint_key(
            molecule
                .avalon_fingerprint(&AvalonFingerprintParams::default())
                .expect("Avalon fingerprint"),
        ),
        fingerprint_key(
            molecule
                .topological_fingerprint(&TopologicalFingerprintParams::default())
                .expect("RDKFingerprint"),
        ),
        fingerprint_key(
            molecule
                .atom_pair_fingerprint(&AtomPairFingerprintParams::default())
                .expect("AtomPair fingerprint"),
        ),
        fingerprint_key(
            topological_torsion_fingerprint(molecule, &TopologicalTorsionFingerprintParams::default())
                .expect("Topological Torsion fingerprint"),
        ),
        fingerprint_key(
            molecule
                .layered_fingerprint(&LayeredFingerprintParams::default())
                .expect("Layered fingerprint"),
        ),
    ]
}

fn configured_layered_key(molecule: &Molecule, params: &LayeredFingerprintParams) -> (usize, Vec<usize>) {
    fingerprint_key(
        molecule
            .layered_fingerprint(params)
            .expect("configured Layered fingerprint"),
    )
}

#[test]
fn layered_calls_compose_with_every_fingerprint_family_without_shared_state_crosstalk() {
    let molecule = Arc::new(
        Molecule::from_smiles("CC[C@H](F)Cl.c1ccncc1O")
            .expect("shared molecule")
            .with_prop("composition", "unchanged"),
    );
    let before = molecule.as_ref().clone();
    let before_smiles = molecule.to_smiles(true).expect("source SMILES");

    let topology_only = LayeredFingerprintParams {
        layers: LayeredFingerprintLayers::TOPOLOGY,
        min_path: 1,
        max_path: 4,
        fp_size: 257,
        atom_counts: Some(vec![7; molecule.num_atoms()]),
        set_only_bits: None,
        branched_paths: false,
        from_atoms: Some(vec![0, 6]),
    };
    let ring_and_aromatic = LayeredFingerprintParams {
        layers: LayeredFingerprintLayers::RING_PRESENCE
            | LayeredFingerprintLayers::RING_SIZE
            | LayeredFingerprintLayers::AROMATICITY,
        min_path: 2,
        max_path: 6,
        fp_size: 509,
        atom_counts: None,
        set_only_bits: None,
        branched_paths: true,
        from_atoms: None,
    };

    let family_baseline = every_family_key(&molecule);
    let topology_baseline = configured_layered_key(&molecule, &topology_only);
    let ring_baseline = configured_layered_key(&molecule, &ring_and_aromatic);

    for _ in 0..4 {
        assert_eq!(configured_layered_key(&molecule, &ring_and_aromatic), ring_baseline);
        assert_eq!(every_family_key(&molecule), family_baseline);
        assert_eq!(configured_layered_key(&molecule, &topology_only), topology_baseline);
        assert_eq!(every_family_key(&molecule), family_baseline);
    }

    let expected = (
        family_baseline.clone(),
        topology_baseline.clone(),
        ring_baseline.clone(),
    );
    let handles = (0..12)
        .map(|_| {
            let molecule = Arc::clone(&molecule);
            let topology_only = topology_only.clone();
            let ring_and_aromatic = ring_and_aromatic.clone();
            std::thread::spawn(move || {
                (
                    every_family_key(&molecule),
                    configured_layered_key(&molecule, &topology_only),
                    configured_layered_key(&molecule, &ring_and_aromatic),
                )
            })
        })
        .collect::<Vec<_>>();
    for handle in handles {
        assert_eq!(handle.join().expect("fingerprint worker"), expected);
    }

    assert_eq!(molecule.as_ref(), &before);
    assert_eq!(molecule.prop("composition"), Some("unchanged"));
    assert_eq!(
        molecule.to_smiles(true).expect("unchanged source SMILES"),
        before_smiles
    );
}
