use std::collections::BTreeMap;

use cosmolkit_core::{
    AtomPairFingerprintGenerator, AtomPairFingerprintParams, FingerprintFuncArguments, Molecule,
};

fn map(entries: &[(u64, i32)]) -> BTreeMap<u64, i32> {
    entries.iter().copied().collect()
}

#[test]
fn modern_core_matches_deprecated_rdkit_default_and_count_simulation_outputs() {
    // Oracle: pinned RDKit 2026.03.1 rdMolDescriptors AtomPair entry points.
    let molecule = Molecule::from_smiles("CCCO").unwrap();
    let default_generator =
        AtomPairFingerprintGenerator::new(&AtomPairFingerprintParams::default()).unwrap();
    let sparse = default_generator
        .sparse_count_fingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .unwrap();
    assert_eq!(
        sparse.nonzero_elements(),
        &map(&[
            (558_113, 1),
            (558_114, 1),
            (558_145, 1),
            (1_590_307, 1),
            (1_590_337, 1),
            (1_590_338, 1),
        ])
    );

    let generator_128 = AtomPairFingerprintGenerator::new(&AtomPairFingerprintParams {
        n_bits: 128,
        ..Default::default()
    })
    .unwrap();
    let count = generator_128
        .count_fingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .unwrap();
    assert_eq!(
        count.nonzero_elements(),
        &map(&[(15, 1), (28, 1), (30, 1), (78, 2), (95, 1)])
    );
    let four_bits = generator_128
        .fingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .unwrap();
    assert_eq!(four_bits.on_bits(), &[56, 57, 60, 112, 120, 124]);

    let generator_96 = AtomPairFingerprintGenerator::new(&AtomPairFingerprintParams {
        n_bits: 96,
        count_bounds: vec![1, 2, 3],
        ..Default::default()
    })
    .unwrap();
    let three_bits = generator_96
        .fingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .unwrap();
    assert_eq!(three_bits.on_bits(), &[42, 43, 45, 84, 90, 93]);
}

#[test]
fn modern_core_matches_deprecated_rdkit_custom_chiral_and_filter_outputs() {
    // Oracle: pinned RDKit 2026.03.1 with the deprecated custom-invariant path.
    let molecule = Molecule::from_smiles("C[C@H](O)F").unwrap();
    let generator = AtomPairFingerprintGenerator::new(&AtomPairFingerprintParams {
        n_bits: 64,
        min_distance: 1,
        max_distance: 2,
        use_chirality: true,
        ..Default::default()
    })
    .unwrap();
    let arguments = || FingerprintFuncArguments {
        from_atoms: Some(vec![1]),
        ignore_atoms: Some(vec![3]),
        custom_atom_invariants: Some(vec![11, 22, 33, 44]),
        ..Default::default()
    };

    let sparse = generator
        .sparse_count_fingerprint(&molecule, &mut arguments())
        .unwrap();
    assert_eq!(
        sparse.nonzero_elements(),
        &map(&[(1_442_145, 1), (2_163_393, 1)])
    );
    let count = generator
        .count_fingerprint(&molecule, &mut arguments())
        .unwrap();
    assert_eq!(count.nonzero_elements(), &map(&[(0, 1), (37, 1)]));
    let explicit = generator.fingerprint(&molecule, &mut arguments()).unwrap();
    assert_eq!(explicit.on_bits(), &[0, 20]);
}

#[test]
fn modern_core_matches_deprecated_rdkit_three_dimensional_outputs() {
    // Oracle: pinned RDKit 2026.03.1 with use2D=false and confId=0.
    let molecule = Molecule::from_mol_block(
        "\n     RDKit          3D\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.6000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\nM  END\n",
    )
    .unwrap();
    let generator = AtomPairFingerprintGenerator::new(&AtomPairFingerprintParams {
        n_bits: 32,
        min_distance: 2,
        max_distance: 3,
        use_2d: false,
        ..Default::default()
    })
    .unwrap();
    let arguments = || FingerprintFuncArguments {
        conf_id: 0,
        ..Default::default()
    };

    let sparse = generator
        .sparse_count_fingerprint(&molecule, &mut arguments())
        .unwrap();
    assert_eq!(
        sparse.nonzero_elements(),
        &map(&[(541_731, 1), (558_114, 1)])
    );
    let count = generator
        .count_fingerprint(&molecule, &mut arguments())
        .unwrap();
    assert_eq!(count.nonzero_elements(), &map(&[(28, 1), (30, 1)]));
    let explicit = generator.fingerprint(&molecule, &mut arguments()).unwrap();
    assert_eq!(explicit.on_bits(), &[16, 24]);
}
