//! Original CCCO conditions mapped to the existing canonical owned values.
#![cfg(feature = "cap-fingerprints")]
use cosmolkit::{AtomPairFingerprintParams, Fingerprint, FingerprintAdditionalOutput, Molecule};
use std::collections::BTreeMap;

const BITS: [u32; 6] = [624, 1144, 1336, 1337, 1404, 1596];
type Snapshot = (
    Option<Vec<u32>>,
    Option<Vec<Vec<u64>>>,
    Option<BTreeMap<u64, Vec<(u32, u32)>>>,
    Option<BTreeMap<u64, Vec<Vec<i32>>>>,
    Option<BTreeMap<u64, Vec<Vec<i32>>>>,
);
fn snapshot(output: &FingerprintAdditionalOutput) -> Snapshot {
    (
        output.atom_counts().map(<[u32]>::to_vec),
        output.atom_to_bits().map(<[Vec<u64>]>::to_vec),
        output.bit_info_map().cloned(),
        output.atoms_per_bit().cloned(),
        output.bit_paths().cloned(),
    )
}
fn expected() -> Snapshot {
    let pairs: BTreeMap<u64, Vec<(u32, u32)>> = [
        (624, vec![(0, 3)]),
        (1144, vec![(0, 2)]),
        (1336, vec![(1, 2), (1, 3)]),
        (1337, vec![(1, 2), (1, 3)]),
        (1404, vec![(0, 1)]),
        (1596, vec![(2, 3)]),
    ]
    .into();
    let atoms = pairs
        .iter()
        .map(|(&bit, pairs)| {
            (
                bit,
                pairs
                    .iter()
                    .map(|&(a, b)| vec![a as i32, b as i32])
                    .collect(),
            )
        })
        .collect();
    (
        Some(vec![3, 3, 3, 3]),
        Some(vec![
            vec![624, 1144, 1404],
            vec![1336, 1337, 1404],
            vec![1144, 1336, 1337, 1596],
            vec![624, 1336, 1337, 1596],
        ]),
        Some(pairs),
        Some(atoms),
        None,
    )
}
fn collect() -> (Fingerprint, FingerprintAdditionalOutput) {
    let mut output = FingerprintAdditionalOutput::new();
    output.allocate_atom_counts();
    output.allocate_atom_to_bits();
    output.allocate_bit_info_map();
    output.allocate_atoms_per_bit();
    let fingerprint = Molecule::from_smiles("CCCO")
        .unwrap()
        .atom_pair_fingerprint_with_params(&AtomPairFingerprintParams::default(), Some(&mut output))
        .unwrap();
    (fingerprint, output)
}
#[test]
fn fingerprint_owned_clone_and_detached_bits_retain_original_conditions() {
    let (first, _) = collect();
    let second = first.clone();
    let mut bits = second.on_bits();
    bits.clear();
    assert_eq!(first.n_bits(), 2048);
    assert_eq!(first.on_bits(), BITS);
    assert_eq!(second.on_bits(), BITS);
}
#[test]
fn full_additional_output_snapshot_mutations_do_not_escape() {
    let (_, output) = collect();
    assert_eq!(snapshot(&output), expected());
    let mut detached = snapshot(&output);
    detached.0.as_mut().unwrap()[0] = 0;
    detached.1.as_mut().unwrap()[0].clear();
    detached.2.as_mut().unwrap().get_mut(&624).unwrap().clear();
    detached.3.as_mut().unwrap().get_mut(&624).unwrap()[0].clear();
    assert_eq!(snapshot(&output), expected());
}
#[test]
fn width_and_collection_information_is_observable_from_canonical_values() {
    let (fingerprint, output) = collect();
    assert_eq!(fingerprint.n_bits(), 2048);
    assert!(output.atom_counts().is_some());
    assert!(output.atom_to_bits().is_some());
    assert!(output.bit_info_map().is_some());
    assert!(output.atoms_per_bit().is_some());
    assert!(output.bit_paths().is_none());
    assert_eq!(
        snapshot(&FingerprintAdditionalOutput::new()),
        (None, None, None, None, None)
    );
    let plain = Molecule::from_smiles("CCCO")
        .unwrap()
        .atom_pair_fingerprint_with_params(&AtomPairFingerprintParams::default(), None)
        .unwrap();
    assert_eq!(plain.on_bits(), BITS);
}
#[test]
fn all_five_optional_fields_and_original_fingerprint_survive_output_reuse() {
    let (first, mut output) = collect();
    let detached = snapshot(&output);
    let second = Molecule::from_smiles("CCO")
        .unwrap()
        .atom_pair_fingerprint_with_params(&AtomPairFingerprintParams::default(), Some(&mut output))
        .unwrap();
    assert_eq!(output.atom_counts(), Some([2, 2, 2].as_slice()));
    assert_eq!(output.atom_to_bits().unwrap().len(), 3);
    assert_eq!(detached, expected());
    assert_eq!(first.on_bits(), BITS);
    assert_ne!(second.on_bits(), first.on_bits());
}
