use std::collections::BTreeSet;

use super::{
    FingerprintError, extend_paths, find_all_paths_of_length_n, find_all_paths_of_lengths_m_to_n,
    rdkit_fp_bond_between_atoms,
};
use crate::{AtomSpec, BondOrder, BondSpec, Element, Molecule};

fn bond_membership(molecule: &Molecule, atom_path: &[usize]) -> Vec<bool> {
    let mut membership = vec![false; molecule.num_bonds()];
    for pair in atom_path.windows(2) {
        let bond = rdkit_fp_bond_between_atoms(molecule, pair[0], pair[1])
            .expect("enumerated path edge must exist");
        membership[bond] = true;
    }
    membership
}

#[test]
fn chain_paths_follow_atom_index_order_and_remove_reverse_bond_duplicates() {
    let molecule = Molecule::from_smiles("CCCC").expect("chain");

    assert_eq!(
        find_all_paths_of_length_n(&molecule, 3, false, false, -1, false).unwrap(),
        vec![vec![0, 1, 2], vec![1, 2, 3]]
    );
    assert_eq!(
        find_all_paths_of_length_n(&molecule, 4, false, false, -1, false).unwrap(),
        vec![vec![0, 1, 2, 3]]
    );
    assert_eq!(
        find_all_paths_of_length_n(&molecule, 3, true, false, -1, false).unwrap(),
        vec![vec![0, 1, 2]]
    );
}

#[test]
fn branch_paths_preserve_source_adjacency_scan_order() {
    let molecule = Molecule::from_smiles("CC(C)C").expect("branched molecule");

    assert_eq!(
        find_all_paths_of_length_n(&molecule, 3, false, false, -1, false).unwrap(),
        vec![vec![0, 1, 2], vec![0, 1, 3], vec![2, 1, 3]]
    );
}

#[test]
fn ring_paths_allow_only_the_final_closure_and_deduplicate_by_bond_set() {
    let molecule = Molecule::from_smiles("C1CCC1").expect("four-membered ring");
    let open_paths = find_all_paths_of_length_n(&molecule, 4, false, false, -1, false).unwrap();
    let closed_paths = find_all_paths_of_length_n(&molecule, 5, false, false, -1, false).unwrap();

    assert_eq!(open_paths.len(), 4);
    assert_eq!(open_paths[0], vec![0, 1, 2, 3]);
    assert_eq!(closed_paths, vec![vec![0, 1, 2, 3, 0]]);
    for path in &open_paths {
        assert_eq!(
            path.iter().copied().collect::<BTreeSet<_>>().len(),
            path.len()
        );
    }
    assert_eq!(
        closed_paths[0][..closed_paths[0].len() - 1]
            .iter()
            .copied()
            .collect::<BTreeSet<_>>()
            .len(),
        closed_paths[0].len() - 1
    );
}

#[test]
fn extend_paths_rejects_internal_repeats_but_accepts_a_final_ring_closure() {
    let triangle_adjacency = [0, 1, 1, 1, 0, 1, 1, 1, 0];
    let path = vec![vec![0, 1, 2]];

    assert!(
        extend_paths(&triangle_adjacency, 3, &path, 5, None)
            .unwrap()
            .is_empty()
    );
    assert_eq!(
        extend_paths(&triangle_adjacency, 3, &path, 4, None).unwrap(),
        vec![vec![0, 1, 2, 0]]
    );
}

#[test]
fn explicit_hydrogens_are_included_only_when_requested() {
    let mut builder = Molecule::builder();
    let carbon = builder.add_atom(AtomSpec::new(Element::C));
    for _ in 0..4 {
        let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
        builder
            .add_bond(BondSpec::new(carbon, hydrogen, BondOrder::Single))
            .expect("methane C-H bond");
    }
    let molecule = builder.build().expect("explicit methane graph");
    assert_eq!(molecule.num_atoms(), 5);

    assert!(
        find_all_paths_of_length_n(&molecule, 2, false, false, -1, false)
            .unwrap()
            .is_empty()
    );
    let with_hydrogens = find_all_paths_of_length_n(&molecule, 2, false, true, -1, false).unwrap();
    assert_eq!(with_hydrogens.len(), 4);
    assert!(with_hydrogens.iter().all(|path| {
        path.iter()
            .any(|atom| molecule.atoms()[*atom].atomic_number() == 1)
    }));
}

#[test]
fn disconnected_and_rooted_paths_keep_source_seed_behavior() {
    let disconnected = Molecule::from_smiles("CC.CC").expect("disconnected molecule");
    assert_eq!(
        find_all_paths_of_length_n(&disconnected, 2, false, false, -1, false).unwrap(),
        vec![vec![0, 1], vec![2, 3]]
    );

    let chain = Molecule::from_smiles("CCCCC").expect("five-atom chain");
    assert_eq!(
        find_all_paths_of_length_n(&chain, 3, false, false, 2, false).unwrap(),
        vec![vec![2, 1, 0], vec![2, 3, 4]]
    );
    assert!(
        find_all_paths_of_length_n(&chain, 3, false, false, 99, false)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn target_length_boundaries_and_range_keys_match_source_behavior() {
    let molecule = Molecule::from_smiles("CCC").expect("chain");
    assert_eq!(
        find_all_paths_of_length_n(&molecule, 1, false, false, -1, false).unwrap(),
        vec![vec![0], vec![1], vec![2]]
    );
    assert!(
        find_all_paths_of_length_n(&molecule, 0, false, false, -1, false)
            .unwrap()
            .is_empty()
    );

    let range = find_all_paths_of_lengths_m_to_n(&molecule, 1, 3, false, false, -1, false).unwrap();
    assert_eq!(range.keys().copied().collect::<Vec<_>>(), vec![1, 2, 3]);
    assert!(matches!(
        find_all_paths_of_lengths_m_to_n(&molecule, 3, 2, false, false, -1, false),
        Err(FingerprintError::InvalidArguments { .. })
    ));
}

#[test]
fn fused_ring_paths_are_legal_deterministic_and_bond_set_unique() {
    let molecule = Molecule::from_smiles("C1CCC2CCCCC2C1").expect("fused rings");
    let paths = find_all_paths_of_length_n(&molecule, 6, false, false, -1, false).unwrap();
    assert!(!paths.is_empty());
    assert_eq!(
        paths,
        find_all_paths_of_length_n(&molecule, 6, false, false, -1, false).unwrap()
    );

    let memberships = paths
        .iter()
        .map(|path| bond_membership(&molecule, path))
        .collect::<BTreeSet<_>>();
    assert_eq!(memberships.len(), paths.len());
    assert!(paths.iter().all(|path| {
        path.windows(2)
            .all(|pair| rdkit_fp_bond_between_atoms(&molecule, pair[0], pair[1]).is_some())
    }));
}

#[test]
fn only_shortest_paths_prunes_longer_routes_between_adjacent_ring_atoms() {
    let molecule = Molecule::from_smiles("C1CC1").expect("triangle");
    assert_eq!(
        find_all_paths_of_length_n(&molecule, 3, false, false, -1, false)
            .unwrap()
            .len(),
        3
    );
    assert!(
        find_all_paths_of_length_n(&molecule, 3, false, false, -1, true)
            .unwrap()
            .is_empty()
    );
}
