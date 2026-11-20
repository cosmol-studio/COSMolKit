#![cfg(all(feature = "cap-fragmentation", feature = "cap-smiles"))]

use cosmolkit::{Molecule, SmilesParseParams, SubstructMatchError};

#[test]
fn ether_matches_retain_order_labels_and_both_propyl_occurrences() {
    let input = Molecule::from_smiles("CCCOCCC").unwrap();
    let sites = input.find_brics_bonds().unwrap();
    assert_eq!(
        sites
            .iter()
            .map(|m| (m.bond().index(), m.ends().map(|(a, e)| (a.index(), e))))
            .collect::<Vec<_>>(),
        [(2, [(3, 3), (2, 4)]), (3, [(3, 3), (4, 4)])]
    );
    let one = input.with_bonds_cut(&[sites[0].to_bond_cut()]).unwrap();
    assert_eq!(one.to_smiles().unwrap().as_bytes(), b"[3*]OCCC.[4*]CCC");
    let cuts: Vec<_> = sites.iter().map(|m| m.to_bond_cut()).collect();
    let cut = input.with_bonds_cut(&cuts).unwrap();
    let mut target = input.clone();
    target.cut_bonds_(&cuts).unwrap();
    assert_eq!(target, cut);
    let components: Vec<_> = cut
        .fragments()
        .unwrap()
        .into_iter()
        .map(|m| m.to_smiles().unwrap().as_bytes().to_owned())
        .collect();
    assert_eq!(
        components,
        [b"[4*]CCC".as_slice(), b"[3*]O[3*]", b"[4*]CCC"]
    );
    assert_eq!(input.to_smiles().unwrap().as_bytes(), b"CCCOCCC");
}

#[test]
fn first_matching_rule_wins_even_if_native_bulk_cut_uses_another_label() {
    let m = Molecule::from_smiles("COC12CCCN1C(=O)CC2").unwrap();
    let hits = m.find_brics_bonds().unwrap();
    assert_eq!(hits.len(), 1);
    assert_eq!(
        hits[0].ends().map(|(a, e)| (a.index(), e)),
        [(1, 3), (2, 13)]
    );
    assert!(
        Molecule::from_smiles("C1CCCCC1")
            .unwrap()
            .find_brics_bonds()
            .unwrap()
            .is_empty()
    );
    let alkene = Molecule::from_smiles("CC=CC").unwrap();
    let hits = alkene.find_brics_bonds().unwrap();
    assert_eq!(hits.len(), 1);
    assert_eq!(hits[0].ends().map(|(_, e)| e), [7, 7]);
    assert_eq!(
        alkene
            .with_bonds_cut(&[hits[0].to_bond_cut()])
            .unwrap()
            .to_smiles()
            .unwrap()
            .as_bytes(),
        b"[7*]=CC.[7*]=CC"
    );
}

#[test]
fn uninitialized_ring_state_is_an_error_even_without_an_eligible_bond() {
    let p = SmilesParseParams {
        sanitize: false,
        remove_hs: false,
        ..Default::default()
    };
    let carbon = Molecule::from_smiles_with_params("C", &p).unwrap();
    assert!(matches!(
        carbon.find_brics_bonds(),
        Err(SubstructMatchError::QueryContext(_))
    ));
}
