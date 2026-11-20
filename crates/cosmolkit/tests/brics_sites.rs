#![cfg(all(feature = "cap-fragmentation", feature = "cap-smiles"))]

use cosmolkit::{Molecule, SmilesParseParams, SubstructMatchError};

#[test]
fn terminal_products_are_unique_and_do_not_mutate_the_receiver() {
    for (smi, expected) in [
        ("C", vec!["C"]),
        ("C/C=C/C", vec!["[7*]CC"]),
        ("CCCOCCC", vec!["[3*]O[3*]", "[4*]CCC"]),
    ] {
        let source = Molecule::from_smiles(smi).unwrap();
        let before = source.to_smiles().unwrap();
        let mut products = source.decompose_brics().unwrap();
        products.sort();
        assert_eq!(products, expected);
        assert_eq!(source.to_smiles().unwrap(), before);
    }
}

#[test]
fn regular_r_brics_keeps_rule_orientation_and_recursive_products_distinct_from_cuts() {
    let decane = Molecule::from_smiles("CCCCCCCCCC").unwrap();
    let sites = decane.find_r_brics_bonds().unwrap();
    assert_eq!(
        sites
            .iter()
            .map(|s| s.ends().map(|(a, l)| (a.index(), l)))
            .collect::<Vec<_>>(),
        [
            [(3, 20), (4, 21)],
            [(6, 20), (5, 21)],
            [(2, 22), (1, 23)],
            [(7, 22), (8, 23)]
        ]
    );
    let cuts: Vec<_> = sites.iter().map(|s| s.to_bond_cut()).collect();
    let product = decane.with_bonds_cut(&cuts).unwrap();
    let mut inplace = decane.copy();
    inplace.cut_bonds_(&cuts).unwrap();
    assert_eq!(inplace.to_smiles().unwrap(), product.to_smiles().unwrap());
    let parts: Vec<_> = product
        .fragments()
        .unwrap()
        .into_iter()
        .map(|m| m.to_smiles().unwrap().into_bytes())
        .collect();
    assert_eq!(
        parts,
        [
            b"[23*]CC".as_slice(),
            b"[20*]CC[22*]",
            b"[21*]CC[21*]",
            b"[20*]CC[22*]",
            b"[23*]CC"
        ]
    );
    let mut terminal = decane.decompose_r_brics().unwrap();
    terminal.sort();
    assert_eq!(
        terminal,
        ["[20*]CCCC", "[20*]CCCC[21*]", "[21*]CC", "[21*]CC[21*]"]
    );
    assert_eq!(decane.decompose_brics().unwrap(), ["CCCCCCCCCC"]);

    let lactam = Molecule::from_smiles("O=C(CN1CCCCC1=O)c1ccco1").unwrap();
    let label = |ends: [(cosmolkit::AtomId, u8); 2]| ends.map(|(a, l)| (a.index(), l));
    assert!(
        lactam
            .find_r_brics_bonds()
            .unwrap()
            .iter()
            .any(|s| label(s.ends()) == [(2, 8), (3, 9)])
    );
    assert!(
        lactam
            .find_brics_bonds()
            .unwrap()
            .iter()
            .any(|s| label(s.ends()) == [(2, 8), (3, 10)])
    );
}

#[test]
fn regular_r_brics_discards_unsanitizable_sequences_instead_of_failing_the_query() {
    let mol = Molecule::from_smiles(
        "O=C([O-])C1=CC=C(N2N=C3[N]4=C2C5=CC=CC=C5[O-][Fe+3]4[O-]C6=CC=CC=C63)C=C1",
    )
    .unwrap();
    assert!(matches!(
        mol.decompose_brics(),
        Err(cosmolkit::BricsDecompositionError::UnkeyedProduct(_))
    ));
    let mut products = mol.decompose_r_brics().unwrap();
    products.sort();
    assert_eq!(
        products,
        [
            "[16*]c1ccc([16*])cc1",
            "[6*]C(=O)[O-]",
            "[9*]n1nc2-c3ccccc3[O-]->[Fe+3]3<-[O-]c4ccccc4-c1[n]->32"
        ]
    );
}

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
