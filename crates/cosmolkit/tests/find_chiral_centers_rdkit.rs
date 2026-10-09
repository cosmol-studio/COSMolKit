#![cfg(all(feature = "cap-stereo", feature = "cap-smiles"))]
//! Fixed source-backed regressions for the modern chiral-center query.
use cosmolkit::{Molecule, SmilesParseParams, StereoReadError};

#[test]
fn modern_chiral_centers_match_pinned_rdkit_candidates_and_labels() {
    // RDKit 2026.03.1 / 351f8f378f8ad6bbd517980c38896e66bf907af8:
    // Chem.FindMolChiralCenters(..., useLegacyImplementation=False).
    // Ordinary fixed expectations: no oracle invocation or generated cache.
    let cases: &[(&str, bool, &[(usize, &str)])] = &[
        ("", false, &[]),
        ("", true, &[]),
        ("C", false, &[]),
        ("C", true, &[]),
        ("CC", false, &[]),
        ("CC", true, &[]),
        ("CCO", false, &[]),
        ("CCO", true, &[]),
        ("C[C@H](O)Cl", false, &[(1, "R")]),
        ("C[C@H](O)Cl", true, &[(1, "R")]),
        ("C[C@@H](O)Cl", false, &[(1, "S")]),
        ("C[C@@H](O)Cl", true, &[(1, "S")]),
        ("CC(O)Cl", false, &[]),
        ("CC(O)Cl", true, &[(1, "?")]),
        ("CC=CC", false, &[]),
        ("CC=CC", true, &[]),
        ("C/C=C/C", false, &[]),
        ("C/C=C/C", true, &[]),
        ("C/C=C\\C", false, &[]),
        ("C/C=C\\C", true, &[]),
        ("CC(=O)[O-]", false, &[]),
        ("CC(=O)[O-]", true, &[]),
        ("C[NH3+]", false, &[]),
        ("C[NH3+]", true, &[]),
        ("CC(=O)[O-].C[NH3+]", false, &[]),
        ("CC(=O)[O-].C[NH3+]", true, &[]),
        ("[CH-]1C=CC=C1", false, &[]),
        ("[CH-]1C=CC=C1", true, &[]),
        ("CN1CCC[C@H]1c1cccnc1", false, &[(5, "S")]),
        ("CN1CCC[C@H]1c1cccnc1", true, &[(5, "S")]),
        ("F[C@H](Cl)Br", false, &[(1, "R")]),
        ("F[C@H](Cl)Br", true, &[(1, "R")]),
        ("F[C@@H](Cl)Br", false, &[(1, "S")]),
        ("F[C@@H](Cl)Br", true, &[(1, "S")]),
        ("F[C@](Cl)(Br)I", false, &[(1, "S")]),
        ("F[C@](Cl)(Br)I", true, &[(1, "S")]),
        ("F[C@@](Cl)(Br)I", false, &[(1, "R")]),
        ("F[C@@](Cl)(Br)I", true, &[(1, "R")]),
        ("C1CC(C)C(C)C(C)C1", false, &[]),
        ("C1CC(C)C(C)C(C)C1", true, &[(2, "?"), (4, "?"), (6, "?")]),
        ("C1C[C@H](C)C(C)[C@H](C)C1", false, &[(2, "S"), (6, "R")]),
        (
            "C1C[C@H](C)C(C)[C@H](C)C1",
            true,
            &[(2, "S"), (4, "?"), (6, "R")],
        ),
        (
            "C1C[C@H](C)[C@H](C)[C@H](C)C1",
            false,
            &[(2, "S"), (4, "s"), (6, "R")],
        ),
        (
            "C1C[C@H](C)[C@H](C)[C@H](C)C1",
            true,
            &[(2, "S"), (4, "s"), (6, "R")],
        ),
        ("[13CH3][C@H]([12CH3])F", false, &[(1, "S")]),
        ("[13CH3][C@H]([12CH3])F", true, &[(1, "S")]),
        ("[13CH3]C([12CH3])F", false, &[]),
        ("[13CH3]C([12CH3])F", true, &[(1, "?")]),
        ("[2H][C@](F)(Cl)Br", false, &[(1, "S")]),
        ("[2H][C@](F)(Cl)Br", true, &[(1, "S")]),
        ("F[Pt@SP1](Cl)(Br)I", false, &[]),
        ("F[Pt@SP1](Cl)(Br)I", true, &[]),
        ("S[As@TB1](F)(Cl)(Br)N", false, &[]),
        ("S[As@TB1](F)(Cl)(Br)N", true, &[]),
        ("F[Co@OH1](Cl)(Br)(I)(N)O", false, &[]),
        ("F[Co@OH1](Cl)(Br)(I)(N)O", true, &[]),
        ("C[C@H](F)F", false, &[]),
        ("C[C@H](F)F", true, &[]),
        ("*[C@](F)(Cl)Br", false, &[(1, "S")]),
        ("*[C@](F)(Cl)Br", true, &[(1, "S")]),
        ("C[C@H](F)C(F)(Cl)Br", false, &[(1, "S")]),
        ("C[C@H](F)C(F)(Cl)Br", true, &[(1, "S"), (3, "?")]),
    ];
    for &(smiles, include, expected) in cases {
        let source = Molecule::from_smiles(smiles).unwrap();
        let peer = source.clone();
        let before = source.clone();
        let expected = expected
            .iter()
            .map(|&(i, label)| (i, label.to_owned()))
            .collect::<Vec<_>>();
        for _ in 0..2 {
            assert_eq!(
                source.find_chiral_centers(include).unwrap(),
                expected,
                "{smiles:?}, include_unassigned={include}"
            );
            assert_eq!(source, before);
            assert_eq!(peer, before);
            assert!(std::ptr::eq(source.topology(), peer.topology()));
            assert!(std::ptr::eq(source.properties(), peer.properties()));
        }
    }
}

#[test]
fn chiral_center_query_retains_typed_errors_and_shared_receivers() {
    use std::error::Error as _;
    let source = Molecule::from_smiles_with_params(
        "C[C@H]1CCCC[C@H]1C |atomProp:1._ringStereochemCand.malformed|",
        &SmilesParseParams {
            sanitize: false,
            remove_hs: false,
            skip_cleanup: true,
            ..SmilesParseParams::default()
        },
    )
    .unwrap();
    let peer = source.clone();
    for include in [false, true] {
        let error = source.find_chiral_centers(include).unwrap_err();
        assert!(
            matches!(
                error,
                StereoReadError::PotentialStereo(
                    cosmolkit::PotentialStereoError::InvalidPropertyKind { .. }
                )
            ),
            "{error:?}"
        );
        assert!(error.source().is_some());
        assert_eq!(source, peer);
        assert!(std::ptr::eq(source.topology(), peer.topology()));
    }
}

#[test]
fn chiral_center_query_recalculates_labels_without_changing_existing_assignment() {
    let source = Molecule::from_smiles("C1C[C@H](C)[C@H](C)[C@H](C)C1")
        .unwrap()
        .with_cip_labels()
        .unwrap();
    let peer = source.clone();
    assert!(source.cip_computed().unwrap());
    assert_eq!(
        source.find_chiral_centers(false).unwrap(),
        [(2, "S".into()), (4, "s".into()), (6, "R".into())]
    );
    assert_eq!(source, peer);
    assert!(source.cip_computed().unwrap());
    assert!(std::ptr::eq(source.topology(), peer.topology()));
}
