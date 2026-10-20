#![cfg(feature = "cap-inchi")]

use cosmolkit::{AtomId, BondStereo, Hybridization, InchiErrorKind, InchiReadParams, Molecule};

#[test]
fn inchi_read_keeps_sanitize_errors_through_hydrogen_removal() {
    // RDKit InchiToMol rethrows MolSanitizeException from both removeHs and
    // sanitizeMol: this graph produces sulfur valence 7 at atom 15 in each.
    let identifier =
        "InChI=1S/C8H16O6S2/c9-5-8(14-16(11,12)13)7(10)6-15-3-1-2-4-15/h7-10H,1-6H2/t7-,8+/m1/s1";
    for remove_hs in [false, true] {
        let error =
            Molecule::from_inchi_with_params(identifier, &InchiReadParams::new(true, remove_hs))
                .unwrap_err();
        assert_eq!(error.kind, InchiErrorKind::SanitizeFailed);
        assert!(
            error
                .detail
                .contains("Explicit valence for atom # 15 S, 7, is greater than permitted")
        );
    }
}

#[test]
fn inchi_read_preserves_sanitized_atom_and_bond_state() {
    // RDKit 2026.09.1 InchiToMol -> sanitizeMol/removeHs -> assignStereochemistry.
    let identifier = "InChI=1S/C4H6/c1-3-4-2/h3-4H,1-2H2";
    for sanitize in [false, true] {
        for remove_hs in [false, true] {
            let mol = Molecule::from_inchi_with_params(
                identifier,
                &InchiReadParams::new(sanitize, remove_hs),
            )
            .unwrap();
            assert_eq!(mol.num_atoms(), 4);
            assert_eq!(mol.num_bonds(), 3);
            let hybridization = if sanitize {
                Hybridization::Sp2
            } else {
                Hybridization::Unspecified
            };
            assert!(
                mol.atoms()
                    .iter()
                    .all(|a| a.hybridization() == hybridization)
            );
            assert!(mol.bonds().iter().all(|b| b.is_conjugated() == sanitize));
        }
    }
    let mol = Molecule::from_inchi(identifier).unwrap();
    assert!(
        mol.atoms()
            .iter()
            .all(|a| a.hybridization() == Hybridization::Sp2)
    );
    assert!(mol.bonds().iter().all(|b| b.is_conjugated()));
}

#[test]
fn inchi_unsanitized_stereo_uses_source_fast_ring_membership() {
    // Fixed boundary regression: legacyStereoPerception uses fastFindRings
    // without sanitization. SSSR puts bond 20 in a six-member ring; DFS puts
    // it in a fourteen-member ring, so only the latter retains its E stereo.
    // Expected source state: RDKit 2026.09.1, InchiToMol and Chirality.cpp.
    let identifier = "InChI=1S/C28H25NO4/c1-17-14-24(31)33-27-21(17)11-12-22-25(27)19(10-13-23(30)18-8-6-5-7-9-18)15-20(26(22)32)16-29-28(2,3)4/h5-16,29H,1-4H3/b13-10+,20-16-";
    for sanitize in [false, true] {
        for remove_hs in [false, true] {
            let mol = Molecule::from_inchi_with_params(
                identifier,
                &InchiReadParams::new(sanitize, remove_hs),
            )
            .unwrap();
            let actual: Vec<_> = mol
                .bonds()
                .iter()
                .enumerate()
                .filter(|(_, b)| b.stereo() != BondStereo::None)
                .map(|(i, b)| (i, b.stereo(), b.stereo_atom_references().to_vec()))
                .collect();
            // removeHs also compacts bond indices; these are the actual source
            // rows after each path, not indices borrowed from the raw graph.
            let expected = if !sanitize {
                vec![
                    (11, BondStereo::E, [18, 22]),
                    (20, BondStereo::E, [19, 24]),
                    (23, BondStereo::Z, [28, 25]),
                ]
            } else if remove_hs {
                vec![(10, BondStereo::E, [18, 22]), (20, BondStereo::Z, [28, 25])]
            } else {
                vec![(11, BondStereo::E, [18, 22]), (23, BondStereo::Z, [28, 25])]
            };
            let expected: Vec<_> = expected
                .into_iter()
                .map(|(i, stereo, atoms)| (i, stereo, atoms.map(AtomId::new).to_vec()))
                .collect();
            assert_eq!(
                actual, expected,
                "sanitize={sanitize}, remove_hs={remove_hs}"
            );
        }
    }
}
