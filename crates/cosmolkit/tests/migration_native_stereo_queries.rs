use cosmolkit::{AtomId, LigandRef, Molecule};

#[test]
fn potential_centers_and_lone_pairs_follow_modern_perception() {
    // RDKit 2026.09.1 FindPotentialStereo: unassigned lactate, redundant
    // tagged substituents, pyramidal sulfur and two-neighbor phosphine.
    for (text, expected) in [
        ("CC(O)C(=O)O", vec![(1, false)]),
        ("C[C@](C)(F)Cl", vec![]),
        ("C[S@](=O)CC", vec![(1, true)]),
        ("C[P@H]CC", vec![(1, true)]),
    ] {
        let mol = Molecule::from_smiles(text).unwrap();
        let before = mol.to_binary().unwrap();
        let values = mol.tetrahedral_stereo(true).unwrap();
        assert_eq!(
            values
                .iter()
                .map(|s| (s.center.index(), s.ligands.is_some()))
                .collect::<Vec<_>>(),
            expected
        );
        assert_eq!(
            mol.tetrahedral_stereo(false).unwrap(),
            values
                .iter()
                .filter(|s| s.ligands.is_some())
                .cloned()
                .collect::<Vec<_>>()
        );
        assert_eq!(
            values.iter().map(|s| s.center.index()).collect::<Vec<_>>(),
            mol.find_chiral_centers(true)
                .unwrap()
                .into_iter()
                .map(|s| s.0)
                .collect::<Vec<_>>()
        );
        assert_eq!(mol.to_binary().unwrap(), before);
    }
    let phosphine = Molecule::from_smiles("C[P@H]CC")
        .unwrap()
        .tetrahedral_stereo(true)
        .unwrap();
    assert_eq!(
        phosphine[0].ligands.unwrap(),
        [
            LigandRef::Atom(AtomId::new(0)),
            LigandRef::Atom(AtomId::new(2)),
            LigandRef::ImplicitHydrogen,
            LigandRef::LonePair,
        ]
    );
}

#[test]
fn tetrahedral_stereo_ligand_order_encodes_smiles_handedness() {
    let ccw = Molecule::from_smiles("F[C@](Cl)(Br)I").expect("parse chiral SMILES");
    let cw = Molecule::from_smiles("F[C@@](Cl)(Br)I").expect("parse chiral SMILES");
    let ccw_stereo = ccw.tetrahedral_stereo(true).expect("tetrahedral stereo");
    let cw_stereo = cw.tetrahedral_stereo(true).expect("tetrahedral stereo");
    assert_eq!(ccw_stereo.len(), 1);
    assert_eq!(cw_stereo.len(), 1);
    assert_eq!(ccw_stereo[0].center, AtomId::new(1));
    assert_eq!(cw_stereo[0].center, AtomId::new(1));
    assert_eq!(
        ccw_stereo[0].ligands.unwrap(),
        [
            LigandRef::Atom(AtomId::new(0)),
            LigandRef::Atom(AtomId::new(2)),
            LigandRef::Atom(AtomId::new(3)),
            LigandRef::Atom(AtomId::new(4))
        ]
    );
    assert_eq!(
        cw_stereo[0].ligands.unwrap(),
        [
            LigandRef::Atom(AtomId::new(0)),
            LigandRef::Atom(AtomId::new(2)),
            LigandRef::Atom(AtomId::new(4)),
            LigandRef::Atom(AtomId::new(3))
        ]
    );
}

#[test]
fn tetrahedral_stereo_places_implicit_hydrogen_as_fourth_ligand() {
    let mol = Molecule::from_smiles("[13CH3:7][C@H](F)Cl").expect("parse chiral SMILES");
    let stereo = mol.tetrahedral_stereo(true).expect("tetrahedral stereo");
    assert_eq!(stereo.len(), 1);
    assert_eq!(stereo[0].center, AtomId::new(1));
    assert_eq!(
        stereo[0].ligands.unwrap(),
        [
            LigandRef::Atom(AtomId::new(0)),
            LigandRef::Atom(AtomId::new(2)),
            LigandRef::Atom(AtomId::new(3)),
            LigandRef::ImplicitHydrogen
        ]
    );
}

#[test]
fn native_perception_and_modern_chiral_centers_retain_read_only_receivers() {
    let molecule = Molecule::from_smiles("F[C@H](Cl)Br").unwrap();
    let before = molecule.to_smiles().unwrap();
    molecule.perceive_stereochemistry().unwrap();
    // Pinned RDKit 2026.03.1 FindMolChiralCenters modern branch returns
    // actual tetrahedral centers and CIP labels, not every unspecified atom.
    let labels = molecule.find_chiral_centers(true).unwrap();
    assert_eq!(labels, vec![(1, "R".to_owned())]);
    assert_eq!(molecule.find_chiral_centers(false).unwrap(), labels);
    assert_eq!(molecule.to_smiles().unwrap(), before);
    assert!(Molecule::new().tetrahedral_stereo(true).unwrap().is_empty());
    Molecule::new().perceive_stereochemistry().unwrap();
}
