use cosmolkit::{AtomId, LigandRef, Molecule};

#[test]
fn tetrahedral_stereo_ligand_order_encodes_smiles_handedness() {
    let ccw = Molecule::from_smiles("F[C@](Cl)(Br)I").expect("parse chiral SMILES");
    let cw = Molecule::from_smiles("F[C@@](Cl)(Br)I").expect("parse chiral SMILES");
    let ccw_stereo = ccw.tetrahedral_stereo().expect("tetrahedral stereo");
    let cw_stereo = cw.tetrahedral_stereo().expect("tetrahedral stereo");
    assert_eq!(ccw_stereo.len(), 1);
    assert_eq!(cw_stereo.len(), 1);
    assert_eq!(ccw_stereo[0].center, AtomId::new(1));
    assert_eq!(cw_stereo[0].center, AtomId::new(1));
    assert_eq!(
        ccw_stereo[0].ligands,
        [
            LigandRef::Atom(AtomId::new(0)),
            LigandRef::Atom(AtomId::new(2)),
            LigandRef::Atom(AtomId::new(3)),
            LigandRef::Atom(AtomId::new(4))
        ]
    );
    assert_eq!(
        cw_stereo[0].ligands,
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
    let stereo = mol.tetrahedral_stereo().expect("tetrahedral stereo");
    assert_eq!(stereo.len(), 1);
    assert_eq!(stereo[0].center, AtomId::new(1));
    assert_eq!(
        stereo[0].ligands,
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
    assert!(Molecule::new().tetrahedral_stereo().unwrap().is_empty());
    Molecule::new().perceive_stereochemistry().unwrap();
}
