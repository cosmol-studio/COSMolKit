#![cfg(feature = "descriptors")]
use cosmolkit::{DescriptorError, DescriptorReadError, Molecule, SmilesParseParams};

fn assert_rows(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for (i, (a, b)) in actual.iter().zip(expected).enumerate() {
        assert!((a - b).abs() < 1e-11, "row {i}: {a} vs {b}");
    }
}

#[test]
fn indices_follow_input_order_and_do_not_cross_components() {
    for (smiles, expected) in [
        ("CCO", &[121.0 / 72.0, 0.25, 545.0 / 72.0][..]),
        ("OCC", &[545.0 / 72.0, 0.25, 121.0 / 72.0]),
        ("[2H]OC", &[137.0 / 24.0, 3.5, 31.0 / 24.0]),
        ("*C", &[0.75, 2.25]),
        ("CO.C", &[1.0, 7.0, 0.0]),
    ] {
        let mol = Molecule::from_smiles(smiles).unwrap();
        assert_rows(&mol.estate_indices().unwrap(), expected);
    }
    let ethanol = Molecule::from_smiles("CCO").unwrap();
    let permuted = ethanol
        .with_atoms_renumbered(&[2, 0, 1].map(cosmolkit::AtomId::new))
        .unwrap()
        .with_assigned_valence()
        .unwrap();
    assert_rows(
        &permuted.estate_indices().unwrap(),
        &[545.0 / 72.0, 121.0 / 72.0, 0.25],
    );
}

#[test]
fn isolated_and_no_implicit_atoms_need_no_valence_assignment() {
    let raw = |smiles| {
        Molecule::from_smiles_with_params(
            smiles,
            &SmilesParseParams {
                sanitize: false,
                remove_hs: false,
                ..Default::default()
            },
        )
        .unwrap()
    };
    assert_rows(&raw("C.C").estate_indices().unwrap(), &[0.0, 0.0]);
    assert_rows(&raw("[CH3][NH2]").estate_indices().unwrap(), &[1.5, 4.5]);
    assert!(matches!(
        raw("CCO").estate_indices(),
        Err(DescriptorReadError::Algorithm {
            source: DescriptorError::Valence { .. }
        })
    ));
}

#[test]
fn signed_and_absolute_extrema_include_interior_atoms() {
    let mol = Molecule::from_smiles("CCC(=O)O").unwrap();
    assert_rows(
        &[
            mol.max_estate_index().unwrap(),
            mol.min_estate_index().unwrap(),
            mol.max_abs_estate_index().unwrap(),
            mol.min_abs_estate_index().unwrap(),
        ],
        &[1349.0 / 144.0, -161.0 / 216.0, 1349.0 / 144.0, 2.0 / 9.0],
    );
    let empty = Molecule::new();
    assert!(empty.estate_indices().unwrap().is_empty());
    for result in [
        empty.max_estate_index(),
        empty.min_estate_index(),
        empty.max_abs_estate_index(),
        empty.min_abs_estate_index(),
    ] {
        assert!(matches!(
            result,
            Err(DescriptorReadError::Algorithm {
                source: DescriptorError::EmptyEStateIndices
            })
        ));
    }
}

#[test]
fn the_two_vsa_vectors_bin_different_atom_properties() {
    let mol = Molecule::from_smiles("CCO").unwrap();
    let expected_areas = [
        0.,
        6.606881964512918,
        0.,
        0.,
        0.,
        6.923737199690624,
        0.,
        0.,
        5.106527394840706,
        0.,
        0.,
    ];
    let expected_indices = [0., 0., 545.0 / 72.0, 0., 0., 0., 0., 139.0 / 72.0, 0., 0.];
    assert_rows(&mol.estate_vsa().unwrap(), &expected_areas);
    assert_rows(&mol.vsa_estate().unwrap(), &expected_indices);

    let raw = Molecule::from_smiles_with_params(
        "CCO",
        &SmilesParseParams {
            sanitize: false,
            remove_hs: false,
            ..Default::default()
        },
    )
    .unwrap();
    let valid = raw.with_assigned_valence().unwrap();
    assert_rows(&valid.estate_vsa().unwrap(), &expected_areas);
    assert_rows(&valid.vsa_estate().unwrap(), &expected_indices);
    assert!(matches!(
        valid.num_rings(),
        Err(DescriptorReadError::MissingInitializedRings)
    ));

    let empty = Molecule::new();
    assert_eq!(empty.estate_vsa().unwrap(), [0.0; 11]);
    assert_eq!(empty.vsa_estate().unwrap(), [0.0; 10]);
}

#[test]
fn numbered_vsa_results_are_per_direction_per_molecule_and_per_graph() {
    let mut mol = Molecule::from_smiles("CCO").unwrap();
    let no_h = mol
        .labute_asa_contributions_with_params(false, true)
        .unwrap();
    let h = mol
        .labute_asa_contributions_with_params(true, true)
        .unwrap();
    assert_ne!(no_h.atom_contributions, h.atom_contributions);
    let labute_rows = |m: &Molecule| m.labute_asa_contributions().unwrap().atom_contributions;

    let held_indices = mol.estate_indices().unwrap();
    let area_bins = mol.estate_vsa().unwrap();
    mol.labute_asa_with_params(false, true).unwrap();
    let copy = mol.copy();
    assert_eq!(mol.estate_vsa_6().unwrap(), area_bins[5]);
    assert_eq!(labute_rows(&mol), no_h.atom_contributions);
    assert_eq!(copy.estate_vsa_6().unwrap(), area_bins[5]);
    assert_eq!(labute_rows(&copy), h.atom_contributions);
    assert_eq!(labute_rows(&mol), no_h.atom_contributions);

    mol.vsa_estate_3().unwrap();
    assert_eq!(labute_rows(&mol), h.atom_contributions);
    mol.labute_asa_with_params(false, true).unwrap();
    mol.vsa_estate_3().unwrap();
    assert_eq!(labute_rows(&mol), no_h.atom_contributions);
    mol.vsa_estate().unwrap();
    assert_eq!(labute_rows(&mol), h.atom_contributions);

    mol.remove_atom_(cosmolkit::AtomId::new(1)).unwrap();
    assert_rows(&held_indices, &[121.0 / 72.0, 0.25, 545.0 / 72.0]);
    assert_eq!(mol.estate_indices().unwrap(), [0.0, 0.0]);
    assert_eq!(mol.estate_vsa_6().unwrap(), 0.0);
}

#[test]
fn kekulization_moves_benzene_across_the_surface_area_cutoff() {
    let mut aromatic = Molecule::from_smiles("c1ccccc1").unwrap();
    assert_rows(&[aromatic.vsa_estate_6().unwrap()], &[12.0]);

    let kekule = aromatic.with_kekulized_bonds().unwrap();
    assert_eq!(kekule.vsa_estate_6().unwrap(), 0.0);
    assert_rows(&[kekule.vsa_estate_7().unwrap()], &[12.0]);
    assert_rows(&[aromatic.vsa_estate_6().unwrap()], &[12.0]);

    aromatic.kekulize_bonds_().unwrap();
    assert_eq!(aromatic.vsa_estate_6().unwrap(), 0.0);
    assert_rows(&[aromatic.vsa_estate_7().unwrap()], &[12.0]);
}
