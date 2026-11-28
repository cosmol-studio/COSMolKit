#![cfg(feature = "descriptors")]

use cosmolkit::{
    DescriptorError, DescriptorReadError, GasteigerChargeError, GasteigerChargeParams, Molecule,
    SmilesParseParams,
};

fn parse(smiles: &str) -> Molecule {
    Molecule::from_smiles(smiles).unwrap()
}

#[test]
fn charge_rows_follow_hydrogen_representation_and_supplied_formal_charge() {
    let ethanol = parse("CCO");
    let rows = ethanol.gasteiger_charges().unwrap();
    assert_eq!(rows[2].atom_charge(), -0.3966637147812598);
    assert_eq!(rows[2].attached_hydrogen_charge(), 0.21002230656786322);

    let graph_hs = Molecule::from_smiles_with_params(
        "[H]OC([H])([H])C([H])([H])[H]",
        &SmilesParseParams {
            remove_hs: false,
            ..Default::default()
        },
    )
    .unwrap()
    .gasteiger_charges()
    .unwrap();
    assert_eq!(graph_hs.len(), 9);
    assert!((graph_hs[0].atom_charge() - rows[2].attached_hydrogen_charge()).abs() < 1e-14);
    assert!((graph_hs[1].atom_charge() - rows[2].atom_charge()).abs() < 1e-14);
    assert_eq!(graph_hs[1].attached_hydrogen_charge(), 0.0);
    assert_eq!(
        parse("CN").gasteiger_charges().unwrap(),
        parse("C[NH2]").gasteiger_charges().unwrap()
    );

    let bracket = parse("C[NH2]").gasteiger_charges().unwrap();
    assert!((bracket[1].atom_charge() - -0.33326905).abs() < 1e-7);
    assert!((bracket[1].attached_hydrogen_charge() - 0.23619516).abs() < 1e-7);
    let total: f64 = rows
        .iter()
        .map(|r| r.atom_charge() + r.attached_hydrogen_charge())
        .sum();
    assert!(total.abs() < 2e-12);
    let acetate = parse("CC(=O)[O-]");
    let zero = GasteigerChargeParams::new(0, false).unwrap();
    let initial = acetate.gasteiger_charges_with_params(&zero).unwrap();
    assert_eq!(
        initial
            .iter()
            .map(|row| row.atom_charge())
            .collect::<Vec<_>>(),
        [0.0, 0.0, -0.5, -0.5]
    );
    assert!(
        initial
            .iter()
            .all(|row| row.attached_hydrogen_charge() == 0.0)
    );
    let total: f64 = acetate
        .gasteiger_charges()
        .unwrap()
        .into_iter()
        .map(|row| row.atom_charge() + row.attached_hydrogen_charge())
        .sum();
    assert!((total + 1.0).abs() < 1e-13);
    assert!(
        (parse("O=S(O)(F)(F)F").gasteiger_charges().unwrap()[1].atom_charge() - 0.212415125).abs()
            < 1e-9
    );
}

#[test]
fn reads_do_not_fill_valence_or_hybridization_and_strict_lookup_stays_strict() {
    let raw = Molecule::from_smiles_with_params(
        "CCO",
        &SmilesParseParams {
            sanitize: false,
            remove_hs: false,
            ..Default::default()
        },
    )
    .unwrap();
    let strict_zero = GasteigerChargeParams::new(0, true).unwrap();
    assert!(matches!(raw.gasteiger_charges_with_params(&strict_zero),
        Err(DescriptorReadError::Algorithm { source: DescriptorError::GasteigerCharge {
            source: GasteigerChargeError::MissingParameters { element, mode }
        } }) if element == "C" && mode.is_empty()));
    let zero = GasteigerChargeParams::new(0, false).unwrap();
    let initial = raw.gasteiger_charges_with_params(&zero).unwrap();
    assert!(
        initial
            .iter()
            .all(|row| row.atom_charge().to_bits() == 0.0_f64.to_bits()
                && row.attached_hydrogen_charge().to_bits() == 0.0_f64.to_bits())
    );
    assert!(matches!(
        raw.gasteiger_charges(),
        Err(DescriptorReadError::Algorithm {
            source: DescriptorError::GasteigerCharge {
                source: GasteigerChargeError::Valence {
                    source: cosmolkit::ValenceError::ImplicitValenceCacheNotInitialized { .. }
                }
            }
        })
    ));

    let only_valence = raw.with_assigned_valence().unwrap();
    assert!(
        only_valence
            .gasteiger_charges()
            .unwrap()
            .iter()
            .all(|row| row.atom_charge().is_nan() && row.attached_hydrogen_charge().is_nan())
    );
    let strict = GasteigerChargeParams::new(12, true).unwrap();
    assert!(
        matches!(only_valence.gasteiger_charges_with_params(&strict), Err(DescriptorReadError::Algorithm { source: DescriptorError::GasteigerCharge { source: GasteigerChargeError::MissingParameters {element, mode} } }) if element == "C" && mode.is_empty())
    );
    let ion = parse("[Na+]");
    assert_eq!(ion.gasteiger_charges().unwrap()[0].atom_charge(), 1.0);
    assert!(ion.gasteiger_charges_with_params(&strict).is_err());
}

#[test]
fn extrema_and_vsa_preserve_finite_components_next_to_nan() {
    let ethanol = parse("CCO");
    let finite = ethanol.gasteiger_charges().unwrap();
    let ethanol_area = ethanol.peoe_vsa().unwrap();
    for (smiles, bad, good) in [("[Fe]C.CCO", 0, 2), ("CCO.[Fe]C", 3, 0)] {
        let mol = parse(smiles);
        let rows = mol.gasteiger_charges().unwrap();
        assert_eq!(&rows[good..good + 3], &finite);
        assert!(rows[bad].atom_charge().is_nan());
        assert_eq!(rows[bad].attached_hydrogen_charge(), 0.0);
        assert!(
            rows[bad + 1].atom_charge().is_nan()
                && rows[bad + 1].attached_hydrogen_charge().is_nan()
        );
        // Descriptors.py folds min(chg, minChg)/max(chg, maxChg):
        // trailing finite atoms clear preceding NaNs; trailing NaNs persist.
        let extrema = [
            mol.min_partial_charge(),
            mol.max_partial_charge(),
            mol.min_abs_partial_charge(),
            mol.max_abs_partial_charge(),
        ]
        .map(Result::unwrap);
        if good > bad {
            let expected = [
                -0.3966637147812598_f64,
                0.04022058174908642,
                0.04022058174908642,
                0.3966637147812598,
            ];
            assert_eq!(extrema.map(f64::to_bits), expected.map(f64::to_bits));
        } else {
            assert!(extrema.into_iter().all(f64::is_nan));
        }
        let bins = mol.peoe_vsa().unwrap();
        assert_eq!(&bins[..13], &ethanol_area[..13]);
        assert!(bins[13] > ethanol_area[13]);
        assert_eq!(mol.peoe_vsa_14().unwrap(), bins[13]);
    }
    let acid = parse("CC(=O)O");
    assert_eq!(
        acid.min_abs_partial_charge().unwrap(),
        acid.max_partial_charge().unwrap().abs()
    );
    assert!((acid.min_partial_charge().unwrap() - -0.481434).abs() < 1e-6);
    assert!((acid.max_partial_charge().unwrap() - 0.299685).abs() < 1e-6);
    assert!((acid.max_abs_partial_charge().unwrap() - 0.481434).abs() < 1e-6);
    let empty = Molecule::new();
    assert_eq!(empty.max_partial_charge().unwrap(), -500.0);
    assert_eq!(empty.min_partial_charge().unwrap(), 500.0);
    assert_eq!(empty.min_abs_partial_charge().unwrap(), 500.0);
    assert_eq!(empty.max_abs_partial_charge().unwrap(), 500.0);
    let max_iterations = GasteigerChargeParams::new(i32::MAX.into(), true).unwrap();
    assert!(
        empty
            .gasteiger_charges_with_params(&max_iterations)
            .unwrap()
            .is_empty()
    );
    assert_eq!(empty.peoe_vsa().unwrap(), [0.0; 14]);
}

#[test]
fn charge_bins_read_existing_labute_atom_rows() {
    let mol = parse("CCO");
    let unringed = mol
        .with_atoms_renumbered(&[0, 1, 2].map(cosmolkit::AtomId::new))
        .unwrap()
        .with_assigned_valence()
        .unwrap();
    assert!(matches!(
        unringed.num_rings(),
        Err(DescriptorReadError::MissingInitializedRings)
    ));
    assert_eq!(
        unringed.gasteiger_charges().unwrap(),
        mol.gasteiger_charges().unwrap()
    );
    assert_eq!(unringed.estate_vsa().unwrap(), mol.estate_vsa().unwrap());
    assert_eq!(unringed.vsa_estate().unwrap(), mol.vsa_estate().unwrap());
    let cold = mol.peoe_vsa().unwrap();
    assert_eq!(cold.len(), 14);
    assert!((cold.iter().sum::<f64>() - 18.63714656).abs() < 1e-7);
    assert_eq!(unringed.peoe_vsa().unwrap(), cold);
    mol.labute_asa_with_params(false, true).unwrap();
    let warm = mol.peoe_vsa().unwrap();
    assert!((warm.iter().sum::<f64>() - 18.68276248).abs() < 1e-7);
    assert_eq!(
        cold.iter()
            .zip(&warm)
            .enumerate()
            .filter_map(|(i, (a, b))| (a != b).then_some(i + 1))
            .collect::<Vec<_>>(),
        [1, 7, 8]
    );
    assert_eq!(warm[6], mol.peoe_vsa_7().unwrap());
    let copy = mol.clone();
    copy.estate_vsa_6().unwrap();
    assert_eq!(copy.peoe_vsa().unwrap(), cold);
    assert_eq!(mol.peoe_vsa().unwrap(), warm);
}
