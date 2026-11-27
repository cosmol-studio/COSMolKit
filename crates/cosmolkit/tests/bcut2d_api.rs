#![cfg(feature = "descriptors")]

use cosmolkit::{
    BondId, BondOrder, DescriptorError, DescriptorReadError, GasteigerChargeError, Molecule,
    SdfReadParams, SmilesParseParams,
};

fn raw(smiles: &str) -> Molecule {
    Molecule::from_smiles_with_params(
        smiles,
        &SmilesParseParams {
            sanitize: false,
            remove_hs: false,
            ..Default::default()
        },
    )
    .unwrap()
}

fn assert_slots(actual: &[f64], reference: &[f64; 8]) {
    assert_eq!(actual.len(), reference.len());
    for (slot, (&a, &r)) in actual.iter().zip(reference).enumerate() {
        if r.is_nan() {
            assert_eq!(a.abs().to_bits(), r.abs().to_bits(), "slot {slot}");
            continue;
        }
        assert!(
            (a - r).abs() <= 1e-12 * r.abs().max(1.0),
            "slot {slot}: {a} vs {r}"
        );
    }
}

#[test]
fn mixed_ring_and_long_path_match_all_eight_pinned_eigenvalues() {
    // RDKit 4ce8a063, BCUT2D, Eigen 3.4.0. Compare each weighting and both extrema.
    for (smiles, reference) in [
        (
            "CCO",
            [
                16.248991329327914,
                10.908168035086371,
                1.3402935186707783,
                -1.52235869276019,
                1.2962426382141046,
                -1.5683635470355741,
                3.8358023355894861,
                0.31840561317391147,
            ],
        ),
        (
            "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
            [
                16.161282315860802,
                10.488187986922153,
                2.1102552164294219,
                -2.0215174970665886,
                2.0118198978857045,
                -2.3696654059207396,
                5.6916045887805531,
                0.7054356050886158,
            ],
        ),
        (
            "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
            [
                14.026648237590239,
                10.021047098698981,
                1.9622819293660956,
                -2.0432875999572162,
                2.1597482375902306,
                -1.8458529013010261,
                4.5186482375902344,
                0.51304709869897724,
            ],
        ),
    ] {
        assert_slots(
            &Molecule::from_smiles(smiles).unwrap().bcut2d().unwrap(),
            &reference,
        );
    }
}

#[test]
fn undefined_charge_eigenvalues_leave_the_other_three_pairs_intact() {
    for (smiles, reference) in [
        (
            "[C+9](F)(F)(F)(F)(F)(F)(F)(F)(F)(F)(F)(F)F",
            [
                20.535089560352173,
                10.485910439647832,
                f64::NAN,
                f64::NAN,
                3.8665627872331734,
                -3.3530727872331729,
                5.94006119412735,
                -1.5770611941273487,
            ],
        ),
        (
            "[C+8](F)(F)(F)(F)(F)(F)(F)(F)F",
            [
                20.116374374372636,
                10.900625625627365,
                f64::NAN,
                f64::NAN,
                3.259755255727092,
                -2.7502652557270917,
                5.3624282508407255,
                -1.0034282508407253,
            ],
        ),
    ] {
        assert_slots(
            &Molecule::from_smiles(smiles).unwrap().bcut2d().unwrap(),
            &reference,
        );
    }
}

#[test]
fn one_atom_eigenvalues_are_the_diagonal_even_for_infinite_charge() {
    for smiles in ["[CH8+6]", "[2H][CH7+6]"] {
        let mol = Molecule::from_smiles(smiles).unwrap();
        assert_eq!(
            mol.bcut2d().unwrap(),
            [
                12.011,
                12.011,
                f64::NEG_INFINITY,
                f64::NEG_INFINITY,
                0.08129,
                0.08129,
                3.243,
                3.243,
            ],
        );
    }
}

#[test]
fn isotope_mass_and_removed_hydrogens_have_exact_identities() {
    for (smiles, mass) in [("C", 12.011), ("[13CH4]", 13.00335484), ("[99CH4]", 99.0)] {
        let values = Molecule::from_smiles(smiles).unwrap().bcut2d().unwrap();
        assert_eq!(
            values,
            [
                mass,
                mass,
                -0.07755789119307553,
                -0.07755789119307553,
                0.1441,
                0.1441,
                2.503,
                2.503
            ]
        );
    }
    for mol in [Molecule::new(), raw("[H][H]")] {
        assert!(mol.bcut2d().unwrap().iter().map(|x| x.to_bits()).eq([0; 8]));
    }
    let heavy_water = raw("[2H]O[2H]");
    assert_eq!(
        heavy_water.bcut2d().unwrap(),
        Molecule::from_smiles("O").unwrap().bcut2d().unwrap()
    );
    assert_eq!(heavy_water.num_atoms(), 3);
}

#[test]
fn private_sanitization_leaves_the_input_and_copy_unprepared() {
    for (smiles, equivalent) in [("C1=CC=CC=C1", "c1ccccc1"), ("CCO", "CCO")] {
        let molecule = raw(smiles);
        let copy = molecule.copy();
        let before = molecule.topology().clone();
        let prepared = Molecule::from_smiles(equivalent).unwrap();
        assert_eq!(molecule.bcut2d().unwrap(), prepared.bcut2d().unwrap());
        for unchanged in [&molecule, &copy] {
            assert_eq!(unchanged.topology(), &before);
            assert!(unchanged.estate_indices().is_err());
        }
    }
}

#[test]
fn substance_group_keeps_its_only_hydrogen() {
    let molblock = "hydroxyl\n  COSMolKit\n\n  2  1  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\nM  STY  1   1 SUP\nM  SAL   1  1   1\nM  SMT   1 Hydrogen\nM  END\n";
    let molecule = Molecule::from_sdf_with_params(
        molblock,
        &SdfReadParams {
            remove_hs: false,
            ..Default::default()
        },
    )
    .unwrap();
    let before = molecule.topology().clone();
    let values = molecule.bcut2d().unwrap();
    let middle = (15.999_f64 + 1.008) / 2.;
    let spread = ((15.999_f64 - 1.008) / 2.).hypot(1.0);
    assert!((values[0] - middle - spread).abs() < 1e-12);
    assert!((values[1] - middle + spread).abs() < 1e-12);
    assert_eq!(molecule.topology(), &before);
}

#[test]
fn failures_follow_chemistry_then_bond_order_without_partial_results() {
    let overvalent = raw("C(F)(F)(F)(F)F");
    let before = overvalent.topology().clone();
    assert!(matches!(
        overvalent.bcut2d(),
        Err(DescriptorReadError::Algorithm {
            source: DescriptorError::Hydrogens { .. }
        })
    ));
    assert_eq!(overvalent.topology(), &before);
    for (smiles, source_index) in [("N->B", 0), ("[2H]C([2H])N->B", 3)] {
        let dative = Molecule::from_smiles(smiles).unwrap();
        assert!(
            matches!(dative.bcut2d(), Err(DescriptorReadError::Algorithm {source: DescriptorError::UnsupportedBcutBond {bond, order: BondOrder::Dative}}) if bond == BondId::new(source_index))
        );
    }
    for smiles in ["[Na+]", "[Na+].N->B"] {
        let mol = Molecule::from_smiles(smiles).unwrap();
        let saved = mol.gasteiger_charges().unwrap();
        let copy = mol.copy();
        for input in [&mol, &copy] {
            assert!(
                matches!(input.bcut2d(), Err(DescriptorReadError::Algorithm {source: DescriptorError::GasteigerCharge {source: GasteigerChargeError::MissingParameters {element, ..}}}) if element == "Na")
            );
            assert_eq!(
                input.gasteiger_charges().unwrap()[0]
                    .atom_charge()
                    .to_bits(),
                saved[0].atom_charge().to_bits()
            );
        }
    }
}
