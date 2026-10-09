//! Public ring-query API regressions for the live ring-state slice.
#![cfg(all(feature = "cap-descriptors", feature = "cap-smiles"))]

use cosmolkit::Molecule;

/// Q1: 72 public literal-table calls (18 frozen cases x two remove-H
/// policies x num_rings/num_heterocycles) through the REAL public
/// constructor with sanitize=true, proving the constructor prepared the
/// rows the queries read.
#[test]
fn ring_live_public_q1_literal_table_calls() {
    let cases: [(&str, usize, usize, [u32; 2]); 18] = [
        ("", 0, 0, [0, 0]),
        ("C", 1, 1, [0, 0]),
        ("C1CCCCC1", 6, 6, [1, 0]),
        ("C1=CCCCC1", 6, 6, [1, 0]),
        ("c1ccccc1", 6, 6, [1, 0]),
        ("n1ccccc1", 6, 6, [1, 1]),
        ("N1CCCCC1", 6, 6, [1, 1]),
        ("C1CCC2(CC1)CCCC2", 10, 10, [2, 0]),
        ("c1ccc2ccccc2c1", 10, 10, [2, 0]),
        ("C12C3C4C1C5C2C3C45", 8, 8, [6, 0]),
        ("[H]C1CCCCC1", 7, 6, [1, 0]),
        ("[2H]C1CCCCC1", 7, 7, [1, 0]),
        ("*1CCCC1", 5, 5, [1, 1]),
        ("C1CCC2CCCCC2C1", 10, 10, [2, 0]),
        ("C1CC2CCC1C2", 7, 7, [2, 0]),
        ("c1ccc2CCCCc2c1", 10, 10, [2, 0]),
        ("c1ccncc1C1CCCCC1", 12, 12, [2, 1]),
        ("N12N3N4N1N5N2N3N45", 8, 8, [6, 6]),
    ];
    let mut calls = 0usize;
    for (smiles, keep_atoms, remove_atoms, literals) in cases {
        for remove_hs in [false, true] {
            let label = format!("{smiles}/rh={remove_hs}");
            let molecule = Molecule::from_smiles_with_params(
                smiles,
                &cosmolkit_smiles::SmilesParseParams {
                    sanitize: true,
                    remove_hs,
                    ..cosmolkit_smiles::SmilesParseParams::default()
                },
            )
            .unwrap();
            // The constructor really prepared the expected topology.
            let expected = if remove_hs { remove_atoms } else { keep_atoms };
            assert_eq!(molecule.num_atoms(), expected, "{label}: atoms");
            calls += 1;
            assert_eq!(molecule.num_rings().unwrap(), literals[0], "{label}");
            calls += 1;
            assert_eq!(molecule.num_heterocycles().unwrap(), literals[1], "{label}");
        }
    }
    assert_eq!(calls, 72, "exact census");
}

/// Q2: 108 unchanged public literal-table calls (18 frozen cases x two
/// remove-H policies x num_aromatic_rings/num_saturated_rings/
/// num_aliphatic_rings) through the REAL public constructor with
/// sanitize=true, checking columns 2/3/4 of the frozen table.
#[test]
fn ring_live_public_q2_literal_table_calls() {
    let cases: [(&str, usize, usize, [u32; 3]); 18] = [
        ("", 0, 0, [0, 0, 0]),
        ("C", 1, 1, [0, 0, 0]),
        ("C1CCCCC1", 6, 6, [0, 1, 1]),
        ("C1=CCCCC1", 6, 6, [0, 0, 1]),
        ("c1ccccc1", 6, 6, [1, 0, 0]),
        ("n1ccccc1", 6, 6, [1, 0, 0]),
        ("N1CCCCC1", 6, 6, [0, 1, 1]),
        ("C1CCC2(CC1)CCCC2", 10, 10, [0, 2, 2]),
        ("c1ccc2ccccc2c1", 10, 10, [2, 0, 0]),
        ("C12C3C4C1C5C2C3C45", 8, 8, [0, 6, 6]),
        ("[H]C1CCCCC1", 7, 6, [0, 1, 1]),
        ("[2H]C1CCCCC1", 7, 7, [0, 1, 1]),
        ("*1CCCC1", 5, 5, [0, 1, 1]),
        ("C1CCC2CCCCC2C1", 10, 10, [0, 2, 2]),
        ("C1CC2CCC1C2", 7, 7, [0, 2, 2]),
        ("c1ccc2CCCCc2c1", 10, 10, [1, 0, 1]),
        ("c1ccncc1C1CCCCC1", 12, 12, [1, 1, 1]),
        ("N12N3N4N1N5N2N3N45", 8, 8, [0, 6, 6]),
    ];
    let mut calls = 0usize;
    for (smiles, keep_atoms, remove_atoms, literals) in cases {
        for remove_hs in [false, true] {
            let label = format!("{smiles}/rh={remove_hs}");
            let molecule = Molecule::from_smiles_with_params(
                smiles,
                &cosmolkit_smiles::SmilesParseParams {
                    sanitize: true,
                    remove_hs,
                    ..cosmolkit_smiles::SmilesParseParams::default()
                },
            )
            .unwrap();
            let expected = if remove_hs { remove_atoms } else { keep_atoms };
            assert_eq!(molecule.num_atoms(), expected, "{label}: atoms");
            calls += 1;
            assert_eq!(
                molecule.num_aromatic_rings().unwrap(),
                literals[0],
                "{label}: aromatic"
            );
            calls += 1;
            assert_eq!(
                molecule.num_saturated_rings().unwrap(),
                literals[1],
                "{label}: saturated"
            );
            calls += 1;
            assert_eq!(
                molecule.num_aliphatic_rings().unwrap(),
                literals[2],
                "{label}: aliphatic"
            );
        }
    }
    assert_eq!(calls, 108, "exact census");
}

/// Q3: 216 unchanged public literal-table calls (18 frozen cases x two
/// remove-H policies x the six combined heterocycle/carbocycle queries),
/// checking columns 5-10 of the frozen table.
#[test]
fn ring_live_public_q3_literal_table_calls() {
    let cases: [(&str, usize, usize, [u32; 6]); 18] = [
        ("", 0, 0, [0, 0, 0, 0, 0, 0]),
        ("C", 1, 1, [0, 0, 0, 0, 0, 0]),
        ("C1CCCCC1", 6, 6, [0, 0, 0, 1, 0, 1]),
        ("C1=CCCCC1", 6, 6, [0, 0, 0, 1, 0, 0]),
        ("c1ccccc1", 6, 6, [0, 1, 0, 0, 0, 0]),
        ("n1ccccc1", 6, 6, [1, 0, 0, 0, 0, 0]),
        ("N1CCCCC1", 6, 6, [0, 0, 1, 0, 1, 0]),
        ("C1CCC2(CC1)CCCC2", 10, 10, [0, 0, 0, 2, 0, 2]),
        ("c1ccc2ccccc2c1", 10, 10, [0, 2, 0, 0, 0, 0]),
        ("C12C3C4C1C5C2C3C45", 8, 8, [0, 0, 0, 6, 0, 6]),
        ("[H]C1CCCCC1", 7, 6, [0, 0, 0, 1, 0, 1]),
        ("[2H]C1CCCCC1", 7, 7, [0, 0, 0, 1, 0, 1]),
        ("*1CCCC1", 5, 5, [0, 0, 1, 0, 1, 0]),
        ("C1CCC2CCCCC2C1", 10, 10, [0, 0, 0, 2, 0, 2]),
        ("C1CC2CCC1C2", 7, 7, [0, 0, 0, 2, 0, 2]),
        ("c1ccc2CCCCc2c1", 10, 10, [0, 1, 0, 1, 0, 0]),
        ("c1ccncc1C1CCCCC1", 12, 12, [1, 0, 0, 1, 0, 1]),
        ("N12N3N4N1N5N2N3N45", 8, 8, [0, 0, 6, 0, 6, 0]),
    ];
    let mut calls = 0usize;
    for (smiles, keep_atoms, remove_atoms, literals) in cases {
        for remove_hs in [false, true] {
            let label = format!("{smiles}/rh={remove_hs}");
            let molecule = Molecule::from_smiles_with_params(
                smiles,
                &cosmolkit_smiles::SmilesParseParams {
                    sanitize: true,
                    remove_hs,
                    ..cosmolkit_smiles::SmilesParseParams::default()
                },
            )
            .unwrap();
            let expected = if remove_hs { remove_atoms } else { keep_atoms };
            assert_eq!(molecule.num_atoms(), expected, "{label}: atoms");
            calls += 1;
            assert_eq!(
                molecule.num_aromatic_heterocycles().unwrap(),
                literals[0],
                "{label}: aromatic hetero"
            );
            calls += 1;
            assert_eq!(
                molecule.num_aromatic_carbocycles().unwrap(),
                literals[1],
                "{label}: aromatic carbo"
            );
            calls += 1;
            assert_eq!(
                molecule.num_aliphatic_heterocycles().unwrap(),
                literals[2],
                "{label}: aliphatic hetero"
            );
            calls += 1;
            assert_eq!(
                molecule.num_aliphatic_carbocycles().unwrap(),
                literals[3],
                "{label}: aliphatic carbo"
            );
            calls += 1;
            assert_eq!(
                molecule.num_saturated_heterocycles().unwrap(),
                literals[4],
                "{label}: saturated hetero"
            );
            calls += 1;
            assert_eq!(
                molecule.num_saturated_carbocycles().unwrap(),
                literals[5],
                "{label}: saturated carbo"
            );
        }
    }
    assert_eq!(calls, 216, "exact census");
}

/// F1 runtime parity: under any feature set that enables both cap-smiles
/// and cap-descriptors (including the reduced selection WITHOUT cap-rings),
/// the constructor prepares the same ring counts the full build produces.
#[test]
fn ring_live_public_feature_parity_counts() {
    let benzene = Molecule::from_smiles("c1ccccc1").unwrap();
    assert_eq!(benzene.num_rings().unwrap(), 1);
    assert_eq!(benzene.num_aromatic_rings().unwrap(), 1);
    assert_eq!(benzene.num_aromatic_carbocycles().unwrap(), 1);
    assert_eq!(benzene.num_heterocycles().unwrap(), 0);
    let pyridine = Molecule::from_smiles("n1ccccc1").unwrap();
    assert_eq!(pyridine.num_rings().unwrap(), 1);
    assert_eq!(pyridine.num_aromatic_heterocycles().unwrap(), 1);
    let raw = Molecule::from_smiles_with_params(
        "c1ccccc1",
        &cosmolkit_smiles::SmilesParseParams {
            sanitize: false,
            remove_hs: false,
            ..cosmolkit_smiles::SmilesParseParams::default()
        },
    )
    .unwrap();
    assert!(
        matches!(
            raw.num_rings(),
            Err(cosmolkit::DescriptorReadError::MissingInitializedRings)
        ),
        "builder-only absence keeps the typed error"
    );
    assert_eq!(raw.num_heterocycles().unwrap(), 0);
}
