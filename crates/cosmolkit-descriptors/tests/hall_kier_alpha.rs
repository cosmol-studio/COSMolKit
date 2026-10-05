//! Fixed owner-boundary regressions; no reference execution or generated data.
//! RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8, BSD license:
//! ConnectivityDescriptors.cpp:71-166,267-295; atomic_data.cpp rB0 column;
//! Descriptors/test.cpp:1332-1362 for the ten molecule examples.

use cosmolkit_descriptors::{DescriptorError, HALL_KIER_ALPHA_VERSION, hall_kier_alpha};
use cosmolkit_model::{Atom, AtomId, AtomSpec, Element, Hybridization};

const HYBRIDIZATIONS: [Hybridization; 9] = [
    Hybridization::Unspecified,
    Hybridization::S,
    Hybridization::Sp,
    Hybridization::Sp2,
    Hybridization::Sp3,
    Hybridization::Sp2d,
    Hybridization::Sp3d,
    Hybridization::Sp3d2,
    Hybridization::Other,
];

fn atom(id: usize, symbol: &str, hybridization: Hybridization) -> Atom {
    Atom::from_spec(
        AtomId::new(id),
        AtomSpec::new(Element::from_symbol(symbol).unwrap()).with_hybridization(hybridization),
    )
}

fn bits(rows: &[f64]) -> Vec<u64> {
    rows.iter().map(|value| value.to_bits()).collect()
}

#[test]
fn hall_kier_special_elements_use_every_stored_hybridization() {
    assert_eq!(HALL_KIER_ALPHA_VERSION, "1.2.0");
    // Fixed literal columns: UNSPECIFIED, S, SP, SP2, SP3, SP2D, SP3D,
    // SP3D2, OTHER. Non-SP/SP2 classifications take the source default.
    let cases: [(&str, [f64; 9]); 10] = [
        ("H", [0.0; 9]),
        ("C", [0.0, 0.0, -0.22, -0.13, 0.0, 0.0, 0.0, 0.0, 0.0]),
        (
            "N",
            [
                -0.04, -0.04, -0.29, -0.20, -0.04, -0.04, -0.04, -0.04, -0.04,
            ],
        ),
        (
            "O",
            [
                -0.04, -0.04, -0.04, -0.20, -0.04, -0.04, -0.04, -0.04, -0.04,
            ],
        ),
        ("F", [-0.07; 9]),
        ("P", [0.43, 0.43, 0.43, 0.30, 0.43, 0.43, 0.43, 0.43, 0.43]),
        ("S", [0.35, 0.35, 0.35, 0.22, 0.35, 0.35, 0.35, 0.35, 0.35]),
        ("Cl", [0.29; 9]),
        ("Br", [0.48; 9]),
        ("I", [0.73; 9]),
    ];
    let mut calls = 0;
    for (symbol, expected) in cases {
        for (hybridization, value) in HYBRIDIZATIONS.into_iter().zip(expected) {
            let atoms = [atom(900, symbol, hybridization)];
            let before = atoms.clone();
            for _ in 0..2 {
                assert_eq!(
                    hall_kier_alpha(&atoms, None).unwrap().to_bits(),
                    value.to_bits()
                );
                let mut sink = [1234.5];
                let result = hall_kier_alpha(&atoms, Some(&mut sink)).unwrap();
                assert_eq!(
                    result.to_bits(),
                    value.to_bits(),
                    "{symbol} {hybridization:?}"
                );
                assert_eq!(sink[0].to_bits(), value.to_bits());
                assert_eq!(atoms, before);
                calls += 2;
            }
        }
    }
    assert_eq!(calls, 360);
}

#[test]
fn hall_kier_fallback_uses_pinned_bond_radii() {
    // Independent pinned rB0 literals, never read from the tested owner or
    // the core table to obtain expectations. 0.77 is source carbon rB0,
    // distinct from its covalent radius (0.76).
    let cases: [(&str, f64); 15] = [
        ("He", 0.7 / 0.77 - 1.0),
        ("Li", 1.23 / 0.77 - 1.0),
        ("B", 0.82 / 0.77 - 1.0),
        ("Ne", 0.7 / 0.77 - 1.0),
        ("Na", 1.54 / 0.77 - 1.0),
        ("Mg", 1.36 / 0.77 - 1.0),
        ("Al", 1.18 / 0.77 - 1.0),
        ("Si", 0.937 / 0.77 - 1.0),
        ("Ar", 1.74 / 0.77 - 1.0),
        ("Fe", 1.17 / 0.77 - 1.0),
        ("Se", 1.17 / 0.77 - 1.0),
        ("Xe", 1.98 / 0.77 - 1.0),
        ("U", 1.58 / 0.77 - 1.0),
        ("Cm", -1.0),
        ("Og", -1.0),
    ];
    let mut calls = 0;
    for (symbol, expected) in cases {
        for hybridization in HYBRIDIZATIONS {
            let atoms = [atom(8, symbol, hybridization)];
            let before = atoms.clone();
            assert_eq!(
                hall_kier_alpha(&atoms, None).unwrap().to_bits(),
                expected.to_bits()
            );
            let mut sink = [f64::NAN, -9876.5];
            let result = hall_kier_alpha(&atoms, Some(&mut sink)).unwrap();
            assert_eq!(
                result.to_bits(),
                expected.to_bits(),
                "{symbol} {hybridization:?}"
            );
            assert_eq!(sink[0].to_bits(), expected.to_bits());
            assert_eq!(sink[1].to_bits(), (-9876.5_f64).to_bits());
            assert_eq!(atoms, before);
            calls += 2;
        }
    }
    assert_eq!(calls, 270);
}

#[test]
fn hall_kier_wildcards_and_oversized_tail_preserve_bits() {
    let atoms = [
        atom(33, "*", Hybridization::Sp2),
        atom(44, "C", Hybridization::Sp),
        atom(55, "*", Hybridization::Other),
        atom(66, "I", Hybridization::Unspecified),
        atom(77, "*", Hybridization::Sp3),
    ];
    let before = atoms.clone();
    let nan = f64::from_bits(0x7ff8_0000_0000_0042);
    let mut sink = [nan, 88.0, f64::NEG_INFINITY, 99.0, -0.0, f64::INFINITY, nan];
    let expected = [
        nan,
        -0.22,
        f64::NEG_INFINITY,
        0.73,
        -0.0,
        f64::INFINITY,
        nan,
    ];
    for _ in 0..2 {
        let scalar = hall_kier_alpha(&atoms, Some(&mut sink)).unwrap();
        assert_eq!(scalar.to_bits(), 0.51_f64.to_bits());
        assert_eq!(bits(&sink), bits(&expected));
        assert_eq!(
            hall_kier_alpha(&atoms, None).unwrap().to_bits(),
            scalar.to_bits()
        );
        assert_eq!(atoms, before);
    }
}

#[test]
fn hall_kier_short_buffer_errors_before_writes() {
    let atoms = [
        atom(0, "Cl", Hybridization::Sp3),
        atom(1, "C", Hybridization::Sp),
        atom(2, "*", Hybridization::Unspecified),
    ];
    let before = atoms.clone();
    for len in 0..atoms.len() {
        let mut sink = vec![f64::from_bits(0x7ff8_0000_0000_0033); len];
        let sink_before = bits(&sink);
        let error = hall_kier_alpha(&atoms, Some(&mut sink)).unwrap_err();
        assert_eq!(
            error,
            DescriptorError::InvalidHallKierContributionRows {
                actual: len,
                minimum: 3
            }
        );
        assert!(error.to_string().contains("hall_kier_alpha"));
        assert_eq!(bits(&sink), sink_before);
        assert_eq!(atoms, before);
    }
    // Wildcards count toward the source size precondition despite skipping
    // all writes; an insufficient all-wildcard buffer still errors.
    let wildcards = [
        atom(0, "*", Hybridization::Unspecified),
        atom(1, "*", Hybridization::Sp),
    ];
    let mut sink = [9.0];
    assert_eq!(
        hall_kier_alpha(&wildcards, Some(&mut sink)),
        Err(DescriptorError::InvalidHallKierContributionRows {
            actual: 1,
            minimum: 2
        })
    );
    assert_eq!(sink, [9.0]);
}

#[test]
fn hall_kier_empty_and_all_wildcard_inputs_return_positive_zero() {
    let mut empty = [];
    let mut sink = [-0.0, f64::from_bits(0x7ff8_0000_0000_0055)];
    let before = bits(&sink);
    assert_eq!(hall_kier_alpha(&[], None).unwrap().to_bits(), 0);
    assert_eq!(hall_kier_alpha(&[], Some(&mut empty)).unwrap().to_bits(), 0);
    assert_eq!(hall_kier_alpha(&[], Some(&mut sink)).unwrap().to_bits(), 0);
    let atoms = [
        atom(700, "*", Hybridization::Other),
        atom(600, "*", Hybridization::Sp2),
    ];
    assert_eq!(
        hall_kier_alpha(&atoms, Some(&mut sink)).unwrap().to_bits(),
        0
    );
    assert_eq!(bits(&sink), before);
}

#[test]
fn hall_kier_contributions_follow_rows_and_sum_in_source_order() {
    // Deliberately unordered and repeated detached IDs must not index sinks.
    let atoms = [
        atom(usize::MAX, "S", Hybridization::Sp2),
        atom(17, "N", Hybridization::Sp),
        atom(17, "Cl", Hybridization::Other),
        atom(2, "H", Hybridization::Sp2),
    ];
    let mut sink = [7.0; 4];
    assert_eq!(
        hall_kier_alpha(&atoms, Some(&mut sink)).unwrap().to_bits(),
        0.22_f64.to_bits()
    );
    assert_eq!(bits(&sink), bits(&[0.22, -0.29, 0.29, 0.0]));
    // This permutation has the same mathematical sum but different binary64
    // rounding. A sorted or reassociated reduction changes source behavior.
    let reordered = [
        atoms[2].clone(),
        atoms[0].clone(),
        atoms[1].clone(),
        atoms[3].clone(),
    ];
    assert_eq!(
        hall_kier_alpha(&reordered, None).unwrap().to_bits(),
        0.22000000000000003_f64.to_bits()
    );
}

#[test]
fn hall_kier_uses_stored_hybridization_and_ignores_other_atom_state() {
    let neutral = atom(0, "C", Hybridization::Unspecified);
    let decorated = Atom::from_spec(
        AtomId::new(1),
        AtomSpec::new(Element::C)
            .with_isotope(13)
            .with_formal_charge(1)
            .with_explicit_hydrogens(3)
            .with_aromatic(true)
            .with_hybridization(Hybridization::Unspecified),
    );
    let atoms = [neutral, decorated];
    let before = atoms.clone();
    let mut sink = [999.0; 2];
    assert_eq!(
        hall_kier_alpha(&atoms, Some(&mut sink)).unwrap().to_bits(),
        0
    );
    assert_eq!(bits(&sink), vec![0, 0]);
    assert_eq!(atoms, before);
    let mut changed = atoms.clone();
    changed[1].set_hybridization(Hybridization::Sp2);
    assert_eq!(
        hall_kier_alpha(&changed, None).unwrap().to_bits(),
        (-0.13_f64).to_bits()
    );
    assert_eq!(atoms, before);
}

#[test]
fn hall_kier_upstream_molecule_examples_with_real_preparation() {
    let cases: [(&str, f64); 10] = [
        ("C=O", -0.33),
        ("CCC1(CC)C(=O)NC(=O)N(C)C1=O", -1.39),
        ("OCC(O)C(O)C(O)C(O)CO", -0.24),
        ("OCC1OC(O)C(O)C(O)C1O", -0.24),
        ("Fc1c[nH]c(=O)[nH]c1=O", -1.39),
        ("OC1CNC(C(=O)O)C1", -0.61),
        ("CCCc1[nH]c(=S)[nH]c(=O)c1", -0.90),
        ("CN(CCCl)CCCl", 0.54),
        ("CBr", 0.48),
        ("CI", 0.73),
    ];
    for (smiles, expected) in cases {
        // Only fixture setup performs real chemistry. The counted owner
        // calls below use the unchanged, prepared borrowed atom rows.
        let record = cosmolkit_smiles::parse_smiles(smiles, &Default::default()).unwrap();
        let prepared =
            cosmolkit_core::sanitize_topology(&record.topology, &Default::default()).unwrap();
        let before = prepared.topology.clone();
        let scalar = hall_kier_alpha(&prepared.topology.atoms, None).unwrap();
        let mut sink = vec![f64::NAN; prepared.topology.atoms.len()];
        let with_sink = hall_kier_alpha(&prepared.topology.atoms, Some(&mut sink)).unwrap();
        assert!(
            (scalar - expected).abs() < 1e-12,
            "{smiles}: {scalar} vs {expected}"
        );
        assert_eq!(scalar.to_bits(), with_sink.to_bits());
        assert_eq!(scalar.to_bits(), sink.iter().sum::<f64>().to_bits());
        assert_eq!(prepared.topology, before);
    }
}
