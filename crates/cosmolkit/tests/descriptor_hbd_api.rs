//! HBD-PUBLIC1-50: general prepared `num_hbd` public regressions.
#![cfg(all(feature = "cap-descriptors", feature = "cap-smiles"))]

use cosmolkit::{DescriptorReadError, Molecule, SmilesParseParams};

/// Frozen 13-case (SMILES, expected) literal table for the general
/// SMARTS-based donor count, independently recorded by ROOT from RDKit
/// 2026.03.1 under BOTH removeHs policies BEFORE any CK run. Isolated S
/// has two H and is 0; ammonium is 1; never derived from the SUT.
const HBD_CASES: [(&str, u32); 13] = [
    ("", 0),
    ("CCO", 1),
    ("NCC(=O)O", 2),
    ("[NH4+]", 1),
    ("c1cc[nH]c1", 1),
    ("c1ccccc1", 0),
    ("CC#CC", 0),
    ("[H]N([H])[H]", 1),
    ("CS", 1),
    ("N", 1),
    ("[OH2]", 0),
    ("S", 0),
    ("COC", 0),
];

fn sanitized(smiles: &str, remove_hs: bool) -> Molecule {
    let params = SmilesParseParams {
        sanitize: true,
        remove_hs,
        ..Default::default()
    };
    Molecule::from_smiles_with_params(smiles, &params).unwrap()
}

/// The frozen 104-call product: 13 cases x 2 remove-H policies x 2
/// receivers (original + Arc-sharing peer clone) x 2 repeats. The census
/// increments only after an ACTUAL invocation; real constructor
/// prerequisites are asserted BEFORE every call and row/bond/property
/// counts stay unchanged after it.
#[test]
fn descriptor_hbd_public_literal_product() {
    let mut calls = 0usize;
    for (smiles, expected) in HBD_CASES {
        for remove_hs in [false, true] {
            let label = format!("{smiles}/rh={remove_hs}");
            let original = sanitized(smiles, remove_hs);
            let peer = original.clone();
            let rows = original.num_atoms();
            let bonds = original.num_bonds();
            let properties = original.properties().clone();
            if smiles.is_empty() {
                assert_eq!(rows, 0, "{label}: empty molecule is a real 0-row fact");
                assert_eq!(bonds, 0, "{label}: empty molecule has no bonds");
            }
            for (who, receiver) in [("original", &original), ("peer", &peer)] {
                for repeat in 0..2 {
                    let call_label = format!("{label}/{who}#{repeat}");
                    assert_eq!(
                        receiver.num_atoms(),
                        rows,
                        "{call_label}: prerequisite rows before the call"
                    );
                    assert_eq!(
                        receiver.num_bonds(),
                        bonds,
                        "{call_label}: prerequisite bonds before the call"
                    );
                    let count = receiver
                        .num_hbd()
                        .unwrap_or_else(|error| panic!("{call_label}: {error:?}"));
                    calls += 1;
                    assert_eq!(count, expected, "{call_label}: literal output");
                    assert_eq!(receiver.num_atoms(), rows, "{call_label}: rows unchanged");
                    assert_eq!(receiver.num_bonds(), bonds, "{call_label}: bonds unchanged");
                    assert_eq!(
                        receiver.properties(),
                        &properties,
                        "{call_label}: properties unchanged"
                    );
                }
            }
        }
    }
    assert_eq!(calls, 104, "exact 104-call census");
}

/// Separately counted source distinctions: the general atom-count pattern
/// and the direct Lipinski donor-hydrogen sum disagree on the SAME real
/// molecules by frozen source definition (never fitted).
#[test]
fn descriptor_hbd_public_supplementary_general_vs_lipinski() {
    let mut calls = 0usize;
    for remove_hs in [false, true] {
        let label = format!("rh={remove_hs}");
        let thioether = sanitized("CS", remove_hs);
        assert_eq!(thioether.num_hbd().unwrap(), 1, "CS general {label}");
        calls += 1;
        assert_eq!(thioether.lipinski_hbd().unwrap(), 0, "CS direct {label}");
        calls += 1;
        let ammonium = sanitized("[NH4+]", remove_hs);
        assert_eq!(ammonium.num_hbd().unwrap(), 1, "ammonium general {label}");
        calls += 1;
        assert_eq!(
            ammonium.lipinski_hbd().unwrap(),
            4,
            "ammonium direct {label}"
        );
        calls += 1;
    }
    assert_eq!(calls, 8, "exact 8-call supplementary census");
}

/// Real raw missing-state regressions on PUBLIC calls: unsanitized
/// molecules report the typed MissingPreparedValence with no fabricated
/// cause, and the whole input survives the failed call. This is the
/// documented prepared-state boundary limitation, not a chemical
/// divergence.
#[test]
fn descriptor_hbd_public_missing_prepared_valence() {
    let mut calls = 0usize;
    for (label, smiles) in [("raw-cco", "CCO"), ("pentavalent", "C(C)(C)(C)(C)C")] {
        let params = SmilesParseParams {
            sanitize: false,
            remove_hs: false,
            ..Default::default()
        };
        let molecule = Molecule::from_smiles_with_params(smiles, &params).unwrap();
        let rows = molecule.num_atoms();
        let bonds = molecule.num_bonds();
        let properties = molecule.properties().clone();
        let error = molecule.num_hbd().unwrap_err();
        calls += 1;
        assert!(
            matches!(error, DescriptorReadError::MissingPreparedValence),
            "{label}: got {error:?}"
        );
        assert!(
            std::error::Error::source(&error).is_none(),
            "{label}: no fabricated cause"
        );
        assert_eq!(molecule.num_atoms(), rows, "{label}: rows preserved");
        assert_eq!(molecule.num_bonds(), bonds, "{label}: bonds preserved");
        assert_eq!(
            molecule.properties(),
            &properties,
            "{label}: properties preserved"
        );
    }
    assert_eq!(calls, 2, "exact raw missing-state census");
}
