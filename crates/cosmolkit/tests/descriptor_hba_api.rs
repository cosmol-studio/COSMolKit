//! HBA-PUBLIC1-42: general prepared `num_hba` public regressions.
#![cfg(all(feature = "cap-descriptors", feature = "cap-smiles"))]

use cosmolkit::{DescriptorReadError, Molecule, SmilesParseParams};

/// Frozen 15-case (SMILES, expected) literal table for the general
/// SMARTS-based acceptor count, fixed BEFORE implementation from the
/// pinned RDKit 351f8f3 behavior (including the test.cpp:2227-2234
/// #8997 discriminators `c1cccn1C` -> 0 and `c1cccc(=O)n1C` -> 1).
/// The same expected value holds under BOTH constructor policies
/// (sanitize=true, remove_hydrogens=false and true); counts never
/// depend on this frozen table being derived from the SUT.
const HBA_CASES: [(&str, u32); 15] = [
    ("", 0),
    ("CCO", 1),
    ("COC", 1),
    ("[OH2]", 0),
    ("CC(=O)O", 1),
    ("CC(=O)[O-]", 2),
    ("CC(=O)N", 1),
    ("CC#N", 1),
    ("CN", 1),
    ("Nc1ccccc1", 1),
    ("c1ccncc1", 1),
    ("c1ccoc1", 1),
    ("c1cccc(=O)n1C", 1),
    ("c1cccn1C", 0),
    ("CC(=O)OCC", 2),
];

fn sanitized(smiles: &str, remove_hydrogens: bool) -> Molecule {
    let params = SmilesParseParams {
        sanitize: true,
        remove_hydrogens,
        ..Default::default()
    };
    Molecule::from_smiles_with_params(smiles, &params).unwrap()
}

/// The frozen 120-call product: 15 cases x 2 remove-H policies x 2
/// receivers (original + Arc-sharing peer clone) x 2 repeats. The census
/// increments only after an ACTUAL invocation; the constructed input
/// identity (row/bond/property facts) is asserted BEFORE every call and
/// preserved after it.
#[test]
fn descriptor_hba_public_literal_product() {
    let mut calls = 0usize;
    for (smiles, expected) in HBA_CASES {
        for remove_hydrogens in [false, true] {
            let label = format!("{smiles}/rh={remove_hydrogens}");
            let original = sanitized(smiles, remove_hydrogens);
            let peer = original.clone();
            // Actual fixture facts captured from the real constructor.
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
                        "{call_label}: fixture rows before the call"
                    );
                    assert_eq!(
                        receiver.num_bonds(),
                        bonds,
                        "{call_label}: fixture bonds before the call"
                    );
                    let count = receiver
                        .num_hba()
                        .unwrap_or_else(|error| panic!("{call_label}: {error:?}"));
                    calls += 1;
                    assert_eq!(count, expected, "{call_label}: literal output");
                    // Whole-input preservation after the same call.
                    assert_eq!(receiver.num_atoms(), rows, "{call_label}: rows preserved");
                    assert_eq!(receiver.num_bonds(), bonds, "{call_label}: bonds preserved");
                    assert_eq!(
                        receiver.properties(),
                        &properties,
                        "{call_label}: properties preserved"
                    );
                }
            }
        }
    }
    assert_eq!(calls, 120, "exact 120-call census");
}

/// Separately counted supplementary distinctions: the general recursive
/// count and the direct Lipinski N/O row sum disagree on the SAME real
/// molecules by design (frozen source-defined distinction, not fitted).
#[test]
fn descriptor_hba_public_supplementary_general_vs_lipinski() {
    let mut calls = 0usize;
    for remove_hydrogens in [false, true] {
        let label = format!("rh={remove_hydrogens}");
        let acid = sanitized("CC(=O)O", remove_hydrogens);
        assert_eq!(acid.num_hba().unwrap(), 1, "acid general {label}");
        calls += 1;
        assert_eq!(acid.lipinski_hba().unwrap(), 2, "acid direct {label}");
        calls += 1;
        let thiophene = sanitized("c1ccsc1", remove_hydrogens);
        assert_eq!(thiophene.num_hba().unwrap(), 1, "thiophene general {label}");
        calls += 1;
        assert_eq!(
            thiophene.lipinski_hba().unwrap(),
            0,
            "thiophene direct {label}"
        );
        calls += 1;
    }
    assert_eq!(calls, 8, "exact 8-call supplementary census");
}

/// Real missing-state regressions on PUBLIC calls: an unsanitized
/// molecule (raw graph, no prepared valence) reports the typed
/// MissingPreparedValence FIRST — before any ring state or chemistry is
/// consulted — with no fabricated cause, and the whole input survives
/// the failed call. This is the documented prepared-state boundary
/// limitation, not a chemical divergence.
#[test]
fn descriptor_hba_public_missing_prepared_valence() {
    let mut calls = 0usize;
    for (label, smiles) in [("raw-cco", "CCO"), ("pentavalent", "C(C)(C)(C)(C)C")] {
        let params = SmilesParseParams {
            sanitize: false,
            remove_hydrogens: false,
            ..Default::default()
        };
        let molecule = Molecule::from_smiles_with_params(smiles, &params).unwrap();
        let rows = molecule.num_atoms();
        let bonds = molecule.num_bonds();
        let properties = molecule.properties().clone();
        let error = molecule.num_hba().unwrap_err();
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
    assert_eq!(calls, 2, "exact missing-state census");
}
