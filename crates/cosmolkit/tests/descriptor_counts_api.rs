//! Public descriptor count-query API regressions (DQ-PUBLIC).

use cosmolkit::DescriptorReadError;
use cosmolkit::{Molecule, SmilesParseParams};

/// Frozen DQ-PUBLIC reference table (installed RDKit 2026.03.1,
/// sanitize=true, both removeHs policies, recorded BEFORE implementation):
/// (SMILES, explicit rows false/true, heavy, total, Lipinski HBA, Lipinski
/// HBD, CSP3 f64 bits). Every constructor receives the EXACT original
/// SMILES under both real parser policies.
const CASES: [(&str, usize, usize, u32, u32, u32, u32, u64); 20] = [
    ("", 0, 0, 0, 0, 0, 0, 0x0000_0000_0000_0000),
    ("C", 1, 1, 1, 5, 0, 0, 0x3ff0_0000_0000_0000),
    ("CCO", 3, 3, 3, 9, 1, 1, 0x3ff0_0000_0000_0000),
    ("[NH4+]", 1, 1, 1, 5, 1, 4, 0x0000_0000_0000_0000),
    ("[O-]", 1, 1, 1, 1, 1, 0, 0x0000_0000_0000_0000),
    ("[H][H]", 2, 2, 0, 2, 0, 0, 0x0000_0000_0000_0000),
    ("[2H]O[2H]", 3, 3, 1, 3, 1, 2, 0x0000_0000_0000_0000),
    ("[13CH4]", 1, 1, 1, 5, 0, 0, 0x3ff0_0000_0000_0000),
    ("N", 1, 1, 1, 4, 1, 3, 0x0000_0000_0000_0000),
    ("O", 1, 1, 1, 3, 1, 2, 0x0000_0000_0000_0000),
    ("C=O", 2, 2, 2, 4, 1, 0, 0x0000_0000_0000_0000),
    ("O=C(N)N", 4, 4, 4, 8, 3, 4, 0x0000_0000_0000_0000),
    ("n1ccccc1", 6, 6, 6, 11, 1, 0, 0x0000_0000_0000_0000),
    ("[nH]1cccc1", 5, 5, 5, 10, 1, 1, 0x0000_0000_0000_0000),
    ("C1CCCCC1", 6, 6, 6, 18, 0, 0, 0x3ff0_0000_0000_0000),
    ("CC(N)C(=O)O", 6, 6, 6, 13, 3, 3, 0x3fe5_5555_5555_5555),
    ("[H]N([H])[H]", 4, 1, 1, 4, 1, 3, 0x0000_0000_0000_0000),
    ("C[C+](C)C", 4, 4, 4, 13, 0, 0, 0x3fe8_0000_0000_0000),
    ("*", 1, 1, 0, 1, 0, 0, 0x0000_0000_0000_0000),
    ("CC#N", 3, 3, 3, 6, 1, 0, 0x3fe0_0000_0000_0000),
];

fn sanitized(smiles: &str, remove_hydrogens: bool) -> Molecule {
    let params = SmilesParseParams {
        remove_hydrogens,
        ..Default::default()
    };
    Molecule::from_smiles_with_params(smiles, &params).unwrap()
}

#[test]
fn descriptor_query_dq04_num_heavy_atoms() {
    // Frozen 160-call product: 20 cases x 2 remove-H policies x 2 receivers
    // (original + Arc-sharing peer clone) x 2 consecutive calls; bitwise
    // results and atom/bond/property public snapshots retained after every
    // call. The census increments only at an actual method invocation.
    let mut calls = 0usize;
    for (smiles, rows_false, rows_true, heavy, _, _, _, _) in CASES {
        for remove_hydrogens in [false, true] {
            let expected_rows = if remove_hydrogens {
                rows_true
            } else {
                rows_false
            };
            let original = sanitized(smiles, remove_hydrogens);
            assert_eq!(
                original.num_atoms(),
                expected_rows,
                "{smiles} policy={remove_hydrogens} row count"
            );
            let peer = original.clone();
            for (label, receiver) in [("original", &original), ("peer", &peer)] {
                for repeat in 0..2 {
                    let atoms = receiver.num_atoms();
                    let bonds = receiver.num_bonds();
                    let properties = receiver.properties().clone();
                    let result = receiver.num_heavy_atoms().unwrap();
                    calls += 1;
                    assert_eq!(
                        result, heavy,
                        "{smiles} policy={remove_hydrogens} {label} #{repeat}"
                    );
                    assert_eq!(receiver.num_atoms(), atoms, "{label} atom snapshot");
                    assert_eq!(receiver.num_bonds(), bonds, "{label} bond snapshot");
                    assert_eq!(
                        receiver.properties(),
                        &properties,
                        "{label} property snapshot"
                    );
                }
            }
        }
    }
    assert_eq!(calls, 160, "exact 160-call census");

    // Raw-input no-valence discriminator (separately counted): the
    // unsanitized constructor leaves no prepared valence, but the
    // topology-only heavy query stays callable with the same literals.
    let mut raw_calls = 0usize;
    for (smiles, _, _, heavy, _, _, _, _) in CASES {
        for remove_hydrogens in [false, true] {
            let params = SmilesParseParams {
                sanitize: false,
                remove_hydrogens,
                ..Default::default()
            };
            let raw = Molecule::from_smiles_with_params(smiles, &params).unwrap();
            assert_eq!(
                raw.num_heavy_atoms().unwrap(),
                heavy,
                "raw {smiles} policy={remove_hydrogens}"
            );
            raw_calls += 1;
        }
    }
    assert_eq!(raw_calls, 40, "raw discriminator separately counted");
}

#[test]
fn descriptor_query_dq05_total_atom_count() {
    // Frozen 160-call product for total_atom_count (rows + attached
    // hydrogens, includeNeighbors=false): explicit-H non-double-counting
    // is retained by the [H]N([H])[H] row (policy false: 4 rows -> total
    // 4, exactly the 3 explicit H atom rows plus the 1 N row; policy
    // true: 1 row -> total 4) and the row-count distinction by asserting
    // Molecule::num_atoms (row length) separately from the total.
    let mut calls = 0usize;
    for (smiles, rows_false, rows_true, _, total, _, _, _) in CASES {
        for remove_hydrogens in [false, true] {
            let expected_rows = if remove_hydrogens {
                rows_true
            } else {
                rows_false
            };
            let original = sanitized(smiles, remove_hydrogens);
            assert_eq!(
                original.num_atoms(),
                expected_rows,
                "{smiles} policy={remove_hydrogens} row count"
            );
            let peer = original.clone();
            for (label, receiver) in [("original", &original), ("peer", &peer)] {
                for repeat in 0..2 {
                    let atoms = receiver.num_atoms();
                    let bonds = receiver.num_bonds();
                    let properties = receiver.properties().clone();
                    let result = receiver.total_atom_count().unwrap();
                    calls += 1;
                    assert_eq!(
                        result, total,
                        "{smiles} policy={remove_hydrogens} {label} #{repeat}"
                    );
                    assert_eq!(
                        receiver.num_atoms(),
                        atoms,
                        "{label} row-length accessor unchanged"
                    );
                    // Row-length accessor keeps its own meaning: wherever
                    // the frozen table has rows != total (C 1->5, CCO
                    // 3->9, ...), the literal total assertion above is the
                    // distinction; rows == total holds only where the
                    // frozen reference says so (empty, [O-], [H][H],
                    // [2H]O[2H]).
                    assert_eq!(
                        u32::try_from(atoms).unwrap(),
                        if remove_hydrogens {
                            rows_true
                        } else {
                            rows_false
                        } as u32,
                        "{label} num_atoms stays the row-length accessor"
                    );
                    assert_eq!(receiver.num_bonds(), bonds, "{label} bond snapshot");
                    assert_eq!(
                        receiver.properties(),
                        &properties,
                        "{label} property snapshot"
                    );
                }
            }
        }
    }
    assert_eq!(calls, 160, "exact 160-call census");

    // All 40 unsanitized constructor states: exactly MissingPreparedValence
    // — no silent cache preparation, even for empty/zero-carbon cases.
    let mut error_calls = 0usize;
    for (smiles, _, _, _, _, _, _, _) in CASES {
        for remove_hydrogens in [false, true] {
            let params = SmilesParseParams {
                sanitize: false,
                remove_hydrogens,
                ..Default::default()
            };
            let raw = Molecule::from_smiles_with_params(smiles, &params).unwrap();
            let err = raw.total_atom_count().unwrap_err();
            error_calls += 1;
            assert!(
                matches!(err, DescriptorReadError::MissingPreparedValence),
                "raw {smiles} policy={remove_hydrogens}: expected MissingPreparedValence, got {err:?}"
            );
        }
    }
    assert_eq!(error_calls, 40, "all 40 unsanitized states typed-error");
}

#[test]
fn descriptor_query_dq06_lipinski_hba() {
    // Frozen 160-call Lipinski HBA product: the DIRECT N/O row count, not
    // general recursive NumHBA. Frozen discriminators retained on every
    // source row: quaternary [NH4+] still counts (general HBA would not),
    // urea O=C(N)N has HBA 3, aromatic [nH]1cccc1 has HBA 1, C=O has 1,
    // and pure hydrocarbons have 0.
    let mut calls = 0usize;
    for (smiles, rows_false, rows_true, _, _, hba, _, _) in CASES {
        for remove_hydrogens in [false, true] {
            let expected_rows = if remove_hydrogens {
                rows_true
            } else {
                rows_false
            };
            let original = sanitized(smiles, remove_hydrogens);
            assert_eq!(
                original.num_atoms(),
                expected_rows,
                "{smiles} policy={remove_hydrogens} row count"
            );
            let peer = original.clone();
            for (label, receiver) in [("original", &original), ("peer", &peer)] {
                for repeat in 0..2 {
                    let atoms = receiver.num_atoms();
                    let bonds = receiver.num_bonds();
                    let properties = receiver.properties().clone();
                    let result = receiver.lipinski_hba().unwrap();
                    calls += 1;
                    assert_eq!(
                        result, hba,
                        "{smiles} policy={remove_hydrogens} {label} #{repeat}"
                    );
                    assert_eq!(receiver.num_atoms(), atoms, "{label} atom snapshot");
                    assert_eq!(receiver.num_bonds(), bonds, "{label} bond snapshot");
                    assert_eq!(
                        receiver.properties(),
                        &properties,
                        "{label} property snapshot"
                    );
                }
            }
        }
    }
    assert_eq!(calls, 160, "exact 160-call census");

    // Topology-only raw discriminator (separately counted): the direct N/O
    // count stays callable on unsanitized states with the same literals.
    let mut raw_calls = 0usize;
    for (smiles, _, _, _, _, hba, _, _) in CASES {
        for remove_hydrogens in [false, true] {
            let params = SmilesParseParams {
                sanitize: false,
                remove_hydrogens,
                ..Default::default()
            };
            let raw = Molecule::from_smiles_with_params(smiles, &params).unwrap();
            assert_eq!(
                raw.lipinski_hba().unwrap(),
                hba,
                "raw {smiles} policy={remove_hydrogens}"
            );
            raw_calls += 1;
        }
    }
    assert_eq!(raw_calls, 40, "raw discriminator separately counted");
}

#[test]
fn descriptor_query_dq07_lipinski_hbd() {
    // Frozen 160-call Lipinski donor-HYDROGEN-sum product (hydrogens on
    // N/O rows, not the donor-atom count). Deuterium and ammonium frozen
    // expectations are retained on every receiver: [2H]O[2H]=2 (isotopic
    // H neighbors count), [NH4+]=4, [H]N([H])[H]=3 under BOTH policies,
    // [nH]1cccc1=1, C=O=0.
    let mut calls = 0usize;
    for (smiles, rows_false, rows_true, _, _, _, hbd, _) in CASES {
        for remove_hydrogens in [false, true] {
            let expected_rows = if remove_hydrogens {
                rows_true
            } else {
                rows_false
            };
            let original = sanitized(smiles, remove_hydrogens);
            assert_eq!(
                original.num_atoms(),
                expected_rows,
                "{smiles} policy={remove_hydrogens} row count"
            );
            let peer = original.clone();
            for (label, receiver) in [("original", &original), ("peer", &peer)] {
                for repeat in 0..2 {
                    let atoms = receiver.num_atoms();
                    let bonds = receiver.num_bonds();
                    let properties = receiver.properties().clone();
                    let result = receiver.lipinski_hbd().unwrap();
                    calls += 1;
                    assert_eq!(
                        result, hbd,
                        "{smiles} policy={remove_hydrogens} {label} #{repeat}"
                    );
                    assert_eq!(receiver.num_atoms(), atoms, "{label} atom snapshot");
                    assert_eq!(receiver.num_bonds(), bonds, "{label} bond snapshot");
                    assert_eq!(
                        receiver.properties(),
                        &properties,
                        "{label} property snapshot"
                    );
                }
            }
        }
    }
    assert_eq!(calls, 160, "exact 160-call census");

    // All 40 unsanitized states: exactly MissingPreparedValence (the
    // hydrogen sum needs the prepared assignment; no silent preparation).
    let mut error_calls = 0usize;
    for (smiles, _, _, _, _, _, _, _) in CASES {
        for remove_hydrogens in [false, true] {
            let params = SmilesParseParams {
                sanitize: false,
                remove_hydrogens,
                ..Default::default()
            };
            let raw = Molecule::from_smiles_with_params(smiles, &params).unwrap();
            let err = raw.lipinski_hbd().unwrap_err();
            error_calls += 1;
            assert!(
                matches!(err, DescriptorReadError::MissingPreparedValence),
                "raw {smiles} policy={remove_hydrogens}: expected MissingPreparedValence, got {err:?}"
            );
        }
    }
    assert_eq!(error_calls, 40, "all 40 unsanitized states typed-error");
}

#[test]
fn descriptor_query_dq08_fraction_csp3() {
    // Frozen 160-call exact-bits CSP3 product. Every constructor receives
    // the EXACT original frozen SMILES under both real parser policies
    // (remove_hydrogens routes through the parser's existing RemoveHs
    // owner; no replacement string). Zero-carbon rows return +0.0 bits
    // and the charged carbon C[C+](C)C yields 0x3fe8000000000000 (3/4).
    let mut calls = 0usize;
    for (smiles, rows_false, rows_true, _, _, _, _, csp3_bits) in CASES {
        for remove_hydrogens in [false, true] {
            let expected_rows = if remove_hydrogens {
                rows_true
            } else {
                rows_false
            };
            let original = sanitized(smiles, remove_hydrogens);
            assert_eq!(
                original.num_atoms(),
                expected_rows,
                "{smiles} policy={remove_hydrogens} row count"
            );
            let peer = original.clone();
            for (label, receiver) in [("original", &original), ("peer", &peer)] {
                for repeat in 0..2 {
                    let atoms = receiver.num_atoms();
                    let bonds = receiver.num_bonds();
                    let properties = receiver.properties().clone();
                    let result = receiver.fraction_csp3().unwrap();
                    calls += 1;
                    assert_eq!(
                        result.to_bits(),
                        csp3_bits,
                        "{smiles} policy={remove_hydrogens} {label} #{repeat}"
                    );
                    assert_eq!(receiver.num_atoms(), atoms, "{label} atom snapshot");
                    assert_eq!(receiver.num_bonds(), bonds, "{label} bond snapshot");
                    assert_eq!(
                        receiver.properties(),
                        &properties,
                        "{label} property snapshot"
                    );
                }
            }
        }
    }
    assert_eq!(calls, 160, "exact 160-call census");

    // All 40 unsanitized states: exactly MissingPreparedValence, including
    // zero-carbon and charged-carbon cases.
    let mut error_calls = 0usize;
    for (smiles, _, _, _, _, _, _, _) in CASES {
        for remove_hydrogens in [false, true] {
            let params = SmilesParseParams {
                sanitize: false,
                remove_hydrogens,
                ..Default::default()
            };
            let raw = Molecule::from_smiles_with_params(smiles, &params).unwrap();
            let err = raw.fraction_csp3().unwrap_err();
            error_calls += 1;
            assert!(
                matches!(err, DescriptorReadError::MissingPreparedValence),
                "raw {smiles} policy={remove_hydrogens}: expected MissingPreparedValence, got {err:?}"
            );
        }
    }
    assert_eq!(error_calls, 40, "all 40 unsanitized states typed-error");
}
