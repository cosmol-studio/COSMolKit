//! Public descriptor count-query API regressions (DQ-PUBLIC).
#![cfg(all(feature = "cap-descriptors", feature = "cap-smiles"))]

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

/// HETERO-REPAIR: the frozen 320-call public product (20 cases x 2
/// remove-H policies x 2 raw/sanitized constructor states x 2
/// original/peer receivers x 2 repeats) for the topology-only
/// num_heteroatoms query. The 12-call four-block persistent-storage
/// proof (3 fixed molecules x 2 receivers x 2 repeats) is relocated
/// into the internal private unit module in crates/cosmolkit/src/
/// descriptors.rs (Arc identities are crate-private); the superseded
/// external value-equality version is removed here.
#[test]
fn descriptor_heteroatoms_public_product() {
    const HETERO_CASES: [(&str, u32); 20] = [
        ("", 0),
        ("C", 0),
        ("CCO", 1),
        ("[NH4+]", 1),
        ("[O-]", 1),
        ("[H][H]", 0),
        ("[2H]O[2H]", 1),
        ("[13CH4]", 0),
        ("N", 1),
        ("O", 1),
        ("C=O", 1),
        ("O=C(N)N", 3),
        ("n1ccccc1", 1),
        ("[nH]1cccc1", 1),
        ("C1CCCCC1", 0),
        ("CC(N)C(=O)O", 3),
        ("[H]N([H])[H]", 1),
        ("C[C+](C)C", 0),
        ("*", 1),
        ("CC#N", 1),
    ];

    // Frozen literal CONSTRUCTOR prerequisites (HETERO-REPAIR, from the
    // frozen SMILES/policy and the DQ literal table; never derived from
    // num_heteroatoms output). Ordered atomic numbers per case under both
    // remove-H policies, ordered isotopes, and explicit-H row
    // specifications. (smiles, ordered atomic numbers keep, ordered
    // isotopes keep, explicit-H rows keep)
    const PREREQS: [(&str, &[u8], &[Option<u16>], &[(usize, u8)]); 20] = [
        ("", &[], &[], &[]),
        ("C", &[6], &[None], &[]),
        ("CCO", &[6, 6, 8], &[None, None, None], &[]),
        ("[NH4+]", &[7], &[None], &[(0, 4)]),
        ("[O-]", &[8], &[None], &[]),
        ("[H][H]", &[1, 1], &[None, None], &[]),
        ("[2H]O[2H]", &[1, 8, 1], &[Some(2), None, Some(2)], &[]),
        ("[13CH4]", &[6], &[Some(13)], &[(0, 4)]),
        ("N", &[7], &[None], &[]),
        ("O", &[8], &[None], &[]),
        ("C=O", &[6, 8], &[None, None], &[]),
        ("O=C(N)N", &[8, 6, 7, 7], &[None, None, None, None], &[]),
        ("n1ccccc1", &[7, 6, 6, 6, 6, 6], &[None; 6], &[]),
        ("[nH]1cccc1", &[7, 6, 6, 6, 6], &[None; 5], &[(0, 1)]),
        ("C1CCCCC1", &[6, 6, 6, 6, 6, 6], &[None; 6], &[]),
        ("CC(N)C(=O)O", &[6, 6, 7, 6, 8, 8], &[None; 6], &[]),
        // Ammonia: H rows survive ONLY remove-H=false ([1,7,1,1]); the
        // remove policy yields just [7] with N explicit-H 3 AFTER
        // removal. All-hydrogen [H][H] keeps both rows under both
        // policies; deuterium rows survive both policies.
        ("[H]N([H])[H]", &[1, 7, 1, 1], &[None; 4], &[(1, 0)]),
        ("C[C+](C)C", &[6, 6, 6, 6], &[None; 4], &[]),
        ("*", &[0], &[None], &[]),
        ("CC#N", &[6, 6, 7], &[None, None, None], &[]),
    ];

    // Public-surface prerequisite check through Molecule::atoms() (rows,
    // atomic numbers, isotopes, explicit H) BEFORE every invocation; a
    // contradictory prerequisite fails visibly. Never SUT-derived.
    fn assert_input_prerequisites(
        label: &str,
        smiles: &str,
        atoms: &[cosmolkit_model::Atom],
        remove_hydrogens: bool,
    ) {
        let index = PREREQS
            .iter()
            .position(|(candidate, _, _, _)| *candidate == smiles)
            .unwrap_or_else(|| panic!("{label}: unknown prerequisite case {smiles:?}"));
        let (_, atomic_numbers, isotopes, explicit_h) = PREREQS[index];
        let all_hydrogen = atomic_numbers.iter().all(|&z| z == 1);
        let kept: Vec<(u8, Option<u16>)> = atomic_numbers
            .iter()
            .zip(isotopes.iter())
            .filter(|(z, isotope)| {
                !remove_hydrogens || **z != 1 || isotope.is_some() || all_hydrogen
            })
            .map(|(z, isotope)| (*z, *isotope))
            .collect();
        assert_eq!(atoms.len(), kept.len(), "{label}: atom row count");
        for (row, (expected_z, expected_isotope)) in atoms.iter().zip(kept.iter()) {
            assert_eq!(row.atomic_number(), *expected_z, "{label}: atomic number");
            assert_eq!(row.isotope(), *expected_isotope, "{label}: isotope");
        }
        if remove_hydrogens && smiles == "[H]N([H])[H]" {
            assert_eq!(atoms[0].atomic_number(), 7, "{label}: N kept");
            assert_eq!(atoms[0].explicit_hydrogens(), 3, "{label}: N H=3");
        } else {
            // All other explicit-H atom specifications remain the frozen
            // literals under both policies, not just remove-H=false.
            for &(row_index, expected_h) in explicit_h {
                assert_eq!(
                    atoms[row_index].explicit_hydrogens(),
                    expected_h,
                    "{label}: explicit H row {row_index}"
                );
            }
        }
    }

    let mut calls = 0usize;
    for (smiles, expected) in HETERO_CASES {
        for remove_hydrogens in [false, true] {
            for sanitized_constructor in [false, true] {
                // Truthful constructor-policy labeling: sanitize=false
                // builds the raw graph EXCEPT that remove-H=true still
                // runs the real RemoveHs pass (which sanitizes on its
                // own); this is the actual constructor policy, not an
                // assertion of a wholly unprepared graph.
                let params = SmilesParseParams {
                    sanitize: sanitized_constructor,
                    remove_hydrogens,
                    ..Default::default()
                };
                let original = Molecule::from_smiles_with_params(smiles, &params).unwrap();
                let peer = original.clone();
                for receiver in [&original, &peer] {
                    for _repeat in 0..2 {
                        let label =
                            format!("{smiles}/rh={remove_hydrogens}/s={sanitized_constructor}");
                        // BEFORE-call constructor/element/H/isotope
                        // prerequisite plus fresh public property clone.
                        assert_input_prerequisites(
                            &label,
                            smiles,
                            receiver.atoms(),
                            remove_hydrogens,
                        );
                        let properties_before = receiver.properties().clone();
                        let atoms = receiver.num_atoms();
                        let bonds = receiver.num_bonds();
                        let count = receiver
                            .num_heteroatoms()
                            .unwrap_or_else(|error| panic!("{label}: {error:?}"));
                        calls += 1;
                        assert_eq!(count, expected, "{label}");
                        assert_eq!(receiver.num_atoms(), atoms, "{label}: atoms");
                        assert_eq!(receiver.num_bonds(), bonds, "{label}: bonds");
                        assert_eq!(
                            receiver.properties(),
                            &properties_before,
                            "{label}: properties after the same call"
                        );
                    }
                }
            }
        }
    }
    assert_eq!(calls, 320, "exact census");
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

    // SmilesParse.cpp runs assignStereochemistry when sanitize OR removeHs
    // is true; Chirality.cpp fills a missing property cache non-strictly.
    // Retain all 40 exact constructor inputs and both policies. Only the
    // false/false branch lacks prepared valence. No query prepares a cache.
    let mut constructor_calls = 0usize;
    let mut prepared_calls = 0usize;
    let mut constructor_error_calls = 0usize;
    let mut error_calls = 0usize;
    for (smiles, _, _, _, expected, _, _, _) in CASES {
        for remove_hydrogens in [false, true] {
            let params = SmilesParseParams {
                sanitize: false,
                remove_hydrogens,
                ..Default::default()
            };
            let raw = Molecule::from_smiles_with_params(smiles, &params).unwrap();
            let original = raw.to_builder();
            let result = raw.total_atom_count();
            constructor_calls += 1;
            if remove_hydrogens {
                assert_eq!(
                    result.unwrap(),
                    expected,
                    "source-prepared raw {smiles} policy={remove_hydrogens}"
                );
                prepared_calls += 1;
            } else {
                let err = result.unwrap_err();
                assert!(
                    matches!(err, DescriptorReadError::MissingPreparedValence),
                    "raw {smiles} policy={remove_hydrogens}: {err:?}"
                );
                constructor_error_calls += 1;
            }
            assert_eq!(
                raw.to_builder(),
                original,
                "query preserves all detached blocks"
            );

            // Preserve all 40 original MissingPreparedValence controls on
            // canonical builder values, which carry the exact semantic blocks
            // and deliberately have no derived-cache installation authority.
            let unprepared = original.clone().build().unwrap();
            assert_eq!(unprepared.to_builder(), original);
            let err = unprepared.total_atom_count().unwrap_err();
            error_calls += 1;
            assert!(
                matches!(err, DescriptorReadError::MissingPreparedValence),
                "uncached {smiles} policy={remove_hydrogens}: {err:?}"
            );
            assert_eq!(
                unprepared.to_builder(),
                original,
                "failed query preserves all detached blocks"
            );
        }
    }
    assert_eq!(
        constructor_calls, 40,
        "all original constructor inputs retained"
    );
    assert_eq!(
        prepared_calls, 20,
        "source removeHs stereo preparation branch"
    );
    assert_eq!(
        constructor_error_calls, 20,
        "source false/false unprepared branch"
    );
    assert_eq!(
        error_calls, 40,
        "all 40 canonical uncached states typed-error"
    );
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

    // SmilesParse.cpp runs assignStereochemistry when sanitize OR removeHs
    // is true; Chirality.cpp fills a missing property cache non-strictly.
    // Retain all 40 exact constructor inputs and both policies. Only the
    // false/false branch lacks prepared valence. No query prepares a cache.
    let mut constructor_calls = 0usize;
    let mut prepared_calls = 0usize;
    let mut constructor_error_calls = 0usize;
    let mut error_calls = 0usize;
    for (smiles, _, _, _, _, _, expected, _) in CASES {
        for remove_hydrogens in [false, true] {
            let params = SmilesParseParams {
                sanitize: false,
                remove_hydrogens,
                ..Default::default()
            };
            let raw = Molecule::from_smiles_with_params(smiles, &params).unwrap();
            let original = raw.to_builder();
            let result = raw.lipinski_hbd();
            constructor_calls += 1;
            if remove_hydrogens {
                assert_eq!(
                    result.unwrap(),
                    expected,
                    "source-prepared raw {smiles} policy={remove_hydrogens}"
                );
                prepared_calls += 1;
            } else {
                let err = result.unwrap_err();
                assert!(
                    matches!(err, DescriptorReadError::MissingPreparedValence),
                    "raw {smiles} policy={remove_hydrogens}: {err:?}"
                );
                constructor_error_calls += 1;
            }
            assert_eq!(
                raw.to_builder(),
                original,
                "query preserves all detached blocks"
            );

            // Preserve all 40 original MissingPreparedValence controls on
            // canonical builder values, which carry the exact semantic blocks
            // and deliberately have no derived-cache installation authority.
            let unprepared = original.clone().build().unwrap();
            assert_eq!(unprepared.to_builder(), original);
            let err = unprepared.lipinski_hbd().unwrap_err();
            error_calls += 1;
            assert!(
                matches!(err, DescriptorReadError::MissingPreparedValence),
                "uncached {smiles} policy={remove_hydrogens}: {err:?}"
            );
            assert_eq!(
                unprepared.to_builder(),
                original,
                "failed query preserves all detached blocks"
            );
        }
    }
    assert_eq!(
        constructor_calls, 40,
        "all original constructor inputs retained"
    );
    assert_eq!(
        prepared_calls, 20,
        "source removeHs stereo preparation branch"
    );
    assert_eq!(
        constructor_error_calls, 20,
        "source false/false unprepared branch"
    );
    assert_eq!(
        error_calls, 40,
        "all 40 canonical uncached states typed-error"
    );
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

    // SmilesParse.cpp runs assignStereochemistry when sanitize OR removeHs
    // is true; Chirality.cpp fills a missing property cache non-strictly.
    // Retain all 40 exact constructor inputs and both policies. Only the
    // false/false branch lacks prepared valence. No query prepares a cache.
    let mut constructor_calls = 0usize;
    let mut prepared_calls = 0usize;
    let mut constructor_error_calls = 0usize;
    let mut error_calls = 0usize;
    for (smiles, _, _, _, _, _, _, expected) in CASES {
        for remove_hydrogens in [false, true] {
            let params = SmilesParseParams {
                sanitize: false,
                remove_hydrogens,
                ..Default::default()
            };
            let raw = Molecule::from_smiles_with_params(smiles, &params).unwrap();
            let original = raw.to_builder();
            let result = raw.fraction_csp3().map(f64::to_bits);
            constructor_calls += 1;
            if remove_hydrogens {
                assert_eq!(
                    result.unwrap(),
                    expected,
                    "source-prepared raw {smiles} policy={remove_hydrogens}"
                );
                prepared_calls += 1;
            } else {
                let err = result.unwrap_err();
                assert!(
                    matches!(err, DescriptorReadError::MissingPreparedValence),
                    "raw {smiles} policy={remove_hydrogens}: {err:?}"
                );
                constructor_error_calls += 1;
            }
            assert_eq!(
                raw.to_builder(),
                original,
                "query preserves all detached blocks"
            );

            // Preserve all 40 original MissingPreparedValence controls on
            // canonical builder values, which carry the exact semantic blocks
            // and deliberately have no derived-cache installation authority.
            let unprepared = original.clone().build().unwrap();
            assert_eq!(unprepared.to_builder(), original);
            let err = unprepared.fraction_csp3().unwrap_err();
            error_calls += 1;
            assert!(
                matches!(err, DescriptorReadError::MissingPreparedValence),
                "uncached {smiles} policy={remove_hydrogens}: {err:?}"
            );
            assert_eq!(
                unprepared.to_builder(),
                original,
                "failed query preserves all detached blocks"
            );
        }
    }
    assert_eq!(
        constructor_calls, 40,
        "all original constructor inputs retained"
    );
    assert_eq!(
        prepared_calls, 20,
        "source removeHs stereo preparation branch"
    );
    assert_eq!(
        constructor_error_calls, 20,
        "source false/false unprepared branch"
    );
    assert_eq!(
        error_calls, 40,
        "all 40 canonical uncached states typed-error"
    );
}
