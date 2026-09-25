use cosmolkit_core::KekulizeError;
use cosmolkit_model::BondOrder;
use cosmolkit_smiles::{
    CxSmilesFields, CxSmilesWriteParams, SmilesParseError, SmilesParseParams, SmilesWriteParams,
    finalize_smiles_stereo, parse_smiles, write_cx_smiles_with_params, write_smiles,
    write_smiles_with_params,
};

fn write(input: &str, params: &SmilesWriteParams) -> String {
    let record = parse_smiles(input, &Default::default()).expect("parse pinned case");
    write_smiles_with_params(&record, params).expect("write pinned case")
}

fn renumbered_record(
    record: &cosmolkit_smiles::SmilesRecord,
    old_atom_order: &[usize],
) -> cosmolkit_smiles::SmilesRecord {
    let old_atom_order = old_atom_order
        .iter()
        .copied()
        .map(cosmolkit_model::AtomId::new)
        .collect::<Vec<_>>();
    let mut renumbered = record.clone();
    renumbered.topology = record
        .topology
        .reordered_atoms(&old_atom_order)
        .expect("renumber source-equivalent graph")
        .0;
    renumbered
}

#[test]
fn writer_kekulization_error_preserves_typed_core_source() {
    let mut record = parse_smiles("C", &Default::default()).expect("parse carbon");
    record.topology.atoms[0].set_aromatic(true);
    let error = write_smiles_with_params(
        &record,
        &SmilesWriteParams {
            do_kekule: true,
            ..Default::default()
        },
    )
    .expect_err("aromatic atom outside a ring cannot be kekulized");

    assert!(matches!(
        &error,
        SmilesParseError::WriterKekulize(KekulizeError::AromaticAtomOutsideRing { atom })
            if atom.index() == 0
    ));
    let source = std::error::Error::source(&error)
        .expect("writer error preserves the core kekulization source");
    assert!(source.downcast_ref::<KekulizeError>().is_some());
}

#[test]
fn writer_disconnected_kekulization_error_preserves_fragment_local_source_and_input() {
    let parser = SmilesParseParams {
        sanitize: false,
        remove_hydrogens: false,
        ..Default::default()
    };
    let mut record = parse_smiles("c1ccccc1.C", &parser).expect("parse unsanitized fragments");
    record.topology.atoms[6].set_aromatic(true);
    let before = record.clone();
    let params = SmilesWriteParams {
        do_isomeric_smiles: true,
        do_kekule: true,
        canonical: true,
        clean_stereo: true,
        ..Default::default()
    };

    let error = write_smiles_with_params(&record, &params)
        .expect_err("the non-ring atom in the second component cannot be kekulized");
    assert!(matches!(
        &error,
        SmilesParseError::WriterKekulize(KekulizeError::AromaticAtomOutsideRing { atom })
            if atom.index() == 0
    ));
    let source = std::error::Error::source(&error)
        .expect("writer error preserves the typed component-local core source");
    assert!(matches!(
        source.downcast_ref::<KekulizeError>(),
        Some(KekulizeError::AromaticAtomOutsideRing { atom }) if atom.index() == 0
    ));
    assert_eq!(
        record, before,
        "failed writing must leave its input unchanged"
    );
}

// Oracle: pinned RDKit 2026.03.1, source revision
// 351f8f378f8ad6bbd517980c38896e66bf907af8c. Expected values were verified
// against that installed version after reading SmilesWrite.cpp and its tests;
// no historical Rust output or generated expected-data cache is used.
#[test]
fn atom_and_bond_tokens_match_pinned_rdkit() {
    let params = SmilesWriteParams {
        canonical: false,
        ..Default::default()
    };
    for (input, expected) in [
        ("CC", "CC"),
        ("C=C", "C=C"),
        ("C#N", "C#N"),
        ("C$C", "C$C"),
        ("c1ccccc1", "c1ccccc1"),
        ("C:C", "C:C"),
        ("N->B", "N->B"),
        ("C~C", "C~C"),
        ("[NH4+]", "[NH4+]"),
        ("[13CH4]", "[13CH4]"),
    ] {
        let record = parse_smiles(input, &Default::default()).expect("parse pinned case");
        assert_eq!(
            write_smiles_with_params(&record, &params).unwrap(),
            expected,
            "{input}"
        );
    }

    let canonical = parse_smiles("N->B", &Default::default()).unwrap();
    assert_eq!(write_smiles(&canonical).unwrap(), "B<-N");
}

#[test]
fn bond_token_endpoints_and_canonical_slash_reversal_match_pinned_rdkit() {
    // Oracle: pinned RDKit 2026.03.1 / revision
    // 351f8f378f8ad6bbd517980c38896e66bf907af8c. Profiles and exact states
    // are frozen in S17_bond_tokens.md.
    let canonical = SmilesWriteParams::default();
    let source_order = SmilesWriteParams {
        canonical: false,
        ..Default::default()
    };
    for (input, expected) in [
        ("N->B", "B<-N"),
        ("B<-N", "B<-N"),
        (r"Br\C=C/F", r"F/C=C\Br"),
    ] {
        assert_eq!(write(input, &canonical), expected, "canonical {input}");
    }
    for (input, expected) in [("N->B", "N->B"), ("B<-N", "B<-N")] {
        assert_eq!(
            write(input, &source_order),
            expected,
            "source order {input}"
        );
    }
    let reversed_slash = parse_smiles(r"Br\C=C/F", &Default::default()).unwrap();
    assert_eq!(
        write_smiles_with_params(&reversed_slash, &source_order).unwrap(),
        r"Br/C=C\F",
        "source-order traversal normalizes the two endpoint directions"
    );

    let mut zero_order = parse_smiles("CC", &Default::default()).unwrap();
    zero_order.topology.bonds[0].set_order(BondOrder::Zero);
    let before = zero_order.clone();
    assert_eq!(
        write_smiles_with_params(&zero_order, &source_order).unwrap(),
        "C~C"
    );
    assert_eq!(zero_order, before, "writer mutated the zero-order input");
}

#[test]
fn atom_bracket_predicate_matches_pinned_rdkit_branches() {
    // Oracle: RDKit 2026.03.1 / revision
    // 351f8f378f8ad6bbd517980c38896e66bf907af8c, MolToSmiles defaults
    // (canonical=true, doIsomericSmiles=true, allHsExplicit=false). The two
    // platinum cases are also pinned in SmilesParse/catch_tests.cpp under
    // "atoms bound to metals should always have Hs specified". The other
    // branch inputs were probed against this installed pinned build after
    // reading SmilesWrite.cpp::atomNeedsBracket.
    let params = SmilesWriteParams::default();
    for (input, expected) in [
        ("C", "C"),
        ("[SiH4]", "[SiH4]"),
        ("[NH4+]", "[NH4+]"),
        ("[13CH4]", "[13CH4]"),
        ("[CH4:7]", "[CH4:7]"),
        ("[CH3]", "[CH3]"),
        ("c1cc[nH]c1", "c1cc[nH]c1"),
        ("Cl[Pt](F)([NH2])[OH]", "[NH2][Pt]([OH])([F])[Cl]"),
        ("Cl[Pt](F)(<-[NH3])[OH]", "[NH3]->[Pt]([OH])([F])[Cl]"),
    ] {
        assert_eq!(write(input, &params), expected, "{input}");
    }
}

#[test]
fn writer_flags_control_isotope_direction_bond_hydrogen_and_kekule_tokens() {
    assert_eq!(
        write(
            "[13CH4]",
            &SmilesWriteParams {
                do_isomeric_smiles: false,
                ..Default::default()
            }
        ),
        "C"
    );
    assert_eq!(
        write(
            "C/C=C/C",
            &SmilesWriteParams {
                do_isomeric_smiles: false,
                ..Default::default()
            }
        ),
        "CC=CC"
    );
    assert_eq!(
        write(
            "CC",
            &SmilesWriteParams {
                all_bonds_explicit: true,
                ..Default::default()
            }
        ),
        "C-C"
    );
    assert_eq!(
        write(
            "N->B",
            &SmilesWriteParams {
                canonical: false,
                include_dative_bonds: false,
                ..Default::default()
            }
        ),
        "[NH3]B"
    );
    assert_eq!(
        write(
            "N->B",
            &SmilesWriteParams {
                canonical: false,
                all_bonds_explicit: true,
                include_dative_bonds: false,
                ..Default::default()
            }
        ),
        "[NH3]-B"
    );
    assert_eq!(
        write(
            "CC",
            &SmilesWriteParams {
                all_hydrogens_explicit: true,
                ..Default::default()
            }
        ),
        "[CH3][CH3]"
    );
    assert_eq!(
        write(
            "c1ccccc1",
            &SmilesWriteParams {
                do_kekule: true,
                ..Default::default()
            }
        ),
        "C1=CC=CC=C1"
    );
    assert_eq!(
        write(
            "c1cc[nH]c1",
            &SmilesWriteParams {
                do_kekule: true,
                ..Default::default()
            }
        ),
        "C1=CNC=C1"
    );
}

#[test]
fn aromatic_endpoint_bond_omission_matches_pinned_rdkit() {
    // Oracle: pinned RDKit 2026.03.1 / revision
    // 351f8f378f8ad6bbd517980c38896e66bf907af8c. The source probe used one
    // unsanitized bond, assigned endpoint/bond aromatic flags, refreshed the
    // property cache, and changed only allBondsExplicit between writes.
    let cases = [
        ("CC", [true, true], BondOrder::Single, false, "c-c", "c-c"),
        ("CC", [true, false], BondOrder::Single, false, "cC", "c-C"),
        ("CC", [false, true], BondOrder::Single, false, "Cc", "C-c"),
        ("CC", [false, false], BondOrder::Single, false, "CC", "C-C"),
        ("CC", [true, true], BondOrder::Single, true, "cc", "c-c"),
        ("CC", [true, true], BondOrder::Double, true, "cc", "c=c"),
        ("CC", [true, true], BondOrder::Double, false, "c=c", "c=c"),
        ("CC", [true, true], BondOrder::Aromatic, true, "cc", "c:c"),
        ("CC", [true, false], BondOrder::Aromatic, true, "c:C", "c:C"),
        ("CC", [false, true], BondOrder::Aromatic, true, "C:c", "C:c"),
        (
            "CC",
            [false, false],
            BondOrder::Aromatic,
            true,
            "C:C",
            "C:C",
        ),
        ("**", [true, true], BondOrder::Aromatic, true, "*:*", "*:*"),
        ("*C", [true, true], BondOrder::Aromatic, true, "*c", "*:c"),
    ];

    let default_params = SmilesWriteParams {
        canonical: false,
        ..Default::default()
    };
    let explicit_params = SmilesWriteParams {
        canonical: false,
        all_bonds_explicit: true,
        ..Default::default()
    };
    for (input, aromatic_atoms, order, aromatic_bond, expected, expected_explicit) in cases {
        let mut record = parse_smiles(input, &Default::default()).expect("parse two-atom case");
        record.topology.atoms[0].set_aromatic(aromatic_atoms[0]);
        record.topology.atoms[1].set_aromatic(aromatic_atoms[1]);
        record.topology.bonds[0].set_order(order);
        record.topology.bonds[0].set_aromatic(aromatic_bond);
        let before = record.clone();

        assert_eq!(
            write_smiles_with_params(&record, &default_params).unwrap(),
            expected,
            "{input}, atoms={aromatic_atoms:?}, order={order:?}, bond_aromatic={aromatic_bond}"
        );
        assert_eq!(
            write_smiles_with_params(&record, &explicit_params).unwrap(),
            expected_explicit,
            "allBondsExplicit: {input}, atoms={aromatic_atoms:?}, order={order:?}, bond_aromatic={aromatic_bond}"
        );
        assert_eq!(record, before, "writer mutated input state for {input}");
    }
}

#[test]
fn atom_hydrogen_charge_and_isotope_fields_match_pinned_rdkit() {
    // Pinned RDKit 2026.03.1 / revision
    // 351f8f378f8ad6bbd517980c38896e66bf907af8c. MolFromSmiles uses its
    // defaults; MolToSmiles uses the default legacy canonical/isomeric profile
    // except where the listed allHsExplicit/doIsomericSmiles option changes.
    let default_params = SmilesWriteParams::default();
    for (input, explicit_hydrogens, implicit_hydrogens) in [("N", 0, 3), ("[NH3]", 3, 0)] {
        let record = parse_smiles(input, &Default::default()).unwrap();
        let before = record.clone();
        assert_eq!(
            record.topology.atoms[0].explicit_hydrogens(),
            explicit_hydrogens,
            "explicit-H source count for {input}"
        );
        let valence = cosmolkit_core::assign_valence_with_options_for_topology(
            &record.topology,
            cosmolkit_core::ValenceModel::RdkitLike,
            false,
        )
        .unwrap();
        assert_eq!(
            valence.implicit_hydrogens[0], implicit_hydrogens,
            "implicit-H source count for {input}"
        );
        assert_eq!(
            write_smiles_with_params(&record, &default_params).unwrap(),
            "N",
            "default H suppression for {input}"
        );
        assert_eq!(
            write_smiles_with_params(
                &record,
                &SmilesWriteParams {
                    all_hydrogens_explicit: true,
                    ..Default::default()
                }
            )
            .unwrap(),
            "[NH3]",
            "all-H bracket count for {input}"
        );
        assert_eq!(record, before, "writing preserves {input}");
    }

    // Charge is emitted regardless of isomeric output; positive and negative
    // magnitudes use the pinned source's distinct +/-1 and +/-N branches.
    for (input, expected) in [
        ("[NH4+]", "[NH4+]"),
        ("[Cl-]", "[Cl-]"),
        ("[Mg+2]", "[Mg+2]"),
        ("[O-2]", "[O-2]"),
    ] {
        let record = parse_smiles(input, &Default::default()).unwrap();
        let before = record.clone();
        for do_isomeric_smiles in [true, false] {
            let params = SmilesWriteParams {
                do_isomeric_smiles,
                ..Default::default()
            };
            assert_eq!(
                write_smiles_with_params(&record, &params).unwrap(),
                expected,
                "charge formatting for {input}, isomeric={do_isomeric_smiles}"
            );
        }
        assert_eq!(record, before, "writing preserves {input}");
    }

    let isotope = parse_smiles("[13CH4]", &Default::default()).unwrap();
    let before = isotope.clone();
    for (params, expected) in [
        (default_params, "[13CH4]"),
        (
            SmilesWriteParams {
                do_isomeric_smiles: false,
                ..Default::default()
            },
            "C",
        ),
        (
            SmilesWriteParams {
                do_isomeric_smiles: false,
                all_hydrogens_explicit: true,
                ..Default::default()
            },
            "[CH4]",
        ),
    ] {
        assert_eq!(
            write_smiles_with_params(&isotope, &params).unwrap(),
            expected,
            "isotope/H formatting with {params:?}"
        );
    }
    assert_eq!(isotope, before, "writing preserves isotope source state");
}

#[test]
fn aromatic_custom_symbols_and_bracket_labels_match_pinned_rdkit() {
    // Pinned RDKit 2026.03.1 / revision
    // 351f8f378f8ad6bbd517980c38896e66bf907af8c, MolToSmiles defaults with
    // only doKekule toggled for the two aromatic custom-symbol cases.
    let mut aromatic = parse_smiles("c1ccccc1", &Default::default()).unwrap();
    aromatic.topology.atoms[0]
        .set_prop("smilesSymbol", "X")
        .unwrap();
    let before = aromatic.clone();
    for (do_kekule, expected) in [(false, "c1c[xH]ccc1"), (true, "C1=C[XH]=CC=C1")] {
        let params = SmilesWriteParams {
            do_kekule,
            ..Default::default()
        };
        assert_eq!(
            write_smiles_with_params(&aromatic, &params).unwrap(),
            expected,
            "do_kekule={do_kekule}"
        );
    }
    assert_eq!(aromatic, before, "writing preserves custom aromatic input");

    let mut bracketed_label = parse_smiles("[13CH4:7]", &Default::default()).unwrap();
    bracketed_label.topology.atoms[0]
        .set_prop("_supplementalSmilesLabel", " ;|{}\\")
        .unwrap();
    assert_eq!(write_smiles(&bracketed_label).unwrap(), "[13CH4:7] ;|{}\\");

    let mut custom_label = parse_smiles("C", &Default::default()).unwrap();
    custom_label.topology.atoms[0]
        .set_prop("smilesSymbol", "X")
        .unwrap();
    custom_label.topology.atoms[0]
        .set_prop("_supplementalSmilesLabel", "TAG")
        .unwrap();
    assert_eq!(write_smiles(&custom_label).unwrap(), "[XH4]TAG");
}

#[test]
fn custom_symbols_and_supplemental_labels_are_appended_without_escaping() {
    let mut custom = parse_smiles("C", &Default::default()).unwrap();
    custom.topology.atoms[0]
        .set_prop("smilesSymbol", "X]\\")
        .unwrap();
    assert_eq!(write_smiles(&custom).unwrap(), "[X]\\H4]");

    let mut supplemental = parse_smiles("C", &Default::default()).unwrap();
    supplemental.topology.atoms[0]
        .set_prop("_supplementalSmilesLabel", " ;|{}\\")
        .unwrap();
    assert_eq!(write_smiles(&supplemental).unwrap(), "C ;|{}\\");
}

#[test]
fn traversal_emits_source_order_branches_and_fragment_order() {
    let noncanonical = SmilesWriteParams {
        canonical: false,
        ..Default::default()
    };
    assert_eq!(write("C(O)(N)F", &noncanonical), "C(O)(N)F");
    assert_eq!(write("O.CC", &noncanonical), "O.CC");
    assert_eq!(write("C1(C2CC2)CC1.O", &noncanonical), "C1(C2CC2)CC1.O");

    // Pinned `detail::MolToSmiles` ranks/traverses canonically and sorts
    // disconnected fragment tuples by their serialized text.
    assert_eq!(write("C(O)(N)F", &SmilesWriteParams::default()), "NC(O)F");
    assert_eq!(write("O.CC", &SmilesWriteParams::default()), "CC.O");
}

#[test]
fn ring_labels_follow_closure_order_and_reuse_the_lowest_display_slot() {
    let noncanonical = SmilesWriteParams {
        canonical: false,
        ..Default::default()
    };
    for (input, expected) in [
        ("C1CC2CCC1C2", "C1CC2CCC1C2"),
        ("C1CCC2(CC1)CCCC2", "C1CCC2(CC1)CCCC2"),
        ("C12C3C4C1C5C2C3C45", "C12C3C4C1C1C2C3C41"),
    ] {
        assert_eq!(write(input, &noncanonical), expected, "{input}");
    }
}

#[test]
fn writer_ring_label_tokens_match_pinned_rdkit_at_9_10_99_100() {
    fn ring_token(index: usize) -> String {
        if index < 10 {
            index.to_string()
        } else if index < 100 {
            format!("%{index}")
        } else {
            format!("%({index})")
        }
    }

    // Pinned RDKit 2026.03.1, explicit legacy source profile, canonical=false,
    // doIsomericSmiles=false, rootedAtAtom=0, sanitize=false, and the other
    // writer defaults. The source emits these root-preserving fixtures
    // byte-for-byte unchanged; see S18_ring_labels.md for the oracle probe.
    let parser = SmilesParseParams {
        sanitize: false,
        remove_hydrogens: false,
        ..Default::default()
    };
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: false,
        rooted_at_atom: Some(cosmolkit_model::AtomId::new(0)),
        ..Default::default()
    };

    for (count, output_bytes) in [(9, 29), (10, 36), (99, 659), (100, 672)] {
        let mut input = String::from("*");
        for index in 1..=count {
            input.push_str(&ring_token(index));
        }
        input.push('*');
        for index in 1..=count {
            input.push('*');
            input.push_str(&ring_token(index));
        }

        let record = parse_smiles(&input, &parser).expect("parse pinned ring-label fixture");
        assert_eq!(record.topology.atoms.len(), count + 2, "count={count}");
        assert_eq!(record.topology.bonds.len(), 2 * count + 1, "count={count}");
        let before = record.clone();
        let output = write_smiles_with_params(&record, &params).unwrap();
        assert_eq!(output, input, "count={count}");
        assert_eq!(output.len(), output_bytes, "count={count}");
        assert_eq!(record, before, "writer mutated count={count} input");
    }
}

#[test]
fn dative_ring_closure_keeps_source_traversal_orientation() {
    assert_eq!(write("N1CCB->1", &SmilesWriteParams::default()), "B1CCN<-1");
    assert_eq!(write("N1CCB<-1", &SmilesWriteParams::default()), "B1CCN->1");
}

#[test]
fn equal_canonical_fragments_sort_by_source_atom_and_bond_order() {
    let mut record = parse_smiles("C(O)(N)F.C(O)(N)F", &Default::default()).unwrap();
    let new_atom_order = [0, 6, 2, 4, 1, 3, 5, 7].map(cosmolkit_model::AtomId::new);
    record.topology = record.topology.reordered_atoms(&new_atom_order).unwrap().0;
    for (index, atom) in record.topology.atoms.iter_mut().enumerate() {
        let label = format!("a{index}");
        atom.set_prop("atomLabel", &label).unwrap();
    }

    let params = CxSmilesWriteParams {
        smiles: SmilesWriteParams::default(),
        fields: CxSmilesFields::ATOM_LABELS,
    };
    assert_eq!(
        write_cx_smiles_with_params(&record, &params).unwrap(),
        "NC(O)F.NC(O)F |$a1;a3;a6;a7;a2;a0;a4;a5$|"
    );
}

#[test]
fn multi_fragment_subset_maps_noncontiguous_atom_and_bond_rows() {
    // Pinned RDKit 2026.03.1, legacy perception, canonical/isomeric output,
    // cleanStereo=true, rootedAtAtom=-1, and CX_ATOM_LABELS plus
    // CX_COORDINATE_BONDS. Four components take the source subset path. The
    // CC=C component retains source atoms 2..=4 and bonds 1..=2, while atom
    // rows 0, 1, and 5..=8 and bond rows 0, 3, and 4 are outside it. Cleaning
    // the orphan ENDUPRIGHT on source bond 1 must not leave a slash there.
    let mut record = parse_smiles("CC.CC=C.CN->B.C", &Default::default()).unwrap();
    record.topology.bonds[1].set_direction(cosmolkit_types::BondDirection::EndUpRight);
    for (index, atom) in record.topology.atoms.iter_mut().enumerate() {
        atom.set_prop("atomLabel", &format!("a{index}")).unwrap();
    }
    assert_eq!(record.topology.atoms.len(), 9);
    assert_eq!(record.topology.bonds.len(), 5);
    assert_eq!(
        record.topology.bonds[1].direction(),
        cosmolkit_types::BondDirection::EndUpRight
    );
    assert_eq!(
        record.topology.bonds[1].stereo(),
        cosmolkit_types::BondStereo::None
    );

    let before = record.clone();
    let params = CxSmilesWriteParams {
        smiles: SmilesWriteParams::default(),
        fields: CxSmilesFields::ATOM_LABELS | CxSmilesFields::COORDINATE_BONDS,
    };
    assert_eq!(
        write_cx_smiles_with_params(&record, &params).unwrap(),
        "B[NH2]C.C.C=CC.CC |$a7;a6;a5;a8;a4;a3;a2;a0;a1$,C:1.0|"
    );
    assert_eq!(
        record, before,
        "subset preparation must preserve the source rows and direction"
    );
}

#[test]
fn writer_subset_remaps_forward_ring_stereo_references_after_renumbering() {
    // Pinned RDKit 2026.03.1, legacy perception, MolToCXSmiles with canonical
    // and isomeric output, cleanStereo=true, rootedAtAtom=-1, and CX_BOND_CFG.
    // The four-component input uses the BONDS_BETWEEN_ATOMS subset path: the
    // STEREOANY ring bond has both stereo references after its endpoint rows,
    // and one reference is row 12 while the selected ring has only 10 atoms.
    // The source writer output is C.C.C.C1=CCCCCCCCC1 |ctu:0|.
    let mut source = parse_smiles("C1CCCCC=CCCC1.C.C.C", &Default::default()).unwrap();
    source.topology.bonds[5].set_stereo_atoms(Some([
        cosmolkit_model::AtomId::new(4),
        cosmolkit_model::AtomId::new(7),
    ]));
    source.topology.bonds[5]
        .set_stereo(cosmolkit_types::BondStereo::Any)
        .unwrap();
    let record = renumbered_record(&source, &[5, 10, 6, 11, 4, 12, 0, 1, 2, 3, 8, 9, 7]);
    let bond = &record.topology.bonds[5];
    assert_eq!(bond.begin(), cosmolkit_model::AtomId::new(0));
    assert_eq!(bond.end(), cosmolkit_model::AtomId::new(2));
    assert_eq!(bond.stereo(), cosmolkit_types::BondStereo::Any);
    assert_eq!(
        bond.stereo_atoms(),
        Some([
            cosmolkit_model::AtomId::new(4),
            cosmolkit_model::AtomId::new(12),
        ])
    );

    let before = record.clone();
    let params = CxSmilesWriteParams {
        smiles: SmilesWriteParams::default(),
        fields: CxSmilesFields::BOND_CFG,
    };
    assert_eq!(
        write_cx_smiles_with_params(&record, &params).unwrap(),
        "C.C.C.C1=CCCCCCCCC1 |ctu:0|"
    );
    assert_eq!(record, before, "fragment preparation preserves the input");
}

#[test]
fn writer_subset_reports_removed_template_attachment_target() {
    // TemplateAttachmentOrder is additional COSMolKit model state with no
    // one-to-one RDKit atom field, so this is an explicit local typed-error
    // regression rather than an RDKit-parity output row. The writer must not
    // retain a stale reference when its selected subset excludes the target.
    let mut record = parse_smiles("CC.C.C.C", &Default::default()).unwrap();
    let order = cosmolkit_model::TemplateAttachmentOrder::new(vec![
        cosmolkit_model::TemplateAttachment::new(cosmolkit_model::AtomId::new(2), "outside"),
    ])
    .unwrap();
    record.topology.atoms[0] = cosmolkit_model::Atom::from_spec(
        cosmolkit_model::AtomId::new(0),
        cosmolkit_model::AtomSpec::new(cosmolkit_types::Element::C)
            .with_template_attachment_order(order),
    );
    let before = record.clone();

    let error = write_smiles_with_params(
        &record,
        &SmilesWriteParams {
            canonical: false,
            ..Default::default()
        },
    )
    .expect_err("a selected atom cannot retain a template target outside its fragment");
    assert!(matches!(
        &error,
        SmilesParseError::WriterStereo(detail)
            if detail.contains("template attachment entry 0 loses referenced atom 2")
    ));
    assert_eq!(record, before, "failed extraction preserves the input");
}

#[test]
fn tetrahedral_writer_options_preserve_source_order_and_hydrogen_tokens() {
    let input = "N[C@@H](C)C(=O)O";
    assert_eq!(
        write(input, &SmilesWriteParams::default()),
        "C[C@H](N)C(=O)O"
    );
    assert_eq!(
        write(
            input,
            &SmilesWriteParams {
                canonical: false,
                ..Default::default()
            }
        ),
        "N[C@@H](C)C(=O)O"
    );
    assert_eq!(
        write(
            input,
            &SmilesWriteParams {
                do_isomeric_smiles: false,
                ..Default::default()
            }
        ),
        "CC(N)C(=O)O"
    );
    assert_eq!(
        write(
            input,
            &SmilesWriteParams {
                all_bonds_explicit: true,
                ..Default::default()
            }
        ),
        "C-[C@H](-N)-C(=O)-O"
    );
    assert_eq!(
        write(
            input,
            &SmilesWriteParams {
                all_hydrogens_explicit: true,
                ..Default::default()
            }
        ),
        "[CH3][C@H]([NH2])[C](=[O])[OH]"
    );
}

#[test]
fn canonical_double_bond_directions_follow_ring_and_isomeric_options() {
    let input = r"C1=CC/C=C2C3=C/CC=CC=CC\3C\2C=C1";
    let record = parse_smiles(input, &Default::default()).expect("parse pinned case");
    let before = record.clone();
    assert_eq!(
        write_smiles(&record).unwrap(),
        r"C1=CC/C=C2\C3=C\CC=CC=CC3C2C=C1"
    );
    assert_eq!(
        record, before,
        "writer must leave the detached input unchanged"
    );
    assert_eq!(
        write_smiles_with_params(
            &record,
            &SmilesWriteParams {
                do_isomeric_smiles: false,
                ..Default::default()
            }
        )
        .unwrap(),
        "C1=CCC=C2C3=CCC=CC=CC3C2C=C1"
    );
    assert_eq!(
        write(
            "C/C=C/C",
            &SmilesWriteParams {
                do_isomeric_smiles: false,
                all_bonds_explicit: true,
                ..Default::default()
            }
        ),
        "C/C=C/C"
    );
}

#[test]
fn non_tetrahedral_writer_emits_sp_tb_and_oh_permutations() {
    for (input, isomeric, nonisomeric) in [
        (
            "[Pt@SP1](F)(Cl)(Br)I",
            "[F][Pt@SP1]([Cl])([Br])[I]",
            "[F][Pt]([Cl])([Br])[I]",
        ),
        (
            "[P@TB1](F)(Cl)(Br)(I)N",
            "N[P@TB8](F)(Cl)(Br)I",
            "NP(F)(Cl)(Br)I",
        ),
        (
            "[Co@OH1](F)(Cl)(Br)(I)(N)O",
            "[NH2][Co@OH9]([OH])([F])([Cl])([Br])[I]",
            "[NH2][Co]([OH])([F])([Cl])([Br])[I]",
        ),
    ] {
        assert_eq!(
            write(input, &SmilesWriteParams::default()),
            isomeric,
            "{input}"
        );
        assert_eq!(
            write(
                input,
                &SmilesWriteParams {
                    do_isomeric_smiles: false,
                    ..Default::default()
                }
            ),
            nonisomeric,
            "{input}"
        );
    }
}

#[test]
fn rooted_writer_starts_requested_component_at_atom_without_reordering_other_components() {
    let record = parse_smiles("CCCO", &Default::default()).unwrap();
    for (canonical, expected) in [(true, "C(O)CC"), (false, "C(CC)O")] {
        let params = SmilesWriteParams {
            canonical,
            rooted_at_atom: Some(cosmolkit_model::AtomId::new(2)),
            ..Default::default()
        };
        assert_eq!(
            write_smiles_with_params(&record, &params).unwrap(),
            expected
        );
    }

    let disconnected = parse_smiles("O.CCC", &Default::default()).unwrap();
    for (canonical, expected) in [(true, "C(C)C.O"), (false, "O.C(C)C")] {
        let params = SmilesWriteParams {
            canonical,
            rooted_at_atom: Some(cosmolkit_model::AtomId::new(2)),
            ..Default::default()
        };
        assert_eq!(
            write_smiles_with_params(&disconnected, &params).unwrap(),
            expected
        );
    }
}

#[test]
fn rooted_writer_matches_source_compact_index_for_interleaved_component_rows() {
    let parsed = parse_smiles("FCN.O", &Default::default()).unwrap();
    let record = renumbered_record(&parsed, &[0, 3, 1, 2]);
    let before = record.clone();

    for (root, expected) in [(0, "FCN.O"), (2, "NCF.O")] {
        let params = SmilesWriteParams {
            canonical: false,
            rooted_at_atom: Some(cosmolkit_model::AtomId::new(root)),
            ..Default::default()
        };
        assert_eq!(
            write_smiles_with_params(&record, &params).unwrap(),
            expected
        );
    }

    let error = write_smiles_with_params(
        &record,
        &SmilesWriteParams {
            canonical: false,
            rooted_at_atom: Some(cosmolkit_model::AtomId::new(3)),
            ..Default::default()
        },
    )
    .expect_err("pinned compact-fragment root index is outside the selected fragment");
    assert!(matches!(
        error,
        SmilesParseError::WriterRootAtomOutOfRange {
            atom_index: 3,
            atom_count: 3
        }
    ));
    assert_eq!(record, before);
}

#[test]
fn rooted_writer_preflight_error_precedes_kekulization_failure() {
    let mut record = parse_smiles("C", &Default::default()).unwrap();
    record.topology.atoms[0].set_aromatic(true);
    let before = record.clone();

    let error = write_smiles_with_params(
        &record,
        &SmilesWriteParams {
            do_kekule: true,
            rooted_at_atom: Some(cosmolkit_model::AtomId::new(1)),
            ..Default::default()
        },
    )
    .expect_err("source root precondition runs before source kekulization");
    assert!(matches!(
        error,
        SmilesParseError::WriterRootAtomOutOfRange {
            atom_index: 1,
            atom_count: 1
        }
    ));
    assert_eq!(record, before);
}

#[test]
fn rooted_writer_rejects_out_of_range_atom_but_keeps_empty_molecule_behavior() {
    let record = parse_smiles("CCCO", &Default::default()).unwrap();
    let error = write_smiles_with_params(
        &record,
        &SmilesWriteParams {
            rooted_at_atom: Some(cosmolkit_model::AtomId::new(4)),
            ..Default::default()
        },
    )
    .expect_err("rootedAtAtom must be in range");
    assert!(matches!(
        error,
        SmilesParseError::WriterRootAtomOutOfRange {
            atom_index: 4,
            atom_count: 4
        }
    ));

    let empty = parse_smiles("", &Default::default()).unwrap();
    assert_eq!(
        write_smiles_with_params(
            &empty,
            &SmilesWriteParams {
                rooted_at_atom: Some(cosmolkit_model::AtomId::new(0)),
                ..Default::default()
            }
        )
        .unwrap(),
        ""
    );
}

#[test]
fn ignored_atom_maps_change_canonical_ranking_but_remain_serialized() {
    let record = parse_smiles("C([CH3:99])([CH3:1])O", &Default::default()).unwrap();
    let source_atom_maps = record
        .topology
        .atoms
        .iter()
        .map(|atom| atom.atom_map())
        .collect::<Vec<_>>();
    assert_eq!(write_smiles(&record).unwrap(), "OC([CH3:1])[CH3:99]");
    assert_eq!(
        record
            .topology
            .atoms
            .iter()
            .map(|atom| atom.atom_map())
            .collect::<Vec<_>>(),
        source_atom_maps
    );
    let ignored = write_smiles_with_params(
        &record,
        &SmilesWriteParams {
            ignore_atom_map_numbers: true,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(ignored, "[CH3:99]C([CH3:1])O");
    assert!(ignored.contains(":99") && ignored.contains(":1"));
    assert_eq!(
        record
            .topology
            .atoms
            .iter()
            .map(|atom| atom.atom_map())
            .collect::<Vec<_>>(),
        source_atom_maps
    );

    let noncanonical = write_smiles_with_params(
        &record,
        &SmilesWriteParams {
            canonical: false,
            ignore_atom_map_numbers: true,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(noncanonical, "C(C)(C)O");
    assert_eq!(
        record
            .topology
            .atoms
            .iter()
            .map(|atom| atom.atom_map())
            .collect::<Vec<_>>(),
        source_atom_maps
    );

    let zero_map = parse_smiles("[CH3:0]O", &Default::default()).unwrap();
    assert_eq!(
        write_smiles_with_params(
            &zero_map,
            &SmilesWriteParams {
                ignore_atom_map_numbers: true,
                ..Default::default()
            }
        )
        .unwrap(),
        "CO"
    );
}

#[test]
fn nonisomeric_canonical_ranking_ignores_suppressed_isotope_and_is_renumbering_stable() {
    let record = parse_smiles("CC(C)C[13CH3]", &Default::default()).unwrap();
    let before = record.clone();
    let reversed = renumbered_record(&record, &[4, 3, 2, 1, 0]);

    for candidate in [&record, &reversed] {
        assert_eq!(
            write_smiles_with_params(
                candidate,
                &SmilesWriteParams {
                    do_isomeric_smiles: false,
                    ..Default::default()
                }
            )
            .unwrap(),
            "CCC(C)C"
        );
        assert_eq!(write_smiles(candidate).unwrap(), "CC(C)C[13CH3]");
    }
    assert_eq!(
        record, before,
        "canonical writing must leave input unchanged"
    );
}

#[test]
fn nonisomeric_canonical_ranking_ignores_suppressed_tetrahedral_stereo() {
    for (input, isomeric) in [
        ("C[C@H](N)C(C)O", "CC(O)[C@H](C)N"),
        ("C[C@@H](N)C(C)O", "CC(O)[C@@H](C)N"),
    ] {
        let record = parse_smiles(input, &Default::default()).unwrap();
        let before = record.clone();
        let reversed = renumbered_record(
            &record,
            &(0..record.topology.atoms.len()).rev().collect::<Vec<_>>(),
        );
        for candidate in [&record, &reversed] {
            assert_eq!(
                write_smiles_with_params(
                    candidate,
                    &SmilesWriteParams {
                        do_isomeric_smiles: false,
                        ..Default::default()
                    }
                )
                .unwrap(),
                "CC(N)C(C)O",
                "{input}"
            );
            assert_eq!(write_smiles(candidate).unwrap(), isomeric, "{input}");
        }
        assert_eq!(
            record, before,
            "canonical writing must leave input unchanged"
        );
    }
}

#[test]
fn nonisomeric_canonical_ranking_ignores_suppressed_double_bond_stereo() {
    for (input, isomeric) in [("C/C=C/C", "C/C=C/C"), (r"C/C=C\C", r"C/C=C\C")] {
        let record = parse_smiles(input, &Default::default()).unwrap();
        let before = record.clone();
        let reversed = renumbered_record(
            &record,
            &(0..record.topology.atoms.len()).rev().collect::<Vec<_>>(),
        );
        for candidate in [&record, &reversed] {
            assert_eq!(
                write_smiles_with_params(
                    candidate,
                    &SmilesWriteParams {
                        do_isomeric_smiles: false,
                        ..Default::default()
                    }
                )
                .unwrap(),
                "CC=CC",
                "{input}"
            );
            assert_eq!(write_smiles(candidate).unwrap(), isomeric, "{input}");
        }
        assert_eq!(
            record, before,
            "canonical writing must leave input unchanged"
        );
    }
}

#[test]
fn nonisomeric_canonical_fallback_keeps_pinned_cx_output_maps() {
    // Pinned RDKit 2026.03.1, legacy stereo perception, with
    // doIsomericSmiles=false, cleanStereo=false, doKekule=false, and
    // CX_ATOM_LABELS | CX_BOND_CFG. The chiral component is paired with a
    // second component in each source order so canonical and noncanonical
    // output maps exercise their distinct source ordering rules. The writer's
    // temporary fallback state is private; assert its observable output maps
    // and that source stereo state is left unchanged.
    for (input, noncanonical, canonical) in [
        (
            "C[C@H](N)C(C)O.CC",
            "CC(N)C(C)O.CC |$a0;a1;a2;a3;a4;a5;a6;a7$|",
            "CC.CC(N)C(C)O |$a6;a7;a0;a1;a2;a3;a4;a5$|",
        ),
        (
            "CC.C[C@H](N)C(C)O",
            "CC.CC(N)C(C)O |$a0;a1;a2;a3;a4;a5;a6;a7$|",
            "CC.CC(N)C(C)O |$a0;a1;a2;a3;a4;a5;a6;a7$|",
        ),
    ] {
        let mut record = parse_smiles(input, &Default::default()).unwrap();
        for (index, atom) in record.topology.atoms.iter_mut().enumerate() {
            atom.set_prop("atomLabel", &format!("a{index}")).unwrap();
        }
        let before = record.clone();

        for (canonical, expected) in [(false, noncanonical), (true, canonical)] {
            let params = CxSmilesWriteParams {
                smiles: SmilesWriteParams {
                    canonical,
                    do_isomeric_smiles: false,
                    clean_stereo: false,
                    do_kekule: false,
                    ..Default::default()
                },
                fields: CxSmilesFields::ATOM_LABELS | CxSmilesFields::BOND_CFG,
            };
            assert_eq!(
                write_cx_smiles_with_params(&record, &params).unwrap(),
                expected,
                "canonical={canonical}, input={input}"
            );
            assert_eq!(
                record, before,
                "writer must not mutate input stereo or source-order state"
            );
        }
    }
}

#[test]
fn ordinary_writer_removes_modeled_stereo_groups_before_canonical_ranking() {
    let input = "C[C@H](N)C[C@@H](N)C |o1:1,4|";
    let record = parse_smiles(input, &Default::default()).unwrap();
    assert_eq!(record.topology.stereo_groups.len(), 1);
    assert_eq!(
        record.topology.stereo_groups[0].kind(),
        cosmolkit_model::StereoGroupKind::Or
    );
    assert_eq!(
        record.topology.stereo_groups[0].atoms(),
        &[
            cosmolkit_model::AtomId::new(1),
            cosmolkit_model::AtomId::new(4)
        ]
    );
    let before = record.clone();
    let reversed = renumbered_record(
        &record,
        &(0..record.topology.atoms.len()).rev().collect::<Vec<_>>(),
    );
    let without_group = parse_smiles("C[C@H](N)C[C@@H](N)C", &Default::default()).unwrap();

    for (isomeric, expected) in [(true, "C[C@H](N)C[C@H](C)N"), (false, "CC(N)CC(C)N")] {
        let params = SmilesWriteParams {
            do_isomeric_smiles: isomeric,
            ..Default::default()
        };
        for candidate in [&record, &reversed] {
            assert_eq!(
                write_smiles_with_params(candidate, &params).unwrap(),
                expected
            );
        }
        assert_eq!(
            write_smiles_with_params(&without_group, &params).unwrap(),
            expected
        );
    }

    let cx = CxSmilesWriteParams {
        smiles: SmilesWriteParams::default(),
        fields: CxSmilesFields::ENHANCED_STEREO,
    };
    assert_eq!(
        write_cx_smiles_with_params(&record, &cx).unwrap(),
        "C[C@H](N)C[C@H](C)N |o1:1,4|"
    );
    assert_eq!(
        record, before,
        "writer must preserve the typed stereo group"
    );
}

#[test]
fn clean_stereo_default_cleans_invalid_tetrahedral_tag() {
    let mut record = parse_smiles("CC(C)(O)F", &Default::default()).unwrap();
    record.topology.atoms[1].set_chiral_tag(cosmolkit_types::ChiralTag::TetrahedralCw);
    record.properties.clear_prop("_StereochemDone");

    assert_eq!(write_smiles(&record).unwrap(), "CC(C)(O)F");
}

#[test]
fn clean_stereo_option_selects_source_cleanup_and_preservation() {
    let mut record = parse_smiles("CC(C)(O)F", &Default::default()).unwrap();
    record.topology.atoms[1].set_chiral_tag(cosmolkit_types::ChiralTag::TetrahedralCw);
    record.properties.clear_prop("_StereochemDone");

    assert_eq!(
        write_smiles_with_params(
            &record,
            &SmilesWriteParams {
                clean_stereo: false,
                ..Default::default()
            }
        )
        .unwrap(),
        "C[C@@](C)(O)F"
    );
    assert_eq!(
        write_smiles_with_params(
            &record,
            &SmilesWriteParams {
                clean_stereo: true,
                ..Default::default()
            }
        )
        .unwrap(),
        "CC(C)(O)F"
    );

    record
        .properties
        .set_computed_prop("_StereochemDone", "1")
        .unwrap();
    assert_eq!(write_smiles(&record).unwrap(), "C[C@@](C)(O)F");
}

#[test]
fn current_stereo_wrapper_candidates_match_pinned_writer_text() {
    let mut record =
        parse_smiles("F[C@H](Cl)C(Cl)(Br)C(F)(F)C/C=C/F", &Default::default()).unwrap();
    record.properties.clear_prop("_StereochemDone");
    let valence = cosmolkit_core::assign_valence_with_options_for_topology(
        &record.topology,
        cosmolkit_core::ValenceModel::RdkitLike,
        false,
    )
    .unwrap();
    let rings = cosmolkit_core::fast_find_rings_from_parts(
        record.topology.atoms.len(),
        &record.topology.bonds,
        &record.topology.adjacency,
    )
    .unwrap();

    let candidates = [
        (
            "existing false/false depiction wrapper",
            cosmolkit_core::assign_legacy_stereochemistry_for_depiction(
                record.topology.clone(),
                &valence,
                &rings,
            )
            .unwrap(),
        ),
        (
            "existing true/true cleanup wrapper",
            cosmolkit_core::assign_legacy_stereochemistry(
                record.topology.clone(),
                &valence,
                &rings,
            )
            .unwrap(),
        ),
    ];

    for (name, topology) in candidates {
        let mut prepared = record.clone();
        prepared.topology = topology;
        prepared
            .properties
            .set_computed_prop("_StereochemDone", "1")
            .unwrap();
        assert_eq!(
            write_smiles(&prepared).unwrap(),
            "F/C=C/CC(F)(F)C(Cl)(Br)[C@H](F)Cl",
            "{name}"
        );
    }
}

#[test]
fn current_stereo_wrapper_candidates_preserve_pinned_large_ring_cx_output() {
    // Pinned RDKit 2026.03.1 / revision
    // 351f8f378f8ad6bbd517980c38896e66bf907af8c, with explicit perception
    // profiles and identical parser-finished state. D3-D4 freezes these source
    // rows independently: legacy gives the first output for
    // (false,false), (true,false), and (true,true); modern gives the second
    // output for each pair. The different profile is observable in both the
    // parsed Z/Cis state and CX suffix, not attributable to either flag.
    const PINNED_LEGACY_CX: &str = "C1=C\\CCCCCCCC/1";
    const PINNED_MODERN_CX: &str = "C1=C\\CCCCCCCC/1 |c:0|";
    assert_ne!(PINNED_LEGACY_CX, PINNED_MODERN_CX);

    // COSMolKit exposes the actual writer clean_stereo option only. Its false
    // path uses the existing legacy false/false wrapper; its true path uses
    // the existing legacy true/true wrapper. There is no CK modern-profile or
    // true/false assignment entry, so neither source-only variant is called
    // or represented by a fabricated wrapper here.
    let parser = SmilesParseParams {
        sanitize: true,
        remove_hydrogens: false,
        ..Default::default()
    };
    let parsed = parse_smiles("C1CCCCC=CCCC1 |c:5|", &parser).unwrap();
    assert_eq!(parsed.properties.prop("_needsDetectBondStereo"), Some("1"));
    assert_eq!(parsed.properties.prop("_StereochemDone"), None);
    assert_eq!(
        parsed.topology.bonds[5].stereo(),
        cosmolkit_types::BondStereo::Cis
    );
    let post_chemistry = cosmolkit_smiles::SmilesRecord {
        topology: cosmolkit_core::sanitize_topology(
            &parsed.topology,
            &cosmolkit_core::SanitizeParams::default(),
        )
        .unwrap()
        .topology,
        coordinates: parsed.coordinates.clone(),
        properties: parsed.properties.clone(),
    };
    let finalized = finalize_smiles_stereo(post_chemistry, &parser).unwrap();
    assert_eq!(finalized.properties.prop("_needsDetectBondStereo"), None);
    assert_eq!(finalized.properties.prop("_StereochemDone"), None);
    assert_eq!(
        finalized.topology.bonds[4].direction(),
        cosmolkit_types::BondDirection::EndUpRight
    );
    assert_eq!(
        finalized.topology.bonds[5].stereo(),
        cosmolkit_types::BondStereo::Z
    );
    assert_eq!(
        finalized.topology.bonds[6].direction(),
        cosmolkit_types::BondDirection::EndDownRight
    );

    // Both available writer options start from this exact, real finalized
    // detached record. The source legacy output is checked for this topology;
    // the separately frozen modern output remains an explicit unavailable CK
    // profile row, not an alternative accepted result.
    let finalized_snapshot = finalized.clone();
    for clean_stereo in [false, true] {
        let output = write_cx_smiles_with_params(
            &finalized,
            &cosmolkit_smiles::CxSmilesWriteParams {
                smiles: SmilesWriteParams {
                    clean_stereo,
                    ..Default::default()
                },
                fields: cosmolkit_smiles::CxSmilesFields::BOND_CFG,
            },
        )
        .unwrap();
        assert_eq!(output, PINNED_LEGACY_CX, "clean_stereo={clean_stereo}");
    }
    assert_eq!(finalized, finalized_snapshot, "writer preserves its input");
}

#[test]
fn pending_cx_direction_phase_is_the_first_large_ring_writer_divergence() {
    // Pinned RDKit 2026.03.1, revision
    // 351f8f378f8ad6bbd517980c38896e66bf907af8c: use the explicit legacy
    // profile represented by the detached finalizer. For sanitize=true and
    // removeHs=false, MolFromSmiles runs chemistry before its pending CX
    // direction phase; MolToCXSmiles receives the resulting finalized state.
    // The modern raw-parser writer output is retained by the separate raw
    // boundary case; it is not an expectation for this legacy finalized path.
    let parser = SmilesParseParams {
        sanitize: true,
        remove_hydrogens: false,
        ..Default::default()
    };
    let parsed = parse_smiles("C1CCCCC=CCCC1 |c:5|", &parser).unwrap();
    let parsed_snapshot = parsed.clone();
    assert_eq!(parsed.properties.prop("_CXSMILES_Data"), Some("|c:5|"));
    assert_eq!(parsed.properties.prop("_needsDetectBondStereo"), Some("1"));
    assert_eq!(parsed.properties.prop("_StereochemDone"), None);
    assert!(parsed.coordinates.conformers_2d.is_empty());
    assert!(parsed.coordinates.conformers_3d.is_empty());
    assert!(
        parsed
            .topology
            .atoms
            .iter()
            .all(|atom| { atom.chiral_tag() == cosmolkit_types::ChiralTag::Unspecified })
    );
    assert_eq!(
        parsed.topology.bonds[4].direction(),
        cosmolkit_types::BondDirection::None
    );
    assert_eq!(
        parsed.topology.bonds[5].stereo(),
        cosmolkit_types::BondStereo::Cis
    );
    assert_eq!(
        parsed.topology.bonds[5].stereo_atoms(),
        Some([
            cosmolkit_model::AtomId::new(4),
            cosmolkit_model::AtomId::new(7)
        ])
    );
    assert_eq!(
        parsed.topology.bonds[6].direction(),
        cosmolkit_types::BondDirection::None
    );

    // Mirror the actual source order: chemistry first, then the one existing
    // canonical SMILES finalizer. Keep each stage snapshot for the boundary.
    let post_chemistry = cosmolkit_smiles::SmilesRecord {
        topology: cosmolkit_core::sanitize_topology(
            &parsed.topology,
            &cosmolkit_core::SanitizeParams::default(),
        )
        .unwrap()
        .topology,
        coordinates: parsed.coordinates.clone(),
        properties: parsed.properties.clone(),
    };
    assert_eq!(
        parsed, parsed_snapshot,
        "chemistry leaves parser output intact"
    );
    assert_eq!(
        post_chemistry.properties.prop("_needsDetectBondStereo"),
        Some("1")
    );
    assert_eq!(
        post_chemistry.topology.bonds[4].direction(),
        cosmolkit_types::BondDirection::None
    );
    assert_eq!(
        post_chemistry.topology.bonds[6].direction(),
        cosmolkit_types::BondDirection::None
    );
    let post_chemistry_snapshot = post_chemistry.clone();
    let finalized = finalize_smiles_stereo(post_chemistry, &parser).unwrap();
    assert_eq!(
        post_chemistry_snapshot
            .properties
            .prop("_needsDetectBondStereo"),
        Some("1"),
        "the pre-finalizer input snapshot retains pending parser work"
    );
    assert_eq!(finalized.properties.prop("_CXSMILES_Data"), Some("|c:5|"));
    assert_eq!(finalized.properties.prop("_needsDetectBondStereo"), None);
    assert_eq!(
        finalized.topology.bonds[4].direction(),
        cosmolkit_types::BondDirection::EndUpRight
    );
    assert_eq!(
        finalized.topology.bonds[5].stereo(),
        cosmolkit_types::BondStereo::Z
    );
    assert_eq!(
        finalized.topology.bonds[5].stereo_atoms(),
        Some([
            cosmolkit_model::AtomId::new(4),
            cosmolkit_model::AtomId::new(7)
        ])
    );
    assert_eq!(
        finalized.topology.bonds[6].direction(),
        cosmolkit_types::BondDirection::EndDownRight
    );

    // This is the exact MolToCXSmiles(CX_BOND_CFG) writer operation and the
    // default canonical/isomeric/cleanStereo=true writer profile. E4's
    // explicit legacy post-finalizer oracle row is C1=C\\CCCCCCCC/1.
    let finalized_snapshot = finalized.clone();
    let output = write_cx_smiles_with_params(
        &finalized,
        &CxSmilesWriteParams {
            smiles: SmilesWriteParams::default(),
            fields: CxSmilesFields::BOND_CFG,
        },
    )
    .unwrap();
    assert_eq!(output, "C1=C\\CCCCCCCC/1");
    assert_eq!(
        finalized, finalized_snapshot,
        "writing preserves finalized input"
    );
}

#[test]
fn pending_cx_stereo_text_is_preserved_when_clean_stereo_is_false() {
    // The detached CX parser leaves the source direction marker pending. This
    // freezes the clean=false output branch; it does not isolate the missing
    // direction phase, which remains covered by the failing case above.
    let record = parse_smiles("C1CCCCC=CCCC1 |c:5|", &Default::default()).unwrap();
    assert_eq!(record.properties.prop("_needsDetectBondStereo"), Some("1"));
    let expected = "C1=C\\CCCCCCCC/1 |c:0|";
    let params = CxSmilesWriteParams {
        smiles: SmilesWriteParams {
            clean_stereo: false,
            ..Default::default()
        },
        fields: CxSmilesFields::BOND_CFG,
    };

    assert_eq!(
        write_cx_smiles_with_params(&record, &params).unwrap(),
        expected
    );
}

#[test]
fn raw_pending_cx_writer_matches_the_unfinalized_source_boundary() {
    // Pinned source matrix: RDKit 2026.03.1, revision
    // 351f8f378f8ad6bbd517980c38896e66bf907af8c. The raw record has no
    // coordinates, sanitize=false, removeHs=false, a pending CX marker, Cis
    // on bond 5 with references [4, 7], and no directions on bonds 4/6.
    // Its state hash is
    // 9b2ea598e1f99968fbe3a00e39358564991c5efae58edbe424bf3dbc6de1775e.
    // Direct MolToCXSmiles(CX_BOND_CFG) with cleanStereo=true and
    // doIsomericSmiles=true returns `C1=CCCCCCCCC1` in the explicitly selected
    // legacy profile. The explicitly selected modern profile is recorded
    // separately and returns `C1=C\\CCCCCCCC/1 |c:0|`; this test follows the
    // fixed legacy source profile used by the detached writer core.
    let parser = SmilesParseParams {
        sanitize: false,
        remove_hydrogens: false,
        ..Default::default()
    };
    let record = parse_smiles("C1CCCCC=CCCC1 |c:5|", &parser).unwrap();
    assert_eq!(record.properties.prop("_needsDetectBondStereo"), Some("1"));
    assert_eq!(record.properties.prop("_StereochemDone"), None);
    assert_eq!(record.properties.prop("_CXSMILES_Data"), Some("|c:5|"));
    assert!(record.coordinates.conformers_2d.is_empty());
    assert!(record.coordinates.conformers_3d.is_empty());
    assert!(
        record
            .topology
            .atoms
            .iter()
            .all(|atom| atom.chiral_tag() == cosmolkit_types::ChiralTag::Unspecified)
    );
    assert_eq!(
        record.topology.bonds[4].direction(),
        cosmolkit_types::BondDirection::None
    );
    assert_eq!(
        record.topology.bonds[5].stereo(),
        cosmolkit_types::BondStereo::Cis
    );
    assert_eq!(
        record.topology.bonds[5].stereo_atoms(),
        Some([
            cosmolkit_model::AtomId::new(4),
            cosmolkit_model::AtomId::new(7)
        ])
    );
    assert_eq!(
        record.topology.bonds[6].direction(),
        cosmolkit_types::BondDirection::None
    );

    let before = record.clone();
    // MolToSmiles passes cleanStereo and defaults force and
    // flagPossibleStereoCenters to false. This record has no done marker, so
    // force=false does not skip assignment. Keep the writer-visible
    // canonical/isomeric options and CX field set fixed across both calls;
    // the detached writer has not yet wired the exact
    // clean=true/possible=false core entrypoint.
    let clean_false_params = CxSmilesWriteParams {
        smiles: SmilesWriteParams {
            do_isomeric_smiles: true,
            canonical: true,
            clean_stereo: false,
            ..Default::default()
        },
        fields: CxSmilesFields::BOND_CFG,
    };
    assert_eq!(
        write_cx_smiles_with_params(&record, &clean_false_params).unwrap(),
        "C1=C\\CCCCCCCC/1 |c:0|"
    );

    let params = CxSmilesWriteParams {
        smiles: SmilesWriteParams {
            do_isomeric_smiles: true,
            canonical: true,
            clean_stereo: true,
            ..Default::default()
        },
        fields: CxSmilesFields::BOND_CFG,
    };
    assert_eq!(
        write_cx_smiles_with_params(&record, &params).unwrap(),
        "C1=CCCCCCCCC1"
    );
    assert_eq!(
        record, before,
        "writing must leave the raw parser record unchanged"
    );
}

#[test]
fn cx_write_after_sanitize_and_smiles_finalization_uses_finalized_state() {
    // This is the pinned legacy profile used by the existing detached core
    // finalizer. The exact input, operation, and CX_BOND_CFG writer profile
    // match E4's post-finalizer source row for sanitize=true/removeHs=false.
    let parser = SmilesParseParams {
        sanitize: true,
        remove_hydrogens: false,
        ..Default::default()
    };
    let mut record = parse_smiles("C1CCCCC=CCCC1 |c:5|", &parser).unwrap();
    assert_eq!(record.properties.prop("_needsDetectBondStereo"), Some("1"));

    // Match Molecule::from_smiles_with_params: requested chemistry completes
    // before the single existing SMILES finalization stage.
    record.topology = cosmolkit_core::sanitize_topology(
        &record.topology,
        &cosmolkit_core::SanitizeParams::default(),
    )
    .unwrap()
    .topology;
    assert_eq!(record.properties.prop("_needsDetectBondStereo"), Some("1"));
    assert_eq!(
        record.topology.bonds[4].direction(),
        cosmolkit_types::BondDirection::None
    );
    assert_eq!(
        record.topology.bonds[6].direction(),
        cosmolkit_types::BondDirection::None
    );

    let finalized = finalize_smiles_stereo(record, &parser).unwrap();
    assert_eq!(finalized.properties.prop("_needsDetectBondStereo"), None);
    assert_eq!(
        finalized.topology.bonds[4].direction(),
        cosmolkit_types::BondDirection::EndUpRight
    );
    assert_eq!(
        finalized.topology.bonds[5].stereo(),
        cosmolkit_types::BondStereo::Z
    );
    assert_eq!(
        finalized.topology.bonds[6].direction(),
        cosmolkit_types::BondDirection::EndDownRight
    );

    let before = finalized.clone();
    let params = CxSmilesWriteParams {
        smiles: SmilesWriteParams::default(),
        fields: CxSmilesFields::BOND_CFG,
    };
    assert_eq!(
        write_cx_smiles_with_params(&finalized, &params).unwrap(),
        "C1=C\\CCCCCCCC/1"
    );
    assert_eq!(
        finalized, before,
        "writing must not finalize the record again"
    );
}

#[test]
fn both_false_smiles_finalization_retains_pending_marker_and_clears_directions() {
    let parser = SmilesParseParams {
        sanitize: false,
        remove_hydrogens: false,
        ..Default::default()
    };
    let raw = parse_smiles("C/C=C/C |c:1|", &parser).unwrap();
    assert_eq!(raw.properties.prop("_needsDetectBondStereo"), Some("1"));
    assert_eq!(
        raw.topology.bonds[0].direction(),
        cosmolkit_types::BondDirection::EndUpRight
    );
    assert_eq!(
        raw.topology.bonds[1].stereo(),
        cosmolkit_types::BondStereo::Cis
    );
    assert_eq!(
        raw.topology.bonds[2].direction(),
        cosmolkit_types::BondDirection::EndUpRight
    );

    let both_false = finalize_smiles_stereo(raw, &parser).unwrap();
    assert_eq!(
        both_false.properties.prop("_needsDetectBondStereo"),
        Some("1")
    );
    assert_eq!(
        both_false.topology.bonds[0].direction(),
        cosmolkit_types::BondDirection::EndUpRight
    );
    assert_eq!(
        both_false.topology.bonds[1].stereo(),
        cosmolkit_types::BondStereo::Cis
    );
    assert_eq!(
        both_false.topology.bonds[1].stereo_atoms(),
        Some([
            cosmolkit_model::AtomId::new(0),
            cosmolkit_model::AtomId::new(3)
        ])
    );
    assert_eq!(
        both_false.topology.bonds[2].direction(),
        cosmolkit_types::BondDirection::EndUpRight
    );
}

#[test]
fn both_false_finalization_clears_begin_wedge_and_dash_flags() {
    let parser = SmilesParseParams {
        sanitize: false,
        remove_hydrogens: false,
        ..Default::default()
    };
    let mut raw = parse_smiles("CCC", &parser).unwrap();
    raw.topology.bonds[0].set_direction(cosmolkit_types::BondDirection::BeginWedge);
    raw.topology.bonds[1].set_direction(cosmolkit_types::BondDirection::BeginDash);

    // Pinned Chirality.cpp::clearSingleBondDirFlags(mol, true) clears
    // BEGINWEDGE/BEGINDASH while retaining ENDDOWNRIGHT/ENDUPRIGHT slash flags.
    let finalized = finalize_smiles_stereo(raw, &parser).unwrap();
    assert_eq!(
        finalized.topology.bonds[0].direction(),
        cosmolkit_types::BondDirection::None
    );
    assert_eq!(
        finalized.topology.bonds[1].direction(),
        cosmolkit_types::BondDirection::None
    );
    assert_eq!(finalized.topology.bonds[0].prop("_UnknownStereo"), None);
    assert_eq!(finalized.topology.bonds[1].prop("_UnknownStereo"), None);
}

#[test]
fn both_false_finalization_clears_unknown_direction_and_records_unknown_stereo() {
    let parser = SmilesParseParams {
        sanitize: false,
        remove_hydrogens: false,
        ..Default::default()
    };
    let mut raw = parse_smiles("CC", &parser).unwrap();
    raw.topology.bonds[0].set_direction(cosmolkit_types::BondDirection::Unknown);

    // The source sets _UnknownStereo before clearing BondDir::UNKNOWN.
    let finalized = finalize_smiles_stereo(raw, &parser).unwrap();
    assert_eq!(
        finalized.topology.bonds[0].direction(),
        cosmolkit_types::BondDirection::None
    );
    assert_eq!(
        finalized.topology.bonds[0].prop("_UnknownStereo"),
        Some("1")
    );
}

#[test]
fn single_fragment_done_marker_observes_ordinary_and_computed_property_presence() {
    // Pinned MolToSmiles uses hasProp, so empty and false-like values remain
    // guards; the detached property model represents ordinary and computed
    // marker values as strings.
    for computed in [false, true] {
        for value in ["1", "", "0", "false"] {
            let mut record = parse_smiles("CC(C)(O)F", &Default::default()).unwrap();
            record.topology.atoms[1].set_chiral_tag(cosmolkit_types::ChiralTag::TetrahedralCw);
            record.properties.clear_prop("_StereochemDone");
            if computed {
                record
                    .properties
                    .set_computed_prop("_StereochemDone", value)
                    .unwrap();
            } else {
                record
                    .properties
                    .set_prop("_StereochemDone", value)
                    .unwrap();
            }
            assert_eq!(record.properties.prop("_StereochemDone"), Some(value));
            assert_eq!(
                record.properties.is_prop_computed("_StereochemDone"),
                computed
            );

            let before = record.clone();
            assert_eq!(
                write_smiles(&record).unwrap(),
                "C[C@@](C)(O)F",
                "presence must guard assignment (computed={computed}, value={value:?})"
            );
            assert_eq!(
                record, before,
                "writing must leave the marker and raw topology unchanged"
            );
        }
    }
}

#[test]
fn single_component_stereo_transport_keeps_source_identity_and_properties() {
    let mut record = parse_smiles("N[C@H](C)O", &Default::default()).unwrap();
    record.properties.clear_prop("_StereochemDone");
    for (index, atom) in record.topology.atoms.iter_mut().enumerate() {
        atom.set_atom_map(Some(10 + index as u32));
        atom.set_prop("source_atom", &format!("a{index}")).unwrap();
        atom.set_computed_prop("_sourceAtom", &format!("ca{index}"))
            .unwrap();
    }
    for (index, bond) in record.topology.bonds.iter_mut().enumerate() {
        bond.set_prop("source_bond", &format!("b{index}")).unwrap();
        bond.set_computed_prop("_sourceBond", &format!("cb{index}"))
            .unwrap();
    }

    let before = record.clone();
    assert_eq!(
        write_smiles(&record).unwrap(),
        "[NH2:10][C@H:11]([CH3:12])[OH:13]"
    );
    for (index, atom) in record.topology.atoms.iter().enumerate() {
        assert_eq!(atom.prop("source_atom"), Some(format!("a{index}").as_str()));
        assert_eq!(
            atom.prop("_sourceAtom"),
            Some(format!("ca{index}").as_str())
        );
        assert!(atom.is_prop_computed("_sourceAtom"));
    }
    for (index, bond) in record.topology.bonds.iter().enumerate() {
        assert_eq!(bond.prop("source_bond"), Some(format!("b{index}").as_str()));
        assert_eq!(
            bond.prop("_sourceBond"),
            Some(format!("cb{index}").as_str())
        );
        assert!(bond.is_prop_computed("_sourceBond"));
    }
    assert_eq!(
        record, before,
        "writer must preserve source state and identity"
    );
}

#[test]
fn multi_fragment_stereo_merge_leaves_unassigned_challenge_component_untouched() {
    // Pinned RDKit 2026.03.1, explicitly using legacy stereo perception and
    // canonical/isomeric output with cleanStereo=true. The leading ring is
    // challenging and keeps an ordinary _StereochemDone marker after
    // clone/prune; the later singleton fragments use subset extraction and
    // assignment. Their fragment-local atom zero must not overwrite the
    // first ring atom when assigned state is merged back.
    let mut record = parse_smiles(
        "[C@H]1(F)CC[C@H](F)CC1.C.C.C.C",
        &SmilesParseParams::default(),
    )
    .unwrap();
    record.properties.clear_prop("_StereochemDone");
    record.properties.set_prop("_StereochemDone", "1").unwrap();
    assert!(!record.properties.is_prop_computed("_StereochemDone"));

    let before = record.clone();
    assert_eq!(
        write_smiles(&record).unwrap(),
        "C.C.C.C.F[C@H]1CC[C@H](F)CC1"
    );
    assert_eq!(
        record, before,
        "fragment stereo merge must not mutate the source record"
    );
}

#[test]
fn multi_fragment_stereo_merge_remaps_signed_ring_ids_after_prefix() {
    // Pinned RDKit 2026.03.1, explicitly using legacy stereo perception and
    // canonical/isomeric output with cleanStereo=true. A computed done marker
    // is cleared on extracted fragments, so stereo assignment writes signed
    // ring-relative IDs in the ring fragment local index space. The leading
    // CC component makes that ring map begin at source atom two.
    let mut record = parse_smiles(
        "CC.[C@H]1(F)CC[C@H](F)CC1.C.C.C",
        &SmilesParseParams::default(),
    )
    .unwrap();
    record.properties.clear_prop("_StereochemDone");
    record
        .properties
        .set_computed_prop("_StereochemDone", "1")
        .unwrap();
    assert!(record.properties.is_prop_computed("_StereochemDone"));

    let before = record.clone();
    assert_eq!(
        write_smiles(&record).unwrap(),
        "C.C.C.CC.F[C@H]1CC[C@H](F)CC1"
    );
    assert_eq!(
        record, before,
        "fragment stereo merge must preserve the computed marker and source topology"
    );
}

#[test]
fn multi_fragment_challenge_clone_clears_computed_done_marker_on_prune() {
    // Pinned RDKit 2026.03.1 / legacy stereo perception. Five components
    // exercise both getTheFrags paths: the invalid chiral center is
    // challenging and is clone/pruned, while `CC` and the singleton atoms
    // take the computed-property-clearing subset path. RWMol::commitBatchEdit
    // clears computed molecule props after clone/prune, so the copied
    // `_StereochemDone` marker does not suppress cleanup on the challenging
    // component.
    let mut record = parse_smiles("CC(C)(O)F.CC.C.C.C", &Default::default()).unwrap();
    record.topology.atoms[1].set_chiral_tag(cosmolkit_types::ChiralTag::TetrahedralCw);
    record
        .properties
        .set_computed_prop("_StereochemDone", "1")
        .unwrap();
    assert_eq!(record.properties.prop("_StereochemDone"), Some("1"));
    assert!(record.properties.is_prop_computed("_StereochemDone"));

    let before = record.clone();
    assert_eq!(write_smiles(&record).unwrap(), "C.C.C.CC.CC(C)(O)F");
    assert_eq!(
        record, before,
        "writer fragment preparation must leave the input marker and invalid source tag unchanged"
    );
}

#[test]
fn multi_fragment_challenge_clone_prune_preserves_ordinary_done_marker() {
    // Pinned RDKit 2026.03.1 / legacy MolToSmiles profile. The ordinary
    // `_StereochemDone` property survives clone/prune, so the copied fragment
    // skips the source writer's stereo cleanup guard.
    let mut record = parse_smiles("CC(C)(O)F.CC.C.C.C", &Default::default()).unwrap();
    record.topology.atoms[1].set_chiral_tag(cosmolkit_types::ChiralTag::TetrahedralCw);
    record.properties.clear_prop("_StereochemDone");
    record.properties.set_prop("_StereochemDone", "1").unwrap();
    assert_eq!(record.properties.prop("_StereochemDone"), Some("1"));
    assert!(!record.properties.is_prop_computed("_StereochemDone"));

    let before = record.clone();
    assert_eq!(write_smiles(&record).unwrap(), "C.C.C.CC.C[C@@](C)(O)F");
    assert_eq!(
        record, before,
        "writer fragment preparation must preserve the ordinary input marker and source tag"
    );
}
