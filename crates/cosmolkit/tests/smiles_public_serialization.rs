#![cfg(feature = "cap-smiles")]

use cosmolkit as ck;
use cosmolkit_smiles as owner;

fn detached(molecule: &ck::Molecule) -> owner::SmilesRecord {
    owner::SmilesRecord {
        topology: molecule.topology().clone(),
        coordinates: ck::CoordinateBlock::default(),
        properties: molecule.properties().clone(),
    }
}

fn options(bits: u16) -> ck::SmilesWriteParams {
    ck::SmilesWriteParams {
        do_isomeric_smiles: bits & 1 != 0,
        do_kekule: bits & 2 != 0,
        canonical: bits & 4 != 0,
        clean_stereo: bits & 8 != 0,
        all_bonds_explicit: bits & 16 != 0,
        all_hydrogens_explicit: bits & 32 != 0,
        include_dative_bonds: bits & 64 != 0,
        ignore_atom_map_numbers: bits & 128 != 0,
        rooted_at_atom: None,
    }
}

fn compare(
    actual: Result<ck::PropertyText, ck::SmilesWriteError>,
    expected: Result<ck::PropertyText, owner::SmilesParseError>,
) {
    match (actual, expected) {
        (Ok(actual), Ok(expected)) => assert_eq!(actual, expected),
        (Err(ck::SmilesWriteError::Write(actual)), Err(expected)) => assert_eq!(actual, expected),
        (actual, expected) => panic!("public/owner mismatch: {actual:?} / {expected:?}"),
    }
}

#[test]
fn whole_writer_forwards_all_eight_boolean_options_and_roots() {
    // Adapter matrix, not an independent chemistry-parity claim. The owner
    // has fixed source regressions; this checks no public option is dropped.
    let mut cases = 0;
    for input in ["[CH3:7][C@H](F)c1ccccc1", "N->[Cu+2]"] {
        let molecule = ck::Molecule::from_smiles(input).unwrap();
        let record = detached(&molecule);
        for bits in 0..256 {
            for rooted_at_atom in [None, Some(ck::AtomId::new(0))] {
                let params = ck::SmilesWriteParams {
                    rooted_at_atom,
                    ..options(bits)
                };
                compare(
                    molecule.to_smiles_with_params(&params),
                    owner::write_smiles_with_params(&record, &params),
                );
                cases += 1;
            }
        }
        assert_eq!(molecule.topology(), &record.topology);
        assert_eq!(molecule.properties(), &record.properties);
    }
    assert_eq!(cases, 1024);
}

#[test]
fn default_strings_errors_and_registry_are_public_contracts() {
    // Fixed RDKit default outputs; no oracle or fixture generation in this test.
    for (input, expected) in [
        ("OCC", "CCO"),
        ("c1ccccc1", "c1ccccc1"),
        ("[13CH4]", "[13CH4]"),
        ("", ""),
    ] {
        let molecule = ck::Molecule::from_smiles(input).unwrap();
        assert_eq!(
            molecule.to_smiles().unwrap().as_bytes(),
            expected.as_bytes()
        );
        assert_eq!(
            molecule.to_smiles().unwrap(),
            molecule.to_smiles_with_params(&Default::default()).unwrap()
        );
        assert_eq!(
            molecule.to_cx_smiles().unwrap(),
            molecule
                .to_cx_smiles_with_params(&Default::default())
                .unwrap()
        );
    }
    let molecule = ck::Molecule::from_smiles("CCO").unwrap();
    let error = molecule
        .to_smiles_with_params(&ck::SmilesWriteParams {
            rooted_at_atom: Some(ck::AtomId::new(99)),
            ..Default::default()
        })
        .unwrap_err();
    assert!(matches!(
        error,
        ck::SmilesWriteError::Write(owner::SmilesParseError::WriterRootAtomOutOfRange {
            atom_index: 99,
            atom_count: 3
        })
    ));
    assert!(std::error::Error::source(&error).is_some());
    for name in [
        "to_smiles",
        "to_smiles_with_params",
        "to_cx_smiles",
        "to_cx_smiles_with_params",
        "to_fragment_smiles",
        "to_fragment_smiles_with_params",
        "to_fragment_cx_smiles",
        "to_fragment_cx_smiles_with_params",
        "to_random_smiles",
        "to_random_smiles_with_params",
    ] {
        let id = format!("Molecule.{name}");
        let entries: Vec<_> = ck::binding_contract::BINDING_CONTRACT
            .iter()
            .filter(|entry| entry.semantic_id == id)
            .collect();
        assert_eq!(entries.len(), 1, "{id}");
        assert_eq!(entries[0].feature, "cap-smiles");
        assert_eq!(
            entries[0].status,
            ck::binding_contract::FunctionStatus::Experimental
        );
    }
}

#[test]
fn fragment_option_matrix_preserves_original_ids_and_typed_failures() {
    let molecule = ck::Molecule::from_smiles("CCO").unwrap();
    let record = detached(&molecule);
    let atoms = vec![ck::AtomId::new(1), ck::AtomId::new(2)];
    assert_eq!(
        molecule.to_fragment_smiles(&atoms).unwrap(),
        ck::PropertyText::from("CO")
    );
    assert_eq!(
        molecule.to_fragment_cx_smiles(&atoms).unwrap(),
        ck::PropertyText::from("CO")
    );
    let mut cases = 0;
    for bits in 0..256 {
        for bonds in [None, Some(vec![]), Some(vec![ck::BondId::new(1)])] {
            for symbols in [None, Some(vec!["C".into(), "N".into(), "O".into()])] {
                let params = ck::FragmentSmilesWriteParams {
                    smiles: options(bits),
                    atoms: atoms.clone(),
                    bonds: bonds.clone(),
                    atom_symbols: symbols,
                    bond_symbols: None,
                };
                let expected = owner::write_fragment_smiles_output(
                    &record,
                    &params.smiles,
                    &params.atoms,
                    params.bonds.as_deref(),
                    params.atom_symbols.as_deref(),
                    None,
                    None,
                    None,
                )
                .map(|value| value.text);
                match (molecule.to_fragment_smiles_with_params(&params), expected) {
                    (Ok(actual), Ok(expected)) => assert_eq!(actual, expected),
                    (Err(ck::SmilesWriteError::Fragment(actual)), Err(expected)) => {
                        assert_eq!(actual, expected)
                    }
                    (actual, expected) => panic!("fragment mismatch {actual:?} / {expected:?}"),
                }
                cases += 1;
            }
        }
    }
    assert_eq!(cases, 1536);
    for atoms in [vec![], vec![ck::AtomId::new(99)]] {
        assert!(matches!(
            molecule.to_fragment_smiles(&atoms),
            Err(ck::SmilesWriteError::Fragment(_))
        ));
    }
    assert_eq!(molecule.topology(), &record.topology);
}

#[test]
fn random_vector_all_boolean_options_counts_and_seeds() {
    // One test owns this binary's random stream: comparisons reseed both
    // paths, preserving ordered duplicate outputs rather than comparing sets.
    let molecule = ck::Molecule::from_smiles("C[C@H](F)c1ccccc1").unwrap();
    let record = detached(&molecule);
    let mut cases = 0;
    for bits in 0..16 {
        let params = ck::RandomSmilesWriteParams {
            do_isomeric_smiles: bits & 1 != 0,
            do_kekule: bits & 2 != 0,
            all_bonds_explicit: bits & 4 != 0,
            all_hydrogens_explicit: bits & 8 != 0,
        };
        for count in [0, 1, 4] {
            for seed in [0, 42, 0x8000_0000] {
                // RDGeneral::getRandomGenerator takes i32 and reseeds only
                // for seed > 0. High-bit u32 seeds and zero continue state.
                // Establish identical prior state for both compared paths.
                owner::write_random_smiles_vector(&record, 0, 42, &params).unwrap();
                let expected =
                    owner::write_random_smiles_vector(&record, count, seed, &params).unwrap();
                molecule
                    .to_random_smiles_with_params(0, 42, &params)
                    .unwrap();
                assert_eq!(
                    molecule
                        .to_random_smiles_with_params(count, seed, &params)
                        .unwrap(),
                    expected
                );
                cases += 1;
            }
        }
    }
    assert_eq!(cases, 144);
    let methane = ck::Molecule::from_smiles("C").unwrap();
    assert_eq!(
        methane.to_random_smiles(4, 42).unwrap(),
        vec![ck::PropertyText::from("C"); 4]
    );
    assert_eq!(
        methane.to_random_smiles(4, 0).unwrap(),
        vec![ck::PropertyText::from("C"); 4]
    );
    assert_eq!(molecule.topology(), &record.topology);
}

#[test]
fn cx_field_selection_and_coordinate_ambiguity_matrix() {
    let base = ck::Molecule::from_smiles("C |$site$|").unwrap();
    let coordinates = ck::CoordinateBlock {
        conformers_2d: vec![ck::Conformer2D::new(7, vec![[1.0, 2.0]])],
        conformers_3d: vec![ck::Conformer3D::new(8, vec![[3.0, 4.0, 5.0]], true)],
        ..Default::default()
    };
    let record = owner::SmilesRecord {
        topology: base.topology().clone(),
        coordinates: coordinates.clone(),
        properties: base.properties().clone(),
    };
    let molecule = ck::Molecule::from_parts(
        record.topology.clone(),
        coordinates,
        record.properties.clone(),
    )
    .unwrap();
    let mut cases = 0;
    for bits in 0..8 {
        let mut fields = ck::CxSmilesFields::NONE;
        for (mask, field) in [
            (1, ck::CxSmilesFields::ATOM_LABELS),
            (2, ck::CxSmilesFields::ATOM_PROPS),
            (4, ck::CxSmilesFields::COORDS),
        ] {
            if bits & mask != 0 {
                fields = fields | field;
            }
        }
        for coordinate_selection in [
            ck::CxCoordinateSelection::Auto,
            ck::CxCoordinateSelection::TwoD { id: 7 },
            ck::CxCoordinateSelection::ThreeD { id: 8 },
            ck::CxCoordinateSelection::ThreeD { id: 99 },
        ] {
            let params = ck::CxSmilesWriteParams {
                fields,
                coordinate_selection,
                ..Default::default()
            };
            compare(
                molecule.to_cx_smiles_with_params(&params),
                owner::write_cx_smiles_with_params(&record, &params),
            );
            cases += 1;
        }
    }
    assert_eq!(cases, 32);
    assert!(molecule.to_cx_smiles().is_err());
    let params = ck::CxSmilesWriteParams {
        fields: ck::CxSmilesFields::COORDS,
        coordinate_selection: ck::CxCoordinateSelection::TwoD { id: 7 },
        ..Default::default()
    };
    assert_eq!(
        molecule.to_cx_smiles_with_params(&params).unwrap(),
        ck::PropertyText::from("C |(1,2,)|")
    );
    let fragment = ck::FragmentCxSmilesWriteParams {
        cx: params,
        atoms: vec![ck::AtomId::new(0)],
        ..Default::default()
    };
    assert_eq!(
        molecule
            .to_fragment_cx_smiles_with_params(&fragment)
            .unwrap(),
        ck::PropertyText::from("C |(1,2,)|")
    );
    assert_eq!(molecule.coordinates_2d(), Some(&[[1.0, 2.0]][..]));
    assert_eq!(molecule.conformers_3d(), &record.coordinates.conformers_3d);
    assert_eq!(molecule.topology(), &record.topology);
}
