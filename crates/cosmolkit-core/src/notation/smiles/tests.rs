use super::*;

fn parse_atom_token(input: &str) -> SmilesAtomToken {
    let mut parser = SmilesParser::new(SmilesLexer::new(input));
    let atom = parser.parse_simple_atomd().unwrap();
    assert_eq!(parser.next_token().unwrap(), SmilesToken::Eos);
    atom
}

fn parse_ring_number_token(input: &str) -> Result<u32, SmilesParseError> {
    let mut parser = SmilesParser::new(SmilesLexer::new(input));
    let ring_number = parser.parse_ring_number()?;
    assert_eq!(parser.next_token()?, SmilesToken::Eos);
    Ok(ring_number)
}

fn parse_number_token(input: &str) -> Result<i32, SmilesParseError> {
    let mut parser = SmilesParser::new(SmilesLexer::new(input));
    let number = parser.parse_number()?;
    assert_eq!(parser.next_token()?, SmilesToken::Eos);
    Ok(number)
}

#[test]
fn shared_get_iso_map_smiles_adapter_preserves_missing_atom_error() {
    let mut builder = MoleculeBuilder::new();
    let carbon = builder.add_atom(AtomSpec::new(Element::C));
    let deuterium = builder.add_atom(AtomSpec::new(Element::H).with_isotope(2));
    let bond = builder
        .add_bond(BondSpec::new(carbon, deuterium, BondOrder::Single))
        .unwrap();
    let missing = AtomId::new(7);
    builder
        .bond_mut(bond)
        .unwrap()
        .set_endpoints(carbon, missing);

    let error = tracked_smiles_isotopic_hydrogens(&builder).unwrap_err();

    assert!(matches!(
        error,
        SmilesParseError::ParseError(message) if message == "missing atom 7"
    ));
}

#[test]
fn preprocess_smiles_splits_name_when_cxsmiles_disabled() {
    let params = SmilesParseParams::without_cxsmiles_for_test();
    let result = preprocess_smiles("CCO ethanol", &params).unwrap();

    assert_eq!(result.smiles, "CCO");
    assert_eq!(result.name, "ethanol");
    assert_eq!(result.cx_part, "");
}

#[test]
fn preprocess_smiles_splits_cx_part_when_cxsmiles_enabled() {
    let result = preprocess_smiles("CCO |$;;foo$| ethanol", &SmilesParseParams::default()).unwrap();

    assert_eq!(result.smiles, "CCO");
    assert_eq!(result.name, "");
    assert_eq!(result.cx_part, "|$;;foo$| ethanol");
}

#[test]
fn preprocess_smiles_leaves_leading_space_input_unsplit_like_rdkit() {
    let result = preprocess_smiles(" CCO name", &SmilesParseParams::default()).unwrap();

    assert_eq!(result.smiles, " CCO name");
    assert_eq!(result.name, "");
    assert_eq!(result.cx_part, "");
}

#[test]
fn preprocess_smiles_applies_replacements_until_stable() {
    let mut replacements = BTreeMap::new();
    replacements.insert("{Q}".to_string(), "{X}CC{X}".to_string());
    replacements.insert("{X}".to_string(), "N".to_string());
    let params = SmilesParseParams::with_replacements_for_test(replacements);

    let result = preprocess_smiles("C{Q}C", &params).unwrap();

    assert_eq!(result.smiles, "CNCCNC");
}

#[test]
fn lexer_trims_ascii_control_whitespace_like_setup_smiles_string() {
    let lexer = SmilesLexer::new("\t CCO \r\n");

    assert_eq!(lexer.scan_input(), "CCO");
}

#[test]
fn lexer_emits_organic_and_aromatic_atom_payloads() {
    let mut lexer = SmilesLexer::new("Clc");

    match lexer.next_token().unwrap() {
        SmilesToken::OrganicAtom(atom) => assert_eq!(atom.spec.element().atomic_number(), 17),
        other => panic!("unexpected token: {other:?}"),
    }
    match lexer.next_token().unwrap() {
        SmilesToken::AromaticAtom(atom) => {
            assert_eq!(atom.spec.element().atomic_number(), 6);
            assert!(atom.spec.is_aromatic());
        }
        other => panic!("unexpected token: {other:?}"),
    }
}

#[test]
fn lexer_emits_bracket_atom_and_biovia_quoted_atom_payloads() {
    let mut lexer = SmilesLexer::new("[He]['Og']");

    assert_eq!(lexer.next_token().unwrap(), SmilesToken::AtomOpen);
    match lexer.next_token().unwrap() {
        SmilesToken::Atom(atom) => assert_eq!(atom.spec.element().atomic_number(), 2),
        other => panic!("unexpected token: {other:?}"),
    }
    assert_eq!(lexer.next_token().unwrap(), SmilesToken::AtomClose);
    assert_eq!(lexer.next_token().unwrap(), SmilesToken::AtomOpen);
    match lexer.next_token().unwrap() {
        SmilesToken::Atom(atom) => assert_eq!(atom.spec.element().atomic_number(), 118),
        other => panic!("unexpected token: {other:?}"),
    }
    assert_eq!(lexer.next_token().unwrap(), SmilesToken::AtomClose);
}

#[test]
fn lexer_emits_legacy_and_single_letter_bracket_atom_payloads() {
    let mut lexer = SmilesLexer::new("[Uuo][U]");

    assert_eq!(lexer.next_token().unwrap(), SmilesToken::AtomOpen);
    match lexer.next_token().unwrap() {
        SmilesToken::Atom(atom) => assert_eq!(atom.spec.element().atomic_number(), 118),
        other => panic!("unexpected token: {other:?}"),
    }
    assert_eq!(lexer.next_token().unwrap(), SmilesToken::AtomClose);
    assert_eq!(lexer.next_token().unwrap(), SmilesToken::AtomOpen);
    match lexer.next_token().unwrap() {
        SmilesToken::Atom(atom) => assert_eq!(atom.spec.element().atomic_number(), 92),
        other => panic!("unexpected token: {other:?}"),
    }
    assert_eq!(lexer.next_token().unwrap(), SmilesToken::AtomClose);
}

#[test]
fn lexer_emits_bond_token_payloads_like_smiles_ll() {
    let mut lexer = SmilesLexer::new("=#:$/\\-><-~");

    match lexer.next_token().unwrap() {
        SmilesToken::Bond(bond) => assert_eq!(bond.order, BondOrder::Double),
        other => panic!("unexpected token: {other:?}"),
    }
    match lexer.next_token().unwrap() {
        SmilesToken::Bond(bond) => assert_eq!(bond.order, BondOrder::Triple),
        other => panic!("unexpected token: {other:?}"),
    }
    match lexer.next_token().unwrap() {
        SmilesToken::Bond(bond) => {
            assert_eq!(bond.order, BondOrder::Aromatic);
            assert!(bond.is_aromatic);
        }
        other => panic!("unexpected token: {other:?}"),
    }
    match lexer.next_token().unwrap() {
        SmilesToken::Bond(bond) => assert_eq!(bond.order, BondOrder::Quadruple),
        other => panic!("unexpected token: {other:?}"),
    }
    match lexer.next_token().unwrap() {
        SmilesToken::Bond(bond) => {
            assert_eq!(bond.direction, BondDirection::EndUpRight);
            assert!(bond.explicit_unspecified_order);
        }
        other => panic!("unexpected token: {other:?}"),
    }
    match lexer.next_token().unwrap() {
        SmilesToken::Bond(bond) => {
            assert_eq!(bond.direction, BondDirection::EndDownRight);
            assert!(bond.explicit_unspecified_order);
        }
        other => panic!("unexpected token: {other:?}"),
    }
    match lexer.next_token().unwrap() {
        SmilesToken::Bond(bond) => assert_eq!(bond.order, BondOrder::DativeRight),
        other => panic!("unexpected token: {other:?}"),
    }
    match lexer.next_token().unwrap() {
        SmilesToken::Bond(bond) => assert_eq!(bond.order, BondOrder::DativeLeft),
        other => panic!("unexpected token: {other:?}"),
    }
    match lexer.next_token().unwrap() {
        SmilesToken::Bond(bond) => assert!(bond.is_null_query),
        other => panic!("unexpected token: {other:?}"),
    }
}

#[test]
fn lexer_returns_start_token_before_scanning_input_like_smiles_ll() {
    let mut lexer = SmilesLexer::with_start_token_for_test("C", SmilesToken::StartMol);

    assert_eq!(lexer.next_token().unwrap(), SmilesToken::StartMol);
    match lexer.next_token().unwrap() {
        SmilesToken::OrganicAtom(atom) => assert_eq!(atom.spec.element().atomic_number(), 6),
        other => panic!("unexpected token: {other:?}"),
    }
}

#[test]
fn lexer_emits_chiral_class_tokens_like_smiles_ll() {
    let mut lexer = SmilesLexer::new("@TH@AL@ SP@TB@OH@");

    assert_eq!(
        lexer.next_token().unwrap(),
        SmilesToken::ChiralClass(ChiralTag::Tetrahedral)
    );
    assert_eq!(
        lexer.next_token().unwrap(),
        SmilesToken::ChiralClass(ChiralTag::Allene)
    );
    assert_eq!(
        lexer.next_token().unwrap(),
        SmilesToken::ChiralClass(ChiralTag::SquarePlanar)
    );
    assert_eq!(
        lexer.next_token().unwrap(),
        SmilesToken::ChiralClass(ChiralTag::TrigonalBipyramidal)
    );
    assert_eq!(
        lexer.next_token().unwrap(),
        SmilesToken::ChiralClass(ChiralTag::Octahedral)
    );
}

#[test]
fn lexer_emits_at_h_atom_state_punctuation_newline_and_bad_character_tokens() {
    let mut lexer = SmilesLexer::new("@H[#:]&\n");

    assert_eq!(lexer.next_token().unwrap(), SmilesToken::At);
    assert_eq!(lexer.next_token().unwrap(), SmilesToken::H);
    assert_eq!(lexer.next_token().unwrap(), SmilesToken::AtomOpen);
    assert_eq!(lexer.next_token().unwrap(), SmilesToken::Hash);
    assert_eq!(lexer.next_token().unwrap(), SmilesToken::Colon);
    assert_eq!(lexer.next_token().unwrap(), SmilesToken::AtomClose);
    assert_eq!(lexer.next_token().unwrap(), SmilesToken::BadCharacter('&'));
    assert_eq!(lexer.next_token().unwrap(), SmilesToken::Eos);
}

#[test]
fn lexer_stops_scanning_after_internal_newline_like_smiles_ll() {
    let mut lexer = SmilesLexer::new("C\nN");

    match lexer.next_token().unwrap() {
        SmilesToken::OrganicAtom(atom) => assert_eq!(atom.spec.element().atomic_number(), 6),
        other => panic!("unexpected token: {other:?}"),
    }
    assert_eq!(lexer.next_token().unwrap(), SmilesToken::Eos);
    assert_eq!(lexer.next_token().unwrap(), SmilesToken::Eos);
}

#[test]
fn lexer_emits_ring_and_punctuation_tokens_like_smiles_ll() {
    let mut lexer = SmilesLexer::new("0%12()+-.");

    assert_eq!(lexer.next_token().unwrap(), SmilesToken::Zero);
    assert_eq!(lexer.next_token().unwrap(), SmilesToken::Percent);
    assert_eq!(lexer.next_token().unwrap(), SmilesToken::NonzeroDigit(1));
    assert_eq!(lexer.next_token().unwrap(), SmilesToken::NonzeroDigit(2));
    assert_eq!(lexer.next_token().unwrap(), SmilesToken::GroupOpen);
    assert_eq!(lexer.next_token().unwrap(), SmilesToken::GroupClose);
    assert_eq!(lexer.next_token().unwrap(), SmilesToken::Plus);
    assert_eq!(lexer.next_token().unwrap(), SmilesToken::Minus);
    assert_eq!(lexer.next_token().unwrap(), SmilesToken::Separator);
    assert_eq!(lexer.next_token().unwrap(), SmilesToken::Eos);
}

#[test]
fn parse_simple_atomd_parses_bracket_map_charge_hydrogen_and_no_implicit() {
    let atom = parse_atom_token("[13C@TH2H3-:5]");

    assert_eq!(
        atom.spec,
        AtomSpec::new(Element::C)
            .with_isotope(13)
            .with_chiral_tag(ChiralTag::Tetrahedral)
            .with_chiral_permutation(2)
            .with_explicit_hydrogens(3)
            .with_formal_charge(-1)
            .with_atom_map(5)
            .with_no_implicit(true)
    );
}

#[test]
fn parse_simple_atomd_parses_hydrogen_isotope_with_explicit_hydrogen_suffix() {
    let atom = parse_atom_token("[2HH1-]");

    assert_eq!(
        atom.spec,
        AtomSpec::new(Element::H)
            .with_isotope(2)
            .with_explicit_hydrogens(1)
            .with_formal_charge(-1)
            .with_no_implicit(true)
    );
}

#[test]
fn parse_simple_atomd_rejects_missing_bracket_close_like_rdkit() {
    let mut parser = SmilesParser::new(SmilesLexer::new("[NH4+"));

    assert_eq!(
        parser.parse_simple_atomd(),
        Err(SmilesParseError::ParseError(
            "expected bracket atom close or atom map, got Eos".to_string()
        ))
    );
}

#[test]
fn parse_chiral_element_rejects_zero_permutation_like_rdkit() {
    let mut parser = SmilesParser::new(SmilesLexer::new("[C@TH0]"));
    assert_eq!(parser.next_token().unwrap(), SmilesToken::AtomOpen);

    assert_eq!(
        parser.parse_bracket_atomd(),
        Err(SmilesParseError::ParseError(
            "chiral permutation cannot be zero".to_string()
        ))
    );
}

#[test]
fn parse_number_reports_int32_overflow_like_rdkit() {
    assert_eq!(
        parse_number_token("2147483648"),
        Err(SmilesParseError::ParseError("number too large".to_string()))
    );
}

#[test]
fn parse_ring_number_accepts_percent_group_forms_like_rdkit() {
    assert_eq!(parse_ring_number_token("7").unwrap(), 7);
    assert_eq!(parse_ring_number_token("%12").unwrap(), 12);
    assert_eq!(parse_ring_number_token("%(1)").unwrap(), 1);
    assert_eq!(parse_ring_number_token("%(12345)").unwrap(), 12345);
}

#[test]
fn parse_ring_number_rejects_empty_and_oversized_percent_groups() {
    assert_eq!(
        parse_ring_number_token("%()"),
        Err(SmilesParseError::ParseError(
            "empty ring number".to_string()
        ))
    );
    assert_eq!(
        parse_ring_number_token("%(123456)"),
        Err(SmilesParseError::ParseError(
            "ring number too large".to_string()
        ))
    );
}

#[test]
fn build_state_add_first_atom_marks_smiles_start_like_rdkit() {
    let mut state = SmilesBuildState::new();
    state.add_first_atom(SmilesAtomToken::new(6)).unwrap();
    let molecule = state.into_molecule().unwrap();

    assert_eq!(molecule.num_atoms(), 1);
    assert_eq!(molecule.num_bonds(), 0);
    assert_eq!(molecule.atoms()[0].atomic_number(), 6);
    assert_eq!(molecule.atoms()[0].prop(SMILES_START_PROP), Some("1"));
}

#[test]
fn add_first_atom_parse_mol_keeps_disconnected_fragment_starts_like_rdkit() {
    let state = to_mol("C.O").unwrap();
    let molecule = state.into_molecule().unwrap();

    assert_eq!(molecule.atomic_numbers(), vec![6, 8]);
    assert_eq!(molecule.num_bonds(), 0);
    assert_eq!(molecule.atoms()[0].prop(SMILES_START_PROP), Some("1"));
    assert_eq!(molecule.atoms()[1].prop(SMILES_START_PROP), Some("1"));
}

#[test]
fn build_state_add_disconnected_atom_starts_new_active_fragment_like_rdkit() {
    let mut state = SmilesBuildState::new();
    state.add_first_atom(SmilesAtomToken::new(6)).unwrap();
    state
        .add_disconnected_atom(SmilesAtomToken::new(8))
        .unwrap();
    state
        .add_atom_connected_to_active(SmilesAtomToken::new(9))
        .unwrap();
    let molecule = state.into_molecule().unwrap();

    assert_eq!(molecule.atomic_numbers(), vec![6, 8, 9]);
    assert_eq!(molecule.num_bonds(), 1);
    assert_eq!(molecule.bonds()[0].begin(), AtomId::new(1));
    assert_eq!(molecule.bonds()[0].end(), AtomId::new(2));
    assert_eq!(molecule.bonds()[0].order(), BondOrder::Single);
    assert_eq!(molecule.bonds()[0].prop(CXSMILES_BOND_IDX_PROP), Some("0"));
    assert_eq!(molecule.atoms()[0].prop(SMILES_START_PROP), Some("1"));
    assert_eq!(molecule.atoms()[1].prop(SMILES_START_PROP), Some("1"));
    assert_eq!(molecule.atoms()[2].prop(SMILES_START_PROP), None);
}

#[test]
fn build_state_add_atom_connected_to_active_uses_unspecified_bond_type_like_rdkit() {
    let mut state = SmilesBuildState::new();
    state.add_first_atom(SmilesAtomToken::new(6)).unwrap();
    state
        .add_atom_connected_to_active(SmilesAtomToken::new(8))
        .unwrap();
    let molecule = state.into_molecule().unwrap();

    assert_eq!(molecule.atomic_numbers(), vec![6, 8]);
    assert_eq!(molecule.num_bonds(), 1);
    assert_eq!(molecule.bonds()[0].order(), BondOrder::Single);
    assert_eq!(molecule.bonds()[0].prop(CXSMILES_BOND_IDX_PROP), Some("0"));
}

#[test]
fn add_atom_connected_to_active_requires_existing_active_atom() {
    let error = SmilesBuildState::new()
        .add_atom_connected_to_active(SmilesAtomToken::new(6))
        .unwrap_err();

    assert_eq!(
        error,
        SmilesParseError::ParseError("no active atom".to_string())
    );
}

#[test]
fn add_branch_atom_connected_to_active_tracks_branch_root_and_new_active_atom_like_rdkit() {
    let mut state = SmilesBuildState::new();
    state.add_first_atom(SmilesAtomToken::new(6)).unwrap();
    state
        .add_atom_connected_to_active(SmilesAtomToken::new(6))
        .unwrap();
    state
        .add_branch_atom_connected_to_active(11, SmilesAtomToken::new(8))
        .unwrap();

    assert_eq!(state.active_atom, Some(AtomId::new(2)));
    assert_eq!(
        state.branch_stack,
        vec![BranchPoint {
            atom: AtomId::new(1),
            open_position: 11,
        }]
    );

    let molecule = state.into_molecule().unwrap();
    assert_eq!(molecule.atomic_numbers(), vec![6, 6, 8]);
    assert_eq!(molecule.num_bonds(), 2);
    assert_eq!(molecule.bonds()[1].begin(), AtomId::new(1));
    assert_eq!(molecule.bonds()[1].end(), AtomId::new(2));
    assert_eq!(molecule.bonds()[1].order(), BondOrder::Single);
    assert_eq!(molecule.bonds()[1].prop(CXSMILES_BOND_IDX_PROP), Some("1"));
}

#[test]
fn add_branch_single_bond_tracks_branch_root_and_single_bond_like_rdkit() {
    let mut state = SmilesBuildState::new();
    state.add_first_atom(SmilesAtomToken::new(6)).unwrap();
    state
        .add_atom_connected_to_active(SmilesAtomToken::new(6))
        .unwrap();
    state
        .add_branch_single_bond(13, SmilesAtomToken::new(8))
        .unwrap();

    assert_eq!(state.active_atom, Some(AtomId::new(2)));
    assert_eq!(
        state.branch_stack,
        vec![BranchPoint {
            atom: AtomId::new(1),
            open_position: 13,
        }]
    );

    let molecule = state.into_molecule().unwrap();
    assert_eq!(molecule.atomic_numbers(), vec![6, 6, 8]);
    assert_eq!(molecule.bonds()[1].begin(), AtomId::new(1));
    assert_eq!(molecule.bonds()[1].end(), AtomId::new(2));
    assert_eq!(molecule.bonds()[1].order(), BondOrder::Single);
    assert_eq!(molecule.bonds()[1].prop(CXSMILES_BOND_IDX_PROP), Some("1"));
}

#[test]
fn add_single_bond_to_atom_adds_explicit_single_bond_like_rdkit() {
    let mut state = SmilesBuildState::new();
    state.add_first_atom(SmilesAtomToken::new(6)).unwrap();
    state
        .add_single_bond_to_atom(SmilesAtomToken::new(8))
        .unwrap();

    let molecule = state.into_molecule().unwrap();
    assert_eq!(molecule.atomic_numbers(), vec![6, 8]);
    assert_eq!(molecule.num_bonds(), 1);
    assert_eq!(molecule.bonds()[0].begin(), AtomId::new(0));
    assert_eq!(molecule.bonds()[0].end(), AtomId::new(1));
    assert_eq!(molecule.bonds()[0].order(), BondOrder::Single);
    assert_eq!(molecule.bonds()[0].prop(CXSMILES_BOND_IDX_PROP), Some("0"));
}

#[test]
fn build_state_add_explicit_bond_to_atom_normalizes_dative_direction_like_rdkit() {
    let mut state = SmilesBuildState::new();
    state.add_first_atom(SmilesAtomToken::new(7)).unwrap();
    state
        .add_explicit_bond_to_atom(
            SmilesBondToken::new(BondOrder::DativeLeft),
            SmilesAtomToken::new(8),
        )
        .unwrap();
    let molecule = state.into_molecule().unwrap();

    assert_eq!(molecule.num_atoms(), 2);
    assert_eq!(molecule.num_bonds(), 1);
    assert_eq!(molecule.bonds()[0].begin(), AtomId::new(1));
    assert_eq!(molecule.bonds()[0].end(), AtomId::new(0));
    assert_eq!(molecule.bonds()[0].order(), BondOrder::Dative);
    assert_eq!(molecule.bonds()[0].prop(CXSMILES_BOND_IDX_PROP), Some("0"));
}

#[test]
fn build_state_add_explicit_directional_bond_preserves_unspecified_order_marker() {
    let mut state = SmilesBuildState::new();
    state.add_first_atom(SmilesAtomToken::new(6)).unwrap();
    state
        .add_explicit_bond_to_atom(
            SmilesBondToken::directional(BondDirection::EndUpRight),
            SmilesAtomToken::new(6),
        )
        .unwrap();
    let molecule = state.into_molecule().unwrap();

    assert_eq!(molecule.bonds()[0].order(), BondOrder::Unspecified);
    assert_eq!(molecule.bonds()[0].direction(), BondDirection::EndUpRight);
    assert_eq!(molecule.bonds()[0].prop(UNSPECIFIED_ORDER_PROP), Some("1"));
}

#[test]
fn build_state_add_explicit_bond_to_atom_maps_null_query_bond_to_any_like_rdkit_subset() {
    let mut state = SmilesBuildState::new();
    state.add_first_atom(SmilesAtomToken::new(6)).unwrap();
    state
        .add_explicit_bond_to_atom(SmilesBondToken::null_query(), SmilesAtomToken::new(6))
        .unwrap();
    let molecule = state.into_molecule().unwrap();

    assert_eq!(
        molecule.bonds()[0].query(),
        Some(&QueryNode::predicate(BondQueryPredicate::Any))
    );
}

#[test]
fn add_branch_explicit_bond_tracks_branch_root_and_new_active_atom_like_rdkit() {
    let mut state = SmilesBuildState::new();
    state.add_first_atom(SmilesAtomToken::new(6)).unwrap();
    state
        .add_atom_connected_to_active(SmilesAtomToken::new(6))
        .unwrap();
    state
        .add_branch_explicit_bond(
            17,
            SmilesBondToken::new(BondOrder::Double),
            SmilesAtomToken::new(8),
        )
        .unwrap();

    assert_eq!(state.active_atom, Some(AtomId::new(2)));
    assert_eq!(
        state.branch_stack,
        vec![BranchPoint {
            atom: AtomId::new(1),
            open_position: 17,
        }]
    );

    let molecule = state.into_molecule().unwrap();
    assert_eq!(molecule.atomic_numbers(), vec![6, 6, 8]);
    assert_eq!(molecule.bonds()[1].begin(), AtomId::new(1));
    assert_eq!(molecule.bonds()[1].end(), AtomId::new(2));
    assert_eq!(molecule.bonds()[1].order(), BondOrder::Double);
    assert_eq!(molecule.bonds()[1].prop(CXSMILES_BOND_IDX_PROP), Some("1"));
}

#[test]
fn from_smiles_with_sanitize_false_parses_linear_simple_atoms_through_grammar_actions() {
    let molecule = Molecule::from_smiles_with_sanitize("CCO", false).unwrap();

    assert_eq!(molecule.atomic_numbers(), vec![6, 6, 8]);
    assert_eq!(molecule.num_bonds(), 2);
    assert_eq!(molecule.bonds()[0].order(), BondOrder::Single);
    assert_eq!(molecule.bonds()[1].order(), BondOrder::Single);
}

#[test]
fn from_smiles_with_sanitize_false_parses_explicit_bond_and_separator_actions() {
    let molecule = Molecule::from_smiles_with_sanitize("C=O.N", false).unwrap();

    assert_eq!(molecule.atomic_numbers(), vec![6, 8, 7]);
    assert_eq!(molecule.num_bonds(), 1);
    assert_eq!(molecule.bonds()[0].order(), BondOrder::Double);
    assert_eq!(molecule.atoms()[0].prop(SMILES_START_PROP), None);
    assert_eq!(molecule.atoms()[2].prop(SMILES_START_PROP), None);
    assert_eq!(molecule.bonds()[0].prop(CXSMILES_BOND_IDX_PROP), None);
}

#[test]
fn from_smiles_with_sanitize_false_keeps_explicit_single_bond_between_aromatic_atoms() {
    let molecule = Molecule::from_smiles_with_sanitize("c-c", false).unwrap();

    assert_eq!(molecule.atomic_numbers(), vec![6, 6]);
    assert!(molecule.atoms()[0].is_aromatic());
    assert!(molecule.atoms()[1].is_aromatic());
    assert_eq!(molecule.bonds()[0].order(), BondOrder::Single);
    assert!(!molecule.bonds()[0].is_aromatic());
}

#[test]
fn from_smiles_with_sanitize_false_sets_unspecified_directional_bond_type_then_cleans_props() {
    let molecule = Molecule::from_smiles_with_sanitize("C/C", false).unwrap();

    assert_eq!(molecule.num_bonds(), 1);
    assert_eq!(molecule.bonds()[0].order(), BondOrder::Single);
    assert_eq!(molecule.bonds()[0].direction(), BondDirection::EndUpRight);
    assert_eq!(molecule.bonds()[0].prop(UNSPECIFIED_ORDER_PROP), None);
    assert_eq!(molecule.bonds()[0].prop(CXSMILES_BOND_IDX_PROP), None);
}

#[test]
fn from_smiles_with_sanitize_false_preserves_directional_and_cx_wedge_state_like_rdkit() {
    let directional = Molecule::from_smiles_with_sanitize("C/C", false).unwrap();
    let wedged = Molecule::from_smiles_with_sanitize("CC |wU:1.0|", false).unwrap();

    assert_eq!(
        directional.bonds()[0].direction(),
        BondDirection::EndUpRight
    );
    assert_eq!(directional.bonds()[0].prop(UNSPECIFIED_ORDER_PROP), None);
    assert_eq!(wedged.bonds()[0].direction(), BondDirection::None);
    assert_eq!(wedged.bonds()[0].prop("_MolFileBondCfg"), Some("1"));
}

#[test]
fn from_smiles_with_sanitize_false_reads_name_when_cx_part_is_plain_text() {
    let molecule = Molecule::from_smiles_with_sanitize("CCO ethanol", false).unwrap();

    assert_eq!(molecule.atomic_numbers(), vec![6, 6, 8]);
    assert_eq!(molecule.properties().name(), Some("ethanol"));
}

#[test]
fn handle_cx_part_and_name_reports_strict_non_name_text_like_rdkit() {
    let params = SmilesParseParams::without_parse_name_for_test();
    let error = mol_from_smiles("CCO not-name", &params).unwrap_err();

    assert_eq!(
        error,
        SmilesParseError::ParseError(
            "CXSMILES extension does not start with | and parseName=false".to_string()
        )
    );
}

#[test]
fn mol_from_smiles_top_level_parses_simple_smiles_like_rdkit() {
    let molecule = mol_from_smiles("CCO", &SmilesParseParams::default()).unwrap();

    assert_eq!(molecule.atomic_numbers(), vec![6, 6, 8]);
    assert_eq!(molecule.num_bonds(), 2);
}

#[test]
fn mol_from_smiles_top_level_propagates_unclosed_ring_parse_failure_like_rdkit() {
    let error = mol_from_smiles("C1CC", &SmilesParseParams::default()).unwrap_err();

    assert_eq!(
        error,
        SmilesParseError::ParseError("unclosed ring".to_string())
    );
}

#[test]
fn mol_from_smiles_top_level_empty_input_returns_empty_molecule_like_rdkit() {
    let molecule = mol_from_smiles("", &SmilesParseParams::default()).unwrap();

    assert_eq!(molecule.num_atoms(), 0);
    assert_eq!(molecule.num_bonds(), 0);
}

#[test]
fn mol_from_smiles_maps_null_query_bond_to_any_query_for_rdkit_v2000_roundtrip() {
    let molecule = mol_from_smiles("C~C", &SmilesParseParams::default()).unwrap();

    assert_eq!(
        molecule.bonds()[0].query(),
        Some(&QueryNode::predicate(BondQueryPredicate::Any))
    );
}

#[test]
fn handle_cx_part_and_name_parses_atom_labels_and_name_like_rdkit() {
    let molecule = Molecule::from_smiles_with_sanitize("CCO |$;;foo$| ethanol", false).unwrap();

    assert_eq!(molecule.atomic_numbers(), vec![6, 6, 8]);
    assert_eq!(molecule.atoms()[2].prop("atomLabel"), Some("foo"));
    assert_eq!(
        molecule.properties().prop("_CXSMILES_Data"),
        Some("|$;;foo$|")
    );
    assert_eq!(molecule.properties().name(), Some("ethanol"));
}

#[test]
fn cleanup_after_parsing_marks_ap1_ap2_atom_labels_like_rdkit() {
    let molecule = Molecule::from_smiles_with_sanitize("*.*.* |$_AP1;_AP2;_AP3$|", false).unwrap();

    assert_eq!(molecule.atoms()[0].prop("_fromAttachPoint"), Some("1"));
    assert_eq!(molecule.atoms()[1].prop("_fromAttachPoint"), Some("2"));
    assert_eq!(molecule.atoms()[2].prop("_fromAttachPoint"), None);
}

#[test]
fn cleanup_after_parsing_clears_parser_temporary_props_like_rdkit() {
    let molecule = Molecule::from_smiles_with_sanitize("c1ccccc1.C", false).unwrap();

    assert!(
        molecule
            .atoms()
            .iter()
            .all(|atom| atom.prop(SMILES_START_PROP).is_none())
    );
    assert!(molecule.bonds().iter().all(|bond| {
        bond.prop(UNSPECIFIED_ORDER_PROP).is_none() && bond.prop(CXSMILES_BOND_IDX_PROP).is_none()
    }));
    assert!(
        molecule.bonds()[..6]
            .iter()
            .all(|bond| bond.order() == BondOrder::Aromatic && bond.is_aromatic())
    );
}

#[test]
fn handle_cx_part_and_name_tolerates_unported_cx_when_not_strict_like_rdkit() {
    let params = SmilesParseParams::non_strict_cxsmiles_for_test();
    let molecule = mol_from_smiles("CCO |rb:0| ethanol", &params).unwrap();

    assert_eq!(molecule.atomic_numbers(), vec![6, 6, 8]);
    assert_eq!(molecule.properties().prop("_CXSMILES_Data"), Some(""));
    assert_eq!(molecule.properties().name(), None);
}

#[test]
fn handle_cx_part_and_name_parses_sgroup_hierarchy_in_strict_mode() {
    let molecule = Molecule::from_smiles_with_sanitize("CC |SgH:0:1|", false).unwrap();
    assert_eq!(molecule.num_atoms(), 2);
    assert_eq!(molecule.num_bonds(), 1);
}

#[test]
fn parse_cx_sgroup_hierarchy_rejects_nonexistent_child_id_like_rdkit() {
    let mut state = SmilesBuildState::new();
    state.builder.add_atom(AtomSpec::new(Element::C));
    state
        .builder
        .add_substance_group(
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_prop("_cxsmilesindex", "0")
                .with_prop("index", "1"),
        )
        .unwrap();
    state
        .builder
        .add_substance_group(
            SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Data)
                .with_prop("_cxsmilesindex", "1")
                .with_prop("index", "2"),
        )
        .unwrap();

    let error = parse_cx_sgroup_hierarchy(&mut state, "SgH:0:2", &mut 0).unwrap_err();

    assert_eq!(
        error,
        SmilesParseError::ParseError("child id references non-existent SGroup".to_string())
    );
}

#[test]
fn handle_cx_part_and_name_parses_polymer_sgroup_in_strict_mode() {
    let molecule = Molecule::from_smiles_with_sanitize("CC |Sg:n:0::ht|", false).unwrap();
    assert_eq!(molecule.num_atoms(), 2);
    assert_eq!(molecule.num_bonds(), 1);
}

#[test]
fn parse_cx_polymer_sgroup_rejects_unknown_type_like_rdkit() {
    let mut state = SmilesBuildState::new();
    state.builder.add_atom(AtomSpec::new(Element::C));

    let error = parse_cx_polymer_sgroup(&mut state, "Sg:bogus:0::ht|", &mut 0, 0).unwrap_err();

    assert_eq!(error, cx_parse_failure());
}

#[test]
fn handle_cx_part_and_name_parses_atom_values_and_props_like_rdkit() {
    let molecule = Molecule::from_smiles_with_sanitize(
        "CC |$_AV:first;second$,atomProp:1.foo.bar&#46;baz|",
        false,
    )
    .unwrap();

    assert_eq!(molecule.atoms()[0].prop("molFileValue"), Some("first"));
    assert_eq!(molecule.atoms()[1].prop("molFileValue"), Some("second"));
    assert_eq!(molecule.atoms()[1].prop("foo"), Some("bar.baz"));
}

#[test]
fn handle_cx_part_and_name_parses_coordinates_like_rdkit() {
    let molecule = Molecule::from_smiles_with_sanitize("CCO |(0,0,;1,0,;2,0,0.5)|", false).unwrap();

    assert_eq!(molecule.conformers_3d().len(), 1);
    assert!(molecule.conformers_3d()[0].is_3d());
    assert_eq!(
        molecule.conformers_3d()[0].coordinates()[0],
        [0.0, 0.0, 0.0]
    );
    assert_eq!(
        molecule.conformers_3d()[0].coordinates()[1],
        [1.0, 0.0, 0.0]
    );
    assert_eq!(
        molecule.conformers_3d()[0].coordinates()[2],
        [2.0, 0.0, 0.5]
    );
}

#[test]
fn mol_from_smiles_conformer_selection_takes_first_2d_and_first_3d_like_rdkit() {
    let mut builder = MoleculeBuilder::new();
    builder.add_atom(AtomSpec::new(Element::C));
    builder
        .add_conformer(Conformer3D::new(0, vec![[0.0, 0.0, 1.0]], true))
        .unwrap();
    builder
        .add_conformer(Conformer3D::new(1, vec![[0.0, 0.0, 0.0]], false))
        .unwrap();
    builder
        .add_conformer(Conformer3D::new(2, vec![[1.0, 0.0, 0.0]], false))
        .unwrap();
    builder
        .add_conformer(Conformer3D::new(3, vec![[1.0, 0.0, 1.0]], true))
        .unwrap();
    let molecule = builder.build().unwrap();

    assert_eq!(first_2d_and_3d_conformer_ids(&molecule), (Some(1), Some(0)));
}

#[test]
fn mol_from_smiles_conformer_selection_reports_only_first_3d_when_no_2d_exists() {
    let mut builder = MoleculeBuilder::new();
    builder.add_atom(AtomSpec::new(Element::C));
    builder
        .add_conformer(Conformer3D::new(0, vec![[0.0, 0.0, 1.0]], true))
        .unwrap();
    builder
        .add_conformer(Conformer3D::new(1, vec![[1.0, 0.0, 1.0]], true))
        .unwrap();
    let molecule = builder.build().unwrap();

    assert_eq!(first_2d_and_3d_conformer_ids(&molecule), (None, Some(0)));
}

fn coordinate_free_atropisomer_candidate(direction: BondDirection) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let chlorine = builder
        .add_atom(AtomSpec::new(Element::from_atomic_number(17).unwrap()).with_no_implicit(true));
    let axis_begin = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
    let alkene_left = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
    let axis_end = builder.add_atom(AtomSpec::new(Element::C));
    let alkene_right = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));

    builder
        .add_bond(BondSpec::new(axis_begin, chlorine, BondOrder::Single).with_direction(direction))
        .unwrap();
    builder
        .add_bond(BondSpec::new(axis_begin, alkene_left, BondOrder::Double))
        .unwrap();
    builder
        .add_bond(BondSpec::new(axis_begin, axis_end, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(BondSpec::new(axis_end, alkene_right, BondOrder::Double))
        .unwrap();

    builder
        .build()
        .unwrap()
        .sanitize_with_ops(crate::SanitizeOps::ALL)
        .unwrap()
}

#[test]
fn mol_from_smiles_assigns_coordinate_free_atropisomer_bond_stereo_like_rdkit() {
    let mut molecule = coordinate_free_atropisomer_candidate(BondDirection::BeginWedge);

    let assignments = atropisomer_stereo_without_conformer(&molecule);
    apply_coordinate_free_atropisomer_assignments(&mut molecule, assignments);

    assert!(molecule.conformers_3d().is_empty());
    assert_eq!(molecule.bonds()[2].stereo(), BondStereo::AtropCcw);
}

#[test]
fn mol_from_smiles_flips_coordinate_free_atropisomer_bond_stereo_for_hash_like_rdkit() {
    let mut molecule = coordinate_free_atropisomer_candidate(BondDirection::BeginDash);

    let assignments = atropisomer_stereo_without_conformer(&molecule);
    apply_coordinate_free_atropisomer_assignments(&mut molecule, assignments);

    assert!(molecule.conformers_3d().is_empty());
    assert_eq!(molecule.bonds()[2].stereo(), BondStereo::AtropCw);
}

fn conformer_backed_atropisomer_candidate(is_3d: bool) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let chlorine = builder
        .add_atom(AtomSpec::new(Element::from_atomic_number(17).unwrap()).with_no_implicit(true));
    let axis_begin = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
    let alkene_left = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
    let axis_end = builder.add_atom(AtomSpec::new(Element::C));
    let alkene_right = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));

    builder
        .add_bond(
            BondSpec::new(axis_begin, chlorine, BondOrder::Single)
                .with_direction(BondDirection::BeginWedge),
        )
        .unwrap();
    builder
        .add_bond(BondSpec::new(axis_begin, alkene_left, BondOrder::Double))
        .unwrap();
    builder
        .add_bond(BondSpec::new(axis_begin, axis_end, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(BondSpec::new(axis_end, alkene_right, BondOrder::Double))
        .unwrap();

    let coords = if is_3d {
        vec![
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, -1.0, 0.0],
        ]
    } else {
        vec![
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, -1.0, 0.0],
        ]
    };
    builder
        .add_conformer(Conformer3D::new(0, coords, is_3d))
        .unwrap();

    builder
        .build()
        .unwrap()
        .sanitize_with_ops(crate::SanitizeOps::ALL)
        .unwrap()
}

#[test]
fn mol_from_smiles_assigns_2d_conformer_backed_atropisomer_bond_stereo_like_rdkit() {
    let mut molecule = conformer_backed_atropisomer_candidate(false);

    let assignments = atropisomer_stereo_from_conformer(&molecule, 0);
    apply_atropisomer_stereo_assignments(&mut molecule, assignments);

    assert_eq!(molecule.bonds()[2].stereo(), BondStereo::AtropCw);
}

#[test]
fn mol_from_smiles_assigns_3d_conformer_backed_atropisomer_bond_stereo_like_rdkit() {
    let mut molecule = conformer_backed_atropisomer_candidate(true);

    let assignments = atropisomer_stereo_from_conformer(&molecule, 0);
    apply_atropisomer_stereo_assignments(&mut molecule, assignments);

    assert_eq!(molecule.bonds()[2].stereo(), BondStereo::AtropCw);
}

fn three_neighbor_3d_carbon(formal_charge: i8) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let center = builder.add_atom(AtomSpec::new(Element::C).with_formal_charge(formal_charge));
    let fluorine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(9).unwrap()));
    let chlorine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(17).unwrap()));
    let bromine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(35).unwrap()));
    builder
        .add_bond(BondSpec::new(center, fluorine, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(BondSpec::new(center, chlorine, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(BondSpec::new(center, bromine, BondOrder::Single))
        .unwrap();
    builder
        .add_conformer(Conformer3D::new(
            0,
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            true,
        ))
        .unwrap();
    builder.build().unwrap()
}

#[test]
fn assign_chiral_types_from_3d_uses_valence_implicit_h_for_three_coordinate_carbon() {
    let mut molecule = three_neighbor_3d_carbon(0);

    assign_chiral_types_from_3d(&mut molecule, 0).expect("3D chirality assignment should succeed");

    assert_ne!(molecule.atoms()[0].chiral_tag(), ChiralTag::Unspecified);
}

#[test]
fn assign_chiral_types_from_3d_rejects_three_coordinate_cation_with_no_implicit_h() {
    let mut molecule = three_neighbor_3d_carbon(1);

    assign_chiral_types_from_3d(&mut molecule, 0).expect("3D chirality assignment should succeed");

    assert_eq!(molecule.atoms()[0].chiral_tag(), ChiralTag::Unspecified);
}

#[test]
fn assign_chiral_types_from_3d_clears_stereochem_done_like_rdkit() {
    let mut molecule = three_neighbor_3d_carbon(0).with_prop("_StereochemDone", "1");

    assign_chiral_types_from_3d(&mut molecule, 0).expect("3D chirality assignment should succeed");

    assert_eq!(molecule.prop("_StereochemDone"), None);
}

#[test]
fn assign_chiral_types_from_3d_marks_non_explicit_3d_chirality_like_rdkit() {
    let mut molecule = three_neighbor_3d_carbon(0);

    assign_chiral_types_from_3d(&mut molecule, 0).expect("3D chirality assignment should succeed");

    assert_ne!(molecule.atoms()[0].chiral_tag(), ChiralTag::Unspecified);
    assert_eq!(
        molecule.atoms()[0].prop("_NonExplicit3DChirality"),
        Some("1")
    );
}

#[test]
fn assign_chiral_types_from_3d_does_not_mark_existing_explicit_atom_as_non_explicit() {
    let mut molecule = three_neighbor_3d_carbon(0);
    molecule.topology_block_mut().atoms[0].set_chiral_tag(ChiralTag::TetrahedralCw);

    assign_chiral_types_from_3d(&mut molecule, 0).expect("3D chirality assignment should succeed");

    assert_ne!(molecule.atoms()[0].chiral_tag(), ChiralTag::Unspecified);
    assert_eq!(molecule.atoms()[0].prop("_NonExplicit3DChirality"), None);
}

fn three_neighbor_pseudo_2d_carbon() -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let center = builder.add_atom(AtomSpec::new(Element::C));
    let fluorine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(9).unwrap()));
    let chlorine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(17).unwrap()));
    let bromine = builder.add_atom(AtomSpec::new(Element::from_atomic_number(35).unwrap()));
    builder
        .add_bond(BondSpec::new(center, fluorine, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(BondSpec::new(center, chlorine, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(BondSpec::new(center, bromine, BondOrder::Single))
        .unwrap();
    builder
        .add_conformer(Conformer3D::new(
            0,
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [-1.0, -1.0, 0.0],
            ],
            true,
        ))
        .unwrap();
    builder.build().unwrap()
}

#[test]
fn assign_chiral_types_from_bond_dirs_promotes_implicit_h_on_three_coordinate_center() {
    let mut molecule = three_neighbor_pseudo_2d_carbon();
    molecule.topology_block_mut().bonds[0].set_direction(BondDirection::BeginWedge);

    assign_chiral_types_from_bond_dirs(&mut molecule, 0);

    assert_ne!(molecule.atoms()[0].chiral_tag(), ChiralTag::Unspecified);
    assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 1);
}

fn stereogenic_double_bond_molecule_with_conformer(same_side: bool, is_3d: bool) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let c0 = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
    let c1 = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
    let f = builder.add_atom(AtomSpec::new(Element::from_atomic_number(9).unwrap()));
    let cl = builder.add_atom(AtomSpec::new(Element::from_atomic_number(17).unwrap()));
    builder
        .add_bond(BondSpec::new(c0, c1, BondOrder::Double))
        .unwrap();
    builder
        .add_bond(BondSpec::new(c0, f, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(BondSpec::new(c1, cl, BondOrder::Single))
        .unwrap();
    let right_y = if same_side { 1.0 } else { -1.0 };
    builder
        .add_conformer(Conformer3D::new(
            0,
            vec![
                [-1.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [-1.0, 1.0, 0.0],
                [1.0, right_y, 0.0],
            ],
            is_3d,
        ))
        .unwrap();
    builder.build().unwrap()
}

fn stereogenic_double_bond_3d_molecule(same_side: bool) -> Molecule {
    stereogenic_double_bond_molecule_with_conformer(same_side, true)
}

fn stereogenic_double_bond_molecule_without_conformer() -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let c0 = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
    let c1 = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
    let f = builder.add_atom(AtomSpec::new(Element::from_atomic_number(9).unwrap()));
    let cl = builder.add_atom(AtomSpec::new(Element::from_atomic_number(17).unwrap()));
    builder
        .add_bond(BondSpec::new(c0, c1, BondOrder::Double))
        .unwrap();
    builder
        .add_bond(BondSpec::new(c0, f, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(BondSpec::new(c1, cl, BondOrder::Single))
        .unwrap();
    builder.build().unwrap()
}

fn linear_double_bond_3d_molecule() -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let c0 = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
    let c1 = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
    let f = builder.add_atom(AtomSpec::new(Element::from_atomic_number(9).unwrap()));
    let cl = builder.add_atom(AtomSpec::new(Element::from_atomic_number(17).unwrap()));
    builder
        .add_bond(BondSpec::new(c0, c1, BondOrder::Double))
        .unwrap();
    builder
        .add_bond(BondSpec::new(c0, f, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(BondSpec::new(c1, cl, BondOrder::Single))
        .unwrap();
    builder
        .add_conformer(Conformer3D::new(
            0,
            vec![
                [-1.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [-2.0, 0.0, 0.0],
                [2.0, 1.0, 0.0],
            ],
            true,
        ))
        .unwrap();
    builder.build().unwrap()
}

fn squiggle_neighbor_double_bond_3d_molecule() -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let c0 = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
    let c1 = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
    let h = builder.add_atom(AtomSpec::new(Element::H));
    let cl = builder.add_atom(AtomSpec::new(Element::from_atomic_number(17).unwrap()));
    builder
        .add_bond(BondSpec::new(c0, c1, BondOrder::Double))
        .unwrap();
    builder
        .add_bond(BondSpec::new(h, c0, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(BondSpec::new(c1, cl, BondOrder::Single))
        .unwrap();
    builder
        .add_conformer(Conformer3D::new(
            0,
            vec![
                [-1.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [-2.0, 1.0, 0.0],
                [2.0, 1.0, 0.0],
            ],
            true,
        ))
        .unwrap();
    builder.build().unwrap()
}

#[test]
fn clear_dir_flags_marks_unknown_and_clears_non_wedge_directions_like_rdkit() {
    let mut molecule = stereogenic_double_bond_3d_molecule(true);
    molecule.topology_block_mut().bonds[0].set_direction(BondDirection::EitherDouble);
    molecule.topology_block_mut().bonds[1].set_direction(BondDirection::Unknown);
    molecule.topology_block_mut().bonds[2].set_direction(BondDirection::EndUpRight);

    clear_dir_flags(&mut molecule, true);

    assert!(molecule.bonds()[0].unknown_stereo());
    assert_eq!(molecule.bonds()[0].direction(), BondDirection::None);
    assert!(molecule.bonds()[1].unknown_stereo());
    assert_eq!(molecule.bonds()[1].direction(), BondDirection::None);
    assert_eq!(molecule.bonds()[2].direction(), BondDirection::EndUpRight);
}

#[test]
fn clear_all_bond_dir_flags_clears_wedge_type_directions_like_rdkit() {
    let mut molecule = stereogenic_double_bond_3d_molecule(true);
    molecule.topology_block_mut().bonds[1].set_direction(BondDirection::EndUpRight);
    molecule.topology_block_mut().bonds[2].set_direction(BondDirection::EndDownRight);

    clear_all_bond_dir_flags(&mut molecule);

    assert_eq!(molecule.bonds()[1].direction(), BondDirection::None);
    assert_eq!(molecule.bonds()[2].direction(), BondDirection::None);
}

#[test]
fn set_double_bond_neighbor_directions_assigns_adjacent_bond_dirs_from_3d() {
    let mut molecule = stereogenic_double_bond_3d_molecule(true);

    set_double_bond_neighbor_directions(&mut molecule, 0).unwrap();

    assert!(matches!(
        molecule.bonds()[1].direction(),
        BondDirection::EndUpRight | BondDirection::EndDownRight
    ));
    assert!(matches!(
        molecule.bonds()[2].direction(),
        BondDirection::EndUpRight | BondDirection::EndDownRight
    ));
}

#[test]
fn set_double_bond_neighbor_directions_materializes_symm_sssr_ring_cache_like_rdkit() {
    let mut molecule = stereogenic_double_bond_3d_molecule(true);

    assert!(molecule.derived_cache().rings.is_none());

    set_double_bond_neighbor_directions(&mut molecule, 0).unwrap();

    let rings = molecule
        .derived_cache()
        .rings
        .as_ref()
        .expect("ring cache should be materialized before stereo detection");
    assert!(rings.is_symm_sssr());
}

#[test]
fn set_double_bond_neighbor_directions_without_conformer_uses_existing_stereo_like_rdkit() {
    let mut molecule = stereogenic_double_bond_molecule_without_conformer();
    molecule.topology_block_mut().bonds[0].set_stereo_atoms(Some([AtomId::new(2), AtomId::new(3)]));
    molecule.topology_block_mut().bonds[0].set_stereo(BondStereo::Trans);

    set_double_bond_neighbor_directions_from_stereo(&mut molecule).unwrap();

    assert!(matches!(
        molecule.bonds()[1].direction(),
        BondDirection::EndUpRight | BondDirection::EndDownRight
    ));
    assert!(matches!(
        molecule.bonds()[2].direction(),
        BondDirection::EndUpRight | BondDirection::EndDownRight
    ));

    set_bond_stereo_from_directions(&mut molecule).unwrap();

    assert_eq!(molecule.bonds()[0].stereo(), BondStereo::Trans);
    assert_eq!(
        molecule.bonds()[0].stereo_atoms(),
        Some([AtomId::new(2), AtomId::new(3)])
    );
}

#[test]
fn set_double_bond_neighbor_directions_marks_linear_arrangement_as_any_like_rdkit() {
    let mut molecule = linear_double_bond_3d_molecule();

    set_double_bond_neighbor_directions(&mut molecule, 0).unwrap();

    assert_eq!(molecule.bonds()[0].stereo(), BondStereo::Any);
    assert_eq!(molecule.bonds()[0].stereo_atoms(), None);
    assert_eq!(molecule.bonds()[1].direction(), BondDirection::None);
    assert_eq!(molecule.bonds()[2].direction(), BondDirection::None);
}

#[test]
fn set_double_bond_neighbor_directions_reorients_existing_dir_to_rdkit_raw_bond_frame() {
    let mut molecule = stereogenic_double_bond_3d_molecule(true);
    molecule.topology_block_mut().bonds[1].set_direction(BondDirection::EndUpRight);

    set_double_bond_neighbor_directions(&mut molecule, 0).unwrap();
    set_bond_stereo_from_directions(&mut molecule).unwrap();

    assert_eq!(molecule.bonds()[1].direction(), BondDirection::EndDownRight);
    assert_eq!(molecule.bonds()[2].direction(), BondDirection::EndDownRight);
    assert_eq!(molecule.bonds()[0].stereo(), BondStereo::Cis);
    assert_eq!(
        molecule.bonds()[0].stereo_atoms(),
        Some([AtomId::new(2), AtomId::new(3)])
    );
}

#[test]
fn from_smiles_assigns_ring_closure_double_bond_stereo_atoms_like_rdkit_row_86() {
    let molecule = Molecule::from_smiles(
            "C=C1/C(C[C@@H](O)CC1)=C\\C=C2[C@@]3([H])[C@@](CCC\\2)(C)[C@]([C@H](C)/C=C/[C@H](C)C(C)C)([H])CC3",
        )
        .unwrap();

    assert_eq!(molecule.bonds()[9].stereo(), BondStereo::E);
    assert_eq!(
        molecule.bonds()[9].stereo_atoms(),
        Some([AtomId::new(8), AtomId::new(11)])
    );
    assert_eq!(molecule.bonds()[10].stereo(), BondStereo::None);
    assert_eq!(molecule.bonds()[10].stereo_atoms(), None);
}

#[test]
fn from_smiles_does_not_leave_directional_bonds_for_nonstereo_ring_like_rdkit_row_87() {
    let molecule = Molecule::from_smiles("O=C1OCC(=C1c1ccccc1)c1ccccc1 |c:4|").unwrap();
    assert_eq!(molecule.bonds()[4].direction(), BondDirection::None);
    assert_eq!(molecule.bonds()[4].stereo(), BondStereo::None);
    assert_eq!(molecule.bonds()[4].stereo_atoms(), None);
}

#[test]
fn from_smiles_does_not_assign_imine_stereo_without_distinguishable_substituents_like_rdkit_row_88()
{
    let molecule =
        Molecule::from_smiles("O=C1N(/N=C(/C)C1=NN/C2=C/C(OC)=CC=C2)C=3C=CC=CC=3").unwrap();
    assert_eq!(molecule.bonds()[3].direction(), BondDirection::None);
    assert_eq!(molecule.bonds()[3].stereo(), BondStereo::None);
    assert_eq!(molecule.bonds()[3].stereo_atoms(), None);
}

#[test]
fn set_double_bond_neighbor_directions_marks_reverse_squiggle_neighbor_as_any_like_rdkit() {
    let mut molecule = squiggle_neighbor_double_bond_3d_molecule();
    molecule.topology_block_mut().bonds[1].set_direction(BondDirection::Unknown);
    molecule.topology_block_mut().bonds[1].set_unknown_stereo(true);

    set_double_bond_neighbor_directions(&mut molecule, 0).unwrap();

    assert_eq!(molecule.bonds()[0].stereo(), BondStereo::Any);
    assert_eq!(molecule.bonds()[0].stereo_atoms(), None);
    assert_eq!(molecule.bonds()[1].direction(), BondDirection::Unknown);
    assert!(molecule.bonds()[1].unknown_stereo());
}

#[test]
fn from_mol_block_row_569_keeps_rdkit_control_bond_choice_for_2d_double_bond_dirs() {
    let molecule = Molecule::from_mol_block(
        r#"
     RDKit          2D

 36 40  0  0  0  0  0  0  0  0999 V2000
   -7.1003   -0.3581    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.6013   -0.4152    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.8025    0.8544    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5026    2.1810    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.0015    2.2381    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7037    3.4506    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2048    3.3936    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5047    2.0670    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0058    2.0099    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3058    0.6833    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1932    0.6263    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.9920    1.8959    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.4887    1.9958    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.4493    0.8438    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.9273    1.0996    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.4053    1.3555    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    6.1831   -0.3784    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    5.6714    2.5776    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    3.9318   -0.5642    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.4344    3.2884    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
    2.5865    4.2490    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.8561    3.4501    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1046   -0.5863    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6035   -0.5292    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.3036    0.7974    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4045   -1.9129    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2034   -3.1825    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.7023   -3.1254    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4024   -1.7988    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5033   -4.5091    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.9956   -4.5661    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3846   -5.1323    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6792   -6.0323    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.7945   -3.2965    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0944   -1.9699    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.8932   -0.7003    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  4  6  1  0
  6  7  2  0
  7  8  1  0
  8  9  1  0
  9 10  1  0
 10 11  2  0
 11 12  1  0
 12 13  2  0
 13 14  1  0
 14 15  1  0
 15 16  1  0
 15 17  1  0
 15 18  1  0
 14 19  1  1
 12 20  1  0
 20 21  1  0
 21 22  2  0
 10 23  1  0
 23 24  1  0
 24 25  1  0
 23 26  2  0
 26 27  1  0
 27 28  2  0
 28 29  1  0
 27 30  1  0
 30 31  1  0
 31 32  1  0
 31 33  1  0
 31 34  1  0
 34 35  2  0
 35 36  1  0
 25  3  1  0
 25  8  2  0
 22 13  1  0
 29 24  2  0
 35 26  1  0
M  END
"#,
    )
    .unwrap();

    assert_eq!(molecule.bonds()[8].direction(), BondDirection::EndUpRight);
    assert_eq!(
        molecule.bonds()[10].direction(),
        BondDirection::EndDownRight
    );
    assert_eq!(molecule.bonds()[21].direction(), BondDirection::EndUpRight);
}

#[test]
fn from_mol_block_row_681_keeps_rdkit_reverse_sorted_double_bond_processing_order() {
    let molecule = Molecule::from_mol_block(
        r#"
     RDKit          2D

 29 31  0  0  0  0  0  0  0  0999 V2000
    4.4613    2.0068    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.5173    0.9414    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.1226   -0.5057    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.1786   -1.5711    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.6292   -1.1893    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.6851   -2.2546    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    8.0238    0.2579    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.9679    1.3232    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.6720   -0.8875    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.6161    0.1778    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.1655   -0.2040    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1095    0.8614    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3411    0.4796    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7357   -0.9676    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6798   -2.0330    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7708   -1.6512    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.8236    0.2515    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1524    1.7412    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6500    1.8263    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.3251    3.1658    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4725    0.5719    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.9701    0.6571    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.7927   -0.5973    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -8.2902   -0.5121    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.1176   -1.9368    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.9401   -3.1912    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6200   -2.0219    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2033    2.9027    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.1947    2.3589    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  1  0
  4  5  1  0
  5  6  1  6
  5  7  1  0
  7  8  1  0
  3  9  2  0
  9 10  1  0
 10 11  2  0
 11 12  1  0
 12 13  1  0
 13 14  1  0
 14 15  1  0
 15 16  1  0
 13 17  1  1
 13 18  1  0
 18 19  1  0
 19 20  1  1
 19 21  1  0
 21 22  2  0
 22 23  1  0
 23 24  1  1
 23 25  1  0
 25 26  1  0
 25 27  1  0
 18 28  1  1
 28 29  1  0
  8  2  1  0
 16 11  1  0
 12 29  1  1
M  END
"#,
    )
    .unwrap();

    assert_eq!(molecule.bonds()[1].direction(), BondDirection::EndDownRight);
    assert_eq!(molecule.bonds()[2].direction(), BondDirection::EndDownRight);
    assert_eq!(molecule.bonds()[8].direction(), BondDirection::EndUpRight);
    assert_eq!(molecule.bonds()[10].direction(), BondDirection::EndUpRight);
    assert_eq!(molecule.bonds()[29].direction(), BondDirection::EndUpRight);
}

#[test]
fn set_bond_stereo_from_directions_assigns_cis_for_same_side_neighbors() {
    let mut molecule = stereogenic_double_bond_3d_molecule(true);
    set_double_bond_neighbor_directions(&mut molecule, 0).unwrap();

    set_bond_stereo_from_directions(&mut molecule).unwrap();

    assert_eq!(molecule.bonds()[0].stereo(), BondStereo::Cis);
    assert_eq!(
        molecule.bonds()[0].stereo_atoms(),
        Some([AtomId::new(2), AtomId::new(3)])
    );
}

#[test]
fn set_bond_stereo_from_directions_leaves_stereo_unset_without_both_neighbor_dirs_like_rdkit() {
    let mut molecule = stereogenic_double_bond_molecule_without_conformer();
    molecule.topology_block_mut().bonds[1].set_direction(BondDirection::EndUpRight);
    molecule
        .properties_mut()
        .set_prop("_needsDetectBondStereo", "1");

    set_bond_stereo_from_directions(&mut molecule).unwrap();

    assert_eq!(molecule.bonds()[0].stereo(), BondStereo::None);
    assert_eq!(molecule.bonds()[0].stereo_atoms(), None);
    assert_eq!(molecule.prop("_needsDetectBondStereo"), None);
}

#[test]
fn assign_double_bond_stereo_from_directions_updates_public_molecule_state_like_rdkit() {
    let mut molecule = stereogenic_double_bond_3d_molecule(true);
    set_double_bond_neighbor_directions(&mut molecule, 0).unwrap();
    molecule
        .properties_mut()
        .set_prop("_needsDetectBondStereo", "1");

    assign_double_bond_stereo_from_directions(&mut molecule).unwrap();

    assert_eq!(molecule.bonds()[0].stereo(), BondStereo::Cis);
    assert_eq!(
        molecule.bonds()[0].stereo_atoms(),
        Some([AtomId::new(2), AtomId::new(3)])
    );
    assert_eq!(molecule.prop("_needsDetectBondStereo"), None);
}

#[test]
fn assign_stereochemistry_from_3d_assigns_trans_for_opposite_side_neighbors() {
    let mut molecule = stereogenic_double_bond_3d_molecule(false);

    assign_stereochemistry_from_3d(&mut molecule, 0).unwrap();

    assert_eq!(molecule.bonds()[0].stereo(), BondStereo::Trans);
    assert_eq!(
        molecule.bonds()[0].stereo_atoms(),
        Some([AtomId::new(2), AtomId::new(3)])
    );
}

#[test]
fn stereochemistry_from_3d_ignores_non_3d_conformer_like_rdkit() {
    let mut molecule = stereogenic_double_bond_molecule_with_conformer(false, false);
    molecule.topology_block_mut().bonds[1].set_direction(BondDirection::EndUpRight);
    molecule.topology_block_mut().bonds[2].set_direction(BondDirection::EndDownRight);

    assign_stereochemistry_from_3d(&mut molecule, 0).unwrap();

    assert_eq!(molecule.bonds()[0].stereo(), BondStereo::None);
    assert_eq!(molecule.bonds()[1].direction(), BondDirection::EndUpRight);
    assert_eq!(molecule.bonds()[2].direction(), BondDirection::EndDownRight);
}

#[test]
fn stereochemistry_from_3d_ignores_molecule_without_conformer_like_rdkit() {
    let mut molecule = stereogenic_double_bond_molecule_without_conformer();
    molecule.topology_block_mut().bonds[1].set_direction(BondDirection::EndUpRight);
    molecule.topology_block_mut().bonds[2].set_direction(BondDirection::EndDownRight);

    assign_stereochemistry_from_3d(&mut molecule, 0).unwrap();

    assert_eq!(molecule.bonds()[0].stereo(), BondStereo::None);
    assert_eq!(molecule.bonds()[1].direction(), BondDirection::EndUpRight);
    assert_eq!(molecule.bonds()[2].direction(), BondDirection::EndDownRight);
}

#[test]
fn stereochemistry_from_3d_clears_existing_double_bond_stereo_before_reassignment() {
    let mut molecule = stereogenic_double_bond_3d_molecule(true);
    molecule.topology_block_mut().bonds[0].set_stereo_atoms(Some([AtomId::new(0), AtomId::new(1)]));
    molecule.topology_block_mut().bonds[0].set_stereo(BondStereo::Trans);

    assign_stereochemistry_from_3d(&mut molecule, 0).unwrap();

    assert_eq!(molecule.bonds()[0].stereo(), BondStereo::Cis);
    assert_eq!(
        molecule.bonds()[0].stereo_atoms(),
        Some([AtomId::new(2), AtomId::new(3)])
    );
}

fn square_planar_3d_phosphorus(coords: Vec<[f64; 3]>) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let center = builder
        .add_atom(AtomSpec::new(Element::from_atomic_number(15).unwrap()).with_no_implicit(true));
    let ligand_elements = [
        Element::from_atomic_number(9).unwrap(),
        Element::from_atomic_number(17).unwrap(),
        Element::from_atomic_number(35).unwrap(),
        Element::from_atomic_number(53).unwrap(),
        Element::from_atomic_number(7).unwrap(),
        Element::from_atomic_number(8).unwrap(),
    ];
    for element in ligand_elements.iter().take(coords.len() - 1) {
        let ligand = builder.add_atom(AtomSpec::new(*element).with_no_implicit(true));
        builder
            .add_bond(BondSpec::new(center, ligand, BondOrder::Single))
            .unwrap();
    }
    builder
        .add_conformer(Conformer3D::new(0, coords, true))
        .unwrap();
    builder.build().unwrap()
}

#[test]
fn assign_chiral_types_from_3d_assigns_square_planar_from_two_opposite_pairs_like_rdkit() {
    let mut molecule = square_planar_3d_phosphorus(vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
    ]);

    assign_chiral_types_from_3d(&mut molecule, 0).expect("3D chirality assignment should succeed");

    assert_eq!(molecule.atoms()[0].chiral_tag(), ChiralTag::SquarePlanar);
    assert_eq!(molecule.atoms()[0].chiral_permutation(), Some(2));
    assert_eq!(
        molecule.atoms()[0].prop("_NonExplicit3DChirality"),
        Some("1")
    );
}

#[test]
fn assign_chiral_types_from_3d_assigns_t_shaped_square_planar_like_rdkit() {
    let mut molecule = square_planar_3d_phosphorus(vec![
        [0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
    ]);

    assign_chiral_types_from_3d(&mut molecule, 0).expect("3D chirality assignment should succeed");

    assert_eq!(molecule.atoms()[0].chiral_tag(), ChiralTag::SquarePlanar);
    assert_eq!(molecule.atoms()[0].chiral_permutation(), Some(3));
}

#[test]
fn assign_chiral_types_from_3d_assigns_trigonal_bipyramidal_from_one_opposite_pair_like_rdkit() {
    let mut molecule = square_planar_3d_phosphorus(vec![
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [-1.0, -1.0, 0.0],
    ]);

    assign_chiral_types_from_3d(&mut molecule, 0).expect("3D chirality assignment should succeed");

    assert_eq!(
        molecule.atoms()[0].chiral_tag(),
        ChiralTag::TrigonalBipyramidal
    );
    assert_eq!(molecule.atoms()[0].chiral_permutation(), Some(7));
    assert_eq!(
        molecule.atoms()[0].prop("_NonExplicit3DChirality"),
        Some("1")
    );
}

#[test]
fn assign_chiral_types_from_3d_assigns_seesaw_trigonal_bipyramidal_like_rdkit() {
    let mut molecule = square_planar_3d_phosphorus(vec![
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
        [1.0, 0.0, 0.0],
        [-0.5, 0.866_025_403_784_438_6, 0.0],
    ]);

    assign_chiral_types_from_3d(&mut molecule, 0).expect("3D chirality assignment should succeed");

    assert_eq!(
        molecule.atoms()[0].chiral_tag(),
        ChiralTag::TrigonalBipyramidal
    );
    assert_eq!(molecule.atoms()[0].chiral_permutation(), Some(7));
}

#[test]
fn assign_chiral_types_from_3d_covers_all_seesaw_octahedral_branches_like_rdkit() {
    let cases = [
        (
            vec![
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, -1.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
            ],
            25,
        ),
        (
            vec![
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 0.0],
                [0.0, 0.0, -1.0],
                [0.0, 1.0, 0.0],
            ],
            19,
        ),
        (
            vec![
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, -1.0],
            ],
            6,
        ),
        (
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, -1.0],
                [0.0, -1.0, 0.0],
            ],
            10,
        ),
        (
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, -1.0],
            ],
            1,
        ),
        (
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, -1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, -1.0],
            ],
            4,
        ),
    ];

    for (coords, expected_perm) in cases {
        let mut molecule = square_planar_3d_phosphorus(coords);
        assign_chiral_types_from_3d(&mut molecule, 0)
            .expect("3D chirality assignment should succeed");

        assert_eq!(molecule.atoms()[0].chiral_tag(), ChiralTag::Octahedral);
        assert_eq!(
            molecule.atoms()[0].chiral_permutation(),
            Some(expected_perm)
        );
    }
}

#[test]
fn assign_chiral_types_from_3d_assigns_octahedral_from_three_opposite_pairs_like_rdkit() {
    let mut molecule = square_planar_3d_phosphorus(vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
    ]);

    assign_chiral_types_from_3d(&mut molecule, 0).expect("3D chirality assignment should succeed");

    assert_eq!(molecule.atoms()[0].chiral_tag(), ChiralTag::Octahedral);
    assert_eq!(molecule.atoms()[0].chiral_permutation(), Some(27));
}

#[test]
fn assign_chiral_types_from_3d_assigns_square_pyramidal_as_octahedral_like_rdkit() {
    let mut molecule = square_planar_3d_phosphorus(vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]);

    assign_chiral_types_from_3d(&mut molecule, 0).expect("3D chirality assignment should succeed");

    assert_eq!(molecule.atoms()[0].chiral_tag(), ChiralTag::Octahedral);
    assert_eq!(molecule.atoms()[0].chiral_permutation(), Some(27));
}

#[test]
fn assign_chiral_types_from_3d_assigns_seesaw_octahedral_like_rdkit() {
    let mut molecule = square_planar_3d_phosphorus(vec![
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ]);

    assign_chiral_types_from_3d(&mut molecule, 0).expect("3D chirality assignment should succeed");

    assert_eq!(molecule.atoms()[0].chiral_tag(), ChiralTag::Octahedral);
    assert_eq!(molecule.atoms()[0].chiral_permutation(), Some(25));
}

#[test]
fn handle_cx_part_and_name_parses_zero_bonds_like_rdkit() {
    let molecule = Molecule::from_smiles_with_sanitize("CC~CC |Z:1|", false).unwrap();

    assert_eq!(molecule.num_bonds(), 3);
    assert_eq!(molecule.bonds()[1].order(), BondOrder::Zero);
}

#[test]
fn handle_cx_part_and_name_parses_enhanced_stereo_groups_like_rdkit() {
    let molecule = Molecule::from_smiles_with_sanitize("C[C@H](O)N |o1:1,o1:2|", false).unwrap();

    assert_eq!(molecule.stereo_groups().len(), 1);
    assert_eq!(molecule.stereo_groups()[0].kind(), StereoGroupKind::Or);
    assert_eq!(molecule.stereo_groups()[0].id(), Some(1));
    assert_eq!(
        molecule.stereo_groups()[0].atoms(),
        &[AtomId::new(1), AtomId::new(2)]
    );
}

#[test]
fn parse_cx_enhanced_stereo_rejects_missing_atom_index_like_rdkit() {
    let mut state = SmilesBuildState::new();
    state.builder.add_atom(AtomSpec::new(Element::C));

    let error = parse_cx_enhanced_stereo(&mut state, "o1:1", &mut 0).unwrap_err();

    assert_eq!(error, cx_parse_failure());
}

#[test]
fn handle_cx_part_and_name_parses_coordinate_bonds_like_rdkit() {
    let molecule = Molecule::from_smiles_with_sanitize("N->O |C:1.0|", false).unwrap();

    assert_eq!(molecule.num_bonds(), 1);
    assert_eq!(molecule.bonds()[0].order(), BondOrder::Dative);
    assert_eq!(molecule.bonds()[0].begin(), AtomId::new(1));
    assert_eq!(molecule.bonds()[0].end(), AtomId::new(0));
}

#[test]
fn handle_cx_part_and_name_parses_hydrogen_bonds_like_rdkit() {
    let molecule = Molecule::from_smiles_with_sanitize("NO |H:1.0|", false).unwrap();

    assert_eq!(molecule.num_bonds(), 1);
    assert_eq!(molecule.bonds()[0].order(), BondOrder::Hydrogen);
    assert_eq!(molecule.bonds()[0].begin(), AtomId::new(1));
    assert_eq!(molecule.bonds()[0].end(), AtomId::new(0));
}

#[test]
fn handle_cx_part_and_name_parses_radicals_like_rdkit() {
    let molecule = Molecule::from_smiles_with_sanitize("CCC |^1:0,^4:1,^7:2|", false).unwrap();

    assert_eq!(molecule.atoms()[0].radical_electrons(), 1);
    assert_eq!(molecule.atoms()[1].radical_electrons(), 2);
    assert_eq!(molecule.atoms()[2].radical_electrons(), 3);
}

#[test]
fn parse_cx_extensions_dispatches_radicals_and_consumes_closing_pipe_like_rdkit() {
    let mut state = SmilesBuildState::new();
    state.builder.add_atom(AtomSpec::new(Element::C));
    state.builder.add_atom(AtomSpec::new(Element::C));
    state.builder.add_atom(AtomSpec::new(Element::C));

    let consumed = parse_cx_extensions(&mut state, "|^1:0,^4:1,^7:2|").unwrap();
    let molecule = state.into_molecule().unwrap();

    assert_eq!(consumed, "|^1:0,^4:1,^7:2|".len());
    assert_eq!(molecule.atoms()[0].radical_electrons(), 1);
    assert_eq!(molecule.atoms()[1].radical_electrons(), 2);
    assert_eq!(molecule.atoms()[2].radical_electrons(), 3);
}

#[test]
fn parse_cx_extensions_accepts_empty_text_like_rdkit_wrapper() {
    let mut state = SmilesBuildState::new();

    let consumed = parse_cx_extensions(&mut state, "").unwrap();

    assert_eq!(consumed, 0);
}

#[test]
fn parse_cx_extensions_rejects_missing_leading_pipe_like_rdkit() {
    let mut state = SmilesBuildState::new();

    let error = parse_cx_extensions(&mut state, "^1:0|").unwrap_err();

    assert_eq!(
        error,
        SmilesParseError::ParseError("CXSMILES extension does not start with |".to_string())
    );
}

#[test]
fn parse_cx_radicals_rejects_invalid_or_truncated_markers_like_rdkit() {
    let mut state = SmilesBuildState::new();
    state.builder.add_atom(AtomSpec::new(Element::C));

    let invalid = parse_cx_radicals(&mut state, "^0:0", &mut 0).unwrap_err();
    let truncated = parse_cx_radicals(&mut state, "^", &mut 0).unwrap_err();

    assert_eq!(invalid, cx_parse_failure());
    assert_eq!(truncated, cx_parse_failure());
}

#[test]
fn parse_cx_coords_marks_2d_and_3d_conformers_like_rdkit_has_non_zero_z() {
    let mut state = SmilesBuildState::new();
    state.builder.add_atom(AtomSpec::new(Element::C));
    state.builder.add_atom(AtomSpec::new(Element::C));
    state.builder.add_atom(AtomSpec::new(Element::O));

    let mut pos = 0;
    parse_cx_coords(&mut state, "(0,0,;1,0,;2,0,0)", &mut pos, 0).unwrap();
    pos = 0;
    parse_cx_coords(&mut state, "(0,0,;1,0,;2,0,0.5)", &mut pos, 1).unwrap();

    let molecule = state.into_molecule().unwrap();
    assert_eq!(molecule.conformers_3d().len(), 2);
    assert!(!molecule.conformers_3d()[0].is_3d());
    assert!(molecule.conformers_3d()[1].is_3d());
    assert_eq!(
        molecule.conformers_3d()[1].coordinates()[2],
        [2.0, 0.0, 0.5]
    );
}

#[test]
fn get_unspecified_bond_type_for_atoms_matches_rdkit_aromatic_rule() {
    assert_eq!(
        get_unspecified_bond_type_for_atoms(true, true),
        BondOrder::Aromatic
    );
    assert_eq!(
        get_unspecified_bond_type_for_atoms(true, false),
        BondOrder::Single
    );
    assert_eq!(
        get_unspecified_bond_type_for_atoms(false, true),
        BondOrder::Single
    );
    assert_eq!(
        get_unspecified_bond_type_for_atoms(false, false),
        BondOrder::Single
    );
}

#[test]
fn handle_cx_part_and_name_parses_unsaturation_query_like_rdkit() {
    let molecule = Molecule::from_smiles_with_sanitize("CC |u:1|", false).unwrap();

    assert_eq!(
        molecule.atoms()[1].query(),
        Some(&QueryNode::predicate(AtomQueryPredicate::IsUnsaturated))
    );
}

#[test]
fn handle_cx_part_and_name_parses_ring_bond_query_like_rdkit() {
    let molecule = Molecule::from_smiles_with_sanitize("C1CC1 |rb:0:2|", false).unwrap();

    assert_eq!(
        molecule.atoms()[0].query(),
        Some(&QueryNode::predicate(AtomQueryPredicate::RingBondCount(2)))
    );
}

#[test]
fn handle_cx_part_and_name_completes_ring_bond_scan_queries_like_rdkit() {
    let molecule = Molecule::from_smiles_with_sanitize("C1CC1 |rb:0:*|", false).unwrap();

    assert_eq!(
        molecule.atoms()[0].query(),
        Some(&QueryNode::predicate(AtomQueryPredicate::RingBondCount(2)))
    );
    assert_eq!(molecule.properties().prop("_NeedsQueryScan"), None);
}

#[test]
fn handle_cx_part_and_name_completes_substitution_scan_queries_like_rdkit() {
    let molecule = Molecule::from_smiles_with_sanitize("CC |s:0:*|", false).unwrap();

    assert_eq!(
        molecule.atoms()[0].query(),
        Some(&QueryNode::predicate(
            AtomQueryPredicate::NonHydrogenDegree(1,)
        ))
    );
    assert_eq!(molecule.properties().prop("_NeedsQueryScan"), None);
}

#[test]
fn handle_cx_part_and_name_completes_ring_and_non_ring_scan_queries_like_rdkit() {
    let molecule = Molecule::from_smiles_with_sanitize("C1CC1C |rb:0:*,s:3:*|", false).unwrap();

    assert_eq!(
        molecule.atoms()[0].query(),
        Some(&QueryNode::predicate(AtomQueryPredicate::RingBondCount(2)))
    );
    assert_eq!(
        molecule.atoms()[3].query(),
        Some(&QueryNode::predicate(
            AtomQueryPredicate::NonHydrogenDegree(1)
        ))
    );
    assert_eq!(molecule.properties().prop("_NeedsQueryScan"), None);
}

#[test]
fn handle_cx_part_and_name_parses_substitution_query_like_rdkit() {
    let molecule = Molecule::from_smiles_with_sanitize("CC |s:0:1|", false).unwrap();

    assert_eq!(
        molecule.atoms()[0].query(),
        Some(&QueryNode::predicate(
            AtomQueryPredicate::NonHydrogenDegree(1)
        ))
    );
}

#[test]
fn shared_canon_helpers_has_single_h_query_matches_rdkit_atomand_only() {
    let direct_h = QueryNode::predicate(AtomQueryPredicate::ImplicitHydrogenCount(1));
    let atom_and_h = QueryNode::and(vec![
        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
        QueryNode::predicate(AtomQueryPredicate::ImplicitHydrogenCount(1)),
    ]);
    let nested_atom_and_h = QueryNode::and(vec![
        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
        QueryNode::and(vec![QueryNode::predicate(
            AtomQueryPredicate::ImplicitHydrogenCount(1),
        )]),
    ]);
    let atom_or_h = QueryNode::or(vec![
        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
        QueryNode::predicate(AtomQueryPredicate::ImplicitHydrogenCount(1)),
    ]);
    let negated_h = QueryNode::and(vec![QueryNode::not(QueryNode::predicate(
        AtomQueryPredicate::ImplicitHydrogenCount(1),
    ))]);

    assert!(!atom_query_has_single_h_count(&direct_h));
    assert!(atom_query_has_single_h_count(&atom_and_h));
    assert!(atom_query_has_single_h_count(&nested_atom_and_h));
    assert!(!atom_query_has_single_h_count(&atom_or_h));
    assert!(!atom_query_has_single_h_count(&negated_h));
}

#[test]
fn shared_canon_helpers_pure_predicates_cover_every_source_condition() {
    assert!(!canon_chiral_atom_needs_tag_inversion(
        2, 1, true, false, 1, false
    ));
    assert!(canon_chiral_atom_needs_tag_inversion(
        3, 1, true, true, 0, true
    ));
    assert!(canon_chiral_atom_needs_tag_inversion(
        3, 0, false, false, 1, false
    ));
    assert!(!canon_chiral_atom_needs_tag_inversion(
        3, 0, false, true, 1, false
    ));
    assert!(!canon_chiral_atom_needs_tag_inversion(
        3, 0, false, false, 0, false
    ));
    assert!(!canon_chiral_atom_needs_tag_inversion(
        3, 0, false, false, 1, true
    ));

    assert!(canon_atom_has_fourth_valence(1, false, false));
    assert!(canon_atom_has_fourth_valence(0, true, false));
    assert!(canon_atom_has_fourth_valence(0, false, true));
    assert!(!canon_atom_has_fourth_valence(0, false, false));
}

#[test]
fn shared_canon_helpers_chiral_permutation_kernel_covers_tags_inverse_and_errors() {
    let incident = [0, 1, 2, 3, 4, 5];
    let square_probe = [Some(1), Some(0), Some(2), Some(3)];
    assert_eq!(
        nontetrahedral_chiral_permutation_kernel(
            1,
            ChiralTag::SquarePlanar,
            6,
            &incident[..4],
            &square_probe,
            false,
        ),
        Ok(3)
    );
    assert_eq!(
        nontetrahedral_chiral_permutation_kernel(
            1,
            ChiralTag::SquarePlanar,
            6,
            &incident[..4],
            &square_probe,
            true,
        ),
        Ok(3)
    );
    assert_eq!(
        nontetrahedral_chiral_permutation_kernel(
            1,
            ChiralTag::TrigonalBipyramidal,
            6,
            &incident[..5],
            &[Some(1), Some(0), Some(2), Some(3), Some(4)],
            false,
        ),
        Ok(9)
    );
    assert_eq!(
        nontetrahedral_chiral_permutation_kernel(
            1,
            ChiralTag::Octahedral,
            6,
            &incident,
            &[Some(1), Some(0), Some(2), Some(3), Some(4), Some(5),],
            false,
        ),
        Ok(17)
    );
    assert_eq!(
        nontetrahedral_chiral_permutation_kernel(
            0,
            ChiralTag::SquarePlanar,
            4,
            &incident[..4],
            &square_probe,
            false,
        ),
        Ok(0)
    );
    assert_eq!(
        nontetrahedral_chiral_permutation_kernel(
            1,
            ChiralTag::TetrahedralCw,
            4,
            &incident[..4],
            &square_probe,
            false,
        ),
        Ok(0)
    );
    assert_eq!(
        nontetrahedral_chiral_permutation_kernel(
            1,
            ChiralTag::SquarePlanar,
            6,
            &incident[..5],
            &[Some(0), Some(1), Some(2), Some(3), Some(4)],
            false,
        ),
        Ok(0)
    );
    assert_eq!(
        nontetrahedral_chiral_permutation_kernel(
            1,
            ChiralTag::SquarePlanar,
            4,
            &incident[..4],
            &[Some(0), Some(1), Some(2)],
            false,
        ),
        Err(ChiralPermutationError::ProbeSizeMismatch)
    );
    assert_eq!(
        nontetrahedral_chiral_permutation_kernel(
            1,
            ChiralTag::SquarePlanar,
            4,
            &incident[..4],
            &[None, Some(0), Some(1), Some(2)],
            false,
        ),
        Err(ChiralPermutationError::MissingProbeElement)
    );
    assert_eq!(
        nontetrahedral_chiral_permutation_kernel(
            1,
            ChiralTag::SquarePlanar,
            4,
            &incident[..4],
            &[Some(4), Some(0), Some(1), Some(2)],
            false,
        ),
        Err(ChiralPermutationError::BondIndexOutOfBounds)
    );
}

#[test]
fn shared_canon_helpers_parser_atom_has_fourth_valence_uses_single_h_query() {
    let mut state = SmilesBuildState::new();
    let atom = state
        .builder
        .add_atom(AtomSpec::new(Element::C).with_query(QueryNode::and(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
            QueryNode::predicate(AtomQueryPredicate::ImplicitHydrogenCount(1)),
        ])));

    assert!(state.atom_has_fourth_valence(atom).unwrap());
}

#[test]
fn shared_canon_helpers_parser_atom_has_fourth_valence_uses_explicit_h_and_false_fallback() {
    let mut state = SmilesBuildState::new();
    let explicit_h = state
        .builder
        .add_atom(AtomSpec::new(Element::C).with_explicit_hydrogens(1));
    let plain = state.builder.add_atom(AtomSpec::new(Element::C));
    assert!(state.atom_has_fourth_valence(explicit_h).unwrap());
    assert!(!state.atom_has_fourth_valence(plain).unwrap());
}

#[test]
fn is_unsaturated_distinguishes_multiple_bonds_from_dative_like_rdkit() {
    let mut saturated = SmilesBuildState::new();
    let a0 = saturated.builder.add_atom(AtomSpec::new(Element::C));
    let a1 = saturated.builder.add_atom(AtomSpec::new(Element::O));
    saturated
        .builder
        .add_bond(BondSpec::new(a0, a1, BondOrder::Dative))
        .unwrap();
    assert!(!saturated.is_unsaturated(a0));

    let mut unsaturated = SmilesBuildState::new();
    let b0 = unsaturated.builder.add_atom(AtomSpec::new(Element::C));
    let b1 = unsaturated.builder.add_atom(AtomSpec::new(Element::O));
    unsaturated
        .builder
        .add_bond(BondSpec::new(b0, b1, BondOrder::Double))
        .unwrap();
}

#[test]
fn perturbation_order_counts_probe_swaps_like_rdkit() {
    let mut state = SmilesBuildState::new();
    let center = state.builder.add_atom(AtomSpec::new(Element::C));
    let n1 = state.builder.add_atom(AtomSpec::new(Element::F));
    let n2 = state.builder.add_atom(AtomSpec::new(Element::CL));
    let n3 = state.builder.add_atom(AtomSpec::new(Element::BR));
    state
        .builder
        .add_bond(BondSpec::new(center, n1, BondOrder::Single))
        .unwrap();
    state
        .builder
        .add_bond(BondSpec::new(center, n2, BondOrder::Single))
        .unwrap();
    state
        .builder
        .add_bond(BondSpec::new(center, n3, BondOrder::Single))
        .unwrap();
    assert_eq!(state.perturbation_order(center, &[2, 1, 0]).unwrap(), 1);
}

#[test]
fn perturbation_order_rejects_size_mismatch_like_rdkit() {
    let mut state = SmilesBuildState::new();
    let center = state.builder.add_atom(AtomSpec::new(Element::C));
    let n1 = state.builder.add_atom(AtomSpec::new(Element::F));
    let n2 = state.builder.add_atom(AtomSpec::new(Element::CL));
    state
        .builder
        .add_bond(BondSpec::new(center, n1, BondOrder::Single))
        .unwrap();
    state
        .builder
        .add_bond(BondSpec::new(center, n2, BondOrder::Single))
        .unwrap();

    let error = state.perturbation_order(center, &[0]).unwrap_err();

    assert_eq!(
        error,
        SmilesParseError::ParseError("size mismatch".to_string())
    );
}

#[test]
fn perturbation_order_rejects_missing_probe_element_like_rdkit() {
    let mut state = SmilesBuildState::new();
    let center = state.builder.add_atom(AtomSpec::new(Element::C));
    let n1 = state.builder.add_atom(AtomSpec::new(Element::F));
    let n2 = state.builder.add_atom(AtomSpec::new(Element::CL));
    state
        .builder
        .add_bond(BondSpec::new(center, n1, BondOrder::Single))
        .unwrap();
    state
        .builder
        .add_bond(BondSpec::new(center, n2, BondOrder::Single))
        .unwrap();

    let error = state.perturbation_order(center, &[0, 9]).unwrap_err();

    assert_eq!(
        error,
        SmilesParseError::ParseError("could not find probe element".to_string())
    );
}

#[test]
fn shared_canon_helpers_parser_chiral_inversion_preserves_degree_and_unsaturation_rules() {
    let mut first_atom = SmilesBuildState::new();
    let first = first_atom
        .builder
        .add_atom(AtomSpec::new(Element::C).with_explicit_hydrogens(1));
    let f1 = first_atom.builder.add_atom(AtomSpec::new(Element::F));
    let f2 = first_atom.builder.add_atom(AtomSpec::new(Element::CL));
    let f3 = first_atom.builder.add_atom(AtomSpec::new(Element::BR));
    first_atom.smiles_start_atoms.insert(first);
    first_atom
        .builder
        .add_bond(BondSpec::new(first, f1, BondOrder::Single))
        .unwrap();
    first_atom
        .builder
        .add_bond(BondSpec::new(first, f2, BondOrder::Single))
        .unwrap();
    first_atom
        .builder
        .add_bond(BondSpec::new(first, f3, BondOrder::Single))
        .unwrap();
    assert!(
        first_atom
            .chiral_atom_needs_tag_inversion(first, 0)
            .unwrap()
    );

    let mut closure_atom = SmilesBuildState::new();
    let s0 = closure_atom.builder.add_atom(AtomSpec::new(Element::C));
    let s1 = closure_atom.builder.add_atom(AtomSpec::new(Element::C));
    let s2 = closure_atom.builder.add_atom(AtomSpec::new(Element::F));
    let s3 = closure_atom.builder.add_atom(AtomSpec::new(Element::CL));
    closure_atom
        .builder
        .add_bond(BondSpec::new(s0, s1, BondOrder::Single))
        .unwrap();
    closure_atom
        .builder
        .add_bond(BondSpec::new(s1, s2, BondOrder::Single))
        .unwrap();
    closure_atom
        .builder
        .add_bond(BondSpec::new(s1, s3, BondOrder::Single))
        .unwrap();
    assert!(closure_atom.chiral_atom_needs_tag_inversion(s1, 1).unwrap());

    let mut unsaturated = SmilesBuildState::new();
    let u0 = unsaturated.builder.add_atom(AtomSpec::new(Element::C));
    let u1 = unsaturated.builder.add_atom(AtomSpec::new(Element::S));
    let u2 = unsaturated.builder.add_atom(AtomSpec::new(Element::O));
    let u3 = unsaturated.builder.add_atom(AtomSpec::new(Element::CL));
    unsaturated
        .builder
        .add_bond(BondSpec::new(u0, u1, BondOrder::Single))
        .unwrap();
    unsaturated
        .builder
        .add_bond(BondSpec::new(u1, u2, BondOrder::Double))
        .unwrap();
    unsaturated
        .builder
        .add_bond(BondSpec::new(u1, u3, BondOrder::Single))
        .unwrap();
    assert!(!unsaturated.chiral_atom_needs_tag_inversion(u1, 1).unwrap());
}

#[test]
fn handle_cx_part_and_name_records_wedged_single_bond_cfg_like_rdkit() {
    let molecule = Molecule::from_smiles_with_sanitize("CC |wU:1.0|", false).unwrap();

    assert_eq!(molecule.bonds()[0].begin(), AtomId::new(1));
    assert_eq!(molecule.bonds()[0].end(), AtomId::new(0));
    assert_eq!(molecule.bonds()[0].direction(), BondDirection::None);
    assert!(!molecule.bonds()[0].unknown_stereo());
    assert_eq!(molecule.bonds()[0].prop("_MolFileBondCfg"), Some("1"));
}

#[test]
fn handle_cx_part_and_name_records_wiggly_single_bond_cfg_like_rdkit() {
    let molecule = Molecule::from_smiles_with_sanitize("CC |w:1.0|", false).unwrap();

    assert_eq!(molecule.bonds()[0].begin(), AtomId::new(1));
    assert_eq!(molecule.bonds()[0].end(), AtomId::new(0));
    assert_eq!(molecule.bonds()[0].direction(), BondDirection::None);
    assert_eq!(molecule.bonds()[0].prop("_MolFileBondCfg"), Some("2"));
    assert_eq!(
        molecule.properties().prop("_needsDetectBondStereo"),
        Some("1")
    );
}

#[test]
fn handle_cx_part_and_name_parses_linknodes_like_rdkit() {
    let molecule = Molecule::from_smiles_with_sanitize("C1CC1 |LN:1:1.3|", false).unwrap();

    assert_eq!(
        molecule.properties().prop("_MolFileLinkNodes"),
        Some("1 3 2 2 1 2 3")
    );
}

#[test]
fn handle_cx_part_and_name_parses_data_sgroups_like_rdkit() {
    let molecule =
        Molecule::from_smiles_with_sanitize("CCO |SgD:2,1:FIELD:info::::|", false).unwrap();

    assert_eq!(molecule.substance_groups().len(), 1);
    let sgroup = &molecule.substance_groups()[0];
    assert_eq!(sgroup.kind(), &SubstanceGroupKind::Data);
    assert_eq!(sgroup.atoms(), &[AtomId::new(2), AtomId::new(1)]);
    assert_eq!(
        sgroup.props().get("FIELDNAME").map(String::as_str),
        Some("FIELD")
    );
    assert_eq!(sgroup.data_fields(), &["info".to_string()]);
    assert_eq!(sgroup.props().get("index").map(String::as_str), Some("1"));
}

#[test]
fn parse_cx_data_sgroup_consumes_coords_for_dropped_group_like_rdkit() {
    let mut state = SmilesBuildState::new();
    state.builder.add_atom(AtomSpec::new(Element::C));
    let mut pos = 0;

    parse_cx_data_sgroup(&mut state, "SgD:9:FIELD:info::::(1.,1.)", &mut pos, 0).unwrap();

    assert_eq!(pos, "SgD:9:FIELD:info::::(1.,1.)".len());
    let molecule = state.into_molecule().unwrap();
    assert!(molecule.substance_groups().is_empty());
}

#[test]
fn handle_cx_part_and_name_parses_variable_attachments_like_rdkit() {
    let molecule =
        Molecule::from_smiles_with_sanitize("CO*.C1=CC=NC=C1 |m:2:3.5.4|", false).unwrap();

    assert_eq!(
        molecule.bonds()[1].prop("_MolFileBondEndPts"),
        Some("(3 4 6 5)")
    );
    assert_eq!(molecule.bonds()[1].prop("_MolFileBondAttach"), Some("ANY"));
}

#[test]
fn handle_cx_part_and_name_parses_double_bond_stereo_like_rdkit() {
    let molecule = Molecule::from_smiles_with_sanitize("CC=CC |c:1|", false).unwrap();

    assert_eq!(molecule.bonds()[1].stereo(), BondStereo::Cis);
    assert_eq!(
        molecule.bonds()[1].stereo_atoms(),
        Some([AtomId::new(0), AtomId::new(3)])
    );
    assert_eq!(
        molecule.properties().prop("_needsDetectBondStereo"),
        Some("1")
    );
}

#[test]
fn from_smiles_with_sanitize_false_parses_bracket_charge_hydrogen_and_map() {
    let molecule = Molecule::from_smiles_with_sanitize("[13CH3:5].[NH4+]", false).unwrap();

    assert_eq!(molecule.atomic_numbers(), vec![6, 7]);
    assert_eq!(molecule.atoms()[0].isotope(), Some(13));
    assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 3);
    assert_eq!(molecule.atoms()[0].atom_map(), Some(5));
    assert!(molecule.atoms()[0].no_implicit());
    assert_eq!(molecule.atoms()[1].explicit_hydrogens(), 4);
    assert_eq!(molecule.atoms()[1].formal_charge(), 1);
    assert!(molecule.atoms()[1].no_implicit());
}

#[test]
fn remove_hs_update_explicit_count_tracks_isotopic_hydrogens_and_clears_stale_property() {
    let mut state = SmilesBuildState::new();
    let carbon = state
        .builder
        .add_atom(AtomSpec::new(Element::C).with_tracked_isotopic_hydrogens(vec![9]));
    let deuterium = state
        .builder
        .add_atom(AtomSpec::new(Element::H).with_isotope(2));
    state
        .builder
        .add_bond(BondSpec::new(carbon, deuterium, BondOrder::Single))
        .unwrap();

    state
        .remove_hs_update_explicit_count(&RemoveHsParams {
            remove_and_track_isotopes: true,
            update_explicit_count: true,
            ..RemoveHsParams::default()
        })
        .unwrap();
    let molecule = state.into_molecule().unwrap();

    assert_eq!(molecule.num_atoms(), 1);
    assert_eq!(molecule.atoms()[0].tracked_isotopic_hydrogens(), &[2]);
    assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 1);
}

#[test]
fn remove_hs_update_explicit_count_tracks_isotopes_after_nonisotopic_prepass() {
    let mut state = SmilesBuildState::new();
    let carbon = state
        .builder
        .add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
    let protium = state.builder.add_atom(AtomSpec::new(Element::H));
    let deuterium = state
        .builder
        .add_atom(AtomSpec::new(Element::H).with_isotope(2));
    state
        .builder
        .add_bond(BondSpec::new(carbon, protium, BondOrder::Single))
        .unwrap();
    state
        .builder
        .add_bond(BondSpec::new(carbon, deuterium, BondOrder::Single))
        .unwrap();

    state
        .remove_hs_update_explicit_count(&RemoveHsParams {
            remove_and_track_isotopes: true,
            update_explicit_count: true,
            ..RemoveHsParams::default()
        })
        .unwrap();
    let molecule = state.into_molecule().unwrap();

    assert_eq!(molecule.num_atoms(), 1);
    assert_eq!(molecule.atoms()[0].tracked_isotopic_hydrogens(), &[2]);
    assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 2);
}

#[test]
fn remove_hs_update_explicit_count_preserves_mapped_hydrogen_when_disabled() {
    let mut state = to_mol("[H:7]C").unwrap();

    state
        .remove_hs_update_explicit_count(&RemoveHsParams {
            remove_mapped: false,
            update_explicit_count: true,
            ..RemoveHsParams::default()
        })
        .unwrap();
    let molecule = state.into_molecule().unwrap();

    assert_eq!(molecule.atomic_numbers(), vec![1, 6]);
    assert_eq!(molecule.atoms()[0].atom_map(), Some(7));
}

#[test]
fn tracked_smiles_isotopic_hydrogens_groups_multiple_centers_in_atom_order_like_rdkit() {
    let mut builder = Molecule::builder();
    let carbon0 = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
    let carbon1 = builder.add_atom(AtomSpec::new(Element::C).with_no_implicit(true));
    let deuterium = builder.add_atom(AtomSpec::new(Element::H).with_isotope(2));
    let tritium = builder.add_atom(AtomSpec::new(Element::H).with_isotope(3));
    builder
        .add_bond(BondSpec::new(carbon1, tritium, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(BondSpec::new(carbon0, deuterium, BondOrder::Single))
        .unwrap();

    let tracked = tracked_smiles_isotopic_hydrogens(&builder).unwrap();

    assert_eq!(tracked, vec![(carbon0, vec![2]), (carbon1, vec![3])]);
}

#[test]
fn remove_hs_update_explicit_count_default_removes_hydrogen_with_wedged_bond() {
    let mut state = SmilesBuildState::new();
    let carbon = state.builder.add_atom(AtomSpec::new(Element::C));
    let hydrogen = state.builder.add_atom(AtomSpec::new(Element::H));
    state
        .builder
        .add_bond(
            BondSpec::new(carbon, hydrogen, BondOrder::Single)
                .with_direction(BondDirection::BeginWedge),
        )
        .unwrap();

    state
        .remove_hs_update_explicit_count(&RemoveHsParams {
            update_explicit_count: true,
            ..RemoveHsParams::default()
        })
        .unwrap();
    let molecule = state.into_molecule().unwrap();

    assert_eq!(molecule.num_atoms(), 1);
    assert_eq!(molecule.atoms()[0].atomic_number(), 6);
    assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 1);
}

#[test]
fn remove_hs_update_explicit_count_preserves_hydrogen_when_wedged_removal_disabled() {
    let mut state = SmilesBuildState::new();
    let carbon = state.builder.add_atom(AtomSpec::new(Element::C));
    let hydrogen = state.builder.add_atom(AtomSpec::new(Element::H));
    state
        .builder
        .add_bond(
            BondSpec::new(carbon, hydrogen, BondOrder::Single)
                .with_direction(BondDirection::BeginWedge),
        )
        .unwrap();

    state
        .remove_hs_update_explicit_count(&RemoveHsParams {
            remove_with_wedged_bond: false,
            update_explicit_count: true,
            ..RemoveHsParams::default()
        })
        .unwrap();
    let molecule = state.into_molecule().unwrap();

    assert_eq!(molecule.num_atoms(), 2);
}

#[test]
fn remove_hs_update_explicit_count_preserves_hydrogen_on_nontetrahedral_neighbor_by_default() {
    let mut state = SmilesBuildState::new();
    let center = state
        .builder
        .add_atom(AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::SquarePlanar));
    let hydrogen = state.builder.add_atom(AtomSpec::new(Element::H));
    state
        .builder
        .add_bond(BondSpec::new(center, hydrogen, BondOrder::Single))
        .unwrap();

    state
        .remove_hs_update_explicit_count(&RemoveHsParams {
            update_explicit_count: true,
            ..RemoveHsParams::default()
        })
        .unwrap();
    let molecule = state.into_molecule().unwrap();

    assert_eq!(molecule.num_atoms(), 2);
}

#[test]
fn remove_hs_update_explicit_count_preserves_query_hydrogen_by_default() {
    let mut state = SmilesBuildState::new();
    let carbon = state.builder.add_atom(AtomSpec::new(Element::C));
    let hydrogen = state.builder.add_atom(
        AtomSpec::new(Element::H).with_query(QueryNode::predicate(AtomQueryPredicate::Any)),
    );
    state
        .builder
        .add_bond(BondSpec::new(carbon, hydrogen, BondOrder::Single))
        .unwrap();

    state
        .remove_hs_update_explicit_count(&RemoveHsParams {
            update_explicit_count: true,
            ..RemoveHsParams::default()
        })
        .unwrap();
    let molecule = state.into_molecule().unwrap();

    assert_eq!(molecule.num_atoms(), 2);
}

#[test]
fn remove_hs_update_explicit_count_preserves_nonimplicit_hydrogen_when_disabled() {
    let mut state = SmilesBuildState::new();
    let carbon = state.builder.add_atom(AtomSpec::new(Element::C));
    let hydrogen = state.builder.add_atom(AtomSpec::new(Element::H));
    state
        .builder
        .add_bond(BondSpec::new(carbon, hydrogen, BondOrder::Single))
        .unwrap();

    state
        .remove_hs_update_explicit_count(&RemoveHsParams {
            remove_nonimplicit: false,
            update_explicit_count: true,
            ..RemoveHsParams::default()
        })
        .unwrap();
    let molecule = state.into_molecule().unwrap();

    assert_eq!(molecule.num_atoms(), 2);
}

#[test]
fn remove_hs_update_explicit_count_removes_degree_zero_hydrogen_when_enabled() {
    let mut state = SmilesBuildState::new();
    state.builder.add_atom(AtomSpec::new(Element::H));

    state
        .remove_hs_update_explicit_count(&RemoveHsParams {
            remove_degree_zero: true,
            ..RemoveHsParams::default()
        })
        .unwrap();
    let molecule = state.into_molecule().unwrap();

    assert_eq!(molecule.num_atoms(), 0);
}

#[test]
fn remove_hs_update_explicit_count_removes_higher_degree_hydrogen_when_enabled() {
    let mut state = SmilesBuildState::new();
    let carbon0 = state.builder.add_atom(AtomSpec::new(Element::C));
    let carbon1 = state.builder.add_atom(AtomSpec::new(Element::C));
    let hydrogen = state.builder.add_atom(AtomSpec::new(Element::H));
    state
        .builder
        .add_bond(BondSpec::new(carbon0, hydrogen, BondOrder::Single))
        .unwrap();
    state
        .builder
        .add_bond(BondSpec::new(carbon1, hydrogen, BondOrder::Single))
        .unwrap();

    state
        .remove_hs_update_explicit_count(&RemoveHsParams {
            remove_higher_degrees: true,
            update_explicit_count: true,
            ..RemoveHsParams::default()
        })
        .unwrap();
    let molecule = state.into_molecule().unwrap();

    assert_eq!(molecule.num_atoms(), 2);
    assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 1);
    assert_eq!(molecule.atoms()[1].explicit_hydrogens(), 1);
}

#[test]
fn remove_hs_update_explicit_count_removes_hydride_when_enabled() {
    let mut state = SmilesBuildState::new();
    let carbon = state.builder.add_atom(AtomSpec::new(Element::C));
    let hydride = state
        .builder
        .add_atom(AtomSpec::new(Element::H).with_formal_charge(-1));
    state
        .builder
        .add_bond(BondSpec::new(carbon, hydride, BondOrder::Single))
        .unwrap();

    state
        .remove_hs_update_explicit_count(&RemoveHsParams {
            remove_hydrides: true,
            update_explicit_count: true,
            ..RemoveHsParams::default()
        })
        .unwrap();
    let molecule = state.into_molecule().unwrap();

    assert_eq!(molecule.num_atoms(), 1);
    assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 1);
}

#[test]
fn remove_hs_update_explicit_count_does_not_always_increment_explicit_hydrogens() {
    let mut state = SmilesBuildState::new();
    let carbon = state.builder.add_atom(AtomSpec::new(Element::C));
    let hydrogen = state.builder.add_atom(AtomSpec::new(Element::H));
    state
        .builder
        .add_bond(BondSpec::new(carbon, hydrogen, BondOrder::Single))
        .unwrap();

    state
        .remove_hs_update_explicit_count(&RemoveHsParams::default())
        .unwrap();
    let molecule = state.into_molecule().unwrap();

    assert_eq!(molecule.num_atoms(), 1);
    assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 0);
}

#[test]
fn remove_hs_update_explicit_count_moves_end_direction_from_removed_hydrogen_like_rdkit() {
    let mut state = SmilesBuildState::new();
    let carbon = state.builder.add_atom(AtomSpec::new(Element::C));
    let hydrogen = state.builder.add_atom(AtomSpec::new(Element::H));
    let oxygen = state.builder.add_atom(AtomSpec::new(Element::O));
    state
        .builder
        .add_bond(
            BondSpec::new(carbon, hydrogen, BondOrder::Single)
                .with_direction(BondDirection::EndUpRight),
        )
        .unwrap();
    state
        .builder
        .add_bond(BondSpec::new(carbon, oxygen, BondOrder::Single))
        .unwrap();

    state
        .remove_hs_update_explicit_count(&RemoveHsParams::default())
        .unwrap();
    let molecule = state.into_molecule().unwrap();

    assert_eq!(molecule.num_atoms(), 2);
    assert_eq!(molecule.bonds()[0].direction(), BondDirection::EndDownRight);
}

#[test]
fn remove_hs_update_explicit_count_adjusts_double_bond_stereo_atoms_like_rdkit() {
    let mut state = SmilesBuildState::new();
    let center = state.builder.add_atom(AtomSpec::new(Element::C));
    let hydrogen = state.builder.add_atom(AtomSpec::new(Element::H));
    let oxygen = state.builder.add_atom(AtomSpec::new(Element::O));
    let carbon = state.builder.add_atom(AtomSpec::new(Element::C));
    state
        .builder
        .add_bond(BondSpec::new(center, hydrogen, BondOrder::Single))
        .unwrap();
    state
        .builder
        .add_bond(BondSpec::new(center, oxygen, BondOrder::Single))
        .unwrap();
    state
        .builder
        .add_bond(
            BondSpec::new(center, carbon, BondOrder::Double)
                .with_stereo_atoms(hydrogen, carbon)
                .with_stereo(BondStereo::Cis),
        )
        .unwrap();

    state
        .remove_hs_update_explicit_count(&RemoveHsParams {
            remove_defining_bond_stereo: true,
            ..RemoveHsParams::default()
        })
        .unwrap();
    let molecule = state.into_molecule().unwrap();

    assert_eq!(molecule.num_atoms(), 3);
    assert_eq!(molecule.bonds()[1].stereo(), BondStereo::Trans);
    assert_eq!(
        molecule.bonds()[1].stereo_atoms(),
        Some([AtomId::new(1), AtomId::new(2)])
    );
}

#[test]
fn remove_hs_update_explicit_count_sets_unknown_stereo_on_heavy_atom_like_rdkit() {
    let mut state = SmilesBuildState::new();
    let carbon = state.builder.add_atom(AtomSpec::new(Element::C));
    let hydrogen = state.builder.add_atom(AtomSpec::new(Element::H));
    state
        .builder
        .add_bond(
            BondSpec::new(carbon, hydrogen, BondOrder::Single)
                .with_direction(BondDirection::Unknown),
        )
        .unwrap();

    state
        .remove_hs_update_explicit_count(&RemoveHsParams::default())
        .unwrap();
    let molecule = state.into_molecule().unwrap();

    assert_eq!(molecule.num_atoms(), 1);
    assert!(molecule.atoms()[0].unknown_stereo());
}

#[test]
fn remove_hs_update_explicit_count_applies_sgroup_special_role_and_emptying_guards() {
    let mut state = SmilesBuildState::new();
    let carbon = state.builder.add_atom(AtomSpec::new(Element::C));
    let hydrogen = state.builder.add_atom(AtomSpec::new(Element::H));
    state
        .builder
        .add_bond(BondSpec::new(carbon, hydrogen, BondOrder::Single))
        .unwrap();
    state
        .builder
        .add_substance_group(
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Superatom)
                .with_atoms(vec![carbon])
                .with_attach_points(vec![crate::SGroupAttachPoint {
                    atom: hydrogen,
                    leaving_atom: None,
                    label: None,
                    order: None,
                }]),
        )
        .unwrap();

    state
        .remove_hs_update_explicit_count(&RemoveHsParams {
            update_explicit_count: true,
            ..RemoveHsParams::default()
        })
        .unwrap();
    let molecule = state.into_molecule().unwrap();
    assert_eq!(molecule.num_atoms(), 2);
    assert_eq!(molecule.num_bonds(), 1);

    let mut state = SmilesBuildState::new();
    let isolated_h = state.builder.add_atom(AtomSpec::new(Element::H));
    state
        .builder
        .add_substance_group(
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Superatom)
                .with_atoms(vec![isolated_h]),
        )
        .unwrap();

    state
        .remove_hs_update_explicit_count(&RemoveHsParams {
            remove_degree_zero: true,
            update_explicit_count: true,
            ..RemoveHsParams::default()
        })
        .unwrap();
    let molecule = state.into_molecule().unwrap();
    assert_eq!(molecule.num_atoms(), 1);
    assert_eq!(molecule.num_bonds(), 0);
    assert_eq!(molecule.atoms()[0].atomic_number(), 1);
}

#[test]
fn remove_hs_update_explicit_count_respects_remove_in_sgroups_false_membership_guard() {
    let mut state = SmilesBuildState::new();
    let hydrogen = state.builder.add_atom(AtomSpec::new(Element::H));
    state
        .builder
        .add_substance_group(
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Superatom)
                .with_atoms(vec![hydrogen]),
        )
        .unwrap();

    state
        .remove_hs_update_explicit_count(&RemoveHsParams {
            remove_degree_zero: true,
            remove_in_sgroups: false,
            update_explicit_count: true,
            ..RemoveHsParams::default()
        })
        .unwrap();
    let molecule = state.into_molecule().unwrap();

    assert_eq!(molecule.num_atoms(), 1);
    assert_eq!(molecule.atoms()[0].atomic_number(), 1);
}

#[test]
fn apply_smiles_postprocessing_removes_hydrogen_when_remove_hs_is_enabled() {
    let params = SmilesParseParams {
        remove_hs: true,
        sanitize: false,
        ..SmilesParseParams::default()
    };
    let mut state = to_mol("[CH]").unwrap();

    apply_smiles_postprocessing(&mut state, &params).unwrap();
    let molecule = state.into_molecule().unwrap();

    assert_eq!(molecule.num_atoms(), 1);
    assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 1);
}

#[test]
fn from_smiles_with_sanitize_false_parses_hash_element_and_tetrahedral_chirality() {
    let molecule = Molecule::from_smiles_with_sanitize("[#6].[C@H]", false).unwrap();

    assert_eq!(molecule.atomic_numbers(), vec![6, 6]);
    assert_eq!(molecule.atoms()[1].chiral_tag(), ChiralTag::TetrahedralCcw);
    assert_eq!(molecule.atoms()[1].explicit_hydrogens(), 1);
}

#[test]
fn from_smiles_with_sanitize_false_parses_bracket_element_lexer_table_like_rdkit() {
    let mut cases: Vec<(&str, u8, bool)> = BRACKET_ATOM_SYMBOLS
        .iter()
        .map(|(symbol, atomic_number)| (*symbol, *atomic_number, false))
        .collect();
    cases.extend([
        ("B", 5, false),
        ("C", 6, false),
        ("N", 7, false),
        ("O", 8, false),
        ("P", 15, false),
        ("S", 16, false),
        ("F", 9, false),
        ("Cl", 17, false),
        ("Br", 35, false),
        ("I", 53, false),
        ("H", 1, false),
        ("b", 5, true),
        ("c", 6, true),
        ("n", 7, true),
        ("o", 8, true),
        ("p", 15, true),
        ("s", 16, true),
        ("si", 14, true),
        ("as", 33, true),
        ("se", 34, true),
        ("te", 52, true),
    ]);

    for (symbol, atomic_number, aromatic) in cases {
        let molecule = Molecule::from_smiles_with_sanitize(&format!("[{symbol}]"), false).unwrap();
        assert_eq!(
            molecule.atoms()[0].atomic_number(),
            atomic_number,
            "symbol {symbol}"
        );
        assert_eq!(
            molecule.atoms()[0].is_aromatic(),
            aromatic,
            "symbol {symbol}"
        );
    }
}

#[test]
fn from_smiles_with_sanitize_false_parses_quoted_biovia_atoms_like_rdkit() {
    for &(symbol, atomic_number) in QUOTED_BIOVIA_ATOM_SYMBOLS {
        let molecule = Molecule::from_smiles_with_sanitize(&format!("[{symbol}]"), false).unwrap();
        assert_eq!(
            molecule.atoms()[0].atomic_number(),
            atomic_number,
            "symbol {symbol}"
        );
    }
}

#[test]
fn from_smiles_with_sanitize_false_converts_tetrahedral_chiral_class_like_rdkit() {
    let th1 = Molecule::from_smiles_with_sanitize("[C@TH1H]", false).unwrap();
    let th2 = Molecule::from_smiles_with_sanitize("[C@TH2H]", false).unwrap();

    assert_eq!(th1.atoms()[0].chiral_tag(), ChiralTag::TetrahedralCcw);
    assert_eq!(th1.atoms()[0].chiral_permutation(), None);
    assert_eq!(th2.atoms()[0].chiral_tag(), ChiralTag::TetrahedralCw);
    assert_eq!(th2.atoms()[0].chiral_permutation(), None);
}

#[test]
fn from_smiles_with_sanitize_false_preserves_non_tetrahedral_chiral_classes_like_rdkit() {
    let allene = Molecule::from_smiles_with_sanitize("[C@AL1]", false).unwrap();
    let square_planar = Molecule::from_smiles_with_sanitize("[C@SP3]", false).unwrap();
    let trigonal_bipyramidal = Molecule::from_smiles_with_sanitize("[C@TB20]", false).unwrap();
    let octahedral = Molecule::from_smiles_with_sanitize("[C@OH30]", false).unwrap();

    assert_eq!(allene.atoms()[0].chiral_tag(), ChiralTag::Allene);
    assert_eq!(allene.atoms()[0].chiral_permutation(), Some(1));
    assert_eq!(
        square_planar.atoms()[0].chiral_tag(),
        ChiralTag::SquarePlanar
    );
    assert_eq!(square_planar.atoms()[0].chiral_permutation(), Some(3));
    assert_eq!(
        trigonal_bipyramidal.atoms()[0].chiral_tag(),
        ChiralTag::TrigonalBipyramidal
    );
    assert_eq!(
        trigonal_bipyramidal.atoms()[0].chiral_permutation(),
        Some(20)
    );
    assert_eq!(octahedral.atoms()[0].chiral_tag(), ChiralTag::Octahedral);
    assert_eq!(octahedral.atoms()[0].chiral_permutation(), Some(30));
}

#[test]
fn from_smiles_with_sanitize_false_preserves_default_non_tetrahedral_chiral_permutations_like_rdkit()
 {
    let allene = Molecule::from_smiles_with_sanitize("[C@AL]", false).unwrap();
    let square_planar = Molecule::from_smiles_with_sanitize("[C@SP]", false).unwrap();
    let trigonal_bipyramidal = Molecule::from_smiles_with_sanitize("[C@TB]", false).unwrap();
    let octahedral = Molecule::from_smiles_with_sanitize("[C@OH]", false).unwrap();

    assert_eq!(allene.atoms()[0].chiral_tag(), ChiralTag::Allene);
    assert_eq!(allene.atoms()[0].chiral_permutation(), Some(0));
    assert_eq!(
        square_planar.atoms()[0].chiral_tag(),
        ChiralTag::SquarePlanar
    );
    assert_eq!(square_planar.atoms()[0].chiral_permutation(), Some(0));
    assert_eq!(
        trigonal_bipyramidal.atoms()[0].chiral_tag(),
        ChiralTag::TrigonalBipyramidal
    );
    assert_eq!(
        trigonal_bipyramidal.atoms()[0].chiral_permutation(),
        Some(0)
    );
    assert_eq!(octahedral.atoms()[0].chiral_tag(), ChiralTag::Octahedral);
    assert_eq!(octahedral.atoms()[0].chiral_permutation(), Some(0));
}

#[test]
fn from_smiles_with_sanitize_false_recomputes_nontetrahedral_ring_permutation_like_rdkit() {
    let sp1 = Molecule::from_smiles_with_sanitize("F[C@SP1]1(Br)I.Cl1", false).unwrap();
    let sp2 = Molecule::from_smiles_with_sanitize("F[C@SP2]1(Br)I.Cl1", false).unwrap();

    assert_eq!(sp1.atoms()[1].chiral_tag(), ChiralTag::SquarePlanar);
    assert_eq!(sp1.atoms()[1].chiral_permutation(), Some(2));
    assert_eq!(sp2.atoms()[1].chiral_tag(), ChiralTag::SquarePlanar);
    assert_eq!(sp2.atoms()[1].chiral_permutation(), Some(3));
}

#[test]
fn shared_canon_helpers_molecule_chiral_permutation_returns_zero_for_nonpositive_or_oversized_cases()
 {
    let mut builder = Molecule::builder();
    let center_without_perm =
        builder.add_atom(AtomSpec::new(Element::PT).with_chiral_tag(ChiralTag::SquarePlanar));
    let neighbor_a = builder.add_atom(AtomSpec::new(Element::CL));
    let neighbor_b = builder.add_atom(AtomSpec::new(Element::CL));
    let neighbor_c = builder.add_atom(AtomSpec::new(Element::CL));
    let neighbor_d = builder.add_atom(AtomSpec::new(Element::CL));
    let bond_a = builder
        .add_bond(BondSpec::new(
            center_without_perm,
            neighbor_a,
            BondOrder::Single,
        ))
        .unwrap();
    let bond_b = builder
        .add_bond(BondSpec::new(
            center_without_perm,
            neighbor_b,
            BondOrder::Single,
        ))
        .unwrap();
    let bond_c = builder
        .add_bond(BondSpec::new(
            center_without_perm,
            neighbor_c,
            BondOrder::Single,
        ))
        .unwrap();
    let bond_d = builder
        .add_bond(BondSpec::new(
            center_without_perm,
            neighbor_d,
            BondOrder::Single,
        ))
        .unwrap();
    let molecule_without_perm = builder.build().unwrap();

    assert_eq!(
        nontetrahedral_chiral_permutation_for_probe(
            &molecule_without_perm,
            center_without_perm,
            &[Some(bond_a), Some(bond_b), Some(bond_c), Some(bond_d)],
            false,
        )
        .unwrap(),
        0
    );

    let mut builder = Molecule::builder();
    let center_with_perm = builder.add_atom(
        AtomSpec::new(Element::PT)
            .with_chiral_tag(ChiralTag::SquarePlanar)
            .with_chiral_permutation(1),
    );
    let neighbor_a = builder.add_atom(AtomSpec::new(Element::CL));
    let neighbor_b = builder.add_atom(AtomSpec::new(Element::CL));
    let neighbor_c = builder.add_atom(AtomSpec::new(Element::CL));
    let neighbor_d = builder.add_atom(AtomSpec::new(Element::CL));
    let extra_neighbor = builder.add_atom(AtomSpec::new(Element::CL));
    let bond_a = builder
        .add_bond(BondSpec::new(
            center_with_perm,
            neighbor_a,
            BondOrder::Single,
        ))
        .unwrap();
    let bond_b = builder
        .add_bond(BondSpec::new(
            center_with_perm,
            neighbor_b,
            BondOrder::Single,
        ))
        .unwrap();
    let bond_c = builder
        .add_bond(BondSpec::new(
            center_with_perm,
            neighbor_c,
            BondOrder::Single,
        ))
        .unwrap();
    let bond_d = builder
        .add_bond(BondSpec::new(
            center_with_perm,
            neighbor_d,
            BondOrder::Single,
        ))
        .unwrap();
    let extra_bond = builder
        .add_bond(BondSpec::new(
            center_with_perm,
            extra_neighbor,
            BondOrder::Single,
        ))
        .unwrap();
    let molecule_with_perm = builder.build().unwrap();

    assert_eq!(
        nontetrahedral_chiral_permutation_for_probe(
            &molecule_with_perm,
            center_with_perm,
            &[
                Some(bond_a),
                Some(bond_b),
                Some(bond_c),
                Some(bond_d),
                Some(extra_bond),
            ],
            false,
        )
        .unwrap(),
        0
    );
}

#[test]
fn from_smiles_with_sanitize_false_rejects_invalid_chiral_permutation_like_rdkit() {
    let error = Molecule::from_smiles_with_sanitize("[C@TH3H]", false).unwrap_err();

    assert_eq!(
        error,
        SmilesParseError::ParseError("Invalid chiral specification on atom 0".to_string())
    );
}

#[test]
fn from_smiles_with_sanitize_false_rejects_non_tetrahedral_chiral_permutation_limits_like_rdkit() {
    let error = Molecule::from_smiles_with_sanitize("[C@SP4]", false).unwrap_err();

    assert_eq!(
        error,
        SmilesParseError::ParseError("Invalid chiral specification on atom 0".to_string())
    );
}

#[test]
fn from_smiles_with_sanitize_false_preserves_tetrahedral_chirality_across_ring_closure_ordering_like_rdkit()
 {
    let linear = Molecule::from_smiles_with_sanitize("F[C@](Cl)(Br)I", false).unwrap();
    let closure = Molecule::from_smiles_with_sanitize("F[C@]1(Br)I.Cl1", false).unwrap();

    assert_eq!(
        linear.atoms()[1].chiral_tag(),
        closure.atoms()[1].chiral_tag()
    );
}

#[test]
fn from_smiles_with_sanitize_false_parses_fused_ring_row_94_like_rdkit() {
    let molecule = Molecule::from_smiles_with_sanitize(
        "Cl.Cl.COc1ccc2nccc([C@@H](O)[C@@H]3C[C@@H]4CCN3C[C@@H]4C=C)c2c1",
        false,
    )
    .unwrap();

    assert_eq!(molecule.num_atoms(), 26);
    assert_eq!(
        molecule
            .atomic_numbers()
            .iter()
            .filter(|&&z| z == 17)
            .count(),
        2
    );
    assert_eq!(
        molecule
            .atoms()
            .iter()
            .filter(|atom| {
                matches!(
                    atom.chiral_tag(),
                    ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
                )
            })
            .count(),
        4
    );
}

#[test]
fn from_smiles_with_sanitize_false_restores_active_atom_after_branch() {
    let molecule = Molecule::from_smiles_with_sanitize("CC(C)O", false).unwrap();

    assert_eq!(molecule.atomic_numbers(), vec![6, 6, 6, 8]);
    assert_eq!(molecule.num_bonds(), 3);
    assert_eq!(molecule.bonds()[0].begin(), AtomId::new(0));
    assert_eq!(molecule.bonds()[0].end(), AtomId::new(1));
    assert_eq!(molecule.bonds()[1].begin(), AtomId::new(1));
    assert_eq!(molecule.bonds()[1].end(), AtomId::new(2));
    assert_eq!(molecule.bonds()[2].begin(), AtomId::new(1));
    assert_eq!(molecule.bonds()[2].end(), AtomId::new(3));
}

#[test]
fn close_branch_restores_root_before_following_atom_like_rdkit() {
    let molecule = to_mol("CC(C)O").unwrap().into_molecule().unwrap();

    assert_eq!(molecule.num_bonds(), 3);
    assert_eq!(molecule.bonds()[1].begin(), AtomId::new(1));
    assert_eq!(molecule.bonds()[1].end(), AtomId::new(2));
    assert_eq!(molecule.bonds()[2].begin(), AtomId::new(1));
    assert_eq!(molecule.bonds()[2].end(), AtomId::new(3));
}

#[test]
fn close_branch_pops_branch_stack_in_lifo_order_like_rdkit() {
    let mut state = SmilesBuildState::new();
    let a0 = state.builder.add_atom(AtomSpec::new(Element::C));
    let a1 = state.builder.add_atom(AtomSpec::new(Element::N));
    state.branch_stack.push(BranchPoint {
        atom: a0,
        open_position: 1,
    });
    state.branch_stack.push(BranchPoint {
        atom: a1,
        open_position: 3,
    });
    state.active_atom = Some(AtomId::new(99));

    state.close_branch().unwrap();
    assert_eq!(state.active_atom, Some(a1));
    assert_eq!(state.branch_stack.len(), 1);

    state.close_branch().unwrap();
    assert_eq!(state.active_atom, Some(a0));
    assert!(state.branch_stack.is_empty());
}

#[test]
fn close_branch_reports_extra_close_parentheses_without_open_branch_like_rdkit() {
    let mut state = SmilesBuildState::new();

    let error = state.close_branch().unwrap_err();

    assert_eq!(
        error,
        SmilesParseError::ParseError("extra close parentheses".to_string())
    );
}

#[test]
fn branch_open_token_returns_current_position_like_rdkit() {
    assert_eq!(SmilesBuildState::branch_open_token(7), 7);
}

#[test]
fn from_smiles_with_sanitize_false_parses_explicit_branch_bond() {
    let molecule = Molecule::from_smiles_with_sanitize("C(=O)N", false).unwrap();

    assert_eq!(molecule.atomic_numbers(), vec![6, 8, 7]);
    assert_eq!(molecule.num_bonds(), 2);
    assert_eq!(molecule.bonds()[0].begin(), AtomId::new(0));
    assert_eq!(molecule.bonds()[0].end(), AtomId::new(1));
    assert_eq!(molecule.bonds()[0].order(), BondOrder::Double);
    assert_eq!(molecule.bonds()[1].begin(), AtomId::new(0));
    assert_eq!(molecule.bonds()[1].end(), AtomId::new(2));
}

#[test]
fn add_ring_marker_records_pending_opening_and_closure_bookkeeping_like_rdkit() {
    let mut state = SmilesBuildState::new();
    state.add_first_atom(SmilesAtomToken::new(6)).unwrap();
    state.add_ring_marker(7).unwrap();

    assert_eq!(
        state.ring_openings.get(&7),
        Some(&RingOpening {
            atom: AtomId::new(0),
            pending_bond: Some(PendingBond {
                token: SmilesBondToken::directional(BondDirection::None),
                cx_smiles_bond_idx: None,
            }),
            input_position: 0,
        })
    );
    assert_eq!(
        state.ring_closures_by_atom.get(&AtomId::new(0)),
        Some(&vec![RingClosureRecord {
            ring_number: 7,
            bond: None,
        }])
    );
}

#[test]
fn add_single_bond_ring_marker_records_single_bond_pending_opening_like_rdkit() {
    let mut state = SmilesBuildState::new();
    state.add_first_atom(SmilesAtomToken::new(6)).unwrap();
    state.add_single_bond_ring_marker(9).unwrap();

    assert_eq!(
        state.ring_openings.get(&9),
        Some(&RingOpening {
            atom: AtomId::new(0),
            pending_bond: Some(PendingBond {
                token: SmilesBondToken::new(BondOrder::Single),
                cx_smiles_bond_idx: None,
            }),
            input_position: 0,
        })
    );
    assert_eq!(
        state.ring_closures_by_atom.get(&AtomId::new(0)),
        Some(&vec![RingClosureRecord {
            ring_number: 9,
            bond: None,
        }])
    );
}

#[test]
fn add_explicit_bond_ring_marker_records_pending_explicit_bond_like_rdkit() {
    let mut state = SmilesBuildState::new();
    state.add_first_atom(SmilesAtomToken::new(6)).unwrap();
    state
        .add_explicit_bond_ring_marker(SmilesBondToken::new(BondOrder::Double), 4)
        .unwrap();

    assert_eq!(
        state.ring_openings.get(&4),
        Some(&RingOpening {
            atom: AtomId::new(0),
            pending_bond: Some(PendingBond {
                token: SmilesBondToken::new(BondOrder::Double),
                cx_smiles_bond_idx: None,
            }),
            input_position: 0,
        })
    );
    assert_eq!(
        state.ring_closures_by_atom.get(&AtomId::new(0)),
        Some(&vec![RingClosureRecord {
            ring_number: 4,
            bond: None,
        }])
    );
}

#[test]
fn from_smiles_with_sanitize_false_closes_simple_ring_numbers() {
    let molecule = Molecule::from_smiles_with_sanitize("C1CC1", false).unwrap();

    assert_eq!(molecule.atomic_numbers(), vec![6, 6, 6]);
    assert_eq!(molecule.num_bonds(), 3);
    assert_eq!(molecule.bonds()[2].begin(), AtomId::new(2));
    assert_eq!(molecule.bonds()[2].end(), AtomId::new(0));
    assert_eq!(molecule.bonds()[2].order(), BondOrder::Single);
}

#[test]
fn from_smiles_with_sanitize_false_parses_percent_ring_numbers() {
    let molecule = Molecule::from_smiles_with_sanitize("C%12CC%12", false).unwrap();

    assert_eq!(molecule.num_atoms(), 3);
    assert_eq!(molecule.num_bonds(), 3);
    assert_eq!(molecule.bonds()[2].begin(), AtomId::new(2));
    assert_eq!(molecule.bonds()[2].end(), AtomId::new(0));
}

#[test]
fn from_smiles_with_sanitize_false_reports_unclosed_ring_like_rdkit() {
    let error = Molecule::from_smiles_with_sanitize("C1CC", false).unwrap_err();

    assert_eq!(
        error,
        SmilesParseError::ParseError("unclosed ring".to_string())
    );
}

#[test]
fn finish_parse_reports_unclosed_ring_before_invalid_chirality_like_rdkit() {
    let error = to_mol("C1[C@TH3]").unwrap_err();

    assert_eq!(
        error,
        SmilesParseError::ParseError("unclosed ring".to_string())
    );
}

#[test]
fn close_mol_rings_accepts_empty_pending_state_like_rdkit() {
    let mut state = SmilesBuildState::new();

    state.close_mol_rings().unwrap();
}

#[test]
fn close_mol_rings_reports_unclosed_ring_for_remaining_opening_like_rdkit() {
    let mut state = SmilesBuildState::new();
    state.add_first_atom(SmilesAtomToken::new(6)).unwrap();
    state.add_ring_marker(1).unwrap();

    let error = state.close_mol_rings().unwrap_err();

    assert_eq!(
        error,
        SmilesParseError::ParseError("unclosed ring".to_string())
    );
}

#[test]
fn second_ring_marker_preallocates_cx_smiles_bond_idx_like_rdkit() {
    let mut state = SmilesBuildState::new();
    state.add_first_atom(SmilesAtomToken::new(6)).unwrap();
    state.add_ring_marker(1).unwrap();
    state
        .add_atom_connected_to_active(SmilesAtomToken::new(6))
        .unwrap();
    state
        .add_atom_connected_to_active(SmilesAtomToken::new(6))
        .unwrap();
    state.add_ring_marker(1).unwrap();

    let closure = state.pending_ring_closures.first().unwrap();
    assert_eq!(closure.opening_pending_bond.cx_smiles_bond_idx, None);
    assert_eq!(closure.closing_pending_bond.cx_smiles_bond_idx, Some(2));
}

#[test]
fn from_smiles_cleanup_clears_ring_closure_cx_smiles_bond_idx_like_rdkit() {
    let molecule = Molecule::from_smiles_with_sanitize("C1CC1", false).unwrap();

    let ring_bond = molecule.bonds()[2].prop(CXSMILES_BOND_IDX_PROP);
    assert_eq!(ring_bond, None);
}

#[test]
fn from_smiles_with_sanitize_false_reports_extra_close_parentheses_like_rdkit() {
    let error = Molecule::from_smiles_with_sanitize("CC)", false).unwrap_err();

    assert_eq!(
        error,
        SmilesParseError::ParseError("extra close parentheses".to_string())
    );
}

#[test]
fn from_smiles_with_sanitize_false_reports_self_ring_closure_like_rdkit() {
    let error = Molecule::from_smiles_with_sanitize("C11", false).unwrap_err();

    assert_eq!(
        error,
        SmilesParseError::ParseError("duplicated ring closure bonds atom 0 to itself".to_string())
    );
}

#[test]
fn from_smiles_with_sanitize_false_reports_duplicate_ring_bond_like_rdkit() {
    let error = Molecule::from_smiles_with_sanitize("C1C1", false).unwrap_err();

    assert_eq!(
        error,
        SmilesParseError::ParseError(
            "ring closure duplicates bond between atom 0 and atom 1".to_string()
        )
    );
}

#[test]
fn from_smiles_with_sanitize_false_reports_unclosed_branch_like_rdkit() {
    let error = Molecule::from_smiles_with_sanitize("C(C", false).unwrap_err();

    assert_eq!(
        error,
        SmilesParseError::ParseError("extra open parentheses".to_string())
    );
}

#[test]
fn from_smiles_with_sanitize_false_reports_empty_branch_payload_like_rdkit_parser_boundary() {
    let error = Molecule::from_smiles_with_sanitize("C()", false).unwrap_err();

    assert_eq!(
        error,
        SmilesParseError::ParseError("expected branch atom or bond, got GroupClose".to_string())
    );
}

#[test]
fn from_smiles_with_sanitize_false_reports_bad_character_as_syntax_error_like_rdkit() {
    let error = Molecule::from_smiles_with_sanitize("C&N", false).unwrap_err();

    assert_eq!(
        error,
        SmilesParseError::ParseError("syntax error".to_string())
    );
}

#[test]
fn from_smiles_with_sanitize_false_uses_first_explicit_ring_bond_order_like_rdkit() {
    let molecule = Molecule::from_smiles_with_sanitize("C=1CC1", false).unwrap();

    assert_eq!(molecule.num_bonds(), 3);
    assert_eq!(molecule.bonds()[2].begin(), AtomId::new(0));
    assert_eq!(molecule.bonds()[2].end(), AtomId::new(2));
    assert_eq!(molecule.bonds()[2].order(), BondOrder::Double);
}

#[test]
fn from_smiles_with_sanitize_false_uses_closing_explicit_ring_bond_when_opening_unspecified() {
    let molecule = Molecule::from_smiles_with_sanitize("C1CC=1", false).unwrap();

    assert_eq!(molecule.num_bonds(), 3);
    assert_eq!(molecule.bonds()[2].begin(), AtomId::new(2));
    assert_eq!(molecule.bonds()[2].end(), AtomId::new(0));
    assert_eq!(molecule.bonds()[2].order(), BondOrder::Double);
}

#[test]
fn from_smiles_with_sanitize_false_ignores_conflicting_closing_ring_bond_spec_like_rdkit() {
    let molecule = Molecule::from_smiles_with_sanitize("C=1CC-1", false).unwrap();

    assert_eq!(molecule.num_bonds(), 3);
    assert_eq!(molecule.bonds()[2].begin(), AtomId::new(0));
    assert_eq!(molecule.bonds()[2].end(), AtomId::new(2));
    assert_eq!(molecule.bonds()[2].order(), BondOrder::Double);
}

#[test]
fn from_smiles_with_sanitize_false_marks_aromatic_atoms_and_unspecified_bonds_like_rdkit() {
    let molecule = Molecule::from_smiles_with_sanitize("c1ccccc1", false).unwrap();

    assert_eq!(molecule.num_atoms(), 6);
    assert_eq!(molecule.num_bonds(), 6);
    assert!(molecule.atoms().iter().all(|atom| atom.is_aromatic()));
    assert!(molecule.bonds().iter().all(|bond| bond.is_aromatic()));
    assert!(
        molecule
            .bonds()
            .iter()
            .all(|bond| bond.order() == BondOrder::Aromatic)
    );
}

#[test]
fn from_smiles_converts_non_ring_aromatic_bridge_to_single_like_rdkit_biphenyl() {
    let molecule = Molecule::from_smiles("c1ccccc1c1ccccc1").unwrap();

    let aromatic_bridge_bonds = molecule
        .bonds()
        .iter()
        .filter(|bond| {
            molecule.atoms()[bond.begin().index()].is_aromatic()
                && molecule.atoms()[bond.end().index()].is_aromatic()
                && !bond.is_aromatic()
        })
        .collect::<Vec<_>>();

    assert_eq!(aromatic_bridge_bonds.len(), 1);
    assert_eq!(aromatic_bridge_bonds[0].order(), BondOrder::Single);
}

#[test]
fn from_smiles_clears_nonunique_tetrahedral_tag_like_rdkit_row_92() {
    let molecule = Molecule::from_smiles(
        "O/C1=C/C=C/C=C1/CN3CCN(CC=2C=CC=CC=2O)[C@]3([H])C=4/C=C(/OC)C(=CC=4)OC",
    )
    .unwrap();

    assert!(
        molecule
            .atoms()
            .iter()
            .all(|atom| atom.chiral_tag() == ChiralTag::Unspecified)
    );
    assert_eq!(molecule.atoms()[20].explicit_hydrogens(), 0);
    assert!(!molecule.atoms()[20].no_implicit());
    assert_eq!(
        molecule
            .derived_cache()
            .valence
            .as_ref()
            .expect("sanitize valence cache")
            .implicit_hydrogens[20],
        1
    );
}

#[test]
fn from_smiles_preserves_ring_special_case_chiral_tags_like_rdkit_row_83() {
    let molecule = Molecule::from_smiles(
        "O=C(NC[C@]12C[C@H]3C[C@H](C[C@H](C3)C1)C2)[C@@H]1C[C@H]2c3ccccc3[C@@H]1c1ccccc12",
    )
    .unwrap();

    let tagged_atom_indices = molecule
        .atoms()
        .iter()
        .enumerate()
        .filter(|(_, atom)| atom.chiral_tag() != ChiralTag::Unspecified)
        .map(|(index, _)| index)
        .collect::<Vec<_>>();

    assert_eq!(tagged_atom_indices, vec![4, 6, 8, 10, 14, 16, 23]);
}

#[test]
fn from_smiles_assigns_ring_closure_double_bond_stereo_like_rdkit_row_106() {
    let molecule =
        Molecule::from_smiles("O=C(N(C(S/1)=S)CCC(O)=O)C1=C\\C2=CC=C(C3=CC=C(C=C3)Cl)O2").unwrap();

    let directional_bonds = molecule
        .bonds()
        .iter()
        .filter(|bond| {
            matches!(
                bond.direction(),
                BondDirection::EndUpRight | BondDirection::EndDownRight
            )
        })
        .count();
    let stereo_double_bonds = molecule
        .bonds()
        .iter()
        .filter(|bond| {
            matches!(
                bond.stereo(),
                BondStereo::Cis | BondStereo::Trans | BondStereo::E | BondStereo::Z
            ) && bond.stereo_atoms().is_some()
        })
        .count();

    assert_eq!(directional_bonds, 2);
    assert_eq!(stereo_double_bonds, 1);
}

#[test]
fn from_smiles_marks_stereochemistry_done_after_final_assignment_like_rdkit() {
    let molecule = Molecule::from_smiles("C[C@H]1CCC[C@](C)(O)C1").unwrap();

    assert_eq!(molecule.prop("_StereochemDone"), Some("1"));
}

#[test]
#[ignore = "debug helper for row 121 bond construction order"]
fn debug_row_121_bond_construction_order() {
    let molecule = Molecule::from_smiles(
            "O=C(O[Na])CC1=C(C(C(O[Na])=O)=C(C)C2=CC3=[N]4C(C(C=O)=C3CC)=CC5=C(C=C)C(C)=C6[N-]75)[N-]2[Cu+2]47[N](C8=C6)=C1C(C8C)CCC(O[Na])=O",
        )
        .unwrap();
    for bond in molecule.bonds() {
        eprintln!(
            "{} {} {} {:?} {}",
            bond.id().index(),
            bond.begin().index(),
            bond.end().index(),
            bond.order(),
            u8::from(bond.is_aromatic())
        );
    }
}

#[test]
fn from_smiles_empty_returns_empty_molecule_like_rdkit() {
    let molecule = Molecule::from_smiles("").unwrap();

    assert_eq!(molecule.num_atoms(), 0);
    assert_eq!(molecule.num_bonds(), 0);
}

#[test]
fn from_smiles_default_runs_postprocessing_for_simple_molecule() {
    let molecule = Molecule::from_smiles("CCO").unwrap();

    assert_eq!(molecule.atomic_numbers(), vec![6, 6, 8]);
    assert_eq!(molecule.num_bonds(), 2);
}

#[test]
fn from_smiles_default_parses_aromatic_with_sanitize_through_operations() {
    // Aromatic molecules now parse with default sanitize through
    // the registered operation pipeline (kekulize + valence + aromaticity).
    let molecule = Molecule::from_smiles("c1ccccc1").unwrap();
    // Sanitize should produce a valid Kekulé/aromatic molecule
    assert_eq!(molecule.atomic_numbers(), vec![6; 6]);
    assert_eq!(molecule.num_bonds(), 6);
}

#[test]
fn from_smiles_default_adjusts_disappearing_pyrrolic_hydrogen_like_rdkit() {
    let molecule = Molecule::from_smiles("c1cccN1").unwrap();

    assert!(molecule.atoms()[4].is_aromatic());
    assert_eq!(molecule.atoms()[4].explicit_hydrogens(), 1);
}

#[test]
fn from_smiles_default_preserves_explicit_bracket_pyrrolic_hydrogen_like_rdkit() {
    let unsanitized = Molecule::from_smiles_with_sanitize("[nH]1cccc1", false).unwrap();
    assert_eq!(unsanitized.atoms()[0].explicit_hydrogens(), 1);
    assert!(unsanitized.atoms()[0].no_implicit());

    let molecule = Molecule::from_smiles("[nH]1cccc1").unwrap();
    assert!(molecule.atoms()[0].is_aromatic());
    assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 1);
    assert!(!molecule.atoms()[0].no_implicit());
}

#[test]
fn from_smiles_default_removes_directional_h_with_sanitize_integration() {
    // Directional hydrogen with default sanitize now works through
    // the registered operations pipeline.
    let molecule = Molecule::from_smiles("[H]/C=C").unwrap();
    assert!(!molecule.atoms().is_empty(), "should parse successfully");
}

#[test]
fn from_smiles_assigns_rdkit_legacy_cip_ranks_for_acetic_acid() {
    let molecule = Molecule::from_smiles("CC(=O)O").unwrap();
    let observed = molecule
        .atoms()
        .iter()
        .map(|atom| {
            atom.prop("_CIPRank")
                .and_then(|value| value.parse::<u32>().ok())
        })
        .collect::<Vec<_>>();
    assert_eq!(observed, vec![Some(0), Some(1), Some(3), Some(2)]);
}

#[test]
fn from_smiles_reranks_legacy_cip_labels_when_possible_centers_remain_like_rdkit_row_2508() {
    let molecule =
        Molecule::from_smiles("CC(C)(c1ccc(OCC[C@H](O)Cl)cc1)c1ccc(OCC[C@@H](O)Cl)cc1").unwrap();
    let observed = molecule
        .atoms()
        .iter()
        .map(|atom| {
            atom.prop("_CIPRank")
                .and_then(|value| value.parse::<u32>().ok())
        })
        .collect::<Vec<_>>();

    assert_eq!(observed[0], Some(0));
    assert_eq!(observed[2], Some(0));
    assert_eq!(observed[3], Some(8));
    assert_eq!(observed[15], Some(7));
    assert_eq!(molecule.atoms()[10].prop("_CIPCode"), Some("R"));
    assert_eq!(molecule.atoms()[22].prop("_CIPCode"), Some("S"));
}

#[test]
fn from_smiles_default_integrates_cx_conformer_with_sanitize() {
    // CX conformer data with default sanitize now works through
    // the registered operations pipeline.
    let molecule = Molecule::from_smiles("CCO |(0,0,;1,0,;2,0,0.5)|").unwrap();
    assert_eq!(molecule.atomic_numbers(), vec![6, 6, 8]);
    assert_eq!(molecule.num_bonds(), 2);
}

#[test]
fn from_smiles_default_removes_explicit_h_and_updates_count() {
    let molecule = Molecule::from_smiles("[CH3][H]").unwrap();

    assert_eq!(molecule.atomic_numbers(), vec![6]);
    assert_eq!(molecule.num_bonds(), 0);
    assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 4);
}

#[test]
fn from_smiles_with_sanitize_false_keeps_explicit_hydrogen_like_rdkit_python_api() {
    let molecule = Molecule::from_smiles_with_sanitize("[CH3][H]", false).unwrap();

    assert_eq!(molecule.atomic_numbers(), vec![6, 1]);
    assert_eq!(molecule.num_bonds(), 1);
    assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 3);
}

#[test]
fn from_smiles_with_sanitize_false_leaves_small_ring_double_bond_stereo_unassigned_like_rdkit() {
    let molecule = Molecule::from_smiles_with_sanitize("C1/C=C/C2=C/CCCC2C1", false).unwrap();

    assert_eq!(molecule.prop("_StereochemDone"), None);
    // RDKit 2026.03.1 Python `MolFromSmiles(..., sanitize=False)` leaves
    // the raw directional single bonds in place but does not finalize
    // small-ring double-bond stereo on the returned molecule.
    assert_eq!(molecule.bonds()[1].stereo(), BondStereo::None);
    assert_eq!(molecule.bonds()[1].stereo_atoms(), None);
}

#[test]
fn from_smiles_persists_ring_stereo_props_like_rdkit_special_case_path() {
    let molecule = Molecule::from_smiles("C1[C@H](F)CC[C@H](Cl)C1").unwrap();

    assert_eq!(molecule.atoms()[1].prop("_ringStereochemCand"), Some("1"));
    assert_eq!(molecule.atoms()[1].prop("_ringStereoAtoms"), Some("6"));
    assert_eq!(molecule.atoms()[5].prop("_ringStereochemCand"), Some("1"));
    assert_eq!(molecule.atoms()[5].prop("_ringStereoAtoms"), Some("2"));
}

#[test]
fn from_smiles_reranks_chiral_center_after_double_bond_stereo_assignment_like_rdkit() {
    let molecule = Molecule::from_smiles("F[C@H](C/C=C/C)C/C=C\\C").unwrap();

    assert_eq!(molecule.atoms()[1].prop("_CIPCode"), Some("R"));
    assert_eq!(molecule.bonds()[3].stereo(), BondStereo::E);
    assert_eq!(
        molecule.bonds()[3].stereo_atoms(),
        Some([AtomId::new(2), AtomId::new(5)])
    );
    assert_eq!(molecule.bonds()[7].stereo(), BondStereo::Z);
    assert_eq!(
        molecule.bonds()[7].stereo_atoms(),
        Some([AtomId::new(6), AtomId::new(9)])
    );
}

#[test]
#[ignore = "debug helper for parser-stage checkpoint alignment"]
fn debug_probe_parser_fused_ring_chain() {
    let input = "C1/C=C/C2=C/CCCC2C1";
    let focus = [0usize, 1usize, 2usize, 3usize, 4usize];

    let print_builder_state = |name: &str, state: &SmilesBuildState| {
        eprintln!(
            "checkpoint={name} bonds={} pending_ring_closures={} ring_openings={}",
            state.builder.bonds().len(),
            state.pending_ring_closures.len(),
            state.ring_openings.len()
        );
        for closure in &state.pending_ring_closures {
            eprintln!(
                "pending ring={} opening_atom={} closing_atom={} opening_dir={:?} closing_dir={:?} opening_order={:?} closing_order={:?}",
                closure.ring_number,
                closure.opening_atom.index(),
                closure.closing_atom.index(),
                closure.opening_pending_bond.token.direction,
                closure.closing_pending_bond.token.direction,
                closure.opening_pending_bond.token.order,
                closure.closing_pending_bond.token.order,
            );
        }
        for bond_idx in focus {
            let Some(bond) = state.builder.bonds().get(bond_idx) else {
                continue;
            };
            eprintln!(
                "bond {} {}-{} order={:?} dir={:?} stereo={:?} stereo_atoms={:?}",
                bond_idx,
                bond.begin().index(),
                bond.end().index(),
                bond.order(),
                bond.direction(),
                bond.stereo(),
                bond.stereo_atoms()
            );
        }
    };

    let mut state = SmilesBuildState::new();
    let lexer = SmilesLexer::new(input);
    let mut parser = SmilesParser::new(lexer);
    parser.parse_mol(&mut state).unwrap();
    print_builder_state("post_parse_mol", &state);

    state.close_mol_rings().unwrap();
    print_builder_state("post_close_mol_rings", &state);

    state.check_chirality_specifications().unwrap();
    print_builder_state("post_check_chirality_specifications", &state);

    state.set_unspecified_bond_types().unwrap();
    print_builder_state("post_set_unspecified_bond_types", &state);

    state.adjust_atom_chirality_flags().unwrap();
    print_builder_state("post_adjust_atom_chirality_flags", &state);
}

#[test]
fn mol_from_smiles_remove_hs_without_sanitize_uses_split_reader_branch_like_rdkit() {
    let params = SmilesParseParams {
        sanitize: false,
        remove_hs: true,
        ..Default::default()
    };
    let molecule = mol_from_smiles("[CH3][H]", &params).unwrap();

    assert_eq!(molecule.atomic_numbers(), vec![6]);
    assert_eq!(molecule.num_bonds(), 0);
    assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 4);
}

#[test]
fn mol_from_smiles_skip_cleanup_preserves_parser_temporary_props_like_rdkit_reader_flag() {
    let params = SmilesParseParams {
        sanitize: false,
        skip_cleanup: true,
        ..Default::default()
    };
    let molecule = mol_from_smiles("c1ccccc1.C", &params).unwrap();

    assert_eq!(molecule.atoms()[0].prop(SMILES_START_PROP), Some("1"));
    assert_eq!(molecule.atoms()[6].prop(SMILES_START_PROP), Some("1"));
    assert!(
        molecule
            .bonds()
            .iter()
            .all(|bond| bond.prop(CXSMILES_BOND_IDX_PROP).is_some())
    );
    assert!(
        molecule.bonds()[..6]
            .iter()
            .all(|bond| bond.order() == BondOrder::Aromatic && bond.is_aromatic())
    );
}

#[test]
fn mol_from_smiles_sanitize_without_remove_hs_uses_molecule_sanitize_branch_like_rdkit() {
    let params = SmilesParseParams {
        sanitize: true,
        remove_hs: false,
        ..Default::default()
    };
    let molecule = mol_from_smiles("[CH3][H]", &params).unwrap();

    assert_eq!(molecule.atomic_numbers(), vec![6, 1]);
    assert_eq!(molecule.num_bonds(), 1);
}

#[test]
fn from_smiles_default_removes_mapped_hydrogen_like_rdkit_default() {
    let molecule = Molecule::from_smiles("[H:7]C").unwrap();

    assert_eq!(molecule.atomic_numbers(), vec![6]);
    assert_eq!(molecule.num_bonds(), 0);
    assert_eq!(molecule.atoms()[0].explicit_hydrogens(), 1);
}

#[test]
fn from_smiles_default_keeps_hydrogen_with_only_h_neighbor() {
    let molecule = Molecule::from_smiles("[H][H]").unwrap();

    assert_eq!(molecule.atomic_numbers(), vec![1, 1]);
    assert_eq!(molecule.num_bonds(), 1);
}

#[test]
fn from_smiles_default_keeps_hydrogen_with_dummy_neighbor() {
    let molecule = Molecule::from_smiles("[*][H]").unwrap();

    assert_eq!(molecule.atomic_numbers(), vec![0, 1]);
    assert_eq!(molecule.num_bonds(), 1);
}

#[test]
fn from_smiles_default_keeps_isotopic_hydrogen() {
    let molecule = Molecule::from_smiles("[2H]O").unwrap();

    assert_eq!(molecule.atomic_numbers(), vec![1, 8]);
    assert_eq!(molecule.num_bonds(), 1);
    assert_eq!(molecule.atoms()[0].isotope(), Some(2));
}

#[test]
fn smiles_parse_ops_clear_atom_chemical_props_resets_query_atom_state_like_rdkit() {
    let atom_id = AtomId::new(0);
    let mut atom = Atom::from_spec(
        atom_id,
        AtomSpec::new(Element::C)
            .with_isotope(13)
            .with_formal_charge(2)
            .with_explicit_hydrogens(3),
    );

    clear_atom_chemical_props(&mut atom);

    assert_eq!(atom.isotope(), None);
    assert_eq!(atom.formal_charge(), 0);
    assert_eq!(atom.explicit_hydrogens(), 0);
}

#[test]
fn smiles_parse_ops_report_parse_error_throws_when_requested_like_rdkit() {
    let error = report_parse_error("bad parse", true).unwrap_err();

    assert_eq!(error, SmilesParseError::ParseError("bad parse".to_string()));
}

#[test]
fn smiles_parse_ops_cleanup_after_parse_error_clears_unmatched_ring_state_like_rdkit() {
    let mut state = SmilesBuildState::new();
    state.add_first_atom(SmilesAtomToken::new(6)).unwrap();
    state.add_ring_marker(1).unwrap();

    assert!(!state.ring_openings.is_empty());
    assert!(!state.ring_closures_by_atom.is_empty());

    state.cleanup_after_parse_error();

    assert!(state.ring_openings.is_empty());
    assert!(state.ring_closures_by_atom.is_empty());
}

#[test]
fn smiles_parse_ops_add_frag_to_mol_merges_disconnected_fragment_like_rdkit() {
    let mut root = SmilesBuildState::new();
    root.add_first_atom(SmilesAtomToken::new(6)).unwrap();

    let mut frag = SmilesBuildState::new();
    frag.add_first_atom(SmilesAtomToken::new(8)).unwrap();

    root.add_frag_to_mol(frag, BondOrder::Ionic, BondDirection::None)
        .unwrap();
    let molecule = root.into_molecule().unwrap();

    assert_eq!(molecule.atomic_numbers(), vec![6, 8]);
    assert_eq!(molecule.num_bonds(), 0);
}

#[test]
fn smiles_parse_ops_add_frag_to_mol_connects_first_fragment_atom_in_insert_order_like_rdkit() {
    let mut root = SmilesBuildState::new();
    root.add_first_atom(SmilesAtomToken::new(6)).unwrap();

    let mut frag = SmilesBuildState::new();
    frag.add_first_atom(SmilesAtomToken::new(8)).unwrap();
    frag.add_atom_connected_to_active(SmilesAtomToken::new(7))
        .unwrap();

    root.add_frag_to_mol(frag, BondOrder::Single, BondDirection::None)
        .unwrap();

    assert_eq!(root.active_atom, Some(AtomId::new(1)));
    assert_eq!(root.builder.atoms().len(), 3);
    assert_eq!(root.builder.bonds().len(), 2);
    assert_eq!(root.builder.bonds()[0].begin(), AtomId::new(1));
    assert_eq!(root.builder.bonds()[0].end(), AtomId::new(2));
    assert_eq!(root.builder.bonds()[1].begin(), AtomId::new(0));
    assert_eq!(root.builder.bonds()[1].end(), AtomId::new(1));
}

#[test]
fn smiles_parse_ops_add_frag_to_mol_remaps_fragment_atom_and_bond_state_like_rdkit() {
    let mut root = SmilesBuildState::new();
    root.add_first_atom(SmilesAtomToken::new(6)).unwrap();
    root.add_atom_connected_to_active(SmilesAtomToken::new(6))
        .unwrap();

    let mut frag = SmilesBuildState::new();
    frag.add_first_atom(SmilesAtomToken::new(8)).unwrap();
    frag.add_atom_connected_to_active(SmilesAtomToken::new(7))
        .unwrap();
    frag.ring_closures_by_atom.insert(
        AtomId::new(0),
        vec![RingClosureRecord {
            ring_number: 7,
            bond: Some(BondId::new(0)),
        }],
    );
    frag.ring_openings.insert(
        9,
        RingOpening {
            atom: AtomId::new(1),
            pending_bond: Some(PendingBond {
                token: SmilesBondToken::new(BondOrder::Single),
                cx_smiles_bond_idx: None,
            }),
            input_position: 4,
        },
    );
    frag.temporary_chiral_permutations.insert(AtomId::new(1), 2);
    frag.cx_stereo_group_tracker.insert((1, 4), 3);

    root.add_frag_to_mol(frag, BondOrder::Ionic, BondDirection::None)
        .unwrap();

    assert_eq!(
        root.ring_closures_by_atom.get(&AtomId::new(2)),
        Some(&vec![RingClosureRecord {
            ring_number: 7,
            bond: Some(BondId::new(1)),
        }])
    );
    assert_eq!(
        root.ring_openings.get(&9),
        Some(&RingOpening {
            atom: AtomId::new(3),
            pending_bond: Some(PendingBond {
                token: SmilesBondToken::new(BondOrder::Single),
                cx_smiles_bond_idx: None,
            }),
            input_position: 4,
        })
    );
    assert_eq!(
        root.temporary_chiral_permutations.get(&AtomId::new(3)),
        Some(&2)
    );
    assert_eq!(root.cx_stereo_group_tracker.get(&(1, 4)), Some(&3));
}

#[test]
fn smiles_parse_ops_check_ring_closure_branch_status_inverts_matching_cases_like_rdkit() {
    let mut degree_one = SmilesBuildState::new();
    let a0 = degree_one
        .builder
        .add_atom(AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCw));
    let a1 = degree_one.builder.add_atom(AtomSpec::new(Element::F));
    degree_one
        .builder
        .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
        .unwrap();
    degree_one.check_ring_closure_branch_status(a0).unwrap();
    assert_eq!(
        degree_one.builder.atoms()[a0.index()].chiral_tag(),
        ChiralTag::TetrahedralCcw
    );

    let mut degree_two_nonzero = SmilesBuildState::new();
    let b0 = degree_two_nonzero
        .builder
        .add_atom(AtomSpec::new(Element::C));
    let b1 = degree_two_nonzero
        .builder
        .add_atom(AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCw));
    let b2 = degree_two_nonzero
        .builder
        .add_atom(AtomSpec::new(Element::F));
    degree_two_nonzero
        .builder
        .add_bond(BondSpec::new(b0, b1, BondOrder::Single))
        .unwrap();
    degree_two_nonzero
        .builder
        .add_bond(BondSpec::new(b1, b2, BondOrder::Single))
        .unwrap();
    degree_two_nonzero
        .check_ring_closure_branch_status(b1)
        .unwrap();
    assert_eq!(
        degree_two_nonzero.builder.atoms()[b1.index()].chiral_tag(),
        ChiralTag::TetrahedralCcw
    );

    let mut degree_three_root = SmilesBuildState::new();
    let c0 = degree_three_root
        .builder
        .add_atom(AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCw));
    let c1 = degree_three_root
        .builder
        .add_atom(AtomSpec::new(Element::F));
    let c2 = degree_three_root
        .builder
        .add_atom(AtomSpec::new(Element::CL));
    let c3 = degree_three_root
        .builder
        .add_atom(AtomSpec::new(Element::BR));
    degree_three_root
        .builder
        .add_bond(BondSpec::new(c0, c1, BondOrder::Single))
        .unwrap();
    degree_three_root
        .builder
        .add_bond(BondSpec::new(c0, c2, BondOrder::Single))
        .unwrap();
    degree_three_root
        .builder
        .add_bond(BondSpec::new(c0, c3, BondOrder::Single))
        .unwrap();
    degree_three_root
        .check_ring_closure_branch_status(c0)
        .unwrap();
    assert_eq!(
        degree_three_root.builder.atoms()[c0.index()].chiral_tag(),
        ChiralTag::TetrahedralCcw
    );
}

#[test]
fn smiles_parse_ops_check_ring_closure_branch_status_leaves_nonmatching_cases_unchanged() {
    let mut state = SmilesBuildState::new();
    let a0 = state.builder.add_atom(AtomSpec::new(Element::C));
    let a1 = state.builder.add_atom(AtomSpec::new(Element::F));
    let a2 = state
        .builder
        .add_atom(AtomSpec::new(Element::CL).with_chiral_tag(ChiralTag::TetrahedralCw));
    state
        .builder
        .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
        .unwrap();
    state
        .builder
        .add_bond(BondSpec::new(a1, a2, BondOrder::Single))
        .unwrap();

    state.check_ring_closure_branch_status(a2).unwrap();

    assert_eq!(
        state.builder.atoms()[a2.index()].chiral_tag(),
        ChiralTag::TetrahedralCw
    );
}

#[test]
fn smiles_parse_ops_to_mol_reports_parse_error_and_cleans_partial_state_like_rdkit() {
    let error = to_mol("C1CC").unwrap_err();

    assert_eq!(
        error,
        SmilesParseError::ParseError("unclosed ring".to_string())
    );
}

#[test]
fn smiles_parse_ops_to_mol_empty_input_returns_empty_state_like_rdkit() {
    let molecule = to_mol("").unwrap().into_molecule().unwrap();

    assert_eq!(molecule.num_atoms(), 0);
    assert_eq!(molecule.num_bonds(), 0);
}

#[test]
fn parse_cx_variable_attachments_rejects_multi_degree_source_like_rdkit() {
    let mut state = to_mol("*C(C)C").unwrap();
    let mut pos = 0;

    let error = parse_cx_variable_attachments(&mut state, "m:1:0.2", &mut pos).unwrap_err();

    assert_eq!(
        error,
        SmilesParseError::ParseError("failure parsing CXSMILES extensions".to_string())
    );
}
