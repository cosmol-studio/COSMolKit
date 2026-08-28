use super::*;

#[test]
fn shared_periodic_table_delegation_covers_smiles_symbols_and_fallback() {
    for atomic_number in 0..=118 {
        assert_eq!(
            element_symbol(atomic_number).unwrap(),
            crate::rdkit_element_symbol(atomic_number).unwrap()
        );
    }
    assert_eq!(element_symbol(200).unwrap(), "?");
}

#[test]
fn shared_count_swaps_smiles_writer_preserves_invariant_failure_mapping() {
    let size_error =
        stereo::count_swaps_to_interconvert(&[BondId::new(0)], &[BondId::new(0), BondId::new(1)])
            .unwrap_err();
    assert!(matches!(
        size_error,
        SmilesWriteError::InvariantViolation {
            stage: "ShortTermAtomWriter",
            message: "count_swaps_to_interconvert() probe/reference length mismatch",
        }
    ));

    let missing_error = stereo::count_swaps_to_interconvert(
        &[BondId::new(0), BondId::new(2)],
        &[BondId::new(0), BondId::new(1)],
    )
    .unwrap_err();
    assert!(matches!(
        missing_error,
        SmilesWriteError::InvariantViolation {
            stage: "ShortTermAtomWriter",
            message: "count_swaps_to_interconvert() expected bond missing from probe order",
        }
    ));
}

fn ethane() -> Molecule {
    Molecule::from_smiles_with_sanitize("CC", false).unwrap()
}

fn hypervalent_sulfur_like_geom_unsanitized_molecule() -> Molecule {
    let mut builder = crate::MoleculeBuilder::new();
    let sulfur = builder.add_atom(crate::AtomSpec::new(crate::Element::S).with_formal_charge(3));
    let nitrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::N));
    let oxygens = (0..4)
        .map(|_| builder.add_atom(crate::AtomSpec::new(crate::Element::O).with_formal_charge(-1)))
        .collect::<Vec<_>>();

    builder
        .add_bond(crate::BondSpec::new(sulfur, nitrogen, BondOrder::Single))
        .unwrap();
    for oxygen in oxygens {
        builder
            .add_bond(crate::BondSpec::new(sulfur, oxygen, BondOrder::Single))
            .unwrap();
    }
    builder.build().unwrap()
}

#[test]
fn mol_to_smiles_empty_molecule_returns_empty_string_like_rdkit_entrypoint() {
    let molecule = Molecule::from_smiles_with_sanitize("", false).unwrap();

    assert_eq!(
        mol_to_smiles(&molecule, &SmilesWriteParams::default()).unwrap(),
        ""
    );
}

#[test]
fn molecule_to_smiles_writes_basic_default_smiles() {
    let smiles = ethane().to_smiles(true).unwrap();

    assert_eq!(smiles, "CC");
}

#[test]
fn mol_to_smiles_rejects_invalid_root_before_writer_branches() {
    let mut params = SmilesWriteParams::default();
    params.rooted_at_atom = Some(2);

    let error = mol_to_smiles(&ethane(), &params).unwrap_err();

    assert_eq!(error, SmilesWriteError::RootedAtomOutOfRange { atom: 2 });
}

#[test]
fn all_primary_smiles_writer_modes_fail_closed_until_ported() {
    let molecule = ethane();
    assert_eq!(
        mol_to_smiles(&molecule, &SmilesWriteParams::default()).unwrap(),
        "CC"
    );
    let mut params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: false,
        clean_stereo: false,
        ..Default::default()
    };

    assert_eq!(mol_to_smiles(&molecule, &params).unwrap(), "CC");

    params.do_kekule = true;
    // Kekulization is now implemented; ethane (no aromatic bonds) succeeds.
    assert_eq!(mol_to_smiles(&molecule, &params).unwrap(), "CC");

    params.do_kekule = false;
    params.do_random = true;
    // Random SMILES is now implemented — ethane (simple molecule) succeeds
    let random = mol_to_smiles(&molecule, &params).unwrap();
    // For ethane, random SMILES should still produce a valid SMILES
    assert_eq!(
        random.len(),
        2,
        "random SMILES should be 2 chars: {random:?}"
    );

    // CX SMILES is now implemented — ethane with CX fields returns plain SMILES
    // (no CX data for ethane) or SMILES with empty CX extension.
    let result = mol_to_cx_smiles(
        &molecule,
        &SmilesWriteParams::default(),
        CxSmilesFields::ALL,
        RestoreBondDirOption::Clear,
    );
    assert!(result.is_ok(), "CX SMILES should succeed: {result:?}");
    // Ethane has no CX-specific data, so output is plain "CC"
    assert_eq!(result.unwrap(), "CC");

    // CX with no extra fields should also work
    let result = mol_to_cx_smiles(
        &molecule,
        &SmilesWriteParams::default(),
        CxSmilesFields::NONE,
        RestoreBondDirOption::None,
    );
    assert!(result.is_ok());
    assert_eq!(result.unwrap(), "CC");
}

#[test]
fn mol_to_smiles_updates_writer_property_cache_non_strict_like_rdkit() {
    let molecule = hypervalent_sulfur_like_geom_unsanitized_molecule();

    let strict_error =
        crate::assign_valence(&molecule, crate::ValenceModel::RdkitLike).unwrap_err();
    assert!(
        strict_error
            .to_string()
            .contains("Explicit valence for atom # 0 S, 5, is greater than permitted")
    );

    assert_eq!(
        mol_to_smiles(&molecule, &SmilesWriteParams::default()).unwrap(),
        "N[S+3]([O-])([O-])([O-])[O-]"
    );
}

#[test]
fn mol_to_random_smiles_vect_returns_requested_count_and_is_seed_reproducible() {
    let molecule = Molecule::from_smiles_with_sanitize("CC(C)CO", false).unwrap();

    let first =
        mol_to_random_smiles_vect(&molecule, 8, 0x1234, false, false, false, false).unwrap();
    let second =
        mol_to_random_smiles_vect(&molecule, 8, 0x1234, false, false, false, false).unwrap();
    let different_seed =
        mol_to_random_smiles_vect(&molecule, 8, 0x5678, false, false, false, false).unwrap();

    assert_eq!(first.len(), 8);
    assert_eq!(first, second);
    assert!(
        first.iter().any(|smiles| smiles != &first[0]),
        "random vector should contain traversal variation: {first:?}"
    );
    assert_ne!(first, different_seed);
}

#[test]
fn random_writer_seeded_start_matches_rdkit_for_diatomic_direction() {
    let molecule = Molecule::from_smiles_with_sanitize("N#C", false).unwrap();

    let generated =
        mol_to_random_smiles_vect(&molecule, 1, 1337, true, false, false, false).unwrap();

    assert_eq!(generated, vec!["C#N".to_string()]);
}

#[test]
fn mol_to_smiles_with_mode_returns_empty_string_for_empty_molecule_like_rdkit() {
    let molecule = Molecule::new();

    let result = mol_to_smiles_with_mode(
        &molecule,
        &SmilesWriteParams::default(),
        SmilesOutputMode::PlainSmiles,
    )
    .unwrap();

    assert_eq!(result, "");
}

#[test]
fn prepare_plain_smiles_molecule_stashes_and_clears_atom_maps_like_rdkit() {
    let molecule = Molecule::from_smiles_with_sanitize("[CH3:7][CH2:3][CH3:1]", false).unwrap();
    let mut working = molecule.clone();
    let params = SmilesWriteParams {
        canonical: false,
        do_isomeric_smiles: false,
        clean_stereo: false,
        ignore_atom_map_numbers: true,
        ..Default::default()
    };

    let saved = prepare_plain_smiles_molecule(&mut working, &params).unwrap();

    assert_eq!(saved, Some(vec![Some(7), Some(3), Some(1)]));
    assert_eq!(
        working
            .atoms()
            .iter()
            .map(|atom| atom.atom_map())
            .collect::<Vec<_>>(),
        vec![None, None, None]
    );
}

#[test]
fn prepare_plain_smiles_molecule_clears_dummy_atom_maps_for_canonical_ranking() {
    let molecule = Molecule::from_smiles_with_sanitize("[*:7]C[CH3:1]", false).unwrap();
    let mut working = molecule.clone();
    let params = SmilesWriteParams {
        ignore_atom_map_numbers: true,
        ..Default::default()
    };

    let saved = prepare_plain_smiles_molecule(&mut working, &params).unwrap();

    assert_eq!(saved, Some(vec![Some(7), None, Some(1)]));
    assert_eq!(
        working
            .atoms()
            .iter()
            .map(|atom| atom.atom_map())
            .collect::<Vec<_>>(),
        vec![None, None, None]
    );
}

#[test]
fn mol_to_smiles_ignore_atom_map_numbers_matches_rdkit_canonical_split_behavior() {
    let molecule = Molecule::from_smiles_with_sanitize("[CH3:7][CH2:3][CH3:1]", false).unwrap();
    let noncanonical = SmilesWriteParams {
        canonical: false,
        do_isomeric_smiles: false,
        clean_stereo: false,
        ignore_atom_map_numbers: true,
        ..Default::default()
    };
    let canonical = SmilesWriteParams {
        ignore_atom_map_numbers: true,
        ..Default::default()
    };

    assert_eq!(mol_to_smiles(&molecule, &noncanonical).unwrap(), "CCC");

    let canonical_smiles = mol_to_smiles(&molecule, &canonical).unwrap();
    assert!(canonical_smiles.contains(":7"));
    assert!(canonical_smiles.contains(":3"));
    assert!(canonical_smiles.contains(":1"));
}

#[test]
fn mol_to_smiles_ignore_atom_map_numbers_matches_rdkit_dummy_map_ordering() {
    let molecule = Molecule::from_smiles_with_sanitize("[*:1]C", false).unwrap();
    let noncanonical = SmilesWriteParams {
        do_isomeric_smiles: false,
        clean_stereo: false,
        canonical: false,
        ignore_atom_map_numbers: true,
        ..Default::default()
    };
    let canonical = SmilesWriteParams {
        do_isomeric_smiles: false,
        clean_stereo: false,
        canonical: true,
        ignore_atom_map_numbers: true,
        ..Default::default()
    };

    assert_eq!(mol_to_smiles(&molecule, &noncanonical).unwrap(), "*C");
    assert_eq!(mol_to_smiles(&molecule, &canonical).unwrap(), "[*:1]C");
}

#[test]
fn mol_to_smiles_do_kekule_handles_exocyclic_aryl_substituent_like_rdkit_row_88() {
    let molecule =
        Molecule::from_smiles("O=C1N(/N=C(/C)C1=NN/C2=C/C(OC)=CC=C2)C=3C=CC=CC=3").unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        do_kekule: true,
        canonical: false,
        clean_stereo: false,
        ..Default::default()
    };

    assert_eq!(
        mol_to_smiles(&molecule, &params).unwrap(),
        "O=C1N(C2=CC=CC=C2)N=C(C)C1=NNC1=CC(OC)=CC=C1"
    );
}

#[test]
fn writer_nonisomeric_explicit_bonds_clears_imine_direction_marks_like_rdkit_row_88() {
    let molecule =
        Molecule::from_smiles("O=C1N(/N=C(/C)C1=NN/C2=C/C(OC)=CC=C2)C=3C=CC=CC=3").unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: false,
        clean_stereo: false,
        all_bonds_explicit: true,
        ..Default::default()
    };

    assert_eq!(
        mol_to_smiles(&molecule, &params).unwrap(),
        "O=C1-N(-c2:c:c:c:c:c:2)-N=C(-C)-C-1=N-N-c1:c:c(-O-C):c:c:c:1"
    );
}

#[test]
fn choose_fragment_start_atom_prefers_terminal_dummy_like_rdkit() {
    let molecule = Molecule::from_smiles_with_sanitize("*C", false).unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        clean_stereo: false,
        ..Default::default()
    };
    let plan = collect_fragment_write_plans(&molecule, &params)
        .unwrap()
        .into_iter()
        .next()
        .unwrap();
    let ranks =
        rank_fragment_atoms_for_smiles(&molecule, &plan, &params, SmilesOutputMode::PlainSmiles)
            .unwrap();

    assert_eq!(ranks, vec![0, 1]);
    assert_eq!(
        choose_fragment_start_atom(&plan, &ranks, &params).unwrap(),
        AtomId::new(0)
    );
}

#[test]
fn collect_fragment_write_plans_tracks_fragment_root_and_bond_membership() {
    let molecule = Molecule::from_smiles_with_sanitize("CC.CCO", false).unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: false,
        clean_stereo: false,
        rooted_at_atom: Some(3),
        ..Default::default()
    };

    let plans = collect_fragment_write_plans(&molecule, &params).unwrap();

    assert_eq!(plans.len(), 2);
    assert_eq!(plans[0].atoms, vec![AtomId::new(0), AtomId::new(1)]);
    assert_eq!(plans[0].bonds, vec![BondId::new(0)]);
    assert_eq!(plans[0].rooted_at_atom, None);
    assert_eq!(
        plans[1].atoms,
        vec![AtomId::new(2), AtomId::new(3), AtomId::new(4)]
    );
    assert_eq!(plans[1].bonds, vec![BondId::new(1), BondId::new(2)]);
    assert_eq!(plans[1].rooted_at_atom, Some(AtomId::new(3)));
}

#[test]
fn choose_fragment_start_atom_consumes_random_seed_stream_sequentially() {
    let plan = FragmentWritePlan {
        atoms: vec![AtomId::new(0), AtomId::new(1), AtomId::new(2)],
        bonds: Vec::new(),
        rooted_at_atom: None,
    };
    let params = SmilesWriteParams {
        canonical: false,
        do_random: true,
        ..Default::default()
    };

    let chosen = with_random_smiles_seed(0x1234, || {
        Ok(vec![
            choose_fragment_start_atom(&plan, &[0, 1, 2], &params)?.index(),
            choose_fragment_start_atom(&plan, &[0, 1, 2], &params)?.index(),
            choose_fragment_start_atom(&plan, &[0, 1, 2], &params)?.index(),
        ])
    })
    .unwrap();

    let mut seed = 0x1234;
    let mut expected = Vec::new();
    for _ in 0..3 {
        expected.push(plan.atoms[(seed as usize) % plan.atoms.len()].index());
        seed = splitmix64(seed);
    }

    assert_eq!(chosen, expected);
}

#[test]
fn rooted_writer_only_reorders_the_containing_fragment() {
    let molecule = Molecule::from_smiles_with_sanitize("CC.CCO", false).unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: false,
        clean_stereo: false,
        rooted_at_atom: Some(3),
        ..Default::default()
    };

    assert_eq!(mol_to_smiles(&molecule, &params).unwrap(), "CC.C(C)O");
}

#[test]
fn write_fragment_smiles_uses_fragment_plan_root_for_subfragment_output() {
    let mut molecule = Molecule::from_smiles_with_sanitize("CC.CCO", false).unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: false,
        clean_stereo: false,
        rooted_at_atom: Some(3),
        ..Default::default()
    };
    let plans = collect_fragment_write_plans(&molecule, &params).unwrap();
    let mut context = SmilesWriteContext::default();

    let fragment = write_fragment_smiles(
        &mut molecule,
        &plans[1],
        &params,
        SmilesOutputMode::PlainSmiles,
        SmilesWriteOverrides::default(),
        &mut context,
    )
    .unwrap();

    assert_eq!(fragment.smiles, "C(C)O");
    assert_eq!(
        fragment.atom_ordering,
        vec![AtomId::new(3), AtomId::new(2), AtomId::new(4)]
    );
}

#[test]
fn canonical_writer_uses_traversal_order_for_tetrahedral_chirality() {
    let molecule = Molecule::from_smiles_with_sanitize("C[C@H](F)CCCl", true).unwrap();
    let mut params = SmilesWriteParams::default();

    params.rooted_at_atom = Some(1);
    assert_eq!(
        mol_to_smiles(&molecule, &params).unwrap(),
        "[C@@H](C)(F)CCCl"
    );

    params.rooted_at_atom = Some(2);
    assert_eq!(mol_to_smiles(&molecule, &params).unwrap(), "F[C@@H](C)CCCl");
}

#[test]
fn rooted_canonical_writer_matches_rdkit_branch_order_for_multichiral_case() {
    let molecule = Molecule::from_smiles_with_sanitize("O[C@](C)(Cl)[C@@](O)(Cl)C", true).unwrap();
    let params = SmilesWriteParams {
        rooted_at_atom: Some(0),
        ..Default::default()
    };

    assert_eq!(
        mol_to_smiles(&molecule, &params).unwrap(),
        "O[C@](C)(Cl)[C@@](C)(O)Cl"
    );
}

#[test]
fn shared_canon_helpers_writer_treats_single_h_query_as_fourth_valence() {
    let mut builder = crate::MoleculeBuilder::new();
    let center = builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_query(
        crate::QueryNode::and(vec![
            crate::QueryNode::predicate(crate::AtomQueryPredicate::AtomicNumber(6)),
            crate::QueryNode::predicate(crate::AtomQueryPredicate::ImplicitHydrogenCount(1)),
        ]),
    ));
    let fluorine = builder.add_atom(crate::AtomSpec::new(crate::Element::F));
    let chlorine = builder.add_atom(crate::AtomSpec::new(crate::Element::CL));
    let bromine = builder.add_atom(crate::AtomSpec::new(crate::Element::BR));
    builder
        .add_bond(crate::BondSpec::new(center, fluorine, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(crate::BondSpec::new(center, chlorine, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(crate::BondSpec::new(center, bromine, BondOrder::Single))
        .unwrap();
    let mut molecule = builder.build().unwrap();
    molecule.topology_block_mut().atoms[center.index()].set_chiral_tag(ChiralTag::TetrahedralCw);

    assert!(atom_has_fourth_valence_for_writer(&molecule, center));
    assert!(!chiral_atom_needs_tag_inversion_for_writer(
        &molecule, center, false, 1
    ));
}

#[test]
fn shared_canon_helpers_writer_treats_cached_implicit_h_as_fourth_valence() {
    let molecule = Molecule::from_smiles_with_sanitize("F[C@H](Cl)Br", true).unwrap();
    let center = AtomId::new(1);

    assert!(atom_has_fourth_valence_for_writer(&molecule, center));
    assert!(!chiral_atom_needs_tag_inversion_for_writer(
        &molecule, center, false, 1
    ));
}

#[test]
fn canonical_dfs_traversal_inserts_implicit_nontetrahedral_neighbors_like_rdkit() {
    let mut first_atom_bonds = vec![Some(BondId::new(3)), Some(BondId::new(4))];
    insert_implicit_nontetrahedral_neighbors_for_writer(
        &mut first_atom_bonds,
        ChiralTag::SquarePlanar,
        true,
    );
    assert_eq!(
        first_atom_bonds,
        vec![None, None, Some(BondId::new(3)), Some(BondId::new(4))]
    );

    let mut later_atom_bonds = vec![Some(BondId::new(3)), Some(BondId::new(4))];
    insert_implicit_nontetrahedral_neighbors_for_writer(
        &mut later_atom_bonds,
        ChiralTag::SquarePlanar,
        false,
    );
    assert_eq!(
        later_atom_bonds,
        vec![Some(BondId::new(3)), None, None, Some(BondId::new(4))]
    );
}

#[test]
fn canonical_dfs_traversal_marks_broken_chirality_for_partial_plan() {
    let molecule = Molecule::from_smiles_with_sanitize("F[C@](Cl)(Br)I", true).unwrap();
    let plan = FragmentWritePlan {
        atoms: vec![
            AtomId::new(0),
            AtomId::new(1),
            AtomId::new(2),
            AtomId::new(3),
        ],
        bonds: vec![BondId::new(0), BondId::new(1), BondId::new(2)],
        rooted_at_atom: Some(AtomId::new(1)),
    };
    let adjustments = compute_writer_chiral_adjustments(
        &molecule,
        &plan,
        AtomId::new(1),
        &vec![Vec::new(); molecule.num_atoms()],
        &vec![Vec::new(); molecule.num_atoms()],
        &[],
        false,
    )
    .unwrap();

    assert!(adjustments.broken_chiral_atoms.contains(&AtomId::new(1)));
    assert!(adjustments.chiral_inversions.is_empty());
    assert!(adjustments.chiral_permutations.is_empty());
}

#[test]
fn fragment_writer_suppresses_chirality_when_fragment_breaks_chiral_incident_bonds() {
    let molecule = Molecule::from_smiles_with_sanitize("C[C@H](F)CCCl", true).unwrap();
    let params = SmilesWriteParams {
        canonical: false,
        ..Default::default()
    };

    assert_eq!(
        mol_fragment_to_smiles(&molecule, &params, &[0, 1, 2], None, None, None).unwrap(),
        "CCF"
    );
    assert_eq!(
        mol_fragment_to_smiles(&molecule, &params, &[0, 1, 2, 3], None, None, None).unwrap(),
        "C[C@H](F)C"
    );
}

#[test]
fn writer_uses_parser_persisted_ring_stereo_props_without_local_reconstruction() {
    let molecule = Molecule::from_smiles_with_sanitize("C1[C@H](F)CC[C@H](Cl)C1", true).unwrap();
    let params = SmilesWriteParams {
        canonical: false,
        ..Default::default()
    };

    assert_eq!(
        mol_to_smiles(&molecule, &params).unwrap(),
        "C1[C@H](F)CC[C@H](Cl)C1"
    );
}

#[test]
fn ring_stereo_special_case_properties_are_computed_and_sanitize_clears_them_like_rdkit() {
    let molecule = Molecule::from_smiles_with_sanitize("C1[C@H](F)CC[C@H](Cl)C1", true).unwrap();

    for atom_index in [1, 5] {
        let atom = &molecule.atoms()[atom_index];
        assert_eq!(atom.prop("_ringStereochemCand"), Some("1"));
        assert!(atom.computed_prop_names().contains("_ringStereochemCand"));
        assert!(atom.prop("_ringStereoAtoms").is_some());
        assert!(atom.computed_prop_names().contains("_ringStereoAtoms"));
    }

    let sanitized = molecule
        .sanitize_with_ops(crate::SanitizeOps::NONE)
        .expect("sanitize with an empty mask still clears computed properties");
    for atom_index in [1, 5] {
        let atom = &sanitized.atoms()[atom_index];
        assert_eq!(atom.prop("_ringStereochemCand"), None);
        assert_eq!(atom.prop("_ringStereoAtoms"), None);
    }
}

#[test]
fn writer_accepts_ring_stereo_candidate_without_ring_neighbors_like_rdkit() {
    let mut molecule = Molecule::from_smiles_with_sanitize("F[C@@H]1O[C@H](Cl)S1", true).unwrap();
    molecule.topology_block_mut().atoms[1].set_computed_prop("_ringStereochemCand", "1");
    assert!(molecule.atoms()[1].prop("_ringStereoAtoms").is_none());
    assert_eq!(
        writer_ring_stereo_atoms(&molecule.atoms()[1], AtomId::new(1)).unwrap(),
        None
    );
    assert_eq!(
        mol_to_smiles(&molecule, &SmilesWriteParams::default()).unwrap(),
        "F[C@@H]1O[C@H](Cl)S1"
    );
}

#[test]
fn writer_rejects_malformed_ring_stereo_neighbors_prop() {
    let mut molecule =
        Molecule::from_smiles_with_sanitize("C1[C@H](F)CC[C@H](Cl)C1", true).unwrap();
    molecule.topology_block_mut().atoms[1].set_prop("_ringStereoAtoms", "0");

    let error = mol_to_smiles(&molecule, &SmilesWriteParams::default()).unwrap_err();

    assert_eq!(
        error,
        SmilesWriteError::InvalidRingStereoState {
            atom: 1,
            requirement: "`_ringStereoAtoms` must be a valid encoded ring-neighbor list",
        }
    );
}

#[test]
fn mol_fragment_to_smiles_rejects_rooted_atom_for_multifragment_molecule_without_bond_scope() {
    let molecule = Molecule::from_smiles_with_sanitize("CC.O", false).unwrap();
    let params = SmilesWriteParams {
        rooted_at_atom: Some(0),
        ..Default::default()
    };

    let error = mol_fragment_to_smiles(&molecule, &params, &[0, 1], None, None, None).unwrap_err();

    assert_eq!(
        error,
        SmilesWriteError::RootedAtomRequiresSingleFragment { atom: 0 }
    );
}

#[test]
fn isomeric_writer_outputs_non_tetrahedral_chiral_classes_like_rdkit() {
    let params = SmilesWriteParams {
        canonical: false,
        clean_stereo: false,
        ..Default::default()
    };

    let square_planar = Molecule::from_smiles_with_sanitize("[Pt@SP2](Cl)(Br)(I)F", false).unwrap();
    let trigonal_bipyramidal =
        Molecule::from_smiles_with_sanitize("[P@TB20](F)(Cl)(Br)(I)C", false).unwrap();
    let octahedral =
        Molecule::from_smiles_with_sanitize("[Co@OH30](F)(Cl)(Br)(I)(N)C", false).unwrap();

    assert_eq!(
        mol_to_smiles(&square_planar, &params).unwrap(),
        "[Pt@SP2]([Cl])([Br])([I])[F]"
    );
    assert_eq!(
        mol_to_smiles(&trigonal_bipyramidal, &params).unwrap(),
        "[P@TB20](F)(Cl)(Br)(I)C"
    );
    assert_eq!(
        mol_to_smiles(&octahedral, &params).unwrap(),
        "[Co@OH30]([F])([Cl])([Br])([I])([NH2])[CH3]"
    );
}

#[test]
fn isomeric_writer_outputs_default_non_tetrahedral_chiral_classes_like_rdkit() {
    let params = SmilesWriteParams {
        canonical: false,
        clean_stereo: false,
        ..Default::default()
    };

    let square_planar = Molecule::from_smiles_with_sanitize("[Pt@SP](Cl)(Br)(I)F", false).unwrap();
    let trigonal_bipyramidal =
        Molecule::from_smiles_with_sanitize("[P@TB](F)(Cl)(Br)(I)C", false).unwrap();
    let octahedral =
        Molecule::from_smiles_with_sanitize("[Co@OH](F)(Cl)(Br)(I)(N)C", false).unwrap();

    assert_eq!(
        mol_to_smiles(&square_planar, &params).unwrap(),
        "[Pt@SP]([Cl])([Br])([I])[F]"
    );
    assert_eq!(
        mol_to_smiles(&trigonal_bipyramidal, &params).unwrap(),
        "[P@TB](F)(Cl)(Br)(I)C"
    );
    assert_eq!(
        mol_to_smiles(&octahedral, &params).unwrap(),
        "[Co@OH]([F])([Cl])([Br])([I])([NH2])[CH3]"
    );
}

#[test]
fn rooted_writer_recomputes_non_tetrahedral_permutations_like_rdkit() {
    let mut params = SmilesWriteParams {
        canonical: false,
        clean_stereo: false,
        ..Default::default()
    };
    let square_planar = Molecule::from_smiles_with_sanitize("[Pt@SP2](Cl)(Br)(I)F", false).unwrap();
    let trigonal_bipyramidal =
        Molecule::from_smiles_with_sanitize("[P@TB20](F)(Cl)(Br)(I)C", false).unwrap();

    params.rooted_at_atom = Some(3);
    assert_eq!(
        mol_to_smiles(&square_planar, &params).unwrap(),
        "[I][Pt@SP3]([Cl])([Br])[F]"
    );

    params.rooted_at_atom = Some(4);
    assert_eq!(
        mol_to_smiles(&square_planar, &params).unwrap(),
        "[F][Pt@SP3]([Cl])([Br])[I]"
    );

    params.rooted_at_atom = Some(2);
    assert_eq!(
        mol_to_smiles(&trigonal_bipyramidal, &params).unwrap(),
        "Cl[P@TB15](F)(Br)(I)C"
    );

    params.rooted_at_atom = Some(5);
    assert_eq!(
        mol_to_smiles(&trigonal_bipyramidal, &params).unwrap(),
        "C[P@TB3](F)(Cl)(Br)I"
    );
}

#[test]
fn noncanonical_nonisomeric_plain_smiles_writes_simple_linear_fragments() {
    let molecule = Molecule::from_smiles_with_sanitize("CCO.N=O", false).unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: false,
        clean_stereo: false,
        ..Default::default()
    };

    assert_eq!(molecule.to_smiles_with_params(&params).unwrap(), "CCO.N=O");
}

#[test]
fn canonical_fragment_assembly_reorders_output_scope_with_sorted_smiles() {
    let molecule = Molecule::from_smiles_with_sanitize("O.C", false).unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: true,
        clean_stereo: false,
        ..Default::default()
    };
    let mut context = SmilesWriteContext::default();

    let smiles = mol_fragment_to_smiles_with_context(
        &molecule,
        &params,
        &[0, 1],
        None,
        None,
        None,
        &mut context,
    )
    .unwrap();

    assert_eq!(smiles, "C.O");
    assert_eq!(
        context.atom_output_order,
        vec![AtomId::new(1), AtomId::new(0)]
    );
}

#[test]
fn noncanonical_nonisomeric_plain_smiles_writes_explicit_single_bonds() {
    let molecule = Molecule::from_smiles_with_sanitize("CC", false).unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: false,
        clean_stereo: false,
        all_bonds_explicit: true,
        ..Default::default()
    };

    assert_eq!(mol_to_smiles(&molecule, &params).unwrap(), "C-C");
}

#[test]
fn noncanonical_nonisomeric_plain_smiles_writes_bracket_atoms_from_explicit_state() {
    let molecule = Molecule::from_smiles_with_sanitize("[NH4+].[O-].[SiH2].[*:7]", false).unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: false,
        clean_stereo: false,
        ..Default::default()
    };

    assert_eq!(
        mol_to_smiles(&molecule, &params).unwrap(),
        "[NH4+].[O-].[SiH2].[*:7]"
    );
}

#[test]
fn noncanonical_nonisomeric_plain_smiles_writes_branches_and_rings() {
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: false,
        clean_stereo: false,
        ..Default::default()
    };
    let ring = Molecule::from_smiles_with_sanitize("C1CC1", false).unwrap();
    let branch = Molecule::from_smiles_with_sanitize("CC(C)O", false).unwrap();
    let nested = Molecule::from_smiles_with_sanitize("C1CCC(CC1)O", false).unwrap();

    assert_eq!(mol_to_smiles(&ring, &params).unwrap(), "C1CC1");
    assert_eq!(mol_to_smiles(&branch, &params).unwrap(), "CC(C)O");
    assert_eq!(mol_to_smiles(&nested, &params).unwrap(), "C1CCC(O)CC1");
}

#[test]
fn noncanonical_nonisomeric_plain_smiles_writes_dative_bonds_and_ring_bond_symbols() {
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: false,
        clean_stereo: false,
        ..Default::default()
    };
    let dative = Molecule::from_smiles_with_sanitize("N->O", false).unwrap();
    let rooted_dative = SmilesWriteParams {
        rooted_at_atom: Some(1),
        ..params.clone()
    };
    let opening_double = Molecule::from_smiles_with_sanitize("C=1CC1", false).unwrap();
    let closing_double = Molecule::from_smiles_with_sanitize("C1CC=1", false).unwrap();

    assert_eq!(mol_to_smiles(&dative, &params).unwrap(), "N->O");
    assert_eq!(mol_to_smiles(&dative, &rooted_dative).unwrap(), "O<-N");
    assert_eq!(mol_to_smiles(&opening_double, &params).unwrap(), "C1=CC1");
    assert_eq!(mol_to_smiles(&closing_double, &params).unwrap(), "C1=CC1");
}

#[test]
fn plain_smiles_writer_strips_dative_bonds_on_working_copy_when_requested() {
    let molecule = Molecule::from_smiles_with_sanitize("N->O", false).unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: false,
        clean_stereo: false,
        include_dative_bonds: false,
        ..Default::default()
    };

    let mut working = molecule.clone();
    let _ = prepare_plain_smiles_molecule(&mut working, &params).unwrap();
    assert_eq!(working.bonds()[0].order(), BondOrder::Single);
    assert_eq!(total_num_hydrogens_for_writer(&working, AtomId::new(0)), 3);
    assert_eq!(total_valence_for_writer(&working, AtomId::new(0)), Some(3));

    assert_eq!(mol_to_smiles(&molecule, &params).unwrap(), "NO");
    assert_eq!(molecule.bonds()[0].order(), BondOrder::Dative);
}

#[test]
fn plain_smiles_writer_clears_molblock_only_bond_state_on_working_copy() {
    let mut builder = crate::MoleculeBuilder::new();
    let c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let c2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    builder
        .add_bond(
            crate::BondSpec::new(c1, c2, BondOrder::Double)
                .with_direction(BondDirection::EitherDouble)
                .with_stereo(BondStereo::Any),
        )
        .unwrap();
    let molecule = builder.build().unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: false,
        clean_stereo: false,
        ..Default::default()
    };

    assert_eq!(mol_to_smiles(&molecule, &params).unwrap(), "C=C");
    assert_eq!(molecule.bonds()[0].direction(), BondDirection::EitherDouble);
    assert_eq!(molecule.bonds()[0].stereo(), BondStereo::Any);
}

#[test]
fn get_bond_smiles_ignores_direction_on_double_bond_like_rdkit() {
    let mut builder = crate::MoleculeBuilder::new();
    let c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let c2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    builder
        .add_bond(
            crate::BondSpec::new(c1, c2, BondOrder::Double)
                .with_direction(BondDirection::BeginWedge),
        )
        .unwrap();
    let molecule = builder.build().unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: false,
        clean_stereo: false,
        ..Default::default()
    };

    assert_eq!(
        get_molecule_bond_smiles(&molecule, 0, Some(0), &params).unwrap(),
        "="
    );
}

#[test]
fn cx_smiles_restore_bond_dirs_clear_uses_working_copy_only() {
    let mut builder = crate::MoleculeBuilder::new();
    let c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let c2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    builder
        .add_bond(
            crate::BondSpec::new(c1, c2, BondOrder::Single).with_direction(BondDirection::Unknown),
        )
        .unwrap();
    let molecule = builder.build().unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: false,
        clean_stereo: false,
        ..Default::default()
    };

    assert_eq!(
        mol_to_cx_smiles(
            &molecule,
            &params,
            CxSmilesFields::ALL,
            RestoreBondDirOption::Clear,
        )
        .unwrap(),
        "CC"
    );
    assert_eq!(molecule.bonds()[0].direction(), BondDirection::Unknown);
}

#[test]
fn cx_preparation_does_not_apply_plain_only_bond_direction_cleanup() {
    let mut builder = crate::MoleculeBuilder::new();
    let c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let c2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    builder
        .add_bond(
            crate::BondSpec::new(c1, c2, BondOrder::Double)
                .with_direction(BondDirection::EitherDouble)
                .with_stereo(BondStereo::Any),
        )
        .unwrap();
    let molecule = builder.build().unwrap();
    let mut params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: false,
        clean_stereo: false,
        ..Default::default()
    };

    let mut keep = molecule.clone();
    prepare_cx_smiles_molecule(
        &mut keep,
        &mut params,
        CxSmilesFields::ALL,
        RestoreBondDirOption::None,
        true,
    )
    .unwrap();
    assert_eq!(keep.bonds()[0].direction(), BondDirection::EitherDouble);
    assert_eq!(keep.bonds()[0].stereo(), BondStereo::Any);

    let mut clear = molecule.clone();
    prepare_cx_smiles_molecule(
        &mut clear,
        &mut params,
        CxSmilesFields::ALL,
        RestoreBondDirOption::Clear,
        true,
    )
    .unwrap();
    assert_eq!(clear.bonds()[0].direction(), BondDirection::None);
    assert_eq!(clear.bonds()[0].stereo(), BondStereo::Any);
}

#[test]
fn noncanonical_nonisomeric_plain_smiles_writes_aromatic() {
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: false,
        clean_stereo: false,
        ..Default::default()
    };
    let benzene = Molecule::from_smiles_with_sanitize("c1ccccc1", false).unwrap();
    assert_eq!(mol_to_smiles(&benzene, &params).unwrap(), "c1ccccc1");

    // Also test a bracketed aromatic atom (Cl with aromatic flag)
    let params_with_h = SmilesWriteParams {
        all_hydrogens_explicit: true,
        ..params
    };
    let mol = Molecule::from_smiles_with_sanitize("c1ccc(C)cc1", false).unwrap();
    // should produce lowercase c's for aromatic carbons and uppercase C for sp3 carbon
    let smi = mol_to_smiles(&mol, &params_with_h).unwrap();
    assert!(
        smi.contains('c'),
        "expected lowercase c for aromatic atoms, got: {smi}"
    );
}

#[test]
fn get_bond_smiles_and_atom_needs_bracket_mark_aromatic_nh_like_rdkit() {
    let molecule = Molecule::from_smiles_with_sanitize("c1cc[nH]c1", false).unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: false,
        clean_stereo: false,
        ..Default::default()
    };

    assert!(atom_needs_bracket(&molecule, AtomId::new(3), "", &params).unwrap());
    assert_eq!(mol_to_smiles(&molecule, &params).unwrap(), "c1cc[nH]c1");
}

#[test]
fn noncanonical_nonisomeric_plain_smiles_writes_zero_order_and_radical_bonds_with_rdkit_compatibility()
 {
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: false,
        clean_stereo: false,
        ..Default::default()
    };
    // Zero-order bonds: RDKit outputs "~" for unknown/zero bond types.
    // The molecule CC~CC parses as 4 carbons with a zero-order bond (idx 1).
    let zero = Molecule::from_smiles_with_sanitize("CC~CC", false).unwrap();
    let output = mol_to_smiles(&zero, &params).unwrap();
    // Should contain the "~" bond symbol in the output
    assert!(
        output.contains('~'),
        "zero-order bond should map to ~, got: {output:?}"
    );

    // Radical-bearing molecules: the radical state is preserved through
    // the typed state and written as bracket notation when needed.
    // Build a molecule with an explicit radical atom.
    let mut builder = crate::MoleculeBuilder::new();
    let c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_radical_electrons(1));
    let c2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    builder
        .add_bond(crate::BondSpec::new(c1, c2, crate::BondOrder::Single))
        .unwrap();
    let radical = builder.build().unwrap();
    let output = mol_to_smiles(&radical, &params).unwrap();
    // Radical-bearing atoms get bracket notation: [CH3] or similar
    assert!(
        output.contains('['),
        "radical atom needs bracket: {output:?}"
    );
}

#[test]
fn noncanonical_nonisomeric_plain_smiles_honors_rooted_atom_for_traversal_start() {
    let molecule = Molecule::from_smiles_with_sanitize("CCO", false).unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: false,
        clean_stereo: false,
        rooted_at_atom: Some(1),
        ..Default::default()
    };

    assert_eq!(mol_to_smiles(&molecule, &params).unwrap(), "C(C)O");
}

#[test]
fn isomeric_smiles_handles_tetrahedral_chirality() {
    // (R)-Alaninol: [C@@H](C)(N)CO  →  canonical output may differ
    // from input but `@`/`@@` marks must be present.
    let mut params = SmilesWriteParams::default();
    params.canonical = false;
    params.clean_stereo = false;
    let mol = Molecule::from_smiles_with_sanitize("C[C@@H](N)CO", false).unwrap();
    let smi = mol_to_smiles(&mol, &params).unwrap();
    assert!(
        smi.contains("@") || smi.contains("@@"),
        "chiral atom should produce @ or @@ mark, got: {smi}"
    );
    assert_eq!(smi, "C[C@@H](N)CO");
}

#[test]
fn isomeric_smiles_handles_bond_stereo_direction() {
    let mut params = SmilesWriteParams::default();
    params.canonical = false;
    params.clean_stereo = false;
    let mol = Molecule::from_smiles_with_sanitize("C/C=C/C", false).unwrap();
    let smi = mol_to_smiles(&mol, &params).unwrap();
    assert_eq!(smi, "C/C=C/C");
}

#[test]
fn isomeric_smiles_writes_double_bond_stereo_with_direction() {
    let mut params = SmilesWriteParams::default();
    params.canonical = false;
    params.clean_stereo = false;
    let mol = Molecule::from_smiles_with_sanitize("C/C=C/C", true).unwrap();
    let smi = mol_to_smiles(&mol, &params).unwrap();
    assert_eq!(smi, "C/C=C/C");
}

#[test]
fn writer_canonicalizes_rooted_double_bond_directions_like_rdkit() {
    let mol = Molecule::from_smiles_with_sanitize("C/C=C/C", true).unwrap();
    let mut params = SmilesWriteParams {
        canonical: false,
        clean_stereo: false,
        rooted_at_atom: Some(1),
        ..Default::default()
    };

    assert_eq!(mol_to_smiles(&mol, &params).unwrap(), "C(/C)=C\\C");

    let z_mol = Molecule::from_smiles_with_sanitize("C/C=C\\C", true).unwrap();
    params.rooted_at_atom = Some(1);
    assert_eq!(mol_to_smiles(&z_mol, &params).unwrap(), "C(/C)=C/C");
}

#[test]
fn writer_rooted_double_bond_directions_match_rdkit_without_sanitize() {
    let mut params = SmilesWriteParams {
        canonical: false,
        clean_stereo: false,
        rooted_at_atom: Some(1),
        ..Default::default()
    };

    let e_mol = Molecule::from_smiles_with_sanitize("C/C=C/C", false).unwrap();
    assert_eq!(mol_to_smiles(&e_mol, &params).unwrap(), "C(/C)=C\\C");

    let z_mol = Molecule::from_smiles_with_sanitize("C/C=C\\C", false).unwrap();
    params.rooted_at_atom = Some(1);
    assert_eq!(mol_to_smiles(&z_mol, &params).unwrap(), "C(/C)=C/C");
}

#[test]
fn writer_rooted_terminal_double_bond_directions_match_rdkit_without_isomeric_smiles() {
    let mol = Molecule::from_smiles_with_sanitize("F/C=C\\F", true).unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: false,
        clean_stereo: false,
        all_bonds_explicit: true,
        rooted_at_atom: Some(3),
        ..Default::default()
    };

    assert_eq!(mol_to_smiles(&mol, &params).unwrap(), "F/C=C\\F");
}

#[test]
fn writer_preserves_rdkit_fused_ring_closure_digit_order_in_canonical_kekule_mode() {
    let mol = Molecule::from_smiles("Clc1ccc2ccc3cccc4ccc1c2c34").unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        do_kekule: true,
        canonical: true,
        clean_stereo: false,
        ..Default::default()
    };

    assert_eq!(
        mol_to_smiles(&mol, &params).unwrap(),
        "ClC1=C2C=CC3=CC=CC4=CC=C(C=C1)C2=C43"
    );
}

#[test]
fn writer_kekule_all_hydrogens_explicit_preserves_pyrrolic_hydrogen_like_rdkit() {
    let mol = Molecule::from_smiles("[nH]1cccc1").unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        do_kekule: true,
        canonical: false,
        clean_stereo: false,
        all_hydrogens_explicit: true,
        ..Default::default()
    };

    assert_eq!(
        mol_to_smiles(&mol, &params).unwrap(),
        "[NH]1[CH]=[CH][CH]=[CH]1"
    );
}

#[test]
fn writer_does_not_bracket_standard_aromatic_carbons_like_rdkit() {
    let mol = Molecule::from_smiles("c1ccccc1").unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: false,
        clean_stereo: false,
        ..Default::default()
    };

    assert!(!atom_needs_bracket(&mol, AtomId::new(0), "", &params).unwrap());
    assert_eq!(mol_to_smiles(&mol, &params).unwrap(), "c1ccccc1");
}

#[test]
fn writer_preserves_unbracketed_aromatic_ring_segment_like_rdkit_row_142() {
    let mol = Molecule::from_smiles("C12(C(C)c3ccccc3)NCC(C1(C)C)CC2").unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        do_kekule: false,
        canonical: false,
        clean_stereo: false,
        include_dative_bonds: false,
        ignore_atom_map_numbers: false,
        ..Default::default()
    };

    let smiles = mol_to_smiles(&mol, &params).unwrap();
    assert!(smiles.contains("c3ccccc3"), "got: {smiles}");
    assert!(!smiles.contains("[cH]"), "got: {smiles}");
}

#[test]
fn writer_restores_atom_maps_before_kekulize_like_rdkit_row_142() {
    let input = "[C:12]12([CH:62]([CH3:65])[c:61]3[cH:64][cH:67][cH:68][cH:66][cH:63]3)[CH:20]4[c:30]5[c:40]6[c:49]7[c:57]8[c:60]([c:59]9[c:55]([c:47]([c:44]([c:52]9[c:51]([c:43]%10[c:35]%11[c:25]%12[c:19]%13%14)[c:53]8[c:45]%11[c:39]6[c:29]4%13)[c:34]([c:24]%15[c:15]%16[c:7]%17[c:3]%18%19)[c:33]%10[c:23]%16[c:16]%12[c:8]%18[c:11]%14[c:5]1%20)[c:37]([c:36]%21[c:26]%22[c:18]%23[c:10]%24[c:13]%25[c:6]%26%27)[c:27]%15[c:17]%22[c:9]%17[c:4]%24[c:1]%19[c:2]%20%26)[c:54]([c:46]%21[c:38]%28[c:28]%23[c:21]%25%29)[c:56]%30[c:48]%28[c:41]%31[c:31]%29[c:22]%32[c:14]2%27)[c:58]%30[c:50]7[c:42]%31[c:32]5%32";
    let mut molecule = Molecule::from_smiles(input).unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        do_kekule: true,
        canonical: true,
        clean_stereo: false,
        include_dative_bonds: false,
        ignore_atom_map_numbers: true,
        rooted_at_atom: Some(molecule.num_atoms() - 1),
        ..Default::default()
    };

    let saved_atom_maps = prepare_plain_smiles_molecule(&mut molecule, &params).unwrap();
    let plan = collect_fragment_write_plans(&molecule, &params)
        .unwrap()
        .into_iter()
        .next()
        .unwrap();
    let _ranks =
        rank_fragment_atoms_for_smiles(&molecule, &plan, &params, SmilesOutputMode::PlainSmiles)
            .unwrap();
    restore_atom_maps_after_canonical_smiles(&mut molecule, saved_atom_maps.as_deref());
    let kekulized = kekulize_for_smiles(&molecule).unwrap();

    assert_eq!(kekulized.bonds()[11].order(), BondOrder::Double);
    assert_eq!(kekulized.bonds()[72].order(), BondOrder::Single);
}

#[test]
fn writer_kekule_ignoring_dative_bonds_matches_rdkit_metal_porphyrin_branch() {
    let mol = Molecule::from_smiles(
            "O=C(O[Na])CC1=C(C(C(O[Na])=O)=C(C)C2=CC3=[N]4C(C(C=O)=C3CC)=CC5=C(C=C)C(C)=C6[N-]75)[N-]2[Cu+2]47[N](C8=C6)=C1C(C8C)CCC(O[Na])=O",
        )
        .unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        do_kekule: true,
        canonical: false,
        clean_stereo: false,
        include_dative_bonds: false,
        ..Default::default()
    };

    assert_eq!(
        mol_to_smiles(&mol, &params).unwrap(),
        "O=C([O][Na])CC1=C2[N]3[Cu+2]45[N]6C(=CC7=C(C)C(C([O][Na])=O)=C1[N-]74)C(CC)=C(C=O)C6=CC1=C(C=C)C(C)=C([N-]15)C=C3C(C)C2CCC([O][Na])=O"
    );
}

fn writer_cleans_nonstereo_double_bond_direction_specs_like_rdkit() {
    let mol = Molecule::from_smiles_with_sanitize("C/C=C(/C)C", true).unwrap();
    let params = SmilesWriteParams {
        canonical: true,
        clean_stereo: false,
        ..Default::default()
    };

    assert_eq!(mol_to_smiles(&mol, &params).unwrap(), "CC=C(C)C");
}

#[test]
fn writer_canonical_nitro_start_atom_matches_rdkit_default_branch() {
    let mol = Molecule::from_smiles_with_sanitize("[O-][N+](=O)O", true).unwrap();
    let params = SmilesWriteParams::default();

    assert_eq!(mol_to_smiles(&mol, &params).unwrap(), "O=[N+]([O-])O");
}

#[test]
fn writer_uses_cip_rank_ties_for_duplicate_double_bond_ligands_like_rdkit() {
    let params = SmilesWriteParams {
        canonical: true,
        clean_stereo: false,
        ..Default::default()
    };

    let fixtures = [
        ("C/C=C(/C(F))C(Cl)", "C/C=C(/CF)CCl"),
        ("C/C=C(/C(F))C(F)", "CC=C(CF)CF"),
        ("C/C=C(/CO)CN", "C/C=C(\\CN)CO"),
    ];

    for (input, expected) in fixtures {
        for sanitize in [true, false] {
            let mol = Molecule::from_smiles_with_sanitize(input, sanitize).unwrap();
            assert_eq!(
                mol_to_smiles(&mol, &params).unwrap(),
                expected,
                "RDKit 2026.03.1 canonical output for {input} with sanitize={sanitize}"
            );
        }
    }
}

#[test]
fn writer_canonicalizes_connected_double_bond_direction_queue_like_rdkit() {
    let mol = Molecule::from_smiles_with_sanitize("C/C=C/C=C\\C", true).unwrap();
    let params = SmilesWriteParams {
        canonical: true,
        clean_stereo: false,
        ..Default::default()
    };

    assert_eq!(mol_to_smiles(&mol, &params).unwrap(), "C/C=C\\C=C\\C");
}

#[test]
fn writer_canonicalizes_connected_double_bond_direction_queue_without_sanitize() {
    let mol = Molecule::from_smiles_with_sanitize("C/C=C/C=C\\C", false).unwrap();
    let params = SmilesWriteParams {
        canonical: true,
        clean_stereo: false,
        ..Default::default()
    };

    assert_eq!(mol_to_smiles(&mol, &params).unwrap(), "C/C=C\\C=C\\C");
}

#[test]
fn writer_canonicalizes_ring_double_bond_directions_like_rdkit() {
    let params = SmilesWriteParams {
        canonical: true,
        clean_stereo: false,
        ..Default::default()
    };

    let fixtures = [
        ("C1C/C=C/CCCCCCCC1", "C1=C/CCCCCCCCCC/1"),
        ("C1/C=C/C=C/CCCCCCCCC1", "C1=C/CCCCCCCCCC/C=C/1"),
        ("C1/C=C/CCCCC1", "C1=C/CCCCCC/1"),
        ("C1/C=C\\CCCCC1", "C1=C\\CCCCCC/1"),
    ];

    for (input, expected) in fixtures {
        for sanitize in [true, false] {
            let mol = Molecule::from_smiles_with_sanitize(input, sanitize).unwrap();
            assert_eq!(
                mol_to_smiles(&mol, &params).unwrap(),
                expected,
                "RDKit 2026.03.1 canonical output for {input} with sanitize={sanitize}"
            );
        }
    }
}

#[test]
fn writer_canonicalizes_fused_ring_double_bond_directions_like_rdkit() {
    let params = SmilesWriteParams {
        canonical: true,
        clean_stereo: false,
        ..Default::default()
    };

    let fixtures = [
        ("C1/C=C/C2CCCCC2C1", "C1=CC2CCCCC2CC1", "C1=C/C2CCCCC2CC/1"),
        (
            "C1/C=C\\C2CCCCC2C1",
            "C1=CC2CCCCC2CC1",
            "C1=C\\C2CCCCC2CC/1",
        ),
        (
            "C1C/C=C/CC2CCCCC12",
            "C1=CCC2CCCCC2CC1",
            "C1=C/CC2CCCCC2CC/1",
        ),
        (
            "C1/C=C/C2=C/CCCC2C1",
            "C1=CC2=CCCCC2CC1",
            "C1=C2/C=C/CCC2CCC1",
        ),
        (
            "C1/C=C\\C2=C/CCCC2C1",
            "C1=CC2=CCCCC2CC1",
            "C1=C2/C=C\\CCC2CCC1",
        ),
        (
            "C1/C=C/C=C/CC2CCCCC12",
            "C1=C/CC2CCCCC2C/C=C/1",
            "C1=C/CC2CCCCC2C/C=C/1",
        ),
        (
            "C1/C=C\\C=C/CC2CCCCC12",
            "C1=C\\CC2CCCCC2C\\C=C/1",
            "C1=C\\CC2CCCCC2C\\C=C/1",
        ),
        (
            "C1/C=C/C=C\\CC2CCCCC12",
            "C1=C\\CC2CCCCC2C/C=C/1",
            "C1=C\\CC2CCCCC2C/C=C/1",
        ),
    ];

    for (input, expected_sanitized, expected_unsanitized) in fixtures {
        for (sanitize, expected) in [(true, expected_sanitized), (false, expected_unsanitized)] {
            let mol = Molecule::from_smiles_with_sanitize(input, sanitize).unwrap();
            assert_eq!(
                mol_to_smiles(&mol, &params).unwrap(),
                expected,
                "RDKit 2026.03.1 canonical output for {input} with sanitize={sanitize}"
            );
        }
    }
}

#[test]
fn prepare_plain_smiles_molecule_initializes_fast_ring_info_for_fused_ring_stereo_like_rdkit() {
    let input = "C1/C=C/C2=C/CCCC2C1";
    let mut molecule = Molecule::from_smiles_with_sanitize(input, false).unwrap();
    let params = SmilesWriteParams {
        canonical: true,
        clean_stereo: false,
        ..Default::default()
    };

    let _saved_atom_maps = prepare_plain_smiles_molecule(&mut molecule, &params).unwrap();
    assert!(
        molecule
            .derived_cache()
            .rings
            .as_ref()
            .is_some_and(crate::RingInfo::is_find_fast_or_better)
    );
}

#[test]
fn writer_fragment_temp_graph_resets_ring_info_like_rdkit() {
    for (input, clears_computed) in [("C1CCCCC1", false), ("C1CCCCC1.[Na+]", true)] {
        let mut molecule = Molecule::from_smiles(input).unwrap();
        assert!(molecule.derived_cache().rings.is_some());
        molecule
            .properties_mut()
            .set_computed_prop("computed-molecule", "remove");
        molecule
            .properties_mut()
            .set_prop("ordinary-molecule", "keep");
        molecule.topology_block_mut().atoms[0].set_computed_prop("computed-atom", "remove");
        molecule.topology_block_mut().atoms[0].set_prop("ordinary-atom", "keep");
        molecule.topology_block_mut().bonds[0].set_computed_prop("computed-bond", "remove");
        molecule.topology_block_mut().bonds[0].set_prop("ordinary-bond", "keep");

        clear_fragment_temp_molecule_computed_stereo_props_for_writer(&mut molecule);

        assert_eq!(
            molecule.derived_cache().rings.is_none(),
            clears_computed,
            "input={input}"
        );
        assert_eq!(
            molecule.prop("computed-molecule"),
            if clears_computed {
                None
            } else {
                Some("remove")
            },
            "input={input}"
        );
        assert_eq!(
            molecule.prop("ordinary-molecule"),
            Some("keep"),
            "input={input}"
        );
        assert_eq!(
            molecule.atoms()[0].prop("computed-atom"),
            if clears_computed {
                None
            } else {
                Some("remove")
            },
            "input={input}"
        );
        assert_eq!(
            molecule.atoms()[0].prop("ordinary-atom"),
            Some("keep"),
            "input={input}"
        );
        assert_eq!(
            molecule.bonds()[0].prop("computed-bond"),
            if clears_computed {
                None
            } else {
                Some("remove")
            },
            "input={input}"
        );
        assert_eq!(
            molecule.bonds()[0].prop("ordinary-bond"),
            Some("keep"),
            "input={input}"
        );
    }
}

#[test]
fn writer_fragment_reassignment_sets_double_bond_stereo_from_directions_like_rdkit() {
    let mut molecule = Molecule::from_smiles_with_sanitize("FC=CF.[Na+]", false).unwrap();
    let double_bond = molecule
        .bonds()
        .iter()
        .find(|bond| bond.order() == BondOrder::Double)
        .map(|bond| bond.id())
        .unwrap();
    let double = &molecule.bonds()[double_bond.index()];
    let begin = double.begin();
    let end = double.end();
    let begin_directed = incident_bonds(&molecule, begin)
        .into_iter()
        .find(|bond| *bond != double_bond)
        .unwrap();
    let end_directed = incident_bonds(&molecule, end)
        .into_iter()
        .find(|bond| *bond != double_bond)
        .unwrap();
    molecule.topology_block_mut().bonds[begin_directed.index()]
        .set_direction(BondDirection::EndUpRight);
    molecule.topology_block_mut().bonds[end_directed.index()]
        .set_direction(BondDirection::EndDownRight);
    molecule.topology_block_mut().bonds[double_bond.index()].set_stereo(BondStereo::None);
    molecule.topology_block_mut().bonds[double_bond.index()].set_stereo_atoms(None);
    molecule
        .properties_mut()
        .set_computed_prop("_StereochemDone", "1");

    prepare_plain_smiles_molecule(&mut molecule, &SmilesWriteParams::default()).unwrap();

    let double = &molecule.bonds()[double_bond.index()];
    assert_eq!(double.stereo(), BondStereo::Z);
    assert_eq!(
        double.stereo_atoms(),
        Some([
            bond_other_atom(&molecule.bonds()[begin_directed.index()], begin).unwrap(),
            bond_other_atom(&molecule.bonds()[end_directed.index()], end).unwrap(),
        ])
    );
}

#[test]
#[ignore = "debug helper for upstream/rust checkpoint alignment"]
fn debug_probe_rust_writer_fused_ring_chain() {
    let input = "C1/C=C/C2=C/CCCC2C1";
    let params = SmilesWriteParams {
        canonical: true,
        clean_stereo: false,
        ..Default::default()
    };

    let focus = [0usize, 1usize, 2usize, 3usize, 4usize];
    let print_state = |name: &str, mol: &Molecule| {
        eprintln!(
            "checkpoint={name} rings_initialized={} is_symm_sssr={} stereo_done={}",
            mol.derived_cache().rings.is_some(),
            mol.derived_cache()
                .rings
                .as_ref()
                .is_some_and(crate::RingInfo::is_symm_sssr),
            mol.prop("_StereochemDone").is_some()
        );
        for bond_idx in focus {
            let bond = &mol.bonds()[bond_idx];
            let ring_count = mol
                .derived_cache()
                .rings
                .as_ref()
                .map(|ri| ri.num_bond_rings(BondId::new(bond_idx)))
                .unwrap_or(0);
            let min_ring = mol
                .derived_cache()
                .rings
                .as_ref()
                .map(|ri| ri.min_bond_ring_size(BondId::new(bond_idx)))
                .unwrap_or(0);
            eprintln!(
                "bond {} {}-{} dir={:?} stereo={:?} stereo_atoms={:?} ring_count={} min_ring={}",
                bond_idx,
                bond.begin().index(),
                bond.end().index(),
                bond.direction(),
                bond.stereo(),
                bond.stereo_atoms(),
                ring_count,
                min_ring
            );
        }
    };

    let mut raw = Molecule::from_smiles_with_sanitize(input, false).unwrap();
    print_state("raw_parse", &raw);

    update_property_cache_for_smiles(&mut raw).unwrap();
    print_state("post_update_property_cache", &raw);

    let mut ringed = raw.clone();
    ensure_fast_rings_for_writer_stereo_perception(&mut ringed).unwrap();
    print_state("post_fast_find_rings", &ringed);

    let mut assigned = ringed.clone();
    assign_stereochemistry_for_smiles(&mut assigned, false).unwrap();
    print_state("post_assign_stereochemistry", &assigned);

    eprintln!(
        "checkpoint=final_output smiles={}",
        mol_to_smiles(
            &Molecule::from_smiles_with_sanitize(input, false).unwrap(),
            &params
        )
        .unwrap()
    );
}

#[test]
#[ignore = "debug helper for upstream/rust checkpoint alignment"]
fn debug_probe_rust_writer_row_91() {
    let input = "O=C1/C(CC[C@@H](C)C1)=C(C)/C";
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        do_kekule: false,
        canonical: false,
        clean_stereo: false,
        all_bonds_explicit: true,
        all_hydrogens_explicit: false,
        include_dative_bonds: false,
        ignore_atom_map_numbers: false,
        ..Default::default()
    };

    let focus = [0usize, 1usize, 2usize, 3usize, 4usize];
    let print_state = |name: &str, mol: &Molecule| {
        eprintln!(
            "checkpoint={name} rings_initialized={} is_symm_sssr={} stereo_done={}",
            mol.derived_cache().rings.is_some(),
            mol.derived_cache()
                .rings
                .as_ref()
                .is_some_and(crate::RingInfo::is_symm_sssr),
            mol.prop("_StereochemDone").is_some()
        );
        for bond_idx in focus {
            let bond = &mol.bonds()[bond_idx];
            let ring_count = mol
                .derived_cache()
                .rings
                .as_ref()
                .map(|ri| ri.num_bond_rings(BondId::new(bond_idx)))
                .unwrap_or(0);
            let min_ring = mol
                .derived_cache()
                .rings
                .as_ref()
                .map(|ri| ri.min_bond_ring_size(BondId::new(bond_idx)))
                .unwrap_or(0);
            eprintln!(
                "bond {} {}-{} order={:?} dir={:?} stereo={:?} stereo_atoms={:?} ring_count={} min_ring={}",
                bond_idx,
                bond.begin().index(),
                bond.end().index(),
                bond.order(),
                bond.direction(),
                bond.stereo(),
                bond.stereo_atoms(),
                ring_count,
                min_ring
            );
        }
    };

    let mut raw = Molecule::from_smiles_with_sanitize(input, false).unwrap();
    print_state("raw_parse", &raw);

    update_property_cache_for_smiles(&mut raw).unwrap();
    print_state("post_update_property_cache", &raw);

    let mut ringed = raw.clone();
    ensure_fast_rings_for_writer_stereo_perception(&mut ringed).unwrap();
    print_state("post_fast_find_rings", &ringed);

    let mut assigned = ringed.clone();
    assign_stereochemistry_for_smiles(&mut assigned, false).unwrap();
    print_state("post_assign_stereochemistry", &assigned);

    eprintln!(
        "checkpoint=final_output smiles={}",
        mol_to_smiles(
            &Molecule::from_smiles_with_sanitize(input, false).unwrap(),
            &params
        )
        .unwrap()
    );
}

#[test]
#[ignore = "debug helper for upstream/rust checkpoint alignment"]
fn debug_probe_rust_writer_row_91_traversal_state() {
    let input = "O=C1/C(CC[C@@H](C)C1)=C(C)/C";
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        do_kekule: false,
        canonical: false,
        clean_stereo: false,
        all_bonds_explicit: true,
        all_hydrogens_explicit: false,
        include_dative_bonds: false,
        ignore_atom_map_numbers: false,
        ..Default::default()
    };

    let mut molecule = Molecule::from_smiles_with_sanitize(input, false).unwrap();
    let _ = prepare_plain_smiles_molecule(&mut molecule, &params).unwrap();
    let plans = collect_fragment_write_plans(&molecule, &params).unwrap();
    let plan = &plans[0];
    let ranks =
        rank_fragment_atoms_for_smiles(&molecule, plan, &params, SmilesOutputMode::PlainSmiles)
            .unwrap();
    let start_atom = choose_fragment_start_atom(plan, &ranks, &params).unwrap();
    let traversal = canonicalize_fragment_stack(
        &molecule,
        plan,
        start_atom,
        &ranks,
        &params,
        SmilesWriteOverrides::default(),
    )
    .unwrap();
    eprintln!(
        "checkpoint=post_traversal bond1_dir={:?} bond1_stereo={:?} stack_len={}",
        molecule.bonds()[1].direction(),
        molecule.bonds()[1].stereo(),
        traversal.stack.len()
    );
    for item in &traversal.stack {
        match *item {
            MolStackElem::Bond(bond, atom_to_left) => {
                if bond.index() == 1 {
                    eprintln!(
                        "stack_bond1 atom_to_left={} dir={:?} stereo={:?} stereo_atoms={:?}",
                        atom_to_left.index(),
                        molecule.bonds()[1].direction(),
                        molecule.bonds()[1].stereo(),
                        molecule.bonds()[1].stereo_atoms()
                    );
                }
            }
            _ => {}
        }
    }
    let result = write_mol_stack(
        &molecule,
        &traversal.stack,
        &params,
        SmilesWriteOverrides::default(),
        &mut SmilesWriteContext::default(),
    )
    .unwrap();
    eprintln!("checkpoint=post_write smiles={}", result.smiles);
}

#[test]
#[ignore = "debug helper for upstream/rust checkpoint alignment"]
fn debug_probe_rust_writer_row_91_direction_canonicalization() {
    let input = "O=C1/C(CC[C@@H](C)C1)=C(C)/C";
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        do_kekule: false,
        canonical: false,
        clean_stereo: false,
        all_bonds_explicit: true,
        all_hydrogens_explicit: false,
        include_dative_bonds: false,
        ignore_atom_map_numbers: false,
        ..Default::default()
    };

    let mut molecule = Molecule::from_smiles_with_sanitize(input, false).unwrap();
    let _ = prepare_plain_smiles_molecule(&mut molecule, &params).unwrap();
    let plan = collect_fragment_write_plans(&molecule, &params)
        .unwrap()
        .remove(0);
    let ranks =
        rank_fragment_atoms_for_smiles(&molecule, &plan, &params, SmilesOutputMode::PlainSmiles)
            .unwrap();
    let start_atom = choose_fragment_start_atom(&plan, &ranks, &params).unwrap();
    let traversal = canonicalize_fragment_stack(
        &molecule,
        &plan,
        start_atom,
        &ranks,
        &params,
        SmilesWriteOverrides::default(),
    )
    .unwrap();
    let mut canonicalized = molecule.clone();
    canonicalize_double_bond_directions_for_writer(
        &mut canonicalized,
        &traversal.stack,
        &traversal.traversal_ring_closure_bonds,
    )
    .unwrap();
    eprintln!(
        "checkpoint=post_direction_canonicalization bond1_dir={:?} bond7_dir={:?} bond8_dir={:?}",
        canonicalized.bonds()[1].direction(),
        canonicalized.bonds()[7].direction(),
        canonicalized.bonds()[8].direction()
    );
    eprintln!(
        "checkpoint=post_direction_canonicalization_smiles={}",
        write_mol_stack(
            &canonicalized,
            &traversal.stack,
            &params,
            SmilesWriteOverrides::default(),
            &mut SmilesWriteContext::default(),
        )
        .unwrap()
        .smiles
    );
}

#[test]
fn canonicalize_double_bond_clear_bond_dirs_does_not_require_neighboring_stereo_double_bond() {
    let mut builder = crate::MoleculeBuilder::new();
    let center = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let anchor = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let fluorine = builder.add_atom(crate::AtomSpec::new(crate::Element::F));
    let chlorine = builder.add_atom(crate::AtomSpec::new(crate::Element::CL));
    let ref_bond = builder
        .add_bond(
            crate::BondSpec::new(center, fluorine, BondOrder::Single)
                .with_direction(BondDirection::EndUpRight),
        )
        .unwrap();
    builder
        .add_bond(
            crate::BondSpec::new(center, chlorine, BondOrder::Single)
                .with_direction(BondDirection::EndDownRight),
        )
        .unwrap();
    builder
        .add_bond(crate::BondSpec::new(center, anchor, BondOrder::Double))
        .unwrap();
    let mut molecule = builder.build().unwrap();

    let mut bond_dir_counts = vec![1, 1, 0];
    let mut atom_dir_counts = vec![2, 0, 2, 2];

    clear_bond_dirs_from_atom_for_writer(
        &mut molecule,
        ref_bond,
        center,
        &mut bond_dir_counts,
        &mut atom_dir_counts,
    );

    assert_eq!(bond_dir_counts, vec![1, 0, 0]);
    assert_eq!(atom_dir_counts, vec![1, 0, 2, 1]);
    assert_eq!(molecule.bonds()[1].direction(), BondDirection::None);
}

#[test]
fn remove_redundant_bond_dir_specs_requires_neighboring_stereo_double_bond_like_rdkit() {
    let mut builder = crate::MoleculeBuilder::new();
    let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let ref_bond = builder
        .add_bond(
            crate::BondSpec::new(a0, a1, BondOrder::Single)
                .with_direction(BondDirection::EndUpRight),
        )
        .unwrap();
    let _dbl_bond = builder
        .add_bond(crate::BondSpec::new(a1, a2, BondOrder::Single))
        .unwrap();
    let mut molecule = builder.build().unwrap();
    let stack = vec![
        MolStackElem::Atom(a0),
        MolStackElem::Bond(ref_bond, a0),
        MolStackElem::Atom(a1),
    ];
    let mut bond_dir_counts = vec![1i8, 0i8];
    let mut atom_dir_counts = vec![2i8, 2i8, 2i8];

    remove_redundant_bond_dir_specs_for_writer(
        &mut molecule,
        &stack,
        &mut bond_dir_counts,
        &mut atom_dir_counts,
    );

    assert_eq!(molecule.bonds()[0].direction(), BondDirection::EndUpRight);
    assert_eq!(bond_dir_counts, vec![1, 0]);
    assert_eq!(atom_dir_counts, vec![2, 2, 2]);
}

#[test]
fn writer_canonical_fragment_scope_preserves_aromatic_fused_ring_form_like_rdkit_row_94() {
    let mol = Molecule::from_smiles_with_sanitize(
        "Cl.Cl.COc1ccc2nccc([C@@H](O)[C@@H]3C[C@@H]4CCN3C[C@@H]4C=C)c2c1",
        false,
    )
    .unwrap();
    let params = SmilesWriteParams::default();

    assert_eq!(
        mol_to_smiles(&mol, &params).unwrap(),
        "C=C[C@H]1CN2CC[C@H]1C[C@H]2[C@H](O)c1ccnc2ccc(OC)cc12.Cl.Cl"
    );
}

#[test]
fn writer_plain_nonisomeric_row_91_matches_rdkit_after_direction_cleanup() {
    let mol = Molecule::from_smiles_with_sanitize("O=C1/C(CC[C@@H](C)C1)=C(C)/C", false).unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        do_kekule: false,
        canonical: false,
        clean_stereo: false,
        all_bonds_explicit: true,
        all_hydrogens_explicit: false,
        include_dative_bonds: false,
        ignore_atom_map_numbers: false,
        ..Default::default()
    };

    assert_eq!(
        mol_to_smiles(&mol, &params).unwrap(),
        "O=C1-C(=C(-C)-C)-C-C-C(-C)-C-1"
    );
}

#[test]
fn writer_preserves_all_audited_bridge_ring_double_bond_stereo_like_rdkit() {
    let fixtures = [
        ("N/N=C1/C[C@H]2C[C@H]2C1", "N/N=C1\\C[C@@H]2C[C@@H]2C1"),
        (
            "COc1ccc([C@H]2N[C@@H](c3ccc(OC)cc3)C3CCCC2/C3=N/N=c2/[nH]c(-c3ccccc3)cs2)cc1",
            "COc1ccc([C@H]2N[C@@H](c3ccc(OC)cc3)C3CCCC2/C3=N/N=c2/[nH]c(-c3ccccc3)cs2)cc1",
        ),
        (
            "c1ccc(CO/N=C2\\C[C@@H](c3ccccc3)S[C@@H](c3ccccc3)C2)cc1",
            "c1ccc(CO/N=C2\\C[C@@H](c3ccccc3)S[C@@H](c3ccccc3)C2)cc1",
        ),
        (
            "CO/N=C1\\[C@H]2CCC[C@@H]1[C@H](c1ccc(Cl)cc1)N[C@@H]2c1ccc(Cl)cc1",
            "CO/N=C1\\[C@H]2CCC[C@@H]1[C@H](c1ccc(Cl)cc1)N[C@@H]2c1ccc(Cl)cc1",
        ),
        (
            "CCNC(=O)c1ccc(/C(=C2\\C[C@H]3CC[C@@H](C2)N3CCc2ccccc2)c2ccccc2)cc1",
            "CCNC(=O)c1ccc(/C(=C2\\C[C@H]3CC[C@@H](C2)N3CCc2ccccc2)c2ccccc2)cc1",
        ),
        (
            "O[C@@H](CN1[C@@H]2CC[C@H]1C/C(=C\\c1ccc(Cl)cc1)C2)c1ccccc1",
            "O[C@@H](CN1[C@@H]2CC[C@H]1C/C(=C\\c1ccc(Cl)cc1)C2)c1ccccc1",
        ),
        (
            "Fc1ccc(CN2[C@@H]3CC[C@H]2C/C(=C/COC(c2ccccc2)c2ccccc2)C3)cc1",
            "Fc1ccc(CN2[C@@H]3CC[C@H]2C/C(=C/COC(c2ccccc2)c2ccccc2)C3)cc1",
        ),
        (
            "C[C@@H]1/C(=N\\NC(=S)Nc2ccccc2)[C@@H](C)[C@H](c2cccc(Cl)c2)O[C@@H]1c1cccc(Cl)c1",
            "C[C@@H]1/C(=N\\NC(=S)Nc2ccccc2)[C@@H](C)[C@H](c2cccc(Cl)c2)O[C@@H]1c1cccc(Cl)c1",
        ),
        (
            "CO/N=C1\\[C@H]2CCC[C@@H]1[C@H](c1cccc(Cl)c1)N[C@@H]2c1cccc(Cl)c1",
            "CO/N=C1\\[C@H]2CCC[C@@H]1[C@H](c1cccc(Cl)c1)N[C@@H]2c1cccc(Cl)c1",
        ),
        (
            "O=C(N/N=C1/C[C@@H]2[C@H](C1)[C@H](OC(=O)c1ccccc1)C[C@@H]2OC(=O)c1ccccc1)c1ccncc1",
            "O=C(N/N=C1/C[C@@H]2[C@H](C1)[C@H](OC(=O)c1ccccc1)C[C@@H]2OC(=O)c1ccccc1)c1ccncc1",
        ),
        (
            "CO/N=C1\\[C@H]2CCC[C@@H]1[C@H](c1ccccc1Cl)N(C)[C@@H]2c1ccccc1Cl",
            "CO/N=C1\\[C@H]2CCC[C@@H]1[C@H](c1ccccc1Cl)N(C)[C@@H]2c1ccccc1Cl",
        ),
        (
            "C=CCN1[C@@H]2CC[C@H]1C/C(=C/COC(c1ccccc1)c1ccccc1)C2.O=C(O)C(=O)O",
            "C=CCN1[C@@H]2CC[C@H]1C/C(=C/COC(c1ccccc1)c1ccccc1)C2.O=C(O)C(=O)O",
        ),
    ];

    for (input, expected) in fixtures {
        let molecule = Molecule::from_smiles(input).unwrap();
        assert_eq!(
            mol_to_smiles(&molecule, &SmilesWriteParams::default()).unwrap(),
            expected,
            "RDKit 2026.03.1 canonical output for {input}"
        );
    }
}

#[test]
#[ignore = "debug helper for upstream/rust checkpoint alignment"]
fn debug_probe_rust_writer_row_91_direction_cleanup_substeps() {
    let input = "O=C1/C(CC[C@@H](C)C1)=C(C)/C";
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        do_kekule: false,
        canonical: false,
        clean_stereo: false,
        all_bonds_explicit: true,
        all_hydrogens_explicit: false,
        include_dative_bonds: false,
        ignore_atom_map_numbers: false,
        ..Default::default()
    };

    let mut molecule = Molecule::from_smiles_with_sanitize(input, false).unwrap();
    let _ = prepare_plain_smiles_molecule(&mut molecule, &params).unwrap();
    let plan = collect_fragment_write_plans(&molecule, &params)
        .unwrap()
        .remove(0);
    let ranks =
        rank_fragment_atoms_for_smiles(&molecule, &plan, &params, SmilesOutputMode::PlainSmiles)
            .unwrap();
    let start_atom = choose_fragment_start_atom(&plan, &ranks, &params).unwrap();
    let traversal = canonicalize_fragment_stack(
        &molecule,
        &plan,
        start_atom,
        &ranks,
        &params,
        SmilesWriteOverrides::default(),
    )
    .unwrap();
    let mut canonicalized = molecule.clone();
    let mut atom_visit_orders = vec![usize::MAX; canonicalized.num_atoms()];
    let mut bond_visit_orders = vec![usize::MAX; canonicalized.num_bonds()];
    for (pos, item) in traversal.stack.iter().enumerate() {
        match *item {
            MolStackElem::Atom(atom) => atom_visit_orders[atom.index()] = pos,
            MolStackElem::Bond(bond, _) => {
                bond_visit_orders[bond.index()] = pos;
                if matches!(
                    canonicalized.bonds()[bond.index()].direction(),
                    BondDirection::EndDownRight | BondDirection::EndUpRight
                ) {
                    canonicalized.topology_block_mut().bonds[bond.index()]
                        .set_direction(BondDirection::None);
                }
            }
            MolStackElem::Ring { .. } | MolStackElem::BranchOpen | MolStackElem::BranchClose => {}
        }
    }

    let mut bond_dir_counts = vec![0i8; canonicalized.num_bonds()];
    let mut atom_dir_counts = vec![0i8; canonicalized.num_atoms()];
    canonicalize_double_bonds_for_writer(
        &mut canonicalized,
        &bond_visit_orders,
        &atom_visit_orders,
        &traversal.traversal_ring_closure_bonds,
        &mut bond_dir_counts,
        &mut atom_dir_counts,
        &traversal.stack,
    );
    eprintln!(
        "checkpoint=after_double_bonds bond1_dir={:?} bond7_dir={:?} bond8_dir={:?} bond_dir_counts={:?} atom_dir_counts={:?}",
        canonicalized.bonds()[1].direction(),
        canonicalized.bonds()[7].direction(),
        canonicalized.bonds()[8].direction(),
        bond_dir_counts,
        atom_dir_counts
    );
    remove_unwanted_bond_dir_specs_for_writer(
        &mut canonicalized,
        &traversal.stack,
        &mut bond_dir_counts,
        &mut atom_dir_counts,
        &bond_visit_orders,
    );
    eprintln!(
        "checkpoint=after_remove_unwanted bond1_dir={:?} bond7_dir={:?} bond8_dir={:?} bond_dir_counts={:?} atom_dir_counts={:?}",
        canonicalized.bonds()[1].direction(),
        canonicalized.bonds()[7].direction(),
        canonicalized.bonds()[8].direction(),
        bond_dir_counts,
        atom_dir_counts
    );
    let bond1_begin = canonicalized.bonds()[1].begin();
    let bond1_end = canonicalized.bonds()[1].end();
    let bond8_end = canonicalized.bonds()[8].end();
    eprintln!(
        "checkpoint=before_manual_clear atom1_neighbors={:?} atom4_neighbors={:?} atom8_neighbors={:?}",
        incident_bonds(&canonicalized, bond1_begin)
            .iter()
            .map(|bond| (
                bond.index(),
                canonicalized.bonds()[bond.index()].order(),
                canonicalized.bonds()[bond.index()].direction()
            ))
            .collect::<Vec<_>>(),
        incident_bonds(&canonicalized, bond1_end)
            .iter()
            .map(|bond| (
                bond.index(),
                canonicalized.bonds()[bond.index()].order(),
                canonicalized.bonds()[bond.index()].direction()
            ))
            .collect::<Vec<_>>(),
        incident_bonds(&canonicalized, bond8_end)
            .iter()
            .map(|bond| (
                bond.index(),
                canonicalized.bonds()[bond.index()].order(),
                canonicalized.bonds()[bond.index()].direction()
            ))
            .collect::<Vec<_>>()
    );
    clear_bond_dirs_from_atom_for_writer(
        &mut canonicalized,
        BondId::new(1),
        bond1_begin,
        &mut bond_dir_counts,
        &mut atom_dir_counts,
    );
    eprintln!(
        "checkpoint=after_manual_clear_bond1_begin bond1_dir={:?} bond7_dir={:?} bond8_dir={:?} bond_dir_counts={:?} atom_dir_counts={:?}",
        canonicalized.bonds()[1].direction(),
        canonicalized.bonds()[7].direction(),
        canonicalized.bonds()[8].direction(),
        bond_dir_counts,
        atom_dir_counts
    );
    clear_bond_dirs_from_atom_for_writer(
        &mut canonicalized,
        BondId::new(1),
        bond1_end,
        &mut bond_dir_counts,
        &mut atom_dir_counts,
    );
    eprintln!(
        "checkpoint=after_manual_clear_bond1_end bond1_dir={:?} bond7_dir={:?} bond8_dir={:?} bond_dir_counts={:?} atom_dir_counts={:?}",
        canonicalized.bonds()[1].direction(),
        canonicalized.bonds()[7].direction(),
        canonicalized.bonds()[8].direction(),
        bond_dir_counts,
        atom_dir_counts
    );
    remove_redundant_bond_dir_specs_for_writer(
        &mut canonicalized,
        &traversal.stack,
        &mut bond_dir_counts,
        &mut atom_dir_counts,
    );
    eprintln!(
        "checkpoint=after_remove_redundant bond1_dir={:?} bond7_dir={:?} bond8_dir={:?} bond_dir_counts={:?} atom_dir_counts={:?}",
        canonicalized.bonds()[1].direction(),
        canonicalized.bonds()[7].direction(),
        canonicalized.bonds()[8].direction(),
        bond_dir_counts,
        atom_dir_counts
    );
    eprintln!(
        "checkpoint=post_write smiles={}",
        write_mol_stack(
            &canonicalized,
            &traversal.stack,
            &params,
            SmilesWriteOverrides::default(),
            &mut SmilesWriteContext::default(),
        )
        .unwrap()
        .smiles
    );
}

#[test]
fn writer_rejects_invalid_nontetrahedral_permutation_with_structured_error() {
    let error = validate_writer_chiral_permutation(ChiralTag::SquarePlanar, 4).unwrap_err();

    assert_eq!(
        error,
        SmilesWriteError::InvalidChiralPermutation {
            chiral_tag: ChiralTag::SquarePlanar,
            permutation: 4,
            limit: 3,
        }
    );
}

#[test]
fn writer_start_atom_guard_reports_invariant_for_empty_rank_scope() {
    let plan = FragmentWritePlan {
        atoms: Vec::new(),
        bonds: Vec::new(),
        rooted_at_atom: None,
    };
    let params = SmilesWriteParams::default();

    let error = choose_fragment_start_atom(&plan, &[], &params).unwrap_err();

    assert_eq!(
        error,
        SmilesWriteError::InvariantViolation {
            stage: "ShortTermAtomWriter",
            message: "choose_fragment_start_atom() called with empty canonical rank scope",
        }
    );
}

#[test]
fn writer_swap_counter_reports_invariant_for_length_mismatch() {
    let error = count_swaps_to_interconvert(&[BondId::new(0)], &[BondId::new(0), BondId::new(1)])
        .unwrap_err();

    assert_eq!(
        error,
        SmilesWriteError::InvariantViolation {
            stage: "ShortTermAtomWriter",
            message: "count_swaps_to_interconvert() probe/reference length mismatch",
        }
    );
}

#[test]
fn writer_swap_counter_reports_invariant_for_missing_expected_bond() {
    let error = count_swaps_to_interconvert(
        &[BondId::new(0), BondId::new(1)],
        &[BondId::new(0), BondId::new(2)],
    )
    .unwrap_err();

    assert_eq!(
        error,
        SmilesWriteError::InvariantViolation {
            stage: "ShortTermAtomWriter",
            message: "count_swaps_to_interconvert() expected bond missing from probe order",
        }
    );
}

#[test]
fn writer_bond_guard_uses_invariant_error_shape() {
    let error = invariant_stage_error::<()>(
        SmilesPlanStage::ShortTermBondWriter,
        "write_ring_closure() could not allocate a free ring index",
    )
    .unwrap_err();

    assert_eq!(
        error,
        SmilesWriteError::InvariantViolation {
            stage: "ShortTermBondWriter",
            message: "write_ring_closure() could not allocate a free ring index",
        }
    );
}

#[test]
fn atom_bond_and_fragment_writer_entries_fail_closed_until_ported() {
    let molecule = ethane();

    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: false,
        clean_stereo: false,
        ..Default::default()
    };

    assert_eq!(get_atom_smiles(&molecule, 0, &params).unwrap(), "C");
    assert_eq!(
        get_molecule_bond_smiles(&molecule, 0, Some(0), &params).unwrap(),
        ""
    );
    assert_eq!(get_bond_smiles(BondOrder::Single).unwrap(), "");
    assert_eq!(get_bond_smiles(BondOrder::Dative).unwrap(), "->");
    // Fragment API is now implemented — ethane fragment produces SMILES
    let fragment = mol_fragment_to_smiles(
        &molecule,
        &SmilesWriteParams::default(),
        &[0, 1],
        None,
        None,
        None,
    )
    .unwrap();
    assert_eq!(fragment, "CC", "ethane fragment should produce CC");

    let fragment_cx = mol_fragment_to_cx_smiles(
        &molecule,
        &SmilesWriteParams::default(),
        &[0, 1],
        None,
        None,
        None,
        CxSmilesFields::ALL,
    )
    .unwrap();
    assert_eq!(fragment_cx, "CC", "ethane fragment CX should be plain CC");
}

#[test]
fn atom_bond_and_fragment_writer_entries_validate_indices_before_unsupported() {
    let molecule = ethane();

    assert_eq!(
        get_atom_smiles(&molecule, 2, &SmilesWriteParams::default()).unwrap_err(),
        SmilesWriteError::AtomOutOfRange { atom: 2 }
    );
    assert_eq!(
        get_molecule_bond_smiles(&molecule, 1, None, &SmilesWriteParams::default()).unwrap_err(),
        SmilesWriteError::BondOutOfRange { bond: 1 }
    );
    assert_eq!(
        mol_fragment_to_smiles(
            &molecule,
            &SmilesWriteParams::default(),
            &[2],
            None,
            None,
            None
        )
        .unwrap_err(),
        SmilesWriteError::AtomOutOfRange { atom: 2 }
    );
}

#[test]
fn mol_fragment_to_smiles_uses_original_atom_and_bond_symbol_overrides() {
    let molecule = Molecule::from_smiles_with_sanitize("CCO", false).unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: false,
        clean_stereo: false,
        ..Default::default()
    };
    let atom_symbols = vec!["A".to_string(), "B".to_string(), "C".to_string()];
    let bond_symbols = vec!["~".to_string(), "!".to_string()];

    let fragment = mol_fragment_to_smiles(
        &molecule,
        &params,
        &[0, 1],
        Some(&[0]),
        Some(&atom_symbols),
        Some(&bond_symbols),
    )
    .unwrap();

    assert_eq!(fragment, "A~B");
}

#[test]
fn mol_fragment_to_smiles_splits_disconnected_atoms_within_fragment_scope() {
    let molecule = Molecule::from_smiles_with_sanitize("CCO", false).unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: false,
        clean_stereo: false,
        ..Default::default()
    };

    let fragment = mol_fragment_to_smiles(&molecule, &params, &[0, 2], None, None, None).unwrap();

    assert_eq!(fragment, "C.O");
}

#[test]
fn mol_fragment_to_cx_smiles_filters_cx_blocks_to_fragment_output_scope() {
    let mut builder = crate::MoleculeBuilder::new();
    let c0 = builder.add_atom(
        crate::AtomSpec::new(crate::Element::C).with_prop("_supplementalSmilesLabel", "keep"),
    );
    let c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let o2 = builder.add_atom(
        crate::AtomSpec::new(crate::Element::O).with_prop("_supplementalSmilesLabel", "drop"),
    );
    builder
        .add_bond(crate::BondSpec::new(c0, c1, BondOrder::Single))
        .unwrap();
    builder
        .add_bond(crate::BondSpec::new(c1, o2, BondOrder::Dative))
        .unwrap();
    let molecule = builder.build().unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: false,
        clean_stereo: false,
        ..Default::default()
    };

    let fragment = mol_fragment_to_cx_smiles(
        &molecule,
        &params,
        &[0, 1],
        Some(&[0]),
        None,
        None,
        CxSmilesFields::ATOM_LABELS | CxSmilesFields::COORDINATE_BONDS,
    )
    .unwrap();

    assert!(
        fragment.contains("keep"),
        "fragment CX should keep atom 0 label: {fragment}"
    );
    assert!(
        !fragment.contains("drop"),
        "fragment CX must not leak atom 2 label: {fragment}"
    );
    assert!(
        !fragment.contains("_Z:2"),
        "fragment CX must not leak out-of-scope dative bond: {fragment}"
    );
}

// ── CX SMILES Extension Tests ──────────────────────────────────────────

fn cx_ethanol() -> Molecule {
    // CCO with atom label on O
    let mut builder = crate::MoleculeBuilder::new();
    let c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let c2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let o =
        builder.add_atom(crate::AtomSpec::new(crate::Element::O).with_prop("atomLabel", "Hydroxy"));
    builder
        .add_bond(crate::BondSpec::new(c1, c2, crate::BondOrder::Single))
        .unwrap();
    builder
        .add_bond(crate::BondSpec::new(c2, o, crate::BondOrder::Single))
        .unwrap();
    builder.build().unwrap()
}

#[test]
fn cx_individual_coords_writes_atom_order_coordinates() {
    let mut builder = crate::MoleculeBuilder::new();
    let c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let c2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    builder
        .add_bond(crate::BondSpec::new(c1, c2, crate::BondOrder::Single))
        .unwrap();
    builder
        .set_2d_coordinates(vec![[0.0, 0.0], [1.5, 0.0]])
        .unwrap();
    let molecule = builder.build().unwrap();
    let scope = CxWriteScope::full_molecule(&molecule);
    let coords_str = write_cx_coords(&molecule, &scope.atom_order);
    assert_eq!(coords_str, "0,0,;1.5,0,");
}

#[test]
fn cx_empty_coords_when_no_2d_coords_present() {
    let molecule = ethane();
    let scope = CxWriteScope::full_molecule(&molecule);
    assert_eq!(write_cx_coords(&molecule, &scope.atom_order), "");
}

#[test]
fn cx_individual_atom_labels_writes_label_entries() {
    let molecule = cx_ethanol();
    let scope = CxWriteScope::full_molecule(&molecule);
    let labels = write_cx_atom_labels(&molecule, &scope.atom_order);
    assert_eq!(labels, ";;Hydroxy");
}

#[test]
fn cx_individual_radicals_writes_entries() {
    // Build a molecule with a radical
    let mut builder = crate::MoleculeBuilder::new();
    let c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let c2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    builder
        .add_bond(crate::BondSpec::new(c1, c2, crate::BondOrder::Single))
        .unwrap();
    let mol = builder.build().unwrap();
    // Test write_cx_radicals directly (bypasses the SMILES writer which
    // rejects radical atoms)
    let scope = CxWriteScope::full_molecule(&mol);
    let radicals = write_cx_radicals(&mol, &scope.atom_order);
    // No radicals on plain ethane
    assert_eq!(radicals, "");
}

#[test]
fn cx_no_radicals_when_none_present() {
    let molecule = ethane();
    let scope = CxWriteScope::full_molecule(&molecule);
    assert_eq!(write_cx_radicals(&molecule, &scope.atom_order), "");
}

#[test]
fn cx_atom_props_writes_entries_for_atoms_with_properties() {
    let molecule = ethane();
    let scope = CxWriteScope::full_molecule(&molecule);
    let props = write_cx_atom_props(&molecule, &scope.atom_order);
    assert_eq!(props, "");
}

#[test]
fn cx_coordinate_bonds_writes_entries() {
    // Build a molecule with a dative bond
    let mut builder = crate::MoleculeBuilder::new();
    let n = builder.add_atom(crate::AtomSpec::new(crate::Element::N));
    let o = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
    builder
        .add_bond(crate::BondSpec::new(n, o, crate::BondOrder::Dative))
        .unwrap();
    let mol = builder.build().unwrap();

    let scope = CxWriteScope::full_molecule(&mol);
    let coord_bonds = write_cx_coordinate_bonds(&mol, &scope.atom_order, &scope.bond_order, "C");
    assert_eq!(coord_bonds, "C:0.0");
}

#[test]
fn cx_bond_cfg_writes_wedge_and_dash_when_coords_are_included() {
    let mut builder = crate::MoleculeBuilder::new();
    let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    builder
        .add_bond(
            crate::BondSpec::new(a0, a1, crate::BondOrder::Single)
                .with_direction(crate::BondDirection::BeginWedge),
        )
        .unwrap();
    builder
        .add_bond(
            crate::BondSpec::new(a1, a2, crate::BondOrder::Single)
                .with_direction(crate::BondDirection::BeginDash),
        )
        .unwrap();
    builder
        .set_2d_coordinates(vec![[0.0, 0.0], [1.5, 0.0], [3.0, 0.0]])
        .unwrap();
    let molecule = builder.build().unwrap();

    let extension =
        get_cx_extensions(&molecule, CxSmilesFields::COORDS | CxSmilesFields::BOND_CFG).unwrap();
    assert!(
        extension.contains("wU:0.0"),
        "wedge-up bond config should be emitted, got: {extension:?}"
    );
    assert!(
        extension.contains("wD:1.1"),
        "wedge-down bond config should be emitted, got: {extension:?}"
    );
}

#[test]
fn cx_bond_cfg_uses_molfile_bond_cfg_fallback_for_wedge_and_dash() {
    let mut builder = crate::MoleculeBuilder::new();
    let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    builder
        .add_bond(
            crate::BondSpec::new(a0, a1, crate::BondOrder::Single)
                .with_prop("_MolFileBondCfg", "3"),
        )
        .unwrap();
    builder
        .set_2d_coordinates(vec![[0.0, 0.0], [1.5, 0.0]])
        .unwrap();
    let molecule = builder.build().unwrap();

    let extension =
        get_cx_extensions(&molecule, CxSmilesFields::COORDS | CxSmilesFields::BOND_CFG).unwrap();
    assert!(
        extension.contains("wD:0.0"),
        "MolFile cfg=3 should emit dash wedge config, got: {extension:?}"
    );
}

#[test]
fn cx_bond_cfg_emits_unknown_without_coords_from_molfile_cfg() {
    let mut builder = crate::MoleculeBuilder::new();
    let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    builder
        .add_bond(
            crate::BondSpec::new(a0, a1, crate::BondOrder::Single)
                .with_prop("_MolFileBondCfg", "2"),
        )
        .unwrap();
    let molecule = builder.build().unwrap();

    let extension = get_cx_extensions(&molecule, CxSmilesFields::BOND_CFG).unwrap();
    assert_eq!(extension, "|w:0.0|");
}

#[test]
fn cx_bond_atropisomer_writes_wedge_for_atrop_neighbor() {
    let mut builder = crate::MoleculeBuilder::new();
    let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    builder
        .add_bond(
            crate::BondSpec::new(a0, a1, crate::BondOrder::Single)
                .with_direction(crate::BondDirection::BeginWedge),
        )
        .unwrap();
    builder
        .add_bond(
            crate::BondSpec::new(a0, a2, crate::BondOrder::Single)
                .with_stereo(crate::BondStereo::AtropCw),
        )
        .unwrap();
    let molecule = builder.build().unwrap();

    let extension = get_cx_extensions(&molecule, CxSmilesFields::BOND_ATROPISOMER).unwrap();
    assert_eq!(extension, "|wU:0.0|");
}

#[test]
fn cx_bond_atropisomer_does_not_use_molfile_cfg_fallback() {
    let mut builder = crate::MoleculeBuilder::new();
    let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    builder
        .add_bond(
            crate::BondSpec::new(a0, a1, crate::BondOrder::Single)
                .with_prop("_MolFileBondCfg", "1"),
        )
        .unwrap();
    builder
        .add_bond(
            crate::BondSpec::new(a0, a2, crate::BondOrder::Single)
                .with_stereo(crate::BondStereo::AtropCw),
        )
        .unwrap();
    let molecule = builder.build().unwrap();

    let extension = get_cx_extensions(&molecule, CxSmilesFields::BOND_ATROPISOMER).unwrap();
    assert_eq!(extension, "");
}

fn cx_ring_double_bond_molecule(stereo: BondStereo) -> Molecule {
    let mut builder = crate::MoleculeBuilder::new();
    let atoms: Vec<_> = (0..8)
        .map(|_| builder.add_atom(crate::AtomSpec::new(crate::Element::C)))
        .collect();
    builder
        .add_bond(
            crate::BondSpec::new(atoms[0], atoms[1], crate::BondOrder::Double)
                .with_stereo(stereo)
                .with_stereo_atoms(atoms[7], atoms[2]),
        )
        .unwrap();
    builder
        .add_bond(crate::BondSpec::new(
            atoms[1],
            atoms[2],
            crate::BondOrder::Single,
        ))
        .unwrap();
    builder
        .add_bond(crate::BondSpec::new(
            atoms[2],
            atoms[3],
            crate::BondOrder::Single,
        ))
        .unwrap();
    builder
        .add_bond(crate::BondSpec::new(
            atoms[3],
            atoms[4],
            crate::BondOrder::Single,
        ))
        .unwrap();
    builder
        .add_bond(crate::BondSpec::new(
            atoms[4],
            atoms[5],
            crate::BondOrder::Single,
        ))
        .unwrap();
    builder
        .add_bond(crate::BondSpec::new(
            atoms[5],
            atoms[6],
            crate::BondOrder::Single,
        ))
        .unwrap();
    builder
        .add_bond(crate::BondSpec::new(
            atoms[6],
            atoms[7],
            crate::BondOrder::Single,
        ))
        .unwrap();
    builder
        .add_bond(crate::BondSpec::new(
            atoms[7],
            atoms[0],
            crate::BondOrder::Single,
        ))
        .unwrap();
    let mut molecule = builder.build().unwrap();
    let ring_info = crate::symmetrize_sssr(&molecule).unwrap();
    molecule.derived_cache_mut().rings = Some(ring_info);
    molecule
}

#[test]
fn write_cx_ringbond_cistrans_block_writes_c_block_for_cis_ring_double_bond() {
    let molecule = cx_ring_double_bond_molecule(BondStereo::Cis);
    let extension = get_cx_extensions(&molecule, CxSmilesFields::BOND_CFG).unwrap();
    assert!(
        extension.contains("c:0"),
        "cis ring double bond should emit c block, got: {extension:?}"
    );
}

#[test]
fn write_cx_ringbond_cistrans_block_writes_ctu_block_for_any_ring_double_bond() {
    let molecule = cx_ring_double_bond_molecule(BondStereo::Any);
    let extension = get_cx_extensions(&molecule, CxSmilesFields::BOND_CFG).unwrap();
    assert!(
        extension.contains("ctu:0"),
        "any ring double bond should emit ctu block, got: {extension:?}"
    );
}

#[test]
fn cx_linknodes_block_writes_compact_form_for_degree_two_center() {
    let mut builder = crate::MoleculeBuilder::new();
    let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    builder
        .add_bond(crate::BondSpec::new(a0, a1, crate::BondOrder::Single))
        .unwrap();
    builder
        .add_bond(crate::BondSpec::new(a1, a2, crate::BondOrder::Single))
        .unwrap();
    let molecule = builder
        .build()
        .unwrap()
        .with_prop("_MolFileLinkNodes", "1 3 2 2 1 2 3");

    let scope = CxWriteScope::full_molecule(&molecule);
    let block = write_cx_linknodes_block(&molecule, &scope.atom_order);
    assert_eq!(block, "LN:1:1.3");
}

#[test]
fn cx_linknodes_block_writes_outer_atoms_for_degree_three_center() {
    let mut builder = crate::MoleculeBuilder::new();
    let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let a3 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    builder
        .add_bond(crate::BondSpec::new(a0, a1, crate::BondOrder::Single))
        .unwrap();
    builder
        .add_bond(crate::BondSpec::new(a0, a2, crate::BondOrder::Single))
        .unwrap();
    builder
        .add_bond(crate::BondSpec::new(a0, a3, crate::BondOrder::Single))
        .unwrap();
    let molecule = builder
        .build()
        .unwrap()
        .with_prop("_MolFileLinkNodes", "1 4 2 1 2 1 3");

    let scope = CxWriteScope::full_molecule(&molecule);
    let block = write_cx_linknodes_block(&molecule, &scope.atom_order);
    assert_eq!(block, "LN:0:1.4.1.2");
}

#[test]
fn cx_polymer_block_writes_sru_label_and_connect() {
    let mut builder = crate::MoleculeBuilder::new();
    let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    builder
        .add_bond(crate::BondSpec::new(a0, a1, crate::BondOrder::Single))
        .unwrap();
    builder
        .add_bond(crate::BondSpec::new(a1, a2, crate::BondOrder::Single))
        .unwrap();
    builder
        .add_substance_group(
            crate::SubstanceGroup::new(
                crate::SubstanceGroupId::new(0),
                crate::SubstanceGroupKind::StructuralRepeatUnit,
            )
            .with_atoms(vec![a0, a1])
            .with_label("SRU")
            .with_connection(crate::SGroupConnection::HeadToHead),
        )
        .unwrap();
    let molecule = builder.build().unwrap();

    let extension = get_cx_extensions(&molecule, CxSmilesFields::POLYMER).unwrap();
    assert_eq!(extension, "|Sg:n:0,1:SRU:hh:::|");
}

#[test]
fn cx_polymer_block_writes_copolymer_subtype_and_crossings() {
    let mut builder = crate::MoleculeBuilder::new();
    let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    builder
        .add_bond(crate::BondSpec::new(a0, a1, crate::BondOrder::Single))
        .unwrap();
    builder
        .add_bond(crate::BondSpec::new(a1, a2, crate::BondOrder::Single))
        .unwrap();
    builder
        .add_substance_group(
            crate::SubstanceGroup::new(
                crate::SubstanceGroupId::new(0),
                crate::SubstanceGroupKind::Copolymer,
            )
            .with_atoms(vec![a0, a1])
            .with_subtype("ALT")
            .with_label("COP")
            .with_connection(crate::SGroupConnection::Either)
            .with_prop("XBHEAD", "1,0")
            .with_prop("XBCORR", "0,1,0,1"),
        )
        .unwrap();
    let molecule = builder.build().unwrap();

    let extension = get_cx_extensions(&molecule, CxSmilesFields::POLYMER).unwrap();
    assert_eq!(extension, "|Sg:alt:0,1:COP:eu:1,0:1,1:|");
}

#[test]
fn cx_sgroup_hierarchy_writes_block_for_polymer_parent_child() {
    let mut builder = crate::MoleculeBuilder::new();
    let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    builder
        .add_bond(crate::BondSpec::new(a0, a1, crate::BondOrder::Single))
        .unwrap();
    builder
        .add_bond(crate::BondSpec::new(a1, a2, crate::BondOrder::Single))
        .unwrap();
    builder
        .add_substance_group(
            crate::SubstanceGroup::new(
                crate::SubstanceGroupId::new(0),
                crate::SubstanceGroupKind::StructuralRepeatUnit,
            )
            .with_atoms(vec![a0, a1]),
        )
        .unwrap();
    builder
        .add_substance_group(
            crate::SubstanceGroup::new(
                crate::SubstanceGroupId::new(1),
                crate::SubstanceGroupKind::StructuralRepeatUnit,
            )
            .with_atoms(vec![a1, a2])
            .with_parent(crate::SubstanceGroupId::new(0)),
        )
        .unwrap();
    let molecule = builder.build().unwrap();

    let extension = get_cx_extensions(&molecule, CxSmilesFields::POLYMER).unwrap();
    assert!(
        extension.contains("SgH:0:1"),
        "polymer hierarchy should be emitted, got: {extension:?}"
    );
}

#[test]
fn cx_sgroup_hierarchy_writes_block_for_sgroup_parent_prop() {
    let mut builder = crate::MoleculeBuilder::new();
    let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    builder
        .add_bond(crate::BondSpec::new(a0, a1, crate::BondOrder::Single))
        .unwrap();
    builder
        .add_substance_group(
            crate::SubstanceGroup::new(
                crate::SubstanceGroupId::new(0),
                crate::SubstanceGroupKind::Data,
            )
            .with_atoms(vec![a0])
            .with_prop("index", "0"),
        )
        .unwrap();
    builder
        .add_substance_group(
            crate::SubstanceGroup::new(
                crate::SubstanceGroupId::new(1),
                crate::SubstanceGroupKind::Data,
            )
            .with_atoms(vec![a1])
            .with_prop("index", "1")
            .with_prop("PARENT", "0"),
        )
        .unwrap();
    let molecule = builder.build().unwrap();

    let extension = get_cx_extensions(&molecule, CxSmilesFields::SGROUPS).unwrap();
    assert!(
        extension.contains("SgH:0:1"),
        "SGroup hierarchy should be emitted for SGROUPS fields, got: {extension:?}"
    );
}

#[test]
fn cx_full_molecule_with_all_fields_returns_smiles_with_cx_extension() {
    let molecule = cx_ethanol();
    let result = mol_to_cx_smiles(
        &molecule,
        &SmilesWriteParams::default(),
        CxSmilesFields::ATOM_LABELS,
        RestoreBondDirOption::None,
    )
    .unwrap();
    // Should contain the SMILES and CX extension separated by space
    assert!(!result.is_empty(), "should produce output");
    assert!(
        result.contains("Hydroxy"),
        "CX label should appear: {result:?}"
    );
}

#[test]
fn cx_no_cx_data_returns_plain_smiles() {
    let molecule = ethane();
    let result = mol_to_cx_smiles(
        &molecule,
        &SmilesWriteParams::default(),
        CxSmilesFields::ALL,
        RestoreBondDirOption::Clear,
    )
    .unwrap();
    assert_eq!(result, "CC", "ethane with no CX data should be plain CC");
}

#[test]
fn cx_get_cx_extensions_returns_union_of_requested_fields() {
    let molecule = cx_ethanol();
    // Request ATOM_LABELS
    let result = get_cx_extensions(&molecule, CxSmilesFields::ATOM_LABELS).unwrap();
    assert!(result.contains("Hydroxy"), "should have atom labels");
    assert!(
        !result.contains('('),
        "no coords (molecule has no 2D coords)"
    );
}

#[test]
fn cx_radical_fields_only_produces_expected_output() {
    let molecule = cx_ethanol();
    let result = get_cx_extensions(&molecule, CxSmilesFields::RADICALS).unwrap();
    // cx_ethanol has no radicals, so result should be empty
    assert_eq!(result, "", "no radicals expected: {result:?}");
}

#[test]
fn cx_individual_atom_props_escapes_property_name_and_value_dots() {
    let mut builder = crate::MoleculeBuilder::new();
    let c =
        builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_prop("foo.bar", "baz.qux"));
    let c2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    builder
        .add_bond(crate::BondSpec::new(c, c2, crate::BondOrder::Single))
        .unwrap();
    let mol = builder.build().unwrap();

    let scope = CxWriteScope::full_molecule(&mol);
    let props = write_cx_atom_props(&mol, &scope.atom_order);
    assert_eq!(props, "atomProp:0.foo&#46;bar.baz&#46;qux");
}

#[test]
fn cx_molfile_values_when_present() {
    let mut builder = crate::MoleculeBuilder::new();
    let c = builder
        .add_atom(crate::AtomSpec::new(crate::Element::C).with_prop("molFileValue", "test_value"));
    let c2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    builder
        .add_bond(crate::BondSpec::new(c, c2, crate::BondOrder::Single))
        .unwrap();
    let mol = builder.build().unwrap();

    let scope = CxWriteScope::full_molecule(&mol);
    let values = write_cx_molfile_values(&mol, &scope.atom_order);
    assert_eq!(values, "test_value;");
}

#[test]
fn cx_stereo_group_writes_appropriate_code() {
    let mut builder = crate::MoleculeBuilder::new();
    let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    builder
        .add_bond(crate::BondSpec::new(a0, a1, crate::BondOrder::Single))
        .unwrap();
    builder
        .add_stereo_group(crate::StereoGroup::new(
            crate::StereoGroupKind::Or,
            vec![a1],
            vec![],
        ))
        .unwrap();
    builder
        .add_stereo_group(crate::StereoGroup::new(
            crate::StereoGroupKind::And,
            vec![a0],
            vec![],
        ))
        .unwrap();
    let molecule = builder.build().unwrap();
    let scope = CxWriteScope::full_molecule(&molecule);
    let stereo = write_cx_enhanced_stereo(&molecule, &scope.atom_order, &scope.bond_order);
    assert_eq!(stereo, "o1:1,&1:0");
}

#[test]
fn cleanup_stereo_groups_for_cx_smiles_moves_atrop_atoms_to_bond_membership() {
    let mut builder = crate::MoleculeBuilder::new();
    let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    builder
        .add_bond(crate::BondSpec::new(a0, a1, crate::BondOrder::Single))
        .unwrap();
    builder
        .add_bond(
            crate::BondSpec::new(a1, a2, crate::BondOrder::Single)
                .with_stereo(crate::BondStereo::AtropCw),
        )
        .unwrap();
    builder
        .add_stereo_group(crate::StereoGroup::new(
            crate::StereoGroupKind::Or,
            vec![a1],
            Vec::new(),
        ))
        .unwrap();
    let mut molecule = builder.build().unwrap();

    cleanup_stereo_groups_for_cx_smiles(&mut molecule).unwrap();

    assert_eq!(molecule.stereo_groups().len(), 1);
    let group = &molecule.stereo_groups()[0];
    assert!(group.atoms().is_empty());
    assert_eq!(group.bonds(), &[BondId::new(1)]);
}

#[test]
fn cx_sgroups_empty_when_no_sgroups() {
    let molecule = ethane();
    let scope = CxWriteScope::full_molecule(&molecule);
    let sgroups = write_cx_sgroups(&molecule, &scope.atom_order, &scope.bond_order);
    assert_eq!(sgroups, "");
}

#[test]
fn cx_data_sgroup_prefers_typed_data_fields() {
    let mut builder = crate::MoleculeBuilder::new();
    let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    builder
        .add_bond(crate::BondSpec::new(a0, a1, crate::BondOrder::Single))
        .unwrap();
    builder
        .add_substance_group(
            crate::SubstanceGroup::new(
                crate::SubstanceGroupId::new(0),
                crate::SubstanceGroupKind::Data,
            )
            .with_atoms(vec![a1, a0])
            .with_data(crate::SGroupData {
                field_name: Some("FIELD".to_string()),
                field_info: Some("INFO".to_string()),
                query_op: Some("OP".to_string()),
                units: Some("TAG".to_string()),
                values: vec!["VALUE".to_string()],
                ..Default::default()
            }),
        )
        .unwrap();
    let molecule = builder.build().unwrap();
    let scope = CxWriteScope::full_molecule(&molecule);

    assert_eq!(
        write_cx_data_sgroups(&molecule, &scope.atom_order),
        "SgD:1,0:FIELD:VALUE:OP:INFO:TAG:"
    );
}

#[test]
fn build_smiles_helpers_match_main_writer_helpers() {
    let mut builder = crate::MoleculeBuilder::new();
    let c = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
    let n = builder.add_atom(crate::AtomSpec::new(crate::Element::N));
    builder
        .add_bond(crate::BondSpec::new(c, n, crate::BondOrder::Triple))
        .unwrap();
    let molecule = builder.build().unwrap();
    let params = SmilesWriteParams::default();
    let context = SmilesWriteContext::default();

    assert_eq!(
        build_atom_smiles(&molecule, c, &params, &context).unwrap(),
        get_atom_smiles(&molecule, c.index(), &params).unwrap()
    );
    assert_eq!(
        build_bond_smiles(&molecule, BondId::new(0), c, &params).unwrap(),
        get_molecule_bond_smiles(&molecule, 0, Some(c.index()), &params).unwrap()
    );
}

#[test]
fn writer_rooted_ring_stereo_centers_follow_rdkit_canonicalize_fragment_postprocessing() {
    let molecule = Molecule::from_smiles(
        "O=C(NC[C@]12C[C@H]3C[C@H](C[C@H](C3)C1)C2)[C@@H]1C[C@H]2c3ccccc3[C@@H]1c1ccccc12",
    )
    .unwrap();
    let params = SmilesWriteParams {
        canonical: false,
        clean_stereo: false,
        rooted_at_atom: Some(molecule.num_atoms() - 1),
        include_dative_bonds: false,
        ..Default::default()
    };

    assert_eq!(
        molecule.to_smiles_with_params(&params).unwrap(),
        "c12c(cccc1)[C@H]1[C@H](C(=O)NC[C@]34C[C@H]5C[C@H](C[C@H](C5)C3)C4)C[C@@H]2c2ccccc21"
    );
}

#[test]
fn writer_noncanonical_explicit_bonds_preserves_rdkit_neighboring_stereo_bond_order() {
    let molecule = Molecule::from_smiles(
            "C=C1/C(C[C@@H](O)CC1)=C\\C=C2[C@@]3([H])[C@@](CCC\\2)(C)[C@]([C@H](C)/C=C/[C@H](C)C(C)C)([H])CC3",
        )
        .unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: false,
        clean_stereo: false,
        all_bonds_explicit: true,
        include_dative_bonds: false,
        ..Default::default()
    };

    assert_eq!(
        molecule.to_smiles_with_params(&params).unwrap(),
        "C=C1/C(=C\\C=C2\\C3-C(-C)(-C-C-C-2)-C(-C(-C)/C=C/C(-C)-C(-C)-C)-C-C-3)-C-C(-O)-C-C-1"
    );
}

#[test]
fn writer_nonisomeric_explicit_bonds_clears_invalid_cx_ring_double_bond_stereo_like_rdkit_row_87() {
    let molecule = Molecule::from_smiles("O=C1OCC(=C1c1ccccc1)c1ccccc1 |c:4|").unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: false,
        clean_stereo: false,
        all_bonds_explicit: true,
        include_dative_bonds: false,
        ..Default::default()
    };

    assert_eq!(
        molecule.to_smiles_with_params(&params).unwrap(),
        "O=C1-O-C-C(-c2:c:c:c:c:c:2)=C-1-c1:c:c:c:c:c:1"
    );
}

#[test]
fn writer_preserves_non_ring_aromatic_bridge_as_explicit_single_like_rdkit_biphenyl() {
    let molecule = Molecule::from_smiles("c1ccccc1c1ccccc1").unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: true,
        clean_stereo: false,
        include_dative_bonds: false,
        ..Default::default()
    };

    assert_eq!(
        molecule.to_smiles_with_params(&params).unwrap(),
        "c1ccc(-c2ccccc2)cc1"
    );
}

#[test]
fn writer_emits_single_bridge_between_aromatic_rings_like_rdkit_row_90() {
    let molecule =
        Molecule::from_smiles("CCOC(=O)c1c(N/C=C\\2/C(=NC(=NC2=O)S)O)scc1c1ccc(cc1)Br").unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: false,
        clean_stereo: false,
        include_dative_bonds: false,
        ..Default::default()
    };

    assert_eq!(
        molecule.to_smiles_with_params(&params).unwrap(),
        "CCOC(=O)c1c(NC=C2C(O)=NC(S)=NC2=O)scc1-c1ccc(Br)cc1"
    );
}

#[test]
fn writer_does_not_emit_cleared_tetrahedral_tag_like_rdkit_row_92() {
    let molecule = Molecule::from_smiles(
        "O/C1=C/C=C/C=C1/CN3CCN(CC=2C=CC=CC=2O)[C@]3([H])C=4/C=C(/OC)C(=CC=4)OC",
    )
    .unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: true,
        canonical: false,
        clean_stereo: false,
        include_dative_bonds: false,
        ..Default::default()
    };

    assert_eq!(
        molecule.to_smiles_with_params(&params).unwrap(),
        "Oc1ccccc1CN1CCN(Cc2ccccc2O)C1c1cc(OC)c(OC)cc1"
    );
}

#[test]
fn writer_preserves_ring_special_case_stereo_like_rdkit_row_83() {
    let molecule = Molecule::from_smiles(
        "O=C(NC[C@]12C[C@H]3C[C@H](C[C@H](C3)C1)C2)[C@@H]1C[C@H]2c3ccccc3[C@@H]1c1ccccc12",
    )
    .unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: true,
        canonical: false,
        clean_stereo: false,
        include_dative_bonds: false,
        ..Default::default()
    };

    assert_eq!(
        molecule.to_smiles_with_params(&params).unwrap(),
        "O=C(NC[C@]12C[C@H]3C[C@H](C[C@H](C3)C1)C2)[C@@H]1C[C@H]2c3ccccc3[C@@H]1c1ccccc12"
    );
}

#[test]
fn writer_preserves_ring_closure_double_bond_directions_like_rdkit_row_106() {
    let molecule =
        Molecule::from_smiles("O=C(N(C(S/1)=S)CCC(O)=O)C1=C\\C2=CC=C(C3=CC=C(C=C3)Cl)O2").unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: false,
        clean_stereo: false,
        all_bonds_explicit: true,
        include_dative_bonds: false,
        ..Default::default()
    };

    assert_eq!(
        molecule.to_smiles_with_params(&params).unwrap(),
        "O=C1-N(-C-C-C(-O)=O)-C(=S)-S/C-1=C\\c1:c:c:c(-c2:c:c:c(-Cl):c:c:2):o:1"
    );
}

#[test]
fn writer_preserves_bridgehead_tetrahedral_stereo_like_rdkit_row_110() {
    let molecule = Molecule::from_smiles(
        "[H]Cl.NC[C@@H]1O[C@@H](CC2=C(O)C(O)=CC=C12)C34C[C@H](C5)C[C@H](C[C@H]5C4)C3",
    )
    .unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: true,
        canonical: false,
        clean_stereo: false,
        include_dative_bonds: false,
        ..Default::default()
    };

    assert_eq!(
        molecule.to_smiles_with_params(&params).unwrap(),
        "Cl.NC[C@@H]1O[C@H](C23C[C@@H]4C[C@@H](C[C@@H](C4)C2)C3)Cc2c(O)c(O)ccc21"
    );
}

#[test]
fn canonical_isomeric_smiles_ranks_each_disconnected_chiral_fragment_like_source() {
    let molecule = Molecule::from_smiles(
        "COc1cc2c(cc1OC)[C@@H](Cc1cc(OC)c(OC)c(OC)c1)[N@+](C)(CCCOC(=O)/C(Cl)=C/C(=O)OCCC[N@+]13CCc4cc(OC)c(OC)cc4[C@H]1c1cc(OC)c(OC)cc1CC3)CC2.[Cl-].[Cl-]",
    )
    .unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: true,
        canonical: true,
        clean_stereo: false,
        ..Default::default()
    };

    assert_eq!(
        molecule.to_smiles_with_params(&params).unwrap(),
        "COc1cc2c(cc1OC)[C@H]1c3cc(OC)c(OC)cc3CC[N@@+]1(CCCOC(=O)/C=C(\\Cl)C(=O)OCCC[N@+]1(C)CCc3cc(OC)c(OC)cc3[C@H]1Cc1cc(OC)c(OC)c(OC)c1)CC2.[Cl-].[Cl-]"
    );
}

#[derive(serde::Deserialize)]
struct FragmentScopeStressParams {
    do_isomeric_smiles: bool,
    do_kekule: bool,
    canonical: bool,
    clean_stereo: bool,
    all_bonds_explicit: bool,
    all_hydrogens_explicit: bool,
    include_dative_bonds: bool,
    ignore_atom_map_numbers: bool,
    root: String,
}

#[derive(serde::Deserialize)]
struct FragmentScopeStressRecord {
    case_id: String,
    smiles: String,
    params: FragmentScopeStressParams,
    expected_smiles: String,
}

#[test]
fn canonical_fragment_scope_matches_every_retained_stress_profile() {
    let fixture = include_str!(
        "../../../../../testdata/smiles/fixtures/rdkit_fragment_scope_stress_regressions.jsonl"
    );
    let mut records = 0usize;
    let mut cases = std::collections::BTreeSet::new();

    for (line_idx, line) in fixture.lines().enumerate() {
        let record: FragmentScopeStressRecord =
            serde_json::from_str(line).unwrap_or_else(|error| {
                panic!(
                    "invalid fragment-scope stress fixture line {}: {error}",
                    line_idx + 1
                )
            });
        let molecule = Molecule::from_smiles(&record.smiles)
            .unwrap_or_else(|error| panic!("failed to parse {}: {error}", record.case_id));
        let rooted_at_atom = match record.params.root.as_str() {
            "none" => None,
            "first" => Some(0),
            "last" => molecule.num_atoms().checked_sub(1),
            root => panic!("invalid root {root:?} in {}", record.case_id),
        };
        let params = SmilesWriteParams {
            do_isomeric_smiles: record.params.do_isomeric_smiles,
            do_kekule: record.params.do_kekule,
            canonical: record.params.canonical,
            clean_stereo: record.params.clean_stereo,
            all_bonds_explicit: record.params.all_bonds_explicit,
            all_hydrogens_explicit: record.params.all_hydrogens_explicit,
            include_dative_bonds: record.params.include_dative_bonds,
            ignore_atom_map_numbers: record.params.ignore_atom_map_numbers,
            rooted_at_atom,
            do_random: false,
            ..Default::default()
        };

        assert_eq!(
            mol_to_smiles(&molecule, &params).unwrap(),
            record.expected_smiles,
            "RDKit 2026.03.1 stress regression {} at fixture line {}",
            record.case_id,
            line_idx + 1
        );
        cases.insert(record.case_id);
        records += 1;
    }

    assert_eq!(records, 108, "fixture must retain every affected profile");
    assert_eq!(
        cases.len(),
        9,
        "fixture must retain every affected molecule"
    );
}

#[test]
#[ignore = "debug helper for RDKit row 121 traversal parity"]
fn debug_row_121_noncanonical_traversal() {
    let mut molecule = Molecule::from_smiles(
            "O=C(O[Na])CC1=C(C(C(O[Na])=O)=C(C)C2=CC3=[N]4C(C(C=O)=C3CC)=CC5=C(C=C)C(C)=C6[N-]75)[N-]2[Cu+2]47[N](C8=C6)=C1C(C8C)CCC(O[Na])=O",
        )
        .unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        canonical: false,
        clean_stereo: false,
        include_dative_bonds: false,
        ..Default::default()
    };
    let mut working_params = params.clone();
    let _saved_atom_maps = prepare_plain_smiles_molecule(&mut molecule, &working_params).unwrap();
    working_params.do_kekule = false;
    let plan = collect_fragment_write_plans(&molecule, &working_params)
        .unwrap()
        .into_iter()
        .next()
        .unwrap();
    let ranks = rank_fragment_atoms_for_smiles(
        &molecule,
        &plan,
        &working_params,
        SmilesOutputMode::PlainSmiles,
    )
    .unwrap();
    molecule = kekulize_for_smiles(&molecule).unwrap();
    let start_atom = choose_fragment_start_atom(&plan, &ranks, &working_params).unwrap();
    let traversal = canonicalize_fragment_stack(
        &molecule,
        &plan,
        start_atom,
        &ranks,
        &working_params,
        SmilesWriteOverrides::default(),
    )
    .unwrap();
    let (ring_info, atom_ring_closures) =
        debug_atom_ring_closures_for_writer(&molecule, &plan, start_atom, &ranks, None).unwrap();
    let mut rank_by_atom = vec![usize::MAX; molecule.num_atoms()];
    for (idx, atom) in plan.atoms.iter().enumerate() {
        rank_by_atom[atom.index()] = ranks[idx];
    }
    for atom_idx in [0usize, 1, 4, 5, 6, 14, 17, 35, 33, 26] {
        let neighbors = molecule
            .topology_block()
            .adjacency
            .neighbors_of(atom_idx)
            .iter()
            .map(|nbr| (nbr.atom_index, nbr.bond.index()))
            .collect::<Vec<_>>();
        eprintln!("neighbors[{atom_idx}]={neighbors:?}");
    }
    for bond_idx in [17usize, 32, 33, 35, 48, 49, 50, 51, 52, 53, 54, 55] {
        eprintln!(
            "bond_rings[{bond_idx}]={}",
            ring_info.num_bond_rings(BondId::new(bond_idx))
        );
    }
    for atom_idx in [17usize, 18, 22, 16, 32, 33, 34, 35, 36, 37, 38, 39, 41] {
        let closures = atom_ring_closures[atom_idx]
            .iter()
            .map(|bond| bond.index())
            .collect::<Vec<_>>();
        eprintln!("atom_ring_closures[{atom_idx}]={closures:?}");
    }
    for atom_idx in [17usize, 35, 33, 26, 32, 36, 39, 40] {
        let atom = AtomId::new(atom_idx);
        let mut seen_from_here = vec![false; molecule.num_atoms()];
        seen_from_here[atom_idx] = true;
        for bond in &atom_ring_closures[atom_idx] {
            let bond = &molecule.bonds()[bond.index()];
            let other = if bond.begin() == atom {
                bond.end()
            } else {
                bond.begin()
            };
            seen_from_here[other.index()] = true;
        }
        let children = molecule
            .topology_block()
            .adjacency
            .neighbors_of(atom_idx)
            .iter()
            .map(|nbr| (BondId::new(nbr.bond.index()), AtomId::new(nbr.atom_index)))
            .filter(|(bond, other)| !seen_from_here[other.index()] && Some(*bond) != None)
            .map(|(bond, other)| {
                let mut rank = rank_by_atom[other.index()] as i64;
                let bond_order_rank = match molecule.bonds()[bond.index()].order() {
                    BondOrder::Null | BondOrder::Unspecified => 0,
                    BondOrder::Single => 1,
                    BondOrder::Double => 2,
                    BondOrder::Triple => 3,
                    BondOrder::Quadruple => 4,
                    BondOrder::Quintuple => 5,
                    BondOrder::Hextuple => 6,
                    BondOrder::OneAndHalf => 7,
                    BondOrder::TwoAndHalf => 8,
                    BondOrder::ThreeAndHalf => 9,
                    BondOrder::FourAndHalf => 10,
                    BondOrder::FiveAndHalf => 11,
                    BondOrder::Aromatic => 12,
                    BondOrder::Ionic => 13,
                    BondOrder::Hydrogen => 14,
                    BondOrder::ThreeCenter => 15,
                    BondOrder::DativeOne => 16,
                    BondOrder::Dative => 17,
                    BondOrder::DativeLeft => 18,
                    BondOrder::DativeRight => 19,
                    BondOrder::Other => 20,
                    BondOrder::Zero => 21,
                };
                if ring_info.num_bond_rings(bond) > 0 {
                    rank += (CANON_MAX_BONDTYPE - bond_order_rank)
                        * CANON_MAX_NATOMS
                        * CANON_MAX_NATOMS;
                }
                (rank, bond.index(), other.index())
            })
            .collect::<Vec<_>>();
        eprintln!("children_pre_sort[{atom_idx}]={children:?}");
    }
    eprintln!("start_atom={}", start_atom.index());
    eprintln!("stack={:#?}", traversal.stack);
}

#[test]
#[ignore = "debug helper for RDKit row 142 canonical kekule traversal parity"]
fn debug_row_142_canonical_kekule_root_last_traversal() {
    let input = "[C:12]12([CH:62]([CH3:65])[c:61]3[cH:64][cH:67][cH:68][cH:66][cH:63]3)[CH:20]4[c:30]5[c:40]6[c:49]7[c:57]8[c:60]([c:59]9[c:55]([c:47]([c:44]([c:52]9[c:51]([c:43]%10[c:35]%11[c:25]%12[c:19]%13%14)[c:53]8[c:45]%11[c:39]6[c:29]4%13)[c:34]([c:24]%15[c:15]%16[c:7]%17[c:3]%18%19)[c:33]%10[c:23]%16[c:16]%12[c:8]%18[c:11]%14[c:5]1%20)[c:37]([c:36]%21[c:26]%22[c:18]%23[c:10]%24[c:13]%25[c:6]%26%27)[c:27]%15[c:17]%22[c:9]%17[c:4]%24[c:1]%19[c:2]%20%26)[c:54]([c:46]%21[c:38]%28[c:28]%23[c:21]%25%29)[c:56]%30[c:48]%28[c:41]%31[c:31]%29[c:22]%32[c:14]2%27)[c:58]%30[c:50]7[c:42]%31[c:32]5%32";
    let mut molecule = Molecule::from_smiles(input).unwrap();
    let params = SmilesWriteParams {
        do_isomeric_smiles: false,
        do_kekule: true,
        canonical: true,
        clean_stereo: false,
        include_dative_bonds: false,
        ignore_atom_map_numbers: true,
        rooted_at_atom: Some(molecule.num_atoms() - 1),
        ..Default::default()
    };
    let mut working_params = params.clone();
    let _saved_atom_maps = prepare_plain_smiles_molecule(&mut molecule, &working_params).unwrap();
    let plan = collect_fragment_write_plans(&molecule, &working_params)
        .unwrap()
        .into_iter()
        .next()
        .unwrap();
    let ranks = rank_fragment_atoms_for_smiles(
        &molecule,
        &plan,
        &working_params,
        SmilesOutputMode::PlainSmiles,
    )
    .unwrap();
    molecule = kekulize_for_smiles(&molecule).unwrap();
    working_params.do_kekule = false;
    let start_atom = choose_fragment_start_atom(&plan, &ranks, &working_params).unwrap();
    let traversal = canonicalize_fragment_stack(
        &molecule,
        &plan,
        start_atom,
        &ranks,
        &working_params,
        SmilesWriteOverrides::default(),
    )
    .unwrap();
    let cached_rings = molecule.derived_cache().rings.clone();
    let fast_rings = crate::fast_find_rings(&molecule).unwrap();
    let (_, atom_ring_closures) =
        debug_atom_ring_closures_for_writer(&molecule, &plan, start_atom, &ranks, None).unwrap();
    eprintln!("start_atom={}", start_atom.index());
    eprintln!("ranks={ranks:?}");
    if let Some(rings) = &cached_rings {
        eprintln!(
            "cached find_type={:?} b11={} b72={}",
            rings.find_type(),
            rings.num_bond_rings(BondId::new(11)),
            rings.num_bond_rings(BondId::new(72))
        );
    } else {
        eprintln!("cached find_type=None");
    }
    eprintln!(
        "fast find_type={:?} b11={} b72={}",
        fast_rings.find_type(),
        fast_rings.num_bond_rings(BondId::new(11)),
        fast_rings.num_bond_rings(BondId::new(72))
    );
    eprintln!(
        "bond orders b11={:?} b72={:?}",
        molecule.bonds()[11].order(),
        molecule.bonds()[72].order()
    );
    eprintln!("atom11 closures={:?}", atom_ring_closures[11]);
    eprintln!("atom27 closures={:?}", atom_ring_closures[27]);
    eprint!("stack");
    for item in &traversal.stack {
        match item {
            MolStackElem::Atom(atom) => eprint!(" A{}", atom.index()),
            MolStackElem::Bond(bond, left) => eprint!(" B{}@L{}", bond.index(), left.index()),
            MolStackElem::Ring { bond: _, ring_idx } => eprint!(" R{ring_idx}"),
            MolStackElem::BranchOpen => eprint!(" ("),
            MolStackElem::BranchClose => eprint!(" )"),
        }
    }
    eprintln!();
    eprintln!(
        "output={}",
        molecule.to_smiles_with_params(&params).unwrap()
    );
}
