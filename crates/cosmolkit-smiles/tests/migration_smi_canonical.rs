#[path = "../src/canonical_rank.rs"]
mod canonical_rank;

use canonical_rank::{CanonicalRankPolicy, rank_component_atoms, rank_component_atoms_with_policy};
use cosmolkit_model::{AtomId, BondId};
use cosmolkit_smiles::{
    SmilesParseError, SmilesParseParams, SmilesRecord, SmilesWriteParams, finalize_smiles_stereo,
    parse_smiles, write_smiles, write_smiles_with_params,
};

fn record(smiles: &str) -> SmilesRecord {
    parse_smiles(smiles, &Default::default())
        .unwrap_or_else(|error| panic!("failed to parse {smiles}: {error}"))
}

fn replace_atom_id(record: &mut SmilesRecord, row: usize, id: AtomId) {
    record.topology.atoms[row] = record.topology.atoms[row].clone().with_id(id);
}

fn replace_bond_id(record: &mut SmilesRecord, row: usize, id: BondId) {
    let original = record.topology.bonds[row].clone();
    record.topology.bonds[row] = original.clone().remapped(
        id,
        original.begin(),
        original.end(),
        original.stereo_atoms(),
    );
}

fn s33_standard_writer_input(input: &str) -> SmilesRecord {
    let parse_params = SmilesParseParams::default();
    let parsed = parse_smiles(input, &parse_params).unwrap();
    let remove_params = cosmolkit_core::RemoveHsParams {
        update_explicit_count: true,
        sanitize: parse_params.sanitize,
        ..cosmolkit_core::RemoveHsParams::default()
    };
    let prepared = cosmolkit_core::remove_hydrogens_with_params(
        parsed.topology,
        parsed.coordinates,
        parsed.properties,
        &remove_params,
    )
    .unwrap();
    finalize_smiles_stereo(
        SmilesRecord {
            topology: prepared.topology,
            coordinates: prepared.coordinates,
            properties: prepared.properties,
        },
        &parse_params,
    )
    .unwrap()
}

fn s33_standard_writer_params(isomeric: bool) -> SmilesWriteParams {
    SmilesWriteParams {
        do_isomeric_smiles: isomeric,
        do_kekule: false,
        canonical: true,
        clean_stereo: true,
        rooted_at_atom: None,
        all_bonds_explicit: false,
        all_hydrogens_explicit: false,
        include_dative_bonds: true,
        ignore_atom_map_numbers: false,
    }
}

fn s33_interleave_oc_co_components(record: &SmilesRecord) -> SmilesRecord {
    assert_eq!(record.topology.atoms.len(), 4);
    assert_eq!(record.topology.bonds.len(), 2);
    let new_to_old = [0, 2, 1, 3];
    let mut old_to_new = [0; 4];
    for (new_index, old_index) in new_to_old.into_iter().enumerate() {
        old_to_new[old_index] = new_index;
    }

    let mut interleaved = record.clone();
    interleaved.topology.atoms = new_to_old
        .into_iter()
        .enumerate()
        .map(|(new_index, old_index)| {
            record.topology.atoms[old_index]
                .clone()
                .with_id(AtomId::new(new_index))
        })
        .collect();
    for bond in &mut interleaved.topology.bonds {
        bond.set_endpoints(
            AtomId::new(old_to_new[bond.begin().index()]),
            AtomId::new(old_to_new[bond.end().index()]),
        );
    }
    interleaved.topology.adjacency = cosmolkit_model::AdjacencyList::from_topology(
        interleaved.topology.atoms.len(),
        &interleaved.topology.bonds,
    );
    interleaved.topology.validate().unwrap();
    interleaved
}

#[test]
fn component_atom_selection_rejects_duplicate_and_out_of_range_indices() {
    let parsed = record("CC.O");
    assert_eq!(
        rank_component_atoms(&parsed.topology, &[0, 0]),
        Err(cosmolkit_core::CanonicalRankError::DuplicateComponentAtom { atom_index: 0 })
    );
    assert_eq!(
        rank_component_atoms(&parsed.topology, &[3]),
        Err(
            cosmolkit_core::CanonicalRankError::ComponentAtomOutOfRange {
                atom_index: 3,
                atom_count: 3,
            }
        )
    );
}

#[test]
fn canonical_writer_rejects_duplicate_and_out_of_range_model_row_ids() {
    let mut duplicate_atom = record("CC");
    replace_atom_id(&mut duplicate_atom, 1, AtomId::new(0));
    assert_eq!(
        write_smiles(&duplicate_atom),
        Err(SmilesParseError::Model(
            "atom at position 1 has id 0, expected 1".to_owned()
        ))
    );

    let mut out_of_range_atom = record("CC");
    replace_atom_id(&mut out_of_range_atom, 1, AtomId::new(4));
    assert_eq!(
        write_smiles(&out_of_range_atom),
        Err(SmilesParseError::Model(
            "atom at position 1 has id 4, expected 1".to_owned()
        ))
    );

    let mut duplicate_bond = record("CCC");
    replace_bond_id(&mut duplicate_bond, 1, BondId::new(0));
    assert_eq!(
        write_smiles(&duplicate_bond),
        Err(SmilesParseError::Model(
            "bond at position 1 has id 0, expected 1".to_owned()
        ))
    );

    let mut out_of_range_bond = record("CC");
    replace_bond_id(&mut out_of_range_bond, 0, BondId::new(4));
    assert_eq!(
        write_smiles(&out_of_range_bond),
        Err(SmilesParseError::Model(
            "bond at position 0 has id 4, expected 0".to_owned()
        ))
    );
}

#[test]
fn canonical_writer_rejects_duplicate_edges_before_fragment_mapping() {
    let mut duplicate_edge = record("CC");
    let original = duplicate_edge.topology.bonds[0].clone();
    duplicate_edge
        .topology
        .bonds
        .push(original.clone().remapped(
            BondId::new(1),
            original.end(),
            original.begin(),
            original.stereo_atoms(),
        ));
    assert_eq!(
        write_smiles(&duplicate_edge),
        Err(SmilesParseError::Model(
            "adjacency does not match topology".to_owned()
        ))
    );
}

#[test]
fn canonical_writer_keeps_pinned_disconnected_component_outputs() {
    for (input, expected) in [("OCC", "CCO"), ("N.CCO", "CCO.N")] {
        assert_eq!(write_smiles(&record(input)).unwrap(), expected, "{input}");
    }
}

fn rank_all_atoms(record: &SmilesRecord, policy: CanonicalRankPolicy) -> Vec<usize> {
    let atoms = (0..record.topology.atoms.len()).collect::<Vec<_>>();
    rank_component_atoms_with_policy(&record.topology, &atoms, policy)
        .expect("canonical rank should accept the parsed topology")
}

#[test]
fn canonical_rank_policy_toggles_isotope_chirality_groups_and_atom_maps() {
    let isotope = record("[13CH3]OC");
    let mut isotope_policy = CanonicalRankPolicy::default();
    isotope_policy.break_ties = false;
    let isotope_ranks = rank_all_atoms(&isotope, isotope_policy);
    assert_ne!(isotope_ranks[0], isotope_ranks[2]);
    isotope_policy.include_isotopes = false;
    let isotope_ranks_without_isotopes = rank_all_atoms(&isotope, isotope_policy);
    assert_eq!(
        isotope_ranks_without_isotopes[0],
        isotope_ranks_without_isotopes[2]
    );

    let chiral = record("F[C@H](Cl)OC(F)Cl");
    let mut chirality_policy = CanonicalRankPolicy::default();
    chirality_policy.break_ties = false;
    let chiral_ranks = rank_all_atoms(&chiral, chirality_policy);
    assert_ne!(chiral_ranks[1], chiral_ranks[4]);
    chirality_policy.include_chirality = false;
    let chiral_ranks_without_chirality = rank_all_atoms(&chiral, chirality_policy);
    assert_eq!(
        chiral_ranks_without_chirality[1],
        chiral_ranks_without_chirality[4]
    );

    let mapped = record("[F:1]C([F:2])O");
    let mut map_policy = CanonicalRankPolicy::default();
    map_policy.break_ties = false;
    let map_ranks = rank_all_atoms(&mapped, map_policy);
    assert_ne!(map_ranks[0], map_ranks[2]);
    map_policy.include_atom_maps = false;
    let map_ranks_without_maps = rank_all_atoms(&mapped, map_policy);
    assert_eq!(map_ranks_without_maps[0], map_ranks_without_maps[2]);

    let grouped_input = "F[C@H](Cl)O[C@H](F)Cl |o1:1|";
    let grouped = record(grouped_input);
    let ungrouped = record("F[C@H](Cl)O[C@H](F)Cl");
    let mut group_policy = CanonicalRankPolicy::default();
    group_policy.break_ties = false;
    group_policy.include_stereo_groups = false;
    assert_eq!(
        rank_all_atoms(&grouped, group_policy),
        rank_all_atoms(&ungrouped, group_policy)
    );
    group_policy.include_stereo_groups = true;
    let group_ranks = rank_all_atoms(&grouped, group_policy);
    assert_ne!(group_ranks[1], group_ranks[4]);
}

#[test]
fn canonical_writer_projects_isomeric_rank_flags_to_pinned_outputs() {
    let mut params = SmilesWriteParams::default();
    for (input, isomeric, nonisomeric) in [
        ("FOCN[15F]", "FOCN[15F]", "FNCOF"),
        ("[15F]OCNF", "FNCO[15F]", "FNCOF"),
        (
            "FC(Cl)OCN[C@H](F)Cl",
            "FC(Cl)OCN[C@H](F)Cl",
            "FC(Cl)NCOC(F)Cl",
        ),
        (
            "FC(Cl)NCO[C@H](F)Cl",
            "FC(Cl)NCO[C@H](F)Cl",
            "FC(Cl)NCOC(F)Cl",
        ),
    ] {
        params.do_isomeric_smiles = true;
        assert_eq!(
            write_smiles_with_params(&record(input), &params).unwrap(),
            isomeric,
            "isomeric {input}"
        );

        params.do_isomeric_smiles = false;
        assert_eq!(
            write_smiles_with_params(&record(input), &params).unwrap(),
            nonisomeric,
            "nonisomeric {input}"
        );
    }
}

#[test]
fn canonical_writer_ignore_atom_maps_changes_traversal_and_preserves_labels() {
    let mapped = record("[NH2:1]c1ccccc1");
    let mut params = SmilesWriteParams::default();
    assert_eq!(
        write_smiles_with_params(&mapped, &params).unwrap(),
        "c1ccc([NH2:1])cc1"
    );

    params.ignore_atom_map_numbers = true;
    assert_eq!(
        write_smiles_with_params(&mapped, &params).unwrap(),
        "[NH2:1]c1ccccc1"
    );
}

#[test]
fn canonical_cycle_discovery_resolves_symmetric_ranks_and_preserves_closures() {
    let symmetric_ring = record("C1CCCCC1");
    let mut unresolved_policy = CanonicalRankPolicy::default();
    unresolved_policy.break_ties = false;
    let unresolved_ranks = rank_all_atoms(&symmetric_ring, unresolved_policy);
    assert!(unresolved_ranks.windows(2).all(|pair| pair[0] == pair[1]));

    let resolved_ranks = rank_all_atoms(&symmetric_ring, CanonicalRankPolicy::default());
    let mut unique_ranks = resolved_ranks.clone();
    unique_ranks.sort_unstable();
    unique_ranks.dedup();
    assert_eq!(unique_ranks.len(), symmetric_ring.topology.atoms.len());
    assert_eq!(write_smiles(&symmetric_ring).unwrap(), "C1CCCCC1");

    // Pinned catch_tests.cpp #5585 asserts this exact output and traversal order.
    assert_eq!(write_smiles(&record("OC1CCCN1")).unwrap(), "OC1CCCN1");

    // The pinned RDKit oracle returns this closure ordering for a three-ring cage.
    assert_eq!(
        write_smiles(&record("C1C2C3C1C2C3")).unwrap(),
        "C1C2C3CC2C13"
    );
}

#[test]
fn canonical_stack_orders_branch_children_from_the_rank_map() {
    for (input, expected) in [
        ("CC(CO)(N)F", "CC(N)(F)CO"),
        ("FC(C)(N)O", "CC(N)(O)F"),
        ("CC(=O)N(C)O", "CC(=O)N(C)O"),
    ] {
        assert_eq!(write_smiles(&record(input)).unwrap(), expected, "{input}");
    }
}

#[test]
fn canonical_writer_preserves_pinned_tetrahedral_ring_stereo_permutations() {
    let mut params = SmilesWriteParams::default();
    params.canonical = true;
    params.do_isomeric_smiles = true;
    params.clean_stereo = true;

    for (input, expected) in [
        (
            "CC1CCC(CC1)NO[C@@H]1CC[C@H](C)CC1",
            "CC1CCC(NO[C@H]2CC[C@@H](C)CC2)CC1",
        ),
        (
            "CC1CCC(CC1)ON[C@@H]1CC[C@H](C)CC1",
            "CC1CCC(ON[C@H]2CC[C@@H](C)CC2)CC1",
        ),
    ] {
        assert_eq!(
            write_smiles_with_params(&record(input), &params).unwrap(),
            expected,
            "{input}"
        );
    }
}

#[test]
fn canonical_writer_decodes_signed_ring_relative_stereo_ids_in_source_order() {
    for (input, first_relation, second_relation, expected) in [
        (
            "C1[C@H](F)CC[C@H](Cl)C1",
            "6",
            "2",
            "F[C@H]1CC[C@@H](Cl)CC1",
        ),
        (
            "C1[C@H](F)CC[C@@H](Cl)C1",
            "-6",
            "-2",
            "F[C@H]1CC[C@H](Cl)CC1",
        ),
    ] {
        let mut parsed = record(input);
        parsed.topology.atoms[1]
            .set_prop("_ringStereoAtoms", first_relation)
            .unwrap();
        parsed.topology.atoms[5]
            .set_prop("_ringStereoAtoms", second_relation)
            .unwrap();

        assert_eq!(write_smiles(&parsed).unwrap(), expected, "{input}");
    }
}

#[test]
fn canonical_writer_reports_source_bad_any_cast_for_malformed_ring_stereo_property() {
    let source_error = SmilesParseError::WriterStereo(
        "`_ringStereoAtoms` is not a signed source INT_VECT value (bad_any_cast)".to_owned(),
    );

    for encoded in ["", "6,,2", "not-an-int"] {
        let mut parsed = record("C1[C@H](F)CC[C@H](Cl)C1");
        parsed.topology.atoms[1]
            .set_prop("_ringStereoAtoms", encoded)
            .unwrap();

        assert_eq!(
            write_smiles(&parsed),
            Err(source_error.clone()),
            "{encoded:?}"
        );
    }
}

#[test]
fn canonical_writer_does_not_propagate_from_a_broken_relative_stereo_center() {
    let mut parsed = record("C1[C@H](F)CC[C@H](Cl)C1");
    parsed.topology.atoms[1]
        .set_computed_prop("_ringStereoAtoms", "6")
        .unwrap();
    parsed.topology.atoms[5]
        .set_computed_prop("_ringStereoAtoms", "2")
        .unwrap();
    parsed.topology.atoms[1]
        .set_prop("_brokenChirality", "1")
        .unwrap();

    let mut params = SmilesWriteParams::default();
    params.canonical = false;
    params.do_isomeric_smiles = true;
    params.clean_stereo = true;
    params.rooted_at_atom = Some(AtomId::new(1));

    // Pinned RDKit MolToSmiles, rootedAtAtom=1, emits this when the source
    // center has _brokenChirality, even though it retains its ring reference.
    assert_eq!(
        write_smiles_with_params(&parsed, &params).unwrap(),
        "C1(F)CC[C@H](Cl)CC1"
    );
}

#[test]
fn canonical_writer_maps_ring_relative_references_across_component_order() {
    for (input, first_center, first_reference, second_center, second_reference, root, expected) in [
        (
            "N.C1[C@H](F)CC[C@H](Cl)C1",
            2,
            "7",
            6,
            "3",
            2,
            "N.[C@H]1(F)CC[C@@H](Cl)CC1",
        ),
        (
            "C1[C@H](F)CC[C@H](Cl)C1.N",
            1,
            "6",
            5,
            "2",
            1,
            "[C@H]1(F)CC[C@@H](Cl)CC1.N",
        ),
    ] {
        let mut parsed = record(input);
        parsed.topology.atoms[first_center]
            .set_computed_prop("_ringStereoAtoms", first_reference)
            .unwrap();
        parsed.topology.atoms[second_center]
            .set_computed_prop("_ringStereoAtoms", second_reference)
            .unwrap();
        assert_eq!(
            parsed.topology.atoms[first_center].prop("_ringStereoAtoms"),
            Some(first_reference),
            "first parsed relation for {input}"
        );
        assert_eq!(
            parsed.topology.atoms[second_center].prop("_ringStereoAtoms"),
            Some(second_reference),
            "second parsed relation for {input}"
        );

        let mut params = SmilesWriteParams::default();
        params.canonical = false;
        params.do_isomeric_smiles = true;
        params.clean_stereo = true;
        params.rooted_at_atom = Some(AtomId::new(root));

        // Pinned RDKit 2026.03.1, same raw input and explicit writer profile:
        // both disconnected component order and source-global signed refs are
        // preserved while the ring is prepared as a selected fragment.
        assert_eq!(
            write_smiles_with_params(&parsed, &params).unwrap(),
            expected,
            "{input}"
        );
    }
}

#[test]
fn canonical_writer_orders_s33_components_with_the_pinned_standard_profile_and_preserves_input() {
    // Pinned RDKit 2026.03.1 standard-writer text rows from
    // S33_component_order_oracle.md; map vectors are observed by the private
    // writer owner test because the public standard API returns text only.
    for (input, isomeric, expected) in [
        ("C[C@H](F)Cl.[13CH3]O", true, "C[C@H](F)Cl.[13CH3]O"),
        ("[13CH3]O.C[C@H](F)Cl", true, "C[C@H](F)Cl.[13CH3]O"),
        ("C[C@H](F)Cl.[13CH3]O", false, "CC(F)Cl.CO"),
        ("[13CH3]O.C[C@H](F)Cl", false, "CC(F)Cl.CO"),
        (
            "[CH3:10][C@H:11]([F:12])[Cl:13].[13CH3:20][O:21]",
            true,
            "[13CH3:20][O:21].[CH3:10][C@H:11]([F:12])[Cl:13]",
        ),
        (
            "[13CH3:20][O:21].[CH3:10][C@H:11]([F:12])[Cl:13]",
            true,
            "[13CH3:20][O:21].[CH3:10][C@H:11]([F:12])[Cl:13]",
        ),
        (
            "[CH3:10][C@H:11]([F:12])[Cl:13].[13CH3:20][O:21]",
            false,
            "[CH3:10][CH:11]([F:12])[Cl:13].[CH3:20][O:21]",
        ),
        (
            "[13CH3:20][O:21].[CH3:10][C@H:11]([F:12])[Cl:13]",
            false,
            "[CH3:10][CH:11]([F:12])[Cl:13].[CH3:20][O:21]",
        ),
    ] {
        let input_record = s33_standard_writer_input(input);
        let before = input_record.clone();
        let params = s33_standard_writer_params(isomeric);
        assert_eq!(
            write_smiles_with_params(&input_record, &params).unwrap(),
            expected,
            "{input}; isomeric={isomeric}"
        );
        assert_eq!(input_record, before, "{input}; isomeric={isomeric}");
    }

    let source_order = s33_standard_writer_input("OC.CO");
    let before = source_order.clone();
    assert_eq!(
        write_smiles_with_params(&source_order, &s33_standard_writer_params(true)).unwrap(),
        "CO.CO",
        "equal-text components in source order"
    );
    assert_eq!(source_order, before, "equal-text source-order input");

    let interleaved = s33_interleave_oc_co_components(&source_order);
    let before = interleaved.clone();
    assert_eq!(
        write_smiles_with_params(&interleaved, &s33_standard_writer_params(true)).unwrap(),
        "CO.CO",
        "equal-text components after [0,2,1,3] input permutation"
    );
    assert_eq!(interleaved, before, "equal-text interleaved input");
}
