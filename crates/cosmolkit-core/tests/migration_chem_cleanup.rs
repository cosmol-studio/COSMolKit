use cosmolkit_core::{
    __migration_cleanup::{CleanupError, CleanupParams, cleanup},
    CanonicalRankError, ValenceError,
};
use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, StereoGroup, StereoGroupKind,
    SubstanceGroup, SubstanceGroupId, SubstanceGroupKind, TopologyBlock, TopologyValidationError,
};
use cosmolkit_types::{BondOrder, BondStereo, Element};

fn atom(atomic_number: u8) -> AtomSpec {
    AtomSpec::new(Element::from_atomic_number(atomic_number).expect("test element"))
}

fn bond(begin: usize, end: usize, order: BondOrder) -> BondSpec {
    BondSpec::new(AtomId::new(begin), AtomId::new(end), order)
}

fn topology(atom_specs: Vec<AtomSpec>, bond_specs: Vec<BondSpec>) -> TopologyBlock {
    let atoms = atom_specs
        .into_iter()
        .enumerate()
        .map(|(index, spec)| Atom::from_spec(AtomId::new(index), spec))
        .collect();
    let bonds = bond_specs
        .into_iter()
        .enumerate()
        .map(|(index, spec)| Bond::from_spec(BondId::new(index), spec))
        .collect();
    TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
}

fn params(charge_normalization: bool, organometallics: bool) -> CleanupParams {
    CleanupParams {
        charge_normalization,
        organometallics,
    }
}

fn nitro() -> TopologyBlock {
    topology(
        vec![atom(6), atom(7), atom(8), atom(8)],
        vec![
            bond(0, 1, BondOrder::Single),
            bond(1, 2, BondOrder::Double),
            bond(1, 3, BondOrder::Double),
        ],
    )
}

fn nitrogen_metal(candidate_atomic_number: u8) -> TopologyBlock {
    topology(
        vec![
            atom(7),
            atom(6),
            atom(6),
            atom(6),
            atom(candidate_atomic_number),
        ],
        vec![
            bond(0, 1, BondOrder::Single),
            bond(0, 2, BondOrder::Single),
            bond(0, 3, BondOrder::Single),
            bond(0, 4, BondOrder::Single),
        ],
    )
}

fn converted_bond_count(topology: &TopologyBlock) -> usize {
    topology
        .bonds
        .iter()
        .filter(|bond| bond.order() == BondOrder::Dative)
        .count()
}

#[test]
fn parameter_matrix_is_independent_deterministic_and_value_style() {
    let input = topology(
        vec![
            atom(6),
            atom(7),
            atom(8),
            atom(8),
            atom(7),
            atom(6),
            atom(6),
            atom(6),
            atom(26),
        ],
        vec![
            bond(0, 1, BondOrder::Single),
            bond(1, 2, BondOrder::Double),
            bond(1, 3, BondOrder::Double),
            bond(4, 5, BondOrder::Single),
            bond(4, 6, BondOrder::Single),
            bond(4, 7, BondOrder::Single),
            bond(4, 8, BondOrder::Single),
        ],
    );
    let snapshot = input.clone();

    for (charge_normalization, organometallics) in
        [(false, false), (true, false), (false, true), (true, true)]
    {
        let selected = params(charge_normalization, organometallics);
        let first = cleanup(&input, &selected).unwrap();
        let second = cleanup(&input, &selected).unwrap();
        assert_eq!(first, second);
        assert_eq!(input, snapshot);

        assert_eq!(
            first.atoms[1].formal_charge(),
            i8::from(charge_normalization)
        );
        assert_eq!(
            first.bonds[1].order(),
            if charge_normalization {
                BondOrder::Single
            } else {
                BondOrder::Double
            }
        );
        assert_eq!(
            first.bonds[6].order(),
            if organometallics {
                BondOrder::Dative
            } else {
                BondOrder::Single
            }
        );
    }

    let empty = TopologyBlock::default();
    assert_eq!(cleanup(&empty, &CleanupParams::default()).unwrap(), empty);
}

#[test]
fn nitrogen_cleanup_uses_first_neighbor_and_preserves_all_unrelated_state() {
    let mut input = nitro();
    input.atoms[0].set_prop("atom-note", "kept").unwrap();
    input.bonds[2].set_prop("bond-note", "kept").unwrap();
    input.substance_groups = vec![
        SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
            .with_atoms(vec![AtomId::new(0), AtomId::new(1)])
            .with_bonds(vec![BondId::new(0)]),
    ];
    input.stereo_groups = vec![
        StereoGroup::new(
            StereoGroupKind::Or,
            vec![AtomId::new(1)],
            vec![BondId::new(1)],
        )
        .with_id(17),
    ];
    input.validate().unwrap();
    let snapshot = input.clone();

    let result = cleanup(&input, &params(true, false)).unwrap();
    let mut expected = input.clone();
    expected.atoms[1].set_formal_charge(1);
    expected.atoms[2].set_formal_charge(-1);
    expected.bonds[1].set_order(BondOrder::Single);
    assert_eq!(result, expected);
    assert_eq!(input, snapshot);
    assert_eq!(result.atoms[3].formal_charge(), 0);
    assert_eq!(result.bonds[2].order(), BondOrder::Double);

    let reordered_neighbors = topology(
        vec![atom(6), atom(7), atom(8), atom(8)],
        vec![
            bond(0, 1, BondOrder::Single),
            bond(1, 3, BondOrder::Double),
            bond(1, 2, BondOrder::Double),
        ],
    );
    let reordered_result = cleanup(&reordered_neighbors, &params(true, false)).unwrap();
    assert_eq!(reordered_result.atoms[3].formal_charge(), -1);
    assert_eq!(reordered_result.atoms[2].formal_charge(), 0);
    assert_eq!(reordered_result.bonds[1].order(), BondOrder::Single);
    assert_eq!(reordered_result.bonds[2].order(), BondOrder::Double);
}

#[test]
fn nitrogen_second_pass_handles_n_triple_n_in_both_row_directions() {
    let forward = topology(
        vec![atom(6), atom(7), atom(7), atom(7)],
        vec![
            bond(0, 1, BondOrder::Single),
            bond(1, 2, BondOrder::Double),
            bond(2, 3, BondOrder::Triple),
        ],
    );
    let forward_result = cleanup(&forward, &params(true, false)).unwrap();
    assert_eq!(forward_result.atoms[2].formal_charge(), 1);
    assert_eq!(forward_result.atoms[3].formal_charge(), -1);
    assert_eq!(forward_result.bonds[2].order(), BondOrder::Double);

    let reversed = topology(
        vec![atom(7), atom(7), atom(7), atom(6)],
        vec![
            bond(1, 0, BondOrder::Triple),
            bond(1, 2, BondOrder::Double),
            bond(2, 3, BondOrder::Single),
        ],
    );
    let reversed_result = cleanup(&reversed, &params(true, false)).unwrap();
    assert_eq!(reversed_result.atoms[1].formal_charge(), 1);
    assert_eq!(reversed_result.atoms[0].formal_charge(), -1);
    assert_eq!(reversed_result.bonds[0].order(), BondOrder::Double);

    let isolated = topology(vec![atom(7), atom(7)], vec![bond(0, 1, BondOrder::Triple)]);
    assert_eq!(cleanup(&isolated, &params(true, false)).unwrap(), isolated);
}

#[test]
fn nitrogen_charge_and_normalized_forms_are_no_ops() {
    let charged = topology(
        vec![atom(6), atom(7).with_formal_charge(1), atom(8), atom(8)],
        vec![
            bond(0, 1, BondOrder::Single),
            bond(1, 2, BondOrder::Double),
            bond(1, 3, BondOrder::Double),
        ],
    );
    assert_eq!(cleanup(&charged, &params(true, false)).unwrap(), charged);

    let normalized = topology(
        vec![
            atom(6),
            atom(7).with_formal_charge(1),
            atom(8),
            atom(8).with_formal_charge(-1),
        ],
        vec![
            bond(0, 1, BondOrder::Single),
            bond(1, 2, BondOrder::Double),
            bond(1, 3, BondOrder::Single),
        ],
    );
    assert_eq!(
        cleanup(&normalized, &params(true, false)).unwrap(),
        normalized
    );
}

#[test]
fn phosphorus_positive_case_and_all_source_gates_are_exact() {
    let positive = topology(
        vec![atom(15), atom(8), atom(6), atom(6), atom(6)],
        vec![
            bond(0, 1, BondOrder::Double),
            bond(0, 2, BondOrder::Double),
            bond(0, 3, BondOrder::Single),
            bond(2, 4, BondOrder::Single),
        ],
    );
    let result = cleanup(&positive, &params(true, false)).unwrap();
    assert_eq!(result.atoms[0].formal_charge(), 1);
    assert_eq!(result.atoms[1].formal_charge(), -1);
    assert_eq!(result.bonds[0].order(), BondOrder::Single);
    assert_eq!(result.bonds[1].order(), BondOrder::Double);

    let missing_carbon_or_nitrogen = topology(
        vec![atom(15), atom(8), atom(16), atom(8)],
        vec![
            bond(0, 1, BondOrder::Double),
            bond(0, 2, BondOrder::Double),
            bond(0, 3, BondOrder::Single),
        ],
    );
    let missing_oxygen = topology(
        vec![atom(15), atom(6), atom(7), atom(6), atom(6)],
        vec![
            bond(0, 1, BondOrder::Double),
            bond(0, 2, BondOrder::Double),
            bond(0, 3, BondOrder::Single),
            bond(1, 4, BondOrder::Single),
        ],
    );
    let low_degree = topology(
        vec![atom(15), atom(8), atom(6), atom(6)],
        vec![
            bond(0, 1, BondOrder::Double),
            bond(0, 2, BondOrder::Double),
            bond(2, 3, BondOrder::Single),
        ],
    );
    let precharged = topology(
        vec![
            atom(15).with_formal_charge(1),
            atom(8),
            atom(6),
            atom(6),
            atom(6),
        ],
        vec![
            bond(0, 1, BondOrder::Double),
            bond(0, 2, BondOrder::Double),
            bond(0, 3, BondOrder::Single),
            bond(2, 4, BondOrder::Single),
        ],
    );
    let two_double_oxygens_fail_the_valence_five_gate = topology(
        vec![atom(15), atom(8), atom(8), atom(6), atom(6)],
        vec![
            bond(0, 1, BondOrder::Double),
            bond(0, 2, BondOrder::Double),
            bond(0, 3, BondOrder::Double),
            bond(3, 4, BondOrder::Single),
        ],
    );
    for unchanged in [
        missing_carbon_or_nitrogen,
        missing_oxygen,
        low_degree,
        precharged,
        two_double_oxygens_fail_the_valence_five_gate,
    ] {
        assert_eq!(
            cleanup(&unchanged, &params(true, false)).unwrap(),
            unchanged
        );
    }
}

#[test]
fn chlorine_bromine_and_iodine_cover_valence_three_five_and_seven() {
    for halogen in [17, 35, 53] {
        for double_oxygen_count in 1..=3 {
            let mut atoms = vec![atom(halogen)];
            atoms.extend((0..=double_oxygen_count).map(|_| atom(8)));
            let mut bonds = (1..=double_oxygen_count)
                .map(|oxygen| bond(0, oxygen, BondOrder::Double))
                .collect::<Vec<_>>();
            bonds.push(bond(0, double_oxygen_count + 1, BondOrder::Single));
            let input = topology(atoms, bonds);
            let result = cleanup(&input, &params(true, false)).unwrap();

            assert_eq!(result.atoms[0].formal_charge(), double_oxygen_count as i8);
            for oxygen in 1..=double_oxygen_count {
                assert_eq!(result.atoms[oxygen].formal_charge(), -1);
                assert_eq!(result.bonds[oxygen - 1].order(), BondOrder::Single);
            }
            assert_eq!(result.atoms[double_oxygen_count + 1].formal_charge(), 0);
            assert_eq!(result.bonds[double_oxygen_count].order(), BondOrder::Single);
        }
    }
}

#[test]
fn halogen_precharge_nonoxygen_neighbor_and_all_single_forms_are_no_ops() {
    for input in [
        topology(
            vec![atom(17).with_formal_charge(1), atom(8), atom(8)],
            vec![bond(0, 1, BondOrder::Double), bond(0, 2, BondOrder::Single)],
        ),
        topology(
            vec![atom(35), atom(8), atom(6)],
            vec![bond(0, 1, BondOrder::Double), bond(0, 2, BondOrder::Single)],
        ),
        topology(
            vec![atom(53), atom(8), atom(8), atom(8), atom(8), atom(8)],
            vec![
                bond(0, 1, BondOrder::Single),
                bond(0, 2, BondOrder::Single),
                bond(0, 3, BondOrder::Single),
                bond(0, 4, BondOrder::Single),
                bond(0, 5, BondOrder::Single),
            ],
        ),
    ] {
        assert_eq!(cleanup(&input, &params(true, false)).unwrap(), input);
    }
}

#[test]
fn metal_query_classification_matches_the_complete_source_exclusion_set() {
    let excluded = [
        0, 1, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 33, 34, 35, 36, 52, 53, 54, 85, 86,
    ];
    for atomic_number in 0..=118 {
        let input = nitrogen_metal(atomic_number);
        let result = cleanup(&input, &params(false, true)).unwrap();
        assert_eq!(
            converted_bond_count(&result),
            usize::from(!excluded.contains(&atomic_number)),
            "atomic number {atomic_number}"
        );
    }
}

#[test]
fn charge_effective_element_aromatic_degree_four_and_no_dative_gates_hold() {
    let n_plus = topology(
        vec![
            atom(7).with_formal_charge(1),
            atom(6),
            atom(6),
            atom(6),
            atom(26),
        ],
        vec![
            bond(0, 1, BondOrder::Single),
            bond(0, 2, BondOrder::Single),
            bond(0, 3, BondOrder::Single),
            bond(0, 4, BondOrder::Single),
        ],
    );
    let non_positive_effective = topology(
        vec![
            atom(6).with_formal_charge(7),
            atom(6),
            atom(6),
            atom(6),
            atom(6),
            atom(26),
        ],
        vec![
            bond(0, 1, BondOrder::Single),
            bond(0, 2, BondOrder::Single),
            bond(0, 3, BondOrder::Single),
            bond(0, 4, BondOrder::Single),
            bond(0, 5, BondOrder::Single),
        ],
    );
    for unchanged in [n_plus, non_positive_effective] {
        assert_eq!(
            cleanup(&unchanged, &params(false, true)).unwrap(),
            unchanged
        );
    }

    let aromatic = topology(
        vec![
            atom(6).with_aromatic(true),
            atom(6),
            atom(6),
            atom(6),
            atom(26),
        ],
        vec![
            bond(0, 1, BondOrder::Single),
            bond(0, 2, BondOrder::Single),
            bond(0, 3, BondOrder::Single),
            bond(0, 4, BondOrder::Single),
        ],
    );
    assert_eq!(
        cleanup(&aromatic, &params(false, true)).unwrap().bonds[3].order(),
        BondOrder::Dative
    );
    let mut non_aromatic = aromatic.clone();
    non_aromatic.atoms[0].set_aromatic(false);
    assert_eq!(
        cleanup(&non_aromatic, &params(false, true)).unwrap(),
        non_aromatic
    );

    for no_dative_element in [1, 2, 9, 10] {
        let input = topology(
            vec![atom(no_dative_element), atom(6), atom(6), atom(6), atom(26)],
            vec![
                bond(0, 1, BondOrder::Single),
                bond(0, 2, BondOrder::Single),
                bond(0, 3, BondOrder::Single),
                bond(0, 4, BondOrder::Single),
            ],
        );
        assert_eq!(cleanup(&input, &params(false, true)).unwrap(), input);
    }
}

fn two_metal_fixture(existing_order: Option<BondOrder>, reverse_metals: bool) -> TopologyBlock {
    let (first_metal, second_metal) = if reverse_metals { (29, 26) } else { (26, 29) };
    let mut atoms = vec![
        atom(7),
        atom(6),
        atom(6),
        atom(6),
        atom(first_metal),
        atom(second_metal),
    ];
    let mut bonds = vec![
        bond(0, 1, BondOrder::Single),
        bond(0, 2, BondOrder::Single),
        bond(0, 3, BondOrder::Single),
        bond(0, 4, BondOrder::Single),
        bond(0, 5, BondOrder::Single),
    ];
    if let Some(order) = existing_order {
        atoms.push(atom(30));
        bonds.push(bond(4, 6, order));
    }
    topology(atoms, bonds)
}

#[test]
fn dative_variant_reachability_matches_count_and_valence_source_boundaries() {
    for order in [BondOrder::Dative, BondOrder::DativeOne] {
        let input = two_metal_fixture(Some(order), false);
        let result = cleanup(&input, &params(false, true)).unwrap();
        assert_eq!(result.bonds[3].order(), BondOrder::Single, "{order:?}");
        assert_eq!(result.bonds[4].order(), BondOrder::Dative, "{order:?}");
        assert_eq!(result.bonds[4].begin(), AtomId::new(0));
        assert_eq!(result.bonds[4].end(), AtomId::new(5));
    }

    for order in [BondOrder::DativeLeft, BondOrder::DativeRight] {
        let input = two_metal_fixture(Some(order), false);
        assert!(matches!(
            cleanup(&input, &params(false, true)),
            Err(CleanupError::CanonicalRank(CanonicalRankError::Valence(
                ValenceError::BadBondType {
                    bond: Some(id),
                    order: actual,
                }
            ))) if id == BondId::new(5) && actual == order
        ));
    }
}

#[test]
fn canonical_tie_break_orientation_and_atom_renumbering_choose_same_metal() {
    let normal = cleanup(&two_metal_fixture(None, false), &params(false, true)).unwrap();
    assert_eq!(normal.bonds[3].order(), BondOrder::Single);
    assert_eq!(normal.bonds[4].order(), BondOrder::Dative);
    assert_eq!(
        normal.atoms[normal.bonds[4].end().index()].atomic_number(),
        29
    );
    assert_eq!(normal.bonds[4].begin(), AtomId::new(0));

    let reversed = cleanup(&two_metal_fixture(None, true), &params(false, true)).unwrap();
    assert_eq!(reversed.bonds[3].order(), BondOrder::Dative);
    assert_eq!(reversed.bonds[4].order(), BondOrder::Single);
    assert_eq!(
        reversed.atoms[reversed.bonds[3].end().index()].atomic_number(),
        29
    );
    assert_eq!(reversed.bonds[3].begin(), AtomId::new(0));
}

#[test]
fn sequential_ligands_observe_dative_counts_from_earlier_conversions() {
    let input = topology(
        vec![
            atom(7),
            atom(6),
            atom(6),
            atom(6),
            atom(7),
            atom(6),
            atom(6),
            atom(6),
            atom(26),
            atom(29),
        ],
        vec![
            bond(0, 1, BondOrder::Single),
            bond(0, 2, BondOrder::Single),
            bond(0, 3, BondOrder::Single),
            bond(0, 8, BondOrder::Single),
            bond(0, 9, BondOrder::Single),
            bond(4, 5, BondOrder::Single),
            bond(4, 6, BondOrder::Single),
            bond(4, 7, BondOrder::Single),
            bond(4, 8, BondOrder::Single),
            bond(4, 9, BondOrder::Single),
        ],
    );
    let result = cleanup(&input, &params(false, true)).unwrap();
    let iron_datives = result
        .bonds
        .iter()
        .filter(|bond| bond.order() == BondOrder::Dative && bond.end() == AtomId::new(8))
        .count();
    let copper_datives = result
        .bonds
        .iter()
        .filter(|bond| bond.order() == BondOrder::Dative && bond.end() == AtomId::new(9))
        .count();
    assert_eq!((iron_datives, copper_datives), (1, 1));
}

#[test]
fn no_work_fast_path_avoids_ranking_but_required_ranking_propagates_protocol_debt() {
    let atrop_spec = bond(0, 1, BondOrder::Single).with_stereo(BondStereo::AtropCw);
    let no_work = topology(vec![atom(6), atom(6)], vec![atrop_spec.clone()]);
    assert_eq!(cleanup(&no_work, &params(false, true)).unwrap(), no_work);

    let mut needs_ranking = nitrogen_metal(26);
    let offset = needs_ranking.atoms.len();
    needs_ranking
        .atoms
        .push(Atom::from_spec(AtomId::new(offset), atom(6)));
    needs_ranking
        .atoms
        .push(Atom::from_spec(AtomId::new(offset + 1), atom(6)));
    needs_ranking.bonds.push(Bond::from_spec(
        BondId::new(needs_ranking.bonds.len()),
        BondSpec::new(
            AtomId::new(offset),
            AtomId::new(offset + 1),
            BondOrder::Single,
        )
        .with_stereo(BondStereo::AtropCw),
    ));
    needs_ranking.adjacency =
        AdjacencyList::from_topology(needs_ranking.atoms.len(), &needs_ranking.bonds);
    needs_ranking.validate().unwrap();

    assert!(matches!(
        cleanup(&needs_ranking, &params(false, true)),
        Err(CleanupError::CanonicalRank(
            CanonicalRankError::ProtocolDebt { .. }
        ))
    ));
}

#[test]
fn malformed_topology_and_bad_valence_bond_return_exact_structured_errors() {
    let valid = nitro();
    let stale_adjacency = TopologyBlock {
        adjacency: AdjacencyList::from_topology(valid.atoms.len(), &[]),
        ..valid
    };
    assert!(matches!(
        cleanup(&stale_adjacency, &CleanupParams::default()),
        Err(CleanupError::InvalidTopology(
            TopologyValidationError::AdjacencyMismatch
        ))
    ));

    let bad_bond = topology(vec![atom(6), atom(6)], vec![bond(0, 1, BondOrder::Other)]);
    assert!(matches!(
        cleanup(&bad_bond, &params(false, true)),
        Err(CleanupError::Valence(ValenceError::BadBondType {
            bond: Some(id),
            order: BondOrder::Other,
        })) if id == BondId::new(0)
    ));
}
