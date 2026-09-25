#![cfg(feature = "op-contracts-strict")]

use cosmolkit_core::{
    __migration_hydrogens::{
        remove_hydrogen_candidates, remove_hydrogen_candidates_with_query_state,
    },
    RemoveHsParams,
};
use cosmolkit_model::{
    Atom, AtomId, AtomQueryPredicate, AtomSpec, Bond, BondDirection, BondId, BondQueryPredicate,
    BondSpec, BondStereo, ChiralTag, QueryAtom, QueryBond, QueryNode, QueryStateRef,
    SGroupAttachPoint, SGroupBondRole, SGroupCState, SubstanceGroup, SubstanceGroupId,
    SubstanceGroupKind, TopologyBlock,
};
use cosmolkit_types::{BondOrder, Element};

fn atom(index: usize) -> AtomId {
    AtomId::new(index)
}

fn bond_id(index: usize) -> BondId {
    BondId::new(index)
}

fn bond(begin: usize, end: usize) -> BondSpec {
    BondSpec::new(atom(begin), atom(end), BondOrder::Single)
}

fn topology(
    atom_specs: Vec<AtomSpec>,
    bond_specs: Vec<BondSpec>,
    groups: Vec<SubstanceGroup>,
) -> TopologyBlock {
    let atoms = atom_specs
        .into_iter()
        .enumerate()
        .map(|(index, spec)| Atom::from_spec(atom(index), spec))
        .collect();
    let bonds = bond_specs
        .into_iter()
        .enumerate()
        .map(|(index, spec)| Bond::from_spec(bond_id(index), spec))
        .collect();
    TopologyBlock::try_from_parts(atoms, bonds, groups, Vec::new()).unwrap()
}

fn selected(topology: &TopologyBlock, params: &RemoveHsParams) -> Vec<usize> {
    remove_hydrogen_candidates(topology, params)
        .into_iter()
        .map(AtomId::index)
        .collect()
}

fn carbon_hydrogen(hydrogen: AtomSpec) -> TopologyBlock {
    topology(
        vec![AtomSpec::new(Element::C), hydrogen],
        vec![bond(0, 1)],
        Vec::new(),
    )
}

#[test]
fn query_overlay_hydrogen_selection_reads_current_isotope_not_old_carrier() {
    // The existing source-backed isotope rule keeps isotope-labelled hydrogen
    // by default. Query identity comes from the overlay, isotope from topology.
    let old = carbon_hydrogen(AtomSpec::new(Element::H));
    let atoms: Vec<_> = old
        .atoms
        .iter()
        .map(|atom| {
            QueryAtom::from_carrier_parts(
                atom.clone(),
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(atom.atomic_number())),
            )
        })
        .collect();
    let bonds = vec![QueryBond::from_carrier_parts(
        old.bonds[0].clone(),
        QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
    )];
    let overlay = QueryStateRef::try_for_topology(&atoms, &bonds, &old).unwrap();
    let current = carbon_hydrogen(AtomSpec::new(Element::H).with_isotope(2));
    overlay.validate_for_topology(&current).unwrap();
    let defaults = RemoveHsParams::default();
    assert_eq!(selected(&old, &defaults), vec![1]);
    assert!(selected(&current, &defaults).is_empty());
    assert!(
        remove_hydrogen_candidates_with_query_state(&current, &defaults, Some(overlay)).is_empty()
    );
    let remove_isotopes = RemoveHsParams {
        remove_isotopes: true,
        ..defaults
    };
    assert_eq!(
        remove_hydrogen_candidates_with_query_state(&current, &remove_isotopes, Some(overlay)),
        vec![atom(1)]
    );
}

#[test]
fn defaults_match_all_remove_hs_parameter_source_fields() {
    let params = RemoveHsParams::default();
    assert!(!params.remove_degree_zero);
    assert!(!params.remove_higher_degrees);
    assert!(!params.remove_only_h_neighbors);
    assert!(!params.remove_isotopes);
    assert!(!params.remove_and_track_isotopes);
    assert!(!params.remove_dummy_neighbors);
    assert!(!params.remove_defining_bond_stereo);
    assert!(params.remove_with_wedged_bond);
    assert!(!params.remove_with_query);
    assert!(params.remove_mapped);
    assert!(params.remove_in_sgroups);
    assert!(params.show_warnings);
    assert!(params.remove_nonimplicit);
    assert!(!params.update_explicit_count);
    assert!(!params.remove_hydrides);
    assert!(!params.remove_nontetrahedral_neighbors);
    assert!(params.sanitize);
}

#[test]
fn selection_excludes_non_h_and_preserves_h_atom_order_without_mutation() {
    let source = topology(
        vec![
            AtomSpec::new(Element::H),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::H),
            AtomSpec::new(Element::H),
        ],
        vec![bond(0, 1), bond(1, 2), bond(1, 3)],
        Vec::new(),
    );
    let snapshot = source.clone();
    assert_eq!(selected(&source, &RemoveHsParams::default()), vec![0, 2, 3]);
    assert_eq!(source, snapshot);
}

#[test]
fn degree_zero_and_higher_degree_flags_have_source_polarity() {
    let isolated = topology(vec![AtomSpec::new(Element::H)], Vec::new(), Vec::new());
    assert!(selected(&isolated, &RemoveHsParams::default()).is_empty());
    assert_eq!(
        selected(
            &isolated,
            &RemoveHsParams {
                remove_degree_zero: true,
                ..Default::default()
            }
        ),
        vec![0]
    );

    let bridging = topology(
        vec![
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::H),
            AtomSpec::new(Element::C),
        ],
        vec![bond(0, 1), bond(1, 2)],
        Vec::new(),
    );
    assert!(selected(&bridging, &RemoveHsParams::default()).is_empty());
    assert_eq!(
        selected(
            &bridging,
            &RemoveHsParams {
                remove_higher_degrees: true,
                ..Default::default()
            }
        ),
        vec![1]
    );
}

#[test]
fn only_h_neighbor_flag_controls_h2_without_becoming_a_neighbor_requirement() {
    let h2 = topology(
        vec![AtomSpec::new(Element::H), AtomSpec::new(Element::H)],
        vec![bond(0, 1)],
        Vec::new(),
    );
    assert!(selected(&h2, &RemoveHsParams::default()).is_empty());
    assert_eq!(
        selected(
            &h2,
            &RemoveHsParams {
                remove_only_h_neighbors: true,
                ..Default::default()
            }
        ),
        vec![0, 1]
    );
    let carbon_h = carbon_hydrogen(AtomSpec::new(Element::H));
    assert_eq!(
        selected(
            &carbon_h,
            &RemoveHsParams {
                remove_only_h_neighbors: true,
                ..Default::default()
            }
        ),
        vec![1]
    );
}

#[test]
fn isotope_flags_are_independent_candidate_permissions() {
    let isotope = carbon_hydrogen(AtomSpec::new(Element::H).with_isotope(2));
    assert!(selected(&isotope, &RemoveHsParams::default()).is_empty());
    for params in [
        RemoveHsParams {
            remove_isotopes: true,
            ..Default::default()
        },
        RemoveHsParams {
            remove_and_track_isotopes: true,
            ..Default::default()
        },
    ] {
        assert_eq!(selected(&isotope, &params), vec![1]);
    }
    let zero = carbon_hydrogen(AtomSpec::new(Element::H).with_isotope(0));
    assert_eq!(selected(&zero, &RemoveHsParams::default()), vec![1]);
}

#[test]
fn query_map_and_source_implicit_flags_follow_source_values() {
    let query = carbon_hydrogen(AtomSpec::new(Element::H));
    let query_atoms = vec![
        QueryAtom::from_carrier_parts(
            query.atoms[0].clone(),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
        ),
        QueryAtom::from_parts(
            query.atoms[1].clone(),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(1)),
        ),
    ];
    let query_bonds = vec![QueryBond::from_carrier_parts(
        query.bonds[0].clone(),
        QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
    )];
    let query_state = QueryStateRef::try_for_topology(&query_atoms, &query_bonds, &query).unwrap();
    assert!(query_state.atom_has_query(atom(1)));
    let selected_with_state = |params: &RemoveHsParams| {
        remove_hydrogen_candidates_with_query_state(&query, params, Some(query_state))
            .into_iter()
            .map(AtomId::index)
            .collect::<Vec<_>>()
    };
    assert!(selected_with_state(&RemoveHsParams::default()).is_empty());
    assert_eq!(
        selected_with_state(&RemoveHsParams {
            remove_with_query: true,
            ..Default::default()
        }),
        vec![1]
    );

    let mapped = carbon_hydrogen(AtomSpec::new(Element::H).with_atom_map(7));
    assert_eq!(selected(&mapped, &RemoveHsParams::default()), vec![1]);
    assert!(
        selected(
            &mapped,
            &RemoveHsParams {
                remove_mapped: false,
                ..Default::default()
            }
        )
        .is_empty()
    );
    let zero_map = carbon_hydrogen(AtomSpec::new(Element::H).with_atom_map(0));
    assert_eq!(
        selected(
            &zero_map,
            &RemoveHsParams {
                remove_mapped: false,
                ..Default::default()
            }
        ),
        vec![1]
    );

    let explicit = carbon_hydrogen(AtomSpec::new(Element::H));
    let implicit = carbon_hydrogen(AtomSpec::new(Element::H).with_implicit_hydrogen(true));
    let implicit_only = RemoveHsParams {
        remove_nonimplicit: false,
        ..Default::default()
    };
    assert!(selected(&explicit, &implicit_only).is_empty());
    assert_eq!(selected(&implicit, &implicit_only), vec![1]);
}

#[test]
fn hydride_dummy_and_nontetrahedral_neighbor_flags_have_source_polarity() {
    let hydride = carbon_hydrogen(AtomSpec::new(Element::H).with_formal_charge(-1));
    assert!(selected(&hydride, &RemoveHsParams::default()).is_empty());
    assert_eq!(
        selected(
            &hydride,
            &RemoveHsParams {
                remove_hydrides: true,
                ..Default::default()
            }
        ),
        vec![1]
    );

    let dummy = topology(
        vec![AtomSpec::new(Element::DUMMY), AtomSpec::new(Element::H)],
        vec![bond(0, 1)],
        Vec::new(),
    );
    assert!(selected(&dummy, &RemoveHsParams::default()).is_empty());
    assert_eq!(
        selected(
            &dummy,
            &RemoveHsParams {
                remove_dummy_neighbors: true,
                ..Default::default()
            }
        ),
        vec![1]
    );

    for tag in [
        ChiralTag::SquarePlanar,
        ChiralTag::TrigonalBipyramidal,
        ChiralTag::Octahedral,
    ] {
        let stereo = topology(
            vec![
                AtomSpec::new(Element::PT).with_chiral_tag(tag),
                AtomSpec::new(Element::H),
            ],
            vec![bond(0, 1)],
            Vec::new(),
        );
        assert!(selected(&stereo, &RemoveHsParams::default()).is_empty());
        assert_eq!(
            selected(
                &stereo,
                &RemoveHsParams {
                    remove_nontetrahedral_neighbors: true,
                    ..Default::default()
                }
            ),
            vec![1]
        );
    }
}

#[test]
fn wedge_gate_distinguishes_source_wedges_from_other_directions() {
    for direction in [BondDirection::BeginWedge, BondDirection::BeginDash] {
        let source = topology(
            vec![AtomSpec::new(Element::C), AtomSpec::new(Element::H)],
            vec![bond(0, 1).with_direction(direction)],
            Vec::new(),
        );
        assert!(
            selected(
                &source,
                &RemoveHsParams {
                    remove_with_wedged_bond: false,
                    ..Default::default()
                }
            )
            .is_empty()
        );
        assert_eq!(selected(&source, &RemoveHsParams::default()), vec![1]);
    }
    let directional = topology(
        vec![AtomSpec::new(Element::C), AtomSpec::new(Element::H)],
        vec![bond(0, 1).with_direction(BondDirection::EndUpRight)],
        Vec::new(),
    );
    assert_eq!(
        selected(
            &directional,
            &RemoveHsParams {
                remove_with_wedged_bond: false,
                ..Default::default()
            }
        ),
        vec![1]
    );
}

#[test]
fn defining_double_bond_stereo_gate_matches_degree_two_neighbor_rule() {
    let defining = topology(
        vec![
            AtomSpec::new(Element::H),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
        ],
        vec![
            bond(0, 1),
            BondSpec::new(atom(1), atom(2), BondOrder::Double).with_stereo(BondStereo::Z),
        ],
        Vec::new(),
    );
    assert!(selected(&defining, &RemoveHsParams::default()).is_empty());
    assert_eq!(
        selected(
            &defining,
            &RemoveHsParams {
                remove_defining_bond_stereo: true,
                ..Default::default()
            }
        ),
        vec![0]
    );

    let any = topology(
        vec![
            AtomSpec::new(Element::H),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
        ],
        vec![
            bond(0, 1),
            BondSpec::new(atom(1), atom(2), BondOrder::Double).with_stereo(BondStereo::Any),
        ],
        Vec::new(),
    );
    assert_eq!(selected(&any, &RemoveHsParams::default()), vec![0]);

    let directed_h = topology(
        vec![
            AtomSpec::new(Element::H),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
        ],
        vec![
            bond(0, 1).with_direction(BondDirection::EndUpRight),
            BondSpec::new(atom(1), atom(2), BondOrder::Double),
        ],
        Vec::new(),
    );
    assert!(selected(&directed_h, &RemoveHsParams::default()).is_empty());
}

#[test]
fn sgroup_true_allows_ordinary_membership_but_protects_special_roles() {
    let ordinary = SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
        .with_atoms(vec![atom(0), atom(1)]);
    let source = topology(
        vec![AtomSpec::new(Element::C), AtomSpec::new(Element::H)],
        vec![bond(0, 1)],
        vec![ordinary],
    );
    assert_eq!(selected(&source, &RemoveHsParams::default()), vec![1]);

    let crossing = SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
        .with_atoms(vec![atom(0), atom(1)])
        .with_bonds(vec![bond_id(0)])
        .with_bond_role(bond_id(0), SGroupBondRole::Crossing);
    let source = topology(
        vec![AtomSpec::new(Element::C), AtomSpec::new(Element::H)],
        vec![bond(0, 1)],
        vec![crossing],
    );
    assert!(selected(&source, &RemoveHsParams::default()).is_empty());

    let contained = SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
        .with_atoms(vec![atom(0), atom(1)])
        .with_bonds(vec![bond_id(0)])
        .with_bond_role(bond_id(0), SGroupBondRole::Contained);
    let source = topology(
        vec![AtomSpec::new(Element::C), AtomSpec::new(Element::H)],
        vec![bond(0, 1)],
        vec![contained],
    );
    assert_eq!(selected(&source, &RemoveHsParams::default()), vec![1]);

    let attach = SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
        .with_attach_points(vec![SGroupAttachPoint {
            atom: atom(1),
            leaving_atom: None,
            label: None,
            order: None,
        }]);
    let source = topology(
        vec![AtomSpec::new(Element::C), AtomSpec::new(Element::H)],
        vec![bond(0, 1)],
        vec![attach],
    );
    assert!(selected(&source, &RemoveHsParams::default()).is_empty());

    let leaving = SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
        .with_attach_points(vec![SGroupAttachPoint {
            atom: atom(0),
            leaving_atom: Some(atom(1)),
            label: None,
            order: None,
        }]);
    let source = topology(
        vec![AtomSpec::new(Element::C), AtomSpec::new(Element::H)],
        vec![bond(0, 1)],
        vec![leaving],
    );
    assert_eq!(selected(&source, &RemoveHsParams::default()), vec![1]);

    let cstate = SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
        .with_cstates(vec![SGroupCState {
            bond: bond_id(0),
            vector: [1.0, 0.0, 0.0],
        }]);
    let source = topology(
        vec![AtomSpec::new(Element::C), AtomSpec::new(Element::H)],
        vec![bond(0, 1)],
        vec![cstate],
    );
    assert!(selected(&source, &RemoveHsParams::default()).is_empty());
}

#[test]
fn sgroup_false_uses_complete_includes_atom_membership() {
    let variants =
        vec![
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_atoms(vec![atom(1)]),
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_parent_atoms(vec![atom(1)]),
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_attach_points(vec![SGroupAttachPoint {
                    atom: atom(1),
                    leaving_atom: None,
                    label: None,
                    order: None,
                }]),
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_attach_points(vec![SGroupAttachPoint {
                    atom: atom(0),
                    leaving_atom: Some(atom(1)),
                    label: None,
                    order: None,
                }]),
        ];
    for group in variants {
        let source = topology(
            vec![AtomSpec::new(Element::C), AtomSpec::new(Element::H)],
            vec![bond(0, 1)],
            vec![group],
        );
        assert!(
            selected(
                &source,
                &RemoveHsParams {
                    remove_in_sgroups: false,
                    ..Default::default()
                }
            )
            .is_empty()
        );
    }
}

#[test]
fn sgroup_emptying_filter_covers_atom_parent_combined_and_empty_groups() {
    let atom_only = topology(
        vec![AtomSpec::new(Element::H)],
        Vec::new(),
        vec![
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_atoms(vec![atom(0)]),
        ],
    );
    let degree_zero = RemoveHsParams {
        remove_degree_zero: true,
        ..Default::default()
    };
    assert!(selected(&atom_only, &degree_zero).is_empty());

    let parent_only = topology(
        vec![AtomSpec::new(Element::H)],
        Vec::new(),
        vec![
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_parent_atoms(vec![atom(0)]),
        ],
    );
    assert!(selected(&parent_only, &degree_zero).is_empty());

    let combined = topology(
        vec![AtomSpec::new(Element::H), AtomSpec::new(Element::H)],
        vec![bond(0, 1)],
        vec![
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_atoms(vec![atom(0)])
                .with_parent_atoms(vec![atom(1)]),
        ],
    );
    assert!(
        selected(
            &combined,
            &RemoveHsParams {
                remove_only_h_neighbors: true,
                ..Default::default()
            }
        )
        .is_empty()
    );

    let empty_group = topology(
        vec![AtomSpec::new(Element::C), AtomSpec::new(Element::H)],
        vec![bond(0, 1)],
        vec![SubstanceGroup::new(
            SubstanceGroupId::new(0),
            SubstanceGroupKind::Data,
        )],
    );
    assert_eq!(selected(&empty_group, &RemoveHsParams::default()), vec![1]);
}

#[test]
fn sgroup_emptying_filter_is_ordered_across_overlapping_groups() {
    let source = topology(
        vec![
            AtomSpec::new(Element::H),
            AtomSpec::new(Element::H),
            AtomSpec::new(Element::C),
        ],
        vec![bond(0, 2), bond(1, 2)],
        vec![
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_atoms(vec![atom(0)]),
            SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Data)
                .with_atoms(vec![atom(0), atom(1), atom(2)]),
        ],
    );
    assert_eq!(selected(&source, &RemoveHsParams::default()), vec![1]);
}

#[test]
fn diagnostic_and_explicit_count_flags_do_not_change_candidates() {
    let source = topology(
        vec![AtomSpec::new(Element::H), AtomSpec::new(Element::H)],
        vec![bond(0, 1)],
        Vec::new(),
    );
    let base = RemoveHsParams {
        remove_only_h_neighbors: true,
        ..Default::default()
    };
    let mut quiet = base.clone();
    quiet.show_warnings = false;
    let mut update = base.clone();
    update.update_explicit_count = true;
    assert_eq!(selected(&source, &base), vec![0, 1]);
    assert_eq!(selected(&source, &quiet), selected(&source, &base));
    assert_eq!(selected(&source, &update), selected(&source, &base));
}

#[test]
fn remove_all_source_flag_combination_selects_otherwise_protected_hydrogen() {
    let query_hydride = AtomSpec::new(Element::H)
        .with_isotope(2)
        .with_atom_map(9)
        .with_formal_charge(-1)
        .with_prop("_MolFileAtomQuery", "1")
        .unwrap();
    let source = topology(vec![query_hydride], Vec::new(), Vec::new());
    let remove_all = RemoveHsParams {
        remove_degree_zero: true,
        remove_higher_degrees: true,
        remove_only_h_neighbors: true,
        remove_isotopes: true,
        remove_and_track_isotopes: false,
        remove_dummy_neighbors: true,
        remove_defining_bond_stereo: true,
        remove_with_wedged_bond: true,
        remove_with_query: true,
        remove_mapped: true,
        remove_in_sgroups: true,
        show_warnings: true,
        remove_nonimplicit: true,
        update_explicit_count: false,
        remove_hydrides: true,
        remove_nontetrahedral_neighbors: true,
        sanitize: true,
    };
    assert_eq!(selected(&source, &remove_all), vec![0]);
}
