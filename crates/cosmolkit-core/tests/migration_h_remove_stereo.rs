#![cfg(feature = "op-contracts-strict")]

use cosmolkit_core::{
    __migration_hydrogens::prepare_hydrogen_removal_stereo, HydrogenError, RemoveHsParams,
    ValenceAssignment,
};
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondDirection, BondId, BondSpec, BondStereo, ChiralTag,
    SGroupAttachPoint, SGroupBondRole, SubstanceGroup, SubstanceGroupId, SubstanceGroupKind,
    TopologyBlock,
};
use cosmolkit_types::{BondOrder, Element};

fn atom(index: usize) -> AtomId {
    AtomId::new(index)
}

fn bond_id(index: usize) -> BondId {
    BondId::new(index)
}

fn bond(begin: usize, end: usize, order: BondOrder) -> BondSpec {
    BondSpec::new(atom(begin), atom(end), order)
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

fn valence(explicit: Vec<i32>, implicit: Vec<i32>) -> ValenceAssignment {
    ValenceAssignment {
        explicit_valence: explicit,
        implicit_hydrogens: implicit,
    }
}

fn zero_valence(atom_count: usize) -> ValenceAssignment {
    valence(vec![0; atom_count], vec![0; atom_count])
}

fn prepare(
    source: TopologyBlock,
    candidates: Vec<usize>,
    valence: &ValenceAssignment,
    params: &RemoveHsParams,
) -> Result<cosmolkit_core::__migration_hydrogens::PreparedHydrogenRemoval, HydrogenError> {
    prepare_hydrogen_removal_stereo(
        source,
        candidates.into_iter().map(atom).collect(),
        valence,
        params,
    )
}

#[test]
fn malformed_candidates_and_valence_lengths_are_structured_and_atomic() {
    let source = topology(
        vec![
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::H),
            AtomSpec::new(Element::H),
        ],
        vec![bond(0, 1, BondOrder::Single), bond(0, 2, BondOrder::Single)],
        Vec::new(),
    );
    let snapshot = source.clone();
    let valid_valence = zero_valence(3);
    let cases = [
        (
            vec![3],
            HydrogenError::InvalidRemovalCandidate {
                position: 0,
                atom: atom(3),
                reason: "atom id is out of range",
            },
        ),
        (
            vec![1, 1],
            HydrogenError::InvalidRemovalCandidate {
                position: 1,
                atom: atom(1),
                reason: "candidate id is duplicated",
            },
        ),
        (
            vec![2, 1],
            HydrogenError::InvalidRemovalCandidate {
                position: 1,
                atom: atom(1),
                reason: "candidate ids are not in strict source order",
            },
        ),
        (
            vec![0],
            HydrogenError::InvalidRemovalCandidate {
                position: 0,
                atom: atom(0),
                reason: "candidate is not hydrogen",
            },
        ),
    ];
    for (candidates, expected) in cases {
        assert_eq!(
            prepare(
                source.clone(),
                candidates,
                &valid_valence,
                &RemoveHsParams::default()
            ),
            Err(expected)
        );
        assert_eq!(source, snapshot);
    }

    for (assignment, field, actual) in [
        (valence(vec![0; 2], vec![0; 3]), "explicit_valence", 2),
        (valence(vec![0; 3], vec![0; 2]), "implicit_hydrogens", 2),
    ] {
        assert_eq!(
            prepare(
                source.clone(),
                vec![1],
                &assignment,
                &RemoveHsParams::default()
            ),
            Err(HydrogenError::ValenceAssignmentLength {
                field,
                expected: 3,
                actual,
            })
        );
        assert_eq!(source, snapshot);
    }
}

#[test]
fn explicit_hydrogen_overflow_is_structured_and_atomic() {
    let source = topology(
        vec![
            AtomSpec::new(Element::C)
                .with_explicit_hydrogens(u8::MAX)
                .with_no_implicit(true),
            AtomSpec::new(Element::H),
        ],
        vec![bond(0, 1, BondOrder::Single)],
        Vec::new(),
    );
    let snapshot = source.clone();
    assert_eq!(
        prepare(
            source.clone(),
            vec![1],
            &zero_valence(2),
            &RemoveHsParams::default()
        ),
        Err(HydrogenError::ExplicitHydrogenOverflow {
            atom: atom(0),
            current: u8::MAX,
        })
    );
    assert_eq!(source, snapshot);
}

#[test]
fn isotope_tracking_disabled_preserves_and_enabled_rebuilds_in_bond_order() {
    let specs = vec![
        AtomSpec::new(Element::C).with_tracked_isotopic_hydrogens(vec![9]),
        AtomSpec::new(Element::H)
            .with_isotope(3)
            .with_tracked_isotopic_hydrogens(vec![8]),
        AtomSpec::new(Element::H).with_isotope(2),
        AtomSpec::new(Element::O).with_tracked_isotopic_hydrogens(vec![7]),
        AtomSpec::new(Element::H).with_isotope(0),
        AtomSpec::new(Element::H).with_isotope(2),
        AtomSpec::new(Element::H).with_isotope(3),
        AtomSpec::new(Element::H).with_isotope(2),
        AtomSpec::new(Element::H).with_isotope(4),
    ];
    let bonds = vec![
        bond(0, 2, BondOrder::Single),
        bond(1, 0, BondOrder::Single),
        bond(3, 4, BondOrder::Single),
        bond(5, 6, BondOrder::Single),
        bond(0, 7, BondOrder::Single),
        bond(8, 3, BondOrder::Single),
    ];
    let source = topology(specs, bonds, Vec::new());
    let preserved = prepare(
        source.clone(),
        Vec::new(),
        &zero_valence(9),
        &RemoveHsParams::default(),
    )
    .unwrap();
    assert_eq!(preserved.topology, source);

    let rebuilt = prepare(
        source,
        Vec::new(),
        &zero_valence(9),
        &RemoveHsParams {
            remove_and_track_isotopes: true,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(
        rebuilt.topology.atoms[0].tracked_isotopic_hydrogens(),
        &[2, 3, 2]
    );
    assert!(rebuilt.topology.atoms[3].tracked_isotopic_hydrogens() == [4]);
    for index in [1, 2, 4, 5, 6, 7, 8] {
        assert!(
            rebuilt.topology.atoms[index]
                .tracked_isotopic_hydrogens()
                .is_empty()
        );
    }
}

#[test]
fn direct_explicit_count_rules_cover_option_no_implicit_and_all_chiral_tags() {
    let direct_specs = [
        AtomSpec::new(Element::C),
        AtomSpec::new(Element::C).with_no_implicit(true),
    ];
    for (spec, update_explicit_count) in [
        (direct_specs[0].clone(), true),
        (direct_specs[1].clone(), false),
    ] {
        let source = topology(
            vec![spec, AtomSpec::new(Element::H)],
            vec![bond(0, 1, BondOrder::Single)],
            Vec::new(),
        );
        let result = prepare(
            source,
            vec![1],
            &zero_valence(2),
            &RemoveHsParams {
                update_explicit_count,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(result.topology.atoms[0].explicit_hydrogens(), 1);
    }

    for tag in [
        ChiralTag::TetrahedralCw,
        ChiralTag::TetrahedralCcw,
        ChiralTag::Other,
        ChiralTag::Tetrahedral,
        ChiralTag::Allene,
        ChiralTag::SquarePlanar,
        ChiralTag::TrigonalBipyramidal,
        ChiralTag::Octahedral,
    ] {
        let source = topology(
            vec![
                AtomSpec::new(Element::C).with_chiral_tag(tag),
                AtomSpec::new(Element::H),
            ],
            vec![bond(0, 1, BondOrder::Single)],
            Vec::new(),
        );
        let result = prepare(
            source,
            vec![1],
            &zero_valence(2),
            &RemoveHsParams::default(),
        )
        .unwrap();
        assert_eq!(result.topology.atoms[0].explicit_hydrogens(), 1, "{tag:?}");
    }
}

#[test]
fn aromatic_and_nondefault_valence_count_rules_cover_acceptance_and_rejection() {
    for element in [Element::N, Element::P] {
        let source = topology(
            vec![
                AtomSpec::new(element).with_aromatic(true),
                AtomSpec::new(Element::H),
            ],
            vec![bond(0, 1, BondOrder::Single)],
            Vec::new(),
        );
        let result = prepare(
            source,
            vec![1],
            &zero_valence(2),
            &RemoveHsParams::default(),
        )
        .unwrap();
        assert_eq!(result.topology.atoms[0].explicit_hydrogens(), 1);
    }

    let aromatic_carbon = topology(
        vec![
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::H),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
        ],
        vec![
            bond(0, 1, BondOrder::Single),
            bond(0, 2, BondOrder::Aromatic).with_aromatic(true),
            bond(0, 3, BondOrder::Aromatic).with_aromatic(true),
        ],
        Vec::new(),
    );
    let result = prepare(
        aromatic_carbon,
        vec![1],
        &valence(vec![3, 0, 0, 0], vec![0; 4]),
        &RemoveHsParams::default(),
    )
    .unwrap();
    assert_eq!(result.topology.atoms[0].explicit_hydrogens(), 1);

    let sulfur = topology(
        vec![AtomSpec::new(Element::S), AtomSpec::new(Element::H)],
        vec![bond(0, 1, BondOrder::Single)],
        Vec::new(),
    );
    let result = prepare(
        sulfur,
        vec![1],
        &valence(vec![4, 0], vec![0; 2]),
        &RemoveHsParams::default(),
    )
    .unwrap();
    assert_eq!(result.topology.atoms[0].explicit_hydrogens(), 1);

    let rejected = topology(
        vec![
            AtomSpec::new(Element::C).with_aromatic(true),
            AtomSpec::new(Element::H),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
        ],
        vec![
            bond(0, 1, BondOrder::Single),
            bond(0, 2, BondOrder::Double),
            bond(0, 3, BondOrder::Aromatic).with_aromatic(true),
        ],
        Vec::new(),
    );
    let result = prepare(
        rejected,
        vec![1],
        &valence(vec![3, 0, 0, 0], vec![0; 4]),
        &RemoveHsParams::default(),
    )
    .unwrap();
    assert_eq!(result.topology.atoms[0].explicit_hydrogens(), 0);
}

fn chiral_center(
    tag: ChiralTag,
    permutation: Option<u32>,
    hydrogen_bond_position: usize,
) -> TopologyBlock {
    let mut center = AtomSpec::new(Element::C).with_chiral_tag(tag);
    if let Some(permutation) = permutation {
        center = center.with_chiral_permutation(permutation);
    }
    let endpoints = [(0, 1), (0, 2), (0, 3), (0, 4)];
    let mut specs = Vec::new();
    for (begin, end) in endpoints {
        specs.push(bond(begin, end, BondOrder::Single));
    }
    let hydrogen_atom = hydrogen_bond_position + 1;
    let mut atoms = vec![center];
    atoms.extend((1..=4).map(|index| {
        if index == hydrogen_atom {
            AtomSpec::new(Element::H)
        } else {
            AtomSpec::new(Element::C)
        }
    }));
    topology(atoms, specs, Vec::new())
}

#[test]
fn perturbation_parity_inverts_only_odd_bond_order() {
    let odd = chiral_center(ChiralTag::TetrahedralCw, None, 0);
    let odd_result = prepare(odd, vec![1], &zero_valence(5), &RemoveHsParams::default()).unwrap();
    assert_eq!(
        odd_result.topology.atoms[0].chiral_tag(),
        ChiralTag::TetrahedralCcw
    );

    // H atom 2 has the second bond-table position, so moving it last takes two swaps.
    let even = chiral_center(ChiralTag::TetrahedralCw, None, 1);
    let even_result = prepare(even, vec![2], &zero_valence(5), &RemoveHsParams::default()).unwrap();
    assert_eq!(
        even_result.topology.atoms[0].chiral_tag(),
        ChiralTag::TetrahedralCw
    );
}

#[test]
fn typed_chiral_inversion_tables_and_source_default_cases_are_exact() {
    for (tag, input, expected) in [
        (ChiralTag::Tetrahedral, 1, 2),
        (ChiralTag::Tetrahedral, 9, 0),
        (ChiralTag::TrigonalBipyramidal, 9, 11),
        (ChiralTag::Octahedral, 3, 16),
        (ChiralTag::SquarePlanar, 2, 2),
        (ChiralTag::Allene, 2, 2),
    ] {
        let source = chiral_center(tag, Some(input), 0);
        let result = prepare(
            source,
            vec![1],
            &zero_valence(5),
            &RemoveHsParams::default(),
        )
        .unwrap();
        assert_eq!(
            result.topology.atoms[0].chiral_permutation(),
            Some(expected),
            "{tag:?}"
        );
    }
}

#[test]
fn multiple_hydrogens_are_processed_descending_against_active_graph_state() {
    let source = topology(
        vec![
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::H),
            AtomSpec::new(Element::H),
        ],
        vec![bond(0, 1, BondOrder::Single), bond(0, 2, BondOrder::Single)],
        Vec::new(),
    );
    let result = prepare(
        source,
        vec![1, 2],
        &zero_valence(3),
        &RemoveHsParams {
            update_explicit_count: true,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(result.atoms_to_remove, vec![atom(1), atom(2)]);
    assert_eq!(result.topology.atoms[0].explicit_hydrogens(), 2);
    assert_eq!(result.topology.atoms.len(), 3);
    assert_eq!(result.topology.bonds.len(), 2);
}

#[test]
fn degree_two_defining_stereo_is_cleared_but_none_and_any_are_preserved() {
    for stereo in [
        BondStereo::E,
        BondStereo::Z,
        BondStereo::Cis,
        BondStereo::Trans,
    ] {
        let source = topology(
            vec![
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::H),
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::C),
            ],
            vec![
                bond(0, 1, BondOrder::Single),
                bond(0, 2, BondOrder::Double)
                    .with_stereo_atoms(atom(1), atom(3))
                    .with_stereo(stereo),
                bond(2, 3, BondOrder::Single),
            ],
            Vec::new(),
        );
        let result = prepare(
            source,
            vec![1],
            &zero_valence(4),
            &RemoveHsParams::default(),
        )
        .unwrap();
        assert_eq!(result.topology.bonds[1].stereo(), BondStereo::None);
        assert_eq!(result.topology.bonds[1].stereo_atoms(), None);
    }

    for stereo in [BondStereo::None, BondStereo::Any] {
        let source = topology(
            vec![
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::H),
                AtomSpec::new(Element::C),
            ],
            vec![
                bond(0, 1, BondOrder::Single),
                bond(0, 2, BondOrder::Double).with_stereo(stereo),
            ],
            Vec::new(),
        );
        let result = prepare(
            source,
            vec![1],
            &zero_valence(3),
            &RemoveHsParams::default(),
        )
        .unwrap();
        assert_eq!(result.topology.bonds[1].stereo(), stereo);
    }
}

#[test]
fn degree_three_stereo_reference_replacement_swaps_only_cis_trans() {
    for (stereo, expected) in [
        (BondStereo::Cis, BondStereo::Trans),
        (BondStereo::Trans, BondStereo::Cis),
        (BondStereo::E, BondStereo::E),
        (BondStereo::Z, BondStereo::Z),
    ] {
        let source = topology(
            vec![
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::H),
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::C),
            ],
            vec![
                bond(0, 1, BondOrder::Single),
                bond(0, 2, BondOrder::Double)
                    .with_stereo_atoms(atom(1), atom(4))
                    .with_stereo(stereo),
                bond(0, 3, BondOrder::Single),
                bond(2, 4, BondOrder::Single),
            ],
            Vec::new(),
        );
        let result = prepare(
            source,
            vec![1],
            &zero_valence(5),
            &RemoveHsParams::default(),
        )
        .unwrap();
        assert_eq!(
            result.topology.bonds[1].stereo_atoms(),
            Some([atom(3), atom(4)])
        );
        assert_eq!(result.topology.bonds[1].stereo(), expected);
    }

    let no_reference = topology(
        vec![
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::H),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
        ],
        vec![
            bond(0, 1, BondOrder::Single),
            bond(0, 2, BondOrder::Double).with_stereo(BondStereo::E),
            bond(0, 3, BondOrder::Single),
        ],
        Vec::new(),
    );
    let result = prepare(
        no_reference,
        vec![1],
        &zero_valence(4),
        &RemoveHsParams::default(),
    )
    .unwrap();
    assert_eq!(result.topology.bonds[1].stereo(), BondStereo::E);
    assert_eq!(result.topology.bonds[1].stereo_atoms(), None);
}

#[test]
fn unknown_direction_is_begin_only_and_direction_propagation_is_exact() {
    for (begin, end, expected) in [(0, 1, true), (1, 0, false)] {
        let source = topology(
            vec![AtomSpec::new(Element::C), AtomSpec::new(Element::H)],
            vec![bond(begin, end, BondOrder::Single).with_direction(BondDirection::Unknown)],
            Vec::new(),
        );
        let result = prepare(
            source,
            vec![1],
            &zero_valence(2),
            &RemoveHsParams::default(),
        )
        .unwrap();
        assert_eq!(result.topology.atoms[0].unknown_stereo(), expected);
    }

    let source = topology(
        vec![
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::H),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
        ],
        vec![
            bond(0, 1, BondOrder::Single).with_direction(BondDirection::EndDownRight),
            bond(0, 2, BondOrder::Single),
            bond(0, 3, BondOrder::Single),
        ],
        Vec::new(),
    );
    let result = prepare(
        source,
        vec![1],
        &zero_valence(4),
        &RemoveHsParams::default(),
    )
    .unwrap();
    assert_eq!(result.topology.bonds[1].direction(), BondDirection::None);
    assert_eq!(
        result.topology.bonds[2].direction(),
        BondDirection::EndUpRight
    );

    let copied = topology(
        vec![
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::H),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
        ],
        vec![
            bond(1, 0, BondOrder::Single).with_direction(BondDirection::EndUpRight),
            bond(0, 2, BondOrder::Single),
            bond(0, 3, BondOrder::Double),
        ],
        Vec::new(),
    );
    let result = prepare(
        copied,
        vec![1],
        &zero_valence(4),
        &RemoveHsParams::default(),
    )
    .unwrap();
    assert_eq!(
        result.topology.bonds[1].direction(),
        BondDirection::EndUpRight
    );
    assert_eq!(result.topology.bonds[2].direction(), BondDirection::None);

    let suppressed = topology(
        vec![
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::H),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
        ],
        vec![
            bond(1, 0, BondOrder::Single).with_direction(BondDirection::EndUpRight),
            bond(0, 2, BondOrder::Single),
            bond(0, 3, BondOrder::Single).with_direction(BondDirection::BeginWedge),
            bond(0, 4, BondOrder::Double),
        ],
        Vec::new(),
    );
    let result = prepare(
        suppressed,
        vec![1],
        &zero_valence(5),
        &RemoveHsParams::default(),
    )
    .unwrap();
    assert_eq!(result.topology.bonds[1].direction(), BondDirection::None);
    assert_eq!(result.topology.bonds[3].direction(), BondDirection::None);
}

#[test]
fn sgroup_references_are_cleaned_one_occurrence_without_metadata_loss() {
    let group = SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
        .with_atoms(vec![atom(1), atom(1), atom(0)])
        .with_bonds(vec![bond_id(0), bond_id(0)])
        .with_bond_role(bond_id(0), SGroupBondRole::Contained)
        .with_parent_atoms(vec![atom(1), atom(1)])
        .with_attach_points(vec![
            SGroupAttachPoint {
                atom: atom(0),
                leaving_atom: Some(atom(1)),
                label: Some("leave".into()),
                order: Some(7),
            },
            SGroupAttachPoint {
                atom: atom(1),
                leaving_atom: None,
                label: Some("primary".into()),
                order: Some(8),
            },
        ])
        .with_label("kept")
        .with_prop("key", "value");
    let source = topology(
        vec![AtomSpec::new(Element::C), AtomSpec::new(Element::H)],
        vec![bond(0, 1, BondOrder::Single)],
        vec![group],
    );
    let result = prepare(
        source,
        vec![1],
        &zero_valence(2),
        &RemoveHsParams::default(),
    )
    .unwrap();
    let group = &result.topology.substance_groups[0];
    assert_eq!(group.atoms(), &[atom(1), atom(0)]);
    assert_eq!(group.parent_atoms(), &[atom(1)]);
    assert_eq!(group.bonds(), &[bond_id(0)]);
    assert_eq!(group.bond_role(bond_id(0)), SGroupBondRole::Contained);
    assert_eq!(group.attach_points()[0].atom, atom(0));
    assert_eq!(group.attach_points()[0].leaving_atom, None);
    assert_eq!(group.attach_points()[0].label.as_deref(), Some("leave"));
    assert_eq!(group.attach_points()[1].atom, atom(1));
    assert_eq!(group.label(), Some("kept"));
    assert_eq!(group.props().get("key").map(String::as_str), Some("value"));
}

#[test]
fn higher_degree_h_updates_every_neighbor_and_preserves_original_id_space() {
    let source = topology(
        vec![
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::H),
            AtomSpec::new(Element::N),
        ],
        vec![bond(0, 1, BondOrder::Single), bond(1, 2, BondOrder::Single)],
        Vec::new(),
    );
    let snapshot = source.clone();
    let result = prepare(
        source,
        vec![1],
        &zero_valence(3),
        &RemoveHsParams {
            update_explicit_count: true,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(result.atoms_to_remove, vec![atom(1)]);
    assert_eq!(result.topology.atoms.len(), snapshot.atoms.len());
    assert_eq!(result.topology.bonds.len(), snapshot.bonds.len());
    assert_eq!(result.topology.adjacency, snapshot.adjacency);
    assert_eq!(result.topology.atoms[0].explicit_hydrogens(), 1);
    assert_eq!(result.topology.atoms[2].explicit_hydrogens(), 1);
    assert_eq!(snapshot.atoms[0].explicit_hydrogens(), 0);
    assert_eq!(snapshot.atoms[2].explicit_hydrogens(), 0);
    result.topology.validate().unwrap();
}
