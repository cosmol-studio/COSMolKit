#![cfg(feature = "op-contracts-strict")]

use cosmolkit_core::{
    __migration_sanitize::{
        AdjustHsError, HybridizationAssignment, adjust_hs, cleanup_chirality,
        cleanup_invalid_atropisomers,
    },
    AtropisomerError, ChemistryProblemError, KekulizeError, KekulizeParams, RingFindType, RingInfo,
    RingSearchParams, SanitizeError, SanitizeOperations, SanitizeParams, SanitizeStage,
    StereoError, ValenceAssignment, ValenceError, ValenceModel, ValencePhase,
    assign_valence_with_options_for_topology, detect_chemistry_problems, find_sssr, kekulize,
    sanitize_topology,
};
use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, StereoGroup, StereoGroupKind,
    TopologyBlock,
};
use cosmolkit_types::{BondOrder, BondStereo, ChiralTag, Element, Hybridization};

fn carbon(id: usize) -> Atom {
    Atom::from_spec(AtomId::new(id), AtomSpec::new(Element::C))
}

fn bond(id: usize, begin: usize, end: usize, stereo: BondStereo) -> Bond {
    Bond::from_spec(
        BondId::new(id),
        BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single).with_stereo(stereo),
    )
}

fn topology(atom_count: usize, bonds: Vec<Bond>, groups: Vec<StereoGroup>) -> TopologyBlock {
    TopologyBlock::try_from_parts(
        (0..atom_count).map(carbon).collect(),
        bonds,
        Vec::new(),
        groups,
    )
    .unwrap()
}

fn topology_from_specs(atom_specs: Vec<AtomSpec>, bond_specs: Vec<BondSpec>) -> TopologyBlock {
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
    TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
}

fn atom_spec(atomic_number: u8) -> AtomSpec {
    AtomSpec::new(Element::from_atomic_number(atomic_number).unwrap())
}

fn bond_spec(begin: usize, end: usize, order: BondOrder) -> BondSpec {
    BondSpec::new(AtomId::new(begin), AtomId::new(end), order)
}

fn alternating_cycle(atom_specs: Vec<AtomSpec>) -> TopologyBlock {
    let count = atom_specs.len();
    topology_from_specs(
        atom_specs,
        (0..count)
            .map(|index| {
                bond_spec(
                    index,
                    (index + 1) % count,
                    if index % 2 == 0 {
                        BondOrder::Double
                    } else {
                        BondOrder::Single
                    },
                )
            })
            .collect(),
    )
}

fn aromatic_cycle(size: usize) -> TopologyBlock {
    topology_from_specs(
        (0..size)
            .map(|_| AtomSpec::new(Element::C).with_aromatic(true))
            .collect(),
        (0..size)
            .map(|index| {
                BondSpec::new(
                    AtomId::new(index),
                    AtomId::new((index + 1) % size),
                    BondOrder::Aromatic,
                )
                .with_aromatic(true)
            })
            .collect(),
    )
}

fn hybridization(
    atom_count: usize,
    overrides: &[(usize, Hybridization)],
) -> HybridizationAssignment {
    let mut values = vec![Hybridization::Sp2; atom_count];
    for &(atom, value) in overrides {
        values[atom] = value;
    }
    HybridizationAssignment { values }
}

fn rings(topology: &TopologyBlock) -> RingInfo {
    find_sssr(topology, &RingSearchParams::default()).unwrap()
}

fn cycle(size: usize, stereo: BondStereo) -> TopologyBlock {
    let bonds = (0..size)
        .map(|id| {
            bond(
                id,
                id,
                (id + 1) % size,
                if id == 0 { stereo } else { BondStereo::None },
            )
        })
        .collect();
    topology(size, bonds, Vec::new())
}

fn chiral_atom(
    id: usize,
    tag: ChiralTag,
    hybridization: Hybridization,
    permutation: Option<u32>,
) -> Atom {
    let mut spec = AtomSpec::new(Element::C)
        .with_chiral_tag(tag)
        .with_hybridization(hybridization);
    if let Some(permutation) = permutation {
        spec = spec
            .with_chiral_permutation(permutation)
            .with_prop("_chiralPermutation", permutation.to_string())
            .unwrap();
    }
    Atom::from_spec(AtomId::new(id), spec)
}

fn star(
    tag: ChiralTag,
    hybridization: Hybridization,
    permutation: Option<u32>,
    neighbor_count: usize,
    implicit_hydrogens: i32,
) -> (TopologyBlock, ValenceAssignment) {
    let mut atoms = vec![chiral_atom(0, tag, hybridization, permutation)];
    atoms.extend((1..=neighbor_count).map(carbon));
    let bonds = (0..neighbor_count)
        .map(|id| bond(id, 0, id + 1, BondStereo::None))
        .collect::<Vec<_>>();
    let topology = TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap();
    let mut implicit = vec![0; topology.atoms.len()];
    implicit[0] = implicit_hydrogens;
    (
        topology,
        ValenceAssignment {
            explicit_valence: vec![0; neighbor_count + 1],
            implicit_hydrogens: implicit,
        },
    )
}

#[test]
fn invalid_atropisomer_cleanup_preserves_acyclic_sp2_bonds_in_both_directions() {
    for stereo in [BondStereo::AtropCw, BondStereo::AtropCcw] {
        let input = topology(2, vec![bond(0, 0, 1, stereo)], Vec::new());
        let output =
            cleanup_invalid_atropisomers(&input, &hybridization(2, &[]), &rings(&input)).unwrap();
        assert_eq!(output, input);
    }
}

#[test]
fn invalid_atropisomer_cleanup_clears_each_non_sp2_endpoint_and_both_directions() {
    for stereo in [BondStereo::AtropCw, BondStereo::AtropCcw] {
        for endpoint in [0, 1] {
            let input = topology(2, vec![bond(0, 0, 1, stereo)], Vec::new());
            let output = cleanup_invalid_atropisomers(
                &input,
                &hybridization(2, &[(endpoint, Hybridization::Sp3)]),
                &rings(&input),
            )
            .unwrap();
            assert_eq!(output.bonds[0].stereo(), BondStereo::None);
            assert_eq!(input.bonds[0].stereo(), stereo);
        }
    }
}

#[test]
fn invalid_atropisomer_cleanup_clears_small_ring_but_preserves_eight_member_macrocycle() {
    let seven = cycle(7, BondStereo::AtropCw);
    let seven_output =
        cleanup_invalid_atropisomers(&seven, &hybridization(7, &[]), &rings(&seven)).unwrap();
    assert_eq!(seven_output.bonds[0].stereo(), BondStereo::None);

    let eight = cycle(8, BondStereo::AtropCcw);
    let eight_output =
        cleanup_invalid_atropisomers(&eight, &hybridization(8, &[]), &rings(&eight)).unwrap();
    assert_eq!(eight_output.bonds[0].stereo(), BondStereo::AtropCcw);
}

#[test]
fn invalid_atropisomer_cleanup_updates_only_affected_group_content_and_preserves_order() {
    let affected = StereoGroup::new(
        StereoGroupKind::Or,
        vec![AtomId::new(0), AtomId::new(2)],
        Vec::new(),
    )
    .with_id(17);
    let unaffected =
        StereoGroup::new(StereoGroupKind::And, vec![AtomId::new(1)], Vec::new()).with_id(19);
    let input = topology(
        4,
        vec![
            bond(0, 0, 1, BondStereo::AtropCw),
            bond(1, 2, 3, BondStereo::AtropCcw),
        ],
        vec![affected, unaffected.clone()],
    );
    let snapshot = input.clone();
    let output = cleanup_invalid_atropisomers(
        &input,
        &hybridization(4, &[(0, Hybridization::Sp3)]),
        &rings(&input),
    )
    .unwrap();

    assert_eq!(input, snapshot);
    assert_eq!(output.atoms, input.atoms);
    assert_eq!(output.adjacency, input.adjacency);
    assert_eq!(
        output.bonds.iter().map(Bond::id).collect::<Vec<_>>(),
        vec![BondId::new(0), BondId::new(1)]
    );
    assert_eq!(output.bonds[0].stereo(), BondStereo::None);
    assert_eq!(output.bonds[1].stereo(), BondStereo::AtropCcw);
    assert_eq!(output.stereo_groups[0].id(), Some(17));
    assert_eq!(output.stereo_groups[0].atoms(), &[AtomId::new(0)]);
    assert_eq!(output.stereo_groups[0].bonds(), &[BondId::new(1)]);
    assert_eq!(output.stereo_groups[1], unaffected);
}

#[test]
fn invalid_atropisomer_cleanup_rejects_mismatched_assignments_before_mutation() {
    let input = topology(2, vec![bond(0, 0, 1, BondStereo::AtropCw)], Vec::new());
    assert_eq!(
        cleanup_invalid_atropisomers(
            &input,
            &HybridizationAssignment {
                values: vec![Hybridization::Sp2],
            },
            &rings(&input),
        ),
        Err(AtropisomerError::HybridizationAssignmentLength {
            actual: 1,
            expected: 2,
        })
    );
    assert_eq!(input.bonds[0].stereo(), BondStereo::AtropCw);

    assert_eq!(
        cleanup_invalid_atropisomers(
            &input,
            &hybridization(2, &[]),
            &RingInfo::new(RingFindType::Fast, 2, 1),
        ),
        Err(AtropisomerError::RingInfoNotSssr)
    );
    assert_eq!(
        cleanup_invalid_atropisomers(
            &input,
            &hybridization(2, &[]),
            &RingInfo::new(RingFindType::Sssr, 1, 1),
        ),
        Err(AtropisomerError::RingAtomRowCount {
            actual: 1,
            expected: 2,
        })
    );
    assert_eq!(
        cleanup_invalid_atropisomers(
            &input,
            &hybridization(2, &[]),
            &RingInfo::new(RingFindType::Sssr, 2, 0),
        ),
        Err(AtropisomerError::RingBondRowCount {
            actual: 0,
            expected: 1,
        })
    );

    let invalid = TopologyBlock {
        adjacency: AdjacencyList::from_topology(0, &[]),
        ..input.clone()
    };
    assert!(matches!(
        cleanup_invalid_atropisomers(&invalid, &hybridization(2, &[]), &rings(&input),),
        Err(AtropisomerError::InvalidTopology { .. })
    ));
}

#[test]
fn chirality_cleanup_applies_tetrahedral_hybridization_and_permutation_rules() {
    for tag in [
        ChiralTag::TetrahedralCw,
        ChiralTag::TetrahedralCcw,
        ChiralTag::Tetrahedral,
    ] {
        let (valid, valence) = star(tag, Hybridization::Sp3, Some(2), 1, 0);
        let output = cleanup_chirality(&valid, &valence).unwrap();
        assert_eq!(output.atoms[0].chiral_tag(), tag);
        assert_eq!(output.atoms[0].chiral_permutation(), Some(2));

        let (invalid, valence) = star(tag, Hybridization::Sp2, Some(2), 1, 0);
        let snapshot = invalid.clone();
        let output = cleanup_chirality(&invalid, &valence).unwrap();
        assert_eq!(output.atoms[0].chiral_tag(), ChiralTag::Unspecified);
        assert_eq!(output.atoms[0].chiral_permutation(), Some(2));
        assert_eq!(invalid, snapshot);
    }

    let (too_large, valence) = star(ChiralTag::Tetrahedral, Hybridization::Sp3, Some(3), 1, 0);
    let output = cleanup_chirality(&too_large, &valence).unwrap();
    assert_eq!(output.atoms[0].chiral_tag(), ChiralTag::Tetrahedral);
    assert_eq!(output.atoms[0].chiral_permutation(), Some(0));
    assert_eq!(output.atoms[0].prop("_chiralPermutation"), Some("0"));
}

#[test]
fn chirality_cleanup_enforces_non_tetrahedral_degree_and_exact_permutation_limits() {
    for (tag, maximum_degree, maximum_permutation) in [
        (ChiralTag::SquarePlanar, 4, 3),
        (ChiralTag::TrigonalBipyramidal, 5, 20),
        (ChiralTag::Octahedral, 6, 30),
    ] {
        let (at_lower_degree, valence) = star(
            tag,
            Hybridization::Unspecified,
            Some(maximum_permutation),
            1,
            1,
        );
        let output = cleanup_chirality(&at_lower_degree, &valence).unwrap();
        assert_eq!(output.atoms[0].chiral_tag(), tag);
        assert_eq!(
            output.atoms[0].chiral_permutation(),
            Some(maximum_permutation)
        );
        assert_eq!(
            output.atoms[0].prop("_chiralPermutation"),
            Some(maximum_permutation.to_string().as_str())
        );

        let (too_small, valence) = star(
            tag,
            Hybridization::Unspecified,
            Some(maximum_permutation),
            1,
            0,
        );
        let output = cleanup_chirality(&too_small, &valence).unwrap();
        assert_eq!(output.atoms[0].chiral_tag(), ChiralTag::Unspecified);
        assert_eq!(
            output.atoms[0].chiral_permutation(),
            Some(maximum_permutation)
        );

        let (too_large_degree, valence) = star(
            tag,
            Hybridization::Unspecified,
            Some(maximum_permutation),
            maximum_degree + 1,
            0,
        );
        let output = cleanup_chirality(&too_large_degree, &valence).unwrap();
        assert_eq!(output.atoms[0].chiral_tag(), ChiralTag::Unspecified);
        assert_eq!(
            output.atoms[0].chiral_permutation(),
            Some(maximum_permutation)
        );

        let (too_large_permutation, valence) = star(
            tag,
            Hybridization::Unspecified,
            Some(maximum_permutation + 1),
            2,
            0,
        );
        let output = cleanup_chirality(&too_large_permutation, &valence).unwrap();
        assert_eq!(output.atoms[0].chiral_tag(), tag);
        assert_eq!(output.atoms[0].chiral_permutation(), Some(0));
        assert_eq!(output.atoms[0].prop("_chiralPermutation"), Some("0"));
    }
}

#[test]
fn chirality_cleanup_selectively_rebuilds_enhanced_stereo_groups_in_source_order() {
    let affected = StereoGroup::new(
        StereoGroupKind::Or,
        vec![AtomId::new(0), AtomId::new(1)],
        vec![BondId::new(0), BondId::new(1)],
    )
    .with_id(7);
    let dropped = StereoGroup::new(
        StereoGroupKind::And,
        vec![AtomId::new(0)],
        vec![BondId::new(1)],
    )
    .with_id(8);
    let unaffected = StereoGroup::new(
        StereoGroupKind::Absolute,
        vec![AtomId::new(1)],
        vec![BondId::new(0)],
    )
    .with_id(9);
    let atoms = vec![
        chiral_atom(0, ChiralTag::TetrahedralCw, Hybridization::Sp2, None),
        chiral_atom(1, ChiralTag::TetrahedralCcw, Hybridization::Sp3, None),
        carbon(2),
    ];
    let bonds = vec![
        bond(0, 1, 2, BondStereo::AtropCw),
        bond(1, 0, 2, BondStereo::None),
    ];
    let input = TopologyBlock::try_from_parts(
        atoms,
        bonds,
        Vec::new(),
        vec![affected, dropped, unaffected.clone()],
    )
    .unwrap();
    let snapshot = input.clone();
    let valence = ValenceAssignment {
        explicit_valence: vec![0; 3],
        implicit_hydrogens: vec![0; 3],
    };
    let output = cleanup_chirality(&input, &valence).unwrap();

    assert_eq!(input, snapshot);
    assert_eq!(
        output.atoms.iter().map(Atom::id).collect::<Vec<_>>(),
        vec![AtomId::new(0), AtomId::new(1), AtomId::new(2)]
    );
    assert_eq!(
        output.bonds.iter().map(Bond::id).collect::<Vec<_>>(),
        vec![BondId::new(0), BondId::new(1)]
    );
    assert_eq!(output.stereo_groups.len(), 2);
    assert_eq!(output.stereo_groups[0].id(), Some(7));
    assert_eq!(output.stereo_groups[0].atoms(), &[AtomId::new(1)]);
    assert_eq!(output.stereo_groups[0].bonds(), &[BondId::new(0)]);
    assert_eq!(output.stereo_groups[1], unaffected);
}

#[test]
fn chirality_cleanup_rejects_invalid_topology_and_valence_rows() {
    let (input, valence) = star(
        ChiralTag::SquarePlanar,
        Hybridization::Unspecified,
        Some(3),
        2,
        0,
    );
    assert_eq!(
        cleanup_chirality(
            &input,
            &ValenceAssignment {
                explicit_valence: vec![0],
                implicit_hydrogens: vec![0; 3],
            },
        ),
        Err(StereoError::InvalidValence {
            field: "explicit_valence",
            actual: 1,
            atom_count: 3,
        })
    );
    assert_eq!(
        cleanup_chirality(
            &input,
            &ValenceAssignment {
                explicit_valence: vec![0; 3],
                implicit_hydrogens: vec![-1, 0, 0],
            },
        ),
        Err(StereoError::InvalidValenceValue {
            field: "implicit_hydrogens",
            atom: AtomId::new(0),
            value: -1,
        })
    );

    let invalid = TopologyBlock {
        adjacency: AdjacencyList::from_topology(0, &[]),
        ..input
    };
    assert!(matches!(
        cleanup_chirality(&invalid, &valence),
        Err(StereoError::InvalidTopology(_))
    ));
}

#[test]
fn adjust_hs_transfers_only_implicit_decreases_and_accumulates_explicit_hydrogens_in_atom_order() {
    let atoms = vec![
        Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C).with_explicit_hydrogens(1),
        ),
        carbon(1),
        carbon(2),
    ];
    let input = TopologyBlock::try_from_parts(atoms, Vec::new(), Vec::new(), Vec::new()).unwrap();
    let snapshot = input.clone();
    let current =
        assign_valence_with_options_for_topology(&input, ValenceModel::RdkitLike, false).unwrap();
    assert_eq!(current.implicit_hydrogens, vec![3, 4, 4]);

    let original = ValenceAssignment {
        explicit_valence: current.explicit_valence.clone(),
        // decrease, equality, and increase relative to the recalculated rows
        implicit_hydrogens: vec![5, 4, 3],
    };
    let output = adjust_hs(&input, &original).unwrap();

    assert_eq!(input, snapshot);
    assert_eq!(
        output
            .topology
            .atoms
            .iter()
            .map(Atom::explicit_hydrogens)
            .collect::<Vec<_>>(),
        vec![3, 0, 0]
    );
    assert_eq!(output.valence.implicit_hydrogens, vec![3, 4, 4]);
    assert_eq!(output.valence.explicit_valence, vec![3, 0, 0]);
    assert_eq!(
        output
            .topology
            .atoms
            .iter()
            .map(Atom::id)
            .collect::<Vec<_>>(),
        vec![AtomId::new(0), AtomId::new(1), AtomId::new(2)]
    );
}

#[test]
fn adjust_hs_preserves_aromatic_post_state_and_all_non_atom_references() {
    let atoms = vec![
        Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::N).with_aromatic(true),
        ),
        Atom::from_spec(
            AtomId::new(1),
            AtomSpec::new(Element::C).with_aromatic(true),
        ),
        Atom::from_spec(
            AtomId::new(2),
            AtomSpec::new(Element::C).with_aromatic(true),
        ),
    ];
    let bonds = (0..3)
        .map(|id| {
            Bond::from_spec(
                BondId::new(id),
                BondSpec::new(
                    AtomId::new(id),
                    AtomId::new((id + 1) % 3),
                    BondOrder::Aromatic,
                )
                .with_aromatic(true),
            )
        })
        .collect::<Vec<_>>();
    let group = StereoGroup::new(
        StereoGroupKind::Or,
        vec![AtomId::new(0)],
        vec![BondId::new(0)],
    )
    .with_id(23);
    let input = TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), vec![group]).unwrap();
    let snapshot = input.clone();
    let mut original =
        assign_valence_with_options_for_topology(&input, ValenceModel::RdkitLike, false).unwrap();
    assert_eq!(original.implicit_hydrogens[0], 0);
    original.implicit_hydrogens[0] = 1;

    let output = adjust_hs(&input, &original).unwrap();

    assert_eq!(input, snapshot);
    assert_eq!(output.topology.atoms[0].explicit_hydrogens(), 1);
    assert!(output.topology.atoms[0].is_aromatic());
    assert_eq!(output.topology.bonds, input.bonds);
    assert_eq!(output.topology.adjacency, input.adjacency);
    assert_eq!(output.topology.stereo_groups, input.stereo_groups);
    assert_eq!(output.topology.substance_groups, input.substance_groups);
}

#[test]
fn adjust_hs_leaves_no_implicit_dummy_and_charged_atoms_unchanged() {
    let atoms = vec![
        Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C).with_no_implicit(true),
        ),
        Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::DUMMY)),
        Atom::from_spec(
            AtomId::new(2),
            AtomSpec::new(Element::H).with_formal_charge(1),
        ),
    ];
    let input = TopologyBlock::try_from_parts(atoms, Vec::new(), Vec::new(), Vec::new()).unwrap();
    let original =
        assign_valence_with_options_for_topology(&input, ValenceModel::RdkitLike, false).unwrap();
    assert_eq!(original.implicit_hydrogens, vec![0, 0, 0]);

    let output = adjust_hs(&input, &original).unwrap();

    assert_eq!(output.topology, input);
    assert_eq!(output.valence, original);
}

#[test]
fn adjust_hs_reports_exact_overflow_and_assignment_errors_without_partial_output() {
    let overflow_input = TopologyBlock::try_from_parts(
        vec![Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C).with_explicit_hydrogens(u8::MAX),
        )],
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .unwrap();
    assert_eq!(
        adjust_hs(
            &overflow_input,
            &ValenceAssignment {
                explicit_valence: vec![i32::from(u8::MAX)],
                implicit_hydrogens: vec![1],
            },
        ),
        Err(AdjustHsError::ExplicitHydrogenOverflow {
            atom: AtomId::new(0),
            original_explicit_hydrogens: u8::MAX,
            original_implicit_valence: 1,
            recalculated_implicit_valence: 0,
        })
    );
    assert_eq!(overflow_input.atoms[0].explicit_hydrogens(), u8::MAX);

    let input = topology(2, Vec::new(), Vec::new());
    assert_eq!(
        adjust_hs(
            &input,
            &ValenceAssignment {
                explicit_valence: vec![0],
                implicit_hydrogens: vec![4, 4],
            },
        ),
        Err(AdjustHsError::ValenceAssignmentLength {
            field: "explicit_valence",
            actual: 1,
            expected: 2,
        })
    );
    assert_eq!(
        adjust_hs(
            &input,
            &ValenceAssignment {
                explicit_valence: vec![0, 0],
                implicit_hydrogens: vec![4, -1],
            },
        ),
        Err(AdjustHsError::InvalidValenceRow {
            atom: AtomId::new(1),
            field: "implicit_hydrogens",
            value: -1,
        })
    );

    let invalid = TopologyBlock {
        adjacency: AdjacencyList::from_topology(0, &[]),
        ..input
    };
    assert!(matches!(
        adjust_hs(
            &invalid,
            &ValenceAssignment {
                explicit_valence: vec![0, 0],
                implicit_hydrogens: vec![4, 4],
            },
        ),
        Err(AdjustHsError::InvalidTopology { .. })
    ));
}

#[test]
fn sanitize_pipeline_flag_vocabulary_default_unknown_bits_and_empty_stage_matrix_are_exact() {
    let cases = [
        (SanitizeOperations::NONE, SanitizeStage::None, 0x000),
        (SanitizeOperations::CLEANUP, SanitizeStage::Cleanup, 0x001),
        (
            SanitizeOperations::PROPERTIES,
            SanitizeStage::Properties,
            0x002,
        ),
        (
            SanitizeOperations::SYMM_RINGS,
            SanitizeStage::SymmRings,
            0x004,
        ),
        (SanitizeOperations::KEKULIZE, SanitizeStage::Kekulize, 0x008),
        (
            SanitizeOperations::FIND_RADICALS,
            SanitizeStage::FindRadicals,
            0x010,
        ),
        (
            SanitizeOperations::SET_AROMATICITY,
            SanitizeStage::SetAromaticity,
            0x020,
        ),
        (
            SanitizeOperations::SET_CONJUGATION,
            SanitizeStage::SetConjugation,
            0x040,
        ),
        (
            SanitizeOperations::SET_HYBRIDIZATION,
            SanitizeStage::SetHybridization,
            0x080,
        ),
        (
            SanitizeOperations::CLEANUP_CHIRALITY,
            SanitizeStage::CleanupChirality,
            0x100,
        ),
        (
            SanitizeOperations::ADJUST_HS,
            SanitizeStage::AdjustHs,
            0x200,
        ),
        (
            SanitizeOperations::CLEANUP_ORGANOMETALLICS,
            SanitizeStage::CleanupOrganometallics,
            0x400,
        ),
        (
            SanitizeOperations::CLEANUP_ATROPISOMERS,
            SanitizeStage::CleanupAtropisomers,
            0x800,
        ),
    ];
    for (operation, stage, bits) in cases {
        assert_eq!(operation.bits(), bits);
        assert_eq!(stage as u32, bits);
        assert_eq!(SanitizeOperations::from_bits(bits).unwrap(), operation);
    }
    assert_eq!(SanitizeOperations::ALL.bits(), 0x0fff_ffff);
    assert_eq!(SanitizeOperations::default(), SanitizeOperations::ALL);
    assert_eq!(
        SanitizeParams::default().operations,
        SanitizeOperations::ALL
    );
    assert!(SanitizeOperations::NONE.is_empty());
    assert!(SanitizeOperations::ALL.contains(SanitizeOperations::CLEANUP_ATROPISOMERS));
    let selected = SanitizeOperations::CLEANUP | SanitizeOperations::KEKULIZE;
    assert_eq!(selected.bits(), 0x009);
    assert!(selected.contains(SanitizeOperations::CLEANUP));
    assert!(!selected.contains(SanitizeOperations::PROPERTIES));
    assert!(matches!(
        SanitizeOperations::from_bits(0x1000),
        Err(SanitizeError::InvalidOperations {
            bits: 0x1000,
            unknown_bits: 0x1000,
        })
    ));
    assert!(matches!(
        SanitizeOperations::from_bits(0x1000_0001),
        Err(SanitizeError::InvalidOperations {
            bits: 0x1000_0001,
            unknown_bits: 0x1000_0000,
        })
    ));
    assert_eq!(
        SanitizeOperations::from_bits(0x0fff_ffff).unwrap(),
        SanitizeOperations::ALL
    );

    let empty = TopologyBlock::default();
    for operations in cases
        .into_iter()
        .map(|(operation, _, _)| operation)
        .chain([SanitizeOperations::ALL])
    {
        assert_eq!(
            sanitize_topology(&empty, &SanitizeParams { operations })
                .unwrap()
                .topology,
            empty,
            "isolated operation bits 0x{:x}",
            operations.bits()
        );
    }
}

#[test]
fn sanitize_pipeline_clears_only_computed_properties_and_reports_property_before_kekulize() {
    let atom = AtomSpec::new(Element::C)
        .with_prop("user_atom", "keep")
        .unwrap()
        .with_computed_prop("computed_atom", "drop")
        .unwrap();
    let other = AtomSpec::new(Element::C);
    let bond = BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
        .with_prop("user_bond", "keep")
        .unwrap()
        .with_computed_prop("computed_bond", "drop")
        .unwrap();
    let input = topology_from_specs(vec![atom, other], vec![bond]);
    let snapshot = input.clone();
    let output = sanitize_topology(
        &input,
        &SanitizeParams {
            operations: SanitizeOperations::NONE,
        },
    )
    .unwrap()
    .topology;
    assert_eq!(input, snapshot);
    assert_eq!(output.atoms[0].prop("user_atom"), Some("keep"));
    assert_eq!(output.atoms[0].prop("computed_atom"), None);
    assert_eq!(output.bonds[0].prop("user_bond"), Some("keep"));
    assert_eq!(output.bonds[0].prop("computed_bond"), None);

    let mut overvalent_specs = vec![AtomSpec::new(Element::C)];
    overvalent_specs.extend((0..5).map(|_| AtomSpec::new(Element::H)));
    overvalent_specs.push(
        AtomSpec::new(Element::C)
            .with_aromatic(true)
            .with_no_implicit(true),
    );
    let overvalent = topology_from_specs(
        overvalent_specs,
        (1..=5)
            .map(|neighbor| bond_spec(0, neighbor, BondOrder::Single))
            .collect(),
    );
    let before = overvalent.clone();

    assert!(
        sanitize_topology(
            &overvalent,
            &SanitizeParams {
                operations: SanitizeOperations::NONE,
            },
        )
        .is_ok()
    );
    assert!(matches!(
        sanitize_topology(
            &overvalent,
            &SanitizeParams {
                operations: SanitizeOperations::PROPERTIES | SanitizeOperations::KEKULIZE,
            },
        ),
        Err(SanitizeError::Properties {
            stage: SanitizeStage::Properties,
            ..
        })
    ));
    assert_eq!(overvalent, before);

    let aromatic_not_in_ring = topology_from_specs(
        vec![
            AtomSpec::new(Element::C)
                .with_aromatic(true)
                .with_no_implicit(true),
        ],
        Vec::new(),
    );
    assert!(matches!(
        sanitize_topology(
            &aromatic_not_in_ring,
            &SanitizeParams {
                operations: SanitizeOperations::KEKULIZE,
            },
        ),
        Err(SanitizeError::Kekulize {
            stage: SanitizeStage::Kekulize,
            ..
        })
    ));

    let invalid = TopologyBlock {
        adjacency: AdjacencyList::from_topology(0, &[]),
        ..input
    };
    assert!(matches!(
        sanitize_topology(&invalid, &SanitizeParams::default()),
        Err(SanitizeError::InvalidTopology {
            stage: SanitizeStage::None,
            ..
        })
    ));
}

#[test]
fn sanitize_pipeline_cleanup_stages_are_independent_and_precede_strict_properties() {
    let combined = topology_from_specs(
        vec![
            atom_spec(6),
            atom_spec(7),
            atom_spec(8),
            atom_spec(8),
            atom_spec(7),
            atom_spec(6),
            atom_spec(6),
            atom_spec(6),
            atom_spec(26),
        ],
        vec![
            bond_spec(0, 1, BondOrder::Single),
            bond_spec(1, 2, BondOrder::Double),
            bond_spec(1, 3, BondOrder::Double),
            bond_spec(4, 5, BondOrder::Single),
            bond_spec(4, 6, BondOrder::Single),
            bond_spec(4, 7, BondOrder::Single),
            bond_spec(4, 8, BondOrder::Single),
        ],
    );
    let snapshot = combined.clone();

    let charge_only = sanitize_topology(
        &combined,
        &SanitizeParams {
            operations: SanitizeOperations::CLEANUP,
        },
    )
    .unwrap()
    .topology;
    assert_eq!(charge_only.atoms[1].formal_charge(), 1);
    assert_eq!(charge_only.bonds[1].order(), BondOrder::Single);
    assert_eq!(charge_only.bonds[6].order(), BondOrder::Single);

    let organometallic_only = sanitize_topology(
        &combined,
        &SanitizeParams {
            operations: SanitizeOperations::CLEANUP_ORGANOMETALLICS,
        },
    )
    .unwrap()
    .topology;
    assert_eq!(organometallic_only.atoms[1].formal_charge(), 0);
    assert_eq!(organometallic_only.bonds[1].order(), BondOrder::Double);
    assert_eq!(organometallic_only.bonds[6].order(), BondOrder::Dative);

    let both_with_strict = sanitize_topology(
        &combined,
        &SanitizeParams {
            operations: SanitizeOperations::CLEANUP
                | SanitizeOperations::CLEANUP_ORGANOMETALLICS
                | SanitizeOperations::PROPERTIES,
        },
    )
    .unwrap()
    .topology;
    assert_eq!(both_with_strict.atoms[1].formal_charge(), 1);
    assert_eq!(both_with_strict.bonds[1].order(), BondOrder::Single);
    assert_eq!(both_with_strict.bonds[6].order(), BondOrder::Dative);
    assert_eq!(combined, snapshot);
}

#[test]
fn sanitize_pipeline_materializes_each_perception_and_cleanup_stage_in_stable_order() {
    let aromatic_input = aromatic_cycle(6);
    let direct_kekulize = kekulize(
        &aromatic_input,
        &KekulizeParams {
            mark_atoms_bonds: true,
            canonical: false,
            max_backtracks: KekulizeParams::default().max_backtracks,
        },
    )
    .unwrap()
    .topology;
    let sanitized_kekulize = sanitize_topology(
        &aromatic_input,
        &SanitizeParams {
            operations: SanitizeOperations::KEKULIZE,
        },
    )
    .unwrap()
    .topology;
    assert_eq!(sanitized_kekulize, direct_kekulize);
    assert!(
        sanitized_kekulize
            .atoms
            .iter()
            .all(|atom| !atom.is_aromatic())
    );
    assert!(
        sanitized_kekulize
            .bonds
            .iter()
            .all(|bond| !bond.is_aromatic())
    );

    let radical_input = topology_from_specs(
        vec![
            AtomSpec::new(Element::C)
                .with_no_implicit(true)
                .with_explicit_hydrogens(3),
        ],
        Vec::new(),
    );
    let radical_output = sanitize_topology(
        &radical_input,
        &SanitizeParams {
            operations: SanitizeOperations::FIND_RADICALS,
        },
    )
    .unwrap()
    .topology;
    assert_eq!(radical_output.atoms[0].radical_electrons(), 1);

    let benzene = alternating_cycle((0..6).map(|_| AtomSpec::new(Element::C)).collect());
    let aromatic_output = sanitize_topology(
        &benzene,
        &SanitizeParams {
            operations: SanitizeOperations::SET_AROMATICITY,
        },
    )
    .unwrap()
    .topology;
    assert!(aromatic_output.atoms.iter().all(Atom::is_aromatic));
    assert!(aromatic_output.bonds.iter().all(Bond::is_aromatic));

    let conjugated_output = sanitize_topology(
        &benzene,
        &SanitizeParams {
            operations: SanitizeOperations::SET_CONJUGATION,
        },
    )
    .unwrap()
    .topology;
    assert!(conjugated_output.bonds.iter().all(Bond::is_conjugated));

    let ethene = topology_from_specs(
        vec![AtomSpec::new(Element::C), AtomSpec::new(Element::C)],
        vec![bond_spec(0, 1, BondOrder::Double)],
    );
    let hybridized_output = sanitize_topology(
        &ethene,
        &SanitizeParams {
            operations: SanitizeOperations::SET_HYBRIDIZATION,
        },
    )
    .unwrap()
    .topology;
    assert_eq!(
        hybridized_output.atoms[0].hybridization(),
        Hybridization::Sp2
    );
    assert_eq!(
        hybridized_output.atoms[1].hybridization(),
        Hybridization::Sp2
    );

    let atrop_input = topology_from_specs(
        vec![
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
        ],
        vec![
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
                .with_stereo(BondStereo::AtropCw),
        ],
    );
    let atrop_output = sanitize_topology(
        &atrop_input,
        &SanitizeParams {
            operations: SanitizeOperations::CLEANUP_ATROPISOMERS,
        },
    )
    .unwrap()
    .topology;
    assert_eq!(atrop_output.bonds[0].stereo(), BondStereo::None);

    let chirality_input = topology_from_specs(
        vec![
            AtomSpec::new(Element::C)
                .with_chiral_tag(ChiralTag::TetrahedralCw)
                .with_hybridization(Hybridization::Sp2),
        ],
        Vec::new(),
    );
    let chirality_output = sanitize_topology(
        &chirality_input,
        &SanitizeParams {
            operations: SanitizeOperations::CLEANUP_CHIRALITY,
        },
    )
    .unwrap()
    .topology;
    assert_eq!(
        chirality_output.atoms[0].chiral_tag(),
        ChiralTag::Unspecified
    );

    let pyrrole = topology_from_specs(
        std::iter::once(AtomSpec::new(Element::N))
            .chain((0..4).map(|_| AtomSpec::new(Element::C)))
            .collect(),
        vec![
            bond_spec(0, 1, BondOrder::Single),
            bond_spec(1, 2, BondOrder::Double),
            bond_spec(2, 3, BondOrder::Single),
            bond_spec(3, 4, BondOrder::Double),
            bond_spec(4, 0, BondOrder::Single),
        ],
    );
    let adjusted = sanitize_topology(
        &pyrrole,
        &SanitizeParams {
            operations: SanitizeOperations::SET_AROMATICITY | SanitizeOperations::ADJUST_HS,
        },
    )
    .unwrap()
    .topology;
    assert!(adjusted.atoms[0].is_aromatic());
    assert_eq!(adjusted.atoms[0].explicit_hydrogens(), 1);
}

#[test]
fn sanitize_pipeline_all_matches_all_named_stages_and_is_deterministic_and_atomic() {
    let benzene = alternating_cycle((0..6).map(|_| AtomSpec::new(Element::C)).collect());
    let snapshot = benzene.clone();
    let named = SanitizeOperations::CLEANUP
        | SanitizeOperations::PROPERTIES
        | SanitizeOperations::SYMM_RINGS
        | SanitizeOperations::KEKULIZE
        | SanitizeOperations::FIND_RADICALS
        | SanitizeOperations::SET_AROMATICITY
        | SanitizeOperations::SET_CONJUGATION
        | SanitizeOperations::SET_HYBRIDIZATION
        | SanitizeOperations::CLEANUP_CHIRALITY
        | SanitizeOperations::ADJUST_HS
        | SanitizeOperations::CLEANUP_ORGANOMETALLICS
        | SanitizeOperations::CLEANUP_ATROPISOMERS;
    assert_eq!(named.bits(), 0x0fff);

    let first = sanitize_topology(
        &benzene,
        &SanitizeParams {
            operations: SanitizeOperations::ALL,
        },
    )
    .unwrap()
    .topology;
    let second = sanitize_topology(
        &benzene,
        &SanitizeParams {
            operations: SanitizeOperations::ALL,
        },
    )
    .unwrap()
    .topology;
    let named_output = sanitize_topology(&benzene, &SanitizeParams { operations: named })
        .unwrap()
        .topology;
    assert_eq!(first, second);
    assert_eq!(first, named_output);
    assert_eq!(benzene, snapshot);
    assert_eq!(
        first.atoms.iter().map(Atom::id).collect::<Vec<_>>(),
        (0..6).map(AtomId::new).collect::<Vec<_>>()
    );
    assert_eq!(
        first.bonds.iter().map(Bond::id).collect::<Vec<_>>(),
        (0..6).map(BondId::new).collect::<Vec<_>>()
    );
    assert_eq!(first.adjacency, benzene.adjacency);
    first.validate().unwrap();
}

#[test]
fn chemistry_problem_detection_returns_zero_for_valid_benzene_and_preserves_input() {
    let input = aromatic_cycle(6);
    let snapshot = input.clone();

    let first = detect_chemistry_problems(&input, &SanitizeParams::default()).unwrap();
    let second = detect_chemistry_problems(&input, &SanitizeParams::default()).unwrap();

    assert!(first.problems.is_empty());
    assert!(second.problems.is_empty());
    assert_eq!(input, snapshot);
}

#[test]
fn chemistry_problem_detection_collects_source_rows_then_kekulize_in_exact_order() {
    // Unsanitized source case: CO(C)CFCc1cc1. Atom 1 is overvalent O, atom 4
    // is overvalent F, and atoms 6/7/8 form the non-kekulizable aromatic ring.
    let input = topology_from_specs(
        vec![
            atom_spec(6),
            atom_spec(8),
            atom_spec(6),
            atom_spec(6),
            atom_spec(9),
            atom_spec(6),
            AtomSpec::new(Element::C).with_aromatic(true),
            AtomSpec::new(Element::C).with_aromatic(true),
            AtomSpec::new(Element::C).with_aromatic(true),
        ],
        vec![
            bond_spec(0, 1, BondOrder::Single),
            bond_spec(1, 2, BondOrder::Single),
            bond_spec(1, 3, BondOrder::Single),
            bond_spec(3, 4, BondOrder::Single),
            bond_spec(4, 5, BondOrder::Single),
            bond_spec(5, 6, BondOrder::Single),
            BondSpec::new(AtomId::new(6), AtomId::new(7), BondOrder::Aromatic).with_aromatic(true),
            BondSpec::new(AtomId::new(7), AtomId::new(8), BondOrder::Aromatic).with_aromatic(true),
            BondSpec::new(AtomId::new(8), AtomId::new(6), BondOrder::Aromatic).with_aromatic(true),
        ],
    );
    let snapshot = input.clone();

    let report = detect_chemistry_problems(&input, &SanitizeParams::default()).unwrap();
    assert_eq!(report.problems.len(), 3);
    assert_eq!(report.problems[0].operation, SanitizeStage::Properties);
    assert!(matches!(
        &report.problems[0].error,
        ChemistryProblemError::Valence(ValenceError::InvalidValence {
            atom,
            atomic_number: 8,
            formal_charge: 0,
            phase: ValencePhase::Explicit,
            calculated: Some(3),
            ..
        }) if *atom == AtomId::new(1)
    ));
    assert_eq!(report.problems[1].operation, SanitizeStage::Properties);
    assert!(matches!(
        &report.problems[1].error,
        ChemistryProblemError::Valence(ValenceError::InvalidValence {
            atom,
            atomic_number: 9,
            formal_charge: 0,
            phase: ValencePhase::Explicit,
            calculated: Some(2),
            ..
        }) if *atom == AtomId::new(4)
    ));
    assert_eq!(report.problems[2].operation, SanitizeStage::Kekulize);
    assert!(matches!(
        &report.problems[2].error,
        ChemistryProblemError::Kekulize(KekulizeError::NotKekulizable {
            problem_atoms,
        }) if problem_atoms == &vec![AtomId::new(6), AtomId::new(7), AtomId::new(8)]
    ));

    let repeated = detect_chemistry_problems(&input, &SanitizeParams::default()).unwrap();
    assert_eq!(format!("{report:?}"), format!("{repeated:?}"));
    assert_eq!(input, snapshot);
}

#[test]
fn chemistry_problem_detection_reports_non_ring_aromatic_atom_as_typed_kekulize_problem() {
    let input = topology_from_specs(
        vec![
            AtomSpec::new(Element::C)
                .with_aromatic(true)
                .with_no_implicit(true),
        ],
        Vec::new(),
    );
    let snapshot = input.clone();
    let report = detect_chemistry_problems(
        &input,
        &SanitizeParams {
            operations: SanitizeOperations::KEKULIZE,
        },
    )
    .unwrap();

    assert_eq!(report.problems.len(), 1);
    assert_eq!(report.problems[0].operation, SanitizeStage::Kekulize);
    assert!(matches!(
        &report.problems[0].error,
        ChemistryProblemError::Kekulize(KekulizeError::AromaticAtomOutsideRing { atom })
            if *atom == AtomId::new(0)
    ));
    assert_eq!(input, snapshot);
}

#[test]
fn chemistry_problem_detection_propagates_cleanup_and_invalid_input_errors() {
    let cleanup_failure = topology_from_specs(
        vec![atom_spec(7), atom_spec(6)],
        vec![bond_spec(0, 1, BondOrder::Other)],
    );
    assert!(matches!(
        detect_chemistry_problems(
            &cleanup_failure,
            &SanitizeParams {
                operations: SanitizeOperations::CLEANUP,
            },
        ),
        Err(SanitizeError::Cleanup {
            stage: SanitizeStage::Cleanup,
            ..
        })
    ));

    let valid = topology_from_specs(vec![atom_spec(6)], Vec::new());
    let invalid = TopologyBlock {
        adjacency: AdjacencyList::from_topology(0, &[]),
        ..valid
    };
    assert!(matches!(
        detect_chemistry_problems(&invalid, &SanitizeParams::default()),
        Err(SanitizeError::InvalidTopology {
            stage: SanitizeStage::None,
            ..
        })
    ));
}

#[test]
fn chemistry_problem_detection_honors_flag_independence_and_non_strict_fallback() {
    let organometallic_candidate = topology_from_specs(
        vec![
            atom_spec(7),
            atom_spec(6),
            atom_spec(6),
            atom_spec(6),
            atom_spec(26),
        ],
        vec![
            bond_spec(0, 1, BondOrder::Single),
            bond_spec(0, 2, BondOrder::Single),
            bond_spec(0, 3, BondOrder::Single),
            bond_spec(0, 4, BondOrder::Single),
        ],
    );
    let report = detect_chemistry_problems(
        &organometallic_candidate,
        &SanitizeParams {
            operations: SanitizeOperations::CLEANUP_ORGANOMETALLICS
                | SanitizeOperations::PROPERTIES,
        },
    )
    .unwrap();
    assert_eq!(report.problems.len(), 1);
    assert!(matches!(
        &report.problems[0].error,
        ChemistryProblemError::Valence(ValenceError::InvalidValence { atom, .. })
            if *atom == AtomId::new(0)
    ));

    let mut atoms = vec![atom_spec(6)];
    atoms.extend((0..5).map(|_| atom_spec(1)));
    let overvalent = topology_from_specs(
        atoms,
        (1..=5)
            .map(|neighbor| bond_spec(0, neighbor, BondOrder::Single))
            .collect(),
    );
    let ignored_stages = SanitizeOperations::SYMM_RINGS
        | SanitizeOperations::FIND_RADICALS
        | SanitizeOperations::SET_AROMATICITY
        | SanitizeOperations::SET_CONJUGATION
        | SanitizeOperations::SET_HYBRIDIZATION
        | SanitizeOperations::CLEANUP_CHIRALITY
        | SanitizeOperations::ADJUST_HS
        | SanitizeOperations::CLEANUP_ATROPISOMERS;
    let report = detect_chemistry_problems(
        &overvalent,
        &SanitizeParams {
            operations: ignored_stages,
        },
    )
    .unwrap();
    assert!(report.problems.is_empty());
}
