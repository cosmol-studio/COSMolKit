use cosmolkit_core::{
    PotentialStereoCenter, PotentialStereoDescriptor, PotentialStereoError, PotentialStereoParams,
    PotentialStereoSpecified, PotentialStereoType, RingFindType, RingInfo, RingSearchParams,
    ValenceAssignment, potential_stereo, symmetrized_sssr,
};
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, TopologyBlock, TopologyValidationError,
};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo, ChiralTag, Element};

fn atom(id: usize, element: Element) -> Atom {
    Atom::from_spec(AtomId::new(id), AtomSpec::new(element))
}

fn atom_with(id: usize, spec: AtomSpec) -> Atom {
    Atom::from_spec(AtomId::new(id), spec)
}

fn bond(id: usize, begin: usize, end: usize, order: BondOrder) -> Bond {
    Bond::from_spec(
        BondId::new(id),
        BondSpec::new(AtomId::new(begin), AtomId::new(end), order),
    )
}

fn topology(atoms: Vec<Atom>, bonds: Vec<Bond>) -> TopologyBlock {
    TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
}

fn valence(explicit: Vec<i32>, implicit: Vec<i32>) -> ValenceAssignment {
    ValenceAssignment {
        explicit_valence: explicit,
        implicit_hydrogens: implicit,
    }
}

fn symmetric_empty_rings(topology: &TopologyBlock) -> RingInfo {
    RingInfo::new(
        RingFindType::SymmSssr,
        topology.atoms.len(),
        topology.bonds.len(),
    )
}

fn four_distinct_ligands(tag: ChiralTag, unknown: bool) -> (TopologyBlock, ValenceAssignment) {
    let center = AtomSpec::new(Element::C)
        .with_chiral_tag(tag)
        .with_unknown_stereo(unknown);
    let topology = topology(
        vec![
            atom_with(0, center),
            atom(1, Element::F),
            atom(2, Element::CL),
            atom(3, Element::BR),
            atom(4, Element::I),
        ],
        vec![
            bond(0, 0, 1, BondOrder::Single),
            bond(1, 0, 2, BondOrder::Single),
            bond(2, 0, 3, BondOrder::Single),
            bond(3, 0, 4, BondOrder::Single),
        ],
    );
    (topology, valence(vec![4, 1, 1, 1, 1], vec![0; 5]))
}

fn alkene(stereo: BondStereo, unknown_direction: bool) -> (TopologyBlock, ValenceAssignment) {
    let mut double =
        BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Double).with_stereo(stereo);
    if matches!(
        stereo,
        BondStereo::Cis | BondStereo::Trans | BondStereo::E | BondStereo::Z
    ) {
        double = double.with_stereo_atoms(AtomId::new(0), AtomId::new(3));
    }
    if unknown_direction {
        double = double.with_direction(BondDirection::EitherDouble);
    }
    let topology = topology(
        (0..4).map(|index| atom(index, Element::C)).collect(),
        vec![
            bond(0, 0, 1, BondOrder::Single),
            Bond::from_spec(BondId::new(1), double),
            bond(2, 2, 3, BondOrder::Single),
        ],
    );
    (topology, valence(vec![1, 3, 3, 1], vec![3, 1, 1, 3]))
}

#[test]
fn defaults_empty_and_disconnected_results_are_total_and_deterministic() {
    assert_eq!(
        PotentialStereoParams::default(),
        PotentialStereoParams {
            clean: false,
            flag_possible: true,
            allow_nontetrahedral: true,
        }
    );
    let empty = TopologyBlock::default();
    let first = potential_stereo(
        &empty,
        &valence(vec![], vec![]),
        &symmetric_empty_rings(&empty),
        &PotentialStereoParams::default(),
    )
    .unwrap();
    assert!(first.stereo.is_empty());
    assert!(first.atom_ranks.is_empty());
    assert!(first.ring_relations.is_empty());
    assert!(first.cleaned_topology.is_none());

    let disconnected = topology(vec![atom(0, Element::C), atom(1, Element::C)], vec![]);
    let assignment = valence(vec![0, 0], vec![4, 4]);
    let rings = symmetric_empty_rings(&disconnected);
    let left = potential_stereo(
        &disconnected,
        &assignment,
        &rings,
        &PotentialStereoParams::default(),
    )
    .unwrap();
    let right = potential_stereo(
        &disconnected,
        &assignment,
        &rings,
        &PotentialStereoParams::default(),
    )
    .unwrap();
    assert_eq!(left, right);
    assert_eq!(left.atom_ranks, vec![0, 0]);
}

#[test]
fn atom_rows_cover_possible_specified_unknown_and_flag_possible_false() {
    let (possible_topology, assignment) = four_distinct_ligands(ChiralTag::Unspecified, false);
    let possible = potential_stereo(
        &possible_topology,
        &assignment,
        &symmetric_empty_rings(&possible_topology),
        &PotentialStereoParams::default(),
    )
    .unwrap();
    assert_eq!(possible.stereo.len(), 1);
    assert_eq!(
        possible.stereo[0].centered_on,
        PotentialStereoCenter::Atom(AtomId::new(0))
    );
    assert_eq!(
        possible.stereo[0].stereo_type,
        PotentialStereoType::AtomTetrahedral
    );
    assert_eq!(
        possible.stereo[0].specified,
        PotentialStereoSpecified::Unspecified
    );
    assert_eq!(possible.stereo[0].controlling_atoms.len(), 4);

    let hidden = potential_stereo(
        &possible_topology,
        &assignment,
        &symmetric_empty_rings(&possible_topology),
        &PotentialStereoParams {
            flag_possible: false,
            ..PotentialStereoParams::default()
        },
    )
    .unwrap();
    assert!(hidden.stereo.is_empty());

    let (specified_topology, specified_valence) =
        four_distinct_ligands(ChiralTag::TetrahedralCw, false);
    let specified = potential_stereo(
        &specified_topology,
        &specified_valence,
        &symmetric_empty_rings(&specified_topology),
        &PotentialStereoParams {
            flag_possible: false,
            ..PotentialStereoParams::default()
        },
    )
    .unwrap();
    assert_eq!(specified.stereo.len(), 1);
    assert_eq!(
        specified.stereo[0].specified,
        PotentialStereoSpecified::Specified
    );
    assert!(matches!(
        specified.stereo[0].descriptor,
        PotentialStereoDescriptor::TetrahedralClockwise
            | PotentialStereoDescriptor::TetrahedralCounterclockwise
    ));

    let (unknown_topology, unknown_valence) = four_distinct_ligands(ChiralTag::Unspecified, true);
    let unknown = potential_stereo(
        &unknown_topology,
        &unknown_valence,
        &symmetric_empty_rings(&unknown_topology),
        &PotentialStereoParams {
            flag_possible: false,
            ..PotentialStereoParams::default()
        },
    )
    .unwrap();
    assert_eq!(
        unknown.stereo[0].specified,
        PotentialStereoSpecified::Unknown
    );
}

#[test]
fn duplicate_ligands_are_rejected_and_cleaning_is_detached_and_complete() {
    let mut center = AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCw);
    center = center.with_explicit_hydrogens(0);
    let mut wedge = BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
        .with_direction(BondDirection::BeginWedge);
    wedge = wedge.with_stereo(BondStereo::None);
    let source = topology(
        vec![
            atom_with(0, center),
            atom(1, Element::F),
            atom(2, Element::F),
            atom(3, Element::CL),
            atom(4, Element::BR),
        ],
        vec![
            Bond::from_spec(BondId::new(0), wedge),
            bond(1, 0, 2, BondOrder::Single),
            bond(2, 0, 3, BondOrder::Single),
            bond(3, 0, 4, BondOrder::Single),
        ],
    );
    let original = source.clone();
    let result = potential_stereo(
        &source,
        &valence(vec![4, 1, 1, 1, 1], vec![0; 5]),
        &symmetric_empty_rings(&source),
        &PotentialStereoParams {
            clean: true,
            flag_possible: true,
            allow_nontetrahedral: true,
        },
    )
    .unwrap();
    assert!(result.stereo.is_empty());
    assert_eq!(source, original);
    let cleaned = result.cleaned_topology.unwrap();
    assert_eq!(cleaned.atoms[0].chiral_tag(), ChiralTag::Unspecified);
    assert_eq!(cleaned.bonds[0].direction(), BondDirection::None);
    assert_eq!(cleaned.atoms.len(), source.atoms.len());
    assert_eq!(cleaned.bonds.len(), source.bonds.len());
    cleaned.validate().unwrap();
}

#[test]
fn all_clean_and_flag_possible_combinations_have_explicit_semantics() {
    let (topology, assignment) = four_distinct_ligands(ChiralTag::Unspecified, false);
    let rings = symmetric_empty_rings(&topology);
    for clean in [false, true] {
        for flag_possible in [false, true] {
            let result = potential_stereo(
                &topology,
                &assignment,
                &rings,
                &PotentialStereoParams {
                    clean,
                    flag_possible,
                    allow_nontetrahedral: true,
                },
            )
            .unwrap();
            assert_eq!(result.cleaned_topology.is_some(), clean);
            assert_eq!(result.stereo.len(), usize::from(flag_possible));
        }
    }
}

#[test]
fn double_bond_rows_cover_possible_specified_and_unknown_in_atom_then_bond_order() {
    for (stereo, unknown_direction, specified, descriptor) in [
        (
            BondStereo::None,
            false,
            PotentialStereoSpecified::Unspecified,
            PotentialStereoDescriptor::None,
        ),
        (
            BondStereo::Cis,
            false,
            PotentialStereoSpecified::Specified,
            PotentialStereoDescriptor::BondCis,
        ),
        (
            BondStereo::Any,
            true,
            PotentialStereoSpecified::Unknown,
            PotentialStereoDescriptor::None,
        ),
    ] {
        let (topology, assignment) = alkene(stereo, unknown_direction);
        let result = potential_stereo(
            &topology,
            &assignment,
            &symmetric_empty_rings(&topology),
            &PotentialStereoParams::default(),
        )
        .unwrap();
        assert_eq!(result.stereo.len(), 1);
        let row = &result.stereo[0];
        assert_eq!(row.centered_on, PotentialStereoCenter::Bond(BondId::new(1)));
        assert_eq!(row.stereo_type, PotentialStereoType::BondDouble);
        assert_eq!(row.specified, specified);
        assert_eq!(row.descriptor, descriptor);
        assert_eq!(row.controlling_atoms.len(), 4);
    }
}

#[test]
fn rings_below_eight_exclude_double_bonds_and_eight_membered_rings_allow_them() {
    for (size, expected_rows) in [(7usize, 0usize), (8, 1)] {
        let mut bonds = Vec::new();
        for index in 0..size {
            bonds.push(bond(
                index,
                index,
                (index + 1) % size,
                if index == 0 {
                    BondOrder::Double
                } else {
                    BondOrder::Single
                },
            ));
        }
        let topology = topology(
            (0..size).map(|index| atom(index, Element::C)).collect(),
            bonds,
        );
        let mut explicit = vec![2; size];
        let mut implicit = vec![2; size];
        explicit[0] = 3;
        explicit[1] = 3;
        implicit[0] = 1;
        implicit[1] = 1;
        let rings = symmetrized_sssr(&topology, &RingSearchParams::default()).unwrap();
        let result = potential_stereo(
            &topology,
            &valence(explicit, implicit),
            &rings,
            &PotentialStereoParams::default(),
        )
        .unwrap();
        assert_eq!(result.stereo.len(), expected_rows);
    }
}

#[test]
fn non_tetrahedral_parameter_controls_typed_interpretation() {
    let center = AtomSpec::new(Element::PT)
        .with_chiral_tag(ChiralTag::SquarePlanar)
        .with_chiral_permutation(1);
    let topology = topology(
        vec![
            atom_with(0, center),
            atom(1, Element::F),
            atom(2, Element::CL),
            atom(3, Element::BR),
            atom(4, Element::I),
        ],
        vec![
            bond(0, 0, 1, BondOrder::Single),
            bond(1, 0, 2, BondOrder::Single),
            bond(2, 0, 3, BondOrder::Single),
            bond(3, 0, 4, BondOrder::Single),
        ],
    );
    let assignment = valence(vec![4, 1, 1, 1, 1], vec![0; 5]);
    let rings = symmetric_empty_rings(&topology);
    let enabled = potential_stereo(
        &topology,
        &assignment,
        &rings,
        &PotentialStereoParams::default(),
    )
    .unwrap();
    assert_eq!(
        enabled.stereo[0].stereo_type,
        PotentialStereoType::AtomSquarePlanar
    );
    assert_eq!(enabled.stereo[0].permutation, 1);
    assert_eq!(
        enabled.stereo[0].specified,
        PotentialStereoSpecified::Specified
    );

    let disabled = potential_stereo(
        &topology,
        &assignment,
        &rings,
        &PotentialStereoParams {
            allow_nontetrahedral: false,
            ..PotentialStereoParams::default()
        },
    )
    .unwrap();
    assert!(
        disabled
            .stereo
            .iter()
            .all(|row| row.stereo_type != PotentialStereoType::AtomSquarePlanar)
    );
}

#[test]
fn ring_special_cases_return_stable_signed_bidirectional_relations() {
    let mut atoms = (0..8)
        .map(|index| atom(index, Element::C))
        .collect::<Vec<_>>();
    atoms[0] = atom_with(
        0,
        AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCw),
    );
    atoms[3] = atom_with(
        3,
        AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCcw),
    );
    let mut bonds = (0..6)
        .map(|index| bond(index, index, (index + 1) % 6, BondOrder::Single))
        .collect::<Vec<_>>();
    bonds.push(bond(6, 0, 6, BondOrder::Single));
    bonds.push(bond(7, 3, 7, BondOrder::Single));
    let topology = topology(atoms, bonds);
    let rings = symmetrized_sssr(&topology, &RingSearchParams::default()).unwrap();
    let result = potential_stereo(
        &topology,
        &valence(vec![3, 2, 2, 3, 2, 2, 1, 1], vec![1, 2, 2, 1, 2, 2, 3, 3]),
        &rings,
        &PotentialStereoParams::default(),
    )
    .unwrap();
    assert_eq!(result.ring_relations.len(), 2);
    assert_eq!(result.ring_relations[0].atom, AtomId::new(0));
    assert_eq!(result.ring_relations[0].other, AtomId::new(3));
    assert!(!result.ring_relations[0].same_orientation);
    assert_eq!(result.ring_relations[1].atom, AtomId::new(3));
    assert_eq!(result.ring_relations[1].other, AtomId::new(0));
    assert!(!result.ring_relations[1].same_orientation);
}

#[test]
fn validation_rejects_topology_valence_and_ring_rows_before_indexing() {
    let (topology, assignment) = four_distinct_ligands(ChiralTag::Unspecified, false);
    let mut invalid = topology.clone();
    invalid.atoms[0] = invalid.atoms[0].clone().with_id(AtomId::new(1));
    assert_eq!(
        potential_stereo(
            &invalid,
            &assignment,
            &symmetric_empty_rings(&topology),
            &PotentialStereoParams::default(),
        ),
        Err(PotentialStereoError::InvalidTopology(
            TopologyValidationError::AtomIdMismatch {
                position: 0,
                id: AtomId::new(1),
            }
        ))
    );
    assert_eq!(
        potential_stereo(
            &topology,
            &valence(vec![0; 4], vec![0; 5]),
            &symmetric_empty_rings(&topology),
            &PotentialStereoParams::default(),
        ),
        Err(PotentialStereoError::InvalidValence {
            field: "explicit_valence",
            actual: 4,
            atom_count: 5,
        })
    );
    assert!(matches!(
        potential_stereo(
            &topology,
            &assignment,
            &RingInfo::new(RingFindType::SymmSssr, 4, topology.bonds.len()),
            &PotentialStereoParams::default(),
        ),
        Err(PotentialStereoError::InvalidRingInfo {
            reason: "atom membership row count mismatch",
            value: 4,
            limit: 5,
            ..
        })
    ));
}

#[test]
fn atropisomer_state_fails_closed_until_its_owner_unit_is_integrated() {
    let mut spec = BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single);
    spec = spec.with_stereo(BondStereo::AtropCw);
    let topology = topology(
        vec![atom(0, Element::C), atom(1, Element::C)],
        vec![Bond::from_spec(BondId::new(0), spec)],
    );
    assert_eq!(
        potential_stereo(
            &topology,
            &valence(vec![1, 1], vec![3, 3]),
            &symmetric_empty_rings(&topology),
            &PotentialStereoParams::default(),
        ),
        Err(PotentialStereoError::AtropisomerDependencyUnavailable {
            bond: BondId::new(0),
        })
    );
}
