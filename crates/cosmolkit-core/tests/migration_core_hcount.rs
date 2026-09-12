use cosmolkit_core::{
    ValenceAssignment, ValenceError, ValenceParams, assign_valence, total_hydrogen_count,
};
use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, TopologyBlock,
    TopologyValidationError,
};
use cosmolkit_types::{BondOrder, Element};

fn atom_spec(element: Element) -> AtomSpec {
    AtomSpec::new(element)
}

fn topology(atom_specs: Vec<AtomSpec>, bond_specs: Vec<BondSpec>) -> TopologyBlock {
    let atoms = atom_specs
        .into_iter()
        .enumerate()
        .map(|(index, spec)| Atom::from_spec(AtomId::new(index), spec))
        .collect::<Vec<_>>();
    let bonds = bond_specs
        .into_iter()
        .enumerate()
        .map(|(index, spec)| Bond::from_spec(BondId::new(index), spec))
        .collect::<Vec<_>>();
    let adjacency = AdjacencyList::from_topology(atoms.len(), &bonds);
    TopologyBlock {
        atoms,
        bonds,
        adjacency,
        ..TopologyBlock::default()
    }
}

fn assignment(implicit_hydrogens: Vec<i32>) -> ValenceAssignment {
    ValenceAssignment {
        explicit_valence: Vec::new(),
        implicit_hydrogens,
    }
}

#[test]
fn source_carbon_graph_distinguishes_default_and_neighbor_inclusive_counts() {
    let input = topology(
        vec![
            atom_spec(Element::C),
            atom_spec(Element::C),
            atom_spec(Element::C),
            atom_spec(Element::C),
            atom_spec(Element::H),
        ],
        vec![
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Single),
            BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single),
            BondSpec::new(AtomId::new(2), AtomId::new(3), BondOrder::Single),
            BondSpec::new(AtomId::new(0), AtomId::new(4), BondOrder::Single),
        ],
    );
    let valence = assign_valence(&input, &ValenceParams::default()).unwrap();

    assert_eq!(
        total_hydrogen_count(&input, &valence, AtomId::new(0), false).unwrap(),
        1
    );
    assert_eq!(
        total_hydrogen_count(&input, &valence, AtomId::new(0), true).unwrap(),
        2
    );
    assert_eq!(
        total_hydrogen_count(&input, &valence, AtomId::new(2), false).unwrap(),
        1
    );
}

#[test]
fn explicit_implicit_and_neighbor_components_are_additive_without_mutation() {
    let input = topology(
        vec![
            atom_spec(Element::N).with_explicit_hydrogens(2),
            atom_spec(Element::H),
        ],
        vec![BondSpec::new(
            AtomId::new(0),
            AtomId::new(1),
            BondOrder::Single,
        )],
    );
    let before = input.clone();
    let valence = assignment(vec![3, 0]);

    assert_eq!(
        total_hydrogen_count(&input, &valence, AtomId::new(0), false).unwrap(),
        5
    );
    assert_eq!(
        total_hydrogen_count(&input, &valence, AtomId::new(0), true).unwrap(),
        6
    );
    assert_eq!(input, before);
}

#[test]
fn no_implicit_bypasses_missing_and_negative_assignment_rows() {
    let input = topology(
        vec![
            atom_spec(Element::C)
                .with_no_implicit(true)
                .with_explicit_hydrogens(2)
                .with_implicit_hydrogen(true),
        ],
        vec![],
    );
    assert_eq!(
        total_hydrogen_count(&input, &assignment(vec![]), AtomId::new(0), false).unwrap(),
        2
    );
    assert_eq!(
        total_hydrogen_count(&input, &assignment(vec![-1]), AtomId::new(0), false).unwrap(),
        2
    );
}

#[test]
fn missing_and_negative_implicit_rows_are_exact_uninitialized_errors() {
    let input = topology(vec![atom_spec(Element::C)], vec![]);
    for valence in [assignment(vec![]), assignment(vec![-1])] {
        assert_eq!(
            total_hydrogen_count(&input, &valence, AtomId::new(0), false),
            Err(ValenceError::ImplicitValenceCacheNotInitialized {
                atom: AtomId::new(0),
            })
        );
    }
    assert_eq!(
        total_hydrogen_count(&input, &assignment(vec![0]), AtomId::new(0), false).unwrap(),
        0
    );
}

#[test]
fn neighbor_mode_counts_hydrogen_isotopes_and_excludes_other_elements_and_dummy() {
    let input = topology(
        vec![
            atom_spec(Element::C),
            atom_spec(Element::H),
            atom_spec(Element::H).with_isotope(2),
            atom_spec(Element::H).with_isotope(3),
            atom_spec(Element::HE),
            atom_spec(Element::DUMMY),
        ],
        (1..=5)
            .map(|neighbor| BondSpec::new(AtomId::new(0), AtomId::new(neighbor), BondOrder::Single))
            .collect(),
    );
    assert_eq!(
        total_hydrogen_count(&input, &assignment(vec![0]), AtomId::new(0), true).unwrap(),
        3
    );
}

#[test]
fn neighbor_mode_is_independent_of_bond_order_and_dative_direction() {
    for (order, hydrogen_is_end) in [
        (BondOrder::Single, true),
        (BondOrder::Zero, true),
        (BondOrder::Aromatic, true),
        (BondOrder::Dative, true),
        (BondOrder::Dative, false),
        (BondOrder::DativeLeft, true),
        (BondOrder::DativeRight, false),
        (BondOrder::ThreeCenter, true),
        (BondOrder::Other, false),
    ] {
        let (atoms, bond, target) = if hydrogen_is_end {
            (
                vec![atom_spec(Element::C), atom_spec(Element::H)],
                BondSpec::new(AtomId::new(0), AtomId::new(1), order),
                AtomId::new(0),
            )
        } else {
            (
                vec![atom_spec(Element::H), atom_spec(Element::C)],
                BondSpec::new(AtomId::new(0), AtomId::new(1), order),
                AtomId::new(1),
            )
        };
        let input = topology(atoms, vec![bond]);
        assert_eq!(
            total_hydrogen_count(&input, &assignment(vec![0, 0]), target, true).unwrap(),
            1,
            "order={order:?}, hydrogen_is_end={hydrogen_is_end}"
        );
    }
}

#[test]
fn helium_neon_source_rows_have_zero_total_hydrogens() {
    let input = topology(
        vec![
            atom_spec(Element::HE).with_no_implicit(true),
            atom_spec(Element::NE).with_no_implicit(true),
            atom_spec(Element::HE)
                .with_no_implicit(true)
                .with_formal_charge(1),
            atom_spec(Element::NE)
                .with_no_implicit(true)
                .with_formal_charge(1),
        ],
        vec![],
    );
    for atom in 0..4 {
        assert_eq!(
            total_hydrogen_count(&input, &assignment(vec![]), AtomId::new(atom), false).unwrap(),
            0
        );
    }
}

#[test]
fn invalid_atom_and_topology_errors_preserve_exact_fields() {
    let valid = topology(vec![atom_spec(Element::C)], vec![]);
    assert_eq!(
        total_hydrogen_count(&valid, &assignment(vec![]), AtomId::new(3), false),
        Err(ValenceError::AtomOutOfRange {
            atom: AtomId::new(3),
            atom_count: 1,
        })
    );

    let invalid = TopologyBlock {
        atoms: vec![Atom::from_spec(AtomId::new(1), atom_spec(Element::C))],
        adjacency: AdjacencyList::from_topology(1, &[]),
        ..TopologyBlock::default()
    };
    assert_eq!(
        total_hydrogen_count(&invalid, &assignment(vec![0]), AtomId::new(0), false),
        Err(ValenceError::InvalidTopology {
            source: TopologyValidationError::AtomIdMismatch {
                position: 0,
                id: AtomId::new(1),
            },
        })
    );
}

#[test]
fn signed_source_accumulator_overflow_is_structured() {
    let input = topology(
        vec![atom_spec(Element::C).with_explicit_hydrogens(1)],
        vec![],
    );
    assert_eq!(
        total_hydrogen_count(&input, &assignment(vec![i32::MAX]), AtomId::new(0), false,),
        Err(ValenceError::HydrogenCountOverflow {
            atom: AtomId::new(0),
            explicit: 1,
            implicit: i32::MAX as u32,
            neighbor_hydrogens: 0,
        })
    );
}
