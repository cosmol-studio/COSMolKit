use cosmolkit_core::{
    ValenceError, ValenceModel, assign_radicals, assign_valence_for_topology,
    assign_valence_with_options_for_topology, atom_has_valence_violation_for_topology,
    bond_valence_contrib, calculate_explicit_valence_for_topology,
    calculate_implicit_valence_for_topology,
};
use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, TopologyBlock,
};
use cosmolkit_types::{BondOrder, Element};

fn detached_topology(atom_specs: Vec<AtomSpec>, bond_specs: Vec<BondSpec>) -> TopologyBlock {
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

fn atom_spec(atomic_number: u8) -> AtomSpec {
    AtomSpec::new(Element::from_atomic_number(atomic_number).expect("test element"))
}

fn aromatic_topology(
    atomic_numbers: &[u8],
    edges: &[(usize, usize)],
    explicit_hydrogens: &[(usize, u8)],
) -> TopologyBlock {
    let mut specs = atomic_numbers
        .iter()
        .map(|&atomic_number| atom_spec(atomic_number).with_aromatic(true))
        .collect::<Vec<_>>();
    for &(index, hydrogens) in explicit_hydrogens {
        specs[index] = specs[index]
            .clone()
            .with_explicit_hydrogens(hydrogens)
            .with_no_implicit(true);
    }
    let bonds = edges
        .iter()
        .map(|&(begin, end)| {
            BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Aromatic)
                .with_aromatic(true)
        })
        .collect();
    detached_topology(specs, bonds)
}

#[test]
fn detached_valence_matches_rdkit_2026_03_1_for_benzene() {
    let topology = aromatic_topology(
        &[6; 6],
        &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)],
        &[],
    );
    let assignment = assign_valence_for_topology(&topology, ValenceModel::RdkitLike).unwrap();

    assert_eq!(assignment.explicit_valence, vec![3; 6]);
    assert_eq!(assignment.implicit_hydrogens, vec![1; 6]);
}

#[test]
fn detached_valence_matches_rdkit_2026_03_1_for_fused_aromatic_junctions() {
    let topology = aromatic_topology(
        &[6; 10],
        &[
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 0),
            (7, 2),
        ],
        &[],
    );
    let assignment = assign_valence_for_topology(&topology, ValenceModel::RdkitLike).unwrap();

    assert_eq!(
        assignment.explicit_valence,
        vec![3, 3, 4, 3, 3, 3, 3, 4, 3, 3]
    );
    assert_eq!(
        assignment.implicit_hydrogens,
        vec![1, 1, 0, 1, 1, 1, 1, 0, 1, 1]
    );
}

#[test]
fn detached_valence_matches_rdkit_2026_03_1_for_aromatic_nh() {
    let topology = aromatic_topology(
        &[6, 6, 6, 7, 6],
        &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)],
        &[(3, 1)],
    );
    let assignment = assign_valence_for_topology(&topology, ValenceModel::RdkitLike).unwrap();

    assert_eq!(assignment.explicit_valence, vec![3, 3, 3, 3, 3]);
    assert_eq!(assignment.implicit_hydrogens, vec![1, 1, 1, 0, 1]);
}

#[test]
fn detached_valence_counts_dative_order_only_at_the_acceptor() {
    let topology = detached_topology(
        vec![atom_spec(7), atom_spec(8)],
        vec![BondSpec::new(
            AtomId::new(0),
            AtomId::new(1),
            BondOrder::Dative,
        )],
    );
    let assignment = assign_valence_for_topology(&topology, ValenceModel::RdkitLike).unwrap();

    assert_eq!(assignment.explicit_valence, vec![0, 1]);
    assert_eq!(assignment.implicit_hydrogens, vec![3, 1]);
    assert_eq!(
        bond_valence_contrib(&topology.bonds[0], AtomId::new(0)).unwrap(),
        0.0
    );
    assert_eq!(
        bond_valence_contrib(&topology.bonds[0], AtomId::new(1)).unwrap(),
        1.0
    );
    assert_eq!(
        bond_valence_contrib(&topology.bonds[0], AtomId::new(2)).unwrap(),
        0.0
    );
}

#[test]
fn detached_valence_rejects_unimplemented_directional_dative_types_like_rdkit() {
    for order in [BondOrder::DativeLeft, BondOrder::DativeRight] {
        let topology = detached_topology(
            vec![atom_spec(7), atom_spec(8)],
            vec![BondSpec::new(AtomId::new(0), AtomId::new(1), order)],
        );
        assert_eq!(
            assign_valence_for_topology(&topology, ValenceModel::RdkitLike)
                .unwrap_err()
                .to_string(),
            "Bad bond type"
        );
    }
}

#[test]
fn detached_valence_matches_rdkit_hypervalent_anion_branches() {
    for atomic_number in [15, 16, 33, 34] {
        let mut atoms = vec![atom_spec(atomic_number).with_formal_charge(-1)];
        atoms.extend((0..5).map(|_| atom_spec(9)));
        let bonds = (1..=5)
            .map(|index| BondSpec::new(AtomId::new(0), AtomId::new(index), BondOrder::Single))
            .collect();
        let topology = detached_topology(atoms, bonds);
        let assignment = assign_valence_for_topology(&topology, ValenceModel::RdkitLike).unwrap();

        assert_eq!(assignment.explicit_valence, vec![5, 1, 1, 1, 1, 1]);
        let central_implicit = if matches!(atomic_number, 15 | 33) {
            1
        } else {
            0
        };
        assert_eq!(
            assignment.implicit_hydrogens,
            vec![central_implicit, 0, 0, 0, 0, 0]
        );
    }
}

#[test]
fn detached_valence_preserves_strict_non_strict_and_check_it_semantics() {
    let topology = detached_topology(
        vec![atom_spec(6), atom_spec(8), atom_spec(8), atom_spec(8)],
        vec![
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Double),
            BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Double),
            BondSpec::new(AtomId::new(0), AtomId::new(3), BondOrder::Double),
        ],
    );

    let error = assign_valence_for_topology(&topology, ValenceModel::RdkitLike).unwrap_err();
    assert_eq!(
        error.to_string(),
        "Explicit valence for atom # 0 C, 6, is greater than permitted"
    );
    let assignment =
        assign_valence_with_options_for_topology(&topology, ValenceModel::RdkitLike, false)
            .unwrap();
    assert_eq!(assignment.explicit_valence, vec![6, 2, 2, 2]);
    assert_eq!(assignment.implicit_hydrogens, vec![0, 0, 0, 0]);
    assert_eq!(
        calculate_explicit_valence_for_topology(&topology, AtomId::new(0), false, true).unwrap(),
        -1
    );
    assert!(atom_has_valence_violation_for_topology(&topology, AtomId::new(0)).unwrap());
}

#[test]
fn detached_implicit_valence_preserves_hydrogen_charge_error_semantics() {
    let topology = detached_topology(vec![atom_spec(1).with_formal_charge(2)], vec![]);

    assert_eq!(
        calculate_implicit_valence_for_topology(&topology, AtomId::new(0), 0, true, false)
            .unwrap_err()
            .to_string(),
        "Unreasonable formal charge on atom # 0."
    );
    assert_eq!(
        calculate_implicit_valence_for_topology(&topology, AtomId::new(0), 0, false, false)
            .unwrap(),
        0
    );
    assert_eq!(
        calculate_implicit_valence_for_topology(&topology, AtomId::new(0), 0, false, true).unwrap(),
        -1
    );
}

#[test]
fn detached_valence_violation_checks_effective_periodic_table_row() {
    let topology = detached_topology(vec![atom_spec(9).with_formal_charge(-2)], vec![]);
    assert!(atom_has_valence_violation_for_topology(&topology, AtomId::new(0)).unwrap());

    let overflow = detached_topology(vec![atom_spec(1).with_formal_charge(-120)], vec![]);
    assert!(atom_has_valence_violation_for_topology(&overflow, AtomId::new(0)).unwrap());
}

#[test]
fn detached_radical_assignment_matches_rdkit_main_group_and_hypervalent_branches() {
    let main_group = detached_topology(
        vec![
            atom_spec(6).with_no_implicit(true),
            atom_spec(6)
                .with_no_implicit(true)
                .with_explicit_hydrogens(3),
            atom_spec(7).with_no_implicit(true),
            atom_spec(7)
                .with_no_implicit(true)
                .with_explicit_hydrogens(4)
                .with_formal_charge(1),
        ],
        vec![],
    );
    assert_eq!(
        assign_radicals(&main_group).unwrap().radical_electrons,
        vec![4, 1, 3, 0]
    );

    let mut atoms = vec![atom_spec(16).with_no_implicit(true).with_formal_charge(-1)];
    atoms.extend((0..4).map(|_| atom_spec(9)));
    let topology = detached_topology(
        atoms,
        (1..=4)
            .map(|index| BondSpec::new(AtomId::new(0), AtomId::new(index), BondOrder::Single))
            .collect(),
    );
    assert_eq!(assign_radicals(&topology).unwrap().radical_electrons[0], 1);
}

#[test]
fn detached_errors_remain_structured() {
    let topology = detached_topology(
        vec![atom_spec(6), atom_spec(8), atom_spec(8), atom_spec(8)],
        vec![
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Double),
            BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Double),
            BondSpec::new(AtomId::new(0), AtomId::new(3), BondOrder::Double),
        ],
    );
    assert!(matches!(
        assign_valence_for_topology(&topology, ValenceModel::RdkitLike),
        Err(ValenceError::InvalidValence {
            atom,
            atomic_number: 6,
            formal_charge: 0,
            ..
        }) if atom == AtomId::new(0)
    ));
}
