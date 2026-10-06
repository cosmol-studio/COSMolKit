use cosmolkit_core::{
    RadicalDiagnostic, RadicalError, ValenceError, assign_radicals, required_valence_list,
};
use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, TopologyBlock,
    TopologyValidationError,
};
use cosmolkit_types::{BondOrder, Element};

fn atom_spec(atomic_number: u8) -> AtomSpec {
    AtomSpec::new(Element::from_atomic_number(atomic_number).expect("test element"))
}

fn assigned_atom(atomic_number: u8) -> AtomSpec {
    atom_spec(atomic_number).with_no_implicit(true)
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

fn isolated(spec: AtomSpec) -> TopologyBlock {
    topology(vec![spec], vec![])
}

fn star(central: AtomSpec, neighbor_atomic_number: u8, count: usize) -> TopologyBlock {
    let mut atoms = vec![central];
    atoms.extend((0..count).map(|_| atom_spec(neighbor_atomic_number)));
    let bonds = (1..=count)
        .map(|neighbor| BondSpec::new(AtomId::new(0), AtomId::new(neighbor), BondOrder::Single))
        .collect();
    topology(atoms, bonds)
}

#[test]
fn preserves_stable_row_order_input_and_skipped_atom_radicals() {
    let input = topology(
        vec![
            atom_spec(6).with_radical_electrons(2),
            assigned_atom(6).with_explicit_hydrogens(3),
            atom_spec(0)
                .with_no_implicit(true)
                .with_radical_electrons(5),
            assigned_atom(7),
        ],
        vec![],
    );
    let before = input.clone();

    let assignment = assign_radicals(&input).unwrap();

    assert_eq!(assignment.radical_electrons, vec![2, 1, 5, 3]);
    assert_eq!(assignment.radical_electrons.len(), input.atoms.len());
    assert!(assignment.diagnostics.is_empty());
    assert_eq!(input, before);
}

#[test]
fn skipped_atoms_do_not_interpret_incident_unsupported_bond_orders() {
    for skipped in [
        atom_spec(6).with_radical_electrons(4),
        atom_spec(0)
            .with_no_implicit(true)
            .with_radical_electrons(6),
    ] {
        let input = topology(
            vec![skipped, atom_spec(8).with_radical_electrons(3)],
            vec![BondSpec::new(
                AtomId::new(0),
                AtomId::new(1),
                BondOrder::ThreeCenter,
            )],
        );
        assert_eq!(
            assign_radicals(&input).unwrap().radical_electrons,
            vec![input.atoms[0].radical_electrons(), 3]
        );
    }
}

#[test]
fn incident_contributions_preserve_dative_direction_zero_orders_and_half_order_rounding() {
    for (order, target_is_end, explicit_hydrogens, expected) in [
        (BondOrder::Dative, false, 0, 4),
        (BondOrder::Dative, true, 0, 3),
        (BondOrder::DativeOne, false, 0, 4),
        (BondOrder::DativeOne, true, 0, 3),
        (BondOrder::Zero, false, 0, 4),
        (BondOrder::Hydrogen, false, 0, 4),
        (BondOrder::OneAndHalf, false, 1, 2),
        (BondOrder::TwoAndHalf, false, 0, 2),
    ] {
        let (atoms, bond, target) = if target_is_end {
            (
                vec![atom_spec(8), assigned_atom(6)],
                BondSpec::new(AtomId::new(0), AtomId::new(1), order),
                1,
            )
        } else {
            (
                vec![
                    assigned_atom(6).with_explicit_hydrogens(explicit_hydrogens),
                    atom_spec(8),
                ],
                BondSpec::new(AtomId::new(0), AtomId::new(1), order),
                0,
            )
        };
        let input = topology(atoms, vec![bond]);
        assert_eq!(
            assign_radicals(&input).unwrap().radical_electrons[target],
            expected,
            "order={order:?}, target_is_end={target_is_end}"
        );
    }
}

#[test]
fn rejects_bad_bond_order_with_exact_nested_error_fields() {
    for order in [
        BondOrder::DativeLeft,
        BondOrder::DativeRight,
        BondOrder::ThreeCenter,
        BondOrder::Other,
    ] {
        let input = topology(
            vec![assigned_atom(6), atom_spec(8)],
            vec![BondSpec::new(AtomId::new(0), AtomId::new(1), order)],
        );
        assert_eq!(
            assign_radicals(&input),
            Err(RadicalError::Valence(ValenceError::BadBondType {
                bond: Some(BondId::new(0)),
                order,
            }))
        );
    }
}

#[test]
fn ordinary_branch_covers_two_electron_base_earlier_element_minimum_and_hypervalence() {
    for (atomic_number, expected) in [(1, 1), (2, 0), (5, 3)] {
        assert_eq!(
            assign_radicals(&isolated(assigned_atom(atomic_number)))
                .unwrap()
                .radical_electrons,
            vec![expected]
        );
    }

    let oxygen_with_three_bonds = star(assigned_atom(8), 1, 3);
    assert_eq!(
        assign_radicals(&oxygen_with_three_bonds)
            .unwrap()
            .radical_electrons[0],
        0,
        "a single allowed-valence row must not use -1 as a radical count"
    );

    let hypervalent_sulfur = star(assigned_atom(16).with_formal_charge(-1), 9, 4);
    assert_eq!(
        assign_radicals(&hypervalent_sulfur)
            .unwrap()
            .radical_electrons[0],
        1,
        "the first nonnegative allowed valence is selected in source order"
    );
}

#[test]
fn unlimited_valence_metals_use_bonded_zero_isolated_parity_and_ordered_diagnostics() {
    let input = topology(
        vec![
            assigned_atom(25).with_formal_charge(2),
            assigned_atom(25).with_formal_charge(1),
            assigned_atom(25),
            assigned_atom(25).with_formal_charge(-1),
            assigned_atom(26),
            assigned_atom(26).with_formal_charge(1),
            assigned_atom(59).with_formal_charge(4),
            assigned_atom(92).with_formal_charge(5),
            assigned_atom(26).with_formal_charge(30),
        ],
        vec![],
    );
    let assignment = assign_radicals(&input).unwrap();
    assert_eq!(
        assignment.radical_electrons,
        vec![1, 0, 1, 0, 0, 1, 0, 0, 0]
    );
    assert_eq!(
        assignment.diagnostics,
        vec![
            RadicalDiagnostic::UnusualChargeClamped {
                atom: AtomId::new(6),
                atomic_number: 59,
                formal_charge: 4,
            },
            RadicalDiagnostic::UnusualChargeClamped {
                atom: AtomId::new(7),
                atomic_number: 92,
                formal_charge: 5,
            },
            RadicalDiagnostic::UnusualChargeClamped {
                atom: AtomId::new(8),
                atomic_number: 26,
                formal_charge: 30,
            },
        ]
    );

    for bonded_metal in [
        assigned_atom(26).with_radical_electrons(7),
        assigned_atom(30)
            .with_formal_charge(1)
            .with_radical_electrons(7),
    ] {
        assert_eq!(
            assign_radicals(&star(bonded_metal, 6, 1))
                .unwrap()
                .radical_electrons[0],
            0
        );
    }
    assert_eq!(
        assign_radicals(&star(assigned_atom(41), 53, 5))
            .unwrap()
            .radical_electrons[0],
        0
    );
}

#[test]
fn rdkit_additional_oxidation_state_rows_assign_zero_radicals() {
    for (atomic_number, ligand_counts) in [
        (84, &[4, 6][..]),
        (54, &[4, 6][..]),
        (53, &[3, 5][..]),
        (85, &[3, 5][..]),
    ] {
        for &ligand_count in ligand_counts {
            assert_eq!(
                assign_radicals(&star(assigned_atom(atomic_number), 9, ligand_count))
                    .unwrap()
                    .radical_electrons[0],
                0,
                "atomic_number={atomic_number}, ligand_count={ligand_count}"
            );
        }
    }
}

#[test]
fn rdkit_helium_neon_and_main_group_charge_rows_match() {
    let input = topology(
        vec![
            assigned_atom(2),
            assigned_atom(10),
            assigned_atom(2).with_formal_charge(1),
            assigned_atom(10).with_formal_charge(1),
            assigned_atom(6),
            assigned_atom(6).with_formal_charge(1),
            assigned_atom(6).with_formal_charge(-1),
        ],
        vec![],
    );
    assert_eq!(
        assign_radicals(&input).unwrap().radical_electrons,
        vec![0, 0, 1, 1, 4, 3, 3]
    );
}

#[test]
fn rdkit_github_6370_overwrites_input_radicals_for_every_single_default_valence_row() {
    let mut exercised = 0;
    for atomic_number in 2..=118 {
        let valences = required_valence_list(atomic_number).unwrap();
        if valences.len() != 1 || valences[0] < 0 {
            continue;
        }
        let default_valence = valences[0] as u8;
        for explicit_hydrogens in 0..=default_valence {
            let input = isolated(
                assigned_atom(atomic_number)
                    .with_explicit_hydrogens(explicit_hydrogens)
                    .with_radical_electrons(7),
            );
            assert_eq!(
                assign_radicals(&input).unwrap().radical_electrons,
                vec![default_valence - explicit_hydrogens],
                "atomic_number={atomic_number}, explicit_hydrogens={explicit_hydrogens}"
            );
            exercised += 1;
        }
    }
    assert!(exercised > 0);

    let ammonium = isolated(
        assigned_atom(7)
            .with_explicit_hydrogens(4)
            .with_formal_charge(1)
            .with_radical_electrons(7),
    );
    assert_eq!(
        assign_radicals(&ammonium).unwrap().radical_electrons,
        vec![0]
    );
}

#[test]
fn malformed_topology_is_rejected_before_assignment() {
    let input = TopologyBlock {
        atoms: vec![Atom::from_spec(AtomId::new(1), assigned_atom(6))],
        adjacency: AdjacencyList::from_topology(1, &[]),
        ..TopologyBlock::default()
    };
    assert_eq!(
        assign_radicals(&input),
        Err(RadicalError::InvalidTopology {
            source: TopologyValidationError::AtomIdMismatch {
                position: 0,
                id: AtomId::new(1),
            },
        })
    );
}

#[test]
fn defensive_unreachable_errors_keep_atom_count_and_lookup_fields() {
    assert_eq!(
        RadicalError::PeriodicTableLookup {
            atomic_number: 119,
            field: "valences",
        },
        RadicalError::PeriodicTableLookup {
            atomic_number: 119,
            field: "valences",
        }
    );
    assert_eq!(
        RadicalError::RadicalCountOutOfRange {
            atom: AtomId::new(3),
            count: 256,
        },
        RadicalError::RadicalCountOutOfRange {
            atom: AtomId::new(3),
            count: 256,
        }
    );
}
