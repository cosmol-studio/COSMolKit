use cosmolkit_core::calculate_explicit_valence_for_topology;
use cosmolkit_model::{Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, TopologyBlock};
use cosmolkit_types::{BondOrder, Element};

// RDKit2026.03.1 first canonical numeric rows with the exact single [-1]
// valence list, directly selected from atomic_data.cpp (not a metal heuristic).
const UNSPECIFIED_ROWS: &[u8] = &[
    0, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 57, 58, 59,
    60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 89, 90,
    91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
    112, 113, 114, 115, 116, 117, 118,
];

#[test]
fn aromatic_unspecified_rows_keep_zero_and_explicit_hydrogens_under_unsigned_default() {
    // Atom.cpp347..453 uses unsigned int dv; -1 becomes UINT_MAX. Thus
    // accum(0..3) > dv is false and std::round(accum+0.1) returns accum.
    for &atomic_number in UNSPECIFIED_ROWS {
        for charge in [-8, 0, 8] {
            for hydrogens in 0..=3 {
                let topology = TopologyBlock::try_from_parts(
                    vec![Atom::from_spec(
                        AtomId::new(0),
                        AtomSpec::new(Element::from_atomic_number(atomic_number).unwrap())
                            .with_aromatic(true)
                            .with_formal_charge(charge)
                            .with_explicit_hydrogens(hydrogens),
                    )],
                    vec![],
                    vec![],
                    vec![],
                )
                .unwrap();
                for (strict, check_it) in
                    [(false, false), (false, true), (true, false), (true, true)]
                {
                    assert_eq!(
                        calculate_explicit_valence_for_topology(
                            &topology,
                            AtomId::new(0),
                            strict,
                            check_it
                        ),
                        Ok(i32::from(hydrogens)),
                        "Z={atomic_number}, charge={charge}, H={hydrogens}, strict={strict}, check={check_it}"
                    );
                }
            }
        }
    }
}

#[test]
fn unspecified_rows_preserve_bond_only_aromaticity_and_half_order_rounding() {
    for &atomic_number in UNSPECIFIED_ROWS {
        for (order, expected) in [
            (BondOrder::Zero, 0),
            (BondOrder::Single, 1),
            (BondOrder::OneAndHalf, 2),
            (BondOrder::Aromatic, 2),
        ] {
            let topology = TopologyBlock::try_from_parts(
                vec![
                    Atom::from_spec(
                        AtomId::new(0),
                        AtomSpec::new(Element::from_atomic_number(atomic_number).unwrap()),
                    ),
                    Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::F)),
                ],
                vec![Bond::from_spec(
                    BondId::new(0),
                    BondSpec::new(AtomId::new(0), AtomId::new(1), order).with_aromatic(true),
                )],
                vec![],
                vec![],
            )
            .unwrap();
            for (strict, check_it) in [(false, false), (false, true), (true, false), (true, true)] {
                assert_eq!(
                    calculate_explicit_valence_for_topology(
                        &topology,
                        AtomId::new(0),
                        strict,
                        check_it
                    ),
                    Ok(expected),
                    "Z={atomic_number}, bond={order:?}, strict={strict}, check={check_it}"
                );
            }
        }
    }
}
