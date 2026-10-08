use cosmolkit_core::{
    ValenceModel, assign_implicit_valence_for_atom_from_parts_with_explicit_valence,
    assign_valence_with_options_for_topology, calculate_implicit_valence_for_topology,
};
use cosmolkit_model::{Atom, AtomId, AtomSpec, TopologyBlock};
use cosmolkit_types::Element;

fn carbon(hydrogens: u8, no_implicit: bool) -> TopologyBlock {
    TopologyBlock::try_from_parts(
        vec![Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_explicit_hydrogens(hydrogens)
                .with_no_implicit(no_implicit),
        )],
        vec![],
        vec![],
        vec![],
    )
    .unwrap()
}

#[test]
fn implicit_assignment_returns_source_signed_byte_cache_separately_from_int_kernel() {
    // Atom.cpp calcImplicitValence first checks/recomputes the explicit cache,
    // then stores calculateImplicitValence into Atom.h's std::int8_t field.
    // Literal expected int and stored values preserve both source boundaries.
    for (h, explicit, kernel, stored) in [
        (0, 0, 4, 4),
        (4, 4, 0, 0),
        (127, 127, 0, 0),
        (128, -128, 132, -124),
        (129, -127, 131, -125),
        (130, -126, 130, -126),
        (131, -125, 129, -127),
        (132, -124, 128, -128),
        (133, -123, 127, 127),
        (254, -2, 6, 6),
        (255, -1, 5, 5),
    ] {
        for no_implicit in [false, true] {
            let topology = carbon(h, no_implicit);
            let before = topology.clone();
            let assigned =
                assign_valence_with_options_for_topology(&topology, ValenceModel::RdkitLike, false)
                    .unwrap();
            assert_eq!(assigned.explicit_valence, vec![explicit]);
            assert_eq!(
                calculate_implicit_valence_for_topology(
                    &topology,
                    AtomId::new(0),
                    explicit,
                    false,
                    false
                ),
                Ok(if no_implicit { 0 } else { kernel })
            );
            assert_eq!(
                assigned.implicit_hydrogens,
                vec![if no_implicit { 0 } else { stored }],
                "H={h}, noImplicit={no_implicit}"
            );
            assert_eq!(topology, before);
        }
    }
}

#[test]
fn implicit_cache_wrapper_recomputes_explicit_sentinel_before_no_implicit_guard() {
    // Native wrapper performs calcExplicitValence(strict) before calling the
    // kernel, even when the kernel would return zero for noImplicit.
    let topology = carbon(5, true);
    let before = topology.clone();
    assert_eq!(
        calculate_implicit_valence_for_topology(&topology, AtomId::new(0), -1, true, false),
        Ok(0)
    );
    assert!(
        assign_implicit_valence_for_atom_from_parts_with_explicit_valence(
            &topology.atoms,
            &topology.bonds,
            &topology.adjacency,
            AtomId::new(0),
            -1,
            true
        )
        .is_err()
    );
    assert_eq!(
        assign_implicit_valence_for_atom_from_parts_with_explicit_valence(
            &topology.atoms,
            &topology.bonds,
            &topology.adjacency,
            AtomId::new(0),
            0,
            true
        ),
        Ok(0)
    );
    assert_eq!(topology, before);
}
