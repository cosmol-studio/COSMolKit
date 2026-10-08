use cosmolkit_core::{
    ValenceError, ValenceModel, ValencePhase, assign_explicit_valence_for_atom_from_parts,
    assign_valence_state_for_atom_from_parts, assign_valence_with_options_for_topology,
    calculate_explicit_valence_for_topology,
};
use cosmolkit_model::{Atom, AtomId, AtomSpec, TopologyBlock};
use cosmolkit_types::Element;

fn carbon(h: u8, no_implicit: bool) -> TopologyBlock {
    TopologyBlock::try_from_parts(
        vec![Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_explicit_hydrogens(h)
                .with_no_implicit(no_implicit),
        )],
        vec![],
        vec![],
        vec![],
    )
    .unwrap()
}

#[test]
fn atom_cache_update_uses_both_source_stored_values_in_same_order_as_whole_assignment() {
    // Atom.cpp updatePropertyCache calls calcExplicitValence then
    // calcImplicitValence, whose returned caches both have signed8 width.
    for (h, explicit, implicit) in [
        (0, 0, 4),
        (4, 4, 0),
        (127, 127, 0),
        (128, -128, -124),
        (129, -127, -125),
        (132, -124, -128),
        (133, -123, 127),
        (254, -2, 6),
        (255, -1, 5),
    ] {
        for no_implicit in [false, true] {
            let t = carbon(h, no_implicit);
            let before = t.clone();
            let expected = (explicit, if no_implicit { 0 } else { implicit });
            assert_eq!(
                calculate_explicit_valence_for_topology(&t, AtomId::new(0), false, false),
                Ok(i32::from(h))
            );
            assert_eq!(
                assign_explicit_valence_for_atom_from_parts(
                    &t.atoms,
                    &t.bonds,
                    &t.adjacency,
                    AtomId::new(0),
                    false
                ),
                Ok(explicit),
                "H={h}, noImplicit={no_implicit}"
            );
            assert_eq!(
                assign_valence_state_for_atom_from_parts(
                    &t.atoms,
                    &t.bonds,
                    &t.adjacency,
                    AtomId::new(0),
                    false
                ),
                Ok(expected)
            );
            let whole =
                assign_valence_with_options_for_topology(&t, ValenceModel::RdkitLike, false)
                    .unwrap();
            assert_eq!(
                (whole.explicit_valence[0], whole.implicit_hydrogens[0]),
                expected
            );
            assert_eq!(t, before);
        }
    }
}

#[test]
fn atom_cache_update_propagates_explicit_failure_before_implicit_no_implicit_success() {
    let t = carbon(5, true);
    let before = t.clone();
    assert!(
        matches!(assign_valence_state_for_atom_from_parts(&t.atoms,&t.bonds,&t.adjacency,AtomId::new(0),true),Err(ValenceError::InvalidValence { atom, atomic_number: 6, phase: ValencePhase::Explicit, calculated: Some(5), .. }) if atom==AtomId::new(0))
    );
    assert!(
        matches!(assign_valence_with_options_for_topology(&t,ValenceModel::RdkitLike,true),Err(ValenceError::InvalidValence { atom, atomic_number: 6, phase: ValencePhase::Explicit, calculated: Some(5), .. }) if atom==AtomId::new(0))
    );
    assert_eq!(t, before);
}
