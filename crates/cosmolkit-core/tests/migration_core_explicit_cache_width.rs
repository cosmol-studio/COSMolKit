use cosmolkit_core::{
    ValenceModel, assign_valence_with_options_for_topology, calculate_explicit_valence_for_topology,
};
use cosmolkit_model::{Atom, AtomId, AtomSpec, TopologyBlock};
use cosmolkit_types::Element;

#[test]
fn explicit_assignment_returns_source_signed_byte_cache_at_both_narrowing_boundaries() {
    // Atom.h declares d_explicitValence std::int8_t. Atom.cpp::calcExplicitValence
    // stores into that field before returning it; calculateExplicitValence
    // itself returns int. Literal expectations preserve the separate boundaries.
    for element in [Element::DUMMY, Element::FE] {
        for (hydrogens, stored) in [
            (0, 0),
            (127, 127),
            (128, -128),
            (129, -127),
            (254, -2),
            (255, -1),
        ] {
            for no_implicit in [false, true] {
                let topology = TopologyBlock::try_from_parts(
                    vec![Atom::from_spec(
                        AtomId::new(0),
                        AtomSpec::new(element)
                            .with_explicit_hydrogens(hydrogens)
                            .with_no_implicit(no_implicit),
                    )],
                    vec![],
                    vec![],
                    vec![],
                )
                .unwrap();
                for strict in [false, true] {
                    assert_eq!(
                        calculate_explicit_valence_for_topology(
                            &topology,
                            AtomId::new(0),
                            strict,
                            false
                        ),
                        Ok(i32::from(hydrogens))
                    );
                    let assigned = assign_valence_with_options_for_topology(
                        &topology,
                        ValenceModel::RdkitLike,
                        strict,
                    )
                    .unwrap();
                    assert_eq!(
                        assigned.explicit_valence,
                        vec![stored],
                        "element={element:?}, H={hydrogens}, noImplicit={no_implicit}, strict={strict}"
                    );
                    assert_eq!(assigned.implicit_hydrogens, vec![0]);
                }
            }
        }
    }
}
