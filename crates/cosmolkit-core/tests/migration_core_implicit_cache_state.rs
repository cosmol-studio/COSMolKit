use cosmolkit_core::{calculate_implicit_valence_for_topology, implicit_valence_for_atom};
use cosmolkit_model::{Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, TopologyBlock};
use cosmolkit_types::{BondOrder, Element};

fn isolated(element: Element, hydrogens: u8, radicals: u8, aromatic: bool) -> TopologyBlock {
    TopologyBlock::try_from_parts(
        vec![Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(element)
                .with_explicit_hydrogens(hydrogens)
                .with_radical_electrons(radicals)
                .with_aromatic(aromatic),
        )],
        vec![],
        vec![],
        vec![],
    )
    .unwrap()
}

#[test]
fn implicit_kernel_recomputed_local_does_not_replace_original_cache_in_radical_sum() {
    // Atom.cpp::calculateImplicitValence recomputes the local auto int8 value
    // for its hydrogen special case; explicitPlusRadV separately reads the
    // original d_explicitValence field. For C cache=-1, that sum is -1+rad.
    for (hydrogens, radicals, expected) in [(0, 0, 5), (2, 0, 5), (3, 0, 5), (0, 1, 4), (3, 2, 3)] {
        for aromatic in [false, true] {
            let topology = isolated(Element::C, hydrogens, radicals, aromatic);
            assert_eq!(
                calculate_implicit_valence_for_topology(
                    &topology,
                    AtomId::new(0),
                    -1,
                    false,
                    false
                ),
                Ok(expected)
            );
        }
    }
    // The higher-level optional query computes the explicit cache first,
    // matching Atom::calcImplicitValence prelude; it does not expose -1 cache.
    assert_eq!(
        implicit_valence_for_atom(
            &isolated(Element::C, 0, 0, false),
            AtomId::new(0),
            None,
            true
        ),
        Ok(4)
    );
}

#[test]
fn implicit_kernel_recomputed_local_retains_signed8_width_before_hydrogen_guard() {
    for (hydrogens, expected) in [(0, 1), (128, 2), (255, 2)] {
        assert_eq!(
            calculate_implicit_valence_for_topology(
                &isolated(Element::H, hydrogens, 0, false),
                AtomId::new(0),
                -1,
                false,
                false
            ),
            Ok(expected)
        );
    }
    // int result256 narrows to auto int8 local0 before the neutralH guard.
    let atoms = (0..=256)
        .map(|i| {
            Atom::from_spec(
                AtomId::new(i),
                AtomSpec::new(if i == 0 { Element::H } else { Element::F }),
            )
        })
        .collect();
    let bonds = (1..=256)
        .map(|i| {
            Bond::from_spec(
                BondId::new(i - 1),
                BondSpec::new(AtomId::new(0), AtomId::new(i), BondOrder::Single),
            )
        })
        .collect();
    let topology = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap();
    assert_eq!(
        calculate_implicit_valence_for_topology(&topology, AtomId::new(0), -1, false, false),
        Ok(1)
    );
}
