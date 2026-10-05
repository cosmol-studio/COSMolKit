use cosmolkit_core::atom_metadata;
use cosmolkit_model::{Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, TopologyBlock};
use cosmolkit_types::{BondOrder, Element};
#[test]
fn metadata_retains_degree_explicit_and_implicit_hydrogen_distinction() {
    let atoms = vec![
        Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
        Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::O)),
        Atom::from_spec(AtomId::new(2), AtomSpec::new(Element::H)),
    ];
    let bonds = vec![
        Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        ),
        Bond::from_spec(
            BondId::new(1),
            BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Single),
        ),
    ];
    let topology = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap();
    let before = topology.clone();
    let rows = atom_metadata(&topology).unwrap();
    assert_eq!(rows.len(), 3);
    assert_eq!(
        (
            rows[0].degree,
            rows[0].explicit_valence,
            rows[0].implicit_hydrogens,
            rows[0].total_hydrogens,
            rows[0].total_valence
        ),
        (1, 1, 3, 3, 4)
    );
    assert_eq!(
        (
            rows[1].degree,
            rows[1].explicit_valence,
            rows[1].implicit_hydrogens,
            rows[1].total_hydrogens,
            rows[1].total_valence
        ),
        (2, 2, 0, 0, 2)
    );
    assert_eq!(topology, before);
}
#[test]
fn metadata_propagates_invalid_valence_and_empty_is_present_empty() {
    assert!(atom_metadata(&TopologyBlock::default()).unwrap().is_empty());
    let atom = Atom::from_spec(
        AtomId::new(0),
        AtomSpec::new(Element::C).with_explicit_hydrogens(5),
    );
    let topology = TopologyBlock::try_from_parts(vec![atom], vec![], vec![], vec![]).unwrap();
    assert!(atom_metadata(&topology).is_err());
}

#[test]
fn cached_metadata_retains_non_strict_cache_instead_of_recalculating_strict_valence() {
    use cosmolkit_core::{ValenceAssignment, atom_metadata_from_assignment};
    let input = TopologyBlock::try_from_parts(
        vec![Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C).with_explicit_hydrogens(5),
        )],
        vec![],
        vec![],
        vec![],
    )
    .unwrap();
    let cache = ValenceAssignment {
        explicit_valence: vec![5],
        implicit_hydrogens: vec![0],
    };
    assert!(atom_metadata(&input).is_err());
    let row = atom_metadata_from_assignment(&input, Some(&cache)).unwrap()[0];
    assert_eq!(
        (
            row.explicit_valence,
            row.implicit_hydrogens,
            row.total_hydrogens,
            row.total_valence
        ),
        (5, 0, 5, 5)
    );
}
#[test]
fn cached_metadata_no_implicit_bypasses_a_stale_or_missing_implicit_entry() {
    use cosmolkit_core::{ValenceAssignment, atom_metadata_from_assignment};
    let input = TopologyBlock::try_from_parts(
        vec![Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_no_implicit(true)
                .with_explicit_hydrogens(2),
        )],
        vec![],
        vec![],
        vec![],
    )
    .unwrap();
    for (explicit, implicit_hydrogens) in [(0, vec![4]), (2, vec![4]), (2, vec![-1]), (2, vec![])] {
        let cache = ValenceAssignment {
            explicit_valence: vec![explicit],
            implicit_hydrogens,
        };
        let before = cache.clone();
        let row = atom_metadata_from_assignment(&input, Some(&cache)).unwrap()[0];
        assert_eq!(
            (
                row.implicit_hydrogens,
                row.total_hydrogens,
                row.total_valence
            ),
            (0, 2, explicit)
        );
        assert_eq!(cache, before);
    }
}
#[test]
fn cached_metadata_preserves_explicit_and_implicit_precondition_errors() {
    use cosmolkit_core::{ValenceAssignment, ValenceError, atom_metadata_from_assignment};
    let input = TopologyBlock::try_from_parts(
        vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C))],
        vec![],
        vec![],
        vec![],
    )
    .unwrap();
    for cache in [
        None,
        Some(ValenceAssignment {
            explicit_valence: vec![-1],
            implicit_hydrogens: vec![4],
        }),
    ] {
        assert!(
            matches!(atom_metadata_from_assignment(&input,cache.as_ref()),Err(ValenceError::ExplicitValenceCacheNotInitialized{atom}) if atom==AtomId::new(0))
        );
    }
    let cache = ValenceAssignment {
        explicit_valence: vec![0],
        implicit_hydrogens: vec![-1],
    };
    assert!(
        matches!(atom_metadata_from_assignment(&input,Some(&cache)),Err(ValenceError::ImplicitValenceCacheNotInitialized{atom}) if atom==AtomId::new(0))
    );
    assert!(
        atom_metadata_from_assignment(&TopologyBlock::default(), None)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn cached_metadata_reproduces_source_int8_explicit_cache_boundary() {
    use cosmolkit_core::{
        ValenceError, ValenceParams, assign_valence, atom_metadata_from_assignment,
    };
    // Source-built RDKit: non-strict UpdatePropertyCache, noImplicit=true.
    // Atom.h int8 cache: 127 is readable; 128/200/255 violate GetValence's
    // narrowed-cache > -1 precondition. Keep the supported non-strict path.
    for hydrogens in [4_u8, 127, 128, 200, 255] {
        let topology = TopologyBlock::try_from_parts(
            vec![Atom::from_spec(
                AtomId::new(0),
                AtomSpec::new(Element::C)
                    .with_explicit_hydrogens(hydrogens)
                    .with_no_implicit(true),
            )],
            vec![],
            vec![],
            vec![],
        )
        .unwrap();
        let cache = assign_valence(
            &topology,
            &ValenceParams {
                strict: false,
                ..ValenceParams::default()
            },
        )
        .unwrap();
        let before_topology = topology.clone();
        let before_cache = cache.clone();
        let result = atom_metadata_from_assignment(&topology, Some(&cache));
        if hydrogens <= 127 {
            let row = result.unwrap()[0];
            assert_eq!(
                (
                    row.explicit_valence,
                    row.implicit_hydrogens,
                    row.total_hydrogens,
                    row.total_valence
                ),
                (
                    i32::from(hydrogens),
                    0,
                    i32::from(hydrogens),
                    i32::from(hydrogens)
                )
            );
        } else {
            assert!(
                matches!(result,
                    Err(ValenceError::ExplicitValenceCacheNotInitialized { atom })
                    if atom == AtomId::new(0)
                ),
                "explicit H {hydrogens}: {result:?}"
            );
        }
        assert_eq!(topology, before_topology);
        assert_eq!(cache, before_cache);
    }
}
