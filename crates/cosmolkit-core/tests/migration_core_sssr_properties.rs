use cosmolkit_core::{
    RingFindingError, RingSearchParams, symmetrized_sssr, symmetrized_sssr_with_properties,
};
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, MoleculeProperties, MoleculePropertyError,
    PropertyValue, TopologyBlock,
};
use cosmolkit_types::{BondOrder, Element};

fn graph(n: usize, edges: &[(usize, usize)]) -> TopologyBlock {
    let atoms = (0..n)
        .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
        .collect();
    let bonds = edges
        .iter()
        .enumerate()
        .map(|(i, &(a, b))| {
            Bond::from_spec(
                BondId::new(i),
                BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
            )
        })
        .collect();
    TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
}

#[test]
fn initial_clear_retains_order_and_removes_only_source_extra_ring_membership() {
    // FindRings.cpp clears extraRings before component processing, even for
    // an empty graph; RDProps.h removes just its first computed-list entry.
    for topology in [graph(0, &[]), graph(3, &[(0, 1), (1, 2), (2, 0)])] {
        let mut props = MoleculeProperties::default();
        props.set_prop("head", 1_i32).unwrap();
        props.set_computed_prop("extraRings", "discarded").unwrap();
        props.set_computed_prop("kept", 2_i32).unwrap();
        props.set_prop("tail", 3_i32).unwrap();
        let (result, props) =
            symmetrized_sssr_with_properties(&topology, props, &RingSearchParams::default())
                .unwrap();
        assert_eq!(
            result,
            symmetrized_sssr(&topology, &RingSearchParams::default()).unwrap()
        );
        assert!(props.prop("extraRings").is_none());
        assert_eq!(
            props.prop("__computedProps"),
            Some(&PropertyValue::StringVector(vec!["kept".into()]))
        );
        assert_eq!(
            props
                .ordered_props()
                .map(|(k, _)| k.clone())
                .collect::<Vec<_>>(),
            vec![
                "head".into(),
                "__computedProps".into(),
                "kept".into(),
                "tail".into()
            ]
        );
    }
}

#[test]
fn malformed_computed_list_fails_even_when_extra_rings_is_absent() {
    // RDProps::clearProp reads the computed list before clearing the key.
    for existing_extra in [false, true] {
        let mut props = MoleculeProperties::default();
        props.set_prop("__computedProps", 3_i32).unwrap();
        if existing_extra {
            props.set_prop("extraRings", "retained-on-error").unwrap();
        }
        let before = props.clone();
        let result = symmetrized_sssr_with_properties(
            &graph(0, &[]),
            props.clone(),
            &RingSearchParams::default(),
        );
        assert!(matches!(
            result,
            Err(RingFindingError::MoleculeProperty(
                MoleculePropertyError::ComputedListKind(_)
            ))
        ));
        assert_eq!(props, before);
    }
}

#[test]
fn symmetrization_consumes_transient_cache_across_disconnected_cubes() {
    // Each cube has five SSSR rows and six equally sized face rings; native
    // removeExtraRings appends extras across fragments, then SymmSSSR erases
    // its computed extraRings property after consuming both fragments. The
    // first computed-list creation survives as an empty vector.
    let mut edges = vec![];
    for offset in [0, 8] {
        for i in 0..8 {
            for bit in [1, 2, 4] {
                let j = i ^ bit;
                if i < j {
                    edges.push((offset + i, offset + j));
                }
            }
        }
    }
    let topology = graph(16, &edges);
    let mut props = MoleculeProperties::default();
    props.set_prop("retained", "value").unwrap();
    let mut expected = props.clone();
    expected
        .set_prop("__computedProps", PropertyValue::StringVector(vec![]))
        .unwrap();
    for _ in 0..2 {
        let (info, prepared_properties) =
            symmetrized_sssr_with_properties(&topology, props, &RingSearchParams::default())
                .unwrap();
        props = prepared_properties;
        assert!(info.is_symm_sssr());
        assert_eq!(info.num_rings(), 12);
        for i in 0..16 {
            assert_eq!(info.num_atom_rings(AtomId::new(i)), 3);
        }
        assert!(props.prop("extraRings").is_none());
        assert_eq!(props, expected);
    }
}
