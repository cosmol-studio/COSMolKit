use cosmolkit_core::{
    RingFindType, RingInfo, RingSearchParams, fast_find_rings, set_double_bond_neighbor_directions,
    symmetrized_sssr,
};
use cosmolkit_model::{Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, TopologyBlock};
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
fn source_guard_returns_ring_effect_even_without_candidate_double_bonds() {
    // Chirality.cpp prepares SymmSSSR before the first candidate loop and
    // before the no-bonds-in-play return. This includes empty/acyclic graphs.
    let mut cube_edges = vec![];
    for i in 0..8 {
        for bit in [1, 2, 4] {
            let j = i ^ bit;
            if i < j {
                cube_edges.push((i, j));
            }
        }
    }
    for (topology, expected_count) in [
        (graph(0, &[]), 0),
        (graph(3, &[(0, 1), (1, 2)]), 0),
        (graph(8, &cube_edges), 6),
    ] {
        let before = topology.clone();
        let fast = fast_find_rings(&topology).unwrap();
        let old_fast = fast.clone();
        let result = set_double_bond_neighbor_directions(topology.clone(), &fast, None).unwrap();
        let updated = result
            .ring_update
            .expect("native required ring preparation must be returned");
        assert!(updated.is_symm_sssr());
        assert_eq!(updated.num_rings(), expected_count);
        assert_eq!(result.topology, before);
        assert!(!result.needs_detect_bond_stereo);
        assert_eq!(fast, old_fast);
        if expected_count == 6 {
            for i in 0..8 {
                assert_eq!(updated.num_atom_rings(AtomId::new(i)), 3);
            }
        }
        let already_prepared = symmetrized_sssr(&topology, &RingSearchParams::default()).unwrap();
        let result =
            set_double_bond_neighbor_directions(topology.clone(), &already_prepared, None).unwrap();
        assert!(result.ring_update.is_none());
        assert_eq!(result.topology, before);
        // Initialized-empty Symm is sufficient even on an empty topology.
        if topology.atoms.is_empty() {
            let sufficient = RingInfo::new(RingFindType::SymmSssr, 0, 0);
            assert!(
                set_double_bond_neighbor_directions(topology, &sufficient, None)
                    .unwrap()
                    .ring_update
                    .is_none()
            );
        }
    }
}
