use cosmolkit_core::{
    RingSearchParams, ValenceModel, assign_legacy_stereochemistry_with_assignments,
    assign_valence_with_options_for_topology, find_sssr,
};
use cosmolkit_model::{Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, TopologyBlock};
use cosmolkit_types::{BondOrder, Element};

#[test]
fn clean_legacy_stereo_upgrades_ring_state_even_without_chiral_atoms() {
    // Pinned Chirality.cpp findChiralAtomSpecialCases upgrades before any chiral guard.
    let atoms = (0..8)
        .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
        .collect();
    let mut bonds = Vec::new();
    for i in 0..8 {
        for bit in [1, 2, 4] {
            let j = i ^ bit;
            if i < j {
                bonds.push(Bond::from_spec(
                    BondId::new(bonds.len()),
                    BondSpec::new(AtomId::new(i), AtomId::new(j), BondOrder::Single),
                ));
            }
        }
    }
    let topology = TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap();
    let rings = find_sssr(&topology, &RingSearchParams::default()).unwrap();
    assert!(!rings.is_symm_sssr());
    let before = (topology.clone(), rings.clone());
    let valence =
        assign_valence_with_options_for_topology(&topology, ValenceModel::RdkitLike, false)
            .unwrap();
    let unchanged = assign_legacy_stereochemistry_with_assignments(
        topology.clone(),
        &valence,
        &rings,
        false,
        false,
    )
    .unwrap();
    assert!(unchanged.ring_update.is_none());
    let cleaned = assign_legacy_stereochemistry_with_assignments(
        topology.clone(),
        &valence,
        &rings,
        true,
        false,
    )
    .unwrap();
    let updated = cleaned
        .ring_update
        .expect("source cleanup must retain SymmSSSR effect");
    assert!(updated.is_symm_sssr());
    for i in 0..8 {
        assert_eq!(updated.num_atom_rings(AtomId::new(i)), 3);
    }
    assert_eq!((topology, rings), before);
    assert_eq!(cleaned.topology, before.0);
}
