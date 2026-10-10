use super::*;
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, BondStereo, Element, Hybridization,
    PropertyValue,
};

fn acid_and_methane() -> TopologyBlock {
    let atoms = [
        Element::O,
        Element::C,
        Element::O,
        Element::C,
        Element::C,
        Element::C,
        Element::O,
        Element::O,
        Element::C,
    ]
    .into_iter()
    .enumerate()
    .map(|(i, element)| {
        Atom::from_spec(
            AtomId::new(i),
            AtomSpec::new(element)
                .with_hybridization(if i == 8 {
                    Hybridization::Sp3
                } else {
                    Hybridization::Sp2
                })
                .with_computed_prop(
                    "_CIPRank",
                    PropertyValue::Int([4, 2, 3, 1, 1, 2, 4, 3, 0][i]),
                )
                .unwrap(),
        )
    })
    .collect();
    let bonds = [
        (0, 1, BondOrder::Double),
        (1, 2, BondOrder::Single),
        (1, 3, BondOrder::Single),
        (3, 4, BondOrder::Double),
        (4, 5, BondOrder::Single),
        (5, 6, BondOrder::Double),
        (5, 7, BondOrder::Single),
    ]
    .into_iter()
    .enumerate()
    .map(|(i, (a, b, order))| {
        let mut spec = BondSpec::new(AtomId::new(a), AtomId::new(b), order);
        if i == 3 {
            spec = spec
                .with_stereo(BondStereo::Z)
                .with_stereo_atoms(AtomId::new(1), AtomId::new(5));
        }
        Bond::from_spec(BondId::new(i), spec)
    })
    .collect();
    TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
}

fn acid_and_linked_rings() -> TopologyBlock {
    let zs = [
        Element::O,
        Element::C,
        Element::O,
        Element::C,
        Element::C,
        Element::C,
        Element::O,
        Element::O,
        Element::C,
        Element::C,
        Element::N,
        Element::N,
        Element::N,
        Element::C,
        Element::C,
        Element::C,
        Element::N,
        Element::C,
        Element::C,
        Element::C,
    ];
    let ranks = [
        13, 7, 12, 2, 2, 7, 13, 12, 1, 5, 10, 11, 9, 3, 0, 4, 8, 4, 0, 6,
    ];
    let atoms = zs
        .into_iter()
        .enumerate()
        .map(|(i, z)| {
            Atom::from_spec(
                AtomId::new(i),
                AtomSpec::new(z)
                    .with_hybridization(Hybridization::Sp2)
                    .with_aromatic(i >= 8 && i != 12)
                    .with_computed_prop("_CIPRank", PropertyValue::Int(ranks[i]))
                    .unwrap(),
            )
        })
        .collect();
    let edges = [
        (0, 1, BondOrder::Double),
        (1, 2, BondOrder::Single),
        (1, 3, BondOrder::Single),
        (3, 4, BondOrder::Double),
        (4, 5, BondOrder::Single),
        (5, 6, BondOrder::Double),
        (5, 7, BondOrder::Single),
        (8, 9, BondOrder::Aromatic),
        (9, 10, BondOrder::Aromatic),
        (10, 11, BondOrder::Aromatic),
        (11, 12, BondOrder::Single),
        (12, 13, BondOrder::Single),
        (13, 14, BondOrder::Aromatic),
        (14, 15, BondOrder::Aromatic),
        (15, 16, BondOrder::Aromatic),
        (16, 17, BondOrder::Aromatic),
        (17, 18, BondOrder::Aromatic),
        (11, 19, BondOrder::Aromatic),
        (19, 8, BondOrder::Aromatic),
        (18, 13, BondOrder::Aromatic),
    ];
    let bonds = edges
        .into_iter()
        .enumerate()
        .map(|(i, (a, b, order))| {
            let mut spec = BondSpec::new(AtomId::new(a), AtomId::new(b), order)
                .with_aromatic(order == BondOrder::Aromatic);
            if i == 3 {
                spec = spec
                    .with_stereo(BondStereo::E)
                    .with_stereo_atoms(AtomId::new(1), AtomId::new(5));
            }
            Bond::from_spec(BondId::new(i), spec)
        })
        .collect();
    TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
}

#[test]
fn merging_an_earlier_seed_preserves_the_selected_fragment_position() {
    let topology = acid_and_linked_rings();
    let rings = symmetrized_sssr(&topology, &RingSearchParams::default()).unwrap();
    let mut fragments = compute_initial_coordinates(&topology, &rings, None, false).unwrap();
    // Pinned RDKit computeInitialCoords / expandEfrag: seed order [5-ring,
    // 6-ring, E-bond], then expand the 6-ring and erase the preceding 5-ring.
    assert_eq!(fragments.len(), 2);
    assert_eq!(
        fragments[0].atoms.keys().copied().collect::<Vec<_>>(),
        (8..20).collect::<Vec<_>>()
    );
    assert_eq!(
        fragments[1].atoms.keys().copied().collect::<Vec<_>>(),
        (0..8).collect::<Vec<_>>()
    );
    for f in &mut fragments {
        f.remove_collisions_bond_and_spiro_flip().unwrap();
    }
    for f in &mut fragments {
        f.remove_collisions_open_angles().unwrap();
        f.remove_collisions_shorten_bonds().unwrap();
    }
    orient_and_shift_fragments(&mut fragments, false, None);
    // Native bit patterns at the packing boundary. Internal embedding matched
    // before the repair; only fragment order and consequent translation differ.
    assert_eq!(
        fragments[0].atoms[&8].loc.map(f64::to_bits),
        [13830796075482922640, 13841274293248062788]
    );
    assert_eq!(
        fragments[1].atoms[&0].loc.map(f64::to_bits),
        [4617034042984890368, 13836404803183458688]
    );
}

#[test]
fn unmerged_stereo_seed_keeps_its_original_position() {
    let topology = acid_and_methane();
    let rings = symmetrized_sssr(&topology, &RingSearchParams::default()).unwrap();
    let mut fragments = compute_initial_coordinates(&topology, &rings, None, false).unwrap();
    assert_eq!(fragments.len(), 2);
    assert!(fragments[0].atoms.contains_key(&0));
    assert!(fragments[1].atoms.contains_key(&8));
    for f in &mut fragments {
        f.remove_collisions_bond_and_spiro_flip().unwrap();
    }
    for f in &mut fragments {
        f.remove_collisions_open_angles().unwrap();
        f.remove_collisions_shorten_bonds().unwrap();
    }
    orient_and_shift_fragments(&mut fragments, false, None);
    assert_eq!(fragments[0].atoms[&3].loc, [0.0, 0.0]);
}
