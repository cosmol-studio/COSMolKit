//! One private S0 observation; source literals frozen before CK execution.
use super::*;
use crate::{AtomId, BondId, BondOrder, BondStereo, Hybridization, SmilesParseParams};
use std::sync::Arc;

#[test]
fn drawing_line22_probe_constructor() {
    let result = Molecule::from_smiles_with_params(
        "O[C@](C)(Cl)[C@@](O)(Cl)C",
        &SmilesParseParams {
            sanitize: true,
            remove_hydrogens: true,
            allow_cxsmiles: true,
            strict_cxsmiles: true,
            parse_name: true,
            skip_cleanup: false,
            debug_parse: false,
            replacements: Default::default(),
        },
    );
    if let Err(error) = &result {
        println!("S0_CONSTRUCTOR_ERROR {error:?} {error}");
    }
    let molecule = result.expect("S0 actual constructor result");
    let topology = molecule.topology_arc_runtime();
    let coordinates = molecule.coordinates_arc_runtime();
    let properties = molecule.properties_arc_runtime();
    let cache = molecule.derived_cache_arc_runtime();
    let baseline = (
        topology.as_ref().clone(),
        coordinates.as_ref().clone(),
        properties.as_ref().clone(),
        cache.as_ref().clone(),
    );
    let valence_identity = cache.valence_assignment().map(|v| v as *const _);
    let ring_identity = cache.valid_ring_info().map(|r| r as *const _);
    let checkpoint = || {
        assert!(Arc::ptr_eq(&topology, &molecule.topology_arc_runtime()));
        assert!(Arc::ptr_eq(
            &coordinates,
            &molecule.coordinates_arc_runtime()
        ));
        assert!(Arc::ptr_eq(&properties, &molecule.properties_arc_runtime()));
        assert!(Arc::ptr_eq(&cache, &molecule.derived_cache_arc_runtime()));
        assert_eq!(molecule.topology(), &baseline.0);
        assert_eq!(molecule.coordinate_block_runtime(), &baseline.1);
        assert_eq!(molecule.properties(), &baseline.2);
        assert_eq!(molecule.derived_cache_runtime(), &baseline.3);
        assert_eq!(
            molecule
                .derived_cache_runtime()
                .valence_assignment()
                .map(|v| v as *const _),
            valence_identity
        );
        assert_eq!(
            molecule
                .derived_cache_runtime()
                .valid_ring_info()
                .map(|r| r as *const _),
            ring_identity
        );
    };
    checkpoint();
    let input = molecule.drawing_input();
    checkpoint();
    println!("S0_TOPOLOGY {:#?}", input.topology);
    println!("S0_COORDINATES {:#?}", input.coordinates);
    println!("S0_PROPERTIES {:#?}", input.properties);
    println!(
        "S0_VALENCE {:#?} S0_RINGS {:#?}",
        input.valence, input.rings
    );
    for atom in &input.topology.atoms {
        println!(
            "S0_ATOM {} z={} isotope={:?} charge={} explicit_h={} no_implicit={} radicals={} aromatic={} hybrid={} chiral={} permutation={:?} map={:?} unknown={} props={:?} computed={:?} query=absent_in_concrete_Atom",
            atom.id().index(),
            atom.atomic_number(),
            atom.isotope(),
            atom.formal_charge(),
            atom.explicit_hydrogens(),
            atom.no_implicit(),
            atom.radical_electrons(),
            atom.is_aromatic(),
            atom.hybridization().rdkit_code(),
            atom.chiral_tag().rdkit_code(),
            atom.chiral_permutation(),
            atom.atom_map(),
            atom.unknown_stereo(),
            atom.props(),
            atom.computed_prop_names()
        );
    }
    for bond in &input.topology.bonds {
        println!(
            "S0_BOND {} endpoints=({}, {}) order={:?} aromatic={} conjugated={} direction={} stereo={:?} stereo_atoms={:?} unknown={} props={:?} computed={:?} query=absent_in_concrete_Bond",
            bond.id().index(),
            bond.begin().index(),
            bond.end().index(),
            bond.order(),
            bond.is_aromatic(),
            bond.is_conjugated(),
            bond.direction().rdkit_code(),
            bond.stereo(),
            bond.stereo_atoms(),
            bond.unknown_stereo(),
            bond.props(),
            bond.computed_prop_names()
        );
    }
    let adjacency = (0..input.topology.atoms.len())
        .map(|i| {
            input
                .topology
                .adjacency
                .neighbors_of(i)
                .iter()
                .map(|n| (n.bond.index(), n.atom_index))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    println!("S0_ADJACENCY {adjacency:?}");
    if let Some(rings) = input.rings {
        println!(
            "S0_RING_ROWS atoms={} bonds={} initialized={} quality={:?} atom_members={:?} bond_members={:?}",
            rings.atom_row_count(),
            rings.bond_row_count(),
            rings.is_initialized(),
            rings.find_type(),
            (0..8)
                .map(|i| rings.atom_members(AtomId::new(i)))
                .collect::<Vec<_>>(),
            (0..7)
                .map(|i| rings.bond_members(BondId::new(i)))
                .collect::<Vec<_>>()
        );
    }
    // All actual consumed-state/property observations above precede equality.
    assert_eq!(input.topology.atoms.len(), 8);
    assert_eq!(input.topology.bonds.len(), 7);
    for (i, atom) in input.topology.atoms.iter().enumerate() {
        assert_eq!(atom.id().index(), i);
        assert_eq!(atom.atomic_number(), [8, 6, 6, 17, 6, 8, 17, 6][i]);
        assert_eq!(
            atom.no_implicit(),
            [false, true, false, false, true, false, false, false][i]
        );
        assert_eq!(atom.chiral_tag().rdkit_code(), [0, 2, 0, 0, 1, 0, 0, 0][i]);
        assert_eq!(atom.hybridization(), Hybridization::Sp3);
        assert_eq!(atom.isotope(), None);
        assert_eq!(atom.atom_map(), None);
        assert_eq!(atom.chiral_permutation(), None);
        assert_eq!(
            (
                atom.formal_charge(),
                atom.explicit_hydrogens(),
                atom.radical_electrons()
            ),
            (0, 0, 0)
        );
        assert!(!atom.is_aromatic() && !atom.unknown_stereo());
    }
    for (i, bond) in input.topology.bonds.iter().enumerate() {
        assert_eq!(bond.id().index(), i);
        assert_eq!(
            (bond.begin().index(), bond.end().index()),
            [(0, 1), (1, 2), (1, 3), (1, 4), (4, 5), (4, 6), (4, 7)][i]
        );
        assert_eq!(bond.order(), BondOrder::Single);
        assert_eq!(bond.direction().rdkit_code(), 0);
        assert_eq!(bond.stereo(), BondStereo::None);
        assert_eq!(bond.stereo_atoms(), None);
        assert!(!bond.is_aromatic() && !bond.is_conjugated() && !bond.unknown_stereo());
    }
    assert_eq!(
        adjacency,
        vec![
            vec![(0, 1)],
            vec![(0, 0), (1, 2), (2, 3), (3, 4)],
            vec![(1, 1)],
            vec![(2, 1)],
            vec![(3, 1), (4, 5), (5, 6), (6, 7)],
            vec![(4, 4)],
            vec![(5, 4)],
            vec![(6, 4)]
        ]
    );
    assert!(input.topology.substance_groups.is_empty() && input.topology.stereo_groups.is_empty());
    assert!(
        input.coordinates.conformers_2d.is_empty() && input.coordinates.conformers_3d.is_empty()
    );
    let valence = input
        .valence
        .expect("actual borrowed source-shaped valence");
    assert_eq!(valence.explicit_valence, [1, 4, 1, 1, 4, 1, 1, 1]);
    assert_eq!(valence.implicit_hydrogens, [1, 0, 3, 0, 0, 1, 0, 3]);
    let rings = input
        .rings
        .expect("actual borrowed initialized-empty ring carrier");
    assert!(rings.is_initialized());
    assert_eq!((rings.atom_row_count(), rings.bond_row_count()), (8, 7));
    assert!(rings.atom_rings().is_empty() && rings.bond_rings().is_empty());
    for i in 0..8 {
        assert_eq!(rings.num_atom_rings(AtomId::new(i)), 0);
        assert!(rings.atom_members(AtomId::new(i)).is_empty());
    }
    for i in 0..7 {
        assert_eq!(rings.num_bond_rings(BondId::new(i)), 0);
        assert!(rings.bond_members(BondId::new(i)).is_empty());
    }
}
