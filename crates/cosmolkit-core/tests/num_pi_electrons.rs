//! Source-backed numPi tests with all13 conditions retained from the approved proposal.
use cosmolkit_core::{ValenceError, num_pi_electrons_for_topology};
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, ChiralTag, Element, Hybridization,
    MoleculeProperties, PropertyValue, TopologyBlock,
};
fn star(
    z: u8,
    degree: usize,
    hybrid: Hybridization,
    aromatic: bool,
    hs: u8,
    tag: ChiralTag,
    order: BondOrder,
) -> TopologyBlock {
    let mut atoms = vec![Atom::from_spec(
        AtomId::new(0),
        AtomSpec::new(Element::from_atomic_number(z).unwrap())
            .with_hybridization(hybrid)
            .with_aromatic(aromatic)
            .with_explicit_hydrogens(hs)
            .with_chiral_tag(tag),
    )];
    for i in 1..=degree {
        atoms.push(Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)));
    }
    let bonds = (1..=degree)
        .map(|i| {
            Bond::from_spec(
                BondId::new(i - 1),
                BondSpec::new(AtomId::new(0), AtomId::new(i), order),
            )
        })
        .collect();
    TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
}
#[test]
fn num_pi_aromatic_precedes_every_hybrid_cache_and_bond_branch() {
    for hybrid in [
        Hybridization::Unspecified,
        Hybridization::S,
        Hybridization::Sp,
        Hybridization::Sp2,
        Hybridization::Sp3,
        Hybridization::Sp2d,
        Hybridization::Sp3d,
        Hybridization::Sp3d2,
        Hybridization::Other,
    ] {
        for cache in [None, Some(-128), Some(-1), Some(0), Some(127)] {
            let t = star(
                6,
                1,
                hybrid,
                true,
                255,
                ChiralTag::Unspecified,
                BondOrder::Other,
            );
            assert_eq!(
                num_pi_electrons_for_topology(&t, AtomId::new(0), cache),
                Ok(1)
            );
        }
    }
}
#[test]
fn num_pi_sp3_bypasses_cache_and_bad_bond() {
    for cache in [None, Some(-128), Some(-1), Some(0), Some(127)] {
        let t = star(
            6,
            1,
            Hybridization::Sp3,
            false,
            255,
            ChiralTag::Unspecified,
            BondOrder::Other,
        );
        assert_eq!(
            num_pi_electrons_for_topology(&t, AtomId::new(0), cache),
            Ok(0)
        );
    }
}
#[test]
fn num_pi_nonsp3_all_hybridizations_get_exact_cache() {
    for hybrid in [
        Hybridization::Unspecified,
        Hybridization::S,
        Hybridization::Sp,
        Hybridization::Sp2,
        Hybridization::Sp2d,
        Hybridization::Sp3d,
        Hybridization::Sp3d2,
        Hybridization::Other,
    ] {
        let t = star(
            6,
            1,
            hybrid,
            false,
            0,
            ChiralTag::Unspecified,
            BondOrder::Single,
        );
        assert_eq!(
            num_pi_electrons_for_topology(&t, AtomId::new(0), Some(3)),
            Ok(2)
        );
    }
}
#[test]
fn num_pi_every_defined_signed8_cache_literal() {
    let t = star(
        6,
        0,
        Hybridization::Unspecified,
        false,
        0,
        ChiralTag::Unspecified,
        BondOrder::Single,
    );
    for &(cache, expected) in CACHE {
        assert_eq!(
            num_pi_electrons_for_topology(&t, AtomId::new(0), Some(cache)),
            Ok(expected)
        );
    }
    assert_eq!(CACHE.len(), 128);
}
#[test]
fn num_pi_missing_and_all_negative_cache_fail_before_bond_read() {
    let t = star(
        6,
        1,
        Hybridization::Unspecified,
        false,
        0,
        ChiralTag::Unspecified,
        BondOrder::Other,
    );
    for cache in std::iter::once(None).chain((-128i8..=-1).map(Some)) {
        let err = num_pi_electrons_for_topology(&t, AtomId::new(0), cache).unwrap_err();
        assert!(matches!(
            &err,
            ValenceError::PiElectronExplicitValenceCacheNotInitialized { .. }
        ));
        assert_eq!(
            err.to_string(),
            "getValence(ValenceType::EXPLICIT) called without call to calcExplicitValence()"
        );
    }
}
#[test]
fn num_pi_all_bond_types_and_both_endpoints_use_nonzero_contributions() {
    let mut observations = 0;
    for &(order, first, second) in BONDS {
        let t = star(
            6,
            1,
            Hybridization::Unspecified,
            false,
            0,
            ChiralTag::Unspecified,
            order,
        );
        for (id, expected) in [(AtomId::new(0), first), (AtomId::new(1), second)] {
            match expected {
                Some(expected) => {
                    assert_eq!(num_pi_electrons_for_topology(&t, id, Some(1)), Ok(expected))
                }
                None => assert!(matches!(
                    num_pi_electrons_for_topology(&t, id, Some(1)),
                    Err(ValenceError::BadBondType { .. })
                )),
            }
            observations += 1;
        }
    }
    assert_eq!(observations, 44);
}
#[test]
fn num_pi_explicit_hydrogens_physical_count_and_invariant() {
    for hs in [0, 1, 2, 126, 127] {
        let t = star(
            6,
            0,
            Hybridization::Unspecified,
            false,
            hs,
            ChiralTag::Unspecified,
            BondOrder::Single,
        );
        assert_eq!(
            num_pi_electrons_for_topology(&t, AtomId::new(0), Some(127)),
            Ok(127 - u32::from(hs))
        );
    }
    let t = star(
        6,
        0,
        Hybridization::Unspecified,
        false,
        255,
        ChiralTag::Unspecified,
        BondOrder::Single,
    );
    let err = num_pi_electrons_for_topology(&t, AtomId::new(0), Some(127)).unwrap_err();
    assert!(matches!(
        err,
        ValenceError::PiElectronInvariant {
            explicit_valence: 127,
            physical_bonds: 255,
            ..
        }
    ));
    assert_eq!(err.to_string(), "explicit valence exceeds atom degree");
}
#[test]
fn num_pi_zero_valence_nonzero_physical_bond_keeps_source_error_text() {
    let t = star(
        6,
        1,
        Hybridization::Unspecified,
        false,
        0,
        ChiralTag::Unspecified,
        BondOrder::Single,
    );
    let err = num_pi_electrons_for_topology(&t, AtomId::new(0), Some(0)).unwrap_err();
    assert!(matches!(
        err,
        ValenceError::PiElectronInvariant {
            explicit_valence: 0,
            physical_bonds: 1,
            ..
        }
    ));
    assert_eq!(err.to_string(), "explicit valence exceeds atom degree");
}
#[test]
fn num_pi_dative_and_dativeone_donor_do_not_count_as_physical() {
    for order in [BondOrder::Dative, BondOrder::DativeOne] {
        let t = star(
            6,
            1,
            Hybridization::Unspecified,
            false,
            0,
            ChiralTag::Unspecified,
            order,
        );
        assert_eq!(
            num_pi_electrons_for_topology(&t, AtomId::new(0), Some(0)),
            Ok(0)
        );
        assert!(matches!(
            num_pi_electrons_for_topology(&t, AtomId::new(1), Some(0)),
            Err(ValenceError::PiElectronInvariant { .. })
        ));
    }
}
#[test]
fn num_pi_degree_is_not_the_physical_bond_count() {
    for order in [
        BondOrder::Zero,
        BondOrder::Ionic,
        BondOrder::Hydrogen,
        BondOrder::Unspecified,
    ] {
        let t = star(
            6,
            1,
            Hybridization::Unspecified,
            false,
            0,
            ChiralTag::Unspecified,
            order,
        );
        assert_eq!(
            num_pi_electrons_for_topology(&t, AtomId::new(0), Some(2)),
            Ok(2)
        );
    }
}
#[test]
fn num_pi_selected_atom_out_of_range_remains_structural() {
    assert!(matches!(
        num_pi_electrons_for_topology(&TopologyBlock::default(), AtomId::new(0), None),
        Err(ValenceError::AtomOutOfRange { .. })
    ));
}
#[test]
fn num_pi_malformed_adjacency_reference_is_not_defaulted() {
    let mut t = star(
        6,
        1,
        Hybridization::Unspecified,
        false,
        0,
        ChiralTag::Unspecified,
        BondOrder::Single,
    );
    t.bonds.clear();
    assert!(matches!(
        num_pi_electrons_for_topology(&t, AtomId::new(0), Some(1)),
        Err(ValenceError::AdjacencyBondOutOfRange { .. })
    ));
}
#[test]
fn num_pi_input_blocks_unchanged_on_success_and_failure() {
    let t = star(
        6,
        1,
        Hybridization::Unspecified,
        false,
        0,
        ChiralTag::Unspecified,
        BondOrder::Single,
    );
    let before = t.clone();
    assert_eq!(
        num_pi_electrons_for_topology(&t, AtomId::new(0), Some(3)),
        Ok(2)
    );
    assert!(num_pi_electrons_for_topology(&t, AtomId::new(0), Some(0)).is_err());
    assert_eq!(t, before);
}
const CACHE: &[(i8, u32)] = &[
    (0_i8, 0_u32),
    (1_i8, 1_u32),
    (2_i8, 2_u32),
    (3_i8, 3_u32),
    (4_i8, 4_u32),
    (5_i8, 5_u32),
    (6_i8, 6_u32),
    (7_i8, 7_u32),
    (8_i8, 8_u32),
    (9_i8, 9_u32),
    (10_i8, 10_u32),
    (11_i8, 11_u32),
    (12_i8, 12_u32),
    (13_i8, 13_u32),
    (14_i8, 14_u32),
    (15_i8, 15_u32),
    (16_i8, 16_u32),
    (17_i8, 17_u32),
    (18_i8, 18_u32),
    (19_i8, 19_u32),
    (20_i8, 20_u32),
    (21_i8, 21_u32),
    (22_i8, 22_u32),
    (23_i8, 23_u32),
    (24_i8, 24_u32),
    (25_i8, 25_u32),
    (26_i8, 26_u32),
    (27_i8, 27_u32),
    (28_i8, 28_u32),
    (29_i8, 29_u32),
    (30_i8, 30_u32),
    (31_i8, 31_u32),
    (32_i8, 32_u32),
    (33_i8, 33_u32),
    (34_i8, 34_u32),
    (35_i8, 35_u32),
    (36_i8, 36_u32),
    (37_i8, 37_u32),
    (38_i8, 38_u32),
    (39_i8, 39_u32),
    (40_i8, 40_u32),
    (41_i8, 41_u32),
    (42_i8, 42_u32),
    (43_i8, 43_u32),
    (44_i8, 44_u32),
    (45_i8, 45_u32),
    (46_i8, 46_u32),
    (47_i8, 47_u32),
    (48_i8, 48_u32),
    (49_i8, 49_u32),
    (50_i8, 50_u32),
    (51_i8, 51_u32),
    (52_i8, 52_u32),
    (53_i8, 53_u32),
    (54_i8, 54_u32),
    (55_i8, 55_u32),
    (56_i8, 56_u32),
    (57_i8, 57_u32),
    (58_i8, 58_u32),
    (59_i8, 59_u32),
    (60_i8, 60_u32),
    (61_i8, 61_u32),
    (62_i8, 62_u32),
    (63_i8, 63_u32),
    (64_i8, 64_u32),
    (65_i8, 65_u32),
    (66_i8, 66_u32),
    (67_i8, 67_u32),
    (68_i8, 68_u32),
    (69_i8, 69_u32),
    (70_i8, 70_u32),
    (71_i8, 71_u32),
    (72_i8, 72_u32),
    (73_i8, 73_u32),
    (74_i8, 74_u32),
    (75_i8, 75_u32),
    (76_i8, 76_u32),
    (77_i8, 77_u32),
    (78_i8, 78_u32),
    (79_i8, 79_u32),
    (80_i8, 80_u32),
    (81_i8, 81_u32),
    (82_i8, 82_u32),
    (83_i8, 83_u32),
    (84_i8, 84_u32),
    (85_i8, 85_u32),
    (86_i8, 86_u32),
    (87_i8, 87_u32),
    (88_i8, 88_u32),
    (89_i8, 89_u32),
    (90_i8, 90_u32),
    (91_i8, 91_u32),
    (92_i8, 92_u32),
    (93_i8, 93_u32),
    (94_i8, 94_u32),
    (95_i8, 95_u32),
    (96_i8, 96_u32),
    (97_i8, 97_u32),
    (98_i8, 98_u32),
    (99_i8, 99_u32),
    (100_i8, 100_u32),
    (101_i8, 101_u32),
    (102_i8, 102_u32),
    (103_i8, 103_u32),
    (104_i8, 104_u32),
    (105_i8, 105_u32),
    (106_i8, 106_u32),
    (107_i8, 107_u32),
    (108_i8, 108_u32),
    (109_i8, 109_u32),
    (110_i8, 110_u32),
    (111_i8, 111_u32),
    (112_i8, 112_u32),
    (113_i8, 113_u32),
    (114_i8, 114_u32),
    (115_i8, 115_u32),
    (116_i8, 116_u32),
    (117_i8, 117_u32),
    (118_i8, 118_u32),
    (119_i8, 119_u32),
    (120_i8, 120_u32),
    (121_i8, 121_u32),
    (122_i8, 122_u32),
    (123_i8, 123_u32),
    (124_i8, 124_u32),
    (125_i8, 125_u32),
    (126_i8, 126_u32),
    (127_i8, 127_u32),
];
const BONDS: &[(BondOrder, Option<u32>, Option<u32>)] = &[
    (BondOrder::Unspecified, Some(1_u32), Some(1_u32)),
    (BondOrder::Ionic, Some(1_u32), Some(1_u32)),
    (BondOrder::Zero, Some(1_u32), Some(1_u32)),
    (BondOrder::Hydrogen, Some(1_u32), Some(1_u32)),
    (BondOrder::Single, Some(0_u32), Some(0_u32)),
    (BondOrder::Double, Some(0_u32), Some(0_u32)),
    (BondOrder::Triple, Some(0_u32), Some(0_u32)),
    (BondOrder::Quadruple, Some(0_u32), Some(0_u32)),
    (BondOrder::Quintuple, Some(0_u32), Some(0_u32)),
    (BondOrder::Hextuple, Some(0_u32), Some(0_u32)),
    (BondOrder::OneAndHalf, Some(0_u32), Some(0_u32)),
    (BondOrder::TwoAndHalf, Some(0_u32), Some(0_u32)),
    (BondOrder::ThreeAndHalf, Some(0_u32), Some(0_u32)),
    (BondOrder::FourAndHalf, Some(0_u32), Some(0_u32)),
    (BondOrder::FiveAndHalf, Some(0_u32), Some(0_u32)),
    (BondOrder::Aromatic, Some(0_u32), Some(0_u32)),
    (BondOrder::Dative, Some(1_u32), Some(0_u32)),
    (BondOrder::DativeOne, Some(1_u32), Some(0_u32)),
    (BondOrder::DativeLeft, None, None),
    (BondOrder::DativeRight, None, None),
    (BondOrder::ThreeCenter, None, None),
    (BondOrder::Other, None, None),
];
