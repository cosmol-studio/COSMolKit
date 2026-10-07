//! Full proposed Chi recomputation arithmetic regressions, NOT installed or run.
//! Source pin351f8f378f8ad6bbd517980c38896e66bf907af8 (BSD).
//! Exports below require explicit ROOT decision; no force/cache semantics promised.
use cosmolkit_core::{RingFindType, RingInfo, ValenceAssignment, ValenceError};
use cosmolkit_descriptors::{
    CHI_0_N_VERSION, CHI_0_V_VERSION, CHI_1_N_VERSION, CHI_1_V_VERSION, CHI_2_N_VERSION,
    CHI_2_V_VERSION, CHI_3_N_VERSION, CHI_3_V_VERSION, CHI_4_N_VERSION, CHI_4_V_VERSION,
    CHI_N_N_VERSION, CHI_N_V_VERSION, ChiInput, DescriptorError, DescriptorResult, chi_0_n,
    chi_0_v, chi_1_n, chi_1_v, chi_2_n, chi_2_v, chi_3_n, chi_3_v, chi_4_n, chi_4_v, chi_n_n,
    chi_n_v,
};
use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, CoordinateBlock,
    Element, MoleculeProperties, TopologyBlock,
};
#[derive(Clone, Debug, PartialEq)]
struct Fixture {
    topology: TopologyBlock,
    coordinates: CoordinateBlock,
    properties: MoleculeProperties,
    valence: ValenceAssignment,
    rings: RingInfo,
}
impl Fixture {
    fn input(&self) -> ChiInput<'_> {
        ChiInput::new(&self.topology, &self.valence)
    }
}
// Explicit prepared-read states, not chemical preparation or degree-derived H.
fn graph(z: &[u8], hs: &[u8], edges: &[(usize, usize)]) -> Fixture {
    assert_eq!(z.len(), hs.len());
    let atoms = z
        .iter()
        .zip(hs)
        .enumerate()
        .map(|(i, (&n, &h))| {
            Atom::from_spec(
                AtomId::new(i),
                AtomSpec::new(Element::from_atomic_number(n).unwrap())
                    .with_no_implicit(true)
                    .with_explicit_hydrogens(h),
            )
        })
        .collect::<Vec<_>>();
    let bonds = edges
        .iter()
        .enumerate()
        .map(|(i, &(a, b))| {
            Bond::from_spec(
                BondId::new(i),
                BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
            )
        })
        .collect::<Vec<_>>();
    let adjacency = AdjacencyList::try_from_topology(z.len(), &bonds).unwrap();
    let topology = TopologyBlock {
        atoms,
        bonds,
        adjacency,
        ..Default::default()
    };
    topology.validate().unwrap();
    Fixture {
        topology,
        coordinates: CoordinateBlock::default(),
        properties: MoleculeProperties::default().with_name("Chi preservation sentinel"),
        valence: ValenceAssignment {
            explicit_valence: vec![0; z.len()],
            implicit_hydrogens: vec![0; z.len()],
        },
        rings: RingInfo::new(RingFindType::OtherOrUnknown, z.len(), edges.len()),
    }
}
type Entry = fn(&ChiInput<'_>) -> DescriptorResult<f64>;
const V: [Entry; 5] = [chi_0_v, chi_1_v, chi_2_v, chi_3_v, chi_4_v];
const N: [Entry; 5] = [chi_0_n, chi_1_n, chi_2_n, chi_3_n, chi_4_n];
fn close(a: f64, e: f64) {
    assert!(a.is_finite() && (a - e).abs() <= 1e-12, "{a:?} vs {e:?}");
    if e == 0.0 {
        assert_eq!(a.to_bits(), e.to_bits());
    }
}
fn check(f: &Fixture, v: [f64; 5], n: [f64; 5]) {
    let before = f.clone();
    for (family, expected) in [(V, v), (N, n)] {
        for (entry, e) in family.into_iter().zip(expected) {
            close(entry(&f.input()).unwrap(), e);
            assert_eq!(f, &before, "all five full borrowed values preserved");
        }
    }
    for k in 2..=4 {
        close(chi_n_v(&f.input(), k).unwrap(), v[k as usize]);
        close(chi_n_n(&f.input(), k).unwrap(), n[k as usize]);
    }
    assert_eq!(f, &before);
}
fn errors_preserve(f: &Fixture) {
    let before = f.clone();
    for entry in V.into_iter().chain(N) {
        assert!(entry(&f.input()).is_err());
        assert_eq!(f, &before);
    }
    assert!(chi_n_v(&f.input(), 2).is_err());
    assert!(chi_n_n(&f.input(), 2).is_err());
    assert_eq!(f, &before);
}
#[test]
fn chi_empty_all12() {
    for version in [
        CHI_0_V_VERSION,
        CHI_1_V_VERSION,
        CHI_2_V_VERSION,
        CHI_3_V_VERSION,
        CHI_4_V_VERSION,
        CHI_N_V_VERSION,
        CHI_0_N_VERSION,
        CHI_1_N_VERSION,
        CHI_2_N_VERSION,
        CHI_3_N_VERSION,
        CHI_4_N_VERSION,
        CHI_N_N_VERSION,
    ] {
        assert_eq!(version, "1.2.0");
    }
    check(&graph(&[], &[], &[]), [0.0; 5], [0.0; 5]);
}
#[test]
fn chi_all119_pinned_outer_and_heavy_denominators() {
    for &(z, hk, nv) in SINGLETONS {
        check(
            &graph(&[z], &[0], &[]),
            [hk, 0.0, 0.0, 0.0, 0.0],
            [nv, 0.0, 0.0, 0.0, 0.0],
        );
    }
}
#[test]
fn chi_first_row_zero_and_u32_subtraction() {
    check(&graph(&[6], &[4], &[]), [0.0; 5], [0.0; 5]);
    check(
        &graph(&[6], &[5], &[]),
        [0.000015258789064276357, 0.0, 0.0, 0.0, 0.0],
        [0.000015258789064276357, 0.0, 0.0, 0.0, 0.0],
    );
}
#[test]
fn chi_wildcard_and_h_are_hk_only_exclusions() {
    check(
        &graph(&[0, 1], &[1, 0], &[]),
        [0.0; 5],
        [1.0000152587890643, 0.0, 0.0, 0.0, 0.0],
    );
}
#[test]
fn chi_heavy_halogen_weight_split() {
    check(
        &graph(&[17], &[0], &[]),
        [1.1338934190276817, 0.0, 0.0, 0.0, 0.0],
        [0.3779644730092272, 0.0, 0.0, 0.0, 0.0],
    );
    check(
        &graph(&[35], &[0], &[]),
        [1.9639610121239313, 0.0, 0.0, 0.0, 0.0],
        [0.3779644730092272, 0.0, 0.0, 0.0, 0.0],
    );
}
#[test]
fn chi_chain_four_fixed_and_generic() {
    let f = graph(&[6; 4], &[3, 2, 2, 3], &[(0, 1), (1, 2), (2, 3)]);
    check(
        &f,
        [
            3.414213562373095,
            1.914213562373095,
            0.9999999999999998,
            0.4999999999999999,
            0.0,
        ],
        [
            3.414213562373095,
            1.914213562373095,
            0.9999999999999998,
            0.4999999999999999,
            0.0,
        ],
    );
}
#[test]
fn chi_branched_star_paths_are_not_subgraphs() {
    check(
        &graph(&[6; 5], &[0, 3, 3, 3, 3], &[(0, 1), (0, 2), (0, 3), (0, 4)]),
        [4.5, 2.0, 3.0, 0.0, 0.0],
        [4.5, 2.0, 3.0, 0.0, 0.0],
    );
}
#[test]
fn chi_triangle_closing_atom_counted_once() {
    check(
        &graph(&[6; 3], &[2; 3], &[(0, 1), (1, 2), (2, 0)]),
        [
            2.1213203435596424,
            1.4999999999999996,
            1.060660171779821,
            0.3535533905932737,
            0.0,
        ],
        [
            2.1213203435596424,
            1.4999999999999996,
            1.060660171779821,
            0.3535533905932737,
            0.0,
        ],
    );
}
#[test]
fn chi_square_bond_sets_not_atom_sets() {
    check(
        &graph(&[6; 4], &[2; 4], &[(0, 1), (1, 2), (2, 3), (3, 0)]),
        [
            2.82842712474619,
            1.9999999999999996,
            1.4142135623730947,
            0.9999999999999997,
            0.24999999999999992,
        ],
        [
            2.82842712474619,
            1.9999999999999996,
            1.4142135623730947,
            0.9999999999999997,
            0.24999999999999992,
        ],
    );
}
#[test]
fn chi_final_revisit_nonfirst_is_multiplied_again() {
    check(
        &graph(&[6; 4], &[3, 1, 2, 2], &[(0, 1), (1, 2), (2, 3), (3, 1)]),
        [
            2.991563831562721,
            1.8938468501173518,
            1.6825219847121646,
            0.8660254037844386,
            0.16666666666666666,
        ],
        [
            2.991563831562721,
            1.8938468501173518,
            1.6825219847121646,
            0.8660254037844386,
            0.16666666666666666,
        ],
    );
}
#[test]
fn chi_explicit_h_bonds_in_fixed_one_but_not_generic_one() {
    let f = graph(&[6, 1, 6], &[3, 0, 3], &[(0, 1), (1, 2)]);
    check(&f, [2.0, 0.0, 0.0, 0.0, 0.0], [3.0, 2.0, 0.0, 0.0, 0.0]);
    close(chi_n_v(&f.input(), 0).unwrap(), 3.0);
    close(chi_n_n(&f.input(), 0).unwrap(), 3.0);
    close(chi_n_v(&f.input(), 1).unwrap(), 0.0);
    close(chi_n_n(&f.input(), 1).unwrap(), 0.0);
}
#[test]
fn chi_generic_zero_one_large_and_uintmax_wrap() {
    let f = graph(&[6], &[0], &[]);
    close(chi_0_v(&f.input()).unwrap(), 0.5);
    close(chi_n_v(&f.input(), 0).unwrap(), 1.0);
    close(chi_n_n(&f.input(), 0).unwrap(), 1.0);
    for order in [1, 2, 3, 4, 5, 12, 31, u32::MAX] {
        close(chi_n_v(&f.input(), order).unwrap(), 0.0);
        close(chi_n_n(&f.input(), order).unwrap(), 0.0);
    }
}
#[test]
fn chi_bond_order_and_aromatic_flags_do_not_reweight() {
    for order in [
        BondOrder::Single,
        BondOrder::Double,
        BondOrder::Triple,
        BondOrder::Aromatic,
        BondOrder::Unspecified,
        BondOrder::Dative,
        BondOrder::Quadruple,
        BondOrder::Quintuple,
        BondOrder::Hextuple,
        BondOrder::OneAndHalf,
        BondOrder::TwoAndHalf,
        BondOrder::ThreeAndHalf,
        BondOrder::FourAndHalf,
        BondOrder::FiveAndHalf,
        BondOrder::Ionic,
        BondOrder::Hydrogen,
        BondOrder::ThreeCenter,
        BondOrder::DativeOne,
        BondOrder::DativeLeft,
        BondOrder::DativeRight,
        BondOrder::Other,
        BondOrder::Zero,
    ] {
        let mut f = graph(&[6; 2], &[3; 2], &[(0, 1)]);
        f.topology.bonds[0] = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), order).with_aromatic(true),
        );
        f.topology.adjacency = AdjacencyList::try_from_topology(2, &f.topology.bonds).unwrap();
        check(&f, [2.0, 1.0, 0.0, 0.0, 0.0], [2.0, 1.0, 0.0, 0.0, 0.0]);
    }
}
#[test]
fn chi_bond_insertion_and_atom_permutations() {
    for edges in [
        &[(0, 1), (1, 2), (2, 3)][..],
        &[(2, 3), (0, 1), (1, 2)][..],
        &[(3, 2), (2, 1), (1, 0)][..],
    ] {
        check(
            &graph(&[6; 4], &[3, 2, 2, 3], edges),
            [
                3.414213562373095,
                1.914213562373095,
                0.9999999999999998,
                0.4999999999999999,
                0.0,
            ],
            [
                3.414213562373095,
                1.914213562373095,
                0.9999999999999998,
                0.4999999999999999,
                0.0,
            ],
        );
    }
    // Heterogeneous triangle, all atom weights multiply once in a closed ring.
    let f = graph(&[6, 7, 8], &[3, 1, 2], &[(0, 1), (1, 2), (2, 0)]);
    let g = graph(&[8, 6, 7], &[2, 3, 1], &[(1, 2), (2, 0), (0, 1)]);
    check(
        &f,
        [2.0, 1.25, 0.75, 0.25, 0.0],
        [2.0, 1.25, 0.75, 0.25, 0.0],
    );
    check(
        &g,
        [2.0, 1.25, 0.75, 0.25, 0.0],
        [2.0, 1.25, 0.75, 0.25, 0.0],
    );

    // Exact pinned issue463 graph/H read rows and renumbering; supplemental source relational invariant.
    let z = [8, 6, 7, 6, 7, 6, 6, 16, 7, 6, 6, 6, 6, 6];
    let hs = [0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 2, 2, 3];
    let edges = [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 3),
        (1, 8),
        (8, 9),
        (9, 10),
        (10, 11),
        (11, 12),
        (12, 10),
        (9, 13),
    ];
    let order = [0, 11, 8, 7, 2, 4, 5, 13, 10, 12, 9, 3, 1, 6];
    let f = graph(&z, &hs, &edges);
    let mut old_to_new = [0usize; 14];
    for (new, &old) in order.iter().enumerate() {
        old_to_new[old] = new;
    }
    let zz = order.iter().map(|&i| z[i]).collect::<Vec<_>>();
    let hh = order.iter().map(|&i| hs[i]).collect::<Vec<_>>();
    let ee = edges
        .iter()
        .map(|&(a, b)| (old_to_new[a], old_to_new[b]))
        .collect::<Vec<_>>();
    let g = graph(&zz, &hh, &ee);
    let before_f = f.clone();
    let before_g = g.clone();
    for entry in V.into_iter().chain(N) {
        close(entry(&f.input()).unwrap(), entry(&g.input()).unwrap());
    }
    assert_eq!(f, before_f);
    assert_eq!(g, before_g);
}
#[test]
fn chi_negative_implicit_cache_and_hk_skipped_getters() {
    let mut f = graph(&[6], &[0], &[]);
    f.topology.atoms[0].set_no_implicit(false);
    f.valence.implicit_hydrogens[0] = -1;
    errors_preserve(&f);
    for entry in [chi_0_v as Entry, chi_0_n] {
        assert!(matches!(entry(&f.input()),
 Err(DescriptorError::Valence { source:ValenceError::ImplicitValenceCacheNotInitialized { atom },.. }) if atom==AtomId::new(0)));
    }
    for z in [0, 1] {
        let mut s = graph(&[z], &[0], &[]);
        s.topology.atoms[0].set_no_implicit(false);
        s.valence.implicit_hydrogens[0] = -1;
        let before = s.clone();
        close(chi_0_v(&s.input()).unwrap(), 0.0);
        assert!(chi_0_n(&s.input()).is_err());
        assert_eq!(s, before);
    }
}
#[test]
fn chi_noimplicit_and_prepared_implicit_plus_explicit() {
    let mut f = graph(&[6], &[3], &[]);
    f.valence.implicit_hydrogens[0] = -1;
    check(&f, [1.0, 0.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0, 0.0]);
    f.topology.atoms[0].set_no_implicit(false);
    f.valence.implicit_hydrogens[0] = 1;
    check(&f, [0.0; 5], [0.0; 5]);
}
#[test]
fn chi_signed_total_h_boundary_preserves_owner_cause() {
    // Atom.h stores d_implicitValence in int8_t. Atom::getNumImplicitHs
    // checks that stored value before getTotalNumHs performs its sum.
    // Retain both original i32::MAX inputs as exact source getter failures.
    let mut f = graph(&[6], &[0], &[]);
    f.topology.atoms[0].set_no_implicit(false);
    let assert_missing = |f: &Fixture| {
        let before = f.clone();
        let names = [
            "chi_0_v", "chi_1_v", "chi_2_v", "chi_3_v", "chi_4_v", "chi_0_n", "chi_1_n", "chi_2_n",
            "chi_3_n", "chi_4_n",
        ];
        for (entry, function) in V.into_iter().chain(N).zip(names) {
            let failure = entry(&f.input()).unwrap_err();
            let cause = ValenceError::ImplicitValenceCacheNotInitialized {
                atom: AtomId::new(0),
            };
            assert_eq!(
                failure,
                DescriptorError::Valence {
                    function,
                    source: cause.clone()
                }
            );
            assert_eq!(
                std::error::Error::source(&failure)
                    .unwrap()
                    .downcast_ref::<ValenceError>(),
                Some(&cause)
            );
            assert_eq!(f, &before);
        }
        for order in 2..=4 {
            for (actual, function) in [
                (chi_n_v(&f.input(), order), "chi_n_v"),
                (chi_n_n(&f.input(), order), "chi_n_n"),
            ] {
                assert_eq!(
                    actual.unwrap_err(),
                    DescriptorError::Valence {
                        function,
                        source: ValenceError::ImplicitValenceCacheNotInitialized {
                            atom: AtomId::new(0)
                        }
                    }
                );
                assert_eq!(f, &before);
            }
        }
    };
    for explicit in [0, 1] {
        f.topology.atoms[0].set_explicit_hydrogens(explicit);
        f.valence.implicit_hydrogens[0] = i32::MAX;
        assert_missing(&f);
    }
    for cached in [128, 255, -1] {
        f.valence.implicit_hydrogens[0] = cached;
        assert_missing(&f);
    }
    // Independent source arithmetic: unsigned (4 - 127), then sqrt/reciprocal.
    // Both calculation-width values have the same stored signed-byte value.
    f.topology.atoms[0].set_explicit_hydrogens(0);
    for cached in [127, i32::MAX - 128] {
        f.valence.implicit_hydrogens[0] = cached;
        check(
            &f,
            [0.000015258789280991895, 0.0, 0.0, 0.0, 0.0],
            [0.000015258789280991895, 0.0, 0.0, 0.0, 0.0],
        );
        for entry in [chi_0_v as Entry, chi_0_n as Entry] {
            assert_eq!(
                entry(&f.input()).unwrap().to_bits(),
                0.000015258789280991895_f64.to_bits()
            );
        }
    }
    f.topology.atoms[0].set_explicit_hydrogens(255);
    check(
        &f,
        [0.00001525878973396293, 0.0, 0.0, 0.0, 0.0],
        [0.00001525878973396293, 0.0, 0.0, 0.0, 0.0],
    );
    for entry in [chi_0_v as Entry, chi_0_n as Entry] {
        assert_eq!(
            entry(&f.input()).unwrap().to_bits(),
            0.00001525878973396293_f64.to_bits()
        );
    }
    f.topology.atoms[0].set_explicit_hydrogens(0);
    for cached in [0, 256] {
        f.valence.implicit_hydrogens[0] = cached;
        check(&f, [0.5, 0.0, 0.0, 0.0, 0.0], [0.5, 0.0, 0.0, 0.0, 0.0]);
    }
}
#[test]
fn chi_prepared_dimensions_and_invalid_topology() {
    for field in ["explicit_valence", "implicit_hydrogens"] {
        let mut f = graph(&[6], &[0], &[]);
        if field == "explicit_valence" {
            f.valence.explicit_valence.clear();
        } else {
            f.valence.implicit_hydrogens.clear();
        }
        errors_preserve(&f);
        assert!(
            matches!(chi_0_v(&f.input()),Err(DescriptorError::InvalidValenceRows { function:"chi_0_v",field:actual,actual:0,expected:1 }) if actual==field)
        );
    }
    let mut f = graph(&[6], &[0], &[]);
    f.topology.atoms[0] = f.topology.atoms[0].clone().with_id(AtomId::new(7));
    errors_preserve(&f);
    assert!(matches!(
        chi_0_v(&f.input()),
        Err(DescriptorError::InvalidTopology {
            function: "chi_0_v",
            ..
        })
    ));
}
#[test]
fn chi_supplied_ring_coordinate_property_state_not_reprepared() {
    let mut f = graph(&[6; 3], &[2; 3], &[(0, 1), (1, 2), (2, 0)]);
    f.rings = RingInfo::new(RingFindType::OtherOrUnknown, 0, 0); // Source Chi does not read rings.
    f.properties.set_computed_prop("unrelated", "42").unwrap();
    f.coordinates
        .conformers_2d
        .push(cosmolkit_model::Conformer2D::new(91, vec![[0.125, 0.5]; 3]));
    f.coordinates
        .conformers_3d
        .push(cosmolkit_model::Conformer3D::new(
            17,
            vec![[1.0, 2.0, 3.0]; 3],
            true,
        ));
    f.valence.explicit_valence = vec![7, 8, 9]; // Source Chi reads only implicit H rows.
    f.rings = RingInfo::new(RingFindType::OtherOrUnknown, 3, 3);
    f.rings.add_ring(&[0, 1, 2], &[0, 1, 2]).unwrap();
    check(
        &f,
        [
            2.1213203435596424,
            1.4999999999999996,
            1.060660171779821,
            0.3535533905932737,
            0.0,
        ],
        [
            2.1213203435596424,
            1.4999999999999996,
            1.060660171779821,
            0.3535533905932737,
            0.0,
        ],
    );
}
#[test]
fn chi_original_upstream97_expected_values() {
    for &(smiles, order, is_v, expected) in UPSTREAM {
        let params = cosmolkit_smiles::SmilesParseParams {
            remove_hydrogens: false,
            ..Default::default()
        };
        let record = cosmolkit_smiles::parse_smiles(smiles, &params).unwrap();
        let valence = cosmolkit_core::assign_valence_with_options_for_topology(
            &record.topology,
            cosmolkit_core::ValenceModel::RdkitLike,
            false,
        )
        .unwrap();
        let rings = RingInfo::new(
            RingFindType::OtherOrUnknown,
            record.topology.atoms.len(),
            record.topology.bonds.len(),
        );
        let f = Fixture {
            topology: record.topology,
            coordinates: CoordinateBlock::default(),
            properties: record.properties,
            valence,
            rings,
        };
        let before = f.clone();
        let value = if is_v {
            V[order](&f.input()).unwrap()
        } else {
            N[order](&f.input()).unwrap()
        };
        assert!(
            (value - expected).abs() < 0.002,
            "source original {smiles}/{order}/{is_v}: {value} vs {expected}"
        );
        assert_eq!(f, before);
    }
}
const SINGLETONS: &[(u8, f64, f64)] = &[
    (0, 0.0, 0.0),
    (1, 0.0, 1.0),
    (2, 0.7071067811865475, 0.7071067811865475),
    (3, 1.0, 1.0),
    (4, 0.7071067811865475, 0.7071067811865475),
    (5, 0.5773502691896258, 0.5773502691896258),
    (6, 0.5, 0.5),
    (7, 0.4472135954999579, 0.4472135954999579),
    (8, 0.4082482904638631, 0.4082482904638631),
    (9, 0.3779644730092272, 0.3779644730092272),
    (10, 0.35355339059327373, 0.35355339059327373),
    (11, 3.0, 1.0),
    (12, 2.1213203435596424, 0.7071067811865475),
    (13, 1.7320508075688774, 0.5773502691896258),
    (14, 1.5, 0.5),
    (15, 1.3416407864998738, 0.4472135954999579),
    (16, 1.224744871391589, 0.4082482904638631),
    (17, 1.1338934190276817, 0.3779644730092272),
    (18, 1.0606601717798212, 0.35355339059327373),
    (19, 4.123105625617661, 1.0),
    (20, 2.91547594742265, 0.7071067811865475),
    (21, 2.3804761428476167, 0.5773502691896258),
    (22, 2.0615528128088303, 0.5),
    (23, 1.8439088914585775, 0.4472135954999579),
    (24, 1.6832508230603465, 0.4082482904638631),
    (25, 1.5583874449479593, 0.3779644730092272),
    (26, 1.457737973711325, 0.35355339059327373),
    (27, 1.3743685418725535, 0.3333333333333333),
    (28, 1.3038404810405297, 0.31622776601683794),
    (29, 1.243163121016122, 0.30151134457776363),
    (30, 3.6742346141747673, 0.7071067811865475),
    (31, 3.0, 0.5773502691896258),
    (32, 2.598076211353316, 0.5),
    (33, 2.32379000772445, 0.4472135954999579),
    (34, 2.1213203435596424, 0.4082482904638631),
    (35, 1.9639610121239313, 0.3779644730092272),
    (36, 1.8371173070873836, 0.35355339059327373),
    (37, 5.916079783099616, 1.0),
    (38, 4.183300132670378, 0.7071067811865475),
    (39, 3.415650255319866, 0.5773502691896258),
    (40, 2.958039891549808, 0.5),
    (41, 2.6457513110645907, 0.4472135954999579),
    (42, 2.41522945769824, 0.4082482904638631),
    (43, 2.23606797749979, 0.3779644730092272),
    (44, 2.091650066335189, 0.35355339059327373),
    (45, 1.9720265943665387, 0.3333333333333333),
    (46, 1.8708286933869707, 0.31622776601683794),
    (47, 1.7837651700316897, 0.30151134457776363),
    (48, 4.743416490252569, 0.7071067811865475),
    (49, 3.8729833462074175, 0.5773502691896258),
    (50, 3.3541019662496843, 0.5),
    (51, 3.0, 0.4472135954999579),
    (52, 2.7386127875258306, 0.4082482904638631),
    (53, 2.5354627641855494, 0.3779644730092272),
    (54, 2.3717082451262845, 0.35355339059327373),
    (55, 7.280109889280518, 1.0),
    (56, 5.1478150704935, 0.7071067811865475),
    (57, 4.203173404306163, 0.5773502691896258),
    (58, 3.640054944640259, 0.5),
    (59, 4.281744192888377, 0.5773502691896258),
    (60, 3.7080992435478315, 0.5),
    (61, 3.3166247903554, 0.4472135954999579),
    (62, 3.0276503540974917, 0.4082482904638631),
    (63, 2.803059552906941, 0.3779644730092272),
    (64, 2.622022120425379, 0.35355339059327373),
    (65, 2.472066162365221, 0.3333333333333333),
    (66, 2.345207879911715, 0.31622776601683794),
    (67, 2.23606797749979, 0.30151134457776363),
    (68, 2.1408720964441885, 0.2886751345948129),
    (69, 2.0568833780186058, 0.2773500981126146),
    (70, 1.9820624179302297, 0.2672612419124244),
    (71, 1.9148542155126762, 0.2581988897471611),
    (72, 4.092676385936225, 0.5),
    (73, 3.6606010435446255, 0.4472135954999579),
    (74, 3.34165627596057, 0.4082482904638631),
    (75, 3.093772546815388, 0.3779644730092272),
    (76, 2.8939592256975564, 0.35355339059327373),
    (77, 2.7284509239574835, 0.3333333333333333),
    (78, 2.588435821108957, 0.31622776601683794),
    (79, 2.467976720090587, 0.30151134457776363),
    (80, 6.2048368229954285, 0.7071067811865475),
    (81, 5.066228051190222, 0.5773502691896258),
    (82, 4.387482193696061, 0.5),
    (83, 3.924283374069717, 0.4472135954999579),
    (84, 3.582364210034113, 0.4082482904638631),
    (85, 3.3166247903554, 0.3779644730092272),
    (86, 3.1024184114977142, 0.35355339059327373),
    (87, 9.219544457292887, 1.0),
    (88, 6.519202405202648, 0.7071067811865475),
    (89, 5.322906474223771, 0.5773502691896258),
    (90, 4.6097722286464435, 0.5),
    (91, 5.385164807134504, 0.5773502691896258),
    (92, 4.663689526544408, 0.5),
    (93, 4.171330722922842, 0.4472135954999579),
    (94, 3.8078865529319543, 0.4082482904638631),
    (95, 3.525417908358019, 0.3779644730092272),
    (96, 3.29772648956823, 0.35355339059327373),
    (97, 3.1091263510296048, 0.3333333333333333),
    (98, 2.9495762407505253, 0.31622776601683794),
    (99, 2.812310599683276, 0.30151134457776363),
    (100, 2.692582403567252, 0.2886751345948129),
    (101, 2.586949495507729, 0.2773500981126146),
    (102, 2.4928469095164494, 0.2672612419124244),
    (103, 2.4083189157584592, 0.2581988897471611),
    (104, 7.1063352017759485, 0.7071067811865475),
    (105, 7.14142842854285, 0.7071067811865475),
    (106, 7.176350047203662, 0.7071067811865475),
    (107, 7.211102550927978, 0.7071067811865475),
    (108, 7.245688373094719, 0.7071067811865475),
    (109, 7.280109889280518, 0.7071067811865475),
    (110, 7.3143694191638975, 0.7071067811865475),
    (111, 7.3484692283495345, 0.7071067811865475),
    (112, 7.382411530116699, 0.7071067811865475),
    (113, 7.416198487095663, 0.7071067811865475),
    (114, 7.44983221287567, 0.7071067811865475),
    (115, 7.483314773547883, 0.7071067811865475),
    (116, 7.516648189186454, 0.7071067811865475),
    (117, 7.54983443527075, 0.7071067811865475),
    (118, 7.582875444051551, 0.7071067811865475),
];
const UPSTREAM: &[(&str, usize, bool, f64)] = &[
    ("CCCCCC", 0, true, 4.828),
    ("CCC(C)CC", 0, true, 4.992),
    ("CC(C)CCC", 0, true, 4.992),
    ("CC(C)C(C)C", 0, true, 5.155),
    ("CC(C)(C)CC", 0, true, 5.207),
    ("CCCCCO", 0, true, 4.276),
    ("CCC(O)CC", 0, true, 4.439),
    ("CC(O)(C)CC", 0, true, 4.654),
    ("c1ccccc1O", 0, true, 3.834),
    ("CCCl", 0, true, 2.841),
    ("CCBr", 0, true, 3.671),
    ("CCI", 0, true, 4.242),
    ("CCCCCC", 1, true, 2.914),
    ("CCC(C)CC", 1, true, 2.808),
    ("CC(C)CCC", 1, true, 2.77),
    ("CC(C)C(C)C", 1, true, 2.643),
    ("CC(C)(C)CC", 1, true, 2.561),
    ("CCCCCO", 1, true, 2.523),
    ("CCC(O)CC", 1, true, 2.489),
    ("CC(O)(C)CC", 1, true, 2.284),
    ("c1ccccc1O", 1, true, 2.134),
    ("CCCCCC", 2, true, 1.707),
    ("CCC(C)CC", 2, true, 1.922),
    ("CC(C)CCC", 2, true, 2.183),
    ("CC(C)C(C)C", 2, true, 2.488),
    ("CC(C)(C)CC", 2, true, 2.914),
    ("CCCCCO", 2, true, 1.431),
    ("CCC(O)CC", 2, true, 1.47),
    ("CC(O)(C)CC", 2, true, 2.166),
    ("c1ccccc1O", 2, true, 1.336),
    ("CCCCCC", 3, true, 0.957),
    ("CCC(C)CC", 3, true, 1.394),
    ("CC(C)CCC", 3, true, 0.866),
    ("CC(C)C(C)C", 3, true, 1.333),
    ("CC(C)(C)CC", 3, true, 1.061),
    ("CCCCCO", 3, true, 0.762),
    ("CCC(O)CC", 3, true, 0.943),
    ("CC(O)(C)CC", 3, true, 0.865),
    ("c1ccccc1O", 3, true, 0.756),
    ("CCCCCC", 4, true, 0.5),
    ("CCC(C)CC", 4, true, 0.289),
    ("CC(C)CCC", 4, true, 0.577),
    ("CC(C)C(C)C", 4, true, 0.0),
    ("CC(C)(C)CC", 4, true, 0.0),
    ("CCCCCO", 4, true, 0.362),
    ("CCC(O)CC", 4, true, 0.289),
    ("CC(O)(C)CC", 4, true, 0.0),
    ("c1ccccc1O", 4, true, 0.428),
    ("CCCCCC", 0, false, 4.828),
    ("CCC(C)CC", 0, false, 4.992),
    ("CC(C)CCC", 0, false, 4.992),
    ("CC(C)C(C)C", 0, false, 5.155),
    ("CC(C)(C)CC", 0, false, 5.207),
    ("CCCCCO", 0, false, 4.276),
    ("CCC(O)CC", 0, false, 4.439),
    ("CC(O)(C)CC", 0, false, 4.654),
    ("c1ccccc1O", 0, false, 3.834),
    ("CCCl", 0, false, 2.085),
    ("CCBr", 0, false, 2.085),
    ("CCI", 0, false, 2.085),
    ("CCCCCC", 1, false, 2.914),
    ("CCC(C)CC", 1, false, 2.808),
    ("CC(C)CCC", 1, false, 2.77),
    ("CC(C)C(C)C", 1, false, 2.643),
    ("CC(C)(C)CC", 1, false, 2.561),
    ("CCCCCO", 1, false, 2.523),
    ("CCC(O)CC", 1, false, 2.489),
    ("CC(O)(C)CC", 1, false, 2.284),
    ("c1ccccc1O", 1, false, 2.134),
    ("C=S", 1, false, 0.289),
    ("CCCCCC", 2, false, 1.707),
    ("CCC(C)CC", 2, false, 1.922),
    ("CC(C)CCC", 2, false, 2.183),
    ("CC(C)C(C)C", 2, false, 2.488),
    ("CC(C)(C)CC", 2, false, 2.914),
    ("CCCCCO", 2, false, 1.431),
    ("CCC(O)CC", 2, false, 1.47),
    ("CC(O)(C)CC", 2, false, 2.166),
    ("c1ccccc1O", 2, false, 1.336),
    ("CCCCCC", 3, false, 0.957),
    ("CCC(C)CC", 3, false, 1.394),
    ("CC(C)CCC", 3, false, 0.866),
    ("CC(C)C(C)C", 3, false, 1.333),
    ("CC(C)(C)CC", 3, false, 1.061),
    ("CCCCCO", 3, false, 0.762),
    ("CCC(O)CC", 3, false, 0.943),
    ("CC(O)(C)CC", 3, false, 0.865),
    ("c1ccccc1O", 3, false, 0.756),
    ("CCCCCC", 4, false, 0.5),
    ("CCC(C)CC", 4, false, 0.289),
    ("CC(C)CCC", 4, false, 0.577),
    ("CC(C)C(C)C", 4, false, 0.0),
    ("CC(C)(C)CC", 4, false, 0.0),
    ("CCCCCO", 4, false, 0.362),
    ("CCC(O)CC", 4, false, 0.289),
    ("CC(O)(C)CC", 4, false, 0.0),
    ("c1ccccc1O", 4, false, 0.428),
];

#[test]
fn chi_narrow_input_borrows_only_topology_and_valence() {
    let fixture = graph(&[6; 4], &[3, 2, 2, 3], &[(0, 1), (1, 2), (2, 3)]);
    let before = fixture.clone();
    let input = ChiInput::new(&fixture.topology, &fixture.valence);
    assert!(std::ptr::eq(input.topology(), &fixture.topology));
    assert!(std::ptr::eq(input.valence(), &fixture.valence));
    // No coordinates/properties/ring carrier can enter this interface.
    for (entry, expected) in V
        .into_iter()
        .zip([
            3.414213562373095,
            1.914213562373095,
            0.9999999999999998,
            0.4999999999999999,
            0.0,
        ])
        .chain(N.into_iter().zip([
            3.414213562373095,
            1.914213562373095,
            0.9999999999999998,
            0.4999999999999999,
            0.0,
        ]))
    {
        close(entry(&input).unwrap(), expected);
    }
    close(chi_n_v(&input, 2).unwrap(), 0.9999999999999998);
    close(chi_n_n(&input, 2).unwrap(), 0.9999999999999998);
    assert_eq!(fixture, before);
}
